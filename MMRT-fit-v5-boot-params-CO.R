library(tidyverse)
library(ggplot2)
library(boot)
library(minpack.lm)
library(purrr)
library(svglite)
library(numDeriv)

# Read in data
df <- read_tsv("./data/CO-bulk.tsv")

# Define which columns are soil rate data
soils <- c("RR1","RR2","RR3",
           "DML1","DML2","DML3",
           "BH1","BH2","BH3")

# Reshape into long format and filter bad data
df_long <- df %>%
  pivot_longer(cols = all_of(soils),
               names_to = "soil",
               values_to = "rate_nmol") %>%
  filter(rate_nmol != 0) %>%
  drop_na(rate_nmol, TempK)

# Constants for MMRT
R  <- 8.314
T0 <- 298.15  # Reference temperature (25°C in Kelvin)

# MMRT equation
MMRT <- function(T, H, S, Cp) {
  log((1.38e-23*T) / 6.626e-34) +
    (S/R) - (H/(R*T)) +
    (Cp/R) * (1 - T0/T + log(T/T0))
}

set.seed(123)
# Function to generate predicted values from fitted model
predict_mmrt <- function(fit, temps) {
  pars <- coef(fit)
  MMRT(temps, H = pars["H"], S = pars["S"], Cp = pars["Cp"])
}

# Function to bootstrap MMRT fits and calculate 95% CI
bootstrap_CI <- function(df, fit_obj, n_boot = 1000) {
  temps <- df$TempK
  # Original fitted coefficients
  coefs <- coef(fit_obj)
  
  # Bootstrap residuals and refit model
  boot_preds <- replicate(n_boot, {
    # Resample residuals
    resampled <- residuals(fit_obj)[sample(1:length(residuals(fit_obj)), replace = TRUE)]
    # Create new "bootstrapped" response
    new_response <- fitted(fit_obj) + resampled
    # Refit model with new response
    tryCatch({
      new_fit <- nlsLM(new_response ~ MMRT(temps, H, S, Cp),
                       start = as.list(coefs),
                       control = nls.lm.control(maxiter = 100))
      predict_mmrt(new_fit, temps)
    }, error = function(e) rep(NA, length(temps)))
  })
  
  # Remove columns with NA (failed fits)
  boot_preds <- boot_preds[, colSums(is.na(boot_preds)) == 0]
  
  # Calculate 2.5% and 97.5% quantiles across replicates
  ci_lower <- apply(boot_preds, 1, quantile, 0.025, na.rm = TRUE)
  ci_upper <- apply(boot_preds, 1, quantile, 0.975, na.rm = TRUE)
  
  tibble(TempK = temps,
         log_rate_lower = ci_lower,
         log_rate_upper = ci_upper)
}

##### Fit MMRT for each soil, extract predictions and CI #####
nested_results <- df_long |>
  nest(.by = soil) |>
  mutate(
    # Fit the MMRT model to each soil
    fit = map(data, \(df) nlsLM(log(rate_nmol) ~ MMRT(TempK, H, S, Cp),
                                data = df,
                                start = list(H = 5e4, S = 100, Cp = -1000))),
    
    # Get fitted values (on log scale)
    log_rate_pred = map2(fit, data, \(mod, df) predict(mod, tempseq = df)),
    
    # Exponentiate to return predictions to original scale
    rate_pred = map(log_rate_pred, exp),
    
    # Get 95% confidence intervals for each fit
    CI = map2(data, fit, bootstrap_CI)
  )

# Merge fitted values + confidence intervals back with original data
# First combine each CI table with its data + predictions
df_long_new <- nested_results |>
  mutate(data_with_pred = pmap(list(data, log_rate_pred, rate_pred, CI), \(df, logp, p, ci) {
    df |>
      mutate(log_rate_pred = logp,
             rate_pred = p) |>
      bind_cols(ci)  # combine CI columns safely (same row order)
  })) |>
  select(soil, data_with_pred) |>
  unnest(data_with_pred) |>
  mutate(rate_lower = exp(log_rate_lower),
         rate_upper = exp(log_rate_upper)) |>
  select(-`TempK...5`) |>
  rename(TempK = `TempK...1`) |>
  mutate(TempC = TempK - 273.15)


# Save the result
write_csv(df_long_new, "./output/MMRTfits_withCI_CO.txt")

##### Plotting #####

#Make the faceted plot
df_long_new |>
  # Create grid row/column for prettier facetting
  mutate(row = case_when(grepl("RR", soil) ~ "RR",
                         grepl("DML", soil) ~ "DML",
                         grepl("BH", soil) ~ "BH"),
         cols = case_when(grepl("1", soil) ~ "1",
                          grepl("2", soil) ~ "2",
                          grepl("3", soil) ~ "3")) |>
  ggplot(aes(TempC, rate_nmol)) +
  geom_point(size = 3) +  # raw data
  geom_ribbon(aes(ymin = rate_lower, ymax = rate_upper), 
              fill = "grey80", alpha = 0.5) +  # shaded 95% CI
  geom_smooth(aes(y = rate_pred), colour = "black", linewidth = 1) +  # fitted line
  scale_y_log10(
    limits = c(10^-5, 10),
    breaks = 10^(-4:1),
    labels = parse(text = paste0("10^", -4:1))
  ) +
  scale_x_continuous(
    limits = c(-40, 80),
    breaks = seq(-40, 80, by = 20),
    name = "Temperature (˚C)"
  ) +
  labs(
    y = expression(CO~oxidation~rate~"("*nmol~g^{-1}~h^{-1}*")"),
    x = expression("Temperature ("*~degree*C*")"),
    title = "MMRT fit for bulk-soil CO oxidation\nwith 95% confidence intervals"
  ) +
  theme_minimal(base_size = 14) +  # clean white background
  theme(
    panel.grid.major = element_blank(),   # remove major gridlines
    panel.grid.minor = element_blank(),   # remove minor gridlines
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5), # add border around each plot
    panel.background = element_blank(),  # empty panel background
    plot.background  = element_blank(),   # empty plot background
    axis.ticks.x = element_line(),
    axis.ticks.y = element_line()
  ) +
  facet_grid(row ~ cols)

ggsave("./output/MMRTfit_facet_CI_CO.svg")

####################################################
##### Calculate Topt and Cp for each soil #####
####################################################

# Function to calculate Topt numerically given H, S, Cp
calc_Topt <- function(H, S, Cp, T0 = 298.15) {
  f <- function(T) {
    S - Cp * log(T / T0) - (H - Cp * (T - T0)) / T
  }
  # Try root finding with error catching
  root <- tryCatch({
    uniroot(f, interval = c(250, 350))$root
  }, error = function(e) NA_real_)
  return(root)
}

nested_results <- nested_results %>%
  mutate(
    # Extract model coefficients as named list
    coef_list = map(fit, coef),
    
    # Calculate Topt for each soil's fit using calc_Topt()
    Topt = map_dbl(coef_list, ~ calc_Topt(H = .x["H"], S = .x["S"], Cp = .x["Cp"], T0 = T0)),
    
    # Extract Cp directly (ΔCp‡)
    Cp = map_dbl(coef_list, ~ .x["Cp"])
  )

# Function to bootstrap one soil dataset (Topt calculated numerically)
bootstrap_Topt_Cp <- function(df, n_boot = 1000, T0 = 298.15) {
  temps <- df$TempK
  rates_log <- log(df$rate_nmol)
  
  # Fit original model
  fit_orig <- nlsLM(rates_log ~ MMRT(temps, H, S, Cp),
                    start = list(H = 5e4, S = 100, Cp = -1000),
                    control = nls.lm.control(maxiter = 100))
  
  residuals_orig <- residuals(fit_orig)
  fitted_orig <- fitted(fit_orig)
  pars_orig <- coef(fit_orig)
  
  # Numerical Topt: find max of predicted curve
  temp_grid <- seq(250, 350, by = 0.1)
  preds_orig <- MMRT(temp_grid, pars_orig["H"], pars_orig["S"], pars_orig["Cp"])
  Topt_numeric <- temp_grid[which.max(preds_orig)]
  
  # Storage
  boot_params <- matrix(NA, nrow = n_boot, ncol = 3)
  colnames(boot_params) <- c("H", "S", "Cp")
  Topt_boot <- numeric(n_boot)
  
  for (i in seq_len(n_boot)) {
    resampled_resid <- sample(residuals_orig, replace = TRUE)
    new_response <- fitted_orig + resampled_resid
    
    try({
      fit_boot <- nlsLM(new_response ~ MMRT(temps, H, S, Cp),
                        start = as.list(pars_orig),
                        control = nls.lm.control(maxiter = 100))
      
      boot_pars <- coef(fit_boot)
      boot_params[i, ] <- boot_pars
      
      preds_boot <- MMRT(temp_grid, boot_pars["H"], boot_pars["S"], boot_pars["Cp"])
      Topt_boot[i] <- temp_grid[which.max(preds_boot)]
    }, silent = TRUE)
  }
  
  # Filter
  valid_rows <- complete.cases(boot_params) & !is.na(Topt_boot)
  boot_params <- boot_params[valid_rows, , drop = FALSE]
  Topt_boot <- Topt_boot[valid_rows]
  
  # Confidence intervals
  Topt_CI <- quantile(Topt_boot, probs = c(0.025, 0.975), na.rm = TRUE)
  Cp_CI <- quantile(boot_params[, "Cp"], probs = c(0.025, 0.975), na.rm = TRUE)
  
  list(
    Topt_mean = Topt_numeric,
    Topt_CI_lower = Topt_CI[1],
    Topt_CI_upper = Topt_CI[2],
    TempC = Topt_numeric - 273.15,
    Cp_mean = pars_orig["Cp"],
    Cp_CI_lower = Cp_CI[1],
    Cp_CI_upper = Cp_CI[2]
  )
}
  
# Apply bootstrap per soil
bootstrap_results <- nested_results %>%
  mutate(
    bootstrap = map(data, ~ bootstrap_Topt_Cp(.x, n_boot = 1000, T0 = T0))
  ) %>%
  unnest_wider(bootstrap)

# Format output
soil_bootstrap_summary <- bootstrap_results %>%
  select(soil, Topt_mean, Topt_CI_lower, Topt_CI_upper, TempC,
         Cp_mean, Cp_CI_lower, Cp_CI_upper)

# Save
write_tsv(soil_bootstrap_summary, "./output/MMRT_Topt_numerical_Cp_CO.tsv")



####################################################
##### Finding mean MMRT fit across all 9 sites #####
####################################################

# Your existing pooled fit:
fit.all <- nlsLM(log(rate_nmol) ~ MMRT(TempK, H, S, Cp),
                 data = df_long,
                 start = list(H = 5e4, S = 100, Cp = -1000))

# Get model parameters and variance-covariance matrix
params <- coef(fit.all)
vcov_mat <- vcov(fit.all)

# Create a prediction grid
tempseq <- tibble(
  TempK = seq(253.15, 353.15, by = 1),  # -20 to 80°C in Kelvin
  TempC = TempK - 273.15
)

# Wrapper to compute fit and confidence intervals
MMRT_pred <- function(par, T) {
  MMRT(T, par["H"], par["S"], par["Cp"])
}

get_confidence_interval <- function(Tval) {
  # Compute gradient of MMRT wrt params at Tval
  grad <- grad(func = function(p) MMRT_pred(setNames(p, names(params)), Tval),
               x = unname(params))
  
  # Variance of prediction using delta method
  pred_var <- t(grad) %*% vcov_mat %*% grad
  se_fit <- sqrt(pred_var)
  
  # Predicted log-rate value
  fit_val <- MMRT_pred(params, Tval)
  
  # Return as 1-row tibble with named columns!
  tibble(
    log_rate_pred  = fit_val,
    log_rate_lower = fit_val - 1.96 * se_fit,
    log_rate_upper = fit_val + 1.96 * se_fit
  )
}


# Apply across all temps, result is a tidy tibble
confint_df <- map_dfr(tempseq$TempK, get_confidence_interval)

# Combine with original tempseq and back-transform
tempseq <- tempseq %>%
  bind_cols(confint_df) %>%
  mutate(
    rate_pred  = exp(log_rate_pred),
    rate_lower = exp(log_rate_lower),
    rate_upper = exp(log_rate_upper)
  )

# -------------------------------------------
# Plot: pooled fit with 95% confidence ribbon
# -------------------------------------------
ggplot(df_long_new, aes(TempC, rate_nmol)) +
  #  geom_point(aes(color = soil), size = 2, alpha = 0.6) +
  geom_ribbon(data = tempseq, aes(x = TempC, ymin = rate_lower, ymax = rate_upper),
              fill = "grey80", alpha = 0.5, inherit.aes = FALSE) +
  geom_line(data = tempseq, aes(x = TempC, y = rate_pred),
            color = "black", linewidth = 0.8, inherit.aes = FALSE) +
  scale_y_log10(
    breaks = 10^(-5:1),
    labels = parse(text = paste0("10^", -5:1)),
    limits = c(1e-5, 10)
  ) +
  scale_x_continuous(
    limits = c(-40, 80),
    breaks = seq(-40, 80, by = 20),
    name = "Temperature (˚C)"
  ) +
  labs(
    title = "Pooled CO MMRT fit with 95% confidence interval",
    x = expression("Temperature ("*degree*C*")"),
    y = expression(CO~oxidation~rate~"("*nmol~g^{-1}~h^{-1}*")")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),   # remove major gridlines
    panel.grid.minor = element_blank(),   # remove minor gridlines
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5), # add border around each plot
    panel.background = element_blank(),  # empty panel background
    plot.background  = element_blank(),   # empty plot background
    axis.ticks.x = element_line(),
    axis.ticks.y = element_line()
  )

ggsave("./output/MMRTfit_Pooled_CO_95CI.svg", width = 11.16, height = 8.65, units = "cm")
