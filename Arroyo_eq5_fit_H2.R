library(tidyverse)

# Import data
df_h2 <- read_tsv("./data/h2-bulk.tsv")

# Define which columns are soil rate data
soils <- c("RR1","RR2","RR3",
           "DML1","DML2","DML3",
           "BH1","BH2","BH3")

RR <- c("RR1","RR2","RR3")
DML <- c("DML1","DML2","DML3")
BH <- c("BH1","BH2","BH3")

# Reshape into long format and filter bad data
df_long <- df_h2 %>%
  pivot_longer(cols = all_of(soils),
               names_to = "soil",
               values_to = "rate_nmol") %>%
  mutate(region = case_when(
    soil %in% RR ~ "RR",
    soil %in% DML ~ "DML",
    soil %in% BH ~ "BH")
  ) %>%
  filter(rate_nmol != 0) %>%
  drop_na(rate_nmol, TempK) %>%
  mutate(
    invT = 1 / TempK,
    log_invT = log(invT),
    log_rate = log(rate_nmol)
  )

#Define the equation 5
eq5 <- function(invT, log_invT, log_y0, b, a) {
  log_y0 - b * invT - a * log_invT
}

# Temperature range for predictions
T_seq <- tibble(TempK = seq(min(df_long$TempK), max(df_long$TempK), length.out = 200)) %>%
  mutate(
    invT = 1 / TempK,
    log_invT = log(invT)
  )

set.seed(123)

# Function to bootstrap one iteration
bootstrap_once <- function(df_long) {
  df_boot <- df_long %>% slice_sample(prop = 1, replace = TRUE)
  fit <- tryCatch({
    nls(
      log_rate ~ eq5(invT, log_invT, log_y0, b, a),
      data = df_boot,
      start = list(log_y0 = -10, b = 5000, a = 1),
      control = nls.control(maxiter = 100, warnOnly = TRUE)
    )
  }, error = function(e) NULL)
  
  if (!is.null(fit)) {
    predict(fit, newdata = T_seq)
  } else {
    rep(NA, nrow(T_seq))
  }
}

# Run 1000 bootstraps
boot_mat <- replicate(1000, bootstrap_once(df_long))

# Clean matrix: remove failed fits (all-NA columns)
boot_mat <- boot_mat[, colSums(is.na(boot_mat)) == 0]

# Summarise bootstrap predictions
boot_summary <- T_seq %>%
  mutate(
    log_rate_median = apply(boot_mat, 1, median, na.rm = TRUE),
    log_rate_lower  = apply(boot_mat, 1, quantile, probs = 0.025, na.rm = TRUE),
    log_rate_upper  = apply(boot_mat, 1, quantile, probs = 0.975, na.rm = TRUE)
  )

#Plot
eq5_plot <- ggplot() +
  geom_ribbon(data = boot_summary, aes(x = 1 / TempK, ymin = log_rate_lower, ymax = log_rate_upper),
              fill = "skyblue", alpha = 0.3) +
  geom_point(data = df_long, aes(x = 1 / TempK, y = log_rate, color = soil, shape = soil), size = 2) +
  scale_color_manual(values = c("BH1" = "#E8AEFB", "BH2" = "#BA55D3", "BH3" = "#800080",
                                "DML1" = "#0F99B2", "DML2" = "#1E90FF", "DML3" = "#0000FF",
                                "RR1" = "#FA98C7", "RR2" = "#F71480", "RR3" = "#FFA500")) +
  scale_shape_manual(values = c(
    "RR1" = 16, "RR2" = 16, "RR3" = 16,
    "DML1" = 17, "DML2" = 17, "DML3" = 17,
    "BH1" = 15, "BH2" = 15, "BH3" = 15)) +
  geom_line(data = boot_summary, aes(x = 1 / TempK, y = log_rate_median), color = "darkblue", linewidth = 1) +
  scale_y_continuous(
    limits = c(-10, 2.5),
    breaks = seq(-10, 2.5, by = 2)
  ) +
  scale_x_continuous(
    limits = c(0.0029, 0.0038),
    breaks = seq(0.0028, 0.0038, by = 0.0002)
  ) +
  labs(
    x = "1/K",
    y = "ln(rate)") +
  guides(
    shape = guide_legend(title = "Soil"),
    color = guide_legend(title = "Soil")
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),   # remove minor gridlines
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5), # add border around each plot
    panel.background = element_blank(),  # empty panel background
    plot.background  = element_blank(),   # empty plot background
    axis.ticks.x = element_line(),
    axis.ticks.y = element_line(),
    legend.position = "right"
  )

ggsave("./output/arroyo_eq5_H2_boot1000_95CI.svg", eq5_plot, units = "cm", width = 12, height = 10)



##########Assessing the goodness of fit

#Simple nonlinear fit
fit_eq5 <- nls(
  log_rate ~ eq5(invT, log_invT, log_y0, b, a),
  data = df_long,
  start = list(log_y0 = -10, b = 5000, a = 1)
)

# Calculate total and residual sum of squares to give an R2-like value (even for non-linear models)
log_rate_obs <- df_long$log_rate
log_rate_pred <- predict(fit_eq5)

SStot <- sum((log_rate_obs - mean(log_rate_obs))^2)
SSres <- sum((log_rate_obs - log_rate_pred)^2)

pseudo_R2 <- 1 - (SSres / SStot)

#Residuals plots
df_long$resid <- residuals(fit_eq5)

res_plot <- ggplot(df_long, aes(x = predict(fit_eq5), y = resid)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Predicted log(rate)", y = "Residual") +
  theme_minimal()

ggsave("./output/arroyo_eq5_H2_residuals.svg", res_plot)

#Check the normality of the residuals
qqnorm(residuals(fit_eq5))
qqline(residuals(fit_eq5))
#Saved this plot direct from the output

#Generate AIC value
aic_eq5 <- AIC(fit_eq5)

#Calculate bootstrapped RMSE
rmse <- sqrt(mean((log_rate_obs - log_rate_pred)^2))

#Create a text file with various goodness of fit measures
sink("./output/arroyo_eq5_H2_fit_check.txt")
print("Pseudo R2 value is... ")
pseudo_R2
print("AIC value")
aic_eq5
print("Bootstrapped RMSE")
rmse
sink()


