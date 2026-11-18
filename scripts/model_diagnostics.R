library(LaMa)
library(dplyr)
library(lubridate)

source("functions/ssm_functions.R")

# Load final model
fit <- readRDS("results/final_model/results.rds")

model_data <- fit$data
grid <- fit$grid
par_hat <- fit$model$estimate

y <- model_data$volume_usd

x <- as.matrix(
  model_data[c("sin_hour", "cos_hour", "weekend")]
)

z <- x

stopifnot(
  length(par_hat) == 12,
  !anyNA(y), !anyNA(x), !anyNA(z),
  all(y >= 0)
)

# Model parameters
mu <- par_hat[2]
sigma <- exp(par_hat[3])
phi <- tanh(par_hat[4])
sigma_s <- exp(par_hat[5])
gamma_zero <- par_hat[12]

eta <- as.vector(mu + x %*% par_hat[6:8])
logit_pi <- as.vector(par_hat[1] + z %*% par_hat[9:11])

gamma <- ar1_transition(phi, sigma_s, grid)
predicted <- as.vector(LaMa::stationary(gamma))

# One-step-ahead quantile residuals
set.seed(123)

residuals <- numeric(length(y))
pred_zero <- numeric(length(y))
zero_residuals <- numeric(length(y))

for (t in seq_along(y)) {
  
  pi_state <- plogis(logit_pi[t] + gamma_zero * grid)
  
  pred_zero[t] <- sum(predicted * pi_state)
  
  if (y[t] == 0) {
    lower <- 0
    upper <- pred_zero[t]
    likelihood <- pi_state
  } else {
    cdf <- plnorm(
      y[t],
      meanlog = eta[t] + grid,
      sdlog = sigma
    )
    
    upper <- pred_zero[t] +
      sum(predicted * (1 - pi_state) * cdf)
    
    lower <- upper
    
    likelihood <- (1 - pi_state) *
      dlnorm(
        y[t],
        meanlog = eta[t] + grid,
        sdlog = sigma
      )
  }
  
  u <- runif(1, lower, upper)
  
  residuals[t] <- qnorm(
    pmin(pmax(u, 1e-10), 1 - 1e-10)
  )
  
  zero_residuals[t] <- as.integer(y[t] == 0) - pred_zero[t]
  
  # Filter latent state
  filtered <- predicted * likelihood
  filtered <- filtered / sum(filtered)
  
  predicted <- as.vector(filtered %*% gamma)
}

# Residual summary
cat("\n--- Residual Diagnostics ---\n")
cat("Observations:", length(y), "\n")
cat("Mean:", round(mean(residuals), 4), "\n")
cat("SD:", round(sd(residuals), 4), "\n")
cat("Median:", round(median(residuals), 4), "\n")
cat("Share |residual| > 2:",
    round(mean(abs(residuals) > 2), 4), "\n")
cat("Share |residual| > 3:",
    round(mean(abs(residuals) > 3), 4), "\n")

cat("\n--- Zero Volume Calibration ---\n")
cat("Observed zero share:",
    round(mean(y == 0), 4), "\n")
cat("Predicted zero share:",
    round(mean(pred_zero), 4), "\n")
cat("Mean zero residual:",
    round(mean(zero_residuals), 4), "\n")

acf_values <- acf(
  residuals,
  lag.max = 72,
  plot = FALSE
)$acf

cat("\n--- Residual Autocorrelation ---\n")
cat("ACF lag 1:", round(acf_values[2], 4), "\n")
cat("ACF lag 24:", round(acf_values[25], 4), "\n")
cat("ACF lag 72:", round(acf_values[73], 4), "\n")


# Output directory
dir.create(
  "results/diagnostics",
  recursive = TRUE,
  showWarnings = FALSE
)

# Helper: show and save each plot
show_and_save <- function(filename, plot_fun,
                          width = 900, height = 550) {
  plot_fun()
  
  png(
    file.path("results/diagnostics", filename),
    width = width,
    height = height
  )
  plot_fun()
  dev.off()
}

# 1. Histogram
show_and_save("residual_histogram.png", function() {
  hist(
    residuals,
    breaks = 50,
    probability = TRUE,
    main = "Randomized Quantile Residuals",
    xlab = "Residual"
  )
  curve(dnorm(x), add = TRUE, col = "red", lwd = 2)
})

# 2. Normal Q-Q plot
show_and_save("residual_qqplot.png", function() {
  qqnorm(residuals, main = "Normal Q-Q Plot")
  qqline(residuals, col = "red", lwd = 2)
})

# 3. Residual ACF
show_and_save("residual_acf.png", function() {
  acf(
    residuals,
    lag.max = 72,
    main = "ACF of Quantile Residuals"
  )
})

# 4. Residuals over time
show_and_save("residual_timeseries.png", function() {
  plot(
    model_data$hour,
    residuals,
    type = "l",
    xlab = "Time",
    ylab = "Quantile Residual",
    main = "Residuals Over Time"
  )
  abline(h = 0, col = "red")
}, width = 1100)

# 5. Zero indicator ACF
show_and_save("zero_acf.png", function() {
  acf(
    zero_residuals,
    lag.max = 72,
    main = "ACF of Zero Indicator Residuals"
  )
})

# Zero probability calibration
calibration <- model_data |>
  mutate(
    pi_zero = pred_zero,
    hour_utc = hour(hour)
  ) |>
  group_by(weekend, hour_utc) |>
  summarise(
    observed = mean(volume_usd == 0),
    estimated = mean(pi_zero),
    .groups = "drop"
  )

calibration_summary <- model_data |>
  mutate(pi_zero = pred_zero) |>
  group_by(weekend) |>
  summarise(
    observed = mean(volume_usd == 0),
    estimated = mean(pi_zero),
    .groups = "drop"
  )

cat("\n--- Zero Calibration by Weekend ---\n")
print(calibration_summary)

cat("\n--- Zero Calibration by Hour ---\n")
print(calibration, n = Inf)

# 6. Zero calibration: weekdays
show_and_save("zero_calibration_weekdays.png", function() {
  dat <- calibration |> filter(weekend == 0)
  
  plot(
    dat$hour_utc,
    dat$estimated,
    type = "l",
    lwd = 2,
    ylim = range(c(dat$observed, dat$estimated)),
    xlab = "Hour (UTC)",
    ylab = "Probability of Zero Volume",
    main = "Zero Calibration - Weekdays"
  )
  
  points(dat$hour_utc, dat$observed, pch = 16)
  
  legend(
    "topright",
    legend = c("Estimated", "Observed"),
    lty = c(1, NA),
    pch = c(NA, 16),
    bty = "n"
  )
})

# 7. Zero calibration: weekend
show_and_save("zero_calibration_weekend.png", function() {
  dat <- calibration |> filter(weekend == 1)
  
  plot(
    dat$hour_utc,
    dat$estimated,
    type = "l",
    lwd = 2,
    ylim = range(c(dat$observed, dat$estimated)),
    xlab = "Hour (UTC)",
    ylab = "Probability of Zero Volume",
    main = "Zero Calibration - Weekend"
  )
  
  points(dat$hour_utc, dat$observed, pch = 16)
  
  legend(
    "topright",
    legend = c("Estimated", "Observed"),
    lty = c(1, NA),
    pch = c(NA, 16),
    bty = "n"
  )
})


# Save results
saveRDS(
  data.frame(
    hour = model_data$hour,
    volume_usd = y,
    residual = residuals,
    zero_residual = zero_residuals,
    pi_zero = pred_zero
  ),
  "results/diagnostics/residuals.rds"
)
