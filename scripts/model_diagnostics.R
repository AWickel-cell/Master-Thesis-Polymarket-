library(LaMa)
library(dplyr)
library(lubridate)

source("functions/ssm_functions.R")

# Load state-dependent zero model
fit <- readRDS("results/state_dependent_zero/results.rds")

model_data <- fit$data
grid <- fit$grid
par_hat <- fit$models[["State-dependent zeros"]]$estimate

y <- model_data$volume_usd

x <- as.matrix(
  model_data[c("price", "sin_hour", "cos_hour", "weekend")]
)

z <- as.matrix(
  model_data[c("sin_hour", "cos_hour", "weekend")]
)

stopifnot(
  length(par_hat) == 13,
  !anyNA(y), !anyNA(x), !anyNA(z),
  all(y >= 0)
)

# Model parameters
mu <- par_hat[2]
sigma <- exp(par_hat[3])
phi <- tanh(par_hat[4])
sigma_s <- exp(par_hat[5])
gamma_zero <- par_hat[13]

eta <- as.vector(mu + x %*% par_hat[6:9])
logit_pi <- as.vector(par_hat[1] + z %*% par_hat[10:12])

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
    cdf <- pnorm(
      log(y[t]),
      mean = eta[t] + grid,
      sd = sigma
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
cat("Mean:", round(mean(residuals), 4), "\n")
cat("SD:", round(sd(residuals), 4), "\n")
cat("Share |residual| > 3:",
    round(mean(abs(residuals) > 3), 4), "\n")

# Diagnostic plots
par(mfrow = c(2, 2))

hist(
  residuals,
  breaks = 50,
  probability = TRUE,
  main = "Randomized Quantile Residuals",
  xlab = "Residual"
)
curve(dnorm(x), add = TRUE, col = "red", lwd = 2)

qqnorm(residuals, main = "Normal Q-Q Plot")
qqline(residuals, col = "red", lwd = 2)

acf(
  residuals,
  lag.max = 72,
  main = "ACF of Quantile Residuals"
)

plot(
  model_data$hour,
  residuals,
  type = "l",
  xlab = "Time",
  ylab = "Quantile Residual",
  main = "Residuals Over Time"
)
abline(h = 0, col = "red")

par(mfrow = c(1, 1))

# Zero indicator diagnostics
acf(
  zero_residuals,
  lag.max = 72,
  main = "ACF of Zero Indicator Residuals"
)

cat("\nLjung-Box tests:\n")

print(Box.test(
  residuals,
  lag = 24,
  type = "Ljung-Box"
))

print(Box.test(
  residuals,
  lag = 72,
  type = "Ljung-Box"
))

print(Box.test(
  zero_residuals,
  lag = 24,
  type = "Ljung-Box"
))

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

par(mfrow = c(1, 2))

for (w in 0:1) {
  dat <- calibration |>
    filter(weekend == w)
  
  plot(
    dat$hour_utc,
    dat$estimated,
    type = "l",
    lwd = 2,
    ylim = range(c(dat$observed, dat$estimated)),
    xlab = "Hour (UTC)",
    ylab = "Probability of Zero Volume",
    main = if (w == 0) "Weekdays" else "Weekend"
  )
  
  points(dat$hour_utc, dat$observed, pch = 16)
  
  legend(
    "topright",
    legend = c("Estimated", "Observed"),
    lty = c(1, NA),
    pch = c(NA, 16),
    bty = "n"
  )
}

par(mfrow = c(1, 1))

# Save results
dir.create(
  "results/state_dependent_zero",
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  data.frame(
    hour = model_data$hour,
    volume_usd = y,
    residual = residuals,
    zero_residual = zero_residuals,
    pi_zero = pred_zero
  ),
  "results/state_dependent_zero/residuals.rds"
)
