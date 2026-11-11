library(LaMa)
library(dplyr)
library(lubridate)

source("functions/ssm_functions.R")

fit <- readRDS("results/model_fit.rds")

model_data <- fit$data
grid <- fit$grid
par_hat <- fit$model$estimate

y <- model_data$volume_usd
x <- as.matrix(model_data[c("price", "sin_hour", "cos_hour")])
z <- as.matrix(model_data[c("sin_hour", "cos_hour")])

mu <- par_hat[2]
sigma <- exp(par_hat[3])
phi <- tanh(par_hat[4])
sigma_s <- exp(par_hat[5])

eta <- as.vector(mu + x %*% par_hat[6:8])
pi_zero <- plogis(as.vector(
  par_hat[1] + z %*% par_hat[9:10]
))

gamma <- ar1_transition(phi, sigma_s, grid)
predicted <- as.vector(LaMa::stationary(gamma))

# One-step-ahead residuals
residuals <- rep(NA_real_, length(y))

for (t in seq_along(y)) {
  state_mean <- sum(predicted * grid)
  state_var <- sum(predicted * (grid - state_mean)^2)
  
  if (y[t] > 0) {
    residuals[t] <- (
      log(y[t]) - eta[t] - state_mean
    ) / sqrt(sigma^2 + state_var)
    
    likelihood <- dnorm(
      log(y[t]),
      mean = eta[t] + grid,
      sd = sigma
    ) / y[t]
    
    filtered <- predicted * likelihood
    filtered <- filtered / sum(filtered)
  } else {
    filtered <- predicted
  }
  
  predicted <- as.vector(filtered %*% gamma)
}

# Residual diagnostics
res <- na.omit(residuals)

cat("Mean:", round(mean(res), 4), "\n")
cat("SD:", round(sd(res), 4), "\n")
cat("Share |residual| > 3:",
    round(mean(abs(res) > 3), 4), "\n")

par(mfrow = c(2, 2))

hist(
  res,
  breaks = 50,
  probability = TRUE,
  main = "Residual Distribution",
  xlab = "Standardized Residuals"
)
curve(dnorm(x), add = TRUE, col = "red", lwd = 2)

qqnorm(res, main = "Normal Q-Q Plot")
qqline(res, col = "red", lwd = 2)

acf(
  res,
  lag.max = 72,
  main = "ACF of Residuals"
)

plot(
  model_data$hour,
  residuals,
  type = "l",
  xlab = "Time",
  ylab = "Standardized Residuals",
  main = "Residuals Over Time"
)
abline(h = 0, col = "red")

par(mfrow = c(1, 1))

# Zero probability by hour
hours <- 0:23

z_hours <- cbind(
  sin(2 * pi * hours / 24),
  cos(2 * pi * hours / 24)
)

estimated <- plogis(as.vector(
  par_hat[1] + z_hours %*% par_hat[9:10]
))

observed <- model_data |>
  mutate(hour_utc = hour(hour)) |>
  group_by(hour_utc) |>
  summarise(
    zero_share = mean(volume_usd == 0),
    .groups = "drop"
  )

plot(
  hours, estimated,
  type = "l",
  lwd = 2,
  ylim = range(c(estimated, observed$zero_share)),
  xlab = "Hour (UTC)",
  ylab = "Probability of Zero Volume",
  main = "Observed vs Estimated Zero Probability"
)

points(
  observed$hour_utc,
  observed$zero_share,
  pch = 16
)

legend(
  "topright",
  legend = c("Estimated", "Observed"),
  lty = c(1, NA),
  pch = c(NA, 16),
  bty = "n"
)

# Save results
saveRDS(
  data.frame(
    hour = model_data$hour,
    volume_usd = y,
    residual = residuals
  ),
  "results/model_residuals.rds"
)
