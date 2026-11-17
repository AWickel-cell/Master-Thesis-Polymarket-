library(dplyr)
library(tidyr)
library(lubridate)
library(readr)
library(ggplot2)
library(LaMa)

source("functions/ssm_functions.R")

# Load final model
final <- readRDS("results/final_model/results.rds")
grid <- final$grid

# Load Brent oil prices
oil <- read_csv(
  "https://fred.stlouisfed.org/graph/fredgraph.csv?id=DCOILBRENTEU",
  show_col_types = FALSE
) |>
  transmute(
    date = as.Date(observation_date),
    price = as.numeric(DCOILBRENTEU)
  ) |>
  filter(
    date >= as.Date("2025-11-01"),
    date <= as.Date("2026-11-10")
  ) |>
  arrange(date) |>
  mutate(
    oil_change = 100 * abs(log(price / lag(price)))
  )

# Use previous available oil price change
oil_daily <- oil |>
  filter(!is.na(oil_change)) |>
  transmute(
    date = date + 1,
    oil_change
  )

# Prepare hourly trading volume
model_data <- final$data |>
  mutate(date = as.Date(hour, tz = "UTC")) |>
  left_join(oil_daily, by = "date") |>
  arrange(hour) |>
  fill(oil_change, .direction = "down") |>
  filter(!is.na(oil_change))

cat("Hours:", nrow(model_data), "\n")
cat("Zero share:", round(mean(model_data$volume_usd == 0), 4), "\n")

# EDA
ggplot(oil, aes(date, oil_change)) +
  geom_line() +
  labs(
    title = "Absolute Daily Brent Oil Price Changes",
    x = NULL,
    y = "Absolute Log Return (%)"
  ) +
  theme_minimal()

ggplot(model_data, aes(oil_change, log1p(volume_usd))) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    title = "Oil Price Changes vs. Trading Volume",
    x = "Last Available Absolute Oil Price Change (%)",
    y = "Log(1 + Hourly Trading Volume)"
  ) +
  theme_minimal()

cat(
  "Correlation:",
  cor(model_data$oil_change, log1p(model_data$volume_usd)),
  "\n"
)

# Model data
y <- model_data$volume_usd

x_base <- as.matrix(
  model_data[c("sin_hour", "cos_hour", "weekend")]
)
z <- x_base

x_oil <- cbind(
  x_base,
  oil_change = model_data$oil_change
)

# Fit reference model
cat("\nFitting final SSM on oil sample...\n")

mod_ssm <- nlm(
  nll_ssm_state_zero,
  p = final$model$estimate,
  y = y,
  x = x_base,
  z = z,
  grid = grid,
  iterlim = 500
)

# Fit model with oil covariate
start_oil <- c(
  mod_ssm$estimate[1:8],
  0,
  mod_ssm$estimate[9:12]
)

cat("\nFitting SSM with oil covariate...\n")

mod_oil <- nlm(
  nll_ssm_state_zero,
  p = start_oil,
  y = y,
  x = x_oil,
  z = z,
  grid = grid,
  hessian = TRUE,
  iterlim = 500
)

# Model comparison
models <- list(
  "Final SSM" = mod_ssm,
  "Final SSM + Oil price change" = mod_oil
)

comparison <- data.frame(
  model = names(models),
  loglik = -sapply(models, function(m) m$minimum),
  parameters = sapply(models, function(m) length(m$estimate)),
  aic = sapply(models, function(m)
    2 * m$minimum + 2 * length(m$estimate)),
  bic = sapply(models, function(m)
    2 * m$minimum + log(length(y)) * length(m$estimate)),
  convergence = sapply(models, function(m) m$code),
  row.names = NULL
) |>
  mutate(
    delta_aic = aic - min(aic),
    delta_bic = bic - min(bic)
  )

# Oil coefficient and confidence interval
oil_beta <- mod_oil$estimate[9]

oil_se <- tryCatch({
  H <- mod_oil$hessian
  if (is.null(H) ||
      any(!is.finite(H)) ||
      min(eigen(H, symmetric = TRUE,
                only.values = TRUE)$values) <= 0) {
    NA_real_
  } else {
    sqrt(solve(H)[9, 9])
  }
}, error = function(e) NA_real_)

oil_effect <- data.frame(
  estimate = oil_beta,
  se = oil_se,
  lower = oil_beta - qnorm(0.975) * oil_se,
  upper = oil_beta + qnorm(0.975) * oil_se
)

cat("\nModel comparison:\n")
print(comparison)

cat("\nOil coefficient and 95% confidence interval:\n")
print(oil_effect)

# Save results
dir.create(
  "results/oil_price",
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  list(
    data = model_data,
    models = models,
    comparison = comparison,
    oil_effect = oil_effect
  ),
  "results/oil_price/results.rds"
)
