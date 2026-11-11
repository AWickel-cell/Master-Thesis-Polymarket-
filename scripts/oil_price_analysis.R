library(dplyr)
library(tidyr)
library(lubridate)
library(readr)
library(ggplot2)

source("functions/ssm_functions.R")

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

# Load hourly trading volume
model_data <- read_csv(
  "data/volume_hourly.csv",
  show_col_types = FALSE
) |>
  mutate(
    hour = as_datetime(hour, tz = "UTC"),
    date = as.Date(hour, tz = "UTC")
  ) |>
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

# Fit models
y <- model_data$volume_usd
x <- as.matrix(model_data["oil_change"])
grid <- seq(-5, 5, length.out = 151)

positive_log_y <- log(y[y > 0])

start_ssm <- c(
  qlogis(mean(y == 0)),
  mean(positive_log_y),
  log(sd(positive_log_y)),
  atanh(0.7),
  log(0.3)
)

mod_ssm <- nlm(
  nll_ssm, start_ssm,
  y = y, grid = grid
)

mod_oil <- nlm(
  nll_ssm_cov,
  c(mod_ssm$estimate, 0),
  y = y, x = x, grid = grid
)

models <- list(
  "Gaussian SSM" = mod_ssm,
  "SSM + Oil price change" = mod_oil
)

comparison <- data.frame(
  model = names(models),
  loglik = -sapply(models, function(m) m$minimum),
  aic = sapply(models, function(m)
    2 * m$minimum + 2 * length(m$estimate)),
  bic = sapply(models, function(m)
    2 * m$minimum + log(length(y)) * length(m$estimate)),
  convergence = sapply(models, function(m) m$code),
  row.names = NULL
)

cat("Oil coefficient:", round(mod_oil$estimate[6], 4), "\n")
print(comparison)

# Save results
dir.create("results/oil_price", recursive = TRUE, showWarnings = FALSE)

saveRDS(
  list(
    data = model_data,
    models = models,
    comparison = comparison
  ),
  "results/oil_price/results.rds"
)
