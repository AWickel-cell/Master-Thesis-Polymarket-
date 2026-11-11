library(dplyr)
library(tidyr)
library(lubridate)
library(readr)

source("functions/ssm_functions.R")

# Prepare data
trades <- read_csv(
  "data/trades_raw_v2.csv",
  col_types = cols(
    proxy_wallet = col_character(),
    token_id = col_character(),
    condition_id = col_character(),
    transaction_hash = col_character(),
    .default = col_guess()
  ),
  show_col_types = FALSE
)

yes_prices <- trades |>
  filter(outcome == "Yes") |>
  mutate(
    datetime = as_datetime(timestamp, tz = "UTC"),
    hour = floor_date(datetime, "hour")
  ) |>
  arrange(datetime) |>
  group_by(hour) |>
  summarise(price_yes = last(price), .groups = "drop")

model_data <- read_csv(
  "data/volume_hourly.csv",
  show_col_types = FALSE
) |>
  left_join(yes_prices, by = "hour") |>
  arrange(hour) |>
  fill(price_yes, .direction = "down") |>
  mutate(
    price_lag = lag(price_yes),
    price_change = abs(price_lag - lag(price_lag)),
    price = 100 * price_lag,
    change = 100 * price_change,
    sin_hour = sin(2 * pi * hour(hour) / 24),
    cos_hour = cos(2 * pi * hour(hour) / 24)
  ) |>
  filter(!is.na(price), !is.na(change))

y <- model_data$volume_usd
grid <- seq(-5, 5, length.out = 151)

x_price <- as.matrix(model_data["price"])
x_change <- as.matrix(model_data["change"])
x_time <- as.matrix(model_data[c("sin_hour", "cos_hour")])
x_both <- as.matrix(model_data[c("price", "sin_hour", "cos_hour")])
z_time <- x_time

# Starting values
positive_log_y <- log(y[y > 0])

start_ssm <- c(
  qlogis(mean(y == 0)),
  mean(positive_log_y),
  log(sd(positive_log_y)),
  atanh(0.7),
  log(0.3)
)

# Fit models
mod_ssm <- nlm(
  nll_ssm, start_ssm,
  y = y, grid = grid
)

mod_price <- nlm(
  nll_ssm_cov, c(mod_ssm$estimate, 0),
  y = y, x = x_price, grid = grid
)

mod_change <- nlm(
  nll_ssm_cov, c(mod_ssm$estimate, 0),
  y = y, x = x_change, grid = grid
)

mod_time <- nlm(
  nll_ssm_cov, c(mod_ssm$estimate, 0, 0),
  y = y, x = x_time, grid = grid
)

mod_both <- nlm(
  nll_ssm_cov, c(mod_price$estimate, 0, 0),
  y = y, x = x_both, grid = grid
)

mod_zero_time <- nlm(
  nll_ssm_zero_cov, c(mod_both$estimate, 0, 0),
  y = y, x = x_both, z = z_time, grid = grid
)

# Model comparison
models <- list(
  "Gaussian SSM" = mod_ssm,
  "SSM + Yes price" = mod_price,
  "SSM + Price change" = mod_change,
  "SSM + Time of day" = mod_time,
  "SSM + Price + Time" = mod_both,
  "SSM + Price + Time + Zero Time" = mod_zero_time
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

cat("Hours:", length(y), "\n")
cat("Zero share:", round(mean(y == 0), 4), "\n")

print(comparison)

saveRDS(
  list(
    model = mod_zero_time,
    data = model_data,
    grid = grid
  ),
  "results/model_fit.rds"
)
