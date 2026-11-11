library(dplyr)
library(lubridate)
library(readr)
library(tidyr)

source("functions/ssm_functions.R")

start_time <- as_datetime("2026-11-01 00:00:00", tz = "UTC")
end_time <- as_datetime("2026-11-10 00:00:00", tz = "UTC")
grid <- seq(-5, 5, length.out = 151)

# Load trades and events
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
) |>
  mutate(
    datetime = as_datetime(timestamp, tz = "UTC"),
    volume_usd = size * price
  ) |>
  filter(datetime >= start_time, datetime < end_time)

events <- read_csv(
  "data/senate_events.csv",
  show_col_types = FALSE
) |>
  mutate(time_utc = as_datetime(time_utc, tz = "UTC"))

# Aggregate volume and CNN calls
prepare_data <- function(minutes) {
  unit <- paste(minutes, "minutes")
  
  volume <- trades |>
    mutate(period = floor_date(datetime, unit)) |>
    group_by(period) |>
    summarise(
      volume_usd = sum(volume_usd),
      n_trades = n(),
      .groups = "drop"
    ) |>
    complete(
      period = seq(
        start_time,
        end_time - minutes(minutes),
        by = paste(minutes, "min")
      ),
      fill = list(volume_usd = 0, n_trades = 0)
    ) |>
    arrange(period)
  
  calls <- events |>
    mutate(period = floor_date(time_utc, unit)) |>
    count(period, name = "n_calls")
  
  volume |>
    left_join(calls, by = "period") |>
    mutate(n_calls = replace_na(n_calls, 0L))
}

# Fit SSMs
fit_midterms <- function(model_data, minutes) {
  y <- model_data$volume_usd
  x_count <- as.matrix(model_data["n_calls"])
  x_dummy <- as.matrix(data.frame(
    call_dummy = as.integer(model_data$n_calls > 0)
  ))
  
  log_y <- log(y[y > 0])
  
  start_par <- c(
    qlogis(mean(y == 0)),
    mean(log_y),
    log(sd(log_y)),
    atanh(0.7),
    log(0.3)
  )
  
  mod_ssm <- nlm(
    nll_ssm, start_par,
    y = y, grid = grid
  )
  
  mod_count <- nlm(
    nll_ssm_cov,
    c(mod_ssm$estimate, 0),
    y = y, x = x_count, grid = grid
  )
  
  mod_dummy <- nlm(
    nll_ssm_cov,
    c(mod_ssm$estimate, 0),
    y = y, x = x_dummy, grid = grid
  )
  
  models <- list(
    "Gaussian SSM" = mod_ssm,
    "SSM + CNN count" = mod_count,
    "SSM + CNN dummy" = mod_dummy
  )
  
  comparison <- data.frame(
    resolution = paste0(minutes, "min"),
    model = names(models),
    loglik = -sapply(models, function(m) m$minimum),
    aic = sapply(models, function(m)
      2 * m$minimum + 2 * length(m$estimate)),
    bic = sapply(models, function(m)
      2 * m$minimum + log(length(y)) * length(m$estimate)),
    convergence = sapply(models, function(m) m$code),
    row.names = NULL
  )
  
  cat("\nResolution:", minutes, "minutes\n")
  cat("Intervals:", length(y), "\n")
  cat("Trades:", sum(model_data$n_trades), "\n")
  cat("Zero share:", round(mean(y == 0), 4), "\n")
  cat("CNN calls:", sum(model_data$n_calls), "\n")
  cat("Intervals with calls:", sum(model_data$n_calls > 0), "\n")
  cat("Count coefficient:", round(mod_count$estimate[6], 4), "\n")
  cat("Dummy coefficient:", round(mod_dummy$estimate[6], 4), "\n")
  
  print(comparison)
  
  list(
    data = model_data,
    models = models,
    comparison = comparison
  )
}

# Run both resolutions
results <- list()

for (minutes_value in c(60, 15)) {
  model_data <- prepare_data(minutes_value)
  
  expected <- as.integer(
    as.numeric(difftime(end_time, start_time, units = "mins")) /
      minutes_value
  )
  
  stopifnot(
    nrow(model_data) == expected,
    !anyNA(model_data$period),
    !anyNA(model_data$n_calls),
    sum(model_data$n_calls) == nrow(events)
  )
  
  results[[paste0(minutes_value, "min")]] <-
    fit_midterms(model_data, minutes_value)
}

comparison <- bind_rows(lapply(results, `[[`, "comparison"))

dir.create("results/midterms", recursive = TRUE, showWarnings = FALSE)

write_csv(comparison, "results/midterms/model_comparison.csv")
saveRDS(results, "results/midterms/results.rds")
