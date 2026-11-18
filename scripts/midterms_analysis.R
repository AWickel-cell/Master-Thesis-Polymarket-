library(dplyr)
library(lubridate)
library(readr)
library(tidyr)
library(LaMa)

source("functions/ssm_functions.R")

start_time <- as_datetime("2026-11-01 00:00:00", tz = "UTC")
end_time <- as_datetime("2026-11-10 00:00:00", tz = "UTC")

final <- readRDS("results/final_model/results.rds")
grid <- final$grid

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

stopifnot(
  !anyNA(trades$datetime),
  !anyNA(trades$volume_usd),
  all(trades$volume_usd >= 0),
  !anyNA(events$time_utc)
)

events <- events |>
  filter(time_utc >= start_time, time_utc < end_time)

# Aggregate volume and CNN calls
prepare_data <- function(minutes_value) {
  unit <- paste(minutes_value, "minutes")
  
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
        end_time - lubridate::minutes(minutes_value),
        by = paste(minutes_value, "min")
      ),
      fill = list(volume_usd = 0, n_trades = 0)
    ) |>
    arrange(period)
  
  calls <- events |>
    mutate(period = floor_date(time_utc, unit)) |>
    count(period, name = "n_calls")
  
  volume |>
    left_join(calls, by = "period") |>
    mutate(
      n_calls = replace_na(n_calls, 0L),
      call_dummy = as.integer(n_calls > 0),
      call_lag1 = lag(call_dummy, default = 0L),
      hour_decimal = hour(period) + minute(period) / 60,
      sin_hour = sin(2 * pi * hour_decimal / 24),
      cos_hour = cos(2 * pi * hour_decimal / 24),
      weekend = as.integer(wday(period, week_start = 1) >= 6)
    )
}

# SSM with optional CNN covariate
nll_midterms <- function(par, y, x, z, grid) {
  if (!all(is.finite(par))) return(1e12)
  
  mu <- par[2]
  sigma <- exp(par[3])
  phi <- tanh(par[4])
  sigma_s <- exp(par[5])
  
  if (!all(is.finite(c(sigma, sigma_s))) ||
      sigma < 0.01 || sigma_s < 0.01 ||
      sigma_s > 5 || abs(phi) > 0.999) {
    return(1e12)
  }
  
  n_beta <- ncol(x)
  n_alpha <- ncol(z)
  
  if (length(par) != 6 + n_beta + n_alpha) {
    return(1e12)
  }
  
  beta <- par[6:(5 + n_beta)]
  alpha <- par[(6 + n_beta):(5 + n_beta + n_alpha)]
  gamma_zero <- par[6 + n_beta + n_alpha]
  
  eta <- as.vector(mu + x %*% beta)
  logit_pi <- as.vector(par[1] + z %*% alpha)
  
  nll <- tryCatch({
    Gamma <- ar1_transition(phi, sigma_s, grid)
    delta <- LaMa::stationary(Gamma)
    
    allprobs <- matrix(0, length(y), length(grid))
    positive <- y > 0
    
    for (j in seq_along(grid)) {
      pi_zero <- plogis(logit_pi + gamma_zero * grid[j])
      
      allprobs[!positive, j] <- pi_zero[!positive]
      
      allprobs[positive, j] <-
        (1 - pi_zero[positive]) *
        dlnorm(
          y[positive],
          meanlog = eta[positive] + grid[j],
          sdlog = sigma
        )
    }
    
    -LaMa::forward(
      delta = delta,
      Gamma = Gamma,
      allprobs = allprobs
    )
  }, error = function(e) 1e12)
  
  if (!is.finite(nll)) return(1e12)
  nll
}

# Calculate 95% confidence interval for CNN coefficient
get_ci <- function(model) {
  beta <- model$estimate[9]
  hessian <- model$hessian
  
  se <- tryCatch({
    if (is.null(hessian) ||
        any(!is.finite(hessian)) ||
        min(eigen(
          hessian,
          symmetric = TRUE,
          only.values = TRUE
        )$values) <= 0) {
      NA_real_
    } else {
      sqrt(solve(hessian)[9, 9])
    }
  }, error = function(e) NA_real_)
  
  data.frame(
    estimate = beta,
    se = se,
    lower = beta - qnorm(0.975) * se,
    upper = beta + qnorm(0.975) * se
  )
}

# Fit models
fit_midterms <- function(model_data, resolution) {
  y <- model_data$volume_usd
  
  x_base <- as.matrix(
    model_data[c("sin_hour", "cos_hour", "weekend")]
  )
  z <- x_base
  
  x_count <- cbind(x_base, n_calls = model_data$n_calls)
  x_dummy <- cbind(x_base, call_dummy = model_data$call_dummy)
  x_lag <- cbind(x_base, call_lag1 = model_data$call_lag1)
  
  stopifnot(
    !anyNA(y),
    !anyNA(x_base),
    !anyNA(z),
    !anyNA(x_count),
    !anyNA(x_dummy),
    !anyNA(x_lag)
  )
  
  start_base <- final$model$estimate
  
  cat("\nFitting baseline:", resolution, "\n")
  
  mod_base <- nlm(
    nll_midterms,
    p = start_base,
    y = y,
    x = x_base,
    z = z,
    grid = grid,
    iterlim = 500
  )
  
  # Insert CNN coefficient before zero-model parameters
  start_event <- c(
    mod_base$estimate[1:8],
    0,
    mod_base$estimate[9:12]
  )
  
  cat("\nFitting CNN count:", resolution, "\n")
  
  mod_count <- nlm(
    nll_midterms,
    p = start_event,
    y = y,
    x = x_count,
    z = z,
    grid = grid,
    hessian = TRUE,
    iterlim = 500
  )
  
  cat("\nFitting CNN dummy:", resolution, "\n")
  
  mod_dummy <- nlm(
    nll_midterms,
    p = start_event,
    y = y,
    x = x_dummy,
    z = z,
    grid = grid,
    hessian = TRUE,
    iterlim = 500
  )
  
  cat("\nFitting CNN lag 1:", resolution, "\n")
  
  mod_lag <- nlm(
    nll_midterms,
    p = start_event,
    y = y,
    x = x_lag,
    z = z,
    grid = grid,
    hessian = TRUE,
    iterlim = 500
  )
  
  models <- list(
    "SSM" = mod_base,
    "SSM + CNN count" = mod_count,
    "SSM + CNN dummy" = mod_dummy,
    "SSM + CNN lag 1" = mod_lag
  )
  
  comparison <- data.frame(
    resolution = resolution,
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
  
  confidence_intervals <- bind_rows(
    mutate(get_ci(mod_count), model = "CNN count"),
    mutate(get_ci(mod_dummy), model = "CNN dummy"),
    mutate(get_ci(mod_lag), model = "CNN lag 1")
  ) |>
    mutate(resolution = resolution) |>
    select(resolution, model, estimate, se, lower, upper)
  
  cat("\nResolution:", resolution, "\n")
  cat("Intervals:", length(y), "\n")
  cat("Trades:", sum(model_data$n_trades), "\n")
  cat("Zero share:", round(mean(y == 0), 4), "\n")
  cat("CNN calls:", sum(model_data$n_calls), "\n")
  
  cat("\nModel comparison:\n")
  print(comparison)
  
  cat("\nCNN coefficients and 95% confidence intervals:\n")
  print(confidence_intervals)
  
  if (any(comparison$convergence != 1)) {
    warning(
      "Some models did not reach nlm convergence code 1: ",
      resolution
    )
  }
  
  list(
    data = model_data,
    models = models,
    comparison = comparison,
    confidence_intervals = confidence_intervals
  )
}

# Run both resolutions
results <- list()

for (minutes_value in c(60, 15)) {
  model_data <- prepare_data(minutes_value)
  
  expected <- as.integer(
    as.numeric(difftime(
      end_time, start_time, units = "mins"
    )) / minutes_value
  )
  
  stopifnot(
    nrow(model_data) == expected,
    !anyNA(model_data$period),
    !anyNA(model_data$n_calls),
    sum(model_data$n_calls) == nrow(events)
  )
  
  results[[paste0(minutes_value, "min")]] <-
    fit_midterms(
      model_data,
      paste0(minutes_value, "min")
    )
}

# Save results
comparison <- bind_rows(
  lapply(results, `[[`, "comparison")
)

confidence_intervals <- bind_rows(
  lapply(results, `[[`, "confidence_intervals")
)

dir.create(
  "results/midterms",
  recursive = TRUE,
  showWarnings = FALSE
)

write_csv(
  comparison,
  "results/midterms/model_comparison.csv"
)

write_csv(
  confidence_intervals,
  "results/midterms/confidence_intervals.csv"
)

saveRDS(
  results,
  "results/midterms/results.rds"
)

