
library(dplyr)
library(lubridate)
library(ggplot2)
library(LaMa)
library(readr)

source("functions/ssm_functions.R")

d_values <- c(1, 6, 12, 24)
grid <- seq(-3, 3, length.out = 95)

data <- read_csv("data/volume_hourly.csv", show_col_types = FALSE) |>
  arrange(hour)

stopifnot(!anyNA(data$hour), !anyDuplicated(data$hour))

dir.create("results/resolution_comparison",
           recursive = TRUE, showWarnings = FALSE)

# Baseline likelihoods
nll_baseline <- function(par, y) {
  pi_zero <- plogis(par[1])
  mu <- par[2]
  sigma <- exp(par[3])
  nu <- 2 + exp(par[4])
  
  zero <- y == 0
  z <- log(y[!zero])
  
  -(sum(zero) * log(pi_zero) +
      sum(log1p(-pi_zero) +
            dt((z - mu) / sigma, df = nu, log = TRUE) -
            log(sigma) - z))
}

nll_gaussian <- function(par, y) {
  pi_zero <- plogis(par[1])
  mu <- par[2]
  sigma <- exp(par[3])
  
  zero <- y == 0
  z <- log(y[!zero])
  
  -(sum(zero) * log(pi_zero) +
      sum(log1p(-pi_zero) +
            dnorm(z, mu, sigma, log = TRUE) - z))
}

# Aggregate volume
aggregate_volume <- function(d) {
  data |>
    mutate(period = floor_date(hour, paste(d, "hours"))) |>
    group_by(period) |>
    summarise(
      volume_usd = sum(volume_usd),
      n_hours = n(),
      .groups = "drop"
    ) |>
    filter(n_hours == d)
}

# Fit models
fit_models <- function(y) {
  log_y <- log(y[y > 0])
  p0 <- pmin(pmax(mean(y == 0), 1e-4), 1 - 1e-4)
  
  start_par <- c(
    qlogis(p0),
    mean(log_y),
    log(sd(log_y)),
    log(8)
  )
  
  mod_baseline <- nlm(nll_baseline, start_par, y = y)
  mod_gaussian <- nlm(nll_gaussian, start_par[1:3], y = y)
  
  mod_ssm <- nlm(
    nll_ssm,
    c(start_par[1:3], atanh(0.7), log(0.3)),
    y = y, grid = grid
  )
  
  models <- list(
    "Log-Student-t" = mod_baseline,
    "Log-Gaussian" = mod_gaussian,
    "Gaussian SSM" = mod_ssm
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
  
  list(models = models, comparison = comparison)
}

# One-step-ahead randomized PIT residuals
get_pit_residuals <- function(y, par) {
  pi_zero <- plogis(par[1])
  mu <- par[2]
  sigma <- exp(par[3])
  
  gamma <- ar1_transition(tanh(par[4]), exp(par[5]), grid)
  predicted <- as.vector(LaMa::stationary(gamma))
  
  pit <- numeric(length(y))
  
  for (t in seq_along(y)) {
    if (y[t] == 0) {
      pit[t] <- runif(1, 0, pi_zero)
      likelihood <- rep(pi_zero, length(grid))
    } else {
      z <- log(y[t])
      
      pit[t] <- pi_zero + (1 - pi_zero) *
        sum(predicted * pnorm(z, mu + grid, sigma))
      
      likelihood <- (1 - pi_zero) *
        dnorm(z, mu + grid, sigma) / y[t]
    }
    
    filtered <- predicted * likelihood
    filtered <- filtered / sum(filtered)
    
    predicted <- as.vector(filtered %*% gamma)
  }
  
  qnorm(pmin(pmax(pit, 1e-10), 1 - 1e-10))
}

# Compare resolutions
run_resolution <- function(d) {
  cat("\nResolution:", d, "hours\n")
  
  volume_data <- aggregate_volume(d)
  y <- volume_data$volume_usd
  
  fit <- fit_models(y)
  
  set.seed(123)
  residuals <- get_pit_residuals(
    y, fit$models[["Gaussian SSM"]]$estimate
  )
  
  folder <- file.path("results/resolution_comparison", paste0(d, "h"))
  dir.create(folder, showWarnings = FALSE)
  
  plot_data <- volume_data |>
    mutate(
      log_volume = if_else(volume_usd > 0, log(volume_usd), NA_real_)
    )
  
  p_volume <- ggplot(plot_data, aes(period, log_volume)) +
    geom_line(linewidth = 0.4, na.rm = TRUE) +
    labs(
      title = paste("Log Trading Volume -", d, "h"),
      x = "Date",
      y = "Log Trading Volume"
    ) +
    theme_minimal()
  
  p_qq <- ggplot(data.frame(residual = residuals),
                 aes(sample = residual)) +
    stat_qq(alpha = 0.5) +
    stat_qq_line() +
    labs(
      title = paste("PIT Q-Q Plot -", d, "h"),
      x = "Theoretical Quantiles",
      y = "Sample Quantiles"
    ) +
    theme_minimal()
  
  ggsave(file.path(folder, "volume.png"),
         p_volume, width = 9, height = 5)
  
  ggsave(file.path(folder, "pit_qq.png"),
         p_qq, width = 7, height = 5)
  
  png(file.path(folder, "pit_acf.png"),
      width = 1000, height = 650)
  
  acf(
    residuals,
    lag.max = min(72, length(y) - 1),
    main = paste("ACF of PIT Residuals -", d, "h")
  )
  
  dev.off()
  
  par_ssm <- fit$models[["Gaussian SSM"]]$estimate
  
  summary <- data.frame(
    d = d,
    n_obs = length(y),
    zero_share = mean(y == 0),
    mean_volume = mean(y),
    phi = tanh(par_ssm[4]),
    sigma_s = exp(par_ssm[5]),
    residual_mean = mean(residuals),
    residual_sd = sd(residuals),
    residual_acf1 = as.numeric(
      acf(residuals, lag.max = 1, plot = FALSE)$acf[2]
    )
  )
  
  comparison <- fit$comparison |>
    mutate(d = d, .before = 1)
  
  print(comparison)
  
  list(
    summary = summary,
    comparison = comparison,
    models = fit$models,
    residuals = residuals
  )
}

results <- lapply(d_values, run_resolution)
names(results) <- paste0(d_values, "h")

resolution_summary <- bind_rows(lapply(results, `[[`, "summary"))
model_comparison <- bind_rows(lapply(results, `[[`, "comparison"))

print(resolution_summary)
print(model_comparison)

write_csv(
  resolution_summary,
  "results/resolution_comparison/summary.csv"
)

write_csv(
  model_comparison,
  "results/resolution_comparison/models.csv"
)

saveRDS(
  results,
  "results/resolution_comparison/results.rds"
)
