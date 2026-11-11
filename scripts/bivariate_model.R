library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(LaMa)

source("functions/ssm_functions.R")

# Prepare hourly Yes/No volume
trades <- read.csv("data/trades_raw_v2.csv")
volume <- read.csv("data/volume_hourly.csv")

volume$hour <- as.POSIXct(
  volume$hour,
  format = "%Y-%m-%dT%H:%M:%SZ",
  tz = "UTC"
)

trades <- trades |>
  mutate(
    hour = floor_date(
      as.POSIXct(timestamp, origin = "1970-01-01", tz = "UTC"),
      "hour"
    ),
    volume_usd = size * price
  )

hourly_outcome <- trades |>
  filter(outcome %in% c("Yes", "No")) |>
  group_by(hour, outcome) |>
  summarise(volume_usd = sum(volume_usd), .groups = "drop") |>
  pivot_wider(
    names_from = outcome,
    values_from = volume_usd,
    values_fill = 0
  )

model_data <- volume |>
  select(hour, volume_usd) |>
  left_join(hourly_outcome, by = "hour") |>
  mutate(
    Yes = replace_na(Yes, 0),
    No = replace_na(No, 0)
  ) |>
  arrange(hour)

stopifnot(
  !anyNA(model_data$hour),
  !anyDuplicated(model_data$hour),
  all(abs(model_data$Yes + model_data$No -
            model_data$volume_usd) < 1e-6)
)

y_yes <- model_data$Yes
y_no <- model_data$No

grid <- seq(-3, 3, length.out = 95)

cat("Observations:", length(y_yes), "\n")
cat("Total volume:", sum(y_yes + y_no), "\n")


# Separate univariate SSMs
fit_uni <- function(y) {
  positive <- log(y[y > 0])
  
  start <- c(
    qlogis(mean(y == 0)),
    mean(positive),
    log(sd(positive)),
    atanh(0.7),
    log(0.3)
  )
  
  nlm(nll_ssm, start, y = y, grid = grid)
}

get_estimates <- function(fit) {
  par <- fit$estimate
  
  c(
    pi_zero = plogis(par[1]),
    mu = par[2],
    sigma = exp(par[3]),
    phi = tanh(par[4]),
    sigma_s = exp(par[5])
  )
}

fit_yes <- fit_uni(y_yes)
fit_no <- fit_uni(y_no)

parameters <- data.frame(
  parameter = names(get_estimates(fit_yes)),
  yes = unname(get_estimates(fit_yes)),
  no = unname(get_estimates(fit_no))
)

print(parameters)


# Joint SSM with shared latent state
nll_joint <- function(par, y_yes, y_no, grid) {
  if (!all(is.finite(par))) return(1e12)
  
  pi_yes <- plogis(par[1])
  pi_no <- plogis(par[2])
  mu_yes <- par[3]
  mu_no <- par[4]
  sigma_yes <- exp(par[5])
  sigma_no <- exp(par[6])
  phi <- tanh(par[7])
  sigma_s <- exp(par[8])
  lambda <- par[9]
  
  if (!all(is.finite(c(sigma_yes, sigma_no, sigma_s))) ||
      min(sigma_yes, sigma_no, sigma_s) < 0.01 ||
      sigma_s > 5 || abs(phi) > 0.999) return(1e12)
  
  nll <- tryCatch({
    gamma <- ar1_transition(phi, sigma_s, grid)
    delta <- LaMa::stationary(gamma)
    
    allprobs <- matrix(0, length(y_yes), length(grid))
    pos_yes <- y_yes > 0
    pos_no <- y_no > 0
    
    for (j in seq_along(grid)) {
      dens_yes <- rep(pi_yes, length(y_yes))
      dens_no <- rep(pi_no, length(y_no))
      
      dens_yes[pos_yes] <- (1 - pi_yes) *
        dlnorm(y_yes[pos_yes], mu_yes + grid[j], sigma_yes)
      
      dens_no[pos_no] <- (1 - pi_no) *
        dlnorm(
          y_no[pos_no],
          mu_no + lambda * grid[j],
          sigma_no
        )
      
      allprobs[, j] <- dens_yes * dens_no
    }
    
    -LaMa::forward(
      delta = delta,
      Gamma = gamma,
      allprobs = allprobs
    )
  }, error = function(e) 1e12)
  
  if (!is.finite(nll)) return(1e12)
  nll
}

start_joint <- c(
  qlogis(mean(y_yes == 0)),
  qlogis(mean(y_no == 0)),
  mean(log(y_yes[y_yes > 0])),
  mean(log(y_no[y_no > 0])),
  log(sd(log(y_yes[y_yes > 0]))),
  log(sd(log(y_no[y_no > 0]))),
  atanh(0.7),
  log(0.3),
  1
)

fit_joint <- nlm(
  nll_joint,
  start_joint,
  y_yes = y_yes,
  y_no = y_no,
  grid = grid
)

joint_parameters <- data.frame(
  parameter = c(
    "pi_yes", "pi_no", "mu_yes", "mu_no",
    "sigma_yes", "sigma_no", "phi", "sigma_s", "lambda"
  ),
  estimate = c(
    plogis(fit_joint$estimate[1:2]),
    fit_joint$estimate[3:4],
    exp(fit_joint$estimate[5:6]),
    tanh(fit_joint$estimate[7]),
    exp(fit_joint$estimate[8]),
    fit_joint$estimate[9]
  )
)

print(joint_parameters)


# Correlated SSM with two latent AR(1) processes
grid_cor <- seq(-3, 3, length.out = 11)

states <- expand.grid(
  s_yes = grid_cor,
  s_no = grid_cor
)

transition_cor <- function(
    phi_yes, phi_no, sigma_s_yes,
    sigma_s_no, rho, states
) {
  sy <- states$s_yes
  sn <- states$s_no
  m <- nrow(states)
  
  gamma <- matrix(0, m, m)
  
  for (i in seq_len(m)) {
    dy <- (sy - phi_yes * sy[i]) / sigma_s_yes
    dn <- (sn - phi_no * sn[i]) / sigma_s_no
    
    z <- (dy^2 - 2 * rho * dy * dn + dn^2) /
      (2 * (1 - rho^2))
    
    weights <- exp(-z)
    total <- sum(weights)
    
    if (!is.finite(total) || total <= 0)
      return(NULL)
    
    gamma[i, ] <- weights / total
  }
  
  gamma
}

# Approximate stationary distribution
initial_distribution_cor <- function(
    phi_yes, phi_no, sigma_s_yes,
    sigma_s_no, rho, states
) {
  var_yes <- sigma_s_yes^2 / (1 - phi_yes^2)
  var_no <- sigma_s_no^2 / (1 - phi_no^2)
  
  cov_yn <- rho * sigma_s_yes * sigma_s_no /
    (1 - phi_yes * phi_no)
  
  corr <- cov_yn / sqrt(var_yes * var_no)
  
  if (!is.finite(corr) || abs(corr) >= 1)
    return(NULL)
  
  zy <- states$s_yes / sqrt(var_yes)
  zn <- states$s_no / sqrt(var_no)
  
  z <- (zy^2 - 2 * corr * zy * zn + zn^2) /
    (2 * (1 - corr^2))
  
  weights <- exp(-z)
  
  if (!all(is.finite(weights)) || sum(weights) <= 0)
    return(NULL)
  
  weights / sum(weights)
}

# Correlated model likelihood
nll_cor <- function(par, y_yes, y_no, states) {
  if (!all(is.finite(par))) return(1e12)
  
  pi_yes <- plogis(par[1])
  pi_no <- plogis(par[2])
  mu_yes <- par[3]
  mu_no <- par[4]
  sigma_yes <- exp(par[5])
  sigma_no <- exp(par[6])
  phi_yes <- tanh(par[7])
  phi_no <- tanh(par[8])
  sigma_s_yes <- exp(par[9])
  sigma_s_no <- exp(par[10])
  rho <- tanh(par[11])
  
  if (!all(is.finite(c(
    sigma_yes, sigma_no, sigma_s_yes, sigma_s_no
  ))) ||
  min(sigma_yes, sigma_no, sigma_s_yes, sigma_s_no) < 0.01 ||
  max(sigma_s_yes, sigma_s_no) > 5 ||
  max(abs(phi_yes), abs(phi_no)) > 0.999 ||
  abs(rho) > 0.99) return(1e12)
  
  nll <- tryCatch({
    gamma <- transition_cor(
      phi_yes, phi_no,
      sigma_s_yes, sigma_s_no,
      rho, states
    )
    
    delta <- initial_distribution_cor(
      phi_yes, phi_no,
      sigma_s_yes, sigma_s_no,
      rho, states
    )
    
    if (is.null(gamma) || is.null(delta))
      return(1e12)
    
    allprobs <- matrix(0, length(y_yes), nrow(states))
    pos_yes <- y_yes > 0
    pos_no <- y_no > 0
    
    for (j in seq_len(nrow(states))) {
      dens_yes <- rep(pi_yes, length(y_yes))
      dens_no <- rep(pi_no, length(y_no))
      
      dens_yes[pos_yes] <- (1 - pi_yes) *
        dlnorm(
          y_yes[pos_yes],
          mu_yes + states$s_yes[j],
          sigma_yes
        )
      
      dens_no[pos_no] <- (1 - pi_no) *
        dlnorm(
          y_no[pos_no],
          mu_no + states$s_no[j],
          sigma_no
        )
      
      allprobs[, j] <- dens_yes * dens_no
    }
    
    -LaMa::forward(
      delta = delta,
      Gamma = gamma,
      allprobs = allprobs
    )
  }, error = function(e) 1e12)
  
  if (!is.finite(nll)) return(1e12)
  nll
}

# Starting values from separate SSMs
est_yes <- get_estimates(fit_yes)
est_no <- get_estimates(fit_no)

start_cor <- unname(c(
  qlogis(est_yes["pi_zero"]),
  qlogis(est_no["pi_zero"]),
  est_yes["mu"],
  est_no["mu"],
  log(est_yes["sigma"]),
  log(est_no["sigma"]),
  atanh(max(-0.98, min(est_yes["phi"], 0.98))),
  atanh(max(-0.98, min(est_no["phi"], 0.98))),
  log(est_yes["sigma_s"]),
  log(est_no["sigma_s"]),
  atanh(0.2)
))

initial_nll <- nll_cor(start_cor, y_yes, y_no, states)

cat("Initial correlated NLL:", initial_nll, "\n")

fit_cor <- if (is.finite(initial_nll) && initial_nll < 1e12) {
  nlm(
    nll_cor,
    start_cor,
    y_yes = y_yes,
    y_no = y_no,
    states = states,
    iterlim = 120
  )
} else {
  NULL
}

cor_parameters <- if (!is.null(fit_cor)) {
  data.frame(
    parameter = c(
      "pi_yes", "pi_no",
      "mu_yes", "mu_no",
      "sigma_yes", "sigma_no",
      "phi_yes", "phi_no",
      "sigma_s_yes", "sigma_s_no",
      "rho"
    ),
    estimate = c(
      plogis(fit_cor$estimate[1:2]),
      fit_cor$estimate[3:4],
      exp(fit_cor$estimate[5:6]),
      tanh(fit_cor$estimate[7:8]),
      exp(fit_cor$estimate[9:10]),
      tanh(fit_cor$estimate[11])
    )
  )
} else {
  NULL
}

print(cor_parameters)


# Model comparison
n <- length(y_yes)

loglik <- c(
  -fit_yes$minimum - fit_no$minimum,
  -fit_joint$minimum,
  if (is.null(fit_cor)) NA_real_ else -fit_cor$minimum
)

k <- c(10, 9, 11)

comparison <- data.frame(
  model = c(
    "Separate SSMs",
    "Joint SSM",
    "Correlated SSM"
  ),
  loglik = loglik,
  aic = -2 * loglik + 2 * k,
  bic = -2 * loglik + log(n) * k,
  convergence = c(
    paste(fit_yes$code, fit_no$code, sep = "/"),
    as.character(fit_joint$code),
    if (is.null(fit_cor)) NA_character_ else as.character(fit_cor$code)
  )
)

valid <- is.finite(comparison$loglik) &
  comparison$loglik > -1e11

comparison$delta_aic <- NA_real_
comparison$delta_bic <- NA_real_

comparison$delta_aic[valid] <-
  comparison$aic[valid] - min(comparison$aic[valid])

comparison$delta_bic[valid] <-
  comparison$bic[valid] - min(comparison$bic[valid])

print(comparison)


# Zero-volume dependence
zero_yes <- y_yes == 0
zero_no <- y_no == 0

zero_summary <- data.frame(
  observed_joint_zero = mean(zero_yes & zero_no),
  expected_independent = mean(zero_yes) * mean(zero_no),
  zero_correlation = cor(
    as.numeric(zero_yes),
    as.numeric(zero_no)
  )
)

print(zero_summary)


# Plot model comparison
plot_data <- comparison |>
  filter(is.finite(delta_aic), is.finite(delta_bic)) |>
  select(model, delta_aic, delta_bic) |>
  pivot_longer(
    cols = c(delta_aic, delta_bic),
    names_to = "criterion",
    values_to = "difference"
  )

p <- ggplot(plot_data, aes(model, difference)) +
  geom_col(fill = "grey50", width = 0.5) +
  facet_wrap(~ criterion, scales = "free_y") +
  labs(
    title = "State-Space Model Comparison",
    x = NULL,
    y = "Difference from Best Model"
  ) +
  theme_minimal()

print(p)


# Save results
dir.create(
  "results/bivariate",
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  list(
    comparison = comparison,
    separate_parameters = parameters,
    joint_parameters = joint_parameters,
    correlated_parameters = cor_parameters,
    zero_summary = zero_summary,
    fits = list(
      yes = fit_yes,
      no = fit_no,
      joint = fit_joint,
      correlated = fit_cor
    )
  ),
  "results/bivariate/results.rds"
)

ggsave(
  "results/bivariate/model_comparison.png",
  plot = p,
  width = 8,
  height = 5,
  dpi = 300
)


# Correlated SSM with dependent zeros

nll_cor_zero <- function(par, y_yes, y_no, states) {
  if (!all(is.finite(par))) return(1e12)
  
  pi_yes <- plogis(par[1])
  pi_no <- plogis(par[2])
  mu_yes <- par[3]
  mu_no <- par[4]
  sigma_yes <- exp(par[5])
  sigma_no <- exp(par[6])
  phi_yes <- tanh(par[7])
  phi_no <- tanh(par[8])
  sigma_s_yes <- exp(par[9])
  sigma_s_no <- exp(par[10])
  rho <- tanh(par[11])
  
  # Dependence parameter for joint zero probability
  lower <- max(0, pi_yes + pi_no - 1)
  upper <- min(pi_yes, pi_no)
  p00 <- lower + plogis(par[12]) * (upper - lower)
  
  p01 <- pi_yes - p00
  p10 <- pi_no - p00
  p11 <- 1 - pi_yes - pi_no + p00
  
  if (!all(is.finite(c(
    sigma_yes, sigma_no, sigma_s_yes, sigma_s_no
  ))) ||
  min(sigma_yes, sigma_no, sigma_s_yes, sigma_s_no) < 0.01 ||
  max(sigma_s_yes, sigma_s_no) > 5 ||
  max(abs(phi_yes), abs(phi_no)) > 0.999 ||
  abs(rho) > 0.99 ||
  min(p00, p01, p10, p11) <= 0) return(1e12)
  
  nll <- tryCatch({
    gamma <- transition_cor(
      phi_yes, phi_no,
      sigma_s_yes, sigma_s_no,
      rho, states
    )
    
    delta <- initial_distribution_cor(
      phi_yes, phi_no,
      sigma_s_yes, sigma_s_no,
      rho, states
    )
    
    if (is.null(gamma) || is.null(delta))
      return(1e12)
    
    n <- length(y_yes)
    m <- nrow(states)
    
    allprobs <- matrix(0, n, m)
    
    zero_yes <- y_yes == 0
    zero_no <- y_no == 0
    
    both_zero <- zero_yes & zero_no
    yes_zero <- zero_yes & !zero_no
    no_zero <- !zero_yes & zero_no
    both_positive <- !zero_yes & !zero_no
    
    for (j in seq_len(m)) {
      dens_yes <- rep(1, n)
      dens_no <- rep(1, n)
      
      pos_yes <- !zero_yes
      pos_no <- !zero_no
      
      dens_yes[pos_yes] <- dlnorm(
        y_yes[pos_yes],
        mu_yes + states$s_yes[j],
        sigma_yes
      )
      
      dens_no[pos_no] <- dlnorm(
        y_no[pos_no],
        mu_no + states$s_no[j],
        sigma_no
      )
      
      allprobs[both_zero, j] <- p00
      
      allprobs[yes_zero, j] <-
        p01 * dens_no[yes_zero]
      
      allprobs[no_zero, j] <-
        p10 * dens_yes[no_zero]
      
      allprobs[both_positive, j] <-
        p11 *
        dens_yes[both_positive] *
        dens_no[both_positive]
    }
    
    -LaMa::forward(
      delta = delta,
      Gamma = gamma,
      allprobs = allprobs
    )
  }, error = function(e) 1e12)
  
  if (!is.finite(nll)) return(1e12)
  nll
}


# Start from correlated SSM estimates
stopifnot(
  exists("fit_cor"),
  !is.null(fit_cor),
  length(fit_cor$estimate) == 11
)

start_zero <- fit_cor$estimate

pi_yes_start <- plogis(start_zero[1])
pi_no_start <- plogis(start_zero[2])

lower <- max(0, pi_yes_start + pi_no_start - 1)
upper <- min(pi_yes_start, pi_no_start)

p00_start <- pi_yes_start * pi_no_start

q_start <- (p00_start - lower) / (upper - lower)

start_zero <- c(
  start_zero,
  qlogis(q_start)
)

cat(
  "Initial NLL:",
  nll_cor_zero(start_zero, y_yes, y_no, states),
  "\n"
)


# Fit extended model
fit_cor_zero <- nlm(
  nll_cor_zero,
  start_zero,
  y_yes = y_yes,
  y_no = y_no,
  states = states,
  iterlim = 120
)

cat("Convergence:", fit_cor_zero$code, "\n")
cat("Final NLL:", fit_cor_zero$minimum, "\n")


# Extract estimates
par_zero <- fit_cor_zero$estimate

pi_yes_hat <- plogis(par_zero[1])
pi_no_hat <- plogis(par_zero[2])

lower <- max(0, pi_yes_hat + pi_no_hat - 1)
upper <- min(pi_yes_hat, pi_no_hat)

p00_hat <- lower +
  plogis(par_zero[12]) * (upper - lower)

zero_parameters <- data.frame(
  parameter = c(
    "pi_yes", "pi_no",
    "p00", "p01", "p10", "p11",
    "phi_yes", "phi_no",
    "sigma_s_yes", "sigma_s_no",
    "rho"
  ),
  estimate = c(
    pi_yes_hat,
    pi_no_hat,
    p00_hat,
    pi_yes_hat - p00_hat,
    pi_no_hat - p00_hat,
    1 - pi_yes_hat - pi_no_hat + p00_hat,
    tanh(par_zero[7]),
    tanh(par_zero[8]),
    exp(par_zero[9]),
    exp(par_zero[10]),
    tanh(par_zero[11])
  )
)

print(zero_parameters)


# Compare all four models
comparison_zero <- comparison |>
  select(model, loglik, aic, bic, convergence)

valid_zero <- is.finite(fit_cor_zero$minimum) &&
  fit_cor_zero$minimum < 1e12

comparison_zero <- bind_rows(
  comparison_zero,
  data.frame(
    model = "Correlated SSM + Dependent Zeros",
    loglik = if (valid_zero) -fit_cor_zero$minimum else NA_real_,
    aic = if (valid_zero) 2 * fit_cor_zero$minimum + 24 else NA_real_,
    bic = if (valid_zero) 2 * fit_cor_zero$minimum +
      log(length(y_yes)) * 12 else NA_real_,
    convergence = as.character(fit_cor_zero$code)
  )
)

valid <- is.finite(comparison_zero$aic) &
  is.finite(comparison_zero$bic)

comparison_zero$delta_aic <- NA_real_
comparison_zero$delta_bic <- NA_real_

comparison_zero$delta_aic[valid] <-
  comparison_zero$aic[valid] -
  min(comparison_zero$aic[valid])

comparison_zero$delta_bic[valid] <-
  comparison_zero$bic[valid] -
  min(comparison_zero$bic[valid])

print(comparison_zero)


# Save extended results
saveRDS(
  list(
    fit = fit_cor_zero,
    parameters = zero_parameters,
    comparison = comparison_zero
  ),
  "results/bivariate/dependent_zeros.rds"
)

