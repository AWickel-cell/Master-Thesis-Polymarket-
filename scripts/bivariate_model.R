library(dplyr)
library(tidyr)
library(lubridate)
library(readr)
library(LaMa)

source("functions/ssm_functions.R")

# Settings
grid <- seq(-1.5, 1.5, length.out = 11)

# Load data
final <- readRDS("results/final_model/results.rds")

model_data <- final$data |>
  select(hour, volume_usd, sin_hour, cos_hour, weekend)

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
    hour = floor_date(as_datetime(timestamp, tz = "UTC"), "hour"),
    volume_usd = size * price
  )

hourly_outcome <- trades |>
  filter(outcome %in% c("Yes", "No")) |>
  group_by(hour, outcome) |>
  summarise(
    volume_usd = sum(volume_usd),
    .groups = "drop"
  ) |>
  pivot_wider(
    names_from = outcome,
    values_from = volume_usd,
    values_fill = 0
  )

model_data <- model_data |>
  left_join(hourly_outcome, by = "hour") |>
  mutate(
    Yes = replace_na(Yes, 0),
    No = replace_na(No, 0)
  ) |>
  arrange(hour)

stopifnot(
  nrow(model_data) == nrow(final$data),
  !anyDuplicated(model_data$hour),
  all(abs(model_data$Yes + model_data$No -
            model_data$volume_usd) < 1e-6)
)

y_yes <- model_data$Yes
y_no <- model_data$No

x <- as.matrix(
  model_data[c("sin_hour", "cos_hour", "weekend")]
)
z <- x

n <- length(y_yes)

cat("Observations:", n, "\n")
cat("Zero share Yes:", mean(y_yes == 0), "\n")
cat("Zero share No:", mean(y_no == 0), "\n")

# Transition matrix for a single latent state
transition_sep <- function(phi, sigma_s, grid) {
  logw <- -outer(
    grid, grid,
    function(a, b) (b - phi * a)^2
  ) / (2 * sigma_s^2)
  
  logw <- logw - apply(logw, 1, max)
  w <- exp(logw)
  
  w / rowSums(w)
}

# Univariate state-dependent zero SSM
nll_uni <- function(par, y) {
  if (!all(is.finite(par))) return(1e12)
  
  sigma <- exp(par[3])
  phi <- tanh(par[4])
  sigma_s <- exp(par[5])
  
  if (!all(is.finite(c(sigma, sigma_s))) ||
      sigma < 0.01 || sigma_s < 0.01 ||
      sigma_s > 5 || abs(phi) > 0.999) {
    return(1e12)
  }
  
  result <- tryCatch({
    Gamma <- transition_sep(phi, sigma_s, grid)
    delta <- LaMa::stationary(Gamma)
    
    eta <- as.vector(par[2] + x %*% par[6:8])
    logit_pi <- as.vector(par[1] + z %*% par[9:11])
    
    allprobs <- matrix(0, length(y), length(grid))
    positive <- y > 0
    
    for (j in seq_along(grid)) {
      pi_zero <- plogis(logit_pi + par[12] * grid[j])
      
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
  
  if (!is.finite(result)) return(1e12)
  result
}

# Fit separate SSMs
fit_uni <- function(y) {
  positive <- log(y[y > 0])
  
  start <- c(
    qlogis(mean(y == 0)),
    mean(positive),
    log(sd(positive)),
    atanh(0.9),
    log(0.1),
    rep(0, 6),
    -1
  )
  
  nlm(
    nll_uni,
    p = start,
    y = y,
    iterlim = 500
  )
}

cat("\nFitting Yes SSM...\n")
fit_yes <- fit_uni(y_yes)

cat("\nFitting No SSM...\n")
fit_no <- fit_uni(y_no)

# Joint SSM with one shared latent state
nll_joint <- function(par) {
  if (!all(is.finite(par))) return(1e12)
  
  sigma_y <- exp(par[5])
  sigma_n <- exp(par[6])
  phi <- tanh(par[7])
  sigma_s <- exp(par[8])
  lambda <- par[9]
  
  if (!all(is.finite(c(sigma_y, sigma_n, sigma_s))) ||
      min(sigma_y, sigma_n, sigma_s) < 0.01 ||
      sigma_s > 5 || abs(phi) > 0.999) {
    return(1e12)
  }
  
  result <- tryCatch({
    Gamma <- transition_sep(phi, sigma_s, grid)
    delta <- LaMa::stationary(Gamma)
    
    eta_y <- as.vector(par[3] + x %*% par[10:12])
    eta_n <- as.vector(par[4] + x %*% par[13:15])
    
    logit_y <- as.vector(par[1] + z %*% par[16:18])
    logit_n <- as.vector(par[2] + z %*% par[19:21])
    
    allprobs <- matrix(0, n, length(grid))
    
    pos_y <- y_yes > 0
    pos_n <- y_no > 0
    
    for (j in seq_along(grid)) {
      s <- grid[j]
      
      pi_y <- plogis(logit_y + par[22] * s)
      pi_n <- plogis(logit_n + par[23] * s)
      
      dens_y <- pi_y
      dens_n <- pi_n
      
      dens_y[pos_y] <- (1 - pi_y[pos_y]) *
        dlnorm(
          y_yes[pos_y],
          meanlog = eta_y[pos_y] + s,
          sdlog = sigma_y
        )
      
      dens_n[pos_n] <- (1 - pi_n[pos_n]) *
        dlnorm(
          y_no[pos_n],
          meanlog = eta_n[pos_n] + lambda * s,
          sdlog = sigma_n
        )
      
      allprobs[, j] <- dens_y * dens_n
    }
    
    -LaMa::forward(
      delta = delta,
      Gamma = Gamma,
      allprobs = allprobs
    )
  }, error = function(e) 1e12)
  
  if (!is.finite(result)) return(1e12)
  result
}

start_joint <- c(
  fit_yes$estimate[1],
  fit_no$estimate[1],
  fit_yes$estimate[2],
  fit_no$estimate[2],
  fit_yes$estimate[3],
  fit_no$estimate[3],
  fit_yes$estimate[4],
  fit_yes$estimate[5],
  1,
  fit_yes$estimate[6:8],
  fit_no$estimate[6:8],
  fit_yes$estimate[9:11],
  fit_no$estimate[9:11],
  fit_yes$estimate[12],
  fit_no$estimate[12]
)

cat("\nFitting joint SSM...\n")

fit_joint <- nlm(
  nll_joint,
  p = start_joint,
  iterlim = 500
)

# Two-dimensional state grid
states <- expand.grid(
  s_yes = grid,
  s_no = grid
)

sy <- states$s_yes
sn <- states$s_no

# Correlated AR(1) transition
transition_cor <- function(phi_y, phi_n, sd_y, sd_n, rho) {
  dy <- outer(sy, sy, function(a, b) b - phi_y * a) / sd_y
  dn <- outer(sn, sn, function(a, b) b - phi_n * a) / sd_n
  
  logw <- -(dy^2 - 2 * rho * dy * dn + dn^2) /
    (2 * (1 - rho^2))
  
  logw <- logw - apply(logw, 1, max)
  w <- exp(logw)
  
  w / rowSums(w)
}

# Stationary initial distribution
initial_cor <- function(phi_y, phi_n, sd_y, sd_n, rho) {
  var_y <- sd_y^2 / (1 - phi_y^2)
  var_n <- sd_n^2 / (1 - phi_n^2)
  
  cov_yn <- rho * sd_y * sd_n /
    (1 - phi_y * phi_n)
  
  corr <- cov_yn / sqrt(var_y * var_n)
  
  zy <- sy / sqrt(var_y)
  zn <- sn / sqrt(var_n)
  
  logw <- -(zy^2 - 2 * corr * zy * zn + zn^2) /
    (2 * (1 - corr^2))
  
  w <- exp(logw - max(logw))
  w / sum(w)
}

# Correlated bivariate SSM
nll_cor <- function(par) {
  if (!all(is.finite(par))) return(1e12)
  
  sd_y <- exp(par[5])
  sd_n <- exp(par[6])
  
  phi_y <- tanh(par[7])
  phi_n <- tanh(par[8])
  
  sd_sy <- exp(par[9])
  sd_sn <- exp(par[10])
  rho <- tanh(par[11])
  
  if (
    any(!is.finite(c(sd_y, sd_n, sd_sy, sd_sn))) ||
    min(sd_y, sd_n, sd_sy, sd_sn) < 0.01 ||
    max(sd_sy, sd_sn) > 5 ||
    max(abs(phi_y), abs(phi_n)) > 0.999 ||
    abs(rho) > 0.98
  ) return(1e12)
  
  result <- tryCatch({
    Gamma <- transition_cor(
      phi_y, phi_n, sd_sy, sd_sn, rho
    )
    
    delta <- LaMa::stationary(Gamma)
    
    eta_y <- as.vector(par[3] + x %*% par[12:14])
    eta_n <- as.vector(par[4] + x %*% par[15:17])
    
    logit_y <- as.vector(par[1] + z %*% par[18:20])
    logit_n <- as.vector(par[2] + z %*% par[21:23])
    
    pi_y <- plogis(outer(logit_y, par[24] * sy, "+"))
    pi_n <- plogis(outer(logit_n, par[25] * sn, "+"))
    
    dens_y <- pi_y
    dens_n <- pi_n
    
    pos_y <- y_yes > 0
    pos_n <- y_no > 0
    
    dens_y[pos_y, ] <-
      (1 - pi_y[pos_y, ]) *
      dlnorm(
        y_yes[pos_y],
        meanlog = outer(eta_y[pos_y], sy, "+"),
        sdlog = sd_y
      )
    
    dens_n[pos_n, ] <-
      (1 - pi_n[pos_n, ]) *
      dlnorm(
        y_no[pos_n],
        meanlog = outer(eta_n[pos_n], sn, "+"),
        sdlog = sd_n
      )
    
    allprobs <- dens_y * dens_n
    
    -LaMa::forward(
      delta = delta,
      Gamma = Gamma,
      allprobs = allprobs
    )
  }, error = function(e) 1e12)
  
  if (!is.finite(result)) return(1e12)
  result
}

# Starting values from separate SSMs
start_cor <- c(
  fit_yes$estimate[1],
  fit_no$estimate[1],
  fit_yes$estimate[2],
  fit_no$estimate[2],
  fit_yes$estimate[3],
  fit_no$estimate[3],
  fit_yes$estimate[4],
  fit_no$estimate[4],
  fit_yes$estimate[5],
  fit_no$estimate[5],
  atanh(0.2),
  fit_yes$estimate[6:8],
  fit_no$estimate[6:8],
  fit_yes$estimate[9:11],
  fit_no$estimate[9:11],
  fit_yes$estimate[12],
  fit_no$estimate[12]
)

cat("\nFitting correlated SSM...\n")

fit_cor <- nlm(
  nll_cor,
  p = start_cor,
  iterlim = 500
)

# Model comparison
models <- list(
  "Separate SSMs" = list(
    loglik = -fit_yes$minimum - fit_no$minimum,
    parameters = 24,
    convergence = paste(fit_yes$code, fit_no$code, sep = "/")
  ),
  "Joint SSM" = list(
    loglik = -fit_joint$minimum,
    parameters = 23,
    convergence = as.character(fit_joint$code)
  ),
  "Correlated SSM" = list(
    loglik = -fit_cor$minimum,
    parameters = 25,
    convergence = as.character(fit_cor$code)
  )
)

comparison <- data.frame(
  model = names(models),
  loglik = sapply(models, function(m) m$loglik),
  parameters = sapply(models, function(m) m$parameters),
  convergence = sapply(models, function(m) m$convergence),
  row.names = NULL
) |>
  mutate(
    aic = -2 * loglik + 2 * parameters,
    bic = -2 * loglik + log(n) * parameters,
    delta_aic = aic - min(aic),
    delta_bic = bic - min(bic)
  )

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

# Joint model parameters
joint_parameters <- data.frame(
  parameter = c(
    "phi", "sigma_s", "lambda",
    "gamma_yes", "gamma_no"
  ),
  estimate = c(
    tanh(fit_joint$estimate[7]),
    exp(fit_joint$estimate[8]),
    fit_joint$estimate[9],
    fit_joint$estimate[22:23]
  )
)

correlated_parameters <- data.frame(
  parameter = c(
    "phi_yes", "phi_no",
    "sigma_s_yes", "sigma_s_no", "rho",
    "gamma_yes", "gamma_no"
  ),
  estimate = c(
    tanh(fit_cor$estimate[7:8]),
    exp(fit_cor$estimate[9:10]),
    tanh(fit_cor$estimate[11]),
    fit_cor$estimate[24:25]
  )
)

cat("\nModel comparison:\n")
print(comparison)

cat("\nZero-volume dependence:\n")
print(zero_summary)

cat("\nJoint model parameters:\n")
print(joint_parameters)

cat("\nCorrelated model parameters:\n")
print(correlated_parameters)

# Save results
dir.create(
  "results/bivariate",
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  list(
    data = model_data,
    grid = grid,
    comparison = comparison,
    zero_summary = zero_summary,
    joint_parameters = joint_parameters,
    correlated_parameters = correlated_parameters,
    fits = list(
      yes = fit_yes,
      no = fit_no,
      joint = fit_joint,
      correlated = fit_cor
    )
  ),
  "results/bivariate/results.rds"
)





