library(dplyr)
library(tidyr)
library(lubridate)
library(LaMa)
library(readr)

source("functions/ssm_functions.R")

# Load data
final <- readRDS("results/final_model/results.rds")
model_data <- final$data
grid <- final$grid

trades <- read_csv("data/trades_raw_v2.csv", show_col_types = FALSE)

trades <- trades |>
  mutate(
    hour = floor_date(
      as.POSIXct(as.numeric(timestamp),
                 origin = "1970-01-01", tz = "UTC"),
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

model_data <- model_data |>
  select(hour, volume_usd, sin_hour, cos_hour, weekend) |>
  left_join(hourly_outcome, by = "hour") |>
  mutate(
    Yes = replace_na(Yes, 0),
    No = replace_na(No, 0)
  ) |>
  arrange(hour)

stopifnot(
  nrow(model_data) == 8854,
  !anyNA(model_data$hour),
  !anyDuplicated(model_data$hour),
  all(abs(model_data$Yes + model_data$No -
            model_data$volume_usd) < 1e-6)
)

y_yes <- model_data$Yes
y_no <- model_data$No

x <- as.matrix(model_data[c("sin_hour", "cos_hour", "weekend")])
z <- x

cat("Observations:", length(y_yes), "\n")
cat("Zero share Yes:", mean(y_yes == 0), "\n")
cat("Zero share No:", mean(y_no == 0), "\n")

# Separate state-dependent SSMs
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
    nll_ssm_state_zero,
    p = start,
    y = y,
    x = x,
    z = z,
    grid = grid,
    iterlim = 500
  )
}

cat("\nFitting Yes SSM...\n")
fit_yes <- fit_uni(y_yes)

cat("\nFitting No SSM...\n")
fit_no <- fit_uni(y_no)

# Joint SSM with shared latent state
nll_joint <- function(par, y_yes, y_no, x, z, grid) {
  if (!all(is.finite(par))) return(1e12)
  
  alpha_yes <- par[1]
  alpha_no <- par[2]
  mu_yes <- par[3]
  mu_no <- par[4]
  sigma_yes <- exp(par[5])
  sigma_no <- exp(par[6])
  
  phi <- tanh(par[7])
  sigma_s <- exp(par[8])
  lambda <- par[9]
  
  beta_yes <- par[10:12]
  beta_no <- par[13:15]
  
  zero_yes <- par[16:18]
  zero_no <- par[19:21]
  
  gamma_yes <- par[22]
  gamma_no <- par[23]
  
  if (!all(is.finite(c(sigma_yes, sigma_no, sigma_s))) ||
      min(sigma_yes, sigma_no, sigma_s) < 0.01 ||
      sigma_s > 5 || abs(phi) > 0.999) return(1e12)
  
  eta_yes <- as.vector(mu_yes + x %*% beta_yes)
  eta_no <- as.vector(mu_no + x %*% beta_no)
  
  logit_yes <- as.vector(alpha_yes + z %*% zero_yes)
  logit_no <- as.vector(alpha_no + z %*% zero_no)
  
  nll <- tryCatch({
    Gamma <- ar1_transition(phi, sigma_s, grid)
    delta <- LaMa::stationary(Gamma)
    
    allprobs <- matrix(0, length(y_yes), length(grid))
    
    pos_yes <- y_yes > 0
    pos_no <- y_no > 0
    
    for (j in seq_along(grid)) {
      s <- grid[j]
      
      pi_yes <- plogis(logit_yes + gamma_yes * s)
      pi_no <- plogis(logit_no + gamma_no * s)
      
      dens_yes <- pi_yes
      dens_no <- pi_no
      
      dens_yes[pos_yes] <- (1 - pi_yes[pos_yes]) *
        dlnorm(
          y_yes[pos_yes],
          meanlog = eta_yes[pos_yes] + s,
          sdlog = sigma_yes
        )
      
      dens_no[pos_no] <- (1 - pi_no[pos_no]) *
        dlnorm(
          y_no[pos_no],
          meanlog = eta_no[pos_no] + lambda * s,
          sdlog = sigma_no
        )
      
      allprobs[, j] <- dens_yes * dens_no
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

# Starting values
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
  y_yes = y_yes,
  y_no = y_no,
  x = x,
  z = z,
  grid = grid,
  iterlim = 500
)

# Model comparison
n <- length(y_yes)

comparison <- data.frame(
  model = c("Separate SSMs", "Joint SSM"),
  loglik = c(
    -fit_yes$minimum - fit_no$minimum,
    -fit_joint$minimum
  ),
  parameters = c(24, 23),
  convergence = c(
    paste(fit_yes$code, fit_no$code, sep = "/"),
    as.character(fit_joint$code)
  )
) |>
  mutate(
    aic = -2 * loglik + 2 * parameters,
    bic = -2 * loglik + log(n) * parameters,
    delta_aic = aic - min(aic),
    delta_bic = bic - min(bic)
  )

cat("\nModel comparison:\n")
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

cat("\nZero-volume dependence:\n")
print(zero_summary)

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

cat("\nJoint model parameters:\n")
print(joint_parameters)

# Save results
dir.create(
  "results/bivariate",
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  list(
    comparison = comparison,
    zero_summary = zero_summary,
    joint_parameters = joint_parameters,
    fits = list(
      yes = fit_yes,
      no = fit_no,
      joint = fit_joint
    )
  ),
  "results/bivariate/results.rds"
)
