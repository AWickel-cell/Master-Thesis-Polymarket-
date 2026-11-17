library(dplyr)
library(readr)
library(LaMa)

source("functions/ssm_functions.R")

# Load existing results
final <- readRDS("results/final_model/results.rds")
previous <- readRDS("results/model_comparison/results.rds")

model_data <- final$data
grid <- final$grid
y <- model_data$volume_usd

x <- as.matrix(model_data[c("sin_hour", "cos_hour", "weekend")])
z <- x

mod_baseline <- previous$models[["Baseline"]]
mod_regression <- previous$models[["Regression without latent process"]]
mod_final <- final$model

stopifnot(
  length(y) == nrow(x),
  !anyNA(y), !anyNA(x),
  all(y >= 0),
  length(mod_final$estimate) == 12
)

# Baseline SSM: latent AR(1), constant zero probability
nll_baseline_ssm <- function(par, y, grid) {
  if (!all(is.finite(par))) return(1e12)
  
  pi_zero <- plogis(par[1])
  mu <- par[2]
  sigma <- exp(par[3])
  phi <- tanh(par[4])
  sigma_s <- exp(par[5])
  
  if (sigma < 0.01 || sigma_s < 0.01 ||
      sigma_s > 5 || abs(phi) > 0.999) return(1e12)
  
  nll <- tryCatch({
    gamma <- ar1_transition(phi, sigma_s, grid)
    delta <- LaMa::stationary(gamma)
    
    allprobs <- matrix(0, length(y), length(grid))
    positive <- y > 0
    
    for (j in seq_along(grid)) {
      allprobs[!positive, j] <- pi_zero
      allprobs[positive, j] <- (1 - pi_zero) *
        dlnorm(y[positive],
               meanlog = mu + grid[j],
               sdlog = sigma)
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

# SSM with covariates and state-independent zeros
nll_cov_ssm <- function(par, y, x, z, grid) {
  if (!all(is.finite(par))) return(1e12)
  
  mu <- par[2]
  sigma <- exp(par[3])
  phi <- tanh(par[4])
  sigma_s <- exp(par[5])
  
  if (sigma < 0.01 || sigma_s < 0.01 ||
      sigma_s > 5 || abs(phi) > 0.999) return(1e12)
  
  beta <- par[6:8]
  alpha <- par[9:11]
  
  eta <- as.vector(mu + x %*% beta)
  pi_zero <- plogis(as.vector(par[1] + z %*% alpha))
  
  nll <- tryCatch({
    gamma <- ar1_transition(phi, sigma_s, grid)
    delta <- LaMa::stationary(gamma)
    
    allprobs <- matrix(0, length(y), length(grid))
    positive <- y > 0
    
    for (j in seq_along(grid)) {
      allprobs[!positive, j] <- pi_zero[!positive]
      allprobs[positive, j] <- (1 - pi_zero[positive]) *
        dlnorm(y[positive],
               meanlog = eta[positive] + grid[j],
               sdlog = sigma)
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

# Fit missing models
cat("\nFitting baseline SSM...\n")

start_baseline_ssm <- c(
  mod_baseline$estimate,
  mod_final$estimate[4:5]
)

mod_baseline_ssm <- nlm(
  nll_baseline_ssm,
  p = start_baseline_ssm,
  y = y,
  grid = grid,
  iterlim = 500
)

cat("\nFitting SSM with state-independent zeros...\n")

mod_cov_ssm <- nlm(
  nll_cov_ssm,
  p = mod_final$estimate[1:11],
  y = y,
  x = x,
  z = z,
  grid = grid,
  iterlim = 500
)

# Model comparison
models <- list(
  "Baseline" = mod_baseline,
  "Baseline SSM" = mod_baseline_ssm,
  "Regression without latent process" = mod_regression,
  "SSM with state-independent zeros" = mod_cov_ssm,
  "Final SSM with state-dependent zeros" = mod_final
)

comparison <- data.frame(
  model = names(models),
  loglik = -sapply(models, function(m) m$minimum),
  parameters = sapply(models, function(m) length(m$estimate)),
  aic = sapply(models, function(m) {
    2 * m$minimum + 2 * length(m$estimate)
  }),
  bic = sapply(models, function(m) {
    2 * m$minimum + log(length(y)) * length(m$estimate)
  }),
  convergence = sapply(models, function(m) m$code),
  row.names = NULL
) |>
  mutate(
    delta_aic = aic - min(aic),
    delta_bic = bic - min(bic)
  )

cat("\nObservations:", length(y), "\n")
cat("Zero share:", round(mean(y == 0), 4), "\n")
cat("\nModel comparison:\n")
print(comparison)

# Latent parameter comparison
latent_parameters <- data.frame(
  model = c(
    "Baseline SSM",
    "SSM with state-independent zeros",
    "Final SSM with state-dependent zeros"
  ),
  phi = sapply(
    list(mod_baseline_ssm, mod_cov_ssm, mod_final),
    function(m) tanh(m$estimate[4])
  ),
  sigma_s = sapply(
    list(mod_baseline_ssm, mod_cov_ssm, mod_final),
    function(m) exp(m$estimate[5])
  )
)

cat("\nLatent parameters:\n")
print(latent_parameters)

# Save results
dir.create(
  "results/model_comparison",
  recursive = TRUE,
  showWarnings = FALSE
)

write_csv(
  comparison,
  "results/model_comparison/comparison.csv"
)

write_csv(
  latent_parameters,
  "results/model_comparison/latent_parameters.csv"
)

saveRDS(
  list(
    models = models,
    comparison = comparison,
    latent_parameters = latent_parameters
  ),
  "results/model_comparison/results.rds"
)

