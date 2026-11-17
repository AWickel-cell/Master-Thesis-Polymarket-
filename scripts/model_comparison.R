library(dplyr)
library(readr)

source("functions/ssm_functions.R")

# Load existing SSM
previous <- readRDS("results/day_of_week/results.rds")

model_data <- previous$data
mod_ssm <- previous$models[["Weekend in both"]]

y <- model_data$volume_usd

x <- as.matrix(
  model_data[c("price", "sin_hour", "cos_hour", "weekend")]
)

z <- as.matrix(
  model_data[c("sin_hour", "cos_hour", "weekend")]
)

stopifnot(
  length(y) == nrow(x),
  length(y) == nrow(z),
  !anyNA(y),
  !anyNA(x),
  !anyNA(z),
  all(y >= 0),
  length(mod_ssm$estimate) == 12
)

# Model without latent process
nll_no_latent <- function(par, y, x, z) {
  if (!all(is.finite(par))) return(1e12)
  
  mu <- par[2]
  sigma <- exp(par[3])
  
  if (!is.finite(sigma) || sigma < 0.01) {
    return(1e12)
  }
  
  beta <- par[4:7]
  alpha <- par[8:10]
  
  eta <- mu + as.vector(x %*% beta)
  pi_zero <- plogis(par[1] + as.vector(z %*% alpha))
  
  positive <- y > 0
  
  log_dens <- numeric(length(y))
  
  log_dens[!positive] <- log(pi_zero[!positive])
  
  log_dens[positive] <-
    log1p(-pi_zero[positive]) +
    dnorm(
      log(y[positive]),
      mean = eta[positive],
      sd = sigma,
      log = TRUE
    ) -
    log(y[positive])
  
  nll <- -sum(log_dens)
  
  if (!is.finite(nll)) return(1e12)
  
  nll
}

# Starting values from existing SSM
start_par <- c(
  mod_ssm$estimate[1:3],
  mod_ssm$estimate[6:12]
)

# Fit model without latent process
mod_no_latent <- nlm(
  nll_no_latent,
  p = start_par,
  y = y,
  x = x,
  z = z,
  iterlim = 500
)

# Model comparison
models <- list(
  "Without latent process" = mod_no_latent,
  "With latent AR(1)" = mod_ssm
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
)

comparison <- comparison |>
  mutate(
    delta_aic = aic - min(aic),
    delta_bic = bic - min(bic)
  )

cat("\nObservations:", length(y), "\n")
cat("Zero share:", round(mean(y == 0), 4), "\n")

cat("\nModel comparison:\n")
print(comparison)

# Parameter estimates
parameters <- data.frame(
  parameter = c(
    "logit_pi", "mu", "log_sigma",
    "price", "sin_hour", "cos_hour", "weekend_mean",
    "sin_zero", "cos_zero", "weekend_zero"
  ),
  no_latent = mod_no_latent$estimate,
  ssm = mod_ssm$estimate[c(1:3, 6:12)]
)

cat("\nParameter estimates:\n")
print(parameters)

cat("\nLatent process parameters:\n")
cat("phi:", round(tanh(mod_ssm$estimate[4]), 4), "\n")
cat("sigma_s:", round(exp(mod_ssm$estimate[5]), 4), "\n")

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
  parameters,
  "results/model_comparison/parameters.csv"
)

saveRDS(
  list(
    models = models,
    comparison = comparison,
    parameters = parameters
  ),
  "results/model_comparison/results.rds"
)
