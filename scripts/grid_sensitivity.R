library(dplyr)
library(readr)

source("functions/ssm_functions.R")

# Load final model
final <- readRDS("results/final_model/results.rds")

model_data <- final$data
mod_ssm <- final$model

y <- model_data$volume_usd

x <- as.matrix(
  model_data[c("sin_hour", "cos_hour", "weekend")]
)

z <- as.matrix(
  model_data[c("sin_hour", "cos_hour", "weekend")]
)

stopifnot(
  length(mod_ssm$estimate) == 12,
  !anyNA(y),
  !anyNA(x),
  !anyNA(z)
)

# Grid configurations
grids <- list(
  "[-5,5] 151" = seq(-5, 5, length.out = 151),
  "[-4,4] 151" = seq(-4, 4, length.out = 151),
  "[-5,5] 201" = seq(-5, 5, length.out = 201)
)

# Fit models
fits <- list()
results <- list()

for (name in names(grids)) {
  
  cat("\nGrid:", name, "\n")
  
  grid <- grids[[name]]
  start_time <- Sys.time()
  
  if (name == "[-5,5] 151") {
    fit <- mod_ssm
  } else {
    fit <- nlm(
      nll_ssm_state_zero,
      p = mod_ssm$estimate,
      y = y,
      x = x,
      z = z,
      grid = grid,
      iterlim = 500
    )
  }
  
  runtime <- as.numeric(
    difftime(Sys.time(), start_time, units = "secs")
  )
  
  fits[[name]] <- fit
  
  results[[name]] <- data.frame(
    grid = name,
    points = length(grid),
    lower = min(grid),
    upper = max(grid),
    loglik = -fit$minimum,
    aic = 2 * fit$minimum + 2 * length(fit$estimate),
    bic = 2 * fit$minimum + log(length(y)) * length(fit$estimate),
    phi = tanh(fit$estimate[4]),
    sigma_s = exp(fit$estimate[5]),
    gamma_zero = fit$estimate[12],
    runtime_sec = runtime,
    convergence = fit$code
  )
}

# Results
comparison <- bind_rows(results) |>
  mutate(
    delta_loglik = loglik - loglik[1],
    delta_aic = aic - aic[1]
  )

print(comparison)

# Compare parameter estimates
parameters <- bind_rows(
  lapply(names(fits), function(name) {
    data.frame(
      grid = name,
      parameter = c(
        "logit_pi", "mu", "log_sigma",
        "atanh_phi", "log_sigma_s",
        "sin_hour", "cos_hour", "weekend_mean",
        "sin_zero", "cos_zero", "weekend_zero",
        "gamma_zero"
      ),
      estimate = fits[[name]]$estimate
    )
  })
)

print(parameters)

# Save results
dir.create(
  "results/grid_sensitivity",
  recursive = TRUE,
  showWarnings = FALSE
)

write_csv(
  comparison,
  "results/grid_sensitivity/comparison.csv"
)

write_csv(
  parameters,
  "results/grid_sensitivity/parameters.csv"
)

saveRDS(
  list(
    fits = fits,
    grids = grids,
    comparison = comparison,
    parameters = parameters
  ),
  "results/grid_sensitivity/results.rds"
)
