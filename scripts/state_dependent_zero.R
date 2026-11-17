library(dplyr)
library(readr)

source("functions/ssm_functions.R")

# Load final univariate model
previous <- readRDS("results/day_of_week/results.rds")

model_data <- previous$data
mod_base <- previous$models[["Weekend in both"]]

y <- model_data$volume_usd

x <- as.matrix(
  model_data[c("price", "sin_hour", "cos_hour", "weekend")]
)

z <- as.matrix(
  model_data[c("sin_hour", "cos_hour", "weekend")]
)

grid <- seq(-5, 5, length.out = 151)

# Fit state-dependent zero model
mod_state_zero <- nlm(
  nll_ssm_state_zero,
  p = c(mod_base$estimate, 0),
  y = y,
  x = x,
  z = z,
  grid = grid,
  iterlim = 500,
  hessian = TRUE
)

models <- list(
  "State-independent zeros" = mod_base,
  "State-dependent zeros" = mod_state_zero
)

# Model comparison
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
) |>
  mutate(
    delta_aic = aic - min(aic),
    delta_bic = bic - min(bic)
  )

print(comparison)

# State-dependent zero effect
gamma_zero <- mod_state_zero$estimate[13]

cat("Gamma zero:", round(gamma_zero, 4), "\n")
cat("Phi:", round(tanh(mod_state_zero$estimate[4]), 4), "\n")
cat("Sigma_s:", round(exp(mod_state_zero$estimate[5]), 4), "\n")

# Confidence interval
hessian <- mod_state_zero$hessian

if (!is.null(hessian) &&
    all(is.finite(hessian)) &&
    min(eigen(hessian, symmetric = TRUE,
              only.values = TRUE)$values) > 0) {
  
  se <- sqrt(solve(hessian)[13, 13])
  
  cat("95% CI for gamma zero:",
      round(gamma_zero - 1.96 * se, 4),
      "to",
      round(gamma_zero + 1.96 * se, 4),
      "\n")
} else {
  cat("Confidence interval unavailable\n")
}

# Save results
dir.create(
  "results/state_dependent_zero",
  recursive = TRUE,
  showWarnings = FALSE
)

write_csv(
  comparison,
  "results/state_dependent_zero/comparison.csv"
)

saveRDS(
  list(
    models = models,
    data = model_data,
    grid = grid,
    comparison = comparison
  ),
  "results/state_dependent_zero/results.rds"
)
