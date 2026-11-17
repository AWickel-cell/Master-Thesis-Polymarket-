library(dplyr)
library(readr)

source("functions/ssm_functions.R")

# Load final model and data
previous <- readRDS("results/state_dependent_zero/results.rds")
price_comparison <- readRDS("results/price_comparison/results.rds")

model_data <- previous$data
grid <- previous$grid

mod_final <- price_comparison$models[["Without price"]]

y <- model_data$volume_usd

x <- as.matrix(
  model_data[c("sin_hour", "cos_hour", "weekend")]
)

z <- as.matrix(
  model_data[c("sin_hour", "cos_hour", "weekend")]
)

stopifnot(
  length(mod_final$estimate) == 12,
  length(y) == nrow(x),
  length(y) == nrow(z)
)

# Model summary
comparison <- data.frame(
  model = "State-dependent zero SSM without price",
  loglik = -mod_final$minimum,
  parameters = length(mod_final$estimate),
  aic = 2 * mod_final$minimum + 2 * length(mod_final$estimate),
  bic = 2 * mod_final$minimum +
    log(length(y)) * length(mod_final$estimate),
  convergence = mod_final$code
)

cat("Hours:", length(y), "\n")
cat("Zero share:", round(mean(y == 0), 4), "\n")

print(comparison)

cat("\nLatent process parameters:\n")
cat("phi:", round(tanh(mod_final$estimate[4]), 4), "\n")
cat("sigma_s:", round(exp(mod_final$estimate[5]), 4), "\n")
cat("gamma_zero:", round(mod_final$estimate[12], 4), "\n")

# Save final model
dir.create(
  "results/final_model",
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  list(
    model = mod_final,
    data = model_data,
    grid = grid,
    comparison = comparison
  ),
  "results/final_model/results.rds"
)
