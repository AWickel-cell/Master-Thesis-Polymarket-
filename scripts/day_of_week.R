library(dplyr)
library(tidyr)
library(lubridate)
library(readr)
library(ggplot2)

source("functions/ssm_functions.R")

# Load existing model data
previous <- readRDS("results/model_fit.rds")

model_data <- previous$data |>
  mutate(
    weekday = wday(hour, week_start = 1),
    weekend = as.integer(weekday >= 6)
  )

y <- model_data$volume_usd
grid <- previous$grid

# Descriptive analysis
weekday_summary <- model_data |>
  group_by(weekday) |>
  summarise(
    intervals = n(),
    mean_volume = mean(volume_usd),
    median_volume = median(volume_usd),
    zero_share = mean(volume_usd == 0),
    .groups = "drop"
  )

weekend_summary <- model_data |>
  group_by(weekend) |>
  summarise(
    intervals = n(),
    mean_volume = mean(volume_usd),
    median_volume = median(volume_usd),
    zero_share = mean(volume_usd == 0),
    .groups = "drop"
  )

print(weekday_summary)
print(weekend_summary)

# Plot weekday patterns
p_volume <- ggplot(
  weekday_summary,
  aes(x = factor(weekday), y = mean_volume)
) +
  geom_col(fill = "steelblue") +
  labs(
    x = "Day of week (1 = Monday)",
    y = "Mean hourly volume (USD)"
  ) +
  theme_minimal()

p_zero <- ggplot(
  weekday_summary,
  aes(x = factor(weekday), y = zero_share)
) +
  geom_col(fill = "steelblue") +
  labs(
    x = "Day of week (1 = Monday)",
    y = "Share of zero-volume hours"
  ) +
  theme_minimal()

print(p_volume)
print(p_zero)

# Covariates
x_base <- as.matrix(
  model_data[c("price", "sin_hour", "cos_hour")]
)

z_base <- as.matrix(
  model_data[c("sin_hour", "cos_hour")]
)

x_weekend <- as.matrix(
  model_data[c("price", "sin_hour", "cos_hour", "weekend")]
)

z_weekend <- as.matrix(
  model_data[c("sin_hour", "cos_hour", "weekend")]
)

# Fit models
mod_base <- previous$model

mod_weekend_mean <- nlm(
  nll_ssm_zero_cov,
  c(mod_base$estimate[1:8], 0, mod_base$estimate[9:10]),
  y = y, x = x_weekend, z = z_base, grid = grid,
  hessian = TRUE
)

mod_weekend_zero <- nlm(
  nll_ssm_zero_cov,
  c(mod_base$estimate, 0),
  y = y, x = x_base, z = z_weekend, grid = grid,
  hessian = TRUE
)

mod_weekend_both <- nlm(
  nll_ssm_zero_cov,
  c(
    mod_weekend_mean$estimate[1:9],
    mod_weekend_mean$estimate[10:11],
    0
  ),
  y = y, x = x_weekend, z = z_weekend, grid = grid,
  hessian = TRUE
)

models <- list(
  "Baseline" = mod_base,
  "Weekend in mean" = mod_weekend_mean,
  "Weekend in zero probability" = mod_weekend_zero,
  "Weekend in both" = mod_weekend_both
)

# Model comparison
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

print(comparison)

# Confidence intervals
get_ci <- function(model, index, name) {
  beta <- model$estimate[index]
  hessian <- model$hessian
  
  se <- tryCatch({
    if (is.null(hessian) ||
        any(!is.finite(hessian)) ||
        min(eigen(hessian, symmetric = TRUE,
                  only.values = TRUE)$values) <= 0) {
      stop("Hessian is not positive definite")
    }
    sqrt(solve(hessian)[index, index])
  }, error = function(e) NA_real_)
  
  data.frame(
    model = name,
    estimate = beta,
    se = se,
    lower = beta - qnorm(0.975) * se,
    upper = beta + qnorm(0.975) * se
  )
}

confidence_intervals <- bind_rows(
  get_ci(mod_weekend_mean, 9, "Weekend in mean"),
  get_ci(mod_weekend_zero, 11, "Weekend in zero probability"),
  get_ci(mod_weekend_both, 9, "Weekend in both: mean"),
  get_ci(mod_weekend_both, 12, "Weekend in both: zero probability")
)

print(confidence_intervals)

# Save results
dir.create(
  "results/day_of_week",
  recursive = TRUE,
  showWarnings = FALSE
)

write_csv(
  weekday_summary,
  "results/day_of_week/weekday_summary.csv"
)

write_csv(
  weekend_summary,
  "results/day_of_week/weekend_summary.csv"
)

write_csv(
  comparison,
  "results/day_of_week/model_comparison.csv"
)

write_csv(
  confidence_intervals,
  "results/day_of_week/confidence_intervals.csv"
)

ggsave(
  "results/day_of_week/weekday_volume.png",
  p_volume, width = 8, height = 5
)

ggsave(
  "results/day_of_week/weekday_zero_share.png",
  p_zero, width = 8, height = 5
)

saveRDS(
  list(
    models = models,
    data = model_data,
    comparison = comparison,
    confidence_intervals = confidence_intervals
  ),
  "results/day_of_week/results.rds"
)
