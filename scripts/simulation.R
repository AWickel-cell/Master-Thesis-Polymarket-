library(dplyr)
library(ggplot2)
library(readr)

source("functions/ssm_functions.R")

set.seed(123)

n_obs <- nrow(read_csv("data/volume_hourly.csv", show_col_types = FALSE))
n_sim <- 5
phi_values <- c(0.3, 0.7, 0.99)

mu <- 4.2975
sigma <- 2.1086
sigma_s <- 0.0953
pi_zero <- 0.217

grid <- seq(-3, 3, length.out = 95)

# Simulate Gaussian hurdle SSM
simulate_ssm <- function(n_obs, phi, sigma_s, mu, sigma, pi_zero) {
  state <- numeric(n_obs)
  state[1] <- rnorm(1, 0, sigma_s / sqrt(1 - phi^2))
  
  for (t in 2:n_obs) {
    state[t] <- phi * state[t - 1] + rnorm(1, 0, sigma_s)
  }
  
  y <- numeric(n_obs)
  positive <- runif(n_obs) > pi_zero
  
  y[positive] <- exp(
    rnorm(sum(positive), mu + state[positive], sigma)
  )
  
  y
}

# Fit simulated data
fit_simulation <- function(y, grid) {
  log_y <- log(y[y > 0])
  
  phi_starts <- c(0.3, 0.7, 0.95)
  sigma_s_starts <- c(0.1, 0.3)
  
  starts <- expand.grid(
    phi = phi_starts,
    sigma_s = sigma_s_starts
  )
  
  fits <- vector("list", nrow(starts))
  
  for (i in seq_len(nrow(starts))) {
    start <- c(
      qlogis(mean(y == 0)),
      mean(log_y),
      log(sd(log_y)),
      atanh(starts$phi[i]),
      log(starts$sigma_s[i])
    )
    
    nll_start <- nll_ssm(start, y, grid)
    
    fit <- tryCatch(
      nlm(
        nll_ssm,
        start,
        y = y,
        grid = grid,
        iterlim = 150,
        gradtol = 1e-4,
        steptol = 1e-5,
        hessian = FALSE
      ),
      error = function(e) NULL
    )
    
    if (is.null(fit)) next
    
    fits[[i]] <- data.frame(
      phi_hat = tanh(fit$estimate[4]),
      sigma_s_hat = exp(fit$estimate[5]),
      nll_start = nll_start,
      nll_fit = fit$minimum,
      convergence = fit$code,
      start_phi = starts$phi[i],
      start_sigma_s = starts$sigma_s[i]
    )
  }
  
  fits <- bind_rows(fits)
  
  if (nrow(fits) == 0) {
    return(data.frame(
      phi_hat = NA_real_,
      sigma_s_hat = NA_real_,
      nll_start = NA_real_,
      nll_fit = NA_real_,
      convergence = NA_integer_,
      start_phi = NA_real_,
      start_sigma_s = NA_real_
    ))
  }
  
  valid <- fits |>
    filter(
      convergence %in% c(1, 2),
      is.finite(nll_fit),
      nll_fit < 1e12
    )
  
  if (nrow(valid) == 0) {
    best <- fits[which.min(fits$nll_fit), ]
    best$phi_hat <- NA_real_
    best$sigma_s_hat <- NA_real_
    return(best)
  }
  
  valid[which.min(valid$nll_fit), ]
}

# Run simulations
results <- list()
k <- 1

for (phi_true in phi_values) {
  for (s in seq_len(n_sim)) {
    cat("Phi:", phi_true, "| Simulation:", s, "/", n_sim, "\n")
    
    y <- simulate_ssm(
      n_obs, phi_true, sigma_s, mu, sigma, pi_zero
    )
    
    fit <- fit_simulation(y, grid)
    
    results[[k]] <- data.frame(
      phi_true = phi_true,
      simulation = s,
      fit
    )
    
    k <- k + 1
  }
}

results <- bind_rows(results)

# Summarize estimates
summary_results <- results |>
  group_by(phi_true) |>
  summarise(
    mean_phi_hat = if (all(is.na(phi_hat))) NA_real_
    else mean(phi_hat, na.rm = TRUE),
    sd_phi_hat = if (sum(!is.na(phi_hat)) < 2) NA_real_
    else sd(phi_hat, na.rm = TRUE),
    mean_sigma_s_hat = if (all(is.na(sigma_s_hat))) NA_real_
    else mean(sigma_s_hat, na.rm = TRUE),
    n_converged = sum(
      convergence %in% c(1, 2) & !is.na(phi_hat),
      na.rm = TRUE
    ),
    .groups = "drop"
  )

print(summary_results)

# Check optimization
print(
  results |>
    select(
      phi_true, simulation,
      start_phi, start_sigma_s,
      phi_hat, sigma_s_hat,
      nll_start, nll_fit, convergence
    )
)

# Plot estimated persistence
p <- ggplot(results, aes(factor(phi_true), phi_hat)) +
  geom_boxplot(fill = "grey75", na.rm = TRUE) +
  geom_point(
    position = position_jitter(width = 0.08),
    alpha = 0.7,
    na.rm = TRUE
  ) +
  geom_point(
    data = data.frame(
      phi_true = phi_values,
      phi_hat = phi_values
    ),
    color = "red",
    size = 3
  ) +
  labs(
    x = "True phi",
    y = "Estimated phi"
  ) +
  theme_minimal()

print(p)

# Save results
dir.create(
  "results/simulation",
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  results,
  "results/simulation/results.rds"
)

saveRDS(
  summary_results,
  "results/simulation/summary.rds"
)

ggsave(
  "results/simulation/phi_estimates.png",
  plot = p,
  width = 7,
  height = 5,
  dpi = 300
)
