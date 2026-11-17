library(dplyr)
library(ggplot2)
library(LaMa)

source("functions/ssm_functions.R")

# Load final model
final <- readRDS("results/final_model/results.rds")

mod <- final$model
model_data <- final$data
grid <- final$grid

y <- model_data$volume_usd
x <- as.matrix(model_data[c("sin_hour", "cos_hour", "weekend")])
z <- x

# Model parameters
par <- mod$estimate

mu <- par[2]
sigma <- exp(par[3])
phi <- tanh(par[4])
sigma_s <- exp(par[5])

beta <- par[6:8]
alpha <- par[9:11]
gamma_zero <- par[12]

eta <- as.vector(mu + x %*% beta)
logit_pi <- as.vector(par[1] + z %*% alpha)

# Transition and observation probabilities
Gamma <- ar1_transition(phi, sigma_s, grid)
delta <- LaMa::stationary(Gamma)

allprobs <- matrix(0, length(y), length(grid))
positive <- y > 0

for (j in seq_along(grid)) {
  pi_zero <- plogis(logit_pi + gamma_zero * grid[j])
  
  allprobs[!positive, j] <- pi_zero[!positive]
  
  allprobs[positive, j] <- (1 - pi_zero[positive]) *
    dlnorm(
      y[positive],
      meanlog = eta[positive] + grid[j],
      sdlog = sigma
    )
}

# Smoothed latent state
state_probs <- LaMa::stateprobs(
  delta = delta,
  Gamma = Gamma,
  allprobs = allprobs
)

latent_state <- as.vector(state_probs %*% grid)

plot_data <- data.frame(
  hour = model_data$hour,
  latent_state = latent_state
)

# Check
head(plot_data)
summary(plot_data$latent_state)

# Political events
events <- data.frame(
  date = as.POSIXct(c(
    "2026-01-03",
    "2026-02-28",
    "2026-04-06",
    "2026-09-15",
    "2026-11-03"
  ), tz = "UTC"),
  event = c(
    "Venezuela strike",
    "Iran strikes",
    "Impeachment resolution",
    "Impeachment vote",
    "Midterms"
  ),
  label_y = c(2.95, 2.35, -1.65, -1.65, 2.75),
  vjust = c(-0.5, -0.5, 1.2, 1.2, -0.5)
)

# Plot
p_events <- ggplot(plot_data, aes(hour, latent_state)) +
  geom_line(color = "steelblue", linewidth = 0.5) +
  geom_hline(
    yintercept = 0,
    color = "grey75",
    linewidth = 0.3
  ) +
  geom_vline(
    data = events,
    aes(xintercept = date),
    color = "grey65",
    linetype = "dashed",
    linewidth = 0.45
  ) +
  geom_text(
    data = events,
    aes(x = date, y = label_y, label = event, vjust = vjust),
    inherit.aes = FALSE,
    size = 3.2,
    color = "grey30"
  ) +
  scale_x_datetime(
    date_breaks = "2 months",
    date_labels = "%b %Y"
  ) +
  coord_cartesian(ylim = c(-2, 3.2)) +
  labs(
    title = "Estimated latent state with political events",
    x = "Time",
    y = "Smoothed latent state"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90"),
    plot.title = element_text(face = "bold"),
    plot.margin = margin(15, 20, 10, 10)
  )

print(p_events)

# Save results
dir.create(
  "results/latent_state",
  recursive = TRUE,
  showWarnings = FALSE
)

ggsave(
  "results/latent_state/latent_events_full.png",
  p_events,
  width = 13,
  height = 5.5,
  dpi = 300
)

saveRDS(
  list(
    plot_data = plot_data,
    events = events
  ),
  "results/latent_state/results.rds"
)
