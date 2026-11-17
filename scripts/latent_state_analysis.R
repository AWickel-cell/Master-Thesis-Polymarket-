library(ggplot2)

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

# Save plot
dir.create("results/latent_state",
           recursive = TRUE, showWarnings = FALSE)

ggsave(
  "results/latent_state/latent_events_full.png",
  p_events,
  width = 13,
  height = 5.5,
  dpi = 300
)
