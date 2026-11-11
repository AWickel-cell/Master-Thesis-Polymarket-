library(dplyr)
library(ggplot2)
library(lubridate)
library(readr)

# Prepare data
data <- read_csv("data/volume_hourly.csv", show_col_types = FALSE) |>
  arrange(hour) |>
  mutate(
    log_volume = if_else(volume_usd > 0, log(volume_usd), NA_real_),
    group = cumsum(is.na(log_volume))
  )

stopifnot(!anyNA(data$hour))

# Trading activity summary
print(
  data |>
    summarise(
      total_hours = n(),
      hours_with_trades = sum(n_trades > 0),
      zero_share = mean(n_trades == 0),
      total_trades = sum(n_trades),
      total_volume = sum(volume_usd)
    )
)

# Hourly trading volume
ggplot(data, aes(hour, log_volume, group = group)) +
  geom_line(linewidth = 0.4, na.rm = TRUE) +
  labs(
    title = "Log Hourly Trading Volume",
    x = "Date",
    y = "Log Trading Volume"
  ) +
  theme_minimal()

# One-week zoom
volume_week <- data |>
  filter(
    hour >= as_datetime("2026-10-01", tz = "UTC"),
    hour < as_datetime("2026-10-08", tz = "UTC")
  )

ggplot(volume_week, aes(hour, volume_usd)) +
  geom_point(size = 1.2) +
  labs(
    title = "Hourly Trading Volume: October 1–7, 2026",
    x = "Date",
    y = "Trading Volume (USD)"
  ) +
  theme_minimal()

# Distribution and autocorrelation
hist(
  na.omit(data$log_volume),
  breaks = 40,
  xlab = "Log Trading Volume",
  main = "Distribution of Positive Log Volumes"
)

acf(
  log1p(data$volume_usd),
  lag.max = 72,
  main = "ACF of Log(1 + Hourly Volume)"
)

# Largest volume peak
peak_hour <- data$hour[which.max(data$volume_usd)]

trades <- read_csv(
  "data/trades_raw_v2.csv",
  col_types = cols(
    proxy_wallet = col_character(),
    token_id = col_character(),
    condition_id = col_character(),
    transaction_hash = col_character(),
    .default = col_guess()
  ),
  show_col_types = FALSE
) |>
  mutate(
    datetime = as_datetime(timestamp, tz = "UTC"),
    hour = floor_date(datetime, "hour"),
    volume_usd = size * price
  )

cat("Peak hour:", format(peak_hour, tz = "UTC"), "\n")

print(
  trades |>
    filter(hour == peak_hour) |>
    summarise(
      n_trades = n(),
      total_volume = sum(volume_usd),
      largest_trade = max(volume_usd),
      share_largest = largest_trade / total_volume
    )
)

# April 2026 volume peaks
volume_peak <- data |>
  filter(
    hour >= as_datetime("2026-04-29", tz = "UTC"),
    hour < as_datetime("2026-05-02", tz = "UTC")
  )

ggplot(volume_peak, aes(hour, volume_usd)) +
  geom_line(linewidth = 0.4) +
  geom_point(size = 1) +
  scale_x_datetime(date_breaks = "1 day", date_labels = "%b %d") +
  labs(
    title = "Trading Volume Peaks: April 29 – May 1, 2026",
    x = "Date",
    y = "Trading Volume (USD)"
  ) +
  theme_minimal()

print(
  trades |>
    filter(
      hour >= as_datetime("2026-04-30 05:00:00", tz = "UTC"),
      hour < as_datetime("2026-04-30 10:00:00", tz = "UTC")
    ) |>
    group_by(hour) |>
    summarise(
      n_trades = n(),
      total_volume = sum(volume_usd),
      largest_trade = max(volume_usd),
      share_largest = largest_trade / total_volume,
      .groups = "drop"
    )
)
