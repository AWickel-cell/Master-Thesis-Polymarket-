library(dplyr)
library(tidyr)
library(lubridate)
library(readr)

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
)

start_time <- as_datetime("2025-11-05 23:00:00", tz = "UTC")
end_time <- as_datetime("2026-11-10 00:00:00", tz = "UTC")

trades <- trades |>
  mutate(
    datetime = as_datetime(timestamp, tz = "UTC"),
    volume_usd = size * price,
    hour = floor_date(datetime, "hour")
  )

stopifnot(
  !anyNA(trades$datetime),
  !anyNA(trades$volume_usd),
  all(trades$volume_usd >= 0)
)

trades <- trades |>
  filter(datetime >= start_time, datetime < end_time)

volume_hourly <- trades |>
  group_by(hour) |>
  summarise(
    volume_usd = sum(volume_usd),
    n_trades = n(),
    .groups = "drop"
  ) |>
  complete(
    hour = seq(start_time, end_time - hours(1), by = "hour"),
    fill = list(volume_usd = 0, n_trades = 0)
  ) |>
  arrange(hour)

write_csv(volume_hourly, "data/volume_hourly.csv")

cat("Hours:", nrow(volume_hourly), "\n")
cat("Trades:", sum(volume_hourly$n_trades), "\n")
cat("Zero-volume hours:", sum(volume_hourly$volume_usd == 0), "\n")
cat("Zero-volume percentage:",
    round(mean(volume_hourly$volume_usd == 0) * 100, 2), "%\n")
cat("Total volume:", round(sum(volume_hourly$volume_usd), 2), "USD\n")
