
library(dplyr)
library(lubridate)
library(readr)

# Load trades
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
    volume_usd = size * price
  )

# Define analysis windows (UTC)
windows <- data.frame(
  period = c(
    "January", "March", "April",
    "June_July", "September", "November"
  ),
  start = as_datetime(c(
    "2025-12-30",
    "2026-02-27",
    "2026-04-04",
    "2026-06-25",
    "2026-09-12",
    "2026-11-01"
  ), tz = "UTC"),
  end = as_datetime(c(
    "2026-01-07",
    "2026-03-08",
    "2026-04-11",
    "2026-07-10",
    "2026-09-19",
    "2026-11-08"
  ), tz = "UTC")
)

# Summarise trading activity
analyse_spike <- function(start, end, period) {
  
  data <- trades |>
    filter(datetime >= start, datetime < end)
  
  if (nrow(data) == 0) return(NULL)
  
  trade_volumes <- sort(data$volume_usd, decreasing = TRUE)
  total_volume <- sum(trade_volumes)
  
  wallet_volumes <- data |>
    group_by(proxy_wallet) |>
    summarise(
      volume = sum(volume_usd),
      .groups = "drop"
    )
  
  data.frame(
    period = period,
    n_trades = nrow(data),
    total_volume = total_volume,
    mean_trade = mean(trade_volumes),
    median_trade = median(trade_volumes),
    max_trade = max(trade_volumes),
    top3_share = sum(head(trade_volumes, 3)) / total_volume,
    n_wallets = n_distinct(data$proxy_wallet),
    top_wallet_share = max(wallet_volumes$volume) / total_volume,
    no_share = sum(
      data$volume_usd[data$outcome == "No"]
    ) / total_volume,
    buy_no_share = sum(
      data$volume_usd[
        data$outcome == "No" & data$side == "BUY"
      ]
    ) / total_volume
  )
}

# Compare all periods
spike_summary <- bind_rows(
  lapply(seq_len(nrow(windows)), function(i) {
    analyse_spike(
      windows$start[i],
      windows$end[i],
      windows$period[i]
    )
  })
)

print(spike_summary)

# Largest trades per period
top_trades <- bind_rows(
  lapply(seq_len(nrow(windows)), function(i) {
    
    trades |>
      filter(
        datetime >= windows$start[i],
        datetime < windows$end[i]
      ) |>
      arrange(desc(volume_usd)) |>
      slice_head(n = 5) |>
      mutate(period = windows$period[i]) |>
      select(
        period, datetime, volume_usd,
        outcome, side, proxy_wallet
      )
  })
)

print(top_trades, n = Inf)

# Volume by outcome and side
volume_split <- bind_rows(
  lapply(seq_len(nrow(windows)), function(i) {
    
    trades |>
      filter(
        datetime >= windows$start[i],
        datetime < windows$end[i]
      ) |>
      group_by(outcome, side) |>
      summarise(
        n_trades = n(),
        total_volume = sum(volume_usd),
        .groups = "drop"
      ) |>
      mutate(period = windows$period[i]) |>
      select(period, everything())
  })
)

print(volume_split, n = Inf)

# Save results
dir.create(
  "results/spike_analysis",
  recursive = TRUE,
  showWarnings = FALSE
)

write_csv(
  spike_summary,
  "results/spike_analysis/spike_summary.csv"
)

write_csv(
  top_trades,
  "results/spike_analysis/top_trades.csv"
)

write_csv(
  volume_split,
  "results/spike_analysis/volume_split.csv"
)

saveRDS(
  list(
    windows = windows,
    spike_summary = spike_summary,
    top_trades = top_trades,
    volume_split = volume_split
  ),
  "results/spike_analysis/results.rds"
)
