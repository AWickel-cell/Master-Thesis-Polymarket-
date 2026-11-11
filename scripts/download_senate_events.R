library(rvest)
library(dplyr)
library(stringr)
library(purrr)
library(readr)
library(tibble)
library(lubridate)

base_url <- "https://allthecalls.app"

# Get Senate race pages
page <- read_html(paste0(base_url, "/senate"))

race_urls <- page |>
  html_elements("a[href^='/races/senate/']") |>
  html_attr("href") |>
  unique()

cat("Senate races found:", length(race_urls), "\n")

# Extract CNN calls
get_cnn_call <- function(path) {
  Sys.sleep(0.3)
  
  tables <- read_html(paste0(base_url, path)) |>
    html_elements("table") |>
    html_table(fill = TRUE)
  
  calls <- keep(tables, ~ all(
    c("Source", "Status and time") %in% names(.x)
  ))
  
  if (length(calls) == 0) return(NULL)
  
  cnn <- calls[[1]] |>
    filter(str_detect(Source, "CNN"))
  
  if (nrow(cnn) == 0) return(NULL)
  
  tibble(
    state = str_remove(basename(path), "-senate$"),
    source = "CNN",
    status = cnn[["Status and time"]]
  )
}

events_raw <- map_dfr(race_urls, get_cnn_call) |>
  mutate(
    party = str_remove_all(
      str_extract(status, "\\([DRI]\\)"), "[()]"
    ),
    time_et = str_extract(
      status, "\\d{1,2}:\\d{2}\\s*[AP]M\\s*EST"
    ),
    called = str_detect(status, "^Called:")
  )

dir.create("data", showWarnings = FALSE)
write_csv(events_raw, "data/senate_calls_raw.csv")

cat("Races extracted:", nrow(events_raw), "\n")
cat("CNN calls:", sum(events_raw$called), "\n")

# Prepare Senate events
senate_events <- events_raw |>
  filter(called) |>
  mutate(
    date_et = if_else(
      str_detect(time_et, "AM"),
      "2026-11-04",
      "2026-11-03"
    ),
    datetime_et = as.POSIXct(
      paste(date_et, time_et),
      format = "%Y-%m-%d %I:%M %p EST",
      tz = "America/New_York"
    ),
    time_utc = with_tz(datetime_et, "UTC")
  ) |>
  select(state, party, time_utc, source) |>
  arrange(time_utc)

cat("Missing timestamps:", sum(is.na(senate_events$time_utc)), "\n")

write_csv(senate_events, "data/senate_events.csv")

print(senate_events, n = Inf)
