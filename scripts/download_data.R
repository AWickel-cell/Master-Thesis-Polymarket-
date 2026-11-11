library(httr2)
library(jsonlite)
library(dplyr)

condition_id <- "0x48b0b0bca515f68fccf95af4793dbd0edbfec1f8ec6e8df2c0f69ba74f8c4722"
api_v2 <- "https://data-api.polymarket.com/v2/trades"

cursor <- NULL
all_pages <- list()

repeat {
  req <- request(api_v2) |>
    req_url_query(condition = condition_id, limit = 100) |>
    req_retry(max_tries = 5)
  
  if (!is.null(cursor)) {
    req <- req |> req_url_query(cursor = cursor)
  }
  
  result <- req |>
    req_perform() |>
    resp_body_string() |>
    fromJSON()
  
  if (is.null(result$data)) {
    stop("Download incomplete: missing trade data.")
  }
  
  all_pages[[length(all_pages) + 1]] <- result$data
  
  if (length(all_pages) %% 20 == 0) {
    cat("Downloaded:",
        sum(vapply(all_pages, nrow, integer(1))),
        "trades\n")
  }
  
  if (identical(result$pagination$has_more, FALSE)) break
  
  next_cursor <- result$pagination$next_cursor
  
  if (is.null(next_cursor) || length(next_cursor) != 1 ||
      is.na(next_cursor) || !nzchar(next_cursor) ||
      identical(next_cursor, cursor)) {
    stop("Download incomplete: invalid pagination cursor.")
  }
  
  cursor <- next_cursor
  Sys.sleep(0.2)
}

trades_v2 <- bind_rows(all_pages)

stopifnot(
  nrow(trades_v2) > 0,
  all(trades_v2$condition_id == condition_id)
)

dir.create("data", showWarnings = FALSE)

write.csv(trades_v2, "data/trades_raw_v2.csv", row.names = FALSE)

cat("Download complete\n")
cat("Pages:", length(all_pages), "\n")
cat("Trades:", nrow(trades_v2), "\n")
cat("Period:",
    format(as.POSIXct(min(trades_v2$timestamp),
                      origin = "1970-01-01", tz = "UTC")),
    "to",
    format(as.POSIXct(max(trades_v2$timestamp),
                      origin = "1970-01-01", tz = "UTC")),
    "\n")
