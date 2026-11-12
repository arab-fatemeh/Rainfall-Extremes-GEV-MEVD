# Moving-window analysis of daily rainfall extremes: GEV vs MEVD
# Usage (from the repository folder):
#   Rscript R/run_analysis.R                 -> runs on the synthetic example data
#   edit STATIONS below to use your own files
# Fatemeh Arab Sadeghabadi, 2026

set.seed(2026)
source("R/functions.R")

# ---- stations: name = list(file, reader) ----
STATIONS <- list(
  "Station A (synthetic)" = list(file = "data/example_station_A.csv", reader = read_daily_csv),
  "Station B (synthetic)" = list(file = "data/example_station_B.csv", reader = read_daily_csv)
)
# Examples for other formats:
#   "My station" = list(file = "data/my_ghcn_export.csv", reader = read_ghcn)
#   "My station" = list(file = "data/my_year_blocks.txt",  reader = read_year_blocks)

if (!all(file.exists(sapply(STATIONS, `[[`, "file"))))
  stop("Data files not found. Run  Rscript example/make_synthetic_data.R  first, or edit STATIONS.")

results <- lapply(STATIONS, function(s) analyse(s$reader(s$file)))

dir.create("results", showWarnings = FALSE)
for (nm in names(results)) {
  out <- file.path("results", paste0(gsub("[^A-Za-z0-9]+", "_", nm), "_windows.csv"))
  write.csv(results[[nm]]$windows, out, row.names = FALSE)
}

dir.create("figures", showWarnings = FALSE)
plot_return_levels(results, "figures/return_levels.png")
plot_interval_width(results, "figures/interval_width.png")

for (nm in names(results)) {
  r <- results[[nm]]$windows
  cat(sprintf("%s: %d windows | GEV %.0f-%.0f mm | MEVD %.0f-%.0f mm | median 90%% width: GEV %.0f, MEVD %.0f mm\n",
              nm, nrow(r), min(r$gev), max(r$gev), min(r$mev), max(r$mev),
              median(r$gev_hi - r$gev_lo), median(r$mev_hi - r$mev_lo)))
}
