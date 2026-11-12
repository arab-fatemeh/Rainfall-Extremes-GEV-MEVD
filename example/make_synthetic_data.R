# Makes two synthetic daily rainfall series (150 years each) to try the code.
# These are NOT real observations.
#   Station A: stationary climate
#   Station B: the Weibull scale of wet-day rainfall grows by 20% over the last 60 years
# Wet days follow a seasonal probability; wet-day amounts follow a Weibull distribution.

set.seed(1)
make_station <- function(years, C_trend = 0) {
  out <- NULL
  for (y in years) {
    days <- seq(as.Date(paste0(y, "-01-01")), as.Date(paste0(y, "-12-31")), by = "day")
    doy  <- as.numeric(format(days, "%j"))
    p_wet <- 0.25 + 0.10 * cos(2 * pi * (doy - 300) / 365)          # wetter in autumn
    wet  <- runif(length(days)) < p_wet
    frac <- max(0, (y - (max(years) - 60)) / 60)                     # 0 before the trend starts
    C    <- 7 * (1 + C_trend * frac)
    rain <- ifelse(wet, 1 + rweibull(length(days), shape = 0.8, scale = C), 0)
    out  <- rbind(out, data.frame(date = format(days), rain_mm = round(rain, 1)))
  }
  out
}

dir.create("data", showWarnings = FALSE)
years <- 1871:2020
write.csv(make_station(years, 0),   "data/example_station_A.csv", row.names = FALSE)
write.csv(make_station(years, 0.2), "data/example_station_B.csv", row.names = FALSE)
cat("Written data/example_station_A.csv and data/example_station_B.csv\n")
