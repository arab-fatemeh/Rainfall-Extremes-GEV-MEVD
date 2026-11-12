# Functions for the moving-window analysis of daily rainfall extremes
# GEV (L-moments) vs MEVD with Weibull ordinary events (PWM)
# Fatemeh Arab Sadeghabadi, 2026

# ---- settings (can be changed before sourcing run_analysis.R) ----
if (!exists("Tr"))     Tr     <- 200   # return period (years)
if (!exists("win"))    win    <- 30    # window length (years)
if (!exists("minyrs")) minyrs <- 20    # minimum valid years in a window
if (!exists("thr"))    thr    <- 1     # wet day: rain > thr (mm)
if (!exists("nboot"))  nboot  <- 500   # bootstrap samples
if (!exists("mindays")) mindays <- 360 # a year is valid with at least this many days
Fr <- 1 - 1 / Tr


# ---- readers: all return a list with one vector of daily rain (mm) per year ----

# generic CSV with columns: date (YYYY-MM-DD) and rain_mm
read_daily_csv <- function(f) {
  d <- read.csv(f)
  yr <- as.integer(substr(d$date, 1, 4))
  split(d$rain_mm, yr)
}

# GHCN-Daily CSV export (columns DATE as YYYYMMDD and PRCP in tenths of mm)
read_ghcn <- function(f) {
  d <- read.csv(f)
  split(d$PRCP / 10, d$DATE %/% 10000)
}

# plain text with blocks: a year number, then the daily values of that year
read_year_blocks <- function(f) {
  v <- scan(f, quiet = TRUE)
  pos <- which(v >= 1700 & v <= 2100 & v %% 1 == 0)
  out <- list()
  for (k in seq_along(pos)) {
    last <- if (k < length(pos)) pos[k + 1] - 1 else length(v)
    out[[as.character(v[pos[k]])]] <- v[(pos[k] + 1):last]
  }
  out
}


# ---- GEV with L-moments (Hosking 1990; Hosking et al. 1985) ----

gev_lmom <- function(x) {
  x <- sort(x); n <- length(x); j <- 1:n
  b0 <- mean(x)
  b1 <- sum((j - 1) / (n - 1) * x) / n
  b2 <- sum((j - 1) * (j - 2) / ((n - 1) * (n - 2)) * x) / n
  l1 <- b0; l2 <- 2 * b1 - b0; l3 <- 6 * b2 - 6 * b1 + b0
  t3 <- l3 / l2
  cc <- 2 / (3 + t3) - log(2) / log(3)
  k  <- 7.8590 * cc + 2.9554 * cc^2
  sc <- l2 * k / ((1 - 2^(-k)) * gamma(1 + k))
  lo <- l1 - sc * (1 - gamma(1 + k)) / k
  c(shape = -k, scale = sc, loc = lo)
}

gev_quant <- function(p) {
  unname(p["loc"] + p["scale"] / p["shape"] * ((-log(Fr))^(-p["shape"]) - 1))
}


# ---- MEVD with Weibull ordinary events (Marani & Ignaccolo 2015; Zorzetto et al. 2016) ----

# Weibull fit with probability weighted moments (Greenwood et al. 1979)
wei_pwm <- function(x) {
  x <- sort(x); n <- length(x); i <- 1:n
  M0 <- mean(x)
  M1 <- sum(x * (n - i)) / (n * (n - 1))
  w  <- log(2) / log(M0 / (2 * M1))
  C  <- M0 / gamma(1 + 1 / w)
  c(C = C, w = w)
}

# quantile of the MEVD: solve mean_j [1 - exp(-(x/C_j)^w_j)]^N_j = F
mev_quant <- function(N, C, w) {
  f <- function(x) mean((1 - exp(-(x / C)^w))^N) - Fr
  uniroot(f, c(1e-6, 5000), tol = 1e-8)$root + thr      # add the threshold back
}


# ---- one station: yearly values, moving windows and bootstrap ----

yearly_values <- function(dat) {
  yy <- NULL
  for (y in as.numeric(names(dat))) {
    r <- dat[[as.character(y)]]
    r <- r[!is.na(r)]
    if (length(r) < mindays) next
    e <- r[r > thr] - thr                                # ordinary events (excesses)
    p <- wei_pwm(e)
    yy <- rbind(yy, data.frame(year = y, amax = max(r), N = length(e), C = p["C"], w = p["w"]))
  }
  rownames(yy) <- NULL
  yy
}

analyse <- function(dat) {
  yy <- yearly_values(dat)
  res <- NULL
  for (end in (min(yy$year) + win - 1):max(yy$year)) {
    s <- yy[yy$year > end - win & yy$year <= end, ]
    if (nrow(s) < minyrs) next
    n <- nrow(s)
    gev <- gev_quant(gev_lmom(s$amax))
    mev <- mev_quant(s$N, s$C, s$w)
    gb <- mb <- numeric(nboot)
    for (b in 1:nboot) {
      id <- sample(n, n, replace = TRUE)                 # same years for both methods
      gb[b] <- gev_quant(gev_lmom(s$amax[id]))
      mb[b] <- mev_quant(s$N[id], s$C[id], s$w[id])
    }
    res <- rbind(res, data.frame(end_year = end, n_years = n,
                                 gev = gev, gev_lo = quantile(gb, 0.05), gev_hi = quantile(gb, 0.95),
                                 mev = mev, mev_lo = quantile(mb, 0.05), mev_hi = quantile(mb, 0.95)))
  }
  rownames(res) <- NULL
  list(windows = res, years = yy)
}


# ---- plots ----

pieces <- function(x) split(seq_along(x), cumsum(c(1, diff(x) != 1)))   # gaps break the lines
band <- function(x, lo, hi, col) {
  for (k in pieces(x)) polygon(c(x[k], rev(x[k])), c(lo[k], rev(hi[k])), col = col, border = NA)
}
gline <- function(x, y, ...) for (k in pieces(x)) lines(x[k], y[k], ...)

RED <- "#b2182b"; BLUE <- "#2166ac"

plot_return_levels <- function(results, file) {
  xl <- range(unlist(lapply(results, function(r) r$windows$end_year)))
  yl <- range(unlist(lapply(results, function(r) unlist(r$windows[, c("gev_lo", "gev_hi", "mev_lo", "mev_hi")]))))
  png(file, width = 2000, height = 800 * length(results), res = 220)
  par(mfrow = c(length(results), 1), mar = c(3.5, 4.5, 2, 1), las = 1)
  for (nm in names(results)) {
    r <- results[[nm]]$windows
    plot(NA, xlim = xl, ylim = yl, xlab = "", ylab = paste0(Tr, "-year daily rainfall (mm)"),
         main = nm, adj = 0, font.main = 2)
    grid(nx = NA, ny = NULL, col = "grey90", lty = 1)
    band(r$end_year, r$gev_lo, r$gev_hi, adjustcolor(RED, 0.18))
    band(r$end_year, r$mev_lo, r$mev_hi, adjustcolor(BLUE, 0.30))
    gline(r$end_year, r$gev, col = RED, lwd = 1.3)
    gline(r$end_year, r$mev, col = BLUE, lwd = 1.6)
    if (nm == names(results)[1])
      legend("topleft", bty = "n", cex = 0.8, legend = c("GEV", "MEVD-Weibull"),
             col = c(RED, BLUE), lwd = c(1.3, 1.6))
  }
  mtext(paste0("Final year of the ", win, "-year window"), side = 1, line = 2.3)
  invisible(dev.off())
}

plot_interval_width <- function(results, file) {
  xl <- range(unlist(lapply(results, function(r) r$windows$end_year)))
  png(file, width = 2000, height = 800, res = 220)
  par(mar = c(4, 4.5, 1, 1), las = 1)
  wmax <- max(unlist(lapply(results, function(r) with(r$windows, c(gev_hi - gev_lo, mev_hi - mev_lo)))))
  plot(NA, xlim = xl, ylim = c(0, wmax), xlab = paste0("Final year of the ", win, "-year window"),
       ylab = "Width of 90% interval (mm)")
  lt <- seq_along(results)
  for (i in lt) {
    r <- results[[i]]$windows
    gline(r$end_year, r$gev_hi - r$gev_lo, col = RED, lty = i)
    gline(r$end_year, r$mev_hi - r$mev_lo, col = BLUE, lty = i, lwd = 1.5)
  }
  legend("topleft", bty = "n", cex = 0.8, ncol = 2,
         legend = c(paste("GEV,", names(results)), paste("MEVD,", names(results))),
         col = rep(c(RED, BLUE), each = length(results)), lty = rep(lt, 2))
  invisible(dev.off())
}
