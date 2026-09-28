# ==============================================================================
# truth.R
# What the simulation studies score against: a plug-in true Rt path built from
# a fitted series, and the trailing mean that defines a sliding-window
# estimator's estimand.
# ==============================================================================

# Interpolate a fitted Rt series onto every day 1..n_days.
# rule = 2 extends the first/last fitted value flat over the leading and
# trailing days that no sliding window reached, so the simulator always has a
# defined Rt.
build_true_rt <- function(index, value, n_days) {
  usable <- is.finite(index) & is.finite(value) & value > 0

  if (sum(usable) < 2L) {
    stop("Not enough finite Rt values to build a plug-in true trajectory.")
  }

  stats::approx(
    x = index[usable],
    y = value[usable],
    xout = seq_len(n_days),
    rule = 2
  )$y
}

# Trailing w-day mean: element t is mean(x[(t - w + 1):t]), NA for t < w.
#
# This is the estimand of a w-day sliding-window estimator. fit_epiestim_full()
# sets t_start = t_end - w + 1 and EpiEstim assumes Rt is constant across that
# window, so the quantity it estimates for day t_end is the AVERAGE Rt over
# [t_end - w + 1, t_end], not the instantaneous Rt on day t_end. The two differ
# whenever the truth is trending, which is exactly what makes the windowed arm
# undercover the instantaneous truth in the coverage study. w = 1 returns x
# unchanged, i.e. the instantaneous target.
trailing_mean <- function(x, w) {
  w <- as.integer(w)
  if (length(w) != 1L || is.na(w) || w < 1L) {
    stop("w must be one positive integer.")
  }
  if (w > length(x)) {
    stop("w must not exceed the length of x.")
  }
  if (w == 1L) return(as.numeric(x))

  cumulative <- cumsum(c(0, x))
  out <- rep(NA_real_, length(x))
  ends <- seq.int(w, length(x))
  out[ends] <- (cumulative[ends + 1L] - cumulative[ends + 1L - w]) / w
  out
}
