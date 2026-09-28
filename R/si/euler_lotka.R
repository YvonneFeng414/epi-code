# ==============================================================================
# euler_lotka.R
# The discrete Euler-Lotka equation, connecting a serial interval to the
# exponential growth rate it implies and back. Ported from
# archive/Epiestim_simu/euler_lotka_functions.R, where it explains the
# asymptotic bias of EpiEstim under a wrong SI; here it serves the Queens
# misspecification study the same way.
#
# Base R only, and SI vectors in the project layout of serial_interval.R:
# si[1] is the lag-0 cell (zero), si[k + 1] = P(SI = k). The lag vector `s` is
# derived from that layout by default, so callers pass the SI as they built it.
#
# The tight uniroot tolerance is load-bearing and is kept from the archive:
# uniroot's default (.Machine$double.eps^0.25, about 1.2e-4) left the true-SI
# point at 1.19999 rather than 1.2, which contaminated every bias number
# downstream. A copy of this function that drifted from that tolerance would do
# the same silently, which is why the archive kept one copy and so does this.
# ==============================================================================

# The growth rate r implied by R under serial interval w:
#
#     sum_s w(s) exp(-r s) = 1 / R
#
# The left side decreases in r, so the root is unique. The archive searched
# (1e-9, 2), which assumes R > 1; the Queens truth path drops below 1 after the
# lockdown, so the bracket is widened to negative r and extended on demand.
# With w[1] = 0 the left side at r = 5 is below exp(-5), which brackets any
# R < 148; at the lower end it exceeds exp(|r| * mean lag), which brackets any
# R the studies can meet from the other side.
#
# Only lags with positive mass enter the sum (0 * Inf is NaN), and the lower
# end of the bracket is pulled in so exp(-r * s) stays finite at the longest
# lag: the Britton check hands in 400-lag vectors, where exp(2 * 399)
# overflows. exp(700) is the largest safe magnitude.
discrete_growth_rate <- function(w, R, s = seq_along(w) - 1L,
                                 interval = c(-2, 5)) {
  if (length(w) != length(s)) stop("w and s must have the same length.")
  if (!is.finite(R) || R <= 0) stop("R must be one positive finite value.")
  keep <- w > 0
  w <- w[keep]
  s <- s[keep]
  lower <- max(interval[1L], -700 / max(s))
  stats::uniroot(
    function(r) sum(w * exp(-r * s)) - 1 / R,
    interval = c(lower, interval[2L]),
    extendInt = "yes",
    tol = .Machine$double.eps^0.9
  )$root
}

# The R that an estimator told the SI is `w` converges on once incidence grows
# at rate r and the start-of-series boundary has been left behind:
#
#     R_EL = 1 / sum_s w(s) exp(-r s)
#
# Under the generating SI this returns the R that produced r, so the difference
# R_EL(assumed) - R_EL(true) is the asymptotic misspecification bias.
# Vectorised over r so a whole path is one call. Zero-mass lags are dropped
# for the same reason as above.
euler_lotka_R <- function(w, r, s = seq_along(w) - 1L) {
  if (length(w) != length(s)) stop("w and s must have the same length.")
  keep <- w > 0
  w <- w[keep]
  s <- s[keep]
  1 / vapply(r, function(ri) sum(w * exp(-ri * s)), numeric(1))
}

# The growth rate on every day of an Rt path, under the true SI. Computed once
# per study and shared across every assumed SI, since it depends only on the
# truth. NA where the path is NA (a trailing-mean estimand's leading days).
implied_growth_rates <- function(R_true, si_true) {
  vapply(R_true, function(R) {
    if (is.na(R)) NA_real_ else discrete_growth_rate(si_true, R)
  }, numeric(1))
}

# The Euler-Lotka prediction for a whole path: the R an estimator would report
# on each day under `si_assumed`, if on each day the epidemic were growing at
# that day's implied rate. This is the asymptotic result applied locally, so it
# is an approximation wherever the truth is trending, and exact only where it
# is flat. At si_assumed == si_true it returns R_true (to the uniroot
# tolerance), which is the check to run before trusting the rest.
euler_lotka_path <- function(si_assumed, growth_rates) {
  out <- rep(NA_real_, length(growth_rates))
  ok <- !is.na(growth_rates)
  out[ok] <- euler_lotka_R(si_assumed, growth_rates[ok])
  out
}
