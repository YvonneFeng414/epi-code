################################################################################
# Unit checks for the SI family constructors in R/si/serial_interval.R.
#
# The three misspecification studies fit under five SI families. Everything
# downstream of make_si_family() assumes the returned vector is a proper pmf in
# the make_si() layout, and the family comparison assumes the non-gamma families
# are discretized the same way as gamma_bin. This script checks both, plus the
# one place the construction is known to be lossy (uniform at the grid extremes).
#
# Run from r-proj/:  Rscript tests/test_si_families.R
# Exits non-zero on the first failure, so it can gate a run.
################################################################################

suppressPackageStartupMessages({
  library(EpiLPS)
  library(EpiEstim)
})

if (!file.exists(file.path("R", "si", "serial_interval.R"))) {
  stop("Run this from the r-proj/ directory.")
}
source(file.path("R", "si", "serial_interval.R"))

mean_si <- 7.5
sd_si   <- 3.4
max_si_lag <- 30L

n_fail <- 0L
check <- function(label, pass, detail = "") {
  cat(sprintf("  %-58s %s%s\n", label, if (pass) "PASS" else "FAIL",
              if (nzchar(detail)) paste0("  ", detail) else ""))
  if (!pass) n_fail <<- n_fail + 1L
  invisible(pass)
}

# ==============================================================================
# 1. Shape of the returned object
# ==============================================================================
# Every consumer indexes si[k + 1] = P(SI = k) and passes si[-1] to the
# estimator, so a family that got the lag-0 cell wrong would shift the whole
# distribution a day without any visible error.

cat("\n===== 1. pmf layout =====\n")
for (family in si_families) {
  si <- make_si_family(family, mean_si, sd_si, max_si_lag)
  check(sprintf("%s: sums to 1, si[1] == 0, no negatives", family),
        abs(sum(si) - 1) < 1e-12 && si[1L] == 0 && all(si >= 0) && all(is.finite(si)))
}

# ==============================================================================
# 2. Moment matching
# ==============================================================================
# Each family is parameterised from (mean, sd), so at the true SI every one of
# them should realise those moments. Binning and truncation cost a little; 0.05
# is the same tolerance the study scripts' own off-by-one guard uses.

cat("\n===== 2. realised moments at the true SI =====\n")
for (family in si_families) {
  si <- make_si_family(family, mean_si, sd_si, max_si_lag)
  m <- si_moments(si)
  check(sprintf("%s: mean %.4f, sd %.4f", family, m[["mean"]], m[["sd"]]),
        abs(m[["mean"]] - mean_si) < 0.05 && abs(m[["sd"]] - sd_si) < 0.05)
}

# ==============================================================================
# 3. Agreement with EpiLPS::Idist
# ==============================================================================
# Idist is the independent implementation of the same +/-0.5 binning, so it is
# the oracle for gamma, lognormal and Weibull (it has no uniform). It picks its
# own support length, and the comparison has to be made at THAT length: at a
# common lag 30 the two differ by whatever mass sits outside the shorter
# support, which is a truncation difference and not a disagreement about the
# discretization. Checking at matched Dmax isolates the formula.

cat("\n===== 3. vs EpiLPS::Idist, at matched support length =====\n")
for (pair in list(c("gamma_bin", "gamma"), c("lnorm", "lognorm"),
                  c("weibull", "weibull"))) {
  theirs <- c(0, EpiLPS::Idist(mean = mean_si, sd = sd_si, dist = pair[2L])$pvec)
  mine <- make_si_family(pair[1L], mean_si, sd_si, length(theirs) - 1L)
  tv <- 0.5 * sum(abs(mine - theirs))
  check(sprintf("%s: TV = %.3e at Dmax = %d", pair[1L], tv, length(theirs) - 1L),
        tv < 1e-5)
}

# ==============================================================================
# 4. gamma_discr is untouched, and the discretization gap is what we think
# ==============================================================================
# gamma_discr must stay bit-identical to make_si(): it is the data-generating
# SI, and the whole paired design rests on the simulated series not moving.
# The gamma_discr/gamma_bin gap is the reason gamma_bin exists as a separate
# reference; compare_si_discretization.R measures the same 0.019.

cat("\n===== 4. gamma_discr and the discretization gap =====\n")
check("gamma_discr identical to make_si()",
      identical(make_si_family("gamma_discr", mean_si, sd_si, max_si_lag),
                make_si(mean_si, sd_si, max_si_lag)))

tv_discr <- 0.5 * sum(abs(make_si_family("gamma_discr", mean_si, sd_si, max_si_lag) -
                          make_si_family("gamma_bin", mean_si, sd_si, max_si_lag)))
check(sprintf("TV(gamma_discr, gamma_bin) = %.4f, near the documented 0.019", tv_discr),
      abs(tv_discr - 0.019) < 0.002)

# ==============================================================================
# 5. The whole grid constructs
# ==============================================================================
# 5 families x 17 (mean, sd) settings is what the studies fit. A family that
# fails to construct at one grid point would otherwise surface as a missing key
# deep inside the fit loop.

cat("\n===== 5. every family x setting constructs =====\n")
relative_error <- c(-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8)
grid <- unique(rbind(
  data.frame(mean = mean_si, sd = round(sd_si * (1 + relative_error), 2)),
  data.frame(mean = round(mean_si * (1 + relative_error), 2), sd = sd_si)
))

bad <- 0L
for (family in si_families) {
  for (i in seq_len(nrow(grid))) {
    si <- tryCatch(make_si_family(family, grid$mean[i], grid$sd[i], max_si_lag),
                   error = function(e) NULL)
    if (is.null(si) || abs(sum(si) - 1) > 1e-10) bad <- bad + 1L
  }
}
check(sprintf("%d family x setting combinations construct",
              length(si_families) * nrow(grid)), bad == 0L)

# ==============================================================================
# 6. Where the families actually differ, and where uniform stops meaning it
# ==============================================================================
# Not pass/fail. These two tables are the reason the experiment is worth running
# and the caveat that has to travel with the uniform arm, so they are printed
# rather than asserted.

cat("\n===== 6. lag-1 mass and support at the true SI =====\n")
cat(sprintf("  %-12s %10s %10s %10s\n", "family", "lag-1", "lag-2", "support"))
for (family in si_families) {
  si <- make_si_family(family, mean_si, sd_si, max_si_lag)
  nz <- which(si > 1e-10) - 1L
  cat(sprintf("  %-12s %10.5f %10.5f %6d-%-3d\n",
              family, si[2L], si[3L], min(nz), max(nz)))
}
cat("  (this spread is the effect the family arm is designed to measure)\n")

cat("\n===== 7. uniform moment drift across the grid =====\n")
cat("  Uniform support is mean +/- sd*sqrt(3); where that falls below lag 0.5 the\n")
cat("  truncated, renormalised pmf no longer has the requested moments.\n\n")
cat(sprintf("  %8s %6s | %9s %8s | %8s %8s\n",
            "req.mean", "req.sd", "real.mean", "real.sd", "d_mean", "d_sd"))
drifted <- 0L
for (i in seq_len(nrow(grid))) {
  si <- make_si_family("unif", grid$mean[i], grid$sd[i], max_si_lag)
  m <- si_moments(si)
  d_mean <- m[["mean"]] - grid$mean[i]
  d_sd <- m[["sd"]] - grid$sd[i]
  flag <- max(abs(d_mean), abs(d_sd)) > 0.25
  if (flag) drifted <- drifted + 1L
  cat(sprintf("  %8.2f %6.2f | %9.3f %8.3f | %+8.3f %+8.3f%s\n",
              grid$mean[i], grid$sd[i], m[["mean"]], m[["sd"]], d_mean, d_sd,
              if (flag) "  <- drifted" else ""))
}
cat(sprintf("\n  %d of %d uniform settings drift by more than 0.25.\n",
            drifted, nrow(grid)))
cat("  They are kept and flagged, not dropped: the study CSVs carry realised_mean,\n")
cat("  realised_sd and drift_ok so the affected points stay visible.\n")

# ==============================================================================
cat("\n==============================================================\n")
if (n_fail > 0L) {
  cat(sprintf("%d CHECK(S) FAILED\n", n_fail))
  cat("==============================================================\n")
  quit(status = 1L)
}
cat("All checks passed.\n")
cat("==============================================================\n")
