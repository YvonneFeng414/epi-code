################################################################################
# Equivalence checks for the helpers that were consolidated out of the study
# scripts. Each check carries the ARCHIVED implementation as its oracle - the
# function body as it stood in archive/epi/queens_rt/ - and asserts that the
# shared function reproduces it: bit-for-bit for the simulators (same seed, same
# draws), to 1e-12 for the scorers, on a synthetic fixture.
#
# Run from r-proj/:  Rscript tests/test_shared_helpers.R
# Exits non-zero on the first failure, so it can gate a run.
################################################################################

source(file.path("R", "core", "init.R"))
require_packages(c("EpiEstim", "EpiLPS", "dplyr", "tidyr", "ggplot2", "purrr"))
source_project()

n_fail <- 0L
check <- function(label, pass, detail = "") {
  cat(sprintf("  %-64s %s%s\n", label, if (pass) "PASS" else "FAIL",
              if (nzchar(detail)) paste0("  ", detail) else ""))
  if (!pass) n_fail <<- n_fail + 1L
  invisible(pass)
}

# Largest elementwise difference; NA in the same positions counts as equal,
# NA in different positions as infinitely different.
max_abs_diff <- function(a, b) {
  a <- as.numeric(a); b <- as.numeric(b)
  if (length(a) != length(b) || !identical(is.na(a), is.na(b))) return(Inf)
  d <- abs(a - b)
  if (all(is.na(d))) 0 else max(d, na.rm = TRUE)
}

# ==============================================================================
# 1. The renewal simulator
# ==============================================================================
# Five copies became one. Under the same seed the unified function must return
# the same integers as each original, because the loop arithmetic and the
# sequence of RNG calls are unchanged.

# archive/epi/queens_rt/R_epi_lsp/simulation.R
simulate_renewal_incidence_old <- function(seed_incidence, R_true, si_full, n_days) {
  max_lag <- length(si_full) - 1L
  n_seed <- length(seed_incidence)
  incidence <- numeric(n_days)
  incidence[seq_len(n_seed)] <- seed_incidence
  for (t in seq.int(n_seed + 1L, n_days)) {
    lags <- seq_len(min(max_lag, t - 1L))
    lambda <- R_true[t] * sum(incidence[t - lags] * si_full[lags + 1L])
    if (!is.finite(lambda) || lambda > 1e9) return(NULL)
    incidence[t] <- stats::rpois(1L, lambda = max(lambda, 0))
  }
  if (anyNA(incidence)) return(NULL)
  as.integer(incidence)
}

# archive/epi/queens_rt/epilps_si_misspec.R section 6 (also the MALA twin)
simulate_series_nbinom_old <- function(R, si, seed_incidence, rho) {
  n  <- length(R)
  ns <- length(seed_incidence)
  I  <- numeric(n)
  I[seq_len(ns)] <- seed_incidence
  for (t in seq.int(ns + 1L, n)) {
    k  <- seq_len(min(length(si) - 1L, t - 1L))
    mu <- R[t] * sum(I[t - k] * si[k + 1L])
    if (!is.finite(mu) || mu > 1e9) return(NULL)
    I[t] <- stats::rnbinom(1L, mu = max(mu, 0), size = rho)
  }
  if (anyNA(I)) return(NULL)
  as.integer(I)
}

# archive/epi/queens_rt/epifilter_si_misspec.R section 6
simulate_series_poisson_old <- function(R, si, seed_incidence) {
  n <- length(R)
  ns <- length(seed_incidence)
  I <- numeric(n)
  I[seq_len(ns)] <- seed_incidence
  for (t in seq.int(ns + 1L, n)) {
    k <- seq_len(min(length(si) - 1L, t - 1L))
    mu <- R[t] * sum(I[t - k] * si[k + 1L])
    if (!is.finite(mu) || mu > 1e9) return(NULL)
    I[t] <- stats::rpois(1L, lambda = max(mu, 0))
  }
  if (anyNA(I)) return(NULL)
  as.integer(I)
}

cat("\n===== 1. simulate_renewal_incidence() vs the archived simulators =====\n")

si <- make_si(7.5, 3.4, 30L)
n_days <- 150L
R_true <- 1 + 0.4 * sin(2 * pi * seq_len(n_days) / 60)
seed_inc <- as.integer(round(200 + 30 * sin(seq_len(30))))
rho <- 17.3

sim_pair <- function(seed, old, new) {
  set.seed(seed); a <- lapply(1:20, function(r) old())
  set.seed(seed); b <- lapply(1:20, function(r) new())
  identical(a, b)
}

check("Poisson: identical to simulation.R's simulate_renewal_incidence()",
      sim_pair(1L,
               function() simulate_renewal_incidence_old(seed_inc, R_true, si, n_days),
               function() simulate_renewal_incidence(seed_inc, R_true, si, n_days)))

check("Poisson: identical to epifilter_si_misspec.R's simulate_series()",
      sim_pair(2L,
               function() simulate_series_poisson_old(R_true, si, seed_inc),
               function() simulate_renewal_incidence(seed_inc, R_true, si,
                                                     obs_model = "poisson")))

check("NegBin: identical to epilps_si_misspec.R's simulate_series()",
      sim_pair(3L,
               function() simulate_series_nbinom_old(R_true, si, seed_inc, rho),
               function() simulate_renewal_incidence(seed_inc, R_true, si,
                                                     obs_model = "nbinom",
                                                     overdispersion = rho)))

check("NegBin without an overdispersion is refused",
      inherits(tryCatch(simulate_renewal_incidence(seed_inc, R_true, si, obs_model = "nbinom"),
                        error = function(e) e), "error"))

# ==============================================================================
# 2. The scorer
# ==============================================================================
# score_setting() (MAP), its HPD-carrying twin (MALA) and score_arm() (EpiFilter)
# were three copies of one computation. On a fixture of synthetic fits the
# shared score_interval_fits() must reproduce every column each of them
# produced, including the day-level table.

# archive/epi/queens_rt/epilps_si_misspec.R section 9, with its closed-over
# globals (score_days, truth_scored) made explicit.
score_setting_old <- function(fits, truth_scored, score_days) {
  ok <- vapply(fits, function(f) is.null(f$error) || is.na(f$error), logical(1))
  fits <- fits[ok]
  if (length(fits) == 0L) return(NULL)
  per_rep <- t(vapply(fits, function(f) {
    idx <- match(score_days, f$time)
    lo <- f$lower[idx]; hi <- f$upper[idx]; est <- f$R[idx]
    good <- is.finite(lo) & is.finite(hi) & is.finite(est)
    c(coverage = mean((lo <= truth_scored & truth_scored <= hi)[good]),
      width    = mean((hi - lo)[good]),
      bias     = mean((est - truth_scored)[good]),
      sq_err   = mean(((est - truth_scored)^2)[good]))
  }, numeric(4)))
  est_mat <- vapply(fits, function(f) f$R[match(score_days, f$time)], numeric(length(score_days)))
  lo_mat  <- vapply(fits, function(f) f$lower[match(score_days, f$time)], numeric(length(score_days)))
  hi_mat  <- vapply(fits, function(f) f$upper[match(score_days, f$time)], numeric(length(score_days)))
  day_tbl <- data.frame(
    day = score_days, truth = truth_scored, mean_R = rowMeans(est_mat),
    coverage = rowMeans(lo_mat <= truth_scored & truth_scored <= hi_mat),
    half = rowMeans(hi_mat - lo_mat) / 2
  )
  day_tbl$bias <- day_tbl$mean_R - day_tbl$truth
  mse <- mean(per_rep[, "sq_err"])
  list(
    metrics = data.frame(
      n_replicates = nrow(per_rep), n_failed_fits = sum(!ok),
      Coverage95 = mean(per_rep[, "coverage"]),
      MCSE_Coverage = stats::sd(per_rep[, "coverage"]) / sqrt(nrow(per_rep)),
      MeanCIWidth = mean(per_rep[, "width"]),
      Bias = mean(per_rep[, "bias"]), MSE = mse, RMSE = sqrt(mse),
      rho_hat = mean(vapply(fits, function(f) f$rho, numeric(1))),
      days_bias_exceeds_half = sum(abs(day_tbl$bias) > day_tbl$half),
      pct_bias_exceeds_half = 100 * mean(abs(day_tbl$bias) > day_tbl$half),
      stringsAsFactors = FALSE
    ),
    by_day = day_tbl
  )
}

# archive/epi/queens_rt/epilps_si_misspec_mala.R section 9
score_setting_mala_old <- function(fits, truth_scored, score_days) {
  ok <- vapply(fits, function(f) is.null(f$error) || is.na(f$error), logical(1))
  fits <- fits[ok]
  if (length(fits) == 0L) return(NULL)
  per_rep <- t(vapply(fits, function(f) {
    idx <- match(score_days, f$time)
    lo <- f$lower[idx]; hi <- f$upper[idx]; est <- f$R[idx]
    hlo <- f$hpd_lower[idx]; hhi <- f$hpd_upper[idx]
    good  <- is.finite(lo) & is.finite(hi) & is.finite(est)
    hgood <- is.finite(hlo) & is.finite(hhi)
    c(coverage     = mean((lo <= truth_scored & truth_scored <= hi)[good]),
      width        = mean((hi - lo)[good]),
      hpd_coverage = mean((hlo <= truth_scored & truth_scored <= hhi)[hgood]),
      hpd_width    = mean((hhi - hlo)[hgood]),
      bias         = mean((est - truth_scored)[good]),
      sq_err       = mean(((est - truth_scored)^2)[good]))
  }, numeric(6)))
  n <- length(score_days)
  est_mat <- vapply(fits, function(f) f$R[match(score_days, f$time)], numeric(n))
  lo_mat  <- vapply(fits, function(f) f$lower[match(score_days, f$time)], numeric(n))
  hi_mat  <- vapply(fits, function(f) f$upper[match(score_days, f$time)], numeric(n))
  hlo_mat <- vapply(fits, function(f) f$hpd_lower[match(score_days, f$time)], numeric(n))
  hhi_mat <- vapply(fits, function(f) f$hpd_upper[match(score_days, f$time)], numeric(n))
  day_tbl <- data.frame(
    day = score_days, truth = truth_scored, mean_R = rowMeans(est_mat),
    coverage = rowMeans(lo_mat <= truth_scored & truth_scored <= hi_mat),
    half = rowMeans(hi_mat - lo_mat) / 2,
    hpd_coverage = rowMeans(hlo_mat <= truth_scored & truth_scored <= hhi_mat),
    hpd_half = rowMeans(hhi_mat - hlo_mat) / 2
  )
  day_tbl$bias <- day_tbl$mean_R - day_tbl$truth
  mse <- mean(per_rep[, "sq_err"])
  list(
    metrics = data.frame(
      n_replicates = nrow(per_rep), n_failed_fits = sum(!ok),
      Coverage95 = mean(per_rep[, "coverage"]),
      MCSE_Coverage = stats::sd(per_rep[, "coverage"]) / sqrt(nrow(per_rep)),
      MeanCIWidth = mean(per_rep[, "width"]),
      CoverageHPD95 = mean(per_rep[, "hpd_coverage"]),
      MCSE_CoverageHPD = stats::sd(per_rep[, "hpd_coverage"]) / sqrt(nrow(per_rep)),
      MeanHPDWidth = mean(per_rep[, "hpd_width"]),
      Bias = mean(per_rep[, "bias"]), MSE = mse, RMSE = sqrt(mse),
      rho_hat = mean(vapply(fits, function(f) f$rho, numeric(1))),
      pct_converged = 100 * mean(vapply(fits, function(f) isTRUE(f$converged), logical(1))),
      mean_fit_secs = mean(vapply(fits, function(f) f$secs, numeric(1))),
      days_bias_exceeds_half = sum(abs(day_tbl$bias) > day_tbl$half),
      pct_bias_exceeds_half = 100 * mean(abs(day_tbl$bias) > day_tbl$half),
      stringsAsFactors = FALSE
    ),
    by_day = day_tbl
  )
}

# archive/epi/queens_rt/epifilter_si_misspec.R section 8. Its fits were already
# aligned to score_days and named lo / hi.
score_arm_old <- function(fits, truth_scored, score_days) {
  ok <- vapply(fits, function(f) is.null(f$error) || is.na(f$error), logical(1))
  fits <- fits[ok]
  if (length(fits) == 0L) return(NULL)
  per_rep <- t(vapply(fits, function(f) {
    lo <- f$lo; hi <- f$hi; est <- f$R
    good <- is.finite(lo) & is.finite(hi) & is.finite(est)
    c(coverage = mean((lo <= truth_scored & truth_scored <= hi)[good]),
      width    = mean((hi - lo)[good]),
      bias     = mean((est - truth_scored)[good]),
      sq_err   = mean(((est - truth_scored)^2)[good]))
  }, numeric(4)))
  est_mat <- vapply(fits, function(f) f$R, numeric(length(truth_scored)))
  lo_mat  <- vapply(fits, function(f) f$lo, numeric(length(truth_scored)))
  hi_mat  <- vapply(fits, function(f) f$hi, numeric(length(truth_scored)))
  day_tbl <- data.frame(
    day = score_days, truth = truth_scored, mean_R = rowMeans(est_mat),
    sd_R = apply(est_mat, 1L, stats::sd),
    coverage = rowMeans(lo_mat <= truth_scored & truth_scored <= hi_mat),
    half = rowMeans(hi_mat - lo_mat) / 2
  )
  day_tbl$bias <- day_tbl$mean_R - day_tbl$truth
  mse <- mean(per_rep[, "sq_err"])
  list(
    metrics = data.frame(
      n_replicates = nrow(per_rep), n_failed = sum(!ok),
      Coverage95 = mean(per_rep[, "coverage"]),
      MCSE_Coverage = stats::sd(per_rep[, "coverage"]) / sqrt(nrow(per_rep)),
      MeanCIWidth = mean(per_rep[, "width"]),
      Bias = mean(per_rep[, "bias"]), MSE = mse, RMSE = sqrt(mse),
      half_over_sd = mean(day_tbl$half) / mean(day_tbl$sd_R),
      days_bias_exceeds_half = sum(abs(day_tbl$bias) > day_tbl$half),
      pct_bias_exceeds_half = 100 * mean(abs(day_tbl$bias) > day_tbl$half),
      stringsAsFactors = FALSE
    ),
    by_day = day_tbl
  )
}

cat("\n===== 2. score_interval_fits() vs the archived scorers =====\n")

# A fixture: 8 replicates of a 60-day fit scored on days 31..60, one replicate
# failed, a few non-finite bounds, HPD intervals slightly narrower and shifted.
set.seed(20260906L)
score_days <- 31:60
truth_scored <- 1 + 0.2 * cos(seq_along(score_days) / 5)
make_fit <- function(r) {
  time <- 1:60
  R <- 1 + 0.2 * cos(time / 5) + stats::rnorm(60, 0, 0.05)
  half <- stats::runif(60, 0.05, 0.15)
  lower <- R - half; upper <- R + half
  hpd_lower <- R - 0.9 * half + 0.01; hpd_upper <- R + 0.9 * half + 0.01
  if (r == 3L) { lower[40] <- NA; upper[40] <- NA; hpd_lower[45] <- NA }
  if (r == 5L) return(list(R = NULL, error = "fit failed"))
  list(R = R, lower = lower, upper = upper, hpd_lower = hpd_lower,
       hpd_upper = hpd_upper, time = time, rho = 15 + r, converged = r != 2L,
       secs = 20 + r, error = NA_character_)
}
fits <- lapply(1:8, make_fit)

same_cols <- function(a, b, cols, tol = 1e-12) {
  all(vapply(cols, function(cc) max_abs_diff(a[[cc]], b[[cc]]) <= tol, logical(1)))
}

rho_hat <- function(f) mean(vapply(f, function(x) x$rho, numeric(1)))

new_map <- score_interval_fits(fits, truth_scored, score_days,
                               extra = list(rho_hat = rho_hat))
old_map <- score_setting_old(fits, truth_scored, score_days)
check("MAP: every archived metrics column reproduced",
      same_cols(new_map$metrics, old_map$metrics, names(old_map$metrics)),
      sprintf("%d columns", ncol(old_map$metrics)))
check("MAP: every archived by-day column reproduced",
      same_cols(new_map$by_day, old_map$by_day, names(old_map$by_day)))

new_mala <- score_interval_fits(
  fits, truth_scored, score_days,
  intervals = list(ci = c("lower", "upper"), hpd = c("hpd_lower", "hpd_upper")),
  extra = list(
    rho_hat = rho_hat,
    pct_converged = function(f) 100 * mean(vapply(f, function(x) isTRUE(x$converged), logical(1))),
    mean_fit_secs = function(f) mean(vapply(f, function(x) x$secs, numeric(1)))
  )
)
old_mala <- score_setting_mala_old(fits, truth_scored, score_days)
check("MALA: every archived metrics column reproduced (incl. HPD)",
      same_cols(new_mala$metrics, old_mala$metrics, names(old_mala$metrics)),
      sprintf("%d columns", ncol(old_mala$metrics)))
check("MALA: every archived by-day column reproduced (incl. HPD)",
      same_cols(new_mala$by_day, old_mala$by_day, names(old_mala$by_day)))

# EpiFilter's fits were pre-aligned to score_days with lo/hi names; build the
# equivalent in the new layout and the old one side by side.
ef_new <- lapply(fits, function(f) {
  if (is.null(f$R)) return(f)
  idx <- match(score_days, f$time)
  list(R = f$R[idx], lower = f$lower[idx], upper = f$upper[idx],
       time = score_days, error = NA_character_)
})
ef_old <- lapply(ef_new, function(f) {
  if (is.null(f$R)) return(f)
  list(R = f$R, lo = f$lower, hi = f$upper, error = NA_character_)
})
new_ef <- score_interval_fits(ef_new, truth_scored, score_days)
old_ef <- score_arm_old(ef_old, truth_scored, score_days)
check("EpiFilter: every archived metrics column reproduced (n_failed -> n_failed_fits)",
      same_cols(new_ef$metrics, old_ef$metrics,
                setdiff(names(old_ef$metrics), "n_failed")) &&
        new_ef$metrics$n_failed_fits == old_ef$metrics$n_failed)
check("EpiFilter: every archived by-day column reproduced (incl. sd_R)",
      same_cols(new_ef$by_day, old_ef$by_day, names(old_ef$by_day)))

check("all fits failed returns NULL",
      is.null(score_interval_fits(list(list(R = NULL, error = "x")), truth_scored, score_days)))

# ==============================================================================
# 3. fit_epiestim_last() vs the last row of fit_epiestim_full()
# ==============================================================================
# main.R Section 7A used to call estimate_R() inline for the window ending at
# each forecast origin. EpiEstim's per-window posteriors are independent of the
# other windows, so the last row of the full fit is the same posterior.

cat("\n===== 3. fit_epiestim_last() vs fit_epiestim_full() =====\n")

inc <- read_incidence_csv(project_paths()$incidence_csv)
for (w in c(1L, 7L)) {
  for (t in c(80L, 150L, nrow(inc))) {
    sub <- inc[seq_len(t), , drop = FALSE]
    last <- fit_epiestim_last(sub, si, w)
    full <- fit_epiestim_full(sub, si, w)
    full <- full[nrow(full), ]
    check(sprintf("window %d, origin %d: mean/sd/quantiles agree", w, t),
          last$status == "OK" &&
            max_abs_diff(c(last$R_mean, last$R_sd, last$q025, last$q50, last$q975),
                         c(full$R_mean, full$R_sd, full$lower, full$R, full$upper)) < 1e-12)
  }
}

# ==============================================================================
# 4. run_epifilter() vs the archived inline version
# ==============================================================================
# The archived run_epifilter() closed over six script globals; the shared one
# takes them as arguments and names the bounds lower/upper. Same recursion, so
# the numbers must match on the observed series.

cat("\n===== 4. run_epifilter() vs epifilter_si_misspec.R's inline version =====\n")

run_epifilter_old <- function(incidence, si, Rgrid, m_grid, eta_value, pR0, alpha,
                              filter_start, days, arms = c("smoother", "filter")) {
  tryCatch({
    w <- si[-1L]; n <- length(incidence); L <- numeric(n)
    for (i in seq.int(2L, n)) { k <- seq_len(min(length(w), i - 1L)); L[i] <- sum(incidence[i - k] * w[k]) }
    td <- seq.int(filter_start, length(incidence)); nd <- length(td)
    Rf <- epiFilter(Rgrid, m_grid, eta_value, pR0, nd, L[td], incidence[td], alpha)
    Rs <- epiSmoother(Rgrid, m_grid, Rf[[4L]], Rf[[5L]], nd, Rf[[6L]], alpha)
    idx <- match(days, td)
    pick <- function(o) list(R = o[[1L]][idx], lo = o[[2L]][1L, idx], hi = o[[2L]][2L, idx])
    out <- list(filter = pick(Rf), smoother = pick(Rs), days = days, error = NA_character_)
    if (anyNA(unlist(lapply(out[arms], unlist)))) return(list(error = "NaN in the recursion"))
    out
  }, error = function(e) list(error = conditionMessage(e)))
}

observed <- inc$I[16:nrow(inc)]
days <- 31:length(observed)
grid <- make_epifilter_grid(0.01, 10, 1000L, 0.025)
new_ef_fit <- run_epifilter(observed, si, grid, 0.1, 21L, days)
old_ef_fit <- run_epifilter_old(observed, si, grid$Rgrid, grid$m, 0.1, grid$pR0,
                                grid$alpha, 21L, days)
for (arm in c("smoother", "filter")) {
  check(sprintf("%s: R / lower / upper identical", arm),
        is.na(new_ef_fit$error) && is.na(old_ef_fit$error) &&
          identical(new_ef_fit[[arm]]$R, old_ef_fit[[arm]]$R) &&
          identical(new_ef_fit[[arm]]$lower, old_ef_fit[[arm]]$lo) &&
          identical(new_ef_fit[[arm]]$upper, old_ef_fit[[arm]]$hi))
}
check("time field is the requested days", identical(new_ef_fit$smoother$time, as.integer(days)))

# ==============================================================================
# 5. fit_epiestim_grid() - the harness adapter - vs fit_epiestim_full()
# ==============================================================================
# The adapter is what lets EpiEstim plug into fit_si_grid() beside EpiLPS and
# EpiFilter. Its whole contract: one arm per window in the common layout, with
# the same numbers fit_epiestim_full() returns and `time` equal to its `index`;
# a fitting failure reported as list(error = ...) rather than thrown.

cat("\n===== 5. fit_epiestim_grid() vs fit_epiestim_full() =====\n")

origin <- as.Date("2020-01-01")
obs_df <- data.frame(dates = origin + seq_along(observed) - 1L, I = observed)
grid_fit <- fit_epiestim_grid(observed, si, windows = c(1L, 7L), origin = origin)
check("adapter returns without error", is.na(grid_fit$error))
check("one arm per window, named w1 / w7",
      identical(setdiff(names(grid_fit), "error"), c("w1", "w7")))
for (w in c(1L, 7L)) {
  arm <- grid_fit[[paste0("w", w)]]
  full <- fit_epiestim_full(obs_df, si, w)
  check(sprintf("w%d: R / lower / upper identical to fit_epiestim_full()", w),
        is.na(arm$error) &&
          identical(arm$R, full$R_median) &&
          identical(arm$lower, full$lower) &&
          identical(arm$upper, full$upper))
  check(sprintf("w%d: time is fit_epiestim_full()'s index", w),
        identical(arm$time, as.integer(full$index)))
  check(sprintf("w%d: scored days all present", w),
        !anyNA(match(days, arm$time)))
}

# Too few observations for the window: fit_epiestim_full() stops; the adapter
# must turn that into an error field, since on a worker a throw would lose the
# whole replicate.
too_short <- fit_epiestim_grid(observed[1:3], si, windows = 7L)
check("failure comes back as list(error = message), not a throw",
      is.list(too_short) && !is.null(too_short$error) && !is.na(too_short$error) &&
        is.null(too_short$w7))

# ==============================================================================
cat("\n==============================================================\n")
if (n_fail > 0L) {
  cat(sprintf("%d CHECK(S) FAILED\n", n_fail))
  cat("==============================================================\n")
  quit(status = 1L)
}
cat("All checks passed.\n")
cat("==============================================================\n")
