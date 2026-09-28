# ==============================================================================
# metrics.R
# Scoring functions: Rt credible-interval coverage (two flavours, for the two
# shapes the studies produce), single-forecast summaries, and backtest
# performance metrics.
# ==============================================================================

# ------------------------------------------------------------------------------
# Long-format Rt coverage (main.R Section 6, sim_main.R, truth_source_main.R)
# ------------------------------------------------------------------------------

# Coverage of the 95% Rt credible interval, scored against TWO targets:
#
#   true_R_instant - the instantaneous plug-in true Rt on day t
#   true_R_window  - its trailing epiestim_window-day average, the estimand of
#                    the primary sliding-window arm
#
# Every arm is scored against both so the table shows where an arm's shortfall
# comes from a target mismatch rather than from miscalibration. true_R_own is
# whichever of the two is the arm's own estimand, and the own-estimand coverage
# is the one that actually tests calibration. Coverage is averaged within a
# replicate first and across replicates second, so the across-replicate SD is a
# valid Monte Carlo standard error (this matches the simulation study in epi.R).
rt_coverage_metrics <- function(data) {
  complete <- data %>%
    filter(
      is.finite(R), is.finite(lower), is.finite(upper),
      is.finite(true_R_instant), is.finite(true_R_window), is.finite(true_R_own)
    )

  if (nrow(complete) == 0L) {
    return(data.frame(
      N_replicates = 0L, N_day_estimates = 0L, Estimand = NA_character_,
      Coverage95_Rt_ownEstimand = NA_real_,
      MCSE_Coverage95_Rt_ownEstimand = NA_real_,
      Coverage95_Rt_instantaneous = NA_real_,
      Coverage95_Rt_windowAvg = NA_real_,
      MeanCIWidth = NA_real_,
      Bias_ownEstimand = NA_real_, RMSE_ownEstimand = NA_real_,
      check.names = FALSE
    ))
  }

  estimand_window <- unique(complete$estimand_window)
  if (length(estimand_window) != 1L) {
    stop("rt_coverage_metrics() expects one estimand window per arm.")
  }

  per_replicate <- complete %>%
    group_by(replicate) %>%
    summarise(
      coverage_own = mean(lower <= true_R_own & true_R_own <= upper),
      coverage_instant = mean(lower <= true_R_instant & true_R_instant <= upper),
      coverage_window = mean(lower <= true_R_window & true_R_window <= upper),
      .groups = "drop"
    )

  data.frame(
    N_replicates = nrow(per_replicate),
    N_day_estimates = nrow(complete),
    Estimand = if (estimand_window == 1L) {
      "instantaneous Rt"
    } else {
      paste0(estimand_window, "-day average Rt")
    },
    Coverage95_Rt_ownEstimand = mean(per_replicate$coverage_own),
    MCSE_Coverage95_Rt_ownEstimand = if (nrow(per_replicate) > 1L) {
      stats::sd(per_replicate$coverage_own) / sqrt(nrow(per_replicate))
    } else {
      NA_real_
    },
    Coverage95_Rt_instantaneous = mean(per_replicate$coverage_instant),
    Coverage95_Rt_windowAvg = mean(per_replicate$coverage_window),
    MeanCIWidth = mean(complete$upper - complete$lower),
    Bias_ownEstimand = mean(complete$R - complete$true_R_own),
    RMSE_ownEstimand = sqrt(mean((complete$R - complete$true_R_own)^2)),
    check.names = FALSE
  )
}

# ------------------------------------------------------------------------------
# Per-replicate fit lists (the SI-misspecification studies, the selfcheck's twin)
# ------------------------------------------------------------------------------

# Score a list of per-replicate fits against one known truth.
#
# `fits` is a list, one element per replicate, each a list with `R` and `time`
# (the day index each estimate belongs to) plus the interval bounds named in
# `intervals`, and `error` - NA on success. Fits whose `error` is set are
# dropped; NULL is returned if none survive.
#
# `intervals` names the interval(s) to score. The first is the PRIMARY: it
# defines Coverage95 / MeanCIWidth and the finite-value mask that Bias / MSE
# are computed on. Any further entry is scored alongside under its own name,
# e.g. list(ci = c("lower", "upper"), hpd = c("hpd_lower", "hpd_upper")) adds
# CoverageHPD95 / MCSE_CoverageHPD / MeanHPDWidth and by-day hpd_coverage /
# hpd_half, so interval SHAPE and interval CALIBRATION can be separated.
#
# `extra` is a named list of function(fits) -> scalar for per-study summaries
# that are not interval metrics (rho_hat, pct_converged, mean_fit_secs).
#
# Coverage is computed PER REPLICATE first, so the MCSE across replicates is
# valid; the day-level table is the transpose view, across replicates per day.
# half_over_sd is the interval half-width over the estimator's actual spread
# across replicates: 1.96 is the calibrated value, below it the interval is too
# narrow.
#
# By default a single non-finite day turns the day-level summaries
# (half_over_sd, days/pct_bias_exceeds_half) into NA, as the archived scorers
# did. finite_days_only = TRUE leaves such days out of them instead, as `good`
# leaves them out of the per-replicate metrics. The misspec studies score from
# day 2 and need it: EpiEstim returns no estimate for a window ending on or
# before the assumed SI's mean.
score_interval_fits <- function(fits, truth_scored, score_days,
                                intervals = list(ci = c("lower", "upper")),
                                extra = list(), finite_days_only = FALSE) {
  ok <- vapply(fits, function(f) is.null(f$error) || is.na(f$error), logical(1))
  fits <- fits[ok]
  if (length(fits) == 0L) return(NULL)

  if (length(intervals) == 0L || is.null(names(intervals)) ||
      any(!nzchar(names(intervals)))) {
    stop("intervals must be a named list of c(lower_field, upper_field).")
  }
  primary <- names(intervals)[1L]
  secondary <- names(intervals)[-1L]
  n_days <- length(score_days)

  aligned <- function(f, field) f[[field]][match(score_days, f$time)]

  # --- per-replicate first, so the MCSE across replicates is valid ---
  per_rep <- t(vapply(fits, function(f) {
    est <- aligned(f, "R")
    lo <- aligned(f, intervals[[primary]][1L])
    hi <- aligned(f, intervals[[primary]][2L])
    good <- is.finite(lo) & is.finite(hi) & is.finite(est)

    out <- c(
      coverage = mean((lo <= truth_scored & truth_scored <= hi)[good]),
      width    = mean((hi - lo)[good])
    )
    for (nm in secondary) {
      slo <- aligned(f, intervals[[nm]][1L])
      shi <- aligned(f, intervals[[nm]][2L])
      sgood <- is.finite(slo) & is.finite(shi)
      out[[paste0(nm, "_coverage")]] <-
        mean((slo <= truth_scored & truth_scored <= shi)[sgood])
      out[[paste0(nm, "_width")]] <- mean((shi - slo)[sgood])
    }
    out[["bias"]]   <- mean((est - truth_scored)[good])
    out[["sq_err"]] <- mean(((est - truth_scored)^2)[good])
    out
  }, numeric(4L + 2L * length(secondary))))

  # --- day-level, across replicates ---
  # Matrices are days x replicates, so recycling truth_scored down the columns
  # compares each replicate's interval against the right day.
  as_matrix <- function(field) {
    vapply(fits, function(f) aligned(f, field), numeric(n_days))
  }
  est_mat <- as_matrix("R")
  lo_mat  <- as_matrix(intervals[[primary]][1L])
  hi_mat  <- as_matrix(intervals[[primary]][2L])

  day_tbl <- data.frame(
    day = score_days,
    truth = truth_scored,
    mean_R = rowMeans(est_mat),
    sd_R = apply(est_mat, 1L, stats::sd),
    coverage = rowMeans(lo_mat <= truth_scored & truth_scored <= hi_mat),
    half = rowMeans(hi_mat - lo_mat) / 2
  )
  for (nm in secondary) {
    slo_mat <- as_matrix(intervals[[nm]][1L])
    shi_mat <- as_matrix(intervals[[nm]][2L])
    day_tbl[[paste0(nm, "_coverage")]] <-
      rowMeans(slo_mat <= truth_scored & truth_scored <= shi_mat)
    day_tbl[[paste0(nm, "_half")]] <- rowMeans(shi_mat - slo_mat) / 2
  }
  day_tbl$bias <- day_tbl$mean_R - day_tbl$truth

  mse <- mean(per_rep[, "sq_err"])
  n_rep <- nrow(per_rep)

  metrics <- data.frame(
    n_replicates = n_rep,
    n_failed_fits = sum(!ok),
    Coverage95 = mean(per_rep[, "coverage"]),
    MCSE_Coverage = stats::sd(per_rep[, "coverage"]) / sqrt(n_rep),
    MeanCIWidth = mean(per_rep[, "width"]),
    stringsAsFactors = FALSE
  )
  for (nm in secondary) {
    tag <- toupper(nm)
    metrics[[paste0("Coverage", tag, "95")]] <- mean(per_rep[, paste0(nm, "_coverage")])
    metrics[[paste0("MCSE_Coverage", tag)]] <-
      stats::sd(per_rep[, paste0(nm, "_coverage")]) / sqrt(n_rep)
    metrics[[paste0("Mean", tag, "Width")]] <- mean(per_rep[, paste0(nm, "_width")])
  }
  metrics$Bias <- mean(per_rep[, "bias"])
  metrics$MSE <- mse
  metrics$RMSE <- sqrt(mse)
  finite_day <- if (isTRUE(finite_days_only)) {
    is.finite(day_tbl$half) & is.finite(day_tbl$sd_R) & is.finite(day_tbl$bias)
  } else {
    rep(TRUE, nrow(day_tbl))
  }
  metrics$half_over_sd <- mean(day_tbl$half[finite_day]) / mean(day_tbl$sd_R[finite_day])
  for (nm in names(extra)) {
    metrics[[nm]] <- extra[[nm]](fits)
  }
  exceeds <- (abs(day_tbl$bias) > day_tbl$half)[finite_day]
  metrics$days_bias_exceeds_half <- sum(exceeds)
  metrics$pct_bias_exceeds_half <- 100 * mean(exceeds)

  list(metrics = metrics, by_day = day_tbl)
}

# ------------------------------------------------------------------------------
# Forecast scoring (main.R Section 7)
# ------------------------------------------------------------------------------

# One forecast summary from posterior Rt draws and a common Poisson renewal model.
make_forecast_summary <- function(R_draws, infectiousness, actual, n_draws) {
  if (
    length(R_draws) != n_draws ||
    any(!is.finite(R_draws)) ||
    !is.finite(infectiousness) ||
    infectiousness < 0
  ) {
    return(c(
      pred = NA_real_, lower = NA_real_, upper = NA_real_,
      log_score = NA_real_, interval_score = NA_real_
    ))
  }

  lambda_draws <- pmax(R_draws * infectiousness, 1e-12)
  incidence_draws <- stats::rpois(n_draws, lambda = lambda_draws)

  lower <- unname(stats::quantile(incidence_draws, 0.025, names = FALSE))
  upper <- unname(stats::quantile(incidence_draws, 0.975, names = FALSE))

  # Posterior predictive density:
  # p(I[t+1] = actual | data through t)
  log_score <- log_mean_exp(
    stats::dpois(actual, lambda = lambda_draws, log = TRUE)
  )

  alpha <- 0.05
  interval_score <- (upper - lower) +
    (2 / alpha) * (lower - actual) * as.numeric(actual < lower) +
    (2 / alpha) * (actual - upper) * as.numeric(actual > upper)

  c(
    pred = mean(lambda_draws),
    lower = lower,
    upper = upper,
    log_score = log_score,
    interval_score = interval_score
  )
}

# All metrics below score the one-step-ahead forecast of INCIDENCE I[t+1],
# not the estimated Rt itself. In particular Coverage95_incidence asks how
# often the observed case count fell inside the 95% posterior predictive
# interval - it is NOT the coverage of an Rt credible interval. Rt coverage
# requires a known truth and is measured in the coverage study on semi-synthetic
# replicates calibrated to this series (and in the simulation study, epi.R).
forecast_metrics <- function(data) {
  complete <- data %>%
    filter(
      is.finite(actual), is.finite(pred), is.finite(lower), is.finite(upper),
      is.finite(log_score), is.finite(interval_score)
    )

  if (nrow(complete) == 0L) {
    return(data.frame(
      N = 0L, MAE = NA_real_, RMSE = NA_real_,
      `Bias (prediction - actual)` = NA_real_,
      MeanLogScore = NA_real_, Coverage95_incidence = NA_real_,
      MeanIntervalWidth = NA_real_, MeanIntervalScore = NA_real_,
      check.names = FALSE
    ))
  }

  data.frame(
    N = nrow(complete),
    MAE = mean(abs(complete$pred - complete$actual)),
    RMSE = sqrt(mean((complete$pred - complete$actual)^2)),
    `Bias (prediction - actual)` = mean(complete$pred - complete$actual),
    MeanLogScore = mean(complete$log_score),
    Coverage95_incidence = mean(
      complete$actual >= complete$lower & complete$actual <= complete$upper
    ),
    MeanIntervalWidth = mean(complete$upper - complete$lower),
    MeanIntervalScore = mean(complete$interval_score),
    check.names = FALSE
  )
}
