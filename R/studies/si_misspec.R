# ==============================================================================
# si_misspec.R
# Shared harness for the serial-interval misspecification studies
# (epiestim_si_misspec.R, epilps_si_misspec.R, epilps_si_misspec_mala.R,
# epifilter_si_misspec.R).
#
# All four run the same experiment with a different estimator plugged in:
# simulate replicate epidemics ONCE from a plug-in truth under the TRUE SI,
# then refit every replicate under every assumed SI on a grid that crosses two
# dimensions -
#
#   MOMENTS - two sweeps sharing the correctly specified point, the assumed
#             mean or sd scaled by the same relative errors (config.R's
#             misspec grid, -60%..+60%)
#   FAMILY  - the SHAPE of the assumed SI at those moments: gamma_discr (the
#             generating SI), gamma_bin, lnorm, weibull, unif
#
# - and score coverage / width / bias / RMSE per setting. The design is PAIRED:
# one set of simulated epidemics is refit under every assumed SI, so differences
# across the grid are pure misspecification effect.
#
# The functions here are the steps of that experiment, in the order a driver
# calls them. What differs between drivers - the fit function, the observation
# model, extra intervals or arms, the figure captions - comes in as arguments.
#
# Requires si/serial_interval.R, simulation/renewal.R, scoring/metrics.R,
# data/incidence_io.R and core/parallel.R to be sourced first.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Inputs
# ------------------------------------------------------------------------------

# The observed series, and the plug-in truth the study scores against.
#
# The truth (plugin_truth_paths.csv, from truth_source_main.R) covers the
# "common range" that pipeline selected, which starts some days into the
# observed series. That offset used to be a hard-coded common_start <- 16L in
# every study; it is derived here from the truth file's first date and, when
# the caller passes the documented value, checked against it.
#
# With truth_required = FALSE no truth column is taken (the EpiFilter study
# builds its own truth from its own fit), and common_start must be supplied
# unless the file is read anyway for the hole days (below).
#
# Seeding is one of two kinds. With `burn_in`, the first burn_in days of the
# common range are the seed the simulator copies verbatim. With `I0`, the seed
# is I0 cases on day 1 only (burn_in becomes length(I0)), as in the reference
# study Epiestim_simu/seed/run_sensitivity_queens_nohole.R. Either way scoring
# starts on the day after the seed.
#
# With `hole_threshold`, the days whose EpiEstim 1-day truth (true_R_epiestim)
# falls below it are removed from R_true - the reference study's rule for the
# 2020-07-08 reporting hole. The EpiEstim column locates them whichever truth
# is in use, so every driver drops the same calendar day. `observed` is NOT
# cut: it is the data a truth is fitted on, and the hole is in the data. So
# with a hole removed, n_days is one less than length(observed), and a truth
# built from `observed` must drop hole_days itself.
load_misspec_inputs <- function(paths, truth_column = "true_R_epilps", burn_in = NULL,
                                truth_scheme = "discr_si", common_start = NULL,
                                truth_required = TRUE, I0 = NULL,
                                hole_threshold = NULL) {
  if (is.null(I0) == is.null(burn_in)) {
    stop("Give exactly one of burn_in (observed seed days) or I0 (day-1 seed).")
  }

  obs <- read_incidence_csv(paths$incidence_csv)

  truth <- NULL
  R_true <- NULL
  if (isTRUE(truth_required) || !is.null(hole_threshold)) {
    truth_file <- paths$truth_paths_csv(truth_scheme)
    if (!file.exists(truth_file)) {
      stop("Cannot find ", normalizePath(truth_file, mustWork = FALSE),
           ". Run truth_source_main.R first.")
    }
    truth <- read.csv(truth_file, stringsAsFactors = FALSE)
    needed <- c(if (isTRUE(truth_required)) truth_column,
                if (!is.null(hole_threshold)) "true_R_epiestim")
    missing_cols <- setdiff(needed, names(truth))
    if (length(missing_cols) > 0L) {
      stop("Column(s) ", paste(missing_cols, collapse = ", "), " not in ", truth_file, ".")
    }
    if (isTRUE(truth_required)) R_true <- truth[[truth_column]]

    derived_start <- match(as.Date(truth$date[1L]), obs$dates)
    if (is.na(derived_start)) {
      stop("The truth file's first date is not in the incidence series; the two ",
           "files come from different runs.")
    }
    if (!is.null(common_start) && derived_start != common_start) {
      stop("common_start derived from the truth file is ", derived_start,
           " but the documented value is ", common_start, ".")
    }
    common_start <- derived_start
  } else if (is.null(common_start)) {
    stop("common_start must be given when no truth file is read.")
  }

  observed <- obs$I[seq.int(common_start, nrow(obs))]
  if (!is.null(truth) && length(observed) != nrow(truth)) {
    stop("The truth has ", nrow(truth), " days but the common range has ",
         length(observed), "; the two files come from different runs.")
  }

  hole_days <- integer(0)
  if (!is.null(hole_threshold)) {
    hole_days <- which(truth$true_R_epiestim < hole_threshold)
  }
  keep_days <- setdiff(seq_along(observed), hole_days)
  if (!is.null(R_true)) R_true <- R_true[keep_days]
  n_days <- length(keep_days)

  if (!is.null(I0)) {
    seed_incidence <- as.integer(I0)
    burn_in <- length(seed_incidence)
  } else {
    seed_incidence <- observed[keep_days][seq_len(burn_in)]
  }
  if (burn_in >= n_days) {
    stop("The seed leaves no days to score.")
  }
  score_from <- burn_in + 1L
  score_days <- seq.int(score_from, n_days)

  list(
    obs = obs,
    observed = observed,
    R_true = R_true,
    truth_column = if (is.null(R_true)) NA_character_ else truth_column,
    n_days = n_days,
    common_start = common_start,
    hole_days = hole_days,
    hole_dates = obs$dates[common_start + hole_days - 1L],
    hole_truth = if (is.null(truth)) numeric(0) else truth$true_R_epiestim[hole_days],
    burn_in = burn_in,
    I0 = I0,
    seed_incidence = seed_incidence,
    score_from = score_from,
    score_days = score_days,
    truth_scored = if (is.null(R_true)) NULL else R_true[score_days]
  )
}

# The days removed from the truth, in the form the reference study prints them.
report_hole_days <- function(inp) {
  if (length(inp$hole_days) == 0L) return(invisible(NULL))
  cat("Days removed from the true Rt trajectory (reporting hole):\n")
  for (k in seq_along(inp$hole_days)) {
    d <- inp$hole_days[k]
    cat(sprintf("  day %d (%s): EpiEstim true Rt %.5f, reported cases that day %d\n",
                d, as.character(inp$hole_dates[k]), inp$hole_truth[k],
                as.integer(inp$observed[d])))
  }
  invisible(NULL)
}

# The seed line of the input report, for either kind of seeding.
seed_label <- function(inp) {
  if (!is.null(inp$I0)) {
    sprintf("I0 = %d cases on day 1 only", inp$I0)
  } else {
    sprintf("first %d observed days of the common range (%d to %d cases/day)",
            inp$burn_in, min(inp$seed_incidence), max(inp$seed_incidence))
  }
}

report_misspec_inputs <- function(inp, n_sim,
                                  note = "shared across all SI settings") {
  report_hole_days(inp)
  if (!is.null(inp$R_true)) {
    cat(sprintf("Truth: %s, %d days, Rt %.3f to %.3f\n",
                inp$truth_column, inp$n_days, min(inp$R_true), max(inp$R_true)))
  }
  cat(sprintf("Seed: %s\n", seed_label(inp)))
  cat(sprintf("Scoring days %d-%d (%d days), %d replicates, %s\n",
              inp$score_from, inp$n_days, length(inp$score_days), n_sim, note))
  invisible(NULL)
}

# "-60%..+60%", for report headings that name the relative-error range.
rel_error_range_label <- function(relative_error) {
  r <- 100 * range(relative_error)
  sprintf("%+.0f%%..%+.0f%%", r[1L], r[2L])
}

# ------------------------------------------------------------------------------
# 2. Observation model calibration
# ------------------------------------------------------------------------------

# The NegBin overdispersion EpiLPS (LPSMAP) estimates on the observed series,
# fitted in a child process. This is the data-generating value for the NegBin
# studies; it is taken from the MAP fit even in the MALA study, because
# changing it would change the simulated series and break the pairing.
estimate_overdispersion <- function(incidence, si, K) {
  fit <- run_in_child(fit_epilps_map, list(incidence = incidence, si = si, K = K))
  if (!is.null(fit$error) && !is.na(fit$error)) {
    stop("EpiLPS failed on the observed series: ", fit$error)
  }
  fit$rho
}

# How much noisier NegBin(rho) is than Poisson at a typical count.
nbinom_noise_line <- function(rho, mu = 1000) {
  sd_nb <- sqrt(mu + mu^2 / rho)
  cat(sprintf("  at mu = %.0f: SD %.0f (NegBin) vs %.0f (Poisson), %.1fx the noise\n",
              mu, sd_nb, sqrt(mu), sd_nb / sqrt(mu)))
  invisible(NULL)
}

# ------------------------------------------------------------------------------
# 3. Replicates
# ------------------------------------------------------------------------------

# Simulate n_sim replicate epidemics from R_true under the TRUE SI. Simulated in
# the PARENT process so the RNG stream is reproducible; only the fits fork.
# Diverged replicates are dropped and counted. The MALA study pins RNGkind()
# before calling this, because its workers use L'Ecuyer-CMRG and this stream
# must not.
#
# With zero_run_threshold, a replicate with a run of MORE than that many
# zero-case days went extinct and is dropped too - the reference study's
# find_cols_with_consecutive_zeros() rule. It is applied after every replicate
# is drawn, so it never shifts the RNG stream: the kept series are the same
# draws whatever the threshold.
simulate_replicates <- function(R_true, si_true, seed_incidence, n_sim,
                                obs_model = c("poisson", "nbinom"),
                                overdispersion = NULL, seed,
                                zero_run_threshold = NULL) {
  obs_model <- match.arg(obs_model)
  n_days <- length(R_true)

  set.seed(seed)
  series <- lapply(seq_len(n_sim), function(r) {
    simulate_renewal_incidence(seed_incidence, R_true, si_true, n_days,
                               obs_model = obs_model,
                               overdispersion = overdispersion)
  })

  kept <- !vapply(series, is.null, logical(1))
  series <- series[kept]

  n_extinct <- 0L
  if (!is.null(zero_run_threshold)) {
    extinct <- vapply(series, function(x) {
      runs <- rle(x == 0)
      any(runs$lengths[runs$values] > zero_run_threshold)
    }, logical(1))
    n_extinct <- sum(extinct)
    series <- series[!extinct]
  }

  if (length(series) == 0L) {
    stop("Every replicate diverged or went extinct; check the plug-in true Rt path.")
  }

  if (is.null(zero_run_threshold)) {
    cat(sprintf("\nsimulated %d / %d replicates (%d diverged)\n",
                length(series), n_sim, sum(!kept)))
  } else {
    cat(sprintf("\nsimulated %d / %d replicates (%d diverged, %d extinct: a run of more than %d zero days)\n",
                length(series), n_sim, sum(!kept), n_extinct, zero_run_threshold))
  }
  cat(sprintf("median final-day incidence %.0f\n",
              stats::median(vapply(series, function(x) x[n_days], numeric(1)))))
  series
}

# ------------------------------------------------------------------------------
# 4. The SI settings to fit
# ------------------------------------------------------------------------------

# Two dimensions, crossed. Reported as two sweeps per family that share
# the correctly specified point (so each scenario is a complete curve); fitted
# as the unique (family, mean, sd) triples, so the shared point is only
# computed once per family.
#
# Three flags, because "correct" means three different things:
#   moments_correct - the assumed moments are the true ones, in any family
#   correct         - the fitting SI IS the data-generating SI. Only gamma_discr
#                     at the true moments qualifies, so this stays the control
#                     the reports and the selfcheck comparison select on.
#   family_ref      - the reference the other families are read against. It is
#                     gamma_bin, not gamma_discr: the non-gamma families are
#                     binned at +/-0.5 and gamma_discr is not, and that
#                     discretisation gap (TV 0.019) is large enough to swamp a
#                     family effect if it were left inside the comparison.
#
# Every setting's REALISED moments are recorded. Two reasons they can miss the
# request: discr_si drifts at the extremes (most visibly at mean 1.5, where the
# realised sd is nearer 2.9 than 3.4), and any binned family loses its left tail
# wherever it reaches below lag 0.5 - at requested mean 1.5 the realised mean
# comes out at 3.9 (gamma_bin and uniform), 3.2 (Weibull), 2.4 (lognormal).
# Settings that drift by more than drift_tol are FITTED AND REPORTED, flagged
# by drift_ok = FALSE; the figures draw them with open symbols.
build_si_misspec_grid <- function(fit_families, mean_si, sd_si,
                                  assumed_mean_grid, assumed_sd_grid,
                                  relative_error, max_si_lag, drift_tol) {
  unknown <- setdiff(fit_families, si_families)
  if (length(unknown) > 0L) {
    stop("Unknown SI family/families: ", paste(unknown, collapse = ", "),
         ". Known: ", paste(si_families, collapse = ", "))
  }

  settings <- do.call(rbind, lapply(fit_families, function(family) {
    s <- rbind(
      data.frame(scenario = "vary_sd",
                 assumed_mean = mean_si,
                 assumed_sd = assumed_sd_grid,
                 grid_value = assumed_sd_grid,
                 stringsAsFactors = FALSE),
      data.frame(scenario = "vary_mean",
                 assumed_mean = assumed_mean_grid,
                 assumed_sd = sd_si,
                 grid_value = assumed_mean_grid,
                 stringsAsFactors = FALSE)
    )
    s$rel_error <- rep(relative_error, 2L)
    s$family <- family
    s
  }))

  settings$moments_correct <- settings$assumed_mean == mean_si &
                              settings$assumed_sd == sd_si
  settings$correct    <- settings$family == "gamma_discr" & settings$moments_correct
  settings$family_ref <- settings$family == "gamma_bin"   & settings$moments_correct
  settings$key <- sprintf("%s_%.2f_%.2f", settings$family,
                          settings$assumed_mean, settings$assumed_sd)

  fit_keys <- unique(settings$key)
  fit_specs <- settings[match(fit_keys, settings$key),
                        c("key", "family", "assumed_mean", "assumed_sd")]

  si_list <- lapply(seq_len(nrow(fit_specs)), function(i) {
    make_si_family(fit_specs$family[i], fit_specs$assumed_mean[i],
                   fit_specs$assumed_sd[i], max_si_lag)
  })
  names(si_list) <- fit_specs$key

  realised <- t(vapply(si_list, si_moments, numeric(2)))
  fit_specs$realised_mean <- realised[, "mean"]
  fit_specs$realised_sd   <- realised[, "sd"]
  fit_specs$moment_drift  <- pmax(abs(fit_specs$realised_mean - fit_specs$assumed_mean),
                                  abs(fit_specs$realised_sd - fit_specs$assumed_sd))
  fit_specs$drift_ok      <- fit_specs$moment_drift <= drift_tol

  list(
    settings = settings,
    fit_specs = fit_specs,
    fit_keys = fit_keys,
    si_list = si_list,
    fit_families = fit_families,
    mean_si = mean_si,
    sd_si = sd_si,
    relative_error = relative_error,
    drift_tol = drift_tol
  )
}

# The size of the grid, and the settings whose realised moments drifted.
report_si_grid <- function(grid, n_series, fit_label = "fits") {
  n_settings <- nrow(grid$fit_specs)
  cat(sprintf("\n%d unique SI settings to fit (%d families x %d moment settings;\n",
              n_settings, length(grid$fit_families),
              n_settings / length(grid$fit_families)))
  cat(sprintf("  %d rows reported, the point the two sweeps share is fitted once per family)\n",
              nrow(grid$settings)))
  cat(sprintf("%d replicates x %d settings = %d %s\n",
              n_series, n_settings, n_series * n_settings, fit_label))

  drifted <- grid$fit_specs[!grid$fit_specs$drift_ok, ]
  if (nrow(drifted) > 0L) {
    cat(sprintf("\n%d setting(s) miss the requested moments by more than %.2f. They are\n",
                nrow(drifted), grid$drift_tol))
    cat("fitted and reported, but grid_value overstates how wrong the SI really is,\n")
    cat("so read them off realised_mean / realised_sd:\n")
    d <- drifted[order(drifted$family, drifted$assumed_mean, drifted$assumed_sd), ]
    print(format(d[, c("family", "assumed_mean", "assumed_sd",
                       "realised_mean", "realised_sd", "moment_drift")],
                 digits = 4), row.names = FALSE)
  }
  invisible(NULL)
}

# ------------------------------------------------------------------------------
# 5. Fit every replicate under every SI setting
# ------------------------------------------------------------------------------

# One worker per replicate; inside it, one fit per SI setting. EpiLPS's
# compiled backend was observed to corrupt memory across repeated calls in one
# session, so the parent process never calls it inside a loop. Each worker
# performs one fit per SI setting (17+ of them) rather than one fit total,
# which cuts the number of forks from n_sim * n_settings to n_sim; that was
# checked first (45 consecutive calls in a single process completed cleanly and
# returned bit-identical results on repeated settings).
#
# fit_fun(incidence, si, ...) is the estimator; fit_args are passed through.
# `seeds`, if given, is one integer per replicate: the worker calls
# set.seed(seeds[r], kind = rng_kind) before its first fit, so a stochastic
# estimator (MALA) is reproducible regardless of dispatch order.
#
# Returns one list per replicate, keyed by setting, after dropping replicates
# that failed at the fork level (mclapply reports those as a try-error rather
# than a list).
fit_si_grid <- function(series, si_list, fit_fun, n_cores, fit_args = list(),
                        seeds = NULL, rng_kind = "L'Ecuyer-CMRG",
                        packages = c("EpiEstim", "EpiLPS")) {
  # Force every argument now. The worker closure below captures this frame,
  # and on a PSOCK cluster the frame is serialised as-is: an argument still
  # held as an unevaluated promise (e.g. `grid$si_list`) would be evaluated on
  # the worker, in an environment where the caller's variables do not exist.
  force(series); force(si_list); force(fit_fun); force(fit_args)
  force(seeds); force(rng_kind)

  if (!is.null(seeds) && length(seeds) != length(series)) {
    stop("seeds must have one entry per replicate.")
  }

  cat("\nfitting")
  t_start <- proc.time()

  fits_by_rep <- run_parallel(
    seq_along(series),
    function(r) {
      if (!is.null(seeds)) set.seed(seeds[r], kind = rng_kind)
      incidence <- series[[r]]
      lapply(si_list, function(si) do.call(fit_fun, c(list(incidence, si), fit_args)))
    },
    n_cores = n_cores,
    packages = packages
  )

  elapsed <- (proc.time() - t_start)[["elapsed"]]
  cat(sprintf(" done in %.0f s (%.1f min)\n", elapsed, elapsed / 60))

  rep_ok <- vapply(fits_by_rep,
                   function(x) is.list(x) && length(x) == length(si_list),
                   logical(1))
  if (any(!rep_ok)) {
    cat(sprintf("%d replicate(s) failed at the fork level and were dropped\n",
                sum(!rep_ok)))
    fits_by_rep <- fits_by_rep[rep_ok]
  }
  if (length(fits_by_rep) == 0L) stop("Every replicate failed.")

  fits_by_rep
}

# ------------------------------------------------------------------------------
# 6. Score
# ------------------------------------------------------------------------------

# Score every setting (and every arm, if the estimator returns several) with
# score_interval_fits(). Returns the per-key metrics and by-day tables, keyed
# so expand_si_metrics() can join them back onto the reported sweeps.
# finite_days_only is passed to score_interval_fits().
score_si_grid <- function(fits_by_rep, grid, truth_scored, score_days,
                          arms = NULL,
                          intervals = list(ci = c("lower", "upper")),
                          extra = list(), finite_days_only = FALSE) {
  score_one <- function(key, arm) {
    fits <- lapply(fits_by_rep, function(rep) {
      f <- rep[[key]]
      if (is.null(arm) || (!is.null(f$error) && !is.na(f$error))) f else f[[arm]]
    })
    sc <- score_interval_fits(fits, truth_scored, score_days, intervals, extra,
                              finite_days_only = finite_days_only)
    if (is.null(sc)) return(NULL)

    id <- if (is.null(arm)) {
      data.frame(key = key, stringsAsFactors = FALSE)
    } else {
      data.frame(arm = arm, key = key, stringsAsFactors = FALSE)
    }
    list(metrics = cbind(id, sc$metrics), by_day = cbind(id, sc$by_day))
  }

  scored <- if (is.null(arms)) {
    lapply(grid$fit_keys, score_one, arm = NULL)
  } else {
    unlist(lapply(arms, function(a) lapply(grid$fit_keys, score_one, arm = a)),
           recursive = FALSE)
  }

  list(
    metrics_by_key = do.call(rbind, lapply(scored, function(x) x$metrics)),
    by_day_by_key  = do.call(rbind, lapply(scored, function(x) x$by_day))
  )
}

# Score every arm against its OWN truth. A w-day sliding-window estimator
# targets the trailing mean of Rt over its window, not the instantaneous value
# (trailing_mean() in simulation/truth.R), so scoring every arm against one
# truth would charge SI misspecification for a smoothing offset that is there
# at the correct SI too. `truth_by_arm` is a named list, one FULL-LENGTH truth
# vector per arm; each arm is scored by score_si_grid() on its own vector and
# the frames bound, which is exactly what one multi-arm call would return had
# the truth been shared. Used by the EpiEstim studies, whose arms are windows.
score_si_grid_per_arm <- function(fits_by_rep, grid, truth_by_arm, score_days, ...) {
  arms <- names(truth_by_arm)
  if (is.null(arms) || any(!nzchar(arms))) {
    stop("truth_by_arm must be a named list, one truth vector per arm.")
  }
  scored_arms <- lapply(arms, function(a) {
    score_si_grid(fits_by_rep, grid, truth_by_arm[[a]][score_days], score_days,
                  arms = a, ...)
  })
  list(
    metrics_by_key = do.call(rbind, lapply(scored_arms, function(s) s$metrics_by_key)),
    by_day_by_key  = do.call(rbind, lapply(scored_arms, function(s) s$by_day_by_key))
  )
}

# Expand the per-key scores back to the reported sweeps (each family's shared
# point appears in both), attach the realised moments, and order for writing.
expand_si_metrics <- function(grid, scored, by_arm = FALSE) {
  settings <- grid$settings
  fit_specs <- grid$fit_specs

  metrics <- merge(
    merge(settings,
          fit_specs[, c("key", "realised_mean", "realised_sd", "moment_drift", "drift_ok")],
          by = "key"),
    scored$metrics_by_key, by = "key"
  )
  by_day <- merge(settings[, c("key", "family", "scenario", "grid_value", "correct")],
                  scored$by_day_by_key, by = "key")

  if (isTRUE(by_arm)) {
    metrics <- metrics[order(metrics$arm, metrics$family, metrics$scenario,
                             metrics$grid_value), ]
    by_day <- by_day[order(by_day$arm, by_day$family, by_day$scenario,
                           by_day$grid_value, by_day$day), ]
  } else {
    metrics <- metrics[order(metrics$family, metrics$scenario, metrics$grid_value), ]
    by_day <- by_day[order(by_day$family, by_day$scenario, by_day$grid_value,
                           by_day$day), ]
  }

  list(metrics = metrics, by_day = by_day)
}

# The family comparison at the true moments: what does assuming the wrong
# SHAPE cost when the mean and sd are exactly right? Read against gamma_bin,
# the family that shares its discretisation with the other three. Also records
# where each family puts the most recent day's infectiousness weight
# (lag1_mass) - the mechanism the coverage differences should be explained by.
si_family_comparison <- function(metrics, grid, arms = NULL) {
  fit_families <- grid$fit_families
  si_list <- grid$si_list

  one_arm <- function(d) {
    d <- d[!duplicated(d$family), ]
    d <- d[match(fit_families, d$family), ]
    d <- d[!is.na(d$key), ]
    ref <- d[d$family == "gamma_bin", ]
    if (nrow(ref) == 1L) {
      d$dCoverage_vs_gamma_bin <- d$Coverage95 - ref$Coverage95
      d$dWidth_vs_gamma_bin    <- d$MeanCIWidth - ref$MeanCIWidth
      d$dRMSE_vs_gamma_bin     <- d$RMSE - ref$RMSE
    }
    d
  }

  fam_cmp <- metrics[metrics$moments_correct, ]
  fam_cmp <- if (is.null(arms)) {
    one_arm(fam_cmp)
  } else {
    do.call(rbind, lapply(arms, function(a) one_arm(fam_cmp[fam_cmp$arm == a, ])))
  }

  fam_cmp$lag1_mass <- vapply(fam_cmp$key, function(k) si_list[[k]][2L], numeric(1))
  fam_cmp$max_lag_used <- vapply(
    fam_cmp$key,
    function(k) as.numeric(max(which(si_list[[k]] > 1e-10)) - 1L),
    numeric(1)
  )
  fam_cmp
}

# ------------------------------------------------------------------------------
# 7. Reports
# ------------------------------------------------------------------------------

scenario_labels <- c(vary_sd = "vary_sd", vary_mean = "vary_mean")

held_moment_label <- function(sc, mean_si, sd_si) {
  if (sc == "vary_sd") sprintf("assumed mean held at the true %.1f", mean_si)
  else sprintf("assumed sd held at the true %.1f", sd_si)
}

# One table per (arm x) family x sweep. gamma_discr first, so the original
# study's two tables print exactly where they always did.
report_si_sweeps <- function(metrics, show_cols, mean_si, sd_si,
                             families, arms = NULL) {
  arm_set <- if (is.null(arms)) list(NULL) else as.list(arms)
  for (arm in arm_set) {
    for (fam in families) {
      for (sc in names(scenario_labels)) {
        d <- metrics[metrics$family == fam & metrics$scenario == sc, ]
        if (!is.null(arm)) d <- d[d$arm == arm, ]
        if (nrow(d) == 0L) next
        cat(sprintf("\n===== %s%s | %s (%s) =====\n",
                    if (is.null(arm)) "" else paste0(arm, " | "),
                    fam, sc, held_moment_label(sc, mean_si, sd_si)))
        out <- d[, show_cols]
        out$flag <- ifelse(d$correct, "  <- generating SI",
                    ifelse(d$family_ref, "  <- family reference",
                    ifelse(!d$drift_ok, "  <- moments drifted", "")))
        print(format(out, digits = 4), row.names = FALSE)
      }
    }
  }
  invisible(NULL)
}

report_family_comparison <- function(fam_cmp, fam_cols, mean_si, sd_si, arms = NULL) {
  cat("\n==============================================================\n")
  cat(sprintf("FAMILY COMPARISON at the true moments (mean %.1f, sd %.1f)\n",
              mean_si, sd_si))
  cat("==============================================================\n")
  if ("dCoverage_vs_gamma_bin" %in% names(fam_cmp)) {
    fam_cols <- c(fam_cols, "dCoverage_vs_gamma_bin", "dWidth_vs_gamma_bin")
  }
  if (is.null(arms)) {
    print(format(fam_cmp[, fam_cols], digits = 4), row.names = FALSE)
  } else {
    for (arm in arms) {
      d <- fam_cmp[fam_cmp$arm == arm, ]
      if (nrow(d) == 0L) next
      cat(sprintf("\n--- %s ---\n", arm))
      print(format(d[, fam_cols], digits = 4), row.names = FALSE)
    }
  }
  invisible(NULL)
}

# The two gamma rows differ ONLY in discretisation, so their gap sets the scale
# against which a family difference has to be judged.
report_discretisation_cost <- function(fam_cmp, arms = NULL) {
  cat("\nDiscretisation cost (gamma_discr - gamma_bin, same family, same moments):\n")
  arm_set <- if (is.null(arms)) list(NULL) else as.list(arms)
  for (arm in arm_set) {
    d <- if (is.null(arm)) fam_cmp else fam_cmp[fam_cmp$arm == arm, ]
    gd <- d[d$family == "gamma_discr", ]
    gb <- d[d$family == "gamma_bin", ]
    if (nrow(gd) == 1L && nrow(gb) == 1L) {
      cat(sprintf("  %scoverage %+.4f | CI width %+.4f | RMSE %+.4f\n",
                  if (is.null(arm)) "" else sprintf("%-8s ", arm),
                  gd$Coverage95 - gb$Coverage95,
                  gd$MeanCIWidth - gb$MeanCIWidth,
                  gd$RMSE - gb$RMSE))
    }
  }
  cat("  A family difference of this size or smaller is not distinguishable from\n")
  cat("  the choice of discretisation alone.\n")
  invisible(NULL)
}

# Which SI moment matters more, over the same relative error range - asked
# within each family, so a family that is more or less forgiving of a wrong
# moment shows as a different range rather than being averaged away.
report_sensitivity_ranges <- function(metrics, families, arms = NULL,
                                      cols = c(coverage = "Coverage95", bias = "Bias")) {
  cat(sprintf("\nSensitivity over the same %s relative error:\n",
              rel_error_range_label(metrics$rel_error)))
  arm_set <- if (is.null(arms)) list(NULL) else as.list(arms)
  for (arm in arm_set) {
    for (fam in families) {
      for (sc in names(scenario_labels)) {
        d <- metrics[metrics$family == fam & metrics$scenario == sc, ]
        if (!is.null(arm)) d <- d[d$arm == arm, ]
        if (nrow(d) == 0L) next
        parts <- vapply(names(cols), function(nm) {
          v <- d[[cols[[nm]]]]
          fmt <- if (nm == "bias") "%s %+.4f to %+.4f (range %.4f)" else "%s %.3f to %.3f (range %.3f)"
          sprintf(fmt, nm, min(v), max(v), diff(range(v)))
        }, character(1))
        cat(sprintf("  %s%-11s %-10s %s\n",
                    if (is.null(arm)) "" else sprintf("%-8s ", arm),
                    fam, sc, paste(parts, collapse = " | ")))
      }
    }
  }
  invisible(NULL)
}

# ------------------------------------------------------------------------------
# 8. Euler-Lotka reference
# ------------------------------------------------------------------------------

# The asymptotic bias the discrete Euler-Lotka equation predicts for each
# assumed SI, next to the bias measured. On each scored day the truth implies a
# growth rate under the TRUE SI (implied_growth_rates, si/euler_lotka.R);
# feeding that rate to each ASSUMED SI gives the R an estimator converges on
# once the start-of-series boundary is behind it, and the mean gap to the truth
# is the predicted bias. It is exact where the truth is flat and an
# approximation where it trends; for a constant-Rt scenario it is a true
# prediction. At the generating SI the prediction reproduces the truth, and
# `selfcheck` is the largest deviation from that. Restricted to one arm and one
# family: EpiEstim's own gamma_discr on the instantaneous (1-day) arm is the
# case the theory describes.
euler_lotka_reference <- function(metrics, grid, truth_scored, si_true,
                                  arm = "w1", family = "gamma_discr") {
  growth <- implied_growth_rates(truth_scored, si_true)
  selfcheck <- max(abs(euler_lotka_path(si_true, growth) - truth_scored))

  rows <- metrics[metrics$family == family, ]
  if (!is.null(arm) && "arm" %in% names(rows)) rows <- rows[rows$arm == arm, ]
  if (nrow(rows) == 0L) stop("No rows for family '", family, "' and arm '", arm, "'.")

  rows$R_euler_lotka <- vapply(rows$key, function(k) {
    mean(euler_lotka_path(grid$si_list[[k]], growth))
  }, numeric(1))
  rows$predicted_bias <- rows$R_euler_lotka - mean(truth_scored)
  rows$measured_bias <- rows$Bias
  rows$unexplained_bias <- rows$measured_bias - rows$predicted_bias
  rows <- rows[order(rows$scenario, rows$grid_value),
               c("scenario", "grid_value", "assumed_mean", "assumed_sd",
                 "realised_mean", "realised_sd", "R_euler_lotka",
                 "predicted_bias", "measured_bias", "unexplained_bias",
                 "Coverage95", "RMSE")]
  rownames(rows) <- NULL

  list(table = rows, growth_rates = growth, selfcheck = selfcheck)
}

report_euler_lotka <- function(el, mean_si, sd_si) {
  cat(sprintf("Self-check at the true SI: max |R_EL - truth| = %.2e over %d days\n",
              el$selfcheck, length(el$growth_rates)))
  cat(sprintf("Implied growth rate over the scored days: %.4f to %.4f per day\n",
              min(el$growth_rates), max(el$growth_rates)))
  for (sc in names(scenario_labels)) {
    d <- el$table[el$table$scenario == sc, ]
    if (nrow(d) == 0L) next
    cat(sprintf("\n--- %s (%s) ---\n", sc, held_moment_label(sc, mean_si, sd_si)))
    print(format(d[, c("grid_value", "predicted_bias", "measured_bias",
                       "unexplained_bias", "Coverage95")], digits = 4),
          row.names = FALSE)
  }
  cat("\n  predicted_bias is the Euler-Lotka asymptote averaged over the scored days;\n")
  cat("  unexplained_bias is what the sweep does to the estimator beyond that asymptote.\n")
  invisible(NULL)
}
