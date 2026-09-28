################################################################################
# Revised real-data analysis:
# Rt estimation and one-step-ahead forecasting for Queens, NY
# EpiEstim vs EpiLPS
#
# Main changes from the original script:
#   1. Dates are sorted before cumulative counts are differenced.
#   2. Both methods receive exactly the same discrete serial-interval distribution.
#   3. EpiEstim uses a 7-day sliding window in the primary analysis.
#   4. A fixed, pre-specified burn-in is used for both methods.
#   5. EpiLPS is run in fresh processes to avoid repeated-call instability.
#   6. Expanding-window forecasts use only information available at forecast time.
#   7. EpiEstim Rt draws use its Gamma posterior.
#   8. EpiLPS Rt draws use the lognormal approximation underlying estimR() output.
#   9. Forecast log scores are calculated from Poisson intensity draws, not from
#      already simulated incidence counts.
#  10. An SI-SD sensitivity analysis is included.
#  11. A semi-synthetic parametric-bootstrap study estimates the 95% credible
#      interval coverage rate of Rt itself.
#
# IMPORTANT INTERPRETATION:
# Two different 95% coverage rates are reported and they must not be confused.
#
#   (a) Coverage95_Rt_* (Section 6) is the coverage of the Rt credible interval.
#       The true Rt of the observed Queens epidemic is unknown, so it cannot be
#       measured on the observed series directly. Section 6 therefore runs a
#       semi-synthetic parametric bootstrap: the full-series Rt fit to the real
#       data is treated as a plug-in "true" trajectory, replicate epidemics are
#       simulated from it through the same Poisson renewal model, every method
#       arm is refitted on each replicate, and coverage is scored against that
#       known-by-construction Rt. It is genuine Rt coverage for an epidemic
#       calibrated to Queens, not for the historical Queens epidemic itself,
#       and it is conditional on the assumed plug-in truth (see
#       coverage_truth_source in R/core/config.R).
#
#       Coverage is reported against TWO targets because the methods do not
#       share an estimand: a w-day sliding window estimates the AVERAGE Rt over
#       [t - w + 1, t], whereas EpiLPS estimates the instantaneous Rt on day t.
#       Coverage95_Rt_ownEstimand is the number that tests calibration;
#       Coverage95_Rt_instantaneous and Coverage95_Rt_windowAvg are the same
#       intervals scored against both targets and are NOT comparable with each
#       other. A windowed arm undercovers the instantaneous truth by a
#       deterministic offset whenever the truth is trending, which is a target
#       mismatch and not a defect of the method.
#
#   (b) Coverage95_incidence (Section 7) is the coverage of the one-step-ahead
#       posterior predictive interval for the observed case count I[t+1]. It is
#       measured on the real series and says nothing about Rt calibration.
#
# PROJECT LAYOUT (see README.md for the full map):
#   main.R                 - this file: the analysis pipeline, run from r-proj/
#   R/core/config.R        - all user settings (data path, SI, methods, study settings)
#   R/core/                - init, paths, palette, worker processes
#   R/data/                - incidence loading and validation
#   R/si/                  - serial-interval constructors
#   R/estimators/          - EpiEstim, EpiLPS and EpiFilter fitting
#   R/simulation/          - renewal-model simulator, plug-in truths, scenarios
#   R/scoring/             - Rt coverage and forecast scoring
#   R/studies/             - the experiment harnesses (sensitivity, sim_study, si_misspec)
#   R/plots/               - all figures
#
# Run from r-proj/:  Rscript main.R
# Output: results/queens_rt/
################################################################################

rm(list = ls())

# ==============================================================================
# 0. Packages and modules
# ==============================================================================

source(file.path("R", "core", "init.R"))
require_packages(c("EpiEstim", "EpiLPS", "dplyr", "tidyr", "ggplot2", "purrr"))
source(file.path("R", "core", "config.R"))
source_project()

# ==============================================================================
# 1. Run setup
# ==============================================================================
# The seed is set here, immediately after the settings are read, which is
# where the original config.R set it: the forecast backtest (Section 7) draws
# from this stream after the coverage study has consumed whatever it consumes.

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(forecast_seed)
n_cores <- default_n_cores(n_cores_cap)

# ==============================================================================
# 2. Load and validate incidence data
# ==============================================================================

incidence_data <- load_incidence_data(
  data_file = data_file,
  state_abbr = state_abbr,
  county_name = county_name,
  burn_in_days = burn_in_days
)
n_obs <- nrow(incidence_data)

write.csv(
  incidence_data,
  file.path(output_dir, "queens_daily_incidence.csv"),
  row.names = FALSE
)

# ==============================================================================
# 3. Construct the shared SI distribution
# ==============================================================================

si_full <- make_si(mean_si, sd_si, max_si_lag)
si_epilps <- si_full[-1L]

si_lags <- seq.int(0L, max_si_lag)
realized_si <- si_moments(si_full)
realized_si_mean <- realized_si[["mean"]]
realized_si_sd <- realized_si[["sd"]]

si_table <- data.frame(
  lag = si_lags,
  probability = si_full
)
write.csv(si_table, file.path(output_dir, "serial_interval.csv"), row.names = FALSE)

cat("\n===== Data and SI summary =====\n")
cat("Location: ", county_name, ", ", state_abbr, "\n", sep = "")
cat("Date range: ", as.character(min(incidence_data$dates)), " to ",
    as.character(max(incidence_data$dates)), "\n", sep = "")
cat("Number of days: ", n_obs, "\n", sep = "")
cat("Requested SI mean/SD: ", mean_si, " / ", sd_si, "\n", sep = "")
cat("Realized truncated SI mean/SD: ",
    round(realized_si_mean, 3), " / ", round(realized_si_sd, 3), "\n", sep = "")

# ==============================================================================
# 4. Primary full-series Rt fits
# ==============================================================================

# EpiEstim primary fit: 7-day sliding window.
epiestim_full <- fit_epiestim_full(
  incidence_df = incidence_data,
  si_full = si_full,
  window_length = epiestim_window
)

# EpiLPS primary fit: run once in a fresh process.
epilps_primary_job <- list(
  job_id = "primary",
  incidence = incidence_data$I,
  si_epilps = si_epilps,
  K = K_epilps,
  dates = as.character(incidence_data$dates),
  return_type = "full"
)

epilps_full_raw <- run_epilps_jobs(
  jobs = list(epilps_primary_job),
  cores = n_cores
)[[1L]]

epilps_full <- tidy_epilps_result(epilps_full_raw)

if (all(is.na(epilps_full$R))) {
  stop("The primary EpiLPS fit failed: ", unique(epilps_full$converged)[1L])
}

# Fixed common evaluation start.
evaluation_start_index <- max(
  burn_in_days + 1L,
  max_si_lag + 1L,
  epiestim_window + 1L
)
evaluation_start_date <- incidence_data$dates[evaluation_start_index]

rt_primary <- bind_rows(
  epiestim_full %>%
    select(index, date, method, R, R_mean, R_sd, lower, upper),
  epilps_full %>%
    select(index, date, method, R, R_mean, R_sd, lower, upper)
) %>%
  filter(date >= evaluation_start_date) %>%
  arrange(method, date)

write.csv(
  rt_primary,
  file.path(output_dir, "primary_Rt_estimates.csv"),
  row.names = FALSE
)

cat("Common evaluation start: ", as.character(evaluation_start_date), "\n", sep = "")
cat("EpiLPS optimizer status: ", unique(epilps_full$converged)[1L], "\n", sep = "")

p_primary <- plot_primary_rt(
  rt_primary = rt_primary,
  mean_si = mean_si,
  sd_si = sd_si,
  evaluation_start_date = evaluation_start_date
)

if (interactive()) print(p_primary)
ggsave(
  filename = file.path(output_dir, "primary_Rt_comparison.png"),
  plot = p_primary,
  width = 10,
  height = 7,
  dpi = 300
)

# ==============================================================================
# 5. SI-SD sensitivity analysis
# ==============================================================================
# Varies only the assumed SI SD; the assumed mean is held at mean_si.

cat("\n===== Running SI-SD sensitivity analysis =====\n")

rt_sensitivity <- run_si_sensitivity(
  vary = "sd",
  grid = assumed_sd_grid,
  incidence_data = incidence_data,
  mean_si = mean_si,
  sd_si = sd_si,
  max_si_lag = max_si_lag,
  epiestim_window = epiestim_window,
  K_epilps = K_epilps,
  n_cores = n_cores,
  evaluation_start_date = evaluation_start_date
)

write.csv(
  rt_sensitivity,
  file.path(output_dir, "SI_SD_sensitivity_Rt_estimates.csv"),
  row.names = FALSE
)

p_sensitivity <- plot_si_sensitivity(
  sens_data = rt_sensitivity,
  vary = "sd",
  mean_si = mean_si,
  sd_si = sd_si
)

if (interactive()) print(p_sensitivity)
ggsave(
  filename = file.path(output_dir, "SI_SD_sensitivity_Rt.png"),
  plot = p_sensitivity,
  width = 10,
  height = 7,
  dpi = 300
)

# ==============================================================================
# 5B. SI-mean sensitivity analysis
# ==============================================================================
# Mirrors Section 5, but holds the assumed SD fixed at its true value (sd_si)
# and varies only the assumed mean.

cat("\n===== Running SI-mean sensitivity analysis =====\n")

rt_mean_sensitivity <- run_si_sensitivity(
  vary = "mean",
  grid = assumed_mean_grid,
  incidence_data = incidence_data,
  mean_si = mean_si,
  sd_si = sd_si,
  max_si_lag = max_si_lag,
  epiestim_window = epiestim_window,
  K_epilps = K_epilps,
  n_cores = n_cores,
  evaluation_start_date = evaluation_start_date
)

write.csv(
  rt_mean_sensitivity,
  file.path(output_dir, "SI_mean_sensitivity_Rt_estimates.csv"),
  row.names = FALSE
)

p_mean_sensitivity <- plot_si_sensitivity(
  sens_data = rt_mean_sensitivity,
  vary = "mean",
  mean_si = mean_si,
  sd_si = sd_si
)

if (interactive()) print(p_mean_sensitivity)
ggsave(
  filename = file.path(output_dir, "SI_mean_sensitivity_Rt.png"),
  plot = p_mean_sensitivity,
  width = 10,
  height = 7,
  dpi = 300
)

# ==============================================================================
# 6. Semi-synthetic Rt coverage study (95% credible interval coverage of Rt)
# ==============================================================================
# The true Rt of the observed Queens epidemic is unknown, so coverage of the Rt
# credible interval cannot be scored on the observed series. Instead:
#
#   1. Take the full-series fit from Section 4 as a plug-in "true" Rt path.
#   2. Simulate n_coverage_reps replicate epidemics from that path with the
#      same Poisson renewal model and the same SI used everywhere else,
#      seeding each replicate with the first max_si_lag observed days.
#   3. Refit EpiLPS and every EpiEstim arm in coverage_epiestim_windows on each
#      replicate, otherwise with the settings used in the primary analysis.
#   4. Score how often the 95% interval contains the known true Rt, starting
#      once the widest arm's estimation window has cleared the copied seed
#      incidence (see coverage_start_index in Section 6C).
#
# THE ESTIMAND MATTERS. A w-day sliding window estimates the AVERAGE Rt over
# [t - w + 1, t]; EpiLPS estimates the instantaneous Rt on day t. Scoring both
# against one target penalizes whichever arm's estimand is not that target, and
# the penalty is a deterministic offset, not a calibration failure. On this
# series the plug-in truth trends from about 0.62 to 1.51, so the 7-day average
# differs from the instantaneous value by about 0.037 on average, while the
# 7-day arm's interval half-width is only about 0.069 (Queens' counts are large,
# so the Poisson-only posterior SD is small). That offset alone pushes its
# coverage of the INSTANTANEOUS truth down to roughly 0.75 even though its
# coverage of its OWN estimand is essentially exactly 0.95. epi.R does not show
# this because its truth is constant at 1.0, where window-averaging introduces
# zero offset.
#
# Every arm is therefore scored against both targets, and the own-estimand
# column marks which one actually tests calibration. Coverage figures for
# different targets are not comparable with each other. The EpiEstim window = 1
# arm exists so that the instantaneous comparison against EpiLPS is like-for-
# like; expect much wider intervals there.
#
# The result is the coverage each method achieves on an epidemic calibrated to
# Queens. It is conditional on the plug-in truth: change coverage_truth_source
# and the numbers move, which is why that setting is exposed.
#
# The simulate / refit / truth-join machinery is sim_study.R's, shared with
# sim_main.R and truth_source_main.R: the plug-in truth is wrapped as a scenario
# and handed to the same functions. Each arm's own estimand is then joined per
# arm, which is also what makes a third window in coverage_epiestim_windows
# score correctly.

if (isTRUE(run_rt_coverage)) {
  cat("\n===== Running semi-synthetic Rt coverage study =====\n")

  if (!coverage_truth_source %in% c("epilps", "epiestim", "average")) {
    stop("coverage_truth_source must be 'epilps', 'epiestim', or 'average'.")
  }

  coverage_epiestim_windows <- sort(unique(as.integer(coverage_epiestim_windows)))
  if (length(coverage_epiestim_windows) == 0L ||
      any(is.na(coverage_epiestim_windows)) ||
      any(coverage_epiestim_windows < 1L)) {
    stop("coverage_epiestim_windows must be one or more positive integers.")
  }

  # ---------------------------------------------------------------------------
  # 6A. Plug-in true Rt trajectory
  # ---------------------------------------------------------------------------
  epilps_truth <- build_true_rt(epilps_full$index, epilps_full$R, n_obs)
  epiestim_truth <- build_true_rt(epiestim_full$index, epiestim_full$R, n_obs)

  true_Rt <- switch(
    coverage_truth_source,
    epilps = epilps_truth,
    epiestim = epiestim_truth,
    average = (epilps_truth + epiestim_truth) / 2
  )

  truth_label <- switch(
    coverage_truth_source,
    epilps = "EpiLPS",
    epiestim = paste0("EpiEstim (", epiestim_window, "-day window)"),
    average = "method-averaged"
  )

  if (any(!is.finite(true_Rt)) || any(true_Rt <= 0)) {
    stop("The plug-in true Rt trajectory contains non-positive or non-finite values.")
  }

  # ---------------------------------------------------------------------------
  # 6B. Simulate replicate epidemics
  # ---------------------------------------------------------------------------
  # Seeding with the observed history means every scored day (which starts at
  # evaluation_start_index > max_si_lag) is simulated, never copied.
  seed_incidence <- incidence_data$I[seq_len(max_si_lag)]

  coverage_scenario <- list(
    label = paste0("Plug-in truth: ", truth_label),
    build = function(time) true_Rt
  )

  replicates <- simulate_scenario_replicates(
    scenario = coverage_scenario,
    time = n_obs,
    n_sim = n_coverage_reps,
    I0 = seed_incidence,
    si_true = si_full,
    seed = coverage_seed
  )

  if (replicates$n_failed > 0L) {
    cat("Discarded ", replicates$n_failed, " diverging replicate(s).\n", sep = "")
  }

  # ---------------------------------------------------------------------------
  # 6C. Refit every arm on every replicate and join each arm's own estimand
  # ---------------------------------------------------------------------------
  # fit_scenario_grid() fits every EpiEstim window and EpiLPS (one fresh
  # process per replicate), joins true_R_instant / true_R_window / true_R_own,
  # and trims to the first index every arm can be scored at.
  #
  # A w-day window at day t reads incidence from [t - w + 1, t]. The first
  # max_si_lag days of every replicate are the copied observed series: identical
  # across replicates and NOT generated from the plug-in truth. An arm whose
  # window still overlaps that seed is therefore scored against a truth that did
  # not generate its data, and its coverage collapses for reasons unrelated to
  # calibration - measured at 0.06 for the 7-day arm and 0.42 for EpiLPS over
  # indices 61-66, against 0.95 and 0.92 once the window clears the seed.
  #
  # Scoring therefore starts once the widest arm's window lies entirely inside
  # the simulated region, and the same start is used for every arm so the arms
  # stay comparable. EpiLPS is a global smoother, so no finite warm-up fully
  # decontaminates it, but its dip is empirically confined to the same days.
  coverage_start_index <- evaluation_start_index + max(coverage_epiestim_windows) - 1L
  if (coverage_start_index >= n_obs) {
    stop("The coverage warm-up leaves no days to score; widen the series or narrow the windows.")
  }

  date_lookup <- data.frame(index = seq_len(n_obs), date = incidence_data$dates)

  rt_coverage_draws <- fit_scenario_grid(
    replicates = replicates,
    si_assumed = si_full,
    dates = incidence_data$dates,
    windows = coverage_epiestim_windows,
    reference_window = epiestim_window,
    K_epilps = K_epilps,
    n_cores = n_cores,
    run_epiestim = TRUE,
    run_epilps = TRUE
  ) %>%
    filter(index >= coverage_start_index) %>%
    left_join(date_lookup, by = "index") %>%
    select(replicate, index, date, method, estimand_window,
           R, R_sd, lower, upper, true_R_instant, true_R_window, true_R_own) %>%
    arrange(method, replicate, index)

  write.csv(
    rt_coverage_draws,
    file.path(output_dir, "rt_coverage_replicate_estimates.csv"),
    row.names = FALSE
  )

  # ---------------------------------------------------------------------------
  # 6D. Score coverage over the common evaluation window
  # ---------------------------------------------------------------------------
  rt_coverage_summary <- rt_coverage_draws %>%
    group_by(method) %>%
    group_modify(~rt_coverage_metrics(.x)) %>%
    ungroup() %>%
    mutate(truth_source = coverage_truth_source, .after = method)

  write.csv(
    rt_coverage_summary,
    file.path(output_dir, "rt_coverage_summary.csv"),
    row.names = FALSE
  )

  # Day-level coverage: the fraction of replicates whose interval covered the
  # target on that date, for both targets. Reveals where coverage fails
  # (typically where the truth trends fastest), which the overall rate hides.
  # sim_coverage_by_day() groups on scenario and si_scenario as well; there is
  # one of each here, so they are constants and dropped again afterwards.
  window_target_label <- paste0(epiestim_window, "-day average")

  rt_coverage_by_day <- rt_coverage_draws %>%
    mutate(scenario = "queens", si_scenario = "correct") %>%
    sim_coverage_by_day(reference_window = epiestim_window) %>%
    left_join(date_lookup, by = "index") %>%
    select(method, target, own_estimand, index, date,
           target_R, n_replicates, coverage, mean_ci_width) %>%
    arrange(method, target, index)

  write.csv(
    rt_coverage_by_day,
    file.path(output_dir, "rt_coverage_by_day.csv"),
    row.names = FALSE
  )

  cat("\n===== 95% credible interval coverage of Rt =====\n")
  cat("Plug-in truth: ", truth_label,
      " full-series fit; replicates requested: ", n_coverage_reps, "\n", sep = "")
  cat("Scoring starts ", as.character(incidence_data$dates[coverage_start_index]),
      " (index ", coverage_start_index, "), once the widest window clears the ",
      "copied seed.\n", sep = "")
  cat("The Estimand column names each arm's own target. Coverage against the\n",
      "other target reflects target mismatch, not calibration.\n", sep = "")
  print(as.data.frame(rt_coverage_summary), row.names = FALSE)

  p_coverage <- plot_coverage_by_day(
    rt_coverage_by_day = rt_coverage_by_day,
    truth_label = truth_label,
    epiestim_window = epiestim_window,
    window_target_label = window_target_label
  )

  if (interactive()) print(p_coverage)
  ggsave(
    filename = file.path(output_dir, "rt_coverage_by_day.png"),
    plot = p_coverage,
    width = 10,
    height = 8,
    dpi = 300
  )
}

# ==============================================================================
# 7. Expanding-window one-step-ahead forecast backtest
# ==============================================================================

if (isTRUE(run_backtest)) {
  cat("\n===== Running expanding-window one-step-ahead backtest =====\n")

  # Forecast origin t predicts incidence on day t + 1.
  first_forecast_origin <- max(
    burn_in_days,
    max_si_lag,
    epiestim_window + 1L
  )

  forecast_origins <- seq.int(
    from = first_forecast_origin,
    to = n_obs - 1L,
    by = forecast_stride
  )

  # ---------------------------------------------------------------------------
  # 7A. EpiEstim fits at each forecast origin
  # ---------------------------------------------------------------------------
  epiestim_origin_fits <- map_dfr(
    forecast_origins,
    function(t) {
      cbind(
        data.frame(t = t),
        fit_epiestim_last(
          incidence_df = incidence_data[seq_len(t), , drop = FALSE],
          si_full = si_full,
          window_length = epiestim_window
        )
      )
    }
  )

  # ---------------------------------------------------------------------------
  # 7B. EpiLPS fits at each forecast origin
  # ---------------------------------------------------------------------------
  epilps_origin_jobs <- map(
    forecast_origins,
    function(t) {
      list(
        job_id = as.character(t),
        incidence = incidence_data$I[seq_len(t)],
        si_epilps = si_epilps,
        K = K_epilps,
        dates = as.character(incidence_data$dates[seq_len(t)]),
        return_type = "last"
      )
    }
  )

  epilps_origin_fits <- bind_rows(
    run_epilps_jobs(epilps_origin_jobs, cores = n_cores)
  ) %>%
    transmute(
      t = as.integer(job_id),
      R_mean = R_mean,
      R_sd = R_sd,
      q025 = q025,
      q50 = q50,
      q975 = q975,
      status = converged
    ) %>%
    arrange(t)

  # ---------------------------------------------------------------------------
  # 7C. Construct posterior predictive forecasts
  # ---------------------------------------------------------------------------
  forecast_rows <- vector("list", length(forecast_origins) * 2L)
  row_counter <- 1L

  for (i in seq_along(forecast_origins)) {
    t <- forecast_origins[i]
    actual <- incidence_data$I[t + 1L]
    infectiousness <- calc_infectiousness(incidence_data$I, t, si_full)

    # EpiEstim Gamma posterior draws.
    ee_fit <- epiestim_origin_fits[epiestim_origin_fits$t == t, , drop = FALSE]
    ee_R_draws <- sample_epiestim_R(
      mean_R = ee_fit$R_mean[1L],
      sd_R = ee_fit$R_sd[1L],
      n_draws = n_posterior_draws
    )
    ee_summary <- make_forecast_summary(
      R_draws = ee_R_draws,
      infectiousness = infectiousness,
      actual = actual,
      n_draws = n_posterior_draws
    )

    forecast_rows[[row_counter]] <- data.frame(
      origin_index = t,
      origin_date = incidence_data$dates[t],
      target_date = incidence_data$dates[t + 1L],
      method = paste0("EpiEstim (", epiestim_window, "-day window)"),
      actual = actual,
      infectiousness = infectiousness,
      pred = ee_summary["pred"],
      lower = ee_summary["lower"],
      upper = ee_summary["upper"],
      log_score = ee_summary["log_score"],
      interval_score = ee_summary["interval_score"],
      fit_status = ee_fit$status[1L],
      row.names = NULL
    )
    row_counter <- row_counter + 1L

    # EpiLPS lognormal/Laplace approximation draws.
    el_fit <- epilps_origin_fits[epilps_origin_fits$t == t, , drop = FALSE]
    el_R_draws <- sample_epilps_R(
      q025 = el_fit$q025[1L],
      q50 = el_fit$q50[1L],
      q975 = el_fit$q975[1L],
      n_draws = n_posterior_draws
    )
    el_summary <- make_forecast_summary(
      R_draws = el_R_draws,
      infectiousness = infectiousness,
      actual = actual,
      n_draws = n_posterior_draws
    )

    forecast_rows[[row_counter]] <- data.frame(
      origin_index = t,
      origin_date = incidence_data$dates[t],
      target_date = incidence_data$dates[t + 1L],
      method = "EpiLPS",
      actual = actual,
      infectiousness = infectiousness,
      pred = el_summary["pred"],
      lower = el_summary["lower"],
      upper = el_summary["upper"],
      log_score = el_summary["log_score"],
      interval_score = el_summary["interval_score"],
      fit_status = el_fit$status[1L],
      row.names = NULL
    )
    row_counter <- row_counter + 1L
  }

  backtest <- bind_rows(forecast_rows) %>%
    mutate(
      across(c(actual, infectiousness, pred, lower, upper,
               log_score, interval_score), as.numeric)
    ) %>%
    arrange(method, target_date)

  write.csv(
    backtest,
    file.path(output_dir, "one_step_ahead_backtest.csv"),
    row.names = FALSE
  )

  performance <- backtest %>%
    group_by(method) %>%
    group_modify(~forecast_metrics(.x)) %>%
    ungroup()

  write.csv(
    performance,
    file.path(output_dir, "one_step_ahead_performance.csv"),
    row.names = FALSE
  )

  cat("\n===== One-step-ahead forecast performance =====\n")
  print(performance)

  p_forecast <- plot_forecasts(
    backtest = backtest,
    forecast_stride = forecast_stride
  )

  if (interactive()) print(p_forecast)
  ggsave(
    filename = file.path(output_dir, "one_step_ahead_forecasts.png"),
    plot = p_forecast,
    width = 10,
    height = 7,
    dpi = 300
  )
}

# ==============================================================================
# 8. Save session information
# ==============================================================================

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "sessionInfo.txt")
)

cat("\nAnalysis complete. Results were written to:\n")
cat(normalizePath(output_dir), "\n")
