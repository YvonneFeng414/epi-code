################################################################################
# Fully synthetic simulation study: EpiEstim vs EpiLPS
#
# The second kind of simulation in this project. main.R Section 6 runs a
# SEMI-SYNTHETIC study: the true Rt of the observed Queens epidemic is unknown,
# so a full-series fit to the real data is used as a plug-in truth and each
# replicate is seeded with observed incidence. Its coverage numbers are therefore
# conditional on that plug-in truth.
#
# This script runs a FULLY SYNTHETIC study instead. The true Rt path is written
# down analytically (scenarios.R), so it is known exactly and depends on no fit;
# the only seed is sim_I0 cases on day 1; no real data is read at all. That makes
# it the right place to ask questions the semi-synthetic study cannot answer
# cleanly, and it varies two things:
#
#   1. The SHAPE of the truth - constant, step, linear ramp, sine (rt_scenarios).
#   2. The SI handed to the ESTIMATORS, while epidemics are always generated with
#      the true SI (si_scenarios). This isolates SI misspecification.
#
# WHAT THIS FIXES RELATIVE TO epi.R. epi.R scored every window (nd = 1, 4, 7, 14)
# against the instantaneous truth R_t[valid_days]. That is correct only because
# its truth is constant at 1.0, where window averaging introduces no offset. A
# w-day sliding window estimates the AVERAGE Rt over [t - w + 1, t], so on any
# non-constant scenario scoring it against the instantaneous truth mixes a
# deterministic target-mismatch offset into what looks like a calibration result.
# Every arm here is scored against its own estimand (see build_truth_lookup() in
# sim_study.R), and against both targets for comparison.
#
# Run from r-proj/:  Rscript sim_main.R
# Output: results/simulation_study/
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
# 1. Settings
# ==============================================================================

if (!isTRUE(run_pure_simulation)) {
  stop("run_pure_simulation is FALSE in config.R; nothing to do.")
}

dir.create(sim_output_dir, showWarnings = FALSE, recursive = TRUE)
n_cores <- default_n_cores(n_cores_cap)

sim_windows <- sort(unique(as.integer(sim_windows)))
if (length(sim_windows) == 0L || any(is.na(sim_windows)) ||
    any(sim_windows < 1L)) {
  stop("sim_windows must be one or more positive integers.")
}

# Fail on a malformed scenario now, in seconds, rather than after the first
# expensive round of fits.
validate_rt_scenarios(rt_scenarios, sim_time)
validate_si_scenarios(si_scenarios, sim_epilps_si_scenarios)

start_index <- sim_score_start_index(sim_windows, sim_reference_window)
if (start_index >= sim_time) {
  stop("The warm-up leaves no days to score; raise sim_time or narrow sim_windows.")
}

# The true SI, used only to generate epidemics.
si_true <- make_si(sim_mean_si, sim_sd_si, sim_max_si_lag)

# One assumed SI per misspecification scenario, built the same way so the only
# difference between cells is the (mean, sd) pair.
si_assumed_list <- lapply(si_scenarios, function(values) {
  make_si(values[["mean"]], values[["sd"]], sim_max_si_lag)
})

# Synthetic calendar dates. The study is indexed by day, but fit_epiestim_full()
# and EpiLPS both expect a date column, so a fixed arbitrary origin is supplied
# and every join is keyed on `index` rather than on the date.
sim_dates <- as.Date("2025-01-01") + seq_len(sim_time) - 1L

cat("\n===== Fully synthetic simulation study =====\n")
cat("Rt scenarios:      ", paste(names(rt_scenarios), collapse = ", "), "\n", sep = "")
cat("SI scenarios:      ", paste(names(si_scenarios), collapse = ", "), "\n", sep = "")
cat("EpiLPS runs under: ", paste(sim_epilps_si_scenarios, collapse = ", "), "\n", sep = "")
cat("Days / replicates: ", sim_time, " / ", sim_n_sim, "\n", sep = "")
cat("EpiEstim windows:  ", paste(sim_windows, collapse = ", "), "\n", sep = "")
cat("Scoring starts at index ", start_index,
    " (fixed from the settings, not read off the fits).\n", sep = "")
cat("True SI mean/SD:   ", sim_mean_si, " / ", sim_sd_si, "\n", sep = "")

# ==============================================================================
# 2. Scenario truth table
# ==============================================================================
# Written and plotted before any fitting so the scenario shapes can be checked
# cheaply, and so the estimand offset is on record independently of any estimate.

scenario_truth <- scenario_truth_table(rt_scenarios, sim_time, sim_windows)

write.csv(
  scenario_truth,
  file.path(sim_output_dir, "scenario_truth.csv"),
  row.names = FALSE
)

p_truth <- plot_scenario_truth(scenario_truth)
ggsave(
  filename = file.path(sim_output_dir, "scenario_truth.png"),
  plot = p_truth,
  width = 10,
  height = 11,
  dpi = 300
)

# How far each window's estimand sits from the instantaneous truth over the
# scored range. This is the offset that makes a windowed arm undercover the
# instantaneous target even when it is perfectly calibrated for its own.
estimand_offset <- scenario_truth %>%
  filter(index >= start_index, estimand_window > 1L) %>%
  group_by(scenario, estimand_window) %>%
  summarise(
    mean_abs_offset = mean(abs(true_R_window - true_R_instant)),
    max_abs_offset = max(abs(true_R_window - true_R_instant)),
    .groups = "drop"
  )

cat("\n----- Mean |windowed estimand - instantaneous truth| over scored days -----\n")
print(as.data.frame(estimand_offset), row.names = FALSE, digits = 4)

# ==============================================================================
# 3. Simulate replicates, one set per Rt scenario
# ==============================================================================
# Epidemics depend only on the Rt scenario, so they are generated once here and
# reused across every assumed SI below.
#
# Each scenario gets its own seed, offset from sim_seed by its position, so the
# scenarios are not driven by the same random stream. The first scenario uses
# sim_seed itself, which keeps S1_constant aligned with epi.R's set.seed(123).

cat("\n===== Simulating replicates =====\n")

replicate_sets <- list()
for (i in seq_along(rt_scenarios)) {
  scenario_id <- names(rt_scenarios)[i]

  replicate_sets[[scenario_id]] <- simulate_scenario_replicates(
    scenario = rt_scenarios[[scenario_id]],
    time = sim_time,
    n_sim = sim_n_sim,
    I0 = sim_I0,
    si_true = si_true,
    seed = sim_seed + (i - 1L)
  )

  set_info <- replicate_sets[[scenario_id]]
  cat(scenario_id, ": kept ", length(set_info$series), " / ",
      set_info$n_requested, " replicates",
      if (set_info$n_failed > 0L) {
        paste0(" (", set_info$n_failed, " diverged and were discarded)")
      } else {
        ""
      },
      "; median final-day incidence ",
      stats::median(vapply(set_info$series, function(x) x[sim_time], numeric(1))),
      "\n", sep = "")
}

# ==============================================================================
# 4. Refit and score every (Rt scenario x assumed SI) cell
# ==============================================================================

cat("\n===== Fitting =====\n")

cell_results <- list()

for (scenario_id in names(rt_scenarios)) {
  for (si_name in names(si_scenarios)) {
    run_epilps <- si_name %in% sim_epilps_si_scenarios

    cat("  ", scenario_id, " x SI '", si_name, "'",
        if (run_epilps) " (EpiEstim + EpiLPS)" else " (EpiEstim only)",
        "\n", sep = "")

    draws <- fit_scenario_grid(
      replicates = replicate_sets[[scenario_id]],
      si_assumed = si_assumed_list[[si_name]],
      dates = sim_dates,
      windows = sim_windows,
      reference_window = sim_reference_window,
      K_epilps = K_epilps,
      n_cores = n_cores,
      run_epilps = run_epilps
    )

    cell_results[[paste(scenario_id, si_name, sep = "|")]] <- draws %>%
      mutate(scenario = scenario_id, si_scenario = si_name, .before = 1L)
  }
}

sim_draws <- bind_rows(cell_results) %>%
  mutate(si_scenario = factor(si_scenario, levels = names(si_scenarios))) %>%
  arrange(scenario, si_scenario, method, replicate, index) %>%
  mutate(si_scenario = as.character(si_scenario))

if (isTRUE(sim_write_replicate_estimates)) {
  write.csv(
    sim_draws,
    file.path(sim_output_dir, "sim_replicate_estimates.csv"),
    row.names = FALSE
  )
} else {
  cat("\nSkipping sim_replicate_estimates.csv (",
      format(nrow(sim_draws), big.mark = ","),
      " rows); set sim_write_replicate_estimates <- TRUE to write it.\n", sep = "")
}

# ==============================================================================
# 5. Score
# ==============================================================================
# rt_coverage_metrics() is reused unchanged from the semi-synthetic study.
# Grouping includes method, and method encodes the window, so each group has
# exactly one estimand_window as that function requires.

sim_coverage_summary <- sim_draws %>%
  group_by(scenario, si_scenario, method) %>%
  group_modify(~rt_coverage_metrics(.x)) %>%
  ungroup() %>%
  arrange(scenario, factor(si_scenario, levels = names(si_scenarios)), method)

write.csv(
  sim_coverage_summary,
  file.path(sim_output_dir, "sim_coverage_summary.csv"),
  row.names = FALSE
)

sim_by_day <- sim_coverage_by_day(sim_draws, sim_reference_window)

write.csv(
  sim_by_day,
  file.path(sim_output_dir, "sim_coverage_by_day.csv"),
  row.names = FALSE
)

cat("\n===== 95% credible interval coverage of Rt, correctly specified SI =====\n")
cat("Own-estimand coverage is the column that tests calibration. Coverage of the\n")
cat("other target reflects target mismatch, not a defect of the method.\n\n")
print(
  as.data.frame(
    sim_coverage_summary %>%
      filter(si_scenario == "correct") %>%
      select(scenario, method, Estimand, Coverage95_Rt_ownEstimand,
             MCSE_Coverage95_Rt_ownEstimand, Coverage95_Rt_instantaneous,
             Coverage95_Rt_windowAvg, MeanCIWidth, Bias_ownEstimand,
             RMSE_ownEstimand)
  ),
  row.names = FALSE,
  digits = 3
)

cat("\n===== Own-estimand coverage and bias by assumed SI =====\n")
print(
  as.data.frame(
    sim_coverage_summary %>%
      select(scenario, si_scenario, method,
             Coverage95_Rt_ownEstimand, Bias_ownEstimand)
  ),
  row.names = FALSE,
  digits = 3
)

# ==============================================================================
# 6. Figures
# ==============================================================================

# The Section 2 truth figure redrawn with every fitted arm on top - the EpiEstim
# windows and EpiLPS. Written from here rather than from Section 2 because it
# needs the fits. sim_draws itself is off by default
# (sim_write_replicate_estimates), so this collapsed table is the only per-day
# record of what the estimators produced.
truth_estimates <- sim_mean_estimate_by_day(sim_draws, si_scenario_shown = "correct")

write.csv(
  truth_estimates,
  file.path(sim_output_dir, "scenario_truth_estimates.csv"),
  row.names = FALSE
)

p_truth_est <- plot_scenario_truth(scenario_truth, estimates = truth_estimates)
ggsave(
  filename = file.path(sim_output_dir, "scenario_truth_with_estimates.png"),
  plot = p_truth_est,
  width = 10,
  height = 11,
  dpi = 300
)

p_by_day <- plot_sim_coverage_by_day(sim_by_day, si_scenario_shown = "correct")
ggsave(
  filename = file.path(sim_output_dir, "sim_coverage_by_day.png"),
  plot = p_by_day,
  width = 11,
  height = 10,
  dpi = 300
)

p_si <- plot_si_misspecification(sim_coverage_summary)
ggsave(
  filename = file.path(sim_output_dir, "sim_si_misspecification.png"),
  plot = p_si,
  width = 12,
  height = 8,
  dpi = 300
)

# ==============================================================================
# 7. Save session information
# ==============================================================================

capture.output(
  sessionInfo(),
  file = file.path(sim_output_dir, "sessionInfo.txt")
)

cat("\nSimulation study complete. Results were written to:\n")
cat(normalizePath(sim_output_dir), "\n")
