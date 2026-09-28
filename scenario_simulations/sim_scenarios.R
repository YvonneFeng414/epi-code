################################################################################
# EpiEstim on fully synthetic epidemics: the Rt-shape x T x I0 study.
#
# The port of the archive's main study (Appendix 2), less its Queens scenario.
# Twenty-seven Rt shapes - three constant levels, six step pairs, the same six
# pairs as 7-, 14- and 20-day gradual transitions - at T = 100, plus the three
# constant shapes at T = 250, each crossed with three seed sizes. Every cell:
# n_sim Poisson epidemics from the true SI, EpiEstim under the CORRECT SI with a
# 1-day and a 7-day window, each window scored against its own estimand
# (trailing_mean(); the 1-day arm's is the instantaneous path).
#
# What this measures is EpiEstim's calibration when nothing is misspecified: how
# coverage, bias, RMSE and width depend on the shape of Rt, the length of the
# series and the number of cases. Everything the misspecification studies add
# is read against these rows.
#
# Differences from the archive, all from harmonising onto the shared harness:
# scoring is score_interval_fits() (per-replicate-first coverage with MCSE,
# RMSE, half_over_sd); scoring starts at day max(window) + 1 = 8 rather than
# day 2, so both arms are scored on the same days; the SI is make_si() at lag
# 30; seeds are one integer per cell. Same shapes, same grid, same n_sim.
#
# Run from r-proj/:  Rscript simulations/sim_scenarios.R
# Output: results/simulations/scenarios/
################################################################################

rm(list = ls())

# ==============================================================================
# 0. Packages and modules
# ==============================================================================

source(file.path("R", "core", "init.R"))
require_packages(c("EpiEstim", "dplyr", "purrr", "ggplot2", "patchwork", "parallel"))
source(file.path("R", "core", "config.R"))          # mean_si, sd_si
source(file.path("simulations", "sim_config.R"))    # the sim_ee_* settings
source_project()

# ==============================================================================
# 1. Settings
# ==============================================================================

# TRUE runs three shapes at one T and one I0 with a handful of replicates, in
# a minute or two, to check the pipeline end to end.
quick <- FALSE

n_sim <- if (quick) sim_ee_quick_n_sim else sim_ee_n_sim
times <- if (quick) sim_ee_time else c(sim_ee_time, sim_ee_time_long)
I0_grid <- if (quick) sim_ee_I0 else sim_ee_I0_grid

# Shapes that run at the long horizon (the archive's stepwise_time_values /
# gradual_time_values = 100 only).
long_shapes <- "constant"

windows <- sim_ee_windows
arms <- paste0("w", windows)
max_si_lag <- sim_ee_max_si_lag

output_dir <- sim_ee_results("scenarios")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

n_cores <- default_n_cores(cap = 4L, reserve = 1L, floor = 2L)

# ==============================================================================
# 2. The true SI
# ==============================================================================

si_true <- make_si(mean_si, sd_si, max_si_lag)
true_moments <- check_si_moments(si_true, mean_si)

cat("\n===== EpiEstim on synthetic epidemics: Rt shape x T x I0 =====\n")
cat(sprintf("True SI: realised mean/sd %.4f / %.4f (requested %.1f / %.1f), lags 0-%d\n",
            true_moments[["mean"]], true_moments[["sd"]], mean_si, sd_si, max_si_lag))
cat(sprintf("Windows %s; %d replicates per cell; seeds from %d\n",
            paste(windows, collapse = ", "), n_sim, sim_ee_seed))
if (quick) cat("QUICK MODE: reduced grid and replicates\n")

# ==============================================================================
# 3. Scenarios and the truth table
# ==============================================================================

scenarios_by_T <- lapply(times, function(t) build_epiestim_scenarios(epiestim_scenario_table, t))
names(scenarios_by_T) <- as.character(times)
for (t in times) validate_rt_scenarios(scenarios_by_T[[as.character(t)]], t)

table <- epiestim_scenario_table
if (quick) table <- table[table$subscenario %in% sim_ee_quick_scenarios, ]

truth <- scenario_truth_table(scenarios_by_T[[as.character(sim_ee_time)]][table$subscenario],
                              sim_ee_time, windows)
write.csv(truth, file.path(output_dir, "scenario_truth.csv"), row.names = FALSE)
ggsave(file.path(output_dir, "scenario_truth.pdf"),
       plot_scenario_shapes(truth, window = max(windows), ncol = 3L),
       width = 11, height = 2.4 * ceiling(nrow(table) / 3) + 1)

# ==============================================================================
# 4. The cells
# ==============================================================================

cells <- do.call(rbind, lapply(times, function(t) {
  shapes <- if (t == sim_ee_time) table else table[table$shape %in% long_shapes, ]
  expand.grid(T = t, I0 = I0_grid, subscenario = shapes$subscenario,
              stringsAsFactors = FALSE)
}))
cells <- cells[, c("T", "I0", "subscenario")]
rownames(cells) <- NULL

cat(sprintf("%d cells (%d shapes at T = %d, %d at T = %d, x %d seed sizes)\n",
            nrow(cells), nrow(table), sim_ee_time,
            sum(table$shape %in% long_shapes) * (length(times) > 1L), sim_ee_time_long,
            length(I0_grid)))

# ==============================================================================
# 5. Simulate, fit, score - one cell at a time
# ==============================================================================
# fit_si_grid() with a one-entry SI list: the correct SI only. Each replicate's
# fit carries both windows, and score_fits_per_arm() scores each against its
# own estimand.

metrics_rows <- list()
by_day_rows <- list()
t_start <- proc.time()

for (i in seq_len(nrow(cells))) {
  cell <- cells[i, ]
  target <- sim_target(cell$subscenario, scenarios_by_T[[as.character(cell$T)]],
                       cell$T, cell$I0, windows)

  cat(sprintf("\n[%d/%d] %s | T = %d | I0 = %d", i, nrow(cells),
              target$subscenario, target$T, target$I0))
  series <- simulate_sim_target(target, si_true, n_sim, seed = sim_ee_seed + i)

  fits_by_rep <- fit_si_grid(series, list(correct = si_true), fit_epiestim_grid,
                             n_cores, fit_args = list(windows = windows),
                             packages = "EpiEstim")
  fits <- lapply(fits_by_rep, function(rep) rep[["correct"]])

  scored <- score_fits_per_arm(fits, target)
  desc <- sim_descriptor(target, mean_si = mean_si, sd_si = sd_si)
  metrics_rows[[i]] <- bind_descriptor(desc, scored$metrics)
  by_day_rows[[i]] <- bind_descriptor(desc[, c("subscenario", "T", "I0")], scored$by_day)

  if (i == 1L) {
    report_cost_projection((proc.time() - t_start)[["elapsed"]], i, nrow(cells))
  }
}

metrics <- do.call(rbind, metrics_rows)
by_day <- do.call(rbind, by_day_rows)

write.csv(metrics, file.path(output_dir, "scenario_metrics.csv"), row.names = FALSE)
write.csv(by_day, file.path(output_dir, "scenario_by_day.csv"), row.names = FALSE)
write.csv(format_scenario_table(metrics),
          file.path(output_dir, "scenario_metrics_formatted.csv"), row.names = FALSE)

# ==============================================================================
# 6. Report
# ==============================================================================

cat("\n\n==============================================================\n")
cat("Coverage / bias / RMSE / width by cell (correct SI)\n")
cat("==============================================================\n")
for (arm in arms) {
  d <- metrics[metrics$arm == arm, ]
  cat(sprintf("\n--- %s (%s window) ---\n", arm, window_label(d$window[1L])))
  print(format(d[, c("subscenario", "T", "I0", "Coverage95", "MCSE_Coverage",
                     "Bias", "RMSE", "MeanCIWidth", "half_over_sd")], digits = 4),
        row.names = FALSE)
}

cat("\nhalf/sd is the interval half-width over the estimator's actual spread across\n")
cat("replicates; 1.96 is calibrated, above that is conservative.\n")

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
cat(sprintf("\nWrote results to %s/\n", output_dir))
