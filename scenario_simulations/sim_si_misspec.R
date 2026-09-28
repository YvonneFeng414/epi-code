################################################################################
# EpiEstim under serial-interval misspecification, on synthetic epidemics.
#
# The port of the archive's three misspecification appendices (SI moments,
# SI family, family x moments) in one run, less their Queens target. Nine
# synthetic targets - every constant level and the six gradual pairs at the
# 7-day transition - at T = 100 and I0 = 2000. Per target the epidemics are
# simulated ONCE from the true SI and refit under every assumed SI on the
# harness's grid (R/studies/si_misspec.R): two 9-point moment sweeps at
# -80%..+80% relative error, crossed with five assumed families, so the design
# is paired and every difference across the grid is pure misspecification.
#
# Two arms, each scored against its own estimand: w1 (1-day window, the
# archive's nd = 1, primary) and w7 (the window main.R uses on the real data).
# The truth is known exactly here, which is what this study adds to the Queens
# one: the same sweeps, but no plug-in truth to argue about.
#
# For each target the Euler-Lotka asymptote is tabulated beside the measured
# bias (R/si/euler_lotka.R). For the constant targets the implied growth rate
# is exact, so the asymptote is a true prediction of where the bias should land
# once the start-of-series boundary is behind the estimator.
#
# THIS IS THE LONG RUN: 9 targets x n_sim x 85 SI settings x 2 windows. The
# cost projection after the first target says how long; `quick <- TRUE` runs
# two targets, one family and a handful of replicates in a couple of minutes.
#
# Run from r-proj/:  Rscript simulations/sim_si_misspec.R
# Output: results/simulations/si_misspec/
################################################################################

rm(list = ls())

# ==============================================================================
# 0. Packages and modules
# ==============================================================================

source(file.path("R", "core", "init.R"))
require_packages(c("EpiEstim", "ggplot2", "patchwork", "parallel"))
source(file.path("R", "core", "config.R"))          # mean_si, sd_si, the SI grids
source(file.path("simulations", "sim_config.R"))
source_project()

# ==============================================================================
# 1. Settings
# ==============================================================================

quick <- FALSE

n_sim <- if (quick) sim_ee_quick_n_sim else sim_ee_n_sim
targets <- if (quick) sim_ee_quick_targets else sim_ee_misspec_targets
fit_families <- if (quick) "gamma_discr" else si_families

time <- sim_ee_time
I0 <- sim_ee_I0
windows <- sim_ee_windows
arms <- paste0("w", windows)
max_si_lag <- sim_ee_max_si_lag

# A setting whose realised moments miss the requested ones by more than this is
# flagged, not dropped (see build_si_misspec_grid()).
drift_tol <- 0.25

output_dir <- sim_ee_results("si_misspec")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

n_cores <- default_n_cores(cap = 4L, reserve = 1L, floor = 2L)

# ==============================================================================
# 2. The true SI and the grid of assumed SIs
# ==============================================================================

si_true <- make_si(mean_si, sd_si, max_si_lag)
true_moments <- check_si_moments(si_true, mean_si)

cat("\n===== EpiEstim under SI misspecification: synthetic targets =====\n")
cat(sprintf("True SI: realised mean/sd %.4f / %.4f (requested %.1f / %.1f), lags 0-%d\n",
            true_moments[["mean"]], true_moments[["sd"]], mean_si, sd_si, max_si_lag))
cat(sprintf("T = %d, I0 = %d, %d replicates per target, windows %s\n",
            time, I0, n_sim, paste(windows, collapse = ", ")))
if (quick) cat("QUICK MODE: two targets, one family, reduced replicates\n")

# The grid does not depend on the target, so it is built once.
grid <- build_si_misspec_grid(fit_families, mean_si, sd_si,
                              assumed_mean_grid, assumed_sd_grid,
                              relative_error, max_si_lag, drift_tol)
report_si_grid(grid, n_sim,
               fit_label = sprintf("EpiEstim runs PER TARGET, each fitting %d windows", length(windows)))

plot_si_family_kernels(grid, si_family_colors, file.path(output_dir, "si_family_kernels.pdf"))

# ==============================================================================
# 3. Targets
# ==============================================================================

scenarios <- build_epiestim_scenarios(epiestim_scenario_table, time)
unknown <- setdiff(targets, names(scenarios))
if (length(unknown) > 0L) stop("Unknown target(s): ", paste(unknown, collapse = ", "))

# ==============================================================================
# 4. Simulate, fit, score, per target
# ==============================================================================

metrics_rows <- list()
by_day_rows <- list()
fam_rows <- list()
el_rows <- list()
t_start <- proc.time()

for (i in seq_along(targets)) {
  target <- sim_target(targets[i], scenarios, time, I0, windows)
  cat(sprintf("\n\n===== [%d/%d] %s =====", i, length(targets), target$label))

  series <- simulate_sim_target(target, si_true, n_sim, seed = sim_ee_seed + i)

  fits_by_rep <- fit_si_grid(series, grid$si_list, fit_epiestim_grid, n_cores,
                             fit_args = list(windows = windows), packages = "EpiEstim")

  scored <- score_si_grid_per_arm(fits_by_rep, grid, target$truth_by_arm, target$score_days)
  res <- expand_si_metrics(grid, scored, by_arm = TRUE)

  desc <- sim_descriptor(target, mean_si = mean_si, sd_si = sd_si)
  metrics_rows[[i]] <- bind_descriptor(desc, res$metrics)
  by_day_rows[[i]] <- bind_descriptor(desc[, c("subscenario", "T", "I0")], res$by_day)

  fam_cmp <- si_family_comparison(res$metrics, grid, arms = arms)
  fam_rows[[i]] <- bind_descriptor(desc[, c("subscenario", "major_scenario", "true_rt")], fam_cmp)

  el <- euler_lotka_reference(res$metrics, grid, target$truth_by_arm[["w1"]][target$score_days],
                              si_true, arm = "w1", family = "gamma_discr")
  el_tbl <- el$table
  el_tbl$growth_rate_mean <- mean(el$growth_rates)
  el_tbl$el_selfcheck <- el$selfcheck
  el_rows[[i]] <- bind_descriptor(desc[, c("subscenario", "major_scenario", "true_rt")], el_tbl)

  # Per target, the compact view: the coverage and width ranges of each sweep
  # for the generating family, both arms. The CSV carries every row.
  report_sensitivity_ranges(res$metrics, "gamma_discr", arms = arms,
                            cols = c(coverage = "Coverage95", width = "MeanCIWidth"))

  if (i == 1L) {
    report_cost_projection((proc.time() - t_start)[["elapsed"]], i, length(targets), "target")
  }
}

metrics <- do.call(rbind, metrics_rows)
by_day <- do.call(rbind, by_day_rows)
fam_cmp_all <- do.call(rbind, fam_rows)
el_all <- do.call(rbind, el_rows)

write.csv(metrics, file.path(output_dir, "si_misspec_metrics.csv"), row.names = FALSE)
write.csv(by_day, file.path(output_dir, "si_misspec_by_day.csv"), row.names = FALSE)
write.csv(fam_cmp_all, file.path(output_dir, "si_family_comparison.csv"), row.names = FALSE)
write.csv(el_all, file.path(output_dir, "euler_lotka_reference.csv"), row.names = FALSE)

# ==============================================================================
# 5. Report across targets
# ==============================================================================

cat("\n\n==============================================================\n")
cat("At the correct SI (control point), per target\n")
cat("==============================================================\n")
# The correct SI is the point the two sweeps share, so it has two rows per
# target and arm; one is enough here.
ctrl <- metrics[metrics$correct, c("subscenario", "arm", "Coverage95", "MCSE_Coverage",
                                   "MeanCIWidth", "Bias", "RMSE", "half_over_sd")]
ctrl <- ctrl[!duplicated(ctrl[, c("subscenario", "arm")]), ]
print(format(ctrl[order(ctrl$arm, ctrl$subscenario), ], digits = 4), row.names = FALSE)

cat("\n==============================================================\n")
cat("Euler-Lotka: predicted vs measured bias at the sweep ends (w1, gamma_discr)\n")
cat("==============================================================\n")
ends <- el_all[abs(el_all$grid_value - ifelse(el_all$scenario == "vary_sd", sd_si, mean_si)) >
                 0.7 * ifelse(el_all$scenario == "vary_sd", sd_si, mean_si), ]
print(format(ends[order(ends$subscenario, ends$scenario, ends$grid_value),
                  c("subscenario", "scenario", "grid_value", "predicted_bias",
                    "measured_bias", "unexplained_bias", "Coverage95")], digits = 4),
      row.names = FALSE)
cat(sprintf("\nEuler-Lotka self-check at the true SI, worst target: %.2e\n",
            max(el_all$el_selfcheck)))

# ==============================================================================
# 6. Figures
# ==============================================================================
# One page per target and arm for the moment sweeps (generating family), one
# page per target for the family x moment cross (w1). Multi-page PDFs.

subtitle_base <- sprintf(
  "Epidemics simulated with the true SI (mean = %.1f, sd = %.1f); only the SI assumed by EpiEstim is misspecified.  I0 = %d, T = %d, n_sim = %d.",
  mean_si, sd_si, I0, time, n_sim)

pdf(file.path(output_dir, "si_misspec_pages.pdf"), width = 12, height = 4.6)
for (tgt in targets) {
  for (arm in arms) {
    d <- metrics[metrics$subscenario == tgt & metrics$arm == arm &
                   metrics$family == "gamma_discr", ]
    if (nrow(d) == 0L) next
    print(plot_misspec_page(
      d, relative_error,
      title = sprintf("%s  |  %s window", scenarios[[tgt]]$label, window_label(d$window[1L])),
      subtitle = paste0(subtitle_base, "  Blue = assumed sd wrong, orange = assumed mean wrong; dotted = correct SI.")
    ))
  }
}
invisible(dev.off())

if (length(fit_families) > 1L) {
  pdf(file.path(output_dir, "si_family_pages.pdf"), width = 12, height = 8)
  for (tgt in targets) {
    d <- metrics[metrics$subscenario == tgt & metrics$arm == "w1", ]
    if (nrow(d) == 0L) next
    print(plot_family_page(
      d, relative_error, si_family_colors,
      title = sprintf("%s  |  1-day window  |  assumed SI family x moments", scenarios[[tgt]]$label),
      subtitle = paste0(subtitle_base, "  gamma_discr is the generating SI; open symbols = realised moments drifted > ",
                        drift_tol, ".")
    ))
  }
  invisible(dev.off())
}

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
cat(sprintf("\nWrote results to %s/\n", output_dir))
