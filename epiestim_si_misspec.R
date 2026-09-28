################################################################################
# EpiEstim under serial-interval misspecification, on Poisson data.
#
# The fourth driver on the shared harness in R/studies/si_misspec.R, alongside
# the EpiLPS (MAP, MALA) and EpiFilter studies. Same simulation design (truth
# range, hole removal, I0 seed, replicate count and seed), same two crossed
# dimensions of misspecification, same scorer - so its rows can be read
# directly against theirs.
#
# THE SIMULATION DESIGN is that of Epiestim_simu/seed/
# run_sensitivity_queens_nohole.R, shared through config.R's misspec_* block:
#   * the 2020-07-08 reporting hole (EpiEstim truth Rt 0.0077, a day with zero
#     reported cases, not a collapse in transmission) is dropped from the
#     truth, leaving 253 days;
#   * each replicate is seeded with I0 = 2000 cases on day 1 only and every
#     later day is simulated and scored, from day 2;
#   * replicates with a run of more than 7 zero-case days are dropped as
#     extinct;
#   * the grid is +/-60%, the study grid without its +/-80% points.
# This driver's gamma_discr rows on the w1 arm are the reference study's
# experiment, at misspec_n_sim replicates rather than its 2000.
#
# The data are Poisson because that is EpiEstim's own observation model: the
# posterior for a window is Gamma-Poisson conjugate on the case counts. Giving
# it the observation model it assumes leaves the serial interval as the ONLY
# misspecification, which is what this study is about.
#
# MOMENTS - two sweeps, on config.R's misspec grids, each holding one at its
# true value:
#   1. vary_sd   - assumed mean fixed at the true 7.5, assumed sd 1.36 ... 5.44
#   2. vary_mean - assumed sd fixed at the true 3.4, assumed mean 3.0 ... 12.0
#
# FAMILY - the shape of the assumed SI at those moments: gamma_discr (the
# generating SI), gamma_bin, lognormal, Weibull or uniform, moment-matched.
#
# Two arms, the natural internal contrast for a sliding-window estimator:
#   * w1 - 1-day window, the instantaneous estimate. Primary: it is the fit that
#          generated this study's truth (see below), so at the correct SI it is
#          scored against its own assumptions.
#   * w7 - 7-day window, the setting main.R uses for the primary Queens fit.
#
# THE TRUTH IS PER ARM. A w-day EpiEstim window assumes Rt is constant over
# [t - w + 1, t] and estimates the AVERAGE Rt over that window, not the value
# on day t. Scoring w7 against the instantaneous path would therefore charge
# SI misspecification for a smoothing bias that is not misspecification at all
# and would be there at the correct SI. So w1 is scored against the plug-in
# path itself and w7 against its 7-day trailing mean (trailing_mean() in
# R/simulation/truth.R, the estimand of a windowed estimator). This is the one
# place the driver departs from its siblings, and it is why section 8 scores
# the arms separately and binds the results.
#
# The truth path is EpiEstim's own 1-day fit of the real Queens series under
# the CORRECT SI, from truth_source_main.R (plugin_truth_paths.csv, column
# true_R_epiestim) - the same self-consistent choice the EpiLPS studies make
# with true_R_epilps. Real data has no known Rt, so some fit has to stand in
# for one; taking the method's own means nothing in the SI sweeps is
# attributable to a mismatch with another method's idea of a smooth epidemic.
# The flip side holds too: a truth drawn from a method's own fit can flatter
# it, so read the SHAPE of each curve rather than the level at any one point.
#
# Section 11 adds a reference the other studies do not have: the asymptotic
# bias the discrete Euler-Lotka equation predicts for each assumed SI, so the
# measured bias curve can be read as "the theory" or "something else".
#
# Run from r-proj/:  Rscript epiestim_si_misspec.R
# Output: results/epiestim_si_misspec_nohole/
################################################################################

rm(list = ls())

# ==============================================================================
# 0. Packages and modules
# ==============================================================================
# Unlike the EpiLPS and EpiFilter drivers, EpiEstim is the ESTIMATOR here, not
# just the source of discr_si().

source(file.path("R", "core", "init.R"))
require_packages(c("EpiEstim", "parallel"))
source(file.path("R", "core", "config.R"))   # mean_si, sd_si, the misspec_* design
source_project()

# ==============================================================================
# 1. Settings - the simulation design comes from config.R's misspec_* block,
#    shared with the sibling studies
# ==============================================================================

paths <- project_paths()
truth_column <- "true_R_epiestim"

# First day of the common range the pipeline selected. Derived from the truth
# file by load_misspec_inputs() and checked against this documented value.
common_start <- 16L

# The TRUE serial interval (mean_si, sd_si from config.R) is used only to
# GENERATE epidemics. config.R's max_si_lag is 60; 30 is kept here because it
# is what every misspec study uses, and matching it is what makes the simulated
# series identical across them. At the true SI (7.5 / 3.4) lag 30 already
# captures ~100% of the mass.
max_si_lag <- 30L

n_sim    <- misspec_n_sim
sim_seed <- misspec_sim_seed

# The EpiEstim windows fitted as arms. w1 is primary (see header).
windows <- c(1L, 7L)
arms <- paste0("w", windows)
arm_lty <- list(w1 = 1, w7 = 2)
arm_pch <- list(w1 = 19, w7 = 1)

output_dir <- paths$results_dir("epiestim_si_misspec_nohole")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

n_cores <- default_n_cores(cap = 4L, reserve = 1L, floor = 2L)

# The assumed SI families to fit under. Narrow this to trade family coverage
# for runtime; "gamma_discr" alone reproduces the moment sweeps of the original
# study's design.
fit_families <- si_families

# A setting whose realised moments miss the requested ones by more than this is
# flagged, not dropped (see build_si_misspec_grid()).
drift_tol <- 0.25

# ==============================================================================
# 2. The true SI, with the same self-test the other scripts use
# ==============================================================================

si_true <- make_si(mean_si, sd_si, max_si_lag)
true_moments <- check_si_moments(si_true, mean_si)

cat("\n===== EpiEstim under SI misspecification (Poisson data) =====\n")
cat(sprintf("True SI: realised mean/sd %.4f / %.4f (requested %.1f / %.1f), lags 0-%d\n",
            true_moments[["mean"]], true_moments[["sd"]], mean_si, sd_si, max_si_lag))
cat(sprintf("Arms: %s (windows %s days); each arm scored against its own estimand\n",
            paste(arms, collapse = ", "), paste(windows, collapse = ", ")))

# ==============================================================================
# 3. Inputs
# ==============================================================================

inp <- load_misspec_inputs(paths, truth_column, I0 = misspec_I0,
                           common_start = common_start,
                           hole_threshold = misspec_hole_threshold)
report_misspec_inputs(inp, n_sim)

# The estimator, with everything but the data and the SI bound.
fit_fun <- fit_epiestim_grid
fit_args <- list(windows = windows)

# ==============================================================================
# 4. Observation model
# ==============================================================================
# Poisson, EpiEstim's own. No overdispersion is estimated: the NegBin studies
# calibrate rho from an EpiLPS fit because that is EpiLPS's likelihood; here a
# NegBin would ADD an observation-model misspecification on top of the SI one.

# ==============================================================================
# 5. Simulate
# ==============================================================================

series <- simulate_replicates(inp$R_true, si_true, inp$seed_incidence, n_sim,
                              obs_model = "poisson", seed = sim_seed,
                              zero_run_threshold = misspec_zero_run_threshold)

# ==============================================================================
# 6. The SI settings
# ==============================================================================

grid <- build_si_misspec_grid(fit_families, mean_si, sd_si,
                              misspec_assumed_mean_grid, misspec_assumed_sd_grid,
                              misspec_relative_error, max_si_lag, drift_tol)
report_si_grid(grid, length(series),
               fit_label = sprintf("EpiEstim runs, each fitting %d windows", length(windows)))

# ==============================================================================
# 7. Fit
# ==============================================================================
# EpiEstim is deterministic, so no per-replicate seeds are handed to the
# workers.

fits_by_rep <- fit_si_grid(series, grid$si_list, fit_fun, n_cores,
                           fit_args = fit_args, packages = "EpiEstim")

# ==============================================================================
# 8. Score - each arm against its own estimand
# ==============================================================================
# trailing_mean(x, 1) returns x, so the same construction gives the w1 arm the
# instantaneous path. score_si_grid_per_arm() scores each arm on its own
# vector and binds the frames - what one multi-arm call would give, had the
# truth been shared.
#
# Scored from day 2, EpiEstim has no estimate on the first days: it returns
# none for a window ending on or before the assumed SI's mean (days 2-7 at the
# true 7.5, days 2-11 at 12.0), and w7's first window ends on day 8. The scorer
# skips days without a finite estimate, as the reference study's
# complete.cases() does, so the early days scored differ along the mean sweep.
# finite_days_only extends that to the day-level summaries (half_over_sd).

truth_by_arm <- lapply(windows, function(w) trailing_mean(inp$R_true, w))
names(truth_by_arm) <- arms

scored <- score_si_grid_per_arm(fits_by_rep, grid, truth_by_arm, inp$score_days,
                                finite_days_only = TRUE)

res <- expand_si_metrics(grid, scored, by_arm = TRUE)
metrics <- res$metrics
by_day <- res$by_day

write.csv(metrics, file.path(output_dir, "si_misspec_metrics.csv"), row.names = FALSE)
write.csv(by_day, file.path(output_dir, "si_misspec_by_day.csv"), row.names = FALSE)

# ==============================================================================
# 9. Report
# ==============================================================================

show_cols <- c("grid_value", "realised_mean", "realised_sd", "Coverage95",
               "MCSE_Coverage", "MeanCIWidth", "Bias", "RMSE", "half_over_sd")

# The per-setting tables, for the generating family only. The family effect is
# reported by the comparison block below; the CSV carries every row.
report_si_sweeps(metrics, show_cols, mean_si, sd_si,
                 families = "gamma_discr", arms = arms)

fam_cmp <- si_family_comparison(metrics, grid, arms = arms)
write.csv(fam_cmp, file.path(output_dir, "si_family_comparison.csv"), row.names = FALSE)

fam_cols <- c("family", "realised_mean", "realised_sd", "lag1_mass", "max_lag_used",
              "Coverage95", "MCSE_Coverage", "MeanCIWidth", "Bias", "RMSE")
report_family_comparison(fam_cmp, fam_cols, mean_si, sd_si, arms = arms)
report_discretisation_cost(fam_cmp, arms = arms)

cat("\n==============================================================\n")
cat(sprintf("SUMMARY: sensitivity over the same %s relative error\n",
            rel_error_range_label(misspec_relative_error)))
cat("==============================================================\n")
report_sensitivity_ranges(metrics, fit_families, arms = arms,
                          cols = c(coverage = "Coverage95", width = "MeanCIWidth"))

cat("\nAt the correct SI (control point):\n")
for (arm in arms) {
  ctrl <- metrics[metrics$arm == arm & metrics$correct, ]
  if (nrow(ctrl) == 0L) next
  ctrl <- ctrl[1L, ]
  cat(sprintf("  %-8s coverage %.4f | CI width %.4f | half/sd %.2f\n",
              arm, ctrl$Coverage95, ctrl$MeanCIWidth, ctrl$half_over_sd))
}
cat("  (half/sd is the interval half-width over the estimator's actual spread\n")
cat("   across replicates; 1.96 is calibrated, above that is conservative)\n")

# ==============================================================================
# 10. Figures
# ==============================================================================

plot_misspec_relative_error(
  metrics, file.path(output_dir, "si_misspec_relative_error.png"),
  panels = misspec_panels_rmse, relative_error = misspec_relative_error,
  scenario_colors = scenario_colors,
  arms = arms, arm_lty = arm_lty, arm_pch = arm_pch
)

plot_misspec_curves(
  metrics, file.path(output_dir, "si_misspec_curves.png"),
  panels = misspec_panels_rmse, scenario_colors = scenario_colors,
  arms = arms, arm_lty = arm_lty, arm_pch = arm_pch, height = 900,
  caption = sprintf("EpiEstim on Poisson data, %d replicates, gamma_discr SI; solid = 1-day window, dashed = 7-day; dotted = correct SI",
                    length(fits_by_rep))
)

plot_si_family_kernels(grid, si_family_colors,
                       file.path(output_dir, "si_family_kernels.png"))

# w1 is the primary arm, so the family figure uses it alone; five families on
# two arms would be ten lines a panel. The w7 rows are in the CSV.
plot_si_family_curves(
  metrics, file.path(output_dir, "si_family_curves.png"),
  panels = misspec_panels_rmse, grid = grid, family_colors = si_family_colors,
  arm = "w1",
  caption = sprintf("EpiEstim 1-day window on Poisson data, %d replicates; open symbols = realised moments drifted > %.2f",
                    length(fits_by_rep), drift_tol)
)

# ==============================================================================
# 11. Euler-Lotka reference, at the generating family
# ==============================================================================
# If incidence grows at rate r and EpiEstim is told the SI is w, it converges on
# R_EL = 1 / sum_s w(s) exp(-r s) once the start-of-series boundary is behind
# it. On each scored day the truth path implies a growth rate under the TRUE SI
# (implied_growth_rates); feeding that rate to each ASSUMED SI predicts the R
# the w1 arm should report that day, and the mean gap to the truth is the
# predicted bias. It is the asymptotic result applied locally, so it is exact
# where the path is flat and approximate where it trends; the point of the
# table is whether the measured bias is that number or something else. At the
# correct SI the prediction reproduces the truth, which is checked first.

cat("\n==============================================================\n")
cat("Euler-Lotka reference: predicted vs measured bias (w1, gamma_discr)\n")
cat("==============================================================\n")

el <- euler_lotka_reference(metrics, grid, truth_by_arm[["w1"]][inp$score_days],
                            si_true, arm = "w1", family = "gamma_discr")
write.csv(el$table, file.path(output_dir, "euler_lotka_reference.csv"), row.names = FALSE)
report_euler_lotka(el, mean_si, sd_si)

cat(sprintf("\nWrote results to %s/\n", output_dir))
