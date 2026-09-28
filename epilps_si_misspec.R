################################################################################
# EpiLPS under serial-interval misspecification, on negative-binomial data.
#
# A follow-on from epilps_selfcheck.R. That script established two things this
# one takes as given:
#
#   * Simulating Poisson understates the noise the Queens series actually shows.
#     EpiLPS estimates a NegBin overdispersion near 17 on the observed data, so
#     the negative binomial is the honest observation model. Only that arm is
#     kept here.
#   * Even with NegBin data and a CORRECTLY specified serial interval, EpiLPS's
#     95% intervals cover the plug-in truth only about 75% of the time.
#
# The question here is what the serial interval contributes to that. The data
# are always generated with the TRUE SI (a discr_si gamma, mean 7.5, sd 3.4);
# what varies is the SI handed to EpiLPS, along two crossed dimensions.
#
# MOMENTS - two sweeps, each holding one moment at its true value:
#
#   1. vary_sd   - assumed mean fixed at the true 7.5, assumed sd swept
#                  over 3.4 * (0.4 ... 1.6) = 1.36 ... 5.44
#   2. vary_mean - assumed sd fixed at the true 3.4, assumed mean swept
#                  over 7.5 * (0.4 ... 1.6) = 3.0 ... 12.0
#
# Both grids are config.R's misspec_assumed_sd_grid / misspec_assumed_mean_grid,
# i.e. the true value scaled by the same -60%..+60% relative errors, so the two
# sweeps are comparable step-for-step. Both pass through the correctly
# specified point, which is the control and is shared between them.
#
# THE SIMULATION DESIGN is that of Epiestim_simu/seed/
# run_sensitivity_queens_nohole.R, shared with the sibling studies through
# config.R's misspec_* block: each replicate is seeded with I0 = 2000 cases on
# day 1 only and scored from day 2; extinct replicates (a run of more than 7
# zero days) are dropped; and the 2020-07-08 reporting hole is removed from the
# truth. EpiLPS's own truth has no hole there (it is smooth through the day),
# but the same calendar day is dropped so every study runs on the same 253
# days.
#
# FAMILY - the shape of the assumed SI at those moments: gamma, lognormal,
# Weibull or uniform, each moment-matched to the requested (mean, sd). Two
# families can agree on both moments and still disagree sharply about where the
# mass sits: at the true moments the lag-1 cell is 0.00023 (lognormal), 0.0041
# (gamma), 0.0159 (Weibull) and exactly 0 (uniform, whose support starts at lag
# 2 and stops at 13). That cell is the most recent day's contribution to
# infectiousness, so this is not a cosmetic difference.
#
# The family dimension has FIVE levels, not four. gamma_discr is EpiEstim's
# shifted-gamma discretisation and is the data-generating SI; the other four are
# binned at +/-0.5, the scheme EpiLPS::Idist uses. Those two discretisations of
# the same gamma differ by TV 0.019, so gamma_bin is carried as the reference
# the non-gamma families are read against. Reading lognormal against
# gamma_discr instead would mix the family effect with a discretisation effect
# of comparable size. gamma_discr - gamma_bin is reported on its own as the
# discretisation cost.
#
# The design is PAIRED: one set of simulated epidemics is generated once and
# refit under every assumed SI. Differences across the grid are therefore pure
# misspecification effect, with no Monte Carlo noise from re-simulation.
#
# Caveat worth keeping in view: the plug-in truth comes from an EpiLPS fit of
# the real series under the correct SI, so it is smooth in exactly the way
# EpiLPS assumes. That flatters the correctly specified point specifically, and
# is why the interesting quantity is the SHAPE of each curve rather than the
# level at any one grid point.
#
# The experiment's steps live in R/studies/si_misspec.R, shared with the MALA
# and EpiFilter twins; this file is the settings and the EpiLPS-specific parts.
#
# Run from r-proj/:  Rscript epilps_si_misspec.R
# Output: results/epilps_si_misspec_nohole/
################################################################################

rm(list = ls())

# ==============================================================================
# 0. Packages and modules
# ==============================================================================
# EpiEstim is used ONLY for discr_si(), via make_si(). No estimation runs
# through it here.

source(file.path("R", "core", "init.R"))
require_packages(c("EpiLPS", "EpiEstim", "parallel"))
source(file.path("R", "core", "config.R"))   # mean_si, sd_si, misspec_*, K_epilps
source_project()

# ==============================================================================
# 1. Settings - the simulation design comes from config.R's misspec_* block
# ==============================================================================

paths <- project_paths()
truth_column <- "true_R_epilps"

# First day of the common range the pipeline selected. Derived from the truth
# file by load_misspec_inputs() and checked against this documented value.
common_start <- 16L

# The TRUE serial interval (mean_si, sd_si from config.R) is used only to
# GENERATE epidemics. config.R's max_si_lag is 60; 30 is kept here because it
# is what every misspec study and the MALA twin use, and matching it is what
# makes the simulated series identical across them. At the true SI
# (7.5 / 3.4) lag 30 already captures ~100% of the mass.
max_si_lag <- 30L

n_sim     <- misspec_n_sim
sim_seed  <- misspec_sim_seed

# NULL means estimate it from the observed series rather than assume a value.
overdispersion <- NULL

# epilps_selfcheck.R's NegBin coverage is no longer a reference point for the
# control: it seeds with 30 observed days and keeps the hole day, so it runs a
# different experiment.

output_dir <- paths$results_dir("epilps_si_misspec_nohole")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

n_cores <- default_n_cores(cap = 4L, reserve = 1L, floor = 2L)

# The assumed SI families to fit under. The epidemics are always generated from a
# discr_si gamma, so every family other than gamma_discr at the true moments is
# misspecified in SHAPE even when its mean and sd are exactly right. Narrow this
# vector to trade family coverage for runtime; "gamma_discr" alone reproduces
# the original study exactly.
fit_families <- si_families

# A setting whose realised moments miss the requested ones by more than this is
# flagged, not dropped (see build_si_misspec_grid()).
drift_tol <- 0.25

# ==============================================================================
# 2. The true SI, with the same self-test epilps_selfcheck.R uses
# ==============================================================================

si_true <- make_si(mean_si, sd_si, max_si_lag)
true_moments <- check_si_moments(si_true, mean_si)

cat("\n===== EpiLPS under SI misspecification (NegBin data) =====\n")
cat(sprintf("True SI: realised mean/sd %.4f / %.4f (requested %.1f / %.1f), lags 0-%d\n",
            true_moments[["mean"]], true_moments[["sd"]], mean_si, sd_si, max_si_lag))

# ==============================================================================
# 3. Inputs
# ==============================================================================

inp <- load_misspec_inputs(paths, truth_column, I0 = misspec_I0,
                           common_start = common_start,
                           hole_threshold = misspec_hole_threshold)
report_misspec_inputs(inp, n_sim)

# ==============================================================================
# 4. Overdispersion, estimated from the observed series
# ==============================================================================

if (is.null(overdispersion)) {
  overdispersion <- estimate_overdispersion(inp$obs$I, si_true, K_epilps)
}

cat(sprintf("\nNegBin overdispersion estimated on the observed series: rho = %.2f\n",
            overdispersion))
nbinom_noise_line(overdispersion)

# ==============================================================================
# 5. Simulate - negative binomial only
# ==============================================================================
# Renewal process driven by the TRUE SI. Day 1 is the I0 seed and is never
# scored; every later day is a NegBin draw.

series <- simulate_replicates(inp$R_true, si_true, inp$seed_incidence, n_sim,
                              obs_model = "nbinom",
                              overdispersion = overdispersion, seed = sim_seed,
                              zero_run_threshold = misspec_zero_run_threshold)

# ==============================================================================
# 6. The SI settings to fit
# ==============================================================================

grid <- build_si_misspec_grid(fit_families, mean_si, sd_si,
                              misspec_assumed_mean_grid, misspec_assumed_sd_grid,
                              misspec_relative_error, max_si_lag, drift_tol)
report_si_grid(grid, length(series), fit_label = "EpiLPS fits")

# ==============================================================================
# 7. Fit every replicate under every SI setting
# ==============================================================================

fits_by_rep <- fit_si_grid(series, grid$si_list, fit_epilps_map, n_cores,
                           fit_args = list(K = K_epilps))

# ==============================================================================
# 8. Score
# ==============================================================================
# Metric definitions are those of epilps_selfcheck.R. rho_hat is the mean
# fitted overdispersion across replicates.

scored <- score_si_grid(
  fits_by_rep, grid, inp$truth_scored, inp$score_days,
  extra = list(rho_hat = function(fits) mean(vapply(fits, function(f) f$rho, numeric(1)))),
  finite_days_only = TRUE
)
res <- expand_si_metrics(grid, scored)
metrics <- res$metrics
by_day <- res$by_day

write.csv(metrics, file.path(output_dir, "si_misspec_metrics.csv"), row.names = FALSE)
write.csv(by_day, file.path(output_dir, "si_misspec_by_day.csv"), row.names = FALSE)

# ==============================================================================
# 9. Report
# ==============================================================================

show_cols <- c("grid_value", "realised_mean", "realised_sd", "Coverage95",
               "MCSE_Coverage", "MeanCIWidth", "Bias", "RMSE")
report_si_sweeps(metrics, show_cols, mean_si, sd_si, fit_families)

fam_cmp <- si_family_comparison(metrics, grid)
write.csv(fam_cmp, file.path(output_dir, "si_family_comparison.csv"), row.names = FALSE)

fam_cols <- c("family", "realised_mean", "realised_sd", "lag1_mass", "max_lag_used",
              "Coverage95", "MCSE_Coverage", "MeanCIWidth", "Bias", "RMSE")
report_family_comparison(fam_cmp, fam_cols, mean_si, sd_si)
report_discretisation_cost(fam_cmp)

ctrl <- metrics[metrics$correct, ][1L, ]
cat("\n==============================================================\n")
cat("CONTROL POINT (correctly specified SI)\n")
cat(sprintf("  coverage %.4f (MCSE %.4f), width %.4f\n",
            ctrl$Coverage95, ctrl$MCSE_Coverage, ctrl$MeanCIWidth))
cat("==============================================================\n")

report_sensitivity_ranges(metrics, fit_families)

# ==============================================================================
# 10. Figures
# ==============================================================================

run_caption <- sprintf("EpiLPS on NegBin data (rho = %.1f), %d replicates",
                       overdispersion, length(fits_by_rep))

plot_misspec_curves(
  metrics, file.path(output_dir, "si_misspec_curves.png"),
  panels = misspec_panels_bias, scenario_colors = scenario_colors,
  caption = paste0(run_caption, ", gamma_discr SI; dotted = correct SI, dashed = nominal 0.95")
)

plot_misspec_relative_error(
  metrics, file.path(output_dir, "si_misspec_relative_error.png"),
  panels = misspec_panels_rmse, relative_error = misspec_relative_error,
  scenario_colors = scenario_colors
)

plot_si_family_kernels(grid, si_family_colors,
                       file.path(output_dir, "si_family_kernels.png"))

plot_si_family_curves(
  metrics, file.path(output_dir, "si_family_curves.png"),
  panels = misspec_panels_rmse, grid = grid, family_colors = si_family_colors,
  caption = sprintf("%s; open symbols = realised moments drifted > %.2f",
                    run_caption, drift_tol)
)

cat(sprintf("\nWrote results to %s/\n", output_dir))
