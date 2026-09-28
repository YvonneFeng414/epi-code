################################################################################
# EpiFilter under serial-interval misspecification, on Poisson data.
#
# The data are generated with rpois because that is EpiFilter's own observation
# model: epiFilter() scores observations with dpois(Iday, Lday * Rgrid). Giving
# it a correctly specified observation model leaves the serial interval as the
# ONLY misspecification, which is what this study is about.
#
# Two crossed dimensions of misspecification, shared with the EpiLPS studies
# (see R/studies/si_misspec.R):
#
# MOMENTS - two sweeps, on config.R's misspec grids, each holding one at its
# true value:
#   1. vary_sd   - assumed mean fixed at the true 7.5, assumed sd 1.36 ... 5.44
#   2. vary_mean - assumed sd fixed at the true 3.4, assumed mean 3.0 ... 12.0
#
# FAMILY - the shape of the assumed SI at those moments: gamma, lognormal,
# Weibull or uniform, each moment-matched to the requested (mean, sd). Matching
# both moments does not make two families the same distribution: at the true
# moments the lag-1 cell is 0.00023 (lognormal), 0.0041 (gamma), 0.0159
# (Weibull) and exactly 0 (uniform, supported only on lags 2..13).
#
# The SI enters EpiFilter only through the total infectiousness
# Lday[i] = sum_k I[i-k] w[k], following vignetteCOVID.R. Nothing else in the
# recursion sees it, so the family enters through that weighted sum alone.
#
# Two arms, which is the natural internal contrast:
#   * smoother - epiSmoother, retrospective, uses the whole series. Primary.
#   * filter   - epiFilter, causal / real-time. Free, since the smoother needs
#                the forward pass anyway.
#
# The truth is a plug-in path: EpiFilter's own smoothed fit of the real Queens
# series under the CORRECT SI. Real data has no known Rt, so some fit has to
# stand in for one, and taking EpiFilter's own is the self-consistent choice --
# it is shaped the way EpiFilter assumes, so nothing in the SI sweeps is
# attributable to a mismatch with some other method's idea of a smooth epidemic.
# Note the flip side: a truth drawn from a method's own fit can flatter it, so
# the level at any one grid point matters less than the SHAPE of each curve.
#
# THE SIMULATION DESIGN is that of Epiestim_simu/seed/
# run_sensitivity_queens_nohole.R, shared with the sibling studies through
# config.R's misspec_* block: the 2020-07-08 reporting hole is dropped from the
# truth (253 days), each replicate is seeded with I0 = 2000 cases on day 1
# only, extinct replicates (a run of more than 7 zero days) are dropped, and
# the grid is +/-60%. The reference scores from day 2; EpiFilter cannot, see
# filter_start below, so it scores from there.
#
# On hyperparameters: the diffusion noise eta is supplied by the user and never
# fitted, so its default could be doing work the SI sweeps would get credit for.
# It is fixed at Parag's default 0.1 for the sweeps, and section 10 sweeps it at
# the correct SI so the two contributions can be separated.
#
# Run from r-proj/:  Rscript epifilter_si_misspec.R
# Output: results/epifilter_si_misspec_nohole/
################################################################################

rm(list = ls())

# ==============================================================================
# 0. Packages and modules
# ==============================================================================
# EpiFilter is not on CRAN; epiFilter() and epiSmoother() are vendored verbatim
# in R/estimators/epifilter.R (Parag 2021, github.com/kpzoo/EpiFilter), so the
# recursions are the author's own. EpiEstim is used ONLY for discr_si().

source(file.path("R", "core", "init.R"))
require_packages(c("EpiEstim", "parallel"))
source(file.path("R", "core", "config.R"))   # mean_si, sd_si, the misspec_* design
source_project()

# ==============================================================================
# 1. Settings
# ==============================================================================

paths <- project_paths()

# First day of the common range the sibling studies score over, kept so they
# all simulate the same days. This study takes no truth column from the truth
# file (it builds its own below) but reads the file to locate the hole day,
# which also derives this value and checks it.
common_start <- 16L

# The TRUE serial interval (mean_si, sd_si from config.R) is used only to
# GENERATE epidemics. At the true SI (7.5 / 3.4) lag 30 captures ~100% of the
# mass. The widest assumed SI in the sweeps (sd 5.44) loses 0.4% here, which
# make_si() renormalises away.
max_si_lag <- 30L

n_sim    <- misspec_n_sim
sim_seed <- misspec_sim_seed

# EpiFilter grid and state model, following vignetteCOVID.R.
R_min <- 0.01
R_max <- 10
m_grid <- 1000L        # spacing 0.010; m = 2000 moves coverage by < 0.004
eta <- 0.1             # diffusion noise, Parag's default. NOT fitted.
alpha <- 0.025         # so Rhat rows 1 and 2 are the 2.5% and 97.5% quantiles

# EpiFilter cannot start at day 1. Lday[i] = sum_k I[i-k] w[k] is truncated at
# the start of the series, so on day 2 the Queens seed has I = 906 against
# Lday = 0.6, an implied Rt near 1550. That is far off a grid capped at 10, the
# Poisson likelihood underflows to zero on every grid point, and NaN propagates
# through the whole recursion. The vignettes handle this the same way, starting
# their series at day 8 (flu) and day 20 (SARS).
#
# The simulated replicates hit the same wall. They are seeded with a single
# I0 on day 1, so on the first days Lday is I0 times one early SI cell - and
# under a wrong SI that cell can be all but zero (the uniform family puts
# exactly 0 on lag 1; gamma_discr puts 1.3e-3 there at the true mean 7.5, but
# 4e-5 at mean 9 and 2e-9 at mean 12), so
# the implied Rt runs off the grid and the recursion fails the same way. The
# reference study scores from day 2; this one runs, and scores, from
# filter_start, the first day the filter can be trusted on every setting.
filter_start <- 21L

# eta values for the side analysis in section 10.
eta_grid <- c(0.05, 0.1, 0.2, 0.35, 0.5)

arms <- epifilter_arms
arm_lty <- list(smoother = 1, filter = 2)
arm_pch <- list(smoother = 19, filter = 1)

output_dir <- paths$results_dir("epifilter_si_misspec_nohole")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

n_cores <- default_n_cores(cap = 4L, reserve = 1L, floor = 2L)

# The assumed SI families to fit under. Narrow this to trade family coverage
# for runtime; "gamma_discr" alone reproduces the original study.
fit_families <- si_families

# A setting whose realised moments miss the requested ones by more than this is
# flagged, not dropped (see build_si_misspec_grid()).
drift_tol <- 0.25

# ==============================================================================
# 2. The true SI, with the same self-test the other scripts use
# ==============================================================================

si_true <- make_si(mean_si, sd_si, max_si_lag)
true_moments <- check_si_moments(si_true, mean_si)

grid_ef <- make_epifilter_grid(R_min, R_max, m_grid, alpha)

cat("\n===== EpiFilter under SI misspecification (Poisson data) =====\n")
cat(sprintf("True SI: realised mean/sd %.4f / %.4f (requested %.1f / %.1f), lags 0-%d\n",
            true_moments[["mean"]], true_moments[["sd"]], mean_si, sd_si, max_si_lag))
cat(sprintf("Grid: R in [%.2f, %.0f], m = %d (spacing %.4f), eta = %.2f\n",
            R_min, R_max, m_grid, (R_max - R_min) / (m_grid - 1L), eta))

# ==============================================================================
# 3. Inputs
# ==============================================================================

inp <- load_misspec_inputs(paths, I0 = misspec_I0, common_start = common_start,
                           truth_required = FALSE,
                           hole_threshold = misspec_hole_threshold)
observed <- inp$observed
n_days <- inp$n_days

# The seed day aside, every day is scored from day 2 - except that EpiFilter
# has no estimate before filter_start.
score_days <- inp$score_days[inp$score_days >= filter_start]

# The estimator, with everything but the data and the SI bound.
fit_fun <- run_epifilter
fit_args <- list(grid = grid_ef, eta = eta, filter_start = filter_start,
                 days = score_days)

# ==============================================================================
# 4. Truth path
# ==============================================================================
# Built here by smoothing the real observed series under the CORRECT SI, so
# EpiFilter is scored against a truth shaped like its own assumptions.

# The observed series still has the hole day in it, as the data the truth is
# fitted on should; the hole is dropped from the fitted path afterwards, the
# way the reference study drops it from EpiEstim's.
obs_td <- seq.int(filter_start, length(observed))
obs_fit <- run_epifilter(observed, si_true, grid_ef, eta, filter_start, days = obs_td)
if (!is.null(obs_fit$error) && !is.na(obs_fit$error)) {
  stop("EpiFilter failed on the observed series: ", obs_fit$error)
}

# The fit only covers filter_start onward. The leading days are held flat at
# the first fitted value. With a single seed day they are NOT idle: days 2 to
# filter_start - 1 drive the simulator's R[t]. EpiFilter never scores them.
R_true_full <- c(rep(obs_fit$smoother$R[1L], filter_start - 1L), obs_fit$smoother$R)
stopifnot(length(R_true_full) == length(observed))
R_true <- if (length(inp$hole_days) > 0L) R_true_full[-inp$hole_days] else R_true_full
stopifnot(length(R_true) == n_days)
truth_scored <- R_true[score_days]

report_hole_days(inp)
cat(sprintf("Truth: EpiFilter smoothed fit, %d days, Rt %.3f to %.3f (scored %.3f to %.3f)\n",
            length(R_true), min(R_true), max(R_true),
            min(truth_scored), max(truth_scored)))
cat(sprintf("Seed: %s\n", seed_label(inp)))
cat(sprintf("Filter runs from day %d; scoring days %d-%d (%d days), %d replicates\n",
            filter_start, min(score_days), n_days, length(score_days), n_sim))

# ==============================================================================
# 5. Simulate - Poisson, matching EpiFilter's own observation model
# ==============================================================================

series <- simulate_replicates(R_true, si_true, inp$seed_incidence, n_sim,
                              obs_model = "poisson", seed = sim_seed,
                              zero_run_threshold = misspec_zero_run_threshold)

# ==============================================================================
# 6. The SI settings
# ==============================================================================

grid <- build_si_misspec_grid(fit_families, mean_si, sd_si,
                              misspec_assumed_mean_grid, misspec_assumed_sd_grid,
                              misspec_relative_error, max_si_lag, drift_tol)
report_si_grid(grid, length(series), fit_label = "EpiFilter runs, each scored on both arms")

# ==============================================================================
# 7. Fit and score
# ==============================================================================

fits_by_rep <- fit_si_grid(series, grid$si_list, fit_fun, n_cores,
                           fit_args = fit_args, packages = "EpiEstim")

scored <- score_si_grid(fits_by_rep, grid, truth_scored, score_days, arms = arms,
                        finite_days_only = TRUE)
res <- expand_si_metrics(grid, scored, by_arm = TRUE)
metrics <- res$metrics
by_day <- res$by_day

write.csv(metrics, file.path(output_dir, "si_misspec_metrics.csv"), row.names = FALSE)
write.csv(by_day, file.path(output_dir, "si_misspec_by_day.csv"), row.names = FALSE)

# ==============================================================================
# 8. Report
# ==============================================================================

show_cols <- c("grid_value", "realised_mean", "realised_sd", "Coverage95",
               "MCSE_Coverage", "MeanCIWidth", "Bias", "RMSE", "half_over_sd")

# The per-setting tables, for the generating family only. Printing all five
# families x two arms x two sweeps would be twenty tables; the family effect is
# reported by the comparison block below instead, and the CSV carries every row.
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
# 9. Figures
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
  caption = sprintf("EpiFilter on Poisson data, %d replicates, eta = %.2f, gamma_discr SI; dotted = correct SI",
                    length(fits_by_rep), eta)
)

plot_si_family_kernels(grid, si_family_colors,
                       file.path(output_dir, "si_family_kernels.png"))

# The smoother is the primary arm, so the family figure uses it alone; carrying
# both arms as well as five families would be ten lines a panel. The filter's
# rows are in the CSV.
plot_si_family_curves(
  metrics, file.path(output_dir, "si_family_curves.png"),
  panels = misspec_panels_rmse, grid = grid, family_colors = si_family_colors,
  arm = "smoother",
  caption = sprintf("EpiFilter smoother on Poisson data, %d replicates, eta = %.2f; open symbols = realised moments drifted > %.2f",
                    length(fits_by_rep), eta, drift_tol)
)

# ==============================================================================
# 10. eta side analysis, at the correct SI
# ==============================================================================
# eta is EpiFilter's smoothing hyperparameter and is NOT fitted, so its default
# could be doing some of the work that the SI sweeps get credit for. This
# separates the two. `series` is the same set already simulated in section 5;
# the SI is correct here, so only eta changes.

cat("\n==============================================================\n")
cat("eta sensitivity at the CORRECT SI\n")
cat("==============================================================\n")

# The fit is built by a factory so its inputs travel with it: a worker on
# Windows does not share this script's global environment. The arguments are
# forced so the closure carries values, not promises pointing back here.
make_eta_fit <- function(si, grid, eta, filter_start, days) {
  force(si); force(grid); force(eta); force(filter_start); force(days)
  function(I) run_epifilter(I, si, grid, eta, filter_start, days = days)
}

eta_rows <- list()
for (e in eta_grid) {
  fits <- run_parallel(
    series,
    make_eta_fit(si_true, grid_ef, e, filter_start, score_days),
    n_cores = n_cores, packages = "EpiEstim"
  )
  for (arm in arms) {
    inner <- lapply(fits, function(f) if (!is.null(f$error) && !is.na(f$error)) f else f[[arm]])
    sc <- score_interval_fits(inner, truth_scored, score_days)
    if (is.null(sc)) next
    eta_rows[[length(eta_rows) + 1L]] <-
      cbind(data.frame(eta = e, arm = arm, stringsAsFactors = FALSE), sc$metrics)
  }
}
eta_tbl <- do.call(rbind, eta_rows)
write.csv(eta_tbl, file.path(output_dir, "eta_sensitivity.csv"), row.names = FALSE)
print(format(eta_tbl[, c("eta", "arm", "Coverage95", "MCSE_Coverage", "MeanCIWidth",
                         "Bias", "RMSE", "half_over_sd")], digits = 4), row.names = FALSE)

cat(sprintf("\nWrote results to %s/\n", output_dir))
