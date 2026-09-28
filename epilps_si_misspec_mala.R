################################################################################
# EpiLPS under serial-interval misspecification, on negative-binomial data,
# with the posterior explored by MALA instead of the Laplace/MAP approximation.
#
# This is the sampler-swapped twin of epilps_si_misspec.R. The experiment is
# unchanged - same truth (hole day removed), same I0 seed, same overdispersion,
# same two SI sweeps, same five families, same scoring - so that any difference
# in the reported curves is attributable to the inference engine and nothing
# else:
#
#   epilps_si_misspec.R      estimR()      LPSMAP  - Laplace approximation to
#                                                    the posterior of the
#                                                    B-spline coefficients,
#                                                    lognormal quantiles for Rt
#   epilps_si_misspec_mala.R estimRmcmc()  LPSMALA - Metropolis-adjusted
#                                                    Langevin sampler over the
#                                                    same posterior, empirical
#                                                    quantiles for Rt
#
# The design is PAIRED TWICE OVER:
#
#   * within this script, one set of simulated epidemics is generated once and
#     refit under every assumed SI, so differences across the grid are pure
#     misspecification effect;
#   * against epilps_si_misspec.R, the simulated epidemics are bit-identical.
#     That is why the data-generating overdispersion is still taken from a
#     LPSMAP fit of the observed series (section 4), and the simulator reads
#     the same config.R misspec_* settings and is seeded with the same sim_seed
#     under Mersenne-Twister: the MALA sampler is used for
#     ESTIMATION, which is the thing under study, and deliberately not for the
#     one calibration constant that would perturb the data. The MALA estimate
#     of the same constant is reported alongside as a diagnostic.
#
# What MALA buys that the Laplace fit cannot give:
#
#   * intervals that do not assume a lognormal posterior for Rt. Both the
#     equal-tailed 95% interval (scored as the primary, directly comparable to
#     the MAP script) and the 95% HPD interval (scored as a secondary) are
#     recorded, so interval SHAPE and interval CALIBRATION can be separated.
#   * an overdispersion estimate that is a posterior mean rather than a mode.
#
# COST. MALA is roughly 400x an estimR() call, and the family dimension
# multiplies the grid by the number of families. At the defaults below one fit of
# a ~250-day series takes ~25 s, and the grid is 100 replicates x 65 settings
# (5 families x 13 moment settings) = 6500 fits, i.e. ~45 CPU-hours, ~5.6 h wall
# clock on 8 cores. Section 4 times a real fit and section 6 prints the
# projection for the grid actually requested - READ IT before letting the run go.
#
# Two safe knobs if that is too slow: narrow fit_families (gamma_discr alone
# reproduces the original study), or cut n_sim, which costs Monte Carlo
# precision and nothing else - the MCSE column says how much. n_sim comes from
# config.R's misspec_n_sim, shared with the LPSMAP twin; cutting it here alone
# keeps the pairing (the series are drawn in order, so fewer replicates are the
# twin's first ones) but compares against the twin's full set. niter and burnin
# are NOT safe knobs; see the note in section 1.
#
# Run from r-proj/:  Rscript epilps_si_misspec_mala.R
# Output: results/epilps_si_misspec_mala_nohole/
################################################################################

rm(list = ls())

# ==============================================================================
# 0. Packages and modules
# ==============================================================================

source(file.path("R", "core", "init.R"))
require_packages(c("EpiLPS", "EpiEstim", "parallel"))
source(file.path("R", "core", "config.R"))   # mean_si, sd_si, misspec_*, K_epilps
source_project()

# ==============================================================================
# 1. Settings
# ==============================================================================
# Everything above the MALA block is the MAP script's. Changing any of it breaks
# the pairing with that script's simulated data.

paths <- project_paths()
truth_column <- "true_R_epilps"
common_start <- 16L

# config.R sets max_si_lag = 60; 30 is what the MAP script uses, and matching
# it is what makes the simulated series identical. At the true SI (7.5 / 3.4)
# lag 30 captures ~100% of the mass; the grid extremes (sd 5.44, mean 12.0)
# lose a little tail to truncation, identically in both scripts.
max_si_lag <- 30L

n_sim     <- misspec_n_sim
sim_seed  <- misspec_sim_seed

# --- MALA block -----------------------------------------------------------
# estimRmcmc()'s own defaults. The sampler is self-tuning (Langevin proposal
# with an adaptive step size), so there is nothing to hand-calibrate, but note
# that the package does NOT return the acceptance rate or the chain, so no
# per-fit convergence diagnostic beyond optimconverged is available.
#
# DO NOT treat these as a cheap runtime knob. The scored quantity is an
# interval width, and an under-burned chain reports a narrower one: at the
# correct SI, 400 iterations / 150 burn-in gave coverage 0.59 and width 0.119,
# against 0.85 and 0.191 at the values below on the same replicates. The whole
# coverage curve moves with them, so a run at reduced niter is not comparable
# to the LPSMAP twin or to another MALA run. Cut n_sim instead.
mcmc_niter  <- 5000L   # total MALA iterations per fit
mcmc_burnin <- 2000L   # discarded, leaving 3000 draws
mcmc_seed   <- 20260826L

# Reproducibility under forking. Each replicate's worker reseeds its own
# L'Ecuyer-CMRG stream from mcmc_seed + replicate index, so the result does not
# depend on how the workers happen to be scheduled. The simulator in section 5
# is deliberately left on Mersenne-Twister with sim_seed: that is what
# reproduces epilps_si_misspec.R's series exactly.
rep_seed <- function(r) mcmc_seed + 1000L * as.integer(r)

# NULL means estimate it from the observed series rather than assume a value.
overdispersion <- NULL

# The LPSMAP twin's output, read back if present so the control point and the
# curves can be compared directly rather than by eye across two terminals.
map_metrics_file <- file.path(paths$results_dir("epilps_si_misspec_nohole"),
                              "si_misspec_metrics.csv")

output_dir <- paths$results_dir("epilps_si_misspec_mala_nohole")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# The MAP script caps at 4 cores. A MALA fit costs ~300x a MAP fit, so the cap
# is raised here; two cores are left to the rest of the machine because each
# worker holds a (niter - burnin) x n_days draw matrix while it summarises.
n_cores <- default_n_cores(cap = 8L, reserve = 2L, floor = 2L)

# The assumed SI families, matching epilps_si_misspec.R so the two runs stay
# comparable cell for cell. NOTE ON COST: each family multiplies the number of
# LPSMALA fits. Section 6 prints the projection before any fitting starts.
fit_families <- si_families

# The family whose curves the LPSMALA-vs-LPSMAP figure compares. Both sides are
# filtered to it: those panels draw one point per rel_error, so leaving five
# families in would overplot five rows at every x.
map_compare_family <- "gamma_bin"

drift_tol <- 0.25

# ==============================================================================
# 2. The true SI, with the same self-test epilps_selfcheck.R uses
# ==============================================================================

si_true <- make_si(mean_si, sd_si, max_si_lag)
true_moments <- check_si_moments(si_true, mean_si)

cat("\n===== EpiLPS-MALA under SI misspecification (NegBin data) =====\n")
cat(sprintf("Sampler: LPSMALA, %d iterations, %d burn-in, %d retained draws per fit\n",
            mcmc_niter, mcmc_burnin, mcmc_niter - mcmc_burnin))
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
# The DATA-GENERATING value is the LPSMAP estimate, exactly as in
# epilps_si_misspec.R, because changing it would change the simulated series
# and destroy the pairing between the two scripts. The MALA estimate of the
# same quantity is computed alongside and reported: it is a posterior mean
# where the other is a mode, and the gap between them is the first piece of
# evidence about how much the Laplace approximation is smoothing.
#
# The MALA fit here doubles as the timing calibration for the grid.

if (is.null(overdispersion)) {
  overdispersion <- estimate_overdispersion(inp$obs$I, si_true, K_epilps)
}

cat("\ncalibrating on the observed series (one MALA fit, also times the grid) ...\n")
obs_fit_mala <- run_in_child(
  # Seed inside the child: the child's own seeding is not reproducible across
  # sessions, an explicit set.seed here is.
  function(incidence, si, K, niter, burnin, seed) {
    set.seed(seed, kind = "L'Ecuyer-CMRG")
    fit_epilps_mala(incidence, si, K, niter, burnin)
  },
  list(incidence = inp$obs$I, si = si_true, K = K_epilps,
       niter = mcmc_niter, burnin = mcmc_burnin, seed = mcmc_seed)
)
if (!is.null(obs_fit_mala$error) && !is.na(obs_fit_mala$error)) {
  stop("EpiLPS (LPSMALA) failed on the observed series: ", obs_fit_mala$error)
}

secs_per_fit <- obs_fit_mala$secs

cat("\nNegBin overdispersion on the observed series:\n")
cat(sprintf("  rho = %.2f  (LPSMAP mode)  <- used to generate the data\n", overdispersion))
cat(sprintf("  rho = %.2f  (LPSMALA posterior mean, diagnostic only)\n", obs_fit_mala$rho))
nbinom_noise_line(overdispersion)

# ==============================================================================
# 5. Simulate - negative binomial only
# ==============================================================================
# Byte-for-byte the MAP script's simulator, under the same RNG kind and seed.
# RNGkind is pinned explicitly rather than left at whatever the session default
# is, because the MALA workers use L'Ecuyer-CMRG and this stream must not.

RNGkind("Mersenne-Twister", "Inversion", "Rejection")
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
report_si_grid(grid, length(series), fit_label = "LPSMALA fits")

n_fits <- length(series) * nrow(grid$fit_specs)
cat(sprintf("projected cost at %.1f s/fit: %.1f CPU-hours, ~%.0f min wall clock on %d cores\n",
            secs_per_fit, n_fits * secs_per_fit / 3600,
            n_fits * secs_per_fit / 60 / n_cores, n_cores))

# ==============================================================================
# 7. Fit every replicate under every SI setting
# ==============================================================================
# Each worker seeds its own L'Ecuyer-CMRG stream from the replicate index, so
# the chains it then runs are reproducible regardless of dispatch order.

fits_by_rep <- fit_si_grid(
  series, grid$si_list, fit_epilps_mala, n_cores,
  fit_args = list(K = K_epilps, niter = mcmc_niter, burnin = mcmc_burnin),
  seeds = vapply(seq_along(series), rep_seed, integer(1))
)

# ==============================================================================
# 8. Score
# ==============================================================================
# Coverage95 / MeanCIWidth are the EQUAL-TAILED interval, which is what
# estimR() reports and therefore the like-for-like comparison with the MAP
# twin. CoverageHPD95 / MeanHPDWidth are the extra thing the sampler makes
# available: same posterior, shortest interval instead of the
# symmetric-in-probability one. Where they disagree, the posterior for Rt is
# skewed on that day.

scored <- score_si_grid(
  fits_by_rep, grid, inp$truth_scored, inp$score_days,
  intervals = list(ci = c("lower", "upper"), hpd = c("hpd_lower", "hpd_upper")),
  extra = list(
    rho_hat = function(fits) mean(vapply(fits, function(f) f$rho, numeric(1))),
    pct_converged = function(fits) 100 * mean(vapply(fits, function(f) isTRUE(f$converged), logical(1))),
    mean_fit_secs = function(fits) mean(vapply(fits, function(f) f$secs, numeric(1)))
  ),
  finite_days_only = TRUE
)
res <- expand_si_metrics(grid, scored)
metrics <- res$metrics
by_day <- res$by_day

write.csv(metrics, file.path(output_dir, "si_misspec_mala_metrics.csv"), row.names = FALSE)
write.csv(by_day, file.path(output_dir, "si_misspec_mala_by_day.csv"), row.names = FALSE)

# ==============================================================================
# 9. Report
# ==============================================================================

show_cols <- c("grid_value", "realised_mean", "realised_sd", "Coverage95",
               "MCSE_Coverage", "MeanCIWidth", "CoverageHPD95", "MeanHPDWidth",
               "Bias", "RMSE")
report_si_sweeps(metrics, show_cols, mean_si, sd_si, fit_families)

cat(sprintf("\nsampler health: %.1f%% of fits reported optimconverged, %.1f s mean fit time\n",
            mean(scored$metrics_by_key$pct_converged),
            mean(scored$metrics_by_key$mean_fit_secs)))

fam_cmp <- si_family_comparison(metrics, grid)
write.csv(fam_cmp, file.path(output_dir, "si_family_comparison.csv"), row.names = FALSE)

fam_cols <- c("family", "realised_mean", "realised_sd", "lag1_mass", "max_lag_used",
              "Coverage95", "MCSE_Coverage", "MeanCIWidth",
              "CoverageHPD95", "MeanHPDWidth", "Bias", "RMSE")
report_family_comparison(fam_cmp, fam_cols, mean_si, sd_si)
report_discretisation_cost(fam_cmp)

# --- the control point, against both reference numbers we have ---------------
ctrl <- metrics[metrics$correct, ][1L, ]
cat("\n==============================================================\n")
cat("CONTROL POINT (correctly specified SI)\n")
cat(sprintf("  LPSMALA, equal-tailed : %.4f (MCSE %.4f), width %.4f\n",
            ctrl$Coverage95, ctrl$MCSE_Coverage, ctrl$MeanCIWidth))
cat(sprintf("  LPSMALA, HPD          : %.4f (MCSE %.4f), width %.4f\n",
            ctrl$CoverageHPD95, ctrl$MCSE_CoverageHPD, ctrl$MeanHPDWidth))

# The LPSMAP twin's output. `correct` means the same thing in both the old and
# the new format (gamma_discr at the true moments), so the control-point line
# below works against either. The per-arm FIGURE needs the family column, which
# only a twin run since the family dimension was added will have - map_has_family
# gates it, so an older CSV degrades to no figure rather than to a silently
# overplotted one.
map_metrics <- NULL
map_has_family <- FALSE
if (file.exists(map_metrics_file)) {
  map_metrics <- read.csv(map_metrics_file, stringsAsFactors = FALSE)
  map_has_family <- "family" %in% names(map_metrics)
  map_ctrl <- map_metrics[map_metrics$correct %in% c(TRUE, "TRUE"), ][1L, ]
  cat(sprintf("  epilps_si_misspec.R   : %.4f (MCSE %.4f), width %.4f  [LPSMAP twin]\n",
              map_ctrl$Coverage95, map_ctrl$MCSE_Coverage, map_ctrl$MeanCIWidth))
  cat(sprintf("  difference vs twin     : %+.4f coverage, %+.4f width\n",
              ctrl$Coverage95 - map_ctrl$Coverage95,
              ctrl$MeanCIWidth - map_ctrl$MeanCIWidth))
  cat("  (same simulated data by construction, so this difference is the sampler)\n")
} else {
  cat(sprintf("  (%s not found - run epilps_si_misspec.R for the LPSMAP comparison)\n",
              map_metrics_file))
}
cat("==============================================================\n")

report_sensitivity_ranges(metrics, fit_families)

# Equal-tailed vs HPD: how much interval shape alone is worth.
cat("\nEqual-tailed vs HPD, across the whole grid:\n")
cat(sprintf("  coverage: HPD - equal-tailed ranges %+.4f to %+.4f (mean %+.4f)\n",
            min(metrics$CoverageHPD95 - metrics$Coverage95),
            max(metrics$CoverageHPD95 - metrics$Coverage95),
            mean(metrics$CoverageHPD95 - metrics$Coverage95)))
cat(sprintf("  width   : HPD is %.1f%% to %.1f%% of the equal-tailed width\n",
            100 * min(metrics$MeanHPDWidth / metrics$MeanCIWidth),
            100 * max(metrics$MeanHPDWidth / metrics$MeanCIWidth)))

# ==============================================================================
# 10. Figures
# ==============================================================================

run_caption <- sprintf("EpiLPS-MALA on NegBin data (rho = %.1f), %d replicates",
                       overdispersion, length(fits_by_rep))

plot_misspec_curves(
  metrics, file.path(output_dir, "si_misspec_mala_curves.png"),
  panels = misspec_panels_bias, scenario_colors = scenario_colors, hpd = TRUE,
  caption = sprintf("%s, %d draws/fit, gamma_discr SI; dotted = correct SI, dashed = nominal 0.95",
                    run_caption, mcmc_niter - mcmc_burnin),
  caption_cex = 0.75
)

plot_misspec_relative_error(
  metrics, file.path(output_dir, "si_misspec_mala_relative_error.png"),
  panels = misspec_panels_rmse, relative_error = misspec_relative_error,
  scenario_colors = scenario_colors
)

# LPSMALA vs LPSMAP on identical data. Only drawn when the twin has been run AND
# its CSV carries the family column: both sides are filtered to one family, since
# these panels plot a single point per rel_error and five families would draw
# five overlapping rows at every x.
if (!is.null(map_metrics) && !map_has_family) {
  cat(sprintf("\nSkipping the LPSMALA-vs-LPSMAP figure: %s predates the family\n",
              map_metrics_file))
  cat("dimension and has no `family` column. Re-run epilps_si_misspec.R to get it.\n")
}

# The twin may have been run under a narrower fit_families than this script, so
# the requested comparison family can be missing from one side. Falling through
# without checking would draw the figure with one of the two curves silently
# absent, which looks like a result rather than a missing input.
if (!is.null(map_metrics) && map_has_family) {
  shared_families <- intersect(unique(metrics$family), unique(map_metrics$family))
  if (!map_compare_family %in% shared_families) {
    fallback <- intersect(c("gamma_bin", "gamma_discr", fit_families), shared_families)
    if (length(fallback) == 0L) {
      cat("\nSkipping the LPSMALA-vs-LPSMAP figure: this run and the LPSMAP twin\n")
      cat("share no SI family, so there is nothing to compare like for like.\n")
      map_has_family <- FALSE
    } else {
      cat(sprintf("\nmap_compare_family '%s' is not in the LPSMAP twin's output;\n",
                  map_compare_family))
      cat(sprintf("falling back to '%s', which both runs have.\n", fallback[1L]))
      map_compare_family <- fallback[1L]
    }
  }
}

if (!is.null(map_metrics) && map_has_family) {
  png(file.path(output_dir, "si_misspec_mala_vs_map.png"),
      width = 1500, height = 950, res = 130)
  op <- par(mfrow = c(2, 3), mar = c(4.2, 4.4, 3, 1), oma = c(0, 0, 2, 0))

  cmp_cols <- c(mala = "#2a78d6", hpd = "#7fb2ec", map = "#333333")
  cmp_panels <- list(
    list(mala = "Coverage95", hpd = "CoverageHPD95", map = "Coverage95",
         lab = "95% CI coverage", main = "coverage", href = 0.95),
    list(mala = "MeanCIWidth", hpd = "MeanHPDWidth", map = "MeanCIWidth",
         lab = "mean CI width", main = "CI width", href = NA),
    list(mala = "RMSE", hpd = NA, map = "RMSE",
         lab = "RMSE of Rt", main = "RMSE", href = NA)
  )

  for (sc in sweep_scenarios) {
    d <- metrics[metrics$scenario == sc & metrics$family == map_compare_family, ]
    d <- d[order(d$rel_error), ]
    m <- map_metrics[map_metrics$scenario == sc &
                     map_metrics$family == map_compare_family, ]
    m <- m[order(m$rel_error), ]

    for (p in cmp_panels) {
      ys <- c(d[[p$mala]], m[[p$map]], p$href)
      if (!is.na(p$hpd)) ys <- c(ys, d[[p$hpd]])
      # Headroom so the legend sits in empty space rather than over a curve;
      # the panels differ in shape (falling, rising, U), so a fixed corner
      # would collide in at least one of them.
      ylim <- range(ys, na.rm = TRUE)
      ylim[2L] <- ylim[2L] + 0.28 * diff(ylim)
      plot(NA, xlim = range(misspec_relative_error) * 100, ylim = ylim,
           xlab = rel_error_xlab, ylab = p$lab, main = paste0(sc, ": ", p$main))
      if (!is.na(p$href)) abline(h = p$href, lty = 2, col = "grey40")
      abline(v = 0, lty = 3, col = "grey25")
      lines(m$rel_error * 100, m[[p$map]], type = "b", pch = 17, lty = 1,
            col = cmp_cols[["map"]])
      lines(d$rel_error * 100, d[[p$mala]], type = "b", pch = 19, lty = 1,
            col = cmp_cols[["mala"]])
      if (!is.na(p$hpd)) {
        lines(d$rel_error * 100, d[[p$hpd]], type = "b", pch = 1, lty = 2,
              col = cmp_cols[["hpd"]])
      }
      legend("top",
             legend = c("LPSMAP", "LPSMALA equal-tailed",
                        if (!is.na(p$hpd)) "LPSMALA HPD"),
             col = cmp_cols[c("map", "mala", if (!is.na(p$hpd)) "hpd")],
             pch = c(17, 19, if (!is.na(p$hpd)) 1),
             lty = c(1, 1, if (!is.na(p$hpd)) 2),
             bty = "o", bg = "white", box.col = NA, cex = 0.7)
    }
  }

  mtext(sprintf("LPSMALA vs LPSMAP on identical simulated data, %s SI - differences are the sampler alone",
                map_compare_family), outer = TRUE, cex = 0.8)
  par(op)
  invisible(dev.off())
}

plot_si_family_kernels(grid, si_family_colors,
                       file.path(output_dir, "si_family_kernels.png"))

plot_si_family_curves(
  metrics, file.path(output_dir, "si_family_curves.png"),
  panels = misspec_panels_rmse, grid = grid, family_colors = si_family_colors,
  caption = sprintf("%s; open symbols = realised moments drifted > %.2f",
                    run_caption, drift_tol)
)

cat(sprintf("\nWrote results to %s/\n", output_dir))
