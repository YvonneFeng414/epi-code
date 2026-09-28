################################################################################
# EpiEstim Rt trajectories on synthetic epidemics, under the correct SI and
# under a wrong one.
#
# The port of the archive's two trajectory figure sets, less their Queens page.
# The metrics tables report one number per series; these pages put the
# trajectory back - which days the bias appears on, which direction it goes,
# whether the credible interval still holds the truth.
#
# Five targets, one shape of each kind, all 1.2 -> 0.8: constant, step, and
# the 7-, 14- and 20-day gradual transitions. T = 100, I0 = 2000, 1-day window
# (the archive's nd = 1). Each target is simulated ONCE and refit under nine
# assumed SIs: the correct one and +/-40%, +/-80% relative error on the mean
# (sd correct) and on the sd (mean correct). The correct-SI fit gives the
# first page set; the rest give the second. The archive ran the misspecified
# fits as a separate 500-replicate study; here they share the epidemics with
# the correct fit, which is both cheaper and paired.
#
# Page set 1 (rt_trajectory_pages.pdf), one page per target: left, one
# representative epidemic (the one whose final size is nearest the median) with
# its posterior median and 95% CrI over the truth; right, the mean estimate
# over all epidemics with the mean CrI as a band and the 2.5 / 97.5 percentiles
# of the estimate across epidemics as dashed lines. A band narrower than the
# dashed lines is an overconfident posterior.
#
# Page set 2 (rt_misspec_trajectory_pages.pdf), one 2 x 2 page per target: the
# mean estimate per assumed SI over the truth (top) and per-day coverage
# (bottom), mean sweep left and sd sweep right.
#
# Run from r-proj/:  Rscript simulations/sim_rt_trajectories.R
# Output: results/simulations/rt_trajectories/
################################################################################

rm(list = ls())

# ==============================================================================
# 0. Packages and modules
# ==============================================================================

source(file.path("R", "core", "init.R"))
require_packages(c("EpiEstim", "ggplot2", "patchwork", "parallel"))
source(file.path("R", "core", "config.R"))
source(file.path("simulations", "sim_config.R"))
source_project()

# ==============================================================================
# 1. Settings
# ==============================================================================

quick <- FALSE

n_sim <- if (quick) sim_ee_quick_n_sim else sim_ee_n_sim
targets <- if (quick) sim_ee_quick_targets else sim_ee_traj_targets

time <- sim_ee_time
I0 <- sim_ee_I0
windows <- 1L                    # the trajectory figures are the 1-day arm
max_si_lag <- sim_ee_max_si_lag
errors <- sim_ee_traj_errors
drift_tol <- 0.25

output_dir <- sim_ee_results("rt_trajectories")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

n_cores <- default_n_cores(cap = 4L, reserve = 1L, floor = 2L)

# ==============================================================================
# 2. The true SI and the assumed SIs
# ==============================================================================
# build_si_misspec_grid() on the five-point error subset: nine unique settings
# (the correct point is shared by both sweeps), keyed and with realised moments,
# in the same vocabulary as the misspecification study.

si_true <- make_si(mean_si, sd_si, max_si_lag)
true_moments <- check_si_moments(si_true, mean_si)

grid <- build_si_misspec_grid("gamma_discr", mean_si, sd_si,
                              round(mean_si * (1 + errors), 2),
                              round(sd_si * (1 + errors), 2),
                              errors, max_si_lag, drift_tol)
correct_key <- grid$settings$key[grid$settings$correct][1L]
palette <- rel_error_palette(errors)

cat("\n===== EpiEstim Rt trajectories on synthetic epidemics =====\n")
cat(sprintf("True SI: realised mean/sd %.4f / %.4f (requested %.1f / %.1f), lags 0-%d\n",
            true_moments[["mean"]], true_moments[["sd"]], mean_si, sd_si, max_si_lag))
cat(sprintf("T = %d, I0 = %d, %d replicates per target, 1-day window\n", time, I0, n_sim))
report_si_grid(grid, n_sim, fit_label = "EpiEstim runs per target")
if (quick) cat("QUICK MODE: two targets, reduced replicates\n")

# ==============================================================================
# 3. Targets
# ==============================================================================

scenarios <- build_epiestim_scenarios(epiestim_scenario_table, time)
unknown <- setdiff(targets, names(scenarios))
if (length(unknown) > 0L) stop("Unknown target(s): ", paste(unknown, collapse = ", "))

# ==============================================================================
# 4. Simulate, fit, summarise, per target
# ==============================================================================

traj_rows <- list()
misspec_rows <- list()
pages_correct <- list()
pages_misspec <- list()
t_start <- proc.time()

subtitle_correct <- function(target, n_used) {
  paste0(
    sprintf("True SI mean = %.1f, sd = %.1f, correctly specified. I0 = %d, T = %d, n_sim = %d, 1-day window.\n",
            mean_si, sd_si, target$I0, target$T, n_used),
    "Left: one epidemic's posterior median and 95% CrI. Right: the mean estimate over all epidemics; the band is the mean CrI, ",
    "the dashed lines the 2.5 / 97.5 percentiles of the estimate across epidemics."
  )
}

subtitle_misspec <- function(target, n_used) {
  paste0(
    sprintf("Epidemics simulated with the true SI (mean = %.1f, sd = %.1f); only the SI assumed by EpiEstim is wrong. I0 = %d, T = %d, n_sim = %d, 1-day window.\n",
            mean_si, sd_si, target$I0, target$T, n_used),
    "Top: mean estimate per assumed SI over the truth (black). Bottom: fraction of epidemics whose 95% CrI holds the truth, per day; dashed = 0.95."
  )
}

for (i in seq_along(targets)) {
  target <- sim_target(targets[i], scenarios, time, I0, windows)
  cat(sprintf("\n\n===== [%d/%d] %s =====", i, length(targets), target$label))

  series <- simulate_sim_target(target, si_true, n_sim, seed = sim_ee_seed + i)
  fits_by_rep <- fit_si_grid(series, grid$si_list, fit_epiestim_grid, n_cores,
                             fit_args = list(windows = windows), packages = "EpiEstim")
  truth <- target$truth_by_arm[["w1"]]
  days <- target$score_days
  desc <- sim_descriptor(target, mean_si = mean_si, sd_si = sd_si)

  # --- correct SI: the representative epidemic and the Monte Carlo summary ---
  fits_correct <- arm_fits(lapply(fits_by_rep, function(rep) rep[[correct_key]]), "w1")
  rep_index <- pick_representative(series)
  single <- fit_trajectory(fits_correct[[rep_index]], truth, days)
  summary <- summarise_trajectories(fits_correct, truth, days)
  n_used <- max(summary$n_used)

  traj_rows[[i]] <- bind_descriptor(
    cbind(desc, data.frame(n_sim_used = n_used, representative_sim = rep_index)),
    summary
  )
  pages_correct[[i]] <- plot_trajectory_page(
    single, summary, n_used,
    title = target$label, subtitle = subtitle_correct(target, n_used)
  )

  # --- every assumed SI: per-setting summaries, expanded to the two sweeps ---
  per_key <- lapply(grid$fit_keys, function(k) {
    s <- summarise_trajectories(arm_fits(lapply(fits_by_rep, function(rep) rep[[k]]), "w1"),
                                truth, days)
    cbind(data.frame(key = k, stringsAsFactors = FALSE), s)
  })
  per_key <- do.call(rbind, per_key)
  settings <- grid$settings[, c("key", "scenario", "rel_error", "assumed_mean", "assumed_sd")]
  traj <- merge(settings, per_key, by = "key")
  traj$error_label <- rel_error_factor(traj$rel_error, errors)
  traj <- traj[order(traj$scenario, traj$rel_error, traj$day), ]

  misspec_rows[[i]] <- bind_descriptor(
    cbind(desc, data.frame(n_sim_used = n_used)),
    traj[, c("scenario", "rel_error", "assumed_mean", "assumed_sd", "day", "true_R",
             "mean_est", "mean_lower", "mean_upper", "mc_lower", "mc_upper",
             "coverage", "n_used", "error_label")]
  )
  pages_misspec[[i]] <- plot_misspec_trajectory_page(
    traj, palette,
    title = paste0(target$label, "  |  assumed SI wrong"),
    subtitle = subtitle_misspec(target, n_used)
  )

  # The scalar view of the same fits, for cross-reading against the tables.
  ctrl <- traj[traj$rel_error == 0 & traj$scenario == "vary_sd", ]
  cat(sprintf("\ncorrect SI: mean bias %+.4f, mean coverage %.4f over days %d-%d (%d replicates)\n",
              mean(ctrl$mean_est - ctrl$true_R, na.rm = TRUE),
              mean(ctrl$coverage, na.rm = TRUE), min(days), max(days), n_used))

  if (i == 1L) {
    report_cost_projection((proc.time() - t_start)[["elapsed"]], i, length(targets), "target")
  }
}

# ==============================================================================
# 5. Write
# ==============================================================================

traj_summary <- do.call(rbind, traj_rows)
misspec_summary <- do.call(rbind, misspec_rows)
misspec_summary$error_label <- as.character(misspec_summary$error_label)

write.csv(traj_summary, file.path(output_dir, "rt_trajectory_summary.csv"), row.names = FALSE)
write.csv(misspec_summary, file.path(output_dir, "rt_misspec_trajectory_summary.csv"),
          row.names = FALSE)

pdf(file.path(output_dir, "rt_trajectory_pages.pdf"), width = 12, height = 4.6)
for (p in pages_correct) print(p)
invisible(dev.off())

pdf(file.path(output_dir, "rt_misspec_trajectory_pages.pdf"), width = 11, height = 8)
for (p in pages_misspec) print(p)
invisible(dev.off())

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
cat(sprintf("\nWrote results to %s/\n", output_dir))
