# ============================================================
# Rt TRAJECTORY: FOUR ESTIMATORS, GRADUAL Rt 1.1 -> 1.3
#
# ENTRY POINT
#
#   Rscript run_trajectories.R          # full run, 100 epidemics
#   Rscript run_trajectories.R smoke    # smoke test, 4 epidemics
#
# Run from simulating_scenarios/.
#
# Simulates epidemics from the gradual20_1.1_to_1.3 truth (flat
# 1.1 to day 49, +0.01/day on days 50-69, flat 1.3 from day 70)
# under each estimator's own observation model: Poisson for
# EpiEstim and EpiFilter, negative binomial (rho =
# epilps_nb_rho) for both EpiLPS arms. Same seed and renewal
# process for both, so EpiEstim/EpiFilter share one paired set
# and MAP/MALA share the other.
#
# The simulator, scenario definitions, trajectory summariser
# and panel functions are sourced unchanged from the seed study
# (~/Desktop/Epiestim_simu), so the EpiEstim row reproduces the
# seed study's own trajectory figure.
#
# Outputs, in results/ (or results_smoke/):
#   rt_trajectory_plots.pdf     - the only figure file
#       page 1: estimator rows x (one epidemic | averaged)
#       page 2: per-day bias, coverage and CrI width overlaid
#   rt_trajectory_summary.csv   - per-day rows, all estimators
#   rt_trajectory_metrics.csv   - scalar metrics per estimator
#   sessionInfo.txt
# ============================================================

suppressPackageStartupMessages({
  library(EpiEstim)
  library(dplyr)
  library(parallel)
})

if (!file.exists("estimators.R")) {
  stop("Run this from the simulating_scenarios/ directory.")
}


# ============================================================
# 1. Code
# ============================================================

source("config.R")

for (f in c("simulation_functions.R", "scenarios.R",
            "sensitivity_plot_functions.R", "rt_trajectory_functions.R")) {
  source(file.path(seed_dir, f))
}

source("estimators.R")
source("simulation_nb.R")

dir.create(output_dir, showWarnings = FALSE)

cat(sprintf("Mode: %s | n_sim = %d | cores = %d | output: %s/\n",
            if (SMOKE) "SMOKE" else "FULL", traj_n_sim, n_cores, output_dir))


# ============================================================
# 2. Target and simulated epidemics
# ============================================================

scenario_nm <- traj_subscenarios[[1L]]
target <- build_scenarios(time = traj_time)[[scenario_nm]]
if (is.null(target)) stop("build_scenarios() does not define ", scenario_nm)

R_t <- target$R_t
w <- EpiEstim::discr_si(seq(0, traj_time), mean_si, sd_si)
dates <- start_date + seq_len(traj_time) - 1

set.seed(sim_seed)
simulations <- simulate_many(R_t = R_t, w = w, time = traj_time,
                             I0 = traj_I0, n_sim = traj_n_sim)

# EpiLPS's own data: the same renewal process, negative-binomial
# draws (simulation_nb.R). Same seed, so it is reproducible and
# identical across run_trajectories.R and run_si_misspec.R.
set.seed(sim_seed)
simulations_nb <- simulate_many_nb(R_t = R_t, w = w, time = traj_time,
                                   I0 = traj_I0, n_sim = traj_n_sim,
                                   rho = epilps_nb_rho)

sim_data <- list(poisson = simulations, negbin = simulations_nb)

for (d in names(sim_data)) {
  x <- sim_data[[d]]
  cat(sprintf("Simulated %d %s epidemics; final-day incidence median %.0f (range %.0f-%.0f)\n",
              ncol(x), if (d == "negbin") sprintf("negative-binomial (rho = %g)", epilps_nb_rho) else "Poisson",
              median(x[traj_time, ]), min(x[traj_time, ]), max(x[traj_time, ])))
}


# ============================================================
# 3. Fit all four
# ============================================================

fits <- list()

cat("\nFitting EpiEstim ...\n")
fits$epiestim <- fit_epiestim(simulations, dates)

cat("Fitting EpiFilter ...\n")
fits$epifilter <- fit_epifilter(simulations)

cat("Fitting EpiLPS (MAP) ...\n")
fits$epilps_map <- fit_epilps_map(simulations_nb)

cat("Fitting EpiLPS (MALA) ...\n")
fits$epilps_mala <- fit_epilps_mala(simulations_nb)

for (k in names(fits)) {
  cat(sprintf("  %-22s %7.1f s   failed replicates: %d\n",
              estimator_labels[[k]], fits[[k]]$secs, fits[[k]]$n_failed))
  if (length(fits[[k]]$errors) > 0L) {
    cat(paste0("    ", head(fits[[k]]$errors, 5L), collapse = "\n"), "\n")
  }
}


# ============================================================
# 4. Per-day summaries
# ============================================================

summaries <- lapply(fits, summarise_rt_trajectories, R_t = R_t)

# One representative epidemic per data set: the Poisson one for
# EpiEstim / EpiFilter, the negative-binomial one for EpiLPS.
sim_index <- vapply(names(fits), function(k) {
  pick_representative_sim(sim_data[[estimator_data[[k]]]])
}, integer(1))

summary_df <- bind_rows(lapply(names(summaries), function(k) {
  data.frame(
    estimator = estimator_labels[[k]],
    data_model = estimator_data[[k]],
    subscenario = scenario_nm,
    true_rt = target$true_rt_display,
    T = traj_time,
    I0 = traj_I0,
    transition_duration = target$transition_duration,
    n_sim_used = ncol(simulations),
    representative_sim = sim_index[[k]],
    summaries[[k]],
    stringsAsFactors = FALSE
  )
}))

write.csv(summary_df, file.path(output_dir, "rt_trajectory_summary.csv"),
          row.names = FALSE)


# ============================================================
# 5. Scalar metrics, over every replicate-day cell
#
# "own days": every day that estimator reports.
# "common":   days common_from_day..T, which all four report.
# ============================================================

cell_metrics <- function(traj, true_R, keep_rows) {
  est <- traj$est[keep_rows, , drop = FALSE]
  lo  <- traj$lower[keep_rows, , drop = FALSE]
  hi  <- traj$upper[keep_rows, , drop = FALSE]
  tr  <- matrix(true_R[keep_rows], nrow = nrow(est), ncol = ncol(est))
  err <- est - tr

  data.frame(
    n_cells = sum(is.finite(est)),
    bias = mean(err, na.rm = TRUE),
    rmse = sqrt(mean(err^2, na.rm = TRUE)),
    coverage = mean(lo <= tr & hi >= tr, na.rm = TRUE),
    ci_width = mean(hi - lo, na.rm = TRUE)
  )
}

metrics_df <- bind_rows(lapply(names(fits), function(k) {
  traj <- fits[[k]]
  true_R <- summaries[[k]]$true_R
  own <- which(rowSums(is.finite(traj$est)) > 0L)
  common <- which(traj$t_end >= common_from_day)

  bind_rows(
    data.frame(window = sprintf("own days (%d-%d)", min(traj$t_end[own]),
                                max(traj$t_end[own])),
               cell_metrics(traj, true_R, own)),
    data.frame(window = sprintf("common (%d-%d)", common_from_day, traj_time),
               cell_metrics(traj, true_R, common))
  ) %>%
    mutate(estimator = estimator_labels[[k]],
           data_model = estimator_data[[k]],
           n_failed = traj$n_failed,
           secs = round(traj$secs, 1),
           mean_rho = if (is.null(traj$rho)) NA_real_ else mean(traj$rho, na.rm = TRUE),
           .before = 1)
}))

write.csv(metrics_df, file.path(output_dir, "rt_trajectory_metrics.csv"),
          row.names = FALSE)

cat("\nMetrics:\n")
print(metrics_df %>% mutate(across(where(is.double), ~ signif(.x, 4))),
      row.names = FALSE)


# ============================================================
# 6. The figure: one PDF, two pages
# ============================================================

page_title <- paste0(
  "Gradual Rt ", gsub("→", "->", target$true_rt_display, fixed = TRUE),
  " (", target$transition_duration, "-day transition)"
)

# One y range for all eight panels, so the rows compare directly.
# Taken from display_from_day on (see config.R): EpiLPS's first
# fitted days run into the hundreds and would flatten every
# panel. Those days are still drawn; coord_cartesian() only
# clips the view.
in_common <- function(k) summaries[[k]]$t_end >= display_from_day

y_limits <- rt_y_limits(
  lowers = unlist(lapply(names(fits), function(k) list(
    fits[[k]]$lower[in_common(k), sim_index[[k]]],
    summaries[[k]]$mean_lower[in_common(k)],
    summaries[[k]]$mc_lower[in_common(k)]))),
  uppers = unlist(lapply(names(fits), function(k) list(
    fits[[k]]$upper[in_common(k), sim_index[[k]]],
    summaries[[k]]$mean_upper[in_common(k)],
    summaries[[k]]$mc_upper[in_common(k)]))),
  true_R = R_t
)

# The truth climbs to the top right, so keep the key top left
legend_top_left <- theme(legend.position.inside = c(0.03, 0.97),
                         legend.justification = c(0, 1))

row_panels <- lapply(names(fits), function(k) {
  lab <- estimator_labels[[k]]
  list(
    plot_rt_single(summaries[[k]], fits[[k]], sim_index[[k]], y_limits,
                   title = paste0(lab, ": one simulated epidemic")) +
      legend_top_left,
    plot_rt_montecarlo(summaries[[k]], y_limits, ncol(simulations),
                       title = paste0(lab, ": averaged over ",
                                      ncol(simulations), " epidemics")) +
      legend_top_left
  )
})

subtitle <- paste0(
  "Renewal data, I0 = ", traj_I0, ", T = ", traj_time,
  ", n_sim = ", ncol(simulations), ", correct serial interval (mean = ",
  mean_si, ", sd = ", sd_si, "). Each estimator gets data from its own ",
  "observation model: EpiEstim and EpiFilter Poisson,\nEpiLPS negative binomial ",
  "(rho = ", epilps_nb_rho, "), same seed. Left column: the representative ",
  "epidemic of each data set (Poisson #", sim_index[["epiestim"]],
  ", negative binomial #", sim_index[["epilps_map"]], ").\n",
  "Band = 95% credible interval; dashed = 2.5-97.5 percentile spread of the ",
  "point estimate across simulations. EpiEstim: 1-day window. EpiFilter: ",
  "smoother, eta = ", epifilter_eta, ", from day ", epifilter_start, ".\n",
  "EpiLPS: K = ", epilps_K, ", fitted to days ", epilps_start, "-", traj_time,
  "; MALA ", mcmc_niter, " iterations, ", mcmc_burnin, " burn-in."
)

page1 <- wrap_plots(unlist(row_panels, recursive = FALSE), ncol = 2) +
  plot_annotation(
    title = page_title, subtitle = subtitle,
    theme = theme(plot.title = element_text(face = "bold", size = 13),
                  plot.subtitle = element_text(size = 8.5))
  )

# Page 2: the four estimators overlaid, per day
est_colours <- setNames(c("#2471A3", "#C0392B", "#1E8449", "#7D3C98"),
                        unname(estimator_labels))

per_day <- summary_df %>%
  mutate(bias = mean_est - true_R, ci_width = mean_upper - mean_lower,
         estimator = factor(estimator, levels = unname(estimator_labels)))

transition_marks <- list(
  geom_vline(xintercept = c(49.5, 69.5), linetype = "dotted", colour = "grey50")
)

overlay_panel <- function(y, ylab, title, hline = NULL) {
  p <- ggplot(per_day, aes(x = t_end, y = .data[[y]], colour = estimator)) +
    transition_marks
  if (!is.null(hline)) {
    p <- p + geom_hline(yintercept = hline, linetype = "dashed", colour = "grey30")
  }
  p +
    geom_line(linewidth = 0.6, na.rm = TRUE) +
    scale_colour_manual(values = est_colours) +
    labs(title = title, x = "day", y = ylab) +
    misspec_theme(c(0.03, 0.03), c(0, 0)) +
    theme(legend.position = "bottom")
}

# Axis ranges from display_from_day, for the same reason as page 1
per_day_common <- filter(per_day, t_end >= display_from_day)

page2 <- wrap_plots(
  overlay_panel("bias", "mean(estimate) - true Rt", "Per-day bias", hline = 0) +
    coord_cartesian(ylim = range(per_day_common$bias, na.rm = TRUE)),
  overlay_panel("coverage", "coverage", "Per-day 95% CrI coverage", hline = 0.95) +
    coord_cartesian(ylim = c(0, 1)),
  overlay_panel("ci_width", "mean CrI width", "Per-day 95% CrI width") +
    coord_cartesian(ylim = c(0, max(per_day_common$ci_width, na.rm = TRUE))),
  ncol = 1
) +
  plot_annotation(
    title = paste0(page_title, ": per-day comparison"),
    subtitle = paste0("Dotted lines mark the transition (days 50-69). Axis ranges ",
                      "are set from days ", display_from_day, "-", traj_time,
                      " (view only; metrics start at day ", common_from_day, "). ",
                      "n_sim = ", ncol(simulations), "; with few replicates ",
                      "coverage moves in coarse steps."),
    theme = theme(plot.title = element_text(face = "bold", size = 13),
                  plot.subtitle = element_text(size = 8.5))
  )

pdf_file <- file.path(output_dir, "rt_trajectory_plots.pdf")
pdf(pdf_file, width = 12, height = 15)
print(page1)
print(page2)
invisible(dev.off())

cat("\nWrote", pdf_file, "\n")

writeLines(capture.output(sessionInfo()),
           file.path(output_dir, "sessionInfo.txt"))


# ============================================================
# 7. Checks - after the outputs are written, so a failure
#    still leaves everything on disk to inspect
# ============================================================

expected_R <- c(rep(1.1, 49), seq(1.11, 1.30, by = 0.01), rep(1.3, 31))

problems <- character(0)

if (!isTRUE(all.equal(R_t, expected_R, tolerance = 1e-9))) {
  problems <- c(problems, "true Rt path does not match the plan's table")
}

for (k in names(summaries)) {
  if (!isTRUE(all.equal(summaries[[k]]$true_R, expected_R[summaries[[k]]$t_end],
                        tolerance = 1e-9))) {
    problems <- c(problems, paste(estimator_labels[[k]], "true_R column is wrong"))
  }
  if (fits[[k]]$n_failed > 0L) {
    problems <- c(problems, sprintf("%s: %d replicate(s) returned no estimate",
                                    estimator_labels[[k]], fits[[k]]$n_failed))
  }
}

if (length(problems) > 0L) {
  stop("Checks failed:\n  ", paste(problems, collapse = "\n  "))
}

cat("All checks passed: true Rt path matches, every estimator fitted every replicate.\n")
