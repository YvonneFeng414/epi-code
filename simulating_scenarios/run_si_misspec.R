# ============================================================
# SI MISSPECIFICATION: FOUR ESTIMATORS, GRADUAL Rt 1.1 -> 1.3
#
# ENTRY POINT
#
#   Rscript run_si_misspec.R          # full run, 100 epidemics
#   Rscript run_si_misspec.R smoke    # smoke test, 4 epidemics
#
# Run from simulating_scenarios/.
#
# The same epidemics as run_trajectories.R (same seed, same
# simulators, always the TRUE serial interval): Poisson data for
# EpiEstim and EpiFilter, negative-binomial data (rho =
# epilps_nb_rho) for both EpiLPS arms. Refitted by every
# estimator under five assumed SIs - si_arms in config.R:
# correct, then mean and sd each moved -/+ si_misspec_rel. The design is paired, so any
# difference between arms is the misspecification alone, and
# the Correct-SI arm reproduces run_trajectories.R exactly.
#
# Common burn-in: everything is scored and drawn from day
# si_misspec_score_from (computed below; see config.R), because EpiEstim's
# own start day moves with the assumed SI mean.
#
# Outputs, in <output_dir>_si_misspec/ (e.g. results_si_misspec/):
#   si_misspec_plots.pdf     - the only figure file
#       page 1: estimator rows x (Rt over time | time-varying bias)
#       page 2: estimator rows x (coverage | CrI width), per day
#   si_misspec_summary.csv   - per-day rows, estimator x arm
#   si_misspec_metrics.csv   - scalar metrics, estimator x arm x window
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

misspec_dir <- paste0(output_dir, "_si_misspec")
dir.create(misspec_dir, showWarnings = FALSE)

cat(sprintf("Mode: %s | I0 = %g | n_sim = %d | cores = %d | output: %s/\n",
            if (SMOKE) "SMOKE" else "FULL", traj_I0, traj_n_sim, n_cores,
            misspec_dir))


# ============================================================
# 2. Target and simulated epidemics - identical to
#    run_trajectories.R, so the two experiments share data
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

# EpiLPS's own data - identical to run_trajectories.R
set.seed(sim_seed)
simulations_nb <- simulate_many_nb(R_t = R_t, w = w, time = traj_time,
                                   I0 = traj_I0, n_sim = traj_n_sim,
                                   rho = epilps_nb_rho)

cat(sprintf("EpiLPS data: negative binomial, rho = %g\n", epilps_nb_rho))
cat(sprintf("Simulated %d epidemics per data set with the TRUE SI (mean %.2f, sd %.2f)\n",
            ncol(simulations), mean_si, sd_si))


# ============================================================
# 3. The arms, and what truncation at lag 30 costs each one
#
# EpiFilter and EpiLPS get make_si(), truncated at lag 30 and
# renormalised; EpiEstim builds its own untruncated SI. The
# mass dropped is what the two groups see differently.
# ============================================================

si_arms$mass_past_lag30 <- mapply(function(m, s) {
  full <- EpiEstim::discr_si(0:300, m, s)
  sum(full[-(1:31)]) / sum(full)
}, si_arms$a_mean, si_arms$a_sd)

cat("\nAssumed SI per arm:\n")
print(transform(si_arms, mass_past_lag30 = signif(mass_past_lag30, 2)),
      row.names = FALSE)


# ============================================================
# 4. Fit every estimator under every arm
# ============================================================

# EpiEstim warns on every early window of every epidemic; that
# is expected (those windows come back NA) and would bury
# anything else in the log.
quiet_early <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("too early", conditionMessage(w))) invokeRestart("muffleWarning")
  })
}

fitters <- list(
  epiestim    = function(m, s) quiet_early(fit_epiestim(simulations, dates, m, s)),
  epifilter   = function(m, s) fit_epifilter(simulations, m, s),
  epilps_map  = function(m, s) fit_epilps_map(simulations_nb, m, s),
  epilps_mala = function(m, s) fit_epilps_mala(simulations_nb, m, s)
)

fits <- list()

for (a in seq_len(nrow(si_arms))) {
  arm <- si_arms$arm[a]
  cat(sprintf("\n--- %s (assumed mean %.2f, sd %.2f) ---\n",
              arm, si_arms$a_mean[a], si_arms$a_sd[a]))

  fits[[arm]] <- list()
  for (k in names(fitters)) {
    fits[[arm]][[k]] <- fitters[[k]](si_arms$a_mean[a], si_arms$a_sd[a])
    cat(sprintf("  %-22s %6.1f s   failed replicates: %d\n",
                estimator_labels[[k]], fits[[arm]][[k]]$secs,
                fits[[arm]][[k]]$n_failed))
  }
}


# ============================================================
# 5. Per-day summaries
# ============================================================

summary_df <- bind_rows(lapply(seq_len(nrow(si_arms)), function(a) {
  arm <- si_arms$arm[a]
  bind_rows(lapply(names(fitters), function(k) {
    data.frame(
      estimator = estimator_labels[[k]],
      data_model = estimator_data[[k]],
      arm = arm,
      a_mean = si_arms$a_mean[a],
      a_sd = si_arms$a_sd[a],
      subscenario = scenario_nm,
      I0 = traj_I0,
      n_sim_used = ncol(simulations),
      summarise_rt_trajectories(fits[[arm]][[k]], R_t = R_t),
      stringsAsFactors = FALSE
    )
  }))
}))

write.csv(summary_df, file.path(misspec_dir, "si_misspec_summary.csv"),
          row.names = FALSE)

# Common burn-in (see config.R): the latest day on which any
# estimator x arm first has an estimate - EpiEstim's own start
# moves with the assumed mean - and never before common_from_day.
first_day <- summary_df %>%
  filter(!is.na(mean_est)) %>%
  group_by(estimator, arm) %>%
  summarise(first = min(t_end), .groups = "drop")

si_misspec_score_from <- as.integer(max(common_from_day, first_day$first))

# The phase windows below start at day 25
if (si_misspec_score_from >= 25L) {
  stop("Common burn-in reaches day ", si_misspec_score_from,
       "; the phase windows assume it ends before day 25.")
}

cat(sprintf("\nCommon burn-in: scored from day %d (latest first-estimated day; set by %s)\n",
            si_misspec_score_from,
            paste(unique(with(first_day[first_day$first == max(first_day$first), ],
                              paste(estimator, "/", arm))), collapse = ", ")))


# ============================================================
# 6. Scalar metrics over replicate-day cells
#
# The scored window is si_misspec_score_from..T, the same days
# for every estimator and arm. The four phases split it: EpiLPS
# start-up (to day 24), flat 1.1 (25-49), the rise (50-69),
# flat 1.3 (70-100).
# ============================================================

scored_window <- paste0(si_misspec_score_from, "-", traj_time)

metric_windows <- setNames(
  list(si_misspec_score_from:traj_time, si_misspec_score_from:24,
       25:49, 50:69, 70:traj_time),
  c(scored_window, paste0(si_misspec_score_from, "-24"),
    "25-49", "50-69", "70-100")
)

cell_metrics <- function(traj, days) {
  rows <- which(traj$t_end %in% days)
  est <- traj$est[rows, , drop = FALSE]
  lo  <- traj$lower[rows, , drop = FALSE]
  hi  <- traj$upper[rows, , drop = FALSE]
  tr  <- matrix(R_t[traj$t_end[rows]], nrow = nrow(est), ncol = ncol(est))
  err <- est - tr

  data.frame(
    n_cells = sum(is.finite(est)),
    bias = mean(err, na.rm = TRUE),
    rmse = sqrt(mean(err^2, na.rm = TRUE)),
    coverage = mean(lo <= tr & hi >= tr, na.rm = TRUE),
    ci_width = mean(hi - lo, na.rm = TRUE)
  )
}

metrics_df <- bind_rows(lapply(seq_len(nrow(si_arms)), function(a) {
  arm <- si_arms$arm[a]
  bind_rows(lapply(names(fitters), function(k) {
    traj <- fits[[arm]][[k]]
    bind_rows(lapply(names(metric_windows), function(win) {
      data.frame(
        estimator = estimator_labels[[k]],
        data_model = estimator_data[[k]],
        arm = arm,
        a_mean = si_arms$a_mean[a],
        a_sd = si_arms$a_sd[a],
        window = win,
        cell_metrics(traj, metric_windows[[win]]),
        n_failed = traj$n_failed,
        secs = round(traj$secs, 1),
        mean_rho = if (is.null(traj$rho)) NA_real_ else mean(traj$rho, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }))
  }))
}))

write.csv(metrics_df, file.path(misspec_dir, "si_misspec_metrics.csv"),
          row.names = FALSE)

cat("\nMetrics over the scored window (days ", scored_window, "):\n", sep = "")
print(
  metrics_df %>%
    filter(window == scored_window) %>%
    select(estimator, arm, bias, rmse, coverage, ci_width, n_failed) %>%
    mutate(across(where(is.double), ~ signif(.x, 3))),
  row.names = FALSE
)


# ============================================================
# 7. The figure: one PDF, two pages
#
# Layout of queens_rt/combined/trajectory_panel.png: estimator
# rows, lettered tinted strips, one line per assumed SI with a
# ribbon for the mean 95% CrI, one shared legend on top.
# ============================================================

arm_levels <- c("True Rt", si_arms$arm)

# In si_arms order: correct, mean -, mean +, sd -, sd +
arm_colours <- setNames(
  c("#000000", "#2471A3", "#E8853A", "#C0392B", "#7F7F7F", "#1A9E5B"),
  arm_levels
)

arm_lty <- setNames(
  c("solid", "solid", "longdash", "longdash", "longdash", "longdash"),
  arm_levels
)

row_fill <- setNames(c("#E2EFDA", "#DEEAF6", "#FBE4D5", "#E2D4F0"),
                     unname(estimator_labels))

per_day <- summary_df %>%
  filter(t_end >= si_misspec_score_from) %>%
  mutate(
    estimator = factor(estimator, levels = unname(estimator_labels)),
    arm = factor(arm, levels = arm_levels),
    bias = mean_est - true_R,
    bias_lo = mean_lower - true_R,
    bias_hi = mean_upper - true_R,
    ci_width = mean_upper - mean_lower
  )

# The truth is known on every day, so it spans the whole axis
truth_df <- data.frame(t_end = seq_len(traj_time), true_R = R_t,
                       arm = factor("True Rt", levels = arm_levels))

# Axis ranges from display_from_day on (view only; see config.R)
view <- filter(per_day, t_end >= display_from_day)
rt_lim    <- range(c(view$mean_lower, view$mean_upper, view$true_R), na.rm = TRUE)
bias_lim  <- range(c(view$bias_lo, view$bias_hi, 0), na.rm = TRUE)
width_lim <- c(0, max(view$ci_width, na.rm = TRUE))

panel_theme <- theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 0.4),
        axis.title = element_text(size = 8.5),
        axis.text = element_text(size = 7.5),
        plot.margin = margin(2, 4, 2, 2))

strip_theme <- function(est) {
  theme(strip.background = element_rect(fill = row_fill[[est]], colour = "grey30",
                                        linewidth = 0.4),
        strip.text = element_text(size = 8.5, face = "bold", hjust = 0,
                                  margin = margin(3, 3, 3, 5)))
}

# key_levels: page 1 keys the black reference as "True Rt";
# page 2 has no truth line, so its key lists the arms only.
arm_scales <- function(key_levels) list(
  scale_colour_manual(values = arm_colours, limits = key_levels,
                      breaks = key_levels, drop = FALSE),
  scale_fill_manual(values = arm_colours, breaks = key_levels, guide = "none"),
  scale_linetype_manual(values = arm_lty, limits = key_levels,
                        breaks = key_levels, drop = FALSE)
)

transition_marks <- geom_vline(xintercept = c(49.5, 69.5), linetype = "dotted",
                               colour = "grey50")

# One panel: `y` is the line, `lo`/`hi` the ribbon (NULL = none),
# `ref` a reference series drawn as the "True Rt" key entry,
# `hline` an unkeyed reference level.
make_panel <- function(est, letter, what, y, lo, hi, ylim, ylab, ref = NULL,
                       hline = NULL, key_levels = arm_levels) {
  d <- filter(per_day, estimator == est)
  d$strip <- sprintf("%s. %s (%s)", letter, what, est)

  p <- ggplot(d, aes(x = t_end)) + transition_marks
  if (!is.null(hline)) {
    p <- p + geom_hline(yintercept = hline, colour = "black", linewidth = 0.7)
  }
  if (!is.null(lo)) {
    p <- p + geom_ribbon(aes(ymin = .data[[lo]], ymax = .data[[hi]], fill = arm),
                         alpha = 0.13, colour = NA, na.rm = TRUE)
  }
  if (!is.null(ref)) {
    p <- p + geom_line(data = ref, aes(y = y, colour = arm, linetype = arm),
                       linewidth = 0.9)
  }
  p +
    geom_line(aes(y = .data[[y]], colour = arm, linetype = arm),
              linewidth = 0.5, na.rm = TRUE) +
    facet_wrap(~ strip) +
    arm_scales(key_levels) +
    coord_cartesian(xlim = c(0, traj_time), ylim = ylim) +
    labs(x = "Day", y = ylab) +
    panel_theme + strip_theme(est)
}

ests <- unname(estimator_labels)
lets <- matrix(LETTERS[seq_len(2 * length(ests))], ncol = 2, byrow = TRUE)

truth_ref <- transform(truth_df, y = true_R)
zero_ref  <- data.frame(t_end = c(0, traj_time), y = 0,
                        arm = factor("True Rt", levels = arm_levels))

collect_legend <- function(rows, title, subtitle) {
  wrap_plots(rows, ncol = 1) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = title, subtitle = subtitle,
      theme = theme(plot.title = element_text(face = "bold", size = 13),
                    plot.subtitle = element_text(size = 8.5))
    ) &
    theme(legend.position = "top", legend.direction = "horizontal",
          legend.title = element_blank(), legend.background = element_blank(),
          legend.key = element_blank(), legend.text = element_text(size = 9),
          legend.margin = margin(0, 0, 2, 0))
}

page_title <- paste0(
  "SI misspecification: gradual Rt ",
  gsub("→", "->", target$true_rt_display, fixed = TRUE),
  " (", target$transition_duration, "-day transition)"
)

subtitle <- paste0(
  "Data simulated with the TRUE SI (mean ", mean_si, ", sd ", sd_si,
  "); each line re-fits the same ", ncol(simulations), " epidemics with the ",
  "assumed SI named in the key. I0 = ", traj_I0, ", T = ", traj_time, ".\n",
  "Lines = mean over epidemics; ribbons = mean 95% CrI. Dotted lines mark the ",
  "transition (days 50-69). Scored and drawn from day ",
  si_misspec_score_from, " (common burn-in); axis ranges from day ",
  display_from_day, " (view only).\n",
  "Data: Poisson for EpiEstim and EpiFilter; negative binomial (rho = ",
  epilps_nb_rho, ") for EpiLPS, which is fitted to days ", epilps_start, "-",
  traj_time, " and runs far off the scale on its first days."
)

page1 <- collect_legend(
  lapply(seq_along(ests), function(i) {
    make_panel(ests[i], lets[i, 1], "Estimated Rt over time", "mean_est",
               "mean_lower", "mean_upper", rt_lim, "Rt", truth_ref) |
      make_panel(ests[i], lets[i, 2], "Time-varying bias", "bias",
                 "bias_lo", "bias_hi", bias_lim, "Bias of Rt", zero_ref)
  }),
  page_title, subtitle
)

page2 <- collect_legend(
  lapply(seq_along(ests), function(i) {
    make_panel(ests[i], lets[i, 1], "95% CrI coverage", "coverage",
               NULL, NULL, c(0, 1), "Coverage", hline = 0.95,
               key_levels = si_arms$arm) |
      make_panel(ests[i], lets[i, 2], "95% CrI width", "ci_width",
                 NULL, NULL, width_lim, "Mean CrI width",
                 key_levels = si_arms$arm)
  }),
  paste0(page_title, ": coverage and interval width"),
  paste0("Black line in the coverage panels = nominal 0.95. n_sim = ",
         ncol(simulations), if (SMOKE) "; with 4 epidemics coverage moves in steps of 0.25." else ".")
)

pdf_file <- file.path(misspec_dir, "si_misspec_plots.pdf")
pdf(pdf_file, width = 11, height = 12.5)
print(page1)
print(page2)
invisible(dev.off())

cat("\nWrote", pdf_file, "\n")

writeLines(capture.output(sessionInfo()),
           file.path(misspec_dir, "sessionInfo.txt"))


# ============================================================
# 8. Checks - after the outputs are written
# ============================================================

problems <- character(0)

expected_R <- c(rep(1.1, 49), seq(1.11, 1.30, by = 0.01), rep(1.3, 31))
if (!isTRUE(all.equal(R_t, expected_R, tolerance = 1e-9))) {
  problems <- c(problems, "true Rt path does not match the plan's table")
}

# The Correct-SI arm must reproduce run_trajectories.R: same
# seed, same epidemics, same fitters, same SI.
ref_file <- file.path(output_dir, "rt_trajectory_summary.csv")
if (file.exists(ref_file)) {
  ref <- read.csv(ref_file, stringsAsFactors = FALSE)
  if (all(ref$I0 == traj_I0) && all(ref$n_sim_used == ncol(simulations))) {
    mine <- filter(summary_df, arm == "Correct SI")
    cols <- c("mean_est", "mean_lower", "mean_upper", "coverage")
    for (e in unique(mine$estimator)) {
      a <- mine[mine$estimator == e, c("t_end", cols)]
      b <- ref[ref$estimator == e, c("t_end", cols)]
      if (!isTRUE(all.equal(a, b, tolerance = 1e-8, check.attributes = FALSE))) {
        problems <- c(problems, paste("Correct-SI arm differs from", ref_file, "for", e))
      }
    }
    if (!any(grepl("differs", problems))) {
      cat("Correct-SI arm reproduces", ref_file, "for all four estimators.\n")
    }
  } else {
    cat("Skipped reproduction check:", ref_file, "has a different I0 or n_sim.\n")
  }
} else {
  cat("Skipped reproduction check: no", ref_file, "\n")
}

# Failed fits under misspecification are a RESULT - reported,
# never a reason to stop.
failed <- metrics_df %>%
  filter(window == scored_window, n_failed > 0) %>%
  select(estimator, arm, n_failed)
if (nrow(failed) > 0L) {
  cat("\n*** Fits that returned no estimate (a result, not a crash):\n")
  print(failed, row.names = FALSE)
} else {
  cat("Every estimator fitted every epidemic under every arm.\n")
}

if (length(problems) > 0L) {
  stop("Checks failed:\n  ", paste(problems, collapse = "\n  "))
}

cat("All checks passed.\n")
