# ============================================================
# Redraw page 1 of si_misspec_plots.pdf (Rt over time | bias)
# with the estimates drawn from day 25 only
#
#   Rscript replot_si_misspec.R          # results_si_misspec/
#   I0=100 Rscript replot_si_misspec.R   # results_I0_100_si_misspec/
#
# Run from simulating_scenarios/ after run_si_misspec.R. No
# fitting: reads the per-day si_misspec_summary.csv and writes
# si_misspec_plots_from_day25.pdf next to the original, which
# is left untouched.
#
# Changes from run_si_misspec.R's figure:
#   - every estimator line and ribbon starts at draw_from_day
#     instead of the common burn-in. 25 is display_from_day, the
#     day the axis ranges were already taken from, and the start
#     of the flat-1.1 phase; before it EpiLPS's renewal sum is
#     still filling up and its lines ran off the scale. The
#     black references - true Rt, zero bias - still span the
#     whole axis;
#   - page 1 only (the coverage / width page is dropped), with
#     no title or subtitle, so the page is the panels and key;
#   - EpiFilter's strips read "EpiFilter", without "(smoother)".
# Metrics are not recomputed; si_misspec_metrics.csv still
# scores from the common burn-in.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

if (!file.exists("run_si_misspec.R")) {
  stop("Run this from the simulating_scenarios/ directory.")
}

source("config.R")
for (f in c("simulation_functions.R", "scenarios.R")) {
  source(file.path(seed_dir, f))
}

draw_from_day <- display_from_day

misspec_dir  <- paste0(output_dir, "_si_misspec")
summary_file <- file.path(misspec_dir, "si_misspec_summary.csv")
if (!file.exists(summary_file)) {
  stop("Cannot find ", summary_file, ". Run run_si_misspec.R first.")
}

summary_df <- read.csv(summary_file, stringsAsFactors = FALSE)
n_sim <- summary_df$n_sim_used[1L]

scenario_nm <- traj_subscenarios[[1L]]
target <- build_scenarios(time = traj_time)[[scenario_nm]]
if (is.null(target)) stop("build_scenarios() does not define ", scenario_nm)
R_t <- target$R_t

# The common burn-in, recomputed exactly as run_si_misspec.R
# does, for the console line only.
first_day <- summary_df %>%
  filter(!is.na(mean_est)) %>%
  group_by(estimator, arm) %>%
  summarise(first = min(t_end), .groups = "drop")
si_misspec_score_from <- as.integer(max(common_from_day, first_day$first))

cat(sprintf("Read %s (I0 = %g, n_sim = %d); scored from day %d, drawn from day %d\n",
            summary_file, summary_df$I0[1L], n_sim, si_misspec_score_from,
            draw_from_day))


# ============================================================
# The figure - run_si_misspec.R section 7, with per_day cut at
# draw_from_day
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
  filter(t_end >= draw_from_day) %>%
  mutate(
    estimator = factor(estimator, levels = unname(estimator_labels)),
    arm = factor(arm, levels = arm_levels),
    bias = mean_est - true_R,
    bias_lo = mean_lower - true_R,
    bias_hi = mean_upper - true_R
  )

# The truth is known on every day, so it spans the whole axis
truth_df <- data.frame(t_end = seq_len(traj_time), true_R = R_t,
                       arm = factor("True Rt", levels = arm_levels))

rt_lim    <- range(c(per_day$mean_lower, per_day$mean_upper, R_t), na.rm = TRUE)
bias_lim  <- range(c(per_day$bias_lo, per_day$bias_hi, 0), na.rm = TRUE)

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

# The key lists the black reference as "True Rt", then the arms.
arm_scales <- list(
  scale_colour_manual(values = arm_colours, limits = arm_levels,
                      breaks = arm_levels, drop = FALSE),
  scale_fill_manual(values = arm_colours, breaks = arm_levels, guide = "none"),
  scale_linetype_manual(values = arm_lty, limits = arm_levels,
                        breaks = arm_levels, drop = FALSE)
)

transition_marks <- geom_vline(xintercept = c(49.5, 69.5), linetype = "dotted",
                               colour = "grey50")

# Strip names drop EpiFilter's "(smoother)"; row_fill and the
# data stay keyed by the full estimator label.
strip_name <- function(est) sub(" (smoother)", "", est, fixed = TRUE)

# One panel: `y` is the line, `lo`/`hi` the ribbon, `ref` a
# reference series drawn as the "True Rt" key entry.
make_panel <- function(est, letter, what, y, lo, hi, ylim, ylab, ref) {
  d <- filter(per_day, estimator == est)
  d$strip <- sprintf("%s. %s (%s)", letter, what, strip_name(est))

  ggplot(d, aes(x = t_end)) + transition_marks +
    geom_ribbon(aes(ymin = .data[[lo]], ymax = .data[[hi]], fill = arm),
                alpha = 0.13, colour = NA, na.rm = TRUE) +
    geom_line(data = ref, aes(y = y, colour = arm, linetype = arm),
              linewidth = 0.9) +
    geom_line(aes(y = .data[[y]], colour = arm, linetype = arm),
              linewidth = 0.5, na.rm = TRUE) +
    facet_wrap(~ strip) +
    arm_scales +
    coord_cartesian(xlim = c(0, traj_time), ylim = ylim) +
    labs(x = "Day", y = ylab) +
    panel_theme + strip_theme(est)
}

ests <- unname(estimator_labels)
lets <- matrix(LETTERS[seq_len(2 * length(ests))], ncol = 2, byrow = TRUE)

truth_ref <- transform(truth_df, y = true_R)
zero_ref  <- data.frame(t_end = c(0, traj_time), y = 0,
                        arm = factor("True Rt", levels = arm_levels))

page <- wrap_plots(
  lapply(seq_along(ests), function(i) {
    make_panel(ests[i], lets[i, 1], "Estimated Rt over time", "mean_est",
               "mean_lower", "mean_upper", rt_lim, "Rt", truth_ref) |
      make_panel(ests[i], lets[i, 2], "Time-varying bias", "bias",
                 "bias_lo", "bias_hi", bias_lim, "Bias of Rt", zero_ref)
  }),
  ncol = 1
) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.title = element_blank(), legend.background = element_blank(),
        legend.key = element_blank(), legend.text = element_text(size = 9),
        legend.margin = margin(0, 0, 2, 0))

pdf_file <- file.path(misspec_dir, sprintf("si_misspec_plots_from_day%d.pdf", draw_from_day))
pdf(pdf_file, width = 11, height = 11.8)
print(page)
invisible(dev.off())

cat("Wrote", pdf_file, "\n")
