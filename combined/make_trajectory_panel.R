################################################################################
# Estimated Rt over time, and time-varying bias, under SI misspecification.
#
# Format follows the reference figure: estimator as rows, "estimated Rt over
# time" and "time-varying bias" as the two columns, one line per assumed SI with
# a ribbon for the mean 95% interval, panels lettered with tinted title bars.
#
# THREE DELIBERATE DEPARTURES from that reference, all forced by the data:
#
#   1. The EpiEstim row is NOT the seed study's. That study reports summary
#      metrics per grid point only, with no per-day output, so its EpiEstim
#      result cannot be drawn here. combined/run_epiestim_arm.R rebuilds an
#      EpiEstim arm inside this repo instead - same incidence file, same common
#      range, same hole removal, same seeding, scored against true_R_epiestim.
#      Its control point lands at coverage 0.9502 against the seed study's
#      0.9503, so the rebuild reproduces the seed result it replaces.
#   2. +/-60%, not +/-50%. The grid these scripts were run on is the true value
#      scaled by -80%..+80% in 20-point steps, so 50% is not on it. 60% is the
#      nearest arm and is what is drawn and labelled.
#   3. The truth is the real Queens plug-in Rt over days 31-253, not a synthetic
#      step from 0.8 to 1.2 over 100 days. It ranges 0.62 to 1.51 and moves
#      throughout, so there is no single changepoint to mark.
#
# Restricted to the gamma_bin assumed-SI family, as combined_si_misspec_panel is.
#
# Run from queens_rt/:  Rscript combined/make_trajectory_panel.R
################################################################################

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(dplyr)
})

out_dir <- "combined"
fam     <- "gamma_bin"
rel     <- 0.6          # the misspecification level drawn, as a fraction
mean_si <- 7.5
sd_si   <- 3.4
epifilter_arm <- "smoother"

sources <- list(
  list(label = "EpiEstim",      file = "combined/epiestim_arm_by_day.csv",                  arm = NULL),
  list(label = "EpiLPS (MAP)",  file = "results_epilps_si_misspec/si_misspec_by_day.csv",      arm = NULL),
  list(label = "EpiLPS (MALA)", file = "results_epilps_si_misspec_mala/si_misspec_mala_by_day.csv", arm = NULL),
  list(label = "EpiFilter",     file = "results_epifilter_si_misspec/si_misspec_by_day.csv",   arm = epifilter_arm)
)

# The five assumed SIs drawn, addressed by the (mean, sd) each key encodes.
si_arms <- data.frame(
  arm   = c("Correct SI", "Mean +60%", "Mean -60%", "SD +60%", "SD -60%"),
  a_mean = c(mean_si, mean_si * (1 + rel), mean_si * (1 - rel), mean_si, mean_si),
  a_sd   = c(sd_si,   sd_si,               sd_si,               sd_si * (1 + rel), sd_si * (1 - rel)),
  stringsAsFactors = FALSE
)
si_arms$key <- sprintf("%s_%.2f_%.2f", fam, round(si_arms$a_mean, 2), round(si_arms$a_sd, 2))

# "True Rt" is carried as a series rather than an unmapped geom_line, so that it
# appears in the key. NOTE it is a DIFFERENT curve in each row: every estimator
# is scored against its own plug-in truth (EpiLPS rows against true_R_epilps,
# 0.62-1.51; EpiFilter against its own smoother fit, 0.44-1.95). The rows are
# therefore not three views of one ground truth.
arm_levels <- c("True Rt", si_arms$arm)
arm_colours <- c("True Rt" = "#000000",
                 "Correct SI" = "#2471A3", "Mean +60%" = "#C0392B",
                 "Mean -60%"  = "#E8853A", "SD +60%"   = "#1A9E5B",
                 "SD -60%"    = "#7F7F7F")
arm_lty <- c("True Rt" = "solid",
             "Correct SI" = "solid", "Mean +60%" = "longdash", "Mean -60%" = "longdash",
             "SD +60%" = "longdash", "SD -60%" = "longdash")


# ---- load ---------------------------------------------------------------------

grab <- function(src) {
  if (!file.exists(src$file)) stop("Cannot find ", src$file)
  d <- read.csv(src$file, stringsAsFactors = FALSE)
  if (!is.null(src$arm)) d <- d[d$arm == src$arm, , drop = FALSE]
  d <- d[d$key %in% si_arms$key, , drop = FALSE]
  missing <- setdiff(si_arms$key, unique(d$key))
  if (length(missing) > 0L) {
    stop(src$label, " is missing key(s): ", paste(missing, collapse = ", "))
  }
  d$estimator <- src$label
  d$arm <- si_arms$arm[match(d$key, si_arms$key)]
  d[, c("estimator", "arm", "day", "truth", "mean_R", "half", "bias")]
}

dat <- bind_rows(lapply(sources, grab))
dat$estimator <- factor(dat$estimator, levels = vapply(sources, `[[`, "", "label"))
dat$arm       <- factor(dat$arm, levels = arm_levels)

# The truth is one series per estimator, not one per arm.
truth_df <- dat %>%
  group_by(estimator, day) %>%
  summarise(truth = dplyr::first(truth), .groups = "drop")

write.csv(dat, file.path(out_dir, "trajectory_by_day.csv"), row.names = FALSE)

# ---- figure -------------------------------------------------------------------

row_fill <- c("EpiEstim"     = "#E2EFDA", "EpiLPS (MAP)"  = "#FBE4D5",
              "EpiLPS (MALA)" = "#E2D4F0", "EpiFilter"     = "#DEEAF6")

panel_theme <- theme_bw(base_size = 10) +
  theme(panel.grid  = element_blank(),
        panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 0.4),
        axis.title  = element_text(size = 8.5),
        axis.text   = element_text(size = 7.5),
        plot.margin = margin(2, 4, 2, 2))

strip_theme <- function(est) {
  theme(strip.background = element_rect(fill = row_fill[[est]], colour = "grey30",
                                        linewidth = 0.4),
        strip.text = element_text(size = 8.5, face = "bold", hjust = 0,
                                  margin = margin(3, 3, 3, 5)))
}

rt_lim   <- range(c(dat$mean_R - dat$half, dat$mean_R + dat$half, dat$truth))
bias_lim <- range(c(dat$bias - dat$half, dat$bias + dat$half, 0))

make_rt_panel <- function(est, letter) {
  d <- dat[dat$estimator == est, ]
  d$strip <- sprintf("%s. Estimated Rt over time (%s)", letter, est)
  ggplot(d, aes(x = day)) +
    geom_ribbon(aes(ymin = mean_R - half, ymax = mean_R + half, fill = arm),
                alpha = 0.13, colour = NA) +
    geom_line(data = transform(truth_df[truth_df$estimator == est, ],
                               arm = factor("True Rt", levels = arm_levels)),
              aes(y = truth, colour = arm, linetype = arm), linewidth = 0.9) +
    geom_line(aes(y = mean_R, colour = arm, linetype = arm), linewidth = 0.5) +
    facet_wrap(~ strip) +
    scale_colour_manual(values = arm_colours, limits = arm_levels,
                        breaks = arm_levels, drop = FALSE) +
    scale_fill_manual(values = arm_colours, breaks = arm_levels, guide = "none") +
    scale_linetype_manual(values = arm_lty, limits = arm_levels,
                          breaks = arm_levels, drop = FALSE) +
    coord_cartesian(ylim = rt_lim) +
    labs(x = "Day", y = expression(R[t])) +
    panel_theme + strip_theme(est)
}

make_bias_panel <- function(est, letter) {
  d <- dat[dat$estimator == est, ]
  d$strip <- sprintf("%s. Time-varying bias (%s)", letter, est)
  ggplot(d, aes(x = day)) +
    geom_ribbon(aes(ymin = bias - half, ymax = bias + half, fill = arm),
                alpha = 0.13, colour = NA) +
    geom_line(data = data.frame(day = range(d$day), bias = 0,
                                arm = factor("True Rt", levels = arm_levels)),
              aes(y = bias, colour = arm, linetype = arm), linewidth = 0.9) +
    geom_line(aes(y = bias, colour = arm, linetype = arm), linewidth = 0.5) +
    facet_wrap(~ strip) +
    scale_colour_manual(values = arm_colours, limits = arm_levels,
                        breaks = arm_levels, drop = FALSE) +
    scale_fill_manual(values = arm_colours, breaks = arm_levels, guide = "none") +
    scale_linetype_manual(values = arm_lty, limits = arm_levels,
                          breaks = arm_levels, drop = FALSE) +
    coord_cartesian(ylim = bias_lim) +
    labs(x = "Day", y = expression("Bias of " * R[t])) +
    panel_theme + strip_theme(est)
}

ests <- levels(dat$estimator)
lets <- matrix(LETTERS[seq_len(2 * length(ests))], ncol = 2, byrow = TRUE)

rows <- lapply(seq_along(ests), function(i)
  make_rt_panel(ests[i], lets[i, 1]) | make_bias_panel(ests[i], lets[i, 2]))

page <- wrap_plots(rows, ncol = 1) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.title = element_blank(), legend.background = element_blank(),
        legend.key = element_blank(), legend.text = element_text(size = 9),
        legend.margin = margin(0, 0, 2, 0))

pdf(file.path(out_dir, "trajectory_panel.pdf"), width = 11, height = 11.5)
print(page); invisible(dev.off())
png(file.path(out_dir, "trajectory_panel.png"), width = 11, height = 11.5,
    units = "in", res = 200)
print(page); invisible(dev.off())

# ---- readout ------------------------------------------------------------------

cat("\n===== Rt trajectory and time-varying bias (", fam, ") =====\n", sep = "")
cat(sprintf("days %d-%d, misspecification +/-%.0f%%\n",
            min(dat$day), max(dat$day), 100 * rel))
cat("each estimator is scored against its OWN plug-in truth:\n")
for (e in levels(dat$estimator)) {
  t <- dat$truth[dat$estimator == e]
  cat(sprintf("  %-14s %.2f to %.2f\n", e, min(t), max(t)))
}
cat("NOTE: the EpiEstim row is rebuilt by run_epiestim_arm.R, not taken from the seed study.\n\n")

dat %>%
  group_by(estimator, arm) %>%
  summarise(mean_bias = mean(bias), max_abs_bias = max(abs(bias)),
            mean_halfwidth = mean(half), .groups = "drop") %>%
  as.data.frame() %>%
  (function(x) print(data.frame(lapply(x, function(v)
    if (is.numeric(v)) round(v, 4) else v)), row.names = FALSE))

cat(sprintf("\nWrote %s/trajectory_panel.pdf, .png and trajectory_by_day.csv\n", out_dir))
