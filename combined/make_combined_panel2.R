################################################################################
# panel2: score from day 2, drop nothing, let na.rm skip what is missing.
#
# The companion to combined_si_misspec_panel1, which scores all three arms over
# one shared window (days 14-253) chosen as the first day every estimator has
# recovered from the single-day I0 = 2000 seed. This figure asks the opposite
# question: keep every day each estimator will give a number for, starting at
# day 2, and let the averages skip the rest.
#
# READ THE BIAS PANELS WITH THIS IN MIND. na.rm only removes what is actually
# NA, and the three arms do not fail the same way:
#
#   EpiEstim  - genuine NA, from EpiEstim's own `t_end > final_mean_si` rule.
#               The number of missing leading days slides with the assumed SI
#               mean: first estimate on day 4 at mean 3.0, day 8 at 7.5, day 12
#               at 12.0. Correctly skipped.
#   EpiFilter - no estimate before filter_start = 8, emitted as per-day NA.
#               Correctly skipped.
#   EpiLPS    - NOT reached by na.rm. Its B-spline has no grid cap and returns a
#               finite number however little history it has - 5.6e11 on day 2 at
#               assumed sd -60%. Those days are wrong rather than missing.
#
# Left alone, that single day moves EpiLPS's 252-day average bias to 2.2e9 and
# its RMSE to 3.5e10, and the Bias and RMSE panels render as blank space with
# one spike. So a PLAUSIBILITY CEILING is applied instead:
# r_max_plausible = 5, per replicate-day, blanking the whole triple (point and
# both bounds). The rule is identical across all three arms; only EpiLPS ever
# trips it (EpiEstim peaks at 2.9, EpiFilter at 7.3 but only on 0.1% of cells).
#
# WHAT THE CEILING COSTS. At 5 it also removes days where the truth itself
# legitimately exceeds 5 - the plug-in truths peak at 6.48 (EpiLPS) and 6.74
# (EpiEstim) during the explosive phase - so roughly 25 EpiLPS cells that are
# hard rather than absurd go with them. A ceiling of 10 would clear the truth
# entirely and remove only the absurd values, at the cost of leaving some
# residual spikes in the RMSE panels. 5 is the chosen setting; the trade is
# recorded here rather than left implicit.
#
# The other side of the trade is real too, and panel1 exists because of it:
# scoring from each estimator's own first available day means the number of days
# averaged slides with the grid point, and the days that slide in and out are
# the explosive early ones where single-day bias reaches -4.
#
# EpiLPS (MALA) is drawn from epilps_mala_niter40000/, produced by
# run_scorefrom_mala.sh with the same score_from = 2 and r_max_plausible = 5,
# fitted under gamma_bin only, at 40000 iterations / 16000 burn-in. Not the
# script's default 5000 / 2000: that is under-mixed here (median ESS of the Rt
# draws ~17) and its intervals came out ~3% NARROWER than LPSMAP's, which is an
# artefact. At 40000 the MALA interval is 1-2% wider than LPSMAP's at every grid
# point, as sampling the hyperparameters should make it.
# It shares the MAP arm's simulated epidemics, and the ceiling blanks its
# equal-tailed and HPD bounds together; the equal-tailed interval is what is
# scored here, like-for-like with MAP.
#
# Run from queens_rt/:  Rscript combined/make_combined_panel2.R
################################################################################

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(dplyr)
})

out_dir <- "combined"
fam     <- "gamma_bin"

# EpiFilter has two arms; the smoother is the script's stated primary.
epifilter_arm <- "smoother"

# All four arms now run the same +/-60% grid, so the intersection is a no-op;
# the check is kept because it is what catches a source drifting onto a
# different sweep, which has happened before.
restrict_to_common_grid <- TRUE

exp_dir <- "combined/experiments/score-from-2-rmax5"

sources <- list(
  list(label = "EpiEstim",      file = file.path(exp_dir, "epiestim_arm_metrics.csv"),               arm = NULL),
  list(label = "EpiLPS (MAP)",  file = file.path(exp_dir, "epilps_map/si_misspec_metrics.csv"),      arm = NULL),
  list(label = "EpiLPS (MALA)", file = file.path(exp_dir, "epilps_mala_niter40000/si_misspec_mala_metrics.csv"), arm = NULL),
  list(label = "EpiFilter",     file = file.path(exp_dir, "epifilter/si_misspec_metrics.csv"),       arm = epifilter_arm)
)

for (src in sources) {
  if (!file.exists(src$file)) stop("Cannot find ", normalizePath(src$file, mustWork = FALSE))
}

# ---- 1. Harmonise the sources onto one schema --------------------------------
# All four now write the same column names, so this is a filter and a rename
# rather than the per-source mapping the seed CSV used to need.

take <- function(src) {
  d <- read.csv(src$file, stringsAsFactors = FALSE)
  d <- d[grepl(paste0("^", fam, "_"), d$key), , drop = FALSE]
  if (!is.null(src$arm)) d <- d[d$arm == src$arm, , drop = FALSE]
  if (nrow(d) == 0L) stop("No ", fam, " rows in ", src$file)
  data.frame(method = src$label, scenario = d$scenario, rel_error = d$rel_error,
             Bias = d$Bias, RMSE = d$RMSE,
             Coverage = d$Coverage95, CIWidth = d$MeanCIWidth,
             stringsAsFactors = FALSE)
}

combined <- do.call(rbind, lapply(sources, take))

# Report any (method, scenario, rel_error) cell that is absent. EpiFilter loses
# assumed_mean = 12.0 (+60%): a single 2000-case seed day makes day 1 a spike and
# day 2 a trough, and at that assumed SI length the implied Rt reaches 198
# against a grid capped at 10, so the Poisson likelihood underflows and the
# recursion returns NaN. The cell is genuinely missing, not zero - say so rather
# than letting the line quietly skip a point.
full <- expand.grid(method = unique(combined$method),
                    scenario = unique(combined$scenario),
                    rel_error = sort(unique(combined$rel_error)),
                    stringsAsFactors = FALSE)
have <- paste(combined$method, combined$scenario, combined$rel_error)
gaps <- full[!paste(full$method, full$scenario, full$rel_error) %in% have, ]
if (nrow(gaps) > 0L) {
  cat("\nMISSING CELLS (drawn as gaps in the line):\n")
  for (i in seq_len(nrow(gaps))) {
    cat(sprintf("  %-14s %-10s rel_error %+.1f\n",
                gaps$method[i], gaps$scenario[i], gaps$rel_error[i]))
  }
}

if (isTRUE(restrict_to_common_grid)) {
  per_method <- split(combined$rel_error, combined$method)
  common <- Reduce(intersect, lapply(per_method, unique))
  dropped <- setdiff(sort(unique(combined$rel_error)), common)
  if (length(dropped) > 0L) {
    cat(sprintf("\nRestricted to the common grid: dropped %s\n",
                paste0(sprintf("%+.0f%%", 100 * dropped), collapse = ", ")))
  }
  combined <- combined[combined$rel_error %in% common, , drop = FALSE]
}

method_levels <- vapply(sources, `[[`, "", "label")
combined$method   <- factor(combined$method, levels = method_levels)
combined$scenario <- factor(combined$scenario, levels = c("vary_sd", "vary_mean"),
                            labels = c("assumed sd wrong", "assumed mean wrong"))
combined$rel_error_pct <- 100 * combined$rel_error
combined <- combined[order(combined$method, combined$scenario, combined$rel_error), ]

write.csv(combined, file.path(out_dir, "combined_si_misspec_metrics2.csv"),
          row.names = FALSE)

# ---- 2. Figure ---------------------------------------------------------------
# Row order is MEAN first: it is the arm that actually moves the metrics, so it
# leads. The factor built in section 1 is releveled here rather than there so the
# CSV keeps the vary_sd / vary_mean order the source scripts use.
combined$scenario <- factor(as.character(combined$scenario),
                            levels = c("assumed mean wrong", "assumed sd wrong"))

# Legend order, and the order the palettes below are subset to.
key_order <- c("EpiEstim", "EpiFilter", "EpiLPS (MAP)", "EpiLPS (MALA)")

method_colours <- c("EpiEstim"      = "#1f77b4",
                    "EpiFilter"     = "#ff7f0e",
                    "EpiLPS (MAP)"  = "#2ca02c",
                    "EpiLPS (MALA)" = "#9467bd")[key_order]

# Shape as well as colour: in the Bias and RMSE panels the two EpiLPS arms
# coincide to the line width (same data, same posterior-median estimator), and
# colour alone would render that as a single curve.
method_shapes <- c("EpiEstim"      = 16,   # filled circle
                   "EpiFilter"     = 15,   # filled square
                   "EpiLPS (MAP)"  = 17,   # filled triangle
                   "EpiLPS (MALA)" = 18)[key_order]

# MALA dashed, everything else solid. In the Bias and RMSE panels the two EpiLPS
# curves coincide to well under the line width, and a solid curve drawn over a
# solid curve reads as one series with the other missing. Dashing one of them
# lets the MAP line show through the gaps, so the coincidence is visible as
# coincidence. The two separate cleanly in the Coverage panels, where the dash
# simply distinguishes them as any other line style would.
method_linetypes <- c("EpiEstim"      = "solid",
                      "EpiFilter"     = "solid",
                      "EpiLPS (MAP)"  = "solid",
                      "EpiLPS (MALA)" = "longdash")

combined$method <- factor(as.character(combined$method), levels = key_order)

# One tint per arm, carried by both the panel title bars and the row label.
row_fill <- c("assumed mean wrong" = "#FBE4D5",
              "assumed sd wrong"   = "#DEEAF6")
row_label <- c("assumed mean wrong" = "SI mean misspecified (SD fixed)",
               "assumed sd wrong"   = "SI SD misspecified (mean fixed)")

metric_specs <- list(
  list(column = "Bias",     short = "Bias",            y = expression("Bias of " * R[t]), hline = 0),
  list(column = "RMSE",     short = "RMSE",            y = expression("RMSE of " * R[t]), hline = NA),
  list(column = "Coverage", short = "95% CI coverage", y = "95% CI coverage",             hline = 0.95)
)

x_breaks <- sort(unique(combined$rel_error_pct))

arm_title <- c("assumed mean wrong" = "SI mean misspecified",
               "assumed sd wrong"   = "SI SD misspecified")

# Shared limits down each column, padded, so the same metric is on one scale in
# both rows and the two arms can be compared by eye.
col_limits <- lapply(metric_specs, function(m) {
  r <- range(combined[[m$column]], na.rm = TRUE)
  if (!is.na(m$hline)) r <- range(c(r, m$hline))
  r + c(-1, 1) * 0.08 * diff(r)
})

panel_theme <- function(fill) {
  theme_bw(base_size = 10) +
    theme(panel.grid        = element_blank(),
          panel.border      = element_rect(colour = "grey30", fill = NA, linewidth = 0.4),
          plot.title        = element_text(size = 9, face = "bold", hjust = 0,
                                           margin = margin(4, 4, 4, 6)),
          plot.title.position = "panel",
          plot.background   = element_rect(fill = NA, colour = NA),
          axis.title        = element_text(size = 8.5),
          axis.text         = element_text(size = 7.5),
          plot.margin       = margin(2, 4, 2, 2))
}

# The tinted title bar is drawn as a strip rather than a plot title, because a
# ggplot title cannot carry a background fill.
make_panel <- function(m, arm, letter, limits, show_x) {
  d <- combined[combined$scenario == arm, ]
  d$strip <- sprintf("%s. %s (%s)", letter, m$short, arm_title[[arm]])

  p <- ggplot(d, aes(x = rel_error_pct, y = .data[[m$column]],
                     colour = method, shape = method, linetype = method,
                     group = method))
  if (!is.na(m$hline)) {
    p <- p + geom_hline(yintercept = m$hline, linetype = "dashed",
                        colour = "grey45", linewidth = 0.35)
  }
  p +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey45",
               linewidth = 0.35) +
    geom_line(linewidth = 0.55) +
    geom_point(size = 1.5) +
    facet_wrap(~ strip) +
    scale_colour_manual(values = method_colours, breaks = key_order) +
    scale_shape_manual(values = method_shapes, breaks = key_order) +
    scale_linetype_manual(values = method_linetypes, breaks = key_order) +
    scale_x_continuous(breaks = x_breaks) +
    coord_cartesian(ylim = limits) +
    labs(x = if (show_x) "Relative error in the assumed SI parameter (%)" else NULL,
         y = m$y) +
    panel_theme(row_fill[[arm]]) +
    theme(strip.background = element_rect(fill = row_fill[[arm]], colour = "grey30",
                                          linewidth = 0.4),
          strip.text = element_text(size = 8.5, face = "bold", hjust = 0,
                                    margin = margin(3, 3, 3, 5)))
}

# The rotated tinted label down the left of each row.
make_row_label <- function(arm) {
  ggplot() +
    annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,
             fill = row_fill[[arm]], colour = "grey30", linewidth = 0.4) +
    annotate("text", x = 0.5, y = 0.5, label = row_label[[arm]],
             angle = 90, size = 3.0, fontface = "bold") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void() + theme(plot.margin = margin(2, 1, 2, 2))
}

arms    <- levels(combined$scenario)
letters6 <- matrix(c("A","B","C","D","E","F"), nrow = 2, byrow = TRUE)

rows <- lapply(seq_along(arms), function(i) {
  panels <- lapply(seq_along(metric_specs), function(j)
    make_panel(metric_specs[[j]], arms[i], letters6[i, j], col_limits[[j]],
               show_x = TRUE))
  wrap_plots(c(list(make_row_label(arms[i])), panels), nrow = 1,
             widths = c(0.055, 1, 1, 1))
})

page <- wrap_plots(rows, ncol = 1) +
  plot_layout(guides = "collect") &
  theme(legend.position  = "top",
        legend.direction = "horizontal",
        legend.title     = element_blank(),
        legend.background = element_blank(),
        legend.key       = element_blank(),
        legend.text      = element_text(size = 9),
        legend.margin    = margin(0, 0, 2, 0))

pdf(file.path(out_dir, "combined_si_misspec_panel2.pdf"), width = 11.5, height = 6.6)
print(page); invisible(dev.off())

png(file.path(out_dir, "combined_si_misspec_panel2.png"), width = 11.5, height = 6.6,
    units = "in", res = 200)
print(page); invisible(dev.off())

# ---- 3. Console summary ------------------------------------------------------
cat("\n===== Combined SI misspecification panel (", fam, ") =====\n", sep = "")
cat(sprintf("EpiFilter arm drawn: %s\n", epifilter_arm))
# The four sources are not guaranteed to share a grid - the seed study is
# maintained outside this repo and has been regenerated on a shorter sweep
# before. Report what each method actually covers rather than a product that
# need not multiply out, and say so loudly when they disagree.
grid_by_method <- tapply(combined$rel_error, combined$method,
                         function(x) sort(unique(x)))
npts <- vapply(grid_by_method, length, integer(1))

cat(sprintf("rows: %d\n", nrow(combined)))
for (m in names(grid_by_method)) {
  cat(sprintf("  %-14s %d grid points, %+.0f%% to %+.0f%%\n", m, npts[[m]],
              100 * min(grid_by_method[[m]]), 100 * max(grid_by_method[[m]])))
}
if (length(unique(npts)) > 1L) {
  cat("\n  NOTE: the sources do not share a grid. Curves are drawn over each\n")
  cat("        method's own range, so the panels are not equally wide in x for\n")
  cat("        every method. Check the source CSVs before comparing extremes.\n")
}
cat("\n")

ctrl <- combined[combined$rel_error == 0, ]
cat("Correctly specified point (rel_error = 0):\n")
print(data.frame(method = ctrl$method, scenario = ctrl$scenario,
                 Coverage = round(ctrl$Coverage, 4), Bias = round(ctrl$Bias, 4),
                 RMSE = round(ctrl$RMSE, 4), CIWidth = round(ctrl$CIWidth, 4)),
      row.names = FALSE)

cat("\nCoverage range over the -80%..+80% grid:\n")
combined %>%
  group_by(method, scenario) %>%
  summarise(min = min(Coverage), max = max(Coverage), .groups = "drop") %>%
  mutate(range = max - min) %>%
  as.data.frame() %>%
  (function(x) print(data.frame(lapply(x, function(v)
    if (is.numeric(v)) round(v, 3) else v)), row.names = FALSE))

cat(sprintf("\nWrote %s/\n", out_dir))
