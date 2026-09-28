# ==============================================================================
# sim_epiestim_figures.R
# ggplot pages for the pure-synthetic EpiEstim studies (simulations/). Every
# function returns a patchwork or ggplot object; the drivers print them into an
# open pdf() device (multi-page files) or hand them to ggsave() (single pages).
# Ported from the archive's sensitivity_plot_functions.R,
# rt_trajectory_functions.R, rt_misspec_functions.R and the figure sections of
# rt_bias_convergence_main.R and britton_validation.R, retargeted at the
# harness's column names.
#
# All figure text is plain ASCII: the pdf() device's default fonts have no
# subscript or arrow glyphs, so "Rt" and "->" rather than their Unicode forms.
#
# Requires ggplot2 and patchwork attached; core/palette.R for scenario_colors
# and si_family_colors; plots/si_misspec_figures.R for sweep_legend and
# rel_error_xlab; plots/ggplot_figures.R for theme_rt().
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Theme, page frame, relative-error scale
# ------------------------------------------------------------------------------

# Base-R look: white panel, full box, no grid, the key INSIDE the panel with no
# title. `legend` / `legend_just` place it (ggplot2 >= 3.5 syntax).
sim_theme <- function(legend = c(0.97, 0.97), legend_just = c(1, 1)) {
  theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
      legend.position = "inside",
      legend.position.inside = legend,
      legend.justification = legend_just,
      legend.title = element_blank(),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.key.size = unit(0.9, "lines"),
      legend.text = element_text(size = 8),
      legend.margin = margin(0, 0, 0, 0),
      axis.title = element_text(size = 9)
    )
}

# The same look with the key OUTSIDE, for pages whose panels share one colour
# scale and let patchwork collect a single legend at the bottom.
sim_theme_bottom <- function() {
  theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
      legend.position = "bottom",
      legend.title = element_text(size = 9),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.key.size = unit(0.9, "lines"),
      legend.text = element_text(size = 8.5),
      axis.title = element_text(size = 9)
    )
}

# Panels on one page under a bold title and a small subtitle. `collect` pools
# the panels' legends into one at the bottom.
sim_page <- function(panels, nrow, title = NULL, subtitle = NULL, collect = FALSE) {
  page <- wrap_plots(panels, nrow = nrow)
  if (collect) page <- page + plot_layout(guides = "collect")
  page <- page + plot_annotation(
    title = title, subtitle = subtitle,
    theme = theme(plot.title = element_text(face = "bold", size = 12),
                  plot.subtitle = element_text(size = 8.5))
  )
  if (collect) page <- page & theme(legend.position = "bottom")
  page
}

# "+40%" / "-80%" labels, as an ordered factor so legends run from the most
# negative to the most positive error rather than alphabetically.
rel_error_labels <- function(errors) {
  paste0(ifelse(errors > 0, "+", ""), round(100 * errors), "%")
}

rel_error_factor <- function(errors, all_errors) {
  factor(rel_error_labels(errors), levels = rel_error_labels(sort(unique(all_errors))))
}

# Diverging palette keyed by label: blue for an SI assumed too SMALL, red for
# one assumed too LARGE, grey for the correct one so it reads as the reference.
# Deliberately not scenario_colors - there blue/orange separate the two SWEEPS;
# here the sweeps are separate panels and colour carries the size of the error.
rel_error_palette <- function(all_errors) {
  errors <- sort(unique(all_errors))
  cols <- character(length(errors))
  neg <- errors < 0
  pos <- errors > 0
  cols[errors == 0] <- "#4D4D4D"
  if (any(neg)) cols[neg] <- grDevices::colorRampPalette(c("#08519C", "#9ECAE1"))(sum(neg))
  if (any(pos)) cols[pos] <- grDevices::colorRampPalette(c("#FCAE91", "#A50F15"))(sum(pos))
  names(cols) <- rel_error_labels(errors)
  cols
}

# ------------------------------------------------------------------------------
# 2. Scenario shapes
# ------------------------------------------------------------------------------

# Every Rt shape on one sheet: the instantaneous truth in black and, where a
# windowed arm exists, its trailing-mean estimand in the vary_sd blue. Input is
# scenario_truth_table(); facets read in registry order.
plot_scenario_shapes <- function(truth_table, window = 7L, ncol = 3L) {
  first_window <- min(truth_table$estimand_window)
  inst <- truth_table[truth_table$estimand_window == first_window, ]
  win <- truth_table[truth_table$estimand_window == window &
                       is.finite(truth_table$true_R_window), ]

  p <- ggplot() +
    geom_hline(yintercept = 1, linetype = 3, colour = "grey50") +
    geom_line(data = inst, aes(x = index, y = true_R_instant),
              colour = "black", linewidth = 0.7)
  if (nrow(win) > 0L && window > 1L) {
    p <- p + geom_line(data = win, aes(x = index, y = true_R_window),
                       colour = scenario_colors[["vary_sd"]], linewidth = 0.5)
  }
  p +
    facet_wrap(~label, ncol = ncol, scales = "free_y") +
    labs(
      title = "True Rt shapes",
      subtitle = if (window > 1L) {
        sprintf("Black = instantaneous Rt (the 1-day arm's estimand); blue = trailing %d-day mean (the %d-day arm's estimand).",
                window, window)
      } else {
        "Black = instantaneous Rt."
      },
      x = "day", y = "Rt"
    ) +
    theme_rt(base_size = 9)
}

# ------------------------------------------------------------------------------
# 3. Moment sweeps and family pages
# ------------------------------------------------------------------------------

# The three metric panels of a misspecification page, in figure order, with
# where each one's inside key goes: coverage descends so the key sits bottom
# right, RMSE rises so top right, width rises from the left so bottom left.
sim_metric_specs <- list(
  list(col = "Coverage95",  title = "Coverage", ylab = "95% CI coverage",
       hline = 0.95, legend = c(0.97, 0.03), just = c(1, 0)),
  list(col = "RMSE",        title = "RMSE",     ylab = "RMSE of Rt",
       hline = NA,   legend = c(0.97, 0.98), just = c(1, 1)),
  list(col = "MeanCIWidth", title = "CI width", ylab = "mean CI width",
       hline = NA,   legend = c(0.03, 0.03), just = c(0, 0))
)

# One metric against relative error, both sweeps overlaid: blue = assumed sd
# wrong, orange = assumed mean wrong. Dashed = nominal (coverage only), dotted
# = the correctly specified point.
make_sweep_panel <- function(d, spec, x_breaks) {
  p <- ggplot(d, aes(x = 100 * rel_error, y = .data[[spec$col]], colour = scenario))
  if (!is.na(spec$hline)) {
    p <- p + geom_hline(yintercept = spec$hline, linetype = "dashed",
                        colour = "grey30", linewidth = 0.4)
  }
  p +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey30", linewidth = 0.4) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1.8, shape = 16) +
    scale_colour_manual(values = scenario_colors, labels = sweep_legend,
                        breaks = names(sweep_legend)) +
    scale_x_continuous(breaks = x_breaks) +
    labs(title = spec$title, x = rel_error_xlab, y = spec$ylab) +
    sim_theme(spec$legend, spec$just)
}

# One page for one target, one family, one arm: coverage / RMSE / width against
# relative error, both sweeps in each panel. `metrics` is already restricted.
plot_misspec_page <- function(metrics, relative_error, title = NULL, subtitle = NULL) {
  panels <- lapply(sim_metric_specs, function(spec) {
    make_sweep_panel(metrics, spec, 100 * relative_error)
  })
  sim_page(panels, nrow = 1, title, subtitle)
}

# One metric for one sweep, one line per assumed family; open symbols mark
# settings whose realised moments drifted past drift_tol (the x-axis is the
# REQUESTED moment, which overstates the error there).
make_family_panel <- function(d, spec, x_breaks, family_colors, sweep_title) {
  d$drift_ok <- factor(d$drift_ok, levels = c(TRUE, FALSE))
  p <- ggplot(d, aes(x = 100 * rel_error, y = .data[[spec$col]], colour = family))
  if (!is.na(spec$hline)) {
    p <- p + geom_hline(yintercept = spec$hline, linetype = "dashed",
                        colour = "grey30", linewidth = 0.4)
  }
  p +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey30", linewidth = 0.4) +
    geom_line(linewidth = 0.55) +
    geom_point(aes(shape = drift_ok), size = 1.8) +
    scale_colour_manual(values = family_colors, name = "assumed SI family",
                        breaks = names(family_colors), drop = FALSE) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), drop = FALSE,
                       labels = c(`TRUE` = "moments as requested",
                                  `FALSE` = "realised moments drifted"),
                       name = NULL) +
    scale_x_continuous(breaks = x_breaks) +
    labs(title = paste0(sweep_title, ": ", spec$title), x = rel_error_xlab, y = spec$ylab) +
    sim_theme_bottom()
}

# One page for one target and one arm: two rows (sd sweep, mean sweep) by the
# three metrics, colour = family, one collected legend. The family x moment
# cross the archive spread over two appendices.
plot_family_page <- function(metrics, relative_error, family_colors,
                             title = NULL, subtitle = NULL) {
  families <- intersect(names(family_colors), unique(metrics$family))
  colors <- family_colors[families]
  panels <- list()
  for (sc in c("vary_sd", "vary_mean")) {
    d <- metrics[metrics$scenario == sc, ]
    for (spec in sim_metric_specs) {
      panels[[length(panels) + 1L]] <-
        make_family_panel(d, spec, 100 * relative_error, colors, sweep_short[[sc]])
    }
  }
  sim_page(panels, nrow = 2, title, subtitle, collect = TRUE)
}

# ------------------------------------------------------------------------------
# 4. Rt trajectories under the correct SI
# ------------------------------------------------------------------------------

# The estimate reuses the sweep blue; the truth is black so it reads as the
# reference; the Monte Carlo spread is the sweep orange.
rt_colours <- c("True Rt" = "#000000", "Estimate" = "#2a78d6", "MC spread" = "#eb6834")

# A common y range for the panels of a page. The early windows have little data
# behind them and produce enormous credible intervals, which left alone flatten
# everything else; clip to the 0.5 / 99.5 percentiles of the bounds instead,
# and let the range follow the data rather than start at zero. Applied with
# coord_cartesian(), which clips the VIEW - ylim() would drop rows before the
# ribbon is drawn and silently deform it.
rt_y_limits <- function(lowers, uppers, true_R) {
  y_lo <- min(stats::quantile(unlist(lowers), 0.005, na.rm = TRUE),
              min(true_R, na.rm = TRUE), na.rm = TRUE)
  y_hi <- max(stats::quantile(unlist(uppers), 0.995, na.rm = TRUE),
              max(true_R, na.rm = TRUE), na.rm = TRUE)
  pad <- 0.05 * (y_hi - y_lo)
  c(max(0, y_lo - pad), y_hi + pad)
}

rt_panel_frame <- function(p, title, y_limits, legend = c(0.97, 0.97),
                           legend_just = c(1, 1)) {
  p +
    scale_colour_manual(values = rt_colours, breaks = names(rt_colours)) +
    # The ribbon shares the estimate's colour; a fill key would repeat it.
    scale_fill_manual(values = rt_colours, breaks = names(rt_colours), guide = "none") +
    coord_cartesian(ylim = y_limits) +
    labs(title = title, x = "day", y = "Rt") +
    sim_theme(legend, legend_just)
}

# Panel (a): one representative epidemic - the posterior median with its 95%
# credible interval, and the truth drawn last so the ribbon never hides it.
plot_rt_single <- function(single, y_limits, title = "One simulated epidemic") {
  p <- ggplot(single, aes(x = day)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = "Estimate"), alpha = 0.25) +
    geom_line(aes(y = est, colour = "Estimate"), linewidth = 0.6) +
    geom_line(aes(y = true_R, colour = "True Rt"), linewidth = 0.7)
  rt_panel_frame(p, title, y_limits)
}

# Panel (b): averaged over every epidemic. The gap between the blue and black
# lines is the bias; the ribbon width is the CI width; the ribbon against the
# dashed orange lines is the coverage - a ribbon narrower than the spread of
# the estimates means the posterior is overconfident.
plot_rt_montecarlo <- function(summary, y_limits, n_sim_used, title = NULL) {
  if (is.null(title)) title <- paste0("Averaged over ", n_sim_used, " simulated epidemics")
  p <- ggplot(summary, aes(x = day)) +
    geom_ribbon(aes(ymin = mean_lower, ymax = mean_upper, fill = "Estimate"), alpha = 0.25) +
    geom_line(aes(y = mc_lower, colour = "MC spread"), linetype = "dashed", linewidth = 0.4) +
    geom_line(aes(y = mc_upper, colour = "MC spread"), linetype = "dashed", linewidth = 0.4) +
    geom_line(aes(y = mean_est, colour = "Estimate"), linewidth = 0.6) +
    geom_line(aes(y = true_R, colour = "True Rt"), linewidth = 0.7)
  rt_panel_frame(p, title, y_limits)
}

# One page: the two panels side by side, sharing a y range. The first days
# after a seed spike fall under EpiEstim's posterior-CV threshold and come back
# NA ("too early in the epidemic"); they are dropped here rather than left for
# ggplot to warn about on every line.
plot_trajectory_page <- function(single, summary, n_sim_used, title = NULL, subtitle = NULL) {
  single <- single[!is.na(single$est), ]
  summary <- summary[!is.na(summary$mean_est), ]
  y_limits <- rt_y_limits(
    lowers = list(single$lower, summary$mean_lower, summary$mc_lower),
    uppers = list(single$upper, summary$mean_upper, summary$mc_upper),
    true_R = summary$true_R
  )
  panels <- list(plot_rt_single(single, y_limits),
                 plot_rt_montecarlo(summary, y_limits, n_sim_used))
  sim_page(panels, nrow = 1, title, subtitle)
}

# ------------------------------------------------------------------------------
# 5. Rt trajectories under a wrong SI
# ------------------------------------------------------------------------------

# Top row: the mean estimate per grid point over the truth in black. The gap
# between a coloured line and the black one IS the misspecification bias, day
# by day. `arm_df` is one sweep's rows with columns day, true_R, mean_est,
# error_label (a rel_error_factor()).
plot_misspec_estimates <- function(arm_df, palette, y_limits, title) {
  truth_df <- unique(arm_df[, c("day", "true_R")])
  ggplot() +
    geom_line(data = truth_df, aes(x = day, y = true_R), colour = "black", linewidth = 0.9) +
    geom_line(data = arm_df, aes(x = day, y = mean_est, colour = error_label), linewidth = 0.55) +
    scale_colour_manual(values = palette, name = "error in assumed SI moment", drop = FALSE) +
    coord_cartesian(ylim = y_limits) +
    labs(title = title, x = "day", y = "Rt") +
    sim_theme_bottom()
}

# Bottom row: per-day coverage. The scalar in a metrics table is one number for
# the whole series and hides that coverage collapses on some days and is
# untouched on others.
plot_misspec_coverage <- function(arm_df, palette, title) {
  ggplot(arm_df) +
    geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey30", linewidth = 0.4) +
    geom_line(aes(x = day, y = coverage, colour = error_label), linewidth = 0.55) +
    scale_colour_manual(values = palette, name = "error in assumed SI moment", drop = FALSE) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = title, x = "day", y = "95% CrI coverage") +
    sim_theme_bottom()
}

# One page, 2 x 2: estimates (top) and coverage (bottom), mean sweep left and
# sd sweep right, one shared y range on the top row and one collected legend.
# `traj` has a `scenario` column (vary_mean / vary_sd) plus the columns above.
# Days EpiEstim declined to estimate (NA, see plot_trajectory_page) are dropped.
plot_misspec_trajectory_page <- function(traj, palette, title = NULL, subtitle = NULL) {
  traj <- traj[!is.na(traj$mean_est), ]
  mean_arm <- traj[traj$scenario == "vary_mean", ]
  sd_arm <- traj[traj$scenario == "vary_sd", ]
  y_limits <- rt_y_limits(lowers = list(traj$mean_est), uppers = list(traj$mean_est),
                          true_R = traj$true_R)
  panels <- list(
    plot_misspec_estimates(mean_arm, palette, y_limits, "Assumed SI MEAN wrong (sd correct)"),
    plot_misspec_estimates(sd_arm, palette, y_limits, "Assumed SI SD wrong (mean correct)"),
    plot_misspec_coverage(mean_arm, palette, "Coverage, assumed mean wrong"),
    plot_misspec_coverage(sd_arm, palette, "Coverage, assumed sd wrong")
  )
  sim_page(panels, nrow = 2, title, subtitle, collect = TRUE)
}

# ------------------------------------------------------------------------------
# 6. Per-day bias convergence
# ------------------------------------------------------------------------------

padded_range <- function(x, pad = 0.08) {
  rng <- range(x, na.rm = TRUE)
  rng + c(-1, 1) * pad * diff(rng)
}

# One panel of per-day bias, one line per assumed sd, with optional dotted
# Euler-Lotka asymptotes in matching colours (`el_df`: el_bias, colour).
make_bias_panel <- function(df, y_col, title, y_limits, y_label, palette, el_df = NULL) {
  p <- ggplot(df, aes(x = day, y = .data[[y_col]], colour = error_label)) +
    geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.4)
  if (!is.null(el_df)) {
    # geom_hline never inherits aesthetics, so the colour is passed outside aes().
    p <- p + geom_hline(data = el_df, aes(yintercept = el_bias),
                        colour = el_df$colour, linetype = "dotted", linewidth = 0.5)
  }
  p +
    geom_line(linewidth = 0.55) +
    scale_colour_manual(values = palette, name = "error in assumed SI sd", drop = FALSE) +
    coord_cartesian(ylim = y_limits) +
    labs(title = title, x = "day", y = y_label) +
    sim_theme_bottom()
}

# Three panels (two if no correctly specified sd is in the grid): the full
# series showing the start-of-series spike decay; from `zoom_from` the raw bias
# against the dotted asymptotes; and the same minus the correctly specified
# curve, which cancels the Monte Carlo error common to all fits and is the panel
# that settles whether the bias lands on Euler-Lotka. `plot_df` has day, bias,
# bias_rel, error_label; `el_df` has el_bias, colour.
plot_bias_convergence_page <- function(plot_df, el_df, palette, zoom_from = 30L,
                                       title = NULL, subtitle = NULL) {
  plot_df <- plot_df[!is.na(plot_df$bias), ]
  zoom_df <- plot_df[plot_df$day >= zoom_from, ]
  raw_label <- "bias in Rt  (mean estimate - true)"
  has_ref <- any(is.finite(plot_df$bias_rel))

  panels <- list(
    make_bias_panel(plot_df, "bias", "Full series: the start-of-series artefact",
                    range(plot_df$bias, na.rm = TRUE), raw_label, palette),
    make_bias_panel(zoom_df, "bias", paste0("From day ", zoom_from, ": raw bias"),
                    padded_range(c(zoom_df$bias, el_df$el_bias)), raw_label, palette, el_df)
  )
  if (has_ref) {
    panels[[3L]] <- make_bias_panel(
      zoom_df, "bias_rel",
      paste0("From day ", zoom_from, ": minus the correctly-specified curve"),
      padded_range(c(zoom_df$bias_rel, el_df$el_bias)),
      "bias relative to correct sd", palette, el_df
    )
  }
  sim_page(panels, nrow = 1, title, subtitle, collect = TRUE)
}

# ------------------------------------------------------------------------------
# 7. Britton & Scalia Tomba validation
# ------------------------------------------------------------------------------

britton_series_colours <- c(EL = "#2a78d6", EPI = "#eb6834")
britton_series_labels <- c(EL = "Euler-Lotka  1 / sum w(s) exp(-rs)",
                           EPI = "EpiEstim R-hat (1-day window)")

# The theoretical curve over the fine sd grid as a line, EpiEstim at the marker
# grid as points; dashed = the truth (or zero bias), dotted = the correct sd.
make_britton_panel <- function(df, df_marks, y_col, mark_col, y_label, hline,
                               sd_true, title, legend, legend_just) {
  ggplot() +
    geom_hline(yintercept = hline, linetype = "dashed", colour = "grey30", linewidth = 0.4) +
    geom_vline(xintercept = sd_true, linetype = "dotted", colour = "grey30", linewidth = 0.4) +
    geom_line(data = df, aes(x = assumed_sd_si, y = .data[[y_col]], colour = "EL"),
              linewidth = 0.6) +
    geom_point(data = df_marks, aes(x = assumed_sd_si, y = .data[[mark_col]], colour = "EPI"),
               size = 2, shape = 16) +
    scale_colour_manual(values = britton_series_colours, labels = britton_series_labels,
                        breaks = names(britton_series_labels)) +
    labs(title = title, x = sprintf("assumed SI sd (true sd = %.1f)", sd_true), y = y_label) +
    sim_theme(legend, legend_just)
}

# 2 x 2: columns = the true R values, row 1 the estimate against the theory,
# row 2 the same as relative bias. `results` is the validation table.
plot_britton_page <- function(results, R_values, sd_true, title = NULL, subtitle = NULL) {
  panels <- vector("list", 2L * length(R_values))
  for (i in seq_along(R_values)) {
    R_true <- R_values[i]
    df <- results[results$true_R == R_true, ]
    df$EL_bias_pct <- 100 * (df$R_EL_discrete - R_true) / R_true
    df_marks <- df[df$is_grid_point, ]

    panels[[i]] <- make_britton_panel(
      df, df_marks, "R_EL_discrete", "R_hat_nd1", "estimated R", R_true, sd_true,
      sprintf("True R = %.1f   (r = %.4f, doubling %.1f d)",
              R_true, df$growth_rate[1L], df$doubling_time[1L]),
      legend = c(0.97, 0.97), legend_just = c(1, 1)
    )
    # Bottom-left key: the bias curve only ever descends, so a top-right key
    # would sit on the zero line.
    panels[[i + length(R_values)]] <- make_britton_panel(
      df, df_marks, "EL_bias_pct", "bias_rel_pct", "relative bias in R (%)", 0, sd_true,
      sprintf("Relative bias, true R = %.1f", R_true),
      legend = c(0.03, 0.03), legend_just = c(0, 0)
    )
  }
  sim_page(panels, nrow = 2, title, subtitle)
}

britton_recon_colours <- c(ALL = "#eb6834", LATE = "#2a78d6", EL = "#1a9e5b", STO = "#7F7F7F")

# The reconciliation: the deterministic replica of the stochastic set-up scored
# over the study's windows (ALL), the same restricted to late windows (LATE),
# the Euler-Lotka prediction (EL), and - when the stochastic study has run -
# its published bias as open grey circles (STO).
plot_britton_reconciliation <- function(recon, sd_true, late_from, scored_from,
                                        title = NULL, subtitle = NULL) {
  long <- rbind(
    data.frame(assumed_sd_si = recon$assumed_sd_si, bias = recon$bias_scored_windows, series = "ALL"),
    data.frame(assumed_sd_si = recon$assumed_sd_si, bias = recon$bias_late_windows, series = "LATE"),
    data.frame(assumed_sd_si = recon$assumed_sd_si, bias = recon$bias_euler_lotka, series = "EL")
  )
  labels <- c(ALL = sprintf("windows from day %d (as the stochastic study scores)", scored_from),
              LATE = paste0("windows after day ", late_from),
              EL = "Euler-Lotka prediction",
              STO = "stochastic study (sim_si_misspec)")

  p <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30", linewidth = 0.4) +
    geom_vline(xintercept = sd_true, linetype = "dotted", colour = "grey30", linewidth = 0.4) +
    geom_line(data = long, aes(x = assumed_sd_si, y = bias, colour = series), linewidth = 0.6) +
    geom_point(data = long, aes(x = assumed_sd_si, y = bias, colour = series), size = 1.8, shape = 16)

  if (!all(is.na(recon$bias_stochastic))) {
    sto <- data.frame(assumed_sd_si = recon$assumed_sd_si, bias = recon$bias_stochastic,
                      series = "STO")
    p <- p + geom_point(data = sto, aes(x = assumed_sd_si, y = bias, colour = series),
                        size = 3.2, shape = 1, stroke = 0.9)
  }

  p <- p +
    scale_colour_manual(values = britton_recon_colours, labels = labels,
                        breaks = names(labels)) +
    labs(x = sprintf("assumed SI sd (true sd = %.1f)", sd_true), y = "bias in R") +
    sim_theme(c(0.03, 0.97), c(0, 1))

  p + plot_annotation(
    title = title, subtitle = subtitle,
    theme = theme(plot.title = element_text(face = "bold", size = 12),
                  plot.subtitle = element_text(size = 8.5))
  )
}
