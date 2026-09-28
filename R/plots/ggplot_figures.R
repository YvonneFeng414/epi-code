# ==============================================================================
# ggplot_figures.R
# All ggplot figures. Each function returns a plot object; the drivers print and
# ggsave()s them. theme_rt() is the shared theme every figure used to repeat
# inline.
# ==============================================================================

# Shared theme: identical to the theme_bw() + theme(...) block that was
# repeated at the end of every ggplot chain in the original script.
theme_rt <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      strip.background = element_rect(fill = "grey95"),
      panel.grid.minor = element_blank()
    )
}

# Section 4: faceted primary Rt comparison.
# Faceted comparison avoids overlapping uncertainty ribbons.
plot_primary_rt <- function(rt_primary, mean_si, sd_si, evaluation_start_date) {
  ggplot(rt_primary, aes(x = date, y = R)) +
    geom_ribbon(
      aes(ymin = lower, ymax = upper, fill = method),
      alpha = 0.20,
      color = NA,
      show.legend = FALSE
    ) +
    geom_line(aes(color = method), linewidth = 0.7, show.legend = FALSE) +
    geom_hline(yintercept = 1, linetype = 3) +
    facet_wrap(~method, ncol = 1) +
    coord_cartesian(ylim = c(0, 5)) +
    labs(
      title = "Estimated reproduction number in Queens, NY",
      subtitle = paste0(
        "Shared SI: mean = ", mean_si, ", SD = ", sd_si,
        "; evaluation begins ", evaluation_start_date
      ),
      x = "Date",
      y = expression(R[t])
    ) +
    theme_rt()
}

# Sections 5/5B: Rt sensitivity to the assumed SI SD ("sd") or mean ("mean").
plot_si_sensitivity <- function(sens_data, vary = c("sd", "mean"), mean_si, sd_si) {
  vary <- match.arg(vary)

  if (vary == "sd") {
    color_values <- sens_data$assumed_sd_si
    plot_title <- "Sensitivity of estimated Rt to assumed serial-interval SD"
    plot_subtitle <- paste0("Assumed SI mean held fixed at ", mean_si)
    color_label <- "Assumed SI SD"
  } else {
    color_values <- sens_data$assumed_mean_si
    plot_title <- "Sensitivity of estimated Rt to assumed serial-interval mean"
    plot_subtitle <- paste0("Assumed SI SD held fixed at its true value, ", sd_si)
    color_label <- "Assumed SI mean"
  }

  ggplot(sens_data, aes(x = date, y = R, color = factor(color_values))) +
    geom_line(linewidth = 0.65) +
    geom_hline(yintercept = 1, linetype = 3) +
    facet_wrap(~method, ncol = 1) +
    coord_cartesian(ylim = c(0, 5)) +
    labs(
      title = plot_title,
      subtitle = plot_subtitle,
      x = "Date",
      y = expression(R[t]),
      color = color_label
    ) +
    theme_rt()
}

# Section 6: day-level coverage across replicates, one panel per target.
plot_coverage_by_day <- function(rt_coverage_by_day,
                                 truth_label,
                                 epiestim_window,
                                 window_target_label) {
  ggplot(
    rt_coverage_by_day,
    aes(x = date, y = coverage, color = method)
  ) +
    geom_line(linewidth = 0.65) +
    geom_hline(yintercept = 0.95, linetype = 2) +
    facet_wrap(~target, ncol = 1) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      title = expression(paste("Day-level 95% credible interval coverage of ", R[t])),
      subtitle = paste0(
        "Scored against both targets (panels); semi-synthetic replicates from the ",
        truth_label, " full-series fit.\nEach arm is calibrated for one panel ",
        "only: EpiLPS and EpiEstim (1-day window) for 'Instantaneous',\nEpiEstim (",
        epiestim_window, "-day window) for '", window_target_label,
        "'. Dashed line = nominal 0.95."
      ),
      x = "Date",
      y = "Coverage across replicates",
      color = "Method"
    ) +
    theme_rt()
}

# ------------------------------------------------------------------------------
# Fully synthetic simulation study (sim_main.R)
# ------------------------------------------------------------------------------

# The true Rt path of every scenario, with each window's trailing mean overlaid.
# Doubles as the cheap sanity check that a scenario has the intended shape and as
# a direct picture of the estimand offset: wherever the coloured lines separate
# from the black one, a windowed arm and an instantaneous arm are targeting
# genuinely different quantities.
#
# `estimates` is optional so sim_main.R can draw this figure twice: once in
# Section 2 before anything has been fitted (the sanity check), and again after
# fitting with sim_mean_estimate_by_day() output, which turns the estimand offset
# from an assertion about window averaging into something visible - each fitted
# arm should sit on the solid line of ITS OWN window, not on the black one.
#
# With estimates shown the two aesthetics carry different things: COLOUR is the
# estimand (which quantity the line is aiming at) and LINETYPE is the estimator.
# EpiLPS therefore shares the window-1 colour with the 1-day EpiEstim arm, which
# is correct rather than a collision - both target the instantaneous Rt - and the
# two are told apart by dash versus dot.
# `estimate_colors` switches that off. Supply a vector keyed by method_family and
# the estimates layer colours by ESTIMATOR instead of by estimand window. That is
# what a study with only one window wants: there the window colour carries no
# information (every line is window 1) while the method does, and it keeps the
# colour meaning consistent with the other figures in the same output directory.
# NULL keeps the colour-by-window default, which is correct for sim_main.R where
# four genuinely different windows are on screen.
#
# `title` / `subtitle` override the defaults below, which describe sim_main.R's
# fully synthetic scenarios. A caller whose truth comes from somewhere else, or
# which has no windowed arm at all, needs to say so rather than inherit prose
# about a layer its figure does not contain.
plot_scenario_truth <- function(truth_table, estimates = NULL,
                                title = NULL, subtitle = NULL,
                                estimate_colors = NULL) {
  # trailing_mean() is undefined for the first w - 1 days; dropping those rows
  # here keeps ggplot from warning about them on every call.
  windowed <- truth_table %>%
    filter(estimand_window > 1L, is.finite(true_R_window))

  show_estimates <- !is.null(estimates) && nrow(estimates) > 0L

  # Shared numeric ordering for the colour scale, taken from the windows actually
  # drawn. Without it the legend runs 4, 7, 14, 1: window 1 has no trailing-mean
  # line of its own, so it first appears in the estimates layer and gets trained
  # onto the scale last.
  window_levels <- sort(unique(c(
    windowed$estimand_window,
    if (show_estimates) estimates$estimand_window
  )))

  color_by_method <- show_estimates && !is.null(estimate_colors)

  # The black instantaneous truth is drawn FIRST so the fitted lines sit on top of
  # it. The reverse order hides an unbiased arm completely under the thicker black
  # stroke, which reads as "not plotted" rather than "lands exactly on the truth".
  p <- ggplot(truth_table, aes(x = index)) +
    geom_line(aes(y = true_R_instant), color = "black", linewidth = 0.7) +
    geom_line(
      data = windowed,
      aes(y = true_R_window, color = factor(estimand_window, window_levels)),
      linewidth = 0.6
    )

  if (show_estimates) {
    # The estimate table is keyed on the scenario id; the facets are keyed on the
    # label, so carry the ordered label factor over from the truth table rather
    # than making the caller know about it. Window 1 is excluded from the truth
    # overlay above because it coincides with the black line, but its ESTIMATE is
    # a distinct series and is drawn here.
    estimates <- estimates %>%
      left_join(distinct(truth_table, scenario, label), by = "scenario")

    estimate_aes <- if (color_by_method) {
      aes(y = mean_R, color = method_family, linetype = method_family)
    } else {
      aes(y = mean_R, color = factor(estimand_window, window_levels),
          linetype = method_family)
    }

    # When colour also maps to the estimator, the linetype scale must carry the
    # SAME legend name or ggplot draws two separate legends listing identical
    # entries. Matching names merges them into one.
    estimator_legend <- if (color_by_method) {
      "Estimator"
    } else {
      "Estimator\n(mean over replicates)"
    }

    p <- p + geom_line(data = estimates, estimate_aes, linewidth = 0.6) +
      scale_linetype_manual(
        values = c(EpiEstim = 2, EpiLPS = 3),
        name = estimator_legend
      )
  }

  if (is.null(title)) {
    title <- if (show_estimates) {
      expression(paste("True ", R[t], " scenarios, their windowed estimands, and the fitted arms"))
    } else {
      expression(paste("True ", R[t], " scenarios and their windowed estimands"))
    }
  }

  if (is.null(subtitle)) {
    subtitle <- if (show_estimates) {
      paste(
        "Black = instantaneous truth. Solid coloured = trailing mean over a",
        "w-day window, i.e. what a w-day\nsliding-window estimator targets.",
        "Dashed = EpiEstim, dotted = EpiLPS, both averaged over replicates under",
        "\nthe correctly specified SI. Colour is the estimand, so EpiLPS shares",
        "window 1 with the 1-day arm;\neach fitted line should track the solid",
        "line of its own colour."
      )
    } else {
      paste(
        "Black = instantaneous truth (the estimand of EpiLPS and the 1-day arm).",
        "Coloured = trailing mean over a w-day\nwindow, which is what a w-day",
        "sliding-window estimator actually targets."
      )
    }
  }

  # The colour scale depends on what colour was mapped to: estimator names when
  # estimate_colors was supplied, otherwise the ordered window levels.
  p <- p + if (color_by_method) {
    scale_color_manual(values = estimate_colors)
  } else {
    scale_color_discrete(limits = as.character(window_levels))
  }

  p +
    geom_hline(yintercept = 1, linetype = 3) +
    facet_wrap(~label, ncol = 1, scales = "free_y") +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Day",
      y = expression(R[t]),
      color = if (color_by_method) "Estimator" else "Window (days)"
    ) +
    theme_rt()
}

# Day-level coverage, faceted by scenario and target so target mismatch is
# visible as a whole panel sitting below 0.95 rather than as one number.
# `subtitle` overrides the default, which describes sim_main.R's fully synthetic
# replicates and its windowed arms. Callers whose replicates are not synthetic,
# or who have no windowed arm, must supply their own rather than ship a figure
# that misdescribes itself.
#
# `method_colors` is a named vector mapping method to colour. Supply it whenever
# the same arms appear in more than one figure: colour has to follow the entity
# across a set of figures, or a reader carries the wrong mapping from one to the
# next. NULL keeps ggplot's default hue scale.
plot_sim_coverage_by_day <- function(sim_by_day, si_scenario_shown,
                                     subtitle = NULL, method_colors = NULL) {
  if (is.null(subtitle)) {
    subtitle <- paste0(
      "Fully synthetic replicates; assumed SI = '", si_scenario_shown,
      "'. Each arm is calibrated for one column only:\nthe 'Instantaneous' ",
      "column for EpiLPS and the 1-day arm, the window column for the ",
      "matching windowed arm.\nDashed line = nominal 0.95."
    )
  }

  p <- ggplot(
    sim_by_day %>% filter(si_scenario == si_scenario_shown),
    aes(x = index, y = coverage, color = method)
  ) +
    geom_line(linewidth = 0.6) +
    geom_hline(yintercept = 0.95, linetype = 2)

  if (!is.null(method_colors)) {
    p <- p + scale_color_manual(values = method_colors)
  }

  p +
    facet_grid(scenario ~ target) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      title = expression(paste("Day-level 95% credible interval coverage of ", R[t])),
      subtitle = subtitle,
      x = "Day",
      y = "Coverage across replicates",
      color = "Method"
    ) +
    theme_rt()
}

# How own-estimand coverage and bias respond to getting the SI wrong. Assumed SI
# on the x axis, ordered so the misspecification direction reads left to right.
plot_si_misspecification <- function(sim_summary) {
  long <- sim_summary %>%
    select(scenario, si_scenario, method,
           Coverage95_Rt_ownEstimand, Bias_ownEstimand) %>%
    pivot_longer(
      cols = c(Coverage95_Rt_ownEstimand, Bias_ownEstimand),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(
      metric = if_else(
        metric == "Coverage95_Rt_ownEstimand",
        "95% coverage of own estimand",
        "Bias (estimate - own estimand)"
      ),
      si_scenario = factor(si_scenario, levels = unique(sim_summary$si_scenario))
    )

  # One reference line per metric: nominal coverage, and zero bias.
  reference <- data.frame(
    metric = c("95% coverage of own estimand", "Bias (estimate - own estimand)"),
    yintercept = c(0.95, 0)
  )

  ggplot(long, aes(x = si_scenario, y = value, color = method, group = method)) +
    geom_hline(data = reference, aes(yintercept = yintercept), linetype = 2) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1.4) +
    facet_grid(metric ~ scenario, scales = "free_y") +
    labs(
      title = "Effect of serial-interval misspecification on Rt recovery",
      subtitle = paste(
        "Epidemics are always generated with the true SI (mean 7.5, SD 3.4);",
        "only the SI given to the estimators\nis wrong. Dashed lines mark",
        "nominal coverage and zero bias."
      ),
      x = "Assumed serial interval",
      y = NULL,
      color = "Method"
    ) +
    theme_rt() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# ------------------------------------------------------------------------------
# Self-consistency study (truth_source_main.R)
# ------------------------------------------------------------------------------

# Day-level bias against the reported interval half-width.
#
# The overlay figure cannot show this. An unbiased arm's fitted line lands within
# a stroke width of the truth on an axis spanning the whole Rt range, so "no bias"
# and "not plotted" look identical. Plotting the DIFFERENCE puts each arm on its
# own scale, where a bias of 0.007 and one of 0.11 are both legible.
#
# The ribbon is +/- half the mean reported CI width, centred on zero, so the
# figure carries its own reading rule: wherever the line leaves the ribbon, the
# systematic bias alone exceeds the interval half-width and no amount of correct
# noise modelling can cover the truth on that day. That is the mechanism behind a
# coverage number, rather than the number itself.
#
# `data` needs columns: scenario, label (facet), method_family (colour), index,
# bias, half.
plot_arm_residuals <- function(data, method_colors = NULL,
                               title = NULL, subtitle = NULL) {
  if (is.null(title)) {
    title <- expression(paste("Day-level bias against the reported interval half-width"))
  }
  if (is.null(subtitle)) {
    subtitle <- paste(
      "Line = mean estimate minus truth, averaged over replicates. Band = +/-",
      "half the mean reported 95% CI\nwidth. Where the line leaves the band the",
      "bias alone exceeds the interval, so coverage on that day is\nimpossible",
      "however well the noise is modelled."
    )
  }

  p <- ggplot(data, aes(x = index)) +
    geom_ribbon(
      aes(ymin = -half, ymax = half, fill = method_family),
      alpha = 0.18,
      color = NA
    ) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4, color = "grey35") +
    geom_line(aes(y = bias, color = method_family), linewidth = 0.6)

  if (!is.null(method_colors)) {
    p <- p +
      scale_color_manual(values = method_colors) +
      scale_fill_manual(values = method_colors)
  }

  # No legend: each facet holds exactly one estimator and the strip names it, so
  # a legend would just repeat the panel headings. What the line and the band
  # mean is prose, and lives in the subtitle.
  p +
    facet_wrap(~label, ncol = 1, scales = "free_y") +
    guides(color = "none", fill = "none") +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Day",
      y = "Estimate - truth"
    ) +
    theme_rt()
}

# Headline metrics for the two matched arms. The four measures are on different
# units, so they get one panel each with its own scale - never a shared axis and
# never two y-scales on one panel.
#
# `panels` is long, one row per (method, metric), with: method, metric (an
# ordered factor giving panel order), value, lower/upper for error bars (NA where
# an interval is not meaningful), and label (pre-formatted, since each metric
# wants its own precision). Reference lines come from `references`, a data frame
# of metric + yintercept, so only the panels that have a meaningful reference
# get one.
#
# Colours are categorical slots 1 and 2 of the validated palette, assigned in
# fixed order and never cycled. Identity is carried by the x-axis label as well
# as the fill, so it is never colour-alone.
#
# `x_var` names the column on the x axis. It defaults to method, one bar per arm.
# Set it to another grouping (e.g. an SI scheme) to compare the same arms across
# a second factor: bars then dodge within each x group and colour keeps meaning
# the method, rather than spending extra categorical hues on the interaction.
plot_arm_metrics <- function(panels, references = NULL,
                             title = NULL, subtitle = NULL,
                             method_colors = NULL, x_var = "method") {
  if (is.null(method_colors)) {
    method_colors <- stats::setNames(
      c("#2a78d6", "#eb6834")[seq_along(levels(panels$method))],
      levels(panels$method)
    )
  }

  dodge <- position_dodge(width = 0.7)

  p <- ggplot(panels, aes(x = .data[[x_var]], y = value, fill = method))

  if (!is.null(references) && nrow(references) > 0L) {
    p <- p + geom_hline(
      data = references,
      aes(yintercept = yintercept),
      linetype = 2,
      linewidth = 0.4,
      color = "grey35"
    )
  }

  p +
    geom_col(width = 0.62, position = dodge) +
    geom_errorbar(
      aes(ymin = lower, ymax = upper),
      width = 0.14,
      linewidth = 0.4,
      color = "grey20",
      na.rm = TRUE,
      position = dodge
    ) +
    geom_text(
      aes(label = label, y = pmax(value, upper, na.rm = TRUE)),
      vjust = -0.6,
      size = 3.4,
      color = "grey15",
      position = dodge
    ) +
    facet_wrap(~metric, scales = "free_y", nrow = 2) +
    scale_fill_manual(values = method_colors) +
    # Headroom for the value labels, and a baseline anchored at zero.
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL, fill = NULL) +
    theme_rt() +
    theme(
      legend.position = "bottom",
      panel.grid.major.x = element_blank()
    )
}

# Section 7: expanding-window one-step-ahead forecasts vs observed incidence.
plot_forecasts <- function(backtest, forecast_stride) {
  ggplot(backtest, aes(x = target_date)) +
    geom_ribbon(
      aes(ymin = lower, ymax = upper, fill = method),
      alpha = 0.18,
      color = NA,
      show.legend = FALSE
    ) +
    geom_line(aes(y = pred, color = method), linewidth = 0.65,
              show.legend = FALSE) +
    geom_point(aes(y = actual), size = 0.7, alpha = 0.65) +
    facet_wrap(~method, ncol = 1) +
    labs(
      title = "Expanding-window one-step-ahead forecasts",
      subtitle = paste0(
        "Common Poisson renewal forecast; fixed evaluation start; stride = ",
        forecast_stride, " day(s)"
      ),
      x = "Target date",
      y = "Daily reported cases"
    ) +
    theme_rt()
}

# ------------------------------------------------------------------------------
# Shared by truth_source_main.R and compare_si_discretization.R
# ------------------------------------------------------------------------------

# The two plug-in truths side by side, before either arm is simulated from.
# `rt_compare` has columns day, EpiEstim, EpiLPS.
plot_plugin_truths <- function(rt_compare, colors = estimator_colors) {
  ggplot(rt_compare, aes(x = day)) +
    geom_line(aes(y = EpiEstim, color = "EpiEstim"), linewidth = 0.7, alpha = 0.8) +
    geom_line(aes(y = EpiLPS, color = "EpiLPS"), linewidth = 1.1) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
    scale_color_manual(values = colors) +
    labs(
      x = "Day",
      y = expression(R[t]),
      color = "Method",
      title = "Plug-in true Rt paths: EpiEstim vs EpiLPS"
    ) +
    theme_rt()
}

# The share of an arm's MSE that is systematic bias, from the day-level tables:
# mean over scored days of (mean estimate - truth)^2. A credible interval is
# sized for sampling noise, so an arm whose error is mostly bias cannot cover
# its own truth no matter how well its noise model is calibrated. Joined on
# `keys` plus index; returns one row per key combination.
arm_bias_share <- function(by_day, estimates, keys = c("scenario", "method")) {
  inner_join(
    by_day %>% select(all_of(keys), index, target_R),
    estimates %>% select(all_of(keys), index, mean_R),
    by = c(keys, "index")
  ) %>%
    group_by(across(all_of(keys))) %>%
    summarise(mean_bias2 = mean((mean_R - target_R)^2), .groups = "drop")
}

# The four headline metrics plot_arm_metrics() draws, in panel order.
arm_metric_levels <- c(
  "95% CI coverage of own truth",
  "Mean CI width",
  "RMSE",
  "Share of MSE that is systematic bias (%)"
)

# The reference line for the coverage panel only.
arm_metric_references <- function() {
  data.frame(
    metric = factor(arm_metric_levels[1], levels = arm_metric_levels),
    yintercept = 0.95
  )
}

# Long-format panel table for plot_arm_metrics(). `panel_source` has one row
# per arm (x `keys`) with coverage, coverage_mcse, width, rmse and bias_share;
# each metric gets its own precision in `label`, and only coverage carries an
# interval (+/- 1.96 MCSE).
build_arm_metric_panels <- function(panel_source, keys = "method") {
  bind_rows(
    panel_source %>% transmute(
      across(all_of(keys)),
      metric = arm_metric_levels[1], value = coverage,
      lower = coverage - 1.96 * coverage_mcse,
      upper = coverage + 1.96 * coverage_mcse,
      label = sprintf("%.3f", coverage)
    ),
    panel_source %>% transmute(
      across(all_of(keys)),
      metric = arm_metric_levels[2], value = width,
      lower = NA_real_, upper = NA_real_, label = sprintf("%.3f", width)
    ),
    panel_source %>% transmute(
      across(all_of(keys)),
      metric = arm_metric_levels[3], value = rmse,
      lower = NA_real_, upper = NA_real_, label = sprintf("%.4f", rmse)
    ),
    panel_source %>% transmute(
      across(all_of(keys)),
      metric = arm_metric_levels[4], value = bias_share,
      lower = NA_real_, upper = NA_real_, label = sprintf("%.0f%%", bias_share)
    )
  ) %>%
    mutate(metric = factor(metric, levels = arm_metric_levels))
}
