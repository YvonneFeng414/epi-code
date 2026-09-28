# ==============================================================================
# scenarios.R
# Scenario definitions for the fully synthetic simulation study (sim_main.R).
#
# Two independent scenario dimensions:
#
#   rt_scenarios - the TRUE Rt trajectory that generates the epidemics. Because
#                  it is written down analytically, the truth is known exactly
#                  and does not depend on any fit. This is what separates the
#                  fully synthetic study from the semi-synthetic coverage study
#                  in main.R Section 6, where the "truth" is a plug-in estimate.
#
#   si_scenarios - the serial interval handed to the ESTIMATORS. Epidemics are
#                  always generated with the true SI, so this dimension isolates
#                  the effect of SI misspecification on Rt recovery.
# ==============================================================================

# ------------------------------------------------------------------------------
# Rt trajectory constructors
# ------------------------------------------------------------------------------
# Each constructor returns a function of `time` that yields a length-`time`
# numeric Rt path. Deferring `time` keeps the registry below independent of
# sim_time, so changing the series length needs no edits to the scenarios.

# Constant Rt. Window averaging introduces no offset here, which is why this is
# the only scenario directly comparable with epi.R (see sim_main.R).
rt_constant <- function(value) {
  function(time) rep(value, time)
}

# Abrupt change on `change_day` (e.g. an intervention taking effect).
# The worst case for a sliding window: the window smears the jump across w days.
rt_step <- function(from, to, change_day) {
  function(time) {
    change_day <- as.integer(change_day)
    if (change_day < 2L || change_day > time) {
      stop("change_day must lie inside the simulated period.")
    }
    c(rep(from, change_day - 1L), rep(to, time - change_day + 1L))
  }
}

# Gradual change: `from` up to change_day - 1, then a linear ramp of
# `transition_days` values that ends AT `to`, then `to`. The ramp's first value
# is one step past `from`, not `from` itself (the archive's seq(...)[-1]), so a
# 7-day transition is exactly seven days at intermediate levels. A step is the
# transition_days = 1 case.
rt_gradual <- function(from, to, change_day, transition_days) {
  function(time) {
    change_day <- as.integer(change_day)
    transition_days <- as.integer(transition_days)
    if (change_day < 2L || transition_days < 1L ||
        change_day + transition_days - 1L > time) {
      stop("The transition must lie inside the simulated period.")
    }
    ramp <- seq(from, to, length.out = transition_days + 1L)[-1L]
    c(rep(from, change_day - 1L), ramp,
      rep(to, time - change_day - transition_days + 1L))
  }
}

# Linear trend across the whole period. Mirrors the shape the real Queens series
# shows in main.R Section 6, where the plug-in truth trends from 0.62 to 1.51.
rt_ramp <- function(from, to) {
  function(time) seq(from, to, length.out = time)
}

# Sinusoidal Rt. Tests whether a smoother can follow the rate of change; EpiLPS
# uses a K-dimensional B-spline basis, so a period short relative to the series
# length is the case where its basis may be too coarse.
rt_sine <- function(centre, amplitude, period) {
  function(time) centre + amplitude * sin(2 * pi * seq_len(time) / period)
}

# ------------------------------------------------------------------------------
# Rt scenario registry
# ------------------------------------------------------------------------------
# The list names become the `scenario` column in every output file, so they are
# kept short and sortable. `label` is the human-readable form used in plots.
rt_scenarios <- list(
  S1_constant = list(
    label = "Constant Rt = 1",
    build = rt_constant(1.0)
  ),
  S2_step = list(
    label = "Step 1.5 -> 0.7 at day 50",
    build = rt_step(1.5, 0.7, 50L)
  ),
  S3_ramp = list(
    label = "Linear ramp 0.7 -> 1.5",
    build = rt_ramp(0.7, 1.5)
  ),
  S4_sine = list(
    label = "Sine: centre 1.0, amplitude 0.3, period 30",
    build = rt_sine(1.0, 0.3, 30)
  )
)

# ------------------------------------------------------------------------------
# SI misspecification registry
# ------------------------------------------------------------------------------
# Assumed (mean, sd) used for FITTING. "correct" matches the SI used to generate
# the epidemics; the others are deliberately wrong by a moderate-to-large margin
# (mean shifted -33%/+33%, sd shifted -56%/+76%). These are the same five
# scenarios as si_misspecification.R, kept identical so the two are comparable.
si_scenarios <- list(
  correct       = c(mean = 7.5, sd = 3.4),
  mean_too_low  = c(mean = 5.0, sd = 3.4),
  mean_too_high = c(mean = 10.0, sd = 3.4),
  sd_too_low    = c(mean = 7.5, sd = 1.5),
  sd_too_high   = c(mean = 7.5, sd = 6.0)
)

# ------------------------------------------------------------------------------
# The EpiEstim simulation study's 27 sub-scenarios (simulations/)
# ------------------------------------------------------------------------------
# The archive's build_scenarios() as a table: three constant levels, six
# step pairs, and the same six pairs as gradual transitions of 7, 14 and 20
# days. Row order is the archive's, which is the order the appendix table
# reads in. `true_rt` is the display label, in ASCII ("->", and "1" not "1.0"
# for the constant level, both as the archive wrote them) because the pdf()
# device's default fonts have no arrow glyph.
epiestim_scenario_table <- local({
  pairs <- data.frame(r_from = c(1.2, 0.8, 1.3, 1.1, 0.9, 0.7),
                      r_to   = c(0.8, 1.2, 1.1, 1.3, 0.7, 0.9))
  constant <- data.frame(shape = "constant", r_from = c(1.0, 0.8, 1.2),
                         r_to = NA_real_, transition_days = NA_integer_)
  step <- data.frame(shape = "step", pairs, transition_days = 1L)
  gradual <- do.call(rbind, lapply(c(7L, 14L, 20L), function(d) {
    data.frame(shape = "gradual", pairs, transition_days = d)
  }))
  tab <- rbind(constant, step, gradual)

  fmt <- function(x) formatC(x, format = "fg")
  tab$major_scenario <- c(constant = "Scenario 1: Constant Rt",
                          step     = "Scenario 2: Stepwise Rt",
                          gradual  = "Scenario 3: Gradual Rt")[tab$shape]
  tab$true_rt <- ifelse(tab$shape == "constant", fmt(tab$r_from),
                        paste0(fmt(tab$r_from), " -> ", fmt(tab$r_to)))
  tab$transition_duration <- ifelse(tab$shape == "constant", "--",
                                    as.character(tab$transition_days))
  prefix <- ifelse(tab$shape == "step", "step", paste0("gradual", tab$transition_days))
  tab$subscenario <- ifelse(tab$shape == "constant",
                            sprintf("constant_%.1f", tab$r_from),
                            paste0(prefix, "_", fmt(tab$r_from), "_to_", fmt(tab$r_to)))

  tab <- tab[, c("subscenario", "major_scenario", "shape", "r_from", "r_to",
                 "transition_days", "true_rt", "transition_duration")]
  rownames(tab) <- NULL
  tab
})

# One registry entry per table row, in the same list(label, build) form as
# rt_scenarios so build_rt_scenario() and scenario_truth_table() accept it, with
# the table's descriptor columns carried along for the output files. The change
# lands on day round(time / 2) - the first day at the new level - which is the
# archive's rt_change = round(time / 2) - 1 (its LAST day at the old level).
build_epiestim_scenarios <- function(table = epiestim_scenario_table, time) {
  change_day <- as.integer(round(time / 2))
  out <- lapply(seq_len(nrow(table)), function(i) {
    row <- table[i, ]
    build <- switch(
      row$shape,
      constant = rt_constant(row$r_from),
      step     = rt_step(row$r_from, row$r_to, change_day),
      gradual  = rt_gradual(row$r_from, row$r_to, change_day, row$transition_days),
      stop("Unknown scenario shape: ", row$shape)
    )
    label <- switch(
      row$shape,
      constant = paste0("Constant Rt = ", row$true_rt),
      step     = paste0("Step ", row$true_rt, " at day ", change_day),
      gradual  = paste0("Gradual ", row$true_rt, " over ", row$transition_days,
                        " days from day ", change_day)
    )
    list(label = label, build = build, shape = row$shape,
         major_scenario = row$major_scenario, true_rt = row$true_rt,
         transition_duration = row$transition_duration)
  })
  names(out) <- table$subscenario
  out
}

# ------------------------------------------------------------------------------
# Construction and validation
# ------------------------------------------------------------------------------

# Evaluate one scenario's `build` closure and check the resulting trajectory.
# A non-positive Rt would make the renewal model meaningless, and a non-finite
# one would silently poison every downstream metric, so both are fatal here
# rather than at scoring time.
build_rt_scenario <- function(scenario, time) {
  if (!is.list(scenario) || !is.function(scenario$build)) {
    stop("A scenario must be a list with a `build` function.")
  }

  R_true <- scenario$build(time)

  if (length(R_true) != time) {
    stop("Scenario '", scenario$label, "' produced ", length(R_true),
         " values for ", time, " days.")
  }
  if (any(!is.finite(R_true)) || any(R_true <= 0)) {
    stop("Scenario '", scenario$label,
         "' produced non-finite or non-positive Rt values.")
  }

  R_true
}

# Build every scenario once, up front, so a malformed scenario fails in seconds
# instead of after the first expensive round of fits.
validate_rt_scenarios <- function(scenarios, time) {
  if (!is.list(scenarios) || length(scenarios) == 0L) {
    stop("rt_scenarios must be a non-empty list.")
  }
  if (is.null(names(scenarios)) || any(!nzchar(names(scenarios)))) {
    stop("Every Rt scenario must be named; the name becomes the scenario column.")
  }

  invisible(lapply(scenarios, build_rt_scenario, time = time))
}

# Check the SI registry and that every requested EpiLPS subset actually exists.
validate_si_scenarios <- function(scenarios, epilps_subset) {
  if (!is.list(scenarios) || length(scenarios) == 0L) {
    stop("si_scenarios must be a non-empty list.")
  }
  if (is.null(names(scenarios)) || any(!nzchar(names(scenarios)))) {
    stop("Every SI scenario must be named.")
  }

  for (name in names(scenarios)) {
    values <- scenarios[[name]]
    if (!all(c("mean", "sd") %in% names(values))) {
      stop("SI scenario '", name, "' must have named 'mean' and 'sd' entries.")
    }
    if (any(!is.finite(values)) || values[["mean"]] <= 1 || values[["sd"]] <= 0) {
      stop("SI scenario '", name, "' has an invalid mean or sd.")
    }
  }

  unknown <- setdiff(epilps_subset, names(scenarios))
  if (length(unknown) > 0L) {
    stop("sim_epilps_si_scenarios names unknown SI scenario(s): ",
         paste(unknown, collapse = ", "))
  }

  invisible(TRUE)
}

# Long-format truth table for reporting and plotting: the instantaneous truth
# plus the trailing average over every estimator window, which is what makes the
# window-averaging offset visible before any fitting happens.
scenario_truth_table <- function(scenarios, time, windows) {
  out <- map_dfr(names(scenarios), function(scenario_id) {
    R_true <- build_rt_scenario(scenarios[[scenario_id]], time)

    map_dfr(sort(unique(as.integer(windows))), function(w) {
      data.frame(
        scenario = scenario_id,
        label = scenarios[[scenario_id]]$label,
        index = seq_len(time),
        estimand_window = w,
        true_R_instant = R_true,
        true_R_window = trailing_mean(R_true, w),
        stringsAsFactors = FALSE
      )
    })
  })

  # Order the label factor by registry position so facets read S1..Sn rather
  # than alphabetically by label text.
  out$label <- factor(
    out$label,
    levels = vapply(scenarios, function(s) s$label, character(1))
  )
  out
}
