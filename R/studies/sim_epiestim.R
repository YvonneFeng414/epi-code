# ==============================================================================
# sim_epiestim.R
# Harness for the pure-synthetic EpiEstim studies in simulations/. Those
# studies generate epidemics from an analytic Rt (simulation/scenarios.R,
# epiestim_scenario_table) seeded with a single I0 and no real data, then run
# EpiEstim under the correct SI or a grid of wrong ones. The heavy lifting is
# the same as the Queens misspecification studies - simulate_replicates(),
# fit_si_grid(), score_interval_fits() - and this file supplies what those
# expect but load_misspec_inputs() would otherwise have read from data: a
# target record with the truth, the seed, the scored days and each window's
# own estimand. It also carries the pieces the archive's trajectory and
# bias-convergence figures need (per-day summaries across replicates, the
# representative replicate, the convergence-day rule) and the appendix-style
# table formatter.
#
# Requires simulation/scenarios.R, simulation/truth.R, scoring/metrics.R and
# studies/si_misspec.R to be sourced first.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Targets
# ------------------------------------------------------------------------------

# One simulation target: a scenario built at horizon `time`, seeded with I0 on
# day 1, scored from the day after the widest window closes so every arm has an
# estimate on every scored day. `truth_by_arm` holds each window's own estimand
# (trailing_mean(); w = 1 is the instantaneous path), full length so callers
# subset by score_days. The descriptor fields come straight from the scenario
# registry entry (build_epiestim_scenarios()).
sim_target <- function(name, scenarios, time, I0, windows) {
  sc <- scenarios[[name]]
  if (is.null(sc)) {
    stop("Unknown scenario '", name, "'. Known: ",
         paste(names(scenarios), collapse = ", "))
  }
  windows <- sort(unique(as.integer(windows)))
  if (any(windows < 1L)) stop("windows must be positive integers.")

  R_true <- build_rt_scenario(sc, time)
  score_from <- max(windows) + 1L
  if (score_from >= time) stop("The widest window leaves no days to score.")

  truth_by_arm <- lapply(windows, function(w) trailing_mean(R_true, w))
  names(truth_by_arm) <- paste0("w", windows)

  list(
    subscenario = name,
    major_scenario = sc$major_scenario,
    true_rt = sc$true_rt,
    transition_duration = sc$transition_duration,
    label = sc$label,
    T = as.integer(time),
    I0 = as.integer(I0),
    windows = windows,
    arms = names(truth_by_arm),
    R_true = R_true,
    seed_incidence = as.integer(I0),
    score_from = score_from,
    score_days = seq.int(score_from, time),
    truth_by_arm = truth_by_arm
  )
}

# Poisson replicates from the target, seeded with its I0. A thin name over
# simulate_replicates() so the drivers read as one vocabulary.
simulate_sim_target <- function(target, si_true, n_sim, seed) {
  simulate_replicates(target$R_true, si_true, target$seed_incidence, n_sim,
                      obs_model = "poisson", seed = seed)
}

# The descriptor columns every output row of these studies starts with.
sim_descriptor <- function(target, arm = NULL, mean_si = NULL, sd_si = NULL) {
  d <- data.frame(
    major_scenario = target$major_scenario,
    subscenario = target$subscenario,
    true_rt = target$true_rt,
    T = target$T,
    I0 = target$I0,
    transition_duration = target$transition_duration,
    stringsAsFactors = FALSE
  )
  if (!is.null(arm)) {
    d$arm <- arm
    d$window <- as.integer(sub("^w", "", arm))
  }
  if (!is.null(mean_si)) d$true_mean_si <- mean_si
  if (!is.null(sd_si)) d$true_sd_si <- sd_si
  d
}

# Repeat a one-row descriptor down the rows of `df` and bind the two.
bind_descriptor <- function(desc, df) {
  if (is.null(df) || nrow(df) == 0L) return(NULL)
  out <- cbind(desc[rep(1L, nrow(df)), , drop = FALSE], df)
  rownames(out) <- NULL
  out
}

# ------------------------------------------------------------------------------
# 2. Scoring
# ------------------------------------------------------------------------------

# Pull one arm out of every multi-arm fit, keeping a failed fit as-is so the
# scorer can drop it.
arm_fits <- function(fits, arm) {
  lapply(fits, function(f) {
    if (!is.null(f$error) && !is.na(f$error)) f else f[[arm]]
  })
}

# Score one setting's fits (one multi-arm fit per replicate) with
# score_interval_fits(), each arm against its own estimand from the target.
# Returns the metrics and by-day frames with `arm` and `window` columns, the
# same shape score_si_grid_per_arm() produces for a whole grid.
score_fits_per_arm <- function(fits, target, ...) {
  rows <- lapply(target$arms, function(arm) {
    sc <- score_interval_fits(arm_fits(fits, arm),
                              target$truth_by_arm[[arm]][target$score_days],
                              target$score_days, ...)
    if (is.null(sc)) return(NULL)
    id <- data.frame(arm = arm, window = as.integer(sub("^w", "", arm)),
                     stringsAsFactors = FALSE)
    list(metrics = cbind(id, sc$metrics), by_day = cbind(id, sc$by_day))
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0L) stop("Every fit failed on every arm.")
  list(
    metrics = do.call(rbind, lapply(rows, function(r) r$metrics)),
    by_day  = do.call(rbind, lapply(rows, function(r) r$by_day))
  )
}

# ------------------------------------------------------------------------------
# 3. Trajectories
# ------------------------------------------------------------------------------

# The per-day view across replicates that the trajectory figures draw: the
# mean estimate, the mean credible interval (the average band an analyst would
# see), the 2.5 / 97.5 percentiles of the point estimate across replicates (how
# far the estimate itself scatters), per-day coverage and the number of
# replicates with an estimate that day. `fits` are one arm's fits, one per
# replicate; `truth` is the FULL-LENGTH estimand for that arm. The archive's
# summarise_rt_trajectories() on the harness's fit lists.
summarise_trajectories <- function(fits, truth, days) {
  ok <- vapply(fits, function(f) is.null(f$error) || is.na(f$error), logical(1))
  fits <- fits[ok]
  if (length(fits) == 0L) stop("Every fit failed.")

  n_days <- length(days)
  aligned <- function(f, field) f[[field]][match(days, f$time)]
  as_matrix <- function(field) {
    matrix(vapply(fits, aligned, numeric(n_days), field = field), nrow = n_days)
  }
  est <- as_matrix("R")
  lo <- as_matrix("lower")
  hi <- as_matrix("upper")
  truth_d <- truth[days]

  # A day with no estimate in any replicate gives NaN from the row means and
  # -Inf/Inf from the quantiles; both are reported as NA.
  as_na <- function(x) { x[!is.finite(x)] <- NA_real_; x }
  q <- function(p) apply(est, 1L, stats::quantile, probs = p, na.rm = TRUE, names = FALSE)

  data.frame(
    day = as.integer(days),
    true_R = truth_d,
    mean_est = as_na(rowMeans(est, na.rm = TRUE)),
    mean_lower = as_na(rowMeans(lo, na.rm = TRUE)),
    mean_upper = as_na(rowMeans(hi, na.rm = TRUE)),
    mc_lower = as_na(suppressWarnings(q(0.025))),
    mc_upper = as_na(suppressWarnings(q(0.975))),
    coverage = as_na(rowMeans(lo <= truth_d & truth_d <= hi, na.rm = TRUE)),
    n_used = rowSums(!is.na(est))
  )
}

# One replicate's own trajectory, for the single-epidemic panel.
fit_trajectory <- function(fit, truth, days) {
  aligned <- function(field) fit[[field]][match(days, fit$time)]
  data.frame(
    day = as.integer(days),
    true_R = truth[days],
    est = aligned("R"),
    lower = aligned("lower"),
    upper = aligned("upper")
  )
}

# The replicate whose final size is closest to the median, so the single-
# epidemic panel shows a typical epidemic rather than whichever one is first.
pick_representative <- function(series) {
  final_size <- vapply(series, sum, numeric(1))
  which.min(abs(final_size - stats::median(final_size)))
}

# ------------------------------------------------------------------------------
# 4. Convergence to the Euler-Lotka asymptote
# ------------------------------------------------------------------------------

# First day from which the per-day bias stays within tol (relative) of the
# Euler-Lotka asymptote for the rest of the series. Taking the day after the
# LAST violation rather than the first success stops a curve that merely
# crosses the band on its way past from being credited. The correctly
# specified SI has an asymptote of exactly zero, where a relative tolerance is
# meaningless, so it reports NA. Ported verbatim from the archive.
convergence_day <- function(day, bias, el, tol = 0.1) {
  if (abs(el) < 1e-10) return(NA_integer_)
  ok <- !is.na(bias) & abs(bias - el) < tol * abs(el)
  usable <- !is.na(bias)
  if (!any(ok)) return(NA_integer_)
  last_bad <- suppressWarnings(max(day[usable & !ok]))
  if (!is.finite(last_bad)) return(as.integer(min(day[usable])))
  after <- day[day > last_bad]
  if (length(after) == 0L) return(NA_integer_)
  as.integer(min(after))
}

# ------------------------------------------------------------------------------
# 5. The appendix-style table
# ------------------------------------------------------------------------------

scenario_section_order <- c("Scenario 1: Constant Rt", "Scenario 2: Stepwise Rt",
                            "Scenario 3: Gradual Rt")
scenario_section_labels <- c("Constant Rt", "Stepwise Rt", "Gradual Rt")

# The archive's display order within a section. Note "0.8 -> 1.2" sorts before
# "1.2 -> 0.8", the reverse of the registry order - kept as the appendix had it.
scenario_rt_order <- c("1", "0.8", "1.2", "0.8 -> 1.2", "1.2 -> 0.8",
                       "1.3 -> 1.1", "1.1 -> 1.3", "0.9 -> 0.7", "0.7 -> 0.9")

window_label <- function(w) {
  ifelse(w == 1L, "Daily", ifelse(w == 7L, "Weekly", paste0(w, "-day")))
}

# The scenario study's metrics laid out as the archive's appendix table: one
# header row per section, rows sorted by the archive's rule, repeated
# True Rt / T / I0 / Transition duration blanked hierarchically, Bias and MSE
# in 3-significant-digit scientific notation, coverage and width to 3 dp.
# Every column is character; the CSV is the deliverable.
format_scenario_table <- function(metrics) {
  m <- metrics
  m$section <- match(m$major_scenario, scenario_section_order)
  if (anyNA(m$section)) stop("Unknown major_scenario in metrics.")
  m$rt_order <- match(m$true_rt, scenario_rt_order)
  m$rt_order[is.na(m$rt_order)] <- 99L
  m$transition_order <- suppressWarnings(as.numeric(m$transition_duration))
  m$transition_order[is.na(m$transition_order)] <- 0
  m <- m[order(m$section, m$rt_order, m$T, m$I0, m$transition_order, m$window), ]

  blank_row <- function(cols) {
    data.frame(as.list(stats::setNames(rep("", length(cols)), cols)),
               check.names = FALSE, stringsAsFactors = FALSE)
  }
  cols <- c("True Rt", "T", "I0", "Transition duration", "Window",
            "Bias", "MSE", "Coverage rate", "CI width")

  sections <- lapply(sort(unique(m$section)), function(s) {
    d <- m[m$section == s, ]
    n <- nrow(d)
    keys <- list(d$true_rt, d$T, d$I0, d$transition_duration)
    # TRUE where keys 1..k all equal the previous row's: blank that level.
    repeated_to <- function(k) {
      same <- c(FALSE, rep(TRUE, n - 1L))
      for (j in seq_len(k)) same <- same & c(FALSE, keys[[j]][-1L] == keys[[j]][-n])
      same
    }
    display <- function(x, k) { x <- as.character(x); x[repeated_to(k)] <- ""; x }

    header <- blank_row(cols)
    header[["True Rt"]] <- scenario_section_labels[s]
    body <- data.frame(
      `True Rt` = display(d$true_rt, 1L),
      T = display(d$T, 2L),
      I0 = display(d$I0, 3L),
      `Transition duration` = display(d$transition_duration, 4L),
      Window = window_label(d$window),
      Bias = sprintf("%.3E", d$Bias),
      MSE = sprintf("%.3E", d$MSE),
      `Coverage rate` = formatC(d$Coverage95, format = "f", digits = 3),
      `CI width` = formatC(d$MeanCIWidth, format = "f", digits = 3),
      check.names = FALSE, stringsAsFactors = FALSE
    )
    rbind(header, body)
  })
  out <- do.call(rbind, sections)
  rownames(out) <- NULL
  out
}

# ------------------------------------------------------------------------------
# 6. Console helpers
# ------------------------------------------------------------------------------

# After the first unit of work, say how long the rest will take, so a run that
# will take hours announces itself before it has taken them.
report_cost_projection <- function(elapsed_secs, done, total, unit = "cell") {
  remaining <- (total - done) * elapsed_secs / done
  cat(sprintf("\n%d of %d %ss done in %.0f s; projected %.1f min for the remaining %d\n",
              done, total, unit, elapsed_secs, remaining / 60, total - done))
  invisible(remaining)
}
