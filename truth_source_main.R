################################################################################
# this file is used for testing the self-consistency of the estimators
# only implemented on eplips and epiestim
# 
# Self-consistency study on the observed Queens series.
#
# The third kind of simulation in this project, and the only one where each
# estimator is paired with a truth it produced itself:
#
#   Observed Queens incidence
#     |-- EpiEstim 1-day full fit -- plug-in truth -- N replicates -- refit
#     |                                                               EpiEstim only
#     `-- EpiLPS      full fit ---- plug-in truth -- N replicates -- refit
#                                                                    EpiLPS only
#
# There is NO cross-fitting: EpiEstim is never scored against the EpiLPS truth
# or the reverse. The two arms are independent.
#
# WHAT THIS DOES AND DOES NOT MEASURE. Each arm is a parametric bootstrap
# calibration check: regenerate data from a method's own fit, and see whether
# that method recovers it with correctly calibrated uncertainty. That is a well
# posed question, and a method can fail it.
#
# It is NOT a head-to-head verdict on which estimator is better. Every arm is
# graded on a truth that is smooth in exactly the way that arm assumes, which is
# maximally favourable to it, so the two arms' numbers are not comparable as a
# ranking. A method can be perfectly self-consistent and still be wrong about the
# real epidemic. Compare CI WIDTHS between the arms rather than coverage: equal
# coverage at very different widths means very different sharpness.
#
# Relation to the other two studies. main.R Section 6 is semi-synthetic but picks
# ONE plug-in truth via coverage_truth_source and is tangled into the full Queens
# analysis. sim_main.R is fully synthetic, with an analytic truth that depends on
# no fit. This script sits between them: real data supplies the truth, but each
# method supplies its own.
#
# Run from r-proj/:  Rscript truth_source_main.R
# Output: results/truth_source_<si_discretization>/
################################################################################

rm(list = ls())

# ==============================================================================
# 0. Packages and modules
# ==============================================================================
# config.R is deliberately NOT sourced: its settings belong to the main.R run
# (a 7-day window, max_si_lag 60) and this study carries its own.

source(file.path("R", "core", "init.R"))
require_packages(c("EpiEstim", "EpiLPS", "dplyr", "tidyr", "ggplot2", "purrr"))
source_project()

# ==============================================================================
# 1. Settings
# ==============================================================================

paths <- project_paths()
incidence_file <- paths$incidence_csv

# The SI is treated as known and correct throughout: it generates the epidemics
# AND is handed to both estimators. SI misspecification is sim_main.R's question,
# not this one.
true_mean_si <- 7.5
true_sd_si   <- 3.4
max_si_lag   <- 30L

# Which discretization turns (mean, sd) into a probability vector. The two are
# NOT the same distribution:
#
#   "discr_si" - EpiEstim::discr_si() via make_si(). A gamma fitted to mean - 1,
#                shifted one day, discretized with the Cori et al. triangular
#                kernel. Support 0..max_si_lag.
#   "idist"    - EpiLPS::Idist(). A gamma fitted to the mean directly, binned at
#                +/- 0.5 and renormalized. Support 1..Dmax, Dmax chosen by the
#                0.9999 quantile (27 here), so max_si_lag does not apply.
#
# Both realise mean/sd ~ 7.5/3.4 but differ by total variation 0.019, and by 3.2x
# in the lag-1 cell - the weight that drives the most recent days' contribution
# to infectiousness. EpiLPS's own documentation uses Idist, so anyone following
# that workflow is on the second one; running both is how we find out whether the
# headline result depends on the convention.
#
# Whichever is chosen, BOTH estimators receive it. That invariant is what keeps
# the two arms comparable and must not be broken.
si_discretization <- "discr_si"

# The only EpiEstim arm. window = 1 both matches the fit that generates the
# EpiEstim truth and makes the arm's estimand the instantaneous Rt, which is what
# EpiLPS targets - so the two arms are answering the same question.
epiestim_window  <- 1L
reference_window <- 1L

n_sim <- 100L

# Seeding. A single I0 on day 1 is epi.R's and sim_main.R's convention, and it is
# fine there because their truth is a synthetic path with no scale attached to
# it. It does NOT transfer to a truth taken from a real epidemic: one seed day
# supplies far too little infectiousness next to the ~30 days of history the
# observed series had behind it, so the replicates collapse to single digits and
# stay there. Measured on this series, I0 = 100 gives a median of 1-12 cases/day
# against an observed 86-786, and every metric then describes a near-extinct
# epidemic rather than Queens - the 1-day EpiEstim arm's mean CI width blows up
# past 3.0 purely from Poisson noise on tiny counts.
#
# Seeding with the first seed_burn_in observed days instead, which is what
# main.R Section 6 does, keeps the replicates at the observed scale. Those copied
# days are never scored: seed_burn_in already excludes them.
seed_from_observed <- TRUE
I0 <- 100L   # only used when seed_from_observed is FALSE

# Drop the leading days where EpiEstim's estimate explodes on the epidemic's tiny
# early counts. Data-dependent, which the rest of this repo avoids; see the note
# in section 5.
truth_max_R <- 5

# Extra unscored days at the head of each replicate: seeding is an artificial
# start the renewal recursion needs roughly a full SI span to forget.
#
# Deliberately a literal, NOT max_si_lag. Idist truncates its support at 27 while
# discr_si runs to 30, so tying the burn-in to the SI length would make the two
# discretizations seed differently and score different days - which would confound
# exactly the comparison this setting exists to enable.
seed_burn_in <- 30L

K_epilps <- 30L
sim_seed <- 20260812L

# Figure colours come from R/core/palette.R (arm_levels, arm_colors,
# estimator_colors), shared with compare_si_discretization.R so a reader moving
# between the two carries the same mapping.

output_dir <- paths$truth_source_dir(si_discretization)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

n_cores <- default_n_cores(cap = 4L)

# ==============================================================================
# 2. Load the observed incidence series
# ==============================================================================
# This file is already the finished dates/I product that main.R writes out
# (see read_incidence_csv() for what is checked).

incidence_data <- read_incidence_csv(incidence_file)
n_obs <- nrow(incidence_data)

# ==============================================================================
# 3. Shared serial interval
# ==============================================================================
# One SI object for everything, so the two estimators cannot differ merely
# because they built the distribution differently. This is why EpiEstim runs
# under non_parametric_si with an explicit si_distr rather than parametric_si:
# parametric_si would have EpiEstim construct its own SI internally while EpiLPS
# received a separately constructed one.

# Both branches return the EpiEstim convention: index 1 is lag 0 and is zero.
# Idist() emits lags 1..Dmax already summing to 1, so prepending the zero cell is
# all that is needed to line the two up.
si_full <- switch(
  si_discretization,
  discr_si = make_si(true_mean_si, true_sd_si, max_si_lag),
  idist = c(0, EpiLPS::Idist(mean = true_mean_si, sd = true_sd_si,
                             dist = "gamma")$pvec),
  stop("si_discretization must be 'discr_si' or 'idist'.")
)

if (any(!is.finite(si_full)) || any(si_full < 0) ||
    abs(sum(si_full) - 1) > 1e-8 || si_full[1L] != 0) {
  stop("The constructed SI is not a valid distribution with a zero lag-0 cell.")
}

si_epilps <- si_full[-1L]

si_lags <- seq.int(0L, length(si_full) - 1L)
realised_si <- si_moments(si_full)

cat("\n===== Self-consistency study: each method against its own truth =====\n")
cat("Series:            ", as.character(min(incidence_data$dates)), " to ",
    as.character(max(incidence_data$dates)), " (", n_obs, " days)\n", sep = "")
cat("SI discretization: ", si_discretization, ", support lag 0-",
    max(si_lags), "\n", sep = "")
cat("SI realised mean/SD: ", sprintf("%.4f / %.4f", realised_si[["mean"]],
    realised_si[["sd"]]), " (requested ", true_mean_si, " / ", true_sd_si, ")\n",
    sep = "")
cat("SI mean/SD:        ", true_mean_si, " / ", true_sd_si,
    " (known and correct everywhere)\n", sep = "")
cat("Replicates per arm:", n_sim, "\n")

# ==============================================================================
# 4. Two plug-in truths, one per method
# ==============================================================================

cat("\n===== Fitting the two plug-in truths on the observed series =====\n")

epiestim_truth_fit <- fit_epiestim_full(
  incidence_df = incidence_data,
  si_full = si_full,
  window_length = epiestim_window
)

epilps_truth_fit <- bind_rows(
  run_epilps_jobs(
    list(list(
      job_id = "observed",
      incidence = incidence_data$I,
      si_epilps = si_epilps,
      K = K_epilps,
      dates = as.character(incidence_data$dates),
      return_type = "full"
    )),
    cores = 1L
  )
)

if (all(is.na(epilps_truth_fit$q50))) {
  stop("EpiLPS failed on the observed series: ", epilps_truth_fit$converged[1L])
}

# ==============================================================================
# 5. Common day range
# ==============================================================================
# The truncation rule is applied to the EPIESTIM path only, because it is the one
# that explodes on the epidemic's tiny early counts; EpiLPS is smooth and does
# not. Both truths are then cut to the SAME range so the two arms' bias, MSE and
# coverage are computed over identical days and remain comparable.
#
# Note this makes the evaluation window depend on the estimates themselves, which
# main.R and sim_main.R deliberately avoid (they fix it from settings). Change
# truth_max_R to a fixed burn-in if that dependence matters for your use.

qualifying <- which(epiestim_truth_fit$R <= truth_max_R)
if (length(qualifying) == 0L) {
  stop("No EpiEstim estimate falls at or below truth_max_R = ", truth_max_R,
       "; the series never settles, so no common range can be chosen.")
}

common_start <- epiestim_truth_fit$index[qualifying[1L]]
common_index <- seq.int(common_start, n_obs)
sim_time <- length(common_index)

if (sim_time <= seed_burn_in + 1L) {
  stop("The common range leaves no days to score after the seed burn-in.")
}

# build_true_rt() interpolates onto every day of the common range and extends
# flat at the edges, so both paths are defined on all sim_time days even though
# the EpiEstim fit starts at index 2 and EpiLPS at index 1.
rt_epiestim <- build_true_rt(
  index = epiestim_truth_fit$index - common_start + 1L,
  value = epiestim_truth_fit$R,
  n_days = sim_time
)

rt_epilps <- build_true_rt(
  index = epilps_truth_fit$index - common_start + 1L,
  value = epilps_truth_fit$q50,
  n_days = sim_time
)

sim_dates <- incidence_data$dates[common_index]

cat("Common range:      day ", common_start, " to ", n_obs, " (", sim_time,
    " days, from ", as.character(min(sim_dates)), ")\n", sep = "")
cat("EpiEstim truth Rt: ", sprintf("%.3f to %.3f, mean %.3f",
    min(rt_epiestim), max(rt_epiestim), mean(rt_epiestim)), "\n", sep = "")
cat("EpiLPS truth Rt:   ", sprintf("%.3f to %.3f, mean %.3f",
    min(rt_epilps), max(rt_epilps), mean(rt_epilps)), "\n", sep = "")

write.csv(
  data.frame(
    index = seq_len(sim_time),
    date = sim_dates,
    true_R_epiestim = rt_epiestim,
    true_R_epilps = rt_epilps
  ),
  file.path(output_dir, "plugin_truth_paths.csv"),
  row.names = FALSE
)

# The two truths side by side, before either arm is simulated from.
ggsave(
  filename = file.path(output_dir, "plugin_truth_comparison.png"),
  plot = plot_plugin_truths(
    data.frame(day = seq_len(sim_time), EpiEstim = rt_epiestim, EpiLPS = rt_epilps)
  ),
  width = 10, height = 6, dpi = 300
)

# ==============================================================================
# 6. Simulate replicates from each truth
# ==============================================================================
# Each fitted truth is wrapped as an rt_scenarios-shaped entry so the whole
# sim_main.R machinery applies unchanged: build_rt_scenario() only asks for a
# build(time) closure returning `time` positive finite values. Written out
# explicitly rather than in a loop so each closure captures the intended vector.

truth_scenarios <- list(
  epiestim_truth = list(
    label = "Plug-in truth: EpiEstim (1-day window)",
    build = function(time) rt_epiestim
  ),
  epilps_truth = list(
    label = "Plug-in truth: EpiLPS",
    build = function(time) rt_epilps
  )
)

validate_rt_scenarios(truth_scenarios, sim_time)

# Both arms are seeded identically, so any difference between them comes from the
# truth path and the estimator, never from the starting conditions.
seed_incidence <- if (isTRUE(seed_from_observed)) {
  incidence_data$I[common_index[seq_len(seed_burn_in)]]
} else {
  I0
}

cat("\n===== Simulating replicates =====\n")
cat("Seed: ",
    if (isTRUE(seed_from_observed)) {
      paste0("the first ", seed_burn_in, " observed days of the common range (",
             paste(range(seed_incidence), collapse = " to "), " cases/day)")
    } else {
      paste0(I0, " cases on day 1 only")
    },
    "\n", sep = "")

replicate_sets <- list()
for (i in seq_along(truth_scenarios)) {
  arm_id <- names(truth_scenarios)[i]

  replicate_sets[[arm_id]] <- simulate_scenario_replicates(
    scenario = truth_scenarios[[arm_id]],
    time = sim_time,
    n_sim = n_sim,
    I0 = seed_incidence,
    si_true = si_full,
    seed = sim_seed + (i - 1L)
  )

  set_info <- replicate_sets[[arm_id]]
  cat(arm_id, ": kept ", length(set_info$series), " / ", set_info$n_requested,
      " replicates",
      if (set_info$n_failed > 0L) {
        paste0(" (", set_info$n_failed, " diverged and were discarded)")
      } else {
        ""
      },
      "; median final-day incidence ",
      stats::median(vapply(set_info$series, function(x) x[sim_time], numeric(1))),
      "\n", sep = "")
}

# ==============================================================================
# 7. Refit - each arm with its own method only
# ==============================================================================
# The run_epiestim / run_epilps toggles are what keep the arms matched. Anything
# else here would be cross-fitting, which this design excludes by construction.

cat("\n===== Fitting =====\n")

cat("  epiestim_truth x EpiEstim (1-day window)\n")
epiestim_arm <- fit_scenario_grid(
  replicates = replicate_sets$epiestim_truth,
  si_assumed = si_full,
  dates = sim_dates,
  windows = epiestim_window,
  reference_window = reference_window,
  K_epilps = K_epilps,
  n_cores = n_cores,
  run_epiestim = TRUE,
  run_epilps = FALSE
)

cat("  epilps_truth x EpiLPS\n")
epilps_arm <- fit_scenario_grid(
  replicates = replicate_sets$epilps_truth,
  si_assumed = si_full,
  dates = sim_dates,
  windows = epiestim_window,
  reference_window = reference_window,
  K_epilps = K_epilps,
  n_cores = n_cores,
  run_epiestim = FALSE,
  run_epilps = TRUE
)

# The arm id goes in the `scenario` column and si_scenario is a constant, because
# that is what the scoring and plotting helpers group on. There is only one SI
# here, and it is the correct one.
draws <- bind_rows(
  epiestim_arm %>% mutate(scenario = "epiestim_truth", .before = 1L),
  epilps_arm %>% mutate(scenario = "epilps_truth", .before = 1L)
) %>%
  mutate(si_scenario = "correct", .after = scenario) %>%
  arrange(scenario, method, replicate, index)

# ==============================================================================
# 8. Score
# ==============================================================================
# With every window at 1, sim_score_start_index() returns only 2, so the seed
# burn-in is what actually governs where scoring begins.

window_start <- sim_score_start_index(epiestim_window, reference_window)
start_index <- max(window_start, seed_burn_in + 1L)

cat("\nScoring starts at index ", start_index, " (",
    if (start_index > window_start) "the seed burn-in binds" else "the window binds",
    "; ", sim_time - start_index + 1L, " days scored).\n", sep = "")

draws <- draws %>% filter(index >= start_index)

metrics_summary <- draws %>%
  group_by(scenario, si_scenario, method) %>%
  group_modify(~rt_coverage_metrics(.x)) %>%
  ungroup() %>%
  mutate(MSE_ownEstimand = RMSE_ownEstimand^2, .after = RMSE_ownEstimand) %>%
  arrange(scenario, method)

# Coverage95_Rt_instantaneous and Coverage95_Rt_windowAvg are duplicates of
# Coverage95_Rt_ownEstimand by construction here: reference_window is 1 and every
# arm's estimand_window is 1, so build_truth_lookup()'s three truth columns are
# the same vector. They are kept in the CSV only because rt_coverage_metrics() is
# shared with main.R and sim_main.R.
write.csv(
  metrics_summary,
  file.path(output_dir, "truth_source_metrics.csv"),
  row.names = FALSE
)

cat("\n===== Each method against its own plug-in truth =====\n")
cat("Both arms target the instantaneous Rt, so Estimand is the same in each row.\n")
cat("Coverage is not a ranking across arms - each is graded on a truth it made.\n")
cat("MeanCIWidth is the comparable column.\n\n")
print(
  as.data.frame(
    metrics_summary %>%
      select(scenario, method, N_replicates, N_day_estimates,
             Coverage95_Rt_ownEstimand, MCSE_Coverage95_Rt_ownEstimand,
             MeanCIWidth, Bias_ownEstimand, RMSE_ownEstimand, MSE_ownEstimand)
  ),
  row.names = FALSE,
  digits = 4
)

# ==============================================================================
# 9. Day-level results and figures
# ==============================================================================

# sim_coverage_by_day() pivots the instantaneous and window targets into two rows
# per estimate. Here they hold identical numbers, so keeping both would give the
# figure two identical facet columns; own_estimand selects the instantaneous one.
by_day <- sim_coverage_by_day(draws, reference_window) %>%
  filter(own_estimand)

write.csv(
  by_day,
  file.path(output_dir, "truth_source_coverage_by_day.csv"),
  row.names = FALSE
)

truth_estimates <- sim_mean_estimate_by_day(draws, si_scenario_shown = "correct")

write.csv(
  truth_estimates,
  file.path(output_dir, "truth_source_estimates.csv"),
  row.names = FALSE
)

# With a single 1-day window the trailing-mean layer inside plot_scenario_truth()
# is empty, so these show the black instantaneous truth and the fitted arms only.
# Both figures therefore need their own captions: the defaults describe
# sim_main.R's synthetic scenarios and its windowed estimands, neither of which
# applies here.
truth_table <- scenario_truth_table(truth_scenarios, sim_time, epiestim_window)

ggsave(
  filename = file.path(output_dir, "truth_paths.png"),
  plot = plot_scenario_truth(
    truth_table,
    title = expression(paste("Plug-in true ", R[t], " paths from the observed Queens series")),
    subtitle = paste(
      "Each panel is one method's full-series fit, used as the truth that",
      "generates that arm's replicates.\nThe first days are shown but never",
      "scored: the seed burn-in excludes them."
    )
  ),
  width = 10, height = 7, dpi = 300
)

# Restricted to the scored range. Over the full range the early spike (Rt near 7)
# compresses the y axis and hides the region the metrics actually describe.
ggsave(
  filename = file.path(output_dir, "truth_with_estimates.png"),
  plot = plot_scenario_truth(
    truth_table %>% filter(index >= start_index),
    estimates = truth_estimates,
    title = expression(paste("Each method against its own plug-in true ", R[t])),
    subtitle = paste(
      "Black = that arm's plug-in truth, drawn underneath. Coloured = the fitted",
      "arm, averaged over replicates.\nScored range only. A fitted line lying on",
      "the black one means no bias - at this scale EpiEstim's bias (max 0.03)\nis",
      "narrower than the stroke, so see the residual figure to actually read it."
    ),
    estimate_colors = estimator_colors
  ),
  width = 10, height = 7, dpi = 300
)

# --- Residual figure ----------------------------------------------------------
# The overlay above cannot resolve a bias of 0.007 against an axis spanning the
# whole Rt range. This puts each arm's bias on its own scale, against the width
# of the interval that is supposed to absorb it.
arm_residuals <- inner_join(
  by_day %>% select(scenario, method, index, target_R, mean_ci_width),
  truth_estimates %>% select(scenario, method, method_family, index, mean_R),
  by = c("scenario", "method", "index")
) %>%
  left_join(distinct(truth_table, scenario, label), by = "scenario") %>%
  mutate(bias = mean_R - target_R, half = mean_ci_width / 2)

write.csv(
  arm_residuals,
  file.path(output_dir, "truth_source_residuals.csv"),
  row.names = FALSE
)

cat("\nDays where |bias| exceeds the reported CI half-width:\n")
print(
  as.data.frame(
    arm_residuals %>%
      group_by(scenario, method) %>%
      summarise(
        days = n(),
        exceeding = sum(abs(bias) > half),
        pct = 100 * mean(abs(bias) > half),
        .groups = "drop"
      )
  ),
  row.names = FALSE,
  digits = 3
)

ggsave(
  filename = file.path(output_dir, "truth_residuals.png"),
  plot = plot_arm_residuals(
    arm_residuals,
    method_colors = estimator_colors,
    subtitle = paste(
      "Line = mean estimate minus that arm's own plug-in truth. Band = +/- half",
      "the mean reported 95% CI\nwidth. Where the line leaves the band, the bias",
      "alone exceeds the interval, so coverage on that\nday is impossible however",
      "well the noise is modelled. Note the two panels' y scales differ."
    )
  ),
  width = 10, height = 7, dpi = 300
)

# --- Headline metrics figure --------------------------------------------------
# The bias share is the interpretive payload: a credible interval is sized for
# sampling noise, so an arm whose error is mostly systematic bias cannot cover
# its own truth no matter how well its noise model is calibrated.
day_bias <- arm_bias_share(by_day, truth_estimates)

panel_source <- metrics_summary %>%
  left_join(day_bias, by = c("scenario", "method")) %>%
  transmute(
    method = factor(method, levels = arm_levels),
    coverage = Coverage95_Rt_ownEstimand,
    coverage_mcse = MCSE_Coverage95_Rt_ownEstimand,
    width = MeanCIWidth,
    rmse = RMSE_ownEstimand,
    bias_share = 100 * mean_bias2 / MSE_ownEstimand
  )

metric_panels <- build_arm_metric_panels(panel_source)

ggsave(
  filename = file.path(output_dir, "arm_metrics.png"),
  plot = plot_arm_metrics(
    metric_panels,
    references = arm_metric_references(),
    title = "Each estimator scored against its own plug-in truth",
    subtitle = paste0(
      "Queens incidence, ", n_sim, " replicates per arm, ",
      sim_time - start_index + 1L, " scored days, SI known and correct.\n",
      "The two arms have DIFFERENT truths, so this is not a ranking: the ",
      "1-day EpiEstim arm does no smoothing,\nso recovering its own fit is ",
      "close to automatic. Error bar on coverage is +/- 1.96 MCSE; dashed ",
      "line is nominal 0.95."
    ),
    method_colors = arm_colors
  ),
  width = 10, height = 7.5, dpi = 300
)

ggsave(
  filename = file.path(output_dir, "coverage_by_day.png"),
  plot = plot_sim_coverage_by_day(
    by_day,
    si_scenario_shown = "correct",
    subtitle = paste(
      "Semi-synthetic replicates: the truth is each method's own fit to the",
      "observed Queens series, and the SI\nis known and correct. Each arm is",
      "refit only on its own truth, so a panel below 0.95 is that method",
      "\nfailing to recover itself. Dashed line = nominal 0.95."
    ),
    method_colors = arm_colors
  ),
  width = 9, height = 7, dpi = 300
)

# ==============================================================================
# 10. Save session information
# ==============================================================================

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "sessionInfo.txt")
)

cat("\nSelf-consistency study complete. Results were written to:\n")
cat(normalizePath(output_dir), "\n")
