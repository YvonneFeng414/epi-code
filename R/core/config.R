# ==============================================================================
# config.R
# User settings for the Queens, NY Rt analysis (EpiEstim vs EpiLPS).
# Sourced by main.R, sim_main.R and make_queens_incidence.R, and by the
# SI-misspecification studies for the shared SI grids. Values are unchanged from
# Section 0 of the original monolithic script.
#
# This file is settings ONLY. It creates no directories, sets no seed and
# probes no hardware; each driver does those things itself, where they can be
# seen. That is what makes it safe for a study with its own settings to source
# this one for the few values it shares.
# ==============================================================================

# ----------------------------
# File and location settings
# ----------------------------
data_file  <- project_paths()$data_file
state_abbr <- "NY"
county_name <- "Queens"

# ----------------------------
# Serial interval settings
# ----------------------------
mean_si <- 7.5
sd_si   <- 3.4
# 60 instead of 30: at sd = 20 a lag-30 truncation would drop ~6.5% of the SI
# mass (shape = (mean/sd)^2 = 0.14 makes the Gamma very heavy-tailed), which
# would distort that scenario after make_si()'s renormalization. 60 keeps
# truncation to ~2.7% for sd = 20 while changing nothing for sd_si = 3.4
# (already ~100% captured within lag 30).
max_si_lag <- 60L

# Sensitivity analysis: both grids are the true value scaled by the same
# set of relative errors (-80% to +80%), so the two analyses are directly
# comparable step-for-step.
relative_error <- c(-0.8, -0.6, -0.4, -0.2, 0,
                      0.2,  0.4,  0.6,  0.8)

assumed_mean_grid <- round(mean_si * (1 + relative_error), 2)
assumed_sd_grid   <- round(sd_si * (1 + relative_error), 2)

# ----------------------------
# SI-misspecification studies (the four *_si_misspec*.R drivers)
# ----------------------------
# The simulation design of Epiestim_simu/seed/run_sensitivity_queens_nohole.R,
# shared here so the four drivers - and the bit-identical series the EpiLPS MAP
# and MALA twins depend on - cannot drift apart. Each driver still keeps its
# own plug-in truth and observation model.
#
# The grid is the one above minus both +/-80% points: at -80% the assumed sd is
# 0.68 days and the assumed mean 1.5, well outside anything a real analysis
# would assume. relative_error and the two grids line up position by position,
# so one logical index subsets all three. main.R's sweeps keep the full grid.
misspec_keep_error        <- !relative_error %in% c(-0.8, 0.8)
misspec_relative_error    <- relative_error[misspec_keep_error]
misspec_assumed_mean_grid <- assumed_mean_grid[misspec_keep_error]
misspec_assumed_sd_grid   <- assumed_sd_grid[misspec_keep_error]

# Seed cases on day 1 only; every later day is simulated and scored.
misspec_I0 <- 2000L

# The observed series has one day (2020-07-08) with zero reported cases because
# the cumulative count was not updated. EpiEstim's 1-day fit turns that into a
# true Rt of 0.0077 - the prior-driven posterior when no cases are seen, not a
# real collapse in transmission. Days whose EpiEstim truth falls below this are
# dropped from every driver's truth path, by calendar day.
misspec_hole_threshold <- 0.2

# A replicate with a run of more than this many zero-case days went extinct and
# is dropped (7 is about one mean serial interval).
misspec_zero_run_threshold <- 7L

misspec_n_sim    <- 100L
misspec_sim_seed <- 123L

# ----------------------------
# Method settings
# ----------------------------
epiestim_window <- 7L   # primary EpiEstim window
K_epilps <- 30L # number of B-spline basis functions for EpiLPS

# The same fixed burn-in/evaluation start is used for all methods.
# It is not selected from the estimated Rt values.
burn_in_days <- 45L

# ----------------------------
# Rt coverage study settings (Section 6)
# ----------------------------
run_rt_coverage <- TRUE

# Number of semi-synthetic replicate epidemics. The Monte Carlo standard error
# of a coverage rate near 0.95 is roughly sqrt(0.95 * 0.05 / n_coverage_reps),
# i.e. about 2.2 percentage points at 100 replicates.
n_coverage_reps <- 100L

# Which full-series fit supplies the plug-in "true" Rt trajectory:
#   "epilps"   - EpiLPS posterior median (smooth, defined on every day)
#   "epiestim" - EpiEstim posterior median from the 7-day sliding window
#   "average"  - pointwise mean of the two medians where both exist
# This is a genuine design choice, not a neutral one: a truth drawn from one
# method's fit is smooth in exactly the way that method assumes, which can
# flatter it. Re-running with a different source is the natural sensitivity
# check on the coverage numbers.
coverage_truth_source <- "epilps"

# EpiEstim arms in the coverage study. A w-day sliding window estimates the
# AVERAGE Rt over [t - w + 1, t], so window = 1 is the only EpiEstim arm whose
# estimand is the instantaneous Rt that EpiLPS targets; it is included so the
# instantaneous comparison is like-for-like, alongside the primary window.
coverage_epiestim_windows <- c(1L, epiestim_window)

# The first max_si_lag days of each replicate are copied from the observed
# series so the renewal recursion starts from a realistic history; they are
# never scored (the evaluation start already excludes them).
coverage_seed <- 20260727L

# ----------------------------
# Forecast settings
# ----------------------------
run_backtest <- TRUE
forecast_stride <- 1L  # use 1 for every day; use 7 for weekly forecast origins
n_posterior_draws <- 5000L

# main.R calls set.seed(forecast_seed) immediately after sourcing this file,
# which is where the original config.R called set.seed() itself. The backtest
# (Section 7) draws from this stream after the coverage study has consumed
# whatever it consumes, so the two sections must stay in that order for the
# forecast numbers to reproduce.
forecast_seed <- 20260726L

# Limit parallel EpiLPS jobs to avoid excessive memory use. Drivers turn this
# into a worker count with default_n_cores(n_cores_cap).
n_cores_cap <- 4L

# ----------------------------
# Output settings
# ----------------------------
output_dir <- project_paths()$results_dir("queens_rt")

# ==============================================================================
# Fully synthetic simulation study (sim_main.R)
# ==============================================================================
# A different kind of study from the coverage study above. There the "true" Rt is
# a plug-in estimate from the real Queens fit and each replicate is seeded with
# observed incidence. Here the true Rt is written down analytically (see
# scenarios.R), the only seed is sim_I0 on day 1, and no real data is read at
# all - so the truth is known exactly rather than assumed.

run_pure_simulation <- TRUE

sim_time  <- 100L   # days per simulated series
sim_n_sim <- 500L   # replicate epidemics per Rt scenario
sim_I0    <- 100L   # seed cases on day 1

# The TRUE serial interval, used only to GENERATE epidemics. The SI handed to the
# estimators comes from si_scenarios in scenarios.R, which is what makes the
# misspecification dimension possible.
sim_mean_si <- 7.5
sim_sd_si   <- 3.4

# 30 is enough: a Gamma with mean 7.5 / sd 3.4 has essentially no mass beyond
# lag 30. epi.R used discr_si(seq(0, time), ...), i.e. lag 100, which is
# numerically equivalent but produces a si_distr of length time + 1 - exactly at
# the boundary EpiEstim allows, for no benefit.
sim_max_si_lag <- 30L

# EpiEstim arms, one per sliding-window length. These are epi.R's nd_values.
# A w-day window estimates the AVERAGE Rt over [t - w + 1, t], so each arm is
# scored against its own trailing-mean target (see sim_study.R).
sim_windows <- c(1L, 4L, 7L, 14L)

# The fixed cross-arm comparison target, held the same for every arm so the
# window-average column is comparable across the table.
sim_reference_window <- 7L

sim_seed <- 123L

# EpiLPS dominates the runtime: its fits depend on the assumed SI, so they cannot
# be shared across SI scenarios the way they are shared across windows, giving
# length(rt_scenarios) * length(si_scenarios) * sim_n_sim fresh-process fits.
# Narrow this to e.g. "correct" to trade misspecification coverage for time.
# Written as a literal so config.R stays independent of scenarios.R.
sim_epilps_si_scenarios <- c(
  "correct", "mean_too_low", "mean_too_high", "sd_too_low", "sd_too_high"
)

sim_output_dir <- project_paths()$results_dir("simulation_study")

# The replicate-level estimate table is one row per
# (Rt scenario, SI scenario, arm, replicate, day). At the full settings above
# that is several million rows / hundreds of MB, so it is off by default; the
# summary and day-level files carry the results. Turn it on for a small n_sim
# when you need to inspect individual fits.
sim_write_replicate_estimates <- FALSE
