# ============================================================
# SETTINGS: four-estimator Rt trajectory experiment
#
# One scenario, gradual20_1.1_to_1.3 (Rt 1.1 -> 1.3 over a
# 20-day transition), fitted by EpiEstim, EpiFilter,
# EpiLPS (MAP) and EpiLPS (MALA) on one shared set of
# simulated epidemics.
#
# Smoke mode:  Rscript run_trajectories.R smoke
#   (or SMOKE=1 in the environment) - 4 replicates, written to
#   results_smoke/. Every estimator setting is kept at its full
#   value; only the replicate count drops.
# ============================================================


# ============================================================
# 0. Mode and paths
# ============================================================

cli_args <- commandArgs(trailingOnly = TRUE)

SMOKE <- "smoke" %in% cli_args ||
  identical(Sys.getenv("SMOKE"), "1")


# EpiFilter is not on CRAN; these are Parag's scripts, vendored
# verbatim in queens_rt/.
epifilter_dir <- file.path("..", "queens_rt", "R_epifilter")



# ============================================================
# 1. Data-generating mechanism (matches the seed study)
# ============================================================

sim_seed <- 123L

traj_subscenarios <- c("gradual20_1.1_to_1.3")

traj_time <- 100L
# Seed size. Override without editing:  I0=100 Rscript run_trajectories.R
# A non-default I0 writes to its own folder, e.g. results_I0_100/.
traj_I0 <- as.numeric(Sys.getenv("I0", "2000"))

output_dir <- paste0(
  if (SMOKE) "results_smoke" else "results",
  if (traj_I0 != 2000) paste0("_I0_", traj_I0) else ""
)

traj_n_sim <- if (SMOKE) 4L else 100L

# TRUE serial interval, used both to simulate and to fit
mean_si <- 7.5
sd_si   <- 3.4

start_date <- as.Date("2025-01-01")


# ============================================================
# 2. Estimator settings
# ============================================================

# EpiEstim: 1-day sliding window, as in the seed figures
epiestim_nd <- 1L

# EpiFilter - grid and state model from vignetteCOVID.R, as in
# queens_rt/epifilter_si_misspec.R. filter_start = 8: with a
# 2000-case seed on day 1, total infectiousness is truncated
# on the first days and the implied Rt runs far off the grid,
# which NaNs the recursion (see that script for the numbers).
epifilter_R_min  <- 0.01
epifilter_R_max  <- 10
epifilter_m      <- 1000L
epifilter_eta    <- 0.1
epifilter_alpha  <- 0.025
epifilter_start  <- 8L

# EpiLPS B-spline basis size, as in the queens_rt studies
epilps_K <- 30L

# EpiLPS burn-in: the fit sees days epilps_start..T only, the
# same start as EpiFilter's recursion. EpiLPS smooths the
# log-incidence itself with the spline (estimR: muhat =
# exp(B theta)), so the day-1 seed spike (2000, then ~2 on
# day 2) is a misfit it can only absorb as overdispersion - on
# the full series rho ~ 8, which inflates every interval. Unlike
# EpiFilter, EpiLPS builds its own total infectiousness from
# what it is given, so days 1..epilps_start-1 are also missing
# from its early renewal sums.
epilps_start <- epifilter_start

# EpiLPS's observation model is negative binomial and it cannot
# fit a Poisson one, so its arms get their own data, simulated
# with rnbinom (simulation_nb.R): variance mu + mu^2 / rho.
# rho is an ASSUMED noise level - the simulated scenario has no
# observed series to estimate it from. 20 sits between the two
# Queens estimates (17.36 on the raw series, 22.49 cleaned).
# EpiEstim and EpiFilter keep the Poisson data.
epilps_nb_rho <- 20

# Which data each estimator is fitted to
estimator_data <- c(
  epiestim    = "poisson",
  epifilter   = "poisson",
  epilps_map  = "negbin",
  epilps_mala = "negbin"
)

# EpiLPS MALA - the queens_rt defaults. Not reduced in smoke
# mode, so the smoke test exercises the real sampler settings.
mcmc_niter  <- 5000L
mcmc_burnin <- 2000L
mcmc_seed   <- 20260826L

# Replicate r's MALA seed, independent of fork scheduling
rep_seed <- function(r) mcmc_seed + 1000L * as.integer(r)

n_cores <- max(1L, parallel::detectCores(logical = TRUE) - 1L)

# Metrics are also reported over the days every estimator
# covers, so the four can be compared on the same window.
common_from_day <- epifilter_start

# Figure axis ranges only - never the metrics. For its first
# ~17 fitted days (99% of the SI mass) EpiLPS's renewal sum is
# still filling up and it returns Rt in the hundreds, which
# would flatten every panel. Those days are still drawn.
display_from_day <- 25L


# ------------------------------------------------------------
# SI misspecification (run_si_misspec.R only)
#
# The data are always simulated with the TRUE SI above; each
# arm changes only the SI the estimators are told. One moment
# is moved at a time, the other held at its true value.
# ------------------------------------------------------------

si_misspec_rel <- 0.6

si_misspec_pct <- round(100 * si_misspec_rel)

si_arms <- data.frame(
  arm = c("Correct SI",
          sprintf("Mean -%d%%", si_misspec_pct), sprintf("Mean +%d%%", si_misspec_pct),
          sprintf("SD -%d%%", si_misspec_pct), sprintf("SD +%d%%", si_misspec_pct)),
  a_mean = mean_si * c(1, 1 - si_misspec_rel, 1 + si_misspec_rel, 1, 1),
  a_sd = sd_si * c(1, 1, 1, 1 - si_misspec_rel, 1 + si_misspec_rel),
  stringsAsFactors = FALSE
)

# Common burn-in for this experiment: every estimator and arm is
# scored (and drawn) from the same day. EpiEstim will not
# estimate before a day set by the ASSUMED mean (day 4 at mean
# -50%, 8 at the true mean, 12 at mean +50%), so the latest
# first-estimated day across all estimators and arms is used,
# as the seed study's sens_burn_in_windows does. It depends on
# the grid, so run_si_misspec.R computes it from the fits
# (si_misspec_score_from) rather than it being fixed here. The
# fits are unchanged: EpiFilter and EpiLPS are still fitted
# from day 8.


# ============================================================
# 3. Display order and labels
# ============================================================

estimator_labels <- c(
  epiestim   = "EpiEstim",
  epifilter  = "EpiFilter (smoother)",
  epilps_map = "EpiLPS (MAP)",
  epilps_mala = "EpiLPS (MALA)"
)
