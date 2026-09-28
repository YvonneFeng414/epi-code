# ==============================================================================
# sim_config.R
# Shared settings for the pure-synthetic EpiEstim studies in this folder. Same
# contract as R/core/config.R: settings only - no directory is created, no seed
# is set, no hardware is probed. Every driver sources R/core/config.R first
# (for the true SI mean_si / sd_si and the relative-error grids) and then this.
#
# These are the archive's config.R values for its synthetic entry points, with
# the Queens-only settings dropped. n_sim = 2000 is the archive's; the quick_*
# values are for smoke runs and are switched on by `quick <- TRUE` in a driver.
# ==============================================================================

# Series length. The main study also runs the constant shapes at the long
# horizon (the archive's time_values = c(100, 250) with stepwise / gradual
# shapes at 100 only).
sim_ee_time <- 100L
sim_ee_time_long <- 250L

# Seed cases on day 1. The main study crosses the shapes with the whole grid;
# every other study runs at the single middle value.
sim_ee_I0_grid <- c(100L, 2000L, 10000L)
sim_ee_I0 <- 2000L

sim_ee_n_sim <- 2000L

# EpiEstim arms. w1 is the archive's nd = 1 and the primary arm; w7 is the
# window main.R uses on the real data. Each is scored against its own estimand.
sim_ee_windows <- c(1L, 7L)

# Lags kept in every SI vector. At the true SI (7.5 / 3.4) lag 30 captures
# ~100% of the mass; the widest assumed sd in the sweeps (6.12) loses 0.8%,
# which make_si() renormalises away. Matched to the Queens studies.
sim_ee_max_si_lag <- 30L

# Base seed; each simulated cell adds its own index.
sim_ee_seed <- 20260907L

# The nine synthetic targets of the misspecification appendices: every
# constant level, and the six gradual pairs at the 7-day transition.
sim_ee_misspec_targets <- c(
  "constant_1.0", "constant_0.8", "constant_1.2",
  "gradual7_1.2_to_0.8", "gradual7_0.8_to_1.2", "gradual7_1.3_to_1.1",
  "gradual7_1.1_to_1.3", "gradual7_0.9_to_0.7", "gradual7_0.7_to_0.9"
)

# The trajectory figures: one shape of each kind, all 1.2 -> 0.8.
sim_ee_traj_targets <- c(
  "constant_1.2", "step_1.2_to_0.8",
  "gradual7_1.2_to_0.8", "gradual14_1.2_to_0.8", "gradual20_1.2_to_0.8"
)

# The relative errors drawn as trajectories: a subset of relative_error, since
# nine lines in one panel is an unreadable fan. Also the sd grid of the bias
# convergence study.
sim_ee_traj_errors <- c(-0.8, -0.4, 0, 0.4, 0.8)

# Where a study's outputs go: results/simulations/<name>/.
sim_ee_results <- function(name) project_paths()$results_dir(file.path("simulations", name))

# Smoke-run settings, used when a driver sets quick <- TRUE.
sim_ee_quick_n_sim <- 20L
sim_ee_quick_scenarios <- c("constant_1.2", "step_1.2_to_0.8", "gradual7_1.2_to_0.8")
sim_ee_quick_targets <- c("constant_1.2", "gradual7_1.2_to_0.8")
