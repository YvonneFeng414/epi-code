################################################################################
# Checks for the pure-synthetic EpiEstim harness (simulations/). Where the
# archive had a function of its own, its body is carried here as the oracle
# and the shared function must reproduce it; the rest are contract checks on
# the new helpers.
#
# Run from r-proj/:  Rscript tests/test_sim_epiestim.R
# Exits non-zero on the first failure, so it can gate a run.
################################################################################

source(file.path("R", "core", "init.R"))
require_packages(c("EpiEstim", "dplyr", "purrr", "ggplot2", "patchwork"))
source(file.path("R", "core", "config.R"))
source_project()

n_fail <- 0L
check <- function(label, pass, detail = "") {
  cat(sprintf("  %-64s %s%s\n", label, if (pass) "PASS" else "FAIL",
              if (nzchar(detail)) paste0("  ", detail) else ""))
  if (!pass) n_fail <<- n_fail + 1L
  invisible(pass)
}

max_abs_diff <- function(a, b) {
  a <- as.numeric(a); b <- as.numeric(b)
  if (length(a) != length(b) || !identical(is.na(a), is.na(b))) return(Inf)
  d <- abs(a - b)
  if (all(is.na(d))) 0 else max(d, na.rm = TRUE)
}

# ==============================================================================
# 1. Rt shape constructors vs the archived ones
# ==============================================================================

cat("\n===== 1. rt_step() / rt_gradual() vs archive/Epiestim_simu/simulation_functions.R =====\n")

# archive/Epiestim_simu/simulation_functions.R
make_step_R_old <- function(R_before, R_after, change_day, time) {
  c(rep(R_before, change_day), rep(R_after, time - change_day))
}
make_gradual_R_old <- function(R_before, R_after, change_day, transition_days, time) {
  transition <- seq(from = R_before, to = R_after, length.out = transition_days + 1)[-1]
  out <- c(rep(R_before, change_day), transition, rep(R_after, time - change_day - transition_days))
  if (length(out) != time) stop("Generated R_t does not have the correct length.")
  out
}

for (time in c(100L, 250L)) {
  rt_change <- round(time / 2) - 1          # the archive: last day at the old level
  change_day <- as.integer(round(time / 2)) # ours: first day at the new level
  check(sprintf("step 1.2 -> 0.8 at T = %d identical", time),
        identical(rt_step(1.2, 0.8, change_day)(time),
                  make_step_R_old(1.2, 0.8, rt_change, time)))
  for (d in c(7L, 14L, 20L)) {
    check(sprintf("gradual 1.2 -> 0.8 over %d days at T = %d identical", d, time),
          max_abs_diff(rt_gradual(1.2, 0.8, change_day, d)(time),
                       make_gradual_R_old(1.2, 0.8, rt_change, d, time)) < 1e-12)
  }
}

# ==============================================================================
# 2. The scenario registry
# ==============================================================================

cat("\n===== 2. epiestim_scenario_table / build_epiestim_scenarios() =====\n")

tab <- epiestim_scenario_table
check("27 sub-scenarios", nrow(tab) == 27L)
check("archive order: constants, steps, gradual 7 / 14 / 20",
      identical(tab$subscenario[c(1, 4, 10, 16, 22)],
                c("constant_1.0", "step_1.2_to_0.8", "gradual7_1.2_to_0.8",
                  "gradual14_1.2_to_0.8", "gradual20_1.2_to_0.8")))
check("constant 1.0 displays as \"1\" (the archive's label)",
      tab$true_rt[tab$subscenario == "constant_1.0"] == "1")
check("ASCII arrow in the pair labels",
      all(grepl(" -> ", tab$true_rt[tab$shape != "constant"], fixed = TRUE)))
check("transition_duration is \"--\" / \"1\" / \"7\" / \"14\" / \"20\"",
      identical(sort(unique(tab$transition_duration)), sort(c("--", "1", "7", "14", "20"))))

scen <- build_epiestim_scenarios(tab, 100L)
check("registry names are the subscenarios", identical(names(scen), tab$subscenario))
paths <- lapply(scen, build_rt_scenario, time = 100L)
check("every path has length 100 and is positive",
      all(vapply(paths, function(p) length(p) == 100L && all(p > 0), logical(1))))
check("step changes on day 50", paths$step_1.2_to_0.8[49] == 1.2 && paths$step_1.2_to_0.8[50] == 0.8)
check("gradual7 reaches 0.8 on day 56 and not before",
      paths$gradual7_1.2_to_0.8[56] == 0.8 && paths$gradual7_1.2_to_0.8[55] > 0.8)
check("scenario_truth_table() accepts the registry",
      nrow(scenario_truth_table(scen[1:3], 100L, c(1L, 7L))) == 3L * 2L * 100L)

# ==============================================================================
# 3. renewal_expectation() vs the archived deterministic_renewal()
# ==============================================================================

cat("\n===== 3. renewal_expectation() vs britton_validation.R's deterministic_renewal() =====\n")

# archive/Epiestim_simu/britton_validation.R
deterministic_renewal_old <- function(R, w, time, I0) {
  I <- numeric(time)
  I[1] <- I0
  for (t in 2:time) {
    max_s <- min(t, length(w))
    I[t] <- R * sum(I[t - seq_len(max_s) + 1] * w[seq_len(max_s)])
  }
  I
}

si <- make_si(mean_si, sd_si, 30L)
for (R in c(0.8, 1.2)) {
  new <- renewal_expectation(2000, R, si, 100L)
  old <- deterministic_renewal_old(R, si, 100L, 2000)
  check(sprintf("R = %.1f, I0 = 2000, T = 100: relative difference < 1e-12", R),
        max(abs(new - old) / pmax(abs(old), 1)) < 1e-12)
}
check("vector R_true accepted", length(renewal_expectation(100, rep(1.1, 60), si)) == 60L)

# ==============================================================================
# 4. Per-arm scoring helpers
# ==============================================================================

cat("\n===== 4. score_si_grid_per_arm() / score_fits_per_arm() / summarise_trajectories() =====\n")

scen100 <- build_epiestim_scenarios(epiestim_scenario_table, 60L)
target <- sim_target("gradual7_1.2_to_0.8", scen100, 60L, 2000L, c(1L, 7L))
check("sim_target scores from max(window) + 1", target$score_from == 8L)
check("sim_target's w1 truth is the path, w7 the trailing mean",
      identical(target$truth_by_arm$w1, target$R_true) &&
        max_abs_diff(target$truth_by_arm$w7, trailing_mean(target$R_true, 7L)) == 0)

series <- simulate_sim_target(target, si, 4L, seed = 1L)
grid <- build_si_misspec_grid("gamma_discr", mean_si, sd_si,
                              round(mean_si * (1 + c(-0.4, 0, 0.4)), 2),
                              round(sd_si * (1 + c(-0.4, 0, 0.4)), 2),
                              c(-0.4, 0, 0.4), 30L, 0.25)
fits_by_rep <- lapply(series, function(s) {
  lapply(grid$si_list, function(si_a) fit_epiestim_grid(s, si_a, windows = c(1L, 7L)))
})

per_arm <- score_si_grid_per_arm(fits_by_rep, grid, target$truth_by_arm, target$score_days)
manual <- lapply(c("w1", "w7"), function(a) {
  score_si_grid(fits_by_rep, grid, target$truth_by_arm[[a]][target$score_days],
                target$score_days, arms = a)
})
check("score_si_grid_per_arm() == two score_si_grid() calls bound",
      identical(per_arm$metrics_by_key, do.call(rbind, lapply(manual, `[[`, "metrics_by_key"))) &&
        identical(per_arm$by_day_by_key, do.call(rbind, lapply(manual, `[[`, "by_day_by_key"))))

correct_key <- grid$settings$key[grid$settings$correct][1L]
fits_correct <- lapply(fits_by_rep, function(rep) rep[[correct_key]])
sc <- score_fits_per_arm(fits_correct, target)
check("score_fits_per_arm(): one metrics row per arm, arm / window columns",
      nrow(sc$metrics) == 2L && identical(sc$metrics$arm, c("w1", "w7")) &&
        identical(sc$metrics$window, c(1L, 7L)))
check("score_fits_per_arm() matches score_si_grid_per_arm() at the correct key",
      max_abs_diff(sc$metrics$Coverage95,
                   per_arm$metrics_by_key$Coverage95[per_arm$metrics_by_key$key == correct_key]) < 1e-12)

traj <- summarise_trajectories(arm_fits(fits_correct, "w1"), target$truth_by_arm$w1, target$score_days)
by_day <- score_interval_fits(arm_fits(fits_correct, "w1"),
                              target$truth_by_arm$w1[target$score_days], target$score_days)$by_day
check("summarise_trajectories(): coverage equals score_interval_fits() by-day",
      max_abs_diff(traj$coverage, by_day$coverage) < 1e-12)
check("summarise_trajectories(): mean_est equals by-day mean_R",
      max_abs_diff(traj$mean_est, by_day$mean_R) < 1e-12)
check("summarise_trajectories(): n_used is the replicate count", all(traj$n_used == 4L))
check("summarise_trajectories(): mc_lower <= mean_est <= mc_upper",
      all(traj$mc_lower <= traj$mean_est + 1e-12 & traj$mean_est <= traj$mc_upper + 1e-12))
single <- fit_trajectory(fits_correct[[1L]]$w1, target$truth_by_arm$w1, target$score_days)
check("fit_trajectory(): one row per scored day, est between bounds",
      nrow(single) == length(target$score_days) &&
        all(single$lower <= single$est & single$est <= single$upper))
check("pick_representative() returns an index into the series",
      pick_representative(series) %in% seq_along(series))

# ==============================================================================
# 5. convergence_day()
# ==============================================================================

cat("\n===== 5. convergence_day() =====\n")

day <- 1:10
check("el = 0 is NA", is.na(convergence_day(day, rep(0.1, 10), 0)))
check("never within tolerance is NA", is.na(convergence_day(day, rep(1, 10), 0.1)))
check("within tolerance throughout returns the first usable day",
      convergence_day(day, rep(0.1, 10), 0.1) == 1L)
check("crosses the band and leaves: day after the LAST violation",
      convergence_day(day, c(1, 0.1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1), 0.1) == 5L)
check("last day violates: NA", is.na(convergence_day(day, c(rep(0.1, 9), 1), 0.1)))
check("NA days are skipped", convergence_day(day, c(NA, NA, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1), 0.1) == 4L)

# ==============================================================================
# 6. open_figure_device()
# ==============================================================================

cat("\n===== 6. open_figure_device() =====\n")

for (ext in c("pdf", "png")) {
  f <- tempfile(fileext = paste0(".", ext))
  open_figure_device(f, 600, 400, 100)
  plot(1:3)
  invisible(dev.off())
  check(sprintf("writes a non-empty .%s", ext), file.exists(f) && file.info(f)$size > 0)
  unlink(f)
}

# ==============================================================================
# 7. format_scenario_table()
# ==============================================================================

cat("\n===== 7. format_scenario_table() =====\n")

fixture <- data.frame(
  major_scenario = c("Scenario 1: Constant Rt", "Scenario 1: Constant Rt",
                     "Scenario 1: Constant Rt", "Scenario 2: Stepwise Rt"),
  true_rt = c("1.2", "1.2", "1.2", "1.2 -> 0.8"),
  T = c(100L, 100L, 100L, 100L), I0 = c(100L, 100L, 2000L, 2000L),
  transition_duration = c("--", "--", "--", "1"),
  window = c(1L, 7L, 1L, 1L),
  Bias = c(0.0123, -0.0004, 0.5, 1e-5), MSE = c(0.01, 0.02, 0.03, 0.04),
  Coverage95 = c(0.951, 0.9, 0.95, 0.94), MeanCIWidth = c(1.1111, 0.5, 0.4, 0.3),
  stringsAsFactors = FALSE
)
tbl <- format_scenario_table(fixture)
check("two section header rows", sum(tbl$T == "" & tbl$Window == "") == 2L)
check("header labels", identical(tbl[["True Rt"]][tbl$Window == ""], c("Constant Rt", "Stepwise Rt")))
body <- tbl[tbl$Window != "", ]
check("True Rt blanked on repeated rows within a section",
      identical(body[["True Rt"]], c("1.2", "", "", "1.2 -> 0.8")))
check("I0 shown when it changes, blanked when it repeats",
      identical(body$I0, c("100", "", "2000", "2000")))
check("Bias in 3-significant-digit scientific notation",
      all(grepl("^-?[0-9]\\.[0-9]{3}E[+-][0-9]{2}$", body$Bias)))
check("Coverage rate to 3 dp", identical(body[["Coverage rate"]], c("0.951", "0.900", "0.950", "0.940")))
check("Window labels", identical(body$Window, c("Daily", "Weekly", "Daily", "Daily")))
check("nine ASCII column names",
      identical(names(tbl), c("True Rt", "T", "I0", "Transition duration", "Window",
                              "Bias", "MSE", "Coverage rate", "CI width")))

# ==============================================================================
cat("\n==============================================================\n")
if (n_fail > 0L) {
  cat(sprintf("%d CHECK(S) FAILED\n", n_fail))
  cat("==============================================================\n")
  quit(status = 1L)
}
cat("All checks passed.\n")
cat("==============================================================\n")
