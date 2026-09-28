################################################################################
# Per-day bias in EpiEstim's Rt when only the assumed SI sd is wrong, run long
# enough to see where it lands.
#
# The port of the archive's bias-convergence diagnostic. One scenario, one
# sweep, one arm: constant Rt = 1.2, the assumed sd at -80%, -40%, 0, +40%,
# +80% of the truth with the assumed mean held correct, 1-day window, T = 200.
#
# The question it answers: the misspecification studies average bias over
# every scored day, and their bias column carries the OPPOSITE sign to the one
# theory predicts above the true sd. The reason is the start of the series -
# lambda_t is truncated there, a wider assumed SI loses more of it, and R-hat
# is pushed up on the early days by far more than the asymptotic effect. This
# study asks how many days that artefact takes to decay, and whether what is
# left converges on the Euler-Lotka value R_EL(assumed sd) - R_true that
# sim_britton_validation.R derives deterministically.
#
# Scoring starts on day 2 here, not day max(window) + 1: the early days are the
# subject. Convergence is judged on the bias DIFFERENCED against the correctly
# specified curve: all five fits share the same epidemics, so their Monte Carlo
# error is common and cancels, leaving the part actually caused by the wrong
# sd. converged_from_day is the first day from which that stays within 10% of
# the asymptote for the rest of the series (convergence_day() in
# R/studies/sim_epiestim.R).
#
# Run from r-proj/:  Rscript simulations/sim_bias_convergence.R
# Output: results/simulations/bias_convergence/
################################################################################

rm(list = ls())

# ==============================================================================
# 0. Packages and modules
# ==============================================================================

source(file.path("R", "core", "init.R"))
require_packages(c("EpiEstim", "ggplot2", "patchwork", "parallel"))
source(file.path("R", "core", "config.R"))
source(file.path("simulations", "sim_config.R"))
source_project()

# ==============================================================================
# 1. Settings
# ==============================================================================

quick <- FALSE

n_sim <- if (quick) sim_ee_quick_n_sim else sim_ee_n_sim

R_const <- 1.2
time <- 200L
I0 <- sim_ee_I0
windows <- 1L
max_si_lag <- sim_ee_max_si_lag

# The sd sweep: the trajectory subset of relative_error, mean held at the truth.
errors <- sim_ee_traj_errors
assumed_sds <- round(sd_si * (1 + errors), 2)

# Days tabulated on the console, and where the zoomed panels start.
report_days <- c(10L, 20L, 30L, 50L, 80L, 100L, 150L, 200L)
zoom_from <- 30L

output_dir <- sim_ee_results("bias_convergence")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

n_cores <- default_n_cores(cap = 4L, reserve = 1L, floor = 2L)

# ==============================================================================
# 2. True SI, assumed SIs, and the Euler-Lotka asymptotes
# ==============================================================================
# r is the growth rate the TRUE SI implies at R = 1.2; the asymptote for each
# assumed sd is the R that SI reports at that same r, minus the truth. The
# correctly specified sd gives exactly zero, which is checked.

si_true <- make_si(mean_si, sd_si, max_si_lag)
true_moments <- check_si_moments(si_true, mean_si)

si_list <- lapply(assumed_sds, function(s) make_si(mean_si, s, max_si_lag))
names(si_list) <- sprintf("sd_%.2f", assumed_sds)

growth_rate <- discrete_growth_rate(si_true, R_const)
el_bias <- vapply(si_list, function(si) euler_lotka_R(si, growth_rate) - R_const, numeric(1))
stopifnot(abs(el_bias[errors == 0]) < 1e-10)

cat("\n===== Per-day bias convergence: constant Rt = 1.2, assumed sd wrong =====\n")
cat(sprintf("True SI: realised mean/sd %.4f / %.4f (requested %.1f / %.1f), lags 0-%d\n",
            true_moments[["mean"]], true_moments[["sd"]], mean_si, sd_si, max_si_lag))
cat(sprintf("T = %d, I0 = %d, %d replicates, 1-day window; r = %.6f per day (doubling %.1f d)\n",
            time, I0, n_sim, growth_rate, log(2) / growth_rate))
cat("Euler-Lotka asymptotes (R_EL - 1.2) by assumed sd:\n")
print(data.frame(relative_error = errors, assumed_sd = assumed_sds, el_bias = el_bias),
      row.names = FALSE, digits = 6)
if (quick) cat("QUICK MODE: reduced replicates\n")

# ==============================================================================
# 3. Simulate and fit
# ==============================================================================

scenarios <- build_epiestim_scenarios(epiestim_scenario_table, time)
target <- sim_target("constant_1.2", scenarios, time, I0, windows)   # scores from day 2

series <- simulate_sim_target(target, si_true, n_sim, seed = sim_ee_seed + 1L)
cat(sprintf("median day-%d incidence %.0f\n", time,
            stats::median(vapply(series, function(x) x[time], numeric(1)))))

fits_by_rep <- fit_si_grid(series, si_list, fit_epiestim_grid, n_cores,
                           fit_args = list(windows = windows), packages = "EpiEstim")

# ==============================================================================
# 4. Per-day bias, differenced against the correct sd, and convergence
# ==============================================================================

truth <- target$truth_by_arm[["w1"]]
days <- target$score_days

rows <- lapply(seq_along(si_list), function(k) {
  s <- summarise_trajectories(arm_fits(lapply(fits_by_rep, function(rep) rep[[k]]), "w1"),
                              truth, days)
  data.frame(relative_error = errors[k], assumed_sd_si = assumed_sds[k],
             t_end = s$day, true_R = s$true_R, mean_est = s$mean_est,
             bias = s$mean_est - s$true_R, el_bias = el_bias[[k]], n_used = s$n_used)
})
conv_df <- do.call(rbind, rows)

ref <- conv_df[conv_df$relative_error == 0, c("t_end", "bias")]
names(ref)[2L] <- "bias_ref"
conv_df <- merge(conv_df, ref, by = "t_end", all.x = TRUE)
conv_df$bias_rel <- conv_df$bias - conv_df$bias_ref

conv_summary <- do.call(rbind, lapply(seq_along(assumed_sds), function(k) {
  d <- conv_df[conv_df$assumed_sd_si == assumed_sds[k], ]
  d <- d[order(d$t_end), ]
  last <- d[d$t_end == time, ]
  data.frame(
    relative_error = errors[k], assumed_sd_si = assumed_sds[k], el_bias = el_bias[[k]],
    raw_bias_last = last$bias, rel_bias_last = last$bias_rel,
    raw_gap_to_EL = last$bias - el_bias[[k]], rel_gap_to_EL = last$bias_rel - el_bias[[k]],
    converged_from_day = convergence_day(d$t_end, d$bias_rel, el_bias[[k]])
  )
}))
conv_df$converged_from_day <- conv_summary$converged_from_day[
  match(conv_df$assumed_sd_si, conv_summary$assumed_sd_si)]

conv_df <- conv_df[order(conv_df$relative_error, conv_df$t_end),
                   c("relative_error", "assumed_sd_si", "t_end", "true_R", "mean_est",
                     "bias", "bias_ref", "bias_rel", "el_bias", "converged_from_day", "n_used")]
rownames(conv_df) <- NULL
write.csv(conv_df, file.path(output_dir, "bias_convergence.csv"), row.names = FALSE)
write.csv(conv_summary, file.path(output_dir, "bias_convergence_summary.csv"), row.names = FALSE)

# ==============================================================================
# 5. Report
# ==============================================================================

day_table <- function(col) {
  d <- conv_df[conv_df$t_end %in% report_days, c("t_end", "assumed_sd_si", col)]
  wide <- stats::reshape(d, idvar = "t_end", timevar = "assumed_sd_si", direction = "wide")
  names(wide) <- sub(paste0("^", col, "\\."), "sd ", names(wide))
  wide[order(wide$t_end), ]
}

cat("\n==============================================================\n")
cat("Raw bias (mean estimate - 1.2) on selected days, by assumed sd\n")
cat("==============================================================\n")
print(format(day_table("bias"), digits = 4), row.names = FALSE)

cat("\n==============================================================\n")
cat("Bias minus the correctly-specified curve, by assumed sd\n")
cat("==============================================================\n")
print(format(day_table("bias_rel"), digits = 4), row.names = FALSE)

cat("\n==============================================================\n")
cat(sprintf("At day %d, against the Euler-Lotka asymptote\n", time))
cat("==============================================================\n")
print(format(conv_summary, digits = 4), row.names = FALSE)
cat("\n  converged_from_day: first day from which the differenced bias stays within\n")
cat("  10% of the asymptote to the end of the series (NA for the correct sd, whose\n")
cat("  asymptote is zero). A late or NA day with a small rel_gap_to_EL means the\n")
cat("  Monte Carlo error alone exceeds the tolerance, not that the bias missed.\n")

# ==============================================================================
# 6. Figure
# ==============================================================================

palette <- rel_error_palette(errors)
plot_df <- data.frame(day = conv_df$t_end, bias = conv_df$bias, bias_rel = conv_df$bias_rel,
                      error_label = rel_error_factor(conv_df$relative_error, errors))
el_df <- data.frame(el_bias = el_bias, colour = palette[rel_error_labels(errors)],
                    stringsAsFactors = FALSE)

page <- plot_bias_convergence_page(
  plot_df, el_df, palette, zoom_from = zoom_from,
  title = sprintf("Per-day bias in Rt when only the assumed SI sd is wrong  (constant Rt = %.1f)", R_const),
  subtitle = paste0(
    sprintf("Epidemics simulated with the TRUE serial interval (mean = %.1f, sd = %.1f); only the sd assumed by EpiEstim is wrong.  I0 = %d, T = %d, n_sim = %d, 1-day window.\n",
            mean_si, sd_si, I0, time, n_sim),
    "Left: the first days carry a large positive artefact - lambda_t is still truncated by the start of the series, and a wider assumed SI loses more of it.\n",
    "Middle: the raw bias once that has decayed; dotted = the Euler-Lotka value. Right: minus the correctly-specified curve, which cancels the Monte Carlo error\n",
    "common to all five fits and leaves the part caused by the wrong sd - which is what Euler-Lotka predicts."
  )
)
ggsave(file.path(output_dir, "bias_convergence.pdf"), page, width = 16, height = 5.6)

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
cat(sprintf("\nWrote results to %s/\n", output_dir))
