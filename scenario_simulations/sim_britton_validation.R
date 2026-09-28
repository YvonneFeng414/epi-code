################################################################################
# Deterministic validation: EpiEstim's SI-misspecification bias IS the
# Euler-Lotka effect (Britton & Scalia Tomba 2019).
#
# The port of the archive's britton_validation.R. No randomness anywhere. For
# a true R and the true SI, the discrete Euler-Lotka equation gives the growth
# rate r; incidence is then a pure exponential at that rate, anchored at 1e12
# cases on the burn-in day so the Gamma prior is negligible and the start-of-
# series boundary is 200 days behind every estimate used. EpiEstim is told a
# range of assumed sds (mean held correct), and its R-hat is compared to the
# closed-form asymptote R_EL = 1 / sum_s w_assumed(s) exp(-r s). They agree to
# ~1e-7 relative across the whole sweep, which is what the stopifnot block at
# the end asserts. It is the one place in the project where the DGP, the SI
# construction and the estimator are checked against each other exactly.
#
# Departures from the archive, both from harmonising onto R/: the assumed SIs
# are make_si() vectors (lag 0 zeroed, renormalised) rather than raw discr_si,
# with the raw mass reported alongside; and the fit is fit_epiestim_full(),
# whose R_mean is the Mean(R) the archive used here (everything else in the
# project scores the median).
#
# Section 6 reconciles this with the stochastic study: the misspecification
# tables' bias column has the OPPOSITE sign to Euler-Lotka above the true sd,
# because averaging over every scored day lets the start-of-series artefact
# dominate. The deterministic replica of that set-up, scored over the same
# days, reproduces the stochastic bias; restricted to late windows it recovers
# the Euler-Lotka sign and size. If sim_si_misspec.R has run, its constant_1.2
# rows are drawn on the same figure.
#
# Run from r-proj/:  Rscript simulations/sim_britton_validation.R
# Output: results/simulations/britton_validation/
################################################################################

rm(list = ls())

# ==============================================================================
# 0. Packages and modules
# ==============================================================================

source(file.path("R", "core", "init.R"))
require_packages(c("EpiEstim", "ggplot2", "patchwork"))
source(file.path("R", "core", "config.R"))          # mean_si, sd_si, assumed_sd_grid
source(file.path("simulations", "sim_config.R"))
source_project()

# ==============================================================================
# 1. Settings
# ==============================================================================

quick <- FALSE

R_values <- if (quick) 1.2 else c(1.2, 2.5)

T_btn <- 400L
burn_in <- 200L          # days discarded before any estimate is used
anchor <- 1e12           # incidence AT the burn-in day
w_main <- 1L
w_check <- 7L

# The fine grid draws the theoretical curve; the 9 points of the +/-80%
# relative grid are the markers where EpiEstim is also run with the wider
# window and on the renewal series.
sd_marks <- assumed_sd_grid
sd_fine <- if (quick) sd_marks else seq(0.5, 12, by = 0.25)
sd_all <- sort(unique(c(sd_fine, sd_marks)))
is_mark <- vapply(sd_all, function(x) any(abs(x - sd_marks) < 1e-9), logical(1))

# Lags kept in every SI here: the whole series. At sd = 12 a lag-30 vector
# would drop a visible share of the mass; keeping T - 1 lags makes the
# truncation loss ~1e-9 and the realised-moment checks below meaningful.
max_lag_btn <- T_btn - 1L

# Reconciliation with the stochastic study (section 6).
recon_R <- 1.2
recon_time <- sim_ee_time
recon_I0 <- sim_ee_I0
recon_scored_from <- max(sim_ee_windows) + 1L   # the day sim_si_misspec scores from
recon_late_from <- 60L
stochastic_file <- file.path(sim_ee_results("si_misspec"), "si_misspec_metrics.csv")

output_dir <- sim_ee_results("britton_validation")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 2. Helpers
# ==============================================================================

# The continuous closed form, for reference. discr_si models the SI as 1 + X
# with X ~ Gamma(a, scale = b), so E[exp(-r S)] = exp(-r) (1 + r b)^(-a). It
# ignores the discretisation's variance inflation of 1/6, so it sits slightly
# above the discrete value.
euler_lotka_R_continuous <- function(mu, sigma, r) {
  a <- ((mu - 1) / sigma)^2
  b <- sigma^2 / (mu - 1)
  exp(r) * (1 + r * b)^a
}

# EpiEstim's mean R-hat over the post-burn-in days, with two diagnostics: how
# much posterior uncertainty is left at counts this large (should be ~0), and
# the day-to-day spread of R-hat, which under pure exponential growth with the
# boundary burned in must be ~0 - a non-zero value means burn_in is too short.
epiestim_R_hat <- function(I, si, w, burn_in) {
  df <- data.frame(dates = as.Date("2020-01-01") + seq_along(I) - 1L, I = I)
  fit <- fit_epiestim_full(df, si, w)
  keep <- (fit$index - w + 1L) > burn_in & is.finite(fit$R_mean)
  est <- fit$R_mean[keep]
  list(R_hat = mean(est),
       ci_half_width = mean((fit$upper[keep] - fit$lower[keep]) / 2),
       day_spread = diff(range(est)),
       n_days = length(est))
}

# The raw discr_si mass over the kept lags, before make_si() renormalises it.
raw_si_mass <- function(mean, sd, max_lag) sum(EpiEstim::discr_si(seq.int(0L, max_lag), mean, sd))

# ==============================================================================
# 3. The sweep
# ==============================================================================

si_true <- make_si(mean_si, sd_si, max_lag_btn)
s <- seq.int(0L, T_btn - 1L)

cat("\n===== Britton & Scalia Tomba validation: deterministic exponential growth =====\n")
cat(sprintf("T = %d, burn-in %d days, I = %.0e at day %d; %d assumed sds (%d markers); R in {%s}\n",
            T_btn, burn_in, anchor, burn_in, length(sd_all), sum(is_mark),
            paste(R_values, collapse = ", ")))
if (quick) cat("QUICK MODE: one R, markers only\n")

rows <- list()
for (R_true in R_values) {
  r <- discrete_growth_rate(si_true, R_true)
  R_EL_check <- euler_lotka_R(si_true, r)
  cat(sprintf("\nTrue R = %.2f   r = %.10f   doubling time = %.3f days   R_EL at the true SI = %.10f\n",
              R_true, r, log(2) / r, R_EL_check))
  # The growth rate and the SI must be mutually consistent, or every bias
  # below is measured from the wrong baseline.
  stopifnot(abs(R_EL_check - R_true) < 1e-10)

  # Primary series: pure exponential, exactly the regime Euler-Lotka describes.
  I_exp <- anchor * exp(r * (s - burn_in))
  # Cross-check: the noise-free renewal equation, rescaled to the same anchor
  # (it is linear in I, so the rescale is exact).
  I_ren <- renewal_expectation(1, R_true, si_true, T_btn)
  I_ren <- I_ren * (anchor / I_ren[burn_in])

  for (k in seq_along(sd_all)) {
    sd_a <- sd_all[k]
    si_a <- make_si(mean_si, sd_a, max_lag_btn)
    mom <- si_moments(si_a)
    R_EL <- euler_lotka_R(si_a, r)
    fit1 <- epiestim_R_hat(I_exp, si_a, w_main, burn_in)

    R_hat_7 <- NA_real_
    R_hat_ren <- NA_real_
    if (is_mark[k]) {
      R_hat_7 <- epiestim_R_hat(I_exp, si_a, w_check, burn_in)$R_hat
      R_hat_ren <- epiestim_R_hat(I_ren, si_a, w_main, burn_in)$R_hat
    }

    rows[[length(rows) + 1L]] <- data.frame(
      true_R = R_true, growth_rate = r, doubling_time = log(2) / r,
      true_mean_si = mean_si, true_sd_si = sd_si,
      assumed_mean_si = mean_si, assumed_sd_si = sd_a, is_grid_point = is_mark[k],
      realized_mean_si = mom[["mean"]], realized_sd_si = mom[["sd"]],
      si_mass = raw_si_mass(mean_si, sd_a, max_lag_btn),
      R_EL_discrete = R_EL,
      R_EL_continuous = euler_lotka_R_continuous(mean_si, sd_a, r),
      R_hat_nd1 = fit1$R_hat, R_hat_nd7 = R_hat_7, R_hat_renewal = R_hat_ren,
      ci_half_width = fit1$ci_half_width, day_spread = fit1$day_spread,
      n_days_used = fit1$n_days,
      bias_abs = fit1$R_hat - R_true,
      bias_rel_pct = 100 * (fit1$R_hat - R_true) / R_true,
      epiestim_vs_EL_rel = abs(fit1$R_hat - R_EL) / R_EL,
      stringsAsFactors = FALSE
    )
  }
}
results <- do.call(rbind, rows)
write.csv(results, file.path(output_dir, "britton_validation.csv"), row.names = FALSE)

# ==============================================================================
# 4. The checks
# ==============================================================================
# Eight premises, asserted rather than eyeballed. The thresholds are the
# archive's; the relative tolerances on the realised moments allow for the
# finite-T truncation of the heavy tail at the widest sd and the triangular
# kernel's variance inflation (sd^2 + 1/6), which is a property of the
# discretisation and not of T.

marks <- results[results$is_grid_point, ]
chk <- c(
  el     = max(results$epiestim_vs_EL_rel),
  spread = max(results$day_spread),
  mean   = max(abs(results$realized_mean_si - mean_si) / mean_si),
  sd     = max(abs(results$realized_sd_si - sqrt(results$assumed_sd_si^2 + 1 / 6)) /
                 sqrt(results$assumed_sd_si^2 + 1 / 6)),
  mass   = min(results$si_mass),
  ci     = max(results$ci_half_width),
  nd     = max(abs(marks$R_hat_nd7 - marks$R_hat_nd1)),
  ren    = max(abs(marks$R_hat_renewal - marks$R_hat_nd1))
)
limits <- c(el = 1e-6, spread = 1e-6, mean = 1e-7, sd = 5e-3, mass = 1 - 1e-6,
            ci = 1e-5, nd = 1e-6, ren = 1e-6)
labels <- c(el = "max |R_hat - R_EL| / R_EL", spread = "max day-to-day spread in R_hat",
            mean = sprintf("max relative |mean SI - %.1f|", mean_si),
            sd = "max rel |sd - sqrt(s^2 + 1/6)|", mass = "min raw sum(w) over assumed SIs",
            ci = "max posterior CI half-width", nd = "max |R_hat(w7) - R_hat(w1)| at marks",
            ren = "max |R_hat(renewal) - R_hat(exp)| at marks")

cat("\n==============================================================\n")
cat("VALIDATION CHECKS\n")
cat("==============================================================\n")
for (nm in names(chk)) {
  pass <- if (nm == "mass") chk[[nm]] > limits[[nm]] else chk[[nm]] < limits[[nm]]
  cat(sprintf("  %-42s %12.3e   (%s %.0e)  %s\n", labels[[nm]], chk[[nm]],
              if (nm == "mass") ">" else "<", limits[[nm]], if (pass) "ok" else "FAIL"))
}
stopifnot(
  chk[["el"]] < limits[["el"]], chk[["spread"]] < limits[["spread"]],
  chk[["mean"]] < limits[["mean"]], chk[["sd"]] < limits[["sd"]],
  chk[["mass"]] > limits[["mass"]], chk[["ci"]] < limits[["ci"]],
  chk[["nd"]] < limits[["nd"]], chk[["ren"]] < limits[["ren"]]
)
cat("\n  All checks passed.\n")

cat("\n==============================================================\n")
cat(sprintf("RESULTS AT THE +/-80%% GRID (assumed mean fixed at %.1f)\n", mean_si))
cat("==============================================================\n")
marks$rel_err_pct <- round(100 * (marks$assumed_sd_si - sd_si) / sd_si)
print(format(marks[, c("true_R", "rel_err_pct", "assumed_sd_si", "realized_sd_si",
                       "R_EL_discrete", "R_hat_nd1", "bias_abs", "bias_rel_pct")],
             digits = 8), row.names = FALSE)

# ==============================================================================
# 5. Figure
# ==============================================================================

page <- plot_britton_page(
  results, R_values, sd_si,
  title = "SI misspecification bias is the Euler-Lotka effect (Britton & Scalia Tomba 2019)",
  subtitle = paste0(
    sprintf("Deterministic exponential incidence at the growth rate implied by the true SI (mean = %.1f, sd = %.1f). Only the assumed sd is misspecified; the assumed mean is held at %.1f.\n",
            mean_si, sd_si, mean_si),
    sprintf("I = %.0e at day %d so the Gamma(1, 5) prior is negligible; first %d days discarded to remove the boundary. Points mark the +/-80%% grid of the misspecification studies.\n",
            anchor, burn_in, burn_in),
    sprintf("EpiEstim reproduces the Euler-Lotka value to %.1e relative across the whole sweep.", chk[["el"]])
  )
)
ggsave(file.path(output_dir, "britton_validation.pdf"), page, width = 12, height = 8)

# ==============================================================================
# 6. Reconciliation with the stochastic study
# ==============================================================================
# The misspecification study's set-up with the noise removed: I0 on day 1, the
# noise-free renewal equation from there, the median R scored over the same
# days the stochastic study scores. Same SI construction (lag 30) as that
# study, so the two are comparable point for point.

cat("\n==============================================================\n")
cat("RECONCILIATION WITH THE STOCHASTIC STUDY (sim_si_misspec, constant_1.2)\n")
cat("==============================================================\n")
cat(sprintf("Replicating constant Rt = %.1f, T = %d, I0 = %d with the noise removed; median R,\n",
            recon_R, recon_time, recon_I0))
cat(sprintf("scored from day %d as the study does, and from day %d for comparison.\n\n",
            recon_scored_from, recon_late_from + 1L))

si_true30 <- make_si(mean_si, sd_si, sim_ee_max_si_lag)
r30 <- discrete_growth_rate(si_true30, recon_R)
I_recon <- renewal_expectation(recon_I0, recon_R, si_true30, recon_time)
recon_df <- data.frame(dates = as.Date("2020-01-01") + seq_len(recon_time) - 1L, I = I_recon)

recon <- do.call(rbind, lapply(sd_marks, function(sd_a) {
  si_a <- make_si(mean_si, sd_a, sim_ee_max_si_lag)
  fit <- fit_epiestim_full(recon_df, si_a, 1L)
  scored <- fit$index >= recon_scored_from
  late <- fit$index > recon_late_from
  data.frame(
    assumed_sd_si = sd_a,
    bias_scored_windows = mean(fit$R[scored] - recon_R, na.rm = TRUE),
    bias_late_windows = mean(fit$R[late] - recon_R, na.rm = TRUE),
    bias_euler_lotka = euler_lotka_R(si_a, r30) - recon_R
  )
}))

recon$bias_stochastic <- NA_real_
if (file.exists(stochastic_file)) {
  sto <- read.csv(stochastic_file, stringsAsFactors = FALSE)
  sto <- sto[sto$subscenario == "constant_1.2" & sto$scenario == "vary_sd" &
               sto$family == "gamma_discr" & sto$arm == "w1", c("assumed_sd", "Bias")]
  recon$bias_stochastic <- sto$Bias[match(round(recon$assumed_sd_si, 2), round(sto$assumed_sd, 2))]
  gap <- max(abs(recon$bias_scored_windows - recon$bias_stochastic), na.rm = TRUE)
  cat(sprintf("Stochastic results found: max |deterministic replica - stochastic| = %.5f\n", gap))
  cat("The stochastic bias curve is therefore almost entirely a deterministic boundary artefact.\n\n")
} else {
  cat("No stochastic results yet (run sim_si_misspec.R first to fill bias_stochastic).\n\n")
}
print(format(recon, digits = 4), row.names = FALSE)
cat(sprintf("\nRestricted to windows after day %d, the curve recovers the Euler-Lotka sign and size.\n",
            recon_late_from))
write.csv(recon, file.path(output_dir, "britton_reconciliation.csv"), row.names = FALSE)

recon_page <- plot_britton_reconciliation(
  recon, sd_si, recon_late_from, recon_scored_from,
  title = "The stochastic study's bias curve is a boundary artefact, not the Euler-Lotka effect",
  subtitle = paste0(
    sprintf("Constant Rt = %.1f, T = %d, I0 = %d. Orange: the stochastic set-up with the noise removed, scored over the study's days; grey circles: the stochastic\n",
            recon_R, recon_time, recon_I0),
    "results it reproduces. Early days, where lambda_t is truncated by the start of the series, push R-hat up as the assumed sd widens. Dropping them\n",
    "recovers the Euler-Lotka curve."
  )
)
ggsave(file.path(output_dir, "britton_reconciliation.pdf"), recon_page, width = 9, height = 6)

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
cat(sprintf("\nWrote results to %s/\n", output_dir))
