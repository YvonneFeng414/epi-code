################################################################################
# Independent self-check of the EpiLPS arm, and the observation-model contrast.
#
# This script is DELIBERATELY SELF-CONTAINED. It does not source R/.
# That is the whole point: truth_source_main.R's EpiLPS arm reports 95% credible
# intervals covering their own plug-in truth only 60.6% of the time, and a check
# built from the same helpers cannot tell you whether that number is real or a
# bug in a shared function. This re-implements the simulator, the fitting loop
# and the scoring from scratch. If the two agree, neither is broken.
#
# Please do not "helpfully" refactor this to reuse R/ - the duplication
# is the experiment.
#
# It answers two questions.
#
# 1. REPLICATION. Under a Poisson observation model, with every setting matched
#    to truth_source_main.R, does an independent implementation reproduce 0.606?
#
# 2. OBSERVATION MODEL. The pipeline simulates Poisson. EpiLPS estimates a
#    negative-binomial overdispersion of about 17 on the real Queens series,
#    which at 1000 cases/day means an SD near 242 against Poisson's 32 - roughly
#    eight times the noise. The 0.606 arises because a FIXED smoothing bias
#    (~0.038) is large next to a NARROW interval (half-width ~0.036). Widen the
#    interval by simulating the overdispersion the data actually shows and the
#    bias may be swallowed. If it is, the headline result is an artefact of an
#    unrealistically quiet observation model rather than a property of EpiLPS.
#
# Run from r-proj/:  Rscript epilps_selfcheck.R
# Output: results/epilps_selfcheck/
################################################################################

rm(list = ls())

# ==============================================================================
# 0. Packages
# ==============================================================================
# EpiEstim is used ONLY for discr_si(). No estimation happens through it here.

required_packages <- c("EpiLPS", "EpiEstim", "parallel")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Install first: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(EpiLPS)
  library(EpiEstim)
  library(parallel)
})

# ==============================================================================
# 1. Settings - every one matched to truth_source_main.R
# ==============================================================================
# The Poisson arm is only a valid replication if these agree exactly with the
# pipeline. Any drift here turns a genuine cross-check into a comparison of two
# different experiments.

truth_file     <- file.path("results", "truth_source_discr_si", "plugin_truth_paths.csv")
truth_column   <- "true_R_epilps"
incidence_file <- file.path("results", "queens_rt", "queens_daily_incidence.csv")

common_start <- 16L    # first day of the common range the pipeline selected

mean_si <- 7.5
sd_si   <- 3.4
max_lag <- 30L

n_sim     <- 100L
K_epilps  <- 30L
burn_in   <- 30L       # seed days, copied from observation and never scored
sim_seed  <- 20260813L # the pipeline's epilps arm seed (sim_seed + 1)

obs_models <- c("poisson", "nbinom")

# NULL means estimate it from the observed series rather than assume a value.
overdispersion <- NULL

# The coverage the pipeline reports for this arm, for the cross-check.
pipeline_coverage <- 0.6059821

output_dir <- file.path("results", "epilps_selfcheck")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

available_cores <- parallel::detectCores(logical = TRUE)
if (is.na(available_cores)) available_cores <- 1L
n_cores <- max(2L, min(4L, available_cores - 1L))

# Validated categorical palette slots 1 and 2.
model_colors <- c(poisson = "#2a78d6", nbinom = "#eb6834")

# ==============================================================================
# 2. Serial interval, with a self-test
# ==============================================================================
# si[k + 1] = P(SI = k), so si[1] is the lag-0 cell and is zero.
#
# The assertion below is the point of this section. Pairing I[t-k] with si[k]
# instead of si[k+1] shifts the whole distribution a day later and moves the
# realised mean from 7.5 to 8.5 - a mistake that is invisible in the output but
# changes the experiment. Checking the realised mean against the requested one
# catches it immediately.

si <- EpiEstim::discr_si(seq.int(0L, max_lag), mean_si, sd_si)
si[1L] <- 0
si <- si / sum(si)

si_lags <- seq.int(0L, max_lag)
realised_mean <- sum(si_lags * si)
realised_sd <- sqrt(sum((si_lags - realised_mean)^2 * si))

if (abs(realised_mean - mean_si) > 0.05) {
  stop("Realised SI mean is ", round(realised_mean, 3), ", expected ", mean_si,
       ". The lag indexing is off by one.")
}

cat("\n===== Independent EpiLPS self-check =====\n")
cat(sprintf("SI: realised mean/sd %.4f / %.4f (requested %.1f / %.1f), lags 0-%d\n",
            realised_mean, realised_sd, mean_si, sd_si, max_lag))

# ==============================================================================
# 3. Inputs
# ==============================================================================

for (f in c(truth_file, incidence_file)) {
  if (!file.exists(f)) {
    stop("Cannot find ", normalizePath(f, mustWork = FALSE),
         ". Run truth_source_main.R first.")
  }
}

truth <- read.csv(truth_file)
R_true <- truth[[truth_column]]
n_days <- length(R_true)

obs <- read.csv(incidence_file)
seed_incidence <- obs$I[seq.int(common_start, nrow(obs))][seq_len(burn_in)]

score_from <- burn_in + 1L
score_days <- seq.int(score_from, n_days)

cat(sprintf("Truth: %s, %d days, Rt %.3f to %.3f\n",
            truth_column, n_days, min(R_true), max(R_true)))
cat(sprintf("Seed: first %d observed days of the common range (%d to %d cases/day)\n",
            burn_in, min(seed_incidence), max(seed_incidence)))
cat(sprintf("Scoring days %d-%d (%d days), %d replicates per observation model\n",
            score_from, n_days, length(score_days), n_sim))

# ==============================================================================
# 4. One EpiLPS fit, in its own process
# ==============================================================================
# EpiLPS's compiled backend corrupts memory when estimR() is called more than
# once in a session, so every call gets a fresh fork. Calling it in a plain loop
# will eventually segfault, which is why this is not optional.

fit_one <- function(incidence, si, K) {
  tryCatch({
    f <- EpiLPS::estimR(incidence = incidence, si = si[-1L], K = K)
    list(
      R     = as.numeric(f$RLPS$Rq0.50),   # median, matching the pipeline
      lower = as.numeric(f$RLPS$Rq0.025),
      upper = as.numeric(f$RLPS$Rq0.975),
      time  = as.integer(f$RLPS$Time),
      rho   = as.numeric(f$NegBinoverdisp),
      error = NA_character_
    )
  }, error = function(e) list(R = NULL, error = conditionMessage(e)))
}

fit_many <- function(series_list, si, K, cores) {
  parallel::mclapply(
    series_list,
    function(x) fit_one(x, si, K),
    mc.cores = cores,
    mc.preschedule = FALSE
  )
}

# ==============================================================================
# 5. Overdispersion, estimated from the observed series
# ==============================================================================

if (is.null(overdispersion)) {
  obs_fit <- fit_many(list(obs$I), si, K_epilps, 2L)[[1L]]
  if (!is.null(obs_fit$error) && !is.na(obs_fit$error)) {
    stop("EpiLPS failed on the observed series: ", obs_fit$error)
  }
  overdispersion <- obs_fit$rho
}

cat(sprintf("\nNegBin overdispersion estimated on the observed series: rho = %.2f\n",
            overdispersion))
cat(sprintf("  at mu = 1000: SD %.0f (NegBin) vs %.0f (Poisson), %.1fx the noise\n",
            sqrt(1000 + 1000^2 / overdispersion), sqrt(1000),
            sqrt(1000 + 1000^2 / overdispersion) / sqrt(1000)))

# ==============================================================================
# 6. Simulator
# ==============================================================================
# Renewal process. The first length(seed_incidence) days are copied verbatim and
# never scored; every later day is drawn from the chosen observation model.
#
# Note si[k + 1L]: element k + 1 of si is P(SI = k), so lag k must index k + 1.

simulate_series <- function(R, si, seed_incidence, model, rho) {
  n  <- length(R)
  ns <- length(seed_incidence)
  I  <- numeric(n)
  I[seq_len(ns)] <- seed_incidence

  for (t in seq.int(ns + 1L, n)) {
    k  <- seq_len(min(length(si) - 1L, t - 1L))
    mu <- R[t] * sum(I[t - k] * si[k + 1L])

    if (!is.finite(mu) || mu > 1e9) return(NULL)

    I[t] <- if (model == "poisson") {
      stats::rpois(1L, lambda = max(mu, 0))
    } else {
      stats::rnbinom(1L, mu = max(mu, 0), size = rho)
    }
  }

  if (anyNA(I)) return(NULL)
  as.integer(I)
}

# ==============================================================================
# 7. Run both observation models
# ==============================================================================

metrics <- list()
by_day  <- list()
examples <- list()

for (model in obs_models) {
  cat(sprintf("\n===== %s =====\n", model))

  # Simulate in the PARENT so the RNG stream is reproducible; only the fits fork.
  set.seed(sim_seed)
  series <- lapply(seq_len(n_sim), function(r) {
    simulate_series(R_true, si, seed_incidence, model, overdispersion)
  })

  kept <- !vapply(series, is.null, logical(1))
  series <- series[kept]
  cat(sprintf("simulated %d / %d replicates (%d diverged)\n",
              length(series), n_sim, sum(!kept)))
  cat(sprintf("median final-day incidence %.0f\n",
              stats::median(vapply(series, function(x) x[n_days], numeric(1)))))

  examples[[model]] <- series[[1L]]

  fits <- fit_many(series, si, K_epilps, n_cores)
  ok <- vapply(fits, function(f) is.null(f$error) || is.na(f$error), logical(1))
  if (any(!ok)) {
    cat(sprintf("EpiLPS failed on %d replicate(s); first: %s\n",
                sum(!ok), fits[[which(!ok)[1L]]]$error))
  }
  fits <- fits[ok]
  if (length(fits) == 0L) stop("Every EpiLPS fit failed for model ", model)

  # --- per-replicate coverage first, so the MCSE across replicates is valid ---
  truth_scored <- R_true[score_days]

  per_rep <- t(vapply(fits, function(f) {
    idx <- match(score_days, f$time)
    lo <- f$lower[idx]; hi <- f$upper[idx]; est <- f$R[idx]
    good <- is.finite(lo) & is.finite(hi) & is.finite(est)
    c(
      coverage = mean((lo <= truth_scored & truth_scored <= hi)[good]),
      width    = mean((hi - lo)[good]),
      bias     = mean((est - truth_scored)[good]),
      sq_err   = mean(((est - truth_scored)^2)[good])
    )
  }, numeric(4)))

  # --- day-level, across replicates ---
  est_mat <- vapply(fits, function(f) f$R[match(score_days, f$time)], numeric(length(score_days)))
  lo_mat  <- vapply(fits, function(f) f$lower[match(score_days, f$time)], numeric(length(score_days)))
  hi_mat  <- vapply(fits, function(f) f$upper[match(score_days, f$time)], numeric(length(score_days)))

  day_tbl <- data.frame(
    obs_model = model,
    day = score_days,
    truth = truth_scored,
    mean_R = rowMeans(est_mat),
    # Matrices are days x replicates, so recycling truth_scored down the columns
    # compares each replicate's interval against the right day.
    coverage = rowMeans(lo_mat <= truth_scored & truth_scored <= hi_mat),
    half = rowMeans(hi_mat - lo_mat) / 2
  )
  day_tbl$bias <- day_tbl$mean_R - day_tbl$truth
  by_day[[model]] <- day_tbl

  coverage <- mean(per_rep[, "coverage"])
  mcse <- stats::sd(per_rep[, "coverage"]) / sqrt(nrow(per_rep))
  mse <- mean(per_rep[, "sq_err"])

  metrics[[model]] <- data.frame(
    obs_model = model,
    overdispersion = if (model == "nbinom") overdispersion else NA_real_,
    n_replicates = nrow(per_rep),
    n_days = length(score_days),
    Coverage95 = coverage,
    MCSE_Coverage = mcse,
    MeanCIWidth = mean(per_rep[, "width"]),
    Bias = mean(per_rep[, "bias"]),
    MSE = mse,
    RMSE = sqrt(mse),
    days_bias_exceeds_half = sum(abs(day_tbl$bias) > day_tbl$half),
    pct_bias_exceeds_half = 100 * mean(abs(day_tbl$bias) > day_tbl$half)
  )

  cat(sprintf("coverage %.4f (MCSE %.4f)  CI width %.4f  bias %+.5f  RMSE %.5f\n",
              coverage, mcse, mean(per_rep[, "width"]), mean(per_rep[, "bias"]), sqrt(mse)))
}

metrics <- do.call(rbind, metrics)
by_day  <- do.call(rbind, by_day)

write.csv(metrics, file.path(output_dir, "selfcheck_metrics.csv"), row.names = FALSE)
write.csv(by_day, file.path(output_dir, "selfcheck_by_day.csv"), row.names = FALSE)

# ==============================================================================
# 8. The cross-check
# ==============================================================================
# This is the test that matters. The Poisson arm shares no code with the
# pipeline, so agreement is evidence that neither implementation is broken. A
# failure means one of them IS broken, and the nbinom result below cannot be
# trusted until that is resolved.

pois <- metrics[metrics$obs_model == "poisson", ]
diff <- pois$Coverage95 - pipeline_coverage
tol <- 2 * pois$MCSE_Coverage
verdict <- if (abs(diff) <= tol) "PASS" else "FAIL"

cat("\n==============================================================\n")
cat("REPLICATION CROSS-CHECK (Poisson arm vs truth_source_main.R)\n")
cat(sprintf("  independent : %.4f\n", pois$Coverage95))
cat(sprintf("  pipeline    : %.4f\n", pipeline_coverage))
cat(sprintf("  difference  : %+.4f   tolerance (2 MCSE) %.4f   -> %s\n",
            diff, tol, verdict))
if (verdict == "FAIL") {
  cat("  One of the two implementations is wrong. Do not interpret the\n")
  cat("  observation-model contrast below until this is resolved.\n")
}
cat("==============================================================\n")

cat("\n===== Observation model contrast =====\n")
print(
  metrics[, c("obs_model", "overdispersion", "Coverage95", "MCSE_Coverage",
              "MeanCIWidth", "Bias", "RMSE", "pct_bias_exceeds_half")],
  row.names = FALSE, digits = 4
)

if (nrow(metrics) == 2L) {
  cat(sprintf("\nCI width ratio nbinom/poisson: %.2fx",
              metrics$MeanCIWidth[2] / metrics$MeanCIWidth[1]))
  cat(sprintf("  (simulated noise ratio at mu=1000: %.1fx)\n",
              sqrt(1000 + 1000^2 / overdispersion) / sqrt(1000)))
}

# ==============================================================================
# 9. Figures
# ==============================================================================
# Base R, in the idiom of the script this replaces. One panel per metric with
# its own scale - never a shared axis across quantities in different units.

png(file.path(output_dir, "selfcheck_by_day.png"),
    width = 2400, height = 2400, res = 200)
op <- par(mfrow = c(3, 1), mar = c(4, 4.5, 3, 1))

plot(NA, xlim = range(by_day$day), ylim = c(0, 1),
     xlab = "Day", ylab = "Coverage across replicates",
     main = "95% CI coverage of the arm's own plug-in truth")
abline(h = 0.95, lty = 2, col = "grey35")
for (m in obs_models) {
  d <- by_day[by_day$obs_model == m, ]
  lines(d$day, d$coverage, col = model_colors[[m]], lwd = 2)
}
legend("bottomleft", legend = obs_models, col = model_colors[obs_models],
       lwd = 2, bty = "n")

plot(NA, xlim = range(by_day$day), ylim = range(c(by_day$bias, by_day$half, -by_day$half)),
     xlab = "Day", ylab = "Estimate - truth",
     main = "Day-level bias (lines) against reported CI half-width (dotted)")
abline(h = 0, lty = 2, col = "grey35")
for (m in obs_models) {
  d <- by_day[by_day$obs_model == m, ]
  lines(d$day, d$bias, col = model_colors[[m]], lwd = 2)
  lines(d$day, d$half, col = model_colors[[m]], lty = 3)
  lines(d$day, -d$half, col = model_colors[[m]], lty = 3)
}
legend("bottomleft", legend = obs_models, col = model_colors[obs_models],
       lwd = 2, bty = "n")

plot(NA, xlim = c(1, n_days), ylim = c(0, max(vapply(examples, max, numeric(1)))),
     xlab = "Day", ylab = "Incidence",
     main = "One simulated replicate per observation model")
for (m in obs_models) lines(seq_len(n_days), examples[[m]], col = model_colors[[m]])
abline(v = burn_in + 0.5, lty = 2, col = "grey35")
legend("topright", legend = c(obs_models, "scoring starts"),
       col = c(model_colors[obs_models], "grey35"), lty = c(1, 1, 2), bty = "n")

par(op)
invisible(dev.off())

cat("\nWritten to:\n"); cat(normalizePath(output_dir), "\n")
