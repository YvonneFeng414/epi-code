################################################################################
# The missing EpiEstim arm, built inside this repo's framework.
#
# WHY THIS EXISTS. The EpiEstim column of combined_si_misspec_panel comes from
# the seed study, which reports summary metrics per grid point and no per-day
# output, so the EpiEstim row of the trajectory figure could not be drawn from
# it. Rather than add a by-day writer to a study maintained outside this repo,
# this script runs EpiEstim the way the other three arms are run here, so it is
# built from the same incidence file, the same common range, the same hole
# removal and the same seeding, and produces the same by-day table.
#
# DESIGN: aligned to the seed study on the three axes that were agreed -
# I0, RNG seed and scoring start - with everything else that those three touch
# brought along, because they do not stand alone:
#
#   I0            2000 on day 1 only, the seed's sens_I0[1]. Not the 30 copied
#                 observed days this repo used before.
#   RNG seed      set.seed(123), the seed's value.
#   scoring       day 14 onward, 240 days. NOT the seed's nominal day 2: EpiLPS
#                 needs until day 14 to recover from the single-day seed (its
#                 estimate is 790 against a truth of 5.74 on day 2), and the
#                 three arms are only comparable over a shared window. EpiEstim
#                 itself returns NA for days 2-7 here, so the seed's own
#                 effective start is day 8, not day 2.
#   truth         combined/seed_aligned_truth.csv - the seed's estimate_real_rt()
#                 output with its own Rt < 0.2 hole rule. Shared by all four arms.
#   grid          +/-60%, the seed's current 7 points.
#   observation   Poisson, EpiEstim's own model and the seed's.
#   window        1 day, matching the seed's nd_values[1].
#   SI            handed in explicitly via non_parametric_si so every arm is
#                 given the same pmf. The seed uses parametric_si and lets
#                 EpiEstim build its own; that is the one residual difference and
#                 it is worth 2.2e-04 on the truth path.
#
# Run from queens_rt/:  Rscript combined/run_epiestim_arm.R
# Output: combined/epiestim_arm_by_day.csv, combined/epiestim_arm_metrics.csv
################################################################################

suppressPackageStartupMessages({ library(EpiEstim); library(parallel) })

if (!file.exists(file.path("R_epi_lsp", "utils.R"))) {
  stop("Run this from the queens_rt/ directory.")
}
source(file.path("R_epi_lsp", "utils.R"))

truth_file <- file.path("combined", "seed_aligned_truth.csv")
out_dir    <- "combined"
mean_si      <- 7.5
sd_si        <- 3.4
max_si_lag   <- 30L
n_sim        <- 100L
I0           <- 2000L      # seed sens_I0[1]: one seeded day, not copied observation
score_from   <- 14L        # earliest day all three estimators have recovered
sim_seed     <- 123L       # the seed study's
# Plausibility ceiling on a single replicate-day estimate. Inf here, so the
# committed scripts are unchanged; the score-from experiment turns it on.
#
# It exists because scoring from day 2 exposes estimates that are not merely
# imprecise but arithmetically absurd: EpiLPS returns 5.6e11 on day 2 at
# assumed sd -60%, because its B-spline has no grid cap and lambda_t is near
# zero that early. Those are finite, so na.rm cannot see them, and one of them
# moves a 252-day average by nine orders of magnitude.
#
# When set, an estimate above it takes its whole triple (point, lower, upper)
# out with it - an interval around an absurd point estimate is not evidence
# about coverage either. The rule is applied identically to all three arms;
# in practice only EpiLPS ever trips it.
#
# NOTE a ceiling of 5 also removes days where the TRUTH legitimately exceeds 5
# (the Queens plug-in truth peaks at 6.74 during the explosive phase), so it
# discards some genuinely hard days along with the absurd ones. 10 would clear
# the truth entirely. 5 is the chosen setting; the trade is recorded here.
r_max_plausible <- Inf

zero_run_threshold <- 7L
epiestim_window    <- 1L
fam <- "gamma_bin"

relative_error <- c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6)

n_cores <- max(1L, min(8L, parallel::detectCores(logical = TRUE) - 2L))

# ---- 1. Truth, shared with every other arm ------------------------------------

if (!file.exists(truth_file)) {
  stop("Cannot find ", truth_file, ". Run combined/build_seed_truth.R first.")
}
R_true       <- read.csv(truth_file, stringsAsFactors = FALSE)$true_R
n_days       <- length(R_true)
score_days   <- seq.int(score_from, n_days)
truth_scored <- R_true[score_days]

cat(sprintf("Truth: seed-aligned, %d days, Rt %.3f to %.3f\n",
            n_days, min(R_true), max(R_true)))
cat(sprintf("Seeding: I0 = %d on day 1; scoring days %d-%d (%d days)\n",
            I0, score_from, n_days, length(score_days)))

si_true <- make_si(mean_si, sd_si, max_si_lag)

cat(sprintf("Scored truth spans %.3f to %.3f\n", min(truth_scored), max(truth_scored)))

# ---- 2. Simulate: Poisson, one seeded day, the seed's RNG stream --------------

simulate_series <- function(R, si, I0) {
  n <- length(R)
  I <- numeric(n); I[1L] <- I0
  for (t in seq.int(2L, n)) {
    k  <- seq_len(min(length(si) - 1L, t - 1L))
    mu <- R[t] * sum(I[t - k] * si[k + 1L])
    if (!is.finite(mu) || mu > 1e9) return(NULL)
    I[t] <- stats::rpois(1L, lambda = max(mu, 0))
  }
  if (anyNA(I)) return(NULL)
  as.integer(I)
}

set.seed(sim_seed)
series  <- lapply(seq_len(n_sim), function(r) simulate_series(R_true, si_true, I0))
kept    <- !vapply(series, is.null, logical(1))
series  <- series[kept]
extinct <- vapply(series, has_extinction_run, logical(1), threshold = zero_run_threshold)
series  <- series[!extinct]
if (length(series) == 0L) stop("no usable replicates")
cat(sprintf("simulated %d / %d replicates (%d diverged, %d went extinct)\n",
            length(series), n_sim, sum(!kept), sum(extinct)))

# ---- 3. The assumed SIs ------------------------------------------------------

settings <- rbind(
  data.frame(scenario = "vary_sd",   assumed_mean = mean_si,
             assumed_sd = round(sd_si * (1 + relative_error), 2), stringsAsFactors = FALSE),
  data.frame(scenario = "vary_mean", assumed_mean = round(mean_si * (1 + relative_error), 2),
             assumed_sd = sd_si, stringsAsFactors = FALSE)
)
settings$rel_error  <- rep(relative_error, 2L)
settings$grid_value <- ifelse(settings$scenario == "vary_sd",
                              settings$assumed_sd, settings$assumed_mean)
settings$key <- sprintf("%s_%.2f_%.2f", fam, settings$assumed_mean, settings$assumed_sd)
settings$correct <- settings$assumed_mean == mean_si & settings$assumed_sd == sd_si

fit_keys  <- unique(settings$key)
fit_specs <- settings[match(fit_keys, settings$key), ]

# ---- 4. One EpiEstim fit -----------------------------------------------------
# t_start begins at 2 (EpiEstim conditions on day 1) and the window is 1 day, so
# t_end indexes the day the estimate refers to.

fit_one <- function(incidence, si) {
  t_start <- seq.int(2L, length(incidence) - epiestim_window + 1L)
  t_end   <- t_start + epiestim_window - 1L
  f <- tryCatch(
    EpiEstim::estimate_R(
      incid = incidence, method = "non_parametric_si",
      config = EpiEstim::make_config(list(t_start = t_start, t_end = t_end,
                                          si_distr = si))),
    error = function(e) NULL)
  if (is.null(f)) return(NULL)
  idx <- match(score_days, f$R$t_end)
  if (anyNA(idx)) return(NULL)
  list(R  = f$R$`Median(R)`[idx],
       lo = f$R$`Quantile.0.025(R)`[idx],
       hi = f$R$`Quantile.0.975(R)`[idx])
}

# ---- 5. Run the grid ---------------------------------------------------------

cat(sprintf("fitting %d settings x %d replicates on %d cores",
            nrow(fit_specs), length(series), n_cores))
t0 <- proc.time()

by_day_list <- list(); metrics_list <- list()

for (i in seq_len(nrow(fit_specs))) {
  si_assumed <- make_si_family(fam, fit_specs$assumed_mean[i],
                               fit_specs$assumed_sd[i], max_si_lag)
  fits <- parallel::mclapply(series, function(I) fit_one(I, si_assumed),
                             mc.cores = n_cores)
  ok <- !vapply(fits, is.null, logical(1))
  fits <- fits[ok]
  if (length(fits) == 0L) stop("every fit failed at key ", fit_specs$key[i])

  est <- vapply(fits, `[[`, numeric(length(truth_scored)), "R")
  lo  <- vapply(fits, `[[`, numeric(length(truth_scored)), "lo")
  hi  <- vapply(fits, `[[`, numeric(length(truth_scored)), "hi")

  # Plausibility ceiling, per replicate-day. Blanks the whole triple so the
  # day leaves the point estimate AND the coverage tally.
  absurd <- !is.na(est) & est > r_max_plausible
  est[absurd] <- NA_real_; lo[absurd] <- NA_real_; hi[absurd] <- NA_real_

  # na.rm throughout: EpiEstim returns NA for every day at or before the
  # realised mean of the ASSUMED SI (its own `t_end > final_mean_si` rule), so
  # how many leading days are missing slides with the grid point - 2 days at
  # assumed mean 3.0, 10 at 12.0. At score_from = 14 nothing is missing and
  # these are no-ops; below that they are what keeps a missing day out of the
  # average instead of turning the whole metric into NA.
  day_tbl <- data.frame(
    key = fit_specs$key[i], day = score_days, truth = truth_scored,
    mean_R = rowMeans(est, na.rm = TRUE),
    sd_R = apply(est, 1L, stats::sd, na.rm = TRUE),
    coverage = rowMeans(lo <= truth_scored & truth_scored <= hi, na.rm = TRUE),
    half = rowMeans(hi - lo, na.rm = TRUE) / 2, stringsAsFactors = FALSE)
  day_tbl$bias <- day_tbl$mean_R - day_tbl$truth
  by_day_list[[i]] <- day_tbl

  per_rep_bias <- colMeans(est - truth_scored, na.rm = TRUE)
  per_rep_cov  <- colMeans(lo <= truth_scored & truth_scored <= hi, na.rm = TRUE)
  per_rep_wid  <- colMeans(hi - lo, na.rm = TRUE)
  per_rep_sq   <- colMeans((est - truth_scored)^2, na.rm = TRUE)
  metrics_list[[i]] <- data.frame(
    key = fit_specs$key[i], n_replicates = length(fits), n_failed = sum(!ok),
    Coverage95 = mean(per_rep_cov),
    MCSE_Coverage = stats::sd(per_rep_cov) / sqrt(length(fits)),
    MeanCIWidth = mean(per_rep_wid), Bias = mean(per_rep_bias),
    MSE = mean(per_rep_sq), RMSE = sqrt(mean(per_rep_sq)), stringsAsFactors = FALSE)
  cat(".")
}
cat(sprintf(" done in %.0f s\n", (proc.time() - t0)[["elapsed"]]))

by_day  <- merge(settings[, c("key", "scenario", "grid_value", "rel_error", "correct")],
                 do.call(rbind, by_day_list), by = "key")
metrics <- merge(settings, do.call(rbind, metrics_list), by = "key")
by_day$family <- fam; metrics$family <- fam

write.csv(by_day,  file.path(out_dir, "epiestim_arm_by_day.csv"),  row.names = FALSE)
write.csv(metrics, file.path(out_dir, "epiestim_arm_metrics.csv"), row.names = FALSE)

cat("\nControl point and the +/-60% arms:\n")
sel <- metrics[abs(metrics$rel_error) %in% c(0, 0.6), ]
print(data.frame(scenario = sel$scenario, rel_error = sel$rel_error,
                 Coverage = round(sel$Coverage95, 4), Bias = round(sel$Bias, 4),
                 RMSE = round(sel$RMSE, 4), Width = round(sel$MeanCIWidth, 4)),
      row.names = FALSE)
cat(sprintf("\nWrote %s/epiestim_arm_by_day.csv and _metrics.csv\n", out_dir))
