# ============================================================
# THE FOUR ESTIMATORS
#
# Every fitter takes the T x n_sim incidence matrix and returns
# the same layout as fit_rt_trajectories() in the seed study -
#
#   list(t_start, t_end, est, lower, upper)
#
# with one row per day 2..T (1-day windows, t_start = t_end),
# so summarise_rt_trajectories() and the seed plot functions
# work on all four unchanged. Days an estimator does not report
# stay NA. Day 1 is the seed and is never estimated.
#
# Extra fields: n_failed (replicates with no estimate at all),
# errors (their messages) and secs (wall time).
#
# est is the posterior MEDIAN for all four, and lower / upper
# the equal-tailed 95% interval.
#
# si_mean / si_sd are the serial interval the estimator is TOLD
# (default: the true one). The data are always simulated with
# the true SI, so passing anything else is a misspecification.
#
# Requires config.R and the seed rt_trajectory_functions.R.
# ============================================================


# ------------------------------------------------------------
# Shared pieces
# ------------------------------------------------------------

traj_days <- function(time) seq.int(2L, time)

empty_traj <- function(time, n_sim) {
  days <- traj_days(time)
  m <- matrix(NA_real_, nrow = length(days), ncol = n_sim)
  list(t_start = days, t_end = days, est = m, lower = m, upper = m)
}

# The SI handed to EpiFilter and EpiLPS: EpiEstim's own
# discretisation (the one the data were simulated with), lag 0
# zeroed, truncated at max_lag and renormalised. Same as
# make_si() in the queens_rt scripts. The mass past lag 30 is
# below 1e-4. si[1] is lag 0; the estimators get si[-1].
make_si <- function(mean_si, sd_si, max_lag = 30L) {
  out <- EpiEstim::discr_si(seq.int(0L, max_lag), mean_si, sd_si)
  out[1L] <- 0
  out / sum(out)
}

# Run one fit per replicate in forked workers. EpiLPS's compiled
# backend has been seen to corrupt memory across repeated calls
# in one process (queens_rt/R_epi_lsp/epilps.R), so the parent
# never fits inside a loop; mc.preschedule = FALSE gives each
# replicate its own fork. mc.cores >= 2 so mclapply really forks.
run_forked <- function(n_sim, fit_one) {
  res <- parallel::mclapply(
    seq_len(n_sim),
    fit_one,
    mc.cores = max(2L, min(n_cores, n_sim)),
    mc.preschedule = FALSE
  )
  # A fork-level failure comes back as a try-error, not a list
  lapply(res, function(x) {
    if (inherits(x, "try-error") || is.null(x)) {
      list(error = paste("fork failed:", as.character(x)))
    } else {
      x
    }
  })
}

# Put per-replicate day-indexed results into the matrices
collect_fits <- function(fits, time, n_sim, secs) {
  out <- empty_traj(time, n_sim)
  errors <- character(0)

  for (i in seq_len(n_sim)) {
    f <- fits[[i]]
    if (!is.null(f$error) && !is.na(f$error)) {
      errors <- c(errors, sprintf("rep %d: %s", i, f$error))
      next
    }
    idx <- match(f$day, out$t_end)
    keep <- !is.na(idx)
    out$est[idx[keep], i]   <- f$R[keep]
    out$lower[idx[keep], i] <- f$lower[keep]
    out$upper[idx[keep], i] <- f$upper[keep]
  }

  out$n_failed <- sum(colSums(is.finite(out$est)) == 0L)
  out$errors <- errors
  out$secs <- secs
  out
}


# ============================================================
# 1. EpiEstim - the seed study's own fitter, unchanged
# ============================================================

fit_epiestim <- function(simulations, dates,
                         si_mean = mean_si, si_sd = sd_si) {
  secs <- system.time(
    traj <- fit_rt_trajectories(
      simulations = simulations,
      dates = dates,
      nd = epiestim_nd,
      mean_si = si_mean,
      sd_si = si_sd
    )
  )[["elapsed"]]

  traj$n_failed <- sum(colSums(is.finite(traj$est)) == 0L)
  traj$errors <- character(0)
  traj$secs <- secs
  traj
}


# ============================================================
# 2. EpiFilter - smoother arm (retrospective, whole series)
#
# Follows run_epifilter() in queens_rt/epifilter_si_misspec.R.
# The recursion runs from epifilter_start; earlier days stay NA.
# ============================================================

total_infectiousness <- function(incidence, si) {
  w <- si[-1L]
  n <- length(incidence)
  L <- numeric(n)
  for (i in seq.int(2L, n)) {
    k <- seq_len(min(length(w), i - 1L))
    L[i] <- sum(incidence[i - k] * w[k])
  }
  L
}

fit_epifilter <- function(simulations, si_mean = mean_si, si_sd = sd_si) {
  source(file.path(epifilter_dir, "epiFilter.R"), local = TRUE)
  source(file.path(epifilter_dir, "epiSmoother.R"), local = TRUE)

  si <- make_si(si_mean, si_sd)
  Rgrid <- seq(epifilter_R_min, epifilter_R_max, length.out = epifilter_m)
  pR0 <- rep(1 / epifilter_m, epifilter_m)
  time <- nrow(simulations)

  fit_one <- function(i) {
    tryCatch({
      incidence <- simulations[, i]
      L <- total_infectiousness(incidence, si)
      td <- seq.int(epifilter_start, time)
      nd <- length(td)

      Rf <- epiFilter(Rgrid, epifilter_m, epifilter_eta, pR0, nd,
                      L[td], incidence[td], epifilter_alpha)
      Rs <- epiSmoother(Rgrid, epifilter_m, Rf[[4L]], Rf[[5L]], nd,
                        Rf[[6L]], epifilter_alpha)

      # Rhat rows 1 and 2 are the a and 1 - a quantiles
      out <- list(day = td, R = Rs[[1L]], lower = Rs[[2L]][1L, ],
                  upper = Rs[[2L]][2L, ], error = NA_character_)

      # NA inside the recursion window means the Poisson
      # likelihood underflowed on every grid point
      if (anyNA(c(out$R, out$lower, out$upper))) {
        return(list(error = "NaN in the recursion"))
      }
      out
    }, error = function(e) list(error = conditionMessage(e)))
  }

  secs <- system.time(
    fits <- run_forked(ncol(simulations), fit_one)
  )[["elapsed"]]

  collect_fits(fits, time, ncol(simulations), secs)
}


# ============================================================
# 3. EpiLPS (MAP) - Laplace approximation, estimR()
#
# Fitted to days epilps_start..T only (burn-in; see config.R).
# RLPS$Time counts from 1 within that window, so it is mapped
# back to the original day number.
# ============================================================

fit_epilps_map <- function(simulations, si_mean = mean_si, si_sd = sd_si) {
  si <- make_si(si_mean, si_sd)
  fit_days <- seq.int(epilps_start, nrow(simulations))  # burn-in, see config.R

  fit_one <- function(i) {
    tryCatch({
      f <- EpiLPS::estimR(incidence = simulations[fit_days, i],
                          si = si[-1L], K = epilps_K)
      list(day = fit_days[as.integer(f$RLPS$Time)],
           R = as.numeric(f$RLPS$Rq0.50),
           lower = as.numeric(f$RLPS$Rq0.025),
           upper = as.numeric(f$RLPS$Rq0.975),
           rho = as.numeric(f$NegBinoverdisp),
           error = NA_character_)
    }, error = function(e) list(error = conditionMessage(e)))
  }

  secs <- system.time(
    fits <- run_forked(ncol(simulations), fit_one)
  )[["elapsed"]]

  out <- collect_fits(fits, nrow(simulations), ncol(simulations), secs)
  out$rho <- vapply(fits, function(f) if (is.null(f$rho)) NA_real_ else f$rho,
                    numeric(1))
  out
}


# ============================================================
# 4. EpiLPS (MALA) - Langevin sampler, estimRmcmc()
#
# Seeded per replicate inside the fork, so the draws do not
# depend on how mclapply schedules the jobs.
# ============================================================

fit_epilps_mala <- function(simulations, si_mean = mean_si, si_sd = sd_si) {
  si <- make_si(si_mean, si_sd)
  fit_days <- seq.int(epilps_start, nrow(simulations))  # burn-in, see config.R

  fit_one <- function(i) {
    tryCatch({
      set.seed(rep_seed(i))
      f <- EpiLPS::estimRmcmc(incidence = simulations[fit_days, i],
                              si = si[-1L], K = epilps_K, niter = mcmc_niter,
                              burnin = mcmc_burnin, progressbar = FALSE)
      list(day = fit_days[as.integer(f$RLPS$Time)],
           R = as.numeric(f$RLPS$Rq0.50),
           lower = as.numeric(f$RLPS$Rq0.025),
           upper = as.numeric(f$RLPS$Rq0.975),
           rho = as.numeric(f$NegBinoverdisp),
           error = NA_character_)
    }, error = function(e) list(error = conditionMessage(e)))
  }

  secs <- system.time(
    fits <- run_forked(ncol(simulations), fit_one)
  )[["elapsed"]]

  out <- collect_fits(fits, nrow(simulations), ncol(simulations), secs)
  out$rho <- vapply(fits, function(f) if (is.null(f$rho)) NA_real_ else f$rho,
                    numeric(1))
  out
}
