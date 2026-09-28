# ==============================================================================
# epilps.R
# EpiLPS-specific fitting helpers, in two styles for two kinds of caller.
#
# 1. The JOB RUNNER (epilps_worker / run_epilps_jobs / tidy_epilps_result):
#    every fit in its own fresh R process, returning a data frame in the common
#    Rt layout. Used by main.R, sensitivity.R and sim_study.R, where fits are
#    few and each one must be isolated.
#
# 2. The IN-PROCESS FITTERS (fit_epilps_map / fit_epilps_mala): a plain call
#    returning a list, for the misspecification studies, which fork once per
#    replicate and fit 17+ SI settings inside that fork. See fit_si_grid() in
#    studies/si_misspec.R for why that is safe.
#
# Both hand EpiLPS si[-1]: the lag-1-onward vector.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Fresh-process job runner
# ------------------------------------------------------------------------------

# Worker used in a fresh child R process for every EpiLPS fit.
# job$return_type is either "full" or "last".
epilps_worker <- function(job) {
  fit <- EpiLPS::estimR(
    incidence = job$incidence,
    si = job$si_epilps,
    K = job$K,
    dates = job$dates
  )

  out <- fit$RLPS
  out$Time <- as.character(out$Time)

  if (identical(job$return_type, "last")) {
    out <- out[nrow(out), , drop = FALSE]
  }

  data.frame(
    job_id = job$job_id,
    index = if (identical(job$return_type, "last")) {
      length(job$incidence)
    } else {
      seq_along(job$incidence)
    },
    date = out$Time,
    R_mean = as.numeric(out$R),
    R_sd = as.numeric(out$Rsd),
    q025 = as.numeric(out$Rq0.025),
    q50 = as.numeric(out$Rq0.50),
    q975 = as.numeric(out$Rq0.975),
    overdispersion = fit$NegBinoverdisp,
    converged = as.character(fit$optimconverged),
    stringsAsFactors = FALSE
  )
}

safe_epilps_worker <- function(job) {
  tryCatch(
    epilps_worker(job),
    error = function(e) {
      data.frame(
        job_id = job$job_id,
        index = if (identical(job$return_type, "last")) {
          length(job$incidence)
        } else {
          seq_along(job$incidence)
        },
        date = if (identical(job$return_type, "last")) {
          tail(job$dates, 1L)
        } else {
          job$dates
        },
        R_mean = NA_real_,
        R_sd = NA_real_,
        q025 = NA_real_,
        q50 = NA_real_,
        q975 = NA_real_,
        overdispersion = NA_real_,
        converged = paste0("ERROR: ", conditionMessage(e)),
        stringsAsFactors = FALSE
      )
    }
  )
}

# Run each EpiLPS job in a fresh process.
# This avoids reusing the compiled backend within a single R process.
run_epilps_jobs <- function(jobs, cores = 1L) {
  if (length(jobs) == 0L) return(list())

  if (.Platform$OS.type != "windows") {
    if (length(jobs) == 1L) {
      child <- parallel::mcparallel(safe_epilps_worker(jobs[[1L]]))
      return(unname(parallel::mccollect(child)))
    }

    # mc.cores must be at least 2 here; otherwise mclapply evaluates jobs in
    # the current process and defeats the fresh-process safeguard.
    fork_cores <- max(2L, min(as.integer(cores), length(jobs)))

    return(
      parallel::mclapply(
        jobs,
        safe_epilps_worker,
        mc.cores = fork_cores,
        mc.preschedule = FALSE,
        mc.set.seed = TRUE
      )
    )
  }

  # Windows does not support forked processes. The function below is fully
  # self-contained so that callr can evaluate it in a clean R session.
  lapply(
    jobs,
    function(job) {
      callr::r(
        func = function(job) {
          tryCatch({
            fit <- EpiLPS::estimR(
              incidence = job$incidence,
              si = job$si_epilps,
              K = job$K,
              dates = job$dates
            )

            out <- fit$RLPS
            out$Time <- as.character(out$Time)
            if (identical(job$return_type, "last")) {
              out <- out[nrow(out), , drop = FALSE]
            }

            data.frame(
              job_id = job$job_id,
              index = if (identical(job$return_type, "last")) {
                length(job$incidence)
              } else {
                seq_along(job$incidence)
              },
              date = out$Time,
              R_mean = as.numeric(out$R),
              R_sd = as.numeric(out$Rsd),
              q025 = as.numeric(out$Rq0.025),
              q50 = as.numeric(out$Rq0.50),
              q975 = as.numeric(out$Rq0.975),
              overdispersion = fit$NegBinoverdisp,
              converged = as.character(fit$optimconverged),
              stringsAsFactors = FALSE
            )
          }, error = function(e) {
            data.frame(
              job_id = job$job_id,
              index = if (identical(job$return_type, "last")) {
                length(job$incidence)
              } else {
                seq_along(job$incidence)
              },
              date = if (identical(job$return_type, "last")) {
                tail(job$dates, 1L)
              } else {
                job$dates
              },
              R_mean = NA_real_, R_sd = NA_real_,
              q025 = NA_real_, q50 = NA_real_, q975 = NA_real_,
              overdispersion = NA_real_,
              converged = paste0("ERROR: ", conditionMessage(e)),
              stringsAsFactors = FALSE
            )
          })
        },
        args = list(job = job),
        spinner = FALSE,
        show = FALSE
      )
    }
  )
}

# Convert the raw result of a return_type = "full" EpiLPS job into the common
# Rt data frame layout used by the primary fit and both sensitivity analyses.
tidy_epilps_result <- function(raw) {
  raw %>%
    transmute(
      index = index,
      date = as.Date(date),
      method = "EpiLPS",
      R = q50,
      R_mean = R_mean,
      R_sd = R_sd,
      lower = q025,
      upper = q975,
      overdispersion = overdispersion,
      converged = converged
    )
}

# ------------------------------------------------------------------------------
# 2. In-process fitters for the misspecification studies
# ------------------------------------------------------------------------------
# Both return a list in the layout score_interval_fits() reads: R (the posterior
# median, matching the pipeline), lower / upper (equal-tailed 95%), time, rho
# (the NegBin overdispersion), and error - NA on success, the message otherwise.
# A failed fit is a list with R = NULL, so a caller can drop it rather than crash.

# LPSMAP: the Laplace approximation to the posterior of the B-spline
# coefficients, lognormal quantiles for Rt.
fit_epilps_map <- function(incidence, si, K) {
  tryCatch({
    f <- EpiLPS::estimR(incidence = incidence, si = si[-1L], K = K)
    list(
      R     = as.numeric(f$RLPS$Rq0.50),
      lower = as.numeric(f$RLPS$Rq0.025),
      upper = as.numeric(f$RLPS$Rq0.975),
      time  = as.integer(f$RLPS$Time),
      rho   = as.numeric(f$NegBinoverdisp),
      error = NA_character_
    )
  }, error = function(e) list(R = NULL, error = conditionMessage(e)))
}

# LPSMALA: a Metropolis-adjusted Langevin sampler over the same posterior,
# empirical quantiles for Rt. Returns the same fields plus the 95% HPD bounds,
# the convergence flag and the fit time, so equal-tailed and HPD intervals can
# be scored side by side. progressbar must be FALSE: in a forked worker the bar
# writes interleaved garbage to stdout. rho here is a posterior mean, not a mode.
fit_epilps_mala <- function(incidence, si, K, niter, burnin) {
  tryCatch({
    f <- EpiLPS::estimRmcmc(incidence = incidence, si = si[-1L], K = K,
                            niter = niter, burnin = burnin,
                            progressbar = FALSE)
    list(
      R         = as.numeric(f$RLPS$Rq0.50),
      lower     = as.numeric(f$RLPS$Rq0.025),
      upper     = as.numeric(f$RLPS$Rq0.975),
      hpd_lower = as.numeric(f$HPD95_Rt[[1L]]),
      hpd_upper = as.numeric(f$HPD95_Rt[[2L]]),
      time      = as.integer(f$RLPS$Time),
      rho       = as.numeric(f$NegBinoverdisp),
      converged = isTRUE(f$optimconverged),
      secs      = as.numeric(f$LPS_elapsed),
      error     = NA_character_
    )
  }, error = function(e) list(R = NULL, error = conditionMessage(e)))
}
