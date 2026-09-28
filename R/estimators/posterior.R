# ==============================================================================
# posterior.R
# Posterior Rt draws reconstructed from each estimator's reported summary, used
# by the one-step-ahead forecast backtest (main.R Section 7), plus the
# log-mean-exp the forecast log score needs.
# ==============================================================================

# Numerically stable log(mean(exp(log_values))).
log_mean_exp <- function(log_values) {
  log_values <- log_values[is.finite(log_values)]
  if (length(log_values) == 0L) return(-Inf)
  max_log <- max(log_values)
  max_log + log(mean(exp(log_values - max_log)))
}

# Reconstruct EpiEstim's Gamma posterior from posterior mean and SD.
sample_epiestim_R <- function(mean_R, sd_R, n_draws) {
  if (!is.finite(mean_R) || !is.finite(sd_R) || mean_R <= 0 || sd_R <= 0) {
    return(rep(NA_real_, n_draws))
  }

  shape <- (mean_R / sd_R)^2
  rate  <- mean_R / (sd_R^2)
  stats::rgamma(n_draws, shape = shape, rate = rate)
}

# EpiLPS estimR() reports lognormal-based quantiles for its MAP/Laplace fit.
# Recover the lognormal parameters from the reported median and 95% quantiles.
sample_epilps_R <- function(q025, q50, q975, n_draws) {
  values <- c(q025, q50, q975)
  if (any(!is.finite(values)) || any(values <= 0) || q025 >= q975) {
    return(rep(NA_real_, n_draws))
  }

  z975 <- stats::qnorm(0.975)
  meanlog <- log(q50)
  sdlog <- (log(q975) - log(q025)) / (2 * z975)

  if (!is.finite(sdlog) || sdlog <= 0) {
    return(rep(NA_real_, n_draws))
  }

  stats::rlnorm(n_draws, meanlog = meanlog, sdlog = sdlog)
}
