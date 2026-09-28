# ==============================================================================
# serial_interval.R
# Discrete serial-interval constructors shared by every estimator and study.
# Layout convention throughout the project: si[1] is the lag-0 cell and is
# zero, si[k + 1] = P(SI = k). EpiEstim receives the whole vector; EpiLPS and
# EpiFilter receive si[-1].
# ==============================================================================

# Construct one discrete SI distribution shared by EpiEstim and EpiLPS.
# EpiEstim requires lag 0 as the first element; EpiLPS receives lag 1 onward.
make_si <- function(mean_si, sd_si, max_lag = 30L) {
  if (!is.numeric(mean_si) || length(mean_si) != 1L || mean_si <= 1) {
    stop("mean_si must be one numeric value greater than 1.")
  }
  if (!is.numeric(sd_si) || length(sd_si) != 1L || sd_si <= 0) {
    stop("sd_si must be one positive numeric value.")
  }

  out <- EpiEstim::discr_si(seq.int(0L, max_lag), mean_si, sd_si)
  out[1L] <- 0

  if (any(!is.finite(out)) || any(out < 0) || sum(out) <= 0) {
    stop("The constructed SI distribution is invalid.")
  }

  # Normalize after truncation at max_lag.
  out / sum(out)
}

# ------------------------------------------------------------------------------
# SI family misspecification
# ------------------------------------------------------------------------------
# make_si() above is Gamma-only: EpiEstim::discr_si() fits a gamma to mean - 1
# and shifts it a day, which is Cori's scheme and cannot be pointed at another
# family. The constructors below cover the families the misspecification studies
# fit under, so a study can vary the SHAPE of the assumed SI and not only its
# moments.
#
# The five levels are:
#
#   gamma_discr - make_si(), i.e. discr_si. The data-generating SI.
#   gamma_bin   - gamma under the +/-0.5 binning below.
#   lnorm       - lognormal, same binning
#   weibull     - Weibull, same binning
#   unif        - uniform, same binning
#
# gamma_bin exists because gamma_discr and the other three do NOT discretize the
# same way, and the gap is not negligible: TV(gamma_discr, gamma_bin) = 0.0188,
# the same 0.019 compare_si_discretization.R measures, and the lag-1 cell differs
# 3.2-fold. Comparing lnorm against gamma_discr would therefore mix the family
# effect with a discretization effect of comparable size. gamma_bin is the
# reference the other families are read against; gamma_discr - gamma_bin is the
# discretization cost, reported separately.
si_families <- c("gamma_discr", "gamma_bin", "lnorm", "weibull", "unif")

# Parameters that give a continuous family the requested mean and sd, so the
# only thing differing across families at one (mean, sd) is the density shape.
si_family_params <- function(family, mean_si, sd_si) {
  switch(
    family,
    gamma = list(shape = mean_si^2 / sd_si^2, rate = mean_si / sd_si^2),
    lnorm = {
      sdlog <- sqrt(log(1 + sd_si^2 / mean_si^2))
      list(meanlog = log(mean_si) - sdlog^2 / 2, sdlog = sdlog)
    },
    weibull = {
      # The coefficient of variation depends on the shape alone, so solve
      # CV(k) = sd/mean for k and read the scale off the mean. No closed form.
      # The bracket spans CV 0.02 to 20, far wider than any grid used here.
      cv <- sd_si / mean_si
      shape <- stats::uniroot(
        function(k) sqrt(gamma(1 + 2 / k) / gamma(1 + 1 / k)^2 - 1) - cv,
        interval = c(0.05, 200)
      )$root
      list(shape = shape, scale = mean_si / gamma(1 + 1 / shape))
    },
    unif = list(min = mean_si - sd_si * sqrt(3), max = mean_si + sd_si * sqrt(3)),
    stop("Unknown SI family: ", family)
  )
}

# One discrete SI from any of the families above, in the same layout make_si()
# returns: si[1] is the lag-0 cell and is zero, si[k + 1] = P(SI = k). That makes
# it a drop-in wherever make_si() is used.
#
# Non-gamma_discr families are binned at +/-0.5, i.e. w_k = F(k + 0.5) -
# F(k - 0.5) over lags 1..max_lag, then renormalized. This is EpiLPS::Idist()'s
# scheme and agrees with it to total variation < 1e-5 once the two are compared
# at the same support length.
#
# Mass below lag 0.5 and beyond max_lag is renormalized away, and at the extremes
# of the misspecification grids that is not a rounding detail: at a requested
# mean of 1.5 the realised mean comes back as 3.9 (binned gamma and uniform),
# 3.2 (Weibull) or 2.4 (lognormal), because the left tail the request implies
# falls off the bottom of the support. gamma_discr survives it, shifting a fitted
# gamma rather than binning one. So callers must check si_moments() on what they
# get rather than assume the request was honoured - the study scripts record
# realised_mean / realised_sd and flag the settings that drifted.
make_si_family <- function(family, mean_si, sd_si, max_lag = 30L) {
  if (!is.character(family) || length(family) != 1L) {
    stop("family must be one character value.")
  }
  if (identical(family, "gamma_discr")) {
    return(make_si(mean_si, sd_si, max_lag))
  }
  if (!is.numeric(mean_si) || length(mean_si) != 1L || mean_si <= 1) {
    stop("mean_si must be one numeric value greater than 1.")
  }
  if (!is.numeric(sd_si) || length(sd_si) != 1L || sd_si <= 0) {
    stop("sd_si must be one positive numeric value.")
  }

  base <- sub("_bin$", "", family)
  cdf <- switch(
    base,
    gamma   = stats::pgamma,
    lnorm   = stats::plnorm,
    weibull = stats::pweibull,
    unif    = stats::punif,
    stop("Unknown SI family: ", family)
  )

  params <- si_family_params(base, mean_si, sd_si)
  lags <- seq_len(max_lag)
  weights <- pmax(
    do.call(cdf, c(list(lags + 0.5), params)) -
      do.call(cdf, c(list(lags - 0.5), params)),
    0
  )

  out <- c(0, weights)
  if (any(!is.finite(out)) || any(out < 0) || sum(out) <= 0) {
    stop("The constructed SI distribution is invalid for family '", family,
         "' at mean ", mean_si, ", sd ", sd_si, ".")
  }

  out / sum(out)
}

# Realised moments of a discrete SI, where si[k + 1] = P(SI = k).
si_moments <- function(si) {
  lags <- seq.int(0L, length(si) - 1L)
  m <- sum(lags * si)
  c(mean = m, sd = sqrt(sum((lags - m)^2 * si)))
}

# The off-by-one guard every simulation study runs on its data-generating SI.
# Pairing I[t-k] with si[k] instead of si[k+1] shifts the whole distribution a
# day later and moves the realised mean from 7.5 to 8.5 - a mistake that is
# invisible in the output but changes the experiment. Checking the realised
# mean against the requested one catches it immediately. Returns the moments.
check_si_moments <- function(si, mean_si, tol = 0.05) {
  m <- si_moments(si)
  if (abs(m[["mean"]] - mean_si) > tol) {
    stop("Realised true SI mean is ", round(m[["mean"]], 3),
         ", expected ", mean_si, ". The lag indexing is off by one.")
  }
  invisible(m)
}
