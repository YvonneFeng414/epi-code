# ==============================================================================
# epifilter.R
# EpiFilter (Parag 2020): Bayesian recursive filtering and smoothing of Rt on a
# grid, with a Poisson renewal observation model and a diffusion state model.
#
# Section 1 is VENDORED. EpiFilter is not on CRAN; epiFilter() and
# epiSmoother() below are taken verbatim from https://github.com/kpzoo/EpiFilter
# so the recursions are the author's own. Do not edit them.
#
# Section 2 is this project's wrapper: the grid, the total-infectiousness
# series the SI enters through, and one forward+backward pass returning both
# arms in the common R / lower / upper / time layout the scorer reads.
# ==============================================================================

# ==============================================================================
# 1. Vendored from kpzoo/EpiFilter - DO NOT EDIT
# ==============================================================================

######################################################################
## Bayesian recursive filtering via EpiFilter
# From: Parag, KV, (2020) “Improved real-time estimation of reproduction numbers
# at low case incidence and between epidemic waves” BioRxiv.
######################################################################

# Assumptions
# - observation model is Poisson renewal equation (as in EpiEstim)
# - reproduction number state space model is a simple diffusion

# Inputs - grid on reproduction numbers (Rgrid), size of grid (m), diffusion noise (eta),
# prior on R (pR0), max time (nday), total infectiousness (Lday), incidence (Iday), confidence (a)

# Output - mean (Rmean), median (Rmed), 50% and 95% quantiles of estimates (Rhat),
# causal posterior over R (pR), pre-update (pRup) and state transition matrix (pstate)

epiFilter <- function(Rgrid, m, eta, pR0, nday, Lday, Iday, a){

  # Probability vector for R and prior
  pR = matrix(0, nday, m); pRup = pR
  pR[1, ] = pR0; pRup[1, ] = pR0

  # Mean and median estimates
  Rmean = rep(0, nday); Rmed = Rmean
  # 50% and 95% (depends on a) confidence on R
  Rhat = matrix(0, 4, nday)

  # Initialise mean
  Rmean[1] = pR[1, ]%*%Rgrid
  # CDF of prior
  Rcdf0 = cumsum(pR0)
  # Initialise quartiles
  idm = which(Rcdf0 >= 0.5, 1); Rmed[1] = Rgrid[idm[1]]
  id1 = which(Rcdf0 >= a, 1); id2 = which(Rcdf0 >= 1-a, 1)
  id3 = which(Rcdf0 >= 0.25, 1); id4 = which(Rcdf0 >= 0.75, 1)
  Rhat[1, 1] = Rgrid[id1[1]]; Rhat[2, 1] = Rgrid[id2[1]]
  Rhat[3, 1] = Rgrid[id3[1]]; Rhat[4, 1] = Rgrid[id4[1]]

  # Precompute state distributions for R transitions
  pstate = matrix(0, m, m);
  for(j in 1:m){
    pstate[j, ] = dnorm(Rgrid[j], Rgrid, sqrt(Rgrid)*eta)
  }

  # Update prior to posterior sequentially
  for(i in 2:nday){
    # Compute mean from Poisson renewal (observation model)
    rate = Lday[i]*Rgrid
    # Probabilities of observations
    pI = dpois(Iday[i], rate)

    # State predictions for R
    pRup[i, ]  = pR[i-1, ]%*%pstate
    # Update to posterior over R
    pR[i, ] = pRup[i, ]*pI
    pR[i, ] = pR[i, ]/sum(pR[i, ])

    # Posterior mean and CDF
    Rmean[i] = pR[i, ]%*%Rgrid
    Rcdf = cumsum(pR[i, ])

    # Quantiles for estimates
    idm = which(Rcdf >= 0.5, 1); Rmed[i] = Rgrid[idm[1]]
    id1 = which(Rcdf >= a, 1); id2 = which(Rcdf >= 1-a, 1)
    id3 = which(Rcdf >= 0.25, 1); id4 = which(Rcdf >= 0.75, 1)
    Rhat[1, i] = Rgrid[id1[1]]; Rhat[2, i] = Rgrid[id2[1]]
    Rhat[3, i] = Rgrid[id3[1]]; Rhat[4, i] = Rgrid[id4[1]]
  }

  # Main outputs: estimates of R and states
  epiFilter = list(Rmed, Rhat, Rmean, pR, pRup, pstate)
}

######################################################################
## Bayesian recursive smoothing via EpiFilter
# From: Parag, KV, (2020) “Improved real-time estimation of reproduction numbers
# at low case incidence and between epidemic waves” BioRxiv.
######################################################################

# Assumptions
# - observation model is Poisson renewal equation (as in EpiEstim)
# - reproduction number state space model is a simple diffusion
# - must have run epiFilter first to obtain forward distribution pR
# - method makes a backward pass to generate qR

# Inputs - grid on reproduction numbers (Rgrid), size of grid (m), filtered posterior (pR),
# update pre-filter (pRup), max time (nday), state transition matrix (pstate), confidence (a)

# Output - mean (Rmean), median (Rmed), lower (Rlow) amd upper (Rhigh) quantiles of estimates,
# smoothed posterior over R (qR) which is backwards and forwards

epiSmoother <- function(Rgrid, m, pR, pRup, nday, pstate, a){

  # Last smoothed distribution same as filtered
  qR = matrix(0, nday, m); qR[nday, ] = pR[nday, ]

  # Main smoothing equation iteratively computed
  for(i in seq(nday-1, 1)){
    # Remove zeros
    pRup[i+1, pRup[i+1, ] == 0] = 10^-8

    # Integral term in smoother
    integ = qR[i+1, ]/pRup[i+1, ]
    integ = integ%*%pstate

    # Smoothed posterior over Rgrid
    qR[i, ] = pR[i, ]*integ
    # Force a normalisation
    qR[i, ] = qR[i, ]/sum(qR[i, ]);
  }

  # Mean, median estimats of R
  Rmean = rep(0, nday); Rmed = Rmean
  # 50% and 95% (depends on a) confidence on R
  Rhat = matrix(0, 4, nday)

  # Compute at every time point
  for (i in 1:nday) {
    # Posterior mean and CDF
    Rmean[i] = qR[i, ]%*%Rgrid
    Rcdf = cumsum(qR[i, ])

    # Quantiles for estimates
    idm = which(Rcdf >= 0.5); Rmed[i] = Rgrid[idm[1]]
    id1 = which(Rcdf >= a, 1); id2 = which(Rcdf >= 1-a, 1)
    id3 = which(Rcdf >= 0.25, 1); id4 = which(Rcdf >= 0.75, 1)
    Rhat[1, i] = Rgrid[id1[1]]; Rhat[2, i] = Rgrid[id2[1]]
    Rhat[3, i] = Rgrid[id3[1]]; Rhat[4, i] = Rgrid[id4[1]]
  }

  # Main outputs: estimates of R and states
  epiSmoother = list(Rmed, Rhat, Rmean, qR)
}

# ==============================================================================
# 2. Project wrapper
# ==============================================================================

epifilter_arms <- c("smoother", "filter")

# The R grid and the uniform prior over it, following vignetteCOVID.R. `alpha`
# is the tail probability, so alpha = 0.025 makes Rhat rows 1 and 2 the 2.5%
# and 97.5% quantiles.
make_epifilter_grid <- function(R_min, R_max, m, alpha) {
  m <- as.integer(m)
  list(
    Rgrid = seq(R_min, R_max, length.out = m),
    m = m,
    pR0 = rep(1 / m, m),
    alpha = alpha
  )
}

# The SI enters EpiFilter only here: Lday[i] = sum_k I[i-k] w[k], following
# vignetteCOVID.R's wdist = dgamma(1:nday, ...) and
# Lday[i] = sum(Iday[(i-1):1] * wdist[1:(i-1)]). EpiFilter wants the
# lag-1-onward vector, so w = si[-1].
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

# One forward+backward pass, returning both arms on `days`:
#   filter   - epiFilter, causal / real-time
#   smoother - epiSmoother, retrospective, uses the whole series
#
# EpiFilter cannot start at day 1 (see epifilter_si_misspec.R on filter_start),
# so the recursion runs over filter_start..length(incidence) and `days` are
# matched into that window. Each arm is a list in the common layout
# (R, lower, upper, time, error); a NaN anywhere in the recursion, or an error,
# comes back as list(error = message) so the caller can drop the replicate.
run_epifilter <- function(incidence, si, grid, eta, filter_start, days,
                          arms = epifilter_arms) {
  tryCatch({
    L <- total_infectiousness(incidence, si)
    td <- seq.int(filter_start, length(incidence))
    nd <- length(td)

    Rf <- epiFilter(grid$Rgrid, grid$m, eta, grid$pR0, nd, L[td], incidence[td], grid$alpha)
    Rs <- epiSmoother(grid$Rgrid, grid$m, Rf[[4L]], Rf[[5L]], nd, Rf[[6L]], grid$alpha)

    idx <- match(days, td)
    pick <- function(o) {
      list(R = o[[1L]][idx], lower = o[[2L]][1L, idx], upper = o[[2L]][2L, idx],
           time = as.integer(days), error = NA_character_)
    }
    out <- list(filter = pick(Rf), smoother = pick(Rs), days = days,
                error = NA_character_)

    numbers <- unlist(lapply(out[arms], function(a) c(a$R, a$lower, a$upper)))
    if (anyNA(numbers)) {
      return(list(error = "NaN in the recursion"))
    }
    out
  }, error = function(e) list(error = conditionMessage(e)))
}
