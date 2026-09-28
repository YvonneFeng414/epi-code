# ==============================================================================
# renewal.R
# The Poisson / negative-binomial renewal equation, in its two uses: simulating
# an epidemic forward from a known Rt path, and computing the infectiousness
# that a one-step-ahead forecast multiplies by Rt.
#
# si_full[1] is lag 0 and is zero; si_full[lag + 1] is P(SI = lag). Element
# k + 1 of si is P(SI = k), so lag k must index k + 1 - see check_si_moments().
# ==============================================================================

# Simulate one replicate epidemic from a known Rt trajectory.
#
# The first length(seed_incidence) days are copied verbatim (from the observed
# series, or a single I0) and are never scored. Every later day is a draw from
# the observation model: Poisson, which is EpiEstim's and EpiFilter's own, or
# negative binomial with the given overdispersion (size), which is what EpiLPS
# estimates on the Queens series and what the EpiLPS misspecification studies
# simulate.
#
# A replicate whose recursion diverges comes back as NULL rather than as NAs:
# rpois() returns NA for a rate beyond the integer range, and a silent NA would
# poison every downstream metric.
simulate_renewal_incidence <- function(seed_incidence, R_true, si_full,
                                       n_days = length(R_true),
                                       obs_model = c("poisson", "nbinom"),
                                       overdispersion = NULL) {
  obs_model <- match.arg(obs_model)
  if (obs_model == "nbinom" &&
      (!is.numeric(overdispersion) || length(overdispersion) != 1L ||
       !is.finite(overdispersion) || overdispersion <= 0)) {
    stop("obs_model = 'nbinom' needs one positive finite overdispersion (size).")
  }

  max_lag <- length(si_full) - 1L
  n_seed <- length(seed_incidence)

  if (n_seed < 1L || n_seed >= n_days) {
    stop("seed_incidence must be shorter than the requested series length.")
  }
  if (length(R_true) != n_days) {
    stop("R_true must have one value per simulated day.")
  }

  incidence <- numeric(n_days)
  incidence[seq_len(n_seed)] <- seed_incidence

  for (t in seq.int(n_seed + 1L, n_days)) {
    lags <- seq_len(min(max_lag, t - 1L))
    lambda <- R_true[t] * sum(incidence[t - lags] * si_full[lags + 1L])

    if (!is.finite(lambda) || lambda > 1e9) return(NULL)

    incidence[t] <- if (obs_model == "poisson") {
      stats::rpois(1L, lambda = max(lambda, 0))
    } else {
      stats::rnbinom(1L, mu = max(lambda, 0), size = overdispersion)
    }
  }

  if (anyNA(incidence)) return(NULL)
  as.integer(incidence)
}

# The same recursion with the observation draw replaced by its expectation:
# the noise-free renewal equation. Used by the Britton validation as the
# cross-check that the renewal transient converges onto the pure exponential
# the Euler-Lotka equation describes. A scalar R_true is held constant. Kept
# apart from simulate_renewal_incidence() so that function stays bit-identical
# to its archived oracle in tests/test_shared_helpers.R.
renewal_expectation <- function(seed_incidence, R_true, si_full,
                                n_days = if (length(R_true) > 1L) length(R_true) else NULL) {
  if (is.null(n_days)) stop("n_days is needed when R_true is a single value.")
  if (length(R_true) == 1L) R_true <- rep(R_true, n_days)

  max_lag <- length(si_full) - 1L
  n_seed <- length(seed_incidence)
  if (n_seed < 1L || n_seed >= n_days) {
    stop("seed_incidence must be shorter than the requested series length.")
  }
  if (length(R_true) != n_days) {
    stop("R_true must have one value per simulated day.")
  }

  incidence <- numeric(n_days)
  incidence[seq_len(n_seed)] <- seed_incidence
  for (t in seq.int(n_seed + 1L, n_days)) {
    lags <- seq_len(min(max_lag, t - 1L))
    incidence[t] <- R_true[t] * sum(incidence[t - lags] * si_full[lags + 1L])
  }
  incidence
}

# Calculate infectiousness Lambda_(t+1) from observed incidence through day t.
calc_infectiousness <- function(incidence, t, si_full) {
  max_lag <- length(si_full) - 1L
  if (t < max_lag) {
    stop("t must be at least the maximum SI lag.")
  }

  lags <- seq_len(max_lag)
  past_indices <- t + 1L - lags
  sum(incidence[past_indices] * si_full[lags + 1L])
}
