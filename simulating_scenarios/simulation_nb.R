# ============================================================
# NEGATIVE-BINOMIAL RENEWAL SIMULATOR (EpiLPS arms only)
#
# EpiLPS fits a negative-binomial observation model and has no
# Poisson option, so its arms are given negative-binomial data:
#
#   I_t ~ NegBin(mean = R_t * Lambda_t, size = epilps_nb_rho)
#
# i.e. variance mu + mu^2 / rho. Everything else is a line-for-
# line copy of simulate_epidemic() / simulate_many() in the seed
# study's simulation_functions.R - the same seed day, the same
# SI indexing (w[1] = lag 0) - so the ONLY difference from the
# Poisson data EpiEstim and EpiFilter get is the draw.
# ============================================================

simulate_epidemic_nb <- function(R_t, w, time, I0, rho) {

  if (length(R_t) != time) {
    stop("Length of R_t must equal time.")
  }

  I <- numeric(time)
  I[1] <- I0

  for (t in 2:time) {
    max_s <- min(t, length(w))
    lambda_t <- sum(I[t - seq_len(max_s) + 1] * w[seq_len(max_s)])
    I[t] <- rnbinom(n = 1, mu = R_t[t] * lambda_t, size = rho)
  }

  I
}

simulate_many_nb <- function(R_t, w, time, I0, n_sim, rho) {

  simulations <- matrix(NA_real_, nrow = time, ncol = n_sim)

  for (i in seq_len(n_sim)) {
    simulations[, i] <- simulate_epidemic_nb(R_t, w, time, I0, rho)
  }

  simulations
}
