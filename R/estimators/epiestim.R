# ==============================================================================
# epiestim.R
# EpiEstim-specific fitting helpers. Both run non_parametric_si with an
# explicit si_distr, so the SI EpiEstim uses is exactly the vector the caller
# built and shared with the other estimators.
# ==============================================================================

# Fit EpiEstim once over all requested sliding windows.
fit_epiestim_full <- function(incidence_df, si_full, window_length) {
  n <- nrow(incidence_df)
  first_end <- window_length + 1L  # makes the first t_start equal to 2

  if (n < first_end) {
    stop("Not enough observations for the selected EpiEstim window.")
  }

  t_end <- seq.int(first_end, n)
  t_start <- t_end - window_length + 1L

  fit <- suppressWarnings(
    EpiEstim::estimate_R(
      incid = incidence_df,
      method = "non_parametric_si",
      config = EpiEstim::make_config(list(
        t_start = t_start,
        t_end = t_end,
        si_distr = si_full
      ))
    )
  )

  data.frame(
    index = t_end,
    date = incidence_df$dates[t_end],
    method = paste0("EpiEstim (", window_length, "-day window)"),
    R = fit$R$`Median(R)`,
    R_mean = fit$R$`Mean(R)`,
    R_median = fit$R$`Median(R)`,
    R_sd = fit$R$`Std(R)`,
    lower = fit$R$`Quantile.0.025(R)`,
    upper = fit$R$`Quantile.0.975(R)`,
    stringsAsFactors = FALSE
  )
}

# Fit ONE window, ending on the last row of incidence_df. This is the forecast
# origin fit: at origin t the backtest hands in the series through t and wants
# the posterior for the window [t - w + 1, t] only. An error becomes a row of
# NAs with the message in `status`, so a failed origin is recorded rather than
# aborting the whole backtest.
fit_epiestim_last <- function(incidence_df, si_full, window_length) {
  t <- nrow(incidence_df)
  t_start <- t - window_length + 1L

  fit <- tryCatch(
    suppressWarnings(
      EpiEstim::estimate_R(
        incid = incidence_df,
        method = "non_parametric_si",
        config = EpiEstim::make_config(list(
          t_start = t_start,
          t_end = t,
          si_distr = si_full
        ))
      )
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(data.frame(
      R_mean = NA_real_,
      R_sd = NA_real_,
      q025 = NA_real_,
      q50 = NA_real_,
      q975 = NA_real_,
      status = paste0("ERROR: ", conditionMessage(fit)),
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    R_mean = fit$R$`Mean(R)`[1L],
    R_sd = fit$R$`Std(R)`[1L],
    q025 = fit$R$`Quantile.0.025(R)`[1L],
    q50 = fit$R$`Median(R)`[1L],
    q975 = fit$R$`Quantile.0.975(R)`[1L],
    status = "OK",
    stringsAsFactors = FALSE
  )
}

# The harness adapter: what makes EpiEstim pluggable into fit_si_grid() and
# run_parallel() alongside fit_epilps_map() and run_epifilter(). It satisfies
# the shared estimator contract fit_fun(incidence, si, ...) -> one arm per
# window, each in the common layout (R, lower, upper, time, error), plus a
# top-level error field. fit_epiestim_full() cannot be handed to the harness as
# it stands: it returns a data frame keyed on `index` (the scorer aligns on
# `time`), carries no `error` field, and throws rather than reporting failure,
# which on a worker would take the whole replicate's batch of fits down with it.
#
# `incidence` is a bare count vector, as the simulator returns it. EpiEstim
# insists on a date column, so an arbitrary fixed origin supplies one; nothing
# downstream reads the dates, every join is on `time`. There is no `days`
# argument: run_epifilter() needs one because its recursion cannot start at day
# 1, but EpiEstim returns every window and score_interval_fits() picks the
# scored days out by match(score_days, f$time).
#
# R is the posterior median, as in fit_epiestim_full()'s R column and as the
# original study scored it. A failure comes back as list(error = message) so the
# caller can drop the replicate, as run_epifilter() does.
fit_epiestim_grid <- function(incidence, si, windows = c(1L, 7L),
                              origin = as.Date("2020-01-01")) {
  tryCatch({
    incidence_df <- data.frame(
      dates = origin + seq_along(incidence) - 1L,
      I = as.numeric(incidence)
    )
    arms <- lapply(windows, function(w) {
      fit <- fit_epiestim_full(incidence_df, si, w)
      list(
        R = fit$R_median,
        lower = fit$lower,
        upper = fit$upper,
        time = as.integer(fit$index),
        error = NA_character_
      )
    })
    names(arms) <- paste0("w", windows)
    c(arms, list(error = NA_character_))
  }, error = function(e) list(error = conditionMessage(e)))
}
