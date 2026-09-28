# ==============================================================================
# sensitivity.R
# Shared driver for the SI sensitivity analyses. The original script ran the
# same code twice (Section 5 varying the assumed SI SD, Section 5B varying the
# assumed SI mean); both are now one function parameterized by `vary`.
# Requires utils.R, epiestim.R and epilps.R to be sourced first.
# ==============================================================================

# vary = "sd":   grid holds assumed SD values, mean is held at mean_si.
# vary = "mean": grid holds assumed mean values, SD is held at sd_si.
# Returns the combined, evaluation-window-filtered Rt sensitivity data frame,
# with the same columns (and label column) as the original per-section output.
run_si_sensitivity <- function(vary = c("sd", "mean"),
                               grid,
                               incidence_data,
                               mean_si,
                               sd_si,
                               max_si_lag,
                               epiestim_window,
                               K_epilps,
                               n_cores,
                               evaluation_start_date) {
  vary <- match.arg(vary)

  # Build the SI for one grid point, holding the other SI moment fixed.
  make_current_si <- function(value) {
    if (vary == "sd") {
      make_si(mean_si, value, max_si_lag)
    } else {
      make_si(value, sd_si, max_si_lag)
    }
  }

  # EpiEstim sensitivity fits.
  epiestim_sens <- map_dfr(
    grid,
    function(value) {
      fit_epiestim_full(
        incidence_df = incidence_data,
        si_full = make_current_si(value),
        window_length = epiestim_window
      ) %>%
        mutate(
          assumed_mean_si = if (vary == "mean") value else mean_si,
          assumed_sd_si = if (vary == "sd") value else sd_si
        )
    }
  )

  # EpiLPS sensitivity fits: one clean process per SI scenario.
  epilps_sens_jobs <- map(
    grid,
    function(value) {
      list(
        job_id = paste0(vary, "_", value),
        incidence = incidence_data$I,
        si_epilps = make_current_si(value)[-1L],
        K = K_epilps,
        dates = as.character(incidence_data$dates),
        return_type = "full"
      )
    }
  )

  epilps_sens_raw <- bind_rows(
    run_epilps_jobs(epilps_sens_jobs, cores = n_cores)
  )

  epilps_sens <- epilps_sens_raw %>%
    mutate(
      assumed_value = as.numeric(sub(paste0("^", vary, "_"), "", job_id)),
      assumed_mean_si = if (vary == "mean") assumed_value else mean_si,
      assumed_sd_si = if (vary == "sd") assumed_value else sd_si
    ) %>%
    transmute(
      index = index,
      date = as.Date(date),
      method = "EpiLPS",
      R = q50,
      R_mean = R_mean,
      R_sd = R_sd,
      lower = q025,
      upper = q975,
      assumed_mean_si = assumed_mean_si,
      assumed_sd_si = assumed_sd_si,
      converged = converged
    )

  out <- bind_rows(
    epiestim_sens %>%
      select(
        index, date, method, R, R_mean, R_sd, lower, upper,
        assumed_mean_si, assumed_sd_si
      ) %>%
      mutate(converged = NA_character_),
    epilps_sens
  ) %>%
    filter(date >= evaluation_start_date)

  # Keep the per-scenario label column the original sections produced.
  if (vary == "sd") {
    out <- out %>%
      mutate(assumed_sd_label = paste0("Assumed SI SD = ", assumed_sd_si)) %>%
      arrange(method, assumed_sd_si, date)
  } else {
    out <- out %>%
      mutate(assumed_mean_label = paste0("Assumed SI mean = ", assumed_mean_si)) %>%
      arrange(method, assumed_mean_si, date)
  }

  out
}



