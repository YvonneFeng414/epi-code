# ==============================================================================
# sim_study.R
# Driver for the fully synthetic simulation study (sim_main.R).
#
# Split into two layers because the epidemics depend ONLY on the Rt scenario:
# they are always generated with the true SI, so one set of replicates per Rt
# scenario is refitted under every assumed SI. That is what keeps the simulation
# cost at length(rt_scenarios) * n_sim instead of multiplying by the SI grid.
#
# Requires utils.R, simulation.R, epiestim.R, epilps.R, metrics.R and
# scenarios.R to be sourced first.
# ==============================================================================

# ------------------------------------------------------------------------------
# Layer 1: simulate replicates for one Rt scenario
# ------------------------------------------------------------------------------
# Always uses si_true. A replicate whose renewal recursion diverges comes back as
# NULL from simulate_renewal_incidence() and is dropped here rather than being
# allowed to become NAs downstream; the count is returned so the caller can
# report it, because a scenario that quietly loses most of its replicates would
# otherwise look like a successful run.
simulate_scenario_replicates <- function(scenario, time, n_sim, I0, si_true,
                                        seed) {
  R_true <- build_rt_scenario(scenario, time)

  set.seed(seed)
  series <- map(
    seq_len(n_sim),
    function(r) {
      simulate_renewal_incidence(
        seed_incidence = I0,
        R_true = R_true,
        si_full = si_true,
        n_days = time
      )
    }
  )

  failed <- vapply(series, is.null, logical(1))
  if (all(failed)) {
    stop("Every replicate diverged for scenario '", scenario$label,
         "'; lower the Rt path, sim_I0, or sim_time.")
  }

  list(
    R_true = R_true,
    series = series[!failed],
    replicate_ids = which(!failed),
    n_requested = n_sim,
    n_failed = sum(failed)
  )
}

# ------------------------------------------------------------------------------
# Per-arm truth table
# ------------------------------------------------------------------------------
# The one piece of genuinely new logic relative to main.R Section 6, which only
# ever handles a single window and so can keep two fixed truth columns.
#
# Here several windows run side by side, and a w-day sliding window estimates the
# AVERAGE Rt over [t - w + 1, t]. Each arm's own estimand is therefore the
# trailing mean over ITS OWN window, which has to be joined per arm:
#
#   true_R_instant - the instantaneous truth; the estimand of EpiLPS and of the
#                    1-day EpiEstim arm
#   true_R_window  - the trailing mean over reference_window, held fixed across
#                    arms so that column is comparable down the table
#   true_R_own     - the arm's own estimand; the only one that tests calibration
#
# Scoring against a target that is not an arm's estimand produces a deterministic
# offset, not evidence of miscalibration. epi.R scored every window against the
# instantaneous truth, which is correct only for a constant Rt.
build_truth_lookup <- function(R_true, windows, reference_window) {
  instant <- as.numeric(R_true)
  reference <- trailing_mean(R_true, reference_window)

  map_dfr(
    sort(unique(as.integer(windows))),
    function(w) {
      data.frame(
        estimand_window = w,
        index = seq_along(instant),
        true_R_instant = instant,
        true_R_window = reference,
        true_R_own = if (w == 1L) instant else trailing_mean(R_true, w),
        stringsAsFactors = FALSE
      )
    }
  )
}

# The first index every arm can be scored at. Fixed in advance from the settings
# rather than read off the fitted values (epi.R used which(!is.na(...))[1], which
# makes the evaluation window depend on the data). A w-day EpiEstim arm's first
# estimate lands on index w + 1; trailing_mean() is defined from index w onward.
sim_score_start_index <- function(windows, reference_window) {
  max(max(as.integer(windows)) + 1L, as.integer(reference_window))
}

# ------------------------------------------------------------------------------
# Layer 2: refit and score one (Rt scenario x assumed SI) cell
# ------------------------------------------------------------------------------
# `replicates` is the output of simulate_scenario_replicates(); `si_assumed` is a
# full SI vector (lag 0 first) as returned by make_si().
#
# The two run_* toggles select which arms this cell contains. sim_main.R sets
# run_epilps = FALSE on the SI scenarios it cannot afford EpiLPS for; a study that
# pairs each estimator with its own plug-in truth (truth_source_main.R) needs the
# mirror image, one cell of each. At least one must be TRUE.
fit_scenario_grid <- function(replicates,
                              si_assumed,
                              dates,
                              windows,
                              reference_window,
                              K_epilps,
                              n_cores,
                              run_epiestim = TRUE,
                              run_epilps = TRUE) {
  if (!isTRUE(run_epiestim) && !isTRUE(run_epilps)) {
    stop("At least one of run_epiestim / run_epilps must be TRUE.")
  }

  windows <- sort(unique(as.integer(windows)))
  replicate_ids <- replicates$replicate_ids

  # --- EpiEstim arms: one fit per (window, replicate) -------------------------
  epiestim_draws <- NULL
  if (isTRUE(run_epiestim)) {
    epiestim_draws <- map_dfr(
      windows,
      function(window_length) {
        map_dfr(
          seq_along(replicate_ids),
          function(j) {
            replicate_incidence <- data.frame(
              dates = dates,
              I = replicates$series[[j]]
            )

            fit <- tryCatch(
              fit_epiestim_full(
                incidence_df = replicate_incidence,
                si_full = si_assumed,
                window_length = window_length
              ),
              error = function(e) NULL
            )

            if (is.null(fit)) return(NULL)

            fit %>%
              mutate(
                replicate = replicate_ids[j],
                estimand_window = window_length
              ) %>%
              select(replicate, index, method, estimand_window,
                     R, R_sd, lower, upper)
          }
        )
      }
    )
  }

  # --- EpiLPS arm: one fresh process per replicate ----------------------------
  # EpiLPS output does not depend on the window (it is not a windowed estimator),
  # so one fit per replicate serves the whole table. It DOES depend on the
  # assumed SI, which is why this cannot be hoisted out of the SI loop.
  epilps_draws <- NULL
  if (isTRUE(run_epilps)) {
    epilps_jobs <- map(
      seq_along(replicate_ids),
      function(j) {
        list(
          job_id = paste0("rep_", replicate_ids[j]),
          incidence = replicates$series[[j]],
          si_epilps = si_assumed[-1L],
          K = K_epilps,
          dates = as.character(dates),
          return_type = "full"
        )
      }
    )

    epilps_draws <- bind_rows(
      run_epilps_jobs(epilps_jobs, cores = n_cores)
    ) %>%
      transmute(
        replicate = as.integer(sub("^rep_", "", job_id)),
        index = index,
        method = "EpiLPS",
        # The LPS smoother targets the instantaneous Rt, so its estimand is the
        # trailing average over a window of one day.
        estimand_window = 1L,
        R = q50,
        R_sd = R_sd,
        lower = q025,
        upper = q975,
        converged = converged
      )

    failed_fits <- epilps_draws %>%
      filter(is.na(R)) %>%
      distinct(replicate, converged)

    if (nrow(failed_fits) > 0L) {
      cat("    EpiLPS failed on ", nrow(failed_fits),
          " replicate(s); first message: ", failed_fits$converged[1L], "\n",
          sep = "")
    }

    epilps_draws <- epilps_draws %>%
      select(replicate, index, method, estimand_window, R, R_sd, lower, upper)
  }

  # --- Join the per-arm truth and trim to the common evaluation window --------
  start_index <- sim_score_start_index(windows, reference_window)

  bind_rows(epiestim_draws, epilps_draws) %>%
    filter(index >= start_index) %>%
    left_join(
      build_truth_lookup(replicates$R_true, windows, reference_window),
      by = c("estimand_window", "index")
    ) %>%
    arrange(method, replicate, index)
}

# ------------------------------------------------------------------------------
# Day-level coverage across replicates
# ------------------------------------------------------------------------------
# Same shape as main.R Section 6's rt_coverage_by_day: the fraction of replicates
# whose interval covered each target on each day. The overall rate hides where
# coverage fails, which for these scenarios is exactly where the truth moves
# fastest (the step at day 50, the steepest part of the sine).
sim_coverage_by_day <- function(draws, reference_window) {
  window_target_label <- paste0(reference_window, "-day average")

  draws %>%
    select(scenario, si_scenario, method, estimand_window, replicate, index,
           lower, upper, true_R_instant, true_R_window) %>%
    pivot_longer(
      cols = c(true_R_instant, true_R_window),
      names_to = "target",
      values_to = "target_R"
    ) %>%
    mutate(
      target = if_else(
        target == "true_R_instant", "Instantaneous", window_target_label
      ),
      own_estimand = (target == "Instantaneous") == (estimand_window == 1L)
    ) %>%
    filter(is.finite(target_R), is.finite(lower), is.finite(upper)) %>%
    group_by(scenario, si_scenario, method, target, own_estimand, index) %>%
    summarise(
      target_R = first(target_R),
      n_replicates = n(),
      coverage = mean(lower <= target_R & target_R <= upper),
      mean_ci_width = mean(upper - lower),
      .groups = "drop"
    ) %>%
    arrange(scenario, si_scenario, method, target, index)
}

# ------------------------------------------------------------------------------
# Day-level mean estimate across replicates
# ------------------------------------------------------------------------------
# Collapses the replicate dimension to one trajectory per arm, which is what
# plot_scenario_truth() overlays on the true Rt paths. Every method is kept by
# default; the correctly specified SI is the default cell so that any gap from an
# arm's own target is window averaging or estimator bias rather than SI
# misspecification. Pass method_pattern to narrow to one estimator.
#
# method_family separates the estimator from the arm: `method` is
# "EpiEstim (7-day window)" or "EpiLPS", and the plot needs the estimator alone
# because EpiLPS and the 1-day EpiEstim arm share an estimand and so share a
# colour, and are told apart by linetype.
#
# Note that `draws` is already trimmed to index >= sim_score_start_index(), so
# these trajectories legitimately start later than the truth they are drawn over.
sim_mean_estimate_by_day <- function(draws,
                                     si_scenario_shown = "correct",
                                     method_pattern = NULL) {
  selected <- draws %>%
    filter(si_scenario == si_scenario_shown)

  if (!is.null(method_pattern)) {
    selected <- selected %>% filter(grepl(method_pattern, method))
  }

  selected %>%
    mutate(
      method_family = if_else(grepl("^EpiLPS", method), "EpiLPS", "EpiEstim")
    ) %>%
    group_by(scenario, method, method_family, estimand_window, index) %>%
    summarise(
      n_replicates = n(),
      mean_R = mean(R, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(scenario, method_family, estimand_window, index)
}
