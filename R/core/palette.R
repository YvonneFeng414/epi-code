# ==============================================================================
# palette.R
# Colour vectors shared across figures, so the same entity carries the same
# colour from one figure to the next. Slots 1 and 2 of the validated
# categorical palette (blue, orange): adjacent-pair CVD separation and
# normal-vision separation both clear their floors on a light surface, and
# both sit above 3:1 contrast against it.
# ==============================================================================

# The two SI moment sweeps of the misspecification studies.
scenario_colors <- c(vary_sd = "#2a78d6", vary_mean = "#eb6834")

# The five assumed-SI families (see si/serial_interval.R). gamma_discr is the
# data-generating SI and is drawn in near-black so it reads as the reference.
si_family_colors <- c(gamma_discr = "#333333", gamma_bin = "#2a78d6",
                      lnorm = "#eb6834", weibull = "#1a9e5b", unif = "#9147c7")

# The two matched arms of the self-consistency study. Two keyings of the same
# two colours because the figures key on different columns: `method` is the
# full arm name, `method_family` is the estimator alone.
arm_levels <- c("EpiEstim (1-day window)", "EpiLPS")
arm_colors <- stats::setNames(c("#2a78d6", "#eb6834"), arm_levels)
estimator_colors <- stats::setNames(c("#2a78d6", "#eb6834"), c("EpiEstim", "EpiLPS"))
