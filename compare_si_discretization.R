################################################################################
# Does the self-consistency result depend on how the serial interval is
# discretized?
#
# make_si() (EpiEstim::discr_si) and EpiLPS::Idist() do NOT discretize the same
# continuous distribution:
#
#   discr_si - gamma fitted to mean - 1, shifted a day, Cori triangular kernel
#   Idist    - gamma fitted to the mean directly, binned at +/- 0.5
#
# Both realise mean/sd ~ 7.5/3.4, but they differ by total variation 0.019 and by
# 3.2x in the lag-1 cell, which is the weight driving the most recent days'
# contribution to infectiousness. EpiLPS's own workflow uses Idist, so if the
# headline number moves between the two, the finding does not transfer to anyone
# following that workflow.
#
# This script does no fitting. Run truth_source_main.R once per scheme first:
#
#   si_discretization <- "discr_si"   -> results/truth_source_discr_si/
#   si_discretization <- "idist"      -> results/truth_source_idist/
#
# Run from r-proj/:  Rscript compare_si_discretization.R
# Output: results/si_discretization/
################################################################################

rm(list = ls())

source(file.path("R", "core", "init.R"))
require_packages(c("EpiEstim", "EpiLPS", "dplyr", "tidyr", "ggplot2"))
source_project()

paths <- project_paths()
schemes <- c(discr_si = paths$truth_source_dir("discr_si"),
             idist    = paths$truth_source_dir("idist"))

output_dir <- paths$results_dir("si_discretization")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 1. Load both runs
# ==============================================================================

missing_runs <- schemes[!file.exists(file.path(schemes, "truth_source_metrics.csv"))]
if (length(missing_runs) > 0L) {
  stop("Missing run output for: ", paste(names(missing_runs), collapse = ", "),
       ". Run truth_source_main.R under each si_discretization first.")
}

metrics <- bind_rows(lapply(names(schemes), function(s) {
  read.csv(file.path(schemes[[s]], "truth_source_metrics.csv")) %>%
    mutate(scheme = s, .before = 1L)
}))

truth_paths <- bind_rows(lapply(names(schemes), function(s) {
  read.csv(file.path(schemes[[s]], "plugin_truth_paths.csv")) %>%
    mutate(scheme = s, .before = 1L)
}))

# The day-level bias share, recomputed per scheme the same way the single-run
# figure does it.
bias_share <- bind_rows(lapply(names(schemes), function(s) {
  arm_bias_share(
    read.csv(file.path(schemes[[s]], "truth_source_coverage_by_day.csv")),
    read.csv(file.path(schemes[[s]], "truth_source_estimates.csv"))
  ) %>%
    mutate(scheme = s, .before = 1L)
}))

combined <- metrics %>%
  left_join(bias_share, by = c("scheme", "scenario", "method")) %>%
  mutate(bias_share = 100 * mean_bias2 / MSE_ownEstimand)

write.csv(
  combined,
  file.path(output_dir, "si_discretization_metrics.csv"),
  row.names = FALSE
)

# ==============================================================================
# 2. The two SI vectors, on the record
# ==============================================================================

si_discr <- make_si(7.5, 3.4, 30L)
si_idist <- c(0, EpiLPS::Idist(mean = 7.5, sd = 3.4, dist = "gamma")$pvec)
max_lag <- max(length(si_discr), length(si_idist)) - 1L
pad <- function(x) c(x, rep(0, max_lag + 1L - length(x)))

si_table <- data.frame(
  lag = 0:max_lag,
  discr_si = pad(si_discr),
  idist = pad(si_idist)
) %>%
  mutate(difference = idist - discr_si)

write.csv(si_table, file.path(output_dir, "si_vectors.csv"), row.names = FALSE)

cat("\n===== The two serial intervals =====\n")
for (nm in c("discr_si", "idist")) {
  v <- si_table[[nm]]
  m <- si_moments(v)
  cat(sprintf("  %-9s support 0-%2d  realised mean/sd %.4f / %.4f  lag-1 mass %.5f\n",
              nm, max(si_table$lag[v > 0]), m[["mean"]], m[["sd"]], v[2L]))
}
cat(sprintf("  total variation distance %.4f, max gap %.4f at lag %d\n",
            0.5 * sum(abs(si_table$difference)),
            max(abs(si_table$difference)),
            si_table$lag[which.max(abs(si_table$difference))]))

# Sanity check: if the two runs did not score the same days, every comparison
# below is confounded and the caveat has to travel with the numbers.
scored <- combined %>% distinct(scheme, N_day_estimates)
if (n_distinct(scored$N_day_estimates) > 1L) {
  cat("\nWARNING: the two schemes scored different numbers of day-estimates:\n")
  print(as.data.frame(scored), row.names = FALSE)
  cat("Their metrics are therefore not computed over identical days.\n")
}

truth_span <- truth_paths %>%
  group_by(scheme) %>%
  summarise(days = n(), first_date = min(as.Date(date)), .groups = "drop")
cat("\nCommon range per scheme:\n")
print(as.data.frame(truth_span), row.names = FALSE)

# ==============================================================================
# 3. Delta table
# ==============================================================================

deltas <- combined %>%
  select(scheme, method, Coverage95_Rt_ownEstimand, MeanCIWidth,
         RMSE_ownEstimand, bias_share) %>%
  pivot_longer(-c(scheme, method), names_to = "metric") %>%
  pivot_wider(names_from = scheme, values_from = value) %>%
  mutate(change = idist - discr_si)

cat("\n===== Effect of the SI discretization on each arm =====\n")
print(as.data.frame(deltas), row.names = FALSE, digits = 4)

# ==============================================================================
# 4. Figure
# ==============================================================================

panel_source <- combined %>%
  transmute(
    scheme = factor(scheme, levels = names(schemes)),
    method = factor(method, levels = arm_levels),
    coverage = Coverage95_Rt_ownEstimand,
    coverage_mcse = MCSE_Coverage95_Rt_ownEstimand,
    width = MeanCIWidth,
    rmse = RMSE_ownEstimand,
    bias_share
  )

metric_panels <- build_arm_metric_panels(panel_source, keys = c("scheme", "method"))

ggsave(
  filename = file.path(output_dir, "si_discretization_metrics.png"),
  plot = plot_arm_metrics(
    metric_panels,
    references = arm_metric_references(),
    title = "Does the serial-interval discretization change the result?",
    subtitle = paste(
      "Same study run twice. 'discr_si' = EpiEstim's shifted-gamma kernel;",
      "'idist' = EpiLPS's +/-0.5 binning.\nWithin each run both estimators share",
      "one SI. Each arm is still scored only against its own plug-in truth, so",
      "\nthe two colours are not a ranking - what matters is whether a colour",
      "moves between the two x groups."
    ),
    method_colors = arm_colors,
    x_var = "scheme"
  ),
  width = 10, height = 7.5, dpi = 300
)

cat("\nComparison written to:\n")
cat(normalizePath(output_dir), "\n")
