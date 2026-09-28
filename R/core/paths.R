# ==============================================================================
# paths.R
# Every file location the project reads or writes, in one place. All paths are
# relative to r-proj/ (see init.R's guard).
#
# Three inputs are produced by one script and consumed by others, so their
# locations are named here rather than spelled out in each consumer:
#
#   incidence_csv        written by make_queens_incidence.R / main.R, read by
#                        truth_source_main.R and every misspecification study
#   truth_paths_csv(s)   written by truth_source_main.R under scheme s, read by
#                        the EpiLPS misspecification studies and the selfcheck
#   results_dir(name)    one output directory per study
# ==============================================================================

project_paths <- function() {
  results_root <- "results"
  list(
    # The raw county-level cumulative-case file lives at the root of r-proj/.
    data_file = "fact_2020-12-01.csv",

    results_root = results_root,
    results_dir = function(name) file.path(results_root, name),

    incidence_csv = file.path(results_root, "queens_rt", "queens_daily_incidence.csv"),

    truth_source_dir = function(scheme) {
      file.path(results_root, paste0("truth_source_", scheme))
    },
    truth_paths_csv = function(scheme = "discr_si") {
      file.path(results_root, paste0("truth_source_", scheme), "plugin_truth_paths.csv")
    }
  )
}
