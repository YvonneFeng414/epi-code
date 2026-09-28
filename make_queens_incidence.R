################################################################################
# make_queens_incidence.R
#
# Build results/queens_rt/queens_daily_incidence.csv on its own.
#
# This is Section 2 of main.R lifted out and nothing else. main.R writes the
# same file on its way through a much longer pipeline (Rt fits, sensitivity
# sweeps, the coverage study, the forecast backtest), so regenerating the
# incidence series used to mean running all of that. This script does only the
# load-validate-write step.
#
# WHY IT USES config.R AND data_loader.R RATHER THAN RE-IMPLEMENTING THEM
#
# Five scripts read this CSV:
#
#   truth_source_main.R        epilps_si_misspec.R
#   epilps_selfcheck.R         epilps_si_misspec_mala.R
#   epifilter_si_misspec.R
#
# All of their results are conditional on its exact contents, and two of them
# (epilps_si_misspec.R and its MALA twin) are a matched pair whose comparison
# only holds because they were driven by identical data. A second, independent
# implementation of the loader would be free to drift from main.R's - a
# different date parse, a different first-positive-day rule, a different
# rounding - and the drift would surface as changed downstream numbers with no
# obvious cause. So the settings and the loader are taken from the same files
# main.R uses, and this script contributes no data logic of its own.
#
# SAFETY: the file is compared against what is already on disk before anything
# is written. If the contents would change, the script says so loudly and
# refuses to overwrite unless you pass --force, because a silent change would
# invalidate every result listed above without touching those scripts.
#
# Run from r-proj/:
#   Rscript make_queens_incidence.R           # write if absent, verify if present
#   Rscript make_queens_incidence.R --force   # overwrite even if contents differ
################################################################################

rm(list = ls())

args  <- commandArgs(trailingOnly = TRUE)
force <- "--force" %in% args

# ==============================================================================
# 0. Packages and modules
# ==============================================================================
# dplyr is all data_loader.R needs; main.R's other packages belong to the
# analysis sections this script does not run.

source(file.path("R", "core", "init.R"))
require_packages("dplyr")
source(file.path("R", "core", "config.R"))   # data_file, state_abbr, county_name,
                                             # burn_in_days, output_dir
source_project()                             # load_incidence_data()

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- project_paths()$incidence_csv

cat("\n===== Building the Queens daily incidence series =====\n")
cat(sprintf("source   : %s\n", normalizePath(data_file, mustWork = FALSE)))
cat(sprintf("location : %s, %s\n", county_name, state_abbr))
cat(sprintf("target   : %s\n", out_file))

# ==============================================================================
# 1. Load and validate
# ==============================================================================
# load_incidence_data() does the work and stops on any of: missing columns,
# unparseable dates, duplicate dates, a gap in the daily sequence, non-finite
# counts, negative counts from cumulative-case revisions, an all-zero series,
# or too few days left for burn_in_days. Reaching the next line means the
# series passed all of them.

incidence_data <- load_incidence_data(
  data_file    = data_file,
  state_abbr   = state_abbr,
  county_name  = county_name,
  burn_in_days = burn_in_days
)

cat(sprintf("\nloaded %d days: %s to %s\n",
            nrow(incidence_data),
            as.character(min(incidence_data$dates)),
            as.character(max(incidence_data$dates))))
cat(sprintf("daily cases: min %d, median %.0f, max %d, total %d\n",
            min(incidence_data$I), stats::median(incidence_data$I),
            max(incidence_data$I), sum(incidence_data$I)))
cat(sprintf("series is trimmed to start at the first positive-incidence day (I = %d)\n",
            incidence_data$I[1L]))

# ==============================================================================
# 2. Compare against what is already on disk
# ==============================================================================
# write.csv() to a temporary file first, so the comparison is on exactly what
# would be written rather than on a re-parse of the data frame. The digest is
# taken over the LINES, not the raw bytes: R writes CRLF on Windows and LF
# elsewhere, and a line-ending difference is not a change of contents.

tmp_file <- tempfile(fileext = ".csv")
write.csv(incidence_data, tmp_file, row.names = FALSE)
on.exit(unlink(tmp_file), add = TRUE)

lines_md5 <- function(path) {
  lines <- readLines(path, warn = FALSE)
  digest_file <- tempfile()
  on.exit(unlink(digest_file), add = TRUE)
  writeBin(charToRaw(paste(lines, collapse = "\n")), digest_file)
  unname(tools::md5sum(digest_file))
}

new_md5 <- lines_md5(tmp_file)

if (!file.exists(out_file)) {
  file.copy(tmp_file, out_file, overwrite = FALSE)
  cat(sprintf("\nWROTE %s (did not exist)\n", out_file))
  cat(sprintf("md5: %s\n", new_md5))
} else {
  old_md5 <- lines_md5(out_file)

  if (identical(new_md5, old_md5)) {
    cat(sprintf("\nUNCHANGED - the file on disk already matches (md5 %s)\n", new_md5))
    cat("Nothing was written. Downstream results stay valid.\n")
  } else {
    old <- read.csv(out_file, stringsAsFactors = FALSE)
    cat("\n**************************************************************\n")
    cat("CONTENTS WOULD CHANGE\n")
    cat(sprintf("  on disk : %d rows, md5 %s\n", nrow(old), old_md5))
    cat(sprintf("  rebuilt : %d rows, md5 %s\n", nrow(incidence_data), new_md5))

    if (nrow(old) == nrow(incidence_data)) {
      differing <- which(old$I != incidence_data$I |
                         as.character(old$dates) != as.character(incidence_data$dates))
      cat(sprintf("  %d of %d rows differ; first at row %d\n",
                  length(differing), nrow(old),
                  if (length(differing)) differing[1L] else NA_integer_))
    }

    cat("\nEvery downstream result was produced from the file on disk:\n")
    cat("  truth_source_main.R    epilps_si_misspec.R\n")
    cat("  epilps_selfcheck.R     epilps_si_misspec_mala.R\n")
    cat("  epifilter_si_misspec.R\n")
    cat("Overwriting invalidates all of them, and the paired LPSMAP/LPSMALA\n")
    cat("comparison in particular depends on both arms sharing this input.\n")
    cat("**************************************************************\n")

    if (!force) {
      stop("Refusing to overwrite. Re-run with --force if the change is intended.")
    }

    file.copy(tmp_file, out_file, overwrite = TRUE)
    cat(sprintf("\n--force given: OVERWROTE %s\n", out_file))
    cat("Re-run the downstream scripts before comparing any results.\n")
  }
}

cat(sprintf("\nDone. %s\n", normalizePath(out_file, mustWork = FALSE)))
