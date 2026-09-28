# ==============================================================================
# incidence_io.R
# Reader for the finished daily incidence file (dates, I) that
# make_queens_incidence.R / main.R write and every downstream study reads.
#
# load_incidence_data() in data_loader.R cannot be used for this: it parses the
# raw fact_*.csv, with cases/COUNTY_NAME columns and a cumulative series to
# difference. This file is already the two-column product of that step.
# ==============================================================================

# Read and validate. Stops on a missing column, an unparseable date, a negative
# or non-finite count, or a gap in the daily sequence. Counts come back as
# integers and dates as Date.
read_incidence_csv <- function(path) {
  if (!file.exists(path)) {
    stop("Cannot find the incidence file: ", normalizePath(path, mustWork = FALSE))
  }

  incidence_data <- read.csv(path, stringsAsFactors = FALSE)

  missing_columns <- setdiff(c("dates", "I"), names(incidence_data))
  if (length(missing_columns) > 0L) {
    stop("The incidence file is missing column(s): ",
         paste(missing_columns, collapse = ", "))
  }

  incidence_data$dates <- as.Date(incidence_data$dates)

  if (anyNA(incidence_data$dates)) {
    stop("At least one date in the incidence file could not be parsed.")
  }
  if (any(!is.finite(incidence_data$I)) || any(incidence_data$I < 0)) {
    stop("Incidence must be non-negative and finite.")
  }

  expected_dates <- seq(min(incidence_data$dates), max(incidence_data$dates),
                        by = "day")
  if (length(incidence_data$dates) != length(expected_dates) ||
      !all(incidence_data$dates == expected_dates)) {
    stop("The incidence file does not contain a complete daily date sequence.")
  }

  incidence_data$I <- as.integer(round(incidence_data$I))
  incidence_data
}
