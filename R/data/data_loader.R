# ==============================================================================
# data_loader.R
# Load and validate the county-level incidence data (Section 2 of the original
# script). Returns a two-column data frame (dates, I) trimmed to start at the
# first positive-incidence day.
# ==============================================================================

load_incidence_data <- function(data_file, state_abbr, county_name, burn_in_days) {
  if (!file.exists(data_file)) {
    stop("Cannot find data file: ", normalizePath(data_file, mustWork = FALSE))
  }

  raw_data <- read.csv(
    data_file,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required_columns <- c(
    "date", "cases", "State", "COUNTY_NAME", "STATE_NAME"
  )
  missing_columns <- setdiff(required_columns, names(raw_data))
  if (length(missing_columns) > 0L) {
    stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
  }

  clean_data <- raw_data %>%
    mutate(date = as.Date(date, format = "%m/%d/%Y")) %>%
    arrange(STATE_NAME, COUNTY_NAME, date) %>%
    group_by(STATE_NAME, COUNTY_NAME) %>%
    mutate(daily_cases = c(first(cases), diff(cases))) %>%
    ungroup()

  if (anyNA(clean_data$date)) {
    stop("At least one date could not be parsed.")
  }

  incidence_data <- clean_data %>%
    filter(State == state_abbr, COUNTY_NAME == county_name) %>%
    arrange(date) %>%
    transmute(
      dates = date,
      I = as.numeric(daily_cases)
    )

  if (nrow(incidence_data) == 0L) {
    stop("No rows found for ", county_name, ", ", state_abbr, ".")
  }

  if (anyDuplicated(incidence_data$dates)) {
    stop("Duplicate dates were found for the selected county.")
  }

  expected_dates <- seq(
    min(incidence_data$dates),
    max(incidence_data$dates),
    by = "day"
  )
  # identical() would compare storage mode too (dplyr's Date column is stored
  # as double, seq.Date()'s result as integer) and spuriously fail even when
  # the dates match exactly - compare values instead.
  if (length(incidence_data$dates) != length(expected_dates) ||
      !all(incidence_data$dates == expected_dates)) {
    stop("The selected county does not have a complete daily date sequence.")
  }

  if (any(!is.finite(incidence_data$I))) {
    stop("Incidence contains non-finite values.")
  }

  if (any(incidence_data$I < 0)) {
    bad_dates <- incidence_data$dates[incidence_data$I < 0]
    stop(
      "Negative daily counts were created by revisions to cumulative cases. ",
      "Resolve them before Rt estimation. First affected date: ", bad_dates[1L]
    )
  }

  first_positive <- which(incidence_data$I > 0)[1L]
  if (is.na(first_positive)) {
    stop("The selected series contains no positive incidence.")
  }

  incidence_data <- incidence_data[first_positive:nrow(incidence_data), , drop = FALSE]
  incidence_data$I <- as.integer(round(incidence_data$I))
  row.names(incidence_data) <- NULL

  n_obs <- nrow(incidence_data)
  if (n_obs <= burn_in_days + 1L) {
    stop("Not enough observations after trimming for the chosen burn-in.")
  }

  incidence_data
}
