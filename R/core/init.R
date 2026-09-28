# ==============================================================================
# init.R
# Bootstrap for every driver script. Source this first, from the project root:
#
#   source(file.path("R", "core", "init.R"))
#   require_packages(c("EpiEstim", "EpiLPS", "dplyr", ...))
#   source(file.path("R", "core", "config.R"))   # only drivers that want main.R's settings
#   source_project()
#
# Nothing here sets a seed, creates a directory or reads data. Those are driver
# decisions, made in the driver where they can be seen.
# ==============================================================================

# Every path in the project is relative to r-proj/, so fail early and clearly
# if the working directory is anything else.
if (!file.exists(file.path("R", "core", "init.R"))) {
  stop("Run this from the r-proj/ directory (R/core/init.R not found).")
}

# Check that packages are installed, then attach them. library() must happen
# BEFORE any module function is called: metrics.R and sim_study.R use bare
# dplyr verbs, so filter() would otherwise resolve to stats::filter at call
# time. callr is added on Windows whenever EpiLPS is requested, because the
# fresh-process job runner needs it there (forking is unavailable).
require_packages <- function(pkgs) {
  if (.Platform$OS.type == "windows" && "EpiLPS" %in% pkgs) {
    pkgs <- union(pkgs, "callr")
  }
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop("Install first: ", paste(missing, collapse = ", "))
  }
  suppressPackageStartupMessages({
    for (p in pkgs) library(p, character.only = TRUE)
  })
  invisible(pkgs)
}

# The module folders, in dependency order. core/ knows nothing about
# epidemiology; si/, estimators/, simulation/ and scoring/ are leaves that use
# only core/ and packages; studies/ composes the leaves into experiments;
# plots/ consumes the data frames studies/ and scoring/ produce.
project_module_dirs <- c("core", "data", "si", "estimators", "simulation",
                         "scoring", "studies", "plots")

# Source every module under R/. config.R is skipped: it is settings, not code,
# and a driver sources it explicitly if it wants main.R's settings. Every
# module is pure function (and registry) definitions, so sourcing all of them
# is harmless for a driver that only needs some.
source_project <- function(root = "R") {
  skip <- c("init.R", "config.R")
  for (d in project_module_dirs) {
    files <- sort(list.files(file.path(root, d), pattern = "\\.[Rr]$",
                             full.names = TRUE))
    for (f in files) {
      if (basename(f) %in% skip) next
      source(f, local = FALSE)
    }
  }
  invisible(TRUE)
}

# paths.R is needed by config.R, which is sourced before source_project(), so
# it is loaded here.
source(file.path("R", "core", "paths.R"), local = FALSE)
