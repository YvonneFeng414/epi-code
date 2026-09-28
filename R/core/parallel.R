# ==============================================================================
# parallel.R
# Worker-process helpers. The EpiLPS studies never call the compiled backend
# repeatedly in the parent process: every fit, or every replicate's batch of
# fits, runs in a child. On Unix that is a fork (parallel::mclapply /
# mcparallel), exactly as the original scripts did it. Windows cannot fork, so
# the same calls dispatch to a PSOCK cluster or a callr child instead.
#
# Workers on Windows are bootstrapped with source_project(), so any project
# function is callable inside FUN. Anything else FUN needs must be passed
# through its arguments or captured in a non-global closure environment - a
# free variable that lives in the parent's global environment does not exist
# on a PSOCK worker.
# ==============================================================================

# Number of worker processes: `cap` at most, leaving `reserve` cores to the rest
# of the machine, never below `floor`. The studies use floor = 2 because
# mclapply with mc.cores = 1 evaluates in the parent and defeats the
# fresh-process safeguard.
default_n_cores <- function(cap = 4L, reserve = 1L, floor = 1L) {
  available <- parallel::detectCores(logical = TRUE)
  if (is.na(available)) available <- 1L
  max(as.integer(floor), min(as.integer(cap), available - as.integer(reserve)))
}

# What a Windows worker runs before its first task.
.bootstrap_worker <- function(wd, packages) {
  setwd(wd)
  source(file.path("R", "core", "init.R"))
  for (p in packages) suppressPackageStartupMessages(library(p, character.only = TRUE))
  source_project()
  invisible(NULL)
}

# lapply(X, FUN) across workers.
#
# Unix: mclapply with mc.preschedule = FALSE - one fork per element, dispatched
# on demand. This is byte-for-byte the call the original studies made.
# Windows: a PSOCK cluster of n_cores workers, load-balanced. Workers persist
# across elements, so a FUN that calls EpiLPS::estimR() will call it repeatedly
# in one process; epilps_si_misspec.R's note (45 consecutive calls clean) is the
# basis for accepting that.
run_parallel <- function(X, FUN, n_cores, packages = c("EpiEstim", "EpiLPS")) {
  if (.Platform$OS.type != "windows") {
    return(parallel::mclapply(X, FUN, mc.cores = n_cores, mc.preschedule = FALSE))
  }

  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterCall(cl, .bootstrap_worker, getwd(), packages)
  parallel::parLapplyLB(cl, X, FUN)
}

# One call of fun(args) in a fresh child process, returning its value. Used for
# the single calibration fits (overdispersion, MALA timing) that the studies
# run before the grid. `fun` must take everything it needs through `args`.
run_in_child <- function(fun, args = list(), packages = c("EpiEstim", "EpiLPS")) {
  if (.Platform$OS.type != "windows") {
    child <- parallel::mcparallel(do.call(fun, args))
    return(parallel::mccollect(child)[[1L]])
  }

  callr::r(
    func = function(fun, args, wd, packages) {
      setwd(wd)
      source(file.path("R", "core", "init.R"))
      for (p in packages) suppressPackageStartupMessages(library(p, character.only = TRUE))
      source_project()
      do.call(fun, args)
    },
    args = list(fun = fun, args = args, wd = getwd(), packages = packages),
    spinner = FALSE,
    show = FALSE
  )
}
