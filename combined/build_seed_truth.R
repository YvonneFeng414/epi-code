################################################################################
# The seed study's true Rt trajectory, written out once for the whole repo.
#
# WHY THIS IS A SEPARATE STEP. Four scripts (the EpiEstim arm, EpiFilter, and the
# two EpiLPS arms) must be driven by a bit-identical truth, or the comparison
# between them measures the truth as much as the estimator. Recomputing it in
# each of them would invite exactly the drift make_queens_incidence.R's header
# warns about - a different SI construction or a different hole rule in one file
# and the four arms silently diverge. So it is built here and read from
# combined/seed_aligned_truth.csv everywhere else.
#
# The pipeline is the seed study's, reproduced exactly:
#
#   1. load_incidence_data()  raw fact_*.csv -> daily cases, validated, trimmed
#                             to the first positive day
#   2. estimate_real_rt()     EpiEstim parametric_si, 1-day windows from t = 2,
#                             everything before the first Median(R) <= 5 dropped
#   3. hole removal           days with true Rt < 0.2 (the seed's rule, not the
#                             incidence == 0 rule this repo used before; on this
#                             series both select 2020-07-08 and the output is the
#                             same, but the seed's rule is what is authoritative
#                             here)
#
# NOTE it is NOT the same as results_truth_source_discr_si/plugin_truth_paths.csv's
# true_R_epiestim column, which differs by up to 2.2e-04: that one is built with
# non_parametric_si over make_si()'s lag-30 truncated-and-renormalised pmf, while
# the seed lets EpiEstim construct its own full-support gamma internally.
#
# Run from queens_rt/:  Rscript combined/build_seed_truth.R
################################################################################

suppressPackageStartupMessages({ library(EpiEstim); library(dplyr) })

seed_dir  <- "/Users/yixuanfeng/Desktop/Epiestim_simu/seed"
data_file <- "/Users/yixuanfeng/Desktop/Epiestim_simu/fact_2020-12-01.csv"
out_file  <- "combined/seed_aligned_truth.csv"

# The seed's own settings, read from its config rather than retyped.
for (f in c(file.path(seed_dir, "config.R"), file.path(seed_dir, "queens_functions.R"))) {
  if (!file.exists(f)) stop("Cannot find ", f)
}
source(file.path(seed_dir, "config.R"))          # mean_si, sd_si, rt_cap, hole
source(file.path(seed_dir, "queens_functions.R"))# load_incidence_data, estimate_real_rt

hole_threshold <- 0.2   # the seed driver's literal

incidence_data <- load_incidence_data(
  data_file    = data_file,
  state_abbr   = queens_state_abbr,
  county_name  = queens_county_name,
  burn_in_days = queens_burn_in_days
)

R_t_full <- estimate_real_rt(
  incidence_data = incidence_data,
  mean_si = mean_si, sd_si = sd_si,
  rt_cap = rt_cap, nd_real = 1
)

hole_days <- which(R_t_full < hole_threshold)

cat("\n===== Seed-aligned true Rt =====\n")
cat(sprintf("incidence      : %d days, %s to %s\n", nrow(incidence_data),
            as.character(min(incidence_data$dates)), as.character(max(incidence_data$dates))))
cat(sprintf("after rt_cap %g : %d days\n", rt_cap, length(R_t_full)))
for (d in hole_days) {
  cat(sprintf("hole removed   : day %d, true Rt %.5f, %d reported cases\n", d, R_t_full[d],
              incidence_data$I[d + (nrow(incidence_data) - length(R_t_full))]))
}

# The other arms identify the hole from the incidence (I == 0) because their
# truths are smooth and never dip below 0.2. Assert the two rules still pick the
# same day - if they ever diverge the arms would silently score different days.
repo_inc  <- read.csv(file.path("results_queens_rt_kimi", "queens_daily_incidence.csv"),
                      stringsAsFactors = FALSE)
repo_hole <- which(repo_inc$I[seq.int(16L, nrow(repo_inc))] == 0)
if (!identical(as.integer(hole_days), as.integer(repo_hole))) {
  stop("Hole rules disagree: seed Rt<0.2 picks ", paste(hole_days, collapse = ","),
       " but incidence==0 picks ", paste(repo_hole, collapse = ","),
       ". The arms would score different days.")
}
cat(sprintf("hole rules agree: seed Rt<0.2 and incidence==0 both select day %s\n",
            paste(hole_days, collapse = ",")))

R_t_queens <- R_t_full[-hole_days]

cat(sprintf("truth          : %d days, Rt %.4f to %.4f\n",
            length(R_t_queens), min(R_t_queens), max(R_t_queens)))

write.csv(data.frame(index = seq_along(R_t_queens), true_R = R_t_queens),
          out_file, row.names = FALSE)
cat(sprintf("\nWrote %s\n", out_file))
