# simulations/ — EpiEstim on fully synthetic epidemics

The pure-synthetic side of `archive/Epiestim_simu/`, ported onto the shared
modules in `R/`. Nothing here reads data: every epidemic is generated from an
analytic Rt shape (`epiestim_scenario_table` in `R/simulation/scenarios.R`)
seeded with a single `I0`, so the truth is known exactly. The Queens-facing
part of the same archive is the top-level `epiestim_si_misspec.R`.

**Run from `r-proj/`**, e.g. `Rscript simulations/sim_scenarios.R`. Paths are
relative to `r-proj/`; the init guard stops any driver started elsewhere.
Settings shared across the folder are in `sim_config.R`; the true SI and the
relative-error grids come from `R/core/config.R`. All figures are PDF.

Every driver has a `quick <- FALSE` switch near the top. `TRUE` runs a small
subset with 20 replicates in a minute or two — use it to check the pipeline
before a full run.

## The drivers

| Script | What it does | Writes to `results/simulations/…` | Full run |
|---|---|---|---|
| `sim_scenarios.R` | 27 Rt shapes at T = 100 (+ 3 constant at T = 250) × I0 ∈ {100, 2000, 10000}, EpiEstim under the correct SI, 1-day and 7-day windows each scored against its own estimand | `scenarios/`: `scenario_metrics.csv`, `scenario_metrics_formatted.csv` (appendix layout), `scenario_by_day.csv`, `scenario_truth.csv`, `scenario_truth.pdf` | ~40 min |
| `sim_si_misspec.R` | Nine synthetic targets × the full SI grid (two moment sweeps × five families), both windows; Euler–Lotka reference per target | `si_misspec/`: `si_misspec_metrics.csv`, `si_misspec_by_day.csv`, `si_family_comparison.csv`, `euler_lotka_reference.csv`, `si_misspec_pages.pdf`, `si_family_pages.pdf`, `si_family_kernels.pdf` | hours — read the projection it prints after the first target |
| `sim_rt_trajectories.R` | Five targets; Rt trajectory pages under the correct SI (one epidemic + Monte Carlo average) and under ±40 / ±80 % errors on the mean and sd | `rt_trajectories/`: `rt_trajectory_summary.csv`, `rt_misspec_trajectory_summary.csv`, `rt_trajectory_pages.pdf`, `rt_misspec_trajectory_pages.pdf` | ~15 min |
| `sim_bias_convergence.R` | Constant Rt = 1.2 to T = 200, assumed sd wrong: per-day bias, differenced against the correct sd, versus the Euler–Lotka asymptote | `bias_convergence/`: `bias_convergence.csv`, `bias_convergence_summary.csv`, `bias_convergence.pdf` | minutes |
| `sim_britton_validation.R` | Deterministic exponential growth: EpiEstim's R-hat equals the Euler–Lotka value to ~1e-7 across 56 assumed sds; eight `stopifnot` premises; reconciliation with the stochastic study | `britton_validation/`: `britton_validation.csv`, `britton_validation.pdf`, `britton_reconciliation.csv`, `britton_reconciliation.pdf` | < 1 min |

The drivers are independent. The one soft dependency: `sim_britton_validation.R`
fills its `bias_stochastic` column from `si_misspec/si_misspec_metrics.csv` when
that exists, so run it (again) after `sim_si_misspec.R` to complete the
reconciliation figure.

## What changed from the archive

Same Rt shapes, same grids, same `n_sim = 2000`, same Euler–Lotka machinery.
Harmonised onto `R/`:

- **DGP and fits**: `simulate_renewal_incidence()` (Poisson, `I0` on day 1),
  `fit_epiestim_grid()` (the harness adapter), `fit_si_grid()` for the grid.
- **Scoring**: `score_interval_fits()` — coverage per replicate first (so the
  MCSE is valid), RMSE, `half_over_sd`. Scoring starts at `max(window) + 1`
  (day 8 with the 7-day arm present) so both arms are scored on the same days;
  the archive scored from day 2. `sim_bias_convergence.R` scores from day 2
  because the early days are its subject.
- **Each window scored against its own estimand** (`trailing_mean()`), which
  the archive did implicitly by fitting daily windows only.
- **SI**: `make_si()` at lag 30 (lag 399 in the Britton check), the project's
  five families (`gamma_discr`, `gamma_bin`, `lnorm`, `weibull`, `unif`) in
  place of the archive's four.
- **Seeds**: one integer per simulated cell, not one stream from `set.seed(123)`.
- **Dropped**: Scenario 4 (Queens real Rt), the extinction filter (dead code for
  synthetic targets), Word/`flextable` output, PNG output.

Numbers will therefore not match the archive's appendices to the digit; the
design does, and the rows are comparable with the Queens misspecification
studies at the top level.

## The archive's own results

The outputs the archive scripts produced (their CSVs, PDFs, PNGs, Word appendices
and run logs) are copied — not moved — into `results/simulations/<study>/archive/`
so each ported study sits next to the numbers it replaces: Appendix 2 under
`scenarios/`, Appendices 3–5 under `si_misspec/`, and the trajectory,
bias-convergence and Britton outputs under their studies. Column names and
scoring differ from the new files (see above), so compare shapes and orderings
rather than digits. `results/` is git-ignored; `archive/Epiestim_simu/` remains
the frozen original.

## Tests

```
Rscript tests/test_sim_epiestim.R
```

Carries the archive's `make_step_R()`, `make_gradual_R()` and
`deterministic_renewal()` as oracles for the new constructors, and checks the
per-arm scoring, trajectory summaries, `convergence_day()`, the PDF device
switch and the appendix table formatter.
