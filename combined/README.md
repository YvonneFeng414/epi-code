# Four Rt estimators under the same SI misspecification grid

`combined_si_misspec_panel.pdf` / `.png` put the seed study's EpiEstim result on
the same axes as the three estimators in this repo. Three metrics as columns
(Bias, RMSE, 95% CI coverage), the two misspecification arms as rows with the
**mean** arm on top, panels lettered A–F, estimator key along the top. Each
estimator carries its own shape as well as its own colour. Y scales are shared
down each column.

**Reading the two EpiLPS curves.** In panels A, B, D and E the MAP and MALA
curves coincide to well under the line width, so three devices are used to stop
that reading as a missing series:

1. distinct marker shapes (▲ MAP, ◆ MALA) as well as distinct colours;
2. MALA drawn as a long dash, so the solid green MAP line shows through the gaps
   instead of being painted over;
3. an in-panel note giving the measured maximum discrepancy.

The discrepancies are 0.0003 (Bias, both arms), 0.0021 and 0.0019 (RMSE) —
between 0.5% and 2.3% of the range each panel spans, i.e. not resolvable at any
print size. This is expected, not a bug: the two arms share the simulated data
*and* the point estimator (both report the posterior median `Rq0.50`), so only
the intervals can differ, and only the Coverage panels can show it.

| file | contents |
|---|---|
| `combined_si_misspec_panel.pdf` / `.png` | the 6-panel figure, 12 × 7.2 in |
| `combined_si_misspec_metrics.csv` | 72 rows: 4 methods × 2 arms × 9 grid points |
| `make_combined_panel.R` | the script that builds both |

Sources, all filtered to the `gamma_bin` assumed-SI family:

| method | source | observation model | replicates |
|---|---|---|---|
| EpiEstim | `seed/sensitivity_results_queens_nohole.csv` | Poisson | 2000 |
| EpiLPS (MAP) | `results_epilps_si_misspec/` | **negative binomial**, ρ = 17.36 | 100 |
| EpiLPS (MALA) | `results_epilps_si_misspec_mala/` | **negative binomial**, ρ = 17.36 | 100 |
| EpiFilter | `results_epifilter_si_misspec/` (smoother arm) | Poisson | 100 |

**Evaluated on the common grid: −60% to +60%, 7 points.** The three repo
estimators were run out to ±80%, but the seed study (maintained outside this
repo) currently stops at ±60%, so the ±80% column is dropped and every method is
averaged over the same settings. `restrict_to_common_grid` in the script controls
this; it intersects the four grids rather than hardcoding a range, so if the seed
sweep is regenerated at ±80% the restriction lifts by itself. Each run prints
what was dropped and each method's resulting range.

---

## Read this before reading the figure

**The four arms do not face the same data.** The two EpiLPS scripts simulate
negative-binomial incidence; EpiEstim and EpiFilter simulate Poisson. At a mean
of 1000 cases/day that is an observation SD of **242 against 32 — 7.7× the
noise**. The vertical separation between the EpiLPS pair and the other two in
the Coverage panels is therefore **confounded with the observation model and is
not an estimator ranking**.

The intuitive correction — "EpiLPS looks worse only because its data is noisier"
— is *wrong*, and measured to be wrong. Refitting EpiLPS on Poisson data drawn
from the same truth and the same seed days (20 replicates, correct SI) makes its
coverage **worse**, not better:

| data EpiLPS is fitted to | recovered ρ | Coverage | CI width | RMSE |
|---|---|---|---|---|
| negative binomial (as in the figure) | 17.9 | **0.727** | 0.1661 | 0.0719 |
| Poisson | 25.2 | **0.599** | 0.0700 | 0.0427 |

Poisson data makes the point estimate sharper (RMSE 0.072 → 0.043) but the
interval collapses faster still (0.166 → 0.070), so coverage falls. Putting all
four arms on Poisson would therefore *widen* the gap in the Coverage panels, not
close it. Note also that ρ comes out at 25 on data that is exactly Poisson, where
the true value is infinite: EpiLPS's "overdispersion" is absorbing B-spline
misfit, not only observation noise.

This cannot be fixed by fitting EpiLPS as a Poisson model, because the package
does not offer one — see *Can EpiLPS fit Poisson?* below.

What is safe to read across methods is the *shape* of each curve — how each
estimator degrades as the assumed SI moves — and what is safe to read within the
EpiLPS pair is the MAP/MALA gap, a controlled comparison on identical data.

Two further caveats:

- **Each estimator is scored against its own plug-in truth.** EpiEstim against an
  EpiEstim fit of the real series, EpiLPS against `true_R_epilps`, EpiFilter
  against its own smoother fit. Every truth is smooth in exactly the way its own
  method assumes, which flatters each method at its control point. This cuts
  *against* the EpiLPS result below, not for it.
- **`gamma_bin` is not exactly EpiEstim's own SI.** EpiEstim uses `discr_si`
  (= `gamma_discr`), a different discretisation of the same gamma; the two differ
  by 0.019 in total variation. EpiEstim's curve has the right shape to compare and
  one discretisation step of slack in its level.

---

## What the data says

### 1. The assumed mean is what matters; the assumed sd barely registers

Mean coverage across each 9-point sweep, against that method's own control point:

| method | sd sweep loss | mean sweep loss | ratio |
|---|---|---|---|
| EpiEstim | 0.028 | **0.177** | 6.4× |
| EpiFilter | 0.034 | **0.212** | 6.3× |
| EpiLPS (MAP) | 0.008 | **0.187** | 22.3× |
| EpiLPS (MALA) | 0.006 | **0.198** | 32.2× |

Every estimator loses **6–32× more coverage to a wrong mean than to a wrong sd**,
and the four agree on the size of the mean-sweep loss to within 0.035 (0.177 to
0.212) despite different observation models, different plug-in truths and a 20×
difference in replicate count. This is the most robust result in the figure. RMSE
says the same: over the sd sweep the worst/best ratio is 1.0 (both EpiLPS arms)
and 1.3 (EpiEstim, EpiFilter); over the mean sweep it is 2.2, 2.6 and 2.9.

The two EpiLPS arms are the most lopsided — essentially indifferent to the
assumed sd (loss 0.006–0.008) while losing as much as anyone to a wrong mean.

The practical reading: effort spent pinning down the serial interval's **mean** is
worth far more than effort spent on its **sd**, whichever estimator you use.

### 2. EpiLPS undercovers at the correct SI — and it is not the Laplace approximation alone

At `rel_error = 0`, where the fitted SI *is* the generating SI:

| method | Coverage | Bias | RMSE | CI width |
|---|---|---|---|---|
| EpiEstim | 0.950 | +0.0058 | 0.1135 | 0.4095 |
| EpiFilter | 0.950 | +0.0017 | 0.0626 | 0.2297 |
| EpiLPS (MALA) | 0.826 | −0.0079 | 0.0733 | 0.1976 |
| EpiLPS (MAP) | 0.753 | −0.0079 | 0.0720 | 0.1716 |

EpiEstim and EpiFilter land on nominal 95%. Both EpiLPS arms miss it *with a
correctly specified SI*, scored against a truth built from an EpiLPS fit — the
arrangement that should flatter them most. Some of the gap is the negative
binomial noise (caveat above), but note that the noise cuts both ways: it widens
the intervals too, and they still fall short.

### 3. The MALA/MAP gap is confined to the intervals

This is the one fully controlled comparison in the figure — identical simulated
data, identical point estimator (both take the posterior median `Rq0.50`):

```
max |Δ Bias|   0.00033
max |Δ RMSE|   0.00213
Δ Coverage     +0.043 to +0.078  (mean +0.069)
CI width ratio  1.145 to 1.161   (mean 1.152)
```

The point estimates are identical to four decimal places — which is why the MAP
and MALA curves coincide in the Bias and RMSE panels, and why MALA is drawn
as a long dash so the solid MAP curve shows through the gaps. The Laplace
approximation's lognormal interval is **~13% narrower** than the sampled
posterior, and buys back 6.9 coverage points when replaced by MALA. But MALA
still reaches only 0.826, so **the Laplace approximation explains roughly a third
of EpiLPS's undercoverage here; the rest is not a sampler artefact.**

### 4. EpiFilter is the most efficient estimator on this series

At the correct SI, EpiFilter matches EpiEstim's 0.950 coverage with **45% less
RMSE** (0.0626 vs 0.1135) and **56% of the interval width** (0.230 vs 0.410).
Both face Poisson data on the same real Rt trajectory, so this comparison is not
confounded by the observation model — only by the different plug-in truths.
EpiEstim's interval is carrying a lot of slack to reach nominal coverage.

### 5. Bias is small everywhere; this is a story about interval width

Largest |bias| anywhere on the grid is 0.0336 (EpiEstim, mean +60%), against a
true Rt that ranges 0.62–7.80. Every coverage collapse in the figure is driven by
interval width, not by point estimates drifting off the truth.

The bias sweeps are, however, **asymmetric in the mean arm**, and differently so:

Bias at the grid extremes, mean arm:

| method | at −60% | at +60% |
|---|---|---|
| EpiEstim | −0.0271 | **+0.0336** |
| EpiFilter | −0.0039 | **+0.0157** |
| EpiLPS (MAP) | −0.0080 | +0.0024 |
| EpiLPS (MALA) | −0.0082 | +0.0020 |

Overstating the serial interval biases EpiFilter upward about 4× harder than
understating it biases it down, and EpiEstim carries the largest bias at both
ends. Both EpiLPS arms are far flatter at the high end — their +60% bias is an
order of magnitude below EpiEstim's.

Coverage does **not** follow the bias asymmetry. In the mean arm EpiEstim runs
0.592 (−60%) against 0.669 (+60%) and EpiLPS 0.293 / 0.527 — the *understated*
end is worse for both, even though that is the end where bias is smaller.
EpiFilter goes the other way, 0.614 / 0.559. Bias and coverage are being driven
by different mechanisms here, which is why both panels are worth reading.

---

## Can EpiLPS fit Poisson?

No. `estimR()` and `estimRmcmc()` take no distribution argument; the likelihood
is hard-coded negative binomial with the overdispersion ρ estimated under a
Gamma(`a_rho`, `b_rho`) prior (defaults 1e-4, 1e-4). The package's *simulator*
`episim()` does expose `dist` and `overdisp`, so EpiLPS can generate Poisson data
but not fit it.

A negative binomial tends to a Poisson as ρ → ∞, so the obvious workaround is to
force ρ upward through the `priors` argument. Measured on Poisson data, it does
not work — ρ will not move, and the fit degenerates instead:

| prior (`a_rho`, `b_rho`) | prior mean of ρ | recovered ρ | Coverage | CI width |
|---|---|---|---|---|
| 1e-4, 1e-4 (default) | 1 | 25.2 | 0.604 | 0.0700 |
| 1e3, 1e-2 | 1e5 | 25.3 | 0.519 | 0.0632 |
| 1e6, 1e-6 | 1e12 | 26.0 | 0.072 | 0.0166 |

ρ stays pinned near 25 across twelve orders of magnitude of prior mean while the
interval collapses to 0.0166 — numerical breakdown of the Laplace step, not a
valid Poisson fit.

**Consequence for this comparison.** An apples-to-apples rerun with one shared
observation model is not available on the EpiLPS side. The options are to move
EpiEstim and EpiFilter onto negative-binomial data instead (both would need their
simulators changed and rerunning), or to keep the observation models as they are
and read the figure the way this README recommends — curve shape across methods,
absolute level only within the EpiLPS pair.

---

## Reproducing

```sh
cd queens_rt
Rscript combined/make_combined_panel.R
```

Reads the three `*_metrics.csv` files already in `results_*/` plus the seed CSV
by absolute path. To draw EpiFilter's filter arm instead of the smoother, change
`epifilter_arm` at the top of the script; to use a different assumed-SI family,
change `fam` (`gamma_discr`, `gamma_bin`, `lnorm`, `weibull`, `unif`).
