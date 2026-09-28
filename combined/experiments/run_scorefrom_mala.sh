#!/bin/zsh
# The EpiLPS-MALA arm of the score-from experiment, run on its own.
#
# run_scorefrom_experiment.sh deliberately skips MALA: a MALA fit costs ~300x a
# MAP fit, and the other three arms finish in minutes. This runner applies the
# same score_from / r_max_plausible patches to epilps_si_misspec_mala.R and
# writes into the same combined/experiments/score-from-<N>[-rmax<R>]/ folder, as
# epilps_mala/, next to the epilps_map/ it is paired with. Run the main runner
# first with the same arguments: the MALA script reads that epilps_map/ back for
# its MALA-vs-MAP comparison (map_metrics_file is patched to point there).
#
# Third argument is the assumed-SI families to fit, space-separated (default
# "gamma_bin", the only family combined_si_misspec_panel2 draws). All five
# families is ~5x the cost. The simulated epidemics are bit-identical either
# way; the MALA chains for a given family are not, because each replicate's
# worker runs one L'Ecuyer stream through every setting in turn, so the draws a
# gamma_bin fit sees depend on how many fits came before it. Same distribution,
# different Monte Carlo realisation.
#
# Fourth and fifth arguments are mcmc_niter and mcmc_burnin (default: the
# script's own 5000 / 2000). The script is right that these are not a cheap
# runtime knob, but 5000 / 2000 is not converged on this design either:
# measured on 8 replicates at the correct SI, the median effective sample size of
# the Rt draws is ~17 of 3000, and the interval comes out narrower than LPSMAP's
# (width 0.205 vs 0.218). 20000 / 8000 gives ESS ~58 and width 0.218, and
# 40000 / 16000 gives ESS ~110 and width 0.224. A non-default length writes to
# epilps_mala_niter<N>/ so runs at different lengths sit side by side.
#
#   ./combined/experiments/run_scorefrom_mala.sh 2 5            # panel2's arm
#   ./combined/experiments/run_scorefrom_mala.sh 2 5 "gamma_discr gamma_bin"
#   ./combined/experiments/run_scorefrom_mala.sh 2 5 gamma_bin 20000 8000
set -e
SF=${1:-2}
RMAX=${2:-Inf}
FAMS=${3:-gamma_bin}
NITER=${4:-5000}
BURNIN=${5:-2000}
OUT="combined/experiments/score-from-$SF"
[ "$RMAX" != "Inf" ] && OUT="$OUT-rmax$RMAX"
MALA_DIR="epilps_mala"
[ "$NITER" != "5000" -o "$BURNIN" != "2000" ] && MALA_DIR="epilps_mala_niter$NITER"
TMP=$(mktemp -d)
mkdir -p "$OUT"

[ -f "$OUT/epilps_map/si_misspec_metrics.csv" ] || \
  echo "WARNING: $OUT/epilps_map/ not found - the MALA-vs-MAP comparison will be skipped"

FAM_R="c(\"${FAMS// /\", \"}\")"

sed -e "s/^score_from <- 14L\$/score_from <- ${SF}L/" \
    -e "s/^r_max_plausible <- Inf\$/r_max_plausible <- $RMAX/" \
    -e "s|^output_dir <- \"results_epilps_si_misspec_mala\"\$|output_dir <- \"$OUT/$MALA_DIR\"|" \
    -e "s|^map_metrics_file <- file.path(\"results_epilps_si_misspec\", |map_metrics_file <- file.path(\"$OUT/epilps_map\", |" \
    -e "s/^fit_families <- c(.*)\$/fit_families <- $FAM_R/" \
    -e "s/^mcmc_niter  <- 5000L .*/mcmc_niter  <- ${NITER}L/" \
    -e "s/^mcmc_burnin <- 2000L .*/mcmc_burnin <- ${BURNIN}L/" \
    epilps_si_misspec_mala.R > "$TMP/mala.R"

grep -qE "^score_from <- ${SF}L\$"               "$TMP/mala.R" || { echo "score_from patch failed"; exit 1; }
grep -qE "^r_max_plausible <- ${RMAX}\$"         "$TMP/mala.R" || { echo "r_max patch failed"; exit 1; }
grep -qF "output_dir <- \"$OUT/$MALA_DIR\""      "$TMP/mala.R" || { echo "output_dir patch failed"; exit 1; }
grep -qF "file.path(\"$OUT/epilps_map\","        "$TMP/mala.R" || { echo "map_metrics_file patch failed"; exit 1; }
grep -qF "fit_families <- $FAM_R"                "$TMP/mala.R" || { echo "fit_families patch failed"; exit 1; }
grep -qE "^mcmc_niter  <- ${NITER}L\$"           "$TMP/mala.R" || { echo "mcmc_niter patch failed"; exit 1; }
grep -qE "^mcmc_burnin <- ${BURNIN}L\$"          "$TMP/mala.R" || { echo "mcmc_burnin patch failed"; exit 1; }

echo "=== MALA: score_from = $SF, r_max_plausible = $RMAX, families = $FAMS, niter = $NITER / $BURNIN -> $OUT/$MALA_DIR ==="
Rscript "$TMP/mala.R" > "$OUT/$MALA_DIR.log" 2>&1 && echo "  EpiLPS MALA done"
rm -rf "$TMP"
