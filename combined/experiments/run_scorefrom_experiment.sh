#!/bin/zsh
# Score from an earlier day and let na.rm handle whatever is missing.
#
# The committed scripts score from day 14 - the first day all three estimators
# have recovered from the single-day I0 = 2000 seed - so that every arm averages
# over the same set of days. This runner answers the complement: drop nothing,
# start at day 2, and let each estimator contribute whatever days it will.
#
# The na.rm = TRUE on the aggregations lives in the committed scripts, not here.
# It is a no-op at score_from = 14 (verified bit-identical) and only bites once
# the window moves earlier, so score_from is the single knob this runner turns.
#
# WHAT na.rm DOES AND DOES NOT REACH:
#   EpiEstim  - genuine NA, from its own `t_end > final_mean_si` rule. The
#               number of missing leading days slides with the assumed SI mean
#               (2 at mean 3.0, 10 at 12.0). Skipped correctly.
#   EpiFilter - no estimate before filter_start = 8, emitted as per-day NA.
#               Skipped correctly.
#   EpiLPS    - NOT reached. Its B-spline has no grid cap and returns a finite
#               number however little history it has: 848 on day 1, 244 on day 2
#               against a truth of 6.5. Those are wrong, not missing, and they
#               enter the average - roughly +1.08 of bias spread over 252 days,
#               against -0.006 for the day-14 window. That asymmetry is the
#               point of the comparison, not a bug to work around.
#
# Writes to combined/experiments/score-from-<N>/ so the day-14 results behind
# combined_si_misspec_panel1 are left intact.
#
# MALA is deliberately not run.
#
# Second argument is the plausibility ceiling r_max_plausible (default Inf, i.e.
# off). Setting it to 5 is what produces the committed panel2: it removes the
# absurd EpiLPS values that na.rm cannot reach, at the cost of also removing
# days where the truth legitimately exceeds 5. See the note in each script.
set -e
SF=${1:-2}
RMAX=${2:-Inf}
OUT="combined/experiments/score-from-$SF"
[ "$RMAX" != "Inf" ] && OUT="$OUT-rmax$RMAX"
TMP=$(mktemp -d)
mkdir -p "$OUT"

sed -e "s/^score_from   <- 14L.*/score_from   <- ${SF}L/" \
    -e "s|^out_dir    <- \"combined\"|out_dir    <- \"$OUT\"|" \
    combined/run_epiestim_arm.R | sed "s/^r_max_plausible <- Inf/r_max_plausible <- $RMAX/" > "$TMP/ee.R"

sed -e "s/^score_from <- 14L/score_from <- ${SF}L/" \
    -e "s|^output_dir <- \"results_epilps_si_misspec\"|output_dir <- \"$OUT/epilps_map\"|" \
    epilps_si_misspec.R | sed "s/^r_max_plausible <- Inf/r_max_plausible <- $RMAX/" > "$TMP/map.R"

sed -e "s/^score_from <- 14L/score_from <- ${SF}L/" \
    -e "s|^output_dir <- \"results_epifilter_si_misspec\"|output_dir <- \"$OUT/epifilter\"|" \
    epifilter_si_misspec.R | sed "s/^r_max_plausible <- Inf/r_max_plausible <- $RMAX/" > "$TMP/ef.R"

for f in ee map ef; do
  grep -qE "^score_from +<- +${SF}L" "$TMP/$f.R" || { echo "score_from patch failed in $f"; exit 1; }
  grep -qE "^r_max_plausible <- ${RMAX}\$" "$TMP/$f.R" || { echo "r_max patch failed in $f"; exit 1; }
done

echo "=== score_from = $SF, r_max_plausible = $RMAX -> $OUT ==="
Rscript "$TMP/ee.R"  > "$OUT/epiestim.log"   2>&1 && echo "  EpiEstim   done"
Rscript "$TMP/map.R" > "$OUT/epilps_map.log" 2>&1 && echo "  EpiLPS MAP done"
Rscript "$TMP/ef.R"  > "$OUT/epifilter.log"  2>&1 && echo "  EpiFilter  done"
rm -rf "$TMP"
