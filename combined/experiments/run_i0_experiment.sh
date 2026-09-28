#!/bin/zsh
# Controlled I0 experiment: run all three arms at a given I0, everything else
# held fixed. The scripts are copied and patched rather than parameterised in
# place so the committed versions keep their own settings and the two arms of
# the experiment cannot drift.
#
# filter_start is forced to 8 for BOTH I0 values. At I0 = 2000 the default 5
# fails on 100/100 replicates at assumed mean +60% (day 3-6 have L ~ 0 against
# a long assumed SI, the implied Rt reaches 8794, and the Poisson likelihood
# underflows on every grid point). 8 starts the recursion past that window and
# fails 0/100. Using a different filter_start per arm would confound the
# comparison, so both use 8.
set -e
I0=$1
OUT="combined/experiments/I0-$I0"
TMP=$(mktemp -d)
mkdir -p "$OUT"

sed -e "s/^I0           <- 2000L/I0           <- ${I0}L/" \
    -e "s|out_dir    <- \"combined\"|out_dir    <- \"$OUT\"|" \
    combined/run_epiestim_arm.R > "$TMP/ee.R"

sed -e "s/^I0         <- 2000L/I0         <- ${I0}L/" \
    -e "s|^output_dir <- \"results_epilps_si_misspec\"|output_dir <- \"$OUT/epilps_map\"|" \
    epilps_si_misspec.R > "$TMP/map.R"

sed -e "s/^I0       <- 2000L/I0       <- ${I0}L/" \
    -e "s/^filter_start <- 5L/filter_start <- 8L/" \
    -e "s|^output_dir <- \"results_epifilter_si_misspec\"|output_dir <- \"$OUT/epifilter\"|" \
    epifilter_si_misspec.R > "$TMP/ef.R"

for f in ee map ef; do
  grep -qE "^I0 +<- +${I0}L" "$TMP/$f.R" || { echo "I0 patch failed in $f"; exit 1; }
done
grep -q "^filter_start <- 8L" "$TMP/ef.R" || { echo "filter_start patch failed"; exit 1; }

echo "=== I0 = $I0 -> $OUT ==="
Rscript "$TMP/ee.R"  > "$OUT/epiestim.log"  2>&1 && echo "  EpiEstim  done"
Rscript "$TMP/map.R" > "$OUT/epilps_map.log" 2>&1 && echo "  EpiLPS MAP done"
Rscript "$TMP/ef.R"  > "$OUT/epifilter.log" 2>&1 && echo "  EpiFilter done"
rm -rf "$TMP"
