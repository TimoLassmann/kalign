#!/bin/sh
# Unattended start of the nucfit fits: wait for prepare -> data check -> identity checks -> fast-fit smoke ->
# launch the six fit pipelines (folds 0-4 + all). Any failed check stops before the fits (logged).
# The accurate confidence check failing only holds step 2 (runs/nucconf_hold), not the fast fits.
#   setsid -f sh scripts/nucfit_driver.sh      (from WS/presets)
set -u
WS=$HOME/tmp/kalign-nucfit
PY=$HOME/tmp/kalign-refit/venv-fit/bin/python
L=$WS/runs/driver.log
export TMPDIR=$WS/tmp KALIGN_REFIT_TMP=$WS/tmp OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "$WS/runs"
cd "$WS/presets" || exit 1
say() { echo "$(date +%FT%T) $*" >> "$L"; }
say "driver start"
while ! grep -q '^exit' "$WS/prepare.log"; do sleep 60; done
grep -q '^exit 0' "$WS/prepare.log" || { say "STOP: prepare failed"; exit 1; }
say "prepare done"
$PY -u -m recipe.nucfit_data check --ws "$WS" >> "$L" 2>&1 || { say "STOP: data check failed"; exit 1; }
say "data check passed"
$PY -u -m recipe.nucfit_identity --ws "$WS" --workers 30 >> "$L" 2>&1 || { say "STOP: identity checks failed"; exit 1; }
say "identity checks passed"
$PY -u -m recipe.nucconf check --ws "$WS" >> "$L" 2>&1 || { touch "$WS/runs/nucconf_hold"; say "confidence check failed: step 2 held"; }
$PY -u -m recipe.nucfast --fold 0 --out "$WS/runs/smoke/nucfast" --limit 40 --pop-size 8 --gens 1,1 >> "$L" 2>&1 \
    || { say "STOP: fast-fit smoke failed"; exit 1; }
say "smoke passed; launching six pipelines"
for f in 0 1 2 3 4 all; do
    setsid -f sh scripts/nucfit_pipeline.sh "$f" < /dev/null > /dev/null 2>&1
done
say "launched"
