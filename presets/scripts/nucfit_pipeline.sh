#!/bin/sh
# One fit of the nucfit plan, end to end: fast fit -> accurate step 1 (ensemble) -> accurate step 2 (confidence).
#   sh nucfit_pipeline.sh FOLD          (0..4 or all; run from WS/presets)
# Every step is resumable; a step whose selection.json exists is skipped. Logs in WS/runs/pipeline/<fit>.log.
set -u
FOLD=$1
WS=$HOME/tmp/kalign-nucfit
PY=$HOME/tmp/kalign-refit/venv-fit/bin/python
NAME=$( [ "$FOLD" = all ] && echo final || echo "fold$FOLD" )
export KALIGN_SCORE_PYTHON=$WS/kalign/venv-score/bin/python TMPDIR=$WS/tmp KALIGN_REFIT_TMP=$WS/tmp OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "$WS/runs/pipeline"
LOG=$WS/runs/pipeline/$NAME.log
cd "$WS/presets" || exit 1
step() {   # step NAME DIR CMD...
    n=$1; d=$2; shift 2
    if [ -f "$d/selection.json" ]; then echo "$(date +%FT%T) $n: done before, skipped" >> "$LOG"; return 0; fi
    echo "$(date +%FT%T) $n: start" >> "$LOG"
    "$@" >> "$LOG" 2>&1
    rc=$?
    echo "$(date +%FT%T) $n: exit $rc" >> "$LOG"
    return $rc
}
step nucfast "$WS/runs/nucfast/$NAME" $PY -u -m recipe.nucfast --fold "$FOLD" --out "$WS/runs/nucfast/$NAME" || exit 1
step nucens "$WS/runs/nucens/$NAME" $PY -u -m recipe.nucens --fold "$FOLD" --out "$WS/runs/nucens/$NAME" || exit 1
while [ -f "$WS/runs/nucconf_hold" ]; do
    echo "$(date +%FT%T) nucconf held (confidence check failed; remove runs/nucconf_hold to continue)" >> "$LOG"
    sleep 3600
done
step nucconf "$WS/runs/nucconf/$NAME" $PY -u -m recipe.nucconf fit --fits "$FOLD" || exit 1
echo "$(date +%FT%T) pipeline $NAME complete" >> "$LOG"
