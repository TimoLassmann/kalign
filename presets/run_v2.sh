#!/bin/bash
# Amendment v2 (docs/PLAN.md): re-select two modes per pool at the knees of every v1.1 fit's own Pareto front,
# test them out of fold, record Freeze 2, score the final selection in sample, pool. Workstation only.
#   presets/run_v2.sh FREEZE_v2.json        (manifest of `python -m recipe.freeze record --rev <tag>`)
set -euo pipefail
FREEZE="$1"
WS=~/tmp/kalign-refit
PY="$WS/venv-fit/bin/python"
SRC="$WS/runs/v1_1"; OUT="$WS/runs/v2"
cd "$(dirname "$0")"
mkdir -p "$OUT"
log() { echo "$(date +%FT%T) $*" | tee -a "$OUT/progress.log"; }
"$PY" -m recipe.freeze check --manifest "$FREEZE" || { log "REFUSED: code differs from the frozen commit"; exit 4; }
cp "$FREEZE" "$OUT/freeze1.json"
for p in protein nucleotide; do
  log "START select $p"
  "$PY" -m recipe.select_v2 --pool $p --src "$SRC/refit" --dst "$OUT/refit" --externals "$SRC/externals_$p.tsv" \
    | tee -a "$OUT/progress.log"
done
for p in protein nucleotide; do
  log "START test $p"
  "$PY" -m recipe.nested test --pool $p --ws "$WS" --out "$OUT/refit" >> "$OUT/run.log" 2>&1
  log "DONE test $p"
done
"$PY" - "$OUT" <<'PYEOF'
import hashlib, json, sys, time
from pathlib import Path
out = Path(sys.argv[1])
h = {str(p.relative_to(out)): hashlib.sha256(p.read_bytes()).hexdigest()
     for p in sorted((out / "refit").glob("*/fold*/oof_rows.tsv"))}
(out / "freeze2.json").write_text(json.dumps({"recorded": time.strftime("%Y-%m-%dT%H:%M:%S"),
                                              "oof_rows_sha256": h}, indent=1) + "\n")
print("freeze 2:", len(h), "files")
PYEOF
log "FREEZE2 recorded"
for p in protein nucleotide; do
  log "START insample $p"
  "$PY" -m recipe.nested insample --pool $p --ws "$WS" --out "$OUT/refit" >> "$OUT/run.log" 2>&1
  log "START pool $p"
  "$PY" -m recipe.pool_results --pool $p --out "$OUT/refit" --externals "$SRC/externals_$p.tsv" \
      "$WS/runs/addendum_a1/${p}_kalign361_"*.tsv >> "$OUT/run.log" 2>&1
done
"$PY" - "$OUT" <<'PYEOF'
import hashlib, json, sys
from pathlib import Path
out = Path(sys.argv[1])
old = json.loads((out / "freeze2.json").read_text())["oof_rows_sha256"]
new = {str(p.relative_to(out)): hashlib.sha256(p.read_bytes()).hexdigest()
       for p in sorted((out / "refit").glob("*/fold*/oof_rows.tsv"))}
sys.exit(0 if old == new else "freeze 2 violated: out-of-fold rows changed")
PYEOF
log "RUN complete"
