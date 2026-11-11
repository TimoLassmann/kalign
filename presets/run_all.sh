#!/usr/bin/env bash
# One command for the whole pooled nested-CV re-fit (Freeze 1 procedure), on the workstation.
#
#   presets/run_all.sh start  --freeze FREEZE1.json        detached, resumable, returns at once
#   presets/run_all.sh run    --freeze FREEZE1.json        same, in the foreground
#   presets/run_all.sh status                              progress log and running processes
#   presets/run_all.sh start  --smoke                      tiny end-to-end smoke test (separate runs dir/caches)
#
# Options: --ws DIR (default $KALIGN_REFIT_WS or ~/tmp/kalign-refit)   --out DIR (default WS/runs/full, smoke:
#          WS/runs/smoke_full)   --stop-after STEP   --freeze FILE (manifest of `python -m recipe.freeze record`)
#
# Steps (each writes OUT/steps/<name>.done and one line to OUT/progress.log; a restart skips finished steps,
# and every step is itself resumable from its checkpoints and from the evaluation cache):
#   1 cv_protein, cv_nucleotide  nested CV: 5 fits per pool (run together, they share the 30 aligner slots),
#                                then the out-of-fold test of each fold on ALL its cases
#   2 pool_oof                   pooled out-of-fold tables. FREEZE 2: sha256 of every oof_rows.tsv is recorded
#                                (OUT/freeze2.json) BEFORE any final fit starts; --stop-after pool_oof stops here
#                                so that the results can be committed and tagged
#   3 final_protein, final_nucleotide   the identical fit on the whole pool + in-sample scores
#   4 externals                  MAFFT L-INS-i, MUSCLE 5, Clustal Omega on all cases (comparison only)
#   5 retime                     4-worker re-timing of finalists and 3.7.0 presets (exclusive use of the machine)
#   6 pool_full                  pooled tables with cpu4, externals and optimism; verifies freeze 2 is intact
#   7 decide                     decision per mode (rule of Freeze 1) and presets/out/presets.json
#   8 report                     OUT/report/report.md and tables
#
# Refuses to start (and to resume) if presets/ differs from the frozen commit (recipe.freeze check), or if the
# data tables changed since the run began. Without --freeze only --smoke runs are possible.
set -euo pipefail

PRESETS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WS="${KALIGN_REFIT_WS:-$HOME/tmp/kalign-refit}"
CMD="${1:-}"; shift || true
FREEZE=""; OUT=""; SMOKE=0; STOP_AFTER=""
while [ $# -gt 0 ]; do
  case "$1" in
    --freeze) FREEZE="$2"; shift 2;;
    --ws) WS="$2"; shift 2;;
    --out) OUT="$2"; shift 2;;
    --smoke) SMOKE=1; shift;;
    --stop-after) STOP_AFTER="$2"; shift 2;;
    *) echo "unknown option: $1" >&2; exit 2;;
  esac
done
export KALIGN_REFIT_WS="$WS"
if [ -z "$OUT" ]; then
  if [ "$SMOKE" = 1 ]; then OUT="$WS/runs/smoke_full"; else OUT="$WS/runs/full"; fi
fi
PY="$WS/venv-fit/bin/python"                    # engine, pooling, decision, report (pymoo, numpy, scipy)
PYA="$WS/kalign/960a4f6/venv/bin/python"        # externals (same wheel; tools are called as subprocesses)
LOG="$OUT/progress.log"

# ---- settings: full run = registered settings (no overrides); smoke = tiny, everything separate
NEST_EXTRA=(); FIT_ARGS=""; DEC_EXTRA=(); POOL_EXTRA=(); EXT_EXTRA=(); DEC_SMOKE=(); FOLDS="0,1,2,3,4"
if [ "$SMOKE" = 1 ]; then
  FOLDS="0,1"
  FIT_ARGS="--pop-size 12 --ref-partitions 2 --stage1-gens 2,2 --stage2-gens 1,1 --limit-core 60"
  NEST_EXTRA=(--allow-smoke --cache "$OUT/cache.sqlite" --test-limit 60 --fit-args "$FIT_ARGS")
  DEC_EXTRA=(--allow-smoke --retime-cache "$OUT/retime4.sqlite" --test-limit 60)
  POOL_EXTRA=(--allow-smoke --allow-partial); DEC_SMOKE=(--allow-smoke)
  EXT_EXTRA=(--cache "$OUT/cache_externals.sqlite" --limit-per-pool 30)
fi

stamp() { date +%Y-%m-%dT%H:%M:%S; }
plog() { echo "$(stamp) $*" | tee -a "$LOG"; }

usage() { sed -n '2,32p' "${BASH_SOURCE[0]}"; }

case "$CMD" in
  status)
    echo "== $LOG"; cat "$LOG" 2>/dev/null || echo "(no progress log)"
    echo "== processes"; pgrep -af "recipe\.|run_all.sh" || echo "(none running)"
    exit 0;;
  start)
    mkdir -p "$OUT"
    args=(run --ws "$WS" --out "$OUT"); [ -n "$FREEZE" ] && args+=(--freeze "$FREEZE")
    [ "$SMOKE" = 1 ] && args+=(--smoke); [ -n "$STOP_AFTER" ] && args+=(--stop-after "$STOP_AFTER")
    setsid nohup "$PRESETS/run_all.sh" "${args[@]}" > "$OUT/run_all.out" 2>&1 < /dev/null &
    echo "started (pid $!); progress: $LOG ; output: $OUT/run_all.out"
    exit 0;;
  run) ;;
  *) usage; exit 2;;
esac

# ------------------------------------------------------------------------------------------ run
mkdir -p "$OUT/steps" "$OUT/logs"
exec 9> "$OUT/run.lock"
flock -n 9 || { echo "another run_all.sh is active on $OUT" >&2; exit 3; }
cd "$PRESETS"

# ---- guards: frozen code, frozen data
if [ -n "$FREEZE" ]; then
  [ -f "$FREEZE" ] || { echo "freeze manifest not found: $FREEZE" >&2; exit 4; }
  if [ -f "$OUT/freeze1.json" ] && ! cmp -s "$FREEZE" "$OUT/freeze1.json"; then
    echo "REFUSED: $OUT was started under a different freeze manifest" >&2; exit 4
  fi
  "$PY" -m recipe.freeze check --manifest "$FREEZE" || { echo "REFUSED: code differs from the frozen commit" >&2; exit 4; }
  [ -f "$OUT/freeze1.json" ] || cp "$FREEZE" "$OUT/freeze1.json"
elif [ "$SMOKE" = 1 ]; then
  echo "smoke run without a freeze manifest: code is NOT verified against a commit" >&2
else
  echo "REFUSED: a real run needs --freeze FILE (python -m recipe.freeze record --rev <tagged commit> --out FILE)" >&2
  exit 4
fi
if [ -f "$OUT/data_hashes.json" ]; then
  "$PY" -m recipe.freeze hash-data --ws "$WS" --verify "$OUT/data_hashes.json" || { echo "REFUSED: data changed" >&2; exit 4; }
else
  "$PY" -m recipe.freeze hash-data --ws "$WS" --out "$OUT/data_hashes.json"
fi
COMMIT="$("$PY" -c "import json,sys;m=json.load(open(sys.argv[1]));print(m['commit'][:12], m.get('tag') or '')" "${FREEZE:-/dev/null}" 2>/dev/null || echo unfrozen)"
plog "RUN  start (commit/tag: $COMMIT; smoke=$SMOKE; out=$OUT; $(nproc) cores)"

# ---- step runner
FAILED=0
done_p() { [ -f "$OUT/steps/$1.done" ]; }
run_step() {          # run_step NAME cmd...   (foreground; log to OUT/logs/NAME.log)
  local name="$1"; shift
  if done_p "$name"; then plog "SKIP $name (already done)"; return 0; fi
  plog "START $name"
  local t0=$SECONDS
  if "$@" >> "$OUT/logs/$name.log" 2>&1; then
    touch "$OUT/steps/$name.done"; plog "DONE $name elapsed=$((SECONDS - t0))s"
  else
    plog "FAIL $name rc=$? (see $OUT/logs/$name.log)"; FAILED=1; return 1
  fi
}
parallel2() {         # parallel2 NAME_A "cmd a" NAME_B "cmd b": both run, both are waited for
  local na="$1" ca="$2" nb="$3" cb="$4" pa pb ra=0 rb=0
  run_step "$na" bash -c "$ca" & pa=$!
  run_step "$nb" bash -c "$cb" & pb=$!
  wait "$pa" || ra=$?; wait "$pb" || rb=$?
  [ "$ra" = 0 ] && [ "$rb" = 0 ]
}
stop_here() { if [ "$STOP_AFTER" = "$1" ]; then plog "STOP after $1 (--stop-after)"; exit 0; fi; }
need() { [ "$FAILED" = 0 ] || { plog "RUN  aborted: a step failed"; exit 1; }; }

NESTED="$PY -m recipe.nested"
Q() { printf '%q ' "$@"; }          # shell-quote the arguments (never called with none)
NEST_COMMON=(--out "$OUT/refit" "${NEST_EXTRA[@]}")

# ---- 1  nested CV
parallel2 cv_protein "$NESTED cv --pool protein --folds $FOLDS $(Q "${NEST_COMMON[@]}")" \
          cv_nucleotide "$NESTED cv --pool nucleotide --folds $FOLDS $(Q "${NEST_COMMON[@]}")" || FAILED=1
need; stop_here cv

# ---- 2  pooled out-of-fold results + Freeze 2 record
pool_oof() {
  for p in protein nucleotide; do
    "$PY" -m recipe.pool_results --pool $p --out "$OUT/refit" "${POOL_EXTRA[@]}" || return 1
  done
  "$PY" - "$OUT" <<'EOF'
import hashlib, json, sys, time
from pathlib import Path
out = Path(sys.argv[1])
f = {str(p.relative_to(out)): hashlib.sha256(p.read_bytes()).hexdigest()
     for p in sorted((out / "refit").glob("*/fold*/oof_rows.tsv"))}
if not f:
    raise SystemExit("no oof_rows.tsv found")
fz = out / "freeze2.json"
if fz.exists():
    old = json.loads(fz.read_text())["oof_rows_sha256"]
    if old != f:
        raise SystemExit("FREEZE 2 VIOLATED: out-of-fold rows changed after they were recorded")
else:
    fz.write_text(json.dumps({"recorded": time.strftime("%Y-%m-%dT%H:%M:%S"), "oof_rows_sha256": f}, indent=1) + "\n")
print("freeze 2 recorded/verified:", len(f), "oof_rows files")
EOF
  [ $? = 0 ] || return 1
  for p in protein nucleotide; do
    rm -rf "$OUT/refit/$p/results_freeze2"; cp -r "$OUT/refit/$p/results" "$OUT/refit/$p/results_freeze2"
  done
}
run_step pool_oof pool_oof || true
need
[ -f "$OUT/freeze2.json" ] && plog "FREEZE2 recorded: $(sha256sum "$OUT/freeze2.json" | cut -c1-16) ($OUT/freeze2.json)"
stop_here pool_oof

# ---- 3  final fits
parallel2 final_protein "$NESTED final --pool protein $(Q "${NEST_COMMON[@]}")" \
          final_nucleotide "$NESTED final --pool nucleotide $(Q "${NEST_COMMON[@]}")" || FAILED=1
need; stop_here final

# ---- 4  external aligners (comparison only), all cases of each pool, 30 workers
externals() {
  for p in protein nucleotide; do
    "$PYA" -m recipe.externals --tool all --cases "$WS/data/cases_$p.final.tsv" --workers 30 \
        --out "$OUT/externals_$p.tsv" "${EXT_EXTRA[@]}" || return 1
  done
}
run_step externals externals || true
need; stop_here externals

# ---- 5  4-worker re-timing (exclusive)
retime() {
  for p in protein nucleotide; do
    "$PY" -m recipe.decide retime --pool $p --out "$OUT/refit" --folds "$FOLDS" "${DEC_EXTRA[@]}" || return 1
  done
}
run_step retime retime || true
need; stop_here retime

# ---- 6  pooled tables with cpu4, externals, optimism
pool_full() {
  for p in protein nucleotide; do
    rt="$OUT/refit/$p/retime/retime_rows.tsv"
    "$PY" -m recipe.pool_results --pool $p --out "$OUT/refit" --retime "$rt" \
        --externals "$OUT/externals_$p.tsv" "${POOL_EXTRA[@]}" || return 1
  done
  "$PY" - "$OUT" <<'EOF'
import hashlib, json, sys
from pathlib import Path
out = Path(sys.argv[1])
old = json.loads((out / "freeze2.json").read_text())["oof_rows_sha256"]
new = {str(p.relative_to(out)): hashlib.sha256(p.read_bytes()).hexdigest()
       for p in sorted((out / "refit").glob("*/fold*/oof_rows.tsv"))}
if old != new:
    raise SystemExit("FREEZE 2 VIOLATED: out-of-fold rows differ from the recorded ones")
print("freeze 2 intact:", len(new), "oof_rows files unchanged since recorded")
EOF
}
run_step pool_full pool_full || true
need; stop_here pool_full

# ---- 7  decision + preset table
decide() {
  for p in protein nucleotide; do
    "$PY" -m recipe.decide decide --pool $p --out "$OUT/refit" "${DEC_SMOKE[@]}" || return 1
  done
  "$PY" -m recipe.decide emit --out "$OUT/refit" --dest "$OUT/presets.json"
}
run_step decide decide || true
need; stop_here decide

# ---- 8  report
run_step report "$PY" -m recipe.report --out "$OUT" || true
need
plog "RUN  complete: report $OUT/report/report.md ; presets $OUT/presets.json"
