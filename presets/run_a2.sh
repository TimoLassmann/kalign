#!/bin/bash
# Addendum A2 (docs/PLAN.md): Kalign 3.1 and 3.3.5 (defaults, single thread) on every case of both pools, same
# harness as the external aligners; pooled with the v2 out-of-fold rows. Descriptive only. Workstation.
#   presets/run_a2.sh /path/to/kalign/git/archives   (v3.1.tar.gz and v3.3.5.tar.gz made with `git archive`)
set -euo pipefail
ARCH="$1"
WS=~/tmp/kalign-refit
PY="$WS/kalign/960a4f6/venv/bin/python"; PYF="$WS/venv-fit/bin/python"
OUT="$WS/runs/addendum_a2"; V2="$WS/runs/v2"; V11="$WS/runs/v1_1"
cd "$(dirname "$0")"
mkdir -p "$OUT"
log() { echo "$(date +%FT%T) $*" | tee -a "$OUT/progress.log"; }

# --- builds (x86_64 with the SIMD each release was written for)
K31="$WS/kalign/v3.1"; K335="$WS/kalign/v3.3.5"
if [ ! -x "$K31/kalign31" ]; then
  mkdir -p "$K31" && tar xzf "$ARCH/v3.1.tar.gz" -C "$K31" --strip-components=1
  ( cd "$K31/src"
    printf '#define HAVE_AVX2 1\n#define HAVE_GETPID 1\n#define VERSION "3.1"\n#define PACKAGE_VERSION "3.1"\n#define PACKAGE_NAME "kalign"\n#define PACKAGE_BUGREPORT ""\n' > config.h
    echo 'char build_config[] = "x86_64 AVX2 build of v3.1 (Addendum A2)";' > bc.c
    gcc -O2 -mavx2 -fcommon -DHAVE_CONFIG_H -I. -w -o ../kalign31 run_kalign.c tldevel.c rng.c parameters.c aln_data.c \
      alphabet.c alignment_parameters.c pick_anchor.c bisectingKmeans.c sequence_distance.c euclidean_dist.c \
      alignment.c weave_alignment.c rwalign.c misc.c bpm.c bc.c -lpthread -lm )
fi
if [ ! -x "$K335/build/src/kalign" ]; then
  mkdir -p "$K335" && tar xzf "$ARCH/v3.3.5.tar.gz" -C "$K335" --strip-components=1
  ( mkdir -p "$K335/build" && cd "$K335/build" && /usr/bin/cmake -DCMAKE_BUILD_TYPE=Release -DUSE_OPENMP=OFF \
      -DCMAKE_POLICY_VERSION_MINIMUM=3.5 .. > cmake.log 2>&1 && make -j16 > make.log 2>&1 )
fi
"$K31/kalign31" --version | head -2 | tee -a "$OUT/progress.log"
"$K335/build/src/kalign" --version | head -1 | tee -a "$OUT/progress.log"

# --- runs (8 workers, as Addendum A1)
for p in protein nucleotide; do
  for t in kalign31 kalign335; do
    log "START $p $t"
    "$PY" -m recipe.externals --tool $t --cases "$WS/data/cases_$p.final.tsv" --workers 8 \
        --cache "$OUT/cache_a2.sqlite" --out "$OUT/${p}_$t.tsv" >> "$OUT/run.log" 2>&1
    log "DONE $p $t"
  done
done

# --- pooling with the v2 out-of-fold rows (fold directories linked read-only)
for p in protein nucleotide; do
  mkdir -p "$OUT/refit/$p"
  for d in fold0 fold1 fold2 fold3 fold4 final; do ln -sfn "$V2/refit/$p/$d" "$OUT/refit/$p/$d"; done
  "$PYF" -m recipe.pool_results --pool $p --out "$OUT/refit" --externals "$V11/externals_$p.tsv" \
      "$WS/runs/addendum_a1/${p}_kalign361_"*.tsv "$OUT/${p}_kalign31.tsv" "$OUT/${p}_kalign335.tsv" >> "$OUT/run.log" 2>&1
  log "pooled $p"
done
log "RUN complete"
