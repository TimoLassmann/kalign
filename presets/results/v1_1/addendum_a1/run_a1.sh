#!/bin/bash
# Addendum A1: Kalign 3.6.1 and 3.7.0 CLI on every case, same harness, 8 workers, alternating builds.
cd ~/tmp/kalign-refit/presets_a1
PY=../kalign/960a4f6/venv/bin/python
OUT=../runs/addendum_a1
for p in protein nucleotide; do
  for m in fast default recall accurate; do
    for b in kalign361 kalign370; do
      echo "$(date +%FT%T) START $p $b $m" >> $OUT/progress.log
      $PY -m recipe.externals --tool ${b}_$m --cases ../data/cases_$p.final.tsv --workers 8 \
          --cache $OUT/cache_a1.sqlite --out $OUT/${p}_${b}_$m.tsv >> $OUT/run.log 2>&1 \
        && echo "$(date +%FT%T) DONE $p $b $m" >> $OUT/progress.log \
        || echo "$(date +%FT%T) FAIL $p $b $m" >> $OUT/progress.log
    done
  done
done
echo "$(date +%FT%T) ALL DONE" >> $OUT/progress.log
