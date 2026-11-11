#!/bin/bash
# Addendum A1 pooling: v1.1 OOF rows unchanged + externals + Kalign 3.6.1/3.7.0 CLI rows as external arms.
cd ~/tmp/kalign-refit/presets_a1
PY=../venv-fit/bin/python
A=../runs/addendum_a1; R=../runs/v1_1
for p in protein nucleotide; do
  $PY -m recipe.pool_results --pool $p --out $A/refit --retime $R/refit/$p/retime/retime_rows.tsv \
      --externals $R/externals_$p.tsv $A/${p}_kalign361_*.tsv $A/${p}_kalign370_*.tsv
done
