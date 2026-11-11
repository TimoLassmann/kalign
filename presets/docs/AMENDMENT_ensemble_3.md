# Amendment 3 to the ensemble plan (PLAN_ensemble.md), 2026-10-16

Requested by the user on 2026-10-16. The frozen plan is not edited; this file changes it on top of
`AMENDMENT_ensemble_1.md` and `AMENDMENT_ensemble_2.md`.

## Change

Change (b) of amendment 2 (only fold 0 and the final fit) is reversed: folds 1-4 are run again, so the out-of-fold
evaluation uses all five folds as in the frozen plan (section 7, criteria 1-5; the lockbox stays unusable, as stated
in amendment 1).

Everything else is as amendment 2. Folds 1-4 run with exactly the code, build and settings of the running fold-0 and
final fits: every `run.json` records ensfit sha256 620c6797..., build 534acefb..., generations 20-30, population 56,
seed 1234, amendments 1 + 2. No result of any fit was used for this change; the reason is completeness of the
out-of-fold evaluation for the paper.

## Record

- Fold 0 and final: running since 2026-10-15 21:46:39 (WS time), unchanged.
- Folds 1-4: started 2026-10-16 13:49:33 (WS time), fresh run directories `runs/ensfit/fold1..4`; their amendment-1
  runs (stopped at generation 2) are in `runs/ensfit/_aborted_amend1/`. Alignments already in the cache are reused.
- Launch log: `runs/ensfit/launch_folds1-4.out`. The six fits share the machine-wide 60-slot pool.
