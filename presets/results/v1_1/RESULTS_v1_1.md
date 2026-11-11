# Re-fit v1.1: results summary (recorded 2026-10-07)

Run: Freeze 1 tag `exp/refit-v1.1-freeze1` (commit 109d2fa), workstation, 2026-10-06 20:44 to 2026-10-07 12:25.
Amendment v1.1 (PLAN.md): identical to v1 except the stage-2 generation cap, 40 -> 140; every v1 fit (5 folds +
final, both pools) continued from its stage-2 state. Stage 1, folds, cores, selection and decision rules unchanged.
Full tables: `report/report.md`, `report/tables/`. Generated presets: `presets.json`. v1 for comparison: `../v1/`.

## Convergence (the reason for v1.1)

Stage 2 now stopped by the hypervolume rule in all 12 fits (protein generations 62-119, nucleotide 58-86); none
hit the 140 cap. v1 hit its 40-generation cap in all 12. Stage-2 hypervolume rose 1-4% over v1.

## Decision (unchanged from v1)

| pool | mode | v1.1 preset |
|---|---|---|
| protein | fast | 3.7.0 (rule not met) |
| protein | default | re-fit v1.1 |
| protein | recall = accurate | re-fit v1.1 |
| nucleotide | fast | 3.7.0 |
| nucleotide | default | 3.7.0 (SP -0.039, TC -0.042) |
| nucleotide | recall | re-fit v1.1 (same selection as v1) |
| nucleotide | accurate | 3.7.0 (SP -0.011, TC -0.012) |

## Protein, out of fold vs 3.7.0 (benchmark-weighted; deltas with 95% CI)

| mode | v1 dSP / dTC / dPrecision | v1.1 dSP / dTC / dPrecision | CPU ratio v1 / v1.1 |
|---|---|---|---|
| default | +0.003 / +0.015 / +0.045 | +0.002 / +0.015 / +0.056 | 1.02 / 1.02 |
| accurate | -0.005 / +0.007 / +0.043 | -0.002 [-0.006, +0.001] / +0.011 [+0.004, +0.019] / +0.043 | 1.05 / 1.16 |

Case-weighted, accurate: v1 dSP -0.010, dTC -0.006; v1.1 dSP -0.007, dTC -0.003. Default case-weighted dSP
-0.005 in both. Absolute (benchmark-weighted) accurate: SP 0.755, TC 0.590, precision 0.760, F1 0.751 (3.7.0:
0.757 / 0.579 / 0.717 / 0.733). External aligners unchanged: MUSCLE 5 0.796 / 0.639, MAFFT L-INS-i 0.784 / 0.627,
Clustal Omega 0.761 / 0.593.

## Overfitting (optimism = final-fit in-sample minus OOF)

Protein: SP and TC within +0.006 in every mode (as v1); precision +0.002 to +0.008. No meaningful overfitting.
Nucleotide: fast and recall as v1 (within +-0.003). Default +0.012 to +0.018 SP/TC, accurate -0.007 to -0.015:
the nucleotide stage-2 selections vary between folds and the final fit (flat front); neither mode ships.

## Reading

The longer search converged and moved accurate closer to the decision margin (dSP -0.005 -> -0.002, dTC +0.007 ->
+0.011) at 1.16x the 3.7.0 CPU (within the 1.2x budget); default is unchanged in SP/TC with higher precision.
The ship/keep decision is the same in every mode. v1 was mildly under-searched; the gain from searching further
is small, so the v1 conclusions stand.
