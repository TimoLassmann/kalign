# Re-fit v1: results summary (recorded 2026-10-07)

Run: Freeze 1 tag `exp/refit-v1-freeze1` (commit dbb8744), workstation, 2026-10-05 12:58 to 2026-10-06 20:31.
Freeze 2 (out-of-fold rows hashed before the final fits): 2026-10-06 12:12. Results tag: `exp/refit-v1-results`.
Full tables: `report/report.md`, `report/tables/`. Generated presets: `presets.json`.

## Procedure as run (order)

1. Nested CV, 5 grouped folds per pool: fit on 4/5 (two-stage NSGA-III), select presets by the fixed rule,
   score them on the held-out 1/5. Pooled out-of-fold (OOF) scores = estimate of the procedure's performance.
2. Final fit: the identical procedure on all cases -> the presets in `presets.json`.
3. Optimism (final-fit in-sample minus OOF), external aligners, CPU re-timing at 4 workers, decision rule.

## Presets after the decision rule

| pool | mode | preset |
|---|---|---|
| protein | fast | 3.7.0 (rule not met) |
| protein | default | re-fit v1 |
| protein | recall = accurate | re-fit v1 |
| nucleotide | fast | 3.7.0 (rule not met) |
| nucleotide | default | 3.7.0 (rule not met: SP -0.034, TC -0.035) |
| nucleotide | recall | re-fit v1 |
| nucleotide | accurate | 3.7.0 (rule not met: SP -0.012, TC -0.013) |

Parameters: `presets.json`. Not compiled into Kalign.

## Out-of-fold vs 3.7.0 (protein, benchmark-weighted)

| mode | 3.7.0 SP / TC / precision | re-fit v1 SP / TC / precision |
|---|---|---|
| fast | 0.740 / 0.561 / 0.682 | 0.741 / 0.562 / 0.680 |
| default | 0.738 / 0.548 / 0.685 | 0.741 / 0.563 / 0.730 |
| accurate | 0.757 / 0.579 / 0.717 | 0.752 / 0.586 / 0.760 |

Gains are mainly precision (+0.04) and benchmark-weighted TC; case-weighted, accurate loses SP -0.010 and TC
-0.006; PREFAB and the most divergent cases (< 20% identity) lose SP (-0.01 to -0.03).

## Overfitting (optimism = final-fit in-sample minus OOF)

Protein: SP and TC +0.000 to +0.006 in every mode (case and benchmark weighting); precision +0.015 for default,
-0.008 for accurate. Nucleotide: within about +-0.008. Conclusion: no meaningful overfitting; the procedure's
in-sample scores are close to its out-of-fold scores. Protein gains hold for test cases with and without a
>= 50% identity homolog in training.

## Underfitting / convergence

Stage 1 (single-run presets) stopped by the hypervolume rule in 4/5 protein folds (fold 3 and the final fit
hit 100 generations) and in all nucleotide fits. Stage 2 (ensembles: default, recall, accurate) hit the
40-generation cap in EVERY fit (10 folds + 2 final fits) while still improving: the ensemble search was
compute-limited, not converged. Fold selections differ in structure but have similar scores (flat landscape).
The precision tilt follows from the pre-stated selection rule (max mean z of SP, TC, precision).

## Bias of the comparison

3.7.0 was not re-fitted per fold and was partly tuned on BAliBASE, so the comparison slightly favours 3.7.0.
