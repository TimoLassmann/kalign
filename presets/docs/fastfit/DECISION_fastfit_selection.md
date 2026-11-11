# Decision: fast-fit selection by the knee rule, one fast mode (2026-10-14)

Status: decided by the user on 2026-10-14, while the fits were still running (final fit at generation 61), before
any official selection, re-scoring on the core data or out-of-fold scoring had run. No held-out data had been seen.
This changes section 6 of `PLAN_fastfit.md` (frozen at tag `fastfit-freeze1`) after launch and is reported in the
paper as a dated deviation from the plan, with the reason below.

## Why the registered rule is replaced

The registered fast rule (section 6) is: feasible, Q_small >= Q_small,ref - 0.002, then the highest Q_large. It
never uses CPU, and a "cut-off, then maximise one objective" rule always lands at the end of the front. On the
final fit's screening data (generation 61) it chose a point at R = 1.45 whose Q_large (0.690) is about 0.005 above
points at R = 0.55: nearly 3x the CPU for a gain inside the noise. It also sat on the Q_small cut-off, with small SP
below 4.0.0 fast (0.770 vs 0.776), which conflicts with gate criterion 1 (small SP only). The front has a clear knee
near R = 0.5-0.6; the rule cannot find it.

## The new rule (fast pick)

Applied identically in every fold and in the final fit, on the core data, to the re-scored front (section 5):

1. Feasible, small SP >= small SP of 4.0.0 fast, and small TC >= small TC of 4.0.0 fast.
2. Among these, Q_large,max = the highest Q_large.
3. Keep the points with Q_large >= Q_large,max - delta.
4. Take the lowest R. Ties: higher Q_large, then key. No eligible point: the gm1 fast preset.

delta (defined 2026-10-14 ~20:15, committed before any knee pick was computed on the core data of any fit, and
before delta itself was computed): "within delta of the best" means statistically indistinguishable from the best.
For each fit, let E be the re-scored candidates that pass step 1 on the core data, and b the one with the highest
core Q_large. For each c in E, SE_c is the bootstrap standard error of Q_large(b) - Q_large(c) on that fit's core
data: extHomFam families resampled with replacement (a family brings all its 1,000 and 10,000 cases), Q_large
recomputed per resample exactly as in `summarise()` (mean over k of the per-k case means), paired (same resample for
b and c), 2,000 resamples, seed 1234. Step 3 keeps c when Q_large(b) - Q_large(c) <= 1.96 x SE_c (b is always kept).
Reported alongside as a sensitivity analysis, not used for the decision: the picks at fixed delta = 0.005 and 0.01.

Out-of-fold scoring evaluates both the registered rule's picks and the knee rule's picks on the same held-out cases,
so the effect of the change is visible. The gate (section 7) is unchanged and is applied to the knee rule's final
pick.

## Disclosure: what had been seen when this rule was chosen

The rule was chosen after inspecting the final fit's screening front (its own training data, generation 61),
including the point the rule selects there and the picks at delta = 0.005 and 0.01; the user had named a preferred
point on that front (the interim candidate below). Screening fronts of folds 0-4 (training data only) had also been
inspected. Not seen: any held-out fold case, any lockbox case, and any core-data selection. To keep this from
steering the result: the official pick is whatever the rule selects on the core data, whether or not it is the
preferred point; delta is defined by the noise rule above, not tuned; out-of-fold scoring reports both the registered
rule and this rule on the same held-out cases; the lockbox is used once.

## One fast mode, no ultra-fast mode

The ultra-fast pick is not shipped. On the final fit's screening data it is only slightly cheaper than the knee pick
on small cases and loses about 0.07 SP on large inputs (10,000 SP 0.678 vs 0.751). It stays in the report as a
result (section 7 already reports both picks on the lockbox, 50,000 and 200,000 inputs included).

Exception, decided now: if the knee pick fails lockbox criterion 4 (a 50,000 or 200,000 input does not complete in
4 h at 8 threads) and the ultra pick completes all of them, the ultra mode is reconsidered.

## Intended use

The knee rule's final pick is intended to become a Kalign mode, and the user intends it as the default mode later
on. Which preset name it replaces (fast or default) is decided after the gate.

## Interim candidate (screening data, final fit, generation 61; not the official pick)

The official pick comes from the core data after the fit stops. Numbers below are on the final fit's own screening
cases, so the Kalign row is optimistic; FAMSA was not tuned on them. CPU: total single-threaded CPU over the
screening cases (small set; extHomFam 1,000 + 10,000), same machine and load, relative to FAMSA.

| | Small SP | Small TC | 1k SP | 1k TC | 10k SP | 10k TC | CPU small | CPU large |
|---|---|---|---|---|---|---|---|---|
| FAMSA (defaults, 1 thread) | 0.780 | 0.637 | 0.753 | 0.595 | 0.787 | 0.593 | 30.2 s (1x) | 230.4 s (1x) |
| Kalign 4.0.0 fast | 0.776 | 0.618 | 0.678 | 0.515 | 0.633 | 0.471 | 16.2 s (0.54x) | 218.5 s (0.95x) |
| Knee pick (delta 0.01) | 0.782 | 0.628 | 0.787 | 0.641 | 0.751 | 0.537 | 7.1 s (0.23x) | 148.4 s (0.64x) |

Knee pick parameters: PFASUM60, gap open 8.12, gap extend 1.71, terminal gap extend 0.77, sequence weights 3.67,
VSM off, no consistency, no refine, realign 0; guide tree k-means + Prim, leaf size 925, LCS-indel distance, reduced
alphabet, k-means on BPM distance.

Knee rule (delta 0.01) on every fit's screening data: R 0.62-0.75 in folds 0-4, 0.54 in the final fit, 0.97 in the
seed-4321 repeat; consistency 0 in every fit but the repeat; leaf sizes 925-3,499.

Figures: manuscript repo `k4/explore/fastfit_front/` (`final_front_screen.pdf`, `final_selection_screen.pdf`).
