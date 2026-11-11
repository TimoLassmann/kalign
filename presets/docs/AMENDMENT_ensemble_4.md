# Amendment 4 to the ensemble plan (PLAN_ensemble.md), 2026-10-18

Approved by the user on 2026-10-18, after the final fit's search (30 generations) and its first core re-scoring, and
before any out-of-fold result exists (none has been computed). The frozen plan is not edited; this file changes it on
top of amendments 1-3. Code: `recipe/ensfit_select.py`, `recipe/ensfit_oof.py`, `recipe/ensfit_finish.py`.

## What was seen before this amendment (disclosed)

- The final fit's own re-scoring (recipe.ensfit, steps 3-5) scored the core with an early stop at the first failed
  case. Almost every 10,000 selection family timed out (1,800 s, 1 thread), so each candidate's core summary covered
  only about 40-56 of the training families at 1,000 and a fraction of the small core. Those summaries are not usable;
  an interim reading of them ("SP below the floor on the core") reported to the user was wrong for that reason.
- Under the plan's eligibility rule no candidate could pass the 10,000 part (every input must complete).
- The user decided (2026-10-16) that accurate mode is documented as not feasible above about 5,000 sequences.

## Changes

(a) **10,000 requirement dropped.** Eligibility (plan section 6 step 1) and the gate (section 7 criteria 4 and the
    10,000 part of 5) no longer use the 10,000 inputs. They are reported: the pick and the best point of every fit on
    the fold's 24 selection families at 10,000, and every arm on the out-of-fold families at 10,000 (1 thread,
    1,800 s per input; completion counts and accuracy on the families completed).

(b) **Core re-scoring without early stop, small + 1,000 only.** Every re-scored point runs every small core case and
    every training family at 1,000 to completion.

(c) **Re-scoring set widened.** The plan's set (front of the screening data, P_small >= fast - 0.01, at most 40) plus
    the 40 archive points with the highest screening P_1k among those that beat the fast setting on screening
    precision (small and 1,000) with screening SP at least the fast setting's (no 0.005 slack). Reason: the search
    drives the front to the edge of the SP floor on the 40 screening families, so front points are expected to fall
    just below it on new families; points with margin are the ones that can pass.

(d) **Declared fallback.** The pre-release design with the fast setting as member 1 (seed 1 of the search, key
    f1db0bfd...) is scored out of fold beside the fold picks. If the fold picks fail the gate and it passes, it
    becomes protein accurate; if neither passes, no accurate setting comes out of this experiment.

Unchanged: the searches (no fit is re-run), the SP floor (fast - 0.005, small and 1,000), the knee rule, the gate
criteria 1-3 and the 1,000 part of 5 (here: every out-of-fold input at 1,000 completes within 1,800 s on 1 thread),
bootstrap settings, all five folds (amendment 3), no lockbox (amendment 1).

## Out-of-fold evaluation (plan section 7, as amended)

Arms on every test case of each fold (small; families at 1,000 and 10,000): the fold's pick, the fast setting, the
fallback (d), the pre-release accurate preset (reported) and FAMSA 2 (default options, 1 thread; reported). A fold
without an eligible point counts as the fast setting for the procedure. Gate on the pooled rows of folds 0-4 vs the
fast setting: (1) small precision lower bound > 0; (2) small SP difference >= -0.005 and lower bound >= -0.015;
(3) precision at 1,000 lower bound > 0 and SP difference >= -0.005; (5) every input at 1,000 completes.

## Automation

`recipe.ensfit_finish` (detached on the workstation) runs the selection for each fit as soon as its fit process has
finished, then the out-of-fold run and the gate once folds 0-4 are selected. Outputs: `<fit>/selection_a4.json`,
`runs/ensfit/oof/rows.jsonl`, `runs/ensfit/oof/gate.json`; log `runs/ensfit/finish.log`.
