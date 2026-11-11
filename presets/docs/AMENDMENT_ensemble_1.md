# Amendment 1 to the ensemble plan (PLAN_ensemble.md), 2026-10-15

Approved by the user on 2026-10-15, after the pilot and before any registered generation was scored. The frozen plan
(tag `ens-plan-freeze`) is not edited; this file changes it.

## Reason

The pilot (fold 0, generation 1, `runs/ensfit/pilot_fold0`) had 56/56 candidates above the 600 s CPU limit of the
10,000 timing guard (fold 0: cyclin; fast setting 4.0 s), including both seeds (the pre-release design with the fast
setting as member 1). Because a candidate stops at its first failed case, every candidate became infeasible with
flat constraint violations, which gives the search nothing to follow.

## Change

a. The 10,000 guard input and its constraint are removed from the search. Screening = the 300 small cases + the 40
   screening families at 1,000.
b. The 10,000 check stays only at selection, on the 24 families per fold, exactly as in plan section 6 (every input
   completes in 1,800 s at 1 thread; precision >= fast setting; SP >= fast setting - 0.005).
c. Generation cap 30 instead of 60 (the rest of the stop rule is unchanged: at least 20 generations, HV gain < 0.1 %
   over 5).
d. Unchanged: the CPU limits (7.1 x small, 53.6 x at 1,000), the recall floors, K = 5, member 1 = the fast setting,
   the nested cross-validation and the knee rule.

## Lockbox

Per plan section 9, a changed plan cannot use the lockbox. Gate criterion 6 is dropped: the gate is out of fold only
(criteria 1-5). The 50,000 and 200,000 inputs move to the release timing block.
