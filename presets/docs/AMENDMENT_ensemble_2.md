# Amendment 2 to the ensemble plan (PLAN_ensemble.md), 2026-10-15

Approved by the user on 2026-10-15, after generations 1-2 of the amendment-1 fits. The frozen plan is not edited;
this file changes it on top of `AMENDMENT_ensemble_1.md`.

## Diagnosis

In generations 1 and 2 every candidate of every fit (112 per fit) was infeasible, and only because of the CPU limits
(7.1 x small, 53.6 x at 1,000). The limits came from a pre-freeze measurement in which the fast setting cost 40.4 s on
the 300 small cases; the fit's evaluator (persistent workers, in-call CPU) measures 4.7 s, so the fast setting's cost
was overstated about 9 x and the limits were about 9 x too tight. Through the evaluator the pre-release design costs
29-40 x the fast setting on small cases and 30-65 x at 1,000. Many candidates already beat the fast setting on
precision while meeting the SP floors (e.g. final fit: P_small 0.757 vs 0.724 with SP 0.798 vs 0.782; P_1k 0.848 vs
0.801).

## Change

a. The CPU limits are removed: no CPU constraint and no CPU-based early stop. CPU stays objective 3 (R), so the knee
   rule's lowest-R choice still prefers cheaper points. A candidate with a real failure (error, timeout at 1,800 s)
   stays infeasible and is still stopped at its first failure.
b. Two fits only: fold 0 (the held-out estimate: gate criteria 1-5 on fold 0's test cases only) and the final fit (the
   shipped setting). Folds 1-4 are dropped.
c. CPU axis of the normaliser: R_hi = 3 is meaningless now (R is about 30-50). Over generations 1-2 of the fold-0 and
   final fits, the candidates with complete screening results had R from 30.98 to 49.75 (the others were stopped by
   the CPU early stop, so their R is unknown and higher). The axis maps log R from R_lo = 30 to R_hi = 50 (rounded
   outward to 10); the hypervolume reference point stays 1.1, so points above R of about 52.6 lie outside the
   hypervolume box (they still take part in the search and the selection).

Everything else is as amendment 1: no 10,000 input in the search, the 10,000 check at selection (24 families,
1,800 s, 1 thread, precision >= fast, SP >= fast - 0.005), at most 30 generations, the SP floors, the 1,800 s timeout,
failures infeasible.
