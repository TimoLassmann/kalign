# Amendment 2 to PLAN_ensconf.md: selection instead of consensus (2026-10-29)

Status: written before any selection-mode alignment of a benchmark case exists. Not blind: the folds, members and
grid picks are those of the ensconf re-fit (amendment 1), whose out-of-fold results we have seen. Stated wherever
this is reported.

## Question

With the same members (each fold's pick of amendment 1: its ensemble members, member 1's guide tree, its tree noise
sigma), is it better to output the member alignment the ensemble supports most, with the ensemble's confidence
projected onto it, than to build the consensus alignment from the POARs?

## The selection mode (kalign exp/ens-select 5d56876, on exp/sl-tree-noise ea02ca9)

`kalign_ensemble_config.select_only = 1` (Python `select_only=1`):
1. the K members are aligned and their POARs collected as before;
2. each member k is scored by the existing POAR score: sum over its aligned residue pairs of
   (support - 1) / (K - 1), support = number of members aligning the pair (the member's own vote excluded);
3. output = the member with the highest score (ties: lower index; no preference for member 1, no 5 % rule);
   no consensus, no refinement of the winner;
4. confidence = the existing per-residue confidence on the output (share of the K members aligning each residue
   with each other residue of its column, own vote included, so at least 1/K for a column with two or more
   residues); column confidence = mean over its residues. Unchanged code.
The selected member index is returned (`selected_member`).

## Arms and data

Out of fold, the test cases of folds 0-4 (small; families at 1,000), build 5d56876:
- fast: the fast setting;
- pick: each fold's amendment-1 pick (consensus at its min_support) - identical parameters to the ensconf
  out-of-fold pick arm;
- select: the same parameters plus select_only = 1.
Fold 2 has no pick: in the procedure arms pick and select both count as fast there (as in amendment 1).

## Known-answer tests (before any result is reported)

- S1: on build 5d56876, fast and pick reproduce the ensconf out-of-fold rows of build ea02ca9 on every test case:
  precision, recall, TC, width ratio and pair Brier equal (|d| <= 1e-12).
- S2: on 10 small cases and 2 families at 1,000 (fold 0 test cases, final pick parameters), the select output is
  byte-identical to the selected member aligned alone (n_runs = 1, the member's parameters, seed 42 + k).
- S3: on the same cases, the returned column confidence equals an independent Python computation from the K member
  alignments (|d| <= 1e-5, float).
- S4: on the same cases, the selected member equals the argmax of an independent Python computation of the POAR
  scores of the K members (exact rational arithmetic; ties to the lower index).

## Decision rule

The selection mode replaces the consensus as protein accurate only if all hold:
1. select passes the eight gate criteria of PLAN_ensconf.md section 5 against fast (procedure, folds 0-4);
2. select minus pick (folds 0, 1, 3, 4, where they differ): SP and TC, small and 1,000, lower bound of the 95 %
   interval >= -0.005;
3. select minus pick: pair Brier, small and 1,000, upper bound <= +0.005;
4. families at 1,000 with SP more than 0.2 below fast: select <= 3.
Otherwise the consensus pick (amendment 1) stays. Bootstrap as the gate (10,000, seed 1234, groups / families).

Reported beside the rule, not used by it: CPU, selected-member counts, pair-level calibration (ECE over test-aligned
scored residue pairs, 10 equal-width confidence bins, pairs pooled within the small pool and within the families
at 1,000; reported for pick and select, and for fast as a hard alignment with confidence 1).
