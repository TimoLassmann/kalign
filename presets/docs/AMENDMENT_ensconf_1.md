# Amendment 1 to the ensconf plan (PLAN_ensconf.md), 2026-10-29

Written after the first complete run (selections, out-of-fold rows and gate of 2026-10-29 06:21 were seen; they are
kept in `runs/ensconf/_pre_amend1/` and reported). Decision by the agent under the user's standing delegation
(DECISIONS.md D9); the user had asked first for the failures to be understood.

## What happened
The plan limits each case to 1,800 s on one core (section 1). During the run the tree benchmark of the paper ran on
the same 64-core machine (about 116 processes), so the limit acted on starved processes: all 192 grid failures are
1,800 s wall-clock timeouts after 1,000-1,550 CPU-s. Re-run alone on the idle machine without a limit
(`~/tmp/ensconf_timeouts/rows.jsonl`, wall = CPU within 1 %): 148 finish within 1,800 s (both min_support 2 cases and
146 of 190 min_support 1 cases), 44 take longer (all min_support 1; up to 7,893 s).

## Change
Each grid row that failed with "timeout" and whose idle re-run finished within 1,800 s is replaced by the idle re-run:
its alignment is scored with the plan's scorer and stored with its CPU and wall time (harness tag `+amend1`); the 44
that exceed 1,800 s stay failures. Every replacement is logged (`runs/ensconf/amend1_rows.jsonl`). Then the
selection (section 4), the out-of-fold run and the gate (section 5) are recomputed unchanged. Nothing else changes.

## Why this and not another fix
It applies the plan's own limit as it was meant (1,800 s of one core), so a case's fate no longer depends on what
else ran on the machine. Keeping the starved timeouts would let scheduling decide the selection; dropping the limit
would change the plan. The 44 genuine overruns are a property of the method (the min_support 1 consensus,
RESULTS.md) and stay failures.
