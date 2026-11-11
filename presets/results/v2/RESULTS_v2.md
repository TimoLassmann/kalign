# Re-fit v2: results and record (2026-10-07)

v2 defines the presets that ship. v1 and v1.1 (`../v1/`, `../v1_1/`) and Addendum A1 (Kalign 3.6.1) are kept and
reported alongside it. Plan and rules: `presets/docs/PLAN.md`, sections "Amendment v1.1", "Addendum A1",
"Amendment v2" and "Deviation v2-1".

## What was done, in order

| step | what | where |
|---|---|---|
| Search | Two-stage NSGA-III on each fit's training core: folds 0 to 4 and the final fit ("all"), both pools. Run once in v1, continued to convergence in v1.1 (stage-2 cap 40 to 140 generations). Not repeated in v2. | tags `exp/refit-v1-freeze1`, `exp/refit-v1.1-freeze1` |
| Addendum A1 | Kalign 3.6.1 and 3.7.0 command-line programs on every case, same harness and scorer as the external aligners. Descriptive only. | `../v1_1/addendum_a1/` |
| Amendment v2 | New picking rule (below), written and tagged before any v2 number: `exp/refit-v2-freeze1`. Every shipped preset comes from the procedure; the v1 "keep 3.7.0 unless beaten" decision rule is dropped. | `recipe/select_v2.py`, `run_v2.sh` |
| Pick | The rule applied to the saved Pareto front of every fit (5 folds + final), using only that fit's training cases. | `*/selection_fold*.json`, `*/selection_final.json` |
| Score held out | Each fold's picks aligned on the cases of that fold, which its fit never saw (`recipe.nested test`; mostly from the cache, 2 minutes). Rows hashed (Freeze 2) before the final picks were scored. | `freeze2.json` |
| Score in sample | The final picks on every case (`recipe.nested insample`), only to estimate optimism. | |
| Pool | `recipe.pool_results`: per-case rows pooled over the five held-out folds; 95% BCa intervals, 10,000 resamples. Deviation v2-1 (a missing field in the pooling code) was fixed under tag `exp/refit-v2-fix1`; Freeze 2 verified intact afterwards. | `protein/`, `nucleotide/` |

## Picking rule (Amendment v2)

Quality q = (SP + TC) / 2 on the fit's training core (each benchmark weighted equally); CPU = total single-thread
CPU on that core. A point is acceptable within a CPU limit L if no point costing more (up to L) gains at least
0.002 in q per doubling of CPU. **accurate** = cheapest acceptable point (no limit). **fast** = cheapest acceptable
point within L = CPU of the fastest other aligner on the same core / 5. default = recall = accurate. Precision is
reported, not used. The two constants were set after looking at the final-fit fronts (in sample); the plan says so.

## What "held out" means here

Every case is scored exactly once, by the presets picked in the fold that never saw it (search and pick both used
only the other four folds). The five folds' picks are different parameter sets; the pooled held-out score
estimates how well the procedure (search, then pick by the rule) does on new data. The shipped presets are one more
run of the same procedure on all data. Related cases share a fold (shared reference sequences; Rfam clans).

## Held-out results (benchmark-weighted; CPU per case at 30 concurrent jobs)

| protein, 2,916 cases | SP | TC | precision | CPU (s) |
|---|---|---|---|---|
| Kalign fast | 0.747 [0.735, 0.758] | 0.569 [0.549, 0.589] | 0.683 | 0.05 |
| Kalign accurate | 0.770 [0.759, 0.780] | 0.602 [0.582, 0.621] | 0.718 | 0.95 |
| Clustal Omega | 0.761 [0.749, 0.772] | 0.593 [0.572, 0.612] | 0.712 | 0.42 |
| MAFFT L-INS-i | 0.784 [0.772, 0.794] | 0.627 [0.606, 0.646] | 0.720 | 3.66 |
| MUSCLE 5 | 0.796 [0.786, 0.805] | 0.639 [0.619, 0.657] | 0.729 | 27.6 |

| nucleotide, 3,225 cases | SP | TC | precision | CPU (s) |
|---|---|---|---|---|
| Kalign fast | 0.833 [0.822, 0.843] | 0.691 [0.676, 0.705] | 0.830 | 0.005 |
| Kalign accurate | 0.841 [0.830, 0.851] | 0.701 [0.686, 0.715] | 0.838 | 0.021 |
| Clustal Omega | 0.808 [0.797, 0.818] | 0.663 [0.649, 0.677] | 0.809 | 0.040 |
| MAFFT L-INS-i | 0.831 [0.819, 0.842] | 0.696 [0.681, 0.710] | 0.831 | 0.33 |
| MUSCLE 5 | 0.864 [0.855, 0.873] | 0.736 [0.723, 0.747] | 0.861 | 1.55 |

On protein, Kalign accurate uses about 2 times the CPU of Clustal Omega and a quarter of MAFFT L-INS-i's; Kalign
fast uses about 1/8 of Clustal Omega's. Against 3.7.0 accurate, protein accurate gains +0.013 SP and +0.023 TC.

Optimism (final presets in sample minus held out, benchmark-weighted): protein fast +0.004 SP / +0.003 TC, protein
accurate −0.002 / +0.006, nucleotide both modes ≤ +0.003. Case-weighted, protein fast reaches +0.009.

## BAliBASE 3 by reference category, held out (SP / TC; BB and BBS together)

| category | n | Kalign fast | Kalign accurate | Clustal Omega | MAFFT L-INS-i | MUSCLE 5 |
|---|---|---|---|---|---|---|
| RV11 | 76 | 0.633 / 0.398 | 0.667 / 0.433 | 0.608 / 0.378 | 0.678 / 0.463 | 0.751 / 0.539 |
| RV12 | 88 | 0.925 / 0.823 | 0.926 / 0.817 | 0.905 / 0.789 | 0.939 / 0.849 | 0.947 / 0.870 |
| RV20 | 82 | 0.919 / 0.413 | 0.921 / 0.456 | 0.924 / 0.474 | 0.935 / 0.499 | 0.934 / 0.504 |
| RV30 | 60 | 0.843 / 0.506 | 0.854 / 0.542 | 0.863 / 0.589 | 0.873 / 0.612 | 0.883 / 0.624 |
| RV40 | 48 | 0.870 / 0.529 | 0.896 / 0.508 | 0.900 / 0.595 | 0.914 / 0.602 | 0.909 / 0.558 |
| RV50 | 31 | 0.814 / 0.434 | 0.866 / 0.546 | 0.858 / 0.517 | 0.904 / 0.594 | 0.913 / 0.629 |

Largest change of Kalign accurate against 3.7.0 accurate: RV11 +0.032 SP / +0.035 TC, RV50 +0.076 TC; RV40 is the
one category where it is worse (−0.034 TC). In sample (final presets; `figures/fig_bali_categories_insample.*`)
the category means differ from the held-out ones by −0.011 to +0.025 SP and −0.022 to +0.055 TC; they are not
uniformly higher, because the final presets are different parameter sets from the fold presets and the final fit
saw only 175 of the 385 BAliBASE cases.

## The shipped presets (final fit; full precision in `presets_v2.json`)

Matrix codes: PFASUM43, PFASUM60, CorBLOSUM66 (protein), Kimura 200 PAM (nucleotide). Refinement: none / all /
confident / inline. Gap penalties: gpo open, gpe extend, tgpe terminal extend. vsm_amax: variable scoring matrix
offset. No tree noise in any run.

**Protein fast** (1 run; (SP+TC)/2 = 0.716 on the final core, 0.044 s per case)

| run | matrix | gpo | gpe | tgpe | vsm_amax | refine |
|---|---|---|---|---|---|---|
| 1 | PFASUM43 | 6.190 | 0.809 | 0.600 | 0.000 | none |

shared: seq_weights 2.801, consistency anchors 3, consistency_weight 2.892, realign 1, min_support 0, fill_in off.

**Protein accurate** (5-run ensemble; 0.741, 0.48 s per case)

| run | matrix | gpo | gpe | tgpe | vsm_amax | refine |
|---|---|---|---|---|---|---|
| 1 | CorBLOSUM66 | 7.246 | 1.271 | 1.102 | 0.471 | confident |
| 2 | PFASUM60 | 8.683 | 1.482 | 1.133 | 0.760 | none |
| 3 | PFASUM60 | 8.667 | 0.623 | 0.598 | 1.488 | none |
| 4 | PFASUM60 | 10.046 | 1.217 | 1.729 | 2.520 | confident |
| 5 | PFASUM60 | 9.938 | 0.645 | 0.696 | 2.510 | none |

shared: seq_weights 3.273, consistency anchors 8, consistency_weight 3.558, realign 2, min_support 2,
consistency_merge off, fill_in realign (fill_anchor_frac 0.080, fill_support 1), vsm_mode 0.

**Nucleotide fast** (1 run; 0.769, 0.0066 s per case)

| run | matrix | gpo | gpe | tgpe | vsm_amax | refine |
|---|---|---|---|---|---|---|
| 1 | Kimura 200 PAM | 15.584 | 0.615 | 2.904 | 0.294 | none |

shared: seq_weights 0.428, consistency anchors 1, consistency_weight 1.198, realign 0.

**Nucleotide accurate** (1 run; 0.778, 0.022 s per case)

| run | matrix | gpo | gpe | tgpe | vsm_amax | refine |
|---|---|---|---|---|---|---|
| 1 | Kimura 200 PAM | 15.884 | 1.362 | 2.601 | 0.519 | confident |

shared: seq_weights 0.560, consistency anchors 5, consistency_weight 3.568, realign 1.

Not yet done: compiling these into `lib/src/aln_wrap.c` (the command line still ships the 3.7.0 presets).

## How the figures were made

All charts are drawn by `presets/recipe/report_figures.py` from the run's own files, in the andamentum report style
(`setup_style("report")`, the four andamentum tones, text kept as text, inline SVG); each figure has a CSV of the
plotted numbers beside it. Run with the andamentum tool's Python:
`~/.local/share/uv/tools/andamentum/bin/python -m recipe.report_figures ...` or by calling the functions below.

| figure | function | inputs | held out? |
|---|---|---|---|
| `figures/fig_tradeoff_v2_{protein,nucleotide}` | `fig_tradeoff_v2` | `<pool>/oof_bw.tsv` (v2 pooling, includes 3.6.1 and externals) | yes |
| `figures/fig_bali_categories` | `fig_bali_categories(scope="oof")` | fold `oof_rows.tsv` (v2) + external rows | yes |
| `figures/fig_bali_categories_insample` | `fig_bali_categories(scope="insample")` | final `insample_rows.tsv` (v2) + external rows | no: optimistic |
| `../v1_1/figures/fig_mode_choice_*` | `fig_mode_choice` | final-fit front (`selection.json`), externals and 3.6.1 on the same core | no: training core |

Intervals: BCa over cases for the overall tables; percentile bootstrap over cases (2,000 resamples, seed 1234) for
the BAliBASE categories. CPU: external aligners and Kalign at 30 concurrent jobs (inflated 5 to 30% against an
idle machine, equally for all); Kalign 3.6.1 at 8 jobs (Addendum A1).

## Proposed paper figures

1. **Procedure** (schematic): nested CV, freezes, amendments v1.1 and v2.
2. **Speed between releases**: 3.7.0 against 3.6.1 per mode, same harness (`../v1_1/figures/fig_speed`); the
   exact speed-ups at equal output. Missing and worth adding: CPU and wall time against the number of sequences
   (e.g. 100 to 10,000, BaliFam1000/10000 or HomFam), and wall time against threads; benchmark cases are small.
3. **Held-out accuracy against CPU** (`fig_tradeoff_v2_*`): the main result.
4. **BAliBASE categories with CPU** (`fig_bali_categories`): where Kalign gains and where it lags (RV40).
5. **Mode choice** (`fig_mode_choice_*`): the final-fit curve and the knee rule, marked as training-core scores.
Supplement: convergence, optimism, v1 / v1.1 / v2 side by side, 3.6.1 comparison, in-sample BAliBASE categories.
