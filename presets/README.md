# presets: how Kalign 4's settings were chosen

This folder is the record of how the four built-in settings of Kalign 4
(protein and DNA/RNA, each `fast` and `accurate`) were obtained: the fitting
code, the plans written before each fit, the decisions and the results. It is
not needed to build or use Kalign.

## The procedure

Parameters were fitted with grouped cross-validation (related cases kept in the
same fold) on two pools of
reference alignments: protein (2,916 cases from BAliBASE 3, PREFAB 4, OXBench,
SABRE and BaliFam100) and nucleotide (3,225 cases from BRAliBASE II, MDSA and
Rfam seed alignments). Each fit followed a plan that was written and frozen
before it ran; changes after a freeze are dated amendments or deviations.

| Step | Plan and record | Result |
|---|---|---|
| Pooled re-fit, v1 and v1.1 (search) | `docs/PLAN.md` | `results/v1/`, `results/v1_1/` |
| v2: new picking rule on the v1.1 search | `docs/PLAN.md` (Amendment v2) | `results/v2/RESULTS_v2.md`, `results/v2/presets_v2.json` |
| Protein `fast` | `docs/PLAN_fastfit.md`, `docs/fastfit/` | `results/fastfit/final_pick.json` |
| Protein `accurate` (ensemble) | `docs/PLAN_ensemble.md`, `docs/AMENDMENT_ensemble_*.md`, `docs/PLAN_ensconf.md`, `docs/AMENDMENT_ensconf_*.md`, `docs/IMPL_NOTES_ensfit.md` | `results/ensfit/`, `results/ensconf/final_pick.json` |
| DNA/RNA `fast` and `accurate` | `docs/nucfit/PLAN_nucfit.md` | `results/nucfit/` |

## From the results to Kalign

`release/emit_c_presets.py` writes `lib/src/mode_presets.h` from the
`final_pick.json` files above and `results/v2/presets_v2.json`; the header
records the sha256 of every file it was made from. To check that the shipped
settings are exactly the fitted ones:

```bash
uv run --no-project python presets/release/emit_c_presets.py --check
```

`tests/preset_golden.cmake` (part of `ctest`) checks that the built program
gives the alignments these settings give, and `release/equiv_check.py`
compares the command line, the Python mode API and the fitting code's own
alignment path byte for byte.

## Data

No benchmark data is stored here. `resources/REGISTRY.yaml` lists every input:
download address, version, licence, retrieval date and sha256.
`recipe/download.py` downloads and verifies them, and `docs/step0/data.md`,
`docs/step0/cases_protein.md` and `docs/step0/cases_nucleotide.md` describe how
the cases are built from the downloads; `docs/step0/features.md` describes the
case features used to form the folds. Data, caches and run outputs go to a
working directory outside the repository (default `~/tmp/kalign-refit`).

## Code

- `recipe/`: the fitting pipeline (case lists, folds, search, scoring, nested
  cross-validation, selection, reports); `tests/`: its unit tests.
- `spaces/`: the two search-space definitions the recipe loads (protein
  `refit/space.py`, nucleotide `nested_nuc/nuc_space.py`), unchanged and
  checked against the sha256 values pinned in `recipe/presets_370.py`. Point
  the recipe at them with `export KALIGN_REFIT_TOOLS=presets/spaces`.
- `release/`: `emit_c_presets.py` and `equiv_check.py` (above).
- `verify/`: reproduction of the BAliBASE results.
- `run_all.sh`, `run_v2.sh`, `run_a2.sh`, `scripts/`: launch scripts for the
  fits, written for a 64-thread workstation.

## Notes on the record

- The fits ran from the Kalign development repository. Commits and tags named
  in the plans and records (for example `exp/refit-v2-freeze1`,
  `fastfit-freeze1`, `ens-plan-freeze`), and paths under `experiments/`, refer
  to that development history, which is not part of this repository.
- Each freeze manifest (`results/*/freeze1.json`) holds the sha256 of every
  file frozen for that fit. Files changed by later fits no longer match the
  earlier manifests; the code here is the state after the last fit.
- References to `k4/...` files are to the manuscript repository.
- Edited for publication (2026-11-11): local home-directory paths were
  replaced by `~/` or a description; nothing else changed. Original sha256,
  for checking against the development history:

  | File | Original sha256 |
  |---|---|
  | `resources/REGISTRY.yaml` | `1c6a17b06ee33a384fbd4af9f7c8d662f8220638d458cc0f2e7284d2ba8398c9` |
  | `results/fastfit/lockbox_gate.json` | `cd8d532041c97f2a7bad14e3cb0a872911fafe30c17737fb5917d9e76ded97d1` |
  | `results/v1/report/report.md` | `b567b5b8607a6f8a595f95840aa02cde195f9fb5e0a9a206ebff41848886a863` |
  | `results/v1/report/tables/steps.tsv` | `8d90d80a9155304ae6e33c4421e0269c4bd52e897ba72a605c84dbb1372f3525` |
  | `results/v1_1/report/report.md` | `6f68518d59a8468d2050333b9281b60775757fe48cd1d2286f0f896bdcdaa07e` |
  | `results/v1_1/report/tables/steps.tsv` | `b4eef1ece26780edb78bcd17c90eb9d3ca211b4a7353cc71316b7f608b94f896` |
  | `results/v1_1/addendum_a1/report_build/build.py` | `7285c9571bedcda706040891c08572ddf8fe02e8a99ab9b63c987dbcad2ce222` |
  | `results/v2/report_build/build.py` | `92c7dbe3e970c90083736212a69ae15488d6e09e2eba4687a6964678b38ac855` |
  | `spaces/nested_nuc/nuc_space.py` | `869c22f5d10f24ef7fea6be4e9b5cb110fc0a2215c5f28bf043f795cd05cc1a3` |
  | `docs/step0/cases_nucleotide.md` | `8856c38be6a2c1a89f526a11a8edac910bc59d491b1708631b37e61f9fb0e534` |
  | `recipe/presets_370.py` | `SPACE_SHA256` updated for the edited `nuc_space.py` |
