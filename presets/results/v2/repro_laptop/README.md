# Laptop reproduction of the v2 presets on BAliBASE 3 (2026-10-07)

Independent of the workstation pipeline: fresh build of Kalign from source (commit in `build_commit.txt`) into a
new venv on an Apple M3 (arm64), no cache, the original BAliBASE 3 release (lbgi.fr, `benchmarks/data/downloads/
bb3_release`, 386 cases, BB and BBS), and the official `bali_score` 3.01 (core blocks from the .xml files; built
with `-D_FORTIFY_SOURCE=0 -fno-stack-protector`, without which it aborts on macOS even on the reference itself).
Script: `presets/verify/repro_balibase.py`; comparison: `presets/verify/compare_repro.py`.

The presets are the final v2 presets, so this checks that the numbers reproduce, not overfitting (BAliBASE cases
were partly in the final fit's training data; the held-out numbers are in `../RESULTS_v2.md`).

Result (`compare_workstation.txt`): on the 385 cases both pipelines share, mean SP agrees to within 0.0003 for
all five arms; 97.7-99.2% of cases agree within 0.001 SP. The largest differences are BB20002 and BB40002, the two
cases whose drive5 copy differs from the original release, and a few cases where ARM and x86 floating point break
a near-tie differently. TC from bali_score is 0.004 lower on average in every arm (scorer definition), so
differences between arms are unaffected.

## Earlier Kalign 3 releases (same pipeline: original BAliBASE 3, official bali_score, laptop)

| BAliBASE 3, 386 cases, SP / TC | all | CPU s/case |
|---|---|---|
| Kalign 3.1 (v3.1, Oct 2019, the Kalign 3 paper release) | 0.775 / 0.381 | 0.098 |
| Kalign 3.3.5 (v3.3.5, Nov 2022) | 0.811 / 0.482 | 0.094 |
| Kalign 3.7.0 fast / accurate | 0.830 / 0.515, 0.848 / 0.537 | 0.042, 1.490 |
| new fast / accurate (v2 final presets) | 0.833 / 0.524, 0.852 / 0.567 | 0.108, 1.363 |

Builds: v3.1 is autotools and x86 AVX2 only; it was compiled natively for arm64 without AVX2 (`k31/build/`:
config.h without HAVE_AVX2, a stand-in header that supplies only `_mm_malloc`/`_mm_free` and the two vector types
named in prototypes). Its scalar distance code was probably never used in releases and sums distances in a
different order, so **the 3.1 numbers are not trustworthy** (RV30 SP 0.70 against 0.82 for 3.3.5); an x86 AVX2
build (workstation or Rosetta) is needed to confirm them. v3.3.5 was built with its own CMake switches
(`-DENABLE_SSE=OFF -DENABLE_AVX=OFF -DENABLE_AVX2=OFF -DUSE_OPENMP=OFF`, plus `-DCMAKE_POLICY_VERSION_MINIMUM=3.5`),
a supported configuration. CPU of the scalar builds is not comparable with the SIMD builds.

**Update (Addendum A2, workstation, x86 AVX2 build):** Kalign 3.1 scores SP 0.831 / TC 0.517 on BAliBASE 3 (385
cases, re-fit scorer), against 0.775 / 0.381 for the laptop's native non-AVX build. The laptop 3.1 numbers are a
build artefact and should not be used; see `../addendum_a2/TABLES.md`.
