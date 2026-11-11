# perf/dp: DP profiling and kernel benchmark harness

Research tooling only. Building it does not change Kalign: the option is off by
default, and with it off the instrumentation hooks in `lib/src` expand to
nothing. All 42 library objects have the same machine code (`otool -t`) as
`release/3.7.0`, and `kalign` produces byte-identical output for fast, default,
accurate and recall on `tests/data/{BB11001,BB12006,BB30014,BB40023}.tfa`,
`small.fa` and `tiny.fa`.

## Size limits (laptop-sized runs)

- Simulated data: at most **500 sequences** and **500 residues per sequence**.
  The grid is N in {50, 200, 500} x L in {150, 500} x 2 divergence levels
  (protein), plus DNA N=200, L=500. `make_datasets.py` and `run_profile.py`
  refuse anything larger. `dssim_gen` truncates sequences to 500.
- Use no more than 4 concurrent processes (1 or 2 is better), run kalign with
  `-n 1`, and build with `make -j4` at most.
- `accurate` on N=500, L=500 takes about 2 min and reaches **3.2 GB RSS**,
  mostly from the ensemble tables.
- Keep data, dumps and results outside the repository.

## Build

```bash
export SDKROOT=/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk   # macOS 27
cmake -S . -B $BUILD -DCMAKE_BUILD_TYPE=Release -DKALIGN_PERF=ON
make -C $BUILD -j4 kalign_prof kalign_plain dssim_gen bench_dp
```

| target | what it is |
|---|---|
| `kalign_perf_static` | a copy of the library compiled with `-DKALIGN_PERF_TIMERS` (region timers, plus the dump hook `KALIGN_DP_DUMP`) |
| `kalign_prof` | driver: read, mode preset, align, write; `-j t.json` writes the timer breakdown. The alignment equals `kalign --mode M -n T` |
| `kalign_plain` | the same driver linked to the normal library (used to measure overhead; output parity) |
| `dssim_gen` | FASTA generator built on Kalign's DSSim (`tests/dssim.c`) |
| `bench_dp` | DP microbenchmark and oracle that runs on dumped problems |

Instrumentation: `perf/dp/instr/kperf.h`. Hooks are one macro on the opening
line of a function (`KPERF_REGION`, `KPERF_PHASE_REGION`, `KPERF_DP`,
`KPERF_RUNNER`), so line numbers stay the same. When `KALIGN_PERF_TIMERS` is
not defined, `lib/src/tldevel.h` defines them as empty. Timer overhead is 1.5
to 2.5% of wall time (`kalign_plain` vs `kalign_prof`: 7.97 vs 8.10 s,
36.1 vs 37.0 s, 0.236 vs 0.240 s).

## 1. Regenerate the data

```bash
G="bash -c"     # or a wrapper of your own, for example a memory guard
$G "uv run --no-project python perf/dp/scripts/make_datasets.py --dssim-gen $BUILD/perf/dp/dssim_gen --out $DATA"
```

Seeds are `20261002 + 100*iN + 10*iL + iDiv` for protein and `20266002` for
DNA, and they are recorded in `$DATA/manifest.json`. Divergence is set by
DSSim's `n_obs`: lo = 20 gives about 81% mean pairwise identity, hi = 2 gives
about 44%; DNA uses n_obs = 8 (about 76%). Use `uv run --no-project`, because
plain `uv run` inside the worktree tries to build the Python package.

## 2. Profile

```bash
$G "uv run --no-project python perf/dp/scripts/run_profile.py --prof $BUILD/perf/dp/kalign_prof \
     --data $DATA --results $RES --jobs 2 --timeout 150 \
     --skip prot_N500_L500_lo__accurate prot_N500_L500_hi__accurate"
uv run --no-project python perf/dp/scripts/analyze_profile.py --results $RES --md $RES/summary.md
```

Each JSON file holds self and inclusive ns, call counts and DP cells per
(phase, region). The phases are tree, consbuild, progressive, refine, realign
and ensemble, so a kernel is counted separately in each phase. `dp_*` regions
are kernels without the consistency bonus and `dpc_*` regions are kernels with
it. Profile with `-n 1`, because phases are not propagated to thread-pool
workers.

## 3. Dump DP problems (opt-in)

The instrumented library writes one `.kdp` file per top-level
`aln_runner` call, that is, one complete Hirschberg alignment, when
`KALIGN_DP_DUMP=<dir>` is set. Optional variables:
`KALIGN_DP_DUMP_EVERY=k` (every k-th call per phase and kind),
`KALIGN_DP_DUMP_MAX=n` (cap per phase and kind),
`KALIGN_DP_DUMP_MIN_CELLS` and `KALIGN_DP_DUMP_TAG` (filename prefix; must be
unique per run). The format is in `instr/kdp_format.h` (v2). A dump holds the
sequences or profiles (64 floats per column), the gap and substitution
parameters, the entry states, the real run's raw path, and the sparse
consistency bonus if the call used one.

The reference sample (216 problems, 21 MB) was produced like this:

```bash
$G "cd $REPO && P='uv run --no-project python perf/dp/scripts/run_profile.py --prof $BUILD/perf/dp/kalign_prof --data $DATA --results $DUMPRUN --jobs 2 --timeout 120 --dump $DUMPS'; \
   \$P --modes fast default --only prot_N200_L500_hi prot_N500_L150_lo prot_N50_L500_lo dna_N200_L500 --dump-every 11 --dump-max 5; \
   \$P --modes accurate --only prot_N200_L150_hi prot_N50_L500_lo --dump-every 499 --dump-max 6"
```

The sample has 90 ss, 63 sp and 63 pp problems from the tree, consbuild,
progressive and refine phases. Sizes run from 142x142 to 737x492 (profiles are alignment columns, so they can be longer than 500), and 90
of the problems carry a consistency bonus.

## 4. Benchmark and check

```bash
B=$BUILD/perf/dp/bench_dp
$G "$B --record $ORACLE --selfcheck --json ref.json -q $DUMPS"                 # baseline timing + oracle (about 35 s)
$G "$B --kernel my_kernel --check $ORACLE --json mine.json -q $DUMPS"          # prototype: timing + check
$G "$B --kernel my_kernel --ops none --check $ORACLE -q $DUMPS"                # correctness only (a few s)
$B --list
```

- Timed operations: `fwd` and `bwd` are one full pass each. `align` is the full
  Hirschberg alignment, about twice the cells of one pass. Timings are medians
  over at least 0.05 s and at least 3 repetitions, single-threaded. ns/cell is
  computed against len_a x len_b for all three operations.
- `--check` compares, bit for bit, the last forward and backward state rows,
  the top-level meetup score with its meet column and transition (the global
  DP score), and the path. It reports the top-level score difference and an
  approximate path rescoring difference (see `bench/kdp_score.c`; treat it as
  a relative measure only).
- `--selfcheck`: the `ref` kernel reproduces the path recorded by the real
  Kalign run on 175/175 problems. Problems that used midpoint flips are
  excluded.
- `--no-bonus` drops the consistency bonus. `--kind`, `--max-files` and
  `--min-time` restrict or speed up a run.

### Adding a prototype kernel

Copy `kernels/example_pp_dense.c` to `kernels/<name>.c`, then fill in any of
`.fwd/.bwd/.meet[KDP_SS|KDP_SP|KDP_PP]`, or a complete `.align`. Slots you
leave out fall back to `ref`. Register with `KDP_REGISTER_KERNEL(var)`. Then run
`cmake $BUILD && make -C $BUILD -j4 bench_dp`; the glob picks up the new
file. Semantics and data layout are documented in `bench/kdp_kernel.h`. Your
kernel must apply `sparse_bonus_lookup()` with the same indices as the
reference if it is to match on bonus problems.

The example kernel replaces the sparse pp dot product with a dense 23-term
loop. Its paths are identical (63/63) but its rows are not bit-identical, and
it is about 3.5x slower: without -ffast-math the float reduction is a serial
chain.

**Code alignment matters.** The reference kernels are compiled into `bench_dp`
with `-falign-functions=64` (`KDP_ALIGN_FLAGS`). Purely from link layout, the
ss forward pass measured 1.84 ns/cell in one build and 2.53 ns/cell in
another. pp timings still vary by about ±5% between alignment variants and
runs. Compare a prototype against `ref` built and run in the same binary, and
repeat the measurement.

## Baseline (MacBook Pro M3 Max, single thread, 2026-10-02)

`bench_dp`, kernel `ref`, 216 dumped problems; one full pass takes 34 s with a
peak RSS of 3.5 MB:

| kind | n | fwd ns/cell | bwd ns/cell | align ns/cell | fwd Gcell/s |
|---|---|---|---|---|---|
| ss (seq-seq) | 90 | 1.84 | 1.85 | 3.76 | 0.54 |
| sp (seq-profile) | 63 | 1.88 | 1.87 | 3.85 | 0.53 |
| pp (profile-profile) | 63 | 2.43 | 2.34 | 4.99 | 0.41 |

With and without the bonus: pp 2.77 vs 2.13 fwd ns/cell; sp 1.91 vs 1.84;
ss 1.86 vs 1.84. Simulated profiles are sparse: a pp prof1 column has 1.45
nonzero residues on average, so real, more diverse data will make pp slower.
Inside real runs (timer self time divided by cells) the numbers match: ss and
sp 1.85, pp 2.3 to 2.8, pp with bonus 3.2 to 3.3 ns/cell.

Where the time goes (protein grid, sum over datasets, % of wall time):

| mode | total s | guide tree (all) | tree via full DP | all DP fwd+bwd+meet | profile build/update | extras (VSM, bonus, ensemble) |
|---|---|---|---|---|---|---|
| fast | 5.0 | 55.9 | 0 | 39.9 | 2.2 | 0 |
| default | 25.3 | 11.6 | 0 | 82.1 | 1.7 | 3.2 |
| accurate | 234.6 | 63.8 | 60.3 | 70.2 (incl. tree DP) | 0.2 | 20.5 |

Log-log slopes (time vs N, vs L): DP 1.0 / 1.9; heuristic tree 1.1 / 1.3;
profile 1.0 / 1.0. `accurate` overall is 1.3 / 1.8, its tree DP is 1.3 / 1.95,
and its ensemble extras are 1.8 to 1.9 / 1.3 to 1.4. `analyze_profile.py`
prints the full per-dataset table.

The consistency-bonus column mismatch noted while writing the scorer was
confirmed: the forward kernels looked the bonus up one column off. It is fixed
in d44c237 together with three other DP correctness fixes (ad3be1b, 2a57d36,
2c9a647); the oracle in this README's commands is recorded on the fixed code.

## Exact speed-ups built with this harness (2026-10-02)

All of these leave Kalign's output byte-identical (checked against the 3.7.0
release build on the 13 simulated sets plus tests/data, all modes, 1 and 4
threads, on Apple M3 Max and on an AMD Threadripper 3970X).

| step | what | modes |
|---|---|---|
| A1 | guide tree: reuse pick_anchor's BPM distances; each leaf pair once | all |
| A2 | BPM match table in one pass | all |
| A3 | ensemble POAR: merge per pair, scan instead of binary search | ensembles |
| A4/A5 | wavefront SIMD DP passes (`-DKALIGN_WAVE=ON`), NEON and AVX2 via `lib/src/wave_simd.h`; `lib/src/aln_wave.c` is generated by `tools/gen_aln_wave.py` from `kernels/{ss,sp,pp}_wave.c` | all |
| A6 | smaller per-call padding in the wave kernels | all |
| A7 | BPM: branch-free block step, constant block count | all |

DP kernels (bench_dp, all dumps; full Hirschberg alignment, ns/cell, ref -> wave):
M3 ss 2.7 -> 1.2, sp 2.9 -> 1.65, pp 3.8 -> 2.0; Threadripper ss 7.9 -> 2.05,
sp 8.1 -> 3.8, pp 10.3 -> 4.4.

Tried and not used (exact, but no gain): Farrar striped layout
(`kernels/ss_striped.c`; the lazy fix-up does as much work as the main loop on
Kalign's global DP), per-strip vector step buffers instead of lane extraction
(slower on x86), a SIMD-filtered meetup (`kernels/meet_wave.c`), a scalar
wavefront BPM (`tools/bpm_proto.c`; no better than branch-free alone).
Not exact, not used: summing the pp substitution before adding it to the cell
(`pp_wave*_sub` in history; ~10% faster, last-bit differences).
SIMD BPM across 64-bit blocks: done (8c82cde).
Checked and not built (2026-10-03):
- Two BPM computations interleaved (`tools/bpm_proto.c`, `bpm_simd_x2`), exact:
  vs the SIMD library BPM, M3 -3% / -10% time (length 150 / 500), x86 -12% /
  +9% (register pressure).  Mixed and small; would need call-site changes in
  pick_anchor and d_estimation.
- Continuous wavefront across strips: the 7-step strip edges are at most
  about 3% of full-alignment time (from the per-width profile, KDP_HIRSCH_PROF);
  the per-call fixed costs dominate narrow sub-calls instead.
