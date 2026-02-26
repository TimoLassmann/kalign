[![CMake](https://github.com/TimoLassmann/kalign/actions/workflows/cmake.yml/badge.svg)](https://github.com/TimoLassmann/kalign/actions/workflows/cmake.yml)
[![Python](https://github.com/TimoLassmann/kalign/actions/workflows/python.yml/badge.svg)](https://github.com/TimoLassmann/kalign/actions/workflows/python.yml)
[![Build Python Wheels](https://github.com/TimoLassmann/kalign/actions/workflows/wheels.yml/badge.svg)](https://github.com/TimoLassmann/kalign/actions/workflows/wheels.yml)
![CodeQL](https://github.com/TimoLassmann/kalign/workflows/CodeQL/badge.svg)
![GitHub stars](https://img.shields.io/github/stars/TimoLassmann/kalign)
![GitHub issues](https://img.shields.io/github/issues/TimoLassmann/kalign)

# Kalign

Kalign is a fast multiple sequence alignment program for biological sequences written in C with Python bindings.

## 🚀 Key Features

- **🔥 High Performance**: Fast multiple sequence alignment with multi-threading support
- **⚡ Smart Threading**: Auto-detects CPU cores and uses N-1 threads by default (max 16) for optimal performance
- **🔧 Cross-Platform**: Works on Linux and macOS with multiple build systems (CMake, Zig)
- **📊 Multiple Formats**: FASTA, MSF, Clustal, Stockholm, PHYLIP support
- **🧬 Sequence Types**: Optimized for protein, DNA, RNA, and divergent sequences
- **⚡ SIMD Optimizations**: Vectorized code for x86_64 systems (SSE4.1, AVX, AVX2)
- **🐍 Python Integration**: Modern Python package with comprehensive bioinformatics ecosystem support

## Installation

### From Source (Primary)

#### Prerequisites

- **C compiler** (GCC, Clang, or MSVC)
- **CMake** (3.18 or higher)
- **OpenMP** (optional, for parallelization)

#### Basic Build

```bash
# Download and extract latest release
tar -zxvf kalign-<version>.tar.gz
cd kalign-<version>

# Build
mkdir build && cd build
cmake ..
make
make test
make install
```

#### macOS with Homebrew

On macOS, install dependencies first:

```bash
# Install dependencies
brew install cmake

# For OpenMP support (recommended)
brew install libomp

# Clone and build
git clone https://github.com/TimoLassmann/kalign.git
cd kalign
mkdir build && cd build
cmake ..
make
make test
make install
```

**Note**: On macOS, Kalign automatically configures OpenMP with Homebrew's libomp installation at `/opt/homebrew/opt/libomp/`.

#### Alternative Build Systems
When using `zig`, version `0.12` is required to build this project. 

**Zig Build** (for cross-compilation):
```bash
zig build
```

**Debug Build**:
```bash
cmake -DCMAKE_BUILD_TYPE=Debug ..
make
```

**Without OpenMP**:
```bash
cmake -DUSE_OPENMP=OFF ..
make
```

### Python Package

For development or latest features, install from source:

```bash
pip install git+https://github.com/TimoLassmann/kalign.git
```

For enhanced bioinformatics ecosystem integration:

```bash
pip install "kalign[biopython] @ git+https://github.com/TimoLassmann/kalign.git"    # + Biopython integration
pip install "kalign[skbio] @ git+https://github.com/TimoLassmann/kalign.git"        # + scikit-bio integration  
pip install "kalign[all] @ git+https://github.com/TimoLassmann/kalign.git"          # Full ecosystem support
```

## Usage

### Command Line Interface

```bash
Usage: kalign [--fast | --precise] -i <seq file> -o <out aln>

Modes:

   (default)          : Consistency anchors + VSM (best general-purpose)
   --fast             : VSM only, no consistency (fastest)
   --precise          : Ensemble(3) + VSM + realign (highest precision)

Options:

   --format           : Output format. [Fasta]
   --type             : Alignment type (rna, dna, protein, divergent, internal)
   --gpo              : Gap open penalty. [auto]
   --gpe              : Gap extension penalty. [auto]
   --tgpe             : Terminal gap extension penalty. [auto]
   --refine           : Refinement mode: none, all, confident. [none]
   -n/--nthreads      : Number of threads. [auto: N-1, max 16]

Ensemble options:

   --ensemble N       : Run N ensemble alignments. [off; 5 if no value given]
   --ensemble-seed    : RNG seed for ensemble. [42]
   --min-support N    : Explicit consensus threshold. [auto]
   --save-poar FILE   : Save POAR table for later re-thresholding.
   --load-poar FILE   : Load POAR table (skip alignment).
   --version (-V/-v)  : Prints version.
```

Most users only need to choose a mode. Explicit parameters override mode defaults
(e.g., `--fast --ensemble 5` uses fast base settings with 5 ensemble runs).

#### Threading Behavior

Kalign automatically detects your system's CPU cores and uses N-1 threads by default (leaving one core free), with a maximum of 16 threads.

- **Auto-detection**: Uses CPU cores - 1 (e.g., 15 threads on a 16-core system)
- **Maximum cap**: Never uses more than 16 threads
- **Manual override**: Use `-n/--nthreads` to specify a custom thread count
- **Single-threaded**: Use `-n 1` to disable parallelization

### Input Formats

Kalign accepts:
- **Unaligned sequences**: FASTA format
- **Pre-aligned sequences**: FASTA, MSF, or Clustal format (gaps will be removed and sequences re-aligned)

### Sequence Types

Kalign automatically detects sequence types but offers manual control via `--type`:

- **`protein`**: Uses CorBLOSUM66_13plus substitution matrix (default for protein)
- **`divergent`**: Uses Gonnet 250 substitution matrix for highly divergent proteins
- **`dna`**: DNA parameters (match: +5, mismatch: -4, gap open: -8, gap ext: -6)
- **`rna`**: Optimized parameters for RNA alignments
- **`internal`**: Like DNA but encourages internal gaps (terminal gap penalty: 8)

Fine-tune with `--gpo` (gap open), `--gpe` (gap extension), and `--tgpe` (terminal gap extension).

### Python API

```python
import kalign

sequences = [
    "ATCGATCGATCG",
    "ATCGTCGATCG",
    "ATCGATCATCG"
]

# Default mode — consistency anchors + VSM (best general-purpose)
aligned = kalign.align(sequences)

# Fast mode — no consistency, fastest
aligned = kalign.align(sequences, mode="fast")

# Precise mode — ensemble + realign, highest precision
aligned = kalign.align(sequences, mode="precise")

# Expert: override specific parameters within a mode
aligned = kalign.align(sequences, mode="fast", ensemble=5)

# Ensemble alignment with per-residue confidence scores
result = kalign.align_from_file("proteins.fasta", mode="precise")
print(result.column_confidence)    # per-column confidence [0..1]
print(result.residue_confidence)   # per-residue confidence [0..1]
```

The Python package also provides a `kalign-py` CLI that uses the Python bindings directly (does not require the C binary):

```bash
kalign-py -i sequences.fasta -o aligned.fasta --format fasta
```

For comprehensive Python documentation, see [README-python.md](README-python.md) and the [python-docs directory](python-docs/).

## Examples

### Basic Usage

**Pass sequences via stdin**:
```bash
cat input.fa | kalign -f fasta > out.afa
```

**Combine multiple input files**:
```bash
kalign seqsA.fa seqsB.fa seqsC.fa -f fasta > combined.afa
```

**Use optimal threading** (auto-detected):
```bash
kalign -i sequences.fa -o aligned.afa  # Uses N-1 threads automatically
```

**Custom threading**:
```bash
kalign -i sequences.fa -o aligned.afa -n 8  # Use exactly 8 threads
```

### Ensemble Alignment

Precise mode runs 3 ensemble alignments with realignment, producing the highest precision. For more control, use `--ensemble` directly:

```bash
# Precise mode (ensemble + realign, highest precision)
kalign --precise -i sequences.fa -o aligned.afa

# Or fine-tune: 5 ensemble runs instead of 3
kalign -i sequences.fa -o aligned.afa --ensemble 5

# Save the POAR consensus table for later re-thresholding
kalign -i sequences.fa -o aligned.afa --ensemble 5 --save-poar consensus.poar

# Instantly re-threshold with a different support level (no re-alignment needed)
kalign -i sequences.fa -o rethresholded.afa --load-poar consensus.poar --min-support 3
```

### Format Conversion

**MSF format**:
```bash
kalign -i BB11001.tfa -f msf -o out.msf
```

**Clustal format**:
```bash
kalign -i BB11001.tfa -f clu -o out.clu
```

**Re-align existing alignment**:
```bash
kalign -i BB11001.msf -o out.afa
```

## Library Integration

### CMake Integration

Link Kalign into your C/C++ projects:

```cmake
find_package(kalign)
target_link_libraries(<target> kalign::kalign)
```

**Direct inclusion**:
```cmake
if (NOT TARGET kalign)
  add_subdirectory(<path_to_kalign>/kalign EXCLUDE_FROM_ALL)
endif ()
target_link_libraries(<target> kalign::kalign)
```

### Python Module Development

**Local development**:
```bash
uv pip install -e .
```

**Build Python module with CMake**:
```bash
mkdir build && cd build
cmake -DBUILD_PYTHON_MODULE=ON ..
make
```

### Cutting a New Python Release (PyPI)

This repo is set up to publish to PyPI from GitHub Actions on version tags (`v*`) via `.github/workflows/wheels.yml`.

1) Bump versions
- Update the Python package version in `pyproject.toml` (`[project].version`).
- Update the C/C++ library version in `CMakeLists.txt` (`KALIGN_LIBRARY_VERSION_{MAJOR,MINOR,PATCH}`).
- (Optional but recommended) Add an entry to `ChangeLog`.

2) Sanity check locally (recommended)
```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -U pip build twine
python -m build
twine check dist/*
```

3) Configure trusted publishing on PyPI
- Add this GitHub repo/workflow as a trusted publisher in the PyPI project settings.
- Configure it to match the workflow/environment used in `.github/workflows/wheels.yml` (environment `pypi`).

4) Tag and push
```bash
git tag vX.Y.Z
git push origin vX.Y.Z
```

That tag push triggers the wheel + sdist build, install tests, and then uploads to PyPI.

#### Testing publishing on TestPyPI (optional)

If you want to test the publishing pipeline without uploading to the real PyPI project, this repo also supports a manual TestPyPI publish:

1) Create a TestPyPI API token and add it as the GitHub secret `TEST_PYPI_API_TOKEN`.
2) Run the GitHub Actions workflow **Build Python Wheels** manually and set `publish_target = testpypi`.

To install from TestPyPI while resolving dependencies from PyPI:
```bash
pip install -i https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple <package-name>
```

## Benchmarks

The repository includes an automated benchmark suite that scores kalign against the BAliBASE reference alignments using both the C binary and the Python API. Results are tracked across commits via GitHub Actions.

```bash
# Run benchmarks locally (requires: pip install -e . and a CMake build)
make benchmark                      # full BAliBASE suite
make benchmark BENCH_MAX_CASES=5    # quick smoke test

# Or via Python directly
python -m benchmarks --dataset balibase --method python_api cli -v
```

### Historical Results

Kalign performs well for both speed and accuracy:

#### Balibase
![Balibase_scores](https://user-images.githubusercontent.com/8110320/198513840-0e08a634-bb41-4826-bd58-7fc66eae1054.jpeg)

#### Bralibase  
![Bralibase_scores](https://user-images.githubusercontent.com/8110320/198513850-00e5037f-355f-45ec-828f-ed8d47497272.jpeg)

### Performance Features

- **Multi-threading**: Automatic CPU core detection with OpenMP parallelization
- **SIMD optimizations**: Vectorized algorithms on x86_64 systems (SSE4.1, AVX, AVX2)
- **Bit-parallel algorithms**: Myers' algorithm for efficient alignment
- **Memory optimization**: Custom allocation strategies for large datasets

### Performance Tips

- **Let auto-threading work**: The default N-1 threading usually provides good performance
- **Large datasets**: Consider using `--type internal` for sequences with many gaps
- **Memory**: For very large alignments, monitor memory usage and consider reducing thread count
- **x86_64 systems**: SIMD optimizations provide additional speedup on Intel/AMD processors

## Contributing

We welcome contributions! See our [Contributing Guide](CONTRIBUTING.md) for details on:

- Reporting bugs and requesting features
- Development environment setup  
- Code style guidelines
- Pull request process

## Community Standards

This project follows the [Contributor Covenant Code of Conduct](CODE_OF_CONDUCT.md). By participating, you agree to uphold this code.

## System Requirements

- **Linux**: GCC 4.8+ or Clang 3.4+
- **macOS**: Xcode 8+ or Homebrew GCC/Clang
- **Memory**: ~1GB RAM per 10,000 sequences (typical)
- **CPU**: Any modern processor (additional optimizations on x86_64)

## Troubleshooting

### Common Issues

**macOS OpenMP**: If you see OpenMP-related errors on macOS:
```bash
brew install libomp
# Kalign automatically finds Homebrew's OpenMP installation
```

**Python module**: For Python installation issues:
```bash
pip install --upgrade pip setuptools wheel
pip install git+https://github.com/TimoLassmann/kalign.git
```

**Threading**: If performance seems slow, check thread detection:
```bash
kalign --help  # Shows current thread default
kalign -i test.fa -n 1 -o out.fa  # Force single-threaded for testing
```

For more troubleshooting, see [python-docs/python-troubleshooting.md](python-docs/python-troubleshooting.md).

## Citation

Please cite Kalign in your publications:

1. **Lassmann, Timo.** *Kalign 3: multiple sequence alignment of large data sets.* **Bioinformatics** (2019). [DOI](https://doi.org/10.1093/bioinformatics/btz795) | [PDF](https://academic.oup.com/bioinformatics/advance-article-pdf/doi/10.1093/bioinformatics/btz795/30314127/btz795.pdf)

### Previous Versions

2. **Lassmann, Timo, Oliver Frings, and Erik LL Sonnhammer.** *Kalign2: high-performance multiple alignment of protein and nucleotide sequences allowing external features.* **Nucleic acids research** 37.3 (2008): 858-865. [PubMed](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2647288/)

3. **Lassmann, Timo, and Erik LL Sonnhammer.** *Kalign: an accurate and fast multiple sequence alignment algorithm.* **BMC bioinformatics** 6.1 (2005): 298. [PubMed](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1325270/)

## License

Kalign is licensed under the Apache License, Version 2.0. See [COPYING](COPYING) for details.
