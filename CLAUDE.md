# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Kalign is a fast multiple sequence alignment (MSA) program for biological sequences (protein, DNA, RNA) written in C. The project uses CMake as the primary build system with Zig as an alternative for cross-compilation.

## Build Commands

### CMake Build (Primary)
```bash
mkdir build
cd build
cmake ..
make
make test
make install
```

### Debug Build
```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Debug ..
make
```

### Address Sanitizer Build
```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=ASAN ..
make
```

### Build Without OpenMP
```bash
mkdir build
cd build
cmake -DUSE_OPENMP=OFF ..
make
```

### Zig Build (Alternative)
```bash
zig build
```

## Testing

### Run All Tests
```bash
make test
# or
ctest
```

### Run Specific Test
```bash
ctest -R <test_name>
```

### Test Executables
- `kalign_io_test` - I/O functionality tests
- `kalign_lib_test` - Library API tests  
- `kalign_cmp_test` - Comparison tests
- `kaligncpp` - C++ API tests
- `dssim` - Distance similarity tests

## Code Architecture

### Directory Structure
- `lib/` - Core library implementation
  - `src/` - Library source files (alignment algorithms, I/O, utilities)
  - `include/kalign/` - Public API headers
- `src/` - Main executable (command-line interface)
- `tests/` - Test suite with biological sequence data

### Key Components

#### Alignment Engine (`lib/src/aln_*.c`)
- `aln_seqseq.c` - Sequence-to-sequence alignment
- `aln_seqprofile.c` - Sequence-to-profile alignment  
- `aln_profileprofile.c` - Profile-to-profile alignment
- `aln_controller.c` - Alignment orchestration
- `bpm.c` - Bit-parallel matching (Myers' algorithm)

#### I/O System (`lib/src/msa_*.c`)
- `msa_io.c` - Multi-format sequence reading/writing (FASTA, MSF, Clustal)
- `msa_alloc.c` - Memory management for sequence data
- `msa_op.c` - Sequence operations and manipulations

#### Distance Calculations
- `sequence_distance.c` - Sequence distance metrics
- `euclidean_dist.c` - Euclidean distance calculations
- `bisectingKmeans.c` - K-means clustering for guide tree construction

#### Utilities
- `tldevel.c` - Development utilities and debugging
- `tlmisc.c` - Miscellaneous helper functions
- `alphabet.c` - Sequence alphabet handling (DNA, RNA, protein)
- `task.c` - Task scheduling and threading support

### Performance Features
- **Multi-threading**: OpenMP parallelization
- **SIMD**: SSE4.1, AVX, AVX2 optimizations (enabled by default)
- **Bit-parallel algorithms**: Myers' algorithm for efficient alignment
- **Memory optimization**: Custom allocation strategies

### API Usage
The library provides a C API (`kalign.h`) with C++ compatibility. Key functions:
- Reading sequences from files
- Running alignments with different parameters
- Writing results in various formats
- Memory management helpers

### Build Configuration
- C11 standard required
- OpenMP support (can be disabled with `-DUSE_OPENMP=OFF`)
- SIMD instruction support with runtime detection
- Multiple build types: Release, Debug, ASAN (Address Sanitizer)
- Cross-compilation support via Zig build system

### Testing Strategy
- Unit tests for core components (BPM, distance calculations, I/O)
- Integration tests using Balibase sequence datasets
- C++ API compatibility tests
- Performance benchmarks
- Format conversion tests (FASTA, MSF, Clustal)

## Python Module

The repository includes a Python package that provides Python bindings for the Kalign library. The Python package is configured at the root level using modern Python packaging standards.

### Python Development Commands

#### Install in Development Mode (Recommended)
```bash
uv pip install -e .
```

#### Build Python Package Locally
```bash
uv run python -m build
```

#### Build with CMake (for debugging C extensions)
```bash
mkdir build
cd build
cmake -DBUILD_PYTHON_MODULE=ON ..
make
```

#### Test Python Package
```bash
uv run python -c "import kalign; print(kalign.__version__); seqs=['ATCG','ATCGG']; print(kalign.align(seqs))"
```

#### Run Python Tests
```bash
uv run pytest tests/python/ -v
```

#### Build Wheels for Distribution
```bash
uv run python -m cibuildwheel --output-dir wheelhouse
```

### Python Package Structure
- `python-kalign/__init__.py` - High-level Python API
- `python-kalign/_core.cpp` - pybind11 C++ bindings
- `pyproject.toml` - Modern Python build configuration (at root level)
- `README-python.md` - Python-specific documentation
- `tests/python/` - Python test suite

### Python API Features
- **Simple interface**: `kalign.align(sequences, seq_type="auto")`
- **File support**: `kalign.align_from_file("file.fasta")`
- **Parameter control**: Gap penalties, threading, sequence types
- **Error handling**: Python exceptions for C library errors
- **Memory safety**: Automatic memory management through pybind11

### GitHub Actions Integration
The `.github/workflows/wheels.yml` workflow:
- Builds wheels for Linux, macOS, Windows
- Tests multiple Python versions (3.9-3.13)
- Uses cibuildwheel for cross-platform compatibility
- Publishes to PyPI on tagged releases