#!/bin/bash
# Setup script for running kalign benchmarks on a Linux server.
# Usage: bash setup_server.sh
set -euo pipefail

echo "=== Kalign benchmark server setup ==="

# Check basics
echo "Checking prerequisites..."
command -v cmake >/dev/null 2>&1 || { echo "ERROR: cmake not found. Install with: sudo apt install cmake (or module load cmake)"; exit 1; }
command -v gcc >/dev/null 2>&1 || command -v cc >/dev/null 2>&1 || { echo "ERROR: C compiler not found."; exit 1; }

# Install uv if not present
if ! command -v uv >/dev/null 2>&1; then
    echo "Installing uv..."
    curl -LsSf https://astral.sh/uv/install.sh | sh
    export PATH="$HOME/.local/bin:$PATH"
fi

echo "uv version: $(uv --version)"

# Build C library + Python extension
echo ""
echo "=== Building kalign Python package ==="
uv pip install -e ".[dev]" 2>/dev/null || uv pip install -e .

# Install benchmark dependencies
echo ""
echo "=== Installing benchmark dependencies ==="
uv pip install -e ".[benchmark]"

# Quick smoke test
echo ""
echo "=== Smoke test ==="
uv run python -c "
import kalign
print(f'kalign version: {kalign.__version__}')
result = kalign.align(['ACDEFGHIK', 'ACDGHIK', 'ACDEFHIK'])
print(f'Alignment test: OK ({len(result)} sequences)')
"

# Show system info
echo ""
echo "=== System info ==="
echo "CPU: $(nproc) threads available"
echo "RAM: $(free -h 2>/dev/null | awk '/^Mem:/{print $2}' || echo 'unknown')"
echo "Python: $(uv run python --version)"

echo ""
echo "=== Ready! ==="
echo ""
echo "Run the benchmark suite (see: uv run python -m benchmarks --help):"
echo ""
echo "  uv run python -m benchmarks --dataset balibase --method python_api --mode fast accurate -v"
