#!/bin/bash
# Setup script for running kalign parameter optimization on a Linux server.
# Usage: bash setup_server.sh
set -euo pipefail

echo "=== Kalign optimization server setup ==="

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

# Install optimizer dependencies
echo ""
echo "=== Installing optimizer dependencies ==="
uv pip install pymoo rich

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
echo "Example runs:"
echo ""
echo "  # Quick test (10 minutes):"
echo "  uv run python -m benchmarks.optimize_params --pop-size 20 --n-gen 5 --n-workers 8 --n-threads 1"
echo ""
echo "  # Single-run optimization (~2-3 hours):"
echo "  uv run python -m benchmarks.optimize_params --pop-size 60 --n-gen 50 --n-workers 56 --n-threads 1"
echo ""
echo "  # Resume after interrupt:"
echo "  uv run python -m benchmarks.optimize_params --resume benchmarks/results/optim/gen_checkpoint.pkl --n-gen 80 --n-workers 56 --n-threads 1"
