# Kalign Benchmark Container
#
# Includes kalign, Clustal Omega, MAFFT, and MUSCLE v5 for comparative
# benchmarking on BAliBASE, BRAliBASE, and BaliFam100 datasets.
#
# Build (installs kalign from the 'main' branch):
#   podman build -t kalign-benchmark .
#
# Build from a specific commit for reproducibility:
#   podman build --build-arg KALIGN_REF=abc1234 -t kalign-benchmark .
#
# Run benchmarks:
#   podman run -it \
#     -v ./benchmarks/data:/kalign/benchmarks/data \
#     kalign-benchmark \
#     python -m benchmarks \
#       --dataset balibase --method cli clustalo mafft muscle \
#       --mode fast default recall accurate -v
#
# Run with no arguments to print the benchmark help:
#   podman run -it kalign-benchmark

FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

# System dependencies + alignment tools (clustalo, mafft from apt)
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential cmake g++ git curl \
    python3 python3-pip python3-venv python3-dev \
    clustalo mafft \
    pkg-config \
    && rm -rf /var/lib/apt/lists/*

# ---------- Build MUSCLE v5 from source ----------
RUN cd /tmp && \
    git clone --depth 1 https://github.com/rcedgar/muscle.git && \
    cd muscle/src && \
    sed -i 's/defined(__arm64__)/defined(__arm64__) || defined(__aarch64__)/' myutils.h && \
    bash build_linux.bash && \
    cp ../bin/muscle /usr/local/bin/ && \
    rm -rf /tmp/muscle

# ---------- Python environment ----------
RUN python3 -m venv /venv
ENV PATH="/venv/bin:$PATH"
RUN pip install --no-cache-dir uv

# ---------- Clone kalign and build ----------
ARG KALIGN_REF=main
RUN git clone --branch ${KALIGN_REF} --depth 1 \
    https://github.com/TimoLassmann/kalign.git /kalign
WORKDIR /kalign

# Build kalign C binary
RUN mkdir cbuild && cd cbuild && \
    cmake -DCMAKE_BUILD_TYPE=Release .. && \
    make -j"$(nproc)" && \
    cp src/kalign /usr/local/bin/kalign

# Install Python package with benchmark dependencies
RUN uv pip install --no-cache -e ".[benchmark]"

# ---------- Verify all tools ----------
RUN kalign --version && clustalo --version && mafft --version && muscle -version

# ---------- Data directory (mount point for BAliBASE etc.) ----------
RUN mkdir -p /kalign/benchmarks/data/downloads /kalign/benchmarks/results

CMD ["python", "-m", "benchmarks", "--help"]
