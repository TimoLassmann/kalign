# Kalign Benchmark Container
#
# Includes kalign, Clustal Omega, MAFFT, and MUSCLE v5 for comparative
# benchmarking on BAliBASE, BRAliBASE, and BaliFam100 datasets.
#
# Build:
#   podman build -t kalign-benchmark .
#
# Run the interactive dashboard:
#   podman run -it -p 8050:8050 \
#     -v ./benchmarks/data:/kalign/benchmarks/data \
#     kalign-benchmark
#
# Run a CLI benchmark:
#   podman run -it \
#     -v ./benchmarks/data:/kalign/benchmarks/data \
#     kalign-benchmark \
#     python -m benchmarks \
#       --dataset balibase --method python_api clustalo mafft muscle -v
#
# View results in the dashboard after a CLI run:
#   podman run -it -p 8050:8050 \
#     -v ./benchmarks/data:/kalign/benchmarks/data \
#     -v ./benchmarks/results:/kalign/benchmarks/results \
#     kalign-benchmark

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
# myutils.h checks __arm64__ (macOS) but not __aarch64__ (Linux); add it
RUN cd /tmp && \
    git clone --depth 1 https://github.com/rcedgar/muscle.git && \
    cd muscle/src && \
    sed -i 's/defined(__arm64__)/defined(__arm64__) || defined(__aarch64__)/' myutils.h && \
    bash build_linux.bash && \
    cp ../bin/muscle /usr/local/bin/ && \
    rm -rf /tmp/muscle

# ---------- Copy kalign source and build ----------
COPY . /kalign
WORKDIR /kalign

RUN mkdir -p build && cd build && \
    cmake -DCMAKE_BUILD_TYPE=Release .. && \
    make -j"$(nproc)"

# ---------- Python environment ----------
RUN python3 -m venv /venv
ENV PATH="/venv/bin:/kalign/build/src:$PATH"

RUN pip install --no-cache-dir uv && \
    uv pip install --no-cache -e ".[benchmark]"

# ---------- Verify tools ----------
RUN which kalign && which clustalo && which mafft && which muscle

# ---------- Data & results directories ----------
RUN mkdir -p /kalign/benchmarks/data/downloads /kalign/benchmarks/results

EXPOSE 8050

CMD ["python", "-m", "benchmarks.app", "--host", "0.0.0.0", "--port", "8050"]
