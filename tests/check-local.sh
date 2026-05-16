#!/usr/bin/env bash
#
# tests/check-local.sh — pre-push verification gate for kalign.
#
# Runs four independent test layers; each catches a different class of
# bug.  If every required phase passes, the working tree is safe to push.
#
# Usage:
#   tests/check-local.sh           # run all four phases (~5 min)
#   tests/check-local.sh --quick   # skip the Linux ASAN container (~30s)
#   tests/check-local.sh --help
#
# Phases:
#   1. zig build      — cross-compile sanity across aarch64-macos,
#                       aarch64-linux, x86_64-linux-{gnu,musl}.
#                       Catches GCC-vs-Clang divergence (e.g. an unused
#                       variable that is in fact used inside an #ifdef).
#   2. cmake + ctest  — native macOS Release build + 15 ctests.
#                       Algorithmic correctness, ABI, integration.
#   3. podman ASAN    — Ubuntu container with kalign built under ASAN
#                       and the full ctest suite. Catches Linux glibc
#                       behaviour that Apple's malloc hides (e.g.
#                       uninitialised MMALLOC fields that read as zero
#                       on macOS but garbage on Linux).
#   4. pytest         — Python bindings, mode presets, ecosystem
#                       integration. ~171 tests.
#
# Exits 0 only if every non-skipped phase passes.

set -u
set -o pipefail

cd "$(dirname "$0")/.."   # repo root

# ---- options ----------------------------------------------------------
QUICK=0
for arg in "$@"; do
    case "$arg" in
        --quick) QUICK=1 ;;
        -h|--help)
            sed -n '1,30p' "$0" | sed -n 's/^# \{0,1\}//p'
            exit 0
            ;;
        *)
            printf 'Unknown option: %s\n' "$arg" >&2
            printf 'Try: %s --help\n' "$0" >&2
            exit 2
            ;;
    esac
done

# ---- output helpers ---------------------------------------------------
if [ -t 1 ]; then
    YELLOW=$'\033[1;33m'
    GREEN=$'\033[0;32m'
    RED=$'\033[0;31m'
    NC=$'\033[0m'
else
    YELLOW=''
    GREEN=''
    RED=''
    NC=''
fi

banner() {
    printf '\n%s═══════════════════════════════════════════════════════════════%s\n' "$YELLOW" "$NC"
    printf '%s  %s%s\n' "$YELLOW" "$1" "$NC"
    printf '%s═══════════════════════════════════════════════════════════════%s\n\n' "$YELLOW" "$NC"
}

# ---- result tracking --------------------------------------------------
PASS_LIST=()
FAIL_LIST=()
SKIP_LIST=()

record_pass() { PASS_LIST+=("$1"); }
record_fail() { FAIL_LIST+=("$1"); }
record_skip() { SKIP_LIST+=("$1 — $2"); }

# ---- Phase 1: zig build ----------------------------------------------
phase_zig() {
    banner "Phase 1/4: zig build (cross-compile sanity)"
    if ! command -v zig >/dev/null 2>&1; then
        printf '%s⊘ zig not installed — skipping%s\n' "$YELLOW" "$NC"
        record_skip "Phase 1 (zig build)" "zig not installed"
        return 0
    fi
    if zig build; then
        record_pass "Phase 1 (zig build)"
        return 0
    fi
    record_fail "Phase 1 (zig build)"
    return 1
}

# ---- Phase 2: native CMake + ctest -----------------------------------
phase_cmake() {
    banner "Phase 2/4: native CMake build + ctest"
    local ncpu
    ncpu="$(sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 4)"
    mkdir -p build
    if ( cd build \
         && cmake .. -DCMAKE_BUILD_TYPE=Release > /dev/null \
         && make -j"$ncpu" \
         && ctest --output-on-failure ); then
        record_pass "Phase 2 (cmake + ctest)"
        return 0
    fi
    record_fail "Phase 2 (cmake + ctest)"
    return 1
}

# ---- Phase 3: Linux ASAN container -----------------------------------
phase_memcheck() {
    banner "Phase 3/4: Linux ASAN ctest (podman memcheck container)"
    if [ "$QUICK" = "1" ]; then
        printf '%s⊘ skipped (--quick)%s\n' "$YELLOW" "$NC"
        record_skip "Phase 3 (Linux ASAN)" "--quick"
        return 0
    fi
    if ! command -v podman >/dev/null 2>&1; then
        printf '%s⊘ podman not installed — skipping%s\n' "$YELLOW" "$NC"
        record_skip "Phase 3 (Linux ASAN)" "podman not installed"
        return 0
    fi
    if ! podman info >/dev/null 2>&1; then
        printf '%s⊘ podman machine not running — try:  podman machine start%s\n' "$YELLOW" "$NC"
        record_skip "Phase 3 (Linux ASAN)" "podman machine not running"
        return 0
    fi
    if ! podman build -f Containerfile.memcheck -t kalign-memcheck . ; then
        record_fail "Phase 3 (Linux ASAN — container build)"
        return 1
    fi
    if podman run --rm kalign-memcheck bash -c \
        "cd /kalign/build-asan && ASAN_OPTIONS='detect_leaks=0:halt_on_error=1:abort_on_error=1' ctest --output-on-failure"; then
        record_pass "Phase 3 (Linux ASAN)"
        return 0
    fi
    record_fail "Phase 3 (Linux ASAN)"
    return 1
}

# ---- Phase 4: pytest -------------------------------------------------
phase_python() {
    banner "Phase 4/4: Python tests (pytest)"
    if ! command -v uv >/dev/null 2>&1; then
        printf '%s⊘ uv not installed — skipping%s\n' "$YELLOW" "$NC"
        record_skip "Phase 4 (pytest)" "uv not installed"
        return 0
    fi
    if ! uv pip install -e . \
            --config-settings cmake.args="-DUSE_OPENMP=OFF;-DUSE_THREADPOOL=ON" \
            --force-reinstall --no-deps --quiet; then
        record_fail "Phase 4 (pytest — install)"
        return 1
    fi
    if uv run pytest tests/python/ -q --no-header; then
        record_pass "Phase 4 (pytest)"
        return 0
    fi
    record_fail "Phase 4 (pytest)"
    return 1
}

# ---- run all phases (each independent; no short-circuit) -------------
phase_zig      || true
phase_cmake    || true
phase_memcheck || true
phase_python   || true

# ---- summary ---------------------------------------------------------
banner "Summary"
if [ "${#PASS_LIST[@]}" -gt 0 ]; then
    for item in "${PASS_LIST[@]}"; do
        printf '  %s✓%s %s\n' "$GREEN" "$NC" "$item"
    done
fi
if [ "${#SKIP_LIST[@]}" -gt 0 ]; then
    for item in "${SKIP_LIST[@]}"; do
        printf '  %s⊘%s %s\n' "$YELLOW" "$NC" "$item"
    done
fi
if [ "${#FAIL_LIST[@]}" -gt 0 ]; then
    for item in "${FAIL_LIST[@]}"; do
        printf '  %s✗%s %s\n' "$RED" "$NC" "$item"
    done
fi
echo

n_fail="${#FAIL_LIST[@]}"
if [ "$n_fail" -eq 0 ]; then
    printf '%sAll required phases passed. Safe to push.%s\n' "$GREEN" "$NC"
    exit 0
else
    printf '%s%d phase(s) failed. Do NOT push.%s\n' "$RED" "$n_fail" "$NC"
    exit 1
fi
