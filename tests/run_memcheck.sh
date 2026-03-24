#!/bin/bash
# run_memcheck.sh — Run memory checking tools on kalign
#
# Usage:
#   ./run_memcheck.sh asan       # ASAN + leak detection (C tests)
#   ./run_memcheck.sh valgrind   # Valgrind (C tests)
#   ./run_memcheck.sh python     # Python stress test (normal build)
#   ./run_memcheck.sh all        # All of the above

set -e

MODE="${1:-all}"
N_ITERS="${2:-10}"

TESTDATA="/kalign/tests/data"
INPUT="$TESTDATA/BB11001.tfa"
REF="$TESTDATA/BB11001.msf"
INPUT2="$TESTDATA/BB30014.tfa"
REF2="$TESTDATA/BB30014.msf"

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

banner() {
    echo ""
    echo -e "${YELLOW}============================================================${NC}"
    echo -e "${YELLOW}  $1${NC}"
    echo -e "${YELLOW}============================================================${NC}"
    echo ""
}

# ---- ASAN stress test ----
run_asan() {
    banner "ASAN + Leak Detection (C stress test, $N_ITERS iters)"

    export ASAN_OPTIONS="detect_leaks=1:halt_on_error=0:print_stats=1"

    echo "--- Testing with BB11001 ---"
    /kalign/build-asan/memcheck_stress "$INPUT" "$REF" "$N_ITERS"

    echo ""
    echo "--- Testing with BB30014 ---"
    /kalign/build-asan/memcheck_stress "$INPUT2" "$REF2" "$N_ITERS"

    echo ""
    echo "--- ASAN built-in tests ---"
    cd /kalign/build-asan && ctest --output-on-failure

    echo -e "\n${GREEN}ASAN tests complete.${NC}"
}

# ---- Valgrind stress test ----
run_valgrind() {
    banner "Valgrind Memcheck (C stress test, $N_ITERS iters)"

    echo "--- Testing with BB11001 ---"
    valgrind --leak-check=full \
             --show-leak-kinds=all \
             --track-origins=yes \
             --error-exitcode=1 \
             --errors-for-leak-kinds=definite \
             /kalign/build-debug/memcheck_stress "$INPUT" "$REF" "$N_ITERS"

    echo ""
    echo "--- Testing with BB30014 ---"
    valgrind --leak-check=full \
             --show-leak-kinds=all \
             --track-origins=yes \
             --error-exitcode=1 \
             --errors-for-leak-kinds=definite \
             /kalign/build-debug/memcheck_stress "$INPUT2" "$REF2" "$N_ITERS"

    echo -e "\n${GREEN}Valgrind tests complete.${NC}"
}

# ---- Python stress test ----
run_python() {
    banner "Python Memory Stress Test ($N_ITERS iters)"

    python /kalign/tests/memcheck_stress.py "$N_ITERS"

    echo -e "\n${GREEN}Python stress tests complete.${NC}"
}

# ---- Main ----
case "$MODE" in
    asan)
        run_asan
        ;;
    valgrind)
        run_valgrind
        ;;
    python)
        run_python
        ;;
    all)
        run_asan
        run_valgrind
        run_python
        banner "ALL CHECKS COMPLETE"
        ;;
    *)
        echo "Usage: $0 {asan|valgrind|python|all} [n_iters]"
        exit 1
        ;;
esac
