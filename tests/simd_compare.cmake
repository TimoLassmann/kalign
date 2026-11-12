# The plain kernels and the vector kernels must give the same alignment.
#
# Two ways to call it:
#  1. with -DGOLDEN_FILE=<preset_golden.cmake>: runs every golden case with
#     KALIGN_SIMD=scalar and requires the golden SHA-256 (so the golden
#     values hold for the plain path, too).
#  2. with -DINPUT=<fasta> -DMODE=<fast|accurate> [-DTHREADS=n]: runs the
#     input normally and with KALIGN_SIMD=scalar and requires identical
#     output.
#     Instead of INPUT, -DNSEQ=<n> -DALPHABET=<protein|dna> generates a
#     deterministic input of n short sequences (gen_many_seqs.cmake).
#     With -DEXPECT=<sha256> both outputs must also have this SHA-256.
# Other arguments: KALIGN (executable), OUT (scratch directory).
# DATA (tests/data) is needed for 1.
if(NOT DEFINED THREADS)
  set(THREADS 1)
endif()
file(MAKE_DIRECTORY "${OUT}")

function(run_kalign scalar mode threads input outfile)
  file(REMOVE "${outfile}")
  if(scalar)
    set(ENV{KALIGN_SIMD} "scalar")
  else()
    unset(ENV{KALIGN_SIMD})
  endif()
  execute_process(
    COMMAND "${KALIGN}" --mode ${mode} -n ${threads} -i "${input}" -o "${outfile}" -f fasta
    RESULT_VARIABLE rc
    OUTPUT_QUIET ERROR_VARIABLE se)
  unset(ENV{KALIGN_SIMD})
  if(NOT rc EQUAL 0 OR NOT EXISTS "${outfile}")
    message(FATAL_ERROR "kalign --mode ${mode} -n ${threads} on ${input} failed (scalar=${scalar}, rc=${rc})\n${se}")
  endif()
endfunction()

if(DEFINED GOLDEN_FILE)
  # reuse the GOLDEN list of preset_golden.cmake (only its set(GOLDEN ...) part)
  file(READ "${GOLDEN_FILE}" text)
  string(REGEX MATCH "set\\(GOLDEN[^)]*\\)" golden_set "${text}")
  if(golden_set STREQUAL "")
    message(FATAL_ERROR "no GOLDEN list found in ${GOLDEN_FILE}")
  endif()
  cmake_language(EVAL CODE "${golden_set}")
  set(n_fail 0)
  foreach(entry IN LISTS GOLDEN)
    string(REPLACE ":" ";" parts "${entry}")
    list(GET parts 0 mode)
    list(GET parts 1 input)
    list(GET parts 2 expected)
    set(out_file "${OUT}/scalar_${mode}_${input}.fasta")
    run_kalign(1 ${mode} 1 "${DATA}/${input}" "${out_file}")
    file(SHA256 "${out_file}" got)
    if(got STREQUAL expected)
      message(STATUS "ok   scalar ${mode} ${input}")
    else()
      message(SEND_ERROR "FAIL scalar ${mode} ${input}: sha256 ${got}, expected ${expected}")
      math(EXPR n_fail "${n_fail} + 1")
    endif()
  endforeach()
  if(n_fail GREATER 0)
    message(FATAL_ERROR "${n_fail} scalar golden check(s) failed")
  endif()
else()
  if(NOT DEFINED MODE)
    message(FATAL_ERROR "simd_compare.cmake: -DMODE=... is required")
  endif()
  if(DEFINED NSEQ)
    set(INPUT "${OUT}/many_${ALPHABET}_${NSEQ}.fa")
    execute_process(
      COMMAND ${CMAKE_COMMAND} -DOUTFILE=${INPUT} -DNSEQ=${NSEQ} -DALPHABET=${ALPHABET}
              -P "${CMAKE_CURRENT_LIST_DIR}/gen_many_seqs.cmake"
      RESULT_VARIABLE rc ERROR_VARIABLE se)
    if(NOT rc EQUAL 0)
      message(FATAL_ERROR "input generator failed\n${se}")
    endif()
  elseif(NOT DEFINED INPUT)
    message(FATAL_ERROR "simd_compare.cmake: -DINPUT=... or -DNSEQ=... is required")
  endif()
  get_filename_component(base "${INPUT}" NAME)
  set(o_simd "${OUT}/simd_${MODE}_${THREADS}_${base}")
  set(o_scalar "${OUT}/scalar_${MODE}_${THREADS}_${base}")
  run_kalign(0 ${MODE} ${THREADS} "${INPUT}" "${o_simd}")
  run_kalign(1 ${MODE} ${THREADS} "${INPUT}" "${o_scalar}")
  file(SHA256 "${o_simd}" h1)
  file(SHA256 "${o_scalar}" h2)
  if(NOT h1 STREQUAL h2)
    message(FATAL_ERROR "FAIL ${MODE} -n ${THREADS} ${base}: vector output ${h1} differs from plain output ${h2}")
  endif()
  if(DEFINED EXPECT AND NOT h1 STREQUAL EXPECT)
    message(FATAL_ERROR "FAIL ${MODE} -n ${THREADS} ${base}: output ${h1}, expected ${EXPECT}")
  endif()
  message(STATUS "ok   ${MODE} -n ${THREADS} ${base}: vector and plain outputs are identical (${h1})")
  file(REMOVE "${o_simd}" "${o_scalar}")
endif()
