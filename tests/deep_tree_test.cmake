# Deep guide-tree regression test: N identical sequences give a guide tree
# that is thousands of levels deep.  Used to crash with "thread stack size
# exceeded" on macOS when the tree walk recursed once per tree level.
# Arguments: KALIGN (executable), WORKDIR, NSEQ, MODE, THREADS
file(MAKE_DIRECTORY "${WORKDIR}")
set(in "${WORKDIR}/deep_${MODE}_${NSEQ}.fa")
set(out "${WORKDIR}/deep_${MODE}_${NSEQ}.out.fa")
set(seq "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVV")
string(SUBSTRING "${seq}" 0 60 seq)
set(content "")
math(EXPR last "${NSEQ} - 1")
foreach(i RANGE 0 ${last})
  string(APPEND content ">s${i}\n${seq}\n")
endforeach()
file(WRITE "${in}" "${content}")
execute_process(
  COMMAND "${KALIGN}" --mode ${MODE} -n ${THREADS} -i "${in}" -o "${out}" -f fasta
  RESULT_VARIABLE rc
  OUTPUT_VARIABLE so
  ERROR_VARIABLE se)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "kalign ${MODE} -n ${THREADS} on ${NSEQ} sequences failed (rc=${rc})\n${se}")
endif()
file(STRINGS "${out}" lines)
list(LENGTH lines nlines)
math(EXPR expect "2 * ${NSEQ}")
if(NOT nlines EQUAL expect)
  message(FATAL_ERROR "expected ${expect} output lines, got ${nlines}")
endif()
file(REMOVE "${in}" "${out}")
