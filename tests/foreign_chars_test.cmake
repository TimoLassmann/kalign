# Forcing --type dna on protein input: letters outside the DNA/RNA alphabet
# must be reported in exactly one summary line, not one warning per
# character (and once per ensemble run in accurate mode).
# Arguments: KALIGN (executable), INPUT (protein FASTA), WORKDIR, MODE
file(MAKE_DIRECTORY "${WORKDIR}")
set(out "${WORKDIR}/foreign_${MODE}.fa")
execute_process(
  COMMAND "${KALIGN}" --mode ${MODE} --type dna -n 1 -i "${INPUT}" -o "${out}"
  RESULT_VARIABLE rc
  OUTPUT_VARIABLE so
  ERROR_VARIABLE se)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "kalign --mode ${MODE} --type dna failed (rc=${rc})\n${se}")
endif()
set(all "${so}${se}")
string(REGEX MATCHALL "not in the DNA/RNA alphabet" summary "${all}")
list(LENGTH summary nsummary)
string(REGEX MATCHALL "offending character" old "${all}")
list(LENGTH old nold)
if(NOT nsummary EQUAL 1 OR NOT nold EQUAL 0)
  message(FATAL_ERROR "expected 1 summary line and no per-character warnings, got ${nsummary} and ${nold}\n${all}")
endif()
file(REMOVE "${out}")
