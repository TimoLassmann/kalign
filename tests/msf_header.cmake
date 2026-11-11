# MSF output header: protein input is labelled AA / Type P, nucleotide input
# NA / Type N.
# Output format names are matched exactly; anything else is rejected.
#
# Invoked by ctest as:
#   cmake -DKALIGN=<kalign binary> -DDATA=<tests/data> -DOUT=<scratch dir> -P msf_header.cmake

foreach(v KALIGN DATA OUT)
  if(NOT DEFINED ${v})
    message(FATAL_ERROR "msf_header.cmake: -D${v}=... is required")
  endif()
endforeach()
file(MAKE_DIRECTORY "${OUT}")

macro(fail msg)
  message(SEND_ERROR "FAIL ${msg}")
endmacro()

function(check_msf file kind letter)
  file(READ "${file}" txt)
  if(NOT txt MATCHES "^!!${kind}_MULTIPLE_ALIGNMENT 1.0")
    fail("${file}: first line is not !!${kind}_MULTIPLE_ALIGNMENT 1.0")
  endif()
  if(NOT txt MATCHES "MSF: [0-9]+  Type: ${letter}  ")
    fail("${file}: header line does not have Type: ${letter}")
  endif()
  message(STATUS "ok   ${file} ${kind} / ${letter}")
endfunction()

set(prot "${DATA}/BB11001.tfa")
set(dna "${DATA}/small.fa")

# 1. alignment written as MSF (fast mode)
execute_process(COMMAND "${KALIGN}" -i "${prot}" -o "${OUT}/prot.msf" -f msf --mode fast -n 1 -q
                RESULT_VARIABLE rc OUTPUT_QUIET ERROR_QUIET)
if(NOT rc EQUAL 0)
  fail("protein alignment to MSF failed (${rc})")
else()
  check_msf("${OUT}/prot.msf" AA P)
endif()

execute_process(COMMAND "${KALIGN}" -i "${dna}" -o "${OUT}/dna.msf" -f msf --mode fast -n 1 -q
                RESULT_VARIABLE rc OUTPUT_QUIET ERROR_QUIET)
if(NOT rc EQUAL 0)
  fail("DNA alignment to MSF failed (${rc})")
else()
  check_msf("${OUT}/dna.msf" NA N)
endif()

# 2. format names: exact match only
foreach(f fastq msff cluster xfa)
  execute_process(COMMAND "${KALIGN}" -i "${dna}" -o "${OUT}/bad.out" -f ${f} --mode fast -n 1 -q
                  RESULT_VARIABLE rc OUTPUT_QUIET ERROR_VARIABLE err)
  if(rc EQUAL 0)
    fail("format '${f}' was accepted")
  else()
    message(STATUS "ok   format '${f}' rejected")
  endif()
endforeach()
# Upper-case names are accepted by the library (msa_io.c); the command line
# pre-check in src/parameters.c is covered by another test.
foreach(f fasta fa afa msf clu clustal)
  execute_process(COMMAND "${KALIGN}" -i "${dna}" -o "${OUT}/good.out" -f ${f} --mode fast -n 1 -q
                  RESULT_VARIABLE rc OUTPUT_QUIET ERROR_QUIET)
  if(NOT rc EQUAL 0)
    fail("format '${f}' was rejected")
  else()
    message(STATUS "ok   format '${f}' accepted")
  endif()
endforeach()

