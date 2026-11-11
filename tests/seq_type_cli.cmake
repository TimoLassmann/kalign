# --type and the "accurate" size limit on the command line.
#
#  1. --type accepts auto/protein/dna/rna only (internal, divergent, ... fail).
#  2. --type matching the detected type changes nothing (byte-identical).
#  3. --type overriding detection warns and selects the other presets, also
#     in accurate (5-member ensemble) mode.
#  4. rna and dna select the same (nucleotide) presets.
#  5. accurate refuses more than 5,000 sequences for DNA as for protein.
#
# Invoked by ctest as:
#   cmake -DKALIGN=<kalign binary> -DDATA=<tests/data> -DOUT=<scratch dir> -P seq_type_cli.cmake

foreach(v KALIGN DATA OUT)
  if(NOT DEFINED ${v})
    message(FATAL_ERROR "seq_type_cli.cmake: -D${v}=... is required")
  endif()
endforeach()
file(MAKE_DIRECTORY "${OUT}")
set(n_fail 0)

macro(fail msg)
  message(SEND_ERROR "FAIL ${msg}")
  math(EXPR n_fail "${n_fail} + 1")
endmacro()

function(run_kalign out_var rc_var err_var)
  execute_process(COMMAND "${KALIGN}" ${ARGN}
                  RESULT_VARIABLE rc OUTPUT_QUIET ERROR_VARIABLE err)
  set(${rc_var} "${rc}" PARENT_SCOPE)
  set(${err_var} "${err}" PARENT_SCOPE)
endfunction()

set(dna "${DATA}/small.fa")
set(prot "${DATA}/BB11001.tfa")

# 1. removed / unknown names
foreach(t internal divergent pfasum60 nonsense)
  run_kalign(o rc err -i "${dna}" -o "${OUT}/x.fa" --type ${t} -n 1)
  if(rc EQUAL 0)
    fail("--type ${t} was accepted")
  else()
    message(STATUS "ok   --type ${t} rejected")
  endif()
endforeach()

# 2. matching type: identical to auto
foreach(case "fast;${dna};dna" "accurate;${dna};rna" "fast;${prot};protein" "accurate;${prot};protein")
  list(GET case 0 mode)
  list(GET case 1 in)
  list(GET case 2 t)
  run_kalign(o rc err -i "${in}" -o "${OUT}/auto.fa" --mode ${mode} -n 1 -q)
  run_kalign(o rc2 err2 -i "${in}" -o "${OUT}/typed.fa" --mode ${mode} -n 1 --type ${t})
  file(SHA256 "${OUT}/auto.fa" h1)
  file(SHA256 "${OUT}/typed.fa" h2)
  if(NOT rc EQUAL 0 OR NOT rc2 EQUAL 0 OR NOT h1 STREQUAL h2 OR err2 MATCHES "WARNING")
    fail("--type ${t} --mode ${mode} on matching input differs from auto or warns")
  else()
    message(STATUS "ok   --type ${t} --mode ${mode} matching input == auto")
  endif()
endforeach()

# 3. override: warns, differs from auto, works in accurate (ensemble) mode
foreach(case "fast;${dna};protein;protein" "accurate;${dna};protein;protein"
             "fast;${prot};dna;DNA/RNA" "accurate;${prot};dna;DNA/RNA")
  list(GET case 0 mode)
  list(GET case 1 in)
  list(GET case 2 t)
  list(GET case 3 as)
  run_kalign(o rc err -i "${in}" -o "${OUT}/auto.fa" --mode ${mode} -n 1 -q)
  run_kalign(o rc2 err2 -i "${in}" -o "${OUT}/typed.fa" --mode ${mode} -n 1 --type ${t})
  file(SHA256 "${OUT}/auto.fa" h1)
  file(SHA256 "${OUT}/typed.fa" h2)
  if(NOT rc2 EQUAL 0)
    fail("--type ${t} --mode ${mode} override exited with ${rc2}")
  elseif(h1 STREQUAL h2)
    fail("--type ${t} --mode ${mode} override gave the auto alignment")
  elseif(NOT err2 MATCHES "aligning as ${as} as requested")
    fail("--type ${t} --mode ${mode} override did not warn")
  else()
    message(STATUS "ok   --type ${t} --mode ${mode} overrides detection")
  endif()
endforeach()

# 4. rna == dna
run_kalign(o rc err -i "${prot}" -o "${OUT}/r.fa" --type rna -n 1 -q)
run_kalign(o rc err -i "${prot}" -o "${OUT}/d.fa" --type dna -n 1 -q)
file(SHA256 "${OUT}/r.fa" h1)
file(SHA256 "${OUT}/d.fa" h2)
if(NOT h1 STREQUAL h2)
  fail("--type rna and --type dna differ")
else()
  message(STATUS "ok   --type rna == --type dna")
endif()

# 5. accurate size limit, nucleotide and protein
foreach(case "dna;ACGTACGTAC" "protein;MKTAYIAKQR")
  list(GET case 0 label)
  list(GET case 1 seq)
  set(big "")
  foreach(i RANGE 1 5001)
    string(APPEND big ">s${i}\n${seq}\n")
  endforeach()
  file(WRITE "${OUT}/big_${label}.fa" "${big}")
  run_kalign(o rc err -i "${OUT}/big_${label}.fa" -o "${OUT}/x.fa" --mode accurate -n 1)
  if(rc EQUAL 0 OR NOT err MATCHES "not feasible above 5000")
    fail("${label}: accurate on 5001 sequences was not refused")
  else()
    message(STATUS "ok   ${label}: accurate on 5001 sequences refused")
  endif()
endforeach()

if(n_fail GREATER 0)
  message(FATAL_ERROR "${n_fail} seq_type CLI check(s) failed")
endif()
