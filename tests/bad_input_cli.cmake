# Malformed input must end with a clear error and exit status 1, never a crash;
# an empty sequence among others is dropped (as in Kalign 3.6) and the rest aligned.
# Arguments: KALIGN (executable), WORKDIR, MODE
file(MAKE_DIRECTORY "${WORKDIR}")
set(d "${WORKDIR}")
file(WRITE "${d}/empty.fa" "")
file(WRITE "${d}/header_only.fa" ">a\n")
file(WRITE "${d}/one.fa" ">a\nACGTACGT\n")
file(WRITE "${d}/no_header.fa" "ACGTACGT\nACGTTCGT\n")
file(WRITE "${d}/empty_seq.fa" ">a\n\n>b\nACGTACGT\n>c\nACGTTCGT\n")
file(WRITE "${d}/ok.fa" ">a\nACGTACGT\n>b\nACGTTCGT\n")

# name | input | expected rc | expected message (regex) | output file (optional)
set(cases
  "empty|empty.fa|1|No sequences were found|"
  "header only|header_only.fa|1|Only 1 sequence was found|"
  "one sequence|one.fa|1|Only 1 sequence was found|"
  "no header|no_header.fa|1|not FASTA, MSF or Clustal format|"
  "missing file|missing.fa|1|does not exist|"
  "unwritable output|ok.fa|1|Could not open output file|${d}/no/such/dir/out.fa"
  "empty sequence|empty_seq.fa|0|Removing 1 sequence with a length of 0|")
foreach(c IN LISTS cases)
  string(REPLACE "|" ";" f "${c}")
  list(GET f 0 name)
  list(GET f 1 input)
  list(GET f 2 want_rc)
  list(GET f 3 want_msg)
  list(LENGTH f nf)
  set(out "${d}/out_${MODE}.fa")
  if(nf GREATER 4)
    list(GET f 4 o)
    if(o)
      set(out "${o}")
    endif()
  endif()
  execute_process(
    COMMAND "${KALIGN}" --mode ${MODE} -n 2 -i "${d}/${input}" -o "${out}"
    RESULT_VARIABLE rc OUTPUT_VARIABLE so ERROR_VARIABLE se)
  if(NOT rc STREQUAL "${want_rc}")
    message(FATAL_ERROR "${name} (${MODE}): exit '${rc}', expected ${want_rc}\n${so}${se}")
  endif()
  if(NOT "${so}${se}" MATCHES "${want_msg}")
    message(FATAL_ERROR "${name} (${MODE}): message '${want_msg}' not found\n${so}${se}")
  endif()
  if(want_rc EQUAL 0)
    file(STRINGS "${out}" hdr REGEX "^>")
    list(LENGTH hdr nh)
    if(NOT nh EQUAL 2)
      message(FATAL_ERROR "${name} (${MODE}): ${nh} sequences written, expected 2")
    endif()
  endif()
endforeach()
