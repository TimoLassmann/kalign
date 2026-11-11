# Command-line option rules of kalign.
#
#  1. Fast-only options are rejected with --mode accurate (and recall).
#  2. Accurate-only options are rejected with --mode fast, default and no mode.
#  3. Numbers are parsed strictly; bad values name the option and the value.
#  4. No input: bare "kalign" prints help (exit 0); options without input fail.
#  5. Output format names are matched exactly.
#  6. --add/--existing cannot be combined with --mode or mode-specific options.
#  7. Valid commands still work in fast and accurate mode.
#
# Invoked by ctest as:
#   cmake -DKALIGN=<kalign binary> -DDATA=<tests/data> -DOUT=<scratch dir> -P cli_options.cmake

foreach(v KALIGN DATA OUT)
  if(NOT DEFINED ${v})
    message(FATAL_ERROR "cli_options.cmake: -D${v}=... is required")
  endif()
endforeach()
file(MAKE_DIRECTORY "${OUT}")
set(n_fail 0)

macro(fail msg)
  message(SEND_ERROR "FAIL ${msg}")
  math(EXPR n_fail "${n_fail} + 1")
endmacro()

function(run_kalign rc_var msg_var)
  execute_process(COMMAND "${KALIGN}" ${ARGN}
                  RESULT_VARIABLE rc OUTPUT_VARIABLE so ERROR_VARIABLE se)
  set(${rc_var} "${rc}" PARENT_SCOPE)
  set(${msg_var} "${so}${se}" PARENT_SCOPE)
endfunction()

# expect_error(<regex> <kalign arguments...>): exit status 1 and message
function(expect_error regex)
  run_kalign(rc msg ${ARGN})
  string(REPLACE ";" " " cmd "${ARGN}")
  if(NOT rc EQUAL 1)
    message(SEND_ERROR "FAIL [${cmd}]: exit status ${rc}, expected 1\n${msg}")
  elseif(NOT msg MATCHES "${regex}")
    message(SEND_ERROR "FAIL [${cmd}]: message '${regex}' not found\n${msg}")
  else()
    message(STATUS "ok   [${cmd}] rejected")
  endif()
endfunction()

function(expect_ok)
  run_kalign(rc msg ${ARGN})
  string(REPLACE ";" " " cmd "${ARGN}")
  if(NOT rc EQUAL 0)
    message(SEND_ERROR "FAIL [${cmd}]: exit status ${rc}\n${msg}")
  else()
    message(STATUS "ok   [${cmd}]")
  endif()
endfunction()

set(prot "${DATA}/BB11001.tfa")
set(dna "${DATA}/small.fa")
set(o "${OUT}/o.fa")
set(base -i "${prot}" -o "${o}" -n 1 -q)

# 1. fast-only options with accurate / recall
set(fast_only "--gpo|5" "--gpe|1" "--tgpe|1" "--vsm-mode|1" "--vsm-kappa|0.5" "--vsm-dzero|0.5" "--guide-dist|0")
foreach(mode accurate recall)
  foreach(opt IN LISTS fast_only)
    string(REPLACE "|" ";" optl "${opt}")
    list(GET optl 0 name)
    expect_error("${name} cannot be used with --mode accurate: the accurate ensemble uses its fitted settings. Use --mode fast to set alignment parameters."
                 ${base} --mode ${mode} ${optl})
  endforeach()
endforeach()

# 2. accurate-only options with fast / default / no mode
set(acc_only "--fill-in|1" "--fill-anchor-frac|0.5" "--fill-support|2")
foreach(opt IN LISTS acc_only)
  string(REPLACE "|" ";" optl "${opt}")
  list(GET optl 0 name)
  foreach(modeargs "--mode;fast" "--mode;default" "")
    expect_error("${name} cannot be used with --mode fast: it only applies to the accurate ensemble. Use --mode accurate."
                 ${base} ${modeargs} ${optl})
  endforeach()
endforeach()

# 3. strict numbers (valid mode, so only the parse can fail)
set(bad_numbers
  "--gpo|abc" "--gpo|-5" "--gpo|1x" "--gpo|nan" "--gpo|inf"
  "--gpe|abc" "--gpe|-1" "--tgpe|-0.5"
  "--vsm-kappa|-1" "--vsm-kappa|x" "--vsm-dzero|-1" "--vsm-dzero|1.5z"
  "--guide-dist|1abc" "--guide-dist|2" "--guide-dist|-1"
  "-n|abc" "-n|0" "-n|-2" "-n|1.5" "-n|2x"
  "--confidence-threshold|abc" "--confidence-threshold|1.5" "--confidence-threshold|-0.1"
  "--confidence-style|remov" "--confidence-style|Lowercase")
foreach(bn IN LISTS bad_numbers)
  string(REPLACE "|" ";" f "${bn}")
  list(GET f 0 name)
  list(GET f 1 value)
  expect_error("${name}.*'${value}'" -i "${prot}" -o "${o}" -q ${name} "${value}")
endforeach()
expect_error("--fill-support must be an integer >= 1, got '0'" ${base} --mode accurate --fill-support 0)
expect_error("--fill-support must be an integer >= 1, got 'x'" ${base} --mode accurate --fill-support x)
expect_error("Unknown mode: 'bogus'" ${base} --mode bogus)

# 4. no input
execute_process(COMMAND "${KALIGN}" RESULT_VARIABLE rc OUTPUT_VARIABLE so ERROR_VARIABLE se)
if(NOT rc EQUAL 0 OR NOT so MATCHES "Usage")
  fail("kalign without arguments: exit ${rc}, expected help and exit 0")
else()
  message(STATUS "ok   no arguments prints help")
endif()
expect_error("No input files" -n 1 -q)
expect_error("No input files" --mode accurate -q -o "${o}")

# 5. format names
foreach(fmt bleah safari fastafoo msfx clux faa)
  expect_error("Format ${fmt} not recognized" ${base} -f "${fmt}")
endforeach()
# Names are matched case-insensitively.
foreach(fmt fasta fa afa msf clu clustal FASTA Fa MSF Clu CLUSTAL)
  expect_ok(${base} -f ${fmt})
endforeach()

# 6. --add / --existing conflicts
run_kalign(rc msg -i "${prot}" -o "${OUT}/base.fa" -n 1 -q)
file(WRITE "${OUT}/new.fa" ">added_one\nMKTAYIAKQRQISFVKSHFSRQ\n")
set(add_base --add "${OUT}/new.fa" --existing "${OUT}/base.fa" -o "${OUT}/added.fa" -n 1 -q)
set(add_conflicts "--mode|fast" "--mode|accurate" "--gpo|5" "--gpe|1" "--tgpe|1" "--vsm-mode|1"
    "--vsm-kappa|0.5" "--vsm-dzero|0.5" "--guide-dist|0" "--fill-in|1" "--fill-anchor-frac|0.5" "--fill-support|2")
foreach(opt IN LISTS add_conflicts)
  string(REPLACE "|" ";" optl "${opt}")
  list(GET optl 0 name)
  expect_error("${name} cannot be used with --add: adding sequences uses its own settings."
               ${add_base} ${optl})
endforeach()
expect_ok(${add_base})
expect_error("Both --add and --existing must be specified together" --add "${OUT}/new.fa" -o "${o}" -q)

# 7. valid commands, fast and accurate
foreach(dataset "${prot}" "${dna}")
  expect_ok(-i "${dataset}" -o "${o}" -n 1 -q)
  expect_ok(-i "${dataset}" -o "${o}" -n 1 -q --mode fast --gpo 8 --gpe 1.5 --tgpe 0.5)
  expect_ok(-i "${dataset}" -o "${o}" -n 1 -q --mode default --vsm-mode 1 --vsm-kappa 0.5 --vsm-dzero 0.5 --guide-dist 0)
  expect_ok(-i "${dataset}" -o "${o}" -n 1 -q --mode accurate)
  expect_ok(-i "${dataset}" -o "${o}" -n 1 -q --mode recall --fill-in 1 --fill-anchor-frac 0.5 --fill-support 2)
  expect_ok(-i "${dataset}" -o "${o}" -n 1 -q --mode accurate --confidence-threshold 0.5 --confidence-style remove)
endforeach()
# confidence options with fast: allowed, warns
run_kalign(rc msg -i "${prot}" -o "${o}" -n 1 --mode fast --confidence-threshold 0.5)
if(NOT rc EQUAL 0 OR NOT msg MATCHES "confidence needs --mode accurate")
  fail("--confidence-threshold with fast: exit ${rc}, expected success with a warning\n${msg}")
else()
  message(STATUS "ok   confidence options with fast warn")
endif()

# mixed-case mode names follow the same rules
expect_error("--gpo cannot be used with --mode accurate" ${base} --mode ACCURATE --gpo 5)
expect_error("--fill-in cannot be used with --mode fast" ${base} --mode Fast --fill-in 1)
# confidence-style alone with fast: allowed
expect_ok(${base} --mode fast --confidence-style remove)
# protein input from stdin
execute_process(COMMAND "${KALIGN}" -i - -o "${o}" -n 1 -q
                INPUT_FILE "${prot}" RESULT_VARIABLE rc OUTPUT_VARIABLE so ERROR_VARIABLE se)
if(NOT rc EQUAL 0)
  fail("stdin input: exit ${rc}\n${so}${se}")
else()
  message(STATUS "ok   stdin input")
endif()

if(n_fail GREATER 0)
  message(FATAL_ERROR "${n_fail} CLI option check(s) failed")
endif()
