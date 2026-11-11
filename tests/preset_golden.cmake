# Preset golden test (Kalign 4).
#
# Runs `kalign --mode <m> -n 1` for both modes (and the "default" alias) on
# small protein and DNA inputs in tests/data and compares the SHA-256 of each
# output with golden hashes. The golden outputs were generated through the
# Python path the re-fit evaluated (kalign._core.ensemble_custom_file_to_file
# with the parameters of presets/results/v2/presets_v2.json; protein fast:
# presets/results/fastfit/final_pick.json with the member gap model; protein
# accurate: presets/results/ensconf/final_pick.json), so this test
# checks that lib/src/mode_presets.h reproduces the evaluated presets byte
# for byte.
#
# To regenerate (only when the presets are deliberately changed):
#     uv sync && uv run python tests/preset_golden_regen.py
# and paste the printed lines into GOLDEN below.
#
# Invoked by ctest as:
#   cmake -DKALIGN=<kalign binary> -DDATA=<tests/data> -DOUT=<scratch dir> -P preset_golden.cmake

set(GOLDEN
  "fast:BB11001.tfa:a59470180769f9dfa6b786a29cd510141c9bd4f1dd3b46f15a63cab83ad429ac"
  "fast:BB12006.tfa:99bb7a20f06c867c310ac44ad9bcc968671b7b97aea3250322de1642ad172a66"
  "fast:BB30014.tfa:7717e4968416950f280e7abb51984002796e98a51d836b994a02078df0c13a89"
  "accurate:BB11001.tfa:cdde988057d32d372b214f46eb049fb8ac926153d803e3ec2bbeb5826cbe43ce"
  "accurate:BB12006.tfa:305a04f487b3f368d207dda1164d623465f3da0593f3485e665df64e4cbbcf9b"
  "accurate:BB30014.tfa:33e576b8d7ef3c026f2b7e8ee6fd6a716c7e1b13b229f7007635b64311d7cb29"
  "default:BB11001.tfa:a59470180769f9dfa6b786a29cd510141c9bd4f1dd3b46f15a63cab83ad429ac"
  "default:BB12006.tfa:99bb7a20f06c867c310ac44ad9bcc968671b7b97aea3250322de1642ad172a66"
  "default:BB30014.tfa:7717e4968416950f280e7abb51984002796e98a51d836b994a02078df0c13a89"
  "fast:small.fa:d3d7aced75726181371f7b45097f2f60bb95684c3d604070e6d479d3c65e1a53"
  "accurate:small.fa:7f491b6b8e7ee190bee5ea7b7fcafb8bb58dcd524aa071c7aaa4fb36bb5c1db8"
)

foreach(v KALIGN DATA OUT)
  if(NOT DEFINED ${v})
    message(FATAL_ERROR "preset_golden.cmake: -D${v}=... is required")
  endif()
endforeach()

file(MAKE_DIRECTORY "${OUT}")
set(n_fail 0)
foreach(entry IN LISTS GOLDEN)
  string(REPLACE ":" ";" parts "${entry}")
  list(GET parts 0 mode)
  list(GET parts 1 input)
  list(GET parts 2 expected)
  set(out_file "${OUT}/${mode}_${input}.fasta")
  file(REMOVE "${out_file}")
  execute_process(
    COMMAND "${KALIGN}" --mode ${mode} -n 1 -i "${DATA}/${input}" -o "${out_file}"
    RESULT_VARIABLE rc
    OUTPUT_QUIET ERROR_QUIET)
  if(NOT rc EQUAL 0 OR NOT EXISTS "${out_file}")
    message(SEND_ERROR "FAIL ${mode} ${input}: kalign exited with '${rc}'")
    math(EXPR n_fail "${n_fail} + 1")
    continue()
  endif()
  file(SHA256 "${out_file}" got)
  if(got STREQUAL expected)
    message(STATUS "ok   ${mode} ${input}")
  else()
    message(SEND_ERROR "FAIL ${mode} ${input}: sha256 ${got}, expected ${expected}")
    math(EXPR n_fail "${n_fail} + 1")
  endif()
endforeach()

if(n_fail GREATER 0)
  message(FATAL_ERROR "${n_fail} preset golden check(s) failed")
endif()
