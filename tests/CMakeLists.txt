cmake_minimum_required(VERSION 3.18)

# Build library in sub-directories
if (NOT TARGET kalign::kalign)
  add_subdirectory(${PROJECT_SOURCE_DIR}/lib build EXCLUDE_FROM_ALL)
endif()

# link to installed library 
if (NOT TARGET kalign::kalign)
  find_package(kalign)
endif()

add_executable(largebench
  large_benchmark.c
  )
target_link_libraries(largebench tldevel   ${PROJECT_NAME}_static)

add_executable(kalign_io_test
  kalign_io_test.c
  )

target_link_libraries(kalign_io_test kalign::kalign)


add_executable(dssim
  dssim_test.c
  dssim.c
  )
target_link_libraries(dssim tldevel kalign::kalign)




add_executable(kalign_lib_test
  kalign_lib_test.c
  )
target_link_libraries(kalign_lib_test kalign::kalign)

add_executable(kalign_cmp_test
  kalign_cmp_test.c
  )
target_link_libraries(kalign_cmp_test kalign::kalign)

add_executable(kaligncpp
   kalign_lib_testCXX.cpp)
target_link_libraries(kaligncpp kalign::kalign)

add_test(
  NAME C++Lib
  COMMAND kaligncpp
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )


#BB11001.tfa     
add_test(
  NAME kalign_itest_BB11001.tfa
  COMMAND kalign-bin ${CMAKE_CURRENT_SOURCE_DIR}/data/BB11001.tfa
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )

#BB11001_EOF.msf
add_test(
  NAME kalign_itest_BB11001_EOF.msf
  COMMAND kalign-bin ${CMAKE_CURRENT_SOURCE_DIR}/data/BB11001_EOF.msf
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )

#BB12006.msf
add_test(
  NAME kalign_itest_BB12006.msf
  COMMAND kalign-bin ${CMAKE_CURRENT_SOURCE_DIR}/data/BB12006.msf
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )

#BB12006.tfa
add_test(
  NAME kalign_itest_BB12006.tfa
  COMMAND kalign-bin ${CMAKE_CURRENT_SOURCE_DIR}/data/BB12006.tfa
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )

#BB30014.msf     
add_test(
  NAME kalign_itest_BB30014.msf
  COMMAND kalign-bin ${CMAKE_CURRENT_SOURCE_DIR}/data/BB30014.msf
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )

#BB30014.tfa
add_test(
  NAME kalign_itest_BB30014.tfa
  COMMAND kalign-bin ${CMAKE_CURRENT_SOURCE_DIR}/data/BB30014.tfa
  WORKING_DIRECTORY 
${CMAKE_CURRENT_SOURCE_DIR}
  )


add_executable(kalign_ensemble_test
  kalign_ensemble_test.c
  )
target_link_libraries(kalign_ensemble_test tldevel ${PROJECT_NAME}_static)

add_test(
  NAME kalign_ensemble_test
  COMMAND kalign_ensemble_test ${CMAKE_CURRENT_SOURCE_DIR}/data/BB11001.tfa
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )

# BB40023 yields consensus alignments wider than 4096 columns, which
# previously overflowed the cycle check and scrambled residue order.
add_test(
  NAME kalign_ensemble_test_long
  COMMAND kalign_ensemble_test ${CMAKE_CURRENT_SOURCE_DIR}/data/BB40023.tfa
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )

# The Kalign 4.0 fast and accurate presets (protein and nucleotide) and the
# "default" alias must reproduce the golden SHA-256 values byte for byte (see
# preset_golden.cmake to regenerate). The "recall" alias is not run here.
add_test(
  NAME kalign_preset_golden
  COMMAND ${CMAKE_COMMAND}
          -DKALIGN=$<TARGET_FILE:kalign-bin>
          -DDATA=${CMAKE_CURRENT_SOURCE_DIR}/data
          -DOUT=${CMAKE_CURRENT_BINARY_DIR}/preset_golden
          -P ${CMAKE_CURRENT_SOURCE_DIR}/preset_golden.cmake
  )

# --type (override of the detected sequence type) and the accurate size limit.
add_test(
  NAME kalign_seq_type_cli
  COMMAND ${CMAKE_COMMAND}
          -DKALIGN=$<TARGET_FILE:kalign-bin>
          -DDATA=${CMAKE_CURRENT_SOURCE_DIR}/data
          -DOUT=${CMAKE_CURRENT_BINARY_DIR}/seq_type_cli
          -P ${CMAKE_CURRENT_SOURCE_DIR}/seq_type_cli.cmake
  )

# Option rules, strict number parsing, format names and --add conflicts.
add_test(
  NAME kalign_cli_options
  COMMAND ${CMAKE_COMMAND}
          -DKALIGN=$<TARGET_FILE:kalign-bin>
          -DDATA=${CMAKE_CURRENT_SOURCE_DIR}/data
          -DOUT=${CMAKE_CURRENT_BINARY_DIR}/cli_options
          -P ${CMAKE_CURRENT_SOURCE_DIR}/cli_options.cmake
  )

add_executable(kalign_api_test
  kalign_api_test.c
  )
target_link_libraries(kalign_api_test tldevel ${PROJECT_NAME}_static)

add_test(
  NAME kalign_api_test
  COMMAND kalign_api_test ${CMAKE_CURRENT_SOURCE_DIR}/data/BB11001.tfa
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )

add_test(
  NAME kalign_api_test_large
  COMMAND kalign_api_test ${CMAKE_CURRENT_SOURCE_DIR}/data/BB30014.tfa
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )

add_test(
  NAME DSSIM_takes_a_min
  COMMAND dssim
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )

add_executable(poar_key_test
  poar_key_test.c
  )
target_link_libraries(poar_key_test tldevel ${PROJECT_NAME}_static)

add_test(
  NAME poar_key_test
  COMMAND poar_key_test ${CMAKE_CURRENT_BINARY_DIR}/poar_key_test.poar
  )

add_executable(threadpool_test
  ${PROJECT_SOURCE_DIR}/lib/src/threadpool/test_threadpool.c
  )
target_include_directories(threadpool_test PRIVATE ${PROJECT_SOURCE_DIR}/lib/src/threadpool)
target_link_libraries(threadpool_test ${PROJECT_NAME}_static)

add_test(
  NAME threadpool_tests
  COMMAND threadpool_test
  )

# Malformed input: a clear error and exit status 1, never a crash.
foreach(mode fast accurate)
  add_test(
    NAME kalign_bad_input_${mode}
    COMMAND ${CMAKE_COMMAND} -DKALIGN=$<TARGET_FILE:kalign-bin> -DWORKDIR=${CMAKE_CURRENT_BINARY_DIR}/bad_input_${mode} -DMODE=${mode} -P ${CMAKE_CURRENT_SOURCE_DIR}/bad_input_cli.cmake
  )
endforeach()

# More threads than cores: capped at the number of cores (before, idle workers spun on
# the oversubscribed cores and a tiny input took seconds).
add_test(
  NAME kalign_threads_capped
  COMMAND kalign-bin -i ${CMAKE_CURRENT_SOURCE_DIR}/data/BB11001.tfa -o ${CMAKE_CURRENT_BINARY_DIR}/threads_capped.fa -n 100000
  )
set_tests_properties(kalign_threads_capped PROPERTIES PASS_REGULAR_EXPRESSION "the number of cores" TIMEOUT 30)

# --type dna on protein: one summary line for the replaced letters, not a flood.
foreach(mode fast accurate)
  add_test(
    NAME kalign_foreign_chars_${mode}
    COMMAND ${CMAKE_COMMAND} -DKALIGN=$<TARGET_FILE:kalign-bin> -DINPUT=${CMAKE_CURRENT_SOURCE_DIR}/data/BB11001.tfa -DWORKDIR=${CMAKE_CURRENT_BINARY_DIR}/foreign_chars -DMODE=${mode} -P ${CMAKE_CURRENT_SOURCE_DIR}/foreign_chars_test.cmake
    )
endforeach()

# Deep guide trees must not depend on the stack size of the worker threads.
add_test(
  NAME kalign_deep_tree_accurate_n8
  COMMAND ${CMAKE_COMMAND} -DKALIGN=$<TARGET_FILE:kalign-bin> -DWORKDIR=${CMAKE_CURRENT_BINARY_DIR}/deep_tree -DNSEQ=2000 -DMODE=accurate -DTHREADS=8 -P ${CMAKE_CURRENT_SOURCE_DIR}/deep_tree_test.cmake
  )
add_test(
  NAME kalign_deep_tree_fast_n8
  COMMAND ${CMAKE_COMMAND} -DKALIGN=$<TARGET_FILE:kalign-bin> -DWORKDIR=${CMAKE_CURRENT_BINARY_DIR}/deep_tree -DNSEQ=20000 -DMODE=fast -DTHREADS=8 -P ${CMAKE_CURRENT_SOURCE_DIR}/deep_tree_test.cmake
  )

# The library (kalign, kalign_align_mode) gives the same alignments as the command line.
add_executable(kalign_unity_test
  kalign_unity_test.c
  )
target_link_libraries(kalign_unity_test tldevel ${PROJECT_NAME}_static)

add_test(
  NAME kalign_unity_test
  COMMAND kalign_unity_test $<TARGET_FILE:kalign-bin>
          ${CMAKE_CURRENT_SOURCE_DIR}/data/BB11001.tfa
          ${CMAKE_CURRENT_SOURCE_DIR}/data/small.fa
          ${CMAKE_CURRENT_BINARY_DIR}/unity_test
  )

# MSF header text (AA/NA, Type P/N) and exact output format names.
add_test(
  NAME kalign_msf_header
  COMMAND ${CMAKE_COMMAND} -DKALIGN=$<TARGET_FILE:kalign-bin> -DDATA=${CMAKE_CURRENT_SOURCE_DIR}/data -DOUT=${CMAKE_CURRENT_BINARY_DIR}/msf_header -P ${CMAKE_CURRENT_SOURCE_DIR}/msf_header.cmake
  )

# Run-time vector switch: the aligned allocator, the plain-kernel switch and what
# kalign_simd_name() reports.
add_executable(kalign_simd_test
  kalign_simd_test.c
  )
target_link_libraries(kalign_simd_test tldevel ${PROJECT_NAME}_static)

add_test(
  NAME kalign_simd_test
  COMMAND kalign_simd_test
  )
add_test(
  NAME kalign_simd_test_scalar
  COMMAND kalign_simd_test
  )
set_tests_properties(kalign_simd_test_scalar PROPERTIES ENVIRONMENT "KALIGN_SIMD=scalar")

# The plain kernels (KALIGN_SIMD=scalar) must reproduce the same golden SHA-256
# values as the vector kernels (kalign_preset_golden above).
add_test(
  NAME kalign_preset_golden_scalar
  COMMAND ${CMAKE_COMMAND}
          -DDATA=${CMAKE_CURRENT_SOURCE_DIR}/data
          -DKALIGN=$<TARGET_FILE:kalign-bin>
          -DOUT=${CMAKE_CURRENT_BINARY_DIR}/preset_golden_scalar
          -DGOLDEN_FILE=${CMAKE_CURRENT_SOURCE_DIR}/preset_golden.cmake
          -P ${CMAKE_CURRENT_SOURCE_DIR}/simd_compare.cmake
  )

# Vector and plain kernels must give identical alignments, with one and with
# several threads. Protein inputs above 1,163 sequences and DNA inputs above
# 1,886 (tree_leaf in lib/src/mode_presets.h) go through the k-means guide
# tree; below that the whole tree is single linkage (sl_tree.c). So the
# 1,300-sequence protein cases reach k-means and its Euclidean distance, the
# 1,300-sequence DNA cases the single-linkage tree on the Myers distance, and
# the 2,000-sequence DNA case k-means. These tests do not notice a change in
# the summation order of the Euclidean distance (the alignments stay the
# same); edist_order_test below is the test for that.
foreach(mode fast accurate)
  foreach(alphabet protein dna)
    add_test(
      NAME kalign_simd_vs_scalar_${mode}_${alphabet}_n1
      COMMAND ${CMAKE_COMMAND}
              -DKALIGN=$<TARGET_FILE:kalign-bin>
              -DOUT=${CMAKE_CURRENT_BINARY_DIR}/simd_vs_scalar_${mode}_${alphabet}
              -DMODE=${mode} -DTHREADS=1 -DNSEQ=1300 -DALPHABET=${alphabet}
              -P ${CMAKE_CURRENT_SOURCE_DIR}/simd_compare.cmake
      )
    set_tests_properties(kalign_simd_vs_scalar_${mode}_${alphabet}_n1 PROPERTIES TIMEOUT 600)
  endforeach()
endforeach()
add_test(
  NAME kalign_simd_vs_scalar_fast_protein_n4
  COMMAND ${CMAKE_COMMAND}
          -DKALIGN=$<TARGET_FILE:kalign-bin>
          -DOUT=${CMAKE_CURRENT_BINARY_DIR}/simd_vs_scalar_fast_protein_n4
          -DMODE=fast -DTHREADS=4 -DNSEQ=1300 -DALPHABET=protein
          -P ${CMAKE_CURRENT_SOURCE_DIR}/simd_compare.cmake
  )
set_tests_properties(kalign_simd_vs_scalar_fast_protein_n4 PROPERTIES TIMEOUT 600)
add_test(
  NAME kalign_simd_vs_scalar_fast_dna2000_n1
  COMMAND ${CMAKE_COMMAND}
          -DKALIGN=$<TARGET_FILE:kalign-bin>
          -DOUT=${CMAKE_CURRENT_BINARY_DIR}/simd_vs_scalar_fast_dna2000
          -DMODE=fast -DTHREADS=1 -DNSEQ=2000 -DALPHABET=dna
          -P ${CMAKE_CURRENT_SOURCE_DIR}/simd_compare.cmake
  )
set_tests_properties(kalign_simd_vs_scalar_fast_dna2000_n1 PROPERTIES TIMEOUT 600)

# Guide-tree ties (tests/data/sl_ties): DNA inputs in which two normalised
# Myers distances of the single-linkage tree are equal as fractions (5/41 and
# 15/123) but differ by one unit in the last place as floating-point numbers.
# The plain, NEON and AVX2 kernels stop at different points once a distance
# cannot beat the current link, so the cut-off must be exact for the
# comparison (sl_bpm_plan_pair in lib/src/sl_tree.c); before that was fixed
# these inputs gave different alignments with and without the vector kernels.
# Each output is also pinned by its SHA-256 (fast, accurate). On one platform
# the plain and vector kernels must give the same alignment. Between platforms
# (Arm and x86-64) an exact tie can still be decided differently, because the
# Arm build lets the compiler fuse multiply-adds in the alignment code; these
# six inputs happen to give the same alignment on both, but if a pinned value
# fails on one platform only, look there first.
set(SL_TIES
  "t1:e17df29502290d68b1bf6f1a247351d81b40606dff1883277177a314b2b66640:0e02b19f5eeb4fa4c2ab9c25067ecdb23d859e523bdf32ba93af99464c9c412f"
  "t2:d98873baeee74b4cc4af4ae6c8af277700390bf0d7fa1db32b7e077ecd85dbda:79b7780807724787b63ab6ed4f15bc7a381454e369381e446489f7d42c2ca31d"
  "d3_1_11:e7507f1e2d394e81c1af20d44454b9d695b805b40b59f1fd1df6dc125f62e931:26b1873421d206c4c5b45ce729a4c90829c4bb00c7504ac9ba9a77fdd96fe4fa"
  "d5_1_9:5f94c3444d7b0bd842e0ac6311882fe3d08570c084007d49bbc4d410cfd8e6b4:e70197183fdb22b4cb2307c9752842ffe2ab9d4d0c79e5b2aaaa2bf7fcf63d2b"
  "d5_1_14:ef0a4e01739a62c7d31d8f49e9a17801f4faa7a52962d8321a5dff558e307325:25125dfa008f372e8972055d95dd64509c17a8021c3e9c3c88164e2660714e53"
  "d7_1_10:aff386b851f9c5faebfa4a45e648b7ee0fb4fb5425f859c93fe0434385cd9789:a928cf797eab74471052cd5f82fd1c2eaeefb588120cbf3cafb6a03a97e1da47"
)
foreach(entry IN LISTS SL_TIES)
  string(REPLACE ":" ";" parts "${entry}")
  list(GET parts 0 tie)
  list(GET parts 1 sha_fast)
  list(GET parts 2 sha_accurate)
  foreach(mode fast accurate)
    add_test(
      NAME kalign_sl_ties_${tie}_${mode}
      COMMAND ${CMAKE_COMMAND}
              -DKALIGN=$<TARGET_FILE:kalign-bin>
              -DOUT=${CMAKE_CURRENT_BINARY_DIR}/sl_ties
              -DMODE=${mode} -DTHREADS=1
              -DINPUT=${CMAKE_CURRENT_SOURCE_DIR}/data/sl_ties/${tie}.fa
              -DEXPECT=${sha_${mode}}
              -P ${CMAKE_CURRENT_SOURCE_DIR}/simd_compare.cmake
      )
  endforeach()
endforeach()

# The plain Euclidean distance must round exactly like the AVX2 one (bitwise),
# which is also compared directly when the CPU has AVX2. This is the only test
# that notices a change in the summation order. Its reference sum must not be
# fused into multiply-adds: GCC ignores the FP_CONTRACT pragma in the file and
# fuses by default on Arm, POWER, s390x and RISC-V, so the flag is needed.
add_executable(edist_order_test
  edist_order_test.c
  )
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
  set_source_files_properties(edist_order_test.c PROPERTIES COMPILE_OPTIONS "-ffp-contract=off")
endif()
target_link_libraries(edist_order_test tldevel ${PROJECT_NAME}_static)

add_test(
  NAME edist_order_test
  COMMAND edist_order_test
  )
add_test(
  NAME edist_order_test_scalar
  COMMAND edist_order_test
  )
set_tests_properties(edist_order_test_scalar PROPERTIES ENVIRONMENT "KALIGN_SIMD=scalar")
