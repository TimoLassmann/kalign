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

# Protein presets default/accurate/recall must reproduce the confirmed 3.7.0
# arms byte for byte (golden SHA-256; see preset_golden.cmake to regenerate).
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
