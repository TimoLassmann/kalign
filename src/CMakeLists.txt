
add_subdirectory(tldevel)
add_subdirectory(interface)
add_subdirectory(msa_io)
add_subdirectory(alphabet)
add_subdirectory(tasks)
add_subdirectory(tree)
add_subdirectory(alignment)
add_subdirectory(tools)


add_executable(kalign
  run_kalign.c
  )
if(OpenMP_C_FOUND)
    target_link_libraries(kalign PUBLIC OpenMP::OpenMP_C)
endif()

target_link_libraries(kalign PRIVATE
  tldevel
  module_param
  module_alphabet
  module_tree
  module_aln)



#BB11001.msf
add_test(
  NAME kalign_itest_BB11001.msf
  COMMAND kalign ${CMAKE_CURRENT_SOURCE_DIR}/test_data/BB11001.msf
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )

#BB11001.tfa     
add_test(
  NAME kalign_itest_BB11001.tfa
  COMMAND kalign ${CMAKE_CURRENT_SOURCE_DIR}/test_data/BB11001.tfa
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )

#BB11001_EOF.msf
add_test(
  NAME kalign_itest_BB11001_EOF.msf
  COMMAND kalign ${CMAKE_CURRENT_SOURCE_DIR}/test_data/BB11001_EOF.msf
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )

#BB12006.msf
add_test(
  NAME kalign_itest_BB12006.msf
  COMMAND kalign ${CMAKE_CURRENT_SOURCE_DIR}/test_data/BB12006.msf
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )

#BB12006.tfa
add_test(
  NAME kalign_itest_BB12006.tfa
  COMMAND kalign ${CMAKE_CURRENT_SOURCE_DIR}/test_data/BB12006.tfa
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )

#BB30014.msf     
add_test(
  NAME kalign_itest_BB30014.msf
  COMMAND kalign ${CMAKE_CURRENT_SOURCE_DIR}/test_data/BB12006.tfa
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )

#BB30014.tfa
add_test(
  NAME kalign_itest_BB30014.tfa
  COMMAND kalign ${CMAKE_CURRENT_SOURCE_DIR}/test_data/BB30014.tfa
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )





