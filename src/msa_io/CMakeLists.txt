add_library(module_msa STATIC
  rwalign.c
  #align_io.c
  misc.c
  )

# if(OpenMP_C_FOUND)
#     target_link_libraries(module_msa PUBLIC OpenMP::OpenMP_C)
# endif()


target_include_directories(module_msa
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
  )
target_link_libraries(module_msa PRIVATE
  tldevel
  module_alphabet
  module_param)

include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_executable(misc_utest
  misc.c  
  )
target_link_libraries(misc_utest PRIVATE tldevel module_alphabet)
set_target_properties(misc_utest PROPERTIES COMPILE_FLAGS "-DUTEST_MISC")
add_test(
  NAME misc_utest
  COMMAND misc_utest
  )

add_executable(rwaln_utest
  rwalign.c
  )
target_link_libraries(rwaln_utest PRIVATE
  tldevel
  module_alphabet
  module_param
  )
set_target_properties(rwaln_utest PROPERTIES COMPILE_FLAGS "-DRWALIGN_TEST")
add_test(
  NAME rwaln_utest_a2m.good.1
  COMMAND rwaln_utest ${CMAKE_CURRENT_SOURCE_DIR}/test_data/a2m.good.1
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )

add_test(
  NAME rwaln_utest_a2m.good.2
  COMMAND rwaln_utest ${CMAKE_CURRENT_SOURCE_DIR}/test_data/a2m.good.2
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
  
add_test(
  NAME rwaln_utest_clustal.good.1
  COMMAND rwaln_utest ${CMAKE_CURRENT_SOURCE_DIR}/test_data/clustal.good.1
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )

add_test(
  NAME rwaln_utest_clustal.good.2
  COMMAND rwaln_utest ${CMAKE_CURRENT_SOURCE_DIR}/test_data/clustal.good.2
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )

   
add_test(
  NAME rwaln_utest_afa.good.1
  COMMAND rwaln_utest ${CMAKE_CURRENT_SOURCE_DIR}/test_data/afa.good.1
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )

add_test(
  NAME rwaln_utest_afa.good.2
  COMMAND rwaln_utest ${CMAKE_CURRENT_SOURCE_DIR}/test_data/afa.good.2
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )


add_test(
  NAME rwaln_utest_afa.good.3
  COMMAND rwaln_utest ${CMAKE_CURRENT_SOURCE_DIR}/test_data/afa.good.3
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )



# add_executable(mcc_history_itest
#   mcc_history_plot_test.c
#   )

# target_link_libraries(mcc_history_itest PRIVATE tld-dev module_mcc)

# add_executable(queue_itest
#   queue.c
#   queue_itest.c
#   )
# target_link_libraries(queue_itest PRIVATE tld-dev)
# add_test(
#     NAME queue_itest
#     COMMAND queue_itest
# )
