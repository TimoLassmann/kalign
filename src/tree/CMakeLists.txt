add_library(module_tree STATIC
  bisectingKmeans.c
  sequence_distance.c
  euclidean_dist.c
  pick_anchor.c
  bpm.c
  )

if(OpenMP_C_FOUND)
    target_link_libraries(module_tree PUBLIC OpenMP::OpenMP_C)
endif()


target_include_directories(module_tree
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
  )
target_link_libraries(module_tree PRIVATE
  tldevel
  module_msa
  module_param
  module_alphabet
  module_task)

include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_executable(edist_utest
  euclidean_dist.c
  )
target_link_libraries(edist_utest PRIVATE tldevel module_param)
set_target_properties(edist_utest PROPERTIES COMPILE_FLAGS "-DUTEST_EDIST")
add_test(
  NAME edist_utest
  COMMAND edist_utest
  )


add_executable(bpm_utest
  bpm.c
  )
target_link_libraries(bpm_utest PRIVATE tldevel module_alphabet)
set_target_properties(bpm_utest PROPERTIES COMPILE_FLAGS "-DBPM_UTEST")
add_test(
  NAME bpm_utest
  COMMAND bpm_utest
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
