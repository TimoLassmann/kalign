add_library(module_alphabet STATIC
  alphabet.c 
  )

target_include_directories(module_alphabet
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
  )
target_link_libraries(module_alphabet PRIVATE tldevel)

include_directories(${CMAKE_CURRENT_BINARY_DIR})


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
