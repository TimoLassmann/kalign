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

add_executable(alphabet_utest
  alphabet.c   
  )
target_link_libraries(alphabet_utest PRIVATE tldevel)
set_target_properties(alphabet_utest PROPERTIES COMPILE_FLAGS "-DUTEST_ALPHABET")
add_test(
  NAME alphabet_utest
  COMMAND alphabet_utest
  )
   
