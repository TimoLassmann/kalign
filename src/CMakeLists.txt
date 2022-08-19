
add_subdirectory(tldevel)
add_subdirectory(interface)
add_subdirectory(msa_io)
add_subdirectory(alphabet)
add_subdirectory(tasks)
add_subdirectory(tree)
add_subdirectory(alignment)


add_executable(kalign run_kalign.c)
if(OpenMP_C_FOUND)
    target_link_libraries(kalign PUBLIC OpenMP::OpenMP_C)
endif()

target_link_libraries(kalign PRIVATE
  tldevel
  module_param
  module_alphabet
  module_tree
  module_aln)

 

# add_executable(create_bench create_test_cases.c)
# target_link_libraries(create_bench PRIVATE tld-dev tools)

# add_executable(bench_setup  benchmark_setup.c )
# target_link_libraries(bench_setup PRIVATE tld-dev)


# add_executable(agp_predict main_predict.c)
# target_link_libraries(agp_predict PRIVATE
#   module_io
#   module_agp_core
#   tld-dev
#   tools
#   Threads::Threads
#   )



# set_target_properties(agp_gp PROPERTIES COMPILE_FLAGS "-DGP")


