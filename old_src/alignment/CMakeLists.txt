add_library(module_aln STATIC
  aln_profileprofile.c
  aln_seqprofile.c
  aln_seqseq.c
  aln_mem.c
  aln_setup.c  
  aln_controller.c
  aln_run.c
  aln_param.c
  weave_alignment.c
  )

target_include_directories(module_aln 
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
  )

target_link_libraries(module_aln PRIVATE
  tldevel
  module_msa
  module_param
  module_task
  module_alphabet)

if(OpenMP_CXX_FOUND)
    target_link_libraries(module_aln PUBLIC OpenMP::OpenMP_C)
endif()

  
# ELSE(OPENMP_FOUND AND WITH_OPENMP)
  
# ENDIF(OPENMP_FOUND AND WITH_OPENMP)
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
