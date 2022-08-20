

add_executable(genruncmds
  create_run_cmds.c
  )
target_link_libraries(genruncmds PRIVATE
  tldevel
  )

add_executable(genkaligndiffcmds
  create_run_kalign_diff_cmds.c
  )
target_link_libraries(genkaligndiffcmds PRIVATE
  tldevel
  )


add_executable(benchrunner 
  benchrunner.c
  )
if(OpenMP_C_FOUND)
  target_link_libraries(benchrunner PUBLIC OpenMP::OpenMP_C)
endif()

target_link_libraries(benchrunner PRIVATE
  tldevel
  module_param
  module_alphabet
  module_msa
  )

add_executable(timescorealn
  time_and_score_aln.c
)

target_link_libraries(timescorealn PRIVATE
  tldevel
  module_param
  module_alphabet
  module_tree
  module_aln)


