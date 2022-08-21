add_library(module_msa STATIC
  msa.c
  )

target_include_directories(module_msa
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
  )
target_link_libraries(module_msa PRIVATE
  tldevel
  module_alphabet)

include_directories(${CMAKE_CURRENT_BINARY_DIR})

