cmake_minimum_required(VERSION 3.18)

# Build library in sub-directories
if (NOT TARGET kalign::kalign)
  add_subdirectory(${PROJECT_SOURCE_DIR}/lib build EXCLUDE_FROM_ALL)
endif()

# link to installed library 
# if (NOT TARGET kalign::kalign)
#   find_package(kalign)
# endif()

if(NOT KALIGN_PACKAGE_NAME)
  set(KALIGN_PACKAGE_NAME ${CMAKE_PROJECT_NAME})
endif()

if(NOT KALIGN_PACKAGE_VERSION)
  set(KALIGN_PACKAGE_VERSION ${KALIGN_LIBRARY_VERSION_STRING})
endif()

configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/version.h.in"
  "version.h" )


add_executable(kalign-bin  
  run_kalign.c
  parameters.c
  
  )
# Link with static library.... 
target_link_libraries(kalign-bin PRIVATE tldevel ${PROJECT_NAME}_static )

set_target_properties(kalign-bin
  PROPERTIES OUTPUT_NAME kalign)



add_executable(kalignfmt
  
  run_reformat.c
  parameters.c
  )

install(TARGETS kalign-bin DESTINATION bin)

target_link_libraries(kalignfmt PRIVATE tldevel kalign::kalign)




