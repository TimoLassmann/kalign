cmake_minimum_required(VERSION 3.18)

# Build library in sub-directories
if (NOT TARGET kalign::kalign)
  add_subdirectory(${PROJECT_SOURCE_DIR}/lib build EXCLUDE_FROM_ALL)
endif()

# link to installed library 
if (NOT TARGET kalign::kalign)
  find_package(kalign)
endif()

add_executable(largebench
  large_benchmark.c
  )
target_link_libraries(largebench tldevel   ${PROJECT_NAME}_static)

add_executable(kalign_io_test
  kalign_io_test.c
  )

target_link_libraries(kalign_io_test kalign::kalign)


add_executable(dssim
  dssim_test.c
  dssim.c
  )
target_link_libraries(dssim tldevel kalign::kalign)




add_executable(kalign_lib_test
  kalign_lib_test.c
  )
target_link_libraries(kalign_lib_test kalign::kalign)

add_executable(kalign_cmp_test
  kalign_cmp_test.c
  )
target_link_libraries(kalign_cmp_test kalign::kalign)

add_executable(kaligncpp
   kalign_lib_testCXX.cpp)
target_link_libraries(kaligncpp kalign::kalign)

add_test(
  NAME C++Lib
  COMMAND kaligncpp
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )


#BB11001.tfa     
add_test(
  NAME kalign_itest_BB11001.tfa
  COMMAND kalign-bin ${CMAKE_CURRENT_SOURCE_DIR}/data/BB11001.tfa
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )

#BB11001_EOF.msf
add_test(
  NAME kalign_itest_BB11001_EOF.msf
  COMMAND kalign-bin ${CMAKE_CURRENT_SOURCE_DIR}/data/BB11001_EOF.msf
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )

#BB12006.msf
add_test(
  NAME kalign_itest_BB12006.msf
  COMMAND kalign-bin ${CMAKE_CURRENT_SOURCE_DIR}/data/BB12006.msf
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )

#BB12006.tfa
add_test(
  NAME kalign_itest_BB12006.tfa
  COMMAND kalign-bin ${CMAKE_CURRENT_SOURCE_DIR}/data/BB12006.tfa
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )

#BB30014.msf     
add_test(
  NAME kalign_itest_BB30014.msf
  COMMAND kalign-bin ${CMAKE_CURRENT_SOURCE_DIR}/data/BB12006.tfa
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )

#BB30014.tfa
add_test(
  NAME kalign_itest_BB30014.tfa
  COMMAND kalign-bin ${CMAKE_CURRENT_SOURCE_DIR}/data/BB30014.tfa
  WORKING_DIRECTORY 
${CMAKE_CURRENT_SOURCE_DIR}
  )


add_test(
  NAME DSSIM_takes_a_min
  COMMAND dssim
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
