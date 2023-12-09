cmake_minimum_required(VERSION 3.18)

if(NOT KALIGN_PACKAGE_NAME)
  set(KALIGN_PACKAGE_NAME ${CMAKE_PROJECT_NAME})
endif()

if(NOT KALIGN_PACKAGE_VERSION)
  set(KALIGN_PACKAGE_VERSION ${KALIGN_LIBRARY_VERSION_STRING})
endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")


configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/src/version.h.in" "version.h")


set(source_files
  src/test.c
  src/tldevel.c
  src/tlmisc.c
  src/tlrng.c
  src/esl_stopwatch.c
  src/msa_alloc.c
  src/msa_op.c
  src/msa_io.c
  src/msa_misc.c
  src/msa_check.c
  src/msa_cmp.c
  src/msa_sort.c
  src/alphabet.c

  src/task.c

  src/bisectingKmeans.c
  src/sequence_distance.c
  src/bpm.c
  src/euclidean_dist.c
  src/pick_anchor.c
  
  src/aln_wrap.c

  src/aln_param.c
  src/aln_run.c
  src/aln_mem.c
  src/aln_setup.c
  src/aln_controller.c
  src/aln_seqseq.c
  src/aln_seqprofile.c
  src/aln_profileprofile.c
  src/weave_alignment.c

  # src/coretralign.c
  # src/test.h
  )

add_library(${PROJECT_NAME}_OBJ OBJECT ${source_files})

if(OpenMP_C_FOUND)  
  target_link_libraries(${PROJECT_NAME}_OBJ PRIVATE OpenMP::OpenMP_C)
endif(OpenMP_C_FOUND)


target_include_directories(${PROJECT_NAME}_OBJ
  PRIVATE
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  )

add_library(${PROJECT_NAME})
add_library(${NAMESPACE_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})



generate_export_header(${PROJECT_NAME}
  EXPORT_FILE_NAME ${PROJECT_NAME}/${PROJECT_NAME}_export.h
  )

set(public_header_files
  ${CMAKE_SOURCE_DIR}/lib/include/kalign/kalign.h
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}/${PROJECT_NAME}_export.h
  #${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}/${PROJECT_NAME}_export.h
  )

set_target_properties(${PROJECT_NAME} PROPERTIES
  LANGUAGES C
  CXX_STANDARD 11
  CXX_STANDARD_REQUIRED ON
  CXX_EXTENSIONS OFF
  VERSION ${KALIGN_LIBRARY_VERSION_STRING}
  SOVERSION ${KALIGN_LIBRARY_VERSION_MAJOR}
  PUBLIC_HEADER "${public_header_files}"
  CXX_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN 1
  )

target_sources(${PROJECT_NAME}
  PRIVATE
  ${public_header_files}
  $<TARGET_OBJECTS:${PROJECT_NAME}_OBJ>        
  )

# Define headers for this library. PUBLIC headers are used for
# compiling the library, and will be added to consumers' build
# paths.
target_include_directories(${PROJECT_NAME}
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>  
  )

target_include_directories(${PROJECT_NAME}
  PRIVATE
  src
  
  )

target_link_libraries(${PROJECT_NAME} PRIVATE m)

if(OpenMP_C_FOUND)
  target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_C)
endif()



add_library(${PROJECT_NAME}_static STATIC
  ${public_header_files}
  $<TARGET_OBJECTS:${PROJECT_NAME}_OBJ>        
  )

set_target_properties(${PROJECT_NAME}_static PROPERTIES
  LANGUAGES C
  CXX_STANDARD 11
  CXX_STANDARD_REQUIRED ON
  CXX_EXTENSIONS OFF
  VERSION ${KALIGN_LIBRARY_VERSION_STRING}
  SOVERSION ${KALIGN_LIBRARY_VERSION_MAJOR}
  PUBLIC_HEADER "${public_header_files}"
  CXX_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN 1
  )

target_include_directories(${PROJECT_NAME}_static
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>  
  )

target_include_directories(${PROJECT_NAME}_static
  PRIVATE
  src
  )

target_link_libraries(${PROJECT_NAME}_static PRIVATE m)

if(OpenMP_C_FOUND)
  target_link_libraries(${PROJECT_NAME}_static PRIVATE OpenMP::OpenMP_C)
endif()


add_library(tldevel STATIC
  src/tldevel.c
  src/tlmisc.c
  src/tlrng.c
  src/esl_stopwatch.c
  )

target_include_directories(tldevel 
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src
  )
target_link_libraries(tldevel PRIVATE m)

# adjust this path depending where cmake searches for the target files
set(ConfigPackageLocation "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")

# 'make install' to the correct locations (provided by GNUInstallDirs).
install(TARGETS ${PROJECT_NAME}
  EXPORT ${PROJECT_NAME}Targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Library
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Library
    
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Library # This is for Windows
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME} COMPONENT Development
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )

install(EXPORT ${PROJECT_NAME}Targets
  FILE
    ${PROJECT_NAME}Targets.cmake
  NAMESPACE
    ${NAMESPACE_NAME}::
  DESTINATION
    ${ConfigPackageLocation}
  COMPONENT
    Development
)

configure_package_config_file(
    ${PROJECT_NAME}Config.cmake.in
    ${PROJECT_NAME}Config.cmake
    INSTALL_DESTINATION "${ConfigPackageLocation}"
    PATH_VARS CMAKE_INSTALL_PREFIX
    )
  
write_basic_package_version_file(
    ${PROJECT_NAME}ConfigVersion.cmake
    VERSION ${KALIGN_LIBRARY_VERSION_STRING}
    COMPATIBILITY AnyNewerVersion
    )

install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
  DESTINATION
    ${ConfigPackageLocation}
  COMPONENT
    Development
)


add_executable(bpm_utest
  src/tldevel.c
  src/tlrng.c
  src/esl_stopwatch.c
  src/alphabet.c
  src/bpm.c
  src/bpm_test.c 
  
  
  )
target_link_libraries(bpm_utest PRIVATE m)

add_test(
  NAME bpm_utest
  COMMAND bpm_utest
  )

add_executable(alphabet_utest
  src/tldevel.c
  src/alphabet.c   
  )

set_target_properties(alphabet_utest PROPERTIES COMPILE_FLAGS "-DUTEST_ALPHABET")
add_test(
  NAME alphabet_utest
  COMMAND alphabet_utest
  )
   
add_executable(edist_utest
  src/tldevel.c
  src/tlrng.c
  src/esl_stopwatch.c
  src/euclidean_dist.c
  )


target_link_libraries(edist_utest PRIVATE m)
set_target_properties(edist_utest PROPERTIES COMPILE_FLAGS "-DUTEST_EDIST")
add_test(
  NAME edist_utest
  COMMAND edist_utest
  )


add_executable(task_utest
  src/tldevel.c
  src/tlrng.c
  src/task.c  
  )


target_link_libraries(task_utest PRIVATE m)
set_target_properties(task_utest PROPERTIES COMPILE_FLAGS "-DTASKWRITETEST")
add_test(
  NAME task_utest
  COMMAND task_utest
  )
 
