cmake_minimum_required(VERSION 3.18)

project(kalign LANGUAGES C CXX)

set(NAMESPACE_NAME "kalign")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(BUILD_SHARED_LIBS "Build the shared library" ON)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(GenerateExportHeader)

set(KALIGN_LIBRARY_VERSION_MAJOR 3)
set(KALIGN_LIBRARY_VERSION_MINOR 4)
set(KALIGN_LIBRARY_VERSION_PATCH 0)
set(KALIGN_LIBRARY_VERSION_STRING ${KALIGN_LIBRARY_VERSION_MAJOR}.${KALIGN_LIBRARY_VERSION_MINOR}.${KALIGN_LIBRARY_VERSION_PATCH})


set (CMAKE_C_STANDARD 11)

# SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pg")
# SET(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -pg")
# SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
# SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")

# to compile without open mp:
# cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DUSE_OPENMP=OFF ..

add_compile_options("$<$<CONFIG:RELEASE>:-W;-Wall;-O3;-pedantic>")
add_compile_options("$<$<CONFIG:DEBUG>:-W;-Wall;-O0;-g;-pedantic>")
add_compile_options("$<$<CONFIG:ASAN>:-W;-Wall;-Wextra;-O0;-g;-DDEBUG;-pedantic;-ffunction-sections;-fdata-sections;-fstack-protector-strong;-fsanitize=address>")

if(CMAKE_BUILD_TYPE MATCHES ASAN)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fsanitize=address")
endif()

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  if (MSVC)
    # warning level 4 and all warnings as errors
    add_compile_options(/W4 /WX)
  else()
    # lots of warnings and all warnings as errors
    add_compile_options(-Wall -Wextra -pedantic )
  endif()
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

include(GNUInstallDirs)

include(CTest)
include(CheckCSourceRuns)

option(USE_OPENMP    "Use OpenMP for parallelization" ON)
option(ENABLE_SSE    "Enable compile-time SSE4.1 support." ON)
option(ENABLE_AVX    "Enable compile-time AVX support."    ON)
option(ENABLE_AVX2   "Enable compile-time AVX2 support."   ON)
# option(ENABLE_FMA    "Enable compile-time FMA support."    ON)
# option(ENABLE_AVX512 "Enable compile-time AVX512 support." ON)


if(USE_OPENMP)
  find_package(OpenMP)
  if(OPENMP_FOUND OR OpenMP_FOUND)
    message(STATUS "OpenMP is enabled.")
     
    add_definitions (-DHAVE_OPENMP)
  else(OPENMP_FOUND OR OpenMP_FOUND)
    message(STATUS "OpenMP not supported")
  endif(OPENMP_FOUND OR OpenMP_FOUND)
endif(USE_OPENMP)


if (ENABLE_SSE)
  #
  # Check compiler for SSE4_1 intrinsics
  #
  if (CMAKE_COMPILER_IS_GNUCC OR (CMAKE_C_COMPILER_ID MATCHES "Clang") OR (CMAKE_CXX_COMPILER_ID MATCHES "Clang"))
    set(CMAKE_REQUIRED_FLAGS "-msse4.1")
    check_c_source_runs("
        #include <emmintrin.h>
        #include <smmintrin.h>
        int main()
        {
        __m128i a = _mm_setzero_si128();
        __m128i b = _mm_minpos_epu16(a);
        return 0;
        }"
      HAVE_SSE)
  endif()

  if (HAVE_SSE)        
    message(STATUS "SSE4.1 is enabled - target CPU must support it")
  endif()
  
  if (ENABLE_AVX)

    #
    # Check compiler for AVX intrinsics
    #
    if (CMAKE_COMPILER_IS_GNUCC OR (CMAKE_C_COMPILER_ID MATCHES "Clang") OR (CMAKE_CXX_COMPILER_ID MATCHES "Clang"))
      set(CMAKE_REQUIRED_FLAGS "-mavx")
      check_c_source_runs("
            #include <immintrin.h>
            int main()
            {
              __m256 a, b, c;
              const float src[8] = { 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f };
              float dst[8];
              a = _mm256_loadu_ps( src );
              b = _mm256_loadu_ps( src );
              c = _mm256_add_ps( a, b );
              _mm256_storeu_ps( dst, c );
              int i = 0;
              for( i = 0; i < 8; i++ ){
                if( ( src[i] + src[i] ) != dst[i] ){
                  return -1;
                }
              }
              return 0;
            }"
        HAVE_AVX)
    endif()

    if (HAVE_AVX)
      message(STATUS "AVX is enabled - target CPU must support it")
    endif()
  endif()
  
  if (ENABLE_AVX2)

    #
    # Check compiler for AVX intrinsics
    #
    if (CMAKE_COMPILER_IS_GNUCC OR (CMAKE_C_COMPILER_ID MATCHES "Clang") OR (CMAKE_CXX_COMPILER_ID MATCHES "Clang"))
      set(CMAKE_REQUIRED_FLAGS "-mavx2")
      check_c_source_runs("
          #include <immintrin.h>
          int main()
          {
            __m256i a, b, c;
            const int src[8] = { 1, 2, 3, 4, 5, 6, 7, 8 };
            int dst[8];
            a =  _mm256_loadu_si256( (__m256i*)src );
            b =  _mm256_loadu_si256( (__m256i*)src );
            c = _mm256_add_epi32( a, b );
            _mm256_storeu_si256( (__m256i*)dst, c );
            int i = 0;
            for( i = 0; i < 8; i++ ){
              if( ( src[i] + src[i] ) != dst[i] ){
                return -1;
              }
            }
            return 0;
          }"
        HAVE_AVX2)
    endif()

    if (HAVE_AVX2)
      message(STATUS "AVX2 is enabled - target CPU must support it")
    endif()
  endif()
endif()

if (HAVE_AVX2)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2 -DHAVE_AVX2")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mavx2 -DHAVE_AVX2")
else(HAVE_AVX2)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNOHAVE_AVX2")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DNOHAVE_AVX2")
endif(HAVE_AVX2)

add_subdirectory(lib)
add_subdirectory(src)
add_subdirectory(tests)


MESSAGE(STATUS "")
MESSAGE(STATUS "Configuration:")
MESSAGE(STATUS "--------------------------------------")
MESSAGE(STATUS "Build type                : " ${CMAKE_BUILD_TYPE})
MESSAGE(STATUS "Compiler flags            : " ${CMAKE_C_COMPILE_FLAGS})
MESSAGE(STATUS "Compiler c debug flags    : " ${CMAKE_C_FLAGS_DEBUG})
MESSAGE(STATUS "Compiler c release flags  : " ${CMAKE_C_FLAGS_RELEASE})
MESSAGE(STATUS "Compiler c min size flags : " ${CMAKE_C_FLAGS_MINSIZEREL})
MESSAGE(STATUS "Compiler c flags          : " ${CMAKE_C_FLAGS})
message(STATUS "OpenMP version            : " ${OpenMP_C_VERSION})
message(STATUS "OpenMP flags              : " ${OpenMP_C_FLAGS})

# Package Generator  #######################################################
set(CPACK_PACKAGE_VENDOR "Timo Lassmann")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Kalign")
set(CPACK_PACKAGE_VERSION_MAJOR ${KALIGN_LIBRARY_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${KALIGN_LIBRARY_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${KALIGN_LIBRARY_VERSION_PATCH})
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")
# set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README.org")
set(CPACK_SOURCE_GENERATOR "TGZ;ZIP")
set(CPACK_SOURCE_IGNORE_FILES
  /.cache
  /.git
  /GPATH
  /GTAGS
  /GRTAGS
  /.*build.*
  /.dir-locals.el
  /\\\\.DS_Store
  )
include (CPack)


