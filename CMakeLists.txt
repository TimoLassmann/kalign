cmake_minimum_required(VERSION 3.18)

project(kalign LANGUAGES C CXX)

set(NAMESPACE_NAME "kalign")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(BUILD_SHARED_LIBS "Build the shared library" ON)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(GenerateExportHeader)

set(KALIGN_LIBRARY_VERSION_MAJOR 4)
set(KALIGN_LIBRARY_VERSION_MINOR 0)
set(KALIGN_LIBRARY_VERSION_PATCH 0)
set(KALIGN_LIBRARY_VERSION_STRING ${KALIGN_LIBRARY_VERSION_MAJOR}.${KALIGN_LIBRARY_VERSION_MINOR}.${KALIGN_LIBRARY_VERSION_PATCH})

# Shared-library ABI version, independent of the project version. Bump it
# whenever the C ABI breaks. 4: kalign 3.7.0 added fields to the public
# kalign_run_config / kalign_ensemble_config structs, which
# kalign_run_config_defaults() / kalign_ensemble_config_defaults() return by
# value (the 3.x releases before 3.7.0 used SOVERSION 3).
set(KALIGN_LIBRARY_SOVERSION 4)


set (CMAKE_C_STANDARD 11)

# SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pg")
# SET(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -pg")
# SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
# SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")

# Compiler-specific flags
add_compile_options("$<$<CONFIG:RELEASE>:-W;-Wall;-O3;-pedantic>")
add_compile_options("$<$<CONFIG:DEBUG>:-W;-Wall;-O0;-g;-pedantic>")
add_compile_options("$<$<CONFIG:ASAN>:-W;-Wall;-Wextra;-O0;-g;-DDEBUG;-pedantic;-ffunction-sections;-fdata-sections;-fstack-protector-strong;-fsanitize=address;-DKALIGN_DEBUG_ALIGNED>")

if(CMAKE_BUILD_TYPE MATCHES ASAN)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fsanitize=address")
endif()

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  # lots of warnings
  add_compile_options(-Wall -Wextra -pedantic )
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

include(GNUInstallDirs)

include(CTest)

# ENABLE_AVX2: include the AVX2 kernels in x86-64 builds. They are chosen at
# run time (kalign_simd_on in lib/src/kalign_simd.h); the program itself is
# compiled for the plain x86-64 instruction set, so no flag here raises the
# instruction set. OFF leaves the AVX2 kernels out (-DKALIGN_NO_AVX2).
option(ENABLE_AVX2   "Include the AVX2 kernels (selected at run time on x86-64)." ON)

# Performance tuning parameters
# Threadpool parallelization thresholds (all settings in one place)
set(KALIGN_ALN_SERIAL_THRESHOLD "500" CACHE STRING "Hirschberg DP width below which to run serial")
set(KALIGN_KMEANS_UPGMA_THRESHOLD "50" CACHE STRING "Sequence count below which to use UPGMA instead of parallel k-means")
set(KALIGN_DIST_MIN_SEQS "50" CACHE STRING "Minimum sequences to parallelize distance/anchor computations")
set(KALIGN_PFOR_MIN_CHUNK "10" CACHE STRING "Minimum iterations per chunk in tp_parallel_for")

# Parallelism: built-in work-stealing threadpool (pthreads). With one thread
# no pool is created and every task runs on the calling thread.
find_package(Threads REQUIRED)

if (NOT ENABLE_AVX2)
  add_compile_options(-DKALIGN_NO_AVX2)
  message(STATUS "AVX2 kernels are left out (ENABLE_AVX2=OFF)")
endif()

# Add performance tuning parameters as compile definitions
add_definitions(-DKALIGN_ALN_SERIAL_THRESHOLD=${KALIGN_ALN_SERIAL_THRESHOLD})
add_definitions(-DKALIGN_KMEANS_UPGMA_THRESHOLD=${KALIGN_KMEANS_UPGMA_THRESHOLD})
add_definitions(-DKALIGN_DIST_MIN_SEQS=${KALIGN_DIST_MIN_SEQS})
add_definitions(-DKALIGN_PFOR_MIN_CHUNK=${KALIGN_PFOR_MIN_CHUNK})
option(KALIGN_WAVE "Wavefront SIMD (NEON/AVX2) DP passes, bit-identical to the scalar kernels (lib/src/aln_wave.c)" ON)
if(KALIGN_WAVE)
  add_definitions(-DKALIGN_WAVE)
endif()

add_subdirectory(lib)
add_subdirectory(src)
add_subdirectory(tests)

# DP performance harness (research only, off by default): builds an
# instrumented library variant plus profiling/benchmark tools in perf/dp.
# The regular kalign library and binary are not affected by this option.
option(KALIGN_PERF "Build the perf/dp profiling and DP benchmark harness" OFF)
if(KALIGN_PERF)
  add_subdirectory(perf/dp)
endif()

# Python module build (optional) ##########################################
option(BUILD_PYTHON_MODULE "Build Python extension module" OFF)

if(BUILD_PYTHON_MODULE)
    # Find pybind11 - scikit-build-core will provide this
    find_package(pybind11 CONFIG)

    if(pybind11_FOUND)
        message(STATUS "Building Python extension module")

        # Create the Python extension module
        pybind11_add_module(_core
            python-kalign/_core.cpp
            tests/dssim.c
        )

        # Link against the static kalign library
        target_link_libraries(_core PRIVATE kalign_static)

        target_link_libraries(_core PRIVATE Threads::Threads)

        # Add include directories for DSSim and lib headers
        target_include_directories(_core PRIVATE
            tests
            lib/src
            lib/include
        )

        # Set properties for the Python module
        set_target_properties(_core PROPERTIES
            CXX_STANDARD 11
            CXX_STANDARD_REQUIRED ON
            CXX_EXTENSIONS OFF
        )

        # Install only the Python module (no shared library needed with static linking)
        install(TARGETS _core
                LIBRARY DESTINATION kalign
                COMPONENT python)

        # Install Python source files
        install(DIRECTORY python-kalign/
                DESTINATION kalign
                COMPONENT python
                FILES_MATCHING PATTERN "*.py" PATTERN "*.typed")

        message(STATUS "Python extension module configured")
    else()
        message(WARNING "pybind11 not found - Python module will not be built")
    endif()
else()
    message(STATUS "Python module build disabled (use -DBUILD_PYTHON_MODULE=ON to enable)")
endif()

# Benchmark target ############################################################
# Runs BAliBASE benchmarks comparing the C binary and the Python API.
# Requires: pip install -e .  (to make the kalign Python package available)
#
# Usage:
#   make benchmark                    # full BAliBASE suite, both methods
#   make benchmark BENCH_MAX_CASES=5  # quick smoke test with 5 cases
#
find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_FOUND)
    set(BENCH_MAX_CASES "0" CACHE STRING "Max benchmark cases (0 = all)")
    add_custom_target(benchmark
        COMMAND ${CMAKE_COMMAND} -E echo "Running kalign benchmarks..."
        COMMAND ${Python3_EXECUTABLE} -m benchmarks
            --dataset balibase
            --method python_api cli
            --binary $<TARGET_FILE:kalign-bin>
            --max-cases ${BENCH_MAX_CASES}
            --output ${CMAKE_SOURCE_DIR}/benchmarks/results/latest.json
            -v
        DEPENDS kalign-bin
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMENT "Benchmarking kalign (C binary vs Python API)"
        VERBATIM
    )
endif()

MESSAGE(STATUS "")
MESSAGE(STATUS "Configuration:")
MESSAGE(STATUS "--------------------------------------")
MESSAGE(STATUS "Build type                : " ${CMAKE_BUILD_TYPE})
MESSAGE(STATUS "Compiler flags            : " ${CMAKE_C_COMPILE_FLAGS})
MESSAGE(STATUS "Compiler c debug flags    : " ${CMAKE_C_FLAGS_DEBUG})
MESSAGE(STATUS "Compiler c release flags  : " ${CMAKE_C_FLAGS_RELEASE})
MESSAGE(STATUS "Compiler c min size flags : " ${CMAKE_C_FLAGS_MINSIZEREL})
MESSAGE(STATUS "Compiler c flags          : " ${CMAKE_C_FLAGS})

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
