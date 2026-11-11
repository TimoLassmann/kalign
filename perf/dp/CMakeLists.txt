# perf/dp -- DP performance harness (built only with -DKALIGN_PERF=ON).
#
# Targets
#   kalign_perf_static  instrumented copy of the kalign library
#                       (KALIGN_PERF_TIMERS: region timers + KALIGN_DP_DUMP)
#   kalign_prof         driver: read -> mode preset -> align -> write, prints
#                       a JSON time breakdown (links kalign_perf_static)
#   kalign_plain        same driver on the regular library (timing/output parity)
#   dssim_gen           FASTA generator using Kalign's DSSim simulator
#   bench_dp            DP microbenchmark + oracle over dumped .kdp problems
#                       (links the regular, uninstrumented kalign_static)
cmake_minimum_required(VERSION 3.18)

set(KPERF_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(KLIB_DIR ${PROJECT_SOURCE_DIR}/lib)

# --- instrumented library variant ---------------------------------------
get_target_property(_kalign_obj_sources ${PROJECT_NAME}_OBJ SOURCES)
set(_kperf_lib_sources "")
foreach(_s ${_kalign_obj_sources})
  if(IS_ABSOLUTE ${_s})
    list(APPEND _kperf_lib_sources ${_s})
  else()
    list(APPEND _kperf_lib_sources ${KLIB_DIR}/${_s})
  endif()
endforeach()
list(APPEND _kperf_lib_sources ${KPERF_DIR}/instr/kperf.c)

add_library(kalign_perf_static STATIC ${_kperf_lib_sources})
target_compile_definitions(kalign_perf_static
  PUBLIC KALIGN_PERF_TIMERS
  PRIVATE KALIGN_PACKAGE_VERSION=\"${KALIGN_LIBRARY_VERSION_STRING}\")
target_include_directories(kalign_perf_static
  PUBLIC
    ${KLIB_DIR}/include
    ${CMAKE_BINARY_DIR}/lib
    ${KLIB_DIR}/src
    ${KPERF_DIR}/instr
  PRIVATE
    ${KLIB_DIR}/src/threadpool)
target_link_libraries(kalign_perf_static PUBLIC m)
target_link_libraries(kalign_perf_static PUBLIC Threads::Threads)
add_dependencies(kalign_perf_static ${PROJECT_NAME}_static)

# --- profiling drivers ---------------------------------------------------
add_executable(kalign_prof tools/kalign_prof.c)
target_link_libraries(kalign_prof PRIVATE kalign_perf_static)

add_executable(kalign_plain tools/kalign_prof.c)
target_include_directories(kalign_plain PRIVATE ${KLIB_DIR}/src ${KLIB_DIR}/include ${CMAKE_BINARY_DIR}/lib)
target_link_libraries(kalign_plain PRIVATE ${PROJECT_NAME}_static)
target_link_libraries(kalign_plain PRIVATE Threads::Threads)

# --- dataset generator -----------------------------------------------------
add_executable(dssim_gen tools/dssim_gen.c ${PROJECT_SOURCE_DIR}/tests/dssim.c)
target_include_directories(dssim_gen PRIVATE ${PROJECT_SOURCE_DIR}/tests ${KLIB_DIR}/src
  ${KLIB_DIR}/include ${CMAKE_BINARY_DIR}/lib)
target_link_libraries(dssim_gen PRIVATE ${PROJECT_NAME}_static m)
target_link_libraries(dssim_gen PRIVATE Threads::Threads)

# --- DP microbenchmark ------------------------------------------------------
# Every .c file in kernels/ is compiled in; kernels self-register with
# KDP_REGISTER_KERNEL (see bench/kdp_kernel.h).  Re-run cmake after adding
# a file (CONFIGURE_DEPENDS makes the glob re-evaluate on build).
file(GLOB KDP_KERNEL_SOURCES CONFIGURE_DEPENDS ${KPERF_DIR}/kernels/*.c)
#
# The reference kernels are compiled into bench_dp from lib/src (same flags as
# the library plus -falign-functions=64) instead of being taken from
# kalign_static: on Apple M-series the speed of the DP inner loops depends on
# their code alignment (a 1.84 vs 2.53 ns/cell swing for the seq-seq forward
# pass was observed purely from link layout), so every kernel, reference and
# prototype, is function-aligned to keep timings comparable across rebuilds.
set(KDP_REF_SOURCES
  ${KLIB_DIR}/src/aln_seqseq.c
  ${KLIB_DIR}/src/aln_seqprofile.c
  ${KLIB_DIR}/src/aln_profileprofile.c)
add_executable(bench_dp
  bench/bench_dp.c
  ${KDP_REF_SOURCES}
  bench/kdp_io.c
  bench/kdp_hirsch.c
  bench/kdp_score.c
  bench/kdp_registry.c
  ${KDP_KERNEL_SOURCES})
target_include_directories(bench_dp PRIVATE
  ${KPERF_DIR}/bench ${KPERF_DIR}/instr
  ${KLIB_DIR}/src ${KLIB_DIR}/include ${CMAKE_BINARY_DIR}/lib)
set(KDP_ALIGN_FLAGS "-falign-functions=64" CACHE STRING "code alignment flags for bench_dp kernels")
separate_arguments(_kdp_align UNIX_COMMAND "${KDP_ALIGN_FLAGS}")
target_compile_options(bench_dp PRIVATE ${_kdp_align})
target_link_libraries(bench_dp PRIVATE ${PROJECT_NAME}_static m)
target_link_libraries(bench_dp PRIVATE Threads::Threads)
