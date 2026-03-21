cmake_minimum_required(VERSION 3.18)

# Support standalone build
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    project(threadpool LANGUAGES C)
    set(CMAKE_C_STANDARD 11)
    set(CMAKE_C_STANDARD_REQUIRED ON)
    enable_testing()

    # Optional: find OpenMP for benchmark comparison
    find_package(OpenMP)
endif()

find_package(Threads REQUIRED)

# ── Thread pool library ───────────────────────────────────────────
add_library(threadpool STATIC threadpool.c)
target_include_directories(threadpool PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(threadpool PUBLIC Threads::Threads)
set_target_properties(threadpool PROPERTIES C_STANDARD 11)

# ── Tests ─────────────────────────────────────────────────────────
add_executable(test_threadpool test_threadpool.c)
target_link_libraries(test_threadpool PRIVATE threadpool)
set_target_properties(test_threadpool PROPERTIES C_STANDARD 11)
add_test(NAME threadpool_tests COMMAND test_threadpool)

add_executable(stress_threadpool stress_threadpool.c)
target_link_libraries(stress_threadpool PRIVATE threadpool)
set_target_properties(stress_threadpool PROPERTIES C_STANDARD 11)
add_test(NAME threadpool_stress COMMAND stress_threadpool)

# ── Benchmarks ────────────────────────────────────────────────────
add_executable(bench_threadpool bench_threadpool.c)
target_link_libraries(bench_threadpool PRIVATE threadpool)
set_target_properties(bench_threadpool PROPERTIES C_STANDARD 11)

if(OpenMP_C_FOUND)
    target_link_libraries(bench_threadpool PRIVATE OpenMP::OpenMP_C)
    target_compile_definitions(bench_threadpool PRIVATE HAVE_OPENMP)
endif()
