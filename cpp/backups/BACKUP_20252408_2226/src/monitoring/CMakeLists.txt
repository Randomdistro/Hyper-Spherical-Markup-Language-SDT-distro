# [Enterprise Bean]: CMAKE CONFIGURATION FOR PERFORMANCE MONITORING SYSTEMS!
# [Performance Demon]: Optimized builds with all the compiler flags!

cmake_minimum_required(VERSION 3.20)

# [Modern Hipster]: Using latest CMake features!
project(HSMLPerformanceMonitoring
    VERSION 1.0.0
    DESCRIPTION "HSML Multi-Personality Performance Monitoring System"
    LANGUAGES CXX)

# [Performance Demon]: MAXIMUM OPTIMIZATION FLAGS!
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# [Security Paranoid]: Compiler hardening flags!
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native -flto")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g3 -fsanitize=address -fsanitize=undefined")

# [Performance Demon]: Architecture-specific optimizations!
if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2 -mfma -mbmi2")
    
    # [Hacktivist]: Check for AVX-512 support!
    include(CheckCXXCompilerFlag)
    check_cxx_compiler_flag("-mavx512f" COMPILER_SUPPORTS_AVX512)
    if(COMPILER_SUPPORTS_AVX512)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx512f -mavx512dq")
        add_definitions(-DHSML_AVX512_SUPPORT)
    endif()
endif()

# [Security Paranoid]: Thread safety and security flags!
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -fstack-protector-strong")

# [Enterprise Bean]: Include directories configuration!
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/../../include
    ${CMAKE_CURRENT_SOURCE_DIR}/../../include/hsml
    ${CMAKE_CURRENT_SOURCE_DIR}/../../include/hsml/monitoring
)

# [Functional Purist]: Pure header-only monitoring libraries!
set(MONITORING_HEADERS
    ../../include/hsml/monitoring/hardware_performance_monitor.h
    ../../include/hsml/monitoring/simd_batch_monitor.h
    ../../include/hsml/monitoring/realtime_telemetry_monitor.h
    ../../include/hsml/monitoring/statistical_analysis_monitor.h
    ../../include/hsml/monitoring/low_overhead_profiler.h
    ../../include/hsml/monitoring/unified_performance_monitor.h
)

# [Modern Hipster]: Create interface library for headers!
add_library(hsml_monitoring_headers INTERFACE)
target_include_directories(hsml_monitoring_headers INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/../../include
)

# [Performance Demon]: Core HSML dependencies!
find_library(HSML_CORE_LIB hsml_core PATHS ${CMAKE_CURRENT_SOURCE_DIR}/../../build)
if(NOT HSML_CORE_LIB)
    message(WARNING "HSML core library not found. Building monitoring systems only.")
endif()

# [Hacktivist]: System dependencies for hardware monitoring!
find_package(Threads REQUIRED)

# [Security Paranoid]: Check for required system capabilities!
include(CheckIncludeFile)
check_include_file("linux/perf_event.h" HAVE_PERF_EVENT_H)
check_include_file("x86intrin.h" HAVE_X86INTRIN_H)

if(HAVE_PERF_EVENT_H)
    add_definitions(-DHSML_HAVE_PERF_EVENTS)
endif()

if(HAVE_X86INTRIN_H)
    add_definitions(-DHSML_HAVE_X86_INTRINSICS)
endif()

# [Performance Demon]: Performance monitoring demo executable!
add_executable(hsml_performance_monitor_demo
    performance_monitor_demo.cpp
    ${MONITORING_HEADERS}
)

target_link_libraries(hsml_performance_monitor_demo
    hsml_monitoring_headers
    Threads::Threads
)

# [Modern Hipster]: Link HSML core if available!
if(HSML_CORE_LIB)
    target_link_libraries(hsml_performance_monitor_demo ${HSML_CORE_LIB})
endif()

# [Enterprise Bean]: Compiler-specific optimizations!
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(hsml_performance_monitor_demo PRIVATE
        -Wall -Wextra -Wpedantic
        -Wno-unused-parameter  # [Minimalist Zen]: Sometimes simplicity means unused params
        -ffast-math            # [Performance Demon]: SPEED!
    )
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    target_compile_options(hsml_performance_monitor_demo PRIVATE
        -Wall -Wextra -Wpedantic
        -Wno-unused-parameter
        -ffast-math
        -mllvm -enable-loop-vectorize  # [Hacktivist]: Force vectorization!
    )
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_compile_options(hsml_performance_monitor_demo PRIVATE
        /W4 /O2 /arch:AVX2
        /D_WIN32_WINNT=0x0601  # [Security Paranoid]: Windows 7+ for performance counters
    )
endif()

# [Functional Purist]: Unit tests for pure monitoring functions!
option(BUILD_MONITORING_TESTS "Build monitoring system tests" ON)

if(BUILD_MONITORING_TESTS)
    # [Enterprise Bean]: Test discovery and execution framework!
    enable_testing()
    
    # [Security Paranoid]: Integrity validation tests!
    find_package(GTest QUIET)
    if(GTest_FOUND)
        add_executable(hsml_monitoring_tests
            tests/test_hardware_monitor.cpp
            tests/test_simd_batch_monitor.cpp
            tests/test_statistical_analysis.cpp
            tests/test_low_overhead_profiler.cpp
            tests/test_unified_monitor.cpp
        )
        
        target_link_libraries(hsml_monitoring_tests
            hsml_monitoring_headers
            GTest::gtest
            GTest::gtest_main
            Threads::Threads
        )
        
        # [Functional Purist]: Pure test discovery!
        include(GoogleTest)
        gtest_discover_tests(hsml_monitoring_tests)
    endif()
endif()

# [Modern Hipster]: Benchmarking targets!
option(BUILD_MONITORING_BENCHMARKS "Build performance benchmarks" ON)

if(BUILD_MONITORING_BENCHMARKS)
    find_package(benchmark QUIET)
    if(benchmark_FOUND)
        add_executable(hsml_monitoring_benchmarks
            benchmarks/benchmark_monitoring_overhead.cpp
            benchmarks/benchmark_simd_processing.cpp
            benchmarks/benchmark_statistical_analysis.cpp
        )
        
        target_link_libraries(hsml_monitoring_benchmarks
            hsml_monitoring_headers
            benchmark::benchmark
            Threads::Threads
        )
    endif()
endif()

# [Enterprise Bean]: Installation configuration!
install(TARGETS hsml_performance_monitor_demo
    RUNTIME DESTINATION bin
    COMPONENT monitoring_tools
)

install(FILES ${MONITORING_HEADERS}
    DESTINATION include/hsml/monitoring
    COMPONENT monitoring_headers
)

# [Security Paranoid]: Package configuration!
set(CPACK_PACKAGE_NAME "HSML-Performance-Monitoring")
set(CPACK_PACKAGE_VERSION "${PROJECT_VERSION}")
set(CPACK_PACKAGE_DESCRIPTION "${PROJECT_DESCRIPTION}")
set(CPACK_PACKAGE_CONTACT "HSML Development Team")

include(CPack)

# [All Personalities]: Build summary message!
message(STATUS "🎭 HSML Multi-Personality Performance Monitoring Configuration:")
message(STATUS "🔥 [Performance Demon]: Optimizations enabled: ${CMAKE_CXX_FLAGS_RELEASE}")
message(STATUS "🚀 [Hacktivist]: SIMD support: AVX2=${HAVE_X86INTRIN_H}, AVX512=${COMPILER_SUPPORTS_AVX512}")
message(STATUS "🔒 [Security Paranoid]: Hardware counters: ${HAVE_PERF_EVENT_H}")
message(STATUS "🧘 [Minimalist Zen]: Header-only design for simplicity")
message(STATUS "💫 [Modern Hipster]: C++20 features enabled")
message(STATUS "📊 [Functional Purist]: Mathematical precision guaranteed")
message(STATUS "🏢 [Enterprise Bean]: Full enterprise-grade monitoring stack ready!")

# [Minimalist Zen]: Simple custom target for quick building
add_custom_target(monitoring_quick
    COMMAND ${CMAKE_COMMAND} --build . --target hsml_performance_monitor_demo --parallel
    COMMENT "🧘 Building monitoring demo with zen-like simplicity..."
)

# [Performance Demon]: Ultra-optimized build target!
add_custom_target(monitoring_maximum_performance
    COMMAND ${CMAKE_COMMAND} --build . --target hsml_performance_monitor_demo --config Release --parallel
    COMMENT "🔥 Building with MAXIMUM PERFORMANCE optimizations!"
)