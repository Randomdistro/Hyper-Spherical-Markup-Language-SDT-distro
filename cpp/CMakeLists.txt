cmake_minimum_required(VERSION 3.20)
project(HSML_SDT 
    VERSION 21.0.0
    DESCRIPTION "Hyper-Spherical Markup Language with Spatial Displacement Theory"
    LANGUAGES CXX)

# C++23 is required for concepts and other modern features
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Compiler-specific optimizations
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native -ffast-math")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra -Wpedantic")
    
    # Enable OpenMP for parallel physics calculations
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(hsml_core PRIVATE OpenMP::OpenMP_CXX)
    endif()
endif()

# Build configuration
option(HSML_BUILD_TESTS "Build HSML-SDT tests" ON)
option(HSML_BUILD_DEMOS "Build HSML-SDT demos" ON)
option(HSML_BUILD_WASM "Build WebAssembly target" OFF)
option(HSML_ENABLE_SIMD "Enable SIMD optimizations" ON)

# Include directories
include_directories(include)

# Core HSML-SDT library
add_library(hsml_core STATIC
    src/core/spherical_math.cpp
    src/core/state_21d.cpp
    src/parser/hsml_parser.cpp
    src/parser/hsml_lexer.cpp
    src/physics/sdt_engine.cpp
    src/physics/sdt_entity.cpp
    src/physics/sdt_field.cpp
    src/viewport/steradian_viewport.cpp
    src/render/spherical_renderer.cpp
    src/render/webgl_backend.cpp
    src/render/webgpu_backend.cpp
    src/render/cpu_backend.cpp
)

target_include_directories(hsml_core PUBLIC include)

# Enable SIMD if requested
if(HSML_ENABLE_SIMD)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(hsml_core PRIVATE -mavx2 -mfma)
    elseif(MSVC)
        target_compile_options(hsml_core PRIVATE /arch:AVX2)
    endif()
    target_compile_definitions(hsml_core PRIVATE HSML_ENABLE_SIMD=1)
endif()

# Runtime executable
add_executable(hsml_runtime
    src/runtime/main.cpp
    src/runtime/cli.cpp
    src/runtime/dev_server.cpp
)

target_link_libraries(hsml_runtime PRIVATE hsml_core)

# Platform-specific libraries
if(WIN32)
    target_link_libraries(hsml_core PRIVATE ws2_32)
    target_link_libraries(hsml_runtime PRIVATE ws2_32)
endif()

if(UNIX AND NOT APPLE)
    target_link_libraries(hsml_core PRIVATE pthread)
    target_link_libraries(hsml_runtime PRIVATE pthread)
endif()

# WebAssembly build configuration
if(HSML_BUILD_WASM)
    if(NOT EMSCRIPTEN)
        message(FATAL_ERROR "HSML_BUILD_WASM requires Emscripten")
    endif()
    
    add_executable(hsml_wasm
        src/wasm/wasm_interface.cpp
        src/wasm/js_bindings.cpp
    )
    
    target_link_libraries(hsml_wasm PRIVATE hsml_core)
    
    set_target_properties(hsml_wasm PROPERTIES
        OUTPUT_NAME "hsml"
        SUFFIX ".js"
        LINK_FLAGS "
            -s WASM=1
            -s USE_WEBGL2=1
            -s USE_WEBGPU=1
            -s ALLOW_MEMORY_GROWTH=1
            -s EXPORTED_FUNCTIONS=['_main','_hsml_init','_hsml_tick','_hsml_render']
            -s EXPORTED_RUNTIME_METHODS=['ccall','cwrap']
            -s MODULARIZE=1
            -s EXPORT_NAME='HSML'
            --pre-js ${CMAKE_SOURCE_DIR}/src/wasm/pre.js
            --post-js ${CMAKE_SOURCE_DIR}/src/wasm/post.js
        "
    )
endif()

# Tests
if(HSML_BUILD_TESTS)
    enable_testing()
    
    # Find or fetch a testing framework (using Catch2)
    find_package(Catch2 3 QUIET)
    if(NOT Catch2_FOUND)
        message(STATUS "Catch2 not found, fetching from GitHub")
        include(FetchContent)
        FetchContent_Declare(
            Catch2
            GIT_REPOSITORY https://github.com/catchorg/Catch2.git
            GIT_TAG v3.4.0
        )
        FetchContent_MakeAvailable(Catch2)
    endif()
    
    add_executable(hsml_tests
        tests/spherical_math_tests.cpp
        tests/sdt_engine_tests.cpp
        tests/parser_tests.cpp
        tests/viewport_tests.cpp
        tests/physics_safety_tests.cpp
    )
    
    target_link_libraries(hsml_tests PRIVATE hsml_core Catch2::Catch2WithMain)
    
    # Register tests with CTest
    add_test(NAME spherical_math COMMAND hsml_tests "[spherical]")
    add_test(NAME sdt_physics COMMAND hsml_tests "[physics]")
    add_test(NAME parser COMMAND hsml_tests "[parser]")
    add_test(NAME viewport COMMAND hsml_tests "[viewport]")
    add_test(NAME safety COMMAND hsml_tests "[safety]")
endif()

# Demo applications
if(HSML_BUILD_DEMOS)
    add_executable(demo_basic
        demos/basic/main.cpp
        demos/basic/basic_scene.cpp
    )
    target_link_libraries(demo_basic PRIVATE hsml_core)
    
    add_executable(demo_physics
        demos/physics/main.cpp
        demos/physics/physics_playground.cpp
    )
    target_link_libraries(demo_physics PRIVATE hsml_core)
    
    add_executable(demo_viewport
        demos/viewport/main.cpp
        demos/viewport/viewport_demo.cpp
    )
    target_link_libraries(demo_viewport PRIVATE hsml_core)
endif()

# Installation
install(TARGETS hsml_core hsml_runtime
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
)

install(DIRECTORY include/hsml DESTINATION include)

# Package configuration
include(CMakePackageConfigHelpers)
configure_package_config_file(
    cmake/HSMLConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/HSMLConfig.cmake
    INSTALL_DESTINATION lib/cmake/HSML
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/HSMLConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/HSMLConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/HSMLConfigVersion.cmake
    DESTINATION lib/cmake/HSML
)

# Build summary
message(STATUS "")
message(STATUS "HSML-SDT Build Configuration:")
message(STATUS "  Version: ${PROJECT_VERSION}")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  C++ standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "  Compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "  Build tests: ${HSML_BUILD_TESTS}")
message(STATUS "  Build demos: ${HSML_BUILD_DEMOS}")
message(STATUS "  Build WASM: ${HSML_BUILD_WASM}")
message(STATUS "  Enable SIMD: ${HSML_ENABLE_SIMD}")
message(STATUS "")