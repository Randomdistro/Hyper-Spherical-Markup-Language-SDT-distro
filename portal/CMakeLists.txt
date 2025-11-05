cmake_minimum_required(VERSION 3.20)
project(Portal VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Som (hsml_core library)
set(SOM_ROOT "${CMAKE_SOURCE_DIR}/../cpp")
include_directories(${SOM_ROOT}/include)
link_directories(${SOM_ROOT}/build)

# Find OpenMP (required by Som)
find_package(OpenMP REQUIRED)

# Portal basic client example
add_executable(portal_example
    examples/basic_client.cpp
)

target_link_libraries(portal_example PRIVATE hsml_core OpenMP::OpenMP_CXX)

# Platform-specific libs
if(UNIX AND NOT APPLE)
    target_link_libraries(portal_example PRIVATE pthread)
endif()

message(STATUS "")
message(STATUS "Portal Build Configuration:")
message(STATUS "  Som root: ${SOM_ROOT}")
message(STATUS "  C++ standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "")
