# HSML Visual Studio GUI - CMake Configuration
# GUI-FIRST METHODOLOGY: Perfect build system for perfect GUI architecture

cmake_minimum_required(VERSION 3.20)

# Find Qt6 - Essential for GUI-first architecture
find_package(Qt6 REQUIRED COMPONENTS 
    Core 
    Widgets 
    OpenGL 
    OpenGLWidgets
    Concurrent
    Network
)

# Enable Qt's meta-object system
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

# GUI source files
set(GUI_SOURCES
    hsml_visual_studio_gui.cpp
    hsml_visual_studio_main.cpp
)

# GUI header files (for MOC processing)
set(GUI_HEADERS
    ../include/hsml/gui/hsml_visual_studio_gui.h
)

# Create GUI library
add_library(hsml_gui STATIC
    ${GUI_SOURCES}
    ${GUI_HEADERS}
)

# Target properties for perfect GUI compilation
target_compile_features(hsml_gui PRIVATE cxx_std_20)
target_compile_definitions(hsml_gui PRIVATE
    QT_DISABLE_DEPRECATED_BEFORE=0x060000
    HSML_GUI_FIRST_ARCHITECTURE=1
)

# Include directories
target_include_directories(hsml_gui PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../../include
    ${Qt6Core_INCLUDE_DIRS}
    ${Qt6Widgets_INCLUDE_DIRS}
    ${Qt6OpenGL_INCLUDE_DIRS}
)

# Link Qt6 libraries
target_link_libraries(hsml_gui
    Qt6::Core
    Qt6::Widgets
    Qt6::OpenGL
    Qt6::OpenGLWidgets
    Qt6::Concurrent
    Qt6::Network
    hsml_core
)

# Platform-specific configurations
if(WIN32)
    target_compile_definitions(hsml_gui PRIVATE
        WIN32_LEAN_AND_MEAN
        NOMINMAX
    )
    # Link Windows-specific libraries
    target_link_libraries(hsml_gui
        opengl32
        gdi32
        user32
        kernel32
    )
elseif(APPLE)
    # macOS-specific frameworks
    target_link_libraries(hsml_gui
        "-framework OpenGL"
        "-framework Cocoa"
    )
    # Enable high DPI support
    target_compile_definitions(hsml_gui PRIVATE
        HSML_MACOS_HIGH_DPI=1
    )
elseif(UNIX)
    # Linux-specific libraries
    find_package(OpenGL REQUIRED)
    target_link_libraries(hsml_gui
        OpenGL::GL
        X11
    )
endif()

# Create main GUI executable
add_executable(hsml_visual_studio
    hsml_visual_studio_main.cpp
)

# Executable properties
target_compile_features(hsml_visual_studio PRIVATE cxx_std_20)
target_include_directories(hsml_visual_studio PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../../include
)

# Link GUI library and dependencies
target_link_libraries(hsml_visual_studio
    hsml_gui
    hsml_core
    Qt6::Core
    Qt6::Widgets
    Qt6::OpenGL
    Qt6::OpenGLWidgets
)

# GUI-specific compiler optimizations
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    # Optimize for GUI performance
    if(MSVC)
        target_compile_options(hstml_gui PRIVATE
            /O2          # Maximum optimization
            /Ob2         # Inline function expansion
            /Ot          # Favor fast code
            /Oy          # Frame pointer omission
            /GL          # Whole program optimization
        )
        target_link_options(hsml_visual_studio PRIVATE
            /LTCG        # Link-time code generation
            /SUBSYSTEM:WINDOWS  # Windows GUI application
        )
    else()
        target_compile_options(hsml_gui PRIVATE
            -O3          # Maximum optimization
            -ffast-math  # Fast math operations
            -funroll-loops  # Loop unrolling
            -fomit-frame-pointer  # Frame pointer omission
        )
    endif()
endif()

# Debug configuration
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_definitions(hsml_gui PRIVATE
        HSML_DEBUG_GUI=1
        QT_QML_DEBUG
    )
endif()

# Resource files (if any)
# qt6_add_resources(hsml_gui "hsml_resources"
#     PREFIX "/"
#     FILES
#         resources/icons/hsml_logo.png
#         resources/themes/dark_theme.qss
#         resources/themes/light_theme.qss
# )

# Installation rules
install(TARGETS hsml_visual_studio
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

# Install Qt6 libraries (for deployment)
if(WIN32)
    # Windows deployment
    find_program(QT_DEPLOYQT_EXECUTABLE windeployqt HINTS ${Qt6_DIR}/../../../bin)
    if(QT_DEPLOYQT_EXECUTABLE)
        add_custom_command(TARGET hsml_visual_studio POST_BUILD
            COMMAND ${QT_DEPLOYQT_EXECUTABLE} $<TARGET_FILE:hsml_visual_studio>
            COMMENT "Deploying Qt libraries"
        )
    endif()
elseif(APPLE)
    # macOS deployment
    find_program(QT_DEPLOYQT_EXECUTABLE macdeployqt HINTS ${Qt6_DIR}/../../../bin)
    if(QT_DEPLOYQT_EXECUTABLE)
        add_custom_command(TARGET hsml_visual_studio POST_BUILD
            COMMAND ${QT_DEPLOYQT_EXECUTABLE} $<TARGET_FILE:hsml_visual_studio>
            COMMENT "Deploying Qt libraries"
        )
    endif()
endif()

# Create desktop entry (Linux)
if(UNIX AND NOT APPLE)
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/resources/hsml_visual_studio.desktop.in
        ${CMAKE_CURRENT_BINARY_DIR}/hsml_visual_studio.desktop
        @ONLY
    )
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/hsml_visual_studio.desktop
        DESTINATION share/applications
    )
endif()

# Documentation generation
find_package(Doxygen QUIET)
if(DOXYGEN_FOUND)
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in
        ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        @ONLY
    )
    add_custom_target(gui_docs
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating GUI API documentation with Doxygen"
        VERBATIM
    )
endif()

# Testing for GUI components
if(BUILD_TESTS)
    find_package(Qt6 COMPONENTS Test QUIET)
    if(Qt6Test_FOUND)
        add_subdirectory(tests)
    endif()
endif()

# Performance profiling build
if(ENABLE_PROFILING)
    target_compile_definitions(hsml_gui PRIVATE
        HSML_ENABLE_PROFILING=1
    )
    if(NOT MSVC)
        target_compile_options(hsml_gui PRIVATE -pg)
        target_link_options(hsml_visual_studio PRIVATE -pg)
    endif()
endif()

# Memory debugging
if(ENABLE_MEMORY_DEBUG)
    target_compile_definitions(hsml_gui PRIVATE
        HSML_MEMORY_DEBUG=1
    )
    if(UNIX AND NOT APPLE)
        target_link_libraries(hsml_visual_studio mcheck)
    endif()
endif()

# Print configuration summary
message(STATUS "HSML Visual Studio GUI Configuration:")
message(STATUS "  Qt6 Version: ${Qt6_VERSION}")
message(STATUS "  Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  GUI-First Architecture: Enabled")
message(STATUS "  Platform: ${CMAKE_SYSTEM_NAME}")
message(STATUS "  Compiler: ${CMAKE_CXX_COMPILER_ID}")

if(Qt6_FOUND)
    message(STATUS "  GUI Status: ✓ Perfect Qt6 integration enabled")
else()
    message(WARNING "  GUI Status: ✗ Qt6 not found - GUI will not be available")
endif()

# Performance warnings
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "  Performance Note: Debug build - GUI performance may be reduced")
    message(STATUS "  Recommendation: Use Release build for optimal GUI experience")
endif()

# GUI-FIRST Architecture validation
if(NOT Qt6Widgets_FOUND)
    message(FATAL_ERROR "GUI-FIRST Architecture requires Qt6 Widgets. Please install Qt6 development packages.")
endif()

if(NOT Qt6OpenGL_FOUND)
    message(FATAL_ERROR "HSML spherical rendering requires Qt6 OpenGL. Please install Qt6 OpenGL development packages.")
endif()

message(STATUS "GUI-FIRST methodology: ✓ All requirements satisfied")
message(STATUS "Ready to build the ultimate HSML development environment!")