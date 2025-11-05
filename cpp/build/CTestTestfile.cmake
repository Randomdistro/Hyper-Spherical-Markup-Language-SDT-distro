# CMake generated Testfile for 
# Source directory: /mnt/c/Users/Jimmi/Downloads/schooltools_complete/Hyper-Spherical-Markup-Language-SDT-distro/cpp
# Build directory: /mnt/c/Users/Jimmi/Downloads/schooltools_complete/Hyper-Spherical-Markup-Language-SDT-distro/cpp/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test([=[spherical_math]=] "/mnt/c/Users/Jimmi/Downloads/schooltools_complete/Hyper-Spherical-Markup-Language-SDT-distro/cpp/build/hsml_tests" "[spherical]")
set_tests_properties([=[spherical_math]=] PROPERTIES  _BACKTRACE_TRIPLES "/mnt/c/Users/Jimmi/Downloads/schooltools_complete/Hyper-Spherical-Markup-Language-SDT-distro/cpp/CMakeLists.txt;141;add_test;/mnt/c/Users/Jimmi/Downloads/schooltools_complete/Hyper-Spherical-Markup-Language-SDT-distro/cpp/CMakeLists.txt;0;")
add_test([=[sdt_physics]=] "/mnt/c/Users/Jimmi/Downloads/schooltools_complete/Hyper-Spherical-Markup-Language-SDT-distro/cpp/build/hsml_tests" "[physics]")
set_tests_properties([=[sdt_physics]=] PROPERTIES  _BACKTRACE_TRIPLES "/mnt/c/Users/Jimmi/Downloads/schooltools_complete/Hyper-Spherical-Markup-Language-SDT-distro/cpp/CMakeLists.txt;142;add_test;/mnt/c/Users/Jimmi/Downloads/schooltools_complete/Hyper-Spherical-Markup-Language-SDT-distro/cpp/CMakeLists.txt;0;")
add_test([=[parser]=] "/mnt/c/Users/Jimmi/Downloads/schooltools_complete/Hyper-Spherical-Markup-Language-SDT-distro/cpp/build/hsml_tests" "[parser]")
set_tests_properties([=[parser]=] PROPERTIES  _BACKTRACE_TRIPLES "/mnt/c/Users/Jimmi/Downloads/schooltools_complete/Hyper-Spherical-Markup-Language-SDT-distro/cpp/CMakeLists.txt;143;add_test;/mnt/c/Users/Jimmi/Downloads/schooltools_complete/Hyper-Spherical-Markup-Language-SDT-distro/cpp/CMakeLists.txt;0;")
add_test([=[viewport]=] "/mnt/c/Users/Jimmi/Downloads/schooltools_complete/Hyper-Spherical-Markup-Language-SDT-distro/cpp/build/hsml_tests" "[viewport]")
set_tests_properties([=[viewport]=] PROPERTIES  _BACKTRACE_TRIPLES "/mnt/c/Users/Jimmi/Downloads/schooltools_complete/Hyper-Spherical-Markup-Language-SDT-distro/cpp/CMakeLists.txt;144;add_test;/mnt/c/Users/Jimmi/Downloads/schooltools_complete/Hyper-Spherical-Markup-Language-SDT-distro/cpp/CMakeLists.txt;0;")
add_test([=[safety]=] "/mnt/c/Users/Jimmi/Downloads/schooltools_complete/Hyper-Spherical-Markup-Language-SDT-distro/cpp/build/hsml_tests" "[safety]")
set_tests_properties([=[safety]=] PROPERTIES  _BACKTRACE_TRIPLES "/mnt/c/Users/Jimmi/Downloads/schooltools_complete/Hyper-Spherical-Markup-Language-SDT-distro/cpp/CMakeLists.txt;145;add_test;/mnt/c/Users/Jimmi/Downloads/schooltools_complete/Hyper-Spherical-Markup-Language-SDT-distro/cpp/CMakeLists.txt;0;")
add_test([=[runtime]=] "/mnt/c/Users/Jimmi/Downloads/schooltools_complete/Hyper-Spherical-Markup-Language-SDT-distro/cpp/build/hsml_tests" "[runtime]")
set_tests_properties([=[runtime]=] PROPERTIES  _BACKTRACE_TRIPLES "/mnt/c/Users/Jimmi/Downloads/schooltools_complete/Hyper-Spherical-Markup-Language-SDT-distro/cpp/CMakeLists.txt;146;add_test;/mnt/c/Users/Jimmi/Downloads/schooltools_complete/Hyper-Spherical-Markup-Language-SDT-distro/cpp/CMakeLists.txt;0;")
subdirs("_deps/catch2-build")
