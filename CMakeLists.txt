cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)
project(diffsel CXX)

# Compilation options
set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra -Wpedantic") # normal mode
if (CMAKE_BUILD_TYPE MATCHES COVERAGE)
  message(STATUS "Compiling in coverage measurement mode")
  set(CMAKE_CXX_FLAGS "-O0 -fprofile-arcs -ftest-coverage")
endif()
set(CXX_FEATURES cxx_auto_type cxx_nullptr cxx_final cxx_override cxx_range_for)
include_directories("src" "utils")

# static lib
set(DIFFSEL_LIB
  src/BranchSitePath.cpp
  src/Chain.cpp
  src/Chrono.cpp
  src/CodonSequenceAlignment.cpp
  src/CodonStateSpace.cpp
  src/CodonSubMatrix.cpp
  src/GTRSubMatrix.cpp
  src/MSCodonSubMatrix.cpp
  src/PhyloProcess.cpp
  src/Random.cpp
  src/Sample.cpp
  src/SequenceAlignment.cpp
  src/SubMatrix.cpp
  src/Tree.cpp
  src/Tree.hpp
  src/linalg.cpp
  )
add_library (diffsel_lib STATIC ${DIFFSEL_LIB})
target_compile_features(diffsel_lib PUBLIC ${CXX_FEATURES})

# Diffsel
add_executable(diffsel "src/DiffSel.cpp")
target_link_libraries(diffsel diffsel_lib)
target_compile_features(diffsel PRIVATE ${CXX_FEATURES})

# ReadDiffSel
add_executable(readdiffsel "src/ReadDiffSel.cpp")
target_link_libraries(readdiffsel diffsel_lib)
target_compile_features(readdiffsel PRIVATE ${CXX_FEATURES})

# SingleOmega
add_executable(singleomega "src/SingleOmega.cpp")
target_link_libraries(singleomega diffsel_lib)
target_compile_features(singleomega PRIVATE ${CXX_FEATURES})

# MultiGeneSingleOmega
add_executable(multigenesingleomega "src/MultiGeneSingleOmega.cpp")
target_link_libraries(multigenesingleomega diffsel_lib)
target_compile_features(multigenesingleomega PRIVATE ${CXX_FEATURES})

# Tracecomp
add_executable(tracecomp "utils/TraceComp.cpp" "utils/Correlation.cpp")
target_link_libraries(tracecomp diffsel_lib)
target_compile_features(tracecomp PRIVATE ${CXX_FEATURES})

# Tests
add_executable(tests "src/test.cpp")
target_link_libraries(tests diffsel_lib)
target_compile_features(tests PRIVATE ${CXX_FEATURES})
