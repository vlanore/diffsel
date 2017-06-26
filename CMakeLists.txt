cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)
project(diffsel CXX)

# Compilation options
set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra -Wpedantic") # normal mode
set(CXX_FEATURES cxx_auto_type cxx_nullptr cxx_final cxx_override cxx_range_for)
include_directories("src")


########################
#     STATIC LIBS      #
########################

# Phylogeny
set(DIFFSEL_LIB
  src/BranchSitePath.cpp
  src/Chain.cpp
  src/Chrono.cpp
  src/CodonSequenceAlignment.cpp
  src/CodonStateSpace.cpp
  src/CodonSubMatrix.cpp
  src/GTRSubMatrix.cpp
  src/MSCodonSubMatrix.cpp
  src/MultiGeneSingleOmega.cpp
  src/PhyloProcess.cpp
  src/Random.cpp
  src/Sample.cpp
  src/SequenceAlignment.cpp
  src/StateSpace.cpp
  src/SubMatrix.cpp
  src/TaxonSet.cpp
  src/Tree.cpp
  src/Tree.hpp
  src/linalg.cpp
  )
add_library (diffsel_lib STATIC ${DIFFSEL_LIB})
target_compile_features(diffsel_lib PUBLIC ${CXX_FEATURES})


########################
#     EXECUTABLES      #
########################

# Diffsel
add_executable(diffsel "src/DiffSel.cpp")
target_link_libraries(diffsel diffsel_lib)
target_compile_features(diffsel PRIVATE ${CXX_FEATURES})

# ReadDiffSel
add_executable(readdiffsel "src/ReadDiffSel.cpp")
target_link_libraries(readdiffsel diffsel_lib)
target_compile_features(readdiffsel PRIVATE ${CXX_FEATURES})
