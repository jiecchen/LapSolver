cmake_minimum_required(VERSION 2.8)
project(LapSolver C CXX)

# Necessary packages

find_package (Threads REQUIRED)

# Error messages

set(UNSUPPORTED_CXX "Sorry, the ${CMAKE_CXX_COMPILER_ID} compiler is unsupported.")
set(UNSUPPORTED_CXX "${UNSUPPORTED_CXX} Please use either GCC 4.9+ or Intel 14+.")

# Options

option (STATIC_LINKING "Set this to link statically to appropriate libraries" OFF)

# Detect the compilers, throw an error if it doesn't match the C compiler.

if (NOT "${CMAKE_C_COMPILER_ID}" STREQUAL "${CMAKE_CXX_COMPILER_ID}")
    message (FATAL_ERROR "You must use the same C compiler as C++ compiler.")
endif()

# Detect the compiler, throw an error if it doesn't support Cilk

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    if ("${CMAKE_CXX_COMPILER_VERSION}" VERSION_LESS "14")
        message (FATAL_ERROR "${UNSUPPORTED_CXX}")
    endif ()
    add_definitions("-fast -xHost")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if ("${CMAKE_CXX_COMPILER_VERSION}" VERSION_LESS "4.9")
        message (FATAL_ERROR "${UNSUPPORTED_CXX}")
    endif ()
    add_definitions("-fcilkplus -lcilkrts -O3 -flto -mtune=native -march=native")
else ()
    message (FATAL_ERROR "${UNSUPPORTED_CXX}")
endif ()

# Set up base set of C flags

set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -fopenmp -std=gnu99")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -std=gnu++11")

# Include the source root for easy (java-like) inclusion.

include_directories ("${PROJECT_SOURCE_DIR}/src")

# Build the library.

add_library (kmp SHARED 
    src/algorithms/PartialCholeskyFactorization.cpp
    src/algorithms/PartialCholeskyOrder.cpp
    src/algorithms/PriorityQueue.cpp
    src/algorithms/ShortestPathTree.cpp
    src/algorithms/Stretch.cpp
    src/algorithms/TarjanLCA.cpp
    src/algorithms/TreeDFS.cpp
    src/algorithms/UnionFind.cpp
    src/solvers/ConjugateGradientSolver.cpp
    src/structures/EdgeList.cpp
    src/structures/Graph.cpp
    src/structures/GraphLoader.cpp
    src/structures/TreeChildren.cpp
)

# Set up MKL
# TODO: Make this smart.

set(MKLROOT "/opt/intel/mkl")
include_directories ("${MKLROOT}/include")

if (STATIC_LINKING)
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_L_FLAGS} -static")
    add_definitions("-static")
    find_library(MKL_LP64 NAMES mkl_intel_lp64.a libmkl_intel_lp64.a HINTS ${MKLROOT}/lib ${MKLROOT}/lib/intel64)
    find_library(MKL_CORE NAMES mkl_core.a libmkl_core.a HINTS ${MKLROOT}/lib ${MKLROOT}/lib/intel64)
    find_library(MKL_THREAD NAMES mkl_intel_thread.a libmkl_intel_thread.a HINTS ${MKLROOT}/lib ${MKLROOT}/lib/intel64)
else()
    find_library(MKL_LP64 NAMES mkl_intel_lp64 libmkl_intel_lp64 HINTS ${MKLROOT}/lib ${MKLROOT}/lib/intel64)
    find_library(MKL_CORE NAMES mkl_core libmkl_core HINTS ${MKLROOT}/lib ${MKLROOT}/lib/intel64)
    find_library(MKL_THREAD NAMES mkl_intel_thread libmkl_intel_thread HINTS ${MKLROOT}/lib ${MKLROOT}/lib/intel64)
endif()

set(MKL_LIBS ${MKL_LP64} ${MKL_CORE} ${MKL_THREAD})

# Add all the dependencies

target_link_libraries (kmp ${MKL_LIBS} pthread m)
