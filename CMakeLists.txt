cmake_minimum_required (VERSION 3.1)
project(SBWT)
include(ExternalProject)

## Set maximum k-mer length
if(NOT MAX_KMER_LENGTH)
  set(MAX_KMER_LENGTH 32)
endif()
if(MAX_KMER_LENGTH GREATER 255)
  message(${MAX_KMER_LENGTH})
  message(FATAL_ERROR "Error: MAX_KMER_LENGTH can be at most 255")
endif()
add_definitions(-DMAX_KMER_LENGTH=${MAX_KMER_LENGTH}) # Define for compiler.
message(STATUS "Maximum k-mer length is: " ${MAX_KMER_LENGTH})
message(STATUS "To change the maximum k-mer length, use the option -DMAX_KMER_LENGTH=n")
message(STATUS "(n can be up to 255)")

## Determine build type and set flags accordingly
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release) ## Build Release by default
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-sign-compare")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -march=native")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native")

set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Our sources, excluding the file with the main function
set(NODEBOSS_SOURCES src/globals.cpp)

## Require c++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Fix for GCC (see https://gcc.gnu.org/bugzilla/show_bug.cgi?id=90050 )
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(CXX_FILESYSTEM_LIBRARIES "stdc++fs")
endif()

## Get version number from git tags
find_package(Git)
if(GIT_FOUND)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --tags
    WORKING_DIRECTORY "${local_dir}"
    OUTPUT_VARIABLE SBWT_BUILD_VERSION
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  message(STATUS "SBWT version: ${SBWT_BUILD_VERSION}")
else()
  message(STATUS "GIT not found, version numbering disabled")
  set(SBWT_BUILD_VERSION "unknown")
endif()
string(TIMESTAMP SBWT_BUILD_TIMESTAMP)

## Generate a version.h file containing build version and timestamp
configure_file(${CMAKE_SOURCE_DIR}/version.h.in ${CMAKE_BINARY_DIR}/include/version.h @ONLY)


## Add local dependencies as targets
add_subdirectory(${CMAKE_SOURCE_DIR}/sdsl-lite
                 ${CMAKE_BINARY_DIR}/external/sdsl-lite/build)

## Set includes
include_directories(${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/include/libwheeler ${CMAKE_SOURCE_DIR}/sdsl-lite/include ${CMAKE_SOURCE_DIR}/build/external/sdsl-lite/build/external/libdivsufsort/include ${CMAKE_BINARY_DIR}/include)

## Libraries & executables to build

add_executable(sbwt_build src/sbwt_build.cpp ${NODEBOSS_SOURCES})
add_executable(sbwt_kmer_search src/kmer_search.cpp ${NODEBOSS_SOURCES})
add_executable(arrange_bits research/arrange_bits.cpp ${NODEBOSS_SOURCES})
add_executable(suffix_group_optimization research/suffix_group_optimization.cpp ${NODEBOSS_SOURCES})
add_executable(dump_kmers src/dump_kmers.cpp ${NODEBOSS_SOURCES})

target_link_libraries(sbwt_build sdsl ${CXX_FILESYSTEM_LIBRARIES})
target_link_libraries(sbwt_kmer_search sdsl ${CXX_FILESYSTEM_LIBRARIES})
target_link_libraries(arrange_bits sdsl ${CXX_FILESYSTEM_LIBRARIES})
target_link_libraries(suffix_group_optimization sdsl ${CXX_FILESYSTEM_LIBRARIES})
target_link_libraries(dump_kmers sdsl ${CXX_FILESYSTEM_LIBRARIES})


## Build tests if doing a Debug build
if (CMAKE_BUILD_TYPE MATCHES Debug)
  find_package(Threads REQUIRED)
  set(THREADS_PREFER_PTHREAD_FLAG ON)
  find_package(OpenMP REQUIRED)
  add_executable(sbwt_tests tests/test_main.cpp ${NODEBOSS_SOURCES})
  find_library(GTEST NAMES gtest PATHS googletest/build/lib/)
  target_link_libraries(sbwt_tests ${GTEST} sdsl Threads::Threads OpenMP::OpenMP_CXX ${CXX_FILESYSTEM_LIBRARIES})
else()
#  set_target_properties(gtest PROPERTIES EXCLUDE_FROM_ALL 1)
#  set_target_properties(gtest_main PROPERTIES EXCLUDE_FROM_ALL 1)
endif()