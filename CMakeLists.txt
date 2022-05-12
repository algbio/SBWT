cmake_minimum_required (VERSION 3.1)
project(SBWT)
include(ExternalProject)

## Set maximum k-mer length
if(NOT MAX_KMER_LENGTH)
  set(MAX_KMER_LENGTH 31) ## Default max k-mer length is 31, so that edges are 32-mers and fit into a single machine word
endif()
if(MAX_KMER_LENGTH GREATER 255)
  message(${MAX_KMER_LENGTH})
  message(FATAL_ERROR "Error: MAX_KMER_LENGTH can be at most 255")
endif()
MATH(EXPR K_PLUS_1_MER_LENGTH "${MAX_KMER_LENGTH}+1") # +1 because the construction requires (k+1)-mers
add_definitions(-DKMER_MAX_LENGTH=${K_PLUS_1_MER_LENGTH}) # Define for compiler. Note that the variable name here is KMER_MAX_LENGTH, not MAX_KMER_LENGTH, but that is correct.
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
  message( STATUS "SBWT version: ${SBWT_BUILD_VERSION}")
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
add_executable(kmer-search src/kmer_search.cpp ${NODEBOSS_SOURCES})
add_executable(arrange_bits research/arrange_bits.cpp ${NODEBOSS_SOURCES})
add_executable(suffix_group_optimization research/suffix_group_optimization.cpp ${NODEBOSS_SOURCES})
add_executable(dump_kmers src/dump_kmers.cpp ${NODEBOSS_SOURCES})
add_executable(query_variant_from_file src/query_variant_from_file.cpp ${NODEBOSS_SOURCES})
add_executable(build_plain_matrixboss src/build_plain_matrixboss.cpp ${NODEBOSS_SOURCES})
add_executable(build_variant_from_matrix src/build_variant_from_matrix.cpp ${NODEBOSS_SOURCES})

target_link_libraries(kmer-search sdsl ${CXX_FILESYSTEM_LIBRARIES})
target_link_libraries(arrange_bits sdsl ${CXX_FILESYSTEM_LIBRARIES})
target_link_libraries(suffix_group_optimization sdsl ${CXX_FILESYSTEM_LIBRARIES})
target_link_libraries(dump_kmers sdsl ${CXX_FILESYSTEM_LIBRARIES})
target_link_libraries(query_variant_from_file sdsl ${CXX_FILESYSTEM_LIBRARIES})
target_link_libraries(build_plain_matrixboss sdsl ${CXX_FILESYSTEM_LIBRARIES})
target_link_libraries(build_variant_from_matrix sdsl ${CXX_FILESYSTEM_LIBRARIES})
