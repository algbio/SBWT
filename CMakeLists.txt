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
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Our sources, excluding the file with the main function
set(SBWT_SOURCES src/globals.cpp src/suffix_group_optimization.cpp src/EM_sort/Block.cpp src/EM_sort/EM_sort.cpp src/SeqIO.cpp)

## Require c++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Fix for GCC (see https://gcc.gnu.org/bugzilla/show_bug.cgi?id=90050 )
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(CXX_FILESYSTEM_LIBRARIES "stdc++fs")
endif()

find_library(ZLIB z PATHS ${CMAKE_SOURCE_DIR}/KMC/kmc_tools/libs REQUIRED)
find_library(BZIP2 bz2 PATHS ${CMAKE_BINARY_DIR}/external/bzip2 REQUIRED)
find_package(Threads REQUIRED)

## Get version number from git tags
find_package(Git)
if(GIT_FOUND)
  option(GIT_SUBMODULE "Check submodules during build" ON)
  if(GIT_SUBMODULE)
    message(STATUS "Submodule update")
    execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      RESULT_VARIABLE GIT_SUBMOD_RESULT)
    if(NOT GIT_SUBMOD_RESULT EQUAL "0")
      message(FATAL_ERROR "git submodule update --init --recursive failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
    endif()
  endif()
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

add_subdirectory(${CMAKE_SOURCE_DIR}/KMC
  ${CMAKE_BINARY_DIR}/external/KMC/build)

add_executable(sbwt src/CLI/sbwt.cpp src/CLI/sbwt_build.cpp src/CLI/sbwt_search.cpp src/CLI/sbwt_build_from_plain_matrix.cpp ${SBWT_SOURCES})
target_include_directories(sbwt PRIVATE ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/include/libwheeler ${CMAKE_SOURCE_DIR}/sdsl-lite/include ${CMAKE_SOURCE_DIR}/build/external/sdsl-lite/build/external/libdivsufsort/include ${CMAKE_BINARY_DIR}/include ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/KMC/include)
target_link_libraries(sbwt sdsl ${BZIP2} ${ZLIB} Threads::Threads ${CXX_FILESYSTEM_LIBRARIES} kmc_tools kmc_core)
target_compile_options(sbwt PRIVATE -march=native)

## Build tests if doing a Debug build
if (BUILD_TESTS)
  find_package(Threads REQUIRED)
  set(THREADS_PREFER_PTHREAD_FLAG ON)
  find_package(OpenMP REQUIRED)
  add_executable(sbwt_tests tests/test_main.cpp src/CLI/sbwt_build.cpp src/CLI/sbwt_search.cpp ${SBWT_SOURCES})
  target_include_directories(sbwt_tests PRIVATE ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/include/libwheeler ${CMAKE_SOURCE_DIR}/sdsl-lite/include ${CMAKE_SOURCE_DIR}/build/external/sdsl-lite/build/external/libdivsufsort/include ${CMAKE_BINARY_DIR}/include ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/googletest/googletest/include ${CMAKE_SOURCE_DIR}/KMC/include)
  find_library(GTEST NAMES gtest PATHS ${CMAKE_SOURCE_DIR}/googletest/build/lib/)
  target_link_libraries(sbwt_tests ${GTEST} ${BZIP2} ${ZLIB} sdsl Threads::Threads OpenMP::OpenMP_CXX ${CXX_FILESYSTEM_LIBRARIES} kmc_tools kmc_core)
  target_compile_options(sbwt_tests PRIVATE -march=native)
else()
#  set_target_properties(gtest PROPERTIES EXCLUDE_FROM_ALL 1)
#  set_target_properties(gtest_main PROPERTIES EXCLUDE_FROM_ALL 1)
endif()

# add_library(sbwt_static STATIC ${SBWT_SOURCES})
