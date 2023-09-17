cmake_minimum_required (VERSION 3.1)
project(SBWT)
include(ExternalProject)

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
configure_file(${PROJECT_SOURCE_DIR}/version.h.in ${PROJECT_BINARY_DIR}/include/version.h @ONLY)

## Require zlib
if(NOT ZLIB)
  if(NOT CMAKE_BUILD_ZLIB)
    find_library(ZLIB NAMES z REQUIRED)
  endif()

  if(CMAKE_BUILD_ZLIB OR ZLIB MATCHES NOTFOUND)
    configure_file(CMakeLists_zlib.txt.in ${CMAKE_CURRENT_BINARY_DIR}/external/zlib-download/CMakeLists.txt)
    execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
      RESULT_VARIABLE result
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/zlib-download )
    if(result)
      message(FATAL_ERROR "CMake step for zlib failed: ${result}")
    endif()
    execute_process(COMMAND ${CMAKE_COMMAND} --build .
      RESULT_VARIABLE result
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/external/zlib-download )
    if(result)
      message(FATAL_ERROR "Build step for zlib failed: ${result}")
    endif()
    set(ZLIB ${CMAKE_CURRENT_BINARY_DIR}/external/zlib/libz.a)
    include_directories(${CMAKE_CURRENT_BINARY_DIR}/external/zlib)
    ExternalProject_Add(zlib
      SOURCE_DIR        "${CMAKE_CURRENT_BINARY_DIR}/external/zlib"
      BUILD_IN_SOURCE   1
      BUILD_COMMAND     $(MAKE) -f ${CMAKE_CURRENT_BINARY_DIR}/external/zlib/Makefile CFLAGS=${CMAKE_C_FLAGS} CXXFLAGS=${CMAKE_CXX_FLAGS}
      CONFIGURE_COMMAND ""
      INSTALL_COMMAND   ""
      TEST_COMMAND      ""
      UPDATE_COMMAND    ""
      )
    include_directories(${CMAKE_CURRENT_BINARY_DIR}/external/zlib)
  endif()
endif()

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

if(CMAKE_BUILD_TYPE MATCHES Release)
  add_compile_options(-O3)
endif()

if(CMAKE_BUILD_TYPE MATCHES Debug)
  add_compile_options(-g -Wall -Wextra -Wno-sign-compare)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

## Require c++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Fix for GCC (see https://gcc.gnu.org/bugzilla/show_bug.cgi?id=90050 )
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(CXX_FILESYSTEM_LIBRARIES "stdc++fs")
endif()

find_package(Threads REQUIRED)

## Add local dependencies as targets
add_subdirectory(${PROJECT_SOURCE_DIR}/sdsl-lite
  ${PROJECT_BINARY_DIR}/external/sdsl-lite/build)

add_subdirectory(${PROJECT_SOURCE_DIR}/KMC
  ${PROJECT_BINARY_DIR}/external/KMC/build)

# SeqIO is header-only so it does not need to be added

# Our sources, excluding the file with the main function
set(SBWT_SOURCES
  src/globals.cpp
  src/suffix_group_optimization.cpp
  src/EM_sort/Block.cpp
  src/EM_sort/EM_sort.cpp
  src/run_kmc.cpp
  src/kmc_construct_helper_classes.cpp)

add_executable(sbwt
  src/CLI/sbwt.cpp
  src/CLI/sbwt_build.cpp
  src/CLI/sbwt_search.cpp
  src/CLI/sbwt_build_from_plain_matrix.cpp
  ${SBWT_SOURCES})
target_include_directories(sbwt PRIVATE
  ${PROJECT_SOURCE_DIR}/include/sbwt/
  ${PROJECT_SOURCE_DIR}/include/sbwt/libwheeler
  ${PROJECT_SOURCE_DIR}/sdsl-lite/include
  ${PROJECT_SOURCE_DIR}/SeqIO/include
  ${PROJECT_BINARY_DIR}/external/sdsl-lite/build/external/libdivsufsort/include
  ${PROJECT_BINARY_DIR}/include
  ${PROJECT_SOURCE_DIR}
  ${PROJECT_SOURCE_DIR}/KMC/include)
target_link_libraries(sbwt PRIVATE
  sdsl
  ${ZLIB}
  Threads::Threads
  ${CXX_FILESYSTEM_LIBRARIES}
  kmc_tools
  kmc_core)

## Build tests if doing a Debug build
if (BUILD_TESTS)
  add_subdirectory(${PROJECT_SOURCE_DIR}/googletest
    ${PROJECT_BINARY_DIR}/external/googletest/build)
  find_package(Threads REQUIRED)
  set(THREADS_PREFER_PTHREAD_FLAG ON)
  find_package(OpenMP REQUIRED)
  find_library(GTEST NAMES gtest PATHS ${PROJECT_SOURCE_DIR}/external/googletest/build/lib/)

  add_executable(sbwt_tests
    tests/test_main.cpp
    src/CLI/sbwt_build.cpp
    src/CLI/sbwt_search.cpp
    ${SBWT_SOURCES})
  target_include_directories(sbwt_tests PRIVATE
    ${PROJECT_SOURCE_DIR}/include/sbwt
    ${PROJECT_SOURCE_DIR}/include/sbwt/libwheeler
    ${PROJECT_SOURCE_DIR}/SeqIO/include
    ${PROJECT_SOURCE_DIR}/sdsl-lite/include
    ${PROJECT_BINARY_DIR}/external/sdsl-lite/build/external/libdivsufsort/include
    ${PROJECT_BINARY_DIR}/include
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/googletest/googletest/include
    ${PROJECT_SOURCE_DIR}/KMC/include)
  target_link_libraries(sbwt_tests PRIVATE
    gtest
    ${ZLIB}
    sdsl
    Threads::Threads
    OpenMP::OpenMP_CXX
    ${CXX_FILESYSTEM_LIBRARIES}
    kmc_tools
    kmc_core)
else()
#  set_target_properties(gtest PROPERTIES EXCLUDE_FROM_ALL 1)
#  set_target_properties(gtest_main PROPERTIES EXCLUDE_FROM_ALL 1)
endif()

add_library(sbwt_static STATIC ${SBWT_SOURCES})
target_include_directories(sbwt_static PRIVATE
  ${PROJECT_SOURCE_DIR}/include/sbwt/
  ${PROJECT_SOURCE_DIR}/include/sbwt/libwheeler
  ${PROJECT_SOURCE_DIR}/sdsl-lite/include
  ${PROJECT_SOURCE_DIR}/SeqIO/include
  ${PROJECT_BINARY_DIR}/external/sdsl-lite/build/external/libdivsufsort/include
  ${PROJECT_BINARY_DIR}/include
  ${PROJECT_SOURCE_DIR}
  ${PROJECT_SOURCE_DIR}/KMC/include)
target_link_libraries(sbwt_static PRIVATE
  sdsl
  ${ZLIB}
  Threads::Threads
  ${CXX_FILESYSTEM_LIBRARIES}
  kmc_tools
  kmc_core)
