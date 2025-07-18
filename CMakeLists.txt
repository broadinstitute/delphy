# Delphy cmake build definition
# =============================

cmake_minimum_required(VERSION 3.17.5)

# To make a release:
# 1. Commit final changes and note down commit hash X
# 2. Bump VERSION and BUILD_NUMBER below
# 3. Set BUILD_PREV_GIT_HASH to X
# 4. Commit these changes to CMakeLists.txt and push (git push)
# 5. Tag Git commit (git tag ${VERSION} && git push origin ${VERSION})
# 6. Compile & distribute the release artifacts
project(delphy VERSION 1.0.3)
set(BUILD_NUMBER 2039)
set(BUILD_PREV_GIT_HASH 5b2b742)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

# Git hash in version string (adapted from https://jonathanhamberg.com/post/cmake-embedding-git-hash/)
include(${CMAKE_SOURCE_DIR}/CheckGit.cmake)
CheckGitSetup()

# Set up common options
add_library(common INTERFACE)

# Set up compiler and linker flags
target_compile_options(
  common INTERFACE -Wall -Wextra
  -Werror                            # Treat all warnings as errors
)
add_compile_options(-fno-omit-frame-pointer -mno-omit-leaf-frame-pointer)  # Sane stack traces
if(NOT DEFINED EMSCRIPTEN)
  add_compile_options(
    -Wno-error=maybe-uninitialized   # Ignore something in abseil's raw_hash_set.hpp that I cannot fix
  )
  add_compile_options(-ffast-math)
  #if (NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")  # detects macOS, but still some Intel out there
  if(${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64")
    add_compile_options(-msse3)
  endif()
endif ()

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_compile_options(
    -fsanitize=address
    -fno-optimize-sibling-calls
    -fsanitize-address-use-after-scope
    #-fsanitize=undefined
    )
  add_link_options(
    -fsanitize=address
    -fno-optimize-sibling-calls
    -fsanitize-address-use-after-scope
    #-fsanitize=undefined
    )
endif ()

if(EMSCRIPTEN)
  # Enable emscripten pthreads support
  set(CMAKE_CXX_FLAGS "-pthread ${CMAKE_CXX_FLAGS}")

  # And enable proper support for exceptions at compile-time
  set(CMAKE_CXX_FLAGS "-fwasm-exceptions ${CMAKE_CXX_FLAGS}")
  
  # The version of Clang used by Emscripten 3.1.36 on my system doesn't expose key parts of the
  # C++20 library unless -fexperimental-library is specified.
  # Hopefully, the need for this option will disappear as compiler support for C++20 continues maturing.
  target_compile_options(common INTERFACE -fexperimental-library)
endif()

# Find required and optional libraries
find_package(Microsoft.GSL REQUIRED)
find_package(Threads REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED)
set(ABSL_PROPAGATE_CXX_STD ON)
add_subdirectory(third-party/abseil-cpp EXCLUDE_FROM_ALL)

# Optional UI libraries
find_package(SDL2 QUIET)  # QUIET == no warning if not found
find_package(SDL2_ttf QUIET)

# Core library
add_subdirectory(core)

# Tests (only useful in non-WASM builds, but we build them for WASM anyway to
# ensure that no test code trips up Clang [the default compiler on Mac])
cmake_policy(SET CMP0135 NEW)
include(FetchContent)
# To update to latest GoogleTest release:
# 1. Lookup latest commit SHA, e.g., at the top of https://github.com/google/googletest/commits/main
# 2. Adapt the url below to use the commit SHA & bump timestamp in comment
# 3. Recalculate and update the SHA256 of the artifact by running:
#      curl -sL --output - https://github.com/google/googletest/archive/<COMMIT_SHA>.zip | sha256sum -
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/d25e625364a2078c940da2cb1fcf098aa184fd9f.zip # 2022-04-18
  URL_HASH SHA256=45085e7a631f9da804723b8c71a07a32b2ea2da6b699a4bf85039d4b6e1fb7dc
  )
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_subdirectory(tests)

# Set up delphy_ui executable if we have OpenGL and GLUT.
if(SDL2_FOUND AND SDL2_ttf_FOUND)
  add_executable(delphy_ui tools/delphy_ui.cpp)
  target_include_directories(delphy_ui PUBLIC tools)
  target_link_libraries(delphy_ui
    delphy_core
    SDL2::SDL2 SDL2_ttf::SDL2_ttf)
endif()

# delphy cli executables (only in non-WASM builds)
if(NOT DEFINED EMSCRIPTEN)
  add_executable(delphy tools/delphy.cpp)
  target_include_directories(delphy PUBLIC tools)
  target_link_libraries(delphy delphy_core)

  add_executable(delphy_mcc tools/delphy_mcc.cpp)
  target_include_directories(delphy_mcc PUBLIC tools)
  target_link_libraries(delphy_mcc delphy_core)

  add_executable(beast_trees_to_dphy tools/beast_trees_to_dphy.cpp)
  target_include_directories(beast_trees_to_dphy PUBLIC tools)
  target_link_libraries(beast_trees_to_dphy delphy_core)
endif()

# Build the WASM shim only if compiling with Emscripten
if(EMSCRIPTEN)
  add_executable(delphy tools/delphy_wasm.cpp)
  target_link_libraries(delphy delphy_core)
  target_link_options(delphy PRIVATE
    --no-entry -fwasm-exceptions
    $<$<CONFIG:DEBUG>:-gsource-map --source-map-base ./>
    -pthread -sPTHREAD_POOL_SIZE=navigator.hardwareConcurrency
    -sEXPORTED_FUNCTIONS=_malloc,_free
    -sEXPORTED_RUNTIME_METHODS=UTF8ToString,stringToUTF8OnStack

    # Don't mess up 64-bit integers; see https://v8.dev/features/wasm-bigint
    -sWASM_BIGINT
    -sINITIAL_MEMORY=16mb
    -sALLOW_MEMORY_GROWTH=1

    # Avoid grief with decoding UTF-8 strings in SharedArrayBuffers (which is what
    # the WASM linear memory looks like from the JS side).  See here for details:
    # https://github.com/emscripten-core/emscripten/issues/18034
    -sTEXTDECODER=0
    )
endif()
