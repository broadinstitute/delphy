# Adapted from https://gitlab.com/jhamberg/cmake-examples/-/blob/master/cmake/CheckGit.cmake
# by Jonathan Hamberg (also MIT licensed)

set(CURRENT_LIST_DIR ${CMAKE_CURRENT_LIST_DIR})
if (NOT DEFINED pre_configure_dir)
    set(pre_configure_dir ${CMAKE_CURRENT_LIST_DIR})
endif ()

if (NOT DEFINED post_configure_dir)
    set(post_configure_dir ${CMAKE_BINARY_DIR}/generated)
endif ()

set(pre_configure_file ${pre_configure_dir}/version.cpp.in)
set(post_configure_file ${post_configure_dir}/version.cpp)

function(CheckGitWrite git_hash prev_git_hash git_uncommited_changes)
  file(WRITE ${CMAKE_BINARY_DIR}/git-state.txt
    ${git_hash} "\n"
    ${prev_git_hash} "\n"
    ${git_uncommited_changes})
endfunction()

function(CheckGitRead git_hash prev_git_hash git_uncommited_changes)
    if (EXISTS ${CMAKE_BINARY_DIR}/git-state.txt)
        file(STRINGS ${CMAKE_BINARY_DIR}/git-state.txt CONTENT)
        LIST(GET CONTENT 0 var)
        set(${git_hash} ${var} PARENT_SCOPE)
        LIST(GET CONTENT 1 var)
        set(${prev_git_hash} ${var} PARENT_SCOPE)
        LIST(GET CONTENT 2 var)
        set(${git_uncommitted_changes} ${var} PARENT_SCOPE)
    endif ()
endfunction()

function(CheckGitVersion)
    # Check for uncommitted changes
    execute_process(
        COMMAND git status --porcelain
        WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
        OUTPUT_VARIABLE GIT_STATUS_UNCOMMITTED_CHANGES
        OUTPUT_STRIP_TRAILING_WHITESPACE
        )
    if ("${GIT_STATUS_UNCOMMITTED_CHANGES}" STREQUAL "")
        set(GIT_UNCOMMITED_CHANGES false)
    else ()
        set(GIT_UNCOMMITED_CHANGES true)
    endif ()

    # Get two most recent abbreviated commit hash of the working branch
    execute_process(
        COMMAND git log -2 --format=%h
        WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
        OUTPUT_VARIABLE TWO_GIT_HASHES_STR
        OUTPUT_STRIP_TRAILING_WHITESPACE
        )
    string(REPLACE "\n" ";" TWO_GIT_HASHES_LIST "${TWO_GIT_HASHES_STR}")
    LIST(GET TWO_GIT_HASHES_LIST 0 GIT_HASH)
    LIST(GET TWO_GIT_HASHES_LIST 1 PREV_GIT_HASH)

    CheckGitRead(GIT_HASH_CACHE PREV_GIT_HASH_CACHE GIT_UNCOMMITED_CHANGES_CACHE)
    if (NOT EXISTS ${post_configure_dir})
        file(MAKE_DIRECTORY ${post_configure_dir})
    endif ()

    if (NOT EXISTS ${post_configure_dir}/version.h)
        file(COPY ${pre_configure_dir}/version.h DESTINATION ${post_configure_dir})
    endif()

    if (NOT DEFINED GIT_HASH_CACHE)
        set(GIT_HASH_CACHE "INVALID")
    endif ()
    if (NOT DEFINED PREV_GIT_HASH_CACHE)
        set(PREV_GIT_HASH_CACHE "INVALID")
    endif ()
    if (NOT DEFINED GIT_UNCOMMITED_CHANGES_CACHE)
        set(GIT_UNCOMMITED_CHANGES_CACHE "INVALID")
    endif ()

    # Only update the version.cpp if the hash has changed. This will
    # prevent us from rebuilding the project more than we need to.
    if (NOT ${GIT_HASH} STREQUAL ${GIT_HASH_CACHE}
        OR NOT ${PREV_GIT_HASH} STREQUAL ${PREV_GIT_HASH_CACHE}
        OR NOT ${GIT_UNCOMMITED_CHANGES} STREQUAL "${GIT_UNCOMMITED_CHANGES_CACHE}"
        OR NOT EXISTS ${post_configure_file})
        # Set the GIT_HASH_CACHE and GIT_UNCOMMITED_CHANGES_CACHE variable the next build won't have
        # to regenerate the source file.
        CheckGitWrite(${GIT_HASH} ${PREV_GIT_HASH} ${GIT_UNCOMMITED_CHANGES})

        configure_file(${pre_configure_file} ${post_configure_file} @ONLY)
    endif ()

endfunction()

function(CheckGitSetup)

    add_custom_target(AlwaysCheckGit COMMAND ${CMAKE_COMMAND}
        -DRUN_CHECK_GIT_VERSION=1
        -Dpre_configure_dir=${pre_configure_dir}
        -Dpost_configure_file=${post_configure_dir}
        -DGIT_HASH_CACHE=${GIT_HASH_CACHE}
        -DPREV_GIT_HASH_CACHE=${PREV_GIT_HASH_CACHE}
        -DGIT_UNCOMMITED_CHANGES_CACHE=${GIT_UNCOMMITED_CHANGES_CACHE}
        -DCMAKE_PROJECT_VERSION=${CMAKE_PROJECT_VERSION}
        -DBUILD_NUMBER=${BUILD_NUMBER}
        -DBUILD_PREV_GIT_HASH=${BUILD_PREV_GIT_HASH}
        -P ${CURRENT_LIST_DIR}/CheckGit.cmake
        BYPRODUCTS ${post_configure_file}
        )

    add_library(version ${CMAKE_BINARY_DIR}/generated/version.cpp)
    target_include_directories(version PUBLIC ${CMAKE_BINARY_DIR}/generated)
    add_dependencies(version AlwaysCheckGit)

    CheckGitVersion()
endfunction()

# This is used to run this function from an external cmake process.
if (RUN_CHECK_GIT_VERSION)
    CheckGitVersion()
endif ()
