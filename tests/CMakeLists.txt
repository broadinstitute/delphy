project(tests)

add_executable(tests
  # Source files in alphabetical order
  ancestral_tree_prober_tests.cpp
  api_tests.cpp
  coalescent_tests.cpp
  dates_tests.cpp
  distributions_tests.cpp
  interval_set_tests.cpp
  io_tests.cpp
  missation_map_tests.cpp
  no_bit_gen_tests.cpp
  phylo_tree_calc_tests.cpp
  phylo_tree_tests.cpp
  pop_model_tests.cpp
  safe_gamma_math_tests.cpp
  scalable_coalescent_tests.cpp
  scratch_space_tests.cpp
  sequence_overlay_tests.cpp
  sequence_tests.cpp
  sequence_utils_tests.cpp
  site_states_tree_prober_tests.cpp
  spr_move_tests.cpp
  spr_study_tests.cpp
  staircase_tests.cpp
  tests_main.cpp
  tree_editing_tests.cpp
  tree_prober_tests.cpp
  tree_tests.cpp
  very_scalable_coalescent_tests.cpp
  )

target_link_libraries(tests
  common
  delphy_core
  GTest::gtest
  GTest::gmock)

if(EMSCRIPTEN)
  target_link_options(tests PRIVATE
    ${WASM_LINKER_FLAGS}
    -sPTHREAD_POOL_SIZE=4
    # Without this, Emscripten's EXIT_RUNTIME defaults to false and the Node.js
    # process always exits with code 0, silently swallowing test failures.
    -sEXIT_RUNTIME=1
  )
  # Create a shell script wrapper so tests can be run as ./tests
  file(GENERATE OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/tests"
    CONTENT "#!/bin/sh\ndir=$(dirname \"$0\")\nexec node \"$dir/$<TARGET_FILE_NAME:tests>\" \"$@\"\n"
    FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                     GROUP_READ GROUP_EXECUTE
                     WORLD_READ WORLD_EXECUTE
  )
endif()

# Uncomment to leave assembler intermediates around
#set_target_properties(tests PROPERTIES COMPILE_FLAGS "-save-temps")
