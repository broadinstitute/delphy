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
  )

target_link_libraries(tests
  common
  delphy_core
  GTest::gtest
  GTest::gmock)

# Uncomment to leave assembler intermediates around
#set_target_properties(tests PROPERTIES COMPILE_FLAGS "-save-temps")
