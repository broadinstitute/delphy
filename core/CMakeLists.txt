add_library(delphy_core
  # Source files in alphabetical order
  ancestral_tree_prober.cpp
  api.cpp  # See top of api.fbs for how to regenerate api_generated.h
  beasty_input.cpp
  beasty_output.cpp
  cmdline.cpp
  coalescent.cpp
  dates.cpp
  delphy_output.cpp
  evo_hky.cpp
  evo_model.cpp
  incremental_mcc_tree.cpp
  io.cpp
  mcc_tree.cpp
  newick.cpp
  phylo_tree.cpp
  phylo_tree_calc.cpp
  pop_model.cpp
  run.cpp
  scalable_coalescent.cpp
  scratch_space.cpp
  sequence_utils.cpp
  site_deltas.cpp
  site_states_tree_prober.cpp
  spr_move.cpp
  spr_study.cpp
  staircase.cpp
  subrun.cpp
  tree_editing.cpp
  tree_partitioning.cpp
  very_scalable_coalescent.cpp
  )

target_link_libraries(delphy_core
  common
  version
  Threads::Threads
  Eigen3::Eigen
  Boost::headers
  absl::check
  absl::flat_hash_map
  absl::log
  absl::log_initialize
  absl::random_random
  absl::random_bit_gen_ref
  absl::str_format
  absl::time
  Microsoft.GSL::GSL
  )
target_include_directories(delphy_core INTERFACE .)
target_include_directories(delphy_core PUBLIC ../third-party/flatbuffers/include)
target_include_directories(delphy_core PUBLIC ../third-party/ctpl)
target_include_directories(delphy_core PUBLIC ../third-party/cppcoro/include)
target_include_directories(delphy_core PUBLIC ../third-party/cxxopts/include)

# Uncomment to leave assembler intermediates around
#set_target_properties(delphy_core PROPERTIES COMPILE_FLAGS "-save-temps")
