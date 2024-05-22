#ifndef DELPHY_SUBRUN_H_
#define DELPHY_SUBRUN_H_

#include "absl/log/check.h"
#include "absl/random/bit_gen_ref.h"

#include "estd.h"
#include "evo_model.h"
#include "phylo_tree.h"
#include "very_scalable_coalescent.h"
#include "sequence_overlay.h"
#include "scratch_space.h"

namespace delphy {

class Subrun {
 public:
  Subrun(absl::BitGenRef bitgen, Phylo_tree tree, bool includes_run_root, Global_evo_model evo);

  auto mcmc_sub_iteration() -> void;

  auto tree() const -> const Phylo_tree& { return tree_; }
  auto includes_run_root() const -> bool { return includes_run_root_; }

  // Parameters
  auto evo() const -> const Global_evo_model& { return evo_; }
  auto set_evo(const Global_evo_model& evo) -> void {
    evo_ = evo, invalidate_derived_quantities(); }
  auto coalescent_prior_part() -> Very_scalable_coalescent_prior_part& { return *coalescent_prior_part_; }
  auto set_coalescent_prior_part(Very_scalable_coalescent_prior_part* coalescent_prior_part) -> void {
    coalescent_prior_part_ = coalescent_prior_part, invalidate_derived_quantities();
  }

  // Move config
  auto only_displacing_inner_nodes() const -> bool { return only_displacing_inner_nodes_; }
  auto set_only_displacing_inner_nodes(bool val) -> void { only_displacing_inner_nodes_ = val; }
  auto topology_moves_enabled() const -> bool { return topology_moves_enabled_; }
  auto set_topology_moves_enabled(bool val) -> void { topology_moves_enabled_ = val; }

  // Derived quantities
  auto log_G() const -> double {
    return validate_derived_quantities(), log_G_; }
  auto state_frequencies_of_ref_sequence_per_partition() const -> const Partition_vector<Seq_vector<int>>& {
    return validate_derived_quantities(), state_frequencies_of_ref_sequence_per_partition_; }
  auto lambda_i() const -> const Node_vector<double>& {
    return validate_derived_quantities(), lambda_i_; }
  auto ref_cum_Q_l() const -> const std::vector<double>& {
    return validate_derived_quantities(), ref_cum_Q_l_;
  }
  auto num_sites_missing_at_every_node() const -> const Node_vector<int>& {
    return validate_derived_quantities(), num_sites_missing_at_every_node_; }
  auto log_augmented_coalescent_prior() const -> double {
    return validate_derived_quantities(), log_augmented_coalescent_prior_; }
  
  auto calc_cur_state_frequencies_of_ref_sequence_per_partition() const -> Partition_vector<Seq_vector<int>>;
  auto calc_cur_ref_cum_Q_l() const -> std::vector<double>;
  auto calc_cur_lambda_i(const std::vector<double>& ref_cum_Q_l) const -> Node_vector<double>;
  auto calc_cur_log_G(
      const Node_vector<double>& lambda_i,
      const Partition_vector<Seq_vector<int>>& state_frequencies_of_ref_sequence_per_partition) const
      -> double;
  auto calc_cur_num_sites_missing_at_every_node() const -> Node_vector<int>;
  auto calc_cur_log_augmented_coalescent_prior() const -> double;
  
  auto invalidate_derived_quantities() -> void { derived_quantities_valid_ = false; }
  auto validate_derived_quantities() const -> void {
    if (not derived_quantities_valid_) {
      recalc_derived_quantities();
      derived_quantities_valid_ = true;
    }
  }

  // Debug support
  auto set_original_sequences(
      Node_vector<Sequence_overlay> original_sequences,
      Node_vector<Interval_set<>> original_missing_sites_all)
      -> void {
    CHECK(estd::is_debug_enabled) << "We only track tip sequences in debug builds";
    original_sequences_ = std::move(original_sequences);
    original_missing_sites_all_ = std::move(original_missing_sites_all);
  }
  auto assert_cur_tip_sequences_compatible_with_original_ones() -> void;
  
 private:
  absl::BitGenRef bitgen_;
  Phylo_tree tree_;
  bool includes_run_root_;
  Scratch_space scratch_;

  // Parameters
  Global_evo_model evo_;
  Very_scalable_coalescent_prior_part* coalescent_prior_part_;

  // Move config
  bool only_displacing_inner_nodes_ = false;
  bool topology_moves_enabled_ = true;

  // Derived quantities
  mutable bool derived_quantities_valid_{false};
  mutable double log_G_{std::numeric_limits<double>::quiet_NaN()};
  mutable Partition_vector<Seq_vector<int>> state_frequencies_of_ref_sequence_per_partition_{};
  mutable Node_vector<double> lambda_i_{};
  mutable std::vector<double> ref_cum_Q_l_{};
  mutable Node_vector<int> num_sites_missing_at_every_node_{};
  mutable double log_augmented_coalescent_prior_{std::numeric_limits<double>::quiet_NaN()};
  
  auto recalc_derived_quantities() const -> void;
  auto check_derived_quantities() const -> void;

  // Debug support
  Node_vector<Sequence_overlay> original_sequences_;  // Only nonempty in debug builds!
  Node_vector<Interval_set<>> original_missing_sites_all_; // Only nonempty in debug builds!

  // Move support
  auto pick_random_node() -> Node_index;
  auto pick_random_inner_node() -> Node_index;

  // Moves
  auto inner_node_displace_move() -> void;
  auto branch_reform_move() -> void;
  auto subtree_slide_move() -> void;
  auto wilson_balding_move() -> void;

  auto spr1_move() -> void;

  auto spr_move_core(
      Node_index X,
      Phylo_tree_loc new_nexus,
      double alpha_new_to_old_over_alpha_old_to_new)
      -> void;
};

}  // namespace delphy

#endif // DELPHY_SUBRUN_H_
