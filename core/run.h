#ifndef DELPHY_RUN_H_
#define DELPHY_RUN_H_

#include <random>

#include "ctpl_stl.h"

#include "evo_hky.h"
#include "phylo_tree.h"
#include "scalable_coalescent.h"
#include "tree_partitioning.h"
#include "subrun.h"
#include "very_scalable_coalescent.h"
#include "sequence_overlay.h"

namespace delphy {

class Run {
 public:
  Run(ctpl::thread_pool& thread_pool, absl::BitGenRef bitgen, Phylo_tree tree);

  auto bitgen() -> absl::BitGenRef { return bitgen_; }
  auto tree() -> Phylo_tree& { return tree_; }
  auto tree() const -> const Phylo_tree& { return tree_; }
  auto partition() const -> const Partition& { return tree_partition_; }

  auto step() const -> int64_t { return step_; }
  auto set_step(int64_t num_steps) -> void;
  auto local_moves_per_global_move() const -> int { return local_moves_per_global_move_; }
  auto set_local_moves_per_global_move(int local_moves_per_global_move) -> void;
  auto num_parts() const -> int { return num_parts_; }
  auto set_num_parts(int num_parts) -> void { num_parts_ = num_parts; force_repartition_ = true; }
  
  auto do_mcmc_steps(int steps) -> void;

  auto tree_modified() -> void { invalidate_derived_quantities(); }
  auto coalescent_prior_t_step() const -> double { return coalescent_prior_.t_step(); }
  auto set_coalescent_prior_t_step(double t_step) -> void;
  auto mu() const -> double { return hky_model_.mu; }
  auto set_mu(double mu) -> void { hky_model_.mu = mu, invalidate_derived_quantities(); }
  auto alpha() const -> double { return alpha_; }
  auto set_alpha(double alpha) -> void { alpha_ = alpha, invalidate_derived_quantities(); }
  auto nu() const -> const std::vector<double>& { return nu_; }
  auto set_nu(const std::vector<double>& nu) -> void { nu_ = nu, invalidate_derived_quantities(); }
  auto hky_kappa() const -> double { return hky_model_.kappa; }
  auto set_hky_kappa(double hky_kappa) -> void { hky_model_.kappa = hky_kappa, invalidate_derived_quantities(); }
  auto hky_pi() const -> Seq_vector<double> { return hky_model_.pi_a; }
  auto set_hky_pi(const Seq_vector<double>& hky_pi) -> void { hky_model_.pi_a = hky_pi, invalidate_derived_quantities(); }
  auto pop_model() const -> const Exp_pop_model& { return pop_model_; }
  auto set_pop_model(const Exp_pop_model& pop_model) -> void { pop_model_ = pop_model, invalidate_derived_quantities(); }

  auto evo() const -> const Global_evo_model& {
    return validate_derived_quantities(), evo_; }
  auto pi_a() const -> const Seq_vector<double>& {
    return validate_derived_quantities(), evo_.partition_evo_model[0].pi_a; }
  auto q_ab() const -> const Seq_matrix<double>& {
    return validate_derived_quantities(), evo_.partition_evo_model[0].q_ab; }
  auto log_posterior() const -> double {
    return validate_derived_quantities(), log_G_ + log_coalescent_prior_ + log_other_priors_; }
  auto log_G() const -> double {
    return validate_derived_quantities(), log_G_; }
  auto Ttwiddle_beta_a() const -> const Partition_vector<Seq_vector<double>>& {
    return validate_derived_quantities(), Ttwiddle_beta_a_; }
  auto num_muts() const -> int {
    return validate_derived_quantities(), num_muts_; }
  auto num_muts_ab() const -> const Seq_matrix<int>& {
    return validate_derived_quantities(), num_muts_ab_; }
  auto state_frequencies_of_ref_sequence() const -> const Seq_vector<int>& {
    return validate_derived_quantities(), state_frequencies_of_ref_sequence_; }
  auto log_coalescent_prior() const -> double {
    return validate_derived_quantities(), log_coalescent_prior_; }
  auto log_other_priors() const -> double {
    return validate_derived_quantities(), log_other_priors_; }

  auto calc_cur_log_G() const -> double;
  auto calc_cur_Ttwiddle_beta_a() const -> Partition_vector<Seq_vector<double>>;
  auto calc_cur_log_coalescent_prior() const -> double;
  auto calc_cur_log_augmented_coalescent_prior() const -> double;
  auto calc_cur_num_muts() const -> int;
  auto calc_cur_num_muts_ab() const -> Seq_matrix<int>;
  auto calc_cur_state_frequencies_of_ref_sequence() const -> Seq_vector<int>;
  auto calc_cur_log_other_priors() const -> double;

  auto invalidate_derived_quantities() -> void { derived_quantities_valid_ = false; }
  auto validate_derived_quantities() const -> void {
    if (not derived_quantities_valid_) {
      recalc_derived_quantities();
      derived_quantities_valid_ = true;
    }
  }

  // Selectively suppress some moves for demo purposes
  auto only_displacing_inner_nodes() const -> bool { return only_displacing_inner_nodes_; }
  auto set_only_displacing_inner_nodes(bool only_displacing_inner_nodes) -> void {
    only_displacing_inner_nodes_ = only_displacing_inner_nodes; }
  auto topology_moves_enabled() const -> bool { return topology_moves_enabled_; }
  auto set_topology_moves_enabled(bool topology_moves_enabled) -> void {
    topology_moves_enabled_ = topology_moves_enabled; }
  auto repartitioning_enabled() const -> bool { return repartitioning_enabled_; }
  auto set_repartitioning_enabled(bool repartitioning_enabled) -> void {
    repartitioning_enabled_ = repartitioning_enabled; }
  auto alpha_move_enabled() const -> bool { return alpha_move_enabled_; }
  auto set_alpha_move_enabled(bool enabled) -> void {
    alpha_move_enabled_ = enabled; }
  auto mu_move_enabled() const -> bool { return mu_move_enabled_; }
  auto set_mu_move_enabled(bool enabled) -> void {
    mu_move_enabled_ = enabled; }

  auto repartition() -> void;

 private:
  ctpl::thread_pool* thread_pool_;
  absl::BitGenRef bitgen_;
  Phylo_tree tree_;

  int64_t step_;
  int local_moves_per_global_move_;
  int64_t next_global_move_step_;
  
  Partition tree_partition_;
  int num_parts_;
  bool force_repartition_ = false;
  std::vector<std::mt19937> subbitgens_;
  std::vector<Subrun> subruns_;

  Real_sequence original_ref_sequence_;                    // Only nonempty in debug builds!
  Node_vector<Sequence_overlay> original_sequences_;       // Only nonempty in debug builds!
  Node_vector<Interval_set<>> original_missing_sites_all_; // Only nonempty in debug builds!
  auto save_original_sequences() -> void;
  auto assert_cur_tip_sequences_compatible_with_original_ones() -> void;

  bool only_displacing_inner_nodes_ = false;
  bool topology_moves_enabled_ = true;
  bool repartitioning_enabled_ = true;
  bool alpha_move_enabled_ = false;
  bool mu_move_enabled_ = true;

  // Parameters
  Exp_pop_model pop_model_;
  double alpha_;
  Site_vector<double> nu_;
  Hky_model hky_model_;

  // Derived quantities
  mutable bool derived_quantities_valid_ = false;
  mutable int64_t last_revalidation_step_ = 0;
  mutable Global_evo_model evo_;
  mutable double log_G_ = std::numeric_limits<double>::quiet_NaN();
  mutable Partition_vector<Seq_vector<double>> Ttwiddle_beta_a_ = {};
  mutable int num_muts_ = -1;
  mutable double log_coalescent_prior_ = std::numeric_limits<double>::quiet_NaN();
  mutable double log_other_priors_ = std::numeric_limits<double>::quiet_NaN();
  mutable Seq_matrix<int> num_muts_ab_ = {};
  mutable Seq_vector<int> state_frequencies_of_ref_sequence_ = {};
  mutable Scalable_coalescent_prior coalescent_prior_;
  mutable std::vector<Very_scalable_coalescent_prior_part> coalescent_prior_parts_;

  auto reassemble() -> void;
  auto normalize_root() -> void;
  auto push_tree_to_subruns() -> void;
  auto push_global_params_to_subruns() -> void;
  auto reset_very_scalable_coalescent_parts() -> void;
  
  auto recalc_derived_quantities() const -> void;
  auto check_derived_quantities() const -> void;
  auto check_global_and_local_totals_match() const -> void;

  auto derive_evo() const -> void;

  auto run_local_moves(int count) -> void;
  auto run_global_moves() -> void;
  
  auto mu_move() -> void;
  auto hky_frequencies_move() -> void;
  auto hky_kappa_move() -> void;
  auto gibbs_sample_all_nus() -> void;
  auto alpha_moves() -> void;
  auto pop_size_move() -> void;
  auto pop_growth_rate_move() -> void;
};

}  // namespace delphy

#endif // DELPHY_RUN_H_
