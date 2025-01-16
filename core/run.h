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
  Run(ctpl::thread_pool& thread_pool, std::mt19937 bitgen, Phylo_tree tree);

  auto bitgen() -> absl::BitGenRef { return bitgen_; }
  auto tree() -> Phylo_tree& { return tree_; }
  auto tree() const -> const Phylo_tree& { return tree_; }
  auto partition() const -> const Partition& { return tree_partition_; }

  auto step() const -> int64_t { return step_; }
  auto set_step(int64_t num_steps) -> void;
  auto local_moves_per_global_move() const -> int { return local_moves_per_global_move_; }
  auto set_local_moves_per_global_move(int local_moves_per_global_move) -> void;
  auto num_parts() const -> int { return num_parts_; }
  auto set_num_parts(int num_parts) -> void { num_parts_ = num_parts; partition_stencils_valid_ = false; }
  auto target_coal_prior_cells() const -> int { return target_coal_prior_cells_; }
  auto set_target_coal_prior_cells(int target_coal_prior_cells) -> void {
    target_coal_prior_cells_ = target_coal_prior_cells; partition_stencils_valid_ = false;
  }
  
  auto do_mcmc_steps(int steps) -> void;

  auto tree_modified() -> void { partition_stencils_valid_ = false; invalidate_derived_quantities(); }
  auto coalescent_prior_t_step() const -> double { return coalescent_prior_.t_step(); }
  auto set_coalescent_prior_t_step(double t_step) -> void;
  auto mu() const -> double { return hky_model_.mu; }
  auto set_mu(double mu) -> void { hky_model_.mu = mpox_mu_ = mu, invalidate_derived_quantities(); }
  auto alpha() const -> double { return alpha_; }
  auto set_alpha(double alpha) -> void { alpha_ = alpha, invalidate_derived_quantities(); }
  auto nu() const -> const std::vector<double>& { return nu_; }
  auto set_nu(const std::vector<double>& nu) -> void { nu_ = nu, invalidate_derived_quantities(); }
  auto hky_kappa() const -> double { return hky_model_.kappa; }
  auto set_hky_kappa(double hky_kappa) -> void { hky_model_.kappa = hky_kappa, invalidate_derived_quantities(); }
  auto hky_pi() const -> Seq_vector<double> { return hky_model_.pi_a; }
  auto set_hky_pi(const Seq_vector<double>& hky_pi) -> void { hky_model_.pi_a = hky_pi, invalidate_derived_quantities(); }
  auto pop_model() const -> const Pop_model& { return *pop_model_; }
  auto set_pop_model(std::shared_ptr<const Pop_model> pop_model) -> void {
    pop_model_ = std::move(pop_model);
    coalescent_prior_.pop_model_changed(pop_model_);
    invalidate_derived_quantities(); }

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

  // Mpox hack.  TODO: We really need something much more generic for supporting different
  // forms of posteriors with different sets of parameters.
  //
  // When we turn on "mpox_hack", the evolution model `evo_` switches to 2 site partitions
  // with the following transition rate matrices:
  //
  //               A    C    G    T                           A    C    G    T
  //
  //           /  -1   1/3  1/3  1/3  \                   /   0    0    0    0  \  .
  //           |                      |                   |                     |  .
  //           |  1/3  -1   1/3  1/3  |                   |   0   -2    0    2  |  .
  //  Q_0 = mu |                      |,  Q_1 = Q_0 + mu* |                     |  .
  //           |  1/3  1/3  -1   1/3  |                   |   2    0   -2    0  |  .
  //           |                      |                   |                     |  .
  //           \  1/3  1/3  1/3  -1   /                   \   0    0    0    0  /  .
  //
  // The sites in partition 1 are those with "APOBEC context", i.e., C or T preceded by
  // a T and G or A followed by an A in the sequence of the first tip.  All other sites
  // are in partition 0.  The idea is to use a simple JC model (given how few mutations
  // we see, no point in trying something more sophisticated), and to approximate
  // the much faster (mu*) APOBEC-mediated mutation of 5'-TC-3' and 5'-GA-3' dimers to TT and AA.
  // To model this without introducing couplings between evolution at different sites,
  // we're assuming that mutations are so rare that we'll never see two mutations
  // in adjacent sites, so it's unlikely for a site to change between having or lacking
  // "APOBEC context".  Merging TC and GA into a single partition
  // keeps things simpler still; you could instead imagine 4 different partitions
  // according to whether a site is (1) preceded or not by a T and/or (2) followed or
  // not by an A.
  //
  // This idea is a slight adaptation of what Andrew Rambaut discussed on a Zoom call with us
  // on 23 May 2023, and which they have implemented in BEAST for their mpox analyses
  // (O'Toole et al, Science, Nov 2023).  We're targetting trees with tips only from after the
  // spillover into humans, so there's no attempt to identify a breakpoint along the tree
  // above which only Q_0 applies.
  //
  // The factors of 2 in Q_1 are chosen to match the convention used by O'Toole et al in
  // reporting APOBEC mutation rates (the effective rate at which APOBEC changes would
  // be observed if 50% of sites with APOBEC context were unmutated and 50% were mutated)
  //
  auto mpox_hack_enabled() const -> bool { return mpox_hack_enabled_; }
  auto set_mpox_hack_enabled(bool mpox_hack_enabled) -> void;
  auto mpox_mu() const -> double { return mpox_mu_; }
  auto set_mpox_mu(double mpox_mu) -> void { mpox_mu_ = mpox_mu, invalidate_derived_quantities(); }
  auto mpox_mu_star() const -> double { return mpox_mu_star_; }
  auto set_mpox_mu_star(double mpox_mu_star) -> void { mpox_mu_star_ = mpox_mu_star, invalidate_derived_quantities(); }
  
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
  auto final_pop_size_move_enabled() const -> bool { return final_pop_size_move_enabled_; }
  auto set_final_pop_size_move_enabled(bool enabled) -> void {
    final_pop_size_move_enabled_ = enabled; }
  auto pop_growth_rate_move_enabled() const -> bool { return pop_growth_rate_move_enabled_; }
  auto set_pop_growth_rate_move_enabled(bool enabled) -> void {
    pop_growth_rate_move_enabled_ = enabled; }

 private:
  ctpl::thread_pool* thread_pool_;
  std::mt19937 bitgen_;
  Phylo_tree tree_;

  int64_t step_;
  int local_moves_per_global_move_;
  int64_t next_global_move_step_;
  
  Partition tree_partition_;
  int num_parts_;
  int target_coal_prior_cells_;
  std::vector<std::vector<Partition_part_info>> partition_stencils_;
  mutable bool partition_stencils_valid_ = false;
  mutable int64_t next_partition_stencil_refresh_step_;
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
  bool final_pop_size_move_enabled_ = true;
  bool pop_growth_rate_move_enabled_ = true;

  // Parameters
  std::shared_ptr<const Pop_model> pop_model_;
  double alpha_;
  Site_vector<double> nu_;
  Hky_model hky_model_;

  bool mpox_hack_enabled_ = false;
  double mpox_mu_ = 0.0;
  double mpox_mu_star_ = 0.0;

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

  auto refresh_partition_stencils() -> void;
  auto repartition() -> void;
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
  auto mpox_hack_moves() -> void;
  auto hky_frequencies_move() -> void;
  auto hky_kappa_move() -> void;
  auto gibbs_sample_all_nus() -> void;
  auto alpha_moves() -> void;
  auto pop_size_move() -> void;
  auto pop_growth_rate_move() -> void;
};

}  // namespace delphy

#endif // DELPHY_RUN_H_
