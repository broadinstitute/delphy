#ifndef DELPHY_PHYLO_TREE_CALC_H_
#define DELPHY_PHYLO_TREE_CALC_H_

#include <span>

#include "evo_model.h"
#include "phylo_tree.h"
#include "sequence_overlay.h"

namespace delphy {

template<typename T>
auto is_stochastic_vector(const Seq_vector<T>& p, T eps = 1e-6) -> bool {
  auto sum = T{};
  for (auto p_a : p) {
    if (p_a < 0 || p_a > 1) { return false; }
    sum += p_a;
  }
  if (std::abs(sum - 1) > eps) { return false; }
  return true;
}

template<typename T>
auto is_transition_rate_matrix(const Seq_matrix<T>& Q, T eps = 1e-6) -> bool {
  for (auto a : k_all_real_seq_letters) {
    auto expected_qaa = T{};
    for (auto b : k_all_real_seq_letters) {
      if (a != b) {
        expected_qaa -= Q[a][b];
      }
    }
    if (std::abs(Q[a][a] - expected_qaa) > eps) { return false; }
  }
  return true;
}

auto count_mutations(const Phylo_tree& tree) -> int;

auto view_of_sequence_at(const Phylo_tree& tree, Node_index X) -> Sequence_overlay;
auto view_of_sequence_at(const Phylo_tree& tree, Phylo_tree_loc x) -> Sequence_overlay;

auto reconstruct_missing_sites_at(const Phylo_tree& tree, Node_index node) -> Interval_set<>;
auto reconstruct_missing_sites_at(
    const Phylo_tree& tree,
    Node_index node,
    Scratch_space& scratch)
    -> Scratch_interval_set;
auto is_site_missing_at(const Phylo_tree& tree, Node_index node, Site_index l) -> bool;
auto calc_num_sites_missing_at_every_node(const Phylo_tree& tree) -> Node_vector<int>;
auto recalc_num_sites_missing_upstream(
    const Phylo_tree& tree,
    Node_index node,
    Node_index ancestor,
    Node_vector<int>& num_sites_missing_at_every_node)
    -> void;

auto calc_state_frequencies_per_partition_of(
    const Real_sequence& seq,
    const Global_evo_model& evo)
    -> Partition_vector<Seq_vector<int>>;

auto calc_site_state_at(const Phylo_tree& tree, Phylo_tree_loc loc, Site_index l) -> Real_seq_letter;

// Primitive quantities about tree branch lengths and mutation rates
// -----------------------------------------------------------------
// T_0 = total branch length of tree (ignoring missations)
// T^(l)_a = time in tree during which site `l` has state `a` (missations at site `l` prune subtrees)
// nu^(l) = relative mutation rate of site `l` (E(nu^(l)) = 1)
// q^(l)_a = relative escape rate of state a at site `l` (E_a(q^(l)_a) = 1)
// beta(l) = partition of site l
// u^[beta], q^[beta]_a, q^[beta]_{ab} = evo model parameters for a site in partition beta
//
// Derived quantities
// ------------------
// T^(l) = sum_a T^(l)_a = time of tree for site l (might be lower than T_0 owing to missations)
// T = sum_{a,l} T^(l)_a = total time tracked for all sites in tree
// Ttwiddle^(l) = sum_a q^(l)_a T^(l)_a = time of tree for site l, adjusted for state escape rates
// Ttwiddle^\beta_a = sum_{l in beta} nu^(l) T^(l)_a
//                  = time of tree spent by sites of partition beta in state a, adjusted for relative mutation rates
// Ttwiddle = sum_{a,l} q^(l)_a nu^(l) T^(l)_a = total time tracked for all sites in tree,
//                                               adjusted for state escape rate and relative mutation rates
// lambda_i = sum_{l in xi_i} mu^(l) nu^(l) q_{s^(l)_i}
//          = effective rate at which new mutations accumulate just above node i

auto calc_T(const Phylo_tree& tree) -> double;
auto calc_T_l_a(const Phylo_tree& tree) -> std::vector<Seq_vector<double>>;
auto calc_Ttwiddle_l(
    const Phylo_tree& tree,
    const Global_evo_model& evo)
    -> std::vector<double>;
auto calc_Ttwiddle_a(const Phylo_tree& tree, std::span<const double> nu_l) -> Seq_vector<double>;
auto calc_Ttwiddle_beta_a(
    const Phylo_tree& tree,
    const Global_evo_model& evo)
    -> Partition_vector<Seq_vector<double>>;

auto calc_dTtwiddle_a_dt_for_sequence(const Real_sequence& seq, std::span<const double> nu_l) -> Seq_vector<double>;

// cum_Q^(l)_k = sum_{l < k} Q^(l).  Note, 0 <= k <= L, so there are L+1 elements in the vector
auto calc_cum_Q_l_for_sequence(const Real_sequence& seq, const Global_evo_model& evo) -> std::vector<double>;
auto calc_lambda_for_sequence(const Real_sequence& seq, const Global_evo_model& evo) -> double;
auto calc_lambda_at_node(const Phylo_tree& tree, Node_index node, const Global_evo_model& evo) -> double;
auto calc_lambda_at_node(
    const Phylo_tree& tree,
    Node_index node,
    const Global_evo_model& evo,
    const std::vector<double>& ref_cum_Q_l)
    -> double;
auto calc_lambda_i(
    const Phylo_tree& tree,
    const Global_evo_model& evo,
    const std::vector<double>& ref_cum_Q_l)
    -> Node_vector<double>;
auto recalc_lambda_i_upstream(
    const Phylo_tree& tree,
    Node_index node,
    Node_index ancestor,
    const Global_evo_model& evo,
    Node_vector<double>& lambda_i,
    const std::vector<double>& ref_cum_Q_l)
    -> void;

template<typename Alloc>
auto calc_delta_lambda_across_missations(
    const Global_evo_model& evo,
    const Real_sequence& ref_sequence,
    const std::vector<double>& ref_cum_Q_l,
    const Missation_map<Alloc>& missations)
    -> double {
  auto result = 0.0;
  for (const auto& [mi_start, mi_end] : missations.intervals) {
    result -= ref_cum_Q_l[mi_end] - ref_cum_Q_l[mi_start];
  }
  for (const auto& [mi_site, mi_from] : missations.from_states) {
    auto l = mi_site;
    auto ref_from = ref_sequence[l];
    result -= evo.mu_l(l) * evo.nu_l[l] * (evo.q_l_a(l, mi_from) - evo.q_l_a(l, ref_from));
  }
  return result;
}

template<typename Alloc1, typename Alloc2>
auto calc_delta_lambda_across_branch(
    const Global_evo_model& evo,
    const Real_sequence& ref_sequence,
    const std::vector<double>& ref_cum_Q_l,
    const Mutation_list<Alloc1>& mutations,
    const Missation_map<Alloc2>& missations)
    -> double {
  auto result = 0.0;
  for (const auto& m : mutations) {
    auto l = m.site;
    result += evo.mu_l(l) * evo.nu_l[l] * (evo.q_l_a(l, m.to) - evo.q_l_a(l, m.from));
  }
  result += calc_delta_lambda_across_missations(evo, ref_sequence, ref_cum_Q_l, missations);
  return result;
}

auto calc_log_root_prior(
    const Phylo_tree& tree,
    const Global_evo_model& evo)
    -> double;
auto calc_log_root_prior(
    const Phylo_tree& tree,
    const Global_evo_model& evo,
    const Partition_vector<Seq_vector<int>>& state_frequencies_of_ref_sequence_per_partition)
    -> double;
auto calc_log_G_below_root(
    const Phylo_tree& tree,
    const Global_evo_model& evo)
    -> double;
auto calc_log_G_below_root(
    const Phylo_tree& tree,
    const Global_evo_model& evo,
    const Node_vector<double>& lambda_i,
    const Partition_vector<Seq_vector<int>>& state_frequencies_of_ref_sequence_per_partition)
    -> double;

auto calc_branch_log_G(
    const Phylo_tree& tree,
    Branch_index X,
    double lambda_X,
    const Global_evo_model& evo,
    const Partition_vector<Seq_vector<int>>& state_frequencies_of_ref_sequence_per_partition)
    -> double;

template<typename Mutation_range>
auto calc_branch_log_G(
    double t_P,
    double t_X,
    double lambda_X,
    const Global_evo_model& evo,
    const Mutation_range& mutations)
    -> double {
  // Branch genetic prior of P-X branch in the absence of mutations
  auto result = -lambda_X * (t_X - t_P);
  
  // Correct for mutations
  for (const auto& m : mutations | std::views::reverse) {
    auto l = m.site;
    // Remove contribution to log(G) from segment [t_P,m.t) when state is `m.to`
    // and replace it with the analogous contribution when state is `m.from`
    result -= evo.mu_l(l) * evo.nu_l[l] * (evo.q_l_a(l, m.from) - evo.q_l_a(l, m.to)) * (m.t - t_P);
    result += std::log(evo.mu_l(l) * evo.nu_l[l] * evo.q_l_ab(l, m.from, m.to));
  }
  
  return result;
}

auto calc_path_log_G(
    const Phylo_tree& tree,
    Node_index A,
    Node_index B,
    const Global_evo_model& evo,
    const Node_vector<double>& lambda_i,
    const Partition_vector<Seq_vector<int>>& state_frequencies_of_ref_sequence_per_partition)
    -> double;

auto calc_num_muts(const Phylo_tree& tree) -> int;
auto calc_num_muts_ab(const Phylo_tree& tree) -> Seq_matrix<int>;
auto calc_num_muts_beta_ab(const Phylo_tree& tree, const Global_evo_model& evo) -> Partition_vector<Seq_matrix<int>>;
auto calc_num_muts_l(const Phylo_tree& tree) -> Node_vector<int>;
auto calc_num_muts_l_ab(const Phylo_tree& tree) -> Node_vector<Seq_matrix<int>>;

auto calc_max_tip_time(const Phylo_tree& tree) -> double;

void assert_tip_sequences_compatible_with_original_ones(
    const Phylo_tree& tree,
    const Node_vector<Sequence_overlay>& original_sequences,
    const Node_vector<Interval_set<>>& original_missing_sites_all);

}  // namespace delphy

#endif // DELPHY_PHYLO_TREE_CALC_H_
