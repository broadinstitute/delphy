#ifndef DELPHY_SPR_MOVE_H_
#define DELPHY_SPR_MOVE_H_

#include "evo_model.h"
#include "phylo_tree.h"
#include "phylo_tree_calc.h"
#include "site_deltas.h"

namespace delphy {

// Ideal usage of the machinery below:
//
// auto spr = Spr_move{tree, ...};
// auto old_graft = spr.analyze_graft(X);
// spr.peel_graft(old_graft);
// ...
// spr.move(X, new_tree_loc);
// ...
// auto new_graft = spr.propose_new_graft(X, bitgen);
// if (...accept...) {
//   spr.apply_graft(new_graft);
// } else {
//   spr.move(X, old_tree_loc);
//   spr.apply_graft(old_graft);
// }

// An analysis of the paths of a tree upstream of X that would disappear if X were pruned
struct Spr_graft {

  Node_index X;
  Node_index S;
  double t_P;
  
  // For "inner" grafts, where X is not a child of the root node,
  // branch_infos[i] describes a branch from A to B, where A = parent(B) = parent^{i}(X).
  //
  // For "rooty" grafts, where X *is* a child of the root node,
  // * branch_infos[k_branch_info_P_X] covers sites present only along P->X
  // * branch_infos[k_branch_info_P_S] covers sites present only along P->S
  // * branch_infos[k_branch_info_S_P_X] covers sites present along the S->P->X path
  enum {
    k_branch_info_P_X = 0,
    k_branch_info_P_S = 1,
    k_branch_info_S_P_X = 2
  };
  
  struct Branch_info {
    Node_index A;
    Node_index B;
    bool is_open;  // open := start state is unconstrained
    double T_to_X;
    double partial_lambda_at_A;
    double partial_lambda_at_X;
    Scratch_interval_set warm_sites;
    Scratch_interval_set hot_sites;
    Scratch_vector<Mutation> hot_muts_to_X;
    Site_deltas hot_deltas_to_X;
  };
  Scratch_vector<Branch_info> branch_infos;

  double delta_log_G;
  double log_alpha_mut;
};

struct Spr_move {
  Spr_move(
      Phylo_tree& tree,
      double mu_proposal,
      bool can_change_root,
      const Global_evo_model& evo,
      Node_vector<double>& lambda_i,
      const std::vector<double>& ref_cum_Q_l,
      Node_vector<int>& num_sites_missing_at_every_node)
      : tree{&tree},
        mu_proposal{mu_proposal},
        can_change_root{can_change_root},
        evo{&evo},
        lambda_i{&lambda_i},
        ref_cum_Q_l{&ref_cum_Q_l},
        num_sites_missing_at_every_node{&num_sites_missing_at_every_node} {}

  auto analyze_graft(Node_index X) const -> Spr_graft;
  auto move(Node_index X, Node_index SS, double new_t_P) -> void;
  auto propose_new_graft(Node_index X, absl::BitGenRef bitgen) const -> Spr_graft;
  auto peel_graft(const Spr_graft& graft) -> void;
  auto apply_graft(const Spr_graft& graft) -> void;

  auto count_min_mutations(const Spr_graft& graft) -> int;
  auto count_closed_mutations(const Spr_graft& graft) -> int;
  auto summarize_closed_mutations(const Spr_graft& graft) -> Site_deltas;
  
  Phylo_tree* tree;
  double mu_proposal;
  bool can_change_root;
  const Global_evo_model* evo;
  Node_vector<double>* lambda_i;
  const std::vector<double>* ref_cum_Q_l;
  Node_vector<int>* num_sites_missing_at_every_node;

  // Results
  Spr_graft graft_before;
  Spr_graft graft_after;
  double delta_log_G;
  double log_hastings_ratio;

 private:
  auto start_graft_analysis(Node_index X) const -> Spr_graft;
  auto propose_new_graft_mutations(Spr_graft& graft, absl::BitGenRef bitgen) const -> void;
  auto finish_graft_analysis(Spr_graft& graft) const -> void;
  
  auto start_rooty_graft_analysis(Node_index X) const -> Spr_graft;
  auto propose_new_rooty_graft_mutations(Spr_graft& graft, absl::BitGenRef bitgen) const -> void;
  auto finish_rooty_graft_analysis(Spr_graft& graft) const -> void;
  auto peel_rooty_graft(const Spr_graft& graft) -> void;
  auto apply_rooty_graft(const Spr_graft& graft) -> void;
  auto count_rooty_min_mutations(const Spr_graft& graft) -> int;
  auto count_rooty_closed_mutations(const Spr_graft& graft) -> int;
  auto summarize_rooty_closed_mutations(const Spr_graft& graft) -> Site_deltas;
  
  auto start_inner_graft_analysis(Node_index X) const -> Spr_graft;
  auto propose_new_inner_graft_mutations(Spr_graft& graft, absl::BitGenRef bitgen) const -> void;
  auto finish_inner_graft_analysis(Spr_graft& graft) const -> void;
  auto peel_inner_graft(const Spr_graft& graft) -> void;
  auto apply_inner_graft(const Spr_graft& graft) -> void;
  auto count_inner_min_mutations(const Spr_graft& graft) -> int;
  auto count_inner_closed_mutations(const Spr_graft& graft) -> int;
  auto summarize_inner_closed_mutations(const Spr_graft& graft) -> Site_deltas;
};

// Sample a JC69 mutational history trajectory for a sequence of `L` sites, ranging from `-T` to `0`,
// with site mutation rate `mu`.  Each entry `{l, {from, to}}` in the map `deltas` specifies unequal initial
// and final states of site `l`.  All sites not present in `deltas` have initial and final state `A`,
// which should be adjusted by the caller if needed.
// Any random numbers needed should be generated with `bitgen`.
//
// Returns the mutational history as a vector of timed mutations (`{t, {site, from, to}}`), sorted by increasing time.
//
// The algorithm is heavily inspired by Nielsen's rejection sampling method,
// and uniformization as described in Lartillot 2006:
//
// * Nielsen, R, "Mapping Mutations on Phylogenies", Syst. Biol. 51(5):729–739 (2002)
//   https://dx.doi.org/10.1080/10635150290102393
// * Lartillot, N, "Conjugate Gibbs sampling for Bayesian phylogenetic models", J. Comput. Bio. 13(10):1701-1722 (2006)
//   https://dx.doi.org/10.1089/cmb.2006.13.1701
//
auto sample_mutational_history(
    Site_index L,
    double T,
    double mu,
    const Site_deltas& deltas,
    absl::BitGenRef bitgen)
    -> Scratch_vector<Mutation>;

// Sample a JC69 mutational history trajectory for a sequence of `L` sites, ranging from `-T` to `0`,
// with site mutation rate `mu`.  For each site l, the start state is unconstrained and the final state is `A`;
// this should be adjusted by the caller if needed.
// Any random numbers needed should
// be generated with `bitgen`.
//
// Returns the mutational history as a vector of timed mutations (`{t, {site, from, to}}`), sorted by increasing time.
//
auto sample_unconstrained_mutational_history(
    Site_index L,
    double T,
    double mu,
    absl::BitGenRef bitgen)
    -> Scratch_vector<Mutation>;

// Given a mutational history as returned by `sample_mutational_history` or `sample_unconstrained_mutational_history`,
// possibly filtered, looks for mutations involving sites not in `site_deltas` and rotates all their states so that the
// start and end state match the state present at `end_loc` in `tree`.
//
// The rationale behind separating this from sample_mutational_history is to allow filtering of its output
// to eliminate mutations on sites that are actually missing at end_loc.
//
auto adjust_mutational_history(
    Scratch_vector<Mutation>& history,
    const Site_deltas site_deltas,
    const Phylo_tree& tree,
    Phylo_tree_loc end_loc)
    -> void;

}  // namespace delphy

#endif // DELPHY_SPR_MOVE_H_
