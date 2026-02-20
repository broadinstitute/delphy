#include "spr_move.h"

#include "distributions.h"
#include "phylo_tree_calc.h"
#include "tree_editing.h"

namespace delphy {

auto Spr_move::analyze_graft(Node_index X) const -> Spr_graft {
  auto graft = start_graft_analysis(X);
  finish_graft_analysis(graft);
  return graft;
}

auto Spr_move::propose_new_graft(Node_index X, absl::BitGenRef bitgen) const -> Spr_graft {
  auto graft = start_graft_analysis(X);
  propose_new_graft_mutations(graft, bitgen);
  finish_graft_analysis(graft);
  return graft;
}

auto Spr_move::start_graft_analysis(Node_index X) const -> Spr_graft {
  if (tree->at(X).parent == tree->root) {
    return start_rooty_graft_analysis(X);
  } else {
    return start_inner_graft_analysis(X);
  }
}

auto Spr_move::propose_new_graft_mutations(Spr_graft& graft, absl::BitGenRef bitgen) const -> void {
  if (tree->at(graft.X).parent == tree->root) {
    return propose_new_rooty_graft_mutations(graft, bitgen);
  } else {
    return propose_new_inner_graft_mutations(graft, bitgen);
  }
}

auto Spr_move::finish_graft_analysis(Spr_graft& graft) const -> void {
  if (tree->at(graft.X).parent == tree->root) {
    return finish_rooty_graft_analysis(graft);
  } else {
    return finish_inner_graft_analysis(graft);
  }
}

auto Spr_move::peel_graft(const Spr_graft& graft) -> void {
  CHECK_NE(graft.X, tree->root);
  if (tree->at(graft.X).parent == tree->root) {
    return peel_rooty_graft(graft);
  } else {
    return peel_inner_graft(graft);
  }
}

auto Spr_move::apply_graft(const Spr_graft& graft) -> void {
  CHECK_NE(graft.X, tree->root);
  if (tree->at(graft.X).parent == tree->root) {
    return apply_rooty_graft(graft);
  } else {
    return apply_inner_graft(graft);
  }
}

auto Spr_move::count_min_mutations(const Spr_graft& graft) -> int {
  CHECK_NE(graft.X, tree->root);
  if (tree->at(graft.X).parent == tree->root) {
    return count_rooty_min_mutations(graft);
  } else {
    return count_inner_min_mutations(graft);
  }
}

auto Spr_move::count_closed_mutations(const Spr_graft& graft) -> int {
  CHECK_NE(graft.X, tree->root);
  if (tree->at(graft.X).parent == tree->root) {
    return count_rooty_closed_mutations(graft);
  } else {
    return count_inner_closed_mutations(graft);
  }
}

auto Spr_move::summarize_closed_mutations(const Spr_graft& graft) -> Site_deltas {
  CHECK_NE(graft.X, tree->root);
  if (tree->at(graft.X).parent == tree->root) {
    return summarize_rooty_closed_mutations(graft);
  } else {
    return summarize_inner_closed_mutations(graft);
  }
}

auto Spr_move::start_rooty_graft_analysis(Node_index X) const -> Spr_graft {
  //
  // A rooty graft looks like this:
  //
  //                  +miss_X---muts_X--- X
  //                  |
  // miss_P-muts_P----+ P   <-- root node
  //                  |
  //                  +miss_S---muts_S--- S
  //
  // Sites along the path S->P->X path fall in one of three categories:
  // * Sites present only along P->X ( == miss_S.intervals)
  // * Sites present only along P->S ( == miss_X.intervals)
  // * Sites present all along ( == ALL - miss_P.intervals - miss_S.intervals - miss_X.intervals)
  //
  // We use branch_infos 0 (k_branch_info_P_X), 1 (k_branch_info_P_S) and 2 (k_branch_info_S_P_X)
  // to describe these paths, respectively.

  CHECK_NE(X, tree->root);
  auto t_X = tree->at(X).t;
  auto P = tree->at(X).parent;
  CHECK_EQ(P, tree->root);
  auto t_P = tree->at(P).t;
  auto S = tree->at(P).sibling_of(X);
  auto t_S = tree->at(S).t;

  CHECK(can_change_root);

  const auto& miss_P = tree->at(P).missations;
  const auto& miss_X = tree->at(X).missations;
  const auto& miss_S = tree->at(S).missations;

  const auto& muts_X = tree->at(X).mutations;
  const auto& muts_S = tree->at(S).mutations;
  
  auto graft = Spr_graft{};
  graft.X = X;
  graft.S = S;
  graft.t_P = t_P;

  graft.branch_infos.resize(3);
  
  auto& P_X_path = graft.branch_infos[Spr_graft::k_branch_info_P_X];
  P_X_path.A = P;
  P_X_path.B = X;
  P_X_path.is_open = true;
  P_X_path.T_to_X = t_X - t_P;
  P_X_path.partial_lambda_at_A =
      -1 * calc_delta_lambda_across_missations(*evo, tree->ref_sequence, *ref_cum_Q_l, miss_S);
  P_X_path.warm_sites = miss_S.intervals;
  P_X_path.hot_sites = P_X_path.warm_sites;
  
  P_X_path.partial_lambda_at_X = P_X_path.partial_lambda_at_A;
  for (const auto& m : muts_X) {
    if (P_X_path.hot_sites.contains(m.site)) {
      P_X_path.hot_muts_to_X.push_back(m);
      //push_back_site_deltas(m, P_X_path.hot_deltas_to_X);  // No hot_deltas since P->X is open
      P_X_path.partial_lambda_at_X += evo->mu_l(m.site) * evo->nu_l.at(m.site) *
          (-evo->q_l_a(m.site, m.from) + evo->q_l_a(m.site, m.to));
    }
  }

  
  auto& P_S_path = graft.branch_infos[Spr_graft::k_branch_info_P_S];
  P_S_path.A = P;
  P_S_path.B = S;
  P_S_path.is_open = true;
  P_S_path.T_to_X = t_S - t_P;
  P_S_path.partial_lambda_at_A =
      -1 * calc_delta_lambda_across_missations(*evo, tree->ref_sequence, *ref_cum_Q_l, miss_X);
  P_S_path.warm_sites = miss_X.intervals;
  P_S_path.hot_sites = P_S_path.warm_sites;
  
  P_S_path.partial_lambda_at_X = P_S_path.partial_lambda_at_A;
  for (const auto& m : muts_S) {
    if (P_S_path.hot_sites.contains(m.site)) {
      P_S_path.hot_muts_to_X.push_back(m);
      //push_back_site_deltas(m, P_S_path.hot_deltas_to_X);  // No hot_deltas since P->S is open
      P_S_path.partial_lambda_at_X += evo->mu_l(m.site) * evo->nu_l.at(m.site) *
          (-evo->q_l_a(m.site, m.from) + evo->q_l_a(m.site, m.to));
    }
  }

  
  auto& S_P_X_path = graft.branch_infos[Spr_graft::k_branch_info_S_P_X];
  S_P_X_path.A = S;
  S_P_X_path.B = P;
  S_P_X_path.is_open = false;
  S_P_X_path.T_to_X = (t_S - t_P) + (t_X - t_P);
  S_P_X_path.partial_lambda_at_X = lambda_i->at(X) - P_X_path.partial_lambda_at_X;
  S_P_X_path.partial_lambda_at_A = lambda_i->at(S) - P_S_path.partial_lambda_at_X;
  
  // Use hot_sites and warm_sites as scratch space during the following claculation
  S_P_X_path.hot_sites.insert({0, tree->num_sites()});
  subtract_interval_sets(S_P_X_path.warm_sites, S_P_X_path.hot_sites, miss_P.intervals);
  subtract_interval_sets(S_P_X_path.hot_sites, S_P_X_path.warm_sites, miss_X.intervals);
  subtract_interval_sets(S_P_X_path.warm_sites, S_P_X_path.hot_sites, miss_S.intervals);
  S_P_X_path.hot_sites = S_P_X_path.warm_sites;
  
  for (const auto& m : muts_S | std::views::reverse) {
    if (S_P_X_path.hot_sites.contains(m.site)) {
      auto reversed_m = Mutation{m.to, m.site, m.from, t_P - (m.t - t_P)};
      S_P_X_path.hot_muts_to_X.push_back(reversed_m);
      push_back_site_deltas(reversed_m, S_P_X_path.hot_deltas_to_X);
    }
  }
  for (const auto& m : muts_X) {
    if (S_P_X_path.hot_sites.contains(m.site)) {
      S_P_X_path.hot_muts_to_X.push_back(m);
      push_back_site_deltas(m, S_P_X_path.hot_deltas_to_X);
    }
  }

  return graft;
}

auto Spr_move::propose_new_rooty_graft_mutations(Spr_graft& graft, absl::BitGenRef bitgen) const -> void {
  auto X = graft.X;
  CHECK_NE(X, tree->root);
  auto P = tree->at(X).parent;
  auto S = tree->at(P).sibling_of(X);
  
  for (auto branch_info_idx = 0; branch_info_idx != std::ssize(graft.branch_infos); ++branch_info_idx) {
    auto& branch_info = graft.branch_infos[branch_info_idx];
    
    CHECK(not branch_info.is_open || branch_info.hot_deltas_to_X.empty());
    
    if (not branch_info.hot_sites.empty()) {
      auto new_mutations =
          branch_info.is_open
          ? sample_unconstrained_mutational_history(
              tree->num_sites(), branch_info.T_to_X, mu_proposal, bitgen)
          : sample_mutational_history(
                tree->num_sites(), branch_info.T_to_X, mu_proposal,
                branch_info.hot_deltas_to_X, bitgen);
      if (not new_mutations.empty()) {
        std::erase_if(new_mutations, [&](const auto& m) { return not branch_info.hot_sites.contains(m.site); });

        auto path_end = (branch_info_idx == Spr_graft::k_branch_info_P_S ? tree->node_loc(S) : tree->node_loc(X));
        adjust_mutational_history(new_mutations, branch_info.hot_deltas_to_X, *tree, path_end);
      }
      branch_info.hot_muts_to_X = std::move(new_mutations);

      // For open paths, partial_lambda_at_A may have changed
      if (branch_info.is_open) {
        branch_info.partial_lambda_at_A = branch_info.partial_lambda_at_X;
        for (const auto& m : branch_info.hot_muts_to_X | std::views::reverse) {
          branch_info.partial_lambda_at_A += evo->mu_l(m.site) * evo->nu_l.at(m.site) *
              (+evo->q_l_a(m.site, m.from) - evo->q_l_a(m.site, m.to));
        }
      }
    }
  }
}

auto Spr_move::finish_rooty_graft_analysis(Spr_graft& graft) const -> void {
  auto X = graft.X;
  auto t_X = tree->at(X).t;
  auto P = tree->at(X).parent;
  auto t_P = tree->at(P).t;
  auto S = tree->at(P).sibling_of(X);
  auto t_S = tree->at(S).t;

  const auto& P_X_path = graft.branch_infos[Spr_graft::k_branch_info_P_X];
  const auto& P_S_path = graft.branch_infos[Spr_graft::k_branch_info_P_S];
  const auto& S_P_X_path = graft.branch_infos[Spr_graft::k_branch_info_S_P_X];
  
  // Contribution of hot paths to log_G
  graft.delta_log_G = 0.0;
  graft.delta_log_G += calc_branch_log_G(t_P, t_X, P_X_path.partial_lambda_at_X, *evo, P_X_path.hot_muts_to_X);
  graft.delta_log_G += calc_branch_log_G(t_P, t_S, P_S_path.partial_lambda_at_X, *evo, P_S_path.hot_muts_to_X);
  
  // S-P-X path is a bit trickier since the substitution model need not be time-reversible
  auto S_P_X_muts_P_S_muts = Scratch_vector<Mutation>{};
  for (const auto& m : S_P_X_path.hot_muts_to_X | std::views::reverse) {
    if (m.t < t_P) {
      S_P_X_muts_P_S_muts.push_back(Mutation{m.to, m.site, m.from, t_P + (t_P - m.t)});
    }
  }
  auto S_P_X_muts_P_X_muts = Scratch_vector<Mutation>{};
  for (const auto& m : S_P_X_path.hot_muts_to_X) {
    if (m.t >= t_P) {
      S_P_X_muts_P_X_muts.push_back(m);
    }
  }
  graft.delta_log_G += calc_branch_log_G(t_P, t_X, S_P_X_path.partial_lambda_at_X, *evo, S_P_X_muts_P_X_muts);
  graft.delta_log_G += calc_branch_log_G(t_P, t_S, S_P_X_path.partial_lambda_at_A, *evo, S_P_X_muts_P_S_muts);
  
  // Final complication: root sequence changes
  // For P->X and P->S paths, imagine slipping in the mutations right-to-left from above the root
  // into to their final locations
  for (const auto& m : P_X_path.hot_muts_to_X) {
    graft.delta_log_G += std::log(evo->pi_l_a(m.site, m.from) / evo->pi_l_a(m.site, m.to));
  }
  for (const auto& m : P_S_path.hot_muts_to_X) {
    graft.delta_log_G += std::log(evo->pi_l_a(m.site, m.from) / evo->pi_l_a(m.site, m.to));
  }
  // For the S->P->X path, imagine the mutations all start on the P-X branch and are slid left-to-right
  // into their final locations; those that pass through the root at the ones that end up in the P->S branch.
  for (const auto& m : S_P_X_muts_P_S_muts) {
    graft.delta_log_G += std::log(evo->pi_l_a(m.site, m.from) / evo->pi_l_a(m.site, m.to));
  }

  // Probability of choosing this particular mutational history (excluding factors common to all histories)
  graft.log_alpha_mut = 0.0;
  for (const auto& branch_info : graft.branch_infos) {
    auto L = branch_info.hot_sites.num_sites();
    auto T = branch_info.T_to_X;
    auto M = std::ssize(branch_info.hot_muts_to_X);

    graft.log_alpha_mut +=
        -mu_proposal * L * T
        + M * std::log(mu_proposal / 3);

    if (not branch_info.is_open) {
      auto d = std::ssize(branch_info.hot_deltas_to_X);
      
      auto P_AC_JC = -0.25 * std::expm1(-4./3.*mu_proposal * T);
      auto log_P_AC_JC = std::log(P_AC_JC);
      //auto P_AA_JC = 1.0 - 3 * P_AC_JC;
      auto log_P_AA_JC = std::log1p(-3 * P_AC_JC);
      
      graft.log_alpha_mut -= (L-d)*log_P_AA_JC + d*log_P_AC_JC;
    }
  }
}

auto Spr_move::peel_rooty_graft(const Spr_graft& graft) -> void {
  // The tricky bit here is to fix up the missations and the root sequence
  // TODO: Now that the basic idea is settled and works, there are more direct ways of implementing the below
  CHECK(can_change_root);
  
  auto X = graft.X;
  CHECK_NE(X, tree->root);
  auto t_X = tree->at(X).t;
  auto P = tree->at(X).parent;
  CHECK_EQ(P, tree->root);
  auto t_P = tree->at(P).t;
  auto S = tree->at(P).sibling_of(X);

  CHECK(can_change_root);

  auto& muts_P = tree->at(P).mutations;
  auto& muts_X = tree->at(X).mutations;
  auto& muts_S = tree->at(S).mutations;

  auto& miss_X = tree->at(X).missations;
  auto& miss_S = tree->at(S).missations;

  const auto& P_X_path = graft.branch_infos[Spr_graft::k_branch_info_P_X];
  const auto& P_S_path = graft.branch_infos[Spr_graft::k_branch_info_P_S];
  const auto& S_P_X_path = graft.branch_infos[Spr_graft::k_branch_info_S_P_X];

  // Remove all existing mutations (which may change the root sequence)
  auto ref_to_root_deltas = Site_deltas{};
  for (const auto& m : muts_P) {
    push_back_site_deltas(m, ref_to_root_deltas);
  }
  
  // a) remove all P-X mutations by effectively sliding them left-to-right upstream of the root
  for (const auto& m : muts_X) {
    if (P_X_path.hot_sites.contains(m.site)) {
      if (estd::is_debug_enabled) {
        if (ref_to_root_deltas.contains(m.site)) {
          DCHECK_EQ(ref_to_root_deltas.at(m.site).to, m.from);
        } else {
          DCHECK_EQ(tree->ref_sequence[m.site], m.from);
        }
      }
      push_back_site_deltas(m, ref_to_root_deltas);
      
      DCHECK(not miss_X.contains(m.site));
      DCHECK(miss_S.contains(m.site));
      DCHECK_EQ(miss_S.get_from_state(m.site, tree->ref_sequence), m.from);
      miss_S.set_from_state(m.site, m.to, tree->ref_sequence);
    }
  }

  // b) remove all P-S mutations by effectively sliding them left-to-right upstream of the root
  for (const auto& m : muts_S) {
    if (P_S_path.hot_sites.contains(m.site)) {
      if (estd::is_debug_enabled) {
        if (ref_to_root_deltas.contains(m.site)) {
          DCHECK_EQ(ref_to_root_deltas.at(m.site).to, m.from);
        } else {
          DCHECK_EQ(tree->ref_sequence[m.site], m.from);
        }
      }
      push_back_site_deltas(m, ref_to_root_deltas);
      
      DCHECK(miss_X.contains(m.site));
      DCHECK(not miss_S.contains(m.site));
      DCHECK_EQ(miss_X.get_from_state(m.site, tree->ref_sequence), m.from);
      miss_X.set_from_state(m.site, m.to, tree->ref_sequence);
    }
  }
  
  // c) remove all S-P-X mutations by effectively sliding them right-to-left onto the P-X branch.
  for (const auto& m : muts_S) {
    if (S_P_X_path.hot_sites.contains(m.site)) {
      if (estd::is_debug_enabled) {
        if (ref_to_root_deltas.contains(m.site)) {
          DCHECK_EQ(ref_to_root_deltas.at(m.site).to, m.from);
        } else {
          DCHECK_EQ(tree->ref_sequence[m.site], m.from);
        }
      }
      push_back_site_deltas(m, ref_to_root_deltas);
      
      DCHECK(not miss_X.contains(m.site));
      DCHECK(not miss_S.contains(m.site));
    }
  }
  muts_X.clear();
  muts_S.clear();
  muts_P.clear();

  // Temporarily place fake mutations on the P-X branch to restore the tree's consistency
  auto t_mut_X = 0.5*(t_P + t_X);
  for (const auto& [l, delta] : S_P_X_path.hot_deltas_to_X) {
    if (estd::is_debug_enabled) {
      if (ref_to_root_deltas.contains(l)) {
        DCHECK_EQ(ref_to_root_deltas.at(l).to, delta.from);
      } else {
        DCHECK_EQ(tree->ref_sequence[l], delta.from);
      }
    }
    muts_X.push_back(Mutation{delta.from, l, delta.to, t_mut_X});
  }
  
  // Fix up mutations above root
  for (const auto& [l, delta] : ref_to_root_deltas) {
    CHECK_EQ(tree->ref_sequence[l], delta.from);
    muts_P.push_back(Mutation{delta.from, l, delta.to, -std::numeric_limits<double>::max()});
  }

  // Recalculate lambda_i at P (might differ from lambda_i at S owing to missations above S)
  lambda_i->at(P) = calc_lambda_at_node(*tree, P, *evo, *ref_cum_Q_l);

  assert_phylo_tree_integrity(*tree);
}

auto Spr_move::apply_rooty_graft(const Spr_graft& graft) -> void {
  // The tricky bit here is to fix up the missations and the root sequence
  CHECK(can_change_root);
  
  auto X = graft.X;
  CHECK_NE(X, tree->root);
  auto t_X = tree->at(X).t;
  auto P = tree->at(X).parent;
  CHECK_EQ(P, tree->root);
  auto t_P = tree->at(P).t;
  auto S = tree->at(P).sibling_of(X);
  auto t_S = tree->at(S).t;

  CHECK(can_change_root);

  auto& muts_P = tree->at(P).mutations;
  auto& muts_X = tree->at(X).mutations;
  auto& muts_S = tree->at(S).mutations;

  auto& miss_X = tree->at(X).missations;
  auto& miss_S = tree->at(S).missations;

  const auto& P_X_path = graft.branch_infos[Spr_graft::k_branch_info_P_X];
  const auto& P_S_path = graft.branch_infos[Spr_graft::k_branch_info_P_S];
  const auto& S_P_X_path = graft.branch_infos[Spr_graft::k_branch_info_S_P_X];

  // We assume the graft was either peeled before, or we're here after a `move` operation,
  // which never inserts spurious mutations outside the P-X branch.
  CHECK(muts_S.empty());
  
  // We're about to put in place all the recorded mutations on the P-X branch, so nuke what's already there
  muts_X.clear();
  
  // Record deltas to the root sequence, since they may change
  auto ref_to_root_deltas = Site_deltas{};
  for (const auto& m : muts_P) {
    push_back_site_deltas(m, ref_to_root_deltas);
  }
  muts_P.clear();  // Refilled appropriately below
  
  // Now slide all the new mutations in place in the opposite order (we sort by time later)
  for (const auto& m : P_X_path.hot_muts_to_X | std::views::reverse) {
    muts_X.push_back(m);

    if (estd::is_debug_enabled) {
      if (ref_to_root_deltas.contains(m.site)) {
        DCHECK_EQ(ref_to_root_deltas.at(m.site).to, m.to);
      } else {
        DCHECK_EQ(tree->ref_sequence[m.site], m.to);
      }
    }
    push_back_site_deltas({m.site, m.to, m.from}, ref_to_root_deltas);
    
    DCHECK(not miss_X.contains(m.site));
    DCHECK(miss_S.contains(m.site));
    DCHECK_EQ(miss_S.get_from_state(m.site, tree->ref_sequence), m.to);
    miss_S.set_from_state(m.site, m.from, tree->ref_sequence);
  }
  
  for (const auto& m : P_S_path.hot_muts_to_X | std::views::reverse) {
    muts_S.push_back(m);

    if (estd::is_debug_enabled) {
      if (ref_to_root_deltas.contains(m.site)) {
        DCHECK_EQ(ref_to_root_deltas.at(m.site).to, m.to);
      } else {
        DCHECK_EQ(tree->ref_sequence[m.site], m.to);
      }
    }
    push_back_site_deltas({m.site, m.to, m.from}, ref_to_root_deltas);
    
    DCHECK(miss_X.contains(m.site));
    DCHECK(not miss_S.contains(m.site));
    DCHECK_EQ(miss_X.get_from_state(m.site, tree->ref_sequence), m.to);
    miss_X.set_from_state(m.site, m.from, tree->ref_sequence);
  }
  
  for (const auto& m : S_P_X_path.hot_muts_to_X) {
    if (m.t > t_P) {
      muts_X.push_back(m);
    } else {
      muts_S.push_back(Mutation{m.to, m.site, m.from, t_P + (t_P - m.t)});

      if (estd::is_debug_enabled) {
        if (ref_to_root_deltas.contains(m.site)) {
          DCHECK_EQ(ref_to_root_deltas.at(m.site).to, m.from);
        } else {
          DCHECK_EQ(tree->ref_sequence[m.site], m.from);
        }
      }
      push_back_site_deltas({m.site, m.from, m.to}, ref_to_root_deltas);
      
      DCHECK(not miss_X.contains(m.site));
      DCHECK(not miss_S.contains(m.site));
    }
  }

  // Fix up ordering and mutations above root
  sort_mutations(muts_X);
  sort_mutations(muts_S);
  muts_P.clear();
  for (const auto& [l, delta] : ref_to_root_deltas) {
    CHECK_EQ(tree->ref_sequence[l], delta.from);
    muts_P.push_back(Mutation{delta.from, l, delta.to, -std::numeric_limits<double>::max()});
  }

  // Fix up any roundoff errors in mutation times
  clamp_mutation_times(muts_X, t_P, t_X);
  clamp_mutation_times(muts_S, t_P, t_S);

  // Recalculate lambda_i at the root
  lambda_i->at(P) = lambda_i->at(X) -
        calc_delta_lambda_across_branch(*evo, tree->ref_sequence, *ref_cum_Q_l,
                                        tree->at(X).mutations, tree->at(X).missations);
}

auto Spr_move::count_rooty_min_mutations(const Spr_graft& graft) -> int {
  const auto& S_P_X_path = graft.branch_infos[Spr_graft::k_branch_info_S_P_X];

  // The mutations on the P-X and P-S paths are "open", meaning their existence is not mandated by
  // the delta between the sequences at S and X
  //const auto& P_X_path = graft.branch_infos[Spr_graft::k_branch_info_P_X];
  //const auto& P_S_path = graft.branch_infos[Spr_graft::k_branch_info_P_S];

  return std::ssize(S_P_X_path.hot_deltas_to_X);
}

auto Spr_move::count_rooty_closed_mutations(const Spr_graft& graft) -> int {
  const auto& S_P_X_path = graft.branch_infos[Spr_graft::k_branch_info_S_P_X];

  // The mutations on the P-X and P-S paths are "open", meaning their existence is not mandated by
  // the delta between the sequences at S and X
  //const auto& P_X_path = graft.branch_infos[Spr_graft::k_branch_info_P_X];
  //const auto& P_S_path = graft.branch_infos[Spr_graft::k_branch_info_P_S];

  return std::ssize(S_P_X_path.hot_muts_to_X);
}

auto Spr_move::summarize_rooty_closed_mutations(const Spr_graft& graft) -> Site_deltas {
  const auto& S_P_X_path = graft.branch_infos[Spr_graft::k_branch_info_S_P_X];

  // The mutations on the P-X and P-S paths are "open", meaning their existence is not mandated by
  // the delta between the sequences at S and X
  //const auto& P_X_path = graft.branch_infos[Spr_graft::k_branch_info_P_X];
  //const auto& P_S_path = graft.branch_infos[Spr_graft::k_branch_info_P_S];

  return S_P_X_path.hot_deltas_to_X;
}

auto Spr_move::start_inner_graft_analysis(Node_index X) const -> Spr_graft {
  //
  // An "inner" graft, as opposed to a "rooty" graft, has hot paths going from X to upstream nodes,
  // possibly all the way to the root, but not a hot path that "rounds" the root.  This makes each hot path
  // a strict extension of its predecessor, with a subset of warm sites.
  //
  // The final path on an inner graft may be open if it reaches the root node and we can change the root node.
  
  CHECK_NE(X, tree->root);
  auto t_X = tree->at(X).t;
  auto P = tree->at(X).parent;
  CHECK_NE(P, tree->root);
  auto t_P = tree->at(P).t;
  auto S = tree->at(P).sibling_of(X);

  auto graft = Spr_graft{};
  graft.X = X;
  graft.S = S;
  graft.t_P = t_P;

  // Map out hot path
  // ----------------
  // P-X path is special because we can't correctly enumerate the warm sites without (expensively) walking
  // all the way up to the root; for all other paths, the warm sites are the interesection of all intervening
  // sibling branches' missations' sites.
  auto& P_X_path = graft.branch_infos.emplace_back();
  P_X_path.A = P;
  P_X_path.B = X;
  P_X_path.is_open = false;
  P_X_path.T_to_X = t_X - t_P;
  P_X_path.warm_sites.insert({0, tree->num_sites()});
  subtract_interval_sets(P_X_path.hot_sites, P_X_path.warm_sites, tree->at(S).missations.intervals);

  auto sliding_missations = tree->at(S).missations;
  P_X_path.partial_lambda_at_A = lambda_i->at(X);
  for (const auto& m : tree->at(X).mutations | std::views::reverse) {
    P_X_path.partial_lambda_at_A += evo->mu_l(m.site) * evo->nu_l.at(m.site) *
        (-evo->q_l_a(m.site, m.to) + evo->q_l_a(m.site, m.from));
  }
  // Remove contribution to partial_lambda missations that will continue sliding
  auto next_partial_lambda_at_B = -1 * calc_delta_lambda_across_missations(*evo, tree->ref_sequence, *ref_cum_Q_l,
                                                                           sliding_missations);
  P_X_path.partial_lambda_at_A -= next_partial_lambda_at_B;
  // partial_lambda_at_X, hot_muts_to_X and hot_deltas_to_X filled in later
  // Beyond this line, do not use P_X_path!  it may become a stale pointer as new branch_infos are added
  
  
  // Paths above P-X
  auto cur = P;
  auto parent = tree->at(cur).parent;
  auto sibling = tree->at(parent).sibling_of(cur);
  auto partial_lambda = next_partial_lambda_at_B;
  auto scratch_sites = Scratch_interval_set{};

  while (not sliding_missations.empty()) {

    // See which sites graduate from "warm" to "hot"
    auto& branch_info = graft.branch_infos.emplace_back();
    branch_info.A = parent;
    branch_info.B = cur;
    branch_info.is_open = false;
    branch_info.T_to_X = t_X - tree->at(parent).t;
    branch_info.warm_sites = sliding_missations.intervals;

    for (const auto& m : tree->at(cur).mutations | std::views::reverse) {
      if (sliding_missations.contains(m.site)) {
        partial_lambda += evo->mu_l(m.site) * evo->nu_l.at(m.site) *
            (-evo->q_l_a(m.site, m.to) + evo->q_l_a(m.site, m.from));
        sliding_missations.set_from_state(m.site, m.from, tree->ref_sequence);
      }
    }
    
    subtract_interval_sets(branch_info.hot_sites, branch_info.warm_sites, tree->at(sibling).missations.intervals);
    
    subtract_interval_sets(sliding_missations.intervals, branch_info.warm_sites, branch_info.hot_sites);
    for (auto it = sliding_missations.from_states.begin(); it != sliding_missations.from_states.end(); /*see below*/) {
      const auto& [l, from] = *it;
      if (not sliding_missations.intervals.contains(l)) {
        it = sliding_missations.from_states.erase(it);
      } else {
        ++it;
      }
    }
    
    next_partial_lambda_at_B = -1 * calc_delta_lambda_across_missations(*evo, tree->ref_sequence, *ref_cum_Q_l,
                                                                        sliding_missations);
    branch_info.partial_lambda_at_A = partial_lambda - next_partial_lambda_at_B;  // *only* includes hot sites
    partial_lambda = next_partial_lambda_at_B;
      
    if (parent != tree->root) {
      // Up one branch 
      cur = parent;
      parent = tree->at(cur).parent;
      sibling = tree->at(parent).sibling_of(cur);
      
    } else {
      if (not can_change_root) {
        // Stop sliding missations instead of going around the root
        branch_info.hot_sites = branch_info.warm_sites;
        branch_info.partial_lambda_at_A += partial_lambda;
      } else {
        if (not sliding_missations.empty()) {
          auto& final_open_branch_info = graft.branch_infos.emplace_back();
          final_open_branch_info.A = k_no_node;
          final_open_branch_info.B = tree->root;
          final_open_branch_info.is_open = true;
          final_open_branch_info.T_to_X = t_X - tree->at(parent).t;
          final_open_branch_info.warm_sites = sliding_missations.intervals;
          final_open_branch_info.hot_sites = final_open_branch_info.warm_sites;
          final_open_branch_info.partial_lambda_at_A = partial_lambda;
        }
      }
      sliding_missations.clear();  // We're done...
    }
  }

  // Distribute hot mutations along hot path
  // (accumulate in reverse order, then reverse)
  for (auto i = 0; i != std::ssize(graft.branch_infos); ++i) {
    auto& branch_info_i = graft.branch_infos[i];

    if (branch_info_i.B == tree->root) { continue; }  // Mutations above the root aren't real mutations
    for (const auto& m : tree->at(branch_info_i.B).mutations | std::views::reverse) {
      if (branch_info_i.warm_sites.contains(m.site)) {
        // Search for upstream branch_info where m.site is hot; almost always, it's this branch_info
        auto found = false;
        for (auto j = i; j != std::ssize(graft.branch_infos); ++j) {
          auto& branch_info_j = graft.branch_infos[j];
          
          if (branch_info_j.hot_sites.contains(m.site)) {
            branch_info_j.hot_muts_to_X.push_back(m);
            found = true;
          }
        }
        CHECK(found);
      }
    }
  }

  // Place mutations in the correct order and derive deltas and partial_lambda_X from them
  for (auto& branch_info : graft.branch_infos) {
    std::ranges::reverse(branch_info.hot_muts_to_X);
    
    branch_info.partial_lambda_at_X = branch_info.partial_lambda_at_A;
    for (const auto& m : branch_info.hot_muts_to_X) {
      if (not branch_info.is_open) {
        push_back_site_deltas(m, branch_info.hot_deltas_to_X);  // No deltas for open paths
      }
      branch_info.partial_lambda_at_X += evo->mu_l(m.site) * evo->nu_l.at(m.site) *
          (-evo->q_l_a(m.site, m.from) + evo->q_l_a(m.site, m.to));
    }
  }

  return graft;
}

auto Spr_move::propose_new_inner_graft_mutations(Spr_graft& graft, absl::BitGenRef bitgen) const -> void {
  auto X = graft.X;
  auto path_end = tree->node_loc(X);
  
  for (auto& branch_info : graft.branch_infos) {
    
    if (branch_info.hot_sites.empty()) {
      CHECK(branch_info.hot_muts_to_X.empty());
      continue;
    }
    
    // TODO: if path is open at root, do not constrain state at start
    auto new_mutations =
        branch_info.is_open
        ? sample_unconstrained_mutational_history(
            tree->num_sites(), branch_info.T_to_X, mu_proposal, bitgen)
        : sample_mutational_history(
              tree->num_sites(), branch_info.T_to_X, mu_proposal,
              branch_info.hot_deltas_to_X, bitgen);
    if (not new_mutations.empty()) {
      std::erase_if(new_mutations, [&](const auto& m) { return not branch_info.hot_sites.contains(m.site); });
      if (branch_info.B == X) {
        // In the very first branch, `hot_sites` may cover sites that are actually missing at X
        // (but we don't know about them because the missations are further above X than we every look).
        // Any mutations that show up in new_mutations but not in deltas_to_X_on_hot_sites are very rare,
        // so we can afford an expensive check for them
        std::erase_if(new_mutations, [&](const auto& m) {
          if (branch_info.hot_deltas_to_X.contains(m.site)) { return false; }
          return is_site_missing_at(*tree, X, m.site);
        });
      }
      
      adjust_mutational_history(new_mutations, branch_info.hot_deltas_to_X, *tree, path_end);
    }
    
    branch_info.hot_muts_to_X = std::move(new_mutations);

    // For open paths, partial_lambda_at_A may have changed
    if (branch_info.is_open) {
      branch_info.partial_lambda_at_A = branch_info.partial_lambda_at_X;
      for (const auto& m : branch_info.hot_muts_to_X | std::views::reverse) {
        branch_info.partial_lambda_at_A += evo->mu_l(m.site) * evo->nu_l.at(m.site) *
            (+evo->q_l_a(m.site, m.from) - evo->q_l_a(m.site, m.to));
      }
    }
  }
}

auto Spr_move::finish_inner_graft_analysis(Spr_graft& graft) const -> void {
  auto X = graft.X;
  auto t_X = tree->at(X).t;
  
  // Contribution of hot paths to log_G
  graft.delta_log_G = 0.0;
  for (const auto& branch_info : graft.branch_infos) {
    graft.delta_log_G += calc_branch_log_G(t_X - branch_info.T_to_X, t_X, branch_info.partial_lambda_at_X, *evo,
                                           branch_info.hot_muts_to_X);
  }

  // Possible complication: root sequence changes from final open path
  if (graft.branch_infos.back().is_open) {
    CHECK(can_change_root);
    // Imagine slipping in the mutations right-to-left from above the root into to their final locations
    for (const auto& m : graft.branch_infos.back().hot_muts_to_X) {
      graft.delta_log_G += std::log(evo->pi_l_a(m.site, m.from) / evo->pi_l_a(m.site, m.to));
    }
  }

  // Probability of choosing this particular mutational history (excluding factors common to all histories)
  graft.log_alpha_mut = 0.0;
  for (const auto& branch_info : graft.branch_infos) {
    auto L = branch_info.hot_sites.num_sites();
    if (branch_info.B == X) {
      // Slight complication: we're not really sure that all the hot sites on the P->X path
      // aren't actually missing at X.  But we can still derive the total number of hot sites on the P->X path.
      L = (tree->num_sites() - num_sites_missing_at_every_node->at(X))  // present at X
          - (branch_info.warm_sites.num_sites() - branch_info.hot_sites.num_sites());  // hot above P
    }

    auto T = branch_info.T_to_X;
    auto M = std::ssize(branch_info.hot_muts_to_X);
    
    graft.log_alpha_mut +=
        -mu_proposal * L * T
        + M * std::log(mu_proposal / 3);

    if (not branch_info.is_open) {
      auto d = std::ssize(branch_info.hot_deltas_to_X);
      
      auto P_AC_JC = -0.25 * std::expm1(-4./3.*mu_proposal * T);
      auto log_P_AC_JC = std::log(P_AC_JC);
      //auto P_AA_JC = 1.0 + 3 * P_AC_JC;
      auto log_P_AA_JC = std::log1p(-3 * P_AC_JC);
      
      graft.log_alpha_mut -=(L-d)*log_P_AA_JC + d*log_P_AC_JC;
    }
  }
}

auto Spr_move::peel_inner_graft(const Spr_graft& graft) -> void {
  // The tricky bit here is to fix up the missations
  auto X = graft.X;
  CHECK_NE(X, tree->root);
  auto t_X = tree->at(X).t;
  auto P = tree->at(X).parent;
  CHECK_NE(P, k_no_node);
  auto t_P = tree->at(P).t;

  auto& muts_X = tree->at(X).mutations;

  // First, remove all the mutations on the hot paths.  Effectively, we'll slide them off right-to-left towards the
  // end of each path first, which explains the effects on the `from_states` on missations in offshoot branches.
  //
  // Exception: for the final path, if open, we slide the mutations off left-to-right to above the root,
  // which may change the root sequence
  auto& final_path = graft.branch_infos.back();
  auto ref_to_root_deltas = Site_deltas{};
  if (final_path.is_open) {
    for (const auto& m : tree->at_root().mutations) {
      push_back_site_deltas(m, ref_to_root_deltas);
    }
  }
  
  for (const auto& branch_info : graft.branch_infos) {
    if (branch_info.B == tree->root) { continue; }  // Mutations above root are not real mutations
    if (branch_info.B == X and not final_path.is_open) {
      // Optimization: we have to clear *all* the mutations on this branch, and no missations are affected
      tree->at(X).mutations.clear();
    } else {

      // Mutations not on a final open path are slid right-to-left downstream towards the P-X branch
      for (auto& m : tree->at(branch_info.B).mutations | std::views::reverse) {
        if (branch_info.warm_sites.contains(m.site)) {
          if (not final_path.is_open || not final_path.hot_sites.contains(m.site)) {
            // Slide this mutations downstream onto P-X branch, affecting the from_state of all
            // sibling missations along the way
            for (auto cur = X; cur != branch_info.B; cur = tree->at(cur).parent) {
              auto parent = tree->at(cur).parent;
              auto sibling = tree->at(parent).sibling_of(cur);
              auto& miss_sibling = tree->at(sibling).missations;
              DCHECK(miss_sibling.contains(m.site));
              DCHECK_EQ(miss_sibling.get_from_state(m.site, tree->ref_sequence), m.to);
              miss_sibling.set_from_state(m.site, m.from, tree->ref_sequence);
            }
            m.site = -1;  // Mark for removal
          }
        }
      }
    }
  }

  if (final_path.is_open) {
    for (const auto& branch_info : graft.branch_infos | std::views::reverse) {
      if (branch_info.B == tree->root) { continue; }  // Mutations above root are not real mutations
      
      // Mutations on a final open path are slid left-to-right upstream above the root
      for (auto& m : tree->at(branch_info.B).mutations) {
        if (final_path.hot_sites.contains(m.site)) {
          // Slide this mutations upstream past the root onto P-X branch, affecting the from_state of all
          // sibling missations along the way and the root sequence
          for (auto cur = branch_info.B; cur != tree->root; cur = tree->at(cur).parent) {
            auto parent = tree->at(cur).parent;
            auto sibling = tree->at(parent).sibling_of(cur);
            auto& miss_sibling = tree->at(sibling).missations;
            DCHECK(miss_sibling.contains(m.site));
            DCHECK_EQ(miss_sibling.get_from_state(m.site, tree->ref_sequence), m.from);
            miss_sibling.set_from_state(m.site, m.to, tree->ref_sequence);
          }
          
          // Change in root sequence
          push_back_site_deltas(m, ref_to_root_deltas);
          
          m.site = -1;  // Mark for removal
        }
      }
    }
  }

  for (const auto& branch_info : graft.branch_infos) {
    if (branch_info.B == tree->root) { continue; }  // Mutations above root are not real mutations
    
    // Remove marked mutations
    std::erase_if(tree->at(branch_info.B).mutations, [&](const auto& m) { return m.site == -1; });
  }

  // Temporarily place fake mutations on the P-X branch to restore the tree's consistency
  auto t_mut_X = 0.5*(t_P + t_X);
  for (const auto& branch_info : graft.branch_infos) {
    if (branch_info.B == tree->root) { continue; }  // Mutations above root are not real mutations
    for (const auto& [l, delta] : branch_info.hot_deltas_to_X) {
      muts_X.push_back(Mutation{delta.from, l, delta.to, t_mut_X});
    }
  }
  
  // Fix up mutations above root
  if (final_path.is_open) {
    auto& muts_root = tree->at_root().mutations;
    muts_root.clear();
    for (const auto& [l, delta] : ref_to_root_deltas) {
      CHECK_EQ(tree->ref_sequence[l], delta.from);
      muts_root.push_back(Mutation{delta.from, l, delta.to, -std::numeric_limits<double>::max()});
    }
  }

  // Recalculate lambda_i along the hot path
  for (const auto& branch_info : graft.branch_infos | std::views::take(std::ssize(graft.branch_infos) - 1)) {
    auto A = branch_info.A;
    auto B = branch_info.B;
    lambda_i->at(A) = lambda_i->at(B) -
        calc_delta_lambda_across_branch(*evo, tree->ref_sequence, *ref_cum_Q_l,
                                        tree->at(B).mutations, tree->at(B).missations);
  }

  assert_phylo_tree_integrity(*tree);
}

auto Spr_move::apply_inner_graft(const Spr_graft& graft) -> void {
  // The tricky bit here is to fix up the missations
  auto X = graft.X;
  CHECK_NE(X, tree->root);
  auto& muts_X = tree->at(X).mutations;

  auto& final_path = graft.branch_infos.back();
  
  // We assume the graft was either peeled before, or we're here after an `move` operation,
  // which never inserts spurious mutations outside the P-X branch.
  
  // We're about to put in place all the recorded mutations on the P-X branch, so nuke what's already there
  muts_X.clear();
  
  for (const auto& branch_info : graft.branch_infos) {
    if (branch_info.B == tree->root) { continue; }  // Mutations above root are not real mutations
    for (const auto& m : tree->at(branch_info.B).mutations) {
      DCHECK(not branch_info.warm_sites.contains(m.site));
    }
  }
  if (final_path.is_open) {
    for (const auto& branch_info : graft.branch_infos) {
      if (branch_info.B == tree->root) { continue; }  // Mutations above root are not real mutations
      for (const auto& m : tree->at(branch_info.B).mutations) {
        DCHECK(not final_path.hot_sites.contains(m.site));
      }
    }
  }
  
  // We may have to change the root sequence if the final path is open
  auto ref_to_root_deltas = Site_deltas{};
  if (final_path.is_open) {
    for (const auto& m : tree->at_root().mutations) {
      push_back_site_deltas(m, ref_to_root_deltas);
    }
  }

  // Now insert all the new mutations.  Effectively, we first place them left-to-right just above X, then slide them
  // down to their final position, which explains the effects on the `from_states` on missations in offshoot branches.
  for (const auto& branch_info : graft.branch_infos) {
    if (branch_info.B == X) {
      // Optimization: we have to set *all* the mutations on this branch, and no missations are affected
      tree->at(X).mutations.assign(branch_info.hot_muts_to_X.begin(), branch_info.hot_muts_to_X.end());
    } else {
      if (not branch_info.is_open) {
        for (const auto& m : branch_info.hot_muts_to_X) {
          // Slide this mutations upstream to its final position, affecting the from_state of all
          // sibling missations along the way
          for (auto cur = X; cur != branch_info.A; cur = tree->at(cur).parent) {
            auto parent = tree->at(cur).parent;
            if (tree->at(parent).t <= m.t && m.t < tree->at(cur).t) {
              tree->at(cur).mutations.push_back(m); // Will sort below
              break;
            } else {
              auto sibling = tree->at(parent).sibling_of(cur);
              auto& miss_sibling = tree->at(sibling).missations;
              DCHECK(miss_sibling.contains(m.site));
              DCHECK_EQ(miss_sibling.get_from_state(m.site, tree->ref_sequence), m.from);
              miss_sibling.set_from_state(m.site, m.to, tree->ref_sequence);
            }
          }
        }
      } else {  // branch_info.is_open == true
        for (const auto& m : branch_info.hot_muts_to_X | std::views::reverse) {
          // Slide this mutations downstream to its final position, affecting the from_state of all
          // sibling missations along the way and changing the root sequence
          for (auto cur = X; cur != tree->root; cur = tree->at(cur).parent) {
            auto parent = tree->at(cur).parent;
            if (tree->at(parent).t <= m.t && m.t < tree->at(cur).t) {
              tree->at(cur).mutations.push_back(m); // Will sort below
            }
            if (tree->at(parent).t <= m.t) {
              auto sibling = tree->at(parent).sibling_of(cur);
              auto& miss_sibling = tree->at(sibling).missations;
              DCHECK(miss_sibling.contains(m.site));
              DCHECK_EQ(miss_sibling.get_from_state(m.site, tree->ref_sequence), m.to);
              miss_sibling.set_from_state(m.site, m.from, tree->ref_sequence);
            }
          }

          push_back_site_deltas({m.site, m.to, m.from}, ref_to_root_deltas);
        }
      }
    }
  }

  // Sort all affected branches' mutations and clean up any roundoff errors
  for (const auto& branch_info : graft.branch_infos) {
    if (not branch_info.is_open) {
      auto t_A = tree->at(branch_info.A).t;
      auto t_B = tree->at(branch_info.B).t;
      sort_mutations(tree->at(branch_info.B).mutations);
      clamp_mutation_times(tree->at(branch_info.B).mutations, t_A, t_B);
    }
  }

  // Adjust root sequence if needed
  if (final_path.is_open) {
    tree->at_root().mutations.clear();
    for (const auto& [l, delta] : ref_to_root_deltas) {
      tree->at_root().mutations.push_back(Mutation{delta.from, l, delta.to, -std::numeric_limits<double>::max()});
    }
  }

  // Recalculate lambda_i along the hot path
  for (const auto& branch_info : graft.branch_infos | std::views::take(std::ssize(graft.branch_infos) - 1)) {
    auto A = branch_info.A;
    auto B = branch_info.B;
    lambda_i->at(A) = lambda_i->at(B) -
        calc_delta_lambda_across_branch(*evo, tree->ref_sequence, *ref_cum_Q_l,
                                        tree->at(B).mutations, tree->at(B).missations);
  }

  assert_phylo_tree_integrity(*tree);
}

auto Spr_move::count_inner_min_mutations(const Spr_graft& graft) -> int {
  auto result = 0;
  for (const auto& branch_info : graft.branch_infos) {
    if (not branch_info.is_open) {
      result += std::ssize(branch_info.hot_deltas_to_X);
    }
  }
  return result;
}

auto Spr_move::count_inner_closed_mutations(const Spr_graft& graft) -> int {
  auto result = 0;
  for (const auto& branch_info : graft.branch_infos) {
    if (not branch_info.is_open) {
      result += std::ssize(branch_info.hot_muts_to_X);
    }
  }
  return result;
}

auto Spr_move::summarize_inner_closed_mutations(const Spr_graft& graft) -> Site_deltas {
  auto result = Site_deltas{};
  for (const auto& branch_info : graft.branch_infos) {
    if (not branch_info.is_open) {
      append_site_deltas(result, branch_info.hot_deltas_to_X); // OK because no site is hot in more than one branch_info
    }
  }
  return result;
}

auto Spr_move::move(Node_index X, Node_index SS, double new_t_P) -> void {
  CHECK_NE(X, tree->root);
  DCHECK(not descends_from(*tree, SS, X));
  auto P = tree->at(X).parent;
  auto G = tree->at(P).parent;  // could be k_no_node
  auto S = tree->at(P).sibling_of(X);
  if (SS == P) {
    // We always want to describe the regrafting in terms of the tree *without* the pruned subtree
    SS = S;
  }

  auto GG = tree->at(SS).parent;  // could be k_no_node
  if (GG == P) {
    GG = G;  // Want GG to be SS's parent in the tree *without* the pruned subtree
  }

  auto A = find_MRCA_of(*tree, G, GG);

  // 1. First, we displace X's attachment point until X and SS are siblings
  auto edit = Tree_editing_session{*tree, X, *evo, *lambda_i, *ref_cum_Q_l, *num_sites_missing_at_every_node};

  // Slide X up just below A
  while (tree->at(P).parent != A) {
    edit.slide_P_along_branch(tree->at_parent_of(P).t);
    edit.hop_up();
  }

  // Flip to other branch if needed
  if (not (descends_from(*tree, S, SS) || descends_from(*tree, SS, S))) {
    CHECK_NE(A, k_no_node);
    edit.slide_P_along_branch(tree->at(A).t);
    edit.flip();
  }

  // Slide X down to become SS's sibling
  DCHECK(descends_from(*tree, SS, P));

  auto branches_to_SS = Scratch_vector<Node_index>{};
  auto Xs_current_sibling = tree->at(P).sibling_of(X);
  for (auto cur = SS; cur != Xs_current_sibling; cur = tree->at(cur).parent) {
    branches_to_SS.push_back(cur);
  }

  for (auto Y : branches_to_SS | std::views::reverse) {
    CHECK_NE(tree->at(SS).parent, P);
    edit.slide_P_along_branch(tree->at_parent_of(Y).t);
    edit.hop_down(Y);
  }

  // X is now attached to the correct branch, just have to slide it into the right place
  CHECK_EQ(tree->at(P).sibling_of(X), SS);
  edit.slide_P_along_branch(new_t_P);

  // Done!
  edit.end();
}

static auto choose_different_state(Real_seq_letter s, absl::BitGenRef bitgen) {
  using enum Real_seq_letter;
  auto delta = absl::Uniform(absl::IntervalClosedOpen, bitgen, 1, k_num_real_seq_letters);
  return real_seq_letter_from_index((index_of(s) + delta) % k_num_real_seq_letters);
}

auto sample_mutational_history(
    Site_index L,
    double T,
    double mu,
    const Site_deltas& deltas,
    absl::BitGenRef bitgen)
    -> Scratch_vector<Mutation> {

  auto result = Scratch_vector<Mutation>{};
  
  // For sites with deltas, sample CTMC trajectories starting with `from` state with at least one mutation
  // until we get a trajectory with state `to` at the end.  This is equivalent to doing straight rejection
  // sampling over unconstrained trajectories, but skipping quickly over all proposed trajectories with 0 mutations.
  auto to_states = Scratch_vector<Real_seq_letter>{};
  auto mut_times = Scratch_vector<double>{};
  auto trajectory = Scratch_vector<Mutation>{};
  auto num_muts_dist_ge1 = K_truncated_poisson_distribution{mu*T, 1};
  for (const auto& [l, delta] : deltas) {

    auto n = 0;
    while (true) {

      // Trajectory has n >= 1 mutations:
      // * from delta.from to to_states[0]
      // * from to_states[0] to to_states[1]
      // ...
      // * from to_states[n-2] to to_states[n-1]
      
      n = num_muts_dist_ge1(bitgen);
      to_states.clear();
      auto s = delta.from;
      for (auto i = 0; i != n; ++i) {
        s = choose_different_state(s, bitgen);
        to_states.push_back(s);
      }

      if (s == delta.to) {
        // Accept!
        break;
      } else {
        // Reject and start over!
        continue;
      }
    }
    
    // Times are uniformly distributed across trajectory
    mut_times.clear();
    for (auto i = 0; i != n; ++i) {
      mut_times.push_back(absl::Uniform(absl::IntervalClosedOpen, bitgen, -T, 0.0));
    }
    std::ranges::sort(mut_times);

    auto prev_s = delta.from;
    for (auto i = 0; i != n; ++i) {
      auto next_s = to_states[i];
      auto t = mut_times[i];

      result.push_back(Mutation{prev_s, l, next_s, t});
      
      prev_s = next_s;
    }
  }

  // Conceptually, for sites without deltas, whose start and end states are equal, we also do rejection sampling as
  // above, one site at a time.  However, for the vast majority of the sites, the first proposal of the first mutation
  // time already has t > 0, so the (empty) trajectory is accepted.  If not, the resulting trajectory will almost
  // certainly be rejected: we'd need at least *two* mutations with times between -T and 0, one that mutates away from the
  // starting state and one to return.
  //
  // Given the above, we conceptually iterate over every site, and perform rejection sampling until either an empty
  // trajectory results, or we realize that we're creating a proposed trajectory with 2+ mutations.  For almost every site,
  // we'll just end up sampling an empty trajectory.  Thus, we can easily skip the next `k` such sites until an incipient
  // proposal with 2+ mutations arises.  Once we know that the proposal will have at least 2 mutations, we continue rejection
  // sampling as usual.  The first proposal is completed and accepted if its end state matches its starting
  // state (this happens about 1/3 of the time).  Otherwise, an *entirely new* trajectory is proposed from scratch,
  // which will almost certainly be an empty trajectory.  Hence, the chosen sites do not necessarily end up having
  // mutations.
  //
  // Almost all of the time, this procedure will immediately skip over all sites in the genome, so it's very efficient.
  //
  // In detail:
  // * let p_0 = e^(-mu T) be the probability that the proposed trajectory has 0 mutations
  // * let p_1 = mu T e^(-mu T) be the probability that the proposed trajectory has 1 mutation
  //
  // The probability of rejection sampling eventually proposing a trajectory with 0 mutations without ever proposing
  // one with 2 or more mutations is:
  //
  //  p_* = sum_{k=0}^infty p_1^k p_0 = p_0 / (1 - p_1).
  //
  // Conversely, the probability that rejection sampling at a particular site hits the point where a 2+-mutation
  // trajectory should be proposed is:
  //
  //  1 - p_* = (1 - p0 - p1) / (1 - p1)
  //
  // During rejection sampling over all sites, the number of sites for which empty trajectories will be accepted without
  // a single 2+-mutation trajectory being proposed is distributed as Geo(p_*).  We thus repeatedly pick how many sites
  // to skip over before slowing down and following through on one iteration of rejection sampling at a specific site.
  // For viral-like parameters in a genomic epi setting (e.g., L = 30,000, mu = 1e-3 / year, T = ~2 weeks), we overwhelmingly
  // immediately skip over all sites.
  //
  // Conceptually, we generate A->A trajectories for all L sites, and then filter out and mutations on sites which
  // have explicit deltas, which we treated above.  Equivalently, if we ever stop at a site that has an explicit delta,
  // we immediately skip it.

  auto p_0 = std::exp(-mu*T);
  auto p_1 = mu * T * p_0;
  //auto p_tricky               = (1 - p_0 - p_1) /          (1-p_1);
  //auto one_minus_p_tricky     =      p_0        /          (1-p_1);
  //auto log_one_minus_p_tricky = std::log(p_0)   - std::log1p(-p_1);
  auto log_one_minus_p_tricky   =     -mu*T       - std::log1p(-p_1);
  
  auto l = 0;
  while (true) {
    // We really want to do `delta = std::geometric_distribution{p_tricky}(bitgen); l += delta`, but since
    // p_tricky is very small, the sample drawn might be large enough to overflow integers.
    //
    // So we break down the sampling of delta into an exponential sample followed by flooring, and exit early
    // if the result of the exponential implies that delta >= L.
    //
    // Concretely, if u ~ Expo(-ln(1-p)), then floor(u) ~ Geo(p).  So if u >= L, then we can do an early exit
    auto u = std::exponential_distribution{-log_one_minus_p_tricky}(bitgen);
    if (u >= L) {
      break;  // when we turn u into a Geo(p_tricky) sample and add it to l, we'd have l >= L
    }
    auto delta = std::floor(u);
    l += delta;
    if (l >= L) {
      break;
    }

    if (deltas.contains(l)) {
      // Would have filtered afterwards
      ++l;  // Restart rejection sampling from the *next* site
      continue;
    }

    // Propose trajectory starting from A with n >= 2 mutations:
    // * from A to to_states[0]
    // * from to_states[0] to to_states[1]
    // ...
    // * from to_states[n-2] to to_states[n-1]
    
    auto n = K_truncated_poisson_distribution{mu*T, 2}(bitgen);
    to_states.clear();
    auto s = Real_seq_letter::A;
    for (auto i = 0; i != n; ++i) {
      s = choose_different_state(s, bitgen);
      to_states.push_back(s);
    }
    
    if (s == Real_seq_letter::A) {
      // Accept!
      // Times are uniformly distributed across trajectory
      mut_times.clear();
      for (auto i = 0; i != n; ++i) {
        mut_times.push_back(absl::Uniform(absl::IntervalClosedOpen, bitgen, -T, 0.0));
      }
      std::ranges::sort(mut_times);
      
      auto prev_s = Real_seq_letter::A;
      for (auto i = 0; i != n; ++i) {
        auto next_s = to_states[i];
        auto t = mut_times[i];
        
        result.push_back(Mutation{prev_s, l, next_s, t});
        
        prev_s = next_s;
      }
      
      // We're done with this site; restart rejection sampling at the next site
      ++l;
    } else {
      // Reject!  Continue rejection sampling from *this* site
      continue;
    }
  }

  // Finally, interleave all the mutations in time
  sort_mutations(result);
  
  return result;
}

auto sample_unconstrained_mutational_history(
    Site_index L,
    double T,
    double mu,
    absl::BitGenRef bitgen)
    -> Scratch_vector<Mutation> {

  auto result = Scratch_vector<Mutation>{};

  // This is nearly trivial with the Gillespie algorithm
  auto cur_state_of_site = Scratch_flat_hash_map<Site_index, Real_seq_letter>{};

  // We accumulate mutations in reverse order...
  auto trajectory = Scratch_vector<Mutation>{};
  auto t = 0.0;
  while (true) {
    t -= absl::Exponential(bitgen, mu*L);
    if (t <= -T) { break; }

    auto l = absl::Uniform(absl::IntervalClosedOpen, bitgen, 0, L);

    auto [it, inserted] = cur_state_of_site.try_emplace(l, Real_seq_letter::A);
    auto s = it->second;
    CHECK(not inserted || s == Real_seq_letter::A);
    auto next_s = choose_different_state(s, bitgen);
    
    trajectory.push_back(Mutation{next_s, l, s, t});  // Remember, we're adding mutations right-to-left
    
    it->second = next_s;
  }

  // ...then reverse the result at the end
  std::ranges::reverse(trajectory);

  return trajectory;
}

auto adjust_mutational_history(
    Scratch_vector<Mutation>& history,
    const Site_deltas site_deltas,
    const Phylo_tree& tree,
    Phylo_tree_loc end_loc)
    -> void {

  auto end_states_of_non_mutated_sites = Scratch_flat_hash_map<Site_index, Real_seq_letter>{};
  
  for (auto& m : history | std::views::reverse) {
    
    m.t += end_loc.t;
    
    if (not site_deltas.contains(m.site)) {
      // This code path is rare
      auto end_state = Real_seq_letter::A;
      auto [it, inserted] = end_states_of_non_mutated_sites.try_emplace(m.site, end_state);
      if (not inserted) {
        end_state = it->second;
      } else {
        // inserted end state is wrong, need to correct it
        end_state = it->second = calc_site_state_at(tree, end_loc, m.site);
      }

      // Adjust mutation
      auto delta = index_of(end_state) - index_of(Real_seq_letter::A);  // A = 0, so delta >= 0 always
      m.from = real_seq_letter_from_index((index_of(m.from) + delta) % k_num_real_seq_letters);
      m.to = real_seq_letter_from_index((index_of(m.to) + delta) % k_num_real_seq_letters);
    }
  }
}

}  // namespace delphy
