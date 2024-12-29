#include "tree_editing.h"

#include <absl/log/check.h>

namespace delphy {

Tree_editing_session::Tree_editing_session(
      Phylo_tree& tree,
      Node_index X,
      const Global_evo_model& evo,
      Node_vector<double>& lambda_i,
      const std::vector<double>& ref_cum_Q_l,
      Node_vector<int>& num_sites_missing_at_every_node)
    : tree{&tree},
      X{X},
      evo{&evo},
      lambda_i{&lambda_i},
      ref_cum_Q_l{&ref_cum_Q_l},
      num_sites_missing_at_every_node{&num_sites_missing_at_every_node},
      deltas_nexus_to_X{} {

  CHECK_NE(X, this->tree->root);
  
  // Remove mutations above X
  for (const auto& m : this->tree->at(X).mutations) {
    push_back_site_deltas(m, deltas_nexus_to_X);
  }
  this->tree->at(X).mutations.clear();
}

auto Tree_editing_session::slide_P_along_branch(double new_t_P) -> void {
  auto P = tree->at(X).parent;
  CHECK(not tree->at(P).is_tip());

  if (P == tree->root) {
    // Sliding the root is a completely different ballgame
    return slide_root(new_t_P);
  }
  
  auto t_P = tree->at(P).t;
  auto t_G = (P == tree->root ? -std::numeric_limits<double>::max() : tree->at_parent_of(P).t);
  auto old_t_P = tree->at(P).t;
  auto S = tree->at(P).sibling_of(X);
  auto t_S = tree->at(S).t;

  DCHECK(are_mutation_times_in_range(tree->at(P).mutations, t_G, t_P))
      << P << "\n"
      << [&](){ std::cerr << *tree << "\n"; return "\n"; }();
  DCHECK(are_mutation_times_in_range(tree->at(S).mutations, t_P, t_S));

  if (new_t_P < old_t_P) {
    // Moving P upstream: move last mutations from G-P into P-S into child branches

    auto G_P_muts_after_new_t_P_begin =
        std::ranges::find_if(tree->at(P).mutations, [&](const auto& m) { return m.t >= new_t_P; });
    auto G_P_muts_after_new_t_P_end =
        tree->at(P).mutations.end();

    if (G_P_muts_after_new_t_P_begin != G_P_muts_after_new_t_P_end) {
      std::ranges::reverse(tree->at(S).mutations);
      for (const auto& m : std::ranges::subrange(G_P_muts_after_new_t_P_begin, G_P_muts_after_new_t_P_end)
               | std::views::reverse) {
        
        if (not tree->at(S).missations.contains(m.site)) {
          tree->at(S).mutations.push_back(m);
        } else {
          // This code path is *never* triggerred in SPR moves because graft peeling removes these mutations
          DCHECK_EQ(tree->at(S).missations.get_from_state(m.site, tree->ref_sequence), m.to);
          tree->at(S).missations.set_from_state(m.site, m.from, tree->ref_sequence);
        }
        
        if (not tree->at(X).missations.contains(m.site)) {
          push_front_site_deltas(m, deltas_nexus_to_X);
        } else {
          DCHECK_EQ(tree->at(X).missations.get_from_state(m.site, tree->ref_sequence), m.to);
          tree->at(X).missations.set_from_state(m.site, m.from, tree->ref_sequence);
        }

        lambda_i->at(P) += evo->mu_l(m.site) * evo->nu_l.at(m.site) *
          (-evo->q_l_a(m.site, m.to) + evo->q_l_a(m.site, m.from));
      }
      std::ranges::reverse(tree->at(S).mutations);
      tree->at(P).mutations.erase(G_P_muts_after_new_t_P_begin, G_P_muts_after_new_t_P_end);
    }
  } else {
    // Moving P downstream: move first mutations from P-S into G-P
    auto P_S_muts_before_new_t_P_begin =
        tree->at(S).mutations.begin();
    auto P_S_muts_before_new_t_P_end =
        std::ranges::find_if(tree->at(S).mutations, [&](const auto& m) { return m.t > new_t_P; });
    if (P_S_muts_before_new_t_P_begin != P_S_muts_before_new_t_P_end) {
      for (const auto& m : std::ranges::subrange(P_S_muts_before_new_t_P_begin, P_S_muts_before_new_t_P_end)) {
        
        DCHECK(not tree->at(S).missations.contains(m.site));
        tree->at(P).mutations.push_back(m);
        
        if (not tree->at(X).missations.contains(m.site)) {
          push_front_site_deltas({m.site, m.to, m.from}, deltas_nexus_to_X);
        } else {
          DCHECK_EQ(tree->at(X).missations.get_from_state(m.site, tree->ref_sequence), m.from);
          tree->at(X).missations.set_from_state(m.site, m.to, tree->ref_sequence);
        }

        lambda_i->at(P) += evo->mu_l(m.site) * evo->nu_l.at(m.site) *
          (-evo->q_l_a(m.site, m.from) + evo->q_l_a(m.site, m.to));
      }
      tree->at(S).mutations.erase(P_S_muts_before_new_t_P_begin, P_S_muts_before_new_t_P_end);
    }
  }
  
  tree->at(P).t = new_t_P;
}

auto Tree_editing_session::slide_root(double new_t_P) -> void {
  auto P = tree->at(X).parent;
  CHECK_EQ(P, tree->root);
  auto old_t_P = tree->at(P).t;
  auto S = tree->at(P).sibling_of(X);
  
  // Only sliding down is hard
  if (new_t_P > old_t_P) {
    // Moving root downstream: move first mutations from P-S into G-P
    auto P_S_muts_before_new_t_P_begin = tree->at(S).mutations.begin();
    auto P_S_muts_before_new_t_P_end =
        std::ranges::find_if(tree->at(S).mutations, [&](const auto& m) { return m.t > new_t_P; });
    if (P_S_muts_before_new_t_P_begin != P_S_muts_before_new_t_P_end) {

      auto ref_to_root_deltas = Site_deltas{};
      for (const auto& m : tree->at(P).mutations) {
        push_back_site_deltas(m, ref_to_root_deltas);
      }
      
      for (const auto& m : std::ranges::subrange(P_S_muts_before_new_t_P_begin, P_S_muts_before_new_t_P_end)) {
        
        DCHECK(not tree->at(S).missations.contains(m.site));
        push_back_site_deltas(m, ref_to_root_deltas);  // "tree->at(P).mutations.push_back(m);"
        
        if (not tree->at(X).missations.contains(m.site)) {
          push_front_site_deltas({m.site, m.to, m.from}, deltas_nexus_to_X);
        } else {
          DCHECK_EQ(tree->at(X).missations.get_from_state(m.site, tree->ref_sequence), m.from);
          tree->at(X).missations.set_from_state(m.site, m.to, tree->ref_sequence);
        }

        lambda_i->at(P) += evo->mu_l(m.site) * evo->nu_l.at(m.site) *
          (-evo->q_l_a(m.site, m.from) + evo->q_l_a(m.site, m.to));
      }

      tree->at(P).mutations.clear();
      for (const auto& [l, delta] : ref_to_root_deltas) {
        tree->at(P).mutations.push_back(Mutation{delta.from, l, delta.to, -std::numeric_limits<double>::max()});
      }
      tree->at(S).mutations.erase(P_S_muts_before_new_t_P_begin, P_S_muts_before_new_t_P_end);
    }
  }
  
  tree->at(P).t = new_t_P;
}

auto Tree_editing_session::hop_up() -> void {
  do_hop_up(this->X);
}

auto Tree_editing_session::do_hop_up(Node_index X) -> void {
  // NOTE: parameter `X` intentionally shadows this->X
  
  CHECK_NE(X, tree->root);
  auto P = tree->at(X).parent;
  CHECK(not tree->at(P).is_tip());
  CHECK_NE(P, tree->root);
  CHECK(tree->at(P).mutations.empty());

  auto G = tree->at(P).parent;
  CHECK_EQ(tree->at(P).t, tree->at(G).t);
  auto U = tree->at(G).sibling_of(P);

  auto S = tree->at(P).sibling_of(X);

  // Shuffle mutations and missations around G & P (see picture in tree_editing.h)
  // -----------------------------------------------------------------------------
  
  // 1. Push P's missations down to children
  if (not tree->at(P).missations.empty()) {
    tree->at(X).missations = merge_missations_nondestructively(tree->at(X).missations, tree->at(P).missations);
    tree->at(S).missations = merge_missations_nondestructively(tree->at(S).missations, tree->at(P).missations);
    tree->at(P).missations.clear();
  }
  
  // 2. Swap roles of P & G
  swap(tree->at(P).mutations, tree->at(G).mutations);
  swap(tree->at(P).missations, tree->at(G).missations);

  // 3. Bubble up common missations above G (in *new* topology!  see below)
  CHECK(tree->at(G).missations.empty());  // because of tree->at(P).missations().clear() above
  if (interval_sets_intersect(tree->at(S).missations.intervals, tree->at(U).missations.intervals)) {
    factor_out_common_missations(tree->at(S).missations, tree->at(U).missations, tree->at(G).missations);
  }

  // Change topology
  // ---------------

  // Replace P by G as GG's child
  if (G == tree->root) {
    // Now P will be the root
    tree->root = P;
    tree->at(P).parent = k_no_node;
  } else {
    auto GG = tree->at(G).parent;
    auto GU = tree->at(GG).sibling_of(G);
    tree->at(GG).children = {P, GU};
    tree->at(P).parent = GG;
    CHECK_EQ(tree->at(GU).parent, GG);
  }

  // Reshuffle children of P & G
  tree->at(P).children = {X, G};
  CHECK_EQ(tree->at(X).parent, P);
  tree->at(G).parent = P;
  tree->at(G).children = {S, U};
  tree->at(S).parent = G;
  CHECK_EQ(tree->at(U).parent, G);

  // Update lambda_i and num_sites_missing_at_every_node for P and G
  lambda_i->at(P) = lambda_i->at(G);
  num_sites_missing_at_every_node->at(P) = num_sites_missing_at_every_node->at(G);
  
  lambda_i->at(G) = lambda_i->at(P) +
      calc_delta_lambda_across_missations(*evo, tree->ref_sequence, *ref_cum_Q_l, tree->at(G).missations);
  num_sites_missing_at_every_node->at(G) = num_sites_missing_at_every_node->at(P) +
      tree->at(G).missations.num_sites();
}

auto Tree_editing_session::flip() -> void {
  CHECK_NE(X, tree->root);
  auto P = tree->at(X).parent;
  CHECK(not tree->at(P).is_tip());
  CHECK_NE(P, tree->root);
  CHECK(tree->at(P).mutations.empty());

  auto G = tree->at(P).parent;
  CHECK_EQ(tree->at(P).t, tree->at(G).t);
  auto U = tree->at(G).sibling_of(P);

  auto S = tree->at(P).sibling_of(X);

  // Shuffle mutations and missations around G & P (see picture in tree_editing.h)
  // -----------------------------------------------------------------------------
  
  // 1. Push P's missations down to children
  if (not tree->at(P).missations.empty()) {
    tree->at(S).missations = merge_missations_nondestructively(tree->at(S).missations, tree->at(P).missations);
    tree->at(X).missations = merge_missations_nondestructively(tree->at(X).missations, tree->at(P).missations);
    tree->at(P).missations.clear();
  }

  // 2. Bubble up common missations above P (in *new* topology!  see below)
  if (interval_sets_intersect(tree->at(X).missations.intervals, tree->at(U).missations.intervals)) {
    factor_out_common_missations(tree->at(X).missations, tree->at(U).missations, tree->at(P).missations);
  }

  // Change topology
  // ---------------

  // Reshuffle children of P & G
  tree->at(G).children = {S, P};
  tree->at(S).parent = G;
  CHECK_EQ(tree->at(P).parent, G);

  tree->at(P).children = {X, U};
  CHECK_EQ(tree->at(X).parent, P);
  tree->at(U).parent = P;
  
  // Update lambda_i and num_sites_missing_at_every_node for P
  lambda_i->at(P) = lambda_i->at(G) +
      calc_delta_lambda_across_missations(*evo, tree->ref_sequence, *ref_cum_Q_l, tree->at(P).missations);
  num_sites_missing_at_every_node->at(P) = num_sites_missing_at_every_node->at(G) +
      tree->at(P).missations.num_sites();
}

auto Tree_editing_session::hop_down(Node_index SS) -> void {
  CHECK_NE(X, tree->root);
  auto P = tree->at(X).parent;
  CHECK(not tree->at(P).is_tip());

  CHECK_NE(SS, tree->root);
  auto U = tree->at(SS).parent;
  CHECK_EQ(tree->at(U).parent, P);
  auto UU = tree->at(U).sibling_of(SS);
  CHECK(tree->at(U).mutations.empty());

  do_hop_up(UU);
}

auto Tree_editing_session::end() -> void {
  CHECK(tree->at(X).mutations.empty());
  if (not deltas_nexus_to_X.empty()) {
    auto mut_t = 0.5*(tree->at(X).t + tree->at_parent_of(X).t);
    for (const auto& [l, delta] : deltas_nexus_to_X) {
      tree->at(X).mutations.push_back(Mutation{delta.from, l, delta.to, mut_t});
    }
  }
}

}  // namespace delphy
