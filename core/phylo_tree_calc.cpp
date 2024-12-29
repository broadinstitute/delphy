#include "phylo_tree_calc.h"

#include <absl/log/check.h>

#include "site_deltas.h"

namespace delphy {

auto count_mutations(const Phylo_tree& tree) -> int {
  auto result = 0;
  for (const auto& node : index_order_traversal(tree)) {
    if (node != tree.root) {
      result += std::ssize(tree.at(node).mutations);
    }
  }
  return result;
}

auto view_of_sequence_at(const Phylo_tree& tree, Phylo_tree_loc x) -> Sequence_overlay {
  auto scope = Local_arena_scope{};
  auto site_deltas = Site_deltas{};
  for (auto cur = x.branch; cur != k_no_node; cur = tree.at(cur).parent) {
    for (const auto& m : tree.at(cur).mutations | std::views::reverse) {
      if (m.t <= x.t) {
        push_front_site_deltas(m, site_deltas);
      }
    }
  }
  auto deltas = absl::flat_hash_map<Site_index, Real_seq_letter>{};
  for (const auto& [l, delta] : site_deltas) {
    CHECK_EQ(delta.from, tree.ref_sequence[l]);
    deltas.try_emplace(l, delta.to);
  }
  return Sequence_overlay{tree.ref_sequence, std::move(deltas)};
}

auto view_of_sequence_at(const Phylo_tree& tree, Node_index X) -> Sequence_overlay {
  return view_of_sequence_at(tree, tree.node_loc(X));
}

auto reconstruct_missing_sites_at(
    const Phylo_tree& tree,
    Node_index node)
    -> Scratch_interval_set {
  
  // Accumulate missations all the way up to (and including!) the root
  auto scratch1 = Scratch_interval_set{};
  auto scratch2 = Scratch_interval_set{};
  auto so_far = &scratch1;
  auto cur_scratch = &scratch2;
  for (auto cur = node; cur != k_no_node; cur = tree.at(cur).parent) {
    merge_interval_sets(*cur_scratch, *so_far, tree.at(cur).missations.intervals);
    std::swap(so_far, cur_scratch);
  }
  return *so_far;
}

auto is_site_missing_at(const Phylo_tree& tree, Node_index node, Site_index l) -> bool {
  for (auto cur = node; cur != k_no_node; cur = tree.at(cur).parent) {
    if (tree.at(cur).missations.contains(l)) {
      return true;
    }
  }
  return false;
}

auto calc_num_sites_missing_at_every_node(const Phylo_tree& tree) -> Node_vector<int> {
  auto result = Node_vector<int>(std::ssize(tree));

  for (const auto& node : pre_order_traversal(tree)) {
    auto count_at_parent = node == tree.root ? 0 : result[tree.at(node).parent];
    result[node] = count_at_parent + tree.at(node).missations.num_sites();
  }
  
  return result;
}

auto recalc_num_sites_missing_upstream(
    const Phylo_tree& tree,
    Node_index node,
    Node_index ancestor,
    Node_vector<int>& num_sites_missing_at_every_node)
    -> void {
  CHECK_NE(node, k_no_node);
  auto missing_at_cur = num_sites_missing_at_every_node[node];
  for (auto cur = node; cur != ancestor; cur = tree.at(cur).parent) {
    if (cur != tree.root) {
      auto missing_at_parent = missing_at_cur - tree.at(cur).missations.num_sites();
      num_sites_missing_at_every_node[tree.at(cur).parent] = missing_at_parent;
      missing_at_cur = missing_at_parent;
    }
  }
}

auto calc_state_frequencies_per_partition_of(
    const Real_sequence& seq,
    const Global_evo_model& evo)
    -> Partition_vector<Seq_vector<int>> {
  auto result = Partition_vector<Seq_vector<int>>(evo.num_partitions(), Seq_vector<int>{});
  for (auto l = Site_index{0}; l != std::ssize(seq); ++l) {
    auto partition = evo.partition_for_site[l];
    auto c = seq[l];
    ++result[partition][c];
  }
  return result;
}

auto calc_site_state_at(const Phylo_tree& tree, Phylo_tree_loc query_loc, Site_index l) -> Real_seq_letter {
  for (auto cur = query_loc.branch; cur != k_no_node; cur = tree.at(cur).parent) {
    for (const auto& m : tree.at(cur).mutations | std::views::reverse) {
      if (m.t > query_loc.t) { continue; }
      if (m.site == l) { return m.to; }  // First mutation on site l upstream of query_loc
    }
  }

  // There are no mutations on site l between the reference location (far above the root) and query_loc
  return tree.ref_sequence[l];
}

auto calc_T(const Phylo_tree& tree) -> double {
  auto T = 0.0;
  for (const auto& node : index_order_traversal(tree)) {
    if (node != tree.root) {
      T += tree.at(node).t - tree.at_parent_of(node).t;
    }
  }
  return T;
}

auto calc_T_l_a(const Phylo_tree& tree) -> std::vector<Seq_vector<double>> {
  // First calculate total branch length below every node
  auto T_below_node = Node_vector<double>(std::ssize(tree), 0.0);
  for (const auto& node : post_order_traversal(tree)) {
    auto T = 0.0;
    if (tree.at(node).is_inner_node()) {
      for (const auto& child : tree.at(node).children) {
        T += (tree.at(child).t - tree.at(node).t) + T_below_node[child];
      }
    }
    T_below_node[node] = T;
  }

  // Make a good initial guess for T_l_a as if there were *no* mutations on the tree
  auto T_l_a = Node_vector<Seq_vector<double>>(tree.num_sites(), Seq_vector{0.0, 0.0, 0.0, 0.0});
  auto T = T_below_node[tree.root];
  for (auto l = Site_index{0}; l != tree.num_sites(); ++l) {
    auto a = tree.ref_sequence[l];
    T_l_a[l][a] = T;
  }

  // Every mutation on the tree replaces the contribution from the total branch length downstream of the mutation
  // from `from` to `to`
  for (const auto& node : index_order_traversal(tree)) {
    for (const auto& m : tree.at(node).mutations) {
      auto T_below_mut = T_below_node[node] + (node == tree.root ? 0.0 : tree.at(node).t - m.t);
      T_l_a[m.site][m.from] -= T_below_mut;
      T_l_a[m.site][m.to] += T_below_mut;
    }

    for (const auto& [mi_site, mi_from] : tree.at(node).missations.slow_elements(tree.ref_sequence)) {
      auto T_below_miss = T_below_node[node] +
          (node == tree.root ? 0.0 : tree.at(node).t - tree.at_parent_of(node).t);
      T_l_a[mi_site][mi_from] -= T_below_miss;
    }
  }

  return T_l_a;
}

auto calc_Ttwiddle_l(const Phylo_tree& tree, const Global_evo_model& evo) -> std::vector<double> {
  // TTwiddle^(l) = sum_a q^(l)_a T^(l)_a, but it's more efficient to accumulate it directly
  
  // First calculate total branch length below every node
  auto T_below_node = Node_vector<double>(std::ssize(tree), 0.0);
  for (const auto& node : post_order_traversal(tree)) {
    auto T = 0.0;
    if (tree.at(node).is_inner_node()) {
      for (const auto& child : tree.at(node).children) {
        T += (tree.at(child).t - tree.at(node).t) + T_below_node[child];
      }
    }
    T_below_node[node] = T;
  }

  // Make a good initial guess for Ttwiddle^(l) as if there were *no* mutations on the tree
  auto Ttwiddle_l = Node_vector<double>(tree.num_sites(), 0.0);
  auto T = T_below_node[tree.root];
  for (auto l = Site_index{0}; l != tree.num_sites(); ++l) {
    auto a = tree.ref_sequence[l];
    Ttwiddle_l[l] = evo.q_l_a(l, a) * T;
  }

  // Every mutation on the tree replaces the contribution from the total branch length downstream of the mutation
  // from `from` to `to`
  for (const auto& node : index_order_traversal(tree)) {
    for (const auto& m : tree.at(node).mutations) {
      auto T_below_mut = T_below_node[node] + (node == tree.root ? 0.0 : tree.at(node).t - m.t);
      Ttwiddle_l[m.site] -= evo.q_l_a(m.site, m.from) * T_below_mut;
      Ttwiddle_l[m.site] += evo.q_l_a(m.site, m.to) * T_below_mut;
    }

    for (const auto& [mi_site, mi_from] : tree.at(node).missations.slow_elements(tree.ref_sequence)) {
      auto T_below_miss = T_below_node[node] +
          (node == tree.root ? 0.0 : tree.at(node).t - tree.at_parent_of(node).t);
      Ttwiddle_l[mi_site] -= evo.q_l_a(mi_site, mi_from) * T_below_miss;
    }
  }

  return Ttwiddle_l;
}

auto calc_Ttwiddle_a(const Phylo_tree& tree, std::span<const double> nu_l) -> Seq_vector<double> {
  // TTwiddle_a = \sum_l nu^(l) T^(l)_a
  auto Ttwiddle_a = Seq_vector{0.0, 0.0, 0.0, 0.0};
  auto dTtwiddle_a_dt = calc_dTtwiddle_a_dt_for_sequence(tree.ref_sequence, nu_l);

  for (const auto& [node, children_so_far] : traversal(tree)) {
    if (children_so_far == 0) {
      // Entering `node`
      // At this point, dTtwiddle_a_dt reflects the state at `parent`.
      auto parent = tree.at(node).parent;
      
      // Update dTtwiddle_a_dt to reflect the state at `node`
      // TODO: Treat gaps in one go
      for (const auto& [mi_site, mi_from] : tree.at(node).missations.slow_elements(tree.ref_sequence)) {
        dTtwiddle_a_dt[mi_from] -= nu_l[mi_site];
      }
      for (const auto& m : tree.at(node).mutations) {
        dTtwiddle_a_dt[m.from] -= nu_l[m.site];
        dTtwiddle_a_dt[m.to] += nu_l[m.site];
      }

      if (node != tree.root) {
        // Accumulate portion of TTwiddle_a from the branch that ends in `node`
        auto branch_length = (tree.at(node).t - tree.at(parent).t);
        for (auto a : k_all_real_seq_letters) {
          Ttwiddle_a[a] += dTtwiddle_a_dt[a] * branch_length;
        }
        for (const auto& m : tree.at(node).mutations) {
          Ttwiddle_a[m.to] -= nu_l[m.site] * (m.t - tree.at(parent).t);
          Ttwiddle_a[m.from] += nu_l[m.site] * (m.t - tree.at(parent).t);
        }
      }
    }
    
    if (children_so_far == std::ssize(tree.at(node).children)) {
      // Exiting `node`
      // At this point, dTtwiddle_a_dt reflects the state at `node`
      
      // Update dTtwiddle_a_dt to reflect the state at `parent`
      for (const auto& m : tree.at(node).mutations) {
        dTtwiddle_a_dt[m.to] -= nu_l[m.site];
        dTtwiddle_a_dt[m.from] += nu_l[m.site];
      }

      // TODO: Treat gaps as one
      for (const auto& [mi_site, mi_from] : tree.at(node).missations.slow_elements(tree.ref_sequence)) {
        dTtwiddle_a_dt[mi_from] += nu_l[mi_site];
      }
    }
  }
  
  return Ttwiddle_a;
}

auto calc_Ttwiddle_beta_a(
    const Phylo_tree& tree,
    const Global_evo_model& evo)
    -> Partition_vector<Seq_vector<double>> {
  
  // TTwiddle^\beta_a = \sum_{l in beta} nu^(l) T^(l)_a
  auto Ttwiddle_beta_a = Partition_vector<Seq_vector<double>>(evo.num_partitions(), Seq_vector<double>{0.0});
  auto ntwiddle_beta_a = Partition_vector<Seq_vector<double>>(evo.num_partitions(), Seq_vector<double>{0.0});
  for (auto l = Site_index{0}; l != tree.num_sites(); ++l) {
    auto beta = evo.partition_for_site[l];
    auto s_r = tree.ref_sequence[l];
    ntwiddle_beta_a[beta][s_r] += evo.nu_l[l];
  }

  for (const auto& [node, children_so_far] : traversal(tree)) {
    if (children_so_far == 0) {
      // Entering `node`
      // At this point, ntwiddle_beta_a reflects the state at `parent`.
      auto parent = tree.at(node).parent;
      
      // Update ntwiddle_beta_a to reflect the state at `node`
      // TODO: Treat gaps in one go
      for (const auto& [mi_site, mi_from] : tree.at(node).missations.slow_elements(tree.ref_sequence)) {
        auto beta = evo.partition_for_site[mi_site];
        ntwiddle_beta_a[beta][mi_from] -= evo.nu_l[mi_site];
      }
      for (const auto& m : tree.at(node).mutations) {
        auto beta = evo.partition_for_site[m.site];
        ntwiddle_beta_a[beta][m.from] -= evo.nu_l[m.site];
        ntwiddle_beta_a[beta][m.to] += evo.nu_l[m.site];
      }

      if (node != tree.root) {
        // Accumulate portion of TTwiddle^beta_a from the branch that ends in `node`
        auto branch_length = (tree.at(node).t - tree.at(parent).t);
        for (auto beta = Partition_index{0}; beta != evo.num_partitions(); ++beta) {
          for (auto a : k_all_real_seq_letters) {
            Ttwiddle_beta_a[beta][a] += ntwiddle_beta_a[beta][a] * branch_length;
          }
        }
        for (const auto& m : tree.at(node).mutations | std::views::reverse) {
          // Remove contribution to Ttwiddle^beta_a from segment [t_P,m.t) when state is `m.to`
          // and replace it with the analogous contribution when state is `m.from`
          auto beta = evo.partition_for_site[m.site];
          Ttwiddle_beta_a[beta][m.to] -= evo.nu_l[m.site] * (m.t - tree.at(parent).t);
          Ttwiddle_beta_a[beta][m.from] += evo.nu_l[m.site] * (m.t - tree.at(parent).t);
        }
      }
    }
    
    if (children_so_far == std::ssize(tree.at(node).children)) {
      // Exiting `node`
      // At this point, ntwiddle_beta_a reflects the state at `node`
      
      // Update ntwiddle_beta_a to reflect the state at `parent`
      for (const auto& m : tree.at(node).mutations) {
        auto beta = evo.partition_for_site[m.site];
        ntwiddle_beta_a[beta][m.to] -= evo.nu_l[m.site];
        ntwiddle_beta_a[beta][m.from] += evo.nu_l[m.site];
      }
      // TODO: Treat gaps in one go
      for (const auto& [mi_site, mi_from] : tree.at(node).missations.slow_elements(tree.ref_sequence)) {
        auto beta = evo.partition_for_site[mi_site];
        ntwiddle_beta_a[beta][mi_from] += evo.nu_l[mi_site];
      }
    }
  }
  
  return Ttwiddle_beta_a;
}

auto calc_dTtwiddle_a_dt_for_sequence(const Real_sequence& seq, std::span<const double> nu_l) -> Seq_vector<double> {
  auto dTtwiddle_a_dt = Seq_vector{0.0, 0.0, 0.0, 0.0};
  for (auto l = Site_index{0}; l != std::ssize(seq); ++l) {
    dTtwiddle_a_dt[seq[l]] += nu_l[l];
  }
  return dTtwiddle_a_dt;
}

auto calc_cum_Q_l_for_sequence(const Real_sequence& seq, const Global_evo_model& evo) -> std::vector<double> {
  auto cum_Q_l = std::vector<double>(seq.size()+1, 0.0);  // Note "+1".  See header comment for details
  auto so_far = 0.0;
  cum_Q_l[0] = 0.0;
  for (auto l = Site_index{0}; l != std::ssize(seq); ++l) {
    so_far += evo.mu_l(l) * evo.nu_l[l] * evo.q_l_a(l, seq[l]);
    cum_Q_l[l+1] = so_far;
  }
  return cum_Q_l;
}

auto calc_lambda_for_sequence(
    const Real_sequence& seq,
    const Global_evo_model& evo)
    -> double {
  auto lambda = 0.0;
  for (auto l = Site_index{0}; l != std::ssize(seq); ++l) {
    lambda += evo.mu_l(l) * evo.nu_l[l] * evo.q_l_a(l, seq[l]);
  }
  return lambda;
}

auto calc_lambda_at_node(const Phylo_tree& tree, Node_index node, const Global_evo_model& evo) -> double {
  auto ref_cum_Q_l = calc_cum_Q_l_for_sequence(tree.ref_sequence, evo);
  return calc_lambda_at_node(tree, node, evo, ref_cum_Q_l);
}

auto calc_lambda_at_node(
    const Phylo_tree& tree,
    Node_index node,
    const Global_evo_model& evo,
    const std::vector<double>& ref_cum_Q_l)
    -> double {
  auto result = ref_cum_Q_l.back();  // == lambda_for_ref_sequence
  for (auto cur = node; cur != k_no_node; cur = tree.at(cur).parent) {
    result += calc_delta_lambda_across_branch(
        evo, tree.ref_sequence, ref_cum_Q_l, tree.at(cur).mutations, tree.at(cur).missations);
  }
  return result;
}

auto calc_lambda_i(
    const Phylo_tree& tree,
    const Global_evo_model& evo,
    const std::vector<double>& ref_cum_Q_l)
    -> Node_vector<double> {
  
  auto result = Node_vector<double>(std::ssize(tree));
  auto lambda_ref = ref_cum_Q_l.back();  // == lambda_for_ref_sequence

  for (const auto& node : pre_order_traversal(tree)) {
    auto lambda_parent = node == tree.root ? lambda_ref : result[tree.at(node).parent];
    result[node] = lambda_parent + calc_delta_lambda_across_branch(
        evo, tree.ref_sequence, ref_cum_Q_l, tree.at(node).mutations, tree.at(node).missations);
  }
  
  return result;
}

auto recalc_lambda_i_upstream(
    const Phylo_tree& tree,
    Node_index node,
    Node_index ancestor,
    const Global_evo_model& evo,
    Node_vector<double>& lambda_i,
    const std::vector<double>& ref_cum_Q_l)
    -> void {
  CHECK_NE(node, k_no_node);
  auto lambda_end = lambda_i[node];
  for (auto cur = node; cur != ancestor; cur = tree.at(cur).parent) {
    if (cur != tree.root) {
      auto lambda_parent = lambda_end - calc_delta_lambda_across_branch(
          evo, tree.ref_sequence, ref_cum_Q_l, tree.at(cur).mutations, tree.at(cur).missations);
      lambda_i[tree.at(cur).parent] = lambda_parent;
      lambda_end = lambda_parent;
    }
  }
}

auto calc_log_root_prior(
    const Phylo_tree& tree,
    const Global_evo_model& evo)
    -> double {
  auto state_frequencies_of_ref_sequence_per_partition =
      calc_state_frequencies_per_partition_of(tree.ref_sequence, evo);
  return calc_log_root_prior(tree, evo, state_frequencies_of_ref_sequence_per_partition);
}

auto calc_log_root_prior(
    const Phylo_tree& tree,
    const Global_evo_model& evo,
    const Partition_vector<Seq_vector<int>>& state_frequencies_of_ref_sequence_per_partition)
    -> double {

  auto state_frequencies_at_root_per_partition = state_frequencies_of_ref_sequence_per_partition;
  for (const auto& m : tree.at_root().mutations) {
    auto partition = evo.partition_for_site[m.site];
    --state_frequencies_at_root_per_partition[partition][m.from];
    ++state_frequencies_at_root_per_partition[partition][m.to];
  }
  // TODO: Treat gaps in one go
  for (const auto& [mi_site, mi_from] : tree.at_root().missations.slow_elements(tree.ref_sequence)) {
    auto partition = evo.partition_for_site[mi_site];
    --state_frequencies_at_root_per_partition[partition][mi_from];
  }

  auto result = 0.0;
  for (auto partition = Partition_index{0}; partition != evo.num_partitions(); ++partition) {
    const auto& pi_a = evo.partition_evo_model[partition].pi_a;
    CHECK(is_stochastic_vector(pi_a)) << partition << " " << pi_a;
    for (auto a : k_all_real_seq_letters) {
      if (pi_a[a] != 0.0) {
        result += state_frequencies_at_root_per_partition[partition][a] * std::log(pi_a[a]);
      } else if (state_frequencies_at_root_per_partition[partition][a] != 0) {
        return -std::numeric_limits<double>::infinity();
      }
    }
  }
  return result;
}

auto calc_log_G_below_root(
    const Phylo_tree& tree,
    const Global_evo_model& evo)
    -> double {
  
  return calc_log_G_below_root(tree, evo, calc_lambda_i(tree, evo, calc_cum_Q_l_for_sequence(tree.ref_sequence, evo)),
                               calc_state_frequencies_per_partition_of(tree.ref_sequence, evo));
}

auto calc_log_G_below_root(
    const Phylo_tree& tree,
    const Global_evo_model& evo,
    const Node_vector<double>& lambda_i,
    const Partition_vector<Seq_vector<int>>& state_frequencies_of_ref_sequence_per_partition)
    -> double {

  if (estd::is_debug_enabled) {
    for (auto i = Node_index{0}; i != std::ssize(tree); ++i) {
      CHECK_GE(lambda_i[i], 0.0) << i;
    }
    for (auto beta = Partition_index{0}; beta != evo.num_partitions(); ++beta) {
      CHECK_GE(evo.partition_evo_model[beta].mu, 0.0) << beta;
      CHECK(is_transition_rate_matrix(evo.partition_evo_model[beta].q_ab))
          << beta << " " << evo.partition_evo_model[beta].q_ab;
    }
    for (auto l = Site_index{0}; l != tree.num_sites(); ++l) {
      CHECK_GE(evo.nu_l[l], 0.0) << l;
    }
  }
  
  auto result = 0.0;
  for (const auto& node : index_order_traversal(tree)) {
    if (node != tree.root) {
      result += calc_branch_log_G(tree, node, lambda_i[node], evo, state_frequencies_of_ref_sequence_per_partition);
    }
  }
  return result;
}

auto calc_branch_log_G(
    const Phylo_tree& tree,
    Branch_index X,
    double lambda_X,
    const Global_evo_model& evo,
    const Partition_vector<Seq_vector<int>>& state_frequencies_of_ref_sequence_per_partition)
    -> double {

  if (X == tree.root) {
    return calc_log_root_prior(tree, evo, state_frequencies_of_ref_sequence_per_partition);
  } else {
    return calc_branch_log_G(tree.at_parent_of(X).t, tree.at(X).t, lambda_X, evo, tree.at(X).mutations);
  }
}

auto calc_path_log_G(
    const Phylo_tree& tree,
    Node_index A,
    Node_index B,
    const Global_evo_model& evo,
    const Node_vector<double>& lambda_i,
    const Partition_vector<Seq_vector<int>>& state_frequencies_of_ref_sequence_per_partition)
    -> double {
  DCHECK(descends_from(tree, B, A));

  auto result = 0.0;
  for (auto cur = B; cur != A; cur = tree.at(cur).parent) {
    result += calc_branch_log_G(tree, cur, lambda_i[cur], evo, state_frequencies_of_ref_sequence_per_partition);
  }
  return result;
}

auto calc_num_muts(const Phylo_tree& tree) -> int {
  auto num_muts = 0;
  for (const auto& node : index_order_traversal(tree)) {
    if (node != tree.root) {  // "Mutations" above root are just deltas from reference sequence
      num_muts += std::ssize(tree.at(node).mutations);
    }
  }
  return num_muts;
}

auto calc_num_muts_ab(const Phylo_tree& tree) -> Seq_matrix<int> {
  auto result = Seq_matrix<int>{};
  for (const auto& node : index_order_traversal(tree)) {
    if (node != tree.root) {  // "Mutations" above root are just deltas from reference sequence
      for (const auto& m : tree.at(node).mutations) {
        ++result[m.from][m.to];
      }
    }
  }
  return result;
}

auto calc_num_muts_beta_ab(const Phylo_tree& tree, const Global_evo_model& evo) -> Partition_vector<Seq_matrix<int>> {
  auto result = Partition_vector<Seq_matrix<int>>(evo.num_partitions(), Seq_matrix<int>{});
  for (const auto& node : index_order_traversal(tree)) {
    if (node != tree.root) {  // "Mutations" above root are just deltas from reference sequence
      for (const auto& m : tree.at(node).mutations) {
        auto beta = evo.partition_for_site[m.site];
        ++result[beta][m.from][m.to];
      }
    }
  }
  return result;
}

auto calc_num_muts_l(const Phylo_tree& tree) -> Node_vector<int> {
  auto result = Node_vector<int>(tree.num_sites(), 0);
  for (const auto& node : index_order_traversal(tree)) {
    if (node != tree.root) {  // "Mutations" above root are just deltas from reference sequence
      for (const auto& m : tree.at(node).mutations) {
        ++result[m.site];
      }
    }
  }
  return result;
}

auto calc_num_muts_l_ab(const Phylo_tree& tree) -> Node_vector<Seq_matrix<int>> {
  auto result = Node_vector<Seq_matrix<int>>(tree.num_sites(), Seq_matrix<int>{});
  for (const auto& node : index_order_traversal(tree)) {
    if (node != tree.root) {  // "Mutations" above root are just deltas from reference sequence
      for (const auto& m : tree.at(node).mutations) {
        ++result[m.site][m.from][m.to];
      }
    }
  }
  return result;
}

auto calc_max_tip_time(const Phylo_tree& tree) -> double {
  auto t_max = -std::numeric_limits<double>::infinity();
  for (const auto& node : index_order_traversal(tree)) {
    if (tree.at(node).is_tip() && (tree.at(node).t_max > t_max)) {
      t_max = tree.at(node).t_max;
    }
  }
  return t_max;
}

auto assert_tip_sequences_compatible_with_original_ones(
    const Phylo_tree& tree,
    const Node_vector<Sequence_overlay>& original_sequences,
    const Node_vector<Interval_set<>>& original_missing_sites_all)
    -> void {
  
  if (estd::is_debug_enabled) {
    for (const auto& node : index_order_traversal(tree)) {
      if (tree.at(node).is_tip()) {
        auto scope = Local_arena_scope{};
        auto current_missing_sites = reconstruct_missing_sites_at(tree, node);
        const auto& original_missing_sites = original_missing_sites_all[node];
        CHECK(std::ranges::equal(current_missing_sites, original_missing_sites))
            << absl::StreamFormat(
                "Tip %d has a different set of missing sites ({%s}) than it originally did ({%s})",
            node,
            absl::FormatStreamed(current_missing_sites),
            absl::FormatStreamed(original_missing_sites));
        
        auto current_sequence = view_of_sequence_at(tree, node).materialize();
        const auto& original_sequence = original_sequences[node];
        for (auto l = Site_index{0}; l != tree.num_sites(); ++l) {
          if (not current_missing_sites.contains(l)) {
            CHECK(current_sequence[l] == original_sequence[l]) << absl::StreamFormat(
                "Tip %d now has a %c at site %d, but used to have a %c",
                node, to_char(current_sequence[l]), l+1, to_char(original_sequence[l]));
          }
        }
      }
    }
  }
}

}  // namespace delphy
