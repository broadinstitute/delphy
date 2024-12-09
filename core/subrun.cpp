#include "subrun.h"

#include "distributions.h"
#include "phylo_tree_calc.h"
#include "spr_move.h"
#include "spr_study.h"

namespace delphy {

Subrun::Subrun(absl::BitGenRef bitgen, Phylo_tree tree, bool includes_run_root, Global_evo_model evo)
    : bitgen_{bitgen},
      tree_{std::move(tree)},
      includes_run_root_{includes_run_root},
      scratch_{},
      evo_{std::move(evo)} {
  
  assert_phylo_tree_integrity(tree_, true);  // check even in release builds
}

auto Subrun::recalc_derived_quantities() const -> void {
  // Derived quantities
  state_frequencies_of_ref_sequence_per_partition_ =
      calc_cur_state_frequencies_of_ref_sequence_per_partition();
  ref_cum_Q_l_ = calc_cur_ref_cum_Q_l();
  lambda_i_ = calc_cur_lambda_i(ref_cum_Q_l_);
  log_G_ = calc_cur_log_G(lambda_i_, state_frequencies_of_ref_sequence_per_partition_);
  num_sites_missing_at_every_node_ = calc_cur_num_sites_missing_at_every_node();
  log_augmented_coalescent_prior_ = calc_cur_log_augmented_coalescent_prior();
}

auto Subrun::check_derived_quantities() const -> void {
  if (estd::is_debug_enabled) {
    // We maintain these incrementally except when we force a recalculation
    { auto expected = calc_cur_log_augmented_coalescent_prior();
      CHECK(std::abs(log_augmented_coalescent_prior_ - expected) < 1e-5)
      << log_augmented_coalescent_prior_ << " != " << expected; }
    { auto expected = calc_cur_state_frequencies_of_ref_sequence_per_partition();
      for (auto beta = Partition_index{0}; beta != evo_.num_partitions(); ++beta) {
        CHECK_EQ(state_frequencies_of_ref_sequence_per_partition_[beta], expected[beta])
            << "Partition " << beta; } }
    { auto expected = calc_cur_ref_cum_Q_l();
      for (auto l = Site_index{0}; l <= tree_.num_sites(); ++l) {
        CHECK_LT(std::abs(ref_cum_Q_l_[l] - expected[l]), 1e-8)
            << "Site " << l << ", "
            << ref_cum_Q_l_[l] << " != "
            << expected[l]; } }
    { auto expected = calc_cur_lambda_i(ref_cum_Q_l_);
      for (const auto& node : index_order_traversal(tree_)) {
        CHECK((std::abs(lambda_i_[node] - expected[node])
               / static_cast<double>(tree_.num_sites())) < 1e-8)
            << "Node " << node << ", "
            << lambda_i_[node] << " != "
            << expected[node]; } }
    { auto expected = calc_cur_log_G(lambda_i_, state_frequencies_of_ref_sequence_per_partition_);
      CHECK(std::abs(log_G_ - expected) < 1e-6)
          << log_G_ << " != " << expected; }
    CHECK(num_sites_missing_at_every_node_ == calc_cur_num_sites_missing_at_every_node());
  }
}

auto Subrun::calc_cur_log_G(
    const Node_vector<double>& lambda_i,
    const Partition_vector<Seq_vector<int>>& state_frequencies_of_ref_sequence_per_partition) const
    -> double {
  auto result = 0.0;
  if (includes_run_root_) {
    result += calc_log_root_prior(tree_, evo_, state_frequencies_of_ref_sequence_per_partition);
  }
  result += calc_log_G_below_root(tree_, evo_, lambda_i, state_frequencies_of_ref_sequence_per_partition);
  return result;
}

auto Subrun::calc_cur_state_frequencies_of_ref_sequence_per_partition() const -> Partition_vector<Seq_vector<int>> {
  return calc_state_frequencies_per_partition_of(tree_.ref_sequence, evo_);
}

auto Subrun::calc_cur_lambda_i(const std::vector<double>& ref_cum_Q_l) const -> Node_vector<double> {
  return calc_lambda_i(tree_, evo_, ref_cum_Q_l);
}

auto Subrun::calc_cur_ref_cum_Q_l() const -> std::vector<double> {
  return calc_cum_Q_l_for_sequence(tree_.ref_sequence, evo_);
}

auto Subrun::calc_cur_num_sites_missing_at_every_node() const -> Node_vector<int> {
  return calc_num_sites_missing_at_every_node(tree_);
}

auto Subrun::calc_cur_log_augmented_coalescent_prior() const -> double {
  if (coalescent_prior_part_ == nullptr) {
    return -1.0;
  } else {
    return coalescent_prior_part_->calc_partial_log_prior();
  }
}

auto Subrun::assert_cur_tip_sequences_compatible_with_original_ones() -> void {
  assert_tip_sequences_compatible_with_original_ones(tree_, original_sequences_, original_missing_sites_all_);
}

auto Subrun::mcmc_sub_iteration() -> void {
  // All moves can allocate quickly in the Scratch_space scratch_, which will all be freed at the end of the move
  
  validate_derived_quantities();

  if (only_displacing_inner_nodes_) {
    inner_node_displace_move();
  } else {
    auto total_weight = 15.0 + 15.0;
    if (topology_moves_enabled()) { total_weight += 1.0 + 1.0; }
    auto r = absl::Uniform(bitgen_, 0.0, total_weight);
    if (r < 7.5) { inner_node_displace_move(); }     // Weight 7.5
    else if (r < 15.0) { tip_displace_move(); }      // Weight 7.5
    else if (r < 30.0) { branch_reform_move(); }     // Weight 15.0
    else if (topology_moves_enabled()) {
      if (r < 31.0) { subtree_slide_move(); }        // Weight 1.0
      else { spr1_move(); }                          // Weight 1.0
    }
  }

  check_derived_quantities();

  scratch_.reset();
}

auto Subrun::pick_random_node() -> Node_index {
  DCHECK_GT(std::ssize(tree_), 0);
  return absl::Uniform<Node_index>(absl::IntervalClosedOpen, bitgen_, 0, std::ssize(tree_));
}

auto Subrun::pick_random_inner_node() -> Node_index {
  DCHECK_GT(std::ssize(tree_), 0);
  while (true) {
    auto result = pick_random_node();
    if (tree_.at(result).is_inner_node()) {
      return result;
    }
  }
}

auto Subrun::pick_random_tip() -> Node_index {
  DCHECK_GT(std::ssize(tree_), 0);
  while (true) {
    auto result = pick_random_node();
    if (tree_.at(result).is_tip()) {
      return result;
    }
  }
}

auto Subrun::inner_node_displace_move() -> void {
  if (std::ssize(tree_) < 1) return;

  // We only move the root of the subtree if it's the root of the whole tree
  auto node = pick_random_inner_node();
  if (node == tree_.root && not includes_run_root_) { return; }

  auto t_min = -std::numeric_limits<double>::infinity();
  if (node != tree_.root) {
    t_min = tree_.at_parent_of(node).t;
    for (const auto& m : tree_.at(node).mutations) {
      t_min = std::max(t_min, m.t);
    }
  }

  auto t_max = +std::numeric_limits<double>::infinity();
  for (const auto& child : tree_.at(node).children) {
    t_max = std::min(t_max, tree_.at(child).t);
    for (const auto& m : tree_.at(child).mutations) {
      t_max = std::min(t_max, m.t);
    }
  }

  auto lambda_at_node = lambda_i_[node];
  auto d_logG_dt = 0.0;
  if (node != tree_.root) {
    // parent branch lengthens
    d_logG_dt += -lambda_at_node;
  }
  for (const auto& child : tree_.at(node).children) {
    // child branches shorten (branch start of child is like branch end of node, adjusted for child missations)
    auto lambda_just_below_node = lambda_at_node + calc_delta_lambda_across_missations(
        evo_, tree_.ref_sequence, ref_cum_Q_l_, tree_.at(child).missations);
    d_logG_dt -= -lambda_just_below_node;
  }

   // Pick tree_.at(node).t as ~ exp(d_logG_dt * node.t), bounded by t_min and t_max.
   auto old_node_t = tree_.at(node).t;
  auto log_alpha_old_to_new_over_new_to_old = 0.0;
  auto new_node_t = old_node_t;
  if (node == tree_.root) {
    auto delta_scale = (1 / lambda_i_.at(node)) / 2;  // 95% of the time, won't hop past more than 1 mutation
    new_node_t = old_node_t + absl::Gaussian(bitgen_, 0.0, delta_scale);
    if (new_node_t < t_min || new_node_t > t_max) { return; }
    log_alpha_old_to_new_over_new_to_old = 0.0;
  } else {
    auto t_dist = Bounded_exponential_distribution{d_logG_dt, t_min, t_max};
    new_node_t = t_dist(bitgen_);
    
    // The proposal probability being proportional to the likelihood, the MH ratio will be 1.0 (before considering the coalescent prior)
    log_alpha_old_to_new_over_new_to_old = d_logG_dt * (new_node_t - old_node_t);
  }

  // Forbid zero-length branches or mutation times that coincide with node times
  if (new_node_t == t_min || new_node_t == t_max) { return; }

  CHECK(t_min <= new_node_t) << lambda_at_node << ","
                             << d_logG_dt << "," << t_min << "," << t_max << "," << new_node_t;
  CHECK(new_node_t <= t_max);

  // See if coalescent prior agrees with us
  auto delta_log_G = d_logG_dt * (new_node_t - old_node_t);
  auto delta_log_prior =
      coalescent_prior_part_->calc_delta_partial_log_prior_after_displace_coalescence(old_node_t, new_node_t);
  auto log_mh = delta_log_G + delta_log_prior - log_alpha_old_to_new_over_new_to_old;
  if (log_mh >= 0.0 || absl::Uniform(bitgen_, 0.0, 1.0) < std::exp(log_mh)) {
    // Accept
    coalescent_prior_part_->coalescence_displaced(old_node_t, new_node_t);
    tree_.at(node).t = new_node_t;

    auto dt = new_node_t - old_node_t;
    log_G_ += d_logG_dt * dt;
    log_augmented_coalescent_prior_ += delta_log_prior;
  }
}

auto Subrun::tip_displace_move() -> void {
  if (std::ssize(tree_) < 1) return;

  auto node = pick_random_tip();
  CHECK_NE(node, tree_.root);
  if (tree_.at(node).t_min == tree_.at(node).t_max) { return; }  // No uncertainty!

  auto t_min = std::max(double{tree_.at(node).t_min}, tree_.at_parent_of(node).t);
  for (const auto& m : tree_.at(node).mutations) {
    t_min = std::max(t_min, m.t);
  }

  auto t_max = double{tree_.at(node).t_max};

  auto lambda_at_node = lambda_i_[node];
  auto d_logG_dt = 0.0;
  
  // parent branch lengthens
  d_logG_dt += -lambda_at_node;

  // Pick tree_.at(node).t as ~ exp(d_logG_dt * node.t), bounded by t_min and t_max.
  auto old_node_t = tree_.at(node).t;
  auto log_alpha_old_to_new_over_new_to_old = 0.0;
  
  auto t_dist = Bounded_exponential_distribution{d_logG_dt, t_min, t_max};
  auto new_node_t = t_dist(bitgen_);
  
  // The proposal probability being proportional to the likelihood, the MH ratio will be 1.0 (before considering the coalescent prior)
  log_alpha_old_to_new_over_new_to_old = d_logG_dt * (new_node_t - old_node_t);

  // Forbid zero-length branches or mutation times that coincide with node times
  if (new_node_t == t_min || new_node_t == t_max) { return; }

  CHECK(t_min <= new_node_t) << lambda_at_node << ","
                             << d_logG_dt << "," << t_min << "," << t_max << "," << new_node_t;
  CHECK(new_node_t <= t_max);

  // See if coalescent prior agrees with us
  auto delta_log_G = d_logG_dt * (new_node_t - old_node_t);
  auto delta_log_prior =
      coalescent_prior_part_->calc_delta_partial_log_prior_after_displace_tip(old_node_t, new_node_t);
  auto log_mh = delta_log_G + 0.0*delta_log_prior - log_alpha_old_to_new_over_new_to_old;
  if (log_mh >= 0.0 || absl::Uniform(bitgen_, 0.0, 1.0) < std::exp(log_mh)) {
    // Accept
    coalescent_prior_part_->tip_displaced(old_node_t, new_node_t);
    tree_.at(node).t = new_node_t;

    auto dt = new_node_t - old_node_t;
    log_G_ += d_logG_dt * dt;
    log_augmented_coalescent_prior_ += delta_log_prior;
  }
}

auto Subrun::branch_reform_move() -> void {
  if (std::ssize(tree_) < 3) return;

  auto X = pick_random_node();
  if (X == tree_.root) { return; }  // The branch above the root never has (meaningful) mutations

  auto P = tree_.at(X).parent;
  auto S = tree_.at(P).sibling_of(X);
  auto t_X = tree_.at(X).t;
  auto t_P = tree_.at(P).t;

  if (P == tree_.root) {
    // Allow the mutations on the P-X and P-S branch to dance;
    // the logic below would artificially only reform the mutations on the P-X branch
    spr_move_core(X, {S, t_P}, 1.0);
    //return;  // Fall through instead
  }
  
  // Proposal
  auto& old_mutations = tree_.at(X).mutations;
  auto new_mutations = randomize_branch_mutation_times(tree_, X, bitgen_, scratch_);

  // Accept/reject
  auto delta_log_G = 
      calc_branch_log_G(t_P, t_X, lambda_i_.at(X), evo_, new_mutations)
      - calc_branch_log_G(t_P, t_X, lambda_i_.at(X), evo_, old_mutations);
  auto log_mh = delta_log_G;
  auto accept = log_mh >= 0.0 || absl::Uniform(bitgen_, 0.0, 1.0) < std::exp(log_mh);
  
  if (accept) {
    tree_.at(X).mutations.assign(new_mutations.begin(), new_mutations.end());
    log_G_ += delta_log_G;
  }
}

// Enumerate all branches that straddle time tree and are at or below the branch ending at P,
// not at or below the branch ending at X.  The output iterator accumulates the Node_index's
// of the endpoint of the relevant branches.
template<typename Output_iterator>
static auto enumerate_descendant_branches_straddling(
    const Phylo_tree& tree,
    Node_index P,
    double t,
    Node_index X,
    Output_iterator out)
    -> Output_iterator {
  
  DCHECK(P == tree.root || tree.at_parent_of(P).t <= t);

  if (P == X) {
    // Skip the branch ending at X and all of its descendants
    return out;
  }

  if (t <= tree.at(P).t) {
    *out = P;
    ++out;
  } else if (tree.at(P).is_inner_node()) {
    for (const auto& child : tree.at(P).children) {
      out = enumerate_descendant_branches_straddling(tree, child, t, X, out);
    }
  }
  return out;
}

auto Subrun::subtree_slide_move() -> void {
  if (std::ssize(tree_) < 2) return;

  // Comments below purposely track those in SubtreeSlideOperator#doOperation() in BEAST

  // 1. choose a random node avoiding root
  auto X = pick_random_node();
  if (X == tree_.root) { return; }

  auto P = tree_.at(X).parent;
  auto S = tree_.at(P).sibling_of(X);  // "CiP"

  // 2. choose a delta to move
  auto delta_scale = (1 / lambda_i_.at(X)) / 2;  // 95% of the time, won't hop past more than 1 mutation
  auto delta_t = absl::Gaussian(bitgen_, 0.0, delta_scale);
  auto old_P_t = tree_.at(P).t;
  auto new_P_t = tree_.at(P).t + delta_t;

  // 3. if the move is up
  if (delta_t < 0.0) {

    // 3.1 if the topology will change
    if (P != tree_.root && new_P_t < tree_.at_parent_of(P).t) {
      // find new parent
      auto GG = tree_.at(P).parent;  // new grandparent after move
      auto SS = P;                  // new sibling after move
      while (new_P_t < tree_.at(GG).t) {
        SS = GG;
        GG = tree_.at(GG).parent;
        if (GG == k_no_node) {
          break;
        }
      }
      CHECK((GG == k_no_node || tree_.at(GG).t <= new_P_t) && new_P_t <= tree_.at(SS).t);

      // 3.1.1 if creating a new root... => tree rewired in spr_move_core
      // 3.1.2 no new root... => tree rewired in spr_move_core

      // 3.1.3 count the hypothetical sources of this destination
      auto branches = Scratch_vector<Node_index>{scratch_};
      enumerate_descendant_branches_straddling(tree_, SS, old_P_t, X, std::back_inserter(branches));

      auto alpha_old_to_new = 1.0;
      auto alpha_new_to_old = 1.0 / std::ssize(branches);

      spr_move_core(X, Phylo_tree_loc{SS, new_P_t}, alpha_new_to_old / alpha_old_to_new);
    } else {
      // just change the node height
      spr_move_core(X, Phylo_tree_loc{S, new_P_t}, 1.0);
    }

  } else {
    // 4 if we are sliding the subtree down.

    // 4.0 is it a valid move?
    if (new_P_t > tree_.at(X).t) {
      return;
    }

    // 4.1 will the move change the topology
    if (new_P_t > tree_.at(S).t) {

      auto branches = Scratch_vector<Node_index>{scratch_};
      enumerate_descendant_branches_straddling(tree_, P, new_P_t, X, std::back_inserter(branches));

      // if no valid destinations then return a failure
      if (branches.empty()) {
        return;
      }

      // pick a random parent/child destination edge uniformly from options
      auto branch_index = absl::Uniform(absl::IntervalClosedOpen, bitgen_, 0, std::ssize(branches));
      auto SS = branches[branch_index];

      // 4.1.1 if iP was root... => tree rewired in spr_move_core
      auto alpha_old_to_new = 1.0 / static_cast<double>(branches.size());
      auto alpha_new_to_old = 1.0;

      spr_move_core(X, Phylo_tree_loc{SS, new_P_t}, alpha_new_to_old / alpha_old_to_new);
    } else {
      // just change the node height
      spr_move_core(X, Phylo_tree_loc{S, new_P_t}, 1.0);
    }
  }
}

auto Subrun::wilson_balding_move() -> void {
  if (std::ssize(tree_) < 2) return;

  // Comments below purposely track those in WilsonBalding#proposeTree() in BEAST

  // choose a random node avoiding root
  auto X = Node_index{};
  do {
    X = pick_random_node();
  } while (tree_.root == X);
  auto P = tree_.at(X).parent;

  // choose another random node to insert X above
  auto SS = pick_random_node();  // new sibling of X after move
  auto GG = tree_.at(SS).parent;  // new grandparent of X after move
  
  // make sure that the target branch <GG, SS> is above the subtree being moved
  while ((GG != k_no_node && tree_.at(GG).t >= tree_.at(X).t) || (X == SS)) {
    SS = pick_random_node();
    GG = tree_.at(SS).parent;
  }

  // disallow moves that change the root
  if (SS == tree_.root || P == tree_.root) {
    return;
  }

  if (GG == P || SS == P || GG == X) {
    return;
  }

  auto S = tree_.at(P).sibling_of(X);
  auto G = tree_.at(P).parent;
  auto new_max_age = std::min(tree_.at(X).t, tree_.at(SS).t);
  auto new_range = new_max_age - tree_.at(GG).t;
  auto new_t_P = absl::Uniform(bitgen_, new_max_age - new_range, new_max_age);
  auto old_max_age = std::min(tree_.at(X).t, tree_.at(S).t);
  auto old_range = old_max_age - tree_.at(G).t;

  spr_move_core(X, Phylo_tree_loc{SS, new_t_P}, new_range / old_range);
}

auto Subrun::spr1_move() -> void {
  if (std::ssize(tree_) < 2) return;

  auto chooser = absl::Uniform(bitgen_, 0.0, 1.0);
  auto limit =
      chooser < 0.01 ? std::numeric_limits<int>::max() :
      //chooser < 0.02 ? 2 :
      1;
  
  // Effective mu for all Jukes-Cantor-like calculations
  auto mu_JC = lambda_i_.at(tree_.root) / (tree_.num_sites() - num_sites_missing_at_every_node_.at(tree_.root));
  //std::cerr << absl::StreamFormat("mu_JC = %.3g x 10^3 / site / year\n",
  //                                mu_JC * 365.0 * 1e3);

  // Instead of picking it exactly in the case of a Jukes-Cantor model, we dampen the proposal
  // probability by an exponent `f` < 1.  See comments in `Spr_Study` constructor for details.
  // That way, we systematically somewhat under-propose the low-mutation regions found in the study,
  // which gives the MH criterion ample wiggle room to add corrections without accidentally
  // rejecting an otherwise great move
  auto annealing_factor = 0.8;

  // choose a random node avoiding root
  auto X = Node_index{};
  do {
    X = pick_random_node();
  } while (tree_.root == X);
  auto t_X = tree_.at(X).t;
  auto P = tree_.at(X).parent;
  auto old_t_P = tree_.at(P).t;
  auto old_S = tree_.at(P).sibling_of(X);
  auto old_G = tree_.at(P).parent;  // may be k_no_node

  auto pruning_changes_root = P == tree_.root;
  
  if (pruning_changes_root && not includes_run_root_) {
    // This move would change the root, and we can't do that here
    return;
  }

  // 0. Kick off SPR move by cataloguing and removing all of X's upstream mutations
  // FIXME: There *must* be a way to share the spr move core here
  auto spr = Spr_move{tree_, mu_JC, includes_run_root_,
    evo_, lambda_i_, ref_cum_Q_l_, num_sites_missing_at_every_node_, scratch_};

  auto old_graft = spr.analyze_graft(X);
  spr.peel_graft(old_graft);
  if (pruning_changes_root) {
    CHECK(tree_.at(old_S).mutations.empty());  // This is part of what `peel_graft` does
    CHECK(tree_.at(P).mutations.empty() || tree_.at(P).mutations.back().t == -std::numeric_limits<double>::max());
  }
  
  // 1. Perform an SPR study around P
  auto old_min_muts = spr.count_min_mutations(old_graft);
  auto old_num_muts = spr.count_closed_mutations(old_graft);
  auto old_deltas_P_to_X = spr.summarize_closed_mutations(old_graft);
  auto missing_at_X = reconstruct_missing_sites_at(tree_, X, scratch_);

  auto pre_builder = Spr_study_builder{tree_, X, t_X, missing_at_X, scratch_};
  pre_builder.max_muts_from_start = limit; // Stop when crossing over more than 1 mutation (makes move local)
  pre_builder.seed_fill_from(old_S, 0, std::move(old_deltas_P_to_X), includes_run_root_);

  auto pre_study = Spr_study{std::move(pre_builder), lambda_i_.at(X), annealing_factor, t_X};

  
  // 2. Use SPR study to pick new nexus
  auto new_nexus_region = pre_study.pick_nexus_region(bitgen_);
  auto new_S = pre_study.candidate_regions[new_nexus_region].branch;
  CHECK_NE(new_S, P);
  auto new_t_P = pre_study.pick_time_in_region(new_nexus_region, bitgen_);
  auto log_alpha_old_to_new = pre_study.log_alpha_in_region(new_nexus_region, new_t_P);

  auto new_nexus_min_muts = pre_study.candidate_regions[new_nexus_region].min_muts;
  
  auto t_new_S = tree_.at(new_S).t;
  auto new_G = tree_.at(new_S).parent;  // may be k_no_node
  if (new_G == P) { new_G = old_G; }  // Sliding along the same branch
  auto t_new_G = (new_G == k_no_node) ? -std::numeric_limits<double>::max() : tree_.at(new_G).t;
  
  // Forbid zero-length branches
  if (new_t_P == t_X || new_t_P == t_new_S || new_t_P == t_new_G) {
    spr.apply_graft(old_graft);
    return;
  }

  
  // 3. Pre-flight SPR move
  spr.move(X, new_S, new_t_P);
  auto new_graft = spr.propose_new_graft(X, bitgen_);

  CHECK_EQ(tree_.at(X).parent, P);

  auto grafting_changes_root = P == tree_.root;
  CHECK(includes_run_root_ || not grafting_changes_root);
  if (grafting_changes_root) {
    CHECK(tree_.at(new_S).mutations.empty());  // No mutations should have been added by spr.move()
    CHECK(tree_.at(P).mutations.empty() || tree_.at(P).mutations.back().t == -std::numeric_limits<double>::max());
  }

  // 4. Perform an SPR study from the new position to calculate probability of reverse proposal
  auto new_min_muts = spr.count_min_mutations(new_graft);
  auto new_num_muts = spr.count_closed_mutations(new_graft);
  auto new_deltas_P_to_X = spr.summarize_closed_mutations(new_graft);
  auto post_builder = Spr_study_builder{tree_, X, t_X, missing_at_X, scratch_};
  post_builder.max_muts_from_start = limit; // Stop when crossing over more than 1 mutation (makes move local)
  post_builder.seed_fill_from(new_S, 0, std::move(new_deltas_P_to_X), includes_run_root_);
  
  auto post_study = Spr_study{std::move(post_builder), lambda_i_.at(X), annealing_factor, t_X};

  auto old_nexus_region = post_study.find_region(old_S, old_t_P);
  CHECK_NE(old_nexus_region, -1);
  auto log_alpha_new_to_old = post_study.log_alpha_in_region(old_nexus_region, old_t_P);


  // 4bis. Sanity checking
  CHECK_EQ(new_min_muts, pre_study.candidate_regions[new_nexus_region].min_muts);
  CHECK_EQ(old_min_muts, post_study.candidate_regions[old_nexus_region].min_muts);
  
  
  // 5. Accept/reject?
  auto delta_log_augmented_coalescent_prior =
      coalescent_prior_part_->calc_delta_partial_log_prior_after_displace_coalescence(
          old_t_P, new_t_P);
  auto log_mh = (new_graft.delta_log_G - new_graft.log_alpha_mut)
      - (old_graft.delta_log_G - old_graft.log_alpha_mut)
      + log_alpha_new_to_old - log_alpha_old_to_new
      + delta_log_augmented_coalescent_prior;

  auto accept = log_mh >= 0.0 || absl::Uniform(bitgen_, 0.0, 1.0) < std::exp(log_mh);

  if (false && (pruning_changes_root || grafting_changes_root)) {
    std::cerr << absl::StreamFormat("old_t_P = %g, new_t_P = %g, accept = %s\n",
                                    old_t_P, new_t_P, accept ? "true" : "false");
  }
  
  if (false && (pruning_changes_root || grafting_changes_root) && not accept) { // && (new_num_muts < old_num_muts-5)) {
    auto mu = evo_.mu_l(0);
    auto L = tree_.num_sites();

    if (pruning_changes_root || grafting_changes_root) {
      std::cerr << "****** ";
    } else if (old_graft.branch_infos.size() > 1 || new_graft.branch_infos.size() > 1) {
      std::cerr << "++++++ ";
    }
    std::cerr << absl::StreamFormat("Odd rejection: old %d, new_min %d, new %d, t_n-t_o=%g, log_mh=%g\n",
                                    old_num_muts, new_nexus_min_muts, new_num_muts, new_t_P - old_t_P, log_mh);
    std::cerr << absl::StreamFormat(" - delta_log_augmented_coalescent_prior: %g\n",
                                    delta_log_augmented_coalescent_prior);
    std::cerr << absl::StreamFormat(" - new_graft.delta_log_G = %g, old_graft.delta_log_G = %g => \n",
                                    new_graft.delta_log_G, old_graft.delta_log_G);
    std::cerr << absl::StreamFormat("   log(G_n/G_0) = %g (cf JC expectation = %g)\n",
                                    new_graft.delta_log_G - old_graft.delta_log_G,
                                    -mu*L * (old_t_P - new_t_P)
                                    + (new_num_muts - old_num_muts)*std::log(mu / 3.0));
    std::cerr << absl::StreamFormat(" - new_graft.log_alpha_mut = %g, old_graft.log_alpha_mut = %g =>\n",
                                    new_graft.log_alpha_mut, old_graft.log_alpha_mut);
    std::cerr << absl::StreamFormat("   log(a_m(n->o)/a_m(o->n)) = %g (cf JC expectation = %g)\n",
                                    old_graft.log_alpha_mut - new_graft.log_alpha_mut,
                                    -old_num_muts*std::log(t_X - old_t_P) + new_num_muts*std::log(t_X - new_t_P));
    std::cerr << absl::StreamFormat(" - log_alpha_new_to_old = %g, log_alpha_old_to_new = %g =>\n",
                                    log_alpha_new_to_old, log_alpha_old_to_new);
    std::cerr << absl::StreamFormat("   log(a_SPR(n->o)/a_SPR(o->n)) = %g (cf JC expectation = %g)\n",
                                    log_alpha_new_to_old - log_alpha_old_to_new,
                                    -mu*L * (new_t_P - old_t_P)
                                    + old_num_muts*std::log(mu*(t_X - old_t_P)/3.0)
                                    - new_num_muts*std::log(mu*(t_X - new_t_P)/3.0));
    std::cerr << absl::StreamFormat("                                     (cf JC expectation with f = %g)\n",
                                    annealing_factor * (
                                        -mu*L * (new_t_P - old_t_P)
                                        + old_num_muts*std::log(mu*(t_X - old_t_P)/3.0)
                                        - new_num_muts*std::log(mu*(t_X - new_t_P)/3.0)));
  }
  
  if (accept) {
    spr.apply_graft(new_graft);
    log_G_ -= old_graft.delta_log_G;
    log_G_ += new_graft.delta_log_G;
    log_augmented_coalescent_prior_ += delta_log_augmented_coalescent_prior;
    coalescent_prior_part_->coalescence_displaced(old_t_P, new_t_P);
  } else {
    spr.move(X, old_S, old_t_P);
    spr.apply_graft(old_graft);
  }
}

// SPR = Subtree Pruning and Regrafting
// Prune the subtree rooted at node X and regraft it to just above SS at time new_P_t
// (which must be between tree.at_parent_of(SS).t and tree.at(SS).t, calculated *after* X is pruned).
// If the SPR move is accepted, then the nodes SS and X become siblings.
// The choice of this configuration has accumulated a Hastings ratio
// of alpha_new_to_old_over_alpha_old_to_new
auto Subrun::spr_move_core(
    Node_index X,
    Phylo_tree_loc new_nexus,
    double alpha_new_to_old_over_alpha_old_to_new)
    -> void {

  if (X == tree_.root) { return; }  // Really can't do this!
  if (not includes_run_root_) {
    if (tree_.at(X).parent == tree_.root || new_nexus.branch == tree_.root) {
      // This move would change the root, and we can't do that here
      return;
    }
  }

  auto t_X = tree_.at(X).t;
  auto P = tree_.at(X).parent;
  auto old_t_P = tree_.at(P).t;
  auto old_S = tree_.at(P).sibling_of(X);
  auto new_t_P = new_nexus.t;

  // Forbid zero-length branches
  if (new_t_P == t_X ||
      (new_t_P == tree_.at(new_nexus.branch).t) ||
      (P != tree_.root && new_t_P == tree_.at_parent_of(P).t)) {  // FIXME: This should look at G', not G!
    return;
  }
  
  auto mu_JC = lambda_i_.at(tree_.root) / (tree_.num_sites() - num_sites_missing_at_every_node_.at(tree_.root));
  //std::cerr << absl::StreamFormat("mu_JC = %.3g x 10^3 / site / year\n",
  //                                mu_JC * 365.0 * 1e3);
  auto spr = Spr_move{tree_, mu_JC, includes_run_root_,
    evo_, lambda_i_, ref_cum_Q_l_, num_sites_missing_at_every_node_, scratch_};

  auto old_graft = spr.analyze_graft(X);
  spr.peel_graft(old_graft);
  spr.move(X, new_nexus.branch, new_nexus.t);
  auto new_graft = spr.propose_new_graft(X, bitgen_);
  
  // Accept/reject?
  auto delta_log_augmented_coalescent_prior =
      coalescent_prior_part_->calc_delta_partial_log_prior_after_displace_coalescence(
          old_t_P, new_nexus.t);
  auto log_mh = (new_graft.delta_log_G - new_graft.log_alpha_mut)
      - (old_graft.delta_log_G - old_graft.log_alpha_mut)
      + std::log(alpha_new_to_old_over_alpha_old_to_new)
      + delta_log_augmented_coalescent_prior;

  auto accept = log_mh >= 0.0 || absl::Uniform(bitgen_, 0.0, 1.0) < std::exp(log_mh);
  
  if (accept) {
    spr.apply_graft(new_graft);
    log_G_ -= old_graft.delta_log_G;
    log_G_ += new_graft.delta_log_G;
    log_augmented_coalescent_prior_ += delta_log_augmented_coalescent_prior;
    coalescent_prior_part_->coalescence_displaced(old_t_P, new_nexus.t);
  } else {
    spr.move(X, old_S, old_t_P);
    spr.apply_graft(old_graft);
  }
}

}  // namespace delphy
