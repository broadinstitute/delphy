#include "run.h"

#include "distributions.h"
#include "tree_partitioning.h"
#include "phylo_tree_calc.h"

#include <boost/math/special_functions/gamma.hpp>

namespace delphy {

Run::Run(ctpl::thread_pool& thread_pool, std::mt19937 bitgen, Phylo_tree tree)
    : thread_pool_{&thread_pool},
      bitgen_{bitgen},
      tree_{std::move(tree)},
      num_parts_{1},
      target_coal_prior_cells_{400},
      pop_model_{std::make_shared<Exp_pop_model>(calc_max_tip_time(tree_), 1000.0, 0.0)},
      skygrid_tau_{1.0},  // Prior mean, Gill et al 2012 Eq. 15
      alpha_{10.0},
      nu_(tree_.num_sites(), 1.0),
      evo_{make_single_partition_global_evo_model(tree_.num_sites())},
      coalescent_prior_{pop_model_, tree_.size(), calc_max_tip_time(tree_), 1.0} {

  // In Debug builds, record input tip sequences so that we can check that MCMC moves don't result
  // in incompatible real tip sequences
  save_original_sequences();
  
  // MCMC control
  step_ = 0;
  set_local_moves_per_global_move(-1);  // use default
  next_global_move_step_ = -1;
  next_partition_stencil_refresh_step_ = -1;
  
  // Initialize coalescent prior
  for (const auto& node : index_order_traversal(tree_)) {
    if (tree_.at(node).is_tip()) {
      coalescent_prior_.mark_as_tip(node);
    } else {
      coalescent_prior_.mark_as_coalescence(node);
    }
  }

  // Pick a reasonable starting point for stationary frequencies of HKY model
  auto state_frequencies_of_ref_sequence = calc_cur_state_frequencies_of_ref_sequence();
  auto total_present_sites = static_cast<double>(estd::ranges::sum(state_frequencies_of_ref_sequence));
  auto est_pi = Seq_vector<double>{};
  for (auto a : k_all_real_seq_letters) {
    est_pi[a] = state_frequencies_of_ref_sequence[a] / total_present_sites;
  }
  auto crazy = false;
  for (auto a : k_all_real_seq_letters) {
    if (est_pi[a] < 0.01 || est_pi[a] > 0.99) {
      crazy = true;
    }
  }
  if (crazy) {
    // Just start with something uniform
    std::ranges::fill(est_pi, 0.25);
  }
  hky_model_.mu = 1e-3 / 365.0;
  hky_model_.kappa = 1.0;
  hky_model_.pi_a = est_pi;

  repartition();

  invalidate_derived_quantities();
}

auto Run::refresh_partition_stencils() -> void {
  // Make 10 partition stencils very occasionally; during the run, we switch between these stencils randomly.
  // It's important for correctness that all partitions are equally likely to be chosen during a run.
  // If we change the set of partitions over which this choice takes place *slowly* (~ less than once
  // per effective sample), then it doesn't matter *how* those partitions are chosen.  But if we do it
  // more quickly, the sampling becomes biased.

  if (partition_stencils_valid_ && step_ < next_partition_stencil_refresh_step_) {
    return;
  }

  //std::cerr << "\n\n*****\nREFRESHING PARTITION STENCILS\n*****\n\n";

  partition_stencils_.clear();
  for (auto i = 0; i != 10; ++i) {
    partition_stencils_.push_back(generate_random_partition_stencil(tree_, num_parts_, bitgen_));
  }

  partition_stencils_valid_ = true;
  next_partition_stencil_refresh_step_ = step_ + 200 * local_moves_per_global_move_;
  force_repartition_ = true;
}

auto Run::repartition() -> void {
  // Ensure we have enough bitgens around
  while (std::ssize(subbitgens_) < num_parts_) {
    subbitgens_.emplace_back(bitgen_());
  }

  // Ensure partition stencils are in up to date
  refresh_partition_stencils();

  //std::cerr << "*** REPARTITION ***" << std::endl;
  auto partition_stencils_index = std::uniform_int_distribution<int>{0, int(std::ssize(partition_stencils_)-1)}(bitgen_);
  auto& partition_stencil = partition_stencils_.at(partition_stencils_index);
  tree_partition_ = partition_tree(tree_, partition_stencil);

  // Before we propagate the reference sequence to subruns, normalize
  CHECK(tree_.at_root().missations.from_states.empty());
  normalize_root();

  // A good time to check for consistency, even in release builds
  assert_phylo_tree_integrity(tree_, true);
  
  // Make subruns
  subruns_.clear();
  for (auto& partition_part : tree_partition_.parts()) {
    auto subroot = partition_part.info().cut_point;
    auto subroot_seq_view = view_of_sequence_at(tree_, subroot);
    const auto& ref_seq = tree_.ref_sequence;

    auto subtree = Phylo_tree{partition_part.size()};
    subtree.ref_sequence = ref_seq;
    copy_topology(partition_part, subtree);

    // Missations in above subroot logically come before mutations above subroot,
    // so all their from_states match the reference sequence
    auto root_missation_map = Missation_map<>{};
    root_missation_map.intervals = reconstruct_missing_sites_at(tree_, subroot);
    root_missation_map.from_states.clear();
    
    auto root_mutations = Mutation_list<>{};
    for (const auto& [l, b] : subroot_seq_view.deltas()) {
      if (not root_missation_map.contains(l)) {
        root_mutations.push_back(Mutation{ref_seq[l], l, b, -std::numeric_limits<double>::max()});
      }
    }
    sort_mutations(root_mutations);

    const auto empty_mutation_list = Mutation_list<>{};
    auto original_sequences = Node_vector<Sequence_overlay>{};
    auto original_missing_sites_all = Node_vector<Interval_set<>>{};
    for (const auto& partition_node : index_order_traversal(partition_part)) {
      auto node = partition_part.at(partition_node).orig_tree_index();
      auto subtree_node = partition_node;
      const auto& mutation_list = (node == subroot) ? root_mutations : tree_.at(node).mutations;
      const auto& missation_map = (node == subroot) ? root_missation_map : tree_.at(node).missations;
      subtree.at(subtree_node).name = tree_.at(node).name;
      subtree.at(subtree_node).t = tree_.at(node).t;
      if (subtree.at(subtree_node).is_tip() && not tree_.at(node).is_tip()) {
        // "tip" in subtree is really just a frozen inner node, which we shouldn't move
        subtree.at(subtree_node).t_min = tree_.at(node).t;
        subtree.at(subtree_node).t_max = tree_.at(node).t;
      } else {
        subtree.at(subtree_node).t_min = tree_.at(node).t_min;
        subtree.at(subtree_node).t_max = tree_.at(node).t_max;
      }
      subtree.at(subtree_node).mutations = mutation_list;
      subtree.at(subtree_node).missations = missation_map;
      if (estd::is_debug_enabled) {
        original_sequences.push_back(original_sequences_[node]);
        original_missing_sites_all.push_back(original_missing_sites_all_[node]);
      }
    }

    auto& subrun = subruns_.emplace_back(
        absl::BitGenRef{subbitgens_[subruns_.size()]}, std::move(subtree), subroot == tree_.root, evo_);
    if (estd::is_debug_enabled) {
      subrun.set_original_sequences(std::move(original_sequences), std::move(original_missing_sites_all));
    }
  }

  // Start very_scalable_coalescent_prior from scratch
  reset_very_scalable_coalescent_parts();

  force_repartition_ = false;
}

auto Run::reassemble() -> void {
  CHECK_EQ(subruns_.size(), tree_partition_.parts().size());
  auto num_partitions = std::ssize(tree_partition_.parts());

  for (auto i = 0; i != num_partitions; ++i) {
    const auto& partition_part = tree_partition_.parts()[i];
    const auto& subrun = subruns_[i];

    assert_phylo_tree_integrity(subrun.tree());
    CHECK(subrun.tree().at_root().missations.from_states.empty());

    // Transpose all the data in the subrun onto the main tree
    const auto& subtree = subrun.tree();
    for (const auto& subnode : post_order_traversal(subtree)) {
      auto node = partition_part.at(subnode).orig_tree_index();
      tree_.at(node).t = subtree.at(subnode).t;
      if (tree_.at(node).is_inner_node()) {
        coalescent_prior_.displace_coalescence(node, tree_.at(node).t);
      } else {
        coalescent_prior_.displace_tip(node, tree_.at(node).t);
      }
      if (subnode != subtree.root) {
        tree_.at(node).mutations = subtree.at(subnode).mutations;
        tree_.at(node).missations = subtree.at(subnode).missations;
      }

      // Transpose topology
      // Note: the test is not `tree_.at(node).is_inner_node()` to only handle cut points that are subroots, not subtips
      if (subtree.at(subnode).is_inner_node()) {
        auto subleft = subtree.at(subnode).left_child();
        auto subright = subtree.at(subnode).right_child();
        auto left = partition_part.at(subleft).orig_tree_index();
        auto right = partition_part.at(subright).orig_tree_index();
        tree_.at(node).children = {left, right};
        tree_.at(left).parent = node;
        tree_.at(right).parent = node;

        if (estd::is_debug_enabled) {
          CHECK_LE(tree_.at(node).t, tree_.at(left).t);
          CHECK_LE(tree_.at(node).t, tree_.at(right).t);
        }
      }
    }

    // If the subrun includes the root node, transpose the root index & deltas from ref sequence back
    if (subrun.includes_run_root()) {
      const auto& subtree = subrun.tree();
      auto subroot = subtree.root;
      auto subroot_in_main = partition_part.at(subroot).orig_tree_index();
      tree_.root = subroot_in_main;
      tree_.at_root().parent = k_no_node;
      tree_.at_root().mutations = subtree.at(subroot).mutations;
      tree_.at_root().missations = subtree.at(subroot).missations;
      DCHECK(tree_.ref_sequence == subtree.ref_sequence);  // subrun ref sequence unchanged
    }
  }

  assert_phylo_tree_integrity(tree_);

  invalidate_derived_quantities();
  check_global_and_local_totals_match();
}

auto Run::normalize_root() -> void {
  // Merge any mutations above the root into the reference sequence
  for (const auto& m : tree_.at_root().mutations) {
    --state_frequencies_of_ref_sequence_[m.from];
    ++state_frequencies_of_ref_sequence_[m.to];
  }
  rereference_to_root_sequence(tree_);
}

auto Run::push_global_params_to_subruns() -> void {
  for (auto& subrun : subruns_) {
    subrun.set_evo(evo_);
    subrun.set_only_displacing_inner_nodes(only_displacing_inner_nodes_);
    subrun.set_topology_moves_enabled(topology_moves_enabled_);
  }
  reset_very_scalable_coalescent_parts();  // TODO: Should only push kbars and k(t)?
}

auto Run::reset_very_scalable_coalescent_parts() -> void {
  // Make very scalable coalescent
  auto phylo_subtrees = std::vector<const Phylo_tree*>{};
  for (const auto& subrun : subruns_) {
    phylo_subtrees.push_back(&subrun.tree());
  }
  coalescent_prior_parts_ = make_very_scalable_coalescent_prior_parts(
      phylo_subtrees,
      tree_partition_.root_part_index(),
      pop_model_,
      subbitgens_,
      coalescent_prior_.t_step()  // NOTE: This picks up any changes to t_step midway
  );
  for (auto i = 0; i != std::ssize(phylo_subtrees); ++i) {
    subruns_[i].set_coalescent_prior_part(&coalescent_prior_parts_[i]);
  }
}

auto Run::set_coalescent_prior_t_step(double t_step) -> void {
  //std::cerr << absl::StreamFormat("\n\n*****\n\nResetting coalescent prior resolution to %g days\n\n*******\n\n", t_step);
  coalescent_prior_.reset(t_step);
  log_coalescent_prior_ = calc_cur_log_coalescent_prior();
  // Very scalable coalescent should pick up new t_step in reset_very_scalable_coalescent_parts()
}

auto Run::recalc_derived_quantities() const -> void {
  // Derived quantities
  derive_evo();
  log_G_ = calc_cur_log_G();
  Ttwiddle_beta_a_ = calc_cur_Ttwiddle_beta_a();
  num_muts_ = calc_cur_num_muts();
  num_muts_ab_ = calc_cur_num_muts_ab();
  log_coalescent_prior_ = calc_cur_log_coalescent_prior();
  log_other_priors_ = calc_cur_log_other_priors();
  state_frequencies_of_ref_sequence_ = calc_cur_state_frequencies_of_ref_sequence();
  last_revalidation_step_ = step_;
}

auto Run::check_derived_quantities() const -> void {
  if (estd::is_debug_enabled) {
    // We maintain these incrementally except when we force a recalculation
    { auto expected = calc_cur_log_G();
      CHECK(std::abs(log_G_ - expected) < 1e-6)
          << log_G_ << " != " << expected; }
    { auto expected = calc_cur_Ttwiddle_beta_a();
      CHECK((taxicab_distance(Ttwiddle_beta_a_, expected)
             / static_cast<double>(tree_.num_sites())) < 1e-8)
          << absl::StrJoin(Ttwiddle_beta_a_, ",", absl::StreamFormatter())
          << " != "
          << absl::StrJoin(expected, ",", absl::StreamFormatter()); }
    CHECK_EQ(num_muts_, calc_cur_num_muts());
    CHECK_EQ(num_muts_ab_, calc_cur_num_muts_ab());
    { auto expected = calc_cur_log_coalescent_prior();
      CHECK(std::abs(log_coalescent_prior_ - expected) < 1e-8)
          << log_coalescent_prior_ << " != " << expected; }
    { auto expected = calc_cur_log_other_priors();
      CHECK(std::abs(log_other_priors_ - expected) < 1e-8)
          << log_other_priors_ << " != " << expected; }
    CHECK_EQ(state_frequencies_of_ref_sequence_, calc_cur_state_frequencies_of_ref_sequence());
  }
}

auto Run::check_global_and_local_totals_match() const -> void {
  if (estd::is_debug_enabled) {
    const auto& subrun_with_root = subruns_[tree_partition_.root_part_index()];
    auto sub_log_G = 0.0;
    auto sub_log_augmented_coalescent_prior = 0.0;
    for (const auto& subrun : subruns_) {
      sub_log_G += subrun.log_G();
      sub_log_augmented_coalescent_prior += subrun.log_augmented_coalescent_prior();
    }
    CHECK(std::abs(log_G() - sub_log_G) < 1e-6)
        << log_G() << " != " << sub_log_G;
    auto log_augmented_coalescent_prior = calc_cur_log_augmented_coalescent_prior();
    CHECK(std::abs(log_augmented_coalescent_prior - sub_log_augmented_coalescent_prior) < 1e-6)
        << log_augmented_coalescent_prior << " != " << sub_log_augmented_coalescent_prior;
    CHECK_EQ(state_frequencies_of_ref_sequence(),
             subrun_with_root.state_frequencies_of_ref_sequence_per_partition()[0]);
  }
}

auto Run::set_mpox_hack_enabled(bool mpox_hack_enabled) -> void {
  if (mpox_hack_enabled == mpox_hack_enabled_) { return; }

  mpox_hack_enabled_ = mpox_hack_enabled;
  if (mpox_hack_enabled) {
    // Set up two partitions as described in run.h
    auto first_tip = 0;
    CHECK(tree_.at(first_tip).is_tip());
    auto seq = view_of_sequence_at(tree_, first_tip);

    auto count_apobec_context = 0;
    for (auto l = Site_index{0}; l != tree_.num_sites(); ++l) {
      auto apobec_context =
          ((l != 0) &&
           (seq[l-1] == Real_seq_letter::T) &&
           ((seq[l] == Real_seq_letter::C) ||      // not yet mutated
            (seq[l] == Real_seq_letter::T))) ||    // already mutated
          ((l+1 != tree_.num_sites()) &&
           (seq[l+1] == Real_seq_letter::A) &&
           ((seq[l] == Real_seq_letter::G) ||      // not yet mutated
            (seq[l] == Real_seq_letter::A)));      // already mutated
      if (apobec_context) { ++count_apobec_context; }
      evo_.partition_for_site[l] = apobec_context ? 1 : 0;
    }
    std::cerr << count_apobec_context << " of " << tree_.num_sites() << " sites have APOBEC context\n";
    
    evo_.partition_evo_model = Partition_vector<Site_evo_model>(2, Site_evo_model{});

    mpox_mu_ = hky_model_.mu;
    mpox_mu_star_ = 0.0;
  } else {
    // Set up single partition
    for (auto l = Site_index{0}; l != tree_.num_sites(); ++l) {
      evo_.partition_for_site[l] = 0;
    }
    evo_.partition_evo_model = Partition_vector<Site_evo_model>(1, Site_evo_model{});
    hky_model_.mu = mpox_mu_;
  }
  invalidate_derived_quantities();
}

auto Run::derive_evo() const -> void {
  if (not mpox_hack_enabled_) {
    CHECK_EQ(evo_.num_partitions(), 1);
    evo_.partition_evo_model[0] = hky_model_.derive_site_evo_model();
  } else {
    CHECK_EQ(evo_.num_partitions(), 2);

    // Partition 0 = Sites without APOBEC context
    auto& site_evo_0 = evo_.partition_evo_model[0];
    site_evo_0.mu = mpox_mu_;
    site_evo_0.pi_a = {0.25, 0.25, 0.25, 0.25};
    for (auto a : k_all_real_seq_letters) {
      for (auto b : k_all_real_seq_letters) {
        // Only polymerase-driven errors, following a simple Jukes-Cantor model
        site_evo_0.q_ab[a][b] = ((a == b) ? -1.0 : 1./3.);
      }
    }

    // Partition 1 = Sites with APOBEC context
    auto& site_evo_1 = evo_.partition_evo_model[1];

    // Polymerase-driven mutations are still there
    site_evo_1 = site_evo_0;
    site_evo_1.mu = mpox_mu_;
    auto rho = mpox_mu_star_ / mpox_mu_;
    
    // APOBEC-driven TC -> TT mutation
    site_evo_1.q_ab[Real_seq_letter::C][Real_seq_letter::T] += 2*rho;
    site_evo_1.q_ab[Real_seq_letter::C][Real_seq_letter::C] -= 2*rho;
    
    // APOBEC-driven GA -> AA mutation (== TC -> TT on the complementary strand)
    site_evo_1.q_ab[Real_seq_letter::G][Real_seq_letter::A] += 2*rho;
    site_evo_1.q_ab[Real_seq_letter::G][Real_seq_letter::G] -= 2*rho;
  }
  evo_.nu_l = nu_;
}

auto Run::calc_cur_log_G() const -> double {
  auto result = 0.0;
  result += calc_log_root_prior(tree_, evo_);
  result += calc_log_G_below_root(tree_, evo_);
  return result;
}

auto Run::calc_cur_Ttwiddle_beta_a() const -> Partition_vector<Seq_vector<double>> {
  return calc_Ttwiddle_beta_a(tree_, evo_);
}

auto Run::calc_cur_num_muts() const -> int {
  return calc_num_muts(tree_);
}

auto Run::calc_cur_num_muts_ab() const -> Seq_matrix<int> {
  return calc_num_muts_ab(tree_);
}

auto Run::calc_cur_log_coalescent_prior() const -> double {
  for (const auto& node : index_order_traversal(tree_)) {
    if (tree_.at(node).is_tip()) {
      coalescent_prior_.displace_tip(node, tree_.at(node).t);
    } else {
      coalescent_prior_.displace_coalescence(node, tree_.at(node).t);
    }
  }
  return coalescent_prior_.calc_log_prior();
}

auto Run::calc_cur_log_augmented_coalescent_prior() const -> double {
  auto result = 0.0;
  CHECK(not coalescent_prior_parts_.empty());
  for (const auto& part : coalescent_prior_parts_) {
    result += part.calc_partial_log_prior();
  }
  return result;
}

auto Run::calc_cur_state_frequencies_of_ref_sequence() const -> Seq_vector<int> {
  return calc_state_frequencies_per_partition_of(tree_.ref_sequence, evo_)[0];
}

auto Run::calc_cur_log_other_priors() const -> double {
  auto log_prior = 0.0;

  if (not mpox_hack_enabled_) {
    // Mu - Uniform prior
  } else {
    // Mu & Mu* - Uniform priors
  }

  // Alpha - Exponential prior with mean 1.0
  const auto mean_alpha = 1.0;
  log_prior += -alpha_ / mean_alpha - std::log(mean_alpha);
  
  // nu_l ~ Gamma(alpha, alpha) = (alpha^alpha / Gamma(alpha)) (\nu_l)^{\alpha - 1} exp(-\alpha nu_l)
  auto L = tree_.num_sites();
  log_prior += L * (alpha_ * std::log(alpha_) - std::lgamma(alpha_));
  auto sum_nu = 0.0;
  auto sum_log_nu = 0.0;
  for (auto nu_l : nu_) {
    DCHECK_NE(nu_l, 0.0);
    sum_nu += nu_l;
    sum_log_nu += std::log(nu_l);
  }
  log_prior += (alpha_ - 1) * sum_log_nu - alpha_ * sum_nu;

  if (not mpox_hack_enabled_) {
    // HKY Freqs - Uniform prior
      
    // HKY Kappa - log-normal prior, so log(kappa) has mean 1 and stddev 1.25
    //  pi(kappa) = exp(-(log(kappa) - mean_log_kappa)^2 / (2 sigma_log_kappa^2)) / (kappa * sqrt(2 pi sigma_log_kappa^2))
    const auto mean_log_kappa = 1.0;
    const auto sigma_log_kappa = 1.25;
    log_prior += -std::pow(std::log(hky_kappa()) - mean_log_kappa, 2) / (2 * sigma_log_kappa * sigma_log_kappa)
        - 0.5 * std::log(2 * M_PI * sigma_log_kappa * sigma_log_kappa) - std::log(hky_kappa());
  }

  if (typeid(*pop_model_) == typeid(Exp_pop_model)) {
    const auto& exp_pop_model = static_cast<const Exp_pop_model&>(*pop_model_);
    
    // pop_n0 - 1/x prior
    log_prior -= std::log(exp_pop_model.pop_at_t0());
    
    // pop_g - Laplace prior on the growth rate, with mu 0.001/365 and scale 30.701135/365
    //      pdf(g; mu, scale) = 1/(2*scale) exp(-|g - mu| / scale)
    const auto mu_g = 0.001 / 365.0;
    const auto scale_g = 30.701135 / 365.0;
    log_prior += -std::abs(exp_pop_model.growth_rate() - mu_g) / scale_g - std::log(2 * scale_g);
    
  } else if (typeid(*pop_model_) == typeid(Skygrid_pop_model)) {
    const auto& skygrid_pop_model = static_cast<const Skygrid_pop_model&>(*pop_model_);

    // tau - gamma prior (Gill et al 2012, Eq. 15)
    const auto prior_tau_alpha = 0.001;
    const auto prior_tau_beta = 0.001;

    const auto log_tau = std::log(skygrid_tau_);
    
    log_prior +=
        (prior_tau_alpha - 1) * log_tau
        - prior_tau_beta * skygrid_tau_;

    // gamma - GMRF prior (Gill et al 2012, Eq. 13)
    //
    // Since time intervals are in principle not all the same, we really should
    // be using something like a diffusion law to condition how much consecutive
    // gamma's should vary, e.g., gamma_{i+1} - gamma_i ~ N(0, 2 D (x_{i+1} - x_i)).
    // We don't do that to not deviate from the standard Skygrid parametrization;
    // assuming equal time intervals of width dt = x_1 - x_0, we'd have
    //
    //    2 D dt = 1/tau => D = 1/(2 tau dt).

    for (auto k = 0; k != skygrid_pop_model.M(); ++k) {
      auto delta_gamma = skygrid_pop_model.gamma(k+1) - skygrid_pop_model.gamma(k);
      log_prior +=
          0.5 * log_tau
          - 0.5 * delta_gamma*delta_gamma * skygrid_tau_;
    }
  }

  return log_prior;
}

static auto subrun_exec(
    int /*thread_id*/,
    Subrun& subrun,
    int64_t subrun_substeps)
    -> void {
  
  for (auto i = 0L; i != subrun_substeps; ++i) {
    subrun.mcmc_sub_iteration();
  }
}

// Do `substeps` local moves spread over the subruns, then one global move
auto Run::do_mcmc_steps(int substeps) -> void {
  auto target_step = step_ + substeps;

  // Ensure that changes to tree, such as from normalize_root() at the end of the last mcmc_iterations call
  // or any unexpected outside manipulation, are propagated to subruns.
  repartition();
  
  while (step_ < target_step) {
    
    // Correct for accumulating round-off errors
    if (step_ >= (last_revalidation_step_ + 1'000'000)) { invalidate_derived_quantities(); }
    validate_derived_quantities();
    
    if (step_ >= next_global_move_step_) {
      auto should_repartition = force_repartition_ || repartitioning_enabled_;
      if (should_repartition) { repartition(); }
      run_global_moves();
      next_global_move_step_ = step_ + local_moves_per_global_move_;
    }

    CHECK_LT(step_, next_global_move_step_);
    auto steps_before_next_global_move = std::min(next_global_move_step_, target_step) - step_;
    if (steps_before_next_global_move > 0) {
      push_global_params_to_subruns();
      check_global_and_local_totals_match();
      run_local_moves(steps_before_next_global_move);
      reassemble();
      step_ += steps_before_next_global_move;
    }
  }
  
  normalize_root();
  
  assert_cur_tip_sequences_compatible_with_original_ones();
  assert_phylo_tree_integrity(tree_, true);
}

auto Run::set_step(int64_t step) -> void {
  if (step < 0) {
    throw std::invalid_argument(absl::StrFormat(
        "Step should not be negative (got %d)", step));
  }
  step_ = step;
  next_global_move_step_ = -1;
  next_partition_stencil_refresh_step_ = -1;
}

auto Run::set_local_moves_per_global_move(int local_moves_per_global_move) -> void {
  if (local_moves_per_global_move == -1) {
    local_moves_per_global_move = 50 * std::ssize(tree_);  // A good default
  }
  if (local_moves_per_global_move < 1) {
    throw std::invalid_argument(absl::StrFormat(
        "There must be at least 1 local move per global move (got %d)", local_moves_per_global_move));
  }
  local_moves_per_global_move_ = local_moves_per_global_move;
  next_global_move_step_ = -1;
  next_partition_stencil_refresh_step_ = -1;
}

auto Run::run_local_moves(int count) -> void {
  auto subcount = count / subruns_.size();
  auto subrun_count = std::ssize(subruns_);
  auto subrun_results = std::vector<std::future<void>>{};
  for (auto i = 1; i < subrun_count; ++i) {  // Note: start at 1
    subrun_results.push_back(thread_pool_->push(subrun_exec, std::ref(subruns_[i]), subcount));
  }
  subrun_exec(0, subruns_[0], count - (subrun_count-1) * subcount);  // Subrun 0 runs on this thread
  for (auto& subrun_result : subrun_results) {
    subrun_result.get();
  }
}

auto Run::run_global_moves() -> void {
  // Global moves
  invalidate_derived_quantities();
  validate_derived_quantities();

  // Each group of global moves below is fully independent of the others and is practically a Gibbs sampling
  // of one of the global parameters, so there's no point in doing lots of these moves

  if (not mpox_hack_enabled_) {
    // 1. Gibbs sampling of mu
    // Depends on Ttwiddle_a, {q_a} and num_muts
    if (mu_move_enabled_) {
      mu_move();
      check_derived_quantities();
    }
    
    // 2. Pseudo-Gibbs sampling of evo model parameters
    // Depends only on num_muts_ab and counts of state in root sequence; after enough MCMC moves, the
    // resulting kappa & pi_a's are practically Gibbs sampled
    for (auto i = 0; i != 10; ++i) {
      hky_frequencies_move();
      check_derived_quantities();
      hky_kappa_move();
      check_derived_quantities();
    }
  } else {
    // 1 & 2. Gibbs sampling of mu & mu_star
    mpox_hack_moves();
    check_derived_quantities();
  }

  // 3. Pseudo-Gibbs sampling of alpha
  // On integrating out the nu_l's, p(alpha) depends only on {M_l} and {TTwiddle_l}.  After lots of little alpha moves,
  // alpha is practically Gibbs sampled.  Then the nu_l's really are Gibbs sampled
  if (alpha_move_enabled_) {
    alpha_moves();  // Gibbs samples all nu's as part of the alpha move
    check_derived_quantities();
  }

  // 4-pre. Adjust resolution of staircases used in coalescent prior if they're really out of whack
  auto min_t = tree_.at_root().t;
  auto max_t = calc_max_tip_time(tree_);
  auto cur_t_step = coalescent_prior_.t_step();
  auto target_t_step = (max_t - min_t) / target_coal_prior_cells_;
  auto min_t_step = (1.0) / target_coal_prior_cells_;  // Avoid degenerate cases when tree collapses
  auto ratio_cur_to_target = cur_t_step / target_t_step;
  auto retarget =  // Don't do this unless really needed
      (cur_t_step > min_t_step) &&
      (ratio_cur_to_target < 2./3. || ratio_cur_to_target > 4./3.);
  if (retarget) {
    auto new_t_step = std::max(min_t_step, 0.5 * (cur_t_step + target_t_step)); // Avoid jumping to an extreme
    set_coalescent_prior_t_step(new_t_step);
  }

  // 4. Pseudo-Gibbs sampling of population parameters
  // Depends only on tree topology; after enough MCMC moves, the
  // resulting population model parameters are practically Gibbs sampled
  for (auto i = 0; i != 50; ++i) {
    if (typeid(*pop_model_) == typeid(Exp_pop_model)) {
      if (final_pop_size_move_enabled_) {
        pop_size_move();
        check_derived_quantities();
      }
      if (pop_growth_rate_move_enabled_) {
        pop_growth_rate_move();
        check_derived_quantities();
      }
      
    } else if (typeid(*pop_model_) == typeid(Skygrid_pop_model)) {
      skygrid_tau_move();
      check_derived_quantities();
      
      skygrid_gammas_move();
      check_derived_quantities();
    }
  }
}

auto Run::mu_move() -> void {
  // We follow LeMieux et al (2021) and use a uniform prior for mu.
  // For concreteness, we assume the mutation rates across all partitions are linked,
  // though the infrastructure is already here to allow for unlinked rates.
  //
  // We can Gibbs sample a new value of mu using a Gamma with
  // shape alpha = M+1 and rate beta = Ttwiddle, where M = num_muts() and
  // Ttwiddle = sum_{beta,a} q^[beta]_a * Ttwiddle^beta_a

  CHECK_GE(evo_.num_partitions(), 1);
  for (auto beta = 0; beta != evo_.num_partitions(); ++beta) {
    CHECK_EQ(evo_.partition_evo_model[beta].mu, hky_model_.mu)
        << "Unlinked mutation rates not yet supported";
  }

  auto Ttwiddle = 0.0;
  for (auto beta = 0; beta != evo_.num_partitions(); ++beta) {
    const auto& evo_beta = evo_.partition_evo_model[beta];
    auto Ttwiddle_beta = 0.0;
    for (auto a : k_all_real_seq_letters) {
      Ttwiddle_beta += evo_beta.q_a(a) * Ttwiddle_beta_a_[beta][a];
    }
    Ttwiddle += Ttwiddle_beta;
  }
  
  // WARNING: C++ calls "beta" the *scale* of the distribution, i.e., the reciprocal of the rate
  auto mu_dist = std::gamma_distribution<double>{
      num_muts_ + 1.0,
      1.0 / Ttwiddle};

  auto old_mu = hky_model_.mu;
  auto new_mu = mu_dist(bitgen_);
  hky_model_.mu = new_mu;
  for (auto beta = 0; beta != evo_.num_partitions(); ++beta) {
    evo_.partition_evo_model[beta] = hky_model_.derive_site_evo_model();
  }
  log_G_ += -(new_mu - old_mu) * Ttwiddle + num_muts_ * std::log(new_mu / old_mu);
  log_other_priors_ += 0.0;
}

auto Run::mpox_hack_moves() -> void {
  // We use a uniform prior for mu and mu*.  The dependence of the posterior on mu and mu*
  // conditional on everything else being fixed is thus:
  //
  // P ~ exp{-mu Ttwiddle + 2 mu* TTwiddle*} * (mu / 3)^(M-M*) * (mu / 3 + 2 mu*)^M*,
  //
  // where
  // - M^beta_ab = num_muts_beta_ab(),  // Number of a->b mutations in partition beta
  // - M = sum_{beta,a,b} M^beta_ab     // Total number of mutations
  // - M* = M^1_CT + M^1_GA             // APOBEC-like mutations in sites with APOBEC context
  // - Ttwiddle^beta_a = sum_{l in beta} nu^(l) T^(l)_a
  //                                    // Site-rate adjusted time spent by sites in partition beta in state a
  // - Ttwiddle = sum_{beta,a} Ttwiddle^beta_a  // Site-rate adjusted total branch length over all (N-pruned) site trees
  // - TTwiddle* = TTwiddle^1_C + Ttwiddle^1_G  // Site-rate adjusted time spent by APOBEC-context sites in C/G states
  //
  // Change variables to mu = mu and rho = mu* / mu.  The Jacobian |d(mu,rho)/d(mu,mu*)| is 1/mu, so
  //
  //  P(mu,rho) ~ e^{-mu (Ttwiddle + 2 rho TTwiddle*)} * (mu / 3)^{M-1} * (1 + 6 rho)^M*.  // M-1 from Jacobian
  //
  //  It follows that
  //
  //         mu|rho ~ Gamma[M, Ttwiddle + 2 rho Ttwiddle*], and
  // (1 + 6 rho)|mu ~ Gamma[M* + 1, (mu / 3) * Ttwiddle*]  (subject to rho >= 0)
  //
  // Hence, we can Gibbs sample both mu and rho.

  // TODO: Always calculate this instead of M_ab
  auto M_beta_ab = calc_num_muts_beta_ab(tree_, evo_);
  auto M_ab = Seq_matrix<int>{0};
  for (auto beta = Partition_index{0}; beta != evo_.num_partitions(); ++beta) {
    for (auto a : k_all_real_seq_letters) {
      for (auto b : k_all_real_seq_letters) {
        M_ab[a][b] += M_beta_ab[beta][a][b];
      }
    }
  }
  auto M = num_muts_;  // == sum_{a,b} M_ab
  auto M_star =
      M_beta_ab[1][Real_seq_letter::C][Real_seq_letter::T] +
      M_beta_ab[1][Real_seq_letter::G][Real_seq_letter::A];

  static auto next_step_output = int64_t{-1};
  if (step_ >= next_step_output) {
    next_step_output = step_ + 1'000'000;
    
    std::cerr << "Mutations by type: \n";
    std::cerr << absl::StreamFormat(
        "- C>T = %d (of which %d = %g%% in APOBEC context)\n",
        M_ab[Real_seq_letter::C][Real_seq_letter::T],
        M_beta_ab[1][Real_seq_letter::C][Real_seq_letter::T],
        M_beta_ab[1][Real_seq_letter::C][Real_seq_letter::T] /
        double(M_ab[Real_seq_letter::C][Real_seq_letter::T]) * 100.0);
    std::cerr << absl::StreamFormat(
        "- G>A = %d (of which %d = %g%% in APOBEC context)\n",
        M_ab[Real_seq_letter::G][Real_seq_letter::A],
        M_beta_ab[1][Real_seq_letter::G][Real_seq_letter::A],
        M_beta_ab[1][Real_seq_letter::G][Real_seq_letter::A] /
        double(M_ab[Real_seq_letter::G][Real_seq_letter::A]) * 100.0);
    std::cerr << absl::StreamFormat("- A>C = %d\n", M_ab[Real_seq_letter::A][Real_seq_letter::C]);
    std::cerr << absl::StreamFormat("- A>G = %d\n", M_ab[Real_seq_letter::A][Real_seq_letter::G]);
    std::cerr << absl::StreamFormat("- A>T = %d\n", M_ab[Real_seq_letter::A][Real_seq_letter::T]);
    std::cerr << absl::StreamFormat("- C>A = %d\n", M_ab[Real_seq_letter::C][Real_seq_letter::A]);
    std::cerr << absl::StreamFormat("- C>G = %d\n", M_ab[Real_seq_letter::C][Real_seq_letter::G]);
    std::cerr << absl::StreamFormat("- G>C = %d\n", M_ab[Real_seq_letter::G][Real_seq_letter::C]);
    std::cerr << absl::StreamFormat("- G>T = %d\n", M_ab[Real_seq_letter::G][Real_seq_letter::T]);
    std::cerr << absl::StreamFormat("- T>A = %d\n", M_ab[Real_seq_letter::T][Real_seq_letter::A]);
    std::cerr << absl::StreamFormat("- T>C = %d\n", M_ab[Real_seq_letter::T][Real_seq_letter::C]);
    std::cerr << absl::StreamFormat("- T>G = %d\n", M_ab[Real_seq_letter::T][Real_seq_letter::G]);
  }
  
  auto Ttwiddle = 0.0;
  for (auto beta = Partition_index{0}; beta != evo_.num_partitions(); ++beta) {
    for (auto a : k_all_real_seq_letters) {
      Ttwiddle += Ttwiddle_beta_a_[beta][a];
    }
  }
  auto Ttwiddle_star =
      Ttwiddle_beta_a_[1][Real_seq_letter::C] +
      Ttwiddle_beta_a_[1][Real_seq_letter::G];

  // Do the above for 10 rounds to pseudo-Gibbs sample the tuple (mu, rho)
  for (auto i = 0; i != 10; ++i) {
    // WARNING: C++ calls "beta" the *scale* of the distribution, i.e., the reciprocal of the rate
    
    auto rho = mpox_mu_star_ / mpox_mu_;
    auto Ttwiddle_eff = Ttwiddle + 2 * rho * Ttwiddle_star;

    // mu
    if (mu_move_enabled_) {
      auto mu_dist = std::gamma_distribution<double>{(double)M, 1.0 / Ttwiddle_eff};

      auto old_mu = mpox_mu_;
      auto new_mu = mu_dist(bitgen_);
      mpox_mu_ = new_mu;
      
      log_G_ += -(new_mu - old_mu) * Ttwiddle_eff + num_muts_ * std::log(new_mu / old_mu);
      log_other_priors_ += 0.0;
    }
      
    // rho
    if (Ttwiddle_star > 0.0) {
      // Use inverse transform sampling on k := 1 + 6 rho:
      //
      // p(k) ~ e^{-(mu Ttwiddle*/3) k} * k^M*, k >= 1
      auto min_Q = 0.0;
      auto a = M_star + 1.0;
      auto z = mpox_mu_*Ttwiddle_star/3;
      auto max_Q = (z < 0.1*a ? 1.0 : boost::math::gamma_q(a, z));   // Wierd edge case in WebAssembly
      auto rand_Q = absl::Uniform(absl::IntervalOpenOpen, bitgen_, min_Q, max_Q);
      auto k = boost::math::gamma_q_inv(M_star + 1, rand_Q) / (mpox_mu_*Ttwiddle_star/3);
      k = std::max(1.0, k);  // Avoid edge cases with round-off errors
      
      auto old_rho = rho;
      auto new_rho = (k - 1) / 6.0;
      mpox_mu_star_ = mpox_mu_ * new_rho;
      
      log_G_ += -2 * mpox_mu_ * (new_rho - old_rho) * Ttwiddle_star
          + M_star * std::log((1 + 6*new_rho) / (1 + 6*old_rho));
      log_other_priors_ += 0.0;
    }
  }

  // Update transition rate matrices
  derive_evo();
}

auto Run::hky_frequencies_move() -> void {
  // Adapting priors & moves from Gire et al 2014 BEAST file
  // - frequencies have a uniform prior
  // - frequencies moved by "delta-exchange", i.e., pick two frequencies, add d to one and -d to the other,
  //   where d is uniformly picked in [0, 0.01] (rejecting proposals of frequencies outside [0,1]).
  //
  // For concreteness, we assume the HKY models across all partitions are linked,
  // though the infrastructure is somewhat there already to allow for unlinked models
  // (e.g., we'd need to keep track of M^beta_ab instead of just M_ab).

  // Make a frequency delta exchange
  auto freq_delta = 0.01;
  auto d = absl::Uniform(bitgen_, 0, freq_delta);
  auto ia = absl::Uniform(absl::IntervalClosedOpen, bitgen_, 0, std::ssize(k_all_real_seq_letters));
  auto ib = ia;
  while (ib == ia) {
    ib = absl::Uniform(absl::IntervalClosedOpen, bitgen_, 0, std::ssize(k_all_real_seq_letters));
  }
  auto move_a = *(k_all_real_seq_letters.begin() + ia);
  auto move_b = *(k_all_real_seq_letters.begin() + ib);

  auto old_hky_model = hky_model_;
  auto new_hky_model = old_hky_model;
  new_hky_model.pi_a[move_a] += d;
  if (new_hky_model.pi_a[move_a] <= 0.0 || new_hky_model.pi_a[move_a] >= 1.0) { return; }
  new_hky_model.pi_a[move_b] -= d;
  if (new_hky_model.pi_a[move_b] <= 0.0 || new_hky_model.pi_a[move_b] >= 1.0) { return; }

  auto old_site_evo = old_hky_model.derive_site_evo_model();  // == evo_.partition_evo_model[0];
  auto new_site_evo = new_hky_model.derive_site_evo_model();

  // Calculate change in log_G and acceptance probability
  auto delta_log_G = 0.0;

  // Branch lengths
  for (auto beta = 0; beta != evo_.num_partitions(); ++beta) {
    for (auto a : k_all_real_seq_letters) {
      delta_log_G -= new_site_evo.mu * (new_site_evo.q_a(a) - old_site_evo.q_a(a)) * Ttwiddle_beta_a_[beta][a];
    }
  }

  // Root priors
  auto force_reject = false;
  auto state_frequencies_at_root = state_frequencies_of_ref_sequence_;
  for (const auto& m : tree_.at_root().mutations) {
    --state_frequencies_at_root[m.from];
    ++state_frequencies_at_root[m.to];
  }
  // TODO: Treat gaps in one go
  for (const auto& [mi_site, mi_from] : tree_.at_root().missations.slow_elements(tree_.ref_sequence)) {
    --state_frequencies_at_root[mi_from];
  }
  for (auto a : k_all_real_seq_letters) {
    if (state_frequencies_at_root[a] > 0) {
      if (new_site_evo.pi_a[a] == 0) { force_reject = true; break; }
      delta_log_G += state_frequencies_at_root[a] * std::log(new_site_evo.pi_a[a] / old_site_evo.pi_a[a]);
    }
  }

  // Mutations
  for (auto a : k_all_real_seq_letters) {
    for (auto b : k_all_real_seq_letters) {
      if (a != b && num_muts_ab_[a][b] > 0) {
        if (num_muts_ab_[a][b] > 0 && new_site_evo.q_ab[a][b] == 0) { force_reject = true; break; }
        delta_log_G += num_muts_ab_[a][b] * std::log(new_site_evo.q_ab[a][b] / old_site_evo.q_ab[a][b]);
      }
    }
  }

  auto log_p_acc = delta_log_G;
  if (not force_reject && (log_p_acc > 0 || absl::Uniform(bitgen_, 0.0, 1.0) < std::exp(log_p_acc))) {
    // Accept
    hky_model_ = new_hky_model;
    for (auto beta = 0; beta != evo_.num_partitions(); ++beta) {
      evo_.partition_evo_model[beta] = new_site_evo;
    }
    log_G_ += delta_log_G;
    log_other_priors_ += 0.0;
  } else {
    // Reject
  }
}

auto Run::hky_kappa_move() -> void {
  // Adapting priors & moves from Gire et al 2014 BEAST file
  // - kappa has a log-normal prior, so log(kappa) has mean 1 and stddev 1.25
  // - kappa move scaled by a uniform factor from 0.75 to 1/0.75
  //    alpha(o -> node) = 1 / (o/0.75 - o*0.75)
  //    alpha(node -> o) = 1 / (node/0.75 - node*0.75)
  //    => alpha(node->o) / alpha(o->node) = o / node
  //
  // For concreteness, we assume the HKY models across all partitions are linked,
  // though the infrastructure is somewhat there already to allow for unlinked models
  // (e.g., we'd need to keep track of M^beta_ab instead of just M_ab).

  // Make a kappa change
  auto scale_factor = 0.75;
  auto scale = absl::Uniform(bitgen_, scale_factor, 1.0 / scale_factor);
  
  auto old_hky_model = hky_model_;
  auto new_hky_model = old_hky_model;
  new_hky_model.kappa *= scale;

  auto old_site_evo = old_hky_model.derive_site_evo_model();  // == evo_.partition_evo_model[0];
  auto new_site_evo = new_hky_model.derive_site_evo_model();

  // log-normal prior for kappa:
  //  pi(kappa) ~ exp(-(log(kappa) - mean_log_kappa)^2 / (2 sigma_log_kappa^2)) / kappa
  auto mean_log_kappa = 1.0;
  auto sigma_log_kappa = 1.25;
  auto alpha_kappa_n_to_o_over_o_to_n = old_hky_model.kappa / new_hky_model.kappa;
  auto prior_new_over_prior_old = std::exp(
      (-std::pow(std::log(new_hky_model.kappa) - mean_log_kappa, 2)
          - (-std::pow(std::log(old_hky_model.kappa) - mean_log_kappa, 2))) / (2 * sigma_log_kappa * sigma_log_kappa)
  ) * old_hky_model.kappa / new_hky_model.kappa;

  // Calculate change in log_G and acceptance probability
  auto delta_log_G = 0.0;

  // Branch lengths
  for (auto beta = 0; beta != evo_.num_partitions(); ++beta) {
    for (auto a : k_all_real_seq_letters) {
      delta_log_G -= new_site_evo.mu * (new_site_evo.q_a(a) - old_site_evo.q_a(a)) * Ttwiddle_beta_a_[beta][a];
    }
  }

  // Mutations
  auto force_reject = false;
  for (auto a : k_all_real_seq_letters) {
    for (auto b : k_all_real_seq_letters) {
      if (a != b && num_muts_ab_[a][b] > 0) {
        if (num_muts_ab_[a][b] > 0 && new_site_evo.q_ab[a][b] == 0) { force_reject = true; break; }
        delta_log_G += num_muts_ab_[a][b] * std::log(new_site_evo.q_ab[a][b] / old_site_evo.q_ab[a][b]);
      }
    }
  }

  auto delta_log_other_priors = std::log(prior_new_over_prior_old);
  auto log_p_acc = delta_log_G + delta_log_other_priors + std::log(alpha_kappa_n_to_o_over_o_to_n);
  if (not force_reject && (log_p_acc > 0 || absl::Uniform(bitgen_, 0.0, 1.0) < std::exp(log_p_acc))) {
    // Accept
    hky_model_ = new_hky_model;
    for (auto beta = 0; beta != evo_.num_partitions(); ++beta) {
      evo_.partition_evo_model[beta] = new_site_evo;
    }
    log_G_ += delta_log_G;
    log_other_priors_ += delta_log_other_priors;
  } else {
    // Reject
  }
}

auto Run::gibbs_sample_all_nus() -> void {
  // We can Gibbs sample a new value of nu_l using a Gamma with
  // shape alpha = M^(l) + alpha and rate beta = mu^(l) * TTwiddle^(l) + alpha

  auto Ttwiddle_l = calc_Ttwiddle_l(tree_, evo_);
  auto M_l = calc_num_muts_l(tree_);
  auto old_sum_nu = 0.0;
  auto new_sum_nu = 0.0;

  for (auto l = Site_index{0}; l != tree_.num_sites(); ++l) {
    // WARNING: C++ calls "beta" the *scale* of the distribution, i.e., the reciprocal of the rate
    auto nu_l_dist = std::gamma_distribution<double>{
      M_l[l] + alpha_,
      1.0 / (evo_.mu_l(l) * Ttwiddle_l[l] + alpha_)};
    
    auto old_nu_l = nu_[l];

    // CAREFUL! When M_l[l] == 0 and alpha is very low, nu_l_dist is very concentrated
    // very close to zero, to the point where random samples can come out as *exactly*
    // zero owing to roundoff.  Those zeros then appear in logs (e.g., below), which
    // wreaks havoc.  If nu_[l] is very close to zero, site l is effectively invariant,
    // and no mutation will be ascribed to it.  Hence, the contribution to log_G is
    // essentially zero, and it doesn't matter what the exact value is.  So we cop out
    // here and put a crazy but nonzero lower bound `eps` on nu_[l].  Technically, we
    // should then adjust the nu prior to distinguish between definite values above `eps`
    // and "anything below `eps`".  Since the distinction has no effect on the sampled
    // trees and remaining parameters, just on the reporting of the posterior values, we
    // purposely keep things simple and don't make the above correction.
    //
    // A more consistent approach would be to allow for some sites to be honest-to-god
    // invariant and *then* demand that the noninvariant sites have a relative site rate
    // above eps.  One day, one day...
    //
    // Thanks to David Pascall at the Cambridge Infectious Diseases centre for
    // reporting and helping debug this issue.
    auto new_nu_l = std::max(1e-50, nu_l_dist(bitgen_));  // eps = 1e-50
    nu_[l] = new_nu_l;
    
    auto log_new_over_old_nu_l = std::log(new_nu_l / old_nu_l);
    log_G_ += -evo_.mu_l(l) * (new_nu_l - old_nu_l) * Ttwiddle_l[l] + M_l[l] * log_new_over_old_nu_l;

    old_sum_nu += old_nu_l;
    new_sum_nu += new_nu_l;
    log_other_priors_ += (alpha_ - 1) * log_new_over_old_nu_l;
  }

  log_other_priors_ += -alpha_ * (new_sum_nu - old_sum_nu);

  evo_.nu_l = nu_;  // TODO: No need to store this twice & copy here
  Ttwiddle_beta_a_ = calc_cur_Ttwiddle_beta_a();
}

static auto calc_log_p_alpha(
    const Global_evo_model& evo,
    double alpha,
    const Node_vector<double>& Ttwiddle_l,
    const Node_vector<int>& M_l)
    -> double {
  
  // Excluding the prior pi(alpha), we have,
  //
  //  log p(alpha) = \sum_{l=0}^L (loggamma(M^(l) + alpha) - (M^(l) + alpha) log(mu^(l) Ttwiddle^(l) + alpha))
  //                          - L (loggamma(        alpha) -          alpha  log(                      alpha))

  auto L = std::ssize(M_l);
  auto result = 0.0;
  auto num_sites_with_muts = 0;
  for (auto l = Site_index{0}; l != L; ++l) {
    if (M_l[l] > 0) {
      ++num_sites_with_muts;
      result += std::lgamma(M_l[l] + alpha);
    }
    result -= (M_l[l] + alpha) * std::log(evo.mu_l(l) * Ttwiddle_l[l] + alpha);
  }
  result -= num_sites_with_muts * std::lgamma(alpha) - L * alpha * std::log(alpha);
  return result;
}

auto Run::alpha_moves() -> void {
  auto Ttwiddle_l = calc_Ttwiddle_l(tree_, evo_);
  auto M_l = calc_num_muts_l(tree_);

  // exponential prior for alpha:
  //  pi(alpha) ~ exp(-alpha / mean_alpha)
  auto mean_alpha = 1.0;
  
  auto alpha_before_moves = alpha_;
  auto cur_log_p_alpha = calc_log_p_alpha(evo_, alpha_, Ttwiddle_l, M_l);
  for (auto sub_move = 0; sub_move != 10; ++sub_move) {
    // We follow LeMieux et al 2021 in using an exponential prior with mean 1.0
    // and scaling moves with scale factor 0.90 (see hky_move kappa moves for a detailed description)
    
    auto scale_factor = 0.90;
    auto old_alpha = alpha_;
    auto scale = absl::Uniform(bitgen_, scale_factor, 1 / scale_factor);
    auto new_alpha = scale * old_alpha;
    auto alpha_n_to_o_over_o_to_n = old_alpha / new_alpha;
    
    auto log_prior_new_over_prior_old = -(new_alpha - old_alpha) / mean_alpha;
    
    // Calculate change in `log p(alpha)` and acceptance probability
    auto new_log_p_alpha = calc_log_p_alpha(evo_, new_alpha, Ttwiddle_l, M_l);
    
    auto delta_log_posterior = log_prior_new_over_prior_old + new_log_p_alpha - cur_log_p_alpha;
    auto log_metropolis = delta_log_posterior + std::log(alpha_n_to_o_over_o_to_n);
    if (log_metropolis > 0 || absl::Uniform(bitgen_, 0.0, 1.0) < std::exp(log_metropolis)) {
      // Accept (no change to log_G)
      alpha_ = new_alpha;
      cur_log_p_alpha = new_log_p_alpha;
    }
  }
  
  // See calc_cur_log_other_priors() for discussion on alpha & nu_l priors
  auto L = tree_.num_sites();
  auto sum_nu = 0.0;
  auto sum_log_nu = 0.0;
  for (auto nu_l : nu_) {
    DCHECK_NE(nu_l, 0.0);
    sum_nu += nu_l;
    sum_log_nu += std::log(nu_l);
  }
  log_other_priors_ += 0.0
      + -(alpha_ - alpha_before_moves) / mean_alpha
      + L * (alpha_ * std::log(alpha_) - alpha_before_moves * std::log(alpha_before_moves))
      - L * (std::lgamma(alpha_) - std::lgamma(alpha_before_moves))
      + (alpha_ - alpha_before_moves) * sum_log_nu
      - (alpha_ - alpha_before_moves) * sum_nu;
  
  // IMPORTANT: Since alpha has changed, all the values of nu_l need to be Gibbs sampled!
  gibbs_sample_all_nus();
}

auto Run::pop_size_move() -> void {
  CHECK(typeid(*pop_model_) == typeid(Exp_pop_model));
  auto exp_pop_model = std::static_pointer_cast<const Exp_pop_model>(pop_model_);
    
  // We follow LeMieux et al (2021) and use
  //
  // - a 1/x prior on the pop size at time 0  (this makes the choice of t = 0 irrelevant)
  // - a scale operator with scale factor [0.75, 1/0.75] to change the pop size
    
  auto scale_factor = 0.75;
  auto old_n0 = exp_pop_model->pop_at_t0();
  auto scale = absl::Uniform(bitgen_, scale_factor, 1 / scale_factor);
  auto new_n0 = scale * old_n0;
  auto alpha_n_to_o_over_o_to_n = old_n0 / new_n0;
    
  // 1/x prior for n0:
  //  pi(n0) ~ 1/n0  => pi(new_n0)/pi(old_n0) = 1/scale
  auto log_prior_new_over_prior_old = -std::log(scale);
    
  auto old_log_coal = log_coalescent_prior_;
  auto old_pop_model = exp_pop_model;
  pop_model_ = std::make_shared<Exp_pop_model>(old_pop_model->t0(), new_n0, old_pop_model->growth_rate());
  coalescent_prior_.pop_model_changed(pop_model_);
  auto new_log_coal = calc_cur_log_coalescent_prior();
    
  auto log_metropolis =
      (new_log_coal - old_log_coal) + log_prior_new_over_prior_old + std::log(alpha_n_to_o_over_o_to_n);
  if (log_metropolis > 0 || absl::Uniform(bitgen_, 0.0, 1.0) < std::exp(log_metropolis)) {
    // Accept
    log_coalescent_prior_ = new_log_coal;
    log_other_priors_ += log_prior_new_over_prior_old;
  } else {
    // Reject
    pop_model_ = old_pop_model;
    coalescent_prior_.pop_model_changed(pop_model_);
  }
}

auto Run::pop_growth_rate_move() -> void {
  CHECK(typeid(*pop_model_) == typeid(Exp_pop_model));
  auto exp_pop_model = std::static_pointer_cast<const Exp_pop_model>(pop_model_);
    
  // We follow LeMieux et al (2021) and use
  //
  // - a Laplace prior on the growth rate, with mu 0.001/365 and scale 30.701135/365
  //     pdf(g; mu, scale) = 1/(2*scale) exp(-|g - mu| / scale)
  // - a random walk operator on the growth rate that adds a uniform perturbation in [-delta, +delta],
  //   with delta = 1.0 / 365.0
    
  auto window_size = 1.0 / 365.0;
  auto mu = 0.001 / 365.0;
  auto scale = 30.701135 / 365.0;
    
  auto delta = absl::Uniform(bitgen_, -window_size, +window_size);
  auto old_g = exp_pop_model->growth_rate();
  auto new_g = old_g + delta;
    
  // Laplace prior for g
  //   pi(g) ~ exp(-|g - mu| / b)
  auto log_prior_new_over_prior_old = (std::abs(old_g - mu) - std::abs(new_g - mu)) / scale;
    
  auto old_log_coal = log_coalescent_prior_;
  auto old_pop_model = exp_pop_model;
  pop_model_ = std::make_shared<Exp_pop_model>(old_pop_model->t0(), old_pop_model->pop_at_t0(), new_g);
  coalescent_prior_.pop_model_changed(pop_model_);
  auto new_log_coal = calc_cur_log_coalescent_prior();
    
  auto log_metropolis = (new_log_coal - old_log_coal) + log_prior_new_over_prior_old;
  if (log_metropolis > 0 || absl::Uniform(bitgen_, 0.0, 1.0) < std::exp(log_metropolis)) {
    // Accept
    log_coalescent_prior_ = new_log_coal;
    log_other_priors_ += log_prior_new_over_prior_old;
  } else {
    // Reject
    pop_model_ = old_pop_model;
    coalescent_prior_.pop_model_changed(pop_model_);
  }
}

auto Run::skygrid_tau_move() -> void {
  CHECK(typeid(*pop_model_) == typeid(Skygrid_pop_model));
  auto skygrid_pop_model = std::static_pointer_cast<const Skygrid_pop_model>(pop_model_);

  // If all other parameters are fixed, tau follows a Gamma distribution:
  //
  // P(tau) \propto tau^(alpha-1) exp(-beta tau)
  //                * tau^(M/2) exp(-tau sum_squared_delta_gammas/2)
  //
  // => tau ~ Gamma(alpha + M/2, beta + sum_squared_delta_gammas/2)
  //
  // We can Gibbs sample that directly

  auto sum_squared_delta_gammas = 0.0;
  for (auto k = 0; k != skygrid_pop_model->M(); ++k) {
    auto delta_gamma = skygrid_pop_model->gamma(k+1) - skygrid_pop_model->gamma(k);
    sum_squared_delta_gammas += delta_gamma*delta_gamma;
  }

  auto prior_tau_alpha = 0.001;  // Gill et al 2012, Eq. 15
  auto prior_tau_beta = 0.001;
  auto M = skygrid_pop_model->M();
  
  auto post_alpha = prior_tau_alpha + 0.5 * M;
  auto post_beta = prior_tau_beta + 0.5 * sum_squared_delta_gammas;

  // WARNING: C++ calls "beta" the *scale* of the distribution, i.e., the reciprocal of the rate
  auto tau_dist = std::gamma_distribution<double>{
      post_alpha,
      1.0 / post_beta};
  
  auto old_tau = skygrid_tau_;
  auto new_tau = tau_dist(bitgen_);
  skygrid_tau_ = new_tau;
  
  log_other_priors_ +=
      (post_alpha - 1) * std::log(new_tau/old_tau)
      - post_beta * (new_tau - old_tau);
}

auto Run::skygrid_gammas_move() -> void {
  CHECK(typeid(*pop_model_) == typeid(Skygrid_pop_model));
  auto skygrid_pop_model = std::static_pointer_cast<const Skygrid_pop_model>(pop_model_);

  // Here, the suggested move in Gill et al 2012 is super-smart but super-complicated.
  // For now, we just displace one of the gamma_k's by a delta ~ N(0, 1/tau).
  // Smarter moves might involve working in Fourier (DCT) space (e.g., the zeroth mode
  // displaces all gamma_k's by the same amount, which costs nothing w.r.t to the GMRF prior),
  // or making concerted displacements (e.g., adding delta to gamma_k for k_min <= k <= k_max)

  auto M = skygrid_pop_model->M();
  auto k = absl::Uniform(absl::IntervalClosedClosed, bitgen_, 0, M);
  auto delta_sd = std::sqrt(1/skygrid_tau_);
  
  auto delta = absl::Gaussian(bitgen_, 0.0, delta_sd);

  auto delta_sum_squared_delta_gammas = 0.0;
  if (k > 0) {
    auto old_delta_gamma = skygrid_pop_model->gamma(k) - skygrid_pop_model->gamma(k-1);
    auto new_delta_gamma = old_delta_gamma + delta;
    
    delta_sum_squared_delta_gammas +=
        new_delta_gamma*new_delta_gamma - old_delta_gamma*old_delta_gamma;
  }
  if (k < M) {
    auto old_delta_gamma = skygrid_pop_model->gamma(k+1) - skygrid_pop_model->gamma(k);
    auto new_delta_gamma = old_delta_gamma - delta;
    
    delta_sum_squared_delta_gammas +=
        new_delta_gamma*new_delta_gamma - old_delta_gamma*old_delta_gamma;
  }

  auto log_prior_new_over_prior_old =
      - (skygrid_tau_/2.0) * delta_sum_squared_delta_gammas;
  auto alpha_n_to_o_over_o_to_n = 1.0;  // Symmetric move
  
  auto old_log_coal = log_coalescent_prior_;
  auto old_pop_model = skygrid_pop_model;
  
  auto new_gamma = std::vector<double>{skygrid_pop_model->gamma()};
  new_gamma.at(k) += delta;
  pop_model_ = std::make_shared<Skygrid_pop_model>(old_pop_model->x(), new_gamma, old_pop_model->type());
  coalescent_prior_.pop_model_changed(pop_model_);
  auto new_log_coal = calc_cur_log_coalescent_prior();
    
  auto log_metropolis =
      (new_log_coal - old_log_coal) + log_prior_new_over_prior_old + std::log(alpha_n_to_o_over_o_to_n);
  if (log_metropolis > 0 || absl::Uniform(bitgen_, 0.0, 1.0) < std::exp(log_metropolis)) {
    // Accept
    log_coalescent_prior_ = new_log_coal;
    log_other_priors_ += log_prior_new_over_prior_old;
  } else {
    // Reject
    pop_model_ = old_pop_model;
    coalescent_prior_.pop_model_changed(pop_model_);
  }
}

auto Run::save_original_sequences() -> void {
  if (estd::is_debug_enabled) {
    // Horribly inefficient, but this is only used for paranoid checks in Debug builds

    // Make a copy of the original ref sequence; all sequences are encoded as overlays atop this copy
    original_ref_sequence_ = tree_.ref_sequence;
    
    for (const auto& node : index_order_traversal(tree_)) {
      original_sequences_.push_back(view_of_sequence_at(tree_, node).rebase_to(original_ref_sequence_));
      auto missing_sites = Interval_set<>{};
      missing_sites = reconstruct_missing_sites_at(tree_, node);
      original_missing_sites_all_.push_back(std::move(missing_sites));
    }
  }
}

auto Run::assert_cur_tip_sequences_compatible_with_original_ones() -> void {
  assert_tip_sequences_compatible_with_original_ones(tree_, original_sequences_, original_missing_sites_all_);
}

}  // namespace delphy
