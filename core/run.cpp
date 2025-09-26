#include "run.h"

#include "distributions.h"
#include "tree_partitioning.h"
#include "phylo_tree_calc.h"
#include "dates.h"

#include <numbers>
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
      skygrid_tau_prior_alpha_{0.001},  // Prior, Gill et al 2012 discussion below Eq. 16
      skygrid_tau_prior_beta_{0.001},
      skygrid_low_gamma_barrier_loc_{std::log(1.0)},  // Pop sizes below 1 day are penalized
      skygrid_low_gamma_barrier_scale_{-std::log(0.70)},  // by 1 nat for every 30% drop
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

  const Pop_model& raw_pop_model = *pop_model_;
  if (typeid(raw_pop_model) == typeid(Exp_pop_model)) {
    const auto& exp_pop_model = static_cast<const Exp_pop_model&>(raw_pop_model);
    
    // pop_n0 - 1/x prior
    log_prior -= std::log(exp_pop_model.pop_at_t0());
    
    // pop_g - Laplace prior on the growth rate, with mu 0.001/365 and scale 30.701135/365
    //      pdf(g; mu, scale) = 1/(2*scale) exp(-|g - mu| / scale)
    const auto mu_g = 0.001 / 365.0;
    const auto scale_g = 30.701135 / 365.0;
    log_prior += -std::abs(exp_pop_model.growth_rate() - mu_g) / scale_g - std::log(2 * scale_g);
    
  } else if (typeid(raw_pop_model) == typeid(Skygrid_pop_model)) {
    const auto& skygrid_pop_model = static_cast<const Skygrid_pop_model&>(raw_pop_model);

    const auto tau = skygrid_tau_;
    const auto log_tau = std::log(tau);
    const auto M = skygrid_pop_model.M();
    
    if (skygrid_tau_move_enabled_) {
      // tau - gamma prior (Gill et al 2012, Eq. 15)
      log_prior +=
          (skygrid_tau_prior_alpha_ - 1) * log_tau
          - skygrid_tau_prior_beta_ * tau;
    }
    
    // gamma - GMRF prior (Gill et al 2012, Eq. 13)
    //
    // Since time intervals are in principle not all the same, we really should
    // be using something like a diffusion law to condition how much consecutive
    // gamma's should vary, e.g., gamma_{i+1} - gamma_i ~ N(0, 2 D (x_{i+1} - x_i)).
    // We don't do that to not deviate from the standard Skygrid parametrization;
    // assuming equal time intervals of width dt = x_1 - x_0, we'd have
    //
    //    2 D dt = 1/tau => D = 1/(2 tau dt).

    for (auto k = 1; k <= M; ++k) {
      auto delta_gamma = skygrid_pop_model.gamma(k) - skygrid_pop_model.gamma(k-1);
      log_prior +=
          0.5 * log_tau
          - 0.5 * std::pow(delta_gamma, 2) * tau;
    }

    if (skygrid_low_gamma_barrier_enabled_) {
      // Strong lower bound on gamma_k's
      //
      // Anything implying a coalescent timescale below N(t) = e^{skygrid_gamma_min} (with
      // default value of 1 day = the resolution of tip timestamps) is probably a numerical
      // instability.  These instabilities tend to drive the population curve to extremes that
      // our moveset finds it hard to recover from.  In realistic datasets, this prior should
      // almost always be exactly zero.
      //
      // The barrier on low gamma_k's is quadratic: for each gamma_k, dropping below
      // skygrid_gamma_min by skygrid_delta_gamma reduces the log-posterior by 1.
      //
      for (auto k = 0; k <= M; ++k) {
        if (skygrid_pop_model.gamma(k) < skygrid_low_gamma_barrier_loc_) {
          auto excess = skygrid_low_gamma_barrier_loc_ - skygrid_pop_model.gamma(k);
          log_prior -= std::pow(excess / skygrid_low_gamma_barrier_scale_, 2);
        }
      }
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
  const Pop_model& raw_pop_model = *pop_model_;
  
  if (typeid(raw_pop_model) == typeid(Exp_pop_model)) {
    for (auto i = 0; i != 50; ++i) {
      if (final_pop_size_move_enabled_) {
        pop_size_move();
        check_derived_quantities();
      }
      if (pop_growth_rate_move_enabled_) {
        pop_growth_rate_move();
        check_derived_quantities();
      }
    }
    
  } else if (typeid(raw_pop_model) == typeid(Skygrid_pop_model)) {
    for (auto i = 0; i != 5; ++i) {  // With HMC, each move likely samples each quantity quite well already
      if (skygrid_tau_move_enabled_) {
        skygrid_tau_move();
        check_derived_quantities();
      }
      
      skygrid_gammas_zero_mode_gibbs_move();
      check_derived_quantities();
      
      skygrid_gammas_hmc_move();
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
  const Pop_model& raw_pop_model = *pop_model_;
  CHECK(typeid(raw_pop_model) == typeid(Exp_pop_model));
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
  const Pop_model& raw_pop_model = *pop_model_;
  CHECK(typeid(raw_pop_model) == typeid(Exp_pop_model));
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
  const Pop_model& raw_pop_model = *pop_model_;
  CHECK(typeid(raw_pop_model) == typeid(Skygrid_pop_model));
  auto skygrid_pop_model = std::static_pointer_cast<const Skygrid_pop_model>(pop_model_);

  // If all other parameters are fixed, tau follows a Gamma distribution:
  //
  // P(tau) \propto tau^(alpha-1) exp(-beta tau)
  //                * tau^(M/2) exp(-tau sum_squared_delta_gammas/2)
  //
  // => tau ~ Gamma(alpha + M/2, beta + sum_squared_delta_gammas/2)
  //
  // We can Gibbs sample that directly

  auto M = skygrid_pop_model->M();
  
  auto sum_squared_delta_gammas = 0.0;
  for (auto k = 1; k <= M; ++k) {
    auto delta_gamma = skygrid_pop_model->gamma(k) - skygrid_pop_model->gamma(k-1);
    sum_squared_delta_gammas += std::pow(delta_gamma, 2);
  }

  auto post_alpha = skygrid_tau_prior_alpha_ + 0.5 * M;
  auto post_beta = skygrid_tau_prior_beta_ + 0.5 * sum_squared_delta_gammas;

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

auto Run::skygrid_gammas_hmc_move() -> void {
  // This move uses Hamiltonian Monte Carlo (HMC) to perform effective collective moves of the
  // Skygrid population curve.  The curve is specified by a set of values gamma_k (0<=k<=M),
  // so that N(x_k) = exp(gamma_k) at predefined "knot" times x_k. At intermediate times, the
  // population curve is interpolated in some reasonable way; and before the first knot and
  // after the last knot, the curve is likewise extrapolated (see pop_model.h for exact
  // details).
  //
  // The dependence of the log-posterior on the {gamma_k} is complicated but tractable:
  //
  //   log P = const                              <-- terms independent of {gamma_k}
  //   
  //           - \sum_{c=0}^{C-1} (1/2) Delta k_c (k_c-1) / N_c      \               |
  //                                                                  +---  -U_coal  |
  //           - \sum_{i inner node} log(N(t_i))                     /               |
  //                                                                                 |
  //           - (tau/2) \sum_{k=1}^M (gamma_k - gamma_{k-1})^2      \               |
  //                                                                  +---  -U_prior |
  //           - \sum_{k=0}^M [[ gamma_k < gamma_min ]]               |              |
  //                 ((gamma_k - gamma_min) / Delta_gamma)^2.        /               |
  //
  // Here, c ranges over the C cells of width Delta in the scalable coalescent,
  //       k_c is the average number of active lineages in cell c,
  //       N_c is the average of N(t) over cell c,
  //       gamma_min is a minimum reasonable value for gamma_k,
  //       Delta_gamma is the scale at which deviations of gamma_k below gamma_min
  //          reduce log P by 1 in log P
  //
  // A straightforward implementation of HMC (e.g., similar to Baele et al 2020
  // doi:10.12688/wellcomeopenres.15770.1, used in BEAST X for Skygrid with HMC) regards the
  // {gamma_k}'s as the "positions" of particles in a dynamical system, augments them with
  // momenta variable {p_k} of some arbitrary mass, then sets up a Hamiltonian and derives
  // equations of motion for those variables.  If every gamma_k gets assigned an equal mass,
  // the free system (U_coal = 0) looks like a vibrating string-of-beads with free ends whose
  // normal modes oscillate at very different frequencies: long modes (which are interesting
  // for N(t)!) are slow, while short modes (which are essentially noise on top of N(t)) are
  // much faster.  This is unfortunate because the HMC time step then needs to be small enough
  // to be a fraction of the period of the fastest mode.  That results in a lot of calculation
  // effort being spent on tracking noisy and irrelevant short-wavelength wiggles on top of
  // N(t).
  //
  // An alternative implementation that first Fourier-transforms the {gamma_k}'s into
  // normal mode amplitudes {A_s} and sets up a dynamical system in _those_ coordinates
  // can be configured with different per-mode "masses" such that all the normal modes
  // of the free system oscillate at the same frequency.  That would allow much larger
  // timesteps in the _free_ system: because the irrelevant short-wavelength vibrations
  // are slowed down while the important long-wavelength vibrations are sped up,
  // computational effort is focused where it matter most.  The details and code for such
  // an HMC setup are in the commit history.
  //
  // Sadly, while the above alternative frequency-space HMC seems ideal in the _free_ system
  // (U_coal = 0), it performs poorly once the effects of the actual coalescent prior are
  // included.  After lots of experimentation, the intuition is that while motion of the
  // {gamma_k}'s is dominated by U_prior in time intervals where there is little information,
  // it is dominated by the effects of U_coal in time intervals where there are lots of
  // coalescences and/or active branches.
  //
  // The above experiments suggested a very effective compromise: a real-space HMC, but with
  // different masses for different gamma_k's.  To get an intuition for what those masses
  // should be in the case of the (standard) staircase Skygrid population curve.  Following
  // the notation from Gill et al 2012, we rewrite log P as follows:
  //
  //   log P = const
  //
  //           - \sum_{k=0}^M SS_k exp(-gamma_k)                 <-- Note k=0 lower bound
  //
  //           - \sum_{k=0}^M c_k gamma_k                        <-- Note k=0 lower bound
  //
  //           - \sum_{k=1}^M (tau/2) (gamma_k - gamma_{k-1})^2. <-- Note k=1 lower bound
  //
  //           (--- terms from low-gamma_k guard omitted ---)
  //
  // Here,
  //
  //   SS_k = int_{y_{k-1}}^{y_k} dt (1/2) k_c(t) (k_c(t) - 1), and
  //    c_k = \sum_{i inner} [[ y_{k-1} < t_i <= y_k ]],
  //
  // while
  //
  //          /   -infinity,         k < 0;      |  y_{k-1} < t <= y_k spans what is
  //    y_k = |   x_k,          0 <= k < M;      |  called "Interval k" in pop_model.h,
  //          \   +infinity,    M <  k.          |  but Interval M is extended to +infinity.
  //
  // Roughly speaking, c_k counts the number of coalescence event in interval k, while SS_k
  // counts something like an ``effective'' elapsed time in interval k.  So we expect the
  // average coalescence timescale (~ N(t)) in interval k to be something like SS_k / c_k
  // (and, in the absence of a GMRF prior, it would be exactly that).  We will see this
  // expectation realized below.
  //
  // Ideally, we'd then maximize `log P` with respect to {gamma_k}, calculate the Hessian at
  // that minimum, and apply a Gaussian approximation to identify the natural normal modes and
  // frequencies of the full system, and then implement HMC in terms of _those_ normal mode
  // amplitudes (indeed, this closely tracks what the ``block update'' move in the original
  // Gill et al 2012 paper does).  Actually realizing this idea is challenging and expensive.
  // Moreover, it would only work exactly for the staircase population curve.  So instead, we
  // recover most of the intuition yielded by the above approach with a much simpler
  // approximation that maintains a diagonal mass matrix.
  //
  // To proceed, consider just the dependence of log P on a particular gamma_k, we get
  //
  //   log P = const
  //           - SS_k exp(-gamma_k)
  //           - c_k gamma_k
  //           - (tau/2) (gamma_k - gamma_{k-1})^2 [[ k > 0 ]]
  //           - (tau/2) (gamma_k - gamma_{k+1})^2 [[ k < M ]].
  //           (--- terms from low-gamma_k guard omitted ---)
  //
  // In the limit of large M, we can suppose that gamma_{k-1} =~ gamma_k =~ gamma_{k+1}, but
  // their small oscillations around these equal values will matter.  In that case, ignoring
  // the low-gamma_k guard and focusing again on a specific k,
  //
  //   log P =~ const
  //           - SS_k exp(-gamma_k)
  //           - c_k gamma_k
  //
  //  while
  //
  //   del(log P)/del gamma_k =~
  //           + SS_k exp(-gamma_k)
  //           - c_k
  //           - tau (gamma_k - gamma_{k-1})  \  these terms are terms
  //           - tau (gamma_k - gamma_{k+1})  /  are essentially zero
  //
  // and
  //
  //   del^2(log P)/del gamma_k^2 =~
  //           - SS_k exp(-gamma_k)
  //           - tau                         <-- independent of {gamma_k} !
  //
  // Neglecting the last two terms in del(log P)/del gamma_k places the maximum of log P at
  // a value gammahat_k satisfying:
  //
  //   exp(-gammahat_k) =~ c_k / SS_k
  //    =>  gammahat_k  =~ -log(c_k / SS_k).
  //
  // (Note: in the staircase model,
  //
  //   N(t in interval k) = exp(gamma_k) =~ SS_k/c_k,
  //
  // as alluded to above).
  //
  // Hence, we have our key result:
  //
  //   del^2(log P)/del gamma_k^2 =~ -(tau + c_k).
  //
  // Compare this to a Gaussian distribution:
  //
  //   log P_G(gamma_k) = -(gamma_k - gammahat_k)^2 / (2 sigma^2),
  //   => del^2(log P_G)/del gamma_k^2 = -1 / sigma^2.
  //
  // Hence, the typical width of the posterior along the gamma_k direction near its minimum is
  // something like 1/sqrt(tau + c_k).  This isn't strictly true owing to the coupling between
  // the neighbouring gamma_k's, and it's very possible that the width of the posterior along
  // other directions is much larger, but this gives us a reasonable starting point with a
  // physically pleasing interpretation:
  //
  //   "In intervals where there is little data (c_k << tau), oscillations of the
  //    population curve are driven by the GMRF prior.  In intervals where there is
  //    a lot of data (c_k >> tau), oscillations of the population curve are driven
  //    by the coalescent prior."
  //
  // As we will see below, setting the mass m_k associated with variable gamma_k to
  //
  //    m_k = tau + c_k
  //
  // results in all oscillations of the dynamical system having very similar angular
  // frequencies close to 1, which permits picking a natural time step and trajectory length.
  // Empirically, the above choice of m_k makes the population curve move slowly in intervals
  // with lots of data while allowing much larger and faster oscillations in intervals with
  // little data.  This behaviour makes the HMC moves quite effective at exploring the
  // posterior distribution of {gamma_k}.
  //
  // -------
  //
  // Now for the nitty-gritty.  We first introduce the (reduced) potential
  // energy U:
  //
  //    U := U_coal + U_prior.
  //
  // We augment this system with momentum variables {p_k} having kinetic energy K:
  //
  //    K := sum_{k=0}^M p_k^2 / (2 m_k),
  //
  // At equilibrium, we have:
  //
  //    p_k ~ N(0, m_k).
  //
  // With these energies in place, we obtain the following Hamiltonian H({gamma_k}, {p_k})
  // and associated equations of motion:
  //
  //    H = U + K,
  //    d gamma_k/dt = + del H / del p_k     = p_k / m_k,
  //    d p_k    /dt = - del H / del gamma_k = - del U / del gamma_k =: f_k;
  //
  // i.e.,
  //
  //    d^2 gamma_k/dt^2 = -(f_k / m_k).
  //
  // We integrate these equations of motion using a slightly nonstandard position Verlet
  // integrator over S steps of size dt:
  //
  //    gamma_k(t + dt/2) = gamma_k(t)        + (dt/2) (p_k(   t    ) / m_k),
  //        p_k(t + dt)   =     p_k(t)        +   dt    f_k(t + dt/2)
  //    gamma_k(t + dt)   = gamma_k(t + dt/2) + (dt/2) (p_k( t + dt ) / m_k).
  //
  // (Baele et al 2020 [and almost everyone else, for that matter] use a standard velocity
  // Verlet integrator, as follows:
  //
  //        p_k(t + dt/2) =     p_k(t)        + (dt/2) f_k(t);
  //    gamma_k(t + dt)   = gamma_k(t)        +   dt  (p_k(t + dt/2) / m_k);
  //        p_k(t + dt)   =     p_k(t + dt/2) + (dt/2) f_k(t + dt).
  //
  // The advantage of using the position Verlet integrator above is there is one fewer force
  // evaluation overall, which is a little bit helpful in the typically short HMC trajectories
  // we use.  Otherwise, velocity and position Verlet are both reversible and have very
  // similar energy conservation properties, with position Verlet possibly a little better
  // for large timestemp [see Tuckerman, Berne and Martyna 1992, doi: 10.1063/1.463137]).
  //
  // The two gaps remaining are:
  // (a) how to calculate f_k?
  // (b) how to choose m_k, dt and S?
  //
  // The second question is easier to answer given the above discussion.  As noted above,
  // near the maximum of log P at {gammahat_k} , we'll have something like
  //
  //   d^2 (gamma_k - gammahat_k) / dt^2 = - (gamma_k - gammahat_k) * (tau + c_k) / m_k.
  //                                       \                                    /
  //                                        +------------ ~ f_k ---------------+
  //
  // Hence, setting
  //
  //   m_k = tau + c_k
  //
  // makes gamma_k behave as a harmonic oscillator about gammahat_k, with angular frequency 1.
  // This is not exactly what happens, but it's close enough to be a useful guide to tuning
  // masses and timesteps.
  //
  // For stability and to compensate for the approximate nature of the above argument, we'd
  // like the position Verlet integrator to trace out a full period in about 100 steps,
  // which suggests using:
  //
  //    dt =? 2 pi / 100.
  //
  // If we ran the position Verlet integrator for exactly one period = 100 steps, then to a
  // first approximation, all {gamma_k}'s would end at the same place where they started,
  // regardless of the starting momenta!  Moreover, while the above estimate for dt should
  // work reasonably well in most cases, we risk getting stuck in a region of parameter space
  // where forces are very high and the timestep needs to be substantially lower to improve
  // integration accuracy.
  //
  // We kill two birds with one stone by implementing the so-called Randomized HMC rule of
  // making the total integration time random with an exponential distribution (Bou-Rabee and
  // Maria Sanz-Serna 2017, doi:10.1214/16-AAP1255).  Concretely, instead of picking dt
  // deterministically, we draw a sample from an exponential distribution.  Then we sample for
  // exactly 25 steps; this is exactly 1/4 of a period if dt assumes its mean value and the
  // approximations above are accurate, where we expect the variables to differ most from
  // their initial values:
  //
  //    dt ~ Expo(100 / 2 pi),
  //     S = 25.
  //
  // Note that because this choice is symmetric between the forward and reverse moves, it does
  // not induce any correction terms in the acceptance probabilities below.
  //
  // Incidentally, note that setting the actual value of the angular frequencies to 1 is
  // irrelevant.  Had we chosen a different value for it and scaled {m_k} and dt accordingly,
  // the generated {gamma_k} trajectories would be unchanged.
  //
  // Now for forces.  We have
  //
  //   f_k := -(del U / del gamma_k),
  //
  // where
  //
  //  U = + \sum_{c=0}^{C-1} (1/2) Delta k_c (k_c-1) / N_c
  //      + \sum_{i inner node} log(N(t_i))
  //      + (tau/2) \sum_{k=1}^M (gamma_k - gamma_{k-1})^2
  //      + \sum_{k=0}^M [[ gamma_k < gamma_min ]] ((gamma_k - gamma_min) / Delta_gamma)^2.
  //
  // Hence,
  //
  //   f_k = + \sum_{c=0}^{C-1} (1/2) Delta (k_c (k_c-1) / N_c^2) (del N_c / del gamma_k)
  //           - \sum_{i inner node} [del log(N(t)) / del gamma_k]_{t = t_i}
  //           - tau (  (gamma_k - gamma_{k-1}) [[ k > 0 ]]
  //                  + (gamma_k - gamma_{k+1}) [[ k < M ]] )
  //           - \sum_{k=0}^M [[ gamma_k < gamma_min ]] 2 (gamma_k - gamma_min)/Delta_gamma^2.
  //
  //       = + \sum_{c=0}^{C-1} (1/2) Delta (k_c (k_c-1) / N_c) (del log(N_c) / del gamma_k)
  //           - \sum_{i inner node} [del log(N(t)) / del gamma_k]_{t = t_i}.
  //           - tau (  (gamma_k - gamma_{k-1}) [[ k > 0 ]]
  //                  + (gamma_k - gamma_{k+1}) [[ k < M ]] )
  //           - \sum_{k=0}^M [[ gamma_k < gamma_min ]] 2 (gamma_k - gamma_min)/Delta_gamma^2.
  //
  // In the first sum, only a few terms are non-zero for a given k: these are the cells
  // that overlap with Pop_model::support_of_d_log_N_d_gamma.  The key factors
  // are then given by Pop_model::d_log_int_N_d_gamma
  // and Pop_model::d_log_N_d_gamma, respectively.

  // TODO: The code below could be significantly optimized:
  // - mutating population models in-place to remove allocations
  // - replacing the loop over all tips in the calculation of f_k[k] for each k
  //   by a way to identify the handful of `k`s affected by any one tip
  // - replacing `sum_{i inner} 1/N(t_i)` in the coalescent prior by `sum_c (1/N_c)^{n_c},
  //   where n_c is the number of tips in interval `k`.
  // - merging the position half-steps of consecutive steps
 
  const Pop_model& raw_pop_model = *pop_model_;
  CHECK(typeid(raw_pop_model) == typeid(Skygrid_pop_model));
  auto old_pop_model = std::static_pointer_cast<const Skygrid_pop_model>(pop_model_);
  auto new_pop_model = old_pop_model;

  auto M = old_pop_model->M();

  auto Delta = coalescent_prior_.t_step();
  const auto& k_c = coalescent_prior_.k_bars();
  const auto& N_c = coalescent_prior_.popsize_bars();
  auto C = static_cast<int>(std::ssize(k_c));
  auto tau = skygrid_tau_;

  // real-space "positions" and "momenta"
  auto gamma_k = std::vector<double>(M+1, 0.0);
  auto p_k = std::vector<double>(M+1, 0.0);

  // Count number of coalescences in each Skygrid interval
  // (intervals with more coalescences get assigned a higher mass, which slows them down)
  auto c_k = std::vector<double>(M+1, 0.0);
  for (const auto& node : index_order_traversal(tree_)) {
    if (not tree_.at(node).is_tip()) {
      auto t = tree_.at(node).t;
      auto node_k = std::clamp(old_pop_model->interval_containing_t(t), 0, M);
      c_k[node_k] += 1.0;  // Should be accurate up to 2^52 coalescences...
    }
  }

  // masses of each gamma_k that roughly equalize angular frequencies normal modes to "1"
  auto m_k = std::vector<double>(M+1, 0.0);
  auto inv_m_k = std::vector<double>(M+1, 0.0);       // inverse masses
  for (auto k = 0; k <= M; ++k) {
    m_k[k] = tau + c_k[k];
    inv_m_k[k] = 1.0 / m_k[k];
  }

  // components of forces
  auto f_k = std::vector<double>(M+1, 0.0);
  
  // Helper functions to calculate intermediate quantitites
  auto calc_K = [&]() -> double {
    auto K = 0.0;
    for (auto k = 0; k <= M; ++k) {
      K += 0.5 * std::pow(p_k[k], 2) * inv_m_k[k];
    }
    return K;
  };
  auto calc_U_prior_from_gamma_k_s = [&]() -> double {
    auto U_prior = 0.0;
    for (auto k = 1; k <= M; ++k) {  // Note: k == 0 excluded!
      U_prior += std::pow(gamma_k[k] - gamma_k[k-1], 2);
    }
    U_prior *= 0.5 * tau;

    if (skygrid_low_gamma_barrier_enabled_) {
      // Strong lower bound on gamma_k's: See comments in calc_cur_log_other_priors().
      for (auto k = 0; k <= M; ++k) {
        if (gamma_k[k] < skygrid_low_gamma_barrier_loc_) {
          auto excess = skygrid_low_gamma_barrier_loc_ - gamma_k[k];
          U_prior += std::pow(excess / skygrid_low_gamma_barrier_scale_, 2);
        }
      }
    }
      
    return U_prior;
  };

  auto update_new_pop_model_from_gamma_k_s = [&]() -> void {
    new_pop_model = std::make_shared<Skygrid_pop_model>(
        old_pop_model->x(), gamma_k, old_pop_model->type());
  };

  auto calc_f_k_s_from_gamma_k_s = [&]() -> void {

    auto t_min_coal = coalescent_prior_.cell_lbound(0);
    auto t_max_coal = coalescent_prior_.cell_ubound(C-1);
    
    for (auto k = 0; k <= M; ++k) {
      f_k[k] = 0.0;

      // Careful: support can go down to -inf for k = 0, and up to +inf for k = M
      // (for staircase model, t_min = y_{k-1} and t_max = y_k, defined in the comment above)
      auto [t_min, t_max] = new_pop_model->support_of_d_log_N_d_gamma(k);

      // And careful: support can be fully outside range of coalescent prior cells
      auto c_min = (t_min < t_min_coal) ? 0   : std::clamp(coalescent_prior_.cell_for(t_min), 0, C-1);
      auto c_max = (t_max > t_max_coal) ? C-1 : std::clamp(coalescent_prior_.cell_for(t_max), 0, C-1);
      
      for (auto c = c_min; c <= c_max; ++c) {
        auto a = coalescent_prior_.cell_lbound(c);
        auto b = coalescent_prior_.cell_ubound(c);
        f_k[k] +=
            0.5 * Delta * k_c[c] * (k_c[c] - 1.0) / N_c[c]
            * new_pop_model->d_log_int_N_d_gamma(a, b, k);
      }

      for (const auto& node : index_order_traversal(tree_)) {
        if (not tree_.at(node).is_tip()) {
          f_k[k] -= new_pop_model->d_log_N_d_gamma(tree_.at(node).t, k);
        }
      }

      if (k > 0) {
        f_k[k] -= tau * (gamma_k[k] - gamma_k[k-1]);
      }
      if (k < M) {
        f_k[k] -= tau * (gamma_k[k] - gamma_k[k+1]);
      }

      if (skygrid_low_gamma_barrier_enabled_) {
        // Strong lower bound on gamma_k's: See comments in calc_cur_log_other_priors().
        if (gamma_k[k] < skygrid_low_gamma_barrier_loc_) {
          auto excess = skygrid_low_gamma_barrier_loc_ - gamma_k[k];
          f_k[k] += 2 * excess / std::pow(skygrid_low_gamma_barrier_scale_, 2);
        }
      }
    }
  };

  
  // Start of Leapfrog integrator
  constexpr auto debug_hmc = false;

  // Initial "positions" and "momenta"
  for (auto k = 0; k <= M; ++k) {
    gamma_k[k] = old_pop_model->gamma(k);
    p_k[k] = absl::Gaussian(bitgen_, 0.0, std::sqrt(m_k[k]));
  }
  
  auto old_K = calc_K();
  auto old_U_prior = calc_U_prior_from_gamma_k_s();
  auto old_U_coal = -log_coalescent_prior_;
  auto old_H = old_K + old_U_prior + old_U_coal;

  if (debug_hmc) {
    // Debugging
    std::cerr << absl::StrFormat(
        "Initial:       K = %10.1f, U_prior = %10.1f, U_coal = %10.1f, H = %10.1f\n",
        old_K, old_U_prior, old_U_coal, old_H);
  }

  // Reject outright if initial momenta are outrageously high (see comments in analogous check mid-loop below)
  if (calc_K() > (100.0 * (M+1))) {
    if (debug_hmc) {
      std::cerr << "Rejecting Skygrid HMC because it's initial K is too high: K = " << calc_K() << "...\n";
    }
    pop_model_ = old_pop_model;
    coalescent_prior_.pop_model_changed(old_pop_model);
    return;
  }
  
  // Position Verlet algorithm:
  //
  //    gamma_k(t + dt/2) = gamma_k(t)        + (dt/2) (p_k(   t    ) / m_k),
  //        p_k(t + dt)   =     p_k(t)        +   dt    f_k(t + dt/2)
  //    gamma_k(t + dt)   = gamma_k(t + dt/2) + (dt/2) (p_k( t + dt ) / m_k).
  //
  // For clarity, we don't try to merge the second position update of one step with
  // the first position update of the next step; this isn't the speed bottleneck
  const auto dt = absl::Exponential(bitgen_, 1.0 / (2 * std::numbers::pi / 100.0));

  // RHMC, mean 25 steps
  const auto num_steps = 25;

  for (auto step = 0; step != num_steps; ++step) {

    // "Position" half-step
    for (auto k = 0; k <= M; ++k) {
      gamma_k[k] += 0.5 * dt * p_k[k] * inv_m_k[k];
    }
    update_new_pop_model_from_gamma_k_s();
    coalescent_prior_.pop_model_changed(new_pop_model);  // Updates N_c

    // Force calculation
    calc_f_k_s_from_gamma_k_s();

    if (debug_hmc) {
      // Check force calculation
      auto d_gamma = 1e-4;
      for (auto k = 0; k <= M; ++k) {
        auto old_gamma_k = gamma_k[k];
        
        gamma_k[k] = old_gamma_k + d_gamma;
        update_new_pop_model_from_gamma_k_s();
        coalescent_prior_.pop_model_changed(new_pop_model);
        auto U_plus = -calc_cur_log_coalescent_prior() + calc_U_prior_from_gamma_k_s();

        gamma_k[k] = old_gamma_k - d_gamma;
        update_new_pop_model_from_gamma_k_s();
        coalescent_prior_.pop_model_changed(new_pop_model);
        auto U_minus = -calc_cur_log_coalescent_prior() + calc_U_prior_from_gamma_k_s();

        auto slow_f_k = -(U_plus - U_minus) / (2 * d_gamma);
        CHECK_LE(std::abs(slow_f_k - f_k[k]), std::max(1e-5, 1e-2 * std::max(std::abs(f_k[k]), std::abs(slow_f_k))))
            << k << " " << old_gamma_k << " " << f_k[k] << " " << slow_f_k;

        gamma_k[k] = old_gamma_k;
      }

      // Return to original gamma_k's
      update_new_pop_model_from_gamma_k_s();
      coalescent_prior_.pop_model_changed(new_pop_model);
    }

    // "Momentum" full-step
    for (auto k = 0; k <= M; ++k) {
      p_k[k] += dt * f_k[k];
    }

    // Sanity check: if kinetic energy starts exploding, reject move prematurely
    // (numerics say this is simply going to crash).
    // At equilibrium, K =~ (M+1) / 2.
    // To ensure detailed balance, we also apply this rejection criteria at the very beginning;
    // thus, we _never_ propose trajectories where K > 10 (M+1) at any point.
    if (calc_K() > (100.0 * (M+1))) {
      if (debug_hmc) {
        std::cerr << "Rejecting Skygrid HMC because it's blowing up: K = " << calc_K() << "...\n";
      }
      pop_model_ = old_pop_model;
      coalescent_prior_.pop_model_changed(old_pop_model);
      return;
    }

    // "Position" half-step
    for (auto k = 0; k <= M; ++k) {
      gamma_k[k] += 0.5 * dt * p_k[k] * inv_m_k[k];
    }

    if (debug_hmc) {
      // Debugging (all of the below can go away if debug trace is not needed)
      update_new_pop_model_from_gamma_k_s();
      coalescent_prior_.pop_model_changed(new_pop_model);
      
      auto cur_K = calc_K();
      auto cur_U_prior = calc_U_prior_from_gamma_k_s();
      auto cur_U_coal = -calc_cur_log_coalescent_prior();
      auto cur_H = cur_K + cur_U_prior + cur_U_coal;
      
      std::cerr << absl::StrFormat(
          "Post step %3d: K = %10.1f, U_prior = %10.1f, U_coal = %10.1f, H = %10.1f\n",
          step, cur_K, cur_U_prior, cur_U_coal, cur_H);
    }
  }

  update_new_pop_model_from_gamma_k_s();
  coalescent_prior_.pop_model_changed(new_pop_model);
  
  auto new_K = calc_K();
  auto new_U_prior = calc_U_prior_from_gamma_k_s();
  auto new_U_coal = -calc_cur_log_coalescent_prior();
  auto new_H = new_K + new_U_prior + new_U_coal;

  auto log_metropolis = old_H - new_H;  // HMC => this should be close to 0.0

  // Defend against horrible NaN's and overflow errors
  auto something_blew_up = false;
  for (auto k = 0; k <= M; ++k) {
    if (std::isnan(gamma_k[k])) {
      something_blew_up = true;
    }
  }
  if (std::isnan(log_metropolis)) {
    something_blew_up = true;
  }
  if (something_blew_up && debug_hmc) {
    std::cerr << "About to reject because a NaN appeared somewhere\n";
  }
  
  if (not something_blew_up && (log_metropolis > 0 || absl::Uniform(bitgen_, 0.0, 1.0) < std::exp(log_metropolis))) {
    // Accept
    if (debug_hmc) {
      std::cerr << "Accept!\n";
      for (auto k = 0; k <= M; ++k) {
        std::cerr << absl::StreamFormat("* gamma_k[%2d] = %10.5f  ( at x_k[%2d] = %s)\n",
                                        k, gamma_k[k], k, to_iso_date(new_pop_model->x(k)));
      }
    }
    pop_model_ = new_pop_model;
    log_coalescent_prior_ = -new_U_coal;
    log_other_priors_ += (-new_U_prior) - (-old_U_prior);
  } else {
    // Reject
    if (debug_hmc) {
      std::cerr << "Reject!\n";
    }
    pop_model_ = old_pop_model;
    coalescent_prior_.pop_model_changed(pop_model_);
  }
}

auto Run::skygrid_gammas_zero_mode_gibbs_move() -> void {
  //
  // Amazingly, the overall size of the population is something that can be Gibbs-sampled in
  // the Skygrid model.  First, we decompose {gamma_k} into a baseline `gamma_bar` plus
  // offsets {gamma'_k} which add up to 0:
  //
  //   gamma_bar := 1/(M+1) sum_{k=0}^M gamma_k,  and
  //   gamma'_k  := gamma_k - gamma_bar
  //
  // yields
  //
  //   gamma_k = gamma_bar + gamma'_k   and   sum_{k=0}^M gamma'_k = 0.
  //
  // Clearly, U_prior does not constrain `gamma_bar` at all.  However, U_coal's dependence on
  // `gamma_bar` is particularly simple.  This is best seen with a change of coordinates:
  //
  //    {gamma_k}  ->  {I_bar, gamma'_k | k < M}  (gamma'_M excluded as it's not independent).
  //
  // Here
  //
  //   I_bar := exp(-gamma_bar).
  //
  // We first factor out I_bar from N(t) and N_c, respectively:
  //
  //   Ntwiddle(t) := N(t) I_bar,
  //   Ntwiddle_c  := N_c  I_bar.
  //
  // Note that this makes Ntwiddle(t) and Ntwiddle_c a function purely of {gamma'_k | k < M},
  // and so they are both independent of I_bar.  Conversely, I_bar is by construction
  // independent of {gamma'_k | k < M}.
  //
  // We then write U_coal in terms of I_bar and {gamma'_k | k < M}:
  //
  //                     +---------------------- B ------------------------+
  //                    /                                                   \ |
  //   U_coal = + I_bar \sum_{c=0}^{C-1} (1/2) Delta k_c (k_c-1) / Ntwiddle_c
  //            - N_inner log(I_bar) + \sum_{i inner node} log(Ntwiddle(t_i)).
  //
  // Hence, the dependence of the posterior (~ exp(-U)) on I_bar is almost a straightforward
  // Gamma distribution!  Explicitly,
  //
  //   exp[-U_coal] = I_bar^{N_inner} exp[- B I_bar] * (factors independent of I_bar)
  //
  // That means that I_bar is amenable to Gibbs sampling.  The caveat is that there are
  // Jacobian factors that modify the base result, but not by much.  First, note that an
  // intermediate change from {gamma_k} to gamma_bar & {gamma'_k | k < M} has a constant
  // Jacobian, so that
  //
  //   P({gamma_k}) d{gamma_k} =
  //       (const.) * P({gamma_bar, gamma'_k | k < M}) d{gamma_bar, gamma'_k | k < M}.
  //
  // Hence, the only difficulty comes from switching from gamma_bar to I_bar.  We have
  //
  //           I_bar = exp[ -gamma_bar ]
  //   =>     dI_bar = - exp[ -gamma_bar ] dgamma_bar
  //   => dgamma_bar = - dI_bar / I_bar
  //
  // Thus, when switching from gamma_bar to I_bar, apart from a constant factor, we gain
  // an extra factor of 1/I_bar in the posterior:
  //
  //   P(I_bar) dI_bar = (const.) * I_bar^{N_inner-1} exp[- B I_bar] dI_bar.
  //
  // Hence, we can Gibbs sample I_bar from Gamma(alpha = N_inner, lambda = B).
  
  const Pop_model& raw_pop_model = *pop_model_;
  CHECK(typeid(raw_pop_model) == typeid(Skygrid_pop_model));
  auto old_pop_model = std::static_pointer_cast<const Skygrid_pop_model>(pop_model_);

  auto M = old_pop_model->M();

  auto Delta = coalescent_prior_.t_step();
  const auto& k_c = coalescent_prior_.k_bars();
  const auto& N_c = coalescent_prior_.popsize_bars();
  auto C = std::ssize(k_c);
  
  auto gamma_bar = 0.0;
  for (auto k = 0; k <= M; ++k) {
    gamma_bar += old_pop_model->gamma(k);
  }
  gamma_bar /= (M+1);
  auto I_bar = std::exp(-gamma_bar);
  
  auto N_inner = (tree_.size() - 1) / 2;
  auto B = 0.0;
  for (auto c = 0; c < C; ++c) {
    B += 0.5 * Delta * k_c[c] * (k_c[c] - 1.0) / N_c[c];
  }
  B /= I_bar;

  // WARNING: C++ calls "beta" the *scale* of the distribution, i.e., 1/rate
  auto I_bar_dist = std::gamma_distribution<double>{
    static_cast<double>(N_inner),
    1.0 / B};
  
  auto old_I_bar = I_bar;
  auto new_I_bar = I_bar_dist(bitgen_);
  auto delta_gamma_bar = std::log(old_I_bar / new_I_bar); // -log(new_I_bar)-(-log(old_I_bar))

  auto new_gamma_k = std::vector<double>(M+1, 0.0);
  for (auto k = 0; k <= M; ++k) {
    new_gamma_k[k] = old_pop_model->gamma(k) + delta_gamma_bar;
  }
  
  // Strong lower bound on gamma_k's: See comments in calc_cur_log_other_priors().
  auto old_log_prior_bound = 0.0;
  auto new_log_prior_bound = 0.0;

  if (skygrid_low_gamma_barrier_enabled_) {
    for (auto k = 0; k <= M; ++k) {
      if (old_pop_model->gamma(k) < skygrid_low_gamma_barrier_loc_) {
        auto excess = skygrid_low_gamma_barrier_loc_ - old_pop_model->gamma(k);
        old_log_prior_bound -= std::pow(excess / skygrid_low_gamma_barrier_scale_, 2);
      }
      if (new_gamma_k[k] < skygrid_low_gamma_barrier_loc_) {
        auto excess = skygrid_low_gamma_barrier_loc_ - new_gamma_k[k];
        new_log_prior_bound -= std::pow(excess / skygrid_low_gamma_barrier_scale_, 2);
      }
    }
  }

  auto delta_log_posterior = new_log_prior_bound - old_log_prior_bound;  // Otherwise, this is a Gibbs move

  // Reject if anything becomes a NaN
  auto something_blew_up = false;
  for (auto k = 0; k <= M; ++k) {
    if (std::isnan(new_gamma_k[k])) {
      something_blew_up = true;
    }
  }
  if (std::isnan(delta_log_posterior)) {
    something_blew_up = true;
  }
  if (something_blew_up) {
    if (false) {
      std::cerr << "Rejecting zero-mode Gibbs move because a NaN appeared somewhere\n";
    }
    return;
  }
  
  auto log_metropolis = delta_log_posterior;
  if (log_metropolis >= 0 || absl::Uniform(bitgen_, 0.0, 1.0) < std::exp(log_metropolis)) {
    // Accept
    pop_model_ = std::make_shared<Skygrid_pop_model>(old_pop_model->x(), new_gamma_k, old_pop_model->type());
    coalescent_prior_.pop_model_changed(pop_model_);
    
    log_coalescent_prior_ = calc_cur_log_coalescent_prior();
    log_other_priors_ += new_log_prior_bound - old_log_prior_bound;
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
