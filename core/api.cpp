#include "api.h"

#include <ranges>

#include "run.h"
#include "sequence.h"
#include "dates.h"
#include "site_states_tree_prober.h"
#include "ancestral_tree_prober.h"
#include "phylo_tree_calc.h"

namespace delphy {

static auto to_real_seq_letter(api::RealSeqLetter c) -> Real_seq_letter {
  switch (c) {
    case api::RealSeqLetter_A: return Real_seq_letter::A;
    case api::RealSeqLetter_C: return Real_seq_letter::C;
    case api::RealSeqLetter_G: return Real_seq_letter::G;
    case api::RealSeqLetter_T: return Real_seq_letter::T;
    default: throw std::invalid_argument("c");
  }
}

static auto to_api_real_seq_letter(Real_seq_letter c) -> api::RealSeqLetter {
  switch (c) {
    case Real_seq_letter::A: return api::RealSeqLetter_A;
    case Real_seq_letter::C: return api::RealSeqLetter_C;
    case Real_seq_letter::G: return api::RealSeqLetter_G;
    case Real_seq_letter::T: return api::RealSeqLetter_T;
    default: throw std::invalid_argument("c");
  }
}

auto phylo_tree_to_api_tree(const Phylo_tree& tree) -> flatbuffers::DetachedBuffer {
  auto num_muts = 0;
  auto num_missation_intervals = 0;
  for (const auto& node : index_order_traversal(tree)) {
    num_muts += std::ssize(tree.at(node).mutations);
    num_missation_intervals += tree.at(node).missations.intervals.num_intervals();
  }

  // Make initial buffer size large enough to avoid reallocations later
  auto fbb = flatbuffers::FlatBufferBuilder{
    1024
    + sizeof(api::Node) * std::ssize(tree)
    + sizeof(api::Mutation) * num_muts
    + sizeof(api::MissationInterval) * num_missation_intervals
    + sizeof(api::RealSeqLetter) * tree.num_sites()
  };
  
  auto api_ref_seq = fbb.CreateVector(tree.ref_sequence).o;  // .o = Work around small typing glitch in gen'd code

  auto out_mutation = (api::Mutation*){nullptr};
  auto api_mutations = fbb.CreateUninitializedVectorOfStructs(num_muts, &out_mutation);

  auto out_missation_interval = (api::MissationInterval*){nullptr};
  auto api_missation_intervals = fbb.CreateUninitializedVectorOfStructs(
      num_missation_intervals, &out_missation_interval);

  auto out_node = (api::Node*){nullptr};
  auto api_nodes = fbb.CreateUninitializedVectorOfStructs(std::ssize(tree), &out_node);

  for (auto node = Node_index{0}; node != std::ssize(tree); ++node) {
    *out_node = api::Node{
      tree.at(node).parent,
      tree.at(node).is_tip() ? -1 : tree.at(node).left_child(),
      tree.at(node).is_tip() ? -1 : tree.at(node).right_child(),
      static_cast<float>(tree.at(node).t)};
    ++out_node;

    for (const auto& m : tree.at(node).mutations) {
      *out_mutation = api::Mutation{
        node,
        m.site,
        to_api_real_seq_letter(m.from),
        to_api_real_seq_letter(m.to),
        static_cast<float>(m.t)};
      ++out_mutation;
    }

    auto by_start = [](const Site_interval& r) { return r.first; };
    CHECK(std::ranges::is_sorted(tree.at(node).missations.intervals, {}, by_start));
    for (const auto& [start, end] : tree.at(node).missations.intervals) {
      *out_missation_interval = api::MissationInterval{
        node,
        start,
        end
      };
      ++out_missation_interval;
    }
  }

  auto api_tree = api::CreateTree(fbb, api_nodes, api_mutations, api_missation_intervals, api_ref_seq, tree.root);

  fbb.FinishSizePrefixed(api_tree);

  return fbb.Release();
}

auto phylo_tree_to_api_tree_info(const Phylo_tree& tree) -> flatbuffers::DetachedBuffer {
  auto fbb = flatbuffers::FlatBufferBuilder{};

  auto node_infos = std::vector<flatbuffers::Offset<api::NodeInfo>>{};

  for (const auto& node : index_order_traversal(tree)) {
    auto api_name = fbb.CreateString(tree.at(node).name);
    node_infos.push_back(api::CreateNodeInfo(fbb, api_name));
  }

  auto api_node_infos = fbb.CreateVector(node_infos);
  auto api_tree_info = api::CreateTreeInfo(fbb, api_node_infos);

  fbb.FinishSizePrefixed(api_tree_info);
  return fbb.Release();
}

auto api_tree_and_tree_info_to_phylo_tree(const uint8_t* tree_fb, const uint8_t* tree_info_fb) -> Phylo_tree {
  auto api_tree = flatbuffers::GetSizePrefixedRoot<api::Tree>(tree_fb);
  auto api_tree_info = flatbuffers::GetSizePrefixedRoot<api::TreeInfo>(tree_info_fb);

  auto num_nodes = static_cast<Node_index>(api_tree->nodes()->size());
  CHECK_EQ(num_nodes, std::ssize(*api_tree_info->node_infos()));
  
  auto tree = Phylo_tree{num_nodes};
  
  tree.ref_sequence.reserve(api_tree->ref_seq()->size());
  for (auto c : *api_tree->ref_seq()) {
    tree.ref_sequence.push_back(to_real_seq_letter(api::EnumValuesRealSeqLetter()[c]));
  }

  for (auto node = Node_index{0}; node != num_nodes; ++node) {
    auto api_node = api_tree->nodes()->Get(node);
    auto api_node_info = api_tree_info->node_infos()->Get(node);
    tree.at(node).name = api_node_info->name()->str();
    tree.at(node).parent = api_node->parent();
    if (api_node->left_child() == k_no_node && api_node->right_child() == k_no_node) {
      tree.at(node).children = {};
    } else {
      tree.at(node).children = {api_node->left_child(), api_node->right_child()};
    }
    tree.at(node).t = api_node->t();
  }

  for (const auto& api_mutation : *api_tree->mutations()) {
    auto branch = api_mutation->branch();
    tree.at(branch).mutations.emplace_back(to_real_seq_letter(api_mutation->from()),
                                           api_mutation->site(),
                                           to_real_seq_letter(api_mutation->to()),
                                           api_mutation->t());
  }

  for (const auto& api_missation_interval : *api_tree->missation_intervals()) {
    auto branch = api_missation_interval->branch();
    tree.at(branch).missations.intervals.insert({
        api_missation_interval->start_site(),
        api_missation_interval->end_site()});
  }

  tree.root = api_tree->root_node();

  fix_up_missations(tree);  // Fill in from_states wherever they differ from the reference sequence

  return tree;
}

auto run_to_api_params(const Run& run) -> flatbuffers::DetachedBuffer {
  const auto& tree = run.tree();
  auto fbb = flatbuffers::FlatBufferBuilder{
    512    // Large enough initial buffer to avoid reallocations later
    + sizeof(api::Params)
    + sizeof(double) * tree.num_sites()};
  
  auto params_builder = api::ParamsBuilder{fbb};

  params_builder.add_step(run.step());
  params_builder.add_num_local_moves_per_global_move(run.local_moves_per_global_move());
  params_builder.add_num_parts(run.num_parts());
  params_builder.add_mu(run.mu());
  params_builder.add_alpha(run.alpha());
  if (not std::ranges::all_of(run.nu(), [](auto nu_l) { return nu_l == 1.0; })) {
    auto api_nu = fbb.CreateVector(run.nu());
    params_builder.add_nu(api_nu);
  }
  params_builder.add_hky_kappa(run.hky_kappa());
  params_builder.add_hky_pi_A(run.hky_pi()[Real_seq_letter::A]);
  params_builder.add_hky_pi_C(run.hky_pi()[Real_seq_letter::C]);
  params_builder.add_hky_pi_G(run.hky_pi()[Real_seq_letter::G]);
  params_builder.add_hky_pi_T(run.hky_pi()[Real_seq_letter::T]);
  params_builder.add_pop_t0(run.pop_model().t0());
  params_builder.add_pop_n0(run.pop_model().pop_at_t0());
  params_builder.add_pop_g(run.pop_model().growth_rate());

  params_builder.add_only_displacing_inner_nodes(run.only_displacing_inner_nodes());
  params_builder.add_topology_moves_enabled(run.topology_moves_enabled());
  params_builder.add_repartitioning_enabled(run.repartitioning_enabled());
  params_builder.add_alpha_move_enabled(run.alpha_move_enabled());
  params_builder.add_mu_move_enabled(run.mu_move_enabled());

  params_builder.add_log_posterior(run.log_G() + run.log_coalescent_prior() + run.log_other_priors());
  params_builder.add_log_other_priors(run.log_other_priors());
  params_builder.add_log_coalescent_prior(run.log_coalescent_prior());
  params_builder.add_log_G(run.log_G());
  params_builder.add_total_branch_length(calc_T(run.tree()));

  auto api_params = params_builder.Finish();

  fbb.FinishSizePrefixed(api_params);
  return fbb.Release();
}

auto apply_api_params_to_run(const uint8_t* params_fb, Run& run) -> void {
  // Not every param is settable, but we'll do our best
  // FIXME: Make more things settable, e.g. coalescent prior staircase resolution
  auto api_params = flatbuffers::GetSizePrefixedRoot<api::Params>(params_fb);

  run.set_step(api_params->step());
  run.set_local_moves_per_global_move(api_params->num_local_moves_per_global_move());
  run.set_num_parts(api_params->num_parts());
  run.set_mu(api_params->mu());
  run.set_alpha(api_params->alpha());
  if (api_params->nu() != nullptr) {
    std::vector<double> new_nu{api_params->nu()->cbegin(), api_params->nu()->cend()};
    run.set_nu(new_nu);
  } else {
    run.set_nu(std::vector<double>(run.tree().num_sites(), 1.0));
  }
  run.set_hky_kappa(api_params->hky_kappa());
  run.set_hky_pi({api_params->hky_pi_A(),
                  api_params->hky_pi_C(),
                  api_params->hky_pi_G(),
                  api_params->hky_pi_T()});
  run.set_pop_model(Exp_pop_model{api_params->pop_t0(), api_params->pop_n0(), api_params->pop_g()});

  run.set_only_displacing_inner_nodes(api_params->only_displacing_inner_nodes());
  run.set_topology_moves_enabled(api_params->topology_moves_enabled());
  run.set_repartitioning_enabled(api_params->repartitioning_enabled());
  run.set_alpha_move_enabled(api_params->alpha_move_enabled());
  run.set_mu_move_enabled(api_params->mu_move_enabled());

  // Derived quantities not settable
}

auto api_probe_site_states_on_tree(
    const Phylo_tree& tree,
    double pop_t0,
    double pop_n0,
    double pop_g,
    int32_t site,
    double t_start,
    double t_end,
    int32_t num_t_cells,
    double* out_values)
    -> void {
  
  auto pop_model = Exp_pop_model{pop_t0, pop_n0, pop_g};

  auto results = probe_site_states_on_tree(tree, pop_model, site, t_start, t_end, num_t_cells);

  for (auto s : k_all_real_seq_letters) {
    auto si = index_of(s);
    for (auto cell = 0; cell != num_t_cells; ++cell) {
      out_values[si * num_t_cells + cell] = results[si].at_cell(cell);
    }
  }
}

auto api_probe_ancestors_on_tree(
    const Phylo_tree& tree,
    double pop_t0,
    double pop_n0,
    double pop_g,
    const int32_t* marked_ancestors,
    int32_t num_marked_ancestors,
    double t_start,
    double t_end,
    int32_t num_t_cells,
    double* out_values)
    -> void {
  
  auto pop_model = Exp_pop_model{pop_t0, pop_n0, pop_g};

  if (num_marked_ancestors < 0) {
    throw std::invalid_argument(absl::StrFormat(
        "Number of marked ancestor cannot be negative (got %d)", num_marked_ancestors));
  }

  auto marked_ancestors_span = std::span{marked_ancestors, static_cast<size_t>(num_marked_ancestors)};
  auto results = probe_ancestors_on_tree(tree, pop_model, marked_ancestors_span, t_start, t_end, num_t_cells);

  auto k = num_marked_ancestors;
  for (auto i = 0; i != (k+1); ++i) {
    for (auto cell = 0; cell != num_t_cells; ++cell) {
      out_values[i * num_t_cells + cell] = results[i].at_cell(cell);
    }
  }
}

auto api_render_population_curve(
    double pop_t0,
    double pop_n0,
    double pop_g,
    double t_start,
    double t_end,
    int32_t num_t_cells,
    double* out_values)
    -> void {
  
  auto pop_model = Exp_pop_model{pop_t0, pop_n0, pop_g};

  auto results = render_population_curve(pop_model, t_start, t_end, num_t_cells);

  for (auto cell = 0; cell != num_t_cells; ++cell) {
    out_values[cell] = results.at_cell(cell);
  }
}

auto api_export_mcc_tree(const Mcc_tree& tree, bool inner_nodes_defined_as_mrcas_of_tips) -> std::string {
  auto t_for = [&](const auto& node) {
    return inner_nodes_defined_as_mrcas_of_tips ? tree.at(node).t_mrca() : tree.at(node).t();
  };
  
  auto ss = std::ostringstream{};
  for (const auto& [node, children_so_far] : traversal(tree)) {
    if (children_so_far == 0) {
      // Entering `node`
      if (tree.at(node).is_inner_node()) {
        ss << "(";
      }
    }

    if (0 < children_so_far && children_so_far < std::ssize(tree.at(node).children)) {
      // In between two children
      ss << ",";
    }

    if (children_so_far == std::ssize(tree.at(node).children)) {
      // Exiting `node`
      if (tree.at(node).is_inner_node()) {
        ss << ")";
      }
      if (tree.at(node).is_tip()) {
        const auto& base_tree = *tree.base_trees().front();
        ss << base_tree.at(node).name;
      } else {
        ss << absl::StreamFormat("NODE_%07d", node);
      }
      auto branch_length_in_years = (node == tree.root ? 0.0 : t_for(node) - t_for(tree.at(node).parent)) / 365.0;
      ss << ":" << absl::StreamFormat("%f", branch_length_in_years);
    }
  }
  ss << ";";
  return std::move(ss).str();  // Idiom to move-construct result instead of copying
}

}  // namespace delphy
