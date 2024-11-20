#ifndef DELPHY_API_H_
#define DELPHY_API_H_

#include <cstdint>
#ifdef __APPLE__
#include <machine/endian.h>
#else
#include <endian.h>
#endif
#include <cstddef>

#include "api_generated.h"
#include "run.h"
#include "mcc_tree.h"

namespace delphy {

auto phylo_tree_to_api_tree(const Phylo_tree& tree) -> flatbuffers::DetachedBuffer;
auto phylo_tree_to_api_tree_info(const Phylo_tree& tree) -> flatbuffers::DetachedBuffer;
auto api_tree_and_tree_info_to_phylo_tree(const uint8_t* tree_fb, const uint8_t* tree_info_fb) -> Phylo_tree;

auto run_to_api_params(const Run& run) -> flatbuffers::DetachedBuffer;
auto apply_api_params_to_run(const uint8_t* params_fb, Run& run) -> void;

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
    -> void;

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
    -> void;

auto api_render_population_curve(
    double pop_t0,
    double pop_n0,
    double pop_g,
    double t_start,
    double t_end,
    int32_t num_t_cells,
    double* out_values)
    -> void;

auto api_export_mcc_tree(const Mcc_tree& tree, bool inner_nodes_defined_as_mrcas_of_tips) -> std::string;

}  // namespace delphy

#endif // DELPHY_API_H_
