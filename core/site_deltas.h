#ifndef DELPHY_SITE_DELTAS_H_
#define DELPHY_SITE_DELTAS_H_

#include "evo_model.h"
#include "phylo_tree.h"
#include "sequence_overlay.h"

namespace delphy {

struct Site_delta {
  Real_seq_letter from;
  Real_seq_letter to;
  auto operator<=>(const Site_delta& that) const = default;
};
inline auto operator<<(std::ostream& os, const Site_delta& delta) -> std::ostream& {
  return os << absl::StreamFormat("%c->%c", to_char(delta.from), to_char(delta.to));
}
using Site_deltas = Scratch_flat_hash_map<Site_index, Site_delta>;
inline auto operator<<(std::ostream& os, const std::pair<const int, Site_delta>& delta_entry) -> std::ostream& {
  const auto& [l, delta] = delta_entry;
  return os << absl::StreamFormat("%c%d%c", to_char(delta.from), l, to_char(delta.to));
}
inline auto operator<<(std::ostream& os, const Site_deltas& deltas) -> std::ostream& {
  return os << absl::StreamFormat(
      "Site_deltas{%s}",
      absl::StrJoin(deltas, ", ", absl::StreamFormatter()));
}
inline auto site_deltas_entry(
    Site_index l,
    Real_seq_letter from,
    Real_seq_letter to)
    -> std::pair<const Site_index, Site_delta> {
  return std::make_pair(l, Site_delta{from, to});
}

// Given the site deltas from a tree loc just after z to a tree loc x, and a seq delta at z,
// mutate the site deltas to go from just before z to x.
auto push_front_site_deltas(
    Seq_delta delta_z,
    Site_deltas& site_deltas_post_z_to_x)
    -> void;

// Given the site deltas from a tree loc just before z to a tree loc x, and a seq delta at z,
// mutate the site deltas to go from just after z to x.
auto pop_front_site_deltas(
    Seq_delta delta_z,
    Site_deltas& site_deltas_pre_z_to_x)
    -> void;

// Given the site deltas from a tree loc x to a tree loc just before z, and a seq delta at z,
// mutate the site deltas to go from x to just after z.
auto push_back_site_deltas(
    Seq_delta delta_z,
    Site_deltas& site_deltas_x_to_pre_z)
    -> void;

// Given the site deltas from a tree loc x to a tree loc just after z, and a seq delta at z,
// mutate the site deltas to go from x to just before z.
auto pop_back_site_deltas(
    Seq_delta delta_z,
    Site_deltas& site_deltas_x_to_post_z)
    -> void;

// Given the site deltas from a tree loc x to another (unspecified) tree loc z, and a tree loc a at or above x,
// mutate the site deltas to go from a to z.
auto displace_site_deltas_start_upwards(
    const Phylo_tree& tree,
    Site_deltas& site_deltas,
    Phylo_tree_loc x,
    Phylo_tree_loc a)
    -> void;

// Given the site deltas from a tree loc a to another (unspecified) tree loc z, and a tree loc x at or below a,
// mutate the site deltas to go from x to z.
auto displace_site_deltas_start_downwards(
    const Phylo_tree& tree,
    Site_deltas& site_deltas,
    Phylo_tree_loc a,
    Phylo_tree_loc x,
    Scratch_space& scratch)
    -> void;

// Given the site deltas from x to y and those from y to z, mutate the first site deltas to go from x to z
auto append_site_deltas(Site_deltas& x_to_y_deltas, const Site_deltas& y_to_z_deltas) -> void;

auto calc_site_deltas_between(const Phylo_tree& tree, Phylo_tree_loc x, Phylo_tree_loc y, Scratch_space& scratch) -> Site_deltas;
auto calc_site_deltas_between(const Phylo_tree& tree, Node_index X, Node_index Y, Scratch_space& scratch) -> Site_deltas;

}  // namespace delphy

#endif // DELPHY_SITE_DELTAS_H_
