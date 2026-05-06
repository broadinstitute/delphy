#ifndef DELPHY_SITE_DELTAS_H_
#define DELPHY_SITE_DELTAS_H_

#include <absl/container/flat_hash_map.h>
#include <absl/log/check.h>

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
using Heap_site_deltas = absl::flat_hash_map<Site_index, Site_delta>;
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
template<typename Site_deltas_map>
auto push_front_site_deltas(
    Seq_delta delta_z,
    Site_deltas_map& site_deltas_post_z_to_x)
    -> void {
  //
  //    delta
  //      |    ___ site_deltas ___                ______ site_deltas ______
  //      V   /                   \              /                         \      |
  // -----X--*----X--X--X---X------*--   =>   --*--X-------X--X--X---X------*--
  //      z  ^                     x            ^  z                        x
  //       post-z                             pre-z
  
  auto [it_post_z_to_x, inserted] = site_deltas_post_z_to_x.try_emplace(delta_z.site, delta_z.from, delta_z.to);
  if (not inserted) {  // Merge delta at z with existing post-z-to-x delta
    auto& [_, delta_post_z_to_x] = *it_post_z_to_x;
    CHECK_EQ(delta_z.to, delta_post_z_to_x.from);
    auto& delta_pre_z_to_x = delta_post_z_to_x;
    delta_pre_z_to_x.from = delta_z.from;
    if (delta_pre_z_to_x.from == delta_pre_z_to_x.to) {
      site_deltas_post_z_to_x.erase(it_post_z_to_x);
    }
  }
}

// Given the site deltas from a tree loc just before z to a tree loc x, and a seq delta at z,
// mutate the site deltas to go from just after z to x.
template<typename Site_deltas_map>
auto pop_front_site_deltas(
    Seq_delta delta_z,
    Site_deltas_map& site_deltas_pre_z_to_x)
    -> void {
  //
  //    delta
  //      |    ___ site_deltas ___                ______ site_deltas ______
  //      V   /                   \              /                         \      |
  // -----X--*----X--X--X---X------*--   <=   --*--X-------X--X--X---X------*--
  //      z  ^                     x            ^  z                        x
  //       post-z                             pre-z
  
  push_front_site_deltas(delta_z.inverse(), site_deltas_pre_z_to_x);
}

// Given the site deltas from a tree loc x to a tree loc just before z, and a seq delta at z,
// mutate the site deltas to go from x to just after z.
template<typename Site_deltas_map>
auto push_back_site_deltas(
    Seq_delta delta_z,
    Site_deltas_map& site_deltas_x_to_pre_z)
    -> void {
  //
  //                           delta                                              |
  //      ___ site_deltas ___    |                ______ site_deltas ______       |
  //     /                   \   V               /                         \      |
  // ---*----X--X--X---X------*--X---   =>   ---*----X--X--X---X---------X--*---  |
  //    x                     ^  z              x                           ^     |
  //                        pre-z                                         post-z  |
  
  auto [it_x_to_pre_z, inserted] = site_deltas_x_to_pre_z.try_emplace(delta_z.site, delta_z.from, delta_z.to);
  if (not inserted) {  // Merge delta at z with existing x-to-pre-z delta
    auto& [_, delta_x_to_pre_z] = *it_x_to_pre_z;
    CHECK_EQ(delta_z.from, delta_x_to_pre_z.to);
    auto& delta_x_to_post_z = delta_x_to_pre_z;
    delta_x_to_post_z.to = delta_z.to;
    if (delta_x_to_post_z.from == delta_x_to_post_z.to) {
      site_deltas_x_to_pre_z.erase(it_x_to_pre_z);
    }
  }
}

// Given the site deltas from a tree loc x to a tree loc just after z, and a seq delta at z,
// mutate the site deltas to go from x to just before z.
template<typename Site_deltas_map>
auto pop_back_site_deltas(
    Seq_delta delta_z,
    Site_deltas_map& site_deltas_x_to_post_z)
    -> void {
  //
  //                                                                      delta   |
  //      ______ site_deltas ______                  ___ site_deltas ___    |     |
  //     /                         \                /                   \   V     |
  // ---*----X--X--X---X---------X--*---   =>   ---*----X--X--X---X------*--X---  |
  //    x                           ^              x                     ^  z     |
  //                              post-z                               pre-z      |
  
  push_back_site_deltas(delta_z.inverse(), site_deltas_x_to_post_z);
}

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
    Phylo_tree_loc x)
    -> void;

// Given the site deltas from x to y and those from y to z, mutate the first site deltas to go from x to z
template<typename Site_deltas_map_dst, typename Site_deltas_map_src>
auto append_site_deltas(Site_deltas_map_dst& x_to_y_deltas, const Site_deltas_map_src& y_to_z_deltas) -> void {
  for (const auto& [l, y_to_z_delta] : y_to_z_deltas) {
    push_back_site_deltas({l, y_to_z_delta.from, y_to_z_delta.to}, x_to_y_deltas);
  }
}

auto calc_site_deltas_between(const Phylo_tree& tree, Phylo_tree_loc x, Phylo_tree_loc y) -> Site_deltas;
auto calc_site_deltas_between(const Phylo_tree& tree, Node_index X, Node_index Y) -> Site_deltas;

}  // namespace delphy

#endif // DELPHY_SITE_DELTAS_H_
