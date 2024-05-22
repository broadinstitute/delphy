#include "site_deltas.h"

#include <absl/log/check.h>

namespace delphy {

auto push_front_site_deltas(
    Seq_delta delta_z,
    Site_deltas& site_deltas_post_z_to_x)
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

auto pop_front_site_deltas(
    Seq_delta delta_z,
    Site_deltas& site_deltas_pre_z_to_x)
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

auto push_back_site_deltas(
    Seq_delta delta_z,
    Site_deltas& site_deltas_x_to_pre_z)
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

auto pop_back_site_deltas(
    Seq_delta delta_z,
    Site_deltas& site_deltas_x_to_post_z)
    -> void {
  //
  //
  //                                                                      delta   |
  //      ______ site_deltas ______                  ___ site_deltas ___    |     |
  //     /                         \                /                   \   V     |
  // ---*----X--X--X---X---------X--*---   =>   ---*----X--X--X---X------*--X---  |
  //    x                           ^              x                     ^  z     |
  //                              post-z                               pre-z      |
  
  push_back_site_deltas(delta_z.inverse(), site_deltas_x_to_post_z);
}

auto displace_site_deltas_start_upwards(
    const Phylo_tree& tree,
    Site_deltas& site_deltas,
    Phylo_tree_loc x,
    Phylo_tree_loc a)
    -> void {
  //
  // site_deltas are the deltas from the sequence at x to the sequence at some other downstream point y.
  // We walk over the mutations z on the path x->a and prepend them to site_deltas
  //
  //
  //                                 +-- ... -- y       +---------- ... --------- y
  //                                 |                  |
  //                             +-<-+-- ...            |            +---*-- ...
  //                  z          V   x                  |            |   x
  //       +-< ... -<-*-<- ... <-+           ===>       |  +-- ... --+
  //       V                                            |  |
  //  --*<-+                                          --+--+
  //    a                                               a
  
  DCHECK(descends_from(tree, x, a));
  CHECK_NE(a.branch, k_no_node);
  auto G = tree.at(a.branch).parent;
  
  for (auto cur = x.branch; cur != G; cur = tree.at(cur).parent) {
    for (const auto& m_z : tree.at(cur).mutations | std::views::reverse) {
      if (a.t <= m_z.t && m_z.t <= x.t) {
        push_front_site_deltas(m_z, site_deltas);
      }
    }
  }
}

auto displace_site_deltas_start_downwards(
    const Phylo_tree& tree,
    Site_deltas& site_deltas,
    Phylo_tree_loc a,
    Phylo_tree_loc x,
    Scratch_space& scratch)
    -> void {
  //
  // site_deltas are the deltas from the sequence at a to the sequence at some other downstream point y.
  // We walk over the mutations z on the path a->x and prepend their inverses to site_deltas
  //
  //
  //    +----------------- ... ---------------- y                +-- ... -- y
  //    |                                                        |
  //    |                        +->-*-- ...                 +---+-- ...
  //    |             z          ^   x                       |   x
  //    |  +-> ... ->-*->- ... >-+                 +-- ... --+
  //    |  ^                                       |
  //  --+>-+                                  --*--+
  //    a                                       a
  //
  // The only awkward thing is that we need to build a temporary list of branches from a to x first,
  // since only the path x->a is easy to walk, but the path a->x is not.
  //

  DCHECK(descends_from(tree, x, a));
  CHECK_NE(a.branch, k_no_node);
  auto G = tree.at(a.branch).parent;
  
  auto branches_x_to_a = Scratch_vector<Branch_index>{scratch};
  for (auto cur = x.branch; cur != G; cur = tree.at(cur).parent) {
    branches_x_to_a.push_back(cur);
  }

  for (const auto& cur : branches_x_to_a | std::views::reverse) {  // branches a->x
    for (const auto& m_z : tree.at(cur).mutations) {
      if (a.t <= m_z.t && m_z.t <= x.t) {
        pop_front_site_deltas(m_z, site_deltas);
      }
    }
  }
}

auto append_site_deltas(Site_deltas& x_to_y_deltas, const Site_deltas& y_to_z_deltas) -> void {
  for (const auto& [l, y_to_z_delta] : y_to_z_deltas) {
    push_back_site_deltas({l, y_to_z_delta.from, y_to_z_delta.to}, x_to_y_deltas);
  }
}

auto calc_site_deltas_between(const Phylo_tree& tree, Phylo_tree_loc x, Phylo_tree_loc y, Scratch_space& scratch) -> Site_deltas {
  auto site_deltas = Site_deltas{scratch};
  auto a = find_MRCA_of(tree, x, y);
  
  // Here: site_deltas goes from y to y
  
  if (y != a) {
    displace_site_deltas_start_upwards(tree, site_deltas, y, a);
  }
  
  // Here: site_deltas goes from a to y
  
  if (x != a) {
    displace_site_deltas_start_downwards(tree, site_deltas, a, x, scratch);
  }

  // Here: site_deltas goes from x to y
  
  return site_deltas;
}

auto calc_site_deltas_between(const Phylo_tree& tree, Node_index X, Node_index Y, Scratch_space& scratch) -> Site_deltas {
  return calc_site_deltas_between(tree, tree.node_loc(X), tree.node_loc(Y), scratch);
}

}  // namespace delphy
