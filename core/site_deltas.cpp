#include "site_deltas.h"

#include <absl/log/check.h>

namespace delphy {

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
    Phylo_tree_loc x)
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

  auto branches_x_to_a = Scratch_vector<Branch_index>{};
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

auto calc_site_deltas_between(const Phylo_tree& tree, Phylo_tree_loc x, Phylo_tree_loc y) -> Site_deltas {
  auto site_deltas = Site_deltas{};
  auto a = find_MRCA_of(tree, x, y);

  // Here: site_deltas goes from y to y

  if (y != a) {
    displace_site_deltas_start_upwards(tree, site_deltas, y, a);
  }

  // Here: site_deltas goes from a to y

  if (x != a) {
    displace_site_deltas_start_downwards(tree, site_deltas, a, x);
  }

  // Here: site_deltas goes from x to y

  return site_deltas;
}

auto calc_site_deltas_between(const Phylo_tree& tree, Node_index X, Node_index Y) -> Site_deltas {
  return calc_site_deltas_between(tree, tree.node_loc(X), tree.node_loc(Y));
}

}  // namespace delphy
