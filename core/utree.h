#ifndef DELPHY_UTREE_H_
#define DELPHY_UTREE_H_

#include <array>
#include <functional>
#include <ranges>
#include <vector>

#include "absl/log/check.h"
#include "absl/random/bit_gen_ref.h"
#include "absl/random/distributions.h"
#include "cppcoro/generator.hpp"

#include "estd.h"
#include "interval_set.h"
#include "phylo_tree.h"
#include "site_deltas.h"

namespace delphy {

// Utree: an unrooted, bifurcating, mutation-annotated tree.
//
// Each undirected edge is stored as two directed arcs (one per direction), linked as mates.
// Each arc stores the site deltas in its direction.  Arcs are allocated in mate pairs at
// consecutive even-odd indices: arc i and arc i^1 are mates (Knuth, TAOCP v4A pp. 21-22).
//
// A distinguished "focus node" gives the tree a temporary root-like orientation.  Each
// non-focus node stores an arc_to_focus: the outgoing arc along the path toward the focus.
// The focus node has arc_to_focus = k_no_arc.  Moving the focus is O(path_length).
//
// `ref_sequence` is the reference nucleotide sequence.  `globally_missing_sites` tracks sites
// missing across ALL tips added so far.  `deltas_ref_to_focus` records the site deltas
// from the reference sequence to the current focus node's sequence: if site l maps to
// {from, to}, then the focus has state `to` at site l while the reference has state `from`.

using Arc_index = int;
inline constexpr auto k_no_arc = Arc_index{-1};

// A directed arc in the Utree.  When allocated, `target` is the destination node.
// When on the free list, `next_free_arc` links to the next free pair.
// `deltas` stores the site deltas from origin to target (cleared when free).
struct Uarc {
  union {
    Node_index target;            // The node this arc points TO (when allocated)
    Arc_index next_free_arc;      // Next free arc pair (when on free list)
  };
  Heap_site_deltas deltas;        // Site deltas FROM origin TO target; cleared when free
};

// A node in the Utree.  Tips have 1 non-k_no_arc arc; inner nodes have 3.
// `arc_to_focus` is the outgoing arc toward the focus node (k_no_arc for the focus itself).
struct Unode {
  std::array<Arc_index, 3> arcs = {k_no_arc, k_no_arc, k_no_arc};
  Arc_index arc_to_focus = k_no_arc;
};

struct Utree {
  Real_sequence ref_sequence;
  Interval_set<> globally_missing_sites;   // Sites missing across ALL tips added so far

  std::vector<Uarc> arcs;            // All arcs (preallocated in mate pairs)
  std::vector<Unode> nodes;          // All nodes (preallocated)
  Arc_index arc_free_list_head = k_no_arc;   // Head of free list of arc pairs
  int num_tips = 0;                  // Tips are at indices [0, num_tips)
  int num_inner_nodes_so_far = 0;    // Inner nodes allocated so far, starting at index num_tips

  Node_index focus = -1;             // Current focus node
  Heap_site_deltas deltas_ref_to_focus;  // Deltas from reference sequence to focus

  // Create an empty Utree pre-allocated for `num_tips` tips.
  // Nodes and arc pairs are sized for a full binary tree, plus 2 scratch pairs for split_edge.
  // All arc pairs start on the free list.
  static auto make_empty(int num_tips) -> Utree {
    auto tree = Utree{};
    tree.num_tips = num_tips;
    auto num_nodes = std::max(1, 2 * num_tips - 1);
    auto num_arc_pairs = std::max(1, 2 * num_tips - 3 + 2);  // +2 scratch for split_edge
    tree.nodes.resize(num_nodes);
    tree.arcs.resize(2 * num_arc_pairs);
    tree.arc_free_list_head = 0;
    for (auto i = 0; i < 2 * num_arc_pairs; i += 2) {
      tree.arcs[i].next_free_arc = (i + 2 < 2 * num_arc_pairs) ? (i + 2) : k_no_arc;
    }
    return tree;
  }

  // Return the mate of an arc (the arc representing the same edge in the opposite direction)
  inline auto mate(Arc_index arc) const -> Arc_index { return arc ^ 1; }

  // Return the origin node of an arc (= the target of its mate)
  inline auto origin(Arc_index arc) const -> Node_index { return arcs[mate(arc)].target; }

  // Return the target node of an arc
  inline auto target(Arc_index arc) const -> Node_index { return arcs[arc].target; }

  // Find the outgoing arc from `node` whose target is `target`, or k_no_arc if none
  auto find_arc(Node_index node, Node_index target) const -> Arc_index {
    for (auto a : nodes[node].arcs) {
      if (a != k_no_arc && arcs[a].target == target) { return a; }
    }
    return k_no_arc;
  }

  // Count of non-k_no_arc arc slots for a node
  auto degree(Node_index node) const -> int {
    auto d = 0;
    for (auto a : nodes[node].arcs) { if (a != k_no_arc) { ++d; } }
    return d;
  }

  // Whether node is a tip (degree-1 node)
  auto is_tip(Node_index node) const -> bool { return degree(node) == 1; }

  // Pick a random tip uniformly from [0, num_tips)
  auto pick_random_tip(absl::BitGenRef bitgen) const -> Node_index {
    return static_cast<Node_index>(absl::Uniform<int>(bitgen, 0, num_tips));
  }

  // Pick a random connected node (tip or inner), skipping degree-0 nodes
  auto pick_random_node(absl::BitGenRef bitgen) const -> Node_index {
    auto num_nodes = num_tips + num_inner_nodes_so_far;
    Node_index node;
    do { node = static_cast<Node_index>(absl::Uniform<int>(bitgen, 0, num_nodes)); }
    while (degree(node) == 0);
    return node;
  }

  // Total site deltas across all edges (each undirected edge counted once)
  auto count_deltas() const -> int {
    auto total = 0;
    for (auto i = Arc_index{0}; i < std::ssize(arcs); i += 2) {
      total += count_arc_deltas(i);
    }
    return total;
  }

  // Number of site deltas on an arc
  auto count_arc_deltas(Arc_index arc) const -> int {
    return static_cast<int>(std::ssize(arcs[arc].deltas));
  }

  // Pop a pair from the free list and return the even (base) index.
  // The two arcs of the pair are at base and base+1.
  auto alloc_arc_pair() -> Arc_index {
    CHECK_NE(arc_free_list_head, k_no_arc);
    auto base = arc_free_list_head;
    arc_free_list_head = arcs[base].next_free_arc;
    arcs[base].target = -1;
    arcs[base + 1].target = -1;
    return base;
  }

  // Push a pair back onto the free list.  Accepts either arc from the pair.
  auto free_arc_pair(Arc_index arc) -> void {
    auto base = arc & ~1;
    arcs[base].deltas.clear();
    arcs[base + 1].deltas.clear();
    arcs[base].next_free_arc = arc_free_list_head;
    arc_free_list_head = base;
  }

  // Add an undirected edge between nodes A and B.  Allocates an arc pair, sets targets,
  // and wires both nodes' arc slots.  Returns the arc A->B; the caller can then set deltas.
  auto add_arc(Node_index A, Node_index B) -> Arc_index {
    auto base = alloc_arc_pair();
    auto arc_AB = base;
    auto arc_BA = base + 1;
    arcs[arc_AB].target = B;
    arcs[arc_BA].target = A;

    auto wire = [&](Node_index node, Arc_index arc) {
      for (auto& a : nodes[node].arcs) {
        if (a == k_no_arc) { a = arc; return; }
      }
      CHECK(false) << "No free arc slot in node " << node;
    };
    wire(A, arc_AB);
    wire(B, arc_BA);

    return arc_AB;
  }

  // Set node F as the focus, recomputing all arc_to_focus fields from scratch via DFS.
  // Unlike move_focus_to, this does NOT update deltas_ref_to_focus.
  auto reset_focus(Node_index F) -> void;

  // Remove the edge between tip X and its neighbor M, leaving M as degree-2.
  // Returns M.  The caller must move the focus away from M (if needed) and call
  // merge_through(M) afterward.
  // Pre: is_tip(X), degree(X) == 1, focus != X, num_tips >= 3.
  auto detach_tip(Node_index X) -> Node_index;

  // Merge the two edges incident to degree-2 inner node M into a single edge, removing M
  // from the tree.  Inverse of split_edge.  Pre: degree(M) == 2, focus != M.
  // Returns the new arc A->B.
  auto merge_through(Node_index M) -> Arc_index;

  // Split edge (A,B) by inserting node M.  Replaces the single edge with two: (A,M) and (M,B).
  // `site_delta_side(Seq_delta delta, Node_index A, Node_index B) -> Node_index` decides which
  // side each site delta goes to: return A for the A-M side, or B for the M-B side.
  // M gets two arc slots filled; the third remains k_no_arc for the caller to attach a tip.
  template<typename Site_delta_side>
  auto split_edge(Arc_index arc_AB, Node_index M, Site_delta_side site_delta_side) -> void;

  // Move the focus from its current position to `target`.
  // Two-pass algorithm: (1) walk target -> old focus, reversing arc_to_focus links;
  // (2) walk old focus -> target, applying arc deltas to deltas_ref_to_focus.
  // Hooks are called at each step: pre_arc_hop(arc) before applying deltas,
  // post_arc_hop(arc) after.
  template<typename Pre_arc_hop, typename Post_arc_hop>
  auto move_focus_to(Node_index target, Pre_arc_hop pre_arc_hop, Post_arc_hop post_arc_hop) -> void;

  // Convenience: move focus with only a pre-hop hook.
  template<typename Pre_arc_hop>
  auto move_focus_to(Node_index target, Pre_arc_hop pre_arc_hop) -> void {
    move_focus_to(target, pre_arc_hop, [](Arc_index){});
  }

  // Convenience: move focus without hooks.
  auto move_focus_to(Node_index target) -> void {
    move_focus_to(target, [](Arc_index){}, [](Arc_index){});
  }
};

// Build a rough "guide tree" by inserting tips one at a time in input order.
// Each tip is attached at the edge where it introduces the fewest new site deltas,
// found via a priority-queue-driven branch-and-bound search.
auto build_guide_tree(
    Real_sequence ref_sequence,
    const std::vector<Tip_desc>& tip_descs,
    absl::BitGenRef bitgen,
    const std::function<void(int,int)>& progress_hook = [](int,int){})
    -> Utree;

// Traverse the guide tree in nearest-first order, calling `callback(tip, closest_prev_tip)`
// for each tip.  `closest_prev_tip` is the nearest already-added tip in the guide tree
// (k_no_node for the first tip).
auto for_each_tip_in_nearest_first_order(
    const Utree& guide_tree,
    absl::BitGenRef bitgen,
    const std::function<void(Node_index tip, Node_index closest_prev_tip)>& callback) -> void;

// Build a refined Utree by adding tips in guide-tree nearest-first order.
auto build_refined_tree(
    const Utree& guide_tree,
    const std::vector<Tip_desc>& tip_descs,
    absl::BitGenRef bitgen,
    const std::function<void(int,int)>& progress_hook = [](int,int){}) -> Utree;

// SPR refinement of tip placements.  For each random tip, detach it, search for a better
// reattachment point starting from a random node, and accept or roll back.
// progress_hook(attempts_so_far, max_attempts, cur_deltas): called after each attempt with
// the current total delta count across all edges.
auto spr_refine_tips(
    Utree& tree, const std::vector<Tip_desc>& tip_descs,
    absl::BitGenRef bitgen,
    const std::function<void(int,int,int)>& progress_hook = [](int,int,int){}) -> void;

enum class Rooting_method { regression, midpoint };

struct Rooting_info {
  Node_index root;
  Rooting_method method;
  double r2;       // R^2 of root-to-tip regression (may be <= 0 if no clock signal)
  double lambda;   // overall mutations per day (not per-site)
  double t_MRCA;   // estimated root date (days since epoch)
};

// Root the Utree at the timed midpoint of a diametral path, then estimate
// the mutation rate and root date via OLS regression of root-to-tip mutation
// counts against tip dates.  Modifies the tree in place: inserts a root node.
// The focus location is undefined after this call.
auto midpoint_root_utree(Utree& tree, const std::vector<Tip_desc>& tip_descs) -> Rooting_info;

// Root the Utree at the position that maximizes R^2 of root-to-tip OLS regression against tip dates.
// Falls back to midpoint rooting if regression is not applicable (e.g., all tips have the same date).
// Modifies the tree in place: inserts a root node.  The focus location is undefined after this call.
auto ols_regression_root_utree(Utree& tree, const std::vector<Tip_desc>& tip_descs,
                               absl::BitGenRef bitgen)
    -> Rooting_info;

// Root the Utree at the position that minimizes chi^2 of root-to-tip GLS regression against tip dates.
// Uses the phylogenetic covariance structure (Sherman-Morrison updates) to properly weight tips.
// Falls back to midpoint rooting if regression is not applicable (e.g., all tips have the same date).
// Modifies the tree in place: inserts a root node.  The focus location is undefined after this call.
auto gls_regression_root_utree(Utree& tree, const std::vector<Tip_desc>& tip_descs,
                               absl::BitGenRef bitgen)
    -> Rooting_info;

// Convert a rooted Utree to a Phylo_tree.
// Moves the focus to `rooting_info.root` if not already there.
auto utree_to_phylo_tree(
    Utree& utree, const Rooting_info& rooting_info, const std::vector<Tip_desc>& tip_descs,
    absl::BitGenRef bitgen) -> Phylo_tree;

// Full pipeline: build guide tree, refine it, root it, estimate rate, convert to Phylo_tree.
auto build_initial_phylo_tree(
    Real_sequence ref_sequence, std::vector<Tip_desc> tip_descs,
    absl::BitGenRef bitgen,
    const std::function<void(int,int)>& guide_tree_progress_hook = [](int,int){},
    const std::function<void(int,int,int)>& refined_tree_progress_hook = [](int,int,int){},
    const std::function<void(int,int,int)>& spr_refine_progress_hook = [](int,int,int){},
    const std::function<void(const Rooting_info&)>& rooting_hook = [](const Rooting_info&){})
    -> Phylo_tree;

enum class Arc_direction { entering, leaving };
struct Annotated_arc {
  Arc_index arc;
  Arc_direction direction;
  auto operator<=>(const Annotated_arc&) const = default;
};

// DFS Euler tour of a Utree starting at `source`, yielding each arc hopped along the way
// annotated with direction: `entering` (away from source) or `leaving` (backtracking).
// Every edge is visited exactly twice (once in each direction).
inline auto annotated_arc_euler_tour(const Utree& tree, Node_index source)
    -> cppcoro::generator<Annotated_arc> {
  auto stack = std::vector<Annotated_arc>{};
  for (auto a : tree.nodes[source].arcs) {
    if (a != k_no_arc) {
      stack.push_back({tree.mate(a), Arc_direction::leaving});
      stack.push_back({a, Arc_direction::entering});
    }
  }
  while (not stack.empty()) {
    auto [arc, direction] = stack.back();
    stack.pop_back();
    co_yield Annotated_arc{arc, direction};
    if (direction == Arc_direction::entering) {
      auto B = tree.target(arc);
      for (auto a : tree.nodes[B].arcs) {
        if (a != k_no_arc && a != tree.mate(arc)) {
          stack.push_back({tree.mate(a), Arc_direction::leaving});
          stack.push_back({a, Arc_direction::entering});
        }
      }
    }
  }
}

// Like annotated_arc_euler_tour, but yields only the Arc_index (discarding direction).
inline auto arc_euler_tour(const Utree& tree, Node_index source)
    -> cppcoro::generator<Arc_index> {
  for (auto [arc, direction] : annotated_arc_euler_tour(tree, source)) {
    co_yield arc;
  }
}

// Iterate over all tip node indices in the tree: [0, num_tips).
inline auto utree_tips(const Utree& tree) {
  return std::views::iota(Node_index{0}, static_cast<Node_index>(tree.num_tips));
}

// CHECK-fails on any structural inconsistency in the tree (node degrees, arc wiring,
// arc_to_focus pointers, free list, path consistency of deltas).
auto assert_utree_integrity(const Utree& tree, bool force = false) -> void;

// CHECK-fails if the tree's tip sequences or globally_missing_sites don't match the input tip_descs.
// Assumes the tree already passes assert_utree_integrity.
auto assert_utree_matches_tip_descs(
    const Utree& tree, const std::vector<Tip_desc>& tip_descs, bool force = false) -> void;

template<typename Pre_arc_hop, typename Post_arc_hop>
auto Utree::move_focus_to(Node_index target, Pre_arc_hop pre_arc_hop, Post_arc_hop post_arc_hop) -> void {
  if (target == focus) {
    return;
  }

  // Pass 1 (Reverse): walk from target toward current focus, flipping arc_to_focus links
  // so they point from old focus toward new focus (target).
  auto cur = target;
  auto prev_saved = nodes[cur].arc_to_focus;
  nodes[cur].arc_to_focus = k_no_arc;

  while (prev_saved != k_no_arc) {
    auto next_node = arcs[prev_saved].target;
    auto next_saved = nodes[next_node].arc_to_focus;
    nodes[next_node].arc_to_focus = mate(prev_saved);
    prev_saved = next_saved;
  }

  // Pass 2 (Forward): follow the now-reversed arc_to_focus links from old focus toward
  // target, applying each arc's deltas to deltas_ref_to_focus and calling hooks.
  cur = focus;
  while (cur != target) {
    auto arc_R = nodes[cur].arc_to_focus;
    pre_arc_hop(arc_R);
    for (const auto& [site, delta] : arcs[arc_R].deltas) {
      push_back_site_deltas({site, delta.from, delta.to}, deltas_ref_to_focus);
    }
    post_arc_hop(arc_R);
    cur = arcs[arc_R].target;
  }

  focus = target;
}

template<typename Site_delta_side>
auto Utree::split_edge(Arc_index arc_AB, Node_index M, Site_delta_side site_delta_side) -> void {
  auto arc_BA = mate(arc_AB);
  auto A = origin(arc_AB);
  auto B = target(arc_AB);

  // Allocate new edges: (A,M) and (M,B)
  auto arc_AM = alloc_arc_pair();
  auto arc_MA = mate(arc_AM);
  auto arc_MB = alloc_arc_pair();
  auto arc_BM = mate(arc_MB);

  arcs[arc_AM].target = M;
  arcs[arc_MA].target = A;
  arcs[arc_MB].target = B;
  arcs[arc_BM].target = M;

  // Distribute deltas from old edge to the two new edges
  for (const auto& [site, delta] : arcs[arc_AB].deltas) {
    auto side = site_delta_side({site, delta.from, delta.to}, A, B);
    CHECK(side == A || side == B) << "site_delta_side returned " << side
        << ", expected A=" << A << " or B=" << B;
    if (side == A) {
      arcs[arc_AM].deltas[site] = {delta.from, delta.to};
      arcs[arc_MA].deltas[site] = {delta.to, delta.from};
    } else {
      arcs[arc_MB].deltas[site] = {delta.from, delta.to};
      arcs[arc_BM].deltas[site] = {delta.to, delta.from};
    }
  }

  // Wire M's arc slots (third slot left empty for caller to attach a tip)
  nodes[M].arcs[0] = arc_MA;
  nodes[M].arcs[1] = arc_MB;
  nodes[M].arcs[2] = k_no_arc;

  // Update A's arc slot: replace old arc_AB with arc_AM
  for (auto& a : nodes[A].arcs) {
    if (a == arc_AB) { a = arc_AM; break; }
  }

  // Update B's arc slot: replace old arc_BA with arc_BM
  for (auto& a : nodes[B].arcs) {
    if (a == arc_BA) { a = arc_BM; break; }
  }

  // Update arc_to_focus for A, B, and M if they pointed at old arcs.
  if (nodes[A].arc_to_focus == arc_AB) {
    nodes[A].arc_to_focus = arc_AM;
    nodes[M].arc_to_focus = arc_MB;
  }
  if (nodes[B].arc_to_focus == arc_BA) {
    nodes[B].arc_to_focus = arc_BM;
    nodes[M].arc_to_focus = arc_MA;
  }

  free_arc_pair(arc_AB);
}

}  // namespace delphy

#endif // DELPHY_UTREE_H_
