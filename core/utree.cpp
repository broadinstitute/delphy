#include "utree.h"

#include <cmath>
#include <queue>
#include <random>

namespace delphy {

auto Utree::reset_focus(Node_index F) -> void {
  focus = F;
  nodes[F].arc_to_focus = k_no_arc;
  for (auto [arc, direction] : annotated_arc_euler_tour(*this, F)) {
    if (direction == Arc_direction::entering) {
      nodes[target(arc)].arc_to_focus = mate(arc);
    }
  }
}

// Utree_builder: incremental tree construction engine.
//
// Builds a Utree by inserting tips one at a time with greedy parsimony placement.
// Each tip is attached approximately at the edge where it introduces the fewest new site
// deltas, found via a priority-queue-driven branch-and-bound search.
//
// The builder maintains a "focus node" whose sequence relative to the reference is tracked
// in deltas_ref_to_focus.  When evaluating a candidate edge for attaching a new tip X, the
// cost equals the number of site deltas on the new M-X edge, computed cheaply from
// focus_to_X_deltas_ and the arc deltas on the candidate edge.
//
// A "focal arc" is an arc whose origin is the current focus node.  Several methods below
// operate on focal arcs, which allows cheap cost evaluation without moving the focus.
//
// This class is an implementation detail — not exposed in the header.  build_guide_tree
// (and future Round 3 methods) create a builder and call add_tip in their chosen order.
class Utree_builder {
 public:
  Utree_builder(Real_sequence ref_sequence,
                const std::vector<Tip_desc>& tip_descs,
                absl::BitGenRef bitgen)
      : tip_descs_{tip_descs}, bitgen_{bitgen} {
    auto N = std::ssize(tip_descs);
    if (N > 0) {
      tree_ = Utree::make_empty(N);
    }
    tree_.ref_sequence = std::move(ref_sequence);
    L_ = static_cast<int>(std::ssize(tree_.ref_sequence));
    sqrt_6L_ = std::sqrt(6.0 * L_);
  }

  // Add tip X to the tree.  Tips can be added in any order, but each tip index
  // in [0, N) must be added exactly once.
  auto add_tip(int X) -> void {
    CHECK_GE(X, 0);
    CHECK_LT(X, std::ssize(tip_descs_));
    CHECK_LT(tips_added_, std::ssize(tip_descs_));
    if (tips_added_ == 1) {
      CHECK_NE(X, tree_.focus) << "Tip " << X << " already added";
    } else if (tips_added_ >= 2) {
      CHECK_EQ(tree_.degree(X), 0) << "Tip " << X << " already added";
    }

    if (tips_added_ == 0) {
      add_first_tip(X);
    } else {
      update_globally_missing_sites(X);
      init_focus_to_X_deltas(X);
      auto best_arc = find_best_attachment_arc(X);
      if (best_arc == k_no_arc) {
        attach_tip_directly_to_isolated_focus(X);
      } else {
        move_focus_updating_focus_to_X_deltas(tree_.origin(best_arc), X);
        attach_tip_to_focal_arc(X, best_arc);
      }
    }
    ++tips_added_;
  }

  // Return the finished tree, transferring ownership.
  auto finish() -> Utree { return std::move(tree_); }

 private:
  static constexpr auto k_min_pruning_threshold = 2;
  static constexpr auto pq_cmp = std::greater<>{};
  using Pq_entry = std::pair<int, Arc_index>;

  // Under JC69, the expected number of same-site blips (mutation+reversal pairs) on
  // a path with `cost` differing sites out of L total is approximately Poisson(sigma^2),
  // where sigma = cost / sqrt(6L).  The threshold is the 5-sigma upper bound on how many
  // blips could be simultaneously open, scaled by 10x to absorb site-rate heterogeneity.
  // See plans/2026-05-07-01-adaptive-pruning-threshold.md for the full derivation.
  auto pruning_threshold(int cost) -> int {
    auto sigma = cost / sqrt_6L_;
    auto threshold = static_cast<int>(std::ceil(10.0 * sigma * (sigma + 5)));
    return std::clamp(threshold, k_min_pruning_threshold, L_);
  }

  // Initialize the tree with the very first tip: set it as focus, copy its seq_deltas
  // into deltas_ref_to_focus, and initialize globally_missing_sites from its missation
  // intervals.
  // Pre: tips_added_ == 0.  Post: tree_ has 1 tip, no edges, focus at X.
  auto add_first_tip(int X) -> void {
    CHECK_EQ(tips_added_, 0);
    const auto& tip_X = tip_descs_[X];

    tree_.focus = X;
    tree_.globally_missing_sites = tip_X.missations.intervals;
    for (const auto& sd : tip_X.seq_deltas) {
      CHECK(not tree_.globally_missing_sites.contains(sd.site));
      tree_.deltas_ref_to_focus[sd.site] = {sd.from, sd.to};
    }
  }

  // Update tree_.globally_missing_sites to reflect the intersection with tip X's missing
  // intervals.  If any sites leave globally_missing_sites and tip X has a seq_delta there,
  // add the corresponding entry to deltas_ref_to_focus (the tree's state at that site is
  // now set to tip X's state).
  // Pre: tips_added_ >= 1.
  auto update_globally_missing_sites(int X) -> void {
    CHECK_GE(tips_added_, 1);
    const auto& tip_X = tip_descs_[X];
    const auto& miss_X = tip_X.missations.intervals;

    if (not interval_set_is_subset_of(tree_.globally_missing_sites, miss_X)) {
      auto old_globally_missing_sites = std::move(tree_.globally_missing_sites);
      tree_.globally_missing_sites = intersect_interval_sets(old_globally_missing_sites, miss_X);

      for (const auto& sd : tip_X.seq_deltas) {
        if (old_globally_missing_sites.contains(sd.site)
            && not tree_.globally_missing_sites.contains(sd.site)) {
          // Site sd.site is leaving globally_missing_sites: the tree's state there is
          // implicitly set to tip X's state, which differs from the reference
          push_back_site_deltas(sd, tree_.deltas_ref_to_focus);
        }
      }
    }
  }

  // Compute focus_to_X_deltas_ = (ref -> tip X's sequence) composed with (focus -> ref).
  // Excludes sites missing at tip X.
  // Pre: tips_added_ >= 1.  Post: focus_to_X_deltas_ is valid for tip X.
  auto init_focus_to_X_deltas(int X) -> void {
    CHECK_GE(tips_added_, 1);
    const auto& tip_X = tip_descs_[X];
    const auto& miss_X = tip_X.missations.intervals;

    focus_to_X_deltas_.clear();
    for (const auto& sd : tip_X.seq_deltas) {
      CHECK(not miss_X.contains(sd.site));
      focus_to_X_deltas_[sd.site] = {sd.from, sd.to};
    }
    for (const auto& [site, delta] : tree_.deltas_ref_to_focus) {
      if (not miss_X.contains(site)) {
        push_front_site_deltas({site, delta.to, delta.from}, focus_to_X_deltas_);
      }
    }
  }

  // Branch-and-bound search for the edge where attaching tip X introduces the fewest new
  // site deltas.  Moves the focus during the search.  Returns the best arc, or k_no_arc if
  // the tree has no edges (tips_added_ == 1).
  // Pre: focus_to_X_deltas_ has been initialized for tip X.
  auto find_best_attachment_arc(int X) -> Arc_index {
    // The search explores edges in order of increasing attachment cost (number of new M-X
    // deltas).  Since the priority queue is a min-heap and costs only increase as we move
    // away from the optimum, once the cheapest remaining entry exceeds best_cost + threshold,
    // no remaining entry can improve the best — so we stop.

    // `<=` below ensures a real arc wins ties over k_no_arc
    auto best_cost = static_cast<int>(std::ssize(focus_to_X_deltas_));
    auto best_arc = k_no_arc;

    pq_.clear();

    // Evaluate all arcs outgoing from focus
    for (auto a : tree_.nodes[tree_.focus].arcs) {
      if (a != k_no_arc) {
        auto cost = eval_focal_arc(a);
        if (cost <= best_cost) {
          best_cost = cost;
          best_arc = a;
        }
        pq_.push_back({cost, a});
      }
    }
    std::make_heap(pq_.begin(), pq_.end(), pq_cmp);

    // Process priority queue
    while (not pq_.empty()) {
      std::pop_heap(pq_.begin(), pq_.end(), pq_cmp);
      auto [priority, arc_R] = pq_.back();
      pq_.pop_back();

      // All remaining entries have priority >= this one; if even this one is too expensive,
      // no further entry can improve the best
      if (priority > best_cost + pruning_threshold(best_cost)) {
        break;
      }

      move_focus_updating_focus_to_X_deltas(tree_.target(arc_R), X);

      // Evaluate arcs outgoing from the new focus (excluding the one we came from)
      for (auto a : tree_.nodes[tree_.focus].arcs) {
        if (a != k_no_arc && a != tree_.mate(arc_R)) {
          auto cost = eval_focal_arc(a);
          if (cost <= best_cost) {
            best_cost = cost;
            best_arc = a;
          }
          pq_.push_back({cost, a});
          std::push_heap(pq_.begin(), pq_.end(), pq_cmp);
        }
      }
    }

    return best_arc;
  }

  // Attach tip X with a direct edge from the current focus (no edge splitting).
  // Used when the focus is an isolated node (degree 0), i.e., the second tip being added.
  // Pre: focus_to_X_deltas_ is valid, tips_added_ == 1, focus has degree 0.
  // Post: tip X is wired into the tree.
  auto attach_tip_directly_to_isolated_focus(int X) -> void {
    CHECK_EQ(tips_added_, 1);
    CHECK_EQ(tree_.degree(tree_.focus), 0);
    auto arc_focus_X = tree_.add_arc(tree_.focus, X);
    for (const auto& [site, delta] : focus_to_X_deltas_) {
      tree_.arcs[arc_focus_X].deltas[site] = delta;
      tree_.arcs[tree_.mate(arc_focus_X)].deltas[site] = {delta.to, delta.from};
    }
    tree_.nodes[X].arc_to_focus = tree_.mate(arc_focus_X);
  }

  // Attach tip X by splitting best_arc's edge, inserting a new inner node M, and connecting
  // tip X to M.  Distributes the old edge's deltas between the two new edges to minimize
  // mutations on the M-X edge.
  // Pre: focus_to_X_deltas_ is valid, best_arc is a focal arc (origin == focus).
  // Post: tip X is wired into the tree, inner node count incremented.
  auto attach_tip_to_focal_arc(int X, Arc_index best_arc) -> void {
    CHECK_NE(best_arc, k_no_arc);
    CHECK_EQ(tree_.origin(best_arc), tree_.focus);
    const auto& miss_X = tip_descs_[X].missations.intervals;

    auto M = tree_.num_tips + tree_.num_inner_nodes_so_far;
    tree_.num_inner_nodes_so_far += 1;

    // Build M-to-X deltas incrementally during split_edge: start from focus_to_X_deltas
    // and adjust for each delta placed on the A-M side
    M_to_X_deltas_ = focus_to_X_deltas_;

    tree_.split_edge(best_arc, M, [&](Seq_delta sd, Node_index A, Node_index B) -> Node_index {
      // Distribute each delta to minimize new mutations on the M-X edge.
      // A side (M gets B's state) vs B side (M gets A's state).
      if (miss_X.contains(sd.site) || tree_.globally_missing_sites.contains(sd.site)) {
        auto side = std::bernoulli_distribution{0.5}(bitgen_) ? A : B;
        if (side == A) {
          pop_front_site_deltas({sd.site, sd.from, sd.to}, M_to_X_deltas_);
        }
        return side;
      }
      auto it = focus_to_X_deltas_.find(sd.site);
      if (it != focus_to_X_deltas_.end() && it->second.to == sd.to) {
        // X matches B's state: put delta on A side so M gets B's (= X's) state
        pop_front_site_deltas({sd.site, sd.from, sd.to}, M_to_X_deltas_);
        return A;
      } else if (it == focus_to_X_deltas_.end()) {
        // X matches A's (focus's) state: put delta on B side so M gets A's (= X's) state
        return B;
      } else {
        // X matches neither A nor B: random side
        auto side = std::bernoulli_distribution{0.5}(bitgen_) ? A : B;
        if (side == A) {
          pop_front_site_deltas({sd.site, sd.from, sd.to}, M_to_X_deltas_);
        }
        return side;
      }
    });

    CHECK_EQ(tree_.target(tree_.nodes[M].arc_to_focus), tree_.focus);

    // Wire M-X edge
    auto arc_MX = tree_.add_arc(M, X);
    for (const auto& [site, delta] : M_to_X_deltas_) {
      tree_.arcs[arc_MX].deltas[site] = delta;
      tree_.arcs[tree_.mate(arc_MX)].deltas[site] = {delta.to, delta.from};
    }
    tree_.nodes[X].arc_to_focus = tree_.mate(arc_MX);
  }

  // Move the focus to `target`, updating focus_to_X_deltas_ in tandem by applying
  // pop_front_site_deltas for each arc delta at a non-missing site.
  auto move_focus_updating_focus_to_X_deltas(Node_index target, Node_index X) -> void {
    const auto& miss_X = tip_descs_[X].missations.intervals;
    tree_.move_focus_to(target, [&](Arc_index a) {
      for (const auto& [site, delta] : tree_.arcs[a].deltas) {
        if (not miss_X.contains(site) && not tree_.globally_missing_sites.contains(site)) {
          pop_front_site_deltas({site, delta.from, delta.to}, focus_to_X_deltas_);
        }
      }
    });
  }

  // Evaluate the cost of attaching the current tip at the edge corresponding to a focal arc
  // (an arc whose origin is the current focus).  Cost = number of site deltas that would
  // appear on the new M-X edge.
  // Pre: origin(a) == tree_.focus, focus_to_X_deltas_ is valid.
  auto eval_focal_arc(Arc_index a) -> int {
    CHECK_EQ(tree_.origin(a), tree_.focus);
    auto savings = 0;
    for (const auto& [site, delta] : tree_.arcs[a].deltas) {
      auto it = focus_to_X_deltas_.find(site);
      if (it != focus_to_X_deltas_.end() && it->second.to == delta.to) {
        ++savings;
      }
    }
    return static_cast<int>(std::ssize(focus_to_X_deltas_)) - savings;
  }

  const std::vector<Tip_desc>& tip_descs_;
  absl::BitGenRef bitgen_;
  Utree tree_;
  int tips_added_ = 0;
  int L_ = 0;
  double sqrt_6L_ = 0.0;

  Heap_site_deltas focus_to_X_deltas_;   // Deltas from focus to tip X being added
  Heap_site_deltas M_to_X_deltas_;       // Deltas from the new inner node M to tip X
  std::vector<Pq_entry> pq_;             // Min-heap for branch-and-bound search
};

// Build a rough "guide tree" by inserting tips one at a time in input order.
// Each tip is attached at the edge where it introduces the fewest new site deltas,
// found via a priority-queue-driven branch-and-bound search.
auto build_guide_tree(
    Real_sequence ref_sequence,
    const std::vector<Tip_desc>& tip_descs,
    absl::BitGenRef bitgen)
    -> Utree {
  auto builder = Utree_builder{std::move(ref_sequence), tip_descs, bitgen};
  for (auto k = 0; k < std::ssize(tip_descs); ++k) {
    builder.add_tip(k);
  }
  return builder.finish();
}

auto assert_utree_integrity(const Utree& tree, bool force) -> void {
  if (not estd::is_debug_enabled && not force) { return; }

  auto num_nodes = tree.num_tips + tree.num_inner_nodes_so_far;

  // Empty tree: nothing to check
  if (num_nodes == 0) {
    CHECK_EQ(tree.focus, k_no_node);
    return;
  }

  // Node degree: a single-node tree has degree 0; otherwise tips have 1, inner nodes have 3
  if (num_nodes == 1) {
    CHECK_EQ(tree.degree(0), 0);
  } else {
    for (auto i = Node_index{0}; i < tree.num_tips; ++i) {
      CHECK_EQ(tree.degree(i), 1) << "Tip " << i;
    }
    for (auto i = tree.num_tips; i < num_nodes; ++i) {
      CHECK_EQ(tree.degree(i), 3) << "Inner node " << i;
    }
  }

  // Arc free list: walk and collect free pair base indices
  auto free_pairs = absl::flat_hash_set<Arc_index>{};
  auto num_arcs = std::ssize(tree.arcs);
  for (auto cur = tree.arc_free_list_head; cur != k_no_arc; cur = tree.arcs[cur].next_free_arc) {
    CHECK_GE(cur, 0) << "Free list index out of range";
    CHECK_LT(cur, num_arcs) << "Free list index out of range";
    auto base = cur & ~1;
    CHECK_EQ(cur, base) << "Free list entry " << cur << " is not the base of its pair";
    CHECK(free_pairs.insert(base).second) << "Cycle in free list at arc " << base;
  }

  // Arc wiring (allocated arcs only)
  for (auto base = Arc_index{0}; base < num_arcs; base += 2) {
    if (free_pairs.contains(base)) { continue; }
    auto a = base;
    auto a_mate = tree.mate(a);

    auto target_a = tree.target(a);
    auto origin_a = tree.origin(a);
    CHECK_GE(target_a, 0) << "Arc " << a << " target out of range";
    CHECK_LT(target_a, num_nodes) << "Arc " << a << " target " << target_a << " out of range";
    CHECK_GE(origin_a, 0) << "Arc " << a << " origin out of range";
    CHECK_LT(origin_a, num_nodes) << "Arc " << a << " origin " << origin_a << " out of range";
    CHECK_NE(origin_a, target_a) << "Self-loop at arc " << a;

    // Target contains mate in its arc slots
    auto found_mate_in_target = false;
    for (auto slot : tree.nodes[target_a].arcs) {
      if (slot == a_mate) { found_mate_in_target = true; break; }
    }
    CHECK(found_mate_in_target) << "Arc " << a << " target " << target_a
        << " does not contain mate " << a_mate;

    // Origin contains a in its arc slots
    auto found_a_in_origin = false;
    for (auto slot : tree.nodes[origin_a].arcs) {
      if (slot == a) { found_a_in_origin = true; break; }
    }
    CHECK(found_a_in_origin) << "Arc " << a << " origin " << origin_a
        << " does not contain arc " << a;

    // Deltas are exact reverses
    CHECK_EQ(std::ssize(tree.arcs[a].deltas), std::ssize(tree.arcs[a_mate].deltas))
        << "Arc " << a << " and mate have different delta counts";
    for (const auto& [site, delta] : tree.arcs[a].deltas) {
      CHECK(not tree.globally_missing_sites.contains(site))
          << "Arc " << a << " has delta at globally missing site " << site;
      auto it = tree.arcs[a_mate].deltas.find(site);
      CHECK(it != tree.arcs[a_mate].deltas.end())
          << "Arc " << a << " has delta at site " << site << " but mate does not";
      CHECK_EQ(it->second.from, delta.to)
          << "Arc " << a << " site " << site << ": mate from != arc to";
      CHECK_EQ(it->second.to, delta.from)
          << "Arc " << a << " site " << site << ": mate to != arc from";
    }
  }

  // Graph is a tree + arc_to_focus via Euler tour from focus.
  // Each non-focus node should be entered exactly once. A cycle or extra edge would cause
  // a node to be entered more than once; a disconnected node would never be entered.
  CHECK_GE(tree.focus, 0);
  CHECK_LT(tree.focus, num_nodes);
  CHECK_EQ(tree.nodes[tree.focus].arc_to_focus, k_no_arc);

  auto enter_count = std::vector<int>(num_nodes, 0);
  enter_count[tree.focus] = 1;  // focus is the DFS source, not yielded
  for (auto [arc, direction] : annotated_arc_euler_tour(tree, tree.focus)) {
    if (direction == Arc_direction::entering) {
      auto node = tree.target(arc);
      CHECK_GE(node, 0) << "Euler tour target out of range";
      CHECK_LT(node, num_nodes) << "Euler tour target " << node << " out of range";
      ++enter_count[node];
      CHECK_LE(enter_count[node], 1)
          << "Node " << node << " entered more than once (not a tree)";
    } else {
      // Leaving arc goes from its origin toward the focus — it is that node's arc_to_focus
      auto node = tree.origin(arc);
      CHECK_GE(node, 0) << "Euler tour origin out of range";
      CHECK_LT(node, num_nodes) << "Euler tour origin " << node << " out of range";
      CHECK_EQ(tree.nodes[node].arc_to_focus, arc)
          << "Node " << node << " arc_to_focus mismatch";
    }
  }
  for (auto i = Node_index{0}; i < num_nodes; ++i) {
    CHECK_EQ(enter_count[i], 1) << "Node " << i << " not reached by Euler tour";
  }

  // deltas_ref_to_focus sanity
  auto L = std::ssize(tree.ref_sequence);
  for (const auto& [site, delta] : tree.deltas_ref_to_focus) {
    CHECK_GE(site, 0) << "deltas_ref_to_focus site out of range";
    CHECK_LT(site, L) << "deltas_ref_to_focus site " << site << " out of range";
    CHECK_EQ(delta.from, tree.ref_sequence[site])
        << "deltas_ref_to_focus site " << site << ": from != ref_sequence";
    CHECK_NE(delta.from, delta.to)
        << "deltas_ref_to_focus site " << site << ": self-delta";
    CHECK(not tree.globally_missing_sites.contains(site))
        << "deltas_ref_to_focus at globally missing site " << site;
  }

  // Arc delta sanity (sites in range, no self-deltas)
  for (auto base = Arc_index{0}; base < num_arcs; base += 2) {
    if (free_pairs.contains(base)) { continue; }
    for (const auto& [site, delta] : tree.arcs[base].deltas) {
      CHECK_GE(site, 0) << "Arc " << base << " delta site out of range";
      CHECK_LT(site, L) << "Arc " << base << " delta site " << site << " out of range";
      CHECK_NE(delta.from, delta.to)
          << "Arc " << base << " site " << site << ": self-delta";
    }
  }

  // Path consistency via Euler tour
  auto ref_to_cur = tree.deltas_ref_to_focus;
  for (auto [arc, direction] : annotated_arc_euler_tour(tree, tree.focus)) {
    for (const auto& [site, delta] : tree.arcs[arc].deltas) {
      auto it = ref_to_cur.find(site);
      if (it != ref_to_cur.end()) {
        CHECK_EQ(it->second.to, delta.from)
            << "Path inconsistency at site " << site << " on arc " << arc;
      } else {
        CHECK_EQ(tree.ref_sequence[site], delta.from)
            << "Path inconsistency at site " << site << " on arc " << arc
            << " (expected ref state)";
      }
      push_back_site_deltas({site, delta.from, delta.to}, ref_to_cur);
    }
  }
  // CHECK not CHECK_EQ: Heap_site_deltas has no operator<<
  CHECK(ref_to_cur == tree.deltas_ref_to_focus) << "Euler tour did not round-trip";
}

auto assert_utree_matches_tip_descs(
    const Utree& tree, const std::vector<Tip_desc>& tip_descs, bool force) -> void {
  if (not estd::is_debug_enabled && not force) { return; }

  auto N = std::ssize(tip_descs);
  CHECK_EQ(tree.num_tips, N);

  auto L = std::ssize(tree.ref_sequence);

  // Validate tip_descs inputs
  for (auto i = Node_index{0}; i < N; ++i) {
    const auto& td = tip_descs[i];

    // Missation intervals: non-empty, in range, sorted, non-overlapping, non-consecutive
    td.missations.intervals.assert_valid(L);

    // seq_deltas: unique sites, in range, from == ref, from != to, not missing
    auto seen_sites = absl::flat_hash_set<Site_index>{};
    for (const auto& sd : td.seq_deltas) {
      CHECK_GE(sd.site, 0) << "Tip " << i << " seq_delta site out of range";
      CHECK_LT(sd.site, L) << "Tip " << i << " seq_delta site " << sd.site << " out of range";
      CHECK_EQ(sd.from, tree.ref_sequence[sd.site])
          << "Tip " << i << " seq_delta site " << sd.site << ": from != ref_sequence";
      CHECK_NE(sd.from, sd.to)
          << "Tip " << i << " seq_delta site " << sd.site << ": self-delta";
      CHECK(not td.missations.intervals.contains(sd.site))
          << "Tip " << i << " seq_delta at missing site " << sd.site;
      CHECK(seen_sites.insert(sd.site).second)
          << "Tip " << i << " duplicate seq_delta at site " << sd.site;
    }
  }

  // Globally missing = intersection of all tip missing intervals
  if (N == 0) {
    CHECK(tree.globally_missing_sites.empty());
    return;
  }
  auto expected_globally_missing_sites = tip_descs[0].missations.intervals;
  for (auto i = Node_index{1}; i < N; ++i) {
    expected_globally_missing_sites =
        intersect_interval_sets(expected_globally_missing_sites, tip_descs[i].missations.intervals);
  }
  CHECK_EQ(tree.globally_missing_sites, expected_globally_missing_sites);

  // Tip sequences via Euler tour
  auto ref_to_cur = tree.deltas_ref_to_focus;
  auto check_tip = [&](Node_index tip) {
    const auto& td = tip_descs[tip];
    const auto& miss = td.missations.intervals;

    // Every seq_delta must appear in ref_to_cur
    for (const auto& sd : td.seq_deltas) {
      auto it = ref_to_cur.find(sd.site);
      CHECK(it != ref_to_cur.end())
          << "Tip " << tip << " site " << sd.site << ": expected delta not in ref_to_cur";
      CHECK_EQ(it->second.from, sd.from)
          << "Tip " << tip << " site " << sd.site << ": from mismatch";
      CHECK_EQ(it->second.to, sd.to)
          << "Tip " << tip << " site " << sd.site << ": to mismatch";
    }

    // Every ref_to_cur entry not in seq_deltas must be at a missing site
    auto seq_delta_sites = absl::flat_hash_set<Site_index>{};
    for (const auto& sd : td.seq_deltas) { seq_delta_sites.insert(sd.site); }
    for (const auto& [site, delta] : ref_to_cur) {
      if (not seq_delta_sites.contains(site)) {
        CHECK(miss.contains(site))
            << "Tip " << tip << " site " << site
            << ": extra ref_to_cur entry at non-missing site";
      }
    }
  };

  // Check focus if it's a tip
  if (tree.focus < N) {
    check_tip(tree.focus);
  }

  for (auto [arc, direction] : annotated_arc_euler_tour(tree, tree.focus)) {
    for (const auto& [site, delta] : tree.arcs[arc].deltas) {
      push_back_site_deltas({site, delta.from, delta.to}, ref_to_cur);
    }
    if (direction == Arc_direction::entering) {
      auto node = tree.target(arc);
      if (node < N) {
        check_tip(node);
      }
    }
  }
}

}  // namespace delphy
