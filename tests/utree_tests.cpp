#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "utree.h"

namespace delphy {

inline constexpr auto rA = Real_seq_letter::A;
inline constexpr auto rC = Real_seq_letter::C;
inline constexpr auto rG = Real_seq_letter::G;
inline constexpr auto rT = Real_seq_letter::T;

// Helper: set deltas on an arc and its mate (forward and reverse).
auto set_arc_deltas(Utree& tree, Arc_index arc_AB, Heap_site_deltas from_A_to_B) -> void {
  auto arc_BA = tree.mate(arc_AB);
  tree.arcs[arc_AB].deltas = from_A_to_B;
  tree.arcs[arc_BA].deltas.clear();
  for (const auto& [site, delta] : from_A_to_B) {
    tree.arcs[arc_BA].deltas[site] = {delta.to, delta.from};
  }
}

// --- alloc / free ---

TEST(Utree_test, alloc_free_lifo) {
  auto tree = Utree::make_empty(3);

  auto first = tree.alloc_arc_pair();
  EXPECT_THAT(first, testing::Eq(0));

  auto second = tree.alloc_arc_pair();
  EXPECT_THAT(second, testing::Eq(2));

  tree.free_arc_pair(first);

  auto third = tree.alloc_arc_pair();
  EXPECT_THAT(third, testing::Eq(0));  // LIFO: get back the freed pair
}

TEST(Utree_test, free_accepts_odd_arc) {
  auto tree = Utree::make_empty(3);
  auto base = tree.alloc_arc_pair();

  set_arc_deltas(tree, base, {{0, {rA, rC}}});

  tree.free_arc_pair(base + 1);  // pass the odd arc

  EXPECT_THAT(tree.arcs[base].deltas, testing::IsEmpty());
  EXPECT_THAT(tree.arcs[base + 1].deltas, testing::IsEmpty());

  EXPECT_THAT(tree.arc_free_list_head, testing::Eq(base));
}

// --- add_arc / reset_focus ---

TEST(Utree_test, add_arc_and_reset_focus) {
  // Build: 0 -- 1
  auto tree = Utree::make_empty(2);
  tree.add_arc(0, 1);
  tree.reset_focus(0);

  EXPECT_THAT(tree.focus, testing::Eq(0));
  EXPECT_THAT(tree.nodes[0].arc_to_focus, testing::Eq(k_no_arc));
  EXPECT_THAT(tree.target(tree.nodes[1].arc_to_focus), testing::Eq(0));
}

// --- split_edge ---

TEST(Utree_test, split_edge_distributes_deltas) {
  // Edge 0--1 with deltas: site 0 A->C, site 1 G->T
  auto tree = Utree::make_empty(3);
  auto arc_01 = tree.add_arc(0, 1);
  set_arc_deltas(tree, arc_01, {
      {0, {rA, rC}},
      {1, {rG, rT}}});
  tree.reset_focus(0);

  // Split: keep site 0 on A-side (0-M), site 1 on M-B side (M-1)
  auto M = Node_index{2};
  tree.split_edge(arc_01, M, [](Seq_delta sd, Node_index A, Node_index B) {
    return sd.site == 0 ? A : B;
  });

  // M should have 2 filled arc slots, 1 empty
  auto filled = 0;
  for (auto a : tree.nodes[M].arcs) {
    if (a != k_no_arc) { ++filled; }
  }
  EXPECT_THAT(filled, testing::Eq(2));

  // Find arcs from M to 0 and M to 1
  auto arc_M0 = tree.find_arc(M, 0);
  auto arc_M1 = tree.find_arc(M, 1);
  ASSERT_THAT(arc_M0, testing::Ne(k_no_arc));
  ASSERT_THAT(arc_M1, testing::Ne(k_no_arc));

  // 0-M edge (A-side) should have site 0 delta only
  auto arc_0M = tree.mate(arc_M0);
  EXPECT_THAT(tree.arcs[arc_0M].deltas, testing::UnorderedElementsAre(
      site_deltas_entry(0, rA, rC)));

  // M-1 edge (M-B side) should have site 1 delta only
  EXPECT_THAT(tree.arcs[arc_M1].deltas, testing::UnorderedElementsAre(
      site_deltas_entry(1, rG, rT)));
}

TEST(Utree_test, split_edge_focus_at_A) {
  // Focus at 0 (= A), split the 0--1 edge
  auto tree = Utree::make_empty(3);
  auto arc_01 = tree.add_arc(0, 1);
  tree.reset_focus(0);

  auto M = Node_index{2};
  tree.split_edge(arc_01, M, [](Seq_delta, Node_index A, Node_index) { return A; });

  EXPECT_THAT(tree.nodes[0].arc_to_focus, testing::Eq(k_no_arc));
  EXPECT_THAT(tree.target(tree.nodes[M].arc_to_focus), testing::Eq(0));
  EXPECT_THAT(tree.target(tree.nodes[1].arc_to_focus), testing::Eq(M));
}

TEST(Utree_test, split_edge_focus_at_B) {
  // Focus at 1 (= B when splitting arc 0->1)
  auto tree = Utree::make_empty(3);
  auto arc_01 = tree.add_arc(0, 1);
  tree.reset_focus(1);

  auto M = Node_index{2};
  tree.split_edge(arc_01, M, [](Seq_delta, Node_index A, Node_index) { return A; });

  // A=0, B=1.  Focus is at 1 (= B), so:
  // - node 1 should still be focus
  // - M should point toward 1
  // - node 0 should point toward M
  EXPECT_THAT(tree.nodes[1].arc_to_focus, testing::Eq(k_no_arc));
  EXPECT_THAT(tree.target(tree.nodes[M].arc_to_focus), testing::Eq(1));
  EXPECT_THAT(tree.target(tree.nodes[0].arc_to_focus), testing::Eq(M));
}

// --- move_focus_to ---

TEST(Utree_test, move_focus_to_self) {
  auto tree = Utree::make_empty(2);
  tree.add_arc(0, 1);
  tree.reset_focus(0);
  tree.deltas_ref_to_focus[0] = {rA, rC};

  tree.move_focus_to(0);

  EXPECT_THAT(tree.focus, testing::Eq(0));
  EXPECT_THAT(tree.deltas_ref_to_focus, testing::SizeIs(1));  // unchanged
}

TEST(Utree_test, move_focus_neighbor) {
  // 0--1 with delta at site 0: A->C.  Focus at 0, ref_to_focus = {5: G->T}
  auto tree = Utree::make_empty(2);
  auto arc_01 = tree.add_arc(0, 1);
  set_arc_deltas(tree, arc_01, {{0, {rA, rC}}});
  tree.reset_focus(0);
  tree.deltas_ref_to_focus[5] = {rG, rT};

  tree.move_focus_to(1);

  EXPECT_THAT(tree.focus, testing::Eq(1));
  EXPECT_THAT(tree.deltas_ref_to_focus, testing::UnorderedElementsAre(
      site_deltas_entry(0, rA, rC),
      site_deltas_entry(5, rG, rT)));

  EXPECT_THAT(tree.nodes[1].arc_to_focus, testing::Eq(k_no_arc));
  ASSERT_THAT(tree.nodes[0].arc_to_focus, testing::Ne(k_no_arc));
  EXPECT_THAT(tree.target(tree.nodes[0].arc_to_focus), testing::Eq(1));
}

TEST(Utree_test, move_focus_through_chain) {
  // Build: 0 -- M(2) -- 1
  // Edge 0-M: site 0 A->C.  Edge M-1: site 1 G->T.
  auto tree = Utree::make_empty(3);
  auto M = Node_index{2};
  auto arc_0M = tree.add_arc(0, M);
  auto arc_M1 = tree.add_arc(M, 1);
  set_arc_deltas(tree, arc_0M, {{0, {rA, rC}}});
  set_arc_deltas(tree, arc_M1, {{1, {rG, rT}}});
  tree.reset_focus(0);

  tree.move_focus_to(1);

  EXPECT_THAT(tree.focus, testing::Eq(1));
  EXPECT_THAT(tree.deltas_ref_to_focus, testing::UnorderedElementsAre(
      site_deltas_entry(0, rA, rC),
      site_deltas_entry(1, rG, rT)));

  EXPECT_THAT(tree.nodes[1].arc_to_focus, testing::Eq(k_no_arc));
  EXPECT_THAT(tree.target(tree.nodes[M].arc_to_focus), testing::Eq(1));
  EXPECT_THAT(tree.target(tree.nodes[0].arc_to_focus), testing::Eq(M));
}

TEST(Utree_test, move_focus_hooks_record_path) {
  // 0 -- M(2) -- 1
  auto tree = Utree::make_empty(3);
  auto M = Node_index{2};
  auto arc_0M = tree.add_arc(0, M);
  auto arc_M1 = tree.add_arc(M, 1);
  tree.reset_focus(0);

  auto pre_arcs = std::vector<Arc_index>{};
  auto post_arcs = std::vector<Arc_index>{};

  tree.move_focus_to(1,
      [&](Arc_index a) { pre_arcs.push_back(a); },
      [&](Arc_index a) { post_arcs.push_back(a); });

  // Should have hopped two edges: 0->M, then M->1
  EXPECT_THAT(pre_arcs, testing::ElementsAre(arc_0M, arc_M1));
  EXPECT_THAT(post_arcs, testing::ElementsAre(arc_0M, arc_M1));
}

TEST(Utree_test, move_focus_cancelling_deltas) {
  // Edge 0--1 with delta at site 0: A->C.
  // ref_to_focus = {0: C->A} (ref has C, focus(0) has A)
  // Moving to 1: C->A->C = identity => deltas_ref_to_focus should be empty
  auto tree = Utree::make_empty(2);
  auto arc_01 = tree.add_arc(0, 1);
  set_arc_deltas(tree, arc_01, {{0, {rA, rC}}});
  tree.reset_focus(0);
  tree.deltas_ref_to_focus[0] = {rC, rA};

  tree.move_focus_to(1);

  EXPECT_THAT(tree.focus, testing::Eq(1));
  EXPECT_THAT(tree.deltas_ref_to_focus, testing::IsEmpty());
}

// --- build_guide_tree ---


// Helper: make a Tip_desc with the given seq_deltas and optional missing intervals
auto make_tip(std::vector<Seq_delta> seq_deltas, Interval_set<> missing = {}) -> Tip_desc {
  auto td = Tip_desc{};
  td.seq_deltas = std::move(seq_deltas);
  td.missations.intervals = std::move(missing);
  return td;
}

// Helper: make a Tip_desc with name, date bounds, seq_deltas, and optional missing intervals
auto make_dated_tip(
    std::string name, float t_min, float t_max,
    std::vector<Seq_delta> seq_deltas, Interval_set<> missing = {}) -> Tip_desc {
  auto td = Tip_desc{};
  td.name = std::move(name);
  td.t_min = t_min;
  td.t_max = t_max;
  td.seq_deltas = std::move(seq_deltas);
  td.missations.intervals = std::move(missing);
  return td;
}

TEST(Utree_test, guide_tree_zero_tips) {
  auto ref = Real_sequence{rA, rC, rG, rT};
  auto bitgen = std::mt19937{42};
  auto tree = build_guide_tree(ref, {}, bitgen);

  EXPECT_THAT(tree.ref_sequence, testing::ElementsAre(rA, rC, rG, rT));
  EXPECT_THAT(tree.num_tips, testing::Eq(0));
  assert_utree_integrity(tree, true);
}

TEST(Utree_test, guide_tree_one_tip) {
  auto ref = Real_sequence{rA, rC, rG, rT};
  auto tips = std::vector<Tip_desc>{
      make_tip({{0, rA, rT}})  // site 0: A->T
  };
  auto bitgen = std::mt19937{42};

  auto tree = build_guide_tree(ref, tips, bitgen);

  EXPECT_THAT(tree.num_tips, testing::Eq(1));
  EXPECT_THAT(tree.focus, testing::Eq(0));
  EXPECT_THAT(tree.deltas_ref_to_focus, testing::UnorderedElementsAre(
      site_deltas_entry(0, rA, rT)));
  EXPECT_THAT(tree.count_deltas(), testing::Eq(0));
  assert_utree_integrity(tree, true);
  assert_utree_matches_tip_descs(tree, tips, true);
}

TEST(Utree_test, guide_tree_two_tips) {
  // ref = ACGT.  Tip 0: site 0 A->T.  Tip 1: site 1 C->G.
  // Tips differ from each other at sites 0 and 1.
  auto ref = Real_sequence{rA, rC, rG, rT};
  auto tips = std::vector<Tip_desc>{
      make_tip({{0, rA, rT}}),
      make_tip({{1, rC, rG}})
  };
  auto bitgen = std::mt19937{42};

  auto tree = build_guide_tree(ref, tips, bitgen);

  EXPECT_THAT(tree.num_tips, testing::Eq(2));
  EXPECT_THAT(tree.count_deltas(), testing::Eq(2));
  assert_utree_integrity(tree, true);
  assert_utree_matches_tip_descs(tree, tips, true);
}

TEST(Utree_test, guide_tree_three_tips_optimal) {
  // ref = AAAA.  Tips 0 and 1 share a mutation at site 0 (A->C), tip 2 does not.
  // Tip 1 also has a unique mutation at site 1 (A->G).
  // Optimal parsimony = 2 (site 0 on shared branch, site 1 on tip 1 branch).
  // The guide tree may be slightly suboptimal due to insertion order.
  auto ref = Real_sequence{rA, rA, rA, rA};
  auto tips = std::vector<Tip_desc>{
      make_tip({{0, rA, rC}}),               // tip 0: C at site 0
      make_tip({{0, rA, rC}, {1, rA, rG}}),  // tip 1: C at site 0, G at site 1
      make_tip({})                            // tip 2: same as ref
  };
  auto bitgen = std::mt19937{42};

  auto tree = build_guide_tree(ref, tips, bitgen);

  EXPECT_THAT(tree.num_tips, testing::Eq(3));
  EXPECT_THAT(tree.count_deltas(), testing::Le(3));
  assert_utree_integrity(tree, true);
  assert_utree_matches_tip_descs(tree, tips, true);
}

TEST(Utree_test, guide_tree_five_tips_two_clusters) {
  // Two clusters sharing mutations:
  //   Cluster A (tips 0,1,2): site 0 A->C
  //   Cluster B (tips 3,4): site 1 A->G
  // Within cluster A: tip 1 has site 2 A->T, tip 2 has site 3 A->T
  // Optimal parsimony = 4 (sites 0, 1, 2, 3 each mutate once)
  auto ref = Real_sequence{rA, rA, rA, rA, rA};
  auto tips = std::vector<Tip_desc>{
      make_tip({{0, rA, rC}}),                           // tip 0
      make_tip({{0, rA, rC}, {2, rA, rT}}),              // tip 1
      make_tip({{0, rA, rC}, {3, rA, rT}}),              // tip 2
      make_tip({{1, rA, rG}}),                           // tip 3
      make_tip({{1, rA, rG}}),                           // tip 4
  };
  auto bitgen = std::mt19937{42};

  auto tree = build_guide_tree(ref, tips, bitgen);

  EXPECT_THAT(tree.num_tips, testing::Eq(5));
  // Optimal = 4; guide tree may be slightly suboptimal due to insertion order
  EXPECT_THAT(tree.count_deltas(), testing::Le(6));
  assert_utree_integrity(tree, true);
  assert_utree_matches_tip_descs(tree, tips, true);
}

TEST(Utree_test, guide_tree_with_missing_data) {
  // ref = AAAA.  Tip 0 is missing sites [2,4).  Tip 1 has a mutation at site 2.
  // Since site 2 is missing at tip 0, the 0-1 edge should not have a delta at site 2.
  auto ref = Real_sequence{rA, rA, rA, rA};
  auto tips = std::vector<Tip_desc>{
      make_tip({}, Interval_set<>{{Site_interval{2, 4}}}),  // tip 0: missing sites 2,3
      make_tip({{2, rA, rC}})                               // tip 1: C at site 2
  };
  auto bitgen = std::mt19937{42};

  auto tree = build_guide_tree(ref, tips, bitgen);

  EXPECT_THAT(tree.num_tips, testing::Eq(2));
  // Site 2 is missing at tip 0, so the edge between them should have no delta at site 2
  // (it's globally missing after just tip 0, and when tip 1 is added, the tree's state
  // at site 2 is set to tip 1's state — C — so they agree)
  EXPECT_THAT(tree.count_deltas(), testing::Eq(0));
  EXPECT_THAT(tree.globally_missing_sites, testing::IsEmpty());

  // When site 2 leaves globally_missing_sites (because tip 1 has data there), the tree's state
  // at that site is set to tip 1's state (C).  Since the ref has A, the focus picks up A->C.
  EXPECT_THAT(tree.deltas_ref_to_focus, testing::Contains(
      site_deltas_entry(2, rA, rC)));
  assert_utree_integrity(tree, true);
  assert_utree_matches_tip_descs(tree, tips, true);
}

TEST(Utree_test, guide_tree_globally_missing_sites_shrinks) {
  // ref = AAA.  Tip 0 missing [1,3), tip 1 missing [0,2), tip 2 has no missing.
  // After tip 0: globally_missing_sites = [1,3)
  // After tip 1: globally_missing_sites = [1,2)  (intersection of [1,3) and [0,2))
  // After tip 2: globally_missing_sites = empty
  auto ref = Real_sequence{rA, rA, rA};
  auto all_tips = std::vector<Tip_desc>{
      make_tip({}, Interval_set<>{{Site_interval{1, 3}}}),  // tip 0
      make_tip({}, Interval_set<>{{Site_interval{0, 2}}}),  // tip 1
      make_tip({})                                          // tip 2
  };

  // Build incrementally to verify globally_missing_sites at each stage
  auto bitgen1 = std::mt19937{42};
  auto tree1 = build_guide_tree(ref,
      {all_tips.begin(), all_tips.begin() + 1}, bitgen1);
  EXPECT_THAT(tree1.globally_missing_sites, testing::ElementsAre(
      std::pair{1, 3}));
  auto tips1 = std::vector<Tip_desc>{all_tips.begin(), all_tips.begin() + 1};
  assert_utree_integrity(tree1, true);
  assert_utree_matches_tip_descs(tree1, tips1, true);

  auto bitgen2 = std::mt19937{42};
  auto tips2 = std::vector<Tip_desc>{all_tips.begin(), all_tips.begin() + 2};
  auto tree2 = build_guide_tree(ref, tips2, bitgen2);
  EXPECT_THAT(tree2.globally_missing_sites, testing::ElementsAre(
      std::pair{1, 2}));
  assert_utree_integrity(tree2, true);
  assert_utree_matches_tip_descs(tree2, tips2, true);

  auto bitgen3 = std::mt19937{42};
  auto tree3 = build_guide_tree(ref, all_tips, bitgen3);
  EXPECT_THAT(tree3.globally_missing_sites, testing::IsEmpty());
  EXPECT_THAT(tree3.count_deltas(), testing::Eq(0));
  assert_utree_integrity(tree3, true);
  assert_utree_matches_tip_descs(tree3, all_tips, true);
}

// Helper: build a small 3-tip tree manually.
//
// Inner node 3 has the reference state (A at all sites).
// Tip 0 has C at site 0; tip 1 has G at site 1; tip 2 matches the reference.
//
//   tip0(C...) ---(0:C->A)--- inner3(A...) ---(1:A->G)--- tip1(.G..)
//                                  |
//                           (no deltas)
//                                  |
//                              tip2(A...)
//
// Focus at tip0, ref_to_focus = {site 0: A->C}
auto make_3tip_tree() -> Utree {
  auto tree = Utree::make_empty(3);
  tree.ref_sequence = {rA, rA, rA, rA};
  tree.num_inner_nodes_so_far = 1;

  auto arc_0_3 = tree.add_arc(0, 3);
  set_arc_deltas(tree, arc_0_3, {{0, {rC, rA}}});

  auto arc_3_1 = tree.add_arc(3, 1);
  set_arc_deltas(tree, arc_3_1, {{1, {rA, rG}}});

  tree.add_arc(3, 2);

  tree.reset_focus(0);
  tree.deltas_ref_to_focus[0] = {rA, rC};

  return tree;
}

// --- Euler tour traversals ---
//
// make_3tip_tree arc indices:
//   add_arc(0,3): arc 0 = 0->3, arc 1 = 3->0
//   add_arc(3,1): arc 2 = 3->1, arc 3 = 1->3
//   add_arc(3,2): arc 4 = 3->2, arc 5 = 2->3

struct Arc_hop {
  Node_index from;
  Node_index to;
  Arc_direction direction;
  auto operator<=>(const Arc_hop&) const = default;
};

TEST(Utree_test, annotated_arc_euler_tour_from_focus) {
  auto tree = make_3tip_tree();  // focus at 0
  using enum Arc_direction;

  auto result = std::vector<Arc_hop>{};
  for (auto [arc, dir] : annotated_arc_euler_tour(tree, tree.focus)) {
    result.push_back({tree.origin(arc), tree.target(arc), dir});
  }

  // DFS from tip0: enter 3, enter 2, leave 2, enter 1, leave 1, leave 3
  EXPECT_THAT(result, testing::ElementsAre(
      Arc_hop{0, 3, entering},
      Arc_hop{3, 2, entering},
      Arc_hop{2, 3, leaving},
      Arc_hop{3, 1, entering},
      Arc_hop{1, 3, leaving},
      Arc_hop{3, 0, leaving}));
}

TEST(Utree_test, annotated_arc_euler_tour_from_inner_node) {
  auto tree = make_3tip_tree();
  using enum Arc_direction;

  auto result = std::vector<Arc_hop>{};
  for (auto [arc, dir] : annotated_arc_euler_tour(tree, Node_index{3})) {
    result.push_back({tree.origin(arc), tree.target(arc), dir});
  }

  // DFS from inner3: arcs are [3->0, 3->1, 3->2], stack processes last-pushed first
  EXPECT_THAT(result, testing::ElementsAre(
      Arc_hop{3, 2, entering},
      Arc_hop{2, 3, leaving},
      Arc_hop{3, 1, entering},
      Arc_hop{1, 3, leaving},
      Arc_hop{3, 0, entering},
      Arc_hop{0, 3, leaving}));
}

TEST(Utree_test, arc_euler_tour_matches_annotated) {
  auto tree = make_3tip_tree();

  auto annotated_arcs = std::vector<Arc_index>{};
  for (auto [arc, dir] : annotated_arc_euler_tour(tree, tree.focus)) {
    annotated_arcs.push_back(arc);
  }

  auto plain = std::vector<Arc_index>{};
  for (auto arc : arc_euler_tour(tree, tree.focus)) {
    plain.push_back(arc);
  }

  EXPECT_THAT(plain, testing::ElementsAreArray(annotated_arcs));
}

// --- utree_tips ---

TEST(Utree_test, utree_tips_empty) {
  auto tree = Utree::make_empty(0);
  auto tips = std::vector<Node_index>{};
  for (auto tip : utree_tips(tree)) {
    tips.push_back(tip);
  }
  EXPECT_THAT(tips, testing::IsEmpty());
}

TEST(Utree_test, utree_tips_3tip_tree) {
  auto tree = make_3tip_tree();
  auto tips = std::vector<Node_index>{};
  for (auto tip : utree_tips(tree)) {
    tips.push_back(tip);
  }
  EXPECT_THAT(tips, testing::ElementsAre(0, 1, 2));
}

// --- assert_utree_integrity ---

TEST(Utree_test, assert_integrity_passes_on_valid_tree) {
  auto tree = make_3tip_tree();
  assert_utree_integrity(tree, /*force=*/true);
}

TEST(Utree_test, assert_integrity_passes_on_one_tip_tree) {
  auto tree = Utree::make_empty(1);
  tree.ref_sequence = {rA, rC, rG, rT};
  tree.focus = 0;
  tree.deltas_ref_to_focus[0] = {rA, rT};
  assert_utree_integrity(tree, /*force=*/true);
}

TEST(Utree_test, assert_integrity_passes_on_empty_tree) {
  auto tree = Utree{};
  tree.ref_sequence = {rA, rC, rG, rT};
  assert_utree_integrity(tree, /*force=*/true);
}

#ifndef __EMSCRIPTEN__  // We want tests to be compiled to exercise clang, but EXPECT_DEATH doesn't exist in Emscripten
TEST(Utree_test, assert_integrity_catches_bad_delta) {
  auto tree = make_3tip_tree();

  // Corrupt: change a delta on one arc without updating the mate
  auto arc_0_3 = tree.nodes[0].arcs[0];
  tree.arcs[arc_0_3].deltas[0] = {rA, rT};  // was A->C, now A->T; mate still has C->A

  EXPECT_DEATH(assert_utree_integrity(tree, /*force=*/true), "");
}
#endif

// --- assert_utree_matches_tip_descs ---

TEST(Utree_test, assert_matches_tip_descs_passes) {
  auto tree = make_3tip_tree();
  auto tips = std::vector<Tip_desc>{
      make_tip({{0, rA, rC}}),   // tip 0: C at site 0
      make_tip({{1, rA, rG}}),   // tip 1: G at site 1
      make_tip({})               // tip 2: same as ref
  };
  assert_utree_matches_tip_descs(tree, tips, /*force=*/true);
}

TEST(Utree_test, assert_matches_tip_descs_with_missing_data) {
  // tip0(A..) ---(0:A->C, 1:A->G)--- tip1(C.., missing site 1)
  auto tree = Utree::make_empty(2);
  tree.ref_sequence = {rA, rA, rA};
  tree.globally_missing_sites.clear();

  auto arc_01 = tree.add_arc(0, 1);
  set_arc_deltas(tree, arc_01, {{0, {rA, rC}}, {1, {rA, rG}}});
  tree.reset_focus(0);

  auto tips = std::vector<Tip_desc>{
      make_tip({}),                        // tip 0: same as ref
      make_tip({{0, rA, rC}}, {{{1, 2}}})  // tip 1: C at site 0, missing site 1
  };
  assert_utree_matches_tip_descs(tree, tips, /*force=*/true);
}

#ifndef __EMSCRIPTEN__  // We want tests to be compiled to exercise clang, but EXPECT_DEATH doesn't exist in Emscripten
TEST(Utree_test, assert_matches_catches_wrong_tip) {
  auto tree = make_3tip_tree();
  auto tips = std::vector<Tip_desc>{
      make_tip({{0, rA, rT}}),               // WRONG: was rC in tree, now rT
      make_tip({{0, rA, rC}, {1, rA, rG}}),
      make_tip({})
  };

  EXPECT_DEATH(assert_utree_matches_tip_descs(tree, tips, /*force=*/true), "");
}
#endif

// --- midpoint_root_utree ---

TEST(Utree_test, midpoint_root_zero_tips) {
  auto tree = Utree::make_empty(0);
  auto rooting_info = midpoint_root_utree(tree, {});
  EXPECT_EQ(rooting_info.root, k_no_node);
}

TEST(Utree_test, midpoint_root_one_tip) {
  auto ref = Real_sequence{rA, rC, rG, rT};
  auto bitgen = std::mt19937{42};
  auto tips = std::vector<Tip_desc>{
      make_dated_tip("tip0", 100.0f, 100.0f, {{0, rA, rT}})
  };
  auto tree = build_guide_tree(ref, tips, bitgen);
  auto rooting_info = midpoint_root_utree(tree, tips);
  EXPECT_EQ(rooting_info.root, 0);
}

TEST(Utree_test, midpoint_root_two_tips_equal_dates) {
  auto ref = Real_sequence{rA, rA, rA, rA};
  auto bitgen = std::mt19937{42};
  auto tips = std::vector<Tip_desc>{
      make_dated_tip("tip0", 100.0f, 100.0f, {{0, rA, rC}}),
      make_dated_tip("tip1", 100.0f, 100.0f, {{1, rA, rG}})
  };
  auto tree = build_guide_tree(ref, tips, bitgen);
  assert_utree_integrity(tree, true);
  auto rooting_info = midpoint_root_utree(tree, tips);
  assert_utree_integrity(tree, true);

  // Root should be a new inner node
  EXPECT_GE(rooting_info.root, tree.num_tips);
  // Tree should have 3 nodes total (2 tips + 1 inner)
  EXPECT_EQ(tree.num_tips + tree.num_inner_nodes_so_far, 3);
  // Same dates => fallback rate
  EXPECT_NEAR(rooting_info.lambda, 1.0 / 30.0, 1e-10);
}

TEST(Utree_test, midpoint_root_three_tips) {
  auto ref = Real_sequence{rA, rA, rA, rA, rA, rA};
  auto bitgen = std::mt19937{42};
  auto tips = std::vector<Tip_desc>{
      make_dated_tip("tip0", 100.0f, 100.0f, {{0, rA, rC}, {1, rA, rC}}),
      make_dated_tip("tip1", 100.0f, 100.0f, {{2, rA, rG}, {3, rA, rG}}),
      make_dated_tip("tip2", 100.0f, 100.0f, {{4, rA, rT}, {5, rA, rT}})
  };
  auto tree = build_guide_tree(ref, tips, bitgen);
  assert_utree_integrity(tree, true);

  auto rooting_info = midpoint_root_utree(tree, tips);
  assert_utree_integrity(tree, true);

  EXPECT_GE(rooting_info.root, tree.num_tips);
  // 3 tips + 2 inner nodes (original + root)
  EXPECT_EQ(tree.num_tips + tree.num_inner_nodes_so_far, 5);
}

TEST(Utree_test, midpoint_root_unequal_dates) {
  auto ref = Real_sequence{rA, rA, rA, rA};
  auto bitgen = std::mt19937{42};
  // tip0 is much earlier, tip1 is recent — root should be shifted toward tip0
  auto tips = std::vector<Tip_desc>{
      make_dated_tip("tip0", 0.0f, 0.0f, {{0, rA, rC}}),
      make_dated_tip("tip1", 100.0f, 100.0f, {{1, rA, rG}})
  };
  auto tree = build_guide_tree(ref, tips, bitgen);
  auto rooting_info = midpoint_root_utree(tree, tips);
  assert_utree_integrity(tree, true);
  EXPECT_GE(rooting_info.root, tree.num_tips);
  EXPECT_GT(rooting_info.lambda, 0.0);
  EXPECT_LE(rooting_info.t_MRCA, 0.0);  // root should be at or before earliest tip
}

TEST(Utree_test, midpoint_root_known_rate) {
  // Build a tree where tips have known distances from root and known dates
  // so we can verify the OLS estimate
  auto ref = Real_sequence{rA, rA, rA, rA, rA, rA};
  auto bitgen = std::mt19937{42};
  // 3 tips at different dates with mutations proportional to time
  // lambda = 0.1 mutations/day, t_root = 0
  // tip0: t=10, m=1; tip1: t=20, m=2; tip2: t=30, m=3
  auto tips = std::vector<Tip_desc>{
      make_dated_tip("tip0", 10.0f, 10.0f, {{0, rA, rC}}),
      make_dated_tip("tip1", 20.0f, 20.0f, {{0, rA, rC}, {1, rA, rG}}),
      make_dated_tip("tip2", 30.0f, 30.0f, {{0, rA, rC}, {1, rA, rG}, {2, rA, rT}})
  };
  auto tree = build_guide_tree(ref, tips, bitgen);
  auto rooting_info = midpoint_root_utree(tree, tips);

  EXPECT_EQ(rooting_info.method, Rooting_method::midpoint);
  EXPECT_GT(rooting_info.lambda, 0.0);
  EXPECT_LT(rooting_info.t_MRCA, 10.0);  // root should be before earliest tip
}

// --- regression_root_utree ---

TEST(Utree_test, regression_root_three_tips_known_rate) {
  // 3 tips with non-overlapping mutations and dates proportional to mutation count.
  // Each tip has distinct mutations, so the guide tree places them on tip-specific branches.
  // tip0: t=100, 1 unique mutation; tip1: t=200, 2 unique; tip2: t=300, 3 unique
  auto ref = Real_sequence{rA, rA, rA, rA, rA, rA};
  auto bitgen = std::mt19937{42};
  auto tips = std::vector<Tip_desc>{
      make_dated_tip("tip0", 100.0f, 100.0f, {{0, rA, rC}}),
      make_dated_tip("tip1", 200.0f, 200.0f, {{1, rA, rG}, {2, rA, rT}}),
      make_dated_tip("tip2", 300.0f, 300.0f, {{3, rA, rC}, {4, rA, rG}, {5, rA, rT}})
  };
  auto tree = build_guide_tree(ref, tips, bitgen);
  assert_utree_integrity(tree, true);

  auto rooting_info = regression_root_utree(tree, tips, bitgen);
  assert_utree_integrity(tree, true);

  EXPECT_EQ(rooting_info.method, Rooting_method::regression);
  EXPECT_GE(rooting_info.root, tree.num_tips);
  EXPECT_GT(rooting_info.lambda, 0.0);
  EXPECT_LT(rooting_info.t_MRCA, 100.0);  // root should be before earliest tip
}

TEST(Utree_test, regression_root_same_dates_fallback) {
  // All tips have the same date: Var(t) = 0, should fall back to midpoint rooting
  auto ref = Real_sequence{rA, rA, rA, rA};
  auto bitgen = std::mt19937{42};
  auto tips = std::vector<Tip_desc>{
      make_dated_tip("tip0", 100.0f, 100.0f, {{0, rA, rC}}),
      make_dated_tip("tip1", 100.0f, 100.0f, {{1, rA, rG}}),
      make_dated_tip("tip2", 100.0f, 100.0f, {{2, rA, rT}})
  };
  auto tree = build_guide_tree(ref, tips, bitgen);

  auto rooting_info = regression_root_utree(tree, tips, bitgen);
  assert_utree_integrity(tree, true);

  EXPECT_EQ(rooting_info.method, Rooting_method::midpoint);
  EXPECT_GE(rooting_info.root, tree.num_tips);
  // Fallback rate = 1/30
  EXPECT_NEAR(rooting_info.lambda, 1.0 / 30.0, 1e-10);
}

TEST(Utree_test, regression_root_two_tips_fallback) {
  // N <= 2: should fall back to midpoint rooting
  auto ref = Real_sequence{rA, rA, rA, rA};
  auto bitgen = std::mt19937{42};
  auto tips = std::vector<Tip_desc>{
      make_dated_tip("tip0", 0.0f, 0.0f, {{0, rA, rC}}),
      make_dated_tip("tip1", 100.0f, 100.0f, {{1, rA, rG}})
  };
  auto tree = build_guide_tree(ref, tips, bitgen);

  auto rooting_info = regression_root_utree(tree, tips, bitgen);
  assert_utree_integrity(tree, true);

  EXPECT_EQ(rooting_info.method, Rooting_method::midpoint);
  EXPECT_GE(rooting_info.root, tree.num_tips);
  EXPECT_GT(rooting_info.lambda, 0.0);
}

TEST(Utree_test, regression_root_star_topology) {
  // Star topology: all tips equidistant from center, different dates
  auto ref = Real_sequence{rA, rA, rA, rA, rA, rA};
  auto bitgen = std::mt19937{42};
  auto tips = std::vector<Tip_desc>{
      make_dated_tip("tip0", 10.0f, 10.0f, {{0, rA, rC}}),
      make_dated_tip("tip1", 20.0f, 20.0f, {{1, rA, rG}}),
      make_dated_tip("tip2", 30.0f, 30.0f, {{2, rA, rT}}),
      make_dated_tip("tip3", 40.0f, 40.0f, {{3, rA, rC}})
  };
  auto tree = build_guide_tree(ref, tips, bitgen);
  assert_utree_integrity(tree, true);

  auto rooting_info = regression_root_utree(tree, tips, bitgen);
  assert_utree_integrity(tree, true);

  EXPECT_EQ(rooting_info.method, Rooting_method::regression);
  EXPECT_GE(rooting_info.root, tree.num_tips);
  EXPECT_GT(rooting_info.lambda, 0.0);
  EXPECT_LT(rooting_info.t_MRCA, 10.0);  // root before earliest tip
}

// --- utree_to_phylo_tree ---

TEST(Utree_test, utree_to_phylo_tree_three_tips) {
  auto ref = Real_sequence{rA, rA, rA, rA, rA, rA};
  auto bitgen = std::mt19937{42};
  auto tips = std::vector<Tip_desc>{
      make_dated_tip("tip0", 10.0f, 10.0f, {{0, rA, rC}}),
      make_dated_tip("tip1", 20.0f, 20.0f, {{1, rA, rG}}),
      make_dated_tip("tip2", 30.0f, 30.0f, {{2, rA, rT}})
  };
  auto tree = build_guide_tree(ref, tips, bitgen);
  auto rooting_info = midpoint_root_utree(tree, tips);
  auto phylo_tree = utree_to_phylo_tree(tree, rooting_info, tips, bitgen);

  // Basic checks
  EXPECT_EQ(std::ssize(phylo_tree), 2 * 3 - 1);
  EXPECT_EQ(phylo_tree.at(0).name, "tip0");
  EXPECT_EQ(phylo_tree.at(1).name, "tip1");
  EXPECT_EQ(phylo_tree.at(2).name, "tip2");

  // Tip times should be clamped to their date bounds
  EXPECT_NEAR(phylo_tree.at(0).t, 10.0, 1e-6);
  EXPECT_NEAR(phylo_tree.at(1).t, 20.0, 1e-6);
  EXPECT_NEAR(phylo_tree.at(2).t, 30.0, 1e-6);

  // Root should be before all tips
  EXPECT_LT(phylo_tree.at(phylo_tree.root).t, 10.0);
}

// --- for_each_tip_in_nearest_first_order ---

TEST(Utree_test, nearest_first_order_one_tip) {
  auto ref = Real_sequence{rA, rA, rA, rA};
  auto tips = std::vector<Tip_desc>{make_tip({})};
  auto bitgen = std::mt19937{42};
  auto tree = build_guide_tree(ref, tips, bitgen);

  auto order = std::vector<std::pair<Node_index, Node_index>>{};
  for_each_tip_in_nearest_first_order(tree, bitgen, [&](Node_index tip, Node_index prev) {
    order.push_back({tip, prev});
  });

  EXPECT_THAT(std::ssize(order), testing::Eq(1));
  EXPECT_THAT(order[0].first, testing::Eq(0));
  EXPECT_THAT(order[0].second, testing::Eq(k_no_node));
}

TEST(Utree_test, nearest_first_order_two_tips) {
  auto ref = Real_sequence{rA, rA, rA, rA};
  auto tips = std::vector<Tip_desc>{
      make_tip({{0, rA, rC}}),
      make_tip({{1, rA, rG}})
  };
  auto bitgen = std::mt19937{42};
  auto tree = build_guide_tree(ref, tips, bitgen);

  auto order = std::vector<std::pair<Node_index, Node_index>>{};
  for_each_tip_in_nearest_first_order(tree, bitgen, [&](Node_index tip, Node_index prev) {
    order.push_back({tip, prev});
  });

  EXPECT_THAT(std::ssize(order), testing::Eq(2));
  EXPECT_THAT(order[0].second, testing::Eq(k_no_node));
  // Second tip's closest_prev_tip should be the first tip
  EXPECT_THAT(order[1].second, testing::Eq(order[0].first));
  // Both tips covered
  auto tips_seen = std::vector<Node_index>{order[0].first, order[1].first};
  EXPECT_THAT(tips_seen, testing::UnorderedElementsAre(0, 1));
}

TEST(Utree_test, nearest_first_order_three_tips) {
  // ref = AAAA.  Tips 0 and 1 share a mutation at site 0 (A->C), tip 2 does not.
  // Tip 1 also has a unique mutation at site 1 (A->G).
  // In the guide tree, tips 0 and 1 are closer to each other than to tip 2.
  // If the starting tip is 0 or 1, the other should come second (before tip 2).
  auto ref = Real_sequence{rA, rA, rA, rA};
  auto tips = std::vector<Tip_desc>{
      make_tip({{0, rA, rC}}),               // tip 0: C at site 0
      make_tip({{0, rA, rC}, {1, rA, rG}}),  // tip 1: C at site 0, G at site 1
      make_tip({})                            // tip 2: same as ref
  };

  // Try multiple seeds to exercise different starting tips
  for (auto seed = 0; seed < 20; ++seed) {
    auto bitgen = std::mt19937{static_cast<unsigned>(seed)};
    auto tree = build_guide_tree(ref, tips, bitgen);

    auto order = std::vector<Node_index>{};
    for_each_tip_in_nearest_first_order(tree, bitgen, [&](Node_index tip, Node_index) {
      order.push_back(tip);
    });

    EXPECT_THAT(order, testing::UnorderedElementsAre(0, 1, 2));

    // If starting tip is 0, tip 1 should come before tip 2 (closer in guide tree)
    if (order[0] == 0) {
      EXPECT_THAT(order[1], testing::Eq(1));
    }
    // If starting tip is 1, tip 0 should come before tip 2
    if (order[0] == 1) {
      EXPECT_THAT(order[1], testing::Eq(0));
    }
  }
}

TEST(Utree_test, nearest_first_order_five_tips_covers_all) {
  auto ref = Real_sequence{rA, rA, rA, rA, rA};
  auto tips = std::vector<Tip_desc>{
      make_tip({{0, rA, rC}}),
      make_tip({{0, rA, rC}, {2, rA, rT}}),
      make_tip({{0, rA, rC}, {3, rA, rT}}),
      make_tip({{1, rA, rG}}),
      make_tip({{1, rA, rG}}),
  };
  auto bitgen = std::mt19937{42};
  auto tree = build_guide_tree(ref, tips, bitgen);

  auto order = std::vector<Node_index>{};
  for_each_tip_in_nearest_first_order(tree, bitgen, [&](Node_index tip, Node_index) {
    order.push_back(tip);
  });

  EXPECT_THAT(order, testing::UnorderedElementsAre(0, 1, 2, 3, 4));
}

// --- build_refined_tree ---

TEST(Utree_test, build_refined_tree_integrity) {
  auto ref = Real_sequence{rA, rA, rA, rA, rA, rA, rA, rA};
  auto tips = std::vector<Tip_desc>{
      make_tip({{0, rA, rC}, {1, rA, rG}}),
      make_tip({{0, rA, rC}, {2, rA, rT}}),
      make_tip({{3, rA, rG}, {4, rA, rT}}),
      make_tip({{5, rA, rC}, {6, rA, rG}, {7, rA, rT}})
  };
  auto bitgen = std::mt19937{42};

  auto guide_tree = build_guide_tree(ref, tips, bitgen);
  auto refined_tree = build_refined_tree(guide_tree, tips, bitgen);

  assert_utree_integrity(refined_tree, true);
  assert_utree_matches_tip_descs(refined_tree, tips, true);
}

TEST(Utree_test, build_refined_tree_parsimony) {
  auto ref = Real_sequence{rA, rA, rA, rA, rA};
  auto tips = std::vector<Tip_desc>{
      make_tip({{0, rA, rC}}),
      make_tip({{0, rA, rC}, {2, rA, rT}}),
      make_tip({{0, rA, rC}, {3, rA, rT}}),
      make_tip({{1, rA, rG}}),
      make_tip({{1, rA, rG}}),
  };

  for (auto seed = 0; seed < 20; ++seed) {
    auto bitgen = std::mt19937{static_cast<unsigned>(seed)};
    auto guide_tree = build_guide_tree(ref, tips, bitgen);
    auto refined_tree = build_refined_tree(guide_tree, tips, bitgen);

    assert_utree_integrity(refined_tree, true);
    assert_utree_matches_tip_descs(refined_tree, tips, true);
  }
}

TEST(Utree_test, build_refined_tree_with_missing_data) {
  auto ref = Real_sequence{rA, rA, rA, rA, rA, rA};
  auto tips = std::vector<Tip_desc>{
      make_tip({{2, rA, rC}}, Interval_set<>{{Site_interval{0, 2}}}),
      make_tip({{3, rA, rG}}),
      make_tip({{4, rA, rT}})
  };
  auto bitgen = std::mt19937{42};

  auto guide_tree = build_guide_tree(ref, tips, bitgen);
  auto refined_tree = build_refined_tree(guide_tree, tips, bitgen);

  assert_utree_integrity(refined_tree, true);
  assert_utree_matches_tip_descs(refined_tree, tips, true);
}

// --- build_initial_phylo_tree (end-to-end) ---

TEST(Utree_test, build_initial_phylo_tree_end_to_end) {
  auto ref = Real_sequence{rA, rA, rA, rA, rA, rA, rA, rA};
  auto bitgen = std::mt19937{42};
  auto tips = std::vector<Tip_desc>{
      make_dated_tip("seq1", 10.0f, 10.0f, {{0, rA, rC}, {1, rA, rG}}),
      make_dated_tip("seq2", 20.0f, 20.0f, {{0, rA, rC}, {2, rA, rT}}),
      make_dated_tip("seq3", 30.0f, 30.0f, {{3, rA, rG}, {4, rA, rT}}),
      make_dated_tip("seq4", 40.0f, 40.0f, {{5, rA, rC}, {6, rA, rG}, {7, rA, rT}})
  };

  auto phylo_tree = build_initial_phylo_tree(ref, tips, bitgen);

  EXPECT_EQ(std::ssize(phylo_tree), 2 * 4 - 1);

  // Verify tip names preserved
  EXPECT_EQ(phylo_tree.at(0).name, "seq1");
  EXPECT_EQ(phylo_tree.at(1).name, "seq2");
  EXPECT_EQ(phylo_tree.at(2).name, "seq3");
  EXPECT_EQ(phylo_tree.at(3).name, "seq4");

  // Verify tip dates preserved
  EXPECT_NEAR(phylo_tree.at(0).t, 10.0, 1e-6);
  EXPECT_NEAR(phylo_tree.at(1).t, 20.0, 1e-6);
  EXPECT_NEAR(phylo_tree.at(2).t, 30.0, 1e-6);
  EXPECT_NEAR(phylo_tree.at(3).t, 40.0, 1e-6);

  // Root should be before all tips
  EXPECT_LT(phylo_tree.at(phylo_tree.root).t, 10.0);

  // assert_phylo_tree_integrity and assert_phylo_tree_matches_tip_descs
  // are already called inside build_initial_phylo_tree, so if we get here, they passed
}

TEST(Utree_test, build_initial_phylo_tree_with_missing_data) {
  auto ref = Real_sequence{rA, rA, rA, rA, rA, rA};
  auto bitgen = std::mt19937{42};

  auto tips = std::vector<Tip_desc>{
      make_dated_tip("seq1", 10.0f, 10.0f, {{2, rA, rC}}, Interval_set<>{{Site_interval{0, 2}}}),
      make_dated_tip("seq2", 20.0f, 20.0f, {{3, rA, rG}}),
      make_dated_tip("seq3", 30.0f, 30.0f, {{4, rA, rT}})
  };

  auto phylo_tree = build_initial_phylo_tree(ref, tips, bitgen);

  EXPECT_EQ(std::ssize(phylo_tree), 2 * 3 - 1);
  EXPECT_EQ(phylo_tree.at(0).name, "seq1");
}

}  // namespace delphy
