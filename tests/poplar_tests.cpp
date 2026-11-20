#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "poplar.h"

namespace delphy {

inline constexpr auto sA = Seq_letters::A;
inline constexpr auto sC = Seq_letters::C;
inline constexpr auto sG = Seq_letters::G;
inline constexpr auto sT = Seq_letters::T;
inline constexpr auto sN = Seq_letters::N;

inline constexpr auto rA = Real_seq_letter::A;
inline constexpr auto rC = Real_seq_letter::C;
inline constexpr auto rG = Real_seq_letter::G;
inline constexpr auto rT = Real_seq_letter::T;

class Poplar_test : public testing::Test {
 protected:
  //
  // The tree that we build:
  //
  //                                     +T0N--T1G--G1T--- a (NTNN)
  //                                     |
  //                     +A2N--A0T-A1T---+ x (TTNN)
  //                     |               |
  // (ACAA) A3N- C1A --+ r (AAAN)        +-------- T1G ------- b (TGNN)
  //                     |
  //                     +A1N--------A0T------- T0G ------------------ c (GNAN)
  // 

  Real_sequence ref_sequence{rA, rC, rA, rA};
  Phylo_tree tree{5};

  static constexpr Node_index r = 0;
  static constexpr Node_index x = 1;
  static constexpr Node_index a = 2;
  static constexpr Node_index b = 3;
  static constexpr Node_index c = 4;

  Poplar_test() {
    tree.root = r;
    tree.ref_sequence = ref_sequence;

    tree.at(r).parent = k_no_node;
    tree.at(r).children = {x, c};
    tree.at(r).name = "r";
    tree.at(r).t_min = -std::numeric_limits<float>::max();
    tree.at(r).t_max = +std::numeric_limits<float>::max();
    tree.at(r).t = -1.0;
    tree.at(r).mutations = {Mutation{rC, 1, rA, -std::numeric_limits<double>::max()}};
    tree.at(r).missations = {{Missation{3, rA}}, ref_sequence};

    tree.at(x).parent = r;
    tree.at(x).children = {a, b};
    tree.at(x).name = "x";
    tree.at(x).t_min = -std::numeric_limits<float>::max();
    tree.at(x).t_max = +std::numeric_limits<float>::max();
    tree.at(x).t = 0.0;
    tree.at(x).mutations = {Mutation{rA, 0, rT, -0.5}, Mutation{rA, 1, rT, -0.4}};
    tree.at(x).missations = {{Missation{2, rA}}, ref_sequence};

    tree.at(a).parent = x;
    tree.at(a).children = {};
    tree.at(a).name = "a";
    tree.at(a).t = tree.at(a).t_min = tree.at(a).t_max = 1.0;
    tree.at(a).mutations = {Mutation{rT, 1, rG, 0.3}, Mutation{rG, 1, rT, 0.6}};
    tree.at(a).missations = {{Missation{0, rT}}, ref_sequence};

    tree.at(b).parent = x;
    tree.at(b).children = {};
    tree.at(b).name = "b";
    tree.at(b).t = tree.at(b).t_min = tree.at(b).t_max = 2.0;
    tree.at(b).mutations = {Mutation{rT, 1, rG, 1.0}};

    tree.at(c).parent = r;
    tree.at(c).children = {};
    tree.at(c).name = "c";
    tree.at(c).t = tree.at(c).t_min = tree.at(c).t_max = 3.0;
    tree.at(c).mutations = {Mutation{rA, 0, rT, 0.0}, Mutation{rT, 0, rG, 1.0}};
    tree.at(c).missations = {{Missation{1, rA}}, ref_sequence};

    assert_phylo_tree_integrity(tree);
  }
};

class Poplar_test2 : public testing::Test {
 protected:
  //
  // The tree that we build here (== example tree in poplar.h top-level comments):
  //
  //                    +---A1G---+ a (CGA)
  //                    |
  //         +----A0C---+ x
  //         |          |
  //         |          +(--------+ b (CNA)
  //         |          A1N
  //    -----+ r (AAA)     
  //         |
  //         |            +--A2T--+ c (AAT)
  //         |            |
  //         +------------+ y
  //                      |
  //                      |
  //                      +-------+ d (AAA)
  // 

  Real_sequence ref_sequence{rA, rA, rA};
  Phylo_tree tree{7};

  static constexpr Node_index r = 0;
  static constexpr Node_index x = 1;
  static constexpr Node_index y = 2;
  static constexpr Node_index a = 3;
  static constexpr Node_index b = 4;
  static constexpr Node_index c = 5;
  static constexpr Node_index d = 6;

  Poplar_test2() {
    tree.root = r;
    tree.ref_sequence = ref_sequence;

    tree.at(r).parent = k_no_node;
    tree.at(r).children = {x, y};
    tree.at(r).name = "r";
    tree.at(r).t_min = -std::numeric_limits<float>::max();
    tree.at(r).t_max = +std::numeric_limits<float>::max();
    tree.at(r).t = -1.0;
    tree.at(r).mutations = {};
    tree.at(r).missations = {};

    tree.at(x).parent = r;
    tree.at(x).children = {a, b};
    tree.at(x).name = "x";
    tree.at(x).t_min = -std::numeric_limits<float>::max();
    tree.at(x).t_max = +std::numeric_limits<float>::max();
    tree.at(x).t = 0.0;
    tree.at(x).mutations = {Mutation{rA, 0, rC, -0.5}};
    tree.at(x).missations = {};

    tree.at(y).parent = r;
    tree.at(y).children = {c, d};
    tree.at(y).name = "y";
    tree.at(y).t_min = -std::numeric_limits<float>::max();
    tree.at(y).t_max = +std::numeric_limits<float>::max();
    tree.at(y).t = 0.2;
    tree.at(y).mutations = {};
    tree.at(y).missations = {};

    tree.at(a).parent = x;
    tree.at(a).children = {};
    tree.at(a).name = "a";
    tree.at(a).t = tree.at(a).t_min = tree.at(a).t_max = 1.0;
    tree.at(a).mutations = {Mutation{rA, 1, rG, 0.5}};
    tree.at(a).missations = {};

    tree.at(b).parent = x;
    tree.at(b).children = {};
    tree.at(b).name = "b";
    tree.at(b).t = tree.at(b).t_min = tree.at(b).t_max = 1.0;
    tree.at(b).mutations = {};
    tree.at(b).missations = {{Missation{1, rA}}, ref_sequence};

    tree.at(c).parent = y;
    tree.at(c).children = {};
    tree.at(c).name = "c";
    tree.at(c).t = tree.at(c).t_min = tree.at(c).t_max = 1.0;
    tree.at(c).mutations = {Mutation{rA, 2, rT, 0.5}};
    tree.at(c).missations = {};

    tree.at(d).parent = y;
    tree.at(d).children = {};
    tree.at(d).name = "d";
    tree.at(d).t = tree.at(d).t_min = tree.at(d).t_max = 1.0;
    tree.at(d).mutations = {};
    tree.at(d).missations = {};

    assert_phylo_tree_integrity(tree);
  }
};

static auto seq_to_str(const Sequence& seq) -> std::string {
  auto result = std::string{};
  for (auto a : seq) {
    result.push_back(to_char(a));
  }
  return result;
}

TEST(Poplar_tests, site_state_fingerprint_randomish) {
  // Run through all possible sequences with 3 letters ((2^4 - 1)^3 = 3375) and ensure that they all have
  // different fingerprints
  auto possible_seq_letters = std::vector<Seq_letter>{};
  for (auto can_A : {false, true}) {
    for (auto can_C : {false, true}) {
      for (auto can_G : {false, true}) {
        for (auto can_T : {false, true}) {
          auto seq_letter = 
              (can_A ? Seq_letters::A : Seq_letters::_) |
              (can_C ? Seq_letters::C : Seq_letters::_) |
              (can_G ? Seq_letters::G : Seq_letters::_) |
              (can_T ? Seq_letters::T : Seq_letters::_) |
              Seq_letters::_;

          // Skip "ambiguous" Seq_letter where no state is possible
          if (seq_letter == Seq_letters::_) { continue; }
          
          possible_seq_letters.push_back(seq_letter);
        }
      }
    }
  }
      
  auto fingerprints_2_seqs = absl::flat_hash_map<uint64_t, Sequence>{};
  for (auto s_0 : possible_seq_letters) {
    for (auto s_1 : possible_seq_letters) {
      for (auto s_2 : possible_seq_letters) {
        
        auto sequence = Sequence{{s_0, s_1, s_2}};
        auto fingerprint =
            site_state_fingerprint(0, s_0) ^ 
            site_state_fingerprint(1, s_1) ^ 
            site_state_fingerprint(2, s_2);

        // An early version of this erroneously sign-extended 32-bit fingerprints
        EXPECT_THAT(fingerprint >> 32, testing::Ne(0x00000000));
        EXPECT_THAT(fingerprint >> 32, testing::Ne(0xFFFFFFFF));
        
        auto [it, inserted] = fingerprints_2_seqs.try_emplace(fingerprint, sequence);

        EXPECT_THAT(inserted, testing::IsTrue()) <<
            absl::StreamFormat("Fingerprint collision %s and %s -> %d",
                               seq_to_str(it->second),
                               seq_to_str(sequence),
                               fingerprint);
      }
    }
  }

  // Visual check of actual fingerprints (revealed 32-bit sign extension bug!)
  // for (const auto& [fingerprint, seq] : fingerprints_2_seqs) {
  //   std::cerr << absl::StreamFormat("%16ul <- %s\n", fingerprint, seq_to_str(seq));
  // }
}

TEST_F(Poplar_test, map_site_tree_branch) {
  auto is = [](Node_index start_node, Node_index end_node) { return testing::Eq(std::make_pair(start_node, end_node)); };

  // If a full "missing" check is enabled in the code (flip full_missing_site_check to
  // `true`), then map_site_tree_branch conclusively determines whether a full-branch is
  // missing in the specified pruned site tree and returns {k_no_node, k_no_node}.
  // Otherwise, this check is best-effort to avoid incurring needless costs during normal
  // operation
  auto maybe_is = [](Node_index start_node, Node_index end_node) {
    return testing::AnyOf(
        testing::Eq(std::make_pair(start_node, end_node)),
        testing::Eq(std::make_pair(k_no_node, k_no_node)));
  };
  
  EXPECT_THAT(map_site_tree_branch(tree, r, 0), is(k_no_node, k_no_node));
  EXPECT_THAT(map_site_tree_branch(tree, x, 0), is(c, b));
  EXPECT_THAT(map_site_tree_branch(tree, a, 0), is(k_no_node, k_no_node));
  EXPECT_THAT(map_site_tree_branch(tree, b, 0), is(c, b));
  EXPECT_THAT(map_site_tree_branch(tree, c, 0), is(b, c));
  
  EXPECT_THAT(map_site_tree_branch(tree, r, 1), is(k_no_node, k_no_node));
  EXPECT_THAT(map_site_tree_branch(tree, x, 1), is(k_no_node, x));
  EXPECT_THAT(map_site_tree_branch(tree, a, 1), maybe_is(x, a));   // Can't tell locally that site 1 is missing
  EXPECT_THAT(map_site_tree_branch(tree, b, 1), maybe_is(x, b));   // Can't tell locally that site 1 is missing
  EXPECT_THAT(map_site_tree_branch(tree, c, 1), is(k_no_node, k_no_node));
  
  EXPECT_THAT(map_site_tree_branch(tree, r, 2), is(k_no_node, k_no_node));
  EXPECT_THAT(map_site_tree_branch(tree, x, 2), is(k_no_node, k_no_node));
  EXPECT_THAT(map_site_tree_branch(tree, a, 2), maybe_is(x, a));  // Can't tell locally that site 2 is missing
  EXPECT_THAT(map_site_tree_branch(tree, b, 2), maybe_is(x, b));  // Can't tell locally that site 2 is missing
  EXPECT_THAT(map_site_tree_branch(tree, c, 2), is(k_no_node, c));
  
  EXPECT_THAT(map_site_tree_branch(tree, r, 3), is(k_no_node, k_no_node));
  EXPECT_THAT(map_site_tree_branch(tree, x, 3), is(k_no_node, k_no_node));
  EXPECT_THAT(map_site_tree_branch(tree, a, 3), maybe_is(x, a));  // Can't tell locally that site 3 is missing
  EXPECT_THAT(map_site_tree_branch(tree, b, 3), maybe_is(x, b));  // Can't tell locally that site 3 is missing
  EXPECT_THAT(map_site_tree_branch(tree, c, 3), is(k_no_node, k_no_node));
}

TEST_F(Poplar_test, map_mutation_extent) {
  EXPECT_THAT(map_mutation_extent(tree, r, 0), testing::Eq(std::nullopt));  // "branch" above root doesn't exist
  EXPECT_THAT(map_mutation_extent(tree, x, 0), testing::Eq(Mutation_extent{0, b, rT, c, rG, r}));  // Wraps around root AND hits Edge case (a) =2+ mutations with a net effect
  EXPECT_THAT(map_mutation_extent(tree, a, 0), testing::Eq(std::nullopt));  // site 0 is missing along x-a
  EXPECT_THAT(map_mutation_extent(tree, b, 0), testing::Eq(Mutation_extent{0, b, rT, c, rG, r}));
  EXPECT_THAT(map_mutation_extent(tree, c, 0), testing::Eq(Mutation_extent{0, b, rT, c, rG, r}));
  
  EXPECT_THAT(map_mutation_extent(tree, r, 1), testing::Eq(std::nullopt));  // "branch" above root doesn't exist
  EXPECT_THAT(map_mutation_extent(tree, x, 1), testing::Eq(Mutation_extent{1, k_no_node, rA, x, rT, k_no_node}));  // Edge case (c): mutation above all informative samples
  EXPECT_THAT(map_mutation_extent(tree, a, 1), testing::Eq(Mutation_extent{1, x, rT, a, rT, x}));  // Edge case (b): 2+ mutations with no net effect => "neutral" extent
  EXPECT_THAT(map_mutation_extent(tree, b, 1), testing::Eq(Mutation_extent{1, x, rT, b, rG, x}));
  EXPECT_THAT(map_mutation_extent(tree, c, 1), testing::Eq(std::nullopt));  // site 1 is missing along r-c
  
  EXPECT_THAT(map_mutation_extent(tree, r, 2), testing::Eq(std::nullopt));  // no mutations on site 2
  EXPECT_THAT(map_mutation_extent(tree, x, 2), testing::Eq(std::nullopt));  // no mutations on site 2
  EXPECT_THAT(map_mutation_extent(tree, a, 2), testing::Eq(std::nullopt));  // no mutations on site 2
  EXPECT_THAT(map_mutation_extent(tree, b, 2), testing::Eq(std::nullopt));  // no mutations on site 2
  EXPECT_THAT(map_mutation_extent(tree, c, 2), testing::Eq(std::nullopt));  // no mutations on site 2
  
  EXPECT_THAT(map_mutation_extent(tree, r, 3), testing::Eq(std::nullopt));  // no mutations on site 3
  EXPECT_THAT(map_mutation_extent(tree, x, 3), testing::Eq(std::nullopt));  // no mutations on site 3
  EXPECT_THAT(map_mutation_extent(tree, a, 3), testing::Eq(std::nullopt));  // no mutations on site 3
  EXPECT_THAT(map_mutation_extent(tree, b, 3), testing::Eq(std::nullopt));  // no mutations on site 3
  EXPECT_THAT(map_mutation_extent(tree, c, 3), testing::Eq(std::nullopt));  // no mutations on site 3
}

TEST_F(Poplar_test, map_tree_mutation_extents_AND_annotate_branches_with_extents) {
  // Mutation extents are sorted in the order implied by the default implementation of `<=>`
  auto mapped = map_tree_mutation_extents(tree);
  EXPECT_THAT(mapped.extents, testing::ElementsAre(
      Mutation_extent{0, b, rT, c, rG, r},    // extent 0
      Mutation_extent{1, x, rT, b, rG, x}));  // extent 1

  // The "half" extent (edge case (c)) on site 1 (A at r, T at x) is not an extent; it only forces the root blob to have T
  EXPECT_THAT(mapped.root_state_overrides, testing::ElementsAre(Root_state_override{1, rT}));

  EXPECT_THAT(annotate_branches_with_extents(tree, mapped.extents), testing::Eq(Branch_vector<Branch_extent_info>{
      {{},  {}},     // *-r branch
      {{},  {}},     // r-x branch
      {{},  {}},     // x-a branch
      {{1}, {0, 1}}, // x-b branch
      {{},  {0}}})); // r-c branch
}

TEST_F(Poplar_test, derive_blob_tree) {
  //
  // Here's the test tree, with extent boundaries marked:
  //
  //          "Half" extent (not a real extent: root state override)
  //          (site 1, A at r, T at x)   +T0N--T1G--G1T--- a (NTNN)
  //                     ===============*|
  //                    /+A2N--A0T-A1T---+ x (TTNN)
  //              ...==/ |  /=========\  |
  // (ACAA) A3N- C1A --+ r H           H +-------- T1G ------- b (TGNN)
  //              (AAAN) | H            \====================>  Extent 0  (site 0, T at b, G at c)
  //                     | H              <==================>  Extent 1  (site 1, T at x, G at b)
  //                     | H
  //                     |  \========================================>
  //                     +A1N--------A0T------- T0G ------------------ c (GNAN)
  //
  // So the tree splits up into blobs like this (note that the "half" extent doesn't create blob boundaries):
  //
  //                                     +T0N--T1G--G1T--- a (NTNN)  [0]
  //                                     |
  //                     +A2N--A0T-A1T---+ x (TTNN)
  //                     |
  // (ACAA) A3N- C1A --+ r
  //              (AAAN) |                      +-------- T1G ------- [1]      b (TGNN)  [2]
  //                     |
  //                     |
  //                     +A1N--------A0T------- T0G ------------------         c (GNAN)  [3]
  //
  //                                                                  
  //
  // More specifically, here are the 4 blobs in their full glory:
  //
  //     Parent   Children     Blob Sequence            Nodes(*)   Branches(*)     Active Mut. Exts.
  //     -------  --------     -------------------      -----      -------------   -----------------
  // [0]    -     [1],[3]      {T,G}  T    A    N       *r,x,a     r-x, x-a, r-c   0
  // [1]   [0]    [2]          {T,G}{T,G}  A    N       -          *x-b            0, 1
  // [2]   [1]    -              T    G    A    N       *b         -               -
  // [3]   [0]    -              G    T    A    N       *c         -               -
  //
  // (*) seed node/branch marked with `*`
  //
  // Note that there's no guarantee about the blob indices in the calculated blob tree,
  // so the test has to be written with care.
  
  auto blob_tree = derive_blob_tree(tree);

  assert_tree_integrity(blob_tree);

  // Match blob tree topology and extract actual indices (e.g., `[0]` above has index `b0`)
  ASSERT_THAT(std::ssize(blob_tree), testing::Eq(4));  // If this isn't true, the below could well crash
  auto b0 = blob_tree.root;
  EXPECT_THAT(b0, testing::Ne(k_no_blob));
  EXPECT_THAT(blob_tree.at(b0).children, testing::SizeIs(2));
  
  auto b1 = k_no_blob;
  for (const auto& child : blob_tree.at(b0).children) {
    for (const auto& branch : blob_tree.at(child).branches) {
      if (branch == b) { // x-b
        b1 = child;
      }
    }
  }
  EXPECT_THAT(b1, testing::Ne(k_no_blob));
  EXPECT_THAT(blob_tree.at(b1).children, testing::SizeIs(1));
  
  auto b2 = blob_tree.at(b1).children[0];
  EXPECT_THAT(blob_tree.at(b2).children, testing::IsEmpty());
  
  auto b3 = k_no_blob;
  for (const auto& child : blob_tree.at(b0).children) {
    for (const auto& node : blob_tree.at(child).nodes) {
      if (node == c) {
        b3 = child;
      }
    }
  }
  EXPECT_THAT(b3, testing::Ne(k_no_blob));
  EXPECT_THAT(blob_tree.at(b3).children, testing::IsEmpty());
  EXPECT_THAT((std::set<Blob_index>{b0, b1, b2, b3}), testing::SizeIs(4));  // test uniqueness

  // Check nodes & branches
  EXPECT_THAT(blob_tree.at(b0).nodes, testing::ElementsAre(r, x, a));
  EXPECT_THAT(blob_tree.at(b0).branches, testing::ElementsAre(x, a, c));
  EXPECT_THAT(blob_tree.at(b1).nodes, testing::ElementsAre());
  EXPECT_THAT(blob_tree.at(b1).branches, testing::ElementsAre(b));
  EXPECT_THAT(blob_tree.at(b2).nodes, testing::ElementsAre(b));
  EXPECT_THAT(blob_tree.at(b2).branches, testing::ElementsAre());
  EXPECT_THAT(blob_tree.at(b3).nodes, testing::ElementsAre(c));
  EXPECT_THAT(blob_tree.at(b3).branches, testing::ElementsAre());

  struct Case { Node_index node; Blob_index blob; };  // Also used for branches
  for (const auto& [node, blob] : {Case{r, b0}, Case{x, b0}, Case{a, b0}, Case{b, b2}, Case{c, b3}}) {
    EXPECT_THAT(blob_tree.node_2_blob[node], testing::Eq(blob)) << node << " " << blob;
  }
  for (const auto& [branch, blob] : {Case{r, k_no_blob}, Case{x, b0}, Case{a, b0}, Case{b, b1}, Case{c, b0}}) {
    EXPECT_THAT(blob_tree.branch_2_blob[branch], testing::Eq(blob)) << branch << " " << blob;
  }

  // Check seed nodes & branches
  EXPECT_THAT(blob_tree.at(b0).seed_node, testing::Eq(r));
  EXPECT_THAT(blob_tree.at(b0).seed_branch, testing::Eq(k_no_branch));
  EXPECT_THAT(blob_tree.at(b1).seed_node, testing::Eq(k_no_node));
  EXPECT_THAT(blob_tree.at(b1).seed_branch, testing::Eq(b));
  EXPECT_THAT(blob_tree.at(b2).seed_node, testing::Eq(b));
  EXPECT_THAT(blob_tree.at(b2).seed_branch, testing::Eq(k_no_branch));
  EXPECT_THAT(blob_tree.at(b3).seed_node, testing::Eq(c));
  EXPECT_THAT(blob_tree.at(b3).seed_branch, testing::Eq(k_no_branch));

  // Check mutation extents
  EXPECT_THAT(blob_tree.extents, testing::ElementsAre(
      Mutation_extent{0, b, rT, c, rG, r},    // extent 0
      Mutation_extent{1, x, rT, b, rG, x}));  // extent 1
  EXPECT_THAT(blob_tree.root_state_overrides, testing::ElementsAre(Root_state_override{1, rT}));  // from "half" extent

  // Sites missing in the whole tree are recorded separately and never show up in sequence deltas.  The root blob's
  // deltas (vs. the reference sequence ACAA) only have the ambiguity at site 0 and the override at site 1
  EXPECT_THAT(blob_tree.missing_sites, testing::ElementsAre(Site_interval{3, 4}));
  EXPECT_THAT(blob_tree.at(b0).sequence_delta, testing::ElementsAre(
      Blob_site_delta{0, sA, sT | sG},
      Blob_site_delta{1, sC, sT}));
  EXPECT_THAT(blob_tree.at(b0).active_mutation_extents, testing::ElementsAre(0));
  EXPECT_THAT(blob_tree.at(b1).active_mutation_extents, testing::ElementsAre(0, 1));
  EXPECT_THAT(blob_tree.at(b2).active_mutation_extents, testing::ElementsAre());
  EXPECT_THAT(blob_tree.at(b3).active_mutation_extents, testing::ElementsAre());

  // Check fingerprints
  auto fp = [](Seq_letter s0, Seq_letter s1, Seq_letter s2, Seq_letter s3) {
    return
        site_state_fingerprint(0, s0) ^
        site_state_fingerprint(1, s1) ^
        site_state_fingerprint(2, s2) ^
        site_state_fingerprint(3, s3);
  };
  EXPECT_THAT(blob_tree.at(b0).fingerprint, testing::Eq(fp(sT | sG,   sT   , sA, sN)));
  EXPECT_THAT(blob_tree.at(b1).fingerprint, testing::Eq(fp(sT | sG, sT | sG, sA, sN)));
  EXPECT_THAT(blob_tree.at(b2).fingerprint, testing::Eq(fp(   sT  ,   sG   , sA, sN)));
  EXPECT_THAT(blob_tree.at(b3).fingerprint, testing::Eq(fp(   sG  ,   sT   , sA, sN)));

  // Check blob sequences & corresponding fingerprints (ref seq is ACAA)
  auto b0_sequence = query_blob_sequence(blob_tree, b0);
  EXPECT_THAT(b0_sequence.missing_sites, testing::ElementsAre(Site_interval{3, 4}));
  EXPECT_THAT(b0_sequence.mutations, testing::ElementsAre(
      Blob_sequence_mutation{1, rC, rT}));
  EXPECT_THAT(b0_sequence.ambiguities, testing::ElementsAre(
      Blob_sequence_ambiguity{0, rG, rT}));
  EXPECT_THAT(calc_blob_fingerprint(tree.ref_sequence, b0_sequence), testing::Eq(blob_tree.at(b0).fingerprint));

  auto b1_sequence = query_blob_sequence(blob_tree, b1);
  EXPECT_THAT(b1_sequence.missing_sites, testing::ElementsAre(Site_interval{3, 4}));
  EXPECT_THAT(b1_sequence.mutations, testing::ElementsAre());
  EXPECT_THAT(b1_sequence.ambiguities, testing::ElementsAre(
      Blob_sequence_ambiguity{0, rG, rT},
      Blob_sequence_ambiguity{1, rG, rT}));
  EXPECT_THAT(calc_blob_fingerprint(tree.ref_sequence, b1_sequence), testing::Eq(blob_tree.at(b1).fingerprint));

  auto b2_sequence = query_blob_sequence(blob_tree, b2);
  EXPECT_THAT(b2_sequence.missing_sites, testing::ElementsAre(Site_interval{3, 4}));
  EXPECT_THAT(b2_sequence.mutations, testing::ElementsAre(
        Blob_sequence_mutation{0, rA, rT},
        Blob_sequence_mutation{1, rC, rG}));
  EXPECT_THAT(b2_sequence.ambiguities, testing::ElementsAre());
  EXPECT_THAT(calc_blob_fingerprint(tree.ref_sequence, b2_sequence), testing::Eq(blob_tree.at(b2).fingerprint));

  auto b3_sequence = query_blob_sequence(blob_tree, b3);
  EXPECT_THAT(b3_sequence.missing_sites, testing::ElementsAre(Site_interval{3, 4}));
  EXPECT_THAT(b3_sequence.mutations, testing::ElementsAre(
        Blob_sequence_mutation{0, rA, rG},
        Blob_sequence_mutation{1, rC, rT}));
  EXPECT_THAT(b3_sequence.ambiguities, testing::ElementsAre());
  EXPECT_THAT(calc_blob_fingerprint(tree.ref_sequence, b3_sequence), testing::Eq(blob_tree.at(b3).fingerprint));
}

TEST_F(Poplar_test2, derive_blob_tree_2) {
  //
  // From the example at the top of poplar.h, here's the test tree with extent boundaries marked:
  //
  //                  /==========>
  //         =========  +---A1G---+ a (CGA)
  //       / ==========>|
  //      / /+----A0C---+ x
  //     H H |          |
  //     H H |          +(--------+ b (CNA)
  //     H H |          A1N
  //    -H-H-+ r (AAA)     
  //     H H |             <=====>   Extent 2  (site 2, A at y, T at c)
  //     H H |            +--A2T--+ c (AAT)
  //     H H |            |
  //      \ \+------------+ y
  //       \ ============>|      Extent 0  (site 0, C at x, A at y)
  //        =============>|      Extent 1  (site 1, G at a, A at y)
  //                      +-------+ d (AAA)
  //
  // So the tree splits up into blobs like this:
  //
  //                                          [2]
  //                          +---A1G---      + a (CGA)
  //                          |
  //         +----A0C---      + x  [1]
  //         |                |
  //     [0] |                |
  //         |               +(--------+ b (CNA)
  //         |                A1N
  //    -----+ r (AAA)     
  //         |                        [4]       [5]
  //         |                      --A2T--      + c (AAT)
  //         |                
  //         +------------     + y
  //                           |
  //                       [3] |
  //                           +-------+ d (AAA)
  //
  // More specifically, here are the 6 blobs, their parents and children,
  // their member branches & nodes, and their "sequence" (see below)
  //
  //     Parent   Children     Blob Sequence        Nodes(*)   Branches(*)   Active Mut. Exts.
  //     -------  -----------  ---------------      --------   -----------   -----------------
  // [0]    -     [1],[3]      {A,C}{A,G}  A        *r         r-x, r-y      0, 1
  // [1]   [0]    [2]            C  {A,G}  A        *x,b       x-b, x-a      1
  // [2]   [1]    -              C    G    A        *a         -             -
  // [3]   [0]    [4]            A    A    A        *y,d       y-d           -
  // [4]   [3]    [5]            A    A  {A,T}      -          *y-c          2
  // [5]   [4]    -              A    A    T        *c         -             -
  //
  // (*) seed node/branch marked with `*`
  //
  // Note that there's no guarantee about the blob indices in the calculated blob tree,
  // so the test has to be written with care.
  
  auto blob_tree = derive_blob_tree(tree);

  assert_tree_integrity(blob_tree);
  
  // Match blob tree topology and extract actual indices (e.g., `[0]` above has index `b0`)
  ASSERT_THAT(std::ssize(blob_tree), testing::Eq(6));  // If this isn't true, the below could well crash
  auto b0 = blob_tree.root;
  EXPECT_THAT(b0, testing::Ne(k_no_blob));
  EXPECT_THAT(blob_tree.at(b0).children, testing::SizeIs(2));
  
  auto b1 = k_no_blob;
  for (const auto& child : blob_tree.at(b0).children) {
    for (const auto& branch : blob_tree.at(child).branches) {
      if (branch == b) { // x-b
        b1 = child;
      }
    }
  }
  EXPECT_THAT(b1, testing::Ne(k_no_blob));
  EXPECT_THAT(blob_tree.at(b1).children, testing::SizeIs(1));
  
  auto b2 = k_no_blob;
  for (const auto& child : blob_tree.at(b1).children) {
    for (const auto& node : blob_tree.at(child).nodes) {
      if (node == a) { // a
        b2 = child;
      }
    }
  }
  EXPECT_THAT(b2, testing::Ne(k_no_blob));
  EXPECT_THAT(blob_tree.at(b2).children, testing::IsEmpty());
  
  auto b3 = k_no_blob;
  for (const auto& child : blob_tree.at(b0).children) {
    for (const auto& branch : blob_tree.at(child).branches) {
      if (branch == d) { // y-d
        b3 = child;
      }
    }
  }
  EXPECT_THAT(b3, testing::Ne(k_no_blob));
  EXPECT_THAT(blob_tree.at(b3).children, testing::SizeIs(1));
  
  auto b4 = blob_tree.at(b3).children[0];
  EXPECT_THAT(blob_tree.at(b4).children, testing::SizeIs(1));
  
  auto b5 = blob_tree.at(b4).children[0];
  EXPECT_THAT(blob_tree.at(b5).children, testing::IsEmpty());
  
  EXPECT_THAT((std::set<Blob_index>{b0, b1, b2, b3, b4, b5}), testing::SizeIs(6));  // test uniqueness

  // Check nodes & branches
  EXPECT_THAT(blob_tree.at(b0).nodes, testing::ElementsAre(r));
  EXPECT_THAT(blob_tree.at(b0).branches, testing::ElementsAre(x, y));
  EXPECT_THAT(blob_tree.at(b1).nodes, testing::ElementsAre(x, b));
  EXPECT_THAT(blob_tree.at(b1).branches, testing::ElementsAre(a, b));
  EXPECT_THAT(blob_tree.at(b2).nodes, testing::ElementsAre(a));
  EXPECT_THAT(blob_tree.at(b2).branches, testing::ElementsAre());
  EXPECT_THAT(blob_tree.at(b3).nodes, testing::ElementsAre(y, d));
  EXPECT_THAT(blob_tree.at(b3).branches, testing::ElementsAre(d));
  EXPECT_THAT(blob_tree.at(b4).nodes, testing::ElementsAre());
  EXPECT_THAT(blob_tree.at(b4).branches, testing::ElementsAre(c));
  EXPECT_THAT(blob_tree.at(b5).nodes, testing::ElementsAre(c));
  EXPECT_THAT(blob_tree.at(b5).branches, testing::ElementsAre());

  struct Case { Node_index node; Blob_index blob; };  // Also used for branches
  for (const auto& [node, blob] : {Case{r, b0}, Case{x, b1}, Case{y, b3}, Case{a, b2}, Case{b, b1}, Case{c, b5}, Case{d, b3}}) {
    EXPECT_THAT(blob_tree.node_2_blob[node], testing::Eq(blob)) << node << " " << blob;
  }
  for (const auto& [branch, blob] : {Case{r, k_no_node}, Case{x, b0}, Case{y, b0}, Case{a, b1}, Case{b, b1}, Case{c, b4}, Case{d, b3}}) {
    EXPECT_THAT(blob_tree.branch_2_blob[branch], testing::Eq(blob)) << branch << " " << blob;
  }

  // Check seed nodes & branches
  EXPECT_THAT(blob_tree.at(b0).seed_node, testing::Eq(r));
  EXPECT_THAT(blob_tree.at(b0).seed_branch, testing::Eq(k_no_branch));
  EXPECT_THAT(blob_tree.at(b1).seed_node, testing::Eq(x));
  EXPECT_THAT(blob_tree.at(b1).seed_branch, testing::Eq(k_no_branch));
  EXPECT_THAT(blob_tree.at(b2).seed_node, testing::Eq(a));
  EXPECT_THAT(blob_tree.at(b2).seed_branch, testing::Eq(k_no_branch));
  EXPECT_THAT(blob_tree.at(b3).seed_node, testing::Eq(y));
  EXPECT_THAT(blob_tree.at(b3).seed_branch, testing::Eq(k_no_branch));
  EXPECT_THAT(blob_tree.at(b4).seed_node, testing::Eq(k_no_node));
  EXPECT_THAT(blob_tree.at(b4).seed_branch, testing::Eq(c));
  EXPECT_THAT(blob_tree.at(b5).seed_node, testing::Eq(c));
  EXPECT_THAT(blob_tree.at(b5).seed_branch, testing::Eq(k_no_branch));

  // Check mutation extents
  EXPECT_THAT(blob_tree.extents, testing::ElementsAre(
      Mutation_extent{0, x, rC, y, rA, r},    // extent 0
      Mutation_extent{1, y, rA, a, rG, r},    // extent 1
      Mutation_extent{2, y, rA, c, rT, y}));  // extent 2
  EXPECT_THAT(blob_tree.at(b0).active_mutation_extents, testing::ElementsAre(0, 1));
  EXPECT_THAT(blob_tree.at(b1).active_mutation_extents, testing::ElementsAre(1));
  EXPECT_THAT(blob_tree.at(b2).active_mutation_extents, testing::ElementsAre());
  EXPECT_THAT(blob_tree.at(b3).active_mutation_extents, testing::ElementsAre());
  EXPECT_THAT(blob_tree.at(b4).active_mutation_extents, testing::ElementsAre(2));
  EXPECT_THAT(blob_tree.at(b5).active_mutation_extents, testing::ElementsAre());

  // Check branch extents
  EXPECT_THAT(blob_tree.branch_extents, testing::Eq(Branch_vector<Branch_extent_info>{
      {{},  {}},     // *-r branch
      {{},  {0}},    // r-x branch
      {{},  {0, 1}}, // r-y branch
      {{},  {1}},    // x-a branch
      {{},  {}},     // x-b branch
      {{2}, {2}},    // y-c branch
      {{},  {}}}));   // y-d branch
  
  // Check fingerprints
  auto fp = [](Seq_letter s0, Seq_letter s1, Seq_letter s2) {
    return
        site_state_fingerprint(0, s0) ^
        site_state_fingerprint(1, s1) ^
        site_state_fingerprint(2, s2);
  };
  EXPECT_THAT(blob_tree.at(b0).fingerprint, testing::Eq(fp(sA | sC, sA | sG,   sA   )));
  EXPECT_THAT(blob_tree.at(b1).fingerprint, testing::Eq(fp(   sC  , sA | sG,   sA   )));
  EXPECT_THAT(blob_tree.at(b2).fingerprint, testing::Eq(fp(   sC  ,   sG,      sA   )));
  EXPECT_THAT(blob_tree.at(b3).fingerprint, testing::Eq(fp(   sA  ,   sA   ,   sA   )));
  EXPECT_THAT(blob_tree.at(b4).fingerprint, testing::Eq(fp(   sA  ,   sA   , sA | sT)));
  EXPECT_THAT(blob_tree.at(b5).fingerprint, testing::Eq(fp(   sA  ,   sA   ,   sT   )));

  // Check blob sequences (ref seq is AAA)
  auto b0_sequence = query_blob_sequence(blob_tree, b0);
  EXPECT_THAT(b0_sequence.missing_sites, testing::ElementsAre());
  EXPECT_THAT(b0_sequence.mutations, testing::ElementsAre());
  EXPECT_THAT(b0_sequence.ambiguities, testing::ElementsAre(
      Blob_sequence_ambiguity{0, rA, rC},
      Blob_sequence_ambiguity{1, rA, rG}));
  EXPECT_THAT(calc_blob_fingerprint(tree.ref_sequence, b0_sequence), testing::Eq(blob_tree.at(b0).fingerprint));
  
  auto b1_sequence = query_blob_sequence(blob_tree, b1);
  EXPECT_THAT(b1_sequence.missing_sites, testing::ElementsAre());
  EXPECT_THAT(b1_sequence.mutations, testing::ElementsAre(
        Blob_sequence_mutation{0, rA, rC}));
  EXPECT_THAT(b1_sequence.ambiguities, testing::ElementsAre(
      Blob_sequence_ambiguity{1, rA, rG}));
  EXPECT_THAT(calc_blob_fingerprint(tree.ref_sequence, b1_sequence), testing::Eq(blob_tree.at(b1).fingerprint));
  
  auto b2_sequence = query_blob_sequence(blob_tree, b2);
  EXPECT_THAT(b2_sequence.missing_sites, testing::ElementsAre());
  EXPECT_THAT(b2_sequence.mutations, testing::ElementsAre(
        Blob_sequence_mutation{0, rA, rC},
        Blob_sequence_mutation{1, rA, rG}));
  EXPECT_THAT(b2_sequence.ambiguities, testing::ElementsAre());
  EXPECT_THAT(calc_blob_fingerprint(tree.ref_sequence, b2_sequence), testing::Eq(blob_tree.at(b2).fingerprint));
  
  auto b3_sequence = query_blob_sequence(blob_tree, b3);
  EXPECT_THAT(b3_sequence.missing_sites, testing::ElementsAre());
  EXPECT_THAT(b3_sequence.mutations, testing::ElementsAre());
  EXPECT_THAT(b3_sequence.ambiguities, testing::ElementsAre());
  EXPECT_THAT(calc_blob_fingerprint(tree.ref_sequence, b3_sequence), testing::Eq(blob_tree.at(b3).fingerprint));
  
  auto b4_sequence = query_blob_sequence(blob_tree, b4);
  EXPECT_THAT(b4_sequence.missing_sites, testing::ElementsAre());
  EXPECT_THAT(b4_sequence.mutations, testing::ElementsAre());
  EXPECT_THAT(b4_sequence.ambiguities, testing::ElementsAre(
      Blob_sequence_ambiguity{2, rA, rT}));
  EXPECT_THAT(calc_blob_fingerprint(tree.ref_sequence, b4_sequence), testing::Eq(blob_tree.at(b4).fingerprint));
  
  auto b5_sequence = query_blob_sequence(blob_tree, b5);
  EXPECT_THAT(b5_sequence.missing_sites, testing::ElementsAre());
  EXPECT_THAT(b5_sequence.mutations, testing::ElementsAre(
        Blob_sequence_mutation{2, rA, rT}));
  EXPECT_THAT(b5_sequence.ambiguities, testing::ElementsAre());
  EXPECT_THAT(calc_blob_fingerprint(tree.ref_sequence, b5_sequence), testing::Eq(blob_tree.at(b5).fingerprint));
}

// Regression tests for query_blob_sequence() when a blob boundary (caused by site 0) lies in the
// middle of a path that is "unconstrained" at site 1.  The blob sequence must not depend on the
// intermediate mutations along such a path.

class Poplar_test3 : public testing::Test {
 protected:
  //
  // Half-extent at site 1 whose path r->x->y->e crosses a blob boundary (site 0 extent on branch x-y):
  //
  //                          +(--------- z (AN)       <- z missing site 1
  //                          |
  //                +---------+ x (AA)
  //                |         |                  +(------ w (TN)   <- w missing site 1
  //                |         +--A0T--A1C--------+ y (TC)
  //     -----------+ r (AA)                     |
  //                |                            +--C1G---+ e (TG)
  //                +(--------- s (AN)                     |
  //                                                       +-- p (TG)
  //                                                       +-- q (TG)
  //
  // Site 1 is missing at s, z and w, so the pruned site-1 tree has the single branch r-e, and
  // all intermediate mutations (A1C on x-y, C1G on y-e) are "unconstrained": the state at e is
  // G, and so is the blob sequence at every blob below the half-extent's end node e.
  //
  Real_sequence ref_sequence{rA, rA};
  Phylo_tree tree{9};

  static constexpr Node_index r = 0;
  static constexpr Node_index x = 1;
  static constexpr Node_index s = 2;
  static constexpr Node_index z = 3;
  static constexpr Node_index y = 4;
  static constexpr Node_index w = 5;
  static constexpr Node_index e = 6;
  static constexpr Node_index p = 7;
  static constexpr Node_index q = 8;

  auto set_inner(
      Node_index n, Node_index parent, std::initializer_list<Node_index> child_list, const char* name, double t)
      -> void {
    CHECK_EQ(std::ssize(child_list), 2);
    tree.at(n).parent = parent;
    tree.at(n).children = {*child_list.begin(), *(child_list.begin() + 1)};
    tree.at(n).name = name;
    tree.at(n).t_min = -std::numeric_limits<float>::max();
    tree.at(n).t_max = +std::numeric_limits<float>::max();
    tree.at(n).t = t;
  }

  auto set_tip(Node_index n, Node_index parent, const char* name, double t) -> void {
    tree.at(n).parent = parent;
    tree.at(n).children = {};
    tree.at(n).name = name;
    tree.at(n).t = tree.at(n).t_min = tree.at(n).t_max = t;
  }

  Poplar_test3() {
    tree.root = r;
    tree.ref_sequence = ref_sequence;

    set_inner(r, k_no_node, {x, s}, "r", -1.0);
    set_inner(x, r, {y, z}, "x", 0.0);
    set_tip(s, r, "s", 1.0);
    tree.at(s).missations = {{Missation{1, rA}}, ref_sequence};
    set_tip(z, x, "z", 1.0);
    tree.at(z).missations = {{Missation{1, rA}}, ref_sequence};
    set_inner(y, x, {e, w}, "y", 0.5);
    tree.at(y).mutations = {Mutation{rA, 0, rT, 0.1}, Mutation{rA, 1, rC, 0.2}};
    set_tip(w, y, "w", 1.0);
    tree.at(w).missations = {{Missation{1, rC}}, ref_sequence};
    set_inner(e, y, {p, q}, "e", 1.0);
    tree.at(e).mutations = {Mutation{rC, 1, rG, 0.7}};
    set_tip(p, e, "p", 2.0);
    set_tip(q, e, "q", 2.0);

    assert_phylo_tree_integrity(tree);
  }
};

class Poplar_test4 : public Poplar_test3 {
 protected:
  //
  // Same idea, but with a net-zero (edge case (b)) path at site 1 that spans 2 full branches
  // (x->y->e, since w is missing at site 1) and has a blob boundary (site 0 extent on branch
  // y) in the middle:
  //
  //                          +-------- z (AA)
  //                +---------+ x (AA)
  //                |         |                  +(------ w (TN)   <- w missing site 1
  //                |         +--A0T--A1C--------+ y (TC)
  //     -----------+ r (AA)                     |
  //                |                            +--C1A---+ e (TA)
  //                +--------- s (AA)                     |
  //                                                      +-- p (TA)
  //                                                      +-- q (TA)
  //
  // The pruned site-1 tree has a single branch x-e, with A at both ends, so there's no mutation
  // extent at site 1 and the blob sequence at every blob must have A at site 1.
  //
  Poplar_test4() {
    tree.at(s).missations = {};
    tree.at(z).missations = {};
    tree.at(e).mutations = {Mutation{rC, 1, rA, 0.7}};
    assert_phylo_tree_integrity(tree);
  }
};

// A deliberately naive, brute-force reference for "what is the sequence at this position of the tree?",
// written straight from the definition of pruned site trees and sharing no logic with poplar.cpp.
//
// For each site l, the "pruned site tree" is the unrooted tree spanned by the tips that aren't missing
// at l.  Its vertices are the informative tips and the nodes with 2 informative children; every other
// node (including the root, whose two sides get fused or, if only one side is informative, left
// dangling) lies in the interior of a pruned branch.  The sequence at a position is then:
//
//  - the actual state there, if the position is a vertex of the pruned site tree;
//  - otherwise, a position in the interior of a pruned branch, whose state is the common state at
//    the branch's two end-point vertices if they agree (whatever happens along the way, e.g.,
//    mutations that cancel out), or the ambiguity {a, b} if they differ;
//  - for the dangling branch above the "first" vertex when the root has only 1 informative child (edge
//    case (c)), the state at that vertex;
//  - positions whose subtree has no informative tips (offshoots with missations) borrow the state of
//    the position they hang from;
//  - N if the site is missing in the entire tree.
struct Naive_position {
  Node_index node;
  bool is_branch;  // false: the node itself; true: the branch above the node
};

struct Naive_site_trees {
  const Phylo_tree& tree;
  Node_vector<Site_vector<int>> num_informative_tips;       // [node][l]
  Node_vector<Real_sequence> state;                         // [node]; actual sequence at each node in the full tree

  explicit Naive_site_trees(const Phylo_tree& tree_in) : tree{tree_in} {
    auto num_sites = tree.num_sites();
    num_informative_tips.assign(std::ssize(tree), Site_vector<int>(num_sites, 0));
    state.assign(std::ssize(tree), Real_sequence(num_sites, Real_seq_letter::A));

    for (auto l = Site_index{0}; l != num_sites; ++l) {
      auto missing = Node_vector<bool>(std::ssize(tree), false);
      for (const auto& node : pre_order_traversal(tree)) {
        auto parent = tree.at(node).parent;
        missing[node] = (parent != k_no_node && missing[parent]) || tree.at(node).missations.contains(l);
        auto s = (parent == k_no_node) ? tree.ref_sequence[l] : state[parent][l];
        for (const auto& m : tree.at(node).mutations) {
          if (m.site == l) { s = m.to; }
        }
        state[node][l] = s;
      }
      for (const auto& node : post_order_traversal(tree)) {
        if (tree.at(node).is_tip()) {
          num_informative_tips[node][l] = missing[node] ? 0 : 1;
        } else {
          auto n = 0;
          for (const auto& child : tree.at(node).children) { n += num_informative_tips[child][l]; }
          num_informative_tips[node][l] = n;
        }
      }
    }
  }

  auto expected_sequence(Naive_position pos) const -> Sequence {
    auto result = Sequence{};
    for (auto l = Site_index{0}; l != tree.num_sites(); ++l) {
      result.push_back(expected_state(l, pos));
    }
    return result;
  }

  // The vertices of the unrooted pruned site tree for site `l` are the informative tips and the nodes with
  // 2 informative children.  In particular, the root is never a vertex.
  auto is_vertex_in_unrooted_pruned_site_tree(Site_index l, Node_index n) const -> bool {
    if (n == tree.root || num_informative_tips[n][l] == 0) { return false; }
    if (tree.at(n).is_tip()) { return true; }
    return
        num_informative_tips[tree.at(n).left_child()][l] > 0 &&
        num_informative_tips[tree.at(n).right_child()][l] > 0;
  }

  // Return an arbitrary child of node `n` that has informative downstream tips at site `l`
  auto any_informative_child(Site_index l, Node_index n) const -> Node_index {
    for (const auto& child : tree.at(n).children) {
      if (num_informative_tips[child][l] > 0) { return child; }
    }
    CHECK(false) << "Node " << n << " has no informative child at site " << l;
  }

  // Walk down from `n` (an informative node that is not the root) to the first vertex of the unrooted pruned
  // site tree.  Every node on the way down that isn't a vertex has exactly 1 informative child.
  auto first_vertex_at_or_below(Site_index l, Node_index n) const -> Node_index {
    CHECK_NE(n, tree.root);
    while (not is_vertex_in_unrooted_pruned_site_tree(l, n)) { n = any_informative_child(l, n); }
    return n;
  }

  // State in the interior of the pruned branch that contains the full-tree branch above informative node `c`
  auto branch_state(Site_index l, Node_index c) const -> Seq_letter {
    auto bottom = first_vertex_at_or_below(l, c);
    auto prev = c;
    auto p = tree.at(c).parent;
    auto top = k_no_node;
    while (true) {
      if (p == tree.root) {
        auto other = tree.at(p).sibling_of(prev);
        if (num_informative_tips[other][l] == 0) {
          return to_seq_letter(state[bottom][l]);  // Dangling branch above the root (edge case (c))
        }
        top = first_vertex_at_or_below(l, other);  // Fused across the root
        break;
      }
      if (is_vertex_in_unrooted_pruned_site_tree(l, p)) {
        top = p;
        break;
      }
      prev = p;
      p = tree.at(p).parent;
    }
    return to_seq_letter(state[top][l]) | to_seq_letter(state[bottom][l]);
  }

  auto expected_state(Site_index l, Naive_position pos) const -> Seq_letter {
    if (num_informative_tips[tree.root][l] == 0) { return Seq_letters::N; }

    auto n = pos.node;
    if (num_informative_tips[n][l] == 0) {
      // `n` is in an offshoot with no informative tips (e.g., because of missations), so it takes on the state of
      // the point where the offshoot hangs off the pruned site tree.  That is the nearest informative ancestor of `n`.
      // It can't be a vertex of the unrooted pruned site tree (the offshoot isn't informative, so it has at most
      // 1 informative child), so it lies in the interior of a pruned branch, and we're looking for the state there.
      n = tree.at(n).parent;
      while (num_informative_tips[n][l] == 0) { n = tree.at(n).parent; }
      pos.is_branch = false;
    }

    if (n == tree.root) {
      // The root is never a vertex.  Its position lies in the interior of the pruned branch above any one of its
      // informative children (if both are informative, they're on the same pruned branch, fused across the root)
      return branch_state(l, any_informative_child(l, n));
    }
    if (not pos.is_branch && is_vertex_in_unrooted_pruned_site_tree(l, n)) {
      return to_seq_letter(state[n][l]);
    }
    return branch_state(l, n);
  }
};

static auto fingerprint_of_sequence(const Sequence& seq) -> Blob_fingerprint {
  auto result = k_empty_blob_fingerprint;
  for (auto l = Site_index{0}; l != std::ssize(seq); ++l) {
    result ^= site_state_fingerprint(l, seq[l]);
  }
  return result;
}

// Check the whole blob tree against the naive reference:
//  - every node and branch in a blob must have the blob's fingerprint, i.e., a blob really is a region
//    over which the (ambiguity-aware) sequence is constant;
//  - the sequence reported by query_blob_sequence() for the blob must also agree.
static auto expect_blob_tree_matches_naive_reference(const Blob_tree& blob_tree, const Phylo_tree& tree) -> void {
  auto reference = Naive_site_trees{tree};

  for (auto blob = Blob_index{0}; blob != std::ssize(blob_tree); ++blob) {
    const auto& b = blob_tree.at(blob);

    auto check_member = [&](Naive_position pos, const char* what) {
      auto expected = reference.expected_sequence(pos);
      EXPECT_THAT(b.fingerprint, testing::Eq(fingerprint_of_sequence(expected)))
          << "Blob " << blob << ", " << what << " " << pos.node << ": naive reference says "
          << seq_to_str(expected);
    };
    for (const auto& node : b.nodes) { check_member({node, false}, "node"); }
    for (const auto& branch : b.branches) { check_member({branch, true}, "branch above"); }

    auto seed_pos =
        (b.seed_node != k_no_node)
        ? Naive_position{b.seed_node, false}
        : Naive_position{b.seed_branch, true};
    auto expected = reference.expected_sequence(seed_pos);
    auto queried = query_blob_sequence(blob_tree, blob);
    EXPECT_THAT(calc_blob_fingerprint(tree.ref_sequence, queried), testing::Eq(fingerprint_of_sequence(expected)))
        << "Blob " << blob << ": query_blob_sequence disagrees with the naive reference, which says "
        << seq_to_str(expected);
  }
}

// For every blob, the fingerprint derived incrementally in derive_blob_tree() must agree with
// the fingerprint recalculated from the blob's queried sequence.
static auto expect_blob_fingerprints_consistent(const Blob_tree& blob_tree, const Phylo_tree& tree) -> void {
  for (auto blob = Blob_index{0}; blob != std::ssize(blob_tree); ++blob) {
    auto sequence = query_blob_sequence(blob_tree, blob);
    EXPECT_THAT(calc_blob_fingerprint(tree.ref_sequence, sequence), testing::Eq(blob_tree.at(blob).fingerprint))
        << "Blob " << blob << ": " << blob_tree.at(blob);
  }
}

TEST_F(Poplar_test3, half_extent_path_crossing_blob_boundary) {
  auto blob_tree = derive_blob_tree(tree);
  assert_tree_integrity(blob_tree);

  EXPECT_THAT(blob_tree.extents, testing::ElementsAre(
      Mutation_extent{0, x, rA, y, rT, x}));                   // extent 0: site 0 on branch x-y
  EXPECT_THAT(blob_tree.root_state_overrides, testing::ElementsAre(
      Root_state_override{1, rG}));                            // "half" extent r->e at site 1: root blob has G

  // The site-0 extent splits off blob 1 (branch y only) and blob 2 (node y and everything below)
  ASSERT_THAT(std::ssize(blob_tree), testing::Eq(3));
  auto b_y = blob_tree.node_2_blob[y];
  EXPECT_THAT(b_y, testing::Ne(blob_tree.root));
  EXPECT_THAT(blob_tree.node_2_blob[e], testing::Eq(b_y));

  // Everything at or below the half-extent's end has G at site 1, regardless of the A1C mutation on x-y
  auto sequence = query_blob_sequence(blob_tree, b_y);
  EXPECT_THAT(sequence.mutations, testing::ElementsAre(
      Blob_sequence_mutation{0, rA, rT},
      Blob_sequence_mutation{1, rA, rG}));

  expect_blob_fingerprints_consistent(blob_tree, tree);
}

TEST_F(Poplar_test4, net_zero_path_crossing_blob_boundary) {
  auto blob_tree = derive_blob_tree(tree);
  assert_tree_integrity(blob_tree);

  EXPECT_THAT(blob_tree.extents, testing::ElementsAre(
      Mutation_extent{0, x, rA, y, rT, x}));  // Only site 0; site 1 is net-zero

  ASSERT_THAT(std::ssize(blob_tree), testing::Eq(3));
  auto b_y = blob_tree.node_2_blob[y];
  EXPECT_THAT(b_y, testing::Ne(blob_tree.root));

  // Site 1 is A at both ends of the pruned site tree branch x->e, so the blob has A there (i.e., no mutation)
  auto sequence = query_blob_sequence(blob_tree, b_y);
  EXPECT_THAT(sequence.mutations, testing::ElementsAre(
      Blob_sequence_mutation{0, rA, rT}));

  expect_blob_fingerprints_consistent(blob_tree, tree);
}

class Poplar_test5 : public testing::Test {
 protected:
  //
  // Parallel mutations on both sides of the root (a net-zero path *through* the root):
  //
  //                 +---A0C--- x (C)
  //                 |
  //      ---------- r (A)
  //                 |
  //                 +---A0C--- y (C)
  //
  // The pruned site tree has a single branch x-y whose end-points agree (C), so there is no mutation
  // extent and no ambiguity, and the root blob has C at site 0, not the root's actual state A.
  //
  Real_sequence ref_sequence{rA};
  Phylo_tree tree{3};

  static constexpr Node_index r = 0;
  static constexpr Node_index x = 1;
  static constexpr Node_index y = 2;

  Poplar_test5() {
    tree.root = r;
    tree.ref_sequence = ref_sequence;

    tree.at(r).parent = k_no_node;
    tree.at(r).children = {x, y};
    tree.at(r).name = "r";
    tree.at(r).t_min = -std::numeric_limits<float>::max();
    tree.at(r).t_max = +std::numeric_limits<float>::max();
    tree.at(r).t = -1.0;

    tree.at(x).parent = r;
    tree.at(x).children = {};
    tree.at(x).name = "x";
    tree.at(x).t = tree.at(x).t_min = tree.at(x).t_max = 1.0;
    tree.at(x).mutations = {Mutation{rA, 0, rC, 0.0}};

    tree.at(y).parent = r;
    tree.at(y).children = {};
    tree.at(y).name = "y";
    tree.at(y).t = tree.at(y).t_min = tree.at(y).t_max = 1.0;
    tree.at(y).mutations = {Mutation{rA, 0, rC, 0.0}};

    assert_phylo_tree_integrity(tree);
  }
};

TEST_F(Poplar_test5, neutral_path_through_root) {
  auto blob_tree = derive_blob_tree(tree);
  assert_tree_integrity(blob_tree);

  ASSERT_THAT(std::ssize(blob_tree), testing::Eq(1));
  auto sequence = query_blob_sequence(blob_tree, blob_tree.root);
  EXPECT_THAT(sequence.mutations, testing::ElementsAre(Blob_sequence_mutation{0, rA, rC}));
  EXPECT_THAT(sequence.ambiguities, testing::ElementsAre());
  expect_blob_tree_matches_naive_reference(blob_tree, tree);
}

TEST_F(Poplar_test, matches_naive_reference) {
  expect_blob_tree_matches_naive_reference(derive_blob_tree(tree), tree);
}

TEST_F(Poplar_test2, matches_naive_reference) {
  expect_blob_tree_matches_naive_reference(derive_blob_tree(tree), tree);
}

TEST_F(Poplar_test3, matches_naive_reference) {
  expect_blob_tree_matches_naive_reference(derive_blob_tree(tree), tree);
}

TEST_F(Poplar_test4, matches_naive_reference) {
  expect_blob_tree_matches_naive_reference(derive_blob_tree(tree), tree);
}

}  // namespace delphy
