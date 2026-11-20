#ifndef DELPHY_POPLAR_H_
#define DELPHY_POPLAR_H_

#include <cstdint>
#include <optional>
#include <ostream>
#include <utility>
#include <vector>

#include "absl/log/check.h"
#include "absl/strings/str_format.h"
#include "absl/strings/str_join.h"

#include "interval_set.h"
#include "phylo_tree.h"
#include "sequence.h"
#include "tree.h"

namespace delphy {

// In POPLAR (POPulation-based LAyout of tRees), we partition a Phylo_tree's branches and
// nodes into "blobs", subtrees of the main tree whose sequence is almost determined, up
// to a few sites with unresolved ambiguity.  By construction, all possible resolutions of
// a blob's sequence have comparable posterior mass.
//
// To identify the blobs, we walk through each tree's mutations.  For each mutation, we
// trace out its "extent", the path along which it can easily move.  In most cases, that's
// just the branch on which the mutation lives, reflecting that the timing of the mutation
// along the branch is very uncertain:
//
//                       +--...                               +--...     <===> = Mutation extent
//          <===========>|                       <===========>|
//         +--A0C--------+                      +--------A0C--+
//         |             |                      |             |
//    ...--+             +--...     <--->  ...--+             +--...
//         |                                    |
//         +--...                               +--...
//
// Sometimes, however, a mutation can extend over multiple branches if the offshoots
// coming off the extent have missations at that site (so the subtree's state at that site
// can be imputed with either the `from` or `to` states with comparable ease):
//
//                           +--...                                   +--...
//                    ======>|                                 ======>|
//                   /+------+                                /+--A0C-+
//          <======== |      |                       <======== |      |
//         +--A0C-----+      +--...                 +----------+      +--...
//         |          |                             |          |
//    ...--+          +(--...           <--->  ...--+          +(--...
//         |          C0N                           |          A0N
//         +--...                                   +--...
//
// Those paths can even (easily) go around the root:
//
//                    +--...                              +--...
//         ==========>|                        ==========>|
//        /+----A0C---+                       /+----------+
//       H |          |                      H |          |
//       H |          +--...                 H |          +--...
//    ---H-+ r (A...)              <--->  ---H-+ r (C...)
//       H |            +--...               H |            +--...
//       H |            |                    H |            |
//        \+------------+                     \+-----C0A----+
//         ============>|                      ============>|
//                      +--...                              +--...
//
// Four noteworthy edge cases:
//
// a) Two mutations when one will do (everything is possible!)
//
//                       +--...                               +--...
//          <===========>|                       <===========>|
//         +-A0T---T0C---+                      +-----A0T-T0C-+
//         |             |                      |             |
//    ...--+             +--...     <--->  ...--+             +--...
//         |                                    |
//         +--...                               +--...
//
//    We treat this as a single A0C mutation covering the whole extent.  In spirit, we
//    normalize every branch at each pruned site-tree to have 0 or 1 mutations, and that
//    stands in for every possible mutation path with equal (resp., unequal) endpoints.
//
// b) Two mutations when zero are needed (see above):
//
//                       +--...                               +--...
//          <===========>|                                    |
//         +-A0T---T0A---+                      +-------------+
//         |             |                      |             |
//    ...--+             +--...     <--->  ...--+             +--...
//         |                                    |
//         +--...                               +--...
//
//    Here, we don't record a mutation extent in the final blob tree (see comment above).
//
// c) Even weirder, because of how we treat the root asymmetrically wrt missations, an
//    extent can be unbounded at the "above-the-root" node if only one child subtree of
//    the root is informative at a given site:
//
//                    +--...                              +--...
//         ==========*|                                   |
//        /+----A0C---+                        +----------+
//       H |          |                        |          |
// ...==/  |          +--...                   |          +--...
//    -----+ r (A...)              <--->  -----+ r (C...)
//         |                                   |
//         +(--...                             +(--...
//         A0N                                 C0N
//
//    Here, we don't record a mutation extent either, since we think of a path from the
//    root to a node whose offshoots all have missations at a given site to stand for "any
//    unconstrained mutational history that ends in a particular state" (above, C at
//    site 0), without specifying the details.  The only effect is that the root blob's
//    sequence has the state at the end of that path (C at site 0 in the right panel, even
//    though the root's actual state is A in the left panel), which we record as a
//    `Root_state_override`.
//
// d) A similar subtlety arises when the pruned site tree's branch goes through the root
//    and has the same state at both end-points even though the root's actual state is
//    different (e.g., the same mutation A0C on both children of the root).  There is no
//    ambiguity and no mutation extent, but the root blob's state is the common state at
//    the end-points (C), which again we record as a `Root_state_override`.
//
//                    +--...                              +--...
//         ==========>|                                   |
//        /+----A0C---+ (C...)                 +----------+ (C...)
//       H |          |                        |          |
//       H |          +--...                   |          +--...
//    ---H-+ r (A...)              <--->  -----+ r (C...)
//       H |            +--...                 |            +--...
//       H |            |                      |            |
//        \+----A0C-----+ (C...)               +------------+ (C...)
//         ============>|                                   |
//                      +--...                              +--...
//
// After the process of mapping out mutation extents is complete, the result is 0+ extent
// "boundaries" being placed at the beginning and end of every branch.  We then partition
// the tree's branches and nodes into "blobs" along these extent boundaries.  Note that
// blobs thus consist of connected branches and nodes: the nodes at the endpoint of a
// branch may be part of the same blob, or they may be in different blobs if a mutation
// extent boundary exists at the beginning/end of the branch.
//
// Example: the tree below has mutation extent boundaries at:
//  - the end of the r-x branch (A0C) - extent 0
//  - the end of the r-y branch (A0C and A1G) - extents 0 and 1
//  - the end of the x-a branch (A1G) - extent 1
//  - the start of the y-c branch (A2T) - extent 2
//  - the end of the y-c branch (A2T) - extent 2
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
// In order to match blobs with equal sequence across posterior trees, we assign each blob
// a `uint64_t` "fingerprint".  These are computed by xor'ing together the fingerprints
// for the possible states at each site.  The fingerprint for a state `a` at site `l` is
// given by `site_state_fingerprint(l, a)`, and just needs to be a good hash of the tuple
// (l, a).  Claude Sonnet 5.5 suggested using SplitMix64, as implemented in
// https://prng.di.unimi.it/splitmix.c; it's a 64-bit-to-64-bit bijection with good
// avalanche [i.e., flipping one bit in the input flips about half the bits in the
// output], and that seemed good enough for this purpose.
//
// When a blob's defining sequence `{s_l}` is completely unambiguous, the blob's
// fingerprint is merely
//
//    blob_fingerprint = xor_{0 <= l < L} site_state_fingerprint(l, s_l).
//
// If a state at site `l` is ambiguous between two states `a` and `b`, then the blob
// fingerprint includes the XOR of the site-state fingerprints for both `a` and `b`.
//
// If a blob's fingerprint is known, and the state of site `l` is `a`, then the
// fingerprint of a related blob with state `b` there is simply:
//
//   blob_fingerprint = old_blob_fingerprint ^ site_state_fingerprint(l, a) ^ site_state_fingerprint(l, b)
//
// Here are some worked examples (ssf = site_state_fingerprint):
//
//  Blob sequence         Fingerprint
//  -------------         -----------
//  {A,C}{A,G}  A         (ssf(0, A) ^ ssf(0, C)) ^ (ssf(1, A) ^ ssf(1, G)) ^ ssf(2, A)
//    C  {A,G}  A                      ssf(0, C)  ^ (ssf(1, A) ^ ssf(1, G)) ^ ssf(2, A) = above ^ ssf(0, A)
//    C    G    A                      ssf(0, C)  ^              ssf(1, G)  ^ ssf(2, A) = above ^ ssf(1, A)
//
// The exact sequence of a blob in a specific Phylo_tree is recreated from the blob tree
// itself: each blob records how its sequence differs from its parent blob's (see
// `Blob_site_delta` below), so a blob's sequence is the tree's reference sequence plus
// the deltas along the path from the root blob to it.  The same deltas are used to update
// the fingerprints as we go down the blob tree.  Each blob also records a "seed" node or
// branch, the earliest point on the Phylo_tree that lies in the blob.
//
using Blob_fingerprint = uint64_t;
inline constexpr auto k_empty_blob_fingerprint = Blob_fingerprint{0};
inline auto site_state_fingerprint(Site_index l, Real_seq_letter a) -> Blob_fingerprint {
  DCHECK_GE(l, 0);
  static_assert(std::is_same_v<Blob_fingerprint, uint64_t>, "If Blob_fingerprint is not uint64_t, site_state_fingerprint needs to change significantly");
  auto x = uint64_t{(static_cast<uint64_t>(l) << 8) | index_of(a)};
  auto z = x + 0x9e3779b97f4a7c15;
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
  z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
  return z ^ (z >> 31);
}

inline auto site_state_fingerprint(Site_index l, Seq_letter a) -> Blob_fingerprint {
  auto result = k_empty_blob_fingerprint;
  if (a & Seq_letters::A) { result ^= site_state_fingerprint(l, Real_seq_letter::A); }
  if (a & Seq_letters::C) { result ^= site_state_fingerprint(l, Real_seq_letter::C); }
  if (a & Seq_letters::G) { result ^= site_state_fingerprint(l, Real_seq_letter::G); }
  if (a & Seq_letters::T) { result ^= site_state_fingerprint(l, Real_seq_letter::T); }
  return result;
}

// Map the end-points of the branch in a pruned site tree for site `l` that contains a
// specific branch in the full Phylo_tree (see Supplementary Figure 11 in the Delphy
// paper's SI).  Most often, these will just be the endpoints of the branch in the full
// tree, but missations and going around the root mean that in the general case, it is a
// path in the full tree ending at two specific nodes (possibly `k_no_node` for an
// unbounded path that goes above the root, see edge case (c) above).  Note that a branch
// in the full tree is only part of a branch in the pruned site tree for site `l` if that
// branch is not downstream of a missation for site `l`.  In such cases, this function may
// return `{k_no_node, k_no_node}`, but only if the absence of `full_branch` in the pruned
// site tree can be determined by local exploration of the tree (optionally, can run an
// expensive check to conclusively decide whether the full branch is in the pruned site
// tree, with a POSSIBLE BUG warning issued if it is not).
auto map_site_tree_branch(
    const Phylo_tree& tree,
    Branch_index full_branch,
    Site_index l)
    -> std::pair<Node_index, Node_index>;

// Map the extent of a mutation, including its state at both endpoints and the mrca of the
// endpoints (one of start_node, end_node or the tree root).  The result is std::nullopt
// if there is no mutation on the pruned site-`l` tree branch containing `branch`.  As
// highlighted in edge case (b), multiple mutations may have no net effect; we still
// return a "neutral" extent in that case, with `state_at_start == state_at_end`.  It is
// up to the caller to decide what to do with those (see `map_tree_mutation_extents`).
// For "half" mutation extents (edge case (c) above), one of the endpoints will be
// k_no_node, and the corresponding state is the state at the root node.  Single mutation
// extents aren't directional; we enforce that start_node < end_node.
struct Mutation_extent {
  Site_index site;
  Node_index start_node;
  Real_seq_letter state_at_start;
  Node_index end_node;
  Real_seq_letter state_at_end;
  Node_index mrca_node;
  
  auto operator<=>(const Mutation_extent& that) const = default;
};
inline auto operator<<(std::ostream& os, const Mutation_extent& e) -> std::ostream& {
  return os << absl::StreamFormat("%c%d%c on %d<-%d->%d", to_char(e.state_at_start), e.site, to_char(e.state_at_end), e.start_node, e.mrca_node, e.end_node);
}
auto map_mutation_extent(
    const Phylo_tree& tree,
    Branch_index branch,
    Site_index l)
    -> std::optional<Mutation_extent>;

// The sequence of the root blob can differ from the root's actual sequence (`ref_sequence
// + root.mutations`) at sites whose pruned site tree branch goes through (or ends at) the
// root, because that branch is unconstrained: all that matters is the state at its
// end-points.  Besides genuinely ambiguous mutation extents through the root (which show
// up as active mutation extents of the root blob), there are two cases where the root
// blob has a *certain* state different from the root's actual state: "half" mutation
// extents (edge case (c) above), and "neutral" mutation extents through the root (edge
// case (d) above).  We don't treat those as mutation extents of the blob tree; we only
// record what the root blob's state is at those sites:
struct Root_state_override {
  Site_index site;
  Real_seq_letter state;

  auto operator<=>(const Root_state_override& that) const = default;
};
inline auto operator<<(std::ostream& os, const Root_state_override& o) -> std::ostream& {
  return os << absl::StreamFormat("%c@%d", to_char(o.state), o.site);
}

// `map_tree_mutation_extents` does what it says on the tin: walk through all the
// mutations in a Phylo_tree and build a list of (unique) mutation extents found.  "Half"
// extents and "neutral" extents through the root are not included in `extents`; they are
// only recorded as `root_state_overrides` (also unique, and sorted by site).
using Mutation_extent_index = int;
struct Mapped_mutation_extents {
  std::vector<Mutation_extent> extents;
  std::vector<Root_state_override> root_state_overrides;
};
auto map_tree_mutation_extents(const Phylo_tree& tree) -> Mapped_mutation_extents;

// `annotate_branches_with_extents` makes a map from each branch to (a) the mutation
// extents that begin at the start of the branch and (b) the mutation extents that end at
// the end of the branch.
struct Branch_extent_info {
  std::vector<Mutation_extent_index> mutation_extents_at_branch_start;
  std::vector<Mutation_extent_index> mutation_extents_at_branch_end;
  
  auto operator<=>(const Branch_extent_info& that) const = default;
};
inline auto operator<<(std::ostream& os, const Branch_extent_info& e) -> std::ostream& {
  return os << absl::StreamFormat("<(%s)---(%s)>",
                                  absl::StrJoin(e.mutation_extents_at_branch_start, ",", absl::StreamFormatter()),
                                  absl::StrJoin(e.mutation_extents_at_branch_end, ",", absl::StreamFormatter()));
}
auto annotate_branches_with_extents(
    const Phylo_tree& tree,
    const std::vector<Mutation_extent>& extents)
    -> Branch_vector<Branch_extent_info>;

// Partition tree into blobs.  Attach a blob number to each node and each branch.  Attach
// to each blob a list of "active mutation extents".  Keep track of which mutation extents
// are active at the root.  Build the root blob's sequence & fingerprint from the root's
// actual sequence, the root state overrides above, the extents active at the root and the
// sites missing in the whole tree.  Then derive the sequence & fingerprint of every other
// blob during the traversal.  Also, record which is the root blob, and for each blob,
// what its parent blob is.  Finally, each blob also records which of its constituent
// nodes & branches is the "start" of the blob (the "seed" node/branch).  Note: At this
// point, we differentiate blobs that have the same fingerprint but recur across the tree.
//
// Each blob also records how its sequence differs from its parent blob's (or from the
// tree's reference sequence, for the root blob) as a list of `Blob_site_delta`s, sorted
// by site.  The state at a site is a `Seq_letter`: a single base if the state is certain,
// or two bases if it is ambiguous (the site is covered by an active mutation extent).  Sites
// missing in the entire tree aren't part of any delta; they are the same for every blob, so
// they are recorded once, in `Blob_tree::missing_sites`.  Each delta records both the state
// at the parent blob (`from_state`) and the state at this blob (`to_state`).  These deltas are
// computed top-down during blob construction (in the same breath as the fingerprints), so
// a blob's sequence is just the reference sequence plus the deltas along the path from
// the root blob to it.  In particular, mutations along the way within a mutation extent,
// or in any path that has no net effect, don't affect the sequence of any blob.
//
// Because deltas record old states, a Blob_tree doesn't need the Phylo_tree or its
// reference sequence to answer questions about blob sequences: the first delta touching a
// site along any path from the root blob has the reference sequence's state at that site
// as its `from_state`.
using Blob_index = Node_index;
inline constexpr auto k_no_blob = k_no_node;
struct Blob_site_delta {
  Site_index site;
  Seq_letter from_state;  // State at the parent blob (for the root blob: the reference sequence's state)
  Seq_letter to_state;    // State at this blob

  auto operator==(const Blob_site_delta& that) const -> bool = default;
};
struct Blob_node : public Nary_node {
  Blob_fingerprint fingerprint;
  std::vector<Blob_site_delta> sequence_delta;  // vs. parent blob (or reference sequence for the root blob)
  Node_index seed_node;      // Exactly one of `seed_node` and `seed_branch` is `k_no_node`;
  Branch_index seed_branch;  // the blob starts either at a node or at a branch, but not both
  std::vector<Node_index> nodes;
  std::vector<Branch_index> branches;
  std::vector<Mutation_extent_index> active_mutation_extents;

  auto operator==(const Blob_node& that) const -> bool = default;
};
inline auto operator<<(std::ostream& os, const Blob_node& node) -> std::ostream& {
  return os << absl::StreamFormat(
      "Blob_node{fingerprint=\"0x%016x\", nodes=[%s], branches=[%s], active_mutation_extents=[%s]; parent=%d, children=[%s]}",
      node.fingerprint,
      absl::StrJoin(node.nodes, ", ", absl::StreamFormatter()),
      absl::StrJoin(node.branches, ", ", absl::StreamFormatter()),
      absl::StrJoin(node.active_mutation_extents, ", ", absl::StreamFormatter()),
      node.parent,
      absl::StrJoin(node.children, ", "));
}

struct Blob_tree : public Tree<Blob_node> {
  using Tree::Tree;

  std::vector<Mutation_extent> extents;
  Branch_vector<Branch_extent_info> branch_extents;
  std::vector<Root_state_override> root_state_overrides;  // See above; already baked into root blob's sequence
  Interval_set<> missing_sites;  // Sites missing in the whole tree (the same for every blob)

  Node_vector<Blob_index> node_2_blob;
  Branch_vector<Blob_index> branch_2_blob;
};
auto derive_blob_tree(const Phylo_tree& tree) -> Blob_tree;

// Query for a given blob what its exact "sequence + ambiguity" is.  Build this as a list
// of set of definite missing sites (only for sites missing on the whole tree), definite
// mutations from the reference sequence + list of ambiguous sites.  Computed in a single
// walk from `blob` up to the root blob, combining the `sequence_delta`s along the way.
// Needs neither the Phylo_tree nor the Blob_tree's mutation extents, so the result is
// self-contained.
struct Blob_sequence_mutation {
  Site_index site;
  Real_seq_letter from_state;  // State at the start of the mutation (here, in the reference sequence)
  Real_seq_letter to_state;    // State at the end of the mutation (here, in the blob)

  auto operator==(const Blob_sequence_mutation& that) const -> bool = default;
};
inline auto operator<<(std::ostream& os, const Blob_sequence_mutation& m) -> std::ostream& {
  return os << absl::StreamFormat(
      "Blob_sequence_mutation{%c%d%c}",
      to_char(m.from_state),
      m.site,
      to_char(m.to_state));
}
struct Blob_sequence_ambiguity {
  Site_index site;
  Real_seq_letter state1;  // The two possible states, with state1 < state2 (in A < C < G < T order)
  Real_seq_letter state2;

  auto operator==(const Blob_sequence_ambiguity& that) const -> bool = default;
};
inline auto operator<<(std::ostream& os, const Blob_sequence_ambiguity& a) -> std::ostream& {
  return os << absl::StreamFormat(
      "Blob_sequence_ambiguity{%c|%c @ %d}",
      to_char(a.state1),
      to_char(a.state2),
      a.site);
}
struct Blob_sequence {
  Interval_set<> missing_sites;
  std::vector<Blob_sequence_mutation> mutations;
  std::vector<Blob_sequence_ambiguity> ambiguities;
};
auto query_blob_sequence(
    const Blob_tree& blob_tree,
    Blob_index blob)
    -> Blob_sequence;

// Although we record all the blob fingerprints directly, as a sanity check, we can also
// calculate a blob fingerprint from scratch (and independently of the sequence deltas)
// from a Blob_sequence and the reference sequence it refers to
auto calc_blob_fingerprint(
    const Real_sequence& ref_sequence,
    const Blob_sequence& blob_sequence)
    -> Blob_fingerprint;

// TODO: Calculate all of this for all the Phylo_trees in a list.  Then pick out the stable blobs (> threshold).
// Create collapsed blob tree for an arbitrary one of those Phylo_trees.

// TODO: Implement root-folding.  It'll be something like: find and mark the smallest set of blob nodes that cover
// the blob tree roots ~99% of the time, then in each collapsed blob tree, find the MST of those root-candidate blobs
// and collapse that MST into a single "rooty" blob.

}  // namespace delphy

#endif // DELPHY_POPLAR_H_
