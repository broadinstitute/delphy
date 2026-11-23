# Summary blob trees: aggregate many blob trees into one

## Summary

A single `Blob_tree` describes one posterior tree.  To say something about the posterior as a
whole, we want to know which blobs show up in (almost) every posterior tree.  This plan adds, to
`core/poplar.{h,cpp}`:

* a **support aggregator**, `Blob_support_aggregator`, that accumulates the **posterior support**
  of every blob fingerprint seen so far (the fraction of the blob trees added to it that contain
  a blob with that fingerprint).  Blob trees are added to it one at a time as the MCMC produces
  new posterior samples;
* a `Summary_blob_tree` data structure, and a function `derive_summary_blob_tree()` that builds
  one from an aggregator and a single `Blob_tree` that donates the topology.  The summary tree
  contains only the blobs of that `Blob_tree` whose fingerprints have posterior support above a
  threshold (99% by default), with all the other blobs "compacted away";
* tests; and
* a WASM interface for all of it, with the matching bindings in `delphy-web`'s `delphy_api.ts`.

Splitting accumulation from summarization means that deriving a summary tree costs the same
whether the aggregator holds 10 blob trees or 10,000: it only depends on the size of the one
topology-donating `Blob_tree`.

This is the infrastructure only.  Showing summary blob trees in `delphy_ui` is the next work
unit, and is out of scope here.

## Background

* **Blob fingerprints** (`Blob_node::fingerprint`) are 64-bit hashes of a blob's sequence (with
  ambiguities), built from absolute per-site states, so equal sequences have equal fingerprints
  in different blob trees.  Aggregation only makes sense for blob trees derived from samples of
  the posterior tree distribution of a single dataset, so all of them have the same number of
  sites (nothing checks this, and `Blob_tree` does not grow any information to allow it).
* **Distinct blobs sharing a fingerprint are rare.**  Within one blob tree it takes something
  like a mutation on one branch that is reversed on the very next branch (the blobs before the
  first mutation and after the reversal have the same sequence).  The summary tree contains at
  most one blob per fingerprint: the earliest one.
* **`Mcc_tree` is the precedent** for aggregating many trees: `derive_mcc_tree()` takes a vector
  of `Phylo_tree*` base trees, and records "posterior support" per node
  (`Mcc_node::posterior_support()`, the fraction of base trees with an exact match).  Its WASM
  interface passes an array of tree pointers plus a count, and `delphy_api.ts` fills a
  `Uint32Array` with the `phyloTreePtr_` values.  Unlike `Mcc_tree`, the support aggregation
  here is incremental.
* **Blob trees are n-ary** (`Blob_node` derives from `Nary_node`), and the generic traversals in
  `core/tree.h` work on them.  A `Blob_tree` is self-contained for sequence queries (it needs
  neither the `Phylo_tree` nor its reference sequence), but its node, branch and seed indices only
  mean something relative to the `Phylo_tree` it was derived from.
* **The root blob oscillates** between posterior trees (observed in `delphy_ui`), so it is quite
  likely to have low posterior support.  This work unit deliberately sidesteps that, see "Out of
  scope" below.

## Design

### 1. Accumulating posterior support (`core/poplar.h`)

```cpp
// Accumulates, over blob trees that are samples from the posterior of one dataset, how many of them contain a blob
// with each fingerprint
struct Blob_support_aggregator {
  int num_blob_trees = 0;
  absl::flat_hash_map<Blob_fingerprint, int> blob_tree_counts;  // Number of blob trees with a blob with each fingerprint

  // Add one more blob tree.  A fingerprint that occurs in several blobs of the same blob tree still only counts once
  // for that blob tree.
  auto add_blob_tree(const Blob_tree& blob_tree) -> void;

  // Fraction of the blob trees added so far that contain a blob with this fingerprint (0 for fingerprints
  // never seen).  Requires at least 1 blob tree.
  auto posterior_support_of(Blob_fingerprint fingerprint) const -> double;
};
```

`add_blob_tree()` gathers the fingerprints of all blobs in the tree, deduplicates them (sort +
unique on a scratch vector), increments `num_blob_trees`, and increments the count of each
distinct fingerprint once.  Cost is linear in the size of the added blob tree.  Blob trees are
added one at a time, which is all that is needed: callers add each new posterior sample as it
arrives.

The aggregator stores the support of *all* the fingerprints it has seen, not only the ones that
make it into a summary tree, so that a UI can use it to highlight (or dim) low-support blobs in
individual blob trees.

### 2. The summary tree (`core/poplar.h`)

```cpp
// Index of a blob in a Summary_blob_tree (not to be confused with a Blob_index into a Blob_tree)
using Summary_blob_index = Node_index;

struct Summary_blob_node : public Nary_node {
  Blob_index blob;  // Corresponding blob in the summary tree's topology tree

  auto operator==(const Summary_blob_node& that) const -> bool = default;
};

struct Summary_blob_tree : public Tree<Summary_blob_node> {
  using Tree::Tree;

  const Blob_tree* topology_tree = nullptr;                      // Donates the topology
  const Blob_support_aggregator* support_aggregator = nullptr;   // Holds the supports of all fingerprints
  double posterior_support_threshold = 0.0;                      // Threshold used to decide which blobs to keep

  // Fingerprint and posterior support of a summary blob (via the blob in the topology tree)
  auto fingerprint_of(Summary_blob_index summary_blob) const -> Blob_fingerprint;
  auto posterior_support_of(Summary_blob_index summary_blob) const -> double;
};
```

Design choices:

* **Not self-contained.**  A `Summary_blob_tree` holds pointers to the `Blob_tree` that gives it
  its topology and to the aggregator whose supports it was built from.  The caller is assumed to
  hold on to both indefinitely (just as `Mcc_tree` holds raw pointers to its base trees).  This
  keeps the summary tree small: a summary blob is just an index into the topology tree, and
  everything else (fingerprint, sequence, member nodes and branches, seeds, ...) is one lookup
  away in that blob tree.
* **Supports are not copied into the nodes.**  `Summary_blob_tree::posterior_support_of()` asks
  the aggregator.  A new summary tree is expected to be derived whenever the aggregator is
  updated: querying an old summary tree against an updated aggregator doesn't make sense, and
  nothing goes out of its way to defend against it.
* `Summary_blob_index` is a distinct alias from `Blob_index`, to keep the two index spaces
  apart in signatures.  (`Blob_vector` / `Blob_map` remain for per-blob data in a `Blob_tree`.)

### 3. Building the summary tree

```cpp
inline constexpr auto k_default_summary_blob_support_threshold = 0.99;

// `topology_tree` donates the topology: a blob B in the summary tree has as children the earliest blobs
// below B in `topology_tree` that are kept.  Requires `support_aggregator.num_blob_trees > 0` and a threshold in [0, 1].
auto derive_summary_blob_tree(
    const Blob_support_aggregator& support_aggregator,
    const Blob_tree& topology_tree,
    double posterior_support_threshold = k_default_summary_blob_support_threshold)
    -> Summary_blob_tree;
```

Typically the topology donor is the most recent blob tree added to the aggregator, and the
summary tree is re-derived each time the aggregator is updated, but any blob tree can be used.

**Which blobs are kept.**  Visit the blobs of `topology_tree` in pre-order (parents before
children, children in their stored order).  A blob is kept if:

* it is the root blob of `topology_tree` (always, regardless of its support, see "Out of scope"), or
* the posterior support of its fingerprint is at least the threshold, **and** no blob with the
  same fingerprint has been kept earlier in the traversal.

So at most one blob per fingerprint survives: the earliest one in the traversal.  Among blobs on
the same path from the root, that is the one closest to the root, which is the earliest in time;
the other blobs with that fingerprint, which are downstream of it, are compacted away like any
low-support blob.  Equal-fingerprint blobs that are *not* downstream of each other (rarer
still) are resolved in traversal order.

**Building the tree without a parallel vector.**  Use `traversal(topology_tree)`, which, like in
`derive_blob_tree()`, produces an event when entering a blob (`children_so_far == 0`) and one
when leaving it (`children_so_far == number of children`; both in one go for a leaf).  Keep one
piece of state: `cur_summary_blob`, the nearest kept ancestor-or-self of the blob being visited
(`k_no_node` before the root):

```cpp
for (const auto& [blob, children_so_far] : traversal(topology_tree)) {
  if (children_so_far == 0) {
    // Entering `blob`
    if (is_kept(blob)) {
      // New summary node, whose parent is `cur_summary_blob`; append it to that node's children (or make it the
      // root); record `blob` in it; remember its fingerprint as taken
      cur_summary_blob = new_node;
    }
  }
  if (children_so_far == std::ssize(topology_tree.at(blob).children)) {
    // Leaving `blob`: if it was kept, `cur_summary_blob` is its summary node, so go back up to that node's parent
    if (cur_summary_blob != k_no_node && summary_blob_tree.at(cur_summary_blob).blob == blob) {
      cur_summary_blob = summary_blob_tree.at(cur_summary_blob).parent;
    }
  }
}
```

(Whether `blob` is kept is decided from the aggregator and a `flat_hash_set` of fingerprints
already taken.)  The test on leaving works because every kept summary node records a distinct
blob.  (The implementation happens to number summary nodes in pre-order, but that is incidental:
the indices of summary nodes are unspecified, and neither callers nor tests may rely on them.)

The net effect is that a summary blob's children are exactly the earliest kept blobs below the
corresponding blob in `topology_tree`, in the same left-to-right order.  Cost is linear in the
size of `topology_tree`; the result is deterministic (the hash set is only ever used for
membership tests).

The threshold comparison is simply `posterior_support >= threshold` everywhere.

### 4. WASM interface (`tools/delphy_wasm.cpp`)

A new `// Blob support aggregator and summary blob tree` section after the blob tree section,
following the shape of the existing blob tree exports.  The work is cheap (linear in the size of
one blob tree for both adding a tree and deriving a summary), so everything is synchronous, with
no thread pool or callbacks.

| Export                                                                | Returns                          |
|-----------------------------------------------------------------------|----------------------------------|
| `delphy_blob_support_aggregator_new(ctx)`                             | `Blob_support_aggregator*`       |
| `delphy_blob_support_aggregator_delete(ctx, aggregator)`              | `void`                           |
| `delphy_blob_support_aggregator_add_blob_tree(ctx, aggregator, blob_tree)` | `void`                      |
| `delphy_blob_support_aggregator_get_posterior_support_of(ctx, aggregator, uint64_t fingerprint)` | `double`  |
| `delphy_derive_summary_blob_tree(ctx, aggregator, topology_blob_tree, double threshold)` | `Summary_blob_tree*` |
| `delphy_summary_blob_tree_delete(ctx, tree)`                          | `void`                           |
| `delphy_summary_blob_tree_get_size(ctx, tree)`                        | `Summary_blob_index`             |
| `delphy_summary_blob_tree_get_root(ctx, tree)`                        | `Summary_blob_index`             |
| `delphy_summary_blob_tree_get_parent_of(ctx, tree, blob)`             | `Summary_blob_index`             |
| `delphy_summary_blob_tree_get_num_children_of(ctx, tree, blob)`       | `int`                            |
| `delphy_summary_blob_tree_get_nth_child_of(ctx, tree, blob, n)`       | `Summary_blob_index`             |
| `delphy_summary_blob_tree_get_blob_of(ctx, tree, blob)`               | `Blob_index` (in the topology tree) |
| `delphy_summary_blob_tree_get_fingerprint_of(ctx, tree, blob)`        | `uint64_t` (a BigInt in JS, as for blob trees) |
| `delphy_summary_blob_tree_get_posterior_support_of(ctx, tree, blob)`  | `double`                         |

This lets JS (a) accumulate supports gradually as posterior samples arrive, (b) query the
support of any fingerprint, not just those in a summary tree, and (c) build a summary tree from
the supports accumulated so far and an arbitrary `Blob_tree`.  There are deliberately no getters
for the number of blob trees aggregated or for the threshold used.

Conventions as for the existing blob tree exports: `const` references for getters, `.at()` for
bounds-checked access, a `static_assert` that `Blob_fingerprint` is `uint64_t` next to the
fingerprint getters, and a section comment saying that a summary blob tree holds pointers to the
blob tree and the aggregator it was derived from, so both must outlive it, and that a summary
tree should be re-derived whenever the aggregator is updated.

### 5. `delphy-web` bindings (`src/ts/pythia/delphy_api.ts`, and `pythia.ts`)

Per the project rule, every change in `delphy_wasm.cpp` comes with the matching change in
`delphy_api.ts`:

* `BlobSupportAggregatorPtr` and `SummaryBlobTreePtr` types (next to `BlobTreePtr`), and the raw
  function type declarations and bindings in `delphyCoreRaw` for each export above;
* a `BlobSupportAggregator` class with `delete()`, `addBlobTree(blobTree: BlobTree)` and
  `getPosteriorSupportOf(fingerprint: bigint): number`, created with
  `Delphy.newBlobSupportAggregator()`;
* `Delphy.deriveSummaryBlobTree(aggregator: BlobSupportAggregator, topologyBlobTree: BlobTree,
  threshold: number = 0.99): SummaryBlobTree`;
* a `SummaryBlobTree` class mirroring `BlobTree`'s getters (`delete`, `getSize`, `getRootIndex`,
  `getParentIndexOf`, `getNumChildrenOf`, `getNthChildIndexOf`, `getChildrenOf`, plus
  `getBlobOf`, `getFingerprintOf` as a `bigint` with `BigInt.asUintN(64, ...)`, and
  `getPosteriorSupportOf`).  The name `SummaryBlobTree` avoids a clash with the existing
  `SummaryTree` interface, which is specific to MCC-like trees (base trees and corresponding
  nodes).  It keeps references to the `BlobSupportAggregator` and the topology `BlobTree` it
  was derived from, since the C++ side only holds raw pointers to them;
* thin `newBlobSupportAggregator()` and `deriveSummaryBlobTree()` pass-throughs in `Pythia`,
  like the one for `deriveBlobTree`.

The bundled `delphy.js` / `delphy.wasm` in `delphy-web` must be rebuilt and copied over when the
`delphy-web` side is committed, as for the earlier blob tree work.  No UI changes in either
`delphy-web` or `delphy_ui` in this work unit.

### 6. Housekeeping in `poplar.h`

The TODO at the end of `poplar.h` ("Calculate all of this for all the Phylo_trees in a list.  Then
pick out the stable blobs (> threshold).  Create collapsed blob tree for an arbitrary one of those
Phylo_trees.") is what this work implements, so it is removed.  The root-folding TODO stays.  The
header comment above the new declarations explains posterior support, the threshold, the
compaction and the earliest-blob-per-fingerprint rule.  `poplar.h` also needs
`#include "absl/container/flat_hash_map.h"` for the aggregator's counts (`poplar.cpp` already
includes `absl/container/flat_hash_set.h`, which the traversal needs for the fingerprints already
taken), and per the include conventions each header includes what it uses directly.

## Tests (`tests/poplar_tests.cpp`)

All blob trees in the tests are derived from `Phylo_tree`s with `derive_blob_tree()`; there are no
hand-built mock `Blob_tree`s, so the tests do not need maintenance whenever the `Blob_tree` data
structure changes.

**Fixtures.**  Reuse the existing hand-built `Phylo_tree` fixtures, and create variants of them
by copying the fixture's tree and editing a few mutations (re-checking
`assert_phylo_tree_integrity()`), so that some blobs keep their fingerprints and others change.
Expected supports are read off the derived blob trees themselves (e.g., the fingerprint of the
blob containing a given `Phylo_tree` node, via `node_2_blob`), never hard-coded hashes.  Finding
variants that exercise each scenario below is part of writing the tests, and every scenario
first asserts its preconditions on the derived blob trees (for example that two fingerprints
really coincide, or that an interior blob's fingerprint differs between two trees).

A dedicated **reversal fixture** gives blob trees with repeated fingerprints: a ladder-like
`Phylo_tree` with a mutation (say A0C) on one branch and its reversal (C0A) on the next, with
all tips informative so that each extent stays within a single branch.  The blobs before the
mutation and after the reversal then share a fingerprint, and so do the two ambiguous
single-branch blobs.

**Aggregator**

1. One blob tree: every fingerprint in it has count 1 and support 1; a fingerprint never seen
   has support 0.
2. Several identical blob trees: supports are 1.
3. A mix of variants: fingerprints shared by all variants have support 1, and others have the
   fraction of trees that contain them.
4. **Repeated fingerprints within a tree** (reversal fixture) count once for that tree.
5. **Incremental accumulation:** the counts after adding blob trees one at a time do not
   depend on the order in which they were added, and each addition changes only the counts of
   the fingerprints in the added tree (and `num_blob_trees`).

**Summary tree**

6. One blob tree (from a fixture without repeated fingerprints): with the default threshold, the
   summary has exactly the topology of the blob tree, every blob is kept, and `get_blob_of` maps
   the summary nodes back correctly.
7. Several identical blob trees: same.
8. **Compaction of low-support blobs:** a topology tree with blobs whose fingerprints have low
   support (because the other trees in the aggregator differ there) gives a summary without
   them, whose children are the earliest kept descendants.  This includes the interesting case
   of a low-support *interior* blob with several high-support descendants (they become siblings,
   in the same left-to-right order as in the topology tree), and a run of several consecutive
   low-support blobs.
9. **The root is kept even if its own support is below the threshold** (aggregate blob trees
   whose root blobs differ), and its low support is reported faithfully by
   `posterior_support_of()`.
10. **Threshold:** threshold 0 keeps every distinct fingerprint of the topology tree; threshold 1
    keeps only the fingerprints present in all the blob trees (plus the root); intermediate
    thresholds cut at the right counts, with `>=` semantics.
11. **Topology comes from the donor tree**: two topology donors with different shapes give
    different summaries from the same aggregator, while the supports reported are the same.
12. **Earliest blob per fingerprint** (reversal fixture): each fingerprint appears at most once in
    the summary tree, as the blob closest to the root, and the compacted-away duplicates'
    descendants attach to the nearest kept ancestor.
13. **Summary tree invariants** for all of the above: `assert_tree_integrity()` holds, the root
    is the topology tree's root blob, every node's blob is a distinct blob of the topology
    tree, the fingerprints of the nodes are distinct, every node's parent is its nearest kept
    ancestor in the topology tree, and `Summary_blob_tree::fingerprint_of()` and
    `posterior_support_of()` agree with the topology tree and the aggregator.  The indices of
    summary nodes (e.g., that the root is 0, or that they follow pre-order) are deliberately
    not tested, so tests list summary blobs in pre-order of the summary tree.

An empty aggregator is a `CHECK` failure; that is not worth a death test (there are death tests
in the suite, but this does not need one).

## Files changed

| File                                  | Action | What |
|---------------------------------------|--------|------|
| `core/poplar.h`                       | Modify | Add `Blob_support_aggregator`, `Summary_blob_index`, `Summary_blob_node`, `Summary_blob_tree` (with its `fingerprint_of()` / `posterior_support_of()` member functions), `k_default_summary_blob_support_threshold` and `derive_summary_blob_tree()` with documentation; remove the TODO that this work implements |
| `core/poplar.cpp`                     | Modify | Implement them |
| `tests/poplar_tests.cpp`              | Modify | Reversal fixture, `Phylo_tree` variants and the tests above |
| `tools/delphy_wasm.cpp`               | Modify | New aggregator and summary blob tree section with the exports above |
| `delphy-web: src/ts/pythia/delphy_api.ts` | Modify | Raw declarations, bindings, `BlobSupportAggregator` and `SummaryBlobTree` classes, `Delphy.newBlobSupportAggregator()`, `Delphy.deriveSummaryBlobTree()` |
| `delphy-web: src/ts/pythia/pythia.ts` | Modify | Pass-throughs |
| `delphy-web: src/ts/delphy/delphy.{js,wasm}` | Rebuild | Bundled WASM matching the new exports |

No changes to the CMake files (`poplar.cpp` and `poplar_tests.cpp` are already built) or to
`delphy_ui` (next work unit).

## Verification

1. Build and run the full unit test suite (debug), and the Poplar tests in the release build.
2. Build the WASM target; confirm the new exports show up in the generated `delphy.js`.
3. Type-check `delphy-web` with `tsc --noEmit` against the new declarations.
4. Convince ourselves that the supports are right on something real, ahead of the `delphy_ui`
   work: for a sequence of blob trees derived from a short MCMC run, compare the aggregator's
   counts with a brute-force recount in a throwaway test or script (kept out of the tree if it
   is not worth keeping), and check that the summary tree at threshold 0 has exactly one blob per
   distinct fingerprint of the topology tree.

## Out of scope: rooting the summary tree

The root of the summary tree is simply the root blob of the topology tree, even if that blob has
low posterior support, so that the summary is a proper tree with a single root.  Given that the
root blob oscillates, this will often be a low-support blob.  Later work will refine this: find a
minimal set of "rooty" blobs that were in root position and together account for, say, more than
99% of the posterior support; then, in the topology tree, build the minimal spanning tree of
those rooty blobs and treat it as a single "pseudo-root" node of the summary tree, whose children
are exactly the closest high-support blobs just outside that spanning tree.  (This is the
root-folding TODO that stays in `poplar.h`.)

## Decisions made during review

* **The "earliest" blob for a repeated fingerprint is the first one in pre-order traversal of
  the topology tree.**  For blobs on the same root-to-leaf path (the case that actually arises
  with a reversal), that is exactly the earliest in time.  For equal-fingerprint blobs in
  different lineages, pre-order is not time order, but being strictly time-based would need the
  underlying `Phylo_tree` (the blob tree has no times) and the case is extremely rare.
* **Summary trees do not freeze supports.**  `posterior_support_of()` always asks the aggregator,
  and a summary tree is simply re-derived whenever the aggregator is updated; querying an old
  summary tree against an updated aggregator makes no sense, so nothing defends against it.
* **Blob trees are added to the aggregator one at a time**; there is no bulk interface over a span
  of `Blob_tree*`.

## Suggested order of work

1. `core/poplar.{h,cpp}`: the aggregator and its tests; then the summary tree, its derivation
   and its tests.  Build and run the Poplar tests as you go.
2. `tools/delphy_wasm.cpp` and the delphy-web bindings (`delphy_api.ts`, `pythia.ts`), together,
   since they must stay in sync.
3. Full verification (see above), then commit (and push, if asked).
