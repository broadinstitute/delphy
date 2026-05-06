# Better Tree Init, Round 1: Utree and Guide Tree

This is the first in a series of plans implementing a new tree initialization pipeline for
Delphy, aimed at producing better starting trees for deep phylogenies (many mutations per
branch).

## Background

Delphy performs Bayesian phylogenetic inference via MCMC on a rooted, timed, explicit
mutation-annotated tree (EMAT). To kickstart the MCMC, we must construct an initial
consistent EMAT from input sequences and dates.

### Input data

The input to tree construction is a **reference sequence** plus N **tip descriptions**
(`Tip_desc`, `core/phylo_tree.h:137`). Each tip description contains:
- A name.
- Date bounds (`t_min`, `t_max`).
- `seq_deltas`: a list of `Seq_delta` (`core/mutations.h:67`) records `{site, from, to}`
  describing how the tip's sequence differs from the reference at non-missing sites.
- `missations`: a `Missation_map<>` (`core/mutations.h:125`) tracking missing-data intervals
  and entry states. In all downstream uses of `Tip_desc`, only the `.intervals` field (an
  `Interval_set`) is used.

The output is a `Phylo_tree` — Delphy's internal representation of an EMAT.
Each branch carries an ordered list of mutations (site, from state, to state, time), and the
tree tracks missing data per branch via "missation maps." The `Phylo_tree` also stores a
reference sequence, with the mutations above the root node listing the deltas between the reference and root sequences.

### The problem

The current initializer (`build_usher_like_tree`, `core/phylo_tree.cpp:729`) builds the EMAT
incrementally by adding tips one at a time. It has three problems for deep trees:

1. **Premature commitment to mutation distribution.** Mutations are placed on branches with
   specific orderings and times as tips are added. If mutation order is imputed incorrectly
   early on, downstream tips introduce spurious reversals instead of revising the original
   order.

2. **No ability to revisit the root.** The tree is built rooted from the start: the first
   two tips define the initial root, and subsequent tips either attach within the tree or
   above the root. If the initial root choice is poor (which is likely when the first two
   tips happen to be unrepresentative), every subsequent attachment decision is biased by
   this choice. The algorithm cannot backtrack to select a better root.

3. **O(N^2) construction time.** Every new tip is evaluated against *every* possible
   attachment point in the tree, with no pruning. For large trees, this dominates runtime.

### The new approach

Build an **unrooted mutation-annotated tree** first — where edges carry **site deltas**
(which sites change state, without ordering or timing) — then refine the topology, root it
using temporal signal, and convert to the full EMAT. Starting with an unrooted tree decouples
topology from timing, defers rooting, and simplifies refinement moves.

The approach mirrors a common phylogenetics workflow: UShER builds a rough MAT
(mutation-annotated tree) via greedy parsimony, matOptimize refines its topology via local
rearrangements, and TempEst roots and time-calibrates it. Our pipeline follows the same three
stages — greedy construction, topology refinement, rooting — implemented as a single
integrated system within Delphy. One simplification: UShER and matOptimize maintain global
parsimony (propagating each change through the full tree), but we use purely local operations
— each inserted tip only introduces new mutations on the branch connecting it to the existing
tree, and refinement moves only update the affected edges. We don't expect this to matter much
for our use cases, but it is a difference.

This work also serves as a testbed for new data structure ideas — the arc-based tree
representation and the focus-node mechanism described below — that may prove useful more
broadly in Delphy's MCMC infrastructure.

## Roadmap

The overall strategy is to get a basic end-to-end pipeline working as early as possible,
then improve each stage incrementally. The unrooted tree (`Utree`) is the central data
structure throughout: tips are attached to it by finding the edge that minimizes the number
of new site deltas (a greedy parsimony criterion), and the topology is later refined by
detaching and reattaching nodes.

1. **Round 1 (this plan):** Build the `Utree` data structure and use it to construct a
   rough "guide tree" by inserting tips one at a time in input order. Each tip is attached
   approximately at the edge where it introduces the fewest new site deltas, found via a
   priority-queue-driven branch-and-bound search. The guide tree captures approximate
   relationships — tips nearby in the guide tree are likely genuinely similar — but its
   topology may be poor because the input order is arbitrary.

2. **Round 2:** Implement midpoint rooting on the `Utree` (find the tree's diameter and root
   at its midpoint) and conversion from the rooted `Utree` to a `Phylo_tree` (assign times
   from root-to-tip mutation counts assuming a rough mutation rate, distribute mutations
   along branches, build missation maps). This yields a basic working pipeline: guide tree
   + midpoint root + conversion = usable initial `Phylo_tree`.

3. **Round 3:** Improve the incremental construction by using the guide tree to determine a
   better insertion order. Traverse the guide tree in priority-queue order (shortest edges
   first) to produce an ordering where each tip is close to some already-inserted tip. Then
   rebuild the `Utree` from scratch using this order. This reuses all the Round 1
   infrastructure but should produce a substantially better topology.

4. **Round 4:** Replace midpoint rooting with root-to-tip regression. For each candidate
   root edge, orient the tree and regress root-to-tip mutation counts against tip dates
   using least-squares. Pick the candidate with the best fit. This gives a better
   root position and also estimates the mutation rate and root date for time assignment.

5. **Round 5:** Topology refinement by SPR moves on tips. Detach each tip and search for a
   better reattachment edge using the same branch-and-bound search as construction. Accept
   moves that reduce total parsimony score. This corrects mistakes from the greedy
   construction.

6. **Round 6:** Generalize refinement to SPR moves on subtrees. Detach an inner node
   (separating one subtree from the rest), search for a better reattachment edge using
   Fitch parsimony to jointly optimize the sequences at both the detached node and the
   subtree root. This handles cases where an entire clade is misplaced.

---

## Goal for Round 1

Build a new unrooted tree data structure (`Utree`) and use it to construct a rough "guide
tree" from input tip sequences via greedy incremental insertion with branch-and-bound
pruning. The guide tree's topology may be poor (since tips are added in arbitrary input
order), but it captures approximate relationships: tips nearby in the guide tree are likely
genuinely similar.

The guide tree is built by adding tips one at a time and attaching each new tip
approximately at the edge where it introduces the fewest new site deltas. A
priority-queue-driven search explores candidate edges in order of increasing cost, with
aggressive pruning to keep construction fast. The existing `build_usher_like_tree`
(`core/phylo_tree.cpp:729`) uses a similar greedy parsimony strategy, but evaluates every
possible attachment point without pruning, making it O(N^2). The new approach uses
approximate branch-and-bound to avoid exploring most of the tree.

## Data structures

### Site deltas (existing)

A site delta on an edge records "at site 5, the state changes from C to T along this edge."
This is a symmetric, directionless fact (equivalently, "from T to C" in the other direction).

A `Site_delta` (`core/site_deltas.h:10`) is a `{from, to}` pair of `Real_seq_letter` at a
single site. A `Site_deltas` object (`core/site_deltas.h:18`) is a
`Scratch_flat_hash_map<Site_index, Site_delta>` — a sparse map storing only the sites that
change state. Existing operations handle composition:

- `push_front_site_deltas(delta_z, deltas)`: given deltas from just-after-z to x, prepend
  the mutation at z, producing deltas from just-before-z to x. Cancels if the mutation
  reverses an existing delta.
- `pop_front_site_deltas(delta_z, deltas)`: inverse — given deltas from just-before-z to x,
  remove the mutation at z, producing deltas from just-after-z to x.
- `append_site_deltas(x_to_y, y_to_z)`: compose two delta paths end-to-end.

These are reused directly for traversal and edge merging.

A `Seq_delta` (`core/mutations.h:67`) is a `{site, from, to}` triple — a `Site_delta` at a
named site. It is used as the argument type for `push_front_site_deltas` and friends.

**Allocation note**: `Site_deltas` uses the thread-local arena allocator, which tree
initialization doesn't need. We introduce `Heap_site_deltas = absl::flat_hash_map<Site_index,
Site_delta>` in `core/site_deltas.h` and use it throughout the `Utree` code. The composition
operations are templatized on the map type and moved to `core/site_deltas.h` so they work
with either variant. Eventually, `Site_deltas` itself should be parameterized by allocator
(defaulting to the STL allocator), with `Scratch_site_deltas` as the arena-allocated alias —
mirroring the pattern used for other container types.

### `Utree`

A `Utree` is stored in a new file `core/utree.h` (and `core/utree.cpp`). It represents an
unrooted, bifurcating, mutation-annotated tree using **arcs** (directed edges). Each
undirected edge between nodes A and B is stored as two arcs, one pointing from A to B and
one from B to A; each arc is the other's **mate**. Each arc stores its target node and the
site deltas in that direction:

```cpp
using Arc_index = int;
inline constexpr auto k_no_arc = Arc_index{-1};

struct Uarc {
  union {
    Node_index target;            // The node this arc points TO (when allocated)
    Arc_index next_free_arc;      // Next free arc pair (when on free list)
  };
  Heap_site_deltas deltas;        // Site deltas FROM origin TO target; cleared when free
};
```

When allocated, `target` identifies the destination node; `next_free_arc` is used for
free-list bookkeeping (see "Pre-allocation" below). `deltas` lives outside the union and is
always a valid object.

The origin of an arc is its mate's target; the "Mate trick" below provides a mechanism for
quickly deducing an arc's mate. Both `origin(arc)` and `target(arc)` are provided as helper
methods on `Utree`; `target(arc)` is syntactic sugar for `arcs[arc].target`, preferred for
readability when used alongside `origin(arc)`.

Each node stores its outgoing arcs in a fixed-size array. The `arc_to_focus` field
supports the **focus node** mechanism (described in the next subsection), which gives the
unrooted tree a temporary orientation by designating one node as the focus:

```cpp
struct Unode {
  std::array<Arc_index, 3> arcs = {k_no_arc, k_no_arc, k_no_arc};
  Arc_index arc_to_focus = k_no_arc;  // Outgoing arc toward the focus; k_no_arc for the focus node
};
```

The tree structure:

```cpp
struct Utree {
  Real_sequence ref_sequence;         // Reference sequence
  Interval_set globally_missing_sites; // Sites missing across all tips so far

  std::vector<Uarc> arcs;             // All arcs (preallocated in mate pairs)
  std::vector<Unode> nodes;           // All nodes (preallocated)
  Arc_index arc_free_list_head;       // Head of free list of arc pairs
  int num_tips;                       // Tips are at indices [0, num_tips)
  int num_inner_nodes_so_far;         // Inner nodes allocated so far

  // Focus state
  Node_index focus;                   // Current focus node
  Heap_site_deltas deltas_ref_to_focus;  // Deltas from reference sequence to focus
};
```

Together with the input `tip_descs`, the `Utree` fully describes the dataset: the reference
sequence, each tip's deltas and missations, the tree topology, and the edge annotations. In
a future round, it may make sense to read input FASTA/MAPLE files directly into a `Utree`
representation instead of going through `Tip_desc`.

The next section describes the "focus state." A few additional notes:

- **Mate trick** (Knuth, TAOCP v4A pp. 21-22): Arcs are allocated in pairs at consecutive
even-odd indices: arc `i` and arc `i ^ 1` (XOR with 1) are **mates** — they represent the
same undirected edge in opposite directions. So, for example, the origin of arc `i` is
`arcs[i ^ 1].target`. Without this trick, finding the mate of an arc R to B would
require looking up node B, scanning its arc slots, looking up each target, and finding which
one has R as an outgoing arc. The mate trick replaces this with a single XOR instruction.
The rest of the design works without it — it is purely an addressing optimization, not
essential for correctness.

- **Node layout**: Tips are at indices `[0, num_tips)`. Inner nodes are allocated starting at
index `num_tips` as needed during construction. An unrooted bifurcating tree with N tips has
`2N - 2` nodes and `2N - 3` undirected edges = `2*(2N - 3)` arcs.

- **Degree**: Tips have 1 non-`k_no_arc` entry in `arcs[]`; inner nodes have 3. During
`split_edge`, a node temporarily has degree 2 (third slot empty) before the caller attaches
a tip or subtree.

- **Pre-allocation**: All node and arc slots are allocated upfront for the maximum tree size,
plus a few extra scratch slots to keep `split_edge` clean and simple. The **free list** is an
intrusive singly-linked list threaded through the arc pairs: the `next_free_arc` field of
the even-indexed arc in each free pair stores the index of the next free pair (or `k_no_arc`
for end-of-list).

### Focus node

The tree maintains a distinguished **focus node** F. The focus stores
`deltas_ref_to_focus`: the site deltas from the reference sequence to the focus node's
implicit sequence. This provides a quick way to compute the deltas from the reference
sequence to any node in the tree: walk from the node toward the focus (following
`arc_to_focus` at each node), composing arc deltas along the way, then compose with
`deltas_ref_to_focus`.

Each node's `arc_to_focus` is the outgoing arc that points toward the focus; the remaining
arcs point away from the focus. The focus node itself has `arc_to_focus = k_no_arc`. This
gives the unrooted tree a temporary orientation — like a rooted tree with F as the root —
without altering the topology.

It is also possible to focus on a particular undirected *edge* rather than a node, by having
two nodes' `arc_to_focus` point at each other. We don't use this in Round 1, but note it for
later rounds.

**Focus movement**: Moving the focus from F to a neighbor G is O(1):
- F gets `arc_to_focus` set to the arc pointing toward G.
- G gets `arc_to_focus = k_no_arc`.
- Apply the F-to-G arc's deltas to `deltas_ref_to_focus`.

Moving to a distant node X costs O(path_length): walk from X toward F, flipping each node's
`arc_to_focus` along the way, then walk back updating `deltas_ref_to_focus`.

A **focal arc** is an outgoing arc from the current focus node. Because the focus's sequence
is tracked in `deltas_ref_to_focus`, any quantity that depends on the focus's sequence can be
evaluated at a focal arc without moving the focus — this is a key property exploited during
the branch-and-bound search.

**Path finding**: To find the path (and compute deltas) from any node X to any node Y: move
the focus to Y, then follow parent arcs from X toward the focus (Y), composing deltas along
the way.

## Core operations

### `Utree::make_empty(int num_tips) -> Utree`

Static factory that creates a `Utree` pre-allocated for `num_tips` tips. Sizes the `nodes`
vector for `max(1, 2*num_tips - 2)` nodes (tips + inner nodes) and the `arcs` vector for
`max(1, 2*num_tips - 3 + 2)` arc pairs (one per edge, plus 2 scratch pairs so `split_edge`
can allocate new edges before freeing the old one). All arc pairs start on the free list.

### Helpers

- `alloc_arc_pair() -> Arc_index` / `free_arc_pair(arc)` — pop/push arc pairs from the free list.
- `add_arc(A, B) -> Arc_index` — allocate an edge between A and B, wire both nodes' arc slots, return the arc A→B.
- `find_arc(node, target) -> Arc_index` — look up the outgoing arc from `node` to `target`.
- `degree(node) -> int` — count of non-`k_no_arc` arc slots.

### `reset_focus(Node_index F)`

Set node F as the focus, recomputing all `arc_to_focus` fields from scratch. Unlike
`move_focus_to`, this does **not** update `deltas_ref_to_focus` — the caller is responsible
for setting it. Useful for initial setup or when `deltas_ref_to_focus` will be rebuilt anyway.
Implemented via `annotated_arc_euler_tour`: on each `entering` arc, set the target's
`arc_to_focus` to the mate of that arc.

### `split_edge(Arc_index arc_AB, Node_index M, site_delta_side)`

Split the edge identified by arc `arc_AB` (pointing from A to B) by inserting node M.
Replaces the undirected edge (A, B) with two edges: (A, M) and (M, B).

- Allocate two new edges (A-M and M-B). For each site delta in the old A-B edge, call
  `site_delta_side(Seq_delta{site, from, to}, A, B) -> Node_index`. The returned value
  must be A (delta goes on the A-M side) or B (delta goes on the M-B side); any other
  value CHECK-fails. Free the old edge.
- M gets two arc slots filled (toward A and toward B); the third slot remains `k_no_arc`
  for the caller to attach a tip.
- Sets M's `arc_to_focus` to the new arc pointing toward the focus side (whichever of A, B
  has a path to the focus). Updates A's and B's arc slots to point to the new arcs. If A or
  B had `arc_to_focus` set to one of the old (now-freed) arcs, updates it to the
  corresponding new arc.

Returns nothing. If a caller needs the new arcs, it can look them up from M's outgoing arcs.

### `move_focus_to(Node_index target, pre_arc_hop, post_arc_hop = noop)`

Move the focus from its current position to `target`. Since `arc_to_focus` links point
*toward* the focus, the path to `target` isn't directly available. Two passes solve this:
first walk from `target` toward the old focus, flipping `arc_to_focus` links along the way
so they point toward `target`; then walk back from the old focus toward `target`, applying
each arc's deltas to `deltas_ref_to_focus`. At each step of the forward pass, caller-supplied
hooks `pre_arc_hop(R)` and `post_arc_hop(R)` fire, allowing the caller to update
focus-associated state (e.g., `focus_to_X_deltas`) in tandem. Cost: O(path_length).

### `annotated_arc_euler_tour(const Utree& tree, Node_index source) -> cppcoro::generator<Annotated_arc>`

A generator that yields every arc hopped along a DFS of the tree, starting at `source`,
annotated with direction:

```cpp
enum class Arc_direction { entering, leaving };
struct Annotated_arc { Arc_index arc; Arc_direction direction; };
```

An Euler tour of a tree visits every edge exactly twice (once entering a subtree, once
leaving it). Each yielded arc carries the deltas in the direction of travel — entering arcs
advance deeper into the tree, leaving arcs undo the forward step. Implemented via a
stack-based DFS.

`arc_euler_tour` is a convenience wrapper that discards the direction annotation, yielding
only the `Arc_index`.

### `assert_utree_integrity(const Utree& tree, bool force = false)`

Verify the standalone structural integrity of a `Utree`, without reference to any external
input data. Modeled after `assert_phylo_tree_integrity` (`core/phylo_tree.h`): a free
function that CHECKs on any violation. Gated by `estd::is_debug_enabled || force` — a
no-op in release builds unless `force` is true.

Checks node degrees (tips have 1, inner nodes have 3), arc wiring (targets, mates, and delta
reversal are consistent), `arc_to_focus` tree (acyclic, rooted at focus), free list integrity,
`deltas_ref_to_focus` sanity (sites in range, `from` matches reference, no globally missing
sites), and path consistency (an Euler tour from the focus should compose deltas consistently
and round-trip back to `deltas_ref_to_focus`).

### `assert_utree_matches_tip_descs(const Utree& tree, const std::vector<Tip_desc>& tip_descs, bool force = false)`

Verify that a `Utree`'s tip sequences and globally missing set are consistent with the
input `Tip_desc`s. Assumes the tree already passes `assert_utree_integrity`. Gated by
`estd::is_debug_enabled || force`.

First validates the `tip_descs` inputs themselves (missation intervals well-formed, seq_deltas
consistent with reference). Then checks that `globally_missing_sites` equals the intersection
of all tips' missing intervals, and that the sequence at each tip (reconstructed by composing
arc deltas from the focus via an Euler tour) matches `tip_descs[i].seq_deltas` at non-missing
sites.

## Guide tree construction

The entry point is `build_guide_tree(ref_sequence, tip_descs, bitgen) -> Utree`. The
incremental construction logic lives in a `Utree_builder` class (defined in `core/utree.cpp`,
not exposed in the header). The builder's `add_tip(X)` method adds a single tip X to the
tree; tips can be added in any order, so `build_guide_tree` adds them in input order
(0, 1, ..., N-1) while a future Round 3 can reuse the same machinery with a different
insertion order. `build_guide_tree` is a thin wrapper: create a builder, add each tip in
turn, call `finish()` to get the tree.

### Algorithm

**Invariant**: After adding the first K tips, the `Utree` is in a valid state representing
exactly those K tips. There is no final "clean-up" step after the last tip — the `Utree` is
immediately usable after each insertion.

**Starting point (first tip):**

Initialize the `Utree` with a single tip X and no edges (`add_first_tip`). This is the one
moment when a node has zero outgoing arcs. Set the focus at X and initialize
`deltas_ref_to_focus` from `tip_descs[X].seq_deltas`. Initialize `globally_missing_sites` to
`tip_descs[X].missations.intervals`.

**Adding subsequent tips:**

The remaining tips are added one at a time. For each new tip X:

1. **Update `globally_missing_sites`** (`update_globally_missing_sites`). Remove from
   `globally_missing_sites` any site where X is informative (not missing). For each such
   newly-informative site that is also in `tip_descs[X].seq_deltas` (i.e., X's state differs
   from the reference at that site), add a corresponding entry to `deltas_ref_to_focus` — this
   reflects the implicit decision that the tree's state at previously-globally-missing sites
   now matches X's state. (Only sites that are both leaving `globally_missing_sites` and in
   `tip_descs[X].seq_deltas` need this update, so it is cheap.) As a fast path, if
   `globally_missing_sites` is already a subset of tip X's missing intervals (checked via
   `interval_set_is_subset_of`), the intersection and `deltas_ref_to_focus` update can be
   skipped entirely.

2. **Initialize `focus_to_X_deltas`** (`init_focus_to_X_deltas`). Both `deltas_ref_to_focus`
   and `tip_descs[X].seq_deltas` describe deltas from the same starting point (the reference
   sequence). The delta from the focus to X is their "difference": start with a copy of one,
   then pop or push all entries from the other. This is analogous to vector subtraction
   (focus_to_X = ref_to_X - ref_to_focus). Build a `Heap_site_deltas` from
   `tip_descs[X].seq_deltas` (CHECKing that none of these entries are at sites missing at
   X — `seq_deltas` should only contain non-missing sites). Then for each `(s, {R, F})` in
   `deltas_ref_to_focus` where s is not missing at X, call
   `push_front_site_deltas({s, F, R}, focus_to_X_deltas)` — this prepends the inverse of
   each ref→focus delta, converting ref→X into focus→X.

3. **Search for the best attachment edge** (`find_best_attachment_arc`) **via approximate
   branch-and-bound.** Use a min-priority queue of arcs, keyed by the cost of attaching X at
   that arc's edge. The search starts at the focus and explores candidate edges in order of
   increasing cost, pruning aggressively.

   The intuition behind the branch-and-bound search: suppose the partial tree is correct and
   the true best attachment point is at some edge P, but the focus is currently at Q. The
   path from Q to P is likely parsimonious — a single mutation for every site where Q and P
   differ and no mutation for sites where they agree. Moreover, the arcs branching off the
   Q-P path have small total evolutionary distance, so the chance of "red herring" mutations
   (at sites where Q-P or P-X differ) along those side branches is low. Following the
   direction of steepest cost decrease at each step leads straight from Q to P. A small
   pruning threshold accommodates the occasional stray mutation, while reflecting the
   exponentially decreasing probability of many such coincidences.

   - **Cost of attaching X at an edge** (`eval_focal_arc`)**:** Attaching X at edge (A, B)
     splits it by inserting a new inner node M between A and B, then connecting X to M. The
     existing site deltas on (A, B) are distributed between A-M and M-B to minimize the
     number of new site deltas on M-X: each delta is placed on whichever side makes M's
     state match X at that site. The cost equals the number of remaining site deltas on M-X
     — i.e., the increase in total tree parsimony score. Evaluating the cost at a focal arc
     is cheap since it uses `focus_to_X_deltas` directly without moving the focus:
     ```
     cost = ssize(focus_to_X_deltas) - savings
     savings = |{sites in A->B arc deltas where X matches B's state}|
     ```
     For each site `s` in the A->B arc's deltas: look up `s` in `focus_to_X_deltas`. If
     present and `focus_to_X_deltas[s].to == delta.to` (i.e., the A->B arc changes to the
     same state that X has), X agrees with B — count one saving, because this delta can be
     placed on the A-M side (giving M the same state as B and X at this site), avoiding a
     delta on M-X. That's it: sites where X agrees with A (the focus) are *not* in
     `focus_to_X_deltas`, so they don't contribute to the baseline
     `ssize(focus_to_X_deltas)` and need no saving. The edge can only help, never hurt.
     Update `best_cost` and `best_arc` whenever `cost <= best_cost` (using `<=`, not `<`,
     so that ties always prefer a real arc over the `k_no_arc` default).

   - **Initialization**: Evaluate all edges incident to the focus. Update `best_cost` and
     `best_arc`. Push all arcs leaving the focus onto the priority queue, keyed by the cost
     of attaching X at the corresponding edge. When only one tip has been added so far, the
     focus has no arcs yet, so the queue starts empty and the search immediately terminates
     with `best_arc = k_no_arc`.

   - **Processing**: Pop the lowest-priority arc `R` from the queue. Its priority is the
     cost of attaching X at R's edge, which was already evaluated when R was pushed. If the
     priority exceeds `best_cost + threshold`, stop the search (everything remaining has
     equal or higher priority). Otherwise, move the focus to R's target (updating
     `deltas_ref_to_focus` and `focus_to_X_deltas` as the focus moves — see below). For
     each outgoing arc from R's target except R's mate, evaluate the cost of attaching X at
     that arc's edge, update `best_cost` and `best_arc`, and push the arc onto the queue
     keyed by its cost. (Pruning happens at pop time, not push time: arcs whose priority
     exceeds the threshold are simply discarded when popped.)

   - **Pruning threshold**: Use a low threshold (1-2) for the guide tree — it's OK to
     occasionally pick a slightly-suboptimal attachment point in exchange for faster
     construction.

   - **Focus movement and `focus_to_X_deltas` updates**
     (`move_focus_updating_focus_to_X_deltas`)**:** When the focus moves from node F to
     neighbor G along arc F->G, update `deltas_ref_to_focus` by applying the arc's deltas.
     Update `focus_to_X_deltas` in tandem: for each non-missing-site entry in the arc's
     deltas, apply `pop_front_site_deltas` to `focus_to_X_deltas`. This keeps
     `focus_to_X_deltas` consistent as the focus moves, without recomputing from scratch.

4. **Attach X.** Two cases:

   - **`best_arc == k_no_arc`** (second tip only — no edges exist yet)
     (`attach_tip_directly_to_isolated_focus`): add a direct edge from the focus to tip X
     with `focus_to_X_deltas` as its deltas.

   - **`best_arc != k_no_arc`** (third tip onward): move the focus to `origin(best_arc)`
     (`move_focus_updating_focus_to_X_deltas`), making `best_arc` a focal arc, then call
     `attach_tip_to_focal_arc`, which calls `split_edge(best_arc, M, site_delta_side)` where
     M is the next available inner node. `site_delta_side(delta, A, B)` distributes each
     site delta:
     - X has state matching A's state: delta goes to M-B side (return B).
     - X has state matching B's state: delta goes to A-M side (return A).
     - X matches neither: delta goes to a randomly chosen side (using `bitgen`).
     - Site missing at X or globally missing: delta goes to a randomly chosen side.

     Then create the M-X edge with deltas for each site where X differs from M's
     (now-determined) sequence.

   In both cases, set `nodes[X].arc_to_focus` to the arc from X toward its neighbor (the
   focus or M, respectively).

5. **Leave the focus** near the attachment point. The next tip's search will start from
   there.

### Globally missing sites

Some sites may be missing across the *entire* partial tree (no tip added so far is
informative at that site). These sites are tracked in `globally_missing_sites` (an
`Interval_set`) on the `Utree` itself.

Globally missing sites are excluded from `focus_to_X_deltas` and from edge cost calculations,
just like sites missing at X itself.

### Summary (`add_tip` pseudocode)

```
add_tip(X):
  if tips_added_ == 0:
      add_first_tip(X)
  else:
      update_globally_missing_sites(X)
      init_focus_to_X_deltas(X)
      best_arc = find_best_attachment_arc(X)
      if best_arc == k_no_arc:
          attach_tip_directly_to_isolated_focus(X)
      else:
          move_focus_updating_focus_to_X_deltas(origin(best_arc), X)
          attach_tip_to_focal_arc(X, best_arc)
  ++tips_added_
```

### Output

A `Utree` with all N tips placed. The topology is rough (no refinement), but nearby tips in
this tree are likely genuinely similar. In Round 2, this will be converted to a `Phylo_tree`
via midpoint rooting.

## Files to create/modify

| File | Action | What |
|------|--------|------|
| `core/site_deltas.h` | Modify | Add `Heap_site_deltas` type alias; templatize composition operations |
| `core/site_deltas.cpp` | Modify | Remove bodies of templatized functions (now in header) |
| `core/interval_set.h` | Modify | Add `interval_set_is_subset_of`; add `Interval_set::assert_valid` member method |
| `core/utree.h` | Create | `Utree`, `Unode`, `Uarc`, `Arc_direction`, `Annotated_arc`, `split_edge`, `move_focus_to` (all overloads), `annotated_arc_euler_tour`, `arc_euler_tour`, `target(arc)` helper |
| `core/utree.cpp` | Create | `Utree::reset_focus`, `Utree_builder` (internal class), `build_guide_tree`, `assert_utree_integrity`, `assert_utree_matches_tip_descs` |
| `core/CMakeLists.txt` | Modify | Add `utree.cpp` to sources |
| `tests/utree_tests.cpp` | Create | Unit tests for `Utree` operations and guide tree construction |
| `tests/interval_set_tests.cpp` | Modify | Add tests for `interval_set_is_subset_of` |
| `tests/CMakeLists.txt` | Modify | Add `utree_tests.cpp` to test sources |

## Reused components

- `Site_deltas`, `Site_delta` (`core/site_deltas.h`) — scratch-allocated edge deltas and
  composition operations.
- `push_front_site_deltas`, `pop_front_site_deltas`, `push_back_site_deltas`,
  `append_site_deltas` (`core/site_deltas.h`) — traversal and edge merging.
- `Tip_desc` (`core/phylo_tree.h:137`) — input format.
- `Seq_delta` (`core/mutations.h:67`) — per-site sequence differences.
- `Interval_set` (`core/interval_set.h`) — for `globally_missing_sites` and missation handling.

## Verification

1. **Test infrastructure**: Set up helper functions to build small `Utree`s from inline
   descriptions (e.g., a list of tip sequences or seq_deltas), so that unit tests are
   concise and easy to write.

2. **Unit tests**: Test `split_edge` and `move_focus_to` on small examples. Test
   `annotated_arc_euler_tour` and `arc_euler_tour` on a manually-built 3-tip tree by
   comparing the exact sequence of `{origin, target, direction}` triples against the
   expected DFS order. Test `assert_utree_integrity` and `assert_utree_matches_tip_descs`
   on valid trees (should pass) and on deliberately corrupted trees (should CHECK-fail via
   `EXPECT_DEATH`). Test guide tree construction on small examples (3-5 tips), verify the
   total parsimony score is optimal or near-optimal, and call `assert_utree_integrity` and
   `assert_utree_matches_tip_descs` on each result.

3. **Build**: `cmake --build build/debug --target tests -j$(nproc) && build/debug/tests/tests --gtest_filter='Utree*'`
