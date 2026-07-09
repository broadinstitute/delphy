# Better Tree Init, Round 6: Unified SPR Refinement

This is the sixth in a series of plans implementing a new tree initialization pipeline for
Delphy.  See the [Round 1 plan](2026-05-04-01-better-tree-init-round1-utree-and-guide-tree.md),
section "Roadmap", for the full multi-round plan.

This round unifies tip and subtree SPR refinement under a single framework based on
local Fitch parsimony.  A compact representation (`Relative_fitch_sets`) enables
O(|arc deltas|) candidate evaluation for both cases.  The main changes are:
(1) a generic search algorithm (Move/Eval/Search), (2) subtree SPR initialization and
reattachment, and (3) refactoring `add_tip` to use the same infrastructure.

## Background

Round 5 added SPR refinement of individual tips: detach a tip, search for the best
reattachment edge, and move it there.  This corrects misplaced tips but cannot fix errors
deeper in the topology.  If an entire subtree is attached at the wrong place, it's unlikely that
a sequence of single-tip moves can fix it efficiently.

## Goal

Add subtree SPR refinement.  When regrafting a subtree, the attachment branch and its
immediate sibling branches are rewired to achieve local parsimony.  In principle, the SPR
refinement of tips also implements local parsimony, so it's possible to unify these two
refinement modalities.  This plan develops that unified framework by abstracting either
side of the attachment branch into a single sequence of per-site Fitch sets; the search
for a better attachment location is then generic, while the precise rewiring of mutations
around the attachment branch remains specific to each modality.


## Subtree SPR: overview

### The 5-edge star

Pick a random internal edge M-X where both endpoints have degree 3.  X's other two neighbors
are D and E:

```
... -- M -- X -- D
       |    |
      ...   E
```

Detach: remove the M-X edge, then merge through M (joining M's other two neighbors P and Q
into a direct P-Q edge) and merge through X (joining D and E into a direct D-E edge).  Both M
and X become free connector nodes.  The tree splits into an M-side component (the
complementary tree, which contains the focus) and an X-side component (D's and E's subtrees
joined by the new D-E edge).  (The detach also re-roots the X-side's internal focus links
onto one child before merging through X; this is a focus-bookkeeping detail covered by the
detailed flow below.)

Throughout, **A** and **B** denote the endpoints of the candidate edge M is reattached into,
**P** and **Q** denote M's two original neighbors, and **D** and **E** denote X's two
children.  The tree is unrooted, so D and E are, strictly, X's two *neighbors* other than M;
we call them "children" because the focus lies on the M side, orienting the M-X edge so that
D and E hang below X.  This plan uses "neighbors" and "children" interchangeably for them.

### Reattachment

Search the M-side for the best edge A-B on which to reattach.  Reattaching creates a 5-edge
star:

```
A -- M -- X -- D
     |    |
     B    E
```

M is inserted into A-B and X is re-inserted into the D-E edge (both via `split_edge`); M and
X are then connected by a new edge.  Splitting D-E chooses X's state and distributes the
existing D-E deltas onto the new X-D and X-E arcs, giving local parsimony for free.

### Comparison with tip SPR

In Round 5, the moved node X is a tip with no other neighbors.  There is no D-E subtree; the
star has 3 edges:

```
A -- M -- X
     |
     B
```

The search problem is the same in both cases: find the best A-B edge on the M-side.


## Local Fitch parsimony at a candidate edge

SPR refinement prunes a subtree rooted at X and searches for where to reattach it.  At each
candidate, we insert a connector node M into some edge A-B and wire M to X.  We assign a
**Fitch set** to each side of the M-X edge: Fitch_X for X's subtree, Fitch_M for the
candidate attachment point.  Each is the set of states its endpoint can take while
maintaining local parsimony within its side.  Fitch_X is computed once when X is detached
and does not change across candidates; Fitch_M varies with each candidate.  The concrete
definitions appear in the tables below.

**The M-X edge is all that matters.**  An M-X mutation occurs at a site exactly when Fitch_M
and Fitch_X don't intersect.  If we pick x from Fitch_X and m from Fitch_M, then everything
on either side of the branch is locally parsimonious by construction.  So the search reduces
to: find the candidate that minimizes
**n_mismatch = |{sites where Fitch_M and Fitch_X do not intersect}|**.

Tips and subtrees share the same search problem: minimize n_mismatch.  They differ only in
the per-site Fitch_X.

**Incremental evaluation.**  Given Fitch_M at one candidate attachment arc, Fitch_M at
neighboring candidates can be computed with only a few operations per arc delta -- the
endpoint states change incrementally as we traverse the tree.  This makes the search
O(|arc deltas|) per candidate rather than O(L).


### Fitch_X (fixed for a given X)

Fitch_X is computed once when X is detached.  At each site, Fitch_X is obtained from X's
children via the standard Fitch rule: Fitch_X = Fitch_D & Fitch_E if the intersection is
nonempty, otherwise Fitch_X = Fitch_D | Fitch_E.  Here, Fitch_D = {d} when D is an inner
node or a tip with data at the site, and Fitch_D = {A,C,G,T} when D is a tip missing the
site; likewise for Fitch_E.  This is what makes the parsimony _local_: Fitch_D depends only
on D's state, not on the subtree rooted at D away from X.  For a tip X with no children,
Fitch_X = {x} at non-missing sites and {A,C,G,T} at missing sites.

Concrete cases:

| Scenario                                                             | Fitch_X                | |Fitch_X| |
|----------------------------------------------------------------------|------------------------|-----------|
| Tip, site has data                                                   | {x}                    | 1         |
| Tip, site is missing                                                 | {A,C,G,T}              | 4         |
| Subtree, d = e, both children inner nodes or tips with data at site  | {d}                    | 1         |
| Subtree, d != e, both children inner nodes or tips with data at site | {d, e}                 | 2         |
| Subtree, exactly one child is a tip missing this site                | {state of other child} | 1         |
| Subtree, both children are tips, both missing site                   | {A,C,G,T}              | 4         |


### Fitch_M (determined by the candidate edge)

At each site, Fitch_M is computed from the states at A and B via the same Fitch rule:
Fitch_A = {a} for inner nodes or tips with data, {A,C,G,T} for tips at missing sites; and
likewise for Fitch_B.

In the common case where both A and B are inner nodes or tips with data at the site:

| Condition | Fitch_M |
|-----------|---------|
| a = b     | {a}     |
| a != b    | {a, b}  |

A site has a != b exactly when the candidate edge A-B carries an arc delta at that site.

When A or B is a tip missing the site, Fitch_M differs: the missing side's {A,C,G,T}
intersected with the other side's singleton always gives the other side's state.  For
example, if A is a tip missing site s, then Fitch_M = {b} regardless of a.  As we will see
below, this does not cause any complications: an existing invariant ensures that tip-adjacent
arcs carry no deltas at sites missing at the tip, so a != b never occurs at a site missing at a tip endpoint.


## Fitch sets relative to a reference sequence

Throughout this plan, we need to represent a per-site Fitch set -- a set of nucleotide
states at each site -- and test whether a given state belongs to it.  We do this relative
to a **reference sequence**: a sequence of states (one per site) that serves as a baseline
for compact storage.  The reference sequence is not stored physically; only deviations from
it appear explicitly.

The key insight is that per-site Fitch sets are almost always singletons.  We exploit this
by partitioning sites into three categories by Fitch set size, each with its own
representation:

1. **`resolved_deltas_`**: `Heap_site_deltas`
   (= `absl::flat_hash_map<Site_index, Site_delta>`).
   Covers **resolved sites** (Fitch set size 1).  Contains an entry `{from: r, to: f}` at
   site s when the reference state `r` differs from the sole element `f`.  Absence means the
   reference state equals f (the common sub-case, so the map is sparse).

2. **`ambiguous_masks_`**: `absl::flat_hash_map<Site_index, Seq_letter>`.
   Covers **ambiguous sites** (Fitch set size 2 or 3).  Value: a `Seq_letter` bitmask
   (e.g., `to_seq_letter(d) | to_seq_letter(e)` for a size-2 set).
   **Static**: does not depend on the reference sequence.

3. **`uninformative_sites_`**: `Interval_set`.
   Covers **uninformative sites** (Fitch set size 4), including globally missing sites.
   Any state belongs to the Fitch set.  Missing data tends to clump into contiguous
   intervals, so a few intervals cover many sites.  **Static**: does not depend on the
   reference sequence.

Most genomic sites are resolved (size-1 Fitch sets), ambiguous sites (size-2/3) are
rare, and uninformative sites (size-4) cluster into intervals---so this hybrid
representation is both space-efficient and fast.

These three components partition all sites: each site belongs to
exactly one.  We encapsulate them in a structure called `Relative_fitch_sets`, which
supports two operations:

**`contains(s, state, ref_state) -> bool`**: Is `state` in the Fitch set at site s?
The caller supplies `ref_state`, the current reference state at site s (in `Utree_builder`,
this is always known from the arc delta being processed: it's the delta's `from` state).
1. If `uninformative_sites_.contains(s)`: yes.
2. Let `mask = ambiguous_masks_[s]`.  If found: return `mask & to_seq_letter(state)`.
3. Otherwise (resolved site): look up `resolved_deltas_[s]`.  If found, return
   `state == .to`.  If absent, the sole element equals `ref_state`; return
   `state == ref_state`.

**`on_ref_change(s, old_state, new_state)`**: adjust the representation when the reference
state changes at s.  For resolved sites, this calls `pop_front_site_deltas` on
`resolved_deltas_`.  For ambiguous and uninformative sites, no adjustment is needed (their
representations are absolute).

### Usage in SPR search

In the SPR search, we build a `Relative_fitch_sets` to represent Fitch_X.  The reference
sequence is the sequence at the focus.  The builder already maintains this sequence
incrementally as `deltas_ref_to_focus` -- an existing `Utree` member holding the running
per-site delta from the (physical) reference sequence to the current focus, so the focus
state at site s is `deltas_ref_to_focus[s].to` if present, else the reference state at s.
As the focus moves over an arc, each arc delta triggers one `on_ref_change` call, keeping
the `Relative_fitch_sets` representation consistent.

This same structure works for both tips (`ambiguous_masks_` empty, the other two populated
from the tip descriptor) and subtrees (all three populated from the children's arcs).  The
search algorithm -- move, eval, branch-and-bound -- queries only `contains` and
`on_ref_change`, so it is generic over the kind of X being moved.  The same structure could
also serve `add_tip` in the existing builder, unifying the Round 5 infrastructure with the
new subtree SPR.


## Incremental evaluation of n_mismatch

To evaluate a candidate arc A->B, we first move the focus to node A, incrementally updating
the Fitch_X representation and a running mismatch counter along the way.  Then we examine
the arc deltas on A->B itself.

### n_mismatches_ (running counter)

As the focus moves through the tree, we maintain:

```
n_mismatches_ = |{sites where the state at the focus is not in Fitch_X}|
```

This is the n_mismatch value for a hypothetical candidate arc with no deltas (Fitch_M =
{state at focus} at every site).  It depends only on the focus position and is updated in
O(|arc deltas|) each time the focus hops along an edge: for each arc delta
`{site, old_state -> new_state}`, call `contains` for the old and new states and adjust the
counter accordingly.

### savings (per candidate arc)

At an arc delta site on A-B with delta a -> b, Fitch_M = {a, b} instead of {a}.  Since
{a, b} is a superset of {a}, expanding Fitch_M can only resolve a mismatch, never create
one.  A savings of 1 occurs at a site where:

1. n_mismatches_ counted it: `not contains(s, a, a)`, AND
2. b resolves it: `contains(s, b, a)`.

```
savings = |{arc delta sites on A-B where not contains(s, a, a) AND contains(s, b, a)}|
```

Computed in O(|arc deltas|) per candidate, with at most two `contains` calls (and hence at
most two hash lookups) per delta site.

### Candidate cost

```
n_mismatch = n_mismatches_ - savings
```

### Tip endpoints with missing data

The savings formula assumes Fitch_M = {a, b} at arc delta sites.  This requires both A and
B to have data at the site.  If A is a tip missing site s, Fitch_A = {A,C,G,T} and
Fitch_M = {b} instead of {a, b} -- a narrower set that could cause the formula to
undercount.

However, this situation cannot arise.  There is an invariant (enforced by the bugfix in
commit eab7a69) that arcs incident on tips carry no deltas at sites missing at the tip.
The only edge where this invariant could be violated is the merged edge created by
`merge_through` during detach -- the composition of A-M and M-B deltas can inherit deltas
at sites missing at A or B.  `strip_missing_deltas` (see below) removes these after every
merge, restoring the invariant on the merged edge; with that step in place, every candidate
arc respects the invariant and the formula is exact as stated.  (Removing them is subtler
than a plain erase -- such a delta can carry the *other* endpoint's real data -- so
`strip_missing_deltas` "slides" it off rather than dropping it; see that helper for details.)


## Core operations

Three primitives -- Move, Eval, Search -- are shared by tip and subtree SPR.  They read only
`fitch_X_` and `n_mismatches_`, so they are generic over how `fitch_X_` was initialized.
Each is O(|arc deltas|) per hop/candidate (at most two `contains` calls per delta site).

- **Move** (`move_focus_updating_fitch_X`): hop the focus via `tree_.move_focus_to(target,
  callback)`.  Per arc delta `{s, old -> new}`, apply the `n_mismatches_` update from
  "Incremental evaluation" (still keyed on reference state `old`), then
  `on_ref_change(s, old, new)`.
- **Eval** (`eval_focal_arc`): for a focal arc (origin == focus), return
  `n_mismatches_ - savings` using the savings count from "Incremental evaluation".
- **Search** (`find_best_attachment_arc`): unchanged from Round 5 -- priority-queue
  branch-and-bound that calls Move to traverse and Eval to score, pruning subtrees whose base
  cost exceeds the best n_mismatch so far.  Blind to tip vs subtree.


## Application: tip SPR

### Fitch_X for tips

At non-missing sites: Fitch_X = {x} (the tip's observed state), size 1.
At missing sites: Fitch_X = {A, C, G, T}, size 4.

`Relative_fitch_sets` after initialization:
- `resolved_deltas_`: entries at non-missing sites where focus != x.
- `ambiguous_masks_`: **empty** (no size-2+ sites).
- `uninformative_sites_`: X's missing-site intervals from the tip descriptor.
- `n_mismatches_` = |resolved_deltas_|.

### Initialization

Initialize `Relative_fitch_sets` **relative to the current focus**, wherever it happens to
be, directly from X's tip descriptor and `deltas_ref_to_focus`:

- `resolved_deltas_` = the focus->X deltas at non-missing sites.  This is the composition
  `(ref -> X) o (focus -> ref)`, exactly the existing `init_focus_to_X_deltas` computation.
- `ambiguous_masks_` = empty.
- `uninformative_sites_` = X's missing-site intervals from the tip descriptor.
- `n_mismatches_` = |resolved_deltas_|.

For tip SPR the focus is at M at this point, so `resolved_deltas_` equals the M-X arc deltas
(by path consistency, `focus -> X = M -> X`, and the M-X arc, by the tip-adjacency invariant,
carries no deltas at sites missing at X).  For `add_tip` the focus is wherever the previous
insertion left it and X has no edges yet; the same computation applies unchanged.  A single
focus-relative init thus serves both.

(Considered and rejected: moving the focus onto X first so that `resolved_deltas_` starts
empty, then letting the X->M focus move populate it.  That produces the identical
`Relative_fitch_sets` but requires an extra focus traversal to X, cannot serve `add_tip`
-- where X is not yet in the tree -- and adds a tip-only init path that diverges from the
subtree init.  The focus-relative init above is the unifying choice.)

Then detach X from M, move focus to one of M's former neighbors (updating
`Relative_fitch_sets` via `on_ref_change` along the way), merge through M (creating a
direct edge between M's other two neighbors), and strip missing-site deltas from the
merged edge.

### Reattachment

After finding the best edge A-B: split A-B to insert M (via `split_edge`).  `split_edge`'s
callback decides, per arc delta, which side of the inserted node the delta stays on: it
returns A to keep the delta on the A-M arc, or B to keep it on the M-B arc.  The inserted
node M takes the state of the *opposite* endpoint -- returning A leaves the delta between A
and M, so M ends up with B's state, and vice versa.  (This is the existing `split_edge`
contract from Round 5; the `return A // M gets b` comments in the pseudocode below follow
directly from it.)  We place each arc delta on the side that gives M a state matching x, to
avoid an M-X mutation.  For each delta `{s, a -> b}` on A-B:

- If `contains(s, b, a)`: delta on A-M side -> M gets B's state = b.
- If `contains(s, a, a)`: delta on M-B side -> M gets A's state = a.
- If neither: either side (M-X mutation unavoidable).

Wire M-X directly.  M-X deltas: sites where M's assigned state != x.

### Reduction to Round 5

With `ambiguous_masks_` empty, `contains(s, state, ref_state)` always falls through to the
resolved path: check `resolved_deltas_[s]`.  The eval reduces to "savings if b == f" per
delta site.  This IS the Round 5 algorithm, plus `uninformative_sites_` for tips with
missing data.


## Application: add_tip

`add_tip(X)` inserts a new tip into an existing tree.  X has no edges yet, so there is no
detach step; it reuses tip SPR's search and reattach and differs only in init.  Since
`init_fitch_X_for_tip` is already focus-relative, it produces the correct
`Relative_fitch_sets` unchanged -- just from the tip descriptor rather than an M-X arc.

This subsumes the existing `add_tip`, which used `focus_to_X_deltas_` and its dedicated
methods directly; those are removed.


## Application: subtree SPR (5-edge star)

### Fitch_X for subtrees

X has children D and E.  Fitch_X is computed from their states at initialization, before X
is merged out.

- **Size-1 sites** (d = e): Fitch_X = {d}.  The vast majority.
- **Size-2 sites** (d != e): Fitch_X = {d, e}.
- **Size-4 sites**: both D and E are tips missing this site.

`Relative_fitch_sets` after initialization:
- `resolved_deltas_`: entries at size-1 sites where focus != f (= d = e).
- `ambiguous_masks_`: one entry per site where d != e.  mask = `to_seq_letter(d) | to_seq_letter(e)`.
- `uninformative_sites_`: intersection of D's and E's missing intervals (if both are tips).
- `n_mismatches_`: count across all sizes.

### Initialization

With focus at M (before detaching), walk the arcs M-X, X-D, and X-E.  For each site with
a delta on any of these arcs, record the state at M (the first `from` encountered for that
site) and propagate through the deltas to get the states at D and E.  The result is a map
of `s -> (m, d, e)` for all sites where not all three states are equal.  (Sites missing at
both D and E, if both are tips, go into `uninformative_sites_` instead.)

Build `Relative_fitch_sets` from this map in a single pass:

- **d = e (resolved):** Fitch_X = {d}.  Add `{from: m, to: d}` to `resolved_deltas_`.
  m != d (otherwise all three would be equal), so increment `n_mismatches_`.
- **d != e (ambiguous):** Fitch_X = {d, e}.  Store mask `to_seq_letter(d) | to_seq_letter(e)` in
  `ambiguous_masks_`.  If m is not in {d, e}, increment `n_mismatches_`.

Then remove the M-X edge, re-root the X-side onto one child so X keeps a focus link into
`merge_through` (see the detach flow and the `merge_through` precondition), and merge through
X, producing a direct D-E edge whose deltas encode exactly the `d != e` (ambiguous) sites --
one delta `{from: d, to: e}` per such site.
This D-E edge is the *ordered* record of which state belongs to D and which to E, which
reattachment needs (`ambiguous_masks_` alone is unordered).  Then move the focus to one of
M's former neighbors (updating `Relative_fitch_sets` via `on_ref_change` along the way),
merge through M, and strip missing-site deltas from both merged edges.

### Reattachment

After finding the best edge A-B:

**Step 1: Split A-B to insert M.**  Identical to tip reattachment above -- same `split_edge`
callback, now querying a subtree's Fitch_X -- placing each delta on the side that gives M a
state in Fitch_X.

**Step 2: Split the D-E edge to re-insert X.**  A second `split_edge` re-inserts X into the
D-E edge.  Its callback chooses X's state per site, preferring M's state when it lies in
Fitch_X (avoiding an M-X mutation) and distributing each D-E delta onto the new X-D or X-E
arc accordingly.  Because the D-E edge already carries exactly the ambiguous-site deltas,
this reuses `split_edge` instead of recomputing X-D and X-E by hand.

**Step 3: Wire M-X and set `arc_to_focus`.**  Connect M and X with a new edge carrying a
delta wherever M's state differs from the chosen X state.  Set X's `arc_to_focus` toward M
(hence toward the focus).  Then repoint both D and E toward X.  After Step 2's `split_edge`,
one child already has `arc_to_focus` pointing at X (set by split_edge's update), while the
other (the X-side sink, which was D) still has `arc_to_focus == k_no_arc`.  Rather than
tracking which child was the sink, Step 3 simply repoints both unconditionally — the
already-correct child is written with the same value, and the sink child is fixed.  Both
children then point at X and X points at M, giving the whole X-side a consistent orientation
toward the focus.


## Delta tracking

The `spr_refine` loop needs to know the change in total tree deltas after each move to abort
early after N consecutive non-improvements.  The change has the same form in both SPR
applications:

```
delta_change = best_cost - old_cost
```

where `best_cost` is the eval result (the n_mismatch value) at the chosen reattachment arc
-- the same `best_cost` returned by the search in the outer-loop pseudocode below -- and
`old_cost` is the net delta cost of X at its old position.

**add_tip:** shown here only to complete the picture -- `add_tip` returns `void` and runs in
the guide-tree / refined-tree builders, *not* in the `spr_refine` abort-after-N loop, so it
wires no delta tracking.  For the record, its change degenerates to `delta_change = best_cost`
because `old_cost = 0` (X has no old edges).

**Tip SPR:** `old_cost = d_MX + d_MP + d_MQ - d_PQ`, where d_MX, d_MP, d_MQ are the delta
counts on the three edges around M before detaching, and d_PQ is the delta count on the
merged P-Q edge after `merge_through`.

**Subtree SPR:** `old_cost = d_MX + d_MP + d_MQ + d_XD + d_XE - d_PQ - d_DE`, where d_PQ and
d_DE are the delta counts on the merged P-Q and D-E edges (after stripping).  Each merged
edge cancels mutations that were needlessly split across the two edges it replaces; d_DE
equals the number of ambiguous sites (one delta per site where d != e, since the optimal x
is always in {d, e}).

The outer refinement loop tracks `total_deltas += delta_change` and counts consecutive
non-improvements (delta_change >= 0).  After N consecutive non-improvements, it stops.


## Detailed implementation plan

### New type: `Relative_fitch_sets` (in `core/utree.cpp`)

An implementation-detail struct next to `Utree_builder` (not in the header):

```cpp
struct Relative_fitch_sets {
  Heap_site_deltas resolved_deltas_;                     // size-1 sites: {from: ref, to: f}
  absl::flat_hash_map<Site_index, Seq_letter> ambiguous_masks_;  // size-2/3 sites: bitmask
  Interval_set<> uninformative_sites_;                   // size-4 sites (incl. globally missing)

  auto clear() -> void;

  // Is `state` in the Fitch set at site s?  `ref_state` is the focus's current state at s.
  auto contains(Site_index s, Real_seq_letter state, Real_seq_letter ref_state) const -> bool {
    if (uninformative_sites_.contains(s)) { return true; }
    if (auto it = ambiguous_masks_.find(s); it != ambiguous_masks_.end()) {
      return it->second & to_seq_letter(state);
    }
    if (auto it = resolved_deltas_.find(s); it != resolved_deltas_.end()) {
      return state == it->second.to;
    }
    return state == ref_state;   // sole element equals the (implicit) reference state
  }

  // The focus's reference state at s changed from old_state to new_state.
  auto on_ref_change(Site_index s, Real_seq_letter old_state, Real_seq_letter new_state) -> void {
    // Ambiguous and uninformative representations are absolute -> nothing to do.
    if (uninformative_sites_.contains(s)) { return; }
    if (ambiguous_masks_.contains(s)) { return; }
    pop_front_site_deltas({s, old_state, new_state}, resolved_deltas_);
  }
};
```

Note: `contains`/`on_ref_change` check `uninformative_sites_` and `ambiguous_masks_` first,
so uninformative and ambiguous sites are handled automatically -- the Move and Eval core
operations need **no** explicit `globally_missing`/`missing-at-X` guards (globally missing
sites never appear in arc deltas anyway).  This is a simplification over Round 5, which
guarded every arc delta against `miss_X` and `globally_missing_sites`.

The correctness of dropping those guards rests on the invariant that arc deltas never fall
on missing sites.  Where cheap and unintrusive (e.g. in the Move/Eval loops and the init
walks), add debug `CHECK`s that each arc-delta site is neither in `globally_missing_sites`
nor in the moved node's missing sites -- turning a silent representation bug into a loud one.

### `Utree_builder` changes

Remove `focus_to_X_deltas_` and its dedicated methods (`init_focus_to_X_deltas`,
`move_focus_updating_focus_to_X_deltas`, tip-specific `eval_focal_arc`).  Add:

New members:
- `Relative_fitch_sets fitch_X_;`           -- Fitch_X for the node being (re)placed
- `int n_mismatches_ = 0;`                  -- running counter (see "Incremental evaluation")
- `absl::flat_hash_map<Site_index, Real_seq_letter> m_overrides_;`  -- scratch: M's state at arc-delta sites where M departs from the focus (filled by `split_best_arc_inserting_M`, read by subtree Step 2)
- `int component_abort_threshold_ = 0;`     -- DFS abort threshold for `pick_random_node_in_component` (see "SPR throughput" section)
- `std::vector<Node_index> component_dfs_stack_;`  -- scratch DFS stack for same
- `std::vector<Node_index> component_nodes_;`      -- scratch node list for same
- keep `M_to_X_deltas_`, `pq_`, `best_arcs_`.

New / reworked methods (all public, per Round 5's convention for this .cpp-only class):

| Method | Purpose |
|--------|---------|
| `init_fitch_X_for_tip(int X)` | Build `fitch_X_` + `n_mismatches_` for tip X, relative to the current focus (subsumes `init_focus_to_X_deltas`; also serves `add_tip`). |
| `init_fitch_X_for_subtree(X)` | Build `fitch_X_` + `n_mismatches_` for subtree X (focus at M); derives X's children D, E as its neighbors other than the focus M. |
| `move_focus_updating_fitch_X(target)` | Move core operation.  No X parameter (uninformative sites live in `fitch_X_`). |
| `eval_focal_arc(arc) -> int` | Eval core operation, generalized to `contains`. |
| `find_best_attachment_arc() -> pair<Arc_index,int>` | Search core operation.  No X parameter; `best_cost` seeds from `n_mismatches_`. |
| `focus_state(s) -> Real_seq_letter` | Accessor: the focus's state at site s (`deltas_ref_to_focus[s].to` if present, else the reference state). |
| `m_state(s) -> Real_seq_letter` | Accessor: M's state at site s during subtree reattach (`m_overrides_[s]` if present, else `focus_state(s)`). |
| `split_best_arc_inserting_M(best_arc, M)` | Shared reattachment helper: split `best_arc` to insert connector M, distributing deltas to minimize M-X cost; builds `M_to_X_deltas_` and `m_overrides_`. |
| `wire_M_X(M, X)` | Shared reattachment helper: add the M-X edge from `M_to_X_deltas_` and orient X toward the focus (via M). |
| `attach_tip_to_focal_arc(X, best_arc, M)` | Tip reattachment: `split_best_arc_inserting_M` then `wire_M_X`. |
| `attach_subtree_to_focal_arc(X, best_arc, M, arc_DE)` | Subtree reattachment: `split_best_arc_inserting_M`, split `arc_DE` to re-insert X, `wire_M_X`, and repoint both D and E at X. |
| `strip_missing_deltas(arc)` | Remove, from both directions of `arc`, deltas at sites missing at a tip endpoint (endpoints derived from the arc), sliding any that fall off the focus into `fitch_X_`/`deltas_ref_to_focus`.  Factors out the inline strip Round 5 performs after `merge_through`. |
| `init_component_picker()` | Compute `component_abort_threshold_` and reserve scratch vectors; called once at the start of `spr_refine`. |
| `pick_random_node_in_component(sink)` | Pick a random node in the connected component whose sink is `sink`, via DFS with rejection-sampling fallback (see "SPR throughput" section). |

### Core operations (pseudocode)

```
move_focus_updating_fitch_X(target):
    tree_.move_focus_to(target, pre_arc_hop = [&](arc a):
        for (s, {old, new}) in arcs[a].deltas:
            n_mismatches_ += fitch_X_.contains(s, old, old) - fitch_X_.contains(s, new, old)
            fitch_X_.on_ref_change(s, old, new))

eval_focal_arc(a):   // pre: origin(a) == focus
    savings = 0
    for (s, {a_st, b_st}) in arcs[a].deltas:
        if not fitch_X_.contains(s, a_st, a_st) and fitch_X_.contains(s, b_st, a_st):
            savings += 1
    return n_mismatches_ - savings

find_best_attachment_arc():   // identical structure to Round 5, but:
    best_cost = n_mismatches_                       // was ssize(focus_to_X_deltas_)
    ... priority-queue branch-and-bound, calling eval_focal_arc / move_focus_updating_fitch_X ...
    // pruning_threshold unchanged; returns {best arc (random among ties), best_cost}
```

### Initialization (pseudocode)

```
init_fitch_X_for_tip(X):                               // focus anywhere
    miss_X = tip_descs_[X].missations.intervals
    fitch_X_.clear(); n_mismatches_ = 0
    fitch_X_.uninformative_sites_ = miss_X
    // resolved_deltas_ = focus->X deltas at non-missing sites (== old init_focus_to_X_deltas)
    for sd in tip_descs_[X].seq_deltas:              // ref->X
        fitch_X_.resolved_deltas_[sd.site] = {sd.from, sd.to}
    for (s, d) in deltas_ref_to_focus:               // compose focus->ref
        if not miss_X.contains(s):
            push_front_site_deltas({s, d.to, d.from}, fitch_X_.resolved_deltas_)
    n_mismatches_ = ssize(fitch_X_.resolved_deltas_)

init_fitch_X_for_subtree(X):                           // focus at M
    fitch_X_.clear(); n_mismatches_ = 0
    D, E = X's two neighbors other than the focus M    // derived here, not passed in
    // Sites missing at BOTH children (only when both are tips) are size-4 uninformative.
    // Set them all up front -- including sites with no delta on any incident arc, which the
    // per-site walk below would never visit but which candidate arcs elsewhere may query.
    if D and E are both tips:
        fitch_X_.uninformative_sites_ = intersect(miss_D, miss_E)
    // Walk M-X, X-D, X-E to recover (m, d, e) at each site carrying a delta (see below).
    for each such site s, with m = M's state and d, e = D's, E's states:
        if fitch_X_.uninformative_sites_.contains(s): continue        // already size-4
        d_missing = D is a tip missing s;  e_missing = E is a tip missing s
        if d_missing or e_missing or d == e:             // size-1 (resolved)
            f = d_missing ? e : d                        // one-child-missing collapses to the other
            if m != f: fitch_X_.resolved_deltas_[s] = {m, f}; ++n_mismatches_
        else:                                            // size-2 (ambiguous)
            fitch_X_.ambiguous_masks_[s] = to_seq_letter(d) | to_seq_letter(e)
            if m not in {d, e}: ++n_mismatches_
```

The ordered `(d, e)` per ambiguous site is not stored here: after init we merge through X,
which composes X-D and X-E into a single D-E edge carrying exactly `{from: d, to: e}` at each
ambiguous site.  That edge is the ordered record reattachment consults.

The (m, d, e) walk: read `arcs[arc_MX].deltas` for m->x per site, then compose with
`arcs[arc_XD].deltas` (x->d) and `arcs[arc_XE].deltas` (x->e).  Sites absent from a given arc
inherit the upstream state.  Handle the missing-at-tip-child collapse exactly as the
conceptual "Fitch_X for subtrees" table prescribes.

### Reattachment (pseudocode)

Both attach methods share two helpers: `split_best_arc_inserting_M` (the `split_edge` on
`best_arc` plus `M_to_X_deltas_`/`m_overrides_` bookkeeping) and `wire_M_X` (adding the M-X
edge and orienting X).  The pseudocode below expands them inline for readability.

Tip -- identical to Round 5 `attach_tip_to_focal_arc`, with the membership tests rephrased
through `contains` (this is the "Reduction to Round 5" already described):

The `place_on_A` helper (shared by tip and subtree) pops the site from `M_to_X_deltas_`
only for resolved sites.  Ambiguous sites have no entry in `resolved_deltas_` (and hence in
`M_to_X_deltas_`), and their M-X cost is handled separately in subtree Step 2, so popping
must be skipped for them.

```
attach_tip_to_focal_arc(X, best_arc, M):             // pre: origin(best_arc) == focus
    M_to_X_deltas_ = fitch_X_.resolved_deltas_
    place_on_A(s, a, b):                             // M gets B's state b; shared with subtree Step 1
        if not fitch_X_.ambiguous_masks_.contains(s):
            pop_front({s,a,b}, M_to_X_deltas_)       // skip for ambiguous: no entry in resolved_deltas_
        m_overrides_[s] = b
    split_edge(best_arc, M, [&](sd={s,a,b}, A, B):
        if fitch_X_.uninformative_sites_.contains(s): return random(A,B)
        if fitch_X_.contains(s, b, a):  place_on_A(s,a,b); return A   // M gets b
        if fitch_X_.contains(s, a, a):  return B                      // M gets a
        side = random(A,B); if side==A: place_on_A(s,a,b); return side)
    wire M-X from M_to_X_deltas_; nodes[X].arc_to_focus = mate(arc_MX)
```

Subtree -- insert M into best_arc, re-insert X into the D-E edge (`arc_DE`, created by
merging through X at detach), then wire M-X:

```
attach_subtree_to_focal_arc(X, best_arc, M, arc_DE):      // pre: origin(best_arc) == focus
    // Step 1: split best_arc to insert M (identical callback to the tip case, using the
    //         same place_on_A helper).  Builds M_to_X_deltas_ = resolved_deltas_ adjusted on
    //         each A-side placement; records M's chosen state at candidate-delta sites in m_overrides_.
    M_to_X_deltas_ = fitch_X_.resolved_deltas_
    split_edge(best_arc, M, <same callback as tip, using place_on_A>)

    // Step 2: split arc_DE to re-insert X, choosing X's state per ambiguous site.
    //   m_state(s) = m_overrides_[s] if present else focus_state(s)   // = M's state at s
    //   (focus_state(s) = deltas_ref_to_focus[s].to if present else ref_sequence[s])
    // orig, dst = arc_DE's endpoints (whatever order merge_through produced); d = state at
    // orig, e = state at dst.  Per the split_edge contract, returning an endpoint keeps the
    // delta on that side and gives X the *opposite* endpoint's state.
    split_edge(arc_DE, X, [&](sd={s, d, e}, orig, dst):
        m = m_state(s)
        if m == e:      return orig                    // x = e, no M-X mutation (delta {d->e} on orig-X)
        else if m == d: return dst                     // x = d, no M-X mutation (delta {d->e} on X-dst)
        else:           M_to_X_deltas_[s] = {m, d}; return dst)  // m outside {d,e}: x = d, add M-X mutation

    // Step 3: wire M-X (deltas = M_to_X_deltas_) and orient the X-side toward the focus.
    arc_MX = add_arc(M, X)
    write arcs[arc_MX] / mate from M_to_X_deltas_
    nodes[X].arc_to_focus = mate(arc_MX)               // X -> M -> focus
    CHECK(nodes[D].arc_to_focus == k_no_arc || nodes[E].arc_to_focus == k_no_arc)
    nodes[D].arc_to_focus = find_arc(D, X)             // D -> X
    nodes[E].arc_to_focus = find_arc(E, X)             // E -> X
```

Step 3 unconditionally repoints both D and E toward X.  After split_edge(arc_DE, X), one
child's arc_to_focus is already pointing toward X (set by split_edge's arc_to_focus update),
but the other child was the X-side sink (arc_to_focus == k_no_arc) and needs repointing.
Rather than tracking which child was the sink, both are repointed — the already-correct one
is simply written again with the same value.  The CHECK asserts the precondition: exactly
one of D, E must be the X-side sink (k_no_arc) at this point.

`m_overrides_` is a small scratch map filled during Step 1 (only candidate-delta sites, where
M departs from the focus state); ambiguous sites rarely coincide with candidate deltas, so
`m_state(s)` is almost always just `focus_state(s)`.  M-X deltas come out exactly right: at resolved
sites `x` is the Fitch singleton and `M_to_X_deltas_` already tracks M vs it; at ambiguous
sites Step 2 adds a correction only when M lies outside `{d, e}`.

### `Utree` new method: `remove_edge` (in `core/utree.h` + `.cpp`)

Remove the arc pair between two nodes, clear the matching slot in each node's `arcs`, and
free the pair.  Both endpoints drop by one degree.  If an endpoint's `arc_to_focus` pointed
through the removed arc, it is cleared to `k_no_arc` (making that node a temporary sink);
otherwise `arc_to_focus` is left unchanged.  Used to cut the M-X edge before merging through
M and X.

```cpp
auto Utree::remove_edge(Node_index u, Node_index v) -> void {
  auto arc_uv = find_arc(u, v);
  CHECK_NE(arc_uv, k_no_arc);
  auto arc_vu = mate(arc_uv);
  for (auto& a : nodes[u].arcs) { if (a == arc_uv) { a = k_no_arc; break; } }
  if (nodes[u].arc_to_focus == arc_uv) { nodes[u].arc_to_focus = k_no_arc; }
  for (auto& a : nodes[v].arcs) { if (a == arc_vu) { a = k_no_arc; break; } }
  if (nodes[v].arc_to_focus == arc_vu) { nodes[v].arc_to_focus = k_no_arc; }
  free_arc_pair(arc_uv);
}
```

### `Utree::merge_through`: add a focal-link precondition

`merge_through(M)` already requires `focus != M` and `degree(M) == 2`.  Add one more
precondition: **M must not be a sink** -- its `arc_to_focus` must point at one of its two
arcs (`M.arc_to_focus != k_no_arc`).  Merging through a sink is ill-defined for focus
routing: both neighbors point their `arc_to_focus` at M, so the repair step would redirect
both onto the new edge and leave a transient 2-cycle instead of a clean sink.  The CHECK
turns that into a loud failure.

```cpp
CHECK_NE(nodes[M].arc_to_focus, k_no_arc);   // M must have a valid focal link (not a sink)
```

All existing callers already satisfy this (they merge an interior node after moving the
focus to a neighbor).  The subtree-SPR detach calls `remove_edge(M, X)`, which transiently
makes X a sink (clearing X's `arc_to_focus`); it immediately re-roots the X-side to give X
a valid focal link toward D before calling `merge_through(X)`.

### Helper: `strip_missing_deltas`

`merge_through` composes the two edges it replaces; if an endpoint is a tip, the composed
edge can inherit deltas at sites missing at that tip (from the inner edge).  These violate
the tip-adjacency invariant and must be removed.  Round 5 does this inline after each
`merge_through`; here we factor it into a `Utree_builder` method (it needs `tip_descs_` and
`fitch_X_`) reused by all three merge sites.  It takes just the arc and derives the endpoints
from it.

Removing such a delta is **not** a plain erase -- it can carry real information.  Picture
sliding the mutation along the edge toward the missing endpoint until it "falls off" that
end.  If the missing endpoint is the one *farther* from the focus, the slide is harmless: the
missing tip just inherits the focus-side state, which `deltas_ref_to_focus` already anchors.
But if the missing endpoint **is the focus**, the slide changes the focus's own state at that
site (from `d.from` to `d.to`) -- one step of a focus move -- so we must record it: update
`n_mismatches_`, call `fitch_X_.on_ref_change`, and extend `deltas_ref_to_focus`.  Skipping
this silently drops the far endpoint's real data at that site.  (Pre: `fitch_X_` is valid, so
`strip` runs only after the relevant `init_fitch_X_for_*`.)

```
strip_missing_deltas(arc):                             // endpoints derived from arc
    for node in {origin(arc), target(arc)}:
        if not tree_.is_tip(node): continue
        miss = tip_descs_[node].missations.intervals
        arc_from_node = (origin(arc) == node) ? arc : mate(arc)   // d reads {from: node, to: far}
        for (s, d) in arcs[arc_from_node].deltas where miss.contains(s):
            if node == focus:                          // mutation slides off the focus
                n_mismatches_ += fitch_X_.contains(s, d.from, d.from) - fitch_X_.contains(s, d.to, d.from)
                fitch_X_.on_ref_change(s, d.from, d.to)
                push_back_site_deltas({s, d.from, d.to}, deltas_ref_to_focus)
            erase s from arcs[arc_from_node].deltas and its mate
```

### Debug helpers: `utree_to_newick` / `utree_render` (in `core/utree.h` + `.cpp`)

Two free functions for inspecting a `Utree` while debugging the SPR moves, temporarily rooted
at a given node (default: the focus).  Both drive off the existing `annotated_arc_euler_tour`
(no recursion); tips are labeled `t<i>`, inner nodes `i<i>`, and each branch is annotated
with its parent->child mutations (e.g. `C3G,A7T`):

- `utree_to_newick(tree, root = focus) -> string`: a Newick string with per-branch mutations
  in NHX-style `[&m=...]` comments.
- `utree_render(tree, root = focus) -> string`: an indented ASCII-art tree (one node per
  line, `+--`/```--``` connectors, `[muts]` per branch) for eyeballing topology + mutations.

These proved essential for diagnosing the missing-data strip bug and are worth keeping as
general Utree debugging infrastructure.

### Outer loop: unified `spr_refine`

```
spr_refine(tree, tip_descs, bitgen, progress_hook):
    N = tree.num_tips
    if N <= 2: return
    builder = Utree_builder(move(tree), tip_descs, bitgen); t = builder.tree_
    max_attempts = 30*N; cur_deltas = t.count_deltas(); consecutive_non_improvements = 0

    for attempt in [0, max_attempts):
        // Pick a random internal edge M-X (M has degree 3).
        // M is rejection-sampled until degree 3; non-internal picks do NOT consume an attempt.
        do: M = t.pick_random_node(bitgen); while t.degree(M) != 3
        arc_MX = random arc from M; X = t.target(arc_MX)

        if t.is_tip(X):
            // ---- TIP SPR (Round 5 path, expressed via the unified core) ----
            // X is a tip: no D-E subtree, and fitch_X_ can be built any time from the
            // tip descriptor, so (unlike the subtree branch) we detach FIRST and init
            // afterwards, at the settled focus.
            P, Q = M's other two neighbors           // Round 5 calls these A, B
            d_MX,d_MP,d_MQ = arc delta counts
            if t.focus == X: t.move_focus_to(M)
            t.detach_tip(X)                          // X deg 0 (isolated); M deg 2 (P,Q)
            if t.focus == M: t.move_focus_to(P)      // move off M (plain: fitch_X_ not built yet)
            arc_PQ = t.merge_through(M)              // M freed (reused later)
            // Init before strip: strip may slide a mutation off the focus (an arc_PQ
            // endpoint), which needs fitch_X_ valid.
            t.move_focus_to(t.origin(arc_PQ))        // plain move: fitch_X_ not built yet
            builder.init_fitch_X_for_tip(X)          // focus-relative; X is disconnected
            strip_missing_deltas(arc_PQ)
            d_PQ = |arcs[arc_PQ].deltas|
            old_cost = d_MX + d_MP + d_MQ - d_PQ
            // rollback eval at P-Q
            best_arc = arc_PQ; best_cost = builder.eval_focal_arc(arc_PQ)
            if best_cost >= old_cost:
                do S = t.pick_random_node(bitgen) while S == X  // X-side is one node; any other is M-side
                builder.move_focus_updating_fitch_X(S)
                (found_arc, found_cost) = builder.find_best_attachment_arc()
                if found_cost < best_cost: best_arc, best_cost = found_arc, found_cost
            builder.move_focus_updating_fitch_X(t.origin(best_arc))
            builder.attach_tip_to_focal_arc(X, best_arc, M)
        else:
            // ---- SUBTREE SPR ----
            if t.degree(X) != 3: continue            // (X internal => degree 3; guard anyway)
            P, Q = M's other two neighbors; D, E = X's other two neighbors  // loop needs these
            d_MX,d_MP,d_MQ,d_XD,d_XE = arc delta counts
            t.move_focus_to(M)
            builder.init_fitch_X_for_subtree(X)        // derives its own D, E internally
            // Detach: cut M-X.  Focus (M) is the M-side sink.  Re-root the X-side onto D so
            // X keeps a valid focal link into merge_through (see Notes).
            t.remove_edge(M, X)                        // M deg 2 (P,Q); X deg 2 (D,E)
            // remove_edge clears X's arc_to_focus (it pointed at the now-gone arc_XM),
            // making X a temporary X-side sink.  Re-root: D becomes the X-side sink, X->D.
            t.nodes[D].arc_to_focus = k_no_arc         // D becomes the X-side sink
            t.nodes[X].arc_to_focus = arc_XD           // X -> D
            // Merge both connectors out
            arc_DE = t.merge_through(X)                // focus at M != X; X has focal link X->D; X freed (reused later)
            strip_missing_deltas(arc_DE)
            builder.move_focus_updating_fitch_X(P)     // move off M, updating fitch_X_
            arc_PQ = t.merge_through(M)                // focus at P != M; M freed (reused later)
            strip_missing_deltas(arc_PQ)
            d_PQ = |arcs[arc_PQ].deltas|; d_DE = |arcs[arc_DE].deltas|
            old_cost = d_MX + d_MP + d_MQ + d_XD + d_XE - d_PQ - d_DE
            // rollback eval at P-Q; settle the focus at P (origin of arc_PQ) so
            // eval_focal_arc's origin==focus precondition holds without relying on
            // merge_through's arc orientation
            builder.move_focus_updating_fitch_X(P)
            best_arc = arc_PQ; best_cost = builder.eval_focal_arc(arc_PQ)
            if best_cost >= old_cost:
                S = builder.pick_random_node_in_component(P)  // P is now the M-side sink
                builder.move_focus_updating_fitch_X(S)
                (found_arc, found_cost) = builder.find_best_attachment_arc()
                if found_cost < best_cost: best_arc, best_cost = found_arc, found_cost
            builder.move_focus_updating_fitch_X(t.origin(best_arc))
            builder.attach_subtree_to_focal_arc(X, best_arc, M, arc_DE)

        delta_change = best_cost - old_cost
        cur_deltas += delta_change
        consecutive_non_improvements = (delta_change < 0) ? 0 : consecutive_non_improvements + 1
        if consecutive_non_improvements >= N: break
        progress_hook(attempt+1, max_attempts, cur_deltas)

    tree = builder.finish()
```

Notes:
- **M-side start (`builder.pick_random_node_in_component`)**: the search must start on
  the M-side.  After both `merge_through` calls, P is the focus and the M-side sink.
  `pick_random_node_in_component(P)` picks a random node in P's component using a
  DFS-with-early-abort hybrid; see "SPR throughput: pick_random_node_in_component
  optimization" below for the algorithm history, the abort-threshold derivation, and the
  detailed final algorithm.  For tip SPR the X-side is a single node, so any node other
  than X is on the M-side and no walk is needed.
- **Focus bookkeeping around detach/merge** mirrors Round 5 (`focus != X` before
  `remove_edge`/`detach_tip`; `focus != M` before `merge_through`).  For the subtree branch,
  `init_fitch_X_for_subtree` runs with focus at M; `remove_edge(M, X)` clears X's
  `arc_to_focus` (it pointed at the now-gone M-X arc), making X a temporary X-side sink;
  we immediately re-root the X-side by
  making D the X-side sink and pointing X -> D, merge X out (focus still at M, off the
  X-side), move the focus to P *through the updating move* so `fitch_X_` stays consistent,
  then merge M out.
- **Init timing differs by branch** (intentional): the subtree branch must build `fitch_X_`
  at M *before* detaching, because the X-D and X-E arcs it reads disappear once X is merged
  through -- so every focus move afterward uses the *updating* `move_focus_updating_fitch_X`.
  The tip branch builds `fitch_X_` from the tip descriptor, which survives detachment, so it
  detaches and merges first (using *plain* `move_focus_to`, since `fitch_X_` is not yet
  built) and calls `init_fitch_X_for_tip` at the settled focus; the moves *after* init use
  the updating variant.  Crucially, `init_fitch_X_for_tip` must precede `strip_missing_deltas`
  in this branch, because the strip may slide a mutation off the focus and update `fitch_X_`
  (the subtree branch already inits before it strips).
- **`old_cost` / `delta_change`** use the formulas from the "Delta tracking" section; verify
  against `count_deltas()` in debug tests (see Verification).

### `add_tip` rework

```
add_tip(X):
    if tips_added_ == 0:                                // first tip: no tree yet
        add_first_tip(X)                                // sets focus=X, globally_missing_sites, deltas_ref_to_focus
    else:
        update_globally_missing_sites(X)
        init_fitch_X_for_tip(X)                         // focus wherever it is
        (best_arc, best_cost) = find_best_attachment_arc()
        if best_arc == k_no_arc: attach_tip_directly_to_isolated_focus(X)
        else:
            move_focus_updating_fitch_X(t.origin(best_arc))
            attach_tip_to_focal_arc(X, best_arc, alloc_inner_node())
```

`add_first_tip` sets the focus to X, initialises `globally_missing_sites` from X's missing
intervals, and seeds `deltas_ref_to_focus` from X's seq_deltas — no search, no edges.
`attach_tip_directly_to_isolated_focus` (reached when `tips_added_ == 1`) reads
`fitch_X_.resolved_deltas_` (formerly `focus_to_X_deltas_`).  The guide-tree and
refined-tree builders are otherwise unchanged.

### Integration into the build pipeline

```
build_initial_phylo_tree:
    1. Guide tree
    2. Nearest-first refinement (up to 5 rounds)
    3. SPR refinement (tips + subtrees), up to 30*N attempts   <-- spr_refine
    4. Regression rooting
    5. Conversion to Phylo_tree
```

`spr_refine` replaces the Round 5 `spr_refine_tips` call.  The hook parameter stays
`spr_refine_progress_hook` with the same `(attempts, max_attempts, cur_deltas)` signature
threaded through `build_rough_initial_tree_from_maple` in `cmdline.h`/`cmdline.cpp`, with the
~10%-increment stderr formatting hook at the call site.  `spr_refine_tips` is removed outright
(not kept as a wrapper): the function is renamed to `spr_refine`, and its tests are renamed to
`spr_refine_*` -- tip attempts exercise the same unified core.

### Files to modify

| File | Action | What |
|------|--------|------|
| `core/utree.h` | Modify | Add `Utree::remove_edge`; declare `utree_to_newick`/`utree_render` debug helpers; rename/declare `spr_refine` (replacing `spr_refine_tips`); keep `build_initial_phylo_tree` signature (hook already present from Round 5). |
| `core/utree.cpp` | Modify | Add `Relative_fitch_sets`; rework `Utree_builder` (new members, `init_fitch_X_for_tip`, `init_fitch_X_for_subtree`, `move_focus_updating_fitch_X`, generalized `eval_focal_arc`/`find_best_attachment_arc`, `split_best_arc_inserting_M`, `wire_M_X`, `attach_tip_to_focal_arc`, `attach_subtree_to_focal_arc`, `strip_missing_deltas`); add the focal-link `CHECK` to `merge_through`; remove `focus_to_X_deltas_` machinery; implement `remove_edge`, `utree_to_newick`/`utree_render`, and `spr_refine`; update `add_tip`; pipeline call. |
| `core/cmdline.h` / `.cpp` | Modify | Keep the SPR progress hook wired (rename to match `spr_refine` if desired). |
| `tests/utree_tests.cpp` | Modify | Tests for `remove_edge`, subtree reattachment, and `spr_refine` (see Verification for the full list). |

### Verification

1. **`Relative_fitch_sets`:** implemented as a `.cpp`-internal struct with no header exposure,
   so it is not unit-tested directly.  Its `contains`/`on_ref_change` behaviour across all
   site categories (resolved/ambiguous/uninformative) is exercised indirectly by the tip- and
   subtree-SPR tests below, which drive many focus moves and reattachments.
2. **`remove_edge` tests:** arc pair freed, both endpoints degree-2, `assert_utree_integrity`.
3. **`init_fitch_X_for_tip`:** the `init_focus_to_X_deltas` it subsumes has been removed, so
   there is no reference implementation left to diff against.  Its correctness is guarded
   end-to-end instead: `add_tip` uses it in the guide-tree / refined-tree builders, whose
   `finish()` runs `assert_utree_matches_tip_descs`, and the tip-SPR tests exercise it directly.
4. **Tip-SPR regression:** the Round 5 tip tests, renamed from `spr_refine_tips_*` to
   `spr_refine_*` (`spr_refine_integrity`, `_reduces_deltas`, `_monotonic`, `_small_trees`),
   still pass -- tip attempts exercise the same unified core -- and their total-delta
   bookkeeping matches `count_deltas()`.
5. **Subtree SPR:** several tests cover the subtree path, each asserting
   `assert_utree_integrity` + `assert_utree_matches_tip_descs` (which together verify X-side
   focus routing after reattachment and that the `merge_through` focal-link `CHECK` never
   fires during a valid move):
   - `spr_refine_relocates_subtree`: a deliberately misplaced subtree that `spr_refine`
     relocates, lowering total deltas.
   - `spr_refine_ambiguous_subtree_site`: a subtree whose two children disagree at a site
     (an ambiguous Fitch set), exercising `ambiguous_masks_` and Step 2.
   - `spr_refine_missing_data_subtree`: the missing-data collapse + strip -- a child tip
     missing a site the other child has data at, i.e. the case the strip's focus-slide fixes.
   - `spr_refine_random_stress`: larger randomized trees (many tips, several mutations each,
     some missing intervals), asserting integrity + tip-descriptor match after each of many
     `spr_refine` passes.  This is the test that surfaced the missing-data bug.
6. **Delta accounting:** the running `cur_deltas` tally (which drives the abort logic and the
   reported score) is not asserted against `count_deltas()` per move.  Tree corruption from a
   bad move is instead caught by the per-pass `assert_utree_matches_tip_descs` in the subtree
   / stress tests, and net improvement by `spr_refine_reduces_deltas`.  (A debug
   `cur_deltas == count_deltas()` assert after each move would additionally guard the
   `old_cost`/`delta_change` arithmetic and is cheap to add.)
7. **Convergence:** a second `spr_refine` call yields no further improvement.
8. **End-to-end:** `build_initial_phylo_tree` produces valid `Phylo_tree`s
   (`assert_phylo_tree_integrity`, `assert_phylo_tree_matches_tip_descs`).
9. **Build/run:** `cmake --build build/debug --target tests -j$(nproc) &&
   build/debug/tests/tests --gtest_filter='Utree*'`; then compare total delta counts on the
   standard demo datasets (expecting gains where tip SPR fell short of the MCMC mean: Ebola,
   Zika, RSV-A).

## Initial results

Total parsimony score (mutations / deltas) at each tree-construction stage on the standard
Delphy demo datasets plus the deep RSV-A phylogeny.  Lower is better.  The columns match the
Round 5 table (`2026-06-10-02-better-tree-init-round5-spr-refinement.md`), whose "SPR refined"
column measured tip-only SPR; here that stage runs the unified tip + subtree refinement of
this round.  `max_attempts` was also raised from `10*N` (Round 5) to `30*N`, which the RSV-A
dataset in particular benefits from.  Values vary somewhat run-to-run.

| Dataset                          |    N    |    L    | Guide tree | Refined tree (rounds) | SPR refined (attempts)  | Approx. MCMC mean |
|----------------------------------|---------|---------|------------|-----------------------|-------------------------|-------------------|
| Ebola (Dudas et al, 2017) (*)    |   1,610 |  18,996 |      3,161 |             2,881 (1) |           2,790 ( 5.9N) |             2,786 |
| Ebola (Gire et al, 2014)         |      81 |  18,959 |         55 |                55 (0) |              55 ( 1.0N) |                55 |
| H3N2 (Rambaut et al, 2008)       |     165 |   1,698 |        402 |               392 (4) |             392 ( 1.0N) |               394 |
| H5N1 (Andersen lab, 2025) (*)    |   3,339 |  13,136 |      4,392 |             4,361 (3) |           4,341 ( 4.0N) |             4,347 |
| mpox (O'Toole et al, 2023)       |      41 | 197,209 |        345 |               345 (0) |             344 ( 1.3N) |               343 |
| mpox (Parker et al, 2025)        |     177 | 197,209 |        950 |               944 (2) |             943 ( 1.5N) |               944 |
| SARS-CoV-2 (Lemieux et al, 2021) |     757 |  29,903 |        664 |               664 (0) |             664 ( 1.0N) |               667 |
| Zika (Metsky et al, 2017)        |     174 |  10,807 |      1,246 |             1,246 (0) |           1,213 ( 3.9N) |             1,210 |
| RSV-A (NCBI, 2012--2025) (**)    |   3,000 |  15,277 |     45,627 |            37,149 (3) |          35,334 (20.9N) |               N/A |

(*) Note: Observed substantial variability run-to-run, but the refined tree is always
    substantially better, and the many refinement rounds produce more stable refined scores.
(**) Note: The RSV-A dataset is contaminated by 3 nearby RSV-B sequences that are ~2,500
mutations from the RSV-A sequences.

Adding subtree moves (vs. Round 5's tip-only SPR) improves the SPR stage across the board:

- **Ebola (Dudas):** 2,971 -> 2,790, essentially reaching the MCMC mean (2,786).
- **RSV-A:** 36,813 -> 35,334, and still improving at ~60% of the 30N budget -- the case that
  motivated raising `max_attempts` to `30*N`.
- **Zika:** 1,224 -> 1,213, approaching the MCMC mean (1,210).
- **H3N2 (392), H5N1 (4,341), SARS-CoV-2 (664)** now match or dip just below the MCMC mean
  (which averages over non-parsimonious sampled trees).
- Datasets already at their tip-SPR optimum (Ebola Gire, both mpox) are essentially unchanged.


## SPR throughput: pick_random_node_in_component optimization

The `exp_100000` simulated dataset (N=100,000 tips, exponential growth;
[delphy-2026-paper-data](https://github.com/broadinstitute/delphy-2026-paper-data/)
`sims/exp_100000/inputs/exp_100000.maple`) was used as a throughput benchmark for
`spr_refine`.  Profiling revealed that `pick_random_node_in_component` dominated subtree SPR time.
Four algorithm generations were evaluated; timings are for subtree attempts only,
measured over 2 minutes on a 22-core machine.

| Algorithm | pick_S avg (us) | Sub attempts / 2 min | Notes |
|-----------|----------------|----------------------|-------|
| Rejection sampling | 122,730 | 1,109 | Original.  O(depth) per try; catastrophic when X-side is large (most random picks land on X-side and are rejected). |
| Coroutine walk + reservoir sampling | 6,748 | 19,150 | Walk M-side via `annotated_arc_euler_tour`; reservoir-sample uniformly.  ~18x faster, but coroutine `co_yield` overhead (~20 ns/yield x 200K yields per call) still dominates. |
| Direct iterative DFS + uniform pick | 1,315 | 71,988 | Replace coroutine with explicit DFS stack; collect M-side into a scratch vector, one RNG pick.  ~5x faster.  Still O(N) walk every call. |
| Direct DFS + 5% abort + rejection fallback | 80 | 133,334 | Abort if M-side > 5% of total nodes; fall back to rejection (acceptance rate >= 5%).  ~16x faster. |
| Direct DFS + optimal abort + rejection fallback | 25 | ~133,334 | Abort threshold = sqrt(N_total * log2(N_total)) (~1% of N).  ~3x faster than 5%. |

After this optimization `pick_S` (25 us) is negligible relative to the search phase (332 us).
Throughput is now limited by `find_best_attachment_arc`.


### Final algorithm: `builder.pick_random_node_in_component(sink)`

Three members of `Utree_builder` support this method, initialized once by
`init_component_picker()` at the start of `spr_refine`:

- `component_abort_threshold_` -- computed from tree size; see below
- `component_dfs_stack_`       -- scratch DFS stack, reserved to `abort_threshold + 1`
- `component_nodes_`           -- scratch node list, reserved to `abort_threshold + 1`

**Abort threshold.**  Walking T nodes costs c*T (one DFS node visit each).  Rejection
sampling when the M-side has more than T nodes has acceptance rate >= T/N_total; each
try costs approximately c*log2(N_total) (one arc_to_focus walk of ~log2(N) hops).
Expected rejection cost: c * log2(N_total) * (N_total/T).  Setting the two costs equal:

```
  T^2 = N_total * log2(N_total)
  T   = sqrt(N_total * log2(N_total))
```

The constant c cancels because a DFS node visit and an arc_to_focus hop are the same
operation (pointer chase + bounds check).  For a full tree with N tips,
N_total = 2*N - 1 and T is about 1% of N_total.

**Pre-condition:** `arc_to_focus[sink] == k_no_arc` (sink is the component's focus/sink).

**Algorithm:**

```
pick_random_node_in_component(sink):
    CHECK arc_to_focus[sink] == k_no_arc
    component_nodes_.clear(); component_dfs_stack_.clear()
    for each arc a from sink: push target(a) onto component_dfs_stack_
    while component_dfs_stack_ not empty:
        V = pop_back(component_dfs_stack_)
        push_back V onto component_nodes_
        if |component_nodes_| > component_abort_threshold_:
            // Component is large: fall back to rejection (acceptance rate >= T/N_total).
            loop:
                S = pick_random_node(bitgen_)
                cur = S
                while arc_to_focus[cur] != k_no_arc: cur = target(arc_to_focus[cur])
                if cur == sink: return S
        for each arc a from V where a != k_no_arc and a != arc_to_focus[V]:
            push target(a) onto component_dfs_stack_
    CHECK component_nodes_ not empty
    return component_nodes_[uniform random index in [0, |component_nodes_|)]
```

**Correctness.**  At each node V, skipping `arc_to_focus[V]` avoids going back toward sink,
so the DFS visits exactly the nodes reachable from sink's neighbors.  The rejection path
walks any randomly chosen node to its sink and accepts iff it matches — nodes in other
components reach a different sink and are rejected automatically, with no need to enumerate
excluded nodes explicitly.  The entry CHECK ensures sink is the focus, which is required for
`arc_to_focus` to identify components correctly.
