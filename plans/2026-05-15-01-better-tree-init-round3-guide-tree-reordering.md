# Better Tree Init, Round 3: Guide Tree Reordering

This is the third in a series of plans implementing a new tree initialization pipeline for
Delphy.  See the [Round 1 plan](2026-05-04-01-better-tree-init-round1-utree-and-guide-tree.md)
section "Roadmap" for the full multi-round plan.

## Background

Round 1 built a **guide tree** --- a `Utree` built by inserting tips in input order
(0, 1, ..., N-1) with greedy parsimony placement.  The guide tree captures approximate
relationships --- nearby tips are likely genuinely similar --- but its topology may be poor
because the input order is arbitrary.

## Goal for Round 3

Use the guide tree to determine a better insertion order, then rebuild a **refined tree**
(also a `Utree`) from scratch.  The key idea: traverse the guide tree starting from a random
tip, always adding the nearest unvisited tip next.  Since the guide tree's topology captures
approximate relationships, this order ensures that each new tip is close to something already
in the tree, dramatically reducing the chance of early misplacements.

## Algorithm

For each directed arc I->J in the guide tree, compute `arc_nearest[I->J]`: the nearest tip
reachable from I in the direction of J, and its distance from I.  This is done in two DFS
passes over the guide tree (Passes 1 and 2), after which a priority-queue-driven traversal
(Pass 3) determines the insertion order and builds the refined tree.

### Worked example

Consider a 4-tip guide tree with tips A, B, C, D and inner nodes X, Y (edge sizes in
brackets):

```
    A --[2]-- X --[1]-- Y --[3]-- C
              |         |
             [4]       [1]
              |         |
              B         D
```

Root the DFS at R = A:

```
       R = A
           |
       [A->X, 2]
           |
           X
          / \
         /   \
 [X->B, 4]   [X->Y, 1]
       /       \
      B         Y
               / \
       [Y->C, 3] [Y->D, 1]
             /     \
            C       D
```

The Euler tour from R visits 10 arcs (E = entering = away from R, L = leaving = toward R):

```
  #1  A->X  E     #2  X->B  E     #3  B->X  L     #4  X->Y  E     #5  Y->C  E
  #6  C->Y  L     #7  Y->D  E     #8  D->Y  L     #9  Y->X  L     #10 X->A  L
```

**Pass 1** processes leaving events (post-order) and annotates the mate arc (pointing away
from R):

```
  #3  leaving B->X  ->  annotate X->B = {B, 4}             (tip: dist = edge size)
  #6  leaving C->Y  ->  annotate Y->C = {C, 3}             (tip)
  #8  leaving D->Y  ->  annotate Y->D = {D, 1}             (tip)
  #9  leaving Y->X  ->  annotate X->Y = {D, 1+min(3,1)=2}  uses Y->C (#6), Y->D (#8) OK
  #10 leaving X->A  ->  annotate A->X = {D, 2+min(4,2)=4}  uses X->B (#3), X->Y (#9) OK
```

**Pass 2** runs a second DFS from R and processes entering events (pre-order), annotating
the mate arc (pointing toward R):

```
  #1  entering A->X  ->  annotate X->A = {A, 2}             (A is a tip, dist to A = 0)
  #2  entering X->B  ->  annotate B->X = {D, 4+min(2,2)=6}  uses X->Y (P1#9), X->A (P2#1) OK
  #4  entering X->Y  ->  annotate Y->X = {A, 1+min(4,2)=3}  uses X->B (P1#3), X->A (P2#1) OK
  #5  entering Y->C  ->  annotate C->Y = {D, 3+min(1,3)=4}  uses Y->D (P1#8), Y->X (P2#4) OK
  #7  entering Y->D  ->  annotate D->Y = {C, 1+min(3,3)=4}  uses Y->C (P1#6), Y->X (P2#4) OK
```

Note that events P2#2 and P2#7 involve an arbitrary tie-break (both options have the same
distance).

At each step, the annotations we need have already been computed: Pass 1 annotations are all
available before Pass 2 starts, and within Pass 2 each node's "toward R" annotation is set
before its children's (pre-order guarantees parent before child).

**Pass 3** (starting at tip D).  In each entry below: T = tip being added, I = arc origin,
H = closest previously added tip to I, N = current node during walk.  `d_X_Y` denotes the
guide-tree distance from X to Y.

```
  Start: add D
    Push D->Y (dist 4, nearest C, closest_prev_tip D at dist 0).

  Pop D->Y (T=C, H=D, dist 4).
    Move refined tree focus to D, then add C.
    Walk from D toward C (d_T_I=4, d_I_H=0):
      At Y (d_N_I=1): push Y->X (dist 3, nearest A).
        d_N_H=1+0=1 vs d_N_T=4-1=3 -> H closer -> closest_prev_tip D at dist 1.
      Arrive at C.

  Pop Y->X (T=A, H=D, dist 3, d_I_H=1).
    Move refined tree focus to D, then add A.
    Walk from Y toward A (d_T_I=3, d_I_H=1):
      At X (d_N_I=1): push X->B (dist 4, nearest B).
        d_N_H=1+1=2 vs d_N_T=3-1=2 -> tie -> say T, closest_prev_tip A at dist 2.
      Arrive at A.

  Pop X->B (T=B, H=A, dist 4, d_I_H=2).
    Move refined tree focus to A, then add B.
    Arc target is B -> done.

  Queue empty -> done.  Insertion order: [D, C, A, B].
```

### Detailed steps

#### Pass 1: Annotate guide tree arcs (post-order)

Root the guide tree at an arbitrary node R (e.g., tip 0).  Run
`annotated_arc_euler_tour(guide_tree, R)`.  On each leaving event for arc X->P (where X is
the node being left and P is its parent toward R), annotate the mate arc P->X:

- If X is a tip: `arc_nearest[P->X] = {X, ssize(P->X.deltas)}`.
- If X is inner with children C_1, ..., C_k (the arcs from X other than X->P):
  ```
  i* = argmin_i arc_nearest[X->C_i].dist
  arc_nearest[P->X] = {arc_nearest[X->C_i*].tip, ssize(P->X.deltas) + arc_nearest[X->C_i*].dist}
  ```
  Each `arc_nearest[X->C_i]` was set on an earlier leaving event (post-order guarantees
  children before parent).

#### Pass 2: Annotate guide tree arcs (pre-order)

Run a second `annotated_arc_euler_tour(guide_tree, R)`.  On each entering event for arc P->X
(where P is the parent toward R and X is the node being entered), annotate the mate arc X->P:

- Compute `min_from_P_not_via_X`: the minimum `arc_nearest[a].dist` over all of P's outgoing
  arcs `a` other than P->X.  All such arcs are already annotated: arcs pointing away from R
  were set in Pass 1; arcs pointing toward R were set earlier in this pass (pre-order
  guarantees parent before child).
  - Special case: if P is a tip (degree 1, so it has no arcs other than P->X): set
    `min_from_P_not_via_X = 0` and `best_tip = P` (P itself is the closest tip not via X).
    In practice, this only happens when P = R and R is a tip.
- `arc_nearest[X->P] = {best_tip, ssize(X->P.deltas) + min_from_P_not_via_X}`, where
  `best_tip` is the tip from the minimizing arc (or P itself for the special case).

#### Pass 3: Nearest-first traversal and refined tree construction

Interleave the insertion order computation with the refined tree construction: as each tip is
identified by the guide tree traversal, immediately add it to the refined tree via the
existing `Utree_builder`.

Alongside the insertion order, we also track `closest_prev_tip` for each tip: the nearest
already-added tip in the guide tree.  When building the refined tree, we use this to position
the builder's focus near the insertion point before each `add_tip` call, so the branch-and-bound
search starts close to where the new tip will land.

Each priority queue entry stores `(dist, arc, closest_prev_tip, d_closest_prev_tip)`:
- `dist` = `arc_nearest[arc].dist` (priority key).
- `arc` = the guide tree arc.
- `closest_prev_tip` = the nearest already-added tip to the arc's origin.
- `d_closest_prev_tip` = the guide-tree distance from the arc's origin to `closest_prev_tip`.

1. Pick a random starting tip S.  Add it to the refined tree (`builder.add_tip(S)`).  Push
   S's single outgoing arc with `closest_prev_tip = S, d_closest_prev_tip = 0`.

2. Loop while the queue is non-empty:
   a. Pop the entry with the smallest `dist`.  Let I be the arc's origin, H be the entry's
      `closest_prev_tip`, `d_I_H` be the entry's `d_closest_prev_tip`, and T be
      `arc_nearest[arc].tip` (the nearest tip to add).  Also let `d_T_I = dist` (the
      guide-tree distance from T to I).
   b. Move the refined tree's focus to H, then add T to the refined tree
      (`builder.add_tip(T)`).
   c. Walk from I toward T along the guide tree, starting with the popped arc I->J.  At I
      itself, do nothing (I is already visited and its branching arcs were pushed in an
      earlier iteration).  At each subsequent node N on the path (tracking cumulative distance
      `d_N_I` from N to I), push all of N's outgoing arcs except the arrival arc and the
      continuation arc (the one with `arc_nearest[arc].tip == T`).  For each pushed arc,
      compute that new entry's `closest_prev_tip` and `d_closest_prev_tip`:
      ```
      d_N_H = d_N_I + d_I_H      (distance from N to H, via I)
      d_N_T = d_T_I - d_N_I      (distance from N to T, the just-added tip)
      ```
      If `d_N_T <= d_N_H`: push with `closest_prev_tip = T, d_closest_prev_tip = d_N_T`.
      Otherwise: push with `closest_prev_tip = H, d_closest_prev_tip = d_N_H`.

**Why `arc_nearest.tip == T` correctly identifies the continuation arc:** if
`arc_nearest[I->J].tip == T`, then T is the nearest tip from I via J.  At J, exactly one of
J's subtrees contains T, and T must be the nearest tip in that subtree (otherwise a closer
tip would have been chosen for `arc_nearest[I->J]`).  By induction, the same holds at each
intermediate node, so `arc_nearest[arc].tip == T` uniquely identifies the correct arc at
every step.

**No visited flags needed.**  The guide tree is a tree, so there is exactly one path between
any two nodes.  The set of walked edges forms a connected subtree rooted at the starting tip,
and every arc pushed onto the priority queue leads into a completely disjoint, unwalked
subtree --- so it is impossible to push an arc that was already walked or already in the queue.
The only arcs to exclude at each node are the arrival and continuation arcs, which is purely
local bookkeeping.

**`closest_prev_tip` correctness.**  When an arc is pushed from node N, the only already-added
tips reachable from N are (a) tips in the visited subtree, reachable via the arrival arc, and
(b) T, reachable via the continuation arc.  The nearest already-added tip via the visited
subtree is H at distance `d_N_I + d_I_H` (since H was the nearest to I, and I is `d_N_I`
away from N).  Since `d_N_I` grows and `d_N_T` shrinks along the walk, there is a single
crossover: before it, H is closer; after it, T is closer.

Moreover, `closest_prev_tip` **does not become stale** between push and pop: no tip added
after the push can be closer to N than the recorded value.  Intuitively, any such tip T'
must branch off the I->T walk path, so d(I,T') >= d(I,T) (otherwise T' would have been the
walk's target instead of T).  Since the extra distance from I to N is shared, T remains at
least as close to N as T'.  The same argument applies recursively for tips added from earlier
walks, bottoming out at the starting tip where the closest_prev_tip is exact.
See Appendix A for a detailed proof (not yet reviewed for correctness by a human).

## Data structures

```cpp
struct Arc_nearest {
  Node_index tip = k_no_node;  // Nearest tip reachable in this arc's direction
  int dist = 0;                // Distance from arc origin to that tip (sum of ssize(deltas))
};
```

Stored per arc in a `std::vector<Arc_nearest>` indexed by `Arc_index`, sized to match the
guide tree's arc vector.

## Integration

Change `build_initial_phylo_tree` to accept separate progress hooks for each stage and insert
a refinement step between guide tree construction and rooting:

```cpp
auto build_initial_phylo_tree(
    Real_sequence ref_sequence, std::vector<Tip_desc> tip_descs,
    absl::BitGenRef bitgen,
    const std::function<void(int,int)>& guide_tree_progress_hook = [](int,int){},
    const std::function<void(int,int,int)>& refined_tree_progress_hook = [](int,int,int){})
    -> Phylo_tree;
```

The `refined_tree_progress_hook` takes `(round, tips_so_far, total_tips)`, where `round`
starts at 1.  `build_initial_phylo_tree` wraps this into a `(tips_so_far, total_tips)` hook
that it passes to `build_refined_tree` (which doesn't know about rounds).

```
1. utree = build_guide_tree(std::move(ref_sequence), tip_descs, bitgen,
                            guide_tree_progress_hook)
2. Iteratively refine (up to 5 rounds):
     refined = build_refined_tree(utree, tip_descs, bitgen,
         [&](tips_so_far, total_tips) { refined_tree_progress_hook(round, ...); })
     if refined.count_deltas() >= utree.count_deltas(): stop (no improvement)
     utree = std::move(refined)
3. root = midpoint_root_utree(utree, tip_descs)
4. rate = estimate_rate_and_root_date(utree, root, tip_descs)
5. return utree_to_phylo_tree(utree, root, tip_descs, rate, bitgen)
```

Refinement is applied iteratively because each round can further improve the topology.  The
loop stops when a round fails to reduce the parsimony score (total delta count), or after 5
rounds, whichever comes first.  If a refinement round produces a tree with *more* deltas than
the input, the refinement is abandoned and the input tree is kept.

`build_refined_tree` internally annotates the guide tree (Passes 1-2), creates a
`Utree_builder`, and runs the nearest-first traversal (Pass 3) calling `builder.add_tip` for
each tip as it is identified.

The caller in `cmdline.cpp` passes distinct hooks for the guide tree and refined tree phases
via `build_rough_initial_tree_from_maple` (which gained a second hook parameter).  Later
rounds of this roadmap will add more hooks (e.g., root-to-tip regression, SPR refinement).

## Functions

All new functions are declared in `core/utree.h` and defined in `core/utree.cpp`.
`Arc_nearest` is a file-local struct in `core/utree.cpp` (implementation detail).

### `Utree` helpers (new)

- `count_deltas()`: total site deltas across all edges (each undirected edge counted once,
  via `count_arc_deltas` on even-indexed arcs).  Used by the iterative refinement loop.
- `pick_random_tip(bitgen)`: returns a uniformly random tip index from `[0, num_tips)`.
  Used by `for_each_tip_in_nearest_first_order` to choose the starting tip.

### `for_each_tip_in_nearest_first_order`

```cpp
auto for_each_tip_in_nearest_first_order(
    const Utree& guide_tree,
    absl::BitGenRef bitgen,
    const std::function<void(Node_index tip, Node_index closest_prev_tip)>& callback) -> void;
```

Annotates the guide tree's arcs (Passes 1-2), picks a random starting tip, and traverses
the guide tree in nearest-first order (Pass 3), calling `callback(tip, closest_prev_tip)` for each
tip.  `closest_prev_tip` is the nearest already-added tip in the guide tree (`k_no_node` for the
first tip).  Declared in the header for testability: tests can materialize the insertion
order by collecting into a vector.

### `build_refined_tree`

```cpp
auto build_refined_tree(
    const Utree& guide_tree,
    const std::vector<Tip_desc>& tip_descs,
    absl::BitGenRef bitgen,
    const std::function<void(int,int)>& progress_hook = [](int,int){})
    -> Utree;
```

Builds a refined `Utree` by adding tips in the order determined by `for_each_tip_in_nearest_first_order`.
Copies `guide_tree.ref_sequence` internally to construct the `Utree_builder`, calls
`for_each_tip_in_nearest_first_order` with a callback that moves the builder's focus to
`closest_prev_tip` (if not `k_no_node`) before invoking `builder.add_tip`, then returns the
finished tree.

**`Utree_builder` change:** add a public `move_focus_to(Node_index target)` method that
calls `tree_.move_focus_to(target)`, allowing `build_refined_tree` to position the builder's
focus before each `add_tip` call.  This method is only valid *between* `add_tip` calls --- it
repositions the tree's focus without maintaining the internal `focus_to_X_deltas_` state that
`add_tip` uses (that state is managed internally by `add_tip` itself via
`init_focus_to_X_deltas` and `move_focus_updating_focus_to_X_deltas`).

## Files to modify

| File | Action | What |
|------|--------|------|
| `core/utree.h` | Modify | Add `count_deltas()` and `pick_random_tip()` to `Utree`; update `build_initial_phylo_tree` to take two progress hooks (second with round number); declare `for_each_tip_in_nearest_first_order` and `build_refined_tree` |
| `core/utree.cpp` | Modify | Add public `move_focus_to` to `Utree_builder` (between-`add_tip` focus hint only); define `Arc_nearest` (file-local); implement `for_each_tip_in_nearest_first_order` and `build_refined_tree` (copies `guide_tree.ref_sequence` internally); update `build_initial_phylo_tree` with iterative refinement loop |
| `core/cmdline.h` | Modify | Add `refined_tree_progress_hook` parameter to `build_rough_initial_tree_from_maple` |
| `core/cmdline.cpp` | Modify | Pass distinct guide tree and refined tree hooks to `build_initial_phylo_tree` |
| `tests/utree_tests.cpp` | Modify | Add tests for `for_each_tip_in_nearest_first_order` and `build_refined_tree` |

## Verification

1. **`for_each_tip_in_nearest_first_order` tests**:
   - 1-tip tree: callback called once with tip 0.
   - 2-tip tree: callback called with both tips.
   - 3-tip tree with known topology: verify the order is consistent with nearest-first
     traversal (e.g., if the starting tip is 0 and tip 1 is closer to tip 0 than tip 2 in
     the guide tree, then tip 1 comes before tip 2).
   - Verify the callback is called exactly N times and covers all tips.

2. **`build_refined_tree` tests**: Verify the resulting Utree passes `assert_utree_integrity`
   and `assert_utree_matches_tip_descs` across multiple random seeds.  (The refined tree
   usually has fewer deltas than the guide tree, but not always --- the iterative refinement
   loop in `build_initial_phylo_tree` handles this by abandoning rounds that don't improve.)

3. **End-to-end**: Verify `build_initial_phylo_tree` still produces valid trees (the existing
   end-to-end tests cover this since `build_initial_phylo_tree` now calls `build_refined_tree`
   internally).

4. **Build**: `cmake --build build/debug --target tests -j$(nproc) && build/debug/tests/tests --gtest_filter='Utree*'`

5. **Manual end-to-end test**: Run `delphy_ui` with `--v0-init mp-plus-timing` against
   standard demo datasets.  Spot-check that the results are sensible.

## Preliminary results

Total parsimony score (mutations or deltas) for each tree construction method on the
standard Delphy demo datasets and one deep-phylogeny dataset.  Lower is better.  Results
can vary between runs owing to randomness.

| Dataset                          |    N    |    L    | Guide tree | Refined tree (rounds) | Approx. MCMC mean |
|----------------------------------|---------|---------|------------|-----------------------|-------------------|
| Ebola (Dudas et al, 2017) (*)    |   1,610 |  18,996 |      3,294 |             2,814 (5) |             2,786 |
| Ebola (Gire et al, 2014)         |      81 |  18,959 |         55 |                55 (2) |                55 |
| H3N2 (Rambaut et al, 2008)       |     165 |   1,698 |        402 |               394 (5) |               394 |
| H5N1 (Andersen lab, 2025) (*)    |   3,339 |  13,136 |      5,933 |             4,392 (5) |             4,347 |
| mpox (O'Toole et al, 2023)       |      41 | 197,209 |        345 |               347 (1) |               343 |
| mpox (Parker et al, 2025)        |     177 | 197,209 |        952 |               946 (2) |               944 |
| SARS-CoV-2 (Lemieux et al, 2021) |     757 |  29,903 |        665 |               665 (1) |               667 |
| Zika (Metsky et al, 2017)        |     174 |  10,807 |      1,296 |             1,300 (1) |             1,210 |
| RSV-A (NCBI, 2012--2025) (**)    |   3,000 |  15,277 |     47,934 |            38,540 (3) |               N/A |

(*) Note: Observed substantial variability run-to-run, but refined tree is always substantially better,
    and the many refinement rounds produce more stable refined tree scores
(**) Note: The RSV-A dataset is contaminated by 3 nearby RSV-B sequences that are ~2,500
mutations from the RSV-A sequences.

---

## Appendix A: Proof that `closest_prev_tip` cannot become stale

**Note:** This proof was generated by Claude and has not been reviewed for correctness by a
human.

**Claim:** When an arc pushed from node N is later popped, its recorded `closest_prev_tip` is
still the nearest previously-added tip to N --- no tip added between the push and pop is closer.

**Setup.** An arc is pushed from node N during a walk from origin I toward target T.  At push
time, `closest_prev_tip` = min(d(N,H), d(N,T)), where H is the nearest previously-added tip
via the arrival direction and T is the just-added continuation tip.  We must show that no tip
T' added after the push satisfies d(N,T') < min(d(N,H), d(N,T)).

Since N is an interior node of the walk (degree 3 in a binary tree), it has exactly three
outgoing arcs: arrival (toward I), continuation (toward T), and branching (the arc being
pushed).  In a tree, the path from N to any other node passes through exactly one of these
three arcs.  We consider each case for T'.

**Case 1: T' reachable via N's branching arc.**  The branching arc is the arc being pushed
itself.  T' would be in this arc's subtree.  T' is only added when this arc is popped, at
which point the `closest_prev_tip` is consumed (used to position the builder's focus), not
compared against.  No staleness.

**Case 2: T' reachable via N's continuation arc (T' branches off the N->T path).**  T' is in a
subtree hanging off some node M between N and T on the walk.  The path from I to T' goes
through M: `d(I,T') = d(I,M) + d(M,T')`.  The path from I to T also goes through M:
`d(I,T) = d(I,M) + d(M,T)`.

The walk targeted T because T was the nearest tip from I in this direction (determined by
`arc_nearest` from the popped entry).  Since both T and T' branch off from M on the same side
of I, `arc_nearest` selected T over T', so `d(I,T) <= d(I,T')`.  This gives
`d(M,T) <= d(M,T')`.

For staleness we need d(N,T') < d(N,T).  Since M is between N and T:
```
d(N,T') = d(N,M) + d(M,T') >= d(N,M) + d(M,T) = d(N,T)
```
**Contradiction.**

**Case 3: T' reachable via N's arrival arc (T' is in the visited subtree's direction).**
The path from N to T' goes back through I: `d(N,T') = d(N,I) + d(I,T')`.  The recorded
closest_prev_tip via the arrival direction is H at distance `d(N,H) = d(N,I) + d(I,H)`.
Staleness requires `d(I,T') < d(I,H)`.

T' was added by popping an arc pushed from some node Q.  Q's arc leads into the subtree
containing T'.  There are two sub-cases:

**Case 3a: Q is on the same walk (between I and N).**  Q is between I and N, so
`d(I,Q) + d(Q,N) = d(I,N)`.  The path from I to T' goes through Q:
`d(I,T') = d(I,Q) + d(Q,T')`.  But Q's branching arc was pushed during the same walk that
targeted T.  From I's perspective, T was the nearest tip in this direction, and Q is on the
I->T path.  T' branches off at Q, so both T and T' are reachable from I through Q:
```
d(I,T)  = d(I,Q) + d(Q, ..., T)
d(I,T') = d(I,Q) + d(Q,T')
```
Since T was selected as nearest from I: `d(I,T) <= d(I,T')`, giving
`d(Q, ..., T) <= d(Q,T')`.

Now, for staleness we need d(I,T') < d(I,H).  But at the time N's arc was pushed,
`closest_prev_tip` was set to min(d(N,H), d(N,T)).  We need to show d(I,T') >= d(I,H),
i.e., T' is no closer to I than H.

H's distance from I: `d(I,H)` was the nearest previously-added tip to I at push time.  But
the walk chose T as the target because T was the nearest tip from I *in the walk direction*.
H might be in a different direction from I entirely (via I's own arrival arc from a previous
walk).  However, the walk from I toward T was initiated by popping an entry whose
`closest_prev_tip` already accounted for H.  The key: T' branches off the I->T path at Q, so
`d(I,T') >= d(I,T)`.  And H was chosen at push time with
`d(N,H) = d(N,I) + d(I,H)`, so staleness requires `d(I,T') < d(I,H)`.  Combined:

```
d(I,T) <= d(I,T') < d(I,H)
```

This would mean T is closer to I than H --- but H was the nearest previously-added tip to I
via the arrival direction, which is a *different* direction from the walk.  T is in the walk
direction.  These are different subtrees from I, so this inequality is not itself a
contradiction.  However, it reduces to exactly the same claim one level up: the entry that
was popped to initiate this walk had its own `closest_prev_tip`, and by the same argument
applied to *that* entry's walk, no intermediate tip could have been closer.

**Case 3b: Q is on a previous walk entirely.**  This reduces to Case 3a applied to the
previous walk.  The recursion terminates at the base case: the starting tip S, whose single
arc is pushed with `closest_prev_tip = S, d = 0`.  Since S is the only added tip at that
point, this is exact and cannot become stale.

**Inductive conclusion.**  At the base level, `closest_prev_tip` is exact.  Each subsequent
level's `closest_prev_tip` is either via the continuation direction (Case 2, proven
non-stale) or via the arrival direction (Case 3, which reduces to the previous level).  By
induction, no level can become stale.  QED
