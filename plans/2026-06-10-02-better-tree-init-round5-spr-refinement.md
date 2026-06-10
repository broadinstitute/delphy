# Better Tree Init, Round 5: SPR Refinement of Tips

This is the fifth in a series of plans implementing a new tree initialization pipeline for
Delphy.  See the [Round 1 plan](2026-05-04-01-better-tree-init-round1-utree-and-guide-tree.md),
section "Roadmap", for the full multi-round plan.

## Background

Rounds 1--3 build a `Utree` via greedy parsimony insertion: a rough guide tree (Round 1),
refined by reinserting tips in nearest-first order (Round 3), iterated until convergence.
Round 4 roots the tree via regression.

The nearest-first refinement is powerful --- it typically reduces total deltas by 10--25% on
deep-phylogeny datasets --- but it rebuilds the tree from scratch each round.  Each insertion
is greedy: once a tip is placed, it stays put even if later insertions reveal that an earlier
choice was suboptimal.  As the tree fills in, the greedy choices become increasingly
constrained by the existing topology, and late-arriving tips have little room to correct
early mistakes.

## Goal for Round 5

Add a post-construction refinement step that directly corrects individual tip placements via
SPR (Subtree Prune and Regraft) moves on the unrooted `Utree`.  For each tip, detach it from
the tree, search for the best reattachment edge starting from a random node, and reattach
it --- rolling back to the original attachment point if the search does not find a better
position.

This fills a different niche from Round 3: Round 3 corrects the insertion *order* (each tip
starts close to its neighbors), while Round 5 corrects individual *placements* after the
full topology is known.  The two complement each other: Round 3 builds a good approximate
tree, Round 5 fine-tunes it.

## Algorithm overview

Each refinement attempt picks a random tip and tries to improve its placement:

1. **Detach X:** Remove X's edge to its parent inner node M (M becomes degree-2), then merge
   M away by composing its two remaining edges into one.
2. **Evaluate the rollback cost:** Compute the cost of reattaching X at the merged A-B edge
   (the original attachment region).
3. **Search from a random start:** Move the focus to a random node S and run the same
   priority-queue branch-and-bound search as the builder's `find_best_attachment_arc`.
4. **Accept or roll back:** If the search found a position with cost < rollback cost, accept
   it (reattach at the best edge).  Otherwise, roll back (reattach at the A-B edge).

Starting the search from a random node rather than from the old attachment point is
deliberate: the old position is very likely already a local optimum, so the branch-and-bound
pruning would prevent exploration of distant alternatives that happen to lie on the other
side of an improbable intermediate mutational configuration.  The random start lets the
search discover improvements in other parts of the tree.

### Iteration structure

Up to `10*N` random attempts.  Each attempt picks a random tip, detaches it, searches, and
accepts or rolls back.  Track total tree deltas: if no attempt improves the total for `N`
consecutive attempts, stop early.

```
spr_refine_tips(in_out_tree, tip_descs, bitgen, progress_hook):
    N = in_out_tree.num_tips
    if N <= 2: return

    builder = Utree_builder(std::move(in_out_tree), tip_descs, bitgen)
    tree = builder.tree_    // alias
    max_attempts = 10 * N
    cur_deltas = tree.count_deltas()
    consecutive_non_improvements = 0

    for attempt = 0 to max_attempts - 1:
        X = tree.pick_random_tip(bitgen)

        // Record old state
        arc_XM = X's only arc
        M = target(arc_XM)
        d_XM = ssize(arcs[arc_XM].deltas)
        A, B = M's other two neighbors
        d_AM = ssize(A-M edge deltas)
        d_MB = ssize(M-B edge deltas)

        // Move focus away from X if needed, then detach
        if tree.focus == X:
            tree.move_focus_to(M)
        tree.detach_tip(X)

        // Move focus away from M if needed, then merge
        if tree.focus == M:
            tree.move_focus_to(A)
        arc_AB = tree.merge_through(M)

        // Strip deltas at sites missing at tip endpoints (see note below)
        for node in {A, B}:
            if tree.is_tip(node):
                strip deltas at node's missing sites from arc_AB and its mate

        d_AB = ssize(arcs[arc_AB].deltas)
        old_cost = d_XM + d_AM + d_MB - d_AB

        // Evaluate rollback cost at the original position (A-B edge)
        tree.move_focus_to(origin(arc_AB))
        builder.init_focus_to_X_deltas(X)
        best_arc = arc_AB
        best_cost = builder.eval_focal_arc(arc_AB)

        // Search from a random start, unless rollback already improves on old_cost
        if best_cost >= old_cost:
            do { S = tree.pick_random_node(bitgen) } while (S == X)
            builder.move_focus_updating_focus_to_X_deltas(S, X)
            auto [found_arc, found_cost] = builder.find_best_attachment_arc(X)
            if found_cost < best_cost:
                best_arc = found_arc
                best_cost = found_cost

        // Reattach at best position (rollback or search result)
        builder.move_focus_updating_focus_to_X_deltas(origin(best_arc), X)
        builder.attach_tip_to_focal_arc(X, best_arc, M)
        delta_change = best_cost - old_cost
        cur_deltas += delta_change
        if delta_change < 0:
            consecutive_non_improvements = 0
        else:
            consecutive_non_improvements += 1

        if consecutive_non_improvements >= N:
            break

        progress_hook(attempt + 1, max_attempts, cur_deltas)

    in_out_tree = builder.finish()
```

Note: even when rolling back, the total can decrease because the merged A-B edge may cancel
mutations that were needlessly split across A-M and M-B.

**Missing-site delta stripping:** `merge_through` composes the A-M and M-B deltas into a
single A-B edge.  If A or B is a tip, the composed edge can inherit deltas from the inner
edge (A-M or M-B) at sites where the tip has no data.  These deltas carry no information and
violate the invariant that tip edges have no deltas at missing sites.  The loop strips them
immediately after merging.  The cost tracking (`d_AB`, `old_cost`, `delta_change`) uses the
post-strip delta count, so everything stays consistent.

### Edge cases

- **N <= 2:** No SPR possible (detaching a tip from a 2-tip tree leaves an isolated node).
  `spr_refine_tips` returns immediately.
- **N = 3:** One inner node M.  Detaching any tip leaves a 2-tip tree with one edge.  The
  search finds only that edge and reattaches.  This always produces the same tree (there is
  only one unrooted topology on 3 tips).

## Adapting `Utree_builder` for SPR refinement

The builder (`Utree_builder`) already contains all the branch-and-bound search and tip
attachment machinery needed for SPR refinement.  Rather than extracting a separate class,
adapt `Utree_builder` to support both incremental construction and tip refinement on an
existing tree.

`Utree_builder` remains an implementation-detail class defined in `core/utree.cpp` (not
exposed in the header).

### Changes to `Utree_builder`

1. **Second constructor from an existing tree.**  Takes an already-built `Utree` by move,
   sets `tips_added_` to `num_tips` (all tips already placed), and initializes `L_` and
   `sqrt_6L_` from `ref_sequence`.  No progress hook is needed (SPR refinement passes its
   own hook to `spr_refine_tips`).

2. **`find_best_attachment_arc` returns `std::pair<Arc_index, int>`** (best arc, best cost)
   instead of returning just the arc.

3. **`attach_tip_to_focal_arc` takes an explicit `Node_index M`** instead of allocating M
   internally.  During SPR refinement, the caller reuses the inner node freed by
   `merge_through`.

4. **New `alloc_inner_node()` helper** on `Utree_builder` that returns
   `tree_.num_tips + tree_.num_inner_nodes_so_far++`.  Used by `add_tip` to allocate M
   before calling `attach_tip_to_focal_arc`.

5. **Make all members public.**  `Utree_builder` is an implementation-detail class in the
   `.cpp` file, so there is no encapsulation boundary to protect.  Drop the `private:` section
   and make everything `public:`.  This lets the free function `spr_refine_tips` access
   `init_focus_to_X_deltas`, `eval_focal_arc`, `find_best_attachment_arc`,
   `move_focus_updating_focus_to_X_deltas`, `attach_tip_to_focal_arc`, and `tree_` directly.

### Updated `add_tip` flow

```
add_tip(X):
    update_globally_missing_sites(X)
    init_focus_to_X_deltas(X)
    auto [best_arc, best_cost] = find_best_attachment_arc(X)
    if best_arc == k_no_arc:
        attach_tip_directly_to_isolated_focus(X)
    else:
        move_focus_updating_focus_to_X_deltas(origin(best_arc), X)
        attach_tip_to_focal_arc(X, best_arc, alloc_inner_node())
```

### Public free function

The public API for SPR refinement is a free function that constructs a `Utree_builder` from
the existing tree and uses its search/attach methods directly to implement the refinement
loop (the pseudocode above).  The `Utree_builder` provides the tools; the free function
orchestrates them.

```cpp
auto spr_refine_tips(
    Utree& tree, const std::vector<Tip_desc>& tip_descs,
    absl::BitGenRef bitgen,
    const std::function<void(int,int,int)>& progress_hook = [](int,int,int){}) -> void;
```

Internally, the function moves the tree into a builder, runs the refinement loop calling
`detach_tip`, `merge_through`, `init_focus_to_X_deltas`, `find_best_attachment_arc`,
`attach_tip_to_focal_arc`, etc., and moves the tree back out via `finish()` at the end.

## Core operations

### `Utree::merge_through(Node_index M) -> Arc_index`

Merge the two edges incident to a degree-2 inner node M into a single edge, removing M from
the tree.  This is the inverse of `split_edge`: where `split_edge` inserts a degree-2 node
by replacing one edge with two, `merge_through` removes a degree-2 node by replacing its two
edges with one.

**Preconditions:**
- `degree(M) == 2`
- `focus != M`

**Steps:**

1. Identify M's two arcs (`arc_MA`, `arc_MB`) and their targets (`A`, `B`).

2. Compose deltas: `A_to_B_deltas` = copy of `arcs[mate(arc_MA)].deltas` (the A->M deltas),
   then call `append_site_deltas(A_to_B_deltas, arcs[arc_MB].deltas)` (the M->B deltas).
   This produces the composed A->B deltas, with cancellation at sites where M was the only
   node differing from both A and B.

3. Build `B_to_A_deltas` by reversing each entry of `A_to_B_deltas`.

4. Allocate a new arc pair for A-B.  Set targets and deltas.

5. Update A's and B's arc slots: replace the old arcs (pointing to M) with the new ones
   (pointing to each other).

6. Update `arc_to_focus`: since `focus != M`, exactly one of A or B has `arc_to_focus`
   pointing through M (the one on the far side of M from the focus).  The other has
   `arc_to_focus` pointing elsewhere (toward the focus, away from M).  Update the affected
   node's `arc_to_focus` to the new arc.

7. Clear M's arc slots and `arc_to_focus`.  Free the two old arc pairs.

8. Return the new arc A->B.

### `Utree::detach_tip(Node_index X) -> Node_index`

Member function.  Remove the edge between tip X and its neighbor M, leaving M as a degree-2
node in the tree.  Returns M.  Does **not** call `merge_through` --- the caller is
responsible for managing the focus and calling `merge_through(M)` when ready.

**Preconditions:**
- `tree.is_tip(X)` and `tree.degree(X) == 1`
- `tree.focus != X` (X becomes disconnected, so it cannot be the focus)
- Tree has at least 3 tips

**Steps:**

1. Find X's sole arc `arc_XM`; let `M = target(arc_XM)`, `arc_MX = mate(arc_XM)`.

2. Remove X's arc from M's arc slots (set the slot containing `arc_MX` to `k_no_arc`).

3. Clear X's arc slots and `arc_to_focus`.  Free the X-M arc pair.

4. Return M.

**Note:** After `detach_tip`, the focus can be at M (M is still in the tree, just degree-2).
This lets the caller update any focus-associated state before moving the focus away and
calling `merge_through(M)`.

## Integration

### Pipeline change

Insert SPR refinement between the nearest-first rebuild loop and rooting:

```
build_initial_phylo_tree:
    1. Guide tree
    2. Nearest-first refinement (up to 5 rounds)
    3. SPR tip refinement (up to 10*N attempts)      <-- NEW
    4. Regression rooting
    5. Conversion to Phylo_tree
```

### Progress hook

Add an `spr_refine_progress_hook` parameter to `build_initial_phylo_tree`:

```cpp
auto build_initial_phylo_tree(
    Real_sequence ref_sequence, std::vector<Tip_desc> tip_descs,
    absl::BitGenRef bitgen,
    const std::function<void(int,int)>& guide_tree_progress_hook = [](int,int){},
    const std::function<void(int,int,int)>& refined_tree_progress_hook = [](int,int,int){},
    const std::function<void(int,int,int)>& spr_refine_progress_hook = [](int,int,int){},
    const std::function<void(const Rooting_info&)>& rooting_hook = [](const Rooting_info&){})
    -> Phylo_tree;
```

The `spr_refine_progress_hook` takes `(attempts_so_far, max_attempts, cur_deltas)`, where
`max_attempts = 10*N` and `cur_deltas` is the current total delta count.  Thread it through `build_rough_initial_tree_from_maple` in
`cmdline.h`/`cmdline.cpp` with the same defaulted parameter pattern as the other hooks.
At the call site in `cmdline.cpp`, supply a hook that prints progress to stderr at ~10%
increments (matching the existing pattern for other stages).

## Files to create/modify

| File | Action | What |
|------|--------|------|
| `core/utree.h` | Modify | Add `detach_tip` and `merge_through` methods to `Utree`; declare `spr_refine_tips`; update `build_initial_phylo_tree` signature |
| `core/utree.cpp` | Modify | Implement `detach_tip` and `merge_through`; adapt `Utree_builder` (second constructor, `attach_tip_to_focal_arc` takes M, `find_best_attachment_arc` returns pair); implement free `spr_refine_tips`; update `build_initial_phylo_tree` |
| `core/cmdline.h` | Modify | Thread `spr_refine_progress_hook` through `build_rough_initial_tree_from_maple` |
| `core/cmdline.cpp` | Modify | Forward `spr_refine_progress_hook`; supply formatting hook at call site |
| `tests/utree_tests.cpp` | Modify | Add tests for `merge_through`, `detach_tip`, and `spr_refine_tips` |

## Verification

1. **`merge_through` tests:**
   - 3-tip tree: merge the inner node, verify the result is a single edge with composed
     deltas.  Check that cancelling deltas (A and B agree at a site where M had a different
     state) reduce the delta count.
   - Verify `assert_utree_integrity` passes after merge.
   - Test with focus on each side of M (both directions of `arc_to_focus` update).

2. **`detach_tip` tests:**
   - 3-tip tree: detach each tip in turn, verify M is returned and has degree 2.
   - 4-tip tree: detach a tip, verify M has degree 2 and the rest of the tree is intact.
   - Verify focus at M is allowed (and that `merge_through` works after moving focus away).

3. **`Utree_builder` refactoring tests:**
   - Verify `Utree_builder` still produces correct trees after the changes (existing tests
     exercise this — they should pass unchanged).

4. **`spr_refine_tips` tests:**
   - Construct a deliberately misbuilt tree (e.g., 4 tips where one is attached at the wrong
     edge).  Verify that SPR refinement moves it to the correct position and reduces total
     deltas.
   - Run on the output of `build_refined_tree` and verify `assert_utree_integrity` and
     `assert_utree_matches_tip_descs` pass after refinement.
   - Verify convergence: after refinement, a second call to `spr_refine_tips` should
     produce no further improvement (total deltas unchanged).

5. **End-to-end:** Verify `build_initial_phylo_tree` still produces valid trees
   (`assert_phylo_tree_integrity` and `assert_phylo_tree_matches_tip_descs` pass).

6. **Build:** `cmake --build build/debug --target tests -j$(nproc) && build/debug/tests/tests --gtest_filter='Utree*'`

7. **Manual test:** Run `delphy_ui` with `--v0-init mp-plus-timing` against standard demo
   datasets.  Compare total delta counts before and after SPR refinement.

## Initial results

Total parsimony score (mutations or deltas) for each tree construction method on the
standard Delphy demo datasets and one deep-phylogeny dataset.  Lower is better.  Similar
numbers appear in the Round 3 plan (`2026-05-15-01-better-tree-init-round3-guide-tree-reordering.md`);
the values here differ slightly owing to run-to-run variability and the missing-site delta
bugfix (`2026-06-10-01-fix-spurious-missing-site-deltas.md`).

| Dataset                          |    N    |    L    | Guide tree | Refined tree (rounds) | SPR refined (attempts) | Approx. MCMC mean |
|----------------------------------|---------|---------|------------|-----------------------|------------------------|-------------------|
| Ebola (Dudas et al, 2017) (*)    |   1,610 |  18,996 |      3,147 |             3,072 (2) |           2,971 (5.3N) |             2,786 |
| Ebola (Gire et al, 2014)         |      81 |  18,959 |         55 |                55 (0) |              55 (1.0N) |                55 |
| H3N2 (Rambaut et al, 2008)       |     165 |   1,698 |        401 |               393 (2) |             393 (1.0N) |               394 |
| H5N1 (Andersen lab, 2025) (*)    |   3,339 |  13,136 |      4,619 |             4,357 (1) |           4,351 (1.8N) |             4,347 |
| mpox (O'Toole et al, 2023)       |      41 | 197,209 |        343 |               343 (0) |             343 (1.0N) |               343 |
| mpox (Parker et al, 2025)        |     177 | 197,209 |        947 |               945 (1) |             944 (1.1N) |               944 |
| SARS-CoV-2 (Lemieux et al, 2021) |     757 |  29,903 |        664 |               664 (0) |             664 (1.0N) |               667 |
| Zika (Metsky et al, 2017)        |     174 |  10,807 |      1,247 |             1,247 (0) |           1,224 (5.8N) |             1,210 |
| RSV-A (NCBI, 2012--2025) (**)    |   3,000 |  15,277 |     47,137 |            37,678 (5) |         36,813 (10.0N) |               N/A |

(*) Note: Observed substantial variability run-to-run, but refined tree is always substantially better,
    and the many refinement rounds produce more stable refined tree scores
(**) Note: The RSV-A dataset is contaminated by 3 nearby RSV-B sequences that are ~2,500
mutations from the RSV-A sequences.
