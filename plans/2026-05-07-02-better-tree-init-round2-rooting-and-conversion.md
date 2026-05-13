# Better Tree Init, Round 2: Midpoint Rooting and Phylo_tree Conversion

This is the second in a series of plans implementing a new tree initialization pipeline for
Delphy, continuing from Round 1 (Utree data structure and guide tree construction). See
the [Round 1 plan](2026-05-04-01-better-tree-init-round1-utree-and-guide-tree.md), section
"Roadmap", for the full multi-round plan.

## Background

Round 1 built the `Utree` data structure and a `build_guide_tree` function that produces an
unrooted, bifurcating, mutation-annotated tree by greedy parsimony insertion. The guide tree
captures approximate relationships between tips but has no root and no timing information.

Delphy's MCMC operates on a `Phylo_tree` — a rooted, timed, explicit mutation-annotated tree
(EMAT). To use the guide tree as an MCMC starting point, we must root it and convert it to a
`Phylo_tree`. This plan implements that conversion, yielding a basic end-to-end pipeline:
input sequences → guide tree → midpoint root → `Phylo_tree`.

## Goal for Round 2

Implement midpoint rooting on the `Utree` and conversion from the rooted `Utree` to a
`Phylo_tree`. The result is a `build_initial_phylo_tree` function that takes the same inputs
as `build_usher_like_tree` and returns a `Phylo_tree` ready for MCMC.

Specifically:

1. **Timed midpoint rooting**: Find a diametral pair `(u, v)` — two tips at maximum
   distance by total site delta count (i.e., whose distance equals the diameter `D`).
   When all tip dates are the same, root at the midpoint of the diametral path `u → v`.
   When tip dates differ, use a rough mutation rate estimate to place the root at a
   position along the diametral path that accounts for the time difference between `u`
   and `v`.

2. **Rate estimation**: Simple OLS regression of root-to-tip mutation counts against tip
   dates, giving a rough mutation rate and root date estimate.

3. **Phylo_tree construction**: Orient the tree from the root, assign times to all nodes,
   convert site deltas to timed mutations, build missation maps, and validate.

Root-to-tip regression for *root selection* (trying many candidate roots and picking the one
with the best R²) is deferred to Round 4. Here we use regression only for *time estimation*
at the fixed midpoint root.

## Midpoint rooting

### Finding a diametral pair

The diameter `D` of an unrooted tree is the maximum distance (by total site delta count)
between any two nodes. For a tree with non-negative edge weights, any diametral pair
(two nodes at distance `D`) consists of tips. The standard two-pass algorithm finds a
diametral pair:

1. **Pass 1**: DFS from an arbitrary tip (tip 0). At each node, track cumulative distance
   (sum of `ssize(arc.deltas)` along the path from the start). The farthest node `u` is one
   endpoint of a diametral pair.

2. **Pass 2**: DFS from `u`. The farthest node `v` is the other endpoint. The distance
   `D = dist(u, v)` is the diameter.

Multiple tip pairs may achieve distance `D`, but no simple path has length greater than `D`.
Any diametral pair `(u, v)` is valid for the rooting step.

### Placing the root

**Intuition**: The MRCA of two tips `x` and `y` lies along the path `x → y`.  Under a
strict molecular clock with rate `lambda`, if both tips have the same date, the MRCA is
roughly at the midpoint: the `M` mutations on the path accumulated over two equal
branches, so each branch has length `M / (2 * lambda)` in time. When tips have different
dates `t_x` and `t_y`, the MRCA time is roughly `t_MRCA = (t_x + t_y) / 2 - M / (2 * lambda)`,
shifting the MRCA position toward the older tip.  The stochastic nature of mutations makes
these formulas true only in the mean.

The more mutations there are between `x` and `y`, the earlier their MRCA. The tree root is
the MRCA of the entire tree — equal to the earliest pairwise MRCA across all tip
pairs. Since earlier MRCAs come from longer paths, a diametral pair (the longest path) is
a natural candidate for the pair whose MRCA is the tree root. This isn't guaranteed — with
unequal tip dates, a non-diametral pair with a large date difference could have an earlier
MRCA — but it's a reasonable heuristic that avoids an O(N²) scan over all tip
pairs. (Round 4's root-to-tip regression over candidate roots will refine the root
choice.)

Assume a rough overall mutation rate `lambda` = 1 mutation per 30 days (a ballpark for most
common viruses; the exact value doesn't matter much since OLS regression refines the rate
estimate later). Let `t_u` and `t_v` be the tip dates (midpoints of `[t_min, t_max]`), and
`M = D` the number of deltas between `u` and `v`. The root time `t_R` satisfies:

```
lambda * [(t_u - t_R) + (t_v - t_R)] = M
t_R = (t_u + t_v) / 2 - M / (2 * lambda)
```

If `t_R >= min(t_u, t_v)` (the root would be younger than one of the endpoints), clamp to
`t_R = min(t_u, t_v) - 14` (a 2-week buffer to avoid zero-length branches).

The root is placed at fraction `c` along the `u → v` path:

```
c = (t_u - t_R) / [(t_u - t_R) + (t_v - t_R)]
```

This fraction determines which edge gets the root and how many of its deltas go to each
side. When `t_u = t_v`, `c = 1/2`, recovering standard midpoint rooting.

**Implementation using the focus mechanism**: Each node's `arc_to_focus` is its outgoing
arc toward the current focus node; following these links from any node traces a path to
the focus. After Pass 2, move the focus to `v` (`move_focus_to(v)`). Now `arc_to_focus`
links point toward `v`, so walking from `u` and following `arc_to_focus` at each node
traces the `u → v` diameter path. Accumulate edge sizes (`ssize(arc.deltas)`) along this
walk. The target cumulative distance from `u` is `n_u_total = round(c * D)`. When the
cumulative distance crosses `n_u_total`, the current edge is the root edge. Let
`cum_before` be the cumulative distance before this edge; the number of deltas from this
edge assigned to the `u`-side is `n_u = n_u_total - cum_before`, with the remaining deltas
going to the `v`-side.

### Rooting the Utree

After identifying the root edge and `n_u`, root the `Utree` in place:

1. Allocate the root node `R` at index `num_tips + num_inner_nodes_so_far` (the next
   available inner node slot). Increment `num_inner_nodes_so_far`.

2. Call `split_edge(root_arc, R, site_delta_side)`. The `site_delta_side` callback assigns
   deltas to the `u`-side or `v`-side to match the split count `n_u`. Since site deltas are
   unordered in the map, the specific assignment of individual deltas to each side is
   arbitrary — iterate and assign the first `n_u` to the `u`-side and the rest to the
   `v`-side.

The focus location is undefined after this call; `utree_to_phylo_tree` moves the focus to
the root when it needs it.

**Preallocation**: The current `Utree::make_empty` allocates `2N - 2` nodes for an unrooted
tree with N tips. After rooting, we have `2N - 1` nodes. Change `make_empty` to allocate
`max(1, 2 * num_tips - 1)` nodes so the root node slot is available. The existing 2 scratch
arc pairs suffice: `split_edge` consumes one net arc pair (allocates 2, frees 1), leaving
1 scratch pair unused.

**Integrity check**: `assert_utree_integrity` must allow inner nodes with degree 2 (not
just 3), since the root node created by `split_edge` has only two neighbors.

### Edge cases

- **N = 0**: No tree, no rooting needed.
- **N = 1**: Single tip, no edges. The root is the tip itself. Diameter is 0. The resulting
  `Phylo_tree` has 1 node (the tip is also the root).
- **N = 2**: One edge. Follows the general algorithm without special handling.

## Rate estimation

After rooting, compute root-to-tip mutation counts and estimate a rough mutation rate via
simple OLS regression.

### Root-to-tip mutation counts

DFS from the root R. At each node, the root-to-node distance is the cumulative
sum of `ssize(arc.deltas)` along the path from R. For tips, this gives `m_i` (root-to-tip
mutation count). This reuses `annotated_arc_euler_tour` from Round 1: on each `entering`
arc, increment the running distance by `ssize(arc.deltas)`; on each `leaving` arc, decrement
by the same amount.

### OLS regression

Fit the linear model `m_i = lambda * t_i + b` where `t_i` is the tip date (midpoint of
`[t_min, t_max]` for tips with date uncertainty) and `lambda` is the overall mutation rate
(mutations per day across the whole genome):

```
lambda = Cov(m, t) / Var(t)
b      = mean_m - lambda * mean_t
t_root = -b / lambda = mean_t - mean_m / lambda
```

This is a standard OLS slope/intercept computation — no matrix algebra needed.

**Numerical stability**: Tip dates can have large offsets (e.g., days since 2020-01-01),
so the naive formula `Var(t) = E[t^2] - mean_t^2` suffers from catastrophic cancellation.
Instead, use a two-pass approach:

1. **Pass 1** (loop over tips via `utree_tips`): compute `mean_t`.
2. **Pass 2** (DFS from root): accumulate `Var(t) = E[(t - mean_t)^2]` and
   `Cov(m, t) = E[m * (t - mean_t)]` using deviations `dt_i = t_i - mean_t`.
   Since `E[dt] = 0`, `E[m * dt] = Cov(m, t)` directly — no subtraction of large terms.

### Fallback for degenerate cases

The regression can fail if:
- **All tips have the same date** (`Var(t) = 0`): no temporal signal.
- **Negative slope** (molecular clock violation, e.g., from a poor root choice).
- **N = 1**: No regression possible.

In all degenerate cases, fall back to a heuristic: use `lambda_fallback = 1.0 / 30.0`
overall mutations per day (i.e., across the entire genome, not per-site). Set `t_root =
mean_t - mean_m / lambda_fallback` (or `t_tip` for N = 1 where `mean_m = 0`). This is the
same rough rate used for the timed midpoint rooting.

A future round (Round 4) will use regression to *select* the root, which largely avoids
these issues. For now, we just need finite, non-degenerate times.

## Conversion to Phylo_tree

### Overview

Create a `Phylo_tree` from the rooted `Utree`, the input `tip_descs`, the estimated
overall mutation rate `lambda`, and the estimated root date `t_root`. The conversion has
six steps: set the reference sequence and place root mutations (Step 1), build topology via
DFS (Step 2), bubble up missations (Step 3), randomize mutation times (Step 4), normalize
the reference sequence to the root (Step 5), and validate (Step 6).

### Step 1: Set the reference sequence and place root mutations

Move the focus to the root (`move_focus_to(root)`) so that `deltas_ref_to_focus` records the
ref→root state changes. Set the `Phylo_tree`'s `ref_sequence` to the **original**
`ref_sequence` (the one used to build `tip_descs`), not the root sequence. This is because
`tip_descs`' missation `from_states` are relative to the original `ref_sequence`. Place
mutations above the root from `deltas_ref_to_focus` (each `{site, from, to}` becomes a
`Mutation{from, site, to, t_root}`). The post-processing step
`rereference_to_root_sequence` (Step 5) will later normalize the `Phylo_tree`'s
`ref_sequence` to the root sequence, emptying root mutations and rewiring all missation maps.

### Step 2: Build topology, assign times, and place mutations

Create `Phylo_tree{2*N - 1}` with tips at `0..N-1` and inner nodes at `N..2N-2`.

**Initialize the root** R before the DFS: set `tree.root = R`,
`tree.at(R).parent = k_no_node`, `tree.at(R).t = t_root`. If R is a tip (N = 1), copy
its name, `t_min`, `t_max`, and missations from `tip_descs[R]`, and clamp
`tree.at(R).t = clamp(t_root, t_min, t_max)`. Otherwise (R is an inner node), set
`tree.at(R).t_min = -max_float`, `tree.at(R).t_max = +max_float`.

**Single DFS** from R (using `annotated_arc_euler_tour`), tracking the root-to-node
mutation count `m_X` (starting at 0 for the root):

**On entering arc P → X:**

1. **Set up topology**: Set `tree.at(X).parent = P` and add X to P's children.

2. **Copy metadata**: For tips, copy name, `t_min`, `t_max`, and missations from
   `tip_descs[X]`. For inner nodes, set `t_min = -max_float`, `t_max = +max_float`.

3. **Estimate time**: Update `m_X += ssize(arc.deltas)`. Then
   `t_X_est = t_root + m_X / lambda`.
   - For tips: clamp to date bounds: `t_X = clamp(t_X_est, t_min, t_max)`.
   - For inner nodes: use `t_X_est` directly (fixed up on exit).

4. **Place mutations**: The arc's site deltas from P to X directly give the branch mutations
   for X. Each site delta `{site, from, to}` on the entering arc records the state change
   from the parent's state to the child's state — exactly the mutation that belongs on the
   corresponding `Phylo_tree` branch.

   For each site delta, create `Mutation{from, site, to, t}` with `t` set to the child's
   time `t_X` as a placeholder. Actual mutation times are assigned later by
   `randomize_mutation_times` (Step 4), which samples each mutation's time uniformly in
   `(t_P, t_X]` and sorts them. This avoids the complication that inner node times are not
   yet finalized during this DFS (they get fixed up on exit).

**On leaving node X** (all children processed):

5. **Fix up inner node times** (inner nodes only): Enforce that X is strictly earlier than
   all its children:
   ```
   t_X = min(t_X, min(t_L, t_R) - min_branch_length)
   ```
   where `min_branch_length = max(0.5, 0.5 / lambda)` — half a day or half a mutation
   interval, whichever is larger. This prevents zero-length branches.

6. **Update `m_X`**: Decrement `m_X -= ssize(arc.deltas)` to restore it for the parent.

**After the DFS completes:**

7. **Fix up root time**: Apply the same constraint to R:
   ```
   t_R = min(t_R, min(t_child1, t_child2) - min_branch_length)
   ```
   The DFS does not yield a "leaving" event for the source node, so the root is not covered
   by step 5. This step handles the case where child fixups cascade and push R's children
   earlier than `t_root`.

### Step 3: Build missation maps

1. **Tips**: Missations were already copied from `tip_descs` in Step 2. Inner nodes start
   with empty missations.

2. **`fix_up_missations(tree)`**: Bubble up common missations from children to parents,
   reconstruct `from_states`, and remove mutations at missing sites. This is the same
   post-processing that `build_random_tree` and `build_usher_like_tree` use.

### Step 4: Randomize mutation times

Call `randomize_mutation_times(tree, bitgen)` to assign proper times to all mutations. This
samples each mutation's time uniformly in `(t_parent, t_child]` and sorts mutations on each
branch by time, matching the existing pattern used by `build_random_tree`.

### Step 5: Normalize reference sequence

Call `rereference_to_root_sequence(tree)`. This normalizes the `Phylo_tree`'s `ref_sequence`
to match the root node's actual sequence (removing all mutations above the root), and
rewires all missation maps via `missations.ref_seq_changed` so their `from_states` are
correct relative to the new `ref_sequence`. This is necessary because Step 1 used the
original `ref_sequence` to keep missation `from_states` correct during construction.

### Step 6: Validate

1. **`assert_phylo_tree_integrity(tree, true)`**: Full structural validation. Called with
   `force = true` so it runs even in release builds, since this is a one-time initialization
   step.

2. **Tip sequence verification** (release-build check): For each tip, reconstruct its
   absolute sequence by composing mutations from root to tip (starting from
   `tree.ref_sequence`), and verify it matches the absolute sequence implied by
   `orig_ref_sequence` + `tip_descs[i].seq_deltas` at all non-missing sites. This catches
   any bug in the construction pipeline that silently produces a tree inconsistent with the
   input data. Implemented as a standalone function `assert_phylo_tree_matches_tip_descs`
   that takes the tree, the original reference sequence, and the tip descs, and walks the
   tree once. Declared in `core/phylo_tree.h`, defined in `core/phylo_tree.cpp`, so it can
   also be added to the end of `build_usher_like_tree` and `build_random_tree` (both of
   which have the original reference sequence available).

## Entry point

```cpp
auto build_initial_phylo_tree(
    Real_sequence ref_sequence,
    std::vector<Tip_desc> tip_descs,
    absl::BitGenRef bitgen,
    const std::function<void(int,int)>& progress_hook = [](int,int){})
    -> Phylo_tree;
```

This function orchestrates the full pipeline:
1. Build the guide tree (`build_guide_tree`).
2. Root it at the timed midpoint: `auto root = midpoint_root_utree(...)`.
3. Estimate mutation rate and root date: `auto rate = estimate_rate_and_root_date(..., root, ...)`.
4. Convert to `Phylo_tree`: `utree_to_phylo_tree(..., root, ..., rate, ...)`.

Each step is a separate function (declared in `core/utree.h`, defined in `core/utree.cpp`)
for testability. The `progress_hook` is called during guide tree construction (which
dominates runtime); this requires augmenting the `Utree_builder` machinery to accept and
invoke a progress hook.

### Internal functions

```cpp
// Root the Utree at the timed midpoint of a diametral path.  Modifies the tree in place:
// inserts a root node.  Returns the root node index.
// The focus location is undefined after this call.
auto midpoint_root_utree(Utree& tree, const std::vector<Tip_desc>& tip_descs) -> Node_index;

// Estimate overall mutation rate and root date from root-to-tip distances and tip dates.
// `root` is the root node index returned by midpoint_root_utree.
struct Rate_estimate {
  double lambda;   // overall mutations per day (not per-site)
  double t_root;   // estimated root date (days since epoch)
};
auto estimate_rate_and_root_date(
    const Utree& tree,
    Node_index root,
    const std::vector<Tip_desc>& tip_descs)
    -> Rate_estimate;

// Convert a rooted Utree to a Phylo_tree.
// `root` is the root node index returned by midpoint_root_utree.
// Moves the focus to `root` if not already there (hence non-const).
auto utree_to_phylo_tree(
    Utree& utree,
    Node_index root,
    const std::vector<Tip_desc>& tip_descs,
    const Rate_estimate& rate,
    absl::BitGenRef bitgen)
    -> Phylo_tree;
```

## CLI and WASM integration

### CLI option

Add a new `--v0-init` command-line option that selects the tree initialization method:

| Value            | Method                                 |
|------------------|----------------------------------------|
| `random`         | `build_random_tree`                    |
| `old-usher-like` | `build_usher_like_tree`                |
| `mp-plus-timing` | `build_initial_phylo_tree` (this plan) |

Default: `old-usher-like` (to be switched to `mp-plus-timing` in a future round once the
pipeline is proven).

Deprecate the existing `--v0-init-heuristic` and `--v0-init-random` flags, but continue to
accept them as synonyms for `--v0-init old-usher-like` and `--v0-init random` respectively.

**Implementation**: In `core/cmdline.h`, declare `enum class Init_method { random,
old_usher_like, mp_plus_timing }` and update the signature of
`build_rough_initial_tree_from_maple` to take `Init_method` instead of `bool random`. In
`core/cmdline.cpp`, add the new case and the `--v0-init` option parsing, mapping the
deprecated flags to it.

### WASM interface

The WASM tree-building interface works as follows:
- `delphy_parse_fasta_into_initial_tree_async` and `delphy_parse_maple_into_initial_tree_async`
  (in `tools/delphy_wasm.cpp`) call `build_rough_initial_tree_from_maple` with a hardcoded
  `init_random = false`.
- The TypeScript API in `delphy-web` (`src/ts/pythia/delphy_api.ts`) wraps these WASM
  functions but does not expose the init method choice — it is controlled entirely by the
  hardcoded C++ value.

**Changes needed**:

1. **`tools/delphy_wasm.cpp`**: Change `auto init_random = false;` (lines 218 and 292) to
   use the new `Init_method` enum, defaulting to `Init_method::old_usher_like`. This is the
   single-line change that switches the WASM build to the new method when ready.

2. **`delphy-web`** (`~/github/fathominfo/delphy-web`): No TypeScript changes needed. The
   init method is controlled by the C++ WASM code, not by the TypeScript API. When we're
   ready to switch, we change one line in `delphy_wasm.cpp`, rebuild the WASM binary, and
   copy it to `delphy-web/src/ts/delphy/delphy.wasm`.

## Files to create/modify

| File | Action | What |
|------|--------|------|
| `core/utree.h` | Modify | Add `Rate_estimate` struct; declare `midpoint_root_utree`, `estimate_rate_and_root_date`, `utree_to_phylo_tree`, `build_initial_phylo_tree`; change `make_empty` node preallocation from `2N-2` to `2N-1` |
| `core/utree.cpp` | Modify | Implement all new functions; add progress hook support to `Utree_builder` |
| `core/phylo_tree.h` | Modify | Declare `assert_phylo_tree_matches_tip_descs(const Phylo_tree&, const Real_sequence& orig_ref_sequence, const std::vector<Tip_desc>&, bool force = false)` |
| `core/phylo_tree.cpp` | Modify | Implement `assert_phylo_tree_matches_tip_descs`; add call at end of `build_usher_like_tree` and `build_random_tree` |
| `core/cmdline.h` | Modify | Add `Init_method` enum; update `build_rough_initial_tree_from_maple` signature (`bool random` → `Init_method`) |
| `core/cmdline.cpp` | Modify | Update `build_rough_initial_tree_from_maple` implementation; add `--v0-init` option; deprecate old flags |
| `tools/delphy_wasm.cpp` | Modify | Change `init_random` to `Init_method` enum (lines 218, 292) |
| `tests/utree_tests.cpp` | Modify | Add tests for midpoint rooting, rate estimation, conversion, and end-to-end pipeline |

## Verification

1. **Midpoint rooting tests**:
   - 2-tip tree: root splits the single edge correctly.
   - 3-tip tree with known edge sizes: verify the root is on the correct edge and `D` is
     correct.
   - 3-tip tree with unequal tip dates: verify the timed midpoint shifts the root toward
     the older tip.
2. **Rate estimation tests**:
   - Synthetic data with known rate: verify `lambda` and `t_root` match.
   - All-same-date tips: verify fallback is used and produces finite values.

3. **Conversion tests**:
   - Small trees (3-5 tips): verify the resulting `Phylo_tree` passes
     `assert_phylo_tree_integrity` and `assert_phylo_tree_matches_tip_descs`.
   - Verify tip names, date bounds, and missations are preserved.
   - Verify parent times are strictly before child times.

4. **Manual end-to-end test**: Run `delphy_ui` with `--v0-init mp-plus-timing` against each
   standard demo dataset. Spot-check that the results are sensible and that the
   `assert_phylo_tree_matches_tip_descs` CHECKs don't trigger.

5. **Build**: `cmake --build build/debug --target tests -j$(nproc) && build/debug/tests/tests --gtest_filter='Utree*'`

## Observations from manual testing

For some datasets (e.g., mpox [Parker et al, 2025]), the resulting trees are great. For
others, two issues arise:

1. **Bad midpoint root**: For datasets like H5N1 (Andersen), the midpoint root produces a
   tMRCA many years before the actual tMRCA. This is inherent to midpoint rooting as a
   heuristic — it doesn't account for rate variation or long stems. Should be largely fixed
   by Round 4's root-to-tip regression over candidate roots.

2. **Ladder-like artifacts for closely related sequences**: For datasets like SARS-CoV-2
   (Lemieux et al, 2021) and Ebola (Gire et al, 2014), closely related tips attached
   sequentially produce degenerate ladder topologies. This happens because when multiple
   attachment points have equal cost, the builder always picks the first one found. Two
   fixes for a small follow-up commit: (a) randomize the attachment point among
   equal-best-cost alternatives in `find_best_attachment_arc`; (b) reduce
   `min_branch_length` in the Utree → Phylo_tree conversion cleanup pass.
