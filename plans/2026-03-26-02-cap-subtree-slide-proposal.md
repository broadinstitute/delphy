# Cap Subtree Slide Proposal to Prevent OOM in Coalescent Prior

## Problem

After fixing the unbounded proposal scale in `inner_node_displace_move()`
(see `2026-03-26-01-cap-root-displace-proposal.md`), a separate code path
in `subtree_slide_move()` exhibits the same bug: when `mu` is tiny, its
Gaussian proposal scale blows up, proposing absurdly distant parent times
that cause the coalescent prior to allocate enormous grids.

This bug is **independent** of both the `sample_mutational_history` infinite
loop fix (`2026-03-25-01-fix-sample-mutational-history-infinite-loop.md`)
and the `inner_node_displace_move` fix, though all three arise from the
same root cause (tiny `mu`).

## Investigation

### How the bug was found

After applying the `inner_node_displace_move` fix, we re-ran `sim_048`
of WCSS `04_free_exp_pop` in `delphy-2026-paper-data`.  Memory spiked
to ~2 GiB, the process froze, and we caught it with SIGSTOP.

### What gdb revealed

Thread 3 was stuck in:

```
#0  std::vector<double>::push_back
#1  Very_scalable_coalescent_prior_part::ensure_space (t=-176968602.06)
#4  Very_scalable_coalescent_prior_part::calc_delta_partial_log_prior_after_displace_coalescence
        (old_t=2165.22, new_t=-176968602.06)
#5  Subrun::spr_move_core (X=47, new_nexus.t=-176968602.06)
#6  Subrun::subtree_slide_move                              at subrun.cpp:397
```

Key state from `subtree_slide_move` (frame 6):

- **`old_P_t = 2165.22`**: the parent's current time
- **`new_P_t = -176,968,602.06`**: the proposed new parent time
  (177 million years in the past!)
- **`mu_proposal = 1.04e-13`**: absurdly small

### Root cause

`subtree_slide_move` picks a random non-root node `X`, then proposes
sliding its parent `P` up or down the tree by a Gaussian displacement.
`S` is the sibling of `X` (the other child of `P`).  If the
displacement moves `P` past a neighboring node, the topology changes
via an SPR (subtree prune-and-regraft) operation handled by
`spr_move_core`.

The Gaussian proposal scale is set on line 365 of `subrun.cpp`:

```cpp
auto delta_scale = (1 / lambda_i_.at(X)) / 2;  // 95% of the time, won't hop past more than 1 mutation
auto delta_t = absl::Gaussian(bitgen_, 0.0, delta_scale);
auto new_P_t = tree_.at(P).t + delta_t;
```

This is the exact same pattern as the `inner_node_displace_move` bug:
with `mu = 1e-13` and `L = 30,000`, `lambda_i ≈ 3e-9`, so
`delta_scale ≈ 1.7e8` (170 million years).  A 1-sigma Gaussian draw
easily proposes a parent time hundreds of millions of years away.

When `P` is at or near the root and the move goes "up" (delta_t < 0),
the while loop at lines 378-384 walks up the tree looking for a new
parent, but when it reaches the root (`GG == k_no_node`), it simply
breaks — there is no bound on `new_P_t`.  The proposed time is then
passed to `spr_move_core`, which calls the coalescent prior's
`calc_delta_partial_log_prior_after_displace_coalescence`, triggering
`ensure_space()` to allocate a grid covering 177 million years.

The proposal will almost certainly be rejected by the coalescent prior,
so this is pure waste.

### Why `spr1_move` is not affected

The other SPR move variant (`spr1_move`) uses an `Spr_study` to pick
the new parent time, which already caps the above-root region to
`s_max = s_min + 20 * tree_span` (spr_study.cpp:447).  So `spr1_move`
is already protected.

## Proposed fix

Two changes:

### 1. Cap `delta_scale` in `subtree_slide_move`

Cap `delta_scale` to `tree_span`, exactly as was done for
`inner_node_displace_move`.  The tree span is computed as
`t_max_tip_ - t_early`, where `t_early` is the earliest node in the
tree excluding P:

- If `P` is **not** the root: `t_early = tree_.at(tree_.root).t`
  (the root time).
- If `P` **is** the root: `t_early = std::min(tree_.at(X).t, tree_.at(S).t)`
  (the earlier of the two children of the root).

Using `t_early` instead of `old_P_t` is important: `delta_scale` must
not depend on the current position of the node being moved, or the
Gaussian proposal would no longer be symmetric and detailed balance
would be broken.  With `t_early` defined as above, its value is the
same before and after the move (the children's times don't change in
this move), so the proposal remains symmetric.

```cpp
auto t_early = (P == tree_.root)
    ? std::min(tree_.at(X).t, tree_.at(S).t)
    : tree_.at(tree_.root).t;
auto tree_span = t_max_tip_ - t_early;
CHECK_GE(tree_span, 0.0);
auto delta_scale = std::min((1 / lambda_i_.at(X)) / 2,  // 95% of the time, won't hop past more than 1 mutation
                            tree_span);
```

This ensures that the Gaussian proposal scale is at most the tree span,
preventing absurd jumps while keeping moves ergodic: it remains
*possible* but extremely unlikely for the parent to be pushed many
factors of `tree_span` into the past.  For normal `mu` values,
`delta_scale` is already much smaller than `tree_span`, so the cap
has no effect.

Note: `t_max_tip_` is already available in the Subrun (line 90 of
`subrun.h`).

### 2. Add a diagnostic in `ensure_space`

Add a warning log when `ensure_space()` allocates more than 10,000
new cells.  This makes it easier to detect similar issues in the
future without having to use gdb.

## Implementation

### Step 1: Cap `delta_scale` (subrun.cpp line 365)

Replace lines 365-368:
```cpp
auto delta_scale = (1 / lambda_i_.at(X)) / 2;  // 95% of the time, won't hop past more than 1 mutation
auto delta_t = absl::Gaussian(bitgen_, 0.0, delta_scale);
auto old_P_t = tree_.at(P).t;
auto new_P_t = tree_.at(P).t + delta_t;
```

with:
```cpp
// Cap the proposal scale to the tree span so that we don't propose
// absurdly distant parent times when mu is tiny (making lambda_i
// negligible and delta_scale enormous).  Such proposals would be
// rejected by the coalescent prior anyway, but evaluating the prior
// can consume GiB of memory.
// See plans/2026-03-26-02-cap-subtree-slide-proposal.md for details.
auto t_early = (P == tree_.root)
    ? std::min(tree_.at(X).t, tree_.at(S).t)
    : tree_.at(tree_.root).t;
auto tree_span = t_max_tip_ - t_early;
CHECK_GE(tree_span, 0.0);
auto delta_scale = std::min((1 / lambda_i_.at(X)) / 2,  // 95% of the time, won't hop past more than 1 mutation
                            tree_span);
auto delta_t = absl::Gaussian(bitgen_, 0.0, delta_scale);
auto old_P_t = tree_.at(P).t;
auto new_P_t = old_P_t + delta_t;
```

### Step 2: Add diagnostic in `ensure_space` (very_scalable_coalescent.cpp)

Add `#include "absl/strings/str_format.h"` to
`very_scalable_coalescent.cpp` (after the existing `#include <iostream>`).

Then, in `Very_scalable_coalescent_prior_part::ensure_space()`, before
the allocation loop at line 267, add a warning when the number of new
cells to allocate exceeds 10,000:

```cpp
auto num_new_cells = max_cell - static_cast<int>(std::ssize(popsize_bar_)) + 1;
if (num_new_cells > 10'000) {
  std::cerr << absl::StreamFormat(
      "WARNING: ensure_space: allocating %d new cells for t=%.6g"
      " (t_ref=%.6g, t_step=%.6g)\n",
      num_new_cells, t, t_ref_, t_step_);
}
```

### Step 3: Build and test

Build and run the test suite.  Then re-run with `sim_048` and confirm
that memory stays stable (no spikes).

## Scope

- Modify `subrun.cpp`: replace `delta_scale` computation on line 365
  with a capped version.
- Modify `very_scalable_coalescent.cpp`: add a diagnostic warning in
  `ensure_space()` for large allocations.
- No new tests needed — the existing test suite covers
  `subtree_slide_move` indirectly, and the fix is a tighter
  proposal scale that cannot break correctness (proposals are just
  smaller, not biased; the Gaussian proposal is symmetric so the
  Hastings ratio remains 1).
