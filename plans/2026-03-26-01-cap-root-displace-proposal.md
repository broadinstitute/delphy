# Cap Root Displace Proposal to Prevent OOM in Coalescent Prior

## Problem

When Delphy runs on a dataset with very few or zero mutations (e.g.,
WCSS `04_free_exp_pop/sims/sim_048` of `delphy-2026-paper-data`), the
Gibbs sampler for `mu` drives it extremely small (observed:
`mu = 4.96e-14`).  This causes `inner_node_displace_move()` in
`core/subrun.cpp` to propose absurdly distant root times, which in
turn causes the coalescent prior's `ensure_space()` to allocate
enormous grid structures, consuming multiple GiB of memory and
freezing the process for several seconds.

This bug is **independent** of the `sample_mutational_history` infinite
loop fixed in `2026-03-25-01-fix-sample-mutational-history-infinite-loop.md`,
though it arises from the same root cause (tiny `mu`).

## Investigation

### How the bug was found

After applying the `sample_mutational_history` fix, we ran `sim_048`
again.  Memory spiked to ~5.5% (~3.4 GiB on a 62 GiB machine),
the process froze, and we caught it with SIGSTOP.

### What gdb revealed

One worker thread was stuck in:

```
__GI___exp (x=-3145561.9)
Exp_pop_model::pop_integral (a=-83495755.1, b=-83495754.0)
Very_scalable_coalescent_prior_part::ensure_space (t=-207318156.8)
  <- calc_delta_partial_log_prior_after_displace_coalescence
       (old_t=1828.99, new_t=-207318156.8)
  <- Subrun::inner_node_displace_move
```

Key state from `inner_node_displace_move` (frame 7):

- **`old_node_t = 1828.99`**: the root's current time
- **`new_node_t = -207,318,156.80`**: the proposed new root time
  (207 million years in the past!)
- **`t_min = -inf`**: no lower bound (root has no parent)
- **`t_max = inf`**: gdb showed `inf`, but this is likely a release
  build artifact (optimized-out or stale register).  The root
  definitely has children, so `t_max` should be finite; the code at
  lines 163-169 iterates over children and takes the min of their
  times and mutation times.
- **`mu = 4.96e-14`**: absurdly small

### Root cause

For the root node, `inner_node_displace_move` uses a Gaussian
proposal (lines 189-190 of `subrun.cpp`):

```cpp
auto delta_scale = (1 / lambda_i_.at(node)) / 2;
new_node_t = old_node_t + absl::Gaussian(bitgen_, 0.0, delta_scale);
```

With `mu = 5e-14` and `L = 30,000`:
- `lambda_i ≈ mu * L ≈ 1.5e-9`
- `delta_scale = 1 / (2 * 1.5e-9) ≈ 3.3e8` (330 million years!)

A 1-sigma Gaussian draw can easily propose a root time hundreds of
millions of years away.  The coalescent prior's `ensure_space()` then
tries to build a grid covering that range:

```cpp
auto max_cell = cell_for(t, t_ref_, t_step_);  // (t_ref - t) / t_step => enormous
for (auto i = std::ssize(popsize_bar_); i <= max_cell; ++i) {
    // allocates one entry per cell
}
```

This allocates millions of cells, consuming GiB of memory.

The proposal will almost certainly be rejected by the coalescent prior
(the prior strongly penalizes such distant roots), so this is pure
waste: enormous allocation followed by immediate deallocation.

### Why the proposal scale is unbounded

The comment on line 189 says `delta_scale` is chosen so that "95% of
the time, won't hop past more than 1 mutation."  This is sensible
when `mu` is normal (say 1e-3/site/year): `lambda ≈ 30`, so
`delta_scale ≈ 0.017 years` (~6 days).  But when `mu` is absurdly
small, `lambda` is tiny and the scale blows up — there are effectively
no mutations to hop past, so the scale becomes the entire age of the
universe.

## Proposed fix

Cap `delta_scale` to `tree_span`, where `tree_span` is the time from
the earliest of the root's children to the latest tip.  The cap is
applied at the declaration site so that the final value of
`delta_scale` is visible in one place:

```cpp
auto tree_span = t_max_tip_ - t_max;
CHECK_GE(tree_span, 0.0);
auto delta_scale = std::min((1 / lambda_i_.at(node)) / 2,  // 95% of the time, won't hop past more than 1 mutation
                            tree_span);
```

Here `t_max` is already computed (lines 163-169) as the minimum of
the root's children's times and their mutation times — i.e., the
earliest point the root could be moved forward to.  So `tree_span`
is the time span of the tree below the root.

This ensures that the Gaussian proposal scale is at most the tree
span, preventing absurd jumps while keeping moves ergodic: it remains
*possible* but extremely unlikely for the root to be pushed many
factors of `tree_span` into the past.  For normal `mu` values,
`delta_scale` is already much smaller than `tree_span`, so the cap
has no effect.

Note: `t_max_tip_` is already available in the Subrun (line 90 of
`subrun.h`).

## Implementation

### Step 1: Cap `delta_scale` (subrun.cpp line 189)

Replace the existing `delta_scale` computation on line 189:
```cpp
auto delta_scale = (1 / lambda_i_.at(node)) / 2;  // 95% of the time, won't hop past more than 1 mutation
```

with:
```cpp
// Cap the proposal scale to the tree span so that we don't propose
// absurdly distant root times when mu is tiny (making lambda_i
// negligible and delta_scale enormous).  Such proposals would be
// rejected by the coalescent prior anyway, but evaluating the prior
// requires allocating a grid covering the proposed time range, which
// can consume GiB of memory.
// See plans/2026-03-26-01-cap-root-displace-proposal.md for details.
auto tree_span = t_max_tip_ - t_max;
CHECK_GE(tree_span, 0.0);
auto delta_scale = std::min((1 / lambda_i_.at(node)) / 2,  // 95% of the time, won't hop past more than 1 mutation
                            tree_span);
```

### Step 2: Build and test

Build and run the test suite.  Then re-run with `sim_048` and confirm
that memory stays stable (no spikes to several GiB).

## Scope

- Modify `subrun.cpp` only: replace `delta_scale` computation on
  line 189 with a capped version.
- No new tests needed — the existing test suite covers
  `inner_node_displace_move` indirectly, and the fix is a tighter
  proposal scale that cannot break correctness (proposals are just
  smaller, not biased; the Gaussian proposal is symmetric so the
  Hastings ratio remains 1).
- Do **not** change `ensure_space` — the root cause is the
  unbounded proposal scale, not the grid allocation logic.
