# Plan: Gate release-build consistency checks with `--v0-paranoid`

## Problem

Several `assert_phylo_tree_integrity(tree, true)` calls force full
tree integrity checks (including `assert_missation_consistency`) even
in release builds.  Profiling a WCSS run (200 tips, L=30,000, ~5%
missing data, site-rate heterogeneity) showed these checks consuming
~23% of cpu_core runtime:

- `assert_missation_consistency`: 16.4%
- `Check_EQImpl` (called from within it): 6.3%

These checks walk every node in the tree, iterate over all missing
sites, and verify sequence consistency — all redundant work for a
correct implementation.  They have not failed for many months of
intense Delphy use across diverse datasets.

Separately, `assert_space` in `very_scalable_coalescent.cpp` has its
`estd::is_debug_enabled` guard commented out, so its `CHECK_GE` /
`CHECK_LT` calls run unconditionally in release builds (1.55% of
cpu_core time).

## Goal

1. Make the forced `assert_phylo_tree_integrity` checks opt-in via a
   new `--v0-paranoid` CLI flag.  Without the flag, release builds
   skip these checks (debug builds always run them regardless).
2. Re-enable the `estd::is_debug_enabled` guard in `assert_space` so
   it only runs in debug builds.  This is not gated by `--v0-paranoid`
   — it's a straightforward bug fix (the guard was accidentally
   commented out).

## Changes

### 1. Re-enable `estd::is_debug_enabled` in `assert_space`

In `very_scalable_coalescent.cpp`, uncomment the
`if (estd::is_debug_enabled)` guard around `assert_space`:

```cpp
// Before:
//if (estd::is_debug_enabled) {
    auto cell = cell_for(t, t_ref, t_step);
    ...
    //}

// After:
if (estd::is_debug_enabled) {
    auto cell = cell_for(t, t_ref, t_step);
    ...
}
```

This is independent of the paranoid flag — `assert_space` should
simply never run in release builds.

### 2. Add `paranoid_` flag to `Run`

In `run.h`, add a `bool paranoid_` member (default `false`) and a
setter/getter:

```cpp
auto set_paranoid(bool v) -> void { paranoid_ = v; }
auto paranoid() const -> bool { return paranoid_; }
```

### 3. Replace `force=true` with `paranoid_` in `Run`

In `run.cpp`, the two hot `assert_phylo_tree_integrity` calls are:

- In `repartition()` (line ~127): called during tree repartitioning
- In `do_mcmc_steps()` (line ~652): called after every batch of local
  moves (every `step_granularity` steps)

Change both from `assert_phylo_tree_integrity(tree_, true)` to
`assert_phylo_tree_integrity(tree_, paranoid_)`.

### 4. Check subtree integrity before `Subrun` construction

In `run.cpp`, in `repartition()`, add
`assert_phylo_tree_integrity(subtree, paranoid_)` just before each
`subruns_.emplace_back(...)` call.  This verifies the tree-splitting
logic is correct without threading `paranoid` into `Subrun` itself.
`Subrun`'s constructor no longer calls `assert_phylo_tree_integrity`.

### 5. Keep `force=true` in tree-building functions

The 3 calls in `phylo_tree.cpp` (lines 256, 721, 987) are in
one-shot tree-construction functions
(`push_all_mutations_and_missations_to_tips`, `build_random_tree`,
`build_usher_like_tree`).  These run once during startup, not in the
MCMC loop, and are not performance-critical.  Leave them as
`force=true`.

### 6. Parse `--v0-paranoid` in `cmdline.cpp`

Add a boolean CLI option `--v0-paranoid` (default false) after
`--v0-seed`.  After constructing the `Run` object, call
`run->set_paranoid(...)` with the parsed value.

## Files changed

- `core/very_scalable_coalescent.cpp` — re-enable `is_debug_enabled`
  guard in `assert_space`
- `core/run.h` / `core/run.cpp` — add `paranoid_` member, use in
  integrity checks, add subtree check before subrun construction
- `core/subrun.cpp` — remove `assert_phylo_tree_integrity` from
  constructor (moved to `run.cpp`)
- `core/cmdline.cpp` — parse `--v0-paranoid`, call
  `run->set_paranoid(...)`

## Validation

1. Build and run tests (`ctest`).
2. Run a short Delphy job **without** `--v0-paranoid` — verify
   integrity checks do not appear in a perf profile.
3. Run a short Delphy job **with** `--v0-paranoid` — verify the
   checks do run (will be visibly slower, ~23% overhead).
