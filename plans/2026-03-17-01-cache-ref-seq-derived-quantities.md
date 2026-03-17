# Plan: Cache ref-sequence-derived quantities in Run, share with Subruns

**Status: SHELVED** — not worth the surgical risk during paper review.

## Problem

`calc_cum_Q_l_for_sequence` and `calc_state_frequencies_per_partition_of`
are O(L) scans over `ref_sequence` that depend only on `ref_sequence`
and `evo` — quantities that are identical across all subruns and
already known at the `Run` level.

All subruns share the same `ref_sequence`: in `repartition()` (line
136 of `core/run.cpp`), every subtree gets `subtree.ref_sequence =
ref_seq` where `ref_seq = tree_.ref_sequence`.  No subrun ever
modifies `ref_sequence` — subruns only modify mutations, missations,
topology, and node times.  (The `DCHECK` at line 246 of `reassemble()`
only checks the root-containing subrun, so it's not strong
confirmation by itself; the real guarantee is that no code path in
`Subrun` writes to `ref_sequence`.)

Each subrun independently recomputes these in
`Subrun::recalc_derived_quantities()` after every global move cycle
(~1,000 times over 20M steps for 200 tips).  With `num_parts = 1`
(default), there is one redundant recomputation per cycle; with N
parts, N redundant recomputations.

`Run::recalc_derived_quantities()` (called at line 629 in
`do_mcmc_steps`) already calls `derive_evo()` and computes
`evo_`-dependent quantities *before* `push_global_params_to_subruns()`
(line 641).  So by the time subruns need these values, the `Run`
already has the current `evo_`.

## Profiling notes

Initial profiling with default `--v0-log-every` (= steps/10000)
exaggerated the cost of these functions because each
`do_mcmc_steps()` call triggers a `repartition()`, and with
`step_granularity = 2000` there were ~10,000 repartitions in a 20M
step run.  In realistic runs, `--v0-log-every` is much larger and
repartitions are far less frequent.

With realistic repartition frequency (`--v0-log-every 200000`,
~100 repartitions over 20M steps, `06_missing_data` sim_000):

- **97B cpu_core cycles, 17s wall time** (vs 255B / 39s with default
  log-every)
- `calc_cum_Q_l_for_sequence`: **10.7%** (was 20.3%)
- `calc_state_frequencies_per_partition_of`: **7.1%** (was 13.9%)
- Combined: **17.8%** (was 34%)
- Expected speedup from caching: **~1.2x** (was ~1.3-1.5x)

For comparison, the top MCMC move hotspots with realistic frequency:
- `inner_node_displace_move`: 9.5%
- `__nextafter`: 4.9%
- `__ieee754_log_fma`: 4.1%
- `calc_Ttwiddle_beta_a`: 2.7%
- `branch_reform_move`: 2.3%

## Goal

Compute `ref_cum_Q_l` and
`state_frequencies_of_ref_sequence_per_partition` once in `Run` after
`evo_` changes, and share them read-only with subruns via const
pointers, so subruns skip the redundant O(L) recomputations entirely.

## Changes

### 1. Cache in `Run::recalc_derived_quantities()`

Add two cached members to `Run`:

```cpp
mutable std::vector<double> ref_cum_Q_l_{};
mutable Partition_vector<Seq_vector<int>>
    state_frequencies_of_ref_sequence_per_partition_{};
```

In `Run::recalc_derived_quantities()`, compute them right after
`derive_evo()`:

```cpp
ref_cum_Q_l_ = calc_cum_Q_l_for_sequence(tree_.ref_sequence, evo_);
state_frequencies_of_ref_sequence_per_partition_ =
    calc_state_frequencies_per_partition_of(tree_.ref_sequence, evo_);
```

Note: `Run` already has `state_frequencies_of_ref_sequence_`
(unpartitioned, `Seq_vector<int>`), computed via
`calc_cur_state_frequencies_of_ref_sequence()` which calls
`calc_state_frequencies_per_partition_of(...)[0]`.  Update it to
derive from the cached partitioned version:

```cpp
state_frequencies_of_ref_sequence_ =
    state_frequencies_of_ref_sequence_per_partition_[0];
```

This works because `state_frequencies_of_ref_sequence_` is only used
in the single-partition case (the mpox hack has 2 partitions but does
not use this member for anything partition-sensitive — verify).

### 2. Share cached values with subruns via const pointers

Add const-pointer members to `Subrun`:

```cpp
const std::vector<double>* shared_ref_cum_Q_l_ = nullptr;
const Partition_vector<Seq_vector<int>>*
    shared_state_frequencies_of_ref_sequence_per_partition_ = nullptr;
```

Add a single setter:

```cpp
auto set_shared_ref_seq_derived_quantities(
    const std::vector<double>& ref_cum_Q_l,
    const Partition_vector<Seq_vector<int>>&
        state_frequencies_of_ref_sequence_per_partition) -> void {
  shared_ref_cum_Q_l_ = &ref_cum_Q_l;
  shared_state_frequencies_of_ref_sequence_per_partition_ =
      &state_frequencies_of_ref_sequence_per_partition;
}
```

In `push_global_params_to_subruns()`, call this after `set_evo()`:

```cpp
subrun.set_shared_ref_seq_derived_quantities(
    ref_cum_Q_l_,
    state_frequencies_of_ref_sequence_per_partition_);
```

### 3. Use shared values in `Subrun::recalc_derived_quantities()`

In `Subrun::recalc_derived_quantities()`, use the shared values
instead of recomputing:

```cpp
// Before:
state_frequencies_of_ref_sequence_per_partition_ =
    calc_cur_state_frequencies_of_ref_sequence_per_partition();
ref_cum_Q_l_ = calc_cur_ref_cum_Q_l();

// After:
if (shared_state_frequencies_of_ref_sequence_per_partition_) {
  state_frequencies_of_ref_sequence_per_partition_ =
      *shared_state_frequencies_of_ref_sequence_per_partition_;
} else {
  state_frequencies_of_ref_sequence_per_partition_ =
      calc_cur_state_frequencies_of_ref_sequence_per_partition();
}
if (shared_ref_cum_Q_l_) {
  ref_cum_Q_l_ = *shared_ref_cum_Q_l_;
} else {
  ref_cum_Q_l_ = calc_cur_ref_cum_Q_l();
}
```

The fallback (when pointers are null) keeps `Subrun` usable
standalone, e.g. in tests.  The `calc_cur_*` methods remain for
`check_derived_quantities()` in debug builds.

### 4. Lifetime safety

The shared pointers point into `Run`'s mutable cached members.  These
are recomputed in `Run::recalc_derived_quantities()`, which runs at
line 629 of `do_mcmc_steps` — *before*
`push_global_params_to_subruns()` (line 641) and `run_local_moves()`
(line 643).  The cached values are stable for the entire local-moves
phase.  The pointers are refreshed every time
`push_global_params_to_subruns()` is called, so they always point to
the current values.

After `repartition()`, `push_global_params_to_subruns()` runs before
any subrun iteration (line 641 in the main loop), so newly
constructed subruns will have their pointers set before first use.

## Files changed

- `core/run.h` — add `ref_cum_Q_l_` and
  `state_frequencies_of_ref_sequence_per_partition_` cached members
- `core/run.cpp` — compute in `recalc_derived_quantities()`, pass in
  `push_global_params_to_subruns()`, derive
  `state_frequencies_of_ref_sequence_` from cached partitioned version
- `core/subrun.h` — add shared const-pointer members and setter
- `core/subrun.cpp` — use shared values in
  `recalc_derived_quantities()` when available

## Validation

1. **Baseline profile** (with realistic log-every): 97B cpu_core
   cycles, 17s wall time on `06_missing_data` sim_000, 20M steps,
   `--v0-log-every 200000`.
2. Run the existing Delphy test suite (`build/release/tests/tests`).
3. **After profile**: re-profile on the same dataset and workload,
   compare cpu_core cycles and wall time to the baseline.  Expect
   `calc_cum_Q_l_for_sequence` and
   `calc_state_frequencies_per_partition_of` to drop from ~17.8%
   combined to near-zero in subrun context.

## Expected impact

With realistic repartition frequency, these two functions account for
~17.8% of cpu_core time.  With sharing, each is computed once per
global move cycle (in `Run`) instead of once per subrun per cycle.
Expected speedup: ~1.2x.
