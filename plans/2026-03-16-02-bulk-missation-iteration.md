# Plan: Replace hot `slow_elements` calls with interval-level iteration

## Problem

`Missation_map::slow_elements()` is a coroutine (`cppcoro::generator`)
that yields individual `{site, from_state}` pairs by iterating
site-by-site through each interval.  With ~1,500 missing sites per tip
across ~400 nodes, this produces ~600K coroutine suspend/resume cycles
per full tree traversal.

Profiling a WCSS run (200 tips, L=30,000, ~5% missing data,
site-rate heterogeneity enabled) showed `slow_elements` consuming
**25.6%** of cpu_core time — the single largest hotspot.  The coroutine
overhead (suspend/resume/frame allocation) dominates the actual work.

## Background: `Missation_map` structure

A `Missation_map` stores missing sites as two public fields:

- `intervals`: an `Interval_set` of `[start, end)` pairs (contiguous
  gaps).
- `from_states`: a `flat_map<Site_index, Real_seq_letter>` containing
  only sites where `from_state != ref_sequence[site]`.  This map is
  almost always empty or tiny (0--2 entries per node), since only a
  small fraction of sites are missing at any one node, and mutations
  between a node and the root sequence are rare at those few sites.

`slow_elements()` iterates site-by-site over `intervals`, calling
`get_from_state(l, ref_sequence)` for each site — which does a
`from_states.find()` lookup on every iteration.

## Affected call sites

There are 7 `slow_elements` call sites in `phylo_tree_calc.cpp`.
Of these, 3 are hot in production, 3 are only used in tests (but
should be fixed for consistency since the test functions are the
canonical implementations that the production variants are based on),
and 1 is production code but not hot.

### Hot production call sites

These account for essentially all of the 25.6% overhead.

**`calc_Ttwiddle_beta_a`** (`core/phylo_tree_calc.cpp`, lines 288, 327) — accumulate into a
partitioned vector by `from_state`.  Called from
`Run::recalc_derived_quantities()` (every ~1M steps) and
`Run::gibbs_sample_all_nus()` (every global move cycle).  Profiling
shows ~18.5% of `slow_elements` time attributed to this function.

1. **entry** (line 288):
   ```cpp
   for (const auto& [mi_site, mi_from] : tree.at(node).missations.slow_elements(tree.ref_sequence)) {
       auto beta = evo.partition_for_site[mi_site];
       ntwiddle_beta_a[beta][mi_from] -= evo.nu_l[mi_site];
   }
   ```

2. **exit** (line 327):
   ```cpp
   for (const auto& [mi_site, mi_from] : tree.at(node).missations.slow_elements(tree.ref_sequence)) {
       auto beta = evo.partition_for_site[mi_site];
       ntwiddle_beta_a[beta][mi_from] += evo.nu_l[mi_site];
   }
   ```

**`calc_Ttwiddle_l`** (`core/phylo_tree_calc.cpp`, line 202) — per-site indexed write.  Called
from `Run::gibbs_sample_all_nus()` and `Run::alpha_moves()` (both
global moves).  Profiling shows ~4.9% of `slow_elements` time here.

3. **line 202**:
   ```cpp
   for (const auto& [mi_site, mi_from] : tree.at(node).missations.slow_elements(tree.ref_sequence)) {
       auto T_below_miss = ...;
       Ttwiddle_l[mi_site] -= evo.q_l_a(mi_site, mi_from) * T_below_miss;
   }
   ```

### Test-only call sites (fix for consistency)

These functions are never called in production code, only from tests
in `phylo_tree_calc_tests.cpp`.  They are the canonical
implementations that `calc_Ttwiddle_beta_a` is a partitioned variant
of, so they should be fixed for consistency and to keep the tests
exercising the same code patterns.

**`calc_T_l_a`** (`core/phylo_tree_calc.cpp`, line 160) — per-site indexed write:

4. **line 160**:
   ```cpp
   for (const auto& [mi_site, mi_from] : tree.at(node).missations.slow_elements(tree.ref_sequence)) {
       auto T_below_miss = ...;
       T_l_a[mi_site][mi_from] -= T_below_miss;
   }
   ```

**`calc_Ttwiddle_a`** (`core/phylo_tree_calc.cpp`, lines 225, 257) — accumulate into a
4-element vector by `from_state`:

5. **entry** (line 225):
   ```cpp
   for (const auto& [mi_site, mi_from] : tree.at(node).missations.slow_elements(tree.ref_sequence)) {
       dTtwiddle_a_dt[mi_from] -= nu_l[mi_site];
   }
   ```

6. **exit** (line 257):
   ```cpp
   for (const auto& [mi_site, mi_from] : tree.at(node).missations.slow_elements(tree.ref_sequence)) {
       dTtwiddle_a_dt[mi_from] += nu_l[mi_site];
   }
   ```

### Not hot but fix for consistency

**`calc_log_root_prior`** (`core/phylo_tree_calc.cpp`, line 446) — accumulate by partition and
state, but only over the root node's missations (typically empty or
near-empty).  Called from `Subrun::calc_cur_log_G()` and
`Run::calc_cur_log_G()`.

7. **line 446**:
   ```cpp
   for (const auto& [mi_site, mi_from] : tree.at_root().missations.slow_elements(tree.ref_sequence)) {
       auto partition = evo.partition_for_site[mi_site];
       --state_frequencies_at_root_per_partition[partition][mi_from];
   }
   ```

## Proposed approach

For each call site, replace the `slow_elements` loop with:

1. A direct loop over `missations.intervals` (the raw `[start, end)`
   pairs), processing sites in bulk assuming `from_state ==
   ref_seq[l]`.
2. A small fixup loop over `missations.from_states` to correct the
   (rare) exceptions where `from_state != ref_seq[l]`.

### Pattern for per-site indexed writes (sites 3, 4)

```cpp
// calc_T_l_a: subtract T_below_miss from T_l_a[site][from]
auto T_below_miss = T_below_node[node] +
    (node == tree.root ? 0.0 : tree.at(node).t - tree.at_parent_of(node).t);
for (const auto& [start, end] : tree.at(node).missations.intervals) {
    for (auto l = start; l != end; ++l) {
        T_l_a[l][tree.ref_sequence[l]] -= T_below_miss;
    }
}
// Fix up exceptions where from_state != ref_seq[l]
for (const auto& [l, from] : tree.at(node).missations.from_states) {
    T_l_a[l][tree.ref_sequence[l]] += T_below_miss;  // undo ref_seq assumption
    T_l_a[l][from] -= T_below_miss;                  // apply correct from_state
}
```

The `calc_Ttwiddle_l` variant (site 3) is similar but the per-site
expression involves `evo.q_l_a`, which returns different values for
different states.  The fixup must subtract and add different scalar
values from the same output location:

```cpp
// calc_Ttwiddle_l: subtract q_l_a(site, from) * T_below_miss from Ttwiddle_l[site]
auto T_below_miss = T_below_node[node] +
    (node == tree.root ? 0.0 : tree.at(node).t - tree.at_parent_of(node).t);
for (const auto& [start, end] : tree.at(node).missations.intervals) {
    for (auto l = start; l != end; ++l) {
        Ttwiddle_l[l] -= evo.q_l_a(l, tree.ref_sequence[l]) * T_below_miss;
    }
}
for (const auto& [l, from] : tree.at(node).missations.from_states) {
    Ttwiddle_l[l] += evo.q_l_a(l, tree.ref_sequence[l]) * T_below_miss;  // undo ref_seq assumption
    Ttwiddle_l[l] -= evo.q_l_a(l, from) * T_below_miss;                  // apply correct from_state
}
```

### Pattern for accumulate by `from_state` (sites 1, 2, 5, 6, 7)

```cpp
// calc_Ttwiddle_beta_a entry: ntwiddle_beta_a[beta][from] -= nu_l[site]
for (const auto& [start, end] : tree.at(node).missations.intervals) {
    for (auto l = start; l != end; ++l) {
        auto beta = evo.partition_for_site[l];
        ntwiddle_beta_a[beta][tree.ref_sequence[l]] -= evo.nu_l[l];
    }
}
for (const auto& [l, from] : tree.at(node).missations.from_states) {
    auto beta = evo.partition_for_site[l];
    ntwiddle_beta_a[beta][tree.ref_sequence[l]] += evo.nu_l[l];  // undo ref_seq assumption
    ntwiddle_beta_a[beta][from] -= evo.nu_l[l];                  // apply correct from_state
}
```

The exit loops use `+=` instead of `-=` (and the fixup signs flip
accordingly).

The unpartitioned variants (sites 5, 6) drop the
`evo.partition_for_site` lookup.  The `calc_log_root_prior` variant
(site 7) uses integer `--`/`++` instead of floating-point `-=`/`+=`.

## Scope

- Change the 7 `slow_elements` call sites in `core/phylo_tree_calc.cpp`
  (lines 160, 202, 225, 257, 288, 327, 446).
- Do **not** change `slow_elements` itself or other callers
  (assertions, tree manipulation, tests) -- they are not hot.
- Do **not** add new methods to `Missation_map` -- `intervals` and
  `from_states` are already public.

## Validation

1. **Baseline profile** (before any changes): run a 20M-step Delphy
   job on the WCSS `06_missing_data` sim_000 dataset (200 tips,
   L=30,000, ~5% missing data, site-rate heterogeneity) and record
   `perf stat` cpu_core cycles + wall time as the baseline.
2. Run the existing Delphy test suite (`build/release/tests/tests`).
   
3. **After profile**: re-profile on the same dataset and workload,
   compare cpu_core cycles and wall time to the baseline.  Expect
   `slow_elements` to drop from ~25.6% to near zero.

## Expected impact

Eliminating ~600K coroutine suspend/resume cycles per tree traversal
should remove most of the 25.6% overhead.  The replacement loops do
the same arithmetic with plain `for` loops and a tiny fixup pass.
The `from_states` fixup is negligible (typically 0--2 entries per
node, since only a small fraction of sites are missing and mutations
at those sites are rare).

## Future: O(num_intervals) processing via prefix sums

The bulk loops proposed above are O(num_sites) per node -- they
iterate over every site in every interval.  For the
accumulate-by-state sites (sites 1, 2, 5, 6, 7 — which accumulate a
per-site value like `nu_l[l]` grouped by `ref_seq[l]` and optionally
by partition), we could reduce this to
O(num_intervals) by precomputing prefix sums.

Specifically, if we precomputed:

- `cum_nu_l_by_state[a][l]` = sum of `nu_l[i]` for `i < l` where
  `ref_seq[i] == a`

then for an interval `[start, end)` the total contribution to
state `a` would be `cum_nu_l_by_state[a][end] -
cum_nu_l_by_state[a][start]`, computed in O(1) per interval.

For the partitioned variants (sites 1, 2, 7), a similar prefix sum
partitioned by both state and partition would be needed.

These prefix sums depend only on `ref_sequence`, `nu_l`, and
`partition_for_site`, all of which change only when the evo model
changes (not per-node).  They could be computed once per
`recalc_derived_quantities()` call and reused across all nodes.

This optimization does **not** apply to the per-site indexed write
sites (3: `calc_Ttwiddle_l`, 4: `calc_T_l_a`), which write to a
per-site output array and must touch each site individually
regardless.

This is not worth doing now -- the coroutine elimination alone should
remove the bulk of the overhead -- but is worth keeping in mind for
datasets with very large gaps (e.g., thousands of consecutive missing
sites).
