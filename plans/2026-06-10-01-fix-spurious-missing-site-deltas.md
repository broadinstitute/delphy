# Fix: spurious missing-site deltas on tip edges

## Problem

In `Utree_builder::attach_tip_to_focal_arc`, when splitting an edge to insert a new tip X,
each delta on the split edge is randomly assigned to the A-M or M-B side.  For deltas at
sites missing at X, the code was calling `pop_front_site_deltas` on `M_to_X_deltas_` when
the delta went to the A side.  But `focus_to_X_deltas_` (which seeds `M_to_X_deltas_`)
excludes missing sites, so `pop_front_site_deltas` inserts a *new* entry instead of
adjusting an existing one.  This creates a spurious delta on the M-X edge at a site where X
has no data.

**Consequence:** `eval_focal_arc` correctly excludes missing sites, so it underpredicts the
actual M-X delta count.  The discrepancy is small per-tip (0 or 1) but accumulates over the
tree, inflating `count_deltas()`.  This can cause slightly suboptimal placement decisions
during tree construction (the cost comparison between candidate edges is biased by the
spurious deltas already present in the tree).

## Fix

For missing-site deltas on the split edge, randomly assign the side without touching
`M_to_X_deltas_`:

```cpp
if (miss_X.contains(sd.site) || tree_.globally_missing_sites.contains(sd.site)) {
    return std::bernoulli_distribution{0.5}(bitgen_) ? A : B;
}
```

The M-X edge should never carry deltas at sites where X has no data.

## Files changed

| File | Action | What |
|------|--------|------|
| `core/utree.cpp` | Modify | Remove `pop_front_site_deltas` call for missing-site deltas in `attach_tip_to_focal_arc`; add integrity and tip-matching checks to `Utree_builder::finish()`; add missing-site delta check to `assert_utree_matches_tip_descs` |
| `tests/utree_tests.cpp` | Modify | Remove `assert_matches_tip_descs_with_missing_data` test (it expected the buggy behavior) |

## Verification

1. Temporarily revert the fix in `attach_tip_to_focal_arc` and confirm the new CHECK in
   `assert_utree_matches_tip_descs` fires on datasets with missing data (Dudas, Zika,
   O'Toole, Parker, H5N1, RSV).
2. Re-apply the fix and confirm the CHECK passes on all datasets.
3. Run all Utree tests.

## Preliminary results

Total parsimony score (site deltas) for trees built with `--v0-init mp-plus-timing`.
Lower is better.  "Before" values are from the Round 3 plan
(`2026-05-15-01-better-tree-init-round3-guide-tree-reordering.md`, "Refined tree" column),
which builds and iteratively refines the tree without this fix.  "After" values are with
this fix applied.  Both guide tree and refined tree are affected: any tip inserted after a
prior tip with missing data can pick up spurious deltas when `split_edge` distributes an
existing missing-site delta.

| Dataset                          |    N    |    L    | Before (Round 3) |  After |
|----------------------------------|---------|---------|------------------|--------|
| Ebola (Dudas et al, 2017) (*)    |   1,610 |  18,996 |            2,814 |  2,808 |
| Ebola (Gire et al, 2014)         |      81 |  18,959 |               55 |     55 |
| H3N2 (Rambaut et al, 2008)       |     165 |   1,698 |              394 |    392 |
| H5N1 (Andersen lab, 2025) (*)    |   3,339 |  13,136 |            4,392 |  4,362 |
| mpox (O'Toole et al, 2023)       |      41 | 197,209 |              347 |    343 |
| mpox (Parker et al, 2025)        |     177 | 197,209 |              946 |    945 |
| SARS-CoV-2 (Lemieux et al, 2021) |     757 |  29,903 |              665 |    664 |
| Zika (Metsky et al, 2017)        |     174 |  10,807 |            1,300 |  1,244 |
| RSV-A (NCBI, 2012--2025) (**)    |   3,000 |  15,277 |           38,540 | 38,133 |

(*) Note: Substantial run-to-run variability observed for these datasets
(**) Note: The RSV-A dataset is contaminated by 3 nearby RSV-B sequences that are ~2,500
mutations from the RSV-A sequences.
