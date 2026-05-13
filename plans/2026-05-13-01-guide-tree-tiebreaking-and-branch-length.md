# Guide Tree Tie-Breaking and Branch Length Fix

Follow-up to the [Round 2 plan](2026-05-07-02-better-tree-init-round2-rooting-and-conversion.md),
addressing two issues observed during manual testing (see "Observations from manual testing"
in that plan).

## Problem

1. **Ladder-like artifacts**: When many tips have identical or very similar sequences (e.g.,
   SARS-CoV-2 [Lemieux et al, 2021], Ebola [Gire et al, 2014]), `find_best_attachment_arc`
   in `Utree_builder` finds many edges with the same best cost. It always picks the first
   one encountered during the search, which produces degenerate ladder topologies — each new
   tip attaches next to the previous one, creating a long chain of degree-3 nodes instead of
   a balanced subtree.

2. **Overly conservative branch lengths**: The `min_branch_length` in `utree_to_phylo_tree`
   is `max(0.5, 0.5 / lambda)`, which can be half a day or more. This pushes inner node
   times unnecessarily far into the past, especially for closely related sequences where
   the true internal branch lengths are very short.

## Fix

### Random tie-breaking among equal-cost arcs

In `find_best_attachment_arc` (`core/utree.cpp`), instead of tracking a single `best_arc`,
accumulate all arcs with equal best cost into a `best_arcs_` vector. When a strictly better
cost is found, clear the vector. At the end, choose uniformly at random from the vector
using `absl::Uniform`. This breaks the systematic bias toward the first-found arc.

Add `best_arcs_` as a member of `Utree_builder` (reused across calls, like `pq_`).

### Reduce minimum branch length

In `utree_to_phylo_tree` (`core/utree.cpp`), change `min_branch_length` from
`max(0.5, 0.5 / lambda)` to a fixed `0.1` days. The old value was overly conservative;
0.1 days (~2.4 hours) is short enough to avoid distorting the tree while still preventing
zero-length branches.

## Files to modify

| File | What |
|------|------|
| `core/utree.cpp` | Change `find_best_attachment_arc` to collect equal-cost arcs and pick randomly; add `best_arcs_` member; change `min_branch_length` to `0.1` |

## Verification

1. Existing tests pass: `cmake --build build/debug --target tests -j$(nproc) && build/debug/tests/tests`
2. Manual: run `delphy_ui --v0-init mp-plus-timing` on SARS-CoV-2 (Lemieux et al, 2021)
   and Ebola (Gire et al, 2014) and verify the ladder artifacts are reduced.
