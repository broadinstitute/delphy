# Better Tree Init, Round 4: Root-to-Tip Regression Rooting

This is the fourth in a series of plans implementing a new tree initialization pipeline for
Delphy.  See the [Round 1 plan](2026-05-04-01-better-tree-init-round1-utree-and-guide-tree.md),
section "Roadmap", for the full multi-round plan.

## Background

Round 2 implemented midpoint rooting: find the tree's diameter, place the root at the
(timed) midpoint of a diametral path, then estimate the mutation rate via OLS regression of
root-to-tip mutation counts against tip dates.  This works well when the molecular clock
holds approximately, but midpoint rooting is purely topological -- it doesn't use temporal
information to *choose* the root position.  Even under a strict molecular clock, midpoint
rooting is only strictly correct when all tip dates are equal; with unequal tip dates, it is
a reasonable heuristic but can be seriously wrong.  More generally, for datasets where the
diametral path doesn't pass through the true root (e.g., rate variation, long stems, or
asymmetric tree shapes), midpoint rooting can place the root far from the correct position.

## Goal for Round 4

Replace midpoint rooting with root-to-tip regression rooting.  The idea is to try every
possible root position along every edge of the unrooted tree.  At each candidate position,
count the mutations from the root to every tip and regress those counts against tip dates.
The position whose regression fits best (highest R^2, subject to the fit implying a positive
mutation rate) is chosen as the root.

The key insight is that the OLS regression quantities for each candidate root can be
expressed purely in terms of per-arc subtree statistics, which are computed incrementally
by two DFS passes over the tree.  Root positions are always placed by splitting an edge,
producing a degree-2 root node consistent with bifurcating Phylo_trees.

## Mathematical framework

### The regression model

Given a rooted tree with root R, we can count the number of mutations m_i on the path from R
to each tip i, and we know each tip's sampling date t_i.  A strict molecular clock predicts
a linear relationship:

    m_i = lambda * t_i + c + noise

where lambda is the mutation rate and c = -lambda * t_MRCA is the intercept.  Equivalently,

    t_MRCA = mean_t - mean_m / lambda

Our goal is to find the root position that makes this linear model fit best (highest R^2),
and then read off lambda and t_MRCA from the OLS fit.

**Caveats:** OLS assumes independent, equal-variance Gaussian noise terms for each tip.
Neither assumption holds here.  The {m_i} are correlated because tips sharing phylogenetic
ancestry have correlated root-to-tip distances, and even under tip-wise independence, the
noise variance scales as ~sqrt(lambda * t_i) with finite-size effects for small m_i.
Consequently, the R^2, confidence intervals, and p-values of the regression are not valid
statistical estimates.  However, OLS root-to-tip regression is empirically a good starting
approximation for root placement and rate estimation, and is the standard approach used by
[TempEst](https://github.com/beast-dev/Tempest)
(Rambaut et al., *Virus Evolution* 2016, doi:10.1093/ve/vew007) and
[Clockor2](https://clockor2.github.io/)
(Featherstone et al., *Systematic Biology* 2024, doi:10.1093/sysbio/syae003).

### Notation and global constants

We have an unrooted tree with N tips.  Each tip's date is the midpoint of its date range:
t_i = (tip_descs[i].t_min + tip_descs[i].t_max) / 2.  Tip dates are always used as centered
deviations: dt_i = t_i - mean_t, where mean_t = (1/N) * sum(t_i).  All formulas below use
dt_i, never raw t_i, to avoid catastrophic cancellation.

Global constants (independent of root position):

    N       = number of tips
    mean_t  = (1/N) * sum(t_i)
    Var_t   = (1/N) * sum(dt_i^2)

If Var_t = 0 (all tips have the same date), skip regression entirely and fall back to
midpoint rooting.

### Per-arc subtree statistics

Evaluating R^2 at every candidate root naively would require an O(N) traversal per candidate.
Instead, we define per-arc subtree statistics that let us compute the OLS quantities at any
candidate position in O(1), after two O(N) DFS passes to populate them.

For each arc `a` in the tree, define the **subtree beyond `a`** = Sub(a) = the connected
component containing target(a) after removing the edge.  We associate with each arc an
`Arc_stats` struct containing five fields, all with distances (mutation counts) measured from
target(a):

    n(a)        = number of tips in Sub(a)
    sum_dt(a)   = sum of dt_i for tips in Sub(a)
    sum_m(a)    = sum of dist(target(a), i) for tips in Sub(a)
    sum_m_dt(a) = sum of dist(target(a), i) * dt_i for tips in Sub(a)
    sum_m2(a)   = sum of dist(target(a), i)^2 for tips in Sub(a)

Here `m` stands for "mutations from target to tip" -- the same quantity that becomes the
root-to-tip mutation count when the arc's target is the root.

### Two operations on subtree stats

All recurrences reduce to two operations on Arc_stats structs:

**prepend(s, D)**: Given stats `s` measured from some node Q, and an edge from X to Q with
D deltas, return stats measured from X.  Each tip's mutation count increases by D:

    prepend(s, D).n        = s.n
    prepend(s, D).sum_dt   = s.sum_dt
    prepend(s, D).sum_m    = D * s.n      + s.sum_m
    prepend(s, D).sum_m_dt = D * s.sum_dt + s.sum_m_dt
    prepend(s, D).sum_m2   = D^2 * s.n    + 2*D * s.sum_m + s.sum_m2

**combine(A, B)**: Given two sets of stats measured from the same node X (disjoint tip
sets), return their union.  Element-wise addition of all five fields.

### OLS quantities at a candidate root position

Consider the edge between nodes A and B, with arcs A_to_B and B_to_A = mate(A_to_B).
Let D = ssize(arcs[A_to_B].deltas).  A candidate root at position k in {0, 1, ..., D}
splits the edge so that k deltas lie between A and the root and D-k deltas lie between
the root and B.  The root-to-tip mutation counts are:

    For tip i in Sub(A_to_B) (B's side): m_i(k) = (D - k) + dist(B, i)
    For tip i in Sub(B_to_A) (A's side): m_i(k) = k       + dist(A, i)

Equivalently, the summary statistics for all tips, measured from the root at position k,
are simply:

    root_stats(k) = combine(prepend(stats(A_to_B), D - k), prepend(stats(B_to_A), k))

This is the central formula.  The OLS quantities read off directly:

    mean_m    = root_stats.sum_m / N
    mean_m_dt = root_stats.sum_m_dt / N
    mean_m2   = root_stats.sum_m2 / N

    Cov(m, t) = mean_m_dt                  [since sum(dt_i) = 0]
    Var(m)    = mean_m2 - mean_m^2

A candidate position is valid only if Cov(m, t) > 0 (positive molecular clock signal).
From a valid position, the OLS fit gives:

    lambda = Cov(m, t) / Var(t)
    t_MRCA = mean_t - mean_m / lambda

The predicted mutation counts are m_hat_i = mean_m + lambda * dt_i, and the regression fit
quality is measured by R^2, the fraction of variance in m explained by the linear fit:

    R^2 = 1 - SS_res / SS_tot = Cov(m, t)^2 / (Var(m) * Var(t))

where SS_res = sum((m_i - m_hat_i)^2) and SS_tot = sum((m_i - mean_m)^2) = N * Var(m).
For reference, the residual mean square RMS = SS_res / (N - 2) is the unbiased estimator
of the noise variance.  TempEst and Clockor2 can optionally select the best root by
minimizing RMS instead of maximizing R^2; in practice the difference is minimal.  We use
R^2 here.  Note that R^2, lambda, and t_MRCA are all ratios of covariances and variances,
so any finite-sample-size correction (1/N vs 1/(N-1)) cancels out and the choice is
immaterial.

## Algorithm

### Pass 0: Global constants

Loop over all tips to compute mean_t (first pass) and Var_t = (1/N) * sum(dt_i^2)
(second pass, using dt_i = t_i - mean_t).  If Var_t = 0, fall back to midpoint rooting
immediately.

### Pass 1: Bottom-up (post-order) DFS

Compute stats for every arc pointing **away from** the DFS root F (i.e., from parent toward
child in the DFS tree).

Pick any node F as the DFS root (e.g., tip 0).  Run `annotated_arc_euler_tour(tree, F)`.
On each **leaving** event, the tour yields the arc X->P that we are backtracking along
(X is the node being left, P is its parent toward F).  At this point, we have finished
visiting X's entire subtree, so we can compute stats for the mate arc P->X (pointing
away from F, into Sub(P->X)):

Base case -- X is a tip:

    stats(P->X) = {n: 1, sum_dt: dt_X, sum_m: 0, sum_m_dt: 0, sum_m2: 0}

Recursive case -- X is an inner node with children C_1, ..., C_d (arcs from X away from F),
with D_1, ..., D_d deltas respectively:

    stats(P->X) = combine(prepend(stats(X->C_1), D_1),
                          ...
                          prepend(stats(X->C_d), D_d))

(In a bifurcating tree, d = 2.)  Each stats(X->C_i) was set on an earlier leaving event
(post-order guarantees children before parent).

### Pass 2: Top-down (pre-order) DFS + R^2 evaluation

Compute stats for every arc pointing **toward** F (i.e., from child toward parent in the
DFS tree).  Evaluate R^2 at every candidate root position as we go.

Run `annotated_arc_euler_tour(tree, F)`.  On each **entering** event, the tour yields the
arc P->X that we are descending along (P is the parent toward F, X is the node being
entered).  Compute stats for the mate arc X->P (pointing toward F, into Sub(X->P)):

Combine the prepended stats of all outgoing arcs from P **other than** P->X.  All such arcs
are already annotated: arcs pointing away from F were set in Pass 1; arcs pointing toward F
were set earlier in this pass (pre-order guarantees parent before child).

    stats(X->P) = combine over all outgoing arcs `a` from P, a != P->X:
                      prepend(stats(a), ssize(arcs[a].deltas))

Special case: if P is a tip (degree 1, so P has no arcs other than P->X), it has no other
arcs to combine.  P itself is a tip, so add its own contribution:

    stats(X->P) = {n: 1, sum_dt: dt_P, sum_m: 0, sum_m_dt: 0, sum_m2: 0}

In practice, this only happens when P = F and F is a tip.

**R^2 evaluation:**  At this point, both stats(P->X) and stats(X->P) are known.  Iterate
k from 0 to D (where D = ssize(arcs[P->X].deltas)):

    root_stats = combine(prepend(stats(X->P), k),
                         prepend(stats(P->X), D - k))

    Cov_mt = root_stats.sum_m_dt / N
    if Cov_mt <= 0: skip    [no valid clock signal at this position]

    Var_m = root_stats.sum_m2 / N - (root_stats.sum_m / N)^2
    if Var_m <= 0: skip    [all tips equidistant from this root — no variance to explain]
    R2 = Cov_mt^2 / (Var_m * Var_t)

    Track the best R^2 with tie-breaking: maintain a vector of (arc, k) pairs that share
    the highest R^2 so far.  If R2 exceeds the current best, clear the vector and start
    fresh; if R2 equals the current best, append.  This follows the same pattern as
    `find_best_attachment_arc` in `Utree_builder`.

**Total cost:**  Two DFS passes over the tree, plus iterating over all k values within each
edge.  The total number of k values across all edges equals the total number of site deltas
in the tree (each delta is on exactly one edge), so the k-iteration cost is O(total deltas).
The overall cost is O(N + total deltas).

## Root selection and rate estimation

After the two passes, pick a winner uniformly at random from the vector of equally-best
(arc, k) pairs.  Call the winner `best_arc` and `best_k`.  Root the tree:

1. Allocate the root node R at index `num_tips + num_inner_nodes_so_far`.  Increment
   `num_inner_nodes_so_far`.

2. Call `split_edge(best_arc, R, site_delta_side)`, where the callback counts deltas and
   assigns the first `best_k` (in iteration order) to the origin side and the rest to the
   target side.  This is the same approach as `midpoint_root_utree`.

3. Compute the rate estimate from root_stats at the best position:

       Cov_mt = root_stats.sum_m_dt / N
       lambda = Cov_mt / Var_t
       mean_m = root_stats.sum_m / N
       t_MRCA = mean_t - mean_m / lambda

   Return the root node R, the best R^2, lambda, and t_MRCA.

## Fallback to midpoint rooting

The regression is skipped entirely and we fall back to `midpoint_root_utree` when:

- Var_t = 0 (all tips have the same date -- detected in Pass 0).
- N <= 2 (not enough data for meaningful regression).
- No candidate position across any edge has Cov(m, t) > 0 (no position gives a positive
  clock signal).

## Integration

### Rooting result type

Replace `Rate_estimate` with a `Rooting_info` struct that bundles root placement, method,
and rate estimation results:

```cpp
enum class Rooting_method { regression, midpoint };

struct Rooting_info {
  Node_index root;
  Rooting_method method;
  double r2;       // R^2 of root-to-tip regression (may be <= 0 if no clock signal)
  double lambda;   // overall mutations per day (not per-site)
  double t_MRCA;   // estimated root date (days since epoch)
};
```

### New function

```cpp
auto regression_root_utree(Utree& tree, const std::vector<Tip_desc>& tip_descs,
                           absl::BitGenRef bitgen)
    -> Rooting_info;
```

Finds the best root position by maximizing R^2 of root-to-tip regression against tip dates.
Roots the tree in place and returns a `Rooting_info` with the root node, rate estimate, and
R^2.  The `bitgen` is used to break ties when multiple candidate positions share the highest
R^2.  Falls back to `midpoint_root_utree` when regression is not applicable (see Fallback
section above).

### Updated midpoint rooting

Merge `estimate_rate_and_root_date` into `midpoint_root_utree` so that it also returns
`Rooting_info` directly (computing lambda, t_MRCA, and R^2 after placing the root).
Remove `estimate_rate_and_root_date`.

```cpp
auto midpoint_root_utree(Utree& tree, const std::vector<Tip_desc>& tip_descs)
    -> Rooting_info;
```

### Updated conversion function

`utree_to_phylo_tree` takes `const Rooting_info&` instead of separate `Node_index root`
and `Rate_estimate` parameters:

```cpp
auto utree_to_phylo_tree(
    Utree& utree, const Rooting_info& rooting_info,
    const std::vector<Tip_desc>& tip_descs,
    absl::BitGenRef bitgen) -> Phylo_tree;
```

### Pipeline change

Add a `rooting_hook` callback to `build_initial_phylo_tree` so callers can observe the
rooting result without `utree.cpp` doing any I/O itself:

```cpp
auto build_initial_phylo_tree(
    Real_sequence ref_sequence, std::vector<Tip_desc> tip_descs,
    absl::BitGenRef bitgen,
    const std::function<void(int,int)>& guide_tree_progress_hook = [](int,int){},
    const std::function<void(int,int,int)>& refined_tree_progress_hook = [](int,int,int){},
    const std::function<void(const Rooting_info&)>& rooting_hook = [](const Rooting_info&){})
    -> Phylo_tree;
```

Replace the current midpoint + estimate_rate sequence:
```cpp
auto rooting_info = regression_root_utree(utree, tip_descs, bitgen);
rooting_hook(rooting_info);
auto tree = utree_to_phylo_tree(utree, rooting_info, tip_descs, bitgen);
```

Thread the `rooting_hook` through `build_rough_initial_tree_from_maple` in `cmdline.h`/
`cmdline.cpp` (same defaulted parameter pattern).  At the call site in `cmdline.cpp`,
supply a hook that formats and prints the rooting diagnostics to stderr:

```cpp
[&](const Rooting_info& rooting_info) {
  auto method_str = (rooting_info.method == Rooting_method::regression)
      ? "root-to-tip regression" : "midpoint";
  std::cerr << absl::StreamFormat("- rooting: %s, R^2 = %.4f, ...", method_str, rooting_info.r2)
            << ...;
}
```

This keeps all I/O out of `utree.cpp`.

## Files to modify

| File | Action | What |
|------|--------|------|
| `core/utree.h` | Modify | `Rooting_info`, `Rooting_method`; declare `regression_root_utree`; update signatures |
| `core/utree.cpp` | Modify | Implement `regression_root_utree`; merge `estimate_rate_and_root_date` into `midpoint_root_utree`; update `build_initial_phylo_tree`; remove I/O |
| `core/cmdline.h` | Modify | Thread `rooting_hook` through `build_rough_initial_tree_from_maple` |
| `core/cmdline.cpp` | Modify | Forward `rooting_hook`; supply formatting hook at call site |
| `tests/utree_tests.cpp` | Modify | Add regression rooting tests; update midpoint and conversion tests for `Rooting_info` |

## Verification

1. **Unit tests on small trees** (3-5 tips with known structure):
   - Construct a tree where the correct root is known.  Verify the regression picks the
     right edge and a reasonable k.
   - Verify the rate estimate matches a hand-computed OLS.

2. **Degenerate cases**:
   - All-same-date tips: verify fallback to midpoint rooting.
   - 2-tip tree: verify fallback.
   - Star topology (all tips equidistant from center): verify a sensible root is chosen.

3. **Comparison with midpoint rooting**:
   - On datasets where midpoint rooting is known to be poor (e.g., H5N1 Andersen), verify
     that regression rooting produces a more reasonable tMRCA.

4. **End-to-end**: Run `delphy_ui` with `--v0-init mp-plus-timing` on standard demo datasets.
   Verify `assert_phylo_tree_matches_tip_descs` passes and results are sensible.

5. **Build**: `cmake --build build/debug --target tests -j$(nproc) && build/debug/tests/tests --gtest_filter='Utree*'`

---

## Appendix A: Closed-form optimal k within an edge

The algorithm above iterates k from 0 to D for each edge, giving O(total deltas) total work.
A closed-form solution exists that finds the optimal k per edge in O(1), reducing the total
to O(N).  We prefer the explicit iteration for clarity and because it makes filtering out
positions with Cov(m,t) <= 0 transparent.  In practice, D is small on almost all arcs, so the
iteration is fast; the only scenario where O(total deltas) matters is trees with very long
branches carrying many deltas.

For arc A_to_B with D deltas, let B_to_A = mate(A_to_B).  The OLS quantities at position k
are polynomial in k:

    sum_m(k)    = (D-k) * n(A_to_B) + sum_m(A_to_B)
                +   k   * n(B_to_A) + sum_m(B_to_A)                              -- linear
    sum_m_dt(k) = (D-k) * sum_dt(A_to_B) + sum_m_dt(A_to_B)
                +   k   * sum_dt(B_to_A) + sum_m_dt(B_to_A)                      -- linear
    sum_m2(k)   = (D-k)^2 * n(A_to_B) + 2*(D-k) * sum_m(A_to_B) + sum_m2(A_to_B)
                +  k^2    * n(B_to_A) +   2*k   * sum_m(B_to_A) + sum_m2(B_to_A) -- quadratic

Write sum_m_dt(k) = p0 + p1*k (linear) and define:

    q(k) = N * sum_m2(k) - sum_m(k)^2    [= N^2 * Var(m)]

which is quadratic: q(k) = q0 + q1*k + q2*k^2.

R^2(k) is proportional to (p0 + p1*k)^2 / (q0 + q1*k + q2*k^2).  Setting the derivative
to zero and factoring out (p0 + p1*k):

    2*p1*(q0 + q1*k + q2*k^2) = (p0 + p1*k)*(q1 + 2*q2*k)

The k^2 terms (2*p1*q2*k^2) cancel on both sides, leaving a linear equation:

    k* = (2*p1*q0 - p0*q1) / (2*p0*q2 - p1*q1)

The best integer k in [0, D] is among floor(k*), ceil(k*), 0, and D (evaluating R^2 at up
to 4 candidates).  This gives O(1) per edge, O(N) total.
