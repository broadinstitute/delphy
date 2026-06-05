# Better Tree Init, Round 4tris: Marginal Likelihood Rooting and Timing

This is another in a series of plans implementing a new tree initialization pipeline for
Delphy.  See the [Round 1 plan](2026-05-04-01-better-tree-init-round1-utree-and-guide-tree.md),
section "Roadmap", for the full multi-round plan.

## Background

[Round 4](2026-05-18-01-better-tree-init-round4-regression-rooting.md) implemented OLS
root-to-tip regression rooting: try every possible root position, regress root-to-tip
mutation counts against tip dates via ordinary least squares (OLS), and pick the position with
the best R^2.  This works well empirically but has two theoretical limitations: (1) OLS
treats root-to-tip mutation counts as independent, ignoring the phylogenetic correlations
between tips that share ancestry; and (2) OLS assumes equal variance across tips, when in
reality the variance scales with the expected number of mutations.

[Round 4bis](2026-05-21-01-better-tree-init-round4bis-gls-regression-rooting.md) attempted
to address this via generalized least squares (GLS), following Richard Neher's
[`TreeRegression`](https://github.com/neherlab/treetime/blob/master/treetime/treeregression.py)
framework in TreeTime.  GLS accounts for both phylogenetic correlation and tip-specific
variance by using the tree's covariance structure as a weight matrix.  However, the GLS
approach collapses the effective sample size in practice (see Round 4bis, "Post-implementation
findings"), and Delphy continues to use OLS for rooting.

This plan implements a different approach based on a Gaussian branch-length model.  After
an email exchange with Richard Neher on GLS vs OLS, he pointed us to an unpublished note
of his on
["Detecting outliers in heterochronous phylogenetic trees"](https://github.com/neherlab/outlier-detection/blob/master/outliers.pdf),
which describes a related scheme.  Rather than regressing root-to-tip distances, we model
each branch length individually as having a Gamma posterior (from a flat prior plus observed
Poisson mutations), approximated as a Gaussian.  The key idea is to **progressively integrate out node times from tips to
root**, obtaining a reduced potential that depends only on the root time and the mutation
rate.  When tip dates are fixed, the reduced potential's dependence on the mutation rate
has a simple closed form, enabling joint O(N)-time optimization of root position and
mutation rate.

This Gaussian model is related to the GLS model (both use Gaussian distributions for
mutations on branches), but different enough that we hope it avoids the effective-sample-size
collapse observed with GLS.  We are working on a separate write-up exploring these ideas
in more depth.

## Goal for Round 4tris

Add a new rooting function `gaussian_root_utree` that jointly optimizes the root position
and mutation rate by maximizing the marginal likelihood of a Gaussian branch-length model.
The function returns not just a root position and rate estimate, but also optimal times for
all inner nodes, enabling `utree_to_phylo_tree` to use them directly instead of computing
node times from a linear model.

## Mathematical framework

### Overview: integrating out node times

Given a rooted tree with mutation rate mu, we assign a time tau_i to each node.  The
Gaussian approximation to the branch-length posteriors gives a joint potential (negative
log-likelihood) over the inner node times:

    U = sum_i ((tau_i - tau_{p_i}) - (d_i+1)/mu)^2 / (2 * (d_i+1)/mu^2)

where the sum is over all non-root nodes i with parent p_i and d_i observed mutations.
(We use d_i+1 rather than d_i in the numerator because the Gamma(d_i+1, mu) posterior
arises from a flat prior on branch length plus d_i observed Poisson mutations.  This means
even a zero-mutation branch has positive expected length 1/mu.  We also note that this
potential neglects any ancestry prior such as a coalescent.)

Our strategy is to **progressively integrate out each inner node time** in postorder
(tips to root), using completion of the square at each step.  Because U is a sum of
Gaussians, integrating out a node time produces another Gaussian, and the reduced potential
retains the same functional form with updated parameters.  When we reach the root node R,
the fully reduced potential depends only on tau_R and mu:

    U_reduced = (tau_R - tau_tilde)^2 / (2 * sigma^2) + Delta

This is minimized at tau_R = tau_tilde, with minimum value Delta.  To find the best root,
we pick the position with the lowest Delta.  Each candidate root position also yields an
optimal mutation rate (see below), enabling joint root-and-rate optimization.

### Per-arc quantities

We associate three quantities with each directed arc P -> X, summarizing the net effect
of integrating out all node times in the subtree rooted at X away from P:

- **tau_tilde_{P->X}**: the "neutral" time for node X, considering only constraints from
  the subtree below.
- **sigma^2_{P->X}**: the variance of this neutral time.
- **Delta_{P->X}**: the residual potential accumulated from integrating out the subtree's
  node times.

These quantities depend only on the subtree of X away from P: the label P identifies
*which* subtree, not properties of P or the P-X edge.

### Recursive computation via "shift" and "combine"

The per-arc quantities are computed by two half-steps at each inner node X with parent P
and children C_1, ..., C_N.  After integrating out everything below X, the reduced
potential has the form (other terms involving tau_X from the P direction) +
(tau_X - tau_tilde)^2 / (2*sigma^2) + Delta, where the last two terms are parametrized
by precisely the three per-arc quantities for the arc P -> X.

**Shift half-step** (one per child): integrates out tau_{C_n} by completing the square
on the joint quadratic from the X-C_n branch and the C_n subtree.  Result: a harmonic
well on tau_X:

    sigma^2_{X|C_n} = sigma^2_{X->C_n} + (d+1)/mu^2            (variances add)
    tau_tilde_{X|C_n} = tau_tilde_{X->C_n} - (d+1)/mu           (shift by branch mean)
    Delta_{X|C_n} = Delta_{X->C_n}                               (log normalization terms dropped)

where d = d_{X,C_n} is the mutation count on the X-C_n branch.  (The exact shift formula
for Delta includes a log normalization term from the Gaussian integration; we drop it as
negligible compared to the quadratic residuals that dominate Delta.)

**Combine half-step** (once per node): merges the harmonic wells from all children via
completion-of-the-square into a single well:

    sigma^2_{P->X} = [sum_n 1/sigma^2_{X|C_n}]^{-1}
    tau_tilde_{P->X} = sum_n (sigma^2_{P->X} / sigma^2_{X|C_n}) * tau_tilde_{X|C_n}
    Delta_{P->X} = sum_n [(tau_tilde_{P->X} - tau_tilde_{X|C_n})^2 / (2*sigma^2_{X|C_n}) + Delta_{X|C_n}]

**Tip base case**: tau_tilde = t_tip, sigma^2 = 0, Delta = 0 (tip date is fixed, no
integration needed).

### A-F coefficient decomposition: making mu-dependence explicit

When all tip dates are fixed (sigma^2_tip = 0), the per-arc quantities have a simple
mu-dependence:

    tau_tilde_{P->X} = A_{P->X} + B_{P->X} / mu
    sigma^2_{P->X}   = C_{P->X} / mu^2
    Delta_{P->X}     = D_{P->X} + F_{P->X} * mu + G_{P->X} * mu^2

This works because all sigma^2 values scale as 1/mu^2, so the weights in the combine step
(ratios of sigma^2 values) are mu-independent.  This cancellation fails when tip dates have
nonzero uncertainty.

(The exact Delta would also include an E * ln(mu) term from log normalization factors in
the shift step, where each of the N-2 non-root inner-node integrations contributes one unit to E.
Since E is the same for every candidate root position, dropping it does
not affect root selection.  It does shift the optimal mu from -F/(2G) to a root of
2G*mu^2 + F*mu + E = 0, introducing a small O(1/N)-relative bias in the rate
estimate, which is acceptable for initialization.)

**Tip base case** (using centered dates, dt_X = t_X - mean_t):

    A = dt_X,   B = C = D = F = G = 0

**Shift half-step** (branch with d mutations, child coefficients X->C_n -> "X|C_n"):

    A_{X|C_n} = A_{X->C_n}
    B_{X|C_n} = B_{X->C_n} - (d + 1)
    C_{X|C_n} = C_{X->C_n} + (d + 1)
    D_{X|C_n} = D_{X->C_n}
    F_{X|C_n} = F_{X->C_n}
    G_{X|C_n} = G_{X->C_n}

No special case is needed when the child is a tip (C_{X->C_n} = 0): the formulas work
uniformly.

**Combine half-step** (children X|C_1, ..., X|C_N -> parent arc P->X):

    C_{P->X} = [sum_n 1/C_{X|C_n}]^{-1}                                   (compute first)
    A_{P->X} = sum_n (C_{P->X}/C_{X|C_n}) * A_{X|C_n}                     (weighted average)
    B_{P->X} = sum_n (C_{P->X}/C_{X|C_n}) * B_{X|C_n}                     (weighted average)

    D_{P->X} = sum_n [D_{X|C_n} + (B_{P->X} - B_{X|C_n})^2 / (2 * C_{X|C_n})]
    F_{P->X} = sum_n [F_{X|C_n} + (A_{P->X} - A_{X|C_n})*(B_{P->X} - B_{X|C_n}) / C_{X|C_n}]
    G_{P->X} = sum_n [G_{X|C_n} + (A_{P->X} - A_{X|C_n})^2 / (2 * C_{X|C_n})]

Note: D, F, G formulas use the already-computed A_{P->X} and B_{P->X}.

### Optimal mu and root selection

At a candidate root R, the fully reduced potential is:

    Delta = D + F * mu + G * mu^2

Setting d(Delta)/d(mu) = F + 2*G*mu = 0 gives:

    mu = -F / (2*G)      (valid only when F < 0, i.e., positive clock signal)

and the residual potential at the optimum is:

    Delta_min = D - F^2 / (4*G)

We select the root position with the lowest Delta_min, filtering out any candidate with
F >= 0 (no clock signal).

### Signs of the coefficients

The recursion relations imply that the B coefficients are negative, while the C and G
coefficients are non-negative.  F < 0 is the condition for a positive clock signal.
The A coefficients are bounded by the range of (centered) tip dates in the subtree.

## Algorithm

### Pass 0: Global setup

Loop over all tips to compute mean_t and verify that not all dates are identical (if they
are, fall back to midpoint rooting -- there is no temporal signal to exploit).  Store
centered dates dt_i = t_i - mean_t for each tip; all tip dates entering the algorithm use
these centered values.  Results (node times) are converted back to absolute dates
(+ mean_t) on output.

We also maintain a separate per-arc array `max_node_time[a]`, which tracks the latest
allowable time for target(a) such that all branches in the subtree beyond `a` have length
at least `min_branch_length` (= 0.1 days).  This is accumulated during Passes 1 and 2, in a separate
array from `Gaussian_stats` (see below).

### Data structure

```cpp
struct Gaussian_stats {
  double A = 0.0;   // tau_tilde = A + B/mu  (using centered dates)
  double B = 0.0;
  double C = 0.0;   // sigma^2 = C/mu^2
  double D = 0.0;   // Delta = D + F*mu + G*mu^2
  double F = 0.0;
  double G = 0.0;
};
```

### Two operations on Gaussian_stats

**shift(s, d)**: Given stats `s` for arc X->C (the subtree below C), and a branch
from X to C with d mutations, return the "X|C" stats (the effect on X of everything
below it in the C direction).

    return {s.A,
            s.B - (d+1),
            s.C + (d+1),
            s.D,
            s.F,
            s.G}

**combine(a, b)**: Given two shifted stats (X|C_1 and X|C_2) at the same node X,
return the combined stats for the parent arc P->X.  This is the completion-of-the-square
merge.

    C = 1.0 / (1.0/a.C + 1.0/b.C)           // = a.C * b.C / (a.C + b.C)
    wA = C / a.C                             // weight for a
    wB = C / b.C                             // weight for b     (wA + wB == 1)
    A = wA * a.A + wB * b.A
    B = wA * a.B + wB * b.B
    dA_a = A - a.A,  dB_a = B - a.B          // residuals for child a
    dA_b = A - b.A,  dB_b = B - b.B          // residuals for child b
    D = (a.D + dB_a*dB_a / (2*a.C)) + (b.D + dB_b*dB_b / (2*b.C))
    F = (a.F + dA_a*dB_a / a.C) + (b.F + dA_b*dB_b / b.C)
    G = (a.G + dA_a*dA_a / (2*a.C)) + (b.G + dA_b*dA_b / (2*b.C))
    return {A, B, C, D, F, G}

Note: `combine` is written for exactly 2 children (binary tree).  For the general
case, compute C, A, B first from all children, then accumulate D, F, G.  Delphy's
Utrees are bifurcating (degree-3 inner nodes), so the 2-child case suffices.

### Pass 1: Bottom-up (post-order) DFS

Pick any node O as the DFS focus.  Compute Gaussian_stats for every arc pointing
**away from** O.  Also accumulate `max_node_time` for every such arc (used for
clamping in Pass 4).

Run `annotated_arc_euler_tour(tree, O)`.
On each **leaving** event for arc X->P, compute `stats[mate(X->P)]` = `stats[P->X]`:

Base case -- X is a tip:

    stats[P->X] = {A: dt_X, B: 0, C: 0, D: 0, F: 0, G: 0}
    max_node_time[P->X] = t_X

where dt_X = t_X - mean_t is the centered tip date.

Recursive case -- X is an inner node with children C_1, C_2 (arcs X->C_1 and X->C_2):

    stats[P->X] = combine(shift(stats[X->C_1], d_1),
                           shift(stats[X->C_2], d_2))
    max_node_time[P->X] = min(max_node_time[X->C_1] - min_branch_length,
                               max_node_time[X->C_2] - min_branch_length)

where d_1, d_2 are the mutation counts on the X-C_1 and X-C_2 branches.  The `- min_branch_length`
at each level ensures that every branch between X and the constraining tip has room for
at least `min_branch_length` days of length.

### Pass 2: Top-down (pre-order) DFS

Compute Gaussian_stats (and max_node_time) for every arc pointing **toward** O.

Run `annotated_arc_euler_tour(tree, O)`.  On each **entering** event for arc P->X,
compute `stats[X->P]`:

If P is a tip (only happens when P = O):

    stats[X->P] = {A: dt_P, B: 0, C: 0, D: 0, F: 0, G: 0}
    max_node_time[X->P] = t_P

Otherwise, P has exactly two other outgoing arcs, a_1 and a_2, with d_1 and d_2
mutations respectively:

    stats[X->P] = combine(shift(stats[a_1], d_1), shift(stats[a_2], d_2))
    max_node_time[X->P] = min(max_node_time[a_1] - min_branch_length, max_node_time[a_2] - min_branch_length)

### Pass 3: Root evaluation

Loop over all edges (arc pairs A->B / B->A), skipping free arc pairs (both sides have
C = 0).  For each edge with D mutations, evaluate every candidate position k = 0, 1, ..., D:

    shifted_A = shift(stats[B->A], k)             // k muts from R to A
    shifted_B = shift(stats[A->B], D - k)         // D-k muts from R to B
    root = combine(shifted_A, shifted_B)

    // Check for valid clock signal
    if root.F >= 0.0 or root.G <= 0.0:  skip      // G <= 0.0 is defensive

    // Residual potential at optimal mu
    delta_min = root.D - root.F * root.F / (4.0 * root.G)

    // Track best (lowest delta_min)
    if delta_min < best_delta:
        best_delta = delta_min
        best_arc = arc
        best_k = k

**Total cost:**  Two DFS passes over the tree (O(N) each), plus iterating over all k
values within each edge.  The total number of k values across all edges equals the total
number of site deltas in the tree, so the k-iteration cost is O(total deltas).  The
overall cost is O(N + total deltas), same as Round 4.

## Root selection, rate estimation, and node timing

After the three passes, recompute `root_stats` at the winning position:

    mu = -root.F / (2.0 * root.G)
    t_MRCA = root.A + root.B / mu + mean_t        // convert back from centered dates
    lambda = mu                                     // mutations per day

Clamp the root time before `split_edge` (which frees `best_arc`).  The original arc
pair's `max_node_time` entries give the latest allowable times for each subtree:

    t_MRCA = min(t_MRCA, min(max_node_time[B->A] - min_branch_length, max_node_time[A->B] - min_branch_length))

This establishes the invariant that ensures the lower bound never exceeds the upper
bound at any node during Pass 4.

Allocate the root node R and split the edge, same as Rounds 4/4bis.  `split_edge`
only invalidates the split edge's arcs; all other per-arc `gaussian_stats` and
`max_node_time` entries remain valid.  Pass 4 uses child-direction stats and
`max_node_time` entries from Passes 1/2 (which are unaffected), and R's time is set
directly, so no stats lookup is needed for R.

### Pass 4: Node timing (preorder from R)

After placing the root, assign optimal times to all inner nodes using a preorder
traversal.  At each inner node X with parent P and children C_1, C_2, the optimal time
is a weighted average of neutral times implied by each direction:

    tau_X = (w_P * (tau_P + (d_{PX}+1)/mu) + w_1 * tau_1 + w_2 * tau_2) / (w_P + w_1 + w_2)

where:
- w_P = 1 / (d_{PX}+1): inverse branch flexibility from parent
- w_n = 1 / C_{X|C_n}: inverse subtree flexibility from child C_n
- tau_n = A_{X|C_n} + B_{X|C_n}/mu + mean_t: neutral time from child direction
  (A_{X|C_n} and B_{X|C_n} come from shift(stats[X->C_n], d_n))

The mu^2 factors in the variance denominators cancel, leaving mu-independent weights.

To ensure all branch lengths are at least min_branch_length (= 0.1 days), clamp:

    tau_X = clamp(tau_X,
                  tau_P + min_branch_length,
                  min(max_node_time[X->C_1] - min_branch_length, max_node_time[X->C_2] - min_branch_length))

The lower bound ensures the P-X branch has length >= min_branch_length.  The upper bound subtracts
min_branch_length from each child's max_node_time because `max_node_time[X->C_n]` is the latest
time **C_n** can have, and X must be at least min_branch_length earlier.  (This is not double-counting:
the `- min_branch_length` applied during Passes 1/2 accounts for branches *below* C_n, while this one
accounts for the X-C_n branch itself.)  The lower bound never exceeds the upper bound:
clamping tau_R at the root (see above) establishes the invariant, and it propagates by
induction since max_node_time[GP->P] <= max_node_time[P->X] - min_branch_length for any child X of P.
Because clamped values feed into downstream nodes, the resulting timed tree is plausible
even when the unconstrained optimum would produce very short branches.

The node times are stored in the result struct so that `utree_to_phylo_tree` can use them
directly instead of computing times from a linear model.

### R^2 computation

For reporting purposes, R^2 is computed at the winning root after node timing, via
a single O(N) pass that computes unweighted root-to-tip mutation counts and evaluates
the standard OLS R^2 formula (same approach as Round 4bis).

## Comparison with OLS and GLS regression rooting

| Property | OLS (Round 4) | GLS (Round 4bis) | Gaussian (this plan) |
|----------|---------------|------------------|----------------------|
| What it optimizes | R^2 of OLS regression | chi^2 of GLS regression | Reduced potential (node times integrated out) |
| Root-rate coupling | Joint (R^2 jointly over root, rate, t_MRCA) | Joint (chi^2 jointly over root, rate, t_MRCA) | Joint (mu and root co-optimized analytically) |
| Variance structure | Homoscedastic (equal-weight tips) | Heteroscedastic (precision-weighted) | Per-branch Gaussian (variance = (d+1)/mu^2) |
| Uses tree topology | Only for root-to-tip distances | Covariance from shared branches | Branch-level model, node times integrated out |
| Zero-mutation branches | Contribute nothing to regression | Regularized (epsilon floor) | Expected length 1/mu (Gamma posterior) |
| Tip-date uncertainty | Uses date midpoints | Uses date midpoints | Uses date midpoints (could be extended) |
| Effective sample size | N (full) | Collapses near zero | TBD (expected to be better than GLS) |
| Per-arc statistics | 5 fields (Ols_stats) | 6 fields (Gls_stats) | 6 fields (Gaussian_stats) |
| Combine operation | Element-wise addition | Element-wise addition | Weighted average + residual costs |
| Node timing output | Linear model (lambda * t + c) | Linear model | Direct optimal times for all nodes |

## Fallback to midpoint rooting

`gaussian_root_utree` falls back to `midpoint_root_utree` when:

- All tips have the same date (no temporal signal).
- N <= 2 (degenerate tree).
- No candidate position across any edge has F < 0 (no position gives a positive
  clock signal).

## Integration

### New function

```cpp
auto gaussian_root_utree(Utree& tree, const std::vector<Tip_desc>& tip_descs,
                         absl::BitGenRef bitgen)
    -> Rooting_info;
```

### Rooting_info augmentation

Add a `node_times` field to `Rooting_info`:

```cpp
struct Rooting_info {
  Node_index root;
  Rooting_method method;
  double r2;
  double lambda;
  double t_MRCA;
  std::vector<double> node_times;  // indexed by Node_index; empty for OLS/midpoint
};
```

When `node_times` is non-empty, `utree_to_phylo_tree` uses the provided times directly
for each inner node.  When empty, it falls back to the existing linear model
(t_node = t_MRCA + mutations_from_root / lambda).

### Rooting_method update

```cpp
enum class Rooting_method { regression, gaussian, midpoint };
```

### Pipeline change

In `build_initial_phylo_tree`, replace:

```cpp
auto rooting_info = ols_regression_root_utree(utree, tip_descs, bitgen);
```

with:

```cpp
auto rooting_info = gaussian_root_utree(utree, tip_descs, bitgen);
```

The existing OLS and GLS functions are retained for comparison.

## Files to modify

| File | Action | What |
|------|--------|------|
| `core/utree.h` | Modify | Add `Rooting_method::gaussian`; add `node_times` to `Rooting_info`; declare `gaussian_root_utree` |
| `core/utree.cpp` | Modify | Implement `gaussian_root_utree`; update `utree_to_phylo_tree` to use `node_times`; update `build_initial_phylo_tree` |
| `core/cmdline.cpp` | Modify | Update rooting method formatting string |
| `tests/utree_tests.cpp` | Modify | Add `gaussian_root_*` tests |

## Verification

1. **Numerical validation on small trees** (3-5 tips with known structure):
   - Construct a tree where the correct root and rate are known analytically.
   - Verify that the A-G coefficients match hand-computed values.
   - Verify that mu, t_MRCA, and Delta_min match direct evaluation.
   - Verify that node times match the weighted-average formula.

2. **Consistency with OLS**: On well-behaved datasets (strong clock signal), verify
   that the winning root position and estimated rate are close to the OLS results.

3. **Degenerate cases**:
   - All-same-date tips: verify fallback to midpoint rooting.
   - 2-tip tree: verify fallback.
   - Zero-mutation branches: verify sensible behavior (expected length 1/mu).
   - All-zero-mutation tree with varied tip dates: verify that rooting still works
     (clock signal comes from branch-count asymmetry, not mutations).
   - All-zero-mutation tree with same-date tips: verify fallback (no signal at all).

4. **Node timing**: Verify that Pass 4 produces non-negative branch lengths and that
   the resulting Phylo_tree passes `assert_phylo_tree_matches_tip_descs`.

5. **End-to-end**: Run `delphy_ui` with Gaussian rooting on standard demo datasets.
   Verify results are sensible and compare with OLS.

6. **Build**: `cmake --build build/debug --target tests -j$(nproc) && build/debug/tests/tests --gtest_filter='Utree*'`

## Post-implementation findings

### Observed behavior

After implementing `gaussian_root_utree` as described above, we tested it on the
h5n1-andersen-2025 dataset (3335 tips), the same dataset used for the Round 4bis
GLS evaluation.

| Method   | lambda (mut/yr) | R^2    | t_MRCA     |
|----------|----------------:|-------:|------------|
| OLS      |           37.84 |  0.92  | 2023-12-16 |
| GLS      |            1.6  | ~0     | 2011-10-30 |
| Gaussian |           12.86 |  0.91  | 2023-05-22 |

The Gaussian model produces a sensible root (t_MRCA in 2023, positive R^2), unlike
GLS which collapsed.  The rate estimate (12.86 mut/yr) is lower than OLS (37.84),
and the root is placed about 7 months earlier.  The remaining discrepancy relative
to OLS is likely due to the lack of an ancestry prior (such as a coalescent): without
it, the Gaussian model assigns longer expected branch lengths to zero-mutation
branches (1/mu each), which pulls the rate estimate down and the root earlier.
This is acceptable for initialization --- the MCMC will refine both.

### Comparison across datasets

All three methods (OLS, GLS, Gaussian) applied to the same refined divergence tree for
each dataset.  MCMC columns are posterior means from Delphy MCMC runs (from
`delphy-2026-paper-data`); MCMC lambda = clock.rate * L (for mpox, site-weighted mean of
APOBEC3 and non-APOBEC3 rates).  Note that repeated runs have some variability, but the
general pattern among the methods is robust.

| Dataset                |  Method  |   R^2  | lambda (mut/yr) |   t_MRCA   |
|------------------------|----------|-------:|----------------:|:----------:|
| Ebola (Dudas)          | MCMC     |        |           21.05 | 2014-01-04 |
|                        | OLS      | 0.6013 |           17.17 | 2013-08-11 |
|                        | Gaussian | 0.5970 |           10.89 | 2013-07-08 |
|                        | GLS      | 0.1434 |            4.57 | 2012-06-30 |
|------------------------|----------|-------:|----------------:|:----------:|
| Ebola (Gire)           | MCMC     |        |           42.05 | 2014-03-02 |
|                        | OLS      | 0.7074 |           45.81 | 2014-02-18 |
|                        | Gaussian | 0.4748 |           12.88 | 2013-10-25 |
|                        | GLS      | 0.3149 |            0.86 | 2010-07-22 |
|------------------------|----------|-------:|----------------:|:----------:|
| H3N2                   | MCMC     |        |            5.89 | 1997-10-19 |
|                        | OLS      | 0.9479 |            9.39 | 1999-04-28 |
|                        | Gaussian | 0.9028 |            5.11 | 1996-04-11 |
|                        | GLS      | 0.8053 |            3.69 | 1995-07-18 |
|------------------------|----------|-------:|----------------:|:----------:|
| H5N1                   | MCMC     |        |           42.40 | 2023-12-24 |
|                        | OLS      | 0.9128 |           37.89 | 2023-12-16 |
|                        | Gaussian | 0.9061 |           13.20 | 2023-03-16 |
|                        | GLS      | 0.0004 |            1.65 | 2011-10-27 |
|------------------------|----------|-------:|----------------:|:----------:|
| mpox (O'Toole)         | MCMC     |        |           10.18 | 2016-05-02 |
|                        | OLS      | 0.8372 |            7.87 | 2016-12-26 |
|                        | Gaussian | 0.8206 |            7.72 | 2016-01-05 |
|                        | GLS      | 0.8098 |            7.11 | 2016-04-30 |
|------------------------|----------|-------:|----------------:|:----------:|
| mpox (Parker)          | MCMC     |        |            9.96 | 2016-03-29 |
|                        | OLS      | 0.9348 |            7.31 | 2016-07-24 |
|                        | Gaussian | 0.9345 |            7.37 | 2015-09-15 |
|                        | GLS      | 0.9277 |            5.77 | 2015-12-15 |
|------------------------|----------|-------:|----------------:|:----------:|
| SARS-CoV-2 (Lemieux)   | MCMC     |        |           30.02 | 2019-12-04 |
|                        | OLS      | 0.1987 |           19.82 | 2019-08-31 |
|                        | Gaussian | 0.1916 |           10.84 | 2019-03-20 |
|                        | GLS      | 0.1766 |            1.82 | 2016-04-29 |
|------------------------|----------|-------:|----------------:|:----------:|
| Zika                   | MCMC     |        |           10.49 | 2013-07-12 |
|                        | OLS      | 0.4449 |           17.50 | 2014-02-05 |
|                        | Gaussian | 0.3644 |            4.76 | 2011-05-16 |
|                        | GLS      | 0.2835 |            3.17 | 2010-09-22 |
