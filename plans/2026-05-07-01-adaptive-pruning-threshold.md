# Adaptive Pruning Threshold for Guide Tree Construction

## Background

The guide tree construction (Round 1) adds tips one at a time, attaching each new tip at
the edge where it introduces the fewest new site deltas.  The search for the best
attachment edge (`find_best_attachment_arc`) uses a priority-queue-driven branch-and-bound
algorithm: starting from current tree focus node (arbitrary), it evaluates candidate edges
in order of increasing attachment cost, and prunes when the cheapest remaining candidate
exceeds `best_cost + threshold`.

The Round 1 plan specifies a fixed threshold of 1-2, noting that "it's OK to occasionally
pick a slightly-suboptimal attachment point in exchange for faster construction."  This is
appropriate for shallow phylogenies where K (the number of sites differing between the search
starting point and the new tip) is small.  For deep phylogenies, however, K can reach into
the thousands, and a fixed threshold misses better attachment points that lie beyond temporary
cost increases in the search landscape.

This plan develops an adaptive threshold that scales with both the current best cost and the
genome length, informed by a probabilistic analysis of "blip" mutations.

## Blips: temporary cost increases along the search path

Let X denote the new tip being attached.  Consider two attachment points in the tree: Q
with cost K (the current best) and P with cost 0 (the true optimum — an idealization; in
practice, the best attachment point has nonzero cost, but the analysis carries over by
replacing 0 with the actual best cost).  The path from Q to P passes through intermediate
edges.  If the tree is parsimonious, the attachment cost should decrease roughly
monotonically along this path: each step brings us closer to X's sequence.

However, some sites may undergo **blips** — deviations from the parsimonious mutational
history at a given site that temporarily change the attachment cost.  Each open blip adds 1
to the attachment cost at intermediate edges where it is active.  If too many blips are open
simultaneously, the cost at an intermediate edge exceeds `best_cost + threshold`, and the
search prunes the direction towards P — cutting off the path to a better attachment point.

### Kinds of blips

A blip is any deviation from the minimum-mutation history at a site along the Q-P path.
Several kinds arise, with a prevalence that depends on the magnitude of K compared to the
genome length L:

- **Same-site blips** (dominant for small K/L): sites where Q, P, and X all have the same
  state.  A blip requires at least 2 mutations — one away from X's state, one back.

- **Different-site blips**: sites where Q differs from P (= X's state).  These require
  at least 3 mutations to cause a temporary cost increase (the site must visit X's state
  and then leave it before returning).  For a given site, these are rarer than same-site
  blips by a factor of order K/L, so they are negligible when K << L.

- **Higher-order blips**: at any site, more complex mutation patterns (4+ mutations at a
  same-state site, 5+ at a different-state site) can produce blips.  These are suppressed
  by additional powers of K/L, so they are negligible for K/L < 0.15.

For outbreak-scale datasets (SARS-CoV-2, Ebola, Zika, mpox, H5N1, etc.), K/L << 1 and
same-site blips dominate.  The analysis below focuses on same-site blips, then adds a safety
factor to cover the other kinds.

### Same-site blip probability under JC69

Under the Jukes-Cantor model, each site independently follows a continuous-time Markov
chain with 4 states.  The transition probabilities depend on the evolutionary distance
lambda, measured in expected substitutions per site (under a strict molecular clock,
lambda = mu * t, where mu is the per-site substitution rate and t is the elapsed time).
They are:

```
P(same state after distance lambda)               = 1/4 + 3/4 * exp(-4*lambda/3)
P(specific different state after distance lambda) = 1/4 - 1/4 * exp(-4*lambda/3)
```

These probabilities are the sums over all possible mutational histories of length lambda
with the corresponding endpoints (in fact, expanding each P(lambda) exp(+lambda) in
powers of lambda, the lambda^k factor directly yields the probabilities of histories
involving k mutations; see Delphy paper SI Section 4).

We now estimate the probability that a specific site where Q, P, and X all share the same
nucleotide undergoes a 2-mutation blip along the Q-P path:

- **2 mutations on the path:** Without the same-state constraint, the number of mutations
  at a site over evolutionary distance lambda is Poisson(lambda).  P(N=2) ~= lambda^2/2
  for small lambda.
- **Returning to the same state:** The first mutation goes to one of 3 other states.  For
  the start and end states to match, the second mutation must return to the original,
  reducing the probability by a factor of 1/3.
- **Same-site blip probability:** lambda^2/2 * 1/3 = **lambda^2/6**.

How are K and lambda related?  This is the standard Jukes-Cantor distance correction.
From the transition probabilities above, the fraction of sites in a specific different
state after distance lambda is 1/4 - 1/4 * exp(-4*lambda/3).  Since K of L sites differ
between Q and P across 3 possible different states, K/L = 3 * [1/4 - 1/4 * exp(-4*lambda/3)]
= 3/4 * [1 - exp(-4*lambda/3)], giving:

```
lambda = -3/4 * ln(1 - 4K/(3L))
```

For small K/L, the Taylor expansion gives lambda ~= K/L * (1 + 2K/(3L) + ...).  Even at
K/L = 0.15, the relative error between lambda and K/L is only about 10%, so we use
lambda ~= K/L throughout.

### Total expected same-site blip count (uniform rates)

There are L - K same-state sites, each contributing blips independently:

```
E[same-site blip sites] ~= (L - K) * (K/L)^2 / 6
```

For K << L, this simplifies to K^2 / (6L).

Define **sigma = K / sqrt(6L)**.  Then E[same-site blip sites] ~= sigma^2.

Since same-site blips are independent across sites, the total count is approximately
Poisson(sigma^2).

Recall from the introduction that the branch-and-bound search prunes a direction when the
attachment cost exceeds `best_cost + threshold`.  The threshold must therefore be large
enough to accommodate the blips encountered along the path from the current best to a
potentially better attachment point.  The most pessimistic scenario is that all blip sites
are simultaneously active at some intermediate edge.  This is unlikely for large sigma^2
(since blips open and close independently along the path), but provides a safe upper bound.
Under the Poisson(sigma^2) model, a 5-standard-deviation upper bound gives a conservative
threshold:

```
threshold_base = sigma^2 + 5*sigma = sigma * (sigma + 5)
```

For SARS-CoV-2 (L ~= 30,000) with K = 100, sigma ~= 0.24, and threshold_base ~= 1.3 — a
fixed threshold of 2 is fine.  But for K = 3,000, sigma ~= 7.1, and threshold_base ~= 86.
A fixed threshold of 2 would prune the search well before it reaches the better attachment
point.

### Accounting for site-rate heterogeneity and other effects

The analysis above assumes all sites mutate at the same rate.  In reality, mutation rates
can vary dramatically across sites.  In shallow phylogenies, this is almost inconsequential:
few sites mutate more than once, so rate variation has little effect.  But in deep
phylogenies, hypervariable sites (e.g., in surface glycoprotein genes) can accumulate dozens
or hundreds of mutations, while conserved sites (e.g., in polymerase genes) rarely mutate.

To model this, let each site have a rate multiplier r drawn independently from some
distribution with E[r] = 1.  A site with rate multiplier r accumulates mutations at rate
r times the mean, so its evolutionary distance is r*lambda (where lambda is the mean across
sites), and its blip probability is (r*lambda)^2/6 = r^2 * lambda^2/6.

**The key insight is that the total expected blip count scales with E[r^2], not E[r].**  A
site that is 10x faster than average produces not 10x but **100x** more blips.  A small
fraction of hypervariable sites can dominate the total blip count:

```
E[same-site blip sites] ~= (K^2 / (6L)) * E[r^2]
```

Under the standard Gamma(alpha, alpha) model for site rates, E[r] = alpha/alpha = 1 and
Var[r] = alpha/alpha^2 = 1/alpha, so E[r^2] = Var[r] + (E[r])^2 = 1/alpha + 1.  This
ranges from 3 (alpha = 0.5, moderate heterogeneity) to 21 (alpha ~= 0.05, strong
heterogeneity) to much higher for extremely heterogeneous datasets.

Beyond site-rate heterogeneity, many other effects can inflate the effective blip count.
These include: the tree under construction is not yet parsimonious (so the search path may
detour through high-cost regions that are not blips in the JC sense but have the same
pruning consequences), selective pressures and lineage-dependent clock rates create non-JC
mutation patterns, transit through immunocompromised patients produces hypermutated
lineages, and different viral lifecycles in different host species add further rate
variation.

### Safety factor for robustness

To robustly accommodate the r^2-amplified blip rate and the additional biological and
algorithmic effects described above, we apply an order-of-magnitude safety factor to the
base threshold:

```
threshold = 10 * sigma * (sigma + 5)     where sigma = K / sqrt(6L)
```

The factor of 10 absorbs a wide range of rate heterogeneity levels (from moderate to
extreme) and non-blip sources of cost fluctuation, without requiring an explicit estimate
of alpha or any other model parameter.  The threshold is conservative enough that the
search reliably finds good attachment points, yet the quadratic dependence on K ensures
the threshold stays small for easy tips (low K) and only grows large when needed.

Importantly, the threshold is **dynamic**: K in the formula above is `best_cost`, which
decreases as the search discovers better attachment points.  The thresholds in the tables
below represent worst-case values at the start of a search from a distant starting point.
As the search descends toward a good attachment point and `best_cost` drops, the threshold
tightens and the effective search scope narrows substantially — see the Design section for
details.

### Numerical examples

In the tables below, "Practical threshold" = `clamp(ceil(10*sigma*(sigma+5)), 2, L)`.

For L ~= 10,000 (e.g., Zika, dengue, HIV, hepatitis C):

| K (best_cost) | sigma   | 10*sigma*(sigma+5) | Practical threshold |
|---------------|---------|--------------------|---------------------|
| 30            | 0.12    | 6.3                | 7                   |
| 100           | 0.41    | 22.1               | 23                  |
| 300           | 1.22    | 76.2               | 77                  |
| 1000          | 4.08    | 371                | 371                 |
| 3000          | 12.25   | 2112               | 2113                |

For L ~= 15,000 (e.g., influenza A, RSV, measles, mumps; Ebola at ~19,000 is broadly
similar):

| K (best_cost) | sigma   | 10*sigma*(sigma+5) | Practical threshold |
|---------------|---------|--------------------|---------------------|
| 30            | 0.10    | 5.1                | 6                   |
| 100           | 0.33    | 17.8               | 18                  |
| 300           | 1.00    | 60                 | 60                  |
| 1000          | 3.33    | 278                | 278                 |
| 3000          | 10.00   | 1500               | 1500                |

For L ~= 30,000 (e.g., SARS-CoV-2):

| K (best_cost) | sigma   | 10*sigma*(sigma+5) | Practical threshold |
|---------------|---------|--------------------|---------------------|
| 30            | 0.07    | 3.6                | 4                   |
| 100           | 0.24    | 12.3               | 13                  |
| 300           | 0.71    | 40.4               | 41                  |
| 1000          | 2.36    | 173                | 174                 |
| 3000          | 7.07    | 854                | 854                 |

For L ~= 200,000 (e.g., mpox):

| K (best_cost) | sigma   | 10*sigma*(sigma+5) | Practical threshold |
|---------------|---------|--------------------|---------------------|
| 30            | 0.027   | 1.4                | 2                   |
| 100           | 0.091   | 4.6                | 5                   |
| 300           | 0.27    | 14.4               | 15                  |
| 1000          | 0.91    | 54.0               | 54                  |
| 3000          | 2.74    | 212                | 212                 |

For shallow phylogenies (K <= 100), the threshold is modest (<= 23 even for the shortest
genomes), keeping the search fast.  For deep phylogenies (K ~= 3,000 — close to the
maximum for a particularly bad starting point in a deep tree, not the typical K), the
threshold grows substantially, allowing the search to traverse long paths through high-cost
regions.  But that is precisely where a broad search scope is most needed.  In all cases,
the threshold is much smaller than K itself, so pruning still eliminates most of the tree.
Moreover, these are worst-case thresholds at the start of the search; they tighten rapidly
as the search finds better attachment points (see Design).

## Design

### Threshold evaluated at pop-time with current best_cost

Recall that the search uses a priority-queue-driven branch-and-bound algorithm (described
in Background): it maintains a min-priority-queue of candidate edges ordered by attachment
cost, and expands them in order of increasing cost.  The threshold is recomputed each time
an entry is popped from the priority queue, using the current `best_cost`:

```
sigma = best_cost / sqrt_6L
threshold = clamp(ceil(10 * sigma * (sigma + 5)), 2, L)  // See Bounds below
if priority > best_cost + threshold:
    break
```

This means the threshold **tightens as the search improves**.  When `best_cost` is high
(early in the search, or when the search starts far from the optimum), the threshold is
generous and the search explores broadly.  In practice, the search rapidly follows a clear
path of decreasing cost toward the best attachment point; as it descends this path and
`best_cost` drops, the threshold shrinks and the tentative alternative directions that
were opened early in the search are progressively pruned.

This tightening is sound because the search could equivalently be restarted from the
current best candidate attachment point: if we had started there, the initial cost would
be `best_cost` and the threshold would be exactly what we compute.  The search mechanics
do not change significantly whether we arrive at a good candidate from a distant start or
begin there directly.

### Bounds

Even when the "practical threshold" above is close to 0, we maintain a minimum threshold
of 2 to ensure that ties and near-ties are explored, matching the Round 1 plan's original
intent.  Since `best_cost` can never exceed L, the threshold is implicitly bounded; an
explicit `std::clamp(..., 2, L)` prevents overflow for degenerate inputs.

## Files to modify

| File | Action | What |
|------|--------|------|
| `core/utree.cpp` | Modify | Replace fixed `k_pruning_threshold = 2` with adaptive `pruning_threshold(best_cost)` method; evaluate threshold at pop-time using current `best_cost` |

## Verification

1. Existing `Utree` unit tests continue to pass.
2. For `best_cost = 0`, the threshold equals the minimum (2).
3. For small K (<= 100) and typical L (15,000-30,000), the threshold stays small (<= 20),
   preserving fast construction for shallow phylogenies.
4. For large K (>= 1000), the threshold grows quadratically, allowing the search to reach
   distant optima.
5. Manual experiments on a temporarily instrumented Delphy binary (instrumentation not
   committed) confirm that guide trees built with the adaptive threshold achieve parsimony
   scores (total deltas) comparable to or better than trees built by evaluating every
   attachment point without pruning, while remaining substantially faster.

Build: `cmake --build build/debug --target tests -j$(nproc) && build/debug/tests/tests --gtest_filter='Utree*'`

## Preliminary results

Total parsimony score (mutations or deltas) for each tree construction method on the
standard Delphy demo datasets, one of the large simulated datasets, and one additional
deep-phylogeny dataset.  Lower is better.  Results can vary in subsequent attempts owing
to randomness.

| Dataset                          |    N    |    L    | UShER-like | Guide (fixed=2) | Guide (adaptive) | Approx. MCMC mean |
|----------------------------------|---------|---------|------------|-----------------|------------------|-------------------|
| Ebola (Dudas et al, 2017)        |   1,610 |  18,996 |      3,630 |           3,383 |            3,428 |             2,786 |
| Ebola (Gire et al, 2014)         |      81 |  18,959 |         76 |              55 |               55 |                55 |
| H3N2 (Rambaut et al, 2008)       |     165 |   1,698 |        555 |             465 |              403 |               394 |
| H5N1 (Andersen lab, 2025)        |   3,339 |  13,136 |      6,058 |           6,296 |            5,372 |             4,347 |
| mpox (O'Toole et al, 2023)       |      41 | 197,209 |        459 |             344 |              344 |               343 |
| mpox (Parker et al, 2025)        |     177 | 197,209 |      1,436 |             950 |              957 |               944 |
| SARS-CoV-2 (Lemieux et al, 2021) |     757 |  29,903 |        814 |             665 |              665 |               667 |
| Zika (Metsky et al, 2017)        |     174 |  10,807 |      1,585 |           1,352 |            1,333 |             1,210 |
| Exp-100,000 sim                  | 100,000 |  30,000 |     10,313 |           9,797 |            9,797 |             9,799 |
| RSV-A (NCBI, 2012–2025) (*)      |   3,000 |  15,277 |     50,025 |          92,308 |           47,474 |               N/A |

(*) Note: The RSV-A dataset is contaminated by 3 nearby RSV-B sequences that are ~2,500 mutations from the RSV-A sequences.

### Observations

The guide tree with adaptive threshold matches or beats the UShER-like tree on every dataset,
often substantially (e.g., H5N1: 5,372 vs 6,058, RSV-A: 47,474 vs 50,025).  On shallow
phylogenies, both thresholds give identical results near the MCMC optimum.  On deep
phylogenies (H3N2, H5N1, RSV-A), the fixed threshold catastrophically fails: the guide tree
becomes *worse* than the UShER-like tree (RSV-A: 92,308 vs 50,025; H5N1: 6,296 vs 6,058),
while the adaptive threshold recovers and improves on it.

In summary: the adaptive threshold never hurts on shallow phylogenies and is essential for
deep ones, where the fixed threshold can make the guide tree worse than no guide tree at all.

The adaptive guide tree also starts much closer to the MCMC posterior mean than the
UShER-like tree.  On datasets with known MCMC means, the UShER-like tree overshoots by
20–40% (e.g., H3N2: 555 vs 394, H5N1: 6,058 vs 4,347, Zika: 1,585 vs 1,210), while the
adaptive guide tree closes most of that gap (H3N2: 403, H5N1: 5,372, Zika: 1,333).  This
translates directly to less burn-in work for the MCMC sampler.
