# Better Tree Init, Round 4bis: GLS Regression Rooting

This plan extends
[Round 4](2026-05-18-01-better-tree-init-round4-regression-rooting.md), which implemented
root-to-tip regression rooting using ordinary least squares (OLS).

Conversations with the authors of [Clockor2](https://github.com/clockor2/clockor2)
(Featherstone et al., *Systematic Biology* 2024,
[doi:10.1093/sysbio/syae003](https://doi.org/10.1093/sysbio/syae003)) pointed out that the
O(N) OLS algorithm in Round 4 is not novel: it is conceptually identical to the
linear-time regression rooting described in the original
[TreeTime](https://github.com/neherlab/treetime) paper (Sagulenko, Puller & Neher, *Virus
Evolution* 2018, [doi:10.1093/ve/vex042](https://doi.org/10.1093/ve/vex042)).  The O(N)
OLS rooting first appeared in the TreeTime repo on 2016-02-17 (commit
[`d7aaa326`](https://github.com/neherlab/treetime/commit/d7aaa326d60180202d183b6c2bacc283420fc274),
"[Feature] best root and regression search in linear time").  I had completely missed
this.

Six months after the original TreeTime paper was published, TreeTime's code evolved.  In
July 2018, a new
[`TreeRegression`](https://github.com/neherlab/treetime/blob/17da461f9299706b99a90e622183772f00c5076d/treetime/treeregression.py)
class was added, replacing the original OLS with a more sophisticated **generalized least
squares** (GLS) regression that accounts for the phylogenetic correlation structure of
root-to-tip distances.  The GLS approach is described in a September 2018 bioRxiv preprint
by Neher, ["Efficient estimation of evolutionary rates by covariance aware
regression"](https://doi.org/10.1101/408005), which credits Ho & Ané for the efficient
covariance inversion technique.  This preprint was never published (Neher's [publications
page](https://neherlab.org/pages/publications.html) lists it as "won't submit").  The
original OLS code (in `treetime.py`'s `find_best_root_and_regression`, using the same 5
per-subtree statistics as our Round 4) has since been removed.  The OLS behavior is now a
special case of `TreeRegression`: passing `covariation=False` to
[`setup_TreeRegression`](https://github.com/neherlab/treetime/blob/17da461f9299706b99a90e622183772f00c5076d/treetime/clock_tree.py#L260)
sets all internal branch variances to zero and tip branch variances to 1 (a uniform
constant rather than the mutation count), almost reducing GLS to OLS (technically,
TreeTime always minimizes chi^2 even when `covariation=True`, whereas OLS maximizes R^2,
which is leads to subtly different comparisons between candidate root positions on
different branches; see details below).

This plan explains how the GLS works, derives it from first principles, and then goes on to
implement it for Delphy's tree initialization.  We will not replace the existing OLS code;
the new GLS function (`gls_regression_root_utree`) will have the same signature as the
existing OLS function (renamed from `regression_root_utree` to `ols_regression_root_utree`),
so switching between them is a one-line change in `build_initial_phylo_tree`.

> **Status (2026-05-27):** The GLS function was implemented and is correct (verified against
> direct matrix inversion and against TreeTime's formulas), but experiments on real data
> showed that GLS rooting collapses the effective sample size and produces unusable results.
> `build_initial_phylo_tree` continues to use OLS.  The GLS code is retained but not called.
> See [Post-implementation findings](#post-implementation-findings-gls-gives-poor-rooting-results)
> at the end of this document for details.

## Notation

Throughout this plan, we distinguish the following quantities associated with tree branches:

- **Time lengths** (tau_c, tau_1', ...): the durations in days of tree branches.  These are
  unknown parameters that we are trying to estimate.
- **Per-branch mutation counts** (Zc / zc, Z1' / z1', ...): Z_b is the Poisson random
  variable modeling the number of mutations generated on branch b under a given rate and
  branch time; z_b (lowercase) is its observed realization.  We assume branches are short
  enough that mutation counts and delta counts are interchangeable (no reversals or
  intermediate states).  The observed counts z_b are treated as **data**, fixed and known.
- **Root-to-tip mutation counts** (M1 / m1, M2 / m2, ...): M_i is the sum of per-branch
  random variables along the path from root to tip i; m_i (lowercase) is its observed
  realization.  These are the quantities we regress against tip dates.

All tip dates are used as centered deviations: dt_i = t_i - mean_t, where
mean_t = (1/N) * sum(t_i).  This prevents catastrophic cancellation in expressions involving
sums and products of dates.

We use lambda for the overall mutation rate (mutations per day, summed across all sites),
following Delphy's convention.  (In the broader phylogenetics literature, mu typically denotes
the per-site rate; lambda = mu * L where L is the number of sites.)

We write Sigma for the covariance matrix of the root-to-tip mutation counts and
W = Sigma^{-1} for the precision matrix, which plays the role of a weight matrix in the
GLS regression.

## Motivation: why OLS is insufficient

OLS regression of root-to-tip mutation counts against tip dates treats all such counts as
independent observations with equal variance.  Neither assumption holds.

Consider a tree with two tips sharing a common ancestor:

```
                  root
                   |
                   |
         branch c  |
       time tau_c  |
          zc muts  |
                   |
                   |
                 split
                /     \
               /       \
   branch 1'  /         \  branch 2'
time tau_1'  /           \  time tau_2'
  z1' muts  /             \  z2' muts
           /               \
          /                 \
       tip 1               tip 2
```

The observed root-to-tip mutation counts are m1 = zc + z1' and m2 = zc + z2'.  Suppose we
knew the branch times tau_c, tau_1', and tau_2' and simulated neutral evolution forward along
the branches of the tree.  Assuming that mutations accumulate via a Poisson process with
overall rate lambda, the number of mutations on each branch would be an independent random
variable:

    Zc  ~ Poisson(lambda * tau_c)
    Z1' ~ Poisson(lambda * tau_1')
    Z2' ~ Poisson(lambda * tau_2')

From these, we can derive the random root-to-tip mutation counts M1 = Zc + Z1' and
M2 = Zc + Z2'.  These are correlated because they share the Zc component (mutations on the
common stem).  It is the joint distribution of (M1, M2) that we will compare against the
observed counts (m1, m2) to estimate the rate and root date.

Since the per-branch mutation counts are independent, the means, variances, and covariances
of the root-to-tip mutation counts follow directly from the Poisson model:

    E[M1] = lambda * tau_1          where tau_1 = tau_c + tau_1'
    E[M2] = lambda * tau_2          where tau_2 = tau_c + tau_2'

    Var(M1) = lambda * tau_1        [= Var(Zc) + Var(Z1') = lambda*tau_c + lambda*tau_1']
    Var(M2) = lambda * tau_2
    Cov(M1, M2) = Var(Zc) = lambda * tau_c

In a general tree with N tips, the covariance matrix of the root-to-tip mutation counts is:

    Sigma_ij = Cov(M_i, M_j) = lambda * (time from root to MRCA of tips i and j)

This matrix has the same hierarchical structure as the tree: tips sharing more ancestry have
higher covariance.

OLS ignores this structure entirely, treating Sigma as proportional to the identity.  GLS
uses Sigma to properly account for two effects that OLS ignores: (1) tips that share more ancestry
have correlated root-to-tip distances, and (2) tips farther from the root have
proportionally higher variance.  By down-weighting correlated and high-variance
observations, GLS produces estimates of the rate and root date that have lower mean
squared error than OLS estimates -- i.e., they are closer to the true values on average,
especially when the tree is unbalanced or has long internal branches.  Figure 4 of Neher
(2018) demonstrates this on simulated data: GLS estimates of both t_MRCA and the mutation
rate have dramatically reduced variance compared to OLS.

## Breaking the circularity: mutation counts as variances

The covariance matrix Sigma depends on lambda and the branch times -- the very quantities we
are trying to estimate.  This apparent circularity breaks because the branch times tau and
the rate lambda always appear as their product lambda * tau, which is the expected number of
mutations on a branch.  Under the Poisson model, this expected count is also the variance:

    E[Zc] = Var(Zc) = lambda * tau_c

We estimate this product directly from the observed mutation counts:

    Var(Zc) = lambda * tau_c ~= zc     (without regularization; see next section)
    Var(Z1') = lambda * tau_1' ~= z1'
    Var(Z2') = lambda * tau_2' ~= z2'

So the covariance matrix can be expressed entirely in terms of observed mutation counts:

    Var(M1) ~= z1' + zc = m1        (total mutations, root to tip 1)
    Var(M2) ~= z2' + zc = m2        (total mutations, root to tip 2)
    Cov(M1, M2) ~= zc               (mutations on common stem)

More generally, for any pair of tips:

    Sigma_ij ~= (total mutations on the path from root to MRCA(i, j))  (without regularization)

Every entry is a sum of observed mutation counts on branches.  No unknown parameters appear.
The mutation counts were already estimated during phylogenetic reconstruction (in our case,
from the refined `Utree` as delta counts), and we reuse them here as variance/covariance estimates.

(This is an approximation: the mutation counts are themselves estimated with uncertainty, and
using the same data for both tree construction and regression is a form of double-dipping.
But it is standard practice -- TempEst, Clockor2, and TreeTime all do this -- and the result
is refined by Delphy's MCMC anyway.  Remember: we're using this GLS procedure to arrive at a better
starting point for the MCMC, not as the final answer; removing all the above approximations and
apparent circularity is conceptually equivalent to performing the actual MCMC.)

### Regularization

The previous section modeled the variance of each branch's mutation count as
Var(Z_branch) ~= z_branch (from the Poisson property).  We refine this by adding a small
regularization term, following TreeTime:

    Var(Z_branch) ~= z_branch + epsilon   (with regularization)

where:

    epsilon = 0.05 * (sum of all branch mutation counts in the tree) / N

The sum of all branch mutation counts is obtained via `tree.count_deltas()`, which iterates
over every other arc index (i += 2) to count each undirected edge exactly once.  epsilon is
proportional to the mean branch mutation count, ensuring that the regularization scales with
the overall evolutionary signal in the data.  (For a bifurcating tree with N tips, there are
~2N branches, so epsilon ~= 0.1 times the mean branch mutation count.)

Without regularization, a branch with zero observed mutations would have Var(Z_branch) = 0,
contributing nothing to the variance of tips in its subtree.  Tips reached through many
zero-mutation branches would have artificially low total variance Sigma_ii, causing the GLS to
give them disproportionate weight in the rate and root-date estimates.  If all branches on
a path have zero mutations, Sigma becomes singular.

With regularization, the covariance matrix entries become:

    Sigma_ij = sum over shared branches of (z_branch + epsilon)
         = (total mutations on the path from root to MRCA(i,j))
           + (number of branches on that shared path) * epsilon    (with regularization)

The first term is the data-driven covariance from above.  The second is a small
tree-topology-dependent correction that prevents any branch from having zero variance.
For a branch with many mutations, epsilon is a small perturbation.  For a zero-mutation
branch, epsilon provides a floor that keeps the precision finite.

## The GLS regression

Under the Poisson model and the molecular clock, the expected root-to-tip mutation count
for tip i is:

    E[M_i] = lambda * tau_i
           = lambda * (t_i - t_MRCA)
           = lambda * dt_i + lambda * (mean_t - t_MRCA)

We now reparametrize this model with two unknowns:

- **alpha** = lambda: the evolutionary rate (mutations per day, summed across all sites).
- **beta** = lambda * (mean_t - t_MRCA): the expected number of mutations from the root
  to a tip sampled at the mean date (i.e., at dt_i = 0).

With these parameters,

    E[M_i] = alpha * dt_i + beta

The observed root-to-tip counts m_i are noisy realizations of M_i; two tips with the same
date have the same E[M_i] under the clock but can have very different m_i values.
We now approximate the Poisson noise as Gaussian, which makes the likelihood tractable in
closed form (the exact Poisson joint distribution of correlated root-to-tip counts has no
clean parametric form, and its log-likelihood cannot be maximized in closed form or
decomposed into fixed-size sufficient statistics).  Under this approximation, the vector
M = (M1, ..., MN) follows a multivariate Gaussian with mean alpha * dt + beta * 1 and
covariance matrix Sigma (constructed from branch mutation counts in the previous section).
Writing m = (m1, ..., mN) for the vector of observed root-to-tip mutation counts, the
joint log-likelihood is:

    ln p(m | alpha, beta) = -1/2 * chi^2 - 1/2 * ln det(Sigma) + const

where:

    chi^2 = (m - alpha * dt - beta * 1)^T  W  (m - alpha * dt - beta * 1)

and W = Sigma^{-1} is the precision matrix.  (Since alpha * dt_i + beta = lambda * tau_i, the
residuals m_i - lambda * tau_i do not depend on whether we use alpha & beta or lambda
& t_MRCA as the fit parameters.  Centering the dates is a numerical convenience that
prevents catastrophic cancellation; it does not change chi^2 or the comparison across
candidate roots.)

Since Sigma is determined entirely by the observed mutation counts (not by alpha or beta),
the ln det(Sigma) term is constant with respect to alpha and beta.  Maximizing the likelihood
therefore reduces to minimizing chi^2.

(Strictly, comparing different candidate root positions via profile likelihood requires the
term -1/2 * ln det(Sigma) in addition to -1/2 * chi^2_min.  Both terms vary with the root
position in GLS.  Following TreeTime, we neglect the ln det(Sigma) term: chi^2_min dominates
for any reasonably sized dataset.  If needed, ln det(Sigma) could be tracked incrementally in
O(N) time using the matrix determinant lemma, the determinant analogue of
Sherman-Morrison: det(A + u*v^T) = (1 + v^T * A^{-1} * u) * det(A).  A Bayesian
treatment that integrates out alpha and beta would add a further small correction
(-1/2 * ln det(H), where H is the Hessian of chi^2), but we do not pursue this here.)

Setting the gradient of chi^2 to zero gives the GLS estimates of alpha and beta.  This
requires computing inner products of the form a^T W b, which naively requires forming and
inverting the N x N matrix Sigma.  The next section shows how to avoid this entirely.

## Sufficient statistics for GLS: the `Gls_stats` struct

Expanding chi^2 and collecting terms, we obtain a quadratic function of alpha and
beta whose coefficients involve exactly six inner products with W:

    chi^2 =            1 * sum_m_W_m
          +      alpha^2 * sum_dt_W_dt
          +       beta^2 * sum_1_W_1
          -      2*alpha * sum_m_W_dt
          -       2*beta * sum_m_W_1
          + 2*alpha*beta * sum_dt_W_1

where:

    sum_1_W_1    = 1^T  W 1       (total weight; OLS limit: n)
    sum_dt_W_1   = dt^T W 1       (weighted sum of centered dates; OLS limit: sum_dt = 0)
    sum_m_W_1    = m^T  W 1       (weighted sum of mutations; OLS limit: sum_m)
    sum_dt_W_dt  = dt^T W dt      (weighted sum of squared dates; OLS limit: N * Var_t, const)
    sum_m_W_dt   = m^T  W dt      (weighted cross-product; OLS limit: sum_m_dt)
    sum_m_W_m    = m^T  W m       (weighted sum of squared mutations; OLS limit: sum_m2)

(In the OLS limit, Sigma = W = I.)  These six quantities form a `Gls_stats` struct, the
sufficient statistics for some subset of tips with distances measured from some reference
point.  The naming convention `sum_a_W_b` directly encodes the inner product `a^T W b`,
making the Sherman-Morrison correction pattern (see below) immediately readable.  In the
OLS limit, these reduce to the Round 4 statistics (with the addition of sum_dt_W_dt, which
is a global constant in OLS and was not tracked in Round 4).  (For a subset of k < N tips,
W means the inverse of the k x k submatrix of Sigma for those tips; the formulas for the
elements of this k x k submatrix carry over unchanged from the sections "Breaking the
circularity" and "Regularization".  The Round 4 OLS implementation called these
`Arc_stats`.)

### GLS solution from `Gls_stats`

Setting d(chi^2)/d(beta) = 0:

    beta = (sum_m_W_1 - alpha * sum_dt_W_1) / sum_1_W_1

Setting d(chi^2)/d(alpha) = 0 and substituting beta:

    alpha = (sum_m_W_dt * sum_1_W_1 - sum_m_W_1 * sum_dt_W_1) / (sum_dt_W_dt * sum_1_W_1 - sum_dt_W_1^2)

These two formulas give the full GLS solution in O(1) from `Gls_stats`.  So: if we can
compute `Gls_stats` for a candidate root, we can solve the GLS regression without ever
forming Sigma or W.

### Residual chi^2 and root comparison

Substituting the optimal alpha and beta back into the chi^2 expression and simplifying using
the normal equations (alpha * sum_dt_W_dt + beta * sum_dt_W_1 = sum_m_W_dt and
alpha * sum_dt_W_1 + beta * sum_1_W_1 = sum_m_W_1) gives:

    chi^2_min = sum_m_W_m - alpha * sum_m_W_dt - beta * sum_m_W_1

This is the weighted sum of squared residuals at the best-fit rate and intercept.  To compare
candidate root positions, we pick the one with the smallest chi^2_min (subject to the
validity check alpha > 0; see below).  Ties are broken uniformly at random, following the
same pattern as Round 4: maintain a vector of equally-best (arc, k) pairs and pick one via
`bitgen`.

A candidate position is **valid** only if alpha > 0 (positive molecular clock signal).
Positions with alpha <= 0 are skipped: they imply a non-physical negative or zero mutation
rate.  This is the GLS analog of the Round 4 OLS check Cov(m, t) > 0.

Once the best root position is selected, the rate and root date for `Rooting_info` are
recovered from alpha and beta:

    lambda = alpha
    t_MRCA = mean_t - beta / alpha

For reporting purposes, R^2 is computed at the winning root position only, using the same
unweighted definition as in Round 4.  Since `Gls_stats` tracks precision-weighted inner
products (a^T W b), it does not directly yield the unweighted sums needed for R^2.  Instead,
after root selection, a single O(N) pass over the tree computes the unweighted root-to-tip
mutation counts {m_i} for the chosen root and evaluates:

    mean_m = (1/N) * sum(m_i)
    Cov_mt = (1/N) * sum(m_i * dt_i)      [unweighted]
    Var_m  = (1/N) * sum(m_i^2) - mean_m^2
    R^2    = Cov_mt^2 / (Var_m * Var_t)

R^2 is not the quantity being minimized (that is chi^2_min), but it provides a quick summary
of whether the clock fit is "good" or "bad".  Note that minimizing chi^2_min and maximizing
R^2 are not equivalent, even in the OLS case: chi^2_min = N * Var(m) * (1 - R^2), and Var(m)
varies with root position.  Minimizing chi^2_min is the statistically correct criterion
(it corresponds to maximum likelihood); Round 4's use of R^2 is a simpler heuristic.

## Computing `Gls_stats` recursively on the tree

As in Round 4, we work on the unrooted tree and associate a `Gls_stats` with each arc.
For an arc A->B, `Gls_stats(A->B)` holds sufficient statistics for the K tips on the B
side of the arc, with distances measured from B; the inner products use the K x K
precision matrix W = Sigma^{-1} for those tips.  To compute `Gls_stats(A->B)`, let B->C1,
..., B->Ck be the other arcs at B (excluding B->A), with k_i tips on Ci's side.  First, we
**shift** each `Gls_stats(B->Ci)` (a k_i x k_i problem), moving its reference point from
Ci to B (this accounts for the arc's mutation count and shared covariance).  Second, we
**combine** the shifted stats by element-wise addition: this is correct because the K x K
covariance matrix at B is block-diagonal in the k_i x k_i submatrices (tips in different
subtrees share no branches beyond B).  Both operations are O(1) on the six statistics.

Illustration for k = 2 (node B with two neighboring arcs B->C1 and B->C2):

```
        A
        |
        |  arc A->B
        |
        B
       / \
      /   \
     /     \
    C1      C2
   (k1     (k2
   tips)   tips)
```

Arc B->C1 has z1 observed mutations.  Arc B->C2 has z2 observed mutations.

Step 1 -- start (reference points at C1, C2):

    Gls_stats(B->C1) uses covariance Sigma_1  [k1 x k1]    distances measured from C1
    Gls_stats(B->C2) uses covariance Sigma_2  [k2 x k2]    distances measured from C2

Step 2 -- shift (reference points shifted to B):

    Let P1 = shift(Gls_stats(B->C1), z1).
    Let P2 = shift(Gls_stats(B->C2), z2).

    P1 uses covariance Sigma_1' = Sigma_1 + Var(Z1) * 1*1^T  [k1 x k1]    distances measured from B
    P2 uses covariance Sigma_2' = Sigma_2 + Var(Z2) * 1*1^T  [k2 x k2]    distances measured from B

    where Var(Zi) = zi + epsilon (see "Regularization" above).

Step 3 -- combine (= Gls_stats(A->B), measured from B):

    Gls_stats(A->B) = P1 + P2

                              / Sigma_1'      0    \
    uses covariance Sigma_B = |                    |    [(k1+k2) x (k1+k2)]
                              \    0      Sigma_2' /

We never form these matrices.  Shifting is O(1) via the Sherman-Morrison formula
(see below); combining is O(1) by element-wise addition of the six statistics.

### Base case: a single tip

For a tip X, measured at X itself, the distance from X to X is zero and the 1x1 covariance
matrix Sigma = [0] (no branches have been traversed).  Sigma cannot be inverted into W, so
the `Gls_stats` fields are strictly undefined.  We use `sum_1_W_1 == -1` as a sentinel for
"unshifted tip" and stash the tip's centered date dt_X in `sum_dt_W_1`:

    Gls_stats(P->X) = { -1, dt_X, 0, 0, 0, 0 }     (before shifting, measured from X)

The tip's contribution enters during shifting through the P->X arc (see below).
`shift()` detects the sentinel via `sum_1_W_1 == -1`, reads dt_X from `sum_dt_W_1`, and
directly produces the post-shift result rather than applying the Sherman-Morrison formula
(which requires an invertible Sigma).  Similarly, `combine()` CHECKs that neither operand
has `sum_1_W_1 < 0`, since an unshifted tip must be shifted before combining.

### Shifting through an arc: distance shift and Sherman-Morrison correction

Shifting through an arc with z observed mutations and variance sigma^2 = z + epsilon has
two cases, distinguished by `sum_1_W_1`:

**General case (sum_1_W_1 > 0).**  Suppose we have computed `Gls_stats` for a subtree of k
tips below node X, using the covariance matrix Sigma_old for those tips.  Now we cross the
arc above X, which has z mutations and variance sigma^2.  Two things happen:

1. Every tip's distance from the measurement point increases by z (the distance shift).
2. Every pair of tips in the subtree gains shared covariance sigma^2 from the common branch
   (the covariance update).

These two effects are independent and can be applied in either order.  It is illuminating to
treat them separately:

**Distance shift (= Round 4 OLS `prepend`).**  Adding z to every tip's distance
(m -> m + z*1) with the covariance matrix unchanged:

    sum_1_W_1'    = sum_1_W_1                                   (weights unchanged)
    sum_dt_W_1'   = sum_dt_W_1                                  (dates unchanged)
    sum_m_W_1'    = sum_m_W_1 + z * sum_1_W_1                   (each distance increases by z)
    sum_dt_W_dt'  = sum_dt_W_dt                                 (dates unchanged)
    sum_m_W_dt'   = sum_m_W_dt + z * sum_dt_W_1                 (cross-product shifts)
    sum_m_W_m'    = sum_m_W_m + 2*z*sum_m_W_1 + z^2 * sum_1_W_1 (squared distances shift)

These are exactly the Round 4 OLS `prepend` formulas, now expressed in the `Gls_stats`
notation.

**Sherman-Morrison covariance correction.**  The branch adds shared variance sigma^2 to all
tips in the subtree, so the covariance matrix gains a rank-1 term:

    Sigma_new = Sigma_old + sigma^2 * 1 * 1^T

We need to update the stats to use W_new = Sigma_new^{-1} instead of W_old = Sigma_old^{-1},
with distances unchanged.

Rather than re-inverting the full matrix, we use the Sherman-Morrison formula, a standard
identity for rank-1 updates to matrix inverses:

    (A + u * v^T)^{-1} = A^{-1} - A^{-1} * u * v^T * A^{-1} / (1 + v^T * A^{-1} * u)

With A = Sigma_old, u = v = sigma * 1:

    W_new = W_old - sigma^2 * (W_old * 1) * (1^T * W_old) / (1 + sigma^2 * 1^T * W_old * 1)

We do not need the full matrix W_old to apply this update.  Each component of `Gls_stats`
is an inner product of the form a^T W b where a, b are drawn from {dt, M, 1}.  For any such
pair, the Sherman-Morrison formula gives:

    a^T W_new b = a^T W_old b - sigma^2 * (a^T W_old 1) * (1^T W_old b) * gamma

where gamma = 1 / (1 + sigma^2 * sum_1_W_1) is a factor in [0, 1] that shrinks the
effective precision as shared covariance is added.  Every term on the right is an existing
`Gls_stats` component.  The correction to each component a^T W b is:

    -sigma^2 * (a^T W 1) * (1^T W b) * gamma

Concretely:

    sum_1_W_1'   = sum_1_W_1   - sigma^2 * sum_1_W_1  * sum_1_W_1  * gamma
    sum_dt_W_1'  = sum_dt_W_1  - sigma^2 * sum_dt_W_1 * sum_1_W_1  * gamma
    sum_m_W_1'   = sum_m_W_1   - sigma^2 * sum_m_W_1  * sum_1_W_1  * gamma
    sum_dt_W_dt' = sum_dt_W_dt - sigma^2 * sum_dt_W_1 * sum_dt_W_1 * gamma
    sum_m_W_dt'  = sum_m_W_dt  - sigma^2 * sum_m_W_1  * sum_dt_W_1 * gamma
    sum_m_W_m'   = sum_m_W_m   - sigma^2 * sum_m_W_1  * sum_m_W_1  * gamma

**Combined shift.**  Since the distance shift and Sherman-Morrison correction are independent,
they can be applied in either order.  Applying the distance shift first, then the covariance
correction (which is the more natural order, since it generalizes the OLS `prepend`), the
combined formulas for shifting through an arc with z mutations and variance sigma^2 are:

    gamma = 1 / (1 + sigma^2 * sum_1_W_1)

    sum_1_W_1'   = sum_1_W_1 * gamma
    sum_dt_W_1'  = sum_dt_W_1 * gamma
    sum_m_W_1'   = (sum_m_W_1 + z * sum_1_W_1) * gamma
    sum_dt_W_dt' = sum_dt_W_dt - sigma^2 * sum_dt_W_1^2 * gamma
    sum_m_W_dt'  = (sum_m_W_dt + z * sum_dt_W_1) - sigma^2 * sum_dt_W_1 * (sum_m_W_1 + z * sum_1_W_1) * gamma
    sum_m_W_m'   = (sum_m_W_m + 2*z*sum_m_W_1 + z^2*sum_1_W_1) - sigma^2 * (sum_m_W_1 + z*sum_1_W_1)^2 * gamma

These formulas match TreeTime's `propagate_averages`
([lines 268-280](https://github.com/neherlab/treetime/blob/17da461f9299706b99a90e622183772f00c5076d/treetime/treeregression.py#L268)).

**Tip special case (sum_1_W_1 == -1).**  When `sum_1_W_1 == -1`, the input represents an
unshifted tip X, with dt_X stashed in `sum_dt_W_1`.  The Sherman-Morrison formula does not
apply (Sigma is singular).  Instead, `shift` reads dt_X from `s.sum_dt_W_1` and directly
returns:

    { sum_1_W_1:    1/sigma^2,
      sum_dt_W_1:   dt_X / sigma^2,
      sum_m_W_1:    z / sigma^2,
      sum_dt_W_dt:  dt_X^2 / sigma^2,
      sum_m_W_dt:   z * dt_X / sigma^2,
      sum_m_W_m:    z^2 / sigma^2 }

This is the 1x1 GLS problem: one tip at distance z from the measurement point, with
Sigma = [sigma^2] and W = [1/sigma^2].  (This matches TreeTime's `propagate_averages`,
which also special-cases terminals.)

### Combining independent subtrees

At an internal node with children C1, C2, ..., the children's subtrees have independent
mutation counts (disjoint branches).  The joint covariance matrix is block-diagonal, so the
precision matrix W is also block-diagonal, and every inner product a^T W b decomposes as a
sum over blocks.  Therefore:

    Gls_stats_parent = Gls_stats_C1 + Gls_stats_C2 + ...

Element-wise addition, identical to the Round 4 OLS `combine` operation.

### Evaluating a candidate root along an arc

Once `Gls_stats` has been computed for every arc in both directions (see Implementation for
the pass structure), each arc A->B has `Gls_stats(A->B)` (B-side tips,
measured from B) and `Gls_stats(B->A)` (A-side tips, measured from A).  To evaluate a
candidate root at position k in {0, 1, ..., z}, where k mutations lie between A and the
root and z-k mutations lie between the root and B (same convention as Round 4), we shift
each side's stats to R:

    Gls_stats_A_side = shift(Gls_stats(B->A), k)        (shift by k muts, variance k + epsilon)
    Gls_stats_B_side = shift(Gls_stats(A->B), z - k)    (shift by z-k muts, variance z-k + epsilon)
    Gls_stats_root   = Gls_stats_A_side + Gls_stats_B_side

Then solve the GLS regression from `Gls_stats_root` in O(1).

Note that splitting the arc introduces two regularization terms (one per sub-arc), so the
total variance across the split is (k + epsilon) + (z-k + epsilon) = z + 2*epsilon, slightly
more than the original arc's z + epsilon.  TreeTime handles this differently: it linearly
interpolates the branch variance, using x*(z + epsilon) and (1-x)*(z + epsilon) for the two
halves, preserving the total at z + epsilon.  The difference is small (one extra epsilon).

In the Round 4 OLS case, chi^2 is quadratic in (continuous) k, giving a closed-form
minimizer along each arc.  In GLS, this no longer holds: the Sherman-Morrison correction
makes gamma depend on k, so chi^2 is not quadratic and there is no closed-form minimizer.
Instead, as in Round 4, we evaluate at every integer position k = 0, 1, ..., z along each
arc and pick the best.  (TreeTime uses a different approach: it treats the split position as
a continuous fraction x in [0, 1] along each branch, evaluates chi^2 on a 6-point grid, then
refines with `scipy.optimize.minimize_scalar`; see
[`_optimal_root_along_branch`](https://github.com/neherlab/treetime/blob/17da461f9299706b99a90e622183772f00c5076d/treetime/treeregression.py#L385).)

### Total cost

A bottom-up pass over the tree applies `shift` at each arc and `combine` (addition)
at each node, visiting each arc once: O(N).  A top-down pass computes the complement
(outgroup) statistics for each node, also O(N).  Given both subtree and complement stats for
any edge, the GLS regression at any candidate root position along that edge is O(1).

**Total: O(N)**, with no matrix inversions and no materialization of Sigma.

## Reduction to OLS

Setting sigma^2 = 0 for every branch (no variance weighting) gives gamma = 1, and every
Sherman-Morrison correction term vanishes.  The shift reduces to just the distance
shift:

    sum_1_W_1'   = sum_1_W_1                                    [= Round 4's n]
    sum_dt_W_1'  = sum_dt_W_1                                   [= Round 4's sum_dt]
    sum_m_W_1'   = sum_m_W_1 + z * sum_1_W_1                    [= Round 4's sum_m + z * n]
    sum_dt_W_dt' = sum_dt_W_dt                                  [= Round 4's sum_dt2 (constant)]
    sum_m_W_dt'  = sum_m_W_dt + z * sum_dt_W_1                  [= Round 4's sum_m_dt + z * sum_dt]
    sum_m_W_m'   = sum_m_W_m + 2*z*sum_m_W_1 + z^2 * sum_1_W_1  [= Round 4's sum_m2 + 2*z*sum_m + z^2*n]

These are exactly the Round 4 OLS `prepend` formulas.  Similarly, `combine` reduces to
element-wise addition (Round 4's `combine`).

Round 4 uses 5 statistics (n, sum_dt, sum_m, sum_m_dt, sum_m2) rather than 6.  The "missing"
statistic is sum_dt_W_dt (weighted sum of squared centered dates).  In the OLS case with
centered dates, sum_dt_W_dt = sum(dt_i^2) = N * Var_t, which is a global constant independent
of root position (because the unweighted sum of centered dates is always zero, regardless
of which tips are in the subtree).  With GLS weighting, sum_dt_W_dt varies as the precision
weights change with the candidate root, so it must be tracked per-subtree.

The Round 4 OLS implementation is therefore a conceptual special case of the GLS framework
developed here, obtained in the limit of zero branch variances.  In practice, the two are
separate functions (`ols_regression_root_utree` and `gls_regression_root_utree`) with
identical signatures, using different internal statistics (`Arc_stats` vs `Gls_stats`) and
different shift logic (distance-only vs distance + Sherman-Morrison).

Note that minimizing chi^2 across candidate root positions is not the same as maximizing R^2.
R^2 normalizes by Var(m), which changes with the root position (different roots yield
different root-to-tip distances).  chi^2 uses the precision-weighted residual sum directly
and is the correct criterion for maximum-likelihood root selection.  R^2 is computed once at
the winning root for reporting purposes only.

## Correspondence with TreeTime

TreeTime's naming convention uses "avg" (for "average") and "sq" (for "squared").  TreeTime
combines the distance shift and Sherman-Morrison correction into a single
`propagate_averages` function; we separate them conceptually above to make the connection to
the OLS `prepend` transparent, but the implementation uses the same combined formulas as
TreeTime.

Note on reference points: TreeTime uses the same convention as us -- subtree stats are
measured from the node itself, not from the parent.  For a node X with parent P,
grandparent G, and children C1, ..., Ck:

- `X.Q` = `Gls_stats(P->X)`: stats for tips below X, measured from X.  Computed bottom-up
  as `sum_i shift(Ci.Q, z_Ci_X)`.  Only set for internal nodes; for tips, `X.Q` is never
  assigned -- the tip's contribution enters via `propagate_averages`'s terminal special case
  ([line
  259](https://github.com/neherlab/treetime/blob/17da461f9299706b99a90e622183772f00c5076d/treetime/treeregression.py#L259)),
  which directly returns the 1x1 GLS result (the same special case we handle via the
  `sum_1_W_1 == -1` sentinel).
- `X.O` = `Gls_stats(X->P)`: stats for tips NOT below X, measured from P.  Let S1, ..., Sl
  be X's siblings (the other children of P).  Computed top-down as the sum of shifted
  sibling contributions `sum_j shift(Sj.Q, z_Sj_P)` plus `shift(P.O, z_PG)` (which brings
  P's outgroup stats from G to P).  When P is the root, the parent outgroup term is omitted
  (the root has no parent).  Set for all non-root nodes, including tips.
- `X.Qtot` = stats for ALL tips, measured from X.  Computed as
  `X.Q + shift(X.O, z_XP)`: the subtree stats (already at X) plus the outgroup stats
  shifted from P to X.  Only set for internal nodes (not tips, not the root).

TreeTime's tree is rooted with an explicit degree-2 root node R.  For R:
`R.Qtot = R.Q` (all tips, measured from R); `R.O` is never set.  For a child X of R in a
binary tree (sibling S): `X.O = shift(S.Q, z_SR)` (S-side tips measured from R, no parent
outgroup term).  In our Utree, R does not exist as a node -- the edge between X and S is a
single arc, and we evaluate candidate root positions along it using direct shifts from the
endpoints.

Note on regularization at the root edge: TreeTime's degree-2 root R splits a single
unrooted edge (X-S) into two branches (X-R and R-S), each receiving its own epsilon
regularization term.  The root edge therefore contributes 2\*epsilon of total variance
instead of 1\*epsilon.  Since the input root placement is arbitrary, this is an artifact
that subtly biases the stats at every candidate root position (the double-epsilon
propagates through the O values).  In our Utree, the X-S edge is a single arc with one
epsilon during the passes; two epsilon terms only appear when evaluating a specific
candidate root position, which is more consistent: every candidate is treated symmetrically,
regardless of the input root.  With TreeTime's setup, two differently rooted but otherwise
identical input trees can produce different "best" roots on output; this cannot happen with
our approach.

TreeTime uses `branch_variance = branch_length + 0.05 * total_bl / N` as the default
variance model -- the same regularization described in the "Regularization" section above.

| This plan        | TreeTime variable | TreeTime index | Line(s) |
|------------------|-------------------|----------------|---------|
| sum_1_W_1        | `Q[sii]`          | 5              | [4](https://github.com/neherlab/treetime/blob/17da461f9299706b99a90e622183772f00c5076d/treetime/treeregression.py#L4), [280](https://github.com/neherlab/treetime/blob/17da461f9299706b99a90e622183772f00c5076d/treetime/treeregression.py#L280) |
| sum_dt_W_1       | `Q[tavgii]`       | 0              | [268](https://github.com/neherlab/treetime/blob/17da461f9299706b99a90e622183772f00c5076d/treetime/treeregression.py#L268) |
| sum_m_W_1        | `Q[davgii]`       | 1              | [269](https://github.com/neherlab/treetime/blob/17da461f9299706b99a90e622183772f00c5076d/treetime/treeregression.py#L269) |
| sum_dt_W_dt      | `Q[tsqii]`        | 2              | [271](https://github.com/neherlab/treetime/blob/17da461f9299706b99a90e622183772f00c5076d/treetime/treeregression.py#L271) |
| sum_m_W_dt       | `Q[dtavgii]`      | 3              | [274](https://github.com/neherlab/treetime/blob/17da461f9299706b99a90e622183772f00c5076d/treetime/treeregression.py#L274) |
| sum_m_W_m        | `Q[dsqii]`        | 4              | [275](https://github.com/neherlab/treetime/blob/17da461f9299706b99a90e622183772f00c5076d/treetime/treeregression.py#L275) |
| `combine`        | `+=`              |                | [203](https://github.com/neherlab/treetime/blob/17da461f9299706b99a90e622183772f00c5076d/treetime/treeregression.py#L203) |
| `shift`          | `propagate_averages` |             | [235](https://github.com/neherlab/treetime/blob/17da461f9299706b99a90e622183772f00c5076d/treetime/treeregression.py#L235) |
| gamma            | `denom`           |                | [268](https://github.com/neherlab/treetime/blob/17da461f9299706b99a90e622183772f00c5076d/treetime/treeregression.py#L268) |
| bottom-up pass   | postorder loop in `_calculate_averages` | | [197](https://github.com/neherlab/treetime/blob/17da461f9299706b99a90e622183772f00c5076d/treetime/treeregression.py#L197) |
| top-down pass    | preorder loop in `_calculate_averages`  | | [206](https://github.com/neherlab/treetime/blob/17da461f9299706b99a90e622183772f00c5076d/treetime/treeregression.py#L206) |

## Fallback to midpoint rooting

The same fallback conditions as Round 4 apply.  `gls_regression_root_utree` falls back to
`midpoint_root_utree` when:

- Var_t = 0 (all tips have the same date -- detected before the tree passes).
- N <= 2 (not enough data for meaningful regression).
- No candidate position across any arc has alpha > 0 (no position gives a positive clock
  signal).

## Implementation

### Algorithm

The structure mirrors `regression_root_utree` (Round 4) closely: same bottom-up and
top-down DFS passes, same root selection and edge-splitting logic.  The differences are:

1. `Gls_stats` (6 fields) replaces `Arc_stats` (5 fields).
2. `shift` (distance shift + Sherman-Morrison correction) replaces `prepend` (distance shift
   only).
3. Root comparison uses chi^2_min instead of R^2.
4. An additional epsilon parameter is computed once and threaded through the passes.
5. Root evaluation is a separate pass (Pass 3) over all edges, rather than being interleaved
   with the top-down DFS.
6. R^2 is computed in a final O(N) pass after root selection, not during evaluation.

**Pass 0: global constants.**  Loop over all tips to compute `mean_t` and `Var_t`.  Compute
`epsilon` from the total mutation count:

```cpp
auto Nd = static_cast<double>(N);
auto epsilon = 0.05 * tree.count_deltas() / Nd;
```

If `Var_t <= 0` or `N <= 2`, fall back to `midpoint_root_utree`.

**`Gls_stats` struct and operations.**  Defined locally inside `gls_regression_root_utree`
(same pattern as `Arc_stats` in Round 4):

```cpp
struct Gls_stats {
  double sum_1_W_1 = 0.0;
  double sum_dt_W_1 = 0.0;
  double sum_m_W_1 = 0.0;
  double sum_dt_W_dt = 0.0;
  double sum_m_W_dt = 0.0;
  double sum_m_W_m = 0.0;
};
```

`combine(a, b)`:  Element-wise addition of all six fields (identical pattern to Round 4).
CHECKs that `a.sum_1_W_1 >= 0` and `b.sum_1_W_1 >= 0` (rejects unshifted tip sentinels,
which have `sum_1_W_1 == -1`; allows the all-zeros identity element).

`shift` is a closure (lambda) capturing `epsilon`.  Its signature is
`shift(Gls_stats s, int z)`:  it computes `sigma_sq = z + epsilon` internally.  In the tip
special case (`s.sum_1_W_1 == -1`), it reads dt_X from `s.sum_dt_W_1` and returns the 1x1
GLS result directly.  Otherwise, it applies the combined distance-shift + Sherman-Morrison
formulas from the "Combined shift" section above.

**Pass 1: bottom-up (post-order) DFS.**  Same loop structure as Round 4.  On each leaving
event for arc X->P, compute `gls_stats[mate(X->P)]`:

- Tip X: `gls_stats[P->X] = {-1, dt_of(X), 0, 0, 0, 0}` (sentinel with dt stashed).
- Inner node X: combine shifted stats of all child arcs:
  ```cpp
  auto combined = Gls_stats{};
  for (auto a : tree.nodes[X].arcs) {
    if (a != k_no_arc && a != arc_X_to_P) {
      combined = combine(combined, shift(gls_stats[a], tree.count_arc_deltas(a)));
    }
  }
  ```

**Pass 2: top-down (pre-order) DFS.**  Same loop structure as Round 4.  On each entering
event for arc P->X, compute `gls_stats[mate(P->X)]`:

- Tip P: `gls_stats[X->P] = {-1, dt_of(P), 0, 0, 0, 0}` (sentinel with dt stashed).  This
  happens at most once, when the DFS starts at a tip.
- Inner node P: combine shifted stats of all outgoing arcs from P other than P->X (same
  pattern as Pass 1).

After Passes 1 and 2, `gls_stats[a]` is set for every arc `a` in the tree.

**Pass 3: root evaluation.**  Loop over all edges (arc pairs A->B / B->A).  For each edge,
evaluate chi^2 at every position k = 0, ..., D:

```cpp
auto D = tree.count_arc_deltas(arc_A_to_B);
for (auto k = 0; k <= D; ++k) {
  auto A_side = shift(gls_stats[arc_B_to_A], k);
  auto B_side = shift(gls_stats[arc_A_to_B], D - k);
  auto root_stats = combine(A_side, B_side);

  // GLS solution
  auto denom = root_stats.sum_dt_W_dt * root_stats.sum_1_W_1
             - root_stats.sum_dt_W_1 * root_stats.sum_dt_W_1;
  if (denom <= 0.0) { continue; }
  auto alpha = (root_stats.sum_m_W_dt * root_stats.sum_1_W_1
              - root_stats.sum_m_W_1 * root_stats.sum_dt_W_1) / denom;
  if (alpha <= 0.0) { continue; }
  auto beta = (root_stats.sum_m_W_1 - alpha * root_stats.sum_dt_W_1)
            / root_stats.sum_1_W_1;

  auto chi2 = root_stats.sum_m_W_m - alpha * root_stats.sum_m_W_dt
            - beta * root_stats.sum_m_W_1;

  // Track best (lowest chi^2), with tie-breaking vector
  if (chi2 < best_chi2) { best_chi2 = chi2; best_candidates.clear(); }
  if (chi2 == best_chi2) { best_candidates.push_back({arc_A_to_B, k}); }
}
```

where `best_chi2` is initialized to `std::numeric_limits<double>::infinity()` before the
pass begins.

**Root selection and edge splitting.**  Pick a winner uniformly at random from
`best_candidates`.  Recompute `root_stats` at the winning position (needed for `alpha`,
`beta`, and lambda/t_MRCA), then allocate the root node and split the edge, same as Round 4:

```cpp
auto [best_arc, best_k] = best_candidates[idx];
auto best_D = tree.count_arc_deltas(best_arc);
auto best_root_stats = combine(
    shift(gls_stats[tree.mate(best_arc)], best_k),
    shift(gls_stats[best_arc], best_D - best_k));
```

**R^2 computation (final pass).**  After the root is placed, a single DFS from the root
computes unweighted root-to-tip mutation counts and evaluates R^2:

```cpp
auto cur_dist = 0;
auto sum_m = 0.0, sum_m2 = 0.0, sum_m_dt = 0.0;
for (auto [arc, direction] : annotated_arc_euler_tour(tree, R)) {
  auto arc_deltas = tree.count_arc_deltas(arc);
  if (direction == Arc_direction::entering) {
    cur_dist += arc_deltas;
    if (tree.is_tip(tree.target(arc))) {
      auto m = static_cast<double>(cur_dist);
      auto dt = dt_of(tree.target(arc));
      sum_m += m;  sum_m2 += m * m;  sum_m_dt += m * dt;
    }
  } else {
    cur_dist -= arc_deltas;
  }
}
auto mean_m = sum_m / Nd;
auto cov_mt = sum_m_dt / Nd;
auto var_m = sum_m2 / Nd - mean_m * mean_m;
auto r2 = (var_m > 0.0 && var_t > 0.0) ? (cov_mt * cov_mt) / (var_m * var_t) : 0.0;
```

Return `Rooting_info{R, Rooting_method::regression, r2, alpha, mean_t - beta / alpha}`.

### Integration

**Rename.**  Rename `regression_root_utree` to `ols_regression_root_utree` in both
`utree.h` and `utree.cpp`.

**New function.**  Declare in `utree.h`:

```cpp
auto gls_regression_root_utree(Utree& tree, const std::vector<Tip_desc>& tip_descs,
                               absl::BitGenRef bitgen)
    -> Rooting_info;
```

Same signature as `ols_regression_root_utree`.

**Pipeline change.**  In `build_initial_phylo_tree`, replace:

```cpp
auto rooting_info = regression_root_utree(utree, tip_descs, bitgen);
```

with:

```cpp
auto rooting_info = gls_regression_root_utree(utree, tip_descs, bitgen);
```

The old OLS function remains available for comparison; switching back is a one-line change.

### OLS cleanup for structural parallelism

To make the OLS and GLS functions easy to compare and maintain side by side, the OLS
function (`ols_regression_root_utree`) is refactored to match the GLS function's structure:

- **Rename `Arc_stats` -> `Ols_stats`** and **`prepend` -> `shift`**, mirroring the GLS
  naming (`Gls_stats`, `shift`).  `shift` is the right name: both operations shift the
  measurement reference point through an arc; the OLS version is the distance-only special
  case of the GLS version.
- **Move struct and lambdas after Pass 0**, so both functions have the same top-level
  ordering: early return, Pass 0 globals, struct, shift, combine, stats vector, passes.
- **Split Pass 2 into Pass 2 (top-down DFS) + Pass 3 (root evaluation)**, matching the
  GLS pass structure.  The original OLS interleaved the top-down stats computation with R^2
  evaluation in a single pass; splitting them is a negligible performance cost (one extra
  loop over arc pairs) and makes the two functions structurally parallel.
- **Use `arc_A_to_B` / `arc_B_to_A` naming** in the evaluation loop (matching GLS), with
  `const auto&` aliases for the two sides' stats.

The tip base case remains structurally different between OLS and GLS (OLS uses a direct
`{.n = 1, .sum_dt = dt_of(X), ...}` initializer; GLS uses the `sum_1_W_1 == -1` sentinel
with dt stashed in `sum_dt_W_1`), since this is a genuine algorithmic difference.

### Files to modify

| File | Action | What |
|------|--------|------|
| `core/utree.h` | Modify | Rename `regression_root_utree` -> `ols_regression_root_utree`; declare `gls_regression_root_utree` |
| `core/utree.cpp` | Modify | Rename existing function; implement `gls_regression_root_utree`; update `build_initial_phylo_tree` call |
| `tests/utree_tests.cpp` | Modify | Rename existing regression tests to `ols_`; add parallel `gls_` tests |

No changes to `cmdline.h` or `cmdline.cpp` -- the `rooting_hook` interface is unchanged.

### Verification

1. **Parallel test suites**: the four existing OLS regression tests (three_tips_known_rate,
   same_dates_fallback, two_tips_fallback, star_topology) are renamed to `ols_regression_root_*`
   and duplicated as `gls_regression_root_*` tests with identical structure.

2. **Degenerate cases**: same as Round 4 (same-date fallback, 2-tip fallback, no-positive-
   clock-signal fallback).

3. **End-to-end**: Run `delphy_ui` with the GLS rooting on standard demo datasets.  Verify
   `assert_phylo_tree_matches_tip_descs` passes and results are sensible.

4. **Build**: `cmake --build build/debug --target tests -j$(nproc) && build/debug/tests/tests --gtest_filter='Utree*'`

## Post-implementation findings: GLS gives poor rooting results

### Observed behavior

After implementing `gls_regression_root_utree` as described above, we tested it on the
h5n1-andersen-2025 dataset (3335 tips) from our delphy-2026-paper-data benchmarks.  The
results were strikingly poor:

| Method                                        | lambda (mut/yr) |   R^2   |   t_MRCA   | sum_1_W_1 |
|-----------------------------------------------|----------------:|--------:|------------|----------:|
| OLS                                           |           37.84 |  0.92   | 2023-12-16 |   3335    |
| GLS (sigma^2 = z + epsilon, epsilon ~= 0.255) |            1.6  | ~0      | 2011-10-30 |     ~0.09 |
| GLS (sigma^2 = z + 1)                         |           10.1  |  0.0004 | 2022-11-22 |     ~0.09 |

The `sum_1_W_1` column is the total precision weight `1^T W 1` (where `W = Sigma^{-1}`).
In OLS with `W = I`, this equals N.  In GLS, the Sherman-Morrison corrections shrink it
dramatically.  When `1^T W 1` is near zero, the precision matrix W nearly annihilates the
constant vector **1** -- meaning the data provides almost no constraint on the regression
intercept.  The slope estimate also degrades because its denominator,
`sum_dt_W_dt * sum_1_W_1 - sum_dt_W_1^2`, shrinks toward zero.

### Root cause

The Sherman-Morrison rank-1 update `W' = W - sigma^2 * (W*1)(1^T*W) / (1 + sigma^2*1^T*W*1)` shrinks
`sum_1_W_1` by a factor `gamma = 1/(1 + sigma^2*sum_1_W_1)` at each branch.  After
traversing many internal branches, the cumulative shrinkage drives `sum_1_W_1` toward zero.
This is mathematically correct -- it accounts for the fact that tips in the same clade are not
independent observations -- but it means GLS regression has almost no power to distinguish
candidate root positions.

We tested two variance models:

1. **sigma^2 = z + epsilon** (epsilon = 0.05 * total_muts / N ~= 0.255): the original plan's model,
   analogous to TreeTime's default variance model.
2. **sigma^2 = z + 1**: derived from the posterior mean of branch length tau given z mutations
   under a Poisson model with flat prior, tau|z ~ Gamma(z+1, lambda_0), giving
   E[Var(z|tau)|z] = lambda_0 * E[tau|z] = lambda_0 * (z+1)/lambda_0 = z + 1.  This is lambda_0-independent
   and regularizes zero-mutation branches more aggressively.

Neither model produced usable rooting results.

### How TreeTime handles this

TreeTime's `--covariation` flag triggers a multi-step process:

1. **First pass** (`run(max_iter=0)`): reroots with `var=0` for all internal branches
   (effectively OLS), finding a good root position.
2. **Time tree fit**: after rerooting, fits node times and computes `clock_length` for each
   branch (the regression-predicted branch length).
3. **Second pass** (`reroot()`): reroots again using GLS with
   `var = (clock_length + tip_slack^2 * om) * om` where `om = 1/L` and `tip_slack = 10`.
   Because the variance model now uses the smoothed `clock_length` rather than raw
   `mutation_length`, and because the root is already near-optimal, this GLS pass makes
   only minor adjustments.

The key insight is that **the GLS regression result is root-position-independent** for a
fixed variance model (the Sherman-Morrison formula gives exact updates, so all candidate
root positions are evaluated regardless of the starting root).  TreeTime's multi-step
approach works not because OLS provides a better starting root, but because the `clock_length`-
based variance model used in later GLS passes produces a less aggressive covariance
correction than the raw `mutation_length`-based model.

### Decision: use OLS for rooting

Given that:

- OLS gives excellent rooting results (R^2 = 0.92, sensible t_MRCA and rate)
- GLS with any fixed variance model we tried collapses the effective sample size
- TreeTime's GLS success relies on an iteratively refined variance model that requires a
  full time-tree fit between rerooting passes, which is far more complex than we need
- The GLS implementation is mathematically correct (verified against direct matrix inversion
  on small trees and against TreeTime's `propagate_averages` formulas)

we decided to **use OLS (`ols_regression_root_utree`) for tree initialization rooting**.

The GLS function `gls_regression_root_utree` is retained in the codebase for potential
future use (e.g., if an iterative variance refinement scheme is implemented), but it is not
called from `build_initial_phylo_tree`.
