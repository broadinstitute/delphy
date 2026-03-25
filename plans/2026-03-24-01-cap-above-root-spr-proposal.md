# Cap Above-Root SPR Proposal Time to Prevent Overflow and OOM

## Problem

When Delphy runs on a dataset with very few or zero mutations (e.g., WCSS
`04_free_exp_pop/sims/sim_048`, which has `num_mutations = 0` and
`mu = 7.4e-7 subst/site/year`), the SPR above-root proposal (implemented in
`Spr_study::pick_time_in_region()` in `core/spr_study.cpp`) can generate
proposed root times astronomically far in the past, triggering two failures:

1. **OOM:** `Very_scalable_coalescent_prior_part::ensure_space()` tries to
   allocate vectors indexed by `cell_for(t, ...)`, which returns a huge cell
   index for extreme times, exhausting memory.

2. **Integer overflow:** `cell_for()` returns `static_cast<int>(std::floor(...))`
   — when the argument exceeds `INT_MAX`, this is undefined behavior in C++.
   In practice it wraps to `INT_MIN` (`-2147483648`), causing:
   ```
   Check failed: cell_start <= cell_end (86 vs. -2147483648)
   ```
   at `very_scalable_coalescent.cpp:397`.

### Root cause chain

1. The simulated dataset has 0 mutations (mu drawn from prior was extremely low).
2. With 0 mutations, the genetic prior has no per-mutation terms but retains
   the factor `exp[-mu L T]`, where `T` is total branch length and `L` is
   genome size.  The posterior's dependence on `mu` is then dominated by this
   factor and the prior on `mu`, giving a Gibbs full conditional that drives
   `mu` very close to 0.
3. `lambda_X = mu * L` (per-lineage total mutation rate) approaches 0.
4. In `pick_time_in_region()` (spr_study.cpp:389-412), the above-root
   proposal samples a time parameter
   `s = (t_X - t) + (t_S - t)` (sum of the two branch lengths from the
   proposed new root at time `t` to its children X and S, so
   `t = (t_X + t_S - s) / 2`) from a truncated Gamma with rate parameter
   `lambda_X * f`, where `f ≤ 1` is the annealing factor (a tempered
   likelihood exponent that leaves slack in the MH ratio for imperfect
   proposals; currently hard-coded to 0.8).  The inverse-CDF step divides
   by `lambda_X * f`:
   ```cpp
   auto rand_s = boost::math::gamma_q_inv(f*m + 1, rand_Q) / (lambda_X * f);
   ```
   When `lambda_X * f ≈ 0`, `rand_s` → ∞, producing `rand_t` → -∞.

   (Note: when `lambda_X f s_min` is very large, `Q_max = Q(fm+1,
   lambda_X f s_min)` can underflow to a tiny subnormal or to exactly
   0.0.  If Q_max is exactly 0.0, `log(Q_max) = -inf` makes the region
   weight `-inf`, so the region is never selected and the sampling code
   is never reached.  If Q_max is a nonzero subnormal,
   `absl::Uniform(OpenOpen, 0.0, Q_max)` correctly returns a nonzero
   value and `gamma_q_inv` handles it — so this path does not crash.
   The only confirmed failure mode is the `lambda_X ≈ 0` overflow above.)

   The fix below addresses this by capping `s_max`, with a power-law
   fallback that avoids `gamma_q_inv()` and the division by `lambda_X*f`
   entirely when `lambda_X` is small.

5. This extreme `new_t_P` is passed to
   `calc_delta_partial_log_prior_after_displace_coalescence(old_t_P, new_t_P)`,
   which calls `ensure_space(new_t_P)` and `cell_for(new_t_P, ...)`,
   triggering OOM or integer overflow.

## Current above-root proposal

The SPR proposal for reattachment above the root samples a time parameter
`s = (t_X - t) + (t_S - t)` (sum of the two branch lengths from the new root
to its children X and S) from a one-sided truncated Gamma:

```
p(s) ∝ s^{fm} exp(-λ_X f s),   s ≥ s_min = |t_X - t_S|
```

with no upper bound on `s`.  The region weight is obtained by integrating the
full unnormalized density `p_root(s) = 0.5 · f · λ_X · (mu · s / 3)^{fm} ·
exp(-λ_X f s)` over `[s_min, ∞)`.  Substituting `x = λ_X f s` converts the
integral into a regularized upper incomplete gamma function `Q(a, z)`:

```
W_R = ∫_{s_min}^{∞} p_root(s) ds
    = 0.5 · f · λ_X · (mu/3)^{fm} · ∫_{s_min}^{∞} s^{fm} exp(-λ_X f s) ds
    = 0.5 · (mu/(3λ_X f))^{fm} · Γ(fm+1) · Q(fm+1, λ_X f s_min)
```

The prefactor `0.5 · f · λ_X · (mu/3)^{fm}` ensures this weight is on the
same scale as the below-root region weights (which include `f · λ_X` and
`(mu/3)^{fm}` factors from the same unnormalized density).

Sampling uses inverse-CDF via `Q`:

```cpp
auto Q_min = 0.0;
auto Q_max = boost::math::gamma_q(f*m + 1, lambda_X*f*s_min);
auto rand_Q = absl::Uniform(absl::IntervalOpenOpen, bitgen, Q_min, Q_max);
auto rand_s = boost::math::gamma_q_inv(f*m + 1, rand_Q) / (lambda_X*f);
auto rand_t = 0.5 * (t_X + t_S - rand_s);
```

When `lambda_X → 0`, the `/ (lambda_X*f)` produces overflow.

## Proposed fix: cap s_max at 10× current tree span

### Design

Add an upper bound `s_max` to the above-root proposal, making it a
**doubly-truncated Gamma** on `[s_min, s_max]`.  The cap limits the proposed
new root time to be at most `10 · tree_span` earlier than
`t_early = min(t_X, t_S)`, the earlier of the two children of the proposed
new root:

```
t_new_root ≥ t_early - 10 · tree_span
```

where `tree_span = t_max_tip - t_early` is the time span from the earlier
child to the latest tip in the tree.

Since `t = (t_X + t_S - s) / 2` and the distance from the earlier child to
the new root is `t_early - t = (s - s_min) / 2`, this translates to:

```
(s - s_min) / 2 ≤ 10 · tree_span
s_max = s_min + 20 · tree_span
```

where `tree_span = t_max_tip - min(t_X, t_S)`.

The factor of 20 in s-space corresponds to 10× in root-time-space, because
`s` is the *sum* of two branch lengths (the `s ↔ t` mapping has a factor
of 2).

### Why 10×?

- **Normal operation (λ_X not tiny):** The Gamma distribution's exponential
  tail `exp(-λ_X f s)` already concentrates nearly all mass near `s_min`.
  The cap at 10× the tree span is far beyond the effective support of the
  distribution and has negligible effect on proposal probabilities.

- **Pathological case (λ_X → 0):** The Gamma's peak moves to
  `s ≈ fm / (λ_X f)`, which diverges as `λ_X → 0`.  The distribution
  still has infinite support (any `s` is possible), but most of its mass
  concentrates around the peak — so sampled values of `s` can be
  astronomically large, causing numerical overflow in `cell_for()` and
  `ensure_space()`.  The cap keeps `rand_s` bounded regardless of how
  small `λ_X` is.

- **MCMC correctness:** A root 10× the tree span into the past is already
  extremely implausible under any reasonable coalescent prior.  The coalescent
  prior penalizes such deep roots heavily, so these proposals would be
  rejected with near-certainty anyway.  The cap simply avoids computing
  (and overflowing on) proposals that have zero practical chance of acceptance.

### Effect on detailed balance

The MH acceptance ratio remains correct as long as:
1. The region weight uses the doubly-truncated integral.
2. The conditional density `p(t | region)` is the doubly-truncated Gamma density.
3. `log_alpha_in_region()` computes the correct density for the doubly-truncated case.

All three use the same `Q(a, x)` functions — the only change is that `Q_min`
becomes `Q(fm+1, λ_X f s_max)` instead of 0, and the normalization changes
accordingly.  Since both the forward and reverse SPR studies apply the same
cap, detailed balance is preserved.

### Numerical stability for small λ_X f

There are two numerical regimes to consider:

- **Both Q values close to 1** (`λ_X f s_max < 0.01`): Both
  `Q(fm+1, λ_X f s_min)` and `Q(fm+1, λ_X f s_max)` are close to 1, and
  their difference suffers from catastrophic cancellation.  In this regime,
  `exp(-λ_X f s) ≈ 1` across `[s_min, s_max]`, and the distribution is
  effectively a pure power law `s^{fm}`.  We handle this with a power-law
  fallback (see below).

- **Both Q values close to 0** (`λ_X f s_min ≫ fm+1`): Both Q values are
  in the exponential tail.  When both are representable (nonzero), Q decays
  roughly as `x^{a-1} exp(-x) / Γ(a)` for `x ≫ a`, and since `x_max > x_min`
  we have `Q(x_min) ≫ Q(x_max)`, so `Q(x_min) - Q(x_max) ≈ Q(x_min)` with
  no cancellation.  Boost's `gamma_q` is accurate in this tail, so no special
  handling is needed.  If both Q values underflow to exactly 0.0, then
  `log(0 - 0) = -inf`, the weight is `-inf`, and the region is never
  selected — also harmless.  Physically, this regime corresponds to `s_min`
  already being well past the Gamma peak — the above-root weight is tiny,
  and the region is rarely (or never) selected.

#### CDF derivation for power-law sampling

Let `α = fm + 1` (always ≥ 1 since `f > 0` and `m ≥ 0`).  The unnormalized
density is `s^{fm} = s^{α-1}` on `[s_min, s_max]`.

The CDF is:

```
F(s) = ∫_{s_min}^{s} u^{α-1} du  /  ∫_{s_min}^{s_max} u^{α-1} du
     = (s^α - s_min^α) / (s_max^α - s_min^α)
```

To sample via inverse-CDF, set `U = F(s)` and solve for `s`:

```
U = (s^α - s_min^α) / (s_max^α - s_min^α)
s^α = s_min^α + U · (s_max^α - s_min^α)
s = (s_min^α + U · (s_max^α - s_min^α))^{1/α}
```

Since `α ≥ 1 > 0` and `s_min ≥ 0`, the function `x^{1/α}` is monotonically
increasing.  When `U = 0`, `s = (s_min^α)^{1/α} = s_min`; when `U = 1`,
`s = (s_max^α)^{1/α} = s_max`.  So `s ∈ [s_min, s_max]` for all
`U ∈ [0, 1]`, as required.

#### Weight and density

The full unnormalized above-root density (including the prefactors that make
it comparable to below-root region weights) is:

```
p_root(s) = 0.5 · f · λ_X · (mu · s / 3)^{fm} · exp(-λ_X f s)
```

In the power-law regime (`exp(-λ_X f s) ≈ 1`), the weight integral becomes:

```
W_R = 0.5 · f · λ_X · (mu/3)^{fm} · ∫_{s_min}^{s_max} s^{fm} ds
    = 0.5 · f · λ_X · (mu/3)^{fm} · (s_max^α - s_min^α) / α
```

The factor `f · λ_X · (mu/3)^{fm}` is essential: it is the same prefactor
that appears in below-root region weights (via the `f*lambda_X*(t_max - t_min)`
and `(mu*(t_X - t')/3)^{fm}` terms), ensuring correct relative weighting
between above-root and below-root proposals.

The conditional density (for computing the MH proposal ratio) is simply the
normalized power law:

```
p(s | R) = α · s^{α-1} / (s_max^α - s_min^α)
p(t | R) = 2 · p(s | R)
```

We switch to the power-law fallback when `λ_X f s_max < 0.01`.

## Implementation

The guiding principle: when `lambda_X` is not tiny, the new code should
produce results that are nearly identical to the current code, differing
only by the negligible truncation at `s_max`.  The power-law fallback
activates only when `lambda_X f s_max < 0.01`.

The comments in `spr_study.cpp` should be updated to reflect the changes:
the extensive block comment (lines 261-299) explaining the above-root weight
derivation, the smaller block comments in `pick_time_in_region()` (lines
389-396) and `log_alpha_in_region()` (lines 440-444).

### Step 1: `core/spr_study.h` — Add `t_max_tip` member

Add `double t_max_tip;` to the `Spr_study` struct inputs and update the
constructor signature:

```cpp
Spr_study(
    Spr_study_builder&& builder,
    double lambda_X,
    double annealing_factor,
    double t_X,
    double t_max_tip);  // NEW: time of latest tip (including tip-date uncertainty)
```

### Step 2: `core/spr_study.cpp` — Weight calculation (constructor, lines 317-326)

In the `else` branch (above-root region), compute `s_max` and use the
doubly-truncated weight:

```cpp
} else {
    auto t_S = tree->at(region.branch).t;
    auto s_min = std::abs(t_X - t_S);

    // Cap: limit proposed root to at most 10× the tree span earlier
    // than min(t_X, t_S).  This prevents overflow in cell_for() and
    // ensure_space() when lambda_X is tiny.  The factor of 20 in
    // s-space corresponds to 10× in root-time-space because
    // s = (t_X - t) + (t_S - t), so (s - s_min)/2 is the distance
    // from the closer child to the new root.
    auto t_early = std::min(t_X, t_S);
    auto tree_span = t_max_tip - t_early;
    CHECK_GE(tree_span, 0.0);
    auto s_max = s_min + 20.0 * tree_span;

    auto x_min = lambda_X * f * s_min;
    auto x_max = lambda_X * f * s_max;

    if (x_max < 0.01) {
        // Power-law regime (lambda_X tiny): exp(-lambda_X f s) ≈ 1, so
        //   p_root(s) ≈ 0.5 f lambda_X (mu s/3)^{fm}
        //   W_R = 0.5 f lambda_X (mu/3)^{fm} (s_max^α - s_min^α) / α
        // For the log of (s_max^α - s_min^α), factor out s_max^α and use
        // log1p for numerical stability when s_min/s_max is close to 1:
        //   log(s_max^α - s_min^α) = α log(s_max) + log1p(-(s_min/s_max)^α)
        auto alpha = f * m + 1;
        region.log_W_over_Wmax =
            -std::numbers::ln2
            + std::log(f * lambda_X)
            + f*m*std::log(mu / 3)
            + alpha*std::log(s_max) + std::log1p(-std::pow(s_min / s_max, alpha))
            - std::log(alpha);
    } else {
        // Standard doubly-truncated Gamma.  When x_max is large (normal case),
        // Q(fm+1, x_max) ≈ 0, and the result is nearly identical to the
        // one-sided truncation used in Delphy 1.3.1 and earlier.
        region.log_W_over_Wmax =
            -std::numbers::ln2
            + f*m*std::log(mu / (3 * lambda_X * f))
            + std::lgamma(f*m + 1)
            + std::log(boost::math::gamma_q(f*m + 1, x_min)
                     - boost::math::gamma_q(f*m + 1, x_max));
    }
}
```

### Step 3: `core/spr_study.cpp` — Time sampling in `pick_time_in_region()` (lines 389-412)

Update to use doubly-truncated sampling with the power-law fallback:

```cpp
auto t_S = tree->at(region.branch).t;
auto s_min = std::abs(t_X - t_S);
auto t_early = std::min(t_X, t_S);
auto tree_span = t_max_tip - t_early;
CHECK_GE(tree_span, 0.0);
auto s_max = s_min + 20.0 * tree_span;
auto x_min = lambda_X * f * s_min;
auto x_max = lambda_X * f * s_max;

double rand_s;
if (x_max < 0.01) {
    // Power-law fallback: sample s from s^{fm} on [s_min, s_max]
    // via inverse CDF: s = (s_min^a + U (s_max^a - s_min^a))^{1/a}
    auto alpha = f * m + 1;
    auto U = absl::Uniform(absl::IntervalOpenOpen, bitgen, 0.0, 1.0);
    auto s_min_a = std::pow(s_min, alpha);
    auto s_max_a = std::pow(s_max, alpha);
    rand_s = std::pow(s_min_a + U * (s_max_a - s_min_a), 1.0 / alpha);
} else {
    // Doubly-truncated Gamma via inverse CDF.  When x_max is large
    // (normal case), Q_min ≈ 0 and this is nearly identical to the
    // one-sided truncation used in Delphy 1.3.1 and earlier.
    auto Q_min = boost::math::gamma_q(f*m + 1, x_max);
    auto Q_max = boost::math::gamma_q(f*m + 1, x_min);
    auto rand_Q = absl::Uniform(absl::IntervalOpenOpen, bitgen, Q_min, Q_max);
    try {
        rand_s = boost::math::gamma_q_inv(f*m + 1, rand_Q) / (lambda_X*f);
    } catch (const std::exception& e) {
        // We've seen gamma_q_inv throw in the wild but have been unable to
        // reproduce it.  Log all inputs to help diagnose if it recurs.
        std::cerr << absl::StreamFormat(
            "gamma_q_inv threw: %s\n"
            "  f*m+1=%.17g, rand_Q=%.17g, lambda_X=%.17g, f=%.6g\n"
            "  Q_min=%.17g, Q_max=%.17g, x_min=%.17g, x_max=%.17g\n"
            "  s_min=%.17g, s_max=%.17g, m=%d, t_X=%.17g, t_S=%.17g\n",
            e.what(), f*m+1, rand_Q, lambda_X, f,
            Q_min, Q_max, x_min, x_max,
            s_min, s_max, m, t_X, t_S);
        throw;
    }
}
auto rand_t = 0.5 * (t_X + t_S - rand_s);
CHECK_GE(rand_t, region.t_min - 1e-6);
CHECK_LE(rand_t, region.t_max + 1e-6);
rand_t = std::max(region.t_min, std::min(region.t_max, rand_t));
return rand_t;
```

### Step 4: `core/spr_study.cpp` — Density calculation in `log_alpha_in_region()` (lines 428-458)

Update to use the doubly-truncated density with the power-law fallback:

```cpp
auto t_S = tree->at(region.branch).t;
auto s_min = std::abs(t_X - t_S);
auto t_early = std::min(t_X, t_S);
auto tree_span = t_max_tip - t_early;
CHECK_GE(tree_span, 0.0);
auto s_max = s_min + 20.0 * tree_span;
auto x_min = lambda_X * f * s_min;
auto x_max = lambda_X * f * s_max;
auto s = t_X - t + t_S - t;
CHECK_GE(s, s_min - 1e-6);
if (s > s_max + 1e-6) {
    // The query point t is beyond the cap.  This can happen when the
    // reverse study evaluates the old root position: the old root may
    // be deeper than s_max allows (e.g., when min(t_X, t_S) is close
    // to t_max_tip, making tree_span small).  The density is 0 here
    // (outside the proposal's support), so log_alpha = -inf and the
    // MH ratio rejects the move.  This is correct — detailed balance
    // is preserved; the proposal is simply wasted.
    return -std::numeric_limits<double>::infinity();
}

if (x_max < 0.01) {
    // Power-law density:
    //   p(s | R) = a s^{a-1} / (s_max^a - s_min^a),  where a = fm+1
    //   p(t | R) = 2 p(s | R)
    // For log(s_max^a - s_min^a), factor out s_max^a and use log1p:
    //   log(s_max^a - s_min^a) = a log(s_max) + log1p(-(s_min/s_max)^a)
    auto alpha = f * m + 1;
    return log_p_region +
        std::numbers::ln2 +
        std::log(alpha) +
        (alpha - 1) * std::log(s) +
        -alpha*std::log(s_max) +
        -std::log1p(-std::pow(s_min / s_max, alpha));
} else {
    // Doubly-truncated Gamma density.  When x_max is large (normal case),
    // Q(fm+1, x_max) ≈ 0, and the result is nearly identical to the
    // one-sided truncation used in Delphy 1.3.1 and earlier.
    return log_p_region +
        std::numbers::ln2 +
        std::log(lambda_X*f) +
        f*m * std::log(lambda_X*f*s) +
        -lambda_X*f*s +
        -std::lgamma(f*m+1)
        -std::log(boost::math::gamma_q(f*m+1, x_min)
                - boost::math::gamma_q(f*m+1, x_max));
}
```

### Step 5: Guard against `lambda_X == 0` in `spr1_move()`

Add an early return in `spr1_move()` (subrun.cpp), right after node X is
selected (we need `X` to look up `lambda_i_.at(X)`), when
`lambda_i_.at(X) == 0.0`.  With `lambda_X = 0`, all region weights contain
`log(f * lambda_X) = -∞`, making `sum_W_over_Wmax = 0` and
`pick_nexus_region()` undefined.  In practice `lambda_X = mu * L` is never
exactly 0, but a defensive guard is cheap.

### Step 6: Plumb `t_max_tip` from `Run` to `Spr_study`

`Run` already computes `calc_max_tip_time(tree_)` (which uses `t_max`,
accounting for tip-date uncertainty) in several places.  Make `t_max_tip_`
a mutable derived quantity of `Run` (like `num_muts_`, `log_G_`, etc.),
computed in `recalc_derived_quantities()` and accessible via a
`t_max_tip()` accessor.  This ensures the value is recomputed if the tree
is modified externally between calls to `do_mcmc_steps()` (callers signal
this via `Run::tree_modified()`).  Replace the explicit
`calc_max_tip_time(tree_)` call in `run_global_moves()` (line 732) with
`t_max_tip()`.  (The two uses in the `Run` constructor initializer list
cannot use the accessor since derived quantities are not yet computed.)

In `Subrun`, add `t_max_tip_` as a member grouped with `includes_run_root_`
(both describe the subrun's relationship to the overall tree), with a
getter `t_max_tip()` and setter `set_t_max_tip()`.  Push it from `Run` to
subruns in `push_global_params_to_subruns()`, before `set_evo()` to match
the field order in `subrun.h`.  Both `Spr_study` constructor calls in
`spr1_move()` pass `t_max_tip_` as the new fifth argument.

### Step 7: Build and test

Run existing test suite.  Additionally, test with the pathological dataset
`04_free_exp_pop/sims/sim_048` to confirm it no longer OOMs or overflows.

## Scope

- Modify `spr_study.cpp` (weight, sampling, density — 3 functions, plus
  updating the block comments at lines 261-299, 389-396, and 440-444).
- Modify `spr_study.h` (add `t_max_tip` member, update constructor).
- Modify `subrun.h` / `subrun.cpp` (plumb `t_max_tip` through, add
  `lambda_X == 0` guard in `spr1_move()`).
- Modify `run.cpp` / `run.h` (add `t_max_tip_` as derived quantity with
  accessor, replace explicit `calc_max_tip_time` call, push to subruns).
- Do **not** change `cell_for()` return type or `ensure_space()` — the cap
  prevents extreme values from reaching these functions.
- Do **not** change the below-root proposal — it is bounded by the tree and
  cannot overflow.
