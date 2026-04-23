# Fix `gamma_q` / `gamma_q_inv` nan on Emscripten/WASM

## Problem

When Delphy is compiled to WebAssembly (via Emscripten) and runs on datasets
with deeply divergent lineages, the MCMC can reach states where calls to
Boost's `gamma_q` and `gamma_q_inv` produce `nan`, crashing the run:

```
spr_study.cpp:484 Check failed: rand_t >= region.t_min - 1e-6 (nan vs -1.79769e+308)
```

### Background: `gamma_q`

`gamma_q(a, x)` is the upper regularized incomplete gamma function: the
probability that a `Gamma(a, 1)` random variable (shape `a`, rate 1)
exceeds `x`.  Its behavior for fixed `a`:

- `gamma_q(a, 0) = 1`
- `gamma_q(a, x) ≈ 1` for `x << a` (almost all mass is above `x`)
- `gamma_q(a, a) ≈ 0.5` (the median of `Gamma(a, 1)` is near `a`)
- `gamma_q(a, x) → 0` for `x >> a` (negligible upper tail)

The transition from ~1 to ~0 happens in the region around `x = a`.  Concretely, since the
standard deviation of `Gamma(a, 1)` is `sqrt(a)`, `gamma_q` goes from ~0.98 at `x = a -
2*sqrt(a)` to ~0.02 at `x = a + 2*sqrt(a)`.  This approximate behaviour is more accurate
for higher values of `a` (around `a >= ~20`).

### Where `gamma_q` is used in Delphy

The above-root mutation-directed SPR (mdSPR) proposal samples a candidate
reattachment time from a truncated Gamma distribution.  The tree topology determines a set of
candidate reattachment regions, one of which is above the root, and each
requiring at least `min_muts` mutations (`m` below).  The proposal assigns
weight to each region proportional to the Poisson likelihood of observing
`m` mutations in the corresponding evolutionary distance, modified by an
empirically effective annealing exponent `f = 0.8`.  Specifically, if `s`
is the sum of the two branch lengths from the proposed new root to its
children X and S after regrafting, and `lambda_X` is the cumulative
mutation rate at X (rate per site times number of non-missing sites), then

```
p_root(s) ∝ (lambda_X * s)^{fm} * exp(-lambda_X * f * s)
```

which is a `Gamma(fm + 1, lambda_X * f)` density.  The variable `s` ranges
over `[s_min, s_max]`, where `s_min = |t_X - t_S|` (the difference in
node times of X and S) is the minimum evolutionary distance (the
children must straddle the new root), and
`s_max` is capped at `s_min + 20 * tree_span` to avoid numerical overflow
(see `plans/2026-03-24-01-cap-above-root-spr-proposal.md`).

After substituting `x = lambda_X * f * s`, the weight integral becomes:

```
W_R ∝ gamma_q(fm+1, x_min) - gamma_q(fm+1, x_max)
```

A specific attachment point within this region is sampled from the
corresponding truncated Gamma distribution.  Three operations depend on
`gamma_q` and `gamma_q_inv`:

1. **Weight:** `log(gamma_q(a, x_min) - gamma_q(a, x_max))` — the log
   probability of a `Gamma(a, 1)` variate falling in `[x_min, x_max]`
2. **Sampling:** inverse-CDF via `gamma_q_inv` — draw `Q` uniformly in
   `(gamma_q(a, x_max), gamma_q(a, x_min))`, then invert
3. **Proposal density (`log_alpha`):** same as weight but evaluated at a
   specific point

A separate call site in `run.cpp` (the mpox rho sampler, line 935) performs
a similar truncated Gamma sample with an existing ad-hoc guard.

### Why `gamma_q` produces nan on WASM

By default, Boost promotes `double` arguments to `long double` before
computing `gamma_q`.  It then selects the Lanczos approximation matching
the precision of `long double` and relies on math functions like `expl` and
`powl` to work correctly across the full exponent range of that type.

On **x86-64**, `long double` is 80-bit extended precision
(max_exponent = 16384, log_max ≈ 11356).  Boost selects `lanczos17m64`
(64-bit precision).  Boost correctly implements underflow and overflow
guards to calculate `gamma_q` and `gamma_q_inv` for all parameter
values, provided that the underlying math library functions (`expl`,
`powl`, etc.) are as accurate as documented — which they are on x86-64.

On **Emscripten/WASM**, `long double` is declared as 128-bit IEEE binary128
(max_exponent = 16384), and `numeric_limits` reports these properties.
Boost selects `lanczos24m113` (113-bit precision) and sets its overflow
guards to `log_max_value = 11356`.  However, Emscripten's `long double`
math functions (`expl`, `powl`, etc.) are thin wrappers around the
`double` versions — they overflow at the `double` threshold (~710), far
below the declared 11356.  When Boost's overflow guards let through an
intermediate value that exceeds the `double` range, `pow` overflows to
`inf`, `exp` underflows to 0, and the product becomes `inf * 0 = nan`.

Empirically confirmed (x86-64 Linux, em++ 3.1.69, node v18.20.3,
g++ 13.3.0):

```
# WASM (via em++ / node):
expl(709)   = 8.22e+307   ← OK (same as double)
expl(710)   = inf         ← should be 2.23e+308 for true 128-bit
expl(11355) = inf         ← should be 2.59e+4931

# x86-64 (g++):
expl(709)   = 8.22e+307
expl(710)   = 2.23e+308   ← OK
expl(11355) = 2.59e+4931  ← OK
```

`gamma_q_inv` is affected indirectly: it calls `gamma_q` internally during
its root-finding iterations, inheriting the same promotion and overflow
behavior.  This causes not only nans but also **silently wrong results**
on WASM (e.g., `gamma_q_inv(500, 1e-300)` returns 1245 instead of the
correct 1836).

### The fix: suppress `double` → `long double` promotion

Calling `gamma_q` and `gamma_q_inv` with `promote_double<false>` policy
keeps the entire computation in `double`.  Boost then selects
`lanczos13m53` (53-bit precision) and sets its overflow guards to
`log_max_value = 709`, which matches the actual capabilities of the math
functions on all platforms.

With `promote_double<false>`, all previously-nan cases return the correct
result on both x86-64 and WASM:

```
gamma_q(271.4, 6601, no_promote{}) = 0      ← was nan with default policy on WASM
gamma_q(1000,  100,  no_promote{}) = 1      ← was nan with default policy on WASM
```

Comprehensive testing across 900+ parameter combinations (see Appendix below):

- **x86-64:** zero mismatches above 1e-12 relative error between default
  and no-promote policies (max relative error ~1.3e-13).
- **WASM:** 43 `gamma_q` nan cases and 14 `gamma_q_inv` nan/mismatch
  cases are all fixed by no-promote; the no-promote results agree with the
  x86-64 results.

The `long double` promotion provides negligible precision benefit for
`gamma_q` / `gamma_q_inv`.


## Proposed fix

### Commit 1: Minimal fix — add `promote_double<false>` to all call sites

Define a no-promote policy alias in each file that calls `gamma_q` or
`gamma_q_inv`:

```cpp
using gamma_policy = boost::math::policies::policy<
    boost::math::policies::promote_double<false>>;
```

Pass `gamma_policy{}` as the third argument to every `gamma_q` and
`gamma_q_inv` call:

**spr_study.cpp — Weight (line 369):**
```cpp
std::log(boost::math::gamma_q(f*m + 1, x_min, gamma_policy{})
       - boost::math::gamma_q(f*m + 1, x_max, gamma_policy{}));
```

**spr_study.cpp — Sampling (line 465):**
```cpp
auto Q_min = boost::math::gamma_q(f*m + 1, x_max, gamma_policy{});
auto Q_max = boost::math::gamma_q(f*m + 1, x_min, gamma_policy{});
...
rand_s = boost::math::gamma_q_inv(f*m + 1, rand_Q, gamma_policy{}) / (lambda_X*f);
```

**spr_study.cpp — log_alpha_in_region (line 622):**
```cpp
-std::log(boost::math::gamma_q(f*m+1, x_min, gamma_policy{})
        - boost::math::gamma_q(f*m+1, x_max, gamma_policy{}));
```

**run.cpp — mpox rho sampler (line 935):**
```cpp
// Before:
auto max_Q = (z < 0.1*a ? 1.0 : boost::math::gamma_q(a, z));
...
auto k = boost::math::gamma_q_inv(M_star + 1, rand_Q) / ...;

// After:
auto max_Q = boost::math::gamma_q(a, z, gamma_policy{});
...
auto k = boost::math::gamma_q_inv(M_star + 1, rand_Q, gamma_policy{}) / ...;
```

This removes the ad-hoc `z < 0.1*a ? 1.0 :` guard in `run.cpp`, which was
a workaround for the same underlying issue.

Also add defensive CHECKs at key intermediate points:

1. **After `gamma_q` calls:** `CHECK(not std::isnan(Q))` for each result,
   and `CHECK_GE(Q_max - Q_min, 0.0)`.

2. **Before `Uniform` calls:** `CHECK_LT(Q_min, Q_max)` before every
   `Uniform(OpenOpen, Q_min, Q_max)`.

3. **After `gamma_q_inv`:** `CHECK(not std::isnan(rand_s))` and
   `CHECK(not std::isinf(rand_s))`.

4. **After computing `rand_t`:** `CHECK(not std::isnan(rand_t))` (already
   present; keep it).

5. **After `log_alpha_in_region`:** `CHECK(not std::isnan(result))` at the
   end of the function.


### Commit 2: Refactor into `safe_gamma_math.h`

Create `core/safe_gamma_math.h` with wrapper functions that hide the
`promote_double<false>` policy and make the correct usage discoverable for
future callers.  This header should include an extended comment explaining:
- the problem with Emscripten's `long double` math library functions
- how Boost's `gamma_q` / `gamma_q_inv` interact badly with this
- why `promote_double<false>` avoids the issue

The wrapper functions:

**`safe_gamma_q(a, x)`:** calls `boost::math::gamma_q(a, x, gamma_policy{})`
with a post-call `CHECK` for nan and inf.

**`safe_gamma_q_inv(a, Q)`:** validates `Q` is in `[0, 1]`, calls
`boost::math::gamma_q_inv(a, Q, gamma_policy{})` inside a try/catch that
logs `a` and `Q` if Boost throws, and post-call `CHECK`s for nan and inf.

**`safe_log_gamma_integral(a, x_min, x_max)`:** computes
`log(gamma_q(a, x_min) - gamma_q(a, x_max))` via `safe_gamma_q`, with a
`CHECK` that `Q_hi >= Q_lo`.

**`safe_sample_truncated_gamma(alpha, beta, lo, hi, bitgen)`:**
samples from `Gamma(alpha, beta)` truncated to `[lo, hi]` (where `hi`
may be infinity) using inverse-CDF.  Transforms to unit-rate `Gamma(alpha, 1)` internally
(`y = beta * x`), calls `safe_gamma_q` at both endpoints, samples `Q`
uniformly in the resulting interval via `Uniform(OpenClosed, Q_lo, Q_hi)`,
inverts via `safe_gamma_q_inv`, and transforms back.  The `OpenClosed`
bounds avoid `Q = 0`, which would make `gamma_q_inv` return `inf`.  At the
upper end, `Q = Q_hi` is safe since `gamma_q_inv(a, Q_hi)` returns a
finite value (the lower truncation bound `y_lo`).  This is an intentional safety improvement over the `OpenOpen` interval
used in the existing code, which could sample `Q` arbitrarily close to 0.

Includes `CHECK`s that:
- The Q range is non-degenerate (`Q_lo < Q_hi`)
- The sampled value is finite
- The result is in `[lo, hi]`

Refactor all four call sites (3 in `spr_study.cpp`, 1 in `run.cpp`) to use
these wrappers.  The mpox rho sampler in `run.cpp` should use
`safe_sample_truncated_gamma` directly rather than manually computing Q
values and calling `gamma_q_inv`.  This safely replaces the ad-hoc
`z < 0.1*a ? 1.0 :` guard, which is no longer needed because
`safe_gamma_q` never underflows or overflows.

Create `tests/safe_gamma_math_tests.cpp` with tests exercising extreme
parameter regimes:

1. **`safe_gamma_q`:**
   - Small `a`, normal `x` (e.g., `a=2, x=5`): matches `gamma_q` exactly
   - Large `a`, `x >> a` (e.g., `a=271.4, x=6601`): returns 0
   - Large `a`, `x << a` (e.g., `a=1000, x=100`): returns 1
   - Near-mode `x ≈ a` (e.g., `a=271.4, x=280`): matches `gamma_q`
   - Large `a`, moderate `x` (e.g., `a=500, x=1000`): finite result

2. **`safe_gamma_q_inv`:**
   - Round-trips: `safe_gamma_q_inv(a, safe_gamma_q(a, x)) ≈ x` for
     various `a` and `x` in the non-degenerate regime
   - Extreme Q values: `safe_gamma_q_inv(a, 1e-300)` returns a large
     but finite value

3. **`safe_log_gamma_integral`:**
   - `(a=5, x_min=1, x_max=10)`: finite, ≈ log(1) = 0 (most mass
     captured)
   - `(a=271.4, x_min=148, x_max=6601)`: finite, ≈ 0
   - `(a=5, x_min=1000, x_max=2000)`: `-inf` (both Q ≈ 0)

4. **`safe_sample_truncated_gamma`:**
   - Normal parameters: samples in `[lo, hi]`
   - Large `alpha`: samples in `[lo, hi]`
   - Various `alpha` and `beta` combinations

5. **Comprehensive sweep:** the 900+ parameter combinations used for
   validation (varying `a` across {0.01, ..., 1000}, `x` across many
   orders of magnitude, and `Q` across {1e-300, ..., 1-1e-15}; see Appendix below).  Verify
   `safe_gamma_q` and `safe_gamma_q_inv` return finite, non-nan results
   for all combinations.


## Correctness argument

The `promote_double<false>` policy changes only the intermediate precision
of the Boost computation, not the mathematical result.  Boost's `double`
code path (`lanczos13m53`) has well-tested overflow guards that handle all
parameter ranges correctly.

No new code paths, fallback sampling strategies, or approximations are
introduced.  The proposal density computed by `log_alpha_in_region` still
exactly matches the sampling distribution used by `pick_time_in_region`
(both functions in `spr_study.cpp`), so the Metropolis-Hastings ratio
remains correct.


## Appendix: Validation test program

The following standalone program was used to compare default vs
`promote_double<false>` across 900+ parameter combinations.  It was
compiled and run on both x86-64 (g++ 13.3.0) and WASM (em++ 3.1.69 /
node v18.20.3):

```
# x86-64:
g++ -std=c++20 -O2 -I<boost-include-path> test_gamma_q_no_promote.cpp -o test && ./test

# WASM:
em++ -std=c++20 -O2 -fno-exceptions -I<boost-include-path> test_gamma_q_no_promote.cpp -o test.js && node test.js
```

```cpp
#include <cstdio>
#include <cmath>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/policies/policy.hpp>

using no_promote = boost::math::policies::policy<
    boost::math::policies::promote_double<false>>;

int main() {
  int mismatches = 0;
  int nan_cases = 0;
  int total = 0;
  double max_reldiff = 0;
  double worst_a = 0, worst_x = 0;

  auto test_q = [&](double a, double x) {
    ++total;
    auto r_default = boost::math::gamma_q(a, x);
    auto r_nopro = boost::math::gamma_q(a, x, no_promote());

    bool nan_d = std::isnan(r_default);
    bool nan_n = std::isnan(r_nopro);

    if (nan_d || nan_n) {
      if (nan_d != nan_n) {
        std::printf("NAN MISMATCH: gamma_q(%g, %g) default=%g no_promote=%g\n",
                    a, x, r_default, r_nopro);
        ++nan_cases;
      }
      return;
    }

    if (r_default != r_nopro) {
      double reldiff = 0;
      if (r_default != 0)
        reldiff = std::fabs((r_nopro - r_default) / r_default);
      else if (r_nopro != 0)
        reldiff = std::fabs(r_nopro);

      if (reldiff > max_reldiff) {
        max_reldiff = reldiff;
        worst_a = a;
        worst_x = x;
      }

      if (reldiff > 1e-12) {
        std::printf("MISMATCH: gamma_q(%g, %g) default=%.17g no_promote=%.17g"
                    " reldiff=%.3e\n", a, x, r_default, r_nopro, reldiff);
        ++mismatches;
      }
    }
  };

  auto test_q_inv = [&](double a, double Q) {
    if (Q <= 0 || Q >= 1) return;
    ++total;
    auto r_default = boost::math::gamma_q_inv(a, Q);
    auto r_nopro = boost::math::gamma_q_inv(a, Q, no_promote());

    bool nan_d = std::isnan(r_default);
    bool nan_n = std::isnan(r_nopro);

    if (nan_d || nan_n) {
      if (nan_d != nan_n) {
        std::printf("NAN MISMATCH: gamma_q_inv(%g, %g) default=%g no_promote=%g\n",
                    a, Q, r_default, r_nopro);
        ++nan_cases;
      }
      return;
    }

    if (r_default != r_nopro) {
      double reldiff = 0;
      if (r_default != 0)
        reldiff = std::fabs((r_nopro - r_default) / r_default);
      else if (r_nopro != 0)
        reldiff = std::fabs(r_nopro);

      if (reldiff > max_reldiff) {
        max_reldiff = reldiff;
        worst_a = a;
        worst_x = Q;
      }

      if (reldiff > 1e-12) {
        std::printf("MISMATCH: gamma_q_inv(%g, %g) default=%.17g no_promote=%.17g"
                    " reldiff=%.3e\n", a, Q, r_default, r_nopro, reldiff);
        ++mismatches;
      }
    }
  };

  // ---- gamma_q tests ----
  std::printf("=== gamma_q: varying a and x ===\n");

  // Small a
  for (double a : {0.01, 0.1, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0}) {
    for (double x = 0.001; x <= 1000; x *= 1.7) {
      test_q(a, x);
    }
  }

  // Medium a
  for (double a : {20.0, 50.0, 100.0, 150.0, 200.0}) {
    for (double x = 0.01; x <= 10000; x *= 1.5) {
      test_q(a, x);
    }
  }

  // Large a (the regime Delphy hits)
  for (double a : {250.0, 271.4, 300.0, 400.0, 500.0, 750.0, 1000.0}) {
    // Near mode (x ~ a)
    for (double ratio : {0.01, 0.1, 0.5, 0.8, 0.9, 0.95, 0.99, 1.0,
                         1.01, 1.05, 1.1, 1.2, 1.5, 2.0, 5.0, 10.0,
                         20.0, 50.0}) {
      test_q(a, a * ratio);
    }
    // Very large x
    for (double x : {1000.0, 5000.0, 10000.0, 50000.0, 100000.0}) {
      test_q(a, x);
    }
    // Very small x
    for (double x : {0.001, 0.01, 0.1, 1.0, 10.0}) {
      test_q(a, x);
    }
  }

  // Extreme edge: x=0
  for (double a : {1.0, 10.0, 100.0, 1000.0}) {
    test_q(a, 0.0);
  }

  // Near-integer and near-half-integer a
  for (double a : {10.0, 10.5, 100.0, 100.5}) {
    for (double x = 0.1; x <= 1000; x *= 2) {
      test_q(a, x);
    }
  }

  std::printf("gamma_q: %d total, %d mismatches (>1e-12), %d nan mismatches,"
              " max reldiff=%.3e at a=%g x=%g\n\n",
              total, mismatches, nan_cases, max_reldiff, worst_a, worst_x);

  // ---- gamma_q_inv tests ----
  total = mismatches = nan_cases = 0;
  max_reldiff = 0;

  std::printf("=== gamma_q_inv: varying a and Q ===\n");

  for (double a : {0.5, 1.0, 2.0, 5.0, 10.0, 50.0, 100.0, 200.0,
                   271.4, 300.0, 500.0, 1000.0}) {
    for (double Q : {1e-300, 1e-200, 1e-100, 1e-50, 1e-20, 1e-10, 1e-5,
                     0.001, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99,
                     0.999, 0.9999, 1.0 - 1e-10, 1.0 - 1e-15}) {
      test_q_inv(a, Q);
    }
  }

  std::printf("gamma_q_inv: %d total, %d mismatches (>1e-12), %d nan mismatches,"
              " max reldiff=%.3e at a=%g Q=%g\n",
              total, mismatches, nan_cases, max_reldiff, worst_a, worst_x);

  return 0;
}
```

### Results

**x86-64 (g++ 13.3.0):**
```
gamma_q:     674 total, 0 mismatches, 0 nan mismatches, max reldiff=1.259e-13
gamma_q_inv: 228 total, 0 mismatches, 0 nan mismatches, max reldiff=4.175e-15
```

**WASM (em++ 3.1.69, node v18.20.3):**
```
gamma_q:     674 total, 2 mismatches, 43 nan mismatches, max reldiff=6.418e-08
gamma_q_inv: 228 total, 12 mismatches, 2 nan mismatches, max reldiff=1.000e+00
```

All nan mismatches are cases where the default policy produces nan but
`no_promote` returns the correct finite value.  The `gamma_q_inv`
mismatches on WASM are cases where the default policy returns silently
wrong values (e.g., `gamma_q_inv(500, 1e-300)` returns 1245 instead of the
correct 1836) due to the broken `long double` root-finding.
