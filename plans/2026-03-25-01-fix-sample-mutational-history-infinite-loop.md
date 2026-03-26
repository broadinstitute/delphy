# Fix Infinite Loop in `sample_mutational_history` When `mu*T` Is Tiny

## Problem

When Delphy runs on a dataset with very few or zero mutations (e.g., WCSS
`04_free_exp_pop/sims/sim_048`, which has `num_mutations = 0` and
`mu = 7.4e-7 subst/site/year`), the Gibbs sampler for `mu` drives it
extremely close to 0 (observed: `mu = 3.02e-13`).  This causes an infinite
loop in `sample_mutational_history()` (`core/spr_move.cpp`), consuming
memory without bound until the process crashes or is killed.

This bug is **independent** of the above-root SPR proposal overflow fixed
in `2026-03-24-01-cap-above-root-spr-proposal.md`.  That fix bounds how
far back the proposed root can go (capping the branch length `T`).  This
bug triggers with a perfectly normal branch length (`T = 0.000175 years`,
about 1.5 hours) — the problem is that `mu` is so small that `mu*T`
underflows to the point of exact floating-point cancellation.

## Investigation

### How the bug was found

We ran the new build (with the s_max cap fix applied) against `sim_048` in
relwithdebinfo mode.  Memory grew from ~1.2 GiB to ~6.2 GiB in
~20 seconds.  We added a temporary diagnostic to
`sample_mutational_history()` to warn when the result vector exceeds
1000 mutations:

```cpp
auto mutation_limit_warned = false;
auto check_mutation_limit = [&]() {
  if (std::ssize(result) > 1000 && !mutation_limit_warned) {
    mutation_limit_warned = true;
    std::cerr << absl::StreamFormat(
        "WARNING: sample_mutational_history generated >1000 mutations! "
        "(L=%d, T=%.17g, mu=%.17g, mu*T=%.17g, deltas=%d)\n",
        L, T, mu, mu*T, std::ssize(deltas));
  }
};
```

This diagnostic printed:

```
WARNING: sample_mutational_history generated >1000 mutations!
  (L=30000, T=0.00017547949983054423, mu=3.0238276209404132e-13,
   mu*T=5.3061975849640821e-17, deltas=0)
```

This confirms: `T` is small and normal; `mu` is absurdly small;
`mu*T ≈ 5e-17`.  We then stopped the process (SIGSTOP) and attached
gdb, finding one worker thread stuck in:

```
sample_mutational_history()          (spr_move.cpp:1295)
  <- propose_new_inner_graft_mutations()
    <- propose_new_graft()
      <- spr1_move()
```

### What gdb revealed

Key state from `info locals` in frame 9 (`sample_mutational_history`):

- **`result`**: vector of length **33,572,188** (33.5 million mutations)
  and growing
- **`mutation_limit_warned = true`**: our diagnostic fired (added during
  this session)
- **`l = -2147483647`** (`INT_MIN + 1`): the site index has overflowed
- **`log_one_minus_p_tricky = 0`**: the geometric skip rate is zero
  (see analysis below)
- **`num_muts_dist_ge1.lambda_ = 5.3e-17`**: confirms `mu*T = 5.3e-17`
- **`deltas`** has `size_ = 0`: no delta sites (consistent with 0
  mutations in the tree)
- Mutation site indices were nonsensical, far outside `[0, L)`,
  including values near `INT_MIN` — confirming the site index has
  overflowed

## Root cause analysis

### The geometric skip mechanism (spr_move.cpp lines 1268-1295)

For sites without explicit deltas (i.e., sites where the start and end
states are the same), `sample_mutational_history` uses a geometric skip
to efficiently skip over the vast majority of sites that will have empty
trajectories.  The key quantities:

```cpp
auto p_0 = std::exp(-mu*T);                              // P(0 mutations)
auto p_1 = mu * T * p_0;                                 // P(1 mutation)
auto log_one_minus_p_tricky = -mu*T - std::log1p(-p_1);  // log(p_0 / (1-p_1))
```

where `p_tricky = 1 - p_0/(1-p_1)` is the probability that rejection
sampling at a single site encounters a 2+-mutation trajectory before
accepting an empty one.  The geometric skip draws from
`Exp(-log_one_minus_p_tricky)` to determine how many sites to skip.

### The cancellation

With `mu*T = 5.3e-17`:

1. `-mu*T = -5.3e-17`
2. `p_1 = mu*T * exp(-mu*T) = 5.3e-17 * 1.0 = 5.3e-17`
3. `log1p(-p_1) = log1p(-5.3e-17) = -5.3e-17` (exact in double
   precision: `5.3e-17` is so small that `log(1+x) = x` to the last
   bit)
4. `log_one_minus_p_tricky = -5.3e-17 - (-5.3e-17) = +0.0` (exact
   cancellation to IEEE 754 positive zero, since `x - x = +0.0` for
   any finite `x`)

### The infinite loop

5. `-log_one_minus_p_tricky = -(+0.0) = -0.0` (IEEE 754 negative
   zero)
6. `std::exponential_distribution{-0.0}` is constructed.  The C++
   standard requires `lambda > 0`, so this is technically undefined
   behavior.  In practice (gcc/libstdc++), the implementation computes
   `-log(1-U) / lambda`, which gives `positive / (-0.0) = -∞`.
7. `u = -∞`
8. `u >= L` is `-∞ >= 30000` → **false** — the early exit does not fire
9. `delta = std::floor(-∞) = -∞`
10. `l += delta` — `l` is `int` and `delta` is `double`, so the `+=`
    promotes `l` to `double`, adds (giving `-∞`), and converts the
    result back to `int`.  Converting `-∞` to `int` is **undefined
    behavior**.  The compiler does not warn about this implicit
    conversion by default (`-Wfloat-conversion` would catch it, but
    it is not included in `-Wall -Wextra`).  Since the behavior is
    undefined, the resulting value of `l` is unpredictable — gdb
    showed `l = INT_MIN + 1`, but other platforms or optimization
    levels could produce different values.
11. Whatever `l` ends up as, it is far below `L` (30,000), so
    `l >= L` is **false** — the loop continues.
12. `deltas.contains(l)` is false (no site has such a negative index).
13. The loop proposes a 2+-mutation trajectory at the nonsensical site
    index `l`.  With `mu*T ≈ 5e-17`, the K-truncated Poisson almost
    always draws `n = 2`.  About 1/3 of the time the trajectory is
    accepted (end state matches start state), adding 2 mutations to
    `result` and incrementing `l`.  Otherwise, rejection sampling
    retries at the same site.
14. Back at the top of the `while(true)` loop, `u = -∞` again,
    `l` is corrupted again by the same UB, and the cycle repeats.

The loop generates approximately 2 mutations per accepted trajectory
(~1/3 acceptance rate), so ~2/3 mutations per iteration on average.
It will never terminate.

## Proposed fix

### Option A: Early return when `mu*T` is negligible

The probability of a single site encountering a 2+-mutation trajectory
that survives rejection sampling is `p_tricky ≈ (mu*T)^2 / 2` for
small `mu*T`.  The expected number of such sites across the whole
genome is `L * (mu*T)^2 / 2`.  When this is negligible (say `< 1e-6`),
the entire no-delta loop can be skipped.

The condition `L * (mu*T)^2 / 2 < 1e-6` scales naturally with genome
size:
- `L = 30,000`: triggers when `mu*T < 8.2e-6`
- `L = 1,000`: triggers when `mu*T < 4.5e-5`
- `L = 1,000,000`: triggers when `mu*T < 1.4e-6`

Skip the entire no-delta loop when `L * muT * muT < 2e-6`.  The
delta-site loop (lines 1181-1225) is unaffected — it correctly handles
tiny `mu*T` via `K_truncated_poisson_distribution` with `min_k = 1`.

### Option B: Fix the floating-point cancellation directly

Replace:
```cpp
auto log_one_minus_p_tricky = -mu*T - std::log1p(-p_1);
```
with a computation that avoids cancellation.  Since
`p_0/(1-p_1) = exp(-mu*T) / (1 - mu*T*exp(-mu*T))`, and we want
`log(p_0/(1-p_1))`, we can write:
```cpp
auto log_one_minus_p_tricky = -mu*T - std::log1p(-mu*T * std::exp(-mu*T));
```
This is what the code already does.  The cancellation occurs because
for tiny `x = mu*T`, both terms evaluate to approximately `-x`:
- `-mu*T = -x`
- `log1p(-x * exp(-x)) ≈ log1p(-x) ≈ -x`

Their difference is `O(x^2)`, which is lost to floating-point rounding
when `x` is small enough that `-x` and `log1p(-x * exp(-x))` round to
the same `double`.  In IEEE 754 double precision (52-bit mantissa),
the `O(x^2)` correction is lost when `x^2 / x = x < 2^{-52} ≈ 2e-16`,
which is the machine epsilon.  In practice, cancellation starts causing
significant relative error around `x ≈ 1e-8` (where the `x^2/2` term
has only ~8 significant digits left).

A Taylor expansion gives:
```
log(p_0/(1-p_1)) = -mu*T - log(1 - mu*T*exp(-mu*T))
                  ≈ -(mu*T)^2/2   for small mu*T
```

So a numerically stable version is:
```cpp
auto muT = mu * T;
auto log_one_minus_p_tricky = (muT < 1e-4)
    ? -0.5 * muT * muT   // Taylor: -(mu*T)^2/2 + O((mu*T)^3)
    : -muT - std::log1p(-muT * std::exp(-muT));
```

### Recommendation

Use **both** fixes:
- **Option B** (Taylor expansion) as the primary fix: it makes the
  geometric skip numerically correct for all `mu*T` values, eliminating
  the cancellation.
- **Option A** (skip loop) as a defense-in-depth guard: even with
  the Taylor fix, the geometric skip with tiny `mu*T` will immediately
  exit the loop (since `u` will be enormous), so it's effectively a
  no-op.  But an explicit skip makes the intent clear and protects
  against any other unforeseen numerical issues.

Additionally, the geometric skip loop itself is hardened with defensive
checks (Step 5): the `u >= L` exit condition is widened to catch NaN,
negative, and infinite values; the `double`-to-`int` conversion is
made explicit; and `CHECK`s verify that site indices stay in range.

## Implementation

### Step 1: Add a regression test (spr_move_tests.cpp)

Add a new test case that calls `sample_mutational_history` with a tiny
`mu*T` that reliably triggers the bug.  The test should pass after the
fix is applied and would hang (or crash) without it.

In `tests/spr_move_tests.cpp`, add a new test after the existing
`sample_mutational_history` test:

```cpp
TEST_F(Spr_move_simple_test, sample_mutational_history_tiny_muT) {
  // Regression test for infinite loop when mu*T is tiny (~5e-17).
  // Observed in WCSS 04_free_exp_pop/sims/sim_048 of delphy-2026-paper-data.
  // See plans/2026-03-25-01-fix-sample-mutational-history-infinite-loop.md for details.
  auto L = 30'000;  // Typical genome size
  auto T = 0.000175;  // ~1.5 hours, perfectly normal branch length
  auto mu = 3e-13;  // Absurdly small mu, as observed in sim_048
  // mu*T ≈ 5.3e-17

  auto deltas = Site_deltas{};  // No deltas, as in the observed crash

  auto bitgen = std::mt19937_64(12345);
  auto result = sample_mutational_history(L, T, mu, deltas, bitgen);

  // With mu*T ≈ 5e-17 and L = 30,000, the expected number of sites with
  // 2+-mutation trajectories is L*(mu*T)^2/2 ≈ 4e-29.  So we expect
  // zero mutations in the result.
  EXPECT_THAT(std::ssize(result), testing::Eq(0));
}
```

This test does not need the full tree fixture — only `L`, `T`, `mu`,
and an empty `deltas`.  Without the fix, it will hang in the geometric
skip loop; with the fix, it should return immediately with an empty
result.

Also add a companion test for `sample_unconstrained_mutational_history`
to confirm it handles tiny `mu*T` correctly (it uses the Gillespie
algorithm, so it should exit immediately — but let's nail that down
with a test):

```cpp
TEST_F(Spr_move_simple_test, sample_unconstrained_mutational_history_tiny_muT) {
  // Confirm that sample_unconstrained_mutational_history handles tiny mu*T
  // without hanging.  It uses the Gillespie algorithm (Exp(mu*L) inter-event
  // times), so when mu is tiny the first draw overshoots T and the loop exits
  // immediately.
  // Same parameters as the sim_048 scenario above.
  // See plans/2026-03-25-01-fix-sample-mutational-history-infinite-loop.md for details.
  auto L = 30'000;
  auto T = 0.000175;
  auto mu = 3e-13;

  auto bitgen = std::mt19937_64(12345);
  auto result = sample_unconstrained_mutational_history(L, T, mu, bitgen);

  EXPECT_THAT(std::ssize(result), testing::Eq(0));
}
```

Build and run the tests at this point.  The first test
(`sample_mutational_history_tiny_muT`) should hang or crash,
reproducing the bug.  The second test
(`sample_unconstrained_mutational_history_tiny_muT`) should pass
immediately, confirming that the Gillespie algorithm is unaffected.
Once both behaviors are confirmed, proceed with the fix steps below.

### Step 2: Fix the floating-point cancellation (spr_move.cpp ~line 1273)

Replace:
```cpp
//auto log_one_minus_p_tricky = std::log(p_0) - std::log1p(-p_1);
auto log_one_minus_p_tricky   =     -mu*T     - std::log1p(-p_1);
```

with:
```cpp
auto log_one_minus_p_tricky = (muT < 1e-4)
    ? -0.5 * muT * muT
    : -muT - std::log1p(-p_1);
```

The non-Taylor branch reuses `p_1` (which remains declared above,
along with `p_0`); it is identical to the original formula.  The Taylor branch
avoids the cancellation entirely.

Keep the declarations of `p_0` and `p_1` (lines 1279-1280) and the
commented-out alternatives for `p_tricky` and `one_minus_p_tricky`
(lines 1281-1282) — they document the derivation and `p_1` is still
used in the non-Taylor branch.  The commented-out
`log_one_minus_p_tricky` (line 1283) is replaced by the new code.

Note: `muT` is defined once in Step 4.  Although Step 4 appears
later in this plan, its code (`auto muT = mu * T;`) is placed
earlier in the function (before `p_0` and `p_1`).  Both the skip
guard (Step 4) and this Taylor/original branch reuse the same
`muT` variable.

The threshold `1e-4` is chosen so that the relative error from
truncating the Taylor series at the quadratic term is at most
`|(mu*T)^3/3| / |(mu*T)^2/2| = 2*mu*T/3 < 7e-5`, i.e., the
approximation is accurate to better than 1 part in 10,000.  Above
`1e-4`, the original formula is accurate (the two terms being
subtracted differ enough to avoid significant cancellation).

### Step 3: Verify the Taylor expansion

For `x = mu*T`:
```
log(p_0/(1-p_1)) = log(e^{-x} / (1 - x e^{-x}))
                  = -x - log(1 - x e^{-x})
```

Taylor expand `x e^{-x} = x - x^2 + x^3/2 - ...` for small `x`:
```
1 - x e^{-x} = 1 - x + x^2 - x^3/2 + ...
log(1 - x e^{-x}) = log(1 - x + x^2 - ...) ≈ -x + x^2/2 + ...
```
So:
```
-x - log(1 - x e^{-x}) ≈ -x - (-x + x^2/2) = -x^2/2
```

This confirms `-0.5 * muT * muT` is correct to leading order.

### Step 4: Skip no-delta loop when `mu*T` is negligible (defense-in-depth)

At the start of the no-delta loop section (~line 1268), before the
declarations of `p_0` and `p_1`, add `muT`.  Then change the
initialization of `l` and the `while` condition:

```cpp
auto muT = mu * T;
```

This `muT` variable is also reused by the Taylor expansion fix in
Step 2.

Also update the existing `p_0` and `p_1` declarations to use `muT`
instead of recomputing `mu*T`:
```cpp
auto p_0 = std::exp(-muT);
auto p_1 = muT * p_0;
```

After the declarations of `p_0`, `p_1`, and `log_one_minus_p_tricky`,
replace (~lines 1286-1287):
```cpp
auto l = 0;
while (true) {
```
with:
```cpp
auto l = 0;
if (L * muT * muT < 2e-6) {
  // The expected number of no-delta sites with 2+-mutation trajectories
  // is L*(mu*T)^2/2.  When this is negligible, set l to L so the
  // while (l < L) loop below is skipped entirely.
  l = L;
}
while (l < L) {
```

The loop now terminates naturally when `l` advances past `L`, whether
that happens via the geometric skip, via `++l` after accepting a
trajectory, or immediately because `l` was set to `L` by the guard
above.  All code paths still reach the single `sort_mutations(result)`
call at the bottom of the function.

This threshold scales with genome size `L`: for `L = 30,000` it
skips when `mu*T < 8.2e-6`; for `L = 1,000,000` when
`mu*T < 1.4e-6`.

### Step 5: Harden the geometric skip loop against numerical issues

Even with the Taylor fix (Step 2) and skip guard (Step 4),
add defensive checks to the geometric skip loop so that any future
numerical issues cannot cause undefined behavior or an infinite loop.

**5a.** After drawing `u`, break if it is not a finite positive number
(~line 1295):

Replace:
```cpp
auto u = std::exponential_distribution{-log_one_minus_p_tricky}(bitgen);
if (u >= L) {
  break;  // when we turn u into a Geo(p_tricky) sample and add it to l, we'd have l >= L
}
```
with:
```cpp
auto u = std::exponential_distribution{-log_one_minus_p_tricky}(bitgen);
if (!(u >= 0 && u < L)) {
  // Normally, u >= L means we skip past all remaining sites.
  // But u can also be negative, infinite, or NaN if log_one_minus_p_tricky
  // is corrupted by floating-point cancellation (e.g., cancellation to 0
  // gives exponential_distribution{-0.0}, which returns -inf).
  // In all such cases, we should exit the loop.
  break;
}
```

This subsumes the old `u >= L` check and additionally catches negative,
infinite, and NaN values of `u`, which would arise from numerical
issues in `log_one_minus_p_tricky`.

**5b.** Add `static_cast<int>` to the `l += delta` line (~line 1300)
to document that the conversion from `double` to `int` is safe:

Replace:
```cpp
auto delta = std::floor(u);
l += delta;
```
with:
```cpp
l += static_cast<int>(std::floor(u));  // safe: u is in [0, L), so floor(u) is in [0, L-1]
```

Since `u` is guaranteed to be in `[0, L)` by the check in 5a (and `L`
is a `Site_index`, i.e., `int`), `std::floor(u)` is in `[0, L-1]` and
fits in `int`.  The sum `l + floor(u)` can be up to `2L-2`, which also
fits in `int`.  The existing `if (l >= L) { break; }` check remains.

**5c.** Add a `CHECK` before any mutation is pushed in the no-delta loop
to verify `l` is in range (~line 1339):

```cpp
CHECK_LE(0, l);
CHECK_LT(l, L);
```

Using separate `CHECK_LE` / `CHECK_LT` ensures failure messages show
the values of both sides.

Also add an equivalent `CHECK` in the delta-site loop (earlier in the
function) before mutations are pushed for delta sites.  There, `l`
comes from iterating over the `deltas` map, so this validates the
caller's input rather than internal computation.

### Step 6: Update the block comment (spr_move.cpp ~lines 1242-1277)

Update the long block comment preceding the no-delta loop to document:
- The floating-point cancellation problem when `mu*T` is tiny
- The Taylor expansion fix for `log_one_minus_p_tricky`
- The skip guard when `L * muT * muT < 2e-6`

### Step 7: Remove diagnostic instrumentation (spr_move.cpp)

Remove the temporary `check_mutation_limit` lambda and its call sites
that were added during the investigation of this bug.

### Step 8: Build and test

Run the existing test suite, including the new regression test from
Step 1.  The new test should pass (returning zero mutations).
Additionally, test with the pathological dataset
`04_free_exp_pop/sims/sim_048` to confirm it no longer hangs.

### Step 9: Update test comments (spr_move_tests.cpp)

After confirming the fix works, update the comment on
`sample_mutational_history_tiny_muT` to note that the fix is in place
and to document the pre-fix behavior (exhausts memory, throws
`std::bad_alloc` after ~4 minutes):

```cpp
  // Regression test for infinite loop when mu*T is tiny (~5e-17), fixed by Taylor expansion
  // for log_one_minus_p_tricky and a skip guard.  Without the fix, this test exhausts memory
  // and throws std::bad_alloc after ~4 minutes.
  // Observed in WCSS 04_free_exp_pop/sims/sim_048 of delphy-2026-paper-data.
  // See plans/2026-03-25-01-fix-sample-mutational-history-infinite-loop.md for details.
```

## Scope

- Modify `spr_move.cpp` (Taylor expansion fix, skip guard,
  defensive loop hardening, block comment update, diagnostic removal).
- Add tests in `spr_move_tests.cpp` (regression test for tiny `mu*T`,
  confirmation test for `sample_unconstrained_mutational_history`).
- Do **not** change the delta-site loop — it handles tiny `mu*T`
  correctly.  (A site-range `CHECK` is added there, but no logic
  changes.)
- Do **not** change `sample_unconstrained_mutational_history` — it uses
  the Gillespie algorithm (not the geometric skip), so it is unaffected
  by this bug.  A test is added in Step 1 to confirm.
