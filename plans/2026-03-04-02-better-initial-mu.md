# Better Initial Mu When Prior Is Specified

## Motivation

When the user specifies a Gamma prior on mu (via `--v0-mu-prior-mean/stddev` or
`--v0-mu-prior-alpha/beta`), the initial value of mu defaults to `1e-3 subst/site/year`
regardless.  If the prior mean is far from `1e-3`, the sampler wastes burn-in time
moving toward the prior's typical region.

A better default: when a proper prior is specified and the user has not explicitly set
`--v0-init-mutation-rate`, use the prior mean (`alpha / beta`) as the initial value.

## When Is the Prior Proper?

The Gamma(alpha, beta) prior has a finite mean iff `beta > 0` (and `alpha > 0`).

- `alpha = 1, beta = 0` (default): improper uniform prior — no finite mean; keep the
  current default of `1e-3`.
- `alpha > 0, beta > 0`: proper Gamma — mean = `alpha / beta` (in per-day units
  internally).
- `alpha = 0, beta > 0`: improper Jeffreys prior (`p(mu) ~ 1/mu`) — the mean is
  undefined (the distribution is improper); keep the current default of `1e-3`.

Decision rule: if `mu_prior_alpha > 0` and `mu_prior_beta > 0` and the user did not
explicitly provide `--v0-init-mutation-rate`, set `init_mu = mu_prior_alpha / mu_prior_beta`
(already in per-day units at this point in the code).

## Precedent

The skygrid tau prior already does this (see `cmdline.cpp` line 590):
```cpp
skygrid_tau = prior_alpha / prior_beta;  // Prior mean as initial value
```

## Implementation Plan

### Step 1: `core/cmdline.cpp` — Use prior mean as default initial mu

Currently, `init_mu` is read unconditionally from `--v0-init-mutation-rate` (default `1e-3`)
and converted to per-day units.  The prior parameters are computed afterward.

**Change:** After the prior parameters are resolved, check whether the user explicitly
provided `--v0-init-mutation-rate`.  If not, and the prior is proper, override `init_mu`
with the prior mean:

```cpp
// If the user didn't explicitly set an initial mu and the prior is proper,
// use the prior mean as the initial value
if (opts.count("v0-init-mutation-rate") == 0 && mu_prior_alpha > 0.0 && mu_prior_beta > 0.0) {
  init_mu = mu_prior_alpha / mu_prior_beta;  // already in per-day units
}
```

Note: `opts.count("v0-init-mutation-rate")` returns 0 when the user did not pass the flag
on the command line, even when the option has a `default_value`.  This is because cxxopts's
`parse_default` does not increment `m_count` (unlike the regular `parse` path which does
`++m_count`).  Verified in `third-party/cxxopts/include/cxxopts.hpp` lines 1421-1434.

### Step 2: Update `--v0-init-mutation-rate` help text

Update the help string to mention the new behavior:

```cpp
("v0-init-mutation-rate",
 "Initial (or fixed) value of mutation rate, in subst / site / year.  "
 "If not specified and a proper Gamma prior is set, defaults to the prior mean; "
 "otherwise defaults to 1e-3.",
 cxxopts::value<double>()->default_value("1e-3"))
```

### Step 3: No other files need changes

This change is entirely within `cmdline.cpp` parsing logic.  No changes are needed to:

- **`run.h` / `run.cpp`**: `set_mu()` is called with whatever value `init_mu` holds;
  no new members or logic needed.
- **`api.fbs` / `api.cpp`**: the initial mu is not persisted separately; the current mu
  is stored via `run.mu()` at save time.
- **`delphy_wasm.cpp`**: the web UI sets mu directly; it doesn't go through CLI parsing.
- **`beasty_input.cpp`**: exports the current mu, not the initial value.

### Step 4: Build and test

Build Release and run with `--v0-steps 1 --v0-log-every 1 --v0-tree-every 1
--v0-delphy-snapshot-every 1` to verify the initial mu from the Step 0 log line.
(`--v0-steps 0` is rejected by the CLI.)

1. No prior flags — should show `mu = 1 * 10^-3` subst/site/year.
2. `--v0-mu-prior-mean 5e-4 --v0-mu-prior-stddev 2e-4` — should show `mu = 0.5 * 10^-3`
   subst/site/year.
3. `--v0-mu-prior-mean 5e-4 --v0-mu-prior-stddev 2e-4 --v0-init-mutation-rate 1e-3` —
   should show `mu = 1 * 10^-3` subst/site/year, honoring the explicit override.
4. `--v0-mu-prior-alpha 0 --v0-mu-prior-beta 1000` (Jeffreys prior) — should show
   `mu = 1 * 10^-3` subst/site/year (improper prior, falls back to default).

The terminal output reports mutation rates in subst/site/year, so values can be compared
directly.

## Behavioral Impact Analysis

- **No behavior change for existing users** who don't set a mu prior: the default prior
  is `alpha=1, beta=0` (improper), so the `1e-3` default is preserved.
- **No behavior change for users who set `--v0-init-mutation-rate` explicitly**: the
  explicit value always wins.
- **New behavior only when** a proper prior is specified without an explicit init: the
  prior mean replaces `1e-3`.  This is strictly better — the initial value is already in
  the prior's typical region, reducing burn-in.
- **Edge case — `--v0-fix-mutation-rate` with a prior**: the user might set a prior and
  also fix the mutation rate, expecting it to be fixed at the prior mean.  This
  combination is unusual but would work correctly: `init_mu` would be the prior mean,
  and `fix_mutation_rate` would prevent it from changing.
