# Flexible Inverse-Gamma Prior on Skygrid Average Population Level

The Skygrid population model parameterizes `log(N(t))` via gamma values at `M+1` knots
`(x_0, gamma_0), ..., (x_M, gamma_M)` where `x_0 < x_1 < ... < x_M`.  Between and
beyond knots, `log(N(t))` is either constant per interval (staircase) or linearly
interpolated (log-linear); see `core/pop_model.h` for full details.

The overall population level is captured by the "zero mode"
`gamma_bar = (1/(M+1)) * sum_{k=0}^M gamma_k`, whose geometric-mean population level
is `N_bar = exp(gamma_bar)`.  Currently, `gamma_bar` has an implicit uniform prior
(no constraint on the overall level).  The `skygrid_gammas_zero_mode_gibbs_move`
Gibbs-samples `I_bar = exp(-gamma_bar)` from a Gamma distribution derived from the
coalescent likelihood alone.

This plan adds a configurable `InvGamma(alpha, beta)` prior on
`N_bar = exp(gamma_bar)`, paralleling the `pop_inv_n0_prior_alpha/beta` parameters for
the exponential model (see `plans/2026-03-05-01-flexible-n0-prior.md`).  The default
parameters `(alpha=0, beta=0)` recover the current uniform-on-`gamma_bar` behavior.

The InvGamma prior on `N_bar` is conjugate with the coalescent likelihood, so the
zero-mode Gibbs move remains exact.  The prior does not constrain `gamma_bar` to any
particular range — `gamma_bar` can be any real number (negative, zero, or positive).

A user-friendly CLI interface allows specifying the desired mean and standard deviation
of `N_bar` (in years), which is converted to `(alpha, beta)`.


## Current Implementation

### Skygrid GMRF prior

The GMRF smoothness prior on the population curve is (Gill et al 2012, Eq. 13):

```
log pi_GMRF({gamma_k}) = sum_{k=1}^M [ (1/2)(log(tau) - log(2*pi))
                                        - (tau/2) * (gamma_k - gamma_{k-1})^2 ]
```

This prior constrains only the *differences* between consecutive `gamma_k` values.
It does not constrain the overall level `gamma_bar`.

### Zero-mode Gibbs move (run.cpp:1971-2119)

Decomposes `{gamma_k}` into `gamma_bar` + zero-sum offsets `{gamma'_k}`.  The GMRF
smoothness prior only constrains differences `(gamma_k - gamma_{k-1})`, so it does not
constrain `gamma_bar`.  The coalescent likelihood, expressed in
`I_bar = exp(-gamma_bar)`, gives:

```
exp[-U_coal] ~ I_bar^{N_inner} * exp[-B * I_bar]
```

where:
- `N_inner = (num_tips - 1)` is the number of internal (coalescent) nodes,
- `c` ranges over the `C` cells of width `Delta` in the scalable coalescent prior,
- `k_c` is the average number of active lineages in cell `c`,
- `N_c` is the average of `N(t)` over cell `c`,
- `Ntwiddle_c = N_c * I_bar` factors out the dependence on `I_bar`,
- `B = sum_{c=0}^{C-1} (1/2) * Delta * k_c * (k_c - 1) / Ntwiddle_c`.

Combined with the Jacobian from `gamma_bar` -> `I_bar` (factor `1/I_bar`), the
posterior is:

```
P(I_bar) ~ I_bar^{N_inner - 1} * exp[-B * I_bar]  =>  Gamma(N_inner, B)
```

This is Gibbs-sampled directly (run.cpp:2061-2066).  After Gibbs sampling, an MH
acceptance criterion corrects for the low-gamma barrier, which is NOT absorbed into
the Gamma distribution (run.cpp:2074-2117).

### HMC move (run.cpp:1350-1969)

The potential energy `U = U_coal + U_prior` drives the HMC dynamics.  Currently,
`U_prior` (computed by `calc_U_prior_from_gamma_k_s`, run.cpp:1706-1724) includes:

1. GMRF smoothness: `(tau/2) * sum_{k=1}^M (gamma_k - gamma_{k-1})^2`
2. Low-gamma barrier (quadratic penalty for `gamma_k` below a threshold)

There is no `U_prior` term constraining the overall level.

The forces `f_k = -dU/d(gamma_k)` (run.cpp:1746-1798) likewise have no zero-mode
component from `U_prior`.

### `calc_cur_log_other_priors` (run.cpp:525-537)

For the Skygrid branch, includes the tau prior (if tau is being inferred) and delegates
to `calc_cur_skygrid_gmrf_prior()` (run.cpp:542-595), which computes the GMRF log-prior
plus the low-gamma barrier.  No prior on the overall level.


## Conjugacy Analysis

### What prior does "uniform on `gamma_bar`" correspond to?

A uniform prior on `gamma_bar` is equivalent to:
- In `N_bar = exp(gamma_bar)`: `pi(N_bar) ~ 1/N_bar` (Jeffreys prior)
- In `I_bar = exp(-gamma_bar)`: `pi(I_bar) ~ 1/I_bar`

This is `InvGamma(alpha=0, beta=0)` on `N_bar`, confirming that `(alpha=0, beta=0)` is
the right default.

Note: the default Jeffreys `1/n0` prior on `n0` in the exponential population model
also corresponds to a uniform prior on `log(n0)` — the exact same kind of prior as the
current uniform prior on `gamma_bar` here.

### Inverse-Gamma prior on `N_bar`

An `InvGamma(alpha, beta)` prior on `N_bar` has density:

```
pi(N_bar) ~ N_bar^{-(alpha+1)} * exp(-beta / N_bar)
```

In `gamma_bar` coordinates (including the Jacobian
`|dN_bar/d(gamma_bar)| = exp(gamma_bar)`):

```
pi_{gamma_bar}(gamma_bar) = pi_{N_bar}(exp(gamma_bar)) * exp(gamma_bar)
    ~ exp(gamma_bar)^{-(alpha+1)} * exp(-beta * exp(-gamma_bar)) * exp(gamma_bar)
    = exp(-alpha * gamma_bar) * exp(-beta * exp(-gamma_bar))
    = I_bar^alpha * exp(-beta * I_bar)
```

So `log pi(gamma_bar) = -alpha * gamma_bar - beta * exp(-gamma_bar) + const`.

With `alpha=0, beta=0`: `log pi = 0` (uniform on `gamma_bar`).  Correct.

**No constraint on `gamma_bar`:** Both terms `exp(-alpha * gamma_bar)` and
`exp(-beta * exp(-gamma_bar))` are well-defined for all real `gamma_bar` (negative,
zero, or positive), since `exp(-gamma_bar) > 0` always.  The InvGamma prior does not
restrict the domain of `gamma_bar`.

### Updated zero-mode Gibbs posterior

The full `I_bar` posterior becomes:

```
P(I_bar) ~ [coalescent] * [Jacobian] * [InvGamma prior on N_bar, in I_bar coords]
         = I_bar^{N_inner} * exp(-B*I_bar) * I_bar^{-1} * I_bar^{alpha} * exp(-beta*I_bar)
         = I_bar^{N_inner + alpha - 1} * exp(-(B + beta) * I_bar)
```

This is **`Gamma(N_inner + alpha, B + beta)`** in `I_bar`.

Check: `alpha=0, beta=0` -> `Gamma(N_inner, B)`.  Matches current code.

The conjugacy is preserved: the zero-mode Gibbs move remains exact.

### Connection to exponential model

For the exponential model with `eta = 1/n0`, the `InvGamma(alpha, beta)` prior on `n0`
yields a `Gamma(C + alpha, B_exp + beta)` posterior on `eta` (see n0 plan).  The Skygrid
result is exactly analogous, with `N_inner` playing the role of `C` and `B` playing the
role of `B_exp`.


## Design Discussion

### Why `N_bar` and not `exp(gamma_M)`?

An alternative would be to put an InvGamma prior on `exp(gamma_M)`, the population at
the last knot (analogous to `n0` in the exponential model, which is the population at
the time of the last tip).  This would also preserve the Gibbs structure of the
zero-mode move: since `gamma_M = gamma'_M + gamma_bar` where `gamma'_M` is fixed
during the move, the prior contributes `I_bar^alpha * exp(-beta * exp(-gamma'_M) * I_bar)`
to the `I_bar` posterior, yielding `Gamma(N_inner + alpha, B + beta * exp(-gamma'_M))`.

However, the HMC forces from such a prior would be concentrated entirely on `gamma_M`
(`f_M = -alpha + beta * exp(-gamma_M)`, `f_k = 0` for `k != M`).  We prefer the prior
on `N_bar` precisely because it distributes small, equal forces across all knots instead
of concentrating all the force on a single knot.

More generally, one could put independent `InvGamma(alpha_k, beta_k)` priors on every
`exp(gamma_k)` simultaneously and the zero-mode Gibbs move would still work: the
posterior on `I_bar` would be
`Gamma(N_inner + sum_k alpha_k, B + sum_k beta_k * exp(-gamma'_k))`.  This could be
interesting for future work but is beyond the scope of this plan.

### Naming

The exponential model uses `pop_inv_n0_prior_alpha/beta` where `inv_n0` signals that
alpha and beta parameterize an *Inverse*-Gamma on `n0`.  For the Skygrid, the analogous
quantity is `N_bar = exp(gamma_bar)` (the geometric mean population across knots).

Proposed names: `skygrid_inv_nbar_prior_alpha` and `skygrid_inv_nbar_prior_beta`.

This follows the existing `skygrid_` prefix convention
(cf. `skygrid_tau_prior_alpha/beta`) and uses `inv_nbar` to signal an Inverse-Gamma on
`N_bar`.

### Impact on HMC forces

The InvGamma prior adds a potential energy term:

```
U_invgamma = alpha * gamma_bar + beta * exp(-gamma_bar)
           = alpha * (1/(M+1)) * sum_k gamma_k
             + beta * exp(-(1/(M+1)) * sum_k gamma_k)
```

The force contribution on each `gamma_k` is identical (the prior only constrains the
zero mode):

```
f_k_invgamma = -dU_invgamma / d(gamma_k)
             = (-alpha + beta * exp(-gamma_bar)) / (M+1)
```

This is uniform across all `k`.  Note that `f_k_invgamma = 0` when
`exp(-gamma_bar) = alpha / beta`, i.e., at the mode of the prior in `gamma_bar`
coordinates.  (This value coincides with the mean of the implied `Gamma(alpha, beta)`
distribution on `1/N_bar = exp(-gamma_bar)`, but is its mode shifted by the Jacobian:
the prior in `gamma_bar` coordinates is proportional to
`eta^alpha * exp(-beta * eta)` where `eta = exp(-gamma_bar)`, whose mode is at
`eta = alpha / beta`, whereas the `Gamma(alpha, beta)` pdf
`eta^{alpha-1} * exp(-beta * eta)` has its mode at `(alpha - 1) / beta`.)

The second derivative
`d²U_invgamma / d(gamma_k)² = beta * exp(-gamma_bar) / (M+1)²`
is assumed small (divided by `(M+1)²`).  In principle, `exp(-gamma_bar)` can grow very
large when `gamma_bar` is very negative (very small population), which could make this
term significant.  However, in practice the low-gamma barrier already dominates in that
regime, preventing `gamma_bar` from reaching values where this matters.  If needed, the
mass matrix could be augmented by adding
`beta * exp(-gamma_bar) / (M+1)²` to each `m_k`, but we do not expect this to be
necessary for realistic datasets.  **No mass matrix changes for now.**

### Impact on `U_prior` tracking

The HMC move computes `old_U_prior` and `new_U_prior` and updates `log_other_priors_`
via `log_other_priors_ += (-new_U_prior) - (-old_U_prior)` (run.cpp:1960).

The zero-mode Gibbs move updates `log_other_priors_` for the barrier change
(run.cpp:2117).

Both must now also account for the InvGamma prior change.  For the HMC move, this means
`U_invgamma` is folded into `calc_U_prior_from_gamma_k_s`.  For the zero-mode move, the
InvGamma change must be computed explicitly and added to `log_other_priors_`.

### Impact on debug force check

The HMC move has a debug mode (run.cpp:1857-1883) that numerically verifies forces via
finite differences using `calc_U_prior_from_gamma_k_s`.  Since the InvGamma term will be
added to that function, the debug check will automatically verify the new force term.
We must ensure `gamma_bar` is recomputed inside `calc_U_prior_from_gamma_k_s` from the
current `gamma_k` values (not cached), so that finite-difference perturbations of
individual `gamma_k`'s produce correct numerical derivatives.


## Implementation Plan

### Step 0: `core/pop_model.h` — Add `gamma_bar()` and `N_bar()` to `Skygrid_pop_model`

Add convenience methods after the existing derived-quantity accessors (~line 191):

```cpp
// gamma_bar = mean of gamma_k's = "zero mode" of the log-population curve
auto gamma_bar() const -> double {
  auto result = 0.0;
  for (auto k = 0; k <= M(); ++k) {
    result += gamma(k);
  }
  return result / (M() + 1);
}
// N_bar = geometric mean population across knots = exp(gamma_bar)
auto N_bar() const -> double { return std::exp(gamma_bar()); }
```

These are used in `calc_cur_skygrid_gmrf_prior()`, the log output in `delphy.cpp` and
`delphy_ui.cpp`, and anywhere else a `Skygrid_pop_model` object is available.

Note: the HMC lambdas and parts of the zero-mode Gibbs move operate on local `gamma_k`
arrays (not yet wrapped in a `Skygrid_pop_model`), so those sites still compute
`gamma_bar` inline from the array.

### Step 1: `core/run.h` — Add accessors, setters, and private members

Add after the `skygrid_low_gamma_barrier_scale` accessors (together with other Skygrid
parameter getters/setters):

```cpp
auto skygrid_inv_nbar_prior_alpha() const -> double { return skygrid_inv_nbar_prior_alpha_; }
auto set_skygrid_inv_nbar_prior_alpha(double skygrid_inv_nbar_prior_alpha) -> void {
  skygrid_inv_nbar_prior_alpha_ = skygrid_inv_nbar_prior_alpha, invalidate_derived_quantities(); }
auto skygrid_inv_nbar_prior_beta() const -> double { return skygrid_inv_nbar_prior_beta_; }
auto set_skygrid_inv_nbar_prior_beta(double skygrid_inv_nbar_prior_beta) -> void {
  skygrid_inv_nbar_prior_beta_ = skygrid_inv_nbar_prior_beta, invalidate_derived_quantities(); }
```

Private members (after `skygrid_low_gamma_barrier_scale_`):

```cpp
double skygrid_inv_nbar_prior_alpha_;
double skygrid_inv_nbar_prior_beta_;
```

### Step 2: `core/run.cpp` — Constructor

Add member initializers (after `skygrid_low_gamma_barrier_scale_`):

```cpp
skygrid_inv_nbar_prior_alpha_{0.0},  // Default: uniform prior on gamma_bar (= Jeffreys 1/N_bar)
skygrid_inv_nbar_prior_beta_{0.0},
```

### Step 3: `core/run.cpp` — `calc_cur_skygrid_gmrf_prior()`

Add at the start of the function, before the GMRF loop (after line 550):

```cpp
// InvGamma(alpha, beta) prior on N_bar = exp(gamma_bar)
//   log pi(gamma_bar) = -alpha * gamma_bar - beta * exp(-gamma_bar) + const
// Default alpha=0, beta=0: no contribution (uniform on gamma_bar).
auto gamma_bar = skygrid_pop_model.gamma_bar();
log_prior += -skygrid_inv_nbar_prior_alpha_ * gamma_bar
             - skygrid_inv_nbar_prior_beta_ * std::exp(-gamma_bar);
```

### Step 4: `core/run.cpp` — HMC move, `calc_U_prior_from_gamma_k_s`

Add the InvGamma contribution to `calc_U_prior_from_gamma_k_s` (after the GMRF and
barrier terms, before `return U_prior`):

```cpp
// InvGamma(alpha, beta) prior on N_bar = exp(gamma_bar)
//   U_invgamma = alpha * gamma_bar + beta * exp(-gamma_bar)
{
  auto gamma_bar = 0.0;
  for (auto k = 0; k <= M; ++k) {
    gamma_bar += gamma_k[k];
  }
  gamma_bar /= (M + 1);
  U_prior += skygrid_inv_nbar_prior_alpha_ * gamma_bar
             + skygrid_inv_nbar_prior_beta_ * std::exp(-gamma_bar);
}
```

### Step 5: `core/run.cpp` — HMC move, `calc_f_k_s_from_gamma_k_s`

Add the InvGamma force contribution.  This must be computed once (not per-k) and then
distributed equally.  Add after the existing per-k force loop:

```cpp
// InvGamma(alpha, beta) prior on N_bar = exp(gamma_bar): force on each gamma_k
//   f_k_invgamma = (-alpha + beta * exp(-gamma_bar)) / (M+1)
{
  auto gamma_bar = 0.0;
  for (auto k = 0; k <= M; ++k) {
    gamma_bar += gamma_k[k];
  }
  gamma_bar /= (M + 1);
  auto f_invgamma = (-skygrid_inv_nbar_prior_alpha_
                     + skygrid_inv_nbar_prior_beta_ * std::exp(-gamma_bar))
                    / (M + 1);
  for (auto k = 0; k <= M; ++k) {
    f_k[k] += f_invgamma;
  }
}
```

Note: adding the same force to all `k` is correct since the prior only constrains the
zero mode.

### Step 6: `core/run.cpp` — Zero-mode Gibbs move

**6a. Update Gibbs distribution** (run.cpp:2061-2063):

Change from:

```cpp
auto I_bar_dist = std::gamma_distribution<double>{
    static_cast<double>(N_inner),
    1.0 / B};
```

to:

```cpp
auto I_bar_dist = std::gamma_distribution<double>{
    static_cast<double>(N_inner) + skygrid_inv_nbar_prior_alpha_,
    1.0 / (B + skygrid_inv_nbar_prior_beta_)};
```

With `alpha=0, beta=0`: unchanged.

**6b. Update `log_other_priors_`** (run.cpp:2117):

Currently only accounts for barrier changes.  Must also account for the InvGamma prior
change.  Replace:

```cpp
log_other_priors_ += new_log_prior_bound - old_log_prior_bound;
```

with:

```cpp
// InvGamma prior on N_bar: log pi = -alpha * gamma_bar - beta * exp(-gamma_bar)
auto old_gamma_bar = 0.0;
auto new_gamma_bar = 0.0;
for (auto k = 0; k <= M; ++k) {
  old_gamma_bar += old_pop_model->gamma(k);
  new_gamma_bar += new_gamma_k[k];
}
old_gamma_bar /= (M + 1);
new_gamma_bar /= (M + 1);
auto old_inv_nbar_prior = -skygrid_inv_nbar_prior_alpha_ * old_gamma_bar
                          - skygrid_inv_nbar_prior_beta_ * std::exp(-old_gamma_bar);
auto new_inv_nbar_prior = -skygrid_inv_nbar_prior_alpha_ * new_gamma_bar
                          - skygrid_inv_nbar_prior_beta_ * std::exp(-new_gamma_bar);

log_other_priors_ += (new_inv_nbar_prior - old_inv_nbar_prior)
                   + (new_log_prior_bound - old_log_prior_bound);
```

**6c. Update `delta_log_posterior` and MH check** (run.cpp:2091):

The zero-mode move is Gibbs for the coalescent + InvGamma prior (those are absorbed
into the Gamma distribution parameters).  The barrier is NOT absorbed, so it remains
a correction.  This is already handled: `delta_log_posterior` tracks barrier changes,
and the MH check applies them.  No change needed here.

### Step 7: `core/api.fbs` — Add fields to FlatBuffers schema

Add to the `Params` table using next free ids 46-47:

```
skygrid_inv_nbar_prior_alpha: float64 (id: 46);  // alpha of InvGamma prior on N_bar [unitless, >= 0]
skygrid_inv_nbar_prior_beta: float64 (id: 47);   // beta of InvGamma prior on N_bar [units: days, >= 0]
```

Update "Next free id" comment to 48.

Both default to 0.0 (implicit FlatBuffers default), matching current behavior.

### Step 8: `core/api_generated.h` — Regenerate

Run `cd core && ../third-party/flatbuffers/build/flatc --cpp api.fbs`.

### Step 9: `core/api.cpp` — Serialization and deserialization

**Serialization** (near existing skygrid fields, ~line 273):

```cpp
params_builder.add_skygrid_inv_nbar_prior_alpha(run.skygrid_inv_nbar_prior_alpha());
params_builder.add_skygrid_inv_nbar_prior_beta(run.skygrid_inv_nbar_prior_beta());
```

**Deserialization** (near existing skygrid fields, ~line 392):

```cpp
run.set_skygrid_inv_nbar_prior_alpha(api_params->skygrid_inv_nbar_prior_alpha());
run.set_skygrid_inv_nbar_prior_beta(api_params->skygrid_inv_nbar_prior_beta());
```

### Step 10: `core/cmdline.cpp` — CLI options

Two alternative interfaces, mutually exclusive with each other.

**Interface A: alpha/beta**

```cpp
("v0-skygrid-inv-nbar-prior-alpha",
 "[pop-model == skygrid] Alpha parameter of the Inverse-Gamma prior on N_bar, the "
 "geometric-mean effective population size across Skygrid knots "
 "(N_bar = exp(mean of log-population values gamma_k)).  "
 "pi(N_bar) ~ N_bar^{-(alpha+1)} exp[-beta/N_bar].  "
 "Default alpha=0, beta=0 gives a uniform prior on log(N_bar) (= Jeffreys 1/N_bar prior).",
 cxxopts::value<double>()->default_value("0.0"))
("v0-skygrid-inv-nbar-prior-beta",
 "[pop-model == skygrid] Beta parameter of the Inverse-Gamma prior on N_bar, the "
 "geometric-mean effective population size across Skygrid knots "
 "(N_bar = exp(mean of log-population values gamma_k)).  "
 "pi(N_bar) ~ N_bar^{-(alpha+1)} exp[-beta/N_bar].  "
 "Default alpha=0, beta=0 gives a uniform prior on log(N_bar) (= Jeffreys 1/N_bar prior).  "
 "Units: years.",
 cxxopts::value<double>()->default_value("0.0"))
```

Note: CLI presents beta in years.  Internally, `beta_internal = beta_cli * 365`.

**Interface B: mean/stddev**

```cpp
("v0-skygrid-nbar-prior-mean",
 "[pop-model == skygrid] Mean of the Inverse-Gamma prior on N_bar, the geometric-mean "
 "effective population size across Skygrid knots, in years "
 "(mean = beta / (alpha - 1)).  "
 "Must be specified together with --v0-skygrid-nbar-prior-stddev.  "
 "Mutually exclusive with --v0-skygrid-inv-nbar-prior-alpha / "
 "--v0-skygrid-inv-nbar-prior-beta.",
 cxxopts::value<double>())
("v0-skygrid-nbar-prior-stddev",
 "[pop-model == skygrid] Standard deviation of the Inverse-Gamma prior on N_bar, the "
 "geometric-mean effective population size across Skygrid knots, in years "
 "(var = beta^2 / ((alpha-1)^2 (alpha-2))).  "
 "Must be specified together with --v0-skygrid-nbar-prior-mean.  "
 "Mutually exclusive with --v0-skygrid-inv-nbar-prior-alpha / "
 "--v0-skygrid-inv-nbar-prior-beta.",
 cxxopts::value<double>())
```

Conversion from mean/stddev (in years) to internal alpha/beta (in days):

```
For InvGamma(alpha, beta):
  mean = beta / (alpha - 1),  var = beta^2 / ((alpha-1)^2 * (alpha-2))
So:
  var / mean^2 = 1 / (alpha - 2)  =>  alpha = 2 + mean^2 / var
  beta = mean * (alpha - 1)
  beta_internal = beta_years * 365
```

Note: the mean/stddev interface requires `alpha > 2` for finite variance, which
constrains how diffuse the prior can be.  The alpha/beta interface is more general.

When using mean/stddev, set initial `gamma_bar` to `log(mean_years * 365)` (prior mean
in internal units).  When using alpha/beta directly with `alpha > 1` and `beta > 0`,
set initial `gamma_bar` to `log(beta / (alpha - 1))` (prior mean); otherwise, keep the
current default (`log(3 * 365) ≈ 6.999`).

**Detection:** add new option names to `has_skygrid_pop_model_parameters`
(cmdline.cpp:602).

**Implementation order:** Process the InvGamma prior parameters *before* creating the
`Skygrid_pop_model`.  Compute `inv_nbar_prior_alpha`, `inv_nbar_prior_beta_days`, and
`init_nbar_days` first, then sample a random gamma_k trajectory with the configured tau
precision, shift its mean to `log(init_nbar_days)`, and create the `Skygrid_pop_model`
only once at the end.  This avoids creating the model and then recreating it.

### Step 11: `tools/delphy_wasm.cpp` — WASM getter/setter pairs

Add four `EMSCRIPTEN_KEEPALIVE extern "C"` functions near the skygrid getters:

```cpp
auto delphy_run_get_skygrid_inv_nbar_prior_alpha(Delphy_context& /*ctx*/, const Run& run) -> double;
auto delphy_run_set_skygrid_inv_nbar_prior_alpha(Delphy_context& /*ctx*/, Run& run, double alpha) -> void;
auto delphy_run_get_skygrid_inv_nbar_prior_beta(Delphy_context& /*ctx*/, const Run& run) -> double;
auto delphy_run_set_skygrid_inv_nbar_prior_beta(Delphy_context& /*ctx*/, Run& run, double beta) -> void;
```

### Step 11b: `~/now/delphy-web/src/ts/pythia/delphy_api.ts` — TypeScript API bindings

Add corresponding entries near the skygrid declarations in:
1. `DelphyCoreRaw` interface
2. Module bindings
3. `DelphyRun` class methods

### Step 11c: `tools/delphy.cpp` and `tools/delphy_ui.cpp` — Output `N_bar` and fix units

When the population model is Skygrid, output `N_bar = exp(gamma_bar)` (in years) in the
periodic log lines.  This makes it straightforward to verify that the InvGamma prior is
being applied properly (e.g., by running with a tight prior whose mean is far from the
unconstrained value, and checking that `N_bar` converges to the requested prior mean).

While we're at it, fix the units of all population-related quantities in the log output.
Currently, `n0`, `g`, and `gamma_k` are displayed in internal units (days) without any
unit labels.  The CLI accepts `n0` in years, `g` in e-foldings/year, and `beta` in years,
then converts internally to days — but the log output displays raw internal values.  This
is confusing.  Change all population-related outputs to use years with explicit unit labels:

- `n0`: days -> years (divide by 365), add `yr` label
- `g`: e-foldings/day -> e-foldings/year (multiply by 365), add `e-fold/yr` label
- `gamma_k`: `log(N_days)` -> `log(N_years)` (subtract `log(365)`), add annotation
- `N_bar`: display in years (divide by 365), add `yr` label

**`tools/delphy.cpp`:**

Change `n0` and `g` output (~lines 79-80) from:

```cpp
std::cerr << absl::StreamFormat("n0 = %.1f, ", exp_pop_model.pop_at_t0())
          << absl::StreamFormat("g = %.4f, ", exp_pop_model.growth_rate());
```

to:

```cpp
std::cerr << absl::StreamFormat("n0 = %.4g yr, ", exp_pop_model.pop_at_t0() / 365.0)
          << absl::StreamFormat("g = %.4f e-fold/yr, ", exp_pop_model.growth_rate() * 365.0);
```

Add a Skygrid branch after the Exp_pop_model branch (~line 81):

```cpp
} else if (typeid(pop_model) == typeid(Skygrid_pop_model)) {
  const auto& skygrid_pop_model = static_cast<const Skygrid_pop_model&>(run.pop_model());
  std::cerr << absl::StreamFormat("N_bar = %.4g yr, ", skygrid_pop_model.N_bar() / 365.0);
}
```

**`tools/delphy_ui.cpp`:**

Change `n0` and `g` output (~lines 615-616) from:

```cpp
std::cerr << absl::StreamFormat("n0 = %.4g, ", exp_pop_model.pop_at_t0())
          << absl::StreamFormat("g = %.4g, ", exp_pop_model.growth_rate());
```

to:

```cpp
std::cerr << absl::StreamFormat("n0 = %.4g yr, ", exp_pop_model.pop_at_t0() / 365.0)
          << absl::StreamFormat("g = %.4g e-fold/yr, ", exp_pop_model.growth_rate() * 365.0);
```

Change `gamma_k` output (~lines 620-625) to display in log-years instead of log-days.
Change from:

```cpp
std::cerr << absl::StreamFormat("tau = %.4g, ", ui_run->skygrid_tau());
std::cerr << absl::StreamFormat("gamma_k = %s---[", to_iso_date(skygrid_pop_model.x_lo()));
for (auto k = 0; k != skygrid_pop_model.M(); ++k) {
  if (k != 0) { std::cerr << ", "; }
  std::cerr << absl::StreamFormat("%.4g", skygrid_pop_model.gamma(k));
}
std::cerr << absl::StreamFormat("]---%s, ", to_iso_date(skygrid_pop_model.x_hi()));
```

to:

```cpp
std::cerr << absl::StreamFormat("tau = %.4g, ", ui_run->skygrid_tau());
auto log_365 = std::log(365.0);
std::cerr << absl::StreamFormat("gamma_k (log yr) = %s---[", to_iso_date(skygrid_pop_model.x_lo()));
for (auto k = 0; k != skygrid_pop_model.M(); ++k) {
  if (k != 0) { std::cerr << ", "; }
  std::cerr << absl::StreamFormat("%.4g", skygrid_pop_model.gamma(k) - log_365);
}
std::cerr << absl::StreamFormat("]---%s, ", to_iso_date(skygrid_pop_model.x_hi()));
```

After the `gamma_k` output, add `N_bar`:

```cpp
std::cerr << absl::StreamFormat("N_bar = %.4g yr, ", skygrid_pop_model.N_bar() / 365.0);
```

### Step 12: `doc/dphy_file_format.md` — Update documentation

Add the two new fields to the `Params` table documentation.  Update the next free
field id to 48.  Keep this file exactly in sync with `api.fbs`.

### Step 13: `core/beasty_input.cpp` — BEAST XML export

BEAST 2 (both 2.6.2 and 2.7.7) does not support Skygrid, so no change needed there.

BEAST X (10.5.0):  The approach is to define `N_bar` as a derived statistic from
`skygrid.logPopSize`, then apply an `invgammaPrior` to it.  Concretely:

1. Define a `meanStatistic` (or `productStatistic`/`sumStatistic` with appropriate
   scaling) to compute `gamma_bar = mean(skygrid.logPopSize)`.
2. Apply `exp()` to get `N_bar = exp(gamma_bar)`.
3. Apply `<invgammaPrior shape="alpha" scale="beta">` to `N_bar`.

The stumbling block is step 2: BEAST X may not provide a built-in way to apply an
exponential function when defining a composite statistic.  If no `exp()` transform is
available, emit an XML comment explaining that the InvGamma prior on `N_bar` cannot be
exactly represented.  Investigate BEAST X's `TransformedParameter`,
`ExponentialTransform`, or similar classes before concluding this is impossible.

When default (`alpha=0, beta=0`): no prior element needed (matches BEAST X's default
of no zero-mode prior in `gmrfSkyGridLikelihood`).

### Step 14: Build and test

**14a.** Build debug, run unit tests.

**14b.** Short MCMC run with default priors (should match current behavior):

```bash
./build/debug/delphy \
  --v0-in-fasta ~/now/delphy-2026-paper-data/ebola-gire-2014/delphy_inputs/ebola.fasta \
  --v0-pop-model skygrid --v0-skygrid-num-parameters 10 --v0-skygrid-cutoff 1.5 \
  --v0-steps 1000
```

**14c.** Short MCMC run with non-default priors to exercise new code paths:

```bash
./build/debug/delphy \
  --v0-in-fasta ~/now/delphy-2026-paper-data/ebola-gire-2014/delphy_inputs/ebola.fasta \
  --v0-pop-model skygrid --v0-skygrid-num-parameters 10 --v0-skygrid-cutoff 1.5 \
  --v0-steps 1000 \
  --v0-skygrid-nbar-prior-mean 100 --v0-skygrid-nbar-prior-stddev 50
```

**14d.** Enable `debug_hmc = true` temporarily to verify that the finite-difference
force check passes with non-default InvGamma prior parameters.  Use a release build
and run for 1000000 steps.  Verify from the log output that the debug check occurred
many times.

**14e.** Generate BEAST X XML with non-default priors and verify the output.
