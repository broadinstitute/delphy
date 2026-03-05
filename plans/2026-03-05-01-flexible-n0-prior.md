# Flexible Inverse-Gamma Prior on Exponential Pop Model n0

The exponential population model in Delphy parameterizes the effective population size as
`N(t) = n0 * exp(g * (t - t0))`, where `t0` is the time of the last tip, `n0 = N(t0)` is
the effective population size at that reference time, and `g` is the exponential growth rate.
The parameter `n0` is stored internally in units of days (matching Delphy's internal time
unit); on the CLI it is presented in years.

## Current Prior Implementation

### n0 — Jeffreys prior (1/x)

**Prior calculation in `calc_cur_log_other_priors()`** (run.cpp:510-511):
```cpp
// pop_n0 - 1/x prior
log_prior -= std::log(exp_pop_model.pop_at_t0());
```

**MH move in `pop_size_move()`** (run.cpp:1221-1258): scale operator with scale factor
in [0.75, 1/0.75].  The prior ratio is `1/scale`:
```cpp
auto log_prior_new_over_prior_old = -std::log(scale);
```
The `log_other_priors_` incremental update is:
```cpp
log_other_priors_ += log_prior_new_over_prior_old;
```


## Design Discussion

### n0: from Jeffreys (1/x) to an Inverse-Gamma prior

#### Conjugacy analysis

The coalescent log-likelihood for a fixed tree (topology, times) with exponential population
N(t) = n0 * exp(g * (t - t0)) is:
```
log L = sum_coalescences log(1/N(t_i)) - integral (k choose 2) / N(t) dt
```

Writing **eta = 1/n0**, we get `1/N(t) = eta * exp(-g*(t - t0))`, so:
```
log L(eta) = C * log(eta) - eta * B
```
where:
- C = number of coalescent events (= number of tips - 1)
- B = integral of (k choose 2) * exp(-g*(t - t0)) dt (depends on g and tree, not on eta)

That's the kernel of a **Gamma distribution in eta**, and equivalently the kernel of an
**Inverse-Gamma distribution in n0**.  The Inverse-Gamma(alpha, beta) density on n0 is:
```
pi(n0) ~ n0^{-(alpha+1)} * exp(-beta / n0)
```
(equivalently, 1/n0 ~ Gamma(alpha, beta)).  Combining with the likelihood:
```
pi(n0 | tree, g) ~ n0^{-(C + alpha + 1)} * exp(-(B + beta) / n0)
                = InvGamma(C + alpha, B + beta)
```
which can be Gibbs sampled (by sampling eta ~ Gamma(C + alpha, B + beta) and setting
n0 = 1/eta) — just like mu and skygrid tau.

#### What does the current 1/x prior correspond to?

The current prior `pi(n0) ~ 1/n0 = n0^{-(0+1)} * exp(0)` is **Inverse-Gamma(alpha=0, beta=0)**,
the improper Jeffreys prior on n0.

#### Connection to Skygrid zero mode

Log-linear Skygrid with two knots reduces to an exponential population curve.  In
`skygrid_gammas_zero_mode_gibbs_move`, the overall population level (zero mode of
the gamma vector) is Gibbs sampled because the coalescent likelihood has this same
conjugate structure in 1/N.  Shifting all gammas by delta is equivalent to scaling
1/n0 by exp(delta).  So choosing an Inverse-Gamma prior on n0 keeps the door open
for Gibbs sampling the n0 posterior in the future.

#### Proposed prior

**Inverse-Gamma(alpha, beta) prior on n0**, with:
- Density: `pi(n0) ~ n0^{-(alpha+1)} * exp(-beta / n0)`
- Defaults: `alpha = 0, beta = 0` (Jeffreys 1/x prior on n0 = current behavior)
- Log-prior contribution:
  ```
  log pi(n0) = -(alpha + 1) * log(n0) - beta / n0
  ```
  With alpha = 0, beta = 0: `-log(n0)`, matching the current code's
  `log_prior -= std::log(exp_pop_model.pop_at_t0())`.
- Conjugate posterior: `InvGamma(C + alpha, B + beta)` on n0, or equivalently
  `Gamma(C + alpha, B + beta)` on 1/n0 — ready for future Gibbs sampling.
- Parallels the mu prior pattern: mu uses Gamma(1, 0) for uniform; here we use
  InvGamma(0, 0) for Jeffreys.


## Implementation Plan

Add `pop_inv_n0_prior_alpha` and `pop_inv_n0_prior_beta` parameters following the same pattern as
the mu prior parameters.  The prior on n0 is `pi(n0) ~ n0^{-(alpha+1)} * exp(-beta / n0)`.
Defaults: `alpha = 0, beta = 0` (Jeffreys 1/x prior, matching current behavior).

The naming uses `inv_n0` (rather than just `n0`) to make clear that alpha and beta are the
parameters of an *Inverse*-Gamma prior on n0, not a Gamma prior on n0.  Without this
distinction, users familiar with the analogous `mu_prior_alpha` / `mu_prior_beta` (which are
Gamma parameters on mu) would naturally assume these are also Gamma parameters, leading to
the wrong density.

### Step 1: `core/run.h` — Add accessors, setters, and private members

Add after `set_pop_model` (together with the other parameter getters/setters):
```cpp
auto pop_inv_n0_prior_alpha() const -> double { return pop_inv_n0_prior_alpha_; }
auto set_pop_inv_n0_prior_alpha(double pop_inv_n0_prior_alpha) -> void {
  pop_inv_n0_prior_alpha_ = pop_inv_n0_prior_alpha, invalidate_derived_quantities(); }
auto pop_inv_n0_prior_beta() const -> double { return pop_inv_n0_prior_beta_; }
auto set_pop_inv_n0_prior_beta(double pop_inv_n0_prior_beta) -> void {
  pop_inv_n0_prior_beta_ = pop_inv_n0_prior_beta, invalidate_derived_quantities(); }
```

Private members (after existing pop model members):
```cpp
double pop_inv_n0_prior_alpha_;
double pop_inv_n0_prior_beta_;
```

### Step 2: `core/run.cpp` — Constructor, prior calculation, and MH move

**Constructor:** add member initializers:
```cpp
pop_inv_n0_prior_alpha_{0.0},  // Default: Jeffreys 1/x prior on n0
pop_inv_n0_prior_beta_{0.0},
```

**`calc_cur_log_other_priors()`** (lines 510-511): replace the 1/x prior block:
```cpp
// pop_n0 - Inverse-Gamma(pop_inv_n0_prior_alpha_, pop_inv_n0_prior_beta_) prior
log_prior += -(pop_inv_n0_prior_alpha_ + 1) * std::log(exp_pop_model.pop_at_t0())
             - pop_inv_n0_prior_beta_ / exp_pop_model.pop_at_t0();
```
With default alpha=0, beta=0, this gives `-log(n0)`, matching current behavior.

**`pop_size_move()`** (lines 1224-1258): update the prior ratio.  With
`new_n0 = scale * old_n0`:
```cpp
// We use an Inverse-Gamma(alpha, beta) prior on n0:
//   pi(n0) ~ n0^{-(alpha+1)} * exp(-beta / n0)
// The default alpha=0, beta=0 gives the 1/x prior used in LeMieux et al (2021).
// A scale operator with scale factor [0.75, 1/0.75] proposes new_n0 = scale * old_n0.
//
//   log pi(new) - log pi(old) = -(alpha + 1) * log(scale) - beta * (1/new_n0 - 1/old_n0)

...

auto log_prior_new_over_prior_old =
    -(pop_inv_n0_prior_alpha_ + 1) * std::log(scale)
    - pop_inv_n0_prior_beta_ * (1.0 / new_n0 - 1.0 / old_n0);
```
With alpha=0, beta=0: `-(0+1)*log(scale) = -log(scale)`, matching current behavior.

### Step 3: `core/api.fbs` — Add fields to FlatBuffers schema

Add to the `Params` table using next free ids 40-41:
```
pop_inv_n0_prior_alpha: float64 (id: 40);  // alpha parameter of InvGamma prior on n0 [unitless, >= 0]
pop_inv_n0_prior_beta: float64 (id: 41);   // beta parameter of InvGamma prior on n0 [units: days, >= 0]
```
Both default to 0.0 (implicit FlatBuffers default), which matches the current
Jeffreys prior.  Old `.dphy` files that lack these fields will deserialize correctly.

Update "Next free id" comment to 42.

### Step 4: `core/api_generated.h` — Regenerate

Run `cd core && ../third-party/flatbuffers/build/flatc --cpp api.fbs`.

### Step 5: `core/api.cpp` — Serialization and deserialization

**Serialization:**
```cpp
params_builder.add_pop_inv_n0_prior_alpha(run.pop_inv_n0_prior_alpha());
params_builder.add_pop_inv_n0_prior_beta(run.pop_inv_n0_prior_beta());
```

**Deserialization:**
```cpp
run.set_pop_inv_n0_prior_alpha(api_params->pop_inv_n0_prior_alpha());
run.set_pop_inv_n0_prior_beta(api_params->pop_inv_n0_prior_beta());
```

### Step 6: `core/cmdline.cpp` — CLI options

Two alternative interfaces, mutually exclusive with each other.

**Interface A: alpha/beta**

```cpp
("v0-pop-inv-n0-prior-alpha",
 "[pop-model == exponential] Alpha parameter of the Inverse-Gamma prior on the effective "
 "population size n0: pi(n0) ~ n0^{-(alpha+1)} exp[-beta/n0].  "
 "Default alpha=0, beta=0 gives the Jeffreys 1/x prior.",
 cxxopts::value<double>()->default_value("0.0"))
("v0-pop-inv-n0-prior-beta",
 "[pop-model == exponential] Beta parameter of the Inverse-Gamma prior on the effective "
 "population size n0: pi(n0) ~ n0^{-(alpha+1)} exp[-beta/n0].  "
 "Default alpha=0, beta=0 gives the Jeffreys 1/x prior.  Units: years.",
 cxxopts::value<double>()->default_value("0.0"))
```

Note: the CLI presents beta in years (matching `--v0-init-final-pop-size`).
Internally, `beta_internal = beta_cli * 365` (since Delphy's n0 is in days).

**Interface B: mean/stddev**

```cpp
("v0-pop-n0-prior-mean",
 "[pop-model == exponential] Mean of the Inverse-Gamma prior on effective population size "
 "n0, in years (mean = beta / (alpha - 1)).  "
 "Must be specified together with --v0-pop-n0-prior-stddev.  "
 "Mutually exclusive with --v0-pop-inv-n0-prior-alpha / --v0-pop-inv-n0-prior-beta.",
 cxxopts::value<double>())
("v0-pop-n0-prior-stddev",
 "[pop-model == exponential] Standard deviation of the Inverse-Gamma prior on effective "
 "population size n0, in years (var = beta^2 / ((alpha-1)^2 (alpha-2))).  "
 "Must be specified together with --v0-pop-n0-prior-mean.  "
 "Mutually exclusive with --v0-pop-inv-n0-prior-alpha / --v0-pop-inv-n0-prior-beta.",
 cxxopts::value<double>())
```

Conversion from mean/stddev (in years) to internal alpha/beta (in days):
```
For InvGamma(alpha, beta): mean = beta/(alpha-1), var = beta^2/((alpha-1)^2*(alpha-2))
So: var/mean^2 = 1/(alpha-2)  =>  alpha = 2 + mean^2/var
    beta = mean * (alpha - 1)
    beta_internal = beta_years * 365
```

Note: the mean/stddev interface requires alpha > 2 for finite variance, which
constrains how diffuse the prior can be.  The alpha/beta interface is more general.

**Parsing:** follows the same mutual-exclusion pattern as the mu prior's alpha/beta
vs mean/stddev interfaces.  When using mean/stddev, set initial n0 to the prior mean
(as was done for the mu prior).  When using alpha/beta directly, set initial n0 to the
prior mean `beta / (alpha - 1)` if alpha > 1 and beta > 0; otherwise, fall back to the
current default initial value.

**Detection:** add the new option names to the list of exponential-only options
checked against wrong pop model.

### Step 7: `tools/delphy_wasm.cpp` — WASM getter/setter pairs

Add four `EMSCRIPTEN_KEEPALIVE extern "C"` functions near the exp pop model getters:
```cpp
auto delphy_run_get_pop_inv_n0_prior_alpha(Delphy_context& /*ctx*/, const Run& run) -> double;
auto delphy_run_set_pop_inv_n0_prior_alpha(Delphy_context& /*ctx*/, Run& run, double alpha) -> void;
auto delphy_run_get_pop_inv_n0_prior_beta(Delphy_context& /*ctx*/, const Run& run) -> double;
auto delphy_run_set_pop_inv_n0_prior_beta(Delphy_context& /*ctx*/, Run& run, double beta) -> void;
```

### Step 7b: `~/now/delphy-web/src/ts/pythia/delphy_api.ts` — TypeScript API bindings

Add corresponding entries near the exp pop model declarations in:
1. `DelphyCoreRaw` interface
2. Module bindings
3. `DelphyRun` class methods

### Step 8: `doc/dphy_file_format.md` — Update documentation

Add the two new fields to the `Params` table documentation.

### Step 9: `core/beasty_input.cpp` — BEAST XML export

When n0 has a non-Jeffreys prior (i.e., alpha != 0 or beta != 0), emit an
Inverse-Gamma prior instead of the current OneOnX/oneOnXPrior.

**BEAST 2 (2.6.2 and 2.7.7):**

When non-default:
```xml
<prior id="ePopSizePrior.t:input_alignment" name="distribution" x="@ePopSize.t:input_alignment">
  <InverseGamma id="InverseGamma.popSize" name="distr">
    <parameter ... name="alpha">alpha</parameter>
    <parameter ... name="beta">beta_years</parameter>
  </InverseGamma>
</prior>
```
When default (alpha=0, beta=0), keep the current `<OneOnX>`.

**BEAST X (10.5.0):**

When non-default, emit an `<invgammaPrior>` (or equivalent).
When default, keep the current `<oneOnXPrior>`.

### Step 10: Build and test

**10a.** Build debug, run unit tests.

**10b.** Short MCMC run with default priors (should match current behavior):
```bash
./build/debug/delphy \
  --v0-in-fasta ~/now/delphy-2026-paper-data/ebola-gire-2014/delphy_inputs/ebola.fasta \
  --v0-steps 10000
```

**10c.** Short MCMC run with non-default priors to exercise new code paths:
```bash
./build/debug/delphy \
  --v0-in-fasta ~/now/delphy-2026-paper-data/ebola-gire-2014/delphy_inputs/ebola.fasta \
  --v0-steps 10000 \
  --v0-pop-n0-prior-mean 100 --v0-pop-n0-prior-stddev 50
```

**10d.** Generate BEAST XML in all three formats with non-default priors and verify output.
Validate that BEAST2 and BEAST X can actually parse the generated XML by running them
briefly (e.g., with `--v0-steps 1` or a similar short-run option).  Runnable BEAST
executables are at:
- `~/github/CompEvol/beast2.6.2/bin/beast` (BEAST 2.6.2)
- `~/now/delphy-2026-paper-data/beast277` (BEAST 2.7.7)
- `~/now/delphy-2026-paper-data/beastX1050` (BEAST X 10.5.0)