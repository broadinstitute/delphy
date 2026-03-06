# Flexible Laplace Prior on Exponential Pop Model Growth Rate g

The exponential population model in Delphy parameterizes the effective population size as
`N(t) = n0 * exp(g * (t - t0))`, where `t0` is the time of the last tip, `n0 = N(t0)` is
the effective population size at that reference time, and `g` is the exponential growth rate.
The parameter `g` is stored internally in units of e-foldings / day (matching Delphy's
internal time unit); on the CLI and in BEAST XML it is presented in e-foldings / year.

Currently, g has a hard-coded Laplace prior with location `mu = 0.001/year` and
scale `30.701135/year`.  This plan makes the prior's location and scale configurable,
and adds optional lower and upper bounds on g.

## Current Prior Implementation

### g — Laplace prior with hard-coded location and scale

**Prior calculation in `calc_cur_log_other_priors()`** (run.cpp:516-520):
```cpp
// pop_g - Laplace prior on the growth rate, with mu 0.001/365 and scale 30.701135/365
//      pdf(g; mu, scale) = 1/(2*scale) exp(-|g - mu| / scale)
const auto mu_g = 0.001 / 365.0;
const auto scale_g = 30.701135 / 365.0;
log_prior += -std::abs(exp_pop_model.growth_rate() - mu_g) / scale_g - std::log(2 * scale_g);
```

**MH move in `pop_growth_rate_move()`** (run.cpp:1265-1304): random walk operator
with window size `1.0 / 365.0`.  Same hard-coded `mu = 0.001/365` and `scale = 30.701135/365`.

**BEAST XML export** (beasty_input.cpp): all three BEAST versions emit a Laplace prior
with `mu = 0.001` (per year) and `scale = 30.701135` (per year).

### Units

Internally, g is in units of e-foldings / day.  On the CLI and in BEAST XML, g is in
e-foldings / year.  Conversion: `g_internal = g_cli / 365`.

## Design Discussion

### No conjugacy for g

The coalescent log-likelihood with exponential population N(t) = n0 exp(g(t - t0)) is:
```
log L = C log(eta) - g sum_i(t_i - t0) - eta sum_j (k_j choose 2) integral exp(-g(t - t0)) dt
```
The dependence on g is both linear (from `log(1/N(t_i))` terms) and nonlinear
(from `exp(-g(t-t0))/g` integrals), so the likelihood is not in an exponential family
in g and **there is no conjugate prior**.

### Proposed changes

Make the Laplace prior's **location** (mu) and **scale** configurable, and add
optional **lower and upper bounds** on g (defaulting to unconstrained in both
directions).

Defaults: `mu = 0.001` (per year), `scale = 30.701135` (per year), matching the current
hard-coded values from BEAUti2 2.6.2 / LeMieux et al (2021).  Note: the origin of
the specific value `30.701135` is unclear — it appears in BEAUti2 2.6.2 as the default
Laplace scale for growth rate priors, but the authors never explain its derivation.
BEAUti2 2.7.7 changed this default to `0.5`, a much more informative prior.  With
`scale = 30.701135`, the standard deviation is `scale * sqrt(2) ≈ 43.4` e-foldings/year,
which is effectively uninformative.

The bounds (`--v0-pop-growth-rate-min`, `--v0-pop-growth-rate-max`) are enforced by
rejecting proposals that violate them in `pop_growth_rate_move()`.  The prior density
is not renormalized (the truncation just acts as a hard constraint on the parameter
space); this is equivalent to multiplying the prior by an indicator function
`I(g_min <= g <= g_max)`.

### Special cases of the truncated Laplace

The Laplace prior with bounds reduces to simpler distributions in some cases:

- **g_min = mu, g_max = +∞**: the truncated Laplace becomes a shifted Exponential
  distribution starting at mu, with mean = mu + scale.  In particular, when **mu = 0
  and g_min = 0**, this is `Exponential(rate = 1/scale)` with mean = scale.

- **g_min = -∞, g_max = mu**: mirror image — a reflected Exponential ending at mu.

- **g_min = -∞, g_max = +∞** (no bounds): standard Laplace(mu, scale) with mean = mu
  and stddev = scale * sqrt(2).

- **g_min = g_max**: effectively fixing g to that value.

A CLI shorthand `--v0-pop-g-prior-exponential-with-mean` configures an Exponential prior on g:
- If mean > 0: sets mu = 0, scale = mean, g_min = 0, g_max = +∞.
  This gives Exponential(rate = 1/mean) on g ≥ 0 with E[g] = mean.
- If mean < 0: sets mu = 0, scale = |mean|, g_min = -∞, g_max = 0.
  This gives a reflected Exponential on g ≤ 0 with E[g] = mean.
- mean = 0 is rejected (would give scale = 0, which is invalid).

This option is mutually exclusive with `--v0-pop-g-prior-mu`, `--v0-pop-g-prior-scale`,
`--v0-pop-growth-rate-min`, and `--v0-pop-growth-rate-max`.

### BEAST XML: expressing bounds

Both BEAST2 and BEAST X support `lower` and `upper` attributes on `<parameter>`
elements, which enforce hard bounds on parameters — proposals outside the bounds are
automatically rejected.  So the bounds can be expressed in BEAST XML by setting
`lower="..."` and/or `upper="..."` on the `growthRate` parameter element, without
needing a compound distribution.

### Naming

Follow the established pattern: `pop_g_prior_mu` and `pop_g_prior_scale` for the
Laplace parameters, `pop_g_min` and `pop_g_max` for the bounds.  The prefix `pop_g_`
parallels `pop_inv_n0_` and makes clear these relate to the growth rate g.


## Implementation Plan

### Step 1: `core/run.h` — Add accessors, setters, and private members

Add after `pop_inv_n0_prior_beta` accessors (together with the other parameter
getters/setters):
```cpp
auto pop_g_prior_mu() const -> double { return pop_g_prior_mu_; }
auto set_pop_g_prior_mu(double pop_g_prior_mu) -> void {
  pop_g_prior_mu_ = pop_g_prior_mu, invalidate_derived_quantities(); }
auto pop_g_prior_scale() const -> double { return pop_g_prior_scale_; }
auto set_pop_g_prior_scale(double pop_g_prior_scale) -> void {
  pop_g_prior_scale_ = pop_g_prior_scale, invalidate_derived_quantities(); }
auto pop_g_min() const -> double { return pop_g_min_; }
auto set_pop_g_min(double pop_g_min) -> void {
  pop_g_min_ = pop_g_min, invalidate_derived_quantities(); }
auto pop_g_max() const -> double { return pop_g_max_; }
auto set_pop_g_max(double pop_g_max) -> void {
  pop_g_max_ = pop_g_max, invalidate_derived_quantities(); }
```

Private members (after `pop_inv_n0_prior_beta_`):
```cpp
double pop_g_prior_mu_;
double pop_g_prior_scale_;
double pop_g_min_;
double pop_g_max_;
```

### Step 2: `core/run.cpp` — Constructor, prior calculation, and MH move

**Constructor:** add member initializers:
```cpp
pop_g_prior_mu_{0.001 / 365.0},       // BEAUti2 2.6.2 default, in e-foldings / day
pop_g_prior_scale_{30.701135 / 365.0}, // BEAUti2 2.6.2 default, in e-foldings / day
pop_g_min_{-std::numeric_limits<double>::infinity()},  // No lower bound by default
pop_g_max_{+std::numeric_limits<double>::infinity()},   // No upper bound by default
```

**`calc_cur_log_other_priors()`** (lines 516-520): replace the hard-coded values:
```cpp
// pop_g - Laplace(pop_g_prior_mu_, pop_g_prior_scale_) prior on g in [pop_g_min_, pop_g_max_]
//   pi(g) ~ exp(-|g - mu| / scale) * I(g_min <= g <= g_max)
log_prior += -std::abs(exp_pop_model.growth_rate() - pop_g_prior_mu_) / pop_g_prior_scale_
             - std::log(2 * pop_g_prior_scale_);
```
Note: we do not include a normalization correction for the truncation.
The log-prior is computed up to an additive constant that doesn't depend on g, which
is all that matters for MH ratios.

**`pop_growth_rate_move()`** (lines 1265-1304): replace hard-coded values and add
bound rejection:
```cpp
// Truncated Laplace prior on g with optional bounds [pop_g_min_, pop_g_max_]:
//   pi(g) ~ exp(-|g - mu| / scale) * I(g_min <= g <= g_max)
// By default (following LeMieux et al (2021)), mu = 0.001/365, scale = 30.701135/365,
// and no bounds.
// Random walk operator: uniform perturbation in [-delta, +delta], delta = 1.0 / 365.0.

auto window_size = 1.0 / 365.0;
auto old_g = exp_pop_model->growth_rate();
CHECK(old_g >= pop_g_min_ && old_g <= pop_g_max_)
    << "Growth rate " << old_g << " outside bounds [" << pop_g_min_ << ", " << pop_g_max_ << "]";

auto delta = absl::Uniform(bitgen_, -window_size, +window_size);
auto new_g = old_g + delta;

// Reject if outside bounds
if (new_g < pop_g_min_ || new_g > pop_g_max_) { return; }

auto log_prior_new_over_prior_old =
    (std::abs(old_g - pop_g_prior_mu_) - std::abs(new_g - pop_g_prior_mu_)) / pop_g_prior_scale_;
```

### Step 3: `core/api.fbs` — Add fields to FlatBuffers schema

Add to the `Params` table using next free ids 42-45:
```
// Laplace prior on exponential pop model growth rate g.
// FlatBuffers defaults to 0.0, but actual defaults are nonzero:
//   mu = 0.001/365, scale = 30.701135/365, min = -inf, max = +inf.
// Since pop_g_prior_scale should never be 0 in valid data,
// a value of 0.0 reliably indicates "field not present in old file".
pop_g_prior_mu: float64 (id: 42);    // mu (location) of Laplace prior on g [units: 1/day]
pop_g_prior_scale: float64 (id: 43); // scale of Laplace prior on g [units: 1/day, > 0]
pop_g_min: float64 (id: 44);         // lower bound on g [units: 1/day]
pop_g_max: float64 (id: 45);         // upper bound on g [units: 1/day]
```
Update "Next free id" comment to 46.

### Step 4: `core/api_generated.h` — Regenerate

Run `cd core && ../third-party/flatbuffers/build/flatc --cpp api.fbs`.

### Step 5: `core/api.cpp` — Serialization and deserialization

**Serialization:**
```cpp
params_builder.add_pop_g_prior_mu(run.pop_g_prior_mu());
params_builder.add_pop_g_prior_scale(run.pop_g_prior_scale());
params_builder.add_pop_g_min(run.pop_g_min());
params_builder.add_pop_g_max(run.pop_g_max());
```

**Deserialization:** handle old files where fields default to 0.0.  Use
`pop_g_prior_scale == 0.0` as the sentinel for "fields not present":
```cpp
auto pop_g_prior_scale = api_params->pop_g_prior_scale();
if (pop_g_prior_scale != 0.0) {
  run.set_pop_g_prior_mu(api_params->pop_g_prior_mu());
  run.set_pop_g_prior_scale(pop_g_prior_scale);
  run.set_pop_g_min(api_params->pop_g_min());
  run.set_pop_g_max(api_params->pop_g_max());
} else {
  // Old file without these fields — scale == 0.0 is never a valid value,
  // so this reliably indicates the fields are absent.
  // Use values hard-coded into earlier versions of Delphy.
  run.set_pop_g_prior_mu(0.001 / 365.0);
  run.set_pop_g_prior_scale(30.701135 / 365.0);
  run.set_pop_g_min(-std::numeric_limits<double>::infinity());
  run.set_pop_g_max(+std::numeric_limits<double>::infinity());
}
```

### Step 6: `core/cmdline.cpp` — CLI options

```cpp
("v0-pop-g-prior-mu",
 "[pop-model == exponential] Location (mu) of the Laplace prior on the growth rate g: "
 "pi(g) ~ exp(-|g - mu| / scale).  Units: e-foldings / year.",
 cxxopts::value<double>()->default_value("0.001"))
("v0-pop-g-prior-scale",
 "[pop-model == exponential] Scale of the Laplace prior on the growth rate g: "
 "pi(g) ~ exp(-|g - mu| / scale).  Units: e-foldings / year.",
 cxxopts::value<double>()->default_value("30.701135"))
("v0-pop-growth-rate-min",
 "[pop-model == exponential] Lower bound on the growth rate g, in e-foldings / year.  "
 "E.g., use 0 to constrain the population to be non-declining.  "
 "When g_min = mu = 0, the prior reduces to Exponential(rate = 1/scale).",
 cxxopts::value<double>())
("v0-pop-growth-rate-max",
 "[pop-model == exponential] Upper bound on the growth rate g, in e-foldings / year.",
 cxxopts::value<double>())
("v0-pop-g-prior-exponential-with-mean",
 "[pop-model == exponential] Shorthand for an Exponential prior on the growth rate g "
 "with the specified mean, in e-foldings / year.  "
 "If positive, enforces g >= 0 (growing population); if negative, enforces g <= 0 "
 "(declining population).  Equivalent to setting mu = 0, scale = |mean|, and the "
 "appropriate bound.  "
 "Mutually exclusive with --v0-pop-g-prior-mu, --v0-pop-g-prior-scale, "
 "--v0-pop-growth-rate-min, and --v0-pop-growth-rate-max.",
 cxxopts::value<double>())
```

**Parsing:** Two mutually exclusive interfaces:
- **Interface A** (direct): `--v0-pop-g-prior-mu`, `--v0-pop-g-prior-scale`,
  `--v0-pop-growth-rate-min`, `--v0-pop-growth-rate-max`.
- **Interface B** (shorthand): `--v0-pop-g-prior-exponential-with-mean`.

For both interfaces, convert all CLI values from per-year to per-day units by dividing
by 365: this applies to mu, scale, g_min, g_max, and the exponential mean.  Validate
that scale > 0, that g_min <= g_max (if both specified), and that the initial growth
rate is within [g_min, g_max].

The `Exp_pop_model` is created after all priors are configured, so `init_pop_growth_rate`
can be adjusted before the model is constructed.  If the user explicitly specified any of
the new g prior options (`--v0-pop-g-prior-mu`, `--v0-pop-g-prior-scale`,
`--v0-pop-growth-rate-min`, `--v0-pop-growth-rate-max`, or
`--v0-pop-g-prior-exponential-with-mean`) but did NOT explicitly set
`--v0-init-pop-growth-rate`, override `init_pop_growth_rate` to the mean of the truncated
prior (see closed-form expression below).  This parallels the mu prior's behavior of
defaulting the initial value to the prior mean.  When no new options are specified, the
current default of init_pop_growth_rate = 0 is preserved.  This is only done when new options are
explicitly specified, to preserve the behavior of older versions of Delphy.

Similarly, the n0 prior override now modifies `init_final_pop_size` directly instead of
creating the `Exp_pop_model` eagerly.

Update the `--v0-init-pop-growth-rate` help text to mention the new default behavior.

**Mean of the truncated Laplace on [a, b] with location mu and scale s** (assuming
a <= mu <= b):

```
E[g] = mu + (s/2) * [(p+1)*exp(-p) - (q+1)*exp(-q)] / [1 - (exp(-p) + exp(-q))/2]
```
where `p = (mu - a) / s` and `q = (b - mu) / s`.

Special cases:
- No bounds (a = -inf, b = +inf): E[g] = mu.
- a = mu, b = +inf (exponential starting at mu): E[g] = mu + s.
- a = -inf, b = mu (reflected exponential ending at mu): E[g] = mu - s.

For the `--v0-pop-g-prior-exponential-with-mean` shorthand, the initial g is simply the
specified mean (which equals the distribution mean by construction).

The mean calculation is implemented as a helper function in `cmdline.cpp`:
```cpp
static auto truncated_laplace_mean(double mu, double s, double a, double b) -> double {
  CHECK_GT(s, 0.0);
  CHECK_LE(a, b);
  CHECK_LE(a, mu);
  CHECK_LE(mu, b);

  auto p = (mu - a) / s;  // may be +inf
  auto q = (b - mu) / s;  // may be +inf

  auto result = 0.0;

  // Common special cases (avoid inf arithmetic)
  if (std::isinf(p) && std::isinf(q)) {
    result = mu;                                               // no bounds
  } else if (std::isinf(p)) {                                 // only upper bound
    auto eq = std::exp(-q);
    result = mu + (s / 2) * (-(q + 1) * eq) / (1 - eq / 2);
  } else if (std::isinf(q)) {                                 // only lower bound
    auto ep = std::exp(-p);
    result = mu + (s / 2) * ((p + 1) * ep) / (1 - ep / 2);
  } else if (p + q < 1e-4) {
    result = (a + b) / 2;                                      // Taylor fallback for tight bounds
  } else {
    auto ep = std::exp(-p);
    auto eq = std::exp(-q);
    result = mu + (s / 2) * ((1 + p) * ep - (1 + q) * eq) / (1 - (ep + eq) / 2);
  }

  CHECK_GE(result, a);
  CHECK_LE(result, b);
  return result;
}
```
The general formula is `E[g] = mu + (s/2) * [(p+1)exp(-p) - (q+1)exp(-q)] / [1 - (exp(-p) + exp(-q))/2]`
where `p = (mu - a)/s` and `q = (b - mu)/s`.  Special cases handle infinite bounds
(avoiding inf arithmetic) and the degenerate case where both bounds are close to mu
relative to s (where the 0/0 limit is `(a+b)/2`, the midpoint).

Add the new option names to the list of exponential-only options checked against wrong
pop model.

### Step 7: `tools/delphy_wasm.cpp` — WASM getter/setter pairs

Add eight `EMSCRIPTEN_KEEPALIVE extern "C"` functions near the exp pop model getters:
```cpp
auto delphy_run_get_pop_g_prior_mu(Delphy_context& /*ctx*/, const Run& run) -> double;
auto delphy_run_set_pop_g_prior_mu(Delphy_context& /*ctx*/, Run& run, double mu) -> void;
auto delphy_run_get_pop_g_prior_scale(Delphy_context& /*ctx*/, const Run& run) -> double;
auto delphy_run_set_pop_g_prior_scale(Delphy_context& /*ctx*/, Run& run, double scale) -> void;
auto delphy_run_get_pop_g_min(Delphy_context& /*ctx*/, const Run& run) -> double;
auto delphy_run_set_pop_g_min(Delphy_context& /*ctx*/, Run& run, double g_min) -> void;
auto delphy_run_get_pop_g_max(Delphy_context& /*ctx*/, const Run& run) -> double;
auto delphy_run_set_pop_g_max(Delphy_context& /*ctx*/, Run& run, double g_max) -> void;
```

### Step 7b: `~/now/delphy-web/src/ts/pythia/delphy_api.ts` — TypeScript API bindings

Add corresponding entries near the exp pop model declarations in:
1. `DelphyCoreRaw` interface
2. Module bindings
3. `DelphyRun` class methods

### Step 8: `doc/dphy_file_format.md` — Update documentation

Add the four new fields to the `Params` table documentation.

### Step 9: `core/beasty_input.cpp` — BEAST XML export

Update all three BEAST versions to use the configurable mu and scale values
instead of hard-coded ones.  Convert internal per-day values back to per-year
for XML output by multiplying by 365.  For the bounds, set `lower` and/or `upper`
attributes on the `growthRate` `<parameter>` element when non-default bounds are
specified:

**BEAST 2 (2.6.2 and 2.7.7):**
```xml
<parameter id="growthRate.t:input_alignment" ... lower="g_min_years" upper="g_max_years" ...>
```

**BEAST X (10.5.0):**
```xml
<parameter id="exponential.growthRate" value="..." lower="g_min_years" upper="g_max_years"/>
```

Only emit `lower`/`upper` when the bound is finite.

### Step 10: Build and test

**10a.** Build debug, run unit tests.

**10b.** Short MCMC run with default priors (should match current behavior):
```bash
./build/debug/delphy \
  --v0-in-fasta ~/now/delphy-2026-paper-data/ebola-gire-2014/delphy_inputs/ebola.fasta \
  --v0-steps 0
```

**10c.** Short MCMC run with non-default priors:
```bash
./build/debug/delphy \
  --v0-in-fasta ~/now/delphy-2026-paper-data/ebola-gire-2014/delphy_inputs/ebola.fasta \
  --v0-steps 100 \
  --v0-pop-g-prior-mu 0.01 --v0-pop-g-prior-scale 1.0 --v0-pop-growth-rate-min 0
```

**10d.** Generate BEAST XML in all three formats with non-default priors and verify
that BEAST can parse them.  Runnable BEAST executables are at:
- `~/github/CompEvol/beast2.6.2/bin/beast` (BEAST 2.6.2)
- `~/now/delphy-2026-paper-data/beast277` (BEAST 2.7.7)
- `~/now/delphy-2026-paper-data/beastX1050` (BEAST X 10.5.0)
