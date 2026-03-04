# Flexible Mu Prior

## Existing Pattern: `skygrid_tau_prior_alpha` / `skygrid_tau_prior_beta`

These parameters define a Gamma(alpha, beta) prior on the Skygrid precision parameter `tau`,
with density `p(tau) ~ tau^{alpha-1} * exp(-beta * tau)`.  They appear in 9 files across the
codebase.  Here is where and why each one matters.

### 1. `core/run.h` — Accessors, setters, and private members

**Accessors/setters** (lines 61-64):
```cpp
auto skygrid_tau_prior_alpha() const -> double { return skygrid_tau_prior_alpha_; }
auto set_skygrid_tau_prior_alpha(double skygrid_tau_prior_alpha) -> void {
  skygrid_tau_prior_alpha_ = skygrid_tau_prior_alpha, invalidate_derived_quantities(); }
auto skygrid_tau_prior_beta() const -> double { return skygrid_tau_prior_beta_; }
auto set_skygrid_tau_prior_beta(double skygrid_tau_prior_beta) -> void {
  skygrid_tau_prior_beta_ = skygrid_tau_prior_beta, invalidate_derived_quantities(); }
```

**Private members** (lines 227-228):
```cpp
double skygrid_tau_prior_alpha_;
double skygrid_tau_prior_beta_;
```

**Why:** Central storage and access point.  Setters call `invalidate_derived_quantities()` so
that cached log-prior values are recomputed on next access.

### 2. `core/run.cpp` — Default values, prior calculation, and Gibbs move

**Constructor defaults** (lines 21-22):
```cpp
skygrid_tau_prior_alpha_{0.001},  // Prior, Gill et al 2012 discussion below Eq. 16
skygrid_tau_prior_beta_{0.001},
```

**Prior calculation in `calc_cur_log_other_priors()`** (lines 519-524):
```cpp
if (skygrid_tau_move_enabled_) {
  // tau - gamma prior (Gill et al 2012, Eq. 15)
  log_prior +=
      (skygrid_tau_prior_alpha_ - 1) * log_tau
      - skygrid_tau_prior_beta_ * tau;
}
```
This adds `log Gamma(tau; alpha, beta)` (up to a constant) to the total prior.  It is only
included when the tau move is enabled (i.e., tau is being inferred, not fixed).

**Gibbs move in `skygrid_tau_move()`** (lines 1289-1326):
The posterior for tau is conjugate Gamma:
```
posterior: Gamma(alpha + M/2, beta + sum_squared_delta_gammas / 2)
```
The code Gibbs-samples directly from this posterior:
```cpp
auto post_alpha = skygrid_tau_prior_alpha_ + 0.5 * M;
auto post_beta = skygrid_tau_prior_beta_ + 0.5 * sum_squared_delta_gammas;
auto tau_dist = std::gamma_distribution<double>{post_alpha, 1.0 / post_beta};
auto new_tau = tau_dist(bitgen_);
```
Then it incrementally updates `log_other_priors_`:
```cpp
log_other_priors_ +=
    (post_alpha - 1) * std::log(new_tau/old_tau)
    - post_beta * (new_tau - old_tau);
```
Note: the incremental update uses `post_alpha` and `post_beta` (the full conditional
parameters), not the prior parameters.  This is because `log_other_priors_` includes the
GMRF prior on the gammas, which also depends on tau.

### 3. `core/cmdline.cpp` — CLI option definitions and parsing

**Option definitions** (lines 199-204):
```cpp
("v0-skygrid-tau-prior-alpha",
 "[pop-model == skygrid] When inferring the log-population curve's prior smoothness `tau`,
  the hyperprior on tau is ~ tau^{alpha - 1} exp[-beta tau]",
 cxxopts::value<double>()->default_value("0.001"))
("v0-skygrid-tau-prior-beta",
 "[pop-model == skygrid] ...",
 cxxopts::value<double>()->default_value("0.001"))
```

**Parsing** (lines 511-522): when `--v0-skygrid-infer-prior-smoothness` is true, the alpha
and beta values are read and validated (must be positive), then passed to the Run via setters:
```cpp
auto prior_alpha = opts["v0-skygrid-tau-prior-alpha"].as<double>();
auto prior_beta = opts["v0-skygrid-tau-prior-beta"].as<double>();
if (prior_alpha <= 0.0 || prior_beta <= 0.0) {
  std::cerr << "ERROR: Skygrid tau prior parameters must be positive\n";
  std::exit(EXIT_FAILURE);
}
run->set_skygrid_tau_prior_alpha(prior_alpha);
run->set_skygrid_tau_prior_beta(prior_beta);
skygrid_tau = prior_alpha / prior_beta;  // Prior mean as initial value
```

**Detection** (lines 438-439): used in the check that determines whether any skygrid-related
CLI options were passed (to error if used with the exponential pop model).

### 4. `tools/delphy_wasm.cpp` — WASM/JavaScript interface

**Getter/setter pairs** (lines 946-966):
```cpp
auto delphy_run_get_skygrid_tau_prior_alpha(Delphy_context& /*ctx*/, const Run& run) -> double;
auto delphy_run_set_skygrid_tau_prior_alpha(Delphy_context& /*ctx*/, Run& run, double v) -> void;
auto delphy_run_get_skygrid_tau_prior_beta(Delphy_context& /*ctx*/, const Run& run) -> double;
auto delphy_run_set_skygrid_tau_prior_beta(Delphy_context& /*ctx*/, Run& run, double v) -> void;
```

**Why:** Allows the web-based UI (compiled via Emscripten) to get/set these parameters at
runtime.  Each is an `EMSCRIPTEN_KEEPALIVE extern "C"` function.

### 5. `core/api.fbs` — FlatBuffers schema

**Fields in `Params` table** (lines 110-111):
```
skygrid_tau_prior_alpha: float64 (id: 32);
skygrid_tau_prior_beta: float64 (id: 33);
```

**Why:** These are persisted in `.dphy` files (Delphy's native save format) so that a
resumed run uses the same prior parameters.  The `id:` tags ensure backward/forward
compatibility.  Next free id is 38 (line 95).

### 6. `core/api_generated.h` — Generated FlatBuffers code

Auto-generated getters, setters, and builder methods for the `Params` table.  Default value
is 0.0 in the flatbuffer (the Run constructor applies its own default of 0.001).

**Why:** This file is regenerated from `api.fbs`; no manual editing needed.

### 7. `core/api.cpp` — Serialization and deserialization

**Serialization** (lines 265-266):
```cpp
params_builder.add_skygrid_tau_prior_alpha(run.skygrid_tau_prior_alpha());
params_builder.add_skygrid_tau_prior_beta(run.skygrid_tau_prior_beta());
```

**Deserialization** (lines 365-366):
```cpp
run.set_skygrid_tau_prior_alpha(api_params->skygrid_tau_prior_alpha());
run.set_skygrid_tau_prior_beta(api_params->skygrid_tau_prior_beta());
```

**Why:** Marshals the prior parameters between the `Run` object and its FlatBuffers
serialized form (for save/load of `.dphy` files and for the WASM API boundary).

### 8. `core/beasty_input.cpp` — BEAST XML export

**Export** (lines 1238-1240):
```cpp
os << absl::StreamFormat(
    "  <gammaPrior id=\"skygrid.precision.prior\" shape=\"%g\" scale=\"%g\" offset=\"0.0\">\n",
    run.skygrid_tau_prior_alpha(),
    1.0 / run.skygrid_tau_prior_beta())
```

**Why:** When exporting to BEAST XML format, the Gamma prior is converted from Delphy's
(alpha, rate=beta) parameterization to BEAST's (shape, scale=1/beta) parameterization.

### 9. `doc/dphy_file_format.md` — Documentation

Lines 162-163: mirrors the `api.fbs` schema definition with the same comments.

**Why:** Documents the `.dphy` file format for external consumers.


## Current Mu Prior Implementation

The current code in `run.cpp` implements a **uniform prior** on mu.

**Prior calculation in `calc_cur_log_other_priors()`** (lines 469-473):
```cpp
if (not mpox_hack_enabled_) {
  // Mu - Uniform prior
} else {
  // Mu & Mu* - Uniform priors
}
```
No term is added to `log_prior` for mu — the prior contributes exactly 0.

**Gibbs move in `mu_move()`** (lines 758-796):
```cpp
// We follow LeMieux et al (2021) and use a uniform prior for mu.
// We can Gibbs sample a new value of mu using a Gamma with
// shape alpha = M+1 and rate beta = Ttwiddle
auto mu_dist = std::gamma_distribution<double>{
    num_muts_ + 1.0,
    1.0 / Ttwiddle};
// ...
log_G_ += -(new_mu - old_mu) * Ttwiddle + num_muts_ * std::log(new_mu / old_mu);
log_other_priors_ += 0.0;
```

The Gibbs sampler shape is `M + 1` (where M = num_muts_), consistent with a uniform prior.
With a `Gamma(alpha, beta)` prior on mu where `p(mu) ~ mu^{alpha-1} exp(-beta*mu)`:
- The full conditional posterior is `Gamma(M + alpha, Ttwiddle + beta)`
- Defaults of `alpha = 1, beta = 0` reproduce the current uniform prior (shape = M + 1,
  rate = Ttwiddle), since `Gamma(1, 0)` gives `p(mu) ~ mu^0 * exp(0) = 1` (uniform).


## Implementation Plan

Add `mu_prior_alpha` and `mu_prior_beta` parameters following the same pattern as the
skygrid tau prior parameters.  The prior on mu is `p(mu) ~ mu^{alpha-1} exp(-beta*mu)`.
Defaults: `alpha = 1, beta = 0` (uniform prior, matching current behavior).

### Step 1: `core/run.h` — Add accessors, setters, and private members

Add after the `mu()`/`set_mu()` accessors (line 45):
```cpp
auto mu_prior_alpha() const -> double { return mu_prior_alpha_; }
auto set_mu_prior_alpha(double mu_prior_alpha) -> void {
  mu_prior_alpha_ = mu_prior_alpha, invalidate_derived_quantities(); }
auto mu_prior_beta() const -> double { return mu_prior_beta_; }
auto set_mu_prior_beta(double mu_prior_beta) -> void {
  mu_prior_beta_ = mu_prior_beta, invalidate_derived_quantities(); }
```

Add private members after `hky_model_` (line 233):
```cpp
double mu_prior_alpha_;
double mu_prior_beta_;
```

### Step 2: `core/run.cpp` — Constructor, prior calculation, and Gibbs moves

**Constructor** (after `nu_` to match member declaration order): add member initializers:
```cpp
mu_prior_alpha_{1.0},   // Default: uniform prior on mu
mu_prior_beta_{0.0},
```

**`calc_cur_log_other_priors()`** (lines 469-473): replace the empty uniform-prior block:
```cpp
if (not mpox_hack_enabled_) {
  // Mu - Gamma(mu_prior_alpha_, mu_prior_beta_) prior
  log_prior += (mu_prior_alpha_ - 1) * std::log(mu()) - mu_prior_beta_ * mu();
} else {
  // Mu - Gamma(mu_prior_alpha_, mu_prior_beta_) prior (mu* has a uniform prior)
  log_prior += (mu_prior_alpha_ - 1) * std::log(mpox_mu_) - mu_prior_beta_ * mpox_mu_;
}
```
With the default `alpha=1, beta=0`, both terms are 0, matching current behavior.

**`mu_move()`** (lines 758-796): update Gibbs sampler to incorporate prior parameters.
The full conditional posterior is `Gamma(M + alpha, Ttwiddle + beta)`:
```cpp
auto mu_dist = std::gamma_distribution<double>{
    num_muts_ + mu_prior_alpha_,
    1.0 / (Ttwiddle + mu_prior_beta_)};
```
Update comments accordingly.  Update `log_other_priors_` incrementally:
```cpp
log_other_priors_ += (mu_prior_alpha_ - 1) * std::log(new_mu / old_mu)
                     - mu_prior_beta_ * (new_mu - old_mu);
```

**`mpox_hack_moves()`** (lines 878-895): analogous change.  The mpox hack changes variables
from `(mu, mu*)` to `(mu, rho = mu*/mu)`, introducing a Jacobian of `1/mu` that shifts the
shape by -1.  The full conditional is `Gamma(M + alpha - 1, Ttwiddle_eff + beta)`:
```cpp
auto mu_dist = std::gamma_distribution<double>{
    M + mu_prior_alpha_ - 1.0,
    1.0 / (Ttwiddle_eff + mu_prior_beta_)};
```
And the same `log_other_priors_` incremental update as the non-mpox case:
```cpp
log_other_priors_ += (mu_prior_alpha_ - 1) * std::log(new_mu / old_mu)
                     - mu_prior_beta_ * (new_mu - old_mu);
```

### Step 3: `core/api.fbs` — Add fields to FlatBuffers schema

Add to the `Params` table, using next free ids 38 and 39.  Set the flatbuffer default for
`mu_prior_alpha` to `1.0` so that old `.dphy` files (which lack this field) deserialize to
the current uniform-prior behavior.  `mu_prior_beta` defaults to `0.0` (implicit flatbuffer
default), which is also correct.
```
mu_prior_alpha: float64 = 1.0 (id: 38);  // alpha parameter of mu prior [unitless, >= 0]
mu_prior_beta: float64 (id: 39);         // beta parameter of mu prior [units: days, >= 0]
```
Update the "Next free id" comment to 40.

### Step 4: `core/api_generated.h` — Regenerate

Run `cd core && ../third-party/flatbuffers/build/flatc --cpp api.fbs` to regenerate.

### Step 5: `core/api.cpp` — Serialization and deserialization

**Serialization** (after the `mu` line):
```cpp
params_builder.add_mu_prior_alpha(run.mu_prior_alpha());
params_builder.add_mu_prior_beta(run.mu_prior_beta());
```

**Deserialization** (after the `set_mu` line):
```cpp
run.set_mu_prior_alpha(api_params->mu_prior_alpha());
run.set_mu_prior_beta(api_params->mu_prior_beta());
```

### Step 6: `core/cmdline.cpp` — CLI options

Two alternative interfaces, mutually exclusive with each other.

**Interface A: shape/rate (alpha/beta)**

Option definitions (near `v0-fix-mutation-rate` / `v0-init-mutation-rate`):
```cpp
("v0-mu-prior-alpha",
 "Shape (alpha) parameter of the Gamma prior on the mutation rate mu: "
 "p(mu) ~ mu^{alpha-1} exp[-beta mu].  Default alpha=1, beta=0 gives a uniform prior.",
 cxxopts::value<double>()->default_value("1.0"))
("v0-mu-prior-beta",
 "Rate (beta) parameter of the Gamma prior on the mutation rate mu: "
 "p(mu) ~ mu^{alpha-1} exp[-beta mu].  Default alpha=1, beta=0 gives a uniform prior.  "
 "Units: years.",
 cxxopts::value<double>()->default_value("0.0"))
```

Note: the CLI presents beta in per-year units (matching `--v0-init-mutation-rate`).
Internally, `beta_internal = beta_cli * 365` (since Delphy's mu is in per-day units).

**Interface B: mean/stddev**

```cpp
("v0-mu-prior-mean",
 "Mean of the Gamma prior on the mutation rate mu, in subst / site / year.  "
 "Must be specified together with --v0-mu-prior-stddev.  "
 "Mutually exclusive with --v0-mu-prior-alpha / --v0-mu-prior-beta.",
 cxxopts::value<double>())
("v0-mu-prior-stddev",
 "Standard deviation of the Gamma prior on the mutation rate mu, in subst / site / year.  "
 "Must be specified together with --v0-mu-prior-mean.  "
 "Mutually exclusive with --v0-mu-prior-alpha / --v0-mu-prior-beta.",
 cxxopts::value<double>())
```

Conversion from mean/stddev (in per-year units) to internal alpha/beta (in per-day units):
```
alpha = (mean / stddev)^2                         // dimensionless
beta_internal = mean * 365 / stddev^2             // units: days
```

(Derivation: for Gamma(alpha, beta), mean = alpha/beta, var = alpha/beta^2.  So
alpha = mean^2/var, beta = mean/var.  Since beta has units of time (it multiplies mu,
which has units of 1/time), convert from years to days by multiplying by 365.)

**Parsing:**

```cpp
auto has_mu_prior_alpha_beta =
    (opts.count("v0-mu-prior-alpha") > 0 || opts.count("v0-mu-prior-beta") > 0);
auto has_mu_prior_mean_stddev =
    (opts.count("v0-mu-prior-mean") > 0 || opts.count("v0-mu-prior-stddev") > 0);

if (has_mu_prior_alpha_beta && has_mu_prior_mean_stddev) {
  std::cerr << "ERROR: --v0-mu-prior-alpha/beta and --v0-mu-prior-mean/stddev "
            << "are mutually exclusive\n";
  std::exit(EXIT_FAILURE);
}

if (has_mu_prior_mean_stddev) {
  // Both mean and stddev must be specified together
  if (opts.count("v0-mu-prior-mean") == 0 || opts.count("v0-mu-prior-stddev") == 0) {
    std::cerr << "ERROR: --v0-mu-prior-mean and --v0-mu-prior-stddev "
              << "must be specified together\n";
    std::exit(EXIT_FAILURE);
  }
  auto mean = opts["v0-mu-prior-mean"].as<double>();   // per year
  auto stddev = opts["v0-mu-prior-stddev"].as<double>(); // per year
  // validate mean > 0, stddev > 0
  auto alpha = std::pow(mean / stddev, 2);
  auto beta_internal = mean * 365.0 / (stddev * stddev);
  run->set_mu_prior_alpha(alpha);
  run->set_mu_prior_beta(beta_internal);
} else {
  auto alpha = opts["v0-mu-prior-alpha"].as<double>();
  auto beta_cli = opts["v0-mu-prior-beta"].as<double>();
  // validate alpha >= 0, beta_cli >= 0
  run->set_mu_prior_alpha(alpha);
  run->set_mu_prior_beta(beta_cli * 365.0);  // convert per-year to per-day
}
```

Note: unlike the skygrid tau prior (which requires strictly positive alpha and beta), the mu
prior allows `alpha >= 0` and `beta >= 0`.  `alpha = 0` gives the Jeffreys prior `1/mu`
(proper as long as there is at least one mutation).  `beta = 0` gives an improper prior that
is nevertheless fine in practice (the likelihood provides the exponential decay).

### Step 7: `tools/delphy_wasm.cpp` — WASM getter/setter pairs

Add four `EMSCRIPTEN_KEEPALIVE extern "C"` functions following the existing pattern:
```cpp
auto delphy_run_get_mu_prior_alpha(Delphy_context& /*ctx*/, const Run& run) -> double;
auto delphy_run_set_mu_prior_alpha(Delphy_context& /*ctx*/, Run& run, double mu_prior_alpha) -> void;
auto delphy_run_get_mu_prior_beta(Delphy_context& /*ctx*/, const Run& run) -> double;
auto delphy_run_set_mu_prior_beta(Delphy_context& /*ctx*/, Run& run, double mu_prior_beta) -> void;
```

### Step 7b: `~/now/delphy-web/src/ts/pythia/delphy_api.ts` — TypeScript API bindings

The Delphy Web project (`~/now/delphy-web`) wraps the WASM API in TypeScript.  The skygrid
tau prior parameters appear in three places in `delphy_api.ts`:

**1. `DelphyCoreRaw` interface** (type declarations, near line 628):
```typescript
run_get_mu_prior_alpha: (ctx: DelphyContextPtr, run: RunPtr) => number,
run_set_mu_prior_alpha: (ctx: DelphyContextPtr, run: RunPtr, mu_prior_alpha: number) => void,
run_get_mu_prior_beta: (ctx: DelphyContextPtr, run: RunPtr) => number,
run_set_mu_prior_beta: (ctx: DelphyContextPtr, run: RunPtr, mu_prior_beta: number) => void,
```

**2. Module bindings** (mapping to WASM exports, near line 872):
```typescript
run_get_mu_prior_alpha: Module['_delphy_run_get_mu_prior_alpha'],
run_set_mu_prior_alpha: Module['_delphy_run_set_mu_prior_alpha'],
run_get_mu_prior_beta: Module['_delphy_run_get_mu_prior_beta'],
run_set_mu_prior_beta: Module['_delphy_run_set_mu_prior_beta'],
```

**3. `DelphyRun` class methods** (typed getter/setter wrappers, near line 1477):
```typescript
getMuPriorAlpha(): number {
  return Delphy.delphyCoreRaw.run_get_mu_prior_alpha(this.delphy.ctx, this.run);
}

setMuPriorAlpha(mu_prior_alpha: number): void {
  Delphy.delphyCoreRaw.run_set_mu_prior_alpha(this.delphy.ctx, this.run, mu_prior_alpha);
}

getMuPriorBeta(): number {
  return Delphy.delphyCoreRaw.run_get_mu_prior_beta(this.delphy.ctx, this.run);
}

setMuPriorBeta(mu_prior_beta: number): void {
  Delphy.delphyCoreRaw.run_set_mu_prior_beta(this.delphy.ctx, this.run, mu_prior_beta);
}
```

### Step 8: `doc/dphy_file_format.md` — Update documentation

Add the two new fields to the `Params` table documentation, mirroring the `api.fbs` comments.

### Step 9: `core/beasty_input.cpp` — BEAST XML clock rate prior

All three BEAST versions currently emit a `Uniform` prior on the clock rate (or no explicit
prior).  When `mu_prior_alpha != 1` or `mu_prior_beta != 0` (i.e., a non-uniform prior),
replace the Uniform with a Gamma prior.  When the prior is uniform, keep the existing output.

**Unit conversion:** Delphy's mu is in subst/site/day; BEAST's clock.rate is in
subst/site/year (factor of 365).  Since beta has units of time (days internally),
`beta_year = beta_day / 365`.  BEAST 2 uses shape/rate (with `mode="ShapeRate"`),
so rate = `beta_day / 365`.  BEAST X uses shape/scale, so scale = `365 / beta_day`.

**BEAST 2 (2.6.2 and 2.7.7)** — Replace the `<Uniform>` prior block (lines 328-331 and
699-703) with:
```cpp
if (run.mu_move_enabled()) {
  if (run.mu_prior_alpha() != 1.0 || run.mu_prior_beta() != 0.0) {
    // Gamma prior on clock rate (in per-year units)
    os << absl::StreamFormat("      <prior id=\"ClockPrior.c:input_alignment\" name=\"distribution\" x=\"@clockRate.c:input_alignment\">\n")
       << absl::StreamFormat("        <Gamma id=\"Gamma.clockRate\" name=\"distr\" mode=\"ShapeRate\">\n")
       << absl::StreamFormat("          <parameter id=\"RealParameter.clockRate.alpha\" spec=\"parameter.RealParameter\" estimate=\"false\" name=\"alpha\">%g</parameter>\n",
                             run.mu_prior_alpha())
       << absl::StreamFormat("          <parameter id=\"RealParameter.clockRate.beta\" spec=\"parameter.RealParameter\" estimate=\"false\" name=\"beta\">%g</parameter>\n",
                             run.mu_prior_beta() / 365.0)
       << absl::StreamFormat("        </Gamma>\n")
       << absl::StreamFormat("      </prior>\n");
  } else {
    // Uniform prior (current behavior)
    os << absl::StreamFormat("      <prior id=\"ClockPrior.c:input_alignment\" name=\"distribution\" x=\"@clockRate.c:input_alignment\">\n")
       << absl::StreamFormat("        <Uniform id=\"Uniform.0\" name=\"distr\" upper=\"Infinity\"/>\n")
       << absl::StreamFormat("      </prior>\n");
  }
}
```

**BEAST X (10.5.0)** — Add a `<gammaPrior>` in the prior block (near line 1474, where it
currently says "No CTMC prior in Delphy"):
```cpp
if (run.mu_move_enabled() && (run.mu_prior_alpha() != 1.0 || run.mu_prior_beta() != 0.0)) {
  os << absl::StreamFormat("        <gammaPrior id=\"clock.rate.prior\" shape=\"%g\" scale=\"%g\" offset=\"0.0\">\n",
                           run.mu_prior_alpha(),
                           365.0 / run.mu_prior_beta())
     << absl::StreamFormat("          <parameter idref=\"clock.rate\"/>\n")
     << absl::StreamFormat("        </gammaPrior>\n");
}
```

### Step 10: Build and test

**10a. Build and run unit tests.**

Build debug and run `./build/debug/tests/tests` to check for compilation errors and basic
regressions.

**10b. Run a short real MCMC run in debug mode.**

Unit tests do not exercise a real MCMC run.  Instead, perform a short real run using the
ebola dataset (`~/now/delphy-2026-paper-data/ebola-gire-2014/delphy_inputs/ebola.fasta`),
adapted from the example in `02a_run_delphy.sh`.  Run with the debug build for a small
number of steps (e.g., 10000), without persisting log/tree/dphy outputs:

```bash
./build/debug/tools/delphy \
  --v0-in-fasta ~/now/delphy-2026-paper-data/ebola-gire-2014/delphy_inputs/ebola.fasta \
  --v0-steps 10000
```

The `check_derived_quantities()` assertions that run after every MCMC move will catch any
inconsistency between the incrementally-maintained `log_other_priors_` and the value
recomputed from scratch by `calc_cur_log_other_priors()`.

Also test with a non-default prior (e.g., `--v0-mu-prior-mean 1e-3 --v0-mu-prior-stddev 5e-4`)
to exercise the new code paths.

**10c. Generate BEAST XML in all three formats and verify.**

```bash
./build/debug/tools/delphy --dry-run \
  --v0-in-fasta ~/now/delphy-2026-paper-data/ebola-gire-2014/delphy_inputs/ebola.fasta \
  --v0-mu-prior-mean 1e-3 --v0-mu-prior-stddev 5e-4 \
  --v0-out-beast-xml /tmp/test_beast2_262.xml --v0-out-beast-version 2.6.2

./build/debug/tools/delphy --dry-run \
  --v0-in-fasta ~/now/delphy-2026-paper-data/ebola-gire-2014/delphy_inputs/ebola.fasta \
  --v0-mu-prior-mean 1e-3 --v0-mu-prior-stddev 5e-4 \
  --v0-out-beast-xml /tmp/test_beast2_277.xml --v0-out-beast-version 2.7.7

./build/debug/tools/delphy --dry-run \
  --v0-in-fasta ~/now/delphy-2026-paper-data/ebola-gire-2014/delphy_inputs/ebola.fasta \
  --v0-mu-prior-mean 1e-3 --v0-mu-prior-stddev 5e-4 \
  --v0-out-beast-xml /tmp/test_beastX.xml --v0-out-beast-version X-10.5.0
```

Inspect each output to verify that the `<Gamma>` / `<gammaPrior>` elements reference the
correct clock rate parameter with properly converted shape and scale/rate values.  Also
verify that the default (uniform) case still produces the existing `<Uniform>` output.