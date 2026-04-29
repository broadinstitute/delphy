# Exponential Population Model Minimum Population Barrier

## Summary

Add a configurable minimum population barrier to the exponential population model to prevent
effective population sizes from dropping below a minimum value (default: 1 day). This mirrors the
existing Skygrid model's low-population barrier and prevents numerical issues in deep trees.

The new population curve is: **N(t) = max(min_pop, n_0 * exp(g * (t - t_0)))**
where t_0 is the reference time (typically the time of the most recent tip) and min_pop is the
minimum effective population size. Note that Delphy currently never has t > t_0 during a run,
but the Exp_pop_model code should not make that assumption.

A min_pop value of 0.0 effectively disables the barrier (since effective population is always positive).

Note that this is a change from the default behaviour, but is almost always an improvement.  
The possibility exists to disable the barrier when needed by setting min_pop to 0.

## Problem Statement

### Current Behavior

The exponential population model implements `N(t) = n_0 * exp(g * (t - t_0))`. In deep trees
stretching far into the past:

- With positive growth rates (g > 0), the effective population can become arbitrarily small as t -> -Infinity
- Tiny effective population sizes (< 1 day) can cause numerical issues leading to crashes
- Effective population sizes below 1 day are biologically implausible given typical sampling date
  uncertainty (~1 day) and viral generation times (several days)

### Skygrid Precedent

The Skygrid model already implements a low-population barrier via a quadratic penalty in the
prior (controlled by `skygrid_low_gamma_barrier_loc_` and `skygrid_low_gamma_barrier_scale_`).
This uses a **soft** barrier (penalizes but doesn't prevent). For the exponential model, we want
a **hard** barrier (strict lower bound) which is simpler and more appropriate for a parametric
model.

## Design

### Barrier Mechanism

**Hard barrier**: N(t) = max(min_pop, n_0 * exp(g * (t - t_0)))

- min_pop = 0.0: barrier disabled (natural since n_0*exp(...) > 0 always)
- min_pop > 0: effective population cannot go below this value
- Default: min_pop = 1.0 day

**Interaction with MCMC moves on n_0:** The barrier acts purely on the population model
used by the coalescent likelihood. The prior on n_0 (Inverse-Gamma) still sees the raw n_0
value. This means the MCMC can propose and accept n_0 < min_pop — when it does, the
coalescent likelihood sees N(t) = min_pop everywhere (for g >= 0), making n_0 unidentifiable
from the coalescent alone. In practice this is harmless: the prior pulls n_0 upward, and the
barrier only activates for deep-past times where the data have no resolving power anyway.

### Crossing Point Analysis

The effective population crosses the minimum when:

```
n_0 * exp(g * (t_c - t_0)) = min_pop

=> t_c = t_0 + log(min_pop / n_0) / g
```

When g = 0, the population is constant and there is no crossing point (handled separately).
When |g| is very small but nonzero, t_c = t_0 + log(min_pop/n_0)/g may be very large in
magnitude; this is mathematically correct and causes no numerical issues since t_c is only
used for interval comparisons, not in exponential arguments.

The crossing point t_c can be precalculated in the constructor since the parameters are immutable
after construction.

**Cases for g != 0:**

| g sign | n_0 vs min_pop | t_c location | Unclamped behavior |
|--------|----------------|--------------|---------------------|
| g > 0  | n_0 > min_pop  | t_c < t_0 (past) | Below min_pop for t < t_c, above for t > t_c |
| g > 0  | n_0 < min_pop  | t_c > t_0 (future) | Below min_pop for t < t_c, above for t > t_c |
| g < 0  | n_0 > min_pop  | t_c > t_0 (future) | Above min_pop for t < t_c, below for t > t_c |
| g < 0  | n_0 < min_pop  | t_c < t_0 (past) | Above min_pop for t < t_c, below for t > t_c |

Note: When n_0 = min_pop, then t_c = t_0, and the effective population equals min_pop at t_0.

## Implementation Plan

### 1. `core/pop_model.h` - Exp_pop_model Class Interface

Add a new member for the minimum population:

```cpp
class Exp_pop_model : public Pop_model {
 public:
  Exp_pop_model(double t0, double pop_at_t0, double growth_rate, double min_pop);

  // Accessors
  auto t0() const -> double { return t0_; }
  auto pop_at_t0() const -> double { return pop_at_t0_; }
  auto growth_rate() const -> double { return growth_rate_; }
  auto min_pop() const -> double { return min_pop_; }

  // Virtual methods - implementations updated
  auto pop_at_time(double t) const -> double override;
  auto pop_integral(double a, double b) const -> double override;
  auto intensity_integral(double a, double b) const -> double override;

 private:
  double t0_;
  double pop_at_t0_;
  double growth_rate_;
  double min_pop_;
  double t_c_;  // time where n0*exp(g*(t-t0)) == min_pop; NaN if no crossing

  auto print_to(std::ostream& os) const -> void override {
    os << absl::StreamFormat("Exp_pop_model{t0=%g, n0=%g, g=%g, min_pop=%g}",
                             t0(), pop_at_t0(), growth_rate(), min_pop());
  }
};
```

### 2. `core/pop_model.cpp` - Exp_pop_model Method Implementations

#### Constructor

```cpp
Exp_pop_model::Exp_pop_model(double t0, double pop_at_t0, double growth_rate, double min_pop)
    : t0_{t0}, pop_at_t0_{pop_at_t0}, growth_rate_{growth_rate}, min_pop_{min_pop} {
  if (pop_at_t0 <= 0.0) {
    throw std::invalid_argument(absl::StrFormat(
        "Initial effective population size should be positive (not %f)", pop_at_t0));
  }
  if (min_pop < 0.0) {
    throw std::invalid_argument(absl::StrFormat(
        "Minimum effective population size should be non-negative (not %f)", min_pop));
  }
  // Precompute crossing point
  if (min_pop > 0.0 && growth_rate != 0.0) {
    t_c_ = t0 + std::log(min_pop / pop_at_t0) / growth_rate;
  } else {
    // min_pop == 0 or growth_rate == 0: no crossing (barrier disabled or constant population)
    t_c_ = std::numeric_limits<double>::quiet_NaN();
  }
}
```

#### `pop_at_time(t)`

```cpp
auto Exp_pop_model::pop_at_time(double t) const -> double {
  return std::max(min_pop_, pop_at_t0_ * std::exp((t - t0_) * growth_rate_));
}
```

#### `pop_integral(a, b)` and `intensity_integral(a, b)`

Both integrals split into clamped and unclamped regions at `t_c`, with the same
clamped region (where `n0*exp(g*(t-t0)) < min_pop`). They differ in what they integrate:

| | Clamped region | Unclamped region |
|---|---|---|
| `pop_integral` | `min_pop` | `n0 * exp(g*(t-t0))` |
| `intensity_integral` | `1/min_pop` | `(1/n0) * exp(-g*(t-t0))` |

Note: a tempting shortcut would be to factor out a shared helper and derive
`intensity_integral` via the substitution `n0 → 1/n0, g → -g, min_pop → 1/min_pop`.
This works without the barrier, but NOT with it: flipping the sign of `g` swaps which
region is clamped vs unclamped (since the `max` branching depends on monotonicity),
while the actual clamped region is the same for both integrals. Formally,
`1/max(A,B) = min(1/A,1/B)`, not `max(1/A,1/B)`.

```cpp
auto Exp_pop_model::pop_integral(double a, double b) const -> double {
  CHECK_LE(a, b);
  if (min_pop_ == 0.0) {
    if (growth_rate_ == 0.0) {
      // no barrier, constant pop
      return (b - a) * pop_at_t0_;
    } else {
      // no barrier, exponential pop
      return pop_at_t0_ / growth_rate_
          * std::exp(growth_rate_ * (a - t0_))
          * std::expm1(growth_rate_ * (b - a));
    }
  } else if (growth_rate_ == 0.0) {
    // barrier active, constant pop
    return (b - a) * std::max(min_pop_, pop_at_t0_);
  } else if (growth_rate_ > 0.0) {
    if (b <= t_c_) {
      // g > 0, [a,b] entirely clamped
      return (b - a) * min_pop_;
    } else if (a >= t_c_) {
      // g > 0, [a,b] entirely unclamped
      return pop_at_t0_ / growth_rate_
          * std::exp(growth_rate_ * (a - t0_))
          * std::expm1(growth_rate_ * (b - a));
    } else {
      // g > 0, [a,b] straddles t_c: clamped in [a,t_c], unclamped in [t_c,b]
      return (t_c_ - a) * min_pop_
          + pop_at_t0_ / growth_rate_
          * std::exp(growth_rate_ * (t_c_ - t0_))
          * std::expm1(growth_rate_ * (b - t_c_));
    }
  } else {
    if (a >= t_c_) {
      // g < 0, [a,b] entirely clamped
      return (b - a) * min_pop_;
    } else if (b <= t_c_) {
      // g < 0, [a,b] entirely unclamped
      return pop_at_t0_ / growth_rate_
          * std::exp(growth_rate_ * (a - t0_))
          * std::expm1(growth_rate_ * (b - a));
    } else {
      // g < 0, [a,b] straddles t_c: unclamped in [a,t_c], clamped in [t_c,b]
      return pop_at_t0_ / growth_rate_
          * std::exp(growth_rate_ * (a - t0_))
          * std::expm1(growth_rate_ * (t_c_ - a))
          + (b - t_c_) * min_pop_;
    }
  }
}

auto Exp_pop_model::intensity_integral(double a, double b) const -> double {
  CHECK_LE(a, b);
  if (min_pop_ == 0.0) {
    if (growth_rate_ == 0.0) {
      // no barrier, constant pop
      return (b - a) / pop_at_t0_;
    } else {
      // no barrier, exponential pop
      return -1.0 / (pop_at_t0_ * growth_rate_)
          * std::exp(-growth_rate_ * (a - t0_))
          * std::expm1(-growth_rate_ * (b - a));
    }
  } else if (growth_rate_ == 0.0) {
    // barrier active, constant pop
    return (b - a) / std::max(min_pop_, pop_at_t0_);
  } else {
    auto inv_min_pop = 1.0 / min_pop_;

    if (growth_rate_ > 0.0) {
      if (b <= t_c_) {
        // g > 0, [a,b] entirely clamped
        return (b - a) * inv_min_pop;
      } else if (a >= t_c_) {
        // g > 0, [a,b] entirely unclamped
        return -1.0 / (pop_at_t0_ * growth_rate_)
            * std::exp(-growth_rate_ * (a - t0_))
            * std::expm1(-growth_rate_ * (b - a));
      } else {
        // g > 0, [a,b] straddles t_c: clamped in [a,t_c], unclamped in [t_c,b]
        return (t_c_ - a) * inv_min_pop
            - 1.0 / (pop_at_t0_ * growth_rate_)
            * std::exp(-growth_rate_ * (t_c_ - t0_))
            * std::expm1(-growth_rate_ * (b - t_c_));
      }
    } else {
      if (a >= t_c_) {
        // g < 0, [a,b] entirely clamped
        return (b - a) * inv_min_pop;
      } else if (b <= t_c_) {
        // g < 0, [a,b] entirely unclamped
        return -1.0 / (pop_at_t0_ * growth_rate_)
            * std::exp(-growth_rate_ * (a - t0_))
            * std::expm1(-growth_rate_ * (b - a));
      } else {
        // g < 0, [a,b] straddles t_c: unclamped in [a,t_c], clamped in [t_c,b]
        return -1.0 / (pop_at_t0_ * growth_rate_)
            * std::exp(-growth_rate_ * (a - t0_))
            * std::expm1(-growth_rate_ * (t_c_ - a))
            + (b - t_c_) * inv_min_pop;
      }
    }
  }
}
```

### 3. `core/run.h` - No New Configuration Members

The minimum population is stored entirely within the `Exp_pop_model` instance. No need to
duplicate it in the `Run` class. If access is needed, forward to `pop_model_`.

No changes to `run.h`.

### 4. `core/run.cpp` - Initialize and Use Minimum Population

#### Constructor Initialization (around line 20):

The default Exp_pop_model in the Run constructor currently uses the 3-arg constructor
(without min_pop). We need to update this to use the 4-arg constructor with the default
min_pop of 1.0. However, since we're NOT adding a member to Run, the constructor needs to
create the Exp_pop_model with min_pop=1.0 directly:

```cpp
Run::Run(ctpl::thread_pool& thread_pool, std::mt19937 bitgen, Phylo_tree tree)
    : thread_pool_{&thread_pool},
      bitgen_{bitgen},
      tree_{std::move(tree)},
      num_parts_{1},
      target_coal_prior_cells_{400},
      // Use 4-arg constructor with min_pop = 1.0
      pop_model_{std::make_shared<Exp_pop_model>(calc_max_tip_time(tree_), 1000.0, 0.0, 1.0)},
      // ... rest of initialization unchanged
{ ... }
```

#### MCMC moves

The `pop_size_move` and `pop_growth_rate_move` functions create new Exp_pop_model instances.
These need to be updated to pass the existing min_pop from the old model:

```cpp
// In pop_size_move():
pop_model_ = std::make_shared<Exp_pop_model>(
    old_pop_model->t0(), new_n0, old_pop_model->growth_rate(), old_pop_model->min_pop());

// In pop_growth_rate_move():
pop_model_ = std::make_shared<Exp_pop_model>(
    old_pop_model->t0(), old_pop_model->pop_at_t0(), new_g, old_pop_model->min_pop());
```

There are no getter/setter methods needed in Run since the minimum population is accessed
via the pop_model_ directly.

### 5. `core/api.fbs` - FlatBuffers Schema

Add new field to `ExpPopModel` table. Default to 0.0 so old files automatically disable the
feature:

Note: the existing `api.fbs` comment says n0 has "units of 1/day (it's obscure...)", but this
is wrong. `pop_at_time` returns `N(t) * rho` where rho is the generation time in days, so
the result (and n0) has units of **days**. Fix this comment on both the new `ExpPopModel`
table and the deprecated `pop_n0` field (around line 134). Also fix the same wrong comment
in `doc/dphy_file_format.md`, which duplicates the schema.

```fbs
table ExpPopModel {
  t0: float64      (id: 0); // Reference time t0 for population model [units: days since 2020-01-01]
  n0: float64      (id: 1); // Effective population size at t = t0 [units: days (= N_e * generation_time), > 0]
  g: float64       (id: 2); // Growth rate [units of 1 / day]
  min_pop: float64 (id: 3); // Minimum effective population size [units: days, >= 0, default 0.0 = disabled]
}
```

### 6. `core/api.cpp` - Serialization/Deserialization

#### Serialization (around line 230):

```cpp
pop_model_offset = api::CreateExpPopModel(
    fbb,
    exp_pop_model.t0(),
    exp_pop_model.pop_at_t0(),
    exp_pop_model.growth_rate(),
    exp_pop_model.min_pop()).o;
```

#### Deserialization (around line 340):

For deprecated pop_model (PopModel_NONE), we create an Exp_pop_model from the old fields.
We must pass a min_pop value. Since old files don't have this, we use 0.0 to disable the barrier
(reproducing old behavior):

```cpp
if (api_params->pop_model_type() == api::PopModel_NONE) {
  // Old format without explicit pop_model union
  run.set_pop_model(std::make_shared<Exp_pop_model>(
      api_params->pop_t0(),
      api_params->pop_n0(),
      api_params->pop_g(),
      0.0));  // OLD FILE: minimum population barrier disabled
} else {
  switch (api_params->pop_model_type()) {
    case api::PopModel_ExpPopModel: {
      auto api_exp_pop_model = api_params->pop_model_as_ExpPopModel();
      // For files created before this field was added, min_pop will be 0 (FlatBuffers default)
      // which correctly disables the barrier
      run.set_pop_model(std::make_shared<Exp_pop_model>(
          api_exp_pop_model->t0(),
          api_exp_pop_model->n0(),
          api_exp_pop_model->g(),
          api_exp_pop_model->min_pop()));
      break;
    }
    // ... Skygrid case unchanged
  }
}
```

No sentinel value needed since FlatBuffers will default to 0.0 for missing fields, which is exactly
what we want for old files (disabled barrier).

### 7. `core/cmdline.cpp` - Command Line Interface

Add new option (in the exponential model options group, around line 228-287):

```cpp
("v0-exp-pop-min-pop",
 "[pop-model == exponential] Minimum effective population size for exponential model, in years "
 "(0 = disabled, default = 1 day = 1/365 years).  Effective population sizes below this minimum "
 "are clamped to this value.  This prevents numerical issues in deep trees.",
 cxxopts::value<double>()->default_value(absl::StrFormat("%.5g", 1.0 / 365.0)))
```

In the exponential model configuration section, update the `Exp_pop_model` construction
(around line 824) to pass min_pop:

```cpp
auto min_pop = opts["v0-exp-pop-min-pop"].as<double>() * 365.0;  // convert to days
run->set_pop_model(std::make_unique<Exp_pop_model>(t0, init_final_pop_size, init_pop_growth_rate, min_pop));
```

### 8. `tools/delphy_wasm.cpp` - WASM Interface

Add a new getter and update the constructor to take a 4th parameter. The constructor
change breaks the WASM ABI, so delphy-web must be updated accordingly (see Section 9).

Add new getter function:

```cpp
EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_exp_pop_model_get_min_pop(Delphy_context& /*ctx*/, const Exp_pop_model& exp_pop_model) -> double {
  return exp_pop_model.min_pop();
}
```

Update constructor (no default - callers must be updated):

```cpp
EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_exp_pop_model_new(
    Delphy_context& /*ctx*/,
    double t0,
    double n0,
    double g,
    double min_pop)
    -> Pop_model* {
  return new Exp_pop_model{t0, n0, g, min_pop};
}
```

All existing callers of `delphy_exp_pop_model_new` need to be updated to pass the min_pop parameter.

### 9. `delphy-web` - TypeScript Interface (Separate Repository)

**Note:** `delphy-web` is a separate repository. The changes below describe companion
updates that must be made there after the Delphy core changes land and WASM is recompiled.

Three files need changes: `delphy_api.ts` (WASM bindings), `pythia.ts` (run parameter
config and model construction), and `runui.ts` (TODO for future UI work).

#### `delphy_api.ts` - Update ExpPopModel class:

```typescript
export class ExpPopModel extends PopModel {
  constructor(public t0: number, public n0: number, public g: number, public minPop: number) {
    super();
  }

  toPopModelPtr(ctx: DelphyContextPtr): PopModelPtr {
    return Delphy.delphyCoreRaw.exp_pop_model_new(
      ctx, this.t0, this.n0, this.g, this.minPop);
  }
}
```

#### `delphy_api.ts` - Update popModelPtrToPopModel:

```typescript
const t0 = Delphy.delphyCoreRaw.exp_pop_model_get_t0(ctx, rawPopModel);
const n0 = Delphy.delphyCoreRaw.exp_pop_model_get_n0(ctx, rawPopModel);
const g = Delphy.delphyCoreRaw.exp_pop_model_get_g(ctx, rawPopModel);
const minPop = Delphy.delphyCoreRaw.exp_pop_model_get_min_pop(ctx, rawPopModel);

return new ExpPopModel(t0, n0, g, minPop);
```

#### `delphy_api.ts` - Update delphyCoreRaw interface:

```typescript
exp_pop_model_new: (ctx: DelphyContextPtr, t0: number, n0: number, g: number, minPop: number) => PopModelPtr;
exp_pop_model_get_min_pop: (ctx: DelphyContextPtr, pop_model: ExpPopModelPtr) => number;
```

#### `pythia.ts` - Add minPop to RunParamConfig and wire it through:

Add `minPop: number` to the `RunParamConfig` type (in the exponential model params section).
Set its default to `1.0` (days) in `getEmptyRunParamConfig()`.

Update the `ExpPopModel` construction to pass `minPop`:

```typescript
run.setPopModel(
  new ExpPopModel(calcMaxDateOfTree(run.getTree()),
    runParams.finalPopSize,
    runParams.popGrowthRate,
    runParams.minPop));
```

Update `extractRunParamsFromRun` to read `minPop` back:

```typescript
result.popGrowthRate = popModel.g;
result.minPop = popModel.minPop;
```

#### `runui.ts` - Add TODO for future UI exposure:

In `getAdvancedFormValues()`, after the existing exponential model form fields (fixed pop
size and growth rate), add a TODO for future work:

```typescript
// TODO: expose minPop (minimum population barrier) in the Advanced Options dialog
```

### 10. `tools/beast_trees_to_dphy.cpp` - BEAST Tree Converter

Line 605 constructs an `Exp_pop_model` directly:

```cpp
fake_run.set_pop_model(std::make_shared<Exp_pop_model>(0.0, 42.0 / 365.0, 42.0 * 365.0));
```

Update to the 4-arg constructor with `min_pop = 0.0` (this is a dummy pop model for file
conversion, so the barrier is irrelevant):

```cpp
fake_run.set_pop_model(std::make_shared<Exp_pop_model>(0.0, 42.0 / 365.0, 42.0 * 365.0, 0.0));
```

### 11. `tools/delphy.cpp` - CLI Tool

`tools/delphy.cpp` references `Exp_pop_model` via a `typeid` check (line 77) and prints
`n0` and `g` for diagnostics (lines 79-80). Add `min_pop` to the diagnostic output:

```cpp
std::cerr << absl::StreamFormat("n0 = %.4g yr, ", exp_pop_model.pop_at_t0() / 365.0)
          << absl::StreamFormat("g = %.4f e-fold/yr, ", exp_pop_model.growth_rate() * 365.0)
          << absl::StreamFormat("min_pop = %.4g yr, ", exp_pop_model.min_pop() / 365.0);
```

### 12. `tools/delphy_ui.cpp` - Interactive UI Tool

`tools/delphy_ui.cpp` both reads and constructs `Exp_pop_model` instances.

Diagnostic output (lines 615-616) — add `min_pop`, same as `tools/delphy.cpp`:

```cpp
std::cerr << absl::StreamFormat("n0 = %.4g yr, ", exp_pop_model.pop_at_t0() / 365.0)
          << absl::StreamFormat("g = %.4g e-fold/yr, ", exp_pop_model.growth_rate() * 365.0)
          << absl::StreamFormat("min_pop = %.4g yr, ", exp_pop_model.min_pop() / 365.0);
```

Reset pop size (line 926) — pass existing `min_pop`:

```cpp
ui_run->set_pop_model(std::make_shared<Exp_pop_model>(
    exp_pop_model.t0(), 3.0 * 365.0, exp_pop_model.growth_rate(), exp_pop_model.min_pop()));
```

Reset growth rate (line 947) — pass existing `min_pop`:

```cpp
ui_run->set_pop_model(std::make_shared<Exp_pop_model>(
    exp_pop_model.t0(), exp_pop_model.pop_at_t0(), 0.0, exp_pop_model.min_pop()));
```

### 13. Files That Reference `Exp_pop_model` Without Constructing It (No Changes Needed)

The following files use `Exp_pop_model` only via `typeid` checks and `static_cast` for
reading accessors. They do not call the constructor and need no changes, but are listed
here for completeness:

- `core/beasty_output.cpp` — typeid checks + accessor reads for BEAST log output
- `core/beasty_input.cpp` — typeid checks + accessor reads for BEAST XML import

### 14. C++ Callers — Full Inventory

The 3-arg constructor is removed (no backwards compatibility). All callers that construct
`Exp_pop_model` instances:

- `core/run.cpp`: Run constructor (Section 4) and MCMC moves (Section 4)
- `core/api.cpp`: serialization/deserialization (Section 6)
- `core/cmdline.cpp`: CLI model creation (Section 7)
- `tools/delphy_wasm.cpp`: WASM constructor (Section 8)
- `tools/delphy_ui.cpp`: diagnostic output and reset commands (Section 12)
- `tools/beast_trees_to_dphy.cpp`: dummy model for file conversion (Section 10)
- `tests/pop_model_tests.cpp`: unit tests (lines 41, 42, 49, 75, 694)
- `tests/coalescent_tests.cpp`: coalescent test (line 83)
- `tests/scalable_coalescent_tests.cpp`: scalable coalescent test (line 64)
- `tests/very_scalable_coalescent_tests.cpp`: very scalable coalescent test (line 103)
- `tests/ancestral_tree_prober_tests.cpp`: prober test (line 43)
- `tests/site_states_tree_prober_tests.cpp`: site states test (line 151)

Use the 4-arg constructor with min_pop = 1.0 for new default scenarios, or 0.0 to
disable (e.g., in tests that need to reproduce old behavior exactly).

## Backwards Compatibility

### C++ Library

- Old code calling `Exp_pop_model(t0, n0, g)` will not compile - this is intentional
  to force finding and updating all callers
- New files: min_pop explicitly set to 1.0 in default scenarios

### FlatBuffers Files

- Old files (without min_pop field): FlatBuffers defaults to 0.0, which disables the barrier
  (reproduces old behavior)
- New files: min_pop explicitly written

### CLI

- Old runs (no `--v0-exp-pop-min-pop` option): default to 1.0 day
- Old command-line tools reading old files: min_pop = 0.0 (disabled)

### WASM/delphy-web

- Old WASM modules: won't have `exp_pop_model_get_min_pop` or the 4-arg constructor
- delphy-web must be recompiled against new Delphy core
- All callers in delphy-web that construct ExpPopModel must be updated
- For old `.dphy` files loaded by new delphy-web: min_pop will be read as 0.0 (disabled)
- For new runs in delphy-web: default min_pop = 1.0

## Testing Strategy

1. **Unit tests** for `Exp_pop_model` in a new test file or expanded existing tests:
   - `pop_at_time()` never returns < min_pop
   - `pop_at_time()` returns min_pop when unclamped value < min_pop
   - `pop_integral()` with:
     - g > 0, interval entirely below crossing point
     - g > 0, interval entirely above crossing point
     - g > 0, interval straddling crossing point
     - g < 0, interval entirely below crossing point
     - g < 0, interval entirely above crossing point
     - g < 0, interval straddling crossing point
     - g = 0, constant population above min_pop
     - g = 0, constant population below min_pop
   - `intensity_integral()` with the same scenarios
   - Edge cases: min_pop = 0.0, min_pop = n_0, very small |g|, min_pop >> n_0

2. **Integration test**: MCMC run with exponential model and min_pop, verify no crashes

3. **FlatBuffers round-trip**: write params with min_pop > 0, read back, verify min_pop preserved

4. **Old file compatibility**: load old `.dphy` file, verify min_pop = 0.0

5. **CLI test**: verify `--v0-exp-pop-min-pop` option works correctly

## Migration Path for delphy-web

1. Update Delphy core with this change
2. Recompile to WASM
3. In delphy-web:
   - Update TypeScript interface as described above
   - Find and update ALL code that creates `ExpPopModel` to pass minPop parameter
   - The Run parameter parsing code that reads from saved files will automatically get
     minPop = 0.0 for old files
   - Add TODO comment to expose minPop in Advanced Options UI
4. Test:
   - New runs work with minPop
   - Old `.dphy` files load (minPop = 0.0)
   - New `.dphy` files save with minPop

## File Changes Summary

| File | Changes |
|------|---------|
| `core/pop_model.h` | Add `min_pop_` and `t_c_` members, new 4-arg constructor, `min_pop()` accessor, update inline `print_to` |
| `core/pop_model.cpp` | Update all methods with min_pop handling, precompute t_c_ in constructor |
| `core/run.cpp` | Update Run constructor and MCMC moves to pass min_pop |
| `core/api.fbs` | Add `min_pop` field to `ExpPopModel` table with default=0.0, fix n0 units comments |
| `doc/dphy_file_format.md` | Add `min_pop` field to `ExpPopModel` table, fix n0 units comments (mirrors `api.fbs`) |
| `core/api.cpp` | Update serialization/deserialization with min_pop handling |
| `core/cmdline.cpp` | Add `--v0-exp-pop-min-pop` CLI option, update exp model creation |
| `tools/delphy_wasm.cpp` | Add `delphy_exp_pop_model_get_min_pop`, update constructor to 4-arg |
| `tools/delphy.cpp` | Add `min_pop` to diagnostic output (lines 79-80) |
| `tools/delphy_ui.cpp` | Add `min_pop` to diagnostic output (lines 615-616), pass `min_pop` in reset commands (lines 926, 947) |
| `tools/beast_trees_to_dphy.cpp` | Update 3-arg constructor call (line 605) to 4-arg with min_pop=0.0 |
| `tests/pop_model_tests.cpp` | Update existing 3-arg constructor calls, add new min_pop unit tests |
| `tests/coalescent_tests.cpp` | Update 3-arg constructor call (line 83) to 4-arg |
| `tests/scalable_coalescent_tests.cpp` | Update 3-arg constructor call (line 64) to 4-arg |
| `tests/very_scalable_coalescent_tests.cpp` | Update 3-arg constructor call (line 103) to 4-arg |
| `tests/ancestral_tree_prober_tests.cpp` | Update 3-arg constructor call (line 43) to 4-arg |
| `tests/site_states_tree_prober_tests.cpp` | Update 3-arg constructor call (line 151) to 4-arg |
| `delphy-web/.../delphy_api.ts` | *(Separate repo)* Update TypeScript ExpPopModel class, add minPop, update WASM bindings |
| `delphy-web/.../pythia.ts` | *(Separate repo)* Add `minPop` to `RunParamConfig`, wire through construction and extraction |
| `delphy-web/.../runui.ts` | *(Separate repo)* Add TODO for exposing minPop in Advanced Options dialog |

The following files reference `Exp_pop_model` only via `typeid` checks and accessor reads
(no constructor calls) and need **no changes**:

| File | Usage |
|------|-------|
| `core/beasty_output.cpp` | typeid checks + accessor reads for BEAST log output |
| `core/beasty_input.cpp` | typeid checks + accessor reads for BEAST XML import |

