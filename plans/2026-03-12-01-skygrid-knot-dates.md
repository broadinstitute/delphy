# Explicit Skygrid Knot Date Options

## Motivation

Currently, the Skygrid knot endpoints are specified indirectly:

- `x_M` (last/most-recent knot) is implicitly `t0 = calc_max_tip_time(tree)`.
- `x_0` (first/oldest knot) is `t0 - cutoff`, where cutoff comes from
  `--v0-skygrid-cutoff` (in years).

This is awkward: the user must mentally compute the desired time range as an offset
from the last tip time.  It would be much more natural to specify the endpoints
directly as calendar dates.

This plan adds two new CLI options:

- `--v0-skygrid-first-knot-date YYYY-MM-DD` — sets `x_0` directly.
- `--v0-skygrid-last-knot-date YYYY-MM-DD` — sets `x_M` directly.


## Current Implementation

### Knot computation (cmdline.cpp:836-848)

```cpp
auto skygrid_cutoff = opts["v0-skygrid-cutoff"].as<double>();  // years
// ... validation ...
skygrid_cutoff *= 365.0;  // convert to days

auto M = num_parameters - 1;
auto x_k = std::vector<double>(M+1, 0.0);
auto dt = (skygrid_cutoff / M);
for (auto k = 0; k <= M; ++k) {
  x_k[k] = (t0 - skygrid_cutoff) + k * dt;
}
```

Where `t0 = calc_max_tip_time(tree)` (cmdline.cpp:608).

### Downstream uses of `skygrid_cutoff`

After computing `x_k`, `skygrid_cutoff` is never used again directly.  It appears
only through `dt = skygrid_cutoff / M`, which is used in the tau computation from
`--v0-skygrid-prior-double-half-time` (cmdline.cpp:898-899):

```cpp
const auto D = std::pow(std::log(2.0), 2) / (2 * double_half_time);
skygrid_tau = 1.0 / (2 * D * dt);
```

Verified: `skygrid_cutoff` appears only on lines 836-847 (definition, validation,
unit conversion, and the `x_k` loop).  After that block, only `dt` is referenced
(lines 889, 899).  Since `dt = skygrid_cutoff / M = (x_M - x_0) / M`, we can
compute `dt` directly from `x_0` and `x_M` and eliminate the `skygrid_cutoff`
variable entirely.

### Date parsing

`core/dates.h` provides `parse_iso_date(std::string_view)` which converts a
`YYYY-MM-DD` string to days since 2020-01-01 (Delphy's internal time unit).  This is
already used throughout the codebase.


## Interaction with `--v0-skygrid-cutoff`

The two ways of specifying knot endpoints are **mutually exclusive**:

1. **Old way:** `--v0-skygrid-cutoff` — sets `x_0 = t0 - cutoff`, `x_M = t0`.
2. **New way:** `--v0-skygrid-first-knot-date` + `--v0-skygrid-last-knot-date` —
   sets both endpoints explicitly.

If either date option is specified, cutoff must not be, and vice versa.  The new date
options must always be specified together (both or neither).

### Summary of valid combinations

| first-knot-date | last-knot-date | cutoff      | x_0                     | x_M                    |
|-----------------|----------------|-------------|-------------------------|------------------------|
| not set         | not set        | set or dflt | `t0 - cutoff`           | `t0`                   |
| set             | set            | must not set| `parse_iso_date(first)` | `parse_iso_date(last)` |

Any other combination is an error.


## Implementation Plan

All changes are in `core/cmdline.cpp` unless otherwise noted.

### Step 1: Add CLI option definitions (~after line 297)

Add near the existing `--v0-skygrid-cutoff` definition:

```cpp
("v0-skygrid-first-knot-date",
 "[pop-model == skygrid] Date (YYYY-MM-DD) of the first (oldest) Skygrid knot x_0.  "
 "Example: --v0-skygrid-first-knot-date 2013-12-01.  "
 "Must be specified together with --v0-skygrid-last-knot-date.  "
 "Mutually exclusive with --v0-skygrid-cutoff.",
 cxxopts::value<std::string>())
("v0-skygrid-last-knot-date",
 "[pop-model == skygrid] Date (YYYY-MM-DD) of the last (most recent) Skygrid knot x_M.  "
 "Example: --v0-skygrid-last-knot-date 2014-10-01.  "
 "Must be specified together with --v0-skygrid-first-knot-date.  "
 "Mutually exclusive with --v0-skygrid-cutoff.",
 cxxopts::value<std::string>())
```

### Step 2: Add to `has_skygrid_pop_model_parameters` detection (~line 633)

Add two new lines to the existing detection block:

```cpp
(opts.count("v0-skygrid-first-knot-date") > 0) ||
(opts.count("v0-skygrid-last-knot-date") > 0) ||
```

### Step 3: Replace knot computation block (~lines 836-848)

Replace the existing block with:

```cpp
auto has_first_knot_date = opts.count("v0-skygrid-first-knot-date") > 0;
auto has_last_knot_date = opts.count("v0-skygrid-last-knot-date") > 0;
auto has_cutoff = opts.count("v0-skygrid-cutoff") > 0;

// Date options and cutoff are mutually exclusive
if ((has_first_knot_date || has_last_knot_date) && has_cutoff) {
  std::cerr << "ERROR: --v0-skygrid-first-knot-date / --v0-skygrid-last-knot-date "
            << "and --v0-skygrid-cutoff are mutually exclusive\n";
  std::exit(EXIT_FAILURE);
}

// Date options must be specified together
if (has_first_knot_date != has_last_knot_date) {
  std::cerr << "ERROR: --v0-skygrid-first-knot-date and --v0-skygrid-last-knot-date "
            << "must be specified together\n";
  std::exit(EXIT_FAILURE);
}

// Determine x_0 and x_M
auto x_0 = 0.0;
auto x_M = 0.0;
if (has_first_knot_date) {
  x_0 = parse_iso_date(opts["v0-skygrid-first-knot-date"].as<std::string>());
  x_M = parse_iso_date(opts["v0-skygrid-last-knot-date"].as<std::string>());
} else {
  auto skygrid_cutoff = opts["v0-skygrid-cutoff"].as<double>();  // years
  if (skygrid_cutoff < 0.0) {
    std::cerr << "ERROR: Skygrid cutoff must be positive\n";
    std::exit(EXIT_FAILURE);
  }
  x_M = t0;
  x_0 = t0 - skygrid_cutoff * 365.0;
}

if (x_0 >= x_M) {
  std::cerr << "ERROR: Skygrid first knot (" << to_iso_date(x_0) << ") must be before "
            << "last knot (" << to_iso_date(x_M) << ")\n";
  std::exit(EXIT_FAILURE);
}

auto M = num_parameters - 1;
auto x_k = std::vector<double>(M+1, 0.0);
auto dt = (x_M - x_0) / M;
for (auto k = 0; k <= M; ++k) {
  x_k[k] = x_0 + k * dt;
}
```

The old `skygrid_cutoff` variable is eliminated.  `dt` is computed directly as
`(x_M - x_0) / M`, which works regardless of how `x_0` and `x_M` were determined.

### Step 4: Add `#include "dates.h"` if not already present

`cmdline.cpp` likely already includes `dates.h` (it calls `to_iso_date` elsewhere),
but verify and add if missing.

### Step 5: Build and test

**5a.** Build debug, run unit tests.

**5b.** Verify existing behavior is preserved (no date options, cutoff only):

```bash
./build/debug/delphy \
  --v0-in-fasta ~/now/delphy-2026-paper-data/ebola-gire-2014/delphy_inputs/ebola.fasta \
  --v0-pop-model skygrid --v0-skygrid-num-parameters 10 --v0-skygrid-cutoff 1.5 \
  --v0-steps 100
```

**5c.** Test new date options (both dates):

```bash
./build/debug/delphy \
  --v0-in-fasta ~/now/delphy-2026-paper-data/ebola-gire-2014/delphy_inputs/ebola.fasta \
  --v0-pop-model skygrid --v0-skygrid-num-parameters 10 \
  --v0-skygrid-first-knot-date 2013-12-01 --v0-skygrid-last-knot-date 2014-10-01 \
  --v0-steps 100
```

**5d.** Test error cases:

- Both dates + cutoff: should error.
- Only one date option: should error.
- first-knot-date >= last-knot-date: should error.


## What Does NOT Change

- `Skygrid_pop_model` class (`core/pop_model.h`): already accepts arbitrary `x_k`.
- `Run` class (`core/run.h`, `core/run.cpp`): no changes.
- Tau computation: uses `dt = (x_M - x_0) / M`, which is available regardless of how
  `x_0` and `x_M` were determined.
- N_bar priors, low-pop barrier: unchanged.
- FlatBuffers schema (`core/api.fbs`): the knot positions `x_k` are already serialized
  as part of the `Skygrid_pop_model`; no new fields needed.
- BEAST XML export: knot positions flow through unchanged.
