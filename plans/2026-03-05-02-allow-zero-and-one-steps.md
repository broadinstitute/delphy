# Allow `--v0-steps 0` and `--v0-steps 1`

## Problem

Running Delphy with `--v0-steps 0` (to inspect initial state without any MCMC) or
`--v0-steps 1` (to run a single MCMC step) both fail:

- `--v0-steps 0`: rejected by validation at cmdline.cpp:384 (`steps <= 0`).
- `--v0-steps 1`: the default `log_every = steps / 10000 = 0` (integer division),
  which is then rejected by validation at cmdline.cpp:398 (`log_every <= 0`).
  Same for `tree_every = steps / 1000` and `delphy_snapshot_every = steps / 1000`.

## Desired behavior

- `--v0-steps 0`: parse inputs, build initial tree, print one status line with
  initial parameter values, then exit.  No MCMC moves performed.
- `--v0-steps 1`: perform exactly one MCMC step, printing a status line before and
  after (if `log_every` permits).

## Analysis

The driver loop in `tools/delphy.cpp:165-202` already handles these cases correctly
once it gets past validation:
- For `steps=0`: `step() == 0 >= 0` breaks immediately after the first status line.
- For `steps=1` with `step_granularity=1`: does one step, then breaks.

The `%` operations for `log_now`, `tree_now`, `delphy_output_now` require their
divisors to be nonzero, but this is guaranteed as long as the `<= 0` checks pass.

So the only changes needed are in `core/cmdline.cpp`:

## Fix

### cmdline.cpp: Allow `steps >= 0` and clamp `*_every` defaults to at least 1

1. **Line 384**: Change `steps <= 0` to `steps < 0`.

2. **Lines 394, 408, 437**: Clamp the default values of `log_every`, `tree_every`,
   and `delphy_snapshot_every` to at least 1:
   ```cpp
   auto log_every = std::max(1L, steps / 10000);
   auto tree_every = std::max(1L, steps / 1000);
   auto delphy_snapshot_every = std::max(1L, steps / 1000);
   ```

No changes needed to `tools/delphy.cpp` or `core/run.cpp`.
