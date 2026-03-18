# Plan: Use `to_linear_year` consistently in BEAST log output

## Background

Delphy measures time internally in days since 2020-01-01.  BEAST uses
"years" as its time unit.  Converting between the two requires care
because of leap years: `to_linear_year(t)` in `core/dates.cpp` does a
proper calendar-aware conversion (`year + day_of_year / days_in_year`),
while the simpler `t / 365.0` or `2020.0 + t / 365.0` does not.

The difference is tiny — about 1-2 days over the 2020-2025 range —
but it showed up as a statistically significant bias in the WCSS
tip-date uncertainty study (study 11), which has high statistical
power due to thousands of tip-date samples.  The root cause: the
BEAST X log output in `beasty_output.cpp` mixed the two conventions,
using `to_linear_year` for `age(root)` but simple `/ 365` division
for `rootHeight` and `age(TIP)`.  This made the downstream formula
`rootHeight + age(root) - age(TIP)` produce posterior tip calendar
years that were systematically ~2 days too low.

## Scope

This fix applies only to BEAST **log file** output (`beasty_output.cpp`).

Newick **tree output** (`output_newick_tree`) and tree **input**
(`newick_to_phylo_tree` in `beasty_input.cpp`) are left unchanged.
Newick trees encode only relative time differences (branch lengths),
not absolute node times, so there is no way to know where leap years
fall when reading a tree back in.  The simple `/ 365` and `* 365`
conversions must stay consistent with each other to ensure lossless
roundtripping (e.g., when `delphy_mcc` reads back trees that Delphy
wrote).  The tiny leap-year error in branch lengths is
inconsequential.

Similarly, Newick output in `api.cpp` and BEAST XML generation in
`beasty_input.cpp` are left unchanged — they are either tree I/O
(must roundtrip) or rate/magnitude unit conversions where `* 365` /
`/ 365` is correct.

## Principle

**Every conversion of an absolute time or time delta to "years" in
BEAST log output should go through `to_linear_year` or
`delta_linear_years`.**  No bare `/ 365` or `/ 365.0` divisions.

For *rates* (mutations/site/year, growth rate) and *magnitudes*
(population size), the `* 365` / `/ 365` conversion remains correct
and is left unchanged.

## Changes to `core/beasty_output.cpp`

### Static helper

Add a static helper at the top of `beasty_output.cpp`:

```cpp
// Convert a time delta from days-since-epoch to linear years.
//
// Bare `/ 365` division ignores leap years, introducing a systematic
// error of ~1-2 days over multi-year spans.  This was detected by the
// WCSS tip-date uncertainty study (study 11), where mixing `/ 365`
// and `to_linear_year` in BEAST log output caused posterior tip
// calendar years to be systematically ~2 days too low.
//
// Use this for all time deltas in BEAST log output (rootHeight, etc.)
// instead of `(t_to - t_from) / 365.0`.
static auto delta_linear_years(double t_from, double t_to) -> double {
  return to_linear_year(t_to) - to_linear_year(t_from);
}
```

### `Beasty_log_output_2_x_x::output_log`

The non-mpox path:

1. **Line 163 — `TreeHeight`**: `(beast_t0 - tree.at_root().t) / 365`
   → `delta_linear_years(tree.at_root().t, beast_t0)`

2. **Line 184 — `ePopSize`**: `pop_at_time(beast_t0)/365` — magnitude
   conversion.  **No change.**

3. **Line 189 — `growthRate`**: `pop_growth_rate*365` — rate
   conversion.  **No change.**

4. **Line 207 — tip logP**: `to_linear_year(t_max) - to_linear_year(t_min)`
   — already correct.  **No change.**

5. **Line 210 — tip height**: `to_linear_year(tree.at(node).t)` —
   already correct.  **No change.**

The mpox path (`mpox_hack_enabled == true`):

6. **Line 218 — `to_beast_date` lambda**: defined as
   `2020.0 + delphy_t / 365.0`, used on line 226 for `age(root)`.
   **Remove the lambda entirely** and replace its usage on line 226
   with `to_linear_year(tree.at_root().t)`.

7. **Line 225 — mpox `treeModel.rootHeight`**:
   `(beast_t0 - tree.at_root().t) / 365`
   → `delta_linear_years(tree.at_root().t, beast_t0)`

8. **Line 227 — mpox `treeLength`**: `calc_T(tree) / 365.0` — sum of
   branch lengths; would require a new function to fix properly.
   **Leave as-is.**

### `Beasty_log_output_X_10_5_0::output_log`

1. **Line 362 — `rootHeight`**: `(beast_t0 - tree.at_root().t) / 365`
   → `delta_linear_years(tree.at_root().t, beast_t0)`

2. **Line 363 — `age(root)`**: `to_linear_year(tree.at_root().t)` —
   already correct.  **No change.**

3. **Line 364 — `treeLength`**: `calc_T(tree) / 365.0` — same as
   above, **leave as-is.**

4. **Line 372 — `exponential.popSize`**: magnitude conversion.
   **No change.**

5. **Line 377 — `exponential.growthRate`**: rate conversion.
   **No change.**

6. **Line 409 — `age(TIP)`**: `(beast_t0 - tree.at(node).t) / 365.0`
   → `delta_linear_years(tree.at(node).t, beast_t0)`.  **This is the
   key fix** — the one that caused the WCSS tip-date coverage problem.

### Not changed

- **`output_skygrid_results` line 43** — `skygrid.cutOff`: constant
  parameter, not tied to a calendar position.  **Leave as-is.**

- **`output_newick_tree` lines 596, 603** — Newick branch lengths and
  mutation times.  Must stay as `/ 365.0` to roundtrip with
  `newick_to_phylo_tree` (`beasty_input.cpp:51`).  **Leave as-is.**

- **All rate/magnitude conversions** (`mu * 365`, `pop_size / 365`,
  `growth_rate * 365`, `log(365.0)` for coalescent prior units).
  **Leave as-is.**

## Summary of actual code changes

| Location | Old | New |
|----------|-----|-----|
| `beasty_output.cpp` | — | Add static `delta_linear_years` helper |
| line 163 (BEAST 2) | `(beast_t0 - root.t) / 365` | `delta_linear_years(root.t, beast_t0)` |
| line 218 (BEAST 2 mpox) | `to_beast_date` lambda | Remove entirely |
| line 225 (BEAST 2 mpox) | `(beast_t0 - root.t) / 365` | `delta_linear_years(root.t, beast_t0)` |
| line 226 (BEAST 2 mpox) | `to_beast_date(root.t)` | `to_linear_year(root.t)` |
| line 362 (BEAST X) | `(beast_t0 - root.t) / 365` | `delta_linear_years(root.t, beast_t0)` |
| line 409 (BEAST X) | `(beast_t0 - tip.t) / 365.0` | `delta_linear_years(tip.t, beast_t0)` |
