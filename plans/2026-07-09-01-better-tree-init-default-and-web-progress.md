# Better Tree Init: Make It the Default and Expose Staged Progress in the Web UI

This plan makes the "Better tree init" pipeline (Rounds 1-6, ending in
[Round 6](2026-06-11-01-better-tree-init-round6-unified-spr-refinement.md)) the default
initialization method (both in the web UI and on the CLI), and expands the web loading
progress UI to reflect the four internal build stages of `mp-plus-timing` (guide tree,
refined tree, SPR refinement, and rooting + timing).

## Background

The new `mp-plus-timing` initialization pipeline (`build_initial_phylo_tree`,
`core/utree.cpp:1857`) builds an initial `Phylo_tree` in several internal stages:

1. **Guide tree** (`build_guide_tree`): add tips one at a time in input order. Progress is
   `(tips_so_far, total_tips)`.
2. **Refined guide trees** (`build_refined_tree`, up to `k_max_refinement_rounds = 5`
   rounds, stopping early when a round no longer reduces the total delta count): each round
   re-inserts all tips in a better order. Progress is `(round, tips_so_far, total_tips)`.
3. **SPR refinement** (`spr_refine`, up to `30*N` attempts): detach and regraft random
   tips/subtrees. Progress is `(attempt, max_attempts, cur_deltas)`, where `cur_deltas` is
   the current total number of mutations on the tree.
4. **Rooting + conversion** (`ols_regression_root_utree` + `utree_to_phylo_tree`): rooting
   does three passes over the tree (Pass 1 bottom-up Euler tour, Pass 2 top-down Euler tour,
   Pass 3 per-edge root evaluation). For large trees this is slow enough to warrant a
   progress bar. To avoid baking the current three-pass structure into the API, progress is
   reported as `(substage_id, substage, num_substages, nodes_in_substage, total_in_substage)`:
   - `substage_id`: a **stable integer identifying which pass** (see the `Rooting_substage`
     enum below). The UI maps this id to human-readable text on its side, so no UI strings
     ever cross the core/JS boundary. Unknown ids degrade gracefully to a generic label.
   - `substage`: the 1-based ordinal of the current pass (1..`num_substages`), used only to
     compute the global fraction.
   - `num_substages`: total number of passes (3 for `ols_regression_root_utree`).
   - `nodes_in_substage` / `total_in_substage`: progress within the current pass.

   For `ols_regression_root_utree`, `num_substages = 3` and each pass reports
   `(id, pass, 3, nodes_visited_this_pass, num_nodes)`. A future rooting/timing algorithm
   with a different set of passes just reports different `substage_id`s and `num_substages`;
   the UI computes a global fraction as `((substage - 1) + nodes/total) / num_substages` and
   needs no change beyond (optionally) adding labels for any new `substage_id`s. Rooting also
   reports a one-shot `Rooting_info` via the existing `rooting_hook` on completion.

   The stable substage ids (defined in `core/utree.h` next to `Rooting_info`):

   ```cpp
   enum class Rooting_substage : int {
     bottom_up_timing = 1,       // Pass 1: bottom-up (post-order) subtree-stats DFS
     top_down_timing = 2,        // Pass 2: top-down (pre-order) subtree-stats DFS
     root_candidate_eval = 3,    // Pass 3: per-edge root-position R^2 evaluation
   };
   ```

   UI labels (owned by the web repo): `bottom_up_timing` -> "bottom-up timing",
   `top_down_timing` -> "top-down timing", `root_candidate_eval` -> "root candidate
   evaluation" (they render as "Rooting and timing: bottom-up timing", etc.).

The core entry point already exposes hooks for build stages 1-3 (`core/cmdline.h:29`,
`core/utree.h:304`):

```cpp
auto build_rough_initial_tree_from_maple(
    Maple_file&& in_maple,
    Init_method init_method,
    absl::BitGenRef bitgen,
    const std::function<void(int,int)>&           progress_hook,              // guide tree
    const std::function<void(int,int,int)>&       refined_tree_progress_hook, // refined rounds
    const std::function<void(int,int,int)>&       spr_refine_progress_hook,   // SPR refine
    const std::function<void(const Rooting_info&)>& rooting_hook)             // rooting (completion only)
    -> Phylo_tree;
```

The build-stage-1-3 hooks already exist, but two core-side changes are still needed, both in
**Part A** below: (1) flipping the CLI default init method, and (2) adding the rooting
progress hook -- `ols_regression_root_utree` currently takes no progress hook, and neither
`build_initial_phylo_tree` nor `build_rough_initial_tree_from_maple` has a rooting-*progress*
parameter (they only have the one-shot completion `rooting_hook`). The remaining work is in
the WASM API (`tools/delphy_wasm.cpp`, **Part B**) and the web UI
(`~/github/fathominfo/delphy-web`, **Part C**).

### Current state of the WASM loading path

The two WASM entry points `delphy_parse_fasta_into_initial_tree_async` and
`delphy_parse_maple_into_initial_tree_async` (`tools/delphy_wasm.cpp:170`, `:257`):

- Hardcode `auto init_method = Init_method::old_usher_like;` (`:218`, `:292`).
- Wire only a single build-progress hook (`initial_build_progress_hook_id`) to the core's
  `progress_hook`, leaving the other three hooks defaulted.

The coarse `stage_progress_hook` reports which top-level stage the loader is in: for FASTA,
`{1: read, 2: analyze, 3: build}`; for MAPLE, `{1: read, 2: build}`. This plan expands the
single "build" stage into four coarse stages (guide tree / refined tree / SPR refinement /
rooting and timing); see "Coarse stage expansion" under Part B for the new numbering and the
consequences.

### Current state of the web loading path

- `src/ts/pythia/delphy_api.ts`: `parseFastaIntoInitialTreeAsync` / `parseMapleIntoInitialTreeAsync`
  register JS hook functions, pass their IDs to the raw WASM calls. Raw signatures are also
  declared here (`:486`, `:496`) and wired via `wrapRawAsyncApi` (`:818`).
- `src/ts/pythia/pythia.ts`: `initRunFromFasta` / `initRunFromMaple` (`:227`, `:250`) thread
  the callbacks through.
- `src/ts/ui/uploadui.ts`: defines the actual callbacks that drive a single progress bar
  (`#uploader--progress`) via `showProgress(label, total, soFar)` (`:49`). Currently:
  `parseProgressCallback` (read), `analysisProgressCallback` (analyze),
  `initTreeProgressCallback` (build). Call sites at `:290`, `:295`, `:509`, `:523`, `:543`.

## Goal

1. Make `mp-plus-timing` the default init method used by the WASM loading entry points,
   while keeping it a one-line change to fall back to `old_usher_like`.
2. Flip the **CLI** default (`--v0-init`) to `mp-plus-timing` too, with matching option-doc
   updates.
3. Add a **rooting progress hook** to the core rooter and thread it through to the UI.
4. Replace the single `initial_build_progress_hook_id` in both WASM entry points with four
   hooks: guide-tree, refined-tree, SPR-refine, and rooting progress.
5. Thread those hooks through `delphy_api.ts`, `pythia.ts`, and `uploadui.ts`, render
   distinct, informative progress-bar labels for each build stage, and expand the coarse
   stage indicator to match.

Non-goals: the web UI cannot select the init algorithm (so no fallback hooks are needed for
`old_usher_like`).

## Part A: Core changes (delphy `core/`)

### A1. Flip the CLI default to `mp-plus-timing`

In `process_args` (`core/cmdline.cpp:434`), the default is currently `old_usher_like`:

```cpp
auto init_method = Init_method::old_usher_like;   // -> Init_method::mp_plus_timing
```

Two coupled changes are required so the deprecated flags keep working:

- The deprecated `--v0-init-heuristic` flag means "old-usher-like". Today it relies on the
  default already being `old_usher_like` and sets nothing. After the flip it must set the
  method explicitly. Add a branch alongside the existing `--v0-init-random` handling
  (`core/cmdline.cpp:448`):
  ```cpp
  } else if (opts.count("v0-init-heuristic") > 0) {
    init_method = Init_method::old_usher_like;
  }
  ```
- Update the `--v0-init` option help string (`core/cmdline.cpp:178`) to end with
  `(default: mp-plus-timing)` instead of `(default: old-usher-like)`.

(Leave the `--v0-init-heuristic` / `--v0-init-random` "[deprecated ...]" help text as-is;
their meaning is unchanged.)

### A2. Rooting progress hook: `Rooting_substage` enum

Add to `core/utree.h`, next to `Rooting_info`:

```cpp
// Stable identifiers for the passes of a rooting/timing algorithm.  The UI maps these to
// human-readable labels on its side, so no UI strings cross the core/JS boundary.  Add new
// ids here (never renumber existing ones) if a future rooter has different passes.
enum class Rooting_substage : int {
  bottom_up_timing = 1,       // Pass 1: bottom-up (post-order) subtree-stats DFS
  top_down_timing = 2,        // Pass 2: top-down (pre-order) subtree-stats DFS
  root_candidate_eval = 3,    // Pass 3: per-edge root-position R^2 evaluation
};
```

### A3. Rooting progress hook: `ols_regression_root_utree`

Add a progress-hook parameter (defaulted, so existing callers/tests are unaffected):

```cpp
auto ols_regression_root_utree(
    Utree& tree, const std::vector<Tip_desc>& tip_descs, absl::BitGenRef bitgen,
    const std::function<void(int,int,int,int,int)>& progress_hook = [](int,int,int,int,int){})
    -> Rooting_info;
```

The five ints are `(substage_id, substage, num_substages, nodes_in_substage,
total_in_substage)` as defined in the Background. Inside the implementation
(`core/utree.cpp:1255`):

- Compute `num_nodes = tree.num_tips + tree.num_inner_nodes_so_far` and
  `num_edges = num_nodes - 1`. Set `num_substages = 3`.
- Throttle: emit at most ~100 updates per pass. `auto emit_every = std::max(1, num_nodes/100);`
  Emit once at the start of each pass (`nodes = 0`, so the UI can set the label immediately),
  every `emit_every` items, and once at pass end (`nodes = total`).
- **Pass 1** (`core/utree.cpp:1319`, the `leaving`-direction branch): count processed nodes;
  emit `(bottom_up_timing, 1, 3, n_done, num_nodes)`.
- **Pass 2** (`:1344`, the `entering`-direction branch): emit
  `(top_down_timing, 2, 3, n_done, num_nodes)`.
- **Pass 3** (`:1371`, the per-arc-pair loop, counting only non-free pairs): emit
  `(root_candidate_eval, 3, 3, n_done, num_edges)`.

The `midpoint_root_utree` fallback (taken for `N <= 2`, `var_t <= 0`, or no viable
candidate) is for tiny/degenerate trees and stays unhooked; the UI's progress bar simply
completes when the next coarse stage fires. Only `ols_regression_root_utree` gets the hook,
since that is what the pipeline uses; the same pattern extends to `gls_regression_root_utree`
if it ever becomes the default.

### A4. Thread the hook through the builders

Both `build_initial_phylo_tree` (declared in `core/utree.h`, defined in `core/utree.cpp`)
and `build_rough_initial_tree_from_maple` (declared in `core/cmdline.h`, defined in
`core/cmdline.cpp`) gain a `rooting_progress_hook` parameter of type
`std::function<void(int,int,int,int,int)>`, defaulted to a no-op, placed **immediately
before** the existing `rooting_hook` parameter. Placing it before `rooting_hook` lets a
caller pass the four progress hooks positionally and leave the completion `rooting_hook`
defaulted (used by the WASM entry points in Part B).

- In `build_initial_phylo_tree`, forward the new parameter into `ols_regression_root_utree`
  (replacing the current no-hook call).
- In `build_rough_initial_tree_from_maple`, forward the new parameter to
  `build_initial_phylo_tree` (only the `mp_plus_timing` branch uses it; the `random` and
  `old_usher_like` branches ignore it).

### A5. Existing CLI caller and an incidental warning fix

Threading the new parameter through touches two more spots in `core/cmdline.cpp`:

- **CLI caller.** `process_args` already calls `build_rough_initial_tree_from_maple` with
  positional progress-hook lambdas. Insert a `rooting_progress_hook` lambda as the new
  7th argument (before the completion `rooting_hook` lambda). For the CLI, print one concise
  line per pass on completion, de-duplicated so the throttled per-node emit that lands
  exactly on `total` does not print the same "done" line twice:
  ```cpp
  [&, root_last_done = 0](int /*substage_id*/, int substage, int num_substages,
                          int nodes, int total) mutable {
    if (nodes >= total && substage != root_last_done) {
      root_last_done = substage;
      std::cerr << absl::StreamFormat("- rooting pass %d / %d done\n", substage, num_substages);
    }
  },
  ```
- **Incidental `-Werror=return-type` fix.** `build_rough_initial_tree_from_maple` ends with a
  `switch (init_method)` whose cases all `return`, followed by `CHECK(false) << "Unknown init
  method";`. On this toolchain `CHECK(false)` is not treated as `[[noreturn]]`, so recompiling
  the (previously stale) object file trips `-Werror=return-type`. Add a dead-code
  `return Phylo_tree{};` after the `CHECK` to silence it. (Pre-existing latent issue, surfaced
  only because this edit forces a recompile.)

## Part B: WASM API changes (`tools/delphy_wasm.cpp`)

### B1. Init method selection

Introduce a single named constant near the top of the file so switching back is trivial:

```cpp
// The web UI cannot choose an init method, so pin it here.  To revert to the previous
// behaviour, set this to Init_method::old_usher_like (the refined-tree, SPR-refine, and
// rooting progress hooks below simply never fire for old_usher_like).
static constexpr auto k_web_init_method = Init_method::mp_plus_timing;
```

Replace both local `auto init_method = Init_method::old_usher_like;` lines with
`k_web_init_method`.

Note on the fallback: `build_rough_initial_tree_from_maple` dispatches `old_usher_like` to
`build_usher_like_tree`, which only calls `progress_hook` (the guide-tree hook). The other
three hooks are passed but never invoked. So flipping the constant Just Works with no other
change -- the guide-tree progress bar shows, the other three build stages are skipped.

### B2. New hook parameters

In `delphy_parse_fasta_into_initial_tree_async`, replace the single
`int initial_build_progress_hook_id` parameter with four:

```cpp
int guide_tree_progress_hook_id,   // (tips_so_far: int, total_tips: int) -> void
int refined_tree_progress_hook_id, // (round: int, tips_so_far: int, total_tips: int) -> void
int spr_refine_progress_hook_id,   // (attempt: int, max_attempts: int, cur_muts: int) -> void
int rooting_progress_hook_id,      // (substage_id, substage, num_substages, nodes, total): 5x int -> void
```

(`guide_tree_progress_hook_id` is the renamed `initial_build_progress_hook_id`; keeping the
"guide tree" name makes the split self-documenting.)

Add all four hook ids to the thread-pool task's capture list, and wire the build hooks to
`build_rough_initial_tree_from_maple`. Each build-hook lambda also fires the coarse
`stage_progress_hook` once, on its first invocation (see B3) -- so the lambdas below capture
`stage_progress_hook_id` and a per-lambda `fired` latch too. The completion `rooting_hook`
logs the same diagnostics the CLI prints (so capture `L`, the reference-sequence length,
*before* `maple_file` is moved):

```cpp
auto L = std::ssize(maple_file.ref_sequence);  // captured before maple_file is moved
auto tree = new Phylo_tree{
  build_rough_initial_tree_from_maple(
      std::move(maple_file), k_web_init_method, bitgen,
      [stage_progress_hook_id, guide_tree_progress_hook_id, fired = false]
      (int tips_so_far, int total_tips) mutable {
        if (!fired) { fired = true;
          MAIN_THREAD_ASYNC_EM_ASM({delphyRunHook($0, 3);}, stage_progress_hook_id); }  // FASTA stage 3
        MAIN_THREAD_ASYNC_EM_ASM({delphyRunHook($0, $1, $2);},
            guide_tree_progress_hook_id, tips_so_far, total_tips);
      },
      [stage_progress_hook_id, refined_tree_progress_hook_id, fired = false]
      (int round, int tips_so_far, int total_tips) mutable {
        if (!fired) { fired = true;
          MAIN_THREAD_ASYNC_EM_ASM({delphyRunHook($0, 4);}, stage_progress_hook_id); }
        MAIN_THREAD_ASYNC_EM_ASM({delphyRunHook($0, $1, $2, $3);},
            refined_tree_progress_hook_id, round, tips_so_far, total_tips);
      },
      [stage_progress_hook_id, spr_refine_progress_hook_id, fired = false]
      (int attempt, int max_attempts, int cur_muts) mutable {
        if (!fired) { fired = true;
          MAIN_THREAD_ASYNC_EM_ASM({delphyRunHook($0, 5);}, stage_progress_hook_id); }
        MAIN_THREAD_ASYNC_EM_ASM({delphyRunHook($0, $1, $2, $3);},
            spr_refine_progress_hook_id, attempt, max_attempts, cur_muts);
      },
      [stage_progress_hook_id, rooting_progress_hook_id, fired = false]
      (int substage_id, int substage, int num_substages, int nodes, int total) mutable {
        if (!fired) { fired = true;
          MAIN_THREAD_ASYNC_EM_ASM({delphyRunHook($0, 6);}, stage_progress_hook_id); }
        MAIN_THREAD_ASYNC_EM_ASM({delphyRunHook($0, $1, $2, $3, $4, $5);},
            rooting_progress_hook_id, substage_id, substage, num_substages, nodes, total);
      },
      [L](const Rooting_info& rooting_info) {  // completion: log diagnostics (cerr -> JS console)
        auto method_str = (rooting_info.method == Rooting_method::regression)
            ? "root-to-tip regression" : "midpoint";
        std::cerr << absl::StreamFormat("- rooting: %s, R^2 = %.4f", method_str, rooting_info.r2)
                  << absl::StreamFormat(", lambda = %.4g mut/yr (%.2g * 10^-3 mut/site/yr)",
                                        rooting_info.lambda * 365.0,
                                        rooting_info.lambda * 365.0 / L * 1000)
                  << absl::StreamFormat(", t_MRCA = %s\n", to_iso_date(rooting_info.t_MRCA));
      })
};
```

Make the identical change to `delphy_parse_maple_into_initial_tree_async`, except the coarse
stage numbers shift down by one (MAPLE has no analyze stage): guide tree fires stage 2,
refined tree stage 3, SPR refinement stage 4, rooting stage 5. Capture `L` from
`in_maple.ref_sequence` before the move.

### B3. Coarse stage expansion (`stage_progress_hook`)

Break the single "build" coarse stage into four. New numbering:

- FASTA: `{1: read, 2: analyze, 3: guide tree, 4: refined tree, 5: SPR refinement, 6: rooting & timing}`
- MAPLE: `{1: read, 2: guide tree, 3: refined tree, 4: SPR refinement, 5: rooting & timing}`

Emit each `stage_progress_hook` call at the *start* of the corresponding build sub-stage.
As shown in the B2 code, this is done inside each build-hook lambda, guarded by a per-lambda
`fired` latch so it fires exactly once, on that hook's first invocation. It is a best-effort
signal: a sub-stage that never calls its progress hook (e.g. skipped for a tiny tree) simply
does not fire its stage number. Since the granular progress hooks already carry each stage's
identity and label, this coarse signal is redundant for the label but gives the UI a single
clean "entering stage N" event -- handy for resetting the progress bar or showing a stage
title before the first fine-grained update arrives.

**Consequences:**
- The web `stageCallback` currently only does `console.log` (`uploadui.ts:100`), so expanding
  the numbering has **no functional UI breakage** today -- nothing switches on the value.
- The one real coupling cost is that the FASTA and MAPLE numberings differ (they already do
  today), so any UI code that switches on the raw stage number must know the file format.
  Recommendation: the UI should **not** switch on the raw number; use it only to reset the
  bar / set a title, driven by the granular hooks for actual content.
- Net: low risk, modest polish. Included here because it pairs naturally with resetting the
  bar between stages, but it is the most droppable part of this plan if we want to minimize
  surface area.

### B4. Update the header comments

Update the block comments above both functions (currently at `:113`-`:133` and `:238`-`:255`)
to describe the four build sub-stages and the new coarse stage numbering instead of the
single "building the initial tree, one tip at a time" stage.

## Part C: delphy-web changes (`~/github/fathominfo/delphy-web`)

### C1. `src/ts/pythia/delphy_api.ts`

For `parseFastaIntoInitialTreeAsync` and `parseMapleIntoInitialTreeAsync`:

- Replace the `initialBuildProgressHook` parameter with four:
  ```ts
  guideTreeProgressHook: (tipsSoFar: number, totalTips: number) => void = () => void 0,
  refinedTreeProgressHook: (round: number, tipsSoFar: number, totalTips: number) => void = () => void 0,
  sprRefineProgressHook: (attempt: number, maxAttempts: number, curMuts: number) => void = () => void 0,
  rootingProgressHook: (substageId: number, substage: number, numSubstages: number,
                        nodes: number, total: number) => void = () => void 0,
  ```
- Register each with `withHookAsync` and pass the four IDs (in place of the single
  `initialBuildProgressHookId`) to the raw call.
- Update the raw `delphyCoreRaw` type declarations at `:486` and `:496` to list the four new
  hook IDs. `wrapRawAsyncApi` is generic over positional args, so `:818`-`:819` need no
  change.

### C2. `src/ts/pythia/pythia.ts`

`initRunFromFasta` (`:227`) and `initRunFromMaple` (`:250`): replace the single
`initTreeProgressCallback` parameter with `guideTreeProgressCallback`,
`refinedTreeProgressCallback`, `sprRefineProgressCallback`, and `rootingProgressCallback`,
and forward them to the `delphy` calls.

### C3. `src/ts/ui/uploadui.ts`

Replace `initTreeProgressCallback` (`:115`) with four callbacks that drive the same progress
bar with distinct labels. The rooting callback owns the `substage_id -> label` mapping and
computes a global fraction across passes:

```ts
guideTreeProgressCallback = (tipsSoFar: number, totalTips: number) => {
  showProgress(`Building guide tree${warningsLabelAddendum()}`, totalTips, tipsSoFar);
},
refinedTreeProgressCallback = (round: number, tipsSoFar: number, totalTips: number) => {
  showProgress(`Refining guide tree (round ${round})${warningsLabelAddendum()}`, totalTips, tipsSoFar);
},
sprRefineProgressCallback = (attempt: number, maxAttempts: number, curMuts: number) => {
  showProgress(`Optimizing tree: ${curMuts} mutations${warningsLabelAddendum()}`, maxAttempts, attempt);
},
// Keep this map in sync with the Rooting_substage enum in core/utree.h.
rootingSubstageLabels: {[id: number]: string} = {
  1: "bottom-up timing",
  2: "top-down timing",
  3: "root candidate evaluation"
},
rootingProgressCallback = (substageId: number, substage: number, numSubstages: number,
                           nodes: number, total: number) => {
  const what = rootingSubstageLabels[substageId] ?? "rooting and timing";
  // Global fraction across all passes so the bar advances monotonically (scaled to `total`
  // so showProgress renders it as a percentage).
  const soFar = Math.round(((substage - 1) + (total > 0 ? nodes / total : 0)) / numSubstages * total);
  showProgress(`Rooting and timing: ${what}${warningsLabelAddendum()}`, total, soFar);
},
```

These callbacks are declared as entries in the existing `const` comma-chain that already
holds `stageCallback`, `parseProgressCallback`, etc. (so `rootingSubstageLabels` is another
entry in that chain, not a separate top-level `const`).

Notes:
- For SPR refine, the bar fill tracks `attempt / maxAttempts` while the human-readable
  `curMuts` count goes in the label (mirrors the CLI's `num_muts` reporting).
- For rooting, `substageId` selects the label; the bar shows the *global* fraction
  `((substage-1) + nodes/total) / numSubstages` so it advances smoothly across the three
  passes rather than resetting each pass. (The exact `showProgress` argument shaping is a
  detail; the point is one monotonic bar plus a per-pass label.)
- `showProgress` hides the bar when `soFar === 0 || soFar === total`. That's fine here: each
  stage starts near 0 and the bar re-activates as the next stage's first callback fires with
  `soFar > 0`.
- Update all five call sites that currently pass `initTreeProgressCallback` (three
  `initRunFromFasta` calls, two `initRunFromMaple` calls) to pass the four new callbacks
  positionally in their place. (They could instead be bundled into a small object to avoid
  the repetition, but that would change the `pythia.ts` signatures; the positional form keeps
  the change local to the call sites.)
- `stageCallback` is left as the existing `console.log`. Extending it to reset the bar / set a
  stage title on each coarse stage entry (B3) is possible follow-up polish, not done here.

## Design decisions and open questions

1. **Rooting progress uses a generic 5-int hook**, `(substage_id, substage, num_substages,
   nodes_in_substage, total_in_substage)`, so future rooting/timing algorithms with different
   passes need no API change -- only new `Rooting_substage` ids and UI labels. `substage_id`
   is a stable enum owned by the core; the human-readable text lives entirely in the UI.

2. **Coarse `stage_progress_hook` expanded** (Part B3) to guide/refined/SPR/rooting. This is
   low-risk today because the web `stageCallback` only logs, but it is the most droppable
   piece if we want to shrink the change: the granular hooks already identify each stage.

3. **Completion `rooting_hook` logs `Rooting_info` to the console.** Both WASM entry points
   wire the completion `rooting_hook` to a `std::cerr << absl::StreamFormat(...)` line
   mirroring the CLI's rooting log: rooting method, R^2, rate (mut/yr and 10^-3
   mut/site/yr), and t_MRCA as an ISO date. Emscripten forwards `std::cerr` to the JS
   console, so no JS-side hook id or `EM_ASM` embedding is needed (see B2).

4. **Fallback ergonomics.** Reverting the web to `old_usher_like` is a one-line flip of
   `k_web_init_method`; the extra hooks stay wired but silent (only the guide-tree hook fires
   for `old_usher_like`). No conditional code paths.

5. **CLI default flip touches deprecated flags.** `--v0-init-heuristic` must now set
   `old_usher_like` explicitly (A1); miss this and the deprecated flag silently becomes a
   no-op that yields `mp-plus-timing`.

## Files to modify

| Repo | File | Change |
|------|------|--------|
| delphy | `core/cmdline.cpp` | Default `--v0-init` -> `mp_plus_timing`; `--v0-init-heuristic` branch; option help text; forward `rooting_progress_hook` (signature + CLI caller lambda); `return Phylo_tree{}` after the init-method switch |
| delphy | `core/utree.h` | `Rooting_substage` enum; `progress_hook` param on `ols_regression_root_utree`; `rooting_progress_hook` param on `build_initial_phylo_tree` |
| delphy | `core/cmdline.h` | `rooting_progress_hook` param on `build_rough_initial_tree_from_maple` |
| delphy | `core/utree.cpp` | Emit rooting progress in the 3 passes; forward hook through `build_initial_phylo_tree` |
| delphy | `tools/delphy_wasm.cpp` | `k_web_init_method`; split build hook into 4 (+ coarse stage fire, + completion log); comments |
| delphy-web | `src/ts/pythia/delphy_api.ts` | 4 hook params + raw type decls for both parse methods |
| delphy-web | `src/ts/pythia/pythia.ts` | 4 callbacks in `initRunFromFasta`/`initRunFromMaple` |
| delphy-web | `src/ts/ui/uploadui.ts` | 4 progress callbacks (+ rooting substage labels); update 5 call sites |
| delphy-web | `src/ts/delphy/delphy.{js,wasm}` | Regenerated WASM bundle (release build) carrying the new entry-point signatures |

## Verification

1. **CLI default**: `delphy --help` shows `--v0-init ... (default: mp-plus-timing)`; running
   with no `--v0-init` uses the pipeline (confirm via `--v0-steps 0` num_muts matching an
   explicit `--v0-init mp-plus-timing`); `--v0-init-heuristic` and `--v0-init old-usher-like`
   still produce the old behaviour; `--v0-init random` still works.
2. **Rooting hook unit-level check**: with a mid-size dataset, run the pipeline (CLI or a
   test) and confirm the rooting progress hook fires for all three substage ids with
   monotonically advancing `nodes_in_substage`, ending at `(total, total)` per pass.
3. **Build WASM** per the delphy-web build instructions (the release preset
   `conan-emscripten-release`), then copy `build/wasm-release/delphy.js` and `delphy.wasm`
   into `delphy-web/src/ts/delphy/`. Type-check/build the web app (`npx tsc --noEmit`,
   `npm run build`) and run `npm run dev`.
4. **Load a FASTA and a MAPLE dataset** (e.g. from `~/github/broadinstitute/delphy-2026-paper-data`, such as
   `ebola-gire-2014/delphy_inputs/*.maple`) and confirm the progress bar walks through:
   read -> (FASTA: analyze) -> "Building guide tree" -> "Refining guide tree (round k)" ->
   "Optimizing tree: N mutations" -> "Rooting and timing: <pass>" -> run starts.
5. **Confirm the mutation count** shown during SPR refinement decreases over attempts and
   ends near the CLI's step-0 `num_muts` for the same dataset.
6. **Fallback check** (optional): flip `k_web_init_method` to `old_usher_like`, rebuild, and
   confirm only the "Building guide tree" label shows and loading still completes.
