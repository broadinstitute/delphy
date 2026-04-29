# Deterministic Hash Containers

## Goal

Currently, two identically configured Delphy runs with identical seeds produce
different results.  After implementing this plan, they will produce byte-identical
results (logs and trees), both single-threaded and multi-threaded.

Reproducing results across machines also requires using the same number of
threads, since different thread counts produce different (but equally valid) subrun
partitions.  The number of threads defaults to the number of CPUs, which varies
machine to machine.  Wherever we output the seed for later reproducibility, we also
output the number of threads used.

Reproducibility does not extend across different platforms (x86-64 vs arm64 vs
WASM), compilers (gcc vs clang), or optimization levels (debug vs release), since
floating-point operations may be reordered and accumulate different roundoff
errors.  For the same Delphy binary on the same platform, however, runs should be
reproducible across machines.

## Root Cause

Delphy uses Abseil's `flat_hash_map` and `flat_hash_set` as its primary hash
containers (including arena-allocated variants like `Scratch_flat_hash_map`).
Abseil is vendored as a git submodule at `third-party/abseil-cpp`.

These containers (internally called "Swiss tables") have three sources of
nondeterminism that cause iteration order to vary across runs:

1. **`absl::Hash::Seed()`** (`hash.h`): returns `reinterpret_cast<uintptr_t>(&kSeed)`,
   which varies with ASLR (Address Space Layout Randomization — an OS security
   feature that loads programs at a different base address on each execution, so
   the address of any static variable differs from run to run).  Every
   `absl::Hash`-based container inherits this nondeterministic seed.

2. **`PerTableSalt()`** (`raw_hash_set.h`): returns
   `reinterpret_cast<uintptr_t>(ctrl) >> 12`, mixing the table's heap address into
   probe sequences that determine element placement.  Varies across runs due to ASLR
   and allocator state.

3. **`RandomSeed()`** (`raw_hash_set.cc`): returns an incrementing counter XOR'd with
   the address of a static variable.  Used by `ShouldInsertBackwardsForDebug` and
   `ShouldRehashForBugDetection` (the latter active in debug/ASan builds via
   `ABSL_SWISSTABLE_ENABLE_GENERATIONS`).  Both the counter and the address are
   nondeterministic: the counter causes iteration order to depend on how many tables
   have been created previously, and the address varies with ASLR.

Since many code paths iterate over these containers while consuming RNG draws (see
Affected Code Locations below), nondeterministic iteration order causes the MCMC
trajectory to diverge across runs.

Abseil officially rejected deterministic iteration
([issue #339](https://github.com/abseil/abseil-cpp/issues/339),
[issue #720](https://github.com/abseil/abseil-cpp/issues/720)).

## Approach: Build-Time Abseil Patch

We patch the vendored Abseil at CMake configure time.  The patch file lives at
`patches/abseil-deterministic-hash-iteration.patch` and neutralizes each of the
three nondeterminism sources listed above:

1. **`absl::Hash::Seed()`** (called `MixingHashState::Seed()` internally) →
   `return 0x243F6A8885A308D3ULL;` (first 64 bits of pi's fractional part — an
   arbitrary fixed constant)
2. **`PerTableSalt()`** → `return 0;`
3. **`RandomSeed()`** → `return 1;`

The CMake block (in the root `CMakeLists.txt`, before `add_subdirectory(abseil-cpp)`)
applies the patch idempotently:

- `git apply --check --reverse` detects if already applied (e.g., on reconfigure)
- If not applied, `git apply` applies it
- If the patch fails (e.g., Abseil upgraded), the build fails with `FATAL_ERROR`

This ensures anyone checking out the repo and initializing submodules gets the patch
automatically.

Because the patch modifies files inside the `third-party/abseil-cpp` submodule,
`git status` would otherwise show it as having "modified content".  We suppress this
by setting `ignore = dirty` on that submodule in `.gitmodules`.  (`.gitmodules` uses
git-config format, which silently drops comments on rewrite, so this plan is the
canonical explanation for that setting.)

## Determinism Test

`tests/deterministic_abseil_containers_tests.cpp` creates two separate instances each
of `absl::flat_hash_set<int>`, `absl::flat_hash_map<int,int>`, and
`Scratch_flat_hash_set<int>` (Delphy's arena-allocated variant of `flat_hash_set`),
each with 20 elements, and asserts that both instances iterate in the same order.
Without the patch, `PerTableSalt` (derived from heap address) and `RandomSeed`
(incrementing counter) would differ between the two instances, producing different
iteration orders.

We verify determinism by comparing two containers rather than hard-coding expected
iteration orders because debug and release builds produce different orders (debug
builds enable `ABSL_SWISSTABLE_ENABLE_GENERATIONS`, which changes the internal
table layout).

These tests catch:

- Missing or outdated Abseil patch
- Regression from an Abseil upgrade

The tests must pass both in isolation and in the full test suite, confirming no
cross-test interference.

## Alternative Considered and Rejected

**`Stable_hash<T>` functor** — a wrapper around `std::hash<T>` threaded through all
container type aliases and direct instantiations.  We ran an experiment to see if
replacing `absl::Hash` with a deterministic hash functor was sufficient, and
discovered the two other sources of nondeterminism described above (`PerTableSalt`
and `RandomSeed`).  The functor approach also requires touching every container
declaration and is easy to miss in new code.  The Abseil patch fixes all three
sources at once with zero changes to Delphy source.

## Existing Workaround to Remove

`build_usher_like_tree()` in `core/phylo_tree.cpp` currently copies `Site_deltas`
(a hash map) into a sorted vector before a loop that draws RNG values per entry.
This was a targeted workaround for the `absl::Hash::Seed()` nondeterminism.  With
the Abseil determinism patch in place, iteration order is deterministic and the
sorting step is unnecessary.  Remove it and iterate directly over the hash map.

## Seed and Thread Count in Log File Headers

Every BEAST-format log file and trees file starts with comment headers produced by
`stamp_version_into_log_file` (in `core/io.cpp`).  We add two lines:

```
# Seed: 42
# Parallelism: 8
```

The seed is stored in the `Run` object as `prng_seed_` (passed to the constructor,
read via `run.prng_seed()`).  The thread count comes from the existing
`run.num_parts()`, which is set to the thread pool size before output begins.

We route both values through the `Run` object rather than global variables so that
multiple `Delphy_context` instances (e.g., in the WASM interface) can coexist without
clobbering shared state.

The `std::cerr` Parallelism line in `cli_main_loop()` (which goes to the console,
not to log files) was also changed from `std::cout` to `std::cerr` to match the
convention used for all other informational output.

## delphy-web Seed Handling

The web frontend (`delphy-web`) generates a random seed at startup and passes it to
the WASM core for both initial tree construction and MCMC runs.

`constants.ts` defines:

```typescript
const RANDOM_SEED_OVERRIDE = 0;

export const RANDOM_SEED = (RANDOM_SEED_OVERRIDE as number) !== 0
  ? RANDOM_SEED_OVERRIDE
  : (crypto.getRandomValues(new Uint32Array(1))[0] || 1);
```

When `RANDOM_SEED_OVERRIDE` is `0` (the default), a cryptographically random 32-bit
seed is generated.  Setting it to a nonzero value forces that exact seed, enabling
reproduction of a specific run.  The `as number` cast prevents TypeScript from
narrowing the literal type when the override is nonzero (which would make
`!== 0` a type error).

The seed is passed at two points:

1. **Context creation**: `new Delphy(RANDOM_SEED)` — seeds the `bitgen_` used to build
   the initial tree.
2. **Run creation**: `createRun(phyloTree, RANDOM_SEED)` — seeds the MCMC PRNG via the
   `Run` constructor.

At startup, `pythia.ts` logs the seed to the console:

```
Delphy core loaded (version 1.3.3, build 42, commit abc1234, seed 786310996)
```

This way, any bug report can include the seed, and a developer can set
`RANDOM_SEED_OVERRIDE` to that value to reproduce the run.

## Files Changed

### Abseil determinism patch

| File | Change |
|------|--------|
| `patches/abseil-deterministic-hash-iteration.patch` | New: the 3-site Abseil patch |
| `CMakeLists.txt` | Patch application block before `add_subdirectory(abseil-cpp)` |
| `tests/deterministic_abseil_containers_tests.cpp` | New: 3 tests verifying deterministic iteration order |
| `tests/CMakeLists.txt` | Add `deterministic_abseil_containers_tests.cpp` to test sources |

### Workaround removal

| File | Change |
|------|--------|
| `core/phylo_tree.cpp` | Remove sort-before-RNG workaround (see above) |

### Seed and thread count in log file headers

| File | Change |
|------|--------|
| `core/run.h` | Add `prng_seed_` member and `prng_seed()` getter; add `prng_seed` parameter to constructor |
| `core/run.cpp` | Accept and store `prng_seed` in constructor |
| `core/io.h` | `stamp_version_into_log_file` now holds a `const Run&` instead of being an empty struct |
| `core/io.cpp` | `operator<<` reads seed from `run.prng_seed()` and thread count from `run.num_parts()` |
| `core/beasty_output.cpp` | Pass `run` to `stamp_version_into_log_file{run}` (3 call sites) |
| `core/beasty_input.cpp` | Pass `run` to `stamp_version_into_log_file{run}` (3 call sites) |
| `core/cmdline.cpp` | Pass `seed` to `Run` constructor |
| `tools/delphy_wasm.cpp` | Pass `real_seed` to `Run` constructor in `delphy_create_run()` |
| `tools/delphy.cpp` | Change `std::cout` to `std::cerr` for the Parallelism line in `cli_main_loop()` |

### delphy-web seed plumbing (separate repo: `delphy-web`)

| File | Change |
|------|--------|
| `src/ts/constants.ts` | Add `RANDOM_SEED_OVERRIDE` / `RANDOM_SEED` constants (see above) |
| `src/ts/pythia/delphy_api.ts` | `Delphy` constructor and `createRun` each take a required `prngSeed: number` parameter (no default) |
| `src/ts/pythia/pythia.ts` | Import `RANDOM_SEED` from constants; pass it to `new Delphy(RANDOM_SEED)` and all 3 `createRun` calls; log seed at startup |

## Verification

Run Delphy twice with seed 42 and compare log + tree output (ignoring the
invocation-line header which contains output filenames).

| Dataset | Data dir | Pop Model | Threads | Steps |
|---------|----------|-----------|---------|-------|
| Ebola (Dudas 2017) | `ebola-dudas-2017` | Skygrid | 1 | 20,000 |
| Ebola (Dudas 2017) | `ebola-dudas-2017` | Skygrid | all cores | 20,000 |
| SARS-CoV-2 (Lemieux) | `sars-cov-2-lemieux` | Exponential | 1 | 20,000 |
| SARS-CoV-2 (Lemieux) | `sars-cov-2-lemieux` | Exponential | all cores | 20,000 |

Data dirs are under `~/github/broadinstitute/delphy-2026-paper-data/`.  See the
`02a_run_delphy.sh` scripts in each directory for the full set of run options.

All four pairs must produce byte-identical output.

## Affected Code Locations

The following locations iterate hash containers in ways that affect the RNG stream or
mutation ordering.  The key Delphy types involved are `Site_deltas` (a
`Scratch_flat_hash_map<Site_index, Site_delta>` mapping sites to state
changes), `Node_set` (a `flat_hash_set` of tree node indices), and `Site_index` /
`Node_index` (integer indices into the sequence and tree, respectively).

Each location produces a valid MCMC trajectory regardless of iteration order, but
different orders lead to different trajectories.  This inventory confirms that the
Abseil determinism patch covers all of them:

### RNG-consuming iterations

1. **`spr_move.cpp:1181`** — `sample_mutations_on_branch()`: iterates `Site_deltas` and
   draws RNG per site.  This is the hottest path.

2. **`phylo_tree.cpp:555`** — `randomize_branch_mutation_times()`: iterates
   `Scratch_flat_hash_set<Site_index>` in the `complicated` (multi-mutation) path,
   drawing RNG per site.  Rare but real.

3. **`phylo_tree.cpp:478`** — `rewire_tree_through_sequential_accretion()`: iterates
   `Node_set` to assign inner nodes.  Only during initial tree construction.

### Mutation-list construction

These build mutation lists by iterating `Site_deltas`.  The lists are typically sorted
afterward, but intermediate order could theoretically affect floating-point accumulation:

4. **`spr_move.cpp:410, 422, 534, 928, 937, 1054`** — `peel_rooty_graft()`,
   `apply_rooty_graft()`, `peel_inner_graft()`, `apply_inner_graft()`
5. **`tree_editing.cpp:150, 298`** — `change_branch_point_time()`, `hop_up()`
6. **`phylo_tree.cpp:277`** — `add_tree_from_sequence_vector()`
7. **`run.cpp:147`** — `update_seq_above_subroot()`

## Scope Boundaries

- Delphy source changes: removing the sort workaround in `phylo_tree.cpp`, adding
  `prng_seed_` to the `Run` class, plumbing seed and thread count into log file
  headers via `stamp_version_into_log_file`, and fixing `stdout` → `stderr` for an
  informational console line.
- delphy-web changes (separate repo): plumbing `RANDOM_SEED` through to context and
  run creation, adding the override mechanism in `constants.ts`, and logging the seed
  at startup.
- Everything else is in the Abseil patch, the CMake build, and tests.
- The `Stable_hash<T>` functor approach was explored and abandoned in favor of patching.
- Multi-threaded determinism works because subruns use independent RNGs seeded
  deterministically from the main RNG; thread scheduling does not affect the RNG stream
  within each subrun.
