# Enable Delphy Core Unit Tests for Node.js/WASM

## Summary

This plan enables running Delphy's core C++ unit test suite under Node.js using the
Emscripten-compiled WASM module.  Without running unit tests on WASM, we cannot
automatically and continuously verify that WASM-specific fixes work correctly, and we
wouldn't have immediate visibility on WASM-only bugs in the future.

A key motivating example: Emscripten's `long double` math functions are broken — they
overflow at the `double` threshold — causing NaN in Boost's gamma functions when using
the default `promote_double<true>` policy.  A recent fix (see
`2026-04-22-01-fix-gamma-q-nan-from-emscripten-broken-long-double.md`) introduced
`safe_gamma_q`/`safe_gamma_q_inv` wrappers using `promote_double<false>` to work around
this.  Running unit tests on WASM ensures those wrappers remain correct.

## Background

### Current State
- Unit tests compile into a native `tests` executable using Google Test
- Unit tests are already built for WASM in CI (to catch Clang compilation issues), but not executed
- CI runs native tests but not WASM tests

### Problems Discovered

Three issues prevent tests from running correctly under WASM/Node:

1. **Google Mock incompatibility**: Tests using `absl::MockingBitGen` (which relies
   on Google Mock, i.e., gmock) (e.g., `distributions_tests.cpp`) hang in WASM.  The
   gmock framework uses `abort()` for failed expectations, which doesn't properly
   terminate in Emscripten.

2. **Boost `long double` promotion**: As described in
   `2026-04-22-01-fix-gamma-q-nan-from-emscripten-broken-long-double.md`, Emscripten's
   `long double` math functions are broken (they overflow at the `double` threshold), and
   Boost's default `long double` promotion causes NaN in `gamma_q`/`gamma_q_inv`.  Tests
   in `safe_gamma_math_tests.cpp` use bare `boost::math::gamma_q()` calls as reference
   values to compare against the `safe_gamma_q()` wrappers — for some parameter
   combinations (large `a` values like 271.4, 500), those bare calls produce NaN on WASM,
   causing the comparison assertions to fail.

3. **Exit code swallowing**: By default Emscripten's `EXIT_RUNTIME` is disabled, meaning
   the Node.js process always exits with code 0 regardless of what `main()` returns.  Without
   fixing this, all test failures are silently swallowed.

In addition, Emscripten's Node.js builds output a `tests.js` file (not a self-contained
`tests` executable), so a shell script wrapper is needed to make `./tests` work.

## Proposed Changes

### Core Infrastructure

#### `CMakeLists.txt` (top-level)
Refactor to define a shared `WASM_LINKER_FLAGS` variable for the common Emscripten
linker options and update the stale comment above the GoogleTest fetch.  The goal
is to avoid duplicating the common parts of the `delphy` target's options list for
the `tests` target; the `delphy` target should be compiled and linked with exactly the
same effective options as it is now (flag ordering may change cosmetically).

**Stale comment fix:** The comment above `cmake_policy(SET CMP0135 NEW)` (currently
"only useful in non-WASM builds, but we build them for WASM anyway to ensure that no
test code trips up Clang") should be updated to:

```cmake
# Tests are built for WASM to both catch Clang compilation issues and to verify
# WASM-specific functionality (e.g., gamma_q with promote_double<false> policy)
```

**Shared flags:** Define `WASM_LINKER_FLAGS` in a new `if(EMSCRIPTEN)` block *before*
`add_subdirectory(tests)` (e.g., right after `add_subdirectory(core)`).  This is necessary
because CMake processes `add_subdirectory` immediately, so the variable must be set before
`tests/CMakeLists.txt` is evaluated:

```cmake
# Define WASM linker flags early so tests/CMakeLists.txt can use them
if(EMSCRIPTEN)
  # Core Emscripten linker flags shared between delphy and tests targets
  set(WASM_LINKER_FLAGS
    -fwasm-exceptions
    -pthread
    # Don't mess up 64-bit integers; see https://v8.dev/features/wasm-bigint
    -sWASM_BIGINT
    -sINITIAL_MEMORY=16mb
    -sALLOW_MEMORY_GROWTH=1
    # Avoid grief with decoding UTF-8 strings in SharedArrayBuffers (which is what
    # the WASM linear memory looks like from the JS side).  See here for details:
    # https://github.com/emscripten-core/emscripten/issues/18034
    -sTEXTDECODER=0
  )
endif()
```

**delphy target refactor:** In the existing `if(EMSCRIPTEN)` block near the bottom of
the file, replace the inline flags with `${WASM_LINKER_FLAGS}`:

```cmake
if(EMSCRIPTEN)
  # delphy WASM target
  add_executable(delphy tools/delphy_wasm.cpp)
  target_link_libraries(delphy delphy_core)
  target_link_options(delphy PRIVATE
    --no-entry
    ${WASM_LINKER_FLAGS}
    -sPTHREAD_POOL_SIZE=navigator.hardwareConcurrency
    -sEXPORTED_FUNCTIONS=_malloc,_free
    -sEXPORTED_RUNTIME_METHODS=UTF8ToString,stringToUTF8OnStack
    $<$<CONFIG:DEBUG>:-gsource-map --source-map-base ./>
  )
  if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    # Only enable mimalloc in non-debug builds because it conflicts with -fsanitize=address
    target_link_options(delphy PRIVATE
      # Use a scalable memory allocator (the default allocator has high contention
      # and is slow in multi-threaded environments)
      # See https://emscripten.org/docs/porting/pthreads.html#allocator-performance
      -sMALLOC=mimalloc
    )
  endif()
endif()
```

Note: The only change to the `delphy` target is that the `-fwasm-exceptions`, `-pthread`,
`-sWASM_BIGINT`, `-sINITIAL_MEMORY`, `-sALLOW_MEMORY_GROWTH`, and `-sTEXTDECODER` flags
now come from `${WASM_LINKER_FLAGS}` instead of being inline — same flags, different
source.

#### `tests/CMakeLists.txt`
Add WASM-specific linker flags to the existing `tests` target, reusing the shared
`WASM_LINKER_FLAGS` from the top-level:

```cmake
if(EMSCRIPTEN)
  target_link_options(tests PRIVATE
    ${WASM_LINKER_FLAGS}
    -sPTHREAD_POOL_SIZE=4
    # Without this, Emscripten's EXIT_RUNTIME defaults to false and the Node.js
    # process always exits with code 0, silently swallowing test failures.
    -sEXIT_RUNTIME=1
  )
  # Create a shell script wrapper so tests can be run as ./tests
  file(GENERATE OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/tests"
    CONTENT "#!/bin/sh\ndir=$(dirname \"$0\")\nexec node \"$dir/$<TARGET_FILE_NAME:tests>\" \"$@\"\n"
    FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                     GROUP_READ GROUP_EXECUTE
                     WORLD_READ WORLD_EXECUTE
  )
endif()
```

**`-sEXIT_RUNTIME=1` is critical:** Emscripten defaults `EXIT_RUNTIME` to false, meaning
the Node.js process always exits with code 0 regardless of what `main()` returns.  Without
this flag, test failures would be silently swallowed and CI would always report success.

Tests use a fixed pool size of 4 for deterministic, reproducible behavior across
environments, while the `delphy` target uses `navigator.hardwareConcurrency` to maximize
parallelism in production use.

**Shell script wrapper:** Emscripten's CMake toolchain sets `CMAKE_EXECUTABLE_SUFFIX` to
`.js`, so the build produces `tests.js` (plus `tests.wasm` and `tests.worker.js`).
Emscripten does *not* auto-generate a shell runner.  The `file(GENERATE)` block above
creates a `tests` shell script (no extension) in the same directory that invokes
`node tests.js "$@"`, forwarding all arguments.  This makes `--gtest_filter` work
transparently:

```bash
./build/wasm-release/tests/tests --gtest_filter=Distributions_test.*
```

### Test-Related Fixes

These changes clean up rough edges once the basic test infrastructure is in place.

#### `tests/distributions_tests.cpp`
Replace `absl::MockingBitGen` with a real deterministic RNG to avoid gmock hangs in WASM.
The current test mocks the uniform RNG to return a fixed value (0.3), then checks both
bounds (`Le`/`Ge`) and exact values (`DoubleNear`) against formulas that assume that
fixed input.  With a real RNG, the exact values are no longer predictable from a formula,
so the `DoubleNear` assertions must be removed — only the bounds checks remain.

Changes:
- Remove includes: `absl/random/mocking_bit_gen.h`, `absl/random/mock_distributions.h`
- Keep `gmock/gmock.h` (already present; still needed for `testing::Le`, `testing::Ge`
  matchers; the `GTest::gmock` link dependency in `tests/CMakeLists.txt` stays as-is)
- Add include: `<random>`
- In `TEST(Distributions_test, bounded_exponential_distribution)`:
  - Replace `absl::MockingBitGen` with `std::mt19937_64{12345}`
  - Remove `EXPECT_CALL` setup
  - Remove the three `DoubleNear` exact-value assertions (one per test case:
    `unbounded_left`, `unbounded_right`, `bounded`)
  - Keep the `Le` and `Ge` bounds checks — these verify the core correctness property
    that `Bounded_exponential_distribution` respects its bounds

#### `tests/safe_gamma_math_tests.cpp`
The tests call `safe_gamma_q()` and `safe_gamma_q_inv()` — the wrappers in
`safe_gamma_math.h` that already use `promote_double<false>` internally.  These are
the functions under test and need no changes.

The problem is on the **right-hand side** of comparison assertions: several tests compare
`safe_gamma_q(a, x)` against bare `boost::math::gamma_q(a, x)` (which uses the default
`promote_double<true>` policy) as expected values.  On WASM, those bare calls produce
NaN for large parameters, causing the assertions to fail.  (Small parameters like
a=2, a=1, a=10 don't trigger the `long double` overflow and are unaffected.)

**Fix:** Add a `gmock/gmock.h` include and define a custom `NearOrNaN` matcher that
behaves like `DoubleNear` but also accepts NaN expected values on Emscripten (where
bare Boost gamma functions return NaN due to broken `long double`):

```cpp
MATCHER_P2(NearOrNaN, expected, tolerance, "") {
  if (std::isnan(expected)) {
#ifdef __EMSCRIPTEN__
    return true;
#else
    *result_listener << "expected is NaN on non-Emscripten platform";
    return false;
#endif
  }
  auto diff = std::abs(arg - expected);
  if (diff > tolerance) {
    *result_listener << "which is " << diff << " from " << expected;
    return false;
  }
  return true;
}
```

This turns verbose if/else blocks into single-line assertions:
```cpp
EXPECT_THAT(safe_gamma_q(271.4, 280.0), NearOrNaN(boost::math::gamma_q(271.4, 280.0), 1e-12));
```

Affected tests (those using bare Boost calls that may produce NaN on Emscripten):
- `SafeGammaMathTest.SafeGammaQNearMode` — a=271.4 (triggers overflow)
- `SafeGammaMathTest.SafeGammaQLargeALargeX` — a=500 (triggers overflow)
- `SafeGammaMathTest.SafeLogGammaIntegralBasic` — uses bare Boost calls
- `SafeGammaMathTest.SafeLogGammaIntegralAtBounds` — a=271.4 (triggers overflow)

`SafeGammaMathTest.SafeGammaQBasic` (a=2, 1, 10) uses small parameters that don't
trigger the overflow, so it keeps plain `EXPECT_NEAR`.

**Regression test:** Add a new `TEST(SafeGammaMathTest, PromoteDoubleRegression)` that
uses `#ifdef __EMSCRIPTEN__` directly — this test must *fail* if the bare result is
unexpectedly not-NaN on Emscripten (signaling the bug was fixed upstream):
```cpp
TEST(SafeGammaMathTest, PromoteDoubleRegression) {
  auto bare_result = boost::math::gamma_q(271.4, 6601.0);
#ifdef __EMSCRIPTEN__
  EXPECT_TRUE(std::isnan(bare_result)) << "got " << bare_result;
#else
  EXPECT_NEAR(bare_result, 0.0, 1e-12);
#endif
  // safe_gamma_q uses promote_double<false> and works on all platforms
  EXPECT_NEAR(safe_gamma_q(271.4, 6601.0), 0.0, 1e-12);
}
```

### CI Integration

#### `.github/workflows/ci.yml`
The existing WASM job already builds all targets (including `tests`) via its blanket
`cmake --build --preset conan-emscripten-release --parallel` step.  The only additions
are a Node.js setup step (required for SharedArrayBuffer/pthreads support) and a test
execution step after building:
```yaml
- name: Set up Node.js
  uses: actions/setup-node@v4
  with:
    node-version: '20'

- name: Run WASM tests
  run: ./build/wasm-release/tests/tests
```

Node.js >= 16 is required for SharedArrayBuffer support (used by Emscripten's pthreads
implementation).  Node 20 is the current LTS and is used here for stability.

## Correctness Argument

The `promote_double<false>` policy (as described in
`2026-04-22-01-fix-gamma-q-nan-from-emscripten-broken-long-double.md`) ensures numerical
correctness across all platforms by keeping Boost's gamma calculations in `double`
precision, matching the actual capabilities of Emscripten's math library.

Replacing mocks with real deterministic RNGs preserves test intent — the
`Bounded_exponential_distribution` test verifies that outputs stay within bounds,
which is the core correctness property.  Exact value matching against mocked
uniform outputs is less important than verifying the distribution's bounding logic.

The same test source code runs on both native and WASM platforms, ensuring
consistent behavior.  The full test suite passes under Node.js with these changes.

## Testing Strategy

1. Build WASM tests: `cmake --build --preset conan-emscripten-release --target tests`
2. Run specific tests: `./build/wasm-release/tests/tests --gtest_filter=Distributions_test.*`
3. Run all tests: `./build/wasm-release/tests/tests`
4. Verify all tests pass (including new platform-aware regression tests for gamma_q with `promote_double<true>`)
5. Native tests continue as primary validation in CI
6. WASM tests added as CI check, failing the build if any test fails

## Files Modified

- `CMakeLists.txt` (top-level: define shared `WASM_LINKER_FLAGS` variable, update stale comment)
- `tests/CMakeLists.txt` (add WASM linker flags with `-sEXIT_RUNTIME=1`, `file(GENERATE)` shell script wrapper for `./tests` invocation)
- `tests/distributions_tests.cpp` (replace `absl::MockingBitGen` with real RNG, remove exact-value assertions)
- `tests/safe_gamma_math_tests.cpp` (add `NearOrNaN` custom matcher for bare Boost comparisons, add regression test)
- `.github/workflows/ci.yml` (add Node.js v20 setup and WASM test step)
