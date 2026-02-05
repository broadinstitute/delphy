# AGENTS.md - Guide for AI Coding Agents

## CRITICAL: Read CONTRIBUTING.md First

**Before writing any code, you MUST read [CONTRIBUTING.md](CONTRIBUTING.md) in full.**

It defines the coding style for this project. Failure to follow these conventions will result in rejected code.

---

## Key Coding Conventions

These are the most important conventions from CONTRIBUTING.md that you MUST follow:

| Convention | Correct | Incorrect |
|------------|---------|-----------|
| Variable/function names | `snake_case` | `camelCase` |
| Type names (classes, structs, aliases) | `Stroustrup_case` | `snake_case` |
| Boolean negation | `not x` | `!x` |
| Natural logarithm | `log(x)` | `ln(x)` |
| Const style | `const char* str` | `char const* str` |
| Pointer/reference | `int* ptr`, `int& ref` | `int *ptr`, `int &ref` |
| Container size | `std::ssize(c)` | `c.size()` |
| Hash maps/sets | `absl::flat_hash_map` | `std::unordered_map` |
| Data containers | structs with public fields | classes with getters/setters |
| Error handling | exceptions | error codes |
| Function style (AAA) | `auto func(int x) -> int;` | `int func(int x);` |

**Compiler flags:** Code must compile cleanly with `-Wall -Wextra -Werror`.

---

## Repository Structure

```
delphy/
├── core/           # Core library (phylogenetic algorithms, ~20K lines)
│   ├── tree.h          # Generic tree data structures
│   ├── phylo_tree.h    # Phylogenetic tree with mutations and times
│   ├── run.h           # MCMC run orchestration
│   ├── spr_move.h      # Subtree Prune and Regraft moves
│   ├── mutations.h     # Mutation and Missation types
│   ├── sequence.h      # DNA sequence types
│   └── api.fbs         # FlatBuffers schema for serialization
├── tests/          # Google Test suite (27 test files)
├── tools/          # Entry point programs
│   ├── delphy.cpp      # CLI tool
│   ├── delphy_ui.cpp   # OpenGL/SDL visualization UI
│   ├── delphy_mcc.cpp  # MCC tree computation utility
│   └── delphy_wasm.cpp # WebAssembly build
├── third-party/    # Git submodules (abseil-cpp, flatbuffers, ctpl, cppcoro, cxxopts)
├── doc/            # Documentation (dphy_file_format.md)
└── tutorials/      # Jupyter notebooks
```

---

## Build System

### Prerequisites
- CMake 3.17.5+ (Ubuntu 22.04 needs manual upgrade - see INSTALL.md)
- Python 3 with venv
- Conan 2.8.1
- GCC or Clang with C++20 support

### Quick Build

```bash
# First-time setup
python3 -m venv delphy-venv
source delphy-venv/bin/activate
pip3 install 'conan==2.8.1'
conan profile detect
git submodule update --init

# Debug build
conan install . --output-folder=build/debug --build=missing --settings=build_type=Debug
cmake --preset conan-debug
cmake --build --preset conan-debug

# Release build
conan install . --output-folder=build/release --build=missing --settings=build_type=Release
cmake --preset conan-release
cmake --build --preset conan-release
```

### Running Tests

```bash
./build/debug/tests/tests              # All tests (debug)
./build/release/tests/tests            # All tests (release)
./build/debug/tests/tests --gtest_filter="Tree_test.*"  # Specific test suite
```

See [INSTALL.md](INSTALL.md) for WebAssembly builds and detailed instructions.

---

## Key Data Structures

### Node Indexing
- Nodes are referenced by `Node_index` (int), not pointers
- `k_no_node = -1` is the sentinel for null/invalid
- Branch `X` connects `parent(X)` to node `X`
- Access nodes via `tree.at(node)` or `tree.at_parent_of(node)`

### Core Types

```cpp
// From sequence.h
Seq_letter      // uint8 bitmask: A=0b0001, C=0b0010, G=0b0100, T=0b1000, N=0b1111
Real_seq_letter // Enum for unambiguous: A, C, G, T only
Site_index      // int alias for genomic position

// From mutations.h
Mutation        // {from, site, to, t} - timed state change
Seq_delta       // {site, from, to} - mutation without time
Missation       // pair<Site_index, Real_seq_letter> - missing data marker

// From tree.h / phylo_tree.h
Node_index      // int alias for node reference
Phylo_node      // Has: name, t, t_min, t_max, mutations, missations
Phylo_tree      // Tree container with ref_sequence
```

### Tree Traversals

```cpp
for (const auto& node : pre_order_traversal(tree)) { ... }
for (const auto& node : post_order_traversal(tree)) { ... }
for (const auto& node : index_order_traversal(tree)) { ... }
```

---

## Common Gotchas

1. **Use `std::ssize()` not `.size()`** - The codebase strictly avoids signed/unsigned mixing.

2. **Use `not` instead of `!`** - Boolean negation uses the keyword, not the operator.

3. **Times are fractional days** - `t` values are days since 2020-01-01.

4. **Reference sequence is the baseline** - `Phylo_tree::ref_sequence` is the root state; mutations are deltas from it.

5. **Branch vs Node** - A branch is identified by its endpoint. Branch X connects parent(X) to X.

6. **Binary trees only** - Most MCMC code assumes exactly 2 children for inner nodes.

7. **Namespace** - All code lives in `namespace delphy { }`.

8. **Tip vs Inner** - Check with `node.is_tip()` or `node.is_inner_node()`.

9. **Header guards** - Use `#ifndef DELPHY_FILENAME_H_` style.

10. **Debug assertions** - Use `assert_tree_integrity()`, `assert_mutation_consistency()` for validation.

---

## Git Conventions

1. **Never amend pushed commits** - Always create fresh commits. Do not use `git commit --amend` or `git rebase` on commits that have already been pushed.

2. **No agent attribution in commits** - Do not add `Co-Authored-By` or similar agent attribution lines in commit messages. However, agents should identify themselves in code reviews, pull request descriptions, and issue comments.

---

## Dependencies

### Via Conan
- boost/1.84.0 (header-only)
- eigen/3.4.0
- ms-gsl/4.0.0

### Via Git Submodules (third-party/)
- abseil-cpp - Hash maps, random numbers, utilities
- flatbuffers - Serialization
- ctpl - Thread pool
- cppcoro - Coroutines/generators
- cxxopts - Command-line parsing

---

## Testing Patterns

Tests use Google Test:

```cpp
#include <gtest/gtest.h>
#include <gmock/gmock.h>

TEST(Tree_test, example) {
    auto tree = Tree<Binary_node>{2};
    EXPECT_EQ(tree.size(), 2);
    EXPECT_THAT(result, testing::ElementsAre(1, 2));
}
```

Test files follow the pattern `{module}_tests.cpp` in the `tests/` directory.
