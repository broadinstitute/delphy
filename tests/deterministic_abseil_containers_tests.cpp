#include <vector>

#include <gtest/gtest.h>

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"

#include "scratch_space.h"

namespace delphy {

// These tests verify that absl::flat_hash_map/set iteration order is
// deterministic across runs.  They rely on the Abseil patches applied at
// build time (see patches/abseil-deterministic-hash-iteration.patch).
//
// We verify determinism by constructing two separate containers with
// identical elements and checking that they iterate in the same order.
// Without the patch, PerTableSalt (derived from heap address) and
// RandomSeed (incrementing counter) would differ between the two
// containers, producing different iteration orders.

static auto collect_set_order() -> std::vector<int> {
  auto s = absl::flat_hash_set<int>{};
  for (auto i = 0; i < 20; ++i) {
    s.insert(i * 7 + 3);
  }
  auto result = std::vector<int>{};
  for (auto x : s) {
    result.push_back(x);
  }
  return result;
}

TEST(Deterministic_hash_test, flat_hash_set_default_hasher) {
  auto order1 = collect_set_order();
  auto order2 = collect_set_order();
  EXPECT_EQ(order1, order2)
      << "absl::flat_hash_set<int> iteration order is nondeterministic — "
         "is the Abseil deterministic-hash patch applied?";
}

TEST(Deterministic_hash_test, flat_hash_map_default_hasher) {
  auto collect = []() {
    auto m = absl::flat_hash_map<int, int>{};
    for (auto i = 0; i < 20; ++i) {
      m[i * 7 + 3] = i;
    }
    auto keys = std::vector<int>{};
    for (const auto& [k, v] : m) {
      keys.push_back(k);
    }
    return keys;
  };

  auto order1 = collect();
  auto order2 = collect();
  EXPECT_EQ(order1, order2)
      << "absl::flat_hash_map<int,int> iteration order is nondeterministic — "
         "is the Abseil deterministic-hash patch applied?";
}

TEST(Deterministic_hash_test, scratch_flat_hash_set) {
  auto collect = []() {
    auto scope = Local_arena_scope{};
    auto s = Scratch_flat_hash_set<int>{};
    for (auto i = 0; i < 20; ++i) {
      s.insert(i * 7 + 3);
    }
    auto result = std::vector<int>{};
    for (auto x : s) {
      result.push_back(x);
    }
    return result;
  };

  auto order1 = collect();
  auto order2 = collect();
  EXPECT_EQ(order1, order2)
      << "Scratch_flat_hash_set<int> iteration order is nondeterministic — "
         "is the Abseil deterministic-hash patch applied?";
}

}  // namespace delphy
