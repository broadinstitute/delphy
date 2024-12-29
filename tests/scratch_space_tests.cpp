#include "scratch_space.h"

#include "gtest/gtest.h"

namespace delphy {

TEST(Scratch_space_test, vector) {
  auto max_n = 100;

  for (auto iter = 0; iter != 10; ++iter) {
    auto scope = Local_arena_scope{};
    auto xs = Scratch_vector<int>{};
    for (auto i = 0; i != max_n; ++i) {
      xs.push_back(i);
    }
    auto sum = 0;
    for (auto x : xs) {
      sum += x;
    }
    EXPECT_EQ(sum, ((max_n - 1) * max_n) / 2);
  }
}

TEST(Scratch_space_test, unordered_set) {
  auto max_n = 100;

  for (auto iter = 0; iter != 10; ++iter) {
    auto scope = Local_arena_scope{};
    auto xs = Scratch_flat_hash_set<int>{};
    for (auto i = 0; i != max_n; ++i) {
      xs.insert((1337 * i) % (2 * max_n));
    }
    auto sum = 0;
    for (auto x : xs) {
      sum += x;
    }
    EXPECT_GT(sum, 0);
  }
}

TEST(Scratch_space_test, unordered_map) {
  auto max_n = 100;

  for (auto iter = 0; iter != 10; ++iter) {
    auto scope = Local_arena_scope{};
    auto m = Scratch_flat_hash_map<int, int>{};
    for (auto i = 0; i != max_n; ++i) {
      m[i] = 2 * i;
    }
    auto sum = 0;
    for (const auto [key, val] : m) {
      sum += val;
    }
    EXPECT_EQ(sum, 2 * ((max_n - 1) * max_n) / 2);
  }
}

}  // namespace delphy
