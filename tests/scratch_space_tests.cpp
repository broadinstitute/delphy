#include "scratch_space.h"

#include "gtest/gtest.h"

namespace delphy {

TEST(Scratch_space_test, vector) {
  auto starting_arena_size = size_t{100};
  auto max_n = 100;

  auto scratch_space = Scratch_space{starting_arena_size};
  for (auto iter = 0; iter != 10; ++iter) {
    scratch_space.reset();
    auto xs = Scratch_vector<int>{scratch_space};
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
  auto starting_arena_size = size_t{100};
  auto max_n = 100;

  auto scratch_space = Scratch_space{starting_arena_size};
  for (auto iter = 0; iter != 10; ++iter) {
    scratch_space.reset();
    auto xs = Scratch_flat_hash_set<int>{scratch_space};
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
  auto starting_arena_size = size_t{100};
  auto max_n = 100;

  auto scratch_space = Scratch_space{starting_arena_size};
  for (auto iter = 0; iter != 10; ++iter) {
    scratch_space.reset();
    auto m = Scratch_flat_hash_map<int, int>{scratch_space};
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
