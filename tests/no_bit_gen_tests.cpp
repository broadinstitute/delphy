#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <absl/random/distributions.h>

#include "no_bit_gen.h"

namespace delphy {

TEST(No_bit_gen_test, simple) {
  auto bitgen = No_bit_gen{};

  EXPECT_THROW(bitgen(), std::runtime_error);
  EXPECT_THROW(absl::Uniform(bitgen, 0.0, 1.0), std::runtime_error);
  EXPECT_THROW(std::bernoulli_distribution{0.05}(bitgen), std::runtime_error);
}

}  // namespace delphy
