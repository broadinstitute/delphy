#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <absl/random/mocking_bit_gen.h>
#include <absl/random/mock_distributions.h>

#include "distributions.h"

namespace delphy {

TEST(Distributions_test, bounded_exponential_distribution) {
  auto bitgen = absl::MockingBitGen{};
  EXPECT_CALL(absl::MockUniform<double>(), Call(absl::IntervalOpenOpen, bitgen, 0.0, 1.0))
      .WillRepeatedly(testing::Return(0.3));

  auto lambda = 2.3;
  auto inf = std::numeric_limits<double>::infinity();
  
  auto unbounded_left = Bounded_exponential_distribution{+lambda, -inf, 5.0}(bitgen);
  EXPECT_THAT(unbounded_left, testing::Le(5.0));
  EXPECT_THAT(unbounded_left, testing::DoubleNear(5.0 + (1/+lambda)*std::log(0.3), 1e-6));

  auto unbounded_right = Bounded_exponential_distribution{-lambda, 3.0, +inf}(bitgen);
  EXPECT_THAT(unbounded_right, testing::Ge(3.0));
  EXPECT_THAT(unbounded_right, testing::DoubleNear(3.0 + (1/-lambda)*std::log(0.3), 1e-6));

  auto bounded = Bounded_exponential_distribution{lambda, 2.0, 5.0}(bitgen);
  EXPECT_THAT(bounded, testing::Ge(2.0));
  EXPECT_THAT(bounded, testing::Le(5.0));
  EXPECT_THAT(bounded, testing::DoubleNear(2.0 + (1/lambda)*std::log1p(0.3*(std::exp(lambda*(5.0-2.0)) - 1)), 1e-6));
}

}  // namespace delphy
