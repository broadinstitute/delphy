#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "absl/random/mocking_bit_gen.h"
#include "absl/random/mock_distributions.h"
#include "absl/strings/str_format.h"

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

TEST(Distributions_test, k_truncated_poisson_distribution) {
  auto bitgen = std::mt19937_64(12345);
  for (auto lambda : {0.01, 0.1, 1.0, 10.0}) {
    SCOPED_TRACE(absl::StrFormat("Lambda = %g", lambda));
    for (auto min_k : {0, 1, 2, 5, 20}) {
      SCOPED_TRACE(absl::StrFormat("min_k = %d", min_k));
      
      // Draw 100,000 random samples and check count ratios for k = K and k = K+1,
      // where K = max(min_k, floor(lambda))
      auto K = std::max(min_k, static_cast<int>(std::floor(lambda)));
      auto count_K = 0;
      auto count_K_p_1 = 0;
      auto distr = K_truncated_poisson_distribution{lambda, min_k};
      for (auto sample = 0; sample != 100'000; ++sample) {
        auto k = distr(bitgen);

        EXPECT_THAT(k, testing::Ge(min_k));
        if (k == K) { ++count_K; }
        else if (k == (K + 1)) { ++count_K_p_1; }
      }

      SCOPED_TRACE(absl::StrFormat("K = %d, count_K = %d, count_K_p_1 = %d", K, count_K, count_K_p_1));
      
      EXPECT_THAT(count_K, testing::Ge(10));
      EXPECT_THAT(count_K_p_1, testing::Ge(10));
      
      //auto actual_ratio = count_K / double(count_K_p_1);
      auto upper_bound_ratio =
          (count_K + 3*std::sqrt(count_K)) / double(count_K_p_1 - 3*std::sqrt(count_K_p_1));
      auto lower_bound_ratio =
          (count_K - 3*std::sqrt(count_K)) / double(count_K_p_1 + 3*std::sqrt(count_K_p_1));
      auto expected_ratio = 1.0 / (lambda / (K + 1));
      
      EXPECT_THAT(lower_bound_ratio, testing::Lt(expected_ratio));
      EXPECT_THAT(upper_bound_ratio, testing::Gt(expected_ratio));
    }
  }
}

}  // namespace delphy
