// safe_gamma_math_tests.cpp

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <cmath>
#include <limits>
#include <random>

#include "safe_gamma_math.h"

namespace delphy {

// Like testing::DoubleNear, but if `expected` is NaN (which happens on Emscripten for
// bare Boost gamma functions due to broken long double), the match succeeds trivially —
// we can't validate against a broken reference value.  On non-Emscripten platforms,
// a NaN expected value is treated as a test failure.
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

TEST(SafeGammaMathTest, SafeGammaQBasic) {
  // Small a, normal x
  EXPECT_NEAR(safe_gamma_q(2.0, 5.0), boost::math::gamma_q(2.0, 5.0), 1e-12);
  EXPECT_NEAR(safe_gamma_q(1.0, 1.0), boost::math::gamma_q(1.0, 1.0), 1e-12);
  EXPECT_NEAR(safe_gamma_q(10.0, 5.0), boost::math::gamma_q(10.0, 5.0), 1e-12);
}

TEST(SafeGammaMathTest, SafeGammaQLargeX) {
  // Large a, x >> a: should return ~0
  EXPECT_NEAR(safe_gamma_q(271.4, 6601.0), 0.0, 1e-12);
  EXPECT_NEAR(safe_gamma_q(1000.0, 100.0), 1.0, 1e-12);  // x << a: should return ~1
}

TEST(SafeGammaMathTest, SafeGammaQNearMode) {
  // Near mode x ≈ a
  EXPECT_THAT(safe_gamma_q(271.4, 280.0), NearOrNaN(boost::math::gamma_q(271.4, 280.0), 1e-12));
}

TEST(SafeGammaMathTest, SafeGammaQLargeALargeX) {
  // Large a, moderate x
  EXPECT_THAT(safe_gamma_q(500.0, 1000.0), NearOrNaN(boost::math::gamma_q(500.0, 1000.0), 1e-12));
}

TEST(SafeGammaMathTest, SafeGammaQAtZero) {
  // At x = 0, Q(a, 0) = 1 for all a > 0
  for (double a : {1.0, 10.0, 100.0, 1000.0}) {
    EXPECT_NEAR(safe_gamma_q(a, 0.0), 1.0, 1e-12);
  }
}

TEST(SafeGammaMathTest, SafeGammaQInvRoundTrip) {
  // Round-trip: gamma_q_inv(a, gamma_q(a, x)) ≈ x
  std::mt19937 bitgen{42};
  for (double a : {0.5, 1.0, 2.0, 5.0, 10.0, 50.0, 100.0, 271.4, 500.0}) {
    for (double x : {0.1, 1.0, 5.0, 10.0, a/2, a, 2*a, 10*a}) {
      if (x > 0) {
        auto Q = safe_gamma_q(a, x);
        // Only test round-trip if Q is not too close to 0 or 1
        if (Q > 1e-10 && Q < 1.0 - 1e-10) {
          auto x_back = safe_gamma_q_inv(a, Q);
          // Allow small relative error
          EXPECT_NEAR(x_back, x, 1e-6 * (x + 1.0)) << "a=" << a << " x=" << x;
        }
      }
    }
  }
}

TEST(SafeGammaMathTest, SafeGammaQInvExtremeQ) {
  // Extreme Q values — verified via round-trip for precision
  std::mt19937 bitgen{42};
  for (double a : {50.0, 100.0, 271.4, 500.0}) {
    // Very small Q: round-trip should recover Q with reasonable relative error
    auto Q_lo = 1e-300;
    auto x_lo = safe_gamma_q_inv(a, Q_lo);
    auto Q_lo_back = safe_gamma_q(a, x_lo);
    // Allow 1% relative error or absolute error of 1e-310, whichever is larger
    EXPECT_NEAR(Q_lo_back, Q_lo, std::max(0.01 * Q_lo, 1e-310)) << "Round-trip failed for Q=1e-300";

    // Q near 1: round-trip should recover Q within tight absolute tolerance
    auto Q_hi = 1.0 - 1e-15;
    auto x_hi = safe_gamma_q_inv(a, Q_hi);
    auto Q_hi_back = safe_gamma_q(a, x_hi);
    EXPECT_NEAR(Q_hi_back, Q_hi, 1e-12) << "Round-trip failed for Q near 1";

    // Boundary cases
    auto x_zero = safe_gamma_q_inv(a, 0.0);
    EXPECT_TRUE(std::isinf(x_zero) && x_zero > 0) << "Q=0 should return +inf";
    auto x_one = safe_gamma_q_inv(a, 1.0);
    EXPECT_NEAR(x_one, 0.0, 1e-10) << "Q=1 should return ~0";
  }
}

TEST(SafeGammaMathTest, SafeLogGammaIntegralBasic) {
  // Basic case
  EXPECT_THAT(
      safe_log_gamma_integral(5.0, 1.0, 10.0),
      NearOrNaN(std::log(boost::math::gamma_q(5.0, 1.0) - boost::math::gamma_q(5.0, 10.0)), 1e-12));
}

TEST(SafeGammaMathTest, SafeLogGammaIntegralAtBounds) {
  // When x_min < x_max are far apart
  EXPECT_THAT(
      safe_log_gamma_integral(271.4, 148.0, 6601.0),
      NearOrNaN(std::log(boost::math::gamma_q(271.4, 148.0) - boost::math::gamma_q(271.4, 6601.0)), 1e-12));
  // When both x values are large, the difference is ~0, so log(0) = -inf
  auto result = safe_log_gamma_integral(5.0, 1000.0, 2000.0);
  EXPECT_TRUE(std::isinf(result) && result < 0);
}

TEST(SafeGammaMathTest, SafeSampleTruncatedGammaBasic) {
  std::mt19937 bitgen{42};
  auto alpha = 5.0;
  auto beta = 2.0;
  auto lo = 1.0;
  auto hi = 10.0;

  for (int i = 0; i < 100; ++i) {
    auto x = safe_sample_truncated_gamma(alpha, beta, lo, hi, bitgen);
    EXPECT_GE(x, lo);
    EXPECT_LE(x, hi);
  }
}

TEST(SafeGammaMathTest, SafeSampleTruncatedGammaLargeAlpha) {
  std::mt19937 bitgen{42};
  auto alpha = 271.4;
  auto beta = 1.0;
  auto lo = 100.0;
  auto hi = 500.0;

  for (int i = 0; i < 100; ++i) {
    auto x = safe_sample_truncated_gamma(alpha, beta, lo, hi, bitgen);
    EXPECT_GE(x, lo);
    EXPECT_LE(x, hi);
  }
}

TEST(SafeGammaMathTest, SafeSampleTruncatedGammaVariousParams) {
  std::mt19937 bitgen{42};
  for (double alpha : {1.0, 5.0, 10.0, 50.0, 100.0}) {
    for (double beta : {0.5, 1.0, 2.0, 5.0}) {
      auto lo = alpha / beta;
      auto hi = lo + 10.0;
      for (int i = 0; i < 10; ++i) {
        auto x = safe_sample_truncated_gamma(alpha, beta, lo, hi, bitgen);
        EXPECT_GE(x, lo);
        EXPECT_LE(x, hi);
      }
    }
  }
}

TEST(SafeGammaMathTest, SafeSampleTruncatedGammaWithInfinity) {
  std::mt19937 bitgen{42};
  auto alpha = 5.0;
  auto beta = 1.0;
  auto lo = 1.0;
  auto hi = std::numeric_limits<double>::infinity();

  for (int i = 0; i < 100; ++i) {
    auto x = safe_sample_truncated_gamma(alpha, beta, lo, hi, bitgen);
    EXPECT_GE(x, lo);
  }
}

TEST(SafeGammaMathTest, ComprehensiveSweepGammaQ) {
  // Small a
  for (double a : {0.01, 0.1, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0}) {
    for (double x = 0.001; x <= 1000; x *= 1.7) {
      auto result = safe_gamma_q(a, x);
      EXPECT_FALSE(std::isnan(result)) << "a=" << a << " x=" << x;
      EXPECT_FALSE(std::isinf(result)) << "a=" << a << " x=" << x;
    }
  }

  // Medium a
  for (double a : {20.0, 50.0, 100.0, 150.0, 200.0}) {
    for (double x = 0.01; x <= 10000; x *= 1.5) {
      auto result = safe_gamma_q(a, x);
      EXPECT_FALSE(std::isnan(result)) << "a=" << a << " x=" << x;
      EXPECT_FALSE(std::isinf(result)) << "a=" << a << " x=" << x;
    }
  }

  // Large a (the regime Delphy hits)
  for (double a : {250.0, 271.4, 300.0, 400.0, 500.0, 750.0, 1000.0}) {
    for (double ratio : {0.01, 0.1, 0.5, 0.8, 0.9, 0.95, 0.99, 1.0,
                         1.01, 1.05, 1.1, 1.2, 1.5, 2.0, 5.0, 10.0,
                         20.0, 50.0}) {
      auto result = safe_gamma_q(a, a * ratio);
      EXPECT_FALSE(std::isnan(result)) << "a=" << a << " x=" << a * ratio;
      EXPECT_FALSE(std::isinf(result)) << "a=" << a << " x=" << a * ratio;
    }
    for (double x : {1000.0, 5000.0, 10000.0, 50000.0, 100000.0}) {
      auto result = safe_gamma_q(a, x);
      EXPECT_FALSE(std::isnan(result)) << "a=" << a << " x=" << x;
      EXPECT_FALSE(std::isinf(result)) << "a=" << a << " x=" << x;
    }
    for (double x : {0.001, 0.01, 0.1, 1.0, 10.0}) {
      auto result = safe_gamma_q(a, x);
      EXPECT_FALSE(std::isnan(result)) << "a=" << a << " x=" << x;
      EXPECT_FALSE(std::isinf(result)) << "a=" << a << " x=" << x;
    }
  }

  // x = 0
  for (double a : {1.0, 10.0, 100.0, 1000.0}) {
    EXPECT_NEAR(safe_gamma_q(a, 0.0), 1.0, 1e-12);
  }

  // Near-integer and near-half-integer a
  for (double a : {10.0, 10.5, 100.0, 100.5}) {
    for (double x = 0.1; x <= 1000; x *= 2) {
      auto result = safe_gamma_q(a, x);
      EXPECT_FALSE(std::isnan(result)) << "a=" << a << " x=" << x;
      EXPECT_FALSE(std::isinf(result)) << "a=" << a << " x=" << x;
    }
  }
}

TEST(SafeGammaMathTest, ComprehensiveSweepGammaQInv) {
  for (double a : {0.5, 1.0, 2.0, 5.0, 10.0, 50.0, 100.0, 200.0,
                   271.4, 300.0, 500.0, 1000.0}) {
    for (double Q : {1e-300, 1e-200, 1e-100, 1e-50, 1e-20, 1e-10, 1e-5,
                     0.001, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99,
                     0.999, 0.9999, 1.0 - 1e-10, 1.0 - 1e-15}) {
      if (Q <= 0 || Q >= 1) continue;
      auto result = safe_gamma_q_inv(a, Q);
      EXPECT_FALSE(std::isnan(result)) << "a=" << a << " Q=" << Q;
      EXPECT_FALSE(std::isinf(result)) << "a=" << a << " Q=" << Q;
      EXPECT_GT(result, 0.0) << "a=" << a << " Q=" << Q;
    }
  }
}

TEST(SafeGammaMathTest, PromoteDoubleRegression) {
  // Bare boost::math::gamma_q with default promote_double<true> policy.
  // On Emscripten, long double math functions overflow at the double threshold,
  // causing NaN.  If this starts passing, Emscripten may have fixed the bug —
  // consider removing the promote_double<false> workaround.
  auto bare_result = boost::math::gamma_q(271.4, 6601.0);
#ifdef __EMSCRIPTEN__
  EXPECT_TRUE(std::isnan(bare_result)) << "got " << bare_result;
#else
  EXPECT_NEAR(bare_result, 0.0, 1e-12);
#endif
  // safe_gamma_q uses promote_double<false> and works on all platforms
  EXPECT_NEAR(safe_gamma_q(271.4, 6601.0), 0.0, 1e-12);
}

}  // namespace delphy
