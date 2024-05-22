#ifndef DELPHY_DISTRIBUTIONS_H_
#define DELPHY_DISTRIBUTIONS_H_

#include <cmath>

#include <absl/log/check.h>

namespace delphy {

// x is distributed exponentially from a to b with the following pdf:
//
//         / c exp(lambda x),      a < x < b;
//  p(x) = |
//         \ 0,                    otherwise.
//
//  The CDF of this is
//
//         /  0,                                             x <= a;
//  P(x) = |  (c / lambda) [exp(lambda x) - exp(lambda a)],  a < x < b;
//         \  1,                                             b <= x.
//
// For normalization, c = lambda / [exp(lambda b) - exp(lambda a)].
//
// Hence, to sample from p(x), we draw a uniform value u in (0,1) and set
//
//  u = P(x) => x = a + (1/lambda) ln[1 + u * (e^{lambda (b - a)} - 1)].
//
// When |lambda (b-a)| >> 1, either because the ratio really is big or because
// one of the bounds is +/- infinity, we can use the following approximations instead:
//
//  x = b + (1/lambda) ln(u)   (if lambda > 0)
//
// or
//
//  x = a + (1/lambda) ln(u)   (if lambda < 0).
//
class Bounded_exponential_distribution {
 public:
  Bounded_exponential_distribution(double lambda, double a, double b)
      : lambda_{lambda}, a_{a}, b_{b} {
    CHECK(a_ <= b_) << a_ << " <=? " << b_;
    CHECK(not (std::isinf(a_) && std::isinf(b_)));
    CHECK(not (lambda > 0.0 && std::isinf(b_)));
    CHECK(not (lambda < 0.0 && std::isinf(a_)));
  }

  template <typename URBG>
  auto operator()(URBG& bitgen) -> double {
    auto u = absl::Uniform(absl::IntervalOpenOpen, bitgen, 0.0, 1.0);
    auto lambda_times_range = lambda_ * (b_ - a_);
    if (lambda_ == 0.0) {
      return bound(a_ + u * (b_ - a_));
    } else if (lambda_ > 0 && lambda_times_range > 100) {
      return bound(b_ + std::log(u) / lambda_);
    } else if (lambda_ < 0 && lambda_times_range < -100) {
      return bound(a_ + std::log(u) / lambda_);
    } else {
      return bound(a_ + std::log1p(u * (std::exp(lambda_times_range) - 1)) / lambda_);  // log1p(x) = log(1 + x)
    }
  }
  
 private:
  double lambda_;
  double a_;
  double b_;

  auto bound(double x) -> double { return std::ranges::clamp(x, a_, b_); }
};

}  // namespace delphy

#endif // DELPHY_DISTRIBUTIONS_H_
