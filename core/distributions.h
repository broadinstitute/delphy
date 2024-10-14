#ifndef DELPHY_DISTRIBUTIONS_H_
#define DELPHY_DISTRIBUTIONS_H_

#include <cmath>

#include <absl/log/check.h>
#include <absl/random/distributions.h>

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

// k ~ Pois(lambda), with the restriction that k >= min_k, i.e.:
//
//         / c lambda^k / k!,      k >= min_k;
//  p(k) = |
//         \ 0,                    otherwise.
//
class K_truncated_poisson_distribution {
 private:
  double lambda_;
  int min_k_;
  double normalization_ = 0.0;
  double term_before_min_k_ = 0.0;
  double max_k_ = 0.0;
  
 public:
  K_truncated_poisson_distribution(double lambda, int min_k)
      : lambda_{lambda}, min_k_{min_k} {
    CHECK_GT(lambda, 0.0);
    CHECK_GE(min_k, 0);
    
    if (min_k_ <= lambda_) {
      // Use dumb rejection sampling
    } else {
      // Use inverse transform sampling.  The trickiest part here is if min_k >> lambda, then
      // the straightforward way of calculating the normalization constant runs into numerical issues!

      // To avoid accumulation of roundoff errors, we cap the max result
      max_k_ = std::max(10.0 * min_k_, 10.0 * lambda_);
      
      DCHECK_NE(min_k_, 0);
      auto last_term = 1.0;
      auto expm1_lambda = std::expm1(lambda_);
      normalization_ = expm1_lambda;  // Note we already exclude k == 0 here!
      for (auto k = 1; k < min_k_; ++k) {
        
        // On iteration entry:
        // * last_term == lambda^{k-1} / (k-1)!
        // * normalization_ == sum_{m=k}^infty (\lambda^m / m!)
        
        last_term *= lambda_ / k;
        normalization_ -= last_term;
        
        // On iteration exit:
        // * last_term == lambda^k / (k-1)!
        // * normalization_ == sum_{m={k+1}}^infty (\lambda^m / m!)
      }

      // Here:
      // * last_term = lambda^{min_k-1} / (min_k-1)!
      // * normalization_ = sum_{m=min_k}^infty (\lambda^m / m!)

      term_before_min_k_ = last_term;

      if (normalization_ <= 0.0 || fabs(normalization_) < 1e-10*expm1_lambda) {
        // Don't trust this calculation owing to roundoff!  Do it by summing up to a large value of k
        // In typical Delphy usecases, this path should never run (we usually have min_k = 1 or 2)

        normalization_ = 0.0;
        auto new_last_term = last_term;
        for (auto k = min_k_; k < max_k_; ++k) {
        
          // On iteration entry:
          // * new_last_term == lambda^{k-1} / (k-1)!
          // * normalization == sum_{m=min_k}^{k-1} (\lambda^m / m!)
          
          new_last_term *= lambda_ / k;
          normalization_ += new_last_term;
          
          // On iteration exit:
          // * new_last_term == lambda^k / (k-1)!
          // * normalization == sum_{m=min_k}^{k}^infty (\lambda^m / m!)
        }
      }
      
      CHECK_GT(normalization_, 0.0);
    }
  }

  template <typename URBG>
  auto operator()(URBG& bitgen) -> int {
    if (normalization_ == 0.0) {
      // Use dumb rejection sampling
      while (true) {
        auto k = std::poisson_distribution<>{lambda_}(bitgen);
        if (k >= min_k_) {
          return k;
        }
      }
    } else {
      // Use inverse transform sampling (see above)
      auto u = absl::Uniform(absl::IntervalClosedOpen, bitgen, 0.0, normalization_);
      
      auto cum_prob = 0.0;
      auto k = min_k_;
      auto term_k = term_before_min_k_;  // term_k == lambda^{min_k-1} / (min_k_ - 1)!
      while (k < max_k_) {
        term_k *= lambda_ / k;  // term_k == lambda^k / k!
        cum_prob += term_k;
        if (cum_prob > u) { break; }
        ++k;
      }
      return k;
    }
  }
};

}  // namespace delphy

#endif // DELPHY_DISTRIBUTIONS_H_
