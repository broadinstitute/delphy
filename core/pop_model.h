#ifndef DELPHY_POP_MODEL_H_
#define DELPHY_POP_MODEL_H_

#include <ostream>

#include <absl/strings/str_format.h>

#include "staircase.h"

namespace delphy {

class Pop_model {
 public:
  virtual ~Pop_model() = default;

  // N(t) * rho.  Time is measured in calendar units and moves forwards in time from an arbitrary epoch
  // The constant rho is the generation time.
  // Should be positive, but can be zero in special cases (e.g., N(t) = {1 - t for t <= 1, 0 otherwise}).
  // Should handle any value of t, not just nonnegative ones
  virtual auto pop_at_time(double t) const -> double = 0;

  // \int_a^b dt' N(t') * rho
  // We use integrals with two bounds to avoid losing precision when subtracting the results of
  // single-bounded integrals (\int_c^b ... for some fixed c)
  virtual auto pop_integral(double a, double b) const -> double = 0;

  // \int_a^b dt' 1/N(t')
  // We use integrals with two bounds to avoid losing precision when subtracting the results of
  // single-bounded integrals (\int_c^b ... for some fixed c)
  virtual auto intensity_integral(double a, double b) const -> double = 0;

  // Debug printing
  friend auto operator<<(std::ostream& os, const Pop_model& pop_model) -> std::ostream& {
    pop_model.print_to(os);
    return os;
  }

 private:
  // Debug printing details
  virtual auto print_to(std::ostream& os) const -> void = 0;
};

class Const_pop_model : public Pop_model {
 public:
  explicit Const_pop_model(double pop);

  // No setters: change by assigning a new model (resets all params at once and consolidates validation in constructor)
  auto pop() const -> double { return pop_; }

  auto pop_at_time(double t) const -> double override;
  auto pop_integral(double a, double b) const -> double override;
  auto intensity_integral(double a, double b) const -> double override;

 private:
  double pop_;
  
  auto print_to(std::ostream& os) const -> void override { os << absl::StreamFormat("Const_pop_model{pop=%g}", pop()); }
};

class Exp_pop_model : public Pop_model {
 public:
  Exp_pop_model(double t0, double pop_at_t0, double growth_rate);

  // No setters: change by assigning a new model (resets all params at once and consolidates validation in constructor)
  auto t0() const -> double { return t0_; }
  auto pop_at_t0() const -> double { return pop_at_t0_; }
  auto growth_rate() const -> double { return growth_rate_; }

  auto pop_at_time(double t) const -> double override;
  auto pop_integral(double a, double b) const -> double override;
  auto intensity_integral(double a, double b) const -> double override;

 private:
  double t0_;
  double pop_at_t0_;
  double growth_rate_;
  
  auto print_to(std::ostream& os) const -> void override {
    os << absl::StreamFormat("Exp_pop_model{t0=%g, n0=%g, g=%g}",
                             t0(), pop_at_t0(), growth_rate());
  }
};

// A Skygrid model has a population defined piecewise over a finite
// number of predefined intervals:
//
//
//                      |         |       |           |       |
//                      |         |       |           |       |
//         Int.0        |  Int.1  | Int.2 |    ...    | Int.M |   Int.{M+1}
//                      |         |       |           |       |
//                      |         |       |           |       |
//    -INFINITY ...-----+---------+-------+----....---+-------+---... +INFINITY  --> t
//                     x_0       x_1     x_2       x_{M-1}   x_M
//
// There time range (-INFINITY, +INFINITY) is partitioned into M+1 intervals,
// the first and last of which are open-ended.  The boundaries bewteen consecutive intervals
// are called "knots", and are specified at x_0 < x_1 < ... < x_M.  The interval
// for a given time t is defined as:
//
//               /   0,                     t <= x_0;
// Interval(t) = |   k,           x_{k-1} < t <= x_k;    (1 <= k <= M)
//               \   M+1,             x_M < t.
//
// The log-population-size log(N(t)) is specified at each of the knots:
//
//     log N(t_k) =: gamma_k     0 <= k <= M.
//
// We support two variants for interpolating the population size at all other times.
//
// In the traditional Skygrid model (Gill et al 2012, BEAST's "gmrfSkyGridLikelihood",
// k_staircase below), the population is constant over the length of a single interval,
// but can vary from interval to interval:
//
//  Staircase
//  ---------
//                 / gamma_0,            t <= x_0;       // Interval 0
//      log N(t) = | gamma_k,  x_{k-1} < t <= x_k;       // Interval 1 <= k <= M
//                 \ gamma_M,      x_M < t.              // Interval M+1
//
// A somewhat more natural choice (to me, at least; k_log_linear below) is
// to have the population to grow exponentially within a single interval,
// be continuous across intervals, but allowing different growth rates
// (possibly negative) in different intervals.  We parametrize this curve as follows:
//
//  Log-linear
//  ----------
//                 / gamma_0,                                  t <= x_0;  // Interval 0
//      log N(t) = | (1-c) gamma_{k-1} + c gamma_k,  x_{k-1} < t <= x_k;  // Interval 1 <= k <= M
//                 \ gamma_M,                            x_M < t.         // Interval M+1
//
//         [N.B.: for intervals 1 <= k <= M, we have
//                  t = (1-c) x_{k-1} + c x_k
//               => c = (t - x_{k-1}) / (x_k - x_{k-1}) ]
//
class Skygrid_pop_model : public Pop_model {
 public:
  enum class Type {

    // Staircase = the traditional Skygrid model from Gill et al 2012 and BEAST.
    // We allow arbitrary knot times satisfying x_0 < x_1 < ... < x_M.
    // The default choice of knot times is x_k = T - (M-k)/M * K,
    // where T is the time of the latest tip (its upper bound if the date is uncertain),
    // and K is the cutoff duration.  The main notational discrepancy is that,
    // unlike in Gill et al 2012 or in BEAST, Delphy's internal time axis
    // increases towards the future and has a fixed epoch at 2020-01-01:
    //
    //  Here                 Gill et al 2012
    //  ----                 ---------------
    //  t                    T-t
    //  M                    M
    //  T                    0
    //  K                    K
    //  x_k                  x_{M-k}
    //  gamma_k              gamma_{M+1-k}
    //  exp(gamma_k)         theta_{M+1-k}
    
    k_staircase = 1,

    // Log-linear = A continuous population curve s.t. log(N(t)) is continuous
    // and linearly interpolates between the values at the knots.
    // 
    // From looking at published SkyGrid population curve estimates, this is what I
    // _thought_ SkyGrid was doing, before reading the technical details. We keep the
    // original model around for benchmarking against BEAST, but this one seems
    // more natural to me.

    k_log_linear = 2
  };
  
  Skygrid_pop_model(
      std::vector<double> x,      // x[k] = time of knot k
      std::vector<double> gamma,  // gamma[k] = log(N(x[k]))
      Type type);

  // No setters: change by assigning a new model (resets all params at once and consolidates validation in constructor)
  auto x() const -> const std::vector<double>& { return x_; }
  auto x_lo() const -> double { return x_.front(); }
  auto x_hi() const -> double { return x_.back(); }
  auto gamma() const -> const std::vector<double>& { return gamma_; }
  auto type() const -> Type { return type_; }

  // Derived quantities in convenient notation (see notation above)
  auto M() const -> int { return std::ssize(x_) - 1; }
  auto T() const -> double { return x_hi(); }
  auto K() const -> double { return x_hi() - x_lo(); }
  auto x(int k) const -> double { return x_.at(k); }
  auto gamma(int k) const -> double { return gamma_.at(k); }
  auto gamma_below_x_lo() const -> double { return gamma_.front(); }
  auto gamma_above_x_hi() const -> double { return gamma_.back(); }

  // The key functions
  auto pop_at_time(double t) const -> double override;
  auto pop_integral(double a, double b) const -> double override;
  auto intensity_integral(double a, double b) const -> double override;

  // Log-space value and Gradient calculation (for supporting HMC)
  auto log_N(double t) const -> double;  // == std::log(pop_at_time(t))
  auto d_log_N_d_gamma(double t, int k) const -> double;
  auto support_of_d_log_N_d_gamma(int k) const -> std::pair<double, double>;  // range of (inclusive) over which d_log_N_d_gamma(t,k) might be non-zero
  auto d_log_int_N_d_gamma(double a, double b, int k) const -> double;

  // Interval(t) in the top-level comment for Skygrid_pop_model
  auto interval_containing_t(double t) const -> int;

 private:
  std::vector<double> x_;
  std::vector<double> gamma_;
  std::vector<double> minus_gamma_;
  Type type_;
  
  auto print_to(std::ostream& os) const -> void override;

  // pop_integral and intensity_integral are identical except for a flip in the sign of gamma_k
  auto log_int_N_core(double a, double b, const std::vector<double>& gamma_eff) const -> double;
  auto d_log_int_N_d_gamma_core(double a, double b, int k, const std::vector<double>& gamma_eff) const -> double;
};

// Utilities
// ---------

// Each cell spanning [t, t+Delta) records 1/Delta int_t^{t+Delta} dt' N(t')
auto render_population_curve(const Pop_model& pop_model, double t_start, double t_end, int num_t_cells) -> Staircase;

}  // namespace delphy

#endif // DELPHY_POP_MODEL_H_
