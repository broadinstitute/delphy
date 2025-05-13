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

class Skygrid_pop_model : public Pop_model {
 public:
  enum class Type {
    
    // The original Skygrid model (Gill et al 2012, BEAST's "gmrfSkyGridLikelihood"),
    // with constant population between M+1 knots at times x_k:
    //
    //        / exp(gamma_0),            t <= x_0;
    // N(t) = | exp(gamma_k),  x_{k-1} < t <= x_k;
    //        \ exp(gamma_M),      x_M < t.
    //
    // Note: we allow arbitrary knot times satisfying t_0 < t_1 < ... < t_M.
    // The default choice of knot times is t_k = T - (M-k)/M * K,
    // where T is the time of the latest tip (its lower bound if the date is uncertain),
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
    //  gamma_k              gamma_{1 + M-k}
    //  exp(gamma_k)         theta_{1 + M-k}
    
    k_staircase = 1,

    // A continuous population curve s.t. log(N(t)) is continuous and linearly
    // interpolates between the values at the knots.
    // 
    // From looking at published SkyGrid population curve estimates, this is what I
    // _thought_ SkyGrid was doing, before reading the technical details. We keep the
    // original model around for benchmarking against BEAST, but this one seems
    // more natural to me.
    //
    //        / exp(gamma_0),                                  t <= x_0;
    //        |
    // N(t) = | exp((1-c) gamma_{k-1} + c gamma_k),  x_{k-1} < t <= x_k;
    //        | [t= (1-c)     x_{k-1} + c     x_k]
    //        |
    //        \ exp(gamma_M),                            x_M < t.

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

 private:
  std::vector<double> x_;
  std::vector<double> gamma_;
  std::vector<double> minus_gamma_;
  Type type_;
  
  auto print_to(std::ostream& os) const -> void override;

  // auto [k,c] = k_and_c(t);
  // Calculates the interval k >= 1 containing time t, such that t = (1-c) x_{k-1} + c x_k
  // If t <= x_0, returns (1, 0.0).  If t >= x_M, returns (M, 1.0)
  auto k_and_c(double t) const -> std::pair<int, double>;

  // pop_integral and intensity_integral are identical except for a flip in the sign of gamma_k
  auto integral_core(double a, double b, const std::vector<double>& gamma_eff) const -> double;
};

// Utilities
// ---------

// Each cell spanning [t, t+Delta) records 1/Delta int_t^{t+Delta} dt' N(t')
auto render_population_curve(const Pop_model& pop_model, double t_start, double t_end, int num_t_cells) -> Staircase;

}  // namespace delphy

#endif // DELPHY_POP_MODEL_H_
