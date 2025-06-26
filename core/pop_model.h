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

// Utilities
// ---------

// Each cell spanning [t, t+Delta) records 1/Delta int_t^{t+Delta} dt' N(t')
auto render_population_curve(const Pop_model& pop_model, double t_start, double t_end, int num_t_cells) -> Staircase;

}  // namespace delphy

#endif // DELPHY_POP_MODEL_H_
