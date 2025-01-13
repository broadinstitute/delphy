#include "pop_model.h"

#include <absl/log/check.h>
#include <absl/strings/str_format.h>

#include <algorithm>
#include <cmath>

namespace delphy {

Const_pop_model::Const_pop_model(double pop) : pop_{pop} {
  if (pop <= 0.0) {
    throw std::invalid_argument(absl::StrFormat(
        "Population size should be positive (not %f)", pop));
  }
}

auto Const_pop_model::pop_at_time(double) const -> double { return pop_; }
auto Const_pop_model::pop_integral(double a, double b) const -> double { return (b-a) * pop_; }
auto Const_pop_model::intensity_integral(double a, double b) const -> double { return (b-a) / pop_; }

Exp_pop_model::Exp_pop_model(double t0, double pop_at_t0, double growth_rate)
    : t0_{t0}, pop_at_t0_{pop_at_t0}, growth_rate_{growth_rate} {
  if (pop_at_t0 <= 0.0) {
    throw std::invalid_argument(absl::StrFormat(
        "Initial population size should be positive (not %f)", pop_at_t0));
  }
}

auto Exp_pop_model::pop_at_time(double t) const -> double {
  return pop_at_t0_ * std::exp((t - t0_) * growth_rate_);
}

auto Exp_pop_model::pop_integral(double a, double b) const -> double {
  // int_a^b dt' n0 exp[g (t' - t0)] = n0/g [exp[g (b - t0)] - exp[g (a - t0)]]
  //                                 = n0/g exp[g (a - t0)] [exp[g (b - a)] - 1]
  if (growth_rate_ == 0.0) {
    return (b-a) * pop_at_t0_;
  } else {
    return pop_at_t0_ / growth_rate_ *
        std::exp(growth_rate_ * (a - t0_)) *
        std::expm1(growth_rate_ * (b - a));
  }
}

auto Exp_pop_model::intensity_integral(double a, double b) const -> double {
  // int_a^b dt' 1/n0 exp[-g (t - t0))] = 1/(g n0) [exp[-g (a - t0)] - exp[-g (b - t0)]]
  //                                    = 1/(g n0) [exp[-g (b - t0)]] [exp[g (b - a)] - 1]
  if (growth_rate_ == 0.0) {
    return (b-a) / pop_at_t0_;
  } else {
    return (1.0 / (growth_rate_ * pop_at_t0_)) *
        std::exp(-growth_rate_ * (b - t0_)) *
        std::expm1(growth_rate_ * (b - a));
  }
}
  }
}

auto render_population_curve(const Pop_model& pop_model, double t_start, double t_end, int num_t_cells) -> Staircase {
  auto result = Staircase{t_start, t_end, num_t_cells};
  
  for (auto cell = 0; cell != result.num_cells(); ++cell) {
    auto t_lo = cell_lbound(result, cell);
    auto t_hi = cell_ubound(result, cell);
    
    result.at_cell(cell) = pop_model.pop_integral(t_lo, t_hi) / result.cell_size();
  }

  return result;
}

}  // namespace delphy
