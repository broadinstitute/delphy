#include "pop_model.h"

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
auto Const_pop_model::cum_pop_at_time(double t) const -> double { return t * pop_; }
auto Const_pop_model::intensity_at_time(double t) const -> double { return t / pop_; }

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

auto Exp_pop_model::cum_pop_at_time(double t) const -> double {
  // int_t0^t dt' n0 exp[g (t' - t0)] = n0/g [exp[g (t - t0)] - 1]
  if (growth_rate_ == 0.0) {
    return t * pop_at_t0_;
  } else {
    return (std::exp(growth_rate_ * (t - t0_)) - 1) * pop_at_t0_ / growth_rate_;
  }
}

auto Exp_pop_model::intensity_at_time(double t) const -> double {
  // int_t0^t dt' 1/n0 exp[-g (t - t0))] = 1/(g n0) [1 - exp[-g (t-t0)]]
  if (growth_rate_ == 0.0) {
    return t / pop_at_t0_;
  } else {
    return (1 - std::exp(-growth_rate_ * (t-t0_))) / (growth_rate_ * pop_at_t0_);
  }
}

auto render_population_curve(const Pop_model& pop_model, double t_start, double t_end, int num_t_cells) -> Staircase {
  auto result = Staircase{t_start, t_end, num_t_cells};

  auto t_lo = cell_lbound(result, 0);
  auto pop_at_t_lo = pop_model.cum_pop_at_time(t_lo);
  
  for (auto cell = 0; cell != result.num_cells(); ++cell) {
    auto t_hi = cell_ubound(result, cell);
    auto pop_at_t_hi = pop_model.cum_pop_at_time(t_hi);
    
    result.at_cell(cell) = (pop_at_t_hi - pop_at_t_lo) / result.cell_size();
    
    // next cell
    t_lo = t_hi;
    pop_at_t_lo = pop_at_t_hi;
  }

  return result;
}

}  // namespace delphy
