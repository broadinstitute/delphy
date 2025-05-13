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

Skygrid_pop_model::Skygrid_pop_model(
    std::vector<double> x,
    std::vector<double> gamma,
    Type type)
    : x_{std::move(x)}, gamma_{std::move(gamma)}, type_{type} {
  
  if (std::ssize(x_) < 2) {
    throw std::invalid_argument(absl::StrFormat("Skygrid_pop_model needs at least two knots"));
  }
  if (std::ssize(x_) != std::ssize(gamma_)) {
    throw std::invalid_argument(absl::StrFormat(
        "Skygrid_pop_model needs knot times array (`x`) and knot values array (`gamma`) to "
        "have equal size, but there are %d x's and %d gamma's",
        std::ssize(x_), std::ssize(gamma_)));
  }

  // Strictly increasing knot times
  auto M = std::ssize(x_) - 1;
  for (auto i = 0; i != M; ++i) {
    if (not (x_[i] < x_[i+1])) {
      throw std::invalid_argument(absl::StrFormat(
          "Skygrid_pop_model needs strictly increasing knot times, but "
          "x[%d] = %g >= x[%d] = %g",
          i, x_[i], i+1, x_[i+1]));
    }
  }

  // Pre-calculate -gamma
  minus_gamma_.assign(std::ssize(gamma_), 0.0);
  for (auto k = 0; k != std::ssize(gamma_); ++k) {
    minus_gamma_[k] = -gamma_[k];
  }
}

auto Skygrid_pop_model::log_N(double t) const -> double {
  auto k = interval_containing_t(t);
  auto M = this->M();
  
  switch (type_) {
    
    case Type::k_staircase:
      if      (k == 0) { return gamma(0); }
      else if (k <= M) { return gamma(k); }
      else             { return gamma(M); }
      
    case Type::k_log_linear:
      if      (k == 0) { return gamma(0); }
      else if (k <= M) { auto c = (t - x(k-1)) / (x(k) - x(k-1)); return (1-c)*gamma(k-1) + c*gamma(k); }
      else             { return gamma(M); }
      
    default:
      CHECK(false) << "unrecognized type " << static_cast<int>(type_);
  }
}

auto Skygrid_pop_model::pop_at_time(double t) const -> double {
  return std::exp(log_N(t));
}

auto Skygrid_pop_model::d_log_N_d_gamma(double t, int k) const -> double {
  auto kk = interval_containing_t(t);
  auto M = this->M();
  
  switch (type_) {
    
    case Type::k_staircase:
      if      (kk == 0) { return (k == 0)  ? 1 : 0; }
      else if (kk <= M) { return (k == kk) ? 1 : 0; }
      else              { return (k == M)  ? 1 : 0; }
      
    case Type::k_log_linear:
      if      (kk == 0) { return (k == 0) ? 1 : 0; }
      else if (kk <= M) { auto c = (t - x(kk-1)) / (x(kk) - x(kk-1)); return (k == kk-1) ? (1-c) : (k == kk) ? c : 0; }
      else              { return (k == M) ? 1 : 0; }
      
    default:
      CHECK(false) << "unrecognized type " << static_cast<int>(type_);
  }
}

auto Skygrid_pop_model::support_of_d_log_N_d_gamma(int k) const -> std::pair<double, double> {
  switch (type_) {
    
    case Type::k_staircase: {
      auto lo = (k ==   0) ? -std::numeric_limits<double>::infinity() : x_[k-1];
      auto hi = (k == M()) ? +std::numeric_limits<double>::infinity() : x_[k];
      return {lo, hi};
    }
      
    case Type::k_log_linear: {
      auto lo = (k ==   0) ? -std::numeric_limits<double>::infinity() : x_[k-1];
      auto hi = (k == M()) ? +std::numeric_limits<double>::infinity() : x_[k+1];
      return {lo, hi};
    }
      
    default:
      CHECK(false) << "unrecognized type " << static_cast<int>(type_);
  }
}

auto Skygrid_pop_model::log_integral_core(double a, double b, const std::vector<double>& gamma_eff) const -> double {
  CHECK_LE(a, b);

  auto M = this->M();
  CHECK_EQ(std::ssize(gamma_eff), M+1);

  auto ka = interval_containing_t(a);
  auto kb = interval_containing_t(b);

  // Bias ln(integral) by maximum gamma that enters the integral
  // (bias is compensated for at the end of the calculation)
  auto kka = std::max(ka-1, 0);
  auto kkb = std::min(kb,   M);
  CHECK_LE(kka, kkb);
  auto bias = -std::numeric_limits<double>::infinity();
  for (auto k = kka; k <= kkb; ++k) {  // k == kkb is _included_
    bias = std::max(bias, gamma_eff.at(k));
  }

  auto result = 0.0; // == integral_core(a,b,gamma_eff) / std::exp(bias)
  
  // Decompose integral into one piece per overlapping interval
  for (auto k = ka; k <= kb; ++k) { // k == kb is _included_

    // Intersect [a,b] with limits of interval k
    auto lo = std::max(a, (k > 0 ) ? x(k-1) : -std::numeric_limits<double>::infinity());
    auto hi = std::min(b, (k <= M) ? x(k)   : +std::numeric_limits<double>::infinity());
    CHECK_LE(lo, hi);
    
    // Integrate N(t) dt over this intersection:
    //
    //    result += int_lo^hi dt * N(t)/exp(bias)

    if (k == 0) {
      result += std::exp(-bias + gamma_eff[0]) * (hi-lo);
      
    } else if (k == M+1) {
      result += std::exp(-bias + gamma_eff[M]) * (hi-lo);

    } else {
      
      switch (type_) {
        
        case Type::k_staircase:
          result += std::exp(-bias + gamma_eff[k]) * (hi-lo);
          break;
          
        case Type::k_log_linear:
          //
          // int_lo^hi dt exp[-bias + gamma(t)]   (where gamma(t) linearly interpolates
          //                                       from gamma_{k-1} at x_{k-1}
          //                                       to gamma_k at x_k)
          //  = int_lo^hi dt exp[-bias + gamma_lo + (t - lo) * dg/(hi-lo)]   { dg := gamma_hi-gamma_lo }
          //  = exp[-bias + gamma_lo] * { exp(dg) - 1 } * (hi-lo)/dg
          //  = exp[-bias + gamma_lo] * (hi-lo) * { expm1(dg) / dg }

          if (gamma_eff[k] == gamma_eff[k-1]) {
            result += std::exp(-bias + gamma_eff[k]) * (hi - lo);
          } else {
            auto gamma_lo = gamma_eff[k-1] + (lo - x(k-1)) * (gamma_eff[k] - gamma_eff[k-1]) / (x(k) - x(k-1));
            auto gamma_hi = gamma_eff[k-1] + (hi - x(k-1)) * (gamma_eff[k] - gamma_eff[k-1]) / (x(k) - x(k-1));
            auto dg = gamma_hi - gamma_lo;
            result += std::exp(-bias + gamma_lo) * (hi - lo) * (std::expm1(dg) / dg);
          }
          break;
          
        default:
          CHECK(false) << "unrecognized type " << static_cast<int>(type_);
      }
    }
  }
  
  return std::log(result) + bias;
}

auto Skygrid_pop_model::pop_integral(double a, double b) const -> double {
  return std::exp(log_integral_core(a, b, gamma_));
}

auto Skygrid_pop_model::intensity_integral(double a, double b) const -> double {
  return std::exp(log_integral_core(a, b, minus_gamma_));
}

auto Skygrid_pop_model::print_to(std::ostream& os) const -> void {
  os << absl::StreamFormat("Skygrid_pop_model{");
  os << absl::StreamFormat("type=%s, ",
                           (type_ == Type::k_staircase ? "k_staircase" :
                            type_ == Type::k_log_linear ? "k_log_linear" :
                            "???"));
  os << absl::StreamFormat("ln N(t)=[");
  for (auto k = 0; k != M()+1; ++k) {
    if (k != 0) { os << "; "; }
    os << absl::StreamFormat("at t=%g, %g", x(k), gamma(k));
  }
  os << absl::StreamFormat("]}");
}

auto Skygrid_pop_model::interval_containing_t(double t) const -> int {
  auto it = std::ranges::lower_bound(x_, t);  // upper_bound means *(it-1) < t <= *it
  if (it == x_.begin()) {
    return 0;
  } else if (it == x_.end()) {
    return M()+1;
  } else {
    return static_cast<int>(std::ranges::distance(x_.begin(), it));
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
