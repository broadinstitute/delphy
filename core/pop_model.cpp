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
  
  if (x_.empty()) {
    throw std::invalid_argument(absl::StrFormat("Skygrid_pop_model needs at least one knot"));
  }
  if (std::ssize(x_) != std::ssize(gamma_)) {
    throw std::invalid_argument(absl::StrFormat(
        "Skygrid_pop_model needs knot times array (`x`) and knot values array (`gamma`) to "
        "have equal size, but there are %d x's and %d gamma's",
        std::ssize(x_), std::ssize(gamma_)));
  }

  // Strictly increasing knot times
  auto M = std::ssize(x_) - 1;
  for (auto i = 0; i < M; ++i) {
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
  auto [k, c] = k_and_c(t);
  switch (type_) {
    
    case Type::k_staircase:
      if (t <= x_[0]) {
        // here, [k,c] == [1,0.0] by design, but we want "[0,1.0]"
        // (but there's no x[-1] here to make that statement make sense)
        return gamma(0);
      } else {
        return gamma(k);
      }
      
    case Type::k_log_linear:
      return (1-c)*gamma(k-1) + c*gamma(k);
      
    default:
      CHECK(false) << "unrecognized type " << static_cast<int>(type_);
  }
}

auto Skygrid_pop_model::pop_at_time(double t) const -> double {
  return std::exp(log_N(t));
}

auto Skygrid_pop_model::d_log_N_d_gamma(double t, int k) const -> double {
  auto [kk, cc] = k_and_c(t);
  switch (type_) {
    
    case Type::k_staircase:
      if (t <= x_[0]) {
        // here, [kk,cc] == [1,0.0] by design, but we want "[0,1.0]"
        // (but there's no x[-1] here to make that statement make sense)
        return (k == 0) ? 1 : 0;
      } else {
        return (k == kk) ? 1 : 0;
      }
      
    case Type::k_log_linear:
      if      (k == kk-1) { return (1-cc); }
      else if (k == kk)   { return cc;     }
      else                { return 0;      }
      
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

auto Skygrid_pop_model::integral_core(double a, double b, const std::vector<double>& gamma_eff) const -> double {
  CHECK_LE(a, b);

  auto M = this->M();
  CHECK_EQ(std::ssize(gamma_eff), M+1);
  
  // Case 1: Integral fully before lower limit
  if (b <= x_lo()) {
    switch (type_) {
      
      case Type::k_staircase:
      case Type::k_log_linear:
        return std::exp(gamma_eff[0]) * (b-a);
        
      default:
        CHECK(false) << "unrecognized type " << static_cast<int>(type_);
    }
  }

  // Case 2: Integral fully above upper limit
  if (a >= x_hi()) {
    switch (type_) {
      
      case Type::k_staircase:
      case Type::k_log_linear:
        return std::exp(gamma_eff[M]) * (b-a);
        
      default:
        CHECK(false) << "unrecognized type " << static_cast<int>(type_);
    }
  }

  // Case 3: integration range overlaps (t_lo(), t_hi())
  auto result = 0.0;

  // Part between `a` and lower limit
  if (a <= x_lo()) {
    switch (type_) {
      
      case Type::k_staircase:
      case Type::k_log_linear:
        result += std::exp(gamma_eff[0]) * (x_lo() - a);
        break;
        
      default:
        CHECK(false) << "unrecognized type " << static_cast<int>(type_);
    }

    a = x_lo();  // Adjust range!
  }

  // Part between upper limit and `b`
  if (b >= x_hi()) {
    switch (type_) {
      
      case Type::k_staircase:
      case Type::k_log_linear:
        result += std::exp(gamma_eff[M]) * (b - x_hi());
        break;
        
      default:
        CHECK(false) << "unrecognized type " << static_cast<int>(type_);
    }

    b = x_hi();  // Adjust range!
  }

  // Remaining processing is strictly within knot limits
  CHECK_LE(x_lo(), a);
  CHECK_LE(a, b);
  CHECK_LE(b, x_hi());

  // First, a quick helper for integrating the population curve within a single cell
  auto single_cell_integral = [&gamma_eff, this](int k, double ca, double cb) -> double {
    switch (type_) {
      
      case Type::k_staircase:
        return
            (x(k) - x(k-1))
            * std::exp(gamma_eff[k])
            * (cb - ca);
        break;
        
      case Type::k_log_linear: {
        // t = (1-c) x_{k-1} + c x_k  => c = (t - x_{k-1}) / (x_k - x_{k-1})  => dc = dt / (x_k - x_{k-1})
        // int_a^b dt' exp((1-c) gamma_{k-1} + c gamma_k)
        //  = int_ca^cb (x_k - x_{k-1}) dc exp(gamma_{k-1}) * exp(c delta_gamma])  [delta_gamma := gamma_k - gamma_{k-1})
        //  = (x_k - x_{k-1}) exp(gamma_{k-1}) * (exp(delta_gamma cb) - exp(delta_gamma ca)) / delta_gamma
        //  = (x_k - x_{k-1}) exp(gamma_{k-1}) * exp(delta_gamma ca) * [exp(delta_gamma (cb-ca)) - 1] / delta_gamma
        if (auto delta_gamma_eff = gamma_eff[k] - gamma_eff[k-1]; delta_gamma_eff == 0.0) {
          return
              (x(k) - x(k-1))
              * std::exp(gamma_eff[k-1])  // == std::exp(gamma_eff[k])
              * (cb - ca);
        } else {
          return
              (x(k) - x(k-1))
              * std::exp(gamma_eff[k-1])
              * (cb - ca)
              * std::exp(delta_gamma_eff * ca)
              * (std::expm1(delta_gamma_eff * (cb-ca)) / (delta_gamma_eff * (cb-ca)));
        }
        break;
        
        default:
          CHECK(false) << "unrecognized type " << static_cast<int>(type_);
      }
    }
  };
  
  auto [ka, ca] = k_and_c(a);
  auto [kb, cb] = k_and_c(b);

  if (ka == kb) {
    // Special case: integration domain is entirely inside a single cell
    result += single_cell_integral(kb, ca, cb);
    
  } else {

    // General case: the domain...
    // 1. ...starts in the middle of cell x_{ka-1} < x <= x_{ka};
    // 2. ...then covers 0+ complete cells;
    // 3. ...then ends in the middle of cell x_{kb-1} < x <= x_{kb};

    // Partial cell ka
    CHECK_LE(x(ka-1), a);
    CHECK_LE(a, x(ka));
    result += single_cell_integral(ka, (a - x(ka-1)) / (x(ka) - x(ka-1)), 1.0);

    // Intermediate cells
    for (int k = ka+1; k != kb; ++k) {
      result += single_cell_integral(k, 0.0, 1.0);
    }

    // Partial cell kb
    CHECK_LE(x(kb-1), b);
    CHECK_LE(b, x(kb));
    result += single_cell_integral(kb, 0.0, (b - x(kb-1)) / (x(kb) - x(kb-1)));
  }
  
  return result;
}

auto Skygrid_pop_model::pop_integral(double a, double b) const -> double {
  return integral_core(a, b, gamma_);
}

auto Skygrid_pop_model::intensity_integral(double a, double b) const -> double {
  return integral_core(a, b, minus_gamma_);
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

auto Skygrid_pop_model::k_and_c(double t) const -> std::pair<int, double> {
  auto it = std::ranges::lower_bound(x_, t);  // upper_bound means *(it-1) < t <= *it
  if (it == x_.begin()) {
    return {1, 0.0};
  } else if (it == x_.end()) {
    return {M(), 1.0};
  } else {
    auto prev = std::ranges::prev(it);
    auto t_km1 = *prev;
    auto t_k = *it;
    auto c = (t - t_km1) / (t_k - t_km1);
    return {std::ranges::distance(x_.begin(), it), c};
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
