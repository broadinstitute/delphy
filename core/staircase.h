#ifndef DELPHY_STAIRCASE_H_
#define DELPHY_STAIRCASE_H_

#include <iostream>
#include <vector>
#include <span>
#include <cmath>
#include <stdexcept>

#include "absl/strings/str_format.h"

namespace delphy {

// Specifies a function f(x) over the interval [x_start, x_end).
// The domain of f(x) is divided into num_cells cells of
// width cell_size = (x_end - x_start) / num_cells.  Each cell i spans the
// interval [i_start, i_end), with i_start = x_start + i * cell_size and i_end = i_start + cell_size.
// The value f[i] determines the value of the function f(x) in that interval.
class Staircase {
private:
  double x_start_;
  double cell_size_;
  std::vector<double> f_;
  
public:
  using value_type = std::vector<double>::value_type;
  using reference = std::vector<double>::reference;
  using const_reference = std::vector<double>::const_reference;
  using iterator = std::vector<double>::iterator;
  using const_iterator = std::vector<double>::const_iterator;
  using difference_type = std::vector<double>::difference_type;
  using size_type = std::vector<double>::size_type;
  
  Staircase(double x_start, double x_end, int num_cells)
    : x_start_{x_start}, cell_size_{std::numeric_limits<double>::signaling_NaN()}, f_(num_cells, 0.0) {
    if (not (x_start < x_end)) {
      throw std::invalid_argument(absl::StrFormat(
          "Invalid domain: need x_start < x_end, but x_start=%f and x_end=%f", x_start, x_end));
    }
    if (num_cells <= 0) {
      throw std::invalid_argument(absl::StrFormat("Number of cells should be positive, not %d", num_cells));
    }
    
    cell_size_ = (x_end - x_start) / num_cells;
  }
  
  auto x_start() const -> double { return x_start_; }
  auto cell_size() const -> double { return cell_size_; }
  auto num_cells() const -> int { return std::ssize(f_); }
  auto x_end() const -> double { return x_start_ + num_cells() * cell_size_; }

  auto at_cell(int cell) -> double& {
    assert(cell >= 0);
    assert(cell < num_cells());
    return f_[cell];
  }
  auto at_cell(int cell) const -> double { return const_cast<Staircase&>(*this).at_cell(cell); }
  
  auto cell_for(double x) const -> int { return static_cast<int>(std::floor((x - x_start()) / cell_size())); }
  auto at(double x) -> double& { return at_cell(cell_for(x)); }
  auto at(double x) const -> double { return at_cell(cell_for(x)); }

  auto begin() -> iterator { return f_.begin(); }
  auto end() -> iterator { return f_.end(); }
  auto begin() const -> const_iterator { return f_.begin(); }
  auto end() const -> const_iterator { return f_.end(); }
};

// Delicate boundary calculations that are careful about open-closed intervals and roundoff
inline auto cell_for_lbound(const Staircase& staircase, double x) -> int {
  assert(staircase.x_start() <= x);
  assert(x < staircase.x_end());
  return static_cast<int>(std::floor((x - staircase.x_start()) / staircase.cell_size()));
}
inline auto cell_for_ubound(const Staircase& staircase, double x) -> int {
  assert(staircase.x_start() < x);
  assert(x <= staircase.x_end());
  auto cell_offset = std::floor((staircase.x_end() - x) / staircase.cell_size());
  return staircase.num_cells()-1 - static_cast<int>(cell_offset);
}

inline auto cell_lbound(const Staircase& staircase, int cell) -> double {
  return staircase.x_start() + cell * staircase.cell_size();
}
inline auto cell_ubound(const Staircase& staircase, int cell) -> double {
  return cell_lbound(staircase, cell) + staircase.cell_size();
}

// A boxcar is a function with value `height` for `left < x < right` and with value 0 elsewhere
auto add_boxcar(Staircase& staircase, double left, double right, double height) -> void;

// A trapezoid is a function that interpolates linearly from `left_height` at `x = left`
// to `right_height` at `x = right`, and is 0 outside the range `left <= x <= right`
auto add_trapezoid(Staircase& staircase, double left, double right, double left_height, double right_height) -> void;


// A family of staircases with the same x_start, x_end and num_cells
class Staircase_family {
private:
  std::vector<Staircase> members_;
  
public:
  using iterator = std::vector<Staircase>::iterator;
  using const_iterator = std::vector<Staircase>::const_iterator;

  Staircase_family(int num_members, double x_start, double x_end, int num_cells) {
    if (num_members <= 0) {
      throw std::invalid_argument(absl::StrFormat("Number of members should be positive, not %d", num_members));
    }
    members_.reserve(num_members);
    for (int i = 0; i != num_members; ++i) {
      members_.emplace_back(x_start, x_end, num_cells);
    }
  }

  auto num_members() const -> int { return static_cast<int>(std::ssize(members_)); }
  auto x_start() const -> double { return members_[0].x_start(); }
  auto cell_size() const -> double { return members_[0].cell_size(); }
  auto num_cells() const -> int { return members_[0].num_cells(); }
  auto x_end() const -> double { return members_[0].x_end(); }

  auto begin() -> iterator { return members_.begin(); }
  auto end() -> iterator { return members_.end(); }
  auto begin() const -> const_iterator { return members_.begin(); }
  auto end() const -> const_iterator { return members_.end(); }

  auto operator[](int member) -> Staircase& {
    assert(member >= 0);
    assert(member < std::ssize(members_));
    return members_[member];
  }
  auto operator[](int member) const -> const Staircase& { return const_cast<Staircase_family&>(*this)[member]; }
};

// Delicate boundary calculations that are careful about open-closed intervals and roundoff
inline auto cell_for_lbound(const Staircase_family& family, double x) -> int { return cell_for_lbound(family[0], x); }
inline auto cell_for_ubound(const Staircase_family& family, double x) -> int { return cell_for_ubound(family[0], x); }
inline auto cell_lbound(const Staircase_family& family, int cell) -> double { return cell_lbound(family[0], cell); }
inline auto cell_ubound(const Staircase_family& family, int cell) -> double { return cell_ubound(family[0], cell); }

}  // namespace delphy

#endif // DELPHY_STAIRCASE_H_
