#include "staircase.h"

namespace delphy {

auto add_boxcar(Staircase& staircase, double left, double right, double height) -> void {
  if (not (left <= right)) {
    throw std::invalid_argument(absl::StrFormat(
        "Invalid domain: need left <= right, but left=%f and right=%f", left, right));
  }

  // Focus on overlap between [left,right] and [start,end]
  if (left > staircase.x_end() || right < staircase.x_start()) return;  // Doesn't overlap domain
  left = std::max(left, staircase.x_start());
  right = std::min(right, staircase.x_end());
  if (left == right) return;  // Empty boxcar
  assert(staircase.x_start() <= left && left < right && right <= staircase.x_end());
  
  auto cell_start = cell_for_lbound(staircase, left);
  auto cell_end = cell_for_ubound(staircase, right);

  if (cell_start == cell_end) {
    staircase.at_cell(cell_start) += height * (right - left) / staircase.cell_size();
  } else {
    // Handle partial overlaps with first and last cell
    auto iter = staircase.begin() + cell_start;
    *iter += height * (cell_ubound(staircase, cell_start) - left) / staircase.cell_size();
    ++iter;

    auto end = staircase.begin() + cell_end;
    *end += height * (right - cell_lbound(staircase, cell_end)) / staircase.cell_size();

    // The middle is easy
    while (iter != end) {
      *iter += height;
      ++iter;
    }
  }
}

auto add_trapezoid(Staircase& staircase, double left, double right, double left_height, double right_height) -> void {
  if (not (left <= right)) {
    throw std::invalid_argument(absl::StrFormat(
        "Invalid domain: need left <= right, but left=%f and right=%f", left, right));
  }

  auto m = (right_height - left_height) / (right - left);
  auto c = (left_height - m * left);
  auto y_at = [m, c](double x) { return m*x + c; };
  
  // Focus on overlap between [left,right] and [start,end]
  if (left > staircase.x_end() || right < staircase.x_start()) return;  // Doesn't overlap domain
  if (left < staircase.x_start()) {
    left = staircase.x_start();
    left_height = y_at(left);
  }
  if (right > staircase.x_end()) {
    right = staircase.x_end();
    right_height = y_at(right);
  }
  if (left == right) return;  // Empty trapezoid
  assert(staircase.x_start() <= left && left < right && right <= staircase.x_end());
  
  auto cell_start = cell_for_lbound(staircase, left);
  auto cell_end = cell_for_ubound(staircase, right);

  if (cell_start == cell_end) {
    staircase.at_cell(cell_start) += 0.5 * (left_height + right_height) * (right - left) / staircase.cell_size();
  } else {

    // Handle partial overlaps with first and last cell
    auto iter = staircase.begin() + cell_start;
    auto first_cell_ubound = cell_ubound(staircase, cell_start);
    auto first_cell_width = first_cell_ubound - left;
    auto first_cell_avg_height = 0.5 * (y_at(left) + y_at(first_cell_ubound));
    *iter += first_cell_avg_height * first_cell_width / staircase.cell_size();
    ++iter;

    auto end = staircase.begin() + cell_end;
    auto last_cell_lbound = cell_lbound(staircase, cell_end);
    auto last_cell_width = right - cell_lbound(staircase, cell_end);
    auto last_cell_avg_height = 0.5 * (y_at(last_cell_lbound) + y_at(right));
    *end += last_cell_avg_height * last_cell_width / staircase.cell_size();

    // The middle is easy
    auto cell_lbound = first_cell_ubound;
    while (iter != end) {
      auto cell_ubound = cell_lbound + staircase.cell_size();
      auto avg_height = 0.5 * (y_at(cell_lbound) + y_at(cell_ubound));  // height at cell midpoint
      *iter += avg_height;
      ++iter;
      cell_lbound = cell_ubound;  // on to the next cell
    }
  }
}

}  // namespace delphy
