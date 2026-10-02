#ifndef DELPHY_ESTD_H_
#define DELPHY_ESTD_H_

#include <algorithm>
#include <charconv>
#include <numeric>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>

#include "absl/strings/str_format.h"

// estd contains extensions to std that probably should have been there
namespace estd {

// Debug support (we try very hard not to use conditional compilation)
#ifndef NDEBUG
inline constexpr bool is_debug_enabled = true;
#else
inline constexpr bool is_debug_enabled = false;
#endif

// See https://herbsutter.com/2013/06/13/gotw-93-solution-auto-variables-part-2/
template<typename T> auto as_signed(T t){ return std::make_signed_t<T>(t); }
template<typename T> auto as_unsigned(T t){ return std::make_unsigned_t<T>(t); }

namespace ranges {

template<typename Range>
auto sum(Range& r) { return std::accumulate(std::ranges::begin(r), std::ranges::end(r), 0.0, std::plus{}); }

template<std::ranges::range R>
auto to_vec(R&& range) {
  auto result = std::vector<std::ranges::range_value_t<R>>{};
  for (const auto& elem : range) {
    result.push_back(elem);
  }
  return result;
}

}  // namespace ranges

// Simple helper class for viewing a 2D array that's been flattened into a row-major-order 1D array as a 2D array again
// (stand-in for C++23's std::mdspan)
template<typename T>
struct View_2d {
  std::span<T> data;
  int rows;
  int cols;

  View_2d(std::span<T> data_in, int rows_in, int cols_in)
      : data{data_in}, rows{rows_in}, cols{cols_in} {
    if (std::ssize(data) != rows * cols) {
      throw std::invalid_argument(absl::StrFormat(
          "View_2d: data has %d elements, but should have %d x %d = %d",
          std::ssize(data), rows, cols, rows * cols));
    }
  }

  auto operator()(int i, int j) const -> T& {
    return data[i * cols + j];  // element
  }

  auto operator()(int i) const -> std::span<T> {
    return data.subspan(i * cols, cols);  // a full row
  }
};

// Deduce element type from a contiguous range, e.g., `View_2d{vec, rows, cols}`
template<std::ranges::contiguous_range R>
View_2d(R&&, int, int) -> View_2d<std::remove_reference_t<std::ranges::range_reference_t<R>>>;

// overloaded template for std::variant (from https://en.cppreference.com/w/cpp/utility/variant/visit)
template<class... Ts>
struct overloaded : Ts... { using Ts::operator()...; };

}  // namespace estd

#endif // DELPHY_ESTD_H_
