#ifndef DELPHY_ESTD_H_
#define DELPHY_ESTD_H_

#include <algorithm>
#include <charconv>
#include <numeric>
#include <stdexcept>
#include <string>

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

}  // namespace estd

#endif // DELPHY_ESTD_H_
