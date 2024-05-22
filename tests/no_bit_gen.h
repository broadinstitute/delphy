#ifndef DELPHY_NO_BIT_GEN_H_
#define DELPHY_NO_BIT_GEN_H_

#include <absl/random/bit_gen_ref.h>

namespace delphy {

struct No_bit_gen {
  using result_type = uint64_t;

  static constexpr auto min() -> result_type { return std::numeric_limits<result_type>::min(); }
  static constexpr auto max() -> result_type { return std::numeric_limits<result_type>::max(); }

  auto operator()() const -> result_type { throw std::runtime_error("No random number generation here!"); }
};

static_assert(absl::random_internal::is_urbg<No_bit_gen>::value);

}  // namespace delphy

#endif // DELPHY_NO_BIT_GEN_H_
