#include "dates.h"

#include <stdexcept>
#include <iomanip>
#include <iostream>
#include <sstream>

#include <absl/time/time.h>

namespace delphy {

// 2020-01-01 is day 0 (used to be day 1, watch out!)
static auto k_epoch = absl::FromCivil(absl::CivilDay(2020, 1, 1), absl::UTCTimeZone());

auto parse_iso_date(std::string_view iso_date_str) -> double {
  auto d = absl::Time{};
  auto err = std::string{};
  if (not absl::ParseTime("%Y-%m-%d", iso_date_str, &d, &err)) {
    throw std::runtime_error(absl::StrFormat(
        "Badly formatted ISO date: %s (error: %s)", iso_date_str, err));
  }
  return std::round(absl::ToDoubleHours(d - k_epoch) / 24.0);
}

auto to_iso_date(double t) -> std::string {
  auto abs_time = absl::Time{k_epoch + absl::Hours(t * 24.0 + 1)};  // +1 to avoid roundoff
  return absl::FormatTime("%Y-%m-%d", abs_time, absl::UTCTimeZone());
}

// Assign dates to leaf and inner nodes of the tree in a very dumb way:
// * For leaf nodes, keep times read in from import
// * For inner nodes, estimate time from times of children & their mutations
//   (asuming ~13 days per mutation), with about ~1 day extra leeway
auto pseudo_date(Phylo_tree& tree, absl::BitGenRef bitgen) -> void {
  for (const auto& node : post_order_traversal(tree)) {
    if (tree.at(node).is_inner_node()) {
      CHECK(not tree.at(node).children.empty());
      auto left_child = tree.at(node).left_child();
      auto right_child = tree.at(node).right_child();

      // ~1 mutation / 13 days for SARS-CoV-2
      auto est_t_left = tree.at(left_child).t - std::ssize(tree.at(left_child).mutations) * 13.0;
      auto est_t_right = tree.at(right_child).t - std::ssize(tree.at(right_child).mutations) * 13.0;
      
      tree.at(node).t = std::min(est_t_left, est_t_right) - absl::Uniform(bitgen, 0.5, 1.5);
    }
  }
}

}  // namespace delphy
