#include "dates.h"

#include <stdexcept>
#include <iomanip>
#include <iostream>
#include <sstream>

#include <absl/time/civil_time.h>

namespace delphy {

// 2020-01-01 is day 0 (used to be day 1, watch out!)
static auto k_epoch_day = absl::CivilDay(2020, 1, 1);

auto parse_iso_date(std::string_view iso_date_str) -> double {
  auto d = absl::CivilDay{};
  if (not absl::ParseCivilTime(iso_date_str, &d)) {
    throw std::runtime_error(absl::StrFormat(
        "Badly formatted ISO date: %s", iso_date_str));
  }
  return std::round(d - k_epoch_day);
}

auto to_iso_date(double t) -> std::string {
  auto d = k_epoch_day + static_cast<int>(std::floor(t));
  return absl::FormatCivilTime(d);
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
