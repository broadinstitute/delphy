#ifndef DELPHY_DATES_H_
#define DELPHY_DATES_H_

#include <string>

#include "phylo_tree.h"

namespace delphy {

// In this whole code, we measure time in units of "days since 1-Jan-2020", i.e.:
//
//  2020-01-01 => 0
//  2020-01-02 => 1
//  2020-01-03 => 2
//   ...

auto parse_iso_date(std::string_view iso_date_str) -> double;
auto to_iso_date(double t) -> std::string;

auto pseudo_date(Phylo_tree& tree, absl::BitGenRef bitgen) -> void;

}  // namespace delphy

#endif // DELPHY_DATES_H_
