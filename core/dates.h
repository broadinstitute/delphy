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

// Convert a time in "days since 1-Jan-2020" (no ambiguities w.r.t. leap years)
// to BEAST's linear years (see https://beast.community/tip_dates)
// Example: to_linear_year(parse_iso_date("2024-02-01")) == ~ 2024.08333 = 2024 + 1./12.
auto to_linear_year(double t) -> double;

// A month is represented as midnight at the beginning of its first day to midnight at the end of its last day
auto parse_iso_month(std::string_view iso_month_str) -> std::pair<double, double>;

// A year is represented as midnight at the beginning of its first day to midnight at the end of its last day
auto parse_iso_year(std::string_view iso_year_str) -> std::pair<double, double>;

auto pseudo_date(Phylo_tree& tree, absl::BitGenRef bitgen) -> void;

}  // namespace delphy

#endif // DELPHY_DATES_H_
