#ifndef DELPHY_TREE_PROBER_H_
#define DELPHY_TREE_PROBER_H_

#include <numeric>

#include "estd.h"
#include "staircase.h"
#include "pop_model.h"

namespace delphy {

// Tree_prober calculates where on a tree a probe sample at time t will coalesce
class Tree_prober {
private:
  const Staircase_family* branch_counts_by_category_;
  const Pop_model* pop_model_;

  Staircase_family p_;
  
public:
  
  Tree_prober(const Staircase_family& branch_counts_by_category,
              int cells_to_skip,
              const Pop_model& pop_model)
      : Tree_prober(branch_counts_by_category, cells_to_skip, pop_model, std::vector(branch_counts_by_category.num_members(), 0.0)) {}
  
  Tree_prober(const Staircase_family& branch_counts_by_category,
              int cells_to_skip,
              const Pop_model& pop_model,
              std::vector<double> p_initial)
      : branch_counts_by_category_{&branch_counts_by_category},
        pop_model_{&pop_model},
        p_{
          branch_counts_by_category.num_members(),
          cell_lbound(branch_counts_by_category, cells_to_skip),
          branch_counts_by_category.x_end(),
          branch_counts_by_category.num_cells() - cells_to_skip} {
    
    if (std::ssize(p_initial) != num_categories()) {
      throw std::invalid_argument(absl::StrFormat(
          "Invalid vector of starting probabilities: there are %d values but %d categories",
          std::ssize(p_initial), num_categories()));
    }
    for (auto p : p_initial) {
      if (p < 0.0) {
        throw std::out_of_range(absl::StrFormat("Starting probabilities can't be negative (%f)", p));
      }
    }
    if (auto total_p_initial = estd::ranges::sum(p_initial);
        total_p_initial > (1.0 + 1e-6)) {  // Tiny margin for roundoff...
      throw std::out_of_range(absl::StrFormat(
          "Starting probabilities can't add up to more than 1 (%f)", total_p_initial));
    }
    
    
    // At every iteration in the following loop, we're calculating the probability that
    // a probe sample at t_ubound will coalescese with a branch of some category in the current cell,
    // which spans [t_lbound, t_ubound).  If it does, then we calculate the probability that the
    // branch it merges into has a given category.  If it doesn't then the problem reduces to that
    // of a probe sample at t_lbound.  At times before t_start(), a probe sample is the only branch, so
    // it has nothing to coalesce into.

    auto p_before = std::vector{std::move(p_initial)};
    for (auto in_cell = 0; in_cell != branch_counts_by_category.num_cells(); ++in_cell) {
      auto t_lbound = cell_lbound(branch_counts_by_category, in_cell);
      auto t_ubound = cell_ubound(branch_counts_by_category, in_cell);
      auto intensity_over_cell = pop_model.intensity_integral(t_lbound, t_ubound);

      // Rate of coalescence: k(t) / N(t) => Prob of coalescense = 1 - exp(-int_cell dt' k(t') / N(t'))
      auto total_branches = 0.0;
      for (const auto& branches_in_cat : branch_counts_by_category) {
        total_branches += branches_in_cat.at_cell(in_cell);
      }
      auto p_coalesce = 1.0 - std::exp(-total_branches * intensity_over_cell);
      
      for (auto cat = 0; cat != num_categories(); ++cat) {
        const auto& branches_in_cat = branch_counts_by_category[cat];

        auto p_coalesce_cat =
            total_branches == 0.0
            ? 0.0
            : p_coalesce * (branches_in_cat.at_cell(in_cell) / total_branches);  // total_branches: double
        
        auto p_lbound = p_before[cat];
        auto p_ubound = p_coalesce_cat + (1.0 - p_coalesce) * p_lbound;

        if (in_cell >= cells_to_skip) {
          auto out_cell = in_cell - cells_to_skip;
          p_[cat].at_cell(out_cell) = p_ubound;
        }
        p_before[cat] = p_ubound;
      }
    }
  }

  auto branch_counts_by_category() const -> const Staircase_family& { return *branch_counts_by_category_; }
  auto pop_model() const -> const Pop_model& { return *pop_model_; }

  auto num_categories() const -> int { return branch_counts_by_category_->num_members(); }
  auto t_start() const -> double { return p_.x_start(); }
  auto t_end() const -> double { return p_.x_end(); }
  auto num_cells() const -> int { return p_.num_cells(); }

  auto p() const -> const Staircase_family& { return p_; }
  auto p(int category) const -> const Staircase& { return p_[category]; }
};

}  // namespace delphy

#endif // DELPHY_TREE_PROBER_H_
