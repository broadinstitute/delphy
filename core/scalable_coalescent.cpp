#include "scalable_coalescent.h"

#include <cmath>

namespace delphy {

Scalable_coalescent_prior::Scalable_coalescent_prior(
    const Pop_model& pop_model,
    Node_index num_nodes,
    double t_ref,
    double t_step)
    : pop_model_{pop_model},
      node_infos_(num_nodes, {.t = t_ref, .is_tip = false}),
      k_bars_{},
      popsize_bars_{},
      t_ref_{t_ref},
      t_step_{t_step} {
  cells_before_t_ref_ = 0;

  // Conceptually, all tips and coalescences start out at t_ref
}

auto Scalable_coalescent_prior::reset(double t_step) -> void {
  // Send all tips and coalescences back to t_ref and clear all accumulated values
  t_step_ = t_step;
  for (auto& node_info : node_infos_) {
    node_info.t = t_ref_;
  }
  k_bars_.clear();
  popsize_bars_.clear();
  cells_before_t_ref_ = 0;
}

auto Scalable_coalescent_prior::cell_for(double t) const -> int {
  return static_cast<int>(std::floor((t - t_ref_) / t_step_)) + cells_before_t_ref_;
}

auto Scalable_coalescent_prior::assert_space(double t) -> void {
  if (estd::is_debug_enabled) {
    auto cell = cell_for(t);
    auto tot_cells = std::ssize(k_bars_);
    
    CHECK_GE(cell, 0);
    CHECK_LT(cell, tot_cells);
  }
}

auto Scalable_coalescent_prior::ensure_space(double t) -> void {
  auto cell = cell_for(t);
  auto tot_cells = std::ssize(k_bars_);

  if (cell < 0) {
    auto num_cells_to_prepend = -cell;
    auto t_max_new = cell_lbound(0);
    auto cum_popsize_max = pop_model_.cum_pop_at_time(t_max_new);
    for (auto i = 0; i != num_cells_to_prepend; ++i) {
      k_bars_.push_front(1.0);   // The ancestors of the root node can go on forever

      auto t_min_new = t_max_new - t_step_;
      auto cum_popsize_min = pop_model_.cum_pop_at_time(t_min_new);
      popsize_bars_.push_front((cum_popsize_max - cum_popsize_min) / t_step_);

      t_max_new = t_min_new;
      cum_popsize_max = cum_popsize_min;
    }
    cells_before_t_ref_ += num_cells_to_prepend;

  } else if (cell >= tot_cells) {
    auto num_cells_to_append = cell - tot_cells + 1;
    auto t_min_new = cell_ubound(tot_cells - 1);
    auto cum_popsize_min = pop_model_.cum_pop_at_time(t_min_new);
    for (auto i = 0; i != num_cells_to_append; ++i) {
      k_bars_.push_back(0.0);

      auto t_max_new = t_min_new + t_step_;
      auto cum_popsize_max = pop_model_.cum_pop_at_time(t_max_new);
      popsize_bars_.push_back((cum_popsize_max - cum_popsize_min) / t_step_);
      t_min_new = t_max_new;
      cum_popsize_min = cum_popsize_max;
    }
  }

  assert_space(t);
}

auto Scalable_coalescent_prior::add_interval(double t_start, double t_end, double delta_k) -> void {
  if (t_start > t_end) {
    std::swap(t_start, t_end);
  }

  ensure_space(t_start);
  ensure_space(t_end);

  auto cell_start = cell_for(t_start);
  auto cell_end = cell_for(t_end);

  if (cell_start == cell_end) {
    k_bars_[cell_start] += delta_k * (t_end - t_start) / t_step_;
  } else {
    // Handle partial overlaps with first and last cell
    auto iter = k_bars_.begin() + cell_start;
    *iter += delta_k * (cell_ubound(cell_start) - t_start) / t_step_;
    ++iter;

    auto end = k_bars_.begin() + cell_end;
    *end += delta_k * (t_end - cell_lbound(cell_end)) / t_step_;

    // The middle is easy
    while (iter != end) {
      *iter += delta_k;
      ++iter;
    }
  }
}

auto Scalable_coalescent_prior::displace_tip(Node_index node, double new_t) -> void {
  DCHECK(node_infos_.at(node).is_tip);
  auto old_t = node_infos_[node].t;
  if (old_t <= new_t) {
    add_interval(old_t, new_t, +1.0);
  } else {
    add_interval(new_t, old_t, -1.0);
  }
  node_infos_[node].t = new_t;
}

auto Scalable_coalescent_prior::displace_coalescence(Node_index node, double new_t) -> void {
  DCHECK(not node_infos_.at(node).is_tip);
  auto old_t = node_infos_[node].t;
  if (old_t <= new_t) {
    add_interval(old_t, new_t, -1.0);
  } else {
    add_interval(new_t, old_t, +1.0);
  }
  node_infos_[node].t = new_t;
}

auto Scalable_coalescent_prior::pop_model_changed() -> void {
  auto tot_cells = std::ssize(k_bars_);
  auto t_min = cell_lbound(0);
  auto cum_popsize_min = pop_model_.cum_pop_at_time(t_min);

  auto iter_popsize_bar = popsize_bars_.begin();
  for (auto cell = 0; cell != tot_cells; ++cell) {
    auto t_max = t_min + t_step_;
    auto cum_popsize_max = pop_model_.cum_pop_at_time(t_max);

    *iter_popsize_bar = (cum_popsize_max - cum_popsize_min) / t_step_;

    t_min = t_max;
    cum_popsize_min = cum_popsize_max;
    ++iter_popsize_bar;
  }
}

auto Scalable_coalescent_prior::calc_log_prior() const -> double {
  // -[int dt kbar(tree) * (kbar(tree) - 1) / (2 * Nbar(tree))] - sum_i log(N(t_i))
  auto result = 0.0;

  auto tot_cells = std::ssize(k_bars_);
  auto t_min = cell_lbound(0);
  auto cum_popsize_min = pop_model_.cum_pop_at_time(t_min);
  auto iter_k_bar = k_bars_.begin();
  for (auto cell = 0; cell != tot_cells; ++cell) {
    auto t_max = t_min + t_step_;
    auto cum_popsize_max = pop_model_.cum_pop_at_time(t_max);

    auto kbar = *iter_k_bar;
    auto popsize_bar = (cum_popsize_max - cum_popsize_min) / t_step_;

    result -= t_step_ * kbar * (kbar - 1) / (2.0 * popsize_bar);

    t_min = t_max;
    cum_popsize_min = cum_popsize_max;
    ++iter_k_bar;
  }

  for (const auto& node_info : node_infos_) {
    if (not node_info.is_tip) {
      result -= std::log(pop_model_.pop_at_time(node_info.t));
    }
  }

  return result;
}

auto Scalable_coalescent_prior::calc_delta_log_prior_after_displace_coalescence(int i, double new_t) -> double {
  // Gist: recap the loop from calc_log_prior, but only run from min(old_t, new_t) to max(old_t, new_t)
  auto delta_log_prior = 0.0;

  assert(not node_infos_.at(i).is_tip);
  auto old_t = node_infos_.at(i).t;

  if (old_t == new_t) { return 0.0; }
  auto adding_lineages = new_t < old_t;

  auto min_t = std::min(old_t, new_t);
  auto max_t = std::max(old_t, new_t);

  ensure_space(old_t);
  ensure_space(new_t);

  auto cell_start = cell_for(min_t);
  auto cell_end = cell_for(max_t);

  if (cell_start == cell_end) {
    auto kbar = k_bars_[cell_start];
    auto popsize_bar = popsize_bars_[cell_start];
    auto delta_kbar = (adding_lineages ? +1.0 : -1.0) * (max_t - min_t) / t_step_;
    delta_log_prior -=
        t_step_ * ((kbar + delta_kbar) * (kbar + delta_kbar - 1) - kbar * (kbar - 1)) / (2 * popsize_bar);

  } else {
    // Handle ends specially
    auto iter_k_bar = k_bars_.begin() + cell_start;
    auto iter_popsize_bar = popsize_bars_.begin() + cell_start;
    auto kbar = *iter_k_bar;
    auto popsize_bar = *iter_popsize_bar;
    auto delta_kbar = (adding_lineages ? +1.0 : -1.0) * (cell_ubound(cell_start) - min_t) / t_step_;

    delta_log_prior -=
        t_step_ * ((kbar + delta_kbar) * (kbar + delta_kbar - 1) - kbar * (kbar - 1)) / (2 * popsize_bar);

    ++iter_k_bar;
    ++iter_popsize_bar;

    for (auto cell = cell_start + 1; cell < cell_end; ++cell) {
      kbar = *iter_k_bar;
      popsize_bar = *iter_popsize_bar;
      delta_kbar = (adding_lineages ? +1.0 : -1.0);
      delta_log_prior -=
          t_step_ * ((kbar + delta_kbar) * (kbar + delta_kbar - 1) - kbar * (kbar - 1)) / (2 * popsize_bar);
      ++iter_k_bar;
      ++iter_popsize_bar;
    }

    kbar = *iter_k_bar;
    popsize_bar = *iter_popsize_bar;
    delta_kbar = (adding_lineages ? +1.0 : -1.0) * (max_t - cell_lbound(cell_end)) / t_step_;

    delta_log_prior -=
        t_step_ * ((kbar + delta_kbar) * (kbar + delta_kbar - 1) - kbar * (kbar - 1)) / (2 * popsize_bar);
  }

  // Add factor from 1/N(tree)
  delta_log_prior -= std::log(pop_model_.pop_at_time(new_t) / pop_model_.pop_at_time(old_t));

  return delta_log_prior;
}

}  // namespace delphy
