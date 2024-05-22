#include "very_scalable_coalescent.h"

#include <iostream>

namespace delphy {

static auto square(double x) -> double {
  return x * x;
}

static auto cell_for(double t, double t_ref, double t_step) -> int {
  return static_cast<int>(std::floor((t_ref - t) / t_step));
}

static auto cell_ubound(int cell, double t_ref, double t_step) -> double {
  return t_ref - t_step * cell;
}

static auto cell_lbound(int cell, double t_ref, double t_step) -> double {
  return cell_ubound(cell, t_ref, t_step) - t_step;
}

template<typename T>
static auto assert_space(double t, double t_ref, double t_step, const std::vector<T>& k) -> void {
  if (estd::is_debug_enabled) {
    auto cell = cell_for(t, t_ref, t_step);
    auto tot_cells = std::ssize(k);
    
    CHECK_GE(cell, 0);
    CHECK_LT(cell, tot_cells);
  }
}

static auto add_interval(
    double t_start,
    double t_end,
    double delta_k,
    std::vector<double>& k,
    double t_ref,
    double t_step)
    -> void {
  
  if (t_start < t_end) {  // Cell indices grow as time decreases
    std::swap(t_start, t_end);
  }

  assert_space(t_start, t_ref, t_step, k);
  auto cell_start = cell_for(t_start, t_ref, t_step);

  assert_space(t_end, t_ref, t_step, k);
  auto cell_end = cell_for(t_end, t_ref, t_step);

  if (cell_start == cell_end) {
    k[cell_start] += delta_k * (t_start - t_end) / t_step;
  } else {
    // Handle partial overlaps with first and last cell
    auto iter = k.begin() + cell_start;
    *iter += delta_k * (t_start - cell_lbound(cell_start, t_ref, t_step)) / t_step;
    ++iter;

    auto end = k.begin() + cell_end;
    *end += delta_k * (cell_ubound(cell_end, t_ref, t_step) - t_end) / t_step;

    // The middle is easy
    while (iter != end) {
      *iter += delta_k;
      ++iter;
    }
  }
}

auto make_very_scalable_coalescent_prior_parts(
    const std::vector<const Phylo_tree*>& subtrees,
    int root_partition_index,
    const Pop_model& pop_model,
    std::vector<std::mt19937>& prngs,
    double t_step)
    -> std::vector<Very_scalable_coalescent_prior_part> {

  if (subtrees.empty()) {
    return {};
  }

  struct Part_info {
    int part_index;
    const Phylo_tree* subtree;
    double t_min;
    double t_max;
    std::vector<double> k_bar_p;
    std::vector<double> k_twiddle_bar_p;
  };
  auto infos = std::vector<Part_info>(std::ssize(subtrees));

  for (auto i = 0; i != std::ssize(subtrees); ++i) {
    infos[i].part_index = i;
    infos[i].subtree = subtrees[i];
    infos[i].t_min = std::numeric_limits<double>::max();
    infos[i].t_max = std::numeric_limits<double>::min();
  }

  // Find t_min and t_max of every partition
  for (auto& info : infos) {
    for (const auto& node : index_order_traversal(*info.subtree)) {
      info.t_min = std::min(info.t_min, info.subtree->at(node).t);
      info.t_max = std::max(info.t_max, info.subtree->at(node).t);
    }
  }

  // Now min and max over the whole tree
  auto all_t_min = std::ranges::min(infos, {}, [](const auto& part_info) { return part_info.t_min; }).t_min;
  auto all_t_max = std::ranges::max(infos, {}, [](const auto& part_info) { return part_info.t_max; }).t_max;

  // Everything uses all_t_max as the reference time and grows time towards the past
  auto t_ref = all_t_max;

  // Count total number of cells
  auto num_cells = cell_for(all_t_min, t_ref, t_step) + 1;

  // Everything uses all_t_max as the reference time and grows time towards the past.
  // All parts keep part-specific per-cell information, but the number of cells in each part varies
  // depending on the minimum time in that partition (and, for the root partition, can grow dynamically
  // whenever the root is moved further back in time)

  // Count for each cell how many parts are active and allocate per-part structures
  auto num_active_parts = std::vector<int>(num_cells, 0);
  for (auto& info : infos) {
    auto first_cell = cell_for(info.t_max, t_ref, t_step);
    auto last_cell = cell_for(info.t_min, t_ref, t_step);

    for (auto cell = first_cell; cell <= last_cell; ++cell) {
      num_active_parts[cell] += 1;
    }

    auto num_cells_p = last_cell + 1;
    assert(num_cells_p <= num_cells);
    info.k_bar_p.assign(num_cells_p, 0.0);
    info.k_twiddle_bar_p.assign(num_cells_p, 0.0);
  }

  // Make sure there are no cells that are inactive at the beginning due to roundoff in cell boundary calculations
  CHECK_NE(num_active_parts.back(), 0);

  // Calculate k_bar_p for each part
  for (auto& info : infos) {
    for (const auto& node : index_order_traversal(*info.subtree)) {
      if (node != info.subtree->root) {
        add_interval(info.subtree->at_parent_of(node).t, info.subtree->at(node).t, +1.0,
                     info.k_bar_p, t_ref, t_step);
      }
    }
  }

  // Add the last interval from the root node to the earliest time we track
  auto& root_info = infos[root_partition_index];
  add_interval(cell_lbound(num_cells - 1, t_ref, t_step) + 0.0001*t_step,
               root_info.subtree->at_root().t,
               +1.0, root_info.k_bar_p, t_ref, t_step);

  // Allocate and calculate k_bar
  auto k_bar = std::vector<double>(num_cells, 0.0);
  for (auto& info : infos) {
    for (auto i = 0; i < std::ssize(info.k_bar_p); ++i) {
      k_bar[i] += info.k_bar_p[i];
    }
  }

  // Allocate and calculate popsize_bar
  auto popsize_bar = std::vector<double>(num_cells, 0.0);
  auto cum_popsize_max = pop_model.cum_pop_at_time(cell_ubound(0, t_ref, t_step));
  for (auto i = 0; i != num_cells; ++i) {
    auto cum_popsize_min = pop_model.cum_pop_at_time(cell_lbound(i, t_ref, t_step));
    popsize_bar[i] = (cum_popsize_max - cum_popsize_min) / t_step;
    cum_popsize_max = cum_popsize_min;
  }

  // Sample k_twiddle_bar_p for each part
  for (auto& info : infos) {
    for (auto i = 0; i < std::ssize(info.k_twiddle_bar_p); ++i) {
      auto mu = info.k_bar_p[i] - k_bar[i] / num_active_parts[i];
      auto sigma = std::sqrt(popsize_bar[i] / (num_active_parts[i] * t_step));
      info.k_twiddle_bar_p[i] = std::normal_distribution<double>{mu, sigma}(prngs[info.part_index]);
    }
  }

  // Allocate and calculate k_twiddle_bar
  auto k_twiddle_bar = std::vector<double>(num_cells, 0.0);
  for (auto& info : infos) {
    for (auto i = 0; i < std::ssize(info.k_twiddle_bar_p); ++i) {
      k_twiddle_bar[i] += info.k_twiddle_bar_p[i];
    }
  }

  // Package it all up into Very_scalable_coalescent_prior_parts
  auto result = std::vector<Very_scalable_coalescent_prior_part>{};
  for (auto& info : infos) {
    auto includes_tree_root = root_partition_index == info.part_index;
    result.emplace_back(pop_model, *info.subtree, prngs[info.part_index],
                        includes_tree_root, t_ref, t_step,
                        std::move(info.k_bar_p), std::move(info.k_twiddle_bar_p),
                        k_twiddle_bar, popsize_bar, num_active_parts);
  }

  return result;
}

Very_scalable_coalescent_prior_part::Very_scalable_coalescent_prior_part(
    const Pop_model& pop_model,
    const Phylo_tree& subtree,
    std::mt19937& prng,
    bool includes_tree_root,
    double t_ref,
    double t_step,
    std::vector<double> k_bar_p,
    std::vector<double> k_twiddle_bar_p,
    std::vector<double> k_twiddle_bar,
    std::vector<double> popsize_bar,
    std::vector<int> num_active_parts)
    : pop_model_{&pop_model},
      subtree_{&subtree},
      prng_{&prng},
      includes_tree_root_{includes_tree_root},
      k_bar_p_{std::move(k_bar_p)},
      k_twiddle_bar_p_{std::move(k_twiddle_bar_p)},
      k_twiddle_bar_{std::move(k_twiddle_bar)},
      popsize_bar_{std::move(popsize_bar)},
      num_active_parts_{std::move(num_active_parts)},
      t_ref_{t_ref},
      t_step_{t_step} {
}

auto Very_scalable_coalescent_prior_part::ensure_space(double t) -> void {
  if (includes_tree_root_) {
    // May need to add extra cells
    auto max_cell = cell_for(t, t_ref_, t_step_);

    for (auto i = std::ssize(popsize_bar_); i <= max_cell; ++i) {
      auto popsize_bar_i = (pop_model_->cum_pop_at_time(cell_ubound(i, t_ref_, t_step_))
                            - pop_model_->cum_pop_at_time(cell_lbound(i, t_ref_, t_step_))) / t_step_;
      popsize_bar_.push_back(popsize_bar_i);
      num_active_parts_.push_back(1);
    }

    for (auto i = std::ssize(k_bar_p_); i <= max_cell; ++i) {
      auto sigma_bar_p_i = std::sqrt(popsize_bar_[i] / t_step_);
      auto k_twiddle_bar_p_i = std::normal_distribution{0.0, sigma_bar_p_i}(*prng_);
      k_bar_p_.push_back(1.0);
      k_twiddle_bar_p_.push_back(k_twiddle_bar_p_i);
      k_twiddle_bar_.push_back(k_twiddle_bar_p_i);
    }
  }

  assert_space(t, t_ref_, t_step_, k_bar_p_);
  assert_space(t, t_ref_, t_step_, k_twiddle_bar_p_);
  assert_space(t, t_ref_, t_step_, k_twiddle_bar_);
  assert_space(t, t_ref_, t_step_, popsize_bar_);
  assert_space(t, t_ref_, t_step_, num_active_parts_);
}

auto Very_scalable_coalescent_prior_part::coalescence_displaced(double old_t, double new_t) -> void {
  ensure_space(new_t);
  if (old_t <= new_t) {
    add_interval(old_t, new_t, -1.0, k_bar_p_, t_ref_, t_step_);
  } else {
    add_interval(old_t, new_t, +1.0, k_bar_p_, t_ref_, t_step_);
  }
}

auto Very_scalable_coalescent_prior_part::calc_partial_log_prior() const -> double {
  // - int dt [(k_twiddle_bar_p(tree) - mu_bar_p(tree))^2 / (2 sigma_bar_p(tree)^2)
  //           + (k_twiddle_bar(tree) - 1/2) * k_bar_p(tree) / N_bar(tree)
  //           - k_twiddle_bar^2(tree) / 2 N_bar(tree)]  <-- normalization to make it indep. of cells where k_bar_p = 0 or 1
  // - sum_i log(N(t_i))
  //
  // where mu_bar_p(tree) = k_bar_p(tree)
  // and sum_{p active at tree} sigma_bar_p(tree)^2 = N_bar(tree)
  //
  auto result = 0.0;

  for (auto i = 0; i != std::ssize(k_bar_p_); ++i) {
    //auto mu_p_i = k_bar_p_[i];
    //auto sigma2_p_i = popsize_bar_[i] / num_active_parts_[i];
    //result -=
    //    t_step_ * (square(k_twiddle_bar_p_[i] - mu_p_i) / (2 * sigma2_p_i)
    //        + (k_twiddle_bar_[i] - 0.5) * k_bar_p_[i] / popsize_bar_[i]);

    result -= t_step_ / popsize_bar_[i] * (
        0.5 * square(k_twiddle_bar_p_[i] - k_bar_p_[i]) * num_active_parts_[i]
        + (k_twiddle_bar_[i] - 0.5) * k_bar_p_[i]
        - 0.5 * square(k_twiddle_bar_p_[i]) * num_active_parts_[i]
    );
  }

  for (const auto& node : index_order_traversal(*subtree_)) {
    if (subtree_->at(node).is_inner_node()) {
      result -= std::log(pop_model_->pop_at_time(subtree_->at(node).t));
    }
  }

  return result;
}

auto Very_scalable_coalescent_prior_part::calc_delta_partial_log_prior_after_displace_coalescence(
    double old_t,
    double new_t)
    -> double {
  
  // Gist: recap the loop from calc_partial_log_prior, but only run from min(old_t, new_t) to max(old_t, new_t)
  ensure_space(new_t);  // This may change the number of cells, which offsets the log-prior

  auto delta_log_prior = 0.0;

  if (old_t == new_t) { return 0.0; }
  auto adding_lineages = new_t < old_t;

  auto min_t = std::min(old_t, new_t);
  auto max_t = std::max(old_t, new_t);

  auto cell_start = cell_for(max_t, t_ref_, t_step_);
  auto cell_end = cell_for(min_t, t_ref_, t_step_);

  if (cell_start == cell_end) {
    auto i = cell_start;
    auto old_k_bar_p_i = k_bar_p_[i];
    auto delta_k_bar_p_i = (adding_lineages ? +1.0 : -1.0) * (max_t - min_t) / t_step_;
    auto new_k_bar_p_i = old_k_bar_p_i + delta_k_bar_p_i;

    //result -= t_step_ / popsize_bar_[i] * (
    //    0.5 * square(k_twiddle_bar_p_[i] - k_bar_p_[i]) * num_active_parts_[i]
    //        + (k_twiddle_bar_[i] - 0.5) * k_bar_p_[i]
    //);

    delta_log_prior -= t_step_ / popsize_bar_[i] * (
        0.5 * num_active_parts_[i] * (
            square(k_twiddle_bar_p_[i] - new_k_bar_p_i) -
            square(k_twiddle_bar_p_[i] - old_k_bar_p_i))
        + (k_twiddle_bar_[i] - 0.5) * (new_k_bar_p_i - old_k_bar_p_i));

  } else {
    // Handle ends specially
    auto i = cell_start;

    auto old_k_bar_p_i = k_bar_p_[i];
    auto delta_k_bar_p_i = (adding_lineages ? +1.0 : -1.0) * (max_t - cell_lbound(i, t_ref_, t_step_)) / t_step_;
    auto new_k_bar_p_i = old_k_bar_p_i + delta_k_bar_p_i;

    delta_log_prior -= t_step_ / popsize_bar_[i] * (
        0.5 * num_active_parts_[i] * (
            square(k_twiddle_bar_p_[i] - new_k_bar_p_i) -
            square(k_twiddle_bar_p_[i] - old_k_bar_p_i))
        + (k_twiddle_bar_[i] - 0.5) * (new_k_bar_p_i - old_k_bar_p_i)
    );

    for(++i; i < cell_end; ++i) {
      old_k_bar_p_i = k_bar_p_[i];
      delta_k_bar_p_i = (adding_lineages ? +1.0 : -1.0);
      new_k_bar_p_i = old_k_bar_p_i + delta_k_bar_p_i;

      delta_log_prior -= t_step_ / popsize_bar_[i] * (
          0.5 * num_active_parts_[i] * (
              square(k_twiddle_bar_p_[i] - new_k_bar_p_i) -
              square(k_twiddle_bar_p_[i] - old_k_bar_p_i))
          + (k_twiddle_bar_[i] - 0.5) * (new_k_bar_p_i - old_k_bar_p_i));
    }

    assert(i == cell_end);

    old_k_bar_p_i = k_bar_p_[i];
    delta_k_bar_p_i = (adding_lineages ? +1.0 : -1.0) * (cell_ubound(i, t_ref_, t_step_) - min_t) / t_step_;
    new_k_bar_p_i = old_k_bar_p_i + delta_k_bar_p_i;
    
    delta_log_prior -= t_step_ / popsize_bar_[i] * (
        0.5 * num_active_parts_[i] * (
            square(k_twiddle_bar_p_[i] - new_k_bar_p_i) -
            square(k_twiddle_bar_p_[i] - old_k_bar_p_i))
        + (k_twiddle_bar_[i] - 0.5) * (new_k_bar_p_i - old_k_bar_p_i));
  }

  // Add factor from 1/N(tree)
  delta_log_prior -= std::log(pop_model_->pop_at_time(new_t) / pop_model_->pop_at_time(old_t));

  return delta_log_prior;
}

}  // namespace delphy
