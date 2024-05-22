#include "coalescent.h"

#include <algorithm>
#include <cmath>

#include "estd.h"

namespace delphy {

Coalescent_prior::Coalescent_prior(const Pop_model *pop_model, int num_tips)
    : pop_model_{pop_model} {
  for (int i = 0; i != num_tips; ++i) {
    tip_times_.insert(0.0);
  }
  for (int i = 0; i != num_tips - 1; ++i) {
    coalescence_times_.insert(0.0);
  }
}

auto Coalescent_prior::displace_tip(double old_t, double new_t) -> void {
  auto it = tip_times_.find(old_t);
  if (it == tip_times_.end()) {
    throw std::out_of_range{"old_t not a tip time"};
  }
  tip_times_.replace(it, new_t);
}

auto Coalescent_prior::displace_coalescence(double old_t, double new_t) -> void {
  auto it = coalescence_times_.find(old_t);
  if (it == coalescence_times_.end()) {
    throw std::out_of_range{"old_t not a coalescence time"};
  }
  coalescence_times_.replace(it, new_t);
}

auto Coalescent_prior::count_tips_between(double t_min_closed, double t_max_open) const -> int {
  auto rank_t_min = estd::as_signed(tip_times_.lower_bound_rank(t_min_closed));
  auto rank_t_max = estd::as_signed(tip_times_.lower_bound_rank(t_max_open));
  return rank_t_max - rank_t_min;
}

auto Coalescent_prior::count_coalescences_between(double t_min_closed, double t_max_open) const -> int {
  auto rank_t_min = estd::as_signed(coalescence_times_.lower_bound_rank(t_min_closed));
  auto rank_t_max = estd::as_signed(coalescence_times_.lower_bound_rank(t_max_open));
  return rank_t_max - rank_t_min;
}

auto Coalescent_prior::calc_log_prior() const -> double {
  auto coal_it = coalescence_times_.cbegin();
  auto coal_end = coalescence_times_.cend();
  auto tip_it = tip_times_.cbegin();
  auto tip_end = tip_times_.cend();

  auto result = 0.0;

  // Walk through merged list of coalescences and tips to build up factors in calc_log_prior

  auto prev_t = -std::numeric_limits<double>::max();
  auto prev_I = 0.0;  // actual value doesn't matter
  auto active_lineages = 1;
  while (coal_it != coal_end || tip_it != tip_end) {

    auto next_is_coal = (coal_it != coal_end && (tip_it == tip_end || *coal_it <= *tip_it));
    //bool next_is_tip = !next_is_coal;  // (coal_it == coal_end || (tip_it != tip_end && *tip_it < *coal_t)

    // next_is_coal => coal_it != coal_end && next_t == *coal_it
    // next_is_tip => tip_it != tip_end && next_t == *tip_it
    auto next_t = next_is_coal ? *coal_it : *tip_it;

    if (prev_t != next_t) {
      auto next_I = pop_model_->intensity_at_time(next_t);
      auto k = static_cast<double>(active_lineages);
      result -= (k * (k - 1)) / 2 * (next_I - prev_I);
      prev_t = next_t;
      prev_I = next_I;
    }

    if (next_is_coal) {
      // Next thing is a coalescence => 1 more active lineage
      ++active_lineages;
      result -= std::log(pop_model_->pop_at_time(next_t));
      ++coal_it;
    } else {
      // Next thing is a tip => 1 fewer active lineage
      if (active_lineages == 0) {
        throw std::out_of_range{"Number of active lineages dropped below 0?"};
      }
      --active_lineages;
      ++tip_it;
    }
  }

  // Done!
  return result;
}

auto Coalescent_prior::calc_delta_log_prior_after_displace_coalescence(double old_t, double new_t) const -> double {
  // Gist: recap the loop from calc_log_prior, but only run from min(old_t, new_t) to max(old_t, new_t)

  if (old_t == new_t) { return 0.0; }
  auto adding_lineages = new_t < old_t;

  auto min_t = std::min(old_t, new_t);
  auto max_t = std::max(old_t, new_t);

  // Recap situation up to but not including min_t
  auto coal_it = coalescence_times_.lower_bound(min_t);
  auto coal_end = coalescence_times_.end();
  auto tip_it = tip_times_.lower_bound(min_t);
  auto tip_end = tip_times_.end();
  auto active_lineages = 1 + coalescence_times_.rank(coal_it) - tip_times_.rank(tip_it);
  auto prev_t = min_t;
  auto prev_I = pop_model_->intensity_at_time(min_t);

  auto delta_log_prior = 0.0;

  while (coal_it != coal_end || tip_it != tip_end) {

    auto next_is_coal = (coal_it != coal_end && (tip_it == tip_end || *coal_it <= *tip_it));
    //auto next_is_tip = !next_is_coal;  // (coal_it == coal_end || (tip_it != tip_end && *tip_it < *coal_t)

    // next_is_coal => coal_it != coal_end && next_t == *coal_it
    // next_is_tip => tip_it != tip_end && next_t == *tip_it
    auto next_t = std::min(max_t, next_is_coal ? *coal_it : *tip_it);

    if (prev_t != next_t) {
      auto next_I = pop_model_->intensity_at_time(next_t);
      auto k_old = static_cast<double>(active_lineages);
      // k_new = k_old + (adding_lineages ? 1 : -1)
      // result -= (k_new * (k_new - 1)) / 2 * (next_I - prev_I)
      //           - (k_old * (k_old - 1)) / 2 * (next_I - prev_I)
      // so
      auto delta_binom = adding_lineages ? k_old : -(k_old - 1);
      delta_log_prior -= delta_binom * (next_I - prev_I);
      prev_t = next_t;
      prev_I = next_I;
    }

    if (next_t == max_t) {
      // Stop!  Don't look at coalescences & tips >= max_t
      break;
    }

    if (next_is_coal) {
      // Next thing is a coalescence => 1 more active lineage
      ++active_lineages;
      ++coal_it;
    } else {
      // Next thing is a tip => 1 fewer active lineage
      if (active_lineages == 0) {
        throw std::out_of_range{"Number of active lineages dropped below 0?"};
      }
      --active_lineages;
      ++tip_it;
    }
  }

  // Add factor from 1/N(tree)
  delta_log_prior -= std::log(pop_model_->pop_at_time(new_t) / pop_model_->pop_at_time(old_t));

  return delta_log_prior;
}

}  // namespace delphy
