#include "coalescent.h"

#include <algorithm>
#include <cmath>

#include <absl/log/check.h>

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
  auto k = 1;  // number of active lineages just before the next event
  while (coal_it != coal_end || tip_it != tip_end) {

    auto next_is_coal = (coal_it != coal_end && (tip_it == tip_end || *coal_it <= *tip_it));
    //bool next_is_tip = !next_is_coal;  // (coal_it == coal_end || (tip_it != tip_end && *tip_it < *coal_t)

    // next_is_coal => coal_it != coal_end && next_t == *coal_it
    // next_is_tip => tip_it != tip_end && next_t == *tip_it
    auto next_t = next_is_coal ? *coal_it : *tip_it;

    if (k >= 2) {
      CHECK_NE(prev_t, -std::numeric_limits<double>::max());  // k == 1 in the first iteration
      result -= (k * (k - 1)) / 2 * pop_model_->intensity_integral(prev_t, next_t);
    }
    prev_t = next_t;

    if (next_is_coal) {
      // Next thing is a coalescence => 1 more active lineage
      ++k;
      result -= std::log(pop_model_->pop_at_time(next_t));
      ++coal_it;
    } else {
      // Next thing is a tip => 1 fewer active lineage
      CHECK_GT(k, 0);
      --k;
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
  auto k = 1 + static_cast<int>(coalescence_times_.rank(coal_it)) - static_cast<int>(tip_times_.rank(tip_it));  // # of active lineages
  auto prev_t = min_t;

  auto delta_log_prior = 0.0;

  while (coal_it != coal_end || tip_it != tip_end) {

    auto next_is_coal = (coal_it != coal_end && (tip_it == tip_end || *coal_it <= *tip_it));
    //auto next_is_tip = !next_is_coal;  // (coal_it == coal_end || (tip_it != tip_end && *tip_it < *coal_t)

    // next_is_coal => coal_it != coal_end && next_t == *coal_it
    // next_is_tip => tip_it != tip_end && next_t == *tip_it
    auto next_t = std::min(max_t, next_is_coal ? *coal_it : *tip_it);

    if (prev_t != next_t) {
      // k_old = k;
      // k_new = k_old + (adding_lineages ? 1 : -1)
      // result -= (k_new * (k_new - 1)) / 2 * I(prev_t, next_t)
      //           - (k_old * (k_old - 1)) / 2 * I(prev_t, next_t)
      // so
      auto delta_binom = adding_lineages ? k : -(k - 1);
      delta_log_prior -= delta_binom * pop_model_->intensity_integral(prev_t, next_t);
      prev_t = next_t;
    }

    if (next_t == max_t) {
      // Stop!  Don't look at coalescences & tips >= max_t
      break;
    }

    if (next_is_coal) {
      // Next thing is a coalescence => 1 more active lineage
      ++k;
      ++coal_it;
    } else {
      // Next thing is a tip => 1 fewer active lineage
      CHECK_GT(k, 0);
      --k;
      ++tip_it;
    }
  }

  // Add factor from 1/N(tree)
  delta_log_prior -= std::log(pop_model_->pop_at_time(new_t) / pop_model_->pop_at_time(old_t));

  return delta_log_prior;
}

}  // namespace delphy
