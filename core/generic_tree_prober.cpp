#include "generic_tree_prober.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <span>
#include <stdexcept>
#include <variant>
#include <vector>

#include <absl/log/check.h>
#include <absl/strings/str_format.h>
#include <absl/strings/str_join.h>

#include "estd.h"
#include "pop_model.h"

namespace delphy {

auto generic_probe_tree(
    std::span<const Bundle_event> bundle_events,
    int num_bundles,
    const Pop_model& pop_model,
    std::span<const double> probe_times,
    std::span<double> out_results)
    -> void {

  // Input validation
  auto num_probe_times = std::ssize(probe_times);
  auto has_finite_probe_time = false;
  for (const auto& t_i : probe_times) {
    if (std::isnan(t_i)) {
      throw std::invalid_argument(
          "Probe times should not be NaN");
    }
    if (not std::isinf(t_i)) {
      has_finite_probe_time = true;
    }
  }
  if (not probe_times.empty() && not has_finite_probe_time) {
    throw std::invalid_argument(
        "At least one probe time should be finite");
  }
  if (not std::ranges::is_sorted(probe_times)) {
    throw std::invalid_argument(absl::StrFormat(
        "Probe times should be sorted (%s)", absl::StrJoin(probe_times, ", ", absl::StreamFormatter())));
  }
  if (num_bundles < 0) {
    throw std::invalid_argument(absl::StrFormat(
        "Number of bundles shouldn't be negative: %d", num_bundles));
  }
  for (const auto& be : bundle_events) {
    if (be.target_bundle < 0 || be.target_bundle >= num_bundles) {
      throw std::invalid_argument(absl::StrFormat(
          "Event with out-of-range bundle index (%d) (expected in [0, %d))", be.target_bundle, num_bundles));
    }
    if (std::isnan(be.t)) {
      throw std::invalid_argument(
          "Event times should not be NaN");
    }
  }
  if (std::ssize(out_results) != (num_bundles * num_probe_times)) {
    throw std::invalid_argument(absl::StrFormat(
        "Incorrectly sized `out_results` span: should have %d x %d = %d elements, but has %d",
        num_bundles, num_probe_times, num_bundles * num_probe_times, std::ssize(out_results)));
  }

  // Early handling of pathological edge cases (no bundles, no probes)
  if (num_bundles * num_probe_times == 0) {
    return;
  }
  
  // We organize this calculation as a piecewise calculation between discrete events seen in order
  // of increasing time.
  //
  // At any one point, we're accumulating info about the interval (t_1, t_2) between the
  // previous (t_1) and the next (t_2) probe times:
  // 
  // * The probability that a probe at t_2 doesn't attach to the tree before reaching t_1
  // * The probability that a probe at t_2 attaches to each branch bundle b between t_1 and t_2
  //
  // Whenever we hit a probe, we use that info for interval (t_1, t_2) to compute the result for
  // (-infty, t_2).
  //
  // Structuring the calculation as above makes the per-event cost scale as the maximum
  // number of bundles that are active somewhere within an inter-probe interval.  As trees
  // grow by adding more data in the future, the per-event cost thus stays roughly
  // bounded.  The alternative would be to keep track of the attachment probability of
  // _all_ bundles at all times, and for a large tree probe where each branch is a
  // separate bundle, that makes the probing cost scale quadratically with number of tips:
  // a disaster.  The current API forces each probe to record `num_bundles` of
  // probabilities, even in cases where most probabilities are essentially 0.0.  Finding a
  // better API that's still ergonomic would be an important improvement.
  //
  // Between two events at times t_i and t_{i+1}, the population size N and the set and composition 
  // of all active bundles is fixed.  The probability that a probe at t_{i+1} attaches to a bundle
  // breaks down into two cases: (1) it attaches between t_i and t_{i+1}; (2) it doesn't, but it
  // attaches before.  In other words, for a bundle b with k_b active branches and k total active branches:
  //
  //  p_attach_b <- p * (k_b / k) + (1-p) * p_attach_b,
  //
  // where p = 1 - exp(-dt * k / N).
  //
  // In reality, we average 1/N (i.e., the coalescence rate) over small cells (cf. `Scalable_coalescent_prior`,
  // which averages N instead) because calculating exact integrals for a time-varying population over thousands
  // of events can get expensive, especially for Skygrid models.
  
  // Possible events as we move forward through time
  struct Add_to_bundle {
    int b;
  };
  struct Remove_from_bundle {
    int b;
  };
  struct Cross_probe {};
  struct Coal_rate_change {
    double new_coal_rate;
  };
  using Event_data = std::variant<Add_to_bundle, Remove_from_bundle, Cross_probe, Coal_rate_change>;
  struct Event {
    double t;
    Event_data data;
  };
  auto events = std::vector<Event>{};

  // Add all the different kinds of fixed events
  for (auto i = 0; i != std::ssize(probe_times); ++i) {
    events.push_back({probe_times[i], Cross_probe{}});
  }
  for (const auto& e : bundle_events) {
    if (e.adding) {
      events.push_back({e.t, Add_to_bundle{e.target_bundle}});
    } else {
      events.push_back({e.t, Remove_from_bundle{e.target_bundle}});
    }
  }
  
  // No point in going beyond last probe
  auto last_probe_time = probe_times.back();
  std::erase_if(events, [last_probe_time](const auto& e) { return e.t > last_probe_time; });
  
  // Approximate the population curve by lots of small flat segments (cf. `Scalable_coalescent_prior`, see below)
  auto min_t = +std::numeric_limits<double>::max();
  auto max_t = -std::numeric_limits<double>::max();
  auto has_finite_t_events = false;
  for (const auto& e : events) {
    if (not std::isinf(e.t)) {
      min_t = std::min(min_t, e.t);
      max_t = std::max(max_t, e.t);
      has_finite_t_events = true;
    }
  }
  CHECK(has_finite_t_events);  // Guaranteed by the has_finite_probe_time conditional above (throws)
  auto coal_rate = 42.0;  // == 1 / N_e(t) (units of 1/time);  value of 42.0 replaced below

  if (min_t == max_t) {
    // Fallback for very degenerate inputs:
    coal_rate = 1.0 / pop_model.pop_at_time(min_t);
  } else {
    // Cell boundaries: a uniform grid over [min_t, max_t], refined by all the probe times.
    // On its own, the uniform grid can be very coarse around the probes if the tree extends far
    // into the past (e.g., root at t = -1000 and probes in [0, 10] => only ~4 cells cover the probes).
    // Adding the probe times as boundaries guarantees that the population curve is resolved at
    // least as finely as the probes themselves.
    auto tot_cells = 400;  // hard-coded heuristic!
    auto cell_dt = (max_t - min_t) / tot_cells;
    auto cell_boundaries = std::vector<double>{};
    for (auto cell = 0; cell != tot_cells; ++cell) {
      cell_boundaries.push_back(max_t - cell * cell_dt);
    }
    cell_boundaries.push_back(min_t);
    for (const auto& t_i : probe_times) {
      if (min_t < t_i && t_i < max_t) {  // Also excludes infinite probe times
        cell_boundaries.push_back(t_i);
      }
    }
    std::ranges::sort(cell_boundaries);
    auto [first_dup, last_dup] = std::ranges::unique(cell_boundaries);
    cell_boundaries.erase(first_dup, last_dup);

    for (auto cell = 0; cell + 1 < std::ssize(cell_boundaries); ++cell) {
      auto t_min = cell_boundaries[cell];
      auto t_max = cell_boundaries[cell + 1];

      // Here we deviate somewhat from Scalable_coalescent_prior in using mean(1/N) instead of 1/mean(N).
      // Since this calculation does not interact with the MCMC in any way, we can afford a different (and slightly
      // better!) approximation here.  The results of these tree probing calls are user-facing analyses, not
      // components of a likelihood or acceptance ratio.
      auto cell_coal_rate = pop_model.intensity_integral(t_min, t_max) / (t_max - t_min);
      CHECK(not std::isnan(cell_coal_rate));

      events.push_back({t_min, Coal_rate_change{cell_coal_rate}});
      if (cell == 0) {
        coal_rate = cell_coal_rate;  // Earliest value of coal_rate =~ 1/N, used for events before min_t (at -infty)
      }
    }
  }

  // Sort events in increasing order of time (on ties, sort Adds before Removes)
  std::sort(events.begin(), events.end(), [](const Event& a, const Event& b) {
    // Earlier time => earlier event
    if (a.t < b.t) { return true; }
    if (a.t > b.t) { return false; }
    
    // Lower variant index = earlier event => Adds before Removes
    // (incidentally, equal-time probes and coalescence rate changes get processed after Removes)
    static_assert(std::is_same_v<std::variant_alternative_t<0, Event_data>, Add_to_bundle>);
    static_assert(std::is_same_v<std::variant_alternative_t<1, Event_data>, Remove_from_bundle>);
    return a.data.index() < b.data.index();
  });

  // Machinery to keep track of active bundles and their attachment probability in the current inter-probe interval
  struct Active_bundle_info {
    int b;
    int k_b;  // number of active branches in bundle
    double p_attach_b; // prob. of a probe at current time t attaching to this bundle _in this inter-probe interval_
  };
  auto p_no_attach = 1.0;
  auto active_bundles = std::vector<Active_bundle_info>{};
  auto bundle_2_active_bundle = std::vector<int>(num_bundles, -1);  // bundle_2_active_bundle[b] == -1 => b is not active
  auto num_active_branches = 0;

  // Main event loop
  auto t = events.front().t;
  auto cur_probe_index = 0;
  for (const auto& next_event : events) {

    // Advance time to next event
    if (next_event.t > t) {
      CHECK_GE(num_active_branches, 0);
      if (num_active_branches > 0) {
        // The probe at t_2 could attach at any point between t and next_event.t
        auto k = static_cast<double>(num_active_branches);
        auto one_over_k = 1.0 / k;
        auto dt = next_event.t - t;
        auto p_attach_dt = -std::expm1(-k * coal_rate * dt);  // 1 - exp(-k * coal_rate * dt)
        p_no_attach -= p_attach_dt * p_no_attach;             // p_no_attach = (1 - p_attach_dt) * p_no_attach
        for (auto& active_bundle : active_bundles) {
          auto k_b = active_bundle.k_b;
          CHECK_GE(k_b, 0);
          active_bundle.p_attach_b += p_attach_dt * (k_b * one_over_k - active_bundle.p_attach_b);
          //  Alt: p_attach_b = p_attach_dt * (k_b / k) + (1-p_attach_dt) * p_attach_b
        }
      }
    }
    
    // Handle event!
    t = next_event.t;
    std::visit(estd::overloaded{
        
        [&](const Add_to_bundle& e) {
          auto active_bundle = bundle_2_active_bundle[e.b];
          if (active_bundle == -1) {
            active_bundle = std::ssize(active_bundles);
            bundle_2_active_bundle[e.b] = active_bundle;
            active_bundles.push_back(Active_bundle_info{ .b = e.b, .k_b = 0, .p_attach_b = 0.0 });
          }
          ++active_bundles[active_bundle].k_b;
          ++num_active_branches;
        },
        [&](const Remove_from_bundle& e) {
          auto active_bundle = bundle_2_active_bundle[e.b];
          CHECK_GE(active_bundle, 0);
          CHECK_LT(active_bundle, std::ssize(active_bundles));
          CHECK_GT(active_bundles[active_bundle].k_b, 0);  // Not even transiently negative, equal-time Adds come before Removes
          --active_bundles[active_bundle].k_b;
          --num_active_branches;

          // If k_b == 0, keep it around so p_attach_b keeps being scaled down until the next probe
        },
        
        [&](const Cross_probe&) {
          // Accumulate results for all bundles
          for (auto b = 0; b != num_bundles; ++b) {
            auto p = cur_probe_index == 0 ? 0.0 : out_results[b * num_probe_times + cur_probe_index-1];
            out_results[b * num_probe_times + cur_probe_index] = p_no_attach * p;
          }
          for (const auto& active_bundle : active_bundles) {
            out_results[active_bundle.b * num_probe_times + cur_probe_index] += active_bundle.p_attach_b;
          }

          // Deactivate bundles with no active branches
          for (const auto& active_bundle : active_bundles) {
            bundle_2_active_bundle[active_bundle.b] = -1;  // will reset below
          }
          std::erase_if(active_bundles, [](const auto& active_bundle) { return active_bundle.k_b == 0; });
          for (auto ab = 0; ab != std::ssize(active_bundles); ++ab) {
            bundle_2_active_bundle[active_bundles[ab].b] = ab;
          }

          // Reset all attachment probabilities for next inter-probe interval
          p_no_attach = 1.0;
          for (auto& active_bundle : active_bundles) {
            active_bundle.p_attach_b = 0.0;
          }
          ++cur_probe_index;
        },
            
        [&](const Coal_rate_change& e) {
          coal_rate = e.new_coal_rate;
        },
    }, next_event.data);
  }

  // Check we crossed all the probes
  CHECK_EQ(cur_probe_index, num_probe_times);
}

auto make_uniform_probe_times(
    double t_start,
    double t_end,
    int num_t_cells)
    -> std::vector<double> {

  if (not (t_start < t_end)) {
    throw std::invalid_argument(absl::StrFormat(
        "Invalid probe times: need t_start < t_end, but t_start=%g and t_end=%g", t_start, t_end));
  }
  if (num_t_cells <= 0) {
    throw std::invalid_argument(absl::StrFormat("Number of probe cells should be positive, not %d", num_t_cells));
  }

  auto dt = (t_end - t_start) / num_t_cells;
  
  auto result = std::vector<double>{};
  result.reserve(num_t_cells);
  for (auto cell = 0; cell != num_t_cells; ++cell) {
    result.push_back(t_start + (cell+1) * dt);
  }

  return result;
}

}  // namespace delphy
