#ifndef DELPHY_GENERIC_TREE_PROBER_H_
#define DELPHY_GENERIC_TREE_PROBER_H_

#include <span>

#include "pop_model.h"

namespace delphy {

struct Bundle_event {
  double t;
  bool adding;         // true => adding one branch to `target_bundle`; false => removing one branch
  int target_bundle;
};

// Calculate for various times t_i the probability that a probe sample at time t_i is
// descended directly from each of several branch bundles in a tree.  Bundles are indexed
// 0, 1, ..., B-1.  The tree structure is provided implicitly as a series of events where
// a single branch is added or removed from a single bundle as time advances.  Events
// should correspond to a coherent history, i.e., at no time should the net number of
// branches in each bundle drop below 0.  If needed, the branch above the root can be
// encoded by an event at t=-infty; lacking such an event, there's a finite probability of
// a probe never attaching to the tree (so the sum of the attachment probabilities over
// all bundles at a fixed time will be less than 1.0).
//
// All manner of specialized tree probers can be built on top of this interface.
//
// To avoid needless repeated copying, the result is written into a preallocated 2D array
// `out_results`, accessed via a 1D span using row-major ordering.  You'd want to use
// `out_results[b][i]` to refer to the results for bundle `b` and probe time `i`; instead
// use `out_results[b * std::ssize(probe_times) + i]`.  One day, when we move to C++23, we
// can use a `std::mdspan` directly instead.  NOTE: every element of `out_results` is
// overwritten, so no need to initialize it before the call.
auto generic_probe_tree(
    std::span<const Bundle_event> bundle_events,  // in any time order, but should be coherent (see description)
    int num_bundles,
    const Pop_model& pop_model,
    std::span<const double> probe_times,          // sorted past-to-future
    std::span<double> out_results)                // pop. fraction at probe_times[i] descended directly from branch bundle b
    -> void;

// Helper for generic probe times as an adapter to older rigid interfaces
auto make_uniform_probe_times(
    double t_start,
    double t_end,
    int num_t_cells)
    -> std::vector<double>;

}  // namespace delphy

#endif // DELPHY_GENERIC_TREE_PROBER_H_
