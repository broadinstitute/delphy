#ifndef DELPHY_COALESCENT_H_
#define DELPHY_COALESCENT_H_

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ranked_index.hpp>

#include "pop_model.h"

namespace delphy {

// A coalescent prior for a tree with N tips and N-1 inner nodes
// Initially, all node times are set to 0.0.  Use displace_tip and displace_coalescence to change them
class Coalescent_prior {
 public:
  Coalescent_prior(const Pop_model* pop_model, int num_tips);

  auto displace_tip(double old_t, double new_t) -> void;
  auto displace_coalescence(double old_t, double new_t) -> void;

  auto count_tips_between(double t_min_closed, double t_max_open) const -> int;
  auto count_coalescences_between(double t_min_closed, double t_max_open) const -> int;

  auto calc_log_prior() const -> double;
  auto calc_delta_log_prior_after_displace_coalescence(double old_t, double new_t) const -> double;

 private:
  const Pop_model* pop_model_;

  // TimeContainer holds times in a ranked index, so we can quickly count items in a certain time range
  using TimeContainer = ::boost::multi_index::multi_index_container<
      double,
      ::boost::multi_index::indexed_by<
          ::boost::multi_index::ranked_non_unique<::boost::multi_index::identity<double>>
      >
  >;

  TimeContainer tip_times_, coalescence_times_;
};

}  // namespace delphy

#endif // DELPHY_COALESCENT_H_
