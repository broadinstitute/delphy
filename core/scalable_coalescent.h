#ifndef DELPHY_SCALABLE_COALESCENT_H_
#define DELPHY_SCALABLE_COALESCENT_H_

#include <deque>

#include "pop_model.h"
#include "tree.h"

namespace delphy {

// A coalescent prior for a tree with N tips and N-1 inner nodes
// Initially, all node times are set to 0.0.  Use displace_tip and displace_coalescence to change them
class Scalable_coalescent_prior {
 public:
  Scalable_coalescent_prior(std::shared_ptr<const Pop_model> pop_model, Node_index num_nodes, double t_ref, double t_step);

  auto t_step() const -> double { return t_step_; }
  auto k_bars() const -> const std::deque<double>& { return k_bars_; }
  auto popsize_bars() const -> const std::deque<double>& { return popsize_bars_; }

  auto cell_for(double t) const -> int;
  auto cell_lbound(int cell) const -> double { return t_ref_ + ((cell - cells_before_t_ref_) * t_step_); }
  auto cell_ubound(int cell) const -> double { return cell_lbound(cell) + t_step_; }
  
  auto reset(double t_step) -> void;

  auto mark_as_tip(Node_index node) -> void { node_infos_.at(node).is_tip = true; }
  auto mark_as_coalescence(Node_index node) -> void { node_infos_.at(node).is_tip = false; }

  auto displace_tip(Node_index node, double new_t) -> void;
  auto displace_coalescence(Node_index node, double new_t) -> void;
  auto pop_model_changed(std::shared_ptr<const Pop_model> new_pop_model) -> void;

  auto calc_log_prior() const -> double;
  auto calc_delta_log_prior_after_displace_coalescence(int i, double new_t) -> double;

 private:
  std::shared_ptr<const Pop_model> pop_model_;
  struct Node_info {
    double t;
    bool is_tip;
  };
  Node_vector<Node_info> node_infos_;
  std::deque<double> k_bars_;
  std::deque<double> popsize_bars_;
  double t_ref_;
  double t_step_;
  int cells_before_t_ref_;

  auto assert_space(double t) -> void;
  auto ensure_space(double t) -> void;
  auto add_interval(double t_start, double t_end, double delta_k) -> void;
};

}  // namespace delphy

#endif // DELPHY_SCALABLE_COALESCENT_H_
