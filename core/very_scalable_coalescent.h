#ifndef DELPHY_VERY_SCALABLE_COALESCENT_H_
#define DELPHY_VERY_SCALABLE_COALESCENT_H_

#include "pop_model.h"
#include "phylo_tree.h"

namespace delphy {

namespace Very_scalable_coalescent_prior {
auto cell_for(double t, double t_ref, double t_step) -> int;
auto cell_ubound(int cell, double t_ref, double t_step) -> double;
auto cell_lbound(int cell, double t_ref, double t_step) -> double;
auto add_interval(
    double t_start,
    double t_end,
    double delta_k,
    std::vector<double>& k,
    double t_ref,
    double t_step)
    -> void;
};

class Very_scalable_coalescent_prior_part {
 public:
  Very_scalable_coalescent_prior_part(
      std::shared_ptr<const Pop_model> pop_model,
      const Phylo_tree& subtree,
      std::mt19937& prng,
      bool includes_tree_root,
      double t_ref,
      double t_step,
      std::vector<double> k_bar_p,
      std::vector<double> k_twiddle_bar_p,
      std::vector<double> k_twiddle_bar,
      std::vector<double> popsize_bar,
      std::vector<int> num_active_parts);

  auto coalescence_displaced(double old_t, double new_t) -> void;
  auto tip_displaced(double old_t, double new_t) -> void;
  auto calc_partial_log_prior() const -> double;
  auto calc_delta_partial_log_prior_after_displace_coalescence(double old_t, double new_t) -> double;
  auto calc_delta_partial_log_prior_after_displace_tip(double old_t, double new_t) -> double;

 private:
  std::shared_ptr<const Pop_model> pop_model_;
  const Phylo_tree* subtree_;
  std::mt19937* prng_;
  bool includes_tree_root_;
  std::vector<double> k_bar_p_;
  std::vector<double> k_twiddle_bar_p_;
  std::vector<double> k_twiddle_bar_;
  std::vector<double> popsize_bar_;
  std::vector<int> num_active_parts_;
  double t_ref_;
  double t_step_;

  void ensure_space(double t);
  auto calc_delta_partial_log_prior_on_add_interval(double min_t, double max_t, double delta_k) -> double;
};

auto make_very_scalable_coalescent_prior_parts(
    const std::vector<const Phylo_tree*>& subtrees,
    int root_partition_index,
    std::shared_ptr<const Pop_model> pop_model,
    std::vector<std::mt19937>& prngs,
    double t_step)
    -> std::vector<Very_scalable_coalescent_prior_part>;

}  // namespace delphy

#endif // DELPHY_VERY_SCALABLE_COALESCENT_H_
