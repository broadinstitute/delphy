#ifndef DELPHY_SPR_STUDY_H_
#define DELPHY_SPR_STUDY_H_

#include "phylo_tree.h"
#include "site_deltas.h"

namespace delphy {

// A candidate region is a portion of a single branch where a subtree rooted at X might be reattached.  The portion
// ranges from [[b,muts[im-1].t]] (exclusive) to [[b,muts[im].t]] (inclusive), where b is the branch index `branch` and
// `im` is the index `mut_idx` of the mutation that delimits the end of this candidate region.  For convenience,
// muts[-1].t is the branch start time and muts[muts.size()].t is the branch end time.  Reattaching X anywhere on this
// region requires a minimum of min_muts mutations.
//
// Note: we cannot attach X at any time on or after t_X.
//
struct Candidate_region {
  // Inputs
  Branch_index branch;
  int mut_idx;
  double t_min;
  double t_max;
  int min_muts;

  // Outputs
  double log_W_over_Wmax = 0.0;
  double W_over_Wmax = 0.0;

  auto is_above_root() const -> bool { return t_min == -std::numeric_limits<double>::max(); }
  
  auto operator<=>(const Candidate_region& that) const = default;
};
inline auto operator<<(std::ostream& os, const Candidate_region& region) -> std::ostream& {
  return os << absl::StreamFormat("Candidate_region{branch=%d, mut_idx=%d, %g<t<=%g, min_muts=%d, W_over_Wmax=%g}",
                                  region.branch, region.mut_idx,
                                  region.t_min, region.t_max, region.min_muts, region.W_over_Wmax);
}

// An Spr_study_builder represents an ongoing restricted depth-first traversal of a tree.  At any moment, the builder is
// at a certain region of the tree: on branch `cur_branch`, between mutation `cur_mut_idx-1` and `cur_mut_idx`; where
// the rogue indices -1 and std::ssize(muts) denote the beginning and end of the branch, respectively.  The builder
// knows the deltas between the sequence at any point in that region and X.  It has a stack of movements to perform
// (advance to the next/previous region on the same branch, or the first/last region of a child/parent branch, with some
// idea of whether we're traversing forward or we are backtracking).  The caller is responsible for initial priming of
// movements.  After executing a forward movement, more work items may be added to the work stack unless a stopping
// criterion has been reached (e.g., too many deltas, we're at a point later than t_X).
//
// The builder will not enter the branch X if it is present, and will generate regions as if the parent node P of X were
// not in the tree if it's there (i.e., the regions on branch P will be reported as on branch S, and the final region of
// P and initial region of S will be merged).  Further, if P is the root of the tree and `include_above_root` is set to
// `true`, then the region (S, std::ssize(S.mutations)) will be extended to -Infinity, removing all the preceding
// regions (as if those paths had disappeared when X was removed).  If instead `include_above_root` is `false`, then the
// region above P will be removed.  Before these adjustments, the regions describe the tree while X is still attached.
// After these adjustments, the regions reflect the status of the tree after X is detached (`includes_run_root`
// controls, among other things, whether such a detachment would collapse the root to X's sibling if X was a child of
// the root).  These adjustments are important for implementing reversible SPR moves when the subtree rooted at X is not
// physically detached from the tree.  These changes are implemented as a post-processing step.
//
// The traversal will be cut short after encountering more than `max_muts_from_start` mutations from the starting
// point.  By default, this limit is disabled, but it is useful for implementing very local SPR moves.
//
// The builder accumulates a set of `Candidate_region`s as it progreses through its work stack.  The caller is free
// to interleave (a) moving the builder, (b) adding new items to the work stack, and (c) executing the work implied
// by the work stack, so long as internal consistency is maintained.  This may make the builder useful in tricky situations.
//
// We expect most users to simply use the convenience function `seed_fill_from` (see spr_study_tests.cpp and usages in
// phylo_tree.cpp for examples).
//
struct Spr_study_builder {
  // The builder is at the point just downstream of [[cur_branch, muts[cur_mut_idx].t]].
  // If cur_mut_idx == -1, then the builder is at the beginning of cur_branch.
  // If cur_mut_idx == muts.size(), then the builder is at the end of cur_branch
  const Phylo_tree* tree;
  Branch_index cur_branch;
  int cur_mut_idx;
  int cur_muts_from_start;

  Site_deltas cur_to_X_deltas;
  const Scratch_interval_set* missing_at_X;
  
  // Stopping criteria
  Node_index X = k_no_node;
  double t_X = std::numeric_limits<double>::max();
  int max_muts_from_start = std::numeric_limits<int>::max();

  auto cur_region_in_scope() const -> bool {
    return cur_branch != X &&
        cur_muts_from_start <= max_muts_from_start;
  }
  auto region_t_min(Branch_index branch, int mut_idx) const -> double {
    if (branch == tree->root) { return -std::numeric_limits<double>::max(); }
    const auto& muts = tree->at(branch).mutations;
    if (mut_idx == 0) { return tree->branch_begin_t(branch); }
    else { return muts.at(mut_idx - 1).t; }
  }
  auto region_t_max(Branch_index branch, int mut_idx) const -> double {
    if (branch == tree->root) { return tree->at(branch).t; }
    const auto& muts = tree->at(branch).mutations;
    if (mut_idx == std::ssize(muts)) { return tree->branch_end_t(branch); }
    else { return muts.at(mut_idx).t; }
  }

  // Work items
  struct Work_item {
    Branch_index target_branch;
    int target_mut_idx;
    bool is_backtracking;
  };
  Scratch_vector<Work_item> work_stack;

  // Results
  Scratch_vector<Candidate_region> result;
  
  auto visit_cur_region() -> void;
  auto seed_neighbors_except(Branch_index old_branch, int old_mut_idx) -> void;
  auto account_for_Xs_detachment(bool can_change_root) -> void;
  auto remove_regions_in_Xs_future() -> void;
  auto move_to_neighbor(
      Branch_index target_branch,
      int target_mut_idx,
      bool is_backtracking)
      -> void;

  Spr_study_builder(
      const Phylo_tree& tree,
      Node_index X,
      double t_X,
      const Scratch_interval_set& missing_at_X)
      : tree{&tree},
        cur_branch{k_no_node},
        cur_mut_idx{-1},
        cur_muts_from_start{0},
        cur_to_X_deltas{},
        missing_at_X{&missing_at_X},
        X{X},
        t_X{t_X},
        work_stack{},
        result{} {}

  auto is_next_to_cur_region(Branch_index target_branch, int target_mut_idx) -> bool {
    if (cur_branch == k_no_node || target_branch == k_no_node) {
      return true;
    }
    if (cur_branch == target_branch) {
      return abs(target_mut_idx - cur_mut_idx) == 1;
    }
    if (target_branch == tree->at(cur_branch).parent) {
      return cur_mut_idx == 0 && target_mut_idx == std::ssize(tree->at(target_branch).mutations);
    }
    if (cur_branch == tree->at(target_branch).parent) {
      return cur_mut_idx == std::ssize(tree->at(cur_branch).mutations) && target_mut_idx == 0;
    }
    return false;
  }
  
  auto add_forward_movement(Branch_index target_branch, int target_mut_idx) -> void {
    // Push backtrack into stack first
    CHECK(is_next_to_cur_region(target_branch, target_mut_idx));
    work_stack.push_back(Work_item{cur_branch, cur_mut_idx, true});
    work_stack.push_back(Work_item{target_branch, target_mut_idx, false});
  }

  auto do_pending_work() -> void;

  auto seed_fill_from(
      Branch_index cur_branch,
      int cur_mut_idx,
      Site_deltas cur_to_X_deltas,
      bool can_change_root)
      -> void;
};

// An Spr_study is a collection of candidate regions, from which we will sample one new trial nexus
struct Spr_study {
  Spr_study(
      Spr_study_builder&& builder,  // candidate regions moved to study, hence require explicit move context
      double lambda_X,
      double annealing_factor,
      double t_X);
  
  // Inputs
  const Phylo_tree* tree;
  double lambda_X;
  double mu;
  double annealing_factor;  // Should be <= 1
  double t_X;
  Scratch_vector<Candidate_region> candidate_regions;

  // Outputs
  double log_Wmax;
  double sum_W_over_Wmax;
  int all_min_muts;
  double total_weight;

  // Calculations
  auto pick_nexus_region(absl::BitGenRef bitgen) const -> int;
  auto pick_time_in_region(int region_idx, absl::BitGenRef bitgen) const -> double;
  auto find_region(Branch_index branch, double t) const -> int;
  auto log_alpha_in_region(int region_idx, double t) const -> double;

  // Debugging
  auto dump() -> void;
};


// struct Spr_study_result {
//   int chosen_candidate_idx;
//   double t;

//   auto phylo_tree_loc(const Spr_study& study) const -> Phylo_tree_loc {
//     return {study.candidates[chosen_candidate_idx].branch, t};
//   }
// };

// auto calc_weights(Spr_study& study) -> void;
// auto pick_insertion_point(const Spr_study& study, absl::BitGenRef bitgen) -> Spr_study_result;
// auto calc_log_pdf_at_insertion_point(const Spr_study& study, const Spr_study_result& result) -> double;
// auto find_candidate_i(const Spr_study& study, const Phylo_tree_loc& query) -> int;

}  // namespace delphy

#endif // DELPHY_SPR_STUDY_H_
