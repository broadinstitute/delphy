#include "spr_study.h"

#include <numbers>

#include <boost/math/special_functions/gamma.hpp>

namespace delphy {

auto Spr_study_builder::seed_fill_from(
    Branch_index init_branch,
    int init_mut_idx,
    Site_deltas init_to_X_deltas,
    bool can_change_root)
    -> void {
  CHECK(work_stack.empty());
  CHECK_EQ(cur_branch, k_no_node);
  
  cur_to_X_deltas = std::move(init_to_X_deltas);
  add_forward_movement(init_branch, init_mut_idx);
  do_pending_work();

  account_for_Xs_detachment(can_change_root);
  remove_regions_in_Xs_future();
}

auto Spr_study_builder::do_pending_work() -> void {
  while (not work_stack.empty()) {
    auto [target_branch, target_mut_idx, is_backtracking] = work_stack.back();
    work_stack.pop_back();
    
    auto old_branch = cur_branch;
    auto old_mut_idx = cur_mut_idx;
    
    move_to_neighbor(target_branch, target_mut_idx, is_backtracking);
    
    if (not is_backtracking && cur_region_in_scope()) {
      visit_cur_region();
      seed_neighbors_except(old_branch, old_mut_idx);
    }
  }
}

auto Spr_study_builder::move_to_neighbor(
    Branch_index target_branch,
    int target_mut_idx,
    bool is_backtracking)
    -> void {
  
  CHECK(is_next_to_cur_region(target_branch, target_mut_idx));

  if (cur_branch != k_no_node) {
    
    // If we're moving across a mutation, update internal state
    if (target_branch == cur_branch) {
        
      const auto& muts = tree->at(cur_branch).mutations;
      CHECK_GE(cur_mut_idx, 0);
      CHECK_LE(cur_mut_idx, std::ssize(muts));
      
      CHECK_GE(target_mut_idx, 0);
      CHECK_LE(target_mut_idx, std::ssize(muts));
        
      if (target_mut_idx == cur_mut_idx + 1) {         // Moving down past cur_mut_idx
        const auto& m = muts.at(cur_mut_idx);
        if (not missing_at_X->contains(m.site)) {
          DCHECK(not cur_to_X_deltas.contains(m.site) || cur_to_X_deltas.at(m.site).from == m.from);
          pop_front_site_deltas(m, cur_to_X_deltas);
          DCHECK(not cur_to_X_deltas.contains(m.site) || cur_to_X_deltas.at(m.site).from == m.to);
          cur_muts_from_start += (not is_backtracking ? +1 : -1);
        }
          
      } else if (target_mut_idx == cur_mut_idx - 1) {  // Moving up past target_mut_idx
        const auto& m = muts.at(target_mut_idx);
        if (not missing_at_X->contains(m.site)) {
          DCHECK(not cur_to_X_deltas.contains(m.site) || cur_to_X_deltas.at(m.site).from == m.to);
          push_front_site_deltas(m, cur_to_X_deltas);
          DCHECK(not cur_to_X_deltas.contains(m.site) || cur_to_X_deltas.at(m.site).from == m.from);
          cur_muts_from_start += (not is_backtracking ? +1 : -1);
        }
          
      } else {
        CHECK(false) << absl::StreamFormat("Inconsistent work item on branch %d: mut %d to mut %d",
                                           cur_branch, cur_mut_idx, target_mut_idx);
      }
    }
  }

  // We're here
  cur_branch = target_branch;
  cur_mut_idx = target_mut_idx;
}

auto Spr_study_builder::visit_cur_region() -> void {
  result.push_back(Candidate_region{
      .branch = cur_branch,
      .mut_idx = cur_mut_idx,
      .t_min = region_t_min(cur_branch, cur_mut_idx),
      .t_max = region_t_max(cur_branch, cur_mut_idx),
      .min_muts = static_cast<int>(std::ssize(cur_to_X_deltas))
    });
}

auto Spr_study_builder::seed_neighbors_except(Branch_index old_branch, int old_mut_idx) -> void {
  // The current region has up to three neighbours.  Explore all except where we just came from
  auto maybe_move_to = [&](int next_branch, int next_mut_idx) {
    if (next_branch == old_branch && next_mut_idx == old_mut_idx) { return; }
    add_forward_movement(next_branch, next_mut_idx);
  };

  // Up to 1 parent region
  if (cur_branch != tree->root) {
    if (cur_mut_idx > 0) {
      maybe_move_to(cur_branch, cur_mut_idx - 1);
    } else {
      auto parent_branch = tree->at(cur_branch).parent;
      maybe_move_to(parent_branch, std::ssize(tree->at(parent_branch).mutations));
    }
  }

  // Up to 2 child regions
  if (cur_mut_idx < std::ssize(tree->at(cur_branch).mutations)) {
    maybe_move_to(cur_branch, cur_mut_idx + 1);
  } else {
    for (const auto& child : tree->at(cur_branch).children) {
      maybe_move_to(child, 0);
    }
  }
}

auto Spr_study_builder::account_for_Xs_detachment(bool can_change_root) -> void {
  if (X == k_no_node) {
    // X is already detached, no need to pretend anything
    // But we may still need to remove any region above the root
    if (not can_change_root) {
      std::erase_if(result, [r = tree->root](const auto& region) { return region.branch == r; });
    }
    return;
  }
  
  // The regions around P are a bit tricky, since we're supposed to pretend that
  // X has been pruned and P is not there.  And if P is the root and we can change the root,
  // we have to pretend as if the whole P-S branch disappears
  CHECK_NE(X, tree->root);
  auto P = tree->at(X).parent;
  auto S = tree->at(P).sibling_of(X);
  auto num_muts_G_to_P = std::ssize(tree->at(P).mutations);

  for (auto& region : result) {
    // May need to remove any region above the root is we can't change the root
    if (region.branch == tree->root && not can_change_root) {
      region.branch = -1;  // Mark for deletion below;
      continue;
    }

    // Otherwise, no regions except those on the G-P and P-S branches are adjusted
    if (region.branch != S && region.branch != P) { continue; }

    if (P != tree->root) {
      if (region.branch == S) {
        if (region.mut_idx == 0) {
          // The first region in S is merged with the last region in P
          region.t_min = region_t_min(P, num_muts_G_to_P);
        }
        
        // All regions have an offset in mut_idx as the mutations in G-P become part of the G-S "branch"
        region.mut_idx += num_muts_G_to_P;
        
      } else if (region.branch == P) {
        if (region.mut_idx == num_muts_G_to_P) {
          // The first region in S is merged with the last region in P
          region.branch = -1;  // Mark for deletion below
        } else {
          // All other regions are part of the G-S "branch"
          region.branch = S;
        }
      }
      
    } else {
      // Much trickier! P is the root

      if (not can_change_root) {
        // Phew!  Actually, just the regions above P disappear
        if (region.branch == P) {
          region.branch = -1;  // Mark for deletion below
        }
      } else {
        // The last region in the P-S branch becomes the "above root" region; every other region disappears
        if (region.branch == S && region.mut_idx == std::ssize(tree->at(S).mutations)) {
          // We get this index slightly wrong if the mutations above the root are partly cancelled out by the P-S mutations,
          // and/or if the P-S mutations include reversals.  But by this point, this index is purely informative, so we
          // can tolerate the error
          region.mut_idx += num_muts_G_to_P;
          region.t_min = -std::numeric_limits<double>::max();
        } else {
          // Every other region on the G-P or P-S branches disappears
          region.branch = -1;  // Mark for deletion below
        }
      }
    }
  }

  // Delete regions marked above
  std::erase_if(result, [](const auto& region) { return region.branch == -1; });
}

auto Spr_study_builder::remove_regions_in_Xs_future() -> void {
  // We need to fully explore the possible regions before truncating on t_X in case
  // account_for_Xs_detachment() above mucks about with region merging / boundaries
  for (auto& region : result) {
    if (region.t_min >= t_X) {
      region.branch = -1;  // Mark for deletion below
    } else if (region.t_max > t_X) {
      region.t_max = t_X;
    }
  }

  // Delete regions marked above
  std::erase_if(result, [](const auto& region) { return region.branch == -1; });
}

Spr_study::Spr_study(
    Spr_study_builder&& builder,  // candidate regions moved to study, hence require explicit move context
    double lambda_X,
    double annealing_factor,
    double t_X)
    : tree{builder.tree},
      lambda_X{lambda_X},
      annealing_factor{annealing_factor},
      t_X{t_X},
      candidate_regions{std::move(builder.result)} {

  mu = lambda_X / (tree->num_sites() - builder.missing_at_X->num_sites());
  
  // While we *can* calculate the weights for each region exactly in the case of a Jukes-Cantor model
  // (see commented code below), it's often just as effective to use a crude approximation that still
  // favors low-mut regions by roughly the right amount
  
  // We approximate the pdf of a point t in a region R as:
  //
  //   p(t) =~ f*lambda_X * { exp[-lambda_X (t_X - t')] * [mu (t_X - t') / 3]^m }^f,
  //
  // so the weight integrated over the whole region is:
  //
  //   W_R =~ f*lambda_X*(t_max - t_min) * { exp[-lambda_X (t_X - t')] * [mu (t_X - t') / 3]^m }^f,
  //
  //   => log(W_R) =~ log(f*lambda_X*(t_max - t_min)) + f * {-lambda_X (t_X - t') + m * log(mu (t_X - t') / 3)}
  //
  // Here, the region spans the time t_min <= t <= t_max, and t' := (t_min + t_max)/2.
  //
  // Since the factor mu*(t_X - t')/3 is typically tiny, the weights will heavily skew
  // towards regions with the minimum number of mutations.
  //
  // The exponent `f` is an annealing factor that should be at most 1.  When f = 1 and the real model is a
  // Jukes-Cantor model with no site rate heterogeneity and there's no missing data, then the acceptance probability of
  // an SPR move is almost exactly 1 in the MH criterion (i.e., the SPR study would choose the new nexus as if it could
  // Gibbs sample it).  By setting f < 1, we decrease the propensity of the SPR study to immediately go for
  // low-mutation regions, effectively leaving some margin of error in the MH criterion.  It is equivalent
  // to applying a "higher temperature" to the posterior, in the manner of simulated annealing; hence the name.
  //
  // Note: For the region above the root, we need to do something different!  There, we do a much more accurate job,
  // because we need to estimate with reasonable accuracy *where* to regraft (it's no longer enough to guess
  // "somewhere uniformly along the branch", because the branch above the root extends to -Infinity).  We use:
  //
  //   p_root(t) =~ f*lambda_X * {exp[-lambda_X (t_X - t + t_S - t)] * [mu (t_X - t + t_S - t) / 3]^m}^f,  with t <= min(t_X, t_S).
  //
  // Here, S is the branch endpoint to which we're attaching (i.e., the tree root if X were detached).
  //
  // We massage this a bit by setting s := (t_X - t + t_S - t), so |dt/ds| = 1/2, s >= s_min = |t_X - t_S|, and
  //
  //   p_root(s) = 0.5 * f*lambda_X * {exp[-lambda_X s] * [mu s/3]^m}^f
  //             = 0.5 * f*lambda_X * (mu s/3)^{fm} exp[-lambda_X f s].
  //
  // This is simply a truncated Gamma distribution with alpha = fm+1 and beta = lambda_X f.  We can express the
  // integrated weight for points above the root in terms of the normalized upper incomplete gamma function Q(a,z),
  // which thankfully is built into Boost (the function and its inverse):
  //
  //   Q(a,z) := (1/Gamma(a)) int_z^infty x^{a-1} exp(-x) dx.
  //
  // Briefly, Q(a,z) drops smoothly from nearly 1 to nearly 0 around z == a +/- a/2.  See:
  // - https://www.boost.org/doc/libs/1_85_0/libs/math/doc/html/math_toolkit/sf_gamma/igamma.html
  // - https://www.boost.org/doc/libs/1_85_0/libs/math/doc/html/math_toolkit/sf_gamma/igamma_inv.html
  //
  // With that, we can change variables to x = lambda_X f s, and define x_min := lambda_X f s_min, to obtain:
  //
  //   W_R = int_{s_min}^infty p_root(s) ds
  //       = int_{x_min}^infty 0.5 * (mu/(3 lambda_X f))^{fm} x^{fm} exp[-x] dx
  //       = 0.5 * (mu/(3 lambda_X f))^fm Gamma(fm+1) Q(fm+1, x_min)
  //
  //  => log(W_R) = -log(2) + fm log(mu (lambda_X^{-1} / f) / 3) + log Gamma(fm+1) + log Q(fm+1, x_min).
  //
  // Unless we are in exceptional scenarios (e.g., a root branch that is very stretched out with respect
  // to the number of mutations), Q(fm+1, x_min) =~ 1, and this weight is again dominated by the term
  // fm log(mu (lambda_X^{-1} / f) / 3), which plays the same role as fm log(mu (t_X - t') / 3) above,
  // but with (lambda_X^{-1} / f) setting the timescale.

  auto f = annealing_factor;
  
  for (auto& region : candidate_regions) {
    
    auto t_min = region.t_min;
    auto t_max = region.t_max;
    auto m = region.min_muts;

    if (not region.is_above_root()) {
      auto t_prime = 0.5*(t_min + t_max);
      CHECK_LE(t_prime, t_X);
      region.log_W_over_Wmax =
          std::log(f*lambda_X*(t_max - t_min)) +
          f * (- lambda_X * (t_X - t_prime) + m * std::log(mu * (t_X - t_prime) / 3));
      // ^^ will be -inf if t_X == t_prime and/or t_min == t_max
      
    } else {
      auto t_S = tree->at(region.branch).t;
      auto s_min = std::abs(t_X - t_S);
      auto x_min = lambda_X * f * s_min;
      region.log_W_over_Wmax =
          -std::numbers::ln2
          + f*m*std::log(mu / (3 * lambda_X * f))
          + std::lgamma(f*m + 1)
          + std::log(boost::math::gamma_q(f*m + 1, x_min));
    }
  }

  CHECK(not candidate_regions.empty());
  log_Wmax = candidate_regions[0].log_W_over_Wmax;
  for (const auto& region : candidate_regions) {
    log_Wmax = std::max(log_Wmax, region.log_W_over_Wmax);
  }
  
  sum_W_over_Wmax = 0.0;
  for (auto& region : candidate_regions) {
    region.log_W_over_Wmax -= log_Wmax;
    region.W_over_Wmax = std::exp(region.log_W_over_Wmax);
    sum_W_over_Wmax += region.W_over_Wmax;
  }
}

auto Spr_study::dump() -> void {
  // Sort by branch and mut_idx
  auto sorted_regions = std::vector{candidate_regions};
  std::ranges::sort(sorted_regions, [](const auto& a, const auto& b) {
    return (a.branch < b.branch) || ((a.branch == b.branch) &&
                                     (a.mut_idx < b.mut_idx)); });

  auto total_T = 0.0;
  for (const auto& region : sorted_regions) {
    auto T = region.t_max - region.t_min;
    total_T += T;
    std::cerr << absl::StreamFormat("At %d,%d (%g < t < %g = %g days), %d muts (%g time so far)\n",
                                    region.branch, region.mut_idx, region.t_min, region.t_max, T,
                                    region.min_muts, total_T);
  }
}

auto Spr_study::pick_nexus_region(absl::BitGenRef bitgen) const -> int {
  // Choose a candidate region in proportion to its weight.
  //
  // We usually only ever pick one insertion point for a given SPR study, so there's
  // no better way to randomly draw a candidate than to scan all the candidates.
  auto r = absl::Uniform<double>(bitgen, 0.0, sum_W_over_Wmax);
  auto chosen_region_idx = 0;
  for (auto i = 0; i != std::ssize(candidate_regions); ++i) {
    const auto& region = candidate_regions[i];
    if (region.W_over_Wmax >= r) {
      chosen_region_idx = i;
      break;
    } else {
      r -= region.W_over_Wmax;
    }
  }

  return chosen_region_idx;
}

auto Spr_study::pick_time_in_region(int region_idx, absl::BitGenRef bitgen) const -> double {
  // Choose a time t within that region, according to pdf above.
  // If the region is below the root, that's uniformly over the region's span.
  // Otherwise, it's a truncated Gamma distribution
  const auto& region = candidate_regions[region_idx];
  auto t_min = region.t_min;
  auto t_max = region.t_max;
  if (not region.is_above_root()) {
    return absl::Uniform(absl::IntervalOpenClosed, bitgen, t_min, t_max);
  } else {
    // s := t_X - t + t_S - t
    //
    // First sample from the following truncated Gamma with inverse transform sampling:
    //
    //    p_root(s) ~ s^{fm} exp[-lambda_X f s}, with s >= |t_X - t_S|.
    //
    // Then derive t = (t_X + t_S - s) / 2
    //
    auto m = region.min_muts;
    auto f = annealing_factor;
    auto t_S = tree->at(region.branch).t;
    auto s_min = std::abs(t_X - t_S);
    auto Q_min = 0.0;
    auto Q_max = boost::math::gamma_q(f*m + 1, lambda_X*f*s_min);
    auto rand_Q = absl::Uniform(absl::IntervalOpenOpen, bitgen, Q_min, Q_max);
    //std::cerr << absl::StreamFormat("%.6g < Q = %.6g < %.6g, m = %d, f = %.3g\n", Q_min, rand_Q, Q_max, m, f);
    auto rand_s = boost::math::gamma_q_inv(f*m + 1, rand_Q)/(lambda_X*f);
    auto rand_t = 0.5 * (t_X + t_S - rand_s);
    CHECK_GE(rand_t, region.t_min - 1e-6);
    CHECK_LE(rand_t, region.t_max + 1e-6);
    rand_t = std::max(region.t_min, std::min(region.t_max, rand_t));
    return rand_t;
  }
}

// Find the index of the candidate region containing the query point.  Returns -1 if none of them do
auto Spr_study::find_region(Branch_index branch, double t) const -> int {
  // For now, nothing better than linear scan.  Maybe there's a way of organizing the data structures
  // differently
  for (auto i = 0; i != std::ssize(candidate_regions); ++i) {
    const auto& region = candidate_regions[i];
    if (region.branch == branch && region.t_min < t && t <= region.t_max) {
      return i;
    }
  }
  return -1;
}

auto Spr_study::log_alpha_in_region(int region_idx, double t) const -> double {
  // This is relatively easy in terms of the relative weights calculated by calc_weights above,
  // and is independent of t for the simplified pdf above everywhere except above the root.
  // The normalization factor is only important when comparing choosing two insertion points
  // from *different* SPR studies

  const auto& region = candidate_regions[region_idx];
  auto log_p_region = region.log_W_over_Wmax - std::log(sum_W_over_Wmax);
  if (not region.is_above_root()) {
    // p(t | R) = 1/(t_max - t_min)
    return log_p_region - std::log(region.t_max - region.t_min);
  } else {
    // p(s) ~ {exp[- lambda_X s] * s^m}^f = s^{fm} exp[-lambda_X f s]
    // p(s | R) = s^{fm} exp[-lambda_X f s] / {(lambda_X f)^{-(fm + 1)} [1 - P(fm+1, lambda_X f s_min)]}
    //          = (lambda_X f) (lambda_X f s)^{fm} exp[-lambda_X f s] / [1 - P(fm+1, lambda_X f s_min)]}
    // p(t | R) = 2 p(s | R)
    //          = 2 (lambda_X f) (lambda_X f s)^{fm} exp[-lambda_X f s] / [1 - P(fm+1, lambda_X f s_min)]
    auto f = annealing_factor;
    auto m = region.min_muts;
    auto t_S = tree->at(region.branch).t;
    auto s_min = std::abs(t_X - t_S);
    auto s = t_X - t + t_S - t;
    CHECK_GE(s, s_min - 1e-6);
    return log_p_region +
        std::numbers::ln2 +
        std::log(lambda_X*f) +
        f*m * std::log(lambda_X*f*s) +
        -lambda_X*f*s +
        -std::log(boost::math::gamma_q(f*m+1, lambda_X*f*s_min));
  }
}

// 
// auto calc_weights(Spr_study& study) -> void {
//   // Assuming L_X > 0, the weight of each region R is proportional to:
//   //
//   //   W_R = int_{t_min}^{t_max} dt exp[-mu L_X (t_X - t)] * [mu (t_X - t) / 3]^m,
//   //
//   // i.e., the probability of reattaching at any point of the tree is proportional
//   // to the branch genetic prior of the optimal P'-X branch, assuming a JC evo model.
//   //
//   // The above integral can be calculated easily:
//   //
//   //   W_R = [1/(mu L_X)] * m! / (3 L_X)^m * [P(m+1, delta_max) - P(m+1, delta_min)],
//   // with
//   //   delta_min := mu L_X (t_X - t_max)
//   //   delta_max := mu L_X (t_X - t_min)
//   //
//   // where P(a,z) is the normalized lower incomplete gamma function of a and z
//   // (smoothly going from nearly 0 to nearly 1 around z == a +/- a/2).  Thankfully,
//   // this function and its inverse are implemented efficiently in Boost:
//   // - https://www.boost.org/doc/libs/1_84_0/libs/math/doc/html/math_toolkit/sf_gamma/igamma.html
//   // - https://www.boost.org/doc/libs/1_84_0/libs/math/doc/html/math_toolkit/sf_gamma/igamma_inv.html
//   //
//   // To keep things numerically stable, we actually calculate a relative weight by
//   // factoring out [1/(mu L_X)] * (m_min)! / (3 L_X)^(m_min), where m_min is the minimum m
//   // across the entire study, as follows:
//   //
//   //   w_R := W_R / {[1/(mu L_X)] * (m_min)! / (3 L_X)^(m_min)}
//   //        = (m * (m-1) * ... * (m - m_min + 1)) / (3 L_X)^(m - m_min)
//   //          * [P(m+1, delta_max) - P(m+1, delta_min)]
//   //
//   // Hence, the regions with the minimal number of mutations have a relative weight of O(1),
//   // while every additional mutation m' beyond the minimum scales down the relative weight
//   // by a factor of m / (3 L_X).  In genomic epi datasets, m_min / L_X << 1, so the candidate
//   // regions will be overwhelmingly skewed towards the parsimonious ones.

//   CHECK(not study.candidates.empty());
//   CHECK_GT(study.L_X, 0);
  
//   study.all_min_muts = std::ranges::min(study.candidates, {}, [](const auto& R){ return R.min_muts; }).min_muts;
//   auto t_X = study.t_X;
//   auto mu_L = study.mu * study.L_X;
  
//   auto total_weight = 0.0;
//   for (auto& R : study.candidates) {
    
//     CHECK_LE(R.t_max, t_X);
//     CHECK_LE(R.t_min, R.t_max);
    
//     auto m = R.min_muts;
//     auto delta_min = mu_L * (t_X - R.t_max);
//     auto delta_max = mu_L * (t_X - R.t_min);
    
//     R.rel_weight = boost::math::gamma_p(m+1, delta_max) - boost::math::gamma_p(m+1, delta_min);
//     for (auto mm = study.all_min_muts+1; mm <= m; ++mm) {
//       R.rel_weight *= mm / (3 * study.L_X);
//     }
//     total_weight += R.rel_weight;
//   }
//   study.total_weight = total_weight;
// }

// auto pick_insertion_point(const Spr_study& study, absl::BitGenRef bitgen) -> Spr_study_result {
//   // Choose a candidate region in proportion to its weight.
//   //
//   // We usually only ever pick one insertion point for a given SPR study, so there's
//   // no better way to randomly draw a candidate than to scan all the candidates.
//   auto r = absl::Uniform<double>(bitgen, 0.0, study.total_weight);
//   auto chosen_candidate_i = 0;
//   for (auto i = 0; i != std::ssize(study.candidates); ++i) {
//     const auto& candidate = study.candidates[i];
//     if (candidate.rel_weight >= r) {
//       chosen_candidate_i = i;
//       break;
//     } else {
//       r -= candidate.rel_weight;
//     }
//   }

//   // Choose a time t within that region
//   //
//   // Time t should be chosen with probability proportional to:
//   //
//   //    p(t) ~ exp[-mu L_X (t_X - t)] * [mu (t_X - t) / 3]^m  (see calc_weights above)
//   //
//   // with t_min <= t <= t.
//   //
//   // Thankfully, the cdf of this function is related to the incomplete gamma functions.  Define
//   //
//   //        W(t) := int_{t_min}^t dt p(t)
//   //              =  [1/(mu L_X)] * m! / (3 L_X)^m * [P(m+1, delta_max) - P(m+1, delta)],
//   // with
//   //       delta := mu L_X (t_X - t)
//   //   delta_max := mu L_X (t_X - t_min)
//   //
//   // Note that W(t_min) = 0.  To draw t, we just draw u ~ Unif(0, W(t_max)), and solve u = W(t).
//   // Equivalently, define
//   //
//   //       W'(t) := [P(m+1, delta_max) - P(m+1, delta)],
//   //
//   // then choose v ~ Unif(0, W'(t_max)) and solve v = W'(t):
//   //
//   //  P(m+1, delta) = P(m+1, delta_max) - v    [w := rhs is also ~ Unif{0, W'(t_max)}]
//   //      => delta  = P^{-1}_{m+1} [ w ].
//   //
//   // The final inverse function is available in Boost as gamma_p_inv(m+1, w):
//   // - https://www.boost.org/doc/libs/1_84_0/libs/math/doc/html/math_toolkit/sf_gamma/igamma_inv.html
//   //
//   const auto& candidate = study.candidates[chosen_candidate_i];
//   auto m = candidate.min_muts;
//   auto mu_L = study.mu * study.L_X;
//   auto delta_max = mu_L * (study.t_X - candidate.t_min);
//   auto delta_min = mu_L * (study.t_X - candidate.t_max);
//   auto u = absl::Uniform(bitgen, 0.0, 1.0);
//   auto w_max = boost::math::gamma_p(m+1, delta_max);
//   auto w_min = boost::math::gamma_p(m+1, delta_min);
//   auto w = u * (w_min - w_max);
//   auto delta = boost::math::gamma_p_inv(m+1, w + w_max);
//   auto t = study.t_X - (delta / mu_L);

//   //std::cerr << absl::StreamFormat(
//   //    "Candidate %d (branch %d, %g <= t <= %g, m = %d), mu_L=%g, t_X=%g, delta_max=%g, delta_min=%g, w_max=%g, w_min=%g, w=%g, delta=%g, t=%g\n",
//   //    chosen_candidate_i, candidate.branch, candidate.t_min, candidate.t_max, m, mu_L, study.t_X, delta_max, delta_min, w_max, w_min, w, delta, t);
  
//   // Poor-man's defense against roundoff
//   CHECK_GT(t, candidate.t_min);
//   CHECK_LT(t, candidate.t_max);
//   //t = std::max(candidate.t_min, std::min(candidate.t_max, t));

//   return {chosen_candidate_i, t};
// }

// auto calc_log_pdf_at_insertion_point(const Spr_study& study, const Spr_study_result& result) -> double {
//   // This is relatively easy in terms of the relative weights calculated by calc_weights above.
//   // The normalization factor is only important when comparing choosing two insertion points
//   // from *different* SPR studies
  
//   auto mu_L = study.mu * study.L_X;
//   const auto& candidate = study.candidates[result.chosen_candidate_idx];
//   auto m = candidate.min_muts;
//   auto m_min = study.all_min_muts;
  
//   return                                                     // Before log:
//       -mu_L * (study.t_X - result.t)                         //   exp[-mu L_X (t_X - t)]
//       + m * std::log(study.mu * (study.t_X - result.t) / 3)  // * [mu (t_X - t) / 3]^m
//       - (std::log(study.total_weight)                        // / [total_weight
//          - std::log(mu_L)                                    //      * [1/(mu L_X)]
//          + std::lgamma(m_min+1)                              //       * (m_min)!
//          - m_min*std::log(3 * study.L_X));                   //       / (3 L_X)^(m_min)]
// }

// // Find the index of the candidate region containing the query point.  Returns -1 if none of them do
// auto find_candidate_i(const Spr_study& study, const Phylo_tree_loc& query) -> int {
//   // For now, nothing better than linear scan.  Maybe there's a way of organizing the data structures
//   // differently
//   for (auto i = 0; i != std::ssize(study.candidates); ++i) {
//     const auto& candidate = study.candidates[i];
//     if (candidate.branch == query.branch && candidate.t_min <= query.t && query.t <= candidate.t_max) {
//       return i;
//     }
//   }
//   return -1;
// }

}  // namespace delphy
