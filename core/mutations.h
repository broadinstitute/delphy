#ifndef DELPHY_MUTATIONS_H_
#define DELPHY_MUTATIONS_H_

#include <cassert>
#include <vector>

#include <absl/container/flat_hash_set.h>
#include <absl/container/flat_hash_map.h>
#include <absl/strings/str_format.h>
#include <absl/strings/str_join.h>
#include <absl/log/check.h>

#include "sequence.h"
#include "scratch_space.h"
#include "interval_set.h"
#include "tree.h"

namespace delphy {

// A mutation at a particular site and time
struct Mutation {
  // TODO: Reorder fields for packing
  Real_seq_letter from;
  Site_index site;
  Real_seq_letter to;
  double t;
  
  auto operator<=>(const Mutation& that) const = default;
};
inline auto operator<<(std::ostream& os, const Mutation& m) -> std::ostream& {
  return os << absl::StreamFormat("%c%d%c @ %g", to_char(m.from), m.site, to_char(m.to), m.t);
}

// TODO: Make a specialized container for mutation lists that handles 0-2 mutations efficiently
template<typename Alloc = std::allocator<Mutation>>
using Mutation_list = std::vector<Mutation, Alloc>;
using Scratch_mutation_list = Mutation_list<Stl_local_arena_allocator<Mutation>>;

// Mutations on a branch should always be sorted in increasing order of time.
// To facilitate testing, we also break ties by sorting according to site.
inline auto mutations_on_branch_less(const Mutation& ma, const Mutation& mb) -> bool {
  return std::tie(ma.t, ma.site) < std::tie(mb.t, mb.site);
}
template<typename Mutations>
auto sort_mutations(Mutations& ms) -> void {
  std::ranges::sort(ms, mutations_on_branch_less);
}
template<typename Mutations>
auto are_mutations_sorted(Mutations& ms) -> bool {
  return std::ranges::is_sorted(ms, mutations_on_branch_less);
}

// Mutations on a branch should be within the start and end time of a branch.
// Sometimes, round-off errors violate this condition.  clamp_mutation_times cleans this up.
template<typename Mutations>
auto clamp_mutation_times(Mutations& ms, double min_t, double max_t) -> void {
  for (auto& m : ms) {
    m.t = std::clamp(m.t, min_t, max_t);
  }
}
template<typename Mutations>
auto are_mutation_times_in_range(Mutations& ms, double min_t, double max_t) -> bool {
  return std::ranges::all_of(ms, [=](const auto& m){ return min_t <= m.t && m.t <= max_t; });
}

// A Seq_delta is like a mutation, but without an associated time
struct Seq_delta {
  Site_index site;
  Real_seq_letter from;
  Real_seq_letter to;

  Seq_delta() {}
  Seq_delta(Site_index site, Real_seq_letter from, Real_seq_letter to) : site{site}, from{from}, to{to} {}
  Seq_delta(const Mutation& m) : Seq_delta(m.site, m.from, m.to) {}

  auto operator<=>(const Seq_delta& that) const = default;

  auto inverse() const -> Seq_delta { return {site, to, from}; }
  auto place_at(double t) const -> Mutation { return Mutation{from, site, to, t}; }
};
inline auto operator<<(std::ostream& os, const Seq_delta& delta) -> std::ostream& {
  return os << absl::StreamFormat("%c%d%c", to_char(delta.from), delta.site, to_char(delta.to));
}

// "Missation": a marker on a branch that says: from the beginning of this branch down,
// there is no data about this site.  By integrating out the parts of the likelihood
// on nodes below this branch, we reduce the inference problem to one where we have
// complete data at every site, albeit with trees at each site being subtrees of the full
// phylo tree.  We call this idea "N-pruning".
using Missation = std::pair<const Site_index, Real_seq_letter>;
inline auto operator<<(std::ostream& os, const Missation& mi) -> std::ostream& {
  const auto& [site, from] = mi;
  return os << absl::StreamFormat("%c%dN", to_char(from), site);
}

// A Missation_map encodes which sites are missing in sequences downstream of a point X on a phylo tree.
// It also encodes the state at those sites at X.
//
// The naive implementation, which we used for a long time, is a simple map of "site index" to "state".
// A site is missing if it's in the map, and its state at X is the value in the map.
//
// However, data in real sequences is rarely missing in arbitary sites.  Instead, it comes in large
// contiguous "gaps", i.e., a series of consecutive sites.  Moreover, the state of any one site at X is
// very unlikely to differ from the state of the reference sequence of the tree:
// a mutation between X and the point above the root is a rare event.
//
// Hence, we can represent missation maps much more efficiently as a combination of two structures:
// * An interval set of gaps, i.e., a vector of [start,end) pairs describing the ranges of site indices that are missing;
// * A map containing only the delta between the states of missing sites at X vs in the reference sequence.
//
// The following operations can then be implemented very efficiently:
// * Query whether a site is missing
// * Query the state of a missing site at X
// * Update the state of a missing site at X
//
// Further, efficient operations are available on the underlying interval set.  See interval_set.h.
//
// We don't store a pointer to the reference sequence in the Missation_map, since it's the
// same sequence for all maps in a tree.  Instead, we take in the reference sequence as a parameter
// in many calls.  It is the callers responsibility to ensure that the right reference sequence is
// passed in, and that changes in the reference sequence are matched by counteracting changes in
// the Missation_map.
//
template<typename Allocator = std::allocator<Site_interval>>
class Missation_map {
 private:
  using From_states_allocator =
      typename std::allocator_traits<Allocator>::template rebind_alloc<Missation>;
  using From_states_map =
      boost::container::flat_map<Site_index, Real_seq_letter, std::less<Site_index>, From_states_allocator>;
  
 public:
  Interval_set<Allocator> intervals;
  From_states_map from_states;  // A site present here *must* also be contained in `intervals`

  Missation_map(const Allocator& alloc = Allocator{})
      : intervals{alloc}, from_states{From_states_allocator{alloc}} {}
  Missation_map(std::initializer_list<Missation> init_missations,
                const Real_sequence& ref_sequence,
                const Allocator& allocator = Allocator{}) : Missation_map{allocator} {
    for (const auto& mi : init_missations) {
      insert(mi, ref_sequence);
    }
  }
  Missation_map(const Missation_map& that) = default;
  auto operator=(const Missation_map& that) -> Missation_map& = default;
  template<typename Alloc2>
  auto operator=(const Missation_map<Alloc2>& that) -> Missation_map& {
    this->intervals = that.intervals;
    this->from_states.clear();
    this->from_states.insert(boost::container::ordered_unique_range_t{},  // == what we're inserting is ordered & unique
                             that.from_states.begin(), that.from_states.end());
    return *this;
  }
  Missation_map(Missation_map&& that) = default;
  auto operator=(Missation_map&& that) -> Missation_map& = default;
  ~Missation_map() = default;
  

  auto operator==(const Missation_map& that) const -> bool {
    return this->intervals == that.intervals &&
        this->from_states == that.from_states;
  }

  // Don't use this inside a tight loop; instead, treat each interval in one go
  auto slow_elements(const Real_sequence& ref_sequence) const -> cppcoro::generator<Missation> {
    for (const auto& [start, end] : intervals) {
      for (auto l = start; l != end; ++l) {
        co_yield Missation{l, get_from_state(l, ref_sequence)};
      }
    }
  }
  
  auto empty() const -> bool { return intervals.empty(); }
  auto clear() -> void {
    intervals.clear();
    from_states.clear();
  }

  auto size() const -> std::size_t { return intervals.num_sites(); }  // CAREFUL! Prefer num_intervals() or num_sites()
  auto num_intervals() const -> int { return intervals.num_intervals(); }
  auto num_sites() const -> int { return intervals.num_sites(); }

  auto contains(Site_index l) const -> bool { return intervals.contains(l); }
  auto insert(Missation mi, const Real_sequence& ref_sequence) -> void {
    auto [mi_site, mi_from] = mi;
    intervals.insert(mi_site);
    set_from_state(mi_site, mi_from, ref_sequence);
  }
  auto get_from_state(Site_index l, const Real_sequence& ref_sequence) const -> Real_seq_letter {
    DCHECK(intervals.contains(l));
    if (auto it = from_states.find(l); it != from_states.end()) {
      return it->second;
    } else {
      return ref_sequence[l];
    }
  }
  auto set_from_state(Site_index l, Real_seq_letter from, const Real_sequence& ref_sequence) -> void {
    if (l < 0 || l >= std::ssize(ref_sequence)) {
      throw std::out_of_range(absl::StrFormat(
          "Missation out of range, site %d outside of range [0,%d)", l, std::ssize(ref_sequence)));
    }
    DCHECK(intervals.contains(l));
    auto ref_from = ref_sequence[l];
    if (from != ref_from) {
      from_states.insert_or_assign(l, from);
    } else {
      from_states.erase(l);  // harmless if not there already
    }
  }
  
  auto ref_seq_changed(Site_index l, Real_seq_letter old_ref_state, Real_seq_letter new_ref_state) -> void {
    if (not intervals.contains(l)) { return; }  // Ignore irrelevant ref seq changes
    
    if (auto it = from_states.find(l); it != from_states.end()) {
      auto from_state = it->second;
      DCHECK_NE(from_state, old_ref_state);
      if (from_state == new_ref_state) {
        from_states.erase(it);
      } else {
        // Still a delta from the reference
      }
    } else {
      // Not there, but maybe it should be
      auto from_state = old_ref_state;
      if (from_state == new_ref_state) {
        // Still not a delta from the reference
      } else {
        from_states.insert(Missation{l, old_ref_state});
      }
    }
  }
};
template<typename Allocator>
inline auto operator<<(std::ostream& os, const Missation_map<Allocator>& mis) -> std::ostream& {
  return os << absl::StreamFormat(
      "Missation_map{intervals=%s, from_states=[%s]}",
      absl::FormatStreamed(mis.intervals),
      absl::StrJoin(mis.from_states, ", ", absl::StreamFormatter()));
}

using Scratch_missation_map = Missation_map<Stl_local_arena_allocator<Site_interval>>;

template<
  typename Alloc_mi_A,
  typename Alloc_mi_B,
  typename Alloc_result_A,
  typename Alloc_result_B,
  typename Alloc_result_common>
auto factor_out_common_missations(
    const Missation_map<Alloc_mi_A>& mi_A,
    const Missation_map<Alloc_mi_B>& mi_B,
    Missation_map<Alloc_result_A>& result_A,
    Missation_map<Alloc_result_B>& result_B,
    Missation_map<Alloc_result_common>& result_common)
    -> void {
  
  DCHECK(result_A.empty());
  DCHECK(result_B.empty());
  DCHECK(result_common.empty());

  // Intervals
  intersect_interval_sets(result_common.intervals, mi_A.intervals, mi_B.intervals);
  subtract_interval_sets(result_A.intervals, mi_A.intervals, result_common.intervals);
  subtract_interval_sets(result_B.intervals, mi_B.intervals, result_common.intervals);

  // `From` states
  auto it_A = mi_A.from_states.begin();
  auto it_B = mi_B.from_states.begin();
  while (it_A != mi_A.from_states.end() && it_B != mi_B.from_states.end()) {
    auto [mi_A_site, mi_A_from] = *it_A;
    auto [mi_B_site, mi_B_from] = *it_B;

    if (mi_A_site < mi_B_site) {
      result_A.from_states.insert(*it_A);
      ++it_A;

    } else if (mi_B_site < mi_A_site) {
      result_B.from_states.insert(*it_B);
      ++it_B;

    } else {
      DCHECK_EQ(mi_A_site, mi_B_site);
      DCHECK_EQ(mi_A_from, mi_B_from);
      result_common.from_states.insert(*it_A);
      ++it_A;
      ++it_B;
    }
  }
  while (it_A != mi_A.from_states.end()) {
    result_A.from_states.insert(*it_A);
    ++it_A;
  }
  while (it_B != mi_B.from_states.end()) {
    result_B.from_states.insert(*it_B);
    ++it_B;
  }
}

template<typename Alloc_src, typename Alloc_common, typename Alloc_scratch = std::allocator<Site_interval>>
auto factor_out_common_missations(
    Missation_map<Alloc_src>& mi_A,
    Missation_map<Alloc_src>& mi_B,
    Missation_map<Alloc_common>& result_common,
    const Alloc_scratch& scratch_allocator = Alloc_scratch{})
    -> void {
  auto result_A = Missation_map<Alloc_scratch>{scratch_allocator};
  auto result_B = Missation_map<Alloc_scratch>{scratch_allocator};
  factor_out_common_missations(mi_A, mi_B, result_A, result_B, result_common);
  mi_A = result_A;
  mi_B = result_B;
}

// mi_dst = (mi_src_A + mi_src_B)   (assuming no conflict of `from` states in overlapping sites)
template<typename Alloc_dst, typename Alloc_A, typename Alloc_B>
auto merge_missations_nondestructively(
    Missation_map<Alloc_dst>& mi_dst,
    const Missation_map<Alloc_A>& mi_src_A,
    const Missation_map<Alloc_B>& mi_src_B)
    -> void {
  merge_interval_sets(mi_dst.intervals, mi_src_A.intervals, mi_src_B.intervals);
  mi_dst.from_states.clear();
  mi_dst.from_states.insert(mi_src_A.from_states.begin(), mi_src_A.from_states.end());
  mi_dst.from_states.insert(mi_src_B.from_states.begin(), mi_src_B.from_states.end());
}

template<typename Alloc_A, typename Alloc_B, typename Alloc_result = std::allocator<Site_interval>>
auto merge_missations_nondestructively(
    const Missation_map<Alloc_A>& mi_src_A,
    const Missation_map<Alloc_B>& mi_src_B,
    const Alloc_result& result_allocator = Alloc_result{})
    -> Missation_map<Alloc_result> {
  auto result = Missation_map<Alloc_result>{result_allocator};
  merge_missations_nondestructively(result, mi_src_A, mi_src_B);
  return result;
}

template<typename Alloc_dst, typename Alloc_A, typename Alloc_B>
auto subtract_missations_nondestructively(
    Missation_map<Alloc_dst>& mi_dst,
    const Missation_map<Alloc_A>& mi_src_A,
    const Missation_map<Alloc_B>& mi_src_B)
    -> void {
  subtract_interval_sets(mi_dst.intervals, mi_src_A.intervals, mi_src_B.intervals);
  mi_dst.from_states.clear();
  mi_dst.from_states.insert(mi_src_A.from_states.begin(), mi_src_A.from_states.end());
  for (const auto& [l, _] : mi_src_B.from_states) {
    mi_dst.from_states.erase(l);
  }
}

}  // namespace delphy

#endif // DELPHY_MUTATIONS_H_
