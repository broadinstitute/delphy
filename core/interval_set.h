#ifndef DELPHY_INTERVAL_SET_H_
#define DELPHY_INTERVAL_SET_H_

#include <boost/container/flat_map.hpp>

#include <absl/log/check.h>
#include <cppcoro/generator.hpp>

#include "sequence.h"
#include "scratch_space.h"

namespace delphy {

// An interval set is an ordered vector of non-overlapping [start, end) pairs describing a subset of a number line.
//
// This data structure allows for efficient single-pass implementations of the following operations:
// * Querying whether a number is in the interval set
// * Computing the intersection of two interval sets
// * Computing the union of two interval sets
// * Computing the difference of two interval sets
//
// Boost has one implementation of interval sets (Boost::icl), but they are impenetrable and are hard-coded
// to use std::allocator.  A lot of our performance (especially on the WASM side) comes from using
// Scratch_space as a temporary memory allocator.  So instead of trying to fight Boost, we've rolled our own here.

using Site_interval = std::pair<Site_index, Site_index>;

template<typename Allocator = std::allocator<Site_interval>>
class Interval_set {
 private:
  using container_type = boost::container::flat_map<Site_index, Site_index, std::less<Site_index>, Allocator>;
  
 public:
  using value_type = Site_interval;
  using iterator = container_type::iterator;
  using const_iterator = container_type::const_iterator;
      
  Interval_set(const Allocator& allocator = Allocator{}) : intervals_{allocator} {}
  Interval_set(std::initializer_list<Site_index> init_sites,
               const Allocator& allocator = Allocator{}) : Interval_set{allocator} {
    for (const auto& l : init_sites) {
      insert(l);
    }
  }
  Interval_set(std::initializer_list<Site_interval> init_intervals,
               const Allocator& allocator = Allocator{}) : Interval_set{allocator} {
    for (const auto& interval : init_intervals) {
      insert(interval);
    }
  }
  Interval_set(const Interval_set& rhs) = default;
  auto operator=(const Interval_set& rhs) -> Interval_set& = default;
  template<typename Alloc2>  // Allow copy assignment from an Interval_set with a different allocator
  auto operator=(const Interval_set<Alloc2>& rhs) -> Interval_set& {
    intervals_.clear();
    intervals_.insert(boost::container::ordered_unique_range_t{},  // Marker = what we're inserting is ordered & unique
                      rhs.begin(), rhs.end());
    return *this;
  }
  Interval_set(Interval_set&& rhs) = default;
  auto operator=(Interval_set&& rhs) -> Interval_set& {
    intervals_.adopt_sequence(boost::container::ordered_unique_range_t{}, rhs.intervals_.extract_sequence());
    return *this;
  }
  ~Interval_set() = default;

  auto operator==(const Interval_set& that) const -> bool {
    return this->intervals_ == that.intervals_;
  }

  auto begin() const -> const_iterator { return intervals_.begin(); }
  auto end() const -> const_iterator { return intervals_.end(); }

  // Don't use this inside a tight loop; instead, treat each interval in one go
  auto slow_elements() const -> cppcoro::generator<Site_index> {
    for (const auto& [start, end] : intervals_) {
      for (auto l = start; l != end; ++l) {
        co_yield l;
      }
    }
  }

  auto size() const -> std::size_t { return num_sites(); }  // CAREFUL! Prefer num_intervals() or num_sites()
  auto num_intervals() const -> int { return std::ssize(intervals_); }
  auto num_sites() const -> int {
    auto result = 0;
    for (const auto& [start, end] : intervals_) {
      result += end - start;
    }
    return result;
  }
  
  auto clear() -> void { intervals_.clear(); }
  auto empty() const -> bool { return intervals_.empty(); }
  
  auto insert(const Site_interval& interval) -> std::pair<const_iterator, bool> {
    auto [overlap_first, overlap_last] = overlapping_intervals_closed(interval);  // NOTE: `_closed`!
    if (overlap_first == overlap_last) {
      // Easy case: no overlaps
      return intervals_.insert(interval);
    } else {
      // Harder!  Either subsume `interval` into existing interval,
      // or replace overlapping intervals with their union with `interval`

      if (overlap_last == std::ranges::next(overlap_first)) {
        // `interval` overlaps a single interval; is it completely contained within it?
        auto [overlap_start, overlap_end] = *overlap_first;
        if (overlap_start <= interval.first && interval.second <= overlap_end) {
          // Completely contained
          return {overlap_first, false};
        }
      }

      // If we're here, we have to replace the overlapping intervals
      
      // Get the *actual* final overlapping interval, which exists because overlap_first != overlap_last
      auto real_overlap_last = std::ranges::prev(overlap_last);
      
      auto new_interval_start = std::min(interval.first, overlap_first->first);
      auto new_interval_end = std::max(interval.second, real_overlap_last->second);

      intervals_.erase(overlap_first, overlap_last);
      return intervals_.insert({new_interval_start, new_interval_end});
    }
  }
  auto insert(Site_index l) -> std::pair<const_iterator, bool> {
    return insert(Site_interval{l, l+1});
  }

  auto contains(Site_index l) const -> bool {
    auto it = intervals_.upper_bound(l);  // First interval with start > l
    if (it == intervals_.begin()) { return false; }  // All intervals start strictly after `l`
    --it;   // Now `it` points to the last interval with start <= l
    auto [start, end] = *it;
    DCHECK_LE(start, l);
    return l < end;
  }

  // Return the range of intervals that overlaps the query interval, treating the query and all target
  // intervals as CLOSED
  auto overlapping_intervals_closed(Site_interval query) -> std::pair<iterator, iterator> {
    auto [qstart, qend] = query;
    DCHECK_LE(qstart, qend);
    
    // One of two possibilities for qstart wrt upper_bound(qstart)
    //
    //  1.             qstart
    //                   |
    //                   V
    //   -----[-----]------[-----]-----
    //                     ^
    //                     |
    //                  it_start
    //
    auto it_start = intervals_.lower_bound(qstart);  // First interval with qstart <= start
    
    // Perhaps the preceding interval also overlaps `query`
    //
    //  2.     qstart
    //           |
    //           V
    //   -----[-----]------[-----]-----
    //        ^            ^
    //        |            |
    //     prec_it      it_start
    //
    if (it_start != intervals_.begin()) {
      auto prec_it = std::ranges::prev(it_start);  // Last interval with qstart > start
      auto [prec_start, prec_end] = *prec_it;
      DCHECK_LT(prec_start, qstart);
      if (qstart <= prec_end) {  // NOTE `<=` !!!
        // We're in case 2: qstart is inside the interval prec_it
        it_start = prec_it;
      }
    }

    // At this point, it_start points to the first interval overlapping with [qstart, infty), if any

    // The end-point is much easier: upper_bound(qend) directly gives us an iterator to the first
    // interval after the query range
    auto it_end = intervals_.upper_bound(qend);      // First interval with qend < start

    DCHECK(it_start <= it_end);
    
    return {it_start, it_end};
  }

  // Easy access to internals if you know what you're doing
  auto raw_intervals() -> container_type& { return intervals_; }  // Direct access if you know what you're doing
  auto raw_intervals() const -> const container_type& { return intervals_; }
  
  auto no_consecutive_intervals() const -> bool {
    if (empty()) { return true; }
    
    auto it = intervals_.begin();
    auto it2 = std::ranges::next(it);
    
    while (it2 != intervals_.end()) {
      auto [prev_start, prev_end] = *it;
      auto [cur_start, cur_end] = *it2;
      if (cur_start <= prev_end) { return false; }
      it = it2;
      ++it2;
    }
    
    return true;
  }
  
 private:
  // Interval container, stored in a sorted vector in the form of [start, end) pairs
  container_type intervals_;
};

using Scratch_interval_set = Interval_set<Stl_scratch_space<Site_interval>>;

template<typename Allocator>
inline auto operator<<(std::ostream& os, const Interval_set<Allocator>& intervals) -> std::ostream& {
  os << "[";
  for (const auto& [start, end] : intervals) {
    os << absl::StreamFormat("[%d,%d),", start, end);
  }
  os << "]";
  return os;
}

template<typename Alloc_dst, typename Alloc_A, typename Alloc_B>
auto merge_interval_sets(
    Interval_set<Alloc_dst>& dst,
    const Interval_set<Alloc_A>& src_A,
    const Interval_set<Alloc_B>& src_B)
    -> void {

  DCHECK(src_A.no_consecutive_intervals());
  DCHECK(src_B.no_consecutive_intervals());
  
  dst.clear();
  auto intervals = dst.raw_intervals().extract_sequence();

  auto inside = false;
  auto cur_start = Site_index{0};
  auto cur_end = cur_start;
  
  auto it_A = src_A.begin();
  auto it_B = src_B.begin();
  while (not (it_A == src_A.end() && it_B == src_B.end())) {
    auto& first_it = (it_A == src_A.end() ? it_B :  // Note type is `auto&`, so advancing will affect underlying iter
                      it_B == src_B.end() ? it_A :
                      it_A->first <= it_B->first ? it_A : it_B);
    auto [first_start, first_end] = *first_it;
    
    if (not inside) {
      cur_start = first_start;
      cur_end = first_end;
      ++first_it;
      inside = true;
    } else {
      // If first_it can extend cur interval, do it
      if (first_start <= cur_end) {
        cur_end = std::max(cur_end, first_end);
        ++first_it;
      } else {
        // Otherwise, all remaining intervals start beyond cur_end, so output [cur_start, cur_end)
        intervals.push_back(Site_interval{cur_start, cur_end});
        inside = false;
      }
    }
  }

  // Output last interval, if we're inside it
  if (inside) {
    intervals.push_back(Site_interval{cur_start, cur_end});
  }

  dst.raw_intervals().adopt_sequence(boost::container::ordered_unique_range_t{}, std::move(intervals));
  DCHECK(dst.no_consecutive_intervals());
}

template<typename Alloc_A, typename Alloc_B, typename Alloc_result = std::allocator<Site_interval>>
auto merge_interval_sets(
    const Interval_set<Alloc_A>& src_A,
    const Interval_set<Alloc_B>& src_B,
    const Alloc_result& result_allocator = Alloc_result{})
    -> Interval_set<Alloc_result> {
  auto result = Interval_set<>{result_allocator};
  merge_interval_sets(result, src_A, src_B);
  return result;
}

template<typename Alloc_dst, typename Alloc_A, typename Alloc_B>
auto intersect_interval_sets(
    Interval_set<Alloc_dst>& dst,
    const Interval_set<Alloc_A>& src_A,
    const Interval_set<Alloc_B>& src_B)
    -> void {

  DCHECK(src_A.no_consecutive_intervals());
  DCHECK(src_B.no_consecutive_intervals());
  
  dst.clear();
  auto intervals = dst.raw_intervals().extract_sequence();

  auto it_A = src_A.begin();
  auto it_B = src_B.begin();
  while (not (it_A == src_A.end() || it_B == src_B.end())) {
    auto [start_A, end_A] = *it_A;
    auto [start_B, end_B] = *it_B;

    auto start_overlap = std::max(start_A, start_B);
    auto end_overlap = std::min(end_A, end_B);
    auto overlap = start_overlap < end_overlap;

    if (overlap) {
      intervals.push_back(Site_interval{start_overlap, end_overlap});
    }
    
    if (end_A <= end_B) {
      ++it_A;
    } else {
      ++it_B;
    }
  }

  dst.raw_intervals().adopt_sequence(boost::container::ordered_unique_range_t{}, std::move(intervals));
  DCHECK(dst.no_consecutive_intervals());
}

template<typename Alloc_A, typename Alloc_B, typename Alloc_result = std::allocator<Site_interval>>
auto intersect_interval_sets(
    const Interval_set<Alloc_A>& src_A,
    const Interval_set<Alloc_B>& src_B,
    const Alloc_result& result_allocator = Alloc_result{})
    -> Interval_set<Alloc_result> {
  auto result = Interval_set<>{result_allocator};
  intersect_interval_sets(result, src_A, src_B);
  return result;
}


template<typename Alloc1, typename Alloc2>
auto interval_sets_intersect(
    const Interval_set<Alloc1>& A,
    const Interval_set<Alloc2>& B)
    -> bool {

  DCHECK(A.no_consecutive_intervals());
  DCHECK(B.no_consecutive_intervals());
  
  auto it_A = A.begin();
  auto it_B = B.begin();
  while (not (it_A == A.end() || it_B == B.end())) {
    auto [start_A, end_A] = *it_A;
    auto [start_B, end_B] = *it_B;

    auto start_overlap = std::max(start_A, start_B);
    auto end_overlap = std::min(end_A, end_B);
    auto overlap = start_overlap < end_overlap;

    if (overlap) {
      return true;
    }
    
    if (end_A <= end_B) {
      ++it_A;
    } else {
      ++it_B;
    }
  }

  return false;
}

template<typename Alloc_dst, typename Alloc_A, typename Alloc_B>
auto subtract_interval_sets(
    Interval_set<Alloc_dst>& dst,         // dst = src_A - src_B
    const Interval_set<Alloc_A>& src_A,
    const Interval_set<Alloc_B>& src_B)
    -> void {

  DCHECK(src_A.no_consecutive_intervals());
  DCHECK(src_B.no_consecutive_intervals());

  dst.clear();

  if (src_A.empty()) { return; }
  
  auto intervals = dst.raw_intervals().extract_sequence();

  auto it_A = src_A.begin();
  auto it_B = src_B.begin();
  auto [cur_start, cur_end] = *it_A;  // src_A.empty() == false
  auto nextA = [&](){
    ++it_A;
    if (it_A != src_A.end()) {
      cur_start = it_A->first;
      cur_end = it_A->second;
    }
  };
  while (it_A != src_A.end()) {
    if (it_B == src_B.end()) {
      intervals.push_back(Site_interval{cur_start, cur_end});
      nextA();
    } else {
      auto [b_start, b_end] = *it_B;
      //
      // Possible relations of `b` interval to `cur` interval:
      //
      // ------[-----)--- cur
      // --[-)----------- b    (1)
      // --[------)------ b    (2)
      // --[-----------)- b    (3)
      // --------[-)----- b    (4)
      // --------[-----)- b    (5)
      // --------------[- b    (6)
      //
      if (b_start < cur_start) {
        
        if (b_end <= cur_start) {  // (1)
          ++it_B;
          
        } else if (b_end < cur_end) {  // (2)
          cur_start = b_end;
          ++it_B;
          
        } else { // (3)
          nextA();
        }
        
      } else if (b_start < cur_end) {

        if (cur_start < b_start) {
          intervals.push_back(Site_interval{cur_start, b_start});
        }
        
        if (b_end < cur_end) {  // (4)
          cur_start = b_end;
          ++it_B;
          
        } else {  // (5)
          nextA();
        }
        
      } else { // (6)
        intervals.push_back(Site_interval{cur_start, cur_end});
        nextA();
      }
    }
  }

  dst.raw_intervals().adopt_sequence(boost::container::ordered_unique_range_t{}, std::move(intervals));
  DCHECK(dst.no_consecutive_intervals());
}

template<typename Alloc_A, typename Alloc_B, typename Alloc_result = std::allocator<Site_interval>>
auto subtract_interval_sets(
    const Interval_set<Alloc_A>& src_A,
    const Interval_set<Alloc_B>& src_B,
    const Alloc_result& result_allocator = Alloc_result{})
    -> Interval_set<Alloc_result> {
  auto result = Interval_set<>{result_allocator};
  subtract_interval_sets(result, src_A, src_B);
  return result;
}

}  // namespace delphy

#endif // DELPHY_INTERVAL_SET_H_
