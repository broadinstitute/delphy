#ifndef DELPHY_SCRATCH_SPACE_H_
#define DELPHY_SCRATCH_SPACE_H_

#include <cstddef>
#include <cstdint>
#include <forward_list>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <stdexcept>
#include <deque>
#include <stack>

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"

namespace delphy {

// A Scratch_space is a bump allocator that only allocates new blocks of memory from a contiguous arena.
// It provides a reset operation that implicitly discards everything by allowing reuse of the whole arena.
// If an allocation overflows the arena, a new larger arena is created.  A reference to any old arenas
// is kept to discard them on reset() or on destruction.  The net effect is that the allocator can be used
// as efficient scratch space for a function that gets called repeatedly (e.g., an MCMC move): it will
// quickly allocate the maximum working memory necessary and thereafter permit very fast allocations and free
// deletions at the end of each call to the function.
//
// Adapted from https://docs.microsoft.com/en-us/cpp/standard-library/allocators?view=msvc-170
// and https://indiegamedev.net/2022/03/27/custom-c20-memory-allocators-for-stl-containers/
//
// Here is the simplest way to use these allocators:
// ```
// auto scratch = Scratch_space{};
// for (auto i = 0; i != 1000; ++i) {
//   auto big_v = Scratch_vector<int>{scratch};
//   auto node_data = Scratch_unordered_map<Node_index, Calc_data>{};
//   auto active_nodes = Scratch_unordered_set<Node_index>{};
//
//   ... do lots of collection manipulation ...
//
//   scratch.reset();  // Free everything in 2-3 instructions; future allocations recycle allocated scratch space
// }
//
class Scratch_space {
 public:
  explicit Scratch_space(size_t starting_arena_size = 1024L * 1024L) {
    arena_begin_ = arena_cur_ = arena_end_ = nullptr;
    next_arena_size_ = starting_arena_size;
  }
  Scratch_space(const Scratch_space&) = delete;
  Scratch_space& operator=(const Scratch_space&) = delete;
  Scratch_space(Scratch_space&& other) noexcept
      : arena_begin_{other.arena_begin_},
        arena_cur_{other.arena_cur_},
        arena_end_{other.arena_end_},
        next_arena_size_{other.next_arena_size_},
        incinerator_{std::move(other.incinerator_)} {
    other.arena_begin_ = other.arena_cur_ = other.arena_end_ = nullptr;
  }
  Scratch_space& operator=(Scratch_space&& other) noexcept {
    incinerate_old_arenas();
    delete[] arena_begin_;
    arena_begin_ = other.arena_begin_;
    arena_cur_ = other.arena_cur_;
    arena_end_ = other.arena_end_;
    next_arena_size_ = other.next_arena_size_;
    incinerator_ = std::move(other.incinerator_);
    other.arena_begin_ = other.arena_cur_ = other.arena_end_ = nullptr;
    return *this;
  }
  ~Scratch_space() noexcept {
    incinerate_old_arenas();
    delete[] arena_begin_;
  }

  auto allocate(size_t num_bytes) -> void* {
    // Super-fast happy path...
    // Ensure arena_cur_ stays aligned to max_align_t
    num_bytes = (alignof(max_align_t) - 1 + num_bytes) & -alignof(max_align_t);
    auto result = arena_cur_;
    arena_cur_ = reinterpret_cast<std::byte*>(reinterpret_cast<uintptr_t>(arena_cur_) + num_bytes);
    if (arena_cur_ <= arena_end_) {
      return result;
    } else {
      // Slow unhappy path (in a different function to allow inlining of allocate in the common case)
      return allocate_in_new_arena(num_bytes);
    }
  }
  auto free(void*) -> void {
    // Do nothing: call reset() to free all the memory ever allocated by this allocator in one go
  }
  auto reset() -> void {
    incinerate_old_arenas();
    arena_cur_ = arena_begin_;
  }

 private:
  // All three arena pointers are always aligned to max_align_t (i.e., typically, at a 16-byte boundary)
  std::byte* arena_begin_;
  std::byte* arena_cur_;
  std::byte* arena_end_;
  size_t next_arena_size_;

  std::forward_list<std::byte*> incinerator_;

  auto incinerate_old_arenas() -> void {
    for (auto old_arena_begin : incinerator_) {
      delete[] old_arena_begin;
    }
    incinerator_.clear();
  }
  auto allocate_in_new_arena(size_t num_bytes) -> void* {
    incinerator_.push_front(arena_begin_);
    if (num_bytes > next_arena_size_) {
      next_arena_size_ = num_bytes;
    }
    // Pointer returned by new is guaranteed aligned to max_align_t (usually, 16 bytes)
    // See https://stackoverflow.com/a/54590867
    arena_begin_ = new std::byte[next_arena_size_];
    arena_end_ = arena_begin_ + next_arena_size_;
    next_arena_size_ *= 2;

    arena_cur_ = arena_begin_ + num_bytes;
    return arena_begin_;
  }
};

template<typename T = int>  // default for ease of use when allocator will just be rebound to a different type
class Stl_scratch_space {
 public:
  using value_type = T;

  // Implicit constructor allows us to pass the scratch_space directly to suitably-typed STL constructors,
  // e.g., instead of `auto v = std::vector<int, Stl_scratch_space<int>>{Stl_scratch_space(scratch_space)};`
  //          can say `auto v = std::vector<int, Stl_scratch_space<int>>{scratch_space};`
  Stl_scratch_space(Scratch_space& scratch_space) noexcept  // NOLINT(google-explicit-constructor)
      : scratch_space_(&scratch_space) {}
  template<typename U>
  explicit Stl_scratch_space(const Stl_scratch_space<U>& rhs) noexcept : scratch_space_(rhs.scratch_space_) {}

  auto allocate(std::size_t n) -> T* {
    return static_cast<T*>(scratch_space_->allocate(n * sizeof(T)));
  }
  auto deallocate(T* p, std::size_t /*n*/) noexcept -> void { scratch_space_->free(p); }

  auto operator==(const Stl_scratch_space<T>& rhs) const noexcept -> bool = default;

 private:
  template<typename U> friend
  class Stl_scratch_space;
  Scratch_space* scratch_space_;
};

template<typename T>
using Scratch_vector = std::vector<T, Stl_scratch_space<T>>;

template<typename T>
using Scratch_stack = std::stack<T, Scratch_vector<T>>;

template<typename T>
using Scratch_deque = std::deque<T, Stl_scratch_space<T>>;

template<typename K,
         typename V,
         typename Hash = typename absl::flat_hash_map<K, V>::hasher,
         typename Eq = typename absl::flat_hash_map<K, V, Hash>::key_equal>
using Scratch_flat_hash_map = absl::flat_hash_map<K, V, Hash, Eq, Stl_scratch_space<std::pair<const K, V>>>;

template<typename T,
         typename Hash = typename absl::flat_hash_set<T>::hasher,
         typename Eq = typename absl::flat_hash_set<T, Hash>::key_equal>
using Scratch_flat_hash_set = absl::flat_hash_set<T, Hash, Eq, Stl_scratch_space<T>>;

}  // namespace delphy

#endif // DELPHY_SCRATCH_SPACE_H_
