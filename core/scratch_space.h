#ifndef DELPHY_SCRATCH_SPACE_H_
#define DELPHY_SCRATCH_SPACE_H_

#include <cstddef>
#include <cstdint>
#include <vector>
#include <stack>
#include <deque>

#include "absl/log/check.h"
#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"

namespace delphy {

// An `Arena` is a bump allocator, useful for fast, temporary allocations inside
// inner hot loops (e.g., inside an MCMC move).  Allocations occur inside nested scopes.
// Within a scope, allocation proceeds by bumping.  Deallocation happens in bulk at the
// end of a scope.
//
// A thread-local `local_arena` exists in all threads, and is the `Arena` you should
// use unless you have a specific reason to do otherwise.  For efficiency reasons,
// you *must* be inside a Local_arena_scope (either directly or indirectly via a caller)
// to use `local_arena` (see below).
//
// Here is a simple idiomatic example of using `local_arena`:
// ```
// {
//   auto scope = Local_arena_scope{};
//   auto p = local_arena->allocate(100);  // ~5-6 instructions in the common case
//   ...
//   // `scope` ends here, memory for `p` and any other allocations from
//   // `local_arena` released automatically in 2-3 instructions
// }
// ```
//
// Also provided are an STL-style allocator backed by `local_arena`, as well as convenient
// aliases for common STL containers using this allocator.  For example,
// ```
// {
//   auto scope = Local_arena_scope{};
//   auto vec = Scratch_vector{1,2,3,4};
//   vec.push_back(5);
//   ...
//   // `scope` ends here, memory for `vec` released automatically
// }
// ```

class Arena {
  
 private:
  // A block is a contiguous chunk of memory allocated from the system,
  // from which smaller allocations will be sourced by the `Arena`.
  struct Block_desc {
    char* begin = nullptr;  // Empty by default
    char* end = nullptr;

    auto empty() const -> bool { return begin == nullptr && end == nullptr; }
    auto size() const -> std::size_t { return static_cast<std::size_t>(end - begin); }
    auto contains(char* p) const -> bool { return begin <= p && p <= end; }
  };

  static auto alloc_block(std::size_t num_bytes) -> Block_desc {
    // Pointer returned by new is guaranteed aligned to max_align_t (usually, 16 bytes)
    // See https://stackoverflow.com/a/54590867
    auto mem = new char[num_bytes];
    return {.begin = mem, .end = mem + num_bytes};
  }
  static auto free_block(Block_desc& block) -> void {
    delete [] block.begin;
    block = {.begin = nullptr, .end = nullptr };
  }

  // Allocations are sourced from `cur_block` by bumping `cur`, unless there
  // is not enough memory left in that block
  Block_desc cur_block;
  char* cur;

  // Blocks that contain live allocations but are full are kept in `full_blocks`.
  // They may be released on the next call to `pop_scope()`...
  std::vector<Block_desc> full_blocks;

  // ...but we try to hold on to the largest block ever allocated for longer
  // in case we need to reuse it later.  Ideally, after a few warm-up cycles of
  // allocations and releases, the arena will converge to a single block that
  // is large enough to cover all the needs within a single cycle.
  Block_desc spare_block;

  // We keep a stack of `marks` within the arena (the value of `cur` at the
  // time of a call to `push_scope()`).  Deallocation proceeds in bulk by
  // calling `pop_scope()`, which pops the topmost mark in this stack
  // into `cur` and releases all intervening blocks.  See `Arena_scope` below.
  std::stack<char*, std::vector<char*>> marks;

  
 public:
  Arena(std::size_t initial_block_size = 1024L * 1024L) {
    // Instead of allocating the initial block as `cur_block`, we allocate it
    // as `spare_block`, and initialize `cur_block` to a dummy, 0-byte block.
    // If we instead did the natural thing and initialized `cur_block` to
    // an empty block of size `initial_block_size`, then that block would never
    // be released by a `pop_scope()` in the future.
    cur_block = {.begin = nullptr, .end = nullptr};  // OK to call `free_block()` on this
    cur = nullptr;

    spare_block = alloc_block(initial_block_size);
  }
  
  // Keep things simple by forbidding copying or moving
  Arena(const Arena&) = delete;
  Arena& operator=(const Arena&) = delete;
  Arena(Arena&&) = delete;
  Arena& operator=(Arena&&) = delete;
  
  ~Arena() noexcept {
    // Must deallocate all blocks
    free_block(cur_block);
    for (auto& block : full_blocks) {
      free_block(block);
    }
    free_block(spare_block);
  }
  
  auto allocate(std::size_t num_bytes) -> void* {
    // Super-fast happy path...
    DCHECK_NE(this, nullptr);
    DCHECK_LE(cur_block.begin, cur);
    DCHECK_LE(cur, cur_block.end);
    
    // Ensure `cur` stays aligned to `max_align_t`
    num_bytes = (num_bytes + alignof(max_align_t) - 1) & -alignof(max_align_t);
    auto result = cur;
    cur += num_bytes;
    if (cur > cur_block.end) {
      // Slow unhappy path (in a different function to allow inlining of `allocate` in the
      // common case)
      result = slow_allocate(num_bytes);
    }
    
    DCHECK_LE(cur_block.begin, result);
    DCHECK_LE(result + num_bytes, cur_block.end);
    return result;
  }

  auto free(void*) noexcept -> void {
    // Do nothing: call `pop_scope()` to free everything allocated since last `push_scope()`
  }

  auto push_scope() -> void {
    DCHECK_NE(this, nullptr);
    marks.push(cur);
  }
  
  auto pop_scope() noexcept -> void {
    // Super-fast happy path if topmost `mark` is within `cur_block`
    DCHECK_NE(this, nullptr);
    CHECK(not marks.empty());
    if (cur_block.contains(marks.top())) {
      cur = marks.top();
      marks.pop();
    } else {
      slow_pop_scope();
    }
  }

 private:
  // Defined in .cpp to prevent inlining
  auto slow_allocate(std::ptrdiff_t num_bytes) -> char*;
  auto slow_pop_scope() noexcept -> void;
};


// Use an `Arena_scope` is a RAII class to push and pop a scope
// in a particular `Arena`
class Arena_scope {
 private:
  Arena* arena;
  
 public:
  Arena_scope(Arena& arena) : arena{&arena} { this->arena->push_scope(); }
  ~Arena_scope() noexcept { this->arena->pop_scope(); }
};


// Each thread has implicit access to a `local_arena`.  This is particularly
// useful for transparently allowing STL containers to use `Arena`s.
//
// If we simply exposed an `Arena` object directly as a global, then every
// reference to `local_arena` would involve a procedure call to check if the
// variable had been initialized before, and if not, create it and add it to
// a list of variables to destroy at thread exit.
//
// Instead, we store a pointer to the allocated arena, and only refer to the
// underlying object when setting up a `Local_arena_scope`.  This ensures that
// access to `local_arena` takes 1 instruction in the common case, at the cost
// of requiring from users that at least one `Local_arena_scope` be active
// before using `local_arena`.
//
extern thread_local Arena real_local_arena;
extern thread_local Arena* local_arena;


// A `Local_arena_scope` is a RAII class to push and pop a scope of the `local_arena`
// (unlike `Arena_scope`, it occupiues no memory).
class Local_arena_scope {
 public:
  Local_arena_scope() {
    if (local_arena == nullptr) {
      // Actually initialize thread-local arena here!  This reference calls the
      // `Arena` constructor if necessary, and arranges for the destructor to be
      // called on thread exit.
      local_arena = &real_local_arena;
    }
    local_arena->push_scope();
  }
  ~Local_arena_scope() noexcept { local_arena->pop_scope(); }
};


// `Stl_local_arena_allocator` is a STL-compatible allocator that implicitly
// sources all memory from the `local_arena`.  Combined with a `Local_arena_scope`,
// this allows hot loops to use STL containers while rarely incurring any
// memory allocation overhead.
//
// Adapted from https://docs.microsoft.com/en-us/cpp/standard-library/allocators?view=msvc-170
// 
// [Note: below, `T = int` by default for ease of use when allocator will just
//  be rebound to a different type]
template<typename T = int>
class Stl_local_arena_allocator {
 public:
  using value_type = T;

  Stl_local_arena_allocator() {}
  template<typename U>
  explicit Stl_local_arena_allocator(const Stl_local_arena_allocator<U>& /*rhs*/) noexcept {}

  auto allocate(std::size_t n) -> T* {
    return static_cast<T*>(local_arena->allocate(n * sizeof(T)));
  }
  auto deallocate(T* p, std::size_t /*n*/) noexcept -> void {
    local_arena->free(p);
  }

  auto operator==(const Stl_local_arena_allocator<T>& rhs) const noexcept -> bool = default;
};

// Aliases for common STL containers that source memory from `local_arena`
template<typename T>
using Scratch_vector = std::vector<T, Stl_local_arena_allocator<T>>;

template<typename T>
using Scratch_stack = std::stack<T, Scratch_vector<T>>;

template<typename T>
using Scratch_deque = std::deque<T, Stl_local_arena_allocator<T>>;

template<typename K,
         typename V,
         typename Hash = typename absl::flat_hash_map<K, V>::hasher,
         typename Eq = typename absl::flat_hash_map<K, V, Hash>::key_equal>
using Scratch_flat_hash_map = absl::flat_hash_map<K, V, Hash, Eq, Stl_local_arena_allocator<std::pair<const K, V>>>;

template<typename T,
         typename Hash = typename absl::flat_hash_set<T>::hasher,
         typename Eq = typename absl::flat_hash_set<T, Hash>::key_equal>
using Scratch_flat_hash_set = absl::flat_hash_set<T, Hash, Eq, Stl_local_arena_allocator<T>>;

}  // namespace delphy

#endif // DELPHY_SCRATCH_SPACE_H_
