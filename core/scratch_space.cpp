#include "scratch_space.h"

#include <chrono>
#include <thread>
#include <iostream>

namespace delphy {

thread_local Arena real_local_arena{};
thread_local Arena* local_arena = nullptr;  // See comments in scratch2.h

auto Arena::slow_allocate(std::ptrdiff_t num_bytes) -> char* {
  // Current block is full
  CHECK_GT(cur + num_bytes, cur_block.end);
  full_blocks.push_back(cur_block);
  
  // Reuse spare block or allocate brand new block?
  if (std::ssize(spare_block) >= num_bytes) {
    cur_block = spare_block;
    spare_block = {};
    
  } else {
    auto new_block_size = std::max(
        num_bytes,
        2 * std::max(std::ssize(cur_block), std::ssize(spare_block)));
    cur_block = alloc_block(new_block_size);
  }

  // `cur_block` is now completely unused and big enough
  cur = cur_block.begin;

  CHECK_LE(cur_block.begin, cur);
  CHECK_LE(cur + num_bytes, cur_block.end);
  
  auto result = cur;
  cur += num_bytes;
  return result;
}

auto Arena::slow_pop_scope() noexcept -> void {
  CHECK(not marks.empty());
  cur = marks.top();
  marks.pop();
  
  // `cur` is not inside `cur_block`, but it's in block `full_blocks[i]` for some `i`,
  // and everything in `full_blocks[i+1]`, ... `full_blocks.back()` can be burned.  We
  // hang on to the largest block instead of deleting it immediately in case we need the
  // scratch space to grow again.
  
  while (not cur_block.contains(cur)) {
    // "Free" `cur_block`; might trigger release of `spare_block` or `cur_block`
    if (std::ssize(cur_block) < std::ssize(spare_block)) {
      free_block(cur_block);
    } else {
      // Discard (small) `spare_block` and make (large) `cur_block` the spare block
      free_block(spare_block);
      spare_block = cur_block;
    }
    
    CHECK(not full_blocks.empty());
    cur_block = full_blocks.back();
    full_blocks.pop_back();
  }
}

}  // namespace delphy
