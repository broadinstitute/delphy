#include <gtest/gtest.h>

#include <absl/log/initialize.h>

#include "scratch_space.h"

class LocalArenaScopePerTest : public testing::EmptyTestEventListener {
  delphy::Local_arena_scope* scope = nullptr;
  
  void OnTestStart(const testing::TestInfo&) override {
    scope = new delphy::Local_arena_scope{};
  }
  void OnTestEnd(const testing::TestInfo&) override {
    delete scope;
  }
};

auto main(int argc, char **argv) -> int {
  absl::InitializeLog();
  testing::InitGoogleTest(&argc, argv);
  testing::TestEventListeners& listeners = testing::UnitTest::GetInstance()->listeners();
  listeners.Append(new LocalArenaScopePerTest);
  return RUN_ALL_TESTS();
}
