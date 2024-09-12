#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include <cmath>
#include <random>

#include "very_scalable_coalescent.h"

namespace delphy {

TEST(Very_scalable_coalescent_test, cells) {
  const auto t_ref = 8.0;
  const auto t_step = 2.0;

  using namespace Very_scalable_coalescent_prior;
  
  EXPECT_THAT(cell_for(7.5, t_ref, t_step), testing::Eq(0));
  EXPECT_THAT(cell_for(6.1, t_ref, t_step), testing::Eq(0));
  EXPECT_THAT(cell_for(5.9, t_ref, t_step), testing::Eq(1));
  EXPECT_THAT(cell_for(2.5, t_ref, t_step), testing::Eq(2));
  EXPECT_THAT(cell_for(3.5, t_ref, t_step), testing::Eq(2));

  EXPECT_THAT(cell_lbound(0, t_ref, t_step), testing::DoubleNear(6.0, 1e-6));
  EXPECT_THAT(cell_ubound(0, t_ref, t_step), testing::DoubleNear(8.0, 1e-6));
  EXPECT_THAT(cell_lbound(1, t_ref, t_step), testing::DoubleNear(4.0, 1e-6));
  EXPECT_THAT(cell_ubound(1, t_ref, t_step), testing::DoubleNear(6.0, 1e-6));
}

TEST(Very_scalable_coalescent_test, add_interval_single_whole_cell) {
  const auto t_ref = 8.0;
  const auto t_step = 2.0;

  using namespace Very_scalable_coalescent_prior;

  auto cells = std::vector<double>(4, 0.0);

  //                                     (6,8] (4,6] (2,4] (0,2]
  EXPECT_THAT(cells, testing::ElementsAre(0.0,  0.0,  0.0,  0.0));
  add_interval(2.0, 4.0, +2.0, cells, t_ref, t_step);
  EXPECT_THAT(cells, testing::ElementsAre(0.0,  0.0,  2.0,  0.0));
}

TEST(Very_scalable_coalescent_test, add_interval_two_whole_cells) {
  const auto t_ref = 8.0;
  const auto t_step = 2.0;

  using namespace Very_scalable_coalescent_prior;

  auto cells = std::vector<double>(4, 0.0);

  //                                     (6,8] (4,6] (2,4] (0,2]
  EXPECT_THAT(cells, testing::ElementsAre(0.0,  0.0,  0.0,  0.0));
  add_interval(2.0, 6.0, +2.0, cells, t_ref, t_step);
  EXPECT_THAT(cells, testing::ElementsAre(0.0,  2.0,  2.0,  0.0));
}

TEST(Very_scalable_coalescent_test, add_interval_single_partial_cell) {
  const auto t_ref = 8.0;
  const auto t_step = 2.0;

  using namespace Very_scalable_coalescent_prior;

  auto cells = std::vector<double>(4, 0.0);

  //                                     (6,8] (4,6] (2,4] (0,2]
  EXPECT_THAT(cells, testing::ElementsAre(0.0,  0.0,  0.0,  0.0));
  add_interval(3.0, 3.5, +2.0, cells, t_ref, t_step);
  EXPECT_THAT(cells, testing::ElementsAre(0.0,  0.0,  0.5, 0.0));
}

TEST(Very_scalable_coalescent_test, add_interval_two_partial_cells) {
  const auto t_ref = 8.0;
  const auto t_step = 2.0;

  using namespace Very_scalable_coalescent_prior;

  auto cells = std::vector<double>(4, 0.0);

  //                                     (6,8] (4,6] (2,4] (0,2]
  EXPECT_THAT(cells, testing::ElementsAre(0.0,  0.0,  0.0,  0.0));
  add_interval(3.5, 5.0, +2.0, cells, t_ref, t_step);
  EXPECT_THAT(cells, testing::ElementsAre(0.0,  1.0,  0.5,  0.0));
}

TEST(Very_scalable_coalescent_test, add_interval_two_partial_cells_and_one_whole_cell) {
  const auto t_ref = 8.0;
  const auto t_step = 2.0;

  using namespace Very_scalable_coalescent_prior;

  auto cells = std::vector<double>(4, 0.0);

  //                                     (6,8] (4,6] (2,4] (0,2]
  EXPECT_THAT(cells, testing::ElementsAre(0.0,  0.0,  0.0,  0.0));
  add_interval(3.5, 7.0, +2.0, cells, t_ref, t_step);
  EXPECT_THAT(cells, testing::ElementsAre(1.0,  2.0,  0.5,  0.0));
}

}  // namespace delphy
