#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "staircase.h"

namespace delphy {

TEST(Staircase_test, is_an_stl_container) {
  static_assert(std::is_same_v<Staircase::value_type, double>);
  // TODO: Fill this out as needed
}

TEST(Staircase_test, invalid_domain) {
  EXPECT_THROW(( Staircase{0.0, 0.0, 1} ), std::invalid_argument);
}

TEST(Staircase_test, invalid_num_cells) {
  EXPECT_THROW(( Staircase{0.0, 1.0, 0} ), std::invalid_argument);
}

TEST(Staircase_test, simple) {
  const auto staircase = Staircase{1.0, 7.0, 3};
  EXPECT_EQ(staircase.x_start(), 1.0);
  EXPECT_EQ(staircase.cell_size(), 2.0);
  EXPECT_EQ(staircase.num_cells(), 3);
  EXPECT_EQ(staircase.x_end(), 7.0);
  EXPECT_EQ(staircase.at_cell(0), 0.0);
  EXPECT_EQ(staircase.cell_for(3.5), 1);
  EXPECT_EQ(staircase.at(3.5), 0.0);
}

TEST(Staircase_test, cell_bounds) {
  const auto staircase = Staircase{0.0, 10.0, 10};
  EXPECT_EQ(staircase.cell_for(0.5), 0);
  EXPECT_EQ(cell_lbound(staircase, 0), 0.0);
  EXPECT_EQ(cell_ubound(staircase, 0), 1.0);
  EXPECT_EQ(cell_lbound(staircase, 5), 5.0);
  EXPECT_EQ(cell_ubound(staircase, 5), 6.0);
  EXPECT_EQ(cell_for_lbound(staircase, 0.0), 0);
  EXPECT_EQ(cell_for_lbound(staircase, 0.5), 0);
  EXPECT_EQ(cell_for_ubound(staircase, 0.5), 0);
  EXPECT_EQ(cell_for_ubound(staircase, 1.0), 0);
  EXPECT_EQ(cell_for_lbound(staircase, 1.0), 1);
  EXPECT_EQ(cell_for_ubound(staircase, 2.0), 1);
}

TEST(Staircase_test, add_boxcar_invalid_domain) {
  auto staircase = Staircase{0.0, 10.0, 10};
  EXPECT_THROW(( add_boxcar(staircase, 0.0, -4.0, 0.0) ), std::invalid_argument);
}

TEST(Staircase_test, add_boxcar_empty_domain) {
  auto staircase = Staircase{0.0, 10.0, 10};
  EXPECT_NO_THROW(( add_boxcar(staircase, 0.0, 0.0, 0.0) ));
}

TEST(Staircase_test, add_boxcar_overlaps_domain) {
  auto staircase = Staircase{0.0, 10.0, 10};

  add_boxcar(staircase, -1.0, 11.0, 1.0);
  add_boxcar(staircase, -1.0,  1.0, 2.0);
  add_boxcar(staircase,  9.0, 11.0, 3.0);
  
  EXPECT_EQ(staircase.at_cell(0), 3.0);
  EXPECT_EQ(staircase.at_cell(4), 1.0);
  EXPECT_EQ(staircase.at_cell(9), 4.0);
}

TEST(Staircase_test, add_boxcar_doesnt_overlap_domain) {
  auto staircase = Staircase{0.0, 10.0, 10};

  add_boxcar(staircase, -2.0, -1.0, 3.0);
  add_boxcar(staircase, 11.0, 12.0, 7.0);
  
  EXPECT_EQ(staircase.at_cell(0), 0.0);
  EXPECT_EQ(staircase.at_cell(4), 0.0);
  EXPECT_EQ(staircase.at_cell(9), 0.0);
}

TEST(Staircase_test, add_boxcar){
  auto staircase = Staircase{0.0, 10.0, 10};

  EXPECT_EQ(staircase.at_cell(0), 0.0);
  EXPECT_EQ(staircase.at_cell(1), 0.0);
  EXPECT_EQ(staircase.at_cell(2), 0.0);

  add_boxcar(staircase, 2.5, 4.25, 10.0);

  EXPECT_EQ(staircase.at_cell(2), 5.0);
  EXPECT_EQ(staircase.at_cell(3), 10.0);
  EXPECT_EQ(staircase.at_cell(4), 2.5);

  add_boxcar(staircase, 0.0, 1.0, 100.0);

  EXPECT_EQ(staircase.at_cell(0), 100.0);
  EXPECT_EQ(staircase.at_cell(1), 0.0);

  add_boxcar(staircase, 9.0, 10.0, 100.0);

  EXPECT_EQ(staircase.at_cell(8), 0.0);
  EXPECT_EQ(staircase.at_cell(9), 100.0);
}

TEST(Staircase_test, add_trapezoid_invalid_domain) {
  auto staircase = Staircase{0.0, 10.0, 10};
  EXPECT_THROW(( add_trapezoid(staircase, 0.0, -4.0, 0.0, 0.0) ), std::invalid_argument);
}

TEST(Staircase_test, add_trapezoid_empty_domain) {
  auto staircase = Staircase{0.0, 10.0, 10};
  EXPECT_NO_THROW(( add_trapezoid(staircase, 0.0, 0.0, 0.0, 0.0) ));
}

TEST(Staircase_test, add_trapezoid_overlaps_domain) {
  auto staircase = Staircase{0.0, 10.0, 10};

  add_trapezoid(staircase, -1.0, 11.0, 1*-1.0, 1*11.0);
  add_trapezoid(staircase, -1.0,  1.0, 2*-1.0, 2* 1.0);
  add_trapezoid(staircase,  9.0, 11.0, 3* 9.0, 3*11.0);
  
  EXPECT_THAT(staircase.at_cell(0), testing::DoubleNear(3*0.5, 1e-6));
  EXPECT_THAT(staircase.at_cell(4), testing::DoubleNear(1*4.5, 1e-6));
  EXPECT_THAT(staircase.at_cell(9), testing::DoubleNear(4*9.5, 1e-6));
}

TEST(Staircase_test, add_trapezoid_doesnt_overlaps_domain) {
  auto staircase = Staircase{0.0, 10.0, 10};

  add_trapezoid(staircase, -2.0, -1.0, -1.0, 11.0);
  add_trapezoid(staircase, 11.0, 12.0, -1.0, 11.0);
  
  EXPECT_THAT(staircase.at_cell(0), testing::DoubleNear(0, 1e-6));
  EXPECT_THAT(staircase.at_cell(4), testing::DoubleNear(0, 1e-6));
  EXPECT_THAT(staircase.at_cell(9), testing::DoubleNear(0, 1e-6));
}

TEST(Staircase_test, add_trapezoid){
  auto staircase = Staircase{0.0, 10.0, 10};

  EXPECT_EQ(staircase.at_cell(0), 0.0);
  EXPECT_EQ(staircase.at_cell(1), 0.0);
  EXPECT_EQ(staircase.at_cell(2), 0.0);

  add_trapezoid(staircase, 2.5, 4.25, 10.0, 20.0);

  EXPECT_THAT(staircase.at_cell(2), testing::DoubleNear(0.5 * (10.0 + (2.5 - 2.5)/(4.25 - 2.5) * (20.0 - 10.0) +
                                                               10.0 + (3.0 - 2.5)/(4.25 - 2.5) * (20.0 - 10.0)) / 2,
                                                        1e-6));
  EXPECT_THAT(staircase.at_cell(3), testing::DoubleNear(1.0 * (10.0 + (3.0 - 2.5)/(4.25 - 2.5) * (20.0 - 10.0) +
                                                               10.0 + (4.0 - 2.5)/(4.25 - 2.5) * (20.0 - 10.0)) / 2,
                                                        1e-6));
  EXPECT_THAT(staircase.at_cell(4), testing::DoubleNear(0.25 * (10.0 + (4.0 - 2.5)/(4.25 - 2.5) * (20.0 - 10.0) +
                                                                10.0 + (4.25 - 2.5)/(4.25 - 2.5) * (20.0 - 10.0)) / 2,
                                                        1e-6));

  add_trapezoid(staircase, 0.0, 1.0, 100.0, 200.0);

  EXPECT_THAT(staircase.at_cell(0), testing::DoubleNear(150.0, 1e-6));
  EXPECT_THAT(staircase.at_cell(1), testing::DoubleNear(0.0, 1e-6));

  add_trapezoid(staircase, 9.0, 10.0, 200.0, 100.0);

  EXPECT_THAT(staircase.at_cell(8), testing::DoubleNear(0.0, 1e-6));
  EXPECT_THAT(staircase.at_cell(9), testing::DoubleNear(150.0, 1e-6));
}

TEST(staircase_family, invalid_num_members) {
  EXPECT_THROW(( Staircase_family{0, 0.0, 10.0, 10} ), std::invalid_argument);
}

TEST(staircase_family, simple) {
  auto family = Staircase_family{3, 0.0, 10.0, 10};

  EXPECT_EQ(family.num_members(), 3);
  EXPECT_EQ(family.x_start(), 0.0);
  EXPECT_EQ(family.x_end(), 10.0);
  EXPECT_EQ(family.num_cells(), 10);
  
  for (const auto& member : std::as_const(family)) {
    EXPECT_EQ(member.x_start(), 0.0);
    EXPECT_EQ(member.x_end(), 10.0);
    EXPECT_EQ(member.num_cells(), 10);
  }
  
  // Check that different members have independent data
  add_boxcar(family[0], 0.0, 10.0, 1.0);
  add_boxcar(family[1], 5.0, 10.0, 2.0);

  EXPECT_EQ(family[0].at(2.5), 1.0);
  EXPECT_EQ(family[0].at(7.5), 1.0);
  EXPECT_EQ(family[1].at(2.5), 0.0);
  EXPECT_EQ(family[1].at(7.5), 2.0);

  // Check that const accessors work
  for (auto i = 0; i != family.num_members(); ++i) {
    EXPECT_EQ(&family[i], &(std::as_const(family))[i]);
  }
}

}  // namespace delphy
