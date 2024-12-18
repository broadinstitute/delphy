#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "tree_prober.h"

using namespace ::testing;

namespace delphy {

TEST(Tree_prober_test, invalid_p_initial) {
  auto pop = 0.2;   // N_e * rho
  auto pop_model = Const_pop_model{pop};

  auto num_cats = 2;
  auto t_start = 0.0;
  auto t_end = 10.0;
  auto num_cells = 10;
  
  auto branch_counts_by_cat = Staircase_family{num_cats, t_start, t_end, num_cells};

  // Improperly sized (3 values, 2 categories)
  EXPECT_THROW((Tree_prober{branch_counts_by_cat, 0, pop_model, {0.0, 0.0, 0.0}}), std::invalid_argument);

  // Negative values
  EXPECT_THROW((Tree_prober{branch_counts_by_cat, 0, pop_model, {-0.3, -0.2}}), std::out_of_range);

  // Add up to more than 1 (adding up to less is fine: it just means there's an implicit "other" category)
  EXPECT_THROW((Tree_prober{branch_counts_by_cat, 0, pop_model, {0.6, 0.7}}), std::out_of_range);
}

TEST(Tree_prober_test, trivial) {
  auto pop = 0.2;   // N_e * rho
  auto pop_model = Const_pop_model{pop};

  auto num_cats = 2;
  auto t_start = 0.0;
  auto t_end = 1e-5;
  auto num_cells = 1;
  
  auto branch_counts_by_cat = Staircase_family{num_cats, t_start, t_end, num_cells};

  auto p_initial = std::vector{0.3, 0.7};
  auto cells_to_skip = 0;
  auto tree_prober = Tree_prober{branch_counts_by_cat, cells_to_skip, pop_model, p_initial};

  // Output should match p_initial if the time domain is empty
  EXPECT_THAT(tree_prober.p(0).at(t_start), DoubleNear(0.3, 1e-6));
  EXPECT_THAT(tree_prober.p(1).at(t_start), DoubleNear(0.7, 1e-6));
}

TEST(Tree_prober_test, simple) {
  auto pop = 0.2;   // N_e * rho
  auto pop_model = Const_pop_model{pop};

  auto num_cats = 2;
  auto t_start = 0.0;
  auto t_end = 10.0;
  auto num_cells = 10;
  
  auto branch_counts_by_cat = Staircase_family{num_cats, t_start, t_end, num_cells};
  auto& branch_counts_0 = branch_counts_by_cat[0];
  auto& branch_counts_1 = branch_counts_by_cat[1];

  add_boxcar(branch_counts_0, 1.0, 5.0, 1.0);
  add_boxcar(branch_counts_1, 3.0, 10.0, 1.0);

  auto cells_to_skip = 0;
  auto tree_prober = Tree_prober{branch_counts_by_cat, cells_to_skip, pop_model};

  EXPECT_EQ(&tree_prober.branch_counts_by_category(), &branch_counts_by_cat);
  EXPECT_EQ(&tree_prober.pop_model(), &pop_model);
  EXPECT_EQ(tree_prober.num_categories(), num_cats);
  EXPECT_EQ(tree_prober.t_start(), t_start);
  EXPECT_EQ(tree_prober.t_end(), t_end);
  
  auto& p_0 = tree_prober.p(0);
  auto& p_1 = tree_prober.p(1);

  EXPECT_THAT(p_0.at(0.9), Eq(0.0));
  EXPECT_THAT(p_1.at(0.9), Eq(0.0));
  EXPECT_THAT(p_0.at(2.9), DoubleNear(1.0, 1e-2));
  EXPECT_THAT(p_1.at(2.9), Eq(0.0));
  EXPECT_THAT(p_1.at(4.9), Gt(0.4));
  EXPECT_THAT(p_0.at(4.9), Gt(p_1.at(4.9)));
  EXPECT_THAT(p_0.at(9.9), DoubleNear(0.0, 1e-2));
  EXPECT_THAT(p_1.at(9.9), DoubleNear(1.0, 1e-2));
}

TEST(Tree_prober_test, skip_one) {
  auto pop = 0.2;   // N_e * rho
  auto pop_model = Const_pop_model{pop};

  auto num_cats = 2;
  auto t_start = 0.0;
  auto t_end = 10.0;
  auto num_cells = 10;
  auto cell_size = (t_end - t_start) / num_cells;
  
  auto branch_counts_by_cat = Staircase_family{num_cats, t_start, t_end, num_cells};
  auto& branch_counts_0 = branch_counts_by_cat[0];
  auto& branch_counts_1 = branch_counts_by_cat[1];

  add_boxcar(branch_counts_0, 1.0, 5.0, 1.0);
  add_boxcar(branch_counts_1, 3.0, 10.0, 1.0);

  auto cells_to_skip = 1;
  auto tree_prober = Tree_prober{branch_counts_by_cat, cells_to_skip, pop_model};

  EXPECT_EQ(&tree_prober.branch_counts_by_category(), &branch_counts_by_cat);
  EXPECT_EQ(&tree_prober.pop_model(), &pop_model);
  EXPECT_EQ(tree_prober.num_categories(), num_cats);
  EXPECT_THAT(tree_prober.t_start(), testing::DoubleNear(t_start + cell_size * cells_to_skip,
                                                         1e-6));
  EXPECT_EQ(tree_prober.t_end(), t_end);
  
  auto& p_0 = tree_prober.p(0);
  auto& p_1 = tree_prober.p(1);

  //EXPECT_THAT(p_0.at(0.9), Eq(0.0));  // First cell skipped!
  //EXPECT_THAT(p_1.at(0.9), Eq(0.0));  // First cell skipped!
  EXPECT_THAT(p_0.at(2.9), DoubleNear(1.0, 1e-2));
  EXPECT_THAT(p_1.at(2.9), Eq(0.0));
  EXPECT_THAT(p_1.at(4.9), Gt(0.4));
  EXPECT_THAT(p_0.at(4.9), Gt(p_1.at(4.9)));
  EXPECT_THAT(p_0.at(9.9), DoubleNear(0.0, 1e-2));
  EXPECT_THAT(p_1.at(9.9), DoubleNear(1.0, 1e-2));
}

}  // namespace delphy
