#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "ancestral_tree_prober.h"

#include "estd.h"
#include "generic_tree_prober.h"

namespace delphy {

inline constexpr auto rA = Real_seq_letter::A;
inline constexpr auto rC = Real_seq_letter::C;
inline constexpr auto rG = Real_seq_letter::G;
inline constexpr auto rT = Real_seq_letter::T;

class Ancestral_tree_prober_test : public testing::Test {
 protected:
  //
  // The tree that we build:
  //
  // Time:    -1.0        0.0        1.0        2.0        3.0
  //          
  //                       +-- T0C -- a
  //                       |
  //            +-- A0T ---+ x
  //            |          |
  //          r +          +-------- A1G ------- b
  //            |
  //            +------------------- A0G ------------------ c
  // 

  Real_sequence ref_sequence{rA, rA};
  Phylo_tree tree{5};

  static constexpr Node_index r = 0;
  static constexpr Node_index x = 1;
  static constexpr Node_index a = 2;
  static constexpr Node_index b = 3;
  static constexpr Node_index c = 4;

  double const_pop{0.2};
  Const_pop_model const_pop_model{const_pop};

  double exp_pop_n0{1.0};
  double exp_pop_g{0.1};
  Exp_pop_model exp_pop_model{0.0, exp_pop_n0, exp_pop_g, 0.0};

  std::map<std::string, std::reference_wrapper<Pop_model>> pop_models{
    {"Constant", const_pop_model},
    {"Exponential", exp_pop_model}
  };

  Ancestral_tree_prober_test() {
    tree.root = r;
    tree.ref_sequence = ref_sequence;

    tree.at(r).parent = k_no_node;
    tree.at(r).children = {x, c};
    tree.at(r).name = "r";
    tree.at(r).t_min = -std::numeric_limits<float>::max();
    tree.at(r).t_max = +std::numeric_limits<float>::max();
    tree.at(r).t = -1.0;

    tree.at(x).parent = r;
    tree.at(x).children = {a, b};
    tree.at(x).name = "x";
    tree.at(x).t_min = -std::numeric_limits<float>::max();
    tree.at(x).t_max = +std::numeric_limits<float>::max();
    tree.at(x).t = 0.0;
    tree.at(x).mutations = {Mutation{rA, 0, rT, 0.5}};

    tree.at(a).parent = x;
    tree.at(a).children = {};
    tree.at(a).name = "a";
    tree.at(a).t = tree.at(a).t_min = tree.at(a).t_max = 1.0;
    tree.at(a).mutations = {Mutation{rT, 0, rC, 0.5}};

    tree.at(b).parent = x;
    tree.at(b).children = {};
    tree.at(b).name = "b";
    tree.at(b).t = tree.at(b).t_min = tree.at(b).t_max = 2.0;
    tree.at(b).mutations = {Mutation{rA, 1, rG, 1.0}};

    tree.at(c).parent = r;
    tree.at(c).children = {};
    tree.at(c).name = "c";
    tree.at(c).t = tree.at(c).t_min = tree.at(c).t_max = 3.0;
    tree.at(c).mutations = {Mutation{rA, 0, rG, 1.0}};
  }
};

TEST_F(Ancestral_tree_prober_test, invalid_marked_node) {
  auto marked_ancestors = std::vector<Node_index>{static_cast<Node_index>(tree.size()+10)};
  auto num_t_cells = 20;
  auto probe_times = make_uniform_probe_times(-1.0, 3.0, num_t_cells);
  auto raw_results = std::vector<double>((std::ssize(marked_ancestors)+1) * num_t_cells, std::numeric_limits<double>::quiet_NaN());
  EXPECT_THROW((probe_ancestors_on_tree(tree, const_pop_model, marked_ancestors, probe_times, raw_results)),
               std::out_of_range);
}


TEST_F(Ancestral_tree_prober_test, invalid_timelines) {
  auto marked_ancestors = std::vector<Node_index>{r};
  auto probe_times = {-3.5, -4.0, -4.5};
  auto num_t_cells = std::ssize(probe_times);
  auto raw_results = std::vector<double>((std::ssize(marked_ancestors)+1) * num_t_cells, std::numeric_limits<double>::quiet_NaN());
  EXPECT_THROW((probe_ancestors_on_tree(tree, const_pop_model, marked_ancestors, probe_times, raw_results)),
               std::invalid_argument);
}

TEST_F(Ancestral_tree_prober_test, empty) {
  auto t_start = -1.2;
  auto t_end = 3.0;
  auto t_step = 0.2;
  auto num_t_cells = static_cast<int>(std::round((t_end - t_start) / t_step));
  auto probe_times = make_uniform_probe_times(t_start, t_end, num_t_cells);

  auto marked_ancestors = std::vector<Node_index>{};
  auto num_rows = static_cast<int>(std::ssize(marked_ancestors)) + 1;
  auto raw_results = std::vector<double>(num_rows * num_t_cells, std::numeric_limits<double>::quiet_NaN());
  probe_ancestors_on_tree(tree, const_pop_model, marked_ancestors, probe_times, raw_results);
  auto results = estd::View_2d{raw_results, num_rows, num_t_cells};

  // If no nodes are marked, then no marked ancestors are ever ancestral to the probe
  EXPECT_THAT(results(0), testing::Each(testing::Ge(1.0 - 1e-6)));
}

TEST_F(Ancestral_tree_prober_test, trivial) {
  auto probe_times = {-1.9};
  auto num_t_cells = static_cast<int>(std::ssize(probe_times));

  ASSERT_THAT(*(probe_times.end() - 1), testing::Lt(tree.at(r).t));

  auto marked_ancestors = std::vector<Node_index>{r, a};
  auto num_rows = static_cast<int>(std::ssize(marked_ancestors)) + 1;
  auto raw_results = std::vector<double>(num_rows * num_t_cells, std::numeric_limits<double>::quiet_NaN());
  probe_ancestors_on_tree(tree, const_pop_model, marked_ancestors, probe_times, raw_results);
  auto results = estd::View_2d{raw_results, num_rows, num_t_cells};

  // If the probe sample is taken before the root's time, no marked ancestor can be ancestral to it
  EXPECT_THAT(results(0, 0), testing::DoubleNear(0.0, 1e-6));
  EXPECT_THAT(results(1, 0), testing::DoubleNear(0.0, 1e-6));
  EXPECT_THAT(results(2, 0), testing::DoubleNear(1.0, 1e-6));
}

TEST_F(Ancestral_tree_prober_test, typical) {
  for (const auto& [pop_model_name, pop_model] : pop_models) {
    SCOPED_TRACE(pop_model_name);
    
    // We'd like to sample every 0.1 time units from the root to the end
    auto t_start = -1.2;
    auto t_end = 3.0;
    auto t_step = 0.2;
    auto num_t_cells = static_cast<int>(std::round((t_end - t_start) / t_step));
    auto probe_times = make_uniform_probe_times(t_start, t_end, num_t_cells);
    
    // Do it!
    auto marked_ancestors = std::vector<Node_index>{x};
    auto num_rows = static_cast<int>(std::ssize(marked_ancestors)) + 1;
    auto raw_results = std::vector<double>(num_rows * num_t_cells, std::numeric_limits<double>::quiet_NaN());
    probe_ancestors_on_tree(tree, pop_model, marked_ancestors, probe_times, raw_results);
    auto results = estd::View_2d{raw_results, num_rows, num_t_cells};

    // Check that everything is sensible (including the "none" row, all probabilities add up to 1)
    for (auto cell = 0; cell != num_t_cells; ++cell) {
      auto tot_p = 0.0;
      for (auto i = 0; i != num_rows; ++i) {
        auto p = results(i, cell);
        EXPECT_THAT(p, testing::Ge(-1e-6));
        EXPECT_THAT(p, testing::Le(1+1e-6));
        tot_p += p;
      }
      EXPECT_THAT(tot_p, testing::DoubleNear(1.0, 1e-6));
    }
    
    // Visual check
    // std::cout << absl::StreamFormat("Population model: %s\n", absl::FormatStreamed(pop_model));
    // for (auto i = 0; i != num_rows; ++i) {
    //   if (i < std::ssize(marked_ancestors)) {
    //     std::cout << tree.at(marked_ancestors[i]).name << ": ";
    //   } else {
    //     std::cout << "-: ";
    //   }
    //   for (auto cell = 0; cell != num_t_cells; ++cell) {
    //     std::cout << absl::StreamFormat("%.1f, ", results(i, cell));
    //   }
    //   std::cout << "\n";
    // }
  }
}

TEST_F(Ancestral_tree_prober_test, skipped_ancestor) {
  auto t_start = -1.2;
  auto t_end = 3.0;
  auto t_step = 0.2;
  auto num_t_cells = static_cast<int>(std::round((t_end - t_start) / t_step));
  auto probe_times = make_uniform_probe_times(t_start, t_end, num_t_cells);
  
  auto marked_ancestors = std::vector<Node_index>{k_no_node, x, k_no_node};
  auto num_rows = static_cast<int>(std::ssize(marked_ancestors)) + 1;
  auto raw_results = std::vector<double>(num_rows * num_t_cells, std::numeric_limits<double>::quiet_NaN());
  
  probe_ancestors_on_tree(tree, const_pop_model, marked_ancestors, probe_times, raw_results);
  
  auto results = estd::View_2d{raw_results, num_rows, num_t_cells};

  EXPECT_THAT(results(0), testing::Each(testing::Eq(0.0)));
  EXPECT_THAT(results(2), testing::Each(testing::Eq(0.0)));
}

TEST_F(Ancestral_tree_prober_test, t_start_after_root) {
  // Probing from after the root should give the same results as the corresponding tail of
  // a probe that starts before the root
  auto t_step = 0.2;
  auto t_end = 3.0;
  auto marked_ancestors = std::vector<Node_index>{x, a};
  auto num_rows = static_cast<int>(std::ssize(marked_ancestors)) + 1;

  auto full_t_start = -1.2;
  auto full_num_t_cells = static_cast<int>(std::round((t_end - full_t_start) / t_step));
  auto raw_full_results = std::vector<double>(num_rows * full_num_t_cells, std::numeric_limits<double>::quiet_NaN());
  auto full_probe_times = make_uniform_probe_times(full_t_start, t_end, full_num_t_cells);
  probe_ancestors_on_tree(tree, const_pop_model, marked_ancestors, full_probe_times, raw_full_results);
  auto full_results = estd::View_2d{raw_full_results, num_rows, full_num_t_cells};

  auto late_t_start = 0.0;
  auto late_num_t_cells = static_cast<int>(std::round((t_end - late_t_start) / t_step));
  ASSERT_THAT(late_t_start, testing::Gt(tree.at(r).t));
  auto late_probe_times = make_uniform_probe_times(late_t_start, t_end, late_num_t_cells);
  auto raw_late_results = std::vector<double>(num_rows * late_num_t_cells, std::numeric_limits<double>::quiet_NaN());
  probe_ancestors_on_tree(tree, const_pop_model, marked_ancestors, late_probe_times, raw_late_results);
  auto late_results = estd::View_2d{raw_late_results, num_rows, late_num_t_cells};

  auto offset = full_num_t_cells - late_num_t_cells;
  for (auto i = 0; i != num_rows; ++i) {
    for (auto cell = 0; cell != late_num_t_cells; ++cell) {
      EXPECT_THAT(late_results(i, cell), testing::DoubleNear(full_results(i, cell + offset), 1e-6))
          << "i = " << i << ", cell = " << cell;
    }
  }
}

}  // namespace delphy
