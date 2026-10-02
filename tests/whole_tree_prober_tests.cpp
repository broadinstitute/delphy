#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "whole_tree_prober.h"
#include "estd.h"

namespace delphy {

inline constexpr auto rA = Real_seq_letter::A;
inline constexpr auto rC = Real_seq_letter::C;
inline constexpr auto rG = Real_seq_letter::G;
inline constexpr auto rT = Real_seq_letter::T;

class Whole_tree_prober_test : public testing::Test {
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

  Whole_tree_prober_test() {
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

TEST_F(Whole_tree_prober_test, invalid_timelines) {
  auto probe_times = {-3.5, -4.0, -4.5};
  auto raw_results = std::vector<double>(tree.size() * std::ssize(probe_times), std::numeric_limits<double>::quiet_NaN());
  
  EXPECT_THROW((probe_whole_tree(tree, const_pop_model, probe_times, false, raw_results)),
               std::invalid_argument);
}

TEST_F(Whole_tree_prober_test, trivial) {
  auto probe_times = {-2.0};
  auto num_probe_times = static_cast<int>(std::ssize(probe_times));

  ASSERT_THAT(*(probe_times.end() - 1), testing::Lt(tree.at(r).t));

  auto raw_results = std::vector<double>(tree.size() * std::ssize(probe_times), std::numeric_limits<double>::quiet_NaN());
  
  probe_whole_tree(tree, const_pop_model, probe_times, false, raw_results);
  auto results = estd::View_2d{raw_results, tree.size(), num_probe_times};

  // If the probe sample is taken before the root's time, no marked ancestor can be ancestral to it
  EXPECT_THAT(results(r, 0), testing::DoubleNear(1.0, 1e-6));
  EXPECT_THAT(results(a, 0), testing::DoubleNear(0.0, 1e-6));
  EXPECT_THAT(results(b, 0), testing::DoubleNear(0.0, 1e-6));
  EXPECT_THAT(results(c, 0), testing::DoubleNear(0.0, 1e-6));
  EXPECT_THAT(results(x, 0), testing::DoubleNear(0.0, 1e-6));
}

TEST_F(Whole_tree_prober_test, typical) {
  for (const auto& [pop_model_name, pop_model] : pop_models) {
    SCOPED_TRACE(pop_model_name);
    
    // We'd like to sample every 0.1 time units from the root to the end
    auto probe_times = std::vector<double>{};
    for (auto i = -12; i <= 30; i += 2) {
      probe_times.push_back(i / 10.0);
    }
    auto num_probe_times = static_cast<int>(std::ssize(probe_times));
    
    // Do it!
    auto raw_results = std::vector<double>(tree.size() * std::ssize(probe_times), std::numeric_limits<double>::quiet_NaN());
    probe_whole_tree(tree, const_pop_model, probe_times, false, raw_results);
    auto results = estd::View_2d{raw_results, tree.size(), num_probe_times};

    // Check that everything is sensible (including the "root" row, all probabilities add up to 1)
    for (auto i = 0; i != num_probe_times; ++i) {
      auto tot_p = 0.0;
      for (const auto& node : index_order_traversal(tree)) {
        auto p = results(node, i);
        EXPECT_THAT(p, testing::Ge(-1e-6));
        EXPECT_THAT(p, testing::Le(1+1e-6));
        tot_p += p;
      }
      EXPECT_THAT(tot_p, testing::DoubleNear(1.0, 1e-6));
    }
    
    // Visual check
    // std::cout << absl::StreamFormat("Population model: %s\n", absl::FormatStreamed(pop_model));
    // for (const auto& node : index_order_traversal(tree)) {
    //   std::cout << tree.at(node).name << ": ";
    //   for (auto i = 0; i != num_probe_times; ++i) {
    //     std::cout << absl::StreamFormat("%.1f, ", results(node, i));
    //   }
    //   std::cout << "\n";
    // }
  }
}

TEST_F(Whole_tree_prober_test, typical_indirect_descendants) {
  for (const auto& [pop_model_name, pop_model] : pop_models) {
    SCOPED_TRACE(pop_model_name);
    
    // We'd like to sample every 0.1 time units from the root to the end
    auto probe_times = std::vector<double>{};
    for (auto i = -12; i <= 30; i += 2) {
      probe_times.push_back(i / 10.0);
    }
    auto num_probe_times = static_cast<int>(std::ssize(probe_times));
    
    // Do it!
    auto raw_results = std::vector<double>(tree.size() * std::ssize(probe_times), std::numeric_limits<double>::quiet_NaN());
    probe_whole_tree(tree, const_pop_model, probe_times, true, raw_results);
    auto results = estd::View_2d{raw_results, tree.size(), num_probe_times};

    // Check that everything is sensible
    for (auto i = 0; i != num_probe_times; ++i) {
      for (const auto& node : index_order_traversal(tree)) {
        auto p = results(node, i);
        EXPECT_THAT(p, testing::Ge(-1e-6));
        EXPECT_THAT(p, testing::Le(1+1e-6));
      }
    }

    // Check that every branch's coalescence probability exceeds that of its children
    for (const auto& node : index_order_traversal(tree)) {
      for (auto i = 0; i != num_probe_times; ++i) {
        auto p_parent = results(node, i);
        auto p_children = 0.0;
        for (const auto& child : tree.at(node).children) {
          p_children += results(child, i);
        }
        EXPECT_THAT(p_parent, testing::Ge(p_children));
      }
    }
    
    // Visual check
    // std::cout << absl::StreamFormat("Population model: %s\n", absl::FormatStreamed(pop_model));
    // for (const auto& node : index_order_traversal(tree)) {
    //   std::cout << tree.at(node).name << ": ";
    //   for (auto i = 0; i != num_probe_times; ++i) {
    //     std::cout << absl::StreamFormat("%.1f, ", results(node, i));
    //   }
    //   std::cout << "\n";
    // }
  }
}

}  // namespace delphy
