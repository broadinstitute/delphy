#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "interval_set.h"

#include "mutations.h"

#include "absl/hash/hash_testing.h"

namespace delphy {

inline constexpr auto rA = Real_seq_letter::A;
inline constexpr auto rC = Real_seq_letter::C;
inline constexpr auto rG = Real_seq_letter::G;
inline constexpr auto rT = Real_seq_letter::T;

TEST(Missation_map_test, empty) {
  auto mis = Scratch_missation_map();

  EXPECT_THAT(mis.empty(), testing::IsTrue());
}

TEST(Missation_map_test, simple) {
  auto ref_seq = Real_sequence{rA, rC, rG, rT};
  auto mis = Scratch_missation_map({Missation{0, rA}, Missation{1, rG}}, ref_seq);

  EXPECT_THAT(mis.intervals, testing::ElementsAre(
      Site_interval{0, 2}));
  EXPECT_THAT(mis.from_states, testing::ElementsAre(
      //Missation{0, rA},   // implicit from ref_seq
      Missation{1, rG}));
}

TEST(Missation_map_test, missation_printing) {
  auto mi = Missation{4, rA};
  auto ss = std::stringstream{};
  ss << mi;
  EXPECT_THAT(ss.str(), testing::StrEq("A4N"));
}

TEST(Missation_map_test, missation_map_basics) {
  auto ref_seq = Real_sequence{rA, rA, rA, rA, rA, rA, rA};
  auto mi_map = Scratch_missation_map{
    {Missation{4, rA}, Missation{6, rC}, Missation{2, rG}, Missation{4, rA}},
    ref_seq
  };

  EXPECT_THAT(mi_map.num_intervals(), testing::Eq(3));
  EXPECT_THAT(mi_map.num_sites(), testing::Eq(3));
  EXPECT_THAT(mi_map.intervals, testing::ElementsAre(
      Site_interval{2, 3},
      Site_interval{4, 5},
      Site_interval{6, 7}));
  EXPECT_THAT(mi_map.from_states, testing::ElementsAre(
      Missation{2, rG},
      //Missation{4, rA},  // implicit from ref_seq
      Missation{6, rC}));

  EXPECT_THAT(mi_map.get_from_state(2, ref_seq), testing::Eq(rG));
  EXPECT_THAT(mi_map.get_from_state(4, ref_seq), testing::Eq(rA));
  EXPECT_THAT(mi_map.get_from_state(6, ref_seq), testing::Eq(rC));

  // No need for erase so far
  //mi_map.erase(6);
  //
  //EXPECT_THAT(mi_map.size(), testing::Eq(2));

  mi_map.insert(Missation{3, rT}, ref_seq);

  EXPECT_THAT(mi_map.num_intervals(), testing::Eq(2));
  EXPECT_THAT(mi_map.num_sites(), testing::Eq(4));
  EXPECT_THAT(mi_map.intervals, testing::ElementsAre(
      Site_interval{2, 5},
      Site_interval{6, 7}));
  EXPECT_THAT(mi_map.from_states, testing::ElementsAre(
      Missation{2, rG},
      Missation{3, rT},
      //Missation{4, rA},  // implicit from ref_seq
      Missation{6, rC}));
  
  EXPECT_THAT(mi_map.contains(3), testing::IsTrue());
  EXPECT_THAT(mi_map.contains(5), testing::IsFalse());

  ref_seq[2] = rG; mi_map.ref_seq_changed(2, rA, rG);
  ref_seq[3] = rC; mi_map.ref_seq_changed(3, rA, rC);
  ref_seq[4] = rT; mi_map.ref_seq_changed(4, rA, rT);
  ref_seq[5] = rC; mi_map.ref_seq_changed(5, rA, rC);
  ref_seq[6] = rC; mi_map.ref_seq_changed(6, rA, rC);
  
  EXPECT_THAT(mi_map.num_intervals(), testing::Eq(2));
  EXPECT_THAT(mi_map.num_sites(), testing::Eq(4));
  EXPECT_THAT(mi_map.intervals, testing::ElementsAre(
      Site_interval{2, 5},
      Site_interval{6, 7}));
  EXPECT_THAT(mi_map.from_states, testing::ElementsAre(
      //Missation{2, rG},      // now implicit
      Missation{3, rT},
      Missation{4, rA}         // no longer implicit
      /*Missation{6, rC}*/));  // now implicit
}

TEST(Missation_map_test, factor_out_common_missations) {
  auto ref_seq = Real_sequence{rT, rT, rT};
  auto mi_A = Missation_map<>{{Missation{0, rA}, Missation{1, rC}}, ref_seq};
  auto mi_B = Missation_map<>{{Missation{0, rA}, Missation{2, rG}}, ref_seq};
  auto result_common = Missation_map<>{};

  factor_out_common_missations(mi_A, mi_B, result_common);

  EXPECT_THAT(estd::ranges::to_vec(mi_A.slow_elements(ref_seq)), testing::ElementsAre(Missation{1, rC}));
  EXPECT_THAT(estd::ranges::to_vec(mi_B.slow_elements(ref_seq)), testing::ElementsAre(Missation{2, rG}));
  EXPECT_THAT(estd::ranges::to_vec(result_common.slow_elements(ref_seq)), testing::ElementsAre(Missation{0, rA}));
}

TEST(Missation_map_test, merge_missations_nondestructively) {
  auto ref_seq = Real_sequence{rA, rA, rA};
  const auto mi_src_A = Missation_map<>{{Missation{0, rA}, Missation{1, rC}}, ref_seq};
  const auto mi_src_B = Missation_map<>{{Missation{2, rG}}, ref_seq};

  auto result = merge_missations_nondestructively(mi_src_A, mi_src_B);

  EXPECT_THAT(estd::ranges::to_vec(result.slow_elements(ref_seq)), testing::ElementsAre(
      Missation{0, rA}, Missation{1, rC}, Missation{2, rG}));
}

TEST(Missation_map_test, subtract_missations_nondestructively) {
  auto ref_seq = Real_sequence{rA, rA, rA};
  const auto mi_src_A = Missation_map<>{{Missation{0, rA}, Missation{1, rC}}, ref_seq};
  const auto mi_src_B = Missation_map<>{{Missation{1, rC}, Missation{2, rG}}, ref_seq};
  auto result = Missation_map<>{};
  
  subtract_missations_nondestructively(result, mi_src_A, mi_src_B);

  EXPECT_THAT(estd::ranges::to_vec(result.slow_elements(ref_seq)), testing::ElementsAre(
      Missation{0, rA}));
}

}  // namespace delphy
