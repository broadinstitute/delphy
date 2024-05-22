#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <absl/log/check.h>

#include "sequence_overlay.h"

#include "sequence.h"

namespace delphy {

inline constexpr auto rA = Real_seq_letter::A;
inline constexpr auto rC = Real_seq_letter::C;
inline constexpr auto rG = Real_seq_letter::G;
inline constexpr auto rT = Real_seq_letter::T;

TEST(Sequence_overlay_test, trivial) {
  const auto base_seq = Real_sequence{};
  auto overlay = Sequence_overlay{base_seq};

  EXPECT_THAT(overlay, testing::IsEmpty());
  EXPECT_THAT(&overlay.base_sequence(), testing::Eq(&base_seq));
  EXPECT_THAT(std::ssize(overlay), testing::Eq(std::ssize(base_seq)));
  EXPECT_THAT(overlay.deltas(), testing::SizeIs(0));
}

TEST(Sequence_overlay_test, no_changes_direct_const_access) {
  const auto base_seq = Real_sequence{rA, rC, rG, rT};
  const auto overlay = Sequence_overlay{base_seq};

  EXPECT_THAT(overlay[0], testing::Eq(rA));
  EXPECT_THAT(overlay[1], testing::Eq(rC));
  EXPECT_THAT(overlay[2], testing::Eq(rG));
  EXPECT_THAT(overlay[3], testing::Eq(rT));
  
  EXPECT_THAT(overlay.deltas(), testing::SizeIs(0));
}

TEST(Sequence_overlay_test, no_changes_direct_non_const_access) {
  const auto base_seq = Real_sequence{rA, rC, rG, rT};
  auto overlay = Sequence_overlay{base_seq};

  EXPECT_THAT(overlay[0], testing::Eq(rA));
  EXPECT_THAT(overlay[1], testing::Eq(rC));
  EXPECT_THAT(overlay[2], testing::Eq(rG));
  EXPECT_THAT(overlay[3], testing::Eq(rT));

  EXPECT_THAT(overlay.deltas(), testing::SizeIs(0));
}

TEST(Sequence_overlay_test, single_change_direct_access) {
  const auto base_seq = Real_sequence{rA, rA};
  auto overlay = Sequence_overlay{base_seq};

  overlay[1] = rC;

  EXPECT_THAT(overlay[0], testing::Eq(rA));
  EXPECT_THAT(overlay[1], testing::Eq(rC));
  EXPECT_THAT(overlay.deltas(), testing::SizeIs(1));
}

TEST(Sequence_overlay_test, single_reverted_change_leads_to_no_deltas) {
  const auto base_seq = Real_sequence{rA, rA};
  auto overlay = Sequence_overlay{base_seq};

  overlay[1] = rC;
  overlay[1] = base_seq[1];  // Revert

  EXPECT_THAT(overlay[0], testing::Eq(rA));
  EXPECT_THAT(overlay[1], testing::Eq(rA));
  EXPECT_THAT(overlay.deltas(), testing::SizeIs(0));
}

TEST(Sequence_overlay_test, no_changes_iterator_access) {
  const auto base_seq = Real_sequence{rA, rC, rG, rT};
  const auto overlay = Sequence_overlay{base_seq};

  EXPECT_THAT(overlay, testing::ElementsAre(rA, rC, rG, rT));
}

TEST(Sequence_overlay_test, const_iterator_madness) {
  const auto base_seq = Real_sequence{rA, rC, rG, rT};
  const auto overlay = Sequence_overlay{base_seq};

  auto it = overlay.begin();
  auto it2 = it;
  EXPECT_THAT(it, testing::Eq(it2));
  auto it3 = it2++;
  EXPECT_THAT(it, testing::Lt(it2));
  EXPECT_THAT(it3, testing::Eq(it));
  --it2;
  EXPECT_THAT(it2, testing::Eq(it));

  auto it_to_0 = overlay.begin();
  auto it_to_1 = 1 + it_to_0;
  auto it_to_2 = it_to_0;  it_to_2 += 2;
  auto it_to_3 = overlay.end();  // almost there...
  auto it_to_end = it_to_3--;
  
  EXPECT_THAT(*it_to_0, testing::Eq(rA));
  EXPECT_THAT(*it_to_1, testing::Eq(rC));
  EXPECT_THAT(*it_to_2, testing::Eq(rG));
  EXPECT_THAT(*it_to_3, testing::Eq(rT));
  EXPECT_THAT(it_to_end, testing::Eq(overlay.end()));

  EXPECT_THAT(it_to_3 - it_to_1, testing::Eq(2));

  auto itt = it_to_0;
  ++itt;  EXPECT_THAT(itt, testing::Eq(it_to_1));
  --itt;  EXPECT_THAT(itt, testing::Eq(it_to_0));
  itt = it_to_0 + 3;
  itt -= 2;
  EXPECT_THAT(itt, testing::Eq(it_to_1));

  EXPECT_THAT(it_to_1[2], testing::Eq(base_seq[3]));
  EXPECT_THAT(it_to_2 - 1, testing::Eq(it_to_1));

  EXPECT_THAT(overlay, testing::ElementsAre(rA, rC, rG, rT));
}

TEST(Sequence_overlay_test, iterator_madness) {
  const auto base_seq = Real_sequence{rA, rC, rG, rT};
  auto overlay = Sequence_overlay{base_seq};

  auto it = overlay.begin();
  auto it2 = it;
  EXPECT_THAT(it, testing::Eq(it2));
  auto it3 = it2++;
  EXPECT_THAT(it, testing::Lt(it2));
  EXPECT_THAT(it3, testing::Eq(it));
  --it2;
  EXPECT_THAT(it2, testing::Eq(it));

  auto it_to_0 = overlay.begin();
  auto it_to_1 = 1 + it_to_0;
  auto it_to_2 = it_to_0;  it_to_2 += 2;
  auto it_to_3 = overlay.end();  // almost there...
  auto it_to_end = it_to_3--;
  
  EXPECT_THAT(*it_to_0, testing::Eq(rA));
  EXPECT_THAT(*it_to_1, testing::Eq(rC));
  EXPECT_THAT(*it_to_2, testing::Eq(rG));
  EXPECT_THAT(*it_to_3, testing::Eq(rT));
  EXPECT_THAT(it_to_end, testing::Eq(overlay.end()));

  EXPECT_THAT(it_to_3 - it_to_1, testing::Eq(2));

  auto itt = it_to_0;
  ++itt;  EXPECT_THAT(itt, testing::Eq(it_to_1));
  --itt;  EXPECT_THAT(itt, testing::Eq(it_to_0));
  itt = it_to_0 + 3;
  itt -= 2;
  EXPECT_THAT(itt, testing::Eq(it_to_1));

  EXPECT_THAT(it_to_1[2], testing::Eq(base_seq[3]));
  EXPECT_THAT(it_to_2 - 1, testing::Eq(it_to_1));

  EXPECT_THAT(overlay, testing::ElementsAre(rA, rC, rG, rT));
}

TEST(Sequence_overlay_test, algorithms_with_non_const_iterators_1) {
  const auto base_seq = Real_sequence{rC, rA, rT, rG};
  auto overlay = Sequence_overlay{base_seq};

  std::ranges::sort(overlay);

  EXPECT_THAT(base_seq, testing::ElementsAre(rC, rA, rT, rG));
  EXPECT_THAT(overlay, testing::ElementsAre(rA, rC, rG, rT));
}

TEST(Sequence_overlay_test, algorithms_with_non_const_iterators_2) {
  const auto base_seq = Real_sequence{rA, rC, rG, rT};
  auto overlay = Sequence_overlay{base_seq};

  std::ranges::rotate(overlay, overlay.begin() + 2);

  EXPECT_THAT(base_seq, testing::ElementsAre(rA, rC, rG, rT));
  EXPECT_THAT(overlay, testing::ElementsAre(rG, rT, rA, rC));
}

TEST(Sequence_overlay_test, printing) {
  const auto base_seq = Real_sequence{rA, rC, rG, rT};
  auto overlay = Sequence_overlay{base_seq};

  auto ss = std::stringstream{};
  ss << *overlay.begin();
  EXPECT_THAT(ss.str(), testing::StrEq("A"));
}

TEST(Sequence_overlay_test, absl_check_eq) {
  // Without an explicit overload of operator<< for Sequence_overlay_reference, abseil's CHECK_EQ gets confused
  const auto base_seq = Real_sequence{rA, rC, rG, rT};
  auto overlay = Sequence_overlay{base_seq};

  CHECK_EQ(*overlay.begin(), rA);
}

TEST(Sequence_overlay_test, materialize) {
  const auto base_seq = Real_sequence{rA, rC, rG, rT};
  auto overlay = Sequence_overlay{base_seq};
  overlay[0] = rT;
  overlay[3] = rA;

  auto materialized = overlay.materialize();
  
  static_assert(std::is_same_v<Real_sequence, decltype(materialized)>);
  EXPECT_THAT(materialized, testing::ElementsAre(rT, rC, rG, rA));
}

}  // namespace delphy
