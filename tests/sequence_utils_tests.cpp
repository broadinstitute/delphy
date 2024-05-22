#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "sequence_utils.h"

namespace delphy {

inline constexpr auto sA = Seq_letters::A;
inline constexpr auto sC = Seq_letters::C;
inline constexpr auto sG = Seq_letters::G;
inline constexpr auto sT = Seq_letters::T;
inline constexpr auto sN = Seq_letters::N;

inline constexpr auto rA = Real_seq_letter::A;
inline constexpr auto rC = Real_seq_letter::C;
inline constexpr auto rG = Real_seq_letter::G;
inline constexpr auto rT = Real_seq_letter::T;

TEST(Sequence_utils_test, deduce_consensus_sequence_empty) {
  const auto seqs = std::vector<Sequence>{};
  EXPECT_THROW(deduce_consensus_sequence(seqs), std::invalid_argument);
}

TEST(Sequence_utils_test, deduce_consensus_sequence_single) {
  auto seqs = std::vector<Sequence>{};
  seqs.push_back(Sequence{sA, sC, sG, sT});
  EXPECT_EQ(deduce_consensus_sequence(seqs), (Real_sequence{rA, rC, rG, rT}));
}

TEST(Sequence_utils_test, deduce_consensus_sequence_multiple) {
  auto seqs = std::vector<Sequence>{};
  seqs.push_back(Sequence{sA, sC, sT, sT});
  seqs.push_back(Sequence{sA, sA, sT, sG});
  seqs.push_back(Sequence{sC, sC, sG, sG});
  EXPECT_EQ(deduce_consensus_sequence(seqs), (Real_sequence{rA, rC, rT, rG}));
}

TEST(Sequence_utils_test, deduce_consensus_sequence_multiple_with_ambs) {
  auto seqs = std::vector<Sequence>{};
  seqs.push_back(Sequence{sA, sN, sT, sT});
  seqs.push_back(Sequence{sA, sN, sN, sN});
  seqs.push_back(Sequence{sN, sN, sG, sN});
  seqs.push_back(Sequence{sC, sC, sG, sN});
  EXPECT_EQ(deduce_consensus_sequence(seqs), (Real_sequence{rA, rC, rG, rT}));
}

TEST(Sequence_utils_test, deduce_consensus_sequence_multiple_with_unresolvable_ambs) {
  // Unresolvable ambs arbitrarily resolved to A
  auto seqs = std::vector<Sequence>{};
  seqs.push_back(Sequence{sA, sN, sN, sN});
  seqs.push_back(Sequence{sA, sN, sN, sN});
  seqs.push_back(Sequence{sN, sN, sN, sN});
  seqs.push_back(Sequence{sC, sC, sN, sN});
  EXPECT_EQ(deduce_consensus_sequence(seqs), (Real_sequence{rA, rC, rA, rA}));
}

TEST(Sequence_utils_test, deduce_consensus_sequence_view) {
  struct Entry {
    int id;
    Sequence sequence;
  };
  auto entries = std::vector{
    Entry{.id = 1, .sequence = Sequence{sA, sC, sG, sT}},
    Entry{.id = 2, .sequence = Sequence{sA, sC, sG, sT}}
  };
    
  auto consensus_sequence = deduce_consensus_sequence(
      entries | std::views::transform([](const auto& e) { return e.sequence; }));

  EXPECT_EQ(consensus_sequence, (Real_sequence{rA, rC, rG, rT}));
}

TEST(Sequence_utils_test, calculate_delta_from_reference_empty) {
  auto ref_seq = Real_sequence{};
  auto seq = Sequence{};

  const auto delta = calculate_delta_from_reference(seq, ref_seq);

  EXPECT_THAT(delta.seq_deltas, testing::IsEmpty());
  EXPECT_THAT(delta.missations, testing::IsEmpty());
}

TEST(Sequence_utils_test, calculate_delta_from_reference_invalid) {
  auto ref_seq = Real_sequence{rA};
  auto seq = Sequence{sC, sG};

  EXPECT_THROW(calculate_delta_from_reference(seq, ref_seq), std::invalid_argument);
}

TEST(Sequence_utils_test, calculate_delta_from_reference_partially_ambiguous) {
  auto ref_seq = Real_sequence{rA};
  auto seq = Sequence{sA | sC};

  const auto delta = calculate_delta_from_reference(seq, ref_seq);

  EXPECT_THAT(delta.seq_deltas, testing::IsEmpty());
  EXPECT_THAT(estd::ranges::to_vec(delta.missations.slow_elements(ref_seq)),
              testing::ElementsAre(Missation{0, rA}));
}

TEST(Sequence_utils_test, calculate_delta_from_reference_muts) {
  auto ref_seq = Real_sequence{rA, rA, rA};
  auto seq = Sequence{sA, sC, sG};

  const auto delta = calculate_delta_from_reference(seq, ref_seq);

  EXPECT_THAT(delta.seq_deltas,
              testing::ElementsAre(Seq_delta{1, rA, rC},
                                   Seq_delta{2, rA, rG}));
  EXPECT_THAT(delta.missations, testing::IsEmpty());
}

TEST(Sequence_utils_test, calculate_delta_from_reference_muts_with_compatible_ambs) {
  auto ref_seq = Real_sequence{rA, rC, rG};
  auto seq = Sequence{sG, sN, sN};

  const auto delta = calculate_delta_from_reference(seq, ref_seq);

  EXPECT_THAT(delta.seq_deltas,
              testing::ElementsAre(Seq_delta{0, rA, rG}));
  EXPECT_THAT(estd::ranges::to_vec(delta.missations.slow_elements(ref_seq)),
              testing::ElementsAre(Missation{1, rC}, Missation{2, rG}));
}

TEST(Sequence_utils_test, calculate_delta_from_reference_muts_with_incompatible_ambs) {
  // incompatible ambs get changed to N, which is always compatible with the reference sequence
  auto ref_seq = Real_sequence{rA, rC, rG};
  auto seq = Sequence{sG, sG | sT, sA | sT};

  const auto delta = calculate_delta_from_reference(seq, ref_seq);

  EXPECT_THAT(delta.seq_deltas,
              testing::ElementsAre(Seq_delta{0, rA, rG}));
  EXPECT_THAT(estd::ranges::to_vec(delta.missations.slow_elements(ref_seq)),
              testing::ElementsAre(Missation{1, rC}, Missation{2, rG}));
}

}  // namespace delphy
