#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "io.h"

#include <sstream>

#include "dates.h"

namespace delphy {

TEST(Io_test, read_maple) {
  using enum Real_seq_letter;

  auto is = std::istringstream{
    
    ">reference\n"
    "ACGTACGT\n"
    
    ">Seq1|2024-07-31\n"   // ACCTNNGT
    "C 3\n"
    "N 5 2\n"
    
    ">Seq2|2024-06-01\n"
    "A 8\n"
  };

  auto maple_file = read_maple(is);

  EXPECT_THAT(maple_file.ref_sequence, testing::ElementsAre(A, C, G, T, A, C, G, T));
  EXPECT_THAT(maple_file.tip_descs, testing::SizeIs(2));
  
  EXPECT_THAT(maple_file.tip_descs[0].name, testing::Eq("Seq1|2024-07-31"));
  EXPECT_THAT(maple_file.tip_descs[0].t_min, testing::Eq(parse_iso_date("2024-07-16")));
  EXPECT_THAT(maple_file.tip_descs[0].t, testing::Eq(parse_iso_date("2024-07-31")));
  EXPECT_THAT(maple_file.tip_descs[0].t_max, testing::Eq(parse_iso_date("2024-08-15")));
  EXPECT_THAT(maple_file.tip_descs[0].seq_deltas, testing::ElementsAre(Seq_delta{2, G, C}));
  EXPECT_THAT(maple_file.tip_descs[0].missations.intervals, testing::ElementsAre(Site_interval{4, 6}));
  EXPECT_THAT(maple_file.tip_descs[0].missations.from_states, testing::IsEmpty());
  
  EXPECT_THAT(maple_file.tip_descs[1].name, testing::Eq("Seq2|2024-06-01"));
  EXPECT_THAT(maple_file.tip_descs[1].t_min, testing::Eq(parse_iso_date("2024-05-17")));
  EXPECT_THAT(maple_file.tip_descs[1].t, testing::Eq(parse_iso_date("2024-06-01")));
  EXPECT_THAT(maple_file.tip_descs[1].t_max, testing::Eq(parse_iso_date("2024-06-16")));
  EXPECT_THAT(maple_file.tip_descs[1].seq_deltas, testing::ElementsAre(Seq_delta{7, T, A}));
  EXPECT_THAT(maple_file.tip_descs[1].missations.intervals, testing::IsEmpty());
  EXPECT_THAT(maple_file.tip_descs[1].missations.from_states, testing::IsEmpty());
}

}  // namespace delphy
