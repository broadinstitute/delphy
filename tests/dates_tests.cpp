#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <absl/time/time.h>

#include "dates.h"

namespace delphy {

#define test_round_trip(str, t) do {\
    EXPECT_THAT(parse_iso_date(str), testing::DoubleNear(t, 1e-6));     \
    EXPECT_THAT(to_iso_date(t), testing::StrEq(str));                   \
  } while(false)

TEST(Dates_test, round_trip) {
  test_round_trip("2020-01-01", 0);
  test_round_trip("2021-01-01", 366);
  test_round_trip("2021-03-19", 366 - 1 + 31 + 28 + 19);
  test_round_trip("2021-09-28", 366 - 1 + 31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 28);
  test_round_trip("2021-11-21", 366 - 1 + 31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31 + 21);

  // Dates outside (2016, 2024) now that we use Abseil Time instead of hand-coding arithmetic
  test_round_trip("2019-01-01", -365);
  test_round_trip("2018-01-01", -365-365);
  test_round_trip("2017-01-01", -365-365-365);
  test_round_trip("2016-01-01", -365-365-365-366);
  test_round_trip("2024-01-01", +366+365+365+365);
  test_round_trip("2025-01-01", +366+365+365+365+366);

  auto month_lengths = std::array{0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};  // 2021
  auto mm_start = 366;  // 2021-01-01
  for (auto mm = 1; mm <= 12; ++mm) {
    for (auto dd = 1; dd <= month_lengths[mm]; ++dd) {
      auto expected = absl::StrFormat("2021-%02d-%02d", mm, dd);
      SCOPED_TRACE(expected);
      test_round_trip(expected, mm_start + dd - 1);
    }
    mm_start += month_lengths[mm];
  }

  test_round_trip("2021-03-19", 366 - 1 + 31 + 28 + 19);
  test_round_trip("2021-09-27", 366 - 1 + 31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 27);
  test_round_trip("2021-09-28", 366 - 1 + 31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 28);
  test_round_trip("2021-11-21", 366 - 1 + 31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31 + 21);
}

TEST(Dates_test, absl_dates) {
  auto my_epoch = absl::CivilDay(2020, 1, 1);
  EXPECT_EQ(absl::CivilDay(2020,  2,  1) - my_epoch, 31);
  EXPECT_EQ(absl::CivilDay(2020, 12, 15) - my_epoch, 365-31+15);
  EXPECT_EQ(absl::CivilDay(2019,  1,  1) - my_epoch, -365-1+1);
  EXPECT_EQ(absl::CivilDay(2019,  2, 14) - my_epoch, -365-1+31+14);

  auto utc = absl::UTCTimeZone();
  auto fmt = "%a, %e %b %Y %H:%M:%S %z (%Z)";
  auto t = absl::FromCivil(absl::CivilHour(2016,  1,  1, 12), utc);
  EXPECT_EQ("Fri,  1 Jan 2016 12:00:00 +0000 (UTC)",
            absl::FormatTime(fmt, t, utc));

  auto t1 = absl::FromCivil(absl::CivilHour(2016,  1,  1, 12), utc);
  auto t2 = absl::FromCivil(absl::CivilHour(2020,  1,  1, 12), utc);
  EXPECT_EQ(absl::ToDoubleHours(t2 - t1) / 24.0, 1461.0);
}

TEST(Dates_test, parse_iso_month) {
  // 31 days
  EXPECT_THAT(parse_iso_month("2024-12"), testing::Eq(std::pair{parse_iso_date("2024-12-01"), parse_iso_date("2025-01-01")}));
  
  // 30 days
  EXPECT_THAT(parse_iso_month("2024-11"), testing::Eq(std::pair{parse_iso_date("2024-11-01"), parse_iso_date("2024-12-01")}));
  
  // 28 days
  EXPECT_THAT(parse_iso_month("2023-02"), testing::Eq(std::pair{parse_iso_date("2023-02-01"), parse_iso_date("2023-03-01")}));
  
  // 29 days on leap year
  EXPECT_THAT(parse_iso_month("2024-02"), testing::Eq(std::pair{parse_iso_date("2024-02-01"), parse_iso_date("2024-03-01")}));
  
  // The year 1900 was not a leap year
  EXPECT_THAT(parse_iso_month("1900-02"), testing::Eq(std::pair{parse_iso_date("1900-02-01"), parse_iso_date("1900-03-01")}));
  
  // But the year 2000 was (multiple of 400)
  EXPECT_THAT(parse_iso_month("2000-02"), testing::Eq(std::pair{parse_iso_date("2000-02-01"), parse_iso_date("2000-03-01")}));
}

TEST(Dates_test, parse_iso_year) {
  // The year 2024 was a leap year
  EXPECT_THAT(parse_iso_year("2024"), testing::Eq(std::pair{parse_iso_date("2024-01-01"), parse_iso_date("2025-01-01")}));
  
  // The year 2023 was not a leap year
  EXPECT_THAT(parse_iso_year("2023"), testing::Eq(std::pair{parse_iso_date("2023-01-01"), parse_iso_date("2024-01-01")}));
  
  // The year 1900 was not a leap year
  EXPECT_THAT(parse_iso_year("1900"), testing::Eq(std::pair{parse_iso_date("1900-01-01"), parse_iso_date("1901-01-01")}));
  
  // But the year 2000 was (multiple of 400)
  EXPECT_THAT(parse_iso_year("2000"), testing::Eq(std::pair{parse_iso_date("2000-01-01"), parse_iso_date("2001-01-01")}));
}

TEST(Dates_test, to_linear_year) {
  // The year 2024 was a leap year
  EXPECT_THAT(to_linear_year(parse_iso_date("2024-02-01")), testing::DoubleNear(2024.0 + 31.0 / 366.0, 1e-6));
  EXPECT_THAT(to_linear_year(parse_iso_date("2024-03-01")), testing::DoubleNear(2024.0 + (31.0 + 29.0) / 366.0, 1e-6));
  EXPECT_THAT(to_linear_year(parse_iso_date("2024-03-15")), testing::DoubleNear(2024.0 + (31.0 + 29.0 + 15-1) / 366.0, 1e-6));

  // The year 2023 was not a leap year
  EXPECT_THAT(to_linear_year(parse_iso_date("2023-02-01")), testing::DoubleNear(2023.0 + 31.0 / 365.0, 1e-6));
  EXPECT_THAT(to_linear_year(parse_iso_date("2023-03-01")), testing::DoubleNear(2023.0 + (31.0 + 28.0) / 365.0, 1e-6));

  // The year 1900 was not a leap year
  EXPECT_THAT(to_linear_year(parse_iso_date("1900-02-01")), testing::DoubleNear(1900.0 + 31.0 / 365.0, 1e-6));
  EXPECT_THAT(to_linear_year(parse_iso_date("1900-03-01")), testing::DoubleNear(1900.0 + (31.0 + 28.0) / 365.0, 1e-6));

  // But the year 2000 was (multiple of 400)
  EXPECT_THAT(to_linear_year(parse_iso_date("2000-02-01")), testing::DoubleNear(2000.0 + 31.0 / 366.0, 1e-6));
  EXPECT_THAT(to_linear_year(parse_iso_date("2000-03-01")), testing::DoubleNear(2000.0 + (31.0 + 29.0) / 366.0, 1e-6));
}

}  // namespace delphy
