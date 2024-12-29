#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "interval_set.h"

namespace delphy {

TEST(Interval_set_test, empty) {
  auto s = Scratch_interval_set{};

  EXPECT_EQ(s.begin(), s.end());
}

TEST(Interval_set_test, nonoverlapping_inserts) {
  auto s = Scratch_interval_set{};

  auto [it, inserted] = s.insert({100, 200});
  EXPECT_THAT(inserted, testing::IsTrue());
  EXPECT_THAT(it->first, testing::Eq(100));
  EXPECT_THAT(it->second, testing::Eq(200));

  auto [it2, inserted2] = s.insert({300, 400});
  EXPECT_THAT(inserted2, testing::IsTrue());
  EXPECT_THAT(it2->first, testing::Eq(300));
  EXPECT_THAT(it2->second, testing::Eq(400));

  EXPECT_THAT(s, testing::ElementsAre(Site_interval{100, 200}, Site_interval{300, 400}));
}

TEST(Interval_set_test, nonoverlapping_site_inserts) {
  auto s = Scratch_interval_set{};

  auto [it, inserted] = s.insert(100);
  EXPECT_THAT(inserted, testing::IsTrue());
  EXPECT_THAT(it->first, testing::Eq(100));
  EXPECT_THAT(it->second, testing::Eq(101));

  auto [it2, inserted2] = s.insert(300);
  EXPECT_THAT(inserted2, testing::IsTrue());
  EXPECT_THAT(it2->first, testing::Eq(300));
  EXPECT_THAT(it2->second, testing::Eq(301));

  EXPECT_THAT(s, testing::ElementsAre(Site_interval{100, 101}, Site_interval{300, 301}));
}

TEST(Interval_set_test, construct_from_initializer_list) {
  auto s = Scratch_interval_set{{Site_interval{100, 200}, Site_interval{300, 400}}};

  EXPECT_THAT(s, testing::ElementsAre(Site_interval{100, 200}, Site_interval{300, 400}));
}

TEST(Interval_set_test, copy_assign_from_different_allocator) {
  auto s = Scratch_interval_set{{Site_interval{100, 200}, Site_interval{300, 400}}};

  auto s2 = Interval_set<>{};  // Not a Scratch_interval_set!
  
  s2 = s;

  EXPECT_THAT(s2, testing::ElementsAre(Site_interval{100, 200}, Site_interval{300, 400}));
}

TEST(Interval_set_test, contains) {
  auto s = Scratch_interval_set{{Site_interval{100, 200}, Site_interval{300, 400}}};

  EXPECT_THAT(s.contains(50), testing::IsFalse());
  
  EXPECT_THAT(s.contains(99), testing::IsFalse());
  EXPECT_THAT(s.contains(100), testing::IsTrue());
  EXPECT_THAT(s.contains(101), testing::IsTrue());
  
  EXPECT_THAT(s.contains(199), testing::IsTrue());
  EXPECT_THAT(s.contains(200), testing::IsFalse());
  EXPECT_THAT(s.contains(201), testing::IsFalse());
  
  EXPECT_THAT(s.contains(299), testing::IsFalse());
  EXPECT_THAT(s.contains(300), testing::IsTrue());
  EXPECT_THAT(s.contains(301), testing::IsTrue());
  
  EXPECT_THAT(s.contains(399), testing::IsTrue());
  EXPECT_THAT(s.contains(400), testing::IsFalse());
  EXPECT_THAT(s.contains(401), testing::IsFalse());
}

TEST(Interval_set_test, overlapping_intervals_closed) {
  auto s = Scratch_interval_set{{Site_interval{100, 200}, Site_interval{300, 400}}};

  // s.begin() == it_A              it_B                it_C == s.end()
  //               |                  |
  //               V                  V
  //       --------[---------]--------[---------]--------
  //              100       200      300       400
  
  auto it_A = s.begin();
  auto it_B = std::next(it_A);
  auto it_C = s.end();

  EXPECT_THAT(s.overlapping_intervals_closed({ 50,  99}), testing::Pair(it_A, it_A));
  
  EXPECT_THAT(s.overlapping_intervals_closed({ 50, 100}), testing::Pair(it_A, it_B));
  
  EXPECT_THAT(s.overlapping_intervals_closed({ 50, 150}), testing::Pair(it_A, it_B));
  EXPECT_THAT(s.overlapping_intervals_closed({ 50, 450}), testing::Pair(it_A, it_C));
  EXPECT_THAT(s.overlapping_intervals_closed({150, 250}), testing::Pair(it_A, it_B));
  
  EXPECT_THAT(s.overlapping_intervals_closed({250, 270}), testing::Pair(it_B, it_B));
  
  EXPECT_THAT(s.overlapping_intervals_closed({250, 350}), testing::Pair(it_B, it_C));
  EXPECT_THAT(s.overlapping_intervals_closed({450, 500}), testing::Pair(it_C, it_C));
}

TEST(Interval_set_test, strictly_overlapping_inserts) {
  auto s = Scratch_interval_set{};

  auto [it, inserted] = s.insert({100, 200});
  EXPECT_THAT(inserted, testing::IsTrue());
  
  auto [it2, inserted2] = s.insert({150, 250});
  EXPECT_THAT(inserted2, testing::IsTrue());
  EXPECT_THAT(it2->first, testing::Eq(100));
  EXPECT_THAT(it2->second, testing::Eq(250));

  EXPECT_THAT(s, testing::ElementsAre(Site_interval{100, 250}));
}

TEST(Interval_set_test, fully_contained_overlapping_inserts) {
  auto s = Scratch_interval_set{};

  auto [it, inserted] = s.insert({100, 200});
  EXPECT_THAT(inserted, testing::IsTrue());
  
  auto [it2, inserted2] = s.insert({125, 175});
  EXPECT_THAT(inserted2, testing::IsFalse());
  EXPECT_THAT(it2->first, testing::Eq(100));
  EXPECT_THAT(it2->second, testing::Eq(200));

  EXPECT_THAT(s, testing::ElementsAre(Site_interval{100, 200}));
}

TEST(Interval_set_test, complex_overlapping_inserts) {
  auto s = Scratch_interval_set{{Site_interval{100, 200}, Site_interval{300, 400}}};
  
  auto [it2, inserted2] = s.insert({150, 450});
  EXPECT_THAT(inserted2, testing::IsTrue());
  EXPECT_THAT(it2->first, testing::Eq(100));
  EXPECT_THAT(it2->second, testing::Eq(450));

  EXPECT_THAT(s, testing::ElementsAre(Site_interval{100, 450}));
}

TEST(Interval_set_test, tough_merging_inserts) {
  auto s = Scratch_interval_set{{Site_interval{100, 200}, Site_interval{300, 400}}};
  
  auto [it2, inserted2] = s.insert({200, 300});
  EXPECT_THAT(inserted2, testing::IsTrue());
  EXPECT_THAT(it2->first, testing::Eq(100));
  EXPECT_THAT(it2->second, testing::Eq(400));

  EXPECT_THAT(s, testing::ElementsAre(Site_interval{100, 400}));
}

TEST(Interval_set_test, consecutive_site_inserts) {
  // This is very inefficient, but it *should* work
  auto s = Scratch_interval_set{};

  s.insert(10);
  s.insert(11);
  s.insert(12);
  
  s.insert(20);
  s.insert(21);
  s.insert(22);

  EXPECT_THAT(s, testing::ElementsAre(Site_interval{10, 13}, Site_interval{20, 23}));
}

TEST(Interval_set_test, merge) {
  //
  //                                                                                                         1111111111
  //               1111111111222222222233333333334444444444555555555566666666667777777777888888888899999999990000000000
  //     01234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
  // s1: ----------[---------)-[----)--[---------)---------[---------)---------[---------)---------[---------)---------
  // s2: --[----)--------------------[-)--------------[---------)-------------------[---------)-[--)---------[-)-------
  //  m: --[----)--[---------)-[----)[-----------)----[--------------)---------[--------------)-[--------------)-------
  //
  auto s1 = Scratch_interval_set{
    {Site_interval{10,  20},
     Site_interval{22,  27},
     Site_interval{30,  40},
     Site_interval{50,  60},
     Site_interval{70,  80},
     Site_interval{90, 100}}};
  auto s2 = Scratch_interval_set{
    {Site_interval{  2,   7},
     Site_interval{ 28,  30},
     Site_interval{ 45,  55},
     Site_interval{ 75,  85},
     Site_interval{ 87,  90},
     Site_interval{100, 102}}};

  EXPECT_THAT(merge_interval_sets(s1, s2), testing::ElementsAre(
      Site_interval{ 2,   7},
      Site_interval{10,  20},
      Site_interval{22,  27},
      Site_interval{28,  40},
      Site_interval{45,  60},
      Site_interval{70,  85},
      Site_interval{87, 102}));
}

TEST(Interval_set_test, intersect) {
  //
  //                                                                                                         1111111111
  //               1111111111222222222233333333334444444444555555555566666666667777777777888888888899999999990000000000
  //     01234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
  // s1: ----------[---------)-[----)--[---------)---------[---------)---------[---------)---------[---------)---------
  // s2: --[----)--------------------[-)--------------[---------)-------------------[---------)-[--)---------[-)-------
  //  o: --------------------------------------------------[----)-------------------[----)-----------------------------
  //
  auto s1 = Scratch_interval_set{
    {Site_interval{10,  20},
     Site_interval{22,  27},
     Site_interval{30,  40},
     Site_interval{50,  60},
     Site_interval{70,  80},
     Site_interval{90, 100}}};
  auto s2 = Scratch_interval_set{
    {Site_interval{  2,   7},
     Site_interval{ 28,  30},
     Site_interval{ 45,  55},
     Site_interval{ 75,  85},
     Site_interval{ 87,  90},
     Site_interval{100, 102}}};

  EXPECT_THAT(intersect_interval_sets(s1, s2), testing::ElementsAre(
      Site_interval{50,  55},
      Site_interval{75,  80}));

  // interval_sets_intersect
  EXPECT_THAT(interval_sets_intersect(Interval_set{{Site_interval{10, 20}}},
                                      Interval_set{{Site_interval{30, 40}}}), testing::IsFalse());
  EXPECT_THAT(interval_sets_intersect(s1, s2), testing::IsTrue());
}

TEST(Interval_set_test, subtract) {
  //
  //                                                                                                         1111111111
  //               1111111111222222222233333333334444444444555555555566666666667777777777888888888899999999990000000000
  //     01234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
  // s1: ----------[---------)-[----)--[---------)---------[---------)---------[---------)---------[---------)---------
  // s2: --[----)--------------------[-)--------------[---------)-------------------[---------)-[--)---------[-)-------
  //1-2: ----------[---------)-[----)--[---------)--------------[----)---------[----)---- ---------[---------)---------
  //
  //                                                                                                         1111111111
  //               1111111111222222222233333333334444444444555555555566666666667777777777888888888899999999990000000000
  //     01234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
  // s2: --[----)--------------------[-)--------------[---------)-------------------[---------)-[--)---------[-)-------
  // s1: ----------[---------)-[----)--[---------)---------[---------)---------[---------)---------[---------)---------
  //2-1: --[----)--------------------[-)--------------[----)-----------------------------[----)-[--)---------[-)-------
  //
  auto s1 = Scratch_interval_set{
    {Site_interval{10,  20},
     Site_interval{22,  27},
     Site_interval{30,  40},
     Site_interval{50,  60},
     Site_interval{70,  80},
     Site_interval{90, 100}}};
  auto s2 = Scratch_interval_set{
    {Site_interval{  2,   7},
     Site_interval{ 28,  30},
     Site_interval{ 45,  55},
     Site_interval{ 75,  85},
     Site_interval{ 87,  90},
     Site_interval{100, 102}}};

  EXPECT_THAT(subtract_interval_sets(s1, s2), testing::ElementsAre(
      Site_interval{10,  20},
      Site_interval{22,  27},
      Site_interval{30,  40},
      Site_interval{55,  60},
      Site_interval{70,  75},
      Site_interval{90, 100}));

  EXPECT_THAT(subtract_interval_sets(s2, s1), testing::ElementsAre(
      Site_interval{  2,   7},
      Site_interval{ 28,  30},
      Site_interval{ 45,  50},
      Site_interval{ 80,  85},
      Site_interval{ 87,  90},
      Site_interval{100, 102}));
}

TEST(Interval_set_test, subtract_edges) {
  //
  //                                                                                                         1111111111
  //               1111111111222222222233333333334444444444555555555566666666667777777777888888888899999999990000000000
  //     01234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
  // s1: ---[--)------[--)------[--)------[--)------[--)------[--)------[--)------[--)------[--)------[--)------[--)---
  // s2: [--)------[-----)---[-------)----[-)-------[--)------[----)-----[)--------[-)-------[--)--------[-)-----------
  //1-2: ---[--)----------------------------[)--------------------------[)[)------[)--------[)--------[--)------[--)---
  //
  auto s1 = Scratch_interval_set{
    {Site_interval{  3,   6},
     Site_interval{ 13,  16},
     Site_interval{ 23,  26},
     Site_interval{ 33,  36},
     Site_interval{ 43,  46},
     Site_interval{ 53,  56},
     Site_interval{ 63,  66},
     Site_interval{ 73,  76},
     Site_interval{ 83,  86},
     Site_interval{ 93,  96},
     Site_interval{103, 106}}};
  auto s2 = Scratch_interval_set{
    {Site_interval{  0,   3},
     Site_interval{ 10,  16},
     Site_interval{ 20,  28},
     Site_interval{ 33,  35},
     Site_interval{ 43,  46},
     Site_interval{ 53,  58},
     Site_interval{ 64,  65},
     Site_interval{ 74,  76},
     Site_interval{ 84,  87},
     Site_interval{ 96,  98}}};

  EXPECT_THAT(subtract_interval_sets(s1, s2), testing::ElementsAre(
      Site_interval{  3,   6},
      Site_interval{ 35,  36},
      Site_interval{ 63,  64},
      Site_interval{ 65,  66},
      Site_interval{ 73,  74},
      Site_interval{ 83,  84},
      Site_interval{ 93,  96},
      Site_interval{103, 106}));
}

TEST(Interval_set_test, slow_elements) {
  auto s = Scratch_interval_set{{Site_interval{2, 5}, Site_interval{7, 8}}};

  EXPECT_THAT(estd::ranges::to_vec(s.slow_elements()), testing::ElementsAre(2, 3, 4, 7));
}


}  // namespace delphy
