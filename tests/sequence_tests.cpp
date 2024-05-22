#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "sequence.h"

namespace delphy {

TEST(Sequence_test, real_seq_letter_from_index_valid) {
  using enum Real_seq_letter;

  EXPECT_EQ(real_seq_letter_from_index(0), A);
  EXPECT_EQ(real_seq_letter_from_index(1), C);
  EXPECT_EQ(real_seq_letter_from_index(2), G);
  EXPECT_EQ(real_seq_letter_from_index(3), T);
}

TEST(Sequence_test, real_seq_letter_from_index_invalid) {
  using enum Real_seq_letter;

  EXPECT_THROW(real_seq_letter_from_index(-1), std::invalid_argument);
  EXPECT_THROW(real_seq_letter_from_index(4), std::invalid_argument);
}

TEST(Sequence_test, to_char_seq_letter_valid) {
  using namespace Seq_letters;
  
  EXPECT_EQ(to_char(A), 'A');
  EXPECT_EQ(to_char(C), 'C');
  EXPECT_EQ(to_char(G), 'G');
  EXPECT_EQ(to_char(T), 'T');
  
  EXPECT_EQ(to_char(A | _ | G | _), 'R');
  EXPECT_EQ(to_char(_ | C | _ | T), 'Y');
  EXPECT_EQ(to_char(_ | C | G | _), 'S');
  EXPECT_EQ(to_char(A | _ | _ | T), 'W');
  EXPECT_EQ(to_char(_ | _ | G | T), 'K');
  EXPECT_EQ(to_char(A | C | _ | _), 'M');
  
  EXPECT_EQ(to_char(_ | C | G | T), 'B');
  EXPECT_EQ(to_char(A | _ | G | T), 'D');
  EXPECT_EQ(to_char(A | C | _ | T), 'H');
  EXPECT_EQ(to_char(A | C | G | _), 'V');
  
  EXPECT_EQ(to_char(A | C | G | T), 'N');
}

TEST(Sequence_test, to_char_seq_letter_invalid) {
  EXPECT_EQ(to_char(Seq_letter{0b1111'1111}), '?');
}

TEST(Sequence_test, to_char_real_seq_letter_valid) {
  using enum Real_seq_letter;
  
  EXPECT_EQ(to_char(A), 'A');
  EXPECT_EQ(to_char(C), 'C');
  EXPECT_EQ(to_char(G), 'G');
  EXPECT_EQ(to_char(T), 'T');
}

TEST(Sequence_test, to_char_real_seq_letter_invalid) {
  EXPECT_EQ(to_char(Real_seq_letter{0b1111'1111}), '?');
}

TEST(Sequence_test, to_real_seq_letter_valid) {
  EXPECT_EQ(to_real_seq_letter(Seq_letters::A), Real_seq_letter::A);
  EXPECT_EQ(to_real_seq_letter(Seq_letters::C), Real_seq_letter::C);
  EXPECT_EQ(to_real_seq_letter(Seq_letters::G), Real_seq_letter::G);
  EXPECT_EQ(to_real_seq_letter(Seq_letters::T), Real_seq_letter::T);
}

TEST(Sequence_test, to_real_seq_letter_ambiguous) {
  EXPECT_THROW(to_real_seq_letter(Seq_letters::A | Seq_letters::C), std::invalid_argument);
  EXPECT_THROW(to_real_seq_letter(Seq_letters::_), std::invalid_argument);
}

TEST(Sequence_test, seq_vector_int_printing) {
  auto v = Seq_vector{1,2,3,4};
  auto ss = std::stringstream{};
  ss << v;
  EXPECT_THAT(ss.str(), testing::StrEq("Seq_vector{1,2,3,4}"));
}

TEST(Sequence_test, seq_vector_double_printing) {
  auto v = Seq_vector{0.0, 0.25, 0.3, 1e25};
  auto ss = std::stringstream{};
  ss << v;
  EXPECT_THAT(ss.str(), testing::StrEq("Seq_vector{0,0.25,0.3,1e+25}"));
}

TEST(Sequence_test, seq_matrix_int_printing) {
  using enum Real_seq_letter;
  
  auto M = Seq_matrix<int>{};
  M[A][A] = 11;  M[A][C] = 12;  M[A][G] = 13;  M[A][T] = 14;
  M[C][A] = 21;  M[C][C] = 22;  M[C][G] = 23;  M[C][T] = 24;
  M[G][A] = 31;  M[G][C] = 32;  M[G][G] = 33;  M[G][T] = 34;
  M[T][A] = 41;  M[T][C] = 42;  M[T][G] = 43;  M[T][T] = 44;
  
  auto ss = std::stringstream{};
  ss << M;
  EXPECT_THAT(ss.str(), testing::StrEq("Seq_matrix{{11,12,13,14}, {21,22,23,24}, {31,32,33,34}, {41,42,43,44}}"));
}

TEST(Sequence_test, seq_matrix_double_printing) {
  using enum Real_seq_letter;
  
  auto M = Seq_matrix<double>{};
  M[A][A] = 11.5;  M[A][C] = 12.5;  M[A][G] = 13.5;  M[A][T] = 14.5;
  M[C][A] = 21.6;  M[C][C] = 22.6;  M[C][G] = 23.6;  M[C][T] = 24.6;
  M[G][A] = 31.7;  M[G][C] = 32.7;  M[G][G] = 33.7;  M[G][T] = 34.7;
  M[T][A] = 41.8;  M[T][C] = 42.8;  M[T][G] = 43.8;  M[T][T] = 44.8;
  
  auto ss = std::stringstream{};
  ss << M;
  EXPECT_THAT(ss.str(), testing::StrEq(
      "Seq_matrix{{11.5,12.5,13.5,14.5}, {21.6,22.6,23.6,24.6}, {31.7,32.7,33.7,34.7}, {41.8,42.8,43.8,44.8}}"));
}

}  // namespace delphy
