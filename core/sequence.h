#ifndef DELPHY_SEQUENCE_H_
#define DELPHY_SEQUENCE_H_

#include <cassert>
#include <array>
#include <algorithm>
#include <cmath>
#include <vector>
#include <random>

#include "absl/log/check.h"
#include "absl/random/bit_gen_ref.h"
#include "absl/random/distributions.h"
#include "absl/strings/str_format.h"

#include "estd.h"

namespace delphy {

using Seq_letter = uint8_t;
namespace Seq_letters {
inline constexpr Seq_letter
    _{0b0000}, none{_},
    A{0b0001},
    C{0b0010},
    G{0b0100},
    T{0b1000}, U{T},
    N{0b1111};
}

using Sequence = std::vector<Seq_letter>;

inline auto is_ambiguous(Seq_letter s) -> bool {
  return (s & (s - 1)) != 0;  // "not a power of 2"
}

inline auto to_char(Seq_letter s) noexcept -> char {
    using namespace Seq_letters;

    switch (s) {

        // One base
    case A: return 'A';
    case C: return 'C';
    case G: return 'G';
    case T: return 'T';

        // One of two bases
    case A | G: return 'R';
    case C | T: return 'Y';
    case G | C: return 'S';
    case A | T: return 'W';
    case G | T: return 'K';
    case A | C: return 'M';

        // One of three bases
    case _ | C | G | T: return 'B';
    case A | _ | G | T: return 'D';
    case A | C | _ | T: return 'H';
    case A | C | G | _: return 'V';

        // Any base
    case A | C | G | T: return 'N';

    default: return '?';
    }
}
inline auto operator<<(std::ostream& os, Seq_letter s) -> std::ostream& {
  auto c = to_char(s);
  if (c != '?') {
    return os << c;
  } else {
    return os << absl::StreamFormat("[Seq_letter: 0x%02x]", s);
  }
}

inline auto to_seq_letter(char c) -> Seq_letter {
  using namespace Seq_letters;
  
  // See https://en.wikipedia.org/wiki/FASTA_format#Sequence_representation
  // and https://www.bioinformatics.org/sms/iupac.html
  switch (std::toupper(c)) {
    
    // One base
    case 'A': return A;
    case 'C': return C;
    case 'G': return G;
    case 'T': return T;
    case 'U': return T;
      
      // One of two bases
    case 'R': return A | G;
    case 'Y': return C | T;
    case 'S': return G | C;
    case 'W': return A | T;
    case 'K': return G | T;
    case 'M': return A | C;
      
      // One of three bases
    case 'B': return _ | C | G | T;
    case 'D': return A | _ | G | T;
    case 'H': return A | C | _ | T;
    case 'V': return A | C | G | _;
      
    case 'N': // Any base
    case '-': // Any base
    case '?': // Any base
    case '.': // Gap
      return A | C | G | T;
      
    default: return _;
  }
}

enum class Real_seq_letter : uint8_t { A = 0, C = 1, G = 2, T = 3 };
inline constexpr int k_num_real_seq_letters = 4;
inline auto index_of(Real_seq_letter c) -> int { return static_cast<int>(c); }
inline auto real_seq_letter_from_index(int index) -> Real_seq_letter {
  if (index < 0 || index > 3) {
    throw std::invalid_argument(absl::StrFormat(
        "Cannot convert index %d to Real_seq_letter", index));
  }
  return Real_seq_letter{static_cast<std::underlying_type_t<Real_seq_letter>>(index)};
}
inline constexpr std::initializer_list<Real_seq_letter> k_all_real_seq_letters{
    Real_seq_letter::A, Real_seq_letter::C, Real_seq_letter::G, Real_seq_letter::T};

inline auto to_char(Real_seq_letter s) -> char {
  using enum Real_seq_letter;
  switch (s) {
    case A: return 'A';
    case C: return 'C';
    case G: return 'G';
    case T: return 'T';
    default: return '?';
  }
}
inline auto operator<<(std::ostream& os, Real_seq_letter s) -> std::ostream& {
  auto c = to_char(s);
  if (c != '?') {
    return os << c;
  } else {
    return os << absl::StreamFormat("[Real_seq_letter: 0x%02x]", s);
  }
}

inline auto char_to_real_seq_letter(char c) -> Real_seq_letter {
  using enum Real_seq_letter;
  switch (std::toupper(c)) {
    case 'A': return A;
    case 'C': return C;
    case 'G': return G;
    case 'T': return T;
    default:
      throw std::runtime_error(
          absl::StrFormat("Cannot convert character '%c' to Real_seq_letter", c));
  }
}
inline auto to_seq_letter(Real_seq_letter s) -> Seq_letter {
  switch (s) {
    case Real_seq_letter::A: return Seq_letters::A;
    case Real_seq_letter::C: return Seq_letters::C;
    case Real_seq_letter::G: return Seq_letters::G;
    case Real_seq_letter::T: return Seq_letters::T;
    default:
      throw std::runtime_error(
          absl::StrFormat("Unrecognized Real_seq_letter %d", index_of(s)));
  }
}

inline auto to_real_seq_letter(Seq_letter s) -> Real_seq_letter {
  switch (s) {
    case Seq_letters::A: return Real_seq_letter::A;
    case Seq_letters::C: return Real_seq_letter::C;
    case Seq_letters::G: return Real_seq_letter::G;
    case Seq_letters::T: return Real_seq_letter::T;
    default:
      throw std::invalid_argument(absl::StrFormat("Ambiguous Seq_letter %c in to_real_seq_letter", to_char(s)));
  }
}

using Real_sequence = std::vector<Real_seq_letter>;

using Site_index = int;
template<typename T>
using Site_vector = std::vector<T>;

// Seq_matrix<T> is a matrix with elements of type T where row and column indices are A, C, G or T.
// Examples include mutation rate matrices, mutation count matrices and evolution matrices.
// See Peter Gottschling, Discovering Modern C++ (2nd ed), Section 6.6.1.2 for the idea behind the
// CRTP proxy object Seq_matrix_row.  It allows for Seq_matrix elements to be accessed using double
// brackets, i.e., m[a][b], while enforcing that indices are of type Real_seq_letter.
template<typename T>
class Seq_vector;

template<typename Matrix>
class Seq_matrix_row {
  using T = typename Matrix::value_type;

 public:
  Seq_matrix_row(Matrix& matrix, Real_seq_letter a) : matrix_{matrix}, a_{a} {}

  auto operator[](Real_seq_letter b) -> auto& { return matrix_.data_[index_of(a_)][index_of(b)]; }

  operator Seq_vector<T>() const;

 private:
  Matrix& matrix_;
  Real_seq_letter a_;
};

template<typename T>
class Seq_matrix {
 public:
  using value_type = T;

  Seq_matrix() : data_{} {}
  explicit Seq_matrix(const T& initial_values) : data_{} {
    for (auto& row : data_) {
      row.fill(initial_values);
    }
  }
  
  auto operator[](Real_seq_letter c) -> Seq_matrix_row<Seq_matrix> { return {*this, c}; }
  auto operator[](Real_seq_letter c) const -> Seq_matrix_row<const Seq_matrix> { return {*this, c}; }

  template<typename U>
  auto operator+=(const Seq_matrix<U>& rhs) -> Seq_matrix<T>& {
    for (auto a : k_all_real_seq_letters) {
      for (auto b : k_all_real_seq_letters) {
        (*this)[a][b] += rhs[a][b];
      }
    }
    return *this;
  }

 private:
  friend class Seq_matrix_row<Seq_matrix>;
  friend class Seq_matrix_row<const Seq_matrix>;
  std::array<std::array<T, k_num_real_seq_letters>, k_num_real_seq_letters> data_ = {};
};

template<typename T>
auto operator<<(std::ostream& os, const Seq_matrix<T>& M) -> std::ostream& {
  using enum Real_seq_letter;
  return os << absl::StreamFormat(
      "Seq_matrix{{%v,%v,%v,%v}, {%v,%v,%v,%v}, {%v,%v,%v,%v}, {%v,%v,%v,%v}}",
      M[A][A], M[A][C], M[A][G], M[A][T],
      M[C][A], M[C][C], M[C][G], M[C][T],
      M[G][A], M[G][C], M[G][G], M[G][T],
      M[T][A], M[T][C], M[T][G], M[T][T]);
}

template<typename T>
auto operator==(const Seq_matrix<T>& M, const Seq_matrix<T>& N) -> bool {
  for (auto a : k_all_real_seq_letters) {
    for (auto b : k_all_real_seq_letters) {
      if (M[a][b] != N[a][b]) {
        return false;
      }
    }
  }
  return true;
}

// Seq_vector<T> is a vector with elements of type T indexed by A, C, G or T.
// Examples include base frequency vectors, base count vectors and escape rate vectors.
template<typename T>
class Seq_vector : private std::array<T, k_num_real_seq_letters> {
  using Base = std::array<T, k_num_real_seq_letters>;
 public:
  Seq_vector() = default;
  Seq_vector(T tA, T tC, T tG, T tT) : std::array<T, k_num_real_seq_letters>{tA, tC, tG, tT} {}
  explicit Seq_vector(T tAll) : Seq_vector(tAll, tAll, tAll, tAll) {}

  auto operator[](Real_seq_letter c) -> T& { return Base::operator[](index_of(c)); }
  auto operator[](Real_seq_letter c) const -> const T& { return Base::operator[](index_of(c)); }

  auto operator+=(const Seq_vector<T>& rhs) -> Seq_vector<T>& {
    for (auto a : k_all_real_seq_letters) {
      (*this)[a] += rhs[a];
    }
    return *this;
  }
  auto operator-=(const Seq_vector<T>& rhs) -> Seq_vector<T>& {
    for (auto a : k_all_real_seq_letters) {
      (*this)[a] -= rhs[a];
    }
    return *this;
  }
  auto operator/=(T denom) -> Seq_vector<T>& {
    for (T& coeff : *this) {
      coeff /= denom;
    }
    return *this;
  }

  // The rest of the interface is just like std::array: add using declarations here as needed
  using Base::begin;
  using Base::end;
  using Base::size;
};

template<typename T>
auto operator<<(std::ostream& os, const Seq_vector<T>& v) -> std::ostream& {
  using enum Real_seq_letter;
  return os << absl::StreamFormat("Seq_vector{%v,%v,%v,%v}", v[A], v[C], v[G], v[T]);
}

template<typename Matrix>
Seq_matrix_row<Matrix>::operator Seq_vector<T>() const {
  auto result = Seq_vector<T>{};
  for (auto b : k_all_real_seq_letters) {
    result[b] = matrix_.data_[index_of(a_)][index_of(b)];
  }
  return result;
}

template<typename T>
auto operator+(const Seq_vector<T>& v, const Seq_vector<T>& w) -> Seq_vector<T> {
  auto result = Seq_vector<T>{};
  for (auto a : k_all_real_seq_letters) {
    result[a] = v[a] + w[a];
  }
  return result;
}

template<typename T>
auto operator-(const Seq_vector<T>& v, const Seq_vector<T>& w) -> Seq_vector<T> {
  auto result = Seq_vector<T>{};
  for (auto a : k_all_real_seq_letters) {
    result[a] = v[a] - w[a];
  }
  return result;
}

template<typename T>
auto operator*(const Seq_vector<T>& v, const Seq_vector<T>& w) -> Seq_vector<T> {
  auto result = Seq_vector<T>{};
  for (auto a : k_all_real_seq_letters) {
    result[a] = v[a] * w[a];
  }
  return result;
}

template<typename T>
auto operator*(const Seq_vector<T>& v, T k) -> Seq_vector<T> {
  auto result = Seq_vector<T>{};
  for (auto a : k_all_real_seq_letters) {
    result[a] = v[a] * k;
  }
  return result;
}

template<typename T>
auto operator*(T k, const Seq_vector<T>& v) -> Seq_vector<T> {
  auto result = Seq_vector<T>{};
  for (auto a : k_all_real_seq_letters) {
    result[a] = k * v[a];
  }
  return result;
}

template<typename T>
auto operator*(const Seq_matrix<T>& M, const Seq_vector<T>& v) -> Seq_vector<T> {
  auto result = Seq_vector<T>{};
  for (auto a : k_all_real_seq_letters) {
    auto sum = T{};
    for (auto b : k_all_real_seq_letters) {
      sum += M[a][b] * v[b];
    }
    result[a] = sum;
  }
  return result;
}

template<typename T>
auto operator-(const Seq_vector<T>& v) -> Seq_vector<T> {
  auto result = Seq_vector<T>{};
  for (auto a : k_all_real_seq_letters) {
    result[a] = -v[a];
  }
  return result;
}

template<typename T>
auto operator&&(const Seq_vector<T>& v, const Seq_vector<T>& w) -> Seq_vector<T> {
  auto result = Seq_vector<T>{};
  for (auto a : k_all_real_seq_letters) {
    result[a] = v[a] && w[a];
  }
  return result;
}

template<typename T>
auto operator||(const Seq_vector<T>& v, const Seq_vector<T>& w) -> Seq_vector<T> {
  auto result = Seq_vector<T>{};
  for (auto a : k_all_real_seq_letters) {
    result[a] = v[a] || w[a];
  }
  return result;
}

template<typename T>
auto operator==(const Seq_vector<T>& v, const Seq_vector<T>& w) -> bool {
  for (auto a : k_all_real_seq_letters) {
    if (v[a] != w[a]) {
      return false;
    }
  }
  return true;
}

template<typename T>
auto taxicab_distance(const Seq_vector<T>& v, const Seq_vector<T>& w) -> T {
  using std::abs;
  auto result = T{};
  for (auto a : k_all_real_seq_letters) {
    result += abs(v[a] - w[a]);
  }
  return result;
}

template<typename T>
auto taxicab_distance(const std::vector<T>& v, const std::vector<T>& w) -> double {
  using std::abs;
  CHECK_EQ(std::ssize(v), std::ssize(w));
  auto result = 0.0;
  for (auto i = 0; i != std::ssize(v); ++i) {
    result += taxicab_distance(v[i], w[i]);
  }
  return result;
}

inline auto pick_state(const Seq_vector<double>& w, absl::BitGenRef bitgen) -> Real_seq_letter {
  auto p_ = 0.0;
  auto p_A = p_ + w[Real_seq_letter::A];
  auto p_AC = p_A + w[Real_seq_letter::C];
  auto p_ACG = p_AC + w[Real_seq_letter::G];
  auto p_ACGT = p_ACG + w[Real_seq_letter::T];

  auto r = absl::Uniform(absl::IntervalClosedOpen, bitgen, 0.0, p_ACGT);
  if (r < p_AC) {
    return (r < p_A) ? Real_seq_letter::A : Real_seq_letter::C;
  } else {
    return (r < p_ACG) ? Real_seq_letter::G : Real_seq_letter::T;
  }
}

}  // namespace delphy

#endif // DELPHY_SEQUENCE_H_
