#include "evo_hky.h"

#include <Eigen/Dense>

namespace delphy {

auto Hky_model::derive_site_evo_model() const -> Site_evo_model {
  using enum Real_seq_letter;

  auto r = Eigen::Matrix4d{
      { 0, 1, kappa, 1 },
      { 1, 0, 1, kappa },
      { kappa, 1, 0, 1 },
      { 1, kappa, 1, 0 }
  };
  auto real_pi = Eigen::Vector4d{ pi_a[A], pi_a[C], pi_a[G], pi_a[T] };

  auto R = real_pi.transpose() * r * real_pi;
  auto q = Eigen::Matrix4d{};
  for (auto a = 0; a != k_num_real_seq_letters; ++a) {
    q(a, a) = 0.0;
    for (auto b = 0; b != k_num_real_seq_letters; ++b) {
      if (a != b) {
        q(a, b) = r(a, b) / R * real_pi(b);
        q(a, a) -= q(a, b);
      }
    }
  };

  auto q_ab = Seq_matrix<double>{};
  
  q_ab[A][A] = q(0, 0);
  q_ab[A][C] = q(0, 1);
  q_ab[A][G] = q(0, 2);
  q_ab[A][T] = q(0, 3);
  q_ab[C][A] = q(1, 0);
  q_ab[C][C] = q(1, 1);
  q_ab[C][G] = q(1, 2);
  q_ab[C][T] = q(1, 3);
  q_ab[G][A] = q(2, 0);
  q_ab[G][C] = q(2, 1);
  q_ab[G][G] = q(2, 2);
  q_ab[G][T] = q(2, 3);
  q_ab[T][A] = q(3, 0);
  q_ab[T][C] = q(3, 1);
  q_ab[T][G] = q(3, 2);
  q_ab[T][T] = q(3, 3);

  return Site_evo_model{.mu = mu, .pi_a = pi_a, .q_ab = q_ab};
}

}  // namespace delphy
