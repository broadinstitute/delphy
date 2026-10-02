#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <span>
#include <vector>

#include <absl/strings/str_join.h>

#include "estd.h"
#include "generic_tree_prober.h"

using namespace ::testing;

namespace delphy {

using estd::View_2d;

TEST(Generic_tree_prober_test, nan_probe_times) {
  auto pop = 0.2;   // N_e * rho
  auto pop_model = Const_pop_model{pop};
  auto num_bundles = 0;
  auto bundle_events = std::vector<Bundle_event>{};
  auto probe_times = std::vector<double>{0.0, 0.5, std::numeric_limits<double>::quiet_NaN(), 1.0};
  auto num_probe_times = static_cast<int>(std::ssize(probe_times));
  auto raw_results = std::vector<double>(num_bundles * num_probe_times, std::numeric_limits<double>::quiet_NaN());
  
  EXPECT_THROW(generic_probe_tree(bundle_events, num_bundles, pop_model, probe_times, raw_results), std::invalid_argument);
}

TEST(Generic_tree_prober_test, all_infinite_probe_times) {
  auto pop = 0.2;   // N_e * rho
  auto pop_model = Const_pop_model{pop};
  auto num_bundles = 0;
  auto bundle_events = std::vector<Bundle_event>{};
  auto probe_times = std::vector<double>{-std::numeric_limits<double>::infinity(), +std::numeric_limits<double>::infinity()};
  auto num_probe_times = static_cast<int>(std::ssize(probe_times));
  auto raw_results = std::vector<double>(num_bundles * num_probe_times, std::numeric_limits<double>::quiet_NaN());
  
  EXPECT_THROW(generic_probe_tree(bundle_events, num_bundles, pop_model, probe_times, raw_results), std::invalid_argument);
}

TEST(Generic_tree_prober_test, unsorted_probe_times) {
  auto pop = 0.2;   // N_e * rho
  auto pop_model = Const_pop_model{pop};
  auto num_bundles = 0;
  auto bundle_events = std::vector<Bundle_event>{};
  auto probe_times = std::vector<double>{1.0, 0.5, 0.0};
  auto num_probe_times = static_cast<int>(std::ssize(probe_times));
  auto raw_results = std::vector<double>(num_bundles * num_probe_times, std::numeric_limits<double>::quiet_NaN());
  
  EXPECT_THROW(generic_probe_tree(bundle_events, num_bundles, pop_model, probe_times, raw_results), std::invalid_argument);
}

TEST(Generic_tree_prober_test, negative_num_bundles) {
  auto pop = 0.2;   // N_e * rho
  auto pop_model = Const_pop_model{pop};
  auto num_bundles = -1;
  auto bundle_events = std::vector<Bundle_event>{};
  auto probe_times = std::vector<double>{0.0, 0.5, 1.0};
  auto raw_results = std::vector<double>{};  // Can't properly size results
  
  EXPECT_THROW(generic_probe_tree(bundle_events, num_bundles, pop_model, probe_times, raw_results), std::invalid_argument);
}

TEST(Generic_tree_prober_test, bundle_event_with_negative_target_bundle) {
  auto pop = 0.2;   // N_e * rho
  auto pop_model = Const_pop_model{pop};
  auto num_bundles = 1;
  auto bundle_events = std::vector<Bundle_event>{
    { .t = 0.0, .adding = true, .target_bundle = -1 }
  };
  auto probe_times = std::vector<double>{0.0, 0.5, 1.0};
  auto num_probe_times = static_cast<int>(std::ssize(probe_times));
  auto raw_results = std::vector<double>(num_bundles * num_probe_times, std::numeric_limits<double>::quiet_NaN());
  
  EXPECT_THROW(generic_probe_tree(bundle_events, num_bundles, pop_model, probe_times, raw_results), std::invalid_argument);
}

TEST(Generic_tree_prober_test, bundle_event_with_too_big_target_bundle) {
  auto pop = 0.2;   // N_e * rho
  auto pop_model = Const_pop_model{pop};
  auto num_bundles = 1;
  auto bundle_events = std::vector<Bundle_event>{
    { .t = 0.0, .adding = true, .target_bundle = 42 }
  };
  auto probe_times = std::vector<double>{0.0, 0.5, 1.0};
  auto num_probe_times = static_cast<int>(std::ssize(probe_times));
  auto raw_results = std::vector<double>(num_bundles * num_probe_times, std::numeric_limits<double>::quiet_NaN());
  
  EXPECT_THROW(generic_probe_tree(bundle_events, num_bundles, pop_model, probe_times, raw_results), std::invalid_argument);
}

TEST(Generic_tree_prober_test, bundle_event_with_NaN_time) {
  auto pop = 0.2;   // N_e * rho
  auto pop_model = Const_pop_model{pop};
  auto num_bundles = 1;
  auto bundle_events = std::vector<Bundle_event>{
    { .t = std::numeric_limits<double>::quiet_NaN(), .adding = true, .target_bundle = 0 }
  };
  auto probe_times = std::vector<double>{0.0, 0.5, 1.0};
  auto num_probe_times = static_cast<int>(std::ssize(probe_times));
  auto raw_results = std::vector<double>(num_bundles * num_probe_times, std::numeric_limits<double>::quiet_NaN());
  
  EXPECT_THROW(generic_probe_tree(bundle_events, num_bundles, pop_model, probe_times, raw_results), std::invalid_argument);
}

TEST(Generic_tree_prober_test, raw_results_too_small) {
  auto pop = 0.2;   // N_e * rho
  auto pop_model = Const_pop_model{pop};
  auto num_bundles = 1;
  auto bundle_events = std::vector<Bundle_event>{
    { .t = 0.0, .adding = true, .target_bundle = 0 }
  };
  auto probe_times = std::vector<double>{0.0, 0.5, 1.0};
  auto num_probe_times = static_cast<int>(std::ssize(probe_times));
  auto raw_results = std::vector<double>(num_bundles * num_probe_times - 1, std::numeric_limits<double>::quiet_NaN());  // -1 !!!
  
  EXPECT_THROW(generic_probe_tree(bundle_events, num_bundles, pop_model, probe_times, raw_results), std::invalid_argument);
}

TEST(Generic_tree_prober_test, raw_results_too_big) {
  auto pop = 0.2;   // N_e * rho
  auto pop_model = Const_pop_model{pop};
  auto num_bundles = 1;
  auto bundle_events = std::vector<Bundle_event>{
    { .t = 0.0, .adding = true, .target_bundle = 0 }
  };
  auto probe_times = std::vector<double>{0.0, 0.5, 1.0};
  auto num_probe_times = static_cast<int>(std::ssize(probe_times));
  auto raw_results = std::vector<double>(num_bundles * num_probe_times + 1, std::numeric_limits<double>::quiet_NaN());  // +1 !!!
  
  EXPECT_THROW(generic_probe_tree(bundle_events, num_bundles, pop_model, probe_times, raw_results), std::invalid_argument);
}

TEST(Generic_tree_prober_test, no_bundles) {
  auto pop = 0.2;   // N_e * rho
  auto pop_model = Const_pop_model{pop};
  auto num_bundles = 0;
  auto bundle_events = std::vector<Bundle_event>{};
  auto probe_times = std::vector<double>{0.0, 0.5, 1.0};
  auto num_probe_times = static_cast<int>(std::ssize(probe_times));
  auto raw_results = std::vector<double>(num_bundles * num_probe_times, std::numeric_limits<double>::quiet_NaN());
  
  EXPECT_NO_THROW(generic_probe_tree(bundle_events, num_bundles, pop_model, probe_times, raw_results));
}

TEST(Generic_tree_prober_test, no_probe_times) {
  auto pop = 0.2;   // N_e * rho
  auto pop_model = Const_pop_model{pop};
  auto num_bundles = 2;
  auto bundle_events = std::vector<Bundle_event>{
    { .t = 0.0, .adding = true, .target_bundle = 0 },
    { .t = 0.5, .adding = true, .target_bundle = 1 },
  };
  auto probe_times = std::vector<double>{};
  auto num_probe_times = static_cast<int>(std::ssize(probe_times));
  auto raw_results = std::vector<double>(num_bundles * num_probe_times, std::numeric_limits<double>::quiet_NaN());
  
  EXPECT_NO_THROW(generic_probe_tree(bundle_events, num_bundles, pop_model, probe_times, raw_results));
}

TEST(Generic_tree_prober_test, trivial) {
  // A single branch that extends to `-infty`
  auto pop = 0.2;   // N_e * rho
  auto pop_model = Const_pop_model{pop};
  auto num_bundles = 1;
  auto bundle_events = std::vector<Bundle_event>{
    { .t = -std::numeric_limits<double>::infinity(), .adding = true, .target_bundle = 0 },
  };
  auto probe_times = std::vector<double>{0.0};
  auto num_probe_times = static_cast<int>(std::ssize(probe_times));
  auto raw_results = std::vector<double>(num_bundles * num_probe_times, std::numeric_limits<double>::quiet_NaN());

  generic_probe_tree(bundle_events, num_bundles, pop_model, probe_times, raw_results);
  auto results = View_2d{raw_results, num_bundles, num_probe_times};

  EXPECT_THAT(results(0, 0), testing::DoubleNear(1.0, 1e-6));
}

TEST(Generic_tree_prober_test, simple_exact_const) {
  // One infinite branch from t = 0.0.  This has an exact solution even when the
  // population curve is discretized into cells:
  //
  //  p(t) = 1 - exp(-t / N)
  //
  auto pop = 2.0;   // N_e * rho
  auto pop_model = Const_pop_model{pop};
  auto num_bundles = 1;
  auto bundle_events = std::vector<Bundle_event>{
    { .t = 0.0, .adding = true, .target_bundle = 0 },
  };
  auto probe_times = std::vector<double>{0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  auto num_probe_times = static_cast<int>(std::ssize(probe_times));
  auto raw_results = std::vector<double>(num_bundles * num_probe_times, std::numeric_limits<double>::quiet_NaN());

  generic_probe_tree(bundle_events, num_bundles, pop_model, probe_times, raw_results);
  auto results = View_2d{raw_results, num_bundles, num_probe_times};

  for (auto i = 0; i != std::ssize(probe_times); ++i) {
    EXPECT_THAT(results(0, i), testing::DoubleNear(-std::expm1(-probe_times[i] / pop), 1e-6));
  }

  //std::cerr << "results[b=0]: [" << absl::StrJoin(results(0), ", ", absl::StreamFormatter()) << "]\n";
}

TEST(Generic_tree_prober_test, simple_exact_const_2) {
  // Two infinite branchs from t = 0.0.  This has an exact solution even when the
  // population curve is discretized into cells:
  //
  //  p_x(t) = (1/2) (1 - exp(-2 t / N))
  //
  auto pop = 2.0;   // N_e * rho
  auto pop_model = Const_pop_model{pop};
  auto num_bundles = 2;
  auto bundle_events = std::vector<Bundle_event>{
    { .t = 0.0, .adding = true, .target_bundle = 0 },
    { .t = 0.0, .adding = true, .target_bundle = 1 },
  };
  auto probe_times = std::vector<double>{0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  auto num_probe_times = static_cast<int>(std::ssize(probe_times));
  auto raw_results = std::vector<double>(num_bundles * num_probe_times, std::numeric_limits<double>::quiet_NaN());

  generic_probe_tree(bundle_events, num_bundles, pop_model, probe_times, raw_results);
  auto results = View_2d{raw_results, num_bundles, num_probe_times};

  for (auto i = 0; i != std::ssize(probe_times); ++i) {
    EXPECT_THAT(results(0, i), testing::DoubleNear(-0.5*std::expm1(-2 * probe_times[i] / pop), 1e-6));
    EXPECT_THAT(results(1, i), testing::DoubleNear(-0.5*std::expm1(-2 * probe_times[i] / pop), 1e-6));
  }

  //std::cerr << "results[0]: [" << absl::StrJoin(results(0), ", ", absl::StreamFormatter()) << "]\n";
  //std::cerr << "results[1]: [" << absl::StrJoin(results(1), ", ", absl::StreamFormatter()) << "]\n";
}

TEST(Generic_tree_prober_test, simple_exact_exp) {
  // One infinite branch from t = 0.0.  This has an exact solution even when the
  // population curve is discretized into cells:
  //
  //  p(t) = 1 - exp(-int_0^t dt' / N(t'))
  //
  auto t_0 = 0.0;   // N_e(t) * rho = n_0 * exp(g * (t - t_0))
  auto n_0 = 1.0;   // N_e * rho at t = t_0
  auto g = 0.5;     // e-foldings per unit time
  auto pop_model = Exp_pop_model{t_0, n_0, g, 1e-10};
  auto num_bundles = 1;
  auto bundle_events = std::vector<Bundle_event>{
    { .t = 0.0, .adding = true, .target_bundle = 0 },
  };
  auto probe_times = std::vector<double>{0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  auto num_probe_times = static_cast<int>(std::ssize(probe_times));
  auto raw_results = std::vector<double>(num_bundles * num_probe_times, std::numeric_limits<double>::quiet_NaN());

  generic_probe_tree(bundle_events, num_bundles, pop_model, probe_times, raw_results);
  auto results = View_2d{raw_results, num_bundles, num_probe_times};

  for (auto i = 0; i != std::ssize(probe_times); ++i) {
    EXPECT_THAT(results(0, i), testing::DoubleNear(-std::expm1(-pop_model.intensity_integral(0.0, probe_times[i])), 1e-6));
  }

  //std::cerr << "results[0]: [" << absl::StrJoin(results(0), ", ", absl::StreamFormatter()) << "]\n";
}

TEST(Generic_tree_prober_test, simple) {
  auto pop = 0.2;   // N_e * rho
  auto pop_model = Const_pop_model{pop};
  auto num_bundles = 2;
  auto bundle_events = std::vector<Bundle_event>{
    // Two simple branches
    { .t = 1.0, .adding = true, .target_bundle = 0 },
    { .t = 5.0, .adding = false, .target_bundle = 0 },
    
    { .t = 3.0, .adding = true, .target_bundle = 1 },
    { .t = 10.0, .adding = false, .target_bundle = 1 },
  };
  auto probe_times = std::vector<double>{0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  auto num_probe_times = static_cast<int>(std::ssize(probe_times));
  auto raw_results = std::vector<double>(num_bundles * num_probe_times, std::numeric_limits<double>::quiet_NaN());

  generic_probe_tree(bundle_events, num_bundles, pop_model, probe_times, raw_results);
  auto results = View_2d{raw_results, num_bundles, num_probe_times};

  EXPECT_THAT(results(0, 0), testing::DoubleNear(0.0, 1e-2));  // t = 0.0
  EXPECT_THAT(results(0, 3), testing::DoubleNear(1.0, 1e-2));  // t = 3.0
  EXPECT_THAT(results(0, 5), testing::Gt(results(1, 5)));      // t = 5.0
  EXPECT_THAT(results(0, 10), testing::DoubleNear(0.0, 1e-2)); // t = 10.0
  
  EXPECT_THAT(results(1, 0), testing::DoubleNear(0.0, 1e-2));  // t = 0.0
  EXPECT_THAT(results(1, 3), testing::DoubleNear(0.0, 1e-2));  // t = 3.0
  EXPECT_THAT(results(1, 5), testing::Gt(0.4));                // t = 5.0
  EXPECT_THAT(results(1, 10), testing::DoubleNear(1.0, 1e-2)); // t = 10.0

  //std::cerr << "results[0]: [" << absl::StrJoin(results(0), ", ", absl::StreamFormatter()) << "]\n";
  //std::cerr << "results[1]: [" << absl::StrJoin(results(1), ", ", absl::StreamFormatter()) << "]\n";

  for (auto i = 4; i != std::ssize(probe_times); ++i) {  // Before t = 4.0, the probe might go past the root and not attach
    auto p_tot = 0.0;
    for (auto b = 0; b != num_bundles; ++b) {
      p_tot += results(b, i);
    }
    EXPECT_THAT(p_tot, testing::DoubleNear(1.0, 1e-6)) << i;
  }
}

TEST(Generic_tree_prober_test, simple_tree_like) {
  // A branch from t = 0.0 that splits into two at t = t_c.  This has an exact solution even when the
  // population curve is discretized into cells:
  //
  //  p(t) = 1 - exp(-int_0^t k(t') / N(t') dt')
  //
  //         /  1 - exp(-t / N),                   0 <= t <= t_c;
  //       = |
  //         \  1 - exp(-(t_c + 2*(t-t_c)) / N), t_c <= t.
  //
  auto pop = 0.5;   // N_e * rho
  auto pop_model = Const_pop_model{pop};
  auto num_bundles = 1;
  auto t_c = 1.0;
  auto bundle_events = std::vector<Bundle_event>{
    // A branch above the root that splits at t = t_c.
    { .t = 0.0, .adding = true, .target_bundle = 0 },
    { .t = t_c, .adding = false, .target_bundle = 0 },  // end of parent branch
    { .t = t_c, .adding = true, .target_bundle = 0 },   // start of child 1
    { .t = t_c, .adding = true, .target_bundle = 0 },   // start of child 2
  };
  auto probe_times = std::vector<double>{0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0};
  auto num_probe_times = static_cast<int>(std::ssize(probe_times));
  auto raw_results = std::vector<double>(num_bundles * num_probe_times, std::numeric_limits<double>::quiet_NaN());

  generic_probe_tree(bundle_events, num_bundles, pop_model, probe_times, raw_results);
  auto results = View_2d{raw_results, num_bundles, num_probe_times};

  for (auto i = 0; i != std::ssize(probe_times); ++i) {
    auto t = probe_times[i];
    auto expected_p_t = t < t_c
        ? -std::expm1(-t / pop)
        : -std::expm1(-(t_c + 2 *(t - t_c)) / pop);
    EXPECT_THAT(results(0, i), testing::DoubleNear(expected_p_t, 1e-6));
  }

  //std::cerr << "results[0]: [" << absl::StrJoin(results(0), ", ", absl::StreamFormatter()) << "]\n";
}

TEST(Generic_tree_prober_test, simple_tree_like_two_bundles) {
  // A branch from t = 0.0 in one bundle that splits into two at t = t_c in a different bundle.
  // This has an exact solution even when the population curve is discretized into cells:
  //
  //           / 0                                        t in (-inf, 0]
  //  p_0(t) = | 1 - exp(-t / N),                         t in [0, t_c]
  //           \ (1 - exp(-t_c/N)) * exp(-2*(t-t_c)/N),   t in [t_c, +inf)
  //
  //           / 0                                        t in (-inf, t_c]
  //  p_1(t) = |
  //           \ 1 - exp(-2*(t-t_c)/N)                    t in [t_c, +inf)
  //
  auto pop = 0.5;   // N_e * rho
  auto pop_model = Const_pop_model{pop};
  auto num_bundles = 2;
  auto t_c = 1.0;
  auto bundle_events = std::vector<Bundle_event>{
    // A branch above the root that splits at t = t_c.
    { .t = 0.0, .adding = true, .target_bundle = 0 },
    { .t = t_c, .adding = false, .target_bundle = 0 },  // end of parent branch
    { .t = t_c, .adding = true, .target_bundle = 1 },   // start of child 1
    { .t = t_c, .adding = true, .target_bundle = 1 },   // start of child 2
  };
  auto probe_times = std::vector<double>{0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0};
  auto num_probe_times = static_cast<int>(std::ssize(probe_times));
  auto raw_results = std::vector<double>(num_bundles * num_probe_times, std::numeric_limits<double>::quiet_NaN());

  generic_probe_tree(bundle_events, num_bundles, pop_model, probe_times, raw_results);
  auto results = View_2d{raw_results, num_bundles, num_probe_times};

  for (auto i = 0; i != std::ssize(probe_times); ++i) {
    auto t = probe_times[i];
    auto expected_p_t = t < t_c
        ? -std::expm1(-t / pop)
        : (-std::expm1(-t_c/pop)) * std::exp(- 2*(t - t_c) / pop);
    EXPECT_THAT(results(0, i), testing::DoubleNear(expected_p_t, 1e-6));
  }
  
  for (auto i = 0; i != std::ssize(probe_times); ++i) {
    auto t = probe_times[i];
    auto expected_p_t = t < t_c
        ? 0
        : -std::expm1(-2*(t-t_c) / pop);
    EXPECT_THAT(results(1, i), testing::DoubleNear(expected_p_t, 1e-6));
  }

  //std::cerr << "results[0]: [" << absl::StrJoin(results(0), ", ", absl::StreamFormatter()) << "]\n";
  //std::cerr << "results[1]: [" << absl::StrJoin(results(1), ", ", absl::StreamFormatter()) << "]\n";
}

// Brute-force reference solution for a constant population N = `pop`.  Splitting (-infty, t]
// into intervals j = 0, 1, ... between consecutive event times (j = 0 is the most recent), we
// evaluate directly
//
//   p_b(t) = sum_j (k_{b,j} / k_j) (1 - exp(-k_j Delta_j / N)) prod_{i<j} exp(-k_i Delta_i / N).
//
// The branch counts in each interval are obtained by summing all adds and removes before the
// interval's midpoint, so they don't depend on how `generic_probe_tree` orders equal-time events.
static auto reference_probe_const_pop(
    std::span<const Bundle_event> bundle_events,
    int num_bundles,
    double pop,
    double t)
    -> std::vector<double> {
  auto result = std::vector<double>(num_bundles, 0.0);

  // Interval boundaries: t itself and all event times before t, most recent first
  auto boundaries = std::vector<double>{t};
  for (const auto& e : bundle_events) {
    if (e.t < t) {
      boundaries.push_back(e.t);
    }
  }
  std::ranges::sort(boundaries, std::ranges::greater{});
  auto [first_dup, last_dup] = std::ranges::unique(boundaries);
  boundaries.erase(first_dup, last_dup);

  auto p_survive = 1.0;  // prob. of not attaching in any interval more recent than the current one
  for (auto j = 0; j + 1 < std::ssize(boundaries); ++j) {
    auto t_hi = boundaries[j];
    auto t_lo = boundaries[j + 1];
    auto t_mid =
        std::isinf(t_lo) && std::isinf(t_hi) ? 0.0 :
        std::isinf(t_lo)                     ? t_hi - 1.0 :
        std::isinf(t_hi)                     ? t_lo + 1.0 :
        0.5 * (t_lo + t_hi);

    auto k_b = std::vector<int>(num_bundles, 0);
    for (const auto& e : bundle_events) {
      if (e.t < t_mid) {
        k_b[e.target_bundle] += e.adding ? +1 : -1;
      }
    }
    auto k = 0;
    for (auto b = 0; b != num_bundles; ++b) {
      k += k_b[b];
    }
    if (k == 0) {
      continue;  // Nothing to attach to, and survival factor is 1
    }

    auto p_attach = -std::expm1(-k * (t_hi - t_lo) / pop);
    for (auto b = 0; b != num_bundles; ++b) {
      result[b] += p_survive * p_attach * k_b[b] / k;
    }
    p_survive *= 1.0 - p_attach;
  }

  return result;
}

// Compare `generic_probe_tree` against `reference_probe_const_pop` at every probe time
static auto expect_matches_reference_const_pop(
    std::span<const double> raw_results,
    std::span<const Bundle_event> bundle_events,
    int num_bundles,
    double pop,
    std::span<const double> probe_times)
    -> void {
  auto num_probe_times = static_cast<int>(std::ssize(probe_times));
  auto results = View_2d{raw_results, num_bundles, num_probe_times};
  for (auto i = 0; i != num_probe_times; ++i) {
    auto expected = reference_probe_const_pop(bundle_events, num_bundles, pop, probe_times[i]);
    for (auto b = 0; b != num_bundles; ++b) {
      EXPECT_THAT(results(b, i), testing::DoubleNear(expected[b], 1e-10))
          << "b = " << b << ", t = " << probe_times[i];
    }
  }
}

TEST(Generic_tree_prober_test, bundle_reactivated_after_probe) {
  // Bundle 0 has two disjoint branches, [0, 1] and [2, +inf); bundle 1 has one branch, [0.5, +inf).
  // The probe at t = 1.5 deactivates bundle 0 (k_0 = 0 at that point), which moves bundle 1 into
  // a different slot in the active list.  Bundle 0 is then reactivated in a new slot at t = 2.
  //
  // With N = 0.5 (coalescence rate 2 per branch), a probe at t = 2.5 sees, going backwards:
  //   [2, 2.5]:  k = 2 (one each),  k * Delta / N = 2
  //   [1, 2]:    k = 1 (bundle 1),  k * Delta / N = 2
  //   [0.5, 1]:  k = 2 (one each),  k * Delta / N = 2
  //   [0, 0.5]:  k = 1 (bundle 0),  k * Delta / N = 1
  // so
  //   p_0(2.5) = (1/2)(1 - e^-2) + e^-4 (1/2)(1 - e^-2) + e^-6 (1 - e^-1)
  //   p_1(2.5) = (1/2)(1 - e^-2) + e^-2 (1 - e^-2) + e^-4 (1/2)(1 - e^-2)
  auto pop = 0.5;   // N_e * rho
  auto pop_model = Const_pop_model{pop};
  auto num_bundles = 2;
  auto bundle_events = std::vector<Bundle_event>{
    { .t = 0.0, .adding = true,  .target_bundle = 0 },
    { .t = 1.0, .adding = false, .target_bundle = 0 },
    { .t = 2.0, .adding = true,  .target_bundle = 0 },
    { .t = 0.5, .adding = true,  .target_bundle = 1 },
  };
  auto probe_times = std::vector<double>{0.25, 0.75, 1.5, 2.5, 3.0};
  auto num_probe_times = static_cast<int>(std::ssize(probe_times));
  auto raw_results = std::vector<double>(num_bundles * num_probe_times, std::numeric_limits<double>::quiet_NaN());

  generic_probe_tree(bundle_events, num_bundles, pop_model, probe_times, raw_results);
  auto results = View_2d{raw_results, num_bundles, num_probe_times};

  // Sanity-check the reference solution itself against the hand calculation above
  auto em1 = std::exp(-1.0);
  auto em2 = std::exp(-2.0);
  auto expected_p0_2_5 = 0.5*(1 - em2) + em2*em2 * 0.5*(1 - em2) + em2*em2*em2 * (1 - em1);
  auto expected_p1_2_5 = 0.5*(1 - em2) + em2 * (1 - em2) + em2*em2 * 0.5*(1 - em2);
  auto reference_2_5 = reference_probe_const_pop(bundle_events, num_bundles, pop, 2.5);
  EXPECT_THAT(reference_2_5[0], testing::DoubleNear(expected_p0_2_5, 1e-12));
  EXPECT_THAT(reference_2_5[1], testing::DoubleNear(expected_p1_2_5, 1e-12));

  // Now the real test
  expect_matches_reference_const_pop(raw_results, bundle_events, num_bundles, pop, probe_times);
  EXPECT_THAT(results(0, 3), testing::DoubleNear(expected_p0_2_5, 1e-10));
  EXPECT_THAT(results(1, 3), testing::DoubleNear(expected_p1_2_5, 1e-10));
}

TEST(Generic_tree_prober_test, bundle_readded_before_probe) {
  // Bundle 0 has two disjoint branches, [0, 1] and [1.5, +inf); bundle 1 has one branch, [0.5, +inf).
  // No probe falls in (1, 1.5), so bundle 0 keeps its (k_0 = 0) slot in the active list and
  // it is reused when bundle 0 is re-added at t = 1.5.
  auto pop = 0.5;   // N_e * rho
  auto pop_model = Const_pop_model{pop};
  auto num_bundles = 2;
  auto bundle_events = std::vector<Bundle_event>{
    { .t = 0.0, .adding = true,  .target_bundle = 0 },
    { .t = 1.0, .adding = false, .target_bundle = 0 },
    { .t = 1.5, .adding = true,  .target_bundle = 0 },
    { .t = 0.5, .adding = true,  .target_bundle = 1 },
  };
  auto probe_times = std::vector<double>{0.25, 2.0, 3.0};
  auto num_probe_times = static_cast<int>(std::ssize(probe_times));
  auto raw_results = std::vector<double>(num_bundles * num_probe_times, std::numeric_limits<double>::quiet_NaN());

  generic_probe_tree(bundle_events, num_bundles, pop_model, probe_times, raw_results);

  expect_matches_reference_const_pop(raw_results, bundle_events, num_bundles, pop, probe_times);
}

TEST(Generic_tree_prober_test, equal_probe_times) {
  // Same tree as `simple_tree_like_two_bundles`, with repeated probe times, including some
  // at exactly the time t_c when the parent branch ends and the children begin
  auto pop = 0.5;   // N_e * rho
  auto pop_model = Const_pop_model{pop};
  auto num_bundles = 2;
  auto t_c = 1.0;
  auto bundle_events = std::vector<Bundle_event>{
    { .t = 0.0, .adding = true,  .target_bundle = 0 },
    { .t = t_c, .adding = false, .target_bundle = 0 },  // end of parent branch
    { .t = t_c, .adding = true,  .target_bundle = 1 },  // start of child 1
    { .t = t_c, .adding = true,  .target_bundle = 1 },  // start of child 2
  };
  auto probe_times = std::vector<double>{0.5, 0.5, t_c, t_c, t_c, 1.5, 1.5, 2.0};
  auto num_probe_times = static_cast<int>(std::ssize(probe_times));
  auto raw_results = std::vector<double>(num_bundles * num_probe_times, std::numeric_limits<double>::quiet_NaN());

  generic_probe_tree(bundle_events, num_bundles, pop_model, probe_times, raw_results);
  auto results = View_2d{raw_results, num_bundles, num_probe_times};

  expect_matches_reference_const_pop(raw_results, bundle_events, num_bundles, pop, probe_times);
  for (auto b = 0; b != num_bundles; ++b) {
    EXPECT_EQ(results(b, 0), results(b, 1)) << "b = " << b;
    EXPECT_EQ(results(b, 2), results(b, 3)) << "b = " << b;
    EXPECT_EQ(results(b, 3), results(b, 4)) << "b = " << b;
    EXPECT_EQ(results(b, 5), results(b, 6)) << "b = " << b;
  }
}

TEST(Generic_tree_prober_test, events_after_last_probe_are_irrelevant) {
  // Events after the last probe time can't affect any result, even with a time-varying population
  // (the population cells are laid out only over the times up to the last probe)
  auto t_0 = 0.0;   // N_e(t) * rho = n_0 * exp(g * (t - t_0))
  auto n_0 = 1.0;   // N_e * rho at t = t_0
  auto g = 0.5;     // e-foldings per unit time
  auto pop_model = Exp_pop_model{t_0, n_0, g, 1e-10};
  auto num_bundles = 2;
  auto bundle_events_without_extras = std::vector<Bundle_event>{
    { .t = 0.0, .adding = true, .target_bundle = 0 },
    { .t = 0.5, .adding = true, .target_bundle = 1 },
  };
  auto bundle_events_with_extras = std::vector<Bundle_event>{
    { .t = 0.0, .adding = true,  .target_bundle = 0 },
    { .t = 0.5, .adding = true,  .target_bundle = 1 },
    { .t = 3.0, .adding = false, .target_bundle = 0 },  // Extra: ends bundle 0's branch
    { .t = 3.5, .adding = true,  .target_bundle = 1 },  // Extra: second branch in bundle 1...
    { .t = 4.0, .adding = false, .target_bundle = 1 },  // Extra: ...which ends later
    { .t = 5.0, .adding = false, .target_bundle = 1 },  // Extra: ends bundle 1's first branch
  };
  auto probe_times = std::vector<double>{0.25, 1.0, 2.0};
  auto num_probe_times = static_cast<int>(std::ssize(probe_times));
  auto raw_results_without_extras = std::vector<double>(num_bundles * num_probe_times, std::numeric_limits<double>::quiet_NaN());
  auto raw_results_with_extras = std::vector<double>(num_bundles * num_probe_times, std::numeric_limits<double>::quiet_NaN());

  generic_probe_tree(bundle_events_without_extras, num_bundles, pop_model, probe_times, raw_results_without_extras);
  generic_probe_tree(bundle_events_with_extras, num_bundles, pop_model, probe_times, raw_results_with_extras);

  EXPECT_EQ(raw_results_with_extras, raw_results_without_extras);  // Exact equality
}

TEST(Generic_tree_prober_test, infinite_probe_times) {
  // Bundle 0: branch [0, 1]; bundle 1: branch [0.5, +inf); bundle 2: branch [1, +inf).
  // A probe at -inf can't attach anywhere.  A probe at +inf is guaranteed to attach somewhere
  // in [1, +inf), where bundles 1 and 2 each have one branch.
  auto pop = 0.5;   // N_e * rho
  auto pop_model = Const_pop_model{pop};
  auto num_bundles = 3;
  auto bundle_events = std::vector<Bundle_event>{
    { .t = 0.0, .adding = true,  .target_bundle = 0 },
    { .t = 1.0, .adding = false, .target_bundle = 0 },
    { .t = 0.5, .adding = true,  .target_bundle = 1 },
    { .t = 1.0, .adding = true,  .target_bundle = 2 },
  };
  auto inf = std::numeric_limits<double>::infinity();
  auto probe_times = std::vector<double>{-inf, 0.5, 1.5, +inf};
  auto num_probe_times = static_cast<int>(std::ssize(probe_times));
  auto raw_results = std::vector<double>(num_bundles * num_probe_times, std::numeric_limits<double>::quiet_NaN());

  generic_probe_tree(bundle_events, num_bundles, pop_model, probe_times, raw_results);
  auto results = View_2d{raw_results, num_bundles, num_probe_times};

  expect_matches_reference_const_pop(raw_results, bundle_events, num_bundles, pop, probe_times);

  EXPECT_EQ(results(0, 0), 0.0);
  EXPECT_EQ(results(1, 0), 0.0);
  EXPECT_EQ(results(2, 0), 0.0);

  EXPECT_THAT(results(0, 3), testing::DoubleNear(0.0, 1e-12));
  EXPECT_THAT(results(1, 3), testing::DoubleNear(0.5, 1e-12));
  EXPECT_THAT(results(2, 3), testing::DoubleNear(0.5, 1e-12));
}

#ifndef __EMSCRIPTEN__  // We want tests to be compiled to exercise clang, but EXPECT_DEATH doesn't exist in Emscripten
TEST(Generic_tree_prober_test, remove_without_add_dies) {
  auto pop_model = Const_pop_model{0.5};
  auto num_bundles = 1;
  auto bundle_events = std::vector<Bundle_event>{
    { .t = 1.0, .adding = false, .target_bundle = 0 },
  };
  auto probe_times = std::vector<double>{2.0};
  auto num_probe_times = static_cast<int>(std::ssize(probe_times));
  auto raw_results = std::vector<double>(num_bundles * num_probe_times, std::numeric_limits<double>::quiet_NaN());

  EXPECT_DEATH(generic_probe_tree(bundle_events, num_bundles, pop_model, probe_times, raw_results), "");
}

TEST(Generic_tree_prober_test, remove_before_add_dies) {
  auto pop_model = Const_pop_model{0.5};
  auto num_bundles = 1;
  auto bundle_events = std::vector<Bundle_event>{
    { .t = 1.0, .adding = true,  .target_bundle = 0 },
    { .t = 0.5, .adding = false, .target_bundle = 0 },
  };
  auto probe_times = std::vector<double>{2.0};
  auto num_probe_times = static_cast<int>(std::ssize(probe_times));
  auto raw_results = std::vector<double>(num_bundles * num_probe_times, std::numeric_limits<double>::quiet_NaN());

  EXPECT_DEATH(generic_probe_tree(bundle_events, num_bundles, pop_model, probe_times, raw_results), "");
}

TEST(Generic_tree_prober_test, double_remove_dies) {
  // No probe between the two removes, so bundle 0 is still in the active list with k_0 = 0
  // when the second remove arrives
  auto pop_model = Const_pop_model{0.5};
  auto num_bundles = 1;
  auto bundle_events = std::vector<Bundle_event>{
    { .t = 0.0, .adding = true,  .target_bundle = 0 },
    { .t = 1.0, .adding = false, .target_bundle = 0 },
    { .t = 1.5, .adding = false, .target_bundle = 0 },
  };
  auto probe_times = std::vector<double>{2.0};
  auto num_probe_times = static_cast<int>(std::ssize(probe_times));
  auto raw_results = std::vector<double>(num_bundles * num_probe_times, std::numeric_limits<double>::quiet_NaN());

  EXPECT_DEATH(generic_probe_tree(bundle_events, num_bundles, pop_model, probe_times, raw_results), "");
}
#endif

}  // namespace delphy
