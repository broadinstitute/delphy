#ifndef DELPHY_SEQ_VECTOR_MATCHERS_H_
#define DELPHY_SEQ_VECTOR_MATCHERS_H_

#include <gtest/gtest.h>
#include <gmock/gmock.h>

namespace delphy {

MATCHER_P2(seq_vector_double_near, expected, tolerance, "") {
  auto result = true;
  for (auto a : k_all_real_seq_letters) {
    result = result && ExplainMatchResult(testing::DoubleNear(expected[a], tolerance), arg[a], result_listener);
  }
  return result;
}

MATCHER_P(pointwise_seq_vector_double_near, tolerance, "") {
  const auto& [v_actual, v_expected] = arg;
  return ExplainMatchResult(seq_vector_double_near(v_expected, tolerance), v_actual, result_listener);
}

MATCHER_P(pointwise_double_near, tolerance, "") {
  const auto& [v_actual, v_expected] = arg;
  return ExplainMatchResult(testing::DoubleNear(v_expected, tolerance), v_actual, result_listener);
}
  
}  // namespace delphy

#endif // DELPHY_SEQ_VECTOR_MATCHERS_H_
