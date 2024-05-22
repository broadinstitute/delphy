#ifndef DELPHY_EVO_HKY_H_
#define DELPHY_EVO_HKY_H_

#include "evo_model.h"
#include "sequence.h"

namespace delphy {

struct Hky_model {
  double mu{0.0};
  double kappa{1.0};
  Seq_vector<double> pi_a{1.0 / k_num_real_seq_letters};

  auto derive_site_evo_model() const -> Site_evo_model;
};

}  // namespace delphy

#endif // DELPHY_EVO_HKY_H_
