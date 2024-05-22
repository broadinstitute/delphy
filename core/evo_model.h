#ifndef DELPHY_EVO_MODEL_H_
#define DELPHY_EVO_MODEL_H_

#include "sequence.h"

namespace delphy {

using Partition_index = int;
template<typename T>
using Partition_vector = std::vector<T>;

struct Site_evo_model {
  double mu{};
  Seq_vector<double> pi_a{};
  Seq_matrix<double> q_ab{};

  auto q_a(Real_seq_letter a) const -> double { return -q_ab[a][a]; }
};

struct Global_evo_model {
  Global_evo_model(Site_vector<Partition_index> partition_for_site,
                   std::vector<double> nu_l,
                   Partition_vector<Site_evo_model> partition_evo_model);
  
  Site_vector<Partition_index> partition_for_site;
  Site_vector<double> nu_l;
  Partition_vector<Site_evo_model> partition_evo_model;

  auto num_partitions() const -> Partition_index {
    return std::ssize(partition_evo_model);
  }
  
  auto mu_l(Site_index l) const -> double {
    return partition_evo_model[partition_for_site[l]].mu;
  }
  
  auto pi_l_a(Site_index l, Real_seq_letter a) const -> double {
    return partition_evo_model[partition_for_site[l]].pi_a[a];
  }

  auto q_l_a(Site_index l, Real_seq_letter a) const -> double {
    return partition_evo_model[partition_for_site[l]].q_a(a);
  }

  auto q_l_ab(Site_index l, Real_seq_letter a, Real_seq_letter b) const -> double {
    return partition_evo_model[partition_for_site[l]].q_ab[a][b];
  }
};

auto make_single_partition_global_evo_model(Site_index num_sites) -> Global_evo_model;
auto make_global_evo_model(Site_vector<Partition_index> partition_for_site) -> Global_evo_model;

}  // namespace delphy

#endif // DELPHY_EVO_MODEL_H_
