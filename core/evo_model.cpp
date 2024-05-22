#include "evo_model.h"

#include "absl/log/check.h"

namespace delphy {

auto make_single_partition_global_evo_model(Site_index num_sites) -> Global_evo_model {
  return Global_evo_model{
    Site_vector<Partition_index>(num_sites, 0),            // all sites in partition 0
    Site_vector<double>(num_sites, 1.0),                   // nu_l = 1.0 for all l
    Partition_vector<Site_evo_model>{{Site_evo_model{}}}   // Site_evo_model for the one and only partition
  };
}

auto make_global_evo_model(Site_vector<Partition_index> partition_for_site) -> Global_evo_model {
  CHECK(not partition_for_site.empty());
  auto num_sites = std::ssize(partition_for_site);
  auto num_partitions = 1 + std::ranges::max(partition_for_site);
  return {
    std::move(partition_for_site), 
    Site_vector<double>(num_sites, 1.0),         // nu_l = 1.0 for all l
    Partition_vector<Site_evo_model>(num_partitions, Site_evo_model{})
  };
}
Global_evo_model::Global_evo_model(
    Site_vector<Partition_index> partition_for_site,
    std::vector<double> nu_l,
    Partition_vector<Site_evo_model> partition_evo_model)
    : partition_for_site{std::move(partition_for_site)},
      nu_l{std::move(nu_l)},
      partition_evo_model{std::move(partition_evo_model)} {
}


}  // namespace delphy
