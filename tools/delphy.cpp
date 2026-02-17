#include <numeric>
#include <iostream>

#include <boost/program_options.hpp>
#include <absl/log/initialize.h>

#include "beasty_output.h"
#include "delphy_output.h"
#include "dates.h"
#include "cmdline.h"
#include "run.h"
#include "phylo_tree_calc.h"

namespace po = boost::program_options;

namespace delphy {

struct Timestamp {
  int64_t step;
  std::chrono::time_point<std::chrono::high_resolution_clock> time;
};
static std::deque<Timestamp> timestamps;
inline constexpr size_t k_max_timestamps = 10;

static auto print_stats_line(const Run& run) -> void {
  using enum Real_seq_letter;

  const auto& tree = run.tree();
  
  auto M_ts = [&M_ab = run.num_muts_ab()]() {
    return M_ab[A][G] + M_ab[G][A] + M_ab[C][T] + M_ab[T][C];
  }();

  auto sum_nu = 0.0;
  auto sum_nu2 = 0.0;
  for (const auto& nu : run.nu()) {
    sum_nu += nu;
    sum_nu2 += nu * nu;
  }
  auto mean_nu = sum_nu / tree.num_sites();
  auto mean_nu2 = sum_nu2 / tree.num_sites();
  auto sigma_nu = std::sqrt(mean_nu2 - mean_nu*mean_nu);

  auto steps_per_s = 0.0;
  if (timestamps.size() > 1) {
    const auto& earliest_timestamp = timestamps.front();
    const auto& latest_timestamp = timestamps.back();
    steps_per_s = static_cast<double>(latest_timestamp.step - earliest_timestamp.step)
        / (1e-9 * static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(
            latest_timestamp.time - earliest_timestamp.time).count()));
  }

  auto log_G = run.log_G();
  auto log_coalescent_prior = run.log_coalescent_prior();
  auto log_other_priors = run.calc_cur_log_other_priors();
  auto log_posterior = log_G + log_coalescent_prior + log_other_priors;

  auto state_frequencies_at_root = run.state_frequencies_of_ref_sequence();
  for (const auto& m : tree.at_root().mutations) {
    --state_frequencies_at_root[m.from];
    ++state_frequencies_at_root[m.to];
  }
  for (const auto& [mi_site, mi_from] : tree.at_root().missations.slow_elements(tree.ref_sequence)) {
    --state_frequencies_at_root[mi_from];
  }
  
  std::cerr << absl::StreamFormat("%.1f M steps/s, ", steps_per_s / 1e6)
            << absl::StreamFormat("Step %d, ", run.step())
            << absl::StreamFormat("log_posterior = %.2f, ", log_posterior)
            << absl::StreamFormat("log_G = %.2f, ", run.log_G())
            << absl::StreamFormat("log_coal = %.2f, ", run.log_coalescent_prior())
            << absl::StreamFormat("log_other_priors = %.2f, ", log_other_priors)
            << absl::StreamFormat("num_muts = %d, ", run.num_muts())
            << absl::StreamFormat("T = %.2f, ", calc_T(tree))
            << absl::StreamFormat("t_MRCA = %s, ", to_iso_date(tree.at_root().t));
  const auto& pop_model = run.pop_model();
  if (typeid(pop_model) == typeid(Exp_pop_model)) {
    const auto& exp_pop_model = static_cast<const Exp_pop_model&>(run.pop_model());
    std::cerr << absl::StreamFormat("n0 = %.1f, ", exp_pop_model.pop_at_t0())
              << absl::StreamFormat("g = %.4f, ", exp_pop_model.growth_rate());
  }
  
  if (not run.mpox_hack_enabled()) {
    std::cerr << absl::StreamFormat("mu = %.2g * 10^-3 subst / site / year, ", run.mu() * 365 * 1000)
              << absl::StreamFormat("M_ts = %d, M_tv = %d, ", M_ts, run.num_muts() - M_ts)
              << absl::StreamFormat("kappa = %.3g, ", run.hky_kappa())
              << absl::StreamFormat("pi = [%.2g, %.2g, %.2g, %.2g], ",
                                    run.hky_pi()[A],
                                    run.hky_pi()[C],
                                    run.hky_pi()[G],
                                    run.hky_pi()[T]);
  } else {
    std::cerr << absl::StreamFormat("mpox_mu = %.2g * 10^-6 subst / site / year, ", run.mpox_mu()*365*1e6)
              << absl::StreamFormat("mpox_mu* = %.2g * 10^-4 subst / site / year, ", run.mpox_mu_star()*365*1e4);
  }

  std::cerr << (run.alpha_move_enabled()
                ? absl::StreamFormat("alpha = %.3g, ", run.alpha())
                : absl::StreamFormat("alpha OFF, "))
            << (run.alpha_move_enabled()
                ? absl::StreamFormat("nu ~ %.3g +/- %.3g, ", mean_nu, sigma_nu)
                : absl::StreamFormat("nu = 1.000, "))
            << absl::StreamFormat("Root Seq counts: [%d, %d, %d, %d], ",
                                  state_frequencies_at_root[A],
                                  state_frequencies_at_root[C],
                                  state_frequencies_at_root[G],
                                  state_frequencies_at_root[T])
            << absl::StreamFormat("num_muts_ab counts: [%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d]",
                                  run.num_muts_ab()[A][C],
                                  run.num_muts_ab()[A][G],
                                  run.num_muts_ab()[A][T],
                                  run.num_muts_ab()[C][A],
                                  run.num_muts_ab()[C][G],
                                  run.num_muts_ab()[C][T],
                                  run.num_muts_ab()[G][A],
                                  run.num_muts_ab()[G][C],
                                  run.num_muts_ab()[G][T],
                                  run.num_muts_ab()[T][A],
                                  run.num_muts_ab()[T][C],
                                  run.num_muts_ab()[T][G])
            << std::endl;
}

auto cli_main_loop(Processed_cmd_line& c) -> int {

  auto parallelism = c.thread_pool->size();
  std::cout << "# Parallelism: " << parallelism << "\n";
  c.run->set_num_parts(parallelism);

  auto beasty_log_output = std::unique_ptr<Beasty_log_output>{};
  if (c.log_filename.has_value()) {
    auto os_ptr = std::make_unique<std::ofstream>(c.log_filename.value());
    if (not *os_ptr) {
      throw std::runtime_error(absl::StrFormat("Could not open log file '%s' for writing", c.log_filename.value()));
    }
    beasty_log_output = std::make_unique<Beasty_log_output>(os_ptr.get(), c.beast_version, true);
    os_ptr.release();
    beasty_log_output->output_headers(*c.run);
  }
  auto beasty_trees_output = std::unique_ptr<Beasty_trees_output>{};
  if (c.trees_filename.has_value()) {
    auto os_ptr = std::make_unique<std::ofstream>(c.trees_filename.value());
    if (not *os_ptr) {
      throw std::runtime_error(absl::StrFormat("Could not open trees file '%s' for writing", c.trees_filename.value()));
    }
    beasty_trees_output = std::make_unique<Beasty_trees_output>(os_ptr.get(), c.beast_version, true);
    os_ptr.release();
    beasty_trees_output->output_headers(*c.run);
  }
  auto delphy_output = std::unique_ptr<Delphy_output>{};
  if (c.delphy_output_filename.has_value()) {
    auto os_ptr = std::make_unique<std::ofstream>(c.delphy_output_filename.value(), std::ios::binary);
    if (not *os_ptr) {
      throw std::runtime_error(absl::StrFormat("Could not open output file '%s' for writing",
                                               c.delphy_output_filename.value()));
    }
    delphy_output = std::make_unique<Delphy_output>(os_ptr.get(), true);
    os_ptr.release();
    if (c.delphy_output_metadata.has_value()) {
      delphy_output->set_dphy_metadata_blob(c.delphy_output_metadata.value());
    }
    delphy_output->output_preamble(*c.run, c.delphy_snapshot_every);
  }

  auto step_granularity = std::gcd(std::gcd(std::gcd(c.steps, c.log_every), c.tree_every), c.delphy_snapshot_every);
  while (true) {
    auto log_now = (c.run->step() % c.log_every) == 0;
    if (log_now) {
      print_stats_line(*c.run);
      if (beasty_log_output) {
        beasty_log_output->output_log(*c.run);
        beasty_log_output->flush();
      }
    }

    auto tree_now = (c.run->step() % c.tree_every) == 0;
    if (tree_now) {
      if (beasty_trees_output) {
        beasty_trees_output->output_tree(*c.run);
        beasty_trees_output->flush();
      }
    }

    auto delphy_output_now = (c.run->step() % c.delphy_snapshot_every) == 0;
    if (delphy_output_now) {
      if (delphy_output) {
        delphy_output->output_state(*c.run);
        delphy_output->flush();
      }
    }

    if (c.run->step() >= c.steps) {
      break;
    } else {
      c.run->do_mcmc_steps(step_granularity);
      
      timestamps.push_back(Timestamp{c.run->step(), std::chrono::high_resolution_clock::now()});
      if (timestamps.size() > k_max_timestamps) {
        timestamps.pop_front();
      }
    }
  }

  if (beasty_log_output) {
    beasty_log_output->output_footers(*c.run);
  }
  if (beasty_trees_output) {
    beasty_trees_output->output_footers(*c.run);
  }
  if (delphy_output) {
    delphy_output->output_epilog();
  }

  return 0;
}

}  // namespace delphy

auto main(int argc, char** argv) -> int {
  using namespace delphy;

  absl::InitializeLog();

  auto scope = Local_arena_scope{};

  auto c = process_args(argc, argv);

  return cli_main_loop(c);
}
