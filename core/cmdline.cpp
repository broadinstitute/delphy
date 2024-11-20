#include "cmdline.h"

#include <iostream>
#include <fstream>
#include <optional>

#include "absl/random/random.h"
#include "cxxopts.hpp"

#include "dates.h"
#include "io.h"
#include "phylo_tree.h"
#include "sequence_utils.h"
#include "version.h"
#include "beasty_input.h"

namespace delphy {

bool delphy_invoked_via_cli{false};
std::vector<std::string> delphy_cli_args{};

auto build_rough_initial_tree_from_fasta(
    const std::vector<Fasta_entry>& in_fasta,
    bool random,
    absl::BitGenRef bitgen,
    const std::function<void(int,int)>& progress_hook,
    const std::function<void(const std::string&)>& warning_hook)
    -> Phylo_tree {

  if (in_fasta.empty()) {
    throw std::runtime_error("FASTA file has no sequences!");
  }
  
  // Build out "consensus" sequence from all samples, which we'll use to initially impute missing data
  // (helps prevent things like a site that's missing in almost all sequences being imputed as something
  // improbable, and then local topo/resampling moves have a very hard time untangling that error)
  auto consensus_sequence = deduce_consensus_sequence(
      in_fasta | std::views::transform([](const auto& e) { return e.sequence; }));

  // Arbitrary pick the consensus sequence as a reference
  auto resolved_ref_seq = consensus_sequence;

  // Set up tips by converting every sequence into a list of mutations from the reference
  // (essentially, convert the FASTA to a MAPLE file)
  auto tip_descs = std::vector<Tip_desc>{};
  for (auto& fasta_entry : in_fasta) {
    // Ignore sequences that we can't date correctly
    if (auto opt_t = extract_date_from_sequence_id(fasta_entry.id); not opt_t.has_value()) {
      warning_hook(absl::StrFormat("Ignoring sequence '%s', could not determine its date\n", fasta_entry.id));
    } else {
      auto delta = calculate_delta_from_reference(
          fasta_entry.sequence, resolved_ref_seq,
          [&fasta_entry, &warning_hook](const std::string& msg) { warning_hook(
              absl::StrFormat("In sequence '%s', %s", fasta_entry.id, msg)); });
      tip_descs.push_back({
          .name = fasta_entry.id,
          .t = opt_t.value(),
          .seq_deltas = std::move(delta.seq_deltas),
          .missations = std::move(delta.missations)});
    }
  }

  // Join all the tips up in a very rough approximation to greedy parsimony
  if (random) {
    return build_random_tree(std::move(resolved_ref_seq), std::move(tip_descs), bitgen, progress_hook);
  } else {
    return build_usher_like_tree(std::move(resolved_ref_seq), std::move(tip_descs), bitgen, progress_hook);
  }
}

auto build_rough_initial_tree_from_maple(
    Maple_file&& in_maple,
    bool random,
    absl::BitGenRef bitgen)
    -> Phylo_tree {

  if (in_maple.tip_descs.empty()) {
    std::cerr << "MAPLE file has no sequences!\n";
    std::exit(EXIT_FAILURE);
  }
  
  // Join all the tips up in a very rough approximation to greedy parsimony
  if (random) {
    return build_random_tree(std::move(in_maple.ref_sequence), std::move(in_maple.tip_descs), bitgen);
  } else {
    return build_usher_like_tree(std::move(in_maple.ref_sequence), std::move(in_maple.tip_descs), bitgen);
  }
}

auto process_args(int argc, char** argv) -> Processed_cmd_line {
  
  // Before anything else, record CLI invocation (e.g., for log files)
  delphy_invoked_via_cli = true;
  delphy_cli_args = std::vector<std::string>{argv, argv+argc};
  
  cxxopts::Options options("delphy", "Delphy - Fast and scalable Bayesian phylogenetic tree inference");

  options.add_options("Generic options")
      ("version", "print version string")
      ("h,help", "print usage")
      ("n,dry-run", "Dry run: process inputs, generate BEAST XMLs (if needed) but don't proceed with run")
      ;

  options.add_options("Version-0 options: unstable, temporary syntax for specifying Delphy's behaviour")
      ("v0-threads", "Number of threads to use (default: use all cores)",
       cxxopts::value<int>())
      ("v0-seed", "Initial random number seed (default: random)",
       cxxopts::value<uint32_t>())
      ("v0-in-fasta", "input FASTA file", cxxopts::value<std::string>())
      ("v0-in-maple", "input MAPLE file", cxxopts::value<std::string>())
      ("v0-init-heuristic", "Build initial tree with rough heuristics instead of randomly (default)",
       cxxopts::value<bool>())
      ("v0-init-random", "Build initial tree with random",
       cxxopts::value<bool>())
      ("v0-steps", "Total number of steps to run (default: ~100,000 per tip)",
       cxxopts::value<int64_t>())
      ("v0-out-log-file", "Filename for BEAST-like log file",
       cxxopts::value<std::string>())
      ("v0-log-every", "Steps between entries in screen & file logs (default: steps / 100)",
       cxxopts::value<int64_t>())
      ("v0-out-trees-file", "Filename for BEAST-like trees file",
       cxxopts::value<std::string>())
      ("v0-tree-every", "Steps between tree snapshots in tree log (default: steps / 100)",
       cxxopts::value<int64_t>())
      ("v0-out-delphy-file", "Filename for .dphy run file",
       cxxopts::value<std::string>())
      ("v0-delphy-snapshot-every", "Steps between run snapshots in .dphy run file (default: steps / 100)",
       cxxopts::value<int64_t>())
      ("v0-site-rate-heterogeneity", "Enable site rate heterogeneity",
       cxxopts::value<bool>()->default_value("false"))
      ("v0-fix-mutation-rate", "Fix mutation rate",
       cxxopts::value<bool>()->default_value("false"))
      ("v0-init-mutation-rate", "Initial (or fixed) value of mutation rate, in subst / site / year",
       cxxopts::value<double>()->default_value("1e-3"))
      ("v0-out-beast-xml", "Filename for XML input file suitable for BEAST 2.6.2",
       cxxopts::value<std::string>())
      ("v0-mpox-hack", "Enable mpox hack (very crude treatment of APOBEC3 mechanism)",
       cxxopts::value<bool>()->default_value("false"))
      ("v0-target-coal-prior-cells", "Target number of cells to use in parallelized coalescent prior (coalescent prior resolution is adjusted if actual number is more than 33% away from target); higher is more accurate but more expensive",
       cxxopts::value<int>()->default_value("400"))
      ("v0-fix-final-pop-size", "Fix effective population size at time of last tip",
       cxxopts::value<bool>()->default_value("false"))
      ("v0-init-final-pop-size", "Initial (or fixed) value of effective population size at time of last time, in years (the effective population size includes a factor of the generation time)",
       cxxopts::value<double>()->default_value("3.0"))
      ("v0-fix-pop-growth-rate", "Fix effective population size growth rate",
       cxxopts::value<bool>()->default_value("false"))
      ("v0-init-pop-growth-rate", "Initial (or fixed) value of effective population size growth rate, in e-foldings / year",
       cxxopts::value<double>()->default_value("0.0"))
      ;

  try {
    auto opts = options.parse(argc, argv);

    if (opts.count("version")) {
      std::cout << absl::StreamFormat("Delphy Version %s (build %d, commit %s)",
                                      k_delphy_version_string,
                                      k_delphy_build_number,
                                      k_delphy_commit_string) << "\n";
      std::exit(EXIT_SUCCESS);
    }
    if (opts.count("help")) {
      std::cout << options.help() << "\n";
      std::exit(EXIT_SUCCESS);
    }

    // Validate options
    auto threads = int{};
    if (opts.count("v0-threads")) {
      threads = opts["v0-threads"].as<int>();
    } else {
      threads = std::max(1u, std::thread::hardware_concurrency());
    }
    if (threads < 1) {
      std::cerr << "ERROR: threads must be positive, got " << threads << "\n";
      std::exit(EXIT_FAILURE);
    }
    auto thread_pool = std::make_unique<ctpl::thread_pool>(threads);
    
    auto seed = uint32_t{};
    if (opts.count("v0-seed")) {
      seed = opts["v0-seed"].as<uint32_t>();
    } else {
      seed = std::random_device{}();
    }
    auto prng = std::make_unique<std::mt19937>(seed);
    std::cerr << "# Seed: " << seed << "\n";
    
    if (opts.count("v0-init-heuristic") > 0 && opts.count("v0-init-random") > 0) {
      std::cerr << "ERROR: The options --v0-init-heuristic and --v0-init-random are mutually exclusive.  Pick one.\n";
      std::exit(EXIT_FAILURE);
    }
    auto init_random =
        opts.count("v0-init-heuristic") > 0 ? false :
        opts.count("v0-init-random") > 0 ? true :
        false;

    auto tree = Phylo_tree{};
    if (opts.count("v0-in-fasta") > 0 && opts.count("v0-in-maple") > 0) {
      std::cerr << "ERROR: The options --v0-in-fasta and --v0-in-maple are mutually exclusive.  Pick one.\n";
      std::exit(EXIT_FAILURE);
    }
    if (opts.count("v0-in-fasta") > 0) {
      auto in_fasta_filename = opts["v0-in-fasta"].as<std::string>();
      auto in_fasta_is = std::ifstream{in_fasta_filename};
      if (not in_fasta_is) {
        std::cerr << "ERROR: Could not read input FASTA file " << in_fasta_filename << "\n";
        std::exit(EXIT_FAILURE);
      }
      std::cerr << "Reading fasta file " << in_fasta_filename << "\n";
      auto in_fasta = read_fasta(in_fasta_is, [](int num_read_so_far) {
        std::cerr << absl::StreamFormat("- read %d so far\n", num_read_so_far);
      });
      in_fasta_is.close();
      
      std::cerr << (init_random ? "Building random initial tree..." : "Building rough initial tree...") << "\n";
      tree = build_rough_initial_tree_from_fasta(
          in_fasta, init_random, *prng,
          [](int tips_so_far, int total_tips) {
            std::cerr << absl::StreamFormat("- added %d / %d tips\n", tips_so_far, total_tips);
          },
          [](const std::string& warning_msg) {
            std::cerr << absl::StreamFormat("WARNING: %s\n", warning_msg);
          });
    }
    if (opts.count("v0-in-maple") > 0) {
      auto in_maple_filename = opts["v0-in-maple"].as<std::string>();
      auto in_maple_is = std::ifstream{in_maple_filename};
      if (not in_maple_is) {
        std::cerr << "ERROR: Could not read input MAPLE file " << in_maple_filename << "\n";
        std::exit(EXIT_FAILURE);
      }
      std::cerr << "Reading MAPLE file " << in_maple_filename << "\n";
      auto in_maple = read_maple(in_maple_is);
      in_maple_is.close();
      
      std::cerr << (init_random ? "Building random initial tree..." : "Building rough initial tree...") << "\n";
      tree = build_rough_initial_tree_from_maple(std::move(in_maple), init_random, *prng);
    }
    
    assert_phylo_tree_integrity(tree);
    
    auto tip_count = (std::ssize(tree) + 1) / 2;
    std::cerr << absl::StreamFormat("Built initial tree with %d tips\n", tip_count);

    auto steps = 100'000 * tip_count;  // Useful default
    if (opts.count("v0-steps")) {
      steps = opts["v0-steps"].as<int64_t>();
    }
    if (steps <= 0) {
      std::cerr << "ERROR: Number of steps must be positive, got " << steps << "\n";
      std::exit(EXIT_FAILURE);
    }

    auto log_filename = std::optional<std::string>{};
    if (opts.count("v0-out-log-file")) {
      log_filename = opts["v0-out-log-file"].as<std::string>();
    }

    auto log_every = steps / 10000;
    if (opts.count("v0-log-every")) {
      log_every = opts["v0-log-every"].as<int64_t>();
    }
    if (log_every <= 0) {
      std::cerr << "ERROR: Number of steps between log entries must be positive, got " << log_every << "\n";
      std::exit(EXIT_FAILURE);
    }

    auto trees_filename = std::optional<std::string>{};
    if (opts.count("v0-out-trees-file")) {
      trees_filename = opts["v0-out-trees-file"].as<std::string>();
    }

    auto tree_every = steps / 1000;
    if (opts.count("v0-tree-every")) {
      tree_every = opts["v0-tree-every"].as<int64_t>();
    }
    if (tree_every <= 0) {
      std::cerr << "ERROR: Number of steps between tree snapshots must be positive, got " << tree_every << "\n";
      std::exit(EXIT_FAILURE);
    }

    auto delphy_output_filename = std::optional<std::string>{};
    if (opts.count("v0-out-delphy-file")) {
      delphy_output_filename = opts["v0-out-delphy-file"].as<std::string>();
    }

    auto delphy_snapshot_every = steps / 1000;
    if (opts.count("v0-delphy-snapshot-every")) {
      delphy_snapshot_every = opts["v0-delphy-snapshot-every"].as<int64_t>();
    }
    if (delphy_snapshot_every <= 0) {
      std::cerr << "ERROR: Number of steps between .dphy run snapshots must be positive, got " << delphy_snapshot_every << "\n";
      std::exit(EXIT_FAILURE);
    }

    auto alpha_move_enabled = opts["v0-site-rate-heterogeneity"].as<bool>();

    auto fix_mutation_rate = opts["v0-fix-mutation-rate"].as<bool>();
    auto init_mu = opts["v0-init-mutation-rate"].as<double>();
    if (init_mu <= 0.0) {
      std::cerr << "ERROR: Initial mutation rate must be positive, got "
                << init_mu << " subst / site / year"  << "\n";
      std::exit(EXIT_FAILURE);
    }
    init_mu /= 365.0;  // convert to subst / site / day

    auto fix_final_pop_size = opts["v0-fix-final-pop-size"].as<bool>();
    auto init_final_pop_size = opts["v0-init-final-pop-size"].as<double>();
    if (init_final_pop_size <= 0.0) {
      std::cerr << "ERROR: Initial effective population size at time of last tip must be positive, got "
                << init_final_pop_size << " years"  << "\n";
      std::exit(EXIT_FAILURE);
    }
    init_final_pop_size *= 365.0;  // convert to days
    
    auto fix_pop_growth_rate = opts["v0-fix-pop-growth-rate"].as<bool>();
    auto init_pop_growth_rate = opts["v0-init-pop-growth-rate"].as<double>();
    init_pop_growth_rate /= 365.0;  // convert to e-foldings / day
    
    auto mpox_hack_enabled = opts["v0-mpox-hack"].as<bool>();

    auto target_coal_prior_cells = opts["v0-target-coal-prior-cells"].as<int>();
    
    // Create and configure initial run
    auto run = std::make_shared<Run>(*thread_pool, *prng, std::move(tree));
    run->set_mpox_hack_enabled(mpox_hack_enabled);
    run->set_alpha_move_enabled(alpha_move_enabled);
    run->set_mu_move_enabled(not fix_mutation_rate);
    run->set_mu(init_mu);
    run->set_target_coal_prior_cells(target_coal_prior_cells);
    run->set_final_pop_size_move_enabled(not fix_final_pop_size);
    run->set_final_pop_size(init_final_pop_size);
    run->set_pop_growth_rate_move_enabled(not fix_pop_growth_rate);
    run->set_pop_growth_rate(init_pop_growth_rate);

    // Output BEAST input XML if needed
    if (opts.count("v0-out-beast-xml")) {
      auto beast_xml_filename = opts["v0-out-beast-xml"].as<std::string>();
      auto beast_xml_os = std::ofstream{beast_xml_filename};
      if (not beast_xml_os) {
        std::cerr << "ERROR: Could not write BEAST XML input file \'" << beast_xml_filename << "\'" << "\n";
        std::exit(EXIT_FAILURE);
      }

      // Very rough rule of thumb: 10 Delphy steps = 1 BEAST2 step
      export_beast_input(*run, beast_xml_os, steps / 10, std::max(1L, log_every / 10), std::max(1L, tree_every / 10));
    }

    // Dry run => stop here
    if (opts.count("dry-run")) {
      std::exit(EXIT_SUCCESS);
    }
    
    return {
      .prng = std::move(prng),
      .thread_pool = std::move(thread_pool),
      .run = std::move(run),
      .steps = steps,
      .log_filename = std::move(log_filename),
      .log_every = log_every,
      .trees_filename = std::move(trees_filename),
      .tree_every = tree_every,
      .delphy_output_filename = std::move(delphy_output_filename),
      .delphy_snapshot_every = delphy_snapshot_every,
      .alpha_move_enabled = alpha_move_enabled,
      .mu_move_enabled = not fix_mutation_rate,
      .init_mu = init_mu,
      .mpox_hack_enabled = mpox_hack_enabled,
      .final_pop_size_move_enabled = not fix_final_pop_size,
      .init_final_pop_size = init_final_pop_size,
      .pop_growth_rate_move_enabled = not fix_pop_growth_rate,
      .init_pop_growth_rate = init_pop_growth_rate,
      .target_coal_prior_cells = target_coal_prior_cells
    };
    
  } catch (cxxopts::exceptions::exception& x) {
    std::cerr << "ERROR: " << x.what() << "\n" << options.help() << "\n";
    std::exit(EXIT_FAILURE);
  }
}

}  // namespace delphy
