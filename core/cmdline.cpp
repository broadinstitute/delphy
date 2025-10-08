#include "cmdline.h"

#include <iostream>
#include <fstream>
#include <optional>

#include "absl/random/random.h"
#include "cxxopts.hpp"

#include "dates.h"
#include "io.h"
#include "phylo_tree.h"
#include "phylo_tree_calc.h"
#include "sequence_utils.h"
#include "version.h"
#include "beasty_input.h"

namespace delphy {

bool delphy_invoked_via_cli{false};
std::vector<std::string> delphy_cli_args{};

// Repackage a set of equal-length sequences as a reference + diffs
auto fasta_to_maple(
    const std::vector<Fasta_entry>& in_fasta,
    const std::function<void(int,int)>& progress_hook,
    const std::function<void(const std::string&, const Sequence_warning&)>& warning_hook)
    -> Maple_file {
  
  auto result = Maple_file{};
  
  if (in_fasta.empty()) {
    throw std::runtime_error("Input has no sequences!");
  }
  auto L = std::ssize(in_fasta[0].sequence);
  for (auto i = 0; i != std::ssize(in_fasta); ++i) {
    if (std::ssize(in_fasta[i].sequence) != L) {
      throw std::runtime_error(absl::StrFormat(
          "Sequences have different lengths: '%s' has %d sites, while '%s' has %d.  Are these sequences aligned?",
          in_fasta[i].id, std::ssize(in_fasta[i].sequence),
          in_fasta[0].id, std::ssize(in_fasta[0].sequence)));
    }
  }

  auto total_seqs = std::ssize(in_fasta);
  progress_hook(0, total_seqs);
  
  // Build out "consensus" sequence from all samples, which we'll use to initially impute missing data
  // (helps prevent things like a site that's missing in almost all sequences being imputed as something
  // improbable, and then local topo/resampling moves have a very hard time untangling that error)
  auto consensus_sequence = deduce_consensus_sequence(
      in_fasta | std::views::transform([](const auto& e) { return e.sequence; }));

  // Arbitrary pick the consensus sequence as a reference
  result.ref_sequence = consensus_sequence;

  // Set up tips by converting every sequence into a list of mutations from the reference
  // (essentially, convert the FASTA to a MAPLE file)
  auto seqs_so_far = 0;
  for (auto& fasta_entry : in_fasta) {
    // Ignore sequences that we can't date correctly
    if (auto opt_t_range = extract_date_range_from_sequence_id(fasta_entry.id); not opt_t_range.has_value()) {
      warning_hook(fasta_entry.id, Sequence_warnings::No_valid_date{});
    } else {
      
      auto [t_min, t_max] = opt_t_range.value();
      auto delta = calculate_delta_from_reference(
          fasta_entry.id, fasta_entry.sequence, result.ref_sequence, warning_hook);

      result.tip_descs.push_back({
          .name = fasta_entry.id,
          .t_min = static_cast<float>(t_min),
          .t_max = static_cast<float>(t_max),
          .seq_deltas = std::move(delta.seq_deltas),
          .missations = std::move(delta.missations)});
    }

    ++seqs_so_far;
    progress_hook(seqs_so_far, total_seqs);
  }
  
  return result;
}


auto build_rough_initial_tree_from_maple(
    Maple_file&& in_maple,
    bool random,
    absl::BitGenRef bitgen,
    const std::function<void(int,int)>& progress_hook)
    -> Phylo_tree {

  if (in_maple.tip_descs.empty()) {
    throw std::runtime_error("Input has no sequences (after discarding those with warnings)!");
  }
  if (std::ssize(in_maple.tip_descs) == 1) {
    throw std::runtime_error("Input has only 1 sequence (after discarding those with warnings)!");
  }
  
  // Join all the tips up in a very rough approximation to greedy parsimony
  if (random) {
    return build_random_tree(std::move(in_maple.ref_sequence), std::move(in_maple.tip_descs), bitgen, progress_hook);
  } else {
    return build_usher_like_tree(std::move(in_maple.ref_sequence), std::move(in_maple.tip_descs), bitgen, progress_hook);
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
      ("v0-out-delphy-metadata-file", "Filename for JSON metadata blob to append to .dphy run file",
       cxxopts::value<std::string>())
      ("v0-delphy-snapshot-every", "Steps between run snapshots in .dphy run file (default: steps / 100)",
       cxxopts::value<int64_t>())
      ("v0-site-rate-heterogeneity", "Enable site rate heterogeneity",
       cxxopts::value<bool>()->default_value("false"))
      ("v0-fix-mutation-rate", "Fix mutation rate",
       cxxopts::value<bool>()->default_value("false"))
      ("v0-init-mutation-rate", "Initial (or fixed) value of mutation rate, in subst / site / year",
       cxxopts::value<double>()->default_value("1e-3"))
      ("v0-out-beast-xml", "Filename for XML input file suitable for BEAST X or BEAST2 (see --v0-out-beast-version)",
       cxxopts::value<std::string>())
      ("v0-out-beast-version", "BEAST version to target for XML and output files.  Currently allowed values are '2.6.2', '2.7.7', and 'X-10.5.0'.",
       cxxopts::value<std::string>()->default_value("X-10.5.0"))
      ("v0-mpox-hack", "Enable mpox hack (very crude treatment of APOBEC3 mechanism)",
       cxxopts::value<bool>()->default_value("false"))
      ("v0-target-coal-prior-cells", "Target number of cells to use in parallelized coalescent prior (coalescent prior resolution is adjusted if actual number is more than 33% away from target); higher is more accurate but more expensive",
       cxxopts::value<int>()->default_value("400"))

      // Popultation model
      ("v0-pop-model", "Type of population model to use - exponential or skygrid",
       cxxopts::value<std::string>()->default_value("exponential"))
      
      // The following only make sense for an exponential population model
      ("v0-fix-final-pop-size", "[pop-model == exponential] Fix effective population size at time of last tip",
       cxxopts::value<bool>()->default_value("false"))
      ("v0-init-final-pop-size", "[pop-model == exponential] Initial (or fixed) value of effective population size at time of last tip, in years (the effective population size includes a factor of the generation time)",
       cxxopts::value<double>()->default_value("3.0"))
      ("v0-fix-pop-growth-rate", "[pop-model == exponential] Fix effective population size growth rate",
       cxxopts::value<bool>()->default_value("false"))
      ("v0-init-pop-growth-rate", "[pop-model == exponential] Initial (or fixed) value of effective population size growth rate, in e-foldings / year",
       cxxopts::value<double>()->default_value("0.0"))

      // The following only make sense for a Skygrid population model
      ("v0-skygrid-type",
       "[pop-model == skygrid] Type of Skygrid model: \n"
       "* 'staircase' (default) = log N(t) is piecewise flat (standard Gill et al 2012)\n"
       "* 'log-linear' = log N(t) is piecewise linear and continuous\n",
       cxxopts::value<std::string>()->default_value("staircase"))
      ("v0-skygrid-num-parameters",
       "[pop-model == skygrid] Number of parameters used to parametrize N(t) (what BEAUTi calls 'Number of parameters' == one more than what Gill et al 2012 call `M`)",
       cxxopts::value<int>()->default_value("50"))
      ("v0-skygrid-cutoff",
       "[pop-model == skygrid] Time before time of last tip for final transition, in years (what BEAUTi calls 'Time of last transition point' == what Gill et al 2012 calls `x_M`)",
       cxxopts::value<double>()->default_value("0.0"))
      ("v0-skygrid-infer-prior-smoothness",
       "[pop-model == skygrid] Whether to fix the log-population curve's prior smoothness (Delphy default; see --v0-skygrid-prior-double-half-time) or infer it (BEAST default; see --v0-skygrid--tau-prior-alpha and --v0-skygrid-tau-prior-beta)",
       cxxopts::value<bool>()->default_value("false"))
      ("v0-skygrid-prior-double-half-time",
       "[pop-model == skygrid] Typical timescale (in years) over which the a priori population curve fluctuates by a factor of 2 (default = 30 days = 0.0822 years).  Applies only if --v0-skygrid-infer-prior-smoothness is `false`.  More precisely, the log-population curve prior is a 1D random walk with diffusion constant D, such that after after a time T, the log-population is offset by a Gaussian random variable with mean zero and variance 2 D T.  Setting sqrt(2 D T) = log(2) at the 'double-half time' T yields D = log^2(2) / (2 T); at this diffusion rate, after one month, there's a 68% chance that the population has changed by a factor of 2 or less.  In Gill et al's 2012 formulation, this fixes the smoothness precision tau = 1 /(D dt), where dt = cutoff / (num_parameters - 1)",
       cxxopts::value<double>()->default_value("0.0822"))
      ("v0-skygrid-tau",
       "[pop-model == skygrid] Fixed value of Skygrid precision parameter `tau`; mutually exclusive with --v0-skygrid-prior-double-half-time",
       cxxopts::value<double>())
      ("v0-skygrid-tau-prior-alpha",
       "[pop-model == skygrid] When inferring the log-population curve's prior smoothness `tau`, the hyperprior on tau is ~ tau^{alpha - 1} exp[-beta tau]",
       cxxopts::value<double>()->default_value("0.001"))
      ("v0-skygrid-tau-prior-beta",
       "[pop-model == skygrid] When inferring the log-population curve's prior smoothness `tau`, the hyperprior on tau is ~ tau^{alpha - 1} exp[-beta tau]",
       cxxopts::value<double>()->default_value("0.001"))
      ("v0-skygrid-disable-low-pop-barrier",
       "[pop-model == skygrid] Disable quadratic penalty against effective population sizes dropping below a certain limit (--v0-skygrid-low-pop-barrier-loc, default=1 day to match resolution of tip dates).  Such low population sizes sometimes trigger numerical difficulties that cause some MCMC moves to grind to a halt.  Moreover, they cannot reliably be estimated from the data.  We recommend disabling this only if attempting to make an exact comparison to a run from another tool, like BEAST.",
       cxxopts::value<bool>()->default_value("false"))
      ("v0-skygrid-low-pop-barrier-loc",
       "[pop-model == skygrid] Minimum value of N(t), in years, below which we penalize smaller populations (default = 1 day = 1/365 years, to match resolution of tip dates)",
       cxxopts::value<double>()->default_value(absl::StrFormat("%.5g", 1.0 / 365.0)))
      ("v0-skygrid-low-pop-barrier-scale",
       "[pop-model == skygrid] Fraction by which N(t) must drop below minimum for the value of the low-pop barrier to reach a value of 1 nat (i.e., (log(N/N_min) / log(1-scale))^2 = 1 at N = N_min * (1-scale)).",
       cxxopts::value<double>()->default_value("0.30"))
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
    auto prng = std::mt19937{seed};
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

    auto maple_file = Maple_file{};  // Filled in directly with a MAPLE file, or indirectly with a FASTA file
    if (opts.count("v0-in-fasta") > 0) {
      auto in_fasta_filename = opts["v0-in-fasta"].as<std::string>();
      auto in_fasta_is = std::ifstream{in_fasta_filename};
      if (not in_fasta_is) {
        std::cerr << "ERROR: Could not read input FASTA file " << in_fasta_filename << "\n";
        std::exit(EXIT_FAILURE);
      }
      in_fasta_is.seekg(0, std::ios_base::end);
      auto total_bytes = in_fasta_is.tellg();
      in_fasta_is.seekg(0, std::ios_base::beg);
      
      std::cerr << "Reading fasta file " << in_fasta_filename << "\n";
      auto in_fasta = read_fasta(
          in_fasta_is,
          [total_bytes](int seqs_so_far, size_t bytes_so_far) {
            std::cerr << absl::StreamFormat("- read %d sequences so far (%d of %d bytes = %.1f%%)\n",
                                            seqs_so_far, bytes_so_far, total_bytes,
                                            100.0 * bytes_so_far / static_cast<double>(total_bytes));
          },
          default_sequence_warning_hook);
      in_fasta_is.close();

      std::cerr << "Analysing fasta file " << in_fasta_filename << "\n";
      maple_file = fasta_to_maple(
          in_fasta, 
          [](int seqs_so_far, int total_seqs) {
            std::cerr << absl::StreamFormat("- analysed %d / %d sequences\n", seqs_so_far, total_seqs);
          },
          default_sequence_warning_hook);
    }
    if (opts.count("v0-in-maple") > 0) {
      auto in_maple_filename = opts["v0-in-maple"].as<std::string>();
      auto in_maple_is = std::ifstream{in_maple_filename};
      if (not in_maple_is) {
        std::cerr << "ERROR: Could not read input MAPLE file " << in_maple_filename << "\n";
        std::exit(EXIT_FAILURE);
      }
      in_maple_is.seekg(0, std::ios_base::end);
      auto total_bytes = in_maple_is.tellg();
      in_maple_is.seekg(0, std::ios_base::beg);
      
      std::cerr << "Reading MAPLE file " << in_maple_filename << "\n";
      maple_file = read_maple(
          in_maple_is,
          [total_bytes](int seqs_so_far, size_t bytes_so_far) {
            std::cerr << absl::StreamFormat("- read %d sequences so far (%d of %d bytes = %.1f%%)\n",
                                            seqs_so_far, bytes_so_far, total_bytes,
                                            100.0 * bytes_so_far / static_cast<double>(total_bytes));
          },
          default_sequence_warning_hook);
      in_maple_is.close();
    }
    
    std::cerr << (init_random ? "Building random initial tree..." : "Building rough initial tree...") << "\n";
    tree = build_rough_initial_tree_from_maple(
        std::move(maple_file), init_random, prng,
        [](int tips_so_far, int total_tips) {
          std::cerr << absl::StreamFormat("- added %d / %d tips\n", tips_so_far, total_tips);
        });
    
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

    auto delphy_output_metadata = std::optional<std::string>{};
    if (opts.count("v0-out-delphy-metadata-file")) {
      auto delphy_output_metadata_filename = opts["v0-out-delphy-metadata-file"].as<std::string>();
      auto fs = std::ifstream{delphy_output_metadata_filename};
      if (not fs) {
        std::cerr << "ERROR: Could not open .dphy metadata file " << delphy_output_metadata_filename << "\n";
        std::exit(EXIT_FAILURE);
      }

      auto sstr = std::ostringstream{};
      sstr << fs.rdbuf();
      delphy_output_metadata = std::move(sstr).str();
      std::cerr << *delphy_output_metadata << "\n";
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

    auto mpox_hack_enabled = opts["v0-mpox-hack"].as<bool>();

    auto target_coal_prior_cells = opts["v0-target-coal-prior-cells"].as<int>();

    // Create and configure initial run
    auto t0 = calc_max_tip_time(tree);
    auto run = std::make_shared<Run>(*thread_pool, prng, std::move(tree));
    run->set_mpox_hack_enabled(mpox_hack_enabled);
    run->set_alpha_move_enabled(alpha_move_enabled);
    run->set_mu_move_enabled(not fix_mutation_rate);
    run->set_mu(init_mu);
    run->set_target_coal_prior_cells(target_coal_prior_cells);

    // Population model
    auto has_exp_pop_model_parameters =
        (opts.count("v0-fix-final-pop-size") > 0) ||
        (opts.count("v0-init-final-pop-size") > 0) ||
        (opts.count("v0-fix-pop-growth-rate") > 0) ||
        (opts.count("v0-init-pop-growth-rate") > 0);
    auto has_skygrid_pop_model_parameters =
        (opts.count("v0-skygrid-type") > 0) ||
        (opts.count("v0-skygrid-num-parameters") > 0) ||
        (opts.count("v0-skygrid-cutoff") > 0) ||
        (opts.count("v0-skygrid-infer-prior-smoothness") > 0) ||
        (opts.count("v0-skygrid-prior-double-half-time") > 0) ||
        (opts.count("v0-skygrid-tau") > 0) ||
        (opts.count("v0-skygrid-tau-prior-alpha") > 0) ||
        (opts.count("v0-skygrid-tau-prior-beta") > 0) ||
        (opts.count("v0-skygrid-disable-low-pop-barrier") > 0) ||
        (opts.count("v0-skygrid-low-pop-barrier-loc") > 0) ||
        (opts.count("v0-skygrid-low-pop-barrier-scale") > 0);

    if (opts["v0-pop-model"].as<std::string>() == "exponential") {
      
      // Exponential population model
      
      if (has_skygrid_pop_model_parameters) {
        std::cerr << "ERROR: Cannot specify parameters for 'skygrid' model when pop-model is 'exponential'\n";
        std::exit(EXIT_FAILURE);
      }
      
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

      // Configure run
      run->set_pop_model(std::make_unique<Exp_pop_model>(t0, init_final_pop_size, init_pop_growth_rate));
      run->set_final_pop_size_move_enabled(not fix_final_pop_size);
      run->set_pop_growth_rate_move_enabled(not fix_pop_growth_rate);
    
    } else if (opts["v0-pop-model"].as<std::string>() == "skygrid") {
      
      // Skygrid population model
      
      if (has_exp_pop_model_parameters) {
        std::cerr << "ERROR: Cannot specify parameters for 'exponential' model when pop-model is 'skygrid'\n";
        std::exit(EXIT_FAILURE);
      }

      auto skygrid_type = Skygrid_pop_model::Type::k_staircase;
      if (opts["v0-skygrid-type"].as<std::string>() == "staircase") {
        skygrid_type = Skygrid_pop_model::Type::k_staircase;
      } else if (opts["v0-skygrid-type"].as<std::string>() == "log-linear") {
        skygrid_type = Skygrid_pop_model::Type::k_log_linear;
      } else {
        std::cerr << "ERROR: Unknown Skygrid model type '" << opts["v0-skygrid-type"].as<std::string>() << "'\n";
        std::exit(EXIT_FAILURE);
      }

      auto num_parameters = opts["v0-skygrid-num-parameters"].as<int>();
      if (num_parameters < 2) {
        std::cerr << "ERROR: Number of Skygrid parameters must be at least 2\n";
        std::exit(EXIT_FAILURE);
      }

      auto skygrid_cutoff = opts["v0-skygrid-cutoff"].as<double>();
      if (skygrid_cutoff < 0.0) {
        std::cerr << "ERROR: Skygrid cutoff must be positive\n";
        std::exit(EXIT_FAILURE);
      }
      skygrid_cutoff *= 365.0;  // convert to days

      auto M = num_parameters - 1;
      auto x_k = std::vector<double>(M+1, 0.0);
      auto dt = (skygrid_cutoff / M);
      for (auto k = 0; k <= M; ++k) {
        x_k[k] = (t0 - skygrid_cutoff) + k * dt;
      }

      auto skygrid_tau = 0.0;
      auto infer_tau = opts["v0-skygrid-infer-prior-smoothness"].as<bool>();
      if (infer_tau) {
        // Infer a priori smoothness of log-population curve
        auto prior_alpha = opts["v0-skygrid-tau-prior-alpha"].as<double>();
        auto prior_beta = opts["v0-skygrid-tau-prior-beta"].as<double>();
        if (prior_alpha <= 0.0 || prior_beta <= 0.0) {
          std::cerr << "ERROR: Skygrid tau prior parameters must be positive\n";
          std::exit(EXIT_FAILURE);
        }
        
        run->set_skygrid_tau_prior_alpha(prior_alpha);
        run->set_skygrid_tau_prior_beta(prior_beta);

        skygrid_tau = prior_alpha / prior_beta;
        
      } else {
        // Fix a priori smoothness of log-population curve
        if (opts.count("v0-skygrid-tau") > 0 &&
            opts.count("v0-skygrid-prior-double-half-time") > 0) {
          std::cerr << "ERROR: Skygrid tau can be fixed either directly (--v0-skygrid-tau) or by specifying an effective 'double-half' time (--v0-skygrid-prior-double-half-time), but not both\n";
          std::exit(EXIT_FAILURE);
        }

        if (opts.count("v0-skygrid-tau") > 0) {
          skygrid_tau = opts["v0-skygrid-tau"].as<double>();
          if (skygrid_tau <= 0.0) {
            std::cerr << "ERROR: Skygrid tau parameter must be positive\n";
            std::exit(EXIT_FAILURE);
          }
          
        } else {
          auto double_half_time = opts["v0-skygrid-prior-double-half-time"].as<double>();
          double_half_time *= 365.0;  // convert to days
          if (double_half_time <= 0.0) {
            std::cerr << "ERROR: Skygrid prior 'double-half' time must be positive\n";
            std::exit(EXIT_FAILURE);
          }
          
          // Setting tau = 1 / (2 D dt), the prior for the log-population curve looks like
          // a 1D random walk with diffusion constant D.  Hence, on average, a starting
          // log-population changes after a time T by a root-mean-square deviation of
          // `sqrt(2 D T)`.  We parametrize D such that after T = double_half_time, the
          // rms deviation is log(2), i.e., population changes by up to a factor of ~2 in
          // the "double-half time" with 68% probability:
          //
          //   sqrt(2 D T) = log(2)  => D = log^2(2) / (2 T).
          //
          const auto D = std::pow(std::log(2.0), 2) / (2 * double_half_time);
          skygrid_tau = 1.0 / (2 * D * dt);
        }
      }

      // At this point, skygrid_tau is set.  Sample a random trajectory with this
      // precision, then reset its mean to "3 years" (the `skygrid_gammas_zero_mode_gibbs_move`
      // Gibbs sample the mean value of gamma, so we don't need to get this right at all here)
      auto gamma_k = std::vector<double>(M+1, 0.0);
      auto mean_gamma_k = 0.0;
      for (auto k = 1; k <= M; ++k) {
        gamma_k[k] = gamma_k[k-1] + absl::Gaussian(prng, 0.0, std::sqrt(1.0/skygrid_tau));
        mean_gamma_k += gamma_k[k];
      }
      mean_gamma_k /= M+1;
      for (auto k = 0; k <= M; ++k) {
        gamma_k[k] += (-mean_gamma_k) + std::log(3.0 * 365.0);
      }

      // Configure run
      run->set_pop_model(std::make_unique<Skygrid_pop_model>(std::move(x_k), std::move(gamma_k), skygrid_type));
      run->set_skygrid_tau(skygrid_tau);
      run->set_skygrid_tau_move_enabled(infer_tau);  // Prior configured above when infer_tau == true

      // Low-pop barrier
      if (opts["v0-skygrid-disable-low-pop-barrier"].as<bool>()) {
        run->set_skygrid_low_gamma_barrier_enabled(false);
      } else {
        auto low_pop_barrier_loc = opts["v0-skygrid-low-pop-barrier-loc"].as<double>();
        low_pop_barrier_loc *= 365.0;  // convert to days
        auto low_gamma_barrier_loc = std::log(low_pop_barrier_loc);

        auto low_pop_barrier_scale = opts["v0-skygrid-low-pop-barrier-scale"].as<double>();
        auto low_gamma_barrier_scale = -std::log(1 - low_pop_barrier_scale);  // Convert to scale in gamma
        
        run->set_skygrid_low_gamma_barrier_enabled(true);
        run->set_skygrid_low_gamma_barrier_loc(low_gamma_barrier_loc);
        run->set_skygrid_low_gamma_barrier_scale(low_gamma_barrier_scale);
      }
      
    } else {
      std::cerr << "ERROR: Unknown population model '" << opts["v0-pop-model"].as<std::string>() << "'\n";
      std::exit(EXIT_FAILURE);
    }

    // Which version of BEAST to target for input XML and output log & trees files
    auto beast_version = opts["v0-out-beast-version"].as<std::string>();
    
    // Output BEAST input XML if needed
    if (opts.count("v0-out-beast-xml")) {
      auto beast_xml_filename = opts["v0-out-beast-xml"].as<std::string>();
      auto beast_xml_os = std::ofstream{beast_xml_filename};
      if (not beast_xml_os) {
        std::cerr << "ERROR: Could not write BEAST XML input file \'" << beast_xml_filename << "\'" << "\n";
        std::exit(EXIT_FAILURE);
      }

      // Very rough rule of thumb: 10 Delphy steps = 1 BEAST2 step
      export_beast_input(*run, beast_version, beast_xml_os, steps / 10, std::max(1L, log_every / 10), std::max(1L, tree_every / 10));
    }

    // Dry run => stop here
    if (opts.count("dry-run")) {
      std::exit(EXIT_SUCCESS);
    }
    
    return {
      .thread_pool = std::move(thread_pool),
      .run = std::move(run),
      .steps = steps,
      .log_filename = std::move(log_filename),
      .log_every = log_every,
      .trees_filename = std::move(trees_filename),
      .tree_every = tree_every,
      .delphy_output_filename = std::move(delphy_output_filename),
      .delphy_output_metadata = std::move(delphy_output_metadata),
      .delphy_snapshot_every = delphy_snapshot_every,
      .beast_version = std::move(beast_version)
    };
    
  } catch (cxxopts::exceptions::exception& x) {
    std::cerr << "ERROR: " << x.what() << "\n" << options.help() << "\n";
    std::exit(EXIT_FAILURE);
  }
}

}  // namespace delphy
