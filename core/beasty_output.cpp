#include "beasty_output.h"

#include <ranges>

#include "phylo_tree_calc.h"
#include "io.h"
#include "dates.h"

namespace delphy {

// Opaque implementer for each BEAST version
class Beasty_log_output_version_impl {
 public:
  virtual ~Beasty_log_output_version_impl() {}

  virtual auto output_headers(std::ostream& /*os*/, const Run& /*run*/) -> void {};
  virtual auto output_log(std::ostream& /*os*/, const Run& /*run*/) -> void {};
  virtual auto output_footers(std::ostream& /*os*/, const Run& /*run*/) -> void {};
};

static auto output_skygrid_headers(std::ostream& os, const Skygrid_pop_model& pop_model) -> void {
  // Even though we don't support BEAST2 XML input files with Skygrid, having the Skygrid results in the
  // log files in the same format as BEAST X produces them is still useful
  os << "skygrid.isloglinear\t"   // Not output by BEAST, but useful for downstream analysis of Delphy results
     << "skygrid.precision\t";     // tau
  for (auto k = 0; k != (pop_model.M()+1); ++k) {
    os << "skygrid.logPopSize" << (k+1) << "\t";
  }
  os << "skygrid.cutOff\t";        // K
}

static auto output_skygrid_results(std::ostream& os, const Run& run, const Skygrid_pop_model& pop_model) -> void {
  switch (pop_model.type()) {
    case Skygrid_pop_model::Type::k_staircase:  os << "0\t"; break;
    case Skygrid_pop_model::Type::k_log_linear: os << "1\t"; break;
    default:                                    os << "-1\t"; break;
  }
  os << run.skygrid_tau() << "\t";
  for (auto k = 0; k != (pop_model.M()+1); ++k) {
    // gamma = ln(N(t)), where N(t) is measured in days in Delphy, in years in BEAST (hence -ln(365.0) and (M-k))
    os << pop_model.gamma(pop_model.M() - k) - std::log(365.0) << "\t";
  }
  os << (pop_model.K() / 365.0) << "\t";         // years
}

class Beasty_log_output_2_6_2 : public Beasty_log_output_version_impl {
 public:
  virtual auto output_headers(std::ostream& os, const Run& run) -> void {
    if (not run.mpox_hack_enabled()) {
      // WARNING: Ensure this matches the logger output from export_beast_2_6_2_input (beasty_input.cpp)
      os << stamp_version_into_log_file{};
      os << "Sample\t"
         << "numMuts\t"    // Not output by BEAST, but useful for downstream analysis of Delphy results
         << "posterior_for_Delphy\t"        // Not really BEAST's posterior
         << "likelihood_really_logG\t"      // Not really BEAST's likelihood
         << "prior_for_Delphy\t"            // Not really BEAST's prior
         << "treeLikelihood_really_logG\t"  // Not really BEAST's likelihood
         << "TreeHeight\t";
      if (run.mu_move_enabled()) {
        os << "clockRate\t";
      }
      if (run.alpha_move_enabled()) {
        os << "gammaShape\t";
      }
      os << "kappa\t";
      os << "Coalescent\t";
      
      const auto& pop_model = run.pop_model();
      if (typeid(pop_model) == typeid(Exp_pop_model)) {
        if (run.final_pop_size_move_enabled()) {
          os << "ePopSize\t";
        }
        if (run.pop_growth_rate_move_enabled()) {
          os << "growthRate\t";
        }
        
      } else if (typeid(pop_model) == typeid(Skygrid_pop_model)) {
        const auto& skygrid_pop_model = static_cast<const Skygrid_pop_model&>(run.pop_model());
        output_skygrid_headers(os, skygrid_pop_model);
        
      } else {
        std::cerr << "ERROR (output_headers): Unrecognized population model type " << typeid(pop_model).name() << '\n';
      }
      
      os << "freqParameter.1\t"
         << "freqParameter.2\t"
         << "freqParameter.3\t"
         << "freqParameter.4\t";
      os << "\n";
    } else {
      os << "state\t"
         << "numMuts\t"
         << "posterior_for_Delphy\t"
         << "prior\t"
         << "likelihood_really_logG\t"
         << "treeModel.rootHeight\t"
         << "age(root)\t"
         << "treeLength\t";
      
      const auto& pop_model = run.pop_model();
      if (typeid(pop_model) == typeid(Exp_pop_model)) {
        os << "exponential.popSize\t"
           << "exponential.growthRate\t";
        
      } else if (typeid(pop_model) == typeid(Skygrid_pop_model)) {
        const auto& skygrid_pop_model = static_cast<const Skygrid_pop_model&>(run.pop_model());
        output_skygrid_headers(os, skygrid_pop_model);
        
      } else {
        std::cerr << "ERROR (output_headers): Unrecognized population model type " << typeid(pop_model).name() << '\n';
      }
      
      os << "apobec3.clock.rate\t"
         << "non_apobec3.clock.rate\t"
         << "coalescent.all\t"
         << "\n";
    }
  }
    
  virtual auto output_log(std::ostream& os, const Run& run) -> void {
    const auto& tree = run.tree();
  
    // WARNING: Ensure this matches the logger output from export_beast_2_6_2_input (beasty_input.cpp)
  
    // BEAST2 seems to measure time backwards from the time of the latest tip
    auto beast_t0 = -std::numeric_limits<double>::max();
    for (const auto& node : index_order_traversal(tree)) {
      if (tree.at(node).is_tip()) {
        beast_t0 = std::max(beast_t0, tree.at(node).t);
      }
    }
    
    auto num_inner_nodes = (std::ssize(tree) - 1) / 2;
  
    // We don't incrementally update or cache the value of the priors other than log_coalescent_prior
    auto log_other_priors = run.calc_cur_log_other_priors();
    auto log_prior = run.calc_cur_log_coalescent_prior() + log_other_priors;
    
    if (not run.mpox_hack_enabled()) {
      // WARNING: Ensure this matches the logger output from export_beast_2_6_2_input (beasty_input.cpp)
    
      os << run.step() << "\t"
         << run.num_muts() << "\t"
         << log_prior + run.log_G() << "\t"
         << run.log_G() << "\t"  // log(G) instead of log_likelihood (which we don't calculate!)
         << log_prior << "\t"
         << run.log_G() << "\t"  // log(G) instead of log_likelihood (which we don't calculate!)
        
          // We measure time in days since 2020; ref BEAST run in years since time of latest tip
         << (beast_t0 - tree.at_root().t) / 365 << "\t";
    
      if (run.mu_move_enabled()) {
        os << (run.mu() * 365.0) << "\t";
      }
      if (run.alpha_move_enabled()) {
        os << run.alpha() << "\t";
      }
      os << run.hky_kappa() << "\t";
    
      // Coalescent prior has units of (1/time)^(# coalescences)
      // We measure time in days since 2020; ref BEAST run in years since time of latest tip
      os << run.log_coalescent_prior() + num_inner_nodes * std::log(365.0) << "\t";

      const auto& pop_model = run.pop_model();
      if (typeid(pop_model) == typeid(Exp_pop_model)) {
        const auto& exp_pop_model = static_cast<const Exp_pop_model&>(run.pop_model());
        auto pop_growth_rate = exp_pop_model.growth_rate();
      
        if (run.final_pop_size_move_enabled()) {
          // We measure time in days since 2020; ref BEAST run in years since time of latest tip
          os << run.pop_model().pop_at_time(beast_t0)/365 << "\t";
        }
      
        if (run.pop_growth_rate_move_enabled()) {
          // We measure time in days since 2020; ref BEAST run in years since time of latest tip
          os << pop_growth_rate*365 << "\t";
        }
      
      } else if (typeid(pop_model) == typeid(Skygrid_pop_model)) {
        const auto& skygrid_pop_model = static_cast<const Skygrid_pop_model&>(run.pop_model());
        output_skygrid_results(os, run, skygrid_pop_model);
      }
    
      os << run.hky_pi()[Real_seq_letter::A] << "\t"
         << run.hky_pi()[Real_seq_letter::C] << "\t"
         << run.hky_pi()[Real_seq_letter::G] << "\t"
         << run.hky_pi()[Real_seq_letter::T] << "\t";
    
      os << "\n";
  
    } else {  // mpox_hack_enabled == true

      auto to_beast_date = [](double delphy_t) { return 2020.0 + delphy_t / 365.0; };
    
      os << run.step() << "\t"
         << run.num_muts() << "\t"
         << run.log_posterior() << "\t"
         << log_prior << "\t"
         << run.log_G() << "\t"  // log(G) instead of log_likelihood (which we don't calculate!)
         << (beast_t0 - tree.at_root().t) / 365 << "\t"
         << to_beast_date(tree.at_root().t) << "\t"
         << calc_T(tree) / 365.0 << "\t";
    
      const auto& pop_model = run.pop_model();
      if (typeid(pop_model) == typeid(Exp_pop_model)) {
        const auto& exp_pop_model = static_cast<const Exp_pop_model&>(run.pop_model());
        auto pop_growth_rate = exp_pop_model.growth_rate();
      
        os 
            // We measure time in days since 2020; ref BEAST run in years since time of latest tip
            << run.pop_model().pop_at_time(beast_t0)/365 << "\t"
            // We measure time in days since 2020; ref BEAST run in years since time of latest tip
            << pop_growth_rate*365 << "\t";
      
      } else if (typeid(pop_model) == typeid(Skygrid_pop_model)) {
        const auto& skygrid_pop_model = static_cast<const Skygrid_pop_model&>(run.pop_model());
        output_skygrid_results(os, run, skygrid_pop_model);
      }
    
      os << run.mpox_mu_star() * 365 << "\t"
         << run.mpox_mu() * 365 << "\t"
          // Coalescent prior has units of (1/time)^(# coalescences)
          // We measure time in days since 2020; ref BEAST run in years since time of latest tip
         << run.log_coalescent_prior() + num_inner_nodes * std::log(365.0) << "\t";
      os << "\n";
    }
  }
  
  virtual auto output_footers(std::ostream& /*os*/, const Run& /*run*/) -> void {
  }
};


class Beasty_log_output_X_10_5_0 : public Beasty_log_output_version_impl {
 public:
  virtual auto output_headers(std::ostream& os, const Run& run) -> void {
    const auto& tree = run.tree();
    const auto& pop_model = run.pop_model();
    
    // WARNING: Ensure this matches the logger output from export_beast_X_10_5_0_input (beasty_input.cpp)
    os << stamp_version_into_log_file{};
    os << "state\t"
       << "numMuts\t"    // Not output by BEAST, but useful for downstream analysis of Delphy results
       << "posterior_for_Delphy\t"        // Not really BEAST's "joint"
       << "prior_for_Delphy\t"            // Not really BEAST's prior
       << "likelihood_really_logG\t"      // Not really BEAST's likelihood
       << "rootHeight\t"
       << "age(root)\t"
       << "treeLength\t";
      
    if (typeid(pop_model) == typeid(Exp_pop_model)) {
      if (run.final_pop_size_move_enabled()) {
        os << "exponential.popSize\t";
      }
      if (run.pop_growth_rate_move_enabled()) {
        os << "exponential.growthRate\t";
      }
        
    } else if (typeid(pop_model) == typeid(Skygrid_pop_model)) {
      const auto& skygrid_pop_model = static_cast<const Skygrid_pop_model&>(run.pop_model());
      output_skygrid_headers(os, skygrid_pop_model);
        
    } else {
      std::cerr << "ERROR (output_headers): Unrecognized population model type " << typeid(pop_model).name() << '\n';
    }

    os << "kappa\t"
       << "frequencies1\t"
       << "frequencies2\t"
       << "frequencies3\t"
       << "frequencies4\t";
      
    if (run.alpha_move_enabled()) {
      os << "alpha\t";
    }
    if (not run.mpox_hack_enabled()) {
      if (run.mu_move_enabled()) {
        os << "clock.rate\t";
      }
      os << "meanRate\t";
    } else {
      // We don't produce BEAST-equivalent XML inputs when the mpox hack is on,
      // but putting these results in a log file facilitates downstream analysis of Delphy runs
      os << "apobec3.clock.rate\t"
         << "non_apobec3.clock.rate\t";
    }
    
    for (const auto& node : index_order_traversal(tree)) {
      if (tree.at(node).is_tip() && tree.at(node).t_min != tree.at(node).t_max) {
        os << absl::StreamFormat("age(%s)\t", tree.at(node).name);
      }
    }

    os << "treeLikelihood_really_logG\t"   // Not really BEAST's likelihood
       << "branchRates\t";

    if (typeid(pop_model) == typeid(Exp_pop_model)) {
      os << "coalescent\t";
    } else if (typeid(pop_model) == typeid(Skygrid_pop_model)) {
      os << "skygrid\t";
    } else {
      std::cerr << "ERROR (output_headers): Unrecognized population model type " << typeid(pop_model).name() << '\n';
    }

    os << "\n";
  }
    
  virtual auto output_log(std::ostream& os, const Run& run) -> void {
    const auto& tree = run.tree();
    const auto& pop_model = run.pop_model();
  
    // WARNING: Ensure this matches the logger output from export_beast_input (beasty_input.cpp)
  
    // BEAST X measures time backwards from the time of the latest tip
    auto beast_t0 = -std::numeric_limits<double>::max();
    for (const auto& node : index_order_traversal(tree)) {
      if (tree.at(node).is_tip()) {
        beast_t0 = std::max(beast_t0, tree.at(node).t);
      }
    }
    
    auto num_inner_nodes = (std::ssize(tree) - 1) / 2;
  
    // We don't incrementally update or cache the value of the priors other than log_coalescent_prior
    auto log_other_priors = run.calc_cur_log_other_priors();
    auto log_prior = run.calc_cur_log_coalescent_prior() + log_other_priors;
    
    // WARNING: Ensure this matches the logger output from export_beast_X_10_5_0_input (beasty_input.cpp)
    
    os << run.step() << "\t"
       << run.num_muts() << "\t"
       << log_prior + run.log_G() << "\t"
       << log_prior << "\t"
       << run.log_G() << "\t"  // log(G) instead of log_likelihood (which we don't calculate!)
        
        // We measure time in days since 2020; ref BEAST run in years since time of latest tip
       << (beast_t0 - tree.at_root().t) / 365 << "\t"
       << to_linear_year(tree.at_root().t) << "\t"
       << calc_T(tree) / 365.0 << "\t";  // years
    
    if (typeid(pop_model) == typeid(Exp_pop_model)) {
      const auto& exp_pop_model = static_cast<const Exp_pop_model&>(run.pop_model());
      auto pop_growth_rate = exp_pop_model.growth_rate();
      
      if (run.final_pop_size_move_enabled()) {
        // We measure time in days since 2020; ref BEAST run in years since time of latest tip
        os << run.pop_model().pop_at_time(beast_t0)/365 << "\t";
      }
      
      if (run.pop_growth_rate_move_enabled()) {
        // We measure time in days since 2020; ref BEAST run in years since time of latest tip
        os << pop_growth_rate*365 << "\t";
      }
      
    } else if (typeid(pop_model) == typeid(Skygrid_pop_model)) {
      const auto& skygrid_pop_model = static_cast<const Skygrid_pop_model&>(run.pop_model());
      output_skygrid_results(os, run, skygrid_pop_model);
    }
    
    os << run.hky_kappa() << "\t";
    
    os << run.hky_pi()[Real_seq_letter::A] << "\t"
       << run.hky_pi()[Real_seq_letter::C] << "\t"
       << run.hky_pi()[Real_seq_letter::G] << "\t"
       << run.hky_pi()[Real_seq_letter::T] << "\t";
    
    if (run.alpha_move_enabled()) {
      os << run.alpha() << "\t";
    }
    if (not run.mpox_hack_enabled()) {
      if (run.mu_move_enabled()) {
        os << (run.mu() * 365.0) << "\t";
      }
      os << (run.mu() * 365.0) << "\t";  // "meanRate", but we only implement a strict clock for now
    } else {
      // We don't produce BEAST-equivalent XML inputs when the mpox hack is on,
      // but putting these results in a log file facilitates downstream analysis of Delphy runs
      os << run.mpox_mu_star() * 365 << "\t"
         << run.mpox_mu() * 365 << "\t";
    }
    
    for (const auto& node : index_order_traversal(tree)) {
      if (tree.at(node).is_tip() && tree.at(node).t_min != tree.at(node).t_max) {
        os << absl::StreamFormat("%g\t", (beast_t0 - tree.at(node).t) / 365.0);
      }
    }

    os << run.log_G() << "\t"  // log(G) instead of log_likelihood (which we don't calculate!)
       << "0.0\t";             // "branchRates": no idea what this is about

    if (typeid(pop_model) == typeid(Exp_pop_model)) {
      // Coalescent prior has units of (1/time)^(# coalescences)
      // We measure time in days since 2020; ref BEAST run in years since time of latest tip
      os << run.log_coalescent_prior() + num_inner_nodes * std::log(365.0) << "\t";
      
    } else if (typeid(pop_model) == typeid(Skygrid_pop_model)) {
      
      // Coalescent prior has units of (1/time)^(# coalescences)
      // We measure time in days since 2020; ref BEAST run in years since time of latest tip
      //
      auto gmrf_prior = run.calc_cur_skygrid_gmrf_prior();
      
      os << (run.log_coalescent_prior() + num_inner_nodes * std::log(365.0)
             + gmrf_prior) << "\t";
      
    } else {
      std::cerr << "ERROR (output_log): Unrecognized population model type " << typeid(pop_model).name() << '\n';
    }
      
    os << "\n";
  }
  
  virtual auto output_footers(std::ostream& /*os*/, const Run& /*run*/) -> void {
  }
};

Beasty_log_output::Beasty_log_output(
    std::ostream* os,
    std::string beast_version,
    bool own_stream)
    : os_{os}, beast_version_{std::move(beast_version)}, own_stream_{own_stream} {
  if (beast_version_ == "2.6.2") {
    impl_ = std::make_unique<Beasty_log_output_2_6_2>();
  } else if (beast_version_ == "X-10.5.0") {
    impl_ = std::make_unique<Beasty_log_output_X_10_5_0>();
  } else {
    std::cerr << "ERROR: BEAST output generation not currently supported for BEAST version '" << beast_version_ << "'\n";
    impl_ = std::make_unique<Beasty_log_output_version_impl>();  // No-op impl
  }
}

Beasty_log_output::~Beasty_log_output() {
  if (own_stream_) {
    delete os_;
  }
}

auto Beasty_log_output::output_headers(const Run& run) -> void {
  impl_->output_headers(*os_, run);
}

auto Beasty_log_output::output_log(const Run& run) -> void {
  impl_->output_log(*os_, run);
}

auto Beasty_log_output::output_footers(const Run& run) -> void {
  impl_->output_footers(*os_, run);
}

auto Beasty_log_output::flush() -> void {
  os_->flush();
}

Beasty_trees_output::Beasty_trees_output(
    std::ostream* os,
    std::string beast_version,
    bool own_stream)
    : os_{os}, beast_version_{std::move(beast_version)}, own_stream_{own_stream} {}

Beasty_trees_output::~Beasty_trees_output() {
  if (own_stream_) {
    delete os_;
  }
}

auto Beasty_trees_output::output_headers(const Run& run) -> void {
  const auto& tree = run.tree();
  auto num_tips = (std::ssize(tree) + 1) / 2;
  auto next_tip = Node_index{1};
  node_to_tip_.assign(std::ssize(tree), k_no_node);
  for (const auto& node : index_order_traversal(tree)) {
    if (tree.at(node).is_tip()) {
      node_to_tip_[node] = next_tip;
      ++next_tip;
    }
  }

  *os_ << "#NEXUS\n"
       << stamp_version_into_log_file{}
       << "\n"
       << "Begin taxa;\n"
       << "    Dimensions ntax=" << num_tips << ";\n"
       << "        Taxlabels\n";
  
  for (const auto& node : index_order_traversal(tree)) {
    if (tree.at(node).is_tip()) {
      *os_ << "            " << tree.at(node).name << "\n";
    }
  }

  *os_ << "            ;\n"
       << "End;\n"
       << "Begin trees;\n"
       << "    Translate\n";

  auto first_tip = true;
  for (const auto& node : index_order_traversal(tree)) {
    if (tree.at(node).is_tip()) {
      if (first_tip) {
        first_tip = false;
      } else {
        *os_ << ",\n";
      }
      *os_ << absl::StrFormat("        %4d %s", node_to_tip_[node], tree.at(node).name);
    }
  }
  *os_ << "\n";

  *os_ << ";\n";
}

auto Beasty_trees_output::output_tree(const Run& run) -> void {
  *os_ << "tree STATE_" << run.step() << " = ";
  output_newick_tree(run.tree());
  *os_ << "\n";
}

auto Beasty_trees_output::output_footers(const Run& /*run*/) -> void {
  *os_ << "End;";
}

auto Beasty_trees_output::flush() -> void {
  os_->flush();
}

auto Beasty_trees_output::output_newick_tree(const Phylo_tree& tree) -> void {
  //
  // BEAST X and BEAST2 differ slightly in tree output:
  //
  // - BEAST X marks the root node with [&R] while BEAST2 does not
  // - BEAST X annotates each state (before the equals sign) with lnP and joint values
  //   (e.g., `tree STATE_0 [&lnP=-10000,joint=-5000] = [&R] (...);`), BEAST2 does not
  // - BEAST X annotates each branch with a local rate (presumably to support
  //   branch-specific rates), while BEAST2 does not
  //
  // None of these differences is sufficient to warrant separate output routines
  // for each BEAST version.  We just output what BEAST2 does.  If it ever proves
  // necessary to differentiate the trees outputs, we'd implement *_impl delegates
  // per BEAST version as we do for the logs (the version is already passed in
  // and recorded in the Beasty_trees_output constructor).
  //
  for (const auto& [node, children_so_far] : traversal(tree)) {
    if (children_so_far == 0) {
      // Entering `node`
      if (tree.at(node).is_inner_node()) {
        *os_ << "(";
      }
    }

    if (children_so_far > 0 && children_so_far < std::ssize(tree.at(node).children)) {
      // In between two children of an inner node
      *os_ << ",";
    }

    if (children_so_far == std::ssize(tree.at(node).children)) {
      // Exiting `node`
      if (tree.at(node).is_inner_node()) {
        *os_ << ")";
      }
      if (tree.at(node).is_tip()) {
        *os_ << node_to_tip_[node];
      }
      *os_ << ":";

      auto t_parent = (node == tree.root ? tree.at_root().t : tree.at_parent_of(node).t);
      
      // Add mutations, if any, as branch attribute
      auto mutations = std::vector<std::string>{};
      for (const auto& m : tree.at(node).mutations) {
        mutations.push_back(absl::StrFormat(
            "%c%d%c,%g", to_char(m.from), m.site+1, to_char(m.to), m.t - t_parent));
      }
      if (not mutations.empty()) {
        *os_ << "[&mutations={" << absl::StrJoin(mutations, ",") << "}]";
      }
      
      // Days (internal) -> Years (output)
      *os_ << (tree.at(node).t - t_parent) / 365.0;
    }
    
  }
  *os_ << ";";
}

}  // namespace delphy
