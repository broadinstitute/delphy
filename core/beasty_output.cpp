#include "beasty_output.h"

#include <ranges>

#include "phylo_tree_calc.h"
#include "io.h"

namespace delphy {

Beasty_log_output::Beasty_log_output(std::ostream* os, bool own_stream)
    : os_{os}, own_stream_{own_stream} {}

Beasty_log_output::~Beasty_log_output() {
  if (own_stream_) {
    delete os_;
  }
}

auto Beasty_log_output::output_headers(const Run& run) -> void {
  // WARNING: Ensure this matches the logger output from export_beast_input (beasty_input.cpp)
  *os_ << stamp_version_into_log_file{};
  *os_ << "Sample\t"
       << "posterior_for_Delphy\t"
       << "likelihood_really_logG\t"
       << "prior\t"
       << "treeLikelihood_really_logG\t"
       << "TreeHeight\t";
  if (run.mu_move_enabled()) {
    *os_ << "clockRate\t";
  }
  if (run.alpha_move_enabled()) {
    *os_ << "gammaShape\t";
  }
  *os_ << "kappa\t"
       << "CoalescentExponential\t"
       << "ePopSize\t"
       << "growthRate\t"
       << "freqParameter.1\t"
       << "freqParameter.2\t"
       << "freqParameter.3\t"
       << "freqParameter.4\t"
       << "\n";
}

auto Beasty_log_output::output_log(const Run& run) -> void {

  const auto& tree = run.tree();
  
  // WARNING: Ensure this matches the logger output from export_beast_input (beasty_input.cpp)
  
  // BEAST2 seems to measure time backwards from the time of the latest tip
  auto beast_t0 = -std::numeric_limits<double>::max();
  for (const auto& node : index_order_traversal(tree)) {
    if (tree.at(node).is_tip()) {
      beast_t0 = std::max(beast_t0, tree.at(node).t);
    }
  }

  // We don't incrementally update or cache the value of the priors other than log_coalescent_prior
  auto log_other_priors = run.calc_cur_log_other_priors();
  auto log_prior = run.calc_cur_log_coalescent_prior() + log_other_priors;

  auto num_inner_nodes = Node_index{(tree.size() - 1) / 2};

  *os_ << run.step() << "\t"
       << log_prior + run.log_G() << "\t"
       << run.log_G() << "\t"  // log(G) instead of log_likelihood (which we don't calculate!)
       << log_prior << "\t"
       << run.log_G() << "\t"  // log(G) instead of log_likelihood (which we don't calculate!)
      
      // We measure time in days since 2020; ref BEAST run in years since time of latest tip
       << (beast_t0 - tree.at_root().t) / 365 << "\t";
  
  if (run.mu_move_enabled()) {
    // We measure time in days since 2020; ref BEAST run in years since time of latest tip
    *os_ << (run.mu() * 365.0) << "\t";
  }
  if (run.alpha_move_enabled()) {
    *os_ << run.alpha() << "\t";
  }
  *os_ << run.hky_kappa() << "\t"
      
      // Coalescent prior has units of (1/time)^(# coalescences)
      // We measure time in days since 2020; ref BEAST run in years since time of latest tip
       << run.log_coalescent_prior() + num_inner_nodes * std::log(365.0) << "\t"
      
      // We measure time in days since 2020; ref BEAST run in years since time of latest tip
       << run.pop_model().pop_at_time(beast_t0)/365 << "\t"
      
      // We measure time in days since 2020; ref BEAST run in years since time of latest tip
       << run.pop_model().growth_rate()*365 << "\t"
      
       << run.hky_pi()[Real_seq_letter::A] << "\t"
       << run.hky_pi()[Real_seq_letter::C] << "\t"
       << run.hky_pi()[Real_seq_letter::G] << "\t"
       << run.hky_pi()[Real_seq_letter::T]
       << "\n";
}

auto Beasty_log_output::output_footers(const Run& /*run*/) -> void {
}

auto Beasty_log_output::flush() -> void {
  os_->flush();
}

Beasty_trees_output::Beasty_trees_output(std::ostream* os, bool own_stream)
    : os_{os}, own_stream_{own_stream} {}

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
