#include "beasty_input.h"

#include <absl/log/check.h>

#include "dates.h"
#include "newick.h"
#include "phylo_tree_calc.h"
#include "io.h"

namespace delphy {

auto newick_to_phylo_tree(const Newick_tree& ntree, const std::vector<std::string>& tip_names) -> Phylo_tree {
  // Newick tree had better be bifurcating!
  // And we need tip indices to line up across trees...
  auto nnode_to_pnode = Node_vector<Node_index>(std::ssize(ntree), k_no_node);
  auto num_tips = Node_index{(ntree.size() + 1) / 2};
  auto next_inner_node = num_tips;
  for (const auto& nnode : index_order_traversal(ntree)) {
    if (ntree.at(nnode).is_tip()) {
      auto tip_index = int{};
      CHECK(absl::SimpleAtoi(ntree.at(nnode).name, &tip_index));
      nnode_to_pnode.at(nnode) = tip_index - 1;  // tip_index is 1-based!
    } else {
      nnode_to_pnode.at(nnode) = next_inner_node;
      ++next_inner_node;
    }
  }
  
  auto ptree = Phylo_tree{ntree.size()};
  for (const auto& nnode : pre_order_traversal(ntree)) {
    auto pnode = nnode_to_pnode.at(nnode);
    if (nnode != ntree.root) {
      ptree.at(pnode).parent = nnode_to_pnode.at(ntree.at(nnode).parent);
    }
    if (ntree.at(nnode).children.empty()) {
      ptree.at(pnode).children = {};
    } else {
      const auto& children = ntree.at(nnode).children;
      CHECK(children.size() == 2);
      ptree.at(pnode).children = {nnode_to_pnode.at(children[0]), nnode_to_pnode.at(children[1])};
    }
    if (ntree.at(nnode).is_tip()) {
      auto tip_index = int{};
      CHECK(absl::SimpleAtoi(ntree.at(nnode).name, &tip_index));
      ptree.at(pnode).name = tip_names.at(tip_index);
    }
    if (nnode == ntree.root) {
      ptree.at(pnode).t = 0.0;  // root time at 0.0
      ptree.root = pnode;
    } else {
      ptree.at(pnode).t = ptree.at_parent_of(pnode).t + ntree.at(nnode).branch_length * 365.0;  // years to days
    }
  }
  assert_phylo_tree_integrity(ptree);
  return ptree;
}

auto read_beasty_trees(
    std::istream& is,
    std::int64_t min_state,
    std::int64_t every)
    -> std::map<std::int64_t, Phylo_tree> {
  
  // This parser is *extremely* fragile and not very efficient, but good enough to parse
  // the trees produced by Beasty_trees_output and at least one actual BEAST run

  // Adapted slightly from https://stackoverflow.com/a/6500499
  auto trim = [](const std::string& str) {
    auto left_pos = str.find_first_not_of(" \t");
    auto right_pos = str.find_last_not_of(" \t");
    if (left_pos == std::string::npos || right_pos == std::string::npos) {
      return std::string_view{""};
    } else {
      return std::string_view{str}.substr(left_pos, right_pos-left_pos+1);
    }
  };
  
  auto line = std::string{};
  CHECK(std::getline(is, line) && trim(line) == "#NEXUS");
  while (std::getline(is, line) && line[0] == '#') { continue; }
  CHECK(trim(line) == "");
  CHECK(std::getline(is, line) && trim(line) == "Begin taxa;");
  CHECK(std::getline(is, line) && trim(line).starts_with("Dimensions ntax="));
  auto ntax_sv = std::string_view{line.begin() + line.find('=')+1, line.begin() + line.find(';')};
  auto ntax = int{};
  CHECK(absl::SimpleAtoi(ntax_sv, &ntax));
  CHECK(std::getline(is, line) && trim(line) == "Taxlabels");
  
  auto tip_names = std::vector<std::string>{};
  tip_names.push_back("Tip indices are 1-based, not 0-based!");
  for (auto tip_index = 0; tip_index != ntax; ++tip_index) {
    CHECK(std::getline(is, line));
    tip_names.push_back(std::string{trim(line)});
  }
  
  CHECK(std::getline(is, line) && trim(line) == ";");
  CHECK(std::getline(is, line) && trim(line) == "End;");
  while (std::getline(is, line) && trim(line) == "") {
    ;
  }
  CHECK(trim(line) == "Begin trees;");
  CHECK(std::getline(is, line) && trim(line) == "Translate");
  // Skip names, we already have them
  for (auto tip_index = 0; tip_index != ntax; ++tip_index) {
    CHECK(std::getline(is, line));
  }
  CHECK(std::getline(is, line) && trim(line) == ";");

  auto result = std::map<std::int64_t, Phylo_tree>{};
  while (std::getline(is, line)) {
    if (line == "End;") { break; }
    CHECK(line.starts_with("tree STATE_")) << line.substr(0, 50);
    auto state_num_first_pos = std::strlen("tree STATE_");
    auto state_num_last_pos = line.find_first_not_of("0123456789", state_num_first_pos);
    CHECK_NE(state_num_last_pos, std::string::npos);
    auto state_sv = std::string_view{line}.substr(state_num_first_pos, state_num_last_pos-state_num_first_pos);
    auto state = std::int64_t{};
    CHECK(absl::SimpleAtoi(state_sv, &state)) << state_sv;

    if (state < min_state) {
      std::cerr << absl::StreamFormat("Ignoring state %d < %d (burn-in)", state, min_state) << std::endl;
    } else if ((state % every) != 0) {
      std::cerr << absl::StreamFormat("Ignoring state %d (not a multiple of %d)", state, every) << std::endl;
    } else {
      std::cerr << absl::StreamFormat("Reading state %d", state) << std::endl;
      
      auto tree_ss = std::stringstream{line.substr(line.find(" = ") + 3)};
      auto parser = Newick_parser{tree_ss};
      result.insert({state, newick_to_phylo_tree(parser.parse_tree(), tip_names)});
    }
  }

  return result;
}

static auto output_sequences(const Run& run, std::ostream& os) -> void {
  const auto& tree = run.tree();
  
  for (const auto& node : index_order_traversal(tree)) {
    if (tree.at(node).is_tip()) {
      os << absl::StreamFormat("    <sequence id=\"seq_%s\" spec=\"Sequence\" taxon=\"%s\" totalcount=\"4\" value=\"",
                               tree.at(node).name, tree.at(node).name);
      
      auto seq = view_of_sequence_at(tree, node);
      auto missing_sites = reconstruct_missing_sites_at(tree, node);
      
      auto missing_it = missing_sites.begin();
      auto in_missing = false;
      
      for (auto l = Site_index{0}; l != tree.num_sites(); ++l) {
        if (missing_it != missing_sites.end()) {
          const auto& [mi_start, mi_end] = *missing_it;
          if (not in_missing && l == mi_start) {
            in_missing = true;
          } else if (in_missing && l == mi_end) {
            in_missing = false;
            ++missing_it;
          }
        }
        
        if (in_missing) {
          os << 'N';
        } else {
          os << to_char(seq[l]);
        }
      }
      
      os << "\"/>\n";
    }
  }
}

auto export_beast_input(
    const Run& run,
    std::ostream& os,
    int64_t chain_length,
    int64_t log_every,
    int64_t tree_every)
    -> void {
  
  const auto& tree = run.tree();
  
  os << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n"
     << "\n"
     << "<!-- \n"
     << stamp_version_into_log_file{}
     << "-->\n"
     << "\n"
     << "<!-- BEAST2 v2.6.2 input file, modelled on run analyzed in \n"
     << "     LeMieux et al, \"Phylogenetic analysis of SARS-CoV-2 in Boston highlights\n"
     << "     the impact of superspreading events\", Science 371, 588 (2021)\n"
     << "     (https://dx.doi.org/10.1126/science.abe3261) -->\n"
     << "\n"
     << "<!-- This template is essentially what comes out of BEAUti2 after loading\n"
     << "     the input FASTA file and choosing the appropriate:\n"
     << "     * date format (yyyy-M-dd, after last '|')\n"
     << "     * site model (Gamma site model with 0 or 4 gamma categories, HKY subst model)\n"
     << "     * clock model (strict)\n"
     << "     * population prior (Coalescent Exponential Population)\n"
     << "     [leaving all parameters in their default settings] -->\n"
     << "\n";
  
  os << "<beast beautitemplate='Standard' beautistatus='' namespace=\"beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood\" required=\"\" version=\"2.6\">\n"
     << "\n";

  // Sequence data
  os << "  <data\n"
     << "   id=\"input_alignment\"\n"
     << "   spec=\"Alignment\"\n"
     << "   name=\"alignment\">\n";
  output_sequences(run, os);
  os << "  </data>\n";
  os << "\n";

  // Name aliases
  os << "  <map name=\"Uniform\" >beast.math.distributions.Uniform</map>\n";
  os << "  <map name=\"Exponential\" >beast.math.distributions.Exponential</map>\n";
  os << "  <map name=\"LogNormal\" >beast.math.distributions.LogNormalDistributionModel</map>\n";
  os << "  <map name=\"Normal\" >beast.math.distributions.Normal</map>\n";
  os << "  <map name=\"Beta\" >beast.math.distributions.Beta</map>\n";
  os << "  <map name=\"Gamma\" >beast.math.distributions.Gamma</map>\n";
  os << "  <map name=\"LaplaceDistribution\" >beast.math.distributions.LaplaceDistribution</map>\n";
  os << "  <map name=\"prior\" >beast.math.distributions.Prior</map>\n";
  os << "  <map name=\"InverseGamma\" >beast.math.distributions.InverseGamma</map>\n";
  os << "  <map name=\"OneOnX\" >beast.math.distributions.OneOnX</map>\n";
  os << "\n";

  os << "  <run id=\"mcmc\" spec=\"MCMC\" chainLength=\"" << chain_length << "\">\n";
  os << "    <state id=\"state\" spec=\"State\" storeEvery=\"5000\">\n";
  os << "      <tree id=\"Tree.t:input_alignment\" spec=\"beast.evolution.tree.Tree\" name=\"stateNode\">\n";
  os << "        <trait id=\"dateTrait.t:input_alignment\" spec=\"beast.evolution.tree.TraitSet\" "
     << "dateFormat=\"yyyy-M-dd\" traitname=\"date\" value=\"";
  auto first = true;
  for (const auto& node : index_order_traversal(tree)) {
    if (tree.at(node).is_tip()) {
      if (first) {
        first = false;
      } else {
        os << ",";
      }
      os << tree.at(node).name << "=" << to_iso_date(tree.at(node).t);
    }
  }
  os << "\">\n";
  os << "        <taxa id=\"TaxonSet.input_alignment\" spec=\"TaxonSet\">\n";
  os << "          <alignment idref=\"input_alignment\"/>\n";
  os << "        </taxa>\n";
  os << "      </trait>\n";
  os << "      <taxonset idref=\"TaxonSet.input_alignment\"/>\n";
  os << "    </tree>\n";
  os << "\n";

  // Parameter definitions and initial values
  if (run.mu_move_enabled()) {
    os << "    <parameter id=\"clockRate.c:input_alignment\" spec=\"parameter.RealParameter\" name=\"stateNode\">"
       << 1.0 << "</parameter>\n";
  }
  if (run.alpha_move_enabled()) {
    os << "    <parameter id=\"gammaShape.s:input_alignment\" spec=\"parameter.RealParameter\" name=\"stateNode\">"
       << 1.0 << "</parameter>\n";
  }
  os << "    <parameter id=\"kappa.s:input_alignment\" spec=\"parameter.RealParameter\" lower=\"0.0\" name=\"stateNode\">"
     << 2.0 << "</parameter>\n";
  if (run.final_pop_size_move_enabled()) {
    os << "    <parameter id=\"ePopSize.t:input_alignment\" spec=\"parameter.RealParameter\" name=\"stateNode\">"
       << 0.3 << "</parameter>\n";
  }
  if (run.pop_growth_rate_move_enabled()) {
    os << "    <parameter id=\"growthRate.t:input_alignment\" spec=\"parameter.RealParameter\" name=\"stateNode\">"
       << 3.0E-4 << "</parameter>\n";
  }
  os << "    <parameter id=\"freqParameter.s:input_alignment\" spec=\"parameter.RealParameter\" dimension=\"4\" lower=\"0.0\" name=\"stateNode\" upper=\"1.0\">"
     << 0.25 << "</parameter>\n";
  os << "  </state>\n";
  os << "\n";

  // Initial tree
  os << "  <init id=\"RandomTree.t:input_alignment\" spec=\"beast.evolution.tree.RandomTree\" estimate=\"false\" initial=\"@Tree.t:input_alignment\" taxa=\"@input_alignment\">\n";
  os << "    <populationModel id=\"ConstantPopulation0.t:input_alignment\" spec=\"ConstantPopulation\">\n";
  os << "      <parameter id=\"randomPopSize.t:input_alignment\" spec=\"parameter.RealParameter\" name=\"popSize\">1.0</parameter>\n";
  os << "    </populationModel>\n";
  os << "  </init>\n";
  os << "\n";

  // Posterior components
  os << "  <distribution id=\"posterior\" spec=\"util.CompoundDistribution\">\n";
  os << "    <distribution id=\"prior\" spec=\"util.CompoundDistribution\">\n";
  
  os << "      <distribution id=\"CoalescentExponential.t:input_alignment\" spec=\"Coalescent\">\n";
  os << "        <populationModel id=\"ExponentialGrowth.t:input_alignment\" spec=\"ExponentialGrowth\" growthRate=\"";
  if (run.pop_growth_rate_move_enabled()) {
    os << "@growthRate.t:input_alignment";
  } else {
      os << absl::StreamFormat("%g", run.pop_growth_rate()*365.0);  // per year!
  }
  os << "\" popSize=\"";
  if (run.final_pop_size_move_enabled()) {
      os << "@ePopSize.t:input_alignment";
  } else {
      os << absl::StreamFormat("%g", run.final_pop_size()/365.0);  // years!
  }
  os << "\"/>\n";
  os << "        <treeIntervals id=\"TreeIntervals.t:input_alignment\" spec=\"TreeIntervals\" tree=\"@Tree.t:input_alignment\"/>\n";
  os << "      </distribution>\n";

  if (run.mu_move_enabled()) {
    os << "      <prior id=\"ClockPrior.c:input_alignment\" name=\"distribution\" x=\"@clockRate.c:input_alignment\">\n";
    os << "        <Uniform id=\"Uniform.0\" name=\"distr\" upper=\"Infinity\"/>\n";
    os << "      </prior>\n";
  }

  if (run.final_pop_size_move_enabled()) {
    os << "      <prior id=\"ePopSizePrior.t:input_alignment\" name=\"distribution\" x=\"@ePopSize.t:input_alignment\">\n";
    os << "        <OneOnX id=\"OneOnX.1\" name=\"distr\"/>\n";
    os << "      </prior>\n";
  }

  os << "      <prior id=\"FrequenciesPrior.s:input_alignment\" name=\"distribution\" x=\"@freqParameter.s:input_alignment\">\n";
  os << "        <Uniform id=\"Uniform.3\" name=\"distr\"/>\n";
  os << "      </prior>\n";

  if (run.alpha_move_enabled()) {
    auto alpha_mean = 1.0;
    os << "      <prior id=\"GammaShapePrior.s:input_alignment\" name=\"distribution\" x=\"@gammaShape.s:input_alignment\">\n";
    os << "        <Exponential id=\"Exponential.0\" name=\"distr\">\n";
    os << "          <parameter id=\"RealParameter.0\" spec=\"parameter.RealParameter\" estimate=\"false\" name=\"mean\">"
       << alpha_mean << "</parameter>\n";
    os << "        </Exponential>\n";
    os << "      </prior>\n";
  }

  if (run.pop_growth_rate_move_enabled()) {
    auto growth_rate_mu = 0.001;  // per year - these defaults come from BEAUti2
    auto growth_rate_scale = 30.701135;  // per year - these defaults come from BEAUti2
    os << "      <prior id=\"GrowthRatePrior.t:input_alignment\" name=\"distribution\" x=\"@growthRate.t:input_alignment\">\n";
    os << "        <LaplaceDistribution id=\"LaplaceDistribution.0\" name=\"distr\">\n";
    os << "          <parameter id=\"RealParameter.3\" spec=\"parameter.RealParameter\" estimate=\"false\" name=\"mu\">"
       << absl::StreamFormat("%f", growth_rate_mu) << "</parameter>\n";
    os << "          <parameter id=\"RealParameter.4\" spec=\"parameter.RealParameter\" estimate=\"false\" name=\"scale\">"
       << absl::StreamFormat("%f", growth_rate_scale) << "</parameter>\n";
    os << "        </LaplaceDistribution>\n";
    os << "      </prior>\n";
  }

  auto log_kappa_mu = 1.0;  // these defaults come from BEAUti2
  auto log_kappa_sigma = 1.25;  // these defaults come from BEAUti2
  os << "      <prior id=\"KappaPrior.s:input_alignment\" name=\"distribution\" x=\"@kappa.s:input_alignment\">\n";
  os << "        <LogNormal id=\"LogNormalDistributionModel.0\" name=\"distr\">\n";
  os << "          <parameter id=\"RealParameter.1\" spec=\"parameter.RealParameter\" estimate=\"false\" name=\"M\">"
     << log_kappa_mu << "</parameter>\n";
  os << "          <parameter id=\"RealParameter.2\" spec=\"parameter.RealParameter\" estimate=\"false\" name=\"S\">"
     << log_kappa_sigma << "</parameter>\n";
  os << "        </LogNormal>\n";
  os << "      </prior>\n";

  os << "    </distribution>\n";

  os << "    <distribution id=\"likelihood\" spec=\"util.CompoundDistribution\" useThreads=\"true\">\n";
  os << "      <distribution id=\"treeLikelihood.input_alignment\" spec=\"ThreadedTreeLikelihood\" data=\"@input_alignment\" tree=\"@Tree.t:input_alignment\">\n";
  if (run.alpha_move_enabled()) {
    os << "        <siteModel id=\"SiteModel.s:input_alignment\" spec=\"SiteModel\" gammaCategoryCount=\"4\" shape=\"@gammaShape.s:input_alignment\">\n";
  } else {
    os << "        <siteModel id=\"SiteModel.s:input_alignment\" spec=\"SiteModel\" gammaCategoryCount=\"0\">\n";
  }
  os << "          <parameter id=\"mutationRate.s:input_alignment\" spec=\"parameter.RealParameter\" estimate=\"false\" name=\"mutationRate\">1.0</parameter>\n";  // Starting value
  os << "          <parameter id=\"proportionInvariant.s:input_alignment\" spec=\"parameter.RealParameter\" estimate=\"false\" lower=\"0.0\" name=\"proportionInvariant\" upper=\"1.0\">0.0</parameter>\n";
  os << "          <substModel id=\"hky.s:input_alignment\" spec=\"HKY\" kappa=\"@kappa.s:input_alignment\">\n";
  os << "            <frequencies id=\"estimatedFreqs.s:input_alignment\" spec=\"Frequencies\" frequencies=\"@freqParameter.s:input_alignment\"/>\n";
  os << "          </substModel>\n";
  os << "        </siteModel>\n";
  if (run.mu_move_enabled()) {
    os << "        <branchRateModel id=\"StrictClock.c:input_alignment\" spec=\"beast.evolution.branchratemodel.StrictClockModel\" clock.rate=\"@clockRate.c:input_alignment\"/>\n";
  } else {
    os << "        <branchRateModel id=\"StrictClock.c:input_alignment\" spec=\"beast.evolution.branchratemodel.StrictClockModel\" clock.rate=\"" << absl::StreamFormat("%g", run.mu()*365.0) << "\"/>\n";  // per year!
  }
  os << "      </distribution>\n";  // treeLikelihood
  os << "    </distribution>\n";  // likelihood

  os << "  </distribution>\n";  // posterior
  os << "\n";

  // Operators
  if (run.mu_move_enabled()) {
    os << "  <operator id=\"StrictClockRateScaler.c:input_alignment\" spec=\"ScaleOperator\" parameter=\"@clockRate.c:input_alignment\" scaleFactor=\"0.75\" weight=\"3.0\"/>\n";
  
    os << "  <operator id=\"strictClockUpDownOperator.c:input_alignment\" spec=\"UpDownOperator\" scaleFactor=\"0.75\" weight=\"3.0\">\n";
    os << "    <up idref=\"clockRate.c:input_alignment\"/>\n";
    os << "    <down idref=\"Tree.t:input_alignment\"/>\n";
    os << "  </operator>\n";
  }

  if (run.alpha_move_enabled()) {
    os << "  <operator id=\"gammaShapeScaler.s:input_alignment\" spec=\"ScaleOperator\" parameter=\"@gammaShape.s:input_alignment\" scaleFactor=\"0.5\" weight=\"0.1\"/>\n";
  }

  os << "  <operator id=\"KappaScaler.s:input_alignment\" spec=\"ScaleOperator\" parameter=\"@kappa.s:input_alignment\" scaleFactor=\"0.5\" weight=\"0.1\"/>\n";

  os << "  <operator id=\"CoalescentExponentialTreeScaler.t:input_alignment\" spec=\"ScaleOperator\" scaleFactor=\"0.5\" tree=\"@Tree.t:input_alignment\" weight=\"3.0\"/>\n";

  os << "  <operator id=\"CoalescentExponentialTreeRootScaler.t:input_alignment\" spec=\"ScaleOperator\" rootOnly=\"true\" scaleFactor=\"0.5\" tree=\"@Tree.t:input_alignment\" weight=\"3.0\"/>\n";

  os << "  <operator id=\"CoalescentExponentialUniformOperator.t:input_alignment\" spec=\"Uniform\" tree=\"@Tree.t:input_alignment\" weight=\"30.0\"/>\n";

  os << "  <operator id=\"CoalescentExponentialSubtreeSlide.t:input_alignment\" spec=\"SubtreeSlide\" tree=\"@Tree.t:input_alignment\" weight=\"15.0\"/>\n";

  os << "  <operator id=\"CoalescentExponentialNarrow.t:input_alignment\" spec=\"Exchange\" tree=\"@Tree.t:input_alignment\" weight=\"15.0\"/>\n";

  os << "  <operator id=\"CoalescentExponentialWide.t:input_alignment\" spec=\"Exchange\" isNarrow=\"false\" tree=\"@Tree.t:input_alignment\" weight=\"3.0\"/>\n";

  os << "  <operator id=\"CoalescentExponentialWilsonBalding.t:input_alignment\" spec=\"WilsonBalding\" tree=\"@Tree.t:input_alignment\" weight=\"3.0\"/>\n";
  
  if (run.final_pop_size_move_enabled()) {
      os << "  <operator id=\"ePopSizeScaler.t:input_alignment\" spec=\"ScaleOperator\" parameter=\"@ePopSize.t:input_alignment\" scaleFactor=\"0.75\" weight=\"3.0\"/>\n";
  }

  if (run.pop_growth_rate_move_enabled()) {
      os << "  <operator id=\"GrowthRateRandomWalk.t:input_alignment\" spec=\"RealRandomWalkOperator\" parameter=\"@growthRate.t:input_alignment\" weight=\"3.0\" windowSize=\"1.0\"/>\n";
  }

  os << "  <operator id=\"FrequenciesExchanger.s:input_alignment\" spec=\"DeltaExchangeOperator\" delta=\"0.01\" weight=\"0.1\">\n";
  os << "    <parameter idref=\"freqParameter.s:input_alignment\"/>\n";
  os << "  </operator>\n";
  os << "\n";

  // Loggers
  os << "  <logger id=\"tracelog\" spec=\"Logger\" fileName=\"output.log\" logEvery=\"" << log_every << "\" model=\"@posterior\">\n";
  os << "    <log idref=\"posterior\"/>\n";
  os << "    <log idref=\"likelihood\"/>\n";
  os << "    <log idref=\"prior\"/>\n";
  os << "    <log idref=\"treeLikelihood.input_alignment\"/>\n";
  os << "    <log id=\"TreeHeight.t:input_alignment\" spec=\"beast.evolution.tree.TreeHeightLogger\" tree=\"@Tree.t:input_alignment\"/>\n";
  if (run.mu_move_enabled()) {
    os << "    <log idref=\"clockRate.c:input_alignment\"/>\n";
  }
  if (run.alpha_move_enabled()) {
    os << "    <log idref=\"gammaShape.s:input_alignment\"/>\n";
  }
  os << "    <log idref=\"kappa.s:input_alignment\"/>\n";
  os << "    <log idref=\"CoalescentExponential.t:input_alignment\"/>\n";
  if (run.final_pop_size_move_enabled()) {
    os << "    <log idref=\"ePopSize.t:input_alignment\"/>\n";
  }
  if (run.pop_growth_rate_move_enabled()) {
    os << "    <log idref=\"growthRate.t:input_alignment\"/>\n";
  }
  os << "    <log idref=\"freqParameter.s:input_alignment\"/>\n";
  os << "  </logger>\n";
  os << "\n";
  
  os << "  <logger id=\"screenlog\" spec=\"Logger\" logEvery=\"1000\">\n";
  os << "    <log idref=\"posterior\"/>\n";
  os << "    <log idref=\"likelihood\"/>\n";
  os << "    <log idref=\"prior\"/>\n";
  os << "  </logger>\n";
  os << "\n";

  os << "  <logger id=\"treelog.t:input_alignment\" spec=\"Logger\" fileName=\"output.trees\" logEvery=\"" << tree_every << "\" mode=\"tree\">\n";
  os << "    <log id=\"TreeWithMetaDataLogger.t:input_alignment\" spec=\"beast.evolution.tree.TreeWithMetaDataLogger\" tree=\"@Tree.t:input_alignment\"/>\n";
  os << "  </logger>\n";
  os << "\n";
  
  os << "  <operatorschedule id=\"OperatorSchedule\" spec=\"OperatorSchedule\"/>\n";
  os << "\n";
  
  os << "</run>\n";
  os << "</beast>\n";
}

}  // namespace delphy
