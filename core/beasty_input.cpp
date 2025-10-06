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

static auto export_beast_2_6_2_input(
    const Run& run,
    std::ostream& os,
    int64_t chain_length,
    int64_t log_every,
    int64_t tree_every)
    -> void {
  
  const auto& tree = run.tree();

  if (run.mpox_hack_enabled()) {
    std::cerr << "ERROR: BEAST2 input XML generation not currently supported when enabling APOBEC3 treatment (mpox)\n";
    os << "ERROR: BEAST2 input XML generation not currently supported when enabling APOBEC3 treatment (mpox)\n";
    return;
  }
  
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
      os << tree.at(node).name << "=" << to_iso_date(0.5*(tree.at(node).t_min + tree.at(node).t_max));
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
  
  os << "      <distribution id=\"Coalescent.t:input_alignment\" spec=\"Coalescent\">\n";

  const auto& pop_model = run.pop_model();
  if (typeid(pop_model) == typeid(Exp_pop_model)) {
    const auto& exp_pop_model = static_cast<const Exp_pop_model&>(run.pop_model());
    auto final_pop_size = exp_pop_model.pop_at_t0();
    auto pop_growth_rate = exp_pop_model.growth_rate();
    
    os << "        <populationModel id=\"ExponentialGrowth.t:input_alignment\" spec=\"ExponentialGrowth\" growthRate=\"";
    if (run.pop_growth_rate_move_enabled()) {
      os << "@growthRate.t:input_alignment";
    } else {
      os << absl::StreamFormat("%g", pop_growth_rate*365.0);  // per year!
    }
    os << "\" popSize=\"";
    if (run.final_pop_size_move_enabled()) {
      os << "@ePopSize.t:input_alignment";
    } else {
      os << absl::StreamFormat("%g", final_pop_size/365.0);  // years!
    }
    os << "\"/>\n";
    
  } else if (typeid(pop_model) == typeid(Skygrid_pop_model)) {
    // Make invalid XML tag on purpose to stop this BEAST2 XML file from running (with God knows what population model...)
    std::cerr << "ERROR: BEAST2 doesn't implement a Skygrid model!\n";
    os << absl::StreamFormat("<ERROR>BEAST2 doesn't implement a Skygrid model</ERROR>");
    
  } else {
    // Make invalid XML tag on purpose to stop this BEAST2 XML file from running (with God knows what population model...)
    std::cerr << "ERROR: Unrecognized population model type " << typeid(pop_model).name() << '\n';
    os << absl::StreamFormat("<ERROR>UNRECOGNIZED POPULATION MODEL TYPE: %s</ERROR>", typeid(pop_model).name());
  }
  
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

  // Tip-date sampling (prior)
  // See https://www.beast2.org/2015/06/09/sampling-tip-dates.html
  auto first_uncertain_tip = k_no_node;
  for (const auto& node : index_order_traversal(tree)) {
    if (tree.at(node).is_tip() && tree.at(node).t_min != tree.at(node).t_max) {
      if (first_uncertain_tip == k_no_node) {
        first_uncertain_tip = node;
      }
      os << "      <distribution id=\"tip-dist." << tree.at(node).name << "\" spec=\"beast.math.distributions.MRCAPrior\" tipsonly=\"true\" tree=\"@Tree.t:input_alignment\">\n";
      os << "        <taxonset id=\"tip-taxonset." << tree.at(node).name << "\" spec=\"TaxonSet\">\n";
      os << "          <taxon id=\"" << tree.at(node).name << "\" spec=\"Taxon\"/>\n";
      os << "        </taxonset>\n";
      os << "        <Uniform id=\"tip-uniform." << tree.at(node).name << "\" name=\"distr\" "
         << "lower=\"" << absl::StreamFormat("%.5f", to_linear_year(tree.at(node).t_min)) << "\" "
         << "upper=\"" << absl::StreamFormat("%.5f", to_linear_year(tree.at(node).t_max)) << "\"/>\n";
      os << "      </distribution>\n";
    }
  }
  
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

  // Tip-date sampling (operators)
  auto num_uncertain_tips = 0;
  for (const auto& node : index_order_traversal(tree)) {
    if (tree.at(node).is_tip() && tree.at(node).t_min != tree.at(node).t_max) {
      ++num_uncertain_tips;
    }
  }
  if (num_uncertain_tips > 0) {
    auto tot_weight_tip_date_sampling = 10.0;
    auto per_tip_weight_tip_date_sampling = tot_weight_tip_date_sampling / num_uncertain_tips;
    auto max_tip_date_sampling_window_size = 1.0 / (tree.num_sites()*run.mu()*365.0);
    // ^^ max chosen so that a branch ending at a tip with a fixed number of mutations is rarely
    //    overstretched or overcompressed
    
    for (const auto& node : index_order_traversal(tree)) {
      if (tree.at(node).is_tip() && tree.at(node).t_min != tree.at(node).t_max) {
        auto window_size = std::min(max_tip_date_sampling_window_size,
                                    double{tree.at(node).t_max - tree.at(node).t_min} / 4);
        os << "  <operator id=\"tip-operator." << tree.at(node).name << "\" "
           << "windowSize=\"" << window_size << "\" "
           << "spec=\"TipDatesRandomWalker\" "
           << "taxonset=\"@tip-taxonset." << tree.at(node).name << "\" "
           << "tree=\"@Tree.t:input_alignment\" "
           << "weight=\"" << per_tip_weight_tip_date_sampling << "\"/>\n";
      }
    }
  }
  
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
  os << "    <log idref=\"Coalescent.t:input_alignment\"/>\n";
  if (typeid(pop_model) == typeid(Exp_pop_model)) {
    if (run.final_pop_size_move_enabled()) {
      os << "    <log idref=\"ePopSize.t:input_alignment\"/>\n";
    }
    if (run.pop_growth_rate_move_enabled()) {
      os << "    <log idref=\"growthRate.t:input_alignment\"/>\n";
    }
  }
  os << "    <log idref=\"freqParameter.s:input_alignment\"/>\n";
  if (first_uncertain_tip != k_no_node) {
    os << "    <log idref=\"tip-dist." << tree.at(first_uncertain_tip).name << "\"/>\n";
  }
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

static auto output_sequences_X_10_5(const Run& run, std::ostream& os) -> void {
  const auto& tree = run.tree();

  auto num_tips = Node_index{(tree.size() + 1) / 2};
  os << absl::StreamFormat("  <!-- The list of taxa to be analysed (can also include dates/ages).          -->\n")
     << absl::StreamFormat("  <!-- ntax=%-56d -->\n", num_tips);
  
  os << absl::StreamFormat("  <taxa id=\"taxa\">\n");
  for (const auto& node : index_order_traversal(tree)) {
    if (tree.at(node).is_tip()) {
      os << absl::StreamFormat("    <taxon id=\"%s\">\n", tree.at(node).name);
      os << absl::StreamFormat("      <date value=\"%g\" direction=\"forwards\" units=\"years\"",  // Partial line
                               0.5 * (to_linear_year(tree.at(node).t_min) + to_linear_year(tree.at(node).t_max)));

      // Tip date uncertainty
      if (tree.at(node).t_min != tree.at(node).t_max) {
        os << absl::StreamFormat(" uncertainty=\"%g\"",
                                 to_linear_year(tree.at(node).t_max) - to_linear_year(tree.at(node).t_min));
      }
      
      os << absl::StreamFormat("/>\n");
      os << absl::StreamFormat("    </taxon>\n");
    }
  }
  os << absl::StreamFormat("  </taxa>\n")
     << absl::StreamFormat("\n")
     << absl::StreamFormat("\n")
     << absl::StreamFormat("\n");

  os << absl::StreamFormat("  <!-- The sequence alignment (each sequence refers to a taxon above).         -->\n")
     << absl::StreamFormat("  <!-- %-61s -->\n", absl::StrFormat("ntax = %d nchar = %d",
                                                                 num_tips, tree.num_sites()));
  os << absl::StreamFormat("  <alignment id=\"alignment\" dataType=\"nucleotide\">\n");
  for (const auto& node : index_order_traversal(tree)) {
    if (tree.at(node).is_tip()) {
      os << absl::StreamFormat("    <sequence>\n");
      os << absl::StreamFormat("      <taxon idref=\"%s\">\n", tree.at(node).name);
      os << absl::StreamFormat("      ");  // No newline here
      
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
      
      os << absl::StreamFormat("    </sequence>\n");
    }
  }
  os << absl::StreamFormat("  </alignment>\n");
}

static auto export_beast_X_10_5_0_input(
    const Run& run,
    std::ostream& os,
    int64_t chain_length,
    int64_t log_every,
    int64_t tree_every)
    -> void {
  
  const auto& tree = run.tree();

  if (run.mpox_hack_enabled()) {
    std::cerr << "ERROR: BEAST input XML generation not currently supported when enabling APOBEC3 treatment (mpox)\n";
    os << "ERROR: BEAST input XML generation not currently supported when enabling APOBEC3 treatment (mpox)\n";
    return;
  }
  
  // Gather basic info about tree
  auto num_tips = Node_index{(tree.size() + 1) / 2};
  
  auto beast_t0 = -std::numeric_limits<double>::max();  // t0 = time of latest tip
  for (const auto& node : index_order_traversal(tree)) {
    if (tree.at(node).is_tip()) {
      beast_t0 = std::max(beast_t0, tree.at(node).t);
    }
  }
  const auto& pop_model = run.pop_model();
  
  auto has_tip_uncertainty = false;
  for (const auto& node : index_order_traversal(tree)) {
    if (tree.at(node).is_tip() && tree.at(node).t_min != tree.at(node).t_max) {
      has_tip_uncertainty = true;
      break;
    }
  }

  // Off we go!
  os << "<?xml version=\"1.0\" standalone=\"yes\"?>\n"
     << "\n"
     << "<!-- \n"
     << stamp_version_into_log_file{}
     << "-->\n"
     << "\n"
     << "<!-- BEAST X 10.5.0 input file, modelled on run analyzed in \n"
     << "     LeMieux et al, \"Phylogenetic analysis of SARS-CoV-2 in Boston highlights\n"
     << "     the impact of superspreading events\", Science 371, 588 (2021)\n"
     << "     (https://dx.doi.org/10.1126/science.abe3261) -->\n"
     << "\n"
     << "<!-- This template is essentially what comes out of BEAUti X after loading\n"
     << "     the input FASTA file and choosing the appropriate:\n"
     << "     * date format (yyyy-M-dd, after last '|')\n"
     << "     * site model (Gamma site model with 0 or 4 gamma categories, HKY subst model)\n"
     << "     * clock model (strict)\n"
     << "     * population prior (Coalescent: Exponential Population or Coalescent: Hamiltonian Monte Carlo SkyGrid)\n"
     << "     [leaving all parameters in their default settings] -->\n"
     << "\n";
  
  os << "<beast version=\"10.5.0-beta5\">\n"
     << "\n"
     << "\n";

  // Sequence data
  output_sequences_X_10_5(run, os);
  os << absl::StreamFormat("\n")
     << absl::StreamFormat("\n")
     << absl::StreamFormat("\n")
     << absl::StreamFormat("  <!-- The unique patterns from 1 to end                                       -->\n")
     << absl::StreamFormat("  <!-- npatterns=<DELPHY DOESN'T COUNT THESE>                                  -->\n")
     << absl::StreamFormat("  <patterns id=\"patterns\" from=\"1\" strip=\"false\">\n")
     << absl::StreamFormat("    <alignment idref=\"alignment\"/>\n")
     << absl::StreamFormat("  </patterns>\n")
     << absl::StreamFormat("\n")
     << absl::StreamFormat("\n");

  // Starting tree initialized from a random coalescent simulation with constant-size or exponentially growing population
  if (typeid(pop_model) == typeid(Exp_pop_model)) {
    const auto& exp_pop_model = static_cast<const Exp_pop_model&>(run.pop_model());
    auto final_pop_size = exp_pop_model.pop_at_t0();
    auto pop_growth_rate = exp_pop_model.growth_rate();
    
    os << absl::StreamFormat("  <!-- A prior assumption that the population size has grown exponentially     -->\n")
       << absl::StreamFormat("  <!-- throughout the time spanned by the genealogy.                           -->\n")
       << absl::StreamFormat("  <exponentialGrowth id=\"exponential\" units=\"years\">\n")
       << absl::StreamFormat("    <populationSize>\n")
       << absl::StreamFormat("      <parameter id=\"exponential.popSize\" value=\"%g\" lower=\"0.0\"/>\n",
                             run.final_pop_size_move_enabled()
                             ? 1.0                    // Will infer...
                             : final_pop_size/365.0)  // ... vs. fixed  (units: years)
       << absl::StreamFormat("    </populationSize>\n")
       << absl::StreamFormat("    <growthRate>\n")
       << absl::StreamFormat("      <parameter id=\"exponential.growthRate\" value=\"%g\"/>\n",
                             run.final_pop_size_move_enabled()
                             ? 0.0                    // Will infer...
                             : pop_growth_rate*365.0)  // ... vs. fixed  (units: e-foldings per year)
       << absl::StreamFormat("    </growthRate>\n")
       << absl::StreamFormat("  </exponentialGrowth>\n")
       << absl::StreamFormat("\n")
       << absl::StreamFormat("\n")
       << absl::StreamFormat("  <!-- Generate a random starting tree under the coalescent process            -->\n")
       << absl::StreamFormat("  <coalescentSimulator id=\"startingTree\">\n")
       << absl::StreamFormat("    <taxa idref=\"taxa\"/>\n")
       << absl::StreamFormat("    <exponentialGrowth idref=\"exponential\"/>\n")
       << absl::StreamFormat("  </coalescentSimulator>\n")
       << absl::StreamFormat("\n")
       << absl::StreamFormat("\n");

  } else {
    os << absl::StreamFormat("  <!-- This is a simple constant population size coalescent model              -->\n")
       << absl::StreamFormat("  <!-- that is used to generate an initial tree for the chain.                 -->\n")
       << absl::StreamFormat("  <constantSize id=\"initialDemo\" units=\"years\">\n")
       << absl::StreamFormat("    <populationSize>\n")
       << absl::StreamFormat("      <parameter id=\"initialDemo.popSize\" value=\"100.0\"/>\n")
       << absl::StreamFormat("    </populationSize>\n")
       << absl::StreamFormat("  </constantSize>\n")
       << absl::StreamFormat("\n")
       << absl::StreamFormat("\n")
       << absl::StreamFormat("  <!-- Generate a random starting tree under the coalescent process            -->\n")
       << absl::StreamFormat("  <coalescentSimulator id=\"startingTree\">\n")
       << absl::StreamFormat("    <taxa idref=\"taxa\"/>\n")
       << absl::StreamFormat("    <constantSize idref=\"initialDemo\"/>\n")
       << absl::StreamFormat("  </coalescentSimulator>\n")
       << absl::StreamFormat("\n")
       << absl::StreamFormat("\n");
  }

  // Tree model and basic statistics
  os << absl::StreamFormat("  <!-- Generate a tree model                                                   -->\n")
     << absl::StreamFormat("  <treeModel id=\"treeModel\">\n")
     << absl::StreamFormat("    <coalescentTree idref=\"startingTree\"/>\n")
     << absl::StreamFormat("    <rootHeight>\n")
     << absl::StreamFormat("      <parameter id=\"treeModel.rootHeight\"/>")
     << absl::StreamFormat("    </rootHeight>\n")
     << absl::StreamFormat("    <nodeHeights internalNodes=\"true\">\n")
     << absl::StreamFormat("      <parameter id=\"treeModel.internalNodeHeights\"/>\n")
     << absl::StreamFormat("    </nodeHeights>\n")
     << absl::StreamFormat("    <nodeHeights internalNodes=\"true\" rootNode=\"true\">\n")
     << absl::StreamFormat("      <parameter id=\"treeModel.allInternalNodeHeights\"/>\n")
     << absl::StreamFormat("    </nodeHeights>\n");

  if (has_tip_uncertainty) {
    os << absl::StreamFormat("\n");
    os << absl::StreamFormat("    <!-- START Tip date sampling                                                 -->\n");
    for (const auto& node : index_order_traversal(tree)) {
      if (tree.at(node).is_tip() && tree.at(node).t_min != tree.at(node).t_max) {
        os << absl::StreamFormat("    <leafHeight taxon=\"%s\">\n",
                                 tree.at(node).name)
           << absl::StreamFormat("      <parameter id=\"age(%s)\"/>\n",
                                 tree.at(node).name)
           << absl::StreamFormat("    </leafHeight>\n");
      }
    }
    os << absl::StreamFormat("\n");
    os << absl::StreamFormat("    <!-- END Tip date sampling                                                   -->\n");
    os << absl::StreamFormat("\n");
  }
  
  os << absl::StreamFormat("  </treeModel>\n")
     << absl::StreamFormat("\n")
     << absl::StreamFormat("  <!-- Statistic for height of the root of the tree                            -->\n")
     << absl::StreamFormat("  <treeHeightStatistic id=\"rootHeight\">\n")
     << absl::StreamFormat("    <treeModel idref=\"treeModel\"/>\n")
     << absl::StreamFormat("  </treeHeightStatistic>\n")
     << absl::StreamFormat("\n")
     << absl::StreamFormat("  <!-- Statistic for sum of the branch lengths of the tree (tree length)       -->\n")
     << absl::StreamFormat("  <treeLengthStatistic id=\"treeLength\">\n")
     << absl::StreamFormat("    <treeModel idref=\"treeModel\"/>\n")
     << absl::StreamFormat("  </treeLengthStatistic>\n")
     << absl::StreamFormat("\n")
     << absl::StreamFormat("  <!-- Statistic for time of most recent common ancestor of tree               -->\n")
     << absl::StreamFormat("  <tmrcaStatistic id=\"age(root)\" absolute=\"true\">\n")
     << absl::StreamFormat("    <treeModel idref=\"treeModel\"/>\n")
     << absl::StreamFormat("  </tmrcaStatistic>\n")
     << absl::StreamFormat("\n")
     << absl::StreamFormat("\n")
     << absl::StreamFormat("\n");

  // Skygrid flexible population curve
  if (typeid(pop_model) == typeid(Skygrid_pop_model)) {
    const auto& skygrid_pop_model = static_cast<const Skygrid_pop_model&>(run.pop_model());
    
    os << absl::StreamFormat("  <!-- Generate a gmrfSkyGridLikelihood for the Bayesian SkyGrid process       -->\n")
       << absl::StreamFormat("  <gmrfSkyGridLikelihood id=\"skygrid\">\n")
       << absl::StreamFormat("    <populationSizes>\n")
       << absl::StreamFormat("      \n")
       << absl::StreamFormat("      <!-- skygrid.logPopSize is in log units unlike other popSize                 -->\n")
       << absl::StreamFormat("      <parameter id=\"skygrid.logPopSize\" dimension=\"%d\" value=\"1.0\"/>\n",
                             skygrid_pop_model.M() + 1)
       << absl::StreamFormat("    </populationSizes>\n")
       << absl::StreamFormat("    <precisionParameter>\n")
       << absl::StreamFormat("      <parameter id=\"skygrid.precision\" value=\"0.1\" lower=\"0.0\"/>\n")
       << absl::StreamFormat("    </precisionParameter>\n")
       << absl::StreamFormat("    <numGridPoints>\n")
       << absl::StreamFormat("      <parameter id=\"skygrid.numGridPoints\" value=\"%d\"/>\n",
                             skygrid_pop_model.M())
       << absl::StreamFormat("    </numGridPoints>\n")
       << absl::StreamFormat("    <cutOff>\n")
       << absl::StreamFormat("      <parameter id=\"skygrid.cutOff\" value=\"%.5f\"/>",
                             (beast_t0 - skygrid_pop_model.x(0)) / 365.0)
       << absl::StreamFormat("    </cutOff>\n")
       << absl::StreamFormat("    <populationTree>\n")
       << absl::StreamFormat("      <treeModel idref=\"treeModel\"/>\n")
       << absl::StreamFormat("    </populationTree>\n")
       << absl::StreamFormat("  </gmrfSkyGridLikelihood>\n");

    os << absl::StreamFormat("  <gammaPrior id=\"skygrid.precision.prior\" shape=\"%.5f\" scale=\"%.5f\" offset=\"0.0\">\n",
                             run.skygrid_tau_prior_alpha(),
                             1.0 / run.skygrid_tau_prior_beta())
       << absl::StreamFormat("    <parameter idref=\"skygrid.precision\"/>\n")
       << absl::StreamFormat("  </gammaPrior>\n")
       << absl::StreamFormat("  <gmrfSkyrideGradient id=\"gmrfGradientPop\" wrtParameter=\"logPopulationSizes\">\n")
       << absl::StreamFormat("    <gmrfSkyGridLikelihood idref=\"skygrid\"/>\n")
       << absl::StreamFormat("  </gmrfSkyrideGradient>\n")
       << absl::StreamFormat("  <compoundParameter id=\"skygrid.parameters\">\n");
    if (run.skygrid_tau_move_enabled()) {
      os << absl::StreamFormat("    <parameter idref=\"skygrid.precision\"/>\n");
    }
    os << absl::StreamFormat("    <parameter idref=\"skygrid.logPopSize\"/>\n")
       << absl::StreamFormat("  </compoundParameter>\n");
    if (run.skygrid_tau_move_enabled()) {
      os << absl::StreamFormat("  <gmrfSkyrideGradient id=\"gmrfGradientPrec\" wrtParameter=\"precision\">\n")
         << absl::StreamFormat("    <gmrfSkyGridLikelihood idref=\"skygrid\"/>\n")
         << absl::StreamFormat("  </gmrfSkyrideGradient>\n")
         << absl::StreamFormat("  <jointGradient id=\"joint.skygrid.precision\">\n")
         << absl::StreamFormat("    <gmrfSkyrideGradient idref=\"gmrfGradientPrec\"/>\n")
         << absl::StreamFormat("    <gradient>\n")
         << absl::StreamFormat("      <gammaPrior idref=\"skygrid.precision.prior\"/>\n")
         << absl::StreamFormat("      <parameter idref=\"skygrid.precision\"/>\n")
         << absl::StreamFormat("    </gradient>\n")
         << absl::StreamFormat("  </jointGradient>\n");
    }
    os << absl::StreamFormat("  <compoundGradient id=\"full.skygrid.gradient\">\n");
    if (run.skygrid_tau_move_enabled()) {
      os << absl::StreamFormat("    <jointGradient idref=\"joint.skygrid.precision\"/>\n");
    }
    os << absl::StreamFormat("    <gmrfSkyrideGradient idref=\"gmrfGradientPop\"/>\n")
       << absl::StreamFormat("  </compoundGradient>\n")
       << absl::StreamFormat("\n");
    
  } else if (typeid(pop_model) == typeid(Exp_pop_model)) {
    
    os << absl::StreamFormat("  <!-- Generate a coalescent likelihood                                        -->\n")
       << absl::StreamFormat("  <coalescentLikelihood id=\"coalescent\">\n")
       << absl::StreamFormat("    <model>\n")
       << absl::StreamFormat("      <exponentialGrowth idref=\"exponential\"/>\n")
       << absl::StreamFormat("    </model>\n")
       << absl::StreamFormat("    <populationTree>\n")
       << absl::StreamFormat("      <treeModel idref=\"treeModel\"/>\n")
       << absl::StreamFormat("    </populationTree>\n")
       << absl::StreamFormat("  </coalescentLikelihood>\n")
       << absl::StreamFormat("\n")
       << absl::StreamFormat("\n");
      
  } else {
    // Make invalid XML tag on purpose to stop this BEAST X XML file from running (with God knows what population model...)
    std::cerr << "ERROR: Unrecognized population model type " << typeid(pop_model).name() << '\n';
    os << absl::StreamFormat("<ERROR>UNRECOGNIZED POPULATION MODEL TYPE: %s</ERROR>", typeid(pop_model).name());
  }

  // Mutation rates
  os << absl::StreamFormat("  <!-- The strict clock (Uniform rates across branches)                        -->\n")
     << absl::StreamFormat("  <strictClockBranchRates id=\"branchRates\">\n")
     << absl::StreamFormat("    <rate>\n")
     << absl::StreamFormat("      <parameter id=\"clock.rate\" value=\"%g\" lower=\"0.0\"/>\n",
                           run.mu_move_enabled()
                           ? 0.001             // Will infer...
                           : run.mu()*365.0)   // ... vs. fixed (units: mutations per site per year)
     << absl::StreamFormat("    </rate>\n")
     << absl::StreamFormat("  </strictClockBranchRates>\n")
     << absl::StreamFormat("\n")
     << absl::StreamFormat("  <rateStatistic id=\"meanRate\" name=\"meanRate\" mode=\"mean\" internal=\"true\" external=\"true\">\n")
     << absl::StreamFormat("    <treeModel idref=\"treeModel\"/>\n")
     << absl::StreamFormat("    <strictClockBranchRates idref=\"branchRates\"/>\n")
     << absl::StreamFormat("  </rateStatistic>\n")
     << absl::StreamFormat("\n")
     << absl::StreamFormat("\n");

  // Substitution model
  os << absl::StreamFormat("  <!-- The HKY substitution model (Hasegawa, Kishino & Yano, 1985)             -->\n")
     << absl::StreamFormat("  <HKYModel id=\"hky\">\n")
     << absl::StreamFormat("    <frequencies>\n")
     << absl::StreamFormat("      <frequencyModel dataType=\"nucleotide\">\n")
     << absl::StreamFormat("        <frequencies>\n")
     << absl::StreamFormat("          <parameter id=\"frequencies\" value=\"0.25 0.25 0.25 0.25\"/>\n")
     << absl::StreamFormat("        </frequencies>\n")
     << absl::StreamFormat("      </frequencyModel>\n")
     << absl::StreamFormat("    </frequencies>\n")
     << absl::StreamFormat("    <kappa>\n")
     << absl::StreamFormat("      <parameter id=\"kappa\" value=\"2.0\" lower=\"0.0\"/>\n")
     << absl::StreamFormat("    </kappa>\n")
     << absl::StreamFormat("  </HKYModel>\n")
     << absl::StreamFormat("\n");


  // Site model
  os << absl::StreamFormat("  <!-- site model                                                              -->\n")
     << absl::StreamFormat("  <siteModel id=\"siteModel\">\n")
     << absl::StreamFormat("    <substitutionModel>\n")
     << absl::StreamFormat("      <HKYModel idref=\"hky\"/>\n")
     << absl::StreamFormat("    </substitutionModel>\n");

  if (run.alpha_move_enabled()) {
    os << absl::StreamFormat("    <gammaShape gammaCategories=\"4\">\n")
       << absl::StreamFormat("      <parameter id=\"alpha\" value=\"1.0\" lower=\"0.0\"/>\n")
       << absl::StreamFormat("    </gammaShape>\n");
  }

  os << absl::StreamFormat("  </siteModel>\n")
     << absl::StreamFormat("\n")
     << absl::StreamFormat("  <!--                                                                         -->\n")
     << absl::StreamFormat("  <statistic id=\"mu\" name=\"mu\">\n")
     << absl::StreamFormat("    <siteModel idref=\"siteModel\"/>\n")
     << absl::StreamFormat("  </statistic>\n")
     << absl::StreamFormat("\n")
     << absl::StreamFormat("\n");
	

  // Tree likelihood
  os << absl::StreamFormat("  <!-- Likelihood for tree given sequence data                                 -->\n")
      // Below, useAmbiguities = false just means treating things like R and Y as N instead of A|G and C|T, respectively.
      // Missing data is still treated as you would expect.
      // See https://groups.google.com/g/beast-users/c/rWqz6FcTdKU/m/c20Q6djYAgAJ .
     << absl::StreamFormat("  <treeDataLikelihood id=\"treeLikelihood\" useAmbiguities=\"false\" usePreOrder=\"false\">\n")
     << absl::StreamFormat("    <partition>\n")
     << absl::StreamFormat("      <patterns idref=\"patterns\"/>\n")
     << absl::StreamFormat("      <siteModel idref=\"siteModel\"/>\n")
     << absl::StreamFormat("    </partition>\n")
     << absl::StreamFormat("    <treeModel idref=\"treeModel\"/>\n")
     << absl::StreamFormat("    <strictClockBranchRates idref=\"branchRates\"/>\n")
     << absl::StreamFormat("  </treeDataLikelihood>\n")
     << absl::StreamFormat("\n")
     << absl::StreamFormat("\n");
  
  // Operators
  os << absl::StreamFormat("  <!-- Define operators                                                        -->\n")
     << absl::StreamFormat("  <operators id=\"operators\" optimizationSchedule=\"log\">\n")
     << absl::StreamFormat("    <scaleOperator scaleFactor=\"0.75\" weight=\"1\">\n")
     << absl::StreamFormat("      <parameter idref=\"kappa\"/>\n")
     << absl::StreamFormat("    </scaleOperator>\n")
     << absl::StreamFormat("    <deltaExchange delta=\"0.01\" weight=\"1\">\n")
     << absl::StreamFormat("      <parameter idref=\"frequencies\"/>\n")
     << absl::StreamFormat("    </deltaExchange>\n");

  if (run.alpha_move_enabled()) {
    os << absl::StreamFormat("    <scaleOperator scaleFactor=\"0.75\" weight=\"1\">\n")
       << absl::StreamFormat("      <parameter idref=\"alpha\"/>\n")
       << absl::StreamFormat("    </scaleOperator>\n");
  }
  if (run.mu_move_enabled()) {
    os << absl::StreamFormat("    <scaleOperator scaleFactor=\"0.75\" weight=\"3\">\n")
       << absl::StreamFormat("      <parameter idref=\"clock.rate\"/>\n")
       << absl::StreamFormat("    </scaleOperator>\n")
       << absl::StreamFormat("    <upDownOperator scaleFactor=\"0.75\" weight=\"3\">\n")
       << absl::StreamFormat("      <up>\n")
       << absl::StreamFormat("        <parameter idref=\"treeModel.allInternalNodeHeights\"/>\n")
       << absl::StreamFormat("      </up>\n")
       << absl::StreamFormat("      <down>\n")
       << absl::StreamFormat("        <parameter idref=\"clock.rate\"/>\n")
       << absl::StreamFormat("      </down>\n")
       << absl::StreamFormat("    </upDownOperator>\n");
  }
  
  os << absl::StreamFormat("    <subtreeLeap size=\"1.0\" weight=\"%d\">\n",
                           num_tips)
     << absl::StreamFormat("      <treeModel idref=\"treeModel\"/>\n")
     << absl::StreamFormat("    </subtreeLeap>\n")
     << absl::StreamFormat("    <fixedHeightSubtreePruneRegraft weight=\"%g\">\n",
                           num_tips / 10.0)
     << absl::StreamFormat("      <treeModel idref=\"treeModel\"/>\n")
     << absl::StreamFormat("    </fixedHeightSubtreePruneRegraft>\n");

  if (typeid(pop_model) == typeid(Exp_pop_model)) {
    os << absl::StreamFormat("    <scaleOperator scaleFactor=\"0.75\" weight=\"3\">\n")
       << absl::StreamFormat("      <parameter idref=\"exponential.popSize\"/>\n")
       << absl::StreamFormat("    </scaleOperator>\n")
       << absl::StreamFormat("    <randomWalkOperator windowSize=\"1.0\" weight=\"3\">\n")
       << absl::StreamFormat("      <parameter idref=\"exponential.growthRate\"/>\n")
       << absl::StreamFormat("    </randomWalkOperator>\n");
    
  } else if (typeid(pop_model) == typeid(Skygrid_pop_model)) {
    os << absl::StreamFormat("    <hamiltonianMonteCarloOperator weight=\"2\" nSteps=\"50\" stepSize=\"0.01\" mode=\"vanilla\" autoOptimize=\"true\" gradientCheckCount=\"0\" gradientCheckTolerance=\"0.1\" preconditioning=\"none\" preconditioningUpdateFrequency=\"100\">\n")
       << absl::StreamFormat("      <compoundGradient idref=\"full.skygrid.gradient\"/>\n")
       << absl::StreamFormat("      <compoundParameter idref=\"skygrid.parameters\"/>\n");
    if (run.skygrid_tau_move_enabled()) {
      os << absl::StreamFormat("      <signTransform start=\"1\" end=\"1\">\n")
         << absl::StreamFormat("        <compoundParameter idref=\"skygrid.parameters\"/>\n")
         << absl::StreamFormat("      </signTransform>\n");
    }
    os << absl::StreamFormat("    </hamiltonianMonteCarloOperator>\n");
    
  } else {
    // Make invalid XML tag on purpose to stop this BEAST X XML file from running (with God knows what population model...)
    std::cerr << "ERROR: Unrecognized population model type " << typeid(pop_model).name() << '\n';
    os << absl::StreamFormat("<ERROR>UNRECOGNIZED POPULATION MODEL TYPE: %s</ERROR>", typeid(pop_model).name());
  }

  if (has_tip_uncertainty) {
    for (const auto& node : index_order_traversal(tree)) {
      if (tree.at(node).is_tip() && tree.at(node).t_min != tree.at(node).t_max) {
        os << absl::StreamFormat("    <uniformOperator weight=\"1\">\n")
           << absl::StreamFormat("      <parameter idref=\"%s\"/>\n",
                                 tree.at(node).name)
           << absl::StreamFormat("    </uniformOperator>\n");
      }
    }
  }
  
  os << absl::StreamFormat("  </operators>\n")
     << absl::StreamFormat("\n")
     << absl::StreamFormat("\n");

  // MCMC
  os << absl::StreamFormat("  <!-- Define MCMC                                                             -->\n")
     << absl::StreamFormat("  <mcmc id=\"mcmc\" chainLength=\"%d\" autoOptimize=\"true\" operatorAnalysis=\"output.ops\">\n",
                           chain_length)
     << absl::StreamFormat("    <joint id=\"joint\">\n")
     << absl::StreamFormat("      <prior id=\"prior\">\n");
  
  auto log_kappa_mu = 1.0;  // these defaults come from BEAUti2
  auto log_kappa_sigma = 1.25;  // these defaults come from BEAUti2
  os << absl::StreamFormat("        <logNormalPrior mu=\"%g\" sigma=\"%g\" offset=\"0.0\">\n",
                           log_kappa_mu, log_kappa_sigma)
     << absl::StreamFormat("          <parameter idref=\"kappa\"/>\n")
     << absl::StreamFormat("        </logNormalPrior>\n");
  
  os << absl::StreamFormat("        <dirichletPrior alpha=\"1.0\" sumsTo=\"1.0\">\n")
     << absl::StreamFormat("          <parameter idref=\"frequencies\"/>\n")
     << absl::StreamFormat("        </dirichletPrior>\n");

  if (run.alpha_move_enabled()) {
    auto alpha_mean = 1.0;
    os << absl::StreamFormat("        <exponentialPrior mean=\"%g\" offset=\"0.0\">\n",
                             alpha_mean)
       << absl::StreamFormat("          <parameter idref=\"alpha\"/>\n")
       << absl::StreamFormat("        </exponentialPrior>\n");
  }
  
  // No CTMC prior in Delphy
  
  if (typeid(pop_model) == typeid(Exp_pop_model)) {
    
    if (run.final_pop_size_move_enabled()) {
      os << absl::StreamFormat("        <oneOnXPrior>\n")
         << absl::StreamFormat("          <parameter idref=\"exponential.popSize\"/>\n")
         << absl::StreamFormat("        </oneOnXPrior>\n");
    }
    
    if (run.pop_growth_rate_move_enabled()) {
      auto growth_rate_mu = 0.001;  // per year - these defaults come from BEAUti2
      auto growth_rate_scale = 30.701135;  // per year - these defaults come from BEAUti2
      
      os << absl::StreamFormat("        <laplacePrior mean=\"%f\" scale=\"%f\">\n",
                               growth_rate_mu, growth_rate_scale)
         << absl::StreamFormat("          <parameter idref=\"exponential.growthRate\"/>\n")
         << absl::StreamFormat("        </laplacePrior>\n");
    }

    os << absl::StreamFormat("        <coalescentLikelihood idref=\"coalescent\"/>\n")
       << absl::StreamFormat("\n")
       << absl::StreamFormat("\n");
    
  } else if (typeid(pop_model) == typeid(Skygrid_pop_model)) {
    os << absl::StreamFormat("        <gammaPrior idref=\"skygrid.precision.prior\"/>\n")
       << absl::StreamFormat("        <gmrfSkyGridLikelihood idref=\"skygrid\"/>\n");
    
  } else {
    // Make invalid XML tag on purpose to stop this BEAST X XML file from running (with God knows what population model...)
    std::cerr << "ERROR: Unrecognized population model type " << typeid(pop_model).name() << '\n';
    os << absl::StreamFormat("<ERROR>UNRECOGNIZED POPULATION MODEL TYPE: %s</ERROR>", typeid(pop_model).name());
  }

  os << absl::StreamFormat("        <strictClockBranchRates idref=\"branchRates\"/>\n")
     << absl::StreamFormat("      </prior>\n")
     << absl::StreamFormat("      <likelihood id=\"likelihood\">\n")
     << absl::StreamFormat("        <treeDataLikelihood idref=\"treeLikelihood\"/>\n")
     << absl::StreamFormat("      </likelihood>\n")
     << absl::StreamFormat("    </joint>\n")
     << absl::StreamFormat("    <operators idref=\"operators\"/>\n")
     << absl::StreamFormat("\n");

  // Screen logger
  os << absl::StreamFormat("    <!-- write log to screen                                                     -->\n")
     << absl::StreamFormat("    <log id=\"screenLog\" logEvery=\"1000\">\n")
     << absl::StreamFormat("      <column label=\"Joint\" dp=\"4\" width=\"12\">\n")
     << absl::StreamFormat("        <joint idref=\"joint\"/>\n")
     << absl::StreamFormat("      </column>\n")
     << absl::StreamFormat("      <column label=\"Prior\" dp=\"4\" width=\"12\">\n")
     << absl::StreamFormat("        <prior idref=\"prior\"/>\n")
     << absl::StreamFormat("      </column>\n")
     << absl::StreamFormat("      <column label=\"Likelihood\" dp=\"4\" width=\"12\">\n")
     << absl::StreamFormat("        <likelihood idref=\"likelihood\"/>\n")
     << absl::StreamFormat("      </column>\n")
     << absl::StreamFormat("      <column label=\"age(root)\" sf=\"6\" width=\"12\">\n")
     << absl::StreamFormat("        <tmrcaStatistic idref=\"age(root)\"/>\n")
     << absl::StreamFormat("      </column>\n");

  if (run.mu_move_enabled()) {
    os << absl::StreamFormat("      <column label=\"clock.rate\" sf=\"6\" width=\"12\">\n")
       << absl::StreamFormat("        <parameter idref=\"clock.rate\"/>\n")
       << absl::StreamFormat("      </column>\n");
  }
  
  os << absl::StreamFormat("    </log>\n")
     << absl::StreamFormat("\n");

  // Trace logger
  os << absl::StreamFormat("    <!-- write log to file                                                       -->\n")
     << absl::StreamFormat("    <log id=\"fileLog\" logEvery=\"%d\" fileName=\"output.log\" overwrite=\"false\">\n",
                           log_every)
     << absl::StreamFormat("      <joint idref=\"joint\"/>\n")
     << absl::StreamFormat("      <prior idref=\"prior\"/>\n")
     << absl::StreamFormat("      <likelihood idref=\"likelihood\"/>\n")
     << absl::StreamFormat("      <treeHeightStatistic idref=\"rootHeight\"/>\n")
     << absl::StreamFormat("      <tmrcaStatistic idref=\"age(root)\"/>\n")
     << absl::StreamFormat("      <treeLengthStatistic idref=\"treeLength\"/>\n");
  
  if (typeid(pop_model) == typeid(Exp_pop_model)) {
    if (run.final_pop_size_move_enabled()) {
      os << absl::StreamFormat("      <parameter idref=\"exponential.popSize\"/>\n");
    }
    if (run.pop_growth_rate_move_enabled()) {
      os << absl::StreamFormat("      <parameter idref=\"exponential.growthRate\"/>\n");
    }
    
  } else if (typeid(pop_model) == typeid(Skygrid_pop_model)) {
    os << absl::StreamFormat("      <parameter idref=\"skygrid.precision\"/>\n")
       << absl::StreamFormat("      <parameter idref=\"skygrid.logPopSize\"/>\n")
       << absl::StreamFormat("      <parameter idref=\"skygrid.cutOff\"/>\n");
    
  } else {
    // Make invalid XML tag on purpose to stop this BEAST X XML file from running (with God knows what population model...)
    std::cerr << "ERROR: Unrecognized population model type " << typeid(pop_model).name() << '\n';
    os << absl::StreamFormat("<ERROR>UNRECOGNIZED POPULATION MODEL TYPE: %s</ERROR>", typeid(pop_model).name());
  }
  
  os << absl::StreamFormat("      <parameter idref=\"kappa\"/>\n")
     << absl::StreamFormat("      <parameter idref=\"frequencies\"/>\n");

  if (run.alpha_move_enabled()) {
    os << absl::StreamFormat("      <parameter idref=\"alpha\"/>\n");
  }
  if (run.mu_move_enabled()) {
    os << absl::StreamFormat("      <parameter idref=\"clock.rate\"/>\n");
  }
  os << absl::StreamFormat("      <rateStatistic idref=\"meanRate\"/>\n");

  if (has_tip_uncertainty) {
    os << absl::StreamFormat("\n");
    os << absl::StreamFormat("      <!-- START Tip date sampling                                                 -->\n");
    for (const auto& node : index_order_traversal(tree)) {
      os << absl::StreamFormat("      <parameter idref=\"age(%s)\"/>\n",
                               tree.at(node).name);
    }
    os << absl::StreamFormat("\n")
       << absl::StreamFormat("      <!-- END Tip date sampling                                                   -->\n")
       << absl::StreamFormat("\n");
  }
  
  os << absl::StreamFormat("      <treeDataLikelihood idref=\"treeLikelihood\"/>\n")
     << absl::StreamFormat("      <strictClockBranchRates idref=\"branchRates\"/>\n");

  if (typeid(pop_model) == typeid(Exp_pop_model)) {
    os << absl::StreamFormat("      <coalescentLikelihood idref=\"coalescent\"/>\n");
    
  } else if (typeid(pop_model) == typeid(Skygrid_pop_model)) {
    os << absl::StreamFormat("      <gmrfSkyGridLikelihood idref=\"skygrid\"/>\n");
    
  } else {
    // Make invalid XML tag on purpose to stop this BEAST X XML file from running (with God knows what population model...)
    std::cerr << "ERROR: Unrecognized population model type " << typeid(pop_model).name() << '\n';
    os << absl::StreamFormat("<ERROR>UNRECOGNIZED POPULATION MODEL TYPE: %s</ERROR>", typeid(pop_model).name());
  }
  
  os << absl::StreamFormat("    </log>\n")
     << absl::StreamFormat("\n");

  // Tree logger
  os << absl::StreamFormat("    <!-- write tree log to file                                                  -->\n")
     << absl::StreamFormat("    <logTree id=\"treeFileLog\" logEvery=\"%d\" nexusFormat=\"true\" fileName=\"output.trees\" sortTranslationTable=\"true\">\n",
                           tree_every)
     << absl::StreamFormat("      <treeModel idref=\"treeModel\"/>\n")
     << absl::StreamFormat("      <trait name=\"rate\" tag=\"rate\">\n")
     << absl::StreamFormat("        <strictClockBranchRates idref=\"branchRates\"/>\n")
     << absl::StreamFormat("      </trait>\n")
     << absl::StreamFormat("      <joint idref=\"joint\"/>\n")
     << absl::StreamFormat("    </logTree>\n")
     << absl::StreamFormat("\n")
     << absl::StreamFormat("    <!-- write state of Markov chain to checkpoint file                          -->\n")
     << absl::StreamFormat("    <logCheckpoint id=\"checkpointFileLog\" checkpointEvery=\"%d\" checkpointFinal=\"%d\" fileName=\"output.chkpt\" overwrite=\"false\"/>\n",
                           chain_length / 100, chain_length)
     << absl::StreamFormat("  </mcmc>\n")
     << absl::StreamFormat("\n")
     << absl::StreamFormat("  <report>\n")
     << absl::StreamFormat("    <property name=\"timer\">\n")
     << absl::StreamFormat("      <mcmc idref=\"mcmc\"/>\n")
     << absl::StreamFormat("    </property>\n")
     << absl::StreamFormat("  </report>\n")
     << absl::StreamFormat("\n");

  // That's it
  os << absl::StreamFormat("</beast>\n");
}

auto export_beast_input(
    const Run& run,
    std::string_view beast_version,
    std::ostream& os,
    int64_t chain_length,
    int64_t log_every,
    int64_t tree_every)
    -> void {

  if (beast_version == "2.6.2") {
    export_beast_2_6_2_input(run, os, chain_length, log_every, tree_every);
  } else if (beast_version == "X-10.5.0") {
    export_beast_X_10_5_0_input(run, os, chain_length, log_every, tree_every);
  } else {
    std::cerr << "ERROR: BEAST input XML generation not currently supported for BEAST version '" << beast_version << "'\n";
    os << "ERROR: BEAST input XML generation not currently supported for BEAST version '" << beast_version << "'\n";
  }
}

}  // namespace delphy
