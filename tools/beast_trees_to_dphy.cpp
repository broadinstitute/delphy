#include "absl/log/initialize.h"
#include "absl/log/check.h"
#include "absl/random/random.h"

#include "cxxopts.hpp"

#include "beasty_input.h"
#include "version.h"
#include "io.h"
#include "api.h"
#include "cmdline.h"
#include "phylo_tree_calc.h"

namespace delphy {

// Shuffle the nodes in a Phylo_tree to have all the tips come first, in the order they appear on the FASTA file
auto renumber_trees(
    std::map<int64_t, Phylo_tree>& trees,
    const std::vector<Fasta_entry>& fasta)
    -> void {
  
  auto fasta_id_to_idx = std::map<std::string, int>{};
  for (auto i = 0; i != std::ssize(fasta); ++i) {
    const auto& entry = fasta[i];
    auto [it, inserted] = fasta_id_to_idx.try_emplace(entry.id, i);
    if (not inserted) {
      std::cerr << "ERROR: Duplicate FASTA ID: " << entry.id << "\n";
      std::exit(EXIT_FAILURE);
    }
  }

  for (auto& [state, tree] : trees) {
    assert_phylo_tree_integrity(tree);
    
    auto old_id_to_new_id = Node_vector<Node_index>(std::ssize(tree), k_no_node);

    auto num_tips = (std::ssize(tree) + 1) / 2;
    auto next_inner_node_id = num_tips;  // Tips numbered [0, num_tips)
    for (const auto& node : index_order_traversal(tree)) {
      if (tree.at(node).is_tip()) {
        auto new_id = fasta_id_to_idx[tree.at(node).name];
        CHECK_GE(new_id, 0);
        CHECK_LT(new_id, num_tips);
        old_id_to_new_id[node] = new_id;
      } else {
        old_id_to_new_id[node] = next_inner_node_id;
        ++next_inner_node_id;
      }
    }

    // Reorder nodes
    auto new_nodes = tree.nodes; // Copy, with the right allocator
    for (auto old_id = Node_index{0}; old_id != std::ssize(tree); ++old_id) {
      auto new_id = old_id_to_new_id[old_id];
      new_nodes[new_id] = tree.at(old_id);
    }
    swap(tree.nodes, new_nodes);

    // Patch up links
    for (const auto& node : index_order_traversal(tree)) {
      auto old_parent = tree.at(node).parent;
      if (old_parent != k_no_node) {
        tree.at(node).parent = old_id_to_new_id[old_parent];
      }
      if (tree.at(node).is_inner_node()) {
        auto old_left_child = tree.at(node).left_child();
        auto old_right_child = tree.at(node).right_child();
        tree.at(node).children = {old_id_to_new_id[old_left_child], old_id_to_new_id[old_right_child]};
      }
    }
    tree.root = old_id_to_new_id[tree.root];

    assert_phylo_tree_integrity(tree);
  }
}

auto redate_trees(
    std::map<int64_t, Phylo_tree>& trees,
    const std::string& pinned_ref)
    -> void {

  auto maybe_date = extract_date_from_sequence_id(pinned_ref);
  if (not maybe_date.has_value()) {
    std::cerr << absl::StreamFormat("ERROR: Could not determine date of pinned reference sequence '%s'\n", pinned_ref);
    std::exit(EXIT_FAILURE);
  }
  auto target_date = maybe_date.value();
  
  for (auto& [state, tree] : trees) {
    auto t_offset = 0.0;

    // Find pinned_ref in tree and calculate t_offset
    for (const auto& node : index_order_traversal(tree)) {
      if (tree.at(node).name == pinned_ref) {
        t_offset = target_date - tree.at(node).t;
        break;
      }
    }

    // Adjust all times
    for (const auto& node : index_order_traversal(tree)) {
      tree.at(node).t += t_offset;
      for (auto& m : tree.at(node).mutations) {
        m.t += t_offset;
      }
    }
  }
}

// Alias C = G and T = A at APOBEC-context site, as in the BEAST runs of O'Toole et al (2023)
auto apobec_alias(Real_seq_letter s) -> Real_seq_letter {
  switch (s) {
    case Real_seq_letter::G: return Real_seq_letter::C;
    case Real_seq_letter::A: return Real_seq_letter::T;
    default: return s;
  }
}

// (the tests below mimics the BEAST run's behaviour, where C and T are aliased, as are G and A)
auto is_apobec_mut_allowed(Real_seq_letter from, Real_seq_letter to) {
  auto aliased_from = apobec_alias(from);
  auto aliased_to = apobec_alias(to);
  return (aliased_from == aliased_to) ||
      ((aliased_from == Real_seq_letter::C) && (aliased_to == Real_seq_letter::T));
}

auto ancestralize(
    std::map<int64_t, Phylo_tree>& trees,
    const std::vector<Fasta_entry>& fasta,
    absl::BitGenRef bitgen,
    const absl::flat_hash_set<Site_index>& apobec_sites,
    bool forbid_bad_apobec_mutations)
    -> void {

  auto fasta_id_to_entry = std::map<std::string, const Fasta_entry*>{};
  for (const auto& entry : fasta) {
    auto [it, inserted] = fasta_id_to_entry.try_emplace(entry.id, &entry);
    if (not inserted) {
      std::cerr << "ERROR: Duplicate FASTA ID: " << entry.id << "\n";
      std::exit(EXIT_FAILURE);
    }
  }

  auto L = std::ssize(fasta_id_to_entry.begin()->second->sequence);
  auto always_missing = std::vector<char>(L, 1);  // Avoid vector<bool>!
  
  for (auto& [state, tree] : trees) {
    // Init ref sequence with correct size
    tree.ref_sequence = Real_sequence(L, Real_seq_letter::A);
  
    // Try to match tips to FASTA IDs
    for (const auto& node : index_order_traversal(tree)) {
      if (tree.at(node).is_tip()) {
        if (auto it = fasta_id_to_entry.find(tree.at(node).name); it == fasta_id_to_entry.end()) {
          std::cerr << "Could not find tip with name \"" << tree.at(node).name << "\" in FASTA file\n";
          std::exit(EXIT_FAILURE);
        } else {
          const auto& seq = it->second->sequence;
          /*std::cerr << "Sequence for " << tree.at(node).name << ":\n";
          for (auto l = Site_index{0}; l != tree.num_sites(); ++l) {
            std::cerr << to_char(seq[l]);
          }
          std::cerr << "\n";*/
          for (auto l = Site_index{0}; l != tree.num_sites(); ++l) {
            if (not is_ambiguous(seq[l])) {
              always_missing[l] = 0;
            }
          }
        }
      }
    }
  }


  // Main Fitch-Sankoff algorithm: FS sets are encoded using a Sequence of Seq_letters
  // NOTE: This isn't a very efficient implementation, but it's "efficient enough"
  for (auto& [state, tree] : trees) {
    std::cerr << "Running Fitch-Sankoff over tree " << state << "\n";
    
    auto fs_sets = Node_vector<Sequence>(std::ssize(tree), Sequence{});

    // Upwards pass
    for (const auto& node : post_order_traversal(tree)) {
      if (tree.at(node).is_tip()) {
        // FS sets come directly from sequence
        fs_sets.at(node) = fasta_id_to_entry[tree.at(node).name]->sequence;
      } else {
        const auto& left_child_fs = fs_sets[tree.at(node).left_child()];
        const auto& right_child_fs = fs_sets[tree.at(node).right_child()];
        auto& this_fs = fs_sets.at(node) = Sequence(tree.num_sites(), Seq_letters::N);
        for (auto l = Site_index{0}; l != tree.num_sites(); ++l) {
          if ((left_child_fs[l] & right_child_fs[l]) != 0) {
            // Any choice in the FS set requires *no* mutations to the child branches
            this_fs[l] = left_child_fs[l] & right_child_fs[l];
          } else {
            // Any choice in the FS set requires *1* mutation in one of the child branches
            this_fs[l] = left_child_fs[l] | right_child_fs[l];

            // But what if one of those mutations is a forbidden APOBEC mutation?
            if (forbid_bad_apobec_mutations && apobec_sites.contains(l)) {
              // Do this the really pedestrian but obviously correct way
              this_fs[l] = Seq_letters::none;
              auto this_fs_two_muts = Seq_letters::none;
              auto this_fs_one_mut = Seq_letters::none;
              auto this_fs_no_muts = Seq_letters::none;
              
              for (const auto& this_state : k_all_real_seq_letters) {
                for (const auto& left_state : k_all_real_seq_letters) {
                  if ((left_child_fs[l] & to_seq_letter(left_state)) == 0) { continue; }  // Can't be
                  for (const auto& right_state : k_all_real_seq_letters) {
                    if ((right_child_fs[l] & to_seq_letter(right_state)) == 0) { continue; }  // Can't be

                    // At this point, left_state is possible for left child, right_state is possible for right child

                    // Check that the implied this-to-left and this-to-right mutations are allowed
                    auto allowed_left = is_apobec_mut_allowed(this_state, left_state);
                    auto allowed_right = is_apobec_mut_allowed(this_state, right_state);

                    if (allowed_left && allowed_right) {
                      // Ok!  this_state can be used at this node
                      auto num_muts = (this_state == left_state ? 0 : 1) + (this_state == right_state ? 0 : 1);
                      switch (num_muts) {
                        case 0: this_fs_no_muts |= to_seq_letter(this_state); break;
                        case 1: this_fs_one_mut |= to_seq_letter(this_state); break;
                        case 2: this_fs_two_muts |= to_seq_letter(this_state); break;
                        default: CHECK(false) << "???";
                      }
                    }
                  }
                }
              }

              // Pick the most parsimonious assignment
              CHECK_EQ(this_fs_no_muts, Seq_letters::none);  // Otherwise, we'd have landed in the other branch
              if (this_fs_one_mut != Seq_letters::none) {
                this_fs[l] = this_fs_one_mut;
              } else {
                this_fs[l] = this_fs_two_muts;
              }
            }
          }
        }
      }
    }

    // Downwards pass - FS sets resolved arbitrarily for parents before children
    // TODO: resolve randomly instead of naively
    auto resolve_naive = [](Seq_letter s) -> Seq_letter {
      if ((s & Seq_letters::A) != 0) { return Seq_letters::A; }
      if ((s & Seq_letters::C) != 0) { return Seq_letters::C; }
      if ((s & Seq_letters::G) != 0) { return Seq_letters::G; }
      if ((s & Seq_letters::T) != 0) { return Seq_letters::T; }
      std::cerr << "ERROR: Can't resolve Seq_letter " << to_char(s) << "!\n";
      std::exit(EXIT_FAILURE);
    };
    auto& root_fs = fs_sets.at(tree.root);
    for (auto l = Site_index{0}; l != tree.num_sites(); ++l) {
      if (is_ambiguous(root_fs[l])) {
        
        // Ambiguous states at the root owing to not knowing where to place a mutation between the
        // children when these have definite states is a nothing-burger
        auto resolved = resolve_naive(root_fs[l]);
        auto simple_root_resolution = false;
        if (tree.at_root().is_inner_node()) {
          auto left_child = tree.at_root().left_child();
          auto right_child = tree.at_root().right_child();
          auto left_child_fs = fs_sets.at(left_child);
          auto right_child_fs = fs_sets.at(right_child);
          if (not is_ambiguous(left_child_fs[l]) && not is_ambiguous(right_child_fs[l])) {
            // Resolve in favor of closest child (probably what's being used to root the tree anyway)
            auto closer_child = tree.at(left_child).t < tree.at(right_child).t ? left_child : right_child;
            resolved = fs_sets.at(closer_child)[l];
            simple_root_resolution = true;
          }
        }
        
        // Only warn if the site is not missing in all sequences and it's not a simple root resolution above
        if (not always_missing[l] && not simple_root_resolution) {
          std::cerr << "State at site " << l << " of root is ambiguous: " << to_char(root_fs[l])
                    << " and is present in some sequences; resolving arbitrarily to " << to_char(resolved) << "\n";
        }
        root_fs[l] = resolved;
      }
    }
      
    for (const auto& node : pre_order_traversal(tree)) {
      if (node == tree.root) continue;  // Handled above
      auto& this_fs = fs_sets.at(node);
      const auto& parent_fs = fs_sets.at(tree.at(node).parent);
      
      for (auto l = Site_index{0}; l != tree.num_sites(); ++l) {
        if (is_ambiguous(this_fs[l])) {
          CHECK(not is_ambiguous(parent_fs[l]));
          if ((parent_fs[l] & this_fs[l]) != 0) {
            // Resolve by choosing parent state (and avoiding needless mutation)
            this_fs[l] = parent_fs[l];
          } else {
            auto resolved = resolve_naive(this_fs[l]);
            // If bad APOBEC mutations are forbidden, some of the naive choices might actually be illegal!
            if (forbid_bad_apobec_mutations && apobec_sites.contains(l)) {
              auto parent_state = to_real_seq_letter(parent_fs[l]);
              auto valid_this_fs = Seq_letters::none;
              for (auto this_state : k_all_real_seq_letters) {
                if ((this_fs[l] & to_seq_letter(this_state)) != 0 && is_apobec_mut_allowed(parent_state, this_state)) {
                  valid_this_fs |= to_seq_letter(this_state);
                }
              }
              if (valid_this_fs != this_fs[l]) {
                std::cerr << "FS set at site " << l << " of node " << node << " narrowed from " << to_char(this_fs[l])
                          << " to " << to_char(valid_this_fs) << " to forbid bad APOBEC mutations\n";
                this_fs[l] = valid_this_fs;
                resolved = resolve_naive(valid_this_fs);
              }
            }
            if (is_ambiguous(this_fs[l])) {  // Ambiguity may have been resolved by APOBEC rules
              std::cerr << "State at site " << l << " of node " << node << " is ambiguous: " << to_char(this_fs[l])
                        << "; resolving arbitrarily to " << to_char(resolved) << "\n";
            }
            this_fs[l] = resolved;
          }
        }
      }
    }

    // Now turn FS sets into real mutations
    for (auto l = Site_index{0}; l != tree.num_sites(); ++l) {
      CHECK(not is_ambiguous(root_fs[l]));
      tree.ref_sequence[l] = to_real_seq_letter(root_fs[l]);
    }
    tree.at_root().mutations.clear();
    
    for (const auto& node : pre_order_traversal(tree)) {
      if (node == tree.root) continue;  // Handled above
      
      auto& this_fs = fs_sets.at(node);
      const auto& parent_fs = fs_sets.at(tree.at(node).parent);
      tree.at(node).mutations.clear();
      
      for (auto l = Site_index{0}; l != tree.num_sites(); ++l) {
        if (this_fs[l] != parent_fs[l]) {
          auto m = Mutation{
            to_real_seq_letter(parent_fs[l]), l, to_real_seq_letter(this_fs[l]),
            absl::Uniform(bitgen, tree.at_parent_of(node).t, tree.at(node).t)};
          tree.at(node).mutations.push_back(m);

          if (apobec_sites.contains(l) && not is_apobec_mut_allowed(m.from, m.to)) {
            std::cerr << "GAHHH!!! Inadmissible APOBEC-site mutation is part of parsimonious reconstruction: " << m << "\n";
          }
        }
      }

      sort_mutations(tree.at(node).mutations);
    }
    
    assert_phylo_tree_integrity(tree);
  }
}

auto find_subtree_roots(
    std::map<int64_t, Phylo_tree>& trees,
    const std::set<std::string>& only_subtree_below_mrca_of_tips)
    -> std::map<int64_t, Node_index> {
  
  auto state_to_mrca_node = std::map<int64_t, Node_index>{};
  
  // Identify MRCA of given tips in every tree and ensure it identifies subtrees with equal numbers of nodes
  auto num_nodes = -1;
  for (const auto& [state, tree] : trees) {
    if (only_subtree_below_mrca_of_tips.empty()) {
      state_to_mrca_node[state] = tree.root;  // Nothing to subset
      continue;
    }
    
    auto tips_here = std::vector<Node_index>{};
    for (const auto& node : index_order_traversal(tree)) {
      if (tree.at(node).is_tip() && only_subtree_below_mrca_of_tips.contains(tree.at(node).name)) {
        tips_here.push_back(node);
      }
    }
    CHECK_EQ(std::ssize(tips_here), std::ssize(only_subtree_below_mrca_of_tips));

    // Run an MRCA fold over tips
    auto mrca_so_far = k_no_node;
    for (const auto& tip : tips_here) {
      if (mrca_so_far == k_no_node) {
        mrca_so_far = tip;
      } else {
        mrca_so_far = find_MRCA_of(tree, mrca_so_far, tip);
      }
    }

    state_to_mrca_node[state] = mrca_so_far;

    auto descendant_counts = count_all_descendants(tree);  // Expensive, but who cares...
    if (num_nodes == -1) {
      num_nodes = 1 + descendant_counts[mrca_so_far];
    } else {
      CHECK_EQ(num_nodes, 1 + descendant_counts[mrca_so_far]) << state;
    }
  }

  return state_to_mrca_node;
}

auto extract_subtrees(
    std::map<int64_t, Phylo_tree>& trees,
    const std::map<int64_t, Node_index>& state_to_mrca_node)
    -> void {

  // Count nodes in subtrees (in find_subtree_roots, we checked that these were the same in all trees)
  const auto& [first_state, first_tree] = *trees.begin();
  auto descendant_counts = count_all_descendants(first_tree);  // Expensive, but who cares...
  auto num_nodes = 1 + descendant_counts[state_to_mrca_node.at(first_state)];

  // Extract each subtree in turn
  CHECK_NE(num_nodes, -1);
  auto num_tips = (num_nodes + 1) / 2;
  for (auto& [state, tree] : trees) {
    auto mrca_node = state_to_mrca_node.at(state);
    if (mrca_node == tree.root) { continue; }  // Nothing to do
    
    auto next_tip_index = 0;                 // Tips: [0, num_tips)
    auto next_inner_node_index = num_tips;   // Inner nodes: [num_tips, num_nodes)
    auto old_id_to_new_id = absl::flat_hash_map<Node_index, Node_index>{};

    auto new_nodes = std::vector<Phylo_node>(num_nodes, Phylo_node{});

    // Extract nodes and assign new IDs
    for (const auto& old_id : index_order_traversal(tree)) {
      // Wildly inefficient way to traverse subtree of mrca_node, but who cares...
      if (descends_from(tree, old_id, mrca_node)) {
        auto new_id = Node_index{};
        if (tree.at(old_id).is_tip()) {
          CHECK_GE(next_tip_index, 0);
          CHECK_LT(next_tip_index, num_tips);
          new_id = next_tip_index;
          ++next_tip_index;
        } else {
          CHECK_GE(next_inner_node_index, num_tips);
          CHECK_LT(next_inner_node_index, num_nodes);
          new_id = next_inner_node_index;
          ++next_inner_node_index;
        }
        old_id_to_new_id[old_id] = new_id;

        new_nodes[new_id] = tree.at(old_id);
      }
    }
    assert_phylo_tree_integrity(tree);

    // Patch up links
    CHECK(old_id_to_new_id.contains(mrca_node));
    auto new_root = old_id_to_new_id[mrca_node];
    for (auto new_id = Node_index{0}; new_id != num_nodes; ++new_id) {
      auto old_parent = new_nodes.at(new_id).parent;
      if (old_parent == k_no_node) {
        new_nodes.at(new_id).parent = k_no_node;
      } else {
        if (old_id_to_new_id.contains(old_parent)) {
          CHECK_NE(new_id, new_root);
          new_nodes.at(new_id).parent = old_id_to_new_id[old_parent];
        } else {
          CHECK_EQ(new_id, new_root);
          new_nodes.at(new_id).parent = k_no_node;
        }
      }
      if (new_nodes.at(new_id).is_inner_node()) {
        auto old_left_child = new_nodes.at(new_id).left_child();
        auto old_right_child = new_nodes.at(new_id).right_child();
        new_nodes.at(new_id).children = {old_id_to_new_id[old_left_child], old_id_to_new_id[old_right_child]};
      }
    }

    // Set up new tree, with reference sequence including all the mutations upstream of mrca_node
    auto new_tree = Phylo_tree{std::move(new_nodes)};
    new_tree.root = new_root;
    new_tree.at_root().mutations.clear();
    new_tree.ref_sequence = view_of_sequence_at(tree, mrca_node).materialize();

    assert_phylo_tree_integrity(new_tree);
    
    tree = std::move(new_tree);
  }
}

auto fix_up_node_names(std::map<int64_t, Phylo_tree>& trees) -> void {

  for (auto& [state, tree] : trees) {
    // Resolve date uncertainties crudely
    // (I thought I needed this for Delphy UI, but I don't)
    // for (const auto& node : index_order_traversal(tree)) {
    //   if (tree.at(node).is_tip()) {
    //     auto last_pipe_pos = tree.at(node).name.find_last_of('|');
    //     auto last_pos = std::ssize(tree.at(node).name);
    //     auto date_chars = last_pos - (last_pipe_pos+1);
    //     switch(date_chars) {
    //       case 4+1+2+1+2: break;
    //       case 4+1+2: tree.at(node).name += "-15"; break;  // Impute to middle of month
    //       case 4: tree.at(node).name += "-07-01"; break;   // Impute to middle of year
    //       default:
    //         std::cerr << "Invalid date string for '" << tree.at(node).name << "'\n";
    //         std::exit(EXIT_FAILURE);
    //     }
    //   }
    // }

    // Change all but the last pipe in tip names to slashes (bad interaction with Delphy UI)
    for (const auto& node : index_order_traversal(tree)) {
      if (tree.at(node).is_tip()) {
        auto& name = tree.at(node).name;
        auto last_pipe_pos = name.find_last_of('|');
        for (auto i = std::string::size_type{0}; i != last_pipe_pos; ++i) {
          if (name[i] == '|') {
            name[i] = '/';
          }
        }
      }
    }
  }
}

auto write_string(std::ostream& os, std::string_view sv) -> void {
  for (const auto& c : sv) {
    os.put(c);
  }
}

auto write_uint32(std::ostream& os, uint32_t num) -> void {
  os
      .put((num >>  0) & 0xff)
      .put((num >>  8) & 0xff)
      .put((num >> 16) & 0xff)
      .put((num >> 24) & 0xff);
}

auto write_uint64(std::ostream& os, uint64_t num) -> void {
  os
      .put((num >>  0) & 0xff)
      .put((num >>  8) & 0xff)
      .put((num >> 16) & 0xff)
      .put((num >> 24) & 0xff)
      .put((num >> 32) & 0xff)
      .put((num >> 40) & 0xff)
      .put((num >> 48) & 0xff)
      .put((num >> 56) & 0xff);
}

auto write_float32(std::ostream& os, float num) -> void {
  write_uint32(os, std::bit_cast<uint32_t>(num));
}

auto get_flatbuffer_size(const flatbuffers::DetachedBuffer& buffer) -> uint32_t {
  return flatbuffers::GetPrefixedSize(buffer.data()) + 4;  // +4 => account for the size prefix!
}

auto write_flatbuffer_size(std::ostream& os, const flatbuffers::DetachedBuffer& buffer) -> void {
  write_uint32(os, get_flatbuffer_size(buffer));
}

auto write_flatbuffer(std::ostream& os, const flatbuffers::DetachedBuffer& buffer) -> void {
  os.write(reinterpret_cast<const char*>(buffer.data()), get_flatbuffer_size(buffer));
}

auto write_delphy_file(
    std::map<int64_t, Phylo_tree>& trees,
    std::string_view out_filename,
    int every,
    std::mt19937& bitgen)
    -> void {

  const auto APP_MAGIC_NUMBER = "DPHY";
  const auto SAVE_FORMAT_VERSION = 3;
  const auto NO_MORE_TREES = 0;
  const auto config = std::string{
    "{"
    "\"confidence\":90,"
    "\"topology\":0,"
    "\"presentation\":0,"
    "\"spacing\":0,"
    "\"colorBy\":0,"
    "\"burnin\":1,"
    "\"metadataPresent\":0,"
    "\"metadataText\":null,"
    "\"metadataFile\":null,"
    "\"metadataDelimiter\":null,"
    "\"selectedMDField\":-1,"
    "\"metadataColors\":{}"
    "}"};
  
  auto num_threads = 1;  // Don't need to really run anything...
  auto thread_pool = ctpl::thread_pool{num_threads};
  auto fake_run = Run{thread_pool, bitgen, trees.begin()->second};

  // Fill in parameters with values that are hopefully obviously fake
  fake_run.set_mu(42e-3 / 365.0);
  fake_run.set_hky_kappa(42.0);
  fake_run.set_hky_pi({0.1, 0.2, 0.3, 0.4});
  fake_run.set_pop_model(Exp_pop_model{0.0, 42.0 / 365.0, 42.0 * 365.0});
  
  auto api_tree_info = phylo_tree_to_api_tree_info(trees.begin()->second);

  // Here we go
  auto os = std::ofstream{std::string{out_filename}, std::ios::binary};
  write_string(os, APP_MAGIC_NUMBER);
  write_uint32(os, SAVE_FORMAT_VERSION);
  write_uint32(os, k_delphy_version_string.length());
  write_string(os, k_delphy_version_string);
  write_uint32(os, k_delphy_build_number);
  write_uint32(os, k_delphy_commit_string.length());
  write_string(os, k_delphy_commit_string);
  write_uint32(os, 1);  // knee_index
  write_uint32(os, every);  // Steps per sample
  write_uint32(os, 0);  // run.is_alpha_move_enabled()
  write_uint32(os, 1);  // run.is_mpox_hack_enabled()
  write_uint32(os, 0);  // run.is_mu_move_enabled()
  write_float32(os, 0.0);  // run.mu() ???
  write_flatbuffer_size(os, api_tree_info);
  write_flatbuffer(os, api_tree_info);
  
  for (const auto& [state, tree] : trees) {
    fake_run.set_step(state);
    auto api_tree = phylo_tree_to_api_tree(tree);
    auto api_params = run_to_api_params(fake_run);
    write_flatbuffer_size(os, api_tree);
    write_flatbuffer_size(os, api_params);
    write_flatbuffer(os, api_tree);
    write_flatbuffer(os, api_params);
  }

  auto tree_end_position = os.tellp();

  write_uint32(os, NO_MORE_TREES);
  write_uint32(os, config.length());
  write_string(os, config);
  write_uint64(os, tree_end_position);
}

}  // namespace delphy

auto main(int argc, char** argv) -> int {
  using namespace delphy;

  absl::InitializeLog();

  cxxopts::Options options("beast_trees_to_dphy", "Infer plausible mutations on BEAST trees (parsimony model) and write out an equivalent .dphy file (we're targetting exploring the trees with the Delphy UI, not the traces, lineage/mutation prevalences, etc.)");

  options.add_options("Generic options")
      ("version", "print version string")
      ("h,help", "print usage")
      ;

  options.add_options("Input/output options")
      ("fasta", "input FASTA file", cxxopts::value<std::string>())
      ("trees", "input BEAST trees file", cxxopts::value<std::string>())
      ("every", "Only process a BEAST tree if the state number is a multiple of this number",
       cxxopts::value<int64_t>()->default_value("1"))
      ("pinned-ref", "FASTA id (including date string) of tip with fixed time (to convert relative to absolute times)",
       cxxopts::value<std::string>())
      ("apobec-sites", "Optional list of 0-based sites that supposedly only undergo APOBEC3-like mutations",
       cxxopts::value<std::string>())
      ("only-subtree-below-mrca-of",
       "Optional comma-separated list of tips whose MRCA defines the root of a subtree to export",
       cxxopts::value<std::string>())
      ("forbid-bad-apobec-mutations", "When reconstructing mutations, forbid anything but C/G>T/A at an APOBEC site",
       cxxopts::value<bool>()->default_value("false"))
      ("out", "Filename for Delphy .dphy output file",
       cxxopts::value<std::string>())
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

    auto fasta_filename = opts["fasta"].as<std::string>();
    auto trees_filename = opts["trees"].as<std::string>();
    auto every = opts["every"].as<int64_t>();
    auto pinned_ref = opts["pinned-ref"].as<std::string>();
    auto out_filename = opts["out"].as<std::string>();

    auto apobec_sites = absl::flat_hash_set<Site_index>{};
    if (opts.count("apobec-sites")) {
      auto apobec_sites_filename = opts["apobec-sites"].as<std::string>();
      auto apobec_sites_is = std::ifstream{apobec_sites_filename};
      auto l = Site_index{0};
      while (apobec_sites_is >> l) {
        apobec_sites.insert(l);
      }
    }
    auto forbid_bad_apobec_mutations = opts["forbid-bad-apobec-mutations"].as<bool>();
    if (forbid_bad_apobec_mutations && apobec_sites.empty()) {
      std::cerr << "Something's wrong: --forbid-bad-apobec-mutations make no sense if no --apobec-sites specified\n";
      std::exit(EXIT_FAILURE);
    }

    auto only_subtree_below_mrca_of_tips = std::set<std::string>{};
    if (opts.count("only-subtree-below-mrca-of")) {
      auto tip_list = opts["only-subtree-below-mrca-of"].as<std::string>();
      auto iss = std::istringstream{tip_list};
      auto tip = std::string{};
      while (std::getline(iss, tip, ',')) {
        only_subtree_below_mrca_of_tips.insert(tip);
      }
    }
    
    auto fasta_is = std::ifstream{fasta_filename};
    if (not fasta_is) {
      std::cerr << "ERROR: Could not read input FASTA file " << fasta_filename << "\n";
      std::exit(EXIT_FAILURE);
    }
    std::cerr << "Reading fasta file " << fasta_filename << "\n";
    auto fasta = read_fasta(fasta_is);
    fasta_is.close();

    std::cerr << "Read " << std::ssize(fasta) << " entries\n";

    auto trees_is = std::ifstream{trees_filename};
    if (not trees_is) {
      std::cerr << "ERROR: Could not read BEAST trees file " << trees_filename << "\n";
      std::exit(EXIT_FAILURE);
    }
    std::cerr << "Reading BEAST trees file " << trees_filename << "\n";
    auto trees = read_beasty_trees(trees_is, 0, every);
    trees_is.close();

    auto bitgen = std::mt19937{std::random_device{}()};
    std::cerr << "Renumbering tree nodes\n";
    renumber_trees(trees, fasta);
    
    std::cerr << "Finding subtree roots\n";
    auto state_to_mrca_node = find_subtree_roots(trees, only_subtree_below_mrca_of_tips);
    std::cerr << "Reconstructing ancestral states & plausible mutations\n";
    std::cerr << "    WARNING: This uses a parsimony model!  BEWARE!\n";
    ancestralize(trees, fasta, bitgen, apobec_sites, forbid_bad_apobec_mutations);
    std::cerr << "Redating trees to absolute time\n";
    redate_trees(trees, pinned_ref);
    std::cerr << "Extracting subtrees\n";
    extract_subtrees(trees, state_to_mrca_node);
    std::cerr << "Fixing up node names\n";
    fix_up_node_names(trees);

    std::cerr << "Writing out Delphy file " << out_filename << "\n";

    // Temp
    write_delphy_file(trees, out_filename, every, bitgen);
    
  } catch (cxxopts::exceptions::exception& x) {
    std::cerr << "ERROR: " << x.what() << "\n" << options.help() << "\n";
    std::exit(EXIT_FAILURE);
  }
  
  return 0;
}
