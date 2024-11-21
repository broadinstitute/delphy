#include "io.h"

#include <iostream>

#include <boost/algorithm/string/replace.hpp>

#include "version.h"
#include "cmdline.h"
#include "phylo_tree_calc.h"

namespace delphy {

auto read_fasta(
    std::istream& is,
    const std::function<void(int, std::size_t)>& progress_hook)
    -> std::vector<Fasta_entry> {
  
  auto in_seq = false;
  auto cur_id = std::string{};
  auto cur_comments = std::string{};
  auto cur_seq = Sequence{};

  auto result = std::vector<Fasta_entry>{};
  auto bytes_so_far = std::size_t{0};

  for (auto line = std::string{}; getline(is, line);) {
    bytes_so_far += std::ssize(line) + 1;  // `+ 1` from end-of-line character
    
    if (line.empty()) {
      continue;  // ignore empty lines
    }

    if (line[0] == ';') {
      continue;  // comment line, ignore
    }

    if (line[0] == '>') {
      if (in_seq) {
        result.push_back(Fasta_entry{std::move(cur_id), std::move(cur_comments), std::move(cur_seq)});
        progress_hook(std::ssize(result), bytes_so_far);
        cur_seq.clear();
      }

      auto id_begin = std::find_if_not(line.begin() + 1, line.end(), [](auto c) { return std::isspace(c); });
      auto id_end = std::find_if(id_begin, line.end(), [](auto c) { return std::isspace(c); });
      cur_id.assign(id_begin, id_end);

      auto comments_begin = (id_end != line.end()) ? (id_end + 1) : id_end;
      auto comments_end = line.end();
      cur_comments.assign(comments_begin, comments_end);

      in_seq = true;
      continue;
    }

    if (!in_seq) {
      // Sequence data preceding an initial identifier line (e.g., '>SeqName') is read
      // into a sequence named ""
      in_seq = true;
    }

    for (auto c : line) {
      if (auto s = to_seq_letter(c); s != Seq_letters::none) {
        cur_seq.push_back(s);
      }
    }
  }
  if (in_seq) {
    result.push_back(Fasta_entry{std::move(cur_id), std::move(cur_comments), std::move(cur_seq)});
    progress_hook(std::ssize(result), bytes_so_far);
  }

  return result;
}

auto read_maple(std::istream& is) -> Maple_file {
  auto result = Maple_file{};

  auto line = std::string{};
  if (not getline(is, line)) {
    throw std::runtime_error("ERROR: Unexpected EOF while reading MAPLE file reference sequence");
  }
  if (line.empty()) {
    throw std::runtime_error("ERROR: Unexpected empty reference id line in MAPLE file");
  }
  if (line[0] != '>') {
    throw std::runtime_error(absl::StrFormat(
        "ERROR: Expected reference sequence id line to start with '>', but saw this instead: %s",
        line));
  }
  // Ignore reference id

  // Now read reference sequence
  if (not getline(is, line)) {
    throw std::runtime_error("ERROR: Unexpected EOF while reading MAPLE file reference sequence");
  }
  do {
    if (line.empty()) {
      continue;
    }
    if (line[0] == '>') {
      break;
    }
    
    for (auto c : line) {
      if (not std::isspace(c)) {
        auto s = char_to_real_seq_letter(c);
        result.ref_sequence.push_back(s);
      }
    }
  } while (getline(is, line));

  // Now read each sequence in turn
  while (is) {
    if (line.empty() || line[0] != '>') {
      throw std::runtime_error(absl::StrFormat(
          "Expected sequence id line to start with '>', but saw this instead: %s",
          line));
    }

    auto tip_desc = Tip_desc{};
    
    auto id_begin = std::find_if_not(line.begin() + 1, line.end(), [](auto c) { return std::isspace(c); });
    auto id_end = std::find_if(id_begin, line.end(), [](auto c) { return std::isspace(c); });
    tip_desc.name.assign(id_begin, id_end);

    auto maybe_t = extract_date_from_sequence_id(tip_desc.name);
    tip_desc.t = maybe_t.has_value() ? maybe_t.value() : 0.0;

    while(getline(is, line)) {
      if (line.empty()) {
        continue;  // ignore empty lines
      }
      if (line[0] == '>') {
        break;  // on to next tip
      }

      auto ss = std::istringstream{line};
      auto c = char{};
      ss >> c;
      auto s = to_seq_letter(line[0]);
      if (is_ambiguous(s)) {
        // Gap or Missing Data
        auto start_site = Site_index{};  ss >> start_site;  --start_site;  // 0-based indexing
        auto gap_len = Site_index{};  ss >> gap_len;
        auto end_site = start_site + gap_len;
        tip_desc.missations.intervals.insert(Site_interval{start_site, end_site});
      } else {
        // Mutation w.r.t. reference
        auto l = Site_index{};  ss >> l;  --l;  // 0-based indexing
        auto from = result.ref_sequence.at(l);
        tip_desc.seq_deltas.push_back({l, from, to_real_seq_letter(s)});
      }
    }

    // Only include sequence if we can date it!
    if (maybe_t.has_value()) {
      result.tip_descs.push_back(std::move(tip_desc));
    } else {
      std::cerr << absl::StreamFormat("WARNING: Ignoring sequence '%s', could not determine its date\n", tip_desc.name);
    }
  }

  return result;
}

auto operator<<(std::ostream& os, stamp_version_into_log_file) -> std::ostream& {
  os << absl::StreamFormat("# Produced by delphy version %s (build %d, commit %s)\n",
                           k_delphy_version_string, k_delphy_build_number, k_delphy_commit_string);
  if (delphy_invoked_via_cli) {
    os << "# Invocation:";
    for (const auto& arg : delphy_cli_args) {
      // Naive w.r.t. quoting; if you do funny things, you're on your own
      os << " " << boost::replace_all_copy(arg, "--", "<DASH-DASH>");
    }
    os << "\n";
  }
  return os;
}

auto output_resolved_fasta(const Phylo_tree& tree, std::ostream& os) -> void {
  for (const auto& node : index_order_traversal(tree)) {
    if (tree.at(node).is_tip()) {
      
      os << '>' << tree.at(node).name << '\n';
      
      auto seq = view_of_sequence_at(tree, node);  // Missing sites resolved by inheriting from state at missations
      for (auto l = Site_index{0}; l != tree.num_sites(); ++l) {
        os << to_char(seq[l]);
      }
      os << '\n';
    }
  }
}

}  // namespace delphy
