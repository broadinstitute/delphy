#include "io.h"

#include <iostream>

#include <boost/algorithm/string/replace.hpp>

#include "version.h"
#include "cmdline.h"
#include "dates.h"
#include "phylo_tree_calc.h"

namespace delphy {

auto read_fasta(
    std::istream& is,
    const std::function<void(int, std::size_t)>& progress_hook,
    const std::function<void(const std::string&, const Sequence_warning&)>& warning_hook)
    -> std::vector<Fasta_entry> {
  
  auto in_seq = false;
  auto cur_id = std::string{};
  auto cur_seq = Sequence{};
  auto ignore_this_seq = false;
  auto seqs_so_far = 0;

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
        if (not ignore_this_seq) {
          result.push_back(Fasta_entry{std::move(cur_id), std::move(cur_seq)});
        }
        ++seqs_so_far;
        progress_hook(seqs_so_far, bytes_so_far);
        cur_seq.clear();
        ignore_this_seq = false;
      }

      // 2025-04-24: We used to find the first space, if any, and treat everything
      // beyond that as a comment (as in https://www.ncbi.nlm.nih.gov/genbank/fastaformat/).
      // However, this has caused no end of complaints as people produce FASTA files with
      // sequence IDs that contain spaces.  Moreover, both BEAST and IQ-Tree
      // read IDs as everything between the '>' and the end of the line.  So from now on,
      // we'll behave like BEAST and IQ-Tree
      auto id_begin = line.begin() + 1;  // + 1 => skip over initial '>'
      auto id_end = line.end();
      while (id_begin < id_end && std::isspace(*id_begin)) { ++id_begin; }
      while (id_begin < id_end && std::isspace(*(id_end-1))) { --id_end; }
      cur_id.assign(id_begin, id_end);

      in_seq = true;
      continue;
    }

    if (!in_seq) {
      // Sequence data preceding an initial identifier line (e.g., '>SeqName') is read
      // into a sequence named ""
      in_seq = true;
    }

    if (not ignore_this_seq) {
      for (auto c : line) {
        if (isspace(c)) {
          continue;  // Ignore spaces and Unix/Mac/DOS line termination nonsense
        }
        
        if (is_valid_seq_letter(c)) {
          cur_seq.push_back(to_seq_letter(c));
        } else {
          ignore_this_seq = true;
          warning_hook(cur_id, Sequence_warnings::Invalid_state{c});
        }
      }
    }
  }
  if (in_seq) {
    if (not ignore_this_seq) {
      result.push_back(Fasta_entry{std::move(cur_id), std::move(cur_seq)});
    }
    ++seqs_so_far;
    progress_hook(seqs_so_far, bytes_so_far);
  }

  return result;
}

auto read_maple(
    std::istream& is,
    const std::function<void(int, std::size_t)>& progress_hook,
    const std::function<void(const std::string&, const Sequence_warning&)>& warning_hook)
    -> Maple_file {
  
  auto result = Maple_file{};
  auto bytes_so_far = std::size_t{0};

  auto line = std::string{};
  if (not getline(is, line)) {
    throw std::runtime_error("ERROR: Unexpected EOF while reading MAPLE file reference sequence");
  }
  bytes_so_far += std::ssize(line) + 1;  // `+ 1` from end-of-line character
  progress_hook(std::ssize(result.tip_descs), bytes_so_far);
  
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
  // Any ambiguous sites in reference are silently turned into A, but we flag any sequence that includes
  // a mutation to those sites
  auto ambiguous_sites = absl::flat_hash_set<Site_index>{};
  if (not getline(is, line)) {
    throw std::runtime_error("ERROR: Unexpected EOF while reading MAPLE file reference sequence");
  }
  do {
    bytes_so_far += std::ssize(line) + 1;  // `+ 1` from end-of-line character
    if (line.empty()) {
      continue;
    }
    if (line[0] == '>') {
      break;
    }
    
    for (auto c : line) {
      if (std::isspace(c)) {
          continue;  // Ignore spaces and Unix/Mac/DOS line termination nonsense
      }

      if (is_valid_seq_letter(c)) {
        if (is_valid_real_seq_letter(c)) {
          auto s = char_to_real_seq_letter(c);
          result.ref_sequence.push_back(s);
        } else {
          ambiguous_sites.insert(std::ssize(result.ref_sequence));
          result.ref_sequence.push_back(Real_seq_letter::A);
        }
      } else {
        throw std::runtime_error(absl::StrFormat("Reference sequence has invalid state '%c' at site %d",
                                                 c, std::ssize(result.ref_sequence) + 1));
      }
    }
  } while (getline(is, line));
  progress_hook(std::ssize(result.tip_descs), bytes_so_far);

  auto L = std::ssize(result.ref_sequence);

  // Now read each sequence in turn
  while (is) {
    if (line.empty() || line[0] != '>') {
      throw std::runtime_error(absl::StrFormat(
          "Expected sequence id line to start with '>', but saw this instead: %s",
          line));
    }

    auto tip_desc = Tip_desc{};
    auto ignore_this_tip = false;

    // 2025-04-24: See analogous comment in read_fasta above wrt ID vs comment in FASTA files
    auto id_begin = line.begin() + 1;  // + 1 => skip over initial '>'
    auto id_end = line.end();
    while (id_begin < id_end && std::isspace(*id_begin)) { ++id_begin; }
    while (id_begin < id_end && std::isspace(*(id_end-1))) { --id_end; }
    tip_desc.name.assign(id_begin, id_end);

    auto maybe_t_range = extract_date_range_from_sequence_id(tip_desc.name);
    if (not maybe_t_range.has_value()) {
      warning_hook(tip_desc.name, Sequence_warnings::No_valid_date{});
      ignore_this_tip = true;
    } else {
      auto [t_min, t_max] = maybe_t_range.value();
      tip_desc.t_min = t_min;
      tip_desc.t_max = t_max;
    }

    while(getline(is, line)) {
      bytes_so_far += std::ssize(line) + 1;  // `+ 1` from end-of-line character
      if (line.empty()) {
        continue;  // ignore empty lines
      }
      if (line[0] == '>') {
        break;  // on to next tip
      }

      auto ss = std::istringstream{line};
      auto c = char{};
      ss >> c;
      if (not (ss && is_valid_seq_letter(c))) {
        ignore_this_tip = true;
        warning_hook(tip_desc.name, Sequence_warnings::Invalid_state{c});
        continue;
      }
      
      auto s = to_seq_letter(c);
      if (is_ambiguous(s)) {
        // Gap or Missing Data
        auto start_site = Site_index{};  ss >> start_site;  --start_site;  // 0-based indexing
        auto valid_start_site = not ss.fail();  // differs from ss.good() in that EOF is not a failing condition
        auto gap_len = Site_index{1};  ss >> gap_len;  // Default to gap size of 1 if none indicated (e.g., for '?')
        auto end_site = start_site + gap_len;

        if (valid_start_site
            && (0 <= start_site && start_site < L) && (0 < end_site && end_site <= L)
            && (start_site < end_site)) {
          tip_desc.missations.intervals.insert(Site_interval{start_site, end_site});
        } else {
          ignore_this_tip = true;
          warning_hook(tip_desc.name, Sequence_warnings::Invalid_gap{start_site, end_site});
          continue;
        }
      } else {
        // Mutation w.r.t. reference
        auto l = Site_index{};  ss >> l;  --l;  // 0-based indexing
        auto a = (0 <= l && l < L) ? result.ref_sequence.at(l) : Real_seq_letter::A;
        auto b = to_real_seq_letter(s);
        if (ss
            && (0 <= l && l < L)
            && ((a != b) || (a == Real_seq_letter::T)) // Sometimes MAPLE files have spurious "T->U" mutations 
            && not ambiguous_sites.contains(l)) {
          if (not (a == Real_seq_letter::T && b == Real_seq_letter::T)) {  // See above
            tip_desc.seq_deltas.push_back({l, a, b});
          }
        } else {
          ignore_this_tip = true;
          warning_hook(tip_desc.name, Sequence_warnings::Invalid_mutation{a, l, b});
          continue;
        }
      }
    }

    // Only include sequence if we haven't ignored it owing to a parse error above
    if (not ignore_this_tip) {
      result.tip_descs.push_back(std::move(tip_desc));
    }
    
    progress_hook(std::ssize(result.tip_descs), bytes_so_far);
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
      
      os << '>' << tree.at(node).name << '|' << to_iso_date(tree.at(node).t) << '\n';
      
      auto seq = view_of_sequence_at(tree, node);  // Missing sites resolved by inheriting from state at missations
      for (auto l = Site_index{0}; l != tree.num_sites(); ++l) {
        os << to_char(seq[l]);
      }
      os << '\n';
    }
  }
}

}  // namespace delphy
