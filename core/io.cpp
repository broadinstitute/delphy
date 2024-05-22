#include "io.h"

#include <iostream>

#include <boost/algorithm/string/replace.hpp>

#include "dates.h"
#include "version.h"
#include "cmdline.h"

namespace delphy {

auto read_fasta(std::istream& is) -> std::vector<Fasta_entry> {
  auto in_seq = false;
  auto cur_id = std::string{};
  auto cur_comments = std::string{};
  auto cur_seq = Sequence{};

  auto result = std::vector<Fasta_entry>{};

  for (auto line = std::string{}; getline(is, line);) {
    if (line.empty()) {
      continue;  // ignore empty lines
    }

    if (line[0] == ';') {
      continue;  // comment line, ignore
    }

    if (line[0] == '>') {
      if (in_seq) {
        result.push_back(Fasta_entry{std::move(cur_id), std::move(cur_comments), std::move(cur_seq)});
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

}  // namespace delphy
