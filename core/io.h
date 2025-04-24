#ifndef DELPHY_IO_H_
#define DELPHY_IO_H_

#include <ostream>

#include "sequence.h"
#include "sequence_utils.h"
#include "phylo_tree.h"

namespace delphy {

struct Fasta_entry {
  std::string id;
  Sequence sequence;
};

auto read_fasta(
    std::istream& is,
    const std::function<void(int, std::size_t)>& progress_hook = [](int, std::size_t){},
    const std::function<void(const std::string&, const Sequence_warning&)>& warning_hook
    = default_sequence_warning_hook)
    -> std::vector<Fasta_entry>;

struct Maple_file {
  Real_sequence ref_sequence;
  std::vector<Tip_desc> tip_descs;
};

auto read_maple(
    std::istream& is,
    const std::function<void(int, std::size_t)>& progress_hook = [](int, std::size_t){},
    const std::function<void(const std::string&, const Sequence_warning&)>& warning_hook
    = default_sequence_warning_hook)
    -> Maple_file;

struct stamp_version_into_log_file{};  // dummy type for stream command
auto operator<<(std::ostream& os, stamp_version_into_log_file) -> std::ostream&;

auto output_resolved_fasta(const Phylo_tree& tree, std::ostream& os) -> void;

}  // namespace delphy

#endif // DELPHY_IO_H_
