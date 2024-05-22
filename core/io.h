#ifndef DELPHY_IO_H_
#define DELPHY_IO_H_

#include <ostream>

#include "sequence.h"

namespace delphy {

struct Fasta_entry {
  std::string id;
  std::string comments;
  Sequence sequence;
};

auto read_fasta(std::istream& is) -> std::vector<Fasta_entry>;

struct stamp_version_into_log_file{};  // dummy type for stream command
auto operator<<(std::ostream& os, stamp_version_into_log_file) -> std::ostream&;

}  // namespace delphy

#endif // DELPHY_IO_H_
