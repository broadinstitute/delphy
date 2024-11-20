#include "sequence_utils.h"

namespace delphy {

auto calculate_delta_from_reference(
    const Sequence& seq,
    const Real_sequence& ref_seq,
    const std::function<void(const std::string&)>& warning_hook)
    -> Delta_from_reference {
  
  if (std::ssize(seq) != std::ssize(ref_seq)) {
    throw std::invalid_argument(absl::StrFormat(
        "calculate_delta_from_reference(): seq (%d bases) and ref_seq (%d bases) differ in size",
        std::ssize(seq), std::ssize(ref_seq)));
  }
  
  auto result = Delta_from_reference{};

  auto L = static_cast<Site_index>(std::ssize(seq));
  for (auto l = Site_index{0}; l != L; ++l) {
    if (is_ambiguous(seq[l])) {
      if (seq[l] != Seq_letters::N) {
        warning_hook(absl::StrFormat("ambiguous state %c at site %d changed to N\n", to_char(seq[l]), l+1));
      }
      result.missations.insert(Missation{l, ref_seq[l]}, ref_seq);
    } else {
      auto from = ref_seq[l];
      auto to = to_real_seq_letter(seq[l]);
      if (from != to) { result.seq_deltas.push_back({l, from, to}); }
    }
  }

  return result;
}

}  // namespace delphy
