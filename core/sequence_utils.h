#ifndef DELPHY_SEQUENCE_UTILS_H_
#define DELPHY_SEQUENCE_UTILS_H_

#include <ranges>

#include "sequence.h"
#include "mutations.h"

namespace delphy {

// Sequence warning objects
namespace Sequence_warnings {
struct No_valid_date {};
struct Ambiguity_precision_loss {
  Seq_letter original_state;
  Site_index site;
};
struct Invalid_state {
  char state_letter;
};
struct Invalid_gap {
  Site_index start_site;
  Site_index end_site;
};
struct Invalid_mutation {
  Real_seq_letter from;
  Site_index site;
  Real_seq_letter to;
};
}
using Sequence_warning = std::variant<
  Sequence_warnings::No_valid_date,
  Sequence_warnings::Ambiguity_precision_loss,
  Sequence_warnings::Invalid_state,
  Sequence_warnings::Invalid_gap,
  Sequence_warnings::Invalid_mutation
>;

auto default_sequence_warning_hook(const std::string& seq_id, const Sequence_warning& warning) -> void;

// Deduces an unambiguous consensus sequence among a range of (possibly ambiguous) base sequences.
// At a site l, the consensus state is the most common unambiguous state at l (ties are broken in favor
// of A, then C, then G, then T).  If the state is always ambiguous, choose A arbtrarily.
template<std::ranges::range Sequences>
requires std::convertible_to<std::ranges::range_reference_t<Sequences>, const Sequence&>
auto deduce_consensus_sequence(const Sequences& sequences) -> Real_sequence {
  auto first = std::ranges::begin(sequences);
  auto last = std::ranges::end(sequences);

  if (first == last) {
    throw std::invalid_argument("deduce_consensus_sequence passed an empty range of sequences");
  }

  auto L = static_cast<Site_index>(std::ssize(*first));
  auto state_counts = std::vector<Seq_vector<int>>(L, Seq_vector<int>(0));
  for (auto it = first; it != last; ++it) {
    const auto& cur_seq = *it;
    for (auto l = Site_index{0}; l != L; ++l) {
      auto s = cur_seq[l];
      if (not is_ambiguous(s)) {
        ++state_counts[l][to_real_seq_letter(s)];
      }
    }
  }
  
  auto consensus_sequence = Real_sequence(L, Real_seq_letter::A);
  for (auto l = Site_index{0}; l != L; ++l) {
    if (auto most_common_freq = std::ranges::max(state_counts[l]);
        most_common_freq != 0) {
      for (auto a : k_all_real_seq_letters) {
        if (state_counts[l][a] == most_common_freq) {
          consensus_sequence[l] = a;
          break;
        }
      }
    }
  }

  return consensus_sequence;
}

// Summarize delta from an unambiguous reference sequence to a (possibly ambiguous) actual sequence.
// Any ambiguity is turned into an N (completely ambiguous).
struct Delta_from_reference {
  std::vector<Seq_delta> seq_deltas;
  Missation_map<> missations;
};
auto calculate_delta_from_reference(
    const std::string& seq_id,
    const Sequence& seq,
    const Real_sequence& ref_seq,
    const std::function<void(const std::string&, const Sequence_warning&)>& warning_hook
    = default_sequence_warning_hook)
    -> Delta_from_reference;

auto extract_date_from_sequence_id(const std::string_view id) -> std::optional<double>;

}  // namespace delphy

#endif // DELPHY_SEQUENCE_UTILS_H_
