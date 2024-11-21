#include "dates.h"
#include "sequence_utils.h"

namespace delphy {

auto default_sequence_warning_hook(const std::string& seq_id, const Sequence_warning& warning) -> void {
  auto msg = std::string{
    std::visit(estd::overloaded{
      [](const Sequence_warnings::No_valid_date&) -> std::string {
        return "Sequence ignored because its exact date could not be determined "
            "(no ambiguity is allowed for now, use YYYY-MM-DD format)";
      },
      [](const Sequence_warnings::Ambiguity_precision_loss& w) -> std::string {
        return absl::StrFormat("ambiguous state %c at site %d changed to N", to_char(w.original_state), w.site+1);
      }
    }, warning)
  };
  std::cerr << absl::StreamFormat("WARNING (sequence '%s'): %s\n", seq_id, msg);
}

auto calculate_delta_from_reference(
    const std::string& seq_id,
    const Sequence& seq,
    const Real_sequence& ref_seq,
    const std::function<void(const std::string&, const Sequence_warning&)>& warning_hook)
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
        warning_hook(seq_id,  // No sequence ID known
                     Sequence_warnings::Ambiguity_precision_loss{
                       .original_state = seq[l],
                       .site = l
                     });
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

auto extract_date_from_sequence_id(const std::string_view id) -> std::optional<double> {
  // We assume that sequence IDs look like this (these are real examples):
  //
  //   hCoV-19/England/PLYM-3258B175/2022|EPI_ISL_15330011|2022-10-01
  //   hRSV-A-England-160340212-2016-EPI_ISL_11428309-2016-01-19
  //
  // So we check that there's something date-like at the end of the string, and parse that.
  // Bork badly if the date isn't there
  if (id.length() < (1 + 4 + 1 + 2 + 1 + 2)) {
    // Not long enough
    return {};
  }
  auto date_plus_bar = id.substr(id.length() - (1 + 4 + 1 + 2 + 1 + 2));
  if (date_plus_bar[0] != '|' && date_plus_bar[0] != '-') {
    // No separator before date
    return {};
  }
  auto date_str = date_plus_bar.substr(1);
  auto valid =
      std::ranges::all_of(date_str.substr(0, 4), isdigit) &&
          date_str[4] == '-' &&
          std::ranges::all_of(date_str.substr(5, 2), isdigit) &&
          date_str[7] == '-' &&
          std::ranges::all_of(date_str.substr(8, 2), isdigit);
  if (not valid) {
    return {};
  }
  return parse_iso_date(date_str);
}

}  // namespace delphy
