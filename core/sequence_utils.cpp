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
      },
      [](const Sequence_warnings::Invalid_state& w) -> std::string {
        return absl::StrFormat("invalid state '%c'", w.state_letter);
      },
      [](const Sequence_warnings::Invalid_gap& w) -> std::string {
        return absl::StrFormat("invalid gap from %d to %d", w.start_site+1, w.end_site+1);
      },
      [](const Sequence_warnings::Invalid_mutation& w) -> std::string {
        return absl::StrFormat("invalid mutation from %c to %c at site %d", to_char(w.from), to_char(w.to), w.site+1);
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

static auto is_valid_year_str(std::string_view sv) -> bool {
  return
      sv.length() == (4) &&
      std::ranges::all_of(sv.substr(0, 4), isdigit);
}

static auto is_valid_month_str(std::string_view sv) -> bool {
  return
      sv.length() == (4+1+2) &&
      std::ranges::all_of(sv.substr(0, 4), isdigit) &&
      sv[4] == '-' &&
      std::ranges::all_of(sv.substr(5, 2), isdigit);
}

static auto is_valid_date_str(std::string_view sv) -> bool {
  return
      sv.length() == (4+1+2+1+2) &&
      std::ranges::all_of(sv.substr(0, 4), isdigit) &&
      sv[4] == '-' &&
      std::ranges::all_of(sv.substr(5, 2), isdigit) &&
      sv[7] == '-' &&
      std::ranges::all_of(sv.substr(8, 2), isdigit);
}

static auto is_valid_date_range_str(std::string_view sv) -> bool {
  return
      sv.length() == (4+1+2+1+2+1+4+1+2+1+2) &&
      is_valid_date_str(sv.substr(0, 4+1+2+1+2)) &&
      sv[4+1+2+1+2] == '/' &&
      is_valid_date_str(sv.substr(4+1+2+1+2+1, 4+1+2+1+2));
}

auto extract_date_range_from_sequence_id(const std::string_view id) -> std::optional<std::pair<double, double>> {
  // We assume that sequence IDs look like this (these are real examples):
  //
  //   hCoV-19/England/PLYM-3258B175/2022|EPI_ISL_15330011|2022-10-01
  //   hRSV-A-England-160340212-2016-EPI_ISL_11428309-2016-01-19
  //
  // So we check that there's something date-range-like at the end of the string, and parse that.
  // Bork badly if the date range isn't there
  //
  // Date ranges can be in the following formats:
  //
  //   2024-12-26             (exact)
  //   2024-12                (the whole month)
  //   2024                   (the whole year)
  //   2024-11-05/2024-12-26  (arbitrary date range, endpoints *must* be YYYY-MM-DD)
  //

  const auto len_arb_date_range = 4+1+2+1+2+1+4+1+2+1+2;
  const auto len_arb_date       = 4+1+2+1+2;
  const auto len_arb_month      = 4+1+2;
  const auto len_arb_year       = 4;
  
  // At any point, parse_iso_date can fail with an exception; we just report that it failed and swallow the error
  try {
    // Try arbitrary date range first
    do {
      if (id.length() < 1 + len_arb_date_range) {
        break;
      }
      
      auto date_range_plus_bar = id.substr(id.length() - (1 + len_arb_date_range));
      if (date_range_plus_bar[0] != '|' && date_range_plus_bar[0] != '-') {
        // No separator before date range
        break;
      }
      if (not is_valid_date_range_str(date_range_plus_bar.substr(1, len_arb_date_range))) {
        break;
      }
      
      auto min_date_str = date_range_plus_bar.substr(1, len_arb_date);
      auto max_date_str = date_range_plus_bar.substr(1+len_arb_date+1, len_arb_date);
      
      auto t_min = parse_iso_date(min_date_str);
      auto t_max = parse_iso_date(max_date_str);
      return {{t_min, t_max}};
    } while(false);

    // Next try exact date
    do {
      if (id.length() < 1 + len_arb_date) {
        break;
      }
      
      auto date_plus_bar = id.substr(id.length() - (1 + len_arb_date));
      if (date_plus_bar[0] != '|' && date_plus_bar[0] != '-') {
        // No separator before date
        break;
      }
      if (not is_valid_date_str(date_plus_bar.substr(1, len_arb_date))) {
        break;
      }
      
      auto date_str = date_plus_bar.substr(1, len_arb_date);
      auto t_exact = parse_iso_date(date_str);
      return {{t_exact, t_exact}};
    } while(false);

    // Next try month
    do {
      if (id.length() < 1 + len_arb_month) {
        break;
      }
      
      auto month_plus_bar = id.substr(id.length() - (1 + len_arb_month));
      if (month_plus_bar[0] != '|' && month_plus_bar[0] != '-') {
        // No separator before month
        break;
      }
      if (not is_valid_month_str(month_plus_bar.substr(1, len_arb_month))) {
        break;
      }
      
      auto month_str = month_plus_bar.substr(1, len_arb_month);
      auto [t_min, t_max] = parse_iso_month(month_str);
      return {{t_min, t_max}};
    } while(false);

    // Next try year
    do {
      if (id.length() < 1 + len_arb_year) {
        break;
      }
      
      auto year_plus_bar = id.substr(id.length() - (1 + len_arb_year));
      if (year_plus_bar[0] != '|' && year_plus_bar[0] != '-') {
        // No separator before year
        break;
      }
      if (not is_valid_year_str(year_plus_bar.substr(1, len_arb_year))) {
        break;
      }
      
      auto year_str = year_plus_bar.substr(1, len_arb_year);
      auto [t_min, t_max] = parse_iso_year(year_str);
      return {{t_min, t_max}};
    } while(false);

  } catch(...) {
    // Horribly formatted date, bail
  }
  
  // Nothing worked
  return {};
}

}  // namespace delphy
