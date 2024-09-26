#include "delphy_output.h"

#include <ranges>

#include "phylo_tree_calc.h"
#include "io.h"
#include "version.h"

namespace delphy {

static const auto k_app_magic_number = std::string{"DPHY"};
static const auto k_save_format_version = 3;
static const auto k_no_more_trees = 0;
static const auto k_default_dphy_config = std::string{
    "{"
    "\"confidence\":90,"
    "\"topology\":0,"
    "\"presentation\":0,"
    "\"spacing\":0,"
    "\"colorBy\":0,"
    "\"burnin\":0,"
    "\"metadataPresent\":0,"
    "\"metadataText\":null,"
    "\"metadataFile\":null,"
    "\"metadataDelimiter\":null,"
    "\"selectedMDField\":-1,"
    "\"metadataColors\":{}"
    "}"};

Delphy_output::Delphy_output(std::ostream* os, bool own_stream)
    : os_{os}, own_stream_{own_stream} {}

Delphy_output::~Delphy_output() {
  if (own_stream_) {
    delete os_;
  }
}

auto
Delphy_output::write_string(std::string_view sv) -> void {
  write_uint32(sv.length());
  for (const auto& c : sv) {
    os_->put(c);
  }
}

auto
Delphy_output::write_uint32(uint32_t num) -> void {
  (*os_)
      .put((num >>  0) & 0xff)
      .put((num >>  8) & 0xff)
      .put((num >> 16) & 0xff)
      .put((num >> 24) & 0xff);
}

auto
Delphy_output::write_uint64(uint64_t num) -> void {
  (*os_)
      .put((num >>  0) & 0xff)
      .put((num >>  8) & 0xff)
      .put((num >> 16) & 0xff)
      .put((num >> 24) & 0xff)
      .put((num >> 32) & 0xff)
      .put((num >> 40) & 0xff)
      .put((num >> 48) & 0xff)
      .put((num >> 56) & 0xff);
}

auto
Delphy_output::write_float32(float num) -> void {
  write_uint32(std::bit_cast<uint32_t>(num));
}

auto
Delphy_output::get_flatbuffer_size(const flatbuffers::DetachedBuffer& buffer) const -> uint32_t {
  return flatbuffers::GetPrefixedSize(buffer.data()) + 4;  // +4 => account for the size prefix!
}

auto
Delphy_output::write_flatbuffer_size(const flatbuffers::DetachedBuffer& buffer) -> void {
  write_uint32(get_flatbuffer_size(buffer));
}

auto
Delphy_output::write_flatbuffer_contents(const flatbuffers::DetachedBuffer& buffer) -> void {
  os_->write(reinterpret_cast<const char*>(buffer.data()), get_flatbuffer_size(buffer));
}

auto Delphy_output::output_preamble(const Run& run, int64_t steps_per_sample) -> void {
  // See docs/dphy_file_format.md
  for (const auto& c : k_app_magic_number) {
    os_->put(c);
  }
  write_uint32(k_save_format_version);
  write_string(k_delphy_version_string);
  write_uint32(k_delphy_build_number);
  write_string(k_delphy_commit_string);
  write_uint32(0);  // knee_index
  write_uint32(steps_per_sample);  // Steps per sample
  write_uint32(run.alpha_move_enabled() ? 1 : 0);
  write_uint32(run.mpox_hack_enabled() ? 1 : 0);
  write_uint32(run.mu_move_enabled() ? 1 : 0);
  // Other flags like run.final_pop_size_move_enabled() and run.pop_growth_rate_move_enabled() derived from first Params
  write_float32(run.mu());
  
  auto api_tree_info = phylo_tree_to_api_tree_info(run.tree());
  write_flatbuffer_size(api_tree_info);
  write_flatbuffer_contents(api_tree_info);
}

auto Delphy_output::output_state(const Run& run) -> void {
  auto api_tree = phylo_tree_to_api_tree(run.tree());
  auto api_params = run_to_api_params(run);
  write_flatbuffer_size(api_tree);
  write_flatbuffer_size(api_params);
  write_flatbuffer_contents(api_tree);
  write_flatbuffer_contents(api_params);
}

auto Delphy_output::output_epilog() -> void {
  auto tree_end_position = os_->tellp();

  write_uint32(k_no_more_trees);
  write_string(k_default_dphy_config);
  write_uint64(tree_end_position);
}

auto Delphy_output::flush() -> void {
  // Ensure partial .dphy file is valid (epilog is always the same size, so writing a new state
  // here or simply ending the file here will fully overwrite this temporary epilog)
  auto save_position = os_->tellp();
  output_epilog();
  os_->seekp(save_position);
  
  os_->flush();
}

}  // namespace delphy
