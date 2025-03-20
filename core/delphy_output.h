#ifndef DELPHY_DELPHY_OUTPUT_H_
#define DELPHY_DELPHY_OUTPUT_H_

#include <fstream>

#include "api.h"
#include "run.h"

namespace delphy {

class Delphy_output {
 public:
  Delphy_output(std::ostream* os, bool own_stream = false);
  ~Delphy_output();

  auto os() -> std::ostream* { return os_; }

  auto set_dphy_metadata_blob(std::string dphy_metadata_blob) -> void;
  
  auto output_preamble(const Run& run, int64_t steps_per_sample) -> void;
  auto output_state(const Run& run) -> void;
  auto output_epilog() -> void;
  auto flush() -> void;

 private:
  std::ostream* os_;
  bool own_stream_;
  std::string dphy_metadata_blob_;

  auto write_string(std::string_view sv) -> void;
  auto write_uint32(uint32_t num) -> void;
  auto write_uint64(uint64_t num) -> void;
  auto write_float32(float num) -> void;
  
  auto get_flatbuffer_size(const flatbuffers::DetachedBuffer& buffer) const -> uint32_t;
  auto write_flatbuffer_size(const flatbuffers::DetachedBuffer& buffer) -> void;
  auto write_flatbuffer_contents(const flatbuffers::DetachedBuffer& buffer) -> void;
};

}  // namespace delphy

#endif // DELPHY_DELPHY_OUTPUT_H_
