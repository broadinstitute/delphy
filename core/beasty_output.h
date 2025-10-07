#ifndef DELPHY_BEASTY_OUTPUT_H_
#define DELPHY_BEASTY_OUTPUT_H_

#include <fstream>

#include "run.h"

namespace delphy {

class Beasty_log_output_version_impl;  // Opaque implementer for each BEAST version

class Beasty_log_output {
 public:
  Beasty_log_output(std::ostream* os, std::string beast_version, bool own_stream = false);
  ~Beasty_log_output();

  auto os() -> std::ostream* { return os_; }
  
  auto output_headers(const Run& run) -> void;
  auto output_log(const Run& run) -> void;
  auto output_footers(const Run& run) -> void;
  auto flush() -> void;

 private:
  std::ostream* os_;
  std::string beast_version_;
  std::unique_ptr<Beasty_log_output_version_impl> impl_;
  bool own_stream_;
};

class Beasty_trees_output {
 public:
  Beasty_trees_output(std::ostream* os, std::string beast_version, bool own_stream = false);
  ~Beasty_trees_output();

  auto os() -> std::ostream* { return os_; }
  
  auto output_headers(const Run& run) -> void;
  auto output_tree(const Run& run) -> void;
  auto output_footers(const Run& run) -> void;
  auto flush() -> void;

 private:
  std::ostream* os_;
  std::string beast_version_;
  bool own_stream_;
  Node_vector<Node_index> node_to_tip_{};

  void output_newick_tree(const Phylo_tree& tree);
};

}  // namespace delphy

#endif // DELPHY_BEASTY_OUTPUT_H_
