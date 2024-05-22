#ifndef DELPHY_NEWICK_H_
#define DELPHY_NEWICK_H_

#include <absl/container/flat_hash_map.h>

#include "tree.h"

namespace delphy {

using Newick_attributes = absl::flat_hash_map<std::string, std::string>;

class Newick_node : public Nary_node {
 public:
  std::string name{};
  Newick_attributes node_attributes{};
  Newick_attributes branch_attributes{};
  double branch_length{0.0};
};

class Newick_tree : public Tree<Newick_node> {
 public:
  Newick_attributes tree_attributes{};
};

class Newick_lexer {
 public:
  enum class Token_kind {
    k_attributes,
    k_number,
    k_unquoted_label,
    k_quoted_label,
    k_left_paren,
    k_comma,
    k_right_paren,
    k_colon,
    k_semicolon,
    k_invalid,
    k_eof
  };
  static auto to_string(Newick_lexer::Token_kind token_kind) -> std::string_view;

  explicit Newick_lexer(std::istream& is) : is_{&is} {}

  auto next_token() -> Token_kind { return token_kind_ = lex_next_token(); }

  // Info on last token read in
  auto token_kind() -> Token_kind { return token_kind_; }
  auto lexeme() -> std::string_view { return lexeme_; }  // Invalidated by subsequent call to next()

 private:
  std::istream* is_;
  Token_kind token_kind_{Token_kind::k_invalid};
  std::string lexeme_{};  // backing string for string_views in NewickToken

  using int_type = std::istream::int_type;
  static constexpr int_type eof = std::istream::traits_type::eof();

  bool started = false;
  int_type peek_c = eof;
  
  auto prime() -> void;
  auto get() -> int_type;

  auto lex_next_token() -> Token_kind;
  auto lex_whitespace() -> bool;
  auto lex_number() -> bool;
  auto lex_quoted_label() -> bool;
  auto lex_unquoted_label() -> bool;
  auto lex_comment_or_attributes() -> bool;
};

class Newick_parser {
 public:
  explicit Newick_parser(std::istream& is) : lexer_{is} { lexer_.next_token(); }

  auto parse_tree() -> Newick_tree;

 private:
  Newick_lexer lexer_;

  auto match(Newick_lexer::Token_kind token_kind) -> void;
  auto maybe_match(Newick_lexer::Token_kind token_kind) -> bool;
  auto parse_node(Newick_tree& tree, Node_index parent) -> Node_index;
  auto parse_attributes() -> Newick_attributes;
};

}  // namespace delphy

#endif // DELPHY_NEWICK_H_
