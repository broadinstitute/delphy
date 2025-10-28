#include "newick.h"

#include <boost/algorithm/string/replace.hpp>
#include <absl/strings/str_format.h>
#include <absl/strings/str_split.h>

namespace delphy {

// Newick tree grammar (ANTLR v4):
//
// Based loosely on https://github.com/antlr/grammars-v4/blob/master/newick/newick.g4.  Original license follows.
//  The closest things to a spec:
//  - https://evolution.genetics.washington.edu/phylip/newick_doc.html
//  - https://beast.community/nexus_metacomments
//
//BSD License
//Copyright (c) 2020, Tom Everett
//All rights reserved.
//Redistribution and use in source and binary forms, with or without
//modification, are permitted provided that the following conditions
//are met:
//1. Redistributions of source code must retain the above copyright
//   notice, this list of conditions and the following disclaimer.
//2. Redistributions in binary form must reproduce the above copyright
//   notice, this list of conditions and the following disclaimer in the
//   documentation and/or other materials provided with the distribution.
//3. Neither the name of Tom Everett nor the names of its contributors
//   may be used to endorse or promote products derived from this software
//   without specific prior written permission.
//THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//'AS IS' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
//A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
//HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
//SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
//LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
//DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
//THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
//OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// grammar Newick;
//
//tree_
//   : treeAttrs=ATTRIBUTES* root=node ';'
//   ;
//
//node
//   : ( '(' children+=node (',' children+=node)* ')' )?
//     name=( UNQUOTED_LABEL | QUOTED_LABEL | NUMBER )?
//     nodeAttrs=ATTRIBUTES*
//     ( ':' branchAttrs=ATTRIBUTES* length=NUMBER? )?
//   ;
//
//WS             : [ \r\n\t]+ -> skip ;
//ATTRIBUTES     : '[&' .*? ']' ;  // before COMMENT to recognize "[&hello]" as an ATTRIBUTES instead of a COMMENT
//COMMENT        : '[' .*? ']' -> skip ;
//NUMBER         : [+-]? [0-9]+ ('.' [0-9]+)? ([eE] [+-]? [0-9]+)? ;
//UNQUOTED_LABEL : (~[ \r\n\t()[\]':;,])+ ;  // Note: unlike spec, an unquoted label cannot begin with a digit
//QUOTED_LABEL   : '\'' ([^'] | '\'\'')* '\'';

auto Newick_lexer::to_string(Newick_lexer::Token_kind token_kind) -> std::string_view {
  switch (token_kind) {
    case Token_kind::k_attributes:
      return "ATTRIBUTES";
    case Token_kind::k_number:
      return "NUMBER";
    case Token_kind::k_unquoted_label:
      return "UNQUOTED_LABEL";
    case Token_kind::k_quoted_label:
      return "QUOTED_LABEL";
    case Token_kind::k_left_paren:
      return "'('";
    case Token_kind::k_comma:
      return "','";
    case Token_kind::k_right_paren:
      return "')'";
    case Token_kind::k_colon:
      return "':'";
    case Token_kind::k_semicolon:
      return "';'";
    case Token_kind::k_invalid:
      return "INVALID";
    case Token_kind::k_eof:
      return "EOF";
    default:
      throw std::logic_error(absl::StrFormat(
          "Unknown token kind: %d", static_cast<std::underlying_type<Token_kind>::type>(token_kind)));
  }
}

auto Newick_lexer::prime() -> void {
  peek_c = is_->get();
  started = true;
}

auto Newick_lexer::get() -> Newick_lexer::int_type {
  auto c = peek_c;
  if (c != eof) {
    lexeme_ += static_cast<char>(peek_c);
    peek_c = is_->get();
  }
  return c;
}

//WS             : [ \r\n\t]+ -> skip ;
//ATTRIBUTES     : '[&' .*? ']' ;  // before COMMENT to recognize "[&hello]" as an ATTRIBUTES instead of a COMMENT
//COMMENT        : '[' .*? ']' -> skip ;
//NUMBER         : [+-]? [0-9]+ ('.' [0-9]+)? ([eE] [+-]? [0-9]+)? ;
//UNQUOTED_LABEL : (~[ \r\n\t()[\]':;,])+ ;  // Note: unlike spec, an unquoted label cannot begin with a digit
//QUOTED_LABEL   : '\'' ([^'] | '\'\'')* '\'';
//
auto Newick_lexer::lex_next_token() -> Newick_lexer::Token_kind {
  if (!started) {
    prime();
  }

  while (true) {
    lexeme_.clear();
    switch (peek_c) {
      case eof:
        return Token_kind::k_eof;

      case ' ':
      case '\r':
      case '\n':
      case '\t':
        if (not lex_whitespace()) { return Token_kind::k_invalid; }
        break;  // Whitespace => loop and lex another token

      case '(':
        return get(), Token_kind::k_left_paren;
      case ')':
        return get(), Token_kind::k_right_paren;
      case ',':
        return get(), Token_kind::k_comma;
      case ':':
        return get(), Token_kind::k_colon;
      case ';':
        return get(), Token_kind::k_semicolon;

      case ']':
        return get(), Token_kind::k_invalid;

      case '[':
        // ATTRIBUTES or COMMENT
        if (not lex_comment_or_attributes()) { return Token_kind::k_invalid; }
        if (lexeme_.size() >= 2 && lexeme_[1] == '&') { return Token_kind::k_attributes; }
        break;  // Comment => loop and lex another token

      case '\'':
        return lex_quoted_label() ? Token_kind::k_quoted_label : Token_kind::k_invalid;

      default:
        if (std::isdigit(peek_c) || peek_c == '-' || peek_c == '+') {
          return lex_number() ? Token_kind::k_number : Token_kind::k_invalid;
        } else {
          return lex_unquoted_label() ? Token_kind::k_unquoted_label : Token_kind::k_invalid;
        }
    }
  }
}

//WS             : [ \r\n\t]+ -> skip ;
auto Newick_lexer::lex_whitespace() -> bool {
  auto empty = true;
  while (true) {
    switch (peek_c) {
      case ' ':
      case '\r':
      case '\n':
      case '\t':
        get();
        empty = false;
        break;
      default:
        return not empty;
    }
  }
}

//NUMBER         : [+-]? [0-9]+ ('.' [0-9]+)? ([eE] [+-]? [0-9]+)? ;
auto Newick_lexer::lex_number() -> bool {
  if (peek_c == '-' || peek_c == '+') { get(); }
  if (!std::isdigit(peek_c)) { return false; }
  while (std::isdigit(peek_c)) { get(); }
  if (peek_c == '.') {
    get();
    if (!std::isdigit(peek_c)) { return false; }
    while (std::isdigit(peek_c)) { get(); }
  }
  if (peek_c == 'e' || peek_c == 'E') {
    get();
    if (peek_c == '-' || peek_c == '+') { get(); }
    if (!std::isdigit(peek_c)) { return false; }
    while (std::isdigit(peek_c)) { get(); }
  }
  return true;
}

//QUOTED_LABEL   : '\'' ([^'] | '\'\'')* '\'';
auto Newick_lexer::lex_quoted_label() -> bool {
  if (peek_c != '\'') { return false; }
  get();
  while (true) {
    switch (get()) {
      case eof:
        return false;
      case '\'':
        if (peek_c == '\'') { get(); }
        else { return true; }
      default:
        continue;
    }
  }
}

//UNQUOTED_LABEL : (~[ \r\n\t()[\]':;,])+ ;  // Note: unlike spec, an unquoted label cannot begin with a digit
auto Newick_lexer::lex_unquoted_label() -> bool {
  auto empty = true;
  while (true) {
    switch (peek_c) {
      case ' ':
      case '\r':
      case '\n':
      case '\t':
      case '(':
      case ')':
      case '[':
      case ']':
      case '\'':
      case ':':
      case ';':
      case ',':
        return not empty;
      default:
        get();
        empty = false;
        break;
    }
  }
}

//ATTRIBUTES     : '[&' .*? ']' ;  // before COMMENT to recognize "[&hello]" as an ATTRIBUTES instead of a COMMENT
//COMMENT        : '[' .*? ']' -> skip ;
auto Newick_lexer::lex_comment_or_attributes() -> bool {
  get();
  while (peek_c != eof && peek_c != ']') { get(); }
  if (peek_c == ']') { return get(), true; }
  else { return false; }
}

static auto unquote(std::string_view quoted_label) -> std::string {
  assert(quoted_label.starts_with('\'') && quoted_label.ends_with('\''));
  auto unquoted_label = std::string{quoted_label.substr(1, quoted_label.size() - 2)};
  boost::replace_all(unquoted_label, "\'\'", "\'");
  return unquoted_label;
}

auto Newick_parser::match(Newick_lexer::Token_kind token_kind) -> void {
  if (lexer_.token_kind() != token_kind) {
    throw std::runtime_error(absl::StrFormat(
        "Expected %s, but got '%s' instead",
        Newick_lexer::to_string(token_kind), lexer_.lexeme()));
  }
  lexer_.next_token();
}

auto Newick_parser::maybe_match(Newick_lexer::Token_kind token_kind) -> bool {
  if (lexer_.token_kind() == token_kind) {
    lexer_.next_token();
    return true;
  } else {
    return false;
  }
}

auto Newick_parser::parse_tree() -> Newick_tree {
  auto tree = Newick_tree{};
  if (lexer_.token_kind() == Newick_lexer::Token_kind::k_attributes) {
    tree.tree_attributes = parse_attributes();
  }
  tree.root = parse_node(tree, k_no_node);
  match(Newick_lexer::Token_kind::k_semicolon);
  return tree;
}

auto Newick_parser::parse_node(Newick_tree& tree, Node_index parent) -> Node_index {
  auto node = tree.add_node();

  // ( '(' children+=node (',' children+=node)* ')' )?
  if (maybe_match(Newick_lexer::Token_kind::k_left_paren)) {
    tree.at(node).children = {};
    auto new_child = parse_node(tree, node);
    tree.at(node).children.push_back(new_child);
    while (true) {
      if (maybe_match(Newick_lexer::Token_kind::k_comma)) {
        new_child = parse_node(tree, node);
        tree.at(node).children.push_back(new_child);
      } else if (maybe_match(Newick_lexer::Token_kind::k_right_paren)) {
        break;
      } else {
        throw std::runtime_error(absl::StrFormat(
            "Expected ',' or ')', but got '%s' instead", lexer_.lexeme()));
      }
    }
  }

  // name=( UNQUOTED_LABEL | QUOTED_LABEL | NUMBER )?
  switch (lexer_.token_kind()) {
    case Newick_lexer::Token_kind::k_unquoted_label:
    case Newick_lexer::Token_kind::k_number: {
      tree.at(node).name = lexer_.lexeme();
      lexer_.next_token();
      break;
    }
    case Newick_lexer::Token_kind::k_quoted_label: {
      tree.at(node).name = unquote(lexer_.lexeme());
      break;
    }
    default: {
      // No name
      break;
    }
  }

  // nodeAttrs=ATTRIBUTES*
  if (lexer_.token_kind() == Newick_lexer::Token_kind::k_attributes) {
    tree.at(node).node_attributes = parse_attributes();
  }

  // ( ':' branchAttrs=ATTRIBUTES* length=NUMBER? )?
  if (maybe_match(Newick_lexer::Token_kind::k_colon)) {
    if (lexer_.token_kind() == Newick_lexer::Token_kind::k_attributes) {
      tree.at(node).branch_attributes = parse_attributes();
    }
    if (lexer_.token_kind() == Newick_lexer::Token_kind::k_number) {
      if (not absl::SimpleAtod(lexer_.lexeme(), &tree.at(node).branch_length)) {
        throw std::runtime_error(absl::StrFormat(
            "Could not parse branch length '%s'", lexer_.lexeme()));
      }
    }
    match(Newick_lexer::Token_kind::k_number);
  } else {
    tree.at(node).branch_length = 0.0;
  }

  tree.at(node).parent = parent;
  return node;
}

auto Newick_parser::parse_attributes() -> Newick_attributes {
  auto result = Newick_attributes{};
  while (lexer_.token_kind() == Newick_lexer::Token_kind::k_attributes) {
    assert(lexer_.lexeme().starts_with("[&") && lexer_.lexeme().ends_with(']'));
    auto attrs_text = lexer_.lexeme().substr(2, lexer_.lexeme().size() - 3);  // "[&" = 2 chars, "[&]" = 3 chars
    for (auto& attr : absl::StrSplit(attrs_text, ';')) {
      auto equal_it = std::ranges::find(attr, '=');
      if (equal_it != attr.end()) {
        result.try_emplace(std::string{attr.begin(), equal_it},
                           std::string{equal_it + 1, attr.end()});
      } else {
        result.try_emplace(std::string{attr}, "");
      }
    }
    lexer_.next_token();
  }
  return result;
}

}  // namespace delphy
