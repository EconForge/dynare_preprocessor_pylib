/*
 * Copyright © 2019 Dynare Team
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _MACRO_DRIVER_HH
#define _MACRO_DRIVER_HH

#ifdef _PARSING_DRIVER_HH
# error Impossible to include both ../ParsingDriver.hh and Driver.hh
#endif

#include "Parser.hh"
#include "Environment.hh"
#include "Expressions.hh"

#include <stack>
#include <filesystem>

// Declare TokenizerFlexLexer class
#ifndef __FLEX_LEXER_H
# define yyFlexLexer TokenizerFlexLexer
# include <FlexLexer.h>
# undef yyFlexLexer
#endif

namespace macro
{
  /* The lexer class
   * It was necessary to subclass the TokenizerFlexLexer class generated by Flex,
   * since the prototype for TokenizerFlexLexer::yylex() was not convenient.
   */
  class TokenizerFlex : public TokenizerFlexLexer
  {
  public:
    TokenizerFlex(istream *in) : TokenizerFlexLexer{in} { }
    TokenizerFlex(const TokenizerFlex &) = delete;
    TokenizerFlex(TokenizerFlex &&) = delete;
    TokenizerFlex & operator=(const TokenizerFlex &) = delete;
    TokenizerFlex & operator=(TokenizerFlex &&) = delete;

    //! The main lexing function
    Tokenizer::parser::token_type lex(Tokenizer::parser::semantic_type *yylval,
                                      Tokenizer::parser::location_type *yylloc,
                                      macro::Driver &driver);
  };

  //! Implements the macro expansion using a Flex scanner and a Bison parser
  class Driver
  {
  public:
    Environment &env;
  private:
    bool no_line_macro;
    vector<DirectivePtr> statements;
    stack<vector<DirectivePtr>> directive_stack;
  public:
    Driver(Environment &env_arg, bool no_line_macro_arg) :
      env{env_arg}, no_line_macro(no_line_macro_arg) { }
    Driver(const Driver &) = delete;
    Driver(Driver &&) = delete;
    Driver & operator=(const Driver &) = delete;
    Driver & operator=(Driver &&) = delete;

    //! Exception thrown when value of an unknown variable is requested
    class UnknownVariable
    {
    public:
      const string name;
      explicit UnknownVariable(string name_arg) : name{move(name_arg)}
      {
      }
    };

    //! Starts parsing a file, returns output in out
    void parse(const string &file_arg, const string &basename_arg, istream &modfile,
               ostream &output, bool debug, const vector<pair<string, string>> &defines,
               vector<filesystem::path> &paths_arg);

    //! Name of main file being parsed
    string file;

    //! Basename of main file being parsed
    string basename;

    //! Reference to the lexer
    unique_ptr<TokenizerFlex> lexer;

    //! Error handler
    void error(const Tokenizer::parser::location_type &location, const string &message) const;

    inline bool inContext() const { return !directive_stack.empty(); }

    inline void pushContext() { directive_stack.emplace(vector<DirectivePtr>()); }

    inline void pushContextTop(DirectivePtr statement) { directive_stack.top().emplace_back(move(statement)); }

    inline void pushStatements(DirectivePtr statement) { statements.emplace_back(move(statement)); }

    inline vector<DirectivePtr> popContext() { auto top = move(directive_stack.top()); directive_stack.pop(); return top; }
  };
}
#endif
