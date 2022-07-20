/*
 * Copyright © 2019-2022 Dynare Team
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
 * along with Dynare.  If not, see <https://www.gnu.org/licenses/>.
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
    TokenizerFlex &operator=(const TokenizerFlex &) = delete;

    //! The main lexing function
    Tokenizer::parser::token_type lex(Tokenizer::parser::semantic_type *yylval,
                                      Tokenizer::parser::location_type *yylloc,
                                      macro::Driver &driver);
    static void location_increment(Tokenizer::parser::location_type *yylloc, const char *yytext);
  };

  //! Implements the macro expansion using a Flex scanner and a Bison parser
  class Driver
  {
  private:
    vector<DirectivePtr> statements;
    stack<vector<DirectivePtr>> directive_stack;
  public:
    Driver() = default;
    Driver(const Driver &) = delete;
    Driver &operator=(const Driver &) = delete;

    //! Exception thrown when value of an unknown variable is requested
    struct UnknownVariable
    {
      const string name;
    };

    //! Starts parsing a file, modifies `env`, `paths` and `output`
    //! as they are modified by various macro directives
    void parse(const string &file, const string &basename, const istream &modfile,
               bool debug, const vector<pair<string, string>> &defines,
               Environment &env, vector<filesystem::path> &paths, ostream &output);

    //! Name of main file being parsed
    string file;

    //! Basename of main file being parsed
    string basename;

    //! Reference to the lexer
    unique_ptr<TokenizerFlex> lexer;

    //! Error handler
    void error(const Tokenizer::parser::location_type &location, const string &message) const;

    bool
    inContext() const
    {
      return !directive_stack.empty();
    }

    void
    pushContext()
    {
      directive_stack.emplace(vector<DirectivePtr>());
    }

    void
    pushContextTop(DirectivePtr statement)
    {
      directive_stack.top().emplace_back(move(statement));
    }

    void
    pushStatements(DirectivePtr statement)
    {
      statements.emplace_back(move(statement));
    }

    vector<DirectivePtr>
    popContext()
    {
      auto top = move(directive_stack.top());
      directive_stack.pop();
      return top;
    }
  };
}
#endif
