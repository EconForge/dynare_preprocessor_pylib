/*
 * Copyright (C) 2008-2018 Dynare Team
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
# error Impossible to include both ParsingDriver.hh and MacroDriver.hh
#endif

#include <string>
#include <iostream>
#include <stack>
#include <map>
#include <set>
#include <utility>

#include "MacroValue.hh"
#include "MacroBison.hh"

using namespace std;

// Declare MacroFlexLexer class
#ifndef __FLEX_LEXER_H
# define yyFlexLexer MacroFlexLexer
# include <FlexLexer.h>
# undef yyFlexLexer
#endif

//! The lexer class
/*! Actually it was necessary to subclass the MacroFlexLexer class generated by Flex,
  since the prototype for MacroFlexLexer::yylex() was not convenient.
*/
class MacroFlex : public MacroFlexLexer
{
private:
  //! Used to backup all the information related to a given scanning context
  class ScanContext
  {
  public:
    istream *input;
    struct yy_buffer_state *buffer;
    const Macro::parser::location_type yylloc;
    const bool is_for_context;
    const string for_body;
    const Macro::parser::location_type for_body_loc;
    const bool is_comprehension_context;
    const string comprehension_clause;
    const Macro::parser::location_type comprehension_clause_loc;
    const int comprehension_start_condition;
    const int comprehension_iter_nb;
    ScanContext(istream *input_arg, struct yy_buffer_state *buffer_arg,
                Macro::parser::location_type &yylloc_arg, bool is_for_context_arg,
                string for_body_arg,
                Macro::parser::location_type &for_body_loc_arg,
                bool is_comprehension_context_arg,
                string comprehension_clause_arg,
                Macro::parser::location_type &comprehension_clause_loc_arg,
                int comprehension_start_condition_arg,
                int comprehension_iter_nb_arg) :
      input(input_arg), buffer(buffer_arg), yylloc(yylloc_arg), is_for_context(is_for_context_arg),
      for_body(move(for_body_arg)), for_body_loc(for_body_loc_arg),
      is_comprehension_context{is_comprehension_context_arg},
      comprehension_clause{comprehension_clause_arg},
      comprehension_clause_loc{comprehension_clause_loc_arg},
      comprehension_start_condition{comprehension_start_condition_arg},
      comprehension_iter_nb{comprehension_iter_nb_arg}
    {
    }
  };

  //! The stack used to keep track of nested scanning contexts
  stack<ScanContext> context_stack;

  //! Input stream used for initialization of current scanning context
  /*! Kept for deletion at end of current scanning buffer */
  istream *input;

  //! Should we omit the @#line statements ?
  const bool no_line_macro;
  //! The paths to search when looking for .mod files
  vector<string> path;
  //! True iff current context is the body of a loop
  bool is_for_context{false};
  //! If current context is the body of a loop, contains the string of the loop body
  string for_body;
  //! If current context is the body of a loop, contains the location of the beginning of the body
  Macro::parser::location_type for_body_loc;

  //! Temporary variable used in FOR_BODY mode
  string for_body_tmp;
  //! Temporary variable used in FOR_BODY mode
  Macro::parser::location_type for_body_loc_tmp;
  //! Temporary variable used in FOR_BODY mode. Keeps track of the location of the @#for statement, for reporting messages
  Macro::parser::location_type for_stmt_loc_tmp;
  //! Temporary variable used in FOR_BODY mode. Keeps track of number of nested @#for/@#endfor
  int nested_for_nb;
  //! Set to true while parsing a FOR statement (only the statement, not the loop body)
  bool reading_for_statement{false};

  //! Temporary variable used in THEN_BODY and ELSE_BODY modes. Keeps track of number of nested @#if
  int nested_if_nb;
  //! Temporary variable used in THEN_BODY mode
  string then_body_tmp;
  //! Temporary variable used in THEN_BODY mode
  Macro::parser::location_type then_body_loc_tmp;
  //! Temporary variable used in THEN_BODY mode. Keeps track of the location of the @#if statement, for reporting messages
  Macro::parser::location_type if_stmt_loc_tmp;
  //! Temporary variable used in ELSE_BODY mode
  string else_body_tmp;
  //! Temporary variable used in ELSE_BODY mode
  Macro::parser::location_type else_body_loc_tmp;
  //! Set to true while parsing an IF statement (only the statement, not the body)
  bool reading_if_statement{false};

  //! Set to true while parsing the comprehension in a new buffer
  bool is_comprehension_context{false};
  //! Counter for the current comprehension
  int comprehension_iter_nb{0};
  //! The lexer start condition (EXPR or STMT) before entering the comprehension
  int comprehension_start_condition;
  //! Number of nested comprehensions, used during the construction of the new lexer buffer
  int nested_comprehension_nb{0};
  //! Temporary stores for the new lexer buffer
  string comprehension_clause, comprehension_clause_tmp;
  //! Stores for the location of the comprehension clause
  Macro::parser::location_type comprehension_clause_loc, comprehension_clause_loc_tmp;

  //! Output the @#line declaration
  void output_line(Macro::parser::location_type *yylloc) const;

  //! Save current scanning context
  void save_context(Macro::parser::location_type *yylloc);

  //! Restore last scanning context
  void restore_context(Macro::parser::location_type *yylloc);

  //! pushes the colon-separated paths passed to @#includepath onto the path vector
  void push_path(string *includepath, Macro::parser::location_type *yylloc,
                 MacroDriver &driver);

  //! Saves current scanning context and create a new context with content of filename
  /*! Filename must be a newly allocated string which will be deleted by the lexer */
  void create_include_context(string *filename, Macro::parser::location_type *yylloc,
                              MacroDriver &driver);

  //! Saves current scanning context and create a new context based on the "then" body
  void create_then_context(Macro::parser::location_type *yylloc);

  //! Saves current scanning context and create a new context based on the "else" body
  void create_else_context(Macro::parser::location_type *yylloc);

  //! Initialise a new flex buffer with the loop body
  void new_loop_body_buffer(Macro::parser::location_type *yylloc);

  //! Initialize a new flex buffer with the comprehension conditional clause
  void new_comprehension_clause_buffer(Macro::parser::location_type *yylloc);

public:
  MacroFlex(istream *in, ostream *out, bool no_line_macro_arg, vector<string> path_arg);

  //! The main lexing function
  Macro::parser::token_type lex(Macro::parser::semantic_type *yylval,
                                Macro::parser::location_type *yylloc,
                                MacroDriver &driver);
};

//! Implements the macro expansion using a Flex scanner and a Bison parser
class MacroDriver
{
  friend class MacroValue;
private:
  // The map defining an environment
  using env_t = map<string, MacroValuePtr>;

  //! Environment: maps macro variables to their values
  env_t env;

  //! Environment for function currently being evaluated
  vector<env_t> func_env;

  //! Stack used to keep track of (possibly nested) loops
  //! First element is loop variable name(s)
  //! Second is the array over which iteration is done
  //! Third is subscript to be used by next call of iter_loop() (beginning with 0) */
  stack<tuple<vector<string>, shared_ptr<ArrayMV>, int>> loop_stack;

  stack<tuple<vector<string>, shared_ptr<ArrayMV>, int>> comprehension_stack;
public:
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
  /*! \param no_line_macro should we omit the @#line statements ? */
  void parse(const string &f, const string &fb, const string &modfiletxt, ostream &out, bool debug, bool no_line_macro_arg,
             map<string, string> defines, vector<string> path);

  //! Name of main file being parsed
  string file;

  //! Basename of main file being parsed
  string basename;

  //! Whether or not to print @#line
  bool no_line_macro;

  //! Reference to the lexer
  unique_ptr<MacroFlex> lexer;

  //! Used to store the value of the last @#if condition
  bool last_if;

  //! Error handler
  void error(const Macro::parser::location_type &l, const string &m) const;

  //! Print variables
  string printvars(const Macro::parser::location_type &l, const bool save) const;

  //! Set a variable
  void set_variable(const string &name, MacroValuePtr value);

  //! Replace "@{x}" with the value of x (if it exists in the environment) in a string
  //! Check for variable existence first in func_env before checking in env
  string replace_vars_in_str(const string &s) const;

  //! Set a function with arguments
  void set_string_function(const string &name, vector<string> &args, const MacroValuePtr &value);

  //! Push function arguments onto func_env stack setting equal to NULL
  void push_args_into_func_env(const vector<string> &args);

  //! Remove last entry in func_env vector
  void pop_func_env();

  //! Evaluate a function with arguments
  MacroValuePtr eval_string_function(const string &name, const vector<MacroValuePtr> &args);

  //! Get a variable
  MacroValuePtr get_variable(const string &name) const noexcept(false);

  //! Initiate a for loop
  /*! Does not set name = value[1]. You must call iter_loop() for that. */
  void init_loop(const string &name, MacroValuePtr value) noexcept(false);
  /*! Same as above but for looping over tuple array */
  void init_loop(const vector<string> &names, MacroValuePtr value) noexcept(false);

  //! Iterate innermost loop
  /*! Returns false if iteration is no more possible (end of loop);
      in that case it destroys the pointer given to init_loop() */
  bool iter_loop() noexcept(false);

  //! Initializes the evaluation of a comprehension
  void init_comprehension(const vector<string> &names, MacroValuePtr value);
  //! Iterates during the evaluation of the comprehension
  void iter_comprehension();
  //! Helper to construct the value corresponding to the comprehension
  void possibly_add_comprehension_element(vector<MacroValuePtr> &v, MacroValuePtr test_expr) const;
  //! Returns the size of the set over which the current comprehension iterates
  int get_comprehension_iter_nb() const;

  //! Begins an @#if statement
  void begin_if(const MacroValuePtr &value) noexcept(false);

  //! Begins an @#ifdef statement
  void begin_ifdef(const string &name);

  //! Begins an @#ifndef statement
  void begin_ifndef(const string &name);

  //! Executes @#echo directive
  void echo(const Macro::parser::location_type &l, const MacroValuePtr &value) const noexcept(false);

  //! Executes @#error directive
  void error(const Macro::parser::location_type &l, const MacroValuePtr &value) const noexcept(false);
};

#endif // ! MACRO_DRIVER_HH
