/*
 * Copyright (C) 2008 Dynare Team
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

%skeleton "lalr1.cc"
%require "2.3"
%defines

%{
using namespace std;

class MacroDriver;
%}

%name-prefix="Macro"

%parse-param { MacroDriver &driver }
%parse-param { ostream &out }
%lex-param { MacroDriver &driver }

%locations
%initial-action
{
  // Initialize the location filenames
  @$.begin.filename = @$.end.filename = &driver.file;
  // Output first @line statement
  out << "@line \"" << driver.file << "\" 1" << endl;
};

%debug
%error-verbose

%union
{
  string *string_val;
  int int_val;
};

%{
#include "MacroDriver.hh"

/* this "connects" the bison parser in the driver to the flex scanner class
 * object. it defines the yylex() function call to pull the next token from the
 * current lexer object of the driver context. */
#undef yylex
#define yylex driver.lexer->lex
%}

%token <int_val> INT_NUMBER
%token <string_val> NAME
%%

%start statement_list_or_nothing;

statement_list_or_nothing : /* empty */;

%%

void
Macro::parser::error(const Macro::parser::location_type &l,
                     const string &m)
{
  driver.error(l, m);
}

/*
  Local variables:
  mode: C++
  End:
*/
