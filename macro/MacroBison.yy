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

#include "MacroValue.hh"

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
  MacroValue *mv;
};

%{
#include <stdlib.h>  // Pour atoi()
#include "MacroDriver.hh"

/* this "connects" the bison parser in the driver to the flex scanner class
 * object. it defines the yylex() function call to pull the next token from the
 * current lexer object of the driver context. */
#undef yylex
#define yylex driver.lexer->lex
%}

%token DEFINE
%token LPAREN RPAREN LBRACKET RBRACKET COLON EQUAL EOL

%token <int_val> INTEGER
%token <string_val> NAME STRING

%left LOGICAL_OR
%left LOGICAL_AND
%left LESS GREATER LESS_EQUAL GREATER_EQUAL EQUAL_EQUAL EXCLAMATION_EQUAL
%left TIMES DIVIDE
%left PLUS MINUS
%left UMINUS UPLUS EXCLAMATION
%left LBRACKET

%type <mv> expr
%%

%start statement_list_or_nothing;

statement_list_or_nothing : /* empty */
                          | statement_list;

statement_list : statement EOL
               | statement_list statement EOL;

statement : expr
            { *driver.out_stream << $1->toString(); delete $1; }
          | DEFINE NAME EQUAL expr
            { driver.env[*$2] = $4; delete $2; }

expr : INTEGER
       { $$ = new IntMV($1); }
     | STRING
       { $$ = new StringMV(*$1); delete $1; }
     | NAME
       { $$ = driver.env[*$1]->clone(); delete $1; }
     | LPAREN expr RPAREN
       { $$ = $2; }
     | expr PLUS expr
       { $$ = *$1 + *$3; delete $1; delete $3; }
     | expr MINUS expr
       { $$ = *$1 - *$3; delete $1; delete $3; }
     | expr TIMES expr
       { $$ = *$1 * *$3; delete $1; delete $3; }
     | expr DIVIDE expr
       { $$ = *$1 / *$3; delete $1; delete $3; }
     | expr LESS expr
       { $$ = *$1 < *$3; delete $1; delete $3; }
     | expr GREATER expr
       { $$ = *$1 > *$3; delete $1; delete $3; }
     | expr LESS_EQUAL expr
       { $$ = *$1 <= *$3; delete $1; delete $3; }
     | expr GREATER_EQUAL expr
       { $$ = *$1 >= *$3; delete $1; delete $3; }
     | expr EQUAL_EQUAL expr
       { $$ = *$1 == *$3; delete $1; delete $3; }
     | expr EXCLAMATION_EQUAL expr
       { $$ = *$1 != *$3; delete $1; delete $3; }
     | expr LOGICAL_OR expr
       { $$ = *$1 || *$3; delete $1; delete $3; }
     | expr LOGICAL_AND expr
       { $$ = *$1 && *$3; delete $1; delete $3; }
     | MINUS expr %prec UMINUS
       { $$ = -*$2; delete $2;}
     | PLUS expr %prec UPLUS
       { $$ = $2; }
     | EXCLAMATION expr
       { $$ = !*$2; delete $2; }
     | expr LBRACKET expr RBRACKET
       { $$ = (*$1)[*$3]; delete $1; delete $3; }
     ;

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
