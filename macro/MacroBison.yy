/*
 * Copyright (C) 2008-2013 Dynare Team
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

/* Prologue:
   In Bison <= 2.3, it is inserted in both the .cc and .hh files.
   In Bison >= 2.3a, it is inserted only in the .cc file.
   Since Bison 2.4, the new %code directives provide a cleaner way of dealing
   with the prologue.
*/
%{
using namespace std;

#include "MacroValue.hh"

class MacroDriver;
%}

%name-prefix "Macro"

%parse-param { MacroDriver &driver }
%parse-param { ostream &out }
%lex-param { MacroDriver &driver }

%locations
%initial-action
{
  // Initialize the location filenames
  @$.begin.filename = @$.end.filename = &driver.file;
};

%debug
%error-verbose

%union
{
  string *string_val;
  int int_val;
  const MacroValue *mv;
};

%{
#include <cstdlib>  // Pour atoi()
#include "MacroDriver.hh"

/* this "connects" the bison parser in the driver to the flex scanner class
 * object. it defines the yylex() function call to pull the next token from the
 * current lexer object of the driver context. */
#undef yylex
#define yylex driver.lexer->lex

#define TYPERR_CATCH(statement, loc) try        \
    {                                           \
      statement;                                \
    }                                           \
  catch(MacroValue::TypeError &e)               \
    {                                           \
      driver.error(loc, e.message);             \
    }

%}

%token DEFINE LINE FOR IN IF ELSE ENDIF ECHO_DIR ERROR IFDEF IFNDEF
%token LPAREN RPAREN LBRACKET RBRACKET EQUAL EOL LENGTH

%token <int_val> INTEGER
%token <string_val> NAME STRING

%left COMMA
%left LOGICAL_OR
%left LOGICAL_AND
%left LESS GREATER LESS_EQUAL GREATER_EQUAL EQUAL_EQUAL EXCLAMATION_EQUAL
%nonassoc IN
%nonassoc COLON
%left PLUS MINUS
%left TIMES DIVIDE
%left UMINUS UPLUS EXCLAMATION
%left LBRACKET

%type <mv> expr array_expr
%%

%start statement_list_or_nothing;

statement_list_or_nothing : /* empty */
                          | statement_list
                          ;

statement_list : statement EOL
               | statement_list statement EOL
               ;

statement : expr
            { out << $1->toString(); }
          | DEFINE NAME EQUAL expr
            { driver.set_variable(*$2, $4); delete $2; }
          | FOR NAME IN expr
            { TYPERR_CATCH(driver.init_loop(*$2, $4), @$); delete $2; }
          | IF expr
            { TYPERR_CATCH(driver.begin_if($2), @$); }
          | IFDEF NAME
            { TYPERR_CATCH(driver.begin_ifdef(*$2), @$); delete $2; }
          | IFNDEF NAME
            { TYPERR_CATCH(driver.begin_ifndef(*$2), @$); delete $2; }
          | ECHO_DIR expr
            { TYPERR_CATCH(driver.echo(@$, $2), @$); }
          | ERROR expr
            { TYPERR_CATCH(driver.error(@$, $2), @$); }
          | LINE STRING INTEGER
            /* Ignore @#line declarations */
          ;

expr : INTEGER
       { $$ = new IntMV(driver, $1); }
     | STRING
       { $$ = new StringMV(driver, *$1); delete $1; }
     | NAME
       {
         try
           {
             $$ = driver.get_variable(*$1);
           }
         catch(MacroDriver::UnknownVariable(&e))
           {
             error(@$, "Unknown variable: " + e.name);
           }
         delete $1;
       }
     | LENGTH LPAREN array_expr RPAREN
       { TYPERR_CATCH($$ = $3->length(), @$); }
     | LPAREN expr RPAREN
       { $$ = $2; }
     | expr PLUS expr
       { TYPERR_CATCH($$ = *$1 + *$3, @$); }
     | expr MINUS expr
       { TYPERR_CATCH($$ = *$1 - *$3, @$); }
     | expr TIMES expr
       { TYPERR_CATCH($$ = *$1 * *$3, @$); }
     | expr DIVIDE expr
       { TYPERR_CATCH($$ = *$1 / *$3, @$); }
     | expr LESS expr
       { TYPERR_CATCH($$ = *$1 < *$3, @$); }
     | expr GREATER expr
       { TYPERR_CATCH($$ = *$1 > *$3, @$); }
     | expr LESS_EQUAL expr
       { TYPERR_CATCH($$ = *$1 <= *$3, @$); }
     | expr GREATER_EQUAL expr
       { TYPERR_CATCH($$ = *$1 >= *$3, @$); }
     | expr EQUAL_EQUAL expr
       { TYPERR_CATCH($$ = *$1 == *$3, @$); }
     | expr EXCLAMATION_EQUAL expr
       { TYPERR_CATCH($$ = *$1 != *$3, @$); }
     | expr LOGICAL_OR expr
       { TYPERR_CATCH($$ = *$1 || *$3, @$); }
     | expr LOGICAL_AND expr
       { TYPERR_CATCH($$ = *$1 && *$3, @$); }
     | MINUS expr %prec UMINUS
       { TYPERR_CATCH($$ = -*$2, @$); }
     | PLUS expr %prec UPLUS
       { TYPERR_CATCH($$ = +(*$2), @$); }
     | EXCLAMATION expr
       { TYPERR_CATCH($$ = !*$2, @$); }
     | expr LBRACKET array_expr RBRACKET
       {
         TYPERR_CATCH($$ = (*$1)[*$3], @$)
         catch(MacroValue::OutOfBoundsError)
           {
             error(@$, "Index out of bounds");
           }
       }
     | LBRACKET array_expr RBRACKET
       { $$ = $2; }
     | expr COLON expr
       { TYPERR_CATCH($$ = IntMV::new_range(driver, $1, $3), @$); }
     | expr IN expr
       { TYPERR_CATCH($$ = $1->in($3), @$); }
     ;

array_expr : expr
             { $$ = $1->toArray(); }
           | array_expr COMMA expr
             { TYPERR_CATCH($$ = $3->append($1), @$); }
           ;

%%

void
Macro::parser::error(const Macro::parser::location_type &l,
                     const string &m)
{
  driver.error(l, m);
}
