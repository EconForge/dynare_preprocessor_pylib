// -*- C++ -*-
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

%language "c++"
%require "3.0"
%defines
%define api.value.type variant
%define api.namespace {Macro}
%define parse.assert
%define parse.error verbose
%define parse.trace

%code top {
class MacroDriver;
}

%parse-param { MacroDriver &driver }
%parse-param { ostream &out }
%lex-param { MacroDriver &driver }

%locations
%initial-action
{
  // Initialize the location filenames
  @$.begin.filename = @$.end.filename = &driver.file;
};

%code requires {
#include "MacroValue.hh"
}

%code {
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

}

%token COMMA DEFINE LINE FOR IN IF ELSE ENDIF ECHO_DIR ERROR IFDEF IFNDEF
%token LPAREN RPAREN LBRACKET RBRACKET EQUAL EOL LENGTH ECHOMACROVARS SAVE

%token <int> INTEGER
%token <string> NAME STRING

%left LOGICAL_OR
%left LOGICAL_AND
%left LESS GREATER LESS_EQUAL GREATER_EQUAL EQUAL_EQUAL EXCLAMATION_EQUAL
%nonassoc IN
%nonassoc COLON
%left PLUS MINUS
%left TIMES DIVIDE
%precedence UMINUS UPLUS EXCLAMATION
%precedence LBRACKET

%type <vector<string>> func_args
%type <MacroValuePtr> expr
%type <vector<MacroValuePtr>> comma_expr
%%

%start statement_list_or_nothing;

statement_list_or_nothing : %empty
                          | statement_list
                          ;

statement_list : statement EOL
               | statement_list statement EOL
               ;

statement : expr
            { out << $1->toString(); }
          | DEFINE NAME EQUAL expr
            { driver.set_variable($2, $4); }
          | FOR NAME IN expr
            { TYPERR_CATCH(driver.init_loop($2, $4), @$); }
          | IF expr
            { TYPERR_CATCH(driver.begin_if($2), @$); }
          | IFDEF NAME
            { TYPERR_CATCH(driver.begin_ifdef($2), @$); }
          | IFNDEF NAME
            { TYPERR_CATCH(driver.begin_ifndef($2), @$); }
          | ECHO_DIR expr
            { TYPERR_CATCH(driver.echo(@$, $2), @$); }
          | ERROR expr
            { TYPERR_CATCH(driver.error(@$, $2), @$); }
          | LINE STRING INTEGER
            /* Ignore @#line declarations */
          | ECHOMACROVARS
            { driver.printvars(@$, true); }
          | ECHOMACROVARS LPAREN SAVE RPAREN
            { out << driver.printvars(@$, false); }
          | DEFINE NAME LPAREN func_args { driver.push_args_into_func_env($4); } RPAREN EQUAL expr
            {
              TYPERR_CATCH(driver.set_string_function($2, $4, $8), @$);
              driver.pop_func_env();
            }
          ;

func_args : NAME
            { $$ = vector<string>{$1}; }
          | func_args COMMA NAME
            { $1.push_back($3); $$ = $1; }
          ;

expr : INTEGER
       { $$ = make_shared<IntMV>($1); }
     | STRING
       { $$ = make_shared<StringMV>(driver.replace_vars_in_str($1)); }
     | NAME
       {
         try
           {
             $$ = driver.get_variable($1);
           }
         catch(MacroDriver::UnknownVariable(&e))
           {
             error(@$, "Unknown variable: " + e.name);
           }
       }
     | NAME LPAREN comma_expr RPAREN
       { TYPERR_CATCH($$ = driver.eval_string_function($1, $3), @$); }
     | LENGTH LPAREN expr RPAREN
       { TYPERR_CATCH($$ = $3->length(), @$); }
     | LPAREN expr RPAREN
       { $$ = $2; }
     | expr PLUS expr
       { TYPERR_CATCH($$ = $1->plus($3), @$); }
     | expr MINUS expr
       { TYPERR_CATCH($$ = $1->minus($3), @$); }
     | expr TIMES expr
       { TYPERR_CATCH($$ = $1->times($3), @$); }
     | expr DIVIDE expr
       {
         TYPERR_CATCH($$ = $1->divide($3), @$)
         catch(MacroValue::DivisionByZeroError)
           {
             error(@$, "Division by zero");
           }
       }
     | expr LESS expr
       { TYPERR_CATCH($$ = $1->is_less($3), @$); }
     | expr GREATER expr
       { TYPERR_CATCH($$ = $1->is_greater($3), @$); }
     | expr LESS_EQUAL expr
       { TYPERR_CATCH($$ = $1->is_less_equal($3), @$); }
     | expr GREATER_EQUAL expr
       { TYPERR_CATCH($$ = $1->is_greater_equal($3), @$); }
     | expr EQUAL_EQUAL expr
       { $$ = $1->is_equal($3); }
     | expr EXCLAMATION_EQUAL expr
       { $$ = $1->is_different($3); }
     | expr LOGICAL_OR expr
       { TYPERR_CATCH($$ = $1->logical_or($3), @$); }
     | expr LOGICAL_AND expr
       { TYPERR_CATCH($$ = $1->logical_and($3), @$); }
     | MINUS expr %prec UMINUS
       { TYPERR_CATCH($$ = $2->unary_minus(), @$); }
     | PLUS expr %prec UPLUS
       { TYPERR_CATCH($$ = $2->unary_plus(), @$); }
     | EXCLAMATION expr
       { TYPERR_CATCH($$ = $2->logical_not(), @$); }
     | expr LBRACKET expr RBRACKET
       {
         TYPERR_CATCH($$ = $1->subscript($3), @$)
         catch(MacroValue::OutOfBoundsError)
           {
             error(@$, "Index out of bounds");
           }
       }
     | LBRACKET comma_expr RBRACKET
       { $$ = make_shared<ArrayMV>($2); }
     | expr COLON expr
       { TYPERR_CATCH($$ = ArrayMV::range($1, $3), @$); }
     | expr IN expr
       { TYPERR_CATCH($$ = $3->in($1), @$); }
     ;

comma_expr : %empty
             { $$ = vector<MacroValuePtr>{}; } // Empty array
           | expr
             { $$ = vector<MacroValuePtr>{$1}; }
           | comma_expr COMMA expr
             { $1.push_back($3); $$ = $1; }
           ;

%%

void
Macro::parser::error(const Macro::parser::location_type &l,
                     const string &m)
{
  driver.error(l, m);
}
