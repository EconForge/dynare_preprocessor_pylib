// -*- C++ -*-
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

%language "c++"
%require "3.0"
%defines
%define api.value.type variant
%define api.namespace {Tokenizer}
%define parse.assert
%define parse.error verbose
%define parse.trace

%code requires {
namespace macro { class Driver; }
}

%param { macro::Driver &driver }

%locations
%initial-action
{
  // Initialize the location filenames
  @$.begin.filename = @$.end.filename = &driver.file;
};

%code requires {
#include "Directives.hh"
using namespace macro;
}

%code {
#include "Driver.hh"

/* this "connects" the bison parser in the driver to the flex scanner class
 * object. it defines the yylex() function call to pull the next token from the
 * current lexer object of the driver context. */
#undef yylex
#define yylex driver.lexer->lex

}

%token FOR ENDFOR IF IFDEF IFNDEF ELSE ENDIF TRUE FALSE
%token INCLUDE INCLUDEPATH DEFINE EQUAL D_ECHO ERROR
%token COMMA LPAREN RPAREN LBRACKET RBRACKET WHEN
%token BEGIN_EVAL END_EVAL ECHOMACROVARS SAVE

%token EXP LOG LN LOG10 SIN COS TAN ASIN ACOS ATAN
%token SQRT CBRT SIGN MAX MIN FLOOR CEIL TRUNC SUM MOD
%token ERF ERFC GAMMA LGAMMA ROUND NORMPDF NORMCDF LENGTH

%token INT

%left OR
%left AND
%left EQUAL_EQUAL NOT_EQUAL
%left LESS GREATER LESS_EQUAL GREATER_EQUAL
%nonassoc IN
%left COLON
%left UNION
%left INTERSECTION
%left PLUS MINUS
%left TIMES DIVIDE
%precedence UMINUS UPLUS NOT
%precedence CAST_INT
%nonassoc POWER

%token <string> NAME TEXT QUOTED_STRING NUMBER EOL

%type <DirectivePtr> statement
%type <DirectivePtr> directive directive_one_line directive_multiline for if ifdef ifndef text
%type <EvalPtr> eval
%type <ExpressionPtr> expr
%type <FunctionPtr> function
%type <VariablePtr> symbol

%type <vector<VariablePtr>> comma_name
%type <vector<ExpressionPtr>> comma_expr function_args tuple_comma_expr colon_expr

%%

%start statements;

statements : statement
             {
               driver.inContext() ? driver.pushContextTop($1) : driver.pushStatements($1);
             }
           | statements statement
             {
               driver.inContext() ? driver.pushContextTop($2) : driver.pushStatements($2);
             }
           ;

statement : directive
            { $$ = $1; }
          | text
            { $$ = $1; }
          | eval
            { $$ = $1; }
          ;

directive : directive_one_line EOL
            { $$ = $1; }
          | directive_multiline EOL
            { $$ = $1; }
          ;

directive_one_line : INCLUDE expr
                     { $$ = make_shared<Include>($2, driver.env, @$); }
                   | INCLUDEPATH expr
                     { $$ = make_shared<IncludePath>($2, driver.env, @$); }
                   | DEFINE symbol EQUAL expr
                     { $$ = make_shared<Define>($2, $4, driver.env, @$); }
                   | DEFINE function EQUAL expr
                     { $$ = make_shared<Define>($2, $4, driver.env, @$); }
                   | D_ECHO expr
                     { $$ = make_shared<Echo>($2, driver.env, @$); }
                   | ERROR expr
                     { $$ = make_shared<Error>($2, driver.env, @$); }
                   | ECHOMACROVARS
                     { $$ = make_shared<EchoMacroVars>(false, driver.env, @$); }
                   | ECHOMACROVARS LPAREN SAVE RPAREN
                     { $$ = make_shared<EchoMacroVars>(true, driver.env, @$); }
                   ;

directive_multiline : for
                      { $$ = $1; }
                    | if
                      { $$ = $1; }
                    | ifdef
                      { $$ = $1; }
                    | ifndef
                      { $$ = $1; }
                    ;

for : FOR { driver.pushContext(); } NAME IN expr EOL statements ENDFOR
      {
        vector<VariablePtr> vvnp = {make_shared<Variable>($3, driver.env, @3)};
        auto vdp = driver.popContext();
        vdp.emplace_back(make_shared<TextNode>("\n", driver.env, @8));
        $$ = make_shared<For>(vvnp, $5, vdp, driver.env, @$);
      }
    ;

for : FOR { driver.pushContext(); } LPAREN comma_name RPAREN IN expr EOL statements ENDFOR
      {
        vector<VariablePtr> vvnp;
        for (auto & it : $4)
          {
            auto vnp = dynamic_pointer_cast<Variable>(it);
            if (!vnp)
              error(@$, "For loop indices must be variables");
            vvnp.push_back(vnp);
          }
        auto vdp = driver.popContext();
        vdp.emplace_back(make_shared<TextNode>("\n", driver.env, @10));
        $$ = make_shared<For>(vvnp, $7, vdp, driver.env, @$);
      }
    ;

comma_name : NAME
             { $$ = vector<VariablePtr>{make_shared<Variable>($1, driver.env, @$)}; }
           | comma_name COMMA NAME
             { $1.emplace_back(make_shared<Variable>($3, driver.env, @3)); $$ = $1; }
           ;

if_begin : IF { driver.pushContext(); }
         ;

if : if_begin expr EOL statements ENDIF
     {
       auto ifContext = driver.popContext();
       ifContext.emplace_back(make_shared<TextNode>("\n", driver.env, @5));
       $$ = make_shared<If>($2, ifContext, driver.env, @$);
     }
   | if_begin expr EOL statements ELSE EOL { driver.pushContext(); } statements ENDIF
     {
       auto elseContext = driver.popContext();
       elseContext.emplace_back(make_shared<TextNode>("\n", driver.env, @9));
       auto ifContext = driver.popContext();
       ifContext.emplace_back(make_shared<TextNode>("\n", driver.env, @5));
       $$ = make_shared<If>($2, ifContext, elseContext, driver.env, @$);
     }
   ;

ifdef_begin : IFDEF { driver.pushContext(); }
            ;

ifdef : ifdef_begin expr EOL statements ENDIF
        {
          auto ifContext = driver.popContext();
          ifContext.emplace_back(make_shared<TextNode>("\n", driver.env, @5));
          $$ = make_shared<Ifdef>($2, ifContext, driver.env, @$);
        }
      | ifdef_begin expr EOL statements ELSE EOL { driver.pushContext(); } statements ENDIF
        {
          auto elseContext = driver.popContext();
          elseContext.emplace_back(make_shared<TextNode>("\n", driver.env, @9));
          auto ifContext = driver.popContext();
          ifContext.emplace_back(make_shared<TextNode>("\n", driver.env, @5));
          $$ = make_shared<Ifdef>($2, ifContext, elseContext, driver.env, @$);
        }
      ;

ifndef_begin : IFNDEF { driver.pushContext(); }
             ;

ifndef : ifndef_begin expr EOL statements ENDIF
        {
          auto ifContext = driver.popContext();
          ifContext.emplace_back(make_shared<TextNode>("\n", driver.env, @5));
          $$ = make_shared<Ifndef>($2, ifContext, driver.env, @$);
        }
       | ifndef_begin expr EOL statements ELSE EOL { driver.pushContext(); } statements ENDIF
         {
           auto elseContext = driver.popContext();
           elseContext.emplace_back(make_shared<TextNode>("\n", driver.env, @9));
           auto ifContext = driver.popContext();
           ifContext.emplace_back(make_shared<TextNode>("\n", driver.env, @5));
           $$ = make_shared<Ifndef>($2, ifContext, elseContext, driver.env, @$);
         }
       ;

text : TEXT
       { $$ = make_shared<TextNode>($1, driver.env, @$); }
     | EOL
       { $$ = make_shared<TextNode>($1, driver.env, @$); }
     ;

eval : BEGIN_EVAL expr END_EVAL
       { $$ = make_shared<Eval>($2, driver.env, @$); }
     ;

symbol : NAME
         { $$ = make_shared<Variable>($1, driver.env, @$); }
       ;

function : NAME LPAREN RPAREN
           { $$ = make_shared<Function>($1, vector<ExpressionPtr>(), driver.env, @$); }
         | NAME LPAREN function_args RPAREN
           { $$ = make_shared<Function>($1, $3, driver.env, @$); }
         ;

function_args : symbol
                { $$ = vector<ExpressionPtr>{$1}; }
              | function_args COMMA symbol
                { $1.emplace_back($3); $$ = $1; }
              ;

comma_expr : %empty
             { $$ = vector<ExpressionPtr>{}; } // Empty array
           | expr
             { $$ = vector<ExpressionPtr>{$1}; }
           | comma_expr COMMA expr
             { $1.emplace_back($3); $$ = $1; }
           ;

tuple_comma_expr : %empty
                   { $$ = vector<ExpressionPtr>{}; } // Empty tuple
                 | expr COMMA
                   { $$ = vector<ExpressionPtr>{$1}; }
                 | expr COMMA expr
                   { $$ = vector<ExpressionPtr>{$1, $3}; }
                 | tuple_comma_expr COMMA expr
                   { $1.emplace_back($3); $$ = $1; }
                 ;

colon_expr : expr COLON expr
             { $$ = vector<ExpressionPtr>{$1, $3}; }
           | colon_expr COLON expr
             { $1.emplace_back($3); $$ = $1; }
           ;

expr : LPAREN expr RPAREN
       { $$ = $2; }
     | symbol
       { $$ = $1; }
     | NAME LPAREN comma_expr RPAREN
       { $$ = make_shared<Function>($1, $3, driver.env, @$); }
     | TRUE
       { $$ = make_shared<Bool>(true, driver.env, @$); }
     | FALSE
       { $$ = make_shared<Bool>(false, driver.env, @$); }
     | NUMBER
       { $$ = make_shared<Double>($1, driver.env, @$); }
     | QUOTED_STRING
       { $$ = make_shared<String>($1, driver.env, @$); }
     | colon_expr
       {
         if ($1.size() == 2)
           $$ = make_shared<Array>($1[0], $1[1], driver.env, @$);
         else if ($1.size() == 3)
           $$ = make_shared<Array>($1[0], $1[1], $1[2], driver.env, @$);
         else
           error(@$, "The colon operator only works with 2 or 3 arguments");
       }
     | LBRACKET comma_expr RBRACKET
       { $$ = make_shared<Array>($2, driver.env, @$); }
     | symbol LBRACKET comma_expr RBRACKET
       { $1->addIndexing($3); $$ = $1; }
     | LPAREN tuple_comma_expr RPAREN
       { $$ = make_shared<Tuple>($2, driver.env, @$); }
     | LBRACKET expr IN expr WHEN expr RBRACKET
       { $$ = make_shared<Comprehension>(true, $2, $4, $6, driver.env, @$); }
     | LBRACKET expr FOR expr IN expr RBRACKET
       { $$ = make_shared<Comprehension>($2, $4, $6, driver.env, @$); }
     | LBRACKET expr FOR expr IN expr WHEN expr RBRACKET
       { $$ = make_shared<Comprehension>($2, $4, $6, $8, driver.env, @$); }
     | LPAREN INT RPAREN expr %prec CAST_INT
       { $$ = make_shared<UnaryOp>(codes::UnaryOp::cast_int, $4, driver.env, @$); }
     | NOT expr
       { $$ = make_shared<UnaryOp>(codes::UnaryOp::logical_not, $2, driver.env, @$); }
     | MINUS expr %prec UMINUS
       { $$ = make_shared<UnaryOp>(codes::UnaryOp::unary_minus, $2, driver.env, @$); }
     | PLUS expr %prec UPLUS
       { $$ = make_shared<UnaryOp>(codes::UnaryOp::unary_plus, $2, driver.env, @$); }
     | LENGTH LPAREN expr RPAREN
       { $$ = make_shared<UnaryOp>(codes::UnaryOp::length, $3, driver.env, @$); }
     | EXP LPAREN expr RPAREN
       { $$ = make_shared<UnaryOp>(codes::UnaryOp::exp, $3, driver.env, @$); }
     | LOG LPAREN expr RPAREN
       { $$ = make_shared<UnaryOp>(codes::UnaryOp::ln, $3, driver.env, @$); }
     | LN LPAREN expr RPAREN
       { $$ = make_shared<UnaryOp>(codes::UnaryOp::ln, $3, driver.env, @$); }
     | LOG10 LPAREN expr RPAREN
       { $$ = make_shared<UnaryOp>(codes::UnaryOp::log10, $3, driver.env, @$); }
     | SIN LPAREN expr RPAREN
       { $$ = make_shared<UnaryOp>(codes::UnaryOp::sin, $3, driver.env, @$); }
     | COS LPAREN expr RPAREN
       { $$ = make_shared<UnaryOp>(codes::UnaryOp::cos, $3, driver.env, @$); }
     | TAN LPAREN expr RPAREN
       { $$ = make_shared<UnaryOp>(codes::UnaryOp::tan, $3, driver.env, @$); }
     | ASIN LPAREN expr RPAREN
       { $$ = make_shared<UnaryOp>(codes::UnaryOp::asin, $3, driver.env, @$); }
     | ACOS LPAREN expr RPAREN
       { $$ = make_shared<UnaryOp>(codes::UnaryOp::acos, $3, driver.env, @$); }
     | ATAN LPAREN expr RPAREN
       { $$ = make_shared<UnaryOp>(codes::UnaryOp::atan, $3, driver.env, @$); }
     | SQRT LPAREN expr RPAREN
       { $$ = make_shared<UnaryOp>(codes::UnaryOp::sqrt, $3, driver.env, @$); }
     | CBRT LPAREN expr RPAREN
       { $$ = make_shared<UnaryOp>(codes::UnaryOp::cbrt, $3, driver.env, @$); }
     | SIGN LPAREN expr RPAREN
       { $$ = make_shared<UnaryOp>(codes::UnaryOp::sign, $3, driver.env, @$); }
     | FLOOR LPAREN expr RPAREN
       { $$ = make_shared<UnaryOp>(codes::UnaryOp::floor, $3, driver.env, @$); }
     | CEIL LPAREN expr RPAREN
       { $$ = make_shared<UnaryOp>(codes::UnaryOp::ceil, $3, driver.env, @$); }
     | TRUNC LPAREN expr RPAREN
       { $$ = make_shared<UnaryOp>(codes::UnaryOp::trunc, $3, driver.env, @$); }
     | SUM LPAREN expr RPAREN
       { $$ = make_shared<UnaryOp>(codes::UnaryOp::sum, $3, driver.env, @$); }
     | ERF LPAREN expr RPAREN
       { $$ = make_shared<UnaryOp>(codes::UnaryOp::erf, $3, driver.env, @$); }
     | ERFC LPAREN expr RPAREN
       { $$ = make_shared<UnaryOp>(codes::UnaryOp::erfc, $3, driver.env, @$); }
     | GAMMA LPAREN expr RPAREN
       { $$ = make_shared<UnaryOp>(codes::UnaryOp::gamma, $3, driver.env, @$); }
     | LGAMMA LPAREN expr RPAREN
       { $$ = make_shared<UnaryOp>(codes::UnaryOp::lgamma, $3, driver.env, @$); }
     | ROUND LPAREN expr RPAREN
       { $$ = make_shared<UnaryOp>(codes::UnaryOp::round, $3, driver.env, @$); }
     | NORMPDF LPAREN expr RPAREN
       { $$ = make_shared<UnaryOp>(codes::UnaryOp::normpdf, $3, driver.env, @$); }
     | NORMCDF LPAREN expr RPAREN
       { $$ = make_shared<UnaryOp>(codes::UnaryOp::normcdf, $3, driver.env, @$); }
     | expr PLUS expr
       { $$ = make_shared<BinaryOp>(codes::BinaryOp::plus, $1, $3, driver.env, @$); }
     | expr MINUS expr
       { $$ = make_shared<BinaryOp>(codes::BinaryOp::minus, $1, $3, driver.env, @$); }
     | expr TIMES expr
       { $$ = make_shared<BinaryOp>(codes::BinaryOp::times, $1, $3, driver.env, @$); }
     | expr DIVIDE expr
       { $$ = make_shared<BinaryOp>(codes::BinaryOp::divide, $1, $3, driver.env, @$); }
     | expr POWER expr
       { $$ = make_shared<BinaryOp>(codes::BinaryOp::power, $1, $3, driver.env, @$); }
     | expr EQUAL_EQUAL expr
       { $$ = make_shared<BinaryOp>(codes::BinaryOp::equal_equal, $1, $3, driver.env, @$); }
     | expr NOT_EQUAL expr
       { $$ = make_shared<BinaryOp>(codes::BinaryOp::not_equal, $1, $3, driver.env, @$); }
     | expr LESS expr
       { $$ = make_shared<BinaryOp>(codes::BinaryOp::less, $1, $3, driver.env, @$); }
     | expr GREATER expr
       { $$ = make_shared<BinaryOp>(codes::BinaryOp::greater, $1, $3, driver.env, @$); }
     | expr LESS_EQUAL expr
       { $$ = make_shared<BinaryOp>(codes::BinaryOp::less_equal, $1, $3, driver.env, @$); }
     | expr GREATER_EQUAL expr
       { $$ = make_shared<BinaryOp>(codes::BinaryOp::greater_equal, $1, $3, driver.env, @$); }
     | expr AND expr
       { $$ = make_shared<BinaryOp>(codes::BinaryOp::logical_and, $1, $3, driver.env, @$); }
     | expr OR expr
       { $$ = make_shared<BinaryOp>(codes::BinaryOp::logical_or, $1, $3, driver.env, @$); }
     | expr IN expr
       { $$ = make_shared<BinaryOp>(codes::BinaryOp::in, $1, $3, driver.env, @$); }
     | expr UNION expr
       { $$ = make_shared<BinaryOp>(codes::BinaryOp::set_union, $1, $3, driver.env, @$); }
     | expr INTERSECTION expr
       { $$ = make_shared<BinaryOp>(codes::BinaryOp::set_intersection, $1, $3, driver.env, @$); }
     | MAX LPAREN expr COMMA expr RPAREN
       { $$ = make_shared<BinaryOp>(codes::BinaryOp::max, $3, $5, driver.env, @$); }
     | MIN LPAREN expr COMMA expr RPAREN
       { $$ = make_shared<BinaryOp>(codes::BinaryOp::min, $3, $5, driver.env, @$); }
     | MOD LPAREN expr COMMA expr RPAREN
       { $$ = make_shared<BinaryOp>(codes::BinaryOp::mod, $3, $5, driver.env, @$); }
     | NORMPDF LPAREN expr COMMA expr COMMA expr RPAREN
       { $$ = make_shared<TrinaryOp>(codes::TrinaryOp::normpdf, $3, $5, $7, driver.env, @$); }
     | NORMCDF LPAREN expr COMMA expr COMMA expr RPAREN
       { $$ = make_shared<TrinaryOp>(codes::TrinaryOp::normcdf, $3, $5, $7, driver.env, @$); }
     ;

%%

void
Tokenizer::parser::error(const Tokenizer::parser::location_type &l, const string &m)
{
  driver.error(l, m);
}
