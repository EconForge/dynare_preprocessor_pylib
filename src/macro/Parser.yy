// -*- C++ -*-
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

%language "c++"
%require "3.2"
%defines
%define api.value.type variant
%define api.namespace {Tokenizer}
%define api.location.file "ParserLocation.hh"
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

%token FOR ENDFOR IF IFDEF IFNDEF ELSEIF ELSE ENDIF TRUE FALSE
%token INCLUDE INCLUDEPATH DEFINE EQUAL D_ECHO ERROR
%token COMMA LPAREN RPAREN LBRACKET RBRACKET WHEN
%token BEGIN_EVAL END_EVAL ECHOMACROVARS SAVE LINE

%token EXP LOG LN LOG10 SIN COS TAN ASIN ACOS ATAN
%token SQRT CBRT SIGN MAX MIN FLOOR CEIL TRUNC SUM MOD
%token ERF ERFC GAMMA LGAMMA ROUND NORMPDF NORMCDF LENGTH

%token ISEMPTY
%token ISBOOLEAN ISREAL ISSTRING ISTUPLE ISARRAY

%token BOOL REAL STRING TUPLE ARRAY

%token DEFINED

%left OR
%left AND
%left EQUAL_EQUAL NOT_EQUAL
%left LESS GREATER LESS_EQUAL GREATER_EQUAL
%nonassoc IN
/* The COLON operator cannot be given a precedence, because it has both a
   binary and a ternary forms. But technically it belongs here given how the
   grammar rules are organized */
%token COLON
%left UNION
%left INTERSECTION
%left PLUS MINUS
%left TIMES DIVIDE
%precedence UNARY NOT
%precedence CAST
%nonassoc POWER

%token EOL
%token <string> NAME TEXT QUOTED_STRING NUMBER

%type <DirectivePtr> statement
%type <DirectivePtr> directive directive_one_line directive_multiline for if ifdef ifndef text eval
%type <vector<pair<ExpressionPtr, vector<DirectivePtr>>>> if_list if_list1
%type <pair<ExpressionPtr, vector<DirectivePtr>>> elseif else
%type <ExpressionPtr> primary_expr oper_expr colon_expr expr for_when
%type <FunctionPtr> function
%type <VariablePtr> symbol

%type <vector<ExpressionPtr>> comma_expr function_args tuple_comma_expr
%type <vector<string>> name_list

%%

%start statements_or_empty_file;

statements_or_empty_file : %empty
                         | statements
                         ;

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
          | text
          | eval
          ;

directive : directive_one_line EOL
          | directive_multiline EOL
          ;

directive_one_line : INCLUDE expr
                     { $$ = make_shared<Include>($2, @$); }
                   | INCLUDEPATH expr
                     { $$ = make_shared<IncludePath>($2, @$); }
                   | DEFINE symbol
                     { $$ = make_shared<Define>($2, make_shared<Bool>(true, @$), @$); }
                   | DEFINE symbol EQUAL expr
                     { $$ = make_shared<Define>($2, $4, @$); }
                   | DEFINE function EQUAL expr
                     { $$ = make_shared<Define>($2, $4, @$); }
                   | D_ECHO expr
                     { $$ = make_shared<Echo>($2, @$); }
                   | ERROR expr
                     { $$ = make_shared<Error>($2, @$); }
                   | ECHOMACROVARS
                     { $$ = make_shared<EchoMacroVars>(false, @$); }
                   | ECHOMACROVARS name_list
                     { $$ = make_shared<EchoMacroVars>(false, $2, @$); }
                   | ECHOMACROVARS LPAREN SAVE RPAREN
                     { $$ = make_shared<EchoMacroVars>(true, @$); }
                   | ECHOMACROVARS LPAREN SAVE RPAREN name_list
                     { $$ = make_shared<EchoMacroVars>(true, $5, @$); }
                   | LINE QUOTED_STRING NUMBER
                     {
                       // `@#line` is ignored; adjust newlines in output to accord
                       auto l = static_cast<Tokenizer::parser::location_type>(@$);
                       $$ = make_shared<TextNode>(string(l.end.line - l.begin.line + 1, '\n'), @$);
                     }
                   ;

name_list : NAME
            { $$ = {$1}; }
          | name_list NAME
            {
              $1.emplace_back($2);
              $$ = $1;
            }
          ;

directive_multiline : for
                    | if
                    | ifdef
                    | ifndef
                    ;

for_when : %empty
           { $$ = shared_ptr<Expression>(); }
         | WHEN expr
           { $$ = $2; }
         ;

for : FOR { driver.pushContext(); } expr IN expr for_when EOL statements ENDFOR
      {
        vector<VariablePtr> vvnp;
        auto tmpt = dynamic_pointer_cast<Tuple>($3);
        auto tmpv = dynamic_pointer_cast<Variable>($3);
        if (tmpv)
          vvnp.emplace_back(tmpv);
        else if (tmpt)
          for (const auto & it : tmpt->getValue())
            {
              auto vnp = dynamic_pointer_cast<Variable>(it);
              if (!vnp)
                error(@$, "For loop indices must be variables");
              vvnp.emplace_back(vnp);
            }
        else
          error(@1, "For loop indices must be a variable or a tuple");

        auto vdp = driver.popContext();
        vdp.emplace_back(make_shared<TextNode>("\n", @9));

        if (!$6)
          $$ = make_shared<For>(vvnp, $5, vdp, @$);
        else
          {
            auto tmpc = make_shared<Comprehension>(true, $3, $5, $6, @6);
            $$ = make_shared<For>(vvnp, tmpc, vdp, @$);
          }
      }
    ;

if : IF { driver.pushContext(); } if_list ENDIF
     { $$ = make_shared<If>($3, @$); }
   ;

ifdef : IFDEF { driver.pushContext(); } if_list ENDIF
        { $$ = make_shared<Ifdef>($3, @$); }
      ;

ifndef : IFNDEF { driver.pushContext(); } if_list ENDIF
         { $$ = make_shared<Ifndef>($3, @$); }
       ;

if_list : if_list1
        | if_list1 else
          {
            $1.emplace_back($2);
            $$ = $1;
          }
        ;

if_list1 : expr EOL
           {
             auto context = driver.popContext();
             context.emplace_back(make_shared<TextNode>("\n", @2));
             $$ = {{$1, context}};
           }
         | expr EOL statements
           {
             auto context = driver.popContext();
             context.emplace_back(make_shared<TextNode>("\n", @3));
             $$ = {{$1, context}};
           }
         | if_list1 elseif
           {
             $1.emplace_back($2);
             $$ = $1;
           }
         ;

elseif_begin : ELSEIF { driver.pushContext(); } ;

elseif : elseif_begin expr EOL
         {
           auto context = driver.popContext();
           context.emplace_back(make_shared<TextNode>("\n", @3));
           $$ = {$2, context};
         }
       | elseif_begin expr EOL statements
         {
           auto context = driver.popContext();
           context.emplace_back(make_shared<TextNode>("\n", @4));
           $$ = {$2, context};
         }
       ;

else_begin : ELSE { driver.pushContext(); } ;

else : else_begin EOL
       {
         auto context = driver.popContext();
         context.emplace_back(make_shared<TextNode>("\n", @2));
         $$ = {make_shared<Bool>(true, @1), context};
       }
     | else_begin EOL statements
       {
         auto context = driver.popContext();
         context.emplace_back(make_shared<TextNode>("\n", @3));
         $$ = {make_shared<Bool>(true, @1), context};
       }
     ;

text : TEXT
       { $$ = make_shared<TextNode>($1, @$); }
     | EOL
       { $$ = make_shared<TextNode>("\n", @$); }
     ;

eval : BEGIN_EVAL expr END_EVAL
       { $$ = make_shared<Eval>($2, @$); }
     ;

symbol : NAME
         { $$ = make_shared<Variable>($1, @$); }
       ;

function : NAME LPAREN RPAREN
           { $$ = make_shared<Function>($1, vector<ExpressionPtr>(), @$); }
         | NAME LPAREN function_args RPAREN
           { $$ = make_shared<Function>($1, $3, @$); }
         ;

function_args : symbol
                { $$ = {$1}; }
              | function_args COMMA symbol
                { $1.emplace_back($3); $$ = $1; }
              ;

comma_expr : %empty
             { $$ = {}; }
           | expr
             { $$ = {$1}; }
           | comma_expr COMMA expr
             { $1.emplace_back($3); $$ = $1; }
           ;

tuple_comma_expr : %empty
                   { $$ = {}; }
                 | expr COMMA
                   { $$ = {$1}; }
                 | expr COMMA expr
                   { $$ = {$1, $3}; }
                 | tuple_comma_expr COMMA expr
                   { $1.emplace_back($3); $$ = $1; }
                 ;

primary_expr : LPAREN expr RPAREN
               { $$ = $2; }
             | symbol
               { $$ = $1; } // Explicit rule needed for type conversion
             | NAME LBRACKET comma_expr RBRACKET
               {
                 $$ = make_shared<Variable>($1, make_shared<Array>($3, @3), @$);
               }
             | NAME LPAREN comma_expr RPAREN
               { $$ = make_shared<Function>($1, $3, @$); }
             | TRUE
               { $$ = make_shared<Bool>(true, @$); }
             | FALSE
               { $$ = make_shared<Bool>(false, @$); }
             | NUMBER
               { $$ = make_shared<Real>($1, @$); }
             | QUOTED_STRING
               { $$ = make_shared<String>($1, @$); }
             | LBRACKET comma_expr RBRACKET
               { $$ = make_shared<Array>($2, @$); }
             | LPAREN tuple_comma_expr RPAREN
               { $$ = make_shared<Tuple>($2, @$); }
             | LBRACKET expr IN expr WHEN expr RBRACKET
               { $$ = make_shared<Comprehension>(true, $2, $4, $6, @$); }
             | LBRACKET expr FOR expr IN expr RBRACKET
               { $$ = make_shared<Comprehension>($2, $4, $6, @$); }
             | LBRACKET expr FOR expr IN expr WHEN expr RBRACKET
               { $$ = make_shared<Comprehension>($2, $4, $6, $8, @$); }
             | LENGTH LPAREN expr RPAREN
               { $$ = make_shared<UnaryOp>(codes::UnaryOp::length, $3, @$); }
             | ISEMPTY LPAREN expr RPAREN
               { $$ = make_shared<UnaryOp>(codes::UnaryOp::isempty, $3, @$); }
             | ISBOOLEAN LPAREN expr RPAREN
               { $$ = make_shared<UnaryOp>(codes::UnaryOp::isboolean, $3, @$); }
             | ISREAL LPAREN expr RPAREN
               { $$ = make_shared<UnaryOp>(codes::UnaryOp::isreal, $3, @$); }
             | ISSTRING LPAREN expr RPAREN
               { $$ = make_shared<UnaryOp>(codes::UnaryOp::isstring, $3, @$); }
             | ISTUPLE LPAREN expr RPAREN
               { $$ = make_shared<UnaryOp>(codes::UnaryOp::istuple, $3, @$); }
             | ISARRAY LPAREN expr RPAREN
               { $$ = make_shared<UnaryOp>(codes::UnaryOp::isarray, $3, @$); }
             | EXP LPAREN expr RPAREN
               { $$ = make_shared<UnaryOp>(codes::UnaryOp::exp, $3, @$); }
             | LOG LPAREN expr RPAREN
               { $$ = make_shared<UnaryOp>(codes::UnaryOp::ln, $3, @$); }
             | LN LPAREN expr RPAREN
               { $$ = make_shared<UnaryOp>(codes::UnaryOp::ln, $3, @$); }
             | LOG10 LPAREN expr RPAREN
               { $$ = make_shared<UnaryOp>(codes::UnaryOp::log10, $3, @$); }
             | SIN LPAREN expr RPAREN
               { $$ = make_shared<UnaryOp>(codes::UnaryOp::sin, $3, @$); }
             | COS LPAREN expr RPAREN
               { $$ = make_shared<UnaryOp>(codes::UnaryOp::cos, $3, @$); }
             | TAN LPAREN expr RPAREN
               { $$ = make_shared<UnaryOp>(codes::UnaryOp::tan, $3, @$); }
             | ASIN LPAREN expr RPAREN
               { $$ = make_shared<UnaryOp>(codes::UnaryOp::asin, $3, @$); }
             | ACOS LPAREN expr RPAREN
               { $$ = make_shared<UnaryOp>(codes::UnaryOp::acos, $3, @$); }
             | ATAN LPAREN expr RPAREN
               { $$ = make_shared<UnaryOp>(codes::UnaryOp::atan, $3, @$); }
             | SQRT LPAREN expr RPAREN
               { $$ = make_shared<UnaryOp>(codes::UnaryOp::sqrt, $3, @$); }
             | CBRT LPAREN expr RPAREN
               { $$ = make_shared<UnaryOp>(codes::UnaryOp::cbrt, $3, @$); }
             | SIGN LPAREN expr RPAREN
               { $$ = make_shared<UnaryOp>(codes::UnaryOp::sign, $3, @$); }
             | FLOOR LPAREN expr RPAREN
               { $$ = make_shared<UnaryOp>(codes::UnaryOp::floor, $3, @$); }
             | CEIL LPAREN expr RPAREN
               { $$ = make_shared<UnaryOp>(codes::UnaryOp::ceil, $3, @$); }
             | TRUNC LPAREN expr RPAREN
               { $$ = make_shared<UnaryOp>(codes::UnaryOp::trunc, $3, @$); }
             | SUM LPAREN expr RPAREN
               { $$ = make_shared<UnaryOp>(codes::UnaryOp::sum, $3, @$); }
             | ERF LPAREN expr RPAREN
               { $$ = make_shared<UnaryOp>(codes::UnaryOp::erf, $3, @$); }
             | ERFC LPAREN expr RPAREN
               { $$ = make_shared<UnaryOp>(codes::UnaryOp::erfc, $3, @$); }
             | GAMMA LPAREN expr RPAREN
               { $$ = make_shared<UnaryOp>(codes::UnaryOp::gamma, $3, @$); }
             | LGAMMA LPAREN expr RPAREN
               { $$ = make_shared<UnaryOp>(codes::UnaryOp::lgamma, $3, @$); }
             | ROUND LPAREN expr RPAREN
               { $$ = make_shared<UnaryOp>(codes::UnaryOp::round, $3, @$); }
             | NORMPDF LPAREN expr RPAREN
               { $$ = make_shared<UnaryOp>(codes::UnaryOp::normpdf, $3, @$); }
             | NORMCDF LPAREN expr RPAREN
               { $$ = make_shared<UnaryOp>(codes::UnaryOp::normcdf, $3, @$); }
             | MAX LPAREN expr COMMA expr RPAREN
               { $$ = make_shared<BinaryOp>(codes::BinaryOp::max, $3, $5, @$); }
             | MIN LPAREN expr COMMA expr RPAREN
               { $$ = make_shared<BinaryOp>(codes::BinaryOp::min, $3, $5, @$); }
             | MOD LPAREN expr COMMA expr RPAREN
               { $$ = make_shared<BinaryOp>(codes::BinaryOp::mod, $3, $5, @$); }
             | NORMPDF LPAREN expr COMMA expr COMMA expr RPAREN
               { $$ = make_shared<TrinaryOp>(codes::TrinaryOp::normpdf, $3, $5, $7, @$); }
             | NORMCDF LPAREN expr COMMA expr COMMA expr RPAREN
               { $$ = make_shared<TrinaryOp>(codes::TrinaryOp::normcdf, $3, $5, $7, @$); }
             | DEFINED LPAREN NAME RPAREN
               { $$ = make_shared<UnaryOp>(codes::UnaryOp::defined, make_shared<String>($3, @3), @$); }
             ;

oper_expr : primary_expr
          | LPAREN BOOL RPAREN oper_expr %prec CAST
            { $$ = make_shared<UnaryOp>(codes::UnaryOp::cast_bool, $4, @$); }
          | LPAREN REAL RPAREN oper_expr %prec CAST
            { $$ = make_shared<UnaryOp>(codes::UnaryOp::cast_real, $4, @$); }
          | LPAREN STRING RPAREN oper_expr %prec CAST
            { $$ = make_shared<UnaryOp>(codes::UnaryOp::cast_string, $4, @$); }
          | LPAREN TUPLE RPAREN oper_expr %prec CAST
            { $$ = make_shared<UnaryOp>(codes::UnaryOp::cast_tuple, $4, @$); }
          | LPAREN ARRAY RPAREN oper_expr %prec CAST
            { $$ = make_shared<UnaryOp>(codes::UnaryOp::cast_array, $4, @$); }
          | NOT oper_expr
            { $$ = make_shared<UnaryOp>(codes::UnaryOp::logical_not, $2, @$); }
          | MINUS oper_expr %prec UNARY
            { $$ = make_shared<UnaryOp>(codes::UnaryOp::unary_minus, $2, @$); }
          | PLUS oper_expr %prec UNARY
            { $$ = make_shared<UnaryOp>(codes::UnaryOp::unary_plus, $2, @$); }
          | oper_expr PLUS oper_expr
            { $$ = make_shared<BinaryOp>(codes::BinaryOp::plus, $1, $3, @$); }
          | oper_expr MINUS oper_expr
            { $$ = make_shared<BinaryOp>(codes::BinaryOp::minus, $1, $3, @$); }
          | oper_expr TIMES oper_expr
            { $$ = make_shared<BinaryOp>(codes::BinaryOp::times, $1, $3, @$); }
          | oper_expr DIVIDE oper_expr
            { $$ = make_shared<BinaryOp>(codes::BinaryOp::divide, $1, $3, @$); }
          | oper_expr POWER oper_expr
            { $$ = make_shared<BinaryOp>(codes::BinaryOp::power, $1, $3, @$); }
          | oper_expr UNION oper_expr
            { $$ = make_shared<BinaryOp>(codes::BinaryOp::set_union, $1, $3, @$); }
          | oper_expr INTERSECTION oper_expr
            { $$ = make_shared<BinaryOp>(codes::BinaryOp::set_intersection, $1, $3, @$); }
          ;

colon_expr : oper_expr COLON oper_expr
             { $$ = make_shared<Range>($1, $3, @$); }
           | oper_expr COLON oper_expr COLON oper_expr
             { $$ = make_shared<Range>($1, $3, $5, @$); }
           ;

expr : oper_expr
     | colon_expr
     | expr EQUAL_EQUAL expr
       { $$ = make_shared<BinaryOp>(codes::BinaryOp::equal_equal, $1, $3, @$); }
     | expr NOT_EQUAL expr
       { $$ = make_shared<BinaryOp>(codes::BinaryOp::not_equal, $1, $3, @$); }
     | expr LESS expr
       { $$ = make_shared<BinaryOp>(codes::BinaryOp::less, $1, $3, @$); }
     | expr GREATER expr
       { $$ = make_shared<BinaryOp>(codes::BinaryOp::greater, $1, $3, @$); }
     | expr LESS_EQUAL expr
       { $$ = make_shared<BinaryOp>(codes::BinaryOp::less_equal, $1, $3, @$); }
     | expr GREATER_EQUAL expr
       { $$ = make_shared<BinaryOp>(codes::BinaryOp::greater_equal, $1, $3, @$); }
     | expr AND expr
       { $$ = make_shared<BinaryOp>(codes::BinaryOp::logical_and, $1, $3, @$); }
     | expr OR expr
       { $$ = make_shared<BinaryOp>(codes::BinaryOp::logical_or, $1, $3, @$); }
     | expr IN expr
       { $$ = make_shared<BinaryOp>(codes::BinaryOp::in, $1, $3, @$); }
     ;

%%

void
Tokenizer::parser::error(const Tokenizer::parser::location_type &l, const string &m)
{
  driver.error(l, m);
}
