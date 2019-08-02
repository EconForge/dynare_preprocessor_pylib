/* -*- C++ -*- */
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

%{
#include "Driver.hh"

// Announce to Flex the prototype we want for lexing function
#define YY_DECL                                                \
  Tokenizer::parser::token_type                                \
  TokenizerFlex::lex(Tokenizer::parser::semantic_type *yylval, \
                     Tokenizer::parser::location_type *yylloc, \
                     macro::Driver &driver)

// Shortcut to access tokens defined by Bison
using token = Tokenizer::parser::token;

/* By default yylex returns int, we use token_type.
   Unfortunately yyterminate by default returns 0, which is
   not of token_type.  */
#define yyterminate() return Tokenizer::parser::token_type (0);
%}

%option c++

%option prefix="Tokenizer"

%option case-insensitive noinput noyywrap nounput batch debug never-interactive noyymore

%x directive
%x eval
%x expr
%x end_line

%{
  // Increments location counter for every token read
  # define YY_USER_ACTION yylloc->columns(yyleng);
%}

SPC  [ \t]+
EOL  (\r)?\n
CONT \\\\{SPC}*

%%
 /* Code put at the beginning of yylex() */
%{
  // Reset location before reading token
  yylloc->step();
%}

<directive>include         { BEGIN(expr); return token::INCLUDE; }
<directive>includepath     { BEGIN(expr); return token::INCLUDEPATH; }
<directive>define          { BEGIN(expr); return token::DEFINE; }
<directive>echo            { BEGIN(expr); return token::D_ECHO; }
<directive>error           { BEGIN(expr); return token::ERROR; }
<directive>if              { BEGIN(expr); return token::IF; }
<directive>ifdef           { BEGIN(expr); return token::IFDEF; }
<directive>ifndef          { BEGIN(expr); return token::IFNDEF; }
<directive>else            { BEGIN(end_line); return token::ELSE; }
<directive>endif           { BEGIN(end_line); return token::ENDIF; }
<directive>for             { BEGIN(expr); return token::FOR; }
<directive>endfor          { BEGIN(end_line); return token::ENDFOR; }
<directive>echomacrovars   { BEGIN(expr); return token::ECHOMACROVARS; }

<expr,eval>\+              { return token::PLUS; }
<expr,eval>-               { return token::MINUS; }
<expr,eval>\*              { return token::TIMES; }
<expr,eval>\/              { return token::DIVIDE; }
<expr,eval>=               { return token::EQUAL; }
<expr,eval>\^              { return token::POWER; }
<expr,eval><               { return token::LESS; }
<expr,eval>>               { return token::GREATER; }
<expr,eval>>=              { return token::GREATER_EQUAL; }
<expr,eval><=              { return token::LESS_EQUAL; }
<expr,eval>==              { return token::EQUAL_EQUAL; }
<expr,eval>!=              { return token::NOT_EQUAL; }

<expr,eval>&&              { return token::AND; }
<expr,eval>"||"            { return token::OR; }
<expr,eval>!               { return token::NOT; }

<expr,eval>\|              { return token::UNION; }
<expr,eval>&               { return token::INTERSECTION; }

<expr,eval>,               { return token::COMMA; }
<expr,eval>:               { return token::COLON; }

<expr,eval>\(              { return token::LPAREN; }
<expr,eval>\)              { return token::RPAREN; }
<expr,eval>\[              { return token::LBRACKET; }
<expr,eval>\]              { return token::RBRACKET; }
<expr,eval>in              { return token::IN; }
<expr,eval>for             { return token::FOR; }
<expr,eval>when            { return token::WHEN; }
<expr,eval>save            { return token::SAVE; }

<expr,eval>true            { return token::TRUE; }
<expr,eval>false           { return token::FALSE; }

<expr,eval>exp             { return token::EXP; }
<expr,eval>log             { return token::LOG; }
<expr,eval>ln              { return token::LN; }
<expr,eval>log10           { return token::LOG10; }
<expr,eval>sin             { return token::SIN; }
<expr,eval>cos             { return token::COS; }
<expr,eval>tan             { return token::TAN; }
<expr,eval>asin            { return token::ATAN; }
<expr,eval>acos            { return token::ACOS; }
<expr,eval>atan            { return token::ATAN; }
<expr,eval>sqrt            { return token::SQRT; }
<expr,eval>cbrt            { return token::CBRT; }
<expr,eval>sign            { return token::SIGN; }
<expr,eval>max             { return token::MAX; }
<expr,eval>min             { return token::MIN; }
<expr,eval>floor           { return token::FLOOR; }
<expr,eval>ceil            { return token::CEIL; }
<expr,eval>trunc           { return token::TRUNC; }
<expr,eval>mod             { return token::MOD; }
<expr,eval>sum             { return token::SUM; }
<expr,eval>erf             { return token::ERF; }
<expr,eval>erfc            { return token::ERFC; }
<expr,eval>gamma           { return token::GAMMA; }
<expr,eval>lgamma          { return token::LGAMMA; }
<expr,eval>round           { return token::ROUND; }
<expr,eval>length          { return token::LENGTH; }
<expr,eval>normpdf         { return token::NORMPDF; }
<expr,eval>normcdf         { return token::NORMCDF; }

<expr,eval>int             { return token::INT; }

<expr,eval>((([0-9]*\.[0-9]+)|([0-9]+\.))([ed][-+]?[0-9]+)?)|([0-9]+([ed][-+]?[0-9]+)?)|nan|inf {
  yylval->build<string>(yytext);
  return token::NUMBER;
}

<expr,eval>[A-Za-z_][A-Za-z0-9_]* {
  yylval->build<string>(yytext);
  return token::NAME;
}

<expr,eval>\"[^\"]+\" {
  yylval->build<string>(yytext + 1).pop_back();
  return token::QUOTED_STRING;
}

<expr,eval>{SPC}+                         { }
<eval>{EOL}+                              { yylloc->lines(yyleng); yylloc->lines(yyleng); }
<eval>\}                                  { BEGIN(INITIAL); return token::END_EVAL; }

<expr,end_line>{CONT}("//".*)?{SPC}*{EOL} { yylloc->lines(1); yylloc->step(); }
<expr,end_line>{SPC}*("//".*)?{EOL}       {
                                            yylval->build<string>("\n");
                                            yylloc->lines(1);
                                            BEGIN(INITIAL);
                                            return token::EOL;
                                          }

<INITIAL>^{SPC}*@#{SPC}*                  { BEGIN(directive); }
<INITIAL>@\{                              { BEGIN(eval); return token::BEGIN_EVAL; }
<INITIAL>{SPC}*{EOL}                      {
                                            yylval->build<string>(yytext);
                                            yylloc->lines(1);
                                            return token::EOL;
                                          }
<INITIAL><<EOF>>                          { yyterminate(); }

<directive,expr,eval,end_line><<EOF>>     { driver.error(*yylloc, "unexpected end of file"); }

<*>.                                      { yylval->build<string>(yytext); return token::TEXT; }
<*>.|{EOL}                                { driver.error(*yylloc, "character unrecognized by lexer"); }

%%

/* This implementation of TokenizerFlexLexer::yylex() is required to fill the
 * vtable of the class TokenizerFlexLexer. We define the scanner's main yylex
 * function via YY_DECL to reside in the TokenizerFlex class instead. */

#ifdef yylex
# undef yylex
#endif

int
TokenizerFlexLexer::yylex()
{
  cerr << "TokenizerFlexLexer::yylex() has been called; shouldn't arrive here." << endl;
  exit(EXIT_FAILURE);
}
