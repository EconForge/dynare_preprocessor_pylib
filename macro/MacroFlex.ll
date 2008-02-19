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

%{
using namespace std;

#include <fstream>

#include "MacroDriver.hh"
#include "MacroBison.hh"

// Announce to Flex the prototype we want for lexing function
#define YY_DECL                                              \
  Macro::parser::token_type                                  \
  MacroFlex::lex(Macro::parser::semantic_type *yylval,       \
                 Macro::parser::location_type *yylloc,       \
                 MacroDriver &driver)

// Shortcut to access tokens defined by Bison
typedef Macro::parser::token token;

/* By default yylex returns int, we use token_type.
   Unfortunately yyterminate by default returns 0, which is
   not of token_type.  */
#define yyterminate() return Macro::parser::token_type (0);
%}

%option c++

%option prefix="Macro"

%option case-insensitive noyywrap nounput batch debug never-interactive

%x INCLUDE
%x END_INCLUDE
%x MACRO

%{
// Increments location counter for every token read
#define YY_USER_ACTION yylloc->columns(yyleng);
%}
%%
 /* Code put at the beginning of yylex() */
%{
  // Reset location before reading token
  yylloc->step();
%}

 /* Ignore @line declarations, replace them by a blank line */
<INITIAL>^@line[^\r\n]*(\r)?\n     { yylloc->lines(1); yylloc->step(); *yyout << endl; }

<INITIAL>^@include[ \t]+\"   BEGIN(INCLUDE);

<INCLUDE>[^\"\r\n]*          {
                               driver.ifs = new ifstream(yytext, ios::binary);
                               if (driver.ifs->fail())
                                 driver.error(*yylloc, "Could not open " + string(yytext));
                               // Save old location
                               yylloc->step();
                               driver.loc_stack.push(*yylloc);
                               // Reset location
                               yylloc->begin.filename = yylloc->end.filename = new string(yytext);
                               yylloc->begin.line = yylloc->end.line = 1;
                               yylloc->begin.column = yylloc->end.column = 0;
                               // Display @line
                               *yyout << "@line \"" << *yylloc->begin.filename << "\" 1" << endl;
                               // Switch to new buffer
                               /* We don't use yypush_buffer_state(), since it doesn't exist in
                                  Flex 2.5.4 (see Flex 2.5.33 info file - section 11 - for code
                                  example with yypush_buffer_state()) */
                               state_stack.push(YY_CURRENT_BUFFER);
                               yy_switch_to_buffer(yy_create_buffer(driver.ifs, YY_BUF_SIZE));
                               BEGIN(INITIAL);
                            }

<END_INCLUDE>\"[^\r\n]*(\r)?\n  {
                                  yylloc->lines(1);
                                  yylloc->step();
                                  *yyout << "@line \"" << *yylloc->begin.filename << "\" "
                                         << yylloc->begin.line << endl;
                                  BEGIN(INITIAL);
                                }

<INITIAL>@                  { BEGIN(MACRO); }


<MACRO>[ \t\r\f]+           { yylloc->step(); }
<MACRO>@                    { BEGIN(INITIAL); return token::EOL; }
<MACRO>\n                   { BEGIN(INITIAL); return token::EOL; }

<MACRO>[0-9]+               {
                              yylval->int_val = atoi(yytext);
                              return token::INTEGER;
                            }
<MACRO>\(                   { return token::LPAREN; }
<MACRO>\)                   { return token::RPAREN; }
<MACRO>\[                   { return token::LBRACKET; }
<MACRO>\]                   { return token::RBRACKET; }
<MACRO>:                    { return token::COLON; }
<MACRO>,                    { return token::COMMA; }
<MACRO>=                    { return token::EQUAL; }
<MACRO>[!]                  { return token::EXCLAMATION; }
<MACRO>"||"                 { return token::LOGICAL_OR; }
<MACRO>&&                   { return token::LOGICAL_AND; }
<MACRO>"<="                 { return token::LESS_EQUAL; }
<MACRO>">="                 { return token::GREATER_EQUAL; }
<MACRO>"<"                  { return token::LESS; }
<MACRO>">"                  { return token::GREATER; }
<MACRO>"=="                 { return token::EQUAL_EQUAL; }
<MACRO>"!="                 { return token::EXCLAMATION_EQUAL; }
<MACRO>[+]                  { return token::PLUS; }
<MACRO>[-]                  { return token::MINUS; }
<MACRO>[*]                  { return token::TIMES; }
<MACRO>[/]                  { return token::DIVIDE; }

<MACRO>\'[^\']*\'           {
                              yylval->string_val = new string(yytext + 1);
                              yylval->string_val->resize(yylval->string_val->length() - 1);
                              return token::STRING;
                            }

<MACRO>define               { return token::DEFINE; }

<MACRO>[A-Za-z_][A-Za-z0-9_]* {
                                yylval->string_val = new string(yytext);
                                return token::NAME;
                              }


<<EOF>>                     {
                              /* We don't use yypop_buffer_state(), since it doesn't exist in
                                 Flex 2.5.4 (see Flex 2.5.33 info file - section 11 - for code
                                 example with yypop_buffer_state()) */
                              // Quit lexer if end of main file
                              if (state_stack.empty())
                                {
                                  yyterminate();
                                }
                              // Else restore old flex buffer
                              yy_delete_buffer(YY_CURRENT_BUFFER);
                              yy_switch_to_buffer(state_stack.top());
                              state_stack.pop();
                              // And restore old location
                              delete yylloc->begin.filename;
                              *yylloc = driver.loc_stack.top();
                              driver.loc_stack.pop();
                              BEGIN(END_INCLUDE);
                            }

 /* Ignore \r, because under Cygwin, outputting \n automatically adds another \r */
<INITIAL>[\r]+              { yylloc->step(); }

 /* Copy everything else to output */
<INITIAL>[\n]+              { yylloc->lines(yyleng); yylloc->step(); ECHO; }
<INITIAL>.                  { yylloc->step(); ECHO; }

<*>.                        { driver.error(*yylloc, "Macro lexer error: '" + string(yytext) + "'"); }
%%

MacroFlex::MacroFlex(istream* in, ostream* out)
  : MacroFlexLexer(in, out)
{
}

/* This implementation of MacroFlexLexer::yylex() is required to fill the
 * vtable of the class MacroFlexLexer. We define the scanner's main yylex
 * function via YY_DECL to reside in the MacroFlex class instead. */

#ifdef yylex
# undef yylex
#endif

int
MacroFlexLexer::yylex()
{
  cerr << "MacroFlexLexer::yylex() has been called, that should never happen!" << endl;
  exit(-1);
}

/*
  Local variables:
  mode: C++
  End:
*/
