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

%x MACRO
%x FOR_BODY

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

<INITIAL>^@include[ \t]+\"[^\"\r\n]*\"[ \t]*(\r)?\n         {
                               yylloc->lines(1);
                               yylloc->step();
                               // Save old buffer state and location
                               context_stack.push(ScanContext(input, YY_CURRENT_BUFFER, *yylloc, for_body, for_body_loc));
                               // Get filename
                               string *filename = new string(yytext);
                               int dblq_idx1 = filename->find('"');
                               int dblq_idx2 = filename->find('"', dblq_idx1 + 1);
                               filename->erase(dblq_idx2);
                               filename->erase(0, dblq_idx1 + 1);
                               // Open new file
                               input = new ifstream(filename->c_str(), ios::binary);
                               if (input->fail())
                                 driver.error(*yylloc, "Could not open " + *filename);
                               // Reset location
                               yylloc->begin.filename = yylloc->end.filename = filename;
                               yylloc->begin.line = yylloc->end.line = 1;
                               yylloc->begin.column = yylloc->end.column = 0;
                               // We are not in a loop body
                               for_body.erase();
                               // Output @line information
                               output_line(yylloc);
                               // Switch to new buffer
                               yy_switch_to_buffer(yy_create_buffer(input, YY_BUF_SIZE));
                               BEGIN(INITIAL);
                            }

<INITIAL>@                  { BEGIN(MACRO); }

<MACRO>[ \t\r\f]+           { yylloc->step(); }
<MACRO>@                    { BEGIN(INITIAL); return token::EOL; }
<MACRO>\n                   {
                              yylloc->lines(1);
                              yylloc->step();
                              if (reading_for_statement)
                                {
                                  reading_for_statement = false;
                                  for_body_tmp.erase();
                                  for_body_loc_tmp = *yylloc;
                                  nested_for_nb = 0;
                                  BEGIN(FOR_BODY);
                                }
                              else
                                BEGIN(INITIAL);
                              *yyout << endl;
                              return token::EOL;
                            }

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

<MACRO>\"[^\"]*\"           {
                              yylval->string_val = new string(yytext + 1);
                              yylval->string_val->resize(yylval->string_val->length() - 1);
                              return token::STRING;
                            }

<MACRO>line                 { return token::LINE; }
<MACRO>define               { return token::DEFINE; }
<MACRO>for                  { reading_for_statement = true; return token::FOR; }
<MACRO>in                   { return token::IN; }
<MACRO>endfor               { driver.error(*yylloc, "@endfor is not matched by a @for statement"); }

<MACRO>[A-Za-z_][A-Za-z0-9_]* {
                                yylval->string_val = new string(yytext);
                                return token::NAME;
                              }

<MACRO><<EOF>>              { driver.error(*yylloc, "Unexpected end of file while parsing a macro expression"); }

<FOR_BODY>[\n]+             { yylloc->lines(yyleng); yylloc->step(); for_body_tmp.append(yytext); }
<FOR_BODY>@for              { nested_for_nb++; for_body_tmp.append(yytext); }
<FOR_BODY>.                 { for_body_tmp.append(yytext); }
<FOR_BODY><<EOF>>           { driver.error(*yylloc, "Unexpected end of file: @for loop not matched by an @endfor"); }
<FOR_BODY>@endfor[ \t]*(\r)?\n {
                                 if (nested_for_nb)
                                   {
                                     nested_for_nb--;
                                     for_body_tmp.append(yytext);
                                   }
                                 else
                                   {
                                     yylloc->lines(1);
                                     yylloc->step();
                                     // Save old buffer state and location
                                     context_stack.push(ScanContext(input, YY_CURRENT_BUFFER, *yylloc, for_body, for_body_loc));

                                     for_body = for_body_tmp;
                                     for_body_loc = for_body_loc_tmp;
                                     
                                     iter_loop(driver, yylloc);

                                     BEGIN(INITIAL);
                                   }
                               }

<INITIAL><<EOF>>            {
                              // Quit lexer if end of main file
                              if (context_stack.empty())
                                {
                                  yyterminate();
                                }
                              // Else clean current scanning context
                              yy_delete_buffer(YY_CURRENT_BUFFER);
                              delete input;
                              delete yylloc->begin.filename;

                              // If we are not in a loop body, or if the loop has terminated, pop a context
                              if (for_body.empty() || !iter_loop(driver, yylloc))
                                {
                                  // Restore old context
                                  input = context_stack.top().input;
                                  yy_switch_to_buffer(context_stack.top().buffer);
                                  *yylloc = context_stack.top().yylloc;
                                  for_body = context_stack.top().for_body;
                                  for_body_loc = context_stack.top().for_body_loc;
                                  // Remove top of stack
                                  context_stack.pop();
                                  // Dump @line instruction
                                  output_line(yylloc);
                                }
                            }

 /* Ignore \r, because under Cygwin, outputting \n automatically adds another \r */
<INITIAL>[\r]+              { yylloc->step(); }

 /* Copy everything else to output */
<INITIAL>[\n]+              { yylloc->lines(yyleng); yylloc->step(); ECHO; }
<INITIAL>.                  { yylloc->step(); ECHO; }

<*>.                        { driver.error(*yylloc, "Macro lexer error: '" + string(yytext) + "'"); }
%%

MacroFlex::MacroFlex(istream* in, ostream* out)
  : MacroFlexLexer(in, out), input(in), reading_for_statement(false)
{
}

void
MacroFlex::output_line(Macro::parser::location_type *yylloc)
{
  *yyout << endl << "@line \"" << *yylloc->begin.filename << "\" "
         << yylloc->begin.line << endl;
}

bool
MacroFlex::iter_loop(MacroDriver &driver, Macro::parser::location_type *yylloc)
{
  if (!driver.iter_loop())
    return false;

  input = new stringstream(for_body);
  *yylloc = for_body_loc;
  yylloc->begin.filename = yylloc->end.filename = new string(*for_body_loc.begin.filename);
  output_line(yylloc);
  yy_switch_to_buffer(yy_create_buffer(input, YY_BUF_SIZE));

  return true;
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
