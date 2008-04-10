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

%x STMT
%x EXPR
%x FOR_BODY
%x THEN_BODY
%x ELSE_BODY

%{
// Increments location counter for every token read
#define YY_USER_ACTION yylloc->columns(yyleng);
%}

SPC  [ \t]+
EOL  (\r)?\n
CONT \\\\

%%
 /* Code put at the beginning of yylex() */
%{
  // Reset location before reading token
  yylloc->step();
%}

<INITIAL>^{SPC}*@#{SPC}*include{SPC}+\"[^\"\r\n]*\"{SPC}*{EOL} {
                              yylloc->lines(1);
                              yylloc->step();

                              // Get filename
                              string *filename = new string(yytext);
                              int dblq_idx1 = filename->find('"');
                              int dblq_idx2 = filename->find('"', dblq_idx1 + 1);
                              filename->erase(dblq_idx2);
                              filename->erase(0, dblq_idx1 + 1);

                              create_include_context(filename, yylloc, driver);

                              BEGIN(INITIAL);
                            }

<INITIAL>^{SPC}*@#          { yylloc->step(); BEGIN(STMT); }
<INITIAL>@\{                { yylloc->step(); BEGIN(EXPR); }

<EXPR>\}                    { BEGIN(INITIAL); return token::EOL; }

<STMT>{CONT}{SPC}*{EOL}     { yylloc->lines(1); yylloc->step(); }
<STMT>{EOL}                 {
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
                              else if (reading_if_statement)
                                {
                                  reading_if_statement = false;
                                  then_body_tmp.erase();
                                  then_body_loc_tmp = *yylloc;
                                  nested_if_nb = 0;
                                  BEGIN(THEN_BODY);
                                }
                              else
                                {
                                  *yyout << endl;
                                  BEGIN(INITIAL);
                                }
                              return token::EOL;
                            }

<STMT,EXPR>{SPC}+           { yylloc->step(); }

<STMT,EXPR>[0-9]+           {
                              yylval->int_val = atoi(yytext);
                              return token::INTEGER;
                            }
<STMT,EXPR>\(               { return token::LPAREN; }
<STMT,EXPR>\)               { return token::RPAREN; }
<STMT,EXPR>\[               { return token::LBRACKET; }
<STMT,EXPR>\]               { return token::RBRACKET; }
<STMT,EXPR>:                { return token::COLON; }
<STMT,EXPR>,                { return token::COMMA; }
<STMT,EXPR>=                { return token::EQUAL; }
<STMT,EXPR>[!]              { return token::EXCLAMATION; }
<STMT,EXPR>"||"             { return token::LOGICAL_OR; }
<STMT,EXPR>&&               { return token::LOGICAL_AND; }
<STMT,EXPR>"<="             { return token::LESS_EQUAL; }
<STMT,EXPR>">="             { return token::GREATER_EQUAL; }
<STMT,EXPR>"<"              { return token::LESS; }
<STMT,EXPR>">"              { return token::GREATER; }
<STMT,EXPR>"=="             { return token::EQUAL_EQUAL; }
<STMT,EXPR>"!="             { return token::EXCLAMATION_EQUAL; }
<STMT,EXPR>[+]              { return token::PLUS; }
<STMT,EXPR>[-]              { return token::MINUS; }
<STMT,EXPR>[*]              { return token::TIMES; }
<STMT,EXPR>[/]              { return token::DIVIDE; }

<STMT,EXPR>\"[^\"]*\"       {
                              yylval->string_val = new string(yytext + 1);
                              yylval->string_val->resize(yylval->string_val->length() - 1);
                              return token::STRING;
                            }

<STMT>line                  { return token::LINE; }
<STMT>define                { return token::DEFINE; }

<STMT>for                   { reading_for_statement = true; return token::FOR; }
<STMT>in                    { return token::IN; }
<STMT>endfor                { driver.error(*yylloc, "@#endfor is not matched by a @#for statement"); }

<STMT>if                    { reading_if_statement = true; return token::IF; }
<STMT>else                  { driver.error(*yylloc, "@#else is not matched by an @#if statement"); }
<STMT>endif                 { driver.error(*yylloc, "@#endif is not matched by an @#if statement"); }

<STMT>echo                  { return token::ECHO_DIR; }
<STMT>error                 { return token::ERROR; }

<STMT,EXPR>[A-Za-z_][A-Za-z0-9_]* {
                              yylval->string_val = new string(yytext);
                              return token::NAME;
                            }

<EXPR><<EOF>>               { driver.error(*yylloc, "Unexpected end of file while parsing a macro expression"); }
<STMT><<EOF>>               { driver.error(*yylloc, "Unexpected end of file while parsing a macro statement"); }

<FOR_BODY>{EOL}             { yylloc->lines(1); yylloc->step(); for_body_tmp.append(yytext); }
<FOR_BODY>^{SPC}*@#{SPC}*for({SPC}|{CONT}) {
                              nested_for_nb++;
                              for_body_tmp.append(yytext);
                              yylloc->step();
                            }
<FOR_BODY>.                 { for_body_tmp.append(yytext); yylloc->step(); }
<FOR_BODY><<EOF>>           { driver.error(*yylloc, "Unexpected end of file: @#for loop not matched by an @#endfor"); }
<FOR_BODY>^{SPC}*@#{SPC}*endfor{SPC}*{EOL} {
                              yylloc->lines(1);
                              yylloc->step();
                              if (nested_for_nb)
                                {
                                  /* This @#endfor is not the end of the loop body,
                                     but only that of a nested @#for loop */
                                  nested_for_nb--;
                                  for_body_tmp.append(yytext);
                                }
                              else
                                {
                                  // Save old buffer state and location
                                  save_context(yylloc);

                                  for_body = for_body_tmp;
                                  for_body_loc = for_body_loc_tmp;

                                  iter_loop(driver, yylloc);

                                  BEGIN(INITIAL);
                                }
                            }

<THEN_BODY>{EOL}            { yylloc->lines(1); yylloc->step(); then_body_tmp.append(yytext); }
<THEN_BODY>^{SPC}*@#{SPC}*if({SPC}|{CONT}) {
                              nested_if_nb++;
                              then_body_tmp.append(yytext);
                              yylloc->step();
                            }
<THEN_BODY>.                { then_body_tmp.append(yytext); yylloc->step(); }
<THEN_BODY><<EOF>>          { driver.error(*yylloc, "Unexpected end of file: @#if not matched by an @#endif"); }
<THEN_BODY>^{SPC}*@#{SPC}*else{SPC}*{EOL} {
                              yylloc->lines(1);
                              yylloc->step();
                              if (nested_if_nb)
                                then_body_tmp.append(yytext);
                              else
                                {
                                  else_body_tmp.erase();
                                  else_body_loc_tmp = *yylloc;
                                  BEGIN(ELSE_BODY);
                                }
                             }

<THEN_BODY>^{SPC}*@#{SPC}*endif{SPC}*{EOL} {
                              yylloc->lines(1);
                              yylloc->step();
                              if (nested_if_nb)
                                {
                                  /* This @#endif is not the end of the @#if we're parsing,
                                     but only that of a nested @#if */
                                  nested_if_nb--;
                                  then_body_tmp.append(yytext);
                                }
                              else
                                {
                                  if (driver.last_if)
                                    create_then_context(yylloc);
                                  else
                                    output_line(yylloc);

                                  BEGIN(INITIAL);
                                }
                            }

<ELSE_BODY>{EOL}            { yylloc->lines(1); yylloc->step(); else_body_tmp.append(yytext); }
<ELSE_BODY>^{SPC}*@#{SPC}*if({SPC}|{CONT}) {
                              nested_if_nb++;
                              else_body_tmp.append(yytext);
                              yylloc->step();
                            }
<ELSE_BODY>.                { else_body_tmp.append(yytext); yylloc->step(); }
<ELSE_BODY><<EOF>>          { driver.error(*yylloc, "Unexpected end of file: @#if not matched by an @#endif"); }

<ELSE_BODY>^{SPC}*@#{SPC}*endif{SPC}*{EOL} {
                              yylloc->lines(1);
                              yylloc->step();
                              if (nested_if_nb)
                                {
                                  /* This @#endif is not the end of the @#if we're parsing,
                                     but only that of a nested @#if */
                                  nested_if_nb--;
                                  else_body_tmp.append(yytext);
                                }
                              else
                                {
                                  if (driver.last_if)
                                    create_then_context(yylloc);
                                  else
                                    create_else_context(yylloc);

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

                              /* If we are not in a loop body, or if the loop has terminated,
                                 pop a context */
                              if (for_body.empty() || !iter_loop(driver, yylloc))
                                restore_context(yylloc);
                            }

 /* We don't use echo, because under Cygwin it will add an extra \r */
<INITIAL>{EOL}              { yylloc->lines(1); yylloc->step(); *yyout << endl; }

 /* Copy everything else to output */
<INITIAL>.                  { yylloc->step(); ECHO; }

<*>.                        { driver.error(*yylloc, "Macro lexer error: '" + string(yytext) + "'"); }
%%

MacroFlex::MacroFlex(istream* in, ostream* out)
  : MacroFlexLexer(in, out), input(in), reading_for_statement(false), reading_if_statement(false)
{
}

void
MacroFlex::output_line(Macro::parser::location_type *yylloc) const
{
  *yyout << endl << "@#line \"" << *yylloc->begin.filename << "\" "
         << yylloc->begin.line << endl;
}

void
MacroFlex::save_context(Macro::parser::location_type *yylloc)
{
  context_stack.push(ScanContext(input, YY_CURRENT_BUFFER, *yylloc, for_body, for_body_loc));
}

void
MacroFlex::restore_context(Macro::parser::location_type *yylloc)
{
  input = context_stack.top().input;
  yy_switch_to_buffer(context_stack.top().buffer);
  *yylloc = context_stack.top().yylloc;
  for_body = context_stack.top().for_body;
  for_body_loc = context_stack.top().for_body_loc;
  // Remove top of stack
  context_stack.pop();
  // Dump @#line instruction
  output_line(yylloc);
}

void
MacroFlex::create_include_context(string *filename, Macro::parser::location_type *yylloc,
                                  MacroDriver &driver)
{
  save_context(yylloc);
  // Open new file
  input = new ifstream(filename->c_str(), ios::binary);
  if (input->fail())
    driver.error(*yylloc, "Could not open " + *filename);
  // Reset location
  yylloc->begin.filename = yylloc->end.filename = filename;
  yylloc->begin.line = yylloc->end.line = 1;
  yylloc->begin.column = yylloc->end.column = 0;
  // We are not in a loop body
  for_body.clear();
  // Output @#line information
  output_line(yylloc);
  // Switch to new buffer
  yy_switch_to_buffer(yy_create_buffer(input, YY_BUF_SIZE));
}

void
MacroFlex::create_then_context(Macro::parser::location_type *yylloc)
{
  save_context(yylloc);
  input = new stringstream(then_body_tmp);
  *yylloc = then_body_loc_tmp;
  yylloc->begin.filename = yylloc->end.filename = new string(*then_body_loc_tmp.begin.filename);
  for_body.clear();
  output_line(yylloc);
  yy_switch_to_buffer(yy_create_buffer(input, YY_BUF_SIZE));
}

void
MacroFlex::create_else_context(Macro::parser::location_type *yylloc)
{
  save_context(yylloc);
  input = new stringstream(else_body_tmp);
  *yylloc = else_body_loc_tmp;
  yylloc->begin.filename = yylloc->end.filename = new string(*else_body_loc_tmp.begin.filename);
  for_body.clear();
  output_line(yylloc);
  yy_switch_to_buffer(yy_create_buffer(input, YY_BUF_SIZE));
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
