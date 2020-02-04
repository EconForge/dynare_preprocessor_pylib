/*
 * Copyright © 2019-2020 Dynare Team
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

#include "Driver.hh"
#include <regex>

using namespace macro;

void
Driver::parse(const string &file_arg, const string &basename_arg, istream &modfile,
              ostream &output, bool debug, const vector<pair<string, string>> &defines,
              vector<filesystem::path> &paths)
{
  file = file_arg;
  basename = basename_arg;

  if (!defines.empty())
    {
      stringstream command_line_defines_with_endl;
      for (const auto & [var, val] : defines)
        command_line_defines_with_endl << "@#define " << var << " = " << val << endl;
      Driver m(env, true);
      istream is(command_line_defines_with_endl.rdbuf());
      m.parse("command_line_defines", "command_line_defines", is, output, debug, {}, paths);
    }

  // Handle empty files
  if (modfile.rdbuf()->in_avail() == 0)
    return;

  stringstream file_with_endl;
  file_with_endl << modfile.rdbuf() << endl;

  lexer = make_unique<TokenizerFlex>(&file_with_endl);
  lexer->set_debug(debug);

  Tokenizer::parser parser(*this);
  parser.set_debug_level(debug);

  // Launch macro-processing
  parser.parse();

  // Interpret parsed statements
  bool printLine = true;
  for (const auto &statement : statements)
    {
      if (printLine)
        {
          statement->printLineInfo(output, no_line_macro);
          printLine = false;
        }
      statement->interpret(output, no_line_macro, paths);
    }
}

void
Driver::error(const Tokenizer::parser::location_type &location, const string &message) const
{
  cerr << "ERROR in macro-processor: " << location << ": " << message << endl;
  exit(EXIT_FAILURE);
}
