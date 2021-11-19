/*
 * Copyright © 2015-2020 Dynare Team
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

#include <regex>
#include <sstream>
#include <fstream>
#include <filesystem>
#include <algorithm>

#include "macro/Driver.hh"

stringstream
macroExpandModFile(const string &filename, const string &basename, const istream &modfile,
                   bool debug, bool save_macro, string save_macro_file, bool line_macro,
                   const vector<pair<string, string>> &defines,
                   vector<filesystem::path> paths)
{
  // Do macro processing
  stringstream macro_output;
  macro::Environment env = macro::Environment();
  macro::Driver m;
  m.parse(filename, basename, modfile, debug, defines, env, paths, macro_output);
  if (save_macro)
    {
      if (save_macro_file.empty())
        save_macro_file = basename + "-macroexp.mod";
      ofstream macro_output_file(save_macro_file);
      if (macro_output_file.fail())
        {
          cerr << "Cannot open " << save_macro_file << " for macro output" << endl;
          exit(EXIT_FAILURE);
        }

      string str(macro_output.str());
      if (!line_macro)
        {
          /* Remove the @#line directives.
             Unfortunately GCC 11 does not yet support std::regex::multiline
             (despite it being in the C++17 standard), so we are forced to use
             a trick to emulate the “usual” behaviour of the caret ^;
             here, the latter only matches the beginning of file.
             This also means that we are forced to remove the EOL before the
             @#line, and not the one after it (matching the EOL before and the
             EOL after in the same regexp does not work). */
          str = regex_replace(str, regex(R"((^|\r?\n)@#line.*)"), "");
          /* Remove the EOLs at the beginning of the output, the first one
             being a remnant of the first @#line directive. */
          str = regex_replace(str, regex(R"(^(\r?\n)+)"), "");
          /* Replace sequences of several newlines by a single newline (in
             both LF and CR+LF conventions). */
          str = regex_replace(str, regex(R"(\n{2,})"), "\n");
          str = regex_replace(str, regex(R"((\r\n){2,})"), "\r\n");
        }
      macro_output_file << str;
      macro_output_file.close();
    }
  return macro_output;
}
