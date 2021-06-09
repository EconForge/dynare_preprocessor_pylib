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
          str = regex_replace(str, regex(R"((^|\n)\s*@#line.*)"), "");
          auto compareNewline = [](char i, char j) { return i == '\n' && j == '\n'; };
          str.erase(0, str.find_first_not_of('\n'));
          str.erase(unique(str.begin(), str.end(), compareNewline), str.end());
        }
      macro_output_file << str;
      macro_output_file.close();
    }
  return macro_output;
}
