/*
 * Copyright © 2015-2019 Dynare Team
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

#include <sstream>
#include <fstream>
#include <filesystem>

#include "macro/Driver.hh"

void
main1(const string &filename, const string &basename, istream &modfile, bool debug, bool save_macro, string &save_macro_file,
      bool no_line_macro_arg, bool no_empty_line_macro, const vector<pair<string, string>> &defines,
      vector<filesystem::path> &paths, stringstream &macro_output)
{
  // Do macro processing
  macro::Environment env = macro::Environment();
  macro::Driver m(env, no_line_macro_arg);
  m.parse(filename, basename, modfile, macro_output, debug, defines, paths);
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

      string str (macro_output.str());
      if (no_empty_line_macro)
        {
          auto compareNewline = [](char i, char j) { return i == '\n' && j == '\n'; };
          str.erase(0, str.find_first_not_of('\n'));
          str.erase(unique(str.begin(), str.end(), compareNewline), str.end());
        }
      macro_output_file << str;
      macro_output_file.close();
    }
}
