/*
 * Copyright (C) 2003-2008 Dynare Team
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

using namespace std;

#include <iostream>
#include <sstream>
#include <fstream>

#include <cctype>    // for tolower()
#include <algorithm> // for transform()

#include "macro/MacroDriver.hh"

/* Prototype for second part of main function
   Splitting main() in two parts was necessary because ParsingDriver.h and MacroDriver.h can't be
   included simultaneously (because of Bison limitations).
*/
void main2(stringstream &in, string &basename, bool debug, bool clear_all, bool no_tmp_terms);

void
usage()
{
  cerr << "Dynare usage: dynare mod_file [debug] [noclearall] [savemacro] [notmpterms]" << endl;
  exit(EXIT_FAILURE);
}

int
main(int argc, char** argv)
{
  if (argc < 2)
    {
      cerr << "Missing model file!" << endl;
      usage();
    }

  bool clear_all = true;
  bool save_macro = false;
  bool debug = false;
  bool no_tmp_terms = false;

  // Parse options
  for (int arg = 2; arg < argc; arg++)
    {
      if (string(argv[arg]) == string("debug"))
        debug = true;
      else if (string(argv[arg]) == string("noclearall"))
        clear_all = false;
      else if (string(argv[arg]) == string("savemacro"))
        save_macro = true;
      else if (string(argv[arg]) == string("notmpterms"))
        no_tmp_terms = true;
      else
        {
          cerr << "Unknown option: " << argv[arg] << endl;
          usage();
        }
    }

  cout << "Starting Dynare ..." << endl
       << "Starting preprocessing of the model file ..." << endl;

  // Construct basename (check file extension is correct then remove it)
  string basename = argv[1];
  string ext = basename.substr(basename.size() - 4);
  transform(ext.begin(), ext.end(), ext.begin(), (int(*)(int)) tolower); // Convert ext to lowercase
  if (ext != string(".mod") && ext != string(".dyn"))
    {
      cerr << "mod_file extension must be .mod or .dyn!" << endl;
      exit(EXIT_FAILURE);
    }
  basename.erase(basename.size() - 4, 4);

  // Do macro processing
  MacroDriver m;

  stringstream macro_output;
  m.parse(argv[1], macro_output, debug);
  if (save_macro)
    {
      ofstream macro_output_file((basename + "-macroexp.mod").c_str());
      macro_output_file << macro_output.str();
      macro_output_file.close();
    }

  // Do the rest
  main2(macro_output, basename, debug, clear_all, no_tmp_terms);

  return EXIT_SUCCESS;
}
