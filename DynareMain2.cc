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

using namespace std;

#include <iostream>

#include "ParsingDriver.hh"
#include "ModFile.hh"

void
main2(stringstream &in, string &basename, bool debug, bool clear_all)
{
  ParsingDriver p;

  // Do parsing and construct internal representation of mod file
  ModFile *mod_file = p.parse(in, debug);

  // Run checking pass
  mod_file->checkPass();

  // Do computations
  mod_file->computingPass();

  // Write outputs
  mod_file->writeOutputFiles(basename, clear_all);

  delete mod_file;

  cout << "Preprocessing completed." << endl
       << "Starting Matlab computing ..." << endl;
}
