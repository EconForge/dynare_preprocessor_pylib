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

#include <iostream>
#include <fstream>

#include "MacroDriver.hh"

MacroDriver::MacroDriver() : trace_scanning(false), trace_parsing(false)
{
}

MacroDriver::~MacroDriver()
{
}

void
MacroDriver::parse(const string &f, ostream &out)
{
  file = f;
  out_stream = &out;

  ifstream in(f.c_str(), ios::binary);

  lexer = new MacroFlex(&in, &out);
  lexer->set_debug(trace_scanning);

  Macro::parser parser(*this, out);
  parser.parse();

  delete lexer;
}

void
MacroDriver::error(const Macro::parser::location_type &l, const string &m)
{
  cerr << "ERROR: " << l << ": " << m << endl;
  exit(-1);
}
