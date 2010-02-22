/*
 * Copyright (C) 2010 Dynare Team
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

#include <cstdlib>
#include <cassert>
#include <cerrno>
#include <cmath>
#include <iostream>

#include "ExternalFunctionsTable.hh"
#include "SymbolTable.hh"

ExternalFunctionsTable::ExternalFunctionsTable()
{
};

void
ExternalFunctionsTable::addExternalFunction(const int symb_id, const external_function_options external_function_options_arg)
{
  assert(symb_id >= 0);

  if (external_function_options_arg.secondDerivSymbID > eExtFunNotSet &&
      external_function_options_arg.firstDerivSymbID == eExtFunNotSet)
    {
      cerr << "If the second derivative is provided to the external_function() command,"
           << "the first derivative must also be provided." << endl;
      exit(EXIT_FAILURE);
    }

  if (external_function_options_arg.nargs <= 0)
    {
      cerr << "The number of arguments passed to an external function must be > 0." << endl;
      exit(EXIT_FAILURE);
    }

  external_function_options external_function_options_chng = external_function_options_arg;
  if (external_function_options_arg.firstDerivSymbID  == eExtFunSetButNoNameProvided)
    external_function_options_chng.firstDerivSymbID = symb_id;

  if (external_function_options_arg.secondDerivSymbID == eExtFunSetButNoNameProvided)
    external_function_options_chng.secondDerivSymbID = symb_id;

  externalFunctionTable[symb_id] = external_function_options_chng;
}
