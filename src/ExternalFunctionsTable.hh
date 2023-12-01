/*
 * Copyright © 2010-2023 Dynare Team
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

#ifndef EXTERNAL_FUNCTIONS_TABLE_HH
#define EXTERNAL_FUNCTIONS_TABLE_HH

#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

//! Handles external functions
class ExternalFunctionsTable
{
public:
  //! Thrown when trying to access an unknown symbol (by id)
  struct UnknownExternalFunctionSymbolIDException
  {
    //! Symbol ID
    int id;
  };

  /* For all arguments, -2 means not set
   * For firstDerivSymbID and secondDerivSymbID, -1 means that the derivatives are
   * provided in the main function given in the by the "name" option to the
   * external_function() statement.
   */
  struct external_function_options
  {
    int nargs, firstDerivSymbID, secondDerivSymbID;
  };
  using external_function_table_type = map<int, external_function_options>;
  //! Symbol ID used when no external function exists that calculates the derivative
  constexpr static int IDNotSet = -1;
  //! Symbol ID used when the derivative is obtained from the top-level function
  constexpr static int IDSetButNoNameProvided = -2;
  //! Default number of arguments when nargs is not specified
  constexpr static int defaultNargs = 1;

private:
  //! Map containing options provided to external_functions()
  external_function_table_type externalFunctionTable;

public:
  //! Adds an external function to the table as well as its derivative functions
  void addExternalFunction(int symb_id,
                           const external_function_options& external_function_options_arg,
                           bool track_nargs);
  //! See if the function exists in the External Functions Table
  [[nodiscard]] inline bool exists(int symb_id) const;
  //! Get the number of arguments for a given external function
  [[nodiscard]] inline int getNargs(int symb_id) const noexcept(false);
  //! Get the symbol_id of the first derivative function
  [[nodiscard]] inline int getFirstDerivSymbID(int symb_id) const noexcept(false);
  //! Get the symbol_id of the second derivative function
  [[nodiscard]] inline int getSecondDerivSymbID(int symb_id) const noexcept(false);
};

inline bool
ExternalFunctionsTable::exists(int symb_id) const
{
  return externalFunctionTable.contains(symb_id);
}

inline int
ExternalFunctionsTable::getNargs(int symb_id) const noexcept(false)
{
  if (auto it = externalFunctionTable.find(symb_id); it != externalFunctionTable.end())
    return it->second.nargs;
  else
    throw UnknownExternalFunctionSymbolIDException {symb_id};
}

inline int
ExternalFunctionsTable::getFirstDerivSymbID(int symb_id) const noexcept(false)
{
  if (auto it = externalFunctionTable.find(symb_id); it != externalFunctionTable.end())
    return it->second.firstDerivSymbID;
  else
    throw UnknownExternalFunctionSymbolIDException {symb_id};
}

inline int
ExternalFunctionsTable::getSecondDerivSymbID(int symb_id) const noexcept(false)
{
  if (auto it = externalFunctionTable.find(symb_id); it != externalFunctionTable.end())
    return it->second.secondDerivSymbID;
  else
    throw UnknownExternalFunctionSymbolIDException {symb_id};
}

#endif
