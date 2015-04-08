/*
 * Copyright (C) 2010-2015 Dynare Team
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

#ifndef _EXTERNALFUNCTIONSTABLE_HH
#define _EXTERNALFUNCTIONSTABLE_HH

using namespace std;

#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "SymbolTable.hh"

enum ExternalFunctionSetOrNot
  {
    eExtFunSetButNoNameProvided = -2, //! Signifies that the derivative is obtained from the top-level function
    eExtFunNotSet = -1,               //! Signifies that no external function exists that calculates the derivative
    eExtFunSetDefaultNargs = 1        //! This is the default number of arguments when nargs is not specified
  };

//! Handles external functions
class ExternalFunctionsTable
{
public:
  //! Thrown when trying to access an unknown symbol (by id)
  class UnknownExternalFunctionSymbolIDException
  {
  public:
    //! Symbol ID
    int id;
    UnknownExternalFunctionSymbolIDException(int id_arg) : id(id_arg)
    {
    }
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
  typedef map<int, external_function_options> external_function_table_type;
private:
  //! Map containing options provided to external_functions()
  external_function_table_type externalFunctionTable;
public:
  ExternalFunctionsTable();
  //! Adds an external function to the table as well as its derivative functions
  void addExternalFunction(int symb_id, const external_function_options &external_function_options_arg, bool track_nargs);
  //! See if the function exists in the External Functions Table
  inline bool exists(int symb_id) const;
  //! Get the number of arguments for a given external function
  inline int getNargs(int symb_id) const throw (UnknownExternalFunctionSymbolIDException);
  //! Get the symbol_id of the first derivative function
  inline int getFirstDerivSymbID(int symb_id) const throw (UnknownExternalFunctionSymbolIDException);
  //! Get the symbol_id of the second derivative function
  inline int getSecondDerivSymbID(int symb_id) const throw (UnknownExternalFunctionSymbolIDException);
  //! Returns the total number of unique external functions declared or used in the .mod file
  inline int get_total_number_of_unique_model_block_external_functions() const;
  //! Reindex external function table (after removal of extra exo)
  void reindex(SymbolTable &new_symbol_table, SymbolTable &orig_symbol_table);
};

inline bool
ExternalFunctionsTable::exists(int symb_id) const
{
  external_function_table_type::const_iterator iter = externalFunctionTable.find(symb_id);
  return (iter != externalFunctionTable.end());
}

inline int
ExternalFunctionsTable::getNargs(int symb_id) const throw (UnknownExternalFunctionSymbolIDException)
{
  if (exists(symb_id))
    return externalFunctionTable.find(symb_id)->second.nargs;
  else
    throw UnknownExternalFunctionSymbolIDException(symb_id);
}

inline int
ExternalFunctionsTable::getFirstDerivSymbID(int symb_id) const throw (UnknownExternalFunctionSymbolIDException)
{
  if (exists(symb_id))
    return externalFunctionTable.find(symb_id)->second.firstDerivSymbID;
  else
    throw UnknownExternalFunctionSymbolIDException(symb_id);
}

inline int
ExternalFunctionsTable::getSecondDerivSymbID(int symb_id) const throw (UnknownExternalFunctionSymbolIDException)
{
  if (exists(symb_id))
    return externalFunctionTable.find(symb_id)->second.secondDerivSymbID;
  else
    throw UnknownExternalFunctionSymbolIDException(symb_id);
}

inline int
ExternalFunctionsTable::get_total_number_of_unique_model_block_external_functions() const
{
  int number_of_unique_model_block_external_functions = 0;
  for (external_function_table_type::const_iterator it = externalFunctionTable.begin();
       it != externalFunctionTable.end(); it++)
    if (it->second.nargs > 0)
      number_of_unique_model_block_external_functions++;

  return number_of_unique_model_block_external_functions;
}

#endif
