/*
 * Copyright © 2003-2022 Dynare Team
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

#ifndef _SYMBOL_LIST_HH
#define _SYMBOL_LIST_HH

#include <string>
#include <vector>
#include <ostream>
#include <algorithm>

#include "WarningConsolidation.hh"
#include "SymbolTable.hh"

using namespace std;

//! Used to store a list of symbols
/*! This class is no more than a vector<string>, with a pretty-printer for Matlab */
class SymbolList
{
private:
  vector<string> symbols;
public:
  SymbolList() = default;
  // This constructor is deliberately not marked explicit, to allow implicit conversion
  SymbolList(vector<string> symbols_arg);

  class SymbolListException
  {
  public:
    const string message;
    SymbolListException(string message_arg) : message{move(message_arg)}
    {
    };
  };
  //! Remove duplicate symbols
  void removeDuplicates(const string &dynare_command, WarningConsolidation &warnings);
  //! Check symbols to ensure variables have been declared and are endogenous
  void checkPass(WarningConsolidation &warnings, const vector<SymbolType> &types, const SymbolTable &symbol_table) const noexcept(false);
  //! Output content in Matlab format
  /*! Creates a string array for Matlab, stored in variable "varname" */
  void writeOutput(const string &varname, ostream &output) const;
  //! Output content in Matlab format without preceding varname of writeOutput
  void write(ostream &output) const;
  //! Write JSON output
  void writeJsonOutput(ostream &output) const;
  //! Is Empty
  bool
  empty() const
  {
    return symbols.empty();
  };
  //! Return the list of symbols
  vector<string> getSymbols() const;
};

#endif
