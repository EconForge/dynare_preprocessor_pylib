/*
 * Copyright © 2003-2019 Dynare Team
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

#ifndef _SYMBOL_LIST_HH
#define _SYMBOL_LIST_HH

#include <string>
#include <vector>
#include <ostream>

#include "WarningConsolidation.hh"

using namespace std;

//! Used to store a list of symbols
/*! This class is no more than a vector<string>, with a pretty-printer for Matlab */
class SymbolList
{
private:
  //! Internal container for symbol list
  vector<string> symbols;
public:
  //! Adds a symbol to the list
  void addSymbol(const string &symbol);
  //! Removed duplicate symbols
  void removeDuplicates(const string &dynare_command, WarningConsolidation &warnings);
  //! Output content in Matlab format
  /*! Creates a string array for Matlab, stored in variable "varname" */
  void writeOutput(const string &varname, ostream &output) const;
  //! Output content in Matlab format without preceding varname of writeOutput
  void write(ostream &output) const;
  //! Write JSON output
  void writeJsonOutput(ostream &output) const;
  //! Clears all content
  void clear();
  //! Get a copy of the string vector
  vector<string>
  get_symbols() const
  {
    return symbols;
  };
  //! Is Empty
  int
  empty() const
  {
    return symbols.empty();
  };
  //! Return the number of Symbols contained in the list
  int getSize() const;
  //! Return the list of symbols
  vector<string> getSymbols() const;
};

#endif
