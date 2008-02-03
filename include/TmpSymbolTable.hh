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

#ifndef _TMPSYMBOLTABLE_HH
#define _TMPSYMBOLTABLE_HH

#include <string>
#include <vector>
#include <ostream>

#include "SymbolTable.hh"

/*!
  \class  TmpSymbolTable
  \brief  Defines temparary symbol table used with computing tasks
*/
class TmpSymbolTable
{
private :
  /*! list of string TempSymbolTable */
  std::vector<std::string> tmpsymboltable;
  /*! List of symbol Values */
  std::vector<std::string> nameTable;
  //! A reference to enclosing symbol table
  const SymbolTable &symbol_table;
public :
  /*! Constrcutor */
  TmpSymbolTable(const SymbolTable &symbol_table_arg);
  /*! Destructor*/
  ~TmpSymbolTable();
  /*! Adds a temp symbol */
  void AddTempSymbol(const std::string &symbol);
  /*! Adds a temp symbol and its value */
  void AddTempSymbol(const std::string &symbol1, const std::string &symbol2);
  /*! Write TempSymbolTable to output string */
  void writeOutput(const std::string &varname, std::ostream &output) const;
  //! Clears all content
  void clear();
};
//------------------------------------------------------------------------------
#endif
