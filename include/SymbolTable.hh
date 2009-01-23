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

#ifndef _SYMBOLTABLE_HH
#define _SYMBOLTABLE_HH

using namespace std;

#include <map>
#include <string>
#include <vector>
#include <ostream>
#include <iostream>
#include "CodeInterpreter.hh"

//! Stores the symbol table
/*!
  A symbol is given by its name, and is internally represented by a pair (type, id).

  There is a distinct sequence of ids for each type, so two symbol of different types can have the same id.

  Also manages a TeX name for each symbol, which by default is an empty string.
*/
class SymbolTable
{
private:
  //! A symbol is represented by a pair (type, id)
  typedef pair<SymbolType, int> named_symbol_type;

  typedef map<string, named_symbol_type> symbol_table_type;
  //! Maps strings to pairs (type,id)
  symbol_table_type symbol_table;

  typedef map<named_symbol_type, string> inv_symbol_table_type;
  //! Maps pairs (type, id) to names
  inv_symbol_table_type name_table;
  //! Maps pairs (type, id) to TeX names
  inv_symbol_table_type tex_name_table;
public:
  SymbolTable();
  //! Thrown when trying to access an unknown symbol (by name)
  class UnknownSymbolNameException
  {
  public:
    //! Symbol name
    string name;
    UnknownSymbolNameException(const string &name_arg) : name(name_arg) {}
  };
  //! Thrown when trying to access an unknown symbol (by type+id pair)
  class UnknownSymbolIDException
  {
  public:
    //! Symbol type
    SymbolType type;
    //! Symbol ID
    int id;
    UnknownSymbolIDException(SymbolType type_arg, int id_arg) : type(type_arg), id(id_arg) {}
  };
  //! Thrown when trying to declare a symbol twice
  class AlreadyDeclaredException
  {
  public:
    //! Symbol name
    string name;
    //! Was the previous declaration done with the same symbol type ?
    bool same_type;
    AlreadyDeclaredException(const string &name_arg, bool same_type_arg) : name(name_arg), same_type(same_type_arg) {}
  };
  //! Number of declared endogenous variables
  int endo_nbr;
  //! Number of declared exogenous variables
  int exo_nbr;
  //! Number of declared deterministic exogenous variables
  int exo_det_nbr;
  //! Number of declared parameters
  int parameter_nbr;
  //! Number of model local variables
  int model_local_variable_nbr;
  //! Number of modfile local variables
  int modfile_local_variable_nbr;
  //! Number of unknown functions
  int unknown_function_nbr;
  //! Add a symbol
  void addSymbol(const string &name, SymbolType type, const string &tex_name = "") throw (AlreadyDeclaredException);
  //! Tests if symbol already exists
  inline bool exists(const string &name) const;
  //! Get symbol name by type and ID
  inline string getNameByID(SymbolType type, int id) const throw (UnknownSymbolIDException);
  //! Get TeX name by type and ID
  inline string getTeXNameByID(SymbolType type, int id) const throw (UnknownSymbolIDException);
  //! Get type by name
  inline SymbolType getType(const string &name) const throw (UnknownSymbolNameException);
  //! Get ID by name
  inline int getID(const string &name) const throw (UnknownSymbolNameException);
  //! Write output of this class
  void writeOutput(ostream &output) const;
};

inline bool
SymbolTable::exists(const string &name) const
{
  symbol_table_type::const_iterator iter = symbol_table.find(name);
  return (iter != symbol_table.end());
}

inline string
SymbolTable::getNameByID(SymbolType type, int id) const throw (UnknownSymbolIDException)
{
  inv_symbol_table_type::const_iterator iter = name_table.find(make_pair(type, id));
  if (iter != name_table.end())
    return iter->second;
  else
    throw UnknownSymbolIDException(type, id);
}

inline string
SymbolTable::getTeXNameByID(SymbolType type, int id) const throw (UnknownSymbolIDException)
{
  inv_symbol_table_type::const_iterator iter = tex_name_table.find(make_pair(type, id));
  if (iter != tex_name_table.end())
    return iter->second;
  else
    throw UnknownSymbolIDException(type, id);
}

inline SymbolType
SymbolTable::getType(const string &name) const throw (UnknownSymbolNameException)
{
  symbol_table_type::const_iterator iter = symbol_table.find(name);
  if (iter != symbol_table.end())
    return iter->second.first;
  else
    throw UnknownSymbolNameException(name);
}

inline int
SymbolTable::getID(const string &name) const throw (UnknownSymbolNameException)
{
  symbol_table_type::const_iterator iter = symbol_table.find(name);
  if (iter != symbol_table.end())
    return iter->second.second;
  else
    throw UnknownSymbolNameException(name);
}

#endif
