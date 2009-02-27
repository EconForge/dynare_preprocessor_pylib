/*
 * Copyright (C) 2003-2009 Dynare Team
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

#include "CodeInterpreter.hh"

//! Stores the symbol table
/*!
  A symbol is given by its name, and is internally represented by a unique integer.

  When method freeze() is called, computes a distinct sequence of IDs for some types
  (endogenous, exogenous, parameters), which are used by the Matlab/Octave functions.
  We call these "type specific IDs".

  Also manages a TeX name for each symbol, which by default is an empty string.
*/
class SymbolTable
{
private:
  //! Has method freeze() been called?
  bool frozen;

  //! Number of symbols contained in the table
  int size;

  typedef map<string, int> symbol_table_type;
  //! Maps strings to symbol IDs
  symbol_table_type symbol_table;

  //! Maps IDs to names
  vector<string> name_table;
  //! Maps IDs to TeX names
  vector<string> tex_name_table;
  //! Maps IDs to types
  vector<SymbolType> type_table;

  //! Maps symbol IDs to type specific IDs
  vector<int> type_specific_ids;

  //! Maps type specific IDs of endogenous to symbol IDs
  vector<int> endo_ids;
  //! Maps type specific IDs of exogenous to symbol IDs
  vector<int> exo_ids;
  //! Maps type specific IDs of exogenous deterministic to symbol IDs
  vector<int> exo_det_ids;
  //! Maps type specific IDs of parameters to symbol IDs
  vector<int> param_ids;
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
  //! Thrown when trying to access an unknown symbol (by id)
  class UnknownSymbolIDException
  {
  public:
    //! Symbol ID
    int id;
    UnknownSymbolIDException(int id_arg) : id(id_arg) {}
  };
  //! Thrown when trying to access an unknown type specific ID
  class UnknownTypeSpecificIDException
  {
  public:
    int tsid;
    SymbolType type;
    UnknownTypeSpecificIDException(int tsid_arg, SymbolType type_arg) : tsid(tsid_arg), type(type_arg) {}
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
  //! Thrown when table is frozen and trying to modify it
  class FrozenException
  {
  };
  //! Thrown when trying to use the result of freeze() while this method has not yet been called
  class NotYetFrozenException
  {
  };
  //! Add a symbol
  void addSymbol(const string &name, SymbolType type, const string &tex_name = "") throw (AlreadyDeclaredException, FrozenException);
  //! Tests if symbol already exists
  inline bool exists(const string &name) const;
  //! Get symbol name (by ID)
  inline string getName(int id) const throw (UnknownSymbolIDException);
  //! Get TeX name
  inline string getTeXName(int id) const throw (UnknownSymbolIDException);
  //! Get type (by ID)
  inline SymbolType getType(int id) const throw (UnknownSymbolIDException);
  //! Get type (by name)
  inline SymbolType getType(const string &name) const throw (UnknownSymbolNameException);
  //! Get ID (by name)
  inline int getID(const string &name) const throw (UnknownSymbolNameException);
  //! Get ID (by type specific ID)
  int getID(SymbolType type, int tsid) const throw (UnknownTypeSpecificIDException, NotYetFrozenException);
  //! Freeze symbol table
  void freeze() throw (FrozenException);
  //! Change the type of a symbol
  void changeType(int id, SymbolType newtype) throw (UnknownSymbolIDException, FrozenException);
  //! Get type specific ID (by symbol ID)
  inline int getTypeSpecificID(int id) const throw (UnknownSymbolIDException, NotYetFrozenException);
  //! Get type specific ID (by symbol name)
  inline int getTypeSpecificID(const string &name) const throw (UnknownSymbolNameException, NotYetFrozenException);
  //! Get number of endogenous variables
  inline int endo_nbr() const throw (NotYetFrozenException);
  //! Get number of exogenous variables
  inline int exo_nbr() const throw (NotYetFrozenException);
  //! Get number of exogenous deterministic variables
  inline int exo_det_nbr() const throw (NotYetFrozenException);
  //! Get number of parameters
  inline int param_nbr() const throw (NotYetFrozenException);
  //! Returns the greatest symbol ID (the smallest is zero)
  inline int maxID();
  //! Write output of this class
  void writeOutput(ostream &output) const throw (NotYetFrozenException);
};

inline bool
SymbolTable::exists(const string &name) const
{
  symbol_table_type::const_iterator iter = symbol_table.find(name);
  return (iter != symbol_table.end());
}

inline string
SymbolTable::getName(int id) const throw (UnknownSymbolIDException)
{
  if (id < 0 || id >= size)
    throw UnknownSymbolIDException(id);
  else
    return name_table[id];
}

inline string
SymbolTable::getTeXName(int id) const throw (UnknownSymbolIDException)
{
  if (id < 0 || id >= size)
    throw UnknownSymbolIDException(id);
  else
    return tex_name_table[id];
}

inline SymbolType
SymbolTable::getType(int id) const throw (UnknownSymbolIDException)
{
  if (id < 0 || id >= size)
    throw UnknownSymbolIDException(id);
  else
    return type_table[id];
}

inline SymbolType
SymbolTable::getType(const string &name) const throw (UnknownSymbolNameException)
{
  return getType(getID(name));
}

inline int
SymbolTable::getID(const string &name) const throw (UnknownSymbolNameException)
{
  symbol_table_type::const_iterator iter = symbol_table.find(name);
  if (iter != symbol_table.end())
    return iter->second;
  else
    throw UnknownSymbolNameException(name);
}

inline int
SymbolTable::getTypeSpecificID(int id) const throw (UnknownSymbolIDException, NotYetFrozenException)
{
  if (!frozen)
    throw NotYetFrozenException();

  if (id < 0 || id >= size)
    throw UnknownSymbolIDException(id);

  return type_specific_ids[id];
}

inline int
SymbolTable::getTypeSpecificID(const string &name) const throw (UnknownSymbolNameException, NotYetFrozenException)
{
  return getTypeSpecificID(getID(name));
}

inline int
SymbolTable::endo_nbr() const throw (NotYetFrozenException)
{
  if (!frozen)
    throw NotYetFrozenException();

  return endo_ids.size();
}

inline int
SymbolTable::exo_nbr() const throw (NotYetFrozenException)
{
  if (!frozen)
    throw NotYetFrozenException();

  return exo_ids.size();
}

inline int
SymbolTable::exo_det_nbr() const throw (NotYetFrozenException)
{
  if (!frozen)
    throw NotYetFrozenException();

  return exo_det_ids.size();
}

inline int
SymbolTable::param_nbr() const throw (NotYetFrozenException)
{
  if (!frozen)
    throw NotYetFrozenException();

  return param_ids.size();
}

inline int
SymbolTable::maxID()
{
  return(size-1);
}

#endif
