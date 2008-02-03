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

using namespace std;

#include "SymbolTable.hh"
#include "TmpSymbolTable.hh"
#include "Interface.hh"

TmpSymbolTable::TmpSymbolTable(const SymbolTable &symbol_table_arg) :
  symbol_table(symbol_table_arg)
{
}

TmpSymbolTable::~TmpSymbolTable()
{
}

void
TmpSymbolTable::AddTempSymbol(const string &symbol)
{
  // FIXME: add check to verify that symbol exists in symbol_table
  // FIXME: add check to verify that symbol doesn't yet exist in the present table
  tmpsymboltable.push_back(symbol);
}

void
TmpSymbolTable::AddTempSymbol(const string &symbol1, const string &symbol2)
{
  // FIXME: add checks to verify that symbol1 and symbol2 exist in symbol_table
  // FIXME: add check to verify that symbol1 doesn't yet exist in the present table
  tmpsymboltable.push_back(symbol1);
  nameTable.push_back(symbol2);
}

void
TmpSymbolTable::writeOutput(const string &varname, ostream &output) const
{
  output << varname << "=[];" << endl;
  for (vector<string>::const_iterator it = tmpsymboltable.begin();
       it != tmpsymboltable.end(); it++)
    {
      output << varname << " = ";
      output << interfaces::strvcat(varname, "'" + *it + "'") << ";" << endl;
    }
}

void
TmpSymbolTable::clear()
{
  tmpsymboltable.clear();
  nameTable.clear();
}
