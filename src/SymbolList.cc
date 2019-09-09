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

#include "SymbolList.hh"

void
SymbolList::addSymbol(const string &symbol)
{
  symbols.push_back(symbol);
}

void
SymbolList::writeOutput(const string &varname, ostream &output) const
{
  output << varname << " = {";
  for (auto it = symbols.begin();
       it != symbols.end(); ++it)
    {
      if (it != symbols.begin())
        output << ";";
      output << "'" << *it << "'";
    }
  output << "};" << endl;
}

void
SymbolList::writeJsonOutput(ostream &output) const
{
  output << R"("symbol_list": [)";
  for (auto it = symbols.begin();
       it != symbols.end(); ++it)
    {
      if (it != symbols.begin())
        output << ",";
      output << R"(")" << *it << R"(")";
    }
  output << "]";
}

void
SymbolList::clear()
{
  symbols.clear();
}

int
SymbolList::getSize() const
{
  return symbols.size();
}

vector<string>
SymbolList::getSymbols() const
{
  return symbols;
}

void
SymbolList::removeDuplicates(const string &dynare_command, WarningConsolidation &warnings)
{
  vector<string> unique_symbols;
  for (auto & it : symbols)
    if (find(unique_symbols.begin(), unique_symbols.end(), it) == unique_symbols.end())
      unique_symbols.push_back(it);
    else
      warnings << "WARNING: In " << dynare_command << ": " << it
               << " found more than once in symbol list. Removing all but first occurence." << endl;
  symbols = unique_symbols;
}
