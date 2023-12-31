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

#include <regex>

#include "SymbolList.hh"

SymbolList::SymbolList(vector<string> symbols_arg) :
  symbols{move(symbols_arg)}
{
}

void
SymbolList::checkPass(WarningConsolidation &warnings, const vector<SymbolType> &types,
                      const SymbolTable &symbol_table) const noexcept(false)
{
  if (types.empty())
    return;

  smatch m;
  string regex_str = "AUX_EXPECT_|MULT_";
  for (auto type : types)
    if (type == SymbolType::endogenous)
      {
        regex_str += "|AUX_ENDO_|LOG_";
        break;
      }
  regex re("^(" + regex_str +")");
  for (const auto &symbol : symbols)
    {
      if (!symbol_table.exists(symbol))
        {
          if (regex_search(symbol, m, re))
            {
              warnings << "WARNING: symbol_list variable " << symbol << " has not yet been declared. "
                       << "This is being ignored because the variable name corresponds to a possible "
                       << "auxiliary variable name." << endl;
              return;
            }
          else
            throw SymbolListException{"Variable " + symbol +  " was not declared."};
        }

      if (none_of(types.begin(), types.end(),
                  [&](SymbolType type) { return symbol_table.getType(symbol) == type; }))
        {
          string valid_types;
          for (auto type : types)
            switch(type)
              {
              case SymbolType::endogenous:
                valid_types += "endogenous, ";
                break;
              case SymbolType::exogenous:
                valid_types += "exogenous, ";
                break;
              case SymbolType::epilogue:
                valid_types += "epilogue, ";
                break;
              case SymbolType::parameter:
                valid_types += "parameter, ";
                break;
              case SymbolType::exogenousDet:
                valid_types += "exogenousDet, ";
                break;
              case SymbolType::trend:
                valid_types += "trend, ";
                break;
              case SymbolType::logTrend:
                valid_types += "logTrend, ";
                break;
              case SymbolType::modFileLocalVariable:
                valid_types += "modFileLocalVariable, ";
                break;
              case SymbolType::modelLocalVariable:
                valid_types += "modelLocalVariable, ";
                break;
              case SymbolType::externalFunction:
                valid_types += "externalFunction, ";
                break;
              case SymbolType::statementDeclaredVariable:
                valid_types += "statementDeclaredVariable, ";
                break;
              case SymbolType::unusedEndogenous:
                valid_types += "unusedEndogenous, ";
                break;
              case SymbolType::excludedVariable:
                valid_types += "excludedVariable, ";
                break;
              }
          valid_types = valid_types.erase(valid_types.size()-2, 2);
          throw SymbolListException{"Variable " + symbol +  " is not one of {" + valid_types + "}"};
        }
    }
}

void
SymbolList::writeOutput(const string &varname, ostream &output) const
{
  output << varname << " = {";
  for (bool printed_something{false};
       const auto &name : symbols)
    {
      if (exchange(printed_something, true))
        output << ";";
      output << "'" << name << "'";
    }
  output << "};" << endl;
}

void
SymbolList::writeJsonOutput(ostream &output) const
{
  output << R"("symbol_list": [)";
  for (bool printed_something{false};
       const auto &name : symbols)
    {
      if (exchange(printed_something, true))
        output << ",";
      output << R"(")" << name << R"(")";
    }
  output << "]";
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
  for (const auto &it : symbols)
    if (find(unique_symbols.begin(), unique_symbols.end(), it) == unique_symbols.end())
      unique_symbols.push_back(it);
    else
      warnings << "WARNING: In " << dynare_command << ": " << it
               << " found more than once in symbol list. Removing all but first occurrence." << endl;
  symbols = unique_symbols;
}
