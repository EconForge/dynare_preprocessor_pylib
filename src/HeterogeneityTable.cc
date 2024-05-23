/*
 * Copyright © 2024 Dynare Team
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

#include <cassert>
#include <utility>

#include "HeterogeneityTable.hh"
#include "SymbolTable.hh"

void
HeterogeneityTable::setSymbolTable(SymbolTable* symbol_table_arg)
{
  symbol_table = symbol_table_arg;
}

int
HeterogeneityTable::addDimension(string name)
{
  if (name_to_id.contains(name))
    throw AlreadyDeclaredDimensionException {move(name)};

  int id {static_cast<int>(id_to_name.size())};
  name_to_id.emplace(name, id);
  id_to_name.push_back(move(name));
  return id;
}

bool
HeterogeneityTable::exists(const string& name) const
{
  return name_to_id.contains(name);
}

int
HeterogeneityTable::getID(const string& name) const
{
  if (auto it = name_to_id.find(name); it != name_to_id.end())
    return it->second;
  else
    throw UnknownDimensionNameException {name};
}

string
HeterogeneityTable::getName(int id) const
{
  if (id < 0 || id >= static_cast<int>(id_to_name.size()))
    throw UnknownDimensionIDException {id};
  else
    return id_to_name[id];
}

bool
HeterogeneityTable::empty() const
{
  return name_to_id.empty();
}

vector<string>
HeterogeneityTable::getDimensions() const
{
  return id_to_name;
}

int
HeterogeneityTable::size() const
{
  return static_cast<int>(name_to_id.size());
}

void
HeterogeneityTable::addSummedHeterogeneousEndogenous(int symb_id)
{
  assert(symbol_table->getType(symb_id) == SymbolType::heterogeneousEndogenous);
  if (summed_het_endo_to_index.contains(symb_id))
    throw AlreadyDeclaredSummedHeterogeneousEndogenousException {symb_id};

  int index {static_cast<int>(index_to_summed_het_endo.size())};
  summed_het_endo_to_index.emplace(symb_id, index);
  index_to_summed_het_endo.push_back(symb_id);
}

int
HeterogeneityTable::getSummedHeterogenousEndogenousIndex(int symb_id) const
{
  if (auto it = summed_het_endo_to_index.find(symb_id); it != summed_het_endo_to_index.end())
    return it->second;
  else
    throw UnknownSummedHeterogeneousEndogenousException {symb_id};
}

int
HeterogeneityTable::aggregateEndoSize() const
{
  return index_to_summed_het_endo.size();
}

void
HeterogeneityTable::writeOutput(ostream& output) const
{
  for (size_t id {0}; id < id_to_name.size(); id++)
    output << "M_.heterogeneity(" << id + 1 << ").dimension_name = '" << id_to_name[id] << "';"
           << endl;

  output << "M_.heterogeneity_aggregates = {" << endl;
  for (int symb_id : index_to_summed_het_endo)
    output << "'sum', " << symbol_table->getHeterogeneityDimension(symb_id) + 1 << ", "
           << symbol_table->getTypeSpecificID(symb_id) + 1 << ";" << endl;
  output << "};" << endl;
}

void
HeterogeneityTable::writeJsonOutput(ostream& output) const
{
  assert(!empty());
  output << R"("heterogeneity_dimension": [)";
  for (bool first_written {false}; const auto& dim : id_to_name)
    {
      if (exchange(first_written, true))
        output << ", ";
      output << '"' << dim << '"';
    }
  output << "]" << endl;
}
