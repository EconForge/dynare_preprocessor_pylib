/*
 * Copyright © 2026 Dynare Team
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
#include <iostream>
#include <utility>

#include "DatabaseTable.hh"

int
DatabaseTable::addDatabase(string name)
{
  if (name_to_id.contains(name))
    throw AlreadyDeclaredDatabaseException {move(name)};

  int id {static_cast<int>(id_to_name.size())};
  name_to_id.emplace(name, id);
  id_to_name.push_back(move(name));
  return id;
}

bool
DatabaseTable::exists(const string& name) const
{
  return name_to_id.contains(name);
}

int
DatabaseTable::getID(const string& name) const
{
  if (auto it = name_to_id.find(name); it != name_to_id.end())
    return it->second;
  else
    throw UnknownDatabaseNameException {name};
}

string
DatabaseTable::getName(int id) const
{
  if (id < 0 || id >= static_cast<int>(id_to_name.size()))
    throw UnknownDatabaseIDException {id};
  else
    return id_to_name[id];
}

bool
DatabaseTable::empty() const
{
  return name_to_id.empty();
}

void
DatabaseTable::writeOutput(ostream& output) const
{
  output << "M_.database = {";
  for (bool first_written {false}; const auto& name : id_to_name)
    {
      if (exchange(first_written, true))
        output << ", ";
      output << "'" << name << "'";
    }
  output << "};" << '\n';
}

void
DatabaseTable::writeJsonOutput(ostream& output) const
{
  assert(!empty());

  output << R"("database": [)";
  for (bool first_written {false}; const auto& name : id_to_name)
    {
      if (exchange(first_written, true))
        output << ", ";
      output << '"' << name << '"';
    }
  output << "]" << '\n';
}
