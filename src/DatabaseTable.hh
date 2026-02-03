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

#ifndef DATABASE_TABLE_HH
#define DATABASE_TABLE_HH

#include <map>
#include <string>
#include <vector>

using namespace std;

// Stores the databases declared with the “database” command
class DatabaseTable
{
public:
  int addDatabase(string name);
  [[nodiscard]] bool exists(const string& name) const;
  [[nodiscard]] int getID(const string& name) const;
  [[nodiscard]] string getName(int id) const;
  [[nodiscard]] bool empty() const;
  void writeOutput(ostream& output) const;
  void writeJsonOutput(ostream& output) const;

  struct AlreadyDeclaredDatabaseException
  {
    const string name;
  };
  struct UnknownDatabaseNameException
  {
    const string name;
  };
  struct UnknownDatabaseIDException
  {
    const int id;
  };

private:
  map<string, int> name_to_id;
  vector<string> id_to_name;
};

#endif
