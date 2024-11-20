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

#ifndef HETEROGENEITY_TABLE_HH
#define HETEROGENEITY_TABLE_HH

#include <map>
#include <ostream>
#include <string>
#include <vector>

using namespace std;

class SymbolTable; // Forward declaration, to avoid circularity

/*
  There is a guarantee that heterogeneity IDs are increasing, i.e. if dimension A is added after
  dimension B, then the ID of A is greater than the ID of B.
  Moreover, the IDs form a contiguous interval starting at 0.
*/
class HeterogeneityTable
{
private:
  // Maps dimension names to IDs
  map<string, int> name_to_id;
  // Maps dimension IDs to names
  vector<string> id_to_name;

  SymbolTable* symbol_table {nullptr}; // Cannot be a reference, because of circularity

  /* Keeps track of the SUM() operator instances.
     Maps a symbol ID that appears inside a SUM() operator into an index in
     M_.heterogeneity_aggregates */
  map<int, int> summed_het_endo_to_index;

  // Maps an index in M_.heterogeneity_aggregates into a symbol ID
  vector<int> index_to_summed_het_endo;

public:
  void setSymbolTable(SymbolTable* symbol_table_arg);

  struct AlreadyDeclaredDimensionException
  {
    // Dimension name
    const string name;
  };
  struct UnknownDimensionNameException
  {
    // Dimension name
    const string name;
  };
  struct UnknownDimensionIDException
  {
    // Dimension ID
    const int id;
  };

  // Returns the dimension ID
  int addDimension(string name);
  [[nodiscard]] bool exists(const string& name) const;
  [[nodiscard]] int getID(const string& name) const;
  [[nodiscard]] string getName(int id) const;
  [[nodiscard]] bool empty() const;
  [[nodiscard]] vector<string> getDimensions() const;
  [[nodiscard]] int size() const;

  struct AlreadyDeclaredSummedHeterogeneousEndogenousException
  {
    const int symb_id;
  };
  struct UnknownSummedHeterogeneousEndogenousException
  {
    const int symb_id;
  };

  void addSummedHeterogeneousEndogenous(int symb_id);
  [[nodiscard]] int getSummedHeterogenousEndogenousIndex(int symb_id) const;
  [[nodiscard]] int aggregateEndoSize() const;

  void writeOutput(ostream& output) const;
  void writeJsonOutput(ostream& output) const;
};

#endif
