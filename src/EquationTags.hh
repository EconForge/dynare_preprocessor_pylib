/*
 * Copyright © 2020-2023 Dynare Team
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

#ifndef _EQUATION_TAGS_HH
#define _EQUATION_TAGS_HH

#include <map>
#include <set>
#include <string>
#include <optional>

using namespace std;

class EquationTags
{
private:
  map<int, map<string, string>> eqn_tags;
public:
  // Add multiple equation tags for the given equation
  void
  add(int eqn, map<string, string> tags)
  {
    if (eqn_tags.contains(eqn))
      eqn_tags[eqn].insert(move_iterator{tags.begin()}, move_iterator{tags.end()});
    else
      eqn_tags[eqn] = move(tags);
  }

  //! Add a single equation tag for the given equation
  void
  add(int eqn, string key, string value)
  {
    eqn_tags[eqn][move(key)] = move(value);
  }

  //! Clear all equation tag information
  void
  clear()
  {
    eqn_tags.clear();
  }

  //! Erase tags for given equations, using old_eqn_num_2_new as the mapping
  //! to use for the remaining equation numbers
  void erase(const set<int> &eqns, const map<int, int> &old_eqn_num_2_new);

  //! Various functions to get info from equation tags
  //! Get equation tags for a given equation
  map<string, string>
  getTagsByEqn(int eqn) const
  {
    if (auto it = eqn_tags.find(eqn); it != eqn_tags.end())
      return it->second;
    return {};
  }

  //! Get equations that have the given key
  set<int> getEqnsByKey(const string &key) const;

  //! Get equations that have the given key and value
  set<int> getEqnsByTag(const string &key, const string &value) const;

  //! Get the first equation that has the given key and value
  optional<int> getEqnByTag(const string &key, const string &value) const;

  //! Get the tag value given the equation number and key
  optional<string>
  getTagValueByEqnAndKey(int eqn, const string &key) const
  {
    if (auto it = eqn_tags.find(eqn); it != eqn_tags.end())
      if (auto it2 = it->second.find(key); it2 != it->second.end())
        return it2->second;
    return nullopt;
  }

  //! Get the equations marked dynamic
  set<int>
  getDynamicEqns() const
  {
    return getEqnsByTag("dynamic", "");
  }

  //! Returns true if equation tag with key and value exists
  bool
  exists(const string &key, const string &value) const
  {
    return getEqnByTag(key, value).has_value();
  }

  bool
  exists(int eqn) const
  {
    return eqn_tags.contains(eqn);
  }

  //! Returns true if equation tag with key exists for a given equation
  bool
  exists(int eqn, const string &key) const
  {
    auto it = eqn_tags.find(eqn);
    return it != eqn_tags.end() && it->second.contains(key);
  }

  //! Various functions to write equation tags
  void writeCheckSumInfo(ostream &output) const;
  void writeOutput(ostream &output) const;
  void writeLatexOutput(ostream &output, int eqn) const;
  void writeJsonAST(ostream &output, int eq) const;
};

#endif
