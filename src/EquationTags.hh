/*
 * Copyright © 2020-2022 Dynare Team
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

using namespace std;

class EquationTags
{
private:
  map<int, map<string, string>> eqn_tags;
public:
  class TagNotFoundException
  {
  public:
    const string key, value;
    explicit TagNotFoundException(string key_arg, string value_arg)
      : key{move(key_arg)}, value{move(value_arg)}
    {
    }
  };

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
  getTagsByEqn(const int eqn) const
  {
    if (auto it = eqn_tags.find(eqn); it != eqn_tags.end())
      return it->second;
    return map<string, string>{};
  }

  //! Get equations that have the given key
  set<int> getEqnsByKey(const string &key) const;

  //! Get equations that have the given key and value
  set<int> getEqnsByTag(const string &key, const string &value) const;

  //! Get the first equation that has the given key and value
  int getEqnByTag(const string &key, const string &value) const;

  //! Get the tag value given the equation number and key
  string
  getTagValueByEqnAndKey(int eqn, const string &key) const
  {
    return exists(eqn, key) ? eqn_tags.at(eqn).at(key) : "";
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
    try
      {
        getEqnByTag(key, value);
      }
    catch (TagNotFoundException &e)
      {
        return false;
      }
    return true;
  }

  bool
  exists(const int eqn) const
  {
    return eqn_tags.contains(eqn);
  }

  //! Returns true if equation tag with key exists for a given equation
  bool
  exists(const int eqn, const string &key) const
  {
    return exists(eqn) && eqn_tags.at(eqn).contains(key);
  }

  //! Returns true if equation tag with key and value exists for a given equation
  bool
  exists(const int eqn, const string &key, const string &value) const
  {
    return exists(eqn, key) && eqn_tags.at(eqn).at(key) == value;
  }

  //! Various functions to write equation tags
  void writeCheckSumInfo(ostream &output) const;
  void writeOutput(ostream &output) const;
  void writeLatexOutput(ostream &output, int eqn) const;
  void writeJsonAST(ostream &output, const int eq) const;
};

#endif
