/*
 * Copyright © 2020-2021 Dynare Team
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

#include "EquationTags.hh"

#include <regex>

set<int>
EquationTags::getEqnsByKey(const string &key) const
{
  set<int> retval;
  for (const auto & [eqn, tags] : eqn_tags)
    if (tags.find(key) != tags.end())
      retval.insert(eqn);
  return retval;
}

set<int>
EquationTags::getEqnsByTag(const string &key, const string &value) const
{
  set<int> retval;
  for (const auto & [eqn, tags] : eqn_tags)
    if (auto tmp = tags.find(key); tmp != tags.end() && tmp->second == value)
      retval.insert(eqn);
  return retval;
}

int
EquationTags::getEqnByTag(const string &key, const string &value) const
{
  for (const auto & [eqn, tags] : eqn_tags)
    if (auto tmp = tags.find(key); tmp != tags.end() && tmp->second == value)
      return eqn;
  return -1;
}

void
EquationTags::erase(const set<int> &eqns, const map<int, int> &old_eqn_num_2_new)
{
  for (const auto &eqn : eqns)
    eqn_tags.erase(eqn);

  for (const auto & [oldeqn, neweqn] : old_eqn_num_2_new)
    for (auto & [eqn, tags] : eqn_tags)
      if (eqn == oldeqn)
        {
          auto tmp = eqn_tags.extract(eqn);
          tmp.key() = neweqn;
          eqn_tags.insert(move(tmp));
        }
}

void
EquationTags::writeCheckSumInfo(ostream &output) const
{
  for (const auto & [eqn, tags] : eqn_tags)
    for (const auto & [key, value] : tags)
      output << "  " << eqn + 1
             << key << " " << value << endl;
}

void
EquationTags::writeOutput(ostream &output) const
{
  output << "M_.equations_tags = {" << endl;
  for (const auto & [eqn, tags] : eqn_tags)
    for (const auto & [key, value] : tags)
      output << "  " << eqn + 1 << " , '"
             << key << "' , '" << value << "' ;" << endl;
  output << "};" << endl;
}

void
EquationTags::writeOccbinOutput(ostream &output) const
{
  map<int, map<string, string>> occbin_options;
  for (const auto & [eqn, tags] : eqn_tags)
    for (const auto & [key, value] : tags)
      if (key == "pswitch"
          || key == "bind"
          || key == "relax"
          || key == "pcrit")
        occbin_options[eqn][key] = value;

  int idx = 0;
  for (const auto & [eqn, tags] : occbin_options)
    {
      idx++;
      for (const auto & [key, value] : tags)
        output << "M_.occbin.constraint(" << idx << ")."
               << key << " = '" << value << "';" << endl;
      output << "M_.occbin.constraint(" << idx << ").equation = "
             << eqn+1 << ";" << endl;
    }
}

void
EquationTags::writeLatexOutput(ostream &output, int eqn) const
{
  if (!exists(eqn))
    return;

  auto escape_special_latex_symbols
    = [](string str)
      {
        const regex special_latex_chars (R"([&%$#_{}])");
        const regex backslash (R"(\\)");
        const regex tilde (R"(~)");
        const regex carrot (R"(\^)");
        const regex textbackslash (R"(\\textbackslash)");
        str = regex_replace(str, backslash, R"(\textbackslash)");
        str = regex_replace(str, special_latex_chars, R"(\$&)");
        str = regex_replace(str, carrot, R"(\^{})");
        str = regex_replace(str, tilde, R"(\textasciitilde{})");
        return regex_replace(str, textbackslash, R"(\textbackslash{})");
      };

  bool wrote_eq_tag = false;
  output << R"(\noindent[)";
  for (const auto & [key, value] : eqn_tags.at(eqn))
    {
      if (wrote_eq_tag)
        output << ", ";
      output << escape_special_latex_symbols(key);

      if (!value.empty())
        output << "= `" << escape_special_latex_symbols(value) << "'";

      wrote_eq_tag = true;
    }
  output << "]" << endl;
}

void
EquationTags::writeJsonAST(ostream &output, const int eqn) const
{
  if (!exists(eqn))
    return;

  output << R"(, "tags": {)";
  bool wroteFirst = false;
  for (const auto &[key, value] : eqn_tags.at(eqn))
    {
      if (wroteFirst)
        output << ", ";
      else
        wroteFirst = true;
      output << R"(")" << key << R"(": ")" << value << R"(")";
    }
  output << "}";
}

bool
EquationTags::hasOccbinTags() const
{
  for (const auto & [eqn, tags] : eqn_tags)
    for (const auto & [key, value] : tags)
      if (key == "pswitch"
          || key == "bind"
          || key == "relax"
          || key == "pcrit")
        return true;
  return false;
}
