/*
 * Copyright (C) 2018 Dynare Team
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

#include <algorithm>

#include "SubModel.hh"

TrendComponentModelTable::TrendComponentModelTable(SymbolTable &symbol_table_arg) :
  symbol_table(symbol_table_arg)
{
}

void
TrendComponentModelTable::addTrendComponentModel(string name_arg,
                                                 vector<string> eqtags_arg,
                                                 vector<string> trend_eqtags_arg)
{
  if (isExistingTrendComponentModelName(name_arg))
    {
      cerr << "Error: a trend component model already exists with the name " << name_arg << endl;
      exit(EXIT_FAILURE);
    }
  eqtags[name_arg] = move(eqtags_arg);
  trend_eqtags[name_arg] = move(trend_eqtags_arg);
  names.insert(move(name_arg));
}

void
TrendComponentModelTable::setEqNums(map<string, vector<int>> eqnums_arg)
{
  eqnums = move(eqnums_arg);
  setUndiffEqnums();
}

void
TrendComponentModelTable::setTrendEqNums(map<string, vector<int>> trend_eqnums_arg)
{
  trend_eqnums = move(trend_eqnums_arg);
  setUndiffEqnums();
}

void
TrendComponentModelTable::setUndiffEqnums()
{
  if (!nontrend_eqnums.empty() || eqnums.empty() || trend_eqnums.empty())
    return;

  for (const auto &it : eqnums)
    {
      vector<int> nontrend_vec;
      for (auto eq : it.second)
        if (find(trend_eqnums[it.first].begin(), trend_eqnums[it.first].end(), eq)
            == trend_eqnums[it.first].end())
          nontrend_vec.push_back(eq);
      nontrend_eqnums[it.first] = nontrend_vec;
    }
}

void
TrendComponentModelTable::setNonstationary(map<string, vector<bool>> nonstationary_arg)
{
  nonstationary = move(nonstationary_arg);
}

void
TrendComponentModelTable::setLhs(map<string, vector<int>> lhs_arg)
{
  lhs = move(lhs_arg);
}

void
TrendComponentModelTable::setRhs(map<string, vector<set<pair<int, int>>>> rhs_arg)
{
  rhs = move(rhs_arg);
}

void
TrendComponentModelTable::setLhsExprT(map<string, vector<expr_t>> lhs_expr_t_arg)
{
  lhs_expr_t = move(lhs_expr_t_arg);
}

void
TrendComponentModelTable::setMaxLags(map<string, vector<int>> max_lags_arg)
{
  max_lags = move(max_lags_arg);
}

void
TrendComponentModelTable::setDiff(map<string, vector<bool>> diff_arg)
{
  diff = move(diff_arg);
}

void
TrendComponentModelTable::setOrigDiffVar(map<string, vector<int>> orig_diff_var_arg)
{
  orig_diff_var = move(orig_diff_var_arg);
}

map<string, vector<string>>
TrendComponentModelTable::getEqTags() const
{
  return eqtags;
}

vector<string>
TrendComponentModelTable::getEqTags(const string &name_arg) const
{
  checkModelName(name_arg);
  return eqtags.find(name_arg)->second;
}

void
TrendComponentModelTable::checkModelName(const string &name_arg) const
{
  if (!isExistingTrendComponentModelName(name_arg))
    {
      cerr << name_arg
           << " is not a recognized equation tag of a trend component model equation" << endl;
      exit(EXIT_FAILURE);
    }
}

vector<bool>
TrendComponentModelTable::getNonstationary(const string &name_arg) const
{
  checkModelName(name_arg);
  return nonstationary.find(name_arg)->second;
}

vector<int>
TrendComponentModelTable::getLhs(const string &name_arg) const
{
  checkModelName(name_arg);
  return lhs.find(name_arg)->second;
}

vector<expr_t>
TrendComponentModelTable::getLhsExprT(const string &name_arg) const
{
  checkModelName(name_arg);
  return lhs_expr_t.find(name_arg)->second;
}

map<string, vector<string>>
TrendComponentModelTable::getTrendEqTags() const
{
  return trend_eqtags;
}

map<string, vector<int>>
TrendComponentModelTable::getEqNums() const
{
  return eqnums;
}

vector<int>
TrendComponentModelTable::getNonTrendEqNums(const string &name_arg) const
{
  checkModelName(name_arg);
  return nontrend_eqnums.find(name_arg)->second;
}

vector<int>
TrendComponentModelTable::getEqNums(const string &name_arg) const
{
  checkModelName(name_arg);
  return eqnums.find(name_arg)->second;
}

vector<int>
TrendComponentModelTable::getMaxLags(const string &name_arg) const
{
  checkModelName(name_arg);
  return max_lags.find(name_arg)->second;
}

int
TrendComponentModelTable::getMaxLag(const string &name_arg) const
{
  int max_lag_int = 0;
  for (auto it : getMaxLags(name_arg))
    max_lag_int = max(max_lag_int, it);
  return max_lag_int;
}

vector<bool>
TrendComponentModelTable::getDiff(const string &name_arg) const
{
  checkModelName(name_arg);
  return diff.find(name_arg)->second;
}

vector<int>
TrendComponentModelTable::getOrigDiffVar(const string &name_arg) const
{
  checkModelName(name_arg);
  return orig_diff_var.find(name_arg)->second;
}

void
TrendComponentModelTable::writeOutput(ostream &output) const
{
  for (const auto &name : names)
    {
      output << "M_.trend_component." << name << ".model_name = '" << name << "';" << endl
             << "M_.trend_component." << name << ".eqtags = {";
      for (const auto &it : eqtags.at(name))
        output << "'" << it << "'; ";
      output << "};" << endl
             << "M_.trend_component." << name << ".eqn = [";
      for (auto it : eqnums.at(name))
        output << it + 1 << " ";
      output << "];" << endl
             << "M_.trend_component." << name << ".trend_eqn = [";
      for (auto it : trend_eqnums.at(name))
        output << it + 1 << " ";
      output << "];" << endl
             << "M_.trend_component." << name << ".trends = [";
      for (auto it : eqnums.at(name))
        if (find(trend_eqnums.at(name).begin(), trend_eqnums.at(name).end(), it)
            == trend_eqnums.at(name).end())
          output << "false ";
        else
          output << "true ";
      output << "];" << endl
             << "M_.trend_component." << name << ".lhs = [";
      for (auto it : lhs.at(name))
        output << symbol_table.getTypeSpecificID(it) + 1 << " ";
      output << "];" << endl
             << "M_.trend_component." << name << ".max_lag = [";
      for (auto it : max_lags.at(name))
        output << it << " ";
      output << "];" << endl
             << "M_.trend_component." << name << ".diff = [";
      for (const auto &it : diff.at(name))
        output << (it ? "true" : "false") << " ";
      output << "];" << endl
             << "M_.trend_component." << name << ".orig_diff_var = [";
      for (auto it : orig_diff_var.at(name))
        output << (it >= 0 ? symbol_table.getTypeSpecificID(it) + 1 : -1) << " ";
      output << "];" << endl;

      int i = 1;
      for (const auto &it : rhs.at(name))
        {
          output << "M_.trend_component." << name << ".rhs.vars_at_eq{" << i << "}.var = [";
          for (const auto &it1 : it)
            output << symbol_table.getTypeSpecificID(it1.first) + 1 << " ";
          output << "];" << endl
                 << "M_.trend_component." << name << ".rhs.vars_at_eq{" << i << "}.lag = [";
          for (const auto &it1 : it)
            output << it1.second << " ";
          output << "];" << endl;

          i++;
        }
    }
}

void
TrendComponentModelTable::writeJsonOutput(ostream &output) const
{
  for (const auto &name : names)
    {
      if (name != *(names.begin()))
        output << ", ";
      output << "{\"statementName\": \"trend_component_model\","
             << "\"model_name\": \"" << name << "\"}";
    }
}
