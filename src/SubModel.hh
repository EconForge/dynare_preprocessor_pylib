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

#ifndef _SUBMODEL_HH
#define _SUBMODEL_HH

#include <set>
#include <map>
#include <vector>
#include <iostream>

#include "ExprNode.hh"
#include "SymbolTable.hh"

using namespace std;

//! A table with all Trend Component Models in the .mod file
/*!
  The unique name of the trend component model is the identifier
*/
class TrendComponentModelTable
{
private:
  SymbolTable &symbol_table;
  set<string> names;
  map<string, vector<string>> eqtags, trend_eqtags;
  map<string, vector<int>> eqnums, trend_eqnums, nontrend_eqnums, max_lags, lhs, orig_diff_var;
  map<string, vector<set<pair<int, int>>>> rhs;
  map<string, vector<bool>> diff, nonstationary;
  map<string, vector<expr_t>> lhs_expr_t;
  map<string, vector<int>> trend_vars;
public:
  TrendComponentModelTable(SymbolTable &symbol_table_arg);

  //! Add a trend component model
  void addTrendComponentModel(string name_arg, vector<string> eqtags_arg,
                              vector<string> trend_eqtags_arg);

  inline bool isExistingTrendComponentModelName(const string &name_arg) const;
  inline bool empty() const;

  map<string, vector<string>> getEqTags() const;
  vector<string> getEqTags(const string &name_arg) const;
  map<string, vector<string>> getTrendEqTags() const;
  map<string, vector<int>> getEqNums() const;
  map<string, vector<int>> getTrendEqNums() const;
  vector<int> getEqNums(const string &name_arg) const;
  vector<int> getMaxLags(const string &name_arg) const;
  int getMaxLag(const string &name_arg) const;
  vector<int> getLhs(const string &name_arg) const;
  vector<expr_t> getLhsExprT(const string &name_arg) const;
  vector<bool> getDiff(const string &name_arg) const;
  vector<int> getOrigDiffVar(const string &name_arg) const;
  vector<int> getNonTrendEqNums(const string &name_arg) const;
  vector<bool> getNonstationary(const string &name_arg) const;

  void setEqNums(map<string, vector<int>> eqnums_arg);
  void setTrendEqNums(map<string, vector<int>> trend_eqnums_arg);
  void setLhs(map<string, vector<int>> lhs_arg);
  void setRhs(map<string, vector<set<pair<int, int>>>> rhs_arg);
  void setLhsExprT(map<string, vector<expr_t>> lhs_expr_t_arg);
  void setMaxLags(map<string, vector<int>> max_lags_arg);
  void setDiff(map<string, vector<bool>> diff_arg);
  void setOrigDiffVar(map<string, vector<int>> orig_diff_var_arg);
  void setNonstationary(map<string, vector<bool>> nonstationary_arg);
  void setTrendVar(map<string, vector<int>> trend_vars_arg);

  //! Write output of this class
  void writeOutput(ostream &output) const;

  //! Write JSON Output
  void writeJsonOutput(ostream &output) const;

private:
  void checkModelName(const string &name_arg) const;
  void setUndiffEqnums();
};

inline bool
TrendComponentModelTable::isExistingTrendComponentModelName(const string &name_arg) const
{
  return names.find(name_arg) == names.end() ? false : true;
}

inline bool
TrendComponentModelTable::empty() const
{
  return names.empty();
}


class VarModelTable
{
private:
  SymbolTable &symbol_table;
  set<string> names;
  map<string, pair<SymbolList, int>> symbol_list_and_order;
  map<string, vector<string>> eqtags;
  map<string, vector<int>> eqnums, max_lags, lhs, orig_diff_var;
  map<string, vector<set<pair<int, int>>>> rhs;
  map<string, vector<bool>> diff, nonstationary;
  map<string, vector<expr_t>> lhs_expr_t;
  map<string, map<tuple<int, int, int>, int>> AR; // AR: name -> (eqn, lag, lhs_symb_id) -> param_symb_id
public:
  VarModelTable(SymbolTable &symbol_table_arg);

  //! Add a trend component model
  void addVarModel(string name, vector<string> eqtags,
                   pair<SymbolList, int> symbol_list_and_order_arg);

  inline bool isExistingVarModelName(const string &name_arg) const;
  inline bool empty() const;

  map<string, vector<string>> getEqTags() const;
  vector<string> getEqTags(const string &name_arg) const;
  map<string, vector<int>> getEqNums() const;
  vector<int> getEqNums(const string &name_arg) const;
  vector<int> getMaxLags(const string &name_arg) const;
  int getMaxLag(const string &name_arg) const;
  vector<int> getLhs(const string &name_arg) const;
  vector<bool> getNonstationary(const string &name_arg) const;
  map<string, pair<SymbolList, int>> getSymbolListAndOrder() const;
  vector<set<pair<int, int>>> getRhs(const string &name_arg) const;
  vector<expr_t> getLhsExprT(const string &name_arg) const;

  void setEqNums(map<string, vector<int>> eqnums_arg);
  void setLhs(map<string, vector<int>> lhs_arg);
  void setRhs(map<string, vector<set<pair<int, int>>>> rhs_arg);
  void setLhsExprT(map<string, vector<expr_t>> lhs_expr_t_arg);
  void setNonstationary(map<string, vector<bool>> nonstationary_arg);
  void setDiff(map<string, vector<bool>> diff_arg);
  void setMaxLags(map<string, vector<int>> max_lags_arg);
  void setOrigDiffVar(map<string, vector<int>> orig_diff_var_arg);
  void setAR(map<string, map<tuple<int, int, int>, int>> AR_arg);

  //! Write output of this class
  void writeOutput(const string &basename, ostream &output) const;

  //! Write JSON Output
  void writeJsonOutput(ostream &output) const;

private:
  void checkModelName(const string &name_arg) const;
};

inline bool
VarModelTable::isExistingVarModelName(const string &name_arg) const
{
  return names.find(name_arg) == names.end() ? false : true;
}

inline bool
VarModelTable::empty() const
{
  return names.empty();
}

#endif
