/*
 * Copyright © 2018-2021 Dynare Team
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
 * GNU General Public License for more details.SS
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef _SUBMODEL_HH
#define _SUBMODEL_HH

#include <set>
#include <map>
#include <vector>
#include <iostream>

#include "ExprNode.hh"
#include "SymbolTable.hh"
#include "SymbolList.hh"
#include "Statement.hh"

// DynamicModel.hh can’t be included here, otherwise it would be a circular dependency
class DynamicModel;

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
  map<string, vector<string>> eqtags, target_eqtags;
  map<string, vector<int>> eqnums, target_eqnums, nontarget_eqnums, max_lags, lhs, target_lhs, nontarget_lhs, orig_diff_var;
  map<string, vector<set<pair<int, int>>>> rhs;
  map<string, vector<bool>> diff;
  map<string, vector<expr_t>> lhs_expr_t;
  map<string, vector<int>> target_vars;
  map<string, map<tuple<int, int, int>, expr_t>> AR; // name -> (eqn, lag, lhs_symb_id) -> expr_t
  /* Note that A0 in the trend-component model context is not the same thing as
     in the structural VAR context. */
  map<string, map<tuple<int, int, int>, expr_t>> A0, A0star; // name -> (eqn, lag, col) -> expr_t
public:
  explicit TrendComponentModelTable(SymbolTable &symbol_table_arg);

  //! Add a trend component model
  void addTrendComponentModel(string name_arg, vector<string> eqtags_arg,
                              vector<string> target_eqtags_arg);

  inline bool isExistingTrendComponentModelName(const string &name_arg) const;
  inline bool empty() const;

  const map<string, vector<string>> &getEqTags() const;
  const vector<string> &getEqTags(const string &name_arg) const;
  const map<string, vector<string>> &getTargetEqTags() const;
  const map<string, vector<int>> &getEqNums() const;
  const map<string, vector<int>> &getTargetEqNums() const;
  const vector<int> &getTargetEqNums(const string &name_arg) const;
  const vector<int> &getEqNums(const string &name_arg) const;
  const vector<int> &getMaxLags(const string &name_arg) const;
  int getMaxLag(const string &name_arg) const;
  const vector<int> &getLhs(const string &name_arg) const;
  const vector<expr_t> &getLhsExprT(const string &name_arg) const;
  const vector<bool> &getDiff(const string &name_arg) const;
  const vector<int> &getOrigDiffVar(const string &name_arg) const;
  const map<string, vector<int>> &getNonTargetEqNums() const;
  const vector<int> &getNonTargetEqNums(const string &name_arg) const;
  const vector<int> &getNonTargetLhs(const string &name_arg) const;
  const vector<int> &getTargetLhs(const string &name_arg) const;

  void setVals(map<string, vector<int>> eqnums_arg, map<string, vector<int>> target_eqnums_arg,
               map<string, vector<int>> lhs_arg,
               map<string, vector<expr_t>> lhs_expr_t_arg);
  void setRhs(map<string, vector<set<pair<int, int>>>> rhs_arg);
  void setMaxLags(map<string, vector<int>> max_lags_arg);
  void setDiff(map<string, vector<bool>> diff_arg);
  void setOrigDiffVar(map<string, vector<int>> orig_diff_var_arg);
  void setTargetVar(map<string, vector<int>> target_vars_arg);
  void setAR(map<string, map<tuple<int, int, int>, expr_t>> AR_arg);
  void setA0(map<string, map<tuple<int, int, int>, expr_t>> A0_arg,
             map<string, map<tuple<int, int, int>, expr_t>> A0star_arg);

  //! Write output of this class
  void writeOutput(const string &basename, ostream &output) const;

  //! Write JSON Output
  void writeJsonOutput(ostream &output) const;

private:
  void checkModelName(const string &name_arg) const;
  void setNonTargetEqnums();
};

inline bool
TrendComponentModelTable::isExistingTrendComponentModelName(const string &name_arg) const
{
  return names.find(name_arg) != names.end();
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
  map<string, bool> structural; // Whether VARs are structural or reduced-form
  map<string, vector<string>> eqtags;
  map<string, vector<int>> eqnums, max_lags, lhs, lhs_orig_symb_ids, orig_diff_var;
  map<string, vector<set<pair<int, int>>>> rhs; // name -> for each equation: set of pairs (var, lag)
  map<string, vector<bool>> diff;
  map<string, vector<expr_t>> lhs_expr_t;
  map<string, map<tuple<int, int, int>, expr_t>> AR; // name -> (eqn, lag, lhs_symb_id) -> param_expr_t
  /* The A0 matrix is mainly for structural VARs. For reduced-form VARs, it
     will be equal to the identity matrix. Also note that A0 in the structural
     VAR context is not the same thing as in the trend-component model
     context. */
  map<string, map<tuple<int, int>, expr_t>> A0; // name -> (eqn, lhs_symb_id) -> param_expr_t
  map<string, map<int, expr_t>> constants; // name -> eqn -> constant
public:
  explicit VarModelTable(SymbolTable &symbol_table_arg);

  //! Add a VAR model
  void addVarModel(string name, bool structural_arg, vector<string> eqtags);

  inline bool isExistingVarModelName(const string &name_arg) const;
  inline bool empty() const;

  const map<string, bool> &getStructural() const;
  const map<string, vector<string>> &getEqTags() const;
  const vector<string> &getEqTags(const string &name_arg) const;
  const map<string, vector<int>> &getEqNums() const;
  const vector<bool> &getDiff(const string &name_arg) const;
  const vector<int> &getEqNums(const string &name_arg) const;
  const vector<int> &getMaxLags(const string &name_arg) const;
  int getMaxLag(const string &name_arg) const;
  const vector<int> &getLhs(const string &name_arg) const;
  const vector<int> &getLhsOrigIds(const string &name_arg) const;
  const vector<set<pair<int, int>>> &getRhs(const string &name_arg) const;
  const vector<expr_t> &getLhsExprT(const string &name_arg) const;

  void setEqNums(map<string, vector<int>> eqnums_arg);
  void setLhs(map<string, vector<int>> lhs_arg);
  void setRhs(map<string, vector<set<pair<int, int>>>> rhs_arg);
  void setLhsExprT(map<string, vector<expr_t>> lhs_expr_t_arg);
  void setDiff(map<string, vector<bool>> diff_arg);
  void setMaxLags(map<string, vector<int>> max_lags_arg);
  void setOrigDiffVar(map<string, vector<int>> orig_diff_var_arg);
  void setAR(map<string, map<tuple<int, int, int>, expr_t>> AR_arg);
  void setA0(map<string, map<tuple<int, int>, expr_t>> A0_arg);
  void setConstants(map<string, map<int, expr_t>> constants_arg);

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
  return names.find(name_arg) != names.end();
}

inline bool
VarModelTable::empty() const
{
  return names.empty();
}

class PacModelTable
{
private:
  SymbolTable &symbol_table;
  set<string> names;
  map<string, string> aux_model_name;
  map<string, string> discount;
  // The growth expressions belong to the main dynamic_model from the ModFile instance
  map<string, expr_t> growth, original_growth;
  map<string, vector<tuple<int, int, int, double>>> growth_info;

  /* Stores the name of the PAC equation associated to the model.
     pac_model_name → eq_name */
  map<string, string> eq_name;

  /* Stores symb_ids for alphas created by DynamicModel::addPacModelConsistentExpectationEquation()
     pac_model_name → mce_alpha_symb_ids */
  map<string, vector<int>> mce_alpha_symb_ids;
  /* Stores symb_ids for z1s created by DynamicModel::addPacModelConsistentExpectationEquation()
     pac_model_name → mce_z1_symb_id */
  map<string, int> mce_z1_symb_ids;
  /* Stores symb_ids for h0, h1 parameters
     pac_model_name → parameter symb_ids */
  map<string, vector<int>> h0_indices, h1_indices;
  /* Stores indices for growth neutrality parameters
     pac_model_name → growth_neutrality_param_index */
  map<string, int> growth_neutrality_params;

  // Stores LHS vars (only for backward PAC models)
  map<string, vector<int>> lhs;

  // Stores auxiliary model type (only for backward PAC models)
  map<string, string> aux_model_type;

public:
  /* Stores info about PAC equations
     pac_model_name →
         (lhs, optim_share_index, ar_params_and_vars, ec_params_and_vars, non_optim_vars_params_and_constants, additive_vars_params_and_constants, optim_additive_vars_params_and_constants)
  */
  using equation_info_t = map<string,
                              tuple<pair<int, int>, int, vector<tuple<int, int, int>>, pair<int, vector<tuple<int, bool, int>>>, vector<tuple<int, int, int, double>>, vector<tuple<int, int, int, double>>, vector<tuple<int, int, int, double>>>>;
private:
  equation_info_t equation_info;

  int pacEquationMaxLag(const string &name_arg) const;

public:
  explicit PacModelTable(SymbolTable &symbol_table_arg);
  void addPacModel(string name_arg, string aux_model_name_arg, string discount_arg, expr_t growth_arg);
  bool isExistingPacModelName(const string &name_arg) const;
  bool empty() const;
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings);
  void findDiffNodesInGrowth(lag_equivalence_table_t &diff_nodes) const;
  // Called by DynamicModel::substituteDiff()
  void substituteDiffNodesInGrowth(const lag_equivalence_table_t &diff_nodes, ExprNode::subst_table_t &diff_subst_table, vector<BinaryOpNode *> &neweqs);
  // Must be called after substituteDiffNodesInGrowth()
  void transformPass(ExprNode::subst_table_t &diff_subst_table,
                     DynamicModel &dynamic_model, const VarModelTable &var_model_table,
                     const TrendComponentModelTable &trend_component_model_table);
  void writeOutput(const string &basename, ostream &output) const;
  void writeJsonOutput(ostream &output) const;
};


#endif
