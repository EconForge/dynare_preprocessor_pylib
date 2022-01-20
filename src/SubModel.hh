/*
 * Copyright © 2018-2022 Dynare Team
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

class VarExpectationModelTable
{
private:
  SymbolTable &symbol_table;
  set<string> names;
  map<string, expr_t> expression;
  map<string, string> aux_model_name;
  map<string, string> horizon;
  map<string, expr_t> discount;
  map<string, int> time_shift;
  // For each model, list of generated auxiliary param ids, in variable-major order
  map<string, vector<int>> aux_param_symb_ids;
  // Decomposition of the expression
  map<string, vector<tuple<int, int, double>>> vars_params_constants;
public:
  explicit VarExpectationModelTable(SymbolTable &symbol_table_arg);
  void addVarExpectationModel(string name_arg, expr_t expression_arg, string aux_model_name_arg,
                              string horizon_arg, expr_t discount_arg, int time_shift_arg);
  bool isExistingVarExpectationModelName(const string &name_arg) const;
  bool empty() const;
  void substituteUnaryOpsInExpression(const lag_equivalence_table_t &nodes, ExprNode::subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs);
  // Called by DynamicModel::substituteDiff()
  void substituteDiffNodesInExpression(const lag_equivalence_table_t &diff_nodes, ExprNode::subst_table_t &diff_subst_table, vector<BinaryOpNode *> &neweqs);
  void transformPass(ExprNode::subst_table_t &diff_subst_table,
                     DynamicModel &dynamic_model, const VarModelTable &var_model_table,
                     const TrendComponentModelTable &trend_component_model_table);
  void writeOutput(const string &basename, ostream &output) const;
  void writeJsonOutput(ostream &output) const;
};

class PacModelTable
{
private:
  SymbolTable &symbol_table;
  set<string> names;
  map<string, string> aux_model_name;
  map<string, string> discount;
  /* The growth expressions belong to the main dynamic_model from the ModFile
     instance. The growth expression is necessarily nullptr for a model with a
     pac_target_info block. */
  map<string, expr_t> growth, original_growth;
  /* Information about the structure of growth expressions (which must be a
     linear combination of variables).
     Each tuple represents a term: (endo_id, lag, param_id, constant) */
  using growth_info_t = vector<tuple<int, int, int, double>>;
  map<string, growth_info_t> growth_info;
  // The “auxname” option of pac_model (empty if not passed)
  map<string, string> auxname;
  // The “kind” option of pac_model (“undefined” if not passed)
  map<string, PacTargetKind> kind;

  /* Stores the name of the PAC equation associated to the model.
     pac_model_name → eq_name */
  map<string, string> eq_name;

  /* Stores symb_ids for auxiliary endogenous created for the expression
     substituted to the pac_expectation operator:
     - in the backward case, this auxiliary contains exactly the
     pac_expectation value
     - in the MCE case, this auxiliary represents Z₁ (i.e. without the growth
     correction term)
     Note that this structure is not used in the presence of the
     pac_target_info block.
     pac_model_name → symb_id */
  map<string, int> aux_var_symb_ids;
  /* Stores symb_ids for auxiliary parameters created for the expression
     substituted to the pac_expectation operator (excluding the growth
     neutrality correction):
     - in the backward case, contains the “h” parameters
     - in the MCE case, contains the “α” parameters
     Note that this structure is not used in the presence of the
     pac_target_info block.
     pac_model_name → symb_ids */
  map<string, vector<int>> aux_param_symb_ids;
  /* Stores indices for growth neutrality parameters
     pac_model_name → growth_neutrality_param_index.
     This map is not used for PAC models with a pac_target_info block. */
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

public:
  /* (component variable/expr, growth, auxname, kind, coeff. in the linear
     combination, growth_param ID possibly equal to -1, vector of h parameters,
     original_growth, growth_info) */
  using target_component_t = tuple<expr_t, expr_t, string, PacTargetKind, expr_t, int, vector<int>, expr_t, growth_info_t>;

private:
  // pac_model_name → (target variable/expr, auxname_target_nonstationary, target components)
  map<string, tuple<expr_t, string, vector<target_component_t>>> target_info;

  int pacEquationMaxLag(const string &name_arg) const;

  // Return a text representation of a kind (but fails on “unspecified” kind value)
  static string kindToString(PacTargetKind kind);

public:
  explicit PacModelTable(SymbolTable &symbol_table_arg);
  void addPacModel(string name_arg, string aux_model_name_arg, string discount_arg, expr_t growth_arg, string auxname_arg, PacTargetKind kind_arg);
  bool isExistingPacModelName(const string &name_arg) const;
  bool empty() const;
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings);
  // Called by DynamicModel::substituteUnaryOps()
  void substituteUnaryOpsInGrowth(const lag_equivalence_table_t &nodes, ExprNode::subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs);
  void findDiffNodesInGrowth(lag_equivalence_table_t &diff_nodes) const;
  // Called by DynamicModel::substituteDiff()
  void substituteDiffNodesInGrowth(const lag_equivalence_table_t &diff_nodes, ExprNode::subst_table_t &diff_subst_table, vector<BinaryOpNode *> &neweqs);
  // Must be called after substituteDiffNodesInGrowth() and substituteUnaryOpsInGrowth()
  void transformPass(const lag_equivalence_table_t &unary_ops_nodes,
                     ExprNode::subst_table_t &unary_ops_subst_table,
                     const lag_equivalence_table_t &diff_nodes,
                     ExprNode::subst_table_t &diff_subst_table,
                     DynamicModel &dynamic_model, const VarModelTable &var_model_table,
                     const TrendComponentModelTable &trend_component_model_table);
  void writeOutput(const string &basename, ostream &output) const;
  void writeJsonOutput(ostream &output) const;
  void setTargetExpr(const string &name_arg, expr_t target);
  void setTargetAuxnameNonstationary(const string &name_arg, string auxname);
  /* Only the first four elements of the tuple are expected to be set by the
     caller. The other ones will be filled by this class. */
  void addTargetComponent(const string &name_arg, target_component_t component);
  void writeTargetCoefficientsFile(const string &basename) const;
};


#endif
