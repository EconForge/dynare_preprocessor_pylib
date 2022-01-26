/*
 * Copyright © 2007-2022 Dynare Team
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

#include <iostream>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <utility>
#include <limits>

#include "ExprNode.hh"
#include "DataTree.hh"
#include "ModFile.hh"

ExprNode::ExprNode(DataTree &datatree_arg, int idx_arg) : datatree{datatree_arg}, idx{idx_arg}
{
}

expr_t
ExprNode::getDerivative(int deriv_id)
{
  if (!preparedForDerivation)
    prepareForDerivation();

  // Return zero if derivative is necessarily null (using symbolic a priori)
  if (auto it = non_null_derivatives.find(deriv_id); it == non_null_derivatives.end())
    return datatree.Zero;

  // If derivative is stored in cache, use the cached value, otherwise compute it (and cache it)
  if (auto it2 = derivatives.find(deriv_id); it2 != derivatives.end())
    return it2->second;
  else
    {
      expr_t d = computeDerivative(deriv_id);
      derivatives[deriv_id] = d;
      return d;
    }
}

int
ExprNode::precedence(ExprNodeOutputType output_type, const temporary_terms_t &temporary_terms) const
{
  // For a constant, a variable, or a unary op, the precedence is maximal
  return 100;
}

int
ExprNode::precedenceJson(const temporary_terms_t &temporary_terms) const
{
  // For a constant, a variable, or a unary op, the precedence is maximal
  return 100;
}

int
ExprNode::cost(int cost, bool is_matlab) const
{
  // For a terminal node, the cost is null
  return 0;
}

int
ExprNode::cost(const vector<vector<temporary_terms_t>> &blocks_temporary_terms, bool is_matlab) const
{
  // For a terminal node, the cost is null
  return 0;
}

int
ExprNode::cost(const map<pair<int, int>, temporary_terms_t> &temp_terms_map, bool is_matlab) const
{
  // For a terminal node, the cost is null
  return 0;
}

bool
ExprNode::checkIfTemporaryTermThenWrite(ostream &output, ExprNodeOutputType output_type,
                                        const temporary_terms_t &temporary_terms,
                                        const temporary_terms_idxs_t &temporary_terms_idxs) const
{
  if (auto it = temporary_terms.find(const_cast<ExprNode *>(this)); it == temporary_terms.end())
    return false;

  auto it2 = temporary_terms_idxs.find(const_cast<ExprNode *>(this));
  // It is the responsibility of the caller to ensure that all temporary terms have their index
  assert(it2 != temporary_terms_idxs.end());
  output << "T" << LEFT_ARRAY_SUBSCRIPT(output_type)
         << it2->second + ARRAY_SUBSCRIPT_OFFSET(output_type)
         << RIGHT_ARRAY_SUBSCRIPT(output_type);

  return true;
}

pair<expr_t, int>
ExprNode::getLagEquivalenceClass() const
{
  int index = maxLead();

  if (index == numeric_limits<int>::min())
    index = 0; // If no variable in the expression, the equivalence class has size 1

  return { decreaseLeadsLags(index), index };
}

void
ExprNode::collectVariables(SymbolType type, set<int> &result) const
{
  set<pair<int, int>> symbs_lags;
  collectDynamicVariables(type, symbs_lags);
  transform(symbs_lags.begin(), symbs_lags.end(), inserter(result, result.begin()),
            [](auto x) { return x.first; });
}

void
ExprNode::collectEndogenous(set<pair<int, int>> &result) const
{
  set<pair<int, int>> symb_ids_and_lags;
  collectDynamicVariables(SymbolType::endogenous, symb_ids_and_lags);
  for (const auto &[symb_id, lag] : symb_ids_and_lags)
    result.emplace(datatree.symbol_table.getTypeSpecificID(symb_id), lag);
}

void
ExprNode::computeTemporaryTerms(const pair<int, int> &derivOrder,
                                map<pair<int, int>, temporary_terms_t> &temp_terms_map,
                                map<expr_t, pair<int, pair<int, int>>> &reference_count,
                                bool is_matlab) const
{
  // Nothing to do for a terminal node
}

void
ExprNode::computeBlockTemporaryTerms(int blk, int eq, vector<vector<temporary_terms_t>> &blocks_temporary_terms,
                                     map<expr_t, tuple<int, int, int>> &reference_count) const
{
  // Nothing to do for a terminal node
}

void
ExprNode::writeOutput(ostream &output) const
{
  writeOutput(output, ExprNodeOutputType::matlabOutsideModel, {}, {});
}

void
ExprNode::writeOutput(ostream &output, ExprNodeOutputType output_type) const
{
  writeOutput(output, output_type, {}, {});
}

void
ExprNode::writeOutput(ostream &output, ExprNodeOutputType output_type, const temporary_terms_t &temporary_terms, const temporary_terms_idxs_t &temporary_terms_idxs) const
{
  writeOutput(output, output_type, temporary_terms, temporary_terms_idxs, {});
}

void
ExprNode::compile(ostream &CompileCode, unsigned int &instruction_number,
                  bool lhs_rhs, const temporary_terms_t &temporary_terms,
                  const temporary_terms_idxs_t &temporary_terms_idxs, bool dynamic, bool steady_dynamic) const
{
  compile(CompileCode, instruction_number, lhs_rhs, temporary_terms, temporary_terms_idxs, dynamic, steady_dynamic, {});
}

void
ExprNode::writeExternalFunctionOutput(ostream &output, ExprNodeOutputType output_type,
                                      const temporary_terms_t &temporary_terms,
                                      const temporary_terms_idxs_t &temporary_terms_idxs,
                                      deriv_node_temp_terms_t &tef_terms) const
{
  // Nothing to do
}

void
ExprNode::writeJsonExternalFunctionOutput(vector<string> &efout,
                                          const temporary_terms_t &temporary_terms,
                                          deriv_node_temp_terms_t &tef_terms,
                                          bool isdynamic) const
{
  // Nothing to do
}

void
ExprNode::compileExternalFunctionOutput(ostream &CompileCode, unsigned int &instruction_number,
                                        bool lhs_rhs, const temporary_terms_t &temporary_terms,
                                        const temporary_terms_idxs_t &temporary_terms_idxs, bool dynamic, bool steady_dynamic,
                                        deriv_node_temp_terms_t &tef_terms) const
{
  // Nothing to do
}

VariableNode *
ExprNode::createEndoLeadAuxiliaryVarForMyself(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  int n = maxEndoLead();
  assert(n >= 2);

  if (auto it = subst_table.find(this);
      it != subst_table.end())
    return const_cast<VariableNode *>(it->second);

  expr_t substexpr = decreaseLeadsLags(n-1);
  int lag = n-2;

  // Each iteration tries to create an auxvar such that auxvar(+1)=expr(-lag)
  // At the beginning (resp. end) of each iteration, substexpr is an expression (possibly an auxvar) equivalent to expr(-lag-1) (resp. expr(-lag))
  while (lag >= 0)
    {
      expr_t orig_expr = decreaseLeadsLags(lag);
      if (auto it = subst_table.find(orig_expr); it == subst_table.end())
        {
          int symb_id = datatree.symbol_table.addEndoLeadAuxiliaryVar(orig_expr->idx, substexpr);
          neweqs.push_back(datatree.AddEqual(datatree.AddVariable(symb_id, 0), substexpr));
          substexpr = datatree.AddVariable(symb_id, +1);
          assert(dynamic_cast<VariableNode *>(substexpr));
          subst_table[orig_expr] = dynamic_cast<VariableNode *>(substexpr);
        }
      else
        substexpr = const_cast<VariableNode *>(it->second);

      lag--;
    }

  return dynamic_cast<VariableNode *>(substexpr);
}

VariableNode *
ExprNode::createExoLeadAuxiliaryVarForMyself(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  int n = maxExoLead();
  assert(n >= 1);

  if (auto it = subst_table.find(this);
      it != subst_table.end())
    return const_cast<VariableNode *>(it->second);

  expr_t substexpr = decreaseLeadsLags(n);
  int lag = n-1;

  // Each iteration tries to create an auxvar such that auxvar(+1)=expr(-lag)
  // At the beginning (resp. end) of each iteration, substexpr is an expression (possibly an auxvar) equivalent to expr(-lag-1) (resp. expr(-lag))
  while (lag >= 0)
    {
      expr_t orig_expr = decreaseLeadsLags(lag);
      if (auto it = subst_table.find(orig_expr); it == subst_table.end())
        {
          int symb_id = datatree.symbol_table.addExoLeadAuxiliaryVar(orig_expr->idx, substexpr);
          neweqs.push_back(datatree.AddEqual(datatree.AddVariable(symb_id, 0), substexpr));
          substexpr = datatree.AddVariable(symb_id, +1);
          assert(dynamic_cast<VariableNode *>(substexpr));
          subst_table[orig_expr] = dynamic_cast<VariableNode *>(substexpr);
        }
      else
        substexpr = const_cast<VariableNode *>(it->second);

      lag--;
    }

  return dynamic_cast<VariableNode *>(substexpr);
}

bool
ExprNode::isNumConstNodeEqualTo(double value) const
{
  return false;
}

bool
ExprNode::isVariableNodeEqualTo(SymbolType type_arg, int variable_id, int lag_arg) const
{
  return false;
}

void
ExprNode::fillErrorCorrectionRow(int eqn,
                                 const vector<int> &nontarget_lhs,
                                 const vector<int> &target_lhs,
                                 map<tuple<int, int, int>, expr_t> &A0,
                                 map<tuple<int, int, int>, expr_t> &A0star) const
{
  vector<pair<expr_t, int>> terms;
  decomposeAdditiveTerms(terms, 1);

  for (const auto &it : terms)
    {
      pair<int, vector<tuple<int, int, int, double>>> m;
      try
        {
          m = it.first->matchParamTimesLinearCombinationOfVariables();
          for (auto &t : m.second)
            get<3>(t) *= it.second; // Update sign of constants
        }
      catch (MatchFailureException &e)
        {
          /* FIXME: we should not just skip them, but rather verify that they are
             autoregressive terms or residuals (probably by merging the two "fill" procedures) */
          continue;
        }

      // Helper function
      auto one_step_orig = [this](int symb_id) {
                             return datatree.symbol_table.isAuxiliaryVariable(symb_id) ?
                               datatree.symbol_table.getOrigSymbIdForDiffAuxVar(symb_id) : symb_id;
                           };

      /* Verify that all variables belong to the error-correction term.
         FIXME: same remark as above about skipping terms. */
      bool not_ec = false;
      for (const auto &t : m.second)
        {
          int vid = one_step_orig(get<0>(t));
          not_ec = not_ec || (find(target_lhs.begin(), target_lhs.end(), vid) == target_lhs.end()
                              && find(nontarget_lhs.begin(), nontarget_lhs.end(), vid) == nontarget_lhs.end());
        }
      if (not_ec)
        continue;

      // Now fill the matrices
      for (auto [var_id, lag, param_id, constant] : m.second)
        {
          int orig_vid = one_step_orig(var_id);
          int orig_lag = datatree.symbol_table.isAuxiliaryVariable(var_id) ? -datatree.symbol_table.getOrigLeadLagForDiffAuxVar(var_id) : lag;
          if (find(target_lhs.begin(), target_lhs.end(), orig_vid) == target_lhs.end())
            {
              // This an LHS variable, so fill A0
              if (constant != 1)
                {
                  cerr << "ERROR in trend component model: LHS variable should not appear with a multiplicative constant in error correction term" << endl;
                  exit(EXIT_FAILURE);
                }
              if (param_id != -1)
                {
                  cerr << "ERROR in trend component model: spurious parameter in error correction term" << endl;
                  exit(EXIT_FAILURE);
                }
              int colidx = static_cast<int>(distance(nontarget_lhs.begin(), find(nontarget_lhs.begin(), nontarget_lhs.end(), orig_vid)));
              if (A0.find({eqn, -orig_lag, colidx}) != A0.end())
                {
                  cerr << "ExprNode::fillErrorCorrection: Error filling A0 matrix: "
                       << "lag/symb_id encountered more than once in equation" << endl;
                  exit(EXIT_FAILURE);
                }
              A0[{eqn, -orig_lag, colidx}] = datatree.AddVariable(m.first);
            }
          else
            {
              // This is a target, so fill A0star
              int colidx = static_cast<int>(distance(target_lhs.begin(), find(target_lhs.begin(), target_lhs.end(), orig_vid)));
              expr_t e = datatree.AddTimes(datatree.AddVariable(m.first), datatree.AddPossiblyNegativeConstant(-constant));
              if (param_id != -1)
                e = datatree.AddTimes(e, datatree.AddVariable(param_id));
              if (auto coor = tuple(eqn, -orig_lag, colidx); A0star.find(coor) == A0star.end())
                A0star[coor] = e;
              else
                A0star[coor] = datatree.AddPlus(e, A0star[coor]);
            }
        }
    }
}

void
ExprNode::matchMatchedMoment(vector<int> &symb_ids, vector<int> &lags, vector<int> &powers) const
{
  throw MatchFailureException{"Unsupported expression"};
}

bool
ExprNode::isConstant() const
{
  set<pair<int, int>> symbs_lags;
  collectDynamicVariables(SymbolType::endogenous, symbs_lags);
  collectDynamicVariables(SymbolType::exogenous, symbs_lags);
  collectDynamicVariables(SymbolType::exogenousDet, symbs_lags);
  return symbs_lags.empty();
}

bool
ExprNode::hasExogenous() const
{
  set<pair<int, int>> symbs_lags;
  collectDynamicVariables(SymbolType::exogenous, symbs_lags);
  collectDynamicVariables(SymbolType::exogenousDet, symbs_lags);
  return !symbs_lags.empty();
}


NumConstNode::NumConstNode(DataTree &datatree_arg, int idx_arg, int id_arg) :
  ExprNode{datatree_arg, idx_arg},
  id{id_arg}
{
}

int
NumConstNode::countDiffs() const
{
  return 0;
}

void
NumConstNode::prepareForDerivation()
{
  preparedForDerivation = true;
  // All derivatives are null, so non_null_derivatives is left empty
}

expr_t
NumConstNode::computeDerivative(int deriv_id)
{
  return datatree.Zero;
}

void
NumConstNode::writeOutput(ostream &output, ExprNodeOutputType output_type,
                          const temporary_terms_t &temporary_terms,
                          const temporary_terms_idxs_t &temporary_terms_idxs,
                          const deriv_node_temp_terms_t &tef_terms) const
{
  if (!checkIfTemporaryTermThenWrite(output, output_type, temporary_terms, temporary_terms_idxs))
    output << datatree.num_constants.get(id);
}

void
NumConstNode::writeJsonAST(ostream &output) const
{
  output << R"({"node_type" : "NumConstNode", "value" : )";
  output << std::stof(datatree.num_constants.get(id)) << "}";
}

void
NumConstNode::writeJsonOutput(ostream &output,
                              const temporary_terms_t &temporary_terms,
                              const deriv_node_temp_terms_t &tef_terms,
                              bool isdynamic) const
{
  output << datatree.num_constants.get(id);
}

bool
NumConstNode::containsExternalFunction() const
{
  return false;
}

double
NumConstNode::eval(const eval_context_t &eval_context) const noexcept(false)
{
  return datatree.num_constants.getDouble(id);
}

void
NumConstNode::compile(ostream &CompileCode, unsigned int &instruction_number,
                      bool lhs_rhs, const temporary_terms_t &temporary_terms,
                      const temporary_terms_idxs_t &temporary_terms_idxs, bool dynamic, bool steady_dynamic,
                      const deriv_node_temp_terms_t &tef_terms) const
{
  FLDC_ fldc(datatree.num_constants.getDouble(id));
  fldc.write(CompileCode, instruction_number);
}

void
NumConstNode::collectVARLHSVariable(set<expr_t> &result) const
{
  cerr << "ERROR: you can only have variables or unary ops on LHS of VAR" << endl;
  exit(EXIT_FAILURE);
}

void
NumConstNode::collectDynamicVariables(SymbolType type_arg, set<pair<int, int>> &result) const
{
}

void
NumConstNode::computeSubExprContainingVariable(int symb_id, int lag, set<expr_t> &contain_var) const
{
}

BinaryOpNode *
NumConstNode::normalizeEquationHelper(const set<expr_t> &contain_var, expr_t rhs) const
{
  cerr << "NumConstNode::normalizeEquation: this should not happen" << endl;
  exit(EXIT_FAILURE);
}

expr_t
NumConstNode::getChainRuleDerivative(int deriv_id, const map<int, BinaryOpNode *> &recursive_variables)
{
  return datatree.Zero;
}

expr_t
NumConstNode::toStatic(DataTree &static_datatree) const
{
  return static_datatree.AddNonNegativeConstant(datatree.num_constants.get(id));
}

void
NumConstNode::computeXrefs(EquationInfo &ei) const
{
}

expr_t
NumConstNode::clone(DataTree &datatree) const
{
  return datatree.AddNonNegativeConstant(datatree.num_constants.get(id));
}

int
NumConstNode::maxEndoLead() const
{
  return 0;
}

int
NumConstNode::maxExoLead() const
{
  return 0;
}

int
NumConstNode::maxEndoLag() const
{
  return 0;
}

int
NumConstNode::maxExoLag() const
{
  return 0;
}

int
NumConstNode::maxLead() const
{
  return numeric_limits<int>::min();
}

int
NumConstNode::maxLag() const
{
  return numeric_limits<int>::min();
}

int
NumConstNode::maxLagWithDiffsExpanded() const
{
  return numeric_limits<int>::min();
}

expr_t
NumConstNode::undiff() const
{
  return const_cast<NumConstNode *>(this);
}

int
NumConstNode::VarMaxLag(const set<expr_t> &lhs_lag_equiv) const
{
  return 0;
}

expr_t
NumConstNode::decreaseLeadsLags(int n) const
{
  return const_cast<NumConstNode *>(this);
}

expr_t
NumConstNode::decreaseLeadsLagsPredeterminedVariables() const
{
  return const_cast<NumConstNode *>(this);
}

expr_t
NumConstNode::substituteEndoLeadGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const
{
  return const_cast<NumConstNode *>(this);
}

expr_t
NumConstNode::substituteEndoLagGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  return const_cast<NumConstNode *>(this);
}

expr_t
NumConstNode::substituteExoLead(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const
{
  return const_cast<NumConstNode *>(this);
}

expr_t
NumConstNode::substituteExoLag(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  return const_cast<NumConstNode *>(this);
}

expr_t
NumConstNode::substituteExpectation(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool partial_information_model) const
{
  return const_cast<NumConstNode *>(this);
}

expr_t
NumConstNode::substituteAdl() const
{
  return const_cast<NumConstNode *>(this);
}

expr_t
NumConstNode::substituteModelLocalVariables() const
{
  return const_cast<NumConstNode *>(this);
}

expr_t
NumConstNode::substituteVarExpectation(const map<string, expr_t> &subst_table) const
{
  return const_cast<NumConstNode *>(this);
}

void
NumConstNode::findDiffNodes(lag_equivalence_table_t &nodes) const
{
}

void
NumConstNode::findUnaryOpNodesForAuxVarCreation(lag_equivalence_table_t &nodes) const
{
}

int
NumConstNode::findTargetVariable(int lhs_symb_id) const
{
  return -1;
}

expr_t
NumConstNode::substituteDiff(const lag_equivalence_table_t &nodes, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  return const_cast<NumConstNode *>(this);
}

expr_t
NumConstNode::substituteUnaryOpNodes(const lag_equivalence_table_t &nodes, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  return const_cast<NumConstNode *>(this);
}

expr_t
NumConstNode::substitutePacExpectation(const string &name, expr_t subexpr)
{
  return const_cast<NumConstNode *>(this);
}

expr_t
NumConstNode::substitutePacTargetNonstationary(const string &name, expr_t subexpr)
{
  return const_cast<NumConstNode *>(this);
}

expr_t
NumConstNode::differentiateForwardVars(const vector<string> &subset, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  return const_cast<NumConstNode *>(this);
}

bool
NumConstNode::isNumConstNodeEqualTo(double value) const
{
  if (datatree.num_constants.getDouble(id) == value)
    return true;
  else
    return false;
}

bool
NumConstNode::isVariableNodeEqualTo(SymbolType type_arg, int variable_id, int lag_arg) const
{
  return false;
}

bool
NumConstNode::containsPacExpectation(const string &pac_model_name) const
{
  return false;
}

bool
NumConstNode::containsPacTargetNonstationary(const string &pac_model_name) const
{
  return false;
}

expr_t
NumConstNode::replaceTrendVar() const
{
  return const_cast<NumConstNode *>(this);
}

expr_t
NumConstNode::detrend(int symb_id, bool log_trend, expr_t trend) const
{
  return const_cast<NumConstNode *>(this);
}

expr_t
NumConstNode::removeTrendLeadLag(const map<int, expr_t> &trend_symbols_map) const
{
  return const_cast<NumConstNode *>(this);
}

bool
NumConstNode::isInStaticForm() const
{
  return true;
}

bool
NumConstNode::isParamTimesEndogExpr() const
{
  return false;
}

expr_t
NumConstNode::substituteStaticAuxiliaryVariable() const
{
  return const_cast<NumConstNode *>(this);
}

expr_t
NumConstNode::replaceVarsInEquation(map<VariableNode *, NumConstNode *> &table) const
{
  return const_cast<NumConstNode *>(this);
}

VariableNode::VariableNode(DataTree &datatree_arg, int idx_arg, int symb_id_arg, int lag_arg) :
  ExprNode{datatree_arg, idx_arg},
  symb_id{symb_id_arg},
  lag{lag_arg}
{
  // It makes sense to allow a lead/lag on parameters: during steady state calibration, endogenous and parameters can be swapped
  assert(get_type() != SymbolType::externalFunction
         && (lag == 0 || (get_type() != SymbolType::modelLocalVariable && get_type() != SymbolType::modFileLocalVariable)));
}

void
VariableNode::prepareForDerivation()
{
  if (preparedForDerivation)
    return;

  preparedForDerivation = true;

  // Fill in non_null_derivatives
  switch (get_type())
    {
    case SymbolType::endogenous:
    case SymbolType::exogenous:
    case SymbolType::exogenousDet:
    case SymbolType::parameter:
    case SymbolType::trend:
    case SymbolType::logTrend:
      // For a variable or a parameter, the only non-null derivative is with respect to itself
      non_null_derivatives.insert(datatree.getDerivID(symb_id, lag));
      break;
    case SymbolType::modelLocalVariable:
      datatree.getLocalVariable(symb_id)->prepareForDerivation();
      // Non null derivatives are those of the value of the local parameter
      non_null_derivatives = datatree.getLocalVariable(symb_id)->non_null_derivatives;
      break;
    case SymbolType::modFileLocalVariable:
    case SymbolType::statementDeclaredVariable:
    case SymbolType::unusedEndogenous:
      // Such a variable is never derived
      break;
    case SymbolType::externalFunction:
    case SymbolType::epilogue:
      cerr << "VariableNode::prepareForDerivation: impossible case" << endl;
      exit(EXIT_FAILURE);
    case SymbolType::excludedVariable:
      cerr << "VariableNode::prepareForDerivation: impossible case: "
           << "You are trying to derive a variable that has been excluded via model_remove/var_remove/include_eqs/exclude_eqs: "
           << datatree.symbol_table.getName(symb_id) << endl;
      exit(EXIT_FAILURE);
    }
}

expr_t
VariableNode::computeDerivative(int deriv_id)
{
  switch (get_type())
    {
    case SymbolType::endogenous:
    case SymbolType::exogenous:
    case SymbolType::exogenousDet:
    case SymbolType::parameter:
    case SymbolType::trend:
    case SymbolType::logTrend:
      if (deriv_id == datatree.getDerivID(symb_id, lag))
        return datatree.One;
      else
        return datatree.Zero;
    case SymbolType::modelLocalVariable:
      return datatree.getLocalVariable(symb_id)->getDerivative(deriv_id);
    case SymbolType::modFileLocalVariable:
      cerr << "modFileLocalVariable is not derivable" << endl;
      exit(EXIT_FAILURE);
    case SymbolType::statementDeclaredVariable:
      cerr << "statementDeclaredVariable is not derivable" << endl;
      exit(EXIT_FAILURE);
    case SymbolType::unusedEndogenous:
      cerr << "unusedEndogenous is not derivable" << endl;
      exit(EXIT_FAILURE);
    case SymbolType::externalFunction:
    case SymbolType::epilogue:
    case SymbolType::excludedVariable:
      cerr << "VariableNode::computeDerivative: Impossible case!" << endl;
      exit(EXIT_FAILURE);
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

bool
VariableNode::containsExternalFunction() const
{
  if (get_type() == SymbolType::modelLocalVariable)
    return datatree.getLocalVariable(symb_id)->containsExternalFunction();

  return false;
}

void
VariableNode::writeJsonAST(ostream &output) const
{
  output << R"({"node_type" : "VariableNode", )"
         << R"("name" : ")" << datatree.symbol_table.getName(symb_id) << R"(", "type" : ")";
  switch (get_type())
    {
    case SymbolType::endogenous:
      output << "endogenous";
      break;
    case SymbolType::exogenous:
      output << "exogenous";
      break;
    case SymbolType::exogenousDet:
      output << "exogenousDet";
      break;
    case SymbolType::parameter:
      output << "parameter";
      break;
    case SymbolType::modelLocalVariable:
      output << "modelLocalVariable";
      break;
    case SymbolType::modFileLocalVariable:
      output << "modFileLocalVariable";
      break;
    case SymbolType::externalFunction:
      output << "externalFunction";
      break;
    case SymbolType::trend:
      output << "trend";
      break;
    case SymbolType::statementDeclaredVariable:
      output << "statementDeclaredVariable";
      break;
    case SymbolType::logTrend:
      output << "logTrend:";
      break;
    case SymbolType::unusedEndogenous:
      output << "unusedEndogenous";
      break;
    case SymbolType::epilogue:
      output << "epilogue";
      break;
    case SymbolType::excludedVariable:
      cerr << "VariableNode::computeDerivative: Impossible case!" << endl;
      exit(EXIT_FAILURE);
    }
  output << R"(", "lag" : )" << lag << "}";
}

void
VariableNode::writeJsonOutput(ostream &output,
                              const temporary_terms_t &temporary_terms,
                              const deriv_node_temp_terms_t &tef_terms,
                              bool isdynamic) const
{
  if (temporary_terms.find(const_cast<VariableNode *>(this)) != temporary_terms.end())
    {
      output << "T" << idx;
      return;
    }

  output << datatree.symbol_table.getName(symb_id);
  if (isdynamic && lag != 0)
    output << "(" << lag << ")";
}

void
VariableNode::writeOutput(ostream &output, ExprNodeOutputType output_type,
                          const temporary_terms_t &temporary_terms,
                          const temporary_terms_idxs_t &temporary_terms_idxs,
                          const deriv_node_temp_terms_t &tef_terms) const
{
  auto type = get_type();
  if (checkIfTemporaryTermThenWrite(output, output_type, temporary_terms, temporary_terms_idxs))
    return;

  if (isLatexOutput(output_type))
    {
      if (output_type == ExprNodeOutputType::latexDynamicSteadyStateOperator)
        output << R"(\bar)";
      output << "{" << datatree.symbol_table.getTeXName(symb_id) << "}";
      if (output_type == ExprNodeOutputType::latexDynamicModel
          && (type == SymbolType::endogenous || type == SymbolType::exogenous || type == SymbolType::exogenousDet || type == SymbolType::trend || type == SymbolType::logTrend))
        {
          output << "_{t";
          if (lag != 0)
            {
              if (lag > 0)
                output << "+";
              output << lag;
            }
          output << "}";
        }
      return;
    }

  int i;
  switch (type)
    {
    case SymbolType::parameter:
      if (int tsid = datatree.symbol_table.getTypeSpecificID(symb_id);
          output_type == ExprNodeOutputType::matlabOutsideModel)
        output << "M_.params" << "(" << tsid + 1 << ")";
      else
        output << "params" << LEFT_ARRAY_SUBSCRIPT(output_type) << tsid + ARRAY_SUBSCRIPT_OFFSET(output_type) << RIGHT_ARRAY_SUBSCRIPT(output_type);
      break;

    case SymbolType::modelLocalVariable:
      if (output_type == ExprNodeOutputType::matlabDynamicSteadyStateOperator
          || output_type == ExprNodeOutputType::CDynamicSteadyStateOperator)
        {
          output << "(";
          datatree.getLocalVariable(symb_id)->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
          output << ")";
        }
      else
        /* We append underscores to avoid name clashes with "g1" or "oo_".
           But we probably never arrive here because MLV are temporary terms… */
        output << datatree.symbol_table.getName(symb_id) << "__";
      break;

    case SymbolType::modFileLocalVariable:
      output << datatree.symbol_table.getName(symb_id);
      break;

    case SymbolType::endogenous:
      switch (int tsid = datatree.symbol_table.getTypeSpecificID(symb_id);
              output_type)
        {
        case ExprNodeOutputType::juliaDynamicModel:
        case ExprNodeOutputType::matlabDynamicModel:
        case ExprNodeOutputType::CDynamicModel:
          i = datatree.getDynJacobianCol(datatree.getDerivID(symb_id, lag)) + ARRAY_SUBSCRIPT_OFFSET(output_type);
          output <<  "y" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << RIGHT_ARRAY_SUBSCRIPT(output_type);
          break;
        case ExprNodeOutputType::CStaticModel:
        case ExprNodeOutputType::juliaStaticModel:
        case ExprNodeOutputType::matlabStaticModel:
          i = tsid + ARRAY_SUBSCRIPT_OFFSET(output_type);
          output <<  "y" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << RIGHT_ARRAY_SUBSCRIPT(output_type);
          break;
        case ExprNodeOutputType::matlabOutsideModel:
          output << "oo_.steady_state(" << tsid + 1 << ")";
          break;
        case ExprNodeOutputType::juliaDynamicSteadyStateOperator:
        case ExprNodeOutputType::matlabDynamicSteadyStateOperator:
          output << "steady_state" << LEFT_ARRAY_SUBSCRIPT(output_type) << tsid + 1 << RIGHT_ARRAY_SUBSCRIPT(output_type);
          break;
        case ExprNodeOutputType::CDynamicSteadyStateOperator:
          output << "steady_state[" << tsid << "]";
          break;
        case ExprNodeOutputType::juliaSteadyStateFile:
        case ExprNodeOutputType::steadyStateFile:
          output << "ys_" << LEFT_ARRAY_SUBSCRIPT(output_type) << tsid + 1 << RIGHT_ARRAY_SUBSCRIPT(output_type);
          break;
        case ExprNodeOutputType::matlabDseries:
          output << "ds." << datatree.symbol_table.getName(symb_id);
          if (lag != 0)
            output << LEFT_ARRAY_SUBSCRIPT(output_type) << lag << RIGHT_ARRAY_SUBSCRIPT(output_type);
          break;
        case ExprNodeOutputType::epilogueFile:
          output << "ds." << datatree.symbol_table.getName(symb_id);
          output << LEFT_ARRAY_SUBSCRIPT(output_type) << "t";
          if (lag != 0)
            output << lag;
          output << RIGHT_ARRAY_SUBSCRIPT(output_type);
          break;
        case ExprNodeOutputType::occbinDifferenceFile:
          output << "zdatalinear(:," << tsid + 1 << ")";
          break;
        default:
          cerr << "VariableNode::writeOutput: should not reach this point" << endl;
          exit(EXIT_FAILURE);
        }
      break;

    case SymbolType::exogenous:
      i = datatree.symbol_table.getTypeSpecificID(symb_id) + ARRAY_SUBSCRIPT_OFFSET(output_type);
      switch (output_type)
        {
        case ExprNodeOutputType::juliaDynamicModel:
        case ExprNodeOutputType::matlabDynamicModel:
          if (lag > 0)
            output <<  "x" << LEFT_ARRAY_SUBSCRIPT(output_type) << "it_+" << lag << ", " << i
                   << RIGHT_ARRAY_SUBSCRIPT(output_type);
          else if (lag < 0)
            output <<  "x" << LEFT_ARRAY_SUBSCRIPT(output_type) << "it_" << lag << ", " << i
                   << RIGHT_ARRAY_SUBSCRIPT(output_type);
          else
            output <<  "x" << LEFT_ARRAY_SUBSCRIPT(output_type) << "it_, " << i
                   << RIGHT_ARRAY_SUBSCRIPT(output_type);
          break;
        case ExprNodeOutputType::CDynamicModel:
          if (lag == 0)
            output <<  "x[it_+" << i << "*nb_row_x]";
          else if (lag > 0)
            output <<  "x[it_+" << lag << "+" << i << "*nb_row_x]";
          else
            output <<  "x[it_" << lag << "+" << i << "*nb_row_x]";
          break;
        case ExprNodeOutputType::CStaticModel:
        case ExprNodeOutputType::juliaStaticModel:
        case ExprNodeOutputType::matlabStaticModel:
          output << "x" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << RIGHT_ARRAY_SUBSCRIPT(output_type);
          break;
        case ExprNodeOutputType::matlabOutsideModel:
          assert(lag == 0);
          output <<  "oo_.exo_steady_state(" << i << ")";
          break;
        case ExprNodeOutputType::matlabDynamicSteadyStateOperator:
          output <<  "oo_.exo_steady_state(" << i << ")";
          break;
        case ExprNodeOutputType::juliaSteadyStateFile:
        case ExprNodeOutputType::steadyStateFile:
          output << "exo_" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << RIGHT_ARRAY_SUBSCRIPT(output_type);
          break;
        case ExprNodeOutputType::matlabDseries:
          output << "ds." << datatree.symbol_table.getName(symb_id);
          if (lag != 0)
            output << LEFT_ARRAY_SUBSCRIPT(output_type) << lag << RIGHT_ARRAY_SUBSCRIPT(output_type);
          break;
        case ExprNodeOutputType::epilogueFile:
          output << "ds." << datatree.symbol_table.getName(symb_id);
          output << LEFT_ARRAY_SUBSCRIPT(output_type) << "t";
          if (lag != 0)
            output << lag;
          output << RIGHT_ARRAY_SUBSCRIPT(output_type);
          break;
        default:
          cerr << "VariableNode::writeOutput: should not reach this point" << endl;
          exit(EXIT_FAILURE);
        }
      break;

    case SymbolType::exogenousDet:
      i = datatree.symbol_table.getTypeSpecificID(symb_id) + datatree.symbol_table.exo_nbr() + ARRAY_SUBSCRIPT_OFFSET(output_type);
      switch (output_type)
        {
        case ExprNodeOutputType::juliaDynamicModel:
        case ExprNodeOutputType::matlabDynamicModel:
          if (lag > 0)
            output <<  "x" << LEFT_ARRAY_SUBSCRIPT(output_type) << "it_+" << lag << ", " << i
                   << RIGHT_ARRAY_SUBSCRIPT(output_type);
          else if (lag < 0)
            output <<  "x" << LEFT_ARRAY_SUBSCRIPT(output_type) << "it_" << lag << ", " << i
                   << RIGHT_ARRAY_SUBSCRIPT(output_type);
          else
            output <<  "x" << LEFT_ARRAY_SUBSCRIPT(output_type) << "it_, " << i
                   << RIGHT_ARRAY_SUBSCRIPT(output_type);
          break;
        case ExprNodeOutputType::CDynamicModel:
          if (lag == 0)
            output <<  "x[it_+" << i << "*nb_row_x]";
          else if (lag > 0)
            output <<  "x[it_+" << lag << "+" << i << "*nb_row_x]";
          else
            output <<  "x[it_" << lag << "+" << i << "*nb_row_x]";
          break;
        case ExprNodeOutputType::CStaticModel:
        case ExprNodeOutputType::juliaStaticModel:
        case ExprNodeOutputType::matlabStaticModel:
          output << "x" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << RIGHT_ARRAY_SUBSCRIPT(output_type);
          break;
        case ExprNodeOutputType::matlabOutsideModel:
          assert(lag == 0);
          output <<  "oo_.exo_det_steady_state(" << datatree.symbol_table.getTypeSpecificID(symb_id) + 1 << ")";
          break;
        case ExprNodeOutputType::matlabDynamicSteadyStateOperator:
          output <<  "oo_.exo_det_steady_state(" << datatree.symbol_table.getTypeSpecificID(symb_id) + 1 << ")";
          break;
        case ExprNodeOutputType::juliaSteadyStateFile:
        case ExprNodeOutputType::steadyStateFile:
          output << "exo_" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << RIGHT_ARRAY_SUBSCRIPT(output_type);
          break;
        case ExprNodeOutputType::matlabDseries:
          output << "ds." << datatree.symbol_table.getName(symb_id);
          if (lag != 0)
            output << LEFT_ARRAY_SUBSCRIPT(output_type) << lag << RIGHT_ARRAY_SUBSCRIPT(output_type);
          break;
        case ExprNodeOutputType::epilogueFile:
          output << "ds." << datatree.symbol_table.getName(symb_id);
          output << LEFT_ARRAY_SUBSCRIPT(output_type) << "t";
          if (lag != 0)
            output << lag;
          output << RIGHT_ARRAY_SUBSCRIPT(output_type);
          break;
        default:
          cerr << "VariableNode::writeOutput: should not reach this point" << endl;
          exit(EXIT_FAILURE);
        }
      break;
    case SymbolType::epilogue:
      if (output_type == ExprNodeOutputType::epilogueFile)
        {
          output << "ds." << datatree.symbol_table.getName(symb_id);
          output << LEFT_ARRAY_SUBSCRIPT(output_type) << "t";
          if (lag != 0)
            output << lag;
          output << RIGHT_ARRAY_SUBSCRIPT(output_type);
        }
      else if (output_type == ExprNodeOutputType::matlabDseries)
        // Only writing dseries for epilogue_static, hence no need to check lag
        output << "ds." << datatree.symbol_table.getName(symb_id);
      else
        {
          cerr << "VariableNode::writeOutput: Impossible case" << endl;
          exit(EXIT_FAILURE);
        }
      break;
    case SymbolType::unusedEndogenous:
      cerr << "ERROR: You cannot use an endogenous variable in an expression if that variable has not been used in the model block." << endl;
      exit(EXIT_FAILURE);
    case SymbolType::externalFunction:
    case SymbolType::trend:
    case SymbolType::logTrend:
    case SymbolType::statementDeclaredVariable:
    case SymbolType::excludedVariable:
      cerr << "VariableNode::writeOutput: Impossible case" << endl;
      exit(EXIT_FAILURE);
    }
}

expr_t
VariableNode::substituteStaticAuxiliaryVariable() const
{
  if (get_type() == SymbolType::endogenous)
    try
      {
        return datatree.symbol_table.getAuxiliaryVarsExprNode(symb_id)->substituteStaticAuxiliaryVariable();
      }
    catch (SymbolTable::SearchFailedException &e)
      {
      }
  return const_cast<VariableNode *>(this);
}

double
VariableNode::eval(const eval_context_t &eval_context) const noexcept(false)
{
  if (get_type() == SymbolType::modelLocalVariable)
    return datatree.getLocalVariable(symb_id)->eval(eval_context);

  auto it = eval_context.find(symb_id);
  if (it == eval_context.end())
    throw EvalException();

  return it->second;
}

void
VariableNode::compile(ostream &CompileCode, unsigned int &instruction_number,
                      bool lhs_rhs, const temporary_terms_t &temporary_terms,
                      const temporary_terms_idxs_t &temporary_terms_idxs, bool dynamic, bool steady_dynamic,
                      const deriv_node_temp_terms_t &tef_terms) const
{
  auto type = get_type();
  if (type == SymbolType::modelLocalVariable || type == SymbolType::modFileLocalVariable)
    datatree.getLocalVariable(symb_id)->compile(CompileCode, instruction_number, lhs_rhs, temporary_terms, temporary_terms_idxs, dynamic, steady_dynamic, tef_terms);
  else
    {
      int tsid = datatree.symbol_table.getTypeSpecificID(symb_id);
      if (type == SymbolType::exogenousDet)
        tsid += datatree.symbol_table.exo_nbr();
      if (!lhs_rhs)
        {
          if (dynamic)
            {
              if (steady_dynamic) // steady state values in a dynamic model
                {
                  FLDVS_ fldvs{static_cast<uint8_t>(type), static_cast<unsigned int>(tsid)};
                  fldvs.write(CompileCode, instruction_number);
                }
              else
                {
                  if (type == SymbolType::parameter)
                    {
                      FLDV_ fldv{static_cast<int>(type), static_cast<unsigned int>(tsid)};
                      fldv.write(CompileCode, instruction_number);
                    }
                  else
                    {
                      FLDV_ fldv{static_cast<int>(type), static_cast<unsigned int>(tsid), lag};
                      fldv.write(CompileCode, instruction_number);
                    }
                }
            }
          else
            {
              FLDSV_ fldsv{static_cast<uint8_t>(type), static_cast<unsigned int>(tsid)};
              fldsv.write(CompileCode, instruction_number);
            }
        }
      else
        {
          if (dynamic)
            {
              if (steady_dynamic) // steady state values in a dynamic model
                {
                  cerr << "Impossible case: steady_state in rhs of equation" << endl;
                  exit(EXIT_FAILURE);
                }
              else
                {
                  if (type == SymbolType::parameter)
                    {
                      FSTPV_ fstpv{static_cast<int>(type), static_cast<unsigned int>(tsid)};
                      fstpv.write(CompileCode, instruction_number);
                    }
                  else
                    {
                      FSTPV_ fstpv{static_cast<int>(type), static_cast<unsigned int>(tsid), lag};
                      fstpv.write(CompileCode, instruction_number);
                    }
                }
            }
          else
            {
              FSTPSV_ fstpsv{static_cast<uint8_t>(type), static_cast<unsigned int>(tsid)};
              fstpsv.write(CompileCode, instruction_number);
            }
        }
    }
}

void
VariableNode::collectVARLHSVariable(set<expr_t> &result) const
{
  if (get_type() == SymbolType::endogenous && lag == 0)
    result.insert(const_cast<VariableNode *>(this));
  else
    {
      cerr << "ERROR: you can only have endogenous variables or unary ops on LHS of VAR" << endl;
      exit(EXIT_FAILURE);
    }
}

void
VariableNode::collectDynamicVariables(SymbolType type_arg, set<pair<int, int>> &result) const
{
  if (get_type() == type_arg)
    result.emplace(symb_id, lag);
  if (get_type() == SymbolType::modelLocalVariable)
    datatree.getLocalVariable(symb_id)->collectDynamicVariables(type_arg, result);
}

void
VariableNode::computeSubExprContainingVariable(int symb_id_arg, int lag_arg, set<expr_t> &contain_var) const
{
  if (symb_id == symb_id_arg && lag == lag_arg)
    contain_var.insert(const_cast<VariableNode*>(this));
  if (get_type() == SymbolType::modelLocalVariable)
    datatree.getLocalVariable(symb_id)->computeSubExprContainingVariable(symb_id_arg, lag_arg, contain_var);
}

BinaryOpNode *
VariableNode::normalizeEquationHelper(const set<expr_t> &contain_var, expr_t rhs) const
{
  assert(contain_var.count(const_cast<VariableNode *>(this)) > 0);

  if (get_type() == SymbolType::modelLocalVariable)
    return datatree.getLocalVariable(symb_id)->normalizeEquationHelper(contain_var, rhs);

  // This the LHS variable: we have finished the normalization
  return datatree.AddEqual(const_cast<VariableNode *>(this), rhs);
}

expr_t
VariableNode::getChainRuleDerivative(int deriv_id, const map<int, BinaryOpNode *> &recursive_variables)
{
  switch (get_type())
    {
    case SymbolType::endogenous:
    case SymbolType::exogenous:
    case SymbolType::exogenousDet:
    case SymbolType::parameter:
    case SymbolType::trend:
    case SymbolType::logTrend:
      if (deriv_id == datatree.getDerivID(symb_id, lag))
        return datatree.One;
      // If there is in the equation a recursive variable we could use a chaine rule derivation
      else if (auto it = recursive_variables.find(datatree.getDerivID(symb_id, lag));
               it != recursive_variables.end())
        return it->second->arg2->getChainRuleDerivative(deriv_id, recursive_variables);
      else
        return datatree.Zero;

    case SymbolType::modelLocalVariable:
      return datatree.getLocalVariable(symb_id)->getChainRuleDerivative(deriv_id, recursive_variables);
    case SymbolType::modFileLocalVariable:
      cerr << "modFileLocalVariable is not derivable" << endl;
      exit(EXIT_FAILURE);
    case SymbolType::statementDeclaredVariable:
      cerr << "statementDeclaredVariable is not derivable" << endl;
      exit(EXIT_FAILURE);
    case SymbolType::unusedEndogenous:
      cerr << "unusedEndogenous is not derivable" << endl;
      exit(EXIT_FAILURE);
    case SymbolType::externalFunction:
    case SymbolType::epilogue:
    case SymbolType::excludedVariable:
      cerr << "VariableNode::getChainRuleDerivative: Impossible case" << endl;
      exit(EXIT_FAILURE);
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

expr_t
VariableNode::toStatic(DataTree &static_datatree) const
{
  return static_datatree.AddVariable(symb_id);
}

void
VariableNode::computeXrefs(EquationInfo &ei) const
{
  switch (get_type())
    {
    case SymbolType::endogenous:
      ei.endo.emplace(symb_id, lag);
      break;
    case SymbolType::exogenous:
      ei.exo.emplace(symb_id, lag);
      break;
    case SymbolType::exogenousDet:
      ei.exo_det.emplace(symb_id, lag);
      break;
    case SymbolType::parameter:
      ei.param.emplace(symb_id, 0);
      break;
    case SymbolType::modFileLocalVariable:
      datatree.getLocalVariable(symb_id)->computeXrefs(ei);
      break;
    case SymbolType::trend:
    case SymbolType::logTrend:
    case SymbolType::modelLocalVariable:
    case SymbolType::statementDeclaredVariable:
    case SymbolType::unusedEndogenous:
    case SymbolType::externalFunction:
    case SymbolType::epilogue:
    case SymbolType::excludedVariable:
      break;
    }
}

SymbolType
VariableNode::get_type() const
{
  return datatree.symbol_table.getType(symb_id);
}

expr_t
VariableNode::clone(DataTree &datatree) const
{
  return datatree.AddVariable(symb_id, lag);
}

int
VariableNode::maxEndoLead() const
{
  switch (get_type())
    {
    case SymbolType::endogenous:
      return max(lag, 0);
    case SymbolType::modelLocalVariable:
      return datatree.getLocalVariable(symb_id)->maxEndoLead();
    default:
      return 0;
    }
}

int
VariableNode::maxExoLead() const
{
  switch (get_type())
    {
    case SymbolType::exogenous:
      return max(lag, 0);
    case SymbolType::modelLocalVariable:
      return datatree.getLocalVariable(symb_id)->maxExoLead();
    default:
      return 0;
    }
}

int
VariableNode::maxEndoLag() const
{
  switch (get_type())
    {
    case SymbolType::endogenous:
      return max(-lag, 0);
    case SymbolType::modelLocalVariable:
      return datatree.getLocalVariable(symb_id)->maxEndoLag();
    default:
      return 0;
    }
}

int
VariableNode::maxExoLag() const
{
  switch (get_type())
    {
    case SymbolType::exogenous:
      return max(-lag, 0);
    case SymbolType::modelLocalVariable:
      return datatree.getLocalVariable(symb_id)->maxExoLag();
    default:
      return 0;
    }
}

int
VariableNode::maxLead() const
{
  switch (get_type())
    {
    case SymbolType::endogenous:
    case SymbolType::exogenous:
    case SymbolType::exogenousDet:
      return lag;
    case SymbolType::modelLocalVariable:
      return datatree.getLocalVariable(symb_id)->maxLead();
    default:
      return 0;
    }
}

int
VariableNode::maxLag() const
{
  switch (get_type())
    {
    case SymbolType::endogenous:
    case SymbolType::exogenous:
    case SymbolType::exogenousDet:
      return -lag;
    case SymbolType::modelLocalVariable:
      return datatree.getLocalVariable(symb_id)->maxLag();
    default:
      return 0;
    }
}

int
VariableNode::maxLagWithDiffsExpanded() const
{
  switch (get_type())
    {
    case SymbolType::endogenous:
    case SymbolType::exogenous:
    case SymbolType::exogenousDet:
    case SymbolType::epilogue:
      return -lag;
    case SymbolType::modelLocalVariable:
      return datatree.getLocalVariable(symb_id)->maxLagWithDiffsExpanded();
    default:
      return 0;
    }
}

expr_t
VariableNode::undiff() const
{
  return const_cast<VariableNode *>(this);
}

int
VariableNode::VarMaxLag(const set<expr_t> &lhs_lag_equiv) const
{
  auto [lag_equiv_repr, index] = getLagEquivalenceClass();
  if (lhs_lag_equiv.find(lag_equiv_repr) == lhs_lag_equiv.end())
    return 0;
  return maxLag();
}

expr_t
VariableNode::substituteAdl() const
{
  /* Do not recurse into model-local variables definition, rather do it at the
     DynamicModel method level (see the comment there) */
  return const_cast<VariableNode *>(this);
}

expr_t
VariableNode::substituteModelLocalVariables() const
{
  if (get_type() == SymbolType::modelLocalVariable)
    return datatree.getLocalVariable(symb_id);

  return const_cast<VariableNode *>(this);
}

expr_t
VariableNode::substituteVarExpectation(const map<string, expr_t> &subst_table) const
{
  if (get_type() == SymbolType::modelLocalVariable)
    return datatree.getLocalVariable(symb_id)->substituteVarExpectation(subst_table);

  return const_cast<VariableNode *>(this);
}

void
VariableNode::findDiffNodes(lag_equivalence_table_t &nodes) const
{
  if (get_type() == SymbolType::modelLocalVariable)
    datatree.getLocalVariable(symb_id)->findDiffNodes(nodes);
}

void
VariableNode::findUnaryOpNodesForAuxVarCreation(lag_equivalence_table_t &nodes) const
{
  if (get_type() == SymbolType::modelLocalVariable)
    datatree.getLocalVariable(symb_id)->findUnaryOpNodesForAuxVarCreation(nodes);
}

int
VariableNode::findTargetVariable(int lhs_symb_id) const
{
  if (get_type() == SymbolType::modelLocalVariable)
    return datatree.getLocalVariable(symb_id)->findTargetVariable(lhs_symb_id);

  return -1;
}

expr_t
VariableNode::substituteDiff(const lag_equivalence_table_t &nodes, subst_table_t &subst_table,
                             vector<BinaryOpNode *> &neweqs) const
{
  if (get_type() == SymbolType::modelLocalVariable)
    return datatree.getLocalVariable(symb_id)->substituteDiff(nodes, subst_table, neweqs);

  return const_cast<VariableNode *>(this);
}

expr_t
VariableNode::substituteUnaryOpNodes(const lag_equivalence_table_t &nodes, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  if (get_type() == SymbolType::modelLocalVariable)
    return datatree.getLocalVariable(symb_id)->substituteUnaryOpNodes(nodes, subst_table, neweqs);

  return const_cast<VariableNode *>(this);
}

expr_t
VariableNode::substitutePacExpectation(const string &name, expr_t subexpr)
{
  if (get_type() == SymbolType::modelLocalVariable)
    return datatree.getLocalVariable(symb_id)->substitutePacExpectation(name, subexpr);

  return const_cast<VariableNode *>(this);
}

expr_t
VariableNode::substitutePacTargetNonstationary(const string &name, expr_t subexpr)
{
  if (get_type() == SymbolType::modelLocalVariable)
    return datatree.getLocalVariable(symb_id)->substitutePacTargetNonstationary(name, subexpr);

  return const_cast<VariableNode *>(this);
}

expr_t
VariableNode::decreaseLeadsLags(int n) const
{
  switch (get_type())
    {
    case SymbolType::endogenous:
    case SymbolType::exogenous:
    case SymbolType::exogenousDet:
    case SymbolType::trend:
    case SymbolType::logTrend:
      return datatree.AddVariable(symb_id, lag-n);
    case SymbolType::modelLocalVariable:
      return datatree.getLocalVariable(symb_id)->decreaseLeadsLags(n);
    default:
      return const_cast<VariableNode *>(this);
    }
}

expr_t
VariableNode::decreaseLeadsLagsPredeterminedVariables() const
{
  /* Do not recurse into model-local variables definitions, since MLVs are
     already handled by DynamicModel::transformPredeterminedVariables().
     This is also necessary because of #65. */
  if (datatree.symbol_table.isPredetermined(symb_id))
    return decreaseLeadsLags(1);
  else
    return const_cast<VariableNode *>(this);
}

expr_t
VariableNode::substituteEndoLeadGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const
{
  switch (get_type())
    {
    case SymbolType::endogenous:
      if (lag <= 1)
        return const_cast<VariableNode *>(this);
      else
        return createEndoLeadAuxiliaryVarForMyself(subst_table, neweqs);
    case SymbolType::modelLocalVariable:
      if (expr_t value = datatree.getLocalVariable(symb_id); value->maxEndoLead() <= 1)
        return const_cast<VariableNode *>(this);
      else
        return value->substituteEndoLeadGreaterThanTwo(subst_table, neweqs, deterministic_model);
    default:
      return const_cast<VariableNode *>(this);
    }
}

expr_t
VariableNode::substituteEndoLagGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  VariableNode *substexpr;
  int cur_lag;
  switch (get_type())
    {
    case SymbolType::endogenous:
      if (lag >= -1)
        return const_cast<VariableNode *>(this);

      if (auto it = subst_table.find(this); it != subst_table.end())
        return const_cast<VariableNode *>(it->second);

      substexpr = datatree.AddVariable(symb_id, -1);
      cur_lag = -2;

      // Each iteration tries to create an auxvar such that auxvar(-1)=curvar(cur_lag)
      // At the beginning (resp. end) of each iteration, substexpr is an expression (possibly an auxvar) equivalent to curvar(cur_lag+1) (resp. curvar(cur_lag))
      while (cur_lag >= lag)
        {
          VariableNode *orig_expr = datatree.AddVariable(symb_id, cur_lag);
          if (auto it = subst_table.find(orig_expr); it == subst_table.end())
            {
              int aux_symb_id = datatree.symbol_table.addEndoLagAuxiliaryVar(symb_id, cur_lag+1, substexpr);
              neweqs.push_back(datatree.AddEqual(datatree.AddVariable(aux_symb_id, 0), substexpr));
              substexpr = datatree.AddVariable(aux_symb_id, -1);
              subst_table[orig_expr] = substexpr;
            }
          else
            substexpr = const_cast<VariableNode *>(it->second);

          cur_lag--;
        }
      return substexpr;

    case SymbolType::modelLocalVariable:
      if (expr_t value = datatree.getLocalVariable(symb_id); value->maxEndoLag() <= 1)
        return const_cast<VariableNode *>(this);
      else
        return value->substituteEndoLagGreaterThanTwo(subst_table, neweqs);
    default:
      return const_cast<VariableNode *>(this);
    }
}

expr_t
VariableNode::substituteExoLead(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const
{
  switch (get_type())
    {
    case SymbolType::exogenous:
      if (lag <= 0)
        return const_cast<VariableNode *>(this);
      else
        return createExoLeadAuxiliaryVarForMyself(subst_table, neweqs);
    case SymbolType::modelLocalVariable:
      if (expr_t value = datatree.getLocalVariable(symb_id); value->maxExoLead() == 0)
        return const_cast<VariableNode *>(this);
      else
        return value->substituteExoLead(subst_table, neweqs, deterministic_model);
    default:
      return const_cast<VariableNode *>(this);
    }
}

expr_t
VariableNode::substituteExoLag(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  VariableNode *substexpr;
  int cur_lag;
  switch (get_type())
    {
    case SymbolType::exogenous:
      if (lag >= 0)
        return const_cast<VariableNode *>(this);

      if (auto it = subst_table.find(this); it != subst_table.end())
        return const_cast<VariableNode *>(it->second);

      substexpr = datatree.AddVariable(symb_id, 0);
      cur_lag = -1;

      // Each iteration tries to create an auxvar such that auxvar(-1)=curvar(cur_lag)
      // At the beginning (resp. end) of each iteration, substexpr is an expression (possibly an auxvar) equivalent to curvar(cur_lag+1) (resp. curvar(cur_lag))
      while (cur_lag >= lag)
        {
          VariableNode *orig_expr = datatree.AddVariable(symb_id, cur_lag);
          if (auto it = subst_table.find(orig_expr); it == subst_table.end())
            {
              int aux_symb_id = datatree.symbol_table.addExoLagAuxiliaryVar(symb_id, cur_lag+1, substexpr);
              neweqs.push_back(datatree.AddEqual(datatree.AddVariable(aux_symb_id, 0), substexpr));
              substexpr = datatree.AddVariable(aux_symb_id, -1);
              subst_table[orig_expr] = substexpr;
            }
          else
            substexpr = const_cast<VariableNode *>(it->second);

          cur_lag--;
        }
      return substexpr;

    case SymbolType::modelLocalVariable:
      if (expr_t value = datatree.getLocalVariable(symb_id); value->maxExoLag() == 0)
        return const_cast<VariableNode *>(this);
      else
        return value->substituteExoLag(subst_table, neweqs);
    default:
      return const_cast<VariableNode *>(this);
    }
}

expr_t
VariableNode::substituteExpectation(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool partial_information_model) const
{
  if (get_type() == SymbolType::modelLocalVariable)
    return datatree.getLocalVariable(symb_id)->substituteExpectation(subst_table, neweqs, partial_information_model);

  return const_cast<VariableNode *>(this);
}

expr_t
VariableNode::differentiateForwardVars(const vector<string> &subset, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  switch (get_type())
    {
    case SymbolType::endogenous:
      assert(lag <= 1);
      if (lag <= 0
          || (subset.size() > 0
              && find(subset.begin(), subset.end(), datatree.symbol_table.getName(symb_id)) == subset.end()))
        return const_cast<VariableNode *>(this);
      else
        {
          VariableNode *diffvar;
          if (auto it = subst_table.find(this); it != subst_table.end())
            diffvar = const_cast<VariableNode *>(it->second);
          else
            {
              int aux_symb_id = datatree.symbol_table.addDiffForwardAuxiliaryVar(symb_id, datatree.AddMinus(datatree.AddVariable(symb_id, 0),
                                                                                                            datatree.AddVariable(symb_id, -1)));
              neweqs.push_back(datatree.AddEqual(datatree.AddVariable(aux_symb_id, 0), datatree.AddMinus(datatree.AddVariable(symb_id, 0),
                                                                                                                                      datatree.AddVariable(symb_id, -1))));
              diffvar = datatree.AddVariable(aux_symb_id, 1);
              subst_table[this] = diffvar;
            }
          return datatree.AddPlus(datatree.AddVariable(symb_id, 0), diffvar);
        }
    case SymbolType::modelLocalVariable:
      if (expr_t value = datatree.getLocalVariable(symb_id); value->maxEndoLead() <= 0)
        return const_cast<VariableNode *>(this);
      else
        return value->differentiateForwardVars(subset, subst_table, neweqs);
    default:
      return const_cast<VariableNode *>(this);
    }
}

bool
VariableNode::isNumConstNodeEqualTo(double value) const
{
  return false;
}

bool
VariableNode::isVariableNodeEqualTo(SymbolType type_arg, int variable_id, int lag_arg) const
{
  if (get_type() == type_arg && datatree.symbol_table.getTypeSpecificID(symb_id) == variable_id && lag == lag_arg)
    return true;
  else
    return false;
}

bool
VariableNode::containsPacExpectation(const string &pac_model_name) const
{
  if (get_type() == SymbolType::modelLocalVariable)
    return datatree.getLocalVariable(symb_id)->containsPacExpectation(pac_model_name);

  return false;
}

bool
VariableNode::containsPacTargetNonstationary(const string &pac_model_name) const
{
  if (get_type() == SymbolType::modelLocalVariable)
    return datatree.getLocalVariable(symb_id)->containsPacTargetNonstationary(pac_model_name);

  return false;
}

expr_t
VariableNode::replaceTrendVar() const
{
  if (get_type() == SymbolType::modelLocalVariable)
    return datatree.getLocalVariable(symb_id)->replaceTrendVar();

  if (get_type() == SymbolType::trend)
    return datatree.One;
  else if (get_type() == SymbolType::logTrend)
    return datatree.Zero;
  else
    return const_cast<VariableNode *>(this);
}

expr_t
VariableNode::detrend(int symb_id, bool log_trend, expr_t trend) const
{
  if (get_type() == SymbolType::modelLocalVariable)
    return datatree.getLocalVariable(symb_id)->detrend(symb_id, log_trend, trend);

  if (this->symb_id != symb_id)
    return const_cast<VariableNode *>(this);

  if (log_trend)
    {
      if (lag == 0)
        return datatree.AddPlus(const_cast<VariableNode *>(this), trend);
      else
        return datatree.AddPlus(const_cast<VariableNode *>(this), trend->decreaseLeadsLags(-lag));
    }
  else
    {
      if (lag == 0)
        return datatree.AddTimes(const_cast<VariableNode *>(this), trend);
      else
        return datatree.AddTimes(const_cast<VariableNode *>(this), trend->decreaseLeadsLags(-lag));
    }
}

int
VariableNode::countDiffs() const
{
  if (get_type() == SymbolType::modelLocalVariable)
    return datatree.getLocalVariable(symb_id)->countDiffs();

  return 0;
}

expr_t
VariableNode::removeTrendLeadLag(const map<int, expr_t> &trend_symbols_map) const
{
  if (get_type() == SymbolType::modelLocalVariable)
    return datatree.getLocalVariable(symb_id)->removeTrendLeadLag(trend_symbols_map);

  if ((get_type() != SymbolType::trend && get_type() != SymbolType::logTrend) || lag == 0)
    return const_cast<VariableNode *>(this);

  auto it = trend_symbols_map.find(symb_id);
  expr_t noTrendLeadLagNode = datatree.AddVariable(it->first);
  bool log_trend = get_type() == SymbolType::logTrend;
  expr_t trend = it->second;

  if (lag > 0)
    {
      expr_t growthFactorSequence = trend->decreaseLeadsLags(-1);
      if (log_trend)
        {
          for (int i = 1; i < lag; i++)
            growthFactorSequence = datatree.AddPlus(growthFactorSequence, trend->decreaseLeadsLags(-1*(i+1)));
          return datatree.AddPlus(noTrendLeadLagNode, growthFactorSequence);
        }
      else
        {
          for (int i = 1; i < lag; i++)
            growthFactorSequence = datatree.AddTimes(growthFactorSequence, trend->decreaseLeadsLags(-1*(i+1)));
          return datatree.AddTimes(noTrendLeadLagNode, growthFactorSequence);
        }
    }
  else //get_lag < 0
    {
      expr_t growthFactorSequence = trend;
      if (log_trend)
        {
          for (int i = 1; i < abs(lag); i++)
            growthFactorSequence = datatree.AddPlus(growthFactorSequence, trend->decreaseLeadsLags(i));
          return datatree.AddMinus(noTrendLeadLagNode, growthFactorSequence);
        }
      else
        {
          for (int i = 1; i < abs(lag); i++)
            growthFactorSequence = datatree.AddTimes(growthFactorSequence, trend->decreaseLeadsLags(i));
          return datatree.AddDivide(noTrendLeadLagNode, growthFactorSequence);
        }
    }
}

bool
VariableNode::isInStaticForm() const
{
  if (get_type() == SymbolType::modelLocalVariable)
    return datatree.getLocalVariable(symb_id)->isInStaticForm();

  return lag == 0;
}

bool
VariableNode::isParamTimesEndogExpr() const
{
  if (get_type() == SymbolType::modelLocalVariable)
    return datatree.getLocalVariable(symb_id)->isParamTimesEndogExpr();

  return false;
}

expr_t
VariableNode::replaceVarsInEquation(map<VariableNode *, NumConstNode *> &table) const
{
  /* Do not recurse into model-local variables definitions, since MLVs are
     already handled by DynamicModel::simplifyEquations().
     This is also necessary because of #65. */
  for (auto &it : table)
    if (it.first->symb_id == symb_id)
      return it.second;
  return const_cast<VariableNode *>(this);
}

void
VariableNode::matchMatchedMoment(vector<int> &symb_ids, vector<int> &lags, vector<int> &powers) const
{
  /* Used for simple expression outside model block, so no need to special-case
     model local variables */

  if (get_type() != SymbolType::endogenous)
    throw MatchFailureException{"Variable " + datatree.symbol_table.getName(symb_id) + " is not an endogenous"};

  symb_ids.push_back(symb_id);
  lags.push_back(lag);
  powers.push_back(1);
}

UnaryOpNode::UnaryOpNode(DataTree &datatree_arg, int idx_arg, UnaryOpcode op_code_arg, const expr_t arg_arg, int expectation_information_set_arg, int param1_symb_id_arg, int param2_symb_id_arg, string adl_param_name_arg, vector<int> adl_lags_arg) :
  ExprNode{datatree_arg, idx_arg},
  arg{arg_arg},
  expectation_information_set{expectation_information_set_arg},
  param1_symb_id{param1_symb_id_arg},
  param2_symb_id{param2_symb_id_arg},
  op_code{op_code_arg},
  adl_param_name{move(adl_param_name_arg)},
  adl_lags{move(adl_lags_arg)}
{
}

void
UnaryOpNode::prepareForDerivation()
{
  if (preparedForDerivation)
    return;

  preparedForDerivation = true;

  arg->prepareForDerivation();

  // Non-null derivatives are those of the argument (except for STEADY_STATE)
  non_null_derivatives = arg->non_null_derivatives;
  if (op_code == UnaryOpcode::steadyState || op_code == UnaryOpcode::steadyStateParamDeriv
      || op_code == UnaryOpcode::steadyStateParam2ndDeriv)
    datatree.addAllParamDerivId(non_null_derivatives);
}

expr_t
UnaryOpNode::composeDerivatives(expr_t darg, int deriv_id)
{
  expr_t t11, t12, t13, t14, t15;

  switch (op_code)
    {
    case UnaryOpcode::uminus:
      return datatree.AddUMinus(darg);
    case UnaryOpcode::exp:
      return datatree.AddTimes(darg, this);
    case UnaryOpcode::log:
      return datatree.AddDivide(darg, arg);
    case UnaryOpcode::log10:
      t11 = datatree.AddExp(datatree.One);
      t12 = datatree.AddLog10(t11);
      t13 = datatree.AddDivide(darg, arg);
      return datatree.AddTimes(t12, t13);
    case UnaryOpcode::cos:
      t11 = datatree.AddSin(arg);
      t12 = datatree.AddUMinus(t11);
      return datatree.AddTimes(darg, t12);
    case UnaryOpcode::sin:
      t11 = datatree.AddCos(arg);
      return datatree.AddTimes(darg, t11);
    case UnaryOpcode::tan:
      t11 = datatree.AddTimes(this, this);
      t12 = datatree.AddPlus(t11, datatree.One);
      return datatree.AddTimes(darg, t12);
    case UnaryOpcode::acos:
      t11 = datatree.AddSin(this);
      t12 = datatree.AddDivide(darg, t11);
      return datatree.AddUMinus(t12);
    case UnaryOpcode::asin:
      t11 = datatree.AddCos(this);
      return datatree.AddDivide(darg, t11);
    case UnaryOpcode::atan:
      t11 = datatree.AddTimes(arg, arg);
      t12 = datatree.AddPlus(datatree.One, t11);
      return datatree.AddDivide(darg, t12);
    case UnaryOpcode::cosh:
      t11 = datatree.AddSinh(arg);
      return datatree.AddTimes(darg, t11);
    case UnaryOpcode::sinh:
      t11 = datatree.AddCosh(arg);
      return datatree.AddTimes(darg, t11);
    case UnaryOpcode::tanh:
      t11 = datatree.AddTimes(this, this);
      t12 = datatree.AddMinus(datatree.One, t11);
      return datatree.AddTimes(darg, t12);
    case UnaryOpcode::acosh:
      t11 = datatree.AddSinh(this);
      return datatree.AddDivide(darg, t11);
    case UnaryOpcode::asinh:
      t11 = datatree.AddCosh(this);
      return datatree.AddDivide(darg, t11);
    case UnaryOpcode::atanh:
      t11 = datatree.AddTimes(arg, arg);
      t12 = datatree.AddMinus(datatree.One, t11);
      return datatree.AddDivide(darg, t12);
    case UnaryOpcode::sqrt:
      t11 = datatree.AddPlus(this, this);
      return datatree.AddDivide(darg, t11);
    case UnaryOpcode::cbrt:
      t11 = datatree.AddPower(arg, datatree.AddDivide(datatree.Two, datatree.Three));
      t12 = datatree.AddTimes(datatree.Three, t11);
      return datatree.AddDivide(darg, t12);
    case UnaryOpcode::abs:
      t11 = datatree.AddSign(arg);
      return datatree.AddTimes(t11, darg);
    case UnaryOpcode::sign:
      return datatree.Zero;
    case UnaryOpcode::steadyState:
      if (datatree.isDynamic())
        {
          if (datatree.getTypeByDerivID(deriv_id) == SymbolType::parameter)
            {
              auto varg = dynamic_cast<VariableNode *>(arg);
              if (!varg)
                {
                  cerr << "UnaryOpNode::composeDerivatives: STEADY_STATE() should only be used on "
                       << "standalone variables (like STEADY_STATE(y)) to be derivable w.r.t. parameters" << endl;
                  exit(EXIT_FAILURE);
                }
              if (datatree.symbol_table.getType(varg->symb_id) == SymbolType::endogenous)
                return datatree.AddSteadyStateParamDeriv(arg, datatree.getSymbIDByDerivID(deriv_id));
              else
                return datatree.Zero;
            }
          else
            return datatree.Zero;
        }
      else
        return darg;
    case UnaryOpcode::steadyStateParamDeriv:
      assert(datatree.isDynamic());
      if (datatree.getTypeByDerivID(deriv_id) == SymbolType::parameter)
        {
          auto varg = dynamic_cast<VariableNode *>(arg);
          assert(varg);
          assert(datatree.symbol_table.getType(varg->symb_id) == SymbolType::endogenous);
          return datatree.AddSteadyStateParam2ndDeriv(arg, param1_symb_id, datatree.getSymbIDByDerivID(deriv_id));
        }
      else
        return datatree.Zero;
    case UnaryOpcode::steadyStateParam2ndDeriv:
      assert(datatree.isDynamic());
      if (datatree.getTypeByDerivID(deriv_id) == SymbolType::parameter)
        {
          cerr << "3rd derivative of STEADY_STATE node w.r.t. three parameters not implemented" << endl;
          exit(EXIT_FAILURE);
        }
      else
        return datatree.Zero;
    case UnaryOpcode::expectation:
      cerr << "UnaryOpNode::composeDerivatives: not implemented on UnaryOpcode::expectation" << endl;
      exit(EXIT_FAILURE);
    case UnaryOpcode::erf:
    case UnaryOpcode::erfc:
      // x^2
      t11 = datatree.AddPower(arg, datatree.Two);
      // exp(x^2)
      t12 = datatree.AddExp(t11);
      // sqrt(pi)
      t11 = datatree.AddSqrt(datatree.Pi);
      // sqrt(pi)*exp(x^2)
      t13 = datatree.AddTimes(t11, t12);
      // 2/(sqrt(pi)*exp(x^2));
      t14 = datatree.AddDivide(datatree.Two, t13);
      // (2/(sqrt(pi)*exp(x^2)))*dx;
      t15 = datatree.AddTimes(t14, darg);
      if (op_code == UnaryOpcode::erf)
        return t15;
      else // erfc
        return datatree.AddUMinus(t15);
    case UnaryOpcode::diff:
      cerr << "UnaryOpNode::composeDerivatives: not implemented on UnaryOpcode::diff" << endl;
      exit(EXIT_FAILURE);
    case UnaryOpcode::adl:
      cerr << "UnaryOpNode::composeDerivatives: not implemented on UnaryOpcode::adl" << endl;
      exit(EXIT_FAILURE);
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

expr_t
UnaryOpNode::computeDerivative(int deriv_id)
{
  expr_t darg = arg->getDerivative(deriv_id);
  return composeDerivatives(darg, deriv_id);
}

int
UnaryOpNode::cost(const map<pair<int, int>, temporary_terms_t> &temp_terms_map, bool is_matlab) const
{
  // For a temporary term, the cost is null
  for (const auto &it : temp_terms_map)
    if (it.second.find(const_cast<UnaryOpNode *>(this)) != it.second.end())
      return 0;

  return cost(arg->cost(temp_terms_map, is_matlab), is_matlab);
}

int
UnaryOpNode::cost(const vector<vector<temporary_terms_t>> &blocks_temporary_terms, bool is_matlab) const
{
  // For a temporary term, the cost is null
  for (const auto &blk_tt : blocks_temporary_terms)
    for (const auto &eq_tt : blk_tt)
      if (eq_tt.find(const_cast<UnaryOpNode *>(this)) != eq_tt.end())
        return 0;

  return cost(arg->cost(blocks_temporary_terms, is_matlab), is_matlab);
}

int
UnaryOpNode::cost(int cost, bool is_matlab) const
{
  if (op_code == UnaryOpcode::uminus && dynamic_cast<NumConstNode *>(arg))
    return 0; // Cost is zero for a negative constant, as for a positive one

  if (is_matlab)
    // Cost for Matlab files
    switch (op_code)
      {
      case UnaryOpcode::uminus:
      case UnaryOpcode::sign:
        return cost + 70;
      case UnaryOpcode::exp:
        return cost + 160;
      case UnaryOpcode::log:
        return cost + 300;
      case UnaryOpcode::log10:
      case UnaryOpcode::erf:
      case UnaryOpcode::erfc:
        return cost + 16000;
      case UnaryOpcode::cos:
      case UnaryOpcode::sin:
      case UnaryOpcode::cosh:
        return cost + 210;
      case UnaryOpcode::tan:
        return cost + 230;
      case UnaryOpcode::acos:
        return cost + 300;
      case UnaryOpcode::asin:
        return cost + 310;
      case UnaryOpcode::atan:
        return cost + 140;
      case UnaryOpcode::sinh:
        return cost + 240;
      case UnaryOpcode::tanh:
        return cost + 190;
      case UnaryOpcode::acosh:
        return cost + 770;
      case UnaryOpcode::asinh:
        return cost + 460;
      case UnaryOpcode::atanh:
        return cost + 350;
      case UnaryOpcode::sqrt:
      case UnaryOpcode::cbrt:
      case UnaryOpcode::abs:
        return cost + 570;
      case UnaryOpcode::steadyState:
      case UnaryOpcode::steadyStateParamDeriv:
      case UnaryOpcode::steadyStateParam2ndDeriv:
      case UnaryOpcode::expectation:
        return cost;
      case UnaryOpcode::diff:
        cerr << "UnaryOpNode::cost: not implemented on UnaryOpcode::diff" << endl;
        exit(EXIT_FAILURE);
      case UnaryOpcode::adl:
        cerr << "UnaryOpNode::cost: not implemented on UnaryOpcode::adl" << endl;
        exit(EXIT_FAILURE);
      }
  else
    // Cost for C files
    switch (op_code)
      {
      case UnaryOpcode::uminus:
      case UnaryOpcode::sign:
        return cost + 3;
      case UnaryOpcode::exp:
      case UnaryOpcode::acosh:
        return cost + 210;
      case UnaryOpcode::log:
        return cost + 137;
      case UnaryOpcode::log10:
        return cost + 139;
      case UnaryOpcode::cos:
      case UnaryOpcode::sin:
        return cost + 160;
      case UnaryOpcode::tan:
        return cost + 170;
      case UnaryOpcode::acos:
      case UnaryOpcode::atan:
        return cost + 190;
      case UnaryOpcode::asin:
        return cost + 180;
      case UnaryOpcode::cosh:
      case UnaryOpcode::sinh:
      case UnaryOpcode::tanh:
      case UnaryOpcode::erf:
      case UnaryOpcode::erfc:
        return cost + 240;
      case UnaryOpcode::asinh:
        return cost + 220;
      case UnaryOpcode::atanh:
        return cost + 150;
      case UnaryOpcode::sqrt:
      case UnaryOpcode::cbrt:
      case UnaryOpcode::abs:
        return cost + 90;
      case UnaryOpcode::steadyState:
      case UnaryOpcode::steadyStateParamDeriv:
      case UnaryOpcode::steadyStateParam2ndDeriv:
      case UnaryOpcode::expectation:
        return cost;
      case UnaryOpcode::diff:
        cerr << "UnaryOpNode::cost: not implemented on UnaryOpcode::diff" << endl;
        exit(EXIT_FAILURE);
      case UnaryOpcode::adl:
        cerr << "UnaryOpNode::cost: not implemented on UnaryOpcode::adl" << endl;
        exit(EXIT_FAILURE);
      }
  exit(EXIT_FAILURE);
}

void
UnaryOpNode::computeTemporaryTerms(const pair<int, int> &derivOrder,
                                   map<pair<int, int>, temporary_terms_t> &temp_terms_map,
                                   map<expr_t, pair<int, pair<int, int>>> &reference_count,
                                   bool is_matlab) const
{
  expr_t this2 = const_cast<UnaryOpNode *>(this);
  if (auto it = reference_count.find(this2);
      it == reference_count.end())
    {
      reference_count[this2] = { 1, derivOrder };
      arg->computeTemporaryTerms(derivOrder, temp_terms_map, reference_count, is_matlab);
    }
  else
    {
      auto &[nref, min_order] = it->second;
      nref++;
      if (nref * cost(temp_terms_map, is_matlab) > min_cost(is_matlab))
        temp_terms_map[min_order].insert(this2);
    }
}

void
UnaryOpNode::computeBlockTemporaryTerms(int blk, int eq, vector<vector<temporary_terms_t>> &blocks_temporary_terms,
                                        map<expr_t, tuple<int, int, int>> &reference_count) const
{
  expr_t this2 = const_cast<UnaryOpNode *>(this);
  if (auto it = reference_count.find(this2);
      it == reference_count.end())
    {
      reference_count[this2] = { 1, blk, eq };
      arg->computeBlockTemporaryTerms(blk, eq, blocks_temporary_terms, reference_count);
    }
  else
    {
      auto &[nref, first_blk, first_eq] = it->second;
      nref++;
      if (nref * cost(blocks_temporary_terms, false) > min_cost_c)
        blocks_temporary_terms[first_blk][first_eq].insert(this2);
    }
}

bool
UnaryOpNode::containsExternalFunction() const
{
  return arg->containsExternalFunction();
}

void
UnaryOpNode::writeJsonAST(ostream &output) const
{
  output << R"({"node_type" : "UnaryOpNode", "op" : ")";
  switch (op_code)
    {
    case UnaryOpcode::uminus:
      output << "uminus";
      break;
    case UnaryOpcode::exp:
      output << "exp";
      break;
    case UnaryOpcode::log:
      output << "log";
      break;
    case UnaryOpcode::log10:
      output << "log10";
      break;
    case UnaryOpcode::cos:
      output << "cos";
      break;
    case UnaryOpcode::sin:
      output << "sin";
      break;
    case UnaryOpcode::tan:
      output << "tan";
      break;
    case UnaryOpcode::acos:
      output << "acos";
      break;
    case UnaryOpcode::asin:
      output << "asin";
      break;
    case UnaryOpcode::atan:
      output << "atan";
      break;
    case UnaryOpcode::cosh:
      output << "cosh";
      break;
    case UnaryOpcode::sinh:
      output << "sinh";
      break;
    case UnaryOpcode::tanh:
      output << "tanh";
      break;
    case UnaryOpcode::acosh:
      output << "acosh";
      break;
    case UnaryOpcode::asinh:
      output << "asinh";
      break;
    case UnaryOpcode::atanh:
      output << "atanh";
      break;
    case UnaryOpcode::sqrt:
      output << "sqrt";
      break;
    case UnaryOpcode::cbrt:
      output << "cbrt";
      break;
    case UnaryOpcode::abs:
      output << "abs";
      break;
    case UnaryOpcode::sign:
      output << "sign";
      break;
    case UnaryOpcode::diff:
      output << "diff";
      break;
    case UnaryOpcode::adl:
      output << "adl";
      break;
    case UnaryOpcode::steadyState:
      output << "steady_state";
    case UnaryOpcode::steadyStateParamDeriv:
      output << "steady_state_param_deriv";
      break;
    case UnaryOpcode::steadyStateParam2ndDeriv:
      output << "steady_state_param_second_deriv";
      break;
    case UnaryOpcode::expectation:
      output << "expectation";
      break;
    case UnaryOpcode::erf:
      output << "erf";
      break;
    case UnaryOpcode::erfc:
      output << "erfc";
      break;
    }
  output << R"(", "arg" : )";
  arg->writeJsonAST(output);
  switch (op_code)
    {
    case UnaryOpcode::adl:
      output << R"(, "adl_param_name" : ")" << adl_param_name << R"(")"
             << R"(, "lags" : [)";
      for (auto it = adl_lags.begin(); it != adl_lags.end(); ++it)
        {
          if (it != adl_lags.begin())
            output << ", ";
          output << *it;
        }
      output << "]";
      break;
    default:
      break;
    }
  output << "}";
}

void
UnaryOpNode::writeJsonOutput(ostream &output,
                             const temporary_terms_t &temporary_terms,
                             const deriv_node_temp_terms_t &tef_terms,
                             bool isdynamic) const
{
  if (temporary_terms.find(const_cast<UnaryOpNode *>(this)) != temporary_terms.end())
    {
      output << "T" << idx;
      return;
    }

  // Always put parenthesis around uminus nodes
  if (op_code == UnaryOpcode::uminus)
    output << "(";

  switch (op_code)
    {
    case UnaryOpcode::uminus:
      output << "-";
      break;
    case UnaryOpcode::exp:
      output << "exp";
      break;
    case UnaryOpcode::log:
      output << "log";
      break;
    case UnaryOpcode::log10:
      output << "log10";
      break;
    case UnaryOpcode::cos:
      output << "cos";
      break;
    case UnaryOpcode::sin:
      output << "sin";
      break;
    case UnaryOpcode::tan:
      output << "tan";
      break;
    case UnaryOpcode::acos:
      output << "acos";
      break;
    case UnaryOpcode::asin:
      output << "asin";
      break;
    case UnaryOpcode::atan:
      output << "atan";
      break;
    case UnaryOpcode::cosh:
      output << "cosh";
      break;
    case UnaryOpcode::sinh:
      output << "sinh";
      break;
    case UnaryOpcode::tanh:
      output << "tanh";
      break;
    case UnaryOpcode::acosh:
      output << "acosh";
      break;
    case UnaryOpcode::asinh:
      output << "asinh";
      break;
    case UnaryOpcode::atanh:
      output << "atanh";
      break;
    case UnaryOpcode::sqrt:
      output << "sqrt";
      break;
    case UnaryOpcode::cbrt:
      output << "cbrt";
      break;
    case UnaryOpcode::abs:
      output << "abs";
      break;
    case UnaryOpcode::sign:
      output << "sign";
      break;
    case UnaryOpcode::diff:
      output << "diff";
      break;
    case UnaryOpcode::adl:
      output << "adl(";
      arg->writeJsonOutput(output, temporary_terms, tef_terms);
      output << ", '" << adl_param_name << "', [";
      for (auto it = adl_lags.begin(); it != adl_lags.end(); ++it)
        {
          if (it != adl_lags.begin())
            output << ", ";
          output << *it;
        }
      output << "])";
      return;
    case UnaryOpcode::steadyState:
      output << "(";
      arg->writeJsonOutput(output, temporary_terms, tef_terms, isdynamic);
      output << ")";
      return;
    case UnaryOpcode::steadyStateParamDeriv:
      {
        auto varg = dynamic_cast<VariableNode *>(arg);
        assert(varg);
        assert(datatree.symbol_table.getType(varg->symb_id) == SymbolType::endogenous);
        assert(datatree.symbol_table.getType(param1_symb_id) == SymbolType::parameter);
        int tsid_endo = datatree.symbol_table.getTypeSpecificID(varg->symb_id);
        int tsid_param = datatree.symbol_table.getTypeSpecificID(param1_symb_id);
        output << "ss_param_deriv(" << tsid_endo+1 << "," << tsid_param+1 << ")";
      }
      return;
    case UnaryOpcode::steadyStateParam2ndDeriv:
      {
        auto varg = dynamic_cast<VariableNode *>(arg);
        assert(varg);
        assert(datatree.symbol_table.getType(varg->symb_id) == SymbolType::endogenous);
        assert(datatree.symbol_table.getType(param1_symb_id) == SymbolType::parameter);
        assert(datatree.symbol_table.getType(param2_symb_id) == SymbolType::parameter);
        int tsid_endo = datatree.symbol_table.getTypeSpecificID(varg->symb_id);
        int tsid_param1 = datatree.symbol_table.getTypeSpecificID(param1_symb_id);
        int tsid_param2 = datatree.symbol_table.getTypeSpecificID(param2_symb_id);
        output << "ss_param_2nd_deriv(" << tsid_endo+1 << "," << tsid_param1+1
               << "," << tsid_param2+1 << ")";
      }
      return;
    case UnaryOpcode::expectation:
      output << "EXPECTATION(" << expectation_information_set << ")";
      break;
    case UnaryOpcode::erf:
      output << "erf";
      break;
    case UnaryOpcode::erfc:
      output << "erfc";
      break;
    }

  bool close_parenthesis = false;

  /* Enclose argument with parentheses if:
     - current opcode is not uminus, or
     - current opcode is uminus and argument has lowest precedence
  */
  if (op_code != UnaryOpcode::uminus
      || (op_code == UnaryOpcode::uminus
          && arg->precedenceJson(temporary_terms) < precedenceJson(temporary_terms)))
    {
      output << "(";
      close_parenthesis = true;
    }

  // Write argument
  arg->writeJsonOutput(output, temporary_terms, tef_terms, isdynamic);

  if (close_parenthesis)
    output << ")";

  // Close parenthesis for uminus
  if (op_code == UnaryOpcode::uminus)
    output << ")";
}

void
UnaryOpNode::writeOutput(ostream &output, ExprNodeOutputType output_type,
                         const temporary_terms_t &temporary_terms,
                         const temporary_terms_idxs_t &temporary_terms_idxs,
                         const deriv_node_temp_terms_t &tef_terms) const
{
  if (checkIfTemporaryTermThenWrite(output, output_type, temporary_terms, temporary_terms_idxs))
    return;

  // Always put parenthesis around uminus nodes
  if (op_code == UnaryOpcode::uminus)
    output << LEFT_PAR(output_type);

  switch (op_code)
    {
    case UnaryOpcode::uminus:
      output << "-";
      break;
    case UnaryOpcode::exp:
      if (isLatexOutput(output_type))
        output << R"(\exp)";
      else
        output << "exp";
      break;
    case UnaryOpcode::log:
      if (isLatexOutput(output_type))
        output << R"(\log)";
      else
        output << "log";
      break;
    case UnaryOpcode::log10:
      if (isLatexOutput(output_type))
        output << R"(\log_{10})";
      else
        output << "log10";
      break;
    case UnaryOpcode::cos:
      if (isLatexOutput(output_type))
        output << R"(\cos)";
      else
        output << "cos";
      break;
    case UnaryOpcode::sin:
      if (isLatexOutput(output_type))
        output << R"(\sin)";
      else
        output << "sin";
      break;
    case UnaryOpcode::tan:
      if (isLatexOutput(output_type))
        output << R"(\tan)";
      else
        output << "tan";
      break;
    case UnaryOpcode::acos:
      output << "acos";
      break;
    case UnaryOpcode::asin:
      output << "asin";
      break;
    case UnaryOpcode::atan:
      output << "atan";
      break;
    case UnaryOpcode::cosh:
      output << "cosh";
      break;
    case UnaryOpcode::sinh:
      output << "sinh";
      break;
    case UnaryOpcode::tanh:
      output << "tanh";
      break;
    case UnaryOpcode::acosh:
      output << "acosh";
      break;
    case UnaryOpcode::asinh:
      output << "asinh";
      break;
    case UnaryOpcode::atanh:
      output << "atanh";
      break;
    case UnaryOpcode::sqrt:
      if (isLatexOutput(output_type))
        {
          output << R"(\sqrt{)";
          arg->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
          output << "}";
          return;
        }
      output << "sqrt";
      break;
    case UnaryOpcode::cbrt:
      if (isMatlabOutput(output_type))
        {
          output << "nthroot(";
          arg->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
          output << ", 3)";
          return;
        }
      else if (isLatexOutput(output_type))
        {
          output << R"(\sqrt[3]{)";
          arg->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
          output << "}";
          return;
        }
      else
        output << "cbrt";
      break;
    case UnaryOpcode::abs:
      output << "abs";
      break;
    case UnaryOpcode::sign:
      if (output_type == ExprNodeOutputType::CDynamicModel || output_type == ExprNodeOutputType::CStaticModel)
        output << "copysign";
      else
        output << "sign";
      break;
    case UnaryOpcode::steadyState:
      ExprNodeOutputType new_output_type;
      switch (output_type)
        {
        case ExprNodeOutputType::matlabDynamicModel:
        case ExprNodeOutputType::occbinDifferenceFile:
          new_output_type = ExprNodeOutputType::matlabDynamicSteadyStateOperator;
          break;
        case ExprNodeOutputType::latexDynamicModel:
          new_output_type = ExprNodeOutputType::latexDynamicSteadyStateOperator;
          break;
        case ExprNodeOutputType::CDynamicModel:
          new_output_type = ExprNodeOutputType::CDynamicSteadyStateOperator;
          break;
        case ExprNodeOutputType::juliaDynamicModel:
          new_output_type = ExprNodeOutputType::juliaDynamicSteadyStateOperator;
          break;
        default:
          new_output_type = output_type;
          break;
        }
      output << "(";
      arg->writeOutput(output, new_output_type, temporary_terms, temporary_terms_idxs, tef_terms);
      output << ")";
      return;
    case UnaryOpcode::steadyStateParamDeriv:
      {
        auto varg = dynamic_cast<VariableNode *>(arg);
        assert(varg);
        assert(datatree.symbol_table.getType(varg->symb_id) == SymbolType::endogenous);
        assert(datatree.symbol_table.getType(param1_symb_id) == SymbolType::parameter);
        int tsid_endo = datatree.symbol_table.getTypeSpecificID(varg->symb_id);
        int tsid_param = datatree.symbol_table.getTypeSpecificID(param1_symb_id);
        assert(isMatlabOutput(output_type));
        output << "ss_param_deriv(" << tsid_endo+1 << "," << tsid_param+1 << ")";
      }
      return;
    case UnaryOpcode::steadyStateParam2ndDeriv:
      {
        auto varg = dynamic_cast<VariableNode *>(arg);
        assert(varg);
        assert(datatree.symbol_table.getType(varg->symb_id) == SymbolType::endogenous);
        assert(datatree.symbol_table.getType(param1_symb_id) == SymbolType::parameter);
        assert(datatree.symbol_table.getType(param2_symb_id) == SymbolType::parameter);
        int tsid_endo = datatree.symbol_table.getTypeSpecificID(varg->symb_id);
        int tsid_param1 = datatree.symbol_table.getTypeSpecificID(param1_symb_id);
        int tsid_param2 = datatree.symbol_table.getTypeSpecificID(param2_symb_id);
        assert(isMatlabOutput(output_type));
        output << "ss_param_2nd_deriv(" << tsid_endo+1 << "," << tsid_param1+1
               << "," << tsid_param2+1 << ")";
      }
      return;
    case UnaryOpcode::expectation:
      if (!isLatexOutput(output_type))
        {
          cerr << "UnaryOpNode::writeOutput: not implemented on UnaryOpcode::expectation" << endl;
          exit(EXIT_FAILURE);
        }
      output << R"(\mathbb{E}_{t)";
      if (expectation_information_set != 0)
        {
          if (expectation_information_set > 0)
            output << "+";
          output << expectation_information_set;
        }
      output << "}";
      break;
    case UnaryOpcode::erf:
      output << "erf";
      break;
    case UnaryOpcode::erfc:
      output << "erfc";
      break;
    case UnaryOpcode::diff:
      output << "diff";
      break;
    case UnaryOpcode::adl:
      output << "adl";
      break;
    }

  bool close_parenthesis = false;

  /* Enclose argument with parentheses if:
     - current opcode is not uminus, or
     - current opcode is uminus and argument has lowest precedence
  */
  if (op_code != UnaryOpcode::uminus
      || (op_code == UnaryOpcode::uminus
          && arg->precedence(output_type, temporary_terms) < precedence(output_type, temporary_terms)))
    {
      output << LEFT_PAR(output_type);
      if (op_code == UnaryOpcode::sign && (output_type == ExprNodeOutputType::CDynamicModel || output_type == ExprNodeOutputType::CStaticModel))
        output << "1.0,";
      close_parenthesis = true;
    }

  // Write argument
  arg->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);

  if (close_parenthesis)
    output << RIGHT_PAR(output_type);

  // Close parenthesis for uminus
  if (op_code == UnaryOpcode::uminus)
    output << RIGHT_PAR(output_type);
}

void
UnaryOpNode::writeExternalFunctionOutput(ostream &output, ExprNodeOutputType output_type,
                                         const temporary_terms_t &temporary_terms,
                                         const temporary_terms_idxs_t &temporary_terms_idxs,
                                         deriv_node_temp_terms_t &tef_terms) const
{
  arg->writeExternalFunctionOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
}

void
UnaryOpNode::writeJsonExternalFunctionOutput(vector<string> &efout,
                                             const temporary_terms_t &temporary_terms,
                                             deriv_node_temp_terms_t &tef_terms,
                                             bool isdynamic) const
{
  arg->writeJsonExternalFunctionOutput(efout, temporary_terms, tef_terms, isdynamic);
}

void
UnaryOpNode::compileExternalFunctionOutput(ostream &CompileCode, unsigned int &instruction_number,
                                           bool lhs_rhs, const temporary_terms_t &temporary_terms,
                                           const temporary_terms_idxs_t &temporary_terms_idxs, bool dynamic, bool steady_dynamic,
                                           deriv_node_temp_terms_t &tef_terms) const
{
  arg->compileExternalFunctionOutput(CompileCode, instruction_number, lhs_rhs, temporary_terms,
                                     temporary_terms_idxs, dynamic, steady_dynamic, tef_terms);
}

double
UnaryOpNode::eval_opcode(UnaryOpcode op_code, double v) noexcept(false)
{
  switch (op_code)
    {
    case UnaryOpcode::uminus:
      return -v;
    case UnaryOpcode::exp:
      return exp(v);
    case UnaryOpcode::log:
      return log(v);
    case UnaryOpcode::log10:
      return log10(v);
    case UnaryOpcode::cos:
      return cos(v);
    case UnaryOpcode::sin:
      return sin(v);
    case UnaryOpcode::tan:
      return tan(v);
    case UnaryOpcode::acos:
      return acos(v);
    case UnaryOpcode::asin:
      return asin(v);
    case UnaryOpcode::atan:
      return atan(v);
    case UnaryOpcode::cosh:
      return cosh(v);
    case UnaryOpcode::sinh:
      return sinh(v);
    case UnaryOpcode::tanh:
      return tanh(v);
    case UnaryOpcode::acosh:
      return acosh(v);
    case UnaryOpcode::asinh:
      return asinh(v);
    case UnaryOpcode::atanh:
      return atanh(v);
    case UnaryOpcode::sqrt:
      return sqrt(v);
    case UnaryOpcode::cbrt:
      return cbrt(v);
    case UnaryOpcode::abs:
      return abs(v);
    case UnaryOpcode::sign:
      return (v > 0) ? 1 : ((v < 0) ? -1 : 0);
    case UnaryOpcode::steadyState:
      return v;
    case UnaryOpcode::steadyStateParamDeriv:
      cerr << "UnaryOpNode::eval_opcode: not implemented on UnaryOpcode::steadyStateParamDeriv" << endl;
      exit(EXIT_FAILURE);
    case UnaryOpcode::steadyStateParam2ndDeriv:
      cerr << "UnaryOpNode::eval_opcode: not implemented on UnaryOpcode::steadyStateParam2ndDeriv" << endl;
      exit(EXIT_FAILURE);
    case UnaryOpcode::expectation:
      cerr << "UnaryOpNode::eval_opcode: not implemented on UnaryOpcode::expectation" << endl;
      exit(EXIT_FAILURE);
    case UnaryOpcode::erf:
      return erf(v);
    case UnaryOpcode::erfc:
      return erfc(v);
    case UnaryOpcode::diff:
      cerr << "UnaryOpNode::eval_opcode: not implemented on UnaryOpcode::diff" << endl;
      exit(EXIT_FAILURE);
    case UnaryOpcode::adl:
      cerr << "UnaryOpNode::eval_opcode: not implemented on UnaryOpcode::adl" << endl;
      exit(EXIT_FAILURE);
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

double
UnaryOpNode::eval(const eval_context_t &eval_context) const noexcept(false)
{
  double v = arg->eval(eval_context);
  return eval_opcode(op_code, v);
}

void
UnaryOpNode::compile(ostream &CompileCode, unsigned int &instruction_number,
                     bool lhs_rhs, const temporary_terms_t &temporary_terms,
                     const temporary_terms_idxs_t &temporary_terms_idxs, bool dynamic, bool steady_dynamic,
                     const deriv_node_temp_terms_t &tef_terms) const
{
  if (auto this2 = const_cast<UnaryOpNode *>(this);
      temporary_terms.find(this2) != temporary_terms.end())
    {
      if (dynamic)
        {
          FLDT_ fldt(temporary_terms_idxs.at(this2));
          fldt.write(CompileCode, instruction_number);
        }
      else
        {
          FLDST_ fldst(temporary_terms_idxs.at(this2));
          fldst.write(CompileCode, instruction_number);
        }
      return;
    }
  if (op_code == UnaryOpcode::steadyState)
    arg->compile(CompileCode, instruction_number, lhs_rhs, temporary_terms, temporary_terms_idxs, dynamic, true, tef_terms);
  else
    {
      arg->compile(CompileCode, instruction_number, lhs_rhs, temporary_terms, temporary_terms_idxs, dynamic, steady_dynamic, tef_terms);
      FUNARY_ funary{static_cast<uint8_t>(op_code)};
      funary.write(CompileCode, instruction_number);
    }
}

void
UnaryOpNode::collectVARLHSVariable(set<expr_t> &result) const
{
  if (op_code == UnaryOpcode::diff)
    result.insert(const_cast<UnaryOpNode *>(this));
  else
    arg->collectVARLHSVariable(result);
}

void
UnaryOpNode::collectDynamicVariables(SymbolType type_arg, set<pair<int, int>> &result) const
{
  arg->collectDynamicVariables(type_arg, result);
}

void
UnaryOpNode::computeSubExprContainingVariable(int symb_id, int lag, set<expr_t> &contain_var) const
{
  arg->computeSubExprContainingVariable(symb_id, lag, contain_var);
  if (contain_var.count(arg) > 0)
    contain_var.insert(const_cast<UnaryOpNode *>(this));
}

BinaryOpNode *
UnaryOpNode::normalizeEquationHelper(const set<expr_t> &contain_var, expr_t rhs) const
{
  assert(contain_var.count(const_cast<UnaryOpNode *>(this)) > 0);

  switch (op_code)
    {
    case UnaryOpcode::uminus:
      rhs = datatree.AddUMinus(rhs);
      break;
    case UnaryOpcode::exp:
      rhs = datatree.AddLog(rhs);
      break;
    case UnaryOpcode::log:
      rhs = datatree.AddExp(rhs);
      break;
    case UnaryOpcode::log10:
      rhs = datatree.AddPower(datatree.AddNonNegativeConstant("10"), rhs);
      break;
    case UnaryOpcode::cos:
      rhs = datatree.AddAcos(rhs);
      break;
    case UnaryOpcode::sin:
      rhs = datatree.AddAsin(rhs);
      break;
    case UnaryOpcode::tan:
      rhs = datatree.AddAtan(rhs);
      break;
    case UnaryOpcode::acos:
      rhs = datatree.AddCos(rhs);
      break;
    case UnaryOpcode::asin:
      rhs = datatree.AddSin(rhs);
      break;
    case UnaryOpcode::atan:
      rhs = datatree.AddTan(rhs);
      break;
    case UnaryOpcode::cosh:
      rhs = datatree.AddAcosh(rhs);
      break;
    case UnaryOpcode::sinh:
      rhs = datatree.AddAsinh(rhs);
      break;
    case UnaryOpcode::tanh:
      rhs = datatree.AddAtanh(rhs);
      break;
    case UnaryOpcode::acosh:
      rhs = datatree.AddCosh(rhs);
      break;
    case UnaryOpcode::asinh:
      rhs = datatree.AddSinh(rhs);
      break;
    case UnaryOpcode::atanh:
      rhs = datatree.AddTanh(rhs);
      break;
    case UnaryOpcode::sqrt:
      rhs = datatree.AddPower(rhs, datatree.Two);
      break;
    case UnaryOpcode::cbrt:
      rhs = datatree.AddPower(rhs, datatree.Three);
      break;
    default:
      throw NormalizationFailed();
    }

  return arg->normalizeEquationHelper(contain_var, rhs);
}

expr_t
UnaryOpNode::getChainRuleDerivative(int deriv_id, const map<int, BinaryOpNode *> &recursive_variables)
{
  expr_t darg = arg->getChainRuleDerivative(deriv_id, recursive_variables);
  return composeDerivatives(darg, deriv_id);
}

expr_t
UnaryOpNode::buildSimilarUnaryOpNode(expr_t alt_arg, DataTree &alt_datatree) const
{
  switch (op_code)
    {
    case UnaryOpcode::uminus:
      return alt_datatree.AddUMinus(alt_arg);
    case UnaryOpcode::exp:
      return alt_datatree.AddExp(alt_arg);
    case UnaryOpcode::log:
      return alt_datatree.AddLog(alt_arg);
    case UnaryOpcode::log10:
      return alt_datatree.AddLog10(alt_arg);
    case UnaryOpcode::cos:
      return alt_datatree.AddCos(alt_arg);
    case UnaryOpcode::sin:
      return alt_datatree.AddSin(alt_arg);
    case UnaryOpcode::tan:
      return alt_datatree.AddTan(alt_arg);
    case UnaryOpcode::acos:
      return alt_datatree.AddAcos(alt_arg);
    case UnaryOpcode::asin:
      return alt_datatree.AddAsin(alt_arg);
    case UnaryOpcode::atan:
      return alt_datatree.AddAtan(alt_arg);
    case UnaryOpcode::cosh:
      return alt_datatree.AddCosh(alt_arg);
    case UnaryOpcode::sinh:
      return alt_datatree.AddSinh(alt_arg);
    case UnaryOpcode::tanh:
      return alt_datatree.AddTanh(alt_arg);
    case UnaryOpcode::acosh:
      return alt_datatree.AddAcosh(alt_arg);
    case UnaryOpcode::asinh:
      return alt_datatree.AddAsinh(alt_arg);
    case UnaryOpcode::atanh:
      return alt_datatree.AddAtanh(alt_arg);
    case UnaryOpcode::sqrt:
      return alt_datatree.AddSqrt(alt_arg);
    case UnaryOpcode::cbrt:
      return alt_datatree.AddCbrt(alt_arg);
    case UnaryOpcode::abs:
      return alt_datatree.AddAbs(alt_arg);
    case UnaryOpcode::sign:
      return alt_datatree.AddSign(alt_arg);
    case UnaryOpcode::steadyState:
      return alt_datatree.AddSteadyState(alt_arg);
    case UnaryOpcode::steadyStateParamDeriv:
      cerr << "UnaryOpNode::buildSimilarUnaryOpNode: UnaryOpcode::steadyStateParamDeriv can't be translated" << endl;
      exit(EXIT_FAILURE);
    case UnaryOpcode::steadyStateParam2ndDeriv:
      cerr << "UnaryOpNode::buildSimilarUnaryOpNode: UnaryOpcode::steadyStateParam2ndDeriv can't be translated" << endl;
      exit(EXIT_FAILURE);
    case UnaryOpcode::expectation:
      return alt_datatree.AddExpectation(expectation_information_set, alt_arg);
    case UnaryOpcode::erf:
      return alt_datatree.AddErf(alt_arg);
    case UnaryOpcode::erfc:
      return alt_datatree.AddErfc(alt_arg);
    case UnaryOpcode::diff:
      return alt_datatree.AddDiff(alt_arg);
    case UnaryOpcode::adl:
      return alt_datatree.AddAdl(alt_arg, adl_param_name, adl_lags);
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

expr_t
UnaryOpNode::toStatic(DataTree &static_datatree) const
{
  expr_t sarg = arg->toStatic(static_datatree);
  return buildSimilarUnaryOpNode(sarg, static_datatree);
}

void
UnaryOpNode::computeXrefs(EquationInfo &ei) const
{
  arg->computeXrefs(ei);
}

expr_t
UnaryOpNode::clone(DataTree &datatree) const
{
  expr_t substarg = arg->clone(datatree);
  return buildSimilarUnaryOpNode(substarg, datatree);
}

int
UnaryOpNode::maxEndoLead() const
{
  return arg->maxEndoLead();
}

int
UnaryOpNode::maxExoLead() const
{
  return arg->maxExoLead();
}

int
UnaryOpNode::maxEndoLag() const
{
  return arg->maxEndoLag();
}

int
UnaryOpNode::maxExoLag() const
{
  return arg->maxExoLag();
}

int
UnaryOpNode::maxLead() const
{
  return arg->maxLead();
}

int
UnaryOpNode::maxLag() const
{
  return arg->maxLag();
}

int
UnaryOpNode::maxLagWithDiffsExpanded() const
{
  if (op_code == UnaryOpcode::diff)
    return arg->maxLagWithDiffsExpanded() + 1;
  return arg->maxLagWithDiffsExpanded();
}

expr_t
UnaryOpNode::undiff() const
{
  if (op_code == UnaryOpcode::diff)
    return arg;
  return arg->undiff();
}

int
UnaryOpNode::VarMaxLag(const set<expr_t> &lhs_lag_equiv) const
{
  auto [lag_equiv_repr, index] = getLagEquivalenceClass();
  if (lhs_lag_equiv.find(lag_equiv_repr) == lhs_lag_equiv.end())
    return 0;
  return arg->maxLag();
}

expr_t
UnaryOpNode::substituteAdl() const
{
  if (op_code != UnaryOpcode::adl)
    {
      expr_t argsubst = arg->substituteAdl();
      return buildSimilarUnaryOpNode(argsubst, datatree);
    }

  expr_t arg1subst = arg->substituteAdl();
  expr_t retval = nullptr;
  ostringstream inttostr;

  for (auto it = adl_lags.begin(); it != adl_lags.end(); ++it)
    if (it == adl_lags.begin())
      {
        inttostr << *it;
        retval = datatree.AddTimes(datatree.AddVariable(datatree.symbol_table.getID(adl_param_name + "_lag_" + inttostr.str()), 0),
                                   arg1subst->decreaseLeadsLags(*it));
      }
    else
      {
        inttostr.clear();
        inttostr.str("");
        inttostr << *it;
        retval = datatree.AddPlus(retval,
                                  datatree.AddTimes(datatree.AddVariable(datatree.symbol_table.getID(adl_param_name + "_lag_"
                                                                                                     + inttostr.str()), 0),
                                                    arg1subst->decreaseLeadsLags(*it)));
      }
  return retval;
}

expr_t
UnaryOpNode::substituteModelLocalVariables() const
{
  expr_t argsubst = arg->substituteModelLocalVariables();
  return buildSimilarUnaryOpNode(argsubst, datatree);
}

expr_t
UnaryOpNode::substituteVarExpectation(const map<string, expr_t> &subst_table) const
{
  expr_t argsubst = arg->substituteVarExpectation(subst_table);
  return buildSimilarUnaryOpNode(argsubst, datatree);
}

int
UnaryOpNode::countDiffs() const
{
  if (op_code == UnaryOpcode::diff)
    return arg->countDiffs() + 1;
  return arg->countDiffs();
}

bool
UnaryOpNode::createAuxVarForUnaryOpNode() const
{
  switch (op_code)
    {
    case UnaryOpcode::exp:
    case UnaryOpcode::log:
    case UnaryOpcode::log10:
    case UnaryOpcode::cos:
    case UnaryOpcode::sin:
    case UnaryOpcode::tan:
    case UnaryOpcode::acos:
    case UnaryOpcode::asin:
    case UnaryOpcode::atan:
    case UnaryOpcode::cosh:
    case UnaryOpcode::sinh:
    case UnaryOpcode::tanh:
    case UnaryOpcode::acosh:
    case UnaryOpcode::asinh:
    case UnaryOpcode::atanh:
    case UnaryOpcode::sqrt:
    case UnaryOpcode::cbrt:
    case UnaryOpcode::abs:
    case UnaryOpcode::sign:
    case UnaryOpcode::erf:
    case UnaryOpcode::erfc:
      return true;
    default:
      return false;
    }
}

void
UnaryOpNode::findUnaryOpNodesForAuxVarCreation(lag_equivalence_table_t &nodes) const
{
  arg->findUnaryOpNodesForAuxVarCreation(nodes);

  if (!this->createAuxVarForUnaryOpNode())
    return;

  auto [lag_equiv_repr, index] = getLagEquivalenceClass();
  nodes[lag_equiv_repr][index] = const_cast<UnaryOpNode *>(this);
}

void
UnaryOpNode::findDiffNodes(lag_equivalence_table_t &nodes) const
{
  arg->findDiffNodes(nodes);

  if (op_code != UnaryOpcode::diff)
    return;

  auto [lag_equiv_repr, index] = getLagEquivalenceClass();
  nodes[lag_equiv_repr][index] = const_cast<UnaryOpNode *>(this);
}

int
UnaryOpNode::findTargetVariable(int lhs_symb_id) const
{
  return arg->findTargetVariable(lhs_symb_id);
}

expr_t
UnaryOpNode::substituteDiff(const lag_equivalence_table_t &nodes, subst_table_t &subst_table,
                            vector<BinaryOpNode *> &neweqs) const
{
  // If this is not a diff node, then substitute recursively and return
  expr_t argsubst = arg->substituteDiff(nodes, subst_table, neweqs);
  if (op_code != UnaryOpcode::diff)
    return buildSimilarUnaryOpNode(argsubst, datatree);

  if (auto sit = subst_table.find(this);
      sit != subst_table.end())
    return const_cast<VariableNode *>(sit->second);

  auto [lag_equiv_repr, index] = getLagEquivalenceClass();
  auto it = nodes.find(lag_equiv_repr);
  if (it == nodes.end() || it->second.find(index) == it->second.end()
      || it->second.at(index) != this)
    {
      /* diff does not appear in VAR equations, so simply create aux var and return.
         Once the comparison of expression nodes works, come back and remove
         this part, folding into the next loop. */
      int symb_id = datatree.symbol_table.addDiffAuxiliaryVar(argsubst->idx, const_cast<UnaryOpNode *>(this));
      VariableNode *aux_var = datatree.AddVariable(symb_id, 0);
      neweqs.push_back(datatree.AddEqual(aux_var,
                                         datatree.AddMinus(argsubst,
                                                           argsubst->decreaseLeadsLags(1))));
      subst_table[this] = dynamic_cast<VariableNode *>(aux_var);
      return const_cast<VariableNode *>(subst_table[this]);
    }

  /* At this point, we know that this node (and its lagged/leaded brothers)
     must be substituted. We create the auxiliary variable and fill the
     substitution table for all those similar nodes, in an iteration going from
     leads to lags. */
  int last_index = 0;
  VariableNode *last_aux_var = nullptr;
  for (auto rit = it->second.rbegin(); rit != it->second.rend(); ++rit)
    {
      expr_t argsubst = dynamic_cast<UnaryOpNode *>(rit->second)->
        arg->substituteDiff(nodes, subst_table, neweqs);
      auto vn = dynamic_cast<VariableNode *>(argsubst);
      int symb_id;
      if (rit == it->second.rbegin())
        {
          if (vn)
            symb_id = datatree.symbol_table.addDiffAuxiliaryVar(argsubst->idx, rit->second, vn->symb_id, vn->lag);
          else
            symb_id = datatree.symbol_table.addDiffAuxiliaryVar(argsubst->idx, rit->second);

          // make originating aux var & equation
          last_index = rit->first;
          last_aux_var = datatree.AddVariable(symb_id, 0);
          //ORIG_AUX_DIFF = argsubst - argsubst(-1)
          neweqs.push_back(datatree.AddEqual(last_aux_var,
                                             datatree.AddMinus(argsubst,
                                                               argsubst->decreaseLeadsLags(1))));
          subst_table[rit->second] = dynamic_cast<VariableNode *>(last_aux_var);
        }
      else
        {
          // just add equation of form: AUX_DIFF = LAST_AUX_VAR(-1)
          VariableNode *new_aux_var = nullptr;
          for (int i = last_index; i > rit->first; i--)
            {
              if (i == last_index)
                symb_id = datatree.symbol_table.addDiffLagAuxiliaryVar(argsubst->idx, rit->second,
                                                                       last_aux_var->symb_id, last_aux_var->lag - 1);
              else
                symb_id = datatree.symbol_table.addDiffLagAuxiliaryVar(new_aux_var->idx, rit->second,
                                                                       last_aux_var->symb_id, last_aux_var->lag - 1);

              new_aux_var = datatree.AddVariable(symb_id, 0);
              neweqs.push_back(datatree.AddEqual(new_aux_var,
                                                 last_aux_var->decreaseLeadsLags(1)));
              last_aux_var = new_aux_var;
            }
          subst_table[rit->second] = dynamic_cast<VariableNode *>(new_aux_var);
          last_index = rit->first;
        }
    }
  return const_cast<VariableNode *>(subst_table[this]);
}

expr_t
UnaryOpNode::substituteUnaryOpNodes(const lag_equivalence_table_t &nodes, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  if (auto sit = subst_table.find(this);
      sit != subst_table.end())
    return const_cast<VariableNode *>(sit->second);

  /* If the equivalence class of this node is not marked for substitution,
     then substitute recursively and return. */
  auto [lag_equiv_repr, index] = getLagEquivalenceClass();
  auto it = nodes.find(lag_equiv_repr);
  expr_t argsubst = arg->substituteUnaryOpNodes(nodes, subst_table, neweqs);
  if (it == nodes.end())
    return buildSimilarUnaryOpNode(argsubst, datatree);

  string unary_op;
  switch (op_code)
    {
    case UnaryOpcode::exp:
      unary_op = "exp";
      break;
    case UnaryOpcode::log:
      unary_op = "log";
      break;
    case UnaryOpcode::log10:
      unary_op = "log10";
      break;
    case UnaryOpcode::cos:
      unary_op = "cos";
      break;
    case UnaryOpcode::sin:
      unary_op = "sin";
      break;
    case UnaryOpcode::tan:
      unary_op = "tan";
      break;
    case UnaryOpcode::acos:
      unary_op = "acos";
      break;
    case UnaryOpcode::asin:
      unary_op = "asin";
      break;
    case UnaryOpcode::atan:
      unary_op = "atan";
      break;
    case UnaryOpcode::cosh:
      unary_op = "cosh";
      break;
    case UnaryOpcode::sinh:
      unary_op = "sinh";
      break;
    case UnaryOpcode::tanh:
      unary_op = "tanh";
      break;
    case UnaryOpcode::acosh:
      unary_op = "acosh";
      break;
    case UnaryOpcode::asinh:
      unary_op = "asinh";
      break;
    case UnaryOpcode::atanh:
      unary_op = "atanh";
      break;
    case UnaryOpcode::sqrt:
      unary_op = "sqrt";
      break;
    case UnaryOpcode::cbrt:
      unary_op = "cbrt";
      break;
    case UnaryOpcode::abs:
      unary_op = "abs";
      break;
    case UnaryOpcode::sign:
      unary_op = "sign";
      break;
    case UnaryOpcode::erf:
      unary_op = "erf";
      break;
    case UnaryOpcode::erfc:
      unary_op = "erfc";
      break;
    default:
      cerr << "UnaryOpNode::substituteUnaryOpNodes: Shouldn't arrive here" << endl;
      exit(EXIT_FAILURE);
    }

  /* At this point, we know that this node (and its lagged/leaded brothers)
     must be substituted. We create the auxiliary variable and fill the
     substitution table for all those similar nodes, in an iteration going from
     leads to lags. */
  int base_index = it->second.rbegin()->first; // Within the equivalence class,
                                               // index of the node that will
                                               // be used as the definition for
                                               // the aux var.
  VariableNode *aux_var = nullptr;
  for (auto rit = it->second.rbegin(); rit != it->second.rend(); ++rit)
    if (rit == it->second.rbegin())
      {
        /* Verify that we’re not operating on a node with leads, since the
           transformation does not take into account the expectation operator. We only
           need to do this for the first iteration of the loop, because we’re
           going from leads to lags. */
        if (rit->second->maxLead() > 0)
          {
            cerr << "Cannot substitute unary operations that contain leads" << endl;
            exit(EXIT_FAILURE);
          }

        auto argsubst_shifted = argsubst->decreaseLeadsLags(index - base_index);
        auto aux_def = buildSimilarUnaryOpNode(argsubst_shifted, datatree);
        int symb_id;
        if (auto vn = dynamic_cast<VariableNode *>(argsubst_shifted); !vn)
          symb_id = datatree.symbol_table.addUnaryOpAuxiliaryVar(this->idx, aux_def, unary_op);
        else
          symb_id = datatree.symbol_table.addUnaryOpAuxiliaryVar(this->idx, aux_def, unary_op,
                                                                 vn->symb_id, vn->lag);
        aux_var = datatree.AddVariable(symb_id, 0);
        neweqs.push_back(datatree.AddEqual(aux_var, aux_def));
        subst_table[rit->second] = dynamic_cast<VariableNode *>(aux_var);
      }
    else
      subst_table[rit->second] = dynamic_cast<VariableNode *>(aux_var->decreaseLeadsLags(base_index - rit->first));

  assert(subst_table.find(this) != subst_table.end());

  return const_cast<VariableNode *>(subst_table.at(this));
}

expr_t
UnaryOpNode::substitutePacExpectation(const string &name, expr_t subexpr)
{
  expr_t argsubst = arg->substitutePacExpectation(name, subexpr);
  return buildSimilarUnaryOpNode(argsubst, datatree);
}

expr_t
UnaryOpNode::substitutePacTargetNonstationary(const string &name, expr_t subexpr)
{
  expr_t argsubst = arg->substitutePacTargetNonstationary(name, subexpr);
  return buildSimilarUnaryOpNode(argsubst, datatree);
}

expr_t
UnaryOpNode::decreaseLeadsLags(int n) const
{
  expr_t argsubst = arg->decreaseLeadsLags(n);
  return buildSimilarUnaryOpNode(argsubst, datatree);
}

expr_t
UnaryOpNode::decreaseLeadsLagsPredeterminedVariables() const
{
  expr_t argsubst = arg->decreaseLeadsLagsPredeterminedVariables();
  return buildSimilarUnaryOpNode(argsubst, datatree);
}

expr_t
UnaryOpNode::substituteEndoLeadGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const
{
  if (op_code == UnaryOpcode::uminus || deterministic_model)
    {
      expr_t argsubst = arg->substituteEndoLeadGreaterThanTwo(subst_table, neweqs, deterministic_model);
      return buildSimilarUnaryOpNode(argsubst, datatree);
    }
  else
    {
      if (maxEndoLead() >= 2)
        return createEndoLeadAuxiliaryVarForMyself(subst_table, neweqs);
      else
        return const_cast<UnaryOpNode *>(this);
    }
}

expr_t
UnaryOpNode::substituteEndoLagGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  expr_t argsubst = arg->substituteEndoLagGreaterThanTwo(subst_table, neweqs);
  return buildSimilarUnaryOpNode(argsubst, datatree);
}

expr_t
UnaryOpNode::substituteExoLead(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const
{
  if (op_code == UnaryOpcode::uminus || deterministic_model)
    {
      expr_t argsubst = arg->substituteExoLead(subst_table, neweqs, deterministic_model);
      return buildSimilarUnaryOpNode(argsubst, datatree);
    }
  else
    {
      if (maxExoLead() >= 1)
        return createExoLeadAuxiliaryVarForMyself(subst_table, neweqs);
      else
        return const_cast<UnaryOpNode *>(this);
    }
}

expr_t
UnaryOpNode::substituteExoLag(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  expr_t argsubst = arg->substituteExoLag(subst_table, neweqs);
  return buildSimilarUnaryOpNode(argsubst, datatree);
}

expr_t
UnaryOpNode::substituteExpectation(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool partial_information_model) const
{
  if (op_code == UnaryOpcode::expectation)
    {
      if (auto it = subst_table.find(const_cast<UnaryOpNode *>(this)); it != subst_table.end())
        return const_cast<VariableNode *>(it->second);

      //Arriving here, we need to create an auxiliary variable for this Expectation Operator:
      //AUX_EXPECT_(LEAD/LAG)_(period)_(arg.idx) OR
      //AUX_EXPECT_(info_set_name)_(arg.idx)
      int symb_id = datatree.symbol_table.addExpectationAuxiliaryVar(expectation_information_set, arg->idx, const_cast<UnaryOpNode *>(this));
      expr_t newAuxE = datatree.AddVariable(symb_id, 0);

      if (partial_information_model && expectation_information_set == 0)
        if (!dynamic_cast<VariableNode *>(arg))
          {
            cerr << "ERROR: In Partial Information models, EXPECTATION(0)(X) "
                 << "can only be used when X is a single variable." << endl;
            exit(EXIT_FAILURE);
          }

      //take care of any nested expectation operators by calling arg->substituteExpectation(.), then decreaseLeadsLags for this UnaryOpcode::expectation operator
      //arg(lag-period) (holds entire subtree of arg(lag-period)
      expr_t substexpr = (arg->substituteExpectation(subst_table, neweqs, partial_information_model))->decreaseLeadsLags(expectation_information_set);
      assert(substexpr);
      neweqs.push_back(datatree.AddEqual(newAuxE, substexpr)); //AUXE_period_arg.idx = arg(lag-period)
      newAuxE = datatree.AddVariable(symb_id, expectation_information_set);

      assert(dynamic_cast<VariableNode *>(newAuxE));
      subst_table[this] = dynamic_cast<VariableNode *>(newAuxE);
      return newAuxE;
    }
  else
    {
      expr_t argsubst = arg->substituteExpectation(subst_table, neweqs, partial_information_model);
      return buildSimilarUnaryOpNode(argsubst, datatree);
    }
}

expr_t
UnaryOpNode::differentiateForwardVars(const vector<string> &subset, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  expr_t argsubst = arg->differentiateForwardVars(subset, subst_table, neweqs);
  return buildSimilarUnaryOpNode(argsubst, datatree);
}

bool
UnaryOpNode::isNumConstNodeEqualTo(double value) const
{
  return false;
}

bool
UnaryOpNode::isVariableNodeEqualTo(SymbolType type_arg, int variable_id, int lag_arg) const
{
  return false;
}

bool
UnaryOpNode::containsPacExpectation(const string &pac_model_name) const
{
  return arg->containsPacExpectation(pac_model_name);
}

bool
UnaryOpNode::containsPacTargetNonstationary(const string &pac_model_name) const
{
  return arg->containsPacTargetNonstationary(pac_model_name);
}

expr_t
UnaryOpNode::replaceTrendVar() const
{
  expr_t argsubst = arg->replaceTrendVar();
  return buildSimilarUnaryOpNode(argsubst, datatree);
}

expr_t
UnaryOpNode::detrend(int symb_id, bool log_trend, expr_t trend) const
{
  expr_t argsubst = arg->detrend(symb_id, log_trend, trend);
  return buildSimilarUnaryOpNode(argsubst, datatree);
}

expr_t
UnaryOpNode::removeTrendLeadLag(const map<int, expr_t> &trend_symbols_map) const
{
  expr_t argsubst = arg->removeTrendLeadLag(trend_symbols_map);
  return buildSimilarUnaryOpNode(argsubst, datatree);
}

bool
UnaryOpNode::isInStaticForm() const
{
  if (op_code == UnaryOpcode::steadyState || op_code == UnaryOpcode::steadyStateParamDeriv
      || op_code == UnaryOpcode::steadyStateParam2ndDeriv
      || op_code == UnaryOpcode::expectation)
    return false;
  else
    return arg->isInStaticForm();
}

bool
UnaryOpNode::isParamTimesEndogExpr() const
{
  return arg->isParamTimesEndogExpr();
}

expr_t
UnaryOpNode::substituteStaticAuxiliaryVariable() const
{
  if (op_code == UnaryOpcode::diff)
    return datatree.Zero;

  expr_t argsubst = arg->substituteStaticAuxiliaryVariable();
  if (op_code == UnaryOpcode::expectation)
    return argsubst;
  else
    return buildSimilarUnaryOpNode(argsubst, datatree);
}

expr_t
UnaryOpNode::replaceVarsInEquation(map<VariableNode *, NumConstNode *> &table) const
{
  expr_t argsubst = arg->replaceVarsInEquation(table);
  return buildSimilarUnaryOpNode(argsubst, datatree);
}

BinaryOpNode::BinaryOpNode(DataTree &datatree_arg, int idx_arg, const expr_t arg1_arg,
                           BinaryOpcode op_code_arg, const expr_t arg2_arg, int powerDerivOrder_arg) :
  ExprNode{datatree_arg, idx_arg},
  arg1{arg1_arg},
  arg2{arg2_arg},
  op_code{op_code_arg},
  powerDerivOrder{powerDerivOrder_arg}
{
  assert(powerDerivOrder >= 0);
}

void
BinaryOpNode::prepareForDerivation()
{
  if (preparedForDerivation)
    return;

  preparedForDerivation = true;

  arg1->prepareForDerivation();
  arg2->prepareForDerivation();

  // Non-null derivatives are the union of those of the arguments
  // Compute set union of arg1->non_null_derivatives and arg2->non_null_derivatives
  set_union(arg1->non_null_derivatives.begin(),
            arg1->non_null_derivatives.end(),
            arg2->non_null_derivatives.begin(),
            arg2->non_null_derivatives.end(),
            inserter(non_null_derivatives, non_null_derivatives.begin()));
}

expr_t
BinaryOpNode::getNonZeroPartofEquation() const
{
  assert(arg1 == datatree.Zero || arg2 == datatree.Zero);
  if (arg1 == datatree.Zero)
    return arg2;
  return arg1;
}

expr_t
BinaryOpNode::composeDerivatives(expr_t darg1, expr_t darg2)
{
  expr_t t11, t12, t13, t14, t15;

  switch (op_code)
    {
    case BinaryOpcode::plus:
      return datatree.AddPlus(darg1, darg2);
    case BinaryOpcode::minus:
    case BinaryOpcode::equal:
      return datatree.AddMinus(darg1, darg2);
    case BinaryOpcode::times:
      t11 = datatree.AddTimes(darg1, arg2);
      t12 = datatree.AddTimes(darg2, arg1);
      return datatree.AddPlus(t11, t12);
    case BinaryOpcode::divide:
      if (darg2 != datatree.Zero)
        {
          t11 = datatree.AddTimes(darg1, arg2);
          t12 = datatree.AddTimes(darg2, arg1);
          t13 = datatree.AddMinus(t11, t12);
          t14 = datatree.AddTimes(arg2, arg2);
          return datatree.AddDivide(t13, t14);
        }
      else
        return datatree.AddDivide(darg1, arg2);
    case BinaryOpcode::less:
    case BinaryOpcode::greater:
    case BinaryOpcode::lessEqual:
    case BinaryOpcode::greaterEqual:
    case BinaryOpcode::equalEqual:
    case BinaryOpcode::different:
      return datatree.Zero;
    case BinaryOpcode::power:
      if (darg2 == datatree.Zero)
        if (darg1 == datatree.Zero)
          return datatree.Zero;
        else
          if (dynamic_cast<NumConstNode *>(arg2))
            {
              t11 = datatree.AddMinus(arg2, datatree.One);
              t12 = datatree.AddPower(arg1, t11);
              t13 = datatree.AddTimes(arg2, t12);
              return datatree.AddTimes(darg1, t13);
            }
          else
            return datatree.AddTimes(darg1, datatree.AddPowerDeriv(arg1, arg2, powerDerivOrder + 1));
      else
        {
          t11 = datatree.AddLog(arg1);
          t12 = datatree.AddTimes(darg2, t11);
          t13 = datatree.AddTimes(darg1, arg2);
          t14 = datatree.AddDivide(t13, arg1);
          t15 = datatree.AddPlus(t12, t14);
          return datatree.AddTimes(t15, this);
        }
    case BinaryOpcode::powerDeriv:
      if (darg2 == datatree.Zero)
        return datatree.AddTimes(darg1, datatree.AddPowerDeriv(arg1, arg2, powerDerivOrder + 1));
      else
        {
          t11 = datatree.AddTimes(darg2, datatree.AddLog(arg1));
          t12 = datatree.AddMinus(arg2, datatree.AddPossiblyNegativeConstant(powerDerivOrder));
          t13 = datatree.AddTimes(darg1, t12);
          t14 = datatree.AddDivide(t13, arg1);
          t15 = datatree.AddPlus(t11, t14);
          expr_t f = datatree.AddPower(arg1, t12);
          expr_t first_part = datatree.AddTimes(f, t15);

          for (int i = 0; i < powerDerivOrder; i++)
            first_part = datatree.AddTimes(first_part, datatree.AddMinus(arg2, datatree.AddPossiblyNegativeConstant(i)));

          t13 = datatree.Zero;
          for (int i = 0; i < powerDerivOrder; i++)
            {
              t11 = datatree.One;
              for (int j = 0; j < powerDerivOrder; j++)
                if (i != j)
                  {
                    t12 = datatree.AddMinus(arg2, datatree.AddPossiblyNegativeConstant(j));
                    t11 = datatree.AddTimes(t11, t12);
                  }
              t13 = datatree.AddPlus(t13, t11);
            }
          t13 = datatree.AddTimes(darg2, t13);
          t14 = datatree.AddTimes(f, t13);
          return datatree.AddPlus(first_part, t14);
        }
    case BinaryOpcode::max:
      t11 = datatree.AddGreater(arg1, arg2);
      t12 = datatree.AddTimes(t11, darg1);
      t13 = datatree.AddMinus(datatree.One, t11);
      t14 = datatree.AddTimes(t13, darg2);
      return datatree.AddPlus(t14, t12);
    case BinaryOpcode::min:
      t11 = datatree.AddGreater(arg2, arg1);
      t12 = datatree.AddTimes(t11, darg1);
      t13 = datatree.AddMinus(datatree.One, t11);
      t14 = datatree.AddTimes(t13, darg2);
      return datatree.AddPlus(t14, t12);
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

expr_t
BinaryOpNode::unpackPowerDeriv() const
{
  if (op_code != BinaryOpcode::powerDeriv)
    return const_cast<BinaryOpNode *>(this);

  expr_t front = datatree.One;
  for (int i = 0; i < powerDerivOrder; i++)
    front = datatree.AddTimes(front,
                              datatree.AddMinus(arg2,
                                                datatree.AddPossiblyNegativeConstant(i)));
  expr_t tmp = datatree.AddPower(arg1,
                                 datatree.AddMinus(arg2,
                                                   datatree.AddPossiblyNegativeConstant(powerDerivOrder)));
  return datatree.AddTimes(front, tmp);
}

expr_t
BinaryOpNode::computeDerivative(int deriv_id)
{
  expr_t darg1 = arg1->getDerivative(deriv_id);
  expr_t darg2 = arg2->getDerivative(deriv_id);
  return composeDerivatives(darg1, darg2);
}

int
BinaryOpNode::precedence(ExprNodeOutputType output_type, const temporary_terms_t &temporary_terms) const
{
  // A temporary term behaves as a variable
  if (temporary_terms.find(const_cast<BinaryOpNode *>(this)) != temporary_terms.end())
    return 100;

  switch (op_code)
    {
    case BinaryOpcode::equal:
      return 0;
    case BinaryOpcode::equalEqual:
    case BinaryOpcode::different:
      return 1;
    case BinaryOpcode::lessEqual:
    case BinaryOpcode::greaterEqual:
    case BinaryOpcode::less:
    case BinaryOpcode::greater:
      return 2;
    case BinaryOpcode::plus:
    case BinaryOpcode::minus:
      return 3;
    case BinaryOpcode::times:
    case BinaryOpcode::divide:
      return 4;
    case BinaryOpcode::power:
    case BinaryOpcode::powerDeriv:
      if (isCOutput(output_type))
        // In C, power operator is of the form pow(a, b)
        return 100;
      else
        return 5;
    case BinaryOpcode::min:
    case BinaryOpcode::max:
      return 100;
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

int
BinaryOpNode::precedenceJson(const temporary_terms_t &temporary_terms) const
{
  // A temporary term behaves as a variable
  if (temporary_terms.find(const_cast<BinaryOpNode *>(this)) != temporary_terms.end())
    return 100;

  switch (op_code)
    {
    case BinaryOpcode::equal:
      return 0;
    case BinaryOpcode::equalEqual:
    case BinaryOpcode::different:
      return 1;
    case BinaryOpcode::lessEqual:
    case BinaryOpcode::greaterEqual:
    case BinaryOpcode::less:
    case BinaryOpcode::greater:
      return 2;
    case BinaryOpcode::plus:
    case BinaryOpcode::minus:
      return 3;
    case BinaryOpcode::times:
    case BinaryOpcode::divide:
      return 4;
    case BinaryOpcode::power:
    case BinaryOpcode::powerDeriv:
      return 5;
    case BinaryOpcode::min:
    case BinaryOpcode::max:
      return 100;
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

int
BinaryOpNode::cost(const map<pair<int, int>, temporary_terms_t> &temp_terms_map, bool is_matlab) const
{
  // For a temporary term, the cost is null
  for (const auto &it : temp_terms_map)
    if (it.second.find(const_cast<BinaryOpNode *>(this)) != it.second.end())
      return 0;

  int arg_cost = arg1->cost(temp_terms_map, is_matlab) + arg2->cost(temp_terms_map, is_matlab);

  return cost(arg_cost, is_matlab);
}

int
BinaryOpNode::cost(const vector<vector<temporary_terms_t>> &blocks_temporary_terms, bool is_matlab) const
{
  // For a temporary term, the cost is null
  for (const auto &blk_tt : blocks_temporary_terms)
    for (const auto &eq_tt : blk_tt)
      if (eq_tt.find(const_cast<BinaryOpNode *>(this)) != eq_tt.end())
        return 0;

  int arg_cost = arg1->cost(blocks_temporary_terms, is_matlab) + arg2->cost(blocks_temporary_terms, is_matlab);

  return cost(arg_cost, is_matlab);
}

int
BinaryOpNode::cost(int cost, bool is_matlab) const
{
  if (is_matlab)
    // Cost for Matlab files
    switch (op_code)
      {
      case BinaryOpcode::less:
      case BinaryOpcode::greater:
      case BinaryOpcode::lessEqual:
      case BinaryOpcode::greaterEqual:
      case BinaryOpcode::equalEqual:
      case BinaryOpcode::different:
        return cost + 60;
      case BinaryOpcode::plus:
      case BinaryOpcode::minus:
      case BinaryOpcode::times:
        return cost + 90;
      case BinaryOpcode::max:
      case BinaryOpcode::min:
        return cost + 110;
      case BinaryOpcode::divide:
        return cost + 990;
      case BinaryOpcode::power:
      case BinaryOpcode::powerDeriv:
        return cost + (min_cost_matlab/2+1);
      case BinaryOpcode::equal:
        return cost;
      }
  else
    // Cost for C files
    switch (op_code)
      {
      case BinaryOpcode::less:
      case BinaryOpcode::greater:
      case BinaryOpcode::lessEqual:
      case BinaryOpcode::greaterEqual:
      case BinaryOpcode::equalEqual:
      case BinaryOpcode::different:
        return cost + 2;
      case BinaryOpcode::plus:
      case BinaryOpcode::minus:
      case BinaryOpcode::times:
        return cost + 4;
      case BinaryOpcode::max:
      case BinaryOpcode::min:
        return cost + 5;
      case BinaryOpcode::divide:
        return cost + 15;
      case BinaryOpcode::power:
        return cost + 520;
      case BinaryOpcode::powerDeriv:
        return cost + (min_cost_c/2+1);
      case BinaryOpcode::equal:
        return cost;
      }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

void
BinaryOpNode::computeTemporaryTerms(const pair<int, int> &derivOrder,
                                    map<pair<int, int>, temporary_terms_t> &temp_terms_map,
                                    map<expr_t, pair<int, pair<int, int>>> &reference_count,
                                    bool is_matlab) const
{
  expr_t this2 = const_cast<BinaryOpNode *>(this);
  if (auto it = reference_count.find(this2);
      it == reference_count.end())
    {
      // If this node has never been encountered, set its ref count to one,
      //  and travel through its children
      reference_count[this2] = { 1, derivOrder };
      arg1->computeTemporaryTerms(derivOrder, temp_terms_map, reference_count, is_matlab);
      arg2->computeTemporaryTerms(derivOrder, temp_terms_map, reference_count, is_matlab);
    }
  else
    {
      /* If the node has already been encountered, increment its ref count
         and declare it as a temporary term if it is too costly (except if it is
         an equal node: we don't want them as temporary terms) */
      auto &[nref, min_order] = it->second;
      nref++;
      if (nref * cost(temp_terms_map, is_matlab) > min_cost(is_matlab)
          && op_code != BinaryOpcode::equal)
        temp_terms_map[min_order].insert(this2);
    }
}

void
BinaryOpNode::computeBlockTemporaryTerms(int blk, int eq, vector<vector<temporary_terms_t>> &blocks_temporary_terms,
                                         map<expr_t, tuple<int, int, int>> &reference_count) const
{
  expr_t this2 = const_cast<BinaryOpNode *>(this);
  if (auto it = reference_count.find(this2);
      it == reference_count.end())
    {
      reference_count[this2] = { 1, blk, eq };
      arg1->computeBlockTemporaryTerms(blk, eq, blocks_temporary_terms, reference_count);
      arg2->computeBlockTemporaryTerms(blk, eq, blocks_temporary_terms, reference_count);
    }
  else
    {
      auto &[nref, first_blk, first_eq] = it->second;
      nref++;
      if (nref * cost(blocks_temporary_terms, false) > min_cost_c
          && op_code != BinaryOpcode::equal)
        blocks_temporary_terms[first_blk][first_eq].insert(this2);
    }
}

double
BinaryOpNode::eval_opcode(double v1, BinaryOpcode op_code, double v2, int derivOrder) noexcept(false)
{
  switch (op_code)
    {
    case BinaryOpcode::plus:
      return v1 + v2;
    case BinaryOpcode::minus:
      return v1 - v2;
    case BinaryOpcode::times:
      return v1 * v2;
    case BinaryOpcode::divide:
      return v1 / v2;
    case BinaryOpcode::power:
      return pow(v1, v2);
    case BinaryOpcode::powerDeriv:
      if (fabs(v1) < near_zero && v2 > 0
          && derivOrder > v2
          && fabs(v2-nearbyint(v2)) < near_zero)
        return 0.0;
      else
        {
          double dxp = pow(v1, v2-derivOrder);
          for (int i = 0; i < derivOrder; i++)
            dxp *= v2--;
          return dxp;
        }
    case BinaryOpcode::max:
      if (v1 < v2)
        return v2;
      else
        return v1;
    case BinaryOpcode::min:
      if (v1 > v2)
        return v2;
      else
        return v1;
    case BinaryOpcode::less:
      return v1 < v2;
    case BinaryOpcode::greater:
      return v1 > v2;
    case BinaryOpcode::lessEqual:
      return v1 <= v2;
    case BinaryOpcode::greaterEqual:
      return v1 >= v2;
    case BinaryOpcode::equalEqual:
      return v1 == v2;
    case BinaryOpcode::different:
      return v1 != v2;
    case BinaryOpcode::equal:
      throw EvalException();
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

double
BinaryOpNode::eval(const eval_context_t &eval_context) const noexcept(false)
{
  double v1 = arg1->eval(eval_context);
  double v2 = arg2->eval(eval_context);
  return eval_opcode(v1, op_code, v2, powerDerivOrder);
}

void
BinaryOpNode::compile(ostream &CompileCode, unsigned int &instruction_number,
                      bool lhs_rhs, const temporary_terms_t &temporary_terms,
                      const temporary_terms_idxs_t &temporary_terms_idxs, bool dynamic, bool steady_dynamic,
                      const deriv_node_temp_terms_t &tef_terms) const
{
  // If current node is a temporary term
  if (auto this2 = const_cast<BinaryOpNode *>(this);
      temporary_terms.find(this2) != temporary_terms.end())
    {
      if (dynamic)
        {
          FLDT_ fldt(temporary_terms_idxs.at(this2));
          fldt.write(CompileCode, instruction_number);
        }
      else
        {
          FLDST_ fldst(temporary_terms_idxs.at(this2));
          fldst.write(CompileCode, instruction_number);
        }
      return;
    }
  if (op_code == BinaryOpcode::powerDeriv)
    {
      FLDC_ fldc(powerDerivOrder);
      fldc.write(CompileCode, instruction_number);
    }
  arg1->compile(CompileCode, instruction_number, lhs_rhs, temporary_terms, temporary_terms_idxs, dynamic, steady_dynamic, tef_terms);
  arg2->compile(CompileCode, instruction_number, lhs_rhs, temporary_terms, temporary_terms_idxs, dynamic, steady_dynamic, tef_terms);
  FBINARY_ fbinary{static_cast<int>(op_code)};
  fbinary.write(CompileCode, instruction_number);
}

bool
BinaryOpNode::containsExternalFunction() const
{
  return arg1->containsExternalFunction()
    || arg2->containsExternalFunction();
}

void
BinaryOpNode::writeJsonAST(ostream &output) const
{
  output << R"({"node_type" : "BinaryOpNode",)"
         << R"( "op" : ")";
  switch (op_code)
    {
    case BinaryOpcode::plus:
      output << "+";
      break;
    case BinaryOpcode::minus:
      output << "-";
      break;
    case BinaryOpcode::times:
      output << "*";
      break;
    case BinaryOpcode::divide:
      output << "/";
      break;
    case BinaryOpcode::power:
      output << "^";
      break;
    case BinaryOpcode::less:
      output << "<";
      break;
    case BinaryOpcode::greater:
      output << ">";
      break;
    case BinaryOpcode::lessEqual:
      output << "<=";
      break;
    case BinaryOpcode::greaterEqual:
      output << ">=";
      break;
    case BinaryOpcode::equalEqual:
      output << "==";
      break;
    case BinaryOpcode::different:
      output << "!=";
      break;
    case BinaryOpcode::equal:
      output << "=";
      break;
    case BinaryOpcode::max:
      output << "max";
      break;
    case BinaryOpcode::min:
      output << "min";
      break;
    case BinaryOpcode::powerDeriv:
      output << "power_deriv";
      break;
    }
  output << R"(", "arg1" : )";
  arg1->writeJsonAST(output);
  output << R"(, "arg2" : )";
  arg2->writeJsonAST(output);
  output << "}";
}

void
BinaryOpNode::writeJsonOutput(ostream &output,
                              const temporary_terms_t &temporary_terms,
                              const deriv_node_temp_terms_t &tef_terms,
                              bool isdynamic) const
{
  // If current node is a temporary term
  if (temporary_terms.find(const_cast<BinaryOpNode *>(this)) != temporary_terms.end())
    {
      output << "T" << idx;
      return;
    }

  if (op_code == BinaryOpcode::max || op_code == BinaryOpcode::min)
    {
      switch (op_code)
        {
        case BinaryOpcode::max:
          output << "max(";
          break;
        case BinaryOpcode::min:
          output << "min(";
          break;
        default:
          ;
        }
      arg1->writeJsonOutput(output, temporary_terms, tef_terms, isdynamic);
      output << ",";
      arg2->writeJsonOutput(output, temporary_terms, tef_terms, isdynamic);
      output << ")";
      return;
    }

  if (op_code == BinaryOpcode::powerDeriv)
    {
      output << "get_power_deriv(";
      arg1->writeJsonOutput(output, temporary_terms, tef_terms, isdynamic);
      output << ",";
      arg2->writeJsonOutput(output, temporary_terms, tef_terms, isdynamic);
      output << "," << powerDerivOrder << ")";
      return;
    }

  int prec = precedenceJson(temporary_terms);

  bool close_parenthesis = false;

  // If left argument has a lower precedence, or if current and left argument are both power operators,
  // add parenthesis around left argument
  if (auto barg1 = dynamic_cast<BinaryOpNode *>(arg1);
      arg1->precedenceJson(temporary_terms) < prec
      || (op_code == BinaryOpcode::power && barg1 && barg1->op_code == BinaryOpcode::power))
    {
      output << "(";
      close_parenthesis = true;
    }

  // Write left argument
  arg1->writeJsonOutput(output, temporary_terms, tef_terms, isdynamic);

  if (close_parenthesis)
    output << ")";

  // Write current operator symbol
  switch (op_code)
    {
    case BinaryOpcode::plus:
      output << "+";
      break;
    case BinaryOpcode::minus:
      output << "-";
      break;
    case BinaryOpcode::times:
      output << "*";
      break;
    case BinaryOpcode::divide:
      output << "/";
      break;
    case BinaryOpcode::power:
      output << "^";
      break;
    case BinaryOpcode::less:
      output << "<";
      break;
    case BinaryOpcode::greater:
      output << ">";
      break;
    case BinaryOpcode::lessEqual:
      output << "<=";
      break;
    case BinaryOpcode::greaterEqual:
      output << ">=";
      break;
    case BinaryOpcode::equalEqual:
      output << "==";
      break;
    case BinaryOpcode::different:
      output << "!=";
      break;
    case BinaryOpcode::equal:
      output << "=";
      break;
    default:
      ;
    }

  close_parenthesis = false;

  /* Add parenthesis around right argument if:
     - its precedence is lower than those of the current node
     - it is a power operator and current operator is also a power operator
     - it is a minus operator with same precedence than current operator
     - it is a divide operator with same precedence than current operator */
  auto barg2 = dynamic_cast<BinaryOpNode *>(arg2);
  if (int arg2_prec = arg2->precedenceJson(temporary_terms); arg2_prec < prec
      || (op_code == BinaryOpcode::power && barg2 && barg2->op_code == BinaryOpcode::power)
      || (op_code == BinaryOpcode::minus && arg2_prec == prec)
      || (op_code == BinaryOpcode::divide && arg2_prec == prec))
    {
      output << "(";
      close_parenthesis = true;
    }

  // Write right argument
  arg2->writeJsonOutput(output, temporary_terms, tef_terms, isdynamic);

  if (close_parenthesis)
    output << ")";
}

void
BinaryOpNode::writeOutput(ostream &output, ExprNodeOutputType output_type,
                          const temporary_terms_t &temporary_terms,
                          const temporary_terms_idxs_t &temporary_terms_idxs,
                          const deriv_node_temp_terms_t &tef_terms) const
{
  if (checkIfTemporaryTermThenWrite(output, output_type, temporary_terms, temporary_terms_idxs))
    return;

  // Treat derivative of Power
  if (op_code == BinaryOpcode::powerDeriv)
    {
      if (isLatexOutput(output_type))
        unpackPowerDeriv()->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
      else
        {
          if (output_type == ExprNodeOutputType::juliaStaticModel || output_type == ExprNodeOutputType::juliaDynamicModel)
            output << "get_power_deriv(";
          else
            output << "getPowerDeriv(";
          arg1->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
          output << ",";
          arg2->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
          output << "," << powerDerivOrder << ")";
        }
      return;
    }

  // Treat special case of power operator in C, and case of max and min operators
  if ((op_code == BinaryOpcode::power && isCOutput(output_type)) || op_code == BinaryOpcode::max || op_code == BinaryOpcode::min)
    {
      switch (op_code)
        {
        case BinaryOpcode::power:
          output << "pow(";
          break;
        case BinaryOpcode::max:
          if (isCOutput(output_type))
            output << "fmax(";
          else
            output << "max(";
          break;
        case BinaryOpcode::min:
          if (isCOutput(output_type))
            output << "fmin(";
          else
            output << "min(";
          break;
        default:
          ;
        }
      arg1->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
      output << ",";
      arg2->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
      output << ")";
      return;
    }

  int prec = precedence(output_type, temporary_terms);

  bool close_parenthesis = false;

  if (isLatexOutput(output_type) && op_code == BinaryOpcode::divide)
    output << R"(\frac{)";
  else
    {
      // If left argument has a lower precedence, or if current and left argument are both power operators, add parenthesis around left argument
      auto barg1 = dynamic_cast<BinaryOpNode *>(arg1);
      if (arg1->precedence(output_type, temporary_terms) < prec
          || (op_code == BinaryOpcode::power && barg1 != nullptr && barg1->op_code == BinaryOpcode::power))
        {
          output << LEFT_PAR(output_type);
          close_parenthesis = true;
        }
    }

  // Write left argument
  arg1->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);

  if (close_parenthesis)
    output << RIGHT_PAR(output_type);

  if (isLatexOutput(output_type) && op_code == BinaryOpcode::divide)
    output << "}";

  // Write current operator symbol
  switch (op_code)
    {
    case BinaryOpcode::plus:
      output << "+";
      break;
    case BinaryOpcode::minus:
      output << "-";
      break;
    case BinaryOpcode::times:
      if (isLatexOutput(output_type))
        output << R"(\, )";
      else if (output_type == ExprNodeOutputType::occbinDifferenceFile)
        output << ".*"; // This file operates on vectors, see dynare#1826
      else
        output << "*";
      break;
    case BinaryOpcode::divide:
      if (!isLatexOutput(output_type))
        {
          if (output_type == ExprNodeOutputType::occbinDifferenceFile)
            output << "./"; // This file operates on vectors, see dynare#1826
          else
            output << "/";
        }
      break;
    case BinaryOpcode::power:
      if (output_type == ExprNodeOutputType::occbinDifferenceFile)
        output << ".^"; // This file operates on vectors, see dynare#1826
      else
        output << "^";
      break;
    case BinaryOpcode::less:
      output << "<";
      break;
    case BinaryOpcode::greater:
      output << ">";
      break;
    case BinaryOpcode::lessEqual:
      if (isLatexOutput(output_type))
        output << R"(\leq )";
      else
        output << "<=";
      break;
    case BinaryOpcode::greaterEqual:
      if (isLatexOutput(output_type))
        output << R"(\geq )";
      else
        output << ">=";
      break;
    case BinaryOpcode::equalEqual:
      output << "==";
      break;
    case BinaryOpcode::different:
      if (isMatlabOutput(output_type))
        output << "~=";
      else
        {
          if (isCOutput(output_type) || isJuliaOutput(output_type))
            output << "!=";
          else
            output << R"(\neq )";
        }
      break;
    case BinaryOpcode::equal:
      output << "=";
      break;
    default:
      ;
    }

  close_parenthesis = false;

  if (isLatexOutput(output_type) && (op_code == BinaryOpcode::power || op_code == BinaryOpcode::divide))
    output << "{";
  else
    {
      /* Add parenthesis around right argument if:
         - its precedence is lower than those of the current node
         - it is a power operator and current operator is also a power operator
         - it is a minus operator with same precedence than current operator
         - it is a divide operator with same precedence than current operator */
      auto barg2 = dynamic_cast<BinaryOpNode *>(arg2);
      if (int arg2_prec = arg2->precedence(output_type, temporary_terms); arg2_prec < prec
          || (op_code == BinaryOpcode::power && barg2 && barg2->op_code == BinaryOpcode::power && !isLatexOutput(output_type))
          || (op_code == BinaryOpcode::minus && arg2_prec == prec)
          || (op_code == BinaryOpcode::divide && arg2_prec == prec && !isLatexOutput(output_type)))
        {
          output << LEFT_PAR(output_type);
          close_parenthesis = true;
        }
    }

  // Write right argument
  arg2->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);

  if (isLatexOutput(output_type) && (op_code == BinaryOpcode::power || op_code == BinaryOpcode::divide))
    output << "}";

  if (close_parenthesis)
    output << RIGHT_PAR(output_type);
}

void
BinaryOpNode::writeExternalFunctionOutput(ostream &output, ExprNodeOutputType output_type,
                                          const temporary_terms_t &temporary_terms,
                                          const temporary_terms_idxs_t &temporary_terms_idxs,
                                          deriv_node_temp_terms_t &tef_terms) const
{
  arg1->writeExternalFunctionOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
  arg2->writeExternalFunctionOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
}

void
BinaryOpNode::writeJsonExternalFunctionOutput(vector<string> &efout,
                                              const temporary_terms_t &temporary_terms,
                                              deriv_node_temp_terms_t &tef_terms,
                                              bool isdynamic) const
{
  arg1->writeJsonExternalFunctionOutput(efout, temporary_terms, tef_terms, isdynamic);
  arg2->writeJsonExternalFunctionOutput(efout, temporary_terms, tef_terms, isdynamic);
}

void
BinaryOpNode::compileExternalFunctionOutput(ostream &CompileCode, unsigned int &instruction_number,
                                            bool lhs_rhs, const temporary_terms_t &temporary_terms,
                                            const temporary_terms_idxs_t &temporary_terms_idxs, bool dynamic, bool steady_dynamic,
                                            deriv_node_temp_terms_t &tef_terms) const
{
  arg1->compileExternalFunctionOutput(CompileCode, instruction_number, lhs_rhs, temporary_terms,
                                      temporary_terms_idxs, dynamic, steady_dynamic, tef_terms);
  arg2->compileExternalFunctionOutput(CompileCode, instruction_number, lhs_rhs, temporary_terms,
                                      temporary_terms_idxs, dynamic, steady_dynamic, tef_terms);
}

int
BinaryOpNode::VarMaxLag(const set<expr_t> &lhs_lag_equiv) const
{
  return max(arg1->VarMaxLag(lhs_lag_equiv),
             arg2->VarMaxLag(lhs_lag_equiv));
}

void
BinaryOpNode::collectVARLHSVariable(set<expr_t> &result) const
{
  cerr << "ERROR: you can only have variables or unary ops on LHS of VAR" << endl;
  exit(EXIT_FAILURE);
}

void
BinaryOpNode::collectDynamicVariables(SymbolType type_arg, set<pair<int, int>> &result) const
{
  arg1->collectDynamicVariables(type_arg, result);
  arg2->collectDynamicVariables(type_arg, result);
}

expr_t
BinaryOpNode::Compute_RHS(expr_t arg1, expr_t arg2, int op, int op_type) const
{
  switch (op_type)
    {
    case 0: /*Unary Operator*/
      switch (static_cast<UnaryOpcode>(op))
        {
        case UnaryOpcode::uminus:
          return datatree.AddUMinus(arg1);
        case UnaryOpcode::exp:
          return datatree.AddExp(arg1);
        case UnaryOpcode::log:
          return datatree.AddLog(arg1);
        case UnaryOpcode::log10:
          return datatree.AddLog10(arg1);
        default:
          cerr << "BinaryOpNode::Compute_RHS: case not handled";
          exit(EXIT_FAILURE);
        }
      break;
    case 1: /*Binary Operator*/
      switch (static_cast<BinaryOpcode>(op))
        {
        case BinaryOpcode::plus:
          return datatree.AddPlus(arg1, arg2);
        case BinaryOpcode::minus:
          return datatree.AddMinus(arg1, arg2);
        case BinaryOpcode::times:
          return datatree.AddTimes(arg1, arg2);
        case BinaryOpcode::divide:
          return datatree.AddDivide(arg1, arg2);
        case BinaryOpcode::power:
          return datatree.AddPower(arg1, arg2);
        default:
          cerr << "BinaryOpNode::Compute_RHS: case not handled";
          exit(EXIT_FAILURE);
        }
      break;
    }
  return nullptr;
}

void
BinaryOpNode::computeSubExprContainingVariable(int symb_id, int lag, set<expr_t> &contain_var) const
{
  arg1->computeSubExprContainingVariable(symb_id, lag, contain_var);
  arg2->computeSubExprContainingVariable(symb_id, lag, contain_var);
  if (contain_var.count(arg1) > 0 || contain_var.count(arg2) > 0)
    contain_var.insert(const_cast<BinaryOpNode *>(this));
}

BinaryOpNode *
BinaryOpNode::normalizeEquationHelper(const set<expr_t> &contain_var, expr_t rhs) const
{
  assert(contain_var.count(const_cast<BinaryOpNode *>(this)) > 0);

  bool arg1_contains_var = contain_var.count(arg1) > 0;
  bool arg2_contains_var = contain_var.count(arg2) > 0;
  assert(arg1_contains_var || arg2_contains_var);

  if (arg1_contains_var && arg2_contains_var)
    throw NormalizationFailed();

  switch (op_code)
    {
    case BinaryOpcode::plus:
      if (arg1_contains_var)
        rhs = datatree.AddMinus(rhs, arg2);
      else
        rhs = datatree.AddMinus(rhs, arg1);
      break;
    case BinaryOpcode::minus:
      if (arg1_contains_var)
        rhs = datatree.AddPlus(rhs, arg2);
      else
        rhs = datatree.AddMinus(arg1, rhs);
      break;
    case BinaryOpcode::times:
      if (arg1_contains_var)
        rhs = datatree.AddDivide(rhs, arg2);
      else
        rhs = datatree.AddDivide(rhs, arg1);
      break;
    case BinaryOpcode::divide:
      if (arg1_contains_var)
        rhs = datatree.AddTimes(rhs, arg2);
      else
        rhs = datatree.AddDivide(arg1, rhs);
      break;
    case BinaryOpcode::power:
      if (arg1_contains_var)
        rhs = datatree.AddPower(rhs, datatree.AddDivide(datatree.One, arg2));
      else
        // a^f(X)=rhs is normalized in f(X)=ln(rhs)/ln(a)
        rhs = datatree.AddDivide(datatree.AddLog(rhs), datatree.AddLog(arg1));
      break;
    case BinaryOpcode::equal:
      cerr << "BinaryOpCode::normalizeEquationHelper: this case should not happen" << endl;
      exit(EXIT_FAILURE);
    default:
      throw NormalizationFailed();
    }

  if (arg1_contains_var)
    return arg1->normalizeEquationHelper(contain_var, rhs);
  else
    return arg2->normalizeEquationHelper(contain_var, rhs);
}

BinaryOpNode *
BinaryOpNode::normalizeEquation(int symb_id, int lag) const
{
  assert(op_code == BinaryOpcode::equal);

  set<expr_t> contain_var;
  computeSubExprContainingVariable(symb_id, lag, contain_var);

  bool arg1_contains_var = contain_var.count(arg1) > 0;
  bool arg2_contains_var = contain_var.count(arg2) > 0;
  assert(arg1_contains_var || arg2_contains_var);

  if (arg1_contains_var && arg2_contains_var)
    throw NormalizationFailed();

  return arg1_contains_var ? arg1->normalizeEquationHelper(contain_var, arg2)
    : arg2->normalizeEquationHelper(contain_var, arg1);
}

expr_t
BinaryOpNode::getChainRuleDerivative(int deriv_id, const map<int, BinaryOpNode *> &recursive_variables)
{
  expr_t darg1 = arg1->getChainRuleDerivative(deriv_id, recursive_variables);
  expr_t darg2 = arg2->getChainRuleDerivative(deriv_id, recursive_variables);
  return composeDerivatives(darg1, darg2);
}

expr_t
BinaryOpNode::buildSimilarBinaryOpNode(expr_t alt_arg1, expr_t alt_arg2, DataTree &alt_datatree) const
{
  switch (op_code)
    {
    case BinaryOpcode::plus:
      return alt_datatree.AddPlus(alt_arg1, alt_arg2);
    case BinaryOpcode::minus:
      return alt_datatree.AddMinus(alt_arg1, alt_arg2);
    case BinaryOpcode::times:
      return alt_datatree.AddTimes(alt_arg1, alt_arg2);
    case BinaryOpcode::divide:
      return alt_datatree.AddDivide(alt_arg1, alt_arg2);
    case BinaryOpcode::power:
      return alt_datatree.AddPower(alt_arg1, alt_arg2);
    case BinaryOpcode::equal:
      return alt_datatree.AddEqual(alt_arg1, alt_arg2);
    case BinaryOpcode::max:
      return alt_datatree.AddMax(alt_arg1, alt_arg2);
    case BinaryOpcode::min:
      return alt_datatree.AddMin(alt_arg1, alt_arg2);
    case BinaryOpcode::less:
      return alt_datatree.AddLess(alt_arg1, alt_arg2);
    case BinaryOpcode::greater:
      return alt_datatree.AddGreater(alt_arg1, alt_arg2);
    case BinaryOpcode::lessEqual:
      return alt_datatree.AddLessEqual(alt_arg1, alt_arg2);
    case BinaryOpcode::greaterEqual:
      return alt_datatree.AddGreaterEqual(alt_arg1, alt_arg2);
    case BinaryOpcode::equalEqual:
      return alt_datatree.AddEqualEqual(alt_arg1, alt_arg2);
    case BinaryOpcode::different:
      return alt_datatree.AddDifferent(alt_arg1, alt_arg2);
    case BinaryOpcode::powerDeriv:
      return alt_datatree.AddPowerDeriv(alt_arg1, alt_arg2, powerDerivOrder);
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

expr_t
BinaryOpNode::toStatic(DataTree &static_datatree) const
{
  expr_t sarg1 = arg1->toStatic(static_datatree);
  expr_t sarg2 = arg2->toStatic(static_datatree);
  return buildSimilarBinaryOpNode(sarg1, sarg2, static_datatree);
}

void
BinaryOpNode::computeXrefs(EquationInfo &ei) const
{
  arg1->computeXrefs(ei);
  arg2->computeXrefs(ei);
}

expr_t
BinaryOpNode::clone(DataTree &datatree) const
{
  expr_t substarg1 = arg1->clone(datatree);
  expr_t substarg2 = arg2->clone(datatree);
  return buildSimilarBinaryOpNode(substarg1, substarg2, datatree);
}

int
BinaryOpNode::maxEndoLead() const
{
  return max(arg1->maxEndoLead(), arg2->maxEndoLead());
}

int
BinaryOpNode::maxExoLead() const
{
  return max(arg1->maxExoLead(), arg2->maxExoLead());
}

int
BinaryOpNode::maxEndoLag() const
{
  return max(arg1->maxEndoLag(), arg2->maxEndoLag());
}

int
BinaryOpNode::maxExoLag() const
{
  return max(arg1->maxExoLag(), arg2->maxExoLag());
}

int
BinaryOpNode::maxLead() const
{
  return max(arg1->maxLead(), arg2->maxLead());
}

int
BinaryOpNode::maxLag() const
{
  return max(arg1->maxLag(), arg2->maxLag());
}

int
BinaryOpNode::maxLagWithDiffsExpanded() const
{
  return max(arg1->maxLagWithDiffsExpanded(), arg2->maxLagWithDiffsExpanded());
}

expr_t
BinaryOpNode::undiff() const
{
  expr_t arg1subst = arg1->undiff();
  expr_t arg2subst = arg2->undiff();
  return buildSimilarBinaryOpNode(arg1subst, arg2subst, datatree);
}

expr_t
BinaryOpNode::decreaseLeadsLags(int n) const
{
  expr_t arg1subst = arg1->decreaseLeadsLags(n);
  expr_t arg2subst = arg2->decreaseLeadsLags(n);
  return buildSimilarBinaryOpNode(arg1subst, arg2subst, datatree);
}

expr_t
BinaryOpNode::decreaseLeadsLagsPredeterminedVariables() const
{
  expr_t arg1subst = arg1->decreaseLeadsLagsPredeterminedVariables();
  expr_t arg2subst = arg2->decreaseLeadsLagsPredeterminedVariables();
  return buildSimilarBinaryOpNode(arg1subst, arg2subst, datatree);
}

expr_t
BinaryOpNode::substituteEndoLeadGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const
{
  expr_t arg1subst, arg2subst;
  int maxendolead1 = arg1->maxEndoLead(), maxendolead2 = arg2->maxEndoLead();

  if (maxendolead1 < 2 && maxendolead2 < 2)
    return const_cast<BinaryOpNode *>(this);
  if (deterministic_model)
    {
      arg1subst = maxendolead1 >= 2 ? arg1->substituteEndoLeadGreaterThanTwo(subst_table, neweqs, deterministic_model) : arg1;
      arg2subst = maxendolead2 >= 2 ? arg2->substituteEndoLeadGreaterThanTwo(subst_table, neweqs, deterministic_model) : arg2;
      return buildSimilarBinaryOpNode(arg1subst, arg2subst, datatree);
    }
  else
    switch (op_code)
      {
      case BinaryOpcode::plus:
      case BinaryOpcode::minus:
      case BinaryOpcode::equal:
        arg1subst = maxendolead1 >= 2 ? arg1->substituteEndoLeadGreaterThanTwo(subst_table, neweqs, deterministic_model) : arg1;
        arg2subst = maxendolead2 >= 2 ? arg2->substituteEndoLeadGreaterThanTwo(subst_table, neweqs, deterministic_model) : arg2;
        return buildSimilarBinaryOpNode(arg1subst, arg2subst, datatree);
      case BinaryOpcode::times:
      case BinaryOpcode::divide:
        if (maxendolead1 >= 2 && maxendolead2 == 0 && arg2->maxExoLead() == 0)
          {
            arg1subst = arg1->substituteEndoLeadGreaterThanTwo(subst_table, neweqs, deterministic_model);
            return buildSimilarBinaryOpNode(arg1subst, arg2, datatree);
          }
        if (maxendolead1 == 0 && arg1->maxExoLead() == 0
            && maxendolead2 >= 2 && op_code == BinaryOpcode::times)
          {
            arg2subst = arg2->substituteEndoLeadGreaterThanTwo(subst_table, neweqs, deterministic_model);
            return buildSimilarBinaryOpNode(arg1, arg2subst, datatree);
          }
        return createEndoLeadAuxiliaryVarForMyself(subst_table, neweqs);
      default:
        return createEndoLeadAuxiliaryVarForMyself(subst_table, neweqs);
      }
}

expr_t
BinaryOpNode::substituteEndoLagGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  expr_t arg1subst = arg1->substituteEndoLagGreaterThanTwo(subst_table, neweqs);
  expr_t arg2subst = arg2->substituteEndoLagGreaterThanTwo(subst_table, neweqs);
  return buildSimilarBinaryOpNode(arg1subst, arg2subst, datatree);
}

expr_t
BinaryOpNode::substituteExoLead(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const
{
  expr_t arg1subst, arg2subst;
  int maxexolead1 = arg1->maxExoLead(), maxexolead2 = arg2->maxExoLead();

  if (maxexolead1 < 1 && maxexolead2 < 1)
    return const_cast<BinaryOpNode *>(this);
  if (deterministic_model)
    {
      arg1subst = maxexolead1 >= 1 ? arg1->substituteExoLead(subst_table, neweqs, deterministic_model) : arg1;
      arg2subst = maxexolead2 >= 1 ? arg2->substituteExoLead(subst_table, neweqs, deterministic_model) : arg2;
      return buildSimilarBinaryOpNode(arg1subst, arg2subst, datatree);
    }
  else
    switch (op_code)
      {
      case BinaryOpcode::plus:
      case BinaryOpcode::minus:
      case BinaryOpcode::equal:
        arg1subst = maxexolead1 >= 1 ? arg1->substituteExoLead(subst_table, neweqs, deterministic_model) : arg1;
        arg2subst = maxexolead2 >= 1 ? arg2->substituteExoLead(subst_table, neweqs, deterministic_model) : arg2;
        return buildSimilarBinaryOpNode(arg1subst, arg2subst, datatree);
      case BinaryOpcode::times:
      case BinaryOpcode::divide:
        if (maxexolead1 >= 1 && maxexolead2 == 0 && arg2->maxEndoLead() == 0)
          {
            arg1subst = arg1->substituteExoLead(subst_table, neweqs, deterministic_model);
            return buildSimilarBinaryOpNode(arg1subst, arg2, datatree);
          }
        if (maxexolead1 == 0 && arg1->maxEndoLead() == 0
            && maxexolead2 >= 1 && op_code == BinaryOpcode::times)
          {
            arg2subst = arg2->substituteExoLead(subst_table, neweqs, deterministic_model);
            return buildSimilarBinaryOpNode(arg1, arg2subst, datatree);
          }
        return createExoLeadAuxiliaryVarForMyself(subst_table, neweqs);
      default:
        return createExoLeadAuxiliaryVarForMyself(subst_table, neweqs);
      }
}

expr_t
BinaryOpNode::substituteExoLag(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  expr_t arg1subst = arg1->substituteExoLag(subst_table, neweqs);
  expr_t arg2subst = arg2->substituteExoLag(subst_table, neweqs);
  return buildSimilarBinaryOpNode(arg1subst, arg2subst, datatree);
}

expr_t
BinaryOpNode::substituteExpectation(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool partial_information_model) const
{
  expr_t arg1subst = arg1->substituteExpectation(subst_table, neweqs, partial_information_model);
  expr_t arg2subst = arg2->substituteExpectation(subst_table, neweqs, partial_information_model);
  return buildSimilarBinaryOpNode(arg1subst, arg2subst, datatree);
}

expr_t
BinaryOpNode::substituteAdl() const
{
  expr_t arg1subst = arg1->substituteAdl();
  expr_t arg2subst = arg2->substituteAdl();
  return buildSimilarBinaryOpNode(arg1subst, arg2subst, datatree);
}

expr_t
BinaryOpNode::substituteModelLocalVariables() const
{
  expr_t arg1subst = arg1->substituteModelLocalVariables();
  expr_t arg2subst = arg2->substituteModelLocalVariables();
  return buildSimilarBinaryOpNode(arg1subst, arg2subst, datatree);
}

expr_t
BinaryOpNode::substituteVarExpectation(const map<string, expr_t> &subst_table) const
{
  expr_t arg1subst = arg1->substituteVarExpectation(subst_table);
  expr_t arg2subst = arg2->substituteVarExpectation(subst_table);
  return buildSimilarBinaryOpNode(arg1subst, arg2subst, datatree);
}

void
BinaryOpNode::findUnaryOpNodesForAuxVarCreation(lag_equivalence_table_t &nodes) const
{
  arg1->findUnaryOpNodesForAuxVarCreation(nodes);
  arg2->findUnaryOpNodesForAuxVarCreation(nodes);
}

void
BinaryOpNode::findDiffNodes(lag_equivalence_table_t &nodes) const
{
  arg1->findDiffNodes(nodes);
  arg2->findDiffNodes(nodes);
}

expr_t
BinaryOpNode::substituteDiff(const lag_equivalence_table_t &nodes, subst_table_t &subst_table,
                             vector<BinaryOpNode *> &neweqs) const
{
  expr_t arg1subst = arg1->substituteDiff(nodes, subst_table, neweqs);
  expr_t arg2subst = arg2->substituteDiff(nodes, subst_table, neweqs);
  return buildSimilarBinaryOpNode(arg1subst, arg2subst, datatree);
}

expr_t
BinaryOpNode::substituteUnaryOpNodes(const lag_equivalence_table_t &nodes, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  expr_t arg1subst = arg1->substituteUnaryOpNodes(nodes, subst_table, neweqs);
  expr_t arg2subst = arg2->substituteUnaryOpNodes(nodes, subst_table, neweqs);
  return buildSimilarBinaryOpNode(arg1subst, arg2subst, datatree);
}

int
BinaryOpNode::countDiffs() const
{
  return max(arg1->countDiffs(), arg2->countDiffs());
}

expr_t
BinaryOpNode::substitutePacExpectation(const string &name, expr_t subexpr)
{
  expr_t arg1subst = arg1->substitutePacExpectation(name, subexpr);
  expr_t arg2subst = arg2->substitutePacExpectation(name, subexpr);
  return buildSimilarBinaryOpNode(arg1subst, arg2subst, datatree);
}

expr_t
BinaryOpNode::substitutePacTargetNonstationary(const string &name, expr_t subexpr)
{
  expr_t arg1subst = arg1->substitutePacTargetNonstationary(name, subexpr);
  expr_t arg2subst = arg2->substitutePacTargetNonstationary(name, subexpr);
  return buildSimilarBinaryOpNode(arg1subst, arg2subst, datatree);
}

expr_t
BinaryOpNode::differentiateForwardVars(const vector<string> &subset, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  expr_t arg1subst = arg1->differentiateForwardVars(subset, subst_table, neweqs);
  expr_t arg2subst = arg2->differentiateForwardVars(subset, subst_table, neweqs);
  return buildSimilarBinaryOpNode(arg1subst, arg2subst, datatree);
}

expr_t
BinaryOpNode::addMultipliersToConstraints(int i)
{
  int symb_id = datatree.symbol_table.addMultiplierAuxiliaryVar(i);
  expr_t newAuxLM = datatree.AddVariable(symb_id, 0);
  return datatree.AddEqual(datatree.AddTimes(newAuxLM, datatree.AddMinus(arg1, arg2)), datatree.Zero);
}

bool
BinaryOpNode::isNumConstNodeEqualTo(double value) const
{
  return false;
}

bool
BinaryOpNode::isVariableNodeEqualTo(SymbolType type_arg, int variable_id, int lag_arg) const
{
  return false;
}

bool
BinaryOpNode::containsPacExpectation(const string &pac_model_name) const
{
  return arg1->containsPacExpectation(pac_model_name) || arg2->containsPacExpectation(pac_model_name);
}

bool
BinaryOpNode::containsPacTargetNonstationary(const string &pac_model_name) const
{
  return arg1->containsPacTargetNonstationary(pac_model_name) || arg2->containsPacTargetNonstationary(pac_model_name);
}

expr_t
BinaryOpNode::replaceTrendVar() const
{
  expr_t arg1subst = arg1->replaceTrendVar();
  expr_t arg2subst = arg2->replaceTrendVar();
  return buildSimilarBinaryOpNode(arg1subst, arg2subst, datatree);
}

expr_t
BinaryOpNode::detrend(int symb_id, bool log_trend, expr_t trend) const
{
  expr_t arg1subst = arg1->detrend(symb_id, log_trend, trend);
  expr_t arg2subst = arg2->detrend(symb_id, log_trend, trend);
  return buildSimilarBinaryOpNode(arg1subst, arg2subst, datatree);
}

expr_t
BinaryOpNode::removeTrendLeadLag(const map<int, expr_t> &trend_symbols_map) const
{
  expr_t arg1subst = arg1->removeTrendLeadLag(trend_symbols_map);
  expr_t arg2subst = arg2->removeTrendLeadLag(trend_symbols_map);
  return buildSimilarBinaryOpNode(arg1subst, arg2subst, datatree);
}

bool
BinaryOpNode::isInStaticForm() const
{
  return arg1->isInStaticForm() && arg2->isInStaticForm();
}

bool
BinaryOpNode::findTargetVariableHelper1(int lhs_symb_id, int rhs_symb_id) const
{
  if (lhs_symb_id == rhs_symb_id)
    return true;

  try
    {
      if (datatree.symbol_table.isAuxiliaryVariable(rhs_symb_id)
          && lhs_symb_id == datatree.symbol_table.getOrigSymbIdForAuxVar(rhs_symb_id))
        return true;
    }
  catch (...)
    {
    }
  return false;
}

int
BinaryOpNode::findTargetVariableHelper(const expr_t arg1, const expr_t arg2,
                                       int lhs_symb_id) const
{
  set<int> params;
  arg1->collectVariables(SymbolType::parameter, params);
  if (params.size() != 1)
    return -1;

  set<pair<int, int>> endogs;
  arg2->collectDynamicVariables(SymbolType::endogenous, endogs);
  if (auto testarg2 = dynamic_cast<BinaryOpNode *>(arg2);
      endogs.size() == 2 && testarg2 && testarg2->op_code == BinaryOpcode::minus
      && dynamic_cast<VariableNode *>(testarg2->arg1)
      && dynamic_cast<VariableNode *>(testarg2->arg2))
    {
      if (findTargetVariableHelper1(lhs_symb_id, endogs.begin()->first))
        return endogs.rbegin()->first;
      else if (findTargetVariableHelper1(lhs_symb_id, endogs.rbegin()->first))
        return endogs.begin()->first;
    }
  return -1;
}

int
BinaryOpNode::findTargetVariable(int lhs_symb_id) const
{
  int retval = findTargetVariableHelper(arg1, arg2, lhs_symb_id);
  if (retval < 0)
    retval = findTargetVariableHelper(arg2, arg1, lhs_symb_id);
  if (retval < 0)
    retval = arg1->findTargetVariable(lhs_symb_id);
  if (retval < 0)
    retval = arg2->findTargetVariable(lhs_symb_id);
  return retval;
}

void
BinaryOpNode::getPacAREC(int lhs_symb_id, int lhs_orig_symb_id,
                         pair<int, vector<tuple<int, bool, int>>> &ec_params_and_vars,
                         vector<tuple<int, int, int>> &ar_params_and_vars,
                         vector<tuple<int, int, int, double>> &additive_vars_params_and_constants) const
{
  ec_params_and_vars.first = -1;

  vector<pair<expr_t, int>> terms;
  decomposeAdditiveTerms(terms, 1);
  for (auto it = terms.begin(); it != terms.end(); ++it)
    if (auto bopn = dynamic_cast<BinaryOpNode *>(it->first); bopn)
      {
        try
          {
            auto [param_id, target_id] = bopn->matchParamTimesTargetMinusVariable(lhs_orig_symb_id);
            ec_params_and_vars = { param_id, { { target_id, true, 1 }, { lhs_orig_symb_id, false, -1 }}};
            terms.erase(it);
            break;
          }
        catch (MatchFailureException &e)
          {
          }
      }

  if (ec_params_and_vars.first < 0)
    {
      cerr << "Error finding EC part of PAC equation" << endl;
      exit(EXIT_FAILURE);
    }

  for (const auto &it : terms)
    {
      if (dynamic_cast<PacExpectationNode *>(it.first))
        continue;

      pair<int, vector<tuple<int, int, int, double>>> m;
      try
        {
          m = {-1, {it.first->matchVariableTimesConstantTimesParam()}};
        }
      catch (MatchFailureException &e)
        {
          try
            {
              m = it.first->matchParamTimesLinearCombinationOfVariables();
            }
          catch (MatchFailureException &e)
            {
              cerr << "Unsupported expression in PAC equation" << endl;
              exit(EXIT_FAILURE);
            }
        }

      for (auto &t : m.second)
        get<3>(t) *= it.second; // Update sign of constants

      int pid = get<0>(m);
      for (auto [vid, vlag, pidtmp, constant] : m.second)
        {
          if (pid == -1)
            pid = pidtmp;
          else if (pidtmp >= 0)
            {
              cerr << "unexpected parameter found in PAC equation" << endl;
              exit(EXIT_FAILURE);
            }

          if (int vidorig = datatree.symbol_table.getUltimateOrigSymbID(vid);
              vidorig == lhs_symb_id || vidorig == lhs_orig_symb_id)
            {
              // This is an autoregressive term
              if (constant != 1 || pid == -1 || !datatree.symbol_table.isDiffAuxiliaryVariable(vid))
                {
                  cerr << "BinaryOpNode::getPacAREC: autoregressive terms must be of the form 'parameter*diff_lagged_variable" << endl;
                  exit(EXIT_FAILURE);
                }
              int ar_lag = datatree.symbol_table.getOrigLeadLagForDiffAuxVar(vid);
              if (static_cast<int>(ar_params_and_vars.size()) < ar_lag)
                ar_params_and_vars.resize(ar_lag, { -1, -1, 0 });
              ar_params_and_vars[ar_lag-1] = { pid, vid, vlag };
            }
          else
            // This is a residual additive term
            additive_vars_params_and_constants.emplace_back(vid, vlag, pid, constant);
        }
    }
}

bool
BinaryOpNode::isParamTimesEndogExpr() const
{
  if (op_code == BinaryOpcode::times)
    {
      set<int> params;
      auto test_arg1 = dynamic_cast<VariableNode *>(arg1);
      auto test_arg2 = dynamic_cast<VariableNode *>(arg2);
      if (test_arg1)
        arg1->collectVariables(SymbolType::parameter, params);
      else if (test_arg2)
        arg2->collectVariables(SymbolType::parameter, params);
      else
        return false;

      if (params.size() != 1)
        return false;

      params.clear();
      set<pair<int, int>> endogs, exogs;
      if (test_arg1)
        {
          arg2->collectDynamicVariables(SymbolType::endogenous, endogs);
          arg2->collectDynamicVariables(SymbolType::exogenous, exogs);
          arg2->collectVariables(SymbolType::parameter, params);
          if (params.size() == 0 && exogs.size() == 0 && endogs.size() >= 1)
            return true;
        }
      else
        {
          arg1->collectDynamicVariables(SymbolType::endogenous, endogs);
          arg1->collectDynamicVariables(SymbolType::exogenous, exogs);
          arg1->collectVariables(SymbolType::parameter, params);
          if (params.size() == 0 && exogs.size() == 0 && endogs.size() >= 1)
            return true;
        }
    }
  else if (op_code == BinaryOpcode::plus)
    return arg1->isParamTimesEndogExpr() || arg2->isParamTimesEndogExpr();
  return false;
}

expr_t
BinaryOpNode::getPacNonOptimizingPart(int optim_share_symb_id) const
{
  vector<pair<expr_t, int>> factors;
  decomposeMultiplicativeFactors(factors);

  // Search for a factor of the form 1-optim_share
  expr_t one_minus_optim_share = nullptr;
  for (auto [factor, exponent] : factors)
    {
      auto bopn = dynamic_cast<BinaryOpNode *>(factor);
      if (exponent != 1 || !bopn || bopn->op_code != BinaryOpcode::minus)
        continue;
      auto arg1 = dynamic_cast<NumConstNode *>(bopn->arg1);
      auto arg2 = dynamic_cast<VariableNode *>(bopn->arg2);
      if (arg1 && arg2 && arg1->eval({}) == 1 && arg2->symb_id == optim_share_symb_id)
        {
          one_minus_optim_share = factor;
          break;
        }
    }

  if (!one_minus_optim_share)
    return nullptr;

  // Construct the product formed by the other factors and return it
  expr_t non_optim_part = datatree.One;
  for (auto [factor, exponent] : factors)
    if (factor != one_minus_optim_share)
      if (exponent == 1)
        non_optim_part = datatree.AddTimes(non_optim_part, factor);
      else
        non_optim_part = datatree.AddDivide(non_optim_part, factor);

  return non_optim_part;
}

pair<int, expr_t>
BinaryOpNode::getPacOptimizingShareAndExprNodesHelper(int lhs_symb_id, int lhs_orig_symb_id) const
{
  int optim_param_symb_id = -1;
  expr_t optim_part = nullptr;
  set<int> endogs;
  collectVariables(SymbolType::endogenous, endogs);
  // Test whether it contains the LHS in level
  if (endogs.count(lhs_orig_symb_id) > 0)
    {
      set<int> params;
      if (arg1->isParamTimesEndogExpr() && !arg2->isParamTimesEndogExpr())
        {
          optim_part = arg1;
          arg2->collectVariables(SymbolType::parameter, params);
          optim_param_symb_id = *(params.begin());
        }
      else if (arg2->isParamTimesEndogExpr() && !arg1->isParamTimesEndogExpr())
        {
          optim_part = arg2;
          arg1->collectVariables(SymbolType::parameter, params);
          optim_param_symb_id = *(params.begin());
        }
    }
  return {optim_param_symb_id, optim_part};
}

tuple<int, expr_t, expr_t, expr_t>
BinaryOpNode::getPacOptimizingShareAndExprNodes(int lhs_symb_id, int lhs_orig_symb_id) const
{
  vector<pair<expr_t, int>> terms;
  decomposeAdditiveTerms(terms, 1);
  for (auto &it : terms)
    if (dynamic_cast<PacExpectationNode *>(it.first))
      // if the pac_expectation operator is additive in the expression
      // there are no optimizing shares
      return {-1, nullptr, nullptr, nullptr};

  int optim_share;
  expr_t optim_part, non_optim_part, additive_part;
  optim_part = non_optim_part = additive_part = nullptr;

  for (auto it = terms.begin(); it != terms.end(); ++it)
    if (auto bopn = dynamic_cast<BinaryOpNode *>(it->first); bopn)
      {
        tie(optim_share, optim_part)
          = bopn->getPacOptimizingShareAndExprNodesHelper(lhs_symb_id, lhs_orig_symb_id);
        if (optim_share >= 0 && optim_part)
          {
            terms.erase(it);
            break;
          }
      }

  if (!optim_part)
    return {-1, nullptr, nullptr, nullptr};

  for (auto it = terms.begin(); it != terms.end(); ++it)
    if (auto bopn = dynamic_cast<BinaryOpNode *>(it->first); bopn)
      {
        non_optim_part = bopn->getPacNonOptimizingPart(optim_share);
        if (non_optim_part)
          {
            terms.erase(it);
            break;
          }
      }

  if (!non_optim_part)
    return {-1, nullptr, nullptr, nullptr};
  else
    {
      additive_part = datatree.Zero;
      for (auto it : terms)
        additive_part = datatree.AddPlus(additive_part, it.first);
      if (additive_part == datatree.Zero)
        additive_part = nullptr;
    }

  return {optim_share, optim_part, non_optim_part, additive_part};
}

void
BinaryOpNode::fillAutoregressiveRow(int eqn, const vector<int> &lhs, map<tuple<int, int, int>, expr_t> &AR) const
{
  vector<pair<expr_t, int>> terms;
  decomposeAdditiveTerms(terms, 1);
  for (const auto &it : terms)
    {
      int vid, lag, param_id;
      double constant;
      try
        {
          tie(vid, lag, param_id, constant) = it.first->matchVariableTimesConstantTimesParam();
          constant *= it.second;
        }
      catch (MatchFailureException &e)
        {
          continue;
        }

      if (datatree.symbol_table.isDiffAuxiliaryVariable(vid))
        {
          lag = -datatree.symbol_table.getOrigLeadLagForDiffAuxVar(vid);
          vid = datatree.symbol_table.getOrigSymbIdForDiffAuxVar(vid);
        }

      if (find(lhs.begin(), lhs.end(), vid) == lhs.end())
        continue;

      if (AR.find({eqn, -lag, vid}) != AR.end())
        {
          cerr << "BinaryOpNode::fillAutoregressiveRow: Error filling AR matrix: "
               << "lag/symb_id encountered more than once in equation" << endl;
          exit(EXIT_FAILURE);
        }
      if (constant != 1 || param_id == -1)
        {
          cerr << "BinaryOpNode::fillAutoregressiveRow: autoregressive terms must be of the form 'parameter*lagged_variable" << endl;
          exit(EXIT_FAILURE);
        }
      AR[{eqn, -lag, vid}] = datatree.AddVariable(param_id);
    }
}

void
BinaryOpNode::findConstantEquations(map<VariableNode *, NumConstNode *> &table) const
{
  if (op_code == BinaryOpcode::equal)
    {
      // The variable must be contemporaneous (see #83)
      if (auto varg1 = dynamic_cast<VariableNode *>(arg1);
          varg1 && varg1->lag == 0 && dynamic_cast<NumConstNode *>(arg2))
        table[varg1] = dynamic_cast<NumConstNode *>(arg2);
      else if (auto varg2 = dynamic_cast<VariableNode *>(arg2);
               varg2 && varg2->lag == 0 && dynamic_cast<NumConstNode *>(arg1))
        table[varg2] = dynamic_cast<NumConstNode *>(arg1);
    }
}

expr_t
BinaryOpNode::replaceVarsInEquation(map<VariableNode *, NumConstNode *> &table) const
{
  if (op_code == BinaryOpcode::equal)
    for (auto &it : table)
      if ((it.first == arg1 && it.second == arg2) || (it.first == arg2 && it.second == arg1))
        return const_cast<BinaryOpNode *>(this);
  expr_t arg1subst = arg1->replaceVarsInEquation(table);
  expr_t arg2subst = arg2->replaceVarsInEquation(table);
  return buildSimilarBinaryOpNode(arg1subst, arg2subst, datatree);
}

expr_t
BinaryOpNode::substituteStaticAuxiliaryVariable() const
{
  expr_t arg1subst = arg1->substituteStaticAuxiliaryVariable();
  expr_t arg2subst = arg2->substituteStaticAuxiliaryVariable();
  return buildSimilarBinaryOpNode(arg1subst, arg2subst, datatree);
}

expr_t
BinaryOpNode::substituteStaticAuxiliaryDefinition() const
{
  expr_t arg2subst = arg2->substituteStaticAuxiliaryVariable();
  return buildSimilarBinaryOpNode(arg1, arg2subst, datatree);
}

void
BinaryOpNode::matchMatchedMoment(vector<int> &symb_ids, vector<int> &lags, vector<int> &powers) const
{
  if (op_code == BinaryOpcode::times)
    {
      arg1->matchMatchedMoment(symb_ids, lags, powers);
      arg2->matchMatchedMoment(symb_ids, lags, powers);
    }
  else if (op_code == BinaryOpcode::power)
    {
      if (!dynamic_cast<const VariableNode *>(arg1))
        throw MatchFailureException("First argument of power expression must be a variable");
      auto ncn = dynamic_cast<const NumConstNode *>(arg2);
      if (!ncn)
        throw MatchFailureException("Second argument of power expression must be a positive integer");
      double c = datatree.num_constants.getDouble(ncn->id);
      if (c <= 0 || round(c) != c)
        throw MatchFailureException("Second argument of power expression must be a positive integer");
      arg1->matchMatchedMoment(symb_ids, lags, powers);
      powers.back() = static_cast<int>(c);
    }
  else
    throw MatchFailureException("Unsupported binary operator");
}


TrinaryOpNode::TrinaryOpNode(DataTree &datatree_arg, int idx_arg, const expr_t arg1_arg,
                             TrinaryOpcode op_code_arg, const expr_t arg2_arg, const expr_t arg3_arg) :
  ExprNode{datatree_arg, idx_arg},
  arg1{arg1_arg},
  arg2{arg2_arg},
  arg3{arg3_arg},
  op_code{op_code_arg}
{
}

void
TrinaryOpNode::prepareForDerivation()
{
  if (preparedForDerivation)
    return;

  preparedForDerivation = true;

  arg1->prepareForDerivation();
  arg2->prepareForDerivation();
  arg3->prepareForDerivation();

  // Non-null derivatives are the union of those of the arguments
  // Compute set union of arg{1,2,3}->non_null_derivatives
  set<int> non_null_derivatives_tmp;
  set_union(arg1->non_null_derivatives.begin(),
            arg1->non_null_derivatives.end(),
            arg2->non_null_derivatives.begin(),
            arg2->non_null_derivatives.end(),
            inserter(non_null_derivatives_tmp, non_null_derivatives_tmp.begin()));
  set_union(non_null_derivatives_tmp.begin(),
            non_null_derivatives_tmp.end(),
            arg3->non_null_derivatives.begin(),
            arg3->non_null_derivatives.end(),
            inserter(non_null_derivatives, non_null_derivatives.begin()));
}

expr_t
TrinaryOpNode::composeDerivatives(expr_t darg1, expr_t darg2, expr_t darg3)
{

  expr_t t11, t12, t13, t14, t15;

  switch (op_code)
    {
    case TrinaryOpcode::normcdf:
      // normal pdf is inlined in the tree
      expr_t y;
      // sqrt(2*pi)
      t14 = datatree.AddSqrt(datatree.AddTimes(datatree.Two, datatree.Pi));
      // x - mu
      t12 = datatree.AddMinus(arg1, arg2);
      // y = (x-mu)/sigma
      y = datatree.AddDivide(t12, arg3);
      // (x-mu)^2/sigma^2
      t12 = datatree.AddTimes(y, y);
      // -(x-mu)^2/sigma^2
      t13 = datatree.AddUMinus(t12);
      // -((x-mu)^2/sigma^2)/2
      t12 = datatree.AddDivide(t13, datatree.Two);
      // exp(-((x-mu)^2/sigma^2)/2)
      t13 = datatree.AddExp(t12);
      // derivative of a standardized normal
      // t15 = (1/sqrt(2*pi))*exp(-y^2/2)
      t15 = datatree.AddDivide(t13, t14);
      // derivatives thru x
      t11 = datatree.AddDivide(darg1, arg3);
      // derivatives thru mu
      t12 = datatree.AddDivide(darg2, arg3);
      // intermediary sum
      t14 = datatree.AddMinus(t11, t12);
      // derivatives thru sigma
      t11 = datatree.AddDivide(y, arg3);
      t12 = datatree.AddTimes(t11, darg3);
      //intermediary sum
      t11 = datatree.AddMinus(t14, t12);
      // total derivative:
      // (darg1/sigma - darg2/sigma - darg3*(x-mu)/sigma^2) * t15
      // where t15 is the derivative of a standardized normal
      return datatree.AddTimes(t11, t15);
    case TrinaryOpcode::normpdf:
      // (x - mu)
      t11 = datatree.AddMinus(arg1, arg2);
      // (x - mu)/sigma
      t12 = datatree.AddDivide(t11, arg3);
      // darg3 * (x - mu)/sigma
      t11 = datatree.AddTimes(darg3, t12);
      // darg2 - darg1
      t13 = datatree.AddMinus(darg2, darg1);
      // darg2 - darg1 + darg3 * (x - mu)/sigma
      t14 = datatree.AddPlus(t13, t11);
      // ((x - mu)/sigma) * (darg2 - darg1 + darg3 * (x - mu)/sigma)
      t11 = datatree.AddTimes(t12, t14);
      // ((x - mu)/sigma) * (darg2 - darg1 + darg3 * (x - mu)/sigma) - darg3
      t12 = datatree.AddMinus(t11, darg3);
      // this / sigma
      t11 = datatree.AddDivide(this, arg3);
      // total derivative:
      // (this / sigma) * (((x - mu)/sigma) * (darg2 - darg1 + darg3 * (x - mu)/sigma) - darg3)
      return datatree.AddTimes(t11, t12);
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

expr_t
TrinaryOpNode::computeDerivative(int deriv_id)
{
  expr_t darg1 = arg1->getDerivative(deriv_id);
  expr_t darg2 = arg2->getDerivative(deriv_id);
  expr_t darg3 = arg3->getDerivative(deriv_id);
  return composeDerivatives(darg1, darg2, darg3);
}

int
TrinaryOpNode::precedence(ExprNodeOutputType output_type, const temporary_terms_t &temporary_terms) const
{
  // A temporary term behaves as a variable
  if (temporary_terms.find(const_cast<TrinaryOpNode *>(this)) != temporary_terms.end())
    return 100;

  switch (op_code)
    {
    case TrinaryOpcode::normcdf:
    case TrinaryOpcode::normpdf:
      return 100;
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

int
TrinaryOpNode::cost(const map<pair<int, int>, temporary_terms_t> &temp_terms_map, bool is_matlab) const
{
  // For a temporary term, the cost is null
  for (const auto &it : temp_terms_map)
    if (it.second.find(const_cast<TrinaryOpNode *>(this)) != it.second.end())
      return 0;

  int arg_cost = arg1->cost(temp_terms_map, is_matlab)
    + arg2->cost(temp_terms_map, is_matlab)
    + arg3->cost(temp_terms_map, is_matlab);

  return cost(arg_cost, is_matlab);
}

int
TrinaryOpNode::cost(const vector<vector<temporary_terms_t>> &blocks_temporary_terms, bool is_matlab) const
{
  // For a temporary term, the cost is null
  for (const auto &blk_tt : blocks_temporary_terms)
    for (const auto &eq_tt : blk_tt)
      if (eq_tt.find(const_cast<TrinaryOpNode *>(this)) != eq_tt.end())
        return 0;

  int arg_cost = arg1->cost(blocks_temporary_terms, is_matlab)
    + arg2->cost(blocks_temporary_terms, is_matlab)
    + arg3->cost(blocks_temporary_terms, is_matlab);

  return cost(arg_cost, is_matlab);
}

int
TrinaryOpNode::cost(int cost, bool is_matlab) const
{
  if (is_matlab)
    // Cost for Matlab files
    switch (op_code)
      {
      case TrinaryOpcode::normcdf:
      case TrinaryOpcode::normpdf:
        return cost+1000;
      }
  else
    // Cost for C files
    switch (op_code)
      {
      case TrinaryOpcode::normcdf:
      case TrinaryOpcode::normpdf:
        return cost+1000;
      }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

void
TrinaryOpNode::computeTemporaryTerms(const pair<int, int> &derivOrder,
                                     map<pair<int, int>, temporary_terms_t> &temp_terms_map,
                                     map<expr_t, pair<int, pair<int, int>>> &reference_count,
                                     bool is_matlab) const
{
  expr_t this2 = const_cast<TrinaryOpNode *>(this);
  if (auto it = reference_count.find(this2);
      it == reference_count.end())
    {
      // If this node has never been encountered, set its ref count to one,
      //  and travel through its children
      reference_count[this2] = { 1, derivOrder };
      arg1->computeTemporaryTerms(derivOrder, temp_terms_map, reference_count, is_matlab);
      arg2->computeTemporaryTerms(derivOrder, temp_terms_map, reference_count, is_matlab);
      arg3->computeTemporaryTerms(derivOrder, temp_terms_map, reference_count, is_matlab);
    }
  else
    {
      // If the node has already been encountered, increment its ref count
      //  and declare it as a temporary term if it is too costly
      auto &[nref, min_order] = it->second;
      nref++;
      if (nref * cost(temp_terms_map, is_matlab) > min_cost(is_matlab))
        temp_terms_map[min_order].insert(this2);
    }
}

void
TrinaryOpNode::computeBlockTemporaryTerms(int blk, int eq, vector<vector<temporary_terms_t>> &blocks_temporary_terms,
                                          map<expr_t, tuple<int, int, int>> &reference_count) const
{
  expr_t this2 = const_cast<TrinaryOpNode *>(this);
  if (auto it = reference_count.find(this2);
      it == reference_count.end())
    {
      reference_count[this2] = { 1, blk, eq };
      arg1->computeBlockTemporaryTerms(blk, eq, blocks_temporary_terms, reference_count);
      arg2->computeBlockTemporaryTerms(blk, eq, blocks_temporary_terms, reference_count);
      arg3->computeBlockTemporaryTerms(blk, eq, blocks_temporary_terms, reference_count);
    }
  else
    {
      auto &[nref, first_blk, first_eq] = it->second;
      nref++;
      if (nref * cost(blocks_temporary_terms, false) > min_cost_c)
        blocks_temporary_terms[first_blk][first_eq].insert(this2);
    }
}

double
TrinaryOpNode::eval_opcode(double v1, TrinaryOpcode op_code, double v2, double v3) noexcept(false)
{
  switch (op_code)
    {
    case TrinaryOpcode::normcdf:
      return (0.5*(1+erf((v1-v2)/v3/M_SQRT2)));
    case TrinaryOpcode::normpdf:
      return (1/(v3*sqrt(2*M_PI)*exp(pow((v1-v2)/v3, 2)/2)));
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

double
TrinaryOpNode::eval(const eval_context_t &eval_context) const noexcept(false)
{
  double v1 = arg1->eval(eval_context);
  double v2 = arg2->eval(eval_context);
  double v3 = arg3->eval(eval_context);

  return eval_opcode(v1, op_code, v2, v3);
}

void
TrinaryOpNode::compile(ostream &CompileCode, unsigned int &instruction_number,
                       bool lhs_rhs, const temporary_terms_t &temporary_terms,
                       const temporary_terms_idxs_t &temporary_terms_idxs, bool dynamic, bool steady_dynamic,
                       const deriv_node_temp_terms_t &tef_terms) const
{
  // If current node is a temporary term
  if (auto this2 = const_cast<TrinaryOpNode *>(this);
      temporary_terms.find(this2) != temporary_terms.end())
    {
      if (dynamic)
        {
          FLDT_ fldt(temporary_terms_idxs.at(this2));
          fldt.write(CompileCode, instruction_number);
        }
      else
        {
          FLDST_ fldst(temporary_terms_idxs.at(this2));
          fldst.write(CompileCode, instruction_number);
        }
      return;
    }
  arg1->compile(CompileCode, instruction_number, lhs_rhs, temporary_terms, temporary_terms_idxs, dynamic, steady_dynamic, tef_terms);
  arg2->compile(CompileCode, instruction_number, lhs_rhs, temporary_terms, temporary_terms_idxs, dynamic, steady_dynamic, tef_terms);
  arg3->compile(CompileCode, instruction_number, lhs_rhs, temporary_terms, temporary_terms_idxs, dynamic, steady_dynamic, tef_terms);
  FTRINARY_ ftrinary{static_cast<int>(op_code)};
  ftrinary.write(CompileCode, instruction_number);
}

bool
TrinaryOpNode::containsExternalFunction() const
{
  return arg1->containsExternalFunction()
    || arg2->containsExternalFunction()
    || arg3->containsExternalFunction();
}

void
TrinaryOpNode::writeJsonAST(ostream &output) const
{
  output << R"({"node_type" : "TrinaryOpNode", )"
         << R"("op" : ")";
  switch (op_code)
    {
    case TrinaryOpcode::normcdf:
      output << "normcdf";
      break;
    case TrinaryOpcode::normpdf:
      output << "normpdf";
      break;
    }
  output << R"(", "arg1" : )";
  arg1->writeJsonAST(output);
  output << R"(, "arg2" : )";
  arg2->writeJsonAST(output);
  output << R"(, "arg2" : )";
  arg3->writeJsonAST(output);
  output << "}";
}

void
TrinaryOpNode::writeJsonOutput(ostream &output,
                               const temporary_terms_t &temporary_terms,
                               const deriv_node_temp_terms_t &tef_terms,
                               bool isdynamic) const
{
  // If current node is a temporary term
  if (temporary_terms.find(const_cast<TrinaryOpNode *>(this)) != temporary_terms.end())
    {
      output << "T" << idx;
      return;
    }

  switch (op_code)
    {
    case TrinaryOpcode::normcdf:
      output << "normcdf(";
      break;
    case TrinaryOpcode::normpdf:
      output << "normpdf(";
      break;
    }

  arg1->writeJsonOutput(output, temporary_terms, tef_terms, isdynamic);
  output << ",";
  arg2->writeJsonOutput(output, temporary_terms, tef_terms, isdynamic);
  output << ",";
  arg3->writeJsonOutput(output, temporary_terms, tef_terms, isdynamic);
  output << ")";
}

void
TrinaryOpNode::writeOutput(ostream &output, ExprNodeOutputType output_type,
                           const temporary_terms_t &temporary_terms,
                           const temporary_terms_idxs_t &temporary_terms_idxs,
                           const deriv_node_temp_terms_t &tef_terms) const
{
  if (checkIfTemporaryTermThenWrite(output, output_type, temporary_terms, temporary_terms_idxs))
    return;

  switch (op_code)
    {
    case TrinaryOpcode::normcdf:
      if (isCOutput(output_type))
        {
          // In C, there is no normcdf() primitive, so use erf()
          output << "(0.5*(1+erf(((";
          arg1->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
          output << ")-(";
          arg2->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
          output << "))/(";
          arg3->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
          output << ")/M_SQRT2)))";
        }
      else if (isJuliaOutput(output_type))
	{
	  // Julia API is normcdf(mu, sigma, x) !
          output << "normcdf(";
          arg2->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
          output << ",";
          arg3->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
          output << ",";
          arg1->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
          output << ")";
        }
      else
        {
          output << "normcdf(";
          arg1->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
          output << ",";
          arg2->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
          output << ",";
          arg3->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
          output << ")";
        }
      break;
    case TrinaryOpcode::normpdf:
      if (isCOutput(output_type))
        {
          //(1/(v3*sqrt(2*M_PI)*exp(pow((v1-v2)/v3,2)/2)))
          output << "(1/(";
          arg3->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
          output << "*sqrt(2*M_PI)*exp(pow((";
          arg1->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
          output << "-";
          arg2->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
          output << ")/";
          arg3->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
          output << ",2)/2)))";
        }
      else if (isJuliaOutput(output_type))
	{
	  // Julia API is normpdf(mu, sigma, x) !
          output << "normpdf(";
          arg2->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
          output << ",";
          arg3->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
          output << ",";
          arg1->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
          output << ")";
        }
      else
        {
          output << "normpdf(";
          arg1->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
          output << ",";
          arg2->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
          output << ",";
          arg3->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
          output << ")";
        }
      break;
    }
}

void
TrinaryOpNode::writeExternalFunctionOutput(ostream &output, ExprNodeOutputType output_type,
                                           const temporary_terms_t &temporary_terms,
                                           const temporary_terms_idxs_t &temporary_terms_idxs,
                                           deriv_node_temp_terms_t &tef_terms) const
{
  arg1->writeExternalFunctionOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
  arg2->writeExternalFunctionOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
  arg3->writeExternalFunctionOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
}

void
TrinaryOpNode::writeJsonExternalFunctionOutput(vector<string> &efout,
                                               const temporary_terms_t &temporary_terms,
                                               deriv_node_temp_terms_t &tef_terms,
                                               bool isdynamic) const
{
  arg1->writeJsonExternalFunctionOutput(efout, temporary_terms, tef_terms, isdynamic);
  arg2->writeJsonExternalFunctionOutput(efout, temporary_terms, tef_terms, isdynamic);
  arg3->writeJsonExternalFunctionOutput(efout, temporary_terms, tef_terms, isdynamic);
}

void
TrinaryOpNode::compileExternalFunctionOutput(ostream &CompileCode, unsigned int &instruction_number,
                                             bool lhs_rhs, const temporary_terms_t &temporary_terms,
                                             const temporary_terms_idxs_t &temporary_terms_idxs, bool dynamic, bool steady_dynamic,
                                             deriv_node_temp_terms_t &tef_terms) const
{
  arg1->compileExternalFunctionOutput(CompileCode, instruction_number, lhs_rhs, temporary_terms,
                                      temporary_terms_idxs, dynamic, steady_dynamic, tef_terms);
  arg2->compileExternalFunctionOutput(CompileCode, instruction_number, lhs_rhs, temporary_terms,
                                      temporary_terms_idxs, dynamic, steady_dynamic, tef_terms);
  arg3->compileExternalFunctionOutput(CompileCode, instruction_number, lhs_rhs, temporary_terms,
                                      temporary_terms_idxs, dynamic, steady_dynamic, tef_terms);
}

void
TrinaryOpNode::collectVARLHSVariable(set<expr_t> &result) const
{
  cerr << "ERROR: you can only have variables or unary ops on LHS of VAR" << endl;
  exit(EXIT_FAILURE);
}

void
TrinaryOpNode::collectDynamicVariables(SymbolType type_arg, set<pair<int, int>> &result) const
{
  arg1->collectDynamicVariables(type_arg, result);
  arg2->collectDynamicVariables(type_arg, result);
  arg3->collectDynamicVariables(type_arg, result);
}

void
TrinaryOpNode::computeSubExprContainingVariable(int symb_id, int lag, set<expr_t> &contain_var) const
{
  arg1->computeSubExprContainingVariable(symb_id, lag, contain_var);
  arg2->computeSubExprContainingVariable(symb_id, lag, contain_var);
  arg3->computeSubExprContainingVariable(symb_id, lag, contain_var);
  if (contain_var.count(arg1) > 0 || contain_var.count(arg2) > 0 || contain_var.count(arg3) > 0)
    contain_var.insert(const_cast<TrinaryOpNode *>(this));
}

BinaryOpNode *
TrinaryOpNode::normalizeEquationHelper(const set<expr_t> &contain_var, expr_t rhs) const
{
  throw NormalizationFailed();
}

expr_t
TrinaryOpNode::getChainRuleDerivative(int deriv_id, const map<int, BinaryOpNode *> &recursive_variables)
{
  expr_t darg1 = arg1->getChainRuleDerivative(deriv_id, recursive_variables);
  expr_t darg2 = arg2->getChainRuleDerivative(deriv_id, recursive_variables);
  expr_t darg3 = arg3->getChainRuleDerivative(deriv_id, recursive_variables);
  return composeDerivatives(darg1, darg2, darg3);
}

expr_t
TrinaryOpNode::buildSimilarTrinaryOpNode(expr_t alt_arg1, expr_t alt_arg2, expr_t alt_arg3, DataTree &alt_datatree) const
{
  switch (op_code)
    {
    case TrinaryOpcode::normcdf:
      return alt_datatree.AddNormcdf(alt_arg1, alt_arg2, alt_arg3);
    case TrinaryOpcode::normpdf:
      return alt_datatree.AddNormpdf(alt_arg1, alt_arg2, alt_arg3);
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

expr_t
TrinaryOpNode::toStatic(DataTree &static_datatree) const
{
  expr_t sarg1 = arg1->toStatic(static_datatree);
  expr_t sarg2 = arg2->toStatic(static_datatree);
  expr_t sarg3 = arg3->toStatic(static_datatree);
  return buildSimilarTrinaryOpNode(sarg1, sarg2, sarg3, static_datatree);
}

void
TrinaryOpNode::computeXrefs(EquationInfo &ei) const
{
  arg1->computeXrefs(ei);
  arg2->computeXrefs(ei);
  arg3->computeXrefs(ei);
}

expr_t
TrinaryOpNode::clone(DataTree &datatree) const
{
  expr_t substarg1 = arg1->clone(datatree);
  expr_t substarg2 = arg2->clone(datatree);
  expr_t substarg3 = arg3->clone(datatree);
  return buildSimilarTrinaryOpNode(substarg1, substarg2, substarg3, datatree);
}

int
TrinaryOpNode::maxEndoLead() const
{
  return max(arg1->maxEndoLead(), max(arg2->maxEndoLead(), arg3->maxEndoLead()));
}

int
TrinaryOpNode::maxExoLead() const
{
  return max(arg1->maxExoLead(), max(arg2->maxExoLead(), arg3->maxExoLead()));
}

int
TrinaryOpNode::maxEndoLag() const
{
  return max(arg1->maxEndoLag(), max(arg2->maxEndoLag(), arg3->maxEndoLag()));
}

int
TrinaryOpNode::maxExoLag() const
{
  return max(arg1->maxExoLag(), max(arg2->maxExoLag(), arg3->maxExoLag()));
}

int
TrinaryOpNode::maxLead() const
{
  return max(arg1->maxLead(), max(arg2->maxLead(), arg3->maxLead()));
}

int
TrinaryOpNode::maxLag() const
{
  return max(arg1->maxLag(), max(arg2->maxLag(), arg3->maxLag()));
}

int
TrinaryOpNode::maxLagWithDiffsExpanded() const
{
  return max(arg1->maxLagWithDiffsExpanded(),
             max(arg2->maxLagWithDiffsExpanded(), arg3->maxLagWithDiffsExpanded()));
}

expr_t
TrinaryOpNode::undiff() const
{
  expr_t arg1subst = arg1->undiff();
  expr_t arg2subst = arg2->undiff();
  expr_t arg3subst = arg3->undiff();
  return buildSimilarTrinaryOpNode(arg1subst, arg2subst, arg3subst, datatree);
}

int
TrinaryOpNode::VarMaxLag(const set<expr_t> &lhs_lag_equiv) const
{
  return max(arg1->VarMaxLag(lhs_lag_equiv),
             max(arg2->VarMaxLag(lhs_lag_equiv),
                 arg3->VarMaxLag(lhs_lag_equiv)));
}

expr_t
TrinaryOpNode::decreaseLeadsLags(int n) const
{
  expr_t arg1subst = arg1->decreaseLeadsLags(n);
  expr_t arg2subst = arg2->decreaseLeadsLags(n);
  expr_t arg3subst = arg3->decreaseLeadsLags(n);
  return buildSimilarTrinaryOpNode(arg1subst, arg2subst, arg3subst, datatree);
}

expr_t
TrinaryOpNode::decreaseLeadsLagsPredeterminedVariables() const
{
  expr_t arg1subst = arg1->decreaseLeadsLagsPredeterminedVariables();
  expr_t arg2subst = arg2->decreaseLeadsLagsPredeterminedVariables();
  expr_t arg3subst = arg3->decreaseLeadsLagsPredeterminedVariables();
  return buildSimilarTrinaryOpNode(arg1subst, arg2subst, arg3subst, datatree);
}

expr_t
TrinaryOpNode::substituteEndoLeadGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const
{
  if (maxEndoLead() < 2)
    return const_cast<TrinaryOpNode *>(this);
  else if (deterministic_model)
    {
      expr_t arg1subst = arg1->substituteEndoLeadGreaterThanTwo(subst_table, neweqs, deterministic_model);
      expr_t arg2subst = arg2->substituteEndoLeadGreaterThanTwo(subst_table, neweqs, deterministic_model);
      expr_t arg3subst = arg3->substituteEndoLeadGreaterThanTwo(subst_table, neweqs, deterministic_model);
      return buildSimilarTrinaryOpNode(arg1subst, arg2subst, arg3subst, datatree);
    }
  else
    return createEndoLeadAuxiliaryVarForMyself(subst_table, neweqs);
}

expr_t
TrinaryOpNode::substituteEndoLagGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  expr_t arg1subst = arg1->substituteEndoLagGreaterThanTwo(subst_table, neweqs);
  expr_t arg2subst = arg2->substituteEndoLagGreaterThanTwo(subst_table, neweqs);
  expr_t arg3subst = arg3->substituteEndoLagGreaterThanTwo(subst_table, neweqs);
  return buildSimilarTrinaryOpNode(arg1subst, arg2subst, arg3subst, datatree);
}

expr_t
TrinaryOpNode::substituteExoLead(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const
{
  if (maxExoLead() == 0)
    return const_cast<TrinaryOpNode *>(this);
  else if (deterministic_model)
    {
      expr_t arg1subst = arg1->substituteExoLead(subst_table, neweqs, deterministic_model);
      expr_t arg2subst = arg2->substituteExoLead(subst_table, neweqs, deterministic_model);
      expr_t arg3subst = arg3->substituteExoLead(subst_table, neweqs, deterministic_model);
      return buildSimilarTrinaryOpNode(arg1subst, arg2subst, arg3subst, datatree);
    }
  else
    return createExoLeadAuxiliaryVarForMyself(subst_table, neweqs);
}

expr_t
TrinaryOpNode::substituteExoLag(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  expr_t arg1subst = arg1->substituteExoLag(subst_table, neweqs);
  expr_t arg2subst = arg2->substituteExoLag(subst_table, neweqs);
  expr_t arg3subst = arg3->substituteExoLag(subst_table, neweqs);
  return buildSimilarTrinaryOpNode(arg1subst, arg2subst, arg3subst, datatree);
}

expr_t
TrinaryOpNode::substituteExpectation(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool partial_information_model) const
{
  expr_t arg1subst = arg1->substituteExpectation(subst_table, neweqs, partial_information_model);
  expr_t arg2subst = arg2->substituteExpectation(subst_table, neweqs, partial_information_model);
  expr_t arg3subst = arg3->substituteExpectation(subst_table, neweqs, partial_information_model);
  return buildSimilarTrinaryOpNode(arg1subst, arg2subst, arg3subst, datatree);
}

expr_t
TrinaryOpNode::substituteAdl() const
{
  expr_t arg1subst = arg1->substituteAdl();
  expr_t arg2subst = arg2->substituteAdl();
  expr_t arg3subst = arg3->substituteAdl();
  return buildSimilarTrinaryOpNode(arg1subst, arg2subst, arg3subst, datatree);
}

expr_t
TrinaryOpNode::substituteModelLocalVariables() const
{
  expr_t arg1subst = arg1->substituteModelLocalVariables();
  expr_t arg2subst = arg2->substituteModelLocalVariables();
  expr_t arg3subst = arg3->substituteModelLocalVariables();
  return buildSimilarTrinaryOpNode(arg1subst, arg2subst, arg3subst, datatree);
}

expr_t
TrinaryOpNode::substituteVarExpectation(const map<string, expr_t> &subst_table) const
{
  expr_t arg1subst = arg1->substituteVarExpectation(subst_table);
  expr_t arg2subst = arg2->substituteVarExpectation(subst_table);
  expr_t arg3subst = arg3->substituteVarExpectation(subst_table);
  return buildSimilarTrinaryOpNode(arg1subst, arg2subst, arg3subst, datatree);
}

void
TrinaryOpNode::findDiffNodes(lag_equivalence_table_t &nodes) const
{
  arg1->findDiffNodes(nodes);
  arg2->findDiffNodes(nodes);
  arg3->findDiffNodes(nodes);
}

void
TrinaryOpNode::findUnaryOpNodesForAuxVarCreation(lag_equivalence_table_t &nodes) const
{
  arg1->findUnaryOpNodesForAuxVarCreation(nodes);
  arg2->findUnaryOpNodesForAuxVarCreation(nodes);
  arg3->findUnaryOpNodesForAuxVarCreation(nodes);
}

int
TrinaryOpNode::findTargetVariable(int lhs_symb_id) const
{
  int retval = arg1->findTargetVariable(lhs_symb_id);
  if (retval < 0)
    retval = arg2->findTargetVariable(lhs_symb_id);
  if (retval < 0)
    retval = arg3->findTargetVariable(lhs_symb_id);
  return retval;
}

expr_t
TrinaryOpNode::substituteDiff(const lag_equivalence_table_t &nodes, subst_table_t &subst_table,
                              vector<BinaryOpNode *> &neweqs) const
{
  expr_t arg1subst = arg1->substituteDiff(nodes, subst_table, neweqs);
  expr_t arg2subst = arg2->substituteDiff(nodes, subst_table, neweqs);
  expr_t arg3subst = arg3->substituteDiff(nodes, subst_table, neweqs);
  return buildSimilarTrinaryOpNode(arg1subst, arg2subst, arg3subst, datatree);
}

expr_t
TrinaryOpNode::substituteUnaryOpNodes(const lag_equivalence_table_t &nodes, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  expr_t arg1subst = arg1->substituteUnaryOpNodes(nodes, subst_table, neweqs);
  expr_t arg2subst = arg2->substituteUnaryOpNodes(nodes, subst_table, neweqs);
  expr_t arg3subst = arg3->substituteUnaryOpNodes(nodes, subst_table, neweqs);
  return buildSimilarTrinaryOpNode(arg1subst, arg2subst, arg3subst, datatree);
}

int
TrinaryOpNode::countDiffs() const
{
  return max(arg1->countDiffs(), max(arg2->countDiffs(), arg3->countDiffs()));
}

expr_t
TrinaryOpNode::substitutePacExpectation(const string &name, expr_t subexpr)
{
  expr_t arg1subst = arg1->substitutePacExpectation(name, subexpr);
  expr_t arg2subst = arg2->substitutePacExpectation(name, subexpr);
  expr_t arg3subst = arg3->substitutePacExpectation(name, subexpr);
  return buildSimilarTrinaryOpNode(arg1subst, arg2subst, arg3subst, datatree);
}

expr_t
TrinaryOpNode::substitutePacTargetNonstationary(const string &name, expr_t subexpr)
{
  expr_t arg1subst = arg1->substitutePacTargetNonstationary(name, subexpr);
  expr_t arg2subst = arg2->substitutePacTargetNonstationary(name, subexpr);
  expr_t arg3subst = arg3->substitutePacTargetNonstationary(name, subexpr);
  return buildSimilarTrinaryOpNode(arg1subst, arg2subst, arg3subst, datatree);
}

expr_t
TrinaryOpNode::differentiateForwardVars(const vector<string> &subset, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  expr_t arg1subst = arg1->differentiateForwardVars(subset, subst_table, neweqs);
  expr_t arg2subst = arg2->differentiateForwardVars(subset, subst_table, neweqs);
  expr_t arg3subst = arg3->differentiateForwardVars(subset, subst_table, neweqs);
  return buildSimilarTrinaryOpNode(arg1subst, arg2subst, arg3subst, datatree);
}

bool
TrinaryOpNode::isNumConstNodeEqualTo(double value) const
{
  return false;
}

bool
TrinaryOpNode::isVariableNodeEqualTo(SymbolType type_arg, int variable_id, int lag_arg) const
{
  return false;
}

bool
TrinaryOpNode::containsPacExpectation(const string &pac_model_name) const
{
  return (arg1->containsPacExpectation(pac_model_name) || arg2->containsPacExpectation(pac_model_name) || arg3->containsPacExpectation(pac_model_name));
}

bool
TrinaryOpNode::containsPacTargetNonstationary(const string &pac_model_name) const
{
  return arg1->containsPacTargetNonstationary(pac_model_name)
    || arg2->containsPacTargetNonstationary(pac_model_name)
    || arg3->containsPacTargetNonstationary(pac_model_name);
}

expr_t
TrinaryOpNode::replaceTrendVar() const
{
  expr_t arg1subst = arg1->replaceTrendVar();
  expr_t arg2subst = arg2->replaceTrendVar();
  expr_t arg3subst = arg3->replaceTrendVar();
  return buildSimilarTrinaryOpNode(arg1subst, arg2subst, arg3subst, datatree);
}

expr_t
TrinaryOpNode::detrend(int symb_id, bool log_trend, expr_t trend) const
{
  expr_t arg1subst = arg1->detrend(symb_id, log_trend, trend);
  expr_t arg2subst = arg2->detrend(symb_id, log_trend, trend);
  expr_t arg3subst = arg3->detrend(symb_id, log_trend, trend);
  return buildSimilarTrinaryOpNode(arg1subst, arg2subst, arg3subst, datatree);
}

expr_t
TrinaryOpNode::removeTrendLeadLag(const map<int, expr_t> &trend_symbols_map) const
{
  expr_t arg1subst = arg1->removeTrendLeadLag(trend_symbols_map);
  expr_t arg2subst = arg2->removeTrendLeadLag(trend_symbols_map);
  expr_t arg3subst = arg3->removeTrendLeadLag(trend_symbols_map);
  return buildSimilarTrinaryOpNode(arg1subst, arg2subst, arg3subst, datatree);
}

bool
TrinaryOpNode::isInStaticForm() const
{
  return arg1->isInStaticForm() && arg2->isInStaticForm() && arg3->isInStaticForm();
}

bool
TrinaryOpNode::isParamTimesEndogExpr() const
{
  return arg1->isParamTimesEndogExpr()
    || arg2->isParamTimesEndogExpr()
    || arg3->isParamTimesEndogExpr();
}

expr_t
TrinaryOpNode::substituteStaticAuxiliaryVariable() const
{
  expr_t arg1subst = arg1->substituteStaticAuxiliaryVariable();
  expr_t arg2subst = arg2->substituteStaticAuxiliaryVariable();
  expr_t arg3subst = arg3->substituteStaticAuxiliaryVariable();
  return buildSimilarTrinaryOpNode(arg1subst, arg2subst, arg3subst, datatree);
}

expr_t
TrinaryOpNode::replaceVarsInEquation(map<VariableNode *, NumConstNode *> &table) const
{
  expr_t arg1subst = arg1->replaceVarsInEquation(table);
  expr_t arg2subst = arg2->replaceVarsInEquation(table);
  expr_t arg3subst = arg3->replaceVarsInEquation(table);
  return buildSimilarTrinaryOpNode(arg1subst, arg2subst, arg3subst, datatree);
}

AbstractExternalFunctionNode::AbstractExternalFunctionNode(DataTree &datatree_arg,
                                                           int idx_arg,
                                                           int symb_id_arg,
                                                           vector<expr_t> arguments_arg) :
  ExprNode{datatree_arg, idx_arg},
  symb_id{symb_id_arg},
  arguments{move(arguments_arg)}
{
}

void
AbstractExternalFunctionNode::prepareForDerivation()
{
  if (preparedForDerivation)
    return;

  for (auto argument : arguments)
    argument->prepareForDerivation();

  non_null_derivatives = arguments.at(0)->non_null_derivatives;
  for (int i = 1; i < static_cast<int>(arguments.size()); i++)
    set_union(non_null_derivatives.begin(),
              non_null_derivatives.end(),
              arguments.at(i)->non_null_derivatives.begin(),
              arguments.at(i)->non_null_derivatives.end(),
              inserter(non_null_derivatives, non_null_derivatives.begin()));

  preparedForDerivation = true;
}

expr_t
AbstractExternalFunctionNode::computeDerivative(int deriv_id)
{
  assert(datatree.external_functions_table.getNargs(symb_id) > 0);
  vector<expr_t> dargs;
  for (auto argument : arguments)
    dargs.push_back(argument->getDerivative(deriv_id));
  return composeDerivatives(dargs);
}

expr_t
AbstractExternalFunctionNode::getChainRuleDerivative(int deriv_id, const map<int, BinaryOpNode *> &recursive_variables)
{
  assert(datatree.external_functions_table.getNargs(symb_id) > 0);
  vector<expr_t> dargs;
  for (auto argument : arguments)
    dargs.push_back(argument->getChainRuleDerivative(deriv_id, recursive_variables));
  return composeDerivatives(dargs);
}

unsigned int
AbstractExternalFunctionNode::compileExternalFunctionArguments(ostream &CompileCode, unsigned int &instruction_number,
                                                               bool lhs_rhs, const temporary_terms_t &temporary_terms,
                                                               const temporary_terms_idxs_t &temporary_terms_idxs, bool dynamic, bool steady_dynamic,
                                                               const deriv_node_temp_terms_t &tef_terms) const
{
  for (auto argument : arguments)
    argument->compile(CompileCode, instruction_number, lhs_rhs, temporary_terms,
                      temporary_terms_idxs, dynamic, steady_dynamic, tef_terms);
  return (arguments.size());
}

void
AbstractExternalFunctionNode::collectVARLHSVariable(set<expr_t> &result) const
{
  cerr << "ERROR: you can only have variables or unary ops on LHS of VAR" << endl;
  exit(EXIT_FAILURE);
}

void
AbstractExternalFunctionNode::collectDynamicVariables(SymbolType type_arg, set<pair<int, int>> &result) const
{
  for (auto argument : arguments)
    argument->collectDynamicVariables(type_arg, result);
}

double
AbstractExternalFunctionNode::eval(const eval_context_t &eval_context) const noexcept(false)
{
  throw EvalExternalFunctionException();
}

int
AbstractExternalFunctionNode::maxEndoLead() const
{
  int val = 0;
  for (auto argument : arguments)
    val = max(val, argument->maxEndoLead());
  return val;
}

int
AbstractExternalFunctionNode::maxExoLead() const
{
  int val = 0;
  for (auto argument : arguments)
    val = max(val, argument->maxExoLead());
  return val;
}

int
AbstractExternalFunctionNode::maxEndoLag() const
{
  int val = 0;
  for (auto argument : arguments)
    val = max(val, argument->maxEndoLag());
  return val;
}

int
AbstractExternalFunctionNode::maxExoLag() const
{
  int val = 0;
  for (auto argument : arguments)
    val = max(val, argument->maxExoLag());
  return val;
}

int
AbstractExternalFunctionNode::maxLead() const
{
  int val = 0;
  for (auto argument : arguments)
    val = max(val, argument->maxLead());
  return val;
}

int
AbstractExternalFunctionNode::maxLag() const
{
  int val = 0;
  for (auto argument : arguments)
    val = max(val, argument->maxLag());
  return val;
}

int
AbstractExternalFunctionNode::maxLagWithDiffsExpanded() const
{
  int val = 0;
  for (auto argument : arguments)
    val = max(val, argument->maxLagWithDiffsExpanded());
  return val;
}

expr_t
AbstractExternalFunctionNode::undiff() const
{
  vector<expr_t> arguments_subst;
  for (auto argument : arguments)
    arguments_subst.push_back(argument->undiff());
  return buildSimilarExternalFunctionNode(arguments_subst, datatree);
}

int
AbstractExternalFunctionNode::VarMaxLag(const set<expr_t> &lhs_lag_equiv) const
{
  int max_lag = 0;
  for (auto argument : arguments)
    max_lag = max(max_lag, argument->VarMaxLag(lhs_lag_equiv));
  return max_lag;
}

expr_t
AbstractExternalFunctionNode::decreaseLeadsLags(int n) const
{
  vector<expr_t> arguments_subst;
  for (auto argument : arguments)
    arguments_subst.push_back(argument->decreaseLeadsLags(n));
  return buildSimilarExternalFunctionNode(arguments_subst, datatree);
}

expr_t
AbstractExternalFunctionNode::decreaseLeadsLagsPredeterminedVariables() const
{
  vector<expr_t> arguments_subst;
  for (auto argument : arguments)
    arguments_subst.push_back(argument->decreaseLeadsLagsPredeterminedVariables());
  return buildSimilarExternalFunctionNode(arguments_subst, datatree);
}

expr_t
AbstractExternalFunctionNode::substituteEndoLeadGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const
{
  vector<expr_t> arguments_subst;
  for (auto argument : arguments)
    arguments_subst.push_back(argument->substituteEndoLeadGreaterThanTwo(subst_table, neweqs, deterministic_model));
  return buildSimilarExternalFunctionNode(arguments_subst, datatree);
}

expr_t
AbstractExternalFunctionNode::substituteEndoLagGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  vector<expr_t> arguments_subst;
  for (auto argument : arguments)
    arguments_subst.push_back(argument->substituteEndoLagGreaterThanTwo(subst_table, neweqs));
  return buildSimilarExternalFunctionNode(arguments_subst, datatree);
}

expr_t
AbstractExternalFunctionNode::substituteExoLead(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const
{
  vector<expr_t> arguments_subst;
  for (auto argument : arguments)
    arguments_subst.push_back(argument->substituteExoLead(subst_table, neweqs, deterministic_model));
  return buildSimilarExternalFunctionNode(arguments_subst, datatree);
}

expr_t
AbstractExternalFunctionNode::substituteExoLag(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  vector<expr_t> arguments_subst;
  for (auto argument : arguments)
    arguments_subst.push_back(argument->substituteExoLag(subst_table, neweqs));
  return buildSimilarExternalFunctionNode(arguments_subst, datatree);
}

expr_t
AbstractExternalFunctionNode::substituteExpectation(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool partial_information_model) const
{
  vector<expr_t> arguments_subst;
  for (auto argument : arguments)
    arguments_subst.push_back(argument->substituteExpectation(subst_table, neweqs, partial_information_model));
  return buildSimilarExternalFunctionNode(arguments_subst, datatree);
}

expr_t
AbstractExternalFunctionNode::substituteAdl() const
{
  vector<expr_t> arguments_subst;
  for (auto argument : arguments)
    arguments_subst.push_back(argument->substituteAdl());
  return buildSimilarExternalFunctionNode(arguments_subst, datatree);
}

expr_t
AbstractExternalFunctionNode::substituteModelLocalVariables() const
{
  vector<expr_t> arguments_subst;
  for (auto argument : arguments)
    arguments_subst.push_back(argument->substituteModelLocalVariables());
  return buildSimilarExternalFunctionNode(arguments_subst, datatree);
}

expr_t
AbstractExternalFunctionNode::substituteVarExpectation(const map<string, expr_t> &subst_table) const
{
  vector<expr_t> arguments_subst;
  for (auto argument : arguments)
    arguments_subst.push_back(argument->substituteVarExpectation(subst_table));
  return buildSimilarExternalFunctionNode(arguments_subst, datatree);
}

void
AbstractExternalFunctionNode::findDiffNodes(lag_equivalence_table_t &nodes) const
{
  for (auto argument : arguments)
    argument->findDiffNodes(nodes);
}

void
AbstractExternalFunctionNode::findUnaryOpNodesForAuxVarCreation(lag_equivalence_table_t &nodes) const
{
  for (auto argument : arguments)
    argument->findUnaryOpNodesForAuxVarCreation(nodes);
}

int
AbstractExternalFunctionNode::findTargetVariable(int lhs_symb_id) const
{
  for (auto argument : arguments)
    if (int retval = argument->findTargetVariable(lhs_symb_id);
        retval >= 0)
      return retval;
  return -1;
}

expr_t
AbstractExternalFunctionNode::substituteDiff(const lag_equivalence_table_t &nodes, subst_table_t &subst_table,
                                             vector<BinaryOpNode *> &neweqs) const
{
  vector<expr_t> arguments_subst;
  for (auto argument : arguments)
    arguments_subst.push_back(argument->substituteDiff(nodes, subst_table, neweqs));
  return buildSimilarExternalFunctionNode(arguments_subst, datatree);
}

expr_t
AbstractExternalFunctionNode::substituteUnaryOpNodes(const lag_equivalence_table_t &nodes, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  vector<expr_t> arguments_subst;
  for (auto argument : arguments)
    arguments_subst.push_back(argument->substituteUnaryOpNodes(nodes, subst_table, neweqs));
  return buildSimilarExternalFunctionNode(arguments_subst, datatree);
}

int
AbstractExternalFunctionNode::countDiffs() const
{
  int ndiffs = 0;
  for (auto argument : arguments)
    ndiffs = max(ndiffs, argument->countDiffs());
  return ndiffs;
}

expr_t
AbstractExternalFunctionNode::substitutePacExpectation(const string &name, expr_t subexpr)
{
  vector<expr_t> arguments_subst;
  for (auto argument : arguments)
    arguments_subst.push_back(argument->substitutePacExpectation(name, subexpr));
  return buildSimilarExternalFunctionNode(arguments_subst, datatree);
}

expr_t
AbstractExternalFunctionNode::substitutePacTargetNonstationary(const string &name, expr_t subexpr)
{
  vector<expr_t> arguments_subst;
  for (auto argument : arguments)
    arguments_subst.push_back(argument->substitutePacTargetNonstationary(name, subexpr));
  return buildSimilarExternalFunctionNode(arguments_subst, datatree);
}

expr_t
AbstractExternalFunctionNode::differentiateForwardVars(const vector<string> &subset, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  vector<expr_t> arguments_subst;
  for (auto argument : arguments)
    arguments_subst.push_back(argument->differentiateForwardVars(subset, subst_table, neweqs));
  return buildSimilarExternalFunctionNode(arguments_subst, datatree);
}

bool
AbstractExternalFunctionNode::alreadyWrittenAsTefTerm(int the_symb_id, const deriv_node_temp_terms_t &tef_terms) const
{
  if (tef_terms.find({ the_symb_id, arguments }) != tef_terms.end())
    return true;
  return false;
}

int
AbstractExternalFunctionNode::getIndxInTefTerms(int the_symb_id, const deriv_node_temp_terms_t &tef_terms) const noexcept(false)
{
  if (auto it = tef_terms.find({ the_symb_id, arguments });
      it != tef_terms.end())
    return it->second;
  throw UnknownFunctionNameAndArgs();
}

void
AbstractExternalFunctionNode::computeTemporaryTerms(const pair<int, int> &derivOrder,
                                                    map<pair<int, int>, temporary_terms_t> &temp_terms_map,
                                                    map<expr_t, pair<int, pair<int, int>>> &reference_count,
                                                    bool is_matlab) const
{
  /* All external function nodes are declared as temporary terms.

     Given that temporary terms are separated in several functions (residuals,
     jacobian, …), we must make sure that all temporary terms derived from a
     given external function call are assigned just after that call.

     As a consequence, we need to “promote” some terms to a previous level (in
     the sense that residuals come before jacobian), if a temporary term
     corresponding to the same external function call is present in that
     previous level. */

  expr_t this2 = const_cast<AbstractExternalFunctionNode *>(this);
  for (auto &tt : temp_terms_map)
    if (auto it = find_if(tt.second.cbegin(), tt.second.cend(), sameTefTermPredicate());
        it != tt.second.cend())
      {
        tt.second.insert(this2);
        return;
      }

  temp_terms_map[derivOrder].insert(this2);
}

void
AbstractExternalFunctionNode::computeBlockTemporaryTerms(int blk, int eq, vector<vector<temporary_terms_t>> &blocks_temporary_terms,
                                                         map<expr_t, tuple<int, int, int>> &reference_count) const
{
  // See comments in computeTemporaryTerms() for the logic
  expr_t this2 = const_cast<AbstractExternalFunctionNode *>(this);
  for (auto &btt : blocks_temporary_terms)
    for (auto &tt : btt)
      if (auto it = find_if(tt.cbegin(), tt.cend(), sameTefTermPredicate());
          it != tt.cend())
        {
          tt.insert(this2);
          return;
        }

  blocks_temporary_terms[blk][eq].insert(this2);
}

bool
AbstractExternalFunctionNode::isNumConstNodeEqualTo(double value) const
{
  return false;
}

bool
AbstractExternalFunctionNode::isVariableNodeEqualTo(SymbolType type_arg, int variable_id, int lag_arg) const
{
  return false;
}

bool
AbstractExternalFunctionNode::containsPacExpectation(const string &pac_model_name) const
{
  for (auto argument : arguments)
    if (argument->containsPacExpectation(pac_model_name))
      return true;
  return false;
}

bool
AbstractExternalFunctionNode::containsPacTargetNonstationary(const string &pac_model_name) const
{
  for (auto argument : arguments)
    if (argument->containsPacTargetNonstationary(pac_model_name))
      return true;
  return false;
}

expr_t
AbstractExternalFunctionNode::replaceTrendVar() const
{
  vector<expr_t> arguments_subst;
  for (auto argument : arguments)
    arguments_subst.push_back(argument->replaceTrendVar());
  return buildSimilarExternalFunctionNode(arguments_subst, datatree);
}

expr_t
AbstractExternalFunctionNode::detrend(int symb_id, bool log_trend, expr_t trend) const
{
  vector<expr_t> arguments_subst;
  for (auto argument : arguments)
    arguments_subst.push_back(argument->detrend(symb_id, log_trend, trend));
  return buildSimilarExternalFunctionNode(arguments_subst, datatree);
}

expr_t
AbstractExternalFunctionNode::removeTrendLeadLag(const map<int, expr_t> &trend_symbols_map) const
{
  vector<expr_t> arguments_subst;
  for (auto argument : arguments)
    arguments_subst.push_back(argument->removeTrendLeadLag(trend_symbols_map));
  return buildSimilarExternalFunctionNode(arguments_subst, datatree);
}

bool
AbstractExternalFunctionNode::isInStaticForm() const
{
  for (auto argument : arguments)
    if (!argument->isInStaticForm())
      return false;
  return true;
}

bool
AbstractExternalFunctionNode::isParamTimesEndogExpr() const
{
  return false;
}

void
AbstractExternalFunctionNode::computeSubExprContainingVariable(int symb_id, int lag, set<expr_t> &contain_var) const
{
  bool var_present = false;
  for (auto arg : arguments)
    {
      arg->computeSubExprContainingVariable(symb_id, lag, contain_var);
      var_present = var_present || contain_var.count(arg) > 0;
    }
  if (var_present)
    contain_var.insert(const_cast<AbstractExternalFunctionNode *>(this));
}

BinaryOpNode *
AbstractExternalFunctionNode::normalizeEquationHelper(const set<expr_t> &contain_var, expr_t rhs) const
{
  throw NormalizationFailed();
}

void
AbstractExternalFunctionNode::writeExternalFunctionArguments(ostream &output, ExprNodeOutputType output_type,
                                                             const temporary_terms_t &temporary_terms,
                                                             const temporary_terms_idxs_t &temporary_terms_idxs,
                                                             const deriv_node_temp_terms_t &tef_terms) const
{
  for (auto it = arguments.begin(); it != arguments.end(); ++it)
    {
      if (it != arguments.begin())
        output << ",";

      (*it)->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
    }
}

void
AbstractExternalFunctionNode::writeJsonASTExternalFunctionArguments(ostream &output) const
{
  int i = 0;
  output << "{";
  for (auto it = arguments.begin(); it != arguments.end(); ++it, i++)
    {
      if (it != arguments.begin())
        output << ",";

      output << R"("arg)" << i << R"(" : )";
      (*it)->writeJsonAST(output);
    }
  output << "}";
}

void
AbstractExternalFunctionNode::writeJsonExternalFunctionArguments(ostream &output,
                                                                 const temporary_terms_t &temporary_terms,
                                                                 const deriv_node_temp_terms_t &tef_terms,
                                                                 bool isdynamic) const
{
  for (auto it = arguments.begin(); it != arguments.end(); ++it)
    {
      if (it != arguments.begin())
        output << ",";

      (*it)->writeJsonOutput(output, temporary_terms, tef_terms, isdynamic);
    }
}

void
AbstractExternalFunctionNode::writePrhs(ostream &output, ExprNodeOutputType output_type,
                                        const temporary_terms_t &temporary_terms,
                                        const temporary_terms_idxs_t &temporary_terms_idxs,
                                        const deriv_node_temp_terms_t &tef_terms) const
{
  int i = 0;
  for (auto argument : arguments)
    {
      output << "  prhs[" << i++ << "] = mxCreateDoubleScalar("; // All external_function arguments are scalars
      argument->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
      output << ");" << endl;
    }
}

bool
AbstractExternalFunctionNode::containsExternalFunction() const
{
  return true;
}

expr_t
AbstractExternalFunctionNode::substituteStaticAuxiliaryVariable() const
{
  vector<expr_t> arguments_subst;
  for (auto argument : arguments)
    arguments_subst.push_back(argument->substituteStaticAuxiliaryVariable());
  return buildSimilarExternalFunctionNode(arguments_subst, datatree);
}

expr_t
AbstractExternalFunctionNode::replaceVarsInEquation(map<VariableNode *, NumConstNode *> &table) const
{
  vector<expr_t> arguments_subst;
  for (auto argument : arguments)
    arguments_subst.push_back(argument->replaceVarsInEquation(table));
  return buildSimilarExternalFunctionNode(arguments_subst, datatree);
}

ExternalFunctionNode::ExternalFunctionNode(DataTree &datatree_arg,
                                           int idx_arg,
                                           int symb_id_arg,
                                           const vector<expr_t> &arguments_arg) :
  AbstractExternalFunctionNode{datatree_arg, idx_arg, symb_id_arg, arguments_arg}
{
}

expr_t
ExternalFunctionNode::composeDerivatives(const vector<expr_t> &dargs)
{
  vector<expr_t> dNodes;
  for (int i = 0; i < static_cast<int>(dargs.size()); i++)
    dNodes.push_back(datatree.AddTimes(dargs.at(i),
                                       datatree.AddFirstDerivExternalFunction(symb_id, arguments, i+1)));

  expr_t theDeriv = datatree.Zero;
  for (auto &dNode : dNodes)
    theDeriv = datatree.AddPlus(theDeriv, dNode);
  return theDeriv;
}

void
ExternalFunctionNode::compile(ostream &CompileCode, unsigned int &instruction_number,
                              bool lhs_rhs, const temporary_terms_t &temporary_terms,
                              const temporary_terms_idxs_t &temporary_terms_idxs, bool dynamic, bool steady_dynamic,
                              const deriv_node_temp_terms_t &tef_terms) const
{
  if (auto this2 = const_cast<ExternalFunctionNode *>(this);
      temporary_terms.find(this2) != temporary_terms.end())
    {
      if (dynamic)
        {
          FLDT_ fldt(temporary_terms_idxs.at(this2));
          fldt.write(CompileCode, instruction_number);
        }
      else
        {
          FLDST_ fldst(temporary_terms_idxs.at(this2));
          fldst.write(CompileCode, instruction_number);
        }
      return;
    }

  if (!lhs_rhs)
    {
      FLDTEF_ fldtef(getIndxInTefTerms(symb_id, tef_terms));
      fldtef.write(CompileCode, instruction_number);
    }
  else
    {
      FSTPTEF_ fstptef(getIndxInTefTerms(symb_id, tef_terms));
      fstptef.write(CompileCode, instruction_number);
    }
}

void
ExternalFunctionNode::compileExternalFunctionOutput(ostream &CompileCode, unsigned int &instruction_number,
                                                    bool lhs_rhs, const temporary_terms_t &temporary_terms,
                                                    const temporary_terms_idxs_t &temporary_terms_idxs, bool dynamic, bool steady_dynamic,
                                                    deriv_node_temp_terms_t &tef_terms) const
{
  int first_deriv_symb_id = datatree.external_functions_table.getFirstDerivSymbID(symb_id);
  assert(first_deriv_symb_id != ExternalFunctionsTable::IDSetButNoNameProvided);

  for (auto argument : arguments)
    argument->compileExternalFunctionOutput(CompileCode, instruction_number, lhs_rhs, temporary_terms,
                                            temporary_terms_idxs, dynamic, steady_dynamic, tef_terms);

  if (!alreadyWrittenAsTefTerm(symb_id, tef_terms))
    {
      tef_terms[{ symb_id, arguments }] = static_cast<int>(tef_terms.size());
      int indx = getIndxInTefTerms(symb_id, tef_terms);
      int second_deriv_symb_id = datatree.external_functions_table.getSecondDerivSymbID(symb_id);
      assert(second_deriv_symb_id != ExternalFunctionsTable::IDSetButNoNameProvided);

      unsigned int nb_output_arguments = 0;
      if (symb_id == first_deriv_symb_id
          && symb_id == second_deriv_symb_id)
        nb_output_arguments = 3;
      else if (symb_id == first_deriv_symb_id)
        nb_output_arguments = 2;
      else
        nb_output_arguments = 1;
      unsigned int nb_input_arguments = compileExternalFunctionArguments(CompileCode, instruction_number, lhs_rhs, temporary_terms,
                                                                         temporary_terms_idxs, dynamic, steady_dynamic, tef_terms);

      FCALL_ fcall(nb_output_arguments, nb_input_arguments, datatree.symbol_table.getName(symb_id), indx);
      switch (nb_output_arguments)
        {
        case 1:
          fcall.set_function_type(ExternalFunctionType::withoutDerivative);
          break;
        case 2:
          fcall.set_function_type(ExternalFunctionType::withFirstDerivative);
          break;
        case 3:
          fcall.set_function_type(ExternalFunctionType::withFirstAndSecondDerivative);
          break;
        }
      fcall.write(CompileCode, instruction_number);
      FSTPTEF_ fstptef(indx);
      fstptef.write(CompileCode, instruction_number);
    }
}

void
ExternalFunctionNode::writeJsonAST(ostream &output) const
{
  output << R"({"node_type" : "ExternalFunctionNode", )"
         << R"("name" : ")" << datatree.symbol_table.getName(symb_id) << R"(", "args" : [)";
  writeJsonASTExternalFunctionArguments(output);
  output << "]}";
}

void
ExternalFunctionNode::writeJsonOutput(ostream &output,
                                      const temporary_terms_t &temporary_terms,
                                      const deriv_node_temp_terms_t &tef_terms,
                                      bool isdynamic) const
{
  if (temporary_terms.find(const_cast<ExternalFunctionNode *>(this)) != temporary_terms.end())
    {
      output << "T" << idx;
      return;
    }

  try
    {
      int tef_idx = getIndxInTefTerms(symb_id, tef_terms);
      output << "TEF_" << tef_idx;
    }
  catch (UnknownFunctionNameAndArgs &)
    {
      // When writing the JSON output at parsing pass, we don’t use TEF terms
      output << datatree.symbol_table.getName(symb_id) << "(";
      writeJsonExternalFunctionArguments(output, temporary_terms, tef_terms, isdynamic);
      output << ")";
    }
}

void
ExternalFunctionNode::writeOutput(ostream &output, ExprNodeOutputType output_type,
                                  const temporary_terms_t &temporary_terms,
                                  const temporary_terms_idxs_t &temporary_terms_idxs,
                                  const deriv_node_temp_terms_t &tef_terms) const
{
  if (output_type == ExprNodeOutputType::matlabOutsideModel || output_type == ExprNodeOutputType::steadyStateFile
      || output_type == ExprNodeOutputType::juliaSteadyStateFile
      || output_type == ExprNodeOutputType::epilogueFile
      || output_type == ExprNodeOutputType::occbinDifferenceFile
      || isLatexOutput(output_type))
    {
      string name = isLatexOutput(output_type) ? datatree.symbol_table.getTeXName(symb_id)
        : datatree.symbol_table.getName(symb_id);
      output << name << "(";
      writeExternalFunctionArguments(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
      output << ")";
      return;
    }

  if (checkIfTemporaryTermThenWrite(output, output_type, temporary_terms, temporary_terms_idxs))
    return;

  if (isCOutput(output_type))
    output << "*";
  output << "TEF_" << getIndxInTefTerms(symb_id, tef_terms);
}

void
ExternalFunctionNode::writeExternalFunctionOutput(ostream &output, ExprNodeOutputType output_type,
                                                  const temporary_terms_t &temporary_terms,
                                                  const temporary_terms_idxs_t &temporary_terms_idxs,
                                                  deriv_node_temp_terms_t &tef_terms) const
{
  int first_deriv_symb_id = datatree.external_functions_table.getFirstDerivSymbID(symb_id);
  assert(first_deriv_symb_id != ExternalFunctionsTable::IDSetButNoNameProvided);

  for (auto argument : arguments)
    argument->writeExternalFunctionOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);

  if (!alreadyWrittenAsTefTerm(symb_id, tef_terms))
    {
      tef_terms[{ symb_id, arguments }] = static_cast<int>(tef_terms.size());
      int indx = getIndxInTefTerms(symb_id, tef_terms);
      int second_deriv_symb_id = datatree.external_functions_table.getSecondDerivSymbID(symb_id);
      assert(second_deriv_symb_id != ExternalFunctionsTable::IDSetButNoNameProvided);

      if (isCOutput(output_type))
        {
          output << "double *TEF_" << indx;
          if (symb_id == first_deriv_symb_id)
            output << ", *TEFD_" << indx;
          if (symb_id == second_deriv_symb_id)
            output << ", *TEFDD_" << indx;
          output << ";" << endl;

          if (symb_id == first_deriv_symb_id && symb_id == second_deriv_symb_id)
            output << "int TEFDD_" << indx << "_nrows;" << endl;

          int nlhs =
            symb_id == first_deriv_symb_id && symb_id == second_deriv_symb_id ? 3
            : symb_id == first_deriv_symb_id ? 2 : 1;
          output << "{" << endl
                 << "  mxArray *plhs[" << nlhs << "], *prhs[" << arguments.size() << "];" << endl;

          writePrhs(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);

          output << "  mexCallMATLAB(" << nlhs << ", plhs, " << arguments.size() << ", prhs, " << R"(")"
                 << datatree.symbol_table.getName(symb_id) << R"(");)" << endl;

          output << "  TEF_" << indx << " = mxGetPr(plhs[0]);" << endl;
          if (symb_id == first_deriv_symb_id)
            {
              output << "  TEFD_" << indx << " = mxGetPr(plhs[1]);" << endl;
              if (symb_id == second_deriv_symb_id)
                output << "  TEFDD_" << indx << " = mxGetPr(plhs[2]);" << endl
                       << "  TEFDD_" << indx << "_nrows = (int)mxGetM(plhs[2]);" << endl;
            }
          output << "}" << endl;
        }
      else
        {
          if (symb_id == first_deriv_symb_id
              && symb_id == second_deriv_symb_id)
            output << "[TEF_" << indx << ", TEFD_"<< indx << ", TEFDD_"<< indx << "] = ";
          else if (symb_id == first_deriv_symb_id)
            output << "[TEF_" << indx << ", TEFD_"<< indx << "] = ";
          else
            output << "TEF_" << indx << " = ";

          output << datatree.symbol_table.getName(symb_id) << "(";
          writeExternalFunctionArguments(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
          output << ");" << endl;
        }
    }
}

void
ExternalFunctionNode::writeJsonExternalFunctionOutput(vector<string> &efout,
                                                      const temporary_terms_t &temporary_terms,
                                                      deriv_node_temp_terms_t &tef_terms,
                                                      bool isdynamic) const
{
  int first_deriv_symb_id = datatree.external_functions_table.getFirstDerivSymbID(symb_id);
  assert(first_deriv_symb_id != ExternalFunctionsTable::IDSetButNoNameProvided);

  for (auto argument : arguments)
    argument->writeJsonExternalFunctionOutput(efout, temporary_terms, tef_terms, isdynamic);

  if (!alreadyWrittenAsTefTerm(symb_id, tef_terms))
    {
      tef_terms[{ symb_id, arguments }] = static_cast<int>(tef_terms.size());
      int indx = getIndxInTefTerms(symb_id, tef_terms);
      int second_deriv_symb_id = datatree.external_functions_table.getSecondDerivSymbID(symb_id);
      assert(second_deriv_symb_id != ExternalFunctionsTable::IDSetButNoNameProvided);

      stringstream ef;
      ef << R"({"external_function": {)"
         << R"("external_function_term": "TEF_)" << indx << R"(")";

      if (symb_id == first_deriv_symb_id)
        ef << R"(, "external_function_term_d": "TEFD_)" << indx << R"(")";

      if (symb_id == second_deriv_symb_id)
        ef << R"(, "external_function_term_dd": "TEFDD_)" << indx << R"(")";

      ef << R"(, "value": ")" << datatree.symbol_table.getName(symb_id) << "(";
      writeJsonExternalFunctionArguments(ef, temporary_terms, tef_terms, isdynamic);
      ef << R"lit()"}})lit";
      efout.push_back(ef.str());
    }
}

expr_t
ExternalFunctionNode::toStatic(DataTree &static_datatree) const
{
  vector<expr_t> static_arguments;
  for (auto argument : arguments)
    static_arguments.push_back(argument->toStatic(static_datatree));
  return static_datatree.AddExternalFunction(symb_id, static_arguments);
}

void
ExternalFunctionNode::computeXrefs(EquationInfo &ei) const
{
  vector<expr_t> dynamic_arguments;
  for (auto argument : arguments)
    argument->computeXrefs(ei);
}

expr_t
ExternalFunctionNode::clone(DataTree &datatree) const
{
  vector<expr_t> dynamic_arguments;
  for (auto argument : arguments)
    dynamic_arguments.push_back(argument->clone(datatree));
  return datatree.AddExternalFunction(symb_id, dynamic_arguments);
}

expr_t
ExternalFunctionNode::buildSimilarExternalFunctionNode(vector<expr_t> &alt_args, DataTree &alt_datatree) const
{
  return alt_datatree.AddExternalFunction(symb_id, alt_args);
}

function<bool (expr_t)>
ExternalFunctionNode::sameTefTermPredicate() const
{
  return [this](expr_t e) {
           auto e2 = dynamic_cast<ExternalFunctionNode *>(e);
           return (e2 != nullptr && e2->symb_id == symb_id);
         };
}

FirstDerivExternalFunctionNode::FirstDerivExternalFunctionNode(DataTree &datatree_arg,
                                                               int idx_arg,
                                                               int top_level_symb_id_arg,
                                                               const vector<expr_t> &arguments_arg,
                                                               int inputIndex_arg) :
  AbstractExternalFunctionNode{datatree_arg, idx_arg, top_level_symb_id_arg, arguments_arg},
  inputIndex{inputIndex_arg}
{
}

expr_t
FirstDerivExternalFunctionNode::composeDerivatives(const vector<expr_t> &dargs)
{
  vector<expr_t> dNodes;
  for (int i = 0; i < static_cast<int>(dargs.size()); i++)
    dNodes.push_back(datatree.AddTimes(dargs.at(i),
                                       datatree.AddSecondDerivExternalFunction(symb_id, arguments, inputIndex, i+1)));
  expr_t theDeriv = datatree.Zero;
  for (auto &dNode : dNodes)
    theDeriv = datatree.AddPlus(theDeriv, dNode);
  return theDeriv;
}

void
FirstDerivExternalFunctionNode::writeJsonAST(ostream &output) const
{
  output << R"({"node_type" : "FirstDerivExternalFunctionNode", )"
         << R"("name" : ")" << datatree.symbol_table.getName(symb_id) << R"(", "args" : [)";
  writeJsonASTExternalFunctionArguments(output);
  output << "]}";
}

void
FirstDerivExternalFunctionNode::writeJsonOutput(ostream &output,
                                                const temporary_terms_t &temporary_terms,
                                                const deriv_node_temp_terms_t &tef_terms,
                                                bool isdynamic) const
{
  // If current node is a temporary term
  if (temporary_terms.find(const_cast<FirstDerivExternalFunctionNode *>(this)) != temporary_terms.end())
    {
      output << "T" << idx;
      return;
    }

  const int first_deriv_symb_id = datatree.external_functions_table.getFirstDerivSymbID(symb_id);
  assert(first_deriv_symb_id != ExternalFunctionsTable::IDSetButNoNameProvided);

  const int tmpIndx = inputIndex - 1;

  if (first_deriv_symb_id == symb_id)
    output << "TEFD_" << getIndxInTefTerms(symb_id, tef_terms)
           << "[" << tmpIndx << "]";
  else if (first_deriv_symb_id == ExternalFunctionsTable::IDNotSet)
    output << "TEFD_fdd_" << getIndxInTefTerms(symb_id, tef_terms) << "_" << inputIndex;
  else
    output << "TEFD_def_" << getIndxInTefTerms(first_deriv_symb_id, tef_terms)
           << "[" << tmpIndx << "]";
}

void
FirstDerivExternalFunctionNode::writeOutput(ostream &output, ExprNodeOutputType output_type,
                                            const temporary_terms_t &temporary_terms,
                                            const temporary_terms_idxs_t &temporary_terms_idxs,
                                            const deriv_node_temp_terms_t &tef_terms) const
{
  assert(output_type != ExprNodeOutputType::matlabOutsideModel
         && output_type != ExprNodeOutputType::occbinDifferenceFile);

  if (isLatexOutput(output_type))
    {
      output << R"(\frac{\partial )" << datatree.symbol_table.getTeXName(symb_id)
             << R"(}{\partial )" << inputIndex << "}(";
      writeExternalFunctionArguments(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
      output << ")";
      return;
    }

  if (checkIfTemporaryTermThenWrite(output, output_type, temporary_terms, temporary_terms_idxs))
    return;

  const int first_deriv_symb_id = datatree.external_functions_table.getFirstDerivSymbID(symb_id);
  assert(first_deriv_symb_id != ExternalFunctionsTable::IDSetButNoNameProvided);

  const int tmpIndx = inputIndex - 1 + ARRAY_SUBSCRIPT_OFFSET(output_type);

  if (first_deriv_symb_id == symb_id)
    output << "TEFD_" << getIndxInTefTerms(symb_id, tef_terms)
           << LEFT_ARRAY_SUBSCRIPT(output_type) << tmpIndx << RIGHT_ARRAY_SUBSCRIPT(output_type);
  else if (first_deriv_symb_id == ExternalFunctionsTable::IDNotSet)
    {
      if (isCOutput(output_type))
        output << "*";
      output << "TEFD_fdd_" << getIndxInTefTerms(symb_id, tef_terms) << "_" << inputIndex;
    }
  else
    output << "TEFD_def_" << getIndxInTefTerms(first_deriv_symb_id, tef_terms)
           << LEFT_ARRAY_SUBSCRIPT(output_type) << tmpIndx << RIGHT_ARRAY_SUBSCRIPT(output_type);
}

void
FirstDerivExternalFunctionNode::compile(ostream &CompileCode, unsigned int &instruction_number,
                                        bool lhs_rhs, const temporary_terms_t &temporary_terms,
                                        const temporary_terms_idxs_t &temporary_terms_idxs, bool dynamic, bool steady_dynamic,
                                        const deriv_node_temp_terms_t &tef_terms) const
{
  if (auto this2 = const_cast<FirstDerivExternalFunctionNode *>(this);
      temporary_terms.find(this2) != temporary_terms.end())
    {
      if (dynamic)
        {
          FLDT_ fldt(temporary_terms_idxs.at(this2));
          fldt.write(CompileCode, instruction_number);
        }
      else
        {
          FLDST_ fldst(temporary_terms_idxs.at(this2));
          fldst.write(CompileCode, instruction_number);
        }
      return;
    }
  int first_deriv_symb_id = datatree.external_functions_table.getFirstDerivSymbID(symb_id);
  assert(first_deriv_symb_id != ExternalFunctionsTable::IDSetButNoNameProvided);

  if (!lhs_rhs)
    {
      FLDTEFD_ fldtefd(getIndxInTefTerms(symb_id, tef_terms), inputIndex);
      fldtefd.write(CompileCode, instruction_number);
    }
  else
    {
      FSTPTEFD_ fstptefd(getIndxInTefTerms(symb_id, tef_terms), inputIndex);
      fstptefd.write(CompileCode, instruction_number);
    }
}

void
FirstDerivExternalFunctionNode::writeExternalFunctionOutput(ostream &output, ExprNodeOutputType output_type,
                                                            const temporary_terms_t &temporary_terms,
                                                            const temporary_terms_idxs_t &temporary_terms_idxs,
                                                            deriv_node_temp_terms_t &tef_terms) const
{
  assert(output_type != ExprNodeOutputType::matlabOutsideModel);
  int first_deriv_symb_id = datatree.external_functions_table.getFirstDerivSymbID(symb_id);
  assert(first_deriv_symb_id != ExternalFunctionsTable::IDSetButNoNameProvided);

  /* For a node with derivs provided by the user function, call the method
     on the non-derived node */
  if (first_deriv_symb_id == symb_id)
    {
      expr_t parent = datatree.AddExternalFunction(symb_id, arguments);
      parent->writeExternalFunctionOutput(output, output_type, temporary_terms, temporary_terms_idxs,
                                          tef_terms);
      return;
    }

  if (alreadyWrittenAsTefTerm(first_deriv_symb_id, tef_terms))
    return;

  if (isCOutput(output_type))
    if (first_deriv_symb_id == ExternalFunctionsTable::IDNotSet)
      {
        output << "double *TEFD_fdd_" <<  getIndxInTefTerms(symb_id, tef_terms) << "_" << inputIndex << ";" << endl
               << "{" << endl
               << "  const mwSize dims[2] = {1, " << arguments.size() << "};" << endl
               << "  mxArray *plhs[1], *prhs[3];" << endl
               << R"(  prhs[0] = mxCreateString(")" << datatree.symbol_table.getName(symb_id) << R"(");)" << endl
               << "  prhs[1] = mxCreateDoubleScalar(" << inputIndex << ");"<< endl
               << "  prhs[2] = mxCreateCellArray(2, dims);"<< endl;

        int i = 0;
        for (auto argument : arguments)
          {
            output << "  mxSetCell(prhs[2], " << i++ << ", "
                   << "mxCreateDoubleScalar("; // All external_function arguments are scalars
            argument->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
            output << "));" << endl;
          }

        output << "  mexCallMATLAB(1, plhs, 3, prhs," << R"("jacob_element");)" << endl
               << "  TEFD_fdd_" <<  getIndxInTefTerms(symb_id, tef_terms) << "_" << inputIndex
               << " = mxGetPr(plhs[0]);" << endl
               << "}" << endl;
      }
    else
      {
        tef_terms[{ first_deriv_symb_id, arguments }] = static_cast<int>(tef_terms.size());
        int indx = getIndxInTefTerms(first_deriv_symb_id, tef_terms);
        output << "double *TEFD_def_" << indx << ";" << endl
               << "{" << endl
               << "  mxArray *plhs[1], *prhs[" << arguments.size() << "];" << endl;

        writePrhs(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);

        output << "  mexCallMATLAB(1, plhs, " << arguments.size() << ", prhs," << R"(")"
               << datatree.symbol_table.getName(first_deriv_symb_id) << R"(");)" << endl
               << "  TEFD_def_" << indx << " = mxGetPr(plhs[0]);" << endl
               << "}" << endl;
      }
  else
    {
      if (first_deriv_symb_id == ExternalFunctionsTable::IDNotSet)
        output << "TEFD_fdd_" << getIndxInTefTerms(symb_id, tef_terms) << "_" << inputIndex << " = jacob_element('"
               << datatree.symbol_table.getName(symb_id) << "'," << inputIndex << ",{";
      else
        {
          tef_terms[{ first_deriv_symb_id, arguments }] = static_cast<int>(tef_terms.size());
          output << "TEFD_def_" << getIndxInTefTerms(first_deriv_symb_id, tef_terms)
                 << " = " << datatree.symbol_table.getName(first_deriv_symb_id) << "(";
        }

      writeExternalFunctionArguments(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);

      if (first_deriv_symb_id == ExternalFunctionsTable::IDNotSet)
        output << "}";
      output << ");" << endl;
    }
}

void
FirstDerivExternalFunctionNode::writeJsonExternalFunctionOutput(vector<string> &efout,
                                                                const temporary_terms_t &temporary_terms,
                                                                deriv_node_temp_terms_t &tef_terms,
                                                                bool isdynamic) const
{
  int first_deriv_symb_id = datatree.external_functions_table.getFirstDerivSymbID(symb_id);
  assert(first_deriv_symb_id != ExternalFunctionsTable::IDSetButNoNameProvided);

  /* For a node with derivs provided by the user function, call the method
     on the non-derived node */
  if (first_deriv_symb_id == symb_id)
    {
      expr_t parent = datatree.AddExternalFunction(symb_id, arguments);
      parent->writeJsonExternalFunctionOutput(efout, temporary_terms, tef_terms, isdynamic);
      return;
    }

  if (alreadyWrittenAsTefTerm(first_deriv_symb_id, tef_terms))
    return;

  stringstream ef;
  if (first_deriv_symb_id == ExternalFunctionsTable::IDNotSet)
    ef << R"({"first_deriv_external_function": {)"
       << R"("external_function_term": "TEFD_fdd_)" << getIndxInTefTerms(symb_id, tef_terms) << "_" << inputIndex << R"(")"
       << R"(, "analytic_derivative": false)"
       << R"(, "wrt": )" << inputIndex
       << R"(, "value": ")" << datatree.symbol_table.getName(symb_id) << "(";
  else
    {
      tef_terms[{ first_deriv_symb_id, arguments }] = static_cast<int>(tef_terms.size());
      ef << R"({"first_deriv_external_function": {)"
         << R"("external_function_term": "TEFD_def_)" << getIndxInTefTerms(first_deriv_symb_id, tef_terms) << R"(")"
         << R"(, "analytic_derivative": true)"
         << R"(, "value": ")" << datatree.symbol_table.getName(first_deriv_symb_id) << "(";
    }

  writeJsonExternalFunctionArguments(ef, temporary_terms, tef_terms, isdynamic);
  ef << R"lit()"}})lit";
  efout.push_back(ef.str());
}

void
FirstDerivExternalFunctionNode::compileExternalFunctionOutput(ostream &CompileCode, unsigned int &instruction_number,
                                                              bool lhs_rhs, const temporary_terms_t &temporary_terms,
                                                              const temporary_terms_idxs_t &temporary_terms_idxs, bool dynamic, bool steady_dynamic,
                                                              deriv_node_temp_terms_t &tef_terms) const
{
  int first_deriv_symb_id = datatree.external_functions_table.getFirstDerivSymbID(symb_id);
  assert(first_deriv_symb_id != ExternalFunctionsTable::IDSetButNoNameProvided);

  if (first_deriv_symb_id == symb_id || alreadyWrittenAsTefTerm(first_deriv_symb_id, tef_terms))
    return;

  unsigned int nb_add_input_arguments = compileExternalFunctionArguments(CompileCode, instruction_number, lhs_rhs, temporary_terms,
                                                                         temporary_terms_idxs, dynamic, steady_dynamic, tef_terms);
  if (first_deriv_symb_id == ExternalFunctionsTable::IDNotSet)
    {
      unsigned int nb_input_arguments = 0;
      unsigned int nb_output_arguments = 1;
      int indx = getIndxInTefTerms(symb_id, tef_terms);
      FCALL_ fcall(nb_output_arguments, nb_input_arguments, "jacob_element", indx);
      fcall.set_arg_func_name(datatree.symbol_table.getName(symb_id));
      fcall.set_row(inputIndex);
      fcall.set_nb_add_input_arguments(nb_add_input_arguments);
      fcall.set_function_type(ExternalFunctionType::numericalFirstDerivative);
      fcall.write(CompileCode, instruction_number);
      FSTPTEFD_ fstptefd(indx, inputIndex);
      fstptefd.write(CompileCode, instruction_number);
    }
  else
    {
      tef_terms[{ first_deriv_symb_id, arguments }] = static_cast<int>(tef_terms.size());
      int indx = getIndxInTefTerms(symb_id, tef_terms);
      int second_deriv_symb_id = datatree.external_functions_table.getSecondDerivSymbID(symb_id);
      assert(second_deriv_symb_id != ExternalFunctionsTable::IDSetButNoNameProvided);

      unsigned int nb_output_arguments = 1;

      FCALL_ fcall(nb_output_arguments, nb_add_input_arguments, datatree.symbol_table.getName(first_deriv_symb_id), indx);
      fcall.set_function_type(ExternalFunctionType::firstDerivative);
      fcall.write(CompileCode, instruction_number);
      FSTPTEFD_ fstptefd(indx, inputIndex);
      fstptefd.write(CompileCode, instruction_number);
    }
}

expr_t
FirstDerivExternalFunctionNode::clone(DataTree &datatree) const
{
  vector<expr_t> dynamic_arguments;
  for (auto argument : arguments)
    dynamic_arguments.push_back(argument->clone(datatree));
  return datatree.AddFirstDerivExternalFunction(symb_id, dynamic_arguments,
                                                inputIndex);
}

expr_t
FirstDerivExternalFunctionNode::buildSimilarExternalFunctionNode(vector<expr_t> &alt_args, DataTree &alt_datatree) const
{
  return alt_datatree.AddFirstDerivExternalFunction(symb_id, alt_args, inputIndex);
}

expr_t
FirstDerivExternalFunctionNode::toStatic(DataTree &static_datatree) const
{
  vector<expr_t> static_arguments;
  for (auto argument : arguments)
    static_arguments.push_back(argument->toStatic(static_datatree));
  return static_datatree.AddFirstDerivExternalFunction(symb_id, static_arguments,
                                                       inputIndex);
}

void
FirstDerivExternalFunctionNode::computeXrefs(EquationInfo &ei) const
{
  vector<expr_t> dynamic_arguments;
  for (auto argument : arguments)
    argument->computeXrefs(ei);
}

function<bool (expr_t)>
FirstDerivExternalFunctionNode::sameTefTermPredicate() const
{
  int first_deriv_symb_id = datatree.external_functions_table.getFirstDerivSymbID(symb_id);
  if (first_deriv_symb_id == symb_id)
    return [this](expr_t e) {
      auto e2 = dynamic_cast<ExternalFunctionNode *>(e);
      return (e2 && e2->symb_id == symb_id);
    };
  else
    return [this](expr_t e) {
      auto e2 = dynamic_cast<FirstDerivExternalFunctionNode *>(e);
      return (e2 && e2->symb_id == symb_id);
    };
}

SecondDerivExternalFunctionNode::SecondDerivExternalFunctionNode(DataTree &datatree_arg,
                                                                 int idx_arg,
                                                                 int top_level_symb_id_arg,
                                                                 const vector<expr_t> &arguments_arg,
                                                                 int inputIndex1_arg,
                                                                 int inputIndex2_arg) :
  AbstractExternalFunctionNode{datatree_arg, idx_arg, top_level_symb_id_arg, arguments_arg},
  inputIndex1{inputIndex1_arg},
  inputIndex2{inputIndex2_arg}
{
}

expr_t
SecondDerivExternalFunctionNode::composeDerivatives(const vector<expr_t> &dargs)
{
  cerr << "ERROR: third order derivatives of external functions are not implemented" << endl;
  exit(EXIT_FAILURE);
}

void
SecondDerivExternalFunctionNode::writeJsonAST(ostream &output) const
{
  output << R"({"node_type" : "SecondDerivExternalFunctionNode", )"
         << R"("name" : ")" << datatree.symbol_table.getName(symb_id) << R"(", "args" : [)";
  writeJsonASTExternalFunctionArguments(output);
  output << "]}";
}

void
SecondDerivExternalFunctionNode::writeJsonOutput(ostream &output,
                                                 const temporary_terms_t &temporary_terms,
                                                 const deriv_node_temp_terms_t &tef_terms,
                                                 bool isdynamic) const
{
  // If current node is a temporary term
  if (temporary_terms.find(const_cast<SecondDerivExternalFunctionNode *>(this)) != temporary_terms.end())
    {
      output << "T" << idx;
      return;
    }

  const int second_deriv_symb_id = datatree.external_functions_table.getSecondDerivSymbID(symb_id);
  assert(second_deriv_symb_id != ExternalFunctionsTable::IDSetButNoNameProvided);

  const int tmpIndex1 = inputIndex1 - 1;
  const int tmpIndex2 = inputIndex2 - 1;

  if (second_deriv_symb_id == symb_id)
    output << "TEFDD_" << getIndxInTefTerms(symb_id, tef_terms)
           << "[" << tmpIndex1 << "," << tmpIndex2 << "]";
  else if (second_deriv_symb_id == ExternalFunctionsTable::IDNotSet)
    output << "TEFDD_fdd_" << getIndxInTefTerms(symb_id, tef_terms) << "_" << inputIndex1 << "_" << inputIndex2;
  else
    output << "TEFDD_def_" << getIndxInTefTerms(second_deriv_symb_id, tef_terms)
           << "[" << tmpIndex1 << "," << tmpIndex2 << "]";
}

void
SecondDerivExternalFunctionNode::writeOutput(ostream &output, ExprNodeOutputType output_type,
                                             const temporary_terms_t &temporary_terms,
                                             const temporary_terms_idxs_t &temporary_terms_idxs,
                                             const deriv_node_temp_terms_t &tef_terms) const
{
  assert(output_type != ExprNodeOutputType::matlabOutsideModel
         && output_type != ExprNodeOutputType::occbinDifferenceFile);

  if (isLatexOutput(output_type))
    {
      output << R"(\frac{\partial^2 )" << datatree.symbol_table.getTeXName(symb_id)
             << R"(}{\partial )" << inputIndex1 << R"(\partial )" << inputIndex2 << "}(";
      writeExternalFunctionArguments(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
      output << ")";
      return;
    }

  if (checkIfTemporaryTermThenWrite(output, output_type, temporary_terms, temporary_terms_idxs))
    return;

  const int second_deriv_symb_id = datatree.external_functions_table.getSecondDerivSymbID(symb_id);
  assert(second_deriv_symb_id != ExternalFunctionsTable::IDSetButNoNameProvided);

  const int tmpIndex1 = inputIndex1 - 1 + ARRAY_SUBSCRIPT_OFFSET(output_type);
  const int tmpIndex2 = inputIndex2 - 1 + ARRAY_SUBSCRIPT_OFFSET(output_type);

  int indx = getIndxInTefTerms(symb_id, tef_terms);
  if (second_deriv_symb_id == symb_id)
    if (isCOutput(output_type))
      output << "TEFDD_" << indx
             << LEFT_ARRAY_SUBSCRIPT(output_type) << tmpIndex1 << " * TEFDD_" << indx << "_nrows + "
             << tmpIndex2 << RIGHT_ARRAY_SUBSCRIPT(output_type);
    else
      output << "TEFDD_" << getIndxInTefTerms(symb_id, tef_terms)
             << LEFT_ARRAY_SUBSCRIPT(output_type) << tmpIndex1 << "," << tmpIndex2 << RIGHT_ARRAY_SUBSCRIPT(output_type);
  else if (second_deriv_symb_id == ExternalFunctionsTable::IDNotSet)
    {
      if (isCOutput(output_type))
        output << "*";
      output << "TEFDD_fdd_" << getIndxInTefTerms(symb_id, tef_terms) << "_" << inputIndex1 << "_" << inputIndex2;
    }
  else
    if (isCOutput(output_type))
      output << "TEFDD_def_" << getIndxInTefTerms(second_deriv_symb_id, tef_terms)
             << LEFT_ARRAY_SUBSCRIPT(output_type) << tmpIndex1 << " * PROBLEM_" << indx << "_nrows"
             << tmpIndex2 << RIGHT_ARRAY_SUBSCRIPT(output_type);
    else
      output << "TEFDD_def_" << getIndxInTefTerms(second_deriv_symb_id, tef_terms)
             << LEFT_ARRAY_SUBSCRIPT(output_type) << tmpIndex1 << "," << tmpIndex2 << RIGHT_ARRAY_SUBSCRIPT(output_type);
}

void
SecondDerivExternalFunctionNode::writeExternalFunctionOutput(ostream &output, ExprNodeOutputType output_type,
                                                             const temporary_terms_t &temporary_terms,
                                                             const temporary_terms_idxs_t &temporary_terms_idxs,
                                                             deriv_node_temp_terms_t &tef_terms) const
{
  assert(output_type != ExprNodeOutputType::matlabOutsideModel);
  int second_deriv_symb_id = datatree.external_functions_table.getSecondDerivSymbID(symb_id);
  assert(second_deriv_symb_id != ExternalFunctionsTable::IDSetButNoNameProvided);

  /* For a node with derivs provided by the user function, call the method
     on the non-derived node */
  if (second_deriv_symb_id == symb_id)
    {
      expr_t parent = datatree.AddExternalFunction(symb_id, arguments);
      parent->writeExternalFunctionOutput(output, output_type, temporary_terms, temporary_terms_idxs,
                                          tef_terms);
      return;
    }

  if (alreadyWrittenAsTefTerm(second_deriv_symb_id, tef_terms))
    return;

  if (isCOutput(output_type))
    if (second_deriv_symb_id == ExternalFunctionsTable::IDNotSet)
      {
        stringstream ending;
        output << "double *TEFDD_fdd_" <<  getIndxInTefTerms(symb_id, tef_terms) << "_" << inputIndex1 << "_" << inputIndex2 << ";" << endl
               << "{" << endl
               << "  const mwSize dims[2]= {1, " << arguments.size() << "};" << endl
               << "  mxArray *plhs[1], *prhs[4];" << endl
               << R"(  prhs[0] = mxCreateString(")" << datatree.symbol_table.getName(symb_id) << R"(");)" << endl
               << "  prhs[1] = mxCreateDoubleScalar(" << inputIndex1 << ");"<< endl
               << "  prhs[2] = mxCreateDoubleScalar(" << inputIndex2 << ");"<< endl
               << "  prhs[3] = mxCreateCellArray(2, dims);"<< endl;

        int i = 0;
        for (auto argument : arguments)
          {
            output << "  mxSetCell(prhs[3], " << i++ << ", "
                   << "  mxCreateDoubleScalar("; // All external_function arguments are scalars
            argument->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
            output << "));" << endl;
          }

        output << "  mexCallMATLAB(1, plhs, 4, prhs, " << R"("hess_element");)" << endl
               << "  TEFDD_fdd_" <<  getIndxInTefTerms(symb_id, tef_terms) << "_" << inputIndex1 << "_" << inputIndex2
               << " = mxGetPr(plhs[0]);" << endl
               << "}" << endl;
      }
    else
      {
        tef_terms[{ second_deriv_symb_id, arguments }] = static_cast<int>(tef_terms.size());
        int indx = getIndxInTefTerms(second_deriv_symb_id, tef_terms);
        output << "double *TEFDD_def_" << indx << ";" << endl
               << "{" << endl
               << "  mxArray *plhs[1], *prhs[" << arguments.size() << "];" << endl;

        writePrhs(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);

        output << "  mexCallMATLAB(1, plhs, " << arguments.size() << ", prhs, " << R"(")"
               << datatree.symbol_table.getName(second_deriv_symb_id) << R"(");)" << endl
               << "  TEFDD_def_" << indx << " = mxGetPr(plhs[0]);" << endl
               << "}" << endl;
      }
  else
    {
      if (second_deriv_symb_id == ExternalFunctionsTable::IDNotSet)
        output << "TEFDD_fdd_" << getIndxInTefTerms(symb_id, tef_terms) << "_" << inputIndex1 << "_" << inputIndex2
               << " = hess_element('" << datatree.symbol_table.getName(symb_id) << "',"
               << inputIndex1 << "," << inputIndex2 << ",{";
      else
        {
          tef_terms[{ second_deriv_symb_id, arguments }] = static_cast<int>(tef_terms.size());
          output << "TEFDD_def_" << getIndxInTefTerms(second_deriv_symb_id, tef_terms)
                 << " = " << datatree.symbol_table.getName(second_deriv_symb_id) << "(";
        }

      writeExternalFunctionArguments(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);

      if (second_deriv_symb_id == ExternalFunctionsTable::IDNotSet)
        output << "}";
      output << ");" << endl;
    }
}

void
SecondDerivExternalFunctionNode::writeJsonExternalFunctionOutput(vector<string> &efout,
                                                                 const temporary_terms_t &temporary_terms,
                                                                 deriv_node_temp_terms_t &tef_terms,
                                                                 bool isdynamic) const
{
  int second_deriv_symb_id = datatree.external_functions_table.getSecondDerivSymbID(symb_id);
  assert(second_deriv_symb_id != ExternalFunctionsTable::IDSetButNoNameProvided);

  /* For a node with derivs provided by the user function, call the method
     on the non-derived node */
  if (second_deriv_symb_id == symb_id)
    {
      expr_t parent = datatree.AddExternalFunction(symb_id, arguments);
      parent->writeJsonExternalFunctionOutput(efout, temporary_terms, tef_terms, isdynamic);
      return;
    }

  if (alreadyWrittenAsTefTerm(second_deriv_symb_id, tef_terms))
    return;

  stringstream ef;
  if (second_deriv_symb_id == ExternalFunctionsTable::IDNotSet)
    ef << R"({"second_deriv_external_function": {)"
       << R"("external_function_term": "TEFDD_fdd_)" << getIndxInTefTerms(symb_id, tef_terms) << "_" << inputIndex1 << "_" << inputIndex2 << R"(")"
       << R"(, "analytic_derivative": false)"
       << R"(, "wrt1": )" << inputIndex1
       << R"(, "wrt2": )" << inputIndex2
       << R"(, "value": ")" << datatree.symbol_table.getName(symb_id) << "(";
  else
    {
      tef_terms[{ second_deriv_symb_id, arguments }] = static_cast<int>(tef_terms.size());
      ef << R"({"second_deriv_external_function": {)"
         << R"("external_function_term": "TEFDD_def_)" << getIndxInTefTerms(second_deriv_symb_id, tef_terms) << R"(")"
         << R"(, "analytic_derivative": true)"
         << R"(, "value": ")" << datatree.symbol_table.getName(second_deriv_symb_id) << "(";
    }

  writeJsonExternalFunctionArguments(ef, temporary_terms, tef_terms, isdynamic);
  ef << R"lit()"}})lit" << endl;
  efout.push_back(ef.str());
}

expr_t
SecondDerivExternalFunctionNode::clone(DataTree &datatree) const
{
  vector<expr_t> dynamic_arguments;
  for (auto argument : arguments)
    dynamic_arguments.push_back(argument->clone(datatree));
  return datatree.AddSecondDerivExternalFunction(symb_id, dynamic_arguments,
                                                 inputIndex1, inputIndex2);
}

expr_t
SecondDerivExternalFunctionNode::buildSimilarExternalFunctionNode(vector<expr_t> &alt_args, DataTree &alt_datatree) const
{
  return alt_datatree.AddSecondDerivExternalFunction(symb_id, alt_args, inputIndex1, inputIndex2);
}

expr_t
SecondDerivExternalFunctionNode::toStatic(DataTree &static_datatree) const
{
  vector<expr_t> static_arguments;
  for (auto argument : arguments)
    static_arguments.push_back(argument->toStatic(static_datatree));
  return static_datatree.AddSecondDerivExternalFunction(symb_id, static_arguments,
                                                        inputIndex1, inputIndex2);
}

void
SecondDerivExternalFunctionNode::computeXrefs(EquationInfo &ei) const
{
  vector<expr_t> dynamic_arguments;
  for (auto argument : arguments)
    argument->computeXrefs(ei);
}

void
SecondDerivExternalFunctionNode::compile(ostream &CompileCode, unsigned int &instruction_number,
                                         bool lhs_rhs, const temporary_terms_t &temporary_terms,
                                         const temporary_terms_idxs_t &temporary_terms_idxs, bool dynamic, bool steady_dynamic,
                                         const deriv_node_temp_terms_t &tef_terms) const
{
  cerr << "SecondDerivExternalFunctionNode::compile: not implemented." << endl;
  exit(EXIT_FAILURE);
}

void
SecondDerivExternalFunctionNode::compileExternalFunctionOutput(ostream &CompileCode, unsigned int &instruction_number,
                                                               bool lhs_rhs, const temporary_terms_t &temporary_terms,
                                                               const temporary_terms_idxs_t &temporary_terms_idxs, bool dynamic, bool steady_dynamic,
                                                               deriv_node_temp_terms_t &tef_terms) const
{
  cerr << "SecondDerivExternalFunctionNode::compileExternalFunctionOutput: not implemented." << endl;
  exit(EXIT_FAILURE);
}

function<bool (expr_t)>
SecondDerivExternalFunctionNode::sameTefTermPredicate() const
{
  int second_deriv_symb_id = datatree.external_functions_table.getSecondDerivSymbID(symb_id);
  if (second_deriv_symb_id == symb_id)
    return [this](expr_t e) {
      auto e2 = dynamic_cast<ExternalFunctionNode *>(e);
      return (e2 && e2->symb_id == symb_id);
    };
  else
    return [this](expr_t e) {
      auto e2 = dynamic_cast<SecondDerivExternalFunctionNode *>(e);
      return (e2 && e2->symb_id == symb_id);
    };
}

SubModelNode::SubModelNode(DataTree &datatree_arg,
                           int idx_arg,
                           string model_name_arg) :
  ExprNode{datatree_arg, idx_arg},
  model_name{move(model_name_arg)}
{
}

void
SubModelNode::computeTemporaryTerms(const pair<int, int> &derivOrder,
                                          map<pair<int, int>, temporary_terms_t> &temp_terms_map,
                                          map<expr_t, pair<int, pair<int, int>>> &reference_count,
                                          bool is_matlab) const
{
  cerr << "SubModelNode::computeTemporaryTerms not implemented." << endl;
  exit(EXIT_FAILURE);
}

void
SubModelNode::computeBlockTemporaryTerms(int blk, int eq, vector<vector<temporary_terms_t>> &blocks_temporary_terms,
                                               map<expr_t, tuple<int, int, int>> &reference_count) const
{
  cerr << "SubModelNode::computeBlocksTemporaryTerms not implemented." << endl;
  exit(EXIT_FAILURE);
}

expr_t
SubModelNode::toStatic(DataTree &static_datatree) const
{
  cerr << "SubModelNode::toStatic not implemented." << endl;
  exit(EXIT_FAILURE);
}

void
SubModelNode::prepareForDerivation()
{
  cerr << "SubModelNode::prepareForDerivation not implemented." << endl;
  exit(EXIT_FAILURE);
}

expr_t
SubModelNode::computeDerivative(int deriv_id)
{
  cerr << "SubModelNode::computeDerivative not implemented." << endl;
  exit(EXIT_FAILURE);
}

expr_t
SubModelNode::getChainRuleDerivative(int deriv_id, const map<int, BinaryOpNode *> &recursive_variables)
{
  cerr << "SubModelNode::getChainRuleDerivative not implemented." << endl;
  exit(EXIT_FAILURE);
}

int
SubModelNode::maxEndoLead() const
{
  cerr << "SubModelNode::maxEndoLead not implemented." << endl;
  exit(EXIT_FAILURE);
}

int
SubModelNode::maxExoLead() const
{
  cerr << "SubModelNode::maxExoLead not implemented." << endl;
  exit(EXIT_FAILURE);
}

int
SubModelNode::maxEndoLag() const
{
  cerr << "SubModelNode::maxEndoLead not implemented." << endl;
  exit(EXIT_FAILURE);
}

int
SubModelNode::maxExoLag() const
{
  cerr << "SubModelNode::maxExoLead not implemented." << endl;
  exit(EXIT_FAILURE);
}

int
SubModelNode::maxLead() const
{
  cerr << "SubModelNode::maxLead not implemented." << endl;
  exit(EXIT_FAILURE);
}

int
SubModelNode::maxLag() const
{
  cerr << "SubModelNode::maxLag not implemented." << endl;
  exit(EXIT_FAILURE);
}

expr_t
SubModelNode::undiff() const
{
  cerr << "SubModelNode::undiff not implemented." << endl;
  exit(EXIT_FAILURE);
}

int
SubModelNode::VarMaxLag(const set<expr_t> &lhs_lag_equiv) const
{
  cerr << "SubModelNode::VarMaxLag not implemented." << endl;
  exit(EXIT_FAILURE);
}

expr_t
SubModelNode::decreaseLeadsLags(int n) const
{
  cerr << "SubModelNode::decreaseLeadsLags not implemented." << endl;
  exit(EXIT_FAILURE);
}

int
SubModelNode::countDiffs() const
{
  cerr << "SubModelNode::countDiffs not implemented." << endl;
  exit(EXIT_FAILURE);
}


expr_t
SubModelNode::substituteEndoLeadGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const
{
  cerr << "SubModelNode::substituteEndoLeadGreaterThanTwo not implemented." << endl;
  exit(EXIT_FAILURE);
}

expr_t
SubModelNode::substituteEndoLagGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  cerr << "SubModelNode::substituteEndoLagGreaterThanTwo not implemented." << endl;
  exit(EXIT_FAILURE);
}

expr_t
SubModelNode::substituteExoLead(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const
{
  cerr << "SubModelNode::substituteExoLead not implemented." << endl;
  exit(EXIT_FAILURE);
}

expr_t
SubModelNode::substituteExoLag(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  cerr << "SubModelNode::substituteExoLag not implemented." << endl;
  exit(EXIT_FAILURE);
}

bool
SubModelNode::containsExternalFunction() const
{
  return false;
}

double
SubModelNode::eval(const eval_context_t &eval_context) const noexcept(false)
{
  throw EvalException();
}

void
SubModelNode::computeXrefs(EquationInfo &ei) const
{
}

void
SubModelNode::collectVARLHSVariable(set<expr_t> &result) const
{
  cerr << "ERROR: you can only have variables or unary ops on LHS of VAR" << endl;
  exit(EXIT_FAILURE);
}

void
SubModelNode::collectDynamicVariables(SymbolType type_arg, set<pair<int, int>> &result) const
{
}

void
SubModelNode::compile(ostream &CompileCode, unsigned int &instruction_number,
                            bool lhs_rhs, const temporary_terms_t &temporary_terms,
                            const temporary_terms_idxs_t &temporary_terms_idxs, bool dynamic, bool steady_dynamic,
                            const deriv_node_temp_terms_t &tef_terms) const
{
  cerr << "SubModelNode::compile not implemented." << endl;
  exit(EXIT_FAILURE);
}

void
SubModelNode::computeSubExprContainingVariable(int symb_id, int lag, set<expr_t> &contain_var) const
{
}

BinaryOpNode *
SubModelNode::normalizeEquationHelper(const set<expr_t> &contain_var, expr_t rhs) const
{
  throw NormalizationFailed();
}

expr_t
SubModelNode::substituteExpectation(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool partial_information_model) const
{
  return const_cast<SubModelNode *>(this);
}

expr_t
SubModelNode::substituteAdl() const
{
  return const_cast<SubModelNode *>(this);
}

expr_t
SubModelNode::substituteModelLocalVariables() const
{
  return const_cast<SubModelNode *>(this);
}

void
SubModelNode::findDiffNodes(lag_equivalence_table_t &nodes) const
{
}

void
SubModelNode::findUnaryOpNodesForAuxVarCreation(lag_equivalence_table_t &nodes) const
{
}

int
SubModelNode::findTargetVariable(int lhs_symb_id) const
{
  return -1;
}

expr_t
SubModelNode::substituteDiff(const lag_equivalence_table_t &nodes, subst_table_t &subst_table,
                                   vector<BinaryOpNode *> &neweqs) const
{
  return const_cast<SubModelNode *>(this);
}

expr_t
SubModelNode::substituteUnaryOpNodes(const lag_equivalence_table_t &nodes, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  return const_cast<SubModelNode *>(this);
}

bool
SubModelNode::isNumConstNodeEqualTo(double value) const
{
  return false;
}

bool
SubModelNode::isVariableNodeEqualTo(SymbolType type_arg, int variable_id, int lag_arg) const
{
  return false;
}

bool
SubModelNode::isInStaticForm() const
{
  return false;
}

bool
SubModelNode::isParamTimesEndogExpr() const
{
  return false;
}

expr_t
SubModelNode::substituteStaticAuxiliaryVariable() const
{
  return const_cast<SubModelNode *>(this);
}

expr_t
SubModelNode::replaceVarsInEquation(map<VariableNode *, NumConstNode *> &table) const
{
  return const_cast<SubModelNode *>(this);
}

expr_t
SubModelNode::differentiateForwardVars(const vector<string> &subset, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  cerr << "SubModelNode::differentiateForwardVars not implemented." << endl;
  exit(EXIT_FAILURE);
}

expr_t
SubModelNode::decreaseLeadsLagsPredeterminedVariables() const
{
  cerr << "SubModelNode::decreaseLeadsLagsPredeterminedVariables not implemented." << endl;
  exit(EXIT_FAILURE);
}

expr_t
SubModelNode::replaceTrendVar() const
{
  cerr << "SubModelNode::replaceTrendVar not implemented." << endl;
  exit(EXIT_FAILURE);
}

expr_t
SubModelNode::detrend(int symb_id, bool log_trend, expr_t trend) const
{
  cerr << "SubModelNode::detrend not implemented." << endl;
  exit(EXIT_FAILURE);
}

expr_t
SubModelNode::removeTrendLeadLag(const map<int, expr_t> &trend_symbols_map) const
{
  cerr << "SubModelNode::removeTrendLeadLag not implemented." << endl;
  exit(EXIT_FAILURE);
}

VarExpectationNode::VarExpectationNode(DataTree &datatree_arg,
                                       int idx_arg,
                                       string model_name_arg) :
  SubModelNode{datatree_arg, idx_arg, move(model_name_arg)}
{
}

expr_t
VarExpectationNode::clone(DataTree &datatree) const
{
  return datatree.AddVarExpectation(model_name);
}


int
VarExpectationNode::maxLagWithDiffsExpanded() const
{
  /* This node will be substituted by lagged variables, so in theory we should
     return a strictly positive value. But from here this value is not easy to
     compute.
     We return 0, because currently this function is only called from
     DynamicModel::setLeadsLagsOrig(), and the maximum lag will nevertheless be
     correctly computed because the maximum lag of the VAR will be taken into
     account via the corresponding equations. */
  return 0;
}

void
VarExpectationNode::writeOutput(ostream &output, ExprNodeOutputType output_type,
                                const temporary_terms_t &temporary_terms,
                                const temporary_terms_idxs_t &temporary_terms_idxs,
                                const deriv_node_temp_terms_t &tef_terms) const
{
  assert(output_type != ExprNodeOutputType::matlabOutsideModel);

  if (isLatexOutput(output_type))
    {
      output << "VAR_EXPECTATION(" << model_name << ')';
      return;
    }

  cerr << "VarExpectationNode::writeOutput not implemented for non-LaTeX." << endl;
  exit(EXIT_FAILURE);
}

expr_t
VarExpectationNode::substituteVarExpectation(const map<string, expr_t> &subst_table) const
{
  auto it = subst_table.find(model_name);
  if (it == subst_table.end())
    {
      cerr << "ERROR: unknown model '" << model_name << "' used in var_expectation expression" << endl;
      exit(EXIT_FAILURE);
    }
  return it->second;
}

expr_t
VarExpectationNode::substitutePacExpectation(const string &name, expr_t subexpr)
{
  return const_cast<VarExpectationNode *>(this);
}

expr_t
VarExpectationNode::substitutePacTargetNonstationary(const string &name, expr_t subexpr)
{
  return const_cast<VarExpectationNode *>(this);
}

bool
VarExpectationNode::containsPacExpectation(const string &pac_model_name) const
{
  return false;
}

bool
VarExpectationNode::containsPacTargetNonstationary(const string &pac_model_name) const
{
  return false;
}

void
VarExpectationNode::writeJsonAST(ostream &output) const
{
  output << R"({"node_type" : "VarExpectationNode", )"
         << R"("name" : ")" << model_name << R"("})";
}

void
VarExpectationNode::writeJsonOutput(ostream &output,
                                    const temporary_terms_t &temporary_terms,
                                    const deriv_node_temp_terms_t &tef_terms,
                                    bool isdynamic) const
{
  output << "var_expectation("
         << "model_name = " << model_name
         << ")";
}

PacExpectationNode::PacExpectationNode(DataTree &datatree_arg,
                                       int idx_arg,
                                       string model_name_arg) :
  SubModelNode{datatree_arg, idx_arg, move(model_name_arg)}
{
}

expr_t
PacExpectationNode::clone(DataTree &datatree) const
{
  return datatree.AddPacExpectation(model_name);
}

void
PacExpectationNode::writeOutput(ostream &output, ExprNodeOutputType output_type,
                                const temporary_terms_t &temporary_terms,
                                const temporary_terms_idxs_t &temporary_terms_idxs,
                                const deriv_node_temp_terms_t &tef_terms) const
{
  assert(output_type != ExprNodeOutputType::matlabOutsideModel);
  if (isLatexOutput(output_type))
    {
      output << "PAC_EXPECTATION" << LEFT_PAR(output_type) << model_name << RIGHT_PAR(output_type);
      return;
    }
}

int
PacExpectationNode::maxLagWithDiffsExpanded() const
{
  // Same comment as in VarExpectationNode::maxLagWithDiffsExpanded()
  return 0;
}

expr_t
PacExpectationNode::substituteVarExpectation(const map<string, expr_t> &subst_table) const
{
  return const_cast<PacExpectationNode *>(this);
}

bool
PacExpectationNode::containsPacExpectation(const string &pac_model_name) const
{
  if (pac_model_name.empty())
    return true;
  else
    return pac_model_name == model_name;
}

bool
PacExpectationNode::containsPacTargetNonstationary(const string &pac_model_name) const
{
  return false;
}

void
PacExpectationNode::writeJsonAST(ostream &output) const
{
  output << R"({"node_type" : "PacExpectationNode", )"
         << R"("name" : ")" << model_name << R"("})";
}

void
PacExpectationNode::writeJsonOutput(ostream &output,
                                    const temporary_terms_t &temporary_terms,
                                    const deriv_node_temp_terms_t &tef_terms,
                                    bool isdynamic) const
{
  output << "pac_expectation("
         << "model_name = " << model_name
         << ")";
}

expr_t
PacExpectationNode::substitutePacExpectation(const string &name, expr_t subexpr)
{
  if (model_name != name)
    return const_cast<PacExpectationNode *>(this);
  return subexpr;
}

expr_t
PacExpectationNode::substitutePacTargetNonstationary(const string &name, expr_t subexpr)
{
  return const_cast<PacExpectationNode *>(this);
}

PacTargetNonstationaryNode::PacTargetNonstationaryNode(DataTree &datatree_arg,
                                                       int idx_arg,
                                                       string model_name_arg) :
  SubModelNode{datatree_arg, idx_arg, move(model_name_arg)}
{
}

expr_t
PacTargetNonstationaryNode::clone(DataTree &datatree) const
{
  return datatree.AddPacTargetNonstationary(model_name);
}

void
PacTargetNonstationaryNode::writeOutput(ostream &output, ExprNodeOutputType output_type,
                                        const temporary_terms_t &temporary_terms,
                                        const temporary_terms_idxs_t &temporary_terms_idxs,
                                        const deriv_node_temp_terms_t &tef_terms) const
{
  assert(output_type != ExprNodeOutputType::matlabOutsideModel);
  if (isLatexOutput(output_type))
    {
      output << "PAC_TARGET_NONSTATIONARY" << LEFT_PAR(output_type) << model_name << RIGHT_PAR(output_type);
      return;
    }
}

int
PacTargetNonstationaryNode::maxLagWithDiffsExpanded() const
{
  // This node will be replaced by the target lagged by one
  return 1;
}

expr_t
PacTargetNonstationaryNode::substituteVarExpectation(const map<string, expr_t> &subst_table) const
{
  return const_cast<PacTargetNonstationaryNode *>(this);
}

bool
PacTargetNonstationaryNode::containsPacExpectation(const string &pac_model_name) const
{
  return false;
}

bool
PacTargetNonstationaryNode::containsPacTargetNonstationary(const string &pac_model_name) const
{
  if (pac_model_name.empty())
    return true;
  else
    return pac_model_name == model_name;
}

void
PacTargetNonstationaryNode::writeJsonAST(ostream &output) const
{
  output << R"({"node_type" : "PacTargetNonstationaryNode", )"
         << R"("name" : ")" << model_name << R"("})";
}

void
PacTargetNonstationaryNode::writeJsonOutput(ostream &output,
                                            const temporary_terms_t &temporary_terms,
                                            const deriv_node_temp_terms_t &tef_terms,
                                            bool isdynamic) const
{
  output << "pac_target_nonstationary("
         << "model_name = " << model_name
         << ")";
}

expr_t
PacTargetNonstationaryNode::substitutePacExpectation(const string &name, expr_t subexpr)
{
  return const_cast<PacTargetNonstationaryNode *>(this);
}

expr_t
PacTargetNonstationaryNode::substitutePacTargetNonstationary(const string &name, expr_t subexpr)
{
  if (model_name != name)
    return const_cast<PacTargetNonstationaryNode *>(this);
  return subexpr;
}

void
ExprNode::decomposeAdditiveTerms(vector<pair<expr_t, int>> &terms, int current_sign) const
{
  terms.emplace_back(const_cast<ExprNode *>(this), current_sign);
}

void
UnaryOpNode::decomposeAdditiveTerms(vector<pair<expr_t, int>> &terms, int current_sign) const
{
  if (op_code == UnaryOpcode::uminus)
    arg->decomposeAdditiveTerms(terms, -current_sign);
  else
    ExprNode::decomposeAdditiveTerms(terms, current_sign);
}

void
BinaryOpNode::decomposeAdditiveTerms(vector<pair<expr_t, int>> &terms, int current_sign) const
{
  if (op_code == BinaryOpcode::plus || op_code == BinaryOpcode::minus)
    {
      arg1->decomposeAdditiveTerms(terms, current_sign);
      if (op_code == BinaryOpcode::plus)
        arg2->decomposeAdditiveTerms(terms, current_sign);
      else
        arg2->decomposeAdditiveTerms(terms, -current_sign);
    }
  else
    ExprNode::decomposeAdditiveTerms(terms, current_sign);
}

void
ExprNode::decomposeMultiplicativeFactors(vector<pair<expr_t, int>> &factors, int current_exponent) const
{
  factors.emplace_back(const_cast<ExprNode *>(this), current_exponent);
}

void
BinaryOpNode::decomposeMultiplicativeFactors(vector<pair<expr_t, int>> &factors, int current_exponent) const
{
  if (op_code == BinaryOpcode::times || op_code == BinaryOpcode::divide)
    {
      arg1->decomposeMultiplicativeFactors(factors, current_exponent);
      if (op_code == BinaryOpcode::times)
        arg2->decomposeMultiplicativeFactors(factors, current_exponent);
      else
        arg2->decomposeMultiplicativeFactors(factors, -current_exponent);
    }
  else
    ExprNode::decomposeMultiplicativeFactors(factors, current_exponent);
}

tuple<int, int, int, double>
ExprNode::matchVariableTimesConstantTimesParam(bool variable_obligatory) const
{
  int variable_id = -1, lag = 0, param_id = -1;
  double constant = 1.0;
  matchVTCTPHelper(variable_id, lag, param_id, constant, false);
  if (variable_obligatory && variable_id == -1)
    throw MatchFailureException{"No variable in this expression"};
  return {variable_id, lag, param_id, constant};
}

void
ExprNode::matchVTCTPHelper(int &var_id, int &lag, int &param_id, double &constant, bool at_denominator) const
{
  throw MatchFailureException{"Expression not allowed in linear combination of variables"};
}

void
NumConstNode::matchVTCTPHelper(int &var_id, int &lag, int &param_id, double &constant, bool at_denominator) const
{
  double myvalue = eval({});
  if (at_denominator)
    constant /= myvalue;
  else
    constant *= myvalue;
}

void
VariableNode::matchVTCTPHelper(int &var_id, int &lag, int &param_id, double &constant, bool at_denominator) const
{
  if (at_denominator)
    throw MatchFailureException{"A variable or parameter cannot appear at denominator"};

  SymbolType type = get_type();
  if (type == SymbolType::endogenous || type == SymbolType::exogenous)
    {
      if (var_id != -1)
        throw MatchFailureException{"More than one variable in this expression"};
      var_id = symb_id;
      lag = this->lag;
    }
  else if (type == SymbolType::parameter)
    {
      if (param_id != -1)
        throw MatchFailureException{"More than one parameter in this expression"};
      param_id = symb_id;
    }
  else
    throw MatchFailureException{"Symbol " + datatree.symbol_table.getName(symb_id) + " not allowed here"};
}

void
UnaryOpNode::matchVTCTPHelper(int &var_id, int &lag, int &param_id, double &constant, bool at_denominator) const
{
  if (op_code == UnaryOpcode::uminus)
    {
      constant = -constant;
      arg->matchVTCTPHelper(var_id, lag, param_id, constant, at_denominator);
    }
  else
    throw MatchFailureException{"Operator not allowed in this expression"};
}

void
BinaryOpNode::matchVTCTPHelper(int &var_id, int &lag, int &param_id, double &constant, bool at_denominator) const
{
  if (op_code == BinaryOpcode::times || op_code == BinaryOpcode::divide)
    {
      arg1->matchVTCTPHelper(var_id, lag, param_id, constant, at_denominator);
      if (op_code == BinaryOpcode::times)
        arg2->matchVTCTPHelper(var_id, lag, param_id, constant, at_denominator);
      else
        arg2->matchVTCTPHelper(var_id, lag, param_id, constant, !at_denominator);
    }
  else
    throw MatchFailureException{"Operator not allowed in this expression"};
}

vector<tuple<int, int, int, double>>
ExprNode::matchLinearCombinationOfVariables(bool variable_obligatory_in_each_term) const
{
  vector<pair<expr_t, int>> terms;
  decomposeAdditiveTerms(terms);

  vector<tuple<int, int, int, double>> result;

  for (auto [term, sign] : terms)
    {
      auto m = term->matchVariableTimesConstantTimesParam(variable_obligatory_in_each_term);
      get<3>(m) *= sign;
      result.push_back(m);
    }
  return result;
}

pair<int, vector<tuple<int, int, int, double>>>
ExprNode::matchParamTimesLinearCombinationOfVariables() const
{
  auto bopn = dynamic_cast<const BinaryOpNode *>(this);
  if (!bopn || bopn->op_code != BinaryOpcode::times)
    throw MatchFailureException{"Not a multiplicative expression"};

  expr_t param = bopn->arg1, lincomb = bopn->arg2;

  auto is_param = [](expr_t e) {
                    auto vn = dynamic_cast<VariableNode *>(e);
                    return vn && vn->get_type() == SymbolType::parameter;
                  };

  if (!is_param(param))
    {
      swap(param, lincomb);
      if (!is_param(param))
        throw MatchFailureException{"No parameter on either side of the multiplication"};
    }

  return { dynamic_cast<VariableNode *>(param)->symb_id, lincomb->matchLinearCombinationOfVariables() };
}

pair<int, int>
ExprNode::matchParamTimesTargetMinusVariable(int symb_id) const
{
  auto bopn = dynamic_cast<const BinaryOpNode *>(this);
  if (!bopn || bopn->op_code != BinaryOpcode::times)
    throw MatchFailureException{"Not a multiplicative expression"};

  expr_t param = bopn->arg1, minus = bopn->arg2;

  auto is_param = [](expr_t e) {
                    auto vn = dynamic_cast<VariableNode *>(e);
                    return vn && vn->get_type() == SymbolType::parameter;
                  };

  if (!is_param(param))
    {
      swap(param, minus);
      if (!is_param(param))
        throw MatchFailureException{"No parameter on either side of the multiplication"};
    }

  auto bminus = dynamic_cast<const BinaryOpNode *>(minus);
  if (!bminus || bminus->op_code != BinaryOpcode::minus)
    throw MatchFailureException{"Neither factor is a minus operator"};

  auto lhs_level = dynamic_cast<const VariableNode *>(bminus->arg2);
  auto target = dynamic_cast<const VariableNode *>(bminus->arg1);

  auto check_target = [&]()
    {
      if (target->get_type() != SymbolType::endogenous
          && target->get_type() != SymbolType::exogenous)
        return false;
      if (datatree.symbol_table.isAuxiliaryVariable(target->symb_id))
        {
          auto avi = datatree.symbol_table.getAuxVarInfo(target->symb_id);
          if (avi.get_type() == AuxVarType::pacTargetNonstationary && target->lag == -1)
            return true;
          return (avi.get_type() == AuxVarType::unaryOp
                  && avi.get_unary_op() == "log"
                  && avi.get_orig_symb_id() != -1
                  && !datatree.symbol_table.isAuxiliaryVariable(avi.get_orig_symb_id())
                  && target->lag + avi.get_orig_lead_lag() == -1);
        }
      else
        return target->lag == -1;
    };

  if (lhs_level && lhs_level->symb_id == symb_id && target && check_target())
    return { dynamic_cast<VariableNode *>(param)->symb_id, target->symb_id };
  else
    throw MatchFailureException{"Neither factor is of the form (target-variable) where target is endo or exo (possibly logged), and has one lag"};
}

pair<int, expr_t>
ExprNode::matchEndogenousTimesConstant() const
{
  throw MatchFailureException{"This expression is not of the form endogenous*constant"};
}

pair<int, expr_t>
VariableNode::matchEndogenousTimesConstant() const
{
  if (get_type() == SymbolType::endogenous)
    return { symb_id, datatree.One };
  else
    throw MatchFailureException{"This expression is not of the form endogenous*constant"};
}

pair<int, expr_t>
BinaryOpNode::matchEndogenousTimesConstant() const
{
  if (op_code == BinaryOpcode::times)
    {
      if (auto varg1 = dynamic_cast<VariableNode *>(arg1);
          varg1 && varg1->get_type() == SymbolType::endogenous && arg2->isConstant())
        return { varg1->symb_id, arg2 };
      if (auto varg2 = dynamic_cast<VariableNode *>(arg2);
          varg2 && varg2->get_type() == SymbolType::endogenous && arg1->isConstant())
        return { varg2->symb_id, arg1 };
    }
  throw MatchFailureException{"This expression is not of the form endogenous*constant"};
}

pair<vector<pair<int, expr_t>>, expr_t>
ExprNode::matchLinearCombinationOfEndogenousWithConstant() const
{
  vector<pair<expr_t, int>> all_terms;
  decomposeAdditiveTerms(all_terms);

  vector<pair<int, expr_t>> endo_terms;
  expr_t intercept = datatree.Zero;
  for (auto [term, sign] : all_terms)
    if (term->isConstant())
      {
        if (sign == -1)
          intercept = datatree.AddMinus(intercept, term);
        else
          intercept = datatree.AddPlus(intercept, term);
      }
    else
      {
        auto [endo_id, constant] = term->matchEndogenousTimesConstant();
        if (sign == -1)
          constant = datatree.AddUMinus(constant);
        endo_terms.emplace_back(endo_id, constant);
      }
  return { endo_terms, intercept };
}

string
ExprNode::toString() const
{
  ostringstream ss;
  writeJsonOutput(ss, {}, {});
  return ss.str();
}
