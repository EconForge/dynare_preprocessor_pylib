/*
 * Copyright (C) 2007-2018 Dynare Team
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

#include <iostream>
#include <iterator>
#include <algorithm>

#include <cassert>
#include <cmath>

#include <boost/bind.hpp>

#include "ExprNode.hh"
#include "DataTree.hh"
#include "ModFile.hh"

ExprNode::ExprNode(DataTree &datatree_arg) : datatree(datatree_arg), preparedForDerivation(false)
{
  // Add myself to datatree
  datatree.node_list.push_back(this);

  // Set my index and increment counter
  idx = datatree.node_counter++;
}

ExprNode::~ExprNode()
{
}

expr_t
ExprNode::getDerivative(int deriv_id)
{
  if (!preparedForDerivation)
    prepareForDerivation();

  // Return zero if derivative is necessarily null (using symbolic a priori)
  set<int>::const_iterator it = non_null_derivatives.find(deriv_id);
  if (it == non_null_derivatives.end())
    return datatree.Zero;

  // If derivative is stored in cache, use the cached value, otherwise compute it (and cache it)
  map<int, expr_t>::const_iterator it2 = derivatives.find(deriv_id);
  if (it2 != derivatives.end())
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
ExprNode::cost(const temporary_terms_t &temp_terms_map, bool is_matlab) const
{
  // For a terminal node, the cost is null
  return 0;
}

int
ExprNode::cost(const map<NodeTreeReference, temporary_terms_t> &temp_terms_map, bool is_matlab) const
{
  // For a terminal node, the cost is null
  return 0;
}

void
ExprNode::collectVariables(SymbolType type, set<int> &result) const
{
  set<pair<int, int> > symbs_lags;
  collectDynamicVariables(type, symbs_lags);
  transform(symbs_lags.begin(), symbs_lags.end(), inserter(result, result.begin()),
            boost::bind(&pair<int, int>::first, _1));
}

void
ExprNode::collectEndogenous(set<pair<int, int> > &result) const
{
  set<pair<int, int> > symb_ids;
  collectDynamicVariables(eEndogenous, symb_ids);
  for (set<pair<int, int> >::const_iterator it = symb_ids.begin();
       it != symb_ids.end(); it++)
    result.insert(make_pair(datatree.symbol_table.getTypeSpecificID(it->first), it->second));
}

void
ExprNode::collectExogenous(set<pair<int, int> > &result) const
{
  set<pair<int, int> > symb_ids;
  collectDynamicVariables(eExogenous, symb_ids);
  for (set<pair<int, int> >::const_iterator it = symb_ids.begin();
       it != symb_ids.end(); it++)
    result.insert(make_pair(datatree.symbol_table.getTypeSpecificID(it->first), it->second));
}

void
ExprNode::computeTemporaryTerms(map<expr_t, pair<int, NodeTreeReference> > &reference_count,
                                map<NodeTreeReference, temporary_terms_t> &temp_terms_map,
                                bool is_matlab, NodeTreeReference tr) const
{
  // Nothing to do for a terminal node
}

void
ExprNode::computeTemporaryTerms(map<expr_t, int> &reference_count,
                                temporary_terms_t &temporary_terms,
                                map<expr_t, pair<int, int> > &first_occurence,
                                int Curr_block,
                                vector<vector<temporary_terms_t> > &v_temporary_terms,
                                int equation) const
{
  // Nothing to do for a terminal node
}

pair<int, expr_t >
ExprNode::normalizeEquation(int var_endo, vector<pair<int, pair<expr_t, expr_t> > > &List_of_Op_RHS) const
{
  /* nothing to do */
  return (make_pair(0, (expr_t) NULL));
}

void
ExprNode::writeOutput(ostream &output) const
{
  writeOutput(output, oMatlabOutsideModel, temporary_terms_t());
}

void
ExprNode::writeOutput(ostream &output, ExprNodeOutputType output_type) const
{
  writeOutput(output, output_type, temporary_terms_t());
}

void
ExprNode::writeOutput(ostream &output, ExprNodeOutputType output_type, const temporary_terms_t &temporary_terms) const
{
  deriv_node_temp_terms_t tef_terms;
  writeOutput(output, output_type, temporary_terms, tef_terms);
}

void
ExprNode::compile(ostream &CompileCode, unsigned int &instruction_number,
                  bool lhs_rhs, const temporary_terms_t &temporary_terms,
                  const map_idx_t &map_idx, bool dynamic, bool steady_dynamic) const
{
  deriv_node_temp_terms_t tef_terms;
  compile(CompileCode, instruction_number, lhs_rhs, temporary_terms, map_idx, dynamic, steady_dynamic, tef_terms);
}

void
ExprNode::writeExternalFunctionOutput(ostream &output, ExprNodeOutputType output_type,
                                      const temporary_terms_t &temporary_terms,
                                      deriv_node_temp_terms_t &tef_terms) const
{
  // Nothing to do
}

void
ExprNode::writeJsonExternalFunctionOutput(vector<string> &efout,
                                          const temporary_terms_t &temporary_terms,
                                          deriv_node_temp_terms_t &tef_terms,
                                          const bool isdynamic) const
{
  // Nothing to do
}

void
ExprNode::compileExternalFunctionOutput(ostream &CompileCode, unsigned int &instruction_number,
                                        bool lhs_rhs, const temporary_terms_t &temporary_terms,
                                        const map_idx_t &map_idx, bool dynamic, bool steady_dynamic,
                                        deriv_node_temp_terms_t &tef_terms) const
{
  // Nothing to do
}

VariableNode *
ExprNode::createEndoLeadAuxiliaryVarForMyself(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  int n = maxEndoLead();
  assert(n >= 2);

  subst_table_t::const_iterator it = subst_table.find(this);
  if (it != subst_table.end())
    return const_cast<VariableNode *>(it->second);

  expr_t substexpr = decreaseLeadsLags(n-1);
  int lag = n-2;

  // Each iteration tries to create an auxvar such that auxvar(+1)=expr(-lag)
  // At the beginning (resp. end) of each iteration, substexpr is an expression (possibly an auxvar) equivalent to expr(-lag-1) (resp. expr(-lag))
  while (lag >= 0)
    {
      expr_t orig_expr = decreaseLeadsLags(lag);
      it = subst_table.find(orig_expr);
      if (it == subst_table.end())
        {
          int symb_id = datatree.symbol_table.addEndoLeadAuxiliaryVar(orig_expr->idx, substexpr);
          neweqs.push_back(dynamic_cast<BinaryOpNode *>(datatree.AddEqual(datatree.AddVariable(symb_id, 0), substexpr)));
          substexpr = datatree.AddVariable(symb_id, +1);
          assert(dynamic_cast<VariableNode *>(substexpr) != NULL);
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

  subst_table_t::const_iterator it = subst_table.find(this);
  if (it != subst_table.end())
    return const_cast<VariableNode *>(it->second);

  expr_t substexpr = decreaseLeadsLags(n);
  int lag = n-1;

  // Each iteration tries to create an auxvar such that auxvar(+1)=expr(-lag)
  // At the beginning (resp. end) of each iteration, substexpr is an expression (possibly an auxvar) equivalent to expr(-lag-1) (resp. expr(-lag))
  while (lag >= 0)
    {
      expr_t orig_expr = decreaseLeadsLags(lag);
      it = subst_table.find(orig_expr);
      if (it == subst_table.end())
        {
          int symb_id = datatree.symbol_table.addExoLeadAuxiliaryVar(orig_expr->idx, substexpr);
          neweqs.push_back(dynamic_cast<BinaryOpNode *>(datatree.AddEqual(datatree.AddVariable(symb_id, 0), substexpr)));
          substexpr = datatree.AddVariable(symb_id, +1);
          assert(dynamic_cast<VariableNode *>(substexpr) != NULL);
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

bool
ExprNode::isDiffPresent() const
{
  return false;
}

void
ExprNode::getEndosAndMaxLags(map<string, int> &model_endos_and_lags) const
{
}

NumConstNode::NumConstNode(DataTree &datatree_arg, int id_arg) :
  ExprNode(datatree_arg),
  id(id_arg)
{
  // Add myself to the num const map
  datatree.num_const_node_map[id] = this;
}

bool
NumConstNode::isDiffPresent() const
{
  return false;
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
NumConstNode::collectTemporary_terms(const temporary_terms_t &temporary_terms, temporary_terms_inuse_t &temporary_terms_inuse, int Curr_Block) const
{
  temporary_terms_t::const_iterator it = temporary_terms.find(const_cast<NumConstNode *>(this));
  if (it != temporary_terms.end())
    temporary_terms_inuse.insert(idx);
}

void
NumConstNode::writeOutput(ostream &output, ExprNodeOutputType output_type,
                          const temporary_terms_t &temporary_terms,
                          deriv_node_temp_terms_t &tef_terms) const
{
  temporary_terms_t::const_iterator it = temporary_terms.find(const_cast<NumConstNode *>(this));
  if (it != temporary_terms.end())
    if (output_type == oMatlabDynamicModelSparse)
      output << "T" << idx << "(it_)";
    else
      output << "T" << idx;
  else
    output << datatree.num_constants.get(id);
}

void
NumConstNode::writeJsonOutput(ostream &output,
                              const temporary_terms_t &temporary_terms,
                              deriv_node_temp_terms_t &tef_terms,
                              const bool isdynamic) const
{
  output << datatree.num_constants.get(id);
}

bool
NumConstNode::containsExternalFunction() const
{
  return false;
}

double
NumConstNode::eval(const eval_context_t &eval_context) const throw (EvalException, EvalExternalFunctionException)
{
  return (datatree.num_constants.getDouble(id));
}

void
NumConstNode::compile(ostream &CompileCode, unsigned int &instruction_number,
                      bool lhs_rhs, const temporary_terms_t &temporary_terms,
                      const map_idx_t &map_idx, bool dynamic, bool steady_dynamic,
                      deriv_node_temp_terms_t &tef_terms) const
{
  FLDC_ fldc(datatree.num_constants.getDouble(id));
  fldc.write(CompileCode, instruction_number);
}

void
NumConstNode::collectVARLHSVariable(set<expr_t> &result) const
{
}

void
NumConstNode::collectDynamicVariables(SymbolType type_arg, set<pair<int, int> > &result) const
{
}

pair<int, expr_t >
NumConstNode::normalizeEquation(int var_endo, vector<pair<int, pair<expr_t, expr_t> > > &List_of_Op_RHS) const
{
  /* return the numercial constant */
  return (make_pair(0, datatree.AddNonNegativeConstant(datatree.num_constants.get(id))));
}

expr_t
NumConstNode::getChainRuleDerivative(int deriv_id, const map<int, expr_t> &recursive_variables)
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
NumConstNode::cloneDynamic(DataTree &dynamic_datatree) const
{
  return dynamic_datatree.AddNonNegativeConstant(datatree.num_constants.get(id));
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
  return 0;
}

int
NumConstNode::maxLag() const
{
  return 0;
}

expr_t
NumConstNode::undiff() const
{
  return const_cast<NumConstNode *>(this);
}

int
NumConstNode::VarMinLag() const
{
  return 1;
}

void
NumConstNode::VarMaxLag(DataTree &static_datatree, set<expr_t> &static_lhs, int &max_lag) const
{
}

int
NumConstNode::PacMaxLag(vector<int> &lhs) const
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

void
NumConstNode::findDiffNodes(DataTree &static_datatree, diff_table_t &diff_table) const
{
}

void
NumConstNode::findDiffUnaryOpNodes(DataTree &static_datatree, set<UnaryOpNode *> &unary_op_nodes_in_diff) const
{
}

expr_t
NumConstNode::substituteDiff(DataTree &static_datatree, diff_table_t &diff_table, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  return const_cast<NumConstNode *>(this);
}

expr_t
NumConstNode::substituteDiffUnaryOpNodes(DataTree &static_datatree, set<UnaryOpNode *> &nodes, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  return const_cast<NumConstNode *>(this);
}

expr_t
NumConstNode::substitutePacExpectation(map<const PacExpectationNode *, const BinaryOpNode *> &subst_table)
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

void
NumConstNode::getEndosAndMaxLags(map<string, int> &model_endos_and_lags) const
{
}

bool
NumConstNode::containsEndogenous(void) const
{
  return false;
}

bool
NumConstNode::containsExogenous() const
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
NumConstNode::removeTrendLeadLag(map<int, expr_t> trend_symbols_map) const
{
  return const_cast<NumConstNode *>(this);
}

bool
NumConstNode::isInStaticForm() const
{
  return true;
}

void
NumConstNode::setVarExpectationIndex(map<string, pair<SymbolList, int> > &var_model_info)
{
}

void
NumConstNode::walkPacParameters(bool &pac_encountered, pair<int, int> &lhs, set<pair<int, pair<int, int> > > &ec_params_and_vars, set<pair<int, pair<int, int> > > &ar_params_and_vars) const
{
}

void
NumConstNode::addParamInfoToPac(pair<int, int> &lhs_arg, set<pair<int, pair<int, int> > > &ec_params_and_vars_arg, set<pair<int, pair<int, int> > > &ar_params_and_vars_arg)
{
}

void
NumConstNode::fillPacExpectationVarInfo(string &model_name_arg, vector<int> &lhs_arg, int max_lag_arg, vector<bool> &nonstationary_arg, int growth_symb_id_arg, int equation_number_arg)
{
}

bool
NumConstNode::isVarModelReferenced(const string &model_info_name) const
{
  return false;
}

expr_t
NumConstNode::substituteStaticAuxiliaryVariable() const
{
  return const_cast<NumConstNode *>(this);
}

VariableNode::VariableNode(DataTree &datatree_arg, int symb_id_arg, int lag_arg) :
  ExprNode(datatree_arg),
  symb_id(symb_id_arg),
  type(datatree.symbol_table.getType(symb_id_arg)),
  lag(lag_arg)
{
  // Add myself to the variable map
  datatree.variable_node_map[make_pair(symb_id, lag)] = this;

  // It makes sense to allow a lead/lag on parameters: during steady state calibration, endogenous and parameters can be swapped
  assert(type != eExternalFunction
         && (lag == 0 || (type != eModelLocalVariable && type != eModFileLocalVariable)));
}

void
VariableNode::prepareForDerivation()
{
  if (preparedForDerivation)
    return;

  preparedForDerivation = true;

  // Fill in non_null_derivatives
  switch (type)
    {
    case eEndogenous:
    case eExogenous:
    case eExogenousDet:
    case eParameter:
    case eTrend:
    case eLogTrend:
      // For a variable or a parameter, the only non-null derivative is with respect to itself
      non_null_derivatives.insert(datatree.getDerivID(symb_id, lag));
      break;
    case eModelLocalVariable:
      datatree.local_variables_table[symb_id]->prepareForDerivation();
      // Non null derivatives are those of the value of the local parameter
      non_null_derivatives = datatree.local_variables_table[symb_id]->non_null_derivatives;
      break;
    case eModFileLocalVariable:
    case eStatementDeclaredVariable:
    case eUnusedEndogenous:
      // Such a variable is never derived
      break;
    case eExternalFunction:
    case eEndogenousVAR:
      cerr << "VariableNode::prepareForDerivation: impossible case" << endl;
      exit(EXIT_FAILURE);
    }
}

expr_t
VariableNode::computeDerivative(int deriv_id)
{
  switch (type)
    {
    case eEndogenous:
    case eExogenous:
    case eExogenousDet:
    case eParameter:
    case eTrend:
    case eLogTrend:
      if (deriv_id == datatree.getDerivID(symb_id, lag))
        return datatree.One;
      else
        return datatree.Zero;
    case eModelLocalVariable:
      return datatree.local_variables_table[symb_id]->getDerivative(deriv_id);
    case eModFileLocalVariable:
      cerr << "ModFileLocalVariable is not derivable" << endl;
      exit(EXIT_FAILURE);
    case eStatementDeclaredVariable:
      cerr << "eStatementDeclaredVariable is not derivable" << endl;
      exit(EXIT_FAILURE);
    case eUnusedEndogenous:
      cerr << "eUnusedEndogenous is not derivable" << endl;
      exit(EXIT_FAILURE);
    case eExternalFunction:
    case eEndogenousVAR:
      cerr << "Impossible case!" << endl;
      exit(EXIT_FAILURE);
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

void
VariableNode::collectTemporary_terms(const temporary_terms_t &temporary_terms, temporary_terms_inuse_t &temporary_terms_inuse, int Curr_Block) const
{
  temporary_terms_t::const_iterator it = temporary_terms.find(const_cast<VariableNode *>(this));
  if (it != temporary_terms.end())
    temporary_terms_inuse.insert(idx);
  if (type == eModelLocalVariable)
    datatree.local_variables_table[symb_id]->collectTemporary_terms(temporary_terms, temporary_terms_inuse, Curr_Block);
}

bool
VariableNode::containsExternalFunction() const
{
  return false;
}

void
VariableNode::writeJsonOutput(ostream &output,
                              const temporary_terms_t &temporary_terms,
                              deriv_node_temp_terms_t &tef_terms,
                              const bool isdynamic) const
{
  temporary_terms_t::const_iterator it = temporary_terms.find(const_cast<VariableNode *>(this));
  if (it != temporary_terms.end())
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
                          deriv_node_temp_terms_t &tef_terms) const
{
  // If node is a temporary term
  temporary_terms_t::const_iterator it = temporary_terms.find(const_cast<VariableNode *>(this));
  if (it != temporary_terms.end())
    {
      if (output_type == oMatlabDynamicModelSparse)
        output << "T" << idx << "(it_)";
      else
        output << "T" << idx;
      return;
    }

  if (IS_LATEX(output_type))
    {
      if (output_type == oLatexDynamicSteadyStateOperator)
        output << "\\bar";
      output << "{" << datatree.symbol_table.getTeXName(symb_id);
      if (output_type == oLatexDynamicModel
          && (type == eEndogenous || type == eExogenous || type == eExogenousDet || type == eModelLocalVariable || type == eTrend || type == eLogTrend))
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
      output << "}";
      return;
    }

  int i;
  int tsid = datatree.symbol_table.getTypeSpecificID(symb_id);
  switch (type)
    {
    case eParameter:
      if (output_type == oMatlabOutsideModel)
        output << "M_.params" << "(" << tsid + 1 << ")";
      else
        output << "params" << LEFT_ARRAY_SUBSCRIPT(output_type) << tsid + ARRAY_SUBSCRIPT_OFFSET(output_type) << RIGHT_ARRAY_SUBSCRIPT(output_type);
      break;

    case eModelLocalVariable:
      if (output_type == oMatlabDynamicModelSparse || output_type == oMatlabStaticModelSparse
          || output_type == oMatlabDynamicSteadyStateOperator || output_type == oMatlabDynamicSparseSteadyStateOperator
          || output_type == oCDynamicSteadyStateOperator)
        {
          output << "(";
          datatree.local_variables_table[symb_id]->writeOutput(output, output_type, temporary_terms, tef_terms);
          output << ")";
        }
      else
        /* We append underscores to avoid name clashes with "g1" or "oo_" (see
           also ModelTree::writeModelLocalVariables) */
        output << datatree.symbol_table.getName(symb_id) << "__";
      break;

    case eModFileLocalVariable:
      output << datatree.symbol_table.getName(symb_id);
      break;

    case eEndogenous:
      switch (output_type)
        {
        case oJuliaDynamicModel:
        case oMatlabDynamicModel:
        case oCDynamicModel:
          i = datatree.getDynJacobianCol(datatree.getDerivID(symb_id, lag)) + ARRAY_SUBSCRIPT_OFFSET(output_type);
          output <<  "y" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << RIGHT_ARRAY_SUBSCRIPT(output_type);
          break;
        case oCDynamic2Model:
          i = tsid + (lag+1)*datatree.symbol_table.endo_nbr() + ARRAY_SUBSCRIPT_OFFSET(output_type);
          output <<  "y" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << RIGHT_ARRAY_SUBSCRIPT(output_type);
          break;
        case oCStaticModel:
        case oJuliaStaticModel:
        case oMatlabStaticModel:
        case oMatlabStaticModelSparse:
          i = tsid + ARRAY_SUBSCRIPT_OFFSET(output_type);
          output <<  "y" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << RIGHT_ARRAY_SUBSCRIPT(output_type);
          break;
        case oMatlabDynamicModelSparse:
          i = tsid + ARRAY_SUBSCRIPT_OFFSET(output_type);
          if (lag > 0)
            output << "y" << LEFT_ARRAY_SUBSCRIPT(output_type) << "it_+" << lag << ", " << i << RIGHT_ARRAY_SUBSCRIPT(output_type);
          else if (lag < 0)
            output << "y" << LEFT_ARRAY_SUBSCRIPT(output_type) << "it_" << lag << ", " << i << RIGHT_ARRAY_SUBSCRIPT(output_type);
          else
            output << "y" << LEFT_ARRAY_SUBSCRIPT(output_type) << "it_, " << i << RIGHT_ARRAY_SUBSCRIPT(output_type);
          break;
        case oMatlabOutsideModel:
          output << "oo_.steady_state(" << tsid + 1 << ")";
          break;
        case oJuliaDynamicSteadyStateOperator:
        case oMatlabDynamicSteadyStateOperator:
        case oMatlabDynamicSparseSteadyStateOperator:
          output << "steady_state" << LEFT_ARRAY_SUBSCRIPT(output_type) << tsid + 1 << RIGHT_ARRAY_SUBSCRIPT(output_type);
          break;
        case oCDynamicSteadyStateOperator:
          output << "steady_state[" << tsid << "]";
          break;
        case oJuliaSteadyStateFile:
        case oSteadyStateFile:
          output << "ys_" << LEFT_ARRAY_SUBSCRIPT(output_type) << tsid + 1 << RIGHT_ARRAY_SUBSCRIPT(output_type);
          break;
        case oCSteadyStateFile:
          output << "ys_[" << tsid << "]";
          break;
        case oMatlabDseries:
          output << "ds." << datatree.symbol_table.getName(symb_id);
          if (lag != 0)
            output << LEFT_ARRAY_SUBSCRIPT(output_type) << lag << RIGHT_ARRAY_SUBSCRIPT(output_type);
          break;
        default:
          cerr << "VariableNode::writeOutput: should not reach this point" << endl;
          exit(EXIT_FAILURE);
        }
      break;

    case eExogenous:
      i = tsid + ARRAY_SUBSCRIPT_OFFSET(output_type);
      switch (output_type)
        {
        case oJuliaDynamicModel:
        case oMatlabDynamicModel:
        case oMatlabDynamicModelSparse:
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
        case oCDynamicModel:
        case oCDynamic2Model:
          if (lag == 0)
            output <<  "x[it_+" << i << "*nb_row_x]";
          else if (lag > 0)
            output <<  "x[it_+" << lag << "+" << i << "*nb_row_x]";
          else
            output <<  "x[it_" << lag << "+" << i << "*nb_row_x]";
          break;
        case oCStaticModel:
        case oJuliaStaticModel:
        case oMatlabStaticModel:
        case oMatlabStaticModelSparse:
          output << "x" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << RIGHT_ARRAY_SUBSCRIPT(output_type);
          break;
        case oMatlabOutsideModel:
          assert(lag == 0);
          output <<  "oo_.exo_steady_state(" << i << ")";
          break;
        case oMatlabDynamicSteadyStateOperator:
          output <<  "oo_.exo_steady_state(" << i << ")";
          break;
        case oJuliaSteadyStateFile:
        case oSteadyStateFile:
          output << "exo_" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << RIGHT_ARRAY_SUBSCRIPT(output_type);
          break;
        case oCSteadyStateFile:
          output << "exo_[" << i - 1 << "]";
          break;
        case oMatlabDseries:
          output << "ds." << datatree.symbol_table.getName(symb_id);
          if (lag != 0)
            output << LEFT_ARRAY_SUBSCRIPT(output_type) << lag << RIGHT_ARRAY_SUBSCRIPT(output_type);
          break;
        default:
          cerr << "VariableNode::writeOutput: should not reach this point" << endl;
          exit(EXIT_FAILURE);
        }
      break;

    case eExogenousDet:
      i = tsid + datatree.symbol_table.exo_nbr() + ARRAY_SUBSCRIPT_OFFSET(output_type);
      switch (output_type)
        {
        case oJuliaDynamicModel:
        case oMatlabDynamicModel:
        case oMatlabDynamicModelSparse:
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
        case oCDynamicModel:
        case oCDynamic2Model:
          if (lag == 0)
            output <<  "x[it_+" << i << "*nb_row_x]";
          else if (lag > 0)
            output <<  "x[it_+" << lag << "+" << i << "*nb_row_x]";
          else
            output <<  "x[it_" << lag << "+" << i << "*nb_row_x]";
          break;
        case oCStaticModel:
        case oJuliaStaticModel:
        case oMatlabStaticModel:
        case oMatlabStaticModelSparse:
          output << "x" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << RIGHT_ARRAY_SUBSCRIPT(output_type);
          break;
        case oMatlabOutsideModel:
          assert(lag == 0);
          output <<  "oo_.exo_det_steady_state(" << tsid + 1 << ")";
          break;
        case oMatlabDynamicSteadyStateOperator:
          output <<  "oo_.exo_det_steady_state(" << tsid + 1 << ")";
          break;
        case oJuliaSteadyStateFile:
        case oSteadyStateFile:
          output << "exo_" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << RIGHT_ARRAY_SUBSCRIPT(output_type);
          break;
        case oCSteadyStateFile:
          output << "exo_[" << i - 1 << "]";
          break;
        case oMatlabDseries:
          output << "ds." << datatree.symbol_table.getName(symb_id);
          if (lag != 0)
            output << LEFT_ARRAY_SUBSCRIPT(output_type) << lag << RIGHT_ARRAY_SUBSCRIPT(output_type);
          break;
        default:
          cerr << "VariableNode::writeOutput: should not reach this point" << endl;
          exit(EXIT_FAILURE);
        }
      break;

    case eExternalFunction:
    case eTrend:
    case eLogTrend:
    case eStatementDeclaredVariable:
    case eUnusedEndogenous:
    case eEndogenousVAR:
      cerr << "Impossible case" << endl;
      exit(EXIT_FAILURE);
    }
}

expr_t
VariableNode::substituteStaticAuxiliaryVariable() const
{
  if (type == eEndogenous)
    {
      try
        {
          return datatree.symbol_table.getAuxiliaryVarsExprNode(symb_id)->substituteStaticAuxiliaryVariable();
        }
      catch (SymbolTable::SearchFailedException &e)
        {
        }
    }
  return const_cast<VariableNode *>(this);
}

double
VariableNode::eval(const eval_context_t &eval_context) const throw (EvalException, EvalExternalFunctionException)
{
  eval_context_t::const_iterator it = eval_context.find(symb_id);
  if (it == eval_context.end())
    throw EvalException();

  return it->second;
}

void
VariableNode::compile(ostream &CompileCode, unsigned int &instruction_number,
                      bool lhs_rhs, const temporary_terms_t &temporary_terms,
                      const map_idx_t &map_idx, bool dynamic, bool steady_dynamic,
                      deriv_node_temp_terms_t &tef_terms) const
{
  if (type == eModelLocalVariable || type == eModFileLocalVariable)
    datatree.local_variables_table[symb_id]->compile(CompileCode, instruction_number, lhs_rhs, temporary_terms, map_idx, dynamic, steady_dynamic, tef_terms);
  else
    {
      int tsid = datatree.symbol_table.getTypeSpecificID(symb_id);
      if (type == eExogenousDet)
        tsid += datatree.symbol_table.exo_nbr();
      if (!lhs_rhs)
        {
          if (dynamic)
            {
              if (steady_dynamic)  // steady state values in a dynamic model
                {
                  FLDVS_ fldvs(type, tsid);
                  fldvs.write(CompileCode, instruction_number);
                }
              else
                {
                  if (type == eParameter)
                    {
                      FLDV_ fldv(type, tsid);
                      fldv.write(CompileCode, instruction_number);
                    }
                  else
                    {
                      FLDV_ fldv(type, tsid, lag);
                      fldv.write(CompileCode, instruction_number);
                    }
                }
            }
          else
            {
              FLDSV_ fldsv(type, tsid);
              fldsv.write(CompileCode, instruction_number);
            }
        }
      else
        {
          if (dynamic)
            {
              if (steady_dynamic)  // steady state values in a dynamic model
                {
                  cerr << "Impossible case: steady_state in rhs of equation" << endl;
                  exit(EXIT_FAILURE);
                }
              else
                {
                  if (type == eParameter)
                    {
                      FSTPV_ fstpv(type, tsid);
                      fstpv.write(CompileCode, instruction_number);
                    }
                  else
                    {
                      FSTPV_ fstpv(type, tsid, lag);
                      fstpv.write(CompileCode, instruction_number);
                    }
                }
            }
          else
            {
              FSTPSV_ fstpsv(type, tsid);
              fstpsv.write(CompileCode, instruction_number);
            }
        }
    }
}

void
VariableNode::computeTemporaryTerms(map<expr_t, int> &reference_count,
                                    temporary_terms_t &temporary_terms,
                                    map<expr_t, pair<int, int> > &first_occurence,
                                    int Curr_block,
                                    vector<vector<temporary_terms_t> > &v_temporary_terms,
                                    int equation) const
{
  if (type == eModelLocalVariable)
    datatree.local_variables_table[symb_id]->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, Curr_block, v_temporary_terms, equation);
}

void
VariableNode::collectVARLHSVariable(set<expr_t> &result) const
{
  if (type == eEndogenous && lag == 0)
    result.insert(const_cast<VariableNode *>(this));
  else
    {
      cerr << "ERROR: A VAR must have one endogenous variable on the LHS." << endl;
      exit(EXIT_FAILURE);
    }
}

void
VariableNode::collectDynamicVariables(SymbolType type_arg, set<pair<int, int> > &result) const
{
  if (type == type_arg)
    result.insert(make_pair(symb_id, lag));
  if (type == eModelLocalVariable)
    datatree.local_variables_table[symb_id]->collectDynamicVariables(type_arg, result);
}

pair<int, expr_t>
VariableNode::normalizeEquation(int var_endo, vector<pair<int, pair<expr_t, expr_t> > > &List_of_Op_RHS) const
{
  /* The equation has to be normalized with respect to the current endogenous variable ascribed to it.
     The two input arguments are :
     - The ID of the endogenous variable associated to the equation.
     - The list of operators and operands needed to normalize the equation*

     The pair returned by NormalizeEquation is composed of
     - a flag indicating if the expression returned contains (flag = 1) or not (flag = 0)
     the endogenous variable related to the equation.
     If the expression contains more than one occurence of the associated endogenous variable,
     the flag is equal to 2.
     - an expression equal to the RHS if flag = 0 and equal to NULL elsewhere
  */
  if (type == eEndogenous)
    {
      if (datatree.symbol_table.getTypeSpecificID(symb_id) == var_endo && lag == 0)
        /* the endogenous variable */
        return (make_pair(1, (expr_t) NULL));
      else
        return (make_pair(0, datatree.AddVariableInternal(symb_id, lag)));
    }
  else
    {
      if (type == eParameter)
        return (make_pair(0, datatree.AddVariableInternal(symb_id, 0)));
      else
        return (make_pair(0, datatree.AddVariableInternal(symb_id, lag)));
    }
}

expr_t
VariableNode::getChainRuleDerivative(int deriv_id, const map<int, expr_t> &recursive_variables)
{
  switch (type)
    {
    case eEndogenous:
    case eExogenous:
    case eExogenousDet:
    case eParameter:
    case eTrend:
    case eLogTrend:
      if (deriv_id == datatree.getDerivID(symb_id, lag))
        return datatree.One;
      else
        {
          //if there is in the equation a recursive variable we could use a chaine rule derivation
          map<int, expr_t>::const_iterator it = recursive_variables.find(datatree.getDerivID(symb_id, lag));
          if (it != recursive_variables.end())
            {
              map<int, expr_t>::const_iterator it2 = derivatives.find(deriv_id);
              if (it2 != derivatives.end())
                return it2->second;
              else
                {
                  map<int, expr_t> recursive_vars2(recursive_variables);
                  recursive_vars2.erase(it->first);
                  //expr_t c = datatree.AddNonNegativeConstant("1");
                  expr_t d = datatree.AddUMinus(it->second->getChainRuleDerivative(deriv_id, recursive_vars2));
                  //d = datatree.AddTimes(c, d);
                  derivatives[deriv_id] = d;
                  return d;
                }
            }
          else
            return datatree.Zero;
        }
    case eModelLocalVariable:
      return datatree.local_variables_table[symb_id]->getChainRuleDerivative(deriv_id, recursive_variables);
    case eModFileLocalVariable:
      cerr << "ModFileLocalVariable is not derivable" << endl;
      exit(EXIT_FAILURE);
    case eStatementDeclaredVariable:
      cerr << "eStatementDeclaredVariable is not derivable" << endl;
      exit(EXIT_FAILURE);
    case eUnusedEndogenous:
      cerr << "eUnusedEndogenous is not derivable" << endl;
      exit(EXIT_FAILURE);
    case eExternalFunction:
    case eEndogenousVAR:
      cerr << "Impossible case!" << endl;
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
  switch (type)
    {
    case eEndogenous:
      ei.endo.insert(make_pair(symb_id, lag));
      break;
    case eExogenous:
      ei.exo.insert(make_pair(symb_id, lag));
      break;
    case eExogenousDet:
      ei.exo_det.insert(make_pair(symb_id, lag));
      break;
    case eParameter:
      ei.param.insert(make_pair(symb_id, 0));
      break;
    case eTrend:
    case eLogTrend:
    case eModelLocalVariable:
    case eModFileLocalVariable:
    case eStatementDeclaredVariable:
    case eUnusedEndogenous:
    case eExternalFunction:
    case eEndogenousVAR:
      break;
    }
}

expr_t
VariableNode::cloneDynamic(DataTree &dynamic_datatree) const
{
  return dynamic_datatree.AddVariable(symb_id, lag);
}

int
VariableNode::maxEndoLead() const
{
  switch (type)
    {
    case eEndogenous:
      return max(lag, 0);
    case eModelLocalVariable:
      return datatree.local_variables_table[symb_id]->maxEndoLead();
    default:
      return 0;
    }
}

int
VariableNode::maxExoLead() const
{
  switch (type)
    {
    case eExogenous:
      return max(lag, 0);
    case eModelLocalVariable:
      return datatree.local_variables_table[symb_id]->maxExoLead();
    default:
      return 0;
    }
}

int
VariableNode::maxEndoLag() const
{
  switch (type)
    {
    case eEndogenous:
      return max(-lag, 0);
    case eModelLocalVariable:
      return datatree.local_variables_table[symb_id]->maxEndoLag();
    default:
      return 0;
    }
}

int
VariableNode::maxExoLag() const
{
  switch (type)
    {
    case eExogenous:
      return max(-lag, 0);
    case eModelLocalVariable:
      return datatree.local_variables_table[symb_id]->maxExoLag();
    default:
      return 0;
    }
}

int
VariableNode::maxLead() const
{
  switch (type)
    {
    case eEndogenous:
      return lag;
    case eExogenous:
      return lag;
    case eModelLocalVariable:
      return datatree.local_variables_table[symb_id]->maxLead();
    default:
      return 0;
    }
}

int
VariableNode::VarMinLag() const
{
  switch (type)
    {
    case eEndogenous:
      return -lag;
    case eExogenous:
      if (lag > 0)
        return -lag;
      else
        return 1; // Can have contemporaneus exog in VAR
    case eModelLocalVariable:
      return datatree.local_variables_table[symb_id]->VarMinLag();
    default:
      return 1;
    }
}

int
VariableNode::maxLag() const
{
  switch (type)
    {
    case eEndogenous:
      return -lag;
    case eExogenous:
      return -lag;
    case eModelLocalVariable:
      return datatree.local_variables_table[symb_id]->maxLag();
    default:
      return 0;
    }
}

expr_t
VariableNode::undiff() const
{
  return const_cast<VariableNode *>(this);
}

void
VariableNode::VarMaxLag(DataTree &static_datatree, set<expr_t> &static_lhs, int &max_lag) const
{
  if (-lag > max_lag)
    max_lag = -lag;
}

int
VariableNode::PacMaxLag(vector<int> &lhs) const
{
  return -lag;
}

expr_t
VariableNode::substituteAdl() const
{
  return const_cast<VariableNode *>(this);
}

void
VariableNode::findDiffNodes(DataTree &static_datatree, diff_table_t &diff_table) const
{
}

void
VariableNode::findDiffUnaryOpNodes(DataTree &static_datatree, set<UnaryOpNode *> &unary_op_nodes_in_diff) const
{
}

expr_t
VariableNode::substituteDiff(DataTree &static_datatree, diff_table_t &diff_table, subst_table_t &subst_table,
                             vector<BinaryOpNode *> &neweqs) const
{
  return const_cast<VariableNode *>(this);
}

expr_t
VariableNode::substituteDiffUnaryOpNodes(DataTree &static_datatree, set<UnaryOpNode *> &nodes, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  return const_cast<VariableNode *>(this);
}

expr_t
VariableNode::substitutePacExpectation(map<const PacExpectationNode *, const BinaryOpNode *> &subst_table)
{
  return const_cast<VariableNode *>(this);
}

expr_t
VariableNode::decreaseLeadsLags(int n) const
{
  switch (type)
    {
    case eEndogenous:
    case eExogenous:
    case eExogenousDet:
    case eTrend:
    case eLogTrend:
      return datatree.AddVariable(symb_id, lag-n);
    case eModelLocalVariable:
      return datatree.local_variables_table[symb_id]->decreaseLeadsLags(n);
    default:
      return const_cast<VariableNode *>(this);
    }
}

expr_t
VariableNode::decreaseLeadsLagsPredeterminedVariables() const
{
  if (datatree.symbol_table.isPredetermined(symb_id))
    return decreaseLeadsLags(1);
  else
    return const_cast<VariableNode *>(this);
}

expr_t
VariableNode::substituteEndoLeadGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const
{
  expr_t value;
  switch (type)
    {
    case eEndogenous:
      if (lag <= 1)
        return const_cast<VariableNode *>(this);
      else
        return createEndoLeadAuxiliaryVarForMyself(subst_table, neweqs);
    case eModelLocalVariable:
      value = datatree.local_variables_table[symb_id];
      if (value->maxEndoLead() <= 1)
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
  expr_t value;
  subst_table_t::const_iterator it;
  int cur_lag;
  switch (type)
    {
    case eEndogenous:
      if (lag >= -1)
        return const_cast<VariableNode *>(this);

      it = subst_table.find(this);
      if (it != subst_table.end())
        return const_cast<VariableNode *>(it->second);

      substexpr = datatree.AddVariable(symb_id, -1);
      cur_lag = -2;

      // Each iteration tries to create an auxvar such that auxvar(-1)=curvar(cur_lag)
      // At the beginning (resp. end) of each iteration, substexpr is an expression (possibly an auxvar) equivalent to curvar(cur_lag+1) (resp. curvar(cur_lag))
      while (cur_lag >= lag)
        {
          VariableNode *orig_expr = datatree.AddVariable(symb_id, cur_lag);
          it = subst_table.find(orig_expr);
          if (it == subst_table.end())
            {
              int aux_symb_id = datatree.symbol_table.addEndoLagAuxiliaryVar(symb_id, cur_lag+1, substexpr);
              neweqs.push_back(dynamic_cast<BinaryOpNode *>(datatree.AddEqual(datatree.AddVariable(aux_symb_id, 0), substexpr)));
              substexpr = datatree.AddVariable(aux_symb_id, -1);
              subst_table[orig_expr] = substexpr;
            }
          else
            substexpr = const_cast<VariableNode *>(it->second);

          cur_lag--;
        }
      return substexpr;

    case eModelLocalVariable:
      value = datatree.local_variables_table[symb_id];
      if (value->maxEndoLag() <= 1)
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
  expr_t value;
  switch (type)
    {
    case eExogenous:
      if (lag <= 0)
        return const_cast<VariableNode *>(this);
      else
        return createExoLeadAuxiliaryVarForMyself(subst_table, neweqs);
    case eModelLocalVariable:
      value = datatree.local_variables_table[symb_id];
      if (value->maxExoLead() == 0)
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
  expr_t value;
  subst_table_t::const_iterator it;
  int cur_lag;
  switch (type)
    {
    case eExogenous:
      if (lag >= 0)
        return const_cast<VariableNode *>(this);

      it = subst_table.find(this);
      if (it != subst_table.end())
        return const_cast<VariableNode *>(it->second);

      substexpr = datatree.AddVariable(symb_id, 0);
      cur_lag = -1;

      // Each iteration tries to create an auxvar such that auxvar(-1)=curvar(cur_lag)
      // At the beginning (resp. end) of each iteration, substexpr is an expression (possibly an auxvar) equivalent to curvar(cur_lag+1) (resp. curvar(cur_lag))
      while (cur_lag >= lag)
        {
          VariableNode *orig_expr = datatree.AddVariable(symb_id, cur_lag);
          it = subst_table.find(orig_expr);
          if (it == subst_table.end())
            {
              int aux_symb_id = datatree.symbol_table.addExoLagAuxiliaryVar(symb_id, cur_lag+1, substexpr);
              neweqs.push_back(dynamic_cast<BinaryOpNode *>(datatree.AddEqual(datatree.AddVariable(aux_symb_id, 0), substexpr)));
              substexpr = datatree.AddVariable(aux_symb_id, -1);
              subst_table[orig_expr] = substexpr;
            }
          else
            substexpr = const_cast<VariableNode *>(it->second);

          cur_lag--;
        }
      return substexpr;

    case eModelLocalVariable:
      value = datatree.local_variables_table[symb_id];
      if (value->maxExoLag() == 0)
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
  return const_cast<VariableNode *>(this);
}

expr_t
VariableNode::differentiateForwardVars(const vector<string> &subset, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  expr_t value;
  switch (type)
    {
    case eEndogenous:
      assert(lag <= 1);
      if (lag <= 0
          || (subset.size() > 0
              && find(subset.begin(), subset.end(), datatree.symbol_table.getName(symb_id)) == subset.end()))
        return const_cast<VariableNode *>(this);
      else
        {
          subst_table_t::iterator it = subst_table.find(this);
          VariableNode *diffvar;
          if (it != subst_table.end())
            diffvar = const_cast<VariableNode *>(it->second);
          else
            {
              int aux_symb_id = datatree.symbol_table.addDiffForwardAuxiliaryVar(symb_id, datatree.AddMinus(datatree.AddVariable(symb_id, 0),
                                                                                                            datatree.AddVariable(symb_id, -1)));
              neweqs.push_back(dynamic_cast<BinaryOpNode *>(datatree.AddEqual(datatree.AddVariable(aux_symb_id, 0), datatree.AddMinus(datatree.AddVariable(symb_id, 0),
                                                                                                                                      datatree.AddVariable(symb_id, -1)))));
              diffvar = datatree.AddVariable(aux_symb_id, 1);
              subst_table[this] = diffvar;
            }
          return datatree.AddPlus(datatree.AddVariable(symb_id, 0), diffvar);
        }
    case eModelLocalVariable:
      value = datatree.local_variables_table[symb_id];
      if (value->maxEndoLead() <= 0)
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
  if (type == type_arg && datatree.symbol_table.getTypeSpecificID(symb_id) == variable_id && lag == lag_arg)
    return true;
  else
    return false;
}

bool
VariableNode::containsEndogenous(void) const
{
  if (type == eEndogenous)
    return true;
  else
    return false;
}

bool
VariableNode::containsExogenous() const
{
  return (type == eExogenous || type == eExogenousDet);
}

expr_t
VariableNode::replaceTrendVar() const
{
  if (get_type() == eTrend)
    return datatree.One;
  else if (get_type() == eLogTrend)
    return datatree.Zero;
  else
    return const_cast<VariableNode *>(this);
}

expr_t
VariableNode::detrend(int symb_id, bool log_trend, expr_t trend) const
{
  if (get_symb_id() != symb_id)
    return const_cast<VariableNode *>(this);

  if (log_trend)
    {
      if (get_lag() == 0)
        return datatree.AddPlus(const_cast<VariableNode *>(this), trend);
      else
        return datatree.AddPlus(const_cast<VariableNode *>(this), trend->decreaseLeadsLags(-1*get_lag()));
    }
  else
    {
      if (get_lag() == 0)
        return datatree.AddTimes(const_cast<VariableNode *>(this), trend);
      else
        return datatree.AddTimes(const_cast<VariableNode *>(this), trend->decreaseLeadsLags(-1*get_lag()));
    }
}

bool
VariableNode::isDiffPresent() const
{
  return false;
}

expr_t
VariableNode::removeTrendLeadLag(map<int, expr_t> trend_symbols_map) const
{
  if ((get_type() != eTrend && get_type() != eLogTrend) || get_lag() == 0)
    return const_cast<VariableNode *>(this);

  map<int, expr_t>::const_iterator it = trend_symbols_map.find(symb_id);
  expr_t noTrendLeadLagNode = new VariableNode(datatree, it->first, 0);
  bool log_trend = get_type() == eLogTrend;
  expr_t trend = it->second;

  if (get_lag() > 0)
    {
      expr_t growthFactorSequence = trend->decreaseLeadsLags(-1);
      if (log_trend)
        {
          for (int i = 1; i < get_lag(); i++)
            growthFactorSequence = datatree.AddPlus(growthFactorSequence, trend->decreaseLeadsLags(-1*(i+1)));
          return datatree.AddPlus(noTrendLeadLagNode, growthFactorSequence);
        }
      else
        {
          for (int i = 1; i < get_lag(); i++)
            growthFactorSequence = datatree.AddTimes(growthFactorSequence, trend->decreaseLeadsLags(-1*(i+1)));
          return datatree.AddTimes(noTrendLeadLagNode, growthFactorSequence);
        }
    }
  else //get_lag < 0
    {
      expr_t growthFactorSequence = trend;
      if (log_trend)
        {
          for (int i = 1; i < abs(get_lag()); i++)
            growthFactorSequence = datatree.AddPlus(growthFactorSequence, trend->decreaseLeadsLags(i));
          return datatree.AddMinus(noTrendLeadLagNode, growthFactorSequence);
        }
      else
        {
          for (int i = 1; i < abs(get_lag()); i++)
            growthFactorSequence = datatree.AddTimes(growthFactorSequence, trend->decreaseLeadsLags(i));
          return datatree.AddDivide(noTrendLeadLagNode, growthFactorSequence);
        }
    }
}

bool
VariableNode::isInStaticForm() const
{
  return lag == 0;
}

void
VariableNode::setVarExpectationIndex(map<string, pair<SymbolList, int> > &var_model_info)
{
}

void
VariableNode::walkPacParameters(bool &pac_encountered, pair<int, int> &lhs, set<pair<int, pair<int, int> > > &ec_params_and_vars, set<pair<int, pair<int, int> > > &ar_params_and_vars) const
{
}

void
VariableNode::addParamInfoToPac(pair<int, int> &lhs_arg, set<pair<int, pair<int, int> > > &ec_params_and_vars_arg, set<pair<int, pair<int, int> > > &ar_params_and_vars_arg)
{
}

void
VariableNode::fillPacExpectationVarInfo(string &model_name_arg, vector<int> &lhs_arg, int max_lag_arg, vector<bool> &nonstationary_arg, int growth_symb_id_arg, int equation_number_arg)
{
}

bool
VariableNode::isVarModelReferenced(const string &model_info_name) const
{
  return false;
}

void
VariableNode::getEndosAndMaxLags(map<string, int> &model_endos_and_lags) const
{
  string varname = datatree.symbol_table.getName(symb_id);
  if (type == eEndogenous)
    if (model_endos_and_lags.find(varname) == model_endos_and_lags.end())
      model_endos_and_lags[varname] = min(model_endos_and_lags[varname], lag);
    else
      model_endos_and_lags[varname] = lag;
}

UnaryOpNode::UnaryOpNode(DataTree &datatree_arg, UnaryOpcode op_code_arg, const expr_t arg_arg, int expectation_information_set_arg, int param1_symb_id_arg, int param2_symb_id_arg, const string &adl_param_name_arg, vector<int> adl_lags_arg) :
  ExprNode(datatree_arg),
  arg(arg_arg),
  expectation_information_set(expectation_information_set_arg),
  param1_symb_id(param1_symb_id_arg),
  param2_symb_id(param2_symb_id_arg),
  op_code(op_code_arg),
  adl_param_name(adl_param_name_arg),
  adl_lags(adl_lags_arg)
{
  // Add myself to the unary op map
  datatree.unary_op_node_map[make_pair(make_pair(arg, op_code),
                                       make_pair(make_pair(expectation_information_set,
                                                           make_pair(param1_symb_id, param2_symb_id)),
                                                 make_pair(adl_param_name, adl_lags)))] = this;
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
  if (op_code == oSteadyState || op_code == oSteadyStateParamDeriv
      || op_code == oSteadyStateParam2ndDeriv)
    datatree.addAllParamDerivId(non_null_derivatives);
}

expr_t
UnaryOpNode::composeDerivatives(expr_t darg, int deriv_id)
{
  expr_t t11, t12, t13, t14;

  switch (op_code)
    {
    case oUminus:
      return datatree.AddUMinus(darg);
    case oExp:
      return datatree.AddTimes(darg, this);
    case oLog:
      return datatree.AddDivide(darg, arg);
    case oLog10:
      t11 = datatree.AddExp(datatree.One);
      t12 = datatree.AddLog10(t11);
      t13 = datatree.AddDivide(darg, arg);
      return datatree.AddTimes(t12, t13);
    case oCos:
      t11 = datatree.AddSin(arg);
      t12 = datatree.AddUMinus(t11);
      return datatree.AddTimes(darg, t12);
    case oSin:
      t11 = datatree.AddCos(arg);
      return datatree.AddTimes(darg, t11);
    case oTan:
      t11 = datatree.AddTimes(this, this);
      t12 = datatree.AddPlus(t11, datatree.One);
      return datatree.AddTimes(darg, t12);
    case oAcos:
      t11 = datatree.AddSin(this);
      t12 = datatree.AddDivide(darg, t11);
      return datatree.AddUMinus(t12);
    case oAsin:
      t11 = datatree.AddCos(this);
      return datatree.AddDivide(darg, t11);
    case oAtan:
      t11 = datatree.AddTimes(arg, arg);
      t12 = datatree.AddPlus(datatree.One, t11);
      return datatree.AddDivide(darg, t12);
    case oCosh:
      t11 = datatree.AddSinh(arg);
      return datatree.AddTimes(darg, t11);
    case oSinh:
      t11 = datatree.AddCosh(arg);
      return datatree.AddTimes(darg, t11);
    case oTanh:
      t11 = datatree.AddTimes(this, this);
      t12 = datatree.AddMinus(datatree.One, t11);
      return datatree.AddTimes(darg, t12);
    case oAcosh:
      t11 = datatree.AddSinh(this);
      return datatree.AddDivide(darg, t11);
    case oAsinh:
      t11 = datatree.AddCosh(this);
      return datatree.AddDivide(darg, t11);
    case oAtanh:
      t11 = datatree.AddTimes(arg, arg);
      t12 = datatree.AddMinus(datatree.One, t11);
      return datatree.AddTimes(darg, t12);
    case oSqrt:
      t11 = datatree.AddPlus(this, this);
      return datatree.AddDivide(darg, t11);
    case oAbs:
      t11 = datatree.AddSign(arg);
      return datatree.AddTimes(t11, darg);
    case oSign:
      return datatree.Zero;
    case oSteadyState:
      if (datatree.isDynamic())
        {
          if (datatree.getTypeByDerivID(deriv_id) == eParameter)
            {
              VariableNode *varg = dynamic_cast<VariableNode *>(arg);
              if (varg == NULL)
                {
                  cerr << "UnaryOpNode::composeDerivatives: STEADY_STATE() should only be used on "
                       << "standalone variables (like STEADY_STATE(y)) to be derivable w.r.t. parameters" << endl;
                  exit(EXIT_FAILURE);
                }
              if (datatree.symbol_table.getType(varg->symb_id) == eEndogenous)
                return datatree.AddSteadyStateParamDeriv(arg, datatree.getSymbIDByDerivID(deriv_id));
              else
                return datatree.Zero;
            }
          else
            return datatree.Zero;
        }
      else
        return darg;
    case oSteadyStateParamDeriv:
      assert(datatree.isDynamic());
      if (datatree.getTypeByDerivID(deriv_id) == eParameter)
        {
          VariableNode *varg = dynamic_cast<VariableNode *>(arg);
          assert(varg != NULL);
          assert(datatree.symbol_table.getType(varg->symb_id) == eEndogenous);
          return datatree.AddSteadyStateParam2ndDeriv(arg, param1_symb_id, datatree.getSymbIDByDerivID(deriv_id));
        }
      else
        return datatree.Zero;
    case oSteadyStateParam2ndDeriv:
      assert(datatree.isDynamic());
      if (datatree.getTypeByDerivID(deriv_id) == eParameter)
        {
          cerr << "3rd derivative of STEADY_STATE node w.r.t. three parameters not implemented" << endl;
          exit(EXIT_FAILURE);
        }
      else
        return datatree.Zero;
    case oExpectation:
      cerr << "UnaryOpNode::composeDerivatives: not implemented on oExpectation" << endl;
      exit(EXIT_FAILURE);
    case oErf:
      // x^2
      t11 = datatree.AddPower(arg, datatree.Two);
      // exp(x^2)
      t12 =  datatree.AddExp(t11);
      // sqrt(pi)
      t11 = datatree.AddSqrt(datatree.Pi);
      // sqrt(pi)*exp(x^2)
      t13 = datatree.AddTimes(t11, t12);
      // 2/(sqrt(pi)*exp(x^2));
      t14 = datatree.AddDivide(datatree.Two, t13);
      // (2/(sqrt(pi)*exp(x^2)))*dx;
      return datatree.AddTimes(t14, darg);
    case oDiff:
      cerr << "UnaryOpNode::composeDerivatives: not implemented on oDiff" << endl;
      exit(EXIT_FAILURE);
    case oAdl:
      cerr << "UnaryOpNode::composeDerivatives: not implemented on oAdl" << endl;
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
UnaryOpNode::cost(const map<NodeTreeReference, temporary_terms_t> &temp_terms_map, bool is_matlab) const
{
  // For a temporary term, the cost is null
  for (map<NodeTreeReference, temporary_terms_t>::const_iterator it = temp_terms_map.begin();
       it != temp_terms_map.end(); it++)
    if (it->second.find(const_cast<UnaryOpNode *>(this)) != it->second.end())
      return 0;

  return cost(arg->cost(temp_terms_map, is_matlab), is_matlab);
}

int
UnaryOpNode::cost(const temporary_terms_t &temporary_terms, bool is_matlab) const
{
  // For a temporary term, the cost is null
  if (temporary_terms.find(const_cast<UnaryOpNode *>(this)) != temporary_terms.end())
    return 0;

  return cost(arg->cost(temporary_terms, is_matlab), is_matlab);
}

int
UnaryOpNode::cost(int cost, bool is_matlab) const
{
  if (is_matlab)
    // Cost for Matlab files
    switch (op_code)
      {
      case oUminus:
      case oSign:
        return cost + 70;
      case oExp:
        return cost + 160;
      case oLog:
        return cost + 300;
      case oLog10:
      case oErf:
        return cost + 16000;
      case oCos:
      case oSin:
      case oCosh:
        return cost + 210;
      case oTan:
        return cost + 230;
      case oAcos:
        return cost + 300;
      case oAsin:
        return cost + 310;
      case oAtan:
        return cost + 140;
      case oSinh:
        return cost + 240;
      case oTanh:
        return cost + 190;
      case oAcosh:
        return cost + 770;
      case oAsinh:
        return cost + 460;
      case oAtanh:
        return cost + 350;
      case oSqrt:
      case oAbs:
        return cost + 570;
      case oSteadyState:
      case oSteadyStateParamDeriv:
      case oSteadyStateParam2ndDeriv:
      case oExpectation:
        return cost;
      case oDiff:
        cerr << "UnaryOpNode::cost: not implemented on oDiff" << endl;
        exit(EXIT_FAILURE);
      case oAdl:
        cerr << "UnaryOpNode::cost: not implemented on oAdl" << endl;
        exit(EXIT_FAILURE);
      }
  else
    // Cost for C files
    switch (op_code)
      {
      case oUminus:
      case oSign:
        return cost + 3;
      case oExp:
      case oAcosh:
        return cost + 210;
      case oLog:
        return cost + 137;
      case oLog10:
        return cost + 139;
      case oCos:
      case oSin:
        return cost + 160;
      case oTan:
        return cost + 170;
      case oAcos:
      case oAtan:
        return cost + 190;
      case oAsin:
        return cost + 180;
      case oCosh:
      case oSinh:
      case oTanh:
      case oErf:
        return cost + 240;
      case oAsinh:
        return cost + 220;
      case oAtanh:
        return cost + 150;
      case oSqrt:
      case oAbs:
        return cost + 90;
      case oSteadyState:
      case oSteadyStateParamDeriv:
      case oSteadyStateParam2ndDeriv:
      case oExpectation:
        return cost;
      case oDiff:
        cerr << "UnaryOpNode::cost: not implemented on oDiff" << endl;
        exit(EXIT_FAILURE);
      case oAdl:
        cerr << "UnaryOpNode::cost: not implemented on oAdl" << endl;
        exit(EXIT_FAILURE);
      }
  exit(EXIT_FAILURE);
}

void
UnaryOpNode::computeTemporaryTerms(map<expr_t, pair<int, NodeTreeReference> > &reference_count,
                                   map<NodeTreeReference, temporary_terms_t> &temp_terms_map,
                                   bool is_matlab, NodeTreeReference tr) const
{
  expr_t this2 = const_cast<UnaryOpNode *>(this);

  map<expr_t, pair<int, NodeTreeReference > >::iterator it = reference_count.find(this2);
  if (it == reference_count.end())
    {
      reference_count[this2] = make_pair(1, tr);
      arg->computeTemporaryTerms(reference_count, temp_terms_map, is_matlab, tr);
    }
  else
    {
      reference_count[this2] = make_pair(it->second.first + 1, it->second.second);
      if (reference_count[this2].first * cost(temp_terms_map, is_matlab) > MIN_COST(is_matlab))
        temp_terms_map[reference_count[this2].second].insert(this2);
    }
}

void
UnaryOpNode::computeTemporaryTerms(map<expr_t, int> &reference_count,
                                   temporary_terms_t &temporary_terms,
                                   map<expr_t, pair<int, int> > &first_occurence,
                                   int Curr_block,
                                   vector< vector<temporary_terms_t> > &v_temporary_terms,
                                   int equation) const
{
  expr_t this2 = const_cast<UnaryOpNode *>(this);
  map<expr_t, int>::iterator it = reference_count.find(this2);
  if (it == reference_count.end())
    {
      reference_count[this2] = 1;
      first_occurence[this2] = make_pair(Curr_block, equation);
      arg->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, Curr_block, v_temporary_terms, equation);
    }
  else
    {
      reference_count[this2]++;
      if (reference_count[this2] * cost(temporary_terms, false) > MIN_COST_C)
        {
          temporary_terms.insert(this2);
          v_temporary_terms[first_occurence[this2].first][first_occurence[this2].second].insert(this2);
        }
    }
}

void
UnaryOpNode::collectTemporary_terms(const temporary_terms_t &temporary_terms, temporary_terms_inuse_t &temporary_terms_inuse, int Curr_Block) const
{
  temporary_terms_t::const_iterator it = temporary_terms.find(const_cast<UnaryOpNode *>(this));
  if (it != temporary_terms.end())
    temporary_terms_inuse.insert(idx);
  else
    arg->collectTemporary_terms(temporary_terms, temporary_terms_inuse, Curr_Block);
}

bool
UnaryOpNode::containsExternalFunction() const
{
  return arg->containsExternalFunction();
}

void
UnaryOpNode::writeJsonOutput(ostream &output,
                             const temporary_terms_t &temporary_terms,
                             deriv_node_temp_terms_t &tef_terms,
                              const bool isdynamic) const
{
  temporary_terms_t::const_iterator it = temporary_terms.find(const_cast<UnaryOpNode *>(this));
  if (it != temporary_terms.end())
    {
      output << "T" << idx;
      return;
    }

  // Always put parenthesis around uminus nodes
  if (op_code == oUminus)
    output << "(";

  switch (op_code)
    {
    case oUminus:
      output << "-";
      break;
    case oExp:
      output << "exp";
      break;
    case oLog:
      output << "log";
      break;
    case oLog10:
      output << "log10";
      break;
    case oCos:
      output << "cos";
      break;
    case oSin:
      output << "sin";
      break;
    case oTan:
      output << "tan";
      break;
    case oAcos:
      output << "acos";
      break;
    case oAsin:
      output << "asin";
      break;
    case oAtan:
      output << "atan";
      break;
    case oCosh:
      output << "cosh";
      break;
    case oSinh:
      output << "sinh";
      break;
    case oTanh:
      output << "tanh";
      break;
    case oAcosh:
      output << "acosh";
      break;
    case oAsinh:
      output << "asinh";
      break;
    case oAtanh:
      output << "atanh";
      break;
    case oSqrt:
      output << "sqrt";
      break;
    case oAbs:
      output << "abs";
      break;
    case oSign:
      output << "sign";
      break;
    case oDiff:
      output << "diff";
      break;
    case oAdl:
      output << "adl(";
      arg->writeJsonOutput(output, temporary_terms, tef_terms);
      output << ", '" << adl_param_name << "', [";
      for (vector<int>::const_iterator it = adl_lags.begin(); it != adl_lags.end(); it++)
        {
          if (it != adl_lags.begin())
            output << ", ";
          output << *it;
        }
      output << "])";
      return;
    case oSteadyState:
      output << "(";
      arg->writeJsonOutput(output, temporary_terms, tef_terms, isdynamic);
      output << ")";
      return;
    case oSteadyStateParamDeriv:
      {
        VariableNode *varg = dynamic_cast<VariableNode *>(arg);
        assert(varg != NULL);
        assert(datatree.symbol_table.getType(varg->symb_id) == eEndogenous);
        assert(datatree.symbol_table.getType(param1_symb_id) == eParameter);
        int tsid_endo = datatree.symbol_table.getTypeSpecificID(varg->symb_id);
        int tsid_param = datatree.symbol_table.getTypeSpecificID(param1_symb_id);
        output << "ss_param_deriv(" << tsid_endo+1 << "," << tsid_param+1 << ")";
      }
      return;
    case oSteadyStateParam2ndDeriv:
      {
        VariableNode *varg = dynamic_cast<VariableNode *>(arg);
        assert(varg != NULL);
        assert(datatree.symbol_table.getType(varg->symb_id) == eEndogenous);
        assert(datatree.symbol_table.getType(param1_symb_id) == eParameter);
        assert(datatree.symbol_table.getType(param2_symb_id) == eParameter);
        int tsid_endo = datatree.symbol_table.getTypeSpecificID(varg->symb_id);
        int tsid_param1 = datatree.symbol_table.getTypeSpecificID(param1_symb_id);
        int tsid_param2 = datatree.symbol_table.getTypeSpecificID(param2_symb_id);
        output << "ss_param_2nd_deriv(" << tsid_endo+1 << "," << tsid_param1+1
               << "," << tsid_param2+1 << ")";
      }
      return;
    case oExpectation:
      output << "EXPECTATION(" << expectation_information_set << ")";
      break;
    case oErf:
      output << "erf";
      break;
    }

  bool close_parenthesis = false;

  /* Enclose argument with parentheses if:
     - current opcode is not uminus, or
     - current opcode is uminus and argument has lowest precedence
  */
  if (op_code != oUminus
      || (op_code == oUminus
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
  if (op_code == oUminus)
    output << ")";
}

void
UnaryOpNode::writeOutput(ostream &output, ExprNodeOutputType output_type,
                         const temporary_terms_t &temporary_terms,
                         deriv_node_temp_terms_t &tef_terms) const
{
  // If node is a temporary term
  temporary_terms_t::const_iterator it = temporary_terms.find(const_cast<UnaryOpNode *>(this));
  if (it != temporary_terms.end())
    {
      if (output_type == oMatlabDynamicModelSparse)
        output << "T" << idx << "(it_)";
      else
        output << "T" << idx;
      return;
    }

  // Always put parenthesis around uminus nodes
  if (op_code == oUminus)
    output << LEFT_PAR(output_type);

  switch (op_code)
    {
    case oUminus:
      output << "-";
      break;
    case oExp:
      output << "exp";
      break;
    case oLog:
      if (IS_LATEX(output_type))
        output << "\\log";
      else
        output << "log";
      break;
    case oLog10:
      if (IS_LATEX(output_type))
        output << "\\log_{10}";
      else
        output << "log10";
      break;
    case oCos:
      output << "cos";
      break;
    case oSin:
      output << "sin";
      break;
    case oTan:
      output << "tan";
      break;
    case oAcos:
      output << "acos";
      break;
    case oAsin:
      output << "asin";
      break;
    case oAtan:
      output << "atan";
      break;
    case oCosh:
      output << "cosh";
      break;
    case oSinh:
      output << "sinh";
      break;
    case oTanh:
      output << "tanh";
      break;
    case oAcosh:
      output << "acosh";
      break;
    case oAsinh:
      output << "asinh";
      break;
    case oAtanh:
      output << "atanh";
      break;
    case oSqrt:
      output << "sqrt";
      break;
    case oAbs:
      output << "abs";
      break;
    case oSign:
      if (output_type == oCDynamicModel || output_type == oCStaticModel)
        output << "copysign";
      else
        output << "sign";
      break;
    case oSteadyState:
      ExprNodeOutputType new_output_type;
      switch (output_type)
        {
        case oMatlabDynamicModel:
          new_output_type = oMatlabDynamicSteadyStateOperator;
          break;
        case oLatexDynamicModel:
          new_output_type = oLatexDynamicSteadyStateOperator;
          break;
        case oCDynamicModel:
          new_output_type = oCDynamicSteadyStateOperator;
          break;
        case oJuliaDynamicModel:
          new_output_type = oJuliaDynamicSteadyStateOperator;
          break;
        case oMatlabDynamicModelSparse:
          new_output_type = oMatlabDynamicSparseSteadyStateOperator;
          break;
        default:
          new_output_type = output_type;
          break;
        }
      output << "(";
      arg->writeOutput(output, new_output_type, temporary_terms, tef_terms);
      output << ")";
      return;
    case oSteadyStateParamDeriv:
      {
        VariableNode *varg = dynamic_cast<VariableNode *>(arg);
        assert(varg != NULL);
        assert(datatree.symbol_table.getType(varg->symb_id) == eEndogenous);
        assert(datatree.symbol_table.getType(param1_symb_id) == eParameter);
        int tsid_endo = datatree.symbol_table.getTypeSpecificID(varg->symb_id);
        int tsid_param = datatree.symbol_table.getTypeSpecificID(param1_symb_id);
        assert(IS_MATLAB(output_type));
        output << "ss_param_deriv(" << tsid_endo+1 << "," << tsid_param+1 << ")";
      }
      return;
    case oSteadyStateParam2ndDeriv:
      {
        VariableNode *varg = dynamic_cast<VariableNode *>(arg);
        assert(varg != NULL);
        assert(datatree.symbol_table.getType(varg->symb_id) == eEndogenous);
        assert(datatree.symbol_table.getType(param1_symb_id) == eParameter);
        assert(datatree.symbol_table.getType(param2_symb_id) == eParameter);
        int tsid_endo = datatree.symbol_table.getTypeSpecificID(varg->symb_id);
        int tsid_param1 = datatree.symbol_table.getTypeSpecificID(param1_symb_id);
        int tsid_param2 = datatree.symbol_table.getTypeSpecificID(param2_symb_id);
        assert(IS_MATLAB(output_type));
        output << "ss_param_2nd_deriv(" << tsid_endo+1 << "," << tsid_param1+1
               << "," << tsid_param2+1 << ")";
      }
      return;
    case oExpectation:
      if (!IS_LATEX(output_type))
        {
          cerr << "UnaryOpNode::writeOutput: not implemented on oExpectation" << endl;
          exit(EXIT_FAILURE);
        }
      output << "\\mathbb{E}_{t";
      if (expectation_information_set != 0)
        {
          if (expectation_information_set > 0)
            output << "+";
          output << expectation_information_set;
        }
      output << "}";
      break;
    case oErf:
      output << "erf";
      break;
    case oDiff:
      output << "diff";
      break;
    case oAdl:
      output << "adl";
      break;
    }

  bool close_parenthesis = false;

  /* Enclose argument with parentheses if:
     - current opcode is not uminus, or
     - current opcode is uminus and argument has lowest precedence
  */
  if (op_code != oUminus
      || (op_code == oUminus
          && arg->precedence(output_type, temporary_terms) < precedence(output_type, temporary_terms)))
    {
      output << LEFT_PAR(output_type);
      if (op_code == oSign && (output_type == oCDynamicModel || output_type == oCStaticModel))
        output << "1.0,";
      close_parenthesis = true;
    }

  // Write argument
  arg->writeOutput(output, output_type, temporary_terms, tef_terms);

  if (close_parenthesis)
    output << RIGHT_PAR(output_type);

  // Close parenthesis for uminus
  if (op_code == oUminus)
    output << RIGHT_PAR(output_type);
}

void
UnaryOpNode::writeExternalFunctionOutput(ostream &output, ExprNodeOutputType output_type,
                                         const temporary_terms_t &temporary_terms,
                                         deriv_node_temp_terms_t &tef_terms) const
{
  arg->writeExternalFunctionOutput(output, output_type, temporary_terms, tef_terms);
}

void
UnaryOpNode::writeJsonExternalFunctionOutput(vector<string> &efout,
                                             const temporary_terms_t &temporary_terms,
                                             deriv_node_temp_terms_t &tef_terms,
                                             const bool isdynamic) const
{
  arg->writeJsonExternalFunctionOutput(efout, temporary_terms, tef_terms, isdynamic);
}

void
UnaryOpNode::compileExternalFunctionOutput(ostream &CompileCode, unsigned int &instruction_number,
                                           bool lhs_rhs, const temporary_terms_t &temporary_terms,
                                           const map_idx_t &map_idx, bool dynamic, bool steady_dynamic,
                                           deriv_node_temp_terms_t &tef_terms) const
{
  arg->compileExternalFunctionOutput(CompileCode, instruction_number, lhs_rhs, temporary_terms, map_idx,
                                     dynamic, steady_dynamic, tef_terms);
}

double
UnaryOpNode::eval_opcode(UnaryOpcode op_code, double v) throw (EvalException, EvalExternalFunctionException)
{
  switch (op_code)
    {
    case oUminus:
      return (-v);
    case oExp:
      return (exp(v));
    case oLog:
      return (log(v));
    case oLog10:
      return (log10(v));
    case oCos:
      return (cos(v));
    case oSin:
      return (sin(v));
    case oTan:
      return (tan(v));
    case oAcos:
      return (acos(v));
    case oAsin:
      return (asin(v));
    case oAtan:
      return (atan(v));
    case oCosh:
      return (cosh(v));
    case oSinh:
      return (sinh(v));
    case oTanh:
      return (tanh(v));
    case oAcosh:
      return (acosh(v));
    case oAsinh:
      return (asinh(v));
    case oAtanh:
      return (atanh(v));
    case oSqrt:
      return (sqrt(v));
    case oAbs:
      return (abs(v));
    case oSign:
      return (v > 0) ? 1 : ((v < 0) ? -1 : 0);
    case oSteadyState:
      return (v);
    case oSteadyStateParamDeriv:
    case oSteadyStateParam2ndDeriv:
    case oExpectation:
    case oErf:
      return (erf(v));
    case oDiff:
      cerr << "UnaryOpNode::eval_opcode: not implemented on oDiff" << endl;
      exit(EXIT_FAILURE);
    case oAdl:
      cerr << "UnaryOpNode::eval_opcode: not implemented on oAdl" << endl;
      exit(EXIT_FAILURE);
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

double
UnaryOpNode::eval(const eval_context_t &eval_context) const throw (EvalException, EvalExternalFunctionException)
{
  double v = arg->eval(eval_context);

  return eval_opcode(op_code, v);
}

void
UnaryOpNode::compile(ostream &CompileCode, unsigned int &instruction_number,
                     bool lhs_rhs, const temporary_terms_t &temporary_terms,
                     const map_idx_t &map_idx, bool dynamic, bool steady_dynamic,
                     deriv_node_temp_terms_t &tef_terms) const
{
  temporary_terms_t::const_iterator it = temporary_terms.find(const_cast<UnaryOpNode *>(this));
  if (it != temporary_terms.end())
    {
      if (dynamic)
        {
          map_idx_t::const_iterator ii = map_idx.find(idx);
          FLDT_ fldt(ii->second);
          fldt.write(CompileCode, instruction_number);
        }
      else
        {
          map_idx_t::const_iterator ii = map_idx.find(idx);
          FLDST_ fldst(ii->second);
          fldst.write(CompileCode, instruction_number);
        }
      return;
    }
  if (op_code == oSteadyState)
    arg->compile(CompileCode, instruction_number, lhs_rhs, temporary_terms, map_idx, dynamic, true, tef_terms);
  else
    {
      arg->compile(CompileCode, instruction_number, lhs_rhs, temporary_terms, map_idx, dynamic, steady_dynamic, tef_terms);
      FUNARY_ funary(op_code);
      funary.write(CompileCode, instruction_number);
    }
}

void
UnaryOpNode::collectVARLHSVariable(set<expr_t> &result) const
{
  if (op_code == oDiff)
    result.insert(const_cast<UnaryOpNode *>(this));
  else
    arg->collectVARLHSVariable(result);
}

void
UnaryOpNode::collectDynamicVariables(SymbolType type_arg, set<pair<int, int> > &result) const
{
  arg->collectDynamicVariables(type_arg, result);
}

pair<int, expr_t>
UnaryOpNode::normalizeEquation(int var_endo, vector<pair<int, pair<expr_t, expr_t> > > &List_of_Op_RHS) const
{
  pair<bool, expr_t > res = arg->normalizeEquation(var_endo, List_of_Op_RHS);
  int is_endogenous_present = res.first;
  expr_t New_expr_t = res.second;

  if (is_endogenous_present == 2) /* The equation could not be normalized and the process is given-up*/
    return (make_pair(2, (expr_t) NULL));
  else if (is_endogenous_present) /* The argument of the function contains the current values of
                                     the endogenous variable associated to the equation.
                                     In order to normalized, we have to apply the invert function to the RHS.*/
    {
      switch (op_code)
        {
        case oUminus:
          List_of_Op_RHS.push_back(make_pair(oUminus, make_pair((expr_t) NULL, (expr_t) NULL)));
          return (make_pair(1, (expr_t) NULL));
        case oExp:
          List_of_Op_RHS.push_back(make_pair(oLog, make_pair((expr_t) NULL, (expr_t) NULL)));
          return (make_pair(1, (expr_t) NULL));
        case oLog:
          List_of_Op_RHS.push_back(make_pair(oExp, make_pair((expr_t) NULL, (expr_t) NULL)));
          return (make_pair(1, (expr_t) NULL));
        case oLog10:
          List_of_Op_RHS.push_back(make_pair(oPower, make_pair((expr_t) NULL, datatree.AddNonNegativeConstant("10"))));
          return (make_pair(1, (expr_t) NULL));
        case oCos:
          List_of_Op_RHS.push_back(make_pair(oAcos, make_pair((expr_t) NULL, (expr_t) NULL)));
          return (make_pair(1, (expr_t) NULL));
        case oSin:
          List_of_Op_RHS.push_back(make_pair(oAsin, make_pair((expr_t) NULL, (expr_t) NULL)));
          return (make_pair(1, (expr_t) NULL));
        case oTan:
          List_of_Op_RHS.push_back(make_pair(oAtan, make_pair((expr_t) NULL, (expr_t) NULL)));
          return (make_pair(1, (expr_t) NULL));
        case oAcos:
          List_of_Op_RHS.push_back(make_pair(oCos, make_pair((expr_t) NULL, (expr_t) NULL)));
          return (make_pair(1, (expr_t) NULL));
        case oAsin:
          List_of_Op_RHS.push_back(make_pair(oSin, make_pair((expr_t) NULL, (expr_t) NULL)));
          return (make_pair(1, (expr_t) NULL));
        case oAtan:
          List_of_Op_RHS.push_back(make_pair(oTan, make_pair((expr_t) NULL, (expr_t) NULL)));
          return (make_pair(1, (expr_t) NULL));
        case oCosh:
          List_of_Op_RHS.push_back(make_pair(oAcosh, make_pair((expr_t) NULL, (expr_t) NULL)));
          return (make_pair(1, (expr_t) NULL));
        case oSinh:
          List_of_Op_RHS.push_back(make_pair(oAsinh, make_pair((expr_t) NULL, (expr_t) NULL)));
          return (make_pair(1, (expr_t) NULL));
        case oTanh:
          List_of_Op_RHS.push_back(make_pair(oAtanh, make_pair((expr_t) NULL, (expr_t) NULL)));
          return (make_pair(1, (expr_t) NULL));
        case oAcosh:
          List_of_Op_RHS.push_back(make_pair(oCosh, make_pair((expr_t) NULL, (expr_t) NULL)));
          return (make_pair(1, (expr_t) NULL));
        case oAsinh:
          List_of_Op_RHS.push_back(make_pair(oSinh, make_pair((expr_t) NULL, (expr_t) NULL)));
          return (make_pair(1, (expr_t) NULL));
        case oAtanh:
          List_of_Op_RHS.push_back(make_pair(oTanh, make_pair((expr_t) NULL, (expr_t) NULL)));
          return (make_pair(1, (expr_t) NULL));
        case oSqrt:
          List_of_Op_RHS.push_back(make_pair(oPower, make_pair((expr_t) NULL, datatree.Two)));
          return (make_pair(1, (expr_t) NULL));
        case oAbs:
          return (make_pair(2, (expr_t) NULL));
        case oSign:
          return (make_pair(2, (expr_t) NULL));
        case oSteadyState:
          return (make_pair(2, (expr_t) NULL));
        case oErf:
          return (make_pair(2, (expr_t) NULL));
        default:
          cerr << "Unary operator not handled during the normalization process" << endl;
          return (make_pair(2, (expr_t) NULL)); // Could not be normalized
        }
    }
  else
    { /* If the argument of the function do not contain the current values of the endogenous variable
         related to the equation, the function with its argument is stored in the RHS*/
      switch (op_code)
        {
        case oUminus:
          return (make_pair(0, datatree.AddUMinus(New_expr_t)));
        case oExp:
          return (make_pair(0, datatree.AddExp(New_expr_t)));
        case oLog:
          return (make_pair(0, datatree.AddLog(New_expr_t)));
        case oLog10:
          return (make_pair(0, datatree.AddLog10(New_expr_t)));
        case oCos:
          return (make_pair(0, datatree.AddCos(New_expr_t)));
        case oSin:
          return (make_pair(0, datatree.AddSin(New_expr_t)));
        case oTan:
          return (make_pair(0, datatree.AddTan(New_expr_t)));
        case oAcos:
          return (make_pair(0, datatree.AddAcos(New_expr_t)));
        case oAsin:
          return (make_pair(0, datatree.AddAsin(New_expr_t)));
        case oAtan:
          return (make_pair(0, datatree.AddAtan(New_expr_t)));
        case oCosh:
          return (make_pair(0, datatree.AddCosh(New_expr_t)));
        case oSinh:
          return (make_pair(0, datatree.AddSinh(New_expr_t)));
        case oTanh:
          return (make_pair(0, datatree.AddTanh(New_expr_t)));
        case oAcosh:
          return (make_pair(0, datatree.AddAcosh(New_expr_t)));
        case oAsinh:
          return (make_pair(0, datatree.AddAsinh(New_expr_t)));
        case oAtanh:
          return (make_pair(0, datatree.AddAtanh(New_expr_t)));
        case oSqrt:
          return (make_pair(0, datatree.AddSqrt(New_expr_t)));
        case oAbs:
          return (make_pair(0, datatree.AddAbs(New_expr_t)));
        case oSign:
          return (make_pair(0, datatree.AddSign(New_expr_t)));
        case oSteadyState:
          return (make_pair(0, datatree.AddSteadyState(New_expr_t)));
        case oErf:
          return (make_pair(0, datatree.AddErf(New_expr_t)));
        default:
          cerr << "Unary operator not handled during the normalization process" << endl;
          return (make_pair(2, (expr_t) NULL)); // Could not be normalized
        }
    }
  cerr << "UnaryOpNode::normalizeEquation: impossible case" << endl;
  exit(EXIT_FAILURE);
}

expr_t
UnaryOpNode::getChainRuleDerivative(int deriv_id, const map<int, expr_t> &recursive_variables)
{
  expr_t darg = arg->getChainRuleDerivative(deriv_id, recursive_variables);
  return composeDerivatives(darg, deriv_id);
}

expr_t
UnaryOpNode::buildSimilarUnaryOpNode(expr_t alt_arg, DataTree &alt_datatree) const
{
  switch (op_code)
    {
    case oUminus:
      return alt_datatree.AddUMinus(alt_arg);
    case oExp:
      return alt_datatree.AddExp(alt_arg);
    case oLog:
      return alt_datatree.AddLog(alt_arg);
    case oLog10:
      return alt_datatree.AddLog10(alt_arg);
    case oCos:
      return alt_datatree.AddCos(alt_arg);
    case oSin:
      return alt_datatree.AddSin(alt_arg);
    case oTan:
      return alt_datatree.AddTan(alt_arg);
    case oAcos:
      return alt_datatree.AddAcos(alt_arg);
    case oAsin:
      return alt_datatree.AddAsin(alt_arg);
    case oAtan:
      return alt_datatree.AddAtan(alt_arg);
    case oCosh:
      return alt_datatree.AddCosh(alt_arg);
    case oSinh:
      return alt_datatree.AddSinh(alt_arg);
    case oTanh:
      return alt_datatree.AddTanh(alt_arg);
    case oAcosh:
      return alt_datatree.AddAcosh(alt_arg);
    case oAsinh:
      return alt_datatree.AddAsinh(alt_arg);
    case oAtanh:
      return alt_datatree.AddAtanh(alt_arg);
    case oSqrt:
      return alt_datatree.AddSqrt(alt_arg);
    case oAbs:
      return alt_datatree.AddAbs(alt_arg);
    case oSign:
      return alt_datatree.AddSign(alt_arg);
    case oSteadyState:
      return alt_datatree.AddSteadyState(alt_arg);
    case oSteadyStateParamDeriv:
      cerr << "UnaryOpNode::buildSimilarUnaryOpNode: oSteadyStateParamDeriv can't be translated" << endl;
      exit(EXIT_FAILURE);
    case oSteadyStateParam2ndDeriv:
      cerr << "UnaryOpNode::buildSimilarUnaryOpNode: oSteadyStateParam2ndDeriv can't be translated" << endl;
      exit(EXIT_FAILURE);
    case oExpectation:
      return alt_datatree.AddExpectation(expectation_information_set, alt_arg);
    case oErf:
      return alt_datatree.AddErf(alt_arg);
    case oDiff:
      return alt_datatree.AddDiff(alt_arg);
    case oAdl:
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
UnaryOpNode::cloneDynamic(DataTree &dynamic_datatree) const
{
  expr_t substarg = arg->cloneDynamic(dynamic_datatree);
  return buildSimilarUnaryOpNode(substarg, dynamic_datatree);
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
  if (op_code == oDiff)
    return arg->maxLag() + 1;
  return arg->maxLag();
}

expr_t
UnaryOpNode::undiff() const
{
  if (op_code == oDiff)
    return arg;
  return arg->undiff();
}

void
UnaryOpNode::VarMaxLag(DataTree &static_datatree, set<expr_t> &static_lhs, int &max_lag) const
{
  if (op_code != oDiff)
    arg->VarMaxLag(static_datatree, static_lhs, max_lag);
  else
    {
      for (set<expr_t>::const_iterator it = static_lhs.begin();
           it != static_lhs.end(); it++)
        if (*it == this->toStatic(static_datatree))
          {
            int max_lag_tmp = arg->maxLag();
            if (max_lag_tmp > max_lag)
              max_lag = max_lag_tmp;
            return;
          }
      int max_lag_tmp = 0;
      arg->VarMaxLag(static_datatree, static_lhs, max_lag_tmp);
      if (max_lag_tmp + 1 > max_lag)
        max_lag = max_lag_tmp + 1;
    }
}

int
UnaryOpNode::VarMinLag() const
{
  return arg->VarMinLag();
}

int
UnaryOpNode::PacMaxLag(vector<int> &lhs) const
{
  //This will never be an oDiff node
  return arg->PacMaxLag(lhs);
}

expr_t
UnaryOpNode::substituteAdl() const
{
  if (op_code != oAdl)
    {
      expr_t argsubst = arg->substituteAdl();
      return buildSimilarUnaryOpNode(argsubst, datatree);
    }

  expr_t arg1subst = arg->substituteAdl();
  expr_t retval = NULL;
  ostringstream inttostr;

  for (vector<int>::const_iterator it = adl_lags.begin(); it != adl_lags.end(); it++)
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

bool
UnaryOpNode::isDiffPresent() const
{
  if (op_code == oDiff)
    return true;
  return arg->isDiffPresent();
}

void
UnaryOpNode::findDiffUnaryOpNodes(DataTree &static_datatree, set<UnaryOpNode *> &unary_op_nodes_in_diff) const
{
  if (op_code != oDiff)
    {
      arg->findDiffUnaryOpNodes(static_datatree, unary_op_nodes_in_diff);
      return;
    }

  VariableNode *vn = dynamic_cast<VariableNode *>(arg);
  if (vn != NULL)
    return;

  UnaryOpNode *sarg = dynamic_cast<UnaryOpNode *>(arg->toStatic(static_datatree));
  if (sarg == NULL || (!sarg->createAuxVarForUnaryOpNodeInDiffOp() && sarg->get_op_code() != oDiff))
    {
      cerr << "Error: diff operator can only contain certain unary operations" << endl;
      exit(EXIT_FAILURE);
    }

  if (unary_op_nodes_in_diff.find(sarg) != unary_op_nodes_in_diff.end())
    return;

  if (sarg->get_op_code() == oDiff)
    {
      arg->findDiffUnaryOpNodes(static_datatree, unary_op_nodes_in_diff);
      return;
    }

  vn = dynamic_cast<VariableNode *>(sarg->get_arg());
  if (vn == NULL)
    {
      cerr << "Error: A unary op inside of a diff can only take a variable as an argument" << endl;
      exit(EXIT_FAILURE);
    }

  unary_op_nodes_in_diff.insert(sarg);
}

bool
UnaryOpNode::createAuxVarForUnaryOpNodeInDiffOp() const
{
  switch (op_code)
    {
    case oExp:
    case oLog:
    case oLog10:
    case oCos:
    case oSin:
    case oTan:
    case oAcos:
    case oAsin:
    case oAtan:
    case oCosh:
    case oSinh:
    case oTanh:
    case oAcosh:
    case oAsinh:
    case oAtanh:
    case oSqrt:
    case oAbs:
    case oSign:
    case oErf:
      return true;
    default:
      return false;
    }
}

void
UnaryOpNode::findDiffNodes(DataTree &static_datatree, diff_table_t &diff_table) const
{
  if (op_code != oDiff)
    return;

  arg->findDiffNodes(static_datatree, diff_table);

  expr_t sthis = this->toStatic(static_datatree);
  int arg_max_lag = -arg->maxLag();
  diff_table_t::iterator it = diff_table.find(sthis);
  if (it != diff_table.end())
    {
      for (map<int, expr_t>::const_iterator it1 = it->second.begin();
           it1 != it->second.end(); it1++)
        if (arg == it1->second)
          return;
      it->second[arg_max_lag] = const_cast<UnaryOpNode *>(this);
    }
  else
    diff_table[sthis][arg_max_lag] = const_cast<UnaryOpNode *>(this);
}

void
UnaryOpNode::getDiffArgUnaryOperatorIfAny(string &op_handle) const
{
  switch (op_code)
    {
    case oExp:
      op_handle = "@exp";
      break;
    case oLog:
      op_handle = "@log";
      break;
    case oLog10:
      op_handle = "@log10";
      break;
    case oCos:
      op_handle = "@cos";
      break;
    case oSin:
      op_handle = "@sin";
      break;
    case oTan:
      op_handle = "@tan";
      break;
    case oAcos:
      op_handle = "@acos";
      break;
    case oAsin:
      op_handle = "@asin";
      break;
    case oAtan:
      op_handle = "@atan";
      break;
    case oCosh:
      op_handle = "@cosh";
      break;
    case oSinh:
      op_handle = "@sinh";
      break;
    case oTanh:
      op_handle = "@tanh";
      break;
    case oAcosh:
      op_handle = "@acosh";
      break;
    case oAsinh:
      op_handle = "@asinh";
      break;
    case oAtanh:
      op_handle = "@atanh";
      break;
    case oSqrt:
      op_handle = "@sqrt";
      break;
    case oAbs:
      op_handle = "@abs";
      break;
    case oSign:
      op_handle = "@sign";
      break;
    case oErf:
      op_handle = "@erf";
      break;
    default:
      op_handle = "";
      break;
    }
}

expr_t
UnaryOpNode::substituteDiff(DataTree &static_datatree, diff_table_t &diff_table, subst_table_t &subst_table,
                            vector<BinaryOpNode *> &neweqs) const
{
  if (op_code != oDiff)
    {
      expr_t argsubst = arg->substituteDiff(static_datatree, diff_table, subst_table, neweqs);
      return buildSimilarUnaryOpNode(argsubst, datatree);
    }

  subst_table_t::const_iterator sit = subst_table.find(this);
  if (sit != subst_table.end())
    return const_cast<VariableNode *>(sit->second);

  expr_t sthis = dynamic_cast<UnaryOpNode *>(this->toStatic(static_datatree));
  diff_table_t::iterator it = diff_table.find(sthis);
  if (it == diff_table.end() || it->second[-arg->maxLag()] != this)
    {
      cerr << "Internal error encountered. Please report" << endl;
      exit(EXIT_FAILURE);
    }

  int last_arg_max_lag = 0;
  VariableNode *last_aux_var = NULL;
  for (map<int, expr_t>::reverse_iterator rit = it->second.rbegin();
       rit != it->second.rend(); rit++)
    {
      expr_t argsubst = dynamic_cast<UnaryOpNode *>(rit->second)->
          get_arg()->substituteDiff(static_datatree, diff_table, subst_table, neweqs);
      int symb_id;
      VariableNode *vn = dynamic_cast<VariableNode *>(argsubst);
      if (rit == it->second.rbegin())
        {
          if (vn != NULL)
            symb_id = datatree.symbol_table.addDiffAuxiliaryVar(argsubst->idx, argsubst, vn->get_symb_id(), vn->get_lag());
          else
            {
              UnaryOpNode *diffarg = dynamic_cast<UnaryOpNode *>(argsubst);
              if (diffarg != NULL)
                {
                  string op;
                  diffarg->getDiffArgUnaryOperatorIfAny(op);
                  VariableNode *vnarg = dynamic_cast<VariableNode *>(diffarg->get_arg());
                  if (vnarg != NULL)
                    symb_id = datatree.symbol_table.addDiffAuxiliaryVar(argsubst->idx, argsubst,
                                                                        vnarg->get_symb_id(), vnarg->get_lag(), op);
                  else
                    {
                      // The case where we have diff(log(exp(x))) for example
                      cerr << "diffs of nested non-diff expressions are not yet supported" << endl;
                      exit(EXIT_FAILURE);
                    }
                }
              else
                {
                  cerr << "diffs of non unary expressions are not yet supported" << endl;
                  exit(EXIT_FAILURE);
                }
            }
          // make originating aux var & equation
          last_arg_max_lag = rit->first;
          last_aux_var = datatree.AddVariable(symb_id, 0);
          //ORIG_AUX_DIFF = argsubst - argsubst(-1)
          neweqs.push_back(dynamic_cast<BinaryOpNode *>(datatree.AddEqual(last_aux_var,
                                                                          datatree.AddMinus(argsubst,
                                                                                            argsubst->decreaseLeadsLags(1)))));
          subst_table[rit->second] = dynamic_cast<VariableNode *>(last_aux_var);
        }
      else
        {
          // just add equation of form: AUX_DIFF = LAST_AUX_VAR(-1)
          VariableNode *new_aux_var = NULL;
          for (int i = last_arg_max_lag; i > rit->first; i--)
            {
              if (i == last_arg_max_lag)
                symb_id = datatree.symbol_table.addDiffLagAuxiliaryVar(argsubst->idx, argsubst,
                                                                       last_aux_var->get_symb_id(), last_aux_var->get_lag());
              else
                symb_id = datatree.symbol_table.addDiffLagAuxiliaryVar(new_aux_var->idx, new_aux_var,
                                                                       last_aux_var->get_symb_id(), last_aux_var->get_lag());

              new_aux_var = datatree.AddVariable(symb_id, 0);
              neweqs.push_back(dynamic_cast<BinaryOpNode *>(datatree.AddEqual(new_aux_var,
                                                                              last_aux_var->decreaseLeadsLags(1))));
              last_aux_var = new_aux_var;
            }
          subst_table[rit->second] = dynamic_cast<VariableNode *>(new_aux_var);
          last_arg_max_lag = rit->first;
        }
    }
  return const_cast<VariableNode *>(subst_table[this]);
}

expr_t
UnaryOpNode::substituteDiffUnaryOpNodes(DataTree &static_datatree, set<UnaryOpNode *> &nodes, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  UnaryOpNode *sarg, *argtouse;
  if (op_code == oDiff)
    {
      sarg = dynamic_cast<UnaryOpNode *>(arg->toStatic(static_datatree));
      argtouse = dynamic_cast<UnaryOpNode *>(arg);
    }
  else
    {
      sarg = dynamic_cast<UnaryOpNode *>(this->toStatic(static_datatree));
      argtouse = const_cast<UnaryOpNode *>(this);
    }

  if (sarg == NULL || nodes.find(sarg) == nodes.end())
    {
      expr_t argsubst = arg->substituteDiffUnaryOpNodes(static_datatree, nodes, subst_table, neweqs);
      return buildSimilarUnaryOpNode(argsubst, datatree);
    }

  subst_table_t::const_iterator sit = subst_table.find(argtouse);
  if (sit != subst_table.end())
    return const_cast<VariableNode *>(sit->second);

  VariableNode *vn = dynamic_cast<VariableNode *>(dynamic_cast<UnaryOpNode *>(argtouse)->get_arg());
  int lag = vn->get_lag();
  int symb_id = datatree.symbol_table.addUnaryOpInsideDiffAuxiliaryVar(argtouse->idx, argtouse,
                                                                       vn->get_symb_id(), lag);
  VariableNode *aux_var = datatree.AddVariable(symb_id, 0);
  neweqs.push_back(dynamic_cast<BinaryOpNode *>(datatree.AddEqual(aux_var, argtouse)));
  subst_table[argtouse] = dynamic_cast<VariableNode *>(aux_var);
  if (op_code == oDiff)
    return buildSimilarUnaryOpNode(aux_var, datatree);
  return const_cast<VariableNode *>(aux_var);
}

expr_t
UnaryOpNode::substitutePacExpectation(map<const PacExpectationNode *, const BinaryOpNode *> &subst_table)
{
  expr_t argsubst = arg->substitutePacExpectation(subst_table);
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
  if (op_code == oUminus || deterministic_model)
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
  if (op_code == oUminus || deterministic_model)
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
  if (op_code == oExpectation)
    {
      subst_table_t::iterator it = subst_table.find(const_cast<UnaryOpNode *>(this));
      if (it != subst_table.end())
        return const_cast<VariableNode *>(it->second);

      //Arriving here, we need to create an auxiliary variable for this Expectation Operator:
      //AUX_EXPECT_(LEAD/LAG)_(period)_(arg.idx) OR
      //AUX_EXPECT_(info_set_name)_(arg.idx)
      int symb_id = datatree.symbol_table.addExpectationAuxiliaryVar(expectation_information_set, arg->idx, arg);
      expr_t newAuxE = datatree.AddVariable(symb_id, 0);

      if (partial_information_model && expectation_information_set == 0)
        if (dynamic_cast<VariableNode *>(arg) == NULL)
          {
            cerr << "ERROR: In Partial Information models, EXPECTATION(0)(X) "
                 << "can only be used when X is a single variable." << endl;
            exit(EXIT_FAILURE);
          }

      //take care of any nested expectation operators by calling arg->substituteExpectation(.), then decreaseLeadsLags for this oExpectation operator
      //arg(lag-period) (holds entire subtree of arg(lag-period)
      expr_t substexpr = (arg->substituteExpectation(subst_table, neweqs, partial_information_model))->decreaseLeadsLags(expectation_information_set);
      assert(substexpr != NULL);
      neweqs.push_back(dynamic_cast<BinaryOpNode *>(datatree.AddEqual(newAuxE, substexpr))); //AUXE_period_arg.idx = arg(lag-period)
      newAuxE = datatree.AddVariable(symb_id, expectation_information_set);

      assert(dynamic_cast<VariableNode *>(newAuxE) != NULL);
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
UnaryOpNode::containsEndogenous(void) const
{
  return arg->containsEndogenous();
}

bool
UnaryOpNode::containsExogenous() const
{
  return arg->containsExogenous();
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
UnaryOpNode::removeTrendLeadLag(map<int, expr_t> trend_symbols_map) const
{
  expr_t argsubst = arg->removeTrendLeadLag(trend_symbols_map);
  return buildSimilarUnaryOpNode(argsubst, datatree);
}

bool
UnaryOpNode::isInStaticForm() const
{
  if (op_code == oSteadyState || op_code == oSteadyStateParamDeriv
      || op_code == oSteadyStateParam2ndDeriv
      || op_code == oExpectation)
    return false;
  else
    return arg->isInStaticForm();
}

void
UnaryOpNode::setVarExpectationIndex(map<string, pair<SymbolList, int> > &var_model_info)
{
  arg->setVarExpectationIndex(var_model_info);
}

void
UnaryOpNode::walkPacParameters(bool &pac_encountered, pair<int, int> &lhs, set<pair<int, pair<int, int> > > &ec_params_and_vars, set<pair<int, pair<int, int> > > &ar_params_and_vars) const
{
  arg->walkPacParameters(pac_encountered, lhs, ec_params_and_vars, ar_params_and_vars);
}

void
UnaryOpNode::addParamInfoToPac(pair<int, int> &lhs_arg, set<pair<int, pair<int, int> > > &ec_params_and_vars_arg, set<pair<int, pair<int, int> > > &ar_params_and_vars_arg)
{
  arg->addParamInfoToPac(lhs_arg, ec_params_and_vars_arg, ar_params_and_vars_arg);
}

void
UnaryOpNode::fillPacExpectationVarInfo(string &model_name_arg, vector<int> &lhs_arg, int max_lag_arg, vector<bool> &nonstationary_arg, int growth_symb_id_arg, int equation_number_arg)
{
  arg->fillPacExpectationVarInfo(model_name_arg, lhs_arg, max_lag_arg, nonstationary_arg, growth_symb_id_arg, equation_number_arg);
}

bool
UnaryOpNode::isVarModelReferenced(const string &model_info_name) const
{
  return arg->isVarModelReferenced(model_info_name);
}

void
UnaryOpNode::getEndosAndMaxLags(map<string, int> &model_endos_and_lags) const
{
  arg->getEndosAndMaxLags(model_endos_and_lags);
}

expr_t
UnaryOpNode::substituteStaticAuxiliaryVariable() const
{
  expr_t argsubst = arg->substituteStaticAuxiliaryVariable();
  if (op_code == oExpectation)
    return argsubst;
  else
    return buildSimilarUnaryOpNode(argsubst, datatree);
}

BinaryOpNode::BinaryOpNode(DataTree &datatree_arg, const expr_t arg1_arg,
                           BinaryOpcode op_code_arg, const expr_t arg2_arg) :
  ExprNode(datatree_arg),
  arg1(arg1_arg),
  arg2(arg2_arg),
  op_code(op_code_arg),
  powerDerivOrder(0)
{
  datatree.binary_op_node_map[make_pair(make_pair(make_pair(arg1, arg2), powerDerivOrder), op_code)] = this;
}

BinaryOpNode::BinaryOpNode(DataTree &datatree_arg, const expr_t arg1_arg,
                           BinaryOpcode op_code_arg, const expr_t arg2_arg, int powerDerivOrder_arg) :
  ExprNode(datatree_arg),
  arg1(arg1_arg),
  arg2(arg2_arg),
  op_code(op_code_arg),
  powerDerivOrder(powerDerivOrder_arg)
{
  assert(powerDerivOrder >= 0);
  datatree.binary_op_node_map[make_pair(make_pair(make_pair(arg1, arg2), powerDerivOrder), op_code)] = this;
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
    case oPlus:
      return datatree.AddPlus(darg1, darg2);
    case oMinus:
      return datatree.AddMinus(darg1, darg2);
    case oTimes:
      t11 = datatree.AddTimes(darg1, arg2);
      t12 = datatree.AddTimes(darg2, arg1);
      return datatree.AddPlus(t11, t12);
    case oDivide:
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
    case oLess:
    case oGreater:
    case oLessEqual:
    case oGreaterEqual:
    case oEqualEqual:
    case oDifferent:
      return datatree.Zero;
    case oPower:
      if (darg2 == datatree.Zero)
        if (darg1 == datatree.Zero)
          return datatree.Zero;
        else
          if (dynamic_cast<NumConstNode *>(arg2) != NULL)
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
    case oPowerDeriv:
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
          expr_t first_part  = datatree.AddTimes(f, t15);

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
    case oMax:
      t11 = datatree.AddGreater(arg1, arg2);
      t12 = datatree.AddTimes(t11, darg1);
      t13 = datatree.AddMinus(datatree.One, t11);
      t14 = datatree.AddTimes(t13, darg2);
      return datatree.AddPlus(t14, t12);
    case oMin:
      t11 = datatree.AddGreater(arg2, arg1);
      t12 = datatree.AddTimes(t11, darg1);
      t13 = datatree.AddMinus(datatree.One, t11);
      t14 = datatree.AddTimes(t13, darg2);
      return datatree.AddPlus(t14, t12);
    case oEqual:
      return datatree.AddMinus(darg1, darg2);
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

expr_t
BinaryOpNode::unpackPowerDeriv() const
{
  if (op_code != oPowerDeriv)
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
  temporary_terms_t::const_iterator it = temporary_terms.find(const_cast<BinaryOpNode *>(this));
  // A temporary term behaves as a variable
  if (it != temporary_terms.end())
    return 100;

  switch (op_code)
    {
    case oEqual:
      return 0;
    case oEqualEqual:
    case oDifferent:
      return 1;
    case oLessEqual:
    case oGreaterEqual:
    case oLess:
    case oGreater:
      return 2;
    case oPlus:
    case oMinus:
      return 3;
    case oTimes:
    case oDivide:
      return 4;
    case oPower:
    case oPowerDeriv:
      if (IS_C(output_type))
        // In C, power operator is of the form pow(a, b)
        return 100;
      else
        return 5;
    case oMin:
    case oMax:
      return 100;
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

int
BinaryOpNode::precedenceJson(const temporary_terms_t &temporary_terms) const
{
  temporary_terms_t::const_iterator it = temporary_terms.find(const_cast<BinaryOpNode *>(this));
  // A temporary term behaves as a variable
  if (it != temporary_terms.end())
    return 100;

  switch (op_code)
    {
    case oEqual:
      return 0;
    case oEqualEqual:
    case oDifferent:
      return 1;
    case oLessEqual:
    case oGreaterEqual:
    case oLess:
    case oGreater:
      return 2;
    case oPlus:
    case oMinus:
      return 3;
    case oTimes:
    case oDivide:
      return 4;
    case oPower:
    case oPowerDeriv:
      return 5;
    case oMin:
    case oMax:
      return 100;
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

int
BinaryOpNode::cost(const map<NodeTreeReference, temporary_terms_t> &temp_terms_map, bool is_matlab) const
{
  // For a temporary term, the cost is null
  for (map<NodeTreeReference, temporary_terms_t>::const_iterator it = temp_terms_map.begin();
       it != temp_terms_map.end(); it++)
    if (it->second.find(const_cast<BinaryOpNode *>(this)) != it->second.end())
      return 0;

  int arg_cost = arg1->cost(temp_terms_map, is_matlab) + arg2->cost(temp_terms_map, is_matlab);

  return cost(arg_cost, is_matlab);
}

int
BinaryOpNode::cost(const temporary_terms_t &temporary_terms, bool is_matlab) const
{
  // For a temporary term, the cost is null
  if (temporary_terms.find(const_cast<BinaryOpNode *>(this)) != temporary_terms.end())
    return 0;

  int arg_cost = arg1->cost(temporary_terms, is_matlab) + arg2->cost(temporary_terms, is_matlab);

  return cost(arg_cost, is_matlab);
}

int
BinaryOpNode::cost(int cost, bool is_matlab) const
{
  if (is_matlab)
    // Cost for Matlab files
    switch (op_code)
      {
      case oLess:
      case oGreater:
      case oLessEqual:
      case oGreaterEqual:
      case oEqualEqual:
      case oDifferent:
        return cost + 60;
      case oPlus:
      case oMinus:
      case oTimes:
        return cost + 90;
      case oMax:
      case oMin:
        return cost + 110;
      case oDivide:
        return cost + 990;
      case oPower:
      case oPowerDeriv:
        return cost + (MIN_COST_MATLAB/2+1);
      case oEqual:
        return cost;
      }
  else
    // Cost for C files
    switch (op_code)
      {
      case oLess:
      case oGreater:
      case oLessEqual:
      case oGreaterEqual:
      case oEqualEqual:
      case oDifferent:
        return cost + 2;
      case oPlus:
      case oMinus:
      case oTimes:
        return cost + 4;
      case oMax:
      case oMin:
        return cost + 5;
      case oDivide:
        return cost + 15;
      case oPower:
        return cost + 520;
      case oPowerDeriv:
        return cost + (MIN_COST_C/2+1);;
      case oEqual:
        return cost;
      }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

void
BinaryOpNode::computeTemporaryTerms(map<expr_t, pair<int, NodeTreeReference> > &reference_count,
                                    map<NodeTreeReference, temporary_terms_t> &temp_terms_map,
                                    bool is_matlab, NodeTreeReference tr) const
{
  expr_t this2 = const_cast<BinaryOpNode *>(this);
  map<expr_t, pair<int, NodeTreeReference > >::iterator it = reference_count.find(this2);
  if (it == reference_count.end())
    {
      // If this node has never been encountered, set its ref count to one,
      //  and travel through its children
      reference_count[this2] = make_pair(1, tr);
      arg1->computeTemporaryTerms(reference_count, temp_terms_map, is_matlab, tr);
      arg2->computeTemporaryTerms(reference_count, temp_terms_map, is_matlab, tr);
    }
  else
    {
      /* If the node has already been encountered, increment its ref count
         and declare it as a temporary term if it is too costly (except if it is
         an equal node: we don't want them as temporary terms) */
      reference_count[this2] = make_pair(it->second.first + 1, it->second.second);;
      if (reference_count[this2].first * cost(temp_terms_map, is_matlab) > MIN_COST(is_matlab)
          && op_code != oEqual)
        temp_terms_map[reference_count[this2].second].insert(this2);
    }
}

void
BinaryOpNode::computeTemporaryTerms(map<expr_t, int> &reference_count,
                                    temporary_terms_t &temporary_terms,
                                    map<expr_t, pair<int, int> > &first_occurence,
                                    int Curr_block,
                                    vector<vector<temporary_terms_t> > &v_temporary_terms,
                                    int equation) const
{
  expr_t this2 = const_cast<BinaryOpNode *>(this);
  map<expr_t, int>::iterator it = reference_count.find(this2);
  if (it == reference_count.end())
    {
      reference_count[this2] = 1;
      first_occurence[this2] = make_pair(Curr_block, equation);
      arg1->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, Curr_block, v_temporary_terms, equation);
      arg2->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, Curr_block, v_temporary_terms, equation);
    }
  else
    {
      reference_count[this2]++;
      if (reference_count[this2] * cost(temporary_terms, false) > MIN_COST_C
          && op_code != oEqual)
        {
          temporary_terms.insert(this2);
          v_temporary_terms[first_occurence[this2].first][first_occurence[this2].second].insert(this2);
        }
    }
}

double
BinaryOpNode::eval_opcode(double v1, BinaryOpcode op_code, double v2, int derivOrder) throw (EvalException, EvalExternalFunctionException)
{
  switch (op_code)
    {
    case oPlus:
      return (v1 + v2);
    case oMinus:
      return (v1 - v2);
    case oTimes:
      return (v1 * v2);
    case oDivide:
      return (v1 / v2);
    case oPower:
      return (pow(v1, v2));
    case oPowerDeriv:
      if (fabs(v1) < NEAR_ZERO && v2 > 0
          && derivOrder > v2
          && fabs(v2-nearbyint(v2)) < NEAR_ZERO)
        return 0.0;
      else
        {
          double dxp = pow(v1, v2-derivOrder);
          for (int i = 0; i < derivOrder; i++)
            dxp *= v2--;
          return dxp;
        }
    case oMax:
      if (v1 < v2)
        return v2;
      else
        return v1;
    case oMin:
      if (v1 > v2)
        return v2;
      else
        return v1;
    case oLess:
      return (v1 < v2);
    case oGreater:
      return (v1 > v2);
    case oLessEqual:
      return (v1 <= v2);
    case oGreaterEqual:
      return (v1 >= v2);
    case oEqualEqual:
      return (v1 == v2);
    case oDifferent:
      return (v1 != v2);
    case oEqual:
      throw EvalException();
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

double
BinaryOpNode::eval(const eval_context_t &eval_context) const throw (EvalException, EvalExternalFunctionException)
{
  double v1 = arg1->eval(eval_context);
  double v2 = arg2->eval(eval_context);

  return eval_opcode(v1, op_code, v2, powerDerivOrder);
}

void
BinaryOpNode::compile(ostream &CompileCode, unsigned int &instruction_number,
                      bool lhs_rhs, const temporary_terms_t &temporary_terms,
                      const map_idx_t &map_idx, bool dynamic, bool steady_dynamic,
                      deriv_node_temp_terms_t &tef_terms) const
{
  // If current node is a temporary term
  temporary_terms_t::const_iterator it = temporary_terms.find(const_cast<BinaryOpNode *>(this));
  if (it != temporary_terms.end())
    {
      if (dynamic)
        {
          map_idx_t::const_iterator ii = map_idx.find(idx);
          FLDT_ fldt(ii->second);
          fldt.write(CompileCode, instruction_number);
        }
      else
        {
          map_idx_t::const_iterator ii = map_idx.find(idx);
          FLDST_ fldst(ii->second);
          fldst.write(CompileCode, instruction_number);
        }
      return;
    }
  if (op_code == oPowerDeriv)
    {
      FLDC_ fldc(powerDerivOrder);
      fldc.write(CompileCode, instruction_number);
    }
  arg1->compile(CompileCode, instruction_number, lhs_rhs, temporary_terms, map_idx, dynamic, steady_dynamic, tef_terms);
  arg2->compile(CompileCode, instruction_number, lhs_rhs, temporary_terms, map_idx, dynamic, steady_dynamic, tef_terms);
  FBINARY_ fbinary(op_code);
  fbinary.write(CompileCode, instruction_number);
}

void
BinaryOpNode::collectTemporary_terms(const temporary_terms_t &temporary_terms, temporary_terms_inuse_t &temporary_terms_inuse, int Curr_Block) const
{
  temporary_terms_t::const_iterator it = temporary_terms.find(const_cast<BinaryOpNode *>(this));
  if (it != temporary_terms.end())
    temporary_terms_inuse.insert(idx);
  else
    {
      arg1->collectTemporary_terms(temporary_terms, temporary_terms_inuse, Curr_Block);
      arg2->collectTemporary_terms(temporary_terms, temporary_terms_inuse, Curr_Block);
    }
}

bool
BinaryOpNode::containsExternalFunction() const
{
  return arg1->containsExternalFunction()
    || arg2->containsExternalFunction();
}

void
BinaryOpNode::writeJsonOutput(ostream &output,
                              const temporary_terms_t &temporary_terms,
                              deriv_node_temp_terms_t &tef_terms,
                              const bool isdynamic) const
{
  // If current node is a temporary term
  temporary_terms_t::const_iterator it = temporary_terms.find(const_cast<BinaryOpNode *>(this));
  if (it != temporary_terms.end())
    {
      output << "T" << idx;
      return;
    }

  if (op_code == oMax || op_code == oMin)
    {
      switch (op_code)
        {
        case oMax:
          output << "max(";
          break;
        case oMin:
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

  if (op_code == oPowerDeriv)
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
  BinaryOpNode *barg1 = dynamic_cast<BinaryOpNode *>(arg1);
  if (arg1->precedenceJson(temporary_terms) < prec
      || (op_code == oPower && barg1 != NULL && barg1->op_code == oPower))
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
    case oPlus:
      output << "+";
      break;
    case oMinus:
      output << "-";
      break;
    case oTimes:
      output << "*";
      break;
    case oDivide:
      output << "/";
      break;
    case oPower:
      output << "^";
      break;
    case oLess:
      output << "<";
      break;
    case oGreater:
      output << ">";
      break;
    case oLessEqual:
      output << "<=";
      break;
    case oGreaterEqual:
      output << ">=";
      break;
    case oEqualEqual:
      output << "==";
      break;
    case oDifferent:
      output << "!=";
      break;
    case oEqual:
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
  BinaryOpNode *barg2 = dynamic_cast<BinaryOpNode *>(arg2);
  int arg2_prec = arg2->precedenceJson(temporary_terms);
  if (arg2_prec < prec
      || (op_code == oPower && barg2 != NULL && barg2->op_code == oPower)
      || (op_code == oMinus && arg2_prec == prec)
      || (op_code == oDivide && arg2_prec == prec))
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
                          deriv_node_temp_terms_t &tef_terms) const
{
  // If current node is a temporary term
  temporary_terms_t::const_iterator it = temporary_terms.find(const_cast<BinaryOpNode *>(this));
  if (it != temporary_terms.end())
    {
      if (output_type == oMatlabDynamicModelSparse)
        output << "T" << idx << "(it_)";
      else
        output << "T" << idx;
      return;
    }

  // Treat derivative of Power
  if (op_code == oPowerDeriv)
    {
      if (IS_LATEX(output_type))
        unpackPowerDeriv()->writeOutput(output, output_type, temporary_terms, tef_terms);
      else
        {
          if (output_type == oJuliaStaticModel || output_type == oJuliaDynamicModel)
            output << "get_power_deriv(";
          else
            output << "getPowerDeriv(";
          arg1->writeOutput(output, output_type, temporary_terms, tef_terms);
          output << ",";
          arg2->writeOutput(output, output_type, temporary_terms, tef_terms);
          output << "," << powerDerivOrder << ")";
        }
      return;
    }

  // Treat special case of power operator in C, and case of max and min operators
  if ((op_code == oPower && IS_C(output_type)) || op_code == oMax || op_code == oMin)
    {
      switch (op_code)
        {
        case oPower:
          output << "pow(";
          break;
        case oMax:
          output << "max(";
          break;
        case oMin:
          output << "min(";
          break;
        default:
          ;
        }
      arg1->writeOutput(output, output_type, temporary_terms, tef_terms);
      output << ",";
      arg2->writeOutput(output, output_type, temporary_terms, tef_terms);
      output << ")";
      return;
    }

  int prec = precedence(output_type, temporary_terms);

  bool close_parenthesis = false;

  if (IS_LATEX(output_type) && op_code == oDivide)
    output << "\\frac{";
  else
    {
      // If left argument has a lower precedence, or if current and left argument are both power operators, add parenthesis around left argument
      BinaryOpNode *barg1 = dynamic_cast<BinaryOpNode *>(arg1);
      if (arg1->precedence(output_type, temporary_terms) < prec
          || (op_code == oPower && barg1 != NULL && barg1->op_code == oPower))
        {
          output << LEFT_PAR(output_type);
          close_parenthesis = true;
        }
    }

  // Write left argument
  arg1->writeOutput(output, output_type, temporary_terms, tef_terms);

  if (close_parenthesis)
    output << RIGHT_PAR(output_type);

  if (IS_LATEX(output_type) && op_code == oDivide)
    output << "}";

  // Write current operator symbol
  switch (op_code)
    {
    case oPlus:
      output << "+";
      break;
    case oMinus:
      output << "-";
      break;
    case oTimes:
      if (IS_LATEX(output_type))
        output << "\\, ";
      else
        output << "*";
      break;
    case oDivide:
      if (!IS_LATEX(output_type))
        output << "/";
      break;
    case oPower:
      output << "^";
      break;
    case oLess:
      output << "<";
      break;
    case oGreater:
      output << ">";
      break;
    case oLessEqual:
      if (IS_LATEX(output_type))
        output << "\\leq ";
      else
        output << "<=";
      break;
    case oGreaterEqual:
      if (IS_LATEX(output_type))
        output << "\\geq ";
      else
        output << ">=";
      break;
    case oEqualEqual:
      output << "==";
      break;
    case oDifferent:
      if (IS_MATLAB(output_type))
        output << "~=";
      else
        {
          if (IS_C(output_type) || IS_JULIA(output_type))
            output << "!=";
          else
            output << "\\neq ";
        }
      break;
    case oEqual:
      output << "=";
      break;
    default:
      ;
    }

  close_parenthesis = false;

  if (IS_LATEX(output_type) && (op_code == oPower || op_code == oDivide))
    output << "{";
  else
    {
      /* Add parenthesis around right argument if:
         - its precedence is lower than those of the current node
         - it is a power operator and current operator is also a power operator
         - it is a minus operator with same precedence than current operator
         - it is a divide operator with same precedence than current operator */
      BinaryOpNode *barg2 = dynamic_cast<BinaryOpNode *>(arg2);
      int arg2_prec = arg2->precedence(output_type, temporary_terms);
      if (arg2_prec < prec
          || (op_code == oPower && barg2 != NULL && barg2->op_code == oPower && !IS_LATEX(output_type))
          || (op_code == oMinus && arg2_prec == prec)
          || (op_code == oDivide && arg2_prec == prec && !IS_LATEX(output_type)))
        {
          output << LEFT_PAR(output_type);
          close_parenthesis = true;
        }
    }

  // Write right argument
  arg2->writeOutput(output, output_type, temporary_terms, tef_terms);

  if (IS_LATEX(output_type) && (op_code == oPower || op_code == oDivide))
    output << "}";

  if (close_parenthesis)
    output << RIGHT_PAR(output_type);
}

void
BinaryOpNode::writeExternalFunctionOutput(ostream &output, ExprNodeOutputType output_type,
                                          const temporary_terms_t &temporary_terms,
                                          deriv_node_temp_terms_t &tef_terms) const
{
  arg1->writeExternalFunctionOutput(output, output_type, temporary_terms, tef_terms);
  arg2->writeExternalFunctionOutput(output, output_type, temporary_terms, tef_terms);
}

void
BinaryOpNode::writeJsonExternalFunctionOutput(vector<string> &efout,
                                              const temporary_terms_t &temporary_terms,
                                              deriv_node_temp_terms_t &tef_terms,
                                              const bool isdynamic) const
{
  arg1->writeJsonExternalFunctionOutput(efout, temporary_terms, tef_terms, isdynamic);
  arg2->writeJsonExternalFunctionOutput(efout, temporary_terms, tef_terms, isdynamic);
}

void
BinaryOpNode::compileExternalFunctionOutput(ostream &CompileCode, unsigned int &instruction_number,
                                            bool lhs_rhs, const temporary_terms_t &temporary_terms,
                                            const map_idx_t &map_idx, bool dynamic, bool steady_dynamic,
                                            deriv_node_temp_terms_t &tef_terms) const
{
  arg1->compileExternalFunctionOutput(CompileCode, instruction_number, lhs_rhs, temporary_terms, map_idx,
                                      dynamic, steady_dynamic, tef_terms);
  arg2->compileExternalFunctionOutput(CompileCode, instruction_number, lhs_rhs, temporary_terms, map_idx,
                                      dynamic, steady_dynamic, tef_terms);
}

int
BinaryOpNode::VarMinLag() const
{
  return min(arg1->VarMinLag(), arg2->VarMinLag());
}

void
BinaryOpNode::VarMaxLag(DataTree &static_datatree, set<expr_t> &static_lhs, int &max_lag) const
{
  arg1->VarMaxLag(static_datatree, static_lhs, max_lag);
  arg2->VarMaxLag(static_datatree, static_lhs, max_lag);
}

void
BinaryOpNode::collectVARLHSVariable(set<expr_t> &result) const
{
  arg1->collectVARLHSVariable(result);
  arg2->collectVARLHSVariable(result);
}

void
BinaryOpNode::collectDynamicVariables(SymbolType type_arg, set<pair<int, int> > &result) const
{
  arg1->collectDynamicVariables(type_arg, result);
  arg2->collectDynamicVariables(type_arg, result);
}

expr_t
BinaryOpNode::Compute_RHS(expr_t arg1, expr_t arg2, int op, int op_type) const
{
  temporary_terms_t temp;
  switch (op_type)
    {
    case 0: /*Unary Operator*/
      switch (op)
        {
        case oUminus:
          return (datatree.AddUMinus(arg1));
          break;
        case oExp:
          return (datatree.AddExp(arg1));
          break;
        case oLog:
          return (datatree.AddLog(arg1));
          break;
        case oLog10:
          return (datatree.AddLog10(arg1));
          break;
        default:
          cerr << "BinaryOpNode::Compute_RHS: case not handled";
          exit(EXIT_FAILURE);
        }
      break;
    case 1: /*Binary Operator*/
      switch (op)
        {
        case oPlus:
          return (datatree.AddPlus(arg1, arg2));
          break;
        case oMinus:
          return (datatree.AddMinus(arg1, arg2));
          break;
        case oTimes:
          return (datatree.AddTimes(arg1, arg2));
          break;
        case oDivide:
          return (datatree.AddDivide(arg1, arg2));
          break;
        case oPower:
          return (datatree.AddPower(arg1, arg2));
          break;
        default:
          cerr << "BinaryOpNode::Compute_RHS: case not handled";
          exit(EXIT_FAILURE);
        }
      break;
    }
  return ((expr_t) NULL);
}

pair<int, expr_t>
BinaryOpNode::normalizeEquation(int var_endo, vector<pair<int, pair<expr_t, expr_t> > > &List_of_Op_RHS) const
{
  /* Checks if the current value of the endogenous variable related to the equation
     is present in the arguments of the binary operator. */
  vector<pair<int, pair<expr_t, expr_t> > > List_of_Op_RHS1, List_of_Op_RHS2;
  int is_endogenous_present_1, is_endogenous_present_2;
  pair<int, expr_t> res;
  expr_t expr_t_1, expr_t_2;
  res = arg1->normalizeEquation(var_endo, List_of_Op_RHS1);
  is_endogenous_present_1 = res.first;
  expr_t_1 = res.second;

  res = arg2->normalizeEquation(var_endo, List_of_Op_RHS2);
  is_endogenous_present_2 = res.first;
  expr_t_2 = res.second;

  /* If the two expressions contains the current value of the endogenous variable associated to the equation
     the equation could not be normalized and the process is given-up.*/
  if (is_endogenous_present_1 == 2 || is_endogenous_present_2 == 2)
    return (make_pair(2, (expr_t) NULL));
  else if (is_endogenous_present_1 && is_endogenous_present_2)
    return (make_pair(2, (expr_t) NULL));
  else if (is_endogenous_present_1) /*If the current values of the endogenous variable associated to the equation
                                      is present only in the first operand of the expression, we try to normalize the equation*/
    {
      if (op_code == oEqual)       /* The end of the normalization process :
                                      All the operations needed to normalize the equation are applied. */
        {
          pair<int, pair<expr_t, expr_t> > it;
          int oo = List_of_Op_RHS1.size();
          for (int i = 0; i < oo; i++)
            {
              it = List_of_Op_RHS1.back();
              List_of_Op_RHS1.pop_back();
              if (it.second.first && !it.second.second) /*Binary operator*/
                expr_t_2 = Compute_RHS(expr_t_2, (BinaryOpNode *) it.second.first, it.first, 1);
              else if (it.second.second && !it.second.first) /*Binary operator*/
                expr_t_2 = Compute_RHS(it.second.second, expr_t_2, it.first, 1);
              else if (it.second.second && it.second.first) /*Binary operator*/
                expr_t_2 = Compute_RHS(it.second.first, it.second.second, it.first, 1);
              else                                                                                 /*Unary operator*/
                expr_t_2 = Compute_RHS((UnaryOpNode *) expr_t_2, (UnaryOpNode *) it.second.first, it.first, 0);
            }
        }
      else
        List_of_Op_RHS = List_of_Op_RHS1;
    }
  else if (is_endogenous_present_2)
    {
      if (op_code == oEqual)
        {
          int oo = List_of_Op_RHS2.size();
          for (int i = 0; i < oo; i++)
            {
              pair<int, pair<expr_t, expr_t> > it;
              it = List_of_Op_RHS2.back();
              List_of_Op_RHS2.pop_back();
              if (it.second.first && !it.second.second) /*Binary operator*/
                expr_t_1 = Compute_RHS((BinaryOpNode *) expr_t_1, (BinaryOpNode *) it.second.first, it.first, 1);
              else if (it.second.second && !it.second.first) /*Binary operator*/
                expr_t_1 = Compute_RHS((BinaryOpNode *) it.second.second, (BinaryOpNode *) expr_t_1, it.first, 1);
              else if (it.second.second && it.second.first) /*Binary operator*/
                expr_t_1 = Compute_RHS(it.second.first, it.second.second, it.first, 1);
              else
                expr_t_1 = Compute_RHS((UnaryOpNode *) expr_t_1, (UnaryOpNode *) it.second.first, it.first, 0);
            }
        }
      else
        List_of_Op_RHS = List_of_Op_RHS2;
    }
  switch (op_code)
    {
    case oPlus:
      if (!is_endogenous_present_1 && !is_endogenous_present_2)
        {
          List_of_Op_RHS.push_back(make_pair(oMinus, make_pair(datatree.AddPlus(expr_t_1, expr_t_2), (expr_t) NULL)));
          return (make_pair(0, datatree.AddPlus(expr_t_1, expr_t_2)));
        }
      else if (is_endogenous_present_1 && is_endogenous_present_2)
        return (make_pair(1, (expr_t) NULL));
      else if (!is_endogenous_present_1 && is_endogenous_present_2)
        {
          List_of_Op_RHS.push_back(make_pair(oMinus, make_pair(expr_t_1, (expr_t) NULL)));
          return (make_pair(1, expr_t_1));
        }
      else if (is_endogenous_present_1 && !is_endogenous_present_2)
        {
          List_of_Op_RHS.push_back(make_pair(oMinus, make_pair(expr_t_2, (expr_t) NULL)));
          return (make_pair(1, expr_t_2));
        }
      break;
    case oMinus:
      if (!is_endogenous_present_1 && !is_endogenous_present_2)
        {
          List_of_Op_RHS.push_back(make_pair(oMinus, make_pair(datatree.AddMinus(expr_t_1, expr_t_2), (expr_t) NULL)));
          return (make_pair(0, datatree.AddMinus(expr_t_1, expr_t_2)));
        }
      else if (is_endogenous_present_1 && is_endogenous_present_2)
        return (make_pair(1, (expr_t) NULL));
      else if (!is_endogenous_present_1 && is_endogenous_present_2)
        {
          List_of_Op_RHS.push_back(make_pair(oUminus, make_pair((expr_t) NULL, (expr_t) NULL)));
          List_of_Op_RHS.push_back(make_pair(oMinus, make_pair(expr_t_1, (expr_t) NULL)));
          return (make_pair(1, expr_t_1));
        }
      else if (is_endogenous_present_1 && !is_endogenous_present_2)
        {
          List_of_Op_RHS.push_back(make_pair(oPlus, make_pair(expr_t_2, (expr_t) NULL)));
          return (make_pair(1, datatree.AddUMinus(expr_t_2)));
        }
      break;
    case oTimes:
      if (!is_endogenous_present_1 && !is_endogenous_present_2)
        return (make_pair(0, datatree.AddTimes(expr_t_1, expr_t_2)));
      else if (!is_endogenous_present_1 && is_endogenous_present_2)
        {
          List_of_Op_RHS.push_back(make_pair(oDivide, make_pair(expr_t_1, (expr_t) NULL)));
          return (make_pair(1, expr_t_1));
        }
      else if (is_endogenous_present_1 && !is_endogenous_present_2)
        {
          List_of_Op_RHS.push_back(make_pair(oDivide, make_pair(expr_t_2, (expr_t) NULL)));
          return (make_pair(1, expr_t_2));
        }
      else
        return (make_pair(1, (expr_t) NULL));
      break;
    case oDivide:
      if (!is_endogenous_present_1 && !is_endogenous_present_2)
        return (make_pair(0, datatree.AddDivide(expr_t_1, expr_t_2)));
      else if (!is_endogenous_present_1 && is_endogenous_present_2)
        {
          List_of_Op_RHS.push_back(make_pair(oDivide, make_pair((expr_t) NULL, expr_t_1)));
          return (make_pair(1, expr_t_1));
        }
      else if (is_endogenous_present_1 && !is_endogenous_present_2)
        {
          List_of_Op_RHS.push_back(make_pair(oTimes, make_pair(expr_t_2, (expr_t) NULL)));
          return (make_pair(1, expr_t_2));
        }
      else
        return (make_pair(1, (expr_t) NULL));
      break;
    case oPower:
      if (!is_endogenous_present_1 && !is_endogenous_present_2)
        return (make_pair(0, datatree.AddPower(expr_t_1, expr_t_2)));
      else if (is_endogenous_present_1 && !is_endogenous_present_2)
        {
          List_of_Op_RHS.push_back(make_pair(oPower, make_pair(datatree.AddDivide(datatree.One, expr_t_2), (expr_t) NULL)));
          return (make_pair(1, (expr_t) NULL));
        }
      else if (!is_endogenous_present_1 && is_endogenous_present_2)
        {
          /* we have to nomalize a^f(X) = RHS */
          /* First computes the ln(RHS)*/
          List_of_Op_RHS.push_back(make_pair(oLog, make_pair((expr_t) NULL, (expr_t) NULL)));
          /* Second  computes f(X) = ln(RHS) / ln(a)*/
          List_of_Op_RHS.push_back(make_pair(oDivide, make_pair((expr_t) NULL, datatree.AddLog(expr_t_1))));
          return (make_pair(1, (expr_t) NULL));
        }
      break;
    case oEqual:
      if (!is_endogenous_present_1 && !is_endogenous_present_2)
        {
          return (make_pair(0,
                            datatree.AddEqual(datatree.AddVariable(datatree.symbol_table.getID(eEndogenous, var_endo), 0), datatree.AddMinus(expr_t_2, expr_t_1))
                            ));
        }
      else if (is_endogenous_present_1 && is_endogenous_present_2)
        {
          return (make_pair(0,
                            datatree.AddEqual(datatree.AddVariable(datatree.symbol_table.getID(eEndogenous, var_endo), 0), datatree.Zero)
                            ));
        }
      else if (!is_endogenous_present_1 && is_endogenous_present_2)
        {
          return (make_pair(0,
                            datatree.AddEqual(datatree.AddVariable(datatree.symbol_table.getID(eEndogenous, var_endo), 0), /*datatree.AddUMinus(expr_t_1)*/ expr_t_1)
                            ));
        }
      else if (is_endogenous_present_1 && !is_endogenous_present_2)
        {
          return (make_pair(0,
                            datatree.AddEqual(datatree.AddVariable(datatree.symbol_table.getID(eEndogenous, var_endo), 0), expr_t_2)
                            ));
        }
      break;
    case oMax:
      if (!is_endogenous_present_1 && !is_endogenous_present_2)
        return (make_pair(0, datatree.AddMax(expr_t_1, expr_t_2)));
      else
        return (make_pair(1, (expr_t) NULL));
      break;
    case oMin:
      if (!is_endogenous_present_1 && !is_endogenous_present_2)
        return (make_pair(0, datatree.AddMin(expr_t_1, expr_t_2)));
      else
        return (make_pair(1, (expr_t) NULL));
      break;
    case oLess:
      if (!is_endogenous_present_1 && !is_endogenous_present_2)
        return (make_pair(0, datatree.AddLess(expr_t_1, expr_t_2)));
      else
        return (make_pair(1, (expr_t) NULL));
      break;
    case oGreater:
      if (!is_endogenous_present_1 && !is_endogenous_present_2)
        return (make_pair(0, datatree.AddGreater(expr_t_1, expr_t_2)));
      else
        return (make_pair(1, (expr_t) NULL));
      break;
    case oLessEqual:
      if (!is_endogenous_present_1 && !is_endogenous_present_2)
        return (make_pair(0, datatree.AddLessEqual(expr_t_1, expr_t_2)));
      else
        return (make_pair(1, (expr_t) NULL));
      break;
    case oGreaterEqual:
      if (!is_endogenous_present_1 && !is_endogenous_present_2)
        return (make_pair(0, datatree.AddGreaterEqual(expr_t_1, expr_t_2)));
      else
        return (make_pair(1, (expr_t) NULL));
      break;
    case oEqualEqual:
      if (!is_endogenous_present_1 && !is_endogenous_present_2)
        return (make_pair(0, datatree.AddEqualEqual(expr_t_1, expr_t_2)));
      else
        return (make_pair(1, (expr_t) NULL));
      break;
    case oDifferent:
      if (!is_endogenous_present_1 && !is_endogenous_present_2)
        return (make_pair(0, datatree.AddDifferent(expr_t_1, expr_t_2)));
      else
        return (make_pair(1, (expr_t) NULL));
      break;
    default:
      cerr << "Binary operator not handled during the normalization process" << endl;
      return (make_pair(2, (expr_t) NULL)); // Could not be normalized
    }
  // Suppress GCC warning
  cerr << "BinaryOpNode::normalizeEquation: impossible case" << endl;
  exit(EXIT_FAILURE);
}

expr_t
BinaryOpNode::getChainRuleDerivative(int deriv_id, const map<int, expr_t> &recursive_variables)
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
    case oPlus:
      return alt_datatree.AddPlus(alt_arg1, alt_arg2);
    case oMinus:
      return alt_datatree.AddMinus(alt_arg1, alt_arg2);
    case oTimes:
      return alt_datatree.AddTimes(alt_arg1, alt_arg2);
    case oDivide:
      return alt_datatree.AddDivide(alt_arg1, alt_arg2);
    case oPower:
      return alt_datatree.AddPower(alt_arg1, alt_arg2);
    case oEqual:
      return alt_datatree.AddEqual(alt_arg1, alt_arg2);
    case oMax:
      return alt_datatree.AddMax(alt_arg1, alt_arg2);
    case oMin:
      return alt_datatree.AddMin(alt_arg1, alt_arg2);
    case oLess:
      return alt_datatree.AddLess(alt_arg1, alt_arg2);
    case oGreater:
      return alt_datatree.AddGreater(alt_arg1, alt_arg2);
    case oLessEqual:
      return alt_datatree.AddLessEqual(alt_arg1, alt_arg2);
    case oGreaterEqual:
      return alt_datatree.AddGreaterEqual(alt_arg1, alt_arg2);
    case oEqualEqual:
      return alt_datatree.AddEqualEqual(alt_arg1, alt_arg2);
    case oDifferent:
      return alt_datatree.AddDifferent(alt_arg1, alt_arg2);
    case oPowerDeriv:
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
BinaryOpNode::cloneDynamic(DataTree &dynamic_datatree) const
{
  expr_t substarg1 = arg1->cloneDynamic(dynamic_datatree);
  expr_t substarg2 = arg2->cloneDynamic(dynamic_datatree);
  return buildSimilarBinaryOpNode(substarg1, substarg2, dynamic_datatree);
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

expr_t
BinaryOpNode::undiff() const
{
  expr_t arg1subst = arg1->undiff();
  expr_t arg2subst = arg2->undiff();
  return buildSimilarBinaryOpNode(arg1subst, arg2subst, datatree);
}

int
BinaryOpNode::PacMaxLag(vector<int> &lhs) const
{
  return max(arg1->PacMaxLag(lhs), arg2->PacMaxLag(lhs));
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
    {
      switch (op_code)
        {
        case oPlus:
        case oMinus:
        case oEqual:
          arg1subst = maxendolead1 >= 2 ? arg1->substituteEndoLeadGreaterThanTwo(subst_table, neweqs, deterministic_model) : arg1;
          arg2subst = maxendolead2 >= 2 ? arg2->substituteEndoLeadGreaterThanTwo(subst_table, neweqs, deterministic_model) : arg2;
          return buildSimilarBinaryOpNode(arg1subst, arg2subst, datatree);
        case oTimes:
        case oDivide:
          if (maxendolead1 >= 2 && maxendolead2 == 0 && arg2->maxExoLead() == 0)
            {
              arg1subst = arg1->substituteEndoLeadGreaterThanTwo(subst_table, neweqs, deterministic_model);
              return buildSimilarBinaryOpNode(arg1subst, arg2, datatree);
            }
          if (maxendolead1 == 0 && arg1->maxExoLead() == 0
              && maxendolead2 >= 2 && op_code == oTimes)
            {
              arg2subst = arg2->substituteEndoLeadGreaterThanTwo(subst_table, neweqs, deterministic_model);
              return buildSimilarBinaryOpNode(arg1, arg2subst, datatree);
            }
          return createEndoLeadAuxiliaryVarForMyself(subst_table, neweqs);
        default:
          return createEndoLeadAuxiliaryVarForMyself(subst_table, neweqs);
        }
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
    {
      switch (op_code)
        {
        case oPlus:
        case oMinus:
        case oEqual:
          arg1subst = maxexolead1 >= 1 ? arg1->substituteExoLead(subst_table, neweqs, deterministic_model) : arg1;
          arg2subst = maxexolead2 >= 1 ? arg2->substituteExoLead(subst_table, neweqs, deterministic_model) : arg2;
          return buildSimilarBinaryOpNode(arg1subst, arg2subst, datatree);
        case oTimes:
        case oDivide:
          if (maxexolead1 >= 1 && maxexolead2 == 0 && arg2->maxEndoLead() == 0)
            {
              arg1subst = arg1->substituteExoLead(subst_table, neweqs, deterministic_model);
              return buildSimilarBinaryOpNode(arg1subst, arg2, datatree);
            }
          if (maxexolead1 == 0 && arg1->maxEndoLead() == 0
              && maxexolead2 >= 1 && op_code == oTimes)
            {
              arg2subst = arg2->substituteExoLead(subst_table, neweqs, deterministic_model);
              return buildSimilarBinaryOpNode(arg1, arg2subst, datatree);
            }
          return createExoLeadAuxiliaryVarForMyself(subst_table, neweqs);
        default:
          return createExoLeadAuxiliaryVarForMyself(subst_table, neweqs);
        }
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

void
BinaryOpNode::findDiffUnaryOpNodes(DataTree &static_datatree, set<UnaryOpNode *> &unary_op_nodes_in_diff) const
{
  arg1->findDiffUnaryOpNodes(static_datatree, unary_op_nodes_in_diff);
  arg2->findDiffUnaryOpNodes(static_datatree, unary_op_nodes_in_diff);
}

void
BinaryOpNode::findDiffNodes(DataTree &static_datatree, diff_table_t &diff_table) const
{
  arg1->findDiffNodes(static_datatree, diff_table);
  arg2->findDiffNodes(static_datatree, diff_table);
}

expr_t
BinaryOpNode::substituteDiff(DataTree &static_datatree, diff_table_t &diff_table, subst_table_t &subst_table,
                             vector<BinaryOpNode *> &neweqs) const
{
  expr_t arg1subst = arg1->substituteDiff(static_datatree, diff_table, subst_table, neweqs);
  expr_t arg2subst = arg2->substituteDiff(static_datatree, diff_table, subst_table, neweqs);
  return buildSimilarBinaryOpNode(arg1subst, arg2subst, datatree);
}

expr_t
BinaryOpNode::substituteDiffUnaryOpNodes(DataTree &static_datatree, set<UnaryOpNode *> &nodes, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  expr_t arg1subst = arg1->substituteDiffUnaryOpNodes(static_datatree, nodes, subst_table, neweqs);
  expr_t arg2subst = arg2->substituteDiffUnaryOpNodes(static_datatree, nodes, subst_table, neweqs);
  return buildSimilarBinaryOpNode(arg1subst, arg2subst, datatree);
}

bool
BinaryOpNode::isDiffPresent() const
{
  return arg1->isDiffPresent() || arg2->isDiffPresent();
}

expr_t
BinaryOpNode::substitutePacExpectation(map<const PacExpectationNode *, const BinaryOpNode *> &subst_table)
{
  expr_t arg1subst = arg1->substitutePacExpectation(subst_table);
  expr_t arg2subst = arg2->substitutePacExpectation(subst_table);
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
BinaryOpNode::containsEndogenous(void) const
{
  return (arg1->containsEndogenous() || arg2->containsEndogenous());
}

bool
BinaryOpNode::containsExogenous() const
{
  return (arg1->containsExogenous() || arg2->containsExogenous());
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
BinaryOpNode::removeTrendLeadLag(map<int, expr_t> trend_symbols_map) const
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

void
BinaryOpNode::setVarExpectationIndex(map<string, pair<SymbolList, int> > &var_model_info)
{
  arg1->setVarExpectationIndex(var_model_info);
  arg2->setVarExpectationIndex(var_model_info);
}

void
BinaryOpNode::walkPacParametersHelper(const expr_t arg1, const expr_t arg2,
                                      pair<int, int> &lhs,
                                      set<pair<int, pair<int, int> > > &ec_params_and_vars,
                                      set<pair<int, pair<int, int> > > &ar_params_and_vars) const
{
  set<int> params;
  arg1->collectVariables(eParameter, params);
  if (params.size() != 1)
    return;

  set<pair<int, int> > endogs;
  arg2->collectDynamicVariables(eEndogenous, endogs);
  if (endogs.size() == 1)
    ar_params_and_vars.insert(make_pair(*(params.begin()), *(endogs.begin())));
  else if (endogs.size() >= 2)
    {
      BinaryOpNode *testarg2 = dynamic_cast<BinaryOpNode *>(arg2);
      if (testarg2 != NULL && testarg2->get_op_code() == oMinus)
        {
          VariableNode *test_arg1 = dynamic_cast<VariableNode *>(testarg2->get_arg1());
          VariableNode *test_arg2 = dynamic_cast<VariableNode *>(testarg2->get_arg2());
          if (test_arg1 != NULL && test_arg2 != NULL && lhs.first != -1)
            {
              test_arg1->collectDynamicVariables(eEndogenous, endogs);
              ec_params_and_vars.insert(make_pair(*(params.begin()), *(endogs.begin())));
              endogs.clear();
              test_arg2->collectDynamicVariables(eEndogenous, endogs);
              ec_params_and_vars.insert(make_pair(*(params.begin()), *(endogs.begin())));
            }
        }
    }
}

void
BinaryOpNode::walkPacParameters(bool &pac_encountered, pair<int, int> &lhs, set<pair<int, pair<int, int> > > &ec_params_and_vars, set<pair<int, pair<int, int> > > &ar_params_and_vars) const
{
  if (op_code == oTimes)
    {
      int orig_ar_params_and_vars_size = ar_params_and_vars.size();
      int orig_ec_params_and_vars_size = ec_params_and_vars.size();
      walkPacParametersHelper(arg1, arg2, lhs, ec_params_and_vars, ar_params_and_vars);
      if ((int)ar_params_and_vars.size() == orig_ar_params_and_vars_size
          && (int)ec_params_and_vars.size() == orig_ec_params_and_vars_size)
        walkPacParametersHelper(arg2, arg1, lhs, ec_params_and_vars, ar_params_and_vars);
    }
  else if (op_code == oEqual)
    {
      set<pair<int, int> > general_lhs;
      arg1->collectDynamicVariables(eEndogenous, general_lhs);
      if (general_lhs.size() == 1)
        lhs = *(general_lhs.begin());
    }

  arg1->walkPacParameters(pac_encountered, lhs, ec_params_and_vars, ar_params_and_vars);
  arg2->walkPacParameters(pac_encountered, lhs, ec_params_and_vars, ar_params_and_vars);
}

void
BinaryOpNode::addParamInfoToPac(pair<int, int> &lhs_arg, set<pair<int, pair<int, int> > > &ec_params_and_vars_arg, set<pair<int, pair<int, int> > > &ar_params_and_vars_arg)
{
  arg1->addParamInfoToPac(lhs_arg, ec_params_and_vars_arg, ar_params_and_vars_arg);
  arg2->addParamInfoToPac(lhs_arg, ec_params_and_vars_arg, ar_params_and_vars_arg);
}

void
BinaryOpNode::fillPacExpectationVarInfo(string &model_name_arg, vector<int> &lhs_arg, int max_lag_arg, vector<bool> &nonstationary_arg, int growth_symb_id_arg, int equation_number_arg)
{
  arg1->fillPacExpectationVarInfo(model_name_arg, lhs_arg, max_lag_arg, nonstationary_arg, growth_symb_id_arg, equation_number_arg);
  arg2->fillPacExpectationVarInfo(model_name_arg, lhs_arg, max_lag_arg, nonstationary_arg, growth_symb_id_arg, equation_number_arg);
}

bool
BinaryOpNode::isVarModelReferenced(const string &model_info_name) const
{
  return arg1->isVarModelReferenced(model_info_name)
    || arg2->isVarModelReferenced(model_info_name);
}

void
BinaryOpNode::getEndosAndMaxLags(map<string, int> &model_endos_and_lags) const
{
  arg1->getEndosAndMaxLags(model_endos_and_lags);
  arg2->getEndosAndMaxLags(model_endos_and_lags);
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

TrinaryOpNode::TrinaryOpNode(DataTree &datatree_arg, const expr_t arg1_arg,
                             TrinaryOpcode op_code_arg, const expr_t arg2_arg, const expr_t arg3_arg) :
  ExprNode(datatree_arg),
  arg1(arg1_arg),
  arg2(arg2_arg),
  arg3(arg3_arg),
  op_code(op_code_arg)
{
  datatree.trinary_op_node_map[make_pair(make_pair(make_pair(arg1, arg2), arg3), op_code)] = this;
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
    case oNormcdf:
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
    case oNormpdf:
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
  temporary_terms_t::const_iterator it = temporary_terms.find(const_cast<TrinaryOpNode *>(this));
  // A temporary term behaves as a variable
  if (it != temporary_terms.end())
    return 100;

  switch (op_code)
    {
    case oNormcdf:
    case oNormpdf:
      return 100;
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

int
TrinaryOpNode::cost(const map<NodeTreeReference, temporary_terms_t> &temp_terms_map, bool is_matlab) const
{
  // For a temporary term, the cost is null
  for (map<NodeTreeReference, temporary_terms_t>::const_iterator it = temp_terms_map.begin();
       it != temp_terms_map.end(); it++)
    if (it->second.find(const_cast<TrinaryOpNode *>(this)) != it->second.end())
      return 0;

  int arg_cost = arg1->cost(temp_terms_map, is_matlab)
    + arg2->cost(temp_terms_map, is_matlab)
    + arg3->cost(temp_terms_map, is_matlab);

  return cost(arg_cost, is_matlab);
}

int
TrinaryOpNode::cost(const temporary_terms_t &temporary_terms, bool is_matlab) const
{
  // For a temporary term, the cost is null
  if (temporary_terms.find(const_cast<TrinaryOpNode *>(this)) != temporary_terms.end())
    return 0;

  int arg_cost = arg1->cost(temporary_terms, is_matlab)
    + arg2->cost(temporary_terms, is_matlab)
    + arg3->cost(temporary_terms, is_matlab);

  return cost(arg_cost, is_matlab);
}

int
TrinaryOpNode::cost(int cost, bool is_matlab) const
{
  if (is_matlab)
    // Cost for Matlab files
    switch (op_code)
      {
      case oNormcdf:
      case oNormpdf:
        return cost+1000;
      }
  else
    // Cost for C files
    switch (op_code)
      {
      case oNormcdf:
      case oNormpdf:
        return cost+1000;
      }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

void
TrinaryOpNode::computeTemporaryTerms(map<expr_t, pair<int, NodeTreeReference> > &reference_count,
                                     map<NodeTreeReference, temporary_terms_t> &temp_terms_map,
                                     bool is_matlab, NodeTreeReference tr) const
{
  expr_t this2 = const_cast<TrinaryOpNode *>(this);
  map<expr_t, pair<int, NodeTreeReference > >::iterator it = reference_count.find(this2);
  if (it == reference_count.end())
    {
      // If this node has never been encountered, set its ref count to one,
      //  and travel through its children
      reference_count[this2] = make_pair(1, tr);
      arg1->computeTemporaryTerms(reference_count, temp_terms_map, is_matlab, tr);
      arg2->computeTemporaryTerms(reference_count, temp_terms_map, is_matlab, tr);
      arg3->computeTemporaryTerms(reference_count, temp_terms_map, is_matlab, tr);
    }
  else
    {
      // If the node has already been encountered, increment its ref count
      //  and declare it as a temporary term if it is too costly
      reference_count[this2] = make_pair(it->second.first + 1, it->second.second);;
      if (reference_count[this2].first * cost(temp_terms_map, is_matlab) > MIN_COST(is_matlab))
        temp_terms_map[reference_count[this2].second].insert(this2);
    }
}

void
TrinaryOpNode::computeTemporaryTerms(map<expr_t, int> &reference_count,
                                     temporary_terms_t &temporary_terms,
                                     map<expr_t, pair<int, int> > &first_occurence,
                                     int Curr_block,
                                     vector<vector<temporary_terms_t> > &v_temporary_terms,
                                     int equation) const
{
  expr_t this2 = const_cast<TrinaryOpNode *>(this);
  map<expr_t, int>::iterator it = reference_count.find(this2);
  if (it == reference_count.end())
    {
      reference_count[this2] = 1;
      first_occurence[this2] = make_pair(Curr_block, equation);
      arg1->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, Curr_block, v_temporary_terms, equation);
      arg2->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, Curr_block, v_temporary_terms, equation);
      arg3->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, Curr_block, v_temporary_terms, equation);
    }
  else
    {
      reference_count[this2]++;
      if (reference_count[this2] * cost(temporary_terms, false) > MIN_COST_C)
        {
          temporary_terms.insert(this2);
          v_temporary_terms[first_occurence[this2].first][first_occurence[this2].second].insert(this2);
        }
    }
}

double
TrinaryOpNode::eval_opcode(double v1, TrinaryOpcode op_code, double v2, double v3) throw (EvalException, EvalExternalFunctionException)
{
  switch (op_code)
    {
    case oNormcdf:
      return (0.5*(1+erf((v1-v2)/v3/M_SQRT2)));
    case oNormpdf:
      return (1/(v3*sqrt(2*M_PI)*exp(pow((v1-v2)/v3, 2)/2)));
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

double
TrinaryOpNode::eval(const eval_context_t &eval_context) const throw (EvalException, EvalExternalFunctionException)
{
  double v1 = arg1->eval(eval_context);
  double v2 = arg2->eval(eval_context);
  double v3 = arg3->eval(eval_context);

  return eval_opcode(v1, op_code, v2, v3);
}

void
TrinaryOpNode::compile(ostream &CompileCode, unsigned int &instruction_number,
                       bool lhs_rhs, const temporary_terms_t &temporary_terms,
                       const map_idx_t &map_idx, bool dynamic, bool steady_dynamic,
                       deriv_node_temp_terms_t &tef_terms) const
{
  // If current node is a temporary term
  temporary_terms_t::const_iterator it = temporary_terms.find(const_cast<TrinaryOpNode *>(this));
  if (it != temporary_terms.end())
    {
      if (dynamic)
        {
          map_idx_t::const_iterator ii = map_idx.find(idx);
          FLDT_ fldt(ii->second);
          fldt.write(CompileCode, instruction_number);
        }
      else
        {
          map_idx_t::const_iterator ii = map_idx.find(idx);
          FLDST_ fldst(ii->second);
          fldst.write(CompileCode, instruction_number);
        }
      return;
    }
  arg1->compile(CompileCode, instruction_number, lhs_rhs, temporary_terms, map_idx, dynamic, steady_dynamic, tef_terms);
  arg2->compile(CompileCode, instruction_number, lhs_rhs, temporary_terms, map_idx, dynamic, steady_dynamic, tef_terms);
  arg3->compile(CompileCode, instruction_number, lhs_rhs, temporary_terms, map_idx, dynamic, steady_dynamic, tef_terms);
  FTRINARY_ ftrinary(op_code);
  ftrinary.write(CompileCode, instruction_number);
}

void
TrinaryOpNode::collectTemporary_terms(const temporary_terms_t &temporary_terms, temporary_terms_inuse_t &temporary_terms_inuse, int Curr_Block) const
{
  temporary_terms_t::const_iterator it = temporary_terms.find(const_cast<TrinaryOpNode *>(this));
  if (it != temporary_terms.end())
    temporary_terms_inuse.insert(idx);
  else
    {
      arg1->collectTemporary_terms(temporary_terms, temporary_terms_inuse, Curr_Block);
      arg2->collectTemporary_terms(temporary_terms, temporary_terms_inuse, Curr_Block);
      arg3->collectTemporary_terms(temporary_terms, temporary_terms_inuse, Curr_Block);
    }
}

bool
TrinaryOpNode::containsExternalFunction() const
{
  return arg1->containsExternalFunction()
    || arg2->containsExternalFunction()
    || arg3->containsExternalFunction();
}

void
TrinaryOpNode::writeJsonOutput(ostream &output,
                               const temporary_terms_t &temporary_terms,
                               deriv_node_temp_terms_t &tef_terms,
                               const bool isdynamic) const
{
  // If current node is a temporary term
  temporary_terms_t::const_iterator it = temporary_terms.find(const_cast<TrinaryOpNode *>(this));
  if (it != temporary_terms.end())
    {
      output << "T" << idx;
      return;
    }

  switch (op_code)
    {
    case oNormcdf:
      output << "normcdf(";
      break;
    case oNormpdf:
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
                           deriv_node_temp_terms_t &tef_terms) const
{
  // If current node is a temporary term
  temporary_terms_t::const_iterator it = temporary_terms.find(const_cast<TrinaryOpNode *>(this));
  if (it != temporary_terms.end())
    {
      output << "T" << idx;
      return;
    }

  switch (op_code)
    {
    case oNormcdf:
      if (IS_C(output_type))
        {
          // In C, there is no normcdf() primitive, so use erf()
          output << "(0.5*(1+erf(((";
          arg1->writeOutput(output, output_type, temporary_terms, tef_terms);
          output << ")-(";
          arg2->writeOutput(output, output_type, temporary_terms, tef_terms);
          output << "))/(";
          arg3->writeOutput(output, output_type, temporary_terms, tef_terms);
          output << ")/M_SQRT2)))";
        }
      else
        {
          output << "normcdf(";
          arg1->writeOutput(output, output_type, temporary_terms, tef_terms);
          output << ",";
          arg2->writeOutput(output, output_type, temporary_terms, tef_terms);
          output << ",";
          arg3->writeOutput(output, output_type, temporary_terms, tef_terms);
          output << ")";
        }
      break;
    case oNormpdf:
      if (IS_C(output_type))
        {
          //(1/(v3*sqrt(2*M_PI)*exp(pow((v1-v2)/v3,2)/2)))
          output << "(1/(";
          arg3->writeOutput(output, output_type, temporary_terms, tef_terms);
          output << "*sqrt(2*M_PI)*exp(pow((";
          arg1->writeOutput(output, output_type, temporary_terms, tef_terms);
          output << "-";
          arg2->writeOutput(output, output_type, temporary_terms, tef_terms);
          output << ")/";
          arg3->writeOutput(output, output_type, temporary_terms, tef_terms);
          output << ",2)/2)))";
        }
      else
        {
          output << "normpdf(";
          arg1->writeOutput(output, output_type, temporary_terms, tef_terms);
          output << ",";
          arg2->writeOutput(output, output_type, temporary_terms, tef_terms);
          output << ",";
          arg3->writeOutput(output, output_type, temporary_terms, tef_terms);
          output << ")";
        }
      break;
    }
}

void
TrinaryOpNode::writeExternalFunctionOutput(ostream &output, ExprNodeOutputType output_type,
                                           const temporary_terms_t &temporary_terms,
                                           deriv_node_temp_terms_t &tef_terms) const
{
  arg1->writeExternalFunctionOutput(output, output_type, temporary_terms, tef_terms);
  arg2->writeExternalFunctionOutput(output, output_type, temporary_terms, tef_terms);
  arg3->writeExternalFunctionOutput(output, output_type, temporary_terms, tef_terms);
}

void
TrinaryOpNode::writeJsonExternalFunctionOutput(vector<string> &efout,
                                               const temporary_terms_t &temporary_terms,
                                               deriv_node_temp_terms_t &tef_terms,
                                               const bool isdynamic) const
{
  arg1->writeJsonExternalFunctionOutput(efout, temporary_terms, tef_terms, isdynamic);
  arg2->writeJsonExternalFunctionOutput(efout, temporary_terms, tef_terms, isdynamic);
  arg3->writeJsonExternalFunctionOutput(efout, temporary_terms, tef_terms, isdynamic);
}

void
TrinaryOpNode::compileExternalFunctionOutput(ostream &CompileCode, unsigned int &instruction_number,
                                             bool lhs_rhs, const temporary_terms_t &temporary_terms,
                                             const map_idx_t &map_idx, bool dynamic, bool steady_dynamic,
                                             deriv_node_temp_terms_t &tef_terms) const
{
  arg1->compileExternalFunctionOutput(CompileCode, instruction_number, lhs_rhs, temporary_terms, map_idx,
                                      dynamic, steady_dynamic, tef_terms);
  arg2->compileExternalFunctionOutput(CompileCode, instruction_number, lhs_rhs, temporary_terms, map_idx,
                                      dynamic, steady_dynamic, tef_terms);
  arg3->compileExternalFunctionOutput(CompileCode, instruction_number, lhs_rhs, temporary_terms, map_idx,
                                      dynamic, steady_dynamic, tef_terms);
}

void
TrinaryOpNode::collectVARLHSVariable(set<expr_t> &result) const
{
  arg1->collectVARLHSVariable(result);
  arg2->collectVARLHSVariable(result);
  arg3->collectVARLHSVariable(result);
}

void
TrinaryOpNode::collectDynamicVariables(SymbolType type_arg, set<pair<int, int> > &result) const
{
  arg1->collectDynamicVariables(type_arg, result);
  arg2->collectDynamicVariables(type_arg, result);
  arg3->collectDynamicVariables(type_arg, result);
}

pair<int, expr_t>
TrinaryOpNode::normalizeEquation(int var_endo, vector<pair<int, pair<expr_t, expr_t> > > &List_of_Op_RHS) const
{
  pair<int, expr_t> res = arg1->normalizeEquation(var_endo, List_of_Op_RHS);
  bool is_endogenous_present_1 = res.first;
  expr_t expr_t_1 = res.second;
  res = arg2->normalizeEquation(var_endo, List_of_Op_RHS);
  bool is_endogenous_present_2 = res.first;
  expr_t expr_t_2 = res.second;
  res = arg3->normalizeEquation(var_endo, List_of_Op_RHS);
  bool is_endogenous_present_3 = res.first;
  expr_t expr_t_3 = res.second;
  if (!is_endogenous_present_1 && !is_endogenous_present_2 && !is_endogenous_present_3)
    return (make_pair(0, datatree.AddNormcdf(expr_t_1, expr_t_2, expr_t_3)));
  else
    return (make_pair(1, (expr_t) NULL));
}

expr_t
TrinaryOpNode::getChainRuleDerivative(int deriv_id, const map<int, expr_t> &recursive_variables)
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
    case oNormcdf:
      return alt_datatree.AddNormcdf(alt_arg1, alt_arg2, alt_arg3);
    case oNormpdf:
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
TrinaryOpNode::cloneDynamic(DataTree &dynamic_datatree) const
{
  expr_t substarg1 = arg1->cloneDynamic(dynamic_datatree);
  expr_t substarg2 = arg2->cloneDynamic(dynamic_datatree);
  expr_t substarg3 = arg3->cloneDynamic(dynamic_datatree);
  return buildSimilarTrinaryOpNode(substarg1, substarg2, substarg3, dynamic_datatree);
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

expr_t
TrinaryOpNode::undiff() const
{
  expr_t arg1subst = arg1->undiff();
  expr_t arg2subst = arg2->undiff();
  expr_t arg3subst = arg3->undiff();
  return buildSimilarTrinaryOpNode(arg1subst, arg2subst, arg3subst, datatree);
}

int
TrinaryOpNode::VarMinLag() const
{
  return min(min(arg1->VarMinLag(), arg2->VarMinLag()), arg3->VarMinLag());
}

void
TrinaryOpNode::VarMaxLag(DataTree &static_datatree, set<expr_t> &static_lhs, int &max_lag) const
{
  arg1->VarMaxLag(static_datatree, static_lhs, max_lag);
  arg2->VarMaxLag(static_datatree, static_lhs, max_lag);
  arg3->VarMaxLag(static_datatree, static_lhs, max_lag);
}

int
TrinaryOpNode::PacMaxLag(vector<int> &lhs) const
{
  return max(arg1->PacMaxLag(lhs), max(arg2->PacMaxLag(lhs), arg3->PacMaxLag(lhs)));
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

void
TrinaryOpNode::findDiffNodes(DataTree &static_datatree, diff_table_t &diff_table) const
{
  arg1->findDiffNodes(static_datatree, diff_table);
  arg2->findDiffNodes(static_datatree, diff_table);
  arg3->findDiffNodes(static_datatree, diff_table);
}

void
TrinaryOpNode::findDiffUnaryOpNodes(DataTree &static_datatree, set<UnaryOpNode *> &unary_op_nodes_in_diff) const
{
  arg1->findDiffUnaryOpNodes(static_datatree, unary_op_nodes_in_diff);
  arg2->findDiffUnaryOpNodes(static_datatree, unary_op_nodes_in_diff);
  arg3->findDiffUnaryOpNodes(static_datatree, unary_op_nodes_in_diff);
}

expr_t
TrinaryOpNode::substituteDiff(DataTree &static_datatree, diff_table_t &diff_table, subst_table_t &subst_table,
                              vector<BinaryOpNode *> &neweqs) const
{
  expr_t arg1subst = arg1->substituteDiff(static_datatree, diff_table, subst_table, neweqs);
  expr_t arg2subst = arg2->substituteDiff(static_datatree, diff_table, subst_table, neweqs);
  expr_t arg3subst = arg3->substituteDiff(static_datatree, diff_table, subst_table, neweqs);
  return buildSimilarTrinaryOpNode(arg1subst, arg2subst, arg3subst, datatree);
}

expr_t
TrinaryOpNode::substituteDiffUnaryOpNodes(DataTree &static_datatree, set<UnaryOpNode *> &nodes, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  expr_t arg1subst = arg1->substituteDiffUnaryOpNodes(static_datatree, nodes, subst_table, neweqs);
  expr_t arg2subst = arg2->substituteDiffUnaryOpNodes(static_datatree, nodes, subst_table, neweqs);
  expr_t arg3subst = arg3->substituteDiffUnaryOpNodes(static_datatree, nodes, subst_table, neweqs);
  return buildSimilarTrinaryOpNode(arg1subst, arg2subst, arg3subst, datatree);
}

bool
TrinaryOpNode::isDiffPresent() const
{
  return arg1->isDiffPresent() || arg2->isDiffPresent() || arg3->isDiffPresent();
}

expr_t
TrinaryOpNode::substitutePacExpectation(map<const PacExpectationNode *, const BinaryOpNode *> &subst_table)
{
  expr_t arg1subst = arg1->substitutePacExpectation(subst_table);
  expr_t arg2subst = arg2->substitutePacExpectation(subst_table);
  expr_t arg3subst = arg3->substitutePacExpectation(subst_table);
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
TrinaryOpNode::containsEndogenous(void) const
{
  return (arg1->containsEndogenous() || arg2->containsEndogenous() || arg3->containsEndogenous());
}

bool
TrinaryOpNode::containsExogenous() const
{
  return (arg1->containsExogenous() || arg2->containsExogenous() || arg3->containsExogenous());
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
TrinaryOpNode::removeTrendLeadLag(map<int, expr_t> trend_symbols_map) const
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

void
TrinaryOpNode::setVarExpectationIndex(map<string, pair<SymbolList, int> > &var_model_info)
{
  arg1->setVarExpectationIndex(var_model_info);
  arg2->setVarExpectationIndex(var_model_info);
  arg3->setVarExpectationIndex(var_model_info);
}

void
TrinaryOpNode::walkPacParameters(bool &pac_encountered, pair<int, int> &lhs, set<pair<int, pair<int, int> > > &ec_params_and_vars, set<pair<int, pair<int, int> > > &ar_params_and_vars) const
{
  arg1->walkPacParameters(pac_encountered, lhs, ec_params_and_vars, ar_params_and_vars);
  arg2->walkPacParameters(pac_encountered, lhs, ec_params_and_vars, ar_params_and_vars);
  arg3->walkPacParameters(pac_encountered, lhs, ec_params_and_vars, ar_params_and_vars);
}

void
TrinaryOpNode::addParamInfoToPac(pair<int, int> &lhs_arg, set<pair<int, pair<int, int> > > &ec_params_and_vars_arg, set<pair<int, pair<int, int> > > &ar_params_and_vars_arg)
{
  arg1->addParamInfoToPac(lhs_arg, ec_params_and_vars_arg, ar_params_and_vars_arg);
  arg2->addParamInfoToPac(lhs_arg, ec_params_and_vars_arg, ar_params_and_vars_arg);
  arg3->addParamInfoToPac(lhs_arg, ec_params_and_vars_arg, ar_params_and_vars_arg);
}

void
TrinaryOpNode::fillPacExpectationVarInfo(string &model_name_arg, vector<int> &lhs_arg, int max_lag_arg, vector<bool> &nonstationary_arg, int growth_symb_id_arg, int equation_number_arg)
{
  arg1->fillPacExpectationVarInfo(model_name_arg, lhs_arg, max_lag_arg, nonstationary_arg, growth_symb_id_arg, equation_number_arg);
  arg2->fillPacExpectationVarInfo(model_name_arg, lhs_arg, max_lag_arg, nonstationary_arg, growth_symb_id_arg, equation_number_arg);
  arg3->fillPacExpectationVarInfo(model_name_arg, lhs_arg, max_lag_arg, nonstationary_arg, growth_symb_id_arg, equation_number_arg);
}

bool
TrinaryOpNode::isVarModelReferenced(const string &model_info_name) const
{
  return arg1->isVarModelReferenced(model_info_name)
    || arg2->isVarModelReferenced(model_info_name)
    || arg3->isVarModelReferenced(model_info_name);
}

void
TrinaryOpNode::getEndosAndMaxLags(map<string, int> &model_endos_and_lags) const
{
  arg1->getEndosAndMaxLags(model_endos_and_lags);
  arg2->getEndosAndMaxLags(model_endos_and_lags);
  arg3->getEndosAndMaxLags(model_endos_and_lags);
}

expr_t
TrinaryOpNode::substituteStaticAuxiliaryVariable() const
{
  expr_t arg1subst = arg1->substituteStaticAuxiliaryVariable();
  expr_t arg2subst = arg2->substituteStaticAuxiliaryVariable();
  expr_t arg3subst = arg3->substituteStaticAuxiliaryVariable();
  return buildSimilarTrinaryOpNode(arg1subst, arg2subst, arg3subst, datatree);
}

AbstractExternalFunctionNode::AbstractExternalFunctionNode(DataTree &datatree_arg,
                                                           int symb_id_arg,
                                                           const vector<expr_t> &arguments_arg) :
  ExprNode(datatree_arg),
  symb_id(symb_id_arg),
  arguments(arguments_arg)
{
}

void
AbstractExternalFunctionNode::prepareForDerivation()
{
  if (preparedForDerivation)
    return;

  for (vector<expr_t>::const_iterator it = arguments.begin(); it != arguments.end(); it++)
    (*it)->prepareForDerivation();

  non_null_derivatives = arguments.at(0)->non_null_derivatives;
  for (int i = 1; i < (int) arguments.size(); i++)
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
  for (vector<expr_t>::const_iterator it = arguments.begin(); it != arguments.end(); it++)
    dargs.push_back((*it)->getDerivative(deriv_id));
  return composeDerivatives(dargs);
}

expr_t
AbstractExternalFunctionNode::getChainRuleDerivative(int deriv_id, const map<int, expr_t> &recursive_variables)
{
  assert(datatree.external_functions_table.getNargs(symb_id) > 0);
  vector<expr_t> dargs;
  for (vector<expr_t>::const_iterator it = arguments.begin();
       it != arguments.end(); it++)
    dargs.push_back((*it)->getChainRuleDerivative(deriv_id, recursive_variables));
  return composeDerivatives(dargs);
}

unsigned int
AbstractExternalFunctionNode::compileExternalFunctionArguments(ostream &CompileCode, unsigned int &instruction_number,
                                                               bool lhs_rhs, const temporary_terms_t &temporary_terms,
                                                               const map_idx_t &map_idx, bool dynamic, bool steady_dynamic,
                                                               deriv_node_temp_terms_t &tef_terms) const
{
  for (vector<expr_t>::const_iterator it = arguments.begin();
       it != arguments.end(); it++)
    (*it)->compile(CompileCode, instruction_number, lhs_rhs, temporary_terms, map_idx,
                   dynamic, steady_dynamic, tef_terms);
  return (arguments.size());
}

void
AbstractExternalFunctionNode::collectVARLHSVariable(set<expr_t> &result) const
{
   for (vector<expr_t>::const_iterator it = arguments.begin();
        it != arguments.end(); it++)
     (*it)->collectVARLHSVariable(result);
}

void
AbstractExternalFunctionNode::collectDynamicVariables(SymbolType type_arg, set<pair<int, int> > &result) const
{
  for (vector<expr_t>::const_iterator it = arguments.begin();
       it != arguments.end(); it++)
    (*it)->collectDynamicVariables(type_arg, result);
}

void
AbstractExternalFunctionNode::collectTemporary_terms(const temporary_terms_t &temporary_terms, temporary_terms_inuse_t &temporary_terms_inuse, int Curr_Block) const
{
  temporary_terms_t::const_iterator it = temporary_terms.find(const_cast<AbstractExternalFunctionNode *>(this));
  if (it != temporary_terms.end())
    temporary_terms_inuse.insert(idx);
  else
    {
      for (vector<expr_t>::const_iterator it = arguments.begin();
           it != arguments.end(); it++)
        (*it)->collectTemporary_terms(temporary_terms, temporary_terms_inuse, Curr_Block);
    }
}

double
AbstractExternalFunctionNode::eval(const eval_context_t &eval_context) const throw (EvalException, EvalExternalFunctionException)
{
  throw EvalExternalFunctionException();
}

int
AbstractExternalFunctionNode::maxEndoLead() const
{
  int val = 0;
  for (vector<expr_t>::const_iterator it = arguments.begin();
       it != arguments.end(); it++)
    val = max(val, (*it)->maxEndoLead());
  return val;
}

int
AbstractExternalFunctionNode::maxExoLead() const
{
  int val = 0;
  for (vector<expr_t>::const_iterator it = arguments.begin();
       it != arguments.end(); it++)
    val = max(val, (*it)->maxExoLead());
  return val;
}

int
AbstractExternalFunctionNode::maxEndoLag() const
{
  int val = 0;
  for (vector<expr_t>::const_iterator it = arguments.begin();
       it != arguments.end(); it++)
    val = max(val, (*it)->maxEndoLag());
  return val;
}

int
AbstractExternalFunctionNode::maxExoLag() const
{
  int val = 0;
  for (vector<expr_t>::const_iterator it = arguments.begin();
       it != arguments.end(); it++)
    val = max(val, (*it)->maxExoLag());
  return val;
}

int
AbstractExternalFunctionNode::maxLead() const
{
  int val = 0;
  for (vector<expr_t>::const_iterator it = arguments.begin();
       it != arguments.end(); it++)
    val = max(val, (*it)->maxLead());
  return val;
}

int
AbstractExternalFunctionNode::maxLag() const
{
  int val = 0;
  for (vector<expr_t>::const_iterator it = arguments.begin();
       it != arguments.end(); it++)
    val = max(val, (*it)->maxLag());
  return val;
}

expr_t
AbstractExternalFunctionNode::undiff() const
{
  vector<expr_t> arguments_subst;
  for (vector<expr_t>::const_iterator it = arguments.begin(); it != arguments.end(); it++)
    arguments_subst.push_back((*it)->undiff());
  return buildSimilarExternalFunctionNode(arguments_subst, datatree);
}

int
AbstractExternalFunctionNode::VarMinLag() const
{
int val = 0;
  for (vector<expr_t>::const_iterator it = arguments.begin();
       it != arguments.end(); it++)
    val = min(val, (*it)->VarMinLag());
  return val;
}

void
AbstractExternalFunctionNode::VarMaxLag(DataTree &static_datatree, set<expr_t> &static_lhs, int &max_lag) const
{
  for (vector<expr_t>::const_iterator it = arguments.begin();
       it != arguments.end(); it++)
    (*it)->VarMaxLag(static_datatree, static_lhs, max_lag);
}

int
AbstractExternalFunctionNode::PacMaxLag(vector<int> &lhs) const
{
  int val = 0;
  for (vector<expr_t>::const_iterator it = arguments.begin();
       it != arguments.end(); it++)
    val = max(val, (*it)->PacMaxLag(lhs));
  return val;
}

expr_t
AbstractExternalFunctionNode::decreaseLeadsLags(int n) const
{
  vector<expr_t> arguments_subst;
  for (vector<expr_t>::const_iterator it = arguments.begin(); it != arguments.end(); it++)
    arguments_subst.push_back((*it)->decreaseLeadsLags(n));
  return buildSimilarExternalFunctionNode(arguments_subst, datatree);
}

expr_t
AbstractExternalFunctionNode::decreaseLeadsLagsPredeterminedVariables() const
{
  vector<expr_t> arguments_subst;
  for (vector<expr_t>::const_iterator it = arguments.begin(); it != arguments.end(); it++)
    arguments_subst.push_back((*it)->decreaseLeadsLagsPredeterminedVariables());
  return buildSimilarExternalFunctionNode(arguments_subst, datatree);
}

expr_t
AbstractExternalFunctionNode::substituteEndoLeadGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const
{
  vector<expr_t> arguments_subst;
  for (vector<expr_t>::const_iterator it = arguments.begin(); it != arguments.end(); it++)
    arguments_subst.push_back((*it)->substituteEndoLeadGreaterThanTwo(subst_table, neweqs, deterministic_model));
  return buildSimilarExternalFunctionNode(arguments_subst, datatree);
}

expr_t
AbstractExternalFunctionNode::substituteEndoLagGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  vector<expr_t> arguments_subst;
  for (vector<expr_t>::const_iterator it = arguments.begin(); it != arguments.end(); it++)
    arguments_subst.push_back((*it)->substituteEndoLagGreaterThanTwo(subst_table, neweqs));
  return buildSimilarExternalFunctionNode(arguments_subst, datatree);
}

expr_t
AbstractExternalFunctionNode::substituteExoLead(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const
{
  vector<expr_t> arguments_subst;
  for (vector<expr_t>::const_iterator it = arguments.begin(); it != arguments.end(); it++)
    arguments_subst.push_back((*it)->substituteExoLead(subst_table, neweqs, deterministic_model));
  return buildSimilarExternalFunctionNode(arguments_subst, datatree);
}

expr_t
AbstractExternalFunctionNode::substituteExoLag(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  vector<expr_t> arguments_subst;
  for (vector<expr_t>::const_iterator it = arguments.begin(); it != arguments.end(); it++)
    arguments_subst.push_back((*it)->substituteExoLag(subst_table, neweqs));
  return buildSimilarExternalFunctionNode(arguments_subst, datatree);
}

expr_t
AbstractExternalFunctionNode::substituteExpectation(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool partial_information_model) const
{
  vector<expr_t> arguments_subst;
  for (vector<expr_t>::const_iterator it = arguments.begin(); it != arguments.end(); it++)
    arguments_subst.push_back((*it)->substituteExpectation(subst_table, neweqs, partial_information_model));
  return buildSimilarExternalFunctionNode(arguments_subst, datatree);
}

expr_t
AbstractExternalFunctionNode::substituteAdl() const
{
  vector<expr_t> arguments_subst;
  for (vector<expr_t>::const_iterator it = arguments.begin(); it != arguments.end(); it++)
    arguments_subst.push_back((*it)->substituteAdl());
  return buildSimilarExternalFunctionNode(arguments_subst, datatree);
}

void
AbstractExternalFunctionNode::findDiffNodes(DataTree &static_datatree, diff_table_t &diff_table) const
{
  for (vector<expr_t>::const_iterator it = arguments.begin(); it != arguments.end(); it++)
    (*it)->findDiffNodes(static_datatree, diff_table);
}

void
AbstractExternalFunctionNode::findDiffUnaryOpNodes(DataTree &static_datatree, set<UnaryOpNode *> &unary_op_nodes_in_diff) const
{
  for (vector<expr_t>::const_iterator it = arguments.begin(); it != arguments.end(); it++)
    (*it)->findDiffUnaryOpNodes(static_datatree, unary_op_nodes_in_diff);
}

expr_t
AbstractExternalFunctionNode::substituteDiff(DataTree &static_datatree, diff_table_t &diff_table, subst_table_t &subst_table,
                                             vector<BinaryOpNode *> &neweqs) const
{
  vector<expr_t> arguments_subst;
  for (vector<expr_t>::const_iterator it = arguments.begin(); it != arguments.end(); it++)
    arguments_subst.push_back((*it)->substituteDiff(static_datatree, diff_table, subst_table, neweqs));
  return buildSimilarExternalFunctionNode(arguments_subst, datatree);
}

expr_t
AbstractExternalFunctionNode::substituteDiffUnaryOpNodes(DataTree &static_datatree, set<UnaryOpNode *> &nodes, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  vector<expr_t> arguments_subst;
  for (vector<expr_t>::const_iterator it = arguments.begin(); it != arguments.end(); it++)
    arguments_subst.push_back((*it)->substituteDiffUnaryOpNodes(static_datatree, nodes, subst_table, neweqs));
  return buildSimilarExternalFunctionNode(arguments_subst, datatree);
}

bool
AbstractExternalFunctionNode::isDiffPresent() const
{
  bool result = false;
  for (vector<expr_t>::const_iterator it = arguments.begin(); it != arguments.end(); it++)
    result = result || (*it)->isDiffPresent();
  return result;
}

expr_t
AbstractExternalFunctionNode::substitutePacExpectation(map<const PacExpectationNode *, const BinaryOpNode *> &subst_table)
{
  vector<expr_t> arguments_subst;
  for (vector<expr_t>::const_iterator it = arguments.begin(); it != arguments.end(); it++)
    arguments_subst.push_back((*it)->substitutePacExpectation(subst_table));
  return buildSimilarExternalFunctionNode(arguments_subst, datatree);
}

expr_t
AbstractExternalFunctionNode::differentiateForwardVars(const vector<string> &subset, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  vector<expr_t> arguments_subst;
  for (vector<expr_t>::const_iterator it = arguments.begin(); it != arguments.end(); it++)
    arguments_subst.push_back((*it)->differentiateForwardVars(subset, subst_table, neweqs));
  return buildSimilarExternalFunctionNode(arguments_subst, datatree);
}

bool
AbstractExternalFunctionNode::alreadyWrittenAsTefTerm(int the_symb_id, deriv_node_temp_terms_t &tef_terms) const
{
  deriv_node_temp_terms_t::const_iterator it = tef_terms.find(make_pair(the_symb_id, arguments));
  if (it != tef_terms.end())
    return true;
  return false;
}

int
AbstractExternalFunctionNode::getIndxInTefTerms(int the_symb_id, deriv_node_temp_terms_t &tef_terms) const throw (UnknownFunctionNameAndArgs)
{
  deriv_node_temp_terms_t::const_iterator it = tef_terms.find(make_pair(the_symb_id, arguments));
  if (it != tef_terms.end())
    return it->second;
  throw UnknownFunctionNameAndArgs();
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
AbstractExternalFunctionNode::containsEndogenous(void) const
{
  bool result = false;
  for (vector<expr_t>::const_iterator it = arguments.begin(); it != arguments.end(); it++)
    result = result || (*it)->containsEndogenous();
  return result;
}

bool
AbstractExternalFunctionNode::containsExogenous() const
{
  for (vector<expr_t>::const_iterator it = arguments.begin();
       it != arguments.end(); it++)
    if ((*it)->containsExogenous())
      return true;
  return false;
}

expr_t
AbstractExternalFunctionNode::replaceTrendVar() const
{
  vector<expr_t> arguments_subst;
  for (vector<expr_t>::const_iterator it = arguments.begin(); it != arguments.end(); it++)
    arguments_subst.push_back((*it)->replaceTrendVar());
  return buildSimilarExternalFunctionNode(arguments_subst, datatree);
}

expr_t
AbstractExternalFunctionNode::detrend(int symb_id, bool log_trend, expr_t trend) const
{
  vector<expr_t> arguments_subst;
  for (vector<expr_t>::const_iterator it = arguments.begin(); it != arguments.end(); it++)
    arguments_subst.push_back((*it)->detrend(symb_id, log_trend, trend));
  return buildSimilarExternalFunctionNode(arguments_subst, datatree);
}

expr_t
AbstractExternalFunctionNode::removeTrendLeadLag(map<int, expr_t> trend_symbols_map) const
{
  vector<expr_t> arguments_subst;
  for (vector<expr_t>::const_iterator it = arguments.begin(); it != arguments.end(); it++)
    arguments_subst.push_back((*it)->removeTrendLeadLag(trend_symbols_map));
  return buildSimilarExternalFunctionNode(arguments_subst, datatree);
}

bool
AbstractExternalFunctionNode::isInStaticForm() const
{
  for (vector<expr_t>::const_iterator it = arguments.begin(); it != arguments.end(); ++it)
    if (!(*it)->isInStaticForm())
      return false;
  return true;
}

void
AbstractExternalFunctionNode::setVarExpectationIndex(map<string, pair<SymbolList, int> > &var_model_info)
{
  for (vector<expr_t>::const_iterator it = arguments.begin(); it != arguments.end(); it++)
    (*it)->setVarExpectationIndex(var_model_info);
}

void
AbstractExternalFunctionNode::walkPacParameters(bool &pac_encountered, pair<int, int> &lhs, set<pair<int, pair<int, int> > > &ec_params_and_vars, set<pair<int, pair<int, int> > > &ar_params_and_vars) const
{
  for (vector<expr_t>::const_iterator it = arguments.begin(); it != arguments.end(); it++)
    (*it)->walkPacParameters(pac_encountered, lhs, ec_params_and_vars, ar_params_and_vars);
}

void
AbstractExternalFunctionNode::addParamInfoToPac(pair<int, int> &lhs_arg, set<pair<int, pair<int, int> > > &ec_params_and_vars_arg, set<pair<int, pair<int, int> > > &ar_params_and_vars_arg)
{
  for (vector<expr_t>::const_iterator it = arguments.begin(); it != arguments.end(); it++)
    (*it)->addParamInfoToPac(lhs_arg, ec_params_and_vars_arg, ar_params_and_vars_arg);
}

void
AbstractExternalFunctionNode::fillPacExpectationVarInfo(string &model_name_arg, vector<int> &lhs_arg, int max_lag_arg, vector<bool> &nonstationary_arg, int growth_symb_id_arg, int equation_number_arg)
{
  for (vector<expr_t>::const_iterator it = arguments.begin(); it != arguments.end(); it++)
    (*it)->fillPacExpectationVarInfo(model_name_arg, lhs_arg, max_lag_arg, nonstationary_arg, growth_symb_id_arg, equation_number_arg);
}

bool
AbstractExternalFunctionNode::isVarModelReferenced(const string &model_info_name) const
{
  for (vector<expr_t>::const_iterator it = arguments.begin(); it != arguments.end(); ++it)
    if (!(*it)->isVarModelReferenced(model_info_name))
      return true;
  return false;
}

void
AbstractExternalFunctionNode::getEndosAndMaxLags(map<string, int> &model_endos_and_lags) const
{
  for (vector<expr_t>::const_iterator it = arguments.begin(); it != arguments.end(); ++it)
    (*it)->getEndosAndMaxLags(model_endos_and_lags);
}

pair<int, expr_t>
AbstractExternalFunctionNode::normalizeEquation(int var_endo, vector<pair<int, pair<expr_t, expr_t> > >  &List_of_Op_RHS) const
{
  vector<pair<bool, expr_t> > V_arguments;
  vector<expr_t> V_expr_t;
  bool present = false;
  for (vector<expr_t>::const_iterator it = arguments.begin();
       it != arguments.end(); it++)
    {
      V_arguments.push_back((*it)->normalizeEquation(var_endo, List_of_Op_RHS));
      present = present || V_arguments[V_arguments.size()-1].first;
      V_expr_t.push_back(V_arguments[V_arguments.size()-1].second);
    }
  if (!present)
    return (make_pair(0, datatree.AddExternalFunction(symb_id, V_expr_t)));
  else
    return (make_pair(1, (expr_t) NULL));
}

void
AbstractExternalFunctionNode::writeExternalFunctionArguments(ostream &output, ExprNodeOutputType output_type,
                                                             const temporary_terms_t &temporary_terms,
                                                             deriv_node_temp_terms_t &tef_terms) const
{
  for (vector<expr_t>::const_iterator it = arguments.begin();
       it != arguments.end(); it++)
    {
      if (it != arguments.begin())
        output << ",";

      (*it)->writeOutput(output, output_type, temporary_terms, tef_terms);
    }
}

void
AbstractExternalFunctionNode::writeJsonExternalFunctionArguments(ostream &output,
                                                                 const temporary_terms_t &temporary_terms,
                                                                 deriv_node_temp_terms_t &tef_terms,
                                                                 const bool isdynamic) const
{
  for (vector<expr_t>::const_iterator it = arguments.begin();
       it != arguments.end(); it++)
    {
      if (it != arguments.begin())
        output << ",";

      (*it)->writeJsonOutput(output, temporary_terms, tef_terms, isdynamic);
    }
}

void
AbstractExternalFunctionNode::writePrhs(ostream &output, ExprNodeOutputType output_type,
                                        const temporary_terms_t &temporary_terms,
                                        deriv_node_temp_terms_t &tef_terms, const string &ending) const
{
  output << "mxArray *prhs"<< ending << "[nrhs"<< ending << "];" << endl;
  int i = 0;
  for (vector<expr_t>::const_iterator it = arguments.begin();
       it != arguments.end(); it++)
    {
      output << "prhs" << ending << "[" << i++ << "] = mxCreateDoubleScalar("; // All external_function arguments are scalars
      (*it)->writeOutput(output, output_type, temporary_terms, tef_terms);
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
  for (vector<expr_t>::const_iterator it = arguments.begin(); it != arguments.end(); it++)
    arguments_subst.push_back((*it)->substituteStaticAuxiliaryVariable());
  return buildSimilarExternalFunctionNode(arguments_subst, datatree);
}

ExternalFunctionNode::ExternalFunctionNode(DataTree &datatree_arg,
                                           int symb_id_arg,
                                           const vector<expr_t> &arguments_arg) :
  AbstractExternalFunctionNode(datatree_arg, symb_id_arg, arguments_arg)
{
  // Add myself to the external function map
  datatree.external_function_node_map[make_pair(arguments, symb_id)] = this;
}

expr_t
ExternalFunctionNode::composeDerivatives(const vector<expr_t> &dargs)
{
  vector<expr_t> dNodes;
  for (int i = 0; i < (int) dargs.size(); i++)
    dNodes.push_back(datatree.AddTimes(dargs.at(i),
                                       datatree.AddFirstDerivExternalFunction(symb_id, arguments, i+1)));

  expr_t theDeriv = datatree.Zero;
  for (vector<expr_t>::const_iterator it = dNodes.begin(); it != dNodes.end(); it++)
    theDeriv = datatree.AddPlus(theDeriv, *it);
  return theDeriv;
}

void
ExternalFunctionNode::computeTemporaryTerms(map<expr_t, pair<int, NodeTreeReference> > &reference_count,
                                            map<NodeTreeReference, temporary_terms_t> &temp_terms_map,
                                            bool is_matlab, NodeTreeReference tr) const
{
  temp_terms_map[tr].insert(const_cast<ExternalFunctionNode *>(this));
}

void
ExternalFunctionNode::computeTemporaryTerms(map<expr_t, int> &reference_count,
                                            temporary_terms_t &temporary_terms,
                                            map<expr_t, pair<int, int> > &first_occurence,
                                            int Curr_block,
                                            vector< vector<temporary_terms_t> > &v_temporary_terms,
                                            int equation) const
{
  expr_t this2 = const_cast<ExternalFunctionNode *>(this);
  temporary_terms.insert(this2);
  first_occurence[this2] = make_pair(Curr_block, equation);
  v_temporary_terms[Curr_block][equation].insert(this2);
}

void
ExternalFunctionNode::compile(ostream &CompileCode, unsigned int &instruction_number,
                              bool lhs_rhs, const temporary_terms_t &temporary_terms,
                              const map_idx_t &map_idx, bool dynamic, bool steady_dynamic,
                              deriv_node_temp_terms_t &tef_terms) const
{
  temporary_terms_t::const_iterator it = temporary_terms.find(const_cast<ExternalFunctionNode *>(this));
  if (it != temporary_terms.end())
    {
      if (dynamic)
        {
          map_idx_t::const_iterator ii = map_idx.find(idx);
          FLDT_ fldt(ii->second);
          fldt.write(CompileCode, instruction_number);
        }
      else
        {
          map_idx_t::const_iterator ii = map_idx.find(idx);
          FLDST_ fldst(ii->second);
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
                                                    const map_idx_t &map_idx, bool dynamic, bool steady_dynamic,
                                                    deriv_node_temp_terms_t &tef_terms) const
{
  int first_deriv_symb_id = datatree.external_functions_table.getFirstDerivSymbID(symb_id);
  assert(first_deriv_symb_id != eExtFunSetButNoNameProvided);

  for (vector<expr_t>::const_iterator it = arguments.begin();
       it != arguments.end(); it++)
    (*it)->compileExternalFunctionOutput(CompileCode, instruction_number, lhs_rhs, temporary_terms,
                                         map_idx, dynamic, steady_dynamic, tef_terms);

  if (!alreadyWrittenAsTefTerm(symb_id, tef_terms))
    {
      tef_terms[make_pair(symb_id, arguments)] = (int) tef_terms.size();
      int indx = getIndxInTefTerms(symb_id, tef_terms);
      int second_deriv_symb_id = datatree.external_functions_table.getSecondDerivSymbID(symb_id);
      assert(second_deriv_symb_id != eExtFunSetButNoNameProvided);

      unsigned int nb_output_arguments = 0;
      if (symb_id == first_deriv_symb_id
          && symb_id == second_deriv_symb_id)
        nb_output_arguments = 3;
      else if (symb_id == first_deriv_symb_id)
        nb_output_arguments = 2;
      else
        nb_output_arguments = 1;
      unsigned int nb_input_arguments = compileExternalFunctionArguments(CompileCode, instruction_number, lhs_rhs, temporary_terms,
                                                                         map_idx, dynamic, steady_dynamic, tef_terms);

      FCALL_ fcall(nb_output_arguments, nb_input_arguments, datatree.symbol_table.getName(symb_id), indx);
      switch (nb_output_arguments)
        {
        case 1:
          fcall.set_function_type(ExternalFunctionWithoutDerivative);
          break;
        case 2:
          fcall.set_function_type(ExternalFunctionWithFirstDerivative);
          break;
        case 3:
          fcall.set_function_type(ExternalFunctionWithFirstandSecondDerivative);
          break;
        }
      fcall.write(CompileCode, instruction_number);
      FSTPTEF_ fstptef(indx);
      fstptef.write(CompileCode, instruction_number);
    }
}

void
ExternalFunctionNode::writeJsonOutput(ostream &output,
                                      const temporary_terms_t &temporary_terms,
                                      deriv_node_temp_terms_t &tef_terms,
                                      const bool isdynamic) const
{
  temporary_terms_t::const_iterator it = temporary_terms.find(const_cast<ExternalFunctionNode *>(this));
  if (it != temporary_terms.end())
    {
      output << "T" << idx;
      return;
    }

  output << datatree.symbol_table.getName(symb_id) << "(";
  writeJsonExternalFunctionArguments(output, temporary_terms, tef_terms, isdynamic);
  output << ")";
}

void
ExternalFunctionNode::writeOutput(ostream &output, ExprNodeOutputType output_type,
                                  const temporary_terms_t &temporary_terms,
                                  deriv_node_temp_terms_t &tef_terms) const
{
  if (output_type == oMatlabOutsideModel || output_type == oSteadyStateFile
      || output_type == oCSteadyStateFile || output_type == oJuliaSteadyStateFile
      || IS_LATEX(output_type))
    {
      string name = IS_LATEX(output_type) ? datatree.symbol_table.getTeXName(symb_id)
        : datatree.symbol_table.getName(symb_id);
      output << name << "(";
      writeExternalFunctionArguments(output, output_type, temporary_terms, tef_terms);
      output << ")";
      return;
    }

  temporary_terms_t::const_iterator it = temporary_terms.find(const_cast<ExternalFunctionNode *>(this));
  if (it != temporary_terms.end())
    {
      if (output_type == oMatlabDynamicModelSparse)
        output << "T" << idx << "(it_)";
      else
        output << "T" << idx;
      return;
    }

  if (IS_C(output_type))
    output << "*";
  output << "TEF_" << getIndxInTefTerms(symb_id, tef_terms);
}

void
ExternalFunctionNode::writeExternalFunctionOutput(ostream &output, ExprNodeOutputType output_type,
                                                  const temporary_terms_t &temporary_terms,
                                                  deriv_node_temp_terms_t &tef_terms) const
{
  int first_deriv_symb_id = datatree.external_functions_table.getFirstDerivSymbID(symb_id);
  assert(first_deriv_symb_id != eExtFunSetButNoNameProvided);

  for (vector<expr_t>::const_iterator it = arguments.begin();
       it != arguments.end(); it++)
    (*it)->writeExternalFunctionOutput(output, output_type, temporary_terms, tef_terms);

  if (!alreadyWrittenAsTefTerm(symb_id, tef_terms))
    {
      tef_terms[make_pair(symb_id, arguments)] = (int) tef_terms.size();
      int indx = getIndxInTefTerms(symb_id, tef_terms);
      int second_deriv_symb_id = datatree.external_functions_table.getSecondDerivSymbID(symb_id);
      assert(second_deriv_symb_id != eExtFunSetButNoNameProvided);

      if (IS_C(output_type))
        {
          stringstream ending;
          ending << "_tef_" << getIndxInTefTerms(symb_id, tef_terms);
          if (symb_id == first_deriv_symb_id
              && symb_id == second_deriv_symb_id)
            output << "int nlhs" << ending.str() << " = 3;" << endl
                   << "double *TEF_" << indx << ", "
                   << "*TEFD_" << indx << ", "
                   << "*TEFDD_" << indx << ";" << endl;
          else if (symb_id == first_deriv_symb_id)
            output << "int nlhs" << ending.str() << " = 2;" << endl
                   << "double *TEF_" << indx << ", "
                   << "*TEFD_" << indx << "; " << endl;
          else
            output << "int nlhs" << ending.str() << " = 1;" << endl
                   << "double *TEF_" << indx << ";" << endl;

          output << "mxArray *plhs" << ending.str()<< "[nlhs"<< ending.str() << "];" << endl;
          output << "int nrhs" << ending.str()<< " = " << arguments.size() << ";" << endl;
          writePrhs(output, output_type, temporary_terms, tef_terms, ending.str());

          output << "mexCallMATLAB("
                 << "nlhs" << ending.str() << ", "
                 << "plhs" << ending.str() << ", "
                 << "nrhs" << ending.str() << ", "
                 << "prhs" << ending.str() << ", \""
                 << datatree.symbol_table.getName(symb_id) << "\");" << endl;

          if (symb_id == first_deriv_symb_id
              && symb_id == second_deriv_symb_id)
            output << "TEF_" << indx << " = mxGetPr(plhs" << ending.str() << "[0]);" << endl
                   << "TEFD_" << indx << " = mxGetPr(plhs" << ending.str() << "[1]);" << endl
                   << "TEFDD_" << indx << " = mxGetPr(plhs" << ending.str() << "[2]);" << endl
                   << "int TEFDD_" << indx << "_nrows = (int)mxGetM(plhs" << ending.str()<< "[2]);" << endl;
          else if (symb_id == first_deriv_symb_id)
            output << "TEF_" << indx << " = mxGetPr(plhs" << ending.str() << "[0]);" << endl
                   << "TEFD_" << indx << " = mxGetPr(plhs" << ending.str() << "[1]);" << endl;
          else
            output << "TEF_" << indx << " = mxGetPr(plhs" << ending.str() << "[0]);" << endl;
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
          writeExternalFunctionArguments(output, output_type, temporary_terms, tef_terms);
          output << ");" << endl;
        }
    }
}

void
ExternalFunctionNode::writeJsonExternalFunctionOutput(vector<string> &efout,
                                                      const temporary_terms_t &temporary_terms,
                                                      deriv_node_temp_terms_t &tef_terms,
                                                      const bool isdynamic) const
{
  int first_deriv_symb_id = datatree.external_functions_table.getFirstDerivSymbID(symb_id);
  assert(first_deriv_symb_id != eExtFunSetButNoNameProvided);

  for (vector<expr_t>::const_iterator it = arguments.begin();
       it != arguments.end(); it++)
    (*it)->writeJsonExternalFunctionOutput(efout, temporary_terms, tef_terms, isdynamic);

  if (!alreadyWrittenAsTefTerm(symb_id, tef_terms))
    {
      tef_terms[make_pair(symb_id, arguments)] = (int) tef_terms.size();
      int indx = getIndxInTefTerms(symb_id, tef_terms);
      int second_deriv_symb_id = datatree.external_functions_table.getSecondDerivSymbID(symb_id);
      assert(second_deriv_symb_id != eExtFunSetButNoNameProvided);

      stringstream ef;
      ef << "{\"external_function\": {"
         << "\"external_function_term\": \"TEF_" << indx << "\"";

      if (symb_id == first_deriv_symb_id)
        ef << ", \"external_function_term_d\": \"TEFD_" << indx << "\"";

      if (symb_id == second_deriv_symb_id)
        ef << ", \"external_function_term_dd\": \"TEFDD_" << indx << "\"";

      ef << ", \"value\": \"" << datatree.symbol_table.getName(symb_id) << "(";
      writeJsonExternalFunctionArguments(ef, temporary_terms, tef_terms, isdynamic);
      ef << ")\"}}";
      efout.push_back(ef.str());
    }
}

expr_t
ExternalFunctionNode::toStatic(DataTree &static_datatree) const
{
  vector<expr_t> static_arguments;
  for (vector<expr_t>::const_iterator it = arguments.begin();
       it != arguments.end(); it++)
    static_arguments.push_back((*it)->toStatic(static_datatree));
  return static_datatree.AddExternalFunction(symb_id, static_arguments);
}

void
ExternalFunctionNode::computeXrefs(EquationInfo &ei) const
{
  vector<expr_t> dynamic_arguments;
  for (vector<expr_t>::const_iterator it = arguments.begin();
       it != arguments.end(); it++)
    (*it)->computeXrefs(ei);
}

expr_t
ExternalFunctionNode::cloneDynamic(DataTree &dynamic_datatree) const
{
  vector<expr_t> dynamic_arguments;
  for (vector<expr_t>::const_iterator it = arguments.begin();
       it != arguments.end(); it++)
    dynamic_arguments.push_back((*it)->cloneDynamic(dynamic_datatree));
  return dynamic_datatree.AddExternalFunction(symb_id, dynamic_arguments);
}

expr_t
ExternalFunctionNode::buildSimilarExternalFunctionNode(vector<expr_t> &alt_args, DataTree &alt_datatree) const
{
  return alt_datatree.AddExternalFunction(symb_id, alt_args);
}

FirstDerivExternalFunctionNode::FirstDerivExternalFunctionNode(DataTree &datatree_arg,
                                                               int top_level_symb_id_arg,
                                                               const vector<expr_t> &arguments_arg,
                                                               int inputIndex_arg) :
  AbstractExternalFunctionNode(datatree_arg, top_level_symb_id_arg, arguments_arg),
  inputIndex(inputIndex_arg)
{
  // Add myself to the first derivative external function map
  datatree.first_deriv_external_function_node_map[make_pair(make_pair(arguments, inputIndex), symb_id)] = this;
}

void
FirstDerivExternalFunctionNode::computeTemporaryTerms(map<expr_t, pair<int, NodeTreeReference> > &reference_count,
                                                      map<NodeTreeReference, temporary_terms_t> &temp_terms_map,
                                                      bool is_matlab, NodeTreeReference tr) const
{
  temp_terms_map[tr].insert(const_cast<FirstDerivExternalFunctionNode *>(this));
}

void
FirstDerivExternalFunctionNode::computeTemporaryTerms(map<expr_t, int> &reference_count,
                                                      temporary_terms_t &temporary_terms,
                                                      map<expr_t, pair<int, int> > &first_occurence,
                                                      int Curr_block,
                                                      vector< vector<temporary_terms_t> > &v_temporary_terms,
                                                      int equation) const
{
  expr_t this2 = const_cast<FirstDerivExternalFunctionNode *>(this);
  temporary_terms.insert(this2);
  first_occurence[this2] = make_pair(Curr_block, equation);
  v_temporary_terms[Curr_block][equation].insert(this2);
}

expr_t
FirstDerivExternalFunctionNode::composeDerivatives(const vector<expr_t> &dargs)
{
  vector<expr_t> dNodes;
  for (int i = 0; i < (int) dargs.size(); i++)
    dNodes.push_back(datatree.AddTimes(dargs.at(i),
                                       datatree.AddSecondDerivExternalFunction(symb_id, arguments, inputIndex, i+1)));
  expr_t theDeriv = datatree.Zero;
  for (vector<expr_t>::const_iterator it = dNodes.begin(); it != dNodes.end(); it++)
    theDeriv = datatree.AddPlus(theDeriv, *it);
  return theDeriv;
}

void
FirstDerivExternalFunctionNode::writeJsonOutput(ostream &output,
                                                const temporary_terms_t &temporary_terms,
                                                deriv_node_temp_terms_t &tef_terms,
                                                const bool isdynamic) const
{
  // If current node is a temporary term
  temporary_terms_t::const_iterator it = temporary_terms.find(const_cast<FirstDerivExternalFunctionNode *>(this));
  if (it != temporary_terms.end())
    {
      output << "T" << idx;
      return;
    }

  const int first_deriv_symb_id = datatree.external_functions_table.getFirstDerivSymbID(symb_id);
  assert(first_deriv_symb_id != eExtFunSetButNoNameProvided);

  const int tmpIndx = inputIndex - 1;

  if (first_deriv_symb_id == symb_id)
    output << "TEFD_" << getIndxInTefTerms(symb_id, tef_terms)
           << "[" << tmpIndx << "]";
  else if (first_deriv_symb_id == eExtFunNotSet)
    output << "TEFD_fdd_" << getIndxInTefTerms(symb_id, tef_terms) << "_" << inputIndex;
  else
    output << "TEFD_def_" << getIndxInTefTerms(first_deriv_symb_id, tef_terms)
           << "[" << tmpIndx << "]";
}

void
FirstDerivExternalFunctionNode::writeOutput(ostream &output, ExprNodeOutputType output_type,
                                            const temporary_terms_t &temporary_terms,
                                            deriv_node_temp_terms_t &tef_terms) const
{
  assert(output_type != oMatlabOutsideModel);

  if (IS_LATEX(output_type))
    {
      output << "\\frac{\\partial " << datatree.symbol_table.getTeXName(symb_id)
             << "}{\\partial " << inputIndex << "}(";
      writeExternalFunctionArguments(output, output_type, temporary_terms, tef_terms);
      output << ")";
      return;
    }

  // If current node is a temporary term
  temporary_terms_t::const_iterator it = temporary_terms.find(const_cast<FirstDerivExternalFunctionNode *>(this));
  if (it != temporary_terms.end())
    {
      if (output_type == oMatlabDynamicModelSparse)
        output << "T" << idx << "(it_)";
      else
        output << "T" << idx;
      return;
    }

  const int first_deriv_symb_id = datatree.external_functions_table.getFirstDerivSymbID(symb_id);
  assert(first_deriv_symb_id != eExtFunSetButNoNameProvided);

  const int tmpIndx = inputIndex - 1 + ARRAY_SUBSCRIPT_OFFSET(output_type);

  if (first_deriv_symb_id == symb_id)
    output << "TEFD_" << getIndxInTefTerms(symb_id, tef_terms)
           << LEFT_ARRAY_SUBSCRIPT(output_type) << tmpIndx << RIGHT_ARRAY_SUBSCRIPT(output_type);
  else if (first_deriv_symb_id == eExtFunNotSet)
    {
      if (IS_C(output_type))
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
                                        const map_idx_t &map_idx, bool dynamic, bool steady_dynamic,
                                        deriv_node_temp_terms_t &tef_terms) const
{
  temporary_terms_t::const_iterator it = temporary_terms.find(const_cast<FirstDerivExternalFunctionNode *>(this));
  if (it != temporary_terms.end())
    {
      if (dynamic)
        {
          map_idx_t::const_iterator ii = map_idx.find(idx);
          FLDT_ fldt(ii->second);
          fldt.write(CompileCode, instruction_number);
        }
      else
        {
          map_idx_t::const_iterator ii = map_idx.find(idx);
          FLDST_ fldst(ii->second);
          fldst.write(CompileCode, instruction_number);
        }
      return;
    }
  int first_deriv_symb_id = datatree.external_functions_table.getFirstDerivSymbID(symb_id);
  assert(first_deriv_symb_id != eExtFunSetButNoNameProvided);

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
                                                            deriv_node_temp_terms_t &tef_terms) const
{
  assert(output_type != oMatlabOutsideModel);
  int first_deriv_symb_id = datatree.external_functions_table.getFirstDerivSymbID(symb_id);
  assert(first_deriv_symb_id != eExtFunSetButNoNameProvided);

  /* For a node with derivs provided by the user function, call the method
     on the non-derived node */
  if (first_deriv_symb_id == symb_id)
    {
      expr_t parent = datatree.AddExternalFunction(symb_id, arguments);
      parent->writeExternalFunctionOutput(output, output_type, temporary_terms,
                                          tef_terms);
      return;
    }

  if (alreadyWrittenAsTefTerm(first_deriv_symb_id, tef_terms))
    return;

  if (IS_C(output_type))
    if (first_deriv_symb_id == eExtFunNotSet)
      {
        stringstream ending;
        ending << "_tefd_fdd_" << getIndxInTefTerms(symb_id, tef_terms) << "_" << inputIndex;
        output << "int nlhs" << ending.str() << " = 1;" << endl
               << "double *TEFD_fdd_" <<  getIndxInTefTerms(symb_id, tef_terms) << "_" << inputIndex << ";" << endl
               << "mxArray *plhs" << ending.str() << "[nlhs"<< ending.str() << "];" << endl
               << "int nrhs" << ending.str() << " = 3;" << endl
               << "mxArray *prhs" << ending.str() << "[nrhs"<< ending.str() << "];" << endl
               << "mwSize dims" << ending.str() << "[2];" << endl;

        output << "dims" << ending.str() << "[0] = 1;" << endl
               << "dims" << ending.str() << "[1] = " << arguments.size() << ";" << endl;

        output << "prhs" << ending.str() << "[0] = mxCreateString(\"" << datatree.symbol_table.getName(symb_id) << "\");" << endl
               << "prhs" << ending.str() << "[1] = mxCreateDoubleScalar(" << inputIndex << ");"<< endl
               << "prhs" << ending.str() << "[2] = mxCreateCellArray(2, dims" << ending.str() << ");"<< endl;

        int i = 0;
        for (vector<expr_t>::const_iterator it = arguments.begin();
             it != arguments.end(); it++)
          {
            output << "mxSetCell(prhs" << ending.str() << "[2], "
                   << i++ << ", "
                   << "mxCreateDoubleScalar("; // All external_function arguments are scalars
            (*it)->writeOutput(output, output_type, temporary_terms, tef_terms);
            output << "));" << endl;
          }

        output << "mexCallMATLAB("
               << "nlhs" << ending.str() << ", "
               << "plhs" << ending.str() << ", "
               << "nrhs" << ending.str() << ", "
               << "prhs" << ending.str() << ", \""
               << "jacob_element\");" << endl;

        output << "TEFD_fdd_" <<  getIndxInTefTerms(symb_id, tef_terms) << "_" << inputIndex
               << " = mxGetPr(plhs" << ending.str() << "[0]);" << endl;
      }
    else
      {
        tef_terms[make_pair(first_deriv_symb_id, arguments)] = (int) tef_terms.size();
        int indx = getIndxInTefTerms(first_deriv_symb_id, tef_terms);
        stringstream ending;
        ending << "_tefd_def_" << indx;
        output << "int nlhs" << ending.str() << " = 1;" << endl
               << "double *TEFD_def_" << indx << ";" << endl
               << "mxArray *plhs" << ending.str() << "[nlhs"<< ending.str() << "];" << endl
               << "int nrhs" << ending.str() << " = " << arguments.size() << ";" << endl;
        writePrhs(output, output_type, temporary_terms, tef_terms, ending.str());

        output << "mexCallMATLAB("
               << "nlhs" << ending.str() << ", "
               << "plhs" << ending.str() << ", "
               << "nrhs" << ending.str() << ", "
               << "prhs" << ending.str() << ", \""
               << datatree.symbol_table.getName(first_deriv_symb_id) << "\");" << endl;

        output << "TEFD_def_" << indx << " = mxGetPr(plhs" << ending.str() << "[0]);" << endl;
      }
  else
    {
      if (first_deriv_symb_id == eExtFunNotSet)
        output << "TEFD_fdd_" << getIndxInTefTerms(symb_id, tef_terms) << "_" << inputIndex << " = jacob_element('"
               << datatree.symbol_table.getName(symb_id) << "'," << inputIndex << ",{";
      else
        {
          tef_terms[make_pair(first_deriv_symb_id, arguments)] = (int) tef_terms.size();
          output << "TEFD_def_" << getIndxInTefTerms(first_deriv_symb_id, tef_terms)
                 << " = " << datatree.symbol_table.getName(first_deriv_symb_id) << "(";
        }

      writeExternalFunctionArguments(output, output_type, temporary_terms, tef_terms);

      if (first_deriv_symb_id == eExtFunNotSet)
        output << "}";
      output << ");" << endl;
    }
}

void
FirstDerivExternalFunctionNode::writeJsonExternalFunctionOutput(vector<string> &efout,
                                                                const temporary_terms_t &temporary_terms,
                                                                deriv_node_temp_terms_t &tef_terms,
                                                                const bool isdynamic) const
{
  int first_deriv_symb_id = datatree.external_functions_table.getFirstDerivSymbID(symb_id);
  assert(first_deriv_symb_id != eExtFunSetButNoNameProvided);

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
  if (first_deriv_symb_id == eExtFunNotSet)
    ef << "{\"first_deriv_external_function\": {"
       << "\"external_function_term\": \"TEFD_fdd_" << getIndxInTefTerms(symb_id, tef_terms) << "_" << inputIndex << "\""
       << ", \"analytic_derivative\": false"
       << ", \"wrt\": " << inputIndex
       << ", \"value\": \"" << datatree.symbol_table.getName(symb_id) << "(";
  else
    {
      tef_terms[make_pair(first_deriv_symb_id, arguments)] = (int) tef_terms.size();
      ef << "{\"first_deriv_external_function\": {"
         << "\"external_function_term\": \"TEFD_def_" << getIndxInTefTerms(first_deriv_symb_id, tef_terms) << "\""
         << ", \"analytic_derivative\": true"
         << ", \"value\": \"" << datatree.symbol_table.getName(first_deriv_symb_id) << "(";
    }

  writeJsonExternalFunctionArguments(ef, temporary_terms, tef_terms, isdynamic);
  ef << ")\"}}";
  efout.push_back(ef.str());
}

void
FirstDerivExternalFunctionNode::compileExternalFunctionOutput(ostream &CompileCode, unsigned int &instruction_number,
                                                              bool lhs_rhs, const temporary_terms_t &temporary_terms,
                                                              const map_idx_t &map_idx, bool dynamic, bool steady_dynamic,
                                                              deriv_node_temp_terms_t &tef_terms) const
{
  int first_deriv_symb_id = datatree.external_functions_table.getFirstDerivSymbID(symb_id);
  assert(first_deriv_symb_id != eExtFunSetButNoNameProvided);

  if (first_deriv_symb_id == symb_id || alreadyWrittenAsTefTerm(first_deriv_symb_id, tef_terms))
    return;

  unsigned int nb_add_input_arguments = compileExternalFunctionArguments(CompileCode, instruction_number, lhs_rhs, temporary_terms,
                                                                         map_idx, dynamic, steady_dynamic, tef_terms);
  if (first_deriv_symb_id == eExtFunNotSet)
    {
      unsigned int nb_input_arguments = 0;
      unsigned int nb_output_arguments = 1;
      unsigned int indx = getIndxInTefTerms(symb_id, tef_terms);
      FCALL_ fcall(nb_output_arguments, nb_input_arguments, "jacob_element", indx);
      fcall.set_arg_func_name(datatree.symbol_table.getName(symb_id));
      fcall.set_row(inputIndex);
      fcall.set_nb_add_input_arguments(nb_add_input_arguments);
      fcall.set_function_type(ExternalFunctionNumericalFirstDerivative);
      fcall.write(CompileCode, instruction_number);
      FSTPTEFD_ fstptefd(indx, inputIndex);
      fstptefd.write(CompileCode, instruction_number);
    }
  else
    {
      tef_terms[make_pair(first_deriv_symb_id, arguments)] = (int) tef_terms.size();
      int indx = getIndxInTefTerms(symb_id, tef_terms);
      int second_deriv_symb_id = datatree.external_functions_table.getSecondDerivSymbID(symb_id);
      assert(second_deriv_symb_id != eExtFunSetButNoNameProvided);

      unsigned int nb_output_arguments = 1;

      FCALL_ fcall(nb_output_arguments, nb_add_input_arguments, datatree.symbol_table.getName(first_deriv_symb_id), indx);
      fcall.set_function_type(ExternalFunctionFirstDerivative);
      fcall.write(CompileCode, instruction_number);
      FSTPTEFD_ fstptefd(indx, inputIndex);
      fstptefd.write(CompileCode, instruction_number);
    }
}

expr_t
FirstDerivExternalFunctionNode::cloneDynamic(DataTree &dynamic_datatree) const
{
  vector<expr_t> dynamic_arguments;
  for (vector<expr_t>::const_iterator it = arguments.begin();
       it != arguments.end(); it++)
    dynamic_arguments.push_back((*it)->cloneDynamic(dynamic_datatree));
  return dynamic_datatree.AddFirstDerivExternalFunction(symb_id, dynamic_arguments,
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
  for (vector<expr_t>::const_iterator it = arguments.begin();
       it != arguments.end(); it++)
    static_arguments.push_back((*it)->toStatic(static_datatree));
  return static_datatree.AddFirstDerivExternalFunction(symb_id, static_arguments,
                                                       inputIndex);
}

void
FirstDerivExternalFunctionNode::computeXrefs(EquationInfo &ei) const
{
  vector<expr_t> dynamic_arguments;
  for (vector<expr_t>::const_iterator it = arguments.begin();
       it != arguments.end(); it++)
    (*it)->computeXrefs(ei);
}

SecondDerivExternalFunctionNode::SecondDerivExternalFunctionNode(DataTree &datatree_arg,
                                                                 int top_level_symb_id_arg,
                                                                 const vector<expr_t> &arguments_arg,
                                                                 int inputIndex1_arg,
                                                                 int inputIndex2_arg) :
  AbstractExternalFunctionNode(datatree_arg, top_level_symb_id_arg, arguments_arg),
  inputIndex1(inputIndex1_arg),
  inputIndex2(inputIndex2_arg)
{
  // Add myself to the second derivative external function map
  datatree.second_deriv_external_function_node_map[make_pair(make_pair(arguments, make_pair(inputIndex1, inputIndex2)), symb_id)] = this;
}

void
SecondDerivExternalFunctionNode::computeTemporaryTerms(map<expr_t, pair<int, NodeTreeReference> > &reference_count,
                                                       map<NodeTreeReference, temporary_terms_t> &temp_terms_map,
                                                       bool is_matlab, NodeTreeReference tr) const
{
  temp_terms_map[tr].insert(const_cast<SecondDerivExternalFunctionNode *>(this));
}

void
SecondDerivExternalFunctionNode::computeTemporaryTerms(map<expr_t, int> &reference_count,
                                                       temporary_terms_t &temporary_terms,
                                                       map<expr_t, pair<int, int> > &first_occurence,
                                                       int Curr_block,
                                                       vector< vector<temporary_terms_t> > &v_temporary_terms,
                                                       int equation) const
{
  expr_t this2 = const_cast<SecondDerivExternalFunctionNode *>(this);
  temporary_terms.insert(this2);
  first_occurence[this2] = make_pair(Curr_block, equation);
  v_temporary_terms[Curr_block][equation].insert(this2);
}

expr_t
SecondDerivExternalFunctionNode::composeDerivatives(const vector<expr_t> &dargs)

{
  cerr << "ERROR: third order derivatives of external functions are not implemented" << endl;
  exit(EXIT_FAILURE);
}

void
SecondDerivExternalFunctionNode::writeJsonOutput(ostream &output,
                                                 const temporary_terms_t &temporary_terms,
                                                 deriv_node_temp_terms_t &tef_terms,
                                                 const bool isdynamic) const
{
  // If current node is a temporary term
  temporary_terms_t::const_iterator it = temporary_terms.find(const_cast<SecondDerivExternalFunctionNode *>(this));
  if (it != temporary_terms.end())
    {
      output << "T" << idx;
      return;
    }

  const int second_deriv_symb_id = datatree.external_functions_table.getSecondDerivSymbID(symb_id);
  assert(second_deriv_symb_id != eExtFunSetButNoNameProvided);

  const int tmpIndex1 = inputIndex1 - 1;
  const int tmpIndex2 = inputIndex2 - 1;

  if (second_deriv_symb_id == symb_id)
    output << "TEFDD_" << getIndxInTefTerms(symb_id, tef_terms)
           << "[" << tmpIndex1 << "," << tmpIndex2 << "]";
  else if (second_deriv_symb_id == eExtFunNotSet)
    output << "TEFDD_fdd_" << getIndxInTefTerms(symb_id, tef_terms) << "_" << inputIndex1 << "_" << inputIndex2;
  else
    output << "TEFDD_def_" << getIndxInTefTerms(second_deriv_symb_id, tef_terms)
           << "[" << tmpIndex1 << "," << tmpIndex2 << "]";
}

void
SecondDerivExternalFunctionNode::writeOutput(ostream &output, ExprNodeOutputType output_type,
                                             const temporary_terms_t &temporary_terms,
                                             deriv_node_temp_terms_t &tef_terms) const
{
  assert(output_type != oMatlabOutsideModel);

  if (IS_LATEX(output_type))
    {
      output << "\\frac{\\partial^2 " << datatree.symbol_table.getTeXName(symb_id)
             << "}{\\partial " << inputIndex1 << "\\partial " << inputIndex2 << "}(";
      writeExternalFunctionArguments(output, output_type, temporary_terms, tef_terms);
      output << ")";
      return;
    }

  // If current node is a temporary term
  temporary_terms_t::const_iterator it = temporary_terms.find(const_cast<SecondDerivExternalFunctionNode *>(this));
  if (it != temporary_terms.end())
    {
      if (output_type == oMatlabDynamicModelSparse)
        output << "T" << idx << "(it_)";
      else
        output << "T" << idx;
      return;
    }

  const int second_deriv_symb_id = datatree.external_functions_table.getSecondDerivSymbID(symb_id);
  assert(second_deriv_symb_id != eExtFunSetButNoNameProvided);

  const int tmpIndex1 = inputIndex1 - 1 + ARRAY_SUBSCRIPT_OFFSET(output_type);
  const int tmpIndex2 = inputIndex2 - 1 + ARRAY_SUBSCRIPT_OFFSET(output_type);

  int indx = getIndxInTefTerms(symb_id, tef_terms);
  if (second_deriv_symb_id == symb_id)
    if (IS_C(output_type))
      output << "TEFDD_" << indx
             << LEFT_ARRAY_SUBSCRIPT(output_type) << tmpIndex1 << " * TEFDD_" << indx << "_nrows + "
             << tmpIndex2 << RIGHT_ARRAY_SUBSCRIPT(output_type);
    else
      output << "TEFDD_" << getIndxInTefTerms(symb_id, tef_terms)
             << LEFT_ARRAY_SUBSCRIPT(output_type) << tmpIndex1 << "," << tmpIndex2 << RIGHT_ARRAY_SUBSCRIPT(output_type);
  else if (second_deriv_symb_id == eExtFunNotSet)
    {
      if (IS_C(output_type))
        output << "*";
      output << "TEFDD_fdd_" << getIndxInTefTerms(symb_id, tef_terms) << "_" << inputIndex1 << "_" << inputIndex2;
    }
  else
    if (IS_C(output_type))
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
                                                             deriv_node_temp_terms_t &tef_terms) const
{
  assert(output_type != oMatlabOutsideModel);
  int second_deriv_symb_id = datatree.external_functions_table.getSecondDerivSymbID(symb_id);
  assert(second_deriv_symb_id != eExtFunSetButNoNameProvided);

  /* For a node with derivs provided by the user function, call the method
     on the non-derived node */
  if (second_deriv_symb_id == symb_id)
    {
      expr_t parent = datatree.AddExternalFunction(symb_id, arguments);
      parent->writeExternalFunctionOutput(output, output_type, temporary_terms,
                                          tef_terms);
      return;
    }

  if (alreadyWrittenAsTefTerm(second_deriv_symb_id, tef_terms))
    return;

  if (IS_C(output_type))
    if (second_deriv_symb_id == eExtFunNotSet)
      {
        stringstream ending;
        ending << "_tefdd_fdd_" << getIndxInTefTerms(symb_id, tef_terms) << "_" << inputIndex1 << "_" << inputIndex2;
        output << "int nlhs" << ending.str() << " = 1;" << endl
               << "double *TEFDD_fdd_" <<  getIndxInTefTerms(symb_id, tef_terms) << "_" << inputIndex1 << "_" << inputIndex2 << ";" << endl
               << "mxArray *plhs" << ending.str() << "[nlhs"<< ending.str() << "];" << endl
               << "int nrhs" << ending.str() << " = 4;" << endl
               << "mxArray *prhs" << ending.str() << "[nrhs"<< ending.str() << "];" << endl
               << "mwSize dims" << ending.str() << "[2];" << endl;

        output << "dims" << ending.str() << "[0] = 1;" << endl
               << "dims" << ending.str() << "[1] = " << arguments.size() << ";" << endl;

        output << "prhs" << ending.str() << "[0] = mxCreateString(\"" << datatree.symbol_table.getName(symb_id) << "\");" << endl
               << "prhs" << ending.str() << "[1] = mxCreateDoubleScalar(" << inputIndex1 << ");"<< endl
               << "prhs" << ending.str() << "[2] = mxCreateDoubleScalar(" << inputIndex2 << ");"<< endl
               << "prhs" << ending.str() << "[3] = mxCreateCellArray(2, dims" << ending.str() << ");"<< endl;

        int i = 0;
        for (vector<expr_t>::const_iterator it = arguments.begin();
             it != arguments.end(); it++)
          {
            output << "mxSetCell(prhs" << ending.str() << "[3], "
                   << i++ << ", "
                   << "mxCreateDoubleScalar("; // All external_function arguments are scalars
            (*it)->writeOutput(output, output_type, temporary_terms, tef_terms);
            output << "));" << endl;
          }

        output << "mexCallMATLAB("
               << "nlhs" << ending.str() << ", "
               << "plhs" << ending.str() << ", "
               << "nrhs" << ending.str() << ", "
               << "prhs" << ending.str() << ", \""
               << "hess_element\");" << endl;

        output << "TEFDD_fdd_" <<  getIndxInTefTerms(symb_id, tef_terms) << "_" << inputIndex1 << "_" << inputIndex2
               << " = mxGetPr(plhs" << ending.str() << "[0]);" << endl;
      }
    else
      {
        tef_terms[make_pair(second_deriv_symb_id, arguments)] = (int) tef_terms.size();
        int indx = getIndxInTefTerms(second_deriv_symb_id, tef_terms);
        stringstream ending;
        ending << "_tefdd_def_" << indx;

        output << "int nlhs" << ending.str() << " = 1;" << endl
               << "double *TEFDD_def_" << indx << ";" << endl
               << "mxArray *plhs" << ending.str() << "[nlhs"<< ending.str() << "];" << endl
               << "int nrhs" << ending.str() << " = " << arguments.size() << ";" << endl;
        writePrhs(output, output_type, temporary_terms, tef_terms, ending.str());

        output << "mexCallMATLAB("
               << "nlhs" << ending.str() << ", "
               << "plhs" << ending.str() << ", "
               << "nrhs" << ending.str() << ", "
               << "prhs" << ending.str() << ", \""
               << datatree.symbol_table.getName(second_deriv_symb_id) << "\");" << endl;

        output << "TEFDD_def_" << indx << " = mxGetPr(plhs" << ending.str() << "[0]);" << endl;
      }
  else
    {
      if (second_deriv_symb_id == eExtFunNotSet)
        output << "TEFDD_fdd_" << getIndxInTefTerms(symb_id, tef_terms) << "_" << inputIndex1 << "_" << inputIndex2
               << " = hess_element('" << datatree.symbol_table.getName(symb_id) << "',"
               << inputIndex1 << "," << inputIndex2 << ",{";
      else
        {
          tef_terms[make_pair(second_deriv_symb_id, arguments)] = (int) tef_terms.size();
          output << "TEFDD_def_" << getIndxInTefTerms(second_deriv_symb_id, tef_terms)
                 << " = " << datatree.symbol_table.getName(second_deriv_symb_id) << "(";
        }

      writeExternalFunctionArguments(output, output_type, temporary_terms, tef_terms);

      if (second_deriv_symb_id == eExtFunNotSet)
        output << "}";
      output << ");" << endl;
    }
}

void
SecondDerivExternalFunctionNode::writeJsonExternalFunctionOutput(vector<string> &efout,
                                                                 const temporary_terms_t &temporary_terms,
                                                                 deriv_node_temp_terms_t &tef_terms,
                                                                 const bool isdynamic) const
{
  int second_deriv_symb_id = datatree.external_functions_table.getSecondDerivSymbID(symb_id);
  assert(second_deriv_symb_id != eExtFunSetButNoNameProvided);

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
  if (second_deriv_symb_id == eExtFunNotSet)
    ef << "{\"second_deriv_external_function\": {"
       << "\"external_function_term\": \"TEFDD_fdd_" << getIndxInTefTerms(symb_id, tef_terms) << "_" << inputIndex1 << "_" << inputIndex2 << "\""
       << ", \"analytic_derivative\": false"
       << ", \"wrt1\": " << inputIndex1
       << ", \"wrt2\": " << inputIndex2
       << ", \"value\": \"" << datatree.symbol_table.getName(symb_id) << "(";
  else
    {
      tef_terms[make_pair(second_deriv_symb_id, arguments)] = (int) tef_terms.size();
      ef << "{\"second_deriv_external_function\": {"
         << "\"external_function_term\": \"TEFDD_def_" << getIndxInTefTerms(second_deriv_symb_id, tef_terms) << "\""
         << ", \"analytic_derivative\": true"
         << ", \"value\": \"" << datatree.symbol_table.getName(second_deriv_symb_id) << "(";
    }

  writeJsonExternalFunctionArguments(ef, temporary_terms, tef_terms, isdynamic);
  ef << ")\"}}" << endl;
  efout.push_back(ef.str());
}

expr_t
SecondDerivExternalFunctionNode::cloneDynamic(DataTree &dynamic_datatree) const
{
  vector<expr_t> dynamic_arguments;
  for (vector<expr_t>::const_iterator it = arguments.begin();
       it != arguments.end(); it++)
    dynamic_arguments.push_back((*it)->cloneDynamic(dynamic_datatree));
  return dynamic_datatree.AddSecondDerivExternalFunction(symb_id, dynamic_arguments,
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
  for (vector<expr_t>::const_iterator it = arguments.begin();
       it != arguments.end(); it++)
    static_arguments.push_back((*it)->toStatic(static_datatree));
  return static_datatree.AddSecondDerivExternalFunction(symb_id, static_arguments,
                                                        inputIndex1, inputIndex2);
}

void
SecondDerivExternalFunctionNode::computeXrefs(EquationInfo &ei) const
{
  vector<expr_t> dynamic_arguments;
  for (vector<expr_t>::const_iterator it = arguments.begin();
       it != arguments.end(); it++)
    (*it)->computeXrefs(ei);
}

void
SecondDerivExternalFunctionNode::compile(ostream &CompileCode, unsigned int &instruction_number,
                                         bool lhs_rhs, const temporary_terms_t &temporary_terms,
                                         const map_idx_t &map_idx, bool dynamic, bool steady_dynamic,
                                         deriv_node_temp_terms_t &tef_terms) const
{
  cerr << "SecondDerivExternalFunctionNode::compile: not implemented." << endl;
  exit(EXIT_FAILURE);
}

void
SecondDerivExternalFunctionNode::compileExternalFunctionOutput(ostream &CompileCode, unsigned int &instruction_number,
                                                               bool lhs_rhs, const temporary_terms_t &temporary_terms,
                                                               const map_idx_t &map_idx, bool dynamic, bool steady_dynamic,
                                                               deriv_node_temp_terms_t &tef_terms) const
{
  cerr << "SecondDerivExternalFunctionNode::compileExternalFunctionOutput: not implemented." << endl;
  exit(EXIT_FAILURE);
}

VarExpectationNode::VarExpectationNode(DataTree &datatree_arg,
                                       int symb_id_arg,
                                       int forecast_horizon_arg,
                                       const string &model_name_arg) :
  ExprNode(datatree_arg),
  symb_id(symb_id_arg),
  forecast_horizon(forecast_horizon_arg),
  model_name(model_name_arg),
  yidx(-1)
{
  datatree.var_expectation_node_map[make_pair(model_name, make_pair(symb_id, forecast_horizon))] = this;
}

void
VarExpectationNode::computeTemporaryTerms(map<expr_t, pair<int, NodeTreeReference> > &reference_count,
                                          map<NodeTreeReference, temporary_terms_t> &temp_terms_map,
                                          bool is_matlab, NodeTreeReference tr) const
{
  temp_terms_map[tr].insert(const_cast<VarExpectationNode *>(this));
}

void
VarExpectationNode::computeTemporaryTerms(map<expr_t, int> &reference_count,
                                          temporary_terms_t &temporary_terms,
                                          map<expr_t, pair<int, int> > &first_occurence,
                                          int Curr_block,
                                          vector< vector<temporary_terms_t> > &v_temporary_terms,
                                          int equation) const
{
  expr_t this2 = const_cast<VarExpectationNode *>(this);
  temporary_terms.insert(this2);
  first_occurence[this2] = make_pair(Curr_block, equation);
  v_temporary_terms[Curr_block][equation].insert(this2);
}

expr_t
VarExpectationNode::toStatic(DataTree &static_datatree) const
{
  return static_datatree.AddVariable(symb_id);
}

expr_t
VarExpectationNode::cloneDynamic(DataTree &dynamic_datatree) const
{
  return dynamic_datatree.AddVarExpectation(symb_id, forecast_horizon, model_name);
}

void
VarExpectationNode::writeOutput(ostream &output, ExprNodeOutputType output_type,
                                const temporary_terms_t &temporary_terms,
                                deriv_node_temp_terms_t &tef_terms) const
{
  assert(output_type != oMatlabOutsideModel);

  if (IS_LATEX(output_type))
    {
      output << "VAR_" << model_name << LEFT_PAR(output_type)
             << datatree.symbol_table.getTeXName(symb_id)
             << "_{t+" << forecast_horizon << "}" << RIGHT_PAR(output_type);
      return;
    }

  // If current node is a temporary term
  temporary_terms_t::const_iterator it = temporary_terms.find(const_cast<VarExpectationNode *>(this));
  if (it != temporary_terms.end())
    {
      if (output_type == oMatlabDynamicModelSparse)
        output << "T" << idx << "(it_)";
      else
        output << "T" << idx;
      return;
    }

  output << "dynamic_var_forecast_" << model_name << "_" << forecast_horizon << "(" << yidx + 1 << ")";
}

int
VarExpectationNode::maxEndoLead() const
{
  return 0;
}

int
VarExpectationNode::maxExoLead() const
{
  return 0;
}

int
VarExpectationNode::maxEndoLag() const
{
  return 0;
}

int
VarExpectationNode::maxExoLag() const
{
  return 0;
}

int
VarExpectationNode::maxLead() const
{
  return 0;
}

int
VarExpectationNode::maxLag() const
{
  return 0;
}

expr_t
VarExpectationNode::undiff() const
{
  return const_cast<VarExpectationNode *>(this);
}

int
VarExpectationNode::VarMinLag() const
{
  return 1;
}

void
VarExpectationNode::VarMaxLag(DataTree &static_datatree, set<expr_t> &static_lhs, int &max_lag) const
{
}

int
VarExpectationNode::PacMaxLag(vector<int> &lhs) const
{
  return 0;
}

expr_t
VarExpectationNode::decreaseLeadsLags(int n) const
{
  return const_cast<VarExpectationNode *>(this);
}

void
VarExpectationNode::prepareForDerivation()
{
  preparedForDerivation = true;
  // Come back
}

expr_t
VarExpectationNode::computeDerivative(int deriv_id)
{
  return datatree.Zero;
}

expr_t
VarExpectationNode::getChainRuleDerivative(int deriv_id, const map<int, expr_t> &recursive_variables)
{
  return datatree.Zero;
}

bool
VarExpectationNode::containsExternalFunction() const
{
  return false;
}

double
VarExpectationNode::eval(const eval_context_t &eval_context) const throw (EvalException, EvalExternalFunctionException)
{
  eval_context_t::const_iterator it = eval_context.find(symb_id);
  if (it == eval_context.end())
    throw EvalException();

  return it->second;
}

bool
VarExpectationNode::isDiffPresent() const
{
  return false;
}

void
VarExpectationNode::computeXrefs(EquationInfo &ei) const
{
}

void
VarExpectationNode::collectVARLHSVariable(set<expr_t> &result) const
{
}

void
VarExpectationNode::collectDynamicVariables(SymbolType type_arg, set<pair<int, int> > &result) const
{
}

void
VarExpectationNode::collectTemporary_terms(const temporary_terms_t &temporary_terms, temporary_terms_inuse_t &temporary_terms_inuse, int Curr_Block) const
{
  temporary_terms_t::const_iterator it = temporary_terms.find(const_cast<VarExpectationNode *>(this));
  if (it != temporary_terms.end())
    temporary_terms_inuse.insert(idx);
}

void
VarExpectationNode::compile(ostream &CompileCode, unsigned int &instruction_number,
                            bool lhs_rhs, const temporary_terms_t &temporary_terms,
                            const map_idx_t &map_idx, bool dynamic, bool steady_dynamic,
                            deriv_node_temp_terms_t &tef_terms) const
{
  cerr << "VarExpectationNode::compile not implemented." << endl;
  exit(EXIT_FAILURE);
}

pair<int, expr_t >
VarExpectationNode::normalizeEquation(int var_endo, vector<pair<int, pair<expr_t, expr_t> > > &List_of_Op_RHS) const
{
  return make_pair(0, datatree.AddVariableInternal(symb_id, 0));
}

expr_t
VarExpectationNode::substituteEndoLeadGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const
{
  return const_cast<VarExpectationNode *>(this);
}

expr_t
VarExpectationNode::substituteEndoLagGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  return const_cast<VarExpectationNode *>(this);
}

expr_t
VarExpectationNode::substituteExoLead(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const
{
  return const_cast<VarExpectationNode *>(this);
}

expr_t
VarExpectationNode::substituteExoLag(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  return const_cast<VarExpectationNode *>(this);
}

expr_t
VarExpectationNode::substituteExpectation(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool partial_information_model) const
{
  return const_cast<VarExpectationNode *>(this);
}

expr_t
VarExpectationNode::substituteAdl() const
{
  return const_cast<VarExpectationNode *>(this);
}

void
VarExpectationNode::findDiffNodes(DataTree &static_datatree, diff_table_t &diff_table) const
{
}

void
VarExpectationNode::findDiffUnaryOpNodes(DataTree &static_datatree, set<UnaryOpNode *> &unary_op_nodes_in_diff) const
{
}

expr_t
VarExpectationNode::substituteDiff(DataTree &static_datatree, diff_table_t &diff_table, subst_table_t &subst_table,
                                   vector<BinaryOpNode *> &neweqs) const
{
  return const_cast<VarExpectationNode *>(this);
}

expr_t
VarExpectationNode::substituteDiffUnaryOpNodes(DataTree &static_datatree, set<UnaryOpNode *> &nodes, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  return const_cast<VarExpectationNode *>(this);
}

expr_t
VarExpectationNode::substitutePacExpectation(map<const PacExpectationNode *, const BinaryOpNode *> &subst_table)
{
  return const_cast<VarExpectationNode *>(this);
}

expr_t
VarExpectationNode::differentiateForwardVars(const vector<string> &subset, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  return const_cast<VarExpectationNode *>(this);
}

bool
VarExpectationNode::containsEndogenous(void) const
{
  return true;
}

bool
VarExpectationNode::containsExogenous() const
{
  return false;
}

bool
VarExpectationNode::isNumConstNodeEqualTo(double value) const
{
  return false;
}

expr_t
VarExpectationNode::decreaseLeadsLagsPredeterminedVariables() const
{
  return const_cast<VarExpectationNode *>(this);
}

bool
VarExpectationNode::isVariableNodeEqualTo(SymbolType type_arg, int variable_id, int lag_arg) const
{
  return false;
}

expr_t
VarExpectationNode::replaceTrendVar() const
{
  return const_cast<VarExpectationNode *>(this);
}

expr_t
VarExpectationNode::detrend(int symb_id, bool log_trend, expr_t trend) const
{
  return const_cast<VarExpectationNode *>(this);
}

expr_t
VarExpectationNode::removeTrendLeadLag(map<int, expr_t> trend_symbols_map) const
{
  return const_cast<VarExpectationNode *>(this);
}

bool
VarExpectationNode::isInStaticForm() const
{
  return false;
}

bool
VarExpectationNode::isVarModelReferenced(const string &model_info_name) const
{
  return model_name == model_info_name;
}

void
VarExpectationNode::getEndosAndMaxLags(map<string, int> &model_endos_and_lags) const
{
}

void
VarExpectationNode::setVarExpectationIndex(map<string, pair<SymbolList, int> > &var_model_info)
{
  vector<string> vs = var_model_info[model_name].first.get_symbols();
  yidx = find(vs.begin(), vs.end(), datatree.symbol_table.getName(symb_id)) - vs.begin();
}

void
VarExpectationNode::walkPacParameters(bool &pac_encountered, pair<int, int> &lhs, set<pair<int, pair<int, int> > > &ec_params_and_vars, set<pair<int, pair<int, int> > > &ar_params_and_vars) const
{
}

void
VarExpectationNode::addParamInfoToPac(pair<int, int> &lhs_arg, set<pair<int, pair<int, int> > > &ec_params_and_vars_arg, set<pair<int, pair<int, int> > > &ar_params_and_vars_arg)
{
}

void
VarExpectationNode::fillPacExpectationVarInfo(string &model_name_arg, vector<int> &lhs_arg, int max_lag_arg, vector<bool> &nonstationary_arg, int growth_symb_id_arg, int equation_number_arg)
{
}

expr_t
VarExpectationNode::substituteStaticAuxiliaryVariable() const
{
  return const_cast<VarExpectationNode *>(this);
}

void
VarExpectationNode::writeJsonOutput(ostream &output,
                                    const temporary_terms_t &temporary_terms,
                                    deriv_node_temp_terms_t &tef_terms,
                                    const bool isdynamic) const
{
  output << "var_expectation("
         << "forecast_horizon = " << forecast_horizon
         << ", name = " << datatree.symbol_table.getName(symb_id)
         << ", model_name = " << model_name
         << ", yindex = " << yidx
         << ")";
}

PacExpectationNode::PacExpectationNode(DataTree &datatree_arg,
                                       const string &model_name_arg) :
  ExprNode(datatree_arg),
  model_name(model_name_arg)
{
  datatree.pac_expectation_node_map[model_name] = this;
}

void
PacExpectationNode::computeTemporaryTerms(map<expr_t, pair<int, NodeTreeReference> > &reference_count,
                                          map<NodeTreeReference, temporary_terms_t> &temp_terms_map,
                                          bool is_matlab, NodeTreeReference tr) const
{
  temp_terms_map[tr].insert(const_cast<PacExpectationNode *>(this));
}

void
PacExpectationNode::computeTemporaryTerms(map<expr_t, int> &reference_count,
                                          temporary_terms_t &temporary_terms,
                                          map<expr_t, pair<int, int> > &first_occurence,
                                          int Curr_block,
                                          vector< vector<temporary_terms_t> > &v_temporary_terms,
                                          int equation) const
{
  expr_t this2 = const_cast<PacExpectationNode *>(this);
  temporary_terms.insert(this2);
  first_occurence[this2] = make_pair(Curr_block, equation);
  v_temporary_terms[Curr_block][equation].insert(this2);
}

expr_t
PacExpectationNode::toStatic(DataTree &static_datatree) const
{
  return static_datatree.AddPacExpectation(string(model_name));
}

expr_t
PacExpectationNode::cloneDynamic(DataTree &dynamic_datatree) const
{
  return dynamic_datatree.AddPacExpectation(string(model_name));
}

void
PacExpectationNode::writeOutput(ostream &output, ExprNodeOutputType output_type,
                                const temporary_terms_t &temporary_terms,
                                deriv_node_temp_terms_t &tef_terms) const
{
  assert(output_type != oMatlabOutsideModel);

  if (IS_LATEX(output_type))
    {
      output << "PAC_EXPECTATION" << LEFT_PAR(output_type) << model_name << RIGHT_PAR(output_type);
      return;
    }

  output << "M_.pac." << model_name << ".lhs_var = "
         << datatree.symbol_table.getTypeSpecificID(lhs_pac_var.first) + 1 << ";" << endl
         << "M_.pac." << model_name << ".max_lag = " << max_lag << ";" << endl;

  if (growth_symb_id >= 0)
    output << "M_.pac." << model_name << ".growth_neutrality_param_index = "
           << datatree.symbol_table.getTypeSpecificID(growth_param_index) + 1 << ";" << endl;

  output << "M_.pac." << model_name << ".ec.params = ";
  set<pair<int, pair<int, int> > >::const_iterator it = ec_params_and_vars.begin();
  output << datatree.symbol_table.getTypeSpecificID(it->first) + 1
         << ";" << endl
         << "M_.pac." << model_name << ".ec.vars = [";
  for (set<pair<int, pair<int, int> > >::const_iterator it = ec_params_and_vars.begin();
       it != ec_params_and_vars.end(); it++)
    {
      if (it != ec_params_and_vars.begin())
        output << " ";
      output << datatree.symbol_table.getTypeSpecificID(it->second.first) + 1;
    }
  output << "];" << endl
         << "M_.pac." << model_name << ".ar.params = [";
  for (set<pair<int, pair<int, int> > >::const_iterator it = ar_params_and_vars.begin();
       it != ar_params_and_vars.end(); it++)
    {
      if (it != ar_params_and_vars.begin())
        output << " ";
      output << datatree.symbol_table.getTypeSpecificID(it->first) + 1;
    }
  output << "];" << endl
         << "M_.pac." << model_name << ".ar.vars = [";
  for (set<pair<int, pair<int, int> > >::const_iterator it = ar_params_and_vars.begin();
       it != ar_params_and_vars.end(); it++)
    {
      if (it != ar_params_and_vars.begin())
        output << " ";
      output << datatree.symbol_table.getTypeSpecificID(it->second.first) + 1;
    }
  output << "];" << endl
         << "M_.pac." << model_name << ".ar.lags = [";
  for (set<pair<int, pair<int, int> > >::const_iterator it = ar_params_and_vars.begin();
       it != ar_params_and_vars.end(); it++)
    {
      if (it != ar_params_and_vars.begin())
        output << " ";
      output << it->second.second;
    }
  output << "];" << endl
         << "M_.pac." << model_name << ".h0_param_indices = [";
  for (vector<int>::const_iterator it = h0_indices.begin();
       it != h0_indices.end(); it++)
    {
      if (it != h0_indices.begin())
        output << " ";
      output << datatree.symbol_table.getTypeSpecificID(*it) + 1;
    }
  output << "];" << endl
         << "M_.pac." << model_name << ".h1_param_indices = [";
  for (vector<int>::const_iterator it = h1_indices.begin();
       it != h1_indices.end(); it++)
    {
      if (it != h1_indices.begin())
        output << " ";
      output << datatree.symbol_table.getTypeSpecificID(*it) + 1;
    }
  output << "];" << endl;
}

int
PacExpectationNode::maxEndoLead() const
{
  return 0;
}

int
PacExpectationNode::maxExoLead() const
{
  return 0;
}

int
PacExpectationNode::maxEndoLag() const
{
  return 0;
}

int
PacExpectationNode::maxExoLag() const
{
  return 0;
}

int
PacExpectationNode::maxLead() const
{
  return 0;
}

int
PacExpectationNode::maxLag() const
{
  return 0;
}

expr_t
PacExpectationNode::undiff() const
{
  return const_cast<PacExpectationNode *>(this);
}

int
PacExpectationNode::VarMinLag() const
{
  return 1;
}

void
PacExpectationNode::VarMaxLag(DataTree &static_datatree, set<expr_t> &static_lhs, int &max_lag) const
{
}

int
PacExpectationNode::PacMaxLag(vector<int> &lhs) const
{
  return 0;
}

expr_t
PacExpectationNode::decreaseLeadsLags(int n) const
{
  return const_cast<PacExpectationNode *>(this);
}

void
PacExpectationNode::prepareForDerivation()
{
  cerr << "PacExpectationNode::prepareForDerivation: shouldn't arrive here." << endl;
  exit(EXIT_FAILURE);
}

expr_t
PacExpectationNode::computeDerivative(int deriv_id)
{
  cerr << "PacExpectationNode::computeDerivative: shouldn't arrive here." << endl;
  exit(EXIT_FAILURE);
}

expr_t
PacExpectationNode::getChainRuleDerivative(int deriv_id, const map<int, expr_t> &recursive_variables)
{
  cerr << "PacExpectationNode::getChainRuleDerivative: shouldn't arrive here." << endl;
  exit(EXIT_FAILURE);
}

bool
PacExpectationNode::containsExternalFunction() const
{
  return false;
}

double
PacExpectationNode::eval(const eval_context_t &eval_context) const throw (EvalException, EvalExternalFunctionException)
{
  throw EvalException();
}

void
PacExpectationNode::computeXrefs(EquationInfo &ei) const
{
}

void
PacExpectationNode::collectVARLHSVariable(set<expr_t> &result) const
{
}

void
PacExpectationNode::collectDynamicVariables(SymbolType type_arg, set<pair<int, int> > &result) const
{
}

void
PacExpectationNode::collectTemporary_terms(const temporary_terms_t &temporary_terms, temporary_terms_inuse_t &temporary_terms_inuse, int Curr_Block) const
{
  temporary_terms_t::const_iterator it = temporary_terms.find(const_cast<PacExpectationNode *>(this));
  if (it != temporary_terms.end())
    temporary_terms_inuse.insert(idx);
}

void
PacExpectationNode::compile(ostream &CompileCode, unsigned int &instruction_number,
                            bool lhs_rhs, const temporary_terms_t &temporary_terms,
                            const map_idx_t &map_idx, bool dynamic, bool steady_dynamic,
                            deriv_node_temp_terms_t &tef_terms) const
{
  cerr << "PacExpectationNode::compile not implemented." << endl;
  exit(EXIT_FAILURE);
}

bool
PacExpectationNode::isDiffPresent() const
{
  return false;
}

pair<int, expr_t >
PacExpectationNode::normalizeEquation(int var_endo, vector<pair<int, pair<expr_t, expr_t> > > &List_of_Op_RHS) const
{
  //COME BACK
  return make_pair(0, const_cast<PacExpectationNode *>(this));
}

expr_t
PacExpectationNode::substituteEndoLeadGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const
{
  return const_cast<PacExpectationNode *>(this);
}

expr_t
PacExpectationNode::substituteEndoLagGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  return const_cast<PacExpectationNode *>(this);
}

expr_t
PacExpectationNode::substituteExoLead(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const
{
  return const_cast<PacExpectationNode *>(this);
}

expr_t
PacExpectationNode::substituteExoLag(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  return const_cast<PacExpectationNode *>(this);
}

expr_t
PacExpectationNode::substituteExpectation(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool partial_information_model) const
{
  return const_cast<PacExpectationNode *>(this);
}

expr_t
PacExpectationNode::substituteAdl() const
{
  return const_cast<PacExpectationNode *>(this);
}

void
PacExpectationNode::findDiffNodes(DataTree &static_datatree, diff_table_t &diff_table) const
{
}

void
PacExpectationNode::findDiffUnaryOpNodes(DataTree &static_datatree, set<UnaryOpNode *> &unary_op_nodes_in_diff) const
{
}

expr_t
PacExpectationNode::substituteDiff(DataTree &static_datatree, diff_table_t &diff_table, subst_table_t &subst_table,
                                   vector<BinaryOpNode *> &neweqs) const
{
  return const_cast<PacExpectationNode *>(this);
}

expr_t
PacExpectationNode::substituteDiffUnaryOpNodes(DataTree &static_datatree, set<UnaryOpNode *> &nodes, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  return const_cast<PacExpectationNode *>(this);
}

expr_t
PacExpectationNode::differentiateForwardVars(const vector<string> &subset, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  return const_cast<PacExpectationNode *>(this);
}

bool
PacExpectationNode::containsEndogenous(void) const
{
  return true;
}

bool
PacExpectationNode::containsExogenous() const
{
  return false;
}

bool
PacExpectationNode::isNumConstNodeEqualTo(double value) const
{
  return false;
}

expr_t
PacExpectationNode::decreaseLeadsLagsPredeterminedVariables() const
{
  return const_cast<PacExpectationNode *>(this);
}

bool
PacExpectationNode::isVariableNodeEqualTo(SymbolType type_arg, int variable_id, int lag_arg) const
{
  return false;
}

expr_t
PacExpectationNode::replaceTrendVar() const
{
  return const_cast<PacExpectationNode *>(this);
}

expr_t
PacExpectationNode::detrend(int symb_id, bool log_trend, expr_t trend) const
{
  return const_cast<PacExpectationNode *>(this);
}

expr_t
PacExpectationNode::removeTrendLeadLag(map<int, expr_t> trend_symbols_map) const
{
  return const_cast<PacExpectationNode *>(this);
}

bool
PacExpectationNode::isInStaticForm() const
{
  return false;
}

bool
PacExpectationNode::isVarModelReferenced(const string &model_info_name) const
{
  return model_name == model_info_name;
}

void
PacExpectationNode::getEndosAndMaxLags(map<string, int> &model_endos_and_lags) const
{
}

void
PacExpectationNode::setVarExpectationIndex(map<string, pair<SymbolList, int> > &var_model_info)
{
}

expr_t
PacExpectationNode::substituteStaticAuxiliaryVariable() const
{
  return const_cast<PacExpectationNode *>(this);
}

void
PacExpectationNode::writeJsonOutput(ostream &output,
                                    const temporary_terms_t &temporary_terms,
                                    deriv_node_temp_terms_t &tef_terms,
                                    const bool isdynamic) const
{
  output << "pac_expectation("
         << "model_name = " << model_name
         << ")";
}

void
PacExpectationNode::walkPacParameters(bool &pac_encountered, pair<int, int> &lhs, set<pair<int, pair<int, int> > > &ec_params_and_vars, set<pair<int, pair<int, int> > > &ar_params_and_vars) const
{
  pac_encountered = true;
}

void
PacExpectationNode::addParamInfoToPac(pair<int, int> &lhs_arg, set<pair<int, pair<int, int> > > &ec_params_and_vars_arg, set<pair<int, pair<int, int> > > &ar_params_and_vars_arg)
{
  if (lhs_arg.first == -1)
    {
      cerr << "Pac Expectation: error in obtaining LHS varibale." << endl;
      exit(EXIT_FAILURE);
    }

  if (ec_params_and_vars_arg.empty() || ar_params_and_vars_arg.empty())
    {
      cerr << "Pac Expectation: error in obtaining RHS parameters." << endl;
      exit(EXIT_FAILURE);
    }

  lhs_pac_var = lhs_arg;
  ar_params_and_vars = ar_params_and_vars_arg;
  ec_params_and_vars = ec_params_and_vars_arg;
}


void
PacExpectationNode::fillPacExpectationVarInfo(string &model_name_arg, vector<int> &lhs_arg, int max_lag_arg, vector<bool> &nonstationary_arg, int growth_symb_id_arg, int equation_number_arg)
{
  if (model_name != model_name_arg)
    return;

  lhs = lhs_arg;
  max_lag = max_lag_arg;
  growth_symb_id = growth_symb_id_arg;
  equation_number = equation_number_arg;

  for (vector<bool>::const_iterator it = nonstationary_arg.begin();
       it != nonstationary_arg.end(); it++)
    {
      if (*it)
        nonstationary_vars_present = true;
      else
        stationary_vars_present = true;
      if (nonstationary_vars_present && stationary_vars_present)
        break;
    }
}

expr_t
PacExpectationNode::substitutePacExpectation(map<const PacExpectationNode *, const BinaryOpNode *> &subst_table)
{
  map<const PacExpectationNode *, const BinaryOpNode *>::const_iterator myit =
    subst_table.find(const_cast<PacExpectationNode *>(this));
  if (myit != subst_table.end())
    return const_cast<BinaryOpNode *>(myit->second);

  expr_t subExpr = datatree.AddNonNegativeConstant("0");
  if (stationary_vars_present)
    for (int i = 1; i < max_lag + 1; i++)
      for (vector<int>::const_iterator it = lhs.begin(); it != lhs.end(); it++)
        {
          stringstream param_name_h0;
          param_name_h0 << "h0_" << model_name
                        << "_var_" << datatree.symbol_table.getName(*it)
                        << "_lag_" << i;
          int new_param_symb_id = datatree.symbol_table.addSymbol(param_name_h0.str(), eParameter);
          h0_indices.push_back(new_param_symb_id);
          subExpr = datatree.AddPlus(subExpr,
                                     datatree.AddTimes(datatree.AddVariable(new_param_symb_id),
                                                       datatree.AddVariable(*it, -i)));
        }

  if (nonstationary_vars_present)
    for (int i = 1; i < max_lag + 1; i++)
      for (vector<int>::const_iterator it = lhs.begin(); it != lhs.end(); it++)
        {
          stringstream param_name_h1;
          param_name_h1 << "h1_" << model_name
                        << "_var_" << datatree.symbol_table.getName(*it)
                        << "_lag_" << i;
          int new_param_symb_id = datatree.symbol_table.addSymbol(param_name_h1.str(), eParameter);
          h1_indices.push_back(new_param_symb_id);
          subExpr = datatree.AddPlus(subExpr,
                                     datatree.AddTimes(datatree.AddVariable(new_param_symb_id),
                                                       datatree.AddVariable(*it, -i)));
        }

  if (growth_symb_id >= 0)
    {
      growth_param_index = datatree.symbol_table.addSymbol(model_name +
                                                           "_pac_growth_neutrality_correction",
                                                           eParameter);
      subExpr = datatree.AddPlus(subExpr,
                                 datatree.AddTimes(datatree.AddVariable(growth_param_index),
                                                   datatree.AddVariable(growth_symb_id)));
    }

  subst_table[const_cast<PacExpectationNode *>(this)] = dynamic_cast<BinaryOpNode *>(subExpr);

  return subExpr;
}
