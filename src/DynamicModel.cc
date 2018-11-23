/*
 * Copyright (C) 2003-2018 Dynare Team
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
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <cstdio>
#include <cerrno>
#include <algorithm>
#include <iterator>
#include <numeric>

#include "DynamicModel.hh"

void
DynamicModel::copyHelper(const DynamicModel &m)
{
  auto f = [this](const ExprNode *e) { return e->clone(*this); };

  for (const auto &it : m.static_only_equations)
    static_only_equations.push_back(dynamic_cast<BinaryOpNode *>(f(it)));

  auto convert_vector_tt = [f](vector<temporary_terms_t> vtt)
    {
      vector<temporary_terms_t> vtt2;
      for (const auto &tt : vtt)
        {
          temporary_terms_t tt2;
          for (const auto &it : tt)
            tt2.insert(f(it));
          vtt2.push_back(tt2);
        }
      return vtt2;
    };

  for (const auto &it : m.v_temporary_terms)
    v_temporary_terms.push_back(convert_vector_tt(it));

  for (const auto &it : m.first_chain_rule_derivatives)
    first_chain_rule_derivatives[it.first] = f(it.second);

  for (const auto &it : m.equation_type_and_normalized_equation)
    equation_type_and_normalized_equation.push_back(make_pair(it.first, f(it.second)));

  for (const auto &it : m.blocks_derivatives)
    {
      block_derivatives_equation_variable_laglead_nodeid_t v;
      for (const auto &it2 : it)
        v.emplace_back(get<0>(it2), get<1>(it2), get<2>(it2), f(get<3>(it2)));
      blocks_derivatives.push_back(v);
    }

  for (const auto &it : m.dynamic_jacobian)
    dynamic_jacobian[it.first] = f(it.second);

  auto convert_derivative_t = [f](derivative_t dt)
    {
      derivative_t dt2;
      for (const auto &it : dt)
        dt2[it.first] = f(it.second);
      return dt2;
    };
  for (const auto &it : m.derivative_endo)
    derivative_endo.push_back(convert_derivative_t(it));
  for (const auto &it : m.derivative_other_endo)
    derivative_other_endo.push_back(convert_derivative_t(it));
  for (const auto &it : m.derivative_exo)
    derivative_exo.push_back(convert_derivative_t(it));
  for (const auto &it : m.derivative_exo_det)
    derivative_exo_det.push_back(convert_derivative_t(it));

  for (const auto &it : m.pac_expectation_info)
    pac_expectation_info.insert(dynamic_cast<const PacExpectationNode *>(f(it)));
}


DynamicModel::DynamicModel(SymbolTable &symbol_table_arg,
                           NumericalConstants &num_constants_arg,
                           ExternalFunctionsTable &external_functions_table_arg,
                           TrendComponentModelTable &trend_component_model_table_arg,
                           VarModelTable &var_model_table_arg) :
  ModelTree {symbol_table_arg, num_constants_arg, external_functions_table_arg, true},
  trend_component_model_table{trend_component_model_table_arg},
  var_model_table{var_model_table_arg}
{
}

DynamicModel::DynamicModel(const DynamicModel &m) :
  ModelTree {m},
  trend_component_model_table {m.trend_component_model_table},
  var_model_table {m.var_model_table},
  static_only_equations_lineno {m.static_only_equations_lineno},
  static_only_equations_equation_tags {m.static_only_equations_equation_tags},
  deriv_id_table {m.deriv_id_table},
  inv_deriv_id_table {m.inv_deriv_id_table},
  dyn_jacobian_cols_table {m.dyn_jacobian_cols_table},
  max_lag {m.max_lag},
  max_lead {m.max_lead},
  max_endo_lag {m.max_endo_lag},
  max_endo_lead {m.max_endo_lead},
  max_exo_lag {m.max_exo_lag},
  max_exo_lead {m.max_exo_lead},
  max_exo_det_lag {m.max_exo_det_lag},
  max_exo_det_lead {m.max_exo_det_lead},
  max_lag_orig {m.max_lag_orig},
  max_lead_orig {m.max_lead_orig},
  max_endo_lag_orig {m.max_endo_lag_orig},
  max_endo_lead_orig {m.max_endo_lead_orig},
  max_exo_lag_orig {m.max_exo_lag_orig},
  max_exo_lead_orig {m.max_exo_lead_orig},
  max_exo_det_lag_orig {m.max_exo_det_lag_orig},
  max_exo_det_lead_orig {m.max_exo_det_lead_orig},
  xrefs {m.xrefs},
  xref_param  {m.xref_param},
  xref_endo {m.xref_endo},
  xref_exo {m.xref_exo},
  xref_exo_det {m.xref_exo_det},
  nonzero_hessian_eqs {m.nonzero_hessian_eqs},
  v_temporary_terms_inuse {m.v_temporary_terms_inuse},
  map_idx {m.map_idx},
  global_temporary_terms {m.global_temporary_terms},
  block_type_firstequation_size_mfs {m.block_type_firstequation_size_mfs},
  blocks_linear {m.blocks_linear},
  other_endo_block {m.other_endo_block},
  exo_block {m.exo_block},
  exo_det_block {m.exo_det_block},
  block_var_exo {m.block_var_exo},
  block_exo_index {m.block_exo_index},
  block_det_exo_index {m.block_det_exo_index},
  block_other_endo_index {m.block_other_endo_index},
  block_col_type {m.block_col_type},
  variable_block_lead_lag {m.variable_block_lead_lag},
  equation_block {m.equation_block},
  var_expectation_functions_to_write {m.var_expectation_functions_to_write},
  endo_max_leadlag_block {m.endo_max_leadlag_block},
  other_endo_max_leadlag_block {m.other_endo_max_leadlag_block},
  exo_max_leadlag_block {m.exo_max_leadlag_block},
  exo_det_max_leadlag_block {m.exo_det_max_leadlag_block},
  max_leadlag_block {m.max_leadlag_block}
{
  copyHelper(m);
}

DynamicModel &
DynamicModel::operator=(const DynamicModel &m)
{
  ModelTree::operator=(m);

  assert(&trend_component_model_table == &m.trend_component_model_table);
  assert(&var_model_table == &m.var_model_table);

  static_only_equations_lineno = m.static_only_equations_lineno;
  static_only_equations_equation_tags = m.static_only_equations_equation_tags;
  deriv_id_table = m.deriv_id_table;
  inv_deriv_id_table = m.inv_deriv_id_table;
  dyn_jacobian_cols_table = m.dyn_jacobian_cols_table;
  max_lag = m.max_lag;
  max_lead = m.max_lead;
  max_endo_lag = m.max_endo_lag;
  max_endo_lead = m.max_endo_lead;
  max_exo_lag = m.max_exo_lag;
  max_exo_lead = m.max_exo_lead;
  max_exo_det_lag = m.max_exo_det_lag;
  max_exo_det_lead = m.max_exo_det_lead;
  max_lag_orig = m.max_lag_orig;
  max_lead_orig = m.max_lead_orig;
  max_endo_lag_orig = m.max_endo_lag_orig;
  max_endo_lead_orig = m.max_endo_lead_orig;
  max_exo_lag_orig = m.max_exo_lag_orig;
  max_exo_lead_orig = m.max_exo_lead_orig;
  max_exo_det_lag_orig = m.max_exo_det_lag_orig;
  max_exo_det_lead_orig = m.max_exo_det_lead_orig;
  xrefs = m.xrefs;
  xref_param  = m.xref_param;
  xref_endo = m.xref_endo;
  xref_exo = m.xref_exo;
  xref_exo_det = m.xref_exo_det;
  nonzero_hessian_eqs = m.nonzero_hessian_eqs;

  v_temporary_terms.clear();

  v_temporary_terms_inuse = m.v_temporary_terms_inuse;

  first_chain_rule_derivatives.clear();

  map_idx = m.map_idx;
  global_temporary_terms = m.global_temporary_terms;

  equation_type_and_normalized_equation.clear();

  block_type_firstequation_size_mfs = m.block_type_firstequation_size_mfs;

  blocks_derivatives.clear();
  dynamic_jacobian.clear();

  blocks_linear = m.blocks_linear;

  derivative_endo.clear();
  derivative_other_endo.clear();
  derivative_exo.clear();
  derivative_exo_det.clear();

  other_endo_block = m.other_endo_block;
  exo_block = m.exo_block;
  exo_det_block = m.exo_det_block;
  block_var_exo = m.block_var_exo;
  block_exo_index = m.block_exo_index;
  block_det_exo_index = m.block_det_exo_index;
  block_other_endo_index = m.block_other_endo_index;
  block_col_type = m.block_col_type;
  variable_block_lead_lag = m.variable_block_lead_lag;
  equation_block = m.equation_block;
  var_expectation_functions_to_write = m.var_expectation_functions_to_write;

  pac_expectation_info.clear();

  endo_max_leadlag_block = m.endo_max_leadlag_block;
  other_endo_max_leadlag_block = m.other_endo_max_leadlag_block;
  exo_max_leadlag_block = m.exo_max_leadlag_block;
  exo_det_max_leadlag_block = m.exo_det_max_leadlag_block;
  max_leadlag_block = m.max_leadlag_block;

  copyHelper(m);

  return *this;
}

void
DynamicModel::compileDerivative(ofstream &code_file, unsigned int &instruction_number, int eq, int symb_id, int lag, const map_idx_t &map_idx) const
{
  auto it = derivatives[1].find({ eq, getDerivID(symbol_table.getID(SymbolType::endogenous, symb_id), lag) });
  if (it != derivatives[1].end())
    (it->second)->compile(code_file, instruction_number, false, temporary_terms, map_idx, true, false);
  else
    {
      FLDZ_ fldz;
      fldz.write(code_file, instruction_number);
    }
}

void
DynamicModel::compileChainRuleDerivative(ofstream &code_file, unsigned int &instruction_number, int eqr, int varr, int lag, const map_idx_t &map_idx) const
{
  auto it = first_chain_rule_derivatives.find({ eqr, varr, lag });
  if (it != first_chain_rule_derivatives.end())
    (it->second)->compile(code_file, instruction_number, false, temporary_terms, map_idx, true, false);
  else
    {
      FLDZ_ fldz;
      fldz.write(code_file, instruction_number);
    }
}

void
DynamicModel::computeTemporaryTermsOrdered()
{
  map<expr_t, pair<int, int>> first_occurence;
  map<expr_t, int> reference_count;
  BinaryOpNode *eq_node;
  first_chain_rule_derivatives_t::const_iterator it_chr;
  ostringstream tmp_s;
  v_temporary_terms.clear();
  map_idx.clear();

  unsigned int nb_blocks = getNbBlocks();
  v_temporary_terms = vector<vector<temporary_terms_t>>(nb_blocks);
  v_temporary_terms_inuse = vector<temporary_terms_inuse_t>(nb_blocks);
  temporary_terms.clear();

  if (!global_temporary_terms)
    {
      for (unsigned int block = 0; block < nb_blocks; block++)
        {
          reference_count.clear();
          temporary_terms.clear();
          unsigned int block_size = getBlockSize(block);
          unsigned int block_nb_mfs = getBlockMfs(block);
          unsigned int block_nb_recursives = block_size - block_nb_mfs;
          v_temporary_terms[block] = vector<temporary_terms_t>(block_size);
          for (unsigned int i = 0; i < block_size; i++)
            {
              if (i < block_nb_recursives && isBlockEquationRenormalized(block, i))
                getBlockEquationRenormalizedExpr(block, i)->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, block, v_temporary_terms,  i);
              else
                {
                  eq_node = (BinaryOpNode *) getBlockEquationExpr(block, i);
                  eq_node->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, block, v_temporary_terms,  i);
                }
            }
          for (const auto &it : blocks_derivatives[block])
            {
              expr_t id = get<3>(it);
              id->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, block, v_temporary_terms,  block_size-1);
            }
          for (const auto &it : derivative_endo[block])
            it.second->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, block, v_temporary_terms,  block_size-1);
          for (const auto &it : derivative_other_endo[block])
            it.second->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, block, v_temporary_terms,  block_size-1);
          v_temporary_terms_inuse[block] = {};
        }
    }
  else
    {
      for (unsigned int block = 0; block < nb_blocks; block++)
        {
          // Compute the temporary terms reordered
          unsigned int block_size = getBlockSize(block);
          unsigned int block_nb_mfs = getBlockMfs(block);
          unsigned int block_nb_recursives = block_size - block_nb_mfs;
          v_temporary_terms[block] = vector<temporary_terms_t>(block_size);
          for (unsigned int i = 0; i < block_size; i++)
            {
              if (i < block_nb_recursives && isBlockEquationRenormalized(block, i))
                getBlockEquationRenormalizedExpr(block, i)->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, block, v_temporary_terms,  i);
              else
                {
                  eq_node = (BinaryOpNode *) getBlockEquationExpr(block, i);
                  eq_node->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, block, v_temporary_terms, i);
                }
            }
          for (const auto &it : blocks_derivatives[block])
            {
              expr_t id = get<3>(it);
              id->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, block, v_temporary_terms, block_size-1);
            }
          for (const auto &it : derivative_endo[block])
            it.second->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, block, v_temporary_terms, block_size-1);
          for (const auto &it : derivative_other_endo[block])
            it.second->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, block, v_temporary_terms, block_size-1);
        }
      for (unsigned int block = 0; block < nb_blocks; block++)
        {
          // Collect the temporary terms reordered
          unsigned int block_size = getBlockSize(block);
          unsigned int block_nb_mfs = getBlockMfs(block);
          unsigned int block_nb_recursives = block_size - block_nb_mfs;
          set<int> temporary_terms_in_use;
          for (unsigned int i = 0; i < block_size; i++)
            {
              if (i < block_nb_recursives && isBlockEquationRenormalized(block, i))
                getBlockEquationRenormalizedExpr(block, i)->collectTemporary_terms(temporary_terms, temporary_terms_in_use, block);
              else
                {
                  eq_node = (BinaryOpNode *) getBlockEquationExpr(block, i);
                  eq_node->collectTemporary_terms(temporary_terms, temporary_terms_in_use, block);
                }
            }
          for (const auto &it : blocks_derivatives[block])
            {
              expr_t id = get<3>(it);
              id->collectTemporary_terms(temporary_terms, temporary_terms_in_use, block);
            }
          for (const auto &it : derivative_endo[block])
            it.second->collectTemporary_terms(temporary_terms, temporary_terms_in_use, block);
          for (const auto &it : derivative_other_endo[block])
            it.second->collectTemporary_terms(temporary_terms, temporary_terms_in_use, block);
          for (const auto &it : derivative_exo[block])
            it.second->collectTemporary_terms(temporary_terms, temporary_terms_in_use, block);
          for (const auto &it : derivative_exo_det[block])
            it.second->collectTemporary_terms(temporary_terms, temporary_terms_in_use, block);
          v_temporary_terms_inuse[block] = temporary_terms_in_use;
        }
      computeTemporaryTermsMapping();
    }
}

void
DynamicModel::computeTemporaryTermsMapping()
{
  // Add a mapping form node ID to temporary terms order
  int j = 0;
  for (auto temporary_term : temporary_terms)
    map_idx[temporary_term->idx] = j++;
}

void
DynamicModel::writeModelEquationsOrdered_M(const string &basename) const
{
  string tmp_s, sps;
  ostringstream tmp_output, tmp1_output, global_output;
  expr_t lhs = nullptr, rhs = nullptr;
  BinaryOpNode *eq_node;
  ostringstream Ufoss;
  vector<string> Uf(symbol_table.endo_nbr(), "");
  map<expr_t, int> reference_count;
  temporary_terms_t local_temporary_terms;
  ofstream  output;
  int nze, nze_exo, nze_exo_det, nze_other_endo;
  vector<int> feedback_variables;
  ExprNodeOutputType local_output_type;
  Ufoss.str("");

  local_output_type = ExprNodeOutputType::matlabDynamicModelSparse;
  if (global_temporary_terms)
    local_temporary_terms = temporary_terms;

  //----------------------------------------------------------------------
  //For each block
  for (unsigned int block = 0; block < getNbBlocks(); block++)
    {

      //recursive_variables.clear();
      feedback_variables.clear();
      //For a block composed of a single equation determines wether we have to evaluate or to solve the equation
      nze = derivative_endo[block].size();
      nze_other_endo = derivative_other_endo[block].size();
      nze_exo = derivative_exo[block].size();
      nze_exo_det = derivative_exo_det[block].size();
      BlockSimulationType simulation_type = getBlockSimulationType(block);
      unsigned int block_size = getBlockSize(block);
      unsigned int block_mfs = getBlockMfs(block);
      unsigned int block_recursive = block_size - block_mfs;
      deriv_node_temp_terms_t tef_terms;
      local_output_type = ExprNodeOutputType::matlabDynamicModelSparse;
      if (global_temporary_terms)
        local_temporary_terms = temporary_terms;

      int prev_lag;
      unsigned int prev_var, count_col, count_col_endo, count_col_exo, count_col_exo_det, count_col_other_endo;
      map<tuple<int, int, int>, expr_t> tmp_block_endo_derivative;
      for (const auto &it : blocks_derivatives[block])
        tmp_block_endo_derivative[{ get<2>(it), get<1>(it), get<0>(it) }] = get<3>(it);
      prev_var = 999999999;
      prev_lag = -9999999;
      count_col_endo = 0;
      for (const auto &it : tmp_block_endo_derivative)
        {
          int lag = get<0>(it.first);
          unsigned int var = get<1>(it.first);
          if (var != prev_var || lag != prev_lag)
            {
              prev_var = var;
              prev_lag = lag;
              count_col_endo++;
            }
        }
      map<tuple<int, int, int>, expr_t> tmp_block_exo_derivative;
      for (const auto &it : derivative_exo[block])
        tmp_block_exo_derivative[{ get<0>(it.first), get<2>(it.first), get<1>(it.first) }] = it.second;
      prev_var = 999999999;
      prev_lag = -9999999;
      count_col_exo = 0;
      for (const auto &it : tmp_block_exo_derivative)
        {
          int lag = get<0>(it.first);
          unsigned int var = get<1>(it.first);
          if (var != prev_var || lag != prev_lag)
            {
              prev_var = var;
              prev_lag = lag;
              count_col_exo++;
            }
        }
      map<tuple<int, int, int>, expr_t> tmp_block_exo_det_derivative;
      for (const auto &it : derivative_exo_det[block])
        tmp_block_exo_det_derivative[{ get<0>(it.first), get<2>(it.first), get<1>(it.first) }] = it.second;
      prev_var = 999999999;
      prev_lag = -9999999;
      count_col_exo_det = 0;
      for (const auto &it : tmp_block_exo_derivative)
        {
          int lag = get<0>(it.first);
          unsigned int var = get<1>(it.first);
          if (var != prev_var || lag != prev_lag)
            {
              prev_var = var;
              prev_lag = lag;
              count_col_exo_det++;
            }
        }
      map<tuple<int, int, int>, expr_t> tmp_block_other_endo_derivative;
      for (const auto &it : derivative_other_endo[block])
        tmp_block_other_endo_derivative[{ get<0>(it.first), get<2>(it.first), get<1>(it.first) }] = it.second;
      prev_var = 999999999;
      prev_lag = -9999999;
      count_col_other_endo = 0;
      for (const auto &it : tmp_block_other_endo_derivative)
        {
          int lag = get<0>(it.first);
          unsigned int var = get<1>(it.first);
          if (var != prev_var || lag != prev_lag)
            {
              prev_var = var;
              prev_lag = lag;
              count_col_other_endo++;
            }
        }

      tmp1_output.str("");
      tmp1_output << packageDir(basename + ".block") << "/dynamic_" << block+1 << ".m";
      output.open(tmp1_output.str(), ios::out | ios::binary);
      output << "%" << endl
             << "% " << tmp1_output.str() << " : Computes dynamic model for Dynare" << endl
             << "%" << endl
             << "% Warning : this file is generated automatically by Dynare" << endl
             << "%           from model file (.mod)" << endl << endl
             << "%/" << endl;
      if (simulation_type == EVALUATE_BACKWARD || simulation_type == EVALUATE_FORWARD)
        {
          output << "function [y, g1, g2, g3, varargout] = dynamic_" << block+1 << "(y, x, params, steady_state, jacobian_eval, y_kmin, periods)" << endl;
        }
      else if (simulation_type == SOLVE_FORWARD_COMPLETE || simulation_type == SOLVE_BACKWARD_COMPLETE)
        output << "function [residual, y, g1, g2, g3, varargout] = dynamic_" << block+1 << "(y, x, params, steady_state, it_, jacobian_eval)" << endl;
      else if (simulation_type == SOLVE_BACKWARD_SIMPLE || simulation_type == SOLVE_FORWARD_SIMPLE)
        output << "function [residual, y, g1, g2, g3, varargout] = dynamic_" << block+1 << "(y, x, params, steady_state, it_, jacobian_eval)" << endl;
      else
        output << "function [residual, y, g1, g2, g3, b, varargout] = dynamic_" << block+1 << "(y, x, params, steady_state, periods, jacobian_eval, y_kmin, y_size, Periods)" << endl;
      BlockType block_type;
      if (simulation_type == SOLVE_TWO_BOUNDARIES_COMPLETE || simulation_type == SOLVE_TWO_BOUNDARIES_SIMPLE)
        block_type = SIMULTAN;
      else if (simulation_type == SOLVE_FORWARD_COMPLETE || simulation_type == SOLVE_BACKWARD_COMPLETE)
        block_type = SIMULTANS;
      else if ((simulation_type == SOLVE_FORWARD_SIMPLE || simulation_type == SOLVE_BACKWARD_SIMPLE
                || simulation_type == EVALUATE_BACKWARD    || simulation_type == EVALUATE_FORWARD)
               && getBlockFirstEquation(block) < prologue)
        block_type = PROLOGUE;
      else if ((simulation_type == SOLVE_FORWARD_SIMPLE || simulation_type == SOLVE_BACKWARD_SIMPLE
                || simulation_type == EVALUATE_BACKWARD    || simulation_type == EVALUATE_FORWARD)
               && getBlockFirstEquation(block) >= equations.size() - epilogue)
        block_type = EPILOGUE;
      else
        block_type = SIMULTANS;
      output << "  % ////////////////////////////////////////////////////////////////////////" << endl
             << "  % //" << string("                     Block ").substr(int (log10(block + 1))) << block + 1 << " " << BlockType0(block_type)
             << "          //" << endl
             << "  % //                     Simulation type "
             << BlockSim(simulation_type) << "  //" << endl
             << "  % ////////////////////////////////////////////////////////////////////////" << endl;
      //The Temporary terms
      if (simulation_type == EVALUATE_BACKWARD || simulation_type == EVALUATE_FORWARD)
        {
          output << "  if(jacobian_eval)" << endl
                 << "    g1 = spalloc(" << block_mfs  << ", " << count_col_endo << ", " << nze << ");" << endl
                 << "    g1_x=spalloc(" << block_size << ", " << count_col_exo  << ", " << nze_exo << ");" << endl
                 << "    g1_xd=spalloc(" << block_size << ", " << count_col_exo_det  << ", " << nze_exo_det << ");" << endl
                 << "    g1_o=spalloc(" << block_size << ", " << count_col_other_endo << ", " << nze_other_endo << ");" << endl
                 << "  end;" << endl;
        }
      else
        {
          output << "  if(jacobian_eval)" << endl
                 << "    g1 = spalloc(" << block_size << ", " << count_col_endo << ", " << nze << ");" << endl
                 << "    g1_x=spalloc(" << block_size << ", " << count_col_exo  << ", " << nze_exo << ");" << endl
                 << "    g1_xd=spalloc(" << block_size << ", " << count_col_exo_det  << ", " << nze_exo_det << ");" << endl
                 << "    g1_o=spalloc(" << block_size << ", " << count_col_other_endo << ", " << nze_other_endo << ");" << endl
                 << "  else" << endl;
          if (simulation_type == SOLVE_TWO_BOUNDARIES_COMPLETE || simulation_type == SOLVE_TWO_BOUNDARIES_SIMPLE)
            output << "    g1 = spalloc(" << block_mfs << "*Periods, "
                   << block_mfs << "*(Periods+" << max_leadlag_block[block].first+max_leadlag_block[block].second+1 << ")"
                   << ", " << nze << "*Periods);" << endl;
          else
            output << "    g1 = spalloc(" << block_mfs
                   << ", " << block_mfs << ", " << nze << ");" << endl;
          output << "  end;" << endl;
        }

      output << "  g2=0;g3=0;" << endl;
      if (v_temporary_terms_inuse[block].size())
        {
          tmp_output.str("");
          for (int it : v_temporary_terms_inuse[block])
            tmp_output << " T" << it;
          output << "  global" << tmp_output.str() << ";" << endl;
        }
      if (simulation_type == SOLVE_TWO_BOUNDARIES_COMPLETE || simulation_type == SOLVE_TWO_BOUNDARIES_SIMPLE)
        {
          temporary_terms_t tt2;
          for (int i = 0; i < (int) block_size; i++)
            {
              if (v_temporary_terms[block][i].size() && global_temporary_terms)
                {
                  output << "  " << "% //Temporary variables initialization" << endl
                         << "  " << "T_zeros = zeros(y_kmin+periods, 1);" << endl;
                  for (auto it : v_temporary_terms[block][i])
                    {
                      output << "  ";
                      // In the following, "Static" is used to avoid getting the "(it_)" subscripting
                      it->writeOutput(output, ExprNodeOutputType::matlabStaticModelSparse, local_temporary_terms, {});
                      output << " = T_zeros;" << endl;
                    }
                }
            }
        }
      if (simulation_type == SOLVE_BACKWARD_SIMPLE || simulation_type == SOLVE_FORWARD_SIMPLE || simulation_type == SOLVE_BACKWARD_COMPLETE || simulation_type == SOLVE_FORWARD_COMPLETE)
        output << "  residual=zeros(" << block_mfs << ",1);" << endl;
      else if (simulation_type == SOLVE_TWO_BOUNDARIES_COMPLETE || simulation_type == SOLVE_TWO_BOUNDARIES_SIMPLE)
        output << "  residual=zeros(" << block_mfs << ",y_kmin+periods);" << endl;
      if (simulation_type == EVALUATE_BACKWARD)
        output << "  for it_ = (y_kmin+periods):y_kmin+1" << endl;
      if (simulation_type == EVALUATE_FORWARD)
        output << "  for it_ = y_kmin+1:(y_kmin+periods)" << endl;

      if (simulation_type == SOLVE_TWO_BOUNDARIES_COMPLETE || simulation_type == SOLVE_TWO_BOUNDARIES_SIMPLE)
        {
          output << "  b = zeros(periods*y_size,1);" << endl
                 << "  for it_ = y_kmin+1:(periods+y_kmin)" << endl
                 << "    Per_y_=it_*y_size;" << endl
                 << "    Per_J_=(it_-y_kmin-1)*y_size;" << endl
                 << "    Per_K_=(it_-1)*y_size;" << endl;
          sps = "  ";
        }
      else
        if (simulation_type == EVALUATE_BACKWARD || simulation_type == EVALUATE_FORWARD)
          sps = "  ";
        else
          sps = "";
      // The equations
      temporary_terms_idxs_t temporary_terms_idxs;
      for (unsigned int i = 0; i < block_size; i++)
        {
          temporary_terms_t tt2;
          if (v_temporary_terms[block].size())
            {
              output << "  " << "% //Temporary variables" << endl;
              for (auto it : v_temporary_terms[block][i])
                {
                  if (dynamic_cast<AbstractExternalFunctionNode *>(it) != nullptr)
                    it->writeExternalFunctionOutput(output, local_output_type, tt2, temporary_terms_idxs, tef_terms);

                  output << "  " <<  sps;
                  it->writeOutput(output, local_output_type, local_temporary_terms, {}, tef_terms);
                  output << " = ";
                  it->writeOutput(output, local_output_type, tt2, {}, tef_terms);
                  // Insert current node into tt2
                  tt2.insert(it);
                  output << ";" << endl;
                }
            }

          int variable_ID = getBlockVariableID(block, i);
          int equation_ID = getBlockEquationID(block, i);
          EquationType equ_type = getBlockEquationType(block, i);
          string sModel = symbol_table.getName(symbol_table.getID(SymbolType::endogenous, variable_ID));
          eq_node = (BinaryOpNode *) getBlockEquationExpr(block, i);
          lhs = eq_node->get_arg1();
          rhs = eq_node->get_arg2();
          tmp_output.str("");
          lhs->writeOutput(tmp_output, local_output_type, local_temporary_terms, {});
          switch (simulation_type)
            {
            case EVALUATE_BACKWARD:
            case EVALUATE_FORWARD:
            evaluation:     if (simulation_type == SOLVE_TWO_BOUNDARIES_COMPLETE || simulation_type == SOLVE_TWO_BOUNDARIES_SIMPLE)
                output << "    % equation " << getBlockEquationID(block, i)+1 << " variable : " << sModel
                       << " (" << variable_ID+1 << ") " << c_Equation_Type(equ_type) << endl;
              output << "    ";
              if (equ_type == E_EVALUATE)
                {
                  output << tmp_output.str();
                  output << " = ";
                  rhs->writeOutput(output, local_output_type, local_temporary_terms, {});
                }
              else if (equ_type == E_EVALUATE_S)
                {
                  output << "%" << tmp_output.str();
                  output << " = ";
                  if (isBlockEquationRenormalized(block, i))
                    {
                      rhs->writeOutput(output, local_output_type, local_temporary_terms, {});
                      output << endl << "    ";
                      tmp_output.str("");
                      eq_node = (BinaryOpNode *) getBlockEquationRenormalizedExpr(block, i);
                      lhs = eq_node->get_arg1();
                      rhs = eq_node->get_arg2();
                      lhs->writeOutput(output, local_output_type, local_temporary_terms, {});
                      output << " = ";
                      rhs->writeOutput(output, local_output_type, local_temporary_terms, {});
                    }
                }
              else
                {
                  cerr << "Type mismatch for equation " << equation_ID+1  << endl;
                  exit(EXIT_FAILURE);
                }
              output << ";" << endl;
              break;
            case SOLVE_BACKWARD_SIMPLE:
            case SOLVE_FORWARD_SIMPLE:
            case SOLVE_BACKWARD_COMPLETE:
            case SOLVE_FORWARD_COMPLETE:
              if (i < block_recursive)
                goto evaluation;
              feedback_variables.push_back(variable_ID);
              output << "  % equation " << equation_ID+1 << " variable : " << sModel
                     << " (" << variable_ID+1 << ") " << c_Equation_Type(equ_type) << " symb_id=" << symbol_table.getID(SymbolType::endogenous, variable_ID) << endl;
              output << "  " << "residual(" << i+1-block_recursive << ") = (";
              goto end;
            case SOLVE_TWO_BOUNDARIES_COMPLETE:
            case SOLVE_TWO_BOUNDARIES_SIMPLE:
              if (i < block_recursive)
                goto evaluation;
              feedback_variables.push_back(variable_ID);
              output << "    % equation " << equation_ID+1 << " variable : " << sModel
                     << " (" << variable_ID+1 << ") " << c_Equation_Type(equ_type) << " symb_id=" << symbol_table.getID(SymbolType::endogenous, variable_ID) << endl;
              Ufoss << "    b(" << i+1-block_recursive << "+Per_J_) = -residual(" << i+1-block_recursive << ", it_)";
              Uf[equation_ID] += Ufoss.str();
              Ufoss.str("");
              output << "    residual(" << i+1-block_recursive << ", it_) = (";
              goto end;
            default:
            end:
              output << tmp_output.str();
              output << ") - (";
              rhs->writeOutput(output, local_output_type, local_temporary_terms, {});
              output << ");" << endl;
#ifdef CONDITION
              if (simulation_type == SOLVE_TWO_BOUNDARIES_COMPLETE || simulation_type == SOLVE_TWO_BOUNDARIES_SIMPLE)
                output << "  condition(" << i+1 << ")=0;" << endl;
#endif
            }
        }
      // The Jacobian if we have to solve the block
      if (simulation_type == SOLVE_TWO_BOUNDARIES_SIMPLE || simulation_type == SOLVE_TWO_BOUNDARIES_COMPLETE)
        output << "  " << sps << "% Jacobian  " << endl << "    if jacobian_eval" << endl;
      else
        if (simulation_type == SOLVE_BACKWARD_SIMPLE   || simulation_type == SOLVE_FORWARD_SIMPLE
            || simulation_type == SOLVE_BACKWARD_COMPLETE || simulation_type == SOLVE_FORWARD_COMPLETE)
          output << "  % Jacobian  " << endl << "  if jacobian_eval" << endl;
        else
          output << "    % Jacobian  " << endl << "    if jacobian_eval" << endl;
      prev_var = 999999999;
      prev_lag = -9999999;
      count_col = 0;
      for (const auto &it : tmp_block_endo_derivative)
        {
          int lag;
          unsigned int var, eq;
          tie(lag, var, eq) = it.first;
          int eqr = getBlockEquationID(block, eq);
          int varr = getBlockVariableID(block, var);
          if (var != prev_var || lag != prev_lag)
            {
              prev_var = var;
              prev_lag = lag;
              count_col++;
            }

          expr_t id = it.second;

          output << "      g1(" << eq+1 << ", " << count_col << ") = ";
          id->writeOutput(output, local_output_type, local_temporary_terms, {});
          output << "; % variable=" << symbol_table.getName(symbol_table.getID(SymbolType::endogenous, varr))
                 << "(" << lag
                 << ") " << varr+1 << ", " << var+1
                 << ", equation=" << eqr+1 << ", " << eq+1 << endl;
        }
      prev_var = 999999999;
      prev_lag = -9999999;
      count_col = 0;
      for (const auto &it : tmp_block_exo_derivative)
        {
          int lag;
          unsigned int var, eq;
          tie(lag, var, eq) = it.first;
          int eqr = getBlockInitialEquationID(block, eq);
          if (var != prev_var || lag != prev_lag)
            {
              prev_var = var;
              prev_lag = lag;
              count_col++;
            }
          expr_t id = it.second;
          output << "      g1_x(" << eqr+1 << ", " << count_col << ") = ";
          id->writeOutput(output, local_output_type, local_temporary_terms, {});
          output << "; % variable=" << symbol_table.getName(symbol_table.getID(SymbolType::exogenous, var))
                 << "(" << lag
                 << ") " << var+1
                 << ", equation=" << eq+1 << endl;
        }
      prev_var = 999999999;
      prev_lag = -9999999;
      count_col = 0;
      for (const auto &it : tmp_block_exo_det_derivative)
        {
          int lag;
          unsigned int var, eq;
          tie(lag, var, eq) = it.first;
          int eqr = getBlockInitialEquationID(block, eq);
          if (var != prev_var || lag != prev_lag)
            {
              prev_var = var;
              prev_lag = lag;
              count_col++;
            }
          expr_t id = it.second;
          output << "      g1_xd(" << eqr+1 << ", " << count_col << ") = ";
          id->writeOutput(output, local_output_type, local_temporary_terms, {});
          output << "; % variable=" << symbol_table.getName(symbol_table.getID(SymbolType::exogenous, var))
                 << "(" << lag
                 << ") " << var+1
                 << ", equation=" << eq+1 << endl;
        }
      prev_var = 999999999;
      prev_lag = -9999999;
      count_col = 0;
      for (const auto &it : tmp_block_other_endo_derivative)
        {
          int lag;
          unsigned int var, eq;
          tie(lag, var, eq) = it.first;
          int eqr = getBlockInitialEquationID(block, eq);
          if (var != prev_var || lag != prev_lag)
            {
              prev_var = var;
              prev_lag = lag;
              count_col++;
            }
          expr_t id = it.second;

          output << "      g1_o(" << eqr+1 << ", " << /*var+1+(lag+block_max_lag)*block_size*/ count_col << ") = ";
          id->writeOutput(output, local_output_type, local_temporary_terms, {});
          output << "; % variable=" << symbol_table.getName(symbol_table.getID(SymbolType::endogenous, var))
                 << "(" << lag
                 << ") " << var+1
                 << ", equation=" << eq+1 << endl;
        }
      output << "      varargout{1}=g1_x;" << endl
             << "      varargout{2}=g1_xd;" << endl
             << "      varargout{3}=g1_o;" << endl;

      switch (simulation_type)
        {
        case EVALUATE_FORWARD:
        case EVALUATE_BACKWARD:
          output << "    end;" << endl
                 << "  end;" << endl;
          break;
        case SOLVE_BACKWARD_SIMPLE:
        case SOLVE_FORWARD_SIMPLE:
        case SOLVE_BACKWARD_COMPLETE:
        case SOLVE_FORWARD_COMPLETE:
          output << "  else" << endl;
          for (const auto &it : blocks_derivatives[block])
            {
              unsigned int eq, var;
              expr_t id;
              int lag;
              tie(eq, var, lag, id) = it;
              unsigned int eqr = getBlockEquationID(block, eq);
              unsigned int varr = getBlockVariableID(block, var);
              if (lag == 0)
                {
                  output << "    g1(" << eq+1 << ", " << var+1-block_recursive << ") = ";
                  id->writeOutput(output, local_output_type, local_temporary_terms, {});
                  output << "; % variable=" << symbol_table.getName(symbol_table.getID(SymbolType::endogenous, varr))
                         << "(" << lag
                         << ") " << varr+1
                         << ", equation=" << eqr+1 << endl;
                }

            }
          output << "  end;" << endl;
          break;
        case SOLVE_TWO_BOUNDARIES_SIMPLE:
        case SOLVE_TWO_BOUNDARIES_COMPLETE:
          output << "    else" << endl;
          for (const auto &it : blocks_derivatives[block])
            {
              unsigned int eq, var;
              int lag;
              expr_t id;
              tie(eq, var, lag, id) = it;
              unsigned int eqr = getBlockEquationID(block, eq);
              unsigned int varr = getBlockVariableID(block, var);
              ostringstream tmp_output;
              if (eq >= block_recursive && var >= block_recursive)
                {
                  if (lag == 0)
                    Ufoss << "+g1(" << eq+1-block_recursive
                          << "+Per_J_, " << var+1-block_recursive
                          << "+Per_K_)*y(it_, " << varr+1 << ")";
                  else if (lag == 1)
                    Ufoss << "+g1(" << eq+1-block_recursive
                          << "+Per_J_, " << var+1-block_recursive
                          << "+Per_y_)*y(it_+1, " << varr+1 << ")";
                  else if (lag > 0)
                    Ufoss << "+g1(" << eq+1-block_recursive
                          << "+Per_J_, " << var+1-block_recursive
                          << "+y_size*(it_+" << lag-1 << "))*y(it_+" << lag << ", " << varr+1 << ")";
                  else
                    Ufoss << "+g1(" << eq+1-block_recursive
                          << "+Per_J_, " << var+1-block_recursive
                          << "+y_size*(it_" << lag-1 << "))*y(it_" << lag << ", " << varr+1 << ")";
                  Uf[eqr] += Ufoss.str();
                  Ufoss.str("");

                  if (lag == 0)
                    tmp_output << "     g1(" << eq+1-block_recursive << "+Per_J_, "
                               << var+1-block_recursive << "+Per_K_) = ";
                  else if (lag == 1)
                    tmp_output << "     g1(" << eq+1-block_recursive << "+Per_J_, "
                               << var+1-block_recursive << "+Per_y_) = ";
                  else if (lag > 0)
                    tmp_output << "     g1(" << eq+1-block_recursive << "+Per_J_, "
                               << var+1-block_recursive << "+y_size*(it_+" << lag-1 << ")) = ";
                  else if (lag < 0)
                    tmp_output << "     g1(" << eq+1-block_recursive << "+Per_J_, "
                               << var+1-block_recursive << "+y_size*(it_" << lag-1 << ")) = ";
                  output << " " << tmp_output.str();
                  id->writeOutput(output, local_output_type, local_temporary_terms, {});
                  output << ";";
                  output << " %2 variable=" << symbol_table.getName(symbol_table.getID(SymbolType::endogenous, varr))
                         << "(" << lag << ") " << varr+1
                         << ", equation=" << eqr+1 << " (" << eq+1 << ")" << endl;
                }

#ifdef CONDITION
              output << "  if (fabs(condition[" << eqr << "])<fabs(u[" << u << "+Per_u_]))" << endl
                     << "    condition(" << eqr << ")=u(" << u << "+Per_u_);" << endl;
#endif
            }
          for (unsigned int i = 0; i < block_size; i++)
            {
              if (i >= block_recursive)
                output << "  " << Uf[getBlockEquationID(block, i)] << ";" << endl;
#ifdef CONDITION
              output << "  if (fabs(condition(" << i+1 << "))<fabs(u(" << i << "+Per_u_)))" << endl
                     << "    condition(" << i+1 << ")=u(" << i+1 << "+Per_u_);" << endl;
#endif
            }
#ifdef CONDITION
          for (m = 0; m <= ModelBlock->Block_List[block].Max_Lead+ModelBlock->Block_List[block].Max_Lag; m++)
            {
              k = m-ModelBlock->Block_List[block].Max_Lag;
              for (i = 0; i < ModelBlock->Block_List[block].IM_lead_lag[m].size; i++)
                {
                  unsigned int eq = ModelBlock->Block_List[block].IM_lead_lag[m].Equ_Index[i];
                  unsigned int var = ModelBlock->Block_List[block].IM_lead_lag[m].Var_Index[i];
                  unsigned int u = ModelBlock->Block_List[block].IM_lead_lag[m].u[i];
                  unsigned int eqr = ModelBlock->Block_List[block].IM_lead_lag[m].Equ[i];
                  output << "  u(" << u+1 << "+Per_u_) = u(" << u+1 << "+Per_u_) / condition(" << eqr+1 << ");" << endl;
                }
            }
          for (i = 0; i < ModelBlock->Block_List[block].Size; i++)
            output << "  u(" << i+1 << "+Per_u_) = u(" << i+1 << "+Per_u_) / condition(" << i+1 << ");" << endl;
#endif
          output << "    end;" << endl
                 << "  end;" << endl;
          break;
        default:
          break;
        }
      output << "end" << endl;
      output.close();
    }
}

void
DynamicModel::writeModelEquationsCode(const string &basename, const map_idx_t &map_idx) const
{

  ostringstream tmp_output;
  ofstream code_file;
  unsigned int instruction_number = 0;
  bool file_open = false;

  boost::filesystem::create_directories(basename + "/model/bytecode");

  string main_name = basename + "/model/bytecode/dynamic.cod";
  code_file.open(main_name, ios::out | ios::binary | ios::ate);
  if (!code_file.is_open())
    {
      cerr << "Error : Can't open file \"" << main_name << "\" for writing" << endl;
      exit(EXIT_FAILURE);
    }

  int count_u;
  int u_count_int = 0;
  BlockSimulationType simulation_type;
  if ((max_endo_lag > 0) && (max_endo_lead > 0))
    simulation_type = SOLVE_TWO_BOUNDARIES_COMPLETE;
  else if ((max_endo_lag >= 0) && (max_endo_lead == 0))
    simulation_type = SOLVE_FORWARD_COMPLETE;
  else
    simulation_type = SOLVE_BACKWARD_COMPLETE;

  Write_Inf_To_Bin_File(basename + "/model/bytecode/dynamic.bin", u_count_int, file_open, simulation_type == SOLVE_TWO_BOUNDARIES_COMPLETE, symbol_table.endo_nbr());
  file_open = true;

  //Temporary variables declaration
  FDIMT_ fdimt(temporary_terms.size());
  fdimt.write(code_file, instruction_number);

  vector<unsigned int> exo, exo_det, other_endo;

  for (int i = 0; i < symbol_table.exo_det_nbr(); i++)
    exo_det.push_back(i);
  for (int i = 0; i < symbol_table.exo_nbr(); i++)
    exo.push_back(i);

  map<tuple<int, int, int>, expr_t> first_derivatives_reordered_endo;
  map<tuple<int, SymbolType, int, int>, expr_t> first_derivatives_reordered_exo;
  for (const auto & first_derivative : derivatives[1])
    {
      int deriv_id = first_derivative.first[1];
      unsigned int eq = first_derivative.first[0];
      int symb = getSymbIDByDerivID(deriv_id);
      unsigned int var = symbol_table.getTypeSpecificID(symb);
      int lag = getLagByDerivID(deriv_id);
      if (getTypeByDerivID(deriv_id) == SymbolType::endogenous)
        first_derivatives_reordered_endo[{ lag, var, eq }] = first_derivative.second;
      else if (getTypeByDerivID(deriv_id) == SymbolType::exogenous || getTypeByDerivID(deriv_id) == SymbolType::exogenousDet)
        first_derivatives_reordered_exo[{ lag, getTypeByDerivID(deriv_id), var, eq }] = first_derivative.second;
    }
  int prev_var = -1;
  int prev_lag = -999999999;
  int count_col_endo = 0;
  for (const auto &it : first_derivatives_reordered_endo)
    {
      int var, lag;
      tie(lag, var, ignore) = it.first;
      if (prev_var != var || prev_lag != lag)
        {
          prev_var = var;
          prev_lag = lag;
          count_col_endo++;
        }
    }
  prev_var = -1;
  prev_lag = -999999999;
  SymbolType prev_type{SymbolType::unusedEndogenous}; // Any non-exogenous type would do here
  int count_col_exo = 0;
  int count_col_det_exo = 0;

  for (const auto & it : first_derivatives_reordered_exo)
    {
      int var, lag;
      SymbolType type;
      tie(lag, type, var, ignore) = it.first;
      if (prev_var != var || prev_lag != lag || prev_type != type)
        {
          prev_var = var;
          prev_lag = lag;
          prev_type = type;
          if (type == SymbolType::exogenous)
            count_col_exo++;
          else if (type == SymbolType::exogenousDet)
            count_col_det_exo++;
        }
    }

  FBEGINBLOCK_ fbeginblock(symbol_table.endo_nbr(),
                           simulation_type,
                           0,
                           symbol_table.endo_nbr(),
                           variable_reordered,
                           equation_reordered,
                           false,
                           symbol_table.endo_nbr(),
                           max_endo_lag,
                           max_endo_lead,
                           u_count_int,
                           count_col_endo,
                           symbol_table.exo_det_nbr(),
                           count_col_det_exo,
                           symbol_table.exo_nbr(),
                           count_col_exo,
                           0,
                           0,
                           exo_det,
                           exo,
                           other_endo
                           );
  fbeginblock.write(code_file, instruction_number);

  compileTemporaryTerms(code_file, instruction_number, temporary_terms, map_idx, true, false);

  compileModelEquations(code_file, instruction_number, temporary_terms, map_idx, true, false);

  FENDEQU_ fendequ;
  fendequ.write(code_file, instruction_number);

  // Get the current code_file position and jump if eval = true
  streampos pos1 = code_file.tellp();
  FJMPIFEVAL_ fjmp_if_eval(0);
  fjmp_if_eval.write(code_file, instruction_number);
  int prev_instruction_number = instruction_number;

  vector<vector<tuple<int, int, int>>> my_derivatives(symbol_table.endo_nbr());;
  count_u = symbol_table.endo_nbr();
  for (const auto & first_derivative : derivatives[1])
    {
      int deriv_id = first_derivative.first[1];
      if (getTypeByDerivID(deriv_id) == SymbolType::endogenous)
        {
          expr_t d1 = first_derivative.second;
          unsigned int eq = first_derivative.first[0];
          int symb = getSymbIDByDerivID(deriv_id);
          unsigned int var = symbol_table.getTypeSpecificID(symb);
          int lag = getLagByDerivID(deriv_id);
          FNUMEXPR_ fnumexpr(FirstEndoDerivative, eq, var, lag);
          fnumexpr.write(code_file, instruction_number);
          if (!my_derivatives[eq].size())
            my_derivatives[eq].clear();
          my_derivatives[eq].emplace_back(var, lag, count_u);
          d1->compile(code_file, instruction_number, false, temporary_terms, map_idx, true, false);

          FSTPU_ fstpu(count_u);
          fstpu.write(code_file, instruction_number);
          count_u++;
        }
    }
  for (int i = 0; i < symbol_table.endo_nbr(); i++)
    {
      FLDR_ fldr(i);
      fldr.write(code_file, instruction_number);
      if (my_derivatives[i].size())
        {
          for (auto it = my_derivatives[i].begin(); it != my_derivatives[i].end(); it++)
            {
              FLDU_ fldu(get<2>(*it));
              fldu.write(code_file, instruction_number);
              FLDV_ fldv{static_cast<int>(SymbolType::endogenous), static_cast<unsigned int>(get<0>(*it)), get<1>(*it)};
              fldv.write(code_file, instruction_number);
              FBINARY_ fbinary{static_cast<int>(BinaryOpcode::times)};
              fbinary.write(code_file, instruction_number);
              if (it != my_derivatives[i].begin())
                {
                  FBINARY_ fbinary{static_cast<int>(BinaryOpcode::plus)};
                  fbinary.write(code_file, instruction_number);
                }
            }
          FBINARY_ fbinary{static_cast<int>(BinaryOpcode::minus)};
          fbinary.write(code_file, instruction_number);
        }
      FSTPU_ fstpu(i);
      fstpu.write(code_file, instruction_number);
    }

  // Get the current code_file position and jump = true
  streampos pos2 = code_file.tellp();
  FJMP_ fjmp(0);
  fjmp.write(code_file, instruction_number);
  // Set code_file position to previous JMPIFEVAL_ and set the number of instructions to jump
  streampos pos3 = code_file.tellp();
  code_file.seekp(pos1);
  FJMPIFEVAL_ fjmp_if_eval1(instruction_number - prev_instruction_number);
  fjmp_if_eval1.write(code_file, instruction_number);
  code_file.seekp(pos3);
  prev_instruction_number = instruction_number;

  // The Jacobian
  prev_var = -1;
  prev_lag = -999999999;
  count_col_endo = 0;
  for (const auto &it : first_derivatives_reordered_endo)
    {
      unsigned int eq;
      int var, lag;
      tie(lag, var, eq) = it.first;
      expr_t d1 = it.second;
      FNUMEXPR_ fnumexpr(FirstEndoDerivative, eq, var, lag);
      fnumexpr.write(code_file, instruction_number);
      if (prev_var != var || prev_lag != lag)
        {
          prev_var = var;
          prev_lag = lag;
          count_col_endo++;
        }
      d1->compile(code_file, instruction_number, false, temporary_terms, map_idx, true, false);
      FSTPG3_ fstpg3(eq, var, lag, count_col_endo-1);
      fstpg3.write(code_file, instruction_number);
    }
  prev_var = -1;
  prev_lag = -999999999;
  count_col_exo = 0;
  for (const auto & it : first_derivatives_reordered_exo)
    {
      int lag, var, eq;
      tie(lag, ignore, var, eq) = it.first;
      expr_t d1 = it.second;
      FNUMEXPR_ fnumexpr(FirstExoDerivative, eq, var, lag);
      fnumexpr.write(code_file, instruction_number);
      if (prev_var != var || prev_lag != lag)
        {
          prev_var = var;
          prev_lag = lag;
          count_col_exo++;
        }
      d1->compile(code_file, instruction_number, false, temporary_terms, map_idx, true, false);
      FSTPG3_ fstpg3(eq, var, lag, count_col_exo-1);
      fstpg3.write(code_file, instruction_number);
    }
  // Set codefile position to previous JMP_ and set the number of instructions to jump
  pos1 = code_file.tellp();
  code_file.seekp(pos2);
  FJMP_ fjmp1(instruction_number - prev_instruction_number);
  fjmp1.write(code_file, instruction_number);
  code_file.seekp(pos1);

  FENDBLOCK_ fendblock;
  fendblock.write(code_file, instruction_number);
  FEND_ fend;
  fend.write(code_file, instruction_number);
  code_file.close();
}

void
DynamicModel::writeModelEquationsCode_Block(const string &basename, const map_idx_t &map_idx, const bool linear_decomposition) const
{
  struct Uff_l
  {
    int u, var, lag;
    Uff_l *pNext;
  };

  struct Uff
  {
    Uff_l *Ufl, *Ufl_First;
  };

  int i, v;
  string tmp_s;
  ostringstream tmp_output;
  ofstream code_file;
  unsigned int instruction_number = 0;
  expr_t lhs = nullptr, rhs = nullptr;
  BinaryOpNode *eq_node;
  Uff Uf[symbol_table.endo_nbr()];
  map<expr_t, int> reference_count;
  deriv_node_temp_terms_t tef_terms;
  vector<int> feedback_variables;
  bool file_open = false;
  string main_name;
  boost::filesystem::create_directories(basename + "/model/bytecode");
  if (linear_decomposition)
    main_name = basename + "/model/bytecode/non_linear.cod";
  else
    main_name = basename + "/model/bytecode/dynamic.cod";
  code_file.open(main_name, ios::out | ios::binary | ios::ate);
  if (!code_file.is_open())
    {
      cerr << "Error : Can't open file \"" << main_name << "\" for writing" << endl;
      exit(EXIT_FAILURE);
    }
  //Temporary variables declaration

  FDIMT_ fdimt(temporary_terms.size());
  fdimt.write(code_file, instruction_number);

  for (unsigned int block = 0; block < getNbBlocks(); block++)
    {
      feedback_variables.clear();
      if (block > 0)
        {
          FENDBLOCK_ fendblock;
          fendblock.write(code_file, instruction_number);
        }
      int count_u;
      int u_count_int = 0;
      BlockSimulationType simulation_type = getBlockSimulationType(block);
      unsigned int block_size = getBlockSize(block);
      unsigned int block_mfs = getBlockMfs(block);
      unsigned int block_recursive = block_size - block_mfs;
      int block_max_lag = max_leadlag_block[block].first;
      int block_max_lead = max_leadlag_block[block].second;

      if (simulation_type == SOLVE_TWO_BOUNDARIES_SIMPLE || simulation_type == SOLVE_TWO_BOUNDARIES_COMPLETE
          || simulation_type == SOLVE_BACKWARD_COMPLETE || simulation_type == SOLVE_FORWARD_COMPLETE)
        {
          Write_Inf_To_Bin_File_Block(basename, block, u_count_int, file_open,
                                      simulation_type == SOLVE_TWO_BOUNDARIES_COMPLETE || simulation_type == SOLVE_TWO_BOUNDARIES_SIMPLE, linear_decomposition);
          file_open = true;
        }
      map<tuple<int, int, int>, expr_t> tmp_block_endo_derivative;
      for (const auto &it : blocks_derivatives[block])
        tmp_block_endo_derivative[{ get<2>(it), get<1>(it), get<0>(it) }] = get<3>(it);
      map<tuple<int, int, int>, expr_t> tmp_exo_derivative;
      for (const auto &it : derivative_exo[block])
        tmp_exo_derivative[{ get<0>(it.first), get<2>(it.first), get<1>(it.first) }] = it.second;
      map<tuple<int, int, int>, expr_t> tmp_exo_det_derivative;
      for (const auto &it : derivative_exo_det[block])
        tmp_exo_det_derivative[{ get<0>(it.first), get<2>(it.first), get<1>(it.first) }] = it.second;
      map<tuple<int, int, int>, expr_t> tmp_other_endo_derivative;
      for (const auto &it : derivative_other_endo[block])
        tmp_other_endo_derivative[{ get<0>(it.first), get<2>(it.first), get<1>(it.first) }] = it.second;
      int prev_var = -1;
      int prev_lag = -999999999;
      int count_col_endo = 0;
      for (const auto &it : tmp_block_endo_derivative)
        {
          int lag, var;
          tie(lag, var, ignore) = it.first;
          if (prev_var != var || prev_lag != lag)
            {
              prev_var = var;
              prev_lag = lag;
              count_col_endo++;
            }
        }
      unsigned int count_col_det_exo = 0;
      vector<unsigned int> exo_det;
      for (const auto & it : exo_det_block[block])
        for (const auto &it1 : it.second)
          {
            count_col_det_exo++;
            if (find(exo_det.begin(), exo_det.end(), it1) == exo_det.end())
              exo_det.push_back(it1);
          }

      unsigned int count_col_exo = 0;
      vector<unsigned int> exo;
      for (const auto & it : exo_block[block])
        for (const auto &it1 : it.second)
          {
            count_col_exo++;
            if (find(exo.begin(), exo.end(), it1) == exo.end())
              exo.push_back(it1);
          }

      vector<unsigned int> other_endo;
      unsigned int count_col_other_endo = 0;
      for (const auto & it : other_endo_block[block])
        for (const auto &it1 : it.second)
          {
            count_col_other_endo++;
            if (find(other_endo.begin(), other_endo.end(), it1) == other_endo.end())
              other_endo.push_back(it1);
          }

      FBEGINBLOCK_ fbeginblock(block_mfs,
                               simulation_type,
                               getBlockFirstEquation(block),
                               block_size,
                               variable_reordered,
                               equation_reordered,
                               blocks_linear[block],
                               symbol_table.endo_nbr(),
                               block_max_lag,
                               block_max_lead,
                               u_count_int,
                               count_col_endo,
                               exo_det.size(),
                               count_col_det_exo,
                               exo.size(),
                               getBlockExoColSize(block),
                               other_endo.size(),
                               count_col_other_endo,
                               exo_det,
                               exo,
                               other_endo
                               );
      fbeginblock.write(code_file, instruction_number);

      if (linear_decomposition)
        compileTemporaryTerms(code_file, instruction_number, temporary_terms, map_idx, true, false);

      // The equations
      for (i = 0; i < (int) block_size; i++)
        {
          //The Temporary terms
          temporary_terms_t tt2;
          if (v_temporary_terms[block][i].size() && !linear_decomposition)
            {
              for (auto it : v_temporary_terms[block][i])
                {
                  if (dynamic_cast<AbstractExternalFunctionNode *>(it) != nullptr)
                    it->compileExternalFunctionOutput(code_file, instruction_number, false, tt2, map_idx, true, false, tef_terms);

                  FNUMEXPR_ fnumexpr(TemporaryTerm, (int)(map_idx.find(it->idx)->second));
                  fnumexpr.write(code_file, instruction_number);
                  it->compile(code_file, instruction_number, false, tt2, map_idx, true, false, tef_terms);
                  FSTPT_ fstpt((int)(map_idx.find(it->idx)->second));
                  fstpt.write(code_file, instruction_number);
                  // Insert current node into tt2
                  tt2.insert(it);
#ifdef DEBUGC
                  cout << "FSTPT " << v << endl;
                  instruction_number++;
                  code_file.write(&FOK, sizeof(FOK));
                  code_file.write(reinterpret_cast<char *>(&k), sizeof(k));
                  ki++;
#endif

                }
            }
#ifdef DEBUGC
          for (const auto &it : v_temporary_terms[block][i])
            {
              auto ii = map_idx.find(it->idx);
              cout << "map_idx[" << it->idx <<"]=" << ii->second << endl;
            }
#endif

          int variable_ID, equation_ID;
          EquationType equ_type;

          switch (simulation_type)
            {
            evaluation:
            case EVALUATE_BACKWARD:
            case EVALUATE_FORWARD:
              equ_type = getBlockEquationType(block, i);
              {
                FNUMEXPR_ fnumexpr(ModelEquation, getBlockEquationID(block, i));
                fnumexpr.write(code_file, instruction_number);
              }
              if (equ_type == E_EVALUATE)
                {
                  eq_node = (BinaryOpNode *) getBlockEquationExpr(block, i);
                  lhs = eq_node->get_arg1();
                  rhs = eq_node->get_arg2();
                  rhs->compile(code_file, instruction_number, false, temporary_terms, map_idx, true, false);
                  lhs->compile(code_file, instruction_number, true, temporary_terms, map_idx, true, false);
                }
              else if (equ_type == E_EVALUATE_S)
                {
                  eq_node = (BinaryOpNode *) getBlockEquationRenormalizedExpr(block, i);
                  lhs = eq_node->get_arg1();
                  rhs = eq_node->get_arg2();
                  rhs->compile(code_file, instruction_number, false, temporary_terms, map_idx, true, false);
                  lhs->compile(code_file, instruction_number, true, temporary_terms, map_idx, true, false);
                }
              break;
            case SOLVE_BACKWARD_COMPLETE:
            case SOLVE_FORWARD_COMPLETE:
            case SOLVE_TWO_BOUNDARIES_COMPLETE:
            case SOLVE_TWO_BOUNDARIES_SIMPLE:
              if (i < (int) block_recursive)
                goto evaluation;
              variable_ID = getBlockVariableID(block, i);
              equation_ID = getBlockEquationID(block, i);
              feedback_variables.push_back(variable_ID);
              Uf[equation_ID].Ufl = nullptr;
              goto end;
            default:
            end:
              FNUMEXPR_ fnumexpr(ModelEquation, getBlockEquationID(block, i));
              fnumexpr.write(code_file, instruction_number);
              eq_node = (BinaryOpNode *) getBlockEquationExpr(block, i);
              lhs = eq_node->get_arg1();
              rhs = eq_node->get_arg2();
              lhs->compile(code_file, instruction_number, false, temporary_terms, map_idx, true, false);
              rhs->compile(code_file, instruction_number, false, temporary_terms, map_idx, true, false);

              FBINARY_ fbinary{static_cast<int>(BinaryOpcode::minus)};
              fbinary.write(code_file, instruction_number);
              FSTPR_ fstpr(i - block_recursive);
              fstpr.write(code_file, instruction_number);
            }
        }
      FENDEQU_ fendequ;
      fendequ.write(code_file, instruction_number);

      // Get the current code_file position and jump if eval = true
      streampos pos1 = code_file.tellp();
      FJMPIFEVAL_ fjmp_if_eval(0);
      fjmp_if_eval.write(code_file, instruction_number);
      int prev_instruction_number = instruction_number;
      // The Jacobian if we have to solve the block determinsitic block
      if    (simulation_type != EVALUATE_BACKWARD
             && simulation_type != EVALUATE_FORWARD)
        {
          switch (simulation_type)
            {
            case SOLVE_BACKWARD_SIMPLE:
            case SOLVE_FORWARD_SIMPLE:
              {
                FNUMEXPR_ fnumexpr(FirstEndoDerivative, getBlockEquationID(block, 0), getBlockVariableID(block, 0), 0);
                fnumexpr.write(code_file, instruction_number);
              }
              compileDerivative(code_file, instruction_number, getBlockEquationID(block, 0), getBlockVariableID(block, 0), 0, map_idx);
              {
                FSTPG_ fstpg(0);
                fstpg.write(code_file, instruction_number);
              }
              break;

            case SOLVE_BACKWARD_COMPLETE:
            case SOLVE_FORWARD_COMPLETE:
            case SOLVE_TWO_BOUNDARIES_COMPLETE:
            case SOLVE_TWO_BOUNDARIES_SIMPLE:
              count_u = feedback_variables.size();
              for (const auto &it : blocks_derivatives[block])
                {
                  unsigned int eq, var;
                  int lag;
                  tie(eq, var, lag, ignore) = it;
                  unsigned int eqr = getBlockEquationID(block, eq);
                  unsigned int varr = getBlockVariableID(block, var);
                  if (eq >= block_recursive and var >= block_recursive)
                    {
                      if (lag != 0 && (simulation_type == SOLVE_FORWARD_COMPLETE || simulation_type == SOLVE_BACKWARD_COMPLETE))
                        continue;
                      if (!Uf[eqr].Ufl)
                        {
                          Uf[eqr].Ufl = (Uff_l *) malloc(sizeof(Uff_l));
                          Uf[eqr].Ufl_First = Uf[eqr].Ufl;
                        }
                      else
                        {
                          Uf[eqr].Ufl->pNext = (Uff_l *) malloc(sizeof(Uff_l));
                          Uf[eqr].Ufl = Uf[eqr].Ufl->pNext;
                        }
                      Uf[eqr].Ufl->pNext = nullptr;
                      Uf[eqr].Ufl->u = count_u;
                      Uf[eqr].Ufl->var = varr;
                      Uf[eqr].Ufl->lag = lag;
                      FNUMEXPR_ fnumexpr(FirstEndoDerivative, eqr, varr, lag);
                      fnumexpr.write(code_file, instruction_number);
                      compileChainRuleDerivative(code_file, instruction_number, eqr, varr, lag, map_idx);
                      FSTPU_ fstpu(count_u);
                      fstpu.write(code_file, instruction_number);
                      count_u++;
                    }
                }
              for (i = 0; i < (int) block_size; i++)
                {
                  if (i >= (int) block_recursive)
                    {
                      FLDR_ fldr(i-block_recursive);
                      fldr.write(code_file, instruction_number);

                      FLDZ_ fldz;
                      fldz.write(code_file, instruction_number);

                      v = getBlockEquationID(block, i);
                      for (Uf[v].Ufl = Uf[v].Ufl_First; Uf[v].Ufl; Uf[v].Ufl = Uf[v].Ufl->pNext)
                        {
                          FLDU_ fldu(Uf[v].Ufl->u);
                          fldu.write(code_file, instruction_number);
                          FLDV_ fldv{static_cast<int>(SymbolType::endogenous), static_cast<unsigned int>(Uf[v].Ufl->var), Uf[v].Ufl->lag};
                          fldv.write(code_file, instruction_number);

                          FBINARY_ fbinary{static_cast<int>(BinaryOpcode::times)};
                          fbinary.write(code_file, instruction_number);

                          FCUML_ fcuml;
                          fcuml.write(code_file, instruction_number);
                        }
                      Uf[v].Ufl = Uf[v].Ufl_First;
                      while (Uf[v].Ufl)
                        {
                          Uf[v].Ufl_First = Uf[v].Ufl->pNext;
                          free(Uf[v].Ufl);
                          Uf[v].Ufl = Uf[v].Ufl_First;
                        }
                      FBINARY_ fbinary{static_cast<int>(BinaryOpcode::minus)};
                      fbinary.write(code_file, instruction_number);

                      FSTPU_ fstpu(i - block_recursive);
                      fstpu.write(code_file, instruction_number);
                    }
                }
              break;
            default:
              break;
            }
        }
      // Get the current code_file position and jump = true
      streampos pos2 = code_file.tellp();
      FJMP_ fjmp(0);
      fjmp.write(code_file, instruction_number);
      // Set code_file position to previous JMPIFEVAL_ and set the number of instructions to jump
      streampos pos3 = code_file.tellp();
      code_file.seekp(pos1);
      FJMPIFEVAL_ fjmp_if_eval1(instruction_number - prev_instruction_number);
      fjmp_if_eval1.write(code_file, instruction_number);
      code_file.seekp(pos3);
      prev_instruction_number = instruction_number;
      // The Jacobian if we have to solve the block determinsitic block

      prev_var = -1;
      prev_lag = -999999999;
      count_col_endo = 0;
      for (const auto &it : tmp_block_endo_derivative)
        {
          int lag, var;
          unsigned int eq;
          tie(lag, var, eq) = it.first;
          unsigned int eqr = getBlockEquationID(block, eq);
          unsigned int varr = getBlockVariableID(block, var);
          if (prev_var != var || prev_lag != lag)
            {
              prev_var = var;
              prev_lag = lag;
              count_col_endo++;
            }
          FNUMEXPR_ fnumexpr(FirstEndoDerivative, eqr, varr, lag);
          fnumexpr.write(code_file, instruction_number);
          compileDerivative(code_file, instruction_number, eqr, varr, lag, map_idx);
          FSTPG3_ fstpg3(eq, var, lag, count_col_endo-1);
          fstpg3.write(code_file, instruction_number);
        }
      prev_var = -1;
      prev_lag = -999999999;
      count_col_exo = 0;
      for (const auto &it : tmp_exo_derivative)
        {
          int lag, eq, var;
          tie(lag, var, eq) = it.first;
          int eqr = getBlockInitialEquationID(block, eq);
          int varr = getBlockInitialExogenousID(block, var);
          if (prev_var != var || prev_lag != lag)
            {
              prev_var = var;
              prev_lag = lag;
              count_col_exo++;
            }
          expr_t id = it.second;

          FNUMEXPR_ fnumexpr(FirstExoDerivative, eqr, varr, lag);
          fnumexpr.write(code_file, instruction_number);
          id->compile(code_file, instruction_number, false, temporary_terms, map_idx, true, false);
          FSTPG3_ fstpg3(eq, var, lag, /*var*/ count_col_exo-1);
          fstpg3.write(code_file, instruction_number);
        }
      prev_var = -1;
      prev_lag = -999999999;
      int count_col_exo_det = 0;
      for (const auto &it : tmp_exo_det_derivative)
        {
          int lag, eq, var;
          tie(lag, var, eq) = it.first;
          int eqr = getBlockInitialEquationID(block, eq);
          int varr = getBlockInitialDetExogenousID(block, var);
          if (prev_var != var || prev_lag != lag)
            {
              prev_var = var;
              prev_lag = lag;
              count_col_exo_det++;
            }
          expr_t id = it.second;

          FNUMEXPR_ fnumexpr(FirstExodetDerivative, eqr, varr, lag);
          fnumexpr.write(code_file, instruction_number);
          id->compile(code_file, instruction_number, false, temporary_terms, map_idx, true, false);
          FSTPG3_ fstpg3(eq, var, lag, count_col_exo_det-1);
          fstpg3.write(code_file, instruction_number);
        }
      prev_var = -1;
      prev_lag = -999999999;
      count_col_other_endo = 0;
      for (const auto &it : tmp_other_endo_derivative)
        {
          int lag, eq, var;
          tie(lag, var, eq) = it.first;
          int eqr = getBlockInitialEquationID(block, eq);
          int varr = getBlockInitialOtherEndogenousID(block, var);;
          if (prev_var != var || prev_lag != lag)
            {
              prev_var = var;
              prev_lag = lag;
              count_col_other_endo++;
            }
          expr_t id = it.second;

          FNUMEXPR_ fnumexpr(FirstOtherEndoDerivative, eqr, varr, lag);
          fnumexpr.write(code_file, instruction_number);
          id->compile(code_file, instruction_number, false, temporary_terms, map_idx, true, false);
          FSTPG3_ fstpg3(eq, var, lag, count_col_other_endo-1);
          fstpg3.write(code_file, instruction_number);
        }

      // Set codefile position to previous JMP_ and set the number of instructions to jump
      pos1 = code_file.tellp();
      code_file.seekp(pos2);
      FJMP_ fjmp1(instruction_number - prev_instruction_number);
      fjmp1.write(code_file, instruction_number);
      code_file.seekp(pos1);
    }
  FENDBLOCK_ fendblock;
  fendblock.write(code_file, instruction_number);
  FEND_ fend;
  fend.write(code_file, instruction_number);
  code_file.close();
}

void
DynamicModel::writeDynamicMFile(const string &basename) const
{
  writeDynamicModel(basename, false, false);
}

void
DynamicModel::writeDynamicJuliaFile(const string &basename) const
{
  writeDynamicModel(basename, false, true);
}

void
DynamicModel::writeDynamicCFile(const string &basename, const int order) const
{
  boost::filesystem::create_directories(basename + "/model/src");
  string filename = basename + "/model/src/dynamic.c";
  string filename_mex = basename + "/model/src/dynamic_mex.c";
  ofstream mDynamicModelFile, mDynamicMexFile;

  int ntt = temporary_terms_mlv.size() + temporary_terms_derivatives[0].size() + temporary_terms_derivatives[1].size() + temporary_terms_derivatives[2].size() + temporary_terms_derivatives[3].size();

  mDynamicModelFile.open(filename, ios::out | ios::binary);
  if (!mDynamicModelFile.is_open())
    {
      cerr << "Error: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }
  mDynamicModelFile << "/*" << endl
                    << " * " << filename << " : Computes dynamic model for Dynare" << endl
                    << " *" << endl
                    << " * Warning : this file is generated automatically by Dynare" << endl
                    << " *           from model file (.mod)" << endl
                    << " */" << endl
                    << "#include <math.h>" << endl;

  if (external_functions_table.get_total_number_of_unique_model_block_external_functions())
    // External Matlab function, implies Dynamic function will call mex
    mDynamicModelFile << "#include \"mex.h\"" << endl;
  else
    mDynamicModelFile << "#include <stdlib.h>" << endl;

  mDynamicModelFile << "#define max(a, b) (((a) > (b)) ? (a) : (b))" << endl
                    << "#define min(a, b) (((a) > (b)) ? (b) : (a))" << endl;

  // Write function definition if BinaryOpcode::powerDeriv is used
  writePowerDerivCHeader(mDynamicModelFile);

  mDynamicModelFile << endl;

  // Writing the function body
  writeDynamicModel(mDynamicModelFile, true, false);

  mDynamicModelFile << endl;

  writePowerDeriv(mDynamicModelFile);
  mDynamicModelFile.close();

  mDynamicMexFile.open(filename_mex, ios::out | ios::binary);
  if (!mDynamicMexFile.is_open())
    {
      cerr << "Error: Can't open file " << filename_mex << " for writing" << endl;
      exit(EXIT_FAILURE);
    }

  // Writing the gateway routine
  mDynamicMexFile << "/*" << endl
                  << " * " << filename_mex << " : The gateway routine used to call the Dynamic function "
                  << "located in " << filename << endl
                  << " *" << endl
                  << " * Warning : this file is generated automatically by Dynare" << endl
                  << " *           from model file (.mod)" << endl
                  << endl
                  << " */" << endl
                  << endl
                  << "#include <stdlib.h>" << endl
                  << "#include \"mex.h\"" << endl
                  << endl
                  << "const int ntt = " << ntt << ";" << endl
                  << "void dynamic_resid_tt(const double *y, const double *x, int nb_row_x, const double *params, const double *steady_state, int it_, double *T);" << endl
                  << "void dynamic_resid(const double *y, const double *x, int nb_row_x, const double *params, const double *steady_state, int it_, const double *T, double *residual);" << endl
                  << "void dynamic_g1_tt(const double *y, const double *x, int nb_row_x, const double *params, const double *steady_state, int it_, double *T);" << endl
                  << "void dynamic_g1(const double *y, const double *x, int nb_row_x, const double *params, const double *steady_state, int it_, const double *T, double *g1);" << endl
                  << "void dynamic_g2_tt(const double *y, const double *x, int nb_row_x, const double *params, const double *steady_state, int it_, double *T);" << endl
                  << "void dynamic_g2(const double *y, const double *x, int nb_row_x, const double *params, const double *steady_state, int it_, const double *T, double *v2);" << endl
                  << "void dynamic_g3_tt(const double *y, const double *x, int nb_row_x, const double *params, const double *steady_state, int it_, double *T);" << endl
                  << "void dynamic_g3(const double *y, const double *x, int nb_row_x, const double *params, const double *steady_state, int it_, const double *T, double *v3);" << endl
                  << "void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])" << endl
                  << "{" << endl
                  << "  /* Check that no derivatives of higher order than computed are being requested */" << endl
                  << "  if (nlhs > " << order + 1 << ")" << endl
                  << "    mexErrMsgTxt(\"Derivatives of higher order than computed have been requested\");" << endl
                  << "  /* Create a pointer to the input matrix y. */" << endl
                  << "  double *y = mxGetPr(prhs[0]);" << endl
                  << endl
                  << "  /* Create a pointer to the input matrix x. */" << endl
                  << "  double *x = mxGetPr(prhs[1]);" << endl
                  << endl
                  << "  /* Create a pointer to the input matrix params. */" << endl
                  << "  double *params = mxGetPr(prhs[2]);" << endl
                  << endl
                  << "  /* Create a pointer to the input matrix steady_state. */" << endl
                  << "  double *steady_state = mxGetPr(prhs[3]);" << endl
                  << endl
                  << "  /* Fetch time index */" << endl
                  << "  int it_ = (int) mxGetScalar(prhs[4]) - 1;" << endl
                  << endl
                  << "  /* Gets number of rows of matrix x. */" << endl
                  << "  int nb_row_x = mxGetM(prhs[1]);" << endl
                  << endl
                  << "  double *T = (double *) malloc(sizeof(double)*ntt);"
                  << endl
                  << "  if (nlhs >= 1)" << endl
                  << "  {" << endl
                  << "     /* Set the output pointer to the output matrix residual. */" << endl
                  << "     plhs[0] = mxCreateDoubleMatrix(" << equations.size() << ",1, mxREAL);" << endl
                  << "     double *residual = mxGetPr(plhs[0]);" << endl
                  << "     dynamic_resid_tt(y, x, nb_row_x, params, steady_state, it_, T);" << endl
                  << "     dynamic_resid(y, x, nb_row_x, params, steady_state, it_, T, residual);" << endl
                  << "  }" << endl
                  << endl
                  << "  if (nlhs >= 2)" << endl
                  << "  {" << endl
                  << "     /* Set the output pointer to the output matrix g1. */" << endl
                  << "     plhs[1] = mxCreateDoubleMatrix(" << equations.size() << ", " << dynJacobianColsNbr << ", mxREAL);" << endl
                  << "     double *g1 = mxGetPr(plhs[1]);" << endl
                  << "     dynamic_g1_tt(y, x, nb_row_x, params, steady_state, it_, T);" << endl
                  << "     dynamic_g1(y, x, nb_row_x, params, steady_state, it_, T, g1);" << endl
                  << "  }" << endl
                  << endl
                  << " if (nlhs >= 3)" << endl
                  << "  {" << endl
                  << "     /* Set the output pointer to the output matrix v2. */" << endl
                  << "     plhs[2] = mxCreateDoubleMatrix(" << NNZDerivatives[2] << ", " << 3
                  << ", mxREAL);" << endl
                  << "     double *v2 = mxGetPr(plhs[2]);" << endl
                  << "     dynamic_g2_tt(y, x, nb_row_x, params, steady_state, it_, T);" << endl
                  << "     dynamic_g2(y, x, nb_row_x, params, steady_state, it_, T, v2);" << endl
                  << "  }" << endl
                  << endl
                  << " if (nlhs >= 4)" << endl
                  << "  {" << endl
                  << "     /* Set the output pointer to the output matrix v3. */" << endl
                  << "     plhs[3] = mxCreateDoubleMatrix(" << NNZDerivatives[3] << ", " << 3 << ", mxREAL);" << endl
                  << "     double *v3 = mxGetPr(plhs[3]);" << endl
                  << "     dynamic_g3_tt(y, x, nb_row_x, params, steady_state, it_, T);" << endl
                  << "     dynamic_g3(y, x, nb_row_x, params, steady_state, it_, T, v3);" << endl
                  << "  }" << endl
                  << endl
                  << " free(T);"
                  << "}" << endl;
  mDynamicMexFile.close();
}

string
DynamicModel::reform(const string name1) const
{
  string name = name1;
  int pos = name.find("\\", 0);
  while (pos >= 0)
    {
      if (name.substr(pos + 1, 1) != "\\")
        {
          name = name.insert(pos, "\\");
          pos++;
        }
      pos++;
      pos = name.find("\\", pos);
    }
  return (name);
}

void
DynamicModel::printNonZeroHessianEquations(ostream &output) const
{
  if (nonzero_hessian_eqs.size() !=  1)
    output << "[";
  for (auto it = nonzero_hessian_eqs.begin();
       it != nonzero_hessian_eqs.end(); it++)
    {
      if (it != nonzero_hessian_eqs.begin())
        output << " ";
      output << it->first;
    }
  if (nonzero_hessian_eqs.size() != 1)
    output << "]";
}

void
DynamicModel::setNonZeroHessianEquations(map<int, string> &eqs)
{
  for (const auto &it : derivatives[2])
    {
      int eq = it.first[0];
      if (nonzero_hessian_eqs.find(eq) == nonzero_hessian_eqs.end())
        {
          nonzero_hessian_eqs[eq] = "";
          for (auto & equation_tag : equation_tags)
            if (equation_tag.first == eq)
              if (equation_tag.second.first == "name")
                {
                  nonzero_hessian_eqs[eq] = equation_tag.second.second;
                  break;
                }
        }
    }
  eqs = nonzero_hessian_eqs;
}

void
DynamicModel::Write_Inf_To_Bin_File_Block(const string &basename, const int &num,
                                          int &u_count_int, bool &file_open, bool is_two_boundaries, const bool linear_decomposition) const
{
  int j;
  std::ofstream SaveCode;
  string filename;

  if (!linear_decomposition)
    filename = basename + "/model/bytecode/dynamic.bin";
  else
    filename = basename + "/model/bytecode/non_linear.bin";

  if (file_open)
    SaveCode.open(filename, ios::out | ios::in | ios::binary | ios::ate);
  else
    SaveCode.open(filename, ios::out | ios::binary);
  if (!SaveCode.is_open())
    {
      cerr << "Error : Can't open file \"" << filename << "\" for writing" << endl;
      exit(EXIT_FAILURE);
    }
  u_count_int = 0;
  unsigned int block_size = getBlockSize(num);
  unsigned int block_mfs = getBlockMfs(num);
  unsigned int block_recursive = block_size - block_mfs;
  for (const auto &it : blocks_derivatives[num])
    {
      unsigned int eq, var;
      int lag;
      tie(eq, var, lag, ignore) = it;
      if (lag != 0 && !is_two_boundaries)
        continue;
      if (eq >= block_recursive && var >= block_recursive)
        {
          int v = eq - block_recursive;
          SaveCode.write(reinterpret_cast<char *>(&v), sizeof(v));
          int varr = var - block_recursive + lag * block_mfs;
          SaveCode.write(reinterpret_cast<char *>(&varr), sizeof(varr));
          SaveCode.write(reinterpret_cast<char *>(&lag), sizeof(lag));
          int u = u_count_int + block_mfs;
          SaveCode.write(reinterpret_cast<char *>(&u), sizeof(u));
          u_count_int++;
        }
    }

  if (is_two_boundaries)
    u_count_int += block_mfs;
  for (j = block_recursive; j < (int) block_size; j++)
    {
      unsigned int varr = getBlockVariableID(num, j);
      SaveCode.write(reinterpret_cast<char *>(&varr), sizeof(varr));
    }
  for (j = block_recursive; j < (int) block_size; j++)
    {
      unsigned int eqr = getBlockEquationID(num, j);
      SaveCode.write(reinterpret_cast<char *>(&eqr), sizeof(eqr));
    }
  SaveCode.close();
}

void
DynamicModel::writeSparseDynamicMFile(const string &basename) const
{
  string sp;
  ofstream mDynamicModelFile;
  ostringstream tmp, tmp1, tmp_eq;
  bool OK;
  string filename = packageDir(basename) + "/dynamic.m";
  mDynamicModelFile.open(filename, ios::out | ios::binary);
  if (!mDynamicModelFile.is_open())
    {
      cerr << "Error: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }
  mDynamicModelFile << "%" << endl
                    << "% " << filename << " : Computes dynamic model for Dynare" << endl
                    << "%" << endl
                    << "% Warning : this file is generated automatically by Dynare" << endl
                    << "%           from model file (.mod)" << endl << endl
                    << "%/" << endl;

  int Nb_SGE = 0;
  bool open_par = false;

  mDynamicModelFile << "function [varargout] = dynamic(options_, M_, oo_, varargin)" << endl
                    << "  g2=[];g3=[];" << endl;
  //Temporary variables declaration
  OK = true;
  ostringstream tmp_output;
  for (auto temporary_term : temporary_terms)
    {
      if (OK)
        OK = false;
      else
        tmp_output << " ";
      // In the following, "Static" is used to avoid getting the "(it_)" subscripting
      temporary_term->writeOutput(tmp_output, ExprNodeOutputType::matlabStaticModelSparse, temporary_terms, {});
    }
  if (tmp_output.str().length() > 0)
    mDynamicModelFile << "  global " << tmp_output.str() << ";" << endl;

  mDynamicModelFile << "  T_init=zeros(1,options_.periods+M_.maximum_lag+M_.maximum_lead);" << endl;
  tmp_output.str("");
  for (auto temporary_term : temporary_terms)
    {
      tmp_output << "  ";
      // In the following, "Static" is used to avoid getting the "(it_)" subscripting
      temporary_term->writeOutput(tmp_output, ExprNodeOutputType::matlabStaticModelSparse, temporary_terms, {});
      tmp_output << "=T_init;" << endl;
    }
  if (tmp_output.str().length() > 0)
    mDynamicModelFile << tmp_output.str();

  mDynamicModelFile << "  y_kmin=M_.maximum_lag;" << endl
                    << "  y_kmax=M_.maximum_lead;" << endl
                    << "  y_size=M_.endo_nbr;" << endl
                    << "  if(length(varargin)>0)" << endl
                    << "    %it is a simple evaluation of the dynamic model for time _it" << endl
                    << "    y=varargin{1};" << endl
                    << "    x=varargin{2};" << endl
                    << "    params=varargin{3};" << endl
                    << "    steady_state=varargin{4};" << endl
                    << "    it_=varargin{5};" << endl
                    << "    dr=varargin{6};" << endl
                    << "    Per_u_=0;" << endl
                    << "    Per_y_=it_*y_size;" << endl
                    << "    ys=y(it_,:);" << endl;
  tmp.str("");
  tmp_eq.str("");
  unsigned int nb_blocks = getNbBlocks();
  unsigned int block = 0;
  for (int count_call = 1; block < nb_blocks; block++, count_call++)
    {
      unsigned int block_size = getBlockSize(block);
      unsigned int block_mfs = getBlockMfs(block);
      unsigned int block_recursive = block_size - block_mfs;
      BlockSimulationType simulation_type = getBlockSimulationType(block);

      if (simulation_type == EVALUATE_FORWARD || simulation_type == EVALUATE_BACKWARD)
        {
          for (unsigned int ik = 0; ik < block_size; ik++)
            {
              tmp << " " << getBlockVariableID(block, ik)+1;
              tmp_eq << " " << getBlockEquationID(block, ik)+1;
            }
        }
      else
        {
          for (unsigned int ik = block_recursive; ik < block_size; ik++)
            {
              tmp << " " << getBlockVariableID(block, ik)+1;
              tmp_eq << " " << getBlockEquationID(block, ik)+1;
            }
        }
      mDynamicModelFile << "    y_index_eq=[" << tmp_eq.str() << "];" << endl
                        << "    y_index=[" << tmp.str() << "];" << endl;

      switch (simulation_type)
        {
        case EVALUATE_FORWARD:
        case EVALUATE_BACKWARD:
          mDynamicModelFile << "    [y, dr(" << count_call << ").g1, dr(" << count_call << ").g2, dr(" << count_call << ").g3, dr(" << count_call << ").g1_x, dr(" << count_call << ").g1_xd, dr(" << count_call << ").g1_o]=" << basename << ".block.dynamic_" << block + 1 << "(y, x, params, steady_state, 1, it_-1, 1);" << endl
                            << "    residual(y_index_eq)=ys(y_index)-y(it_, y_index);" << endl;
          break;
        case SOLVE_FORWARD_SIMPLE:
        case SOLVE_BACKWARD_SIMPLE:
          mDynamicModelFile << "    [r, y, dr(" << count_call << ").g1, dr(" << count_call << ").g2, dr(" << count_call << ").g3, dr(" << count_call << ").g1_x, dr(" << count_call << ").g1_xd, dr(" << count_call << ").g1_o]=" << basename << ".block.dynamic_" << block + 1 << "(y, x, params, steady_state, it_, 1);" << endl
                            << "    residual(y_index_eq)=r;" << endl;
          break;
        case SOLVE_FORWARD_COMPLETE:
        case SOLVE_BACKWARD_COMPLETE:
          mDynamicModelFile << "    [r, y, dr(" << count_call << ").g1, dr(" << count_call << ").g2, dr(" << count_call << ").g3, dr(" << count_call << ").g1_x, dr(" << count_call << ").g1_xd, dr(" << count_call << ").g1_o]=" << basename << ".block.dynamic_" << block + 1 << "(y, x, params, steady_state, it_, 1);" << endl
                            << "    residual(y_index_eq)=r;" << endl;
          break;
        case SOLVE_TWO_BOUNDARIES_COMPLETE:
        case SOLVE_TWO_BOUNDARIES_SIMPLE:
          mDynamicModelFile << "    [r, y, dr(" << count_call << ").g1, dr(" << count_call << ").g2, dr(" << count_call << ").g3, b, dr(" << count_call << ").g1_x, dr(" << count_call << ").g1_xd, dr(" << count_call << ").g1_o]=" << basename << ".block.dynamic_" <<  block + 1 << "(y, x, params, steady_state, it_-" << max_lag << ", 1, " << max_lag << ", " << block_recursive << "," << "options_.periods" << ");" << endl
                            << "    residual(y_index_eq)=r(:,M_.maximum_lag+1);" << endl;
          break;
        default:
          break;
        }
      tmp_eq.str("");
      tmp.str("");
    }
  if (tmp1.str().length())
    {
      mDynamicModelFile << tmp1.str();
      tmp1.str("");
    }
  mDynamicModelFile << "    varargout{1}=residual;" << endl
                    << "    varargout{2}=dr;" << endl
                    << "    return;" << endl
                    << "  end;" << endl
                    << "  %it is the deterministic simulation of the block decomposed dynamic model" << endl
                    << "  if(options_.stack_solve_algo==0)" << endl
                    << "    mthd='Sparse LU';" << endl
                    << "  elseif(options_.stack_solve_algo==1)" << endl
                    << "    mthd='Relaxation';" << endl
                    << "  elseif(options_.stack_solve_algo==2)" << endl
                    << "    mthd='GMRES';" << endl
                    << "  elseif(options_.stack_solve_algo==3)" << endl
                    << "    mthd='BICGSTAB';" << endl
                    << "  elseif(options_.stack_solve_algo==4)" << endl
                    << "    mthd='OPTIMPATH';" << endl
                    << "  else" << endl
                    << "    mthd='UNKNOWN';" << endl
                    << "  end;" << endl
                    << "  if options_.verbosity" << endl
                    << "    printline(41)" << endl
                    << "    disp(sprintf('MODEL SIMULATION (method=%s):',mthd))" << endl
                    << "    skipline()" << endl
                    << "  end" << endl
                    << "  periods=options_.periods;" << endl
                    << "  maxit_=options_.simul.maxit;" << endl
                    << "  solve_tolf=options_.solve_tolf;" << endl
                    << "  y=oo_.endo_simul';" << endl
                    << "  x=oo_.exo_simul;" << endl
                    << "  params=M_.params;" << endl
                    << "  steady_state=oo_.steady_state;" << endl
                    << "  oo_.deterministic_simulation.status = 0;" << endl;
  for (block = 0; block < nb_blocks; block++)
    {
      unsigned int block_size = getBlockSize(block);
      unsigned int block_mfs = getBlockMfs(block);
      unsigned int block_recursive = block_size - block_mfs;
      BlockSimulationType simulation_type = getBlockSimulationType(block);

      if ((simulation_type == EVALUATE_FORWARD) && (block_size))
        {
          if (open_par)
            mDynamicModelFile << "  end" << endl;
          mDynamicModelFile << "  oo_.deterministic_simulation.status = 1;" << endl
                            << "  oo_.deterministic_simulation.error = 0;" << endl
                            << "  oo_.deterministic_simulation.iterations = 0;" << endl
                            << "  if(isfield(oo_.deterministic_simulation,'block'))" << endl
                            << "    blck_num = length(oo_.deterministic_simulation.block)+1;" << endl
                            << "  else" << endl
                            << "    blck_num = 1;" << endl
                            << "  end;" << endl
                            << "  oo_.deterministic_simulation.block(blck_num).status = 1;" << endl
                            << "  oo_.deterministic_simulation.block(blck_num).error = 0;" << endl
                            << "  oo_.deterministic_simulation.block(blck_num).iterations = 0;" << endl
                            << "  g1=[];g2=[];g3=[];" << endl
                            << "  y=" << basename << ".block.dynamic_" << block + 1 << "(y, x, params, steady_state, 0, y_kmin, periods);" << endl
                            << "  tmp = y(:,M_.block_structure.block(" << block + 1 << ").variable);" << endl
                            << "  if any(isnan(tmp) | isinf(tmp))" << endl
                            << "    disp(['Inf or Nan value during the evaluation of block " << block <<"']);" << endl
                            << "    oo_.deterministic_simulation.status = 0;" << endl
                            << "    oo_.deterministic_simulation.error = 100;" << endl
                            << "    varargout{1} = oo_;" << endl
                            << "    return;" << endl
                            << "  end;" << endl;
        }
      else if ((simulation_type == EVALUATE_BACKWARD) && (block_size))
        {
          if (open_par)
            mDynamicModelFile << "  end" << endl;
          mDynamicModelFile << "  oo_.deterministic_simulation.status = 1;" << endl
                            << "  oo_.deterministic_simulation.error = 0;" << endl
                            << "  oo_.deterministic_simulation.iterations = 0;" << endl
                            << "  if(isfield(oo_.deterministic_simulation,'block'))" << endl
                            << "    blck_num = length(oo_.deterministic_simulation.block)+1;" << endl
                            << "  else" << endl
                            << "    blck_num = 1;" << endl
                            << "  end;" << endl
                            << "  oo_.deterministic_simulation.block(blck_num).status = 1;" << endl
                            << "  oo_.deterministic_simulation.block(blck_num).error = 0;" << endl
                            << "  oo_.deterministic_simulation.block(blck_num).iterations = 0;" << endl
                            << "  g1=[];g2=[];g3=[];" << endl
                            << "  " << basename << ".block.dynamic_" << block + 1 << "(y, x, params, steady_state, 0, y_kmin, periods);" << endl
                            << "  tmp = y(:,M_.block_structure.block(" << block + 1 << ").variable);" << endl
                            << "  if any(isnan(tmp) | isinf(tmp))" << endl
                            << "    disp(['Inf or Nan value during the evaluation of block " << block <<"']);" << endl
                            << "    oo_.deterministic_simulation.status = 0;" << endl
                            << "    oo_.deterministic_simulation.error = 100;" << endl
                            << "    varargout{1} = oo_;" << endl
                            << "    return;" << endl
                            << "  end;" << endl;
        }
      else if ((simulation_type == SOLVE_FORWARD_COMPLETE || simulation_type == SOLVE_FORWARD_SIMPLE) && (block_size))
        {
          if (open_par)
            mDynamicModelFile << "  end" << endl;
          open_par = false;
          mDynamicModelFile << "  g1=0;" << endl
                            << "  r=0;" << endl;
          tmp.str("");
          for (unsigned int ik = block_recursive; ik < block_size; ik++)
            tmp << " " << getBlockVariableID(block, ik)+1;
          mDynamicModelFile << "  y_index = [" << tmp.str() << "];" << endl;
          int nze = blocks_derivatives[block].size();
          mDynamicModelFile << "  if(isfield(oo_.deterministic_simulation,'block'))" << endl
                            << "    blck_num = length(oo_.deterministic_simulation.block)+1;" << endl
                            << "  else" << endl
                            << "    blck_num = 1;" << endl
                            << "  end;" << endl
                            << "  y = solve_one_boundary('" << basename << ".block.dynamic_" <<  block + 1 << "'"
                            << ", y, x, params, steady_state, y_index, " << nze
                            << ", options_.periods, " << blocks_linear[block]
                            << ", blck_num, y_kmin, options_.simul.maxit, options_.solve_tolf, options_.slowc, " << cutoff << ", options_.stack_solve_algo, 1, 1, 0);" << endl
                            << "  tmp = y(:,M_.block_structure.block(" << block + 1 << ").variable);" << endl
                            << "  if any(isnan(tmp) | isinf(tmp))" << endl
                            << "    disp(['Inf or Nan value during the resolution of block " << block <<"']);" << endl
                            << "    oo_.deterministic_simulation.status = 0;" << endl
                            << "    oo_.deterministic_simulation.error = 100;" << endl
                            << "    varargout{1} = oo_;" << endl
                            << "    return;" << endl
                            << "  end;" << endl;
        }
      else if ((simulation_type == SOLVE_BACKWARD_COMPLETE || simulation_type == SOLVE_BACKWARD_SIMPLE) && (block_size))
        {
          if (open_par)
            mDynamicModelFile << "  end" << endl;
          open_par = false;
          mDynamicModelFile << "  g1=0;" << endl
                            << "  r=0;" << endl;
          tmp.str("");
          for (unsigned int ik = block_recursive; ik < block_size; ik++)
            tmp << " " << getBlockVariableID(block, ik)+1;
          mDynamicModelFile << "  y_index = [" << tmp.str() << "];" << endl;
          int nze = blocks_derivatives[block].size();

          mDynamicModelFile << "  if(isfield(oo_.deterministic_simulation,'block'))" << endl
                            << "    blck_num = length(oo_.deterministic_simulation.block)+1;" << endl
                            << "  else" << endl
                            << "    blck_num = 1;" << endl
                            << "  end;" << endl
                            << "  y = solve_one_boundary('" << basename << ".block.dynamic_" <<  block + 1 << "'"
                            <<", y, x, params, steady_state, y_index, " << nze
                            <<", options_.periods, " << blocks_linear[block]
                            <<", blck_num, y_kmin, options_.simul.maxit, options_.solve_tolf, options_.slowc, " << cutoff << ", options_.stack_solve_algo, 1, 1, 0);" << endl
                            << "  tmp = y(:,M_.block_structure.block(" << block + 1 << ").variable);" << endl
                            << "  if any(isnan(tmp) | isinf(tmp))" << endl
                            << "    disp(['Inf or Nan value during the resolution of block " << block <<"']);" << endl
                            << "    oo_.deterministic_simulation.status = 0;" << endl
                            << "    oo_.deterministic_simulation.error = 100;" << endl
                            << "    varargout{1} = oo_;" << endl
                            << "    return;" << endl
                            << "  end;" << endl;
        }
      else if ((simulation_type == SOLVE_TWO_BOUNDARIES_COMPLETE || simulation_type == SOLVE_TWO_BOUNDARIES_SIMPLE) && (block_size))
        {
          if (open_par)
            mDynamicModelFile << "  end" << endl;
          open_par = false;
          Nb_SGE++;
          int nze = blocks_derivatives[block].size();
          mDynamicModelFile << "  y_index=[";
          for (unsigned int ik = block_recursive; ik < block_size; ik++)
            mDynamicModelFile << " " << getBlockVariableID(block, ik)+1;
          mDynamicModelFile << "  ];" << endl
                            << "  if(isfield(oo_.deterministic_simulation,'block'))" << endl
                            << "    blck_num = length(oo_.deterministic_simulation.block)+1;" << endl
                            << "  else" << endl
                            << "    blck_num = 1;" << endl
                            << "  end;" << endl
                            << "  [y oo_] = solve_two_boundaries('" << basename << ".block.dynamic_" <<  block + 1 << "'"
                            <<", y, x, params, steady_state, y_index, " << nze
                            <<", options_.periods, " << max_leadlag_block[block].first
                            <<", " << max_leadlag_block[block].second
                            <<", " << blocks_linear[block]
                            <<", blck_num, y_kmin, options_.simul.maxit, options_.solve_tolf, options_.slowc, " << cutoff << ", options_.stack_solve_algo, options_, M_, oo_);" << endl
                            << "  tmp = y(:,M_.block_structure.block(" << block + 1 << ").variable);" << endl
                            << "  if any(isnan(tmp) | isinf(tmp))" << endl
                            << "    disp(['Inf or Nan value during the resolution of block " << block <<"']);" << endl
                            << "    oo_.deterministic_simulation.status = 0;" << endl
                            << "    oo_.deterministic_simulation.error = 100;" << endl
                            << "    varargout{1} = oo_;" << endl
                            << "    return;" << endl
                            << "  end;" << endl;
        }
    }
  if (open_par)
    mDynamicModelFile << "  end;" << endl;
  open_par = false;
  mDynamicModelFile << "  oo_.endo_simul = y';" << endl
                    << "  varargout{1} = oo_;" << endl
                    << "return;" << endl
                    << "end" << endl;

  mDynamicModelFile.close();

  writeModelEquationsOrdered_M(basename);
}

void
DynamicModel::writeWrapperFunctions(const string &basename, const string &ending) const
{
  string name;
  if (ending == "g1")
    name = "dynamic_resid_g1";
  else if (ending == "g2")
    name= "dynamic_resid_g1_g2";
  else if (ending == "g3")
    name = "dynamic_resid_g1_g2_g3";

  string filename = packageDir(basename) + "/" + name + ".m";
  ofstream output;
  output.open(filename, ios::out | ios::binary);
  if (!output.is_open())
    {
      cerr << "Error: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }

  if (ending == "g1")
    output << "function [residual, g1] = " << name << "(T, y, x, params, steady_state, it_, T_flag)" << endl
           << "% function [residual, g1] = " << name << "(T, y, x, params, steady_state, it_, T_flag)" << endl;
  else if (ending == "g2")
    output << "function [residual, g1, g2] = " << name << "(T, y, x, params, steady_state, it_, T_flag)" << endl
           << "% function [residual, g1, g2] = " << name << "(T, y, x, params, steady_state, it_, T_flag)" << endl;
  else if (ending == "g3")
    output << "function [residual, g1, g2, g3] = " << name << "(T, y, x, params, steady_state, it_, T_flag)" << endl
           << "% function [residual, g1, g2, g3] = " << name << "(T, y, x, params, steady_state, it_, T_flag)" << endl;

  output << "%" << endl
         << "% Wrapper function automatically created by Dynare" << endl
         << "%" << endl
         << endl
         << "    if T_flag" << endl
         << "        T = " << basename << ".dynamic_" << ending << "_tt(T, y, x, params, steady_state, it_);" << endl
         << "    end" << endl;

  if (ending == "g1")
    output << "    residual = " << basename << ".dynamic_resid(T, y, x, params, steady_state, it_, false);" << endl
           << "    g1       = " << basename << ".dynamic_g1(T, y, x, params, steady_state, it_, false);" << endl;
  else if (ending == "g2")
    output << "    [residual, g1] = " << basename << ".dynamic_resid_g1(T, y, x, params, steady_state, it_, false);" << endl
           << "    g2       = " << basename << ".dynamic_g2(T, y, x, params, steady_state, it_, false);" << endl;
  else if (ending == "g3")
    output << "    [residual, g1, g2] = " << basename << ".dynamic_resid_g1_g2(T, y, x, params, steady_state, it_, false);" << endl
           << "    g3       = " << basename << ".dynamic_g3(T, y, x, params, steady_state, it_, false);" << endl;

  output << endl << "end" << endl;
  output.close();
}

void
DynamicModel::writeDynamicModelHelper(const string &basename,
                                      const string &name, const string &retvalname,
                                      const string &name_tt, size_t ttlen,
                                      const string &previous_tt_name,
                                      const ostringstream &init_s,
                                      const ostringstream &end_s,
                                      const ostringstream &s, const ostringstream &s_tt) const
{
  string filename = packageDir(basename) + "/" + name_tt + ".m";
  ofstream output;
  output.open(filename, ios::out | ios::binary);
  if (!output.is_open())
    {
      cerr << "Error: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }

  output << "function T = " << name_tt << "(T, y, x, params, steady_state, it_)" << endl
         << "% function T = " << name_tt << "(T, y, x, params, steady_state, it_)" << endl
         << "%" << endl
         << "% File created by Dynare Preprocessor from .mod file" << endl
         << "%" << endl
         << "% Inputs:" << endl
         << "%   T             [#temp variables by 1]     double  vector of temporary terms to be filled by function" << endl
         << "%   y             [#dynamic variables by 1]  double  vector of endogenous variables in the order stored" << endl
         << "%                                                    in M_.lead_lag_incidence; see the Manual" << endl
         << "%   x             [nperiods by M_.exo_nbr]   double  matrix of exogenous variables (in declaration order)" << endl
         << "%                                                    for all simulation periods" << endl
         << "%   steady_state  [M_.endo_nbr by 1]         double  vector of steady state values" << endl
         << "%   params        [M_.param_nbr by 1]        double  vector of parameter values in declaration order" << endl
         << "%   it_           scalar                     double  time period for exogenous variables for which" << endl
         << "%                                                    to evaluate the model" << endl
         << "%" << endl
         << "% Output:" << endl
         << "%   T           [#temp variables by 1]       double  vector of temporary terms" << endl
         << "%" << endl << endl
         << "assert(length(T) >= " << ttlen << ");" << endl
         << endl;

  if (!previous_tt_name.empty())
    output << "T = " << basename << "." << previous_tt_name << "(T, y, x, params, steady_state, it_);" << endl << endl;

  output << s_tt.str() << endl
         << "end" << endl;
  output.close();

  filename = packageDir(basename) + "/" + name + ".m";
  output.open(filename, ios::out | ios::binary);
  if (!output.is_open())
    {
      cerr << "Error: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }

  output << "function " << retvalname << " = " << name << "(T, y, x, params, steady_state, it_, T_flag)" << endl
         << "% function " << retvalname << " = " << name << "(T, y, x, params, steady_state, it_, T_flag)" << endl
         << "%" << endl
         << "% File created by Dynare Preprocessor from .mod file" << endl
         << "%" << endl
         << "% Inputs:" << endl
         << "%   T             [#temp variables by 1]     double   vector of temporary terms to be filled by function" << endl
         << "%   y             [#dynamic variables by 1]  double   vector of endogenous variables in the order stored" << endl
         << "%                                                     in M_.lead_lag_incidence; see the Manual" << endl
         << "%   x             [nperiods by M_.exo_nbr]   double   matrix of exogenous variables (in declaration order)" << endl
         << "%                                                     for all simulation periods" << endl
         << "%   steady_state  [M_.endo_nbr by 1]         double   vector of steady state values" << endl
         << "%   params        [M_.param_nbr by 1]        double   vector of parameter values in declaration order" << endl
         << "%   it_           scalar                     double   time period for exogenous variables for which" << endl
         << "%                                                     to evaluate the model" << endl
         << "%   T_flag        boolean                    boolean  flag saying whether or not to calculate temporary terms" << endl
         << "%" << endl
         << "% Output:" << endl
         << "%   " << retvalname << endl
         << "%" << endl << endl;

  if (!name_tt.empty())
    output << "if T_flag" << endl
           << "    T = " << basename << "." << name_tt << "(T, y, x, params, steady_state, it_);" << endl
           << "end" << endl;

  output << init_s.str() << endl
         << s.str()
         << end_s.str() << endl
         << "end" << endl;
  output.close();
}

void
DynamicModel::writeDynamicMatlabCompatLayer(const string &basename) const
{
  string filename = packageDir(basename) + "/dynamic.m";
  ofstream output;
  output.open(filename, ios::out | ios::binary);
  if (!output.is_open())
    {
      cerr << "Error: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }
  int ntt = temporary_terms_mlv.size() + temporary_terms_derivatives[0].size() + temporary_terms_derivatives[1].size() + temporary_terms_derivatives[2].size() + temporary_terms_derivatives[3].size();

  output << "function [residual, g1, g2, g3] = dynamic(y, x, params, steady_state, it_)" << endl
         << "    T = NaN(" << ntt << ", 1);" << endl
         << "    if nargout <= 1" << endl
         << "        residual = " << basename << ".dynamic_resid(T, y, x, params, steady_state, it_, true);" << endl
         << "    elseif nargout == 2" << endl
         << "        [residual, g1] = " << basename << ".dynamic_resid_g1(T, y, x, params, steady_state, it_, true);" << endl
         << "    elseif nargout == 3" << endl
         << "        [residual, g1, g2] = " << basename << ".dynamic_resid_g1_g2(T, y, x, params, steady_state, it_, true);" << endl
         << "    else" << endl
         << "        [residual, g1, g2, g3] = " << basename << ".dynamic_resid_g1_g2_g3(T, y, x, params, steady_state, it_, true);" << endl
         << "    end" << endl
         << "end" << endl;

  output.close();
}

void
DynamicModel::writeDynamicModel(ostream &DynamicOutput, bool use_dll, bool julia) const
{
  writeDynamicModel("", DynamicOutput, use_dll, julia);
}

void
DynamicModel::writeDynamicModel(const string &basename, bool use_dll, bool julia) const
{
  ofstream DynamicOutput;
  writeDynamicModel(basename, DynamicOutput, use_dll, julia);
}

void
DynamicModel::writeDynamicModel(const string &basename, ostream &DynamicOutput, bool use_dll, bool julia) const
{
  ostringstream model_tt_output;             // Used for storing model temp vars
  ostringstream model_output;                // Used for storing model equations
  ostringstream jacobian_tt_output;          // Used for storing jacobian temp vars
  ostringstream jacobian_output;             // Used for storing jacobian equations
  ostringstream hessian_tt_output;           // Used for storing Hessian temp vars
  ostringstream hessian_output;              // Used for storing Hessian equations
  ostringstream third_derivatives_tt_output; // Used for storing third order derivatives temp terms
  ostringstream third_derivatives_output;    // Used for storing third order derivatives equations

  ExprNodeOutputType output_type = (use_dll ? ExprNodeOutputType::CDynamicModel :
                                    julia ? ExprNodeOutputType::juliaDynamicModel : ExprNodeOutputType::matlabDynamicModel);

  deriv_node_temp_terms_t tef_terms;
  temporary_terms_t temp_term_union;

  writeModelLocalVariableTemporaryTerms(temp_term_union, temporary_terms_idxs,
                                        model_tt_output, output_type, tef_terms);

  writeTemporaryTerms(temporary_terms_derivatives[0],
                      temp_term_union,
                      temporary_terms_idxs,
                      model_tt_output, output_type, tef_terms);

  writeModelEquations(model_output, output_type, temp_term_union);

  int nrows = equations.size();
  int hessianColsNbr = dynJacobianColsNbr * dynJacobianColsNbr;

  // Writing Jacobian
  if (!derivatives[1].empty())
    {
      writeTemporaryTerms(temporary_terms_derivatives[1],
                          temp_term_union,
                          temporary_terms_idxs,
                          jacobian_tt_output, output_type, tef_terms);

      for (const auto & first_derivative : derivatives[1])
        {
          int eq, var;
          tie(eq, var) = vectorToTuple<2>(first_derivative.first);
          expr_t d1 = first_derivative.second;

          jacobianHelper(jacobian_output, eq, getDynJacobianCol(var), output_type);
          jacobian_output << "=";
          d1->writeOutput(jacobian_output, output_type,
                          temp_term_union, temporary_terms_idxs, tef_terms);
          jacobian_output << ";" << endl;
        }
    }

  // Writing Hessian
  if (!derivatives[2].empty())
    {
      writeTemporaryTerms(temporary_terms_derivatives[2],
                          temp_term_union,
                          temporary_terms_idxs,
                          hessian_tt_output, output_type, tef_terms);

      /* When creating the sparse matrix (in MATLAB or C mode), since storage
         is in column-major order, output the first column, then the second,
         then the third. This gives a significant performance boost in use_dll
         mode (at both compilation and runtime), because it facilitates memory
         accesses and expression reusage. */
      ostringstream col0_output, col1_output, col2_output;

      int k = 0; // Keep the line of a 2nd derivative in v2
      for (const auto & second_derivative : derivatives[2])
        {
          int eq, var1, var2;
          tie(eq, var1, var2) = vectorToTuple<3>(second_derivative.first);
          expr_t d2 = second_derivative.second;

          int id1 = getDynJacobianCol(var1);
          int id2 = getDynJacobianCol(var2);

          int col_nb = id1 * dynJacobianColsNbr + id2;
          int col_nb_sym = id2 * dynJacobianColsNbr + id1;

          ostringstream for_sym;
          if (output_type == ExprNodeOutputType::juliaDynamicModel)
            {
              for_sym << "g2[" << eq + 1 << "," << col_nb + 1 << "]";
              hessian_output << "    @inbounds " << for_sym.str() << " = ";
              d2->writeOutput(hessian_output, output_type, temp_term_union, temporary_terms_idxs, tef_terms);
              hessian_output << endl;
            }
          else
            {
              sparseHelper(2, col0_output, k, 0, output_type);
              col0_output << "=" << eq + 1 << ";" << endl;

              sparseHelper(2, col1_output, k, 1, output_type);
              col1_output << "=" << col_nb + 1 << ";" << endl;

              sparseHelper(2, col2_output, k, 2, output_type);
              col2_output << "=";
              d2->writeOutput(col2_output, output_type, temp_term_union, temporary_terms_idxs, tef_terms);
              col2_output << ";" << endl;

              k++;
            }

          // Treating symetric elements
          if (id1 != id2)
            if (output_type == ExprNodeOutputType::juliaDynamicModel)
              hessian_output << "    @inbounds g2[" << eq + 1 << "," << col_nb_sym + 1 << "] = "
                             << for_sym.str() << endl;
            else
              {
                sparseHelper(2, col0_output, k, 0, output_type);
                col0_output << "=" << eq + 1 << ";" << endl;

                sparseHelper(2, col1_output, k, 1, output_type);
                col1_output << "=" << col_nb_sym + 1 << ";" << endl;

                sparseHelper(2, col2_output, k, 2, output_type);
                col2_output << "=";
                sparseHelper(2, col2_output, k-1, 2, output_type);
                col2_output << ";" << endl;

                k++;
              }
        }

      if (output_type != ExprNodeOutputType::juliaDynamicModel)
        hessian_output << col0_output.str() << col1_output.str() << col2_output.str();
    }

  // Writing third derivatives
  if (!derivatives[3].empty())
    {
      writeTemporaryTerms(temporary_terms_derivatives[3],
                          temp_term_union,
                          temporary_terms_idxs,
                          third_derivatives_tt_output, output_type, tef_terms);

      // See comment above for 2nd order
      ostringstream col0_output, col1_output, col2_output;

      int k = 0; // Keep the line of a 3rd derivative in v3
      for (const auto & third_derivative : derivatives[3])
        {
          int eq, var1, var2, var3;
          tie(eq, var1, var2, var3) = vectorToTuple<4>(third_derivative.first);
          expr_t d3 = third_derivative.second;

          int id1 = getDynJacobianCol(var1);
          int id2 = getDynJacobianCol(var2);
          int id3 = getDynJacobianCol(var3);

          // Reference column number for the g3 matrix
          int ref_col = id1 * hessianColsNbr + id2 * dynJacobianColsNbr + id3;

          ostringstream for_sym;
          if (output_type == ExprNodeOutputType::juliaDynamicModel)
            {
              for_sym << "g3[" << eq + 1 << "," << ref_col + 1 << "]";
              third_derivatives_output << "    @inbounds " << for_sym.str() << " = ";
              d3->writeOutput(third_derivatives_output, output_type, temp_term_union, temporary_terms_idxs, tef_terms);
              third_derivatives_output << endl;
            }
          else
            {
              sparseHelper(3, col0_output, k, 0, output_type);
              col0_output << "=" << eq + 1 << ";" << endl;

              sparseHelper(3, col1_output, k, 1, output_type);
              col1_output << "=" << ref_col + 1 << ";" << endl;

              sparseHelper(3, col2_output, k, 2, output_type);
              col2_output << "=";
              d3->writeOutput(col2_output, output_type, temp_term_union, temporary_terms_idxs, tef_terms);
              col2_output << ";" << endl;
            }

          // Compute the column numbers for the 5 other permutations of (id1,id2,id3)
          // and store them in a set (to avoid duplicates if two indexes are equal)
          set<int> cols;
          cols.insert(id1 * hessianColsNbr + id3 * dynJacobianColsNbr + id2);
          cols.insert(id2 * hessianColsNbr + id1 * dynJacobianColsNbr + id3);
          cols.insert(id2 * hessianColsNbr + id3 * dynJacobianColsNbr + id1);
          cols.insert(id3 * hessianColsNbr + id1 * dynJacobianColsNbr + id2);
          cols.insert(id3 * hessianColsNbr + id2 * dynJacobianColsNbr + id1);

          int k2 = 1; // Keeps the offset of the permutation relative to k
          for (int col : cols)
            if (col != ref_col)
              if (output_type == ExprNodeOutputType::juliaDynamicModel)
                third_derivatives_output << "    @inbounds g3[" << eq + 1 << "," << col + 1 << "] = "
                                         << for_sym.str() << endl;
              else
                {
                  sparseHelper(3, col0_output, k+k2, 0, output_type);
                  col0_output << "=" << eq + 1 << ";" << endl;

                  sparseHelper(3, col1_output, k+k2, 1, output_type);
                  col1_output << "=" << col + 1 << ";" << endl;

                  sparseHelper(3, col2_output, k+k2, 2, output_type);
                  col2_output << "=";
                  sparseHelper(3, col2_output, k, 2, output_type);
                  col2_output << ";" << endl;

                  k2++;
                }
          k += k2;
        }

      if (output_type != ExprNodeOutputType::juliaDynamicModel)
        third_derivatives_output << col0_output.str() << col1_output.str() << col2_output.str();
    }

  if (output_type == ExprNodeOutputType::matlabDynamicModel)
    {
      // Check that we don't have more than 32 nested parenthesis because Matlab does not suppor this. See Issue #1201
      map<string, string> tmp_paren_vars;
      bool message_printed = false;
      fixNestedParenthesis(model_output, tmp_paren_vars, message_printed);
      fixNestedParenthesis(model_tt_output, tmp_paren_vars, message_printed);
      fixNestedParenthesis(jacobian_output, tmp_paren_vars, message_printed);
      fixNestedParenthesis(jacobian_tt_output, tmp_paren_vars, message_printed);
      fixNestedParenthesis(hessian_output, tmp_paren_vars, message_printed);
      fixNestedParenthesis(hessian_tt_output, tmp_paren_vars, message_printed);
      fixNestedParenthesis(third_derivatives_output, tmp_paren_vars, message_printed);
      fixNestedParenthesis(third_derivatives_tt_output, tmp_paren_vars, message_printed);

      ostringstream init_output, end_output;
      init_output << "residual = zeros(" << nrows << ", 1);";
      writeDynamicModelHelper(basename, "dynamic_resid", "residual",
                              "dynamic_resid_tt",
                              temporary_terms_mlv.size() + temporary_terms_derivatives[0].size(),
                              "", init_output, end_output,
                              model_output, model_tt_output);

      init_output.str(string());
      init_output.clear();
      init_output << "g1 = zeros(" << nrows << ", " << dynJacobianColsNbr << ");";
      writeDynamicModelHelper(basename, "dynamic_g1", "g1",
                              "dynamic_g1_tt",
                              temporary_terms_mlv.size() + temporary_terms_derivatives[0].size() + temporary_terms_derivatives[1].size(),
                              "dynamic_resid_tt",
                              init_output, end_output,
                              jacobian_output, jacobian_tt_output);
      writeWrapperFunctions(basename, "g1");

      init_output.str(string());
      init_output.clear();
      if (derivatives[2].size())
        {
          init_output << "v2 = zeros(" << NNZDerivatives[2] << ",3);";
          end_output << "g2 = sparse(v2(:,1),v2(:,2),v2(:,3)," << nrows << "," << hessianColsNbr << ");";
        }
      else
        init_output << "g2 = sparse([],[],[]," << nrows << "," << hessianColsNbr << ");";
      writeDynamicModelHelper(basename, "dynamic_g2", "g2",
                              "dynamic_g2_tt",
                              temporary_terms_mlv.size() + temporary_terms_derivatives[0].size() + temporary_terms_derivatives[1].size()
                              + temporary_terms_derivatives[2].size(),
                              "dynamic_g1_tt",
                              init_output, end_output,
                              hessian_output, hessian_tt_output);
      writeWrapperFunctions(basename, "g2");

      init_output.str(string());
      init_output.clear();
      end_output.str(string());
      end_output.clear();
      int ncols = hessianColsNbr * dynJacobianColsNbr;
      if (derivatives[3].size())
        {
          init_output << "v3 = zeros(" << NNZDerivatives[3] << ",3);";
          end_output << "g3 = sparse(v3(:,1),v3(:,2),v3(:,3)," << nrows << "," << ncols << ");";
        }
      else
        init_output << "g3 = sparse([],[],[]," << nrows << "," << ncols << ");";
      writeDynamicModelHelper(basename, "dynamic_g3", "g3",
                              "dynamic_g3_tt",
                              temporary_terms_mlv.size() + temporary_terms_derivatives[0].size() + temporary_terms_derivatives[1].size()
                              + temporary_terms_derivatives[2].size() + temporary_terms_derivatives[3].size(),
                              "dynamic_g2_tt",
                              init_output, end_output,
                              third_derivatives_output, third_derivatives_tt_output);
      writeWrapperFunctions(basename, "g3");

      writeDynamicMatlabCompatLayer(basename);
    }
  else if (output_type == ExprNodeOutputType::CDynamicModel)
    {
      DynamicOutput << "void dynamic_resid_tt(const double *y, const double *x, int nb_row_x, const double *params, const double *steady_state, int it_, double *T)" << endl
                    << "{" << endl
                    << model_tt_output.str()
                    << "}" << endl
                    << endl
                    << "void dynamic_resid(const double *y, const double *x, int nb_row_x, const double *params, const double *steady_state, int it_, const double *T, double *residual)" << endl
                    << "{" << endl
                    << "  double lhs, rhs;" << endl
                    << model_output.str()
                    << "}" << endl
                    << endl
                    << "void dynamic_g1_tt(const double *y, const double *x, int nb_row_x, const double *params, const double *steady_state, int it_, double *T)" << endl
                    << "{" << endl
                    << jacobian_tt_output.str()
                    << "}" << endl
                    << endl
                    << "void dynamic_g1(const double *y, const double *x, int nb_row_x, const double *params, const double *steady_state, int it_, const double *T, double *g1)" << endl
                    << "{" << endl
                    << jacobian_output.str()
                    << "}" << endl
                    << endl
                    << "void dynamic_g2_tt(const double *y, const double *x, int nb_row_x, const double *params, const double *steady_state, int it_, double *T)" << endl
                    << "{" << endl
                    << hessian_tt_output.str()
                    << "}" << endl
                    << endl
                    << "void dynamic_g2(const double *y, const double *x, int nb_row_x, const double *params, const double *steady_state, int it_, const double *T, double *v2)" << endl
                    << "{" << endl
                    << hessian_output.str()
                    << "}" << endl
                    << endl
                    << "void dynamic_g3_tt(const double *y, const double *x, int nb_row_x, const double *params, const double *steady_state, int it_, double *T)" << endl
                    << "{" << endl
                    << third_derivatives_tt_output.str()
                    << "}" << endl
                    << endl
                    << "void dynamic_g3(const double *y, const double *x, int nb_row_x, const double *params, const double *steady_state, int it_, const double *T, double *v3)" << endl
                    << "{" << endl
                    << third_derivatives_output.str()
                    << "}" << endl;
    }
  else
    {
      string filename = basename + "Dynamic.jl";
      ofstream output;
      output.open(filename, ios::out | ios::binary);
      if (!output.is_open())
        {
          cerr << "Error: Can't open file " << filename << " for writing" << endl;
          exit(EXIT_FAILURE);
        }

      output << "module " << basename << "Dynamic" << endl
             << "#" << endl
             << "# NB: this file was automatically generated by Dynare" << endl
             << "#     from " << basename << ".mod" << endl
             << "#" << endl
             << "using Utils" << endl << endl
             << "export tmp_nbr, dynamic!, dynamicResid!, dynamicG1!, dynamicG2!, dynamicG3!" << endl << endl
             << "#=" << endl
             << "# The comments below apply to all functions contained in this module #" << endl
             << "  NB: The arguments contained on the first line of the function" << endl
             << "      definition are those that are modified in place" << endl << endl
             << "## Exported Functions ##" << endl
             << "  dynamic!      : Wrapper function; computes residuals, Jacobian, Hessian," << endl
             << "                  and third derivatives depending on the arguments provided" << endl
             << "  dynamicResid! : Computes the dynamic model residuals" << endl
             << "  dynamicG1!    : Computes the dynamic model Jacobian" << endl
             << "  dynamicG2!    : Computes the dynamic model Hessian" << endl
             << "  dynamicG3!    : Computes the dynamic model third derivatives" << endl << endl
             << "## Exported Variables ##" << endl
             << "  tmp_nbr       : Vector{Int}(4) respectively the number of temporary variables" << endl
             << "                  for the residuals, g1, g2 and g3." << endl << endl
             << "## Local Functions ##" << endl
             << "  dynamicResidTT! : Computes the dynamic model temporary terms for the residuals" << endl
             << "  dynamicG1TT!    : Computes the dynamic model temporary terms for the Jacobian" << endl
             << "  dynamicG2TT!    : Computes the dynamic model temporary terms for the Hessian" << endl
             << "  dynamicG3TT!    : Computes the dynamic model temporary terms for the third derivatives" << endl << endl
             << "## Function Arguments ##" << endl
             << "  T            : Vector{Float64}(num_temp_terms), temporary terms" << endl
             << "  y            : Vector{Float64}(num_dynamic_vars), endogenous variables in the order stored model_.lead_lag_incidence; see the manual" << endl
             << "  x            : Matrix{Float64}(nperiods,model_.exo_nbr), exogenous variables (in declaration order) for all simulation periods" << endl
             << "  params       : Vector{Float64}(model_.param_nbr), parameter values in declaration order" << endl
             << "  steady_state : Vector{Float64}(model_endo_nbr)" << endl
             << "  it_          : Int, time period for exogenous variables for which to evaluate the model" << endl
             << "  residual     : Vector{Float64}(model_.eq_nbr), residuals of the dynamic model equations in order of declaration of the equations." << endl
             << "  g1           : Matrix{Float64}(model_.eq_nbr, num_dynamic_vars), Jacobian matrix of the dynamic model equations" << endl
             << "                 The rows and columns respectively correspond to equations in order of declaration and variables in order" << endl
             << "                 stored in model_.lead_lag_incidence" << endl
             << "  g2           : spzeros(model_.eq_nbr, (num_dynamic_vars)^2) Hessian matrix of the dynamic model equations" << endl
             << "                 The rows and columns respectively correspond to equations in order of declaration and variables in order" << endl
             << "                 stored in model_.lead_lag_incidence" << endl
             << "  g3           : spzeros(model_.eq_nbr, (num_dynamic_vars)^3) Third order derivative matrix of the dynamic model equations;" << endl
             << "                 The rows and columns respectively correspond to equations in order of declaration and variables in order" << endl
             << "                 stored in model_.lead_lag_incidence" << endl << endl
             << "## Remarks ##" << endl
             << "  [1] `num_dynamic_vars` is the number of non zero entries in the lead lag incidence matrix, `model_.lead_lag_incidence.`" << endl
             << "  [2] The size of `T`, ie the value of `num_temp_terms`, depends on the version of the dynamic model called. The number of temporary variables" << endl    
             << "      used for the different returned objects (residuals, jacobian, hessian or third order derivatives) is given by the elements in `tmp_nbr`" << endl
             << "      exported vector. The first element is the number of temporaries used for the computation of the residuals, the second element is the" << endl
             << "      number of temporaries used for the evaluation of the jacobian matrix, etc. If one calls the version of the dynamic model computing the" << endl
             << "      residuals, the jacobian and hessian matrices, then `T` must have at least `sum(tmp_nbr[1:3])` elements." << endl
             << "=#" << endl << endl;

      // Write the number of temporary terms
      output << "tmp_nbr = zeros(Int,4)" << endl
             << "tmp_nbr[1] = " << temporary_terms_mlv.size() + temporary_terms_derivatives[0].size() << "# Number of temporary terms for the residuals" << endl
             << "tmp_nbr[2] = " << temporary_terms_derivatives[1].size() << "# Number of temporary terms for g1 (jacobian)" << endl
             << "tmp_nbr[3] = " << temporary_terms_derivatives[2].size() << "# Number of temporary terms for g2 (hessian)" << endl
             << "tmp_nbr[4] = " << temporary_terms_derivatives[3].size() << "# Number of temporary terms for g3 (third order derivates)" << endl << endl;

      // dynamicResidTT!
      output << "function dynamicResidTT!(T::Vector{Float64}," << endl
             << "                         y::Vector{Float64}, x::Matrix{Float64}, "
             << "params::Vector{Float64}, steady_state::Vector{Float64}, it_::Int)" << endl
             << model_tt_output.str()
             << "    return nothing" << endl
             << "end" << endl << endl;

      // dynamic!
      output << "function dynamicResid!(T::Vector{Float64}, residual::Vector{Float64}," << endl
             << "                       y::Vector{Float64}, x::Matrix{Float64}, "
             << "params::Vector{Float64}, steady_state::Vector{Float64}, it_::Int, T_flag::Bool)" << endl
             << "    @assert length(T) >= " << temporary_terms_mlv.size() + temporary_terms_derivatives[0].size() << endl
             << "    @assert length(residual) == " << nrows << endl
             << "    @assert length(y)+size(x, 2) == " << dynJacobianColsNbr << endl
             << "    @assert length(params) == " << symbol_table.param_nbr() << endl
             << "    if T_flag" << endl
             << "        dynamicResidTT!(T, y, x, params, steady_state, it_)" << endl
             << "    end" << endl
             << model_output.str()
             << "    return nothing" << endl
             << "end" << endl << endl;

      // dynamicG1TT!
      output << "function dynamicG1TT!(T::Vector{Float64}," << endl
             << "                      y::Vector{Float64}, x::Matrix{Float64}, "
             << "params::Vector{Float64}, steady_state::Vector{Float64}, it_::Int)" << endl
             << "    dynamicResidTT!(T, y, x, params, steady_state, it_)" << endl
             << jacobian_tt_output.str()
             << "    return nothing" << endl
             << "end" << endl << endl;

      // dynamicG1!
      output << "function dynamicG1!(T::Vector{Float64}, g1::Matrix{Float64}," << endl
             << "                    y::Vector{Float64}, x::Matrix{Float64}, "
             << "params::Vector{Float64}, steady_state::Vector{Float64}, it_::Int, T_flag::Bool)" << endl
             << "    @assert length(T) >= "
             << temporary_terms_mlv.size() + temporary_terms_derivatives[0].size() + temporary_terms_derivatives[1].size() << endl
             << "    @assert size(g1) == (" << nrows << ", " << dynJacobianColsNbr << ")" << endl
             << "    @assert length(y)+size(x, 2) == " << dynJacobianColsNbr << endl
             << "    @assert length(params) == " << symbol_table.param_nbr() << endl
             << "    if T_flag" << endl
             << "        dynamicG1TT!(T, y, x, params, steady_state, it_)" << endl
             << "    end" << endl
             << "    fill!(g1, 0.0)" << endl
             << jacobian_output.str()
             << "    return nothing" << endl
             << "end" << endl << endl;

      // dynamicG2TT!
      output << "function dynamicG2TT!(T::Vector{Float64}," << endl
             << "                      y::Vector{Float64}, x::Matrix{Float64}, "
             << "params::Vector{Float64}, steady_state::Vector{Float64}, it_::Int)" << endl
             << "    dynamicG1TT!(T, y, x, params, steady_state, it_)" << endl
             << hessian_tt_output.str()
             << "    return nothing" << endl
             << "end" << endl << endl;

      // dynamicG2!
      output << "function dynamicG2!(T::Vector{Float64}, g2::Matrix{Float64}," << endl
             << "                    y::Vector{Float64}, x::Matrix{Float64}, "
             << "params::Vector{Float64}, steady_state::Vector{Float64}, it_::Int, T_flag::Bool)" << endl
             << "    @assert length(T) >= " << temporary_terms_mlv.size() + temporary_terms_derivatives[0].size() + temporary_terms_derivatives[1].size() + temporary_terms_derivatives[2].size() << endl
             << "    @assert size(g2) == (" << nrows << ", " << hessianColsNbr << ")" << endl
             << "    @assert length(y)+size(x, 2) == " << dynJacobianColsNbr << endl
             << "    @assert length(params) == " << symbol_table.param_nbr() << endl
             << "    if T_flag" << endl
             << "        dynamicG2TT!(T, y, x, params, steady_state, it_)" << endl
             << "    end" << endl
             << "    fill!(g2, 0.0)" << endl
             << hessian_output.str()
             << "    return nothing" << endl
             << "end" << endl << endl;

      // dynamicG3TT!
      output << "function dynamicG3TT!(T::Vector{Float64}," << endl
             << "                      y::Vector{Float64}, x::Matrix{Float64}, "
             << "params::Vector{Float64}, steady_state::Vector{Float64}, it_::Int)" << endl
             << "    dynamicG2TT!(T, y, x, params, steady_state, it_)" << endl
             << third_derivatives_tt_output.str()
             << "    return nothing" << endl
             << "end" << endl << endl;

      // dynamicG3!
      int ncols = hessianColsNbr * dynJacobianColsNbr;
      output << "function dynamicG3!(T::Vector{Float64}, g3::Matrix{Float64}," << endl
             << "                    y::Vector{Float64}, x::Matrix{Float64}, "
             << "params::Vector{Float64}, steady_state::Vector{Float64}, it_::Int, T_flag::Bool)" << endl
             << "    @assert length(T) >= "
             << temporary_terms_mlv.size() + temporary_terms_derivatives[0].size() + temporary_terms_derivatives[1].size() + temporary_terms_derivatives[2].size() + temporary_terms_derivatives[3].size() << endl
             << "    @assert size(g3) == (" << nrows << ", " << ncols << ")" << endl
             << "    @assert length(y)+size(x, 2) == " << dynJacobianColsNbr << endl
             << "    @assert length(params) == " << symbol_table.param_nbr() << endl
             << "    if T_flag" << endl
             << "      dynamicG3TT!(T, y, x, params, steady_state, it_)" << endl
             << "    end" << endl
             << "    fill!(g3, 0.0)" << endl
             << third_derivatives_output.str()
             << "    return nothing" << endl
             << "end" << endl << endl;

      // dynamic!
      output << "function dynamic!(T::Vector{Float64}, residual::Vector{Float64}," << endl
             << "                  y::Vector{Float64}, x::Matrix{Float64}, "
             << "params::Vector{Float64}, steady_state::Vector{Float64}, it_::Int)" << endl
             << "    dynamicResid!(T, residual, y, x, params, steady_state, it_, true)" << endl
             << "    return nothing" << endl
             << "end" << endl
             << endl
             << "function dynamic!(T::Vector{Float64}, residual::Vector{Float64}, g1::Matrix{Float64}," << endl
             << "                  y::Vector{Float64}, x::Matrix{Float64}, "
             << "params::Vector{Float64}, steady_state::Vector{Float64}, it_::Int)" << endl
             << "    dynamicG1!(T, g1, y, x, params, steady_state, it_, true)" << endl
             << "    dynamicResid!(T, residual, y, x, params, steady_state, it_, false)" << endl
             << "    return nothing" << endl
             << "end" << endl
             << endl
             << "function dynamic!(T::Vector{Float64}, residual::Vector{Float64}, g1::Matrix{Float64}, g2::Matrix{Float64}," << endl
             << "                  y::Vector{Float64}, x::Matrix{Float64}, "
             << "params::Vector{Float64}, steady_state::Vector{Float64}, it_::Int)" << endl
             << "    dynamicG2!(T, g2, y, x, params, steady_state, it_, true)" << endl
             << "    dynamicG1!(T, g1, y, x, params, steady_state, it_, false)" << endl
             << "    dynamicResid!(T, residual, y, x, params, steady_state, it_, false)" << endl
             << "    return nothing" << endl
             << "end" << endl
             << endl
             << "function dynamic!(T::Vector{Float64}, residual::Vector{Float64}, g1::Matrix{Float64}, g2::Matrix{Float64}, g3::Matrix{Float64}," << endl
             << "                  y::Vector{Float64}, x::Matrix{Float64}, "
             << "params::Vector{Float64}, steady_state::Vector{Float64}, it_::Int)" << endl
             << "    dynamicG3!(T, g3, y, x, params, steady_state, it_, true)" << endl
             << "    dynamicG2!(T, g2, y, x, params, steady_state, it_, false)" << endl
             << "    dynamicG1!(T, g1, y, x, params, steady_state, it_, false)" << endl
             << "    dynamicResid!(T, residual, y, x, params, steady_state, it_, false)" << endl
             << "    return nothing" << endl
             << "end" << endl
             << "end" << endl;
      output.close();
    }
}

void
DynamicModel::writeOutput(ostream &output, const string &basename, bool block_decomposition, bool linear_decomposition, bool byte_code, bool use_dll, int order, bool estimation_present, bool compute_xrefs, bool julia) const
{
  /* Writing initialisation for M_.lead_lag_incidence matrix
     M_.lead_lag_incidence is a matrix with as many columns as there are
     endogenous variables and as many rows as there are periods in the
     models (nbr of rows = M_.max_lag+M_.max_lead+1)

     The matrix elements are equal to zero if a variable isn't present in the
     model at a given period.
  */

  string modstruct;
  string outstruct;
  if (julia)
    {
      modstruct = "model_.";
      outstruct = "oo_.";
    }
  else
    {
      modstruct = "M_.";
      outstruct = "oo_.";
    }

  output << modstruct << "orig_maximum_endo_lag = " << max_endo_lag_orig << ";" << endl
         << modstruct << "orig_maximum_endo_lead = " << max_endo_lead_orig << ";" << endl
         << modstruct << "orig_maximum_exo_lag = " << max_exo_lag_orig << ";" << endl
         << modstruct << "orig_maximum_exo_lead = " << max_exo_lead_orig << ";" << endl
         << modstruct << "orig_maximum_exo_det_lag = " << max_exo_det_lag_orig << ";" << endl
         << modstruct << "orig_maximum_exo_det_lead = " << max_exo_det_lead_orig << ";" << endl
         << modstruct << "orig_maximum_lag = " << max_lag_orig << ";" << endl
         << modstruct << "orig_maximum_lead = " << max_lead_orig << ";" << endl
         << modstruct << "lead_lag_incidence = [";
  // Loop on endogenous variables
  int nstatic = 0,
    nfwrd   = 0,
    npred   = 0,
    nboth   = 0;
  for (int endoID = 0; endoID < symbol_table.endo_nbr(); endoID++)
    {
      output << endl;
      int sstatic = 1,
        sfwrd   = 0,
        spred   = 0,
        sboth   = 0;
      // Loop on periods
      for (int lag = -max_endo_lag; lag <= max_endo_lead; lag++)
        {
          // Print variableID if exists with current period, otherwise print 0
          try
            {
              int varID = getDerivID(symbol_table.getID(SymbolType::endogenous, endoID), lag);
              output << " " << getDynJacobianCol(varID) + 1;
              if (lag == -1)
                {
                  sstatic = 0;
                  spred = 1;
                }
              else if (lag == 1)
                {
                  if (spred == 1)
                    {
                      sboth = 1;
                      spred = 0;
                    }
                  else
                    {
                      sstatic = 0;
                      sfwrd = 1;
                    }
                }
            }
          catch (UnknownDerivIDException &e)
            {
              output << " 0";
            }
        }
      nstatic += sstatic;
      nfwrd   += sfwrd;
      npred   += spred;
      nboth   += sboth;
      output << ";";
    }
  output << "]';" << endl;
  output << modstruct << "nstatic = " << nstatic << ";" << endl
         << modstruct << "nfwrd   = " << nfwrd   << ";" << endl
         << modstruct << "npred   = " << npred   << ";" << endl
         << modstruct << "nboth   = " << nboth   << ";" << endl
         << modstruct << "nsfwrd   = " << nfwrd+nboth   << ";" << endl
         << modstruct << "nspred   = " << npred+nboth   << ";" << endl
         << modstruct << "ndynamic   = " << npred+nboth+nfwrd << ";" << endl;
  if (!julia)
    output << modstruct << "dynamic_tmp_nbr = zeros(4,1); % Number of temporaries used for the dynamic model" << endl
           << modstruct << "dynamic_tmp_nbr(1) = " << temporary_terms_derivatives[0].size() << "; % Number of temporaries used for the evaluation of the residuals" << endl
           << modstruct << "dynamic_tmp_nbr(2) = " << temporary_terms_derivatives[1].size() << "; % Number of temporaries used for the evaluation of g1 (jacobian)" << endl
           << modstruct << "dynamic_tmp_nbr(3) = " << temporary_terms_derivatives[2].size() << "; % Number of temporaries used for the evaluation of g2 (hessian)" << endl
           << modstruct << "dynamic_tmp_nbr(4) = " << temporary_terms_derivatives[3].size() << "; % Number of temporaries used for the evaluation of g3 (third order derivatives)" << endl;

  // Write equation tags
  if (julia)
    {
      output << modstruct << "equation_tags = [" << endl;
      for (const auto & equation_tag : equation_tags)
        output << "                       EquationTag("
               << equation_tag.first + 1 << " , \""
               << equation_tag.second.first << "\" , \""
               << equation_tag.second.second << "\")" << endl;
      output << "                      ]" << endl;
    }
  else
    {
      output << modstruct << "equations_tags = {" << endl;
      for (const auto & equation_tag : equation_tags)
        output << "  " << equation_tag.first + 1 << " , '"
               << equation_tag.second.first << "' , '"
               << equation_tag.second.second << "' ;" << endl;
      output << "};" << endl;
    }

  /* Say if static and dynamic models differ (because of [static] and [dynamic]
     equation tags) */
  output << modstruct << "static_and_dynamic_models_differ = "
         << (static_only_equations.size() > 0 ?
             (julia ? "true" : "1") :
             (julia ? "false" : "0"))
         << ";" << endl;

  vector<int> state_var;
  for (int endoID = 0; endoID < symbol_table.endo_nbr(); endoID++)
    // Loop on periods
    for (int lag = -max_endo_lag; lag < 0; lag++)
      try
        {
          getDerivID(symbol_table.getID(SymbolType::endogenous, variable_reordered[endoID]), lag);
          if (lag < 0 && find(state_var.begin(), state_var.end(), variable_reordered[endoID]+1) == state_var.end())
            state_var.push_back(variable_reordered[endoID]+1);
        }
      catch (UnknownDerivIDException &e)
        {
        }

  //In case of sparse model, writes the block_decomposition structure of the model
  if (block_decomposition || linear_decomposition)
    {
      vector<int> state_equ;
      int count_lead_lag_incidence = 0;
      int max_lead, max_lag, max_lag_endo, max_lead_endo, max_lag_exo, max_lead_exo, max_lag_exo_det, max_lead_exo_det;
      unsigned int nb_blocks = getNbBlocks();
      for (unsigned int block = 0; block < nb_blocks; block++)
        {
          //For a block composed of a single equation determines wether we have to evaluate or to solve the equation
          count_lead_lag_incidence = 0;
          BlockSimulationType simulation_type = getBlockSimulationType(block);
          int block_size = getBlockSize(block);
          max_lag  = max_leadlag_block[block].first;
          max_lead = max_leadlag_block[block].second;
          max_lag_endo = endo_max_leadlag_block[block].first;
          max_lead_endo = endo_max_leadlag_block[block].second;
          max_lag_exo = exo_max_leadlag_block[block].first;
          max_lead_exo = exo_max_leadlag_block[block].second;
          max_lag_exo_det = exo_det_max_leadlag_block[block].first;
          max_lead_exo_det = exo_det_max_leadlag_block[block].second;
          ostringstream tmp_s, tmp_s_eq;
          tmp_s.str("");
          tmp_s_eq.str("");
          for (int i = 0; i < block_size; i++)
            {
              tmp_s << " " << getBlockVariableID(block, i)+1;
              tmp_s_eq << " " << getBlockEquationID(block, i)+1;
            }
          set<int> exogenous;
          for (const auto & it : exo_block[block])
            for (int it1 : it.second)
              exogenous.insert(it1);
          set<int> exogenous_det;
          for (const auto & it : exo_det_block[block])
            for (int it1 : it.second)
              exogenous_det.insert(it1);
          set<int> other_endogenous;
          for (const auto & it : other_endo_block[block])
            for (int it1 : it.second)
              other_endogenous.insert(it1);
          output << "block_structure.block(" << block+1 << ").Simulation_Type = " << simulation_type << ";" << endl
                 << "block_structure.block(" << block+1 << ").maximum_lag = " << max_lag << ";" << endl
                 << "block_structure.block(" << block+1 << ").maximum_lead = " << max_lead << ";" << endl
                 << "block_structure.block(" << block+1 << ").maximum_endo_lag = " << max_lag_endo << ";" << endl
                 << "block_structure.block(" << block+1 << ").maximum_endo_lead = " << max_lead_endo << ";" << endl
                 << "block_structure.block(" << block+1 << ").maximum_exo_lag = " << max_lag_exo << ";" << endl
                 << "block_structure.block(" << block+1 << ").maximum_exo_lead = " << max_lead_exo << ";" << endl
                 << "block_structure.block(" << block+1 << ").maximum_exo_det_lag = " << max_lag_exo_det << ";" << endl
                 << "block_structure.block(" << block+1 << ").maximum_exo_det_lead = " << max_lead_exo_det << ";" << endl
                 << "block_structure.block(" << block+1 << ").endo_nbr = " << block_size << ";" << endl
                 << "block_structure.block(" << block+1 << ").mfs = " << getBlockMfs(block) << ";" << endl
                 << "block_structure.block(" << block+1 << ").equation = [" << tmp_s_eq.str() << "];" << endl
                 << "block_structure.block(" << block+1 << ").variable = [" << tmp_s.str() << "];" << endl
                 << "block_structure.block(" << block+1 << ").exo_nbr = " << getBlockExoSize(block) << ";" << endl
                 << "block_structure.block(" << block+1 << ").exogenous = [";
          int i = 0;
          for (int exogenou : exogenous)
            if (exogenou >= 0)
              {
                output << " " << exogenou+1;
                i++;
              }
          output << "];" << endl
                 << "block_structure.block(" << block+1 << ").exogenous_det = [";
          i = 0;
          for (int it_exogenous_det : exogenous_det)
            if (it_exogenous_det >= 0)
              {
                output << " " << it_exogenous_det+1;
                i++;
              }
          output << "];" << endl
                 << "block_structure.block(" << block+1 << ").exo_det_nbr = " << i << ";" << endl
                 << "block_structure.block(" << block+1 << ").other_endogenous = [";
          i = 0;
          for (int other_endogenou : other_endogenous)
            if (other_endogenou >= 0)
              {
                output << " " << other_endogenou+1;
                i++;
              }
          output << "];" << endl
                 << "block_structure.block(" << block+1 << ").other_endogenous_block = [";
          i = 0;
          for (int other_endogenou : other_endogenous)
            if (other_endogenou >= 0)
              {
                bool OK = true;
                unsigned int j;
                for (j = 0; j < block && OK; j++)
                  for (unsigned int k = 0; k < getBlockSize(j) && OK; k++)
                    {
                      //printf("*it_other_endogenous=%d, getBlockVariableID(%d, %d)=%d\n",*it_other_endogenous, j, k, getBlockVariableID(j, k));
                      OK = other_endogenou != getBlockVariableID(j, k);
                    }
                if (!OK)
                  output << " " << j;
                i++;
              }
          output << "];" << endl;

          //vector<int> inter_state_var;
          output << "block_structure.block(" << block+1 << ").tm1 = zeros(" << i << ", " << state_var.size() << ");" << endl;
          int count_other_endogenous = 1;
          for (int other_endogenou : other_endogenous)
            {
              for (auto it = state_var.begin(); it != state_var.end(); ++it)
                {
                  //cout << "block = " << block+1 << " state_var = " << *it << " it_other_endogenous=" << *it_other_endogenous + 1 << "\n";
                  if (*it == other_endogenou + 1)
                    {
                      output << "block_structure.block(" << block+1 << ").tm1("
                             << count_other_endogenous << ", "
                             << it - state_var.begin()+1 << ") = 1;" << endl;
                      /*output << "block_structure.block(" << block+1 << ").tm1("
                        << it - state_var.begin()+1 << ", "
                        << count_other_endogenous << ") = 1;\n";*/
                      //cout << "=>\n";
                    }
                }
              count_other_endogenous++;
            }

          output << "block_structure.block(" << block+1 << ").other_endo_nbr = " << i << ";" << endl;

          tmp_s.str("");
          count_lead_lag_incidence = 0;
          dynamic_jacob_map_t reordered_dynamic_jacobian;
          for (const auto & it : blocks_derivatives[block])
            reordered_dynamic_jacobian[{ get<2>(it), get<1>(it), get<0>(it) }] = get<3>(it);
          output << "block_structure.block(" << block+1 << ").lead_lag_incidence = [];" << endl;
          int last_var = -1;
          vector<int> local_state_var;
          vector<int> local_stat_var;
          int n_static = 0, n_backward = 0, n_forward = 0, n_mixed = 0;
          for (int lag = -1; lag < 1+1; lag++)
            {
              last_var = -1;
              for (const auto &it : reordered_dynamic_jacobian)
                {
                  if (lag == get<0>(it.first) && last_var != get<1>(it.first))
                    {
                      if (lag == -1)
                        {
                          local_state_var.push_back(getBlockVariableID(block, get<1>(it.first))+1);
                          n_backward++;
                        }
                      else if (lag == 0)
                        {
                          if (find(local_state_var.begin(), local_state_var.end(), getBlockVariableID(block, get<1>(it.first))+1) == local_state_var.end())
                            {
                              local_stat_var.push_back(getBlockVariableID(block, get<1>(it.first))+1);
                              n_static++;
                            }
                        }
                      else
                        {
                          if (find(local_state_var.begin(), local_state_var.end(), getBlockVariableID(block, get<1>(it.first))+1) != local_state_var.end())
                            {
                              n_backward--;
                              n_mixed++;
                            }
                          else
                            {
                              if (find(local_stat_var.begin(), local_stat_var.end(), getBlockVariableID(block, get<1>(it.first))+1) != local_stat_var.end())
                                n_static--;
                              n_forward++;
                            }
                        }
                      count_lead_lag_incidence++;
                      for (int i = last_var; i < get<1>(it.first)-1; i++)
                        tmp_s << " 0";
                      if (tmp_s.str().length())
                        tmp_s << " ";
                      tmp_s << count_lead_lag_incidence;
                      last_var = get<1>(it.first);
                    }
                }
              for (int i = last_var + 1; i < block_size; i++)
                tmp_s << " 0";
              output << "block_structure.block(" << block+1 << ").lead_lag_incidence = [ block_structure.block(" << block+1 << ").lead_lag_incidence; " << tmp_s.str() << "]; %lag = " << lag << endl;
              tmp_s.str("");
            }
          vector<int> inter_state_var;
          for (auto it_l = local_state_var.begin(); it_l != local_state_var.end(); ++it_l)
            for (auto it = state_var.begin(); it != state_var.end(); ++it)
              if (*it == *it_l)
                inter_state_var.push_back(it - state_var.begin()+1);
          output << "block_structure.block(" << block+1 << ").sorted_col_dr_ghx = [";
          for (int it : inter_state_var)
            output << it << " ";
          output << "];" << endl;
          count_lead_lag_incidence = 0;
          output << "block_structure.block(" << block+1 << ").lead_lag_incidence_other = [];" << endl;
          for (int lag = -1; lag <= 1; lag++)
            {
              tmp_s.str("");
              for (int other_endogenou : other_endogenous)
                {
                  bool done = false;
                  for (int i = 0; i < block_size; i++)
                    {
                      unsigned int eq = getBlockEquationID(block, i);
                      auto it = derivative_other_endo[block].find({ lag, eq, other_endogenou });
                      if (it != derivative_other_endo[block].end())
                        {
                          count_lead_lag_incidence++;
                          tmp_s << " " << count_lead_lag_incidence;
                          done = true;
                          break;
                        }
                    }
                  if (!done)
                    tmp_s << " 0";
                }
              output << "block_structure.block(" << block+1 << ").lead_lag_incidence_other = [ block_structure.block(" << block+1 << ").lead_lag_incidence_other; " << tmp_s.str() << "]; %lag = " << lag << endl;
            }
          output << "block_structure.block(" << block+1 << ").n_static = " << n_static << ";" << endl
                 << "block_structure.block(" << block+1 << ").n_forward = " << n_forward << ";" << endl
                 << "block_structure.block(" << block+1 << ").n_backward = " << n_backward << ";" << endl
                 << "block_structure.block(" << block+1 << ").n_mixed = " << n_mixed << ";" << endl;
        }
      output << modstruct << "block_structure.block = block_structure.block;" << endl;
      string cst_s;
      int nb_endo = symbol_table.endo_nbr();
      output << modstruct << "block_structure.variable_reordered = [";
      for (int i = 0; i < nb_endo; i++)
        output << " " << variable_reordered[i]+1;
      output << "];" << endl;
      output << modstruct << "block_structure.equation_reordered = [";
      for (int i = 0; i < nb_endo; i++)
        output << " " << equation_reordered[i]+1;
      output << "];" << endl;
      vector<int> variable_inv_reordered(nb_endo);

      for (int i = 0; i < nb_endo; i++)
        variable_inv_reordered[variable_reordered[i]] = i;

      for (int it : state_var)
        state_equ.push_back(equation_reordered[variable_inv_reordered[it - 1]]+1);

      map<tuple<int, int, int>,  int> lag_row_incidence;
      for (const auto & first_derivative : derivatives[1])
        {
          int deriv_id = first_derivative.first[1];
          if (getTypeByDerivID(deriv_id) == SymbolType::endogenous)
            {
              int eq = first_derivative.first[0];
              int symb = getSymbIDByDerivID(deriv_id);
              int var = symbol_table.getTypeSpecificID(symb);
              int lag = getLagByDerivID(deriv_id);
              lag_row_incidence[{ lag, eq, var }] = 1;
            }
        }
      int prev_lag = -1000000;
      for (const auto &it : lag_row_incidence)
        {
          if (prev_lag != get<0>(it.first))
            {
              if (prev_lag != -1000000)
                output << "];" << endl;
              prev_lag = get<0>(it.first);
              output << modstruct << "block_structure.incidence(" << max_endo_lag+get<0>(it.first)+1 << ").lead_lag = " << prev_lag << ";" << endl
                     << modstruct << "block_structure.incidence(" << max_endo_lag+get<0>(it.first)+1 << ").sparse_IM = [";
            }
          output << get<1>(it.first)+1 << " " << get<2>(it.first)+1 << ";" << endl;
        }
      output << "];" << endl;
      if (estimation_present)
        {
          ofstream KF_index_file;
          boost::filesystem::create_directories(basename + "/model/bytecode");
          string main_name = basename + "/model/bytecode/kfi";
          KF_index_file.open(main_name, ios::out | ios::binary | ios::ate);
          int n_obs = symbol_table.observedVariablesNbr();
          int n_state = state_var.size();
          for (int it : state_var)
            if (symbol_table.isObservedVariable(symbol_table.getID(SymbolType::endogenous, it-1)))
              n_obs--;

          int n = n_obs + n_state;
          output << modstruct << "nobs_non_statevar = " << n_obs << ";" << endl;
          int nb_diag = 0;
          //map<pair<int,int>, int>::const_iterator  row_state_var_incidence_it = row_state_var_incidence.begin();

          vector<int> i_nz_state_var(n);
          for (int i = 0; i < n_obs; i++)
            i_nz_state_var[i] = n;
          unsigned int lp = n_obs;

          for (unsigned int block = 0; block < nb_blocks; block++)
            {
              int block_size = getBlockSize(block);
              int nze = 0;

              for (int i = 0; i < block_size; i++)
                {
                  int var = getBlockVariableID(block, i);
                  auto it_state_var = find(state_var.begin(), state_var.end(), var+1);
                  if (it_state_var != state_var.end())
                    nze++;
                }
              if (block == 0)
                {
                  set<pair<int, int>> row_state_var_incidence;
                  for (const auto &it : blocks_derivatives[block])
                    {
                      auto it_state_var = find(state_var.begin(), state_var.end(), getBlockVariableID(block, get<1>(it))+1);
                      if (it_state_var != state_var.end())
                        {
                          auto it_state_equ = find(state_equ.begin(), state_equ.end(), getBlockEquationID(block, get<0>(it))+1);
                          if (it_state_equ != state_equ.end())
                            row_state_var_incidence.emplace(it_state_equ - state_equ.begin(), it_state_var - state_var.begin());
                        }

                    }
                  /*tmp_block_endo_derivative[make_pair(it->second.first, make_pair(it->first.second, it->first.first))] = it->second.second;
                    if (block == 0)
                    {

                    vector<int>::const_iterator it_state_equ = find(state_equ.begin(), state_equ.end(), getBlockEquationID(block, i)+1);
                    if (it_state_equ != state_equ.end())
                    {
                    cout << "row_state_var_incidence[make_pair([" << *it_state_equ << "] " << it_state_equ - state_equ.begin() << ", [" << *it_state_var << "] " << it_state_var - state_var.begin() << ")] =  1;\n";
                    row_state_var_incidence.insert(make_pair(it_state_equ - state_equ.begin(), it_state_var - state_var.begin()));
                    }
                    }*/
                  auto  row_state_var_incidence_it = row_state_var_incidence.begin();
                  bool diag = true;
                  int nb_diag_r = 0;
                  while (row_state_var_incidence_it != row_state_var_incidence.end() && diag)
                    {
                      diag = (row_state_var_incidence_it->first == row_state_var_incidence_it->second);
                      if (diag)
                        {
                          int equ = row_state_var_incidence_it->first;
                          row_state_var_incidence_it++;
                          if (equ != row_state_var_incidence_it->first)
                            nb_diag_r++;
                        }

                    }
                  set<pair<int, int>>  col_state_var_incidence;
                  for (const auto & row_state_var_incidence_it : row_state_var_incidence)
                    col_state_var_incidence.emplace(row_state_var_incidence_it.second, row_state_var_incidence_it.first);
                  auto  col_state_var_incidence_it = col_state_var_incidence.begin();
                  diag = true;
                  int nb_diag_c = 0;
                  while (col_state_var_incidence_it != col_state_var_incidence.end() && diag)
                    {
                      diag = (col_state_var_incidence_it->first == col_state_var_incidence_it->second);
                      if (diag)
                        {
                          int var = col_state_var_incidence_it->first;
                          col_state_var_incidence_it++;
                          if (var != col_state_var_incidence_it->first)
                            nb_diag_c++;
                        }
                    }
                  nb_diag = min(nb_diag_r, nb_diag_c);
                  row_state_var_incidence.clear();
                  col_state_var_incidence.clear();
                }
              for (int i = 0; i < nze; i++)
                i_nz_state_var[lp + i] = lp + nze;
              lp += nze;
            }
          output << modstruct << "nz_state_var = [";
          for (unsigned int i = 0; i < lp; i++)
            output << i_nz_state_var[i] << " ";
          output << "];" << endl;
          output << modstruct << "n_diag = " << nb_diag << ";" << endl;
          KF_index_file.write(reinterpret_cast<char *>(&nb_diag), sizeof(nb_diag));

          using index_KF = pair<int, pair<int, int >>;
          vector<index_KF> v_index_KF;
          for (int i = 0; i < n; i++)
            //int i = 0;
            for (int j = n_obs; j < n; j++)
              {
                int j1 = j - n_obs;
                int j1_n_state = j1 * n_state - n_obs;
                if ((i < n_obs) || (i >= nb_diag + n_obs) || (j1 >= nb_diag))
                  for (int k = n_obs; k < i_nz_state_var[i]; k++)
                    {
                      v_index_KF.emplace_back(i + j1 * n, make_pair(i + k * n, k + j1_n_state));
                    }
              }
          int size_v_index_KF = v_index_KF.size();

          KF_index_file.write(reinterpret_cast<char *>(&size_v_index_KF), sizeof(size_v_index_KF));
          for (auto & it : v_index_KF)
            KF_index_file.write(reinterpret_cast<char *>(&it), sizeof(index_KF));

          vector<index_KF> v_index_KF_2;
          int n_n_obs = n * n_obs;
          for (int i = 0; i < n; i++)
            //i = 0;
            for (int j = i; j < n; j++)
              {
                if ((i < n_obs) || (i >= nb_diag + n_obs) || (j < n_obs) || (j >= nb_diag + n_obs))
                  for (int k = n_obs; k < i_nz_state_var[j]; k++)
                    {
                      int k_n = k * n;
                      v_index_KF_2.emplace_back(i * n + j,  make_pair(i + k_n - n_n_obs, j + k_n));
                    }
              }
          int size_v_index_KF_2 = v_index_KF_2.size();

          KF_index_file.write(reinterpret_cast<char *>(&size_v_index_KF_2), sizeof(size_v_index_KF_2));
          for (auto & it : v_index_KF_2)
            KF_index_file.write(reinterpret_cast<char *>(&it), sizeof(index_KF));
          KF_index_file.close();
        }
    }

  output << modstruct << "state_var = [";
  for (int it : state_var)
    output << it << (julia ? "," : " ");
  output << "];" << endl;

  // Writing initialization for some other variables
  if (!julia)
    output << modstruct << "exo_names_orig_ord = [1:" << symbol_table.exo_nbr() << "];" << endl;
  else
    output << modstruct << "exo_names_orig_ord = collect(1:" << symbol_table.exo_nbr() << ");" << endl;

  output << modstruct << "maximum_lag = " << max_lag << ";" << endl
         << modstruct << "maximum_lead = " << max_lead << ";" << endl;

  output << modstruct << "maximum_endo_lag = " << max_endo_lag << ";" << endl
         << modstruct << "maximum_endo_lead = " << max_endo_lead << ";" << endl
         << outstruct << "steady_state = zeros(" << symbol_table.endo_nbr() << (julia ? ")" : ", 1);") << endl;

  output << modstruct << "maximum_exo_lag = " << max_exo_lag << ";" << endl
         << modstruct << "maximum_exo_lead = " << max_exo_lead << ";" << endl
         << outstruct << "exo_steady_state = zeros(" << symbol_table.exo_nbr() <<  (julia ? ")" : ", 1);")   << endl;

  if (symbol_table.exo_det_nbr())
    {
      output << modstruct << "maximum_exo_det_lag = " << max_exo_det_lag << ";" << endl
             << modstruct << "maximum_exo_det_lead = " << max_exo_det_lead << ";" << endl
             << outstruct << "exo_det_steady_state = zeros(" << symbol_table.exo_det_nbr() << (julia ? ")" : ", 1);") << endl;
    }

  output << modstruct << "params = " << (julia ? "fill(NaN, " : "NaN(")
         << symbol_table.param_nbr() << (julia ? ")" : ", 1);") << endl;

  if (compute_xrefs)
    writeXrefs(output);

  // Write number of non-zero derivatives
  // Use -1 if the derivatives have not been computed
  output << modstruct << (julia ? "nnzderivatives" : "NNZDerivatives")
         << " = [" << NNZDerivatives[1] << "; ";
  if (order > 1)
    output << NNZDerivatives[2] << "; ";
  else
    output << "-1; ";

  if (order > 2)
    output << NNZDerivatives[3];
  else
    output << "-1";
  output << "];" << endl;

  // Write PacExpectationInfo
  for (auto it : pac_expectation_info)
    it->ExprNode::writeOutput(output, ExprNodeOutputType::matlabDynamicModel);
}

map<tuple<int, int, int>, expr_t>
DynamicModel::collect_first_order_derivatives_endogenous()
{
  map<tuple<int, int, int>, expr_t> endo_derivatives;
  for (auto & first_derivative : derivatives[1])
    {
      if (getTypeByDerivID(first_derivative.first[1]) == SymbolType::endogenous)
        {
          int eq = first_derivative.first[0];
          int var = symbol_table.getTypeSpecificID(getSymbIDByDerivID(first_derivative.first[1]));
          int lag = getLagByDerivID(first_derivative.first[1]);
          endo_derivatives[{ eq, var, lag }] = first_derivative.second;
        }
    }
  return endo_derivatives;
}

void
DynamicModel::runTrendTest(const eval_context_t &eval_context)
{
  computeDerivIDs();
  testTrendDerivativesEqualToZero(eval_context);
}

void
DynamicModel::updateVarAndTrendModel() const
{
  for (int i = 0; i < 2; i++)
    {
      map<string, vector<int>> eqnums, trend_eqnums;
      if (i == 0)
        eqnums = var_model_table.getEqNums();
      else if (i == 1)
        {
          eqnums = trend_component_model_table.getEqNums();
          trend_eqnums = trend_component_model_table.getTargetEqNums();
        }

      map<string, vector<int>> trend_varr;
      map<string, vector<set<pair<int, int>>>> rhsr;
      for (const auto & it : eqnums)
        {
          vector<int> lhs, trend_var, trend_lhs;
          vector<set<pair<int, int>>> rhs;

          if (i == 1)
            {
              lhs = trend_component_model_table.getLhs(it.first);
              for (auto teqn : trend_eqnums.at(it.first))
                {
                  int eqnidx = 0;
                  for (auto eqn : it.second)
                    {
                      if (eqn == teqn)
                        trend_lhs.push_back(lhs[eqnidx]);
                      eqnidx++;
                    }
                }
            }

          int lhs_idx = 0;
          for (auto eqn : it.second)
            {
              set<pair<int, int>> rhs_set;
              equations[eqn]->get_arg2()->collectDynamicVariables(SymbolType::endogenous, rhs_set);
              rhs.push_back(rhs_set);

              if (i == 1)
                {
                  int lhs_symb_id = lhs[lhs_idx++];
                  if (symbol_table.isAuxiliaryVariable(lhs_symb_id))
                    try
                      {
                        lhs_symb_id = symbol_table.getOrigSymbIdForAuxVar(lhs_symb_id);
                      }
                    catch (...)
                      {
                      }
                  int trend_var_symb_id = equations[eqn]->get_arg2()->findTargetVariable(lhs_symb_id);
                  if (trend_var_symb_id >= 0)
                    {
                      if (symbol_table.isAuxiliaryVariable(trend_var_symb_id))
                        try
                          {
                            trend_var_symb_id = symbol_table.getOrigSymbIdForAuxVar(trend_var_symb_id);
                          }
                        catch (...)
                          {
                          }
                      if (find(trend_lhs.begin(), trend_lhs.end(), trend_var_symb_id) == trend_lhs.end())
                        {
                          cerr << "ERROR: trend found in trend_component equation #" << eqn << " ("
                               << symbol_table.getName(trend_var_symb_id) << ") does not correspond to a trend equation" << endl;
                          exit(EXIT_FAILURE);
                        }
                    }
                  trend_var.push_back(trend_var_symb_id);
                }
            }

          rhsr[it.first] = rhs;
          if (i == 1)
            trend_varr[it.first] = trend_var;
        }

      if (i == 0)
        var_model_table.setRhs(rhsr);
      else if (i == 1)
        {
          trend_component_model_table.setRhs(rhsr);
          trend_component_model_table.setTargetVar(trend_varr);
        }
    }
}

void
DynamicModel::fillVarModelTable() const
{
  map<string, vector<int>> eqnums, lhsr;
  map<string, vector<expr_t>> lhs_expr_tr;
  map<string, vector<set<pair<int, int>>>> rhsr;
  map<string, vector<string>> eqtags = var_model_table.getEqTags();

  for (const auto & it : eqtags)
    {
      vector<int> eqnumber, lhs;
      vector<expr_t> lhs_expr_t;
      vector<set<pair<int, int>>> rhs;

      for (const auto & eqtag : it.second)
        {
          int eqn = -1;
          set<pair<int, int>> lhs_set, lhs_tmp_set, rhs_set;
          for (const auto & equation_tag : equation_tags)
            if (equation_tag.second.first == "name"
                && equation_tag.second.second == eqtag)
              {
                eqn = equation_tag.first;
                break;
              }

          if (eqn == -1)
            {
              cerr << "ERROR: equation tag '" << eqtag << "' not found" << endl;
              exit(EXIT_FAILURE);
            }

          equations[eqn]->get_arg1()->collectDynamicVariables(SymbolType::endogenous, lhs_set);
          equations[eqn]->get_arg1()->collectDynamicVariables(SymbolType::exogenous, lhs_tmp_set);
          equations[eqn]->get_arg1()->collectDynamicVariables(SymbolType::parameter, lhs_tmp_set);

          if (lhs_set.size() != 1 || !lhs_tmp_set.empty())
            {
              cerr << "ERROR: in Equation " << eqtag
                   << ". A VAR may only have one endogenous variable on the LHS. " << endl;
              exit(EXIT_FAILURE);
            }

          auto itlhs = lhs_set.begin();
          if (itlhs->second != 0)
            {
              cerr << "ERROR: in Equation " << eqtag
                   << ". The variable on the LHS of a VAR may not appear with a lead or a lag. "
                   << endl;
              exit(EXIT_FAILURE);
            }

          eqnumber.push_back(eqn);
          lhs.push_back(itlhs->first);
          lhs_set.clear();
          set<expr_t> lhs_expr_t_set;
          equations[eqn]->get_arg1()->collectVARLHSVariable(lhs_expr_t_set);
          lhs_expr_t.push_back(*(lhs_expr_t_set.begin()));

          equations[eqn]->get_arg2()->collectDynamicVariables(SymbolType::endogenous, rhs_set);
          for (const auto & itrhs : rhs_set)
            if (itrhs.second > 0)
              {
                cerr << "ERROR: in Equation " << eqtag
                     << ". A VAR may not have leaded or contemporaneous variables on the RHS. " << endl;
                exit(EXIT_FAILURE);
              }
          rhs.push_back(rhs_set);
        }
      eqnums[it.first] = eqnumber;
      lhsr[it.first] = lhs;
      lhs_expr_tr[it.first] = lhs_expr_t;
      rhsr[it.first] = rhs;
    }
  var_model_table.setEqNums(eqnums);
  var_model_table.setLhs(lhsr);
  var_model_table.setRhs(rhsr);
  var_model_table.setLhsExprT(lhs_expr_tr);

  // Fill AR Matrix
  map<string, map<tuple<int, int, int>, expr_t>> ARr;
  fillAutoregressiveMatrix(ARr, false);
  var_model_table.setAR(ARr);
}

void
DynamicModel::fillVarModelTableFromOrigModel(StaticModel &static_model) const
{
  map<string, vector<int>> lags, orig_diff_var;
  map<string, vector<bool>> diff;
  for (const auto & it : var_model_table.getEqNums())
    {
      set<expr_t> lhs;
      vector<int> orig_diff_var_vec;
      vector<bool> diff_vec;
      for (auto eqn : it.second)
        {
          // ensure no leads in equations
          if (equations[eqn]->get_arg2()->VarMinLag() <= 0)
            {
              cerr << "ERROR in VAR model Equation (#" << eqn << "). "
                   << "Leaded exogenous variables "
                   << "and leaded or contemporaneous endogenous variables not allowed in VAR"
                   << endl;
              exit(EXIT_FAILURE);
            }

          // save lhs variables
          equations[eqn]->get_arg1()->collectVARLHSVariable(lhs);

          equations[eqn]->get_arg1()->countDiffs() > 0 ?
            diff_vec.push_back(true) : diff_vec.push_back(false);
          if (diff_vec.back())
            {
              set<pair<int, int>> diff_set;
              equations[eqn]->get_arg1()->collectDynamicVariables(SymbolType::endogenous, diff_set);

              if (diff_set.size() != 1)
                {
                  cerr << "ERROR: problem getting variable for LHS diff operator in equation "
                       << eqn << endl;
                  exit(EXIT_FAILURE);
                }
              orig_diff_var_vec.push_back(diff_set.begin()->first);
            }
          else
            orig_diff_var_vec.push_back(-1);

        }

      if (it.second.size() != lhs.size())
        {
          cerr << "ERROR: The LHS variables of the VAR model are not unique" << endl;
          exit(EXIT_FAILURE);
        }

      set<expr_t> lhs_static;
      for(const auto & lh : lhs)
        lhs_static.insert(lh->toStatic(static_model));

      vector<int> max_lag;
      for (auto eqn : it.second)
        max_lag.push_back(equations[eqn]->get_arg2()->VarMaxLag(static_model, lhs_static));
      lags[it.first] = max_lag;
      diff[it.first] = diff_vec;
      orig_diff_var[it.first] = orig_diff_var_vec;
    }
  var_model_table.setDiff(diff);
  var_model_table.setMaxLags(lags);
  var_model_table.setOrigDiffVar(orig_diff_var);
}

void
DynamicModel::fillAutoregressiveMatrix(map<string, map<tuple<int, int, int>, expr_t>> &ARr, bool is_trend_component_model) const
{
  auto eqnums = is_trend_component_model ?
    trend_component_model_table.getNonTargetEqNums() : var_model_table.getEqNums();
  for (const auto & it : eqnums)
    {
      int i = 0;
      map<tuple<int, int, int>, expr_t> AR;
      vector<int> lhs = is_trend_component_model ?
        trend_component_model_table.getLhs(it.first) : var_model_table.getLhs(it.first);
      for (auto eqn : it.second)
        equations[eqn]->get_arg2()->fillAutoregressiveRow(i++, lhs, AR);
      ARr[it.first] = AR;
    }
}

void
DynamicModel::fillTrendComponentModelTable() const
{
  map<string, vector<int>> eqnums, trend_eqnums, lhsr;
  map<string, vector<expr_t>> lhs_expr_tr;
  map<string, vector<set<pair<int, int>>>> rhsr;
  map<string, vector<string>> eqtags = trend_component_model_table.getEqTags();
  map<string, vector<string>> trend_eqtags = trend_component_model_table.getTargetEqTags();
  for (const auto & it : trend_eqtags)
    {
      vector<int> trend_eqnumber;
      for (const auto & eqtag : it.second)
        {
          int eqn = -1;
          for (const auto & equation_tag : equation_tags)
            if (equation_tag.second.first == "name"
                && equation_tag.second.second == eqtag)
              {
                eqn = equation_tag.first;
                break;
              }

          if (eqn == -1)
            {
              cerr << "ERROR: trend equation tag '" << eqtag << "' not found" << endl;
              exit(EXIT_FAILURE);
            }
          trend_eqnumber.push_back(eqn);
        }
      trend_eqnums[it.first] = trend_eqnumber;
    }

  for (const auto & it : eqtags)
    {
      vector<int> eqnumber, lhs;
      vector<expr_t> lhs_expr_t;
      vector<set<pair<int, int>>> rhs;

      for (const auto & eqtag : it.second)
        {
          int eqn = -1;
          set<pair<int, int>> lhs_set, lhs_tmp_set, rhs_set;
          for (const auto & equation_tag : equation_tags)
            if (equation_tag.second.first == "name"
                && equation_tag.second.second == eqtag)
              {
                eqn = equation_tag.first;
                break;
              }

          if (eqn == -1)
            {
              cerr << "ERROR: equation tag '" << eqtag << "' not found" << endl;
              exit(EXIT_FAILURE);
            }

          equations[eqn]->get_arg1()->collectDynamicVariables(SymbolType::endogenous, lhs_set);
          equations[eqn]->get_arg1()->collectDynamicVariables(SymbolType::exogenous, lhs_tmp_set);
          equations[eqn]->get_arg1()->collectDynamicVariables(SymbolType::parameter, lhs_tmp_set);

          if (lhs_set.size() != 1 || !lhs_tmp_set.empty())
            {
              cerr << "ERROR: in Equation " << eqtag
                   << ". A trend component model may only have one endogenous variable on the LHS. " << endl;
              exit(EXIT_FAILURE);
            }

          auto itlhs = lhs_set.begin();
          if (itlhs->second != 0)
            {
              cerr << "ERROR: in Equation " << eqtag
                   << ". The variable on the LHS of a trend component model may not appear with a lead or a lag. "
                   << endl;
              exit(EXIT_FAILURE);
            }

          eqnumber.push_back(eqn);
          lhs.push_back(itlhs->first);
          lhs_set.clear();
          set<expr_t> lhs_expr_t_set;
          equations[eqn]->get_arg1()->collectVARLHSVariable(lhs_expr_t_set);
          lhs_expr_t.push_back(*(lhs_expr_t_set.begin()));

          equations[eqn]->get_arg2()->collectDynamicVariables(SymbolType::endogenous, rhs_set);
          for (const auto & itrhs : rhs_set)
            if (itrhs.second > 0)
              {
                cerr << "ERROR: in Equation " << eqtag
                     << ". A trend component model may not have leaded or contemporaneous variables on the RHS. " << endl;
                exit(EXIT_FAILURE);
              }
          rhs.push_back(rhs_set);
        }
      eqnums[it.first] = eqnumber;
      lhsr[it.first] = lhs;
      lhs_expr_tr[it.first] = lhs_expr_t;
      rhsr[it.first] = rhs;
    }
  trend_component_model_table.setRhs(rhsr);
  trend_component_model_table.setVals(eqnums, trend_eqnums, lhsr, lhs_expr_tr);
}

void
DynamicModel::fillErrorComponentMatrix(map<string, map<tuple<int, int, int>, expr_t>> &ECr,
                                       ExprNode::subst_table_t &diff_subst_table) const
{
  for (const auto & it : trend_component_model_table.getEqNums())
    {
      int i = 0;
      map<tuple<int, int, int>, expr_t> EC;
      vector<int> trend_lhs = trend_component_model_table.getTargetLhs(it.first);
      vector<int> nontrend_eqnums = trend_component_model_table.getNonTargetEqNums(it.first);
      vector<int> undiff_nontrend_lhs = getUndiffLHSForPac(it.first, diff_subst_table);
      vector<int> parsed_undiff_nontrend_lhs;

      for (auto eqn : it.second)
        {
          if (find(nontrend_eqnums.begin(), nontrend_eqnums.end(), eqn) != nontrend_eqnums.end())
            parsed_undiff_nontrend_lhs.push_back(undiff_nontrend_lhs.at(i));
          i++;
        }

      i = 0;
      for (auto eqn : it.second)
        if (find(nontrend_eqnums.begin(), nontrend_eqnums.end(), eqn) != nontrend_eqnums.end())
          equations[eqn]->get_arg2()->fillErrorCorrectionRow(i++, parsed_undiff_nontrend_lhs, trend_lhs, EC);
      ECr[it.first] = EC;
    }
}

void
DynamicModel::fillTrendComponentModelTableFromOrigModel(StaticModel &static_model) const
{
  map<string, vector<int>> lags, orig_diff_var;
  map<string, vector<bool>> diff;
  for (const auto & it : trend_component_model_table.getEqNums())
    {
      set<expr_t> lhs;
      vector<int> orig_diff_var_vec;
      vector<bool> diff_vec;
      for (auto eqn : it.second)
        {
          // ensure no leads in equations
          if (equations[eqn]->get_arg2()->VarMinLag() <= 0)
            {
              cerr << "ERROR in trend component model Equation (#" << eqn << "). "
                   << "Leaded exogenous variables "
                   << "and leaded or contemporaneous endogenous variables not allowed in VAR"
                   << endl;
              exit(EXIT_FAILURE);
            }

          // save lhs variables
          equations[eqn]->get_arg1()->collectVARLHSVariable(lhs);

          equations[eqn]->get_arg1()->countDiffs() > 0 ?
            diff_vec.push_back(true) : diff_vec.push_back(false);
          if (diff_vec.back())
            {
              set<pair<int, int>> diff_set;
              equations[eqn]->get_arg1()->collectDynamicVariables(SymbolType::endogenous, diff_set);

              if (diff_set.size() != 1)
                {
                  cerr << "ERROR: problem getting variable for LHS diff operator in equation "
                       << eqn << endl;
                  exit(EXIT_FAILURE);
                }
              orig_diff_var_vec.push_back(diff_set.begin()->first);
            }
          else
            orig_diff_var_vec.push_back(-1);

        }

      if (it.second.size() != lhs.size())
        {
          cerr << "ERROR: The LHS variables of the trend component model are not unique" << endl;
          exit(EXIT_FAILURE);
        }

      set<expr_t> lhs_static;
      for(const auto & lh : lhs)
        lhs_static.insert(lh->toStatic(static_model));

      vector<int> max_lag;
      for (auto eqn : it.second)
        max_lag.push_back(equations[eqn]->get_arg2()->VarMaxLag(static_model, lhs_static));
      lags[it.first] = max_lag;
      diff[it.first] = diff_vec;
      orig_diff_var[it.first] = orig_diff_var_vec;
    }
  trend_component_model_table.setDiff(diff);
  trend_component_model_table.setMaxLags(lags);
  trend_component_model_table.setOrigDiffVar(orig_diff_var);
}

void
DynamicModel::fillTrendComponentmodelTableAREC(ExprNode::subst_table_t &diff_subst_table) const
{
  map<string, map<tuple<int, int, int>, expr_t>> ARr, ECr;
  fillAutoregressiveMatrix(ARr, true);
  trend_component_model_table.setAR(ARr);
  fillErrorComponentMatrix(ECr, diff_subst_table);
  trend_component_model_table.setEC(ECr);
}

void
DynamicModel::addEquationsForVar()
{
  if (var_model_table.empty())
    return;
  map<string, pair<SymbolList, int>> var_symbol_list_and_order =
    var_model_table.getSymbolListAndOrder();

  // List of endogenous variables and the minimum lag value that must exist in the model equations
  map<string, int> var_endos_and_lags, model_endos_and_lags;
  for (const auto & it : var_symbol_list_and_order)
    for (auto & equation : equations)
      if (equation->isVarModelReferenced(it.first))
        {
          vector<string> symbol_list = it.second.first.get_symbols();
          int order = it.second.second;
          for (vector<string>::const_iterator it1 = symbol_list.begin();
               it1 != symbol_list.end(); it1++)
            if (order > 2)
              if (var_endos_and_lags.find(*it1) != var_endos_and_lags.end())
                var_endos_and_lags[*it1] = min(var_endos_and_lags[*it1], -1*order);
              else
                var_endos_and_lags[*it1] = -1*order;
          break;
        }

  if (var_endos_and_lags.empty())
    return;

  // Ensure that the minimum lag value exists in the model equations.
  // If not, add an equation for it
  for (auto & equation : equations)
    equation->getEndosAndMaxLags(model_endos_and_lags);

  int count = 0;
  for (map<string, int>::const_iterator it = var_endos_and_lags.begin();
       it != var_endos_and_lags.end(); it++)
    {
      map<string, int>::const_iterator it1 = model_endos_and_lags.find(it->first);
      if (it1 == model_endos_and_lags.end())
        cerr << "WARNING: Variable used in VAR that is not used in the model: " << it->first << endl;
      else
        if (it->second < it1->second)
          {
            int symb_id = symbol_table.getID(it->first);
            expr_t newvar = AddVariable(symb_id, it->second);
            expr_t auxvar = AddVariable(symbol_table.addVarModelEndoLagAuxiliaryVar(symb_id, it->second, newvar), 0);
            addEquation(AddEqual(newvar, auxvar), -1);
            addAuxEquation(AddEqual(newvar, auxvar));
            count++;
          }
    }

  if (count > 0)
    cout << "Accounting for var_model lags not in model block: added "
         << count << " auxiliary variables and equations." << endl;
}

vector<int>
DynamicModel::getUndiffLHSForPac(const string &aux_model_name,
                                 ExprNode::subst_table_t &diff_subst_table) const
{
  vector<expr_t> lhs_expr_t = trend_component_model_table.getLhsExprT(aux_model_name);
  vector<int> lhs = trend_component_model_table.getLhs(aux_model_name);
  vector<bool> diff = trend_component_model_table.getDiff(aux_model_name);
  vector<int> orig_diff_var = trend_component_model_table.getOrigDiffVar(aux_model_name);
  vector<int> eqnumber = trend_component_model_table.getEqNums(aux_model_name);
  vector<int> nontrend_eqnums = trend_component_model_table.getNonTargetEqNums(aux_model_name);

  for (auto eqn : nontrend_eqnums)
    {
      int i = 0;
      for (vector<int>::const_iterator it1 = eqnumber.begin();
           it1 != eqnumber.end(); it1++, i++)
        if (*it1 == eqn)
          break;

      if (eqnumber[i] != eqn)
        {
          cerr << "ERROR: equation " << eqn << " not found in VAR" << endl;
          exit(EXIT_FAILURE);
        }

      if (diff.at(i) != true)
        {
          cerr << "ERROR: the variable on the LHS of equation #" << eqn
               << " does not have the diff operator applied to it yet you are trying to undiff it."
               << endl;
          exit(EXIT_FAILURE);
        }

      bool printerr = false;
      ExprNode::subst_table_t::const_iterator it1;
      expr_t node = nullptr;
      expr_t aux_var = lhs_expr_t.at(i);
      for (it1 = diff_subst_table.begin(); it1 != diff_subst_table.end(); it1++)
        if (it1->second == aux_var)
          {
            node = const_cast<expr_t>(it1->first);
            break;
          }

      if (node == nullptr)
        {
          cerr << "Unexpected error encountered." << endl;
          exit(EXIT_FAILURE);
        }

      node = node->undiff();
      it1 = diff_subst_table.find(node);
      if (it1 == diff_subst_table.end())
        printerr = true;

      if (printerr)
        { // we have undiffed something like diff(x), hence x is not in diff_subst_table
          lhs_expr_t.at(i) = node;
          lhs.at(i) = dynamic_cast<VariableNode *>(node)->get_symb_id();
        }
      else
        {
          lhs_expr_t.at(i) = const_cast<expr_t>(it1->first);
          lhs.at(i) = const_cast<VariableNode *>(it1->second)->get_symb_id();
        }
    }
  return lhs;
}

void
DynamicModel::walkPacParameters()
{
  for (auto & equation : equations)
    {
      pair<int, int> lhs (-1, -1);
      pair<int, pair<vector<int>, vector<bool>>> ec_params_and_vars;
      set<pair<int, pair<int, int>>> ar_params_and_vars;
      set<pair<int, pair<pair<int, int>, double>>> non_optim_params_vars_and_scaling_factor;

      if (equation->containsPacExpectation())
        {
          int optim_share_index = -1;
          set<int> optim_share;
          expr_t optim_part = nullptr;
          expr_t non_optim_part = nullptr;
          equation->getPacLHS(lhs);
          int lhs_orig_symb_id = lhs.first;
          if (symbol_table.isAuxiliaryVariable(lhs_orig_symb_id))
            try
              {
                lhs_orig_symb_id = symbol_table.getOrigSymbIdForAuxVar(lhs_orig_symb_id);
              }
            catch (...)
              {
              }

          equation->get_arg2()->getPacOptimizingShareAndExprNodes(optim_share,
                                                                  optim_part,
                                                                  non_optim_part);
          if (optim_part == nullptr)
            equation->get_arg2()->getPacOptimizingPart(lhs_orig_symb_id, ec_params_and_vars, ar_params_and_vars);
          else
            {
              optim_share_index = *(optim_share.begin());
              optim_part->getPacOptimizingPart(lhs_orig_symb_id, ec_params_and_vars, ar_params_and_vars);
              non_optim_part->getPacNonOptimizingPart(non_optim_params_vars_and_scaling_factor);
            }
          equation->addParamInfoToPac(lhs,
                                      optim_share_index,
                                      ec_params_and_vars, ar_params_and_vars,
                                      non_optim_params_vars_and_scaling_factor);
        }
    }
}

int
DynamicModel::getPacMaxLag(const string &pac_model_name) const
{
  for (auto & equation : equations)
    if (equation->containsPacExpectation(pac_model_name))
      {
        set<pair<int, int>> endogs;
        equation->get_arg1()->collectDynamicVariables(SymbolType::endogenous, endogs);
        if (endogs.size() != 1)
          {
            cerr << "The LHS of the PAC equation may only be comprised of one endogenous variable"
                 << endl;
            exit(EXIT_FAILURE);
          }
        return equation->PacMaxLag((*(endogs.begin())).first);
      }
  return 0;
}

void
DynamicModel::fillPacExpectationVarInfo(string &pac_model_name,
                                        vector<int> &lhs,
                                        int max_lag,
                                        int pac_max_lag,
                                        vector<bool> &nonstationary,
                                        int growth_symb_id)
{
  for (size_t i = 0; i < equations.size(); i++)
    equations[i]->fillPacExpectationVarInfo(pac_model_name, lhs, max_lag,
                                            pac_max_lag, nonstationary, growth_symb_id, i);
}

void
DynamicModel::substitutePacExpectation()
{
  map<const PacExpectationNode *, const BinaryOpNode *> subst_table;
  for (auto & it : local_variables_table)
    it.second = it.second->substitutePacExpectation(subst_table);

  for (auto & equation : equations)
    {
      auto *substeq = dynamic_cast<BinaryOpNode *>(equation->substitutePacExpectation(subst_table));
      assert(substeq != nullptr);
      equation = substeq;
    }

  for (map<const PacExpectationNode *, const BinaryOpNode *>::const_iterator it = subst_table.begin();
       it != subst_table.end(); it++)
    pac_expectation_info.insert(const_cast<PacExpectationNode *>(it->first));
}

void
DynamicModel::computingPass(bool jacobianExo, int derivsOrder, int paramsDerivsOrder,
                            const eval_context_t &eval_context, bool no_tmp_terms, bool block, bool use_dll,
                            bool bytecode, bool nopreprocessoroutput, bool linear_decomposition)
{
  assert(jacobianExo || (derivsOrder < 2 && paramsDerivsOrder == 0));

  initializeVariablesAndEquations();

  // Prepare for derivation
  computeDerivIDs();

  // Computes dynamic jacobian columns, must be done after computeDerivIDs()
  computeDynJacobianCols(jacobianExo);

  // Compute derivatives w.r. to all endogenous, and possibly exogenous and exogenous deterministic
  set<int> vars;
  for (deriv_id_table_t::const_iterator it = deriv_id_table.begin();
       it != deriv_id_table.end(); it++)
    {
      SymbolType type = symbol_table.getType(it->first.first);
      if (type == SymbolType::endogenous || (jacobianExo && (type == SymbolType::exogenous || type == SymbolType::exogenousDet)))
        vars.insert(it->second);
    }

  // Launch computations
  if (!nopreprocessoroutput)
    cout << "Computing " << (linear_decomposition ? "nonlinear " : "")
         << "dynamic model derivatives (order " << derivsOrder << ")." << endl;

  computeDerivatives(derivsOrder, vars);

  if (paramsDerivsOrder > 0)
    {
      if (!nopreprocessoroutput)
        cout << "Computing dynamic model derivatives w.r.t. parameters (order " << paramsDerivsOrder << ")." << endl;
      computeParamsDerivatives(paramsDerivsOrder);
    }

  jacob_map_t contemporaneous_jacobian, static_jacobian;
  map<tuple<int, int, int>, expr_t> first_order_endo_derivatives;
  // for each block contains pair<Size, Feddback_variable>
  vector<pair<int, int> > blocks;
  vector<unsigned int> n_static, n_forward, n_backward, n_mixed;

  if (linear_decomposition)
    {
      first_order_endo_derivatives = collect_first_order_derivatives_endogenous();
      is_equation_linear = equationLinear(first_order_endo_derivatives);

      evaluateAndReduceJacobian(eval_context, contemporaneous_jacobian, static_jacobian, dynamic_jacobian, cutoff, false);

      if (!computeNaturalNormalization())
        computeNonSingularNormalization(contemporaneous_jacobian, cutoff, static_jacobian, dynamic_jacobian);

      lag_lead_vector_t equation_lag_lead, variable_lag_lead;

      blocks = select_non_linear_equations_and_variables(is_equation_linear, dynamic_jacobian, equation_reordered, variable_reordered,
                                                         inv_equation_reordered, inv_variable_reordered,
                                                         equation_lag_lead, variable_lag_lead,
                                                         n_static, n_forward, n_backward, n_mixed);

      equation_type_and_normalized_equation = equationTypeDetermination(first_order_endo_derivatives, variable_reordered, equation_reordered, 0);
      prologue = 0;
      epilogue = 0;

      block_type_firstequation_size_mfs = reduceBlocksAndTypeDetermination(dynamic_jacobian, blocks, equation_type_and_normalized_equation, variable_reordered, equation_reordered, n_static, n_forward, n_backward, n_mixed, block_col_type, linear_decomposition);

      computeChainRuleJacobian(blocks_derivatives);

      blocks_linear = BlockLinear(blocks_derivatives, variable_reordered);

      collect_block_first_order_derivatives();

      collectBlockVariables();

      global_temporary_terms = true;
      if (!no_tmp_terms)
        computeTemporaryTermsOrdered();

    }

  if (block)
    {

      evaluateAndReduceJacobian(eval_context, contemporaneous_jacobian, static_jacobian, dynamic_jacobian, cutoff, false);

      computeNonSingularNormalization(contemporaneous_jacobian, cutoff, static_jacobian, dynamic_jacobian);

      computePrologueAndEpilogue(static_jacobian, equation_reordered, variable_reordered);

      first_order_endo_derivatives = collect_first_order_derivatives_endogenous();

      equation_type_and_normalized_equation = equationTypeDetermination(first_order_endo_derivatives, variable_reordered, equation_reordered, mfs);

      if (!nopreprocessoroutput)
        cout << "Finding the optimal block decomposition of the model ..." << endl;

      lag_lead_vector_t equation_lag_lead, variable_lag_lead;

      computeBlockDecompositionAndFeedbackVariablesForEachBlock(static_jacobian, dynamic_jacobian, equation_reordered, variable_reordered, blocks, equation_type_and_normalized_equation, false, true, mfs, inv_equation_reordered, inv_variable_reordered, equation_lag_lead, variable_lag_lead, n_static, n_forward, n_backward, n_mixed);

      block_type_firstequation_size_mfs = reduceBlocksAndTypeDetermination(dynamic_jacobian, blocks, equation_type_and_normalized_equation, variable_reordered, equation_reordered, n_static, n_forward, n_backward, n_mixed, block_col_type, linear_decomposition);

      printBlockDecomposition(blocks);

      computeChainRuleJacobian(blocks_derivatives);

      blocks_linear = BlockLinear(blocks_derivatives, variable_reordered);

      collect_block_first_order_derivatives();

      collectBlockVariables();

      global_temporary_terms = true;
      if (!no_tmp_terms)
        computeTemporaryTermsOrdered();
      int k = 0;
      equation_block.resize(equations.size());
      variable_block_lead_lag = vector<tuple<int, int, int>>(equations.size());
      for (unsigned int i = 0; i < getNbBlocks(); i++)
        {
          for (unsigned int j = 0; j < getBlockSize(i); j++)
            {
              equation_block[equation_reordered[k]] = i;
              int l = variable_reordered[k];
              variable_block_lead_lag[l] = { i, variable_lag_lead[l].first, variable_lag_lead[l].second };
              k++;
            }
        }
    }
  else
    {
      computeTemporaryTerms(!use_dll, no_tmp_terms);
      if (bytecode && !no_tmp_terms)
        computeTemporaryTermsMapping();

      /* Must be called after computeTemporaryTerms(), because it depends on
         temporary_terms_mlv to be filled */
      if (paramsDerivsOrder > 0 && !no_tmp_terms)
        computeParamsDerivativesTemporaryTerms();
    }
}

void
DynamicModel::computeXrefs()
{
  int i = 0;
  for (auto & equation : equations)
    {
      ExprNode::EquationInfo ei;
      equation->computeXrefs(ei);
      xrefs[i++] = ei;
    }

  i = 0;
  for (map<int, ExprNode::EquationInfo>::const_iterator it = xrefs.begin();
       it != xrefs.end(); it++, i++)
    {
      computeRevXref(xref_param, it->second.param, i);
      computeRevXref(xref_endo, it->second.endo, i);
      computeRevXref(xref_exo, it->second.exo, i);
      computeRevXref(xref_exo_det, it->second.exo_det, i);
    }
}

void
DynamicModel::computeRevXref(map<pair<int, int>, set<int>> &xrefset, const set<pair<int, int>> &eiref, int eqn)
{
  for (const auto & it : eiref)
    {
      set<int> eq;
      if (xrefset.find(it) != xrefset.end())
        eq = xrefset[it];
      eq.insert(eqn);
      xrefset[it] = eq;
    }
}

void
DynamicModel::writeXrefs(ostream &output) const
{
  output << "M_.xref1.param = cell(1, M_.eq_nbr);" << endl
         << "M_.xref1.endo = cell(1, M_.eq_nbr);" << endl
         << "M_.xref1.exo = cell(1, M_.eq_nbr);" << endl
         << "M_.xref1.exo_det = cell(1, M_.eq_nbr);" << endl;
  int i = 1;
  for (auto it = xrefs.begin();
       it != xrefs.end(); it++, i++)
    {
      output << "M_.xref1.param{" << i << "} = [ ";
      for (const auto & it1 : it->second.param)
        output << symbol_table.getTypeSpecificID(it1.first) + 1 << " ";
      output << "];" << endl;

      output << "M_.xref1.endo{" << i << "} = [ ";
      for (const auto & it1 : it->second.endo)
        output << "struct('id', " << symbol_table.getTypeSpecificID(it1.first) + 1 << ", 'shift', " << it1.second << ");";
      output << "];" << endl;

      output << "M_.xref1.exo{" << i << "} = [ ";
      for (const auto & it1 : it->second.exo)
        output << "struct('id', " << symbol_table.getTypeSpecificID(it1.first) + 1 << ", 'shift', " << it1.second << ");";
      output << "];" << endl;

      output << "M_.xref1.exo_det{" << i << "} = [ ";
      for (const auto & it1 : it->second.exo_det)
        output << "struct('id', " << symbol_table.getTypeSpecificID(it1.first) + 1 << ", 'shift', " << it1.second << ");";
      output << "];" << endl;
    }

  output << "M_.xref2.param = cell(1, M_.param_nbr);" << endl
         << "M_.xref2.endo = cell(1, M_.endo_nbr);" << endl
         << "M_.xref2.exo = cell(1, M_.exo_nbr);" << endl
         << "M_.xref2.exo_det = cell(1, M_.exo_det_nbr);" << endl;
  writeRevXrefs(output, xref_param, "param");
  writeRevXrefs(output, xref_endo, "endo");
  writeRevXrefs(output, xref_exo, "exo");
  writeRevXrefs(output, xref_exo_det, "exo_det");
}

void
DynamicModel::writeRevXrefs(ostream &output, const map<pair<int, int>, set<int>> &xrefmap, const string &type) const
{
  int last_tsid = -1;
  for (const auto & it : xrefmap)
    {
      int tsid = symbol_table.getTypeSpecificID(it.first.first) + 1;
      output << "M_.xref2." << type << "{" << tsid << "} = [ ";
      if (last_tsid == tsid)
        output << "M_.xref2." << type << "{" << tsid << "}; ";
      else
        last_tsid = tsid;

      for (auto it1 = it.second.begin();
           it1 != it.second.end(); it1++)
        if (type == "param")
          output << *it1 + 1 << " ";
        else
          output << "struct('shift', " << it.first.second << ", 'eq', " << *it1+1 << ");";
      output << "];" << endl;
    }
}

map<tuple<int, int, int, int, int>, int>
DynamicModel::get_Derivatives(int block)
{
  int max_lag, max_lead;
  map<tuple<int, int, int, int, int>, int> Derivatives;
  BlockSimulationType simulation_type = getBlockSimulationType(block);
  if (simulation_type == EVALUATE_BACKWARD || simulation_type == EVALUATE_FORWARD)
    {
      max_lag  = 1;
      max_lead = 1;
      setBlockLeadLag(block, max_lag, max_lead);
    }
  else
    {
      max_lag  = getBlockMaxLag(block);
      max_lead = getBlockMaxLead(block);
    }
  int block_size = getBlockSize(block);
  int block_nb_recursive = block_size - getBlockMfs(block);
  for (int lag = -max_lag; lag <= max_lead; lag++)
    {
      for (int eq = 0; eq < block_size; eq++)
        {
          int eqr = getBlockEquationID(block, eq);
          for (int var = 0; var < block_size; var++)
            {
              int varr = getBlockVariableID(block, var);
              if (dynamic_jacobian.find({ lag, eqr, varr }) != dynamic_jacobian.end())
                {
                  bool OK = true;
                  auto its = Derivatives.find({ lag, eq, var, eqr, varr });
                  if (its != Derivatives.end())
                    {
                      if (its->second == 2)
                        OK = false;
                    }

                  if (OK)
                    {
                      if (getBlockEquationType(block, eq) == E_EVALUATE_S && eq < block_nb_recursive)
                        //It's a normalized equation, we have to recompute the derivative using chain rule derivative function
                        Derivatives[{ lag, eq, var, eqr, varr }] = 1;
                      else
                        //It's a feedback equation we can use the derivatives
                        Derivatives[{ lag, eq, var, eqr, varr }] = 0;
                    }
                  if (var < block_nb_recursive)
                    {
                      int eqs = getBlockEquationID(block, var);
                      for (int vars = block_nb_recursive; vars < block_size; vars++)
                        {
                          int varrs = getBlockVariableID(block, vars);
                          //A new derivative needs to be computed using the chain rule derivative function (a feedback variable appears in a recursive equation)
                          if (Derivatives.find({ lag, var, vars, eqs, varrs }) != Derivatives.end())
                            Derivatives[{ lag, eq, vars, eqr, varrs }] = 2;
                        }
                    }
                }
            }
        }
    }
  return Derivatives;
}

void
DynamicModel::computeChainRuleJacobian(blocks_derivatives_t &blocks_endo_derivatives)
{
  map<int, expr_t> recursive_variables;
  unsigned int nb_blocks = getNbBlocks();
  blocks_endo_derivatives = blocks_derivatives_t(nb_blocks);
  for (unsigned int block = 0; block < nb_blocks; block++)
    {
      block_derivatives_equation_variable_laglead_nodeid_t tmp_derivatives;
      recursive_variables.clear();
      int block_size = getBlockSize(block);
      int block_nb_mfs = getBlockMfs(block);
      int block_nb_recursives = block_size - block_nb_mfs;
      blocks_endo_derivatives.push_back(block_derivatives_equation_variable_laglead_nodeid_t(0));
      for (int i = 0; i < block_nb_recursives; i++)
        {
          if (getBlockEquationType(block, i) == E_EVALUATE_S)
            recursive_variables[getDerivID(symbol_table.getID(SymbolType::endogenous, getBlockVariableID(block, i)), 0)] = getBlockEquationRenormalizedExpr(block, i);
          else
            recursive_variables[getDerivID(symbol_table.getID(SymbolType::endogenous, getBlockVariableID(block, i)), 0)] = getBlockEquationExpr(block, i);
        }
      auto Derivatives = get_Derivatives(block);
      for (const auto &it : Derivatives)
        {
          int Deriv_type = it.second;
          int lag, eq, var, eqr, varr;
          tie(lag, eq, var, eqr, varr) = it.first;
          if (Deriv_type == 0)
            first_chain_rule_derivatives[{ eqr, varr, lag }] = derivatives[1][{ eqr, getDerivID(symbol_table.getID(SymbolType::endogenous, varr), lag) }];
          else if (Deriv_type == 1)
            first_chain_rule_derivatives[{ eqr, varr, lag }] = (equation_type_and_normalized_equation[eqr].second)->getChainRuleDerivative(getDerivID(symbol_table.getID(SymbolType::endogenous, varr), lag), recursive_variables);
          else if (Deriv_type == 2)
            {
              if (getBlockEquationType(block, eq) == E_EVALUATE_S && eq < block_nb_recursives)
                first_chain_rule_derivatives[{ eqr, varr, lag }] = (equation_type_and_normalized_equation[eqr].second)->getChainRuleDerivative(getDerivID(symbol_table.getID(SymbolType::endogenous, varr), lag), recursive_variables);
              else
                first_chain_rule_derivatives[{ eqr, varr, lag }] = equations[eqr]->getChainRuleDerivative(getDerivID(symbol_table.getID(SymbolType::endogenous, varr), lag), recursive_variables);
            }
          tmp_derivatives.emplace_back(eq, var, lag, first_chain_rule_derivatives[{ eqr, varr, lag }]);
        }
      blocks_endo_derivatives[block] = tmp_derivatives;
    }
}

void
DynamicModel::collect_block_first_order_derivatives()
{
  //! vector for an equation or a variable indicates the block number
  vector<int> equation_2_block(equation_reordered.size()), variable_2_block(variable_reordered.size());
  unsigned int nb_blocks = getNbBlocks();
  for (unsigned int block = 0; block < nb_blocks; block++)
    {
      unsigned int block_size = getBlockSize(block);
      for (unsigned int i = 0; i < block_size; i++)
        {
          equation_2_block[getBlockEquationID(block, i)] = block;
          variable_2_block[getBlockVariableID(block, i)] = block;
        }
    }
  other_endo_block = vector<lag_var_t>(nb_blocks);
  exo_block = vector<lag_var_t>(nb_blocks);
  exo_det_block = vector<lag_var_t>(nb_blocks);
  derivative_endo = vector<derivative_t>(nb_blocks);
  derivative_other_endo = vector<derivative_t>(nb_blocks);
  derivative_exo = vector<derivative_t>(nb_blocks);
  derivative_exo_det = vector<derivative_t>(nb_blocks);
  endo_max_leadlag_block = vector<pair<int, int>>(nb_blocks, { 0, 0 });
  other_endo_max_leadlag_block = vector<pair<int, int>>(nb_blocks, { 0, 0 });
  exo_max_leadlag_block = vector<pair<int, int>>(nb_blocks, { 0, 0 });
  exo_det_max_leadlag_block = vector<pair<int, int>>(nb_blocks, { 0, 0 });
  max_leadlag_block = vector<pair<int, int>>(nb_blocks, { 0, 0 });
  for (auto & first_derivative : derivatives[1])
    {
      int eq = first_derivative.first[0];
      int var = symbol_table.getTypeSpecificID(getSymbIDByDerivID(first_derivative.first[1]));
      int lag = getLagByDerivID(first_derivative.first[1]);
      int block_eq = equation_2_block[eq];
      int block_var = 0;
      derivative_t tmp_derivative;
      lag_var_t lag_var;
      switch (getTypeByDerivID(first_derivative.first[1]))
        {
        case SymbolType::endogenous:
          block_var = variable_2_block[var];
          if (block_eq == block_var)
            {
              if (lag < 0 && lag < -endo_max_leadlag_block[block_eq].first)
                endo_max_leadlag_block[block_eq] = { -lag, endo_max_leadlag_block[block_eq].second };
              if (lag > 0 && lag > endo_max_leadlag_block[block_eq].second)
                endo_max_leadlag_block[block_eq] = { endo_max_leadlag_block[block_eq].first, lag };
              tmp_derivative = derivative_endo[block_eq];
              tmp_derivative[{ lag, eq, var }] = derivatives[1][{ eq, getDerivID(symbol_table.getID(SymbolType::endogenous, var), lag) }];
              derivative_endo[block_eq] = tmp_derivative;
            }
          else
            {
              if (lag < 0 && lag < -other_endo_max_leadlag_block[block_eq].first)
                other_endo_max_leadlag_block[block_eq] = { -lag, other_endo_max_leadlag_block[block_eq].second };
              if (lag > 0 && lag > other_endo_max_leadlag_block[block_eq].second)
                other_endo_max_leadlag_block[block_eq] = { other_endo_max_leadlag_block[block_eq].first, lag };
              tmp_derivative = derivative_other_endo[block_eq];
              {
                auto it = block_other_endo_index.find(block_eq);
                if (it == block_other_endo_index.end())
                  block_other_endo_index[block_eq][var] = 0;
                else
                  {
                    auto it1 = it->second.find(var);
                    if (it1 == it->second.end())
                      {
                        int size = block_other_endo_index[block_eq].size();
                        block_other_endo_index[block_eq][var] = size;
                      }
                  }
              }
              tmp_derivative[{ lag, eq, var }] = derivatives[1][{ eq, getDerivID(symbol_table.getID(SymbolType::endogenous, var), lag) }];
              derivative_other_endo[block_eq] = tmp_derivative;
              lag_var = other_endo_block[block_eq];
              if (lag_var.find(lag) == lag_var.end())
                lag_var[lag].clear();
              lag_var[lag].insert(var);
              other_endo_block[block_eq] = lag_var;
            }
          break;
        case SymbolType::exogenous:
          if (lag < 0 && lag < -exo_max_leadlag_block[block_eq].first)
            exo_max_leadlag_block[block_eq] = { -lag, exo_max_leadlag_block[block_eq].second };
          if (lag > 0 && lag > exo_max_leadlag_block[block_eq].second)
            exo_max_leadlag_block[block_eq] = { exo_max_leadlag_block[block_eq].first, lag };
          tmp_derivative = derivative_exo[block_eq];
          {
            auto it = block_exo_index.find(block_eq);
            if (it == block_exo_index.end())
              block_exo_index[block_eq][var] = 0;
            else
              {
                auto it1 = it->second.find(var);
                if (it1 == it->second.end())
                  {
                    int size = block_exo_index[block_eq].size();
                    block_exo_index[block_eq][var] = size;
                  }
              }
          }
          tmp_derivative[{ lag, eq, var }] = derivatives[1][{ eq, getDerivID(symbol_table.getID(SymbolType::exogenous, var), lag) }];
          derivative_exo[block_eq] = tmp_derivative;
          lag_var = exo_block[block_eq];
          if (lag_var.find(lag) == lag_var.end())
            lag_var[lag].clear();
          lag_var[lag].insert(var);
          exo_block[block_eq] = lag_var;
          break;
        case SymbolType::exogenousDet:
          if (lag < 0 && lag < -exo_det_max_leadlag_block[block_eq].first)
            exo_det_max_leadlag_block[block_eq] = { -lag, exo_det_max_leadlag_block[block_eq].second };
          if (lag > 0 && lag > exo_det_max_leadlag_block[block_eq].second)
            exo_det_max_leadlag_block[block_eq] = { exo_det_max_leadlag_block[block_eq].first, lag };
          tmp_derivative = derivative_exo_det[block_eq];
          {
            auto it = block_det_exo_index.find(block_eq);
            if (it == block_det_exo_index.end())
              block_det_exo_index[block_eq][var] = 0;
            else
              {
                auto it1 = it->second.find(var);
                if (it1 == it->second.end())
                  {
                    int size = block_det_exo_index[block_eq].size();
                    block_det_exo_index[block_eq][var] = size;
                  }
              }
          }
          tmp_derivative[{ lag, eq, var }] = derivatives[1][{ eq, getDerivID(symbol_table.getID(SymbolType::exogenous, var), lag) }];
          derivative_exo_det[block_eq] = tmp_derivative;
          lag_var = exo_det_block[block_eq];
          if (lag_var.find(lag) == lag_var.end())
            lag_var[lag].clear();
          lag_var[lag].insert(var);
          exo_det_block[block_eq] = lag_var;
          break;
        default:
          break;
        }
      if (lag < 0 && lag < -max_leadlag_block[block_eq].first)
        max_leadlag_block[block_eq] = { -lag, max_leadlag_block[block_eq].second };
      if (lag > 0 && lag > max_leadlag_block[block_eq].second)
        max_leadlag_block[block_eq] = { max_leadlag_block[block_eq].first, lag };
    }

}

void
DynamicModel::collectBlockVariables()
{
  for (unsigned int block = 0; block < getNbBlocks(); block++)
    {
      int prev_var = -1;
      int prev_lag = -999999999;
      int count_col_exo = 0;
      var_t tmp_var_exo;
      for (const auto &it : exo_block[block])
        {
          int lag = it.first;
          for (int var : it.second)
            {
              tmp_var_exo.insert(var);
              if (prev_var != var || prev_lag != lag)
                {
                  prev_var = var;
                  prev_lag = lag;
                  count_col_exo++;
                }
            }
        }
      block_var_exo.emplace_back(tmp_var_exo, count_col_exo);
    }
}

void
DynamicModel::writeDynamicFile(const string &basename, bool block, bool linear_decomposition, bool bytecode, bool use_dll, const string &mexext, const boost::filesystem::path &matlabroot, const boost::filesystem::path &dynareroot, int order, bool julia) const
{
  if (block && bytecode)
    writeModelEquationsCode_Block(basename, map_idx, linear_decomposition);
  else if (!block && bytecode)
    {
      if (linear_decomposition)
        writeModelEquationsCode_Block(basename, map_idx, linear_decomposition);
      writeModelEquationsCode(basename, map_idx);
    }
  else if (block && !bytecode)
    writeSparseDynamicMFile(basename);
  else if (use_dll)
    {
      writeDynamicCFile(basename, order);
      compileDll(basename, "dynamic", mexext, matlabroot, dynareroot);
    }
  else if (julia)
    writeDynamicJuliaFile(basename);
  else
    {
      writeDynamicMFile(basename);
      writeSetAuxiliaryVariables(basename, julia);
    }
}

void
DynamicModel::writeSetAuxiliaryVariables(const string &basename, const bool julia) const
{
  ostringstream output_func_body;
  writeAuxVarRecursiveDefinitions(output_func_body, ExprNodeOutputType::matlabDseries);

  if (output_func_body.str().empty())
    return;

  string func_name = julia ? basename + "_dynamic_set_auxiliary_series" : "dynamic_set_auxiliary_series";
  string filename = julia ? func_name + ".jl" : packageDir(basename) + "/" + func_name + ".m";
  string comment = julia ? "#" : "%";

  ofstream output;
  output.open(filename, ios::out | ios::binary);
  if (!output.is_open())
    {
      cerr << "ERROR: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }

  output << "function ds = " << func_name + "(ds, params)" << endl
         << comment << endl
         << comment << " Status : Computes Auxiliary variables of the dynamic model and returns a dseries" << endl
         << comment << endl
         << comment << " Warning : this file is generated automatically by Dynare" << endl
         << comment << "           from model file (.mod)" << endl << endl
         << output_func_body.str();

  output.close();
}

void
DynamicModel::writeAuxVarRecursiveDefinitions(ostream &output, ExprNodeOutputType output_type) const
{
  deriv_node_temp_terms_t tef_terms;
  temporary_terms_t temporary_terms;
  temporary_terms_idxs_t temporary_terms_idxs;
  for (auto aux_equation : aux_equations)
    if (dynamic_cast<ExprNode *>(aux_equation)->containsExternalFunction())
      dynamic_cast<ExprNode *>(aux_equation)->writeExternalFunctionOutput(output, output_type,
                                                                              temporary_terms,
                                                                              temporary_terms_idxs,
                                                                              tef_terms);
  for (auto aux_equation : aux_equations)
    {
      dynamic_cast<ExprNode *>(aux_equation)->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
      output << ";" << endl;
    }
}

void
DynamicModel::replaceMyEquations(DynamicModel &dynamic_model) const
{
  dynamic_model.equations.clear();
  for (size_t i = 0; i < equations.size(); i++)
    dynamic_model.addEquation(equations[i]->clone(dynamic_model),
                              equations_lineno[i]);
}

void
DynamicModel::computeRamseyPolicyFOCs(const StaticModel &static_model, const bool nopreprocessoroutput)
{
  // Add aux LM to constraints in equations
  // equation[i]->lhs = rhs becomes equation[i]->MULT_(i+1)*(lhs-rhs) = 0
  int i;
  for (i = 0; i < (int) equations.size(); i++)
    {
      auto *substeq = dynamic_cast<BinaryOpNode *>(equations[i]->addMultipliersToConstraints(i));
      assert(substeq != nullptr);
      equations[i] = substeq;
    }
  if (!nopreprocessoroutput)
    cout << "Ramsey Problem: added " << i << " Multipliers." << endl;

  // Add Planner Objective to equations to include in computeDerivIDs
  assert(static_model.equations.size() == 1);
  addEquation(static_model.equations[0]->clone(*this), static_model.equations_lineno[0]);

  // Get max endo lead and max endo lag
  set<pair<int, int>> dynvars;
  int max_eq_lead = 0;
  int max_eq_lag = 0;
  for (auto & equation : equations)
    equation->collectDynamicVariables(SymbolType::endogenous, dynvars);

  for (const auto & dynvar : dynvars)
    {
      int lag = dynvar.second;
      if (max_eq_lead < lag)
        max_eq_lead = lag;
      else if (-max_eq_lag > lag)
        max_eq_lag = -lag;
    }

  // Get Discount Factor
  assert(symbol_table.exists("optimal_policy_discount_factor"));
  int symb_id = symbol_table.getID("optimal_policy_discount_factor");
  assert(symbol_table.getType(symb_id) == SymbolType::parameter);
  expr_t discount_factor_node = AddVariable(symb_id, 0);

  // Create (modified) Lagrangian (so that we can take the derivative once at time t)
  expr_t lagrangian = Zero;
  for (i = 0; i < (int) equations.size(); i++)
    for (int lag = -max_eq_lag; lag <= max_eq_lead; lag++)
      {
        expr_t dfpower = nullptr;
        std::stringstream lagstream;
        lagstream << abs(lag);
        if (lag < 0)
          dfpower = AddNonNegativeConstant(lagstream.str());
        else if (lag == 0)
          dfpower = Zero;
        else
          dfpower = AddMinus(Zero, AddNonNegativeConstant(lagstream.str()));

        lagrangian = AddPlus(AddTimes(AddPower(discount_factor_node, dfpower),
                                      equations[i]->getNonZeroPartofEquation()->decreaseLeadsLags(lag)), lagrangian);
      }

  equations.clear();
  addEquation(AddEqual(lagrangian, Zero), -1);
  computeDerivIDs();

  //Compute derivatives and overwrite equations
  vector<expr_t> neweqs;
  for (deriv_id_table_t::const_iterator it = deriv_id_table.begin();
       it != deriv_id_table.end(); it++)
    // For all endogenous variables with zero lag
    if (symbol_table.getType(it->first.first)  == SymbolType::endogenous && it->first.second == 0)
      neweqs.push_back(AddEqual(equations[0]->getNonZeroPartofEquation()->getDerivative(it->second), Zero));

  // Add new equations
  equations.clear();
  for (auto & neweq : neweqs)
    addEquation(neweq, -1);
}

void
DynamicModel::toNonlinearPart(DynamicModel &non_linear_equations_dynamic_model) const
{
  // Convert model local variables (need to be done first)
  for (const auto & it : local_variables_table)
    non_linear_equations_dynamic_model.AddLocalVariable(it.first, it.second);
}

bool
DynamicModel::ParamUsedWithLeadLag() const
{
  return ParamUsedWithLeadLagInternal();
}

set<int>
DynamicModel::findUnusedEndogenous()
{
  set<int> usedEndo, unusedEndo;
  for (auto & equation : equations)
    equation->collectVariables(SymbolType::endogenous, usedEndo);
  set<int> allEndo = symbol_table.getEndogenous();
  set_difference(allEndo.begin(), allEndo.end(),
                 usedEndo.begin(), usedEndo.end(),
                 inserter(unusedEndo, unusedEndo.begin()));
  return unusedEndo;
}

set<int>
DynamicModel::findUnusedExogenous()
{
  set<int> usedExo, unusedExo, unobservedExo;
  for (auto & equation : equations)
    equation->collectVariables(SymbolType::exogenous, usedExo);
  set<int> observedExo = symbol_table.getObservedExogenous();
  set<int> allExo = symbol_table.getExogenous();
  set_difference(allExo.begin(), allExo.end(),
                 observedExo.begin(), observedExo.end(),
                 inserter(unobservedExo, unobservedExo.begin()));
  set_difference(unobservedExo.begin(), unobservedExo.end(),
                 usedExo.begin(), usedExo.end(),
                 inserter(unusedExo, unusedExo.begin()));
  return unusedExo;
}

void
DynamicModel::setLeadsLagsOrig()
{
  set<pair<int, int>> dynvars;

  for (auto & equation : equations)
    {
      equation->collectDynamicVariables(SymbolType::endogenous, dynvars);
      equation->collectDynamicVariables(SymbolType::exogenous, dynvars);
      equation->collectDynamicVariables(SymbolType::exogenousDet, dynvars);
    }

    for (const auto & dynvar : dynvars)
    {
      int lag = dynvar.second;
      SymbolType type = symbol_table.getType(dynvar.first);

      if (max_lead_orig < lag)
        max_lead_orig= lag;
      else if (-max_lag_orig > lag)
        max_lag_orig = -lag;

      switch (type)
        {
        case SymbolType::endogenous:
          if (max_endo_lead_orig < lag)
            max_endo_lead_orig = lag;
          else if (-max_endo_lag_orig > lag)
            max_endo_lag_orig = -lag;
          break;
        case SymbolType::exogenous:
          if (max_exo_lead_orig < lag)
            max_exo_lead_orig = lag;
          else if (-max_exo_lag_orig > lag)
            max_exo_lag_orig = -lag;
          break;
        case SymbolType::exogenousDet:
          if (max_exo_det_lead_orig < lag)
            max_exo_det_lead_orig = lag;
          else if (-max_exo_det_lag_orig > lag)
            max_exo_det_lag_orig = -lag;
          break;
        default:
          break;
        }
    }
}

void
DynamicModel::computeDerivIDs()
{
  set<pair<int, int>> dynvars;

  for (auto & equation : equations)
    equation->collectDynamicVariables(SymbolType::endogenous, dynvars);

  dynJacobianColsNbr = dynvars.size();

  for (auto & equation : equations)
    {
      equation->collectDynamicVariables(SymbolType::exogenous, dynvars);
      equation->collectDynamicVariables(SymbolType::exogenousDet, dynvars);
      equation->collectDynamicVariables(SymbolType::parameter, dynvars);
      equation->collectDynamicVariables(SymbolType::trend, dynvars);
      equation->collectDynamicVariables(SymbolType::logTrend, dynvars);
    }

  for (const auto & dynvar : dynvars)
    {
      int lag = dynvar.second;
      SymbolType type = symbol_table.getType(dynvar.first);

      /* Setting maximum and minimum lags.

         We don't want these to be affected by lead/lags on parameters: they
         are accepted for facilitating variable flipping, but are simply
         ignored. */
      if (max_lead < lag && type != SymbolType::parameter)
        max_lead = lag;
      else if (-max_lag > lag && type != SymbolType::parameter)
        max_lag = -lag;

      switch (type)
        {
        case SymbolType::endogenous:
          if (max_endo_lead < lag)
            max_endo_lead = lag;
          else if (-max_endo_lag > lag)
            max_endo_lag = -lag;
          break;
        case SymbolType::exogenous:
          if (max_exo_lead < lag)
            max_exo_lead = lag;
          else if (-max_exo_lag > lag)
            max_exo_lag = -lag;
          break;
        case SymbolType::exogenousDet:
          if (max_exo_det_lead < lag)
            max_exo_det_lead = lag;
          else if (-max_exo_det_lag > lag)
            max_exo_det_lag = -lag;
          break;
        default:
          break;
        }

      // Create a new deriv_id
      int deriv_id = deriv_id_table.size();

      deriv_id_table[dynvar] = deriv_id;
      inv_deriv_id_table.push_back(dynvar);
    }
}

SymbolType
DynamicModel::getTypeByDerivID(int deriv_id) const noexcept(false)
{
  return symbol_table.getType(getSymbIDByDerivID(deriv_id));
}

int
DynamicModel::getLagByDerivID(int deriv_id) const noexcept(false)
{
  if (deriv_id < 0 || deriv_id >= (int) inv_deriv_id_table.size())
    throw UnknownDerivIDException();

  return inv_deriv_id_table[deriv_id].second;
}

int
DynamicModel::getSymbIDByDerivID(int deriv_id) const noexcept(false)
{
  if (deriv_id < 0 || deriv_id >= (int) inv_deriv_id_table.size())
    throw UnknownDerivIDException();

  return inv_deriv_id_table[deriv_id].first;
}

int
DynamicModel::getDerivID(int symb_id, int lag) const noexcept(false)
{
  auto it = deriv_id_table.find({ symb_id, lag });
  if (it == deriv_id_table.end())
    throw UnknownDerivIDException();
  else
    return it->second;
}

void
DynamicModel::addAllParamDerivId(set<int> &deriv_id_set)
{
  for (size_t i = 0; i < inv_deriv_id_table.size(); i++)
    if (symbol_table.getType(inv_deriv_id_table[i].first) == SymbolType::parameter)
      deriv_id_set.insert(i);
}

void
DynamicModel::computeDynJacobianCols(bool jacobianExo)
{
  /* Sort the dynamic endogenous variables by lexicographic order over (lag, type_specific_symbol_id)
     and fill the dynamic columns for exogenous and exogenous deterministic */
  map<pair<int, int>, int> ordered_dyn_endo;

  for (deriv_id_table_t::const_iterator it = deriv_id_table.begin();
       it != deriv_id_table.end(); it++)
    {
      const int &symb_id = it->first.first;
      const int &lag = it->first.second;
      const int &deriv_id = it->second;
      SymbolType type = symbol_table.getType(symb_id);
      int tsid = symbol_table.getTypeSpecificID(symb_id);

      switch (type)
        {
        case SymbolType::endogenous:
          ordered_dyn_endo[{ lag, tsid }] = deriv_id;
          break;
        case SymbolType::exogenous:
          // At this point, dynJacobianColsNbr contains the number of dynamic endogenous
          if (jacobianExo)
            dyn_jacobian_cols_table[deriv_id] = dynJacobianColsNbr + tsid;
          break;
        case SymbolType::exogenousDet:
          // At this point, dynJacobianColsNbr contains the number of dynamic endogenous
          if (jacobianExo)
            dyn_jacobian_cols_table[deriv_id] = dynJacobianColsNbr + symbol_table.exo_nbr() + tsid;
          break;
        case SymbolType::parameter:
        case SymbolType::trend:
        case SymbolType::logTrend:
          // We don't assign a dynamic jacobian column to parameters or trend variables
          break;
        default:
          // Shut up GCC
          cerr << "DynamicModel::computeDynJacobianCols: impossible case" << endl;
          exit(EXIT_FAILURE);
        }
    }

  // Fill in dynamic jacobian columns for endogenous
  int sorted_id = 0;
  for (map<pair<int, int>, int>::const_iterator it = ordered_dyn_endo.begin();
       it != ordered_dyn_endo.end(); it++)
    dyn_jacobian_cols_table[it->second] = sorted_id++;

  // Set final value for dynJacobianColsNbr
  if (jacobianExo)
    dynJacobianColsNbr += symbol_table.exo_nbr() + symbol_table.exo_det_nbr();
}

int
DynamicModel::getDynJacobianCol(int deriv_id) const noexcept(false)
{
  auto it = dyn_jacobian_cols_table.find(deriv_id);
  if (it == dyn_jacobian_cols_table.end())
    throw UnknownDerivIDException();
  else
    return it->second;
}

void
DynamicModel::testTrendDerivativesEqualToZero(const eval_context_t &eval_context)
{
  for (deriv_id_table_t::const_iterator it = deriv_id_table.begin();
       it != deriv_id_table.end(); it++)
    if (symbol_table.getType(it->first.first) == SymbolType::trend
        || symbol_table.getType(it->first.first) == SymbolType::logTrend)
      for (int eq = 0; eq < (int) equations.size(); eq++)
        {
          expr_t homogeneq = AddMinus(equations[eq]->get_arg1(),
                                      equations[eq]->get_arg2());

          // Do not run the test if the term inside the log is zero
          if (fabs(homogeneq->eval(eval_context)) > zero_band)
            {
              expr_t testeq = AddLog(homogeneq); // F = log(lhs-rhs)
              testeq = testeq->getDerivative(it->second); // d F / d Trend
              for (deriv_id_table_t::const_iterator endogit = deriv_id_table.begin();
                   endogit != deriv_id_table.end(); endogit++)
                if (symbol_table.getType(endogit->first.first) == SymbolType::endogenous)
                  {
                    double nearZero = testeq->getDerivative(endogit->second)->eval(eval_context); // eval d F / d Trend d Endog
                    if (fabs(nearZero) > zero_band)
                      {
                        cerr << "WARNING: trends not compatible with balanced growth path; the second-order cross partial of equation " << eq + 1 << " (line "
                             << equations_lineno[eq] << ") w.r.t. trend variable "
                             << symbol_table.getName(it->first.first) << " and endogenous variable "
                             << symbol_table.getName(endogit->first.first) << " is not null. " << endl;
                        // Changed to warning. See discussion in #1389
                      }
                  }
            }
        }
}

void
DynamicModel::writeParamsDerivativesFile(const string &basename, bool julia) const
{
  if (!params_derivatives.size())
    return;

  ExprNodeOutputType output_type = (julia ? ExprNodeOutputType::juliaDynamicModel : ExprNodeOutputType::matlabDynamicModel);
  ostringstream tt_output;              // Used for storing model temp vars and equations
  ostringstream jacobian_output;           // Used for storing jacobian equations
  ostringstream hessian_output;            // Used for storing Hessian equations
  ostringstream hessian1_output;           // Used for storing Hessian equations
  ostringstream third_derivs_output;       // Used for storing third order derivatives equations
  ostringstream third_derivs1_output;      // Used for storing third order derivatives equations

  temporary_terms_t temp_term_union;
  deriv_node_temp_terms_t tef_terms;

  writeModelLocalVariableTemporaryTerms(temp_term_union, params_derivs_temporary_terms_idxs, tt_output, output_type, tef_terms);
  for (auto it : { make_pair(0,1), make_pair(1,1), make_pair(0,2), make_pair(1,2), make_pair(2,1) })
    writeTemporaryTerms(params_derivs_temporary_terms.find(it)->second, temp_term_union, params_derivs_temporary_terms_idxs, tt_output, output_type, tef_terms);

  for (const auto & residuals_params_derivative : params_derivatives.find({ 0, 1 })->second)
    {
      int eq, param;
      tie(eq, param) = vectorToTuple<2>(residuals_params_derivative.first);
      expr_t d1 = residuals_params_derivative.second;

      int param_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param)) + 1;

      jacobian_output << "rp" << LEFT_ARRAY_SUBSCRIPT(output_type) << eq+1 << ", " << param_col
                      << RIGHT_ARRAY_SUBSCRIPT(output_type) << " = ";
      d1->writeOutput(jacobian_output, output_type, temp_term_union, params_derivs_temporary_terms_idxs, tef_terms);
      jacobian_output << ";" << endl;
    }

  for (const auto & jacobian_params_derivative : params_derivatives.find({ 1, 1 })->second)
    {
      int eq, var, param;
      tie(eq, var, param) = vectorToTuple<3>(jacobian_params_derivative.first);
      expr_t d2 = jacobian_params_derivative.second;

      int var_col = getDynJacobianCol(var) + 1;
      int param_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param)) + 1;

      hessian_output << "gp" << LEFT_ARRAY_SUBSCRIPT(output_type) << eq+1 << ", " << var_col
                     << ", " << param_col << RIGHT_ARRAY_SUBSCRIPT(output_type) << " = ";
      d2->writeOutput(hessian_output, output_type, temp_term_union, params_derivs_temporary_terms_idxs, tef_terms);
      hessian_output << ";" << endl;
    }

  int i = 1;
  for (const auto &it : params_derivatives.find({ 0, 2 })->second)
    {
      int eq, param1, param2;
      tie(eq, param1, param2) = vectorToTuple<3>(it.first);
      expr_t d2 = it.second;

      int param1_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param1)) + 1;
      int param2_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param2)) + 1;

      hessian1_output << "rpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",1"
                      << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << eq+1 << ";" << endl
                      << "rpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",2"
                      << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << param1_col << ";" << endl
                      << "rpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",3"
                      << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << param2_col << ";" << endl
                      << "rpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",4"
                      << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=";
      d2->writeOutput(hessian1_output, output_type, temp_term_union, params_derivs_temporary_terms_idxs, tef_terms);
      hessian1_output << ";" << endl;

      i++;

      if (param1 != param2)
        {
          // Treat symmetric elements
          hessian1_output << "rpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",1"
                          << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << eq+1 << ";" << endl
                          << "rpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",2"
                          << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << param2_col << ";" << endl
                          << "rpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",3"
                          << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << param1_col << ";" << endl
                          << "rpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",4"
                          << RIGHT_ARRAY_SUBSCRIPT(output_type)
                          << "=rpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i-1 << ",4"
                          << RIGHT_ARRAY_SUBSCRIPT(output_type) << ";" << endl;
          i++;
        }
    }

  i = 1;
  for (const auto &it : params_derivatives.find({ 1, 2 })->second)
    {
      int eq, var, param1, param2;
      tie(eq, var, param1, param2) = vectorToTuple<4>(it.first);
      expr_t d2 = it.second;

      int var_col = getDynJacobianCol(var) + 1;
      int param1_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param1)) + 1;
      int param2_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param2)) + 1;

      third_derivs_output << "gpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",1"
                          << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << eq+1 << ";" << endl
                          << "gpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",2"
                          << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << var_col << ";" << endl
                          << "gpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",3"
                          << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << param1_col << ";" << endl
                          << "gpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",4"
                          << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << param2_col << ";" << endl
                          << "gpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",5"
                          << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=";
      d2->writeOutput(third_derivs_output, output_type, temp_term_union, params_derivs_temporary_terms_idxs, tef_terms);
      third_derivs_output << ";" << endl;

      i++;

      if (param1 != param2)
        {
          // Treat symmetric elements
          third_derivs_output << "gpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",1"
                              << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << eq+1 << ";" << endl
                              << "gpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",2"
                              << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << var_col << ";" << endl
                              << "gpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",3"
                              << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << param2_col << ";" << endl
                              << "gpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",4"
                              << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << param1_col << ";" << endl
                              << "gpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",5"
                              << RIGHT_ARRAY_SUBSCRIPT(output_type)
                              << "=gpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i-1 << ",5"
                              << RIGHT_ARRAY_SUBSCRIPT(output_type) << ";" << endl;
          i++;
        }
    }

  i = 1;
  for (const auto &it : params_derivatives.find({ 2, 1 })->second)
    {
      int eq, var1, var2, param;
      tie(eq, var1, var2, param) = vectorToTuple<4>(it.first);
      expr_t d2 = it.second;

      int var1_col = getDynJacobianCol(var1) + 1;
      int var2_col = getDynJacobianCol(var2) + 1;
      int param_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param)) + 1;

      third_derivs1_output << "hp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",1"
                           << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << eq+1 << ";" << endl
                           << "hp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",2"
                           << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << var1_col << ";" << endl
                           << "hp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",3"
                           << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << var2_col << ";" << endl
                           << "hp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",4"
                           << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << param_col << ";" << endl
                           << "hp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",5"
                           << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=";
      d2->writeOutput(third_derivs1_output, output_type, temp_term_union, params_derivs_temporary_terms_idxs, tef_terms);
      third_derivs1_output << ";" << endl;

      i++;

      if (var1 != var2)
        {
          // Treat symmetric elements
          third_derivs1_output << "hp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",1"
                               << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << eq+1 << ";" << endl
                               << "hp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",2"
                               << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << var2_col << ";" << endl
                               << "hp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",3"
                               << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << var1_col << ";" << endl
                               << "hp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",4"
                               << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << param_col << ";" << endl
                               << "hp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",5"
                               << RIGHT_ARRAY_SUBSCRIPT(output_type)
                               << "=hp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i-1 << ",5"
                               << RIGHT_ARRAY_SUBSCRIPT(output_type) << ";" << endl;
          i++;
        }
    }

  string filename = julia ? basename + "DynamicParamsDerivs.jl" : packageDir(basename) + "/dynamic_params_derivs.m";
  ofstream paramsDerivsFile;
  paramsDerivsFile.open(filename, ios::out | ios::binary);
  if (!paramsDerivsFile.is_open())
    {
      cerr << "ERROR: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }

  if (!julia)
    {
      // Check that we don't have more than 32 nested parenthesis because Matlab does not suppor this. See Issue #1201
      map<string, string> tmp_paren_vars;
      bool message_printed = false;
      fixNestedParenthesis(tt_output, tmp_paren_vars, message_printed);
      fixNestedParenthesis(jacobian_output, tmp_paren_vars, message_printed);
      fixNestedParenthesis(hessian_output, tmp_paren_vars, message_printed);
      fixNestedParenthesis(hessian1_output, tmp_paren_vars, message_printed);
      fixNestedParenthesis(third_derivs_output, tmp_paren_vars, message_printed);
      fixNestedParenthesis(third_derivs1_output, tmp_paren_vars, message_printed);
      paramsDerivsFile << "function [rp, gp, rpp, gpp, hp] = dynamic_params_derivs(y, x, params, steady_state, it_, ss_param_deriv, ss_param_2nd_deriv)" << endl
                       << "%" << endl
                       << "% Compute the derivatives of the dynamic model with respect to the parameters" << endl
                       << "% Inputs :" << endl
                       << "%   y         [#dynamic variables by 1] double    vector of endogenous variables in the order stored" << endl
                       << "%                                                 in M_.lead_lag_incidence; see the Manual" << endl
                       << "%   x         [nperiods by M_.exo_nbr] double     matrix of exogenous variables (in declaration order)" << endl
                       << "%                                                 for all simulation periods" << endl
                       << "%   params    [M_.param_nbr by 1] double          vector of parameter values in declaration order" << endl
                       << "%   steady_state  [M_.endo_nbr by 1] double       vector of steady state values" << endl
                       << "%   it_       scalar double                       time period for exogenous variables for which to evaluate the model" << endl
                       << "%   ss_param_deriv     [M_.eq_nbr by #params]     Jacobian matrix of the steady states values with respect to the parameters" << endl
                       << "%   ss_param_2nd_deriv [M_.eq_nbr by #params by #params] Hessian matrix of the steady states values with respect to the parameters" << endl
                       << "%" << endl
                       << "% Outputs:" << endl
                       << "%   rp        [M_.eq_nbr by #params] double    Jacobian matrix of dynamic model equations with respect to parameters " << endl
                       << "%                                              Dynare may prepend or append auxiliary equations, see M_.aux_vars" << endl
                       << "%   gp        [M_.endo_nbr by #dynamic variables by #params] double    Derivative of the Jacobian matrix of the dynamic model equations with respect to the parameters" << endl
                       << "%                                                           rows: equations in order of declaration" << endl
                       << "%                                                           columns: variables in order stored in M_.lead_lag_incidence" << endl
                       << "%   rpp       [#second_order_residual_terms by 4] double   Hessian matrix of second derivatives of residuals with respect to parameters;" << endl
                       << "%                                                              rows: respective derivative term" << endl
                       << "%                                                              1st column: equation number of the term appearing" << endl
                       << "%                                                              2nd column: number of the first parameter in derivative" << endl
                       << "%                                                              3rd column: number of the second parameter in derivative" << endl
                       << "%                                                              4th column: value of the Hessian term" << endl
                       << "%   gpp      [#second_order_Jacobian_terms by 5] double   Hessian matrix of second derivatives of the Jacobian with respect to the parameters;" << endl
                       << "%                                                              rows: respective derivative term" << endl
                       << "%                                                              1st column: equation number of the term appearing" << endl
                       << "%                                                              2nd column: column number of variable in Jacobian of the dynamic model" << endl
                       << "%                                                              3rd column: number of the first parameter in derivative" << endl
                       << "%                                                              4th column: number of the second parameter in derivative" << endl
                       << "%                                                              5th column: value of the Hessian term" << endl
                       << "%   hp      [#first_order_Hessian_terms by 5] double   Jacobian matrix of derivatives of the dynamic Hessian with respect to the parameters;" << endl
                       << "%                                                              rows: respective derivative term" << endl
                       << "%                                                              1st column: equation number of the term appearing" << endl
                       << "%                                                              2nd column: column number of first variable in Hessian of the dynamic model" << endl
                       << "%                                                              3rd column: column number of second variable in Hessian of the dynamic model" << endl
                       << "%                                                              4th column: number of the parameter in derivative" << endl
                       << "%                                                              5th column: value of the Hessian term" << endl
                       << "%" << endl
                       << "%" << endl
                       << "% Warning : this file is generated automatically by Dynare" << endl
                       << "%           from model file (.mod)" << endl << endl
                       << "T = NaN(" << params_derivs_temporary_terms_idxs.size() << ",1);" << endl
                       << tt_output.str()
                       << "rp = zeros(" << equations.size() << ", "
                       << symbol_table.param_nbr() << ");" << endl
                       << jacobian_output.str()
                       << "gp = zeros(" << equations.size() << ", " << dynJacobianColsNbr << ", " << symbol_table.param_nbr() << ");" << endl
                       << hessian_output.str()
                       << "if nargout >= 3" << endl
                       << "rpp = zeros(" << params_derivatives.find({ 0, 2 })->second.size() << ",4);" << endl
                       << hessian1_output.str()
                       << "gpp = zeros(" << params_derivatives.find({ 1, 2 })->second.size() << ",5);" << endl
                       << third_derivs_output.str()
                       << "end" << endl
                       << "if nargout >= 5" << endl
                       << "hp = zeros(" << params_derivatives.find({ 2, 1 })->second.size() << ",5);" << endl
                       << third_derivs1_output.str()
                       << "end" << endl
                       << "end" << endl;
    }
  else
    paramsDerivsFile << "module " << basename << "DynamicParamsDerivs" << endl
                     << "#" << endl
                     << "# NB: this file was automatically generated by Dynare" << endl
                     << "#     from " << basename << ".mod" << endl
                     << "#" << endl
                     << "export params_derivs" << endl << endl
                     << "function params_derivs(y, x, paramssteady_state, it_, "
                     << "ss_param_deriv, ss_param_2nd_deriv)" << endl
                     << tt_output.str()
                     << "rp = zeros(" << equations.size() << ", "
                     << symbol_table.param_nbr() << ");" << endl
                     << jacobian_output.str()
                     << "gp = zeros(" << equations.size() << ", " << dynJacobianColsNbr << ", " << symbol_table.param_nbr() << ");" << endl
                     << hessian_output.str()
                     << "rpp = zeros(" << params_derivatives.find({ 0, 2 })->second.size() << ",4);" << endl
                     << hessian1_output.str()
                     << "gpp = zeros(" << params_derivatives.find({ 1, 2 })->second.size() << ",5);" << endl
                     << third_derivs_output.str()
                     << "hp = zeros(" << params_derivatives.find({ 2, 1 })->second.size() << ",5);" << endl
                     << third_derivs1_output.str()
                     << "(rp, gp, rpp, gpp, hp)" << endl
                     << "end" << endl
                     << "end" << endl;

  paramsDerivsFile.close();
}

void
DynamicModel::writeLatexFile(const string &basename, const bool write_equation_tags) const
{
  writeLatexModelFile(basename + "_dynamic", ExprNodeOutputType::latexDynamicModel, write_equation_tags);
}

void
DynamicModel::writeLatexOriginalFile(const string &basename, const bool write_equation_tags) const
{
  writeLatexModelFile(basename + "_original", ExprNodeOutputType::latexDynamicModel, write_equation_tags);
}

void
DynamicModel::substituteEndoLeadGreaterThanTwo(bool deterministic_model)
{
  substituteLeadLagInternal(AuxVarType::endoLead, deterministic_model, vector<string>());
}

void
DynamicModel::substituteEndoLagGreaterThanTwo(bool deterministic_model)
{
  substituteLeadLagInternal(AuxVarType::endoLag, deterministic_model, vector<string>());
}

void
DynamicModel::substituteExoLead(bool deterministic_model)
{
  substituteLeadLagInternal(AuxVarType::exoLead, deterministic_model, vector<string>());
}

void
DynamicModel::substituteExoLag(bool deterministic_model)
{
  substituteLeadLagInternal(AuxVarType::exoLag, deterministic_model, vector<string>());
}

void
DynamicModel::substituteLeadLagInternal(AuxVarType type, bool deterministic_model, const vector<string> &subset)
{
  ExprNode::subst_table_t subst_table;
  vector<BinaryOpNode *> neweqs;

  // Substitute in used model local variables
  set<int> used_local_vars;
  for (auto & equation : equations)
    equation->collectVariables(SymbolType::modelLocalVariable, used_local_vars);

  for (int used_local_var : used_local_vars)
    {
      const expr_t value = local_variables_table.find(used_local_var)->second;
      expr_t subst;
      switch (type)
        {
        case AuxVarType::endoLead:
          subst = value->substituteEndoLeadGreaterThanTwo(subst_table, neweqs, deterministic_model);
          break;
        case AuxVarType::endoLag:
          subst = value->substituteEndoLagGreaterThanTwo(subst_table, neweqs);
          break;
        case AuxVarType::exoLead:
          subst = value->substituteExoLead(subst_table, neweqs, deterministic_model);
          break;
        case AuxVarType::exoLag:
          subst = value->substituteExoLag(subst_table, neweqs);
          break;
        case AuxVarType::diffForward:
          subst = value->differentiateForwardVars(subset, subst_table, neweqs);
          break;
        default:
          cerr << "DynamicModel::substituteLeadLagInternal: impossible case" << endl;
          exit(EXIT_FAILURE);
        }
      local_variables_table[used_local_var] = subst;
    }

  // Substitute in equations
  for (auto & equation : equations)
    {
      expr_t subst;
      switch (type)
        {
        case AuxVarType::endoLead:
          subst = equation->substituteEndoLeadGreaterThanTwo(subst_table, neweqs, deterministic_model);
          break;
        case AuxVarType::endoLag:
          subst = equation->substituteEndoLagGreaterThanTwo(subst_table, neweqs);
          break;
        case AuxVarType::exoLead:
          subst = equation->substituteExoLead(subst_table, neweqs, deterministic_model);
          break;
        case AuxVarType::exoLag:
          subst = equation->substituteExoLag(subst_table, neweqs);
          break;
        case AuxVarType::diffForward:
          subst = equation->differentiateForwardVars(subset, subst_table, neweqs);
          break;
        default:
          cerr << "DynamicModel::substituteLeadLagInternal: impossible case" << endl;
          exit(EXIT_FAILURE);
        }
      auto *substeq = dynamic_cast<BinaryOpNode *>(subst);
      assert(substeq != nullptr);
      equation = substeq;
    }

  // Add new equations
  for (auto & neweq : neweqs)
    addEquation(neweq, -1);

  // Order of auxiliary variable definition equations:
  //  - expectation (entered before this function is called)
  //  - lead variables from lower lead to higher lead
  //  - lag variables from lower lag to higher lag
  copy(neweqs.begin(), neweqs.end(), back_inserter(aux_equations));

  if (neweqs.size() > 0)
    {
      cout << "Substitution of ";
      switch (type)
        {
        case AuxVarType::endoLead:
          cout << "endo leads >= 2";
          break;
        case AuxVarType::endoLag:
          cout << "endo lags >= 2";
          break;
        case AuxVarType::exoLead:
          cout << "exo leads";
          break;
        case AuxVarType::exoLag:
          cout << "exo lags";
          break;
        case AuxVarType::expectation:
          cout << "expectation";
          break;
        case AuxVarType::diffForward:
          cout << "forward vars";
          break;
        default:
          cerr << "DynamicModel::substituteLeadLagInternal: impossible case" << endl;
          exit(EXIT_FAILURE);
        }
      cout << ": added " << neweqs.size() << " auxiliary variables and equations." << endl;
    }
}

void
DynamicModel::substituteAdl()
{
  for (auto & equation : equations)
    equation = dynamic_cast<BinaryOpNode *>(equation->substituteAdl());
}

void
DynamicModel::getEquationNumbersFromTags(vector<int> &eqnumbers, set<string> &eqtags) const
{
  for (auto & eqtag : eqtags)
    for (const auto & equation_tag : equation_tags)
      if (equation_tag.second.first == "name"
          && equation_tag.second.second == eqtag)
        {
          eqnumbers.push_back(equation_tag.first);
          break;
        }
}

void
DynamicModel::findPacExpectationEquationNumbers(vector<int> &eqnumbers) const
{
  int i = 0;
  for (auto & equation : equations)
    {
      if (equation->containsPacExpectation())
        if (find(eqnumbers.begin(), eqnumbers.end(), i) == eqnumbers.end())
          eqnumbers.push_back(i);
      i++;
    }
}

void
DynamicModel::substituteUnaryOps(StaticModel &static_model)
{
  vector<int> eqnumbers(equations.size());
  iota(eqnumbers.begin(), eqnumbers.end(), 0);
  substituteUnaryOps(static_model, eqnumbers);
}

void
DynamicModel::substituteUnaryOps(StaticModel &static_model, set<string> &var_model_eqtags)
{
  vector<int> eqnumbers;
  getEquationNumbersFromTags(eqnumbers, var_model_eqtags);
  findPacExpectationEquationNumbers(eqnumbers);
  substituteUnaryOps(static_model, eqnumbers);
}

void
DynamicModel::substituteUnaryOps(StaticModel &static_model, vector<int> &eqnumbers)
{
  diff_table_t nodes;

  // Find matching unary ops that may be outside of diffs (i.e., those with different lags)
  set<int> used_local_vars;
  for (int eqnumber : eqnumbers)
    equations[eqnumber]->collectVariables(SymbolType::modelLocalVariable, used_local_vars);

  // Only substitute unary ops in model local variables that appear in VAR equations
  for (auto & it : local_variables_table)
    if (used_local_vars.find(it.first) != used_local_vars.end())
      it.second->findUnaryOpNodesForAuxVarCreation(static_model, nodes);

  for (int eqnumber : eqnumbers)
    equations[eqnumber]->findUnaryOpNodesForAuxVarCreation(static_model, nodes);

  // Substitute in model local variables
  ExprNode::subst_table_t subst_table;
  vector<BinaryOpNode *> neweqs;
  for (auto & it : local_variables_table)
    it.second = it.second->substituteUnaryOpNodes(static_model, nodes, subst_table, neweqs);

  // Substitute in equations
  for (auto & equation : equations)
    {
      auto *substeq = dynamic_cast<BinaryOpNode *>(equation->
                                                   substituteUnaryOpNodes(static_model, nodes, subst_table, neweqs));
      assert(substeq != nullptr);
      equation = substeq;
    }

  // Add new equations
  for (auto & neweq : neweqs)
    addEquation(neweq, -1);

  copy(neweqs.begin(), neweqs.end(), back_inserter(aux_equations));

  if (subst_table.size() > 0)
    cout << "Substitution of Unary Ops: added " << neweqs.size() << " auxiliary variables and equations." << endl;
}

void
DynamicModel::substituteDiff(StaticModel &static_model, ExprNode::subst_table_t &diff_subst_table)
{
  set<int> used_local_vars;
  for (const auto & equation : equations)
    equation->collectVariables(SymbolType::modelLocalVariable, used_local_vars);

  // Only substitute diffs in model local variables that appear in VAR equations
  diff_table_t diff_table;
  for (auto & it : local_variables_table)
    if (used_local_vars.find(it.first) != used_local_vars.end())
      it.second->findDiffNodes(static_model, diff_table);

  for (const auto & equation : equations)
    equation->findDiffNodes(static_model, diff_table);

  /* Ensure that all diff operators appear once with their argument at current
     period (i.e. maxLag=0).
     If it is not the case, generate the corresponding expressions.
     This is necessary to avoid lags of more than one in the auxiliary
     equation, which would then be modified by subsequent transformations
     (removing lags > 1), which in turn would break the recursive ordering
     of auxiliary equations. See McModelTeam/McModelProject/issues/95 */
  for (auto &it : diff_table)
    {
      auto iterator_arg_max_lag = it.second.rbegin();
      int arg_max_lag = iterator_arg_max_lag->first;
      expr_t arg_max_expr = iterator_arg_max_lag->second;

      /* We compare arg_max_lag with the result of countDiffs(), in order to
         properly handle nested diffs. See McModelTeam/McModelProject/issues/97 */
      while (arg_max_lag < 1 - it.first->countDiffs())
        {
          arg_max_lag++;
          arg_max_expr = arg_max_expr->decreaseLeadsLags(-1);
          it.second[arg_max_lag] = arg_max_expr;
        }
    }

  // Substitute in model local variables
  vector<BinaryOpNode *> neweqs;
  for (auto & it : local_variables_table)
    it.second = it.second->substituteDiff(static_model, diff_table, diff_subst_table, neweqs);

  // Substitute in equations
  for (auto & equation : equations)
    {
      auto *substeq = dynamic_cast<BinaryOpNode *>(equation->
                                                   substituteDiff(static_model, diff_table, diff_subst_table, neweqs));
      assert(substeq != nullptr);
      equation = substeq;
    }

  // Add new equations
  for (auto & neweq : neweqs)
    addEquation(neweq, -1);

  copy(neweqs.begin(), neweqs.end(), back_inserter(aux_equations));

  if (diff_subst_table.size() > 0)
    cout << "Substitution of Diff operator: added " << neweqs.size() << " auxiliary variables and equations." << endl;
}

void
DynamicModel::substituteExpectation(bool partial_information_model)
{
  ExprNode::subst_table_t subst_table;
  vector<BinaryOpNode *> neweqs;

  // Substitute in model local variables
  for (auto & it : local_variables_table)
    it.second = it.second->substituteExpectation(subst_table, neweqs, partial_information_model);

  // Substitute in equations
  for (auto & equation : equations)
    {
      auto *substeq = dynamic_cast<BinaryOpNode *>(equation->substituteExpectation(subst_table, neweqs, partial_information_model));
      assert(substeq != nullptr);
      equation = substeq;
    }

  // Add new equations
  for (auto & neweq : neweqs)
    addEquation(neweq, -1);

  // Add the new set of equations at the *beginning* of aux_equations
  copy(neweqs.rbegin(), neweqs.rend(), front_inserter(aux_equations));

  if (subst_table.size() > 0)
    {
      if (partial_information_model)
        cout << "Substitution of Expectation operator: added " << subst_table.size() << " auxiliary variables and " << neweqs.size() << " auxiliary equations." << endl;
      else
        cout << "Substitution of Expectation operator: added " << neweqs.size() << " auxiliary variables and equations." << endl;
    }
}

void
DynamicModel::transformPredeterminedVariables()
{
  for (auto & it : local_variables_table)
    it.second = it.second->decreaseLeadsLagsPredeterminedVariables();

  for (auto & equation : equations)
    {
      auto *substeq = dynamic_cast<BinaryOpNode *>(equation->decreaseLeadsLagsPredeterminedVariables());
      assert(substeq != nullptr);
      equation = substeq;
    }
}

void
DynamicModel::detrendEquations()
{
  // We go backwards in the list of trend_vars, to deal correctly with I(2) processes
  for (nonstationary_symbols_map_t::const_reverse_iterator it = nonstationary_symbols_map.rbegin();
       it != nonstationary_symbols_map.rend(); ++it)
    for (auto & equation : equations)
      {
        auto *substeq = dynamic_cast<BinaryOpNode *>(equation->detrend(it->first, it->second.first, it->second.second));
        assert(substeq != nullptr);
        equation = dynamic_cast<BinaryOpNode *>(substeq);
      }

  for (auto & equation : equations)
    {
      BinaryOpNode *substeq = dynamic_cast<BinaryOpNode *>(equation->removeTrendLeadLag(trend_symbols_map));
      assert(substeq != nullptr);
      equation = dynamic_cast<BinaryOpNode *>(substeq);
    }
}

void
DynamicModel::removeTrendVariableFromEquations()
{
  for (auto & equation : equations)
    {
      auto *substeq = dynamic_cast<BinaryOpNode *>(equation->replaceTrendVar());
      assert(substeq != nullptr);
      equation = dynamic_cast<BinaryOpNode *>(substeq);
    }
}

void
DynamicModel::differentiateForwardVars(const vector<string> &subset)
{
  substituteLeadLagInternal(AuxVarType::diffForward, true, subset);
}

void
DynamicModel::fillEvalContext(eval_context_t &eval_context) const
{
  // First, auxiliary variables
  for (auto aux_equation : aux_equations)
    {
      assert(aux_equation->get_op_code() == BinaryOpcode::equal);
      auto *auxvar = dynamic_cast<VariableNode *>(aux_equation->get_arg1());
      assert(auxvar != nullptr);
      try
        {
          double val = aux_equation->get_arg2()->eval(eval_context);
          eval_context[auxvar->get_symb_id()] = val;
        }
      catch (ExprNode::EvalException &e)
        {
          // Do nothing
        }
    }

  // Second, model local variables
  for (auto it : local_variables_table)
    {
      try
        {
          const expr_t expression = it.second;
          double val = expression->eval(eval_context);
          eval_context[it.first] = val;
        }
      catch (ExprNode::EvalException &e)
        {
          // Do nothing
        }
    }

  //Third, trend variables
  vector <int> trendVars = symbol_table.getTrendVarIds();
  for (vector <int>::const_iterator it = trendVars.begin();
       it != trendVars.end(); it++)
    eval_context[*it] = 2;                               //not <= 0 bc of log, not 1 bc of powers
}

bool
DynamicModel::isModelLocalVariableUsed() const
{
  set<int> used_local_vars;
  size_t i = 0;
  while (i < equations.size() && used_local_vars.size() == 0)
    {
      equations[i]->collectVariables(SymbolType::modelLocalVariable, used_local_vars);
      i++;
    }
  return used_local_vars.size() > 0;
}

void
DynamicModel::addStaticOnlyEquation(expr_t eq, int lineno, const vector<pair<string, string>> &eq_tags)
{
  auto *beq = dynamic_cast<BinaryOpNode *>(eq);
  assert(beq != nullptr && beq->get_op_code() == BinaryOpcode::equal);

  vector<pair<string, string>> soe_eq_tags;
  for (const auto & eq_tag : eq_tags)
    soe_eq_tags.push_back(eq_tag);

  static_only_equations.push_back(beq);
  static_only_equations_lineno.push_back(lineno);
  static_only_equations_equation_tags.push_back(soe_eq_tags);
}

size_t
DynamicModel::staticOnlyEquationsNbr() const
{
  return static_only_equations.size();
}

size_t
DynamicModel::dynamicOnlyEquationsNbr() const
{
  set<int> eqs;

  for (const auto & equation_tag : equation_tags)
    if (equation_tag.second.first == "dynamic")
      eqs.insert(equation_tag.first);

  return eqs.size();
}

bool
DynamicModel::isChecksumMatching(const string &basename, bool block) const
{
  boost::crc_32_type result;

  std::stringstream buffer;

  // Write equation tags
  for (const auto & equation_tag : equation_tags)
    buffer << "  " << equation_tag.first + 1
           << equation_tag.second.first
           << equation_tag.second.second;

  ExprNodeOutputType buffer_type = block ? ExprNodeOutputType::matlabDynamicModelSparse : ExprNodeOutputType::CDynamicModel;

  for (int eq = 0; eq < (int) equations.size(); eq++)
    {
      BinaryOpNode *eq_node = equations[eq];
      expr_t lhs = eq_node->get_arg1();
      expr_t rhs = eq_node->get_arg2();

      // Test if the right hand side of the equation is empty.
      double vrhs = 1.0;
      try
        {
          vrhs = rhs->eval(eval_context_t());
        }
      catch (ExprNode::EvalException &e)
        {
        }

      if (vrhs != 0) // The right hand side of the equation is not empty ==> residual=lhs-rhs;
        {
          buffer << "lhs =";
          lhs->writeOutput(buffer, buffer_type, temporary_terms, temporary_terms_idxs);
          buffer << ";" << endl;

          buffer << "rhs =";
          rhs->writeOutput(buffer, buffer_type, temporary_terms, temporary_terms_idxs);
          buffer << ";" << endl;

          buffer << "residual" << LEFT_ARRAY_SUBSCRIPT(buffer_type)
                 << eq + ARRAY_SUBSCRIPT_OFFSET(buffer_type)
                 << RIGHT_ARRAY_SUBSCRIPT(buffer_type)
                 << "= lhs-rhs;" << endl;
        }
      else // The right hand side of the equation is empty ==> residual=lhs;
        {
          buffer << "residual" << LEFT_ARRAY_SUBSCRIPT(buffer_type)
                 << eq + ARRAY_SUBSCRIPT_OFFSET(buffer_type)
                 << RIGHT_ARRAY_SUBSCRIPT(buffer_type)
                 << " = ";
          lhs->writeOutput(buffer, buffer_type, temporary_terms, temporary_terms_idxs);
          buffer << ";" << endl;
        }
    }

  const size_t private_buffer_size{1024};
  char private_buffer[private_buffer_size];
  while (buffer)
    {
      buffer.get(private_buffer, private_buffer_size);
      result.process_bytes(private_buffer, strlen(private_buffer));
    }

  bool basename_dir_exists = !boost::filesystem::create_directory(basename);

  // check whether basename directory exist. If not, create it.
  // If it does, read old checksum if it exist
  fstream checksum_file;
  string filename = basename + "/checksum";
  unsigned int old_checksum = 0;
  // read old checksum if it exists
  if (basename_dir_exists)
    {
      checksum_file.open(filename, ios::in | ios::binary);
      if (checksum_file.is_open())
        {
          checksum_file >> old_checksum;
          checksum_file.close();
        }
    }
  // write new checksum file if none or different from old checksum
  if (old_checksum != result.checksum())
    {
      checksum_file.open(filename, ios::out | ios::binary);
      if (!checksum_file.is_open())
        {
          cerr << "ERROR: Can't open file " << filename << endl;
          exit(EXIT_FAILURE);
        }
      checksum_file << result.checksum();
      checksum_file.close();
      return false;
    }

  return true;
}

void
DynamicModel::writeJsonOutput(ostream &output) const
{
  writeJsonModelEquations(output, false);
  output << ", ";
  writeJsonXrefs(output);
  output << ", ";
  writeJsonAST(output);
}

void
DynamicModel::writeJsonAST(ostream &output) const
{
  vector<pair<string, string>> eqtags;
  output << "\"abstract_syntax_tree\":[" << endl;
  for (int eq = 0; eq < (int) equations.size(); eq++)
    {
      if (eq != 0)
        output << ", ";

      output << "{ \"number\":" << eq
             << ", \"line\":" << equations_lineno[eq];

      for (const auto & equation_tag : equation_tags)
        if (equation_tag.first == eq)
          eqtags.push_back(equation_tag.second);

      if (!eqtags.empty())
        {
          output << ", \"tags\": {";
          int i = 0;
          for (vector<pair<string, string>>::const_iterator it = eqtags.begin(); it != eqtags.end(); it++, i++)
            {
              if (i != 0)
                output << ", ";
              output << "\"" << it->first << "\": \"" << it->second << "\"";
            }
          output << "}";
          eqtags.clear();
        }

      output << ", \"AST\": ";
      equations[eq]->writeJsonAST(output);
      output << "}";
    }
  output << "]";
}

void
DynamicModel::writeJsonXrefsHelper(ostream &output, const map<pair<int, int>, set<int>> &xrefs) const
{
  for (auto it = xrefs.begin();
       it != xrefs.end(); it++)
    {
      if (it != xrefs.begin())
        output << ", ";
      output << "{\"name\": \"" << symbol_table.getName(it->first.first) << "\""
             << ", \"shift\": " << it->first.second
             << ", \"equations\": [";
      for (auto it1 = it->second.begin();
           it1 != it->second.end(); it1++)
        {
          if (it1 != it->second.begin())
            output << ", ";
          output << *it1 + 1;
        }
      output << "]}";
    }
}

void
DynamicModel::writeJsonXrefs(ostream &output) const
{
  output << "\"xrefs\": {"
         << "\"parameters\": [";
  writeJsonXrefsHelper(output, xref_param);
  output << "]"
         << ", \"endogenous\": [";
  writeJsonXrefsHelper(output, xref_endo);
  output << "]"
         << ", \"exogenous\": [";
    writeJsonXrefsHelper(output, xref_exo);
  output << "]"
         << ", \"exogenous_deterministic\": [";
  writeJsonXrefsHelper(output, xref_exo_det);
  output << "]}" << endl;
}

void
DynamicModel::writeJsonOriginalModelOutput(ostream &output) const
{
  writeJsonModelEquations(output, false);
  output << ", ";
  writeJsonAST(output);
}

void
DynamicModel::writeJsonDynamicModelInfo(ostream &output) const
{
  output << "\"model_info\": {"
         << "\"lead_lag_incidence\": [";
  // Loop on endogenous variables
  int nstatic = 0,
    nfwrd   = 0,
    npred   = 0,
    nboth   = 0;
  for (int endoID = 0; endoID < symbol_table.endo_nbr(); endoID++)
    {
      if (endoID != 0)
        output << ",";
      output << "[";
      int sstatic = 1,
        sfwrd   = 0,
        spred   = 0,
        sboth   = 0;
      // Loop on periods
      for (int lag = -max_endo_lag; lag <= max_endo_lead; lag++)
        {
          // Print variableID if exists with current period, otherwise print 0
          try
            {
              if (lag != -max_endo_lag)
                output << ",";
              int varID = getDerivID(symbol_table.getID(SymbolType::endogenous, endoID), lag);
              output << " " << getDynJacobianCol(varID) + 1;
              if (lag == -1)
                {
                  sstatic = 0;
                  spred = 1;
                }
              else if (lag == 1)
                {
                  if (spred == 1)
                    {
                      sboth = 1;
                      spred = 0;
                    }
                  else
                    {
                      sstatic = 0;
                      sfwrd = 1;
                    }
                }
            }
          catch (UnknownDerivIDException &e)
            {
              output << " 0";
            }
        }
      nstatic += sstatic;
      nfwrd   += sfwrd;
      npred   += spred;
      nboth   += sboth;
      output << "]";
    }
  output << "], "
         << "\"nstatic\": " << nstatic << ", "
         << "\"nfwrd\": " << nfwrd << ", "
         << "\"npred\": " << npred << ", "
         << "\"nboth\": " << nboth << ", "
         << "\"nsfwrd\": " << nfwrd+nboth << ", "
         << "\"nspred\": " << npred+nboth << ", "
         << "\"ndynamic\": " << npred+nboth+nfwrd << endl;
  output << "}";
}

void
DynamicModel::writeJsonComputingPassOutput(ostream &output, bool writeDetails) const
{
  ostringstream model_local_vars_output;  // Used for storing model local vars
  ostringstream model_output;             // Used for storing model temp vars and equations
  ostringstream jacobian_output;          // Used for storing jacobian equations
  ostringstream hessian_output;           // Used for storing Hessian equations
  ostringstream third_derivatives_output; // Used for storing third order derivatives equations

  deriv_node_temp_terms_t tef_terms;
  temporary_terms_t temp_term_empty;
  temporary_terms_t temp_term_union = temporary_terms_derivatives[0];
  temporary_terms_t temp_term_union_m_1;

  string concat = "";
  int hessianColsNbr = dynJacobianColsNbr * dynJacobianColsNbr;

  writeJsonModelLocalVariables(model_local_vars_output, tef_terms);

  writeJsonTemporaryTerms(temporary_terms_derivatives[0], temp_term_union_m_1, model_output, tef_terms, concat);
  model_output << ", ";
  writeJsonModelEquations(model_output, true);

  // Writing Jacobian
  temp_term_union_m_1 = temp_term_union;
  temp_term_union.insert(temporary_terms_derivatives[1].begin(), temporary_terms_derivatives[1].end());
  concat = "jacobian";
  writeJsonTemporaryTerms(temp_term_union, temp_term_union_m_1, jacobian_output, tef_terms, concat);
  jacobian_output << ", \"jacobian\": {"
                  << "  \"nrows\": " << equations.size()
                  << ", \"ncols\": " << dynJacobianColsNbr
                  << ", \"entries\": [";
  for (auto it = derivatives[1].begin();
       it != derivatives[1].end(); it++)
    {
      if (it != derivatives[1].begin())
        jacobian_output << ", ";

      int eq, var;
      tie(eq, var) = vectorToTuple<2>(it->first);
      int col =  getDynJacobianCol(var);
      expr_t d1 = it->second;

      if (writeDetails)
        jacobian_output << "{\"eq\": " << eq + 1;
      else
        jacobian_output << "{\"row\": " << eq + 1;

      jacobian_output << ", \"col\": " << col + 1;

      if (writeDetails)
        jacobian_output << ", \"var\": \"" << symbol_table.getName(getSymbIDByDerivID(var)) << "\""
                        << ", \"shift\": " << getLagByDerivID(var);

      jacobian_output << ", \"val\": \"";
      d1->writeJsonOutput(jacobian_output, temp_term_union, tef_terms);
      jacobian_output << "\"}" << endl;
    }
  jacobian_output << "]}";

  // Writing Hessian
  temp_term_union_m_1 = temp_term_union;
  temp_term_union.insert(temporary_terms_derivatives[2].begin(), temporary_terms_derivatives[2].end());
  concat = "hessian";
  writeJsonTemporaryTerms(temp_term_union, temp_term_union_m_1, hessian_output, tef_terms, concat);
  hessian_output << ", \"hessian\": {"
                 << "  \"nrows\": " << equations.size()
                 << ", \"ncols\": " << hessianColsNbr
                 << ", \"entries\": [";
  for (auto it = derivatives[2].begin();
       it != derivatives[2].end(); it++)
    {
      if (it != derivatives[2].begin())
        hessian_output << ", ";

      int eq, var1, var2;
      tie(eq, var1, var2) = vectorToTuple<3>(it->first);
      expr_t d2 = it->second;
      int id1 = getDynJacobianCol(var1);
      int id2 = getDynJacobianCol(var2);
      int col_nb = id1 * dynJacobianColsNbr + id2;
      int col_nb_sym = id2 * dynJacobianColsNbr + id1;

      if (writeDetails)
        hessian_output << "{\"eq\": " << eq + 1;
      else
        hessian_output << "{\"row\": " << eq + 1;

      hessian_output << ", \"col\": [" << col_nb + 1;
      if (id1 != id2)
        hessian_output << ", " << col_nb_sym + 1;
      hessian_output << "]";

      if (writeDetails)
        hessian_output << ", \"var1\": \"" << symbol_table.getName(getSymbIDByDerivID(var1)) << "\""
                       << ", \"shift1\": " << getLagByDerivID(var1)
                       << ", \"var2\": \"" << symbol_table.getName(getSymbIDByDerivID(var2)) << "\""
                       << ", \"shift2\": " << getLagByDerivID(var2);

      hessian_output << ", \"val\": \"";
      d2->writeJsonOutput(hessian_output, temp_term_union, tef_terms);
      hessian_output << "\"}" << endl;
    }
  hessian_output << "]}";

  // Writing third derivatives
  temp_term_union_m_1 = temp_term_union;
  temp_term_union.insert(temporary_terms_derivatives[3].begin(), temporary_terms_derivatives[3].end());
  concat = "third_derivatives";
  writeJsonTemporaryTerms(temp_term_union, temp_term_union_m_1, third_derivatives_output, tef_terms, concat);
  third_derivatives_output << ", \"third_derivative\": {"
                           << "  \"nrows\": " << equations.size()
                           << ", \"ncols\": " << hessianColsNbr * dynJacobianColsNbr
                           << ", \"entries\": [";
  for (auto it = derivatives[3].begin();
       it != derivatives[3].end(); it++)
    {
      if (it != derivatives[3].begin())
        third_derivatives_output << ", ";

      int eq, var1, var2, var3;
      tie(eq, var1, var2, var3) = vectorToTuple<4>(it->first);
      expr_t d3 = it->second;

      if (writeDetails)
        third_derivatives_output << "{\"eq\": " << eq + 1;
      else
        third_derivatives_output << "{\"row\": " << eq + 1;

      int id1 = getDynJacobianCol(var1);
      int id2 = getDynJacobianCol(var2);
      int id3 = getDynJacobianCol(var3);
      set<int> cols;
      cols.insert(id1 * hessianColsNbr + id2 * dynJacobianColsNbr + id3);
      cols.insert(id1 * hessianColsNbr + id3 * dynJacobianColsNbr + id2);
      cols.insert(id2 * hessianColsNbr + id1 * dynJacobianColsNbr + id3);
      cols.insert(id2 * hessianColsNbr + id3 * dynJacobianColsNbr + id1);
      cols.insert(id3 * hessianColsNbr + id1 * dynJacobianColsNbr + id2);
      cols.insert(id3 * hessianColsNbr + id2 * dynJacobianColsNbr + id1);

      third_derivatives_output << ", \"col\": [";
      for (auto it2 = cols.begin(); it2 != cols.end(); it2++)
        {
          if (it2 != cols.begin())
            third_derivatives_output << ", ";
          third_derivatives_output << *it2 + 1;
        }
      third_derivatives_output << "]";

      if (writeDetails)
        third_derivatives_output << ", \"var1\": \"" << symbol_table.getName(getSymbIDByDerivID(var1)) << "\""
                                 << ", \"shift1\": " << getLagByDerivID(var1)
                                 << ", \"var2\": \"" << symbol_table.getName(getSymbIDByDerivID(var2)) << "\""
                                 << ", \"shift2\": " << getLagByDerivID(var2)
                                 << ", \"var3\": \"" << symbol_table.getName(getSymbIDByDerivID(var3)) << "\""
                                 << ", \"shift3\": " << getLagByDerivID(var3);

      third_derivatives_output << ", \"val\": \"";
      d3->writeJsonOutput(third_derivatives_output, temp_term_union, tef_terms);
      third_derivatives_output << "\"}" << endl;
    }
  third_derivatives_output << "]}";

  if (writeDetails)
    output << "\"dynamic_model\": {";
  else
    output << "\"dynamic_model_simple\": {";
  output << model_local_vars_output.str()
         << ", " << model_output.str()
         << ", " << jacobian_output.str()
         << ", " << hessian_output.str()
         << ", " << third_derivatives_output.str()
         << "}";
}

void
DynamicModel::writeJsonParamsDerivativesFile(ostream &output, bool writeDetails) const
{
  if (!params_derivatives.size())
    return;

  ostringstream model_local_vars_output;   // Used for storing model local vars
  ostringstream model_output;              // Used for storing model temp vars and equations
  ostringstream jacobian_output;           // Used for storing jacobian equations
  ostringstream hessian_output;            // Used for storing Hessian equations
  ostringstream hessian1_output;           // Used for storing Hessian equations
  ostringstream third_derivs_output;       // Used for storing third order derivatives equations
  ostringstream third_derivs1_output;      // Used for storing third order derivatives equations

  deriv_node_temp_terms_t tef_terms;
  writeJsonModelLocalVariables(model_local_vars_output, tef_terms);

  temporary_terms_t temp_term_union;
  string concat = "all";
  for (auto it : { make_pair(0,1), make_pair(1,1), make_pair(0,2), make_pair(1,2), make_pair(2,1) })
    writeJsonTemporaryTerms(params_derivs_temporary_terms.find(it)->second, temp_term_union, model_output, tef_terms, concat);

  jacobian_output << "\"deriv_wrt_params\": {"
                  << "  \"neqs\": " << equations.size()
                  << ", \"nparamcols\": " << symbol_table.param_nbr()
                  << ", \"entries\": [";
  auto &rp = params_derivatives.find({ 0, 1 })->second;
  for (auto it = rp.begin(); it != rp.end(); it++)
    {
      if (it != rp.begin())
        jacobian_output << ", ";

      int eq, param;
      tie(eq, param) = vectorToTuple<2>(it->first);
      expr_t d1 = it->second;

      int param_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param)) + 1;

      if (writeDetails)
        jacobian_output << "{\"eq\": " << eq + 1;
      else
        jacobian_output << "{\"row\": " << eq + 1;

      jacobian_output << ", \"param_col\": " << param_col + 1;

      if (writeDetails)
        jacobian_output << ", \"param\": \"" << symbol_table.getName(getSymbIDByDerivID(param)) << "\"";

      jacobian_output << ", \"val\": \"";
      d1->writeJsonOutput(jacobian_output, temp_term_union, tef_terms);
      jacobian_output << "\"}" << endl;
    }
  jacobian_output << "]}";

  hessian_output << "\"deriv_jacobian_wrt_params\": {"
                 << "  \"neqs\": " << equations.size()
                 << ", \"nvarcols\": " << dynJacobianColsNbr
                 << ", \"nparamcols\": " << symbol_table.param_nbr()
                 << ", \"entries\": [";
  auto &gp = params_derivatives.find({ 1, 1 })->second;
  for (auto it = gp.begin(); it != gp.end(); it++)
    {
      if (it != gp.begin())
        hessian_output << ", ";

      int eq, var, param;
      tie(eq, var, param) = vectorToTuple<3>(it->first);
      expr_t d2 = it->second;

      int var_col = getDynJacobianCol(var) + 1;
      int param_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param)) + 1;

      if (writeDetails)
        hessian_output << "{\"eq\": " << eq + 1;
      else
        hessian_output << "{\"row\": " << eq + 1;

      hessian_output << ", \"var_col\": " << var_col + 1
                     << ", \"param_col\": " << param_col + 1;

      if (writeDetails)
      hessian_output << ", \"var\": \"" << symbol_table.getName(getSymbIDByDerivID(var)) << "\""
                     << ", \"lag\": " << getLagByDerivID(var)
                     << ", \"param\": \"" << symbol_table.getName(getSymbIDByDerivID(param)) << "\"";

      hessian_output << ", \"val\": \"";
      d2->writeJsonOutput(hessian_output, temp_term_union, tef_terms);
      hessian_output << "\"}" << endl;
    }
  hessian_output << "]}";

  hessian1_output << "\"second_deriv_residuals_wrt_params\": {"
                  << "  \"nrows\": " << equations.size()
                  << ", \"nparam1cols\": " << symbol_table.param_nbr()
                  << ", \"nparam2cols\": " << symbol_table.param_nbr()
                  << ", \"entries\": [";
  auto &rpp = params_derivatives.find({ 0, 2 })->second;
  for (auto it = rpp.begin(); it != rpp.end(); ++it)
    {
      if (it != rpp.begin())
        hessian1_output << ", ";

      int eq, param1, param2;
      tie(eq, param1, param2) = vectorToTuple<3>(it->first);
      expr_t d2 = it->second;

      int param1_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param1)) + 1;
      int param2_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param2)) + 1;

      if (writeDetails)
        hessian1_output << "{\"eq\": " << eq + 1;
      else
        hessian1_output << "{\"row\": " << eq + 1;
      hessian1_output << ", \"param1_col\": " << param1_col + 1
                      << ", \"param2_col\": " << param2_col + 1;

      if (writeDetails)
        hessian1_output << ", \"param1\": \"" << symbol_table.getName(getSymbIDByDerivID(param1)) << "\""
                        << ", \"param2\": \"" << symbol_table.getName(getSymbIDByDerivID(param2)) << "\"";

      hessian1_output << ", \"val\": \"";
      d2->writeJsonOutput(hessian1_output, temp_term_union, tef_terms);
      hessian1_output << "\"}" << endl;
    }
  hessian1_output << "]}";

  third_derivs_output << "\"second_deriv_jacobian_wrt_params\": {"
                      << "  \"neqs\": " << equations.size()
                      << ", \"nvarcols\": " << dynJacobianColsNbr
                      << ", \"nparam1cols\": " << symbol_table.param_nbr()
                      << ", \"nparam2cols\": " << symbol_table.param_nbr()
                      << ", \"entries\": [";
  auto &gpp = params_derivatives.find({ 1, 2 })->second;
  for (auto it = gpp.begin(); it != gpp.end(); ++it)
    {
      if (it != gpp.begin())
        third_derivs_output << ", ";

      int eq, var, param1, param2;
      tie(eq, var, param1, param2) = vectorToTuple<4>(it->first);
      expr_t d2 = it->second;

      int var_col = getDynJacobianCol(var) + 1;
      int param1_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param1)) + 1;
      int param2_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param2)) + 1;

      if (writeDetails)
        third_derivs_output << "{\"eq\": " << eq + 1;
      else
        third_derivs_output << "{\"row\": " << eq + 1;

      third_derivs_output << ", \"var_col\": " << var_col + 1
                          << ", \"param1_col\": " << param1_col + 1
                          << ", \"param2_col\": " << param2_col + 1;

      if (writeDetails)
        third_derivs_output << ", \"var\": \"" << symbol_table.getName(var) << "\""
                            << ", \"lag\": " << getLagByDerivID(var)
                            << ", \"param1\": \"" << symbol_table.getName(getSymbIDByDerivID(param1)) << "\""
                            << ", \"param2\": \"" << symbol_table.getName(getSymbIDByDerivID(param2)) << "\"";

      third_derivs_output << ", \"val\": \"";
      d2->writeJsonOutput(third_derivs_output, temp_term_union, tef_terms);
      third_derivs_output << "\"}" << endl;
    }
  third_derivs_output << "]}" << endl;

  third_derivs1_output << "\"derivative_hessian_wrt_params\": {"
                       << "  \"neqs\": " << equations.size()
                       << ", \"nvar1cols\": " << dynJacobianColsNbr
                       << ", \"nvar2cols\": " << dynJacobianColsNbr
                       << ", \"nparamcols\": " << symbol_table.param_nbr()
                       << ", \"entries\": [";
  auto &hp = params_derivatives.find({ 2, 1 })->second;
  for (auto it = hp.begin(); it != hp.end(); ++it)
    {
      if (it != hp.begin())
        third_derivs1_output << ", ";

      int eq, var1, var2, param;
      tie(eq, var1, var2, param) = vectorToTuple<4>(it->first);
      expr_t d2 = it->second;

      int var1_col = getDynJacobianCol(var1) + 1;
      int var2_col = getDynJacobianCol(var2) + 1;
      int param_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param)) + 1;

      if (writeDetails)
        third_derivs1_output << "{\"eq\": " << eq + 1;
      else
        third_derivs1_output << "{\"row\": " << eq + 1;

      third_derivs1_output << ", \"var1_col\": " << var1_col + 1
                           << ", \"var2_col\": " << var2_col + 1
                           << ", \"param_col\": " << param_col + 1;

      if (writeDetails)
        third_derivs1_output << ", \"var1\": \"" << symbol_table.getName(getSymbIDByDerivID(var1)) << "\""
                             << ", \"lag1\": " << getLagByDerivID(var1)
                             << ", \"var2\": \"" << symbol_table.getName(getSymbIDByDerivID(var2)) << "\""
                             << ", \"lag2\": " << getLagByDerivID(var2)
                             << ", \"param\": \"" << symbol_table.getName(getSymbIDByDerivID(param)) << "\"";

      third_derivs1_output << ", \"val\": \"";
      d2->writeJsonOutput(third_derivs1_output, temp_term_union, tef_terms);
      third_derivs1_output << "\"}" << endl;
    }
  third_derivs1_output << "]}" << endl;

  if (writeDetails)
    output << "\"dynamic_model_params_derivative\": {";
  else
    output << "\"dynamic_model_params_derivatives_simple\": {";
  output << model_local_vars_output.str()
         << ", " << model_output.str()
         << ", " << jacobian_output.str()
         << ", " << hessian_output.str()
         << ", " << hessian1_output.str()
         << ", " << third_derivs_output.str()
         << ", " << third_derivs1_output.str()
         << "}";
}

void
DynamicModel::substituteVarExpectation(const map<string, expr_t> &subst_table)
{
  for (auto & equation : equations)
    equation = dynamic_cast<BinaryOpNode *>(equation->substituteVarExpectation(subst_table));
}
