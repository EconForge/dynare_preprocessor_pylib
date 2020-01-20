/*
 * Copyright © 2003-2019 Dynare Team
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
#include <algorithm>
#include <numeric>
#include <regex>

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
    equation_type_and_normalized_equation.emplace_back(it.first, f(it.second));

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
}

DynamicModel::DynamicModel(SymbolTable &symbol_table_arg,
                           NumericalConstants &num_constants_arg,
                           ExternalFunctionsTable &external_functions_table_arg,
                           TrendComponentModelTable &trend_component_model_table_arg,
                           VarModelTable &var_model_table_arg) :
  ModelTree{symbol_table_arg, num_constants_arg, external_functions_table_arg, true},
  trend_component_model_table{trend_component_model_table_arg},
  var_model_table{var_model_table_arg}
{
}

DynamicModel::DynamicModel(const DynamicModel &m) :
  ModelTree{m},
  trend_component_model_table{m.trend_component_model_table},
  var_model_table{m.var_model_table},
  balanced_growth_test_tol{m.balanced_growth_test_tol},
  static_only_equations_lineno{m.static_only_equations_lineno},
  static_only_equations_equation_tags{m.static_only_equations_equation_tags},
  static_only_equation_tags_xref{m.static_only_equation_tags_xref},
  deriv_id_table{m.deriv_id_table},
  inv_deriv_id_table{m.inv_deriv_id_table},
  dyn_jacobian_cols_table{m.dyn_jacobian_cols_table},
  max_lag{m.max_lag},
  max_lead{m.max_lead},
  max_endo_lag{m.max_endo_lag},
  max_endo_lead{m.max_endo_lead},
  max_exo_lag{m.max_exo_lag},
  max_exo_lead{m.max_exo_lead},
  max_exo_det_lag{m.max_exo_det_lag},
  max_exo_det_lead{m.max_exo_det_lead},
  max_lag_orig{m.max_lag_orig},
  max_lead_orig{m.max_lead_orig},
  max_lag_with_diffs_expanded_orig{m.max_lag_with_diffs_expanded_orig},
  max_endo_lag_orig{m.max_endo_lag_orig},
  max_endo_lead_orig{m.max_endo_lead_orig},
  max_exo_lag_orig{m.max_exo_lag_orig},
  max_exo_lead_orig{m.max_exo_lead_orig},
  max_exo_det_lag_orig{m.max_exo_det_lag_orig},
  max_exo_det_lead_orig{m.max_exo_det_lead_orig},
  xrefs{m.xrefs},
  xref_param{m.xref_param},
  xref_endo{m.xref_endo},
  xref_exo{m.xref_exo},
  xref_exo_det{m.xref_exo_det},
  nonzero_hessian_eqs{m.nonzero_hessian_eqs},
  v_temporary_terms_inuse{m.v_temporary_terms_inuse},
  map_idx{m.map_idx},
  global_temporary_terms{m.global_temporary_terms},
  block_type_firstequation_size_mfs{m.block_type_firstequation_size_mfs},
  blocks_linear{m.blocks_linear},
  other_endo_block{m.other_endo_block},
  exo_block{m.exo_block},
  exo_det_block{m.exo_det_block},
  block_var_exo{m.block_var_exo},
  block_exo_index{m.block_exo_index},
  block_det_exo_index{m.block_det_exo_index},
  block_other_endo_index{m.block_other_endo_index},
  block_col_type{m.block_col_type},
  variable_block_lead_lag{m.variable_block_lead_lag},
  equation_block{m.equation_block},
  var_expectation_functions_to_write{m.var_expectation_functions_to_write},
  endo_max_leadlag_block{m.endo_max_leadlag_block},
  other_endo_max_leadlag_block{m.other_endo_max_leadlag_block},
  exo_max_leadlag_block{m.exo_max_leadlag_block},
  exo_det_max_leadlag_block{m.exo_det_max_leadlag_block},
  max_leadlag_block{m.max_leadlag_block}
{
  copyHelper(m);
}

DynamicModel &
DynamicModel::operator=(const DynamicModel &m)
{
  ModelTree::operator=(m);

  assert(&trend_component_model_table == &m.trend_component_model_table);
  assert(&var_model_table == &m.var_model_table);
  balanced_growth_test_tol = m.balanced_growth_test_tol;

  static_only_equations_lineno = m.static_only_equations_lineno;
  static_only_equations_equation_tags = m.static_only_equations_equation_tags;
  static_only_equation_tags_xref = m.static_only_equation_tags_xref;
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
  max_lag_with_diffs_expanded_orig = m.max_lag_with_diffs_expanded_orig;
  max_endo_lag_orig = m.max_endo_lag_orig;
  max_endo_lead_orig = m.max_endo_lead_orig;
  max_exo_lag_orig = m.max_exo_lag_orig;
  max_exo_lead_orig = m.max_exo_lead_orig;
  max_exo_det_lag_orig = m.max_exo_det_lag_orig;
  max_exo_det_lead_orig = m.max_exo_det_lead_orig;
  xrefs = m.xrefs;
  xref_param = m.xref_param;
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
  if (auto it = derivatives[1].find({ eq, getDerivID(symbol_table.getID(SymbolType::endogenous, symb_id), lag) });
      it != derivatives[1].end())
    it->second->compile(code_file, instruction_number, false, temporary_terms, map_idx, true, false);
  else
    {
      FLDZ_ fldz;
      fldz.write(code_file, instruction_number);
    }
}

void
DynamicModel::compileChainRuleDerivative(ofstream &code_file, unsigned int &instruction_number, int eqr, int varr, int lag, const map_idx_t &map_idx) const
{
  if (auto it = first_chain_rule_derivatives.find({ eqr, varr, lag }); it != first_chain_rule_derivatives.end())
    it->second->compile(code_file, instruction_number, false, temporary_terms, map_idx, true, false);
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
                getBlockEquationRenormalizedExpr(block, i)->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, block, v_temporary_terms, i);
              else
                {
                  eq_node = static_cast<BinaryOpNode *>(getBlockEquationExpr(block, i));
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
                getBlockEquationRenormalizedExpr(block, i)->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, block, v_temporary_terms, i);
              else
                {
                  eq_node = static_cast<BinaryOpNode *>(getBlockEquationExpr(block, i));
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
                  eq_node = static_cast<BinaryOpNode *>(getBlockEquationExpr(block, i));
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
  ofstream output;
  int nze, nze_exo, nze_exo_det, nze_other_endo;
  vector<int> feedback_variables;
  ExprNodeOutputType local_output_type;

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
      for (const auto &it : tmp_block_exo_det_derivative)
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
                || simulation_type == EVALUATE_BACKWARD || simulation_type == EVALUATE_FORWARD)
               && getBlockFirstEquation(block) < prologue)
        block_type = PROLOGUE;
      else if ((simulation_type == SOLVE_FORWARD_SIMPLE || simulation_type == SOLVE_BACKWARD_SIMPLE
                || simulation_type == EVALUATE_BACKWARD || simulation_type == EVALUATE_FORWARD)
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
          for (int i = 0; i < static_cast<int>(block_size); i++)
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
          eq_node = static_cast<BinaryOpNode *>(getBlockEquationExpr(block, i));
          lhs = eq_node->arg1;
          rhs = eq_node->arg2;
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
                      eq_node = static_cast<BinaryOpNode *>(getBlockEquationRenormalizedExpr(block, i));
                      lhs = eq_node->arg1;
                      rhs = eq_node->arg2;
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
        if (simulation_type == SOLVE_BACKWARD_SIMPLE || simulation_type == SOLVE_FORWARD_SIMPLE
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

  filesystem::create_directories(basename + "/model/bytecode");

  string main_name = basename + "/model/bytecode/dynamic.cod";
  code_file.open(main_name, ios::out | ios::binary | ios::ate);
  if (!code_file.is_open())
    {
      cerr << R"(Error : Can't open file ")" << main_name << R"(" for writing)" << endl;
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
  for (const auto & [indices, d1] : derivatives[1])
    {
      int deriv_id = indices[1];
      unsigned int eq = indices[0];
      int symb = getSymbIDByDerivID(deriv_id);
      unsigned int var = symbol_table.getTypeSpecificID(symb);
      int lag = getLagByDerivID(deriv_id);
      if (getTypeByDerivID(deriv_id) == SymbolType::endogenous)
        first_derivatives_reordered_endo[{ lag, var, eq }] = d1;
      else if (getTypeByDerivID(deriv_id) == SymbolType::exogenous || getTypeByDerivID(deriv_id) == SymbolType::exogenousDet)
        first_derivatives_reordered_exo[{ lag, getTypeByDerivID(deriv_id), var, eq }] = d1;
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

  for (const auto &it : first_derivatives_reordered_exo)
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
  for (const auto & [indices, d1] : derivatives[1])
    {
      int deriv_id = indices[1];
      if (getTypeByDerivID(deriv_id) == SymbolType::endogenous)
        {
          unsigned int eq = indices[0];
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
          for (auto it = my_derivatives[i].begin(); it != my_derivatives[i].end(); ++it)
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
  for (const auto &it : first_derivatives_reordered_exo)
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
DynamicModel::writeModelEquationsCode_Block(const string &basename, const map_idx_t &map_idx, bool linear_decomposition) const
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
  filesystem::create_directories(basename + "/model/bytecode");
  if (linear_decomposition)
    main_name = basename + "/model/bytecode/non_linear.cod";
  else
    main_name = basename + "/model/bytecode/dynamic.cod";
  code_file.open(main_name, ios::out | ios::binary | ios::ate);
  if (!code_file.is_open())
    {
      cerr << R"(Error : Can't open file ")" << main_name << R"(" for writing)" << endl;
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
      for (const auto &it : exo_det_block[block])
        for (const auto &it1 : it.second)
          {
            count_col_det_exo++;
            if (find(exo_det.begin(), exo_det.end(), it1) == exo_det.end())
              exo_det.push_back(it1);
          }

      unsigned int count_col_exo = 0;
      vector<unsigned int> exo;
      for (const auto &it : exo_block[block])
        for (const auto &it1 : it.second)
          {
            count_col_exo++;
            if (find(exo.begin(), exo.end(), it1) == exo.end())
              exo.push_back(it1);
          }

      vector<unsigned int> other_endo;
      unsigned int count_col_other_endo = 0;
      for (const auto &it : other_endo_block[block])
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
      for (i = 0; i < static_cast<int>(block_size); i++)
        {
          //The Temporary terms
          temporary_terms_t tt2;
          if (v_temporary_terms[block][i].size() && !linear_decomposition)
            {
              for (auto it : v_temporary_terms[block][i])
                {
                  if (dynamic_cast<AbstractExternalFunctionNode *>(it) != nullptr)
                    it->compileExternalFunctionOutput(code_file, instruction_number, false, tt2, map_idx, true, false, tef_terms);

                  FNUMEXPR_ fnumexpr(TemporaryTerm, static_cast<int>(map_idx.find(it->idx)->second));
                  fnumexpr.write(code_file, instruction_number);
                  it->compile(code_file, instruction_number, false, tt2, map_idx, true, false, tef_terms);
                  FSTPT_ fstpt(static_cast<int>(map_idx.find(it->idx)->second));
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
                  eq_node = static_cast<BinaryOpNode *>(getBlockEquationExpr(block, i));
                  lhs = eq_node->arg1;
                  rhs = eq_node->arg2;
                  rhs->compile(code_file, instruction_number, false, temporary_terms, map_idx, true, false);
                  lhs->compile(code_file, instruction_number, true, temporary_terms, map_idx, true, false);
                }
              else if (equ_type == E_EVALUATE_S)
                {
                  eq_node = static_cast<BinaryOpNode *>(getBlockEquationRenormalizedExpr(block, i));
                  lhs = eq_node->arg1;
                  rhs = eq_node->arg2;
                  rhs->compile(code_file, instruction_number, false, temporary_terms, map_idx, true, false);
                  lhs->compile(code_file, instruction_number, true, temporary_terms, map_idx, true, false);
                }
              break;
            case SOLVE_BACKWARD_COMPLETE:
            case SOLVE_FORWARD_COMPLETE:
            case SOLVE_TWO_BOUNDARIES_COMPLETE:
            case SOLVE_TWO_BOUNDARIES_SIMPLE:
              if (i < static_cast<int>(block_recursive))
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
              eq_node = static_cast<BinaryOpNode *>(getBlockEquationExpr(block, i));
              lhs = eq_node->arg1;
              rhs = eq_node->arg2;
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
      if (simulation_type != EVALUATE_BACKWARD
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
                          Uf[eqr].Ufl = static_cast<Uff_l *>(malloc(sizeof(Uff_l)));
                          Uf[eqr].Ufl_First = Uf[eqr].Ufl;
                        }
                      else
                        {
                          Uf[eqr].Ufl->pNext = static_cast<Uff_l *>(malloc(sizeof(Uff_l)));
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
              for (i = 0; i < static_cast<int>(block_size); i++)
                {
                  if (i >= static_cast<int>(block_recursive))
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
          FSTPG3_ fstpg3(eq, var, lag, count_col_exo-1);
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
DynamicModel::writeDynamicCFile(const string &basename) const
{
  filesystem::create_directories(basename + "/model/src");
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

  mDynamicModelFile << "#include <stdlib.h>" << endl;

  if (external_functions_table.get_total_number_of_unique_model_block_external_functions())
    // External Matlab function, implies Dynamic function will call mex
    mDynamicModelFile
#ifndef __APPLE__
      << "#include <uchar.h>" << endl // For MATLAB ≤ R2011a
#else
      << "typedef uint_least16_t char16_t;" << endl
      << "typedef uint_least32_t char32_t;" << endl // uchar.h does not exist on macOS
#endif
      << R"(#include "mex.h")" << endl;

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
#ifndef __APPLE__
                  << "#include <uchar.h>" << endl // For MATLAB ≤ R2011a
#else
                  << "typedef uint_least16_t char16_t;" << endl
                  << "typedef uint_least32_t char32_t;" << endl // uchar.h does not exist on macOS
#endif
                  << R"(#include "mex.h")" << endl
                  << endl
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
                  << "  if (nlhs > " << computed_derivs_order << ")" << endl
                  << R"(    mexErrMsgTxt("Derivatives of higher order than computed have been requested");)" << endl
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
                  << "  double *T = (double *) malloc(sizeof(double)*" << ntt << ");"
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
DynamicModel::reform(const string &name1) const
{
  string name = name1;
  int pos = name.find(R"(\)", 0);
  while (pos >= 0)
    {
      if (name.substr(pos + 1, 1) != R"(\)")
        {
          name = name.insert(pos, R"(\)");
          pos++;
        }
      pos++;
      pos = name.find(R"(\)", pos);
    }
  return name;
}

void
DynamicModel::printNonZeroHessianEquations(ostream &output) const
{
  if (nonzero_hessian_eqs.size() != 1)
    output << "[";
  for (auto it = nonzero_hessian_eqs.begin();
       it != nonzero_hessian_eqs.end(); ++it)
    {
      if (it != nonzero_hessian_eqs.begin())
        output << " ";
      output << *it + 1;
    }
  if (nonzero_hessian_eqs.size() != 1)
    output << "]";
}

void
DynamicModel::Write_Inf_To_Bin_File_Block(const string &basename, int num,
                                          int &u_count_int, bool &file_open, bool is_two_boundaries, bool linear_decomposition) const
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
      cerr << R"(Error : Can't open file ")" << filename << R"(" for writing)" << endl;
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
  for (j = block_recursive; j < static_cast<int>(block_size); j++)
    {
      unsigned int varr = getBlockVariableID(num, j);
      SaveCode.write(reinterpret_cast<char *>(&varr), sizeof(varr));
    }
  for (j = block_recursive; j < static_cast<int>(block_size); j++)
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
  bool OK = true;
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
      unsigned int block_size = getBlockSize(block),
        block_mfs = getBlockMfs(block),
        block_recursive = block_size - block_mfs;
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

      if (simulation_type == EVALUATE_FORWARD && block_size)
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
      else if (simulation_type == EVALUATE_BACKWARD && block_size)
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
      else if ((simulation_type == SOLVE_FORWARD_COMPLETE || simulation_type == SOLVE_FORWARD_SIMPLE) && block_size)
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
      else if ((simulation_type == SOLVE_BACKWARD_COMPLETE || simulation_type == SOLVE_BACKWARD_SIMPLE) && block_size)
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
      else if ((simulation_type == SOLVE_TWO_BOUNDARIES_COMPLETE || simulation_type == SOLVE_TWO_BOUNDARIES_SIMPLE) && block_size)
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
    name = "dynamic_resid_g1_g2";
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
  vector<ostringstream> d_output(derivatives.size()); // Derivatives output (at all orders, including 0=residual)
  vector<ostringstream> tt_output(derivatives.size()); // Temp terms output (at all orders)

  ExprNodeOutputType output_type = (use_dll ? ExprNodeOutputType::CDynamicModel :
                                    julia ? ExprNodeOutputType::juliaDynamicModel : ExprNodeOutputType::matlabDynamicModel);

  deriv_node_temp_terms_t tef_terms;
  temporary_terms_t temp_term_union;

  writeModelLocalVariableTemporaryTerms(temp_term_union, temporary_terms_idxs,
                                        tt_output[0], output_type, tef_terms);

  writeTemporaryTerms(temporary_terms_derivatives[0],
                      temp_term_union,
                      temporary_terms_idxs,
                      tt_output[0], output_type, tef_terms);

  writeModelEquations(d_output[0], output_type, temp_term_union);

  int nrows = equations.size();
  int hessianColsNbr = dynJacobianColsNbr * dynJacobianColsNbr;

  // Writing Jacobian
  if (!derivatives[1].empty())
    {
      writeTemporaryTerms(temporary_terms_derivatives[1],
                          temp_term_union,
                          temporary_terms_idxs,
                          tt_output[1], output_type, tef_terms);

      for (const auto &first_derivative : derivatives[1])
        {
          auto [eq, var] = vectorToTuple<2>(first_derivative.first);
          expr_t d1 = first_derivative.second;

          jacobianHelper(d_output[1], eq, getDynJacobianCol(var), output_type);
          d_output[1] << "=";
          d1->writeOutput(d_output[1], output_type,
                          temp_term_union, temporary_terms_idxs, tef_terms);
          d_output[1] << ";" << endl;
        }
    }

  // Write derivatives for order ≥ 2
  for (size_t i = 2; i < derivatives.size(); i++)
    if (!derivatives[i].empty())
      {
        writeTemporaryTerms(temporary_terms_derivatives[i],
                            temp_term_union,
                            temporary_terms_idxs,
                            tt_output[i], output_type, tef_terms);

        /* When creating the sparse matrix (in MATLAB or C mode), since storage
           is in column-major order, output the first column, then the second,
           then the third. This gives a significant performance boost in use_dll
           mode (at both compilation and runtime), because it facilitates memory
           accesses and expression reusage. */
        ostringstream col0_output, col1_output, col2_output;

        int k = 0; // Current line index in the 3-column matrix
        for (const auto &[vidx, d] : derivatives[i])
          {
            int eq = vidx[0];

            int col_idx = 0;
            for (size_t j = 1; j < vidx.size(); j++)
              {
                col_idx *= dynJacobianColsNbr;
                col_idx += getDynJacobianCol(vidx[j]);
              }

            if (output_type == ExprNodeOutputType::juliaDynamicModel)
              {
                d_output[i] << "    @inbounds " << "g" << i << "[" << eq + 1 << "," << col_idx + 1 << "] = ";
                d->writeOutput(d_output[i], output_type, temp_term_union, temporary_terms_idxs, tef_terms);
                d_output[i] << endl;
              }
            else
              {
                sparseHelper(i, col0_output, k, 0, output_type);
                col0_output << "=" << eq + 1 << ";" << endl;

                sparseHelper(i, col1_output, k, 1, output_type);
                col1_output << "=" << col_idx + 1 << ";" << endl;

                sparseHelper(i, col2_output, k, 2, output_type);
                col2_output << "=";
                d->writeOutput(col2_output, output_type, temp_term_union, temporary_terms_idxs, tef_terms);
                col2_output << ";" << endl;

                k++;
              }

            // Output symetric elements at order 2
            if (i == 2 && vidx[1] != vidx[2])
              {
                int col_idx_sym = getDynJacobianCol(vidx[2]) * dynJacobianColsNbr + getDynJacobianCol(vidx[1]);

                if (output_type == ExprNodeOutputType::juliaDynamicModel)
                  d_output[2] << "    @inbounds g2[" << eq + 1 << "," << col_idx_sym + 1 << "] = "
                              << "g2[" << eq + 1 << "," << col_idx + 1 << "]" << endl;
                else
                  {
                    sparseHelper(2, col0_output, k, 0, output_type);
                    col0_output << "=" << eq + 1 << ";" << endl;

                    sparseHelper(2, col1_output, k, 1, output_type);
                    col1_output << "=" << col_idx_sym + 1 << ";" << endl;

                    sparseHelper(2, col2_output, k, 2, output_type);
                    col2_output << "=";
                    sparseHelper(2, col2_output, k-1, 2, output_type);
                    col2_output << ";" << endl;

                    k++;
                  }
              }
          }
        if (output_type != ExprNodeOutputType::juliaDynamicModel)
          d_output[i] << col0_output.str() << col1_output.str() << col2_output.str();
      }

  if (output_type == ExprNodeOutputType::matlabDynamicModel)
    {
      // Check that we don't have more than 32 nested parenthesis because Matlab does not suppor this. See Issue #1201
      map<string, string> tmp_paren_vars;
      bool message_printed = false;
      for (auto &it : tt_output)
        fixNestedParenthesis(it, tmp_paren_vars, message_printed);
      for (auto &it : d_output)
        fixNestedParenthesis(it, tmp_paren_vars, message_printed);

      ostringstream init_output, end_output;
      init_output << "residual = zeros(" << nrows << ", 1);";
      writeDynamicModelHelper(basename, "dynamic_resid", "residual",
                              "dynamic_resid_tt",
                              temporary_terms_mlv.size() + temporary_terms_derivatives[0].size(),
                              "", init_output, end_output,
                              d_output[0], tt_output[0]);

      init_output.str("");
      init_output << "g1 = zeros(" << nrows << ", " << dynJacobianColsNbr << ");";
      writeDynamicModelHelper(basename, "dynamic_g1", "g1",
                              "dynamic_g1_tt",
                              temporary_terms_mlv.size() + temporary_terms_derivatives[0].size() + temporary_terms_derivatives[1].size(),
                              "dynamic_resid_tt",
                              init_output, end_output,
                              d_output[1], tt_output[1]);
      writeWrapperFunctions(basename, "g1");

      // For order ≥ 2
      int ncols = dynJacobianColsNbr;
      int ntt = temporary_terms_mlv.size() + temporary_terms_derivatives[0].size() + temporary_terms_derivatives[1].size();
      for (size_t i = 2; i < derivatives.size(); i++)
        {
          ncols *= dynJacobianColsNbr;
          ntt += temporary_terms_derivatives[i].size();
          string gname = "g" + to_string(i);
          string vname = "v" + to_string(i);
          string gprevname = "g" + to_string(i-1);

          init_output.str("");
          end_output.str("");
          if (derivatives[i].size())
            {
              init_output << vname << " = zeros(" << NNZDerivatives[i] << ",3);";
              end_output << gname << " = sparse("
                         << vname << "(:,1),"
                         << vname << "(:,2),"
                         << vname << "(:,3),"
                         << nrows << "," << ncols << ");";
            }
          else
            init_output << gname << " = sparse([],[],[]," << nrows << "," << ncols << ");";
          writeDynamicModelHelper(basename, "dynamic_" + gname, gname,
                                  "dynamic_" + gname + "_tt",
                                  ntt,
                                  "dynamic_" + gprevname + "_tt",
                                  init_output, end_output,
                                  d_output[i], tt_output[i]);
          if (i <= 3)
            writeWrapperFunctions(basename, gname);
        }

      writeDynamicMatlabCompatLayer(basename);
    }
  else if (output_type == ExprNodeOutputType::CDynamicModel)
    {
      for (size_t i = 0; i < d_output.size(); i++)
        {
          string funcname = i == 0 ? "resid" : "g" + to_string(i);
          string argname = i == 0 ? "residual" : i == 1 ? "g1" : "v" + to_string(i);
          DynamicOutput << "void dynamic_" << funcname << "_tt(const double *y, const double *x, int nb_row_x, const double *params, const double *steady_state, int it_, double *T)" << endl
                        << "{" << endl
                        << tt_output[i].str()
                        << "}" << endl
                        << endl
                        << "void dynamic_" << funcname << "(const double *y, const double *x, int nb_row_x, const double *params, const double *steady_state, int it_, const double *T, double *" << argname << ")" << endl
                        << "{" << endl;
          if (i == 0)
            DynamicOutput << "  double lhs, rhs;" << endl;
          DynamicOutput << d_output[i].str()
                        << "}" << endl
                        << endl;
        }
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
             << tt_output[0].str()
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
             << d_output[0].str()
             << "    return nothing" << endl
             << "end" << endl << endl;

      // dynamicG1TT!
      output << "function dynamicG1TT!(T::Vector{Float64}," << endl
             << "                      y::Vector{Float64}, x::Matrix{Float64}, "
             << "params::Vector{Float64}, steady_state::Vector{Float64}, it_::Int)" << endl
             << "    dynamicResidTT!(T, y, x, params, steady_state, it_)" << endl
             << tt_output[1].str()
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
             << d_output[1].str()
             << "    return nothing" << endl
             << "end" << endl << endl;

      // dynamicG2TT!
      output << "function dynamicG2TT!(T::Vector{Float64}," << endl
             << "                      y::Vector{Float64}, x::Matrix{Float64}, "
             << "params::Vector{Float64}, steady_state::Vector{Float64}, it_::Int)" << endl
             << "    dynamicG1TT!(T, y, x, params, steady_state, it_)" << endl
             << tt_output[2].str()
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
             << d_output[2].str()
             << "    return nothing" << endl
             << "end" << endl << endl;

      // dynamicG3TT!
      output << "function dynamicG3TT!(T::Vector{Float64}," << endl
             << "                      y::Vector{Float64}, x::Matrix{Float64}, "
             << "params::Vector{Float64}, steady_state::Vector{Float64}, it_::Int)" << endl
             << "    dynamicG2TT!(T, y, x, params, steady_state, it_)" << endl
             << tt_output[3].str()
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
             << d_output[3].str()
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
DynamicModel::writeDynamicJacobianNonZeroElts(const string &basename) const
{
  vector<pair<int, int>> nzij_pred, nzij_current, nzij_fwrd; // pairs (tsid, equation)
  for (const auto &[indices, d1] : derivatives[1])
    {
      if (symbol_table.getType(getSymbIDByDerivID(indices[1])) != SymbolType::endogenous)
        continue;
      int tsid = symbol_table.getTypeSpecificID(getSymbIDByDerivID(indices[1]));
      int lag = getLagByDerivID(indices[1]);
      if (lag == -1)
        nzij_pred.emplace_back(tsid, indices[0]);
      else if (lag == 0)
        nzij_current.emplace_back(tsid, indices[0]);
      else
        nzij_fwrd.emplace_back(tsid, indices[0]);
    }
  sort(nzij_pred.begin(), nzij_pred.end());
  sort(nzij_current.begin(), nzij_current.end());
  sort(nzij_fwrd.begin(), nzij_fwrd.end());

  ofstream output{"+" + basename + "/dynamic_g1_nz.m", ios::out | ios::binary};
  output << "function [nzij_pred, nzij_current, nzij_fwrd] = dynamic_g1_nz()" << endl
         << "% Returns the coordinates of non-zero elements in the Jacobian, in column-major order, for each lead/lag (only for endogenous)" << endl;
  auto print_nzij = [&output](const vector<pair<int, int>> &nzij, const string &name) {
                      output << "  " << name << " = zeros(" << nzij.size() << ", 2, 'int32');" << endl;
                      int idx = 1;
                      for (const auto &it : nzij)
                        {
                          output << "  " << name << "(" << idx << ",1)=" << it.second+1 << ';'
                                 << " " << name << "(" << idx << ",2)=" << it.first+1 << ';' << endl;
                          idx++;
                        }
                    };
  print_nzij(nzij_pred, "nzij_pred");
  print_nzij(nzij_current, "nzij_current");
  print_nzij(nzij_fwrd, "nzij_fwrd");
  output << "end" << endl;
  output.close();
}

void
DynamicModel::parseIncludeExcludeEquations(const string &inc_exc_eq_tags, set<pair<string, string>> &eq_tag_set, bool exclude_eqs)
{
  string tags;
  if (filesystem::exists(inc_exc_eq_tags))
    {
      ifstream exclude_file;
      exclude_file.open(inc_exc_eq_tags, ifstream::in);
      if (!exclude_file.is_open())
        {
          cerr << "ERROR: Could not open " << inc_exc_eq_tags << endl;
          exit(EXIT_FAILURE);
        }

      string line;
      bool tagname_on_first_line = false;
      while (getline(exclude_file, line))
        {
          removeLeadingTrailingWhitespace(line);
          if (!line.empty())
            if (tags.empty() && line.find("=") != string::npos)
              {
                tagname_on_first_line = true;
                tags += line + "(";
              }
            else
              if (line.find("'") != string::npos)
                tags += line + ",";
              else
                tags += "'" + line + "',";
        }

      if (!tags.empty())
        {
          tags = tags.substr(0, tags.size()-1);
          if (tagname_on_first_line)
            tags += ")";
        }
    }
  else
    tags = inc_exc_eq_tags;
  removeLeadingTrailingWhitespace(tags);

  if (tags.front() == '[' && tags.back() != ']')
    {
      cerr << "Error: " << (exclude_eqs ? "exclude_eqs" : "include_eqs")
           << ": if the first character is '[' the last must be ']'" << endl;
      exit(EXIT_FAILURE);
    }

  if (tags.front() == '[' && tags.back() == ']')
    tags = tags.substr(1, tags.length() - 2);
  removeLeadingTrailingWhitespace(tags);

  regex q(R"(^\w+\s*=)");
  smatch matches;
  string tagname = "name";
  if (regex_search(tags, matches, q))
    {
      tagname = matches[0].str();
      tags = tags.substr(tagname.size(), tags.length() - tagname.size() + 1);
      removeLeadingTrailingWhitespace(tags);
      if (tags.front() == '(' && tags.back() == ')')
        {
          tags = tags.substr(1, tags.length() - 2);
          removeLeadingTrailingWhitespace(tags);
        }
      tagname = tagname.substr(0, tagname.size()-1);
      removeLeadingTrailingWhitespace(tagname);
    }

  string quote_regex = "'[^']+'";
  string non_quote_regex = R"([^,\s]+)";
  regex r(R"((\s*)" + quote_regex + "|" + non_quote_regex + R"(\s*)(,\s*()" + quote_regex + "|" + non_quote_regex + R"()\s*)*)");
  if (!regex_match(tags, r))
    {
      cerr << "Error: " << (exclude_eqs ? "exclude_eqs" : "include_eqs")
           << ": argument is of incorrect format." << endl;
      exit(EXIT_FAILURE);
    }

  regex s(quote_regex + "|" + non_quote_regex);
  for (auto it = sregex_iterator(tags.begin(), tags.end(), s);
       it != sregex_iterator(); ++it)
    {
      auto str = it->str();
      if (str[0] == '\'' && str[str.size()-1] == '\'')
        str = str.substr(1, str.size()-2);
      eq_tag_set.insert({tagname, str});
    }
}

void
DynamicModel::includeExcludeEquations(const string &eqs, bool exclude_eqs)
{
  if (eqs.empty())
    return;

  set<pair<string, string>> eq_tag_set;
  parseIncludeExcludeEquations(eqs, eq_tag_set, exclude_eqs);

  vector<int> excluded_vars
    = ModelTree::includeExcludeEquations(eq_tag_set, exclude_eqs,
                                         equations, equations_lineno,
                                         equation_tags, equation_tags_xref, false);

  // `static_only_equation_tags` is `vector<vector<pair<string, string>>>`
  // while `equation_tags` is `vector<pair<int, pair<string, string>>>`
  // so convert former structure to latter to conform with function call
  int n = 0;
  vector<pair<int, pair<string, string>>> tmp_static_only_equation_tags;
  for (auto &eqn_tags : static_only_equations_equation_tags)
    {
      for (auto &eqn_tag : eqn_tags)
        tmp_static_only_equation_tags.emplace_back(make_pair(n, eqn_tag));
      n++;
    }
  // Ignore output because variables are not excluded when equations marked 'static' are excluded
  ModelTree::includeExcludeEquations(eq_tag_set, exclude_eqs,
                                     static_only_equations, static_only_equations_lineno,
                                     tmp_static_only_equation_tags,
                                     static_only_equation_tags_xref, true);
  if (!eq_tag_set.empty())
    {
      cerr << "ERROR: " << (exclude_eqs ? "exclude_eqs" : "include_eqs") << ": The equations specified by `";
      cerr << eq_tag_set.begin()->first << "= ";
      for (auto &it : eq_tag_set)
        cerr << it.second << ", ";
      cerr << "` were not found." << endl;
      exit(EXIT_FAILURE);
    }

  if (staticOnlyEquationsNbr() != dynamicOnlyEquationsNbr())
    {
      cerr << "ERROR: " << (exclude_eqs ? "exclude_eqs" : "include_eqs")
           << ": You must remove the same number of equations marked `static` as equations marked `dynamic`." << endl;
      exit(EXIT_FAILURE);
    }

  // convert back static equation info
  if (static_only_equations.empty())
    static_only_equations_equation_tags.clear();
  else
    {
      static_only_equations_equation_tags.resize(static_only_equations.size());
      fill(static_only_equations_equation_tags.begin(), static_only_equations_equation_tags.end(), vector<pair<string, string>>());
      for (auto &it : tmp_static_only_equation_tags)
        static_only_equations_equation_tags.at(it.first).emplace_back(it.second);
    }

  // Collect list of used variables in updated list of equations
  set<pair<int, int>> eqn_vars;
  for (const auto &eqn : equations)
    eqn->collectDynamicVariables(SymbolType::endogenous, eqn_vars);
  for (const auto &eqn : static_only_equations)
    eqn->collectDynamicVariables(SymbolType::endogenous, eqn_vars);

  // Change LHS variable type of excluded equation if it is used in an eqution that has been kept
  for (auto ev : excluded_vars)
    {
      bool found = false;
      for (const auto &it : eqn_vars)
        if (it.first == ev)
          {
            symbol_table.changeType(ev, SymbolType::exogenous);
            found = true;
            break;
          }
      if (!found)
        symbol_table.changeType(ev, SymbolType::excludedVariable);
    }
}

void
DynamicModel::writeOutput(ostream &output, const string &basename, bool block_decomposition, bool linear_decomposition, bool byte_code, bool use_dll, bool estimation_present, bool compute_xrefs, bool julia) const
{
  /* Writing initialisation for M_.lead_lag_incidence matrix
     M_.lead_lag_incidence is a matrix with as many columns as there are
     endogenous variables and as many rows as there are periods in the
     models (nbr of rows = M_.max_lag+M_.max_lead+1)

     The matrix elements are equal to zero if a variable isn't present in the
     model at a given period.
  */

  string modstruct, outstruct;
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
         << modstruct << "orig_maximum_lag_with_diffs_expanded = " << max_lag_with_diffs_expanded_orig << ";" << endl
         << modstruct << "lead_lag_incidence = [";
  // Loop on endogenous variables
  int nstatic = 0,
    nfwrd = 0,
    npred = 0,
    nboth = 0;
  for (int endoID = 0; endoID < symbol_table.endo_nbr(); endoID++)
    {
      output << endl;
      int sstatic = 1,
        sfwrd = 0,
        spred = 0,
        sboth = 0;
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
      nfwrd += sfwrd;
      npred += spred;
      nboth += sboth;
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
    {
      output << modstruct << "dynamic_tmp_nbr = [";
      for (size_t i = 0; i < temporary_terms_derivatives.size(); i++)
        output << temporary_terms_derivatives[i].size() + (i == 0 ? temporary_terms_mlv.size() : 0) << "; ";
      output << "];" << endl;
    }

  // Write equation tags
  if (julia)
    {
      output << modstruct << "equation_tags = [" << endl;
      for (const auto &equation_tag : equation_tags)
        output << "                       EquationTag("
               << equation_tag.first + 1 << R"( , ")"
               << equation_tag.second.first << R"(" , ")"
               << equation_tag.second.second << R"("))" << endl;
      output << "                      ]" << endl;
    }
  else
    {
      output << modstruct << "equations_tags = {" << endl;
      for (const auto &equation_tag : equation_tags)
        output << "  " << equation_tag.first + 1 << " , '"
               << equation_tag.second.first << "' , '"
               << equation_tag.second.second << "' ;" << endl;
      output << "};" << endl;
    }

  // Write mapping for variables and equations they are present in
  for (const auto &variable : variableMapping)
    {
      output << modstruct << "mapping." << symbol_table.getName(variable.first) << ".eqidx = [";
      for (auto equation : variable.second)
        output << equation + 1 << " ";
      output << "];" << endl;
    }

  /* Say if static and dynamic models differ (because of [static] and [dynamic]
     equation tags) */
  output << modstruct << "static_and_dynamic_models_differ = "
         << (static_only_equations.size() > 0 ? "true" : "false")
         << ";" << endl;

  // Say if model contains an external function call
  bool has_external_function = false;
  for (auto equation : equations)
    if (equation->containsExternalFunction())
      {
        has_external_function = true;
        break;
      }
  output << modstruct << "has_external_function = "
         << (has_external_function ? "true" : "false")
         << ';' << endl;

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
          max_lag = max_leadlag_block[block].first;
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
          for (const auto &it : exo_block[block])
            for (int it1 : it.second)
              exogenous.insert(it1);
          set<int> exogenous_det;
          for (const auto &it : exo_det_block[block])
            for (int it1 : it.second)
              exogenous_det.insert(it1);
          set<int> other_endogenous;
          for (const auto &it : other_endo_block[block])
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
          for (int exo : exogenous)
            if (exo >= 0)
              {
                output << " " << exo+1;
                i++;
              }
          output << "];" << endl
                 << "block_structure.block(" << block+1 << ").exogenous_det = [";
          i = 0;
          for (int exo_det : exogenous_det)
            if (exo_det >= 0)
              {
                output << " " << exo_det+1;
                i++;
              }
          output << "];" << endl
                 << "block_structure.block(" << block+1 << ").exo_det_nbr = " << i << ";" << endl
                 << "block_structure.block(" << block+1 << ").other_endogenous = [";
          i = 0;
          for (int other_endo : other_endogenous)
            if (other_endo >= 0)
              {
                output << " " << other_endo+1;
                i++;
              }
          output << "];" << endl
                 << "block_structure.block(" << block+1 << ").other_endogenous_block = [";
          i = 0;
          for (int other_endo : other_endogenous)
            if (other_endo >= 0)
              {
                bool OK = true;
                unsigned int j;
                for (j = 0; j < block && OK; j++)
                  for (unsigned int k = 0; k < getBlockSize(j) && OK; k++)
                    OK = other_endo != getBlockVariableID(j, k);
                if (!OK)
                  output << " " << j;
                i++;
              }
          output << "];" << endl;

          output << "block_structure.block(" << block+1 << ").tm1 = zeros(" << i << ", " << state_var.size() << ");" << endl;
          int count_other_endogenous = 1;
          for (int other_endo : other_endogenous)
            {
              for (auto it = state_var.begin(); it != state_var.end(); ++it)
                if (*it == other_endo + 1)
                  output << "block_structure.block(" << block+1 << ").tm1("
                         << count_other_endogenous << ", "
                         << it - state_var.begin()+1 << ") = 1;" << endl;
              count_other_endogenous++;
            }

          output << "block_structure.block(" << block+1 << ").other_endo_nbr = " << i << ";" << endl;

          tmp_s.str("");
          count_lead_lag_incidence = 0;
          dynamic_jacob_map_t reordered_dynamic_jacobian;
          for (const auto &it : blocks_derivatives[block])
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
          for (int &it_l : local_state_var)
            for (auto it = state_var.begin(); it != state_var.end(); ++it)
              if (*it == it_l)
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
              for (int other_endo : other_endogenous)
                {
                  bool done = false;
                  for (int i = 0; i < block_size; i++)
                    {
                      unsigned int eq = getBlockEquationID(block, i);
                      if (derivative_other_endo[block].find({ lag, eq, other_endo })
                          != derivative_other_endo[block].end())
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

      map<tuple<int, int, int>, int> lag_row_incidence;
      for (const auto & [indices, d1] : derivatives[1])
        {
          int deriv_id = indices[1];
          if (getTypeByDerivID(deriv_id) == SymbolType::endogenous)
            {
              int eq = indices[0];
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
          filesystem::create_directories(basename + "/model/bytecode");
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
                  if (find(state_var.begin(), state_var.end(), var+1) != state_var.end())
                    nze++;
                }
              if (block == 0)
                {
                  set<pair<int, int>> row_state_var_incidence;
                  for (const auto &it : blocks_derivatives[block])
                    if (auto it_state_var = find(state_var.begin(), state_var.end(), getBlockVariableID(block, get<1>(it))+1);
                        it_state_var != state_var.end())
                      if (auto it_state_equ = find(state_equ.begin(), state_equ.end(), getBlockEquationID(block, get<0>(it))+1);
                          it_state_equ != state_equ.end())
                        row_state_var_incidence.emplace(it_state_equ - state_equ.begin(), it_state_var - state_var.begin());
                  auto row_state_var_incidence_it = row_state_var_incidence.begin();
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
                  set<pair<int, int>> col_state_var_incidence;
                  for (const auto &it : row_state_var_incidence)
                    col_state_var_incidence.emplace(it.second, it.first);
                  auto col_state_var_incidence_it = col_state_var_incidence.begin();
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
            for (int j = n_obs; j < n; j++)
              {
                int j1 = j - n_obs;
                int j1_n_state = j1 * n_state - n_obs;
                if ((i < n_obs) || (i >= nb_diag + n_obs) || (j1 >= nb_diag))
                  for (int k = n_obs; k < i_nz_state_var[i]; k++)
                    v_index_KF.emplace_back(i + j1 * n, pair(i + k * n, k + j1_n_state));
              }
          int size_v_index_KF = v_index_KF.size();

          KF_index_file.write(reinterpret_cast<char *>(&size_v_index_KF), sizeof(size_v_index_KF));
          for (auto &it : v_index_KF)
            KF_index_file.write(reinterpret_cast<char *>(&it), sizeof(index_KF));

          vector<index_KF> v_index_KF_2;
          int n_n_obs = n * n_obs;
          for (int i = 0; i < n; i++)
            for (int j = i; j < n; j++)
              {
                if ((i < n_obs) || (i >= nb_diag + n_obs) || (j < n_obs) || (j >= nb_diag + n_obs))
                  for (int k = n_obs; k < i_nz_state_var[j]; k++)
                    {
                      int k_n = k * n;
                      v_index_KF_2.emplace_back(i * n + j, pair(i + k_n - n_n_obs, j + k_n));
                    }
              }
          int size_v_index_KF_2 = v_index_KF_2.size();

          KF_index_file.write(reinterpret_cast<char *>(&size_v_index_KF_2), sizeof(size_v_index_KF_2));
          for (auto &it : v_index_KF_2)
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

  // FIXME: implement this for Julia
  if (!julia)
    {
      string empty_cell = "cell(" + to_string(symbol_table.endo_nbr()) + ", 1)";
      output << modstruct << "endo_trends = struct('deflator', " << empty_cell
             << ", 'log_deflator', " << empty_cell << ", 'growth_factor', " << empty_cell
             << ", 'log_growth_factor', " << empty_cell << ");" << endl;
      for (int i = 0; i < symbol_table.endo_nbr(); i++)
        {
          int symb_id = symbol_table.getID(SymbolType::endogenous, i);
          if (auto it = nonstationary_symbols_map.find(symb_id); it != nonstationary_symbols_map.end())
            {
              auto [is_log, deflator] = it->second;
              output << modstruct << "endo_trends(" << i << ")."
                     << (is_log ? "log_deflator" : "deflator") << " = '";
              deflator->writeJsonOutput(output, {}, {});
              output << "';" << endl;

              auto growth_factor = const_cast<DynamicModel *>(this)->AddDivide(deflator, deflator->decreaseLeadsLags(1))->removeTrendLeadLag(trend_symbols_map)->replaceTrendVar();
              output << modstruct << "endo_trends(" << i << ")."
                     << (is_log ? "log_growth_factor" : "growth_factor") << " = '";
              growth_factor->writeJsonOutput(output, {}, {});
              output << "';" << endl;
            }
        }
    }

  if (compute_xrefs)
    writeXrefs(output);

  // Write number of non-zero derivatives
  // Use -1 if the derivatives have not been computed
  output << modstruct << (julia ? "nnzderivatives" : "NNZDerivatives") << " = [";
  for (int i = 1; i < static_cast<int>(NNZDerivatives.size()); i++)
    output << (i > computed_derivs_order ? -1 : NNZDerivatives[i]) << "; ";
  output << "];" << endl;

  // Write Pac Model Consistent Expectation parameter info
  for (auto &it : pac_mce_alpha_symb_ids)
    {
      output << modstruct << "pac." << it.first.first << ".equations." << it.first.second << ".mce.alpha = [";
      for (auto it : it.second)
        output << symbol_table.getTypeSpecificID(it) + 1 << " ";
      output << "];" << endl;
    }

  // Write Pac Model Consistent Expectation Z1 info
  for (auto &it : pac_mce_z1_symb_ids)
    output << modstruct << "pac." << it.first.first << ".equations." << it.first.second << ".mce.z1 = "
           << symbol_table.getTypeSpecificID(it.second) + 1 << ";" << endl;

  // Write Pac lag info
  for (auto &it : pac_eqtag_and_lag)
    output << modstruct << "pac." << it.first.first << ".equations." << it.second.first << ".max_lag = " << it.second.second << ";" << endl;

  // Write Pac equation tag info
  map<string, vector<pair<string, string>>> for_writing;
  for (auto &it : pac_eqtag_and_lag)
    for_writing[it.first.first].emplace_back(it.first.second, it.second.first);

  for (auto &it : for_writing)
    {
      output << modstruct << "pac." << it.first << ".tag_map = [";
      for (auto &it1 : it.second)
        output << "{'" << it1.first << "', '" << it1.second << "'};";
      output << "];" << endl;
    }

  for (auto &it : pac_model_info)
    {
      vector<int> lhs = get<0>(it.second);
      output << modstruct << "pac." << it.first << ".lhs = [";
      for (auto it : lhs)
        output << it + 1 << " ";
      output << "];" << endl;

      if (int growth_param_index = get<1>(it.second);
          growth_param_index >= 0)
        output << modstruct << "pac." << it.first << ".growth_neutrality_param_index = "
               << symbol_table.getTypeSpecificID(growth_param_index) + 1 << ";" << endl;

      output << modstruct << "pac." << it.first << ".auxiliary_model_type = '" << get<2>(it.second) << "';" << endl;
    }

  for (auto &pit : pac_equation_info)
    {
      auto [lhs_pac_var, optim_share_index, ar_params_and_vars, ec_params_and_vars, non_optim_vars_params_and_constants, additive_vars_params_and_constants, optim_additive_vars_params_and_constants] = pit.second;
      string substruct = pit.first.first + ".equations." + pit.first.second + ".";

      output << modstruct << "pac." << substruct << "lhs_var = "
             << symbol_table.getTypeSpecificID(lhs_pac_var.first) + 1 << ";" << endl;

      if (optim_share_index >= 0)
        output << modstruct << "pac." << substruct << "share_of_optimizing_agents_index = "
               << symbol_table.getTypeSpecificID(optim_share_index) + 1 << ";" << endl;

      output << modstruct << "pac." << substruct << "ec.params = "
             << symbol_table.getTypeSpecificID(ec_params_and_vars.first) + 1 << ";" << endl
             << modstruct << "pac." << substruct << "ec.vars = [";
      for (auto it : ec_params_and_vars.second)
        output << symbol_table.getTypeSpecificID(get<0>(it)) + 1 << " ";
      output << "];" << endl
             << modstruct << "pac." << substruct << "ec.istarget = [";
      for (auto it : ec_params_and_vars.second)
        output << (get<1>(it) ? "true " : "false ");
      output << "];" << endl
             << modstruct << "pac." << substruct << "ec.scale = [";
      for (auto it : ec_params_and_vars.second)
        output << get<2>(it) << " ";
      output << "];" << endl
             << modstruct << "pac." << substruct << "ec.isendo = [";
      for (auto it : ec_params_and_vars.second)
        switch (symbol_table.getType(get<0>(it)))
          {
          case SymbolType::endogenous:
            output << "true ";
            break;
          case SymbolType::exogenous:
            output << "false ";
            break;
          default:
            cerr << "expecting endogenous or exogenous" << endl;
            exit(EXIT_FAILURE);
          }
      output << "];" << endl
             << modstruct << "pac." << substruct << "ar.params = [";
      for (auto &it : ar_params_and_vars)
        output << symbol_table.getTypeSpecificID(it.first) + 1 << " ";
      output << "];" << endl
             << modstruct << "pac." << substruct << "ar.vars = [";
      for (auto &it : ar_params_and_vars)
        output << symbol_table.getTypeSpecificID(it.second.first) + 1 << " ";
      output << "];" << endl
             << modstruct << "pac." << substruct << "ar.lags = [";
      for (auto &it : ar_params_and_vars)
        output << it.second.second << " ";
      output << "];" << endl;
      if (!non_optim_vars_params_and_constants.empty())
        {
          output << modstruct << "pac." << substruct << "non_optimizing_behaviour.params = [";
          for (auto &it : non_optim_vars_params_and_constants)
            if (get<2>(it) >= 0)
              output << symbol_table.getTypeSpecificID(get<2>(it)) + 1 << " ";
            else
              output << "NaN ";
          output << "];" << endl
                 << modstruct << "pac." << substruct << "non_optimizing_behaviour.vars = [";
          for (auto &it : non_optim_vars_params_and_constants)
            output << symbol_table.getTypeSpecificID(get<0>(it)) + 1 << " ";
          output << "];" << endl
                 << modstruct << "pac." << substruct << "non_optimizing_behaviour.isendo = [";
          for (auto &it : non_optim_vars_params_and_constants)
            switch (symbol_table.getType(get<0>(it)))
              {
              case SymbolType::endogenous:
                output << "true ";
                break;
              case SymbolType::exogenous:
                output << "false ";
                break;
              default:
                cerr << "expecting endogenous or exogenous" << endl;
                exit(EXIT_FAILURE);
              }
          output << "];" << endl
                 << modstruct << "pac." << substruct << "non_optimizing_behaviour.lags = [";
          for (auto &it : non_optim_vars_params_and_constants)
            output << get<1>(it) << " ";
          output << "];" << endl
                 << modstruct << "pac." << substruct << "non_optimizing_behaviour.scaling_factor = [";
          for (auto &it : non_optim_vars_params_and_constants)
            output << get<3>(it) << " ";
          output << "];" << endl;
        }
      if (!additive_vars_params_and_constants.empty())
        {
          output << modstruct << "pac." << substruct << "additive.params = [";
          for (auto &it : additive_vars_params_and_constants)
            if (get<2>(it) >= 0)
              output << symbol_table.getTypeSpecificID(get<2>(it)) + 1 << " ";
            else
              output << "NaN ";
          output << "];" << endl
                 << modstruct << "pac." << substruct << "additive.vars = [";
          for (auto &it : additive_vars_params_and_constants)
            output << symbol_table.getTypeSpecificID(get<0>(it)) + 1 << " ";
          output << "];" << endl
                 << modstruct << "pac." << substruct << "additive.isendo = [";
          for (auto &it : additive_vars_params_and_constants)
            switch (symbol_table.getType(get<0>(it)))
              {
              case SymbolType::endogenous:
                output << "true ";
                break;
              case SymbolType::exogenous:
                output << "false ";
                break;
              default:
                cerr << "expecting endogenous or exogenous" << endl;
                exit(EXIT_FAILURE);
              }
          output << "];" << endl
                 << modstruct << "pac." << substruct << "additive.lags = [";
          for (auto &it : additive_vars_params_and_constants)
            output << get<1>(it) << " ";
          output << "];" << endl
                 << modstruct << "pac." << substruct << "additive.scaling_factor = [";
          for (auto &it : additive_vars_params_and_constants)
            output << get<3>(it) << " ";
          output << "];" << endl;
        }
      if (!optim_additive_vars_params_and_constants.empty())
        {
          output << modstruct << "pac." << substruct << "optim_additive.params = [";
          for (auto &it : optim_additive_vars_params_and_constants)
            if (get<2>(it) >= 0)
              output << symbol_table.getTypeSpecificID(get<2>(it)) + 1 << " ";
            else
              output << "NaN ";
          output << "];" << endl
                 << modstruct << "pac." << substruct << "optim_additive.vars = [";
          for (auto &it : optim_additive_vars_params_and_constants)
            output << symbol_table.getTypeSpecificID(get<0>(it)) + 1 << " ";
          output << "];" << endl
                 << modstruct << "pac." << substruct << "optim_additive.isendo = [";
          for (auto &it : optim_additive_vars_params_and_constants)
            switch (symbol_table.getType(get<0>(it)))
              {
              case SymbolType::endogenous:
                output << "true ";
                break;
              case SymbolType::exogenous:
                output << "false ";
                break;
              default:
                cerr << "expecting endogenous or exogenous" << endl;
                exit(EXIT_FAILURE);
              }
          output << "];" << endl
                 << modstruct << "pac." << substruct << "optim_additive.lags = [";
          for (auto &it : optim_additive_vars_params_and_constants)
            output << get<1>(it) << " ";
          output << "];" << endl
                 << modstruct << "pac." << substruct << "optim_additive.scaling_factor = [";
          for (auto &it : optim_additive_vars_params_and_constants)
            output << get<3>(it) << " ";
          output << "];" << endl;
        }
      // Create empty h0 and h1 substructures that will be overwritten later if not empty
      output << modstruct << "pac." << substruct << "h0_param_indices = [];" << endl
             << modstruct << "pac." << substruct << "h1_param_indices = [];" << endl;
    }

  for (auto &it : pac_h0_indices)
    {
      output << modstruct << "pac." << it.first.first << ".equations." << it.first.second << ".h0_param_indices = [";
      for (auto it1 : it.second)
        output << symbol_table.getTypeSpecificID(it1) + 1 << " ";
      output << "];" << endl;
    }

  for (auto &it : pac_h1_indices)
    {
      output << modstruct << "pac." << it.first.first << ".equations." << it.first.second << ".h1_param_indices = [";
      for (auto it1 : it.second)
        output << symbol_table.getTypeSpecificID(it1) + 1 << " ";
      output << "];" << endl;
    }
}

map<tuple<int, int, int>, expr_t>
DynamicModel::collect_first_order_derivatives_endogenous()
{
  map<tuple<int, int, int>, expr_t> endo_derivatives;
  for (auto & [indices, d1] : derivatives[1])
    if (getTypeByDerivID(indices[1]) == SymbolType::endogenous)
      {
        int eq = indices[0];
        int var = symbol_table.getTypeSpecificID(getSymbIDByDerivID(indices[1]));
        int lag = getLagByDerivID(indices[1]);
        endo_derivatives[{ eq, var, lag }] = d1;
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
      for (const auto &it : eqnums)
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
              equations[eqn]->arg2->collectDynamicVariables(SymbolType::endogenous, rhs_set);
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
                  int trend_var_symb_id = equations[eqn]->arg2->findTargetVariable(lhs_symb_id);
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

  for (const auto &it : eqtags)
    {
      vector<int> eqnumber, lhs;
      vector<expr_t> lhs_expr_t;
      vector<set<pair<int, int>>> rhs;

      for (const auto &eqtag : it.second)
        {
          int eqn = -1;
          set<pair<int, int>> lhs_set, lhs_tmp_set, rhs_set;
          for (const auto &equation_tag : equation_tags)
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

          equations[eqn]->arg1->collectDynamicVariables(SymbolType::endogenous, lhs_set);
          equations[eqn]->arg1->collectDynamicVariables(SymbolType::exogenous, lhs_tmp_set);
          equations[eqn]->arg1->collectDynamicVariables(SymbolType::parameter, lhs_tmp_set);

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
          equations[eqn]->arg1->collectVARLHSVariable(lhs_expr_t_set);
          lhs_expr_t.push_back(*(lhs_expr_t_set.begin()));

          equations[eqn]->arg2->collectDynamicVariables(SymbolType::endogenous, rhs_set);
          for (const auto &itrhs : rhs_set)
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
  var_model_table.setAR(fillAutoregressiveMatrix(true));
}

void
DynamicModel::fillVarModelTableFromOrigModel() const
{
  map<string, vector<int>> lags, orig_diff_var;
  map<string, vector<bool>> diff;
  for (const auto &it : var_model_table.getEqNums())
    {
      set<expr_t> lhs;
      vector<int> orig_diff_var_vec;
      vector<bool> diff_vec;
      for (auto eqn : it.second)
        {
          // ensure no leads in equations
          if (equations[eqn]->arg2->VarMinLag() <= 0)
            {
              cerr << "ERROR in VAR model Equation (#" << eqn << "). "
                   << "Leaded exogenous variables "
                   << "and leaded or contemporaneous endogenous variables not allowed in VAR"
                   << endl;
              exit(EXIT_FAILURE);
            }

          // save lhs variables
          equations[eqn]->arg1->collectVARLHSVariable(lhs);

          equations[eqn]->arg1->countDiffs() > 0 ?
            diff_vec.push_back(true) : diff_vec.push_back(false);
          if (diff_vec.back())
            {
              set<pair<int, int>> diff_set;
              equations[eqn]->arg1->collectDynamicVariables(SymbolType::endogenous, diff_set);

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

      set<expr_t> lhs_lag_equiv;
      for (const auto &lh : lhs)
        {
          auto [lag_equiv_repr, index] = lh->getLagEquivalenceClass();
          lhs_lag_equiv.insert(lag_equiv_repr);
        }

      vector<int> max_lag;
      for (auto eqn : it.second)
        max_lag.push_back(equations[eqn]->arg2->VarMaxLag(lhs_lag_equiv));
      lags[it.first] = max_lag;
      diff[it.first] = diff_vec;
      orig_diff_var[it.first] = orig_diff_var_vec;
    }
  var_model_table.setDiff(diff);
  var_model_table.setMaxLags(lags);
  var_model_table.setOrigDiffVar(orig_diff_var);
}

map<string, map<tuple<int, int, int>, expr_t>>
DynamicModel::fillAutoregressiveMatrix(bool is_var) const
{
  map<string, map<tuple<int, int, int>, expr_t>> ARr;
  auto eqnums = is_var ?
    var_model_table.getEqNums() : trend_component_model_table.getNonTargetEqNums();
  for (const auto &it : eqnums)
    {
      int i = 0;
      map<tuple<int, int, int>, expr_t> AR;
      vector<int> lhs = is_var ?
        var_model_table.getLhsOrigIds(it.first) : trend_component_model_table.getNonTargetLhs(it.first);
      for (auto eqn : it.second)
        {
          auto bopn = dynamic_cast<BinaryOpNode *>(equations[eqn]->arg2);
          bopn->fillAutoregressiveRow(i++, lhs, AR);
        }
      ARr[it.first] = AR;
    }
  return ARr;
}

void
DynamicModel::fillTrendComponentModelTable() const
{
  map<string, vector<int>> eqnums, trend_eqnums, lhsr;
  map<string, vector<expr_t>> lhs_expr_tr;
  map<string, vector<set<pair<int, int>>>> rhsr;
  map<string, vector<string>> eqtags = trend_component_model_table.getEqTags();
  map<string, vector<string>> trend_eqtags = trend_component_model_table.getTargetEqTags();
  for (const auto &it : trend_eqtags)
    {
      vector<int> trend_eqnumber;
      for (const auto &eqtag : it.second)
        {
          int eqn = -1;
          for (const auto &equation_tag : equation_tags)
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

  for (const auto &it : eqtags)
    {
      vector<int> eqnumber, lhs;
      vector<expr_t> lhs_expr_t;
      vector<set<pair<int, int>>> rhs;

      for (const auto &eqtag : it.second)
        {
          int eqn = -1;
          set<pair<int, int>> lhs_set, lhs_tmp_set, rhs_set;
          for (const auto &equation_tag : equation_tags)
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

          equations[eqn]->arg1->collectDynamicVariables(SymbolType::endogenous, lhs_set);
          equations[eqn]->arg1->collectDynamicVariables(SymbolType::exogenous, lhs_tmp_set);
          equations[eqn]->arg1->collectDynamicVariables(SymbolType::parameter, lhs_tmp_set);

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
          equations[eqn]->arg1->collectVARLHSVariable(lhs_expr_t_set);
          lhs_expr_t.push_back(*(lhs_expr_t_set.begin()));

          equations[eqn]->arg2->collectDynamicVariables(SymbolType::endogenous, rhs_set);
          for (const auto &itrhs : rhs_set)
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

pair<map<string, map<tuple<int, int, int>, expr_t>>, map<string, map<tuple<int, int, int>, expr_t>>>
DynamicModel::fillErrorComponentMatrix(const ExprNode::subst_table_t &diff_subst_table) const
{
  map<string, map<tuple<int, int, int>, expr_t>> A0r, A0starr;

  for (const auto &it : trend_component_model_table.getEqNums())
    {
      int i = 0;
      map<tuple<int, int, int>, expr_t> A0, A0star;
      vector<int> target_lhs = trend_component_model_table.getTargetLhs(it.first);
      vector<int> nontarget_eqnums = trend_component_model_table.getNonTargetEqNums(it.first);
      vector<int> undiff_nontarget_lhs = getUndiffLHSForPac(it.first, diff_subst_table);
      vector<int> parsed_undiff_nontarget_lhs;

      for (auto eqn : it.second)
        {
          if (find(nontarget_eqnums.begin(), nontarget_eqnums.end(), eqn) != nontarget_eqnums.end())
            parsed_undiff_nontarget_lhs.push_back(undiff_nontarget_lhs.at(i));
          i++;
        }

      i = 0;
      for (auto eqn : it.second)
        if (find(nontarget_eqnums.begin(), nontarget_eqnums.end(), eqn) != nontarget_eqnums.end())
          equations[eqn]->arg2->fillErrorCorrectionRow(i++, parsed_undiff_nontarget_lhs, target_lhs, A0, A0star);
      A0r[it.first] = A0;
      A0starr[it.first] = A0star;
    }

  return { A0r, A0starr };
}

void
DynamicModel::fillTrendComponentModelTableFromOrigModel() const
{
  map<string, vector<int>> lags, orig_diff_var;
  map<string, vector<bool>> diff;
  for (const auto &it : trend_component_model_table.getEqNums())
    {
      set<expr_t> lhs;
      vector<int> orig_diff_var_vec;
      vector<bool> diff_vec;
      for (auto eqn : it.second)
        {
          // ensure no leads in equations
          if (equations[eqn]->arg2->VarMinLag() <= 0)
            {
              cerr << "ERROR in trend component model Equation (#" << eqn << "). "
                   << "Leaded exogenous variables "
                   << "and leaded or contemporaneous endogenous variables not allowed in VAR"
                   << endl;
              exit(EXIT_FAILURE);
            }

          // save lhs variables
          equations[eqn]->arg1->collectVARLHSVariable(lhs);

          if (equations[eqn]->arg1->countDiffs() > 0)
            diff_vec.push_back(true);
          else
            diff_vec.push_back(false);
          if (diff_vec.back())
            {
              set<pair<int, int>> diff_set;
              equations[eqn]->arg1->collectDynamicVariables(SymbolType::endogenous, diff_set);

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

      set<expr_t> lhs_lag_equiv;
      for (const auto &lh : lhs)
        {
          auto [lag_equiv_repr, index] = lh->getLagEquivalenceClass();
          lhs_lag_equiv.insert(lag_equiv_repr);
        }

      vector<int> max_lag;
      for (auto eqn : it.second)
        max_lag.push_back(equations[eqn]->arg2->VarMaxLag(lhs_lag_equiv));
      lags[it.first] = max_lag;
      diff[it.first] = diff_vec;
      orig_diff_var[it.first] = orig_diff_var_vec;
    }
  trend_component_model_table.setDiff(diff);
  trend_component_model_table.setMaxLags(lags);
  trend_component_model_table.setOrigDiffVar(orig_diff_var);
}

void
DynamicModel::fillTrendComponentmodelTableAREC(const ExprNode::subst_table_t &diff_subst_table) const
{
  auto ARr = fillAutoregressiveMatrix(false);
  trend_component_model_table.setAR(ARr);
  auto [A0r, A0starr] = fillErrorComponentMatrix(diff_subst_table);
  trend_component_model_table.setA0(A0r, A0starr);
}

void
DynamicModel::addEquationsForVar()
{
  if (var_model_table.empty())
    return;
  auto var_symbol_list_and_order = var_model_table.getSymbolListAndOrder();

  // List of endogenous variables and the minimum lag value that must exist in the model equations
  map<string, int> var_endos_and_lags, model_endos_and_lags;
  for (const auto &it : var_symbol_list_and_order)
    for (auto &equation : equations)
      if (equation->isVarModelReferenced(it.first))
        {
          vector<string> symbol_list = it.second.first.get_symbols();
          int order = it.second.second;
          for (auto &it1 : symbol_list)
            if (order > 2)
              if (var_endos_and_lags.find(it1) != var_endos_and_lags.end())
                var_endos_and_lags[it1] = min(var_endos_and_lags[it1], -order);
              else
                var_endos_and_lags[it1] = -order;
          break;
        }

  if (var_endos_and_lags.empty())
    return;

  // Ensure that the minimum lag value exists in the model equations.
  // If not, add an equation for it
  for (auto &equation : equations)
    equation->getEndosAndMaxLags(model_endos_and_lags);

  int count = 0;
  for (auto &it : var_endos_and_lags)
    if (auto it2 = model_endos_and_lags.find(it.first);
        it2 == model_endos_and_lags.end())
      cerr << "WARNING: Variable used in VAR that is not used in the model: " << it.first << endl;
    else
      if (it.second < it2->second)
        {
          int symb_id = symbol_table.getID(it.first);
          expr_t newvar = AddVariable(symb_id, it.second);
          expr_t auxvar = AddVariable(symbol_table.addVarModelEndoLagAuxiliaryVar(symb_id, it.second, newvar), 0);
          addEquation(AddEqual(newvar, auxvar), -1);
          addAuxEquation(AddEqual(newvar, auxvar));
          count++;
        }

  if (count > 0)
    cout << "Accounting for var_model lags not in model block: added "
         << count << " auxiliary variables and equations." << endl;
}

vector<int>
DynamicModel::getUndiffLHSForPac(const string &aux_model_name,
                                 const ExprNode::subst_table_t &diff_subst_table) const
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
      for (auto it1 = eqnumber.begin(); it1 != eqnumber.end(); ++it1, i++)
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
      expr_t node = nullptr;
      expr_t aux_var = lhs_expr_t.at(i);
      for (const auto &it : diff_subst_table)
        if (it.second == aux_var)
          {
            node = const_cast<expr_t>(it.first);
            break;
          }

      if (!node)
        {
          cerr << "Unexpected error encountered." << endl;
          exit(EXIT_FAILURE);
        }

      node = node->undiff();
      auto it1 = diff_subst_table.find(node);
      if (it1 == diff_subst_table.end())
        printerr = true;

      if (printerr)
        { // we have undiffed something like diff(x), hence x is not in diff_subst_table
          lhs_expr_t.at(i) = node;
          lhs.at(i) = dynamic_cast<VariableNode *>(node)->symb_id;
        }
      else
        {
          lhs_expr_t.at(i) = const_cast<expr_t>(it1->first);
          lhs.at(i) = const_cast<VariableNode *>(it1->second)->symb_id;
        }
    }
  return lhs;
}

map<pair<string, string>, pair<string, int>>
DynamicModel::walkPacParameters(const string &name)
{
  map<pair<string, string>, pair<string, int>> eqtag_and_lag;

  int i = 0;
  for (auto &equation : equations)
    {
      pair<int, int> lhs(-1, -1);
      pair<int, vector<tuple<int, bool, int>>> ec_params_and_vars;
      set<pair<int, pair<int, int>>> ar_params_and_vars;
      vector<tuple<int, int, int, double>> non_optim_vars_params_and_constants, optim_additive_vars_params_and_constants, additive_vars_params_and_constants;

      if (equation->containsPacExpectation())
        {
          set<pair<int, int>> lhss;
          equation->arg1->collectDynamicVariables(SymbolType::endogenous, lhss);
          lhs = *lhss.begin();
          int lhs_symb_id = lhs.first;
          int lhs_orig_symb_id = lhs_symb_id;
          if (symbol_table.isAuxiliaryVariable(lhs_orig_symb_id))
            try
              {
                lhs_orig_symb_id = symbol_table.getOrigSymbIdForAuxVar(lhs_orig_symb_id);
              }
            catch (...)
              {
              }

          auto arg2 = dynamic_cast<BinaryOpNode *>(equation->arg2);
          if (!arg2)
            {
              cerr << "Pac equation in incorrect format" << endl;
              exit(EXIT_FAILURE);
            }
          auto [optim_share_index, optim_part, non_optim_part, additive_part]
            = arg2->getPacOptimizingShareAndExprNodes(lhs_symb_id, lhs_orig_symb_id);

          if (!optim_part)
            {
              auto bopn = dynamic_cast<BinaryOpNode *>(equation->arg2);
              if (!bopn)
                {
                  cerr << "Error in PAC equation" << endl;
                  exit(EXIT_FAILURE);
                }
              bopn->getPacAREC(lhs_symb_id, lhs_orig_symb_id, ec_params_and_vars, ar_params_and_vars, additive_vars_params_and_constants);
            }
          else
            {
              auto bopn = dynamic_cast<BinaryOpNode *>(optim_part);
              if (!bopn)
                {
                  cerr << "Error in PAC equation" << endl;
                  exit(EXIT_FAILURE);
                }
              bopn->getPacAREC(lhs_symb_id, lhs_orig_symb_id, ec_params_and_vars, ar_params_and_vars, optim_additive_vars_params_and_constants);
              try
                {
                  non_optim_vars_params_and_constants = non_optim_part->matchLinearCombinationOfVariables();
                  if (additive_part)
                    additive_vars_params_and_constants = additive_part->matchLinearCombinationOfVariables();
                }
              catch (ExprNode::MatchFailureException &e)
                {
                  cerr << "Error in parsing non-optimizing agents or additive part of PAC equation: "
                       << e.message << endl;
                  exit(EXIT_FAILURE);
                }
            }

          string eqtag;
          for (auto &tag : equation_tags)
            if (tag.first == (&equation - &equations[0]))
              if (tag.second.first == "name")
                {
                  eqtag = tag.second.second;
                  break;
                }
          if (eqtag.empty())
            {
              cerr << "Every equation with a pac expectation must have been assigned an equation tag name" << endl;
              exit(EXIT_FAILURE);
            }
          if (lhs.first == -1)
            {
              cerr << "walkPacParameters: error obtaining LHS varibale." << endl;
              exit(EXIT_FAILURE);
            }
          if (ec_params_and_vars.second.empty() || ar_params_and_vars.empty())
            {
              cerr << "walkPacParameters: error obtaining RHS parameters." << endl;
              exit(EXIT_FAILURE);
            }
          string eq = "eq" + to_string(i++);
          pac_equation_info[{name, eq}] = {lhs, optim_share_index,
                                           ar_params_and_vars, ec_params_and_vars,
                                           non_optim_vars_params_and_constants,
                                           additive_vars_params_and_constants,
                                           optim_additive_vars_params_and_constants};
          eqtag_and_lag[{name, eqtag}] = {eq, 0};
        }
    }
  return eqtag_and_lag;
}

void
DynamicModel::getPacMaxLag(const string &pac_model_name, map<pair<string, string>, pair<string, int>> &eqtag_and_lag) const
{
  for (auto &equation : equations)
    if (equation->containsPacExpectation(pac_model_name))
      {
        set<pair<int, int>> endogs;
        equation->arg1->collectDynamicVariables(SymbolType::endogenous, endogs);
        if (endogs.size() != 1)
          {
            cerr << "The LHS of the PAC equation may only be comprised of one endogenous variable"
                 << endl;
            exit(EXIT_FAILURE);
          }

        string eqtag;
        for (auto &tag : equation_tags)
          if (tag.first == (&equation - &equations[0]))
            if (tag.second.first == "name")
              {
                eqtag = tag.second.second;
                break;
              }
        string eq = eqtag_and_lag[{pac_model_name, eqtag}].first;
        eqtag_and_lag[{pac_model_name, eqtag}] = {eq, equation->PacMaxLag(endogs.begin()->first)};
      }
}

int
DynamicModel::getPacTargetSymbId(const string &pac_model_name) const
{
  for (auto &equation : equations)
    if (equation->containsPacExpectation(pac_model_name))
      {
        pair<int, int> lhs(-1, -1);
        set<pair<int, int>> lhss;
        equation->arg1->collectDynamicVariables(SymbolType::endogenous, lhss);
        lhs = *lhss.begin();
        int lhs_symb_id = lhs.first;
        int lhs_orig_symb_id = lhs_symb_id;
        if (symbol_table.isAuxiliaryVariable(lhs_symb_id))
          try
            {
              lhs_orig_symb_id = symbol_table.getOrigSymbIdForAuxVar(lhs_symb_id);
            }
          catch (...)
            {
            }
        return equation->arg2->getPacTargetSymbId(lhs_symb_id, lhs_orig_symb_id);
      }
  return -1;
}

void
DynamicModel::declarePacModelConsistentExpectationEndogs(const string &name)
{
  int i = 0;
  for (auto &equation : equations)
    if (equation->containsPacExpectation())
      {
        string eqtag;
        for (auto &tag : equation_tags)
          if (tag.first == (&equation - &equations[0]))
            if (tag.second.first == "name")
              {
                eqtag = tag.second.second;
                break;
              }
        if (eqtag.empty())
          {
            cerr << "Every equation with a pac expectation must have been assigned an equation tag name" << endl;
            exit(EXIT_FAILURE);
          }
        string standard_eqtag = "eq" + to_string(i++);
        try
          {
            pac_mce_z1_symb_ids[{name, standard_eqtag}]
              = symbol_table.addSymbol("mce_Z1_" + name + "_" + standard_eqtag, SymbolType::endogenous);
          }
        catch (SymbolTable::AlreadyDeclaredException &e)
          {
            cerr << "Variable name needed by PAC (mce_Z1_" << name << "_" << standard_eqtag << endl;
            exit(EXIT_FAILURE);
          }
      }
}

void
DynamicModel::addPacModelConsistentExpectationEquation(const string &name, int discount_symb_id,
                                                       const map<pair<string, string>, pair<string, int>> &eqtag_and_lag,
                                                       ExprNode::subst_table_t &diff_subst_table)
{
  int pac_target_symb_id = getPacTargetSymbId(name);
  pac_eqtag_and_lag.insert(eqtag_and_lag.begin(), eqtag_and_lag.end());
  int neqs = 0;
  for (auto &it : eqtag_and_lag)
    {
      string eqtag = it.first.second;
      string standard_eqtag = it.second.first;
      int pac_max_lag_m = it.second.second + 1;
      string append_to_name = name + "_" + standard_eqtag;
      if (pac_mce_z1_symb_ids.find({name, standard_eqtag}) == pac_mce_z1_symb_ids.end())
        {
          cerr << "Error finding pac MCE Z1 symb id" << endl;
          exit(EXIT_FAILURE);
        }
      int mce_z1_symb_id = pac_mce_z1_symb_ids[{name, standard_eqtag}];

      expr_t A = One;
      expr_t fp = Zero;
      expr_t beta = AddVariable(discount_symb_id);
      for (int i = 1; i <= pac_max_lag_m; i++)
        try
          {
            int alpha_i_symb_id = symbol_table.addSymbol("mce_alpha_" + append_to_name + "_" + to_string(i),
                                                         SymbolType::parameter);
            pac_mce_alpha_symb_ids[{name, standard_eqtag}].push_back(alpha_i_symb_id);
            A = AddPlus(A, AddVariable(alpha_i_symb_id));
            fp = AddPlus(fp,
                         AddTimes(AddTimes(AddVariable(alpha_i_symb_id),
                                           AddPower(beta, AddPossiblyNegativeConstant(i))),
                                  AddVariable(mce_z1_symb_id, i)));

          }
        catch (SymbolTable::AlreadyDeclaredException &e)
          {
            cerr << "Variable name needed by PAC (mce_alpha_" << append_to_name << "_" << i << ")" << endl;
            exit(EXIT_FAILURE);
          }

      // Add diff nodes and eqs for pac_target_symb_id
      const VariableNode *target_base_diff_node;
      expr_t diff_node_to_search = AddDiff(AddVariable(pac_target_symb_id));
      if (auto sit = diff_subst_table.find(diff_node_to_search);
          sit != diff_subst_table.end())
        target_base_diff_node = sit->second;
      else
        {
          int symb_id = symbol_table.addDiffAuxiliaryVar(diff_node_to_search->idx, diff_node_to_search);
          target_base_diff_node = AddVariable(symb_id);
          addEquation(dynamic_cast<BinaryOpNode *>(AddEqual(const_cast<VariableNode *>(target_base_diff_node),
                                                            AddMinus(AddVariable(pac_target_symb_id),
                                                                     AddVariable(pac_target_symb_id, -1)))), -1);
          neqs++;
        }

      map<int, VariableNode *> target_aux_var_to_add;
      const VariableNode *last_aux_var = target_base_diff_node;
      for (int i = 1; i <= pac_max_lag_m - 1; i++, neqs++)
        {
          expr_t this_diff_node = AddDiff(AddVariable(pac_target_symb_id, i));
          int symb_id = symbol_table.addDiffLeadAuxiliaryVar(this_diff_node->idx, this_diff_node,
                                                             last_aux_var->symb_id, last_aux_var->lag);
          VariableNode *current_aux_var = AddVariable(symb_id);
          addEquation(dynamic_cast<BinaryOpNode *>(AddEqual(current_aux_var,
                                                            AddVariable(last_aux_var->symb_id, 1))), -1);
          last_aux_var = current_aux_var;
          target_aux_var_to_add[i] = current_aux_var;
        }

      expr_t fs = Zero;
      for (int k = 1; k <= pac_max_lag_m - 1; k++)
        {
          expr_t ssum = Zero;
          for (int j = k+1; j <= pac_max_lag_m; j++)
            {
              int alpha_j_symb_id = -1;
              string varname = "mce_alpha_" + append_to_name + "_" + to_string(j);
              try
                {
                  alpha_j_symb_id = symbol_table.getID(varname);
                }
              catch (SymbolTable::UnknownSymbolNameException &e)
                {
                  alpha_j_symb_id = symbol_table.addSymbol(varname, SymbolType::parameter);
                }
              ssum = AddPlus(ssum,
                             AddTimes(AddVariable(alpha_j_symb_id), AddPower(beta, AddPossiblyNegativeConstant(j))));
            }
          fs = AddPlus(fs, AddTimes(ssum, target_aux_var_to_add[k]));
        }
      addEquation(AddEqual(AddVariable(mce_z1_symb_id),
                           AddMinus(AddTimes(A, AddMinus(const_cast<VariableNode *>(target_base_diff_node), fs)), fp)), -1);
      neqs++;
      pac_expectation_substitution[{name, eqtag}] = AddVariable(mce_z1_symb_id);
    }
  cout << "Pac Model Consistent Expectation: added " << neqs << " auxiliary variables and equations." << endl;
}

void
DynamicModel::fillPacModelInfo(const string &pac_model_name,
                               vector<int> lhs,
                               int max_lag,
                               string aux_model_type,
                               const map<pair<string, string>, pair<string, int>> &eqtag_and_lag,
                               const vector<bool> &nonstationary,
                               expr_t growth)
{
  pac_eqtag_and_lag.insert(eqtag_and_lag.begin(), eqtag_and_lag.end());

  bool stationary_vars_present = any_of(nonstationary.begin(), nonstationary.end(), logical_not<bool>());
  bool nonstationary_vars_present = any_of(nonstationary.begin(), nonstationary.end(), [](bool b) { return b; }); // FIXME: use std::identity instead of an anonymous function when we upgrade to C++20

  int growth_param_index = -1;
  if (growth)
    growth_param_index = symbol_table.addSymbol(pac_model_name
                                                +"_pac_growth_neutrality_correction",
                                                SymbolType::parameter);

  for (auto pac_models_and_eqtags : pac_eqtag_and_lag)
    {
      if (pac_models_and_eqtags.first.first != pac_model_name)
        continue;
      string eqtag = pac_models_and_eqtags.first.second;
      string standard_eqtag = pac_models_and_eqtags.second.first;
      expr_t subExpr = Zero;
      if (stationary_vars_present)
        for (int i = 1; i < max_lag + 1; i++)
          for (auto lhsit : lhs)
            {
              stringstream param_name_h0;
              param_name_h0 << "h0_" << pac_model_name
                            << "_" << standard_eqtag
                            << "_var_" << symbol_table.getName(lhsit)
                            << "_lag_" << i;
              int new_param_symb_id = symbol_table.addSymbol(param_name_h0.str(), SymbolType::parameter);
              pac_h0_indices[{pac_model_name, standard_eqtag}].push_back(new_param_symb_id);
              subExpr = AddPlus(subExpr,
                                AddTimes(AddVariable(new_param_symb_id),
                                         AddVariable(lhsit, -i)));
            }

      if (nonstationary_vars_present)
        for (int i = 1; i < max_lag + 1; i++)
          for (auto lhsit : lhs)
            {
              stringstream param_name_h1;
              param_name_h1 << "h1_" << pac_model_name
                            << "_" << standard_eqtag
                            << "_var_" << symbol_table.getName(lhsit)
                            << "_lag_" << i;
              int new_param_symb_id = symbol_table.addSymbol(param_name_h1.str(), SymbolType::parameter);
              pac_h1_indices[{pac_model_name, standard_eqtag}].push_back(new_param_symb_id);
              subExpr = AddPlus(subExpr,
                                AddTimes(AddVariable(new_param_symb_id),
                                         AddVariable(lhsit, -i)));
            }

      if (growth)
        subExpr = AddPlus(subExpr,
                          AddTimes(AddVariable(growth_param_index), growth));

      pac_expectation_substitution[{pac_model_name, eqtag}] = subExpr;
    }
  pac_model_info[pac_model_name] = {move(lhs), growth_param_index, move(aux_model_type)};
}

void
DynamicModel::substitutePacExpectation(const string &pac_model_name)
{
  for (auto &it : pac_expectation_substitution)
    if (it.first.first == pac_model_name)
      for (auto &equation : equations)
        for (auto & [tagged_eq, tag_pair] : equation_tags)
          if (tagged_eq == (&equation - &equations[0])
              && tag_pair.first == "name" && tag_pair.second == it.first.second)
            {
              auto substeq = dynamic_cast<BinaryOpNode *>(equation->substitutePacExpectation(pac_model_name, it.second));
              assert(substeq);
              equation = substeq;
              break;
            }
}

void
DynamicModel::computingPass(bool jacobianExo, int derivsOrder, int paramsDerivsOrder,
                            const eval_context_t &eval_context, bool no_tmp_terms, bool block, bool use_dll,
                            bool bytecode, bool linear_decomposition)
{
  assert(jacobianExo || (derivsOrder < 2 && paramsDerivsOrder == 0));

  initializeVariablesAndEquations();

  // Prepare for derivation
  computeDerivIDs();

  // Computes dynamic jacobian columns, must be done after computeDerivIDs()
  computeDynJacobianCols(jacobianExo);

  // Compute derivatives w.r. to all endogenous, and possibly exogenous and exogenous deterministic
  set<int> vars;
  for (auto &it : deriv_id_table)
    {
      SymbolType type = symbol_table.getType(it.first.first);
      if (type == SymbolType::endogenous || (jacobianExo && (type == SymbolType::exogenous || type == SymbolType::exogenousDet)))
        vars.insert(it.second);
    }

  // Launch computations
  cout << "Computing " << (linear_decomposition ? "nonlinear " : "")
       << "dynamic model derivatives (order " << derivsOrder << ")." << endl;

  computeDerivatives(derivsOrder, vars);

  if (derivsOrder > 1)
    for (const auto &[indices, d2] : derivatives[2])
      nonzero_hessian_eqs.insert(indices[0]);

  if (paramsDerivsOrder > 0)
    {
      cout << "Computing dynamic model derivatives w.r.t. parameters (order " << paramsDerivsOrder << ")." << endl;
      computeParamsDerivatives(paramsDerivsOrder);
    }

  jacob_map_t contemporaneous_jacobian, static_jacobian;
  map<tuple<int, int, int>, expr_t> first_order_endo_derivatives;
  // for each block contains pair<Size, Feddback_variable>
  vector<pair<int, int>> blocks;
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
  for (auto &equation : equations)
    {
      ExprNode::EquationInfo ei;
      equation->computeXrefs(ei);
      xrefs[i++] = ei;
    }

  i = 0;
  for (auto it = xrefs.begin(); it != xrefs.end(); ++it, i++)
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
  for (const auto &it : eiref)
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
  for (auto it = xrefs.begin(); it != xrefs.end(); ++it, i++)
    {
      output << "M_.xref1.param{" << i << "} = [ ";
      for (const auto &it1 : it->second.param)
        output << symbol_table.getTypeSpecificID(it1.first) + 1 << " ";
      output << "];" << endl;

      output << "M_.xref1.endo{" << i << "} = [ ";
      for (const auto &it1 : it->second.endo)
        output << "struct('id', " << symbol_table.getTypeSpecificID(it1.first) + 1 << ", 'shift', " << it1.second << ");";
      output << "];" << endl;

      output << "M_.xref1.exo{" << i << "} = [ ";
      for (const auto &it1 : it->second.exo)
        output << "struct('id', " << symbol_table.getTypeSpecificID(it1.first) + 1 << ", 'shift', " << it1.second << ");";
      output << "];" << endl;

      output << "M_.xref1.exo_det{" << i << "} = [ ";
      for (const auto &it1 : it->second.exo_det)
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
  for (const auto &it : xrefmap)
    {
      int tsid = symbol_table.getTypeSpecificID(it.first.first) + 1;
      output << "M_.xref2." << type << "{" << tsid << "} = [ ";
      if (last_tsid == tsid)
        output << "M_.xref2." << type << "{" << tsid << "}; ";
      else
        last_tsid = tsid;

      for (const auto &it1 : it.second)
        if (type == "param")
          output << it1 + 1 << " ";
        else
          output << "struct('shift', " << it.first.second << ", 'eq', " << it1+1 << ");";
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
      max_lag = 1;
      max_lead = 1;
      setBlockLeadLag(block, max_lag, max_lead);
    }
  else
    {
      max_lag = getBlockMaxLag(block);
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
                  if (auto its = Derivatives.find({ lag, eq, var, eqr, varr });
                      its != Derivatives.end() && its->second == 2)
                    OK = false;

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
          auto [lag, eq, var, eqr, varr] = it.first;
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
  for (auto & [indices, d1] : derivatives[1])
    {
      int eq = indices[0];
      int var = symbol_table.getTypeSpecificID(getSymbIDByDerivID(indices[1]));
      int lag = getLagByDerivID(indices[1]);
      int block_eq = equation_2_block[eq];
      int block_var = 0;
      derivative_t tmp_derivative;
      lag_var_t lag_var;
      switch (getTypeByDerivID(indices[1]))
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

              if (auto it = block_other_endo_index.find(block_eq);
                  it == block_other_endo_index.end())
                block_other_endo_index[block_eq][var] = 0;
              else
                if (auto it1 = it->second.find(var);
                    it1 == it->second.end())
                  {
                    int size = block_other_endo_index[block_eq].size();
                    block_other_endo_index[block_eq][var] = size;
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

          if (auto it = block_exo_index.find(block_eq);
              it == block_exo_index.end())
            block_exo_index[block_eq][var] = 0;
          else
            if (auto it1 = it->second.find(var);
                it1 == it->second.end())
              {
                int size = block_exo_index[block_eq].size();
                block_exo_index[block_eq][var] = size;
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

          if (auto it = block_det_exo_index.find(block_eq);
              it == block_det_exo_index.end())
            block_det_exo_index[block_eq][var] = 0;
          else
            if (auto it1 = it->second.find(var);
                it1 == it->second.end())
              {
                int size = block_det_exo_index[block_eq].size();
                block_det_exo_index[block_eq][var] = size;
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
DynamicModel::writeDynamicFile(const string &basename, bool block, bool linear_decomposition, bool bytecode, bool use_dll, const string &mexext, const filesystem::path &matlabroot, const filesystem::path &dynareroot, bool julia) const
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
      writeDynamicCFile(basename);
      compileDll(basename, "dynamic", mexext, matlabroot, dynareroot);
    }
  else if (julia)
    writeDynamicJuliaFile(basename);
  else
    writeDynamicMFile(basename);
  writeSetAuxiliaryVariables(basename, julia);
}

void
DynamicModel::writeSetAuxiliaryVariables(const string &basename, bool julia) const
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
  for (auto aux_eq : aux_equations)
    if (auto aux_eq2 = dynamic_cast<ExprNode *>(aux_eq);
        aux_eq2->containsExternalFunction())
      aux_eq2->writeExternalFunctionOutput(output, output_type, temporary_terms,
                                           temporary_terms_idxs, tef_terms);
  for (auto aux_eq : aux_equations)
    {
      dynamic_cast<ExprNode *>(aux_eq)->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
      output << ";" << endl;
    }
}

void
DynamicModel::clearEquations()
{
  equations.clear();
  equations_lineno.clear();
  equation_tags.clear();
  equation_tags_xref.clear();
}

void
DynamicModel::replaceMyEquations(DynamicModel &dynamic_model) const
{
  dynamic_model.clearEquations();

  for (size_t i = 0; i < equations.size(); i++)
    dynamic_model.addEquation(equations[i]->clone(dynamic_model), equations_lineno[i]);

  dynamic_model.equation_tags = equation_tags;
  dynamic_model.equation_tags_xref = equation_tags_xref;
}

void
DynamicModel::computeRamseyPolicyFOCs(const StaticModel &static_model)
{
  // Add aux LM to constraints in equations
  // equation[i]->lhs = rhs becomes equation[i]->MULT_(i+1)*(lhs-rhs) = 0
  int i;
  for (i = 0; i < static_cast<int>(equations.size()); i++)
    {
      auto substeq = dynamic_cast<BinaryOpNode *>(equations[i]->addMultipliersToConstraints(i));
      assert(substeq);
      equations[i] = substeq;
    }
  cout << "Ramsey Problem: added " << i << " Multipliers." << endl;

  // Add Planner Objective to equations so that it appears in Lagrangian
  assert(static_model.equations.size() == 1);
  addEquation(static_model.equations[0]->clone(*this), -1);

  // Get max endo lead and max endo lag
  set<pair<int, int>> dynvars;
  int max_eq_lead = 0;
  int max_eq_lag = 0;
  for (auto &equation : equations)
    equation->collectDynamicVariables(SymbolType::endogenous, dynvars);

  for (const auto &[symb_id, lag] : dynvars)
    {
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
  for (i = 0; i < static_cast<int>(equations.size()); i++)
    for (int lag = -max_eq_lag; lag <= max_eq_lead; lag++)
      {
        expr_t dfpower = nullptr;
        stringstream lagstream;
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

  // Save line numbers and tags, see below
  auto old_equations_lineno = equations_lineno;
  auto old_equation_tags = equation_tags;

  // Prepare derivation of the Lagrangian
  clearEquations();
  addEquation(AddEqual(lagrangian, Zero), -1);
  computeDerivIDs();

  /* Compute Lagrangian derivatives.
     Also restore line numbers and tags for FOCs w.r.t. a Lagrange multiplier
     (i.e. a FOC identical to an equation of the original model) */
  vector<expr_t> neweqs;
  vector<int> neweqs_lineno;
  map<int, vector<pair<string, string>>> neweqs_tags;
  for (auto &[symb_id_and_lag, deriv_id] : deriv_id_table)
    {
      auto &[symb_id, lag] = symb_id_and_lag;
      if (symbol_table.getType(symb_id) == SymbolType::endogenous && lag == 0)
        {
          neweqs.push_back(AddEqual(equations[0]->getNonZeroPartofEquation()->getDerivative(deriv_id), Zero));
          if (int i = symbol_table.getEquationNumberForMultiplier(symb_id);
              i != -1)
            {
              // This is a derivative w.r.t. a Lagrange multiplier
              neweqs_lineno.push_back(old_equations_lineno[i]);
              vector<pair<string, string>> tags;
              for (auto &[j, tagpair] : old_equation_tags)
                if (j == i)
                  tags.emplace_back(tagpair);
              neweqs_tags[neweqs.size()-1] = tags;
            }
          else
            neweqs_lineno.push_back(-1);
        }
    }

  // Overwrite equations with the Lagrangian derivatives
  clearEquations();
  for (size_t i = 0; i < neweqs.size(); i++)
    addEquation(neweqs[i], neweqs_lineno[i], neweqs_tags[i]);
}

void
DynamicModel::toNonlinearPart(DynamicModel &non_linear_equations_dynamic_model) const
{
  // Convert model local variables (need to be done first)
  for (const auto &it : local_variables_table)
    non_linear_equations_dynamic_model.AddLocalVariable(it.first, it.second);
}

bool
DynamicModel::ParamUsedWithLeadLag() const
{
  return ParamUsedWithLeadLagInternal();
}

void
DynamicModel::createVariableMapping(int orig_eq_nbr)
{
  for (int ii = 0; ii < orig_eq_nbr; ii++)
    {
      set<int> eqvars;
      equations[ii]->collectVariables(SymbolType::endogenous, eqvars);
      equations[ii]->collectVariables(SymbolType::exogenous, eqvars);
      for (auto eqvar : eqvars)
        {
          eqvar = symbol_table.getUltimateOrigSymbID(eqvar);
          if (eqvar >= 0 && !symbol_table.isAuxiliaryVariable(eqvar))
            variableMapping[eqvar].emplace(ii);
        }
    }
}

void
DynamicModel::expandEqTags()
{
  set<int> existing_tags;
  for (const auto &eqn : equation_tags)
    if (eqn.second.first == "name")
      existing_tags.insert(eqn.first);

  for (int eq = 0; eq < static_cast<int>(equations.size()); eq++)
    if (existing_tags.find(eq) == existing_tags.end())
      if (auto lhs_expr = dynamic_cast<VariableNode *>(equations[eq]->arg1); lhs_expr && equation_tags_xref.find({ "name", symbol_table.getName(lhs_expr->symb_id)}) == equation_tags_xref.end())
        {
          equation_tags.emplace_back(eq, pair("name", symbol_table.getName(lhs_expr->symb_id)));
          equation_tags_xref.emplace(pair("name", symbol_table.getName(lhs_expr->symb_id)), eq);
        }
      else if (equation_tags_xref.find({ "name", to_string(eq+1) }) == equation_tags_xref.end())
        {
          equation_tags.emplace_back(eq, pair("name", to_string(eq+1)));
          equation_tags_xref.emplace(pair("name", to_string(eq+1)), eq);
        }
      else
        {
          cerr << "Error creating default equation tag: cannot assign default tag to equation number " << eq+1 << " because it is already in use" << endl;
          exit(EXIT_FAILURE);
        }

  sort(equation_tags.begin(), equation_tags.end());
}

set<int>
DynamicModel::findUnusedEndogenous()
{
  set<int> usedEndo, unusedEndo;
  for (auto &equation : equations)
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
  for (auto &equation : equations)
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

  for (auto &equation : equations)
    {
      equation->collectDynamicVariables(SymbolType::endogenous, dynvars);
      equation->collectDynamicVariables(SymbolType::exogenous, dynvars);
      equation->collectDynamicVariables(SymbolType::exogenousDet, dynvars);

      max_lag_with_diffs_expanded_orig = max(equation->maxLagWithDiffsExpanded(),
                                             max_lag_with_diffs_expanded_orig);
    }

  for (const auto &dynvar : dynvars)
    {
      int lag = dynvar.second;
      SymbolType type = symbol_table.getType(dynvar.first);

      max_lead_orig = max(lag, max_lead_orig);
      max_lag_orig = max(-lag, max_lag_orig);

      switch (type)
        {
        case SymbolType::endogenous:
          max_endo_lead_orig = max(lag, max_endo_lead_orig);
          max_endo_lag_orig = max(-lag, max_endo_lag_orig);
          break;
        case SymbolType::exogenous:
          max_exo_lead_orig = max(lag, max_exo_lead_orig);
          max_exo_lag_orig = max(-lag, max_exo_lag_orig);
          break;
        case SymbolType::exogenousDet:
          max_exo_det_lead_orig = max(lag, max_exo_det_lead_orig);
          max_exo_det_lag_orig = max(-lag, max_exo_det_lag_orig);
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

  for (auto &equation : equations)
    equation->collectDynamicVariables(SymbolType::endogenous, dynvars);

  dynJacobianColsNbr = dynvars.size();

  for (auto &equation : equations)
    {
      equation->collectDynamicVariables(SymbolType::exogenous, dynvars);
      equation->collectDynamicVariables(SymbolType::exogenousDet, dynvars);
      equation->collectDynamicVariables(SymbolType::parameter, dynvars);
      equation->collectDynamicVariables(SymbolType::trend, dynvars);
      equation->collectDynamicVariables(SymbolType::logTrend, dynvars);
    }

  for (const auto &dynvar : dynvars)
    {
      int lag = dynvar.second;
      SymbolType type = symbol_table.getType(dynvar.first);

      /* Setting maximum and minimum lags.

         We don't want these to be affected by lead/lags on parameters: they
         are accepted for facilitating variable flipping, but are simply
         ignored. */
      if (type != SymbolType::parameter)
        {
          max_lead = max(lag, max_lead);
          max_lag = max(-lag, max_lag);
        }

      switch (type)
        {
        case SymbolType::endogenous:
          max_endo_lead = max(lag, max_endo_lead);
          max_endo_lag = max(-lag, max_endo_lag);
          break;
        case SymbolType::exogenous:
          max_exo_lead = max(lag, max_exo_lead);
          max_exo_lag = max(-lag, max_exo_lag);
          break;
        case SymbolType::exogenousDet:
          max_exo_det_lead = max(lag, max_exo_det_lead);
          max_exo_det_lag = max(-lag, max_exo_det_lag);
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
  if (deriv_id < 0 || deriv_id >= static_cast<int>(inv_deriv_id_table.size()))
    throw UnknownDerivIDException();

  return inv_deriv_id_table[deriv_id].second;
}

int
DynamicModel::getSymbIDByDerivID(int deriv_id) const noexcept(false)
{
  if (deriv_id < 0 || deriv_id >= static_cast<int>(inv_deriv_id_table.size()))
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

  for (auto &it : deriv_id_table)
    {
      int symb_id = it.first.first;
      int lag = it.first.second;
      int deriv_id = it.second;
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
  for (auto &it : ordered_dyn_endo)
    dyn_jacobian_cols_table[it.second] = sorted_id++;

  // Set final value for dynJacobianColsNbr
  if (jacobianExo)
    dynJacobianColsNbr += symbol_table.exo_nbr() + symbol_table.exo_det_nbr();
}

int
DynamicModel::getDynJacobianCol(int deriv_id) const noexcept(false)
{
  if (auto it = dyn_jacobian_cols_table.find(deriv_id);
      it == dyn_jacobian_cols_table.end())
    throw UnknownDerivIDException();
  else
    return it->second;
}

void
DynamicModel::testTrendDerivativesEqualToZero(const eval_context_t &eval_context)
{
  for (auto &it : deriv_id_table)
    if (symbol_table.getType(it.first.first) == SymbolType::trend
        || symbol_table.getType(it.first.first) == SymbolType::logTrend)
      for (int eq = 0; eq < static_cast<int>(equations.size()); eq++)
        {
          expr_t homogeneq = AddMinus(equations[eq]->arg1,
                                      equations[eq]->arg2);

          // Do not run the test if the term inside the log is zero
          if (fabs(homogeneq->eval(eval_context)) > zero_band)
            {
              expr_t testeq = AddLog(homogeneq); // F = log(lhs-rhs)
              testeq = testeq->getDerivative(it.second); // d F / d Trend
              for (auto &endogit : deriv_id_table)
                if (symbol_table.getType(endogit.first.first) == SymbolType::endogenous)
                  {
                    double nearZero = testeq->getDerivative(endogit.second)->eval(eval_context); // eval d F / d Trend d Endog
                    if (fabs(nearZero) > balanced_growth_test_tol)
                      {
                        cerr << "ERROR: trends not compatible with balanced growth path; the second-order cross partial of equation " << eq + 1 << " (line "
                             << equations_lineno[eq] << ") w.r.t. trend variable "
                             << symbol_table.getName(it.first.first) << " and endogenous variable "
                             << symbol_table.getName(endogit.first.first) << " is not null (abs. value = "
                             << fabs(nearZero) << "). If you are confident that your trends are correctly specified, you can raise the value of option 'balanced_growth_test_tol' in the 'model' block." << endl;
                        exit(EXIT_FAILURE);
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
  ostringstream tt_output; // Used for storing model temp vars and equations
  ostringstream rp_output; // 1st deriv. of residuals w.r.t. parameters
  ostringstream gp_output; // 1st deriv. of Jacobian w.r.t. parameters
  ostringstream rpp_output; // 2nd deriv of residuals w.r.t. parameters
  ostringstream gpp_output; // 2nd deriv of Jacobian w.r.t. parameters
  ostringstream hp_output; // 1st deriv. of Hessian w.r.t. parameters
  ostringstream g3p_output; // 1st deriv. of 3rd deriv. matrix w.r.t. parameters

  temporary_terms_t temp_term_union;
  deriv_node_temp_terms_t tef_terms;

  writeModelLocalVariableTemporaryTerms(temp_term_union, params_derivs_temporary_terms_idxs, tt_output, output_type, tef_terms);
  for (const auto &it : params_derivs_temporary_terms)
    writeTemporaryTerms(it.second, temp_term_union, params_derivs_temporary_terms_idxs, tt_output, output_type, tef_terms);

  for (const auto & [indices, d1] : params_derivatives.find({ 0, 1 })->second)
    {
      auto [eq, param] = vectorToTuple<2>(indices);

      int param_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param)) + 1;

      rp_output << "rp" << LEFT_ARRAY_SUBSCRIPT(output_type) << eq+1 << ", " << param_col
                << RIGHT_ARRAY_SUBSCRIPT(output_type) << " = ";
      d1->writeOutput(rp_output, output_type, temp_term_union, params_derivs_temporary_terms_idxs, tef_terms);
      rp_output << ";" << endl;
    }

  for (const auto & [indices, d2] : params_derivatives.find({ 1, 1 })->second)
    {
      auto [eq, var, param] = vectorToTuple<3>(indices);

      int var_col = getDynJacobianCol(var) + 1;
      int param_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param)) + 1;

      gp_output << "gp" << LEFT_ARRAY_SUBSCRIPT(output_type) << eq+1 << ", " << var_col
                << ", " << param_col << RIGHT_ARRAY_SUBSCRIPT(output_type) << " = ";
      d2->writeOutput(gp_output, output_type, temp_term_union, params_derivs_temporary_terms_idxs, tef_terms);
      gp_output << ";" << endl;
    }

  int i = 1;
  for (const auto &[indices, d2] : params_derivatives.find({ 0, 2 })->second)
    {
      auto [eq, param1, param2] = vectorToTuple<3>(indices);

      int param1_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param1)) + 1;
      int param2_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param2)) + 1;

      rpp_output << "rpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",1"
                 << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << eq+1 << ";" << endl
                 << "rpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",2"
                 << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << param1_col << ";" << endl
                 << "rpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",3"
                 << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << param2_col << ";" << endl
                 << "rpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",4"
                 << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=";
      d2->writeOutput(rpp_output, output_type, temp_term_union, params_derivs_temporary_terms_idxs, tef_terms);
      rpp_output << ";" << endl;

      i++;

      if (param1 != param2)
        {
          // Treat symmetric elements
          rpp_output << "rpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",1"
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
  for (const auto &[indices, d2] : params_derivatives.find({ 1, 2 })->second)
    {
      auto [eq, var, param1, param2] = vectorToTuple<4>(indices);

      int var_col = getDynJacobianCol(var) + 1;
      int param1_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param1)) + 1;
      int param2_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param2)) + 1;

      gpp_output << "gpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",1"
                 << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << eq+1 << ";" << endl
                 << "gpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",2"
                 << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << var_col << ";" << endl
                 << "gpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",3"
                 << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << param1_col << ";" << endl
                 << "gpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",4"
                 << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << param2_col << ";" << endl
                 << "gpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",5"
                 << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=";
      d2->writeOutput(gpp_output, output_type, temp_term_union, params_derivs_temporary_terms_idxs, tef_terms);
      gpp_output << ";" << endl;

      i++;

      if (param1 != param2)
        {
          // Treat symmetric elements
          gpp_output << "gpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",1"
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
  for (const auto &[indices, d2] : params_derivatives.find({ 2, 1 })->second)
    {
      auto [eq, var1, var2, param] = vectorToTuple<4>(indices);

      int var1_col = getDynJacobianCol(var1) + 1;
      int var2_col = getDynJacobianCol(var2) + 1;
      int param_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param)) + 1;

      hp_output << "hp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",1"
                << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << eq+1 << ";" << endl
                << "hp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",2"
                << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << var1_col << ";" << endl
                << "hp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",3"
                << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << var2_col << ";" << endl
                << "hp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",4"
                << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << param_col << ";" << endl
                << "hp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",5"
                << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=";
      d2->writeOutput(hp_output, output_type, temp_term_union, params_derivs_temporary_terms_idxs, tef_terms);
      hp_output << ";" << endl;

      i++;

      if (var1 != var2)
        {
          // Treat symmetric elements
          hp_output << "hp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",1"
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

  i = 1;
  for (const auto &[indices, d2] : params_derivatives.find({ 3, 1 })->second)
    {
      auto [eq, var1, var2, var3, param] = vectorToTuple<5>(indices);

      int var1_col = getDynJacobianCol(var1) + 1;
      int var2_col = getDynJacobianCol(var2) + 1;
      int var3_col = getDynJacobianCol(var3) + 1;
      int param_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param)) + 1;

      g3p_output << "g3p" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",1"
                 << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << eq+1 << ";" << endl
                 << "g3p" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",2"
                 << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << var1_col << ";" << endl
                 << "g3p" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",3"
                 << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << var2_col << ";" << endl
                 << "g3p" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",4"
                 << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << var3_col << ";" << endl
                 << "g3p" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",5"
                 << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << param_col << ";" << endl
                 << "g3p" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",6"
                 << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=";
      d2->writeOutput(g3p_output, output_type, temp_term_union, params_derivs_temporary_terms_idxs, tef_terms);
      g3p_output << ";" << endl;

      i++;
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
      fixNestedParenthesis(rp_output, tmp_paren_vars, message_printed);
      fixNestedParenthesis(gp_output, tmp_paren_vars, message_printed);
      fixNestedParenthesis(rpp_output, tmp_paren_vars, message_printed);
      fixNestedParenthesis(gpp_output, tmp_paren_vars, message_printed);
      fixNestedParenthesis(hp_output, tmp_paren_vars, message_printed);
      fixNestedParenthesis(g3p_output, tmp_paren_vars, message_printed);
      paramsDerivsFile << "function [rp, gp, rpp, gpp, hp, g3p] = dynamic_params_derivs(y, x, params, steady_state, it_, ss_param_deriv, ss_param_2nd_deriv)" << endl
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
                       << "%   g3p      [#first_order_g3_terms by 6] double   Jacobian matrix of derivatives of g3 (dynamic 3rd derivs) with respect to the parameters;" << endl
                       << "%                                                              rows: respective derivative term" << endl
                       << "%                                                              1st column: equation number of the term appearing" << endl
                       << "%                                                              2nd column: column number of first variable in g3 of the dynamic model" << endl
                       << "%                                                              3rd column: column number of second variable in g3 of the dynamic model" << endl
                       << "%                                                              4th column: column number of third variable in g3 of the dynamic model" << endl
                       << "%                                                              5th column: number of the parameter in derivative" << endl
                       << "%                                                              6th column: value of the Hessian term" << endl
                       << "%" << endl
                       << "%" << endl
                       << "% Warning : this file is generated automatically by Dynare" << endl
                       << "%           from model file (.mod)" << endl << endl
                       << "T = NaN(" << params_derivs_temporary_terms_idxs.size() << ",1);" << endl
                       << tt_output.str()
                       << "rp = zeros(" << equations.size() << ", "
                       << symbol_table.param_nbr() << ");" << endl
                       << rp_output.str()
                       << "gp = zeros(" << equations.size() << ", " << dynJacobianColsNbr << ", " << symbol_table.param_nbr() << ");" << endl
                       << gp_output.str()
                       << "if nargout >= 3" << endl
                       << "rpp = zeros(" << params_derivatives.find({ 0, 2 })->second.size() << ",4);" << endl
                       << rpp_output.str()
                       << "gpp = zeros(" << params_derivatives.find({ 1, 2 })->second.size() << ",5);" << endl
                       << gpp_output.str()
                       << "end" << endl
                       << "if nargout >= 5" << endl
                       << "hp = zeros(" << params_derivatives.find({ 2, 1 })->second.size() << ",5);" << endl
                       << hp_output.str()
                       << "end" << endl
                       << "if nargout >= 6" << endl
                       << "g3p = zeros(" << params_derivatives.find({ 3, 1 })->second.size() << ",6);" << endl
                       << g3p_output.str()
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
                     << rp_output.str()
                     << "gp = zeros(" << equations.size() << ", " << dynJacobianColsNbr << ", " << symbol_table.param_nbr() << ");" << endl
                     << gp_output.str()
                     << "rpp = zeros(" << params_derivatives.find({ 0, 2 })->second.size() << ",4);" << endl
                     << rpp_output.str()
                     << "gpp = zeros(" << params_derivatives.find({ 1, 2 })->second.size() << ",5);" << endl
                     << gpp_output.str()
                     << "hp = zeros(" << params_derivatives.find({ 2, 1 })->second.size() << ",5);" << endl
                     << hp_output.str()
                     << "g3p = zeros(" << params_derivatives.find({ 3, 1 })->second.size() << ",6);" << endl
                     << g3p_output.str()
                     << "(rp, gp, rpp, gpp, hp, g3p)" << endl
                     << "end" << endl
                     << "end" << endl;

  paramsDerivsFile.close();
}

void
DynamicModel::writeLatexFile(const string &basename, bool write_equation_tags) const
{
  writeLatexModelFile(basename, "dynamic", ExprNodeOutputType::latexDynamicModel, write_equation_tags);
}

void
DynamicModel::writeLatexOriginalFile(const string &basename, bool write_equation_tags) const
{
  writeLatexModelFile(basename, "original", ExprNodeOutputType::latexDynamicModel, write_equation_tags);
}

void
DynamicModel::substituteEndoLeadGreaterThanTwo(bool deterministic_model)
{
  substituteLeadLagInternal(AuxVarType::endoLead, deterministic_model, {});
}

void
DynamicModel::substituteEndoLagGreaterThanTwo(bool deterministic_model)
{
  substituteLeadLagInternal(AuxVarType::endoLag, deterministic_model, {});
}

void
DynamicModel::substituteExoLead(bool deterministic_model)
{
  substituteLeadLagInternal(AuxVarType::exoLead, deterministic_model, {});
}

void
DynamicModel::substituteExoLag(bool deterministic_model)
{
  substituteLeadLagInternal(AuxVarType::exoLag, deterministic_model, {});
}

void
DynamicModel::substituteLeadLagInternal(AuxVarType type, bool deterministic_model, const vector<string> &subset)
{
  ExprNode::subst_table_t subst_table;
  vector<BinaryOpNode *> neweqs;

  // Substitute in used model local variables
  set<int> used_local_vars;
  for (auto &equation : equations)
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
  for (auto &equation : equations)
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
      auto substeq = dynamic_cast<BinaryOpNode *>(subst);
      assert(substeq);
      equation = substeq;
    }

  // Add new equations
  for (auto &neweq : neweqs)
    {
      addEquation(neweq, -1);
      aux_equations.push_back(neweq);
    }

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
  for (auto &equation : equations)
    equation = dynamic_cast<BinaryOpNode *>(equation->substituteAdl());
}

vector<int>
DynamicModel::getEquationNumbersFromTags(const set<string> &eqtags) const
{
  vector<int> eqnumbers;
  for (auto &eqtag : eqtags)
    {
      bool found = false;
      for (const auto &equation_tag : equation_tags)
        if (equation_tag.second.first == "name"
            && equation_tag.second.second == eqtag)
          {
            found = true;
            eqnumbers.push_back(equation_tag.first);
            break;
          }
      if (!found)
        {
          cerr << "ERROR: looking for equation tag " << eqtag << " failed." << endl;
          exit(EXIT_FAILURE);
        }
    }
  return eqnumbers;
}

void
DynamicModel::findPacExpectationEquationNumbers(vector<int> &eqnumbers) const
{
  int i = 0;
  for (auto &equation : equations)
    {
      if (equation->containsPacExpectation()
          && find(eqnumbers.begin(), eqnumbers.end(), i) == eqnumbers.end())
        eqnumbers.push_back(i);
      i++;
    }
}

pair<lag_equivalence_table_t, ExprNode::subst_table_t>
DynamicModel::substituteUnaryOps()
{
  vector<int> eqnumbers(equations.size());
  iota(eqnumbers.begin(), eqnumbers.end(), 0);
  return substituteUnaryOps(eqnumbers);
}

pair<lag_equivalence_table_t, ExprNode::subst_table_t>
DynamicModel::substituteUnaryOps(const set<string> &var_model_eqtags)
{
  vector<int> eqnumbers = getEquationNumbersFromTags(var_model_eqtags);
  findPacExpectationEquationNumbers(eqnumbers);
  return substituteUnaryOps(eqnumbers);
}

pair<lag_equivalence_table_t, ExprNode::subst_table_t>
DynamicModel::substituteUnaryOps(const vector<int> &eqnumbers)
{
  lag_equivalence_table_t nodes;
  ExprNode::subst_table_t subst_table;

  // Mark unary ops to be substituted in model local variables that appear in selected equations
  set<int> used_local_vars;
  for (int eqnumber : eqnumbers)
    equations[eqnumber]->collectVariables(SymbolType::modelLocalVariable, used_local_vars);
  for (auto &it : local_variables_table)
    if (used_local_vars.find(it.first) != used_local_vars.end())
      it.second->findUnaryOpNodesForAuxVarCreation(nodes);

  // Mark unary ops to be substituted in selected equations
  for (int eqnumber : eqnumbers)
    equations[eqnumber]->findUnaryOpNodesForAuxVarCreation(nodes);

  // Substitute in model local variables
  vector<BinaryOpNode *> neweqs;
  for (auto &it : local_variables_table)
    it.second = it.second->substituteUnaryOpNodes(nodes, subst_table, neweqs);

  // Substitute in equations
  for (auto &equation : equations)
    {
      auto substeq = dynamic_cast<BinaryOpNode *>(equation->
                                                  substituteUnaryOpNodes(nodes, subst_table, neweqs));
      assert(substeq);
      equation = substeq;
    }

  // Add new equations
  for (auto &neweq : neweqs)
    {
      addEquation(neweq, -1);
      aux_equations.push_back(neweq);
    }

  if (subst_table.size() > 0)
    cout << "Substitution of Unary Ops: added " << neweqs.size() << " auxiliary variables and equations." << endl;

  return { nodes, subst_table };
}

pair<lag_equivalence_table_t, ExprNode::subst_table_t>
DynamicModel::substituteDiff(vector<expr_t> &pac_growth)
{
  /* Note: at this point, we know that there is no diff operator with a lead,
     because they have been expanded by DataTree::AddDiff().
     Hence we can go forward with the substitution without worrying about the
     expectation operator. */

  lag_equivalence_table_t diff_nodes;
  ExprNode::subst_table_t diff_subst_table;

  // Mark diff operators to be substituted in model local variables
  set<int> used_local_vars;
  for (const auto &equation : equations)
    equation->collectVariables(SymbolType::modelLocalVariable, used_local_vars);
  for (auto &it : local_variables_table)
    if (used_local_vars.find(it.first) != used_local_vars.end())
      it.second->findDiffNodes(diff_nodes);

  // Mark diff operators to be substituted in equations
  for (const auto &equation : equations)
    equation->findDiffNodes(diff_nodes);

  for (const auto &gv : pac_growth)
    if (gv)
      gv->findDiffNodes(diff_nodes);

  // Substitute in model local variables
  vector<BinaryOpNode *> neweqs;
  for (auto &it : local_variables_table)
    it.second = it.second->substituteDiff(diff_nodes, diff_subst_table, neweqs);

  // Substitute in equations
  for (auto &equation : equations)
    {
      auto substeq = dynamic_cast<BinaryOpNode *>(equation->
                                                  substituteDiff(diff_nodes, diff_subst_table, neweqs));
      assert(substeq);
      equation = substeq;
    }

  for (auto &it : pac_growth)
    if (it)
      it = it->substituteDiff(diff_nodes, diff_subst_table, neweqs);

  // Add new equations
  for (auto &neweq : neweqs)
    {
      addEquation(neweq, -1);
      aux_equations.push_back(neweq);
    }

  if (diff_subst_table.size() > 0)
    cout << "Substitution of Diff operator: added " << neweqs.size() << " auxiliary variables and equations." << endl;

  return { diff_nodes, diff_subst_table };
}

void
DynamicModel::substituteExpectation(bool partial_information_model)
{
  ExprNode::subst_table_t subst_table;
  vector<BinaryOpNode *> neweqs;

  // Substitute in model local variables
  for (auto &it : local_variables_table)
    it.second = it.second->substituteExpectation(subst_table, neweqs, partial_information_model);

  // Substitute in equations
  for (auto &equation : equations)
    {
      auto substeq = dynamic_cast<BinaryOpNode *>(equation->substituteExpectation(subst_table, neweqs, partial_information_model));
      assert(substeq);
      equation = substeq;
    }

  // Add new equations
  for (auto &neweq : neweqs)
    {
      addEquation(neweq, -1);
      aux_equations.push_back(neweq);
    }

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
  for (auto &it : local_variables_table)
    it.second = it.second->decreaseLeadsLagsPredeterminedVariables();

  for (auto &equation : equations)
    {
      auto substeq = dynamic_cast<BinaryOpNode *>(equation->decreaseLeadsLagsPredeterminedVariables());
      assert(substeq);
      equation = substeq;
    }
}

void
DynamicModel::detrendEquations()
{
  // We go backwards in the list of trend_vars, to deal correctly with I(2) processes
  for (auto it = nonstationary_symbols_map.crbegin();
       it != nonstationary_symbols_map.crend(); ++it)
    for (auto &equation : equations)
      {
        auto substeq = dynamic_cast<BinaryOpNode *>(equation->detrend(it->first, it->second.first, it->second.second));
        assert(substeq);
        equation = dynamic_cast<BinaryOpNode *>(substeq);
      }

  for (auto &equation : equations)
    {
      auto substeq = dynamic_cast<BinaryOpNode *>(equation->removeTrendLeadLag(trend_symbols_map));
      assert(substeq);
      equation = dynamic_cast<BinaryOpNode *>(substeq);
    }
}

void
DynamicModel::removeTrendVariableFromEquations()
{
  for (auto &equation : equations)
    {
      auto substeq = dynamic_cast<BinaryOpNode *>(equation->replaceTrendVar());
      assert(substeq);
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
      assert(aux_equation->op_code == BinaryOpcode::equal);
      auto auxvar = dynamic_cast<VariableNode *>(aux_equation->arg1);
      assert(auxvar);
      try
        {
          double val = aux_equation->arg2->eval(eval_context);
          eval_context[auxvar->symb_id] = val;
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
  for (int &trendVar : trendVars)
    eval_context[trendVar] = 2; //not <= 0 bc of log, not 1 bc of powers
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
  auto beq = dynamic_cast<BinaryOpNode *>(eq);
  assert(beq && beq->op_code == BinaryOpcode::equal);

  vector<pair<string, string>> soe_eq_tags;
  for (const auto &eq_tag : eq_tags)
    soe_eq_tags.push_back(eq_tag);

  int n = static_only_equations.size();
  static_only_equations.push_back(beq);
  static_only_equations_lineno.push_back(lineno);
  static_only_equations_equation_tags.push_back(soe_eq_tags);
  for (auto &it : soe_eq_tags)
    static_only_equation_tags_xref.emplace(it, n);
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

  for (const auto &equation_tag : equation_tags)
    if (equation_tag.second.first == "dynamic")
      eqs.insert(equation_tag.first);

  return eqs.size();
}

bool
DynamicModel::isChecksumMatching(const string &basename, bool block) const
{
  stringstream buffer;

  // Write equation tags
  for (const auto &equation_tag : equation_tags)
    buffer << "  " << equation_tag.first + 1
           << equation_tag.second.first
           << equation_tag.second.second << endl;

  ExprNodeOutputType buffer_type = block ? ExprNodeOutputType::matlabDynamicModelSparse : ExprNodeOutputType::CDynamicModel;

  deriv_node_temp_terms_t tef_terms;
  temporary_terms_t temp_term_union;
  writeModelLocalVariableTemporaryTerms(temp_term_union, temporary_terms_idxs,
                                        buffer, buffer_type, tef_terms);

  writeTemporaryTerms(temporary_terms_derivatives[0],
                      temp_term_union, temporary_terms_idxs,
                      buffer, buffer_type, tef_terms);

  writeModelEquations(buffer, buffer_type, temp_term_union);

  size_t result = hash<string>{}(buffer.str());

  // check whether basename directory exist. If not, create it.
  // If it does, read old checksum if it exists, return if equal to result
  fstream checksum_file;
  auto filename = filesystem::path{basename} / "checksum";
  if (!filesystem::create_directory(basename))
    {
      checksum_file.open(filename, ios::in | ios::binary);
      if (checksum_file.is_open())
        {
          size_t old_checksum;
          checksum_file >> old_checksum;
          checksum_file.close();
          if (old_checksum == result)
            return true;
        }
    }

  // write new checksum file if none or different from old checksum
  checksum_file.open(filename, ios::out | ios::binary);
  if (!checksum_file.is_open())
    {
      cerr << "ERROR: Can't open file " << filename << endl;
      exit(EXIT_FAILURE);
    }
  checksum_file << result;
  checksum_file.close();
  return false;
}

void
DynamicModel::writeJsonOutput(ostream &output) const
{
  writeJsonModelEquations(output, false);
  output << ", ";
  writeJsonXrefs(output);
  output << ", ";
  writeJsonAST(output);
  output << ", ";
  writeJsonVariableMapping(output);
}

void
DynamicModel::writeJsonAST(ostream &output) const
{
  vector<pair<string, string>> eqtags;
  output << R"("abstract_syntax_tree":[)" << endl;
  for (int eq = 0; eq < static_cast<int>(equations.size()); eq++)
    {
      if (eq != 0)
        output << ", ";

      output << R"({ "number":)" << eq
             << R"(, "line":)" << equations_lineno[eq];

      for (const auto &equation_tag : equation_tags)
        if (equation_tag.first == eq)
          eqtags.push_back(equation_tag.second);

      if (!eqtags.empty())
        {
          output << R"(, "tags": {)";
          int i = 0;
          for (const auto &[name, value] : eqtags)
            {
              if (i != 0)
                output << ", ";
              output << R"(")" << name << R"(": ")" << value << R"(")";
              i++;
            }
          output << "}";
          eqtags.clear();
        }

      output << R"(, "AST": )";
      equations[eq]->writeJsonAST(output);
      output << "}";
    }
  output << "]";
}

void
DynamicModel::writeJsonVariableMapping(ostream &output) const
{
  output << R"("variable_mapping":[)" << endl;
  int ii = 0;
  int end_idx_map = static_cast<int>(variableMapping.size()-1);
  for (const auto &variable : variableMapping)
    {
      output << R"({"name": ")" << symbol_table.getName(variable.first) << R"(", "equations":[)";
      int it = 0;
      int end_idx_eq = static_cast<int>(variable.second.size())-1;
      for (const auto &equation : variable.second)
        for (const auto &equation_tag : equation_tags)
          if (equation_tag.first == equation && equation_tag.second.first == "name")
            output << R"(")" << equation_tag.second.second << (it++ == end_idx_eq ? R"("])" : R"(", )");
      output << (ii++ == end_idx_map ? R"(})" : R"(},)") << endl;
    }
  output << "]";
}

void
DynamicModel::writeJsonXrefsHelper(ostream &output, const map<pair<int, int>, set<int>> &xrefs) const
{
  for (auto it = xrefs.begin(); it != xrefs.end(); ++it)
    {
      if (it != xrefs.begin())
        output << ", ";
      output << R"({"name": ")" << symbol_table.getName(it->first.first) << R"(")"
             << R"(, "shift": )" << it->first.second
             << R"(, "equations": [)";
      for (auto it1 = it->second.begin(); it1 != it->second.end(); ++it1)
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
  output << R"("xrefs": {)"
         << R"("parameters": [)";
  writeJsonXrefsHelper(output, xref_param);
  output << "]"
         << R"(, "endogenous": [)";
  writeJsonXrefsHelper(output, xref_endo);
  output << "]"
         << R"(, "exogenous": [)";
  writeJsonXrefsHelper(output, xref_exo);
  output << "]"
         << R"(, "exogenous_deterministic": [)";
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
  output << R"("model_info": {)"
         << R"("lead_lag_incidence": [)";
  // Loop on endogenous variables
  int nstatic = 0,
    nfwrd = 0,
    npred = 0,
    nboth = 0;
  for (int endoID = 0; endoID < symbol_table.endo_nbr(); endoID++)
    {
      if (endoID != 0)
        output << ",";
      output << "[";
      int sstatic = 1,
        sfwrd = 0,
        spred = 0,
        sboth = 0;
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
      nfwrd += sfwrd;
      npred += spred;
      nboth += sboth;
      output << "]";
    }
  output << "], "
         << R"("nstatic": )" << nstatic << ", "
         << R"("nfwrd": )" << nfwrd << ", "
         << R"("npred": )" << npred << ", "
         << R"("nboth": )" << nboth << ", "
         << R"("nsfwrd": )" << nfwrd+nboth << ", "
         << R"("nspred": )" << npred+nboth << ", "
         << R"("ndynamic": )" << npred+nboth+nfwrd << endl;
  output << "}";
}

void
DynamicModel::writeJsonComputingPassOutput(ostream &output, bool writeDetails) const
{
  ostringstream model_local_vars_output; // Used for storing model local vars
  vector<ostringstream> d_output(derivatives.size()); // Derivatives output (at all orders, including 0=residual)

  deriv_node_temp_terms_t tef_terms;
  temporary_terms_t temp_term_union;

  writeJsonModelLocalVariables(model_local_vars_output, tef_terms);

  writeJsonTemporaryTerms(temporary_terms_derivatives[0], temp_term_union, d_output[0], tef_terms, "");
  d_output[0] << ", ";
  writeJsonModelEquations(d_output[0], true);

  int ncols = dynJacobianColsNbr;
  for (size_t i = 1; i < derivatives.size(); i++)
    {
      string matrix_name = i == 1 ? "jacobian" : i == 2 ? "hessian" : i == 3 ? "third_derivative" : to_string(i) + "th_derivative";
      writeJsonTemporaryTerms(temporary_terms_derivatives[i], temp_term_union, d_output[i], tef_terms, matrix_name);
      temp_term_union.insert(temporary_terms_derivatives[i].begin(), temporary_terms_derivatives[i].end());
      d_output[i] << R"(, ")" << matrix_name  << R"(": {)"
                  << R"(  "nrows": )" << equations.size()
                  << R"(, "ncols": )" << ncols
                  << R"(, "entries": [)";

      for (auto it = derivatives[i].begin(); it != derivatives[i].end(); ++it)
        {
          if (it != derivatives[i].begin())
            d_output[i] << ", ";

          const vector<int> &vidx = it->first;
          expr_t d = it->second;
          int eq = vidx[0];

          int col_idx = 0;
          for (size_t j = 1; j < vidx.size(); j++)
            {
              col_idx *= dynJacobianColsNbr;
              col_idx += getDynJacobianCol(vidx[j]);
            }

          if (writeDetails)
            d_output[i] << R"({"eq": )" << eq + 1;
          else
            d_output[i] << R"({"row": )" << eq + 1;

          d_output[i] << R"(, "col": )" << (i > 1 ? "[" : "") << col_idx + 1;

          if (i == 2 && vidx[1] != vidx[2]) // Symmetric elements in hessian
            {
              int col_idx_sym = getDynJacobianCol(vidx[2]) * dynJacobianColsNbr + getDynJacobianCol(vidx[1]);
              d_output[i] << ", " << col_idx_sym + 1;
            }
          if (i > 1)
            d_output[i] << "]";

          if (writeDetails)
            for (size_t j = 1; j < vidx.size(); j++)
              d_output[i] << R"(, "var)" << (i > 1 ? to_string(j) : "") << R"(": ")" << symbol_table.getName(getSymbIDByDerivID(vidx[j])) << R"(")"
                          << R"(, "shift)" << (i > 1 ? to_string(j) : "") << R"(": )" << getLagByDerivID(vidx[j]);

          d_output[i] << R"(, "val": ")";
          d->writeJsonOutput(d_output[i], temp_term_union, tef_terms);
          d_output[i] << R"("})" << endl;
        }
      d_output[i] << "]}";

      ncols *= dynJacobianColsNbr;
    }

  if (writeDetails)
    output << R"("dynamic_model": {)";
  else
    output << R"("dynamic_model_simple": {)";
  output << model_local_vars_output.str();
  for (const auto &it : d_output)
    output << ", " << it.str();
  output << "}";
}

void
DynamicModel::writeJsonParamsDerivativesFile(ostream &output, bool writeDetails) const
{
  if (!params_derivatives.size())
    return;

  ostringstream model_local_vars_output; // Used for storing model local vars
  ostringstream model_output; // Used for storing model temp vars and equations
  ostringstream rp_output; // 1st deriv. of residuals w.r.t. parameters
  ostringstream gp_output; // 1st deriv. of Jacobian w.r.t. parameters
  ostringstream rpp_output; // 2nd deriv of residuals w.r.t. parameters
  ostringstream gpp_output; // 2nd deriv of Jacobian w.r.t. parameters
  ostringstream hp_output; // 1st deriv. of Hessian w.r.t. parameters
  ostringstream g3p_output; // 1st deriv. of 3rd deriv. matrix w.r.t. parameters

  deriv_node_temp_terms_t tef_terms;
  writeJsonModelLocalVariables(model_local_vars_output, tef_terms);

  temporary_terms_t temp_term_union;
  for (const auto &it : params_derivs_temporary_terms)
    writeJsonTemporaryTerms(it.second, temp_term_union, model_output, tef_terms, "all");

  rp_output << R"("deriv_wrt_params": {)"
            << R"(  "neqs": )" << equations.size()
            << R"(, "nparamcols": )" << symbol_table.param_nbr()
            << R"(, "entries": [)";
  auto &rp = params_derivatives.find({ 0, 1 })->second;
  for (auto it = rp.begin(); it != rp.end(); ++it)
    {
      if (it != rp.begin())
        rp_output << ", ";

      auto [eq, param] = vectorToTuple<2>(it->first);
      expr_t d1 = it->second;

      int param_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param)) + 1;

      if (writeDetails)
        rp_output << R"({"eq": )" << eq + 1;
      else
        rp_output << R"({"row": )" << eq + 1;

      rp_output << R"(, "param_col": )" << param_col + 1;

      if (writeDetails)
        rp_output << R"(, "param": ")" << symbol_table.getName(getSymbIDByDerivID(param)) << R"(")";

      rp_output << R"(, "val": ")";
      d1->writeJsonOutput(rp_output, temp_term_union, tef_terms);
      rp_output << R"("})" << endl;
    }
  rp_output << "]}";

  gp_output << R"("deriv_jacobian_wrt_params": {)"
            << R"(  "neqs": )" << equations.size()
            << R"(, "nvarcols": )" << dynJacobianColsNbr
            << R"(, "nparamcols": )" << symbol_table.param_nbr()
            << R"(, "entries": [)";
  auto &gp = params_derivatives.find({ 1, 1 })->second;
  for (auto it = gp.begin(); it != gp.end(); ++it)
    {
      if (it != gp.begin())
        gp_output << ", ";

      auto [eq, var, param] = vectorToTuple<3>(it->first);
      expr_t d2 = it->second;

      int var_col = getDynJacobianCol(var) + 1;
      int param_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param)) + 1;

      if (writeDetails)
        gp_output << R"({"eq": )" << eq + 1;
      else
        gp_output << R"({"row": )" << eq + 1;

      gp_output << R"(, "var_col": )" << var_col + 1
                << R"(, "param_col": )" << param_col + 1;

      if (writeDetails)
        gp_output << R"(, "var": ")" << symbol_table.getName(getSymbIDByDerivID(var)) << R"(")"
                  << R"(, "lag": )" << getLagByDerivID(var)
                  << R"(, "param": ")" << symbol_table.getName(getSymbIDByDerivID(param)) << R"(")";

      gp_output << R"(, "val": ")";
      d2->writeJsonOutput(gp_output, temp_term_union, tef_terms);
      gp_output << R"("})" << endl;
    }
  gp_output << "]}";

  rpp_output << R"("second_deriv_residuals_wrt_params": {)"
             << R"(  "nrows": )" << equations.size()
             << R"(, "nparam1cols": )" << symbol_table.param_nbr()
             << R"(, "nparam2cols": )" << symbol_table.param_nbr()
             << R"(, "entries": [)";
  auto &rpp = params_derivatives.find({ 0, 2 })->second;
  for (auto it = rpp.begin(); it != rpp.end(); ++it)
    {
      if (it != rpp.begin())
        rpp_output << ", ";

      auto [eq, param1, param2] = vectorToTuple<3>(it->first);
      expr_t d2 = it->second;

      int param1_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param1)) + 1;
      int param2_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param2)) + 1;

      if (writeDetails)
        rpp_output << R"({"eq": )" << eq + 1;
      else
        rpp_output << R"({"row": )" << eq + 1;
      rpp_output << R"(, "param1_col": )" << param1_col + 1
                 << R"(, "param2_col": )" << param2_col + 1;

      if (writeDetails)
        rpp_output << R"(, "param1": ")" << symbol_table.getName(getSymbIDByDerivID(param1)) << R"(")"
                   << R"(, "param2": ")" << symbol_table.getName(getSymbIDByDerivID(param2)) << R"(")";

      rpp_output << R"(, "val": ")";
      d2->writeJsonOutput(rpp_output, temp_term_union, tef_terms);
      rpp_output << R"("})" << endl;
    }
  rpp_output << "]}";

  gpp_output << R"("second_deriv_jacobian_wrt_params": {)"
             << R"(  "neqs": )" << equations.size()
             << R"(, "nvarcols": )" << dynJacobianColsNbr
             << R"(, "nparam1cols": )" << symbol_table.param_nbr()
             << R"(, "nparam2cols": )" << symbol_table.param_nbr()
             << R"(, "entries": [)";
  auto &gpp = params_derivatives.find({ 1, 2 })->second;
  for (auto it = gpp.begin(); it != gpp.end(); ++it)
    {
      if (it != gpp.begin())
        gpp_output << ", ";

      auto [eq, var, param1, param2] = vectorToTuple<4>(it->first);
      expr_t d2 = it->second;

      int var_col = getDynJacobianCol(var) + 1;
      int param1_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param1)) + 1;
      int param2_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param2)) + 1;

      if (writeDetails)
        gpp_output << R"({"eq": )" << eq + 1;
      else
        gpp_output << R"({"row": )" << eq + 1;

      gpp_output << R"(, "var_col": )" << var_col + 1
                 << R"(, "param1_col": )" << param1_col + 1
                 << R"(, "param2_col": )" << param2_col + 1;

      if (writeDetails)
        gpp_output << R"(, "var": ")" << symbol_table.getName(var) << R"(")"
                   << R"(, "lag": )" << getLagByDerivID(var)
                   << R"(, "param1": ")" << symbol_table.getName(getSymbIDByDerivID(param1)) << R"(")"
                   << R"(, "param2": ")" << symbol_table.getName(getSymbIDByDerivID(param2)) << R"(")";

      gpp_output << R"(, "val": ")";
      d2->writeJsonOutput(gpp_output, temp_term_union, tef_terms);
      gpp_output << R"("})" << endl;
    }
  gpp_output << "]}" << endl;

  hp_output << R"("derivative_hessian_wrt_params": {)"
            << R"(  "neqs": )" << equations.size()
            << R"(, "nvar1cols": )" << dynJacobianColsNbr
            << R"(, "nvar2cols": )" << dynJacobianColsNbr
            << R"(, "nparamcols": )" << symbol_table.param_nbr()
            << R"(, "entries": [)";
  auto &hp = params_derivatives.find({ 2, 1 })->second;
  for (auto it = hp.begin(); it != hp.end(); ++it)
    {
      if (it != hp.begin())
        hp_output << ", ";

      auto [eq, var1, var2, param] = vectorToTuple<4>(it->first);
      expr_t d2 = it->second;

      int var1_col = getDynJacobianCol(var1) + 1;
      int var2_col = getDynJacobianCol(var2) + 1;
      int param_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param)) + 1;

      if (writeDetails)
        hp_output << R"({"eq": )" << eq + 1;
      else
        hp_output << R"({"row": )" << eq + 1;

      hp_output << R"(, "var1_col": )" << var1_col + 1
                << R"(, "var2_col": )" << var2_col + 1
                << R"(, "param_col": )" << param_col + 1;

      if (writeDetails)
        hp_output << R"(, "var1": ")" << symbol_table.getName(getSymbIDByDerivID(var1)) << R"(")"
                  << R"(, "lag1": )" << getLagByDerivID(var1)
                  << R"(, "var2": ")" << symbol_table.getName(getSymbIDByDerivID(var2)) << R"(")"
                  << R"(, "lag2": )" << getLagByDerivID(var2)
                  << R"(, "param": ")" << symbol_table.getName(getSymbIDByDerivID(param)) << R"(")";

      hp_output << R"(, "val": ")";
      d2->writeJsonOutput(hp_output, temp_term_union, tef_terms);
      hp_output << R"("})" << endl;
    }
  hp_output << "]}" << endl;

  g3p_output << R"("derivative_g3_wrt_params": {)"
             << R"(  "neqs": )" << equations.size()
             << R"(, "nvar1cols": )" << dynJacobianColsNbr
             << R"(, "nvar2cols": )" << dynJacobianColsNbr
             << R"(, "nvar3cols": )" << dynJacobianColsNbr
             << R"(, "nparamcols": )" << symbol_table.param_nbr()
             << R"(, "entries": [)";
  auto &g3p = params_derivatives.find({ 3, 1 })->second;
  for (auto it = g3p.begin(); it != g3p.end(); ++it)
    {
      if (it != g3p.begin())
        g3p_output << ", ";

      auto [eq, var1, var2, var3, param] = vectorToTuple<5>(it->first);
      expr_t d2 = it->second;

      int var1_col = getDynJacobianCol(var1) + 1;
      int var2_col = getDynJacobianCol(var2) + 1;
      int var3_col = getDynJacobianCol(var3) + 1;
      int param_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param)) + 1;

      if (writeDetails)
        g3p_output << R"({"eq": )" << eq + 1;
      else
        g3p_output << R"({"row": )" << eq + 1;

      g3p_output << R"(, "var1_col": )" << var1_col + 1
                 << R"(, "var2_col": )" << var2_col + 1
                 << R"(, "var3_col": )" << var3_col + 1
                 << R"(, "param_col": )" << param_col + 1;

      if (writeDetails)
        g3p_output << R"(, "var1": ")" << symbol_table.getName(getSymbIDByDerivID(var1)) << R"(")"
                   << R"(, "lag1": )" << getLagByDerivID(var1)
                   << R"(, "var2": ")" << symbol_table.getName(getSymbIDByDerivID(var2)) << R"(")"
                   << R"(, "lag2": )" << getLagByDerivID(var2)
                   << R"(, "var3": ")" << symbol_table.getName(getSymbIDByDerivID(var3)) << R"(")"
                   << R"(, "lag3": )" << getLagByDerivID(var3)
                   << R"(, "param": ")" << symbol_table.getName(getSymbIDByDerivID(param)) << R"(")";

      g3p_output << R"(, "val": ")";
      d2->writeJsonOutput(g3p_output, temp_term_union, tef_terms);
      g3p_output << R"("})" << endl;
    }
  g3p_output << "]}" << endl;

  if (writeDetails)
    output << R"("dynamic_model_params_derivative": {)";
  else
    output << R"("dynamic_model_params_derivatives_simple": {)";
  output << model_local_vars_output.str()
         << ", " << model_output.str()
         << ", " << rp_output.str()
         << ", " << gp_output.str()
         << ", " << rpp_output.str()
         << ", " << gpp_output.str()
         << ", " << hp_output.str()
         << ", " << g3p_output.str()
         << "}";
}

void
DynamicModel::substituteVarExpectation(const map<string, expr_t> &subst_table)
{
  for (auto &equation : equations)
    equation = dynamic_cast<BinaryOpNode *>(equation->substituteVarExpectation(subst_table));
}
