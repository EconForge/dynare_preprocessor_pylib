/*
 * Copyright © 2003-2020 Dynare Team
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

    auto convert_block_derivative = [f](const map<tuple<int, int, int>, expr_t> &dt)
                                    {
                                      map<tuple<int, int, int>, expr_t> dt2;
                                      for (const auto &it : dt)
                                        dt2[it.first] = f(it.second);
                                      return dt2;
                                    };
  for (const auto &it : m.blocks_derivatives_other_endo)
    blocks_derivatives_other_endo.emplace_back(convert_block_derivative(it));
  for (const auto &it : m.blocks_derivatives_exo)
    blocks_derivatives_exo.emplace_back(convert_block_derivative(it));
  for (const auto &it : m.blocks_derivatives_exo_det)
    blocks_derivatives_exo_det.emplace_back(convert_block_derivative(it));

  for (const auto &[key, expr] : m.pac_expectation_substitution)
    pac_expectation_substitution.emplace(key, f(expr));
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
  dynJacobianColsNbr{m.dynJacobianColsNbr},
  variableMapping{m.variableMapping},
  blocks_other_endo{m.blocks_other_endo},
  blocks_exo{m.blocks_exo},
  blocks_exo_det{m.blocks_exo_det},
  blocks_jacob_cols_endo{m.blocks_jacob_cols_endo},
  blocks_jacob_cols_other_endo{m.blocks_jacob_cols_other_endo},
  blocks_jacob_cols_exo{m.blocks_jacob_cols_exo},
  blocks_jacob_cols_exo_det{m.blocks_jacob_cols_exo_det},
  var_expectation_functions_to_write{m.var_expectation_functions_to_write},
  pac_mce_alpha_symb_ids{m.pac_mce_alpha_symb_ids},
  pac_h0_indices{m.pac_h0_indices},
  pac_h1_indices{m.pac_h1_indices},
  pac_mce_z1_symb_ids{m.pac_mce_z1_symb_ids},
  pac_eqtag_and_lag{m.pac_eqtag_and_lag},
  pac_model_info{m.pac_model_info},
  pac_equation_info{m.pac_equation_info}
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
  dynJacobianColsNbr = m.dynJacobianColsNbr;
  variableMapping = m.variableMapping;
  blocks_derivatives_other_endo.clear();
  blocks_derivatives_exo.clear();
  blocks_derivatives_exo_det.clear();
  blocks_other_endo = m.blocks_other_endo;
  blocks_exo = m.blocks_exo;
  blocks_exo_det = m.blocks_exo_det;
  blocks_jacob_cols_endo = m.blocks_jacob_cols_endo;
  blocks_jacob_cols_other_endo = m.blocks_jacob_cols_other_endo;
  blocks_jacob_cols_exo = m.blocks_jacob_cols_exo;
  blocks_jacob_cols_exo_det = m.blocks_jacob_cols_exo_det;

  var_expectation_functions_to_write = m.var_expectation_functions_to_write;

  pac_mce_alpha_symb_ids = m.pac_mce_alpha_symb_ids;
  pac_h0_indices = m.pac_h0_indices;
  pac_h1_indices = m.pac_h1_indices;
  pac_mce_z1_symb_ids = m.pac_mce_z1_symb_ids;
  pac_eqtag_and_lag = m.pac_eqtag_and_lag;
  pac_expectation_substitution.clear();
  pac_model_info = m.pac_model_info;
  pac_equation_info = m.pac_equation_info;

  copyHelper(m);

  return *this;
}

void
DynamicModel::compileDerivative(ofstream &code_file, unsigned int &instruction_number, int eq, int symb_id, int lag, const temporary_terms_t &temporary_terms, const temporary_terms_idxs_t &temporary_terms_idxs) const
{
  if (auto it = derivatives[1].find({ eq, getDerivID(symbol_table.getID(SymbolType::endogenous, symb_id), lag) });
      it != derivatives[1].end())
    it->second->compile(code_file, instruction_number, false, temporary_terms, temporary_terms_idxs, true, false);
  else
    {
      FLDZ_ fldz;
      fldz.write(code_file, instruction_number);
    }
}

void
DynamicModel::compileChainRuleDerivative(ofstream &code_file, unsigned int &instruction_number, int blk, int eq, int var, int lag, const temporary_terms_t &temporary_terms, const temporary_terms_idxs_t &temporary_terms_idxs) const
{
  if (auto it = blocks_derivatives[blk].find({ eq, var, lag });
      it != blocks_derivatives[blk].end())
    it->second->compile(code_file, instruction_number, false, temporary_terms, temporary_terms_idxs, true, false);
  else
    {
      FLDZ_ fldz;
      fldz.write(code_file, instruction_number);
    }
}

void
DynamicModel::additionalBlockTemporaryTerms(int blk,
                                            vector<vector<temporary_terms_t>> &blocks_temporary_terms,
                                            map<expr_t, tuple<int, int, int>> &reference_count) const
{
  for (const auto &[ignore, d] : blocks_derivatives_exo[blk])
    d->computeBlockTemporaryTerms(blk, blocks[blk].size, blocks_temporary_terms, reference_count);
  for (const auto &[ignore, d] : blocks_derivatives_exo_det[blk])
    d->computeBlockTemporaryTerms(blk, blocks[blk].size, blocks_temporary_terms, reference_count);
  for (const auto &[ignore, d] : blocks_derivatives_other_endo[blk])
    d->computeBlockTemporaryTerms(blk, blocks[blk].size, blocks_temporary_terms, reference_count);
}

void
DynamicModel::writeDynamicPerBlockHelper(int blk, ostream &output, ExprNodeOutputType output_type, temporary_terms_t &temporary_terms, int nze_stochastic, int nze_deterministic, int nze_exo, int nze_exo_det, int nze_other_endo) const
{
  BlockSimulationType simulation_type = blocks[blk].simulation_type;
  int block_size = blocks[blk].size;
  int block_mfs_size = blocks[blk].mfs_size;
  int block_recursive_size = blocks[blk].getRecursiveSize();

  deriv_node_temp_terms_t tef_terms;

  auto write_eq_tt = [&](int eq)
                     {
                       for (auto it : blocks_temporary_terms[blk][eq])
                         {
                           if (dynamic_cast<AbstractExternalFunctionNode *>(it))
                             it->writeExternalFunctionOutput(output, output_type, temporary_terms, blocks_temporary_terms_idxs, tef_terms);

                           output << "  ";
                           it->writeOutput(output, output_type, blocks_temporary_terms[blk][eq], blocks_temporary_terms_idxs, tef_terms);
                           output << '=';
                           it->writeOutput(output, output_type, temporary_terms, blocks_temporary_terms_idxs, tef_terms);
                           temporary_terms.insert(it);
                           output << ';' << endl;
                         }
                     };

  // The equations
  for (int eq = 0; eq < block_size; eq++)
    {
      write_eq_tt(eq);

      EquationType equ_type = getBlockEquationType(blk, eq);
      BinaryOpNode *e = getBlockEquationExpr(blk, eq);
      expr_t lhs = e->arg1, rhs = e->arg2;
      switch (simulation_type)
        {
        case BlockSimulationType::evaluateBackward:
        case BlockSimulationType::evaluateForward:
          evaluation:
          if (equ_type == EquationType::evaluateRenormalized)
            {
              e = getBlockEquationRenormalizedExpr(blk, eq);
              lhs = e->arg1;
              rhs = e->arg2;
            }
          else if (equ_type != EquationType::evaluate)
            {
              cerr << "Type mismatch for equation " << getBlockEquationID(blk, eq)+1  << endl;
              exit(EXIT_FAILURE);
            }
          output << "  ";
          lhs->writeOutput(output, output_type, temporary_terms, blocks_temporary_terms_idxs);
          output << '=';
          rhs->writeOutput(output, output_type, temporary_terms, blocks_temporary_terms_idxs);
          output << ';' << endl;
          break;
        case BlockSimulationType::solveBackwardSimple:
        case BlockSimulationType::solveForwardSimple:
        case BlockSimulationType::solveBackwardComplete:
        case BlockSimulationType::solveForwardComplete:
        case BlockSimulationType::solveTwoBoundariesComplete:
        case BlockSimulationType::solveTwoBoundariesSimple:
          if (eq < block_recursive_size)
            goto evaluation;
          output << "  residual" << LEFT_ARRAY_SUBSCRIPT(output_type)
                 << eq-block_recursive_size+ARRAY_SUBSCRIPT_OFFSET(output_type)
                 << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=(";
          goto end;
        default:
        end:
          lhs->writeOutput(output, output_type, temporary_terms, blocks_temporary_terms_idxs);
          output << ")-(";
          rhs->writeOutput(output, output_type, temporary_terms, blocks_temporary_terms_idxs);
          output << ");" << endl;
        }
    }

  // The Jacobian if we have to solve the block

  // Write temporary terms for derivatives
  write_eq_tt(blocks[blk].size);

  if (isCOutput(output_type))
    output << "  if (stochastic_mode) {" << endl;
  else
    output << "  if stochastic_mode" << endl;

  ostringstream i_output, j_output, v_output;
  int line_counter = ARRAY_SUBSCRIPT_OFFSET(output_type);
  for (const auto &[indices, d] : blocks_derivatives[blk])
    {
      auto [eq, var, lag] = indices;
      int jacob_col = blocks_jacob_cols_endo[blk].at({ var, lag });
      i_output << "    g1_i" << LEFT_ARRAY_SUBSCRIPT(output_type) << line_counter
               << RIGHT_ARRAY_SUBSCRIPT(output_type) << '=' << eq+1 << ';' << endl;
      j_output << "    g1_j" << LEFT_ARRAY_SUBSCRIPT(output_type) << line_counter
               << RIGHT_ARRAY_SUBSCRIPT(output_type) << '=' << jacob_col+1 << ';' << endl;
      v_output << "    g1_v" << LEFT_ARRAY_SUBSCRIPT(output_type) << line_counter
               << RIGHT_ARRAY_SUBSCRIPT(output_type) << '=';
      d->writeOutput(v_output, output_type, temporary_terms, blocks_temporary_terms_idxs);
      v_output << ';' << endl;
      line_counter++;
    }
  assert(line_counter == nze_stochastic+ARRAY_SUBSCRIPT_OFFSET(output_type));
  output << i_output.str() << j_output.str() << v_output.str();

  i_output.str("");
  j_output.str("");
  v_output.str("");
  line_counter = ARRAY_SUBSCRIPT_OFFSET(output_type);
  for (const auto &[indices, d] : blocks_derivatives_exo[blk])
    {
      auto [eq, var, lag] = indices;
      int jacob_col = blocks_jacob_cols_exo[blk].at({ var, lag });
      i_output << "    g1_x_i" << LEFT_ARRAY_SUBSCRIPT(output_type) << line_counter
               << RIGHT_ARRAY_SUBSCRIPT(output_type) << '=' << eq+1 << ';' << endl;
      j_output << "    g1_x_j" << LEFT_ARRAY_SUBSCRIPT(output_type) << line_counter
               << RIGHT_ARRAY_SUBSCRIPT(output_type) << '=' << jacob_col+1 << ';' << endl;
      v_output << "    g1_x_v" << LEFT_ARRAY_SUBSCRIPT(output_type) << line_counter
               << RIGHT_ARRAY_SUBSCRIPT(output_type) << '=';
      d->writeOutput(v_output, output_type, temporary_terms, blocks_temporary_terms_idxs);
      v_output << ';' << endl;
      line_counter++;
    }
  assert(line_counter == nze_exo+ARRAY_SUBSCRIPT_OFFSET(output_type));
  output << i_output.str() << j_output.str() << v_output.str();

  i_output.str("");
  j_output.str("");
  v_output.str("");
  line_counter = ARRAY_SUBSCRIPT_OFFSET(output_type);
  for (const auto &[indices, d] : blocks_derivatives_exo_det[blk])
    {
      auto [eq, var, lag] = indices;
      int jacob_col = blocks_jacob_cols_exo_det[blk].at({ var, lag });
      i_output << "    g1_xd_i" << LEFT_ARRAY_SUBSCRIPT(output_type) << line_counter
               << RIGHT_ARRAY_SUBSCRIPT(output_type) << '=' << eq+1 << ';' << endl;
      j_output << "    g1_xd_j" << LEFT_ARRAY_SUBSCRIPT(output_type) << line_counter
               << RIGHT_ARRAY_SUBSCRIPT(output_type) << '=' << jacob_col+1 << ';' << endl;
      v_output << "    g1_xd_v" << LEFT_ARRAY_SUBSCRIPT(output_type) << line_counter
               << RIGHT_ARRAY_SUBSCRIPT(output_type) << '=';
      d->writeOutput(v_output, output_type, temporary_terms, blocks_temporary_terms_idxs);
      v_output << ';' << endl;
      line_counter++;
    }
  assert(line_counter == nze_exo_det+ARRAY_SUBSCRIPT_OFFSET(output_type));
  output << i_output.str() << j_output.str() << v_output.str();

  i_output.str("");
  j_output.str("");
  v_output.str("");
  line_counter = ARRAY_SUBSCRIPT_OFFSET(output_type);
  for (const auto &[indices, d] : blocks_derivatives_other_endo[blk])
    {
      auto [eq, var, lag] = indices;
      int jacob_col = blocks_jacob_cols_other_endo[blk].at({ var, lag });
      i_output << "    g1_o_i" << LEFT_ARRAY_SUBSCRIPT(output_type) << line_counter
               << RIGHT_ARRAY_SUBSCRIPT(output_type) << '=' << eq+1 << ';' << endl;
      j_output << "    g1_o_j" << LEFT_ARRAY_SUBSCRIPT(output_type) << line_counter
               << RIGHT_ARRAY_SUBSCRIPT(output_type) << '=' << jacob_col+1 << ';' << endl;
      v_output << "    g1_o_v" << LEFT_ARRAY_SUBSCRIPT(output_type) << line_counter
               << RIGHT_ARRAY_SUBSCRIPT(output_type) << '=';
      d->writeOutput(v_output, output_type, temporary_terms, blocks_temporary_terms_idxs);
      v_output << ';' << endl;
      line_counter++;
    }
  assert(line_counter == nze_other_endo+ARRAY_SUBSCRIPT_OFFSET(output_type));
  output << i_output.str() << j_output.str() << v_output.str();

  // Deterministic mode
  if (simulation_type != BlockSimulationType::evaluateForward
      && simulation_type != BlockSimulationType::evaluateBackward)
    {
      if (isCOutput(output_type))
        output << "  } else {" << endl;
      else
        output << "  else" << endl;
      i_output.str("");
      j_output.str("");
      v_output.str("");
      line_counter = ARRAY_SUBSCRIPT_OFFSET(output_type);
      if (simulation_type == BlockSimulationType::solveBackwardSimple
          || simulation_type == BlockSimulationType::solveForwardSimple
          || simulation_type == BlockSimulationType::solveBackwardComplete
          || simulation_type == BlockSimulationType::solveForwardComplete)
        for (const auto &[indices, d] : blocks_derivatives[blk])
          {
            auto [eq, var, lag] = indices;
            if (lag == 0)
              {
                i_output << "    g1_i" << LEFT_ARRAY_SUBSCRIPT(output_type) << line_counter
                         << RIGHT_ARRAY_SUBSCRIPT(output_type) << '=' << eq+1 << ';' << endl;
                j_output << "    g1_j" << LEFT_ARRAY_SUBSCRIPT(output_type) << line_counter
                         << RIGHT_ARRAY_SUBSCRIPT(output_type) << '='
                         << var+1-block_recursive_size << ';' << endl;
                v_output << "    g1_v" << LEFT_ARRAY_SUBSCRIPT(output_type) << line_counter
                         << RIGHT_ARRAY_SUBSCRIPT(output_type) << '=';
                d->writeOutput(v_output, output_type, temporary_terms, blocks_temporary_terms_idxs);
                v_output << ';' << endl;
                line_counter++;
              }
          }
      else // solveTwoBoundariesSimple || solveTwoBoundariesComplete
        for (const auto &[indices, d] : blocks_derivatives[blk])
        {
          auto [eq, var, lag] = indices;
          assert(lag >= -1 && lag <= 1);
          if (eq >= block_recursive_size && var >= block_recursive_size)
            {
              i_output << "    g1_i" << LEFT_ARRAY_SUBSCRIPT(output_type) << line_counter
                       << RIGHT_ARRAY_SUBSCRIPT(output_type) << '='
                       << eq+1-block_recursive_size << ';' << endl;
              j_output << "    g1_j" << LEFT_ARRAY_SUBSCRIPT(output_type) << line_counter
                       << RIGHT_ARRAY_SUBSCRIPT(output_type) << '='
                       << var+1-block_recursive_size+block_mfs_size*(lag+1) << ';' << endl;
              v_output << "    g1_v" << LEFT_ARRAY_SUBSCRIPT(output_type) << line_counter
                       << RIGHT_ARRAY_SUBSCRIPT(output_type) << '=';
              d->writeOutput(v_output, output_type, temporary_terms, blocks_temporary_terms_idxs);
              v_output << ';' << endl;
              line_counter++;
            }
        }
      assert(line_counter == nze_deterministic+ARRAY_SUBSCRIPT_OFFSET(output_type));
      output << i_output.str() << j_output.str() << v_output.str();
    }
  if (isCOutput(output_type))
    output << "  }" << endl;
  else
    output << "  end" << endl;
}

int
DynamicModel::nzeDeterministicJacobianForBlock(int blk) const
{
  BlockSimulationType simulation_type = blocks[blk].simulation_type;
  int block_recursive_size = blocks[blk].getRecursiveSize();

  int nze_deterministic = 0;
  if (simulation_type == BlockSimulationType::solveTwoBoundariesComplete
      || simulation_type == BlockSimulationType::solveTwoBoundariesSimple)
    nze_deterministic = count_if(blocks_derivatives[blk].begin(), blocks_derivatives[blk].end(),
                                 [=](const auto &kv) {
                                   auto [eq, var, lag] = kv.first;
                                   return eq >= block_recursive_size && var >= block_recursive_size;
                                 });
  else if (simulation_type == BlockSimulationType::solveBackwardSimple
           || simulation_type == BlockSimulationType::solveForwardSimple
           || simulation_type == BlockSimulationType::solveBackwardComplete
           || simulation_type == BlockSimulationType::solveForwardComplete)
    nze_deterministic = count_if(blocks_derivatives[blk].begin(), blocks_derivatives[blk].end(),
                                 [](const auto &kv) {
                                   auto [eq, var, lag] = kv.first;
                                   return lag == 0;
                                 });
  return nze_deterministic;
}

void
DynamicModel::writeDynamicPerBlockMFiles(const string &basename) const
{
  temporary_terms_t temporary_terms; // Temp terms written so far

  for (int blk = 0; blk < static_cast<int>(blocks.size()); blk++)
    {
      BlockSimulationType simulation_type = blocks[blk].simulation_type;
      int block_size = blocks[blk].size;
      int block_mfs_size = blocks[blk].mfs_size;

      // Number of nonzero derivatives for the various Jacobians
      int nze_stochastic = blocks_derivatives[blk].size();
      int nze_deterministic = nzeDeterministicJacobianForBlock(blk);
      int nze_other_endo = blocks_derivatives_other_endo[blk].size();
      int nze_exo = blocks_derivatives_exo[blk].size();
      int nze_exo_det = blocks_derivatives_exo_det[blk].size();

      string filename = packageDir(basename + ".block") + "/dynamic_" + to_string(blk+1) + ".m";
      ofstream output;
      output.open(filename, ios::out | ios::binary);
      if (!output.is_open())
        {
          cerr << "ERROR: Can't open file " << filename << " for writing" << endl;
          exit(EXIT_FAILURE);
        }

      output << "%" << endl
             << "% " << filename << " : Computes dynamic version of one block" << endl
             << "%" << endl
             << "% Warning : this file is generated automatically by Dynare" << endl
             << "%           from model file (.mod)" << endl << endl
             << "%" << endl;
      if (simulation_type == BlockSimulationType::evaluateBackward
          || simulation_type == BlockSimulationType::evaluateForward)
        output << "function [y, T, g1, varargout] = dynamic_" << blk+1 << "(y, x, params, steady_state, T, it_, stochastic_mode)" << endl;
      else
        output << "function [residual, y, T, g1, varargout] = dynamic_" << blk+1 << "(y, x, params, steady_state, T, it_, stochastic_mode)" << endl;

      output << "  % ////////////////////////////////////////////////////////////////////////" << endl
             << "  % //" << string("                     Block ").substr(static_cast<int>(log10(blk + 1))) << blk+1
             << "                                        //" << endl
             << "  % //                     Simulation type "
             << BlockSim(simulation_type) << "  //" << endl
             << "  % ////////////////////////////////////////////////////////////////////////" << endl;

      if (simulation_type != BlockSimulationType::evaluateForward
          && simulation_type != BlockSimulationType::evaluateBackward)
        output << "  residual=zeros(" << block_mfs_size << ",1);" << endl;

      output << "  if stochastic_mode" << endl
             << "    g1_i=zeros(" << nze_stochastic << ",1);" << endl
             << "    g1_j=zeros(" << nze_stochastic << ",1);" << endl
             << "    g1_v=zeros(" << nze_stochastic << ",1);" << endl
             << "    g1_x_i=zeros(" << nze_exo << ",1);" << endl
             << "    g1_x_j=zeros(" << nze_exo << ",1);" << endl
             << "    g1_x_v=zeros(" << nze_exo << ",1);" << endl
             << "    g1_xd_i=zeros(" << nze_exo_det << ",1);" << endl
             << "    g1_xd_j=zeros(" << nze_exo_det << ",1);" << endl
             << "    g1_xd_v=zeros(" << nze_exo_det << ",1);" << endl
             << "    g1_o_i=zeros(" << nze_other_endo << ",1);" << endl
             << "    g1_o_j=zeros(" << nze_other_endo << ",1);" << endl
             << "    g1_o_v=zeros(" << nze_other_endo << ",1);" << endl;
      if (simulation_type != BlockSimulationType::evaluateForward
          && simulation_type != BlockSimulationType::evaluateBackward)
        output << "  else" << endl
               << "    g1_i=zeros(" << nze_deterministic << ",1);" << endl
               << "    g1_j=zeros(" << nze_deterministic << ",1);" << endl
               << "    g1_v=zeros(" << nze_deterministic << ",1);" << endl;
      output << "  end" << endl
             << endl;

      writeDynamicPerBlockHelper(blk, output, ExprNodeOutputType::matlabDynamicModel, temporary_terms,
                                 nze_stochastic, nze_deterministic, nze_exo, nze_exo_det, nze_other_endo);

      output << endl
             << "  if stochastic_mode" << endl
             << "    g1=sparse(g1_i, g1_j, g1_v, " << block_size << ", " << blocks_jacob_cols_endo[blk].size() << ");" << endl
             << "    varargout{1}=sparse(g1_x_i, g1_x_j, g1_x_v, " << block_size << ", " << blocks_jacob_cols_exo[blk].size() << ");" << endl
             << "    varargout{2}=sparse(g1_xd_i, g1_xd_j, g1_xd_v, " << block_size << ", " << blocks_jacob_cols_exo_det[blk].size() << ");" << endl
             << "    varargout{3}=sparse(g1_o_i, g1_o_j, g1_o_v, " << block_size << ", " << blocks_jacob_cols_other_endo[blk].size() << ");" << endl
             << "  else" << endl;
      switch (simulation_type)
        {
        case BlockSimulationType::evaluateForward:
        case BlockSimulationType::evaluateBackward:
          output << "    g1=[];" << endl;
          break;
        case BlockSimulationType::solveBackwardSimple:
        case BlockSimulationType::solveForwardSimple:
        case BlockSimulationType::solveBackwardComplete:
        case BlockSimulationType::solveForwardComplete:
          output << "    g1=sparse(g1_i, g1_j, g1_v, " << block_mfs_size
                 << ", " << block_mfs_size << ");" << endl;
          break;
        case BlockSimulationType::solveTwoBoundariesSimple:
        case BlockSimulationType::solveTwoBoundariesComplete:
          output << "    g1=sparse(g1_i, g1_j, g1_v, " << block_mfs_size
                 << ", " << 3*block_mfs_size << ");" << endl;
          break;
        default:
          break;
        }
      output << "  end" << endl
             << "end" << endl;
      output.close();
    }
}

void
DynamicModel::writeDynamicPerBlockCFiles(const string &basename) const
{
  temporary_terms_t temporary_terms; // Temp terms written so far

  for (int blk = 0; blk < static_cast<int>(blocks.size()); blk++)
    {
      BlockSimulationType simulation_type = blocks[blk].simulation_type;
      int block_size = blocks[blk].size;
      int block_mfs_size = blocks[blk].mfs_size;

      // Number of nonzero derivatives for the various Jacobians
      int nze_stochastic = blocks_derivatives[blk].size();
      int nze_deterministic = nzeDeterministicJacobianForBlock(blk);
      int nze_other_endo = blocks_derivatives_other_endo[blk].size();
      int nze_exo = blocks_derivatives_exo[blk].size();
      int nze_exo_det = blocks_derivatives_exo_det[blk].size();

      string filename = basename + "/model/src/dynamic_" + to_string(blk+1) + ".c";
      ofstream output;
      output.open(filename, ios::out | ios::binary);
      if (!output.is_open())
        {
          cerr << "ERROR: Can't open file " << filename << " for writing" << endl;
          exit(EXIT_FAILURE);
        }

      output << "/* Block " << blk+1 << endl
             << "   " << BlockSim(simulation_type) << " */" << endl
             << endl
             << "#include <math.h>" << endl
             << "#include <stdlib.h>" << endl
             << "#include <stdbool.h>" << endl
             << R"(#include "mex.h")" << endl
             << endl;

      // Write function definition if BinaryOpcode::powerDeriv is used
      writePowerDerivHeader(output);

      output << endl;

      if (simulation_type == BlockSimulationType::evaluateBackward
          || simulation_type == BlockSimulationType::evaluateForward)
        output << "void dynamic_" << blk+1 << "(double *restrict y, const double *restrict x, int nb_row_x, const double *restrict params, const double *restrict steady_state, double *restrict T, int it_, bool stochastic_mode, double *restrict g1_i, double *restrict g1_j, double *restrict g1_v, double *restrict g1_x_i, double *restrict g1_x_j, double *restrict g1_x_v, double *restrict g1_xd_i, double *restrict g1_xd_j, double *restrict g1_xd_v, double *restrict g1_o_i, double *restrict g1_o_j, double *restrict g1_o_v)" << endl;
      else
        output << "void dynamic_" << blk+1 << "(double *restrict y, const double *restrict x, int nb_row_x, const double *restrict params, const double *restrict steady_state, double *restrict T, int it_, bool stochastic_mode, double *restrict residual, double *restrict g1_i, double *restrict g1_j, double *restrict g1_v, double *restrict g1_x_i, double *restrict g1_x_j, double *restrict g1_x_v, double *restrict g1_xd_i, double *restrict g1_xd_j, double *restrict g1_xd_v, double *restrict g1_o_i, double *restrict g1_o_j, double *restrict g1_o_v)" << endl;
      output << '{' << endl;

      writeDynamicPerBlockHelper(blk, output, ExprNodeOutputType::CDynamicModel, temporary_terms,
                                 nze_stochastic, nze_deterministic, nze_exo, nze_exo_det, nze_other_endo);

      output << '}' << endl
             << endl;

      ostringstream header;
      if (simulation_type == BlockSimulationType::evaluateBackward
          || simulation_type == BlockSimulationType::evaluateForward)
        header << "void dynamic_" << blk+1 << "_mx(mxArray *y, const mxArray *x, const mxArray *params, const mxArray *steady_state, mxArray *T, const mxArray *it_, const mxArray *stochastic_mode, mxArray **g1, mxArray **g1_x, mxArray **g1_xd, mxArray **g1_o)";
      else
        header << "void dynamic_" << blk+1 << "_mx(mxArray *y, const mxArray *x, const mxArray *params, const mxArray *steady_state, mxArray *T, const mxArray *it_, const mxArray *stochastic_mode, mxArray **residual, mxArray **g1, mxArray **g1_x, mxArray **g1_xd, mxArray **g1_o)";
      output << header.str() << endl
             << '{' << endl
             << "  int nb_row_x = mxGetM(x);" << endl;

      if (simulation_type != BlockSimulationType::evaluateForward
          && simulation_type != BlockSimulationType::evaluateBackward)
        output << "  *residual = mxCreateDoubleMatrix(" << block_mfs_size << ",1,mxREAL);" << endl;

      output << "  mxArray *g1_i = NULL, *g1_j = NULL, *g1_v = NULL;" << endl
             << "  mxArray *g1_x_i = NULL, *g1_x_j = NULL, *g1_x_v = NULL;" << endl
             << "  mxArray *g1_xd_i = NULL, *g1_xd_j = NULL, *g1_xd_v = NULL;" << endl
             << "  mxArray *g1_o_i = NULL, *g1_o_j = NULL, *g1_o_v = NULL;" << endl
             << "  if (mxGetScalar(stochastic_mode)) {" << endl
             << "    g1_i=mxCreateDoubleMatrix(" << nze_stochastic << ",1,mxREAL);" << endl
             << "    g1_j=mxCreateDoubleMatrix(" << nze_stochastic << ",1,mxREAL);" << endl
             << "    g1_v=mxCreateDoubleMatrix(" << nze_stochastic << ",1,mxREAL);" << endl
             << "    g1_x_i=mxCreateDoubleMatrix(" << nze_exo << ",1,mxREAL);" << endl
             << "    g1_x_j=mxCreateDoubleMatrix(" << nze_exo << ",1,mxREAL);" << endl
             << "    g1_x_v=mxCreateDoubleMatrix(" << nze_exo << ",1,mxREAL);" << endl
             << "    g1_xd_i=mxCreateDoubleMatrix(" << nze_exo_det << ",1,mxREAL);" << endl
             << "    g1_xd_j=mxCreateDoubleMatrix(" << nze_exo_det << ",1,mxREAL);" << endl
             << "    g1_xd_v=mxCreateDoubleMatrix(" << nze_exo_det << ",1,mxREAL);" << endl
             << "    g1_o_i=mxCreateDoubleMatrix(" << nze_other_endo << ",1,mxREAL);" << endl
             << "    g1_o_j=mxCreateDoubleMatrix(" << nze_other_endo << ",1,mxREAL);" << endl
             << "    g1_o_v=mxCreateDoubleMatrix(" << nze_other_endo << ",1,mxREAL);" << endl;
      if (simulation_type != BlockSimulationType::evaluateForward
          && simulation_type != BlockSimulationType::evaluateBackward)
        output << "  } else {" << endl
               << "    g1_i=mxCreateDoubleMatrix(" << nze_deterministic << ",1,mxREAL);" << endl
               << "    g1_j=mxCreateDoubleMatrix(" << nze_deterministic << ",1,mxREAL);" << endl
               << "    g1_v=mxCreateDoubleMatrix(" << nze_deterministic << ",1,mxREAL);" << endl;
      output << "  }" << endl
             << endl;

      // N.B.: In the following, it_ is decreased by 1, to follow C convention
      if (simulation_type == BlockSimulationType::evaluateBackward
          || simulation_type == BlockSimulationType::evaluateForward)
        output << "  dynamic_" << blk+1 << "(mxGetPr(y), mxGetPr(x), nb_row_x, mxGetPr(params), mxGetPr(steady_state), mxGetPr(T), mxGetScalar(it_)-1, mxGetScalar(stochastic_mode), g1_i ? mxGetPr(g1_i) : NULL, g1_j ? mxGetPr(g1_j) : NULL, g1_v ? mxGetPr(g1_v) : NULL, g1_x_i ? mxGetPr(g1_x_i) : NULL, g1_x_j ? mxGetPr(g1_x_j) : NULL, g1_x_v ? mxGetPr(g1_x_v) : NULL, g1_xd_i ? mxGetPr(g1_xd_i) : NULL, g1_xd_j ? mxGetPr(g1_xd_j) : NULL, g1_xd_v ? mxGetPr(g1_xd_v) : NULL, g1_o_i ? mxGetPr(g1_o_i) : NULL, g1_o_j ? mxGetPr(g1_o_j) : NULL, g1_o_v ? mxGetPr(g1_o_v) : NULL);" << endl;
      else
        output << "  dynamic_" << blk+1 << "(mxGetPr(y), mxGetPr(x), nb_row_x, mxGetPr(params), mxGetPr(steady_state), mxGetPr(T), mxGetScalar(it_)-1, mxGetScalar(stochastic_mode), mxGetPr(*residual), g1_i ? mxGetPr(g1_i) : NULL, g1_j ? mxGetPr(g1_j) : NULL, g1_v ? mxGetPr(g1_v) : NULL, g1_x_i ? mxGetPr(g1_x_i) : NULL, g1_x_j ? mxGetPr(g1_x_j) : NULL, g1_x_v ? mxGetPr(g1_x_v) : NULL, g1_xd_i ? mxGetPr(g1_xd_i) : NULL, g1_xd_j ? mxGetPr(g1_xd_j) : NULL, g1_xd_v ? mxGetPr(g1_xd_v) : NULL, g1_o_i ? mxGetPr(g1_o_i) : NULL, g1_o_j ? mxGetPr(g1_o_j) : NULL, g1_o_v ? mxGetPr(g1_o_v) : NULL);" << endl;

      output << endl
             << "  if (mxGetScalar(stochastic_mode)) {" << endl
             << "    mxArray *m = mxCreateDoubleScalar(" << block_size << ");" << endl
             << "    mxArray *n = mxCreateDoubleScalar(" << blocks_jacob_cols_endo[blk].size() << ");" << endl
             << "    mxArray *plhs[1];" << endl
             << "    mxArray *prhs[5] = { g1_i, g1_j, g1_v, m, n };" << endl
             << R"(    mexCallMATLAB(1, plhs, 5, prhs, "sparse");)" << endl
             << "    *g1=plhs[0];" << endl
             << "    mxDestroyArray(g1_i);" << endl
             << "    mxDestroyArray(g1_j);" << endl
             << "    mxDestroyArray(g1_v);" << endl
             << "    mxDestroyArray(n);" << endl
             << "    n = mxCreateDoubleScalar(" << blocks_jacob_cols_exo[blk].size() << ");" << endl
             << "    mxArray *prhs_x[5] = { g1_x_i, g1_x_j, g1_x_v, m, n };" << endl
             << R"(    mexCallMATLAB(1, plhs, 5, prhs_x, "sparse");)" << endl
             << "    *g1_x=plhs[0];" << endl
             << "    mxDestroyArray(g1_x_i);" << endl
             << "    mxDestroyArray(g1_x_j);" << endl
             << "    mxDestroyArray(g1_x_v);" << endl
             << "    mxDestroyArray(n);" << endl
             << "    n = mxCreateDoubleScalar(" << blocks_jacob_cols_exo_det[blk].size() << ");" << endl
             << "    mxArray *prhs_xd[5] = { g1_xd_i, g1_xd_j, g1_xd_v, m, n };" << endl
             << R"(    mexCallMATLAB(1, plhs, 5, prhs_xd, "sparse");)" << endl
             << "    *g1_xd=plhs[0];" << endl
             << "    mxDestroyArray(g1_xd_i);" << endl
             << "    mxDestroyArray(g1_xd_j);" << endl
             << "    mxDestroyArray(g1_xd_v);" << endl
             << "    mxDestroyArray(n);" << endl
             << "    n = mxCreateDoubleScalar(" << blocks_jacob_cols_other_endo[blk].size() << ");" << endl
             << "    mxArray *prhs_o[5] = { g1_o_i, g1_o_j, g1_o_v, m, n };" << endl
             << R"(    mexCallMATLAB(1, plhs, 5, prhs_o, "sparse");)" << endl
             << "    *g1_o=plhs[0];" << endl
             << "    mxDestroyArray(g1_o_i);" << endl
             << "    mxDestroyArray(g1_o_j);" << endl
             << "    mxDestroyArray(g1_o_v);" << endl
             << "    mxDestroyArray(n);" << endl
             << "    mxDestroyArray(m);" << endl
             << "  } else {" << endl;
      switch (simulation_type)
        {
        case BlockSimulationType::evaluateForward:
        case BlockSimulationType::evaluateBackward:
          output << "    *g1=mxCreateDoubleMatrix(0,0,mxREAL);" << endl;
          break;
        case BlockSimulationType::solveBackwardSimple:
        case BlockSimulationType::solveForwardSimple:
        case BlockSimulationType::solveBackwardComplete:
        case BlockSimulationType::solveForwardComplete:
          output << "    mxArray *m = mxCreateDoubleScalar(" << block_mfs_size << ");" << endl
                 << "    mxArray *n = mxCreateDoubleScalar(" << block_mfs_size << ");" << endl
                 << "    mxArray *plhs[1];" << endl
                 << "    mxArray *prhs[5] = { g1_i, g1_j, g1_v, m, n };" << endl
                 << R"(    mexCallMATLAB(1, plhs, 5, prhs, "sparse");)" << endl
                 << "    *g1=plhs[0];" << endl
                 << "    mxDestroyArray(g1_i);" << endl
                 << "    mxDestroyArray(g1_j);" << endl
                 << "    mxDestroyArray(g1_v);" << endl
                 << "    mxDestroyArray(n);" << endl
                 << "    mxDestroyArray(m);" << endl;
          break;
        case BlockSimulationType::solveTwoBoundariesSimple:
        case BlockSimulationType::solveTwoBoundariesComplete:
          output << "    mxArray *m = mxCreateDoubleScalar(" << block_mfs_size << ");" << endl
                 << "    mxArray *n = mxCreateDoubleScalar(" << 3*block_mfs_size << ");" << endl
                 << "    mxArray *plhs[1];" << endl
                 << "    mxArray *prhs[5] = { g1_i, g1_j, g1_v, m, n };" << endl
                 << R"(    mexCallMATLAB(1, plhs, 5, prhs, "sparse");)" << endl
                 << "    *g1=plhs[0];" << endl
                 << "    mxDestroyArray(g1_i);" << endl
                 << "    mxDestroyArray(g1_j);" << endl
                 << "    mxDestroyArray(g1_v);" << endl
                 << "    mxDestroyArray(n);" << endl
                 << "    mxDestroyArray(m);" << endl;
          break;
        default:
          break;
        }
      output << "    *g1_x=mxCreateDoubleMatrix(0,0,mxREAL);" << endl
             << "    *g1_xd=mxCreateDoubleMatrix(0,0,mxREAL);" << endl
             << "    *g1_o=mxCreateDoubleMatrix(0,0,mxREAL);" << endl
             << "  }" << endl
             << "}" << endl;
      output.close();

      filename = basename + "/model/src/dynamic_" + to_string(blk+1) + ".h";
      ofstream header_output;
      header_output.open(filename, ios::out | ios::binary);
      if (!header_output.is_open())
        {
          cerr << "ERROR: Can't open file " << filename << " for writing" << endl;
          exit(EXIT_FAILURE);
        }
      header_output << header.str() << ';' << endl;
      header_output.close();
    }
}

void
DynamicModel::writeDynamicBytecode(const string &basename) const
{
  ostringstream tmp_output;
  ofstream code_file;
  unsigned int instruction_number = 0;
  bool file_open = false;

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
    simulation_type = BlockSimulationType::solveTwoBoundariesComplete;
  else if ((max_endo_lag >= 0) && (max_endo_lead == 0))
    simulation_type = BlockSimulationType::solveForwardComplete;
  else
    simulation_type = BlockSimulationType::solveBackwardComplete;

  writeBytecodeBinFile(basename + "/model/bytecode/dynamic.bin", u_count_int, file_open, simulation_type == BlockSimulationType::solveTwoBoundariesComplete);
  file_open = true;

  //Temporary variables declaration
  FDIMT_ fdimt(temporary_terms_idxs.size());
  fdimt.write(code_file, instruction_number);

  vector<int> exo, exo_det, other_endo;

  for (int i = 0; i < symbol_table.exo_det_nbr(); i++)
    exo_det.push_back(i);
  for (int i = 0; i < symbol_table.exo_nbr(); i++)
    exo.push_back(i);

  map<tuple<int, int, int>, expr_t> first_derivatives_reordered_endo;
  map<tuple<int, SymbolType, int, int>, expr_t> first_derivatives_reordered_exo;
  for (const auto & [indices, d1] : derivatives[1])
    {
      int deriv_id = indices[1];
      int eq = indices[0];
      int symb = getSymbIDByDerivID(deriv_id);
      int var = symbol_table.getTypeSpecificID(symb);
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
                           endo_idx_block2orig,
                           eq_idx_block2orig,
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
                           other_endo);
  fbeginblock.write(code_file, instruction_number);

  temporary_terms_t temporary_terms_union;
  compileTemporaryTerms(code_file, instruction_number, true, false, temporary_terms_union, temporary_terms_idxs);

  compileModelEquations(code_file, instruction_number, true, false, temporary_terms_union, temporary_terms_idxs);

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
          int eq = indices[0];
          int symb = getSymbIDByDerivID(deriv_id);
          int var = symbol_table.getTypeSpecificID(symb);
          int lag = getLagByDerivID(deriv_id);
          FNUMEXPR_ fnumexpr(FirstEndoDerivative, eq, var, lag);
          fnumexpr.write(code_file, instruction_number);
          if (!my_derivatives[eq].size())
            my_derivatives[eq].clear();
          my_derivatives[eq].emplace_back(var, lag, count_u);
          d1->compile(code_file, instruction_number, false, temporary_terms_union, temporary_terms_idxs, true, false);

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
      auto [lag, var, eq] = it.first;
      expr_t d1 = it.second;
      FNUMEXPR_ fnumexpr(FirstEndoDerivative, eq, var, lag);
      fnumexpr.write(code_file, instruction_number);
      if (prev_var != var || prev_lag != lag)
        {
          prev_var = var;
          prev_lag = lag;
          count_col_endo++;
        }
      d1->compile(code_file, instruction_number, false, temporary_terms_union, temporary_terms_idxs, true, false);
      FSTPG3_ fstpg3(eq, var, lag, count_col_endo-1);
      fstpg3.write(code_file, instruction_number);
    }
  prev_var = -1;
  prev_lag = -999999999;
  count_col_exo = 0;
  for (const auto &it : first_derivatives_reordered_exo)
    {
      auto [lag, ignore, var, eq] = it.first;
      expr_t d1 = it.second;
      FNUMEXPR_ fnumexpr(FirstExoDerivative, eq, var, lag);
      fnumexpr.write(code_file, instruction_number);
      if (prev_var != var || prev_lag != lag)
        {
          prev_var = var;
          prev_lag = lag;
          count_col_exo++;
        }
      d1->compile(code_file, instruction_number, false, temporary_terms_union, temporary_terms_idxs, true, false);
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
DynamicModel::writeDynamicBlockBytecode(const string &basename, bool linear_decomposition) const
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
  vector<int> feedback_variables;
  bool file_open = false;
  string main_name;
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

  FDIMT_ fdimt(blocks_temporary_terms_idxs.size());
  fdimt.write(code_file, instruction_number);

  for (int block = 0; block < static_cast<int>(blocks.size()); block++)
    {
      feedback_variables.clear();
      if (block > 0)
        {
          FENDBLOCK_ fendblock;
          fendblock.write(code_file, instruction_number);
        }
      int count_u;
      int u_count_int = 0;
      BlockSimulationType simulation_type = blocks[block].simulation_type;
      int block_size = blocks[block].size;
      int block_mfs = blocks[block].mfs_size;
      int block_recursive = blocks[block].getRecursiveSize();
      int block_max_lag = blocks[block].max_lag;
      int block_max_lead = blocks[block].max_lead;

      if (simulation_type == BlockSimulationType::solveTwoBoundariesSimple
          || simulation_type == BlockSimulationType::solveTwoBoundariesComplete
          || simulation_type == BlockSimulationType::solveBackwardComplete
          || simulation_type == BlockSimulationType::solveForwardComplete)
        {
          writeBlockBytecodeBinFile(basename, block, u_count_int, file_open,
                                    simulation_type == BlockSimulationType::solveTwoBoundariesComplete || simulation_type == BlockSimulationType::solveTwoBoundariesSimple, linear_decomposition);
          file_open = true;
        }

      FBEGINBLOCK_ fbeginblock(block_mfs,
                               simulation_type,
                               blocks[block].first_equation,
                               block_size,
                               endo_idx_block2orig,
                               eq_idx_block2orig,
                               blocks[block].linear,
                               symbol_table.endo_nbr(),
                               block_max_lag,
                               block_max_lead,
                               u_count_int,
                               blocks_jacob_cols_endo[block].size(),
                               blocks_exo_det[block].size(),
                               blocks_jacob_cols_exo_det[block].size(),
                               blocks_exo[block].size(),
                               blocks_jacob_cols_exo[block].size(),
                               blocks_other_endo[block].size(),
                               blocks_jacob_cols_other_endo[block].size(),
                               vector<int>(blocks_exo_det[block].begin(), blocks_exo_det[block].end()),
                               vector<int>(blocks_exo[block].begin(), blocks_exo[block].end()),
                               vector<int>(blocks_other_endo[block].begin(), blocks_other_endo[block].end()));
      fbeginblock.write(code_file, instruction_number);

      temporary_terms_t temporary_terms_union;
      if (linear_decomposition)
        compileTemporaryTerms(code_file, instruction_number, true, false, temporary_terms_union, blocks_temporary_terms_idxs);

      //The Temporary terms
      deriv_node_temp_terms_t tef_terms;

      auto write_eq_tt = [&](int eq)
                         {
                           for (auto it : blocks_temporary_terms[block][eq])
                             {
                               if (dynamic_cast<AbstractExternalFunctionNode *>(it))
                                 it->compileExternalFunctionOutput(code_file, instruction_number, false, temporary_terms_union, blocks_temporary_terms_idxs, true, false, tef_terms);

                               FNUMEXPR_ fnumexpr(TemporaryTerm, static_cast<int>(blocks_temporary_terms_idxs.at(it)));
                               fnumexpr.write(code_file, instruction_number);
                               it->compile(code_file, instruction_number, false, temporary_terms_union, blocks_temporary_terms_idxs, true, false, tef_terms);
                               FSTPT_ fstpt(static_cast<int>(blocks_temporary_terms_idxs.at(it)));
                               fstpt.write(code_file, instruction_number);
                               temporary_terms_union.insert(it);
#ifdef DEBUGC
                               cout << "FSTPT " << v << endl;
                               instruction_number++;
                               code_file.write(&FOK, sizeof(FOK));
                               code_file.write(reinterpret_cast<char *>(&k), sizeof(k));
                               ki++;
#endif
                             }
                         };

      // The equations
      for (i = 0; i < block_size; i++)
        {
          if (!linear_decomposition)
            write_eq_tt(i);

          int variable_ID, equation_ID;
          EquationType equ_type;

          switch (simulation_type)
            {
            evaluation:
            case BlockSimulationType::evaluateBackward:
            case BlockSimulationType::evaluateForward:
              equ_type = getBlockEquationType(block, i);
              {
                FNUMEXPR_ fnumexpr(ModelEquation, getBlockEquationID(block, i));
                fnumexpr.write(code_file, instruction_number);
              }
              if (equ_type == EquationType::evaluate)
                {
                  eq_node = getBlockEquationExpr(block, i);
                  lhs = eq_node->arg1;
                  rhs = eq_node->arg2;
                  rhs->compile(code_file, instruction_number, false, temporary_terms_union, blocks_temporary_terms_idxs, true, false);
                  lhs->compile(code_file, instruction_number, true, temporary_terms_union, blocks_temporary_terms_idxs, true, false);
                }
              else if (equ_type == EquationType::evaluateRenormalized)
                {
                  eq_node = getBlockEquationRenormalizedExpr(block, i);
                  lhs = eq_node->arg1;
                  rhs = eq_node->arg2;
                  rhs->compile(code_file, instruction_number, false, temporary_terms_union, blocks_temporary_terms_idxs, true, false);
                  lhs->compile(code_file, instruction_number, true, temporary_terms_union, blocks_temporary_terms_idxs, true, false);
                }
              break;
            case BlockSimulationType::solveBackwardComplete:
            case BlockSimulationType::solveForwardComplete:
            case BlockSimulationType::solveTwoBoundariesComplete:
            case BlockSimulationType::solveTwoBoundariesSimple:
              if (i < block_recursive)
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
              eq_node = getBlockEquationExpr(block, i);
              lhs = eq_node->arg1;
              rhs = eq_node->arg2;
              lhs->compile(code_file, instruction_number, false, temporary_terms_union, blocks_temporary_terms_idxs, true, false);
              rhs->compile(code_file, instruction_number, false, temporary_terms_union, blocks_temporary_terms_idxs, true, false);

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
      if (simulation_type != BlockSimulationType::evaluateBackward
          && simulation_type != BlockSimulationType::evaluateForward)
        {
          // Write temporary terms for derivatives
          if (!linear_decomposition)
            write_eq_tt(blocks[block].size);

          switch (simulation_type)
            {
            case BlockSimulationType::solveBackwardSimple:
            case BlockSimulationType::solveForwardSimple:
              {
                FNUMEXPR_ fnumexpr(FirstEndoDerivative, getBlockEquationID(block, 0), getBlockVariableID(block, 0), 0);
                fnumexpr.write(code_file, instruction_number);
              }
              compileDerivative(code_file, instruction_number, getBlockEquationID(block, 0), getBlockVariableID(block, 0), 0, temporary_terms_union, blocks_temporary_terms_idxs);
              {
                FSTPG_ fstpg(0);
                fstpg.write(code_file, instruction_number);
              }
              break;

            case BlockSimulationType::solveBackwardComplete:
            case BlockSimulationType::solveForwardComplete:
            case BlockSimulationType::solveTwoBoundariesComplete:
            case BlockSimulationType::solveTwoBoundariesSimple:
              count_u = feedback_variables.size();
              for (const auto &[indices, ignore] : blocks_derivatives[block])
                {
                  auto [eq, var, lag] = indices;
                  int eqr = getBlockEquationID(block, eq);
                  int varr = getBlockVariableID(block, var);
                  if (eq >= block_recursive and var >= block_recursive)
                    {
                      if (lag != 0
                          && (simulation_type == BlockSimulationType::solveForwardComplete
                              || simulation_type == BlockSimulationType::solveBackwardComplete))
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
                      compileChainRuleDerivative(code_file, instruction_number, block, eq, var, lag, temporary_terms_union, blocks_temporary_terms_idxs);
                      FSTPU_ fstpu(count_u);
                      fstpu.write(code_file, instruction_number);
                      count_u++;
                    }
                }
              for (i = 0; i < block_size; i++)
                {
                  if (i >= block_recursive)
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

      for (const auto &[indices, d] : blocks_derivatives[block])
        {
          auto [eq, var, lag] = indices;
          int eqr = getBlockEquationID(block, eq);
          int varr = getBlockVariableID(block, var);
          FNUMEXPR_ fnumexpr(FirstEndoDerivative, eqr, varr, lag);
          fnumexpr.write(code_file, instruction_number);
          compileDerivative(code_file, instruction_number, eqr, varr, lag, temporary_terms_union, blocks_temporary_terms_idxs);
          FSTPG3_ fstpg3(eq, var, lag, blocks_jacob_cols_endo[block].at({ var, lag }));
          fstpg3.write(code_file, instruction_number);
        }
      for (const auto &[indices, d] : blocks_derivatives_exo[block])
        {
          auto [eqr, var, lag] = indices;
          int eq = getBlockEquationID(block, eqr);
          int varr = 0; // Dummy value, actually unused by the bytecode MEX
          FNUMEXPR_ fnumexpr(FirstExoDerivative, eqr, varr, lag);
          fnumexpr.write(code_file, instruction_number);
          d->compile(code_file, instruction_number, false, temporary_terms_union, blocks_temporary_terms_idxs, true, false);
          FSTPG3_ fstpg3(eq, var, lag, blocks_jacob_cols_exo[block].at({ var, lag }));
          fstpg3.write(code_file, instruction_number);
        }
      for (const auto &[indices, d] : blocks_derivatives_exo_det[block])
        {
          auto [eqr, var, lag] = indices;
          int eq = getBlockEquationID(block, eqr);
          int varr = 0; // Dummy value, actually unused by the bytecode MEX
          FNUMEXPR_ fnumexpr(FirstExodetDerivative, eqr, varr, lag);
          fnumexpr.write(code_file, instruction_number);
          d->compile(code_file, instruction_number, false, temporary_terms_union, blocks_temporary_terms_idxs, true, false);
          FSTPG3_ fstpg3(eq, var, lag, blocks_jacob_cols_exo_det[block].at({ var, lag }));
          fstpg3.write(code_file, instruction_number);
        }
      for (const auto &[indices, d] : blocks_derivatives_other_endo[block])
        {
          auto [eqr, var, lag] = indices;
          int eq = getBlockEquationID(block, eqr);
          int varr = 0; // Dummy value, actually unused by the bytecode MEX
          FNUMEXPR_ fnumexpr(FirstOtherEndoDerivative, eqr, varr, lag);
          fnumexpr.write(code_file, instruction_number);
          d->compile(code_file, instruction_number, false, temporary_terms_union, blocks_temporary_terms_idxs, true, false);
          FSTPG3_ fstpg3(eq, var, lag, blocks_jacob_cols_other_endo[block].at({ var, lag }));
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
  string filename = basename + "/model/src/dynamic.c";

  int ntt = temporary_terms_mlv.size() + temporary_terms_derivatives[0].size() + temporary_terms_derivatives[1].size() + temporary_terms_derivatives[2].size() + temporary_terms_derivatives[3].size();

  ofstream output;
  output.open(filename, ios::out | ios::binary);
  if (!output.is_open())
    {
      cerr << "Error: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }
  output << "/*" << endl
         << " * " << filename << " : Computes dynamic model for Dynare" << endl
         << " *" << endl
         << " * Warning : this file is generated automatically by Dynare" << endl
         << " *           from model file (.mod)" << endl
         << " */" << endl
         << endl
         << "#include <math.h>" << endl
         << "#include <stdlib.h>" << endl
         << R"(#include "mex.h")" << endl
         << endl;

  // Write function definition if BinaryOpcode::powerDeriv is used
  writePowerDeriv(output);

  output << endl;

  writeDynamicModel(output, true, false);

  output << "void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])" << endl
         << "{" << endl
         << "  if (nlhs > " << min(computed_derivs_order + 1, 4) << ")" << endl
         << R"(    mexErrMsgTxt("Derivatives of higher order than computed have been requested");)" << endl
         << "  if (nrhs != 5)" << endl
         << R"(    mexErrMsgTxt("Requires exactly 5 input arguments");)" << endl
         << endl
         << "  double *y = mxGetPr(prhs[0]);" << endl
         << "  double *x = mxGetPr(prhs[1]);" << endl
         << "  double *params = mxGetPr(prhs[2]);" << endl
         << "  double *steady_state = mxGetPr(prhs[3]);" << endl
         << "  int it_ = (int) mxGetScalar(prhs[4]) - 1;" << endl
         << "  int nb_row_x = mxGetM(prhs[1]);" << endl
         << endl
         << "  double *T = (double *) malloc(sizeof(double)*" << ntt << ");" << endl
         << endl
         << "  if (nlhs >= 1)" << endl
         << "    {" << endl
         << "       plhs[0] = mxCreateDoubleMatrix(" << equations.size() << ",1, mxREAL);" << endl
         << "       double *residual = mxGetPr(plhs[0]);" << endl
         << "       dynamic_resid_tt(y, x, nb_row_x, params, steady_state, it_, T);" << endl
         << "       dynamic_resid(y, x, nb_row_x, params, steady_state, it_, T, residual);" << endl
         << "    }" << endl
         << endl
         << "  if (nlhs >= 2)" << endl
         << "    {" << endl
         << "       plhs[1] = mxCreateDoubleMatrix(" << equations.size() << ", " << dynJacobianColsNbr << ", mxREAL);" << endl
         << "       double *g1 = mxGetPr(plhs[1]);" << endl
         << "       dynamic_g1_tt(y, x, nb_row_x, params, steady_state, it_, T);" << endl
         << "       dynamic_g1(y, x, nb_row_x, params, steady_state, it_, T, g1);" << endl
         << "    }" << endl
         << endl
         << "  if (nlhs >= 3)" << endl
         << "    {" << endl
         << "      mxArray *g2_i = mxCreateDoubleMatrix(" << NNZDerivatives[2] << ", " << 1 << ", mxREAL);" << endl
         << "      mxArray *g2_j = mxCreateDoubleMatrix(" << NNZDerivatives[2] << ", " << 1 << ", mxREAL);" << endl
         << "      mxArray *g2_v = mxCreateDoubleMatrix(" << NNZDerivatives[2] << ", " << 1 << ", mxREAL);" << endl
         << "      dynamic_g2_tt(y, x, nb_row_x, params, steady_state, it_, T);" << endl
         << "      dynamic_g2(y, x, nb_row_x, params, steady_state, it_, T, mxGetPr(g2_i), mxGetPr(g2_j), mxGetPr(g2_v));" << endl
         << "      mxArray *m = mxCreateDoubleScalar(" << equations.size() << ");" << endl
         << "      mxArray *n = mxCreateDoubleScalar(" << dynJacobianColsNbr*dynJacobianColsNbr << ");" << endl
         << "      mxArray *plhs_sparse[1], *prhs_sparse[5] = { g2_i, g2_j, g2_v, m, n };" << endl
         << R"(      mexCallMATLAB(1, plhs_sparse, 5, prhs_sparse, "sparse");)" << endl
         << "      plhs[2] = plhs_sparse[0];" << endl
         << "      mxDestroyArray(g2_i);" << endl
         << "      mxDestroyArray(g2_j);" << endl
         << "      mxDestroyArray(g2_v);" << endl
         << "      mxDestroyArray(m);" << endl
         << "      mxDestroyArray(n);" << endl
         << "    }" << endl
         << endl
         << "  if (nlhs >= 4)" << endl
         << "    {" << endl
         << "      mxArray *g3_i = mxCreateDoubleMatrix(" << NNZDerivatives[3] << ", " << 1 << ", mxREAL);" << endl
         << "      mxArray *g3_j = mxCreateDoubleMatrix(" << NNZDerivatives[3] << ", " << 1 << ", mxREAL);" << endl
         << "      mxArray *g3_v = mxCreateDoubleMatrix(" << NNZDerivatives[3] << ", " << 1 << ", mxREAL);" << endl
         << "      dynamic_g3_tt(y, x, nb_row_x, params, steady_state, it_, T);" << endl
         << "      dynamic_g3(y, x, nb_row_x, params, steady_state, it_, T, mxGetPr(g3_i), mxGetPr(g3_j), mxGetPr(g3_v));" << endl
         << "      mxArray *m = mxCreateDoubleScalar(" << equations.size() << ");" << endl
         << "      mxArray *n = mxCreateDoubleScalar(" << dynJacobianColsNbr*dynJacobianColsNbr*dynJacobianColsNbr << ");" << endl
         << "      mxArray *plhs_sparse[1], *prhs_sparse[5] = { g3_i, g3_j, g3_v, m, n };" << endl
         << R"(      mexCallMATLAB(1, plhs_sparse, 5, prhs_sparse, "sparse");)" << endl
         << "      plhs[3] = plhs_sparse[0];" << endl
         << "      mxDestroyArray(g3_i);" << endl
         << "      mxDestroyArray(g3_j);" << endl
         << "      mxDestroyArray(g3_v);" << endl
         << "      mxDestroyArray(m);" << endl
         << "      mxDestroyArray(n);" << endl
         << "    }" << endl
         << endl
         << "  free(T);" << endl
         << "}" << endl;

  output.close();
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
DynamicModel::writeBlockBytecodeBinFile(const string &basename, int num, int &u_count_int,
                                        bool &file_open, bool is_two_boundaries, bool linear_decomposition) const
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
  int block_size = blocks[num].size;
  int block_mfs = blocks[num].mfs_size;
  int block_recursive = blocks[num].getRecursiveSize();
  for (const auto &[indices, ignore] : blocks_derivatives[num])
    {
      auto [eq, var, lag] = indices;
      if (lag != 0 && !is_two_boundaries)
        continue;
      if (eq >= block_recursive && var >= block_recursive)
        {
          int v = eq - block_recursive;
          SaveCode.write(reinterpret_cast<char *>(&v), sizeof(v));
          int varr = var - block_recursive + lag * block_mfs;
          SaveCode.write(reinterpret_cast<char *>(&varr), sizeof(varr));
          SaveCode.write(reinterpret_cast<const char *>(&lag), sizeof(lag));
          int u = u_count_int + block_mfs;
          SaveCode.write(reinterpret_cast<char *>(&u), sizeof(u));
          u_count_int++;
        }
    }

  if (is_two_boundaries)
    u_count_int += block_mfs;
  for (j = block_recursive; j < block_size; j++)
    {
      int varr = getBlockVariableID(num, j);
      SaveCode.write(reinterpret_cast<char *>(&varr), sizeof(varr));
    }
  for (j = block_recursive; j < block_size; j++)
    {
      int eqr = getBlockEquationID(num, j);
      SaveCode.write(reinterpret_cast<char *>(&eqr), sizeof(eqr));
    }
  SaveCode.close();
}

void
DynamicModel::writeDynamicBlockMFile(const string &basename) const
{
  ofstream output;
  string filename = packageDir(basename) + "/dynamic.m";
  output.open(filename, ios::out | ios::binary);
  if (!output.is_open())
    {
      cerr << "Error: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }

  output << "function [residual, y, T, g1, varargout] = dynamic(nblock, y, x, params, steady_state, T, it_, stochastic_mode)" << endl
         << "  switch nblock" << endl;

  for (int blk = 0; blk < static_cast<int>(blocks.size()); blk++)
    {
      output << "    case " << blk+1 << endl;

      BlockSimulationType simulation_type = blocks[blk].simulation_type;

      if (simulation_type == BlockSimulationType::evaluateBackward
          || simulation_type == BlockSimulationType::evaluateForward)
        output << "      [y, T, g1, varargout{1:nargout-4}] = " << basename << ".block.dynamic_" << blk+1 << "(y, x, params, steady_state, T, it_, stochastic_mode);" << endl
               << "      residual = [];" << endl;
      else
        output << "      [residual, y, T, g1, varargout{1:nargout-4}] = " << basename << ".block.dynamic_" << blk+1 << "(y, x, params, steady_state, T, it_, stochastic_mode);" << endl;
    }
  output << "  end" << endl
         << "end" << endl;

  output.close();
}

void
DynamicModel::writeDynamicBlockCFile(const string &basename) const
{
  string filename = basename + "/model/src/dynamic.c";

  ofstream output;
  output.open(filename, ios::out | ios::binary);
  if (!output.is_open())
    {
      cerr << "Error: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }

  output << "#include <math.h>" << endl
         << R"(#include "mex.h")" << endl;

  for (int blk = 0; blk < static_cast<int>(blocks.size()); blk++)
    output << R"(#include "dynamic_)" << blk+1 << R"(.h")" << endl;

  output << endl;
  writePowerDeriv(output);

  output << endl
         << "void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])" << endl
         << "{" << endl
         << "  if (nrhs != 8)" << endl
         << R"(    mexErrMsgTxt("Requires exactly 8 input arguments");)" << endl
         << "  if (nlhs > 7)" << endl
         << R"(    mexErrMsgTxt("Accepts at most 7 output arguments");)" << endl
         << "  int nblock = (int) mxGetScalar(prhs[0]);" << endl
         << "  const mxArray *y = prhs[1], *x = prhs[2], *params = prhs[3], *steady_state = prhs[4], *T = prhs[5], *it_ = prhs[6], *stochastic_mode = prhs[7];" << endl
         << "  mxArray *T_new = mxDuplicateArray(T);" << endl
         << "  mxArray *y_new = mxDuplicateArray(y);" << endl
         << "  mxArray *residual, *g1, *g1_x, *g1_xd, *g1_o;" << endl
         << "  switch (nblock)" << endl
         << "    {" << endl;

  for (int blk = 0; blk < static_cast<int>(blocks.size()); blk++)
    {
      output << "    case " << blk+1 << ':' << endl;

      BlockSimulationType simulation_type = blocks[blk].simulation_type;

      if (simulation_type == BlockSimulationType::evaluateBackward
          || simulation_type == BlockSimulationType::evaluateForward)
        output << "      dynamic_" << blk+1 << "_mx(y_new, x, params, steady_state, T_new, it_, stochastic_mode, &g1, &g1_x, &g1_xd, &g1_o);" << endl
               << "      residual = mxCreateDoubleMatrix(0,0,mxREAL);" << endl;
      else
        output << "      dynamic_" << blk+1 << "_mx(y_new, x, params, steady_state, T_new, it_, stochastic_mode, &residual, &g1, &g1_x, &g1_xd, &g1_o);" << endl;
      output << "      break;" << endl;
    }
  output << "    }" << endl
         << endl
         << "  if (nlhs >= 1)" << endl
         << "    plhs[0] = residual;" << endl
         << "  else" << endl
         << "    mxDestroyArray(residual);" << endl
         << "  if (nlhs >= 2)" << endl
         << "    plhs[1] = y_new;" << endl
         << "  else" << endl
         << "    mxDestroyArray(y_new);" << endl
         << "  if (nlhs >= 3)" << endl
         << "    plhs[2] = T_new;" << endl
         << "  else" << endl
         << "    mxDestroyArray(T_new);" << endl
         << "  if (nlhs >= 4)" << endl
         << "    plhs[3] = g1;" << endl
         << "  else" << endl
         << "    mxDestroyArray(g1);" << endl
         << "  if (nlhs >= 5)" << endl
         << "    plhs[4] = g1_x;" << endl
         << "  else" << endl
         << "    mxDestroyArray(g1_x);" << endl
         << "  if (nlhs >= 6)" << endl
         << "    plhs[5] = g1_xd;" << endl
         << "  else" << endl
         << "    mxDestroyArray(g1_xd);" << endl
         << "  if (nlhs >= 7)" << endl
         << "    plhs[6] = g1_o;" << endl
         << "  else" << endl
         << "    mxDestroyArray(g1_o);" << endl
         << "}" << endl;

  output.close();
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
        ostringstream i_output, j_output, v_output;

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
                i_output << "g" << i << "_i" << LEFT_ARRAY_SUBSCRIPT(output_type)
                         << k + ARRAY_SUBSCRIPT_OFFSET(output_type)
                         << RIGHT_ARRAY_SUBSCRIPT(output_type)
                         << "=" << eq + 1 << ";" << endl;
                j_output << "g" << i << "_j" << LEFT_ARRAY_SUBSCRIPT(output_type)
                         << k + ARRAY_SUBSCRIPT_OFFSET(output_type)
                         << RIGHT_ARRAY_SUBSCRIPT(output_type)
                         << "=" << col_idx + 1 << ";" << endl;
                v_output << "g" << i << "_v" << LEFT_ARRAY_SUBSCRIPT(output_type)
                         << k + ARRAY_SUBSCRIPT_OFFSET(output_type)
                         << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=";
                d->writeOutput(v_output, output_type, temp_term_union, temporary_terms_idxs, tef_terms);
                v_output << ";" << endl;

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
                    i_output << "g" << i << "_i" << LEFT_ARRAY_SUBSCRIPT(output_type)
                             << k + ARRAY_SUBSCRIPT_OFFSET(output_type)
                             << RIGHT_ARRAY_SUBSCRIPT(output_type)
                             << "=" << eq + 1 << ";" << endl;
                    j_output << "g" << i << "_j" << LEFT_ARRAY_SUBSCRIPT(output_type)
                             << k + ARRAY_SUBSCRIPT_OFFSET(output_type)
                             << RIGHT_ARRAY_SUBSCRIPT(output_type)
                             << "=" << col_idx_sym + 1 << ";" << endl;
                    v_output << "g" << i << "_v" << LEFT_ARRAY_SUBSCRIPT(output_type)
                             << k + ARRAY_SUBSCRIPT_OFFSET(output_type)
                             << RIGHT_ARRAY_SUBSCRIPT(output_type) << "="
                             << "g" << i << "_v" << LEFT_ARRAY_SUBSCRIPT(output_type)
                             << k-1 + ARRAY_SUBSCRIPT_OFFSET(output_type)
                             << RIGHT_ARRAY_SUBSCRIPT(output_type) << ";" << endl;

                    k++;
                  }
              }
          }
        if (output_type != ExprNodeOutputType::juliaDynamicModel)
          d_output[i] << i_output.str() << j_output.str() << v_output.str();
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
          string gprevname = "g" + to_string(i-1);

          init_output.str("");
          end_output.str("");
          if (derivatives[i].size())
            {
              init_output << gname << "_i = zeros(" << NNZDerivatives[i] << ",1);" << endl
                          << gname << "_j = zeros(" << NNZDerivatives[i] << ",1);" << endl
                          << gname << "_v = zeros(" << NNZDerivatives[i] << ",1);" << endl;
              end_output << gname << " = sparse("
                         << gname << "_i," << gname << "_j," << gname << "_v,"
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
          DynamicOutput << "void dynamic_" << funcname << "_tt(const double *restrict y, const double *restrict x, int nb_row_x, const double *restrict params, const double *restrict steady_state, int it_, double *restrict T)" << endl
                        << "{" << endl
                        << tt_output[i].str()
                        << "}" << endl
                        << endl
                        << "void dynamic_" << funcname << "(const double *restrict y, const double *restrict x, int nb_row_x, const double *restrict params, const double *restrict steady_state, int it_, const double *restrict T, ";
          if (i == 0)
            DynamicOutput << "double *restrict residual";
          else if (i == 1)
            DynamicOutput << "double *restrict g1";
          else
            DynamicOutput << "double *restrict " << funcname << "_i, double *restrict " << funcname << "_j, double *restrict " << funcname << "_v";
          DynamicOutput << ")" << endl
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
DynamicModel::parseIncludeExcludeEquations(const string &inc_exc_eq_tags,
                                           set<pair<string, string>> &eq_tag_set, bool exclude_eqs)
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
                                         equation_tags, false);

  // Ignore output because variables are not excluded when equations marked 'static' are excluded
  ModelTree::includeExcludeEquations(eq_tag_set, exclude_eqs,
                                     static_only_equations, static_only_equations_lineno,
                                     static_only_equations_equation_tags, true);

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
DynamicModel::writeBlockDriverOutput(ostream &output, const string &basename, const string &modstruct,
                                     const vector<int> &state_var, bool estimation_present) const
{
  for (int blk = 0; blk < static_cast<int>(blocks.size()); blk++)
    {
      int block_size = blocks[blk].size;
      output << modstruct << "block_structure.block(" << blk+1 << ").Simulation_Type = " << static_cast<int>(blocks[blk].simulation_type) << ";" << endl
             << modstruct << "block_structure.block(" << blk+1 << ").maximum_lag = " << blocks[blk].max_lag << ";" << endl
             << modstruct << "block_structure.block(" << blk+1 << ").maximum_lead = " << blocks[blk].max_lead << ";" << endl
             << modstruct << "block_structure.block(" << blk+1 << ").maximum_endo_lag = " << blocks[blk].max_endo_lag << ";" << endl
             << modstruct << "block_structure.block(" << blk+1 << ").maximum_endo_lead = " << blocks[blk].max_endo_lead << ";" << endl
             << modstruct << "block_structure.block(" << blk+1 << ").maximum_exo_lag = " << blocks[blk].max_exo_lag << ";" << endl
             << modstruct << "block_structure.block(" << blk+1 << ").maximum_exo_lead = " << blocks[blk].max_exo_lead << ";" << endl
             << modstruct << "block_structure.block(" << blk+1 << ").maximum_exo_det_lag = " << blocks[blk].max_exo_det_lag << ";" << endl
             << modstruct << "block_structure.block(" << blk+1 << ").maximum_exo_det_lead = " << blocks[blk].max_exo_det_lead << ";" << endl
             << modstruct << "block_structure.block(" << blk+1 << ").endo_nbr = " << block_size << ";" << endl
             << modstruct << "block_structure.block(" << blk+1 << ").mfs = " << blocks[blk].mfs_size << ";" << endl
             << modstruct << "block_structure.block(" << blk+1 << ").equation = [";
      for (int eq = 0; eq < block_size; eq++)
        output << " " << getBlockEquationID(blk, eq)+1;
      output << "];" << endl
             << modstruct << "block_structure.block(" << blk+1 << ").variable = [";
      for (int var = 0; var < block_size; var++)
        output << " " << getBlockVariableID(blk, var)+1;
      output << "];" << endl
             << modstruct << "block_structure.block(" << blk+1 << ").exogenous = [";
      for (int exo : blocks_exo[blk])
        output << " " << exo+1;
      output << "];" << endl
             << modstruct << "block_structure.block(" << blk+1 << ").exo_nbr = " << blocks_exo[blk].size() << ";" << endl
             << modstruct << "block_structure.block(" << blk+1 << ").exogenous_det = [";
      for (int exo_det : blocks_exo_det[blk])
        output << " " << exo_det+1;
      output << "];" << endl
             << modstruct << "block_structure.block(" << blk+1 << ").exo_det_nbr = " << blocks_exo_det[blk].size() << ";" << endl
             << modstruct << "block_structure.block(" << blk+1 << ").other_endogenous = [";
      for (int other_endo : blocks_other_endo[blk])
        output << " " << other_endo+1;
      output << "];" << endl
             << modstruct << "block_structure.block(" << blk+1 << ").other_endogenous_block = [";
      for (int other_endo : blocks_other_endo[blk])
        output << " " << endo2block[other_endo]+1;
      output << "];" << endl;

      output << modstruct << "block_structure.block(" << blk+1 << ").tm1 = zeros(" << blocks_other_endo[blk].size() << ", " << state_var.size() << ");" << endl;
      int line = 1;
      for (auto other_endo : blocks_other_endo[blk])
        {
          if (auto it = find(state_var.begin(), state_var.end(), other_endo);
              it != state_var.end())
            output << modstruct << "block_structure.block(" << blk+1 << ").tm1("
                   << line << ", "
                   << distance(state_var.begin(), it)+1 << ") = 1;" << endl;
          line++;
        }

      output << modstruct << "block_structure.block(" << blk+1 << ").other_endo_nbr = " << blocks_other_endo[blk].size() << ";" << endl;

      int count_lead_lag_incidence = 0;
      vector<int> local_state_var;
      output << modstruct << "block_structure.block(" << blk+1 << ").lead_lag_incidence = [" << endl;
      for (int lag = -1; lag <= 1; lag++)
        {
          for (int var = 0; var < block_size; var++)
            {
              for (int eq = 0; eq < block_size; eq++)
                if (blocks_derivatives[blk].find({ eq, var, lag })
                    != blocks_derivatives[blk].end())
                  {
                    if (lag == -1)
                      local_state_var.push_back(getBlockVariableID(blk, var));
                    output << " " << ++count_lead_lag_incidence;
                    goto var_found;
                  }
              output << " 0";
            var_found:
              ;
            }
          output << ";" << endl;
        }
      output << "];" << endl;

      output << modstruct << "block_structure.block(" << blk+1 << ").sorted_col_dr_ghx = [";
      for (int lsv : local_state_var)
        output << distance(state_var.begin(), find(state_var.begin(), state_var.end(), lsv))+1 << " ";
      output << "];" << endl;

      count_lead_lag_incidence = 0;
      output << modstruct << "block_structure.block(" << blk+1 << ").lead_lag_incidence_other = [" << endl;
      for (int lag = -1; lag <= 1; lag++)
        {
          for (int other_endo : blocks_other_endo[blk])
            {
              for (int eq = 0; eq < block_size; eq++)
                if (blocks_derivatives_other_endo[blk].find({ eq, other_endo, lag })
                    != blocks_derivatives_other_endo[blk].end())
                  {
                    output << " " << ++count_lead_lag_incidence;
                    goto other_endo_found;
                  }
              output << " 0";
            other_endo_found:
              ;
            }
          output << ";" << endl;
        }
      output << "];" << endl;

      output << modstruct << "block_structure.block(" << blk+1 << ").n_static = " << blocks[blk].n_static << ";" << endl
             << modstruct << "block_structure.block(" << blk+1 << ").n_forward = " << blocks[blk].n_forward << ";" << endl
             << modstruct << "block_structure.block(" << blk+1 << ").n_backward = " << blocks[blk].n_backward << ";" << endl
             << modstruct << "block_structure.block(" << blk+1 << ").n_mixed = " << blocks[blk].n_mixed << ";" << endl
             << modstruct << "block_structure.block(" << blk+1 << ").is_linear = " << (blocks[blk].linear ? "true" : "false" ) << ';' << endl
             << modstruct << "block_structure.block(" << blk+1 << ").NNZDerivatives = " << blocks_derivatives[blk].size() << ';' << endl;
    }

  output << modstruct << "block_structure.variable_reordered = [";
  for (int i = 0; i < symbol_table.endo_nbr(); i++)
    output << " " << endo_idx_block2orig[i]+1;
  output << "];" << endl
         << modstruct << "block_structure.equation_reordered = [";
  for (int i = 0; i < symbol_table.endo_nbr(); i++)
    output << " " << eq_idx_block2orig[i]+1;
  output << "];" << endl;

  map<int, set<pair<int, int>>> lag_row_incidence;
  for (const auto &[indices, d1] : derivatives[1])
    if (int deriv_id = indices[1];
        getTypeByDerivID(deriv_id) == SymbolType::endogenous)
      {
        int eq = indices[0];
        int var = symbol_table.getTypeSpecificID(getSymbIDByDerivID(deriv_id));
        int lag = getLagByDerivID(deriv_id);
        lag_row_incidence[lag].insert({ eq, var });
      }
  for (auto [lag, eq_var_set] : lag_row_incidence)
    {
      output << modstruct << "block_structure.incidence(" << max_endo_lag+lag+1 << ").lead_lag = " << lag << ";" << endl
             << modstruct << "block_structure.incidence(" << max_endo_lag+lag+1 << ").sparse_IM = [" << endl;
      for (auto [eq, var] : eq_var_set)
        output << " " << eq+1 << " " << var+1 << ";" << endl;
      output << "];" << endl;
    }
  output << modstruct << "block_structure.dyn_tmp_nbr = " << blocks_temporary_terms_idxs.size() << ';' << endl;

  if (estimation_present)
    {
      filesystem::create_directories(basename + "/model/bytecode");
      string main_name = basename + "/model/bytecode/kfi";
      ofstream KF_index_file;
      KF_index_file.open(main_name, ios::out | ios::binary | ios::ate);
      int n_obs = symbol_table.observedVariablesNbr();
      int n_state = state_var.size();
      for (int it : state_var)
        if (symbol_table.isObservedVariable(symbol_table.getID(SymbolType::endogenous, it)))
          n_obs--;

      int n = n_obs + n_state;
      output << modstruct << "nobs_non_statevar = " << n_obs << ";" << endl;
      int nb_diag = 0;

      vector<int> i_nz_state_var(n);
      for (int i = 0; i < n_obs; i++)
        i_nz_state_var[i] = n;
      int lp = n_obs;

      vector<int> state_equ;
      for (int it : state_var)
        state_equ.push_back(eq_idx_block2orig[endo_idx_orig2block[it]]);

      for (int blk = 0; blk < static_cast<int>(blocks.size()); blk++)
        {
          int nze = 0;
          for (int i = 0; i < blocks[blk].size; i++)
            if (int var = getBlockVariableID(blk, i);
                find(state_var.begin(), state_var.end(), var) != state_var.end())
              nze++;

          if (blk == 0)
            {
              set<pair<int, int>> row_state_var_incidence;
              for (const auto &[idx, ignore] : blocks_derivatives[blk])
                if (auto it_state_var = find(state_var.begin(), state_var.end(), getBlockVariableID(blk, get<1>(idx)));
                    it_state_var != state_var.end())
                  if (auto it_state_equ = find(state_equ.begin(), state_equ.end(), getBlockEquationID(blk, get<0>(idx)));
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
              for (auto [equ, var] : row_state_var_incidence)
                col_state_var_incidence.emplace(var, equ);
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
      for (int i = 0; i < lp; i++)
        output << i_nz_state_var[i] << " ";
      output << "];" << endl
             << modstruct << "n_diag = " << nb_diag << ";" << endl;
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
          if ((i < n_obs) || (i >= nb_diag + n_obs) || (j < n_obs) || (j >= nb_diag + n_obs))
            for (int k = n_obs; k < i_nz_state_var[j]; k++)
              {
                int k_n = k * n;
                v_index_KF_2.emplace_back(i * n + j, pair(i + k_n - n_n_obs, j + k_n));
              }
      int size_v_index_KF_2 = v_index_KF_2.size();

      KF_index_file.write(reinterpret_cast<char *>(&size_v_index_KF_2), sizeof(size_v_index_KF_2));
      for (auto &it : v_index_KF_2)
        KF_index_file.write(reinterpret_cast<char *>(&it), sizeof(index_KF));
      KF_index_file.close();
    }
}

void
DynamicModel::writeDriverOutput(ostream &output, const string &basename, bool block_decomposition, bool linear_decomposition, bool use_dll, bool estimation_present, bool compute_xrefs, bool julia) const
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

      /* Write mapping between model local variables and indices in the temporary
         terms vector (dynare#1722) */
      output << modstruct << "model_local_variables_dynamic_tt_idxs = {" << endl;
      for (auto [mlv, value] : temporary_terms_mlv)
        output << "  '" << symbol_table.getName(mlv->symb_id) << "', "
               << temporary_terms_idxs.at(mlv)+1 << ';' << endl;
      output << "};" << endl;
    }

  // Write equation tags
  equation_tags.writeOutput(output, modstruct, julia);

  // Write Occbin tags
  equation_tags.writeOccbinOutput(output, modstruct, julia);

  // Write mapping for variables and equations they are present in
  if (!julia)
    for (const auto &variable : variableMapping)
      {
	output << modstruct << "mapping." << symbol_table.getName(variable.first) << ".eqidx = [";
	for (auto equation : variable.second)
	  output << equation + 1 << " ";
	output << "];" << endl;
      }
  else
    {
      output << modstruct << "mapping.eqidx = Dict(\n";
      for (const auto &variable : variableMapping)
	{
	  output << "        \""
		 << symbol_table.getName(variable.first)
		 << "\" => [";
	  for (auto equation : variable.second)
	    output << equation + 1 << ", ";
	  output << "]," << endl;
	}
      output << ")" << endl;
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

  // Compute list of state variables, ordered in block-order
  vector<int> state_var;
  for (int endoID = 0; endoID < symbol_table.endo_nbr(); endoID++)
    // Loop on negative lags
    for (int lag = -max_endo_lag; lag < 0; lag++)
      try
        {
          getDerivID(symbol_table.getID(SymbolType::endogenous, endo_idx_block2orig[endoID]), lag);
          if (find(state_var.begin(), state_var.end(), endo_idx_block2orig[endoID]) == state_var.end())
            state_var.push_back(endo_idx_block2orig[endoID]);
        }
      catch (UnknownDerivIDException &e)
        {
        }

  // Write the block structure of the model
  if (block_decomposition || linear_decomposition)
    writeBlockDriverOutput(output, basename, modstruct, state_var, estimation_present);

  output << modstruct << "state_var = [";
  for (int it : state_var)
    output << it+1 << (julia ? "," : " ");
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
          set<pair<int, int>> lhs_set, lhs_tmp_set, rhs_set;
          int eqn = equation_tags.getEqnByTag("name", eqtag);
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
          int eqn = equation_tags.getEqnByTag("name", eqtag);
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
          set<pair<int, int>> lhs_set, lhs_tmp_set, rhs_set;
          int eqn = equation_tags.getEqnByTag("name", eqtag);
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

      if (equation->containsPacExpectation(name))
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

          string eqtag = equation_tags.getTagValueByEqnAndKey(&equation - &equations[0], "name");
          if (eqtag.empty())
            {
              cerr << "Every equation with a pac expectation must have been assigned an equation tag name" << endl;
              exit(EXIT_FAILURE);
            }
          if (lhs.first == -1)
            {
              cerr << "walkPacParameters: error obtaining LHS variable." << endl;
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

        string eqtag = equation_tags.getTagValueByEqnAndKey(&equation - &equations[0], "name");
        string eq = eqtag_and_lag[{pac_model_name, eqtag}].first;
        eqtag_and_lag[{pac_model_name, eqtag}] = {eq, equation->PacMaxLag(endogs.begin()->first)};
      }
}

int
DynamicModel::getPacTargetSymbId(const string &pac_model_name) const
{
  for (auto equation : equations)
    if (equation->containsPacExpectation(pac_model_name))
      {
        set<pair<int, int>> lhss;
        equation->arg1->collectDynamicVariables(SymbolType::endogenous, lhss);
        if (lhss.size() != 1)
          throw PacTargetNotIdentifiedException{pac_model_name, "LHS must contain a single endogenous"};
        int lhs_symb_id = lhss.begin()->first;
        if (!symbol_table.isDiffAuxiliaryVariable(lhs_symb_id))
          throw PacTargetNotIdentifiedException{pac_model_name, "LHS must be a diff operator"};
        int undiff_lhs_symb_id = symbol_table.getOrigSymbIdForAuxVar(lhs_symb_id);
        auto barg2 = dynamic_cast<BinaryOpNode *>(equation->arg2);
        if (!barg2)
          throw PacTargetNotIdentifiedException{pac_model_name, "RHS must be a binary operator"};
        auto [optim_share_index, optim_part, non_optim_part, additive_part]
          = barg2->getPacOptimizingShareAndExprNodes(lhs_symb_id, undiff_lhs_symb_id);
        /* If there is an optimization part, restrict the search to that part,
           since it contains the MCE . */
        expr_t mce = optim_part ? optim_part : equation->arg2;

        vector<pair<expr_t, int>> terms;
        mce->decomposeAdditiveTerms(terms);
        for (auto [term, sign] : terms)
          try
            {
              auto [param_id, target_id] = term->matchParamTimesTargetMinusVariable(undiff_lhs_symb_id);
              return target_id;
            }
          catch (ExprNode::MatchFailureException &)
            {
            }
        throw PacTargetNotIdentifiedException{pac_model_name, "No term of the form parameter*(target-LHS_level)"};
      }
  throw PacTargetNotIdentifiedException{pac_model_name, "No equation with the corresponding pac_expectation operator"};
}

void
DynamicModel::declarePacModelConsistentExpectationEndogs(const string &name)
{
  int i = 0;
  for (auto &equation : equations)
    if (equation->containsPacExpectation(name))
      {
        if (!equation_tags.exists(&equation - &equations[0], "name"))
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
  int pac_target_symb_id;
  try
    {
      pac_target_symb_id = getPacTargetSymbId(name);
    }
  catch (PacTargetNotIdentifiedException &e)
    {
      cerr << "Can't identify target for PAC model " << name << ": " << e.message;
      exit(EXIT_FAILURE);
    }
  pac_eqtag_and_lag.insert(eqtag_and_lag.begin(), eqtag_and_lag.end());
  int neqs = 0;
  for (auto &it : eqtag_and_lag)
    {
      assert(it.first.first == name);
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
      auto create_target_lag = [&](int lag)
        {
          if (symbol_table.isAuxiliaryVariable(pac_target_symb_id))
            {
              // We know it is a log, see ExprNode::matchParamTimesTargetMinusVariable()
              /* We don’t use SymbolTable::getOrigSymbIdForAuxVar(), because it
                 does not work for unary ops, and changing this behaviour might
                 break stuff that relies on an exception in this case. */
              auto avi = symbol_table.getAuxVarInfo(pac_target_symb_id);
              return AddLog(AddVariable(avi.get_orig_symb_id(), lag));
            }
          else
            return dynamic_cast<ExprNode *>(AddVariable(pac_target_symb_id, lag));
        };

      expr_t diff_node_to_search = AddDiff(create_target_lag(0));
      if (auto sit = diff_subst_table.find(diff_node_to_search);
          sit != diff_subst_table.end())
        target_base_diff_node = sit->second;
      else
        {
          int symb_id = symbol_table.addDiffAuxiliaryVar(diff_node_to_search->idx, diff_node_to_search);
          target_base_diff_node = AddVariable(symb_id);
          auto neweq = AddEqual(const_cast<VariableNode *>(target_base_diff_node),
                                AddMinus(create_target_lag(0),
                                         create_target_lag(-1)));
          addEquation(neweq, -1);
          addAuxEquation(neweq);
          neqs++;
        }

      map<int, VariableNode *> target_aux_var_to_add;
      const VariableNode *last_aux_var = target_base_diff_node;
      for (int i = 1; i <= pac_max_lag_m - 1; i++, neqs++)
        {
          expr_t this_diff_node = AddDiff(create_target_lag(i));
          int symb_id = symbol_table.addDiffLeadAuxiliaryVar(this_diff_node->idx, this_diff_node,
                                                             last_aux_var->symb_id, last_aux_var->lag);
          VariableNode *current_aux_var = AddVariable(symb_id);
          auto neweq = AddEqual(current_aux_var, AddVariable(last_aux_var->symb_id, 1));
          addEquation(neweq, -1);
          addAuxEquation(neweq);
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
      auto neweq = AddEqual(AddVariable(mce_z1_symb_id),
                            AddMinus(AddTimes(A, AddMinus(const_cast<VariableNode *>(target_base_diff_node), fs)), fp));
      addEquation(neweq, -1);
      neqs++;
      pac_expectation_substitution[{name, eqtag}] = AddVariable(mce_z1_symb_id);
    }
  cout << "Pac Model Consistent Expectation: added " << neqs << " auxiliary variables and equations for model " << name << "." << endl;
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
        if (equation_tags.exists(&equation - &equations[0], "name", it.first.second))
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


  if (linear_decomposition)
    {
      auto first_order_endo_derivatives = collectFirstOrderDerivativesEndogenous();
      equationLinear(first_order_endo_derivatives);

      auto contemporaneous_jacobian = evaluateAndReduceJacobian(eval_context);

      if (!computeNaturalNormalization())
        computeNonSingularNormalization(contemporaneous_jacobian);

      select_non_linear_equations_and_variables();

      equationTypeDetermination(first_order_endo_derivatives, 0);

      reduceBlockDecomposition();

      computeChainRuleJacobian();

      determineLinearBlocks();

      computeBlockDynJacobianCols();

      if (!no_tmp_terms)
        computeBlockTemporaryTerms();
    }
  else if (block)
    {
      auto contemporaneous_jacobian = evaluateAndReduceJacobian(eval_context);

      computeNonSingularNormalization(contemporaneous_jacobian);

      auto [prologue, epilogue] = computePrologueAndEpilogue();

      auto first_order_endo_derivatives = collectFirstOrderDerivativesEndogenous();

      equationTypeDetermination(first_order_endo_derivatives, mfs);

      cout << "Finding the optimal block decomposition of the model ..." << endl;

      computeBlockDecomposition(prologue, epilogue);

      reduceBlockDecomposition();

      printBlockDecomposition();

      computeChainRuleJacobian();

      determineLinearBlocks();

      computeBlockDynJacobianCols();

      if (!no_tmp_terms)
        computeBlockTemporaryTerms();
    }
  else
    {
      computeTemporaryTerms(!use_dll, no_tmp_terms);

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


map<tuple<int, int, int>, DynamicModel::BlockDerivativeType>
DynamicModel::determineBlockDerivativesType(int blk)
{
  map<tuple<int, int, int>, BlockDerivativeType> derivType;
  int size = blocks[blk].size;
  int nb_recursive = blocks[blk].getRecursiveSize();
  for (int lag = -blocks[blk].max_lag; lag <= blocks[blk].max_lead; lag++)
    for (int eq = 0; eq < size; eq++)
      {
        set<pair<int, int>> endos_and_lags;
        int eq_orig = getBlockEquationID(blk, eq);
        equations[eq_orig]->collectEndogenous(endos_and_lags);
        for (int var = 0; var < size; var++)
          if (int var_orig = getBlockVariableID(blk, var);
              endos_and_lags.find({ var_orig, lag }) != endos_and_lags.end())
            {
              if (getBlockEquationType(blk, eq) == EquationType::evaluateRenormalized
                  && eq < nb_recursive)
                /* It’s a normalized recursive equation, we have to recompute
                   the derivative using the chain rule */
                derivType[{ lag, eq, var }] = BlockDerivativeType::normalizedChainRule;
              else if (derivType.find({ lag, eq, var }) == derivType.end())
                derivType[{ lag, eq, var }] = BlockDerivativeType::standard;

              if (var < nb_recursive)
                for (int feedback_var = nb_recursive; feedback_var < size; feedback_var++)
                  if (derivType.find({ lag, var, feedback_var }) != derivType.end())
                    /* A new derivative needs to be computed using the chain rule
                       (a feedback variable appears in the recursive equation
                       defining the current variable) */
                    derivType[{ lag, eq, feedback_var }] = BlockDerivativeType::chainRule;
            }
      }
  return derivType;
}

void
DynamicModel::computeChainRuleJacobian()
{
  int nb_blocks = blocks.size();
  blocks_derivatives.resize(nb_blocks);
  for (int blk = 0; blk < nb_blocks; blk++)
    {
      int nb_recursives = blocks[blk].getRecursiveSize();

      // Create a map from recursive vars to their defining (normalized) equation
      map<int, BinaryOpNode *> recursive_vars;
      for (int i = 0; i < nb_recursives; i++)
        {
          int deriv_id = getDerivID(symbol_table.getID(SymbolType::endogenous, getBlockVariableID(blk, i)), 0);
          if (getBlockEquationType(blk, i) == EquationType::evaluateRenormalized)
            recursive_vars[deriv_id] = getBlockEquationRenormalizedExpr(blk, i);
          else
            recursive_vars[deriv_id] = getBlockEquationExpr(blk, i);
        }

      // Compute the block derivatives
      for (const auto &[indices, derivType] : determineBlockDerivativesType(blk))
        {
          auto [lag, eq, var] = indices;
          int eq_orig = getBlockEquationID(blk, eq), var_orig = getBlockVariableID(blk, var);
          int deriv_id = getDerivID(symbol_table.getID(SymbolType::endogenous, var_orig), lag);
          expr_t d{nullptr};
          switch (derivType)
            {
            case BlockDerivativeType::standard:
              d = derivatives[1][{ eq_orig, deriv_id }];
              break;
            case BlockDerivativeType::chainRule:
              d = equations[eq_orig]->getChainRuleDerivative(deriv_id, recursive_vars);
              break;
            case BlockDerivativeType::normalizedChainRule:
              d = equation_type_and_normalized_equation[eq_orig].second->getChainRuleDerivative(deriv_id, recursive_vars);
              break;
            }

          if (d != Zero)
            blocks_derivatives[blk][{ eq, var, lag }] = d;
        }
    }
}

void
DynamicModel::computeBlockDynJacobianCols()
{
  int nb_blocks = blocks.size();
  blocks_derivatives_other_endo.resize(nb_blocks);
  blocks_derivatives_exo.resize(nb_blocks);
  blocks_derivatives_exo_det.resize(nb_blocks);
  blocks_other_endo.resize(nb_blocks);
  blocks_exo.resize(nb_blocks);
  blocks_exo_det.resize(nb_blocks);
  // Structures used for lexicographic ordering over (lag, var ID)
  vector<set<pair<int, int>>> dynamic_endo(nb_blocks), dynamic_other_endo(nb_blocks),
    dynamic_exo(nb_blocks), dynamic_exo_det(nb_blocks);

  for (auto & [indices, d1] : derivatives[1])
    {
      int eq_orig = indices[0];
      int block_eq = eq2block[eq_orig];
      int eq = getBlockInitialEquationID(block_eq, eq_orig);
      int var = symbol_table.getTypeSpecificID(getSymbIDByDerivID(indices[1]));
      int lag = getLagByDerivID(indices[1]);
      switch (getTypeByDerivID(indices[1]))
        {
        case SymbolType::endogenous:
          if (block_eq == endo2block[var])
            {
              int var_in_block = getBlockInitialVariableID(block_eq, var);
              dynamic_endo[block_eq].emplace(lag, var_in_block);
            }
          else
            {
              blocks_derivatives_other_endo[block_eq][{ eq, var, lag }] = derivatives[1][{ eq_orig, getDerivID(symbol_table.getID(SymbolType::endogenous, var), lag) }];
              blocks_other_endo[block_eq].insert(var);
              dynamic_other_endo[block_eq].emplace(lag, var);
            }
          break;
        case SymbolType::exogenous:
          blocks_derivatives_exo[block_eq][{ eq, var, lag }] = derivatives[1][{ eq_orig, getDerivID(symbol_table.getID(SymbolType::exogenous, var), lag) }];
          blocks_exo[block_eq].insert(var);
          dynamic_exo[block_eq].emplace(lag, var);
          break;
        case SymbolType::exogenousDet:
          blocks_derivatives_exo_det[block_eq][{ eq, var, lag }] = derivatives[1][{ eq_orig, getDerivID(symbol_table.getID(SymbolType::exogenous, var), lag) }];
          blocks_exo_det[block_eq].insert(var);
          dynamic_exo_det[block_eq].emplace(lag, var);
          break;
        default:
          break;
        }
    }

  // Compute Jacobian column indices
  blocks_jacob_cols_endo.resize(nb_blocks);
  blocks_jacob_cols_other_endo.resize(nb_blocks);
  blocks_jacob_cols_exo.resize(nb_blocks);
  blocks_jacob_cols_exo_det.resize(nb_blocks);
  for (int blk = 0; blk < nb_blocks; blk++)
    {
      int index = 0;
      for (auto [lag, var] : dynamic_endo[blk])
        blocks_jacob_cols_endo[blk][{ var, lag }] = index++;

      index = 0;
      for (auto [lag, var] : dynamic_other_endo[blk])
        blocks_jacob_cols_other_endo[blk][{ var, lag }] = index++;

      index = 0;
      for (auto [lag, var] : dynamic_exo[blk])
        blocks_jacob_cols_exo[blk][{ var, lag }] = index++;

      index = 0;
      for (auto [lag, var] : dynamic_exo_det[blk])
        blocks_jacob_cols_exo_det[blk][{ var, lag }] = index++;
    }
}

void
DynamicModel::writeDynamicFile(const string &basename, bool block, bool linear_decomposition, bool bytecode, bool use_dll, const string &mexext, const filesystem::path &matlabroot, const filesystem::path &dynareroot, bool julia) const
{
  filesystem::path model_dir{basename};
  model_dir /= "model";
  if (use_dll)
    filesystem::create_directories(model_dir / "src");
  if (bytecode)
    filesystem::create_directories(model_dir / "bytecode");

  if (linear_decomposition)
    {
      if (bytecode)
        writeDynamicBlockBytecode(basename, linear_decomposition);
      else
        {
          cerr << "'linear_decomposition' option requires the 'bytecode' option" << endl;
          exit(EXIT_FAILURE);
        }
    }
  else if (block)
    {
      if (bytecode)
        writeDynamicBlockBytecode(basename, linear_decomposition);
      else if (use_dll)
        {
          writeDynamicPerBlockCFiles(basename);
          writeDynamicBlockCFile(basename);
          vector<filesystem::path> src_files{blocks.size() + 1};
          for (int blk = 0; blk < static_cast<int>(blocks.size()); blk++)
            src_files[blk] = model_dir / "src" / ("dynamic_" + to_string(blk+1) + ".c");
          src_files[blocks.size()] = model_dir / "src" / "dynamic.c";
          compileMEX(basename, "dynamic", mexext, src_files, matlabroot, dynareroot);
        }
      else if (julia)
        {
          cerr << "'block' option is not available with Julia" << endl;
          exit(EXIT_FAILURE);
        }
      else // M-files
        {
          writeDynamicPerBlockMFiles(basename);
          writeDynamicBlockMFile(basename);
        }
    }
  else
    {
      if (bytecode)
        writeDynamicBytecode(basename);
      else if (use_dll)
        {
          writeDynamicCFile(basename);
          compileMEX(basename, "dynamic", mexext, { model_dir / "src" / "dynamic.c" },
                     matlabroot, dynareroot);
        }
      else if (julia)
        writeDynamicJuliaFile(basename);
      else
        writeDynamicMFile(basename);
    }

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
  for (auto aux_eq : aux_equations)
    if (aux_eq->containsExternalFunction())
      aux_eq->writeExternalFunctionOutput(output, output_type, {}, {}, tef_terms);
  for (auto aux_eq : aux_equations)
    {
      aux_eq->writeOutput(output, output_type, {}, {}, tef_terms);
      output << ";" << endl;
    }
}

void
DynamicModel::clearEquations()
{
  equations.clear();
  equations_lineno.clear();
  equation_tags.clear();
}

void
DynamicModel::replaceMyEquations(DynamicModel &dynamic_model) const
{
  dynamic_model.clearEquations();

  for (size_t i = 0; i < equations.size(); i++)
    dynamic_model.addEquation(equations[i]->clone(dynamic_model), equations_lineno[i]);

  dynamic_model.equation_tags = equation_tags;
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
  map<int, map<string, string>> neweqs_tags;
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
              map<string, string> tags;
              auto tmp = old_equation_tags.getTagsByEqn(i);
              for (const auto &[key, value] : tmp)
                tags[key] = value;
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
  set<int> existing_tags = equation_tags.getEqnsByKey("name");
  for (int eq = 0; eq < static_cast<int>(equations.size()); eq++)
    if (existing_tags.find(eq) == existing_tags.end())
      if (auto lhs_expr = dynamic_cast<VariableNode *>(equations[eq]->arg1);
          lhs_expr
          && !equation_tags.exists("name", symbol_table.getName(lhs_expr->symb_id)))
        equation_tags.add(eq, "name", symbol_table.getName(lhs_expr->symb_id));
      else if (!equation_tags.exists("name", to_string(eq+1)))
        equation_tags.add(eq, "name", to_string(eq+1));
      else
        {
          cerr << "Error creating default equation tag: cannot assign default tag to equation number " << eq+1 << " because it is already in use" << endl;
          exit(EXIT_FAILURE);
        }
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

set<int>
DynamicModel::getEquationNumbersFromTags(const set<string> &eqtags) const
{
  set<int> eqnumbers;
  for (auto &eqtag : eqtags)
    {
      set<int> tmp = equation_tags.getEqnsByTag("name", eqtag);
      if (tmp.empty())
        {
          cerr << "ERROR: looking for equation tag " << eqtag << " failed." << endl;
          exit(EXIT_FAILURE);
        }
      eqnumbers.insert(tmp.begin(), tmp.end());
    }
  return eqnumbers;
}

void
DynamicModel::findPacExpectationEquationNumbers(set<int> &eqnumbers) const
{
  int i = 0;
  for (auto &equation : equations)
    {
      if (equation->containsPacExpectation()
          && find(eqnumbers.begin(), eqnumbers.end(), i) == eqnumbers.end())
        eqnumbers.insert(i);
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
  set<int> eqnumbers = getEquationNumbersFromTags(var_model_eqtags);
  findPacExpectationEquationNumbers(eqnumbers);
  vector<int> eqnumbers_vec(eqnumbers.begin(), eqnumbers.end());
  return substituteUnaryOps(eqnumbers_vec);
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
  for (int mlv : used_local_vars)
    local_variables_table[mlv]->findUnaryOpNodesForAuxVarCreation(nodes);

  // Mark unary ops to be substituted in selected equations
  for (int eqnumber : eqnumbers)
    equations[eqnumber]->findUnaryOpNodesForAuxVarCreation(nodes);

  // Substitute in model local variables
  vector<BinaryOpNode *> neweqs;
  for (int mlv : used_local_vars)
    local_variables_table[mlv] = local_variables_table[mlv]->substituteUnaryOpNodes(nodes, subst_table, neweqs);

  // Substitute in equations
  for (int eq : eqnumbers)
    {
      auto substeq = dynamic_cast<BinaryOpNode *>(equations[eq]->
                                                  substituteUnaryOpNodes(nodes, subst_table, neweqs));
      assert(substeq);
      equations[eq] = substeq;
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
  for (int trendVar : symbol_table.getTrendVarIds())
    eval_context[trendVar] = 2; //not <= 0 bc of log, not 1 bc of powers
}

bool
DynamicModel::isModelLocalVariableUsed() const
{
  set<int> used_local_vars;
  for (size_t i = 0; i < equations.size() && used_local_vars.empty(); i++)
    equations[i]->collectVariables(SymbolType::modelLocalVariable, used_local_vars);
  return !used_local_vars.empty();
}

void
DynamicModel::addStaticOnlyEquation(expr_t eq, int lineno, const map<string, string> &eq_tags)
{
  auto beq = dynamic_cast<BinaryOpNode *>(eq);
  assert(beq && beq->op_code == BinaryOpcode::equal);

  static_only_equations_equation_tags.add(static_only_equations.size(), eq_tags);
  static_only_equations.push_back(beq);
  static_only_equations_lineno.push_back(lineno);
}

size_t
DynamicModel::staticOnlyEquationsNbr() const
{
  return static_only_equations.size();
}

size_t
DynamicModel::dynamicOnlyEquationsNbr() const
{
  return equation_tags.getDynamicEqns().size();
}

bool
DynamicModel::isChecksumMatching(const string &basename, bool block) const
{
  stringstream buffer;

  // Write equation tags
  equation_tags.writeCheckSumInfo(buffer);

  ExprNodeOutputType buffer_type = ExprNodeOutputType::CDynamicModel;

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
  deriv_node_temp_terms_t tef_terms;
  writeJsonModelLocalVariables(output, false, tef_terms);
  output << ", ";
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

      equation_tags.writeJsonAST(output, eq);

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
        if (auto tmp = equation_tags.getTagValueByEqnAndKey(equation, "name"); !tmp.empty())
          output << R"(")" << tmp << (it++ == end_idx_eq ? R"("])" : R"(", )");
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

  writeJsonModelLocalVariables(model_local_vars_output, true, tef_terms);

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
  writeJsonModelLocalVariables(model_local_vars_output, true, tef_terms);

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

void
DynamicModel::checkNoRemainingPacExpectation() const
{
  for (size_t eq = 0; eq < equations.size(); eq++)
    if (equations[eq]->containsPacExpectation())
      {
        cerr << "ERROR: in equation " << equation_tags.getTagValueByEqnAndKey(eq, "name")
             << ", the pac_expectation operator references an unknown pac_model" << endl;
        exit(EXIT_FAILURE);
      }
}
