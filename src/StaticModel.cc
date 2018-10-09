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

#include <boost/filesystem.hpp>

#include "StaticModel.hh"

void
StaticModel::copyHelper(const StaticModel &m)
{
  auto f = [this](expr_t e) { return e->cloneDynamic(*this); };

  auto convert_vector_tt = [this,f](vector<temporary_terms_t> vtt)
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
  for (const auto &it : m.v_temporary_terms_local)
    v_temporary_terms_local.push_back(convert_vector_tt(it));

  for (const auto &it : m.first_chain_rule_derivatives)
    first_chain_rule_derivatives[it.first] = f(it.second);

  for (const auto &it : m.equation_type_and_normalized_equation)
    equation_type_and_normalized_equation.push_back(make_pair(it.first, f(it.second)));

  for (const auto &it : m.blocks_derivatives)
    {
      block_derivatives_equation_variable_laglead_nodeid_t v;
      for (const auto &it2 : it)
        v.push_back(make_pair(it2.first, make_pair(it2.second.first, f(it2.second.second))));
      blocks_derivatives.push_back(v);
    }

  for (const auto &it : m.dynamic_jacobian)
    dynamic_jacobian[it.first] = f(it.second);

  auto convert_derivative_t = [this,f](derivative_t dt)
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

StaticModel::StaticModel(SymbolTable &symbol_table_arg,
                         NumericalConstants &num_constants_arg,
                         ExternalFunctionsTable &external_functions_table_arg) :
  ModelTree{symbol_table_arg, num_constants_arg, external_functions_table_arg}
{
}

StaticModel::StaticModel(const StaticModel &m) :
  ModelTree {m},
  v_temporary_terms_inuse {m.v_temporary_terms_inuse},
  map_idx {m.map_idx},
  map_idx2 {m.map_idx2},
  global_temporary_terms {m.global_temporary_terms},
  block_type_firstequation_size_mfs {m.block_type_firstequation_size_mfs},
  blocks_linear {m.blocks_linear},
  other_endo_block {m.other_endo_block},
  exo_block {m.exo_block},
  exo_det_block {m.exo_det_block},
  block_col_type {m.block_col_type},
  variable_block_lead_lag {m.variable_block_lead_lag},
  equation_block {m.equation_block},
  endo_max_leadlag_block {m.endo_max_leadlag_block},
  other_endo_max_leadlag_block {m.other_endo_max_leadlag_block},
  exo_max_leadlag_block {m.exo_max_leadlag_block},
  exo_det_max_leadlag_block {m.exo_det_max_leadlag_block},
  max_leadlag_block {m.max_leadlag_block}
{
  copyHelper(m);
}

StaticModel &
StaticModel::operator=(const StaticModel &m)
{
  ModelTree::operator=(m);

  v_temporary_terms.clear();
  v_temporary_terms_local.clear();

  v_temporary_terms_inuse = m.v_temporary_terms_inuse;

  first_chain_rule_derivatives.clear();

  map_idx = m.map_idx;
  map_idx2 = m.map_idx2;
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
  block_col_type = m.block_col_type;
  variable_block_lead_lag = m.variable_block_lead_lag;
  equation_block = m.equation_block;
  endo_max_leadlag_block = m.endo_max_leadlag_block;
  other_endo_max_leadlag_block = m.other_endo_max_leadlag_block;
  exo_max_leadlag_block = m.exo_max_leadlag_block;
  exo_det_max_leadlag_block = m.exo_det_max_leadlag_block;
  max_leadlag_block = m.max_leadlag_block;

  copyHelper(m);

  return *this;
}

void
StaticModel::compileDerivative(ofstream &code_file, unsigned int &instruction_number, int eq, int symb_id, map_idx_t &map_idx, temporary_terms_t temporary_terms) const
{
  auto it = first_derivatives.find({ eq, getDerivID(symbol_table.getID(SymbolType::endogenous, symb_id), 0) });
  if (it != first_derivatives.end())
    (it->second)->compile(code_file, instruction_number, false, temporary_terms, map_idx, false, false);
  else
    {
      FLDZ_ fldz;
      fldz.write(code_file, instruction_number);
    }
}

void
StaticModel::compileChainRuleDerivative(ofstream &code_file, unsigned int &instruction_number, int eqr, int varr, int lag, map_idx_t &map_idx, temporary_terms_t temporary_terms) const
{
  auto it = first_chain_rule_derivatives.find({ eqr, { varr, lag } });
  if (it != first_chain_rule_derivatives.end())
    (it->second)->compile(code_file, instruction_number, false, temporary_terms, map_idx, false, false);
  else
    {
      FLDZ_ fldz;
      fldz.write(code_file, instruction_number);
    }
}

void
StaticModel::computeTemporaryTermsOrdered()
{
  map<expr_t, pair<int, int>> first_occurence;
  map<expr_t, int> reference_count;
  BinaryOpNode *eq_node;
  first_derivatives_t::const_iterator it;
  first_chain_rule_derivatives_t::const_iterator it_chr;
  ostringstream tmp_s;
  v_temporary_terms.clear();
  map_idx.clear();

  unsigned int nb_blocks = getNbBlocks();
  v_temporary_terms = vector< vector<temporary_terms_t>>(nb_blocks);
  v_temporary_terms_local = vector< vector<temporary_terms_t>>(nb_blocks);

  v_temporary_terms_inuse = vector<temporary_terms_inuse_t>(nb_blocks);

  map_idx2 = vector<map_idx_t>(nb_blocks);

  temporary_terms.clear();

  //local temporay terms
  for (unsigned int block = 0; block < nb_blocks; block++)
    {
      map<expr_t, int> reference_count_local;
      reference_count_local.clear();
      map<expr_t, pair<int, int>> first_occurence_local;
      first_occurence_local.clear();
      temporary_terms_t temporary_terms_l;
      temporary_terms_l.clear();

      unsigned int block_size = getBlockSize(block);
      unsigned int block_nb_mfs = getBlockMfs(block);
      unsigned int block_nb_recursives = block_size - block_nb_mfs;
      v_temporary_terms_local[block] = vector<temporary_terms_t>(block_size);

      for (unsigned int i = 0; i < block_size; i++)
        {
          if (i < block_nb_recursives && isBlockEquationRenormalized(block, i))
            getBlockEquationRenormalizedExpr(block, i)->computeTemporaryTerms(reference_count_local, temporary_terms_l, first_occurence_local, block, v_temporary_terms_local,  i);
          else
            {
              eq_node = (BinaryOpNode *) getBlockEquationExpr(block, i);
              eq_node->computeTemporaryTerms(reference_count_local, temporary_terms_l, first_occurence_local, block, v_temporary_terms_local,  i);
            }
        }
      for (block_derivatives_equation_variable_laglead_nodeid_t::const_iterator it = blocks_derivatives[block].begin(); it != (blocks_derivatives[block]).end(); it++)
        {
          expr_t id = it->second.second;
          id->computeTemporaryTerms(reference_count_local, temporary_terms_l, first_occurence_local, block, v_temporary_terms_local,  block_size-1);
        }
      set<int> temporary_terms_in_use;
      temporary_terms_in_use.clear();
      v_temporary_terms_inuse[block] = temporary_terms_in_use;
      computeTemporaryTermsMapping(temporary_terms_l, map_idx2[block]);
    }

  // global temporay terms
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
      for (block_derivatives_equation_variable_laglead_nodeid_t::const_iterator it = blocks_derivatives[block].begin(); it != (blocks_derivatives[block]).end(); it++)
        {
          expr_t id = it->second.second;
          id->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, block, v_temporary_terms, block_size-1);
        }
    }

  for (unsigned int block = 0; block < nb_blocks; block++)
    {
      // Collecte the temporary terms reordered
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
      for (block_derivatives_equation_variable_laglead_nodeid_t::const_iterator it = blocks_derivatives[block].begin(); it != (blocks_derivatives[block]).end(); it++)
        {
          expr_t id = it->second.second;
          id->collectTemporary_terms(temporary_terms, temporary_terms_in_use, block);
        }
      for (int i = 0; i < (int) getBlockSize(block); i++)
        for (auto it = v_temporary_terms[block][i].begin();
             it != v_temporary_terms[block][i].end(); it++)
          (*it)->collectTemporary_terms(temporary_terms, temporary_terms_in_use, block);
      v_temporary_terms_inuse[block] = temporary_terms_in_use;
    }
  computeTemporaryTermsMapping(temporary_terms, map_idx);
}

void
StaticModel::computeTemporaryTermsMapping(temporary_terms_t &temporary_terms, map_idx_t &map_idx)
{
  // Add a mapping form node ID to temporary terms order
  int j = 0;
  for (auto temporary_term : temporary_terms)
    map_idx[temporary_term->idx] = j++;
}

void
StaticModel::writeModelEquationsOrdered_M(const string &basename) const
{
  string tmp_s, sps;
  ostringstream tmp_output, tmp1_output, global_output;
  expr_t lhs = nullptr, rhs = nullptr;
  BinaryOpNode *eq_node;
  map<expr_t, int> reference_count;
  temporary_terms_t local_temporary_terms;
  ofstream  output;
  vector<int> feedback_variables;
  deriv_node_temp_terms_t tef_terms;
  ExprNodeOutputType local_output_type;

  local_output_type = ExprNodeOutputType::matlabStaticModelSparse;
  if (global_temporary_terms)
    local_temporary_terms = temporary_terms;

  //----------------------------------------------------------------------
  //For each block
  for (unsigned int block = 0; block < getNbBlocks(); block++)
    {
      //recursive_variables.clear();
      feedback_variables.clear();
      //For a block composed of a single equation determines wether we have to evaluate or to solve the equation
      BlockSimulationType simulation_type = getBlockSimulationType(block);
      unsigned int block_size = getBlockSize(block);
      unsigned int block_mfs = getBlockMfs(block);
      unsigned int block_recursive = block_size - block_mfs;

      tmp1_output.str("");
      tmp1_output << packageDir(basename + ".block") << "/static_" << block+1 << ".m";
      output.open(tmp1_output.str(), ios::out | ios::binary);
      output << "%\n";
      output << "% " << tmp1_output.str() << " : Computes static model for Dynare\n";
      output << "%\n";
      output << "% Warning : this file is generated automatically by Dynare\n";
      output << "%           from model file (.mod)\n\n";
      output << "%/\n";
      if (simulation_type == EVALUATE_BACKWARD || simulation_type == EVALUATE_FORWARD)
        output << "function y = static_" << block+1 << "(y, x, params)\n";
      else
        output << "function [residual, y, g1] = static_" << block+1 << "(y, x, params)\n";

      BlockType block_type;
      if (simulation_type == SOLVE_FORWARD_COMPLETE || simulation_type == SOLVE_BACKWARD_COMPLETE)
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
      output << "  global options_;" << endl;
      //The Temporary terms
      if (simulation_type != EVALUATE_BACKWARD  && simulation_type != EVALUATE_FORWARD)
        output << " g1 = spalloc("  << block_mfs << ", " << block_mfs << ", " << derivative_endo[block].size() << ");" << endl;

      if (v_temporary_terms_inuse[block].size())
        {
          tmp_output.str("");
          for (int it : v_temporary_terms_inuse[block])
            tmp_output << " T" << it;
          output << "  global" << tmp_output.str() << ";\n";
        }

      if (simulation_type != EVALUATE_BACKWARD && simulation_type != EVALUATE_FORWARD)
        output << "  residual=zeros(" << block_mfs << ",1);\n";

      // The equations
      temporary_terms_idxs_t temporary_terms_idxs;
      for (unsigned int i = 0; i < block_size; i++)
        {
          if (!global_temporary_terms)
            local_temporary_terms = v_temporary_terms[block][i];
          temporary_terms_t tt2;
          tt2.clear();
          if (v_temporary_terms[block].size())
            {
              output << "  " << "% //Temporary variables" << endl;
              for (auto it : v_temporary_terms[block][i])
                {
                  if (dynamic_cast<AbstractExternalFunctionNode *>(it) != nullptr)
                    it->writeExternalFunctionOutput(output, local_output_type, tt2, {}, tef_terms);

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
            evaluation:
              output << "  % equation " << getBlockEquationID(block, i)+1 << " variable : " << sModel
                     << " (" << variable_ID+1 << ") " << c_Equation_Type(equ_type) << endl;
              output << "  ";
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
                      output << "\n  ";
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
                  cerr << "Type mismatch for equation " << equation_ID+1  << "\n";
                  exit(EXIT_FAILURE);
                }
              output << ";\n";
              break;
            case SOLVE_BACKWARD_SIMPLE:
            case SOLVE_FORWARD_SIMPLE:
            case SOLVE_BACKWARD_COMPLETE:
            case SOLVE_FORWARD_COMPLETE:
              if (i < block_recursive)
                goto evaluation;
              feedback_variables.push_back(variable_ID);
              output << "  % equation " << equation_ID+1 << " variable : " << sModel
                     << " (" << variable_ID+1 << ") " << c_Equation_Type(equ_type) << endl;
              output << "  " << "residual(" << i+1-block_recursive << ") = (";
              goto end;
            default:
            end:
              output << tmp_output.str();
              output << ") - (";
              rhs->writeOutput(output, local_output_type, local_temporary_terms, {});
              output << ");\n";
            }
        }
      // The Jacobian if we have to solve the block
      if (simulation_type == SOLVE_BACKWARD_SIMPLE   || simulation_type == SOLVE_FORWARD_SIMPLE
          || simulation_type == SOLVE_BACKWARD_COMPLETE || simulation_type == SOLVE_FORWARD_COMPLETE)
        output << "  " << sps << "% Jacobian  " << endl;
      switch (simulation_type)
        {
        case SOLVE_BACKWARD_SIMPLE:
        case SOLVE_FORWARD_SIMPLE:
        case SOLVE_BACKWARD_COMPLETE:
        case SOLVE_FORWARD_COMPLETE:
          for (auto it = blocks_derivatives[block].begin(); it != (blocks_derivatives[block]).end(); it++)
            {
              unsigned int eq = it->first.first;
              unsigned int var = it->first.second;
              unsigned int eqr = getBlockEquationID(block, eq);
              unsigned int varr = getBlockVariableID(block, var);
              expr_t id = it->second.second;
              output << "    g1(" << eq+1-block_recursive << ", " << var+1-block_recursive << ") = ";
              id->writeOutput(output, local_output_type, local_temporary_terms, {});
              output << "; % variable=" << symbol_table.getName(symbol_table.getID(SymbolType::endogenous, varr))
                     << "(" << 0
                     << ") " << varr+1
                     << ", equation=" << eqr+1 << endl;
            }
          break;
        default:
          break;
        }
      output << "end" << endl;
      output.close();
    }
}

void
StaticModel::writeModelEquationsCode(const string &basename, map_idx_t map_idx) const
{

  ostringstream tmp_output;
  ofstream code_file;
  unsigned int instruction_number = 0;
  bool file_open = false;

  boost::filesystem::create_directories(basename + "/model/bytecode");

  string main_name = basename + "/model/bytecode/static.cod";
  code_file.open(main_name, ios::out | ios::binary | ios::ate);
  if (!code_file.is_open())
    {
      cerr << "Error : Can't open file \"" << main_name << "\" for writing" << endl;
      exit(EXIT_FAILURE);
    }
  int count_u;
  int u_count_int = 0;

  Write_Inf_To_Bin_File(basename + "/model/bytecode/static.bin", u_count_int, file_open, false, symbol_table.endo_nbr());
  file_open = true;

  //Temporary variables declaration
  FDIMST_ fdimst(temporary_terms.size());
  fdimst.write(code_file, instruction_number);
  FBEGINBLOCK_ fbeginblock(symbol_table.endo_nbr(),
                           SOLVE_FORWARD_COMPLETE,
                           0,
                           symbol_table.endo_nbr(),
                           variable_reordered,
                           equation_reordered,
                           false,
                           symbol_table.endo_nbr(),
                           0,
                           0,
                           u_count_int,
                           symbol_table.endo_nbr()
                           );
  fbeginblock.write(code_file, instruction_number);

  // Add a mapping form node ID to temporary terms order
  int j = 0;
  for (auto temporary_term : temporary_terms)
    map_idx[temporary_term->idx] = j++;
  compileTemporaryTerms(code_file, instruction_number, temporary_terms, map_idx, false, false);

  compileModelEquations(code_file, instruction_number, temporary_terms, map_idx, false, false);

  FENDEQU_ fendequ;
  fendequ.write(code_file, instruction_number);

  // Get the current code_file position and jump if eval = true
  streampos pos1 = code_file.tellp();
  FJMPIFEVAL_ fjmp_if_eval(0);
  fjmp_if_eval.write(code_file, instruction_number);
  int prev_instruction_number = instruction_number;

  vector<vector<pair<int, int>>> derivatives;
  derivatives.resize(symbol_table.endo_nbr());
  count_u = symbol_table.endo_nbr();
  for (const auto & first_derivative : first_derivatives)
    {
      int deriv_id = first_derivative.first.second;
      if (getTypeByDerivID(deriv_id) == SymbolType::endogenous)
        {
          expr_t d1 = first_derivative.second;
          unsigned int eq = first_derivative.first.first;
          int symb = getSymbIDByDerivID(deriv_id);
          unsigned int var = symbol_table.getTypeSpecificID(symb);
          FNUMEXPR_ fnumexpr(FirstEndoDerivative, eq, var);
          fnumexpr.write(code_file, instruction_number);
          if (!derivatives[eq].size())
            derivatives[eq].clear();
          derivatives[eq].emplace_back(var, count_u);

          d1->compile(code_file, instruction_number, false, temporary_terms, map_idx, false, false);

          FSTPSU_ fstpsu(count_u);
          fstpsu.write(code_file, instruction_number);
          count_u++;
        }
    }
  for (int i = 0; i < symbol_table.endo_nbr(); i++)
    {
      FLDR_ fldr(i);
      fldr.write(code_file, instruction_number);
      if (derivatives[i].size())
        {
          for (vector<pair<int, int>>::const_iterator it = derivatives[i].begin();
               it != derivatives[i].end(); it++)
            {
              FLDSU_ fldsu(it->second);
              fldsu.write(code_file, instruction_number);
              FLDSV_ fldsv{static_cast<int>(SymbolType::endogenous), static_cast<unsigned int>(it->first)};
              fldsv.write(code_file, instruction_number);
              FBINARY_ fbinary{static_cast<int>(BinaryOpcode::times)};
              fbinary.write(code_file, instruction_number);
              if (it != derivatives[i].begin())
                {
                  FBINARY_ fbinary{static_cast<int>(BinaryOpcode::plus)};
                  fbinary.write(code_file, instruction_number);
                }
            }
          FBINARY_ fbinary{static_cast<int>(BinaryOpcode::minus)};
          fbinary.write(code_file, instruction_number);
        }
      FSTPSU_ fstpsu(i);
      fstpsu.write(code_file, instruction_number);
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

  temporary_terms_t tt2;
  tt2.clear();
  temporary_terms_t tt3;
  tt3.clear();

  // The Jacobian if we have to solve the block determinsitic bloc
  for (const auto & first_derivative : first_derivatives)
    {
      int deriv_id = first_derivative.first.second;
      if (getTypeByDerivID(deriv_id) == SymbolType::endogenous)
        {
          expr_t d1 = first_derivative.second;
          unsigned int eq = first_derivative.first.first;
          int symb = getSymbIDByDerivID(deriv_id);
          unsigned int var = symbol_table.getTypeSpecificID(symb);
          FNUMEXPR_ fnumexpr(FirstEndoDerivative, eq, var);
          fnumexpr.write(code_file, instruction_number);
          if (!derivatives[eq].size())
            derivatives[eq].clear();
          derivatives[eq].emplace_back(var, count_u);

          d1->compile(code_file, instruction_number, false, temporary_terms, map_idx, false, false);
          FSTPG2_ fstpg2(eq, var);
          fstpg2.write(code_file, instruction_number);
        }
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
StaticModel::writeModelEquationsCode_Block(const string &basename, map_idx_t map_idx, vector<map_idx_t> map_idx2) const
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
  deriv_node_temp_terms_t tef_terms;
  bool file_open = false;

  boost::filesystem::create_directories(basename + "/model/bytecode");

  string main_name = basename + "/model/bytecode/static.cod";
  code_file.open(main_name, ios::out | ios::binary | ios::ate);
  if (!code_file.is_open())
    {
      cerr << "Error : Can't open file \"" << main_name << "\" for writing" << endl;
      exit(EXIT_FAILURE);
    }
  //Temporary variables declaration

  FDIMST_ fdimst(temporary_terms.size());
  fdimst.write(code_file, instruction_number);

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

      if (simulation_type == SOLVE_TWO_BOUNDARIES_SIMPLE || simulation_type == SOLVE_TWO_BOUNDARIES_COMPLETE
          || simulation_type == SOLVE_BACKWARD_COMPLETE || simulation_type == SOLVE_FORWARD_COMPLETE)
        {
          Write_Inf_To_Bin_File_Block(basename, block, u_count_int, file_open);
          file_open = true;
        }

      FBEGINBLOCK_ fbeginblock(block_mfs,
                               simulation_type,
                               getBlockFirstEquation(block),
                               block_size,
                               variable_reordered,
                               equation_reordered,
                               blocks_linear[block],
                               symbol_table.endo_nbr(),
                               0,
                               0,
                               u_count_int,
                               /*symbol_table.endo_nbr()*/ block_size
                               );

      fbeginblock.write(code_file, instruction_number);

      // Get the current code_file position and jump if eval = true
      streampos pos1 = code_file.tellp();
      FJMPIFEVAL_ fjmp_if_eval(0);
      fjmp_if_eval.write(code_file, instruction_number);
      int prev_instruction_number = instruction_number;

      for (i = 0; i < (int) block_size; i++)
        {
          //The Temporary terms
          temporary_terms_t tt2;
          tt2.clear();
          if (v_temporary_terms[block].size())
            {
              for (auto it : v_temporary_terms[block][i])
                {
                  if (dynamic_cast<AbstractExternalFunctionNode *>(it) != nullptr)
                    it->compileExternalFunctionOutput(code_file, instruction_number, false, tt2, map_idx, false, false, tef_terms);

                  FNUMEXPR_ fnumexpr(TemporaryTerm, (int)(map_idx.find(it->idx)->second));
                  fnumexpr.write(code_file, instruction_number);
                  it->compile(code_file, instruction_number, false, tt2, map_idx, false, false, tef_terms);
                  FSTPST_ fstpst((int)(map_idx.find(it->idx)->second));
                  fstpst.write(code_file, instruction_number);
                  // Insert current node into tt2
                  tt2.insert(it);
                }
            }

          // The equations
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
                  rhs->compile(code_file, instruction_number, false, temporary_terms, map_idx, false, false);
                  lhs->compile(code_file, instruction_number, true, temporary_terms, map_idx, false, false);
                }
              else if (equ_type == E_EVALUATE_S)
                {
                  eq_node = (BinaryOpNode *) getBlockEquationRenormalizedExpr(block, i);
                  lhs = eq_node->get_arg1();
                  rhs = eq_node->get_arg2();
                  rhs->compile(code_file, instruction_number, false, temporary_terms, map_idx, false, false);
                  lhs->compile(code_file, instruction_number, true, temporary_terms, map_idx, false, false);
                }
              break;
            case SOLVE_BACKWARD_COMPLETE:
            case SOLVE_FORWARD_COMPLETE:
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
              lhs->compile(code_file, instruction_number, false, temporary_terms, map_idx, false, false);
              rhs->compile(code_file, instruction_number, false, temporary_terms, map_idx, false, false);

              FBINARY_ fbinary{static_cast<int>(BinaryOpcode::minus)};
              fbinary.write(code_file, instruction_number);

              FSTPR_ fstpr(i - block_recursive);
              fstpr.write(code_file, instruction_number);
            }
        }
      FENDEQU_ fendequ;
      fendequ.write(code_file, instruction_number);

      // The Jacobian if we have to solve the block
      if    (simulation_type != EVALUATE_BACKWARD
             && simulation_type != EVALUATE_FORWARD)
        {
          switch (simulation_type)
            {
            case SOLVE_BACKWARD_SIMPLE:
            case SOLVE_FORWARD_SIMPLE:
              {
                FNUMEXPR_ fnumexpr(FirstEndoDerivative, 0, 0);
                fnumexpr.write(code_file, instruction_number);
              }
              compileDerivative(code_file, instruction_number, getBlockEquationID(block, 0), getBlockVariableID(block, 0), map_idx, temporary_terms);
              {
                FSTPG_ fstpg(0);
                fstpg.write(code_file, instruction_number);
              }
              break;

            case SOLVE_BACKWARD_COMPLETE:
            case SOLVE_FORWARD_COMPLETE:
              count_u = feedback_variables.size();
              for (auto it = blocks_derivatives[block].begin(); it != (blocks_derivatives[block]).end(); it++)
                {
                  unsigned int eq = it->first.first;
                  unsigned int var = it->first.second;
                  unsigned int eqr = getBlockEquationID(block, eq);
                  unsigned int varr = getBlockVariableID(block, var);
                  if (eq >= block_recursive && var >= block_recursive)
                    {
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
                      FNUMEXPR_ fnumexpr(FirstEndoDerivative, eqr, varr);
                      fnumexpr.write(code_file, instruction_number);
                      compileChainRuleDerivative(code_file, instruction_number, eqr, varr, 0, map_idx, temporary_terms);
                      FSTPSU_ fstpsu(count_u);
                      fstpsu.write(code_file, instruction_number);
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
                          FLDSU_ fldsu(Uf[v].Ufl->u);
                          fldsu.write(code_file, instruction_number);
                          FLDSV_ fldsv{static_cast<int>(SymbolType::endogenous), static_cast<unsigned int>(Uf[v].Ufl->var)};
                          fldsv.write(code_file, instruction_number);

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

                      FSTPSU_ fstpsu(i - block_recursive);
                      fstpsu.write(code_file, instruction_number);

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

      temporary_terms_t tt2;
      tt2.clear();
      temporary_terms_t tt3;
      tt3.clear();
      deriv_node_temp_terms_t tef_terms2;

      for (i = 0; i < (int) block_size; i++)
        {
          if (v_temporary_terms_local[block].size())
            {
              for (auto it = v_temporary_terms_local[block][i].begin();
                   it != v_temporary_terms_local[block][i].end(); it++)
                {
                  if (dynamic_cast<AbstractExternalFunctionNode *>(*it) != nullptr)
                    (*it)->compileExternalFunctionOutput(code_file, instruction_number, false, tt3, map_idx2[block], false, false, tef_terms2);

                  FNUMEXPR_ fnumexpr(TemporaryTerm, (int)(map_idx2[block].find((*it)->idx)->second));
                  fnumexpr.write(code_file, instruction_number);

                  (*it)->compile(code_file, instruction_number, false, tt3, map_idx2[block], false, false, tef_terms);

                  FSTPST_ fstpst((int)(map_idx2[block].find((*it)->idx)->second));
                  fstpst.write(code_file, instruction_number);
                  // Insert current node into tt2
                  tt3.insert(*it);
                  tt2.insert(*it);
                }
            }

          // The equations
          int variable_ID, equation_ID;
          EquationType equ_type;
          switch (simulation_type)
            {
            evaluation_l:
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
                  rhs->compile(code_file, instruction_number, false, tt2, map_idx2[block], false, false);
                  lhs->compile(code_file, instruction_number, true, tt2, map_idx2[block], false, false);
                }
              else if (equ_type == E_EVALUATE_S)
                {
                  eq_node = (BinaryOpNode *) getBlockEquationRenormalizedExpr(block, i);
                  lhs = eq_node->get_arg1();
                  rhs = eq_node->get_arg2();
                  rhs->compile(code_file, instruction_number, false, tt2, map_idx2[block], false, false);
                  lhs->compile(code_file, instruction_number, true, tt2, map_idx2[block], false, false);
                }
              break;
            case SOLVE_BACKWARD_COMPLETE:
            case SOLVE_FORWARD_COMPLETE:
              if (i < (int) block_recursive)
                goto evaluation_l;
              variable_ID = getBlockVariableID(block, i);
              equation_ID = getBlockEquationID(block, i);
              feedback_variables.push_back(variable_ID);
              Uf[equation_ID].Ufl = nullptr;
              goto end_l;
            default:
            end_l:
              FNUMEXPR_ fnumexpr(ModelEquation, getBlockEquationID(block, i));
              fnumexpr.write(code_file, instruction_number);
              eq_node = (BinaryOpNode *) getBlockEquationExpr(block, i);
              lhs = eq_node->get_arg1();
              rhs = eq_node->get_arg2();
              lhs->compile(code_file, instruction_number, false, tt2, map_idx2[block], false, false);
              rhs->compile(code_file, instruction_number, false, tt2, map_idx2[block], false, false);

              FBINARY_ fbinary{static_cast<int>(BinaryOpcode::minus)};
              fbinary.write(code_file, instruction_number);

              FSTPR_ fstpr(i - block_recursive);
              fstpr.write(code_file, instruction_number);
            }
        }
      FENDEQU_ fendequ_l;
      fendequ_l.write(code_file, instruction_number);

      // The Jacobian if we have to solve the block determinsitic bloc
      switch (simulation_type)
        {
        case SOLVE_BACKWARD_SIMPLE:
        case SOLVE_FORWARD_SIMPLE:
          {
            FNUMEXPR_ fnumexpr(FirstEndoDerivative, 0, 0);
            fnumexpr.write(code_file, instruction_number);
          }
          compileDerivative(code_file, instruction_number, getBlockEquationID(block, 0), getBlockVariableID(block, 0), map_idx2[block], tt2 /*temporary_terms*/);
          {
            FSTPG2_ fstpg2(0, 0);
            fstpg2.write(code_file, instruction_number);
          }
          break;
        case EVALUATE_BACKWARD:
        case EVALUATE_FORWARD:
        case SOLVE_BACKWARD_COMPLETE:
        case SOLVE_FORWARD_COMPLETE:
          count_u = feedback_variables.size();
          for (auto it = blocks_derivatives[block].begin(); it != (blocks_derivatives[block]).end(); it++)
            {
              unsigned int eq = it->first.first;
              unsigned int var = it->first.second;
              unsigned int eqr = getBlockEquationID(block, eq);
              unsigned int varr = getBlockVariableID(block, var);
              FNUMEXPR_ fnumexpr(FirstEndoDerivative, eqr, varr, 0);
              fnumexpr.write(code_file, instruction_number);

              compileChainRuleDerivative(code_file, instruction_number, eqr, varr, 0, map_idx2[block], tt2 /*temporary_terms*/);

              FSTPG2_ fstpg2(eq, var);
              fstpg2.write(code_file, instruction_number);
            }
          break;
        default:
          break;
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
StaticModel::Write_Inf_To_Bin_File_Block(const string &basename, const int &num,
                                         int &u_count_int, bool &file_open) const
{
  int j;
  std::ofstream SaveCode;
  string filename = basename + "/model/bytecode/static.bin";
  if (file_open)
    SaveCode.open(filename, ios::out | ios::in | ios::binary | ios::ate);
  else
    SaveCode.open(filename, ios::out | ios::binary);
  if (!SaveCode.is_open())
    {
      cerr << "Error : Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }
  u_count_int = 0;
  unsigned int block_size = getBlockSize(num);
  unsigned int block_mfs = getBlockMfs(num);
  unsigned int block_recursive = block_size - block_mfs;
  for (auto it = blocks_derivatives[num].begin(); it != (blocks_derivatives[num]).end(); it++)
    {
      unsigned int eq = it->first.first;
      unsigned int var = it->first.second;
      int lag = 0;
      if (eq >= block_recursive && var >= block_recursive)
        {
          int v = eq - block_recursive;
          SaveCode.write(reinterpret_cast<char *>(&v), sizeof(v));
          int varr = var - block_recursive;
          SaveCode.write(reinterpret_cast<char *>(&varr), sizeof(varr));
          SaveCode.write(reinterpret_cast<char *>(&lag), sizeof(lag));
          int u = u_count_int + block_mfs;
          SaveCode.write(reinterpret_cast<char *>(&u), sizeof(u));
          u_count_int++;
        }
    }

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

map<pair<int, pair<int, int >>, expr_t>
StaticModel::collect_first_order_derivatives_endogenous()
{
  map<pair<int, pair<int, int >>, expr_t> endo_derivatives;
  for (auto & first_derivative : first_derivatives)
    {
      if (getTypeByDerivID(first_derivative.first.second) == SymbolType::endogenous)
        {
          int eq = first_derivative.first.first;
          int var = symbol_table.getTypeSpecificID(getSymbIDByDerivID(first_derivative.first.second));
          int lag = 0;
          endo_derivatives[{ eq, { var, lag } }] = first_derivative.second;
        }
    }
  return endo_derivatives;
}

void
StaticModel::computingPass(const eval_context_t &eval_context, bool no_tmp_terms, bool hessian, bool thirdDerivatives, int paramsDerivsOrder, bool block, bool bytecode, const bool nopreprocessoroutput)
{
  initializeVariablesAndEquations();

  vector<BinaryOpNode *> neweqs;
  for (unsigned int eq = 0; eq < equations.size() - aux_equations.size(); eq++)
    {
      expr_t eq_tmp = equations[eq]->substituteStaticAuxiliaryVariable();
      neweqs.push_back(dynamic_cast<BinaryOpNode *>(eq_tmp->toStatic(*this)));
    }

  for (auto & aux_equation : aux_equations)
    {
      expr_t eq_tmp = aux_equation->substituteStaticAuxiliaryDefinition();
      neweqs.push_back(dynamic_cast<BinaryOpNode *>(eq_tmp->toStatic(*this)));
    }

  equations.clear();
  copy(neweqs.begin(), neweqs.end(), back_inserter(equations));
  // Compute derivatives w.r. to all endogenous, and possibly exogenous and exogenous deterministic
  set<int> vars;

  for (int i = 0; i < symbol_table.endo_nbr(); i++)
    {
      int id = symbol_table.getID(SymbolType::endogenous, i);
      //      if (!symbol_table.isAuxiliaryVariableButNotMultiplier(id))
      vars.insert(getDerivID(id, 0));
    }

  // Launch computations
  if (!nopreprocessoroutput)
    cout << "Computing static model derivatives:" << endl
         << " - order 1" << endl;
  first_derivatives.clear();

  computeJacobian(vars);

  if (hessian)
    {
      if (!nopreprocessoroutput)
        cout << " - order 2" << endl;
      computeHessian(vars);
    }

  if (thirdDerivatives)
    {
      if (!nopreprocessoroutput)
        cout << " - order 3" << endl;
      computeThirdDerivatives(vars);
    }

  if (paramsDerivsOrder > 0)
    {
      if (!nopreprocessoroutput)
        cout << " - derivatives of Jacobian/Hessian w.r. to parameters" << endl;
      computeParamsDerivatives(paramsDerivsOrder);

      if (!no_tmp_terms)
        computeParamsDerivativesTemporaryTerms();
    }

  if (block)
    {
      jacob_map_t contemporaneous_jacobian, static_jacobian;
      vector<unsigned int> n_static, n_forward, n_backward, n_mixed;

      // for each block contains pair<Size, Feddback_variable>
      vector<pair<int, int>> blocks;

      evaluateAndReduceJacobian(eval_context, contemporaneous_jacobian, static_jacobian, dynamic_jacobian, cutoff, false);

      computeNonSingularNormalization(contemporaneous_jacobian, cutoff, static_jacobian, dynamic_jacobian);

      computePrologueAndEpilogue(static_jacobian, equation_reordered, variable_reordered);

      map<pair<int, pair<int, int>>, expr_t> first_order_endo_derivatives = collect_first_order_derivatives_endogenous();

      equation_type_and_normalized_equation = equationTypeDetermination(first_order_endo_derivatives, variable_reordered, equation_reordered, mfs);

      if (!nopreprocessoroutput)
        cout << "Finding the optimal block decomposition of the model ...\n";

      lag_lead_vector_t equation_lag_lead, variable_lag_lead;

      computeBlockDecompositionAndFeedbackVariablesForEachBlock(static_jacobian, dynamic_jacobian, equation_reordered, variable_reordered, blocks, equation_type_and_normalized_equation, false, false, mfs, inv_equation_reordered, inv_variable_reordered, equation_lag_lead, variable_lag_lead, n_static, n_forward, n_backward, n_mixed);

      block_type_firstequation_size_mfs = reduceBlocksAndTypeDetermination(dynamic_jacobian, blocks, equation_type_and_normalized_equation, variable_reordered, equation_reordered, n_static, n_forward, n_backward, n_mixed, block_col_type, false);

      printBlockDecomposition(blocks);

      computeChainRuleJacobian(blocks_derivatives);

      blocks_linear = BlockLinear(blocks_derivatives, variable_reordered);

      collect_block_first_order_derivatives();

      global_temporary_terms = true;
      if (!no_tmp_terms)
        computeTemporaryTermsOrdered();
    }
  else
    {
      computeTemporaryTerms(true, no_tmp_terms);
      if (bytecode && !no_tmp_terms)
        computeTemporaryTermsMapping(temporary_terms, map_idx);
    }
}

void
StaticModel::writeStaticMFile(const string &basename) const
{
  writeStaticModel(basename, false, false);
}

void
StaticModel::writeWrapperFunctions(const string &basename, const string &ending) const
{
  string name;
  if (ending == "g1")
    name = "static_resid_g1";
  else if (ending == "g2")
    name = "static_resid_g1_g2";
  else if (ending == "g3")
    name = "static_resid_g1_g2_g3";

  string filename = packageDir(basename) + "/" + name + ".m";
  ofstream output;
  output.open(filename, ios::out | ios::binary);
  if (!output.is_open())
    {
      cerr << "Error: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }

  if (ending == "g1")
    output << "function [residual, g1] = " << name << "(T, y, x, params, T_flag)" << endl
           << "% function [residual, g1] = " << name << "(T, y, x, params, T_flag)" << endl;
  else if (ending == "g2")
    output << "function [residual, g1, g2] = " << name << "(T, y, x, params, T_flag)" << endl
           << "% function [residual, g1, g2] = " << name << "(T, y, x, params, T_flag)" << endl;
  else if (ending == "g3")
    output << "function [residual, g1, g2, g3] = " << name << "(T, y, x, params, T_flag)" << endl
           << "% function [residual, g1, g2, g3] = " << name << "(T, y, x, params, T_flag)" << endl;

  output << "%" << endl
         << "% Wrapper function automatically created by Dynare" << endl
         << "%" << endl
         << endl
         << "    if T_flag" << endl
         << "        T = " << basename << ".static_" << ending << "_tt(T, y, x, params);" << endl
         << "    end" << endl;

  if (ending == "g1")
    output << "    residual = " << basename << ".static_resid(T, y, x, params, false);" << endl
           << "    g1       = " << basename << ".static_g1(T, y, x, params, false);" << endl;
  else if (ending == "g2")
    output << "    [residual, g1] = " << basename << ".static_resid_g1(T, y, x, params, false);" << endl
           << "    g2       = " << basename << ".static_g2(T, y, x, params, false);" << endl;
  else if (ending == "g3")
    output << "    [residual, g1, g2] = " << basename << ".static_resid_g1_g2(T, y, x, params, false);" << endl
           << "    g3       = " << basename << ".static_g3(T, y, x, params, false);" << endl;

  output << endl << "end" << endl;
  output.close();
}

void
StaticModel::writeStaticModelHelper(const string &basename,
                                    const string &name, const string &retvalname,
                                    const string &name_tt, size_t ttlen,
                                    const string &previous_tt_name,
                                    const ostringstream &init_s, const ostringstream &end_s,
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

  output << "function T = " << name_tt << "(T, y, x, params)" << endl
         << "% function T = " << name_tt << "(T, y, x, params)" << endl
         << "%" << endl
         << "% File created by Dynare Preprocessor from .mod file" << endl
         << "%" << endl
         << "% Inputs:" << endl
         << "%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function" << endl
         << "%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order" << endl
         << "%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order" << endl
         << "%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order" << endl
         << "%" << endl
         << "% Output:" << endl
         << "%   T         [#temp variables by 1]  double   vector of temporary terms" << endl
         << "%" << endl << endl
         << "assert(length(T) >= " << ttlen << ");" << endl
         << endl;

  if (!previous_tt_name.empty())
    output << "T = " << basename << "." << previous_tt_name << "(T, y, x, params);" << endl << endl;

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

  output << "function " << retvalname << " = " << name << "(T, y, x, params, T_flag)" << endl
         << "% function " << retvalname << " = " << name << "(T, y, x, params, T_flag)" << endl
         << "%" << endl
         << "% File created by Dynare Preprocessor from .mod file" << endl
         << "%" << endl
         << "% Inputs:" << endl
         << "%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function" << endl
         << "%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order" << endl
         << "%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order" << endl
         << "%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order" << endl
         << "%                                              to evaluate the model" << endl
         << "%   T_flag    boolean                 boolean  flag saying whether or not to calculate temporary terms" << endl
         << "%" << endl
         << "% Output:" << endl
         << "%   " << retvalname << endl
         << "%" << endl << endl;

  if (!name_tt.empty())
    output << "if T_flag" << endl
           << "    T = " << basename << "."  << name_tt << "(T, y, x, params);" << endl
           << "end" << endl;

  output << init_s.str() << endl
         << s.str()
         << end_s.str() << endl
         << "end" << endl;
  output.close();
}

void
StaticModel::writeStaticMatlabCompatLayer(const string &basename) const
{
  string filename = packageDir(basename) + "/static.m";
  ofstream output;
  output.open(filename, ios::out | ios::binary);
  if (!output.is_open())
    {
      cerr << "Error: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }
  int ntt = temporary_terms_mlv.size() + temporary_terms_res.size() + temporary_terms_g1.size() + temporary_terms_g2.size() + temporary_terms_g3.size();

  output << "function [residual, g1, g2, g3] = static(y, x, params)" << endl
         << "    T = NaN(" << ntt << ", 1);" << endl
         << "    if nargout <= 1" << endl
         << "        residual = " << basename << ".static_resid(T, y, x, params, true);" << endl
         << "    elseif nargout == 2" << endl
         << "        [residual, g1] = " << basename << ".static_resid_g1(T, y, x, params, true);" << endl
         << "    elseif nargout == 3" << endl
         << "        [residual, g1, g2] = " << basename << ".static_resid_g1_g2(T, y, x, params, true);" << endl
         << "    else" << endl
         << "        [residual, g1, g2, g3] = " << basename << ".static_resid_g1_g2_g3(T, y, x, params, true);" << endl
         << "    end" << endl
         << "end" << endl;

  output.close();
}

void
StaticModel::writeStaticModel(ostream &StaticOutput, bool use_dll, bool julia) const
{
  writeStaticModel("", StaticOutput, use_dll, julia);
}

void
StaticModel::writeStaticModel(const string &basename, bool use_dll, bool julia) const
{
  ofstream StaticOutput;
  writeStaticModel(basename, StaticOutput, use_dll, julia);
}

void
StaticModel::writeStaticModel(const string &basename,
                              ostream &StaticOutput, bool use_dll, bool julia) const
{
  ostringstream model_tt_output;             // Used for storing model temp vars
  ostringstream model_output;                // Used for storing model equations
  ostringstream jacobian_tt_output;          // Used for storing jacobian temp vars
  ostringstream jacobian_output;             // Used for storing jacobian equations
  ostringstream hessian_tt_output;           // Used for storing Hessian temp vars
  ostringstream hessian_output;              // Used for storing Hessian equations
  ostringstream third_derivatives_tt_output; // Used for storing third order derivatives temp terms
  ostringstream third_derivatives_output;    // Used for storing third order derivatives equations
  ostringstream for_sym;

  ExprNodeOutputType output_type = (use_dll ? ExprNodeOutputType::CStaticModel :
                                    julia ? ExprNodeOutputType::juliaStaticModel : ExprNodeOutputType::matlabStaticModel);

  deriv_node_temp_terms_t tef_terms;
  temporary_terms_t temp_term_union;

  for (auto it : temporary_terms_mlv)
    temp_term_union.insert(it.first);
  writeModelLocalVariableTemporaryTerms(temp_term_union, temporary_terms_mlv,
                                        model_tt_output, output_type, tef_terms);

  writeTemporaryTerms(temporary_terms_res,
                      temp_term_union,
                      temporary_terms_idxs,
                      model_tt_output, output_type, tef_terms);
  temp_term_union.insert(temporary_terms_res.begin(), temporary_terms_res.end());

  writeModelEquations(model_output, output_type, temp_term_union);

  int nrows = equations.size();
  int JacobianColsNbr = symbol_table.endo_nbr();
  int hessianColsNbr = JacobianColsNbr*JacobianColsNbr;

  // Write Jacobian w.r. to endogenous only
  if (!first_derivatives.empty())
    {
      writeTemporaryTerms(temporary_terms_g1,
                          temp_term_union,
                          temporary_terms_idxs,
                          jacobian_tt_output, output_type, tef_terms);
      temp_term_union.insert(temporary_terms_g1.begin(), temporary_terms_g1.end());
    }
  for (const auto & first_derivative : first_derivatives)
    {
      int eq, var;
      tie(eq, var) = first_derivative.first;
      expr_t d1 = first_derivative.second;
      int symb_id = getSymbIDByDerivID(var);

      jacobianHelper(jacobian_output, eq, symbol_table.getTypeSpecificID(symb_id), output_type);
      jacobian_output << "=";
      d1->writeOutput(jacobian_output, output_type,
                      temp_term_union, temporary_terms_idxs, tef_terms);
      jacobian_output << ";" << endl;
    }

  int g2ncols = symbol_table.endo_nbr() * symbol_table.endo_nbr();
  // Write Hessian w.r. to endogenous only (only if 2nd order derivatives have been computed)
  if (!second_derivatives.empty())
    {
      writeTemporaryTerms(temporary_terms_g2,
                          temp_term_union,
                          temporary_terms_idxs,
                          hessian_tt_output, output_type, tef_terms);
      temp_term_union.insert(temporary_terms_g2.begin(), temporary_terms_g2.end());

      int k = 0; // Keep the line of a 2nd derivative in v2
      for (const auto & second_derivative : second_derivatives)
        {
          int eq, var1, var2;
          tie(eq, var1, var2) = second_derivative.first;
          expr_t d2 = second_derivative.second;

          int symb_id1 = getSymbIDByDerivID(var1);
          int symb_id2 = getSymbIDByDerivID(var2);

          int tsid1 = symbol_table.getTypeSpecificID(symb_id1);
          int tsid2 = symbol_table.getTypeSpecificID(symb_id2);

          int col_nb = tsid1*symbol_table.endo_nbr()+tsid2;
          int col_nb_sym = tsid2*symbol_table.endo_nbr()+tsid1;

          if (output_type == ExprNodeOutputType::juliaStaticModel)
            {
              for_sym << "g2[" << eq + 1 << "," << col_nb + 1 << "]";
              hessian_output << "  @inbounds " << for_sym.str() << " = ";
              d2->writeOutput(hessian_output, output_type, temp_term_union, temporary_terms_idxs, tef_terms);
              hessian_output << endl;
            }
          else
            {
              sparseHelper(2, hessian_output, k, 0, output_type);
              hessian_output << "=" << eq + 1 << ";" << endl;

              sparseHelper(2, hessian_output, k, 1, output_type);
              hessian_output << "=" << col_nb + 1 << ";" << endl;

              sparseHelper(2, hessian_output, k, 2, output_type);
              hessian_output << "=";
              d2->writeOutput(hessian_output, output_type, temp_term_union, temporary_terms_idxs, tef_terms);
              hessian_output << ";" << endl;

              k++;
            }

          // Treating symetric elements
          if (symb_id1 != symb_id2)
            if (output_type == ExprNodeOutputType::juliaStaticModel)
              hessian_output << "  @inbounds g2[" << eq + 1 << "," << col_nb_sym + 1 << "] = "
                             << for_sym.str() << endl;
            else
              {
                sparseHelper(2, hessian_output, k, 0, output_type);
                hessian_output << "=" << eq + 1 << ";" << endl;

                sparseHelper(2, hessian_output, k, 1, output_type);
                hessian_output << "=" << col_nb_sym + 1 << ";" << endl;

                sparseHelper(2, hessian_output, k, 2, output_type);
                hessian_output << "=";
                sparseHelper(2, hessian_output, k-1, 2, output_type);
                hessian_output << ";" << endl;

                k++;
              }
        }
    }

  // Writing third derivatives
  if (!third_derivatives.empty())
    {
      writeTemporaryTerms(temporary_terms_g3,
                          temp_term_union,
                          temporary_terms_idxs,
                          third_derivatives_tt_output, output_type, tef_terms);
      temp_term_union.insert(temporary_terms_g3.begin(), temporary_terms_g3.end());

      int k = 0; // Keep the line of a 3rd derivative in v3
      for (const auto & third_derivative : third_derivatives)
        {
          int eq, var1, var2, var3;
          tie(eq, var1, var2, var3) = third_derivative.first;
          expr_t d3 = third_derivative.second;

          int id1 = getSymbIDByDerivID(var1);
          int id2 = getSymbIDByDerivID(var2);
          int id3 = getSymbIDByDerivID(var3);

          // Reference column number for the g3 matrix
          int ref_col = id1 * hessianColsNbr + id2 * JacobianColsNbr + id3;

          if (output_type == ExprNodeOutputType::juliaStaticModel)
            {
              for_sym << "g3[" << eq + 1 << "," << ref_col + 1 << "]";
              third_derivatives_output << "  @inbounds " << for_sym.str() << " = ";
              d3->writeOutput(third_derivatives_output, output_type, temp_term_union, temporary_terms_idxs, tef_terms);
              third_derivatives_output << endl;
            }
          else
            {
              sparseHelper(3, third_derivatives_output, k, 0, output_type);
              third_derivatives_output << "=" << eq + 1 << ";" << endl;

              sparseHelper(3, third_derivatives_output, k, 1, output_type);
              third_derivatives_output << "=" << ref_col + 1 << ";" << endl;

              sparseHelper(3, third_derivatives_output, k, 2, output_type);
              third_derivatives_output << "=";
              d3->writeOutput(third_derivatives_output, output_type, temp_term_union, temporary_terms_idxs, tef_terms);
              third_derivatives_output << ";" << endl;
            }

          // Compute the column numbers for the 5 other permutations of (id1,id2,id3)
          // and store them in a set (to avoid duplicates if two indexes are equal)
          set<int> cols;
          cols.insert(id1 * hessianColsNbr + id3 * JacobianColsNbr + id2);
          cols.insert(id2 * hessianColsNbr + id1 * JacobianColsNbr + id3);
          cols.insert(id2 * hessianColsNbr + id3 * JacobianColsNbr + id1);
          cols.insert(id3 * hessianColsNbr + id1 * JacobianColsNbr + id2);
          cols.insert(id3 * hessianColsNbr + id2 * JacobianColsNbr + id1);

          int k2 = 1; // Keeps the offset of the permutation relative to k
          for (int col : cols)
            if (col != ref_col)
              if (output_type == ExprNodeOutputType::juliaStaticModel)
                third_derivatives_output << "  @inbounds g3[" << eq + 1 << "," << col + 1 << "] = "
                                         << for_sym.str() << endl;
              else
                {
                  sparseHelper(3, third_derivatives_output, k+k2, 0, output_type);
                  third_derivatives_output << "=" << eq + 1 << ";" << endl;

                  sparseHelper(3, third_derivatives_output, k+k2, 1, output_type);
                  third_derivatives_output << "=" << col + 1 << ";" << endl;

                  sparseHelper(3, third_derivatives_output, k+k2, 2, output_type);
                  third_derivatives_output << "=";
                  sparseHelper(3, third_derivatives_output, k, 2, output_type);
                  third_derivatives_output << ";" << endl;

                  k2++;
                }
          k += k2;
        }
    }

  if (output_type == ExprNodeOutputType::matlabStaticModel)
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
      init_output << "residual = zeros(" << equations.size() << ", 1);";
      end_output << "if ~isreal(residual)" << endl
                 << "  residual = real(residual)+imag(residual).^2;" << endl
                 << "end";
      writeStaticModelHelper(basename, "static_resid", "residual", "static_resid_tt",
                             temporary_terms_mlv.size() + temporary_terms_res.size(),
                             "", init_output, end_output,
                             model_output, model_tt_output);

      init_output.str(string());
      init_output.clear();
      end_output.str(string());
      end_output.clear();
      init_output << "g1 = zeros(" << equations.size() << ", " << symbol_table.endo_nbr() << ");";
      end_output << "if ~isreal(g1)" << endl
                 << "    g1 = real(g1)+2*imag(g1);" << endl
                 << "end";
      writeStaticModelHelper(basename, "static_g1", "g1", "static_g1_tt",
                             temporary_terms_mlv.size() + temporary_terms_res.size() + temporary_terms_g1.size(),
                             "static_resid_tt",
                             init_output, end_output,
                             jacobian_output, jacobian_tt_output);
      writeWrapperFunctions(basename, "g1");

      init_output.str(string());
      init_output.clear();
      end_output.str(string());
      end_output.clear();
      if (second_derivatives.size())
        {
          init_output << "v2 = zeros(" << NNZDerivatives[1] << ",3);";
          end_output << "g2 = sparse(v2(:,1),v2(:,2),v2(:,3)," << equations.size() << "," << g2ncols << ");";
        }
      else
        init_output << "g2 = sparse([],[],[]," << equations.size() << "," << g2ncols << ");";
      writeStaticModelHelper(basename, "static_g2", "g2", "static_g2_tt",
                             temporary_terms_mlv.size() + temporary_terms_res.size() + temporary_terms_g1.size()
                             + temporary_terms_g2.size(),
                             "static_g1_tt",
                             init_output, end_output,
                             hessian_output, hessian_tt_output);
      writeWrapperFunctions(basename, "g2");

      init_output.str(string());
      init_output.clear();
      end_output.str(string());
      end_output.clear();
      int ncols = hessianColsNbr * JacobianColsNbr;
      if (third_derivatives.size())
        {
          init_output << "v3 = zeros(" << NNZDerivatives[2] << ",3);";
          end_output << "g3 = sparse(v3(:,1),v3(:,2),v3(:,3)," << nrows << "," << ncols << ");";
        }
      else
        init_output << "g3 = sparse([],[],[]," << nrows << "," << ncols << ");";
      writeStaticModelHelper(basename, "static_g3", "g3", "static_g3_tt",
                             temporary_terms_mlv.size() + temporary_terms_res.size() + temporary_terms_g1.size()
                             + temporary_terms_g2.size() + temporary_terms_g3.size(),
                             "static_g2_tt",
                             init_output, end_output,
                             third_derivatives_output, third_derivatives_tt_output);
      writeWrapperFunctions(basename, "g3");

      writeStaticMatlabCompatLayer(basename);
    }
  else if (output_type == ExprNodeOutputType::CStaticModel)
    {
      StaticOutput << "void static_resid_tt(const double *y, const double *x, int nb_row_x, const double *params, double *T)" << endl
                   << "{" << endl
                   << model_tt_output.str()
                   << "}" << endl
                   << endl
                   << "void static_resid(const double *y, const double *x, int nb_row_x, const double *params, const double *T, double *residual)" << endl
                   << "{" << endl
                   << "  double lhs, rhs;" << endl
                   << model_output.str()
                   << "}" << endl
                   << endl
                   << "void static_g1_tt(const double *y, const double *x, int nb_row_x, const double *params, double *T)" << endl
                   << "{" << endl
                   << jacobian_tt_output.str()
                   << "}" << endl
                   << endl
                   << "void static_g1(const double *y, const double *x, int nb_row_x, const double *params, const double *T, double *g1)" << endl
                   << "{" << endl
                   << jacobian_output.str()
                   << "}" << endl
                   << endl
                   << "void static_g2_tt(const double *y, const double *x, int nb_row_x, const double *params, double *T)" << endl
                   << "{" << endl
                   << hessian_tt_output.str()
                   << "}" << endl
                   << endl
                   << "void static_g2(const double *y, const double *x, int nb_row_x, const double *params, const double *T, double *v2)" << endl
                   << "{" << endl
                   << hessian_output.str()
                   << "}" << endl;
    }
  else
    {
      string filename = basename + "Static.jl";
      ofstream output;
      output.open(filename, ios::out | ios::binary);
      if (!output.is_open())
        {
          cerr << "Error: Can't open file " << filename << " for writing" << endl;
          exit(EXIT_FAILURE);
        }

      output << "module " << basename << "Static" << endl
             << "#" << endl
             << "# NB: this file was automatically generated by Dynare" << endl
             << "#     from " << basename << ".mod" << endl
             << "#" << endl
             << "using Utils" << endl << endl
             << "export tmp_nbr, static!, staticResid!, staticG1!, staticG2!, staticG3!" << endl << endl
             << "#=" << endl
             << "# The comments below apply to all functions contained in this module #" << endl
             << "  NB: The arguments contained on the first line of the function" << endl
             << "      definition are those that are modified in place" << endl << endl
             << "## Exported Functions ##" << endl
             << "  static!      : Wrapper function; computes residuals, Jacobian, Hessian," << endl
             << "                 and third order derivatives matroces depending on the arguments provided" << endl
             << "  staticResid! : Computes the static model residuals" << endl
             << "  staticG1!    : Computes the static model Jacobian" << endl
             << "  staticG2!    : Computes the static model Hessian" << endl
             << "  staticG3!    : Computes the static model third derivatives" << endl << endl
             << "## Exported Variables ##" << endl
             << "  tmp_nbr      : Vector{Int}(4) respectively the number of temporary variables" << endl
             << "                 for the residuals, g1, g2 and g3." << endl << endl
             << "## Local Functions ##" << endl
             << "  staticResidTT! : Computes the static model temporary terms for the residuals" << endl
             << "  staticG1TT!    : Computes the static model temporary terms for the Jacobian" << endl
             << "  staticG2TT!    : Computes the static model temporary terms for the Hessian" << endl
             << "  staticG3TT!    : Computes the static model temporary terms for the third derivatives" << endl << endl
             << "## Function Arguments ##" << endl
             << "  T        : Vector{Float64}(num_temp_terms) temporary terms" << endl
             << "  y        : Vector{Float64}(model_.endo_nbr) endogenous variables in declaration order" << endl
             << "  x        : Vector{Float64}(model_.exo_nbr) exogenous variables in declaration order" << endl
             << "  params   : Vector{Float64}(model_.param) parameter values in declaration order" << endl
             << "  residual : Vector{Float64}(model_.eq_nbr) residuals of the static model equations" << endl
             << "             in order of declaration of the equations. Dynare may prepend auxiliary equations," << endl
             << "             see model.aux_vars" << endl
             << "  g1       : Matrix{Float64}(model.eq_nbr,model_.endo_nbr) Jacobian matrix of the static model equations" << endl
             << "             The columns and rows respectively correspond to the variables in declaration order and the" << endl
             << "             equations in order of declaration" << endl
             << "  g2       : spzeros(model.eq_nbr, model_.endo^2) Hessian matrix of the static model equations" << endl
             << "             The columns and rows respectively correspond to the variables in declaration order and the" << endl
             << "             equations in order of declaration" << endl
             << "  g3       : spzeros(model.eq_nbr, model_.endo^3) Third order derivatives matrix of the static model equations" << endl
             << "             The columns and rows respectively correspond to the variables in declaration order and the" << endl
             << "             equations in order of declaration" << endl << endl
             << "## Remarks ##" << endl
             << "  [1] The size of `T`, ie the value of `num_temp_terms`, depends on the version of the static model called. The number of temporary variables" << endl    
             << "      used for the different returned objects (residuals, jacobian, hessian or third order derivatives) is given by the elements in `tmp_nbr`" << endl
             << "      exported vector. The first element is the number of temporaries used for the computation of the residuals, the second element is the" << endl
             << "      number of temporaries used for the evaluation of the jacobian matrix, etc. If one calls the version of the static model computing the" << endl
             << "      residuals, and the jacobian and hessian matrices, then `T` must have at least `sum(tmp_nbr[1:3])` elements." << endl
             << "=#" << endl << endl;

      // Write the number of temporary terms
      output << "tmp_nbr = zeros(Int,4)" << endl
             << "tmp_nbr[1] = " << temporary_terms_mlv.size() + temporary_terms_res.size() << "# Number of temporary terms for the residuals" << endl
             << "tmp_nbr[2] = " << temporary_terms_g1.size() << "# Number of temporary terms for g1 (jacobian)" << endl
             << "tmp_nbr[3] = " << temporary_terms_g2.size() << "# Number of temporary terms for g2 (hessian)" << endl
             << "tmp_nbr[4] = " << temporary_terms_g3.size() << "# Number of temporary terms for g3 (third order derivates)" << endl << endl;

      // staticResidTT!
      output << "function staticResidTT!(T::Vector{Float64}," << endl
             << "                        y::Vector{Float64}, x::Vector{Float64}, params::Vector{Float64})" << endl
             << "    @assert length(T) >= " << temporary_terms_mlv.size() + temporary_terms_res.size()  << endl
             << model_tt_output.str()
             << "    return nothing" << endl
             << "end" << endl << endl;

      // static!
      output << "function staticResid!(T::Vector{Float64}, residual::Vector{Float64}," << endl
             << "                      y::Vector{Float64}, x::Vector{Float64}, params::Vector{Float64}, T0_flag::Bool)" << endl
             << "    @assert length(y) == " << symbol_table.endo_nbr() << endl
             << "    @assert length(x) == " << symbol_table.exo_nbr() << endl
             << "    @assert length(params) == " << symbol_table.param_nbr() << endl
             << "    @assert length(residual) == " << equations.size() << endl
             << "    if T0_flag" << endl
             << "        staticResidTT!(T, y, x, params)" << endl
             << "    end" << endl
             << model_output.str()
             << "    if ~isreal(residual)" << endl
             << "        residual = real(residual)+imag(residual).^2;" << endl
             << "    end" << endl
             << "    return nothing" << endl
             << "end" << endl << endl;

      // staticG1TT!
      output << "function staticG1TT!(T::Vector{Float64}," << endl
             << "                     y::Vector{Float64}, x::Vector{Float64}, params::Vector{Float64}, T0_flag::Bool)" << endl
             << "    if T0_flag" << endl
             << "        staticResidTT!(T, y, x, params)" << endl
             << "    end" << endl
             << jacobian_tt_output.str()
             << "    return nothing" << endl
             << "end" << endl << endl;

      // staticG1!
      output << "function staticG1!(T::Vector{Float64}, g1::Matrix{Float64}," << endl
             << "                   y::Vector{Float64}, x::Vector{Float64}, params::Vector{Float64}, T1_flag::Bool, T0_flag::Bool)" << endl
             << "    @assert length(T) >= "
             << temporary_terms_mlv.size() + temporary_terms_res.size() + temporary_terms_g1.size() << endl
             << "    @assert size(g1) == (" << equations.size() << ", " << symbol_table.endo_nbr() << ")" << endl
             << "    @assert length(y) == " << symbol_table.endo_nbr() << endl
             << "    @assert length(x) == " << symbol_table.exo_nbr() << endl
             << "    @assert length(params) == " << symbol_table.param_nbr() << endl
             << "    if T1_flag" << endl
             << "        staticG1TT!(T, y, x, params, T0_flag)" << endl
             << "    end" << endl
             << "    fill!(g1, 0.0)" << endl
             << jacobian_output.str()
             << "    if ~isreal(g1)" << endl
             << "        g1 = real(g1)+2*imag(g1);" << endl
             << "    end" << endl
             << "    return nothing" << endl
             << "end" << endl << endl;

      // staticG2TT!
      output << "function staticG2TT!(T::Vector{Float64}," << endl
             << "                     y::Vector{Float64}, x::Vector{Float64}, params::Vector{Float64}, T1_flag::Bool, T0_flag::Bool)" << endl
             << "    if T1_flag" << endl
             << "        staticG1TT!(T, y, x, params, TO_flag)" << endl
             << "    end" << endl
             << hessian_tt_output.str()
             << "    return nothing" << endl
             << "end" << endl << endl;

      // staticG2!
      output << "function staticG2!(T::Vector{Float64}, g2::Matrix{Float64}," << endl
             << "                   y::Vector{Float64}, x::Vector{Float64}, params::Vector{Float64}, T2_flag::Bool, T1_flag::Bool, T0_flag::Bool)" << endl
             << "    @assert length(T) >= "
             << temporary_terms_mlv.size() + temporary_terms_res.size() + temporary_terms_g1.size() + temporary_terms_g2.size() << endl
             << "    @assert size(g2) == (" << equations.size() << ", " << g2ncols << ")" << endl
             << "    @assert length(y) == " << symbol_table.endo_nbr() << endl
             << "    @assert length(x) == " << symbol_table.exo_nbr() << endl
             << "    @assert length(params) == " << symbol_table.param_nbr() << endl
             << "    if T2_flag" << endl
             << "        staticG2TT!(T, y, x, params, T1_flag, T0_flag)" << endl
             << "    end" << endl
             << "    fill!(g2, 0.0)" << endl
             << hessian_output.str()
             << "    return nothing" << endl
             << "end" << endl << endl;

      // staticG3TT!
      output << "function staticG3TT!(T::Vector{Float64}," << endl
             << "                     y::Vector{Float64}, x::Vector{Float64}, params::Vector{Float64}, T2_flag::Bool, T1_flag::Bool, T0_flag::Bool)" << endl
             << "    if T2_flag" << endl
             << "        staticG2TT!(T, y, x, params, T1_flag, T0_flag)" << endl
             << "    end" << endl
             << third_derivatives_tt_output.str()
             << "    return nothing" << endl
             << "end" << endl << endl;

      // staticG3!
      int ncols = hessianColsNbr * JacobianColsNbr;
      output << "function staticG3!(T::Vector{Float64}, g3::Matrix{Float64}," << endl
             << "                   y::Vector{Float64}, x::Vector{Float64}, params::Vector{Float64}, T3_flag::Bool, T2_flag::Bool, T1_flag::Bool, T0_flag::Bool)" << endl
             << "    @assert length(T) >= "
             << temporary_terms_mlv.size() + temporary_terms_res.size() + temporary_terms_g1.size() + temporary_terms_g2.size() + temporary_terms_g3.size() << endl
             << "    @assert size(g3) == (" << nrows << ", " << ncols << ")" << endl
             << "    @assert length(y) == " << symbol_table.endo_nbr() << endl
             << "    @assert length(x) == " << symbol_table.exo_nbr() << endl
             << "    @assert length(params) == " << symbol_table.param_nbr() << endl
             << "    if T3_flag" << endl
             << "        staticG3TT!(T, y, x, params, T2_flag, T1_flag, T0_flag)" << endl
             << "    end" << endl
             << "    fill!(g3, 0.0)" << endl
             << third_derivatives_output.str()
             << "    return nothing" << endl
             << "end" << endl << endl;

      // static!
      output << "function static!(T::Vector{Float64}, residual::Vector{Float64}," << endl
             << "                  y::Vector{Float64}, x::Vector{Float64}, params::Vector{Float64})" << endl
             << "    staticResid!(T, residual, y, x, params, true)" << endl
             << "    return nothing" << endl
             << "end" << endl
             << endl
             << "function static!(T::Vector{Float64}, residual::Vector{Float64}, g1::Matrix{Float64}," << endl
             << "                 y::Vector{Float64}, x::Vector{Float64}, params::Vector{Float64})" << endl
             << "    staticG1!(T, g1, y, x, params, true)" << endl
             << "    staticResid!(T, residual, y, x, params, false)" << endl
             << "    return nothing" << endl
             << "end" << endl
             << endl
             << "function static!(T::Vector{Float64}, g1::Matrix{Float64}," << endl
             << "                 y::Vector{Float64}, x::Vector{Float64}, params::Vector{Float64})" << endl
             << "    staticG1!(T, g1, y, x, params, true, false)" << endl
             << "    return nothing" << endl
             << "end" << endl
             << endl
             << "function static!(T::Vector{Float64}, residual::Vector{Float64}, g1::Matrix{Float64}, g2::Matrix{Float64}," << endl
             << "                 y::Vector{Float64}, x::Vector{Float64}, params::Vector{Float64})" << endl
             << "    staticG2!(T, g2, y, x, params, true)" << endl
             << "    staticG1!(T, g1, y, x, params, false)" << endl
             << "    staticResid!(T, residual, y, x, params, false)" << endl
             << "    return nothing" << endl
             << "end" << endl << endl
             << "end" << endl;
      output.close();
    }
}

void
StaticModel::writeStaticCFile(const string &basename) const
{
  // Writing comments and function definition command
  boost::filesystem::create_directories(basename + "/model/src");
  string filename = basename + "/model/src/static.c";
  string filename_mex = basename + "/model/src/static_mex.c";

  int ntt = temporary_terms_mlv.size() + temporary_terms_res.size() + temporary_terms_g1.size() + temporary_terms_g2.size() + temporary_terms_g3.size();

  ofstream output;
  output.open(filename, ios::out | ios::binary);
  if (!output.is_open())
    {
      cerr << "ERROR: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }

  output << "/*" << endl
         << " * " << filename << " : Computes static model for Dynare" << endl
         << " *" << endl
         << " * Warning : this file is generated automatically by Dynare" << endl
         << " *           from model file (.mod)" << endl << endl
         << " */" << endl
#if defined(_WIN32) || defined(__CYGWIN32__)
         << "#ifdef _MSC_VER" << endl
         << "#define _USE_MATH_DEFINES" << endl
         << "#endif" << endl
#endif
         << "#include <math.h>" << endl;

  if (external_functions_table.get_total_number_of_unique_model_block_external_functions())
    // External Matlab function, implies Static function will call mex
    output << "#include \"mex.h\"" << endl;
  else
    output << "#include <stdlib.h>" << endl;

  output << "#define max(a, b) (((a) > (b)) ? (a) : (b))" << endl
         << "#define min(a, b) (((a) > (b)) ? (b) : (a))" << endl;

  // Write function definition if BinaryOpcode::powerDeriv is used
  writePowerDerivCHeader(output);
  writeNormcdfCHeader(output);

  output << endl;

  // Writing the function body
  writeStaticModel(output, true, false);

  output << endl;

  writePowerDeriv(output);
  writeNormcdf(output);
  output.close();

  output.open(filename_mex, ios::out | ios::binary);
  if (!output.is_open())
    {
      cerr << "ERROR: Can't open file " << filename_mex << " for writing" << endl;
      exit(EXIT_FAILURE);
    }

  // Writing the gateway routine
  output << "/*" << endl
         << " * " << filename_mex << " : The gateway routine used to call the Static function "
         << "located in " << filename << endl
         << " *" << endl
         << " * Warning : this file is generated automatically by Dynare" << endl
         << " *           from model file (.mod)" << endl << endl
         << " */" << endl
         << endl
         << "#include <stdlib.h>" << endl
         << "#include \"mex.h\"" << endl
         << endl
         << "const int ntt = " << ntt << ";" << endl
         << "void static_resid_tt(const double *y, const double *x, int nb_row_x, const double *params, double *T);" << endl
         << "void static_resid(const double *y, const double *x, int nb_row_x, const double *params, const double *T, double *residual);" << endl
         << "void static_g1_tt(const double *y, const double *x, int nb_row_x, const double *params, double *T);" << endl
         << "void static_g1(const double *y, const double *x, int nb_row_x, const double *params, const double *T, double *g1);" << endl
         << "void static_g2_tt(const double *y, const double *x, int nb_row_x, const double *params, double *T);" << endl
         << "void static_g2(const double *y, const double *x, int nb_row_x, const double *params, const double *T, double *v2);" << endl
         << "void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])" << endl
         << "{" << endl
         << "  /* Create a pointer to the input matrix y. */" << endl
         << "  double *y = mxGetPr(prhs[0]);" << endl
         << endl
         << "  /* Create a pointer to the input matrix x. */" << endl
         << "  double *x = mxGetPr(prhs[1]);" << endl
         << endl
         << "  /* Create a pointer to the input matrix params. */" << endl
         << "  double *params = mxGetPr(prhs[2]);" << endl
         << endl
         << "  /* Gets number of rows of matrix x. */" << endl
         << "  int nb_row_x = mxGetM(prhs[1]);" << endl
         << endl
         << "  double *T = (double *) malloc(sizeof(double)*ntt);"
         << endl
         << "  if (nlhs >= 1)" << endl
         << "    {" << endl
         << "      /* Set the output pointer to the output matrix residual. */" << endl
         << "      plhs[0] = mxCreateDoubleMatrix(" << equations.size() << ",1, mxREAL);" << endl
         << "      double *residual = mxGetPr(plhs[0]);" << endl
         << "      static_resid_tt(y, x, nb_row_x, params, T);" << endl
         << "      static_resid(y, x, nb_row_x, params, T, residual);" << endl
         << "    }" << endl
         << endl
         << "  if (nlhs >= 2)" << endl
         << "    {" << endl
         << "      /* Set the output pointer to the output matrix g1. */" << endl
         << "      plhs[1] = mxCreateDoubleMatrix(" << equations.size() << ", " << symbol_table.endo_nbr() << ", mxREAL);" << endl
         << "      double *g1 = mxGetPr(plhs[1]);" << endl
         << "      static_g1_tt(y, x, nb_row_x, params, T);" << endl
         << "      static_g1(y, x, nb_row_x, params, T, g1);" << endl
         << "    }" << endl
         << endl
         << "  if (nlhs >= 3)" << endl
         << "    {" << endl
         << "      /* Set the output pointer to the output matrix v2. */" << endl
         << "      plhs[2] = mxCreateDoubleMatrix(" << NNZDerivatives[1] << ", " << 3
         << ", mxREAL);" << endl
         << "      double *v2 = mxGetPr(plhs[2]);" << endl
         << "      static_g2_tt(y, x, nb_row_x, params, T);" << endl
         << "      static_g2(y, x, nb_row_x, params, T, v2);" << endl
         << "    }" << endl
         << endl
         << "  free(T);" << endl
         << "}" << endl;
  output.close();
}

void
StaticModel::writeStaticJuliaFile(const string &basename) const
{
  writeStaticModel(basename, false, true);
}

void
StaticModel::writeStaticFile(const string &basename, bool block, bool bytecode, bool use_dll, bool julia) const
{
  if (block && bytecode)
    writeModelEquationsCode_Block(basename, map_idx, map_idx2);
  else if (!block && bytecode)
    writeModelEquationsCode(basename, map_idx);
  else if (block && !bytecode)
    {
      writeModelEquationsOrdered_M(basename);
      writeStaticBlockMFSFile(basename);
    }
  else if (use_dll)
    writeStaticCFile(basename);
  else if (julia)
    writeStaticJuliaFile(basename);
  else
    writeStaticMFile(basename);
  writeSetAuxiliaryVariables(basename, julia);
}

bool
StaticModel::exoPresentInEqs() const
{
  for (auto equation : equations)
    if (equation->containsExogenous())
      return true;
  return false;
}

void
StaticModel::writeStaticBlockMFSFile(const string &basename) const
{
  string filename = packageDir(basename) + "/static.m";

  ofstream output;
  output.open(filename, ios::out | ios::binary);
  if (!output.is_open())
    {
      cerr << "ERROR: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }

  output << "function [residual, g1, y, var_index] = static(nblock, y, x, params)" << endl
         << "  residual = [];" << endl
         << "  g1 = [];" << endl
         << "  var_index = [];\n" << endl
         << "  switch nblock" << endl;

  unsigned int nb_blocks = getNbBlocks();

  for (int b = 0; b < (int) nb_blocks; b++)
    {

      set<int> local_var;

      output << "    case " << b+1 << endl;

      BlockSimulationType simulation_type = getBlockSimulationType(b);

      if (simulation_type == EVALUATE_BACKWARD || simulation_type == EVALUATE_FORWARD)
        {
          output << "      y_tmp = " << basename << ".block.static_" << b+1 << "(y, x, params);\n";
          ostringstream tmp;
          for (int i = 0; i < (int) getBlockSize(b); i++)
            tmp << " " << getBlockVariableID(b, i)+1;
          output << "      var_index = [" << tmp.str() << "];\n";
          output << "      residual  = y(var_index) - y_tmp(var_index);\n";
          output << "      y = y_tmp;\n";
        }
      else
        output << "      [residual, y, g1] = " << basename << ".block.static_" << b+1 << "(y, x, params);\n";

    }
  output << "  end" << endl
         << "end" << endl;
  output.close();
}

void
StaticModel::writeOutput(ostream &output, bool block) const
{
  output << "M_.static_tmp_nbr = zeros(4,1); % Number of temporaries used for the static model" <<endl
         << "M_.static_tmp_nbr(1) = " << temporary_terms_res.size() << "; % Number of temporaries used for the evaluation of the residuals" << endl
         << "M_.static_tmp_nbr(2) = " << temporary_terms_g1.size() << "; % Number of temporaries used for the evaluation of g1 (jacobian)" << endl
         << "M_.static_tmp_nbr(3) = " << temporary_terms_g2.size() << "; % Number of temporaries used for the evaluation of g2 (hessian)" << endl
         << "M_.static_tmp_nbr(4) = " << temporary_terms_g3.size() << "; % Number of temporaries used for the evaluation of g3 (third order derivatives)" << endl;

  if (!block)
    return;

  unsigned int nb_blocks = getNbBlocks();
  for (int b = 0; b < (int) nb_blocks; b++)
    {
      BlockSimulationType simulation_type = getBlockSimulationType(b);
      unsigned int block_size = getBlockSize(b);
      ostringstream tmp_s, tmp_s_eq;
      tmp_s.str("");
      tmp_s_eq.str("");
      for (unsigned int i = 0; i < block_size; i++)
        {
          tmp_s << " " << getBlockVariableID(b, i)+1;
          tmp_s_eq << " " << getBlockEquationID(b, i)+1;
        }
      output << "block_structure_stat.block(" << b+1 << ").Simulation_Type = " << simulation_type << ";\n";
      output << "block_structure_stat.block(" << b+1 << ").endo_nbr = " << block_size << ";\n";
      output << "block_structure_stat.block(" << b+1 << ").mfs = " << getBlockMfs(block) << ";\n";
      output << "block_structure_stat.block(" << b+1 << ").equation = [" << tmp_s_eq.str() << "];\n";
      output << "block_structure_stat.block(" << b+1 << ").variable = [" << tmp_s.str() << "];\n";
    }
  output << "M_.block_structure_stat.block = block_structure_stat.block;\n";
  string cst_s;
  int nb_endo = symbol_table.endo_nbr();
  output << "M_.block_structure_stat.variable_reordered = [";
  for (int i = 0; i < nb_endo; i++)
    output << " " << variable_reordered[i]+1;
  output << "];\n";
  output << "M_.block_structure_stat.equation_reordered = [";
  for (int i = 0; i < nb_endo; i++)
    output << " " << equation_reordered[i]+1;
  output << "];\n";

  map<pair<int, int>,  int>  row_incidence;
  for (const auto & first_derivative : first_derivatives)
    {
      int deriv_id = first_derivative.first.second;
      if (getTypeByDerivID(deriv_id) == SymbolType::endogenous)
        {
          int eq = first_derivative.first.first;
          int symb = getSymbIDByDerivID(deriv_id);
          int var = symbol_table.getTypeSpecificID(symb);
          //int lag = getLagByDerivID(deriv_id);
          row_incidence[{ eq, var }] = 1;
        }
    }
  output << "M_.block_structure_stat.incidence.sparse_IM = [";
  for (map<pair< int, int >,  int>::const_iterator it = row_incidence.begin(); it != row_incidence.end(); it++)
    {
      output << it->first.first+1 << " " << it->first.second+1 << ";\n";
    }
  output << "];\n";
}

SymbolType
StaticModel::getTypeByDerivID(int deriv_id) const noexcept(false)
{
  if (deriv_id < symbol_table.endo_nbr())
    return SymbolType::endogenous;
  else if (deriv_id < symbol_table.endo_nbr() + symbol_table.param_nbr())
    return SymbolType::parameter;
  else
    throw UnknownDerivIDException();
}

int
StaticModel::getLagByDerivID(int deriv_id) const noexcept(false)
{
  return 0;
}

int
StaticModel::getSymbIDByDerivID(int deriv_id) const noexcept(false)
{
  if (deriv_id < symbol_table.endo_nbr())
    return symbol_table.getID(SymbolType::endogenous, deriv_id);
  else if (deriv_id < symbol_table.endo_nbr() + symbol_table.param_nbr())
    return symbol_table.getID(SymbolType::parameter, deriv_id - symbol_table.endo_nbr());
  else
    throw UnknownDerivIDException();
}

int
StaticModel::getDerivID(int symb_id, int lag) const noexcept(false)
{
  if (symbol_table.getType(symb_id) == SymbolType::endogenous)
    return symbol_table.getTypeSpecificID(symb_id);
  else if (symbol_table.getType(symb_id) == SymbolType::parameter)
    return symbol_table.getTypeSpecificID(symb_id) + symbol_table.endo_nbr();
  else
    return -1;
}

void
StaticModel::addAllParamDerivId(set<int> &deriv_id_set)
{
  for (int i = 0; i < symbol_table.param_nbr(); i++)
    deriv_id_set.insert(i + symbol_table.endo_nbr());
}

map<pair<pair<int, pair<int, int>>, pair<int, int>>, int>
StaticModel::get_Derivatives(int block)
{
  map<pair<pair<int, pair<int, int>>, pair<int, int>>, int> Derivatives;
  Derivatives.clear();
  int block_size = getBlockSize(block);
  int block_nb_recursive = block_size - getBlockMfs(block);
  int lag = 0;
  for (int eq = 0; eq < block_size; eq++)
    {
      int eqr = getBlockEquationID(block, eq);
      for (int var = 0; var < block_size; var++)
        {
          int varr = getBlockVariableID(block, var);
          if (dynamic_jacobian.find({ lag, { eqr, varr } }) != dynamic_jacobian.end())
            {
              bool OK = true;
              map<pair<pair<int, pair<int, int>>, pair<int, int>>, int>::const_iterator its = Derivatives.find({ { lag, { eq, var } }, { eqr, varr } });
              if (its != Derivatives.end())
                {
                  if (its->second == 2)
                    OK = false;
                }

              if (OK)
                {
                  if (getBlockEquationType(block, eq) == E_EVALUATE_S && eq < block_nb_recursive)
                    //It's a normalized equation, we have to recompute the derivative using chain rule derivative function
                    Derivatives[{ { lag, { eq, var } }, { eqr, varr } }] = 1;
                  else
                    //It's a feedback equation we can use the derivatives
                    Derivatives[{ { lag, { eq, var } }, { eqr, varr } }] = 0;
                }
              if (var < block_nb_recursive)
                {
                  int eqs = getBlockEquationID(block, var);
                  for (int vars = block_nb_recursive; vars < block_size; vars++)
                    {
                      int varrs = getBlockVariableID(block, vars);
                      //A new derivative needs to be computed using the chain rule derivative function (a feedback variable appears in a recursive equation)
                      if (Derivatives.find({ { lag, { var, vars } }, { eqs, varrs } }) != Derivatives.end())
                        Derivatives[{ { lag, { eq, vars } }, { eqr, varrs } }] = 2;
                    }
                }
            }
        }
    }

  return (Derivatives);
}

void
StaticModel::computeChainRuleJacobian(blocks_derivatives_t &blocks_derivatives)
{
  map<int, expr_t> recursive_variables;
  unsigned int nb_blocks = getNbBlocks();
  blocks_derivatives = blocks_derivatives_t(nb_blocks);
  for (unsigned int block = 0; block < nb_blocks; block++)
    {
      block_derivatives_equation_variable_laglead_nodeid_t tmp_derivatives;
      recursive_variables.clear();
      BlockSimulationType simulation_type = getBlockSimulationType(block);
      int block_size = getBlockSize(block);
      int block_nb_mfs = getBlockMfs(block);
      int block_nb_recursives = block_size - block_nb_mfs;
      if (simulation_type == SOLVE_TWO_BOUNDARIES_COMPLETE || simulation_type == SOLVE_TWO_BOUNDARIES_SIMPLE)
        {
          blocks_derivatives.push_back(block_derivatives_equation_variable_laglead_nodeid_t(0));
          for (int i = 0; i < block_nb_recursives; i++)
            {
              if (getBlockEquationType(block, i) == E_EVALUATE_S)
                recursive_variables[getDerivID(symbol_table.getID(SymbolType::endogenous, getBlockVariableID(block, i)), 0)] = getBlockEquationRenormalizedExpr(block, i);
              else
                recursive_variables[getDerivID(symbol_table.getID(SymbolType::endogenous, getBlockVariableID(block, i)), 0)] = getBlockEquationExpr(block, i);
            }
          map<pair<pair<int, pair<int, int>>, pair<int, int>>, int> Derivatives = get_Derivatives(block);
          map<pair<pair<int, pair<int, int>>, pair<int, int>>, int>::const_iterator it = Derivatives.begin();
          for (int i = 0; i < (int) Derivatives.size(); i++)
            {
              int Deriv_type = it->second;
              pair<pair<int, pair<int, int>>, pair<int, int>> it_l(it->first);
              it++;
              int lag = it_l.first.first;
              int eq = it_l.first.second.first;
              int var = it_l.first.second.second;
              int eqr = it_l.second.first;
              int varr = it_l.second.second;
              if (Deriv_type == 0)
                first_chain_rule_derivatives[{ eqr, { varr, lag } }] = first_derivatives[{ eqr, getDerivID(symbol_table.getID(SymbolType::endogenous, varr), lag) }];
              else if (Deriv_type == 1)
                first_chain_rule_derivatives[{ eqr, { varr, lag } }] = (equation_type_and_normalized_equation[eqr].second)->getChainRuleDerivative(getDerivID(symbol_table.getID(SymbolType::endogenous, varr), lag), recursive_variables);
              else if (Deriv_type == 2)
                {
                  if (getBlockEquationType(block, eq) == E_EVALUATE_S && eq < block_nb_recursives)
                    first_chain_rule_derivatives[{ eqr, { varr, lag } }] = (equation_type_and_normalized_equation[eqr].second)->getChainRuleDerivative(getDerivID(symbol_table.getID(SymbolType::endogenous, varr), lag), recursive_variables);
                  else
                    first_chain_rule_derivatives[{ eqr, { varr, lag } }] = equations[eqr]->getChainRuleDerivative(getDerivID(symbol_table.getID(SymbolType::endogenous, varr), lag), recursive_variables);
                }
              tmp_derivatives.emplace_back(make_pair(eq, var), make_pair(lag, first_chain_rule_derivatives[make_pair(eqr, make_pair(varr, lag))]));
            }
        }
      else
        {
          blocks_derivatives.push_back(block_derivatives_equation_variable_laglead_nodeid_t(0));
          for (int i = 0; i < block_nb_recursives; i++)
            {
              if (getBlockEquationType(block, i) == E_EVALUATE_S)
                recursive_variables[getDerivID(symbol_table.getID(SymbolType::endogenous, getBlockVariableID(block, i)), 0)] = getBlockEquationRenormalizedExpr(block, i);
              else
                recursive_variables[getDerivID(symbol_table.getID(SymbolType::endogenous, getBlockVariableID(block, i)), 0)] = getBlockEquationExpr(block, i);
            }
          for (int eq = block_nb_recursives; eq < block_size; eq++)
            {
              int eqr = getBlockEquationID(block, eq);
              for (int var = block_nb_recursives; var < block_size; var++)
                {
                  int varr = getBlockVariableID(block, var);
                  expr_t d1 = equations[eqr]->getChainRuleDerivative(getDerivID(symbol_table.getID(SymbolType::endogenous, varr), 0), recursive_variables);
                  if (d1 == Zero)
                    continue;
                  first_chain_rule_derivatives[{ eqr, { varr, 0 } }] = d1;
                  tmp_derivatives.emplace_back(make_pair(eq, var), make_pair(0, first_chain_rule_derivatives[make_pair(eqr, make_pair(varr, 0))]));
                }
            }
        }
      blocks_derivatives[block] = tmp_derivatives;
    }
}

void
StaticModel::collect_block_first_order_derivatives()
{
  //! vector for an equation or a variable indicates the block number
  vector<int> equation_2_block, variable_2_block;
  unsigned int nb_blocks = getNbBlocks();
  equation_2_block = vector<int>(equation_reordered.size());
  variable_2_block = vector<int>(variable_reordered.size());
  for (unsigned int block = 0; block < nb_blocks; block++)
    {
      unsigned int block_size = getBlockSize(block);
      for (unsigned int i = 0; i < block_size; i++)
        {
          equation_2_block[getBlockEquationID(block, i)] = block;
          variable_2_block[getBlockVariableID(block, i)] = block;
        }
    }
  derivative_endo = vector<derivative_t>(nb_blocks);
  endo_max_leadlag_block = vector<pair<int, int>>(nb_blocks, { 0, 0 });
  max_leadlag_block = vector<pair<int, int>>(nb_blocks, { 0, 0 });
  for (auto & first_derivative : first_derivatives)
    {
      int eq = first_derivative.first.first;
      int var = symbol_table.getTypeSpecificID(getSymbIDByDerivID(first_derivative.first.second));
      int lag = 0;
      int block_eq = equation_2_block[eq];
      int block_var = variable_2_block[var];
      max_leadlag_block[block_eq] = { 0, 0 };
      max_leadlag_block[block_eq] = { 0, 0 };
      endo_max_leadlag_block[block_eq] = { 0, 0 };
      endo_max_leadlag_block[block_eq] = { 0, 0 };
      derivative_t tmp_derivative;
      lag_var_t lag_var;
      if (getTypeByDerivID(first_derivative.first.second) == SymbolType::endogenous && block_eq == block_var)
        {
          tmp_derivative = derivative_endo[block_eq];
          tmp_derivative[{ lag, { eq, var } }] = first_derivatives[{ eq, getDerivID(symbol_table.getID(SymbolType::endogenous, var), lag) }];
          derivative_endo[block_eq] = tmp_derivative;
        }
    }
}

void
StaticModel::writeLatexFile(const string &basename, bool write_equation_tags) const
{
  writeLatexModelFile(basename + "_static", ExprNodeOutputType::latexStaticModel, write_equation_tags);
}

void
StaticModel::writeAuxVarInitval(ostream &output, ExprNodeOutputType output_type) const
{
  for (auto aux_equation : aux_equations)
    {
      dynamic_cast<ExprNode *>(aux_equation)->writeOutput(output, output_type);
      output << ";" << endl;
    }
}

void
StaticModel::writeSetAuxiliaryVariables(const string &basename, const bool julia) const
{
  ostringstream output_func_body;
  writeAuxVarRecursiveDefinitions(output_func_body, ExprNodeOutputType::matlabStaticModel);

  if (output_func_body.str().empty())
    return;

  string func_name = julia ? basename + "_set_auxiliary_variables" : "set_auxiliary_variables";
  string filename = julia ? func_name + ".jl" : packageDir(basename) + "/" + func_name + ".m";
  string comment = julia ? "#" : "%";

  ofstream output;
  output.open(filename, ios::out | ios::binary);
  if (!output.is_open())
    {
      cerr << "ERROR: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }

  output << "function y = " << func_name + "(y, x, params)" << endl
         << comment << endl
         << comment << " Status : Computes static model for Dynare" << endl
         << comment << endl
         << comment << " Warning : this file is generated automatically by Dynare" << endl
         << comment << "           from model file (.mod)" << endl << endl
         << output_func_body.str();

  output.close();
}

void
StaticModel::writeAuxVarRecursiveDefinitions(ostream &output, ExprNodeOutputType output_type) const
{
  deriv_node_temp_terms_t tef_terms;
  for (auto aux_equation : aux_equations)
    if (dynamic_cast<ExprNode *>(aux_equation)->containsExternalFunction())
      dynamic_cast<ExprNode *>(aux_equation)->writeExternalFunctionOutput(output, ExprNodeOutputType::matlabStaticModel, {}, {}, tef_terms);
  for (auto aux_equation : aux_equations)
    {
      dynamic_cast<ExprNode *>(aux_equation->substituteStaticAuxiliaryDefinition())->writeOutput(output, output_type);
      output << ";" << endl;
    }
}

void
StaticModel::writeLatexAuxVarRecursiveDefinitions(ostream &output) const
{
  deriv_node_temp_terms_t tef_terms;
  temporary_terms_t temporary_terms;
    temporary_terms_idxs_t temporary_terms_idxs;
  for (auto aux_equation : aux_equations)
    if (dynamic_cast<ExprNode *>(aux_equation)->containsExternalFunction())
      dynamic_cast<ExprNode *>(aux_equation)->writeExternalFunctionOutput(output, ExprNodeOutputType::latexStaticModel,
                                                                              temporary_terms, temporary_terms_idxs, tef_terms);
  for (auto aux_equation : aux_equations)
    {
      output << "\\begin{dmath}" << endl;
      dynamic_cast<ExprNode *>(aux_equation->substituteStaticAuxiliaryDefinition())->writeOutput(output, ExprNodeOutputType::latexStaticModel);
      output << endl << "\\end{dmath}" << endl;
    }
}

void
StaticModel::writeJsonAuxVarRecursiveDefinitions(ostream &output) const
{
  deriv_node_temp_terms_t tef_terms;
  temporary_terms_t temporary_terms;

  for (auto aux_equation : aux_equations)
    if (dynamic_cast<ExprNode *>(aux_equation)->containsExternalFunction())
      {
        vector<string> efout;
        dynamic_cast<ExprNode *>(aux_equation)->writeJsonExternalFunctionOutput(efout,
                                                                                    temporary_terms,
                                                                                    tef_terms,
                                                                                    false);
        for (vector<string>::const_iterator it = efout.begin(); it != efout.end(); it++)
          {
            if (it != efout.begin())
              output << ", ";
            output << *it;
          }
      }

  for (auto aux_equation : aux_equations)
    {
      output << ", {\"lhs\": \"";
      aux_equation->get_arg1()->writeJsonOutput(output, temporary_terms, tef_terms, false);
      output << "\", \"rhs\": \"";
      dynamic_cast<BinaryOpNode *>(aux_equation->substituteStaticAuxiliaryDefinition())->get_arg2()->writeJsonOutput(output, temporary_terms, tef_terms, false);
      output << "\"}";
    }
}

void
StaticModel::writeParamsDerivativesFile(const string &basename, bool julia) const
{
  if (!residuals_params_derivatives.size()
      && !residuals_params_second_derivatives.size()
      && !jacobian_params_derivatives.size()
      && !jacobian_params_second_derivatives.size()
      && !hessian_params_derivatives.size())
    return;

  ExprNodeOutputType output_type = (julia ? ExprNodeOutputType::juliaStaticModel : ExprNodeOutputType::matlabStaticModel);

  ostringstream model_local_vars_output;   // Used for storing model local vars
  ostringstream model_output;              // Used for storing model
  ostringstream jacobian_output;           // Used for storing jacobian equations
  ostringstream hessian_output;            // Used for storing Hessian equations
  ostringstream hessian1_output;           // Used for storing Hessian equations
  ostringstream third_derivs_output;       // Used for storing third order derivatives equations
  ostringstream third_derivs1_output;      // Used for storing third order derivatives equations

  deriv_node_temp_terms_t tef_terms;

  writeTemporaryTerms(params_derivs_temporary_terms, {}, params_derivs_temporary_terms_idxs, model_output, output_type, tef_terms);

  for (const auto & residuals_params_derivative : residuals_params_derivatives)
    {
      int eq, param;
      tie(eq, param) = residuals_params_derivative.first;
      expr_t d1 = residuals_params_derivative.second;

      int param_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param)) + 1;

      jacobian_output << "rp" << LEFT_ARRAY_SUBSCRIPT(output_type)
                      <<  eq+1 << ", " << param_col
                      << RIGHT_ARRAY_SUBSCRIPT(output_type) << " = ";
      d1->writeOutput(jacobian_output, output_type, params_derivs_temporary_terms, params_derivs_temporary_terms_idxs, tef_terms);
      jacobian_output << ";" << endl;
    }

  for (const auto & jacobian_params_derivative : jacobian_params_derivatives)
    {
      int eq, var, param;
      tie(eq, var, param) = jacobian_params_derivative.first;
      expr_t d2 = jacobian_params_derivative.second;

      int var_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(var)) + 1;
      int param_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param)) + 1;

      hessian_output << "gp" << LEFT_ARRAY_SUBSCRIPT(output_type)
                     << eq+1 << ", " << var_col << ", " << param_col
                     << RIGHT_ARRAY_SUBSCRIPT(output_type) << " = ";
      d2->writeOutput(hessian_output, output_type, params_derivs_temporary_terms, params_derivs_temporary_terms_idxs, tef_terms);
      hessian_output << ";" << endl;
    }

  int i = 1;
  for (const auto &it : residuals_params_second_derivatives)
    {
      int eq, param1, param2;
      tie(eq, param1, param2) = it.first;
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
      d2->writeOutput(hessian1_output, output_type, params_derivs_temporary_terms, params_derivs_temporary_terms_idxs, tef_terms);
      hessian1_output << ";" << endl;

      i++;
    }

  i = 1;
  for (const auto &it : jacobian_params_second_derivatives)
    {
      int eq, var, param1, param2;
      tie(eq, var, param1, param2) = it.first;
      expr_t d2 = it.second;

      int var_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(var)) + 1;
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
      d2->writeOutput(third_derivs_output, output_type, params_derivs_temporary_terms, params_derivs_temporary_terms_idxs, tef_terms);
      third_derivs_output << ";" << endl;

      i++;
    }

  i = 1;
  for (const auto &it : hessian_params_derivatives)
    {
      int eq, var1, var2, param;
      tie(eq, var1, var2, param) = it.first;
      expr_t d2 = it.second;

      int var1_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(var1)) + 1;
      int var2_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(var2)) + 1;
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
      d2->writeOutput(third_derivs1_output, output_type, params_derivs_temporary_terms, params_derivs_temporary_terms_idxs, tef_terms);
      third_derivs1_output << ";" << endl;

      i++;
    }

  ofstream paramsDerivsFile;
  string filename = julia ? basename + "StaticParamsDerivs.jl" : packageDir(basename) + "/static_params_derivs.m";
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
      fixNestedParenthesis(model_output, tmp_paren_vars, message_printed);
      fixNestedParenthesis(model_local_vars_output, tmp_paren_vars, message_printed);
      fixNestedParenthesis(jacobian_output, tmp_paren_vars, message_printed);
      fixNestedParenthesis(hessian_output, tmp_paren_vars, message_printed);
      fixNestedParenthesis(hessian1_output, tmp_paren_vars, message_printed);
      fixNestedParenthesis(third_derivs_output, tmp_paren_vars, message_printed);
      fixNestedParenthesis(third_derivs1_output, tmp_paren_vars, message_printed);

      paramsDerivsFile << "function [rp, gp, rpp, gpp, hp] = static_params_derivs(y, x, params)" << endl
                       << "%" << endl
                       << "% Status : Computes derivatives of the static model with respect to the parameters" << endl
                       << "%" << endl
                       << "% Inputs : " << endl
                       << "%   y         [M_.endo_nbr by 1] double    vector of endogenous variables in declaration order" << endl
                       << "%   x         [M_.exo_nbr by 1] double     vector of exogenous variables in declaration order" << endl
                       << "%   params    [M_.param_nbr by 1] double   vector of parameter values in declaration order" << endl
                       << "%" << endl
                       << "% Outputs:" << endl
                       << "%   rp        [M_.eq_nbr by #params] double    Jacobian matrix of static model equations with respect to parameters " << endl
                       << "%                                              Dynare may prepend or append auxiliary equations, see M_.aux_vars" << endl
                       << "%   gp        [M_.endo_nbr by M_.endo_nbr by #params] double    Derivative of the Jacobian matrix of the static model equations with respect to the parameters" << endl
                       << "%                                                           rows: variables in declaration order" << endl
                       << "%                                                           rows: equations in order of declaration" << endl
                       << "%   rpp       [#second_order_residual_terms by 4] double   Hessian matrix of second derivatives of residuals with respect to parameters;" << endl
                       << "%                                                              rows: respective derivative term" << endl
                       << "%                                                              1st column: equation number of the term appearing" << endl
                       << "%                                                              2nd column: number of the first parameter in derivative" << endl
                       << "%                                                              3rd column: number of the second parameter in derivative" << endl
                       << "%                                                              4th column: value of the Hessian term" << endl
                       << "%   gpp      [#second_order_Jacobian_terms by 5] double   Hessian matrix of second derivatives of the Jacobian with respect to the parameters;" << endl
                       << "%                                                              rows: respective derivative term" << endl
                       << "%                                                              1st column: equation number of the term appearing" << endl
                       << "%                                                              2nd column: column number of variable in Jacobian of the static model" << endl
                       << "%                                                              3rd column: number of the first parameter in derivative" << endl
                       << "%                                                              4th column: number of the second parameter in derivative" << endl
                       << "%                                                              5th column: value of the Hessian term" << endl
                       << "%" << endl
                       << "%" << endl
                       << "% Warning : this file is generated automatically by Dynare" << endl
                       << "%           from model file (.mod)" << endl << endl
                       << "T = NaN(" << params_derivs_temporary_terms_idxs.size() << ",1);" << endl
                       << model_local_vars_output.str()
                       << model_output.str()
                       << "rp = zeros(" << equations.size() << ", "
                       << symbol_table.param_nbr() << ");" << endl
                       << jacobian_output.str()
                       << "gp = zeros(" << equations.size() << ", " << symbol_table.endo_nbr() << ", "
                       << symbol_table.param_nbr() << ");" << endl
                       << hessian_output.str()
                       << "if nargout >= 3" << endl
                       << "rpp = zeros(" << residuals_params_second_derivatives.size() << ",4);" << endl
                       << hessian1_output.str()
                       << "gpp = zeros(" << jacobian_params_second_derivatives.size() << ",5);" << endl
                       << third_derivs_output.str()
                       << "end" << endl
                       << "if nargout >= 5" << endl
                       << "hp = zeros(" << hessian_params_derivatives.size() << ",5);" << endl
                       << third_derivs1_output.str()
                       << "end" << endl
                       << "end" << endl;
    }
  else
    paramsDerivsFile << "module " << basename << "StaticParamsDerivs" << endl
                     << "#" << endl
                     << "# NB: this file was automatically generated by Dynare" << endl
                     << "#     from " << basename << ".mod" << endl
                     << "#" << endl
                     << "export params_derivs" << endl << endl
                     << "function params_derivs(y, x, params)" << endl
                     << model_local_vars_output.str()
                     << model_output.str()
                     << "rp = zeros(" << equations.size() << ", "
                     << symbol_table.param_nbr() << ");" << endl
                     << jacobian_output.str()
                     << "gp = zeros(" << equations.size() << ", " << symbol_table.endo_nbr() << ", "
                     << symbol_table.param_nbr() << ");" << endl
                     << hessian_output.str()
                     << "rpp = zeros(" << residuals_params_second_derivatives.size() << ",4);" << endl
                     << hessian1_output.str()
                     << "gpp = zeros(" << jacobian_params_second_derivatives.size() << ",5);" << endl
                     << third_derivs_output.str()
                     << "hp = zeros(" << hessian_params_derivatives.size() << ",5);" << endl
                     << third_derivs1_output.str()
                     << "(rp, gp, rpp, gpp, hp)" << endl
                     << "end" << endl
                     << "end" << endl;

  paramsDerivsFile.close();
}

void
StaticModel::writeJsonOutput(ostream &output) const
{
  writeJsonModelEquations(output, false);
}

void
StaticModel::writeJsonComputingPassOutput(ostream &output, bool writeDetails) const
{
  ostringstream model_local_vars_output;   // Used for storing model local vars
  ostringstream model_output;              // Used for storing model
  ostringstream jacobian_output;           // Used for storing jacobian equations
  ostringstream hessian_output;            // Used for storing Hessian equations
  ostringstream third_derivatives_output;  // Used for storing third order derivatives equations

  deriv_node_temp_terms_t tef_terms;
  temporary_terms_t temp_term_union = temporary_terms_res;
  temporary_terms_t temp_term_union_m_1;

  string concat = "";

  writeJsonModelLocalVariables(model_local_vars_output, tef_terms);

  writeJsonTemporaryTerms(temporary_terms_res, temp_term_union_m_1, model_output, tef_terms, concat);
  model_output << ", ";
  writeJsonModelEquations(model_output, true);

  int nrows = equations.size();
  int JacobianColsNbr = symbol_table.endo_nbr();
  int hessianColsNbr = JacobianColsNbr*JacobianColsNbr;

  // Write Jacobian w.r. to endogenous only
  temp_term_union_m_1 = temp_term_union;
  temp_term_union.insert(temporary_terms_g1.begin(), temporary_terms_g1.end());
  concat = "jacobian";
  writeJsonTemporaryTerms(temp_term_union, temp_term_union_m_1, jacobian_output, tef_terms, concat);
  jacobian_output << ", \"jacobian\": {"
                  << "  \"nrows\": " << nrows
                  << ", \"ncols\": " << JacobianColsNbr
                  << ", \"entries\": [";
  for (auto it = first_derivatives.begin();
       it != first_derivatives.end(); it++)
    {
      if (it != first_derivatives.begin())
        jacobian_output << ", ";

      int eq, var;
      tie(eq, var) = it->first;
      int symb_id = getSymbIDByDerivID(var);
      int col = symbol_table.getTypeSpecificID(symb_id);
      expr_t d1 = it->second;

      if (writeDetails)
        jacobian_output << "{\"eq\": " << eq + 1;
      else
        jacobian_output << "{\"row\": " << eq + 1;

      jacobian_output << ", \"col\": " << col + 1;

      if (writeDetails)
        jacobian_output << ", \"var\": \"" << symbol_table.getName(symb_id) << "\"";

      jacobian_output << ", \"val\": \"";
      d1->writeJsonOutput(jacobian_output, temp_term_union, tef_terms);
      jacobian_output << "\"}" << endl;
    }
  jacobian_output << "]}";

  int g2ncols = symbol_table.endo_nbr() * symbol_table.endo_nbr();
  // Write Hessian w.r. to endogenous only (only if 2nd order derivatives have been computed)
  temp_term_union_m_1 = temp_term_union;
  temp_term_union.insert(temporary_terms_g2.begin(), temporary_terms_g2.end());
  concat = "hessian";
  writeJsonTemporaryTerms(temp_term_union, temp_term_union_m_1, hessian_output, tef_terms, concat);
  hessian_output << ", \"hessian\": {"
                 << "  \"nrows\": " << equations.size()
                 << ", \"ncols\": " << g2ncols
                 << ", \"entries\": [";
  for (auto it = second_derivatives.begin();
       it != second_derivatives.end(); it++)
    {
      if (it != second_derivatives.begin())
        hessian_output << ", ";

      int eq, var1, var2;
      tie(eq, var1, var2) = it->first;
      int symb_id1 = getSymbIDByDerivID(var1);
      int symb_id2 = getSymbIDByDerivID(var2);
      expr_t d2 = it->second;

      int tsid1 = symbol_table.getTypeSpecificID(symb_id1);
      int tsid2 = symbol_table.getTypeSpecificID(symb_id2);

      int col = tsid1*symbol_table.endo_nbr()+tsid2;
      int col_sym = tsid2*symbol_table.endo_nbr()+tsid1;

      if (writeDetails)
        hessian_output << "{\"eq\": " << eq + 1;
      else
        hessian_output << "{\"row\": " << eq + 1;

      hessian_output << ", \"col\": [" << col + 1;

      if (writeDetails)
        hessian_output << ", \"var1\": \"" << symbol_table.getName(symb_id1) << "\""
                       << ", \"var2\": \"" << symbol_table.getName(symb_id2) << "\"";

      if (symb_id1 != symb_id2)
        hessian_output << ", " <<  col_sym + 1;
      hessian_output << "]"
                     << ", \"val\": \"";
      d2->writeJsonOutput(hessian_output, temp_term_union, tef_terms);
      hessian_output << "\"}" << endl;
    }
  hessian_output << "]}";

  // Writing third derivatives
  temp_term_union_m_1 = temp_term_union;
  temp_term_union.insert(temporary_terms_g3.begin(), temporary_terms_g3.end());
  concat = "third_derivatives";
  writeJsonTemporaryTerms(temp_term_union, temp_term_union_m_1, third_derivatives_output, tef_terms, concat);
  third_derivatives_output << ", \"third_derivative\": {"
                           << "  \"nrows\": " << equations.size()
                           << ", \"ncols\": " << hessianColsNbr * JacobianColsNbr
                           << ", \"entries\": [";
  for (auto it = third_derivatives.begin();
       it != third_derivatives.end(); it++)
    {
      if (it != third_derivatives.begin())
        third_derivatives_output << ", ";

      int eq, var1, var2, var3;
      tie(eq, var1, var2, var3) = it->first;
      expr_t d3 = it->second;

      if (writeDetails)
        third_derivatives_output << "{\"eq\": " << eq + 1;
      else
        third_derivatives_output << "{\"row\": " << eq + 1;

      int id1 = getSymbIDByDerivID(var1);
      int id2 = getSymbIDByDerivID(var2);
      int id3 = getSymbIDByDerivID(var3);
      set<int> cols;
      cols.insert(id1 * hessianColsNbr + id2 * JacobianColsNbr + id3);
      cols.insert(id1 * hessianColsNbr + id3 * JacobianColsNbr + id2);
      cols.insert(id2 * hessianColsNbr + id1 * JacobianColsNbr + id3);
      cols.insert(id2 * hessianColsNbr + id3 * JacobianColsNbr + id1);
      cols.insert(id3 * hessianColsNbr + id1 * JacobianColsNbr + id2);
      cols.insert(id3 * hessianColsNbr + id2 * JacobianColsNbr + id1);

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
                                 << ", \"var2\": \"" << symbol_table.getName(getSymbIDByDerivID(var2)) << "\""
                                 << ", \"var3\": \"" << symbol_table.getName(getSymbIDByDerivID(var3)) << "\"";

      third_derivatives_output << ", \"val\": \"";
      d3->writeJsonOutput(third_derivatives_output, temp_term_union, tef_terms);
      third_derivatives_output << "\"}" << endl;
    }
  third_derivatives_output << "]}";

  if (writeDetails)
    output << "\"static_model\": {";
  else
    output << "\"static_model_simple\": {";
  output << model_local_vars_output.str()
         << ", " << model_output.str()
         << ", " << jacobian_output.str()
         << ", " << hessian_output.str()
         << ", " << third_derivatives_output.str()
         << "}";
}

void
StaticModel::writeJsonParamsDerivativesFile(ostream &output, bool writeDetails) const
{
  if (!residuals_params_derivatives.size()
      && !residuals_params_second_derivatives.size()
      && !jacobian_params_derivatives.size()
      && !jacobian_params_second_derivatives.size()
      && !hessian_params_derivatives.size())
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

  temporary_terms_t temp_terms_empty;
  string concat = "all";
  writeJsonTemporaryTerms(params_derivs_temporary_terms, temp_terms_empty, model_output, tef_terms, concat);
  jacobian_output << "\"deriv_wrt_params\": {"
                  << "  \"neqs\": " << equations.size()
                  << ", \"nparamcols\": " << symbol_table.param_nbr()
                  << ", \"entries\": [";
  for (auto it = residuals_params_derivatives.begin();
       it != residuals_params_derivatives.end(); it++)
    {
      if (it != residuals_params_derivatives.begin())
        jacobian_output << ", ";

      int eq, param;
      tie(eq, param) = it->first;
      expr_t d1 = it->second;

      int param_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param)) + 1;

      if (writeDetails)
        jacobian_output << "{\"eq\": " << eq + 1;
      else
        jacobian_output << "{\"row\": " << eq + 1;

      if (writeDetails)
        jacobian_output << ", \"param_col\": " << param_col;

      jacobian_output << ", \"param\": \"" << symbol_table.getName(getSymbIDByDerivID(param)) << "\"";

      jacobian_output << ", \"val\": \"";
      d1->writeJsonOutput(jacobian_output, params_derivs_temporary_terms, tef_terms);
      jacobian_output << "\"}" << endl;
    }
  jacobian_output << "]}";
  hessian_output << "\"deriv_jacobian_wrt_params\": {"
                 << "  \"neqs\": " << equations.size()
                 << ", \"nvarcols\": " << symbol_table.endo_nbr()
                 << ", \"nparamcols\": " << symbol_table.param_nbr()
                 << ", \"entries\": [";
  for (auto it = jacobian_params_derivatives.begin();
       it != jacobian_params_derivatives.end(); it++)
    {
      if (it != jacobian_params_derivatives.begin())
        hessian_output << ", ";

      int eq, var, param;
      tie(eq, var, param) = it->first;
      expr_t d2 = it->second;

      int var_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(var)) + 1;
      int param_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param)) + 1;

      if (writeDetails)
        hessian_output << "{\"eq\": " << eq + 1;
      else
        hessian_output << "{\"row\": " << eq + 1;

      if (writeDetails)
        hessian_output << ", \"var\": \"" << symbol_table.getName(getSymbIDByDerivID(var)) << "\""
                       << ", \"param\": \"" << symbol_table.getName(getSymbIDByDerivID(param)) << "\"";

      hessian_output << ", \"var_col\": " << var_col
                     << ", \"param_col\": " << param_col
                     << ", \"val\": \"";
      d2->writeJsonOutput(hessian_output, params_derivs_temporary_terms, tef_terms);
      hessian_output << "\"}" << endl;
    }
  hessian_output << "]}";

  hessian1_output << "\"second_deriv_residuals_wrt_params\": {"
                  << "  \"nrows\": " << equations.size()
                  << ", \"nparam1cols\": " << symbol_table.param_nbr()
                  << ", \"nparam2cols\": " << symbol_table.param_nbr()
                  << ", \"entries\": [";
  for (auto it = residuals_params_second_derivatives.begin();
       it != residuals_params_second_derivatives.end(); ++it)
    {
      if (it != residuals_params_second_derivatives.begin())
        hessian1_output << ", ";

      int eq, param1, param2;
      tie(eq, param1, param2) = it->first;
      expr_t d2 = it->second;

      int param1_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param1)) + 1;
      int param2_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param2)) + 1;

      if (writeDetails)
        hessian1_output << "{\"eq\": " << eq + 1;
      else
        hessian1_output << "{\"row\": " << eq + 1;

      hessian1_output << ", \"param1_col\": " << param1_col
                      << ", \"param2_col\": " << param2_col;

      if (writeDetails)
        hessian1_output << ", \"param1\": \"" << symbol_table.getName(getSymbIDByDerivID(param1)) << "\""
                        << ", \"param2\": \"" << symbol_table.getName(getSymbIDByDerivID(param2)) << "\"";

      hessian1_output << ", \"val\": \"";
      d2->writeJsonOutput(hessian1_output, params_derivs_temporary_terms, tef_terms);
      hessian1_output << "\"}" << endl;
    }
  hessian1_output << "]}";
  third_derivs_output << "\"second_deriv_jacobian_wrt_params\": {"
                      << "  \"neqs\": " << equations.size()
                      << ", \"nvarcols\": " << symbol_table.endo_nbr()
                      << ", \"nparam1cols\": " << symbol_table.param_nbr()
                      << ", \"nparam2cols\": " << symbol_table.param_nbr()
                      << ", \"entries\": [";
  for (auto it = jacobian_params_second_derivatives.begin();
       it != jacobian_params_second_derivatives.end(); ++it)
    {
      if (it != jacobian_params_second_derivatives.begin())
        third_derivs_output << ", ";

      int eq, var, param1, param2;
      tie(eq, var, param1, param2) = it->first;
      expr_t d2 = it->second;

      int var_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(var)) + 1;
      int param1_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param1)) + 1;
      int param2_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param2)) + 1;

      if (writeDetails)
        third_derivs_output << "{\"eq\": " << eq + 1;
      else
        third_derivs_output << "{\"row\": " << eq + 1;
      third_derivs_output << ", \"var_col\": " << var_col
                          << ", \"param1_col\": " << param1_col
                          << ", \"param2_col\": " << param2_col;

      if (writeDetails)
        third_derivs_output << ", \"var\": \"" << symbol_table.getName(var) << "\""
                            << ", \"param1\": \"" << symbol_table.getName(getSymbIDByDerivID(param1)) << "\""
                            << ", \"param2\": \"" << symbol_table.getName(getSymbIDByDerivID(param2)) << "\"";

      third_derivs_output << ", \"val\": \"";
      d2->writeJsonOutput(third_derivs_output, params_derivs_temporary_terms, tef_terms);
      third_derivs_output << "\"}" << endl;
    }
  third_derivs_output << "]}" << endl;

  third_derivs1_output << "\"derivative_hessian_wrt_params\": {"
                       << "  \"neqs\": " << equations.size()
                       << ", \"nvar1cols\": " << symbol_table.endo_nbr()
                       << ", \"nvar2cols\": " << symbol_table.endo_nbr()
                       << ", \"nparamcols\": " << symbol_table.param_nbr()
                       << ", \"entries\": [";
  for (auto it = hessian_params_derivatives.begin();
       it != hessian_params_derivatives.end(); ++it)
    {
      if (it != hessian_params_derivatives.begin())
        third_derivs1_output << ", ";

      int eq, var1, var2, param;
      tie(eq, var1, var2, param) = it->first;
      expr_t d2 = it->second;

      int var1_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(var1)) + 1;
      int var2_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(var2)) + 1;
      int param_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param)) + 1;

      if (writeDetails)
        third_derivs1_output << "{\"eq\": " << eq + 1;
      else
        third_derivs1_output << "{\"row\": " << eq + 1;

      third_derivs1_output << ", \"var1_col\": " << var1_col
                           << ", \"var2_col\": " << var2_col
                           << ", \"param_col\": " << param_col;

      if (writeDetails)
        third_derivs1_output << ", \"var1\": \"" << symbol_table.getName(getSymbIDByDerivID(var1)) << "\""
                             << ", \"var2\": \"" << symbol_table.getName(getSymbIDByDerivID(var2)) << "\""
                             << ", \"param1\": \"" << symbol_table.getName(getSymbIDByDerivID(param)) << "\"";

      third_derivs1_output << ", \"val\": \"";
      d2->writeJsonOutput(third_derivs1_output, params_derivs_temporary_terms, tef_terms);
      third_derivs1_output << "\"}" << endl;
    }
  third_derivs1_output << "]}" << endl;

  if (writeDetails)
    output << "\"static_model_params_derivative\": {";
  else
    output << "\"static_model_params_derivatives_simple\": {";
  output << model_local_vars_output.str()
         << ", " << model_output.str()
         << ", " << jacobian_output.str()
         << ", " << hessian_output.str()
         << ", " << hessian1_output.str()
         << ", " << third_derivs_output.str()
         << ", " << third_derivs1_output.str()
         << "}";
}
