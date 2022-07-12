/*
 * Copyright © 2003-2022 Dynare Team
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
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <algorithm>
#include <sstream>

#include "StaticModel.hh"
#include "DynamicModel.hh"

StaticModel::StaticModel(SymbolTable &symbol_table_arg,
                         NumericalConstants &num_constants_arg,
                         ExternalFunctionsTable &external_functions_table_arg) :
  ModelTree{symbol_table_arg, num_constants_arg, external_functions_table_arg}
{
}

StaticModel::StaticModel(const StaticModel &m) :
  ModelTree{m}
{
}

StaticModel &
StaticModel::operator=(const StaticModel &m)
{
  ModelTree::operator=(m);

  return *this;
}

StaticModel::StaticModel(const DynamicModel &m) :
  ModelTree{m.symbol_table, m.num_constants, m.external_functions_table}
{
  // Convert model local variables (need to be done first)
  for (int it : m.local_variables_vector)
    AddLocalVariable(it, m.local_variables_table.find(it)->second->toStatic(*this));

  // Convert equations
  int static_only_index = 0;
  set<int> dynamic_equations = m.equation_tags.getDynamicEqns();
  for (int i = 0; i < static_cast<int>(m.equations.size()); i++)
    try
      {
        // If equation is dynamic, replace it by an equation marked [static]
        if (dynamic_equations.contains(i))
          {
            auto [static_only_equations,
                  static_only_equations_lineno,
                  static_only_equations_equation_tags] = m.getStaticOnlyEquationsInfo();

            addEquation(static_only_equations[static_only_index]->toStatic(*this),
                        static_only_equations_lineno[static_only_index],
                        static_only_equations_equation_tags.getTagsByEqn(static_only_index));
            static_only_index++;
          }
        else
          addEquation(m.equations[i]->toStatic(*this),
                      m.equations_lineno[i],
                      m.equation_tags.getTagsByEqn(i));
      }
    catch (DataTree::DivisionByZeroException)
      {
        cerr << "...division by zero error encountered when converting equation " << i << " to static" << endl;
        exit(EXIT_FAILURE);
      }

  // Convert auxiliary equations
  for (auto aux_eq : m.aux_equations)
    addAuxEquation(aux_eq->toStatic(*this));

  user_set_add_flags = m.user_set_add_flags;
  user_set_subst_flags = m.user_set_subst_flags;
  user_set_add_libs = m.user_set_add_libs;
  user_set_subst_libs = m.user_set_subst_libs;
  user_set_compiler = m.user_set_compiler;
}

void
StaticModel::writeBytecodeDerivative(BytecodeWriter &code_file, int eq, int symb_id, const temporary_terms_t &temporary_terms, const temporary_terms_idxs_t &temporary_terms_idxs, const deriv_node_temp_terms_t &tef_terms) const
{
  if (auto it = derivatives[1].find({ eq, getDerivID(symbol_table.getID(SymbolType::endogenous, symb_id), 0) });
      it != derivatives[1].end())
    it->second->writeBytecodeOutput(code_file, ExprNodeBytecodeOutputType::staticModel, temporary_terms, temporary_terms_idxs, tef_terms);
  else
    code_file << FLDZ_{};
}

void
StaticModel::writeBytecodeChainRuleDerivative(BytecodeWriter &code_file, int blk, int eq, int var, int lag, const temporary_terms_t &temporary_terms, const temporary_terms_idxs_t &temporary_terms_idxs, const deriv_node_temp_terms_t &tef_terms) const
{
  if (auto it = blocks_derivatives[blk].find({ eq, var, lag });
      it != blocks_derivatives[blk].end())
    it->second->writeBytecodeOutput(code_file, ExprNodeBytecodeOutputType::staticModel, temporary_terms, temporary_terms_idxs, tef_terms);
  else
    code_file << FLDZ_{};
}

void
StaticModel::writeStaticPerBlockMFiles(const string &basename) const
{
  temporary_terms_t temporary_terms; // Temp terms written so far

  for (int blk = 0; blk < static_cast<int>(blocks.size()); blk++)
    {
      BlockSimulationType simulation_type = blocks[blk].simulation_type;

      string filename = packageDir(basename + ".block") + "/static_" + to_string(blk+1) + ".m";
      ofstream output{filename, ios::out | ios::binary};
      if (!output.is_open())
        {
          cerr << "ERROR: Can't open file " << filename << " for writing" << endl;
          exit(EXIT_FAILURE);
        }
      output << "%" << endl
             << "% " << filename << " : Computes static version of one block" << endl
             << "%" << endl
             << "% Warning : this file is generated automatically by Dynare" << endl
             << "%           from model file (.mod)" << endl << endl
             << "%" << endl;
      if (simulation_type == BlockSimulationType::evaluateBackward
          || simulation_type == BlockSimulationType::evaluateForward)
        output << "function [y, T] = static_" << blk+1 << "(y, x, params, T)" << endl;
      else
        output << "function [residual, y, T, g1] = static_" << blk+1 << "(y, x, params, T)" << endl;

      output << "  % ////////////////////////////////////////////////////////////////////////" << endl
             << "  % //" << "                     Block "s.substr(static_cast<int>(log10(blk + 1))) << blk+1
             << "                                        //" << endl
             << "  % //                     Simulation type "
             << BlockSim(simulation_type) << "  //" << endl
             << "  % ////////////////////////////////////////////////////////////////////////" << endl;

      if (simulation_type != BlockSimulationType::evaluateBackward
          && simulation_type != BlockSimulationType::evaluateForward)
        output << "  residual=zeros(" << blocks[blk].mfs_size << ",1);" << endl
               << "  g1_i=zeros(" << blocks_derivatives[blk].size() << ",1);" << endl
               << "  g1_j=zeros(" << blocks_derivatives[blk].size() << ",1);" << endl
               << "  g1_v=zeros(" << blocks_derivatives[blk].size() << ",1);" << endl
               << endl;

      writeStaticPerBlockHelper<ExprNodeOutputType::matlabStaticModel>(blk, output, temporary_terms);

      if (simulation_type != BlockSimulationType::evaluateBackward
          && simulation_type != BlockSimulationType::evaluateForward)
        output << endl
               << "  g1=sparse(g1_i, g1_j, g1_v, "  << blocks[blk].mfs_size << "," << blocks[blk].mfs_size << ");" << endl;

      output << "end" << endl;
      output.close();
    }
}

void
StaticModel::writeStaticPerBlockCFiles(const string &basename) const
{
  temporary_terms_t temporary_terms; // Temp terms written so far

  for (int blk = 0; blk < static_cast<int>(blocks.size()); blk++)
    {
      BlockSimulationType simulation_type = blocks[blk].simulation_type;

      string filename = basename + "/model/src/static_" + to_string(blk+1) + ".c";
      ofstream output{filename, ios::out | ios::binary};
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
             << R"(#include "mex.h")" << endl
             << endl;

      // Write function definition if BinaryOpcode::powerDeriv is used
      writePowerDerivHeader(output);

      output << endl;

      if (simulation_type == BlockSimulationType::evaluateBackward
          || simulation_type == BlockSimulationType::evaluateForward)
        output << "void static_" << blk+1 << "(double *restrict y, const double *restrict x, const double *restrict params, double *restrict T)" << endl;
      else
        output << "void static_" << blk+1 << "(double *restrict y, const double *restrict x, const double *restrict params, double *restrict T, double *restrict residual, double *restrict g1_i, double *restrict g1_j, double *restrict g1_v)" << endl;
      output << '{' << endl;

      writeStaticPerBlockHelper<ExprNodeOutputType::CStaticModel>(blk, output, temporary_terms);

      output << '}' << endl
             << endl;

      ostringstream header;
      if (simulation_type == BlockSimulationType::evaluateBackward
          || simulation_type == BlockSimulationType::evaluateForward)
        {
          header << "void static_" << blk+1 << "_mx(mxArray *y, const mxArray *x, const mxArray *params, mxArray *T)";
          output << header.str() << endl
                 << '{' << endl
                 << "  static_" << blk+1 << "(mxGetPr(y), mxGetPr(x), mxGetPr(params), mxGetPr(T));" << endl
                 << '}' << endl;
        }
      else
        {
          header << "void static_" << blk+1 << "_mx(mxArray *y, const mxArray *x, const mxArray *params, mxArray *T, mxArray **residual, mxArray **g1)";
          output << header.str() << endl
                 << '{' << endl
                 << "  *residual = mxCreateDoubleMatrix(" << blocks[blk].mfs_size << ",1,mxREAL);" << endl
                 << "  mxArray *g1_i = mxCreateDoubleMatrix(" << blocks_derivatives[blk].size() << ",1,mxREAL);" << endl
                 << "  mxArray *g1_j = mxCreateDoubleMatrix(" << blocks_derivatives[blk].size() << ",1,mxREAL);" << endl
                 << "  mxArray *g1_v = mxCreateDoubleMatrix(" << blocks_derivatives[blk].size() << ",1,mxREAL);" << endl
                 << "  static_" << blk+1 << "(mxGetPr(y), mxGetPr(x), mxGetPr(params), mxGetPr(T), mxGetPr(*residual), mxGetPr(g1_i), mxGetPr(g1_j), mxGetPr(g1_v));" << endl
                 << "  mxArray *plhs[1];" << endl
                 << "  mxArray *m = mxCreateDoubleScalar(" << blocks[blk].mfs_size << ");" << endl
                 << "  mxArray *n = mxCreateDoubleScalar(" << blocks[blk].mfs_size << ");" << endl
                 << "  mxArray *prhs[5] = { g1_i, g1_j, g1_v, m, n };" << endl
                 << R"(  mexCallMATLAB(1, plhs, 5, prhs, "sparse");)" << endl
                 << "  *g1 = plhs[0];" << endl
                 << "  mxDestroyArray(g1_i);" << endl
                 << "  mxDestroyArray(g1_j);" << endl
                 << "  mxDestroyArray(g1_v);" << endl
                 << "  mxDestroyArray(m);" << endl
                 << "  mxDestroyArray(n);" << endl
                 << '}' << endl;
        }

      output.close();

      filename = basename + "/model/src/static_" + to_string(blk+1) + ".h";
      ofstream header_output{filename, ios::out | ios::binary};
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
StaticModel::writeStaticBytecode(const string &basename) const
{
  ostringstream tmp_output;
  bool file_open = false;

  BytecodeWriter code_file{basename + "/model/bytecode/static.cod"};
  int count_u;
  int u_count_int = 0;

  writeBytecodeBinFile(basename + "/model/bytecode/static.bin", u_count_int, file_open, false);
  file_open = true;

  // Compute the union of temporary terms from residuals and 1st derivatives
  temporary_terms_t temporary_terms = temporary_terms_derivatives[0];
  temporary_terms.insert(temporary_terms_derivatives[1].begin(), temporary_terms_derivatives[1].end());

  //Temporary variables declaration
  code_file << FDIMST_{static_cast<int>(temporary_terms.size())}
    << FBEGINBLOCK_{symbol_table.endo_nbr(),
      BlockSimulationType::solveForwardComplete,
      0,
      symbol_table.endo_nbr(),
      endo_idx_block2orig,
      eq_idx_block2orig,
      false,
      symbol_table.endo_nbr(),
      0,
      0,
      u_count_int,
      symbol_table.endo_nbr()};

  temporary_terms_t temporary_terms_union;
  deriv_node_temp_terms_t tef_terms;

  writeBytecodeTemporaryTerms(code_file, ExprNodeBytecodeOutputType::staticModel, temporary_terms_union, temporary_terms_idxs, tef_terms);

  writeBytecodeModelEquations(code_file, ExprNodeBytecodeOutputType::staticModel, temporary_terms_union, temporary_terms_idxs, tef_terms);

  code_file << FENDEQU_{};

  // Get the current code_file position and jump if evaluating
  int pos_jmpifeval = code_file.getInstructionCounter();
  code_file << FJMPIFEVAL_{0}; // Use 0 as jump offset for the time being

  vector<vector<pair<int, int>>> my_derivatives(symbol_table.endo_nbr());
  count_u = symbol_table.endo_nbr();
  for (const auto & [indices, d1] : derivatives[1])
    {
      int deriv_id = indices[1];
      if (getTypeByDerivID(deriv_id) == SymbolType::endogenous)
        {
          int eq = indices[0];
          int var { getTypeSpecificIDByDerivID(deriv_id) };
          code_file << FNUMEXPR_{ExpressionType::FirstEndoDerivative, eq, var};
          if (!my_derivatives[eq].size())
            my_derivatives[eq].clear();
          my_derivatives[eq].emplace_back(var, count_u);

          d1->writeBytecodeOutput(code_file, ExprNodeBytecodeOutputType::staticModel, temporary_terms_union, temporary_terms_idxs, tef_terms);

          code_file << FSTPSU_{count_u};
          count_u++;
        }
    }
  for (int i = 0; i < symbol_table.endo_nbr(); i++)
    {
      code_file << FLDR_{i};
      if (my_derivatives[i].size())
        {
          for (bool printed_something{false};
               const auto &it : my_derivatives[i])
            {
              code_file << FLDSU_{it.second}
                << FLDSV_{SymbolType::endogenous, it.first}
                << FBINARY_{BinaryOpcode::times};
              if (exchange(printed_something, true))
                code_file << FBINARY_{BinaryOpcode::plus};
            }
          code_file << FBINARY_{BinaryOpcode::minus};
        }
      code_file << FSTPSU_{i};
    }

  // Jump unconditionally after the block
  int pos_jmp = code_file.getInstructionCounter();
  code_file << FJMP_{0}; // Use 0 as jump offset for the time being
  // Update jump offset for previous JMPIFEVAL
  code_file.overwriteInstruction(pos_jmpifeval, FJMPIFEVAL_{pos_jmp-pos_jmpifeval});

  temporary_terms_t tt2, tt3;

  // The Jacobian if we have to solve the block determinsitic bloc
  for (const auto & [indices, d1] : derivatives[1])
    {
      int deriv_id = indices[1];
      if (getTypeByDerivID(deriv_id) == SymbolType::endogenous)
        {
          int eq = indices[0];
          int var { getTypeSpecificIDByDerivID(deriv_id) };
          code_file << FNUMEXPR_{ExpressionType::FirstEndoDerivative, eq, var};
          if (!my_derivatives[eq].size())
            my_derivatives[eq].clear();
          my_derivatives[eq].emplace_back(var, count_u);

          d1->writeBytecodeOutput(code_file, ExprNodeBytecodeOutputType::staticModel, temporary_terms_union, temporary_terms_idxs, tef_terms);
          code_file << FSTPG2_{eq, var};
        }
    }

  // Update jump offset for previous JMP
  int pos_end_block = code_file.getInstructionCounter();
  code_file.overwriteInstruction(pos_jmp, FJMP_{pos_end_block-pos_jmp-1});

  code_file << FENDBLOCK_{} << FEND_{};
}

void
StaticModel::writeStaticBlockBytecode(const string &basename) const
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
  expr_t lhs = nullptr, rhs = nullptr;
  BinaryOpNode *eq_node;
  Uff Uf[symbol_table.endo_nbr()];
  map<expr_t, int> reference_count;
  vector<int> feedback_variables;
  bool file_open = false;

  BytecodeWriter code_file{basename + "/model/bytecode/static.cod"};

  //Temporary variables declaration
  code_file << FDIMST_{static_cast<int>(blocks_temporary_terms_idxs.size())};

  temporary_terms_t temporary_terms_union;

  for (int block = 0; block < static_cast<int>(blocks.size()); block++)
    {
      feedback_variables.clear();
      if (block > 0)
        code_file << FENDBLOCK_{};
      int count_u;
      int u_count_int = 0;
      BlockSimulationType simulation_type = blocks[block].simulation_type;
      int block_size = blocks[block].size;
      int block_mfs = blocks[block].mfs_size;
      int block_recursive = blocks[block].getRecursiveSize();

      if (simulation_type == BlockSimulationType::solveTwoBoundariesSimple
          || simulation_type == BlockSimulationType::solveTwoBoundariesComplete
          || simulation_type == BlockSimulationType::solveBackwardComplete
          || simulation_type == BlockSimulationType::solveForwardComplete)
        {
          writeBlockBytecodeBinFile(basename, block, u_count_int, file_open);
          file_open = true;
        }

      code_file << FBEGINBLOCK_{block_mfs,
          simulation_type,
          blocks[block].first_equation,
          block_size,
          endo_idx_block2orig,
          eq_idx_block2orig,
          blocks[block].linear,
          symbol_table.endo_nbr(),
          0,
          0,
          u_count_int,
          block_size};

      // Get the current code_file position and jump if evaluating
      int pos_jmpifeval = code_file.getInstructionCounter();
      code_file << FJMPIFEVAL_{0}; // Use 0 as jump offset for the time being

      //The Temporary terms
      deriv_node_temp_terms_t tef_terms;
      /* Keep a backup of temporary_terms_union here, since temp. terms are
         written a second time below. This is probably unwanted… */
      temporary_terms_t ttu_old = temporary_terms_union;

      auto write_eq_tt = [&](int eq)
                         {
                           for (auto it : blocks_temporary_terms[block][eq])
                             {
                               if (dynamic_cast<AbstractExternalFunctionNode *>(it))
                                 it->writeBytecodeExternalFunctionOutput(code_file, ExprNodeBytecodeOutputType::staticModel, temporary_terms_union, blocks_temporary_terms_idxs, tef_terms);

                               code_file << FNUMEXPR_{ExpressionType::TemporaryTerm, blocks_temporary_terms_idxs.at(it)};
                               it->writeBytecodeOutput(code_file, ExprNodeBytecodeOutputType::staticModel, temporary_terms_union, blocks_temporary_terms_idxs, tef_terms);
                               code_file << FSTPST_{blocks_temporary_terms_idxs.at(it)};
                               temporary_terms_union.insert(it);
                             }
                         };

      for (i = 0; i < block_size; i++)
        {
          write_eq_tt(i);

          // The equations
          int variable_ID, equation_ID;
          EquationType equ_type;
          switch (simulation_type)
            {
            evaluation:
            case BlockSimulationType::evaluateBackward:
            case BlockSimulationType::evaluateForward:
              equ_type = getBlockEquationType(block, i);
              code_file << FNUMEXPR_{ExpressionType::ModelEquation, getBlockEquationID(block, i)};
              if (equ_type == EquationType::evaluate)
                {
                  eq_node = getBlockEquationExpr(block, i);
                  lhs = eq_node->arg1;
                  rhs = eq_node->arg2;
                  rhs->writeBytecodeOutput(code_file, ExprNodeBytecodeOutputType::staticModel, temporary_terms_union, blocks_temporary_terms_idxs, tef_terms);
                  lhs->writeBytecodeOutput(code_file, ExprNodeBytecodeOutputType::staticAssignmentLHS, temporary_terms_union, blocks_temporary_terms_idxs, tef_terms);
                }
              else if (equ_type == EquationType::evaluateRenormalized)
                {
                  eq_node = getBlockEquationRenormalizedExpr(block, i);
                  lhs = eq_node->arg1;
                  rhs = eq_node->arg2;
                  rhs->writeBytecodeOutput(code_file, ExprNodeBytecodeOutputType::staticModel, temporary_terms_union, blocks_temporary_terms_idxs, tef_terms);
                  lhs->writeBytecodeOutput(code_file, ExprNodeBytecodeOutputType::staticAssignmentLHS, temporary_terms_union, blocks_temporary_terms_idxs, tef_terms);
                }
              break;
            case BlockSimulationType::solveBackwardComplete:
            case BlockSimulationType::solveForwardComplete:
              if (i < block_recursive)
                goto evaluation;
              variable_ID = getBlockVariableID(block, i);
              equation_ID = getBlockEquationID(block, i);
              feedback_variables.push_back(variable_ID);
              Uf[equation_ID].Ufl = nullptr;
              goto end;
            default:
            end:
              code_file << FNUMEXPR_{ExpressionType::ModelEquation, getBlockEquationID(block, i)};
              eq_node = getBlockEquationExpr(block, i);
              lhs = eq_node->arg1;
              rhs = eq_node->arg2;
              lhs->writeBytecodeOutput(code_file, ExprNodeBytecodeOutputType::staticModel, temporary_terms_union, blocks_temporary_terms_idxs, tef_terms);
              rhs->writeBytecodeOutput(code_file, ExprNodeBytecodeOutputType::staticModel, temporary_terms_union, blocks_temporary_terms_idxs, tef_terms);

              code_file << FBINARY_{BinaryOpcode::minus} << FSTPR_{i - block_recursive};
            }
        }
      code_file << FENDEQU_{};

      // The Jacobian if we have to solve the block
      if (simulation_type != BlockSimulationType::evaluateBackward
          && simulation_type != BlockSimulationType::evaluateForward)
        {
          // Write temporary terms for derivatives
          write_eq_tt(blocks[block].size);

          switch (simulation_type)
            {
            case BlockSimulationType::solveBackwardSimple:
            case BlockSimulationType::solveForwardSimple:
              code_file << FNUMEXPR_{ExpressionType::FirstEndoDerivative, 0, 0};
              writeBytecodeDerivative(code_file, getBlockEquationID(block, 0), getBlockVariableID(block, 0), temporary_terms_union, blocks_temporary_terms_idxs, tef_terms);
              code_file << FSTPG_{0};
              break;

            case BlockSimulationType::solveBackwardComplete:
            case BlockSimulationType::solveForwardComplete:
              count_u = feedback_variables.size();
              for (const auto &[indices, ignore2] : blocks_derivatives[block])
                {
                  auto [eq, var, ignore] = indices;
                  int eqr = getBlockEquationID(block, eq);
                  int varr = getBlockVariableID(block, var);
                  if (eq >= block_recursive && var >= block_recursive)
                    {
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
                      code_file << FNUMEXPR_{ExpressionType::FirstEndoDerivative, eqr, varr};
                      writeBytecodeChainRuleDerivative(code_file, block, eq, var, 0, temporary_terms_union, blocks_temporary_terms_idxs, tef_terms);
                      code_file << FSTPSU_{count_u};
                      count_u++;
                    }
                }
              for (i = 0; i < block_size; i++)
                if (i >= block_recursive)
                  {
                    code_file << FLDR_{i-block_recursive} << FLDZ_{};

                    v = getBlockEquationID(block, i);
                    for (Uf[v].Ufl = Uf[v].Ufl_First; Uf[v].Ufl; Uf[v].Ufl = Uf[v].Ufl->pNext)
                      code_file << FLDSU_{Uf[v].Ufl->u}
                        << FLDSV_{SymbolType::endogenous, Uf[v].Ufl->var}
                        << FBINARY_{BinaryOpcode::times}
                        << FCUML_{};
                    Uf[v].Ufl = Uf[v].Ufl_First;
                    while (Uf[v].Ufl)
                      {
                        Uf[v].Ufl_First = Uf[v].Ufl->pNext;
                        free(Uf[v].Ufl);
                        Uf[v].Ufl = Uf[v].Ufl_First;
                      }
                    code_file << FBINARY_{BinaryOpcode::minus}
                      << FSTPSU_{i - block_recursive};
                  }
              break;
            default:
              break;
            }
        }

      // Jump unconditionally after the block
      int pos_jmp = code_file.getInstructionCounter();
      code_file << FJMP_{0}; // Use 0 as jump offset for the time being
      // Update jump offset for previous JMPIFEVAL
      code_file.overwriteInstruction(pos_jmpifeval, FJMPIFEVAL_{pos_jmp-pos_jmpifeval});

      tef_terms.clear();
      temporary_terms_union = ttu_old;

      for (i = 0; i < block_size; i++)
        {
          write_eq_tt(i);

          // The equations
          int variable_ID, equation_ID;
          EquationType equ_type;
          switch (simulation_type)
            {
            evaluation_l:
            case BlockSimulationType::evaluateBackward:
            case BlockSimulationType::evaluateForward:
              equ_type = getBlockEquationType(block, i);
              code_file << FNUMEXPR_{ExpressionType::ModelEquation, getBlockEquationID(block, i)};
              if (equ_type == EquationType::evaluate)
                {
                  eq_node = getBlockEquationExpr(block, i);
                  lhs = eq_node->arg1;
                  rhs = eq_node->arg2;
                  rhs->writeBytecodeOutput(code_file, ExprNodeBytecodeOutputType::staticModel, temporary_terms_union, blocks_temporary_terms_idxs, tef_terms);
                  lhs->writeBytecodeOutput(code_file, ExprNodeBytecodeOutputType::staticAssignmentLHS, temporary_terms_union, blocks_temporary_terms_idxs, tef_terms);
                }
              else if (equ_type == EquationType::evaluateRenormalized)
                {
                  eq_node = getBlockEquationRenormalizedExpr(block, i);
                  lhs = eq_node->arg1;
                  rhs = eq_node->arg2;
                  rhs->writeBytecodeOutput(code_file, ExprNodeBytecodeOutputType::staticModel, temporary_terms_union, blocks_temporary_terms_idxs, tef_terms);
                  lhs->writeBytecodeOutput(code_file, ExprNodeBytecodeOutputType::staticAssignmentLHS, temporary_terms_union, blocks_temporary_terms_idxs, tef_terms);
                }
              break;
            case BlockSimulationType::solveBackwardComplete:
            case BlockSimulationType::solveForwardComplete:
              if (i < block_recursive)
                goto evaluation_l;
              variable_ID = getBlockVariableID(block, i);
              equation_ID = getBlockEquationID(block, i);
              feedback_variables.push_back(variable_ID);
              Uf[equation_ID].Ufl = nullptr;
              goto end_l;
            default:
            end_l:
              code_file << FNUMEXPR_{ExpressionType::ModelEquation, getBlockEquationID(block, i)};
              eq_node = getBlockEquationExpr(block, i);
              lhs = eq_node->arg1;
              rhs = eq_node->arg2;
              lhs->writeBytecodeOutput(code_file, ExprNodeBytecodeOutputType::staticModel, temporary_terms_union, blocks_temporary_terms_idxs, tef_terms);
              rhs->writeBytecodeOutput(code_file, ExprNodeBytecodeOutputType::staticModel, temporary_terms_union, blocks_temporary_terms_idxs, tef_terms);

              code_file << FBINARY_{BinaryOpcode::minus} << FSTPR_{i - block_recursive};
            }
        }
      code_file << FENDEQU_{};

      // The Jacobian if we have to solve the block determinsitic bloc

      // Write temporary terms for derivatives
      write_eq_tt(blocks[block].size);

      switch (simulation_type)
        {
        case BlockSimulationType::solveBackwardSimple:
        case BlockSimulationType::solveForwardSimple:
          code_file << FNUMEXPR_{ExpressionType::FirstEndoDerivative, 0, 0};
          writeBytecodeDerivative(code_file, getBlockEquationID(block, 0), getBlockVariableID(block, 0), temporary_terms_union, blocks_temporary_terms_idxs, tef_terms);
          code_file << FSTPG2_{0, 0};
          break;
        case BlockSimulationType::evaluateBackward:
        case BlockSimulationType::evaluateForward:
        case BlockSimulationType::solveBackwardComplete:
        case BlockSimulationType::solveForwardComplete:
          count_u = feedback_variables.size();
          for (const auto &[indices, ignore2] : blocks_derivatives[block])
            {
              auto &[eq, var, ignore] = indices;
              int eqr = getBlockEquationID(block, eq);
              int varr = getBlockVariableID(block, var);
              code_file << FNUMEXPR_{ExpressionType::FirstEndoDerivative, eqr, varr, 0};
              writeBytecodeChainRuleDerivative(code_file, block, eq, var, 0, temporary_terms_union, blocks_temporary_terms_idxs, tef_terms);
              code_file << FSTPG2_{eq, var};
            }
          break;
        default:
          break;
        }
      // Set codefile position to previous JMP_ and set the number of instructions to jump
      // Update jump offset for previous JMP
      int pos_end_block = code_file.getInstructionCounter();
      code_file.overwriteInstruction(pos_jmp, FJMP_{pos_end_block-pos_jmp-1});
    }
  code_file << FENDBLOCK_{} << FEND_{};
}

void
StaticModel::writeBlockBytecodeBinFile(const string &basename, int num,
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
  int block_size = blocks[num].size;
  int block_mfs = blocks[num].mfs_size;
  int block_recursive = blocks[num].getRecursiveSize();
  for (const auto &[indices, ignore2] : blocks_derivatives[num])
    {
      auto [eq, var, ignore] = indices;
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
StaticModel::computingPass(int derivsOrder, int paramsDerivsOrder, const eval_context_t &eval_context, bool no_tmp_terms, bool block)
{
  initializeVariablesAndEquations();

  vector<BinaryOpNode *> neweqs;
  for (int eq = 0; eq < static_cast<int>(equations.size() - aux_equations.size()); eq++)
    {
      expr_t eq_tmp = equations[eq]->substituteStaticAuxiliaryVariable();
      neweqs.push_back(dynamic_cast<BinaryOpNode *>(eq_tmp->toStatic(*this)));
    }

  for (auto &aux_equation : aux_equations)
    {
      expr_t eq_tmp = aux_equation->substituteStaticAuxiliaryDefinition();
      neweqs.push_back(dynamic_cast<BinaryOpNode *>(eq_tmp->toStatic(*this)));
    }

  equations.clear();
  copy(neweqs.begin(), neweqs.end(), back_inserter(equations));

  /* In both MATLAB and Julia, tensors for higher-order derivatives are stored
     in matrices whose columns correspond to variable multi-indices. Since we
     currently are limited to 32-bit signed integers (hence 31 bits) for matrix
     indices, check that we will not overflow (see #89). Note that such a check
     is not needed for parameter derivatives, since tensors for those are not
     stored as matrices. This check is implemented at this place for symmetry
     with DynamicModel::computingPass(). */
  if (log2(symbol_table.endo_nbr())*derivsOrder >= numeric_limits<int>::digits)
    {
      cerr << "ERROR: The static derivatives matrix is too large. Please decrease the approximation order." << endl;
      exit(EXIT_FAILURE);
    }

  // Compute derivatives w.r. to all endogenous
  set<int> vars;
  for (int i = 0; i < symbol_table.endo_nbr(); i++)
    {
      int id = symbol_table.getID(SymbolType::endogenous, i);
      vars.insert(getDerivID(id, 0));
    }

  // Launch computations
  cout << "Computing static model derivatives (order " << derivsOrder << ")." << endl;

  computeDerivatives(derivsOrder, vars);

  if (paramsDerivsOrder > 0)
    {
      cout << "Computing static model derivatives w.r.t. parameters (order " << paramsDerivsOrder << ")." << endl;
      computeParamsDerivatives(paramsDerivsOrder);
    }

  if (block)
    {
      auto contemporaneous_jacobian = evaluateAndReduceJacobian(eval_context);

      computeNonSingularNormalization(contemporaneous_jacobian, false);

      auto [prologue, epilogue] = computePrologueAndEpilogue();

      auto first_order_endo_derivatives = collectFirstOrderDerivativesEndogenous();

      equationTypeDetermination(first_order_endo_derivatives, mfs);

      cout << "Finding the optimal block decomposition of the model ..." << endl;

      computeBlockDecomposition(prologue, epilogue);

      reduceBlockDecomposition();

      printBlockDecomposition();

      computeChainRuleJacobian();

      determineLinearBlocks();

      if (!no_tmp_terms)
        computeBlockTemporaryTerms();
    }
  else
    {
      computeTemporaryTerms(true, no_tmp_terms);

      /* Must be called after computeTemporaryTerms(), because it depends on
         temporary_terms_mlv to be filled */
      if (paramsDerivsOrder > 0 && !no_tmp_terms)
        computeParamsDerivativesTemporaryTerms();
    }
}

void
StaticModel::writeStaticMFile(const string &basename) const
{
  auto [d_output, tt_output] = writeModelFileHelper<ExprNodeOutputType::matlabStaticModel>();

  ostringstream init_output, end_output;
  init_output << "residual = zeros(" << equations.size() << ", 1);";
  end_output << "if ~isreal(residual)" << endl
             << "  residual = real(residual)+imag(residual).^2;" << endl
             << "end";
  writeStaticMFileHelper(basename, "static_resid", "residual", "static_resid_tt",
                         temporary_terms_mlv.size() + temporary_terms_derivatives[0].size(),
                         "", init_output, end_output,
                         d_output[0], tt_output[0]);

  init_output.str("");
  end_output.str("");
  init_output << "g1 = zeros(" << equations.size() << ", " << symbol_table.endo_nbr() << ");";
  end_output << "if ~isreal(g1)" << endl
             << "    g1 = real(g1)+2*imag(g1);" << endl
             << "end";
  writeStaticMFileHelper(basename, "static_g1", "g1", "static_g1_tt",
                         temporary_terms_mlv.size() + temporary_terms_derivatives[0].size() + temporary_terms_derivatives[1].size(),
                         "static_resid_tt",
                         init_output, end_output,
                         d_output[1], tt_output[1]);
  writeStaticMWrapperFunction(basename, "g1");

  // For order ≥ 2
  int ncols{symbol_table.endo_nbr()};
  int ntt{static_cast<int>(temporary_terms_mlv.size() + temporary_terms_derivatives[0].size() + temporary_terms_derivatives[1].size())};
  for (size_t i{2}; i < derivatives.size(); i++)
    {
      ncols *= symbol_table.endo_nbr();
      ntt += temporary_terms_derivatives[i].size();
      string gname{"g" + to_string(i)};
      string gprevname{"g" + to_string(i-1)};

      init_output.str("");
      end_output.str("");
      if (derivatives[i].size())
        {
          init_output << gname << "_i = zeros(" << NNZDerivatives[i] << ",1);" << endl
                      << gname << "_j = zeros(" << NNZDerivatives[i] << ",1);" << endl
                      << gname << "_v = zeros(" << NNZDerivatives[i] << ",1);" << endl;
          end_output << gname << " = sparse("
                     << gname << "_i," << gname << "_j," << gname << "_v,"
                     << equations.size() << "," << ncols << ");";
        }
      else
        init_output << gname << " = sparse([],[],[]," << equations.size() << "," << ncols << ");";
      writeStaticMFileHelper(basename, "static_" + gname, gname,
                             "static_" + gname + "_tt",
                             ntt,
                             "static_" + gprevname + "_tt",
                             init_output, end_output,
                             d_output[i], tt_output[i]);
      if (i <= 3)
        writeStaticMWrapperFunction(basename, gname);
    }

  writeStaticMCompatFile(basename);
}

void
StaticModel::writeStaticMWrapperFunction(const string &basename, const string &ending) const
{
  string name;
  if (ending == "g1")
    name = "static_resid_g1";
  else if (ending == "g2")
    name = "static_resid_g1_g2";
  else if (ending == "g3")
    name = "static_resid_g1_g2_g3";

  string filename = packageDir(basename) + "/" + name + ".m";
  ofstream output{filename, ios::out | ios::binary};
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
StaticModel::writeStaticMFileHelper(const string &basename,
                                    const string &name, const string &retvalname,
                                    const string &name_tt, size_t ttlen,
                                    const string &previous_tt_name,
                                    const ostringstream &init_s, const ostringstream &end_s,
                                    const ostringstream &s, const ostringstream &s_tt) const
{
  string filename = packageDir(basename) + "/" + name_tt + ".m";
  ofstream output{filename, ios::out | ios::binary};
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
StaticModel::writeStaticMCompatFile(const string &basename) const
{
  string filename = packageDir(basename) + "/static.m";
  ofstream output{filename, ios::out | ios::binary};
  if (!output.is_open())
    {
      cerr << "Error: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }
  int ntt = temporary_terms_mlv.size() + temporary_terms_derivatives[0].size() + temporary_terms_derivatives[1].size() + temporary_terms_derivatives[2].size() + temporary_terms_derivatives[3].size();

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
StaticModel::writeStaticCFile(const string &basename) const
{
  // Writing comments and function definition command
  string filename{basename + "/model/src/static.c"};

  int ntt{static_cast<int>(temporary_terms_mlv.size() + temporary_terms_derivatives[0].size() + temporary_terms_derivatives[1].size() + temporary_terms_derivatives[2].size() + temporary_terms_derivatives[3].size())};

  ofstream output{filename, ios::out | ios::binary};
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
         << endl
         << "#include <math.h>" << endl
         << "#include <stdlib.h>" << endl
         << R"(#include "mex.h")" << endl
         << endl;

  // Write function definition if BinaryOpcode::powerDeriv is used
  writePowerDeriv(output);

  output << endl;

  auto [d_output, tt_output] = writeModelFileHelper<ExprNodeOutputType::CStaticModel>();

  for (size_t i = 0; i < d_output.size(); i++)
    {
      string funcname{i == 0 ? "resid" : "g" + to_string(i)};
      output << "void static_" << funcname << "_tt(const double *restrict y, const double *restrict x, const double *restrict params, double *restrict T)" << endl
             << "{" << endl
             << tt_output[i].str()
             << "}" << endl
             << endl
             << "void static_" << funcname << "(const double *restrict y, const double *restrict x, const double *restrict params, const double *restrict T, ";
      if (i == 0)
        output << "double *restrict residual";
      else if (i == 1)
        output << "double *restrict g1";
      else
        output << "double *restrict " << funcname << "_i, double *restrict " << funcname << "_j, double *restrict " << funcname << "_v";
      output << ")" << endl
                   << "{" << endl;
      if (i == 0)
        output << "  double lhs, rhs;" << endl;
      output << d_output[i].str()
             << "}" << endl
             << endl;
    }

  output << "void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])" << endl
         << "{" << endl
         << "  if (nrhs > 3)" << endl
         << R"(    mexErrMsgTxt("Accepts at most 3 output arguments");)" << endl
         << "  if (nrhs != 3)" << endl
         << R"(    mexErrMsgTxt("Requires exactly 3 input arguments");)" << endl
         << "  double *y = mxGetPr(prhs[0]);" << endl
         << "  double *x = mxGetPr(prhs[1]);" << endl
         << "  double *params = mxGetPr(prhs[2]);" << endl
         << endl
         << "  double *T = (double *) malloc(sizeof(double)*" << ntt << ");" << endl
         << endl
         << "  if (nlhs >= 1)" << endl
         << "    {" << endl
         << "      plhs[0] = mxCreateDoubleMatrix(" << equations.size() << ",1, mxREAL);" << endl
         << "      double *residual = mxGetPr(plhs[0]);" << endl
         << "      static_resid_tt(y, x, params, T);" << endl
         << "      static_resid(y, x, params, T, residual);" << endl
         << "    }" << endl
         << endl
         << "  if (nlhs >= 2)" << endl
         << "    {" << endl
         << "      plhs[1] = mxCreateDoubleMatrix(" << equations.size() << ", " << symbol_table.endo_nbr() << ", mxREAL);" << endl
         << "      double *g1 = mxGetPr(plhs[1]);" << endl
         << "      static_g1_tt(y, x, params, T);" << endl
         << "      static_g1(y, x, params, T, g1);" << endl
         << "    }" << endl
         << endl
         << "  if (nlhs >= 3)" << endl
         << "    {" << endl
         << "      mxArray *g2_i = mxCreateDoubleMatrix(" << NNZDerivatives[2] << ", " << 1 << ", mxREAL);" << endl
         << "      mxArray *g2_j = mxCreateDoubleMatrix(" << NNZDerivatives[2] << ", " << 1 << ", mxREAL);" << endl
         << "      mxArray *g2_v = mxCreateDoubleMatrix(" << NNZDerivatives[2] << ", " << 1 << ", mxREAL);" << endl
         << "      static_g2_tt(y, x, params, T);" << endl
         << "      static_g2(y, x, params, T, mxGetPr(g2_i), mxGetPr(g2_j), mxGetPr(g2_v));" << endl
         << "      mxArray *m = mxCreateDoubleScalar(" << equations.size() << ");" << endl
         << "      mxArray *n = mxCreateDoubleScalar(" << symbol_table.endo_nbr()*symbol_table.endo_nbr() << ");" << endl
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
         << "  free(T);" << endl
         << "}" << endl;

  output.close();
}

void
StaticModel::writeStaticJuliaFile(const string &basename) const
{
  auto [d_output, tt_output] = writeModelFileHelper<ExprNodeOutputType::juliaStaticModel>();

  stringstream output;
  output << "module " << basename << "Static" << endl
         << "#" << endl
         << "# NB: this file was automatically generated by Dynare" << endl
         << "#     from " << basename << ".mod" << endl
         << "#" << endl
         << "using StatsFuns" << endl << endl
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
         << "  T        : Vector{<: Real}(num_temp_terms) temporary terms" << endl
         << "  y        : Vector{<: Real}(model_.endo_nbr) endogenous variables in declaration order" << endl
         << "  x        : Vector{<: Real}(model_.exo_nbr) exogenous variables in declaration order" << endl
         << "  params   : Vector{<: Real}(model_.param) parameter values in declaration order" << endl
         << "  residual : Vector{<: Real}(model_.eq_nbr) residuals of the static model equations" << endl
         << "             in order of declaration of the equations. Dynare may prepend auxiliary equations," << endl
         << "             see model.aux_vars" << endl
         << "  g1       : Matrix{<: Real}(model.eq_nbr,model_.endo_nbr) Jacobian matrix of the static model equations" << endl
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
         << "tmp_nbr[1] = " << temporary_terms_mlv.size() + temporary_terms_derivatives[0].size() << "# Number of temporary terms for the residuals" << endl
         << "tmp_nbr[2] = " << temporary_terms_derivatives[1].size() << "# Number of temporary terms for g1 (jacobian)" << endl
         << "tmp_nbr[3] = " << temporary_terms_derivatives[2].size() << "# Number of temporary terms for g2 (hessian)" << endl
         << "tmp_nbr[4] = " << temporary_terms_derivatives[3].size() << "# Number of temporary terms for g3 (third order derivates)" << endl << endl;

  // staticResidTT!
  output << "function staticResidTT!(T::Vector{<: Real}," << endl
         << "                        y::Vector{<: Real}, x::Vector{<: Real}, params::Vector{<: Real})" << endl
         << "    @assert length(T) >= " << temporary_terms_mlv.size() + temporary_terms_derivatives[0].size()  << endl
         << "    @inbounds begin" << endl
         << tt_output[0].str()
	 << "    end" << endl
         << "    return nothing" << endl
         << "end" << endl << endl;

  // static!
  output << "function staticResid!(T::Vector{<: Real}, residual::Vector{<: Real}," << endl
         << "                      y::Vector{<: Real}, x::Vector{<: Real}, params::Vector{<: Real}, T0_flag::Bool)" << endl
         << "    @assert length(y) == " << symbol_table.endo_nbr() << endl
         << "    @assert length(x) == " << symbol_table.exo_nbr() << endl
         << "    @assert length(params) == " << symbol_table.param_nbr() << endl
         << "    @assert length(residual) == " << equations.size() << endl
         << "    if T0_flag" << endl
         << "        staticResidTT!(T, y, x, params)" << endl
         << "    end" << endl
         << "    @inbounds begin" << endl
         << d_output[0].str()
	 << "    end" << endl
         << "    if ~isreal(residual)" << endl
         << "        residual = real(residual)+imag(residual).^2;" << endl
         << "    end" << endl
         << "    return nothing" << endl
         << "end" << endl << endl;

  // staticG1TT!
  output << "function staticG1TT!(T::Vector{<: Real}," << endl
         << "                     y::Vector{<: Real}, x::Vector{<: Real}, params::Vector{<: Real}, T0_flag::Bool)" << endl
         << "    if T0_flag" << endl
         << "        staticResidTT!(T, y, x, params)" << endl
         << "    end" << endl
         << "    @inbounds begin" << endl
         << tt_output[1].str()
	 << "    end" << endl
         << "    return nothing" << endl
         << "end" << endl << endl;

  // staticG1!
  output << "function staticG1!(T::Vector{<: Real}, g1::Matrix{<: Real}," << endl
         << "                   y::Vector{<: Real}, x::Vector{<: Real}, params::Vector{<: Real}, T1_flag::Bool, T0_flag::Bool)" << endl
         << "    @assert length(T) >= "
         << temporary_terms_mlv.size() + temporary_terms_derivatives[0].size() + temporary_terms_derivatives[1].size() << endl
         << "    @assert size(g1) == (" << equations.size() << ", " << symbol_table.endo_nbr() << ")" << endl
         << "    @assert length(y) == " << symbol_table.endo_nbr() << endl
         << "    @assert length(x) == " << symbol_table.exo_nbr() << endl
         << "    @assert length(params) == " << symbol_table.param_nbr() << endl
         << "    if T1_flag" << endl
         << "        staticG1TT!(T, y, x, params, T0_flag)" << endl
         << "    end" << endl
         << "    fill!(g1, 0.0)" << endl
         << "     @inbounds begin" << endl
         << d_output[1].str()
	 << "    end" << endl
         << "    if ~isreal(g1)" << endl
         << "        g1 = real(g1)+2*imag(g1);" << endl
         << "    end" << endl
         << "    return nothing" << endl
         << "end" << endl << endl;

  // staticG2TT!
  output << "function staticG2TT!(T::Vector{<: Real}," << endl
         << "                     y::Vector{<: Real}, x::Vector{<: Real}, params::Vector{<: Real}, T1_flag::Bool, T0_flag::Bool)" << endl
         << "    if T1_flag" << endl
         << "        staticG1TT!(T, y, x, params, TO_flag)" << endl
         << "    end" << endl
         << "    @inbounds begin" << endl
         << tt_output[2].str()
	 << "    end" << endl
         << "    return nothing" << endl
         << "end" << endl << endl;

  // staticG2!
  int hessianColsNbr{symbol_table.endo_nbr() * symbol_table.endo_nbr()};
  output << "function staticG2!(T::Vector{<: Real}, g2::Matrix{<: Real}," << endl
         << "                   y::Vector{<: Real}, x::Vector{<: Real}, params::Vector{<: Real}, T2_flag::Bool, T1_flag::Bool, T0_flag::Bool)" << endl
         << "    @assert length(T) >= "
         << temporary_terms_mlv.size() + temporary_terms_derivatives[0].size() + temporary_terms_derivatives[1].size() + temporary_terms_derivatives[2].size() << endl
         << "    @assert size(g2) == (" << equations.size() << ", " << hessianColsNbr << ")" << endl
         << "    @assert length(y) == " << symbol_table.endo_nbr() << endl
         << "    @assert length(x) == " << symbol_table.exo_nbr() << endl
         << "    @assert length(params) == " << symbol_table.param_nbr() << endl
         << "    if T2_flag" << endl
         << "        staticG2TT!(T, y, x, params, T1_flag, T0_flag)" << endl
         << "    end" << endl
         << "    fill!(g2, 0.0)" << endl
         << "    @inbounds begin" << endl
         << d_output[2].str()
	 << "    end" << endl
         << "    return nothing" << endl
         << "end" << endl << endl;

  // staticG3TT!
  output << "function staticG3TT!(T::Vector{<: Real}," << endl
         << "                     y::Vector{<: Real}, x::Vector{<: Real}, params::Vector{<: Real}, T2_flag::Bool, T1_flag::Bool, T0_flag::Bool)" << endl
         << "    if T2_flag" << endl
         << "        staticG2TT!(T, y, x, params, T1_flag, T0_flag)" << endl
         << "    end" << endl
         << "    @inbounds begin" << endl
         << tt_output[3].str()
	 << "    end" << endl
         << "    return nothing" << endl
         << "end" << endl << endl;

  // staticG3!
  int ncols{hessianColsNbr * symbol_table.endo_nbr()};
  output << "function staticG3!(T::Vector{<: Real}, g3::Matrix{<: Real}," << endl
         << "                   y::Vector{<: Real}, x::Vector{<: Real}, params::Vector{<: Real}, T3_flag::Bool, T2_flag::Bool, T1_flag::Bool, T0_flag::Bool)" << endl
         << "    @assert length(T) >= "
         << temporary_terms_mlv.size() + temporary_terms_derivatives[0].size() + temporary_terms_derivatives[1].size() + temporary_terms_derivatives[2].size() + temporary_terms_derivatives[3].size() << endl
         << "    @assert size(g3) == (" << equations.size() << ", " << ncols << ")" << endl
         << "    @assert length(y) == " << symbol_table.endo_nbr() << endl
         << "    @assert length(x) == " << symbol_table.exo_nbr() << endl
         << "    @assert length(params) == " << symbol_table.param_nbr() << endl
         << "    if T3_flag" << endl
         << "        staticG3TT!(T, y, x, params, T2_flag, T1_flag, T0_flag)" << endl
         << "    end" << endl
         << "    fill!(g3, 0.0)" << endl
         << "    @inbounds begin" << endl
         << d_output[3].str()
	 << "    end" << endl
         << "    return nothing" << endl
         << "end" << endl << endl;

  // static!
  output << "function static!(T::Vector{<: Real}, residual::Vector{<: Real}," << endl
         << "                  y::Vector{<: Real}, x::Vector{<: Real}, params::Vector{<: Real})" << endl
         << "    staticResid!(T, residual, y, x, params, true)" << endl
         << "    return nothing" << endl
         << "end" << endl
         << endl
         << "function static!(T::Vector{<: Real}, residual::Vector{<: Real}, g1::Matrix{<: Real}," << endl
         << "                 y::Vector{<: Real}, x::Vector{<: Real}, params::Vector{<: Real})" << endl
         << "    staticG1!(T, g1, y, x, params, true, true)" << endl
         << "    staticResid!(T, residual, y, x, params, false)" << endl
         << "    return nothing" << endl
         << "end" << endl
         << endl
         << "function static!(T::Vector{<: Real}, g1::Matrix{<: Real}," << endl
         << "                 y::Vector{<: Real}, x::Vector{<: Real}, params::Vector{<: Real})" << endl
         << "    staticG1!(T, g1, y, x, params, true, false)" << endl
         << "    return nothing" << endl
         << "end" << endl
         << endl
         << "function static!(T::Vector{<: Real}, residual::Vector{<: Real}, g1::Matrix{<: Real}, g2::Matrix{<: Real}," << endl
         << "                 y::Vector{<: Real}, x::Vector{<: Real}, params::Vector{<: Real})" << endl
         << "    staticG2!(T, g2, y, x, params, true, true, true)" << endl
         << "    staticG1!(T, g1, y, x, params, false, false)" << endl
         << "    staticResid!(T, residual, y, x, params, false)" << endl
         << "    return nothing" << endl
         << "end" << endl
         << endl;

  // Write function definition if BinaryOpcode::powerDeriv is used
  writePowerDerivJulia(output);

  output << "end" << endl;

  writeToFileIfModified(output, basename + "Static.jl");
}

void
StaticModel::writeStaticFile(const string &basename, bool block, bool use_dll, const string &mexext, const filesystem::path &matlabroot, const filesystem::path &dynareroot, bool julia) const
{
  filesystem::path model_dir{basename};
  model_dir /= "model";
  if (use_dll)
    filesystem::create_directories(model_dir / "src");
  filesystem::create_directories(model_dir / "bytecode");

  if (block)
    {
      writeStaticBlockBytecode(basename);

      if (use_dll)
        {
          writeStaticPerBlockCFiles(basename);
          writeStaticBlockCFile(basename);
          vector<filesystem::path> src_files(blocks.size() + 1);
          for (int blk = 0; blk < static_cast<int>(blocks.size()); blk++)
            src_files[blk] = model_dir / "src" / ("static_" + to_string(blk+1) + ".c");
          src_files[blocks.size()] = model_dir / "src" / "static.c";
          compileMEX(basename, "static", mexext, src_files, matlabroot, dynareroot);
        }
      else if (julia)
        {
          cerr << "'block' option is not available with Julia" << endl;
          exit(EXIT_FAILURE);
        }
      else // M-files
        {
          writeStaticPerBlockMFiles(basename);
          writeStaticBlockMFile(basename);
        }
    }
  else
    {
      writeStaticBytecode(basename);

      if (use_dll)
        {
          writeStaticCFile(basename);
          compileMEX(basename, "static", mexext, { model_dir / "src" / "static.c" },
                     matlabroot, dynareroot);
        }
      else if (julia)
        writeStaticJuliaFile(basename);
      else // M-files
        writeStaticMFile(basename);
    }

  writeSetAuxiliaryVariables(basename, julia);
}

bool
StaticModel::exoPresentInEqs() const
{
  for (auto equation : equations)
    if (equation->hasExogenous())
      return true;
  return false;
}

void
StaticModel::writeStaticBlockMFile(const string &basename) const
{
  string filename = packageDir(basename) + "/static.m";

  ofstream output{filename, ios::out | ios::binary};
  if (!output.is_open())
    {
      cerr << "ERROR: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }

  output << "function [residual, y, T, g1] = static(nblock, y, x, params, T)" << endl
         << "  switch nblock" << endl;

  for (int blk = 0; blk < static_cast<int>(blocks.size()); blk++)
    {
      output << "    case " << blk+1 << endl;

      BlockSimulationType simulation_type = blocks[blk].simulation_type;

      if (simulation_type == BlockSimulationType::evaluateBackward
          || simulation_type == BlockSimulationType::evaluateForward)
        output << "      [y, T] = " << basename << ".block.static_" << blk+1 << "(y, x, params, T);" << endl
               << "      residual = [];" << endl
               << "      g1 = [];" << endl;
      else
        output << "      [residual, y, T, g1] = " << basename << ".block.static_" << blk+1 << "(y, x, params, T);" << endl;

    }
  output << "  end" << endl
         << "end" << endl;
  output.close();
}

void
StaticModel::writeStaticBlockCFile(const string &basename) const
{
  string filename = basename + "/model/src/static.c";

  ofstream output{filename, ios::out | ios::binary};
  if (!output.is_open())
    {
      cerr << "ERROR: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }

  output << "#include <math.h>" << endl
         << R"(#include "mex.h")" << endl;

  for (int blk = 0; blk < static_cast<int>(blocks.size()); blk++)
    output << R"(#include "static_)" << blk+1 << R"(.h")" << endl;

  output << endl;
  writePowerDeriv(output);

  output << endl
         << "void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])" << endl
         << "{" << endl
         << "  if (nrhs != 5)" << endl
         << R"(    mexErrMsgTxt("Requires exactly 5 input arguments");)" << endl
         << "  if (nlhs > 4)" << endl
         << R"(    mexErrMsgTxt("Accepts at most 4 output arguments");)" << endl
         << "  int nblock = (int) mxGetScalar(prhs[0]);" << endl
         << "  const mxArray *y = prhs[1], *x = prhs[2], *params = prhs[3], *T = prhs[4];" << endl
         << "  mxArray *T_new = mxDuplicateArray(T);" << endl
         << "  mxArray *y_new = mxDuplicateArray(y);" << endl
         << "  mxArray *residual, *g1;" << endl
         << "  switch (nblock)" << endl
         << "    {" << endl;

  for (int blk = 0; blk < static_cast<int>(blocks.size()); blk++)
    {
      output << "    case " << blk+1 << ':' << endl;

      BlockSimulationType simulation_type = blocks[blk].simulation_type;

      if (simulation_type == BlockSimulationType::evaluateBackward
          || simulation_type == BlockSimulationType::evaluateForward)
        output << "      static_" << blk+1 << "_mx(y_new, x, params, T_new);" << endl
               << "      residual = mxCreateDoubleMatrix(0,0,mxREAL);" << endl
               << "      g1 = mxCreateDoubleMatrix(0,0,mxREAL);" << endl;
      else
        output << "      static_" << blk+1 << "_mx(y_new, x, params, T_new, &residual, &g1);" << endl;
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
         << "}" << endl;
  output.close();
}

void
StaticModel::writeDriverOutput(ostream &output, bool block) const
{
  output << "M_.static_tmp_nbr = [";
  for (const auto &temporary_terms_derivative : temporary_terms_derivatives)
    output << temporary_terms_derivative.size() << "; ";
  output << "];" << endl;

  /* Write mapping between model local variables and indices in the temporary
     terms vector (dynare#1722) */
  output << "M_.model_local_variables_static_tt_idxs = {" << endl;
  for (auto [mlv, value] : temporary_terms_mlv)
    output << "  '" << symbol_table.getName(mlv->symb_id) << "', "
           << temporary_terms_idxs.at(mlv)+1 << ';' << endl;
  output << "};" << endl;

  if (!block)
    return;

  for (int blk = 0; blk < static_cast<int>(blocks.size()); blk++)
    {
      output << "block_structure_stat.block(" << blk+1 << ").Simulation_Type = " << static_cast<int>(blocks[blk].simulation_type) << ";" << endl
             << "block_structure_stat.block(" << blk+1 << ").endo_nbr = " << blocks[blk].size << ";" << endl
             << "block_structure_stat.block(" << blk+1 << ").mfs = " << blocks[blk].mfs_size << ";" << endl
             << "block_structure_stat.block(" << blk+1 << ").equation = [";
      for (int eq = 0; eq < blocks[blk].size; eq++)
        output << " " << getBlockEquationID(blk, eq)+1;
      output << "];" << endl
             << "block_structure_stat.block(" << blk+1 << ").variable = [";
      for (int var = 0; var < blocks[blk].size; var++)
        output << " " << getBlockVariableID(blk, var)+1;
      output << "];" << endl;
    }
  output << "M_.block_structure_stat.block = block_structure_stat.block;" << endl
         << "M_.block_structure_stat.variable_reordered = [";
  for (int i = 0; i < symbol_table.endo_nbr(); i++)
    output << " " << endo_idx_block2orig[i]+1;
  output << "];" << endl
         << "M_.block_structure_stat.equation_reordered = [";
  for (int i = 0; i < symbol_table.endo_nbr(); i++)
    output << " " << eq_idx_block2orig[i]+1;
  output << "];" << endl;

  set<pair<int, int>> row_incidence;
  for (const auto &[indices, d1] : derivatives[1])
    if (int deriv_id = indices[1];
        getTypeByDerivID(deriv_id) == SymbolType::endogenous)
      {
        int eq = indices[0];
        int var { getTypeSpecificIDByDerivID(deriv_id) };
        row_incidence.emplace(eq, var);
      }
  output << "M_.block_structure_stat.incidence.sparse_IM = [" << endl;
  for (auto [eq, var] : row_incidence)
    output << " " << eq+1 << " " << var+1 << ";" << endl;
  output << "];" << endl
         << "M_.block_structure_stat.tmp_nbr = " << blocks_temporary_terms_idxs.size()
         << ";" << endl;
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
StaticModel::getLagByDerivID([[maybe_unused]] int deriv_id) const noexcept(false)
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
StaticModel::getTypeSpecificIDByDerivID(int deriv_id) const
{
  if (deriv_id < symbol_table.endo_nbr())
    return deriv_id;
  else if (deriv_id < symbol_table.endo_nbr() + symbol_table.param_nbr())
    return deriv_id - symbol_table.endo_nbr();
  else
    throw UnknownDerivIDException();
}

int
StaticModel::getDerivID(int symb_id, [[maybe_unused]] int lag) const noexcept(false)
{
  if (symbol_table.getType(symb_id) == SymbolType::endogenous)
    return symbol_table.getTypeSpecificID(symb_id);
  else if (symbol_table.getType(symb_id) == SymbolType::parameter)
    return symbol_table.getTypeSpecificID(symb_id) + symbol_table.endo_nbr();
  else
    /* See the special treatment in VariableNode::prepareForDerivation(),
       VariableNode::computeDerivative() and VariableNode::getChainRuleDerivative() */
    throw UnknownDerivIDException{};
}

void
StaticModel::addAllParamDerivId(set<int> &deriv_id_set)
{
  for (int i = 0; i < symbol_table.param_nbr(); i++)
    deriv_id_set.insert(i + symbol_table.endo_nbr());
}

void
StaticModel::computeChainRuleJacobian()
{
  int nb_blocks = blocks.size();
  blocks_derivatives.resize(nb_blocks);
  for (int blk = 0; blk < nb_blocks; blk++)
    {
      int nb_recursives = blocks[blk].getRecursiveSize();

      map<int, BinaryOpNode *> recursive_vars;
      for (int i = 0; i < nb_recursives; i++)
        {
          int deriv_id = getDerivID(symbol_table.getID(SymbolType::endogenous, getBlockVariableID(blk, i)), 0);
          if (getBlockEquationType(blk, i) == EquationType::evaluateRenormalized)
            recursive_vars[deriv_id] = getBlockEquationRenormalizedExpr(blk, i);
          else
            recursive_vars[deriv_id] = getBlockEquationExpr(blk, i);
        }

      assert(blocks[blk].simulation_type != BlockSimulationType::solveTwoBoundariesSimple
             && blocks[blk].simulation_type != BlockSimulationType::solveTwoBoundariesComplete);

      int size = blocks[blk].size;

      for (int eq = nb_recursives; eq < size; eq++)
        {
          int eq_orig = getBlockEquationID(blk, eq);
          for (int var = nb_recursives; var < size; var++)
            {
              int var_orig = getBlockVariableID(blk, var);
              expr_t d1 = equations[eq_orig]->getChainRuleDerivative(getDerivID(symbol_table.getID(SymbolType::endogenous, var_orig), 0), recursive_vars);
              if (d1 != Zero)
                blocks_derivatives[blk][{ eq, var, 0 }] = d1;
            }
        }
    }
}

void
StaticModel::writeLatexFile(const string &basename, bool write_equation_tags) const
{
  writeLatexModelFile(basename, "static", ExprNodeOutputType::latexStaticModel, write_equation_tags);
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
StaticModel::writeSetAuxiliaryVariables(const string &basename, bool julia) const
{
  ostringstream output_func_body;
  ExprNodeOutputType output_type = julia ? ExprNodeOutputType::juliaStaticModel : ExprNodeOutputType::matlabStaticModel;
  writeAuxVarRecursiveDefinitions(output_func_body, output_type);

  if (output_func_body.str().empty())
    return;

  string func_name = julia ? "set_auxiliary_variables!" : "set_auxiliary_variables";
  string comment = julia ? "#" : "%";

  stringstream output;
  if (julia)
    output << "module " << basename << "SetAuxiliaryVariables" << endl
           << "export " << func_name << endl;
  output << "function ";
  if (!julia)
    output << "y = ";
  output << func_name << "(y, x, params)" << endl
         << comment << endl
         << comment << " Status : Computes static model for Dynare" << endl
         << comment << endl
         << comment << " Warning : this file is generated automatically by Dynare" << endl
         << comment << "           from model file (.mod)" << endl << endl;
  if (julia)
    output << "@inbounds begin" << endl;
  output << output_func_body.str()
         << "end" << endl;
  if (julia)
    output << "end" << endl
           << "end" << endl;

  writeToFileIfModified(output, julia ? basename + "SetAuxiliaryVariables.jl" : packageDir(basename) + "/" + func_name + ".m");
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
      output << R"(\begin{dmath})" << endl;
      dynamic_cast<ExprNode *>(aux_equation->substituteStaticAuxiliaryDefinition())->writeOutput(output, ExprNodeOutputType::latexStaticModel);
      output << endl << R"(\end{dmath})" << endl;
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
        for (bool printed_something{false};
             const auto &it : efout)
          {
            if (exchange(printed_something, true))
              output << ", ";
            output << it;
          }
      }

  for (auto aux_equation : aux_equations)
    {
      output << R"(, {"lhs": ")";
      aux_equation->arg1->writeJsonOutput(output, temporary_terms, tef_terms, false);
      output << R"(", "rhs": ")";
      dynamic_cast<BinaryOpNode *>(aux_equation->substituteStaticAuxiliaryDefinition())->arg2->writeJsonOutput(output, temporary_terms, tef_terms, false);
      output << R"("})";
    }
}

void
StaticModel::writeJsonOutput(ostream &output) const
{
  deriv_node_temp_terms_t tef_terms;
  writeJsonModelLocalVariables(output, false, tef_terms);
  output << ", ";
  writeJsonModelEquations(output, false);
}

void
StaticModel::writeJsonComputingPassOutput(ostream &output, bool writeDetails) const
{
  ostringstream model_local_vars_output; // Used for storing model local vars
  vector<ostringstream> d_output(derivatives.size()); // Derivatives output (at all orders, including 0=residual)

  deriv_node_temp_terms_t tef_terms;
  temporary_terms_t temp_term_union;

  writeJsonModelLocalVariables(model_local_vars_output, true, tef_terms);

  writeJsonTemporaryTerms(temporary_terms_derivatives[0], temp_term_union, d_output[0], tef_terms, "");
  d_output[0] << ", ";
  writeJsonModelEquations(d_output[0], true);

  int ncols = symbol_table.endo_nbr();
  for (size_t i = 1; i < derivatives.size(); i++)
    {
      string matrix_name = i == 1 ? "jacobian" : i == 2 ? "hessian" : i == 3 ? "third_derivative" : to_string(i) + "th_derivative";
      writeJsonTemporaryTerms(temporary_terms_derivatives[i], temp_term_union, d_output[i], tef_terms, matrix_name);
      temp_term_union.insert(temporary_terms_derivatives[i].begin(), temporary_terms_derivatives[i].end());
      d_output[i] << R"(, ")" << matrix_name  << R"(": {)"
                  << R"(  "nrows": )" << equations.size()
                  << R"(, "ncols": )" << ncols
                  << R"(, "entries": [)";

      for (bool printed_something{false};
           const auto &[vidx, d] : derivatives[i])
        {
          if (exchange(printed_something, true))
            d_output[i] << ", ";

          int eq = vidx[0];

          int col_idx = 0;
          for (size_t j = 1; j < vidx.size(); j++)
            {
              col_idx *= symbol_table.endo_nbr();
              col_idx += getJacobianCol(vidx[j]);
            }

          if (writeDetails)
            d_output[i] << R"({"eq": )" << eq + 1;
          else
            d_output[i] << R"({"row": )" << eq + 1;

          d_output[i] << R"(, "col": )" << (i > 1 ? "[" : "") << col_idx + 1;

          if (i == 2 && vidx[1] != vidx[2]) // Symmetric elements in hessian
            {
              int col_idx_sym = getJacobianCol(vidx[2]) * symbol_table.endo_nbr() + getJacobianCol(vidx[1]);
              d_output[i] << ", " << col_idx_sym + 1;
            }
          if (i > 1)
            d_output[i] << "]";

          if (writeDetails)
            for (size_t j = 1; j < vidx.size(); j++)
              d_output[i] << R"(, "var)" << (i > 1 ? to_string(j) : "") << R"(": ")" << getNameByDerivID(vidx[j]) << R"(")";

          d_output[i] << R"(, "val": ")";
          d->writeJsonOutput(d_output[i], temp_term_union, tef_terms);
          d_output[i] << R"("})" << endl;
        }
      d_output[i] << "]}";

      ncols *= symbol_table.endo_nbr();
    }

  if (writeDetails)
    output << R"("static_model": {)";
  else
    output << R"("static_model_simple": {)";
  output << model_local_vars_output.str();
  for (const auto &it : d_output)
    output << ", " << it.str();
  output << "}";
}

void
StaticModel::writeJsonParamsDerivativesFile(ostream &output, bool writeDetails) const
{
  if (!params_derivatives.size())
    return;

  ostringstream model_local_vars_output; // Used for storing model local vars
  ostringstream model_output; // Used for storing model temp vars and equations
  ostringstream jacobian_output; // Used for storing jacobian equations
  ostringstream hessian_output; // Used for storing Hessian equations
  ostringstream hessian1_output; // Used for storing Hessian equations
  ostringstream third_derivs_output; // Used for storing third order derivatives equations
  ostringstream third_derivs1_output; // Used for storing third order derivatives equations

  deriv_node_temp_terms_t tef_terms;
  writeJsonModelLocalVariables(model_local_vars_output, true, tef_terms);

  temporary_terms_t temp_term_union;
  for (const auto &it : params_derivs_temporary_terms)
    writeJsonTemporaryTerms(it.second, temp_term_union, model_output, tef_terms, "all");

  jacobian_output << R"("deriv_wrt_params": {)"
                  << R"(  "neqs": )" << equations.size()
                  << R"(, "nparamcols": )" << symbol_table.param_nbr()
                  << R"(, "entries": [)";
  for (bool printed_something{false};
       const auto &[vidx, d] : params_derivatives.find({ 0, 1 })->second)
    {
      if (exchange(printed_something, true))
        jacobian_output << ", ";

      auto [eq, param] = vectorToTuple<2>(vidx);

      int param_col { getTypeSpecificIDByDerivID(param) + 1 };

      if (writeDetails)
        jacobian_output << R"({"eq": )" << eq + 1;
      else
        jacobian_output << R"({"row": )" << eq + 1;

      if (writeDetails)
        jacobian_output << R"(, "param_col": )" << param_col;

      jacobian_output << R"(, "param": ")" << getNameByDerivID(param) << R"(")";

      jacobian_output << R"(, "val": ")";
      d->writeJsonOutput(jacobian_output, temp_term_union, tef_terms);
      jacobian_output << R"("})" << endl;
    }
  jacobian_output << "]}";

  hessian_output << R"("deriv_jacobian_wrt_params": {)"
                 << R"(  "neqs": )" << equations.size()
                 << R"(, "nvarcols": )" << symbol_table.endo_nbr()
                 << R"(, "nparamcols": )" << symbol_table.param_nbr()
                 << R"(, "entries": [)";
  for (bool printed_something{false};
       const auto &[vidx, d] : params_derivatives.find({ 1, 1 })->second)
    {
      if (exchange(printed_something, true))
        hessian_output << ", ";

      auto [eq, var, param] = vectorToTuple<3>(vidx);

      int var_col { getTypeSpecificIDByDerivID(var) + 1 };
      int param_col { getTypeSpecificIDByDerivID(param) + 1 };

      if (writeDetails)
        hessian_output << R"({"eq": )" << eq + 1;
      else
        hessian_output << R"({"row": )" << eq + 1;

      if (writeDetails)
        hessian_output << R"(, "var": ")" << getNameByDerivID(var) << R"(")"
                       << R"(, "param": ")" << getNameByDerivID(param) << R"(")";

      hessian_output << R"(, "var_col": )" << var_col
                     << R"(, "param_col": )" << param_col
                     << R"(, "val": ")";
      d->writeJsonOutput(hessian_output, temp_term_union, tef_terms);
      hessian_output << R"("})" << endl;
    }
  hessian_output << "]}";

  hessian1_output << R"("second_deriv_residuals_wrt_params": {)"
                  << R"(  "nrows": )" << equations.size()
                  << R"(, "nparam1cols": )" << symbol_table.param_nbr()
                  << R"(, "nparam2cols": )" << symbol_table.param_nbr()
                  << R"(, "entries": [)";
  for (bool printed_something{false};
       const auto &[vidx, d] : params_derivatives.find({ 0, 2 })->second)
    {
      if (exchange(printed_something, true))
        hessian1_output << ", ";

      auto [eq, param1, param2] = vectorToTuple<3>(vidx);

      int param1_col { getTypeSpecificIDByDerivID(param1) + 1 };
      int param2_col { getTypeSpecificIDByDerivID(param2) + 1 };

      if (writeDetails)
        hessian1_output << R"({"eq": )" << eq + 1;
      else
        hessian1_output << R"({"row": )" << eq + 1;

      hessian1_output << R"(, "param1_col": )" << param1_col
                      << R"(, "param2_col": )" << param2_col;

      if (writeDetails)
        hessian1_output << R"(, "param1": ")" << getNameByDerivID(param1) << R"(")"
                        << R"(, "param2": ")" << getNameByDerivID(param2) << R"(")";

      hessian1_output << R"(, "val": ")";
      d->writeJsonOutput(hessian1_output, temp_term_union, tef_terms);
      hessian1_output << R"("})" << endl;
    }
  hessian1_output << "]}";

  third_derivs_output << R"("second_deriv_jacobian_wrt_params": {)"
                      << R"(  "neqs": )" << equations.size()
                      << R"(, "nvarcols": )" << symbol_table.endo_nbr()
                      << R"(, "nparam1cols": )" << symbol_table.param_nbr()
                      << R"(, "nparam2cols": )" << symbol_table.param_nbr()
                      << R"(, "entries": [)";
  for (bool printed_something{false};
       const auto &[vidx, d] : params_derivatives.find({ 1, 2 })->second)
    {
      if (exchange(printed_something, true))
        third_derivs_output << ", ";

      auto [eq, var, param1, param2] = vectorToTuple<4>(vidx);

      int var_col { getTypeSpecificIDByDerivID(var) + 1 };
      int param1_col { getTypeSpecificIDByDerivID(param1) + 1 };
      int param2_col { getTypeSpecificIDByDerivID(param2) + 1 };

      if (writeDetails)
        third_derivs_output << R"({"eq": )" << eq + 1;
      else
        third_derivs_output << R"({"row": )" << eq + 1;
      third_derivs_output << R"(, "var_col": )" << var_col
                          << R"(, "param1_col": )" << param1_col
                          << R"(, "param2_col": )" << param2_col;

      if (writeDetails)
        third_derivs_output << R"(, "var": ")" << getNameByDerivID(var) << R"(")"
                            << R"(, "param1": ")" << getNameByDerivID(param1) << R"(")"
                            << R"(, "param2": ")" << getNameByDerivID(param2) << R"(")";

      third_derivs_output << R"(, "val": ")";
      d->writeJsonOutput(third_derivs_output, temp_term_union, tef_terms);
      third_derivs_output << R"("})" << endl;
    }
  third_derivs_output << "]}" << endl;

  third_derivs1_output << R"("derivative_hessian_wrt_params": {)"
                       << R"(  "neqs": )" << equations.size()
                       << R"(, "nvar1cols": )" << symbol_table.endo_nbr()
                       << R"(, "nvar2cols": )" << symbol_table.endo_nbr()
                       << R"(, "nparamcols": )" << symbol_table.param_nbr()
                       << R"(, "entries": [)";
  for (bool printed_something{false};
       const auto &[vidx, d] : params_derivatives.find({ 2, 1 })->second)
    {
      if (exchange(printed_something, true))
        third_derivs1_output << ", ";

      auto [eq, var1, var2, param] = vectorToTuple<4>(vidx);

      int var1_col { getTypeSpecificIDByDerivID(var1) + 1 };
      int var2_col { getTypeSpecificIDByDerivID(var2) + 1 };
      int param_col { getTypeSpecificIDByDerivID(param) + 1 };

      if (writeDetails)
        third_derivs1_output << R"({"eq": )" << eq + 1;
      else
        third_derivs1_output << R"({"row": )" << eq + 1;

      third_derivs1_output << R"(, "var1_col": )" << var1_col
                           << R"(, "var2_col": )" << var2_col
                           << R"(, "param_col": )" << param_col;

      if (writeDetails)
        third_derivs1_output << R"(, "var1": ")" << getNameByDerivID(var1) << R"(")"
                             << R"(, "var2": ")" << getNameByDerivID(var2) << R"(")"
                             << R"(, "param1": ")" << getNameByDerivID(param) << R"(")";

      third_derivs1_output << R"(, "val": ")";
      d->writeJsonOutput(third_derivs1_output, temp_term_union, tef_terms);
      third_derivs1_output << R"("})" << endl;
    }
  third_derivs1_output << "]}" << endl;

  if (writeDetails)
    output << R"("static_model_params_derivative": {)";
  else
    output << R"("static_model_params_derivatives_simple": {)";
  output << model_local_vars_output.str()
         << ", " << model_output.str()
         << ", " << jacobian_output.str()
         << ", " << hessian_output.str()
         << ", " << hessian1_output.str()
         << ", " << third_derivs_output.str()
         << ", " << third_derivs1_output.str()
         << "}";
}
