/*
 * Copyright © 2003-2021 Dynare Team
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

void
StaticModel::copyHelper(const StaticModel &m)
{
}

StaticModel::StaticModel(SymbolTable &symbol_table_arg,
                         NumericalConstants &num_constants_arg,
                         ExternalFunctionsTable &external_functions_table_arg) :
  ModelTree{symbol_table_arg, num_constants_arg, external_functions_table_arg}
{
}

StaticModel::StaticModel(const StaticModel &m) :
  ModelTree{m}
{
  copyHelper(m);
}

StaticModel &
StaticModel::operator=(const StaticModel &m)
{
  ModelTree::operator=(m);

  copyHelper(m);

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
        if (dynamic_equations.find(i) != dynamic_equations.end())
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
        cerr << "...division by zero error encountred when converting equation " << i << " to static" << endl;
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
StaticModel::compileDerivative(ofstream &code_file, unsigned int &instruction_number, int eq, int symb_id, const temporary_terms_t &temporary_terms, const temporary_terms_idxs_t &temporary_terms_idxs) const
{
  if (auto it = derivatives[1].find({ eq, getDerivID(symbol_table.getID(SymbolType::endogenous, symb_id), 0) });
      it != derivatives[1].end())
    it->second->compile(code_file, instruction_number, false, temporary_terms, temporary_terms_idxs, false, false);
  else
    {
      FLDZ_ fldz;
      fldz.write(code_file, instruction_number);
    }
}

void
StaticModel::compileChainRuleDerivative(ofstream &code_file, unsigned int &instruction_number, int blk, int eq, int var, int lag, const temporary_terms_t &temporary_terms, const temporary_terms_idxs_t &temporary_terms_idxs) const
{
  if (auto it = blocks_derivatives[blk].find({ eq, var, lag });
      it != blocks_derivatives[blk].end())
    it->second->compile(code_file, instruction_number, false, temporary_terms, temporary_terms_idxs, false, false);
  else
    {
      FLDZ_ fldz;
      fldz.write(code_file, instruction_number);
    }
}

void
StaticModel::writeStaticPerBlockHelper(int blk, ostream &output, ExprNodeOutputType output_type, temporary_terms_t &temporary_terms) const
{
  BlockSimulationType simulation_type = blocks[blk].simulation_type;
  int block_recursive_size = blocks[blk].getRecursiveSize();

  // The equations
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

  for (int eq = 0; eq < blocks[blk].size; eq++)
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
          if (eq < block_recursive_size)
            goto evaluation;
          output << "  residual" << LEFT_ARRAY_SUBSCRIPT(output_type)
                 << eq-block_recursive_size+ARRAY_SUBSCRIPT_OFFSET(output_type)
                 << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=(";
          lhs->writeOutput(output, output_type, temporary_terms, blocks_temporary_terms_idxs);
          output << ")-(";
          rhs->writeOutput(output, output_type, temporary_terms, blocks_temporary_terms_idxs);
          output << ");" << endl;
          break;
        default:
          cerr << "Incorrect type for block " << blk+1 << endl;
          exit(EXIT_FAILURE);
        }
    }
  // The Jacobian if we have to solve the block
  if (simulation_type != BlockSimulationType::evaluateBackward
      && simulation_type != BlockSimulationType::evaluateForward)
    {
      // Write temporary terms for derivatives
      write_eq_tt(blocks[blk].size);

      ostringstream i_output, j_output, v_output;
      int line_counter = ARRAY_SUBSCRIPT_OFFSET(output_type);
      for (const auto &[indices, d] : blocks_derivatives[blk])
        {
          auto [eq, var, ignore] = indices;
          i_output << "  g1_i" << LEFT_ARRAY_SUBSCRIPT(output_type) << line_counter
                   << RIGHT_ARRAY_SUBSCRIPT(output_type) << '=' << eq+1-block_recursive_size
                   << ';' << endl;
          j_output << "  g1_j" << LEFT_ARRAY_SUBSCRIPT(output_type) << line_counter
                   << RIGHT_ARRAY_SUBSCRIPT(output_type) << '=' << var+1-block_recursive_size
                   << ';' << endl;
          v_output << "  g1_v" << LEFT_ARRAY_SUBSCRIPT(output_type) << line_counter
                   << RIGHT_ARRAY_SUBSCRIPT(output_type) << '=';
          d->writeOutput(v_output, output_type, temporary_terms, blocks_temporary_terms_idxs);
          v_output << ';' << endl;
          line_counter++;
        }
      output << i_output.str() << j_output.str() << v_output.str();
    }
}

void
StaticModel::writeStaticPerBlockMFiles(const string &basename) const
{
  temporary_terms_t temporary_terms; // Temp terms written so far

  for (int blk = 0; blk < static_cast<int>(blocks.size()); blk++)
    {
      BlockSimulationType simulation_type = blocks[blk].simulation_type;

      string filename = packageDir(basename + ".block") + "/static_" + to_string(blk+1) + ".m";
      ofstream output;
      output.open(filename, ios::out | ios::binary);
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
             << "  % //" << string("                     Block ").substr(static_cast<int>(log10(blk + 1))) << blk+1
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

      writeStaticPerBlockHelper(blk, output, ExprNodeOutputType::matlabStaticModel, temporary_terms);

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

      writeStaticPerBlockHelper(blk, output, ExprNodeOutputType::CStaticModel, temporary_terms);

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
StaticModel::writeStaticBytecode(const string &basename) const
{
  ostringstream tmp_output;
  ofstream code_file;
  unsigned int instruction_number = 0;
  bool file_open = false;

  string main_name = basename + "/model/bytecode/static.cod";
  code_file.open(main_name, ios::out | ios::binary | ios::ate);
  if (!code_file.is_open())
    {
      cerr << R"(Error : Can't open file ")" << main_name << R"(" for writing)" << endl;
      exit(EXIT_FAILURE);
    }
  int count_u;
  int u_count_int = 0;

  writeBytecodeBinFile(basename + "/model/bytecode/static.bin", u_count_int, file_open, false);
  file_open = true;

  // Compute the union of temporary terms from residuals and 1st derivatives
  temporary_terms_t temporary_terms = temporary_terms_derivatives[0];
  copy(temporary_terms_derivatives[1].begin(), temporary_terms_derivatives[1].end(),
       inserter(temporary_terms, temporary_terms.end()));

  //Temporary variables declaration
  FDIMST_ fdimst(temporary_terms.size());
  fdimst.write(code_file, instruction_number);
  FBEGINBLOCK_ fbeginblock(symbol_table.endo_nbr(),
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
                           symbol_table.endo_nbr());
  fbeginblock.write(code_file, instruction_number);

  temporary_terms_t temporary_terms_union;
  compileTemporaryTerms(code_file, instruction_number, false, false, temporary_terms_union, temporary_terms_idxs);

  compileModelEquations(code_file, instruction_number, false, false, temporary_terms_union, temporary_terms_idxs);

  FENDEQU_ fendequ;
  fendequ.write(code_file, instruction_number);

  // Get the current code_file position and jump if eval = true
  streampos pos1 = code_file.tellp();
  FJMPIFEVAL_ fjmp_if_eval(0);
  fjmp_if_eval.write(code_file, instruction_number);
  int prev_instruction_number = instruction_number;

  vector<vector<pair<int, int>>> my_derivatives(symbol_table.endo_nbr());
  count_u = symbol_table.endo_nbr();
  for (const auto & [indices, d1] : derivatives[1])
    {
      int deriv_id = indices[1];
      if (getTypeByDerivID(deriv_id) == SymbolType::endogenous)
        {
          int eq = indices[0];
          int symb = getSymbIDByDerivID(deriv_id);
          int var = symbol_table.getTypeSpecificID(symb);
          FNUMEXPR_ fnumexpr(ExpressionType::FirstEndoDerivative, eq, var);
          fnumexpr.write(code_file, instruction_number);
          if (!my_derivatives[eq].size())
            my_derivatives[eq].clear();
          my_derivatives[eq].emplace_back(var, count_u);

          d1->compile(code_file, instruction_number, false, temporary_terms_union, temporary_terms_idxs, false, false);

          FSTPSU_ fstpsu(count_u);
          fstpsu.write(code_file, instruction_number);
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
              FLDSU_ fldsu(it->second);
              fldsu.write(code_file, instruction_number);
              FLDSV_ fldsv{static_cast<int>(SymbolType::endogenous), static_cast<unsigned int>(it->first)};
              fldsv.write(code_file, instruction_number);
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

  temporary_terms_t tt2, tt3;

  // The Jacobian if we have to solve the block determinsitic bloc
  for (const auto & [indices, d1] : derivatives[1])
    {
      int deriv_id = indices[1];
      if (getTypeByDerivID(deriv_id) == SymbolType::endogenous)
        {
          int eq = indices[0];
          int symb = getSymbIDByDerivID(deriv_id);
          int var = symbol_table.getTypeSpecificID(symb);
          FNUMEXPR_ fnumexpr(ExpressionType::FirstEndoDerivative, eq, var);
          fnumexpr.write(code_file, instruction_number);
          if (!my_derivatives[eq].size())
            my_derivatives[eq].clear();
          my_derivatives[eq].emplace_back(var, count_u);

          d1->compile(code_file, instruction_number, false, temporary_terms_union, temporary_terms_idxs, false, false);
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
  ofstream code_file;
  unsigned int instruction_number = 0;
  expr_t lhs = nullptr, rhs = nullptr;
  BinaryOpNode *eq_node;
  Uff Uf[symbol_table.endo_nbr()];
  map<expr_t, int> reference_count;
  vector<int> feedback_variables;
  bool file_open = false;

  string main_name = basename + "/model/bytecode/static.cod";
  code_file.open(main_name, ios::out | ios::binary | ios::ate);
  if (!code_file.is_open())
    {
      cerr << R"(Error : Can't open file ")" << main_name << R"(" for writing)" << endl;
      exit(EXIT_FAILURE);
    }
  //Temporary variables declaration

  FDIMST_ fdimst(blocks_temporary_terms_idxs.size());
  fdimst.write(code_file, instruction_number);

  temporary_terms_t temporary_terms_union;

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

      if (simulation_type == BlockSimulationType::solveTwoBoundariesSimple
          || simulation_type == BlockSimulationType::solveTwoBoundariesComplete
          || simulation_type == BlockSimulationType::solveBackwardComplete
          || simulation_type == BlockSimulationType::solveForwardComplete)
        {
          writeBlockBytecodeBinFile(basename, block, u_count_int, file_open);
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
                               0,
                               0,
                               u_count_int,
                               block_size);

      fbeginblock.write(code_file, instruction_number);

      // Get the current code_file position and jump if eval = true
      streampos pos1 = code_file.tellp();
      FJMPIFEVAL_ fjmp_if_eval(0);
      fjmp_if_eval.write(code_file, instruction_number);
      int prev_instruction_number = instruction_number;

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
                                 it->compileExternalFunctionOutput(code_file, instruction_number, false, temporary_terms_union, blocks_temporary_terms_idxs, false, false, tef_terms);

                               FNUMEXPR_ fnumexpr(ExpressionType::TemporaryTerm, static_cast<int>(blocks_temporary_terms_idxs.at(it)));
                               fnumexpr.write(code_file, instruction_number);
                               it->compile(code_file, instruction_number, false, temporary_terms_union, blocks_temporary_terms_idxs, false, false, tef_terms);
                               FSTPST_ fstpst(static_cast<int>(blocks_temporary_terms_idxs.at(it)));
                               fstpst.write(code_file, instruction_number);
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
              {
                FNUMEXPR_ fnumexpr(ExpressionType::ModelEquation, getBlockEquationID(block, i));
                fnumexpr.write(code_file, instruction_number);
              }
              if (equ_type == EquationType::evaluate)
                {
                  eq_node = getBlockEquationExpr(block, i);
                  lhs = eq_node->arg1;
                  rhs = eq_node->arg2;
                  rhs->compile(code_file, instruction_number, false, temporary_terms_union, blocks_temporary_terms_idxs, false, false);
                  lhs->compile(code_file, instruction_number, true, temporary_terms_union, blocks_temporary_terms_idxs, false, false);
                }
              else if (equ_type == EquationType::evaluateRenormalized)
                {
                  eq_node = getBlockEquationRenormalizedExpr(block, i);
                  lhs = eq_node->arg1;
                  rhs = eq_node->arg2;
                  rhs->compile(code_file, instruction_number, false, temporary_terms_union, blocks_temporary_terms_idxs, false, false);
                  lhs->compile(code_file, instruction_number, true, temporary_terms_union, blocks_temporary_terms_idxs, false, false);
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
              FNUMEXPR_ fnumexpr(ExpressionType::ModelEquation, getBlockEquationID(block, i));
              fnumexpr.write(code_file, instruction_number);
              eq_node = getBlockEquationExpr(block, i);
              lhs = eq_node->arg1;
              rhs = eq_node->arg2;
              lhs->compile(code_file, instruction_number, false, temporary_terms_union, blocks_temporary_terms_idxs, false, false);
              rhs->compile(code_file, instruction_number, false, temporary_terms_union, blocks_temporary_terms_idxs, false, false);

              FBINARY_ fbinary{static_cast<int>(BinaryOpcode::minus)};
              fbinary.write(code_file, instruction_number);

              FSTPR_ fstpr(i - block_recursive);
              fstpr.write(code_file, instruction_number);
            }
        }
      FENDEQU_ fendequ;
      fendequ.write(code_file, instruction_number);

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
              {
                FNUMEXPR_ fnumexpr(ExpressionType::FirstEndoDerivative, 0, 0);
                fnumexpr.write(code_file, instruction_number);
              }
              compileDerivative(code_file, instruction_number, getBlockEquationID(block, 0), getBlockVariableID(block, 0), temporary_terms_union, blocks_temporary_terms_idxs);
              {
                FSTPG_ fstpg(0);
                fstpg.write(code_file, instruction_number);
              }
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
                      FNUMEXPR_ fnumexpr(ExpressionType::FirstEndoDerivative, eqr, varr);
                      fnumexpr.write(code_file, instruction_number);
                      compileChainRuleDerivative(code_file, instruction_number, block, eq, var, 0, temporary_terms_union, blocks_temporary_terms_idxs);
                      FSTPSU_ fstpsu(count_u);
                      fstpsu.write(code_file, instruction_number);
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
              {
                FNUMEXPR_ fnumexpr(ExpressionType::ModelEquation, getBlockEquationID(block, i));
                fnumexpr.write(code_file, instruction_number);
              }
              if (equ_type == EquationType::evaluate)
                {
                  eq_node = getBlockEquationExpr(block, i);
                  lhs = eq_node->arg1;
                  rhs = eq_node->arg2;
                  rhs->compile(code_file, instruction_number, false, temporary_terms_union, blocks_temporary_terms_idxs, false, false);
                  lhs->compile(code_file, instruction_number, true, temporary_terms_union, blocks_temporary_terms_idxs, false, false);
                }
              else if (equ_type == EquationType::evaluateRenormalized)
                {
                  eq_node = getBlockEquationRenormalizedExpr(block, i);
                  lhs = eq_node->arg1;
                  rhs = eq_node->arg2;
                  rhs->compile(code_file, instruction_number, false, temporary_terms_union, blocks_temporary_terms_idxs, false, false);
                  lhs->compile(code_file, instruction_number, true, temporary_terms_union, blocks_temporary_terms_idxs, false, false);
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
              FNUMEXPR_ fnumexpr(ExpressionType::ModelEquation, getBlockEquationID(block, i));
              fnumexpr.write(code_file, instruction_number);
              eq_node = getBlockEquationExpr(block, i);
              lhs = eq_node->arg1;
              rhs = eq_node->arg2;
              lhs->compile(code_file, instruction_number, false, temporary_terms_union, blocks_temporary_terms_idxs, false, false);
              rhs->compile(code_file, instruction_number, false, temporary_terms_union, blocks_temporary_terms_idxs, false, false);

              FBINARY_ fbinary{static_cast<int>(BinaryOpcode::minus)};
              fbinary.write(code_file, instruction_number);

              FSTPR_ fstpr(i - block_recursive);
              fstpr.write(code_file, instruction_number);
            }
        }
      FENDEQU_ fendequ_l;
      fendequ_l.write(code_file, instruction_number);

      // The Jacobian if we have to solve the block determinsitic bloc

      // Write temporary terms for derivatives
      write_eq_tt(blocks[block].size);

      switch (simulation_type)
        {
        case BlockSimulationType::solveBackwardSimple:
        case BlockSimulationType::solveForwardSimple:
          {
            FNUMEXPR_ fnumexpr(ExpressionType::FirstEndoDerivative, 0, 0);
            fnumexpr.write(code_file, instruction_number);
          }
          compileDerivative(code_file, instruction_number, getBlockEquationID(block, 0), getBlockVariableID(block, 0), temporary_terms_union, blocks_temporary_terms_idxs);
          {
            FSTPG2_ fstpg2(0, 0);
            fstpg2.write(code_file, instruction_number);
          }
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
              FNUMEXPR_ fnumexpr(ExpressionType::FirstEndoDerivative, eqr, varr, 0);
              fnumexpr.write(code_file, instruction_number);

              compileChainRuleDerivative(code_file, instruction_number, block, eq, var, 0, temporary_terms_union, blocks_temporary_terms_idxs);

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
StaticModel::computingPass(int derivsOrder, int paramsDerivsOrder, const eval_context_t &eval_context, bool no_tmp_terms, bool block, bool bytecode)
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

  // Compute derivatives w.r. to all endogenous
  set<int> vars;
  for (int i = 0; i < symbol_table.endo_nbr(); i++)
    {
      int id = symbol_table.getID(SymbolType::endogenous, i);
      //      if (!symbol_table.isAuxiliaryVariableButNotMultiplier(id))
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
  vector<ostringstream> d_output(derivatives.size()); // Derivatives output (at all orders, including 0=residual)
  vector<ostringstream> tt_output(derivatives.size()); // Temp terms output (at all orders)

  ExprNodeOutputType output_type = (use_dll ? ExprNodeOutputType::CStaticModel :
                                    julia ? ExprNodeOutputType::juliaStaticModel : ExprNodeOutputType::matlabStaticModel);

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
  int JacobianColsNbr = symbol_table.endo_nbr();
  int hessianColsNbr = JacobianColsNbr*JacobianColsNbr;

  auto getJacobCol = [this](int var) { return symbol_table.getTypeSpecificID(getSymbIDByDerivID(var)); };

  // Write Jacobian w.r. to endogenous only
  if (!derivatives[1].empty())
    {
      writeTemporaryTerms(temporary_terms_derivatives[1],
                          temp_term_union,
                          temporary_terms_idxs,
                          tt_output[1], output_type, tef_terms);

      for (const auto & [indices, d1] : derivatives[1])
        {
          auto [eq, var] = vectorToTuple<2>(indices);

          jacobianHelper(d_output[1], eq, getJacobCol(var), output_type);
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
                col_idx *= JacobianColsNbr;
                col_idx += getJacobCol(vidx[j]);
              }

            if (output_type == ExprNodeOutputType::juliaStaticModel)
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
                int col_idx_sym = getJacobCol(vidx[2]) * JacobianColsNbr + getJacobCol(vidx[1]);

                if (output_type == ExprNodeOutputType::juliaStaticModel)
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
        if (output_type != ExprNodeOutputType::juliaStaticModel)
          d_output[i] << i_output.str() << j_output.str() << v_output.str();
      }

  if (output_type == ExprNodeOutputType::matlabStaticModel)
    {
      // Check that we don't have more than 32 nested parenthesis because Matlab does not suppor this. See Issue #1201
      map<string, string> tmp_paren_vars;
      bool message_printed = false;
      for (auto &it : tt_output)
        fixNestedParenthesis(it, tmp_paren_vars, message_printed);
      for (auto &it : d_output)
        fixNestedParenthesis(it, tmp_paren_vars, message_printed);

      ostringstream init_output, end_output;
      init_output << "residual = zeros(" << equations.size() << ", 1);";
      end_output << "if ~isreal(residual)" << endl
                 << "  residual = real(residual)+imag(residual).^2;" << endl
                 << "end";
      writeStaticModelHelper(basename, "static_resid", "residual", "static_resid_tt",
                             temporary_terms_mlv.size() + temporary_terms_derivatives[0].size(),
                             "", init_output, end_output,
                             d_output[0], tt_output[0]);

      init_output.str("");
      end_output.str("");
      init_output << "g1 = zeros(" << equations.size() << ", " << symbol_table.endo_nbr() << ");";
      end_output << "if ~isreal(g1)" << endl
                 << "    g1 = real(g1)+2*imag(g1);" << endl
                 << "end";
      writeStaticModelHelper(basename, "static_g1", "g1", "static_g1_tt",
                             temporary_terms_mlv.size() + temporary_terms_derivatives[0].size() + temporary_terms_derivatives[1].size(),
                             "static_resid_tt",
                             init_output, end_output,
                             d_output[1], tt_output[1]);
      writeWrapperFunctions(basename, "g1");

      // For order ≥ 2
      int ncols = JacobianColsNbr;
      int ntt = temporary_terms_mlv.size() + temporary_terms_derivatives[0].size() + temporary_terms_derivatives[1].size();
      for (size_t i = 2; i < derivatives.size(); i++)
        {
          ncols *= JacobianColsNbr;
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
          writeStaticModelHelper(basename, "static_" + gname, gname,
                                 "static_" + gname + "_tt",
                                 ntt,
                                 "static_" + gprevname + "_tt",
                                 init_output, end_output,
                                 d_output[i], tt_output[i]);
          if (i <= 3)
            writeWrapperFunctions(basename, gname);
        }

      writeStaticMatlabCompatLayer(basename);
    }
  else if (output_type == ExprNodeOutputType::CStaticModel)
    {
      for (size_t i = 0; i < d_output.size(); i++)
        {
          string funcname = i == 0 ? "resid" : "g" + to_string(i);
          StaticOutput << "void static_" << funcname << "_tt(const double *restrict y, const double *restrict x, const double *restrict params, double *restrict T)" << endl
                       << "{" << endl
                       << tt_output[i].str()
                       << "}" << endl
                       << endl
                       << "void static_" << funcname << "(const double *restrict y, const double *restrict x, const double *restrict params, const double *restrict T, ";
          if (i == 0)
            StaticOutput << "double *restrict residual";
          else if (i == 1)
            StaticOutput << "double *restrict g1";
          else
            StaticOutput << "double *restrict " << funcname << "_i, double *restrict " << funcname << "_j, double *restrict " << funcname << "_v";
          StaticOutput << ")" << endl
                       << "{" << endl;
          if (i == 0)
            StaticOutput << "  double lhs, rhs;" << endl;
          StaticOutput << d_output[i].str()
                       << "}" << endl
                       << endl;
        }
    }
  else
    {
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
             << "tmp_nbr[1] = " << temporary_terms_mlv.size() + temporary_terms_derivatives[0].size() << "# Number of temporary terms for the residuals" << endl
             << "tmp_nbr[2] = " << temporary_terms_derivatives[1].size() << "# Number of temporary terms for g1 (jacobian)" << endl
             << "tmp_nbr[3] = " << temporary_terms_derivatives[2].size() << "# Number of temporary terms for g2 (hessian)" << endl
             << "tmp_nbr[4] = " << temporary_terms_derivatives[3].size() << "# Number of temporary terms for g3 (third order derivates)" << endl << endl;

      // staticResidTT!
      output << "function staticResidTT!(T::Vector{Float64}," << endl
             << "                        y::Vector{Float64}, x::Vector{Float64}, params::Vector{Float64})" << endl
             << "    @assert length(T) >= " << temporary_terms_mlv.size() + temporary_terms_derivatives[0].size()  << endl
             << tt_output[0].str()
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
             << d_output[0].str()
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
             << tt_output[1].str()
             << "    return nothing" << endl
             << "end" << endl << endl;

      // staticG1!
      output << "function staticG1!(T::Vector{Float64}, g1::Matrix{Float64}," << endl
             << "                   y::Vector{Float64}, x::Vector{Float64}, params::Vector{Float64}, T1_flag::Bool, T0_flag::Bool)" << endl
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
             << d_output[1].str()
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
             << tt_output[2].str()
             << "    return nothing" << endl
             << "end" << endl << endl;

      // staticG2!
      output << "function staticG2!(T::Vector{Float64}, g2::Matrix{Float64}," << endl
             << "                   y::Vector{Float64}, x::Vector{Float64}, params::Vector{Float64}, T2_flag::Bool, T1_flag::Bool, T0_flag::Bool)" << endl
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
             << d_output[2].str()
             << "    return nothing" << endl
             << "end" << endl << endl;

      // staticG3TT!
      output << "function staticG3TT!(T::Vector{Float64}," << endl
             << "                     y::Vector{Float64}, x::Vector{Float64}, params::Vector{Float64}, T2_flag::Bool, T1_flag::Bool, T0_flag::Bool)" << endl
             << "    if T2_flag" << endl
             << "        staticG2TT!(T, y, x, params, T1_flag, T0_flag)" << endl
             << "    end" << endl
             << tt_output[3].str()
             << "    return nothing" << endl
             << "end" << endl << endl;

      // staticG3!
      int ncols = hessianColsNbr * JacobianColsNbr;
      output << "function staticG3!(T::Vector{Float64}, g3::Matrix{Float64}," << endl
             << "                   y::Vector{Float64}, x::Vector{Float64}, params::Vector{Float64}, T3_flag::Bool, T2_flag::Bool, T1_flag::Bool, T0_flag::Bool)" << endl
             << "    @assert length(T) >= "
             << temporary_terms_mlv.size() + temporary_terms_derivatives[0].size() + temporary_terms_derivatives[1].size() + temporary_terms_derivatives[2].size() + temporary_terms_derivatives[3].size() << endl
             << "    @assert size(g3) == (" << nrows << ", " << ncols << ")" << endl
             << "    @assert length(y) == " << symbol_table.endo_nbr() << endl
             << "    @assert length(x) == " << symbol_table.exo_nbr() << endl
             << "    @assert length(params) == " << symbol_table.param_nbr() << endl
             << "    if T3_flag" << endl
             << "        staticG3TT!(T, y, x, params, T2_flag, T1_flag, T0_flag)" << endl
             << "    end" << endl
             << "    fill!(g3, 0.0)" << endl
             << d_output[3].str()
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
             << "    staticG1!(T, g1, y, x, params, true, true)" << endl
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
}

void
StaticModel::writeStaticCFile(const string &basename) const
{
  // Writing comments and function definition command
  string filename = basename + "/model/src/static.c";

  int ntt = temporary_terms_mlv.size() + temporary_terms_derivatives[0].size() + temporary_terms_derivatives[1].size() + temporary_terms_derivatives[2].size() + temporary_terms_derivatives[3].size();

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
         << endl
         << "#include <math.h>" << endl
         << "#include <stdlib.h>" << endl
         << R"(#include "mex.h")" << endl
         << endl;

  // Write function definition if BinaryOpcode::powerDeriv is used
  writePowerDeriv(output);

  output << endl;

  writeStaticModel(output, true, false);

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
  writeStaticModel(basename, false, true);
}

void
StaticModel::writeStaticFile(const string &basename, bool block, bool bytecode, bool use_dll, const string &mexext, const filesystem::path &matlabroot, const filesystem::path &dynareroot, bool julia) const
{
  filesystem::path model_dir{basename};
  model_dir /= "model";
  if (use_dll)
    filesystem::create_directories(model_dir / "src");
  if (bytecode)
    filesystem::create_directories(model_dir / "bytecode");

  if (block)
    {
      if (bytecode)
        writeStaticBlockBytecode(basename);
      else if (use_dll)
        {
          writeStaticPerBlockCFiles(basename);
          writeStaticBlockCFile(basename);
          vector<filesystem::path> src_files{blocks.size() + 1};
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
      if (bytecode)
        writeStaticBytecode(basename);
      else if (use_dll)
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
    {
      set<int> result;
      equation->collectVariables(SymbolType::exogenous, result);
      equation->collectVariables(SymbolType::exogenousDet, result);
      if (!result.empty())
        return true;
    }
  return false;
}

void
StaticModel::writeStaticBlockMFile(const string &basename) const
{
  string filename = packageDir(basename) + "/static.m";

  ofstream output;
  output.open(filename, ios::out | ios::binary);
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

  ofstream output;
  output.open(filename, ios::out | ios::binary);
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
        int var = symbol_table.getTypeSpecificID(getSymbIDByDerivID(deriv_id));
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

  string func_name = julia ? basename + "_set_auxiliary_variables!" : "set_auxiliary_variables";
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
         << comment << "           from model file (.mod)" << endl << endl
         << output_func_body.str()
         << "end" << endl;
  if (julia)
    output << "end" << endl;

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
        for (auto it = efout.begin(); it != efout.end(); ++it)
          {
            if (it != efout.begin())
              output << ", ";
            output << *it;
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
StaticModel::writeParamsDerivativesFile(const string &basename, bool julia) const
{
  if (!params_derivatives.size())
    return;

  ExprNodeOutputType output_type = (julia ? ExprNodeOutputType::juliaStaticModel : ExprNodeOutputType::matlabStaticModel);

  ostringstream tt_output; // Used for storing temporary terms
  ostringstream jacobian_output; // Used for storing jacobian equations
  ostringstream hessian_output; // Used for storing Hessian equations
  ostringstream hessian1_output; // Used for storing Hessian equations
  ostringstream third_derivs_output; // Used for storing third order derivatives equations
  ostringstream third_derivs1_output; // Used for storing third order derivatives equations

  temporary_terms_t temp_term_union;
  deriv_node_temp_terms_t tef_terms;

  writeModelLocalVariableTemporaryTerms(temp_term_union, params_derivs_temporary_terms_idxs, tt_output, output_type, tef_terms);
  for (const auto &it : params_derivs_temporary_terms)
    writeTemporaryTerms(it.second, temp_term_union, params_derivs_temporary_terms_idxs, tt_output, output_type, tef_terms);

  for (const auto & [indices, d1] : params_derivatives.find({ 0, 1 })->second)
    {
      auto [eq, param] = vectorToTuple<2>(indices);

      int param_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param)) + 1;

      jacobian_output << "rp" << LEFT_ARRAY_SUBSCRIPT(output_type)
                      <<  eq+1 << ", " << param_col
                      << RIGHT_ARRAY_SUBSCRIPT(output_type) << " = ";
      d1->writeOutput(jacobian_output, output_type, temp_term_union, params_derivs_temporary_terms_idxs, tef_terms);
      jacobian_output << ";" << endl;
    }

  for (const auto & [indices, d2] : params_derivatives.find({ 1, 1 })->second)
    {
      auto [eq, var, param] = vectorToTuple<3>(indices);

      int var_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(var)) + 1;
      int param_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param)) + 1;

      hessian_output << "gp" << LEFT_ARRAY_SUBSCRIPT(output_type)
                     << eq+1 << ", " << var_col << ", " << param_col
                     << RIGHT_ARRAY_SUBSCRIPT(output_type) << " = ";
      d2->writeOutput(hessian_output, output_type, temp_term_union, params_derivs_temporary_terms_idxs, tef_terms);
      hessian_output << ";" << endl;
    }

  int i = 1;
  for (const auto &[indices, d2] : params_derivatives.find({ 0, 2 })->second)
    {
      auto [eq, param1, param2] = vectorToTuple<3>(indices);

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
  for (const auto &[indices, d2] : params_derivatives.find({ 1, 2 })->second)
    {
      auto [eq, var, param1, param2] = vectorToTuple<4>(indices);

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
  for (const auto &[indices, d2] : params_derivatives.find({ 2, 1 })->second)
    {
      auto [eq, var1, var2, param] = vectorToTuple<4>(indices);

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
      fixNestedParenthesis(tt_output, tmp_paren_vars, message_printed);
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
                       << tt_output.str()
                       << "rp = zeros(" << equations.size() << ", "
                       << symbol_table.param_nbr() << ");" << endl
                       << jacobian_output.str()
                       << "gp = zeros(" << equations.size() << ", " << symbol_table.endo_nbr() << ", "
                       << symbol_table.param_nbr() << ");" << endl
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
    paramsDerivsFile << "module " << basename << "StaticParamsDerivs" << endl
                     << "#" << endl
                     << "# NB: this file was automatically generated by Dynare" << endl
                     << "#     from " << basename << ".mod" << endl
                     << "#" << endl
                     << "export params_derivs" << endl << endl
                     << "function params_derivs(y, x, params)" << endl
                     << tt_output.str()
                     << "rp = zeros(" << equations.size() << ", "
                     << symbol_table.param_nbr() << ");" << endl
                     << jacobian_output.str()
                     << "gp = zeros(" << equations.size() << ", " << symbol_table.endo_nbr() << ", "
                     << symbol_table.param_nbr() << ");" << endl
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

  auto getJacobCol = [this](int var) { return symbol_table.getTypeSpecificID(getSymbIDByDerivID(var)); };
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
              col_idx *= symbol_table.endo_nbr();
              col_idx += getJacobCol(vidx[j]);
            }

          if (writeDetails)
            d_output[i] << R"({"eq": )" << eq + 1;
          else
            d_output[i] << R"({"row": )" << eq + 1;

          d_output[i] << R"(, "col": )" << (i > 1 ? "[" : "") << col_idx + 1;

          if (i == 2 && vidx[1] != vidx[2]) // Symmetric elements in hessian
            {
              int col_idx_sym = getJacobCol(vidx[2]) * symbol_table.endo_nbr() + getJacobCol(vidx[1]);
              d_output[i] << ", " << col_idx_sym + 1;
            }
          if (i > 1)
            d_output[i] << "]";

          if (writeDetails)
            for (size_t j = 1; j < vidx.size(); j++)
              d_output[i] << R"(, "var)" << (i > 1 ? to_string(j) : "") << R"(": ")" << symbol_table.getName(getSymbIDByDerivID(vidx[j])) << R"(")";

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
  auto &rp = params_derivatives.find({ 0, 1 })->second;
  for (auto it = rp.begin(); it != rp.end(); ++it)
    {
      if (it != rp.begin())
        jacobian_output << ", ";

      auto [eq, param] = vectorToTuple<2>(it->first);
      expr_t d1 = it->second;

      int param_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param)) + 1;

      if (writeDetails)
        jacobian_output << R"({"eq": )" << eq + 1;
      else
        jacobian_output << R"({"row": )" << eq + 1;

      if (writeDetails)
        jacobian_output << R"(, "param_col": )" << param_col;

      jacobian_output << R"(, "param": ")" << symbol_table.getName(getSymbIDByDerivID(param)) << R"(")";

      jacobian_output << R"(, "val": ")";
      d1->writeJsonOutput(jacobian_output, temp_term_union, tef_terms);
      jacobian_output << R"("})" << endl;
    }
  jacobian_output << "]}";

  hessian_output << R"("deriv_jacobian_wrt_params": {)"
                 << R"(  "neqs": )" << equations.size()
                 << R"(, "nvarcols": )" << symbol_table.endo_nbr()
                 << R"(, "nparamcols": )" << symbol_table.param_nbr()
                 << R"(, "entries": [)";
  auto &gp = params_derivatives.find({ 1, 1 })->second;
  for (auto it = gp.begin(); it != gp.end(); ++it)
    {
      if (it != gp.begin())
        hessian_output << ", ";

      auto [eq, var, param] = vectorToTuple<3>(it->first);
      expr_t d2 = it->second;

      int var_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(var)) + 1;
      int param_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param)) + 1;

      if (writeDetails)
        hessian_output << R"({"eq": )" << eq + 1;
      else
        hessian_output << R"({"row": )" << eq + 1;

      if (writeDetails)
        hessian_output << R"(, "var": ")" << symbol_table.getName(getSymbIDByDerivID(var)) << R"(")"
                       << R"(, "param": ")" << symbol_table.getName(getSymbIDByDerivID(param)) << R"(")";

      hessian_output << R"(, "var_col": )" << var_col
                     << R"(, "param_col": )" << param_col
                     << R"(, "val": ")";
      d2->writeJsonOutput(hessian_output, temp_term_union, tef_terms);
      hessian_output << R"("})" << endl;
    }
  hessian_output << "]}";

  hessian1_output << R"("second_deriv_residuals_wrt_params": {)"
                  << R"(  "nrows": )" << equations.size()
                  << R"(, "nparam1cols": )" << symbol_table.param_nbr()
                  << R"(, "nparam2cols": )" << symbol_table.param_nbr()
                  << R"(, "entries": [)";
  auto &rpp = params_derivatives.find({ 0, 2 })->second;
  for (auto it = rpp.begin(); it != rpp.end(); ++it)
    {
      if (it != rpp.begin())
        hessian1_output << ", ";

      auto [eq, param1, param2] = vectorToTuple<3>(it->first);
      expr_t d2 = it->second;

      int param1_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param1)) + 1;
      int param2_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param2)) + 1;

      if (writeDetails)
        hessian1_output << R"({"eq": )" << eq + 1;
      else
        hessian1_output << R"({"row": )" << eq + 1;

      hessian1_output << R"(, "param1_col": )" << param1_col
                      << R"(, "param2_col": )" << param2_col;

      if (writeDetails)
        hessian1_output << R"(, "param1": ")" << symbol_table.getName(getSymbIDByDerivID(param1)) << R"(")"
                        << R"(, "param2": ")" << symbol_table.getName(getSymbIDByDerivID(param2)) << R"(")";

      hessian1_output << R"(, "val": ")";
      d2->writeJsonOutput(hessian1_output, temp_term_union, tef_terms);
      hessian1_output << R"("})" << endl;
    }
  hessian1_output << "]}";

  third_derivs_output << R"("second_deriv_jacobian_wrt_params": {)"
                      << R"(  "neqs": )" << equations.size()
                      << R"(, "nvarcols": )" << symbol_table.endo_nbr()
                      << R"(, "nparam1cols": )" << symbol_table.param_nbr()
                      << R"(, "nparam2cols": )" << symbol_table.param_nbr()
                      << R"(, "entries": [)";
  auto &gpp = params_derivatives.find({ 1, 2 })->second;
  for (auto it = gpp.begin(); it != gpp.end(); ++it)
    {
      if (it != gpp.begin())
        third_derivs_output << ", ";

      auto [eq, var, param1, param2] = vectorToTuple<4>(it->first);
      expr_t d2 = it->second;

      int var_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(var)) + 1;
      int param1_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param1)) + 1;
      int param2_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param2)) + 1;

      if (writeDetails)
        third_derivs_output << R"({"eq": )" << eq + 1;
      else
        third_derivs_output << R"({"row": )" << eq + 1;
      third_derivs_output << R"(, "var_col": )" << var_col
                          << R"(, "param1_col": )" << param1_col
                          << R"(, "param2_col": )" << param2_col;

      if (writeDetails)
        third_derivs_output << R"(, "var": ")" << symbol_table.getName(var) << R"(")"
                            << R"(, "param1": ")" << symbol_table.getName(getSymbIDByDerivID(param1)) << R"(")"
                            << R"(, "param2": ")" << symbol_table.getName(getSymbIDByDerivID(param2)) << R"(")";

      third_derivs_output << R"(, "val": ")";
      d2->writeJsonOutput(third_derivs_output, temp_term_union, tef_terms);
      third_derivs_output << R"("})" << endl;
    }
  third_derivs_output << "]}" << endl;

  third_derivs1_output << R"("derivative_hessian_wrt_params": {)"
                       << R"(  "neqs": )" << equations.size()
                       << R"(, "nvar1cols": )" << symbol_table.endo_nbr()
                       << R"(, "nvar2cols": )" << symbol_table.endo_nbr()
                       << R"(, "nparamcols": )" << symbol_table.param_nbr()
                       << R"(, "entries": [)";
  auto &hp = params_derivatives.find({ 2, 1 })->second;
  for (auto it = hp.begin(); it != hp.end(); ++it)
    {
      if (it != hp.begin())
        third_derivs1_output << ", ";

      auto [eq, var1, var2, param] = vectorToTuple<4>(it->first);
      expr_t d2 = it->second;

      int var1_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(var1)) + 1;
      int var2_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(var2)) + 1;
      int param_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param)) + 1;

      if (writeDetails)
        third_derivs1_output << R"({"eq": )" << eq + 1;
      else
        third_derivs1_output << R"({"row": )" << eq + 1;

      third_derivs1_output << R"(, "var1_col": )" << var1_col
                           << R"(, "var2_col": )" << var2_col
                           << R"(, "param_col": )" << param_col;

      if (writeDetails)
        third_derivs1_output << R"(, "var1": ")" << symbol_table.getName(getSymbIDByDerivID(var1)) << R"(")"
                             << R"(, "var2": ")" << symbol_table.getName(getSymbIDByDerivID(var2)) << R"(")"
                             << R"(, "param1": ")" << symbol_table.getName(getSymbIDByDerivID(param)) << R"(")";

      third_derivs1_output << R"(, "val": ")";
      d2->writeJsonOutput(third_derivs1_output, temp_term_union, tef_terms);
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
