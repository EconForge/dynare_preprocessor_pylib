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
#include <numeric>

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
    AddLocalVariable(it, m.local_variables_table.at(it)->toStatic(*this));

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
StaticModel::writeStaticPerBlockMFiles(const string &basename) const
{
  temporary_terms_t temporary_terms; // Temp terms written so far

  for (int blk = 0; blk < static_cast<int>(blocks.size()); blk++)
    {
      BlockSimulationType simulation_type = blocks[blk].simulation_type;

      filesystem::path filename {packageDir(basename) / "+block" / ("static_" + to_string(blk+1) + ".m")};
      ofstream output{filename, ios::out | ios::binary};
      if (!output.is_open())
        {
          cerr << "ERROR: Can't open file " << filename.string() << " for writing" << endl;
          exit(EXIT_FAILURE);
        }
      output << "%" << endl
             << "% " << filename.string() << " : Computes static version of one block" << endl
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
             << "  % //" << "                     Block "sv.substr(static_cast<int>(log10(blk + 1))) << blk+1
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

vector<filesystem::path>
StaticModel::writeStaticPerBlockCFiles(const string &basename, const string &mexext,
                                       const filesystem::path &matlabroot,
                                       const filesystem::path &dynareroot) const
{
  temporary_terms_t temporary_terms; // Temp terms written so far
  const filesystem::path model_src_dir { filesystem::path{basename} / "model" / "src" };
  vector<filesystem::path> compiled_object_files;

  for (int blk = 0; blk < static_cast<int>(blocks.size()); blk++)
    {
      BlockSimulationType simulation_type = blocks[blk].simulation_type;

      filesystem::path filename { model_src_dir / ("static_" + to_string(blk+1) + ".c") };
      ofstream output{filename, ios::out | ios::binary};
      if (!output.is_open())
        {
          cerr << "ERROR: Can't open file " << filename.string() << " for writing" << endl;
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

      // Compile intermediary object under <MODFILE>/model/src/
      compiled_object_files.emplace_back(compileMEX(model_src_dir, "static_" + to_string(blk+1),
                                                    mexext, { filename }, matlabroot, dynareroot,
                                                    false));

      filename = model_src_dir / ("static_" + to_string(blk+1) + ".h");
      ofstream header_output{filename, ios::out | ios::binary};
      if (!header_output.is_open())
        {
          cerr << "ERROR: Can't open file " << filename << " for writing" << endl;
          exit(EXIT_FAILURE);
        }
      header_output << header.str() << ';' << endl;
      header_output.close();
    }
  return compiled_object_files;
}

void
StaticModel::writeStaticBytecode(const string &basename) const
{
  // First write the .bin file
  int u_count_int { writeBytecodeBinFile(basename + "/model/bytecode/static.bin", false) };

  BytecodeWriter code_file {basename + "/model/bytecode/static.cod"};
  vector<int> eq_idx(equations.size());
  iota(eq_idx.begin(), eq_idx.end(), 0);
  vector<int> endo_idx(symbol_table.endo_nbr());
  iota(endo_idx.begin(), endo_idx.end(), 0);

  // Declare temporary terms and the (single) block
  code_file << FDIMST_{static_cast<int>(temporary_terms_derivatives[0].size()
                                        + temporary_terms_derivatives[1].size())}
            << FBEGINBLOCK_{symbol_table.endo_nbr(),
                            BlockSimulationType::solveForwardComplete,
                            0,
                            symbol_table.endo_nbr(),
                            endo_idx,
                            eq_idx,
                            false,
                            symbol_table.endo_nbr(),
                            0,
                            0,
                            u_count_int,
                            symbol_table.endo_nbr()};

  writeBytecodeHelper<false>(code_file);
}

void
StaticModel::writeStaticBlockBytecode(const string &basename) const
{
  BytecodeWriter code_file {basename + "/model/bytecode/static.cod"};

  const string bin_filename {basename + "/model/bytecode/static.bin"};
  ofstream bin_file {bin_filename, ios::out | ios::binary};
  if (!bin_file.is_open())
    {
      cerr << R"(Error : Can't open file ")" << bin_filename << R"(" for writing)" << endl;
      exit(EXIT_FAILURE);
    }

  // Temporary variables declaration
  code_file << FDIMST_{static_cast<int>(blocks_temporary_terms_idxs.size())};

  for (int block {0}; block < static_cast<int>(blocks.size()); block++)
    {
      const BlockSimulationType simulation_type {blocks[block].simulation_type};
      const int block_size {blocks[block].size};

      const int u_count {simulation_type == BlockSimulationType::solveBackwardComplete
                         || simulation_type == BlockSimulationType::solveForwardComplete
                         ? writeBlockBytecodeBinFile(bin_file, block)
                         : 0};

      code_file << FBEGINBLOCK_{blocks[block].mfs_size,
                                simulation_type,
                                blocks[block].first_equation,
                                block_size,
                                endo_idx_block2orig,
                                eq_idx_block2orig,
                                blocks[block].linear,
                                symbol_table.endo_nbr(),
                                0,
                                0,
                                u_count,
                                block_size};

      writeBlockBytecodeHelper<false>(code_file, block);
    }
  code_file << FEND_{};
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
      cerr << "ERROR: The derivatives matrix of the " << modelClassName() << " is too large. Please decrease the approximation order." << endl;
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
  cout << "Computing " << modelClassName() << " derivatives (order " << derivsOrder << ")." << endl;

  computeDerivatives(derivsOrder, vars);

  if (paramsDerivsOrder > 0)
    {
      cout << "Computing " << modelClassName() << " derivatives w.r.t. parameters (order " << paramsDerivsOrder << ")." << endl;
      computeParamsDerivatives(paramsDerivsOrder);
    }

  computeTemporaryTerms(true, no_tmp_terms);

  if (paramsDerivsOrder > 0 && !no_tmp_terms)
    computeParamsDerivativesTemporaryTerms();

  computingPassBlock(eval_context, no_tmp_terms);
  if (!block_decomposed && block)
    {
      cerr << "ERROR: Block decomposition requested but failed. If your model does not have a steady state, you may want to try the 'no_static' option of the 'model' block." << endl;
      exit(EXIT_FAILURE);
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
                         temporary_terms_derivatives[0].size(),
                         "", init_output, end_output,
                         d_output[0], tt_output[0]);

  init_output.str("");
  end_output.str("");
  init_output << "g1 = zeros(" << equations.size() << ", " << symbol_table.endo_nbr() << ");";
  end_output << "if ~isreal(g1)" << endl
             << "    g1 = real(g1)+2*imag(g1);" << endl
             << "end";
  writeStaticMFileHelper(basename, "static_g1", "g1", "static_g1_tt",
                         temporary_terms_derivatives[0].size() + temporary_terms_derivatives[1].size(),
                         "static_resid_tt",
                         init_output, end_output,
                         d_output[1], tt_output[1]);
  writeStaticMWrapperFunction(basename, "g1");

  // For order ≥ 2
  int ncols{symbol_table.endo_nbr()};
  int ntt { static_cast<int>(temporary_terms_derivatives[0].size() + temporary_terms_derivatives[1].size()) };
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

  filesystem::path filename {packageDir(basename) / (name + ".m")};
  ofstream output{filename, ios::out | ios::binary};
  if (!output.is_open())
    {
      cerr << "ERROR: Can't open file " << filename.string() << " for writing" << endl;
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
  filesystem::path filename {packageDir(basename) / (name_tt + ".m")};
  ofstream output{filename, ios::out | ios::binary};
  if (!output.is_open())
    {
      cerr << "ERROR: Can't open file " << filename.string() << " for writing" << endl;
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

  filename = packageDir(basename) / (name + ".m");
  output.open(filename, ios::out | ios::binary);
  if (!output.is_open())
    {
      cerr << "ERROR: Can't open file " << filename.string() << " for writing" << endl;
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
  filesystem::path filename {packageDir(basename) / "static.m"};
  ofstream output{filename, ios::out | ios::binary};
  if (!output.is_open())
    {
      cerr << "ERROR: Can't open file " << filename.string() << " for writing" << endl;
      exit(EXIT_FAILURE);
    }
  int ntt { static_cast<int>(temporary_terms_derivatives[0].size() + temporary_terms_derivatives[1].size() + temporary_terms_derivatives[2].size() + temporary_terms_derivatives[3].size()) };

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
         << "export static!, staticResid!, staticG1!, staticG2!, staticG3!" << endl << endl
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

  // staticResidTT!
  output << "function staticResidTT!(T::Vector{<: Real}," << endl
         << "                        y::Vector{<: Real}, x::Vector{<: Real}, params::Vector{<: Real})" << endl
         << "    @assert length(T) >= " << temporary_terms_derivatives[0].size()  << endl
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
         << temporary_terms_derivatives[0].size() + temporary_terms_derivatives[1].size() << endl
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
         << temporary_terms_derivatives[0].size() + temporary_terms_derivatives[1].size() + temporary_terms_derivatives[2].size() << endl
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
         << temporary_terms_derivatives[0].size() + temporary_terms_derivatives[1].size() + temporary_terms_derivatives[2].size() + temporary_terms_derivatives[3].size() << endl
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
         << "end" << endl;

  writeToFileIfModified(output, basename + "Static.jl");
}

void
StaticModel::writeStaticFile(const string &basename, bool block, bool use_dll, const string &mexext, const filesystem::path &matlabroot, const filesystem::path &dynareroot, bool julia) const
{
  filesystem::path model_dir{basename};
  model_dir /= "model";
  if (use_dll)
    create_directories(model_dir / "src");
  if (julia)
    create_directories(model_dir / "julia");
  else
    {
      auto plusfolder {packageDir(basename)};
      /* The following is not a duplicate of the same call from
         ModFile::writeMOutput(), because of planner_objective which needs its
         +objective subdirectory */
      create_directories(plusfolder);
      if (block && !use_dll)
        create_directories(plusfolder / "+block");
    }
  create_directories(model_dir / "bytecode");

  // Legacy representation
  if (block)
    {
      writeStaticBlockBytecode(basename);

      if (use_dll)
        {
          auto per_block_object_files { writeStaticPerBlockCFiles(basename, mexext, matlabroot, dynareroot) };
          writeStaticBlockCFile(basename, move(per_block_object_files), mexext, matlabroot, dynareroot);
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
        writeModelCFile<false>(basename, mexext, matlabroot, dynareroot);
      else if (julia)
        writeStaticJuliaFile(basename);
      else // M-files
        writeStaticMFile(basename);
    }

  // Sparse representation
  if (julia)
    writeSparseModelJuliaFiles<false>(basename);

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
  filesystem::path filename {packageDir(basename) / "static.m"};

  ofstream output{filename, ios::out | ios::binary};
  if (!output.is_open())
    {
      cerr << "ERROR: Can't open file " << filename.string() << " for writing" << endl;
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
StaticModel::writeStaticBlockCFile(const string &basename, vector<filesystem::path> per_block_object_files, const string &mexext, const filesystem::path &matlabroot, const filesystem::path &dynareroot) const
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

  per_block_object_files.push_back(filename);
  compileMEX(packageDir(basename), "static", mexext, per_block_object_files, matlabroot, dynareroot);
}

void
StaticModel::writeDriverOutput(ostream &output, bool block) const
{
  output << "M_.static_tmp_nbr = [";
  for (const auto &temporary_terms_derivative : temporary_terms_derivatives)
    output << temporary_terms_derivative.size() << "; ";
  output << "];" << endl;

  if (block)
    writeBlockDriverOutput(output);

  writeDriverSparseIndicesHelper<false>(output);
}

void
StaticModel::writeBlockDriverOutput(ostream &output) const
{
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
  output << "function ";
  if (!julia)
    output << "y = ";
  output << func_name << "(y, x, params)" << endl
         << comment << endl
         << comment << " Status : Computes Auxiliary variables of the " << modelClassName() << endl
         << comment << endl
         << comment << " Warning : this file is generated automatically by Dynare" << endl
         << comment << "           from model file (.mod)" << endl << endl;
  if (julia)
    output << "@inbounds begin" << endl;
  output << output_func_body.str()
         << "end" << endl;
  if (julia)
    output << "end" << endl;

  if (julia)
    writeToFileIfModified(output, basename + "SetAuxiliaryVariables.jl");
  else
    {
      /* Calling writeToFileIfModified() is useless here since we write inside
         a subdirectory deleted at each preprocessor run. */
      filesystem::path filename {packageDir(basename) / (func_name + ".m")};
      ofstream output_file{filename, ios::out | ios::binary};
      if (!output_file.is_open())
        {
          cerr << "ERROR: Can't open file " << filename.string() << " for writing" << endl;
          exit(EXIT_FAILURE);
        }
      output_file << output.str();
      output_file.close();
    }
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
  output << R"("static_tmp_nbr": [)";
  for (bool printed_something {false};
       const auto &tts : temporary_terms_derivatives)
    {
      if (exchange(printed_something, true))
        output << ", ";
      output << tts.size();
    }
  output << "], ";
  writeJsonSparseIndicesHelper<false>(output);
}

void
StaticModel::writeJsonComputingPassOutput(ostream &output, bool writeDetails) const
{
  auto [mlv_output, d_output] { writeJsonComputingPassOutputHelper<false>(writeDetails) };

  if (writeDetails)
    output << R"("static_model": {)";
  else
    output << R"("static_model_simple": {)";
  output << mlv_output.str();
  for (const auto &it : d_output)
    output << ", " << it.str();
  output << "}";
}

void
StaticModel::writeJsonParamsDerivatives(ostream &output, bool writeDetails) const
{
  if (!params_derivatives.size())
    return;

  auto [mlv_output, tt_output, rp_output, gp_output, rpp_output, gpp_output, hp_output, g3p_output]
    { writeJsonParamsDerivativesHelper<false>(writeDetails) };
  // g3p_output is ignored

  if (writeDetails)
    output << R"("static_model_params_derivative": {)";
  else
    output << R"("static_model_params_derivatives_simple": {)";
  output << mlv_output.str()
         << ", " << tt_output.str()
         << ", " << rp_output.str()
         << ", " << gp_output.str()
         << ", " << rpp_output.str()
         << ", " << gpp_output.str()
         << ", " << hp_output.str()
         << "}";
}
