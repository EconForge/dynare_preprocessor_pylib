/*
 * Copyright © 2003-2023 Dynare Team
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
#include <numeric>
#include <regex>
#include <sstream>
#include <string_view>

#include "DynamicModel.hh"
#include "ParsingDriver.hh"

void
DynamicModel::copyHelper(const DynamicModel &m)
{
  auto f = [this](const ExprNode *e) { return e->clone(*this); };

  for (const auto &it : m.static_only_equations)
    static_only_equations.push_back(dynamic_cast<BinaryOpNode *>(f(it)));
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
  dyn_jacobian_ncols{m.dyn_jacobian_ncols},
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
  variableMapping{m.variableMapping},
  blocks_jacob_cols_endo{m.blocks_jacob_cols_endo},
  var_expectation_functions_to_write{m.var_expectation_functions_to_write}
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
  dyn_jacobian_ncols = m.dyn_jacobian_ncols;
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
  variableMapping = m.variableMapping;
  blocks_jacob_cols_endo = m.blocks_jacob_cols_endo;

  var_expectation_functions_to_write = m.var_expectation_functions_to_write;

  copyHelper(m);

  return *this;
}

void
DynamicModel::writeDynamicBytecode(const string &basename) const
{
  // Determine the type of model (used for typing the single block)
  BlockSimulationType simulation_type;
  if (max_endo_lag > 0 && max_endo_lead > 0)
    simulation_type = BlockSimulationType::solveTwoBoundariesComplete;
  else if (max_endo_lag >= 0 && max_endo_lead == 0)
    simulation_type = BlockSimulationType::solveForwardComplete;
  else
    simulation_type = BlockSimulationType::solveBackwardComplete;

  // First write the .bin file
  int u_count_int { writeBytecodeBinFile(basename + "/model/bytecode/dynamic.bin",
                                         simulation_type == BlockSimulationType::solveTwoBoundariesComplete) };

  BytecodeWriter code_file {basename + "/model/bytecode/dynamic.cod"};

  // Declare temporary terms
  code_file << FDIMT_{static_cast<int>(temporary_terms_derivatives[0].size()
                                       + temporary_terms_derivatives[1].size())};

  // Declare the (single) block
  vector<int> exo(symbol_table.exo_nbr()), exo_det(symbol_table.exo_det_nbr());
  iota(exo.begin(), exo.end(), 0);
  iota(exo_det.begin(), exo_det.end(), 0);

  int jacobian_ncols_endo
    { static_cast<int>(count_if(dyn_jacobian_cols_table.begin(), dyn_jacobian_cols_table.end(),
                                [this](const auto &v)
                                { return getTypeByDerivID(v.first) == SymbolType::endogenous; }))
    };
  vector<int> eq_idx(equations.size());
  iota(eq_idx.begin(), eq_idx.end(), 0);
  vector<int> endo_idx(symbol_table.endo_nbr());
  iota(endo_idx.begin(), endo_idx.end(), 0);

  code_file << FBEGINBLOCK_{symbol_table.endo_nbr(),
                            simulation_type,
                            0,
                            symbol_table.endo_nbr(),
                            endo_idx,
                            eq_idx,
                            false,
                            symbol_table.endo_nbr(),
                            u_count_int,
                            jacobian_ncols_endo,
                            symbol_table.exo_det_nbr(),
                            symbol_table.exo_nbr(),
                            exo_det,
                            exo};

  writeBytecodeHelper<true>(code_file);
}

void
DynamicModel::writeDynamicBlockBytecode(const string &basename) const
{
  BytecodeWriter code_file {basename + "/model/bytecode/block/dynamic.cod"};

  const filesystem::path bin_filename {basename + "/model/bytecode/block/dynamic.bin"};
  ofstream bin_file {bin_filename, ios::out | ios::binary};
  if (!bin_file.is_open())
    {
      cerr << R"(Error : Can't open file ")" << bin_filename.string() << R"(" for writing)" << endl;
      exit(EXIT_FAILURE);
    }

  // Temporary variables declaration
  code_file << FDIMT_{static_cast<int>(blocks_temporary_terms_idxs.size())};

  for (int block {0}; block < static_cast<int>(blocks.size()); block++)
    {
      const BlockSimulationType simulation_type {blocks[block].simulation_type};

      // Write section of .bin file except for evaluate blocks and solve simple blocks
      const int u_count {simulation_type == BlockSimulationType::solveTwoBoundariesSimple
                         || simulation_type == BlockSimulationType::solveTwoBoundariesComplete
                         || simulation_type == BlockSimulationType::solveBackwardComplete
                         || simulation_type == BlockSimulationType::solveForwardComplete
                         ? writeBlockBytecodeBinFile(bin_file, block)
                         : 0};

      code_file << FBEGINBLOCK_{blocks[block].mfs_size,
                                simulation_type,
                                blocks[block].first_equation,
                                blocks[block].size,
                                endo_idx_block2orig,
                                eq_idx_block2orig,
                                blocks[block].linear,
                                symbol_table.endo_nbr(),
                                u_count,
                                static_cast<int>(blocks_jacob_cols_endo[block].size())};

      writeBlockBytecodeHelper<true>(code_file, block);
    }
  code_file << FEND_{};
}

void
DynamicModel::writeDynamicMFile(const string &basename) const
{
  auto [d_output, tt_output] = writeModelFileHelper<ExprNodeOutputType::matlabDynamicModel>();

  ostringstream init_output, end_output;
  init_output << "residual = zeros(" << equations.size() << ", 1);";
  writeDynamicMFileHelper(basename, "dynamic_resid", "residual", "dynamic_resid_tt",
                          temporary_terms_derivatives[0].size(),
                          "", init_output, end_output, d_output[0], tt_output[0]);

  init_output.str("");
  init_output << "g1 = zeros(" << equations.size() << ", " << getJacobianColsNbr(false) << ");";
  writeDynamicMFileHelper(basename, "dynamic_g1", "g1", "dynamic_g1_tt",
                          temporary_terms_derivatives[0].size() + temporary_terms_derivatives[1].size(),
                          "dynamic_resid_tt", init_output, end_output, d_output[1], tt_output[1]);
  writeDynamicMWrapperFunction(basename, "g1");

  // For order ≥ 2
  int ncols{getJacobianColsNbr(false)};
  int ntt { static_cast<int>(temporary_terms_derivatives[0].size() + temporary_terms_derivatives[1].size())};
  for (size_t i{2}; i < derivatives.size(); i++)
    {
      ncols *= getJacobianColsNbr(false);
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
      writeDynamicMFileHelper(basename, "dynamic_" + gname, gname, "dynamic_" + gname + "_tt", ntt,
                              "dynamic_" + gprevname + "_tt", init_output, end_output,
                              d_output[i], tt_output[i]);
      if (i <= 3)
        writeDynamicMWrapperFunction(basename, gname);
    }

  writeDynamicMCompatFile(basename);
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
  for (bool printed_something{false};
       int it : nonzero_hessian_eqs)
    {
      if (exchange(printed_something, true))
        output << " ";
      output << it + 1;
    }
  if (nonzero_hessian_eqs.size() != 1)
    output << "]";
}

void
DynamicModel::writeDynamicMWrapperFunction(const string &basename, const string &ending) const
{
  string name;
  if (ending == "g1")
    name = "dynamic_resid_g1";
  else if (ending == "g2")
    name = "dynamic_resid_g1_g2";
  else if (ending == "g3")
    name = "dynamic_resid_g1_g2_g3";

  filesystem::path filename {packageDir(basename) / (name + ".m")};
  ofstream output{filename, ios::out | ios::binary};
  if (!output.is_open())
    {
      cerr << "ERROR: Can't open file " << filename.string() << " for writing" << endl;
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
DynamicModel::writeDynamicMFileHelper(const string &basename,
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

  filename = packageDir(basename) / (name + ".m");
  output.open(filename, ios::out | ios::binary);
  if (!output.is_open())
    {
      cerr << "ERROR: Can't open file " << filename.string() << " for writing" << endl;
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
DynamicModel::writeDynamicMCompatFile(const string &basename) const
{
  filesystem::path filename {packageDir(basename) / "dynamic.m"};
  ofstream output{filename, ios::out | ios::binary};
  if (!output.is_open())
    {
      cerr << "ERROR: Can't open file " << filename.string() << " for writing" << endl;
      exit(EXIT_FAILURE);
    }
  int ntt { static_cast<int>(temporary_terms_derivatives[0].size() + temporary_terms_derivatives[1].size() + temporary_terms_derivatives[2].size() + temporary_terms_derivatives[3].size()) };

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
DynamicModel::writeDynamicJacobianNonZeroEltsFile(const string &basename) const
{
  vector<pair<int, int>> nzij_pred, nzij_current, nzij_fwrd; // pairs (tsid, equation)
  for (const auto &[indices, d1] : derivatives[1])
    {
      if (getTypeByDerivID(indices[1]) != SymbolType::endogenous)
        continue;
      int tsid { getTypeSpecificIDByDerivID(indices[1]) };
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

  const filesystem::path filename {packageDir(basename) / "dynamic_g1_nz.m"};
  ofstream output{filename, ios::out | ios::binary};
  if (!output.is_open())
    {
      cerr << "ERROR: Can't open file " << filename.string() << " for writing" << endl;
      exit(EXIT_FAILURE);
    }
  output << "function [nzij_pred, nzij_current, nzij_fwrd] = dynamic_g1_nz()" << endl
         << "% Returns the coordinates of non-zero elements in the Jacobian, in column-major order, for each lead/lag (only for endogenous)" << endl;
  auto print_nzij = [&output](const vector<pair<int, int>> &nzij, const string &name) {
                      output << "  " << name << " = zeros(" << nzij.size() << ", 2, 'int32');" << endl;
                      for (int idx{1};
                           const auto &it : nzij)
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

vector<pair<string, string>>
DynamicModel::parseIncludeExcludeEquations(const string &inc_exc_option_value, bool exclude_eqs)
{
  auto removeLeadingTrailingWhitespace = [](string &str)
  {
    str.erase(0, str.find_first_not_of("\t\n\v\f\r "));
    str.erase(str.find_last_not_of("\t\n\v\f\r ") + 1);
  };

  string tags;
  if (filesystem::exists(inc_exc_option_value))
    {
      ifstream exclude_file;
      exclude_file.open(inc_exc_option_value, ifstream::in);
      if (!exclude_file.is_open())
        {
          cerr << "ERROR: Could not open " << inc_exc_option_value << endl;
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
    tags = inc_exc_option_value;
  removeLeadingTrailingWhitespace(tags);

  if (tags.front() == '[' && tags.back() != ']')
    {
      cerr << "ERROR: " << (exclude_eqs ? "exclude_eqs" : "include_eqs")
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
      cerr << "ERROR: " << (exclude_eqs ? "exclude_eqs" : "include_eqs")
           << ": argument is of incorrect format." << endl;
      exit(EXIT_FAILURE);
    }

  vector<pair<string, string>> eq_tag_set;
  regex s(quote_regex + "|" + non_quote_regex);
  for (auto it = sregex_iterator(tags.begin(), tags.end(), s);
       it != sregex_iterator(); ++it)
    {
      string_view str {it->str()};
      if (str.front() == '\'' && str.back() == '\'')
        {
          str.remove_prefix(1);
          str.remove_suffix(1);
        }
      eq_tag_set.emplace_back(tagname, str);
    }
  return eq_tag_set;
}

vector<int>
DynamicModel::removeEquationsHelper(set<pair<string, string>> &listed_eqs_by_tag, bool exclude_eqs,
                                    bool excluded_vars_change_type,
                                    vector<BinaryOpNode *> &all_equations,
                                    vector<optional<int>> &all_equations_lineno,
                                    EquationTags &all_equation_tags, bool static_equations) const
{
  if (all_equations.empty())
    return {};

  /* Try to convert the list of equations by tags into a list of equation
     numbers.
     The tag pairs that match an equation are removed from the list, so that
     the caller knows which tag pairs have not been handled. */
  set<int> listed_eqs_by_number;
  for (auto it = listed_eqs_by_tag.begin(); it != listed_eqs_by_tag.end();)
    if (auto tmp = all_equation_tags.getEqnsByTag(it->first, it->second);
        !tmp.empty())
      {
        listed_eqs_by_number.insert(tmp.begin(), tmp.end());
        it = listed_eqs_by_tag.erase(it);
      }
    else
      ++it;

  // Compute the indices of equations to be actually deleted
  set<int> eqs_to_delete_by_number;
  if (exclude_eqs)
    eqs_to_delete_by_number = listed_eqs_by_number;
  else
    for (size_t i = 0; i < all_equations.size(); i++)
      if (!listed_eqs_by_number.contains(i))
        eqs_to_delete_by_number.insert(i);

  // remove from equations, equations_lineno, equation_tags
  vector<BinaryOpNode *> new_equations;
  vector<optional<int>> new_equations_lineno;
  map<int, int> old_eqn_num_2_new;
  vector<int> excluded_vars;
  for (size_t i = 0; i < all_equations.size(); i++)
    if (eqs_to_delete_by_number.contains(i))
      {
        if (excluded_vars_change_type)
          if (auto tmp = all_equation_tags.getTagValueByEqnAndKey(i, "endogenous"); tmp)
            excluded_vars.push_back(symbol_table.getID(*tmp));
          else
            {
              set<int> result;
              all_equations[i]->arg1->collectVariables(SymbolType::endogenous, result);
              if (result.size() == 1)
                excluded_vars.push_back(*result.begin());
              else
                {
                  cerr << "ERROR: Equation " << i+1
                       << " has been excluded but it does not have a single variable on its left-hand side or an `endogenous` tag" << endl;
                  exit(EXIT_FAILURE);
                }
            }
      }
    else
      {
        new_equations.emplace_back(all_equations[i]);
        old_eqn_num_2_new[i] = new_equations.size() - 1;
        new_equations_lineno.emplace_back(all_equations_lineno[i]);
      }
  int n_excl = all_equations.size() - new_equations.size();

  all_equations = new_equations;
  all_equations_lineno = new_equations_lineno;

  all_equation_tags.erase(eqs_to_delete_by_number, old_eqn_num_2_new);

  if (!static_equations)
    for (size_t i = 0; i < excluded_vars.size(); i++)
      for (size_t j = i+1; j < excluded_vars.size(); j++)
        if (excluded_vars[i] == excluded_vars[j])
          {
            cerr << "ERROR: Variable " << symbol_table.getName(i) << " was excluded twice"
                 << " via a model_remove or model_replace statement, or via the include_eqs or exclude_eqs option" << endl;
            exit(EXIT_FAILURE);
          }

  cout << "Excluded " << n_excl << (static_equations ? " static " : " dynamic ")
       << "equation" << (n_excl > 1 ? "s" : "") << " via model_remove or model_replace statement, or via include_eqs or exclude_eqs option" << endl;

  return excluded_vars;
}

void
DynamicModel::removeEquations(const vector<pair<string, string>> &listed_eqs_by_tag, bool exclude_eqs,
                              bool excluded_vars_change_type)
{
  /* Convert the const vector to a (mutable) set */
  set listed_eqs_by_tag2(listed_eqs_by_tag.begin(), listed_eqs_by_tag.end());

  vector<int> excluded_vars = removeEquationsHelper(listed_eqs_by_tag2, exclude_eqs,
                                                    excluded_vars_change_type,
                                                    equations, equations_lineno,
                                                    equation_tags, false);

  // Ignore output because variables are not excluded when equations marked 'static' are excluded
  removeEquationsHelper(listed_eqs_by_tag2, exclude_eqs, excluded_vars_change_type,
                        static_only_equations, static_only_equations_lineno,
                        static_only_equations_equation_tags, true);

  if (!listed_eqs_by_tag2.empty())
    {
      cerr << "ERROR: model_remove/model_replace/exclude_eqs/include_eqs: The equations specified by" << endl;
      for (const auto &[tagname, tagvalue] : listed_eqs_by_tag)
        cerr << " " << tagname << "=" << tagvalue << endl;
      cerr << "were not found." << endl;
      exit(EXIT_FAILURE);
    }

  if (excluded_vars_change_type)
    {
      // Collect list of used variables in updated list of equations
      set<int> eqn_vars;
      for (auto eqn : equations)
        eqn->collectVariables(SymbolType::endogenous, eqn_vars);
      for (auto eqn : static_only_equations)
        eqn->collectVariables(SymbolType::endogenous, eqn_vars);

      /* Change type of endogenous variables determined by excluded equations.
         They become exogenous if they are still used somewhere, otherwise they are
         completely excluded from the model. */
      for (auto ev : excluded_vars)
        if (eqn_vars.contains(ev))
          {
            symbol_table.changeType(ev, SymbolType::exogenous);
            cerr << "Variable '" << symbol_table.getName(ev) << "' turned into an exogenous, as its defining equation has been removed (but it still appears in an equation)" << endl;
          }
        else
          {
            symbol_table.changeType(ev, SymbolType::excludedVariable);
            cerr << "Variable '" << symbol_table.getName(ev) << "' has been excluded from the model, as its defining equation has been removed and it appears nowhere else" << endl;
          }
    }
}

void
DynamicModel::includeExcludeEquations(const string &inc_exc_option_value, bool exclude_eqs)
{
  if (inc_exc_option_value.empty())
    return;

  auto listed_eqs_by_tag = parseIncludeExcludeEquations(inc_exc_option_value, exclude_eqs);

  removeEquations(listed_eqs_by_tag, exclude_eqs, true);

  /* There is already a check about #static and #dynamic in
     ModFile::checkPass(), but the present method is called from
     ModFile::transformPass(), so we must do the check again */
  if (staticOnlyEquationsNbr() != dynamicOnlyEquationsNbr())
    {
      cerr << "ERROR: exclude_eqs/include_eqs: You must remove the same number of equations marked `static` as equations marked `dynamic`." << endl;
      exit(EXIT_FAILURE);
    }
}

void
DynamicModel::writeBlockDriverOutput(ostream &output) const
{
  output << "M_.block_structure.time_recursive = " << boolalpha << time_recursive_block_decomposition << ";" << endl;

  for (int blk = 0; blk < static_cast<int>(blocks.size()); blk++)
    {
      int block_size = blocks[blk].size;
      output << "M_.block_structure.block(" << blk+1 << ").Simulation_Type = " << static_cast<int>(blocks[blk].simulation_type) << ";" << endl
             << "M_.block_structure.block(" << blk+1 << ").endo_nbr = " << block_size << ";" << endl
             << "M_.block_structure.block(" << blk+1 << ").mfs = " << blocks[blk].mfs_size << ";" << endl
             << "M_.block_structure.block(" << blk+1 << ").equation = [";
      for (int eq = 0; eq < block_size; eq++)
        output << " " << getBlockEquationID(blk, eq)+1;
      output << "];" << endl
             << "M_.block_structure.block(" << blk+1 << ").variable = [";
      for (int var = 0; var < block_size; var++)
        output << " " << getBlockVariableID(blk, var)+1;
      output << "];" << endl
             << "M_.block_structure.block(" << blk+1 << ").is_linear = " << boolalpha << blocks[blk].linear << ';' << endl
             << "M_.block_structure.block(" << blk+1 << ").NNZDerivatives = " << blocks_derivatives[blk].size() << ';' << endl;
    }

  writeBlockDriverSparseIndicesHelper<true>(output);

  output << "M_.block_structure.variable_reordered = [";
  for (int i = 0; i < symbol_table.endo_nbr(); i++)
    output << " " << endo_idx_block2orig[i]+1;
  output << "];" << endl
         << "M_.block_structure.equation_reordered = [";
  for (int i = 0; i < symbol_table.endo_nbr(); i++)
    output << " " << eq_idx_block2orig[i]+1;
  output << "];" << endl;

  map<int, set<pair<int, int>>> lag_row_incidence;
  for (const auto &[indices, d1] : derivatives[1])
    if (int deriv_id = indices[1];
        getTypeByDerivID(deriv_id) == SymbolType::endogenous)
      {
        int eq = indices[0];
        int var { getTypeSpecificIDByDerivID(deriv_id) };
        int lag = getLagByDerivID(deriv_id);
        lag_row_incidence[lag].insert({ eq, var });
      }
  for (auto [lag, eq_var_set] : lag_row_incidence)
    {
      output << "M_.block_structure.incidence(" << max_endo_lag+lag+1 << ").lead_lag = " << lag << ";" << endl
             << "M_.block_structure.incidence(" << max_endo_lag+lag+1 << ").sparse_IM = [" << endl;
      for (auto [eq, var] : eq_var_set)
        output << " " << eq+1 << " " << var+1 << ";" << endl;
      output << "];" << endl;
    }
  output << "M_.block_structure.dyn_tmp_nbr = " << blocks_temporary_terms_idxs.size() << ';' << endl;
}

void
DynamicModel::writeDriverOutput(ostream &output, bool compute_xrefs) const
{
  /* Writing initialisation for M_.lead_lag_incidence matrix
     M_.lead_lag_incidence is a matrix with as many columns as there are
     endogenous variables and as many rows as there are periods in the
     models (nbr of rows = M_.max_lag+M_.max_lead+1)

     The matrix elements are equal to zero if a variable isn't present in the
     model at a given period.
  */

  output << "M_.orig_maximum_endo_lag = " << max_endo_lag_orig << ";" << endl
         << "M_.orig_maximum_endo_lead = " << max_endo_lead_orig << ";" << endl
         << "M_.orig_maximum_exo_lag = " << max_exo_lag_orig << ";" << endl
         << "M_.orig_maximum_exo_lead = " << max_exo_lead_orig << ";" << endl
         << "M_.orig_maximum_exo_det_lag = " << max_exo_det_lag_orig << ";" << endl
         << "M_.orig_maximum_exo_det_lead = " << max_exo_det_lead_orig << ";" << endl
         << "M_.orig_maximum_lag = " << max_lag_orig << ";" << endl
         << "M_.orig_maximum_lead = " << max_lead_orig << ";" << endl
         << "M_.orig_maximum_lag_with_diffs_expanded = " << max_lag_with_diffs_expanded_orig << ";" << endl
         << "M_.lead_lag_incidence = [";
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
              output << " " << getJacobianCol(varID, false) + 1;
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
  output << "]';" << endl
         << "M_.nstatic = " << nstatic << ";" << endl
         << "M_.nfwrd   = " << nfwrd   << ";" << endl
         << "M_.npred   = " << npred   << ";" << endl
         << "M_.nboth   = " << nboth   << ";" << endl
         << "M_.nsfwrd   = " << nfwrd+nboth   << ";" << endl
         << "M_.nspred   = " << npred+nboth   << ";" << endl
         << "M_.ndynamic   = " << npred+nboth+nfwrd << ";" << endl
         << "M_.dynamic_tmp_nbr = [";
  for (const auto &tts : temporary_terms_derivatives)
    output << tts.size() << "; ";
  output << "];" << endl;

  // Write equation tags
  equation_tags.writeOutput(output);

  // Write mapping for variables and equations they are present in
  for (const auto &variable : variableMapping)
    {
      output << "M_.mapping." << symbol_table.getName(variable.first) << ".eqidx = [";
      for (auto equation : variable.second)
        output << equation + 1 << " ";
      output << "];" << endl;
    }

  /* Say if static and dynamic models differ (because of [static] and [dynamic]
     equation tags) */
  output << "M_.static_and_dynamic_models_differ = "
         << boolalpha << (static_only_equations.size() > 0)
         << ";" << endl;

  // Say if model contains an external function call
  bool has_external_function = false;
  for (auto equation : equations)
    if (equation->containsExternalFunction())
      {
        has_external_function = true;
        break;
      }
  output << "M_.has_external_function = " << boolalpha << has_external_function
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
  if (block_decomposed)
    writeBlockDriverOutput(output);

  output << "M_.state_var = [";
  for (int it : state_var)
    output << it+1 << " ";
  output << "];" << endl;

  // Writing initialization for some other variables
  output << "M_.exo_names_orig_ord = [1:" << symbol_table.exo_nbr() << "];" << endl;

  output << "M_.maximum_lag = " << max_lag << ";" << endl
         << "M_.maximum_lead = " << max_lead << ";" << endl;

  output << "M_.maximum_endo_lag = " << max_endo_lag << ";" << endl
         << "M_.maximum_endo_lead = " << max_endo_lead << ";" << endl
         << "oo_.steady_state = zeros(" << symbol_table.endo_nbr() << ", 1);" << endl;

  output << "M_.maximum_exo_lag = " << max_exo_lag << ";" << endl
         << "M_.maximum_exo_lead = " << max_exo_lead << ";" << endl
         << "oo_.exo_steady_state = zeros(" << symbol_table.exo_nbr() << ", 1);" << endl;

  if (symbol_table.exo_det_nbr())
    {
      output << "M_.maximum_exo_det_lag = " << max_exo_det_lag << ";" << endl
             << "M_.maximum_exo_det_lead = " << max_exo_det_lead << ";" << endl
             << "oo_.exo_det_steady_state = zeros(" << symbol_table.exo_det_nbr() << ", 1);" << endl;
    }

  output << "M_.params = " << "NaN(" << symbol_table.param_nbr() << ", 1);" << endl;

  string empty_cell = "cell(" + to_string(symbol_table.endo_nbr()) + ", 1)";
  output << "M_.endo_trends = struct('deflator', " << empty_cell
         << ", 'log_deflator', " << empty_cell << ", 'growth_factor', " << empty_cell
         << ", 'log_growth_factor', " << empty_cell << ");" << endl;
  for (int i = 0; i < symbol_table.endo_nbr(); i++)
    {
      int symb_id = symbol_table.getID(SymbolType::endogenous, i);
      if (auto it = nonstationary_symbols_map.find(symb_id); it != nonstationary_symbols_map.end())
        {
          auto [is_log, deflator] = it->second;
          output << "M_.endo_trends(" << i+1 << ")."
                 << (is_log ? "log_deflator" : "deflator") << " = '";
          deflator->writeJsonOutput(output, {}, {});
          output << "';" << endl;

          auto growth_factor = const_cast<DynamicModel *>(this)->AddDivide(deflator, deflator->decreaseLeadsLags(1))->removeTrendLeadLag(trend_symbols_map)->replaceTrendVar();
          output << "M_.endo_trends(" << i+1 << ")."
                 << (is_log ? "log_growth_factor" : "growth_factor") << " = '";
          growth_factor->writeJsonOutput(output, {}, {});
          output << "';" << endl;
        }
    }

  if (compute_xrefs)
    writeXrefs(output);

  // Write number of non-zero derivatives
  // Use -1 if the derivatives have not been computed
  output << "M_.NNZDerivatives = [";
  for (int i = 1; i < static_cast<int>(NNZDerivatives.size()); i++)
    output << (i > computed_derivs_order ? -1 : NNZDerivatives[i]) << "; ";
  output << "];" << endl;

  writeDriverSparseIndicesHelper<true>(output);

  // Write LHS of each equation in text form
  output << "M_.lhs = {" << endl;
  for (auto eq : equations)
    {
      output << "'";
      eq->arg1->writeJsonOutput(output, {}, {});
      output << "'; " << endl;
    }
  output << "};" << endl;
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
  for (bool var : { true, false })
    {
      map<string, vector<optional<int>>> trend_varr;
      map<string, vector<set<pair<int, int>>>> rhsr;
      for (const auto &[model_name, eqns] : (var ? var_model_table.getEqNums()
                                             : trend_component_model_table.getEqNums()))
        {
          vector<int> lhs, trend_lhs;
          vector<optional<int>> trend_var;
          vector<set<pair<int, int>>> rhs;

          if (!var)
            {
              lhs = trend_component_model_table.getLhs(model_name);
              for (auto teqn : trend_component_model_table.getTargetEqNums().at(model_name))
                {
                  int eqnidx = 0;
                  for (auto eqn : eqns)
                    {
                      if (eqn == teqn)
                        trend_lhs.push_back(lhs[eqnidx]);
                      eqnidx++;
                    }
                }
            }

          int lhs_idx = 0;
          for (auto eqn : eqns)
            {
              set<pair<int, int>> rhs_set;
              equations[eqn]->arg2->collectDynamicVariables(SymbolType::endogenous, rhs_set);
              rhs.push_back(rhs_set);

              if (!var)
                {
                  int lhs_symb_id = lhs[lhs_idx++];
                  if (symbol_table.isDiffAuxiliaryVariable(lhs_symb_id))
                    try
                      {
                        lhs_symb_id = symbol_table.getOrigSymbIdForAuxVar(lhs_symb_id);
                      }
                    catch (...)
                      {
                      }
                  optional<int> trend_var_symb_id = equations[eqn]->arg2->findTargetVariable(lhs_symb_id);
                  if (trend_var_symb_id)
                    {
                      if (symbol_table.isDiffAuxiliaryVariable(*trend_var_symb_id))
                        try
                          {
                            trend_var_symb_id = symbol_table.getOrigSymbIdForAuxVar(*trend_var_symb_id);
                          }
                        catch (...)
                          {
                          }
                      if (find(trend_lhs.begin(), trend_lhs.end(), *trend_var_symb_id) == trend_lhs.end())
                        {
                          cerr << "ERROR: trend found in trend_component equation #" << eqn << " ("
                               << symbol_table.getName(*trend_var_symb_id) << ") does not correspond to a trend equation" << endl;
                          exit(EXIT_FAILURE);
                        }
                    }
                  trend_var.push_back(move(trend_var_symb_id));
                }
            }

          rhsr[model_name] = rhs;
          if (!var)
            trend_varr[model_name] = trend_var;
        }

      if (var)
        var_model_table.setRhs(move(rhsr));
      else
        {
          trend_component_model_table.setRhs(move(rhsr));
          trend_component_model_table.setTargetVar(move(trend_varr));
        }
    }
}

void
DynamicModel::fillVarModelTable() const
{
  map<string, vector<int>> eqnums, lhsr;
  map<string, vector<expr_t>> lhs_expr_tr;
  map<string, vector<set<pair<int, int>>>> rhsr;

  for (const auto &[model_name, eqtags] : var_model_table.getEqTags())
    {
      vector<int> eqnumber, lhs;
      vector<expr_t> lhs_expr_t;
      vector<set<pair<int, int>>> rhs;

      for (const auto &eqtag : eqtags)
        {
          set<pair<int, int>> lhs_set, lhs_tmp_set, rhs_set;
          optional<int> eqn { equation_tags.getEqnByTag("name", eqtag) };
          if (!eqn)
            {
              cerr << "ERROR: no equation is named '" << eqtag << "'" << endl;
              exit(EXIT_FAILURE);
            }

          equations[*eqn]->arg1->collectDynamicVariables(SymbolType::endogenous, lhs_set);
          equations[*eqn]->arg1->collectDynamicVariables(SymbolType::exogenous, lhs_tmp_set);
          equations[*eqn]->arg1->collectDynamicVariables(SymbolType::parameter, lhs_tmp_set);

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

          eqnumber.push_back(*eqn);
          lhs.push_back(itlhs->first);
          lhs_set.clear();
          set<expr_t> lhs_expr_t_set;
          equations[*eqn]->arg1->collectVARLHSVariable(lhs_expr_t_set);
          lhs_expr_t.push_back(*(lhs_expr_t_set.begin()));

          equations[*eqn]->arg2->collectDynamicVariables(SymbolType::endogenous, rhs_set);
          rhs.push_back(rhs_set);
        }
      eqnums[model_name] = eqnumber;
      lhsr[model_name] = lhs;
      lhs_expr_tr[model_name] = lhs_expr_t;
      rhsr[model_name] = rhs;
    }
  var_model_table.setEqNums(move(eqnums));
  var_model_table.setLhs(move(lhsr));
  var_model_table.setRhs(move(rhsr));
  var_model_table.setLhsExprT(move(lhs_expr_tr));
}

void
DynamicModel::fillVarModelTableFromOrigModel() const
{
  map<string, vector<int>> lags;
  map<string, vector<optional<int>>> orig_diff_var;
  map<string, vector<bool>> diff;
  for (const auto &[model_name, eqns] : var_model_table.getEqNums())
    {
      set<expr_t> lhs;
      vector<optional<int>> orig_diff_var_vec;
      vector<bool> diff_vec;
      for (auto eqn : eqns)
        {
          // Perform some sanity checks on the RHS
          optional<string> eqtag { equation_tags.getTagValueByEqnAndKey(eqn, "name") };
          set<pair<int, int>> rhs_endo_set, rhs_exo_set;
          equations[eqn]->arg2->collectDynamicVariables(SymbolType::endogenous, rhs_endo_set);
          for (const auto &[symb_id, lag] : rhs_endo_set)
            if (lag > 0)
              {
                cerr << "ERROR: in Equation " << eqtag.value_or(to_string(eqn+1))
                     << ". A VAR model may not have leaded endogenous variables on the RHS. " << endl;
                exit(EXIT_FAILURE);
              }
            else if (!var_model_table.getStructural().at(model_name) && lag == 0)
              {
                cerr << "ERROR: in Equation " << eqtag.value_or(to_string(eqn+1))
                     << ". A non-structural VAR model may not have contemporaneous endogenous variables on the RHS. " << endl;
                exit(EXIT_FAILURE);
              }

          equations[eqn]->arg2->collectDynamicVariables(SymbolType::exogenous, rhs_exo_set);
          for (const auto &[symb_id, lag] : rhs_exo_set)
            if (lag != 0)
              {
                cerr << "ERROR: in Equation " << eqtag.value_or(to_string(eqn+1))
                     << ". A VAR model may not have lagged or leaded exogenous variables on the RHS. " << endl;
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
            orig_diff_var_vec.push_back(nullopt);

        }

      if (eqns.size() != lhs.size())
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
      for (auto eqn : eqns)
        max_lag.push_back(equations[eqn]->arg2->VarMaxLag(lhs_lag_equiv));
      lags[model_name] = max_lag;
      diff[model_name] = diff_vec;
      orig_diff_var[model_name] = orig_diff_var_vec;
    }
  var_model_table.setDiff(move(diff));
  var_model_table.setMaxLags(move(lags));
  var_model_table.setOrigDiffVar(move(orig_diff_var));
}

vector<int>
DynamicModel::getVARDerivIDs(int lhs_symb_id, int lead_lag) const
{
  vector<int> deriv_ids;

  // First directly look for the variable itself
  if (auto it = deriv_id_table.find({ lhs_symb_id, lead_lag });
      it != deriv_id_table.end())
    deriv_ids.push_back(it->second);

  // Then go through auxiliary variables
  for (auto &[key, deriv_id2] : deriv_id_table)
    {
      auto [symb_id2, lead_lag2] = key;
      const AuxVarInfo *avi;
      try
        {
          avi = &symbol_table.getAuxVarInfo(symb_id2);
        }
      catch (SymbolTable::UnknownSymbolIDException)
        {
          continue;
        }

      if (avi->type == AuxVarType::endoLag && avi->orig_symb_id.value() == lhs_symb_id
          && avi->orig_lead_lag.value() + lead_lag2 == lead_lag)
        deriv_ids.push_back(deriv_id2);

      // Handle diff lag auxvar, possibly nested several times
      int diff_lag_depth = 0;
      while (avi->type == AuxVarType::diffLag)
        {
          diff_lag_depth++;
          if (avi->orig_symb_id == lhs_symb_id && lead_lag2 - diff_lag_depth == lead_lag)
            {
              deriv_ids.push_back(deriv_id2);
              break;
            }
          try
            {
              avi = &symbol_table.getAuxVarInfo(avi->orig_symb_id.value());
            }
          catch (SymbolTable::UnknownSymbolIDException)
            {
              break;
            }
        }
    }

  return deriv_ids;
}

void
DynamicModel::fillVarModelTableMatrices()
{
  map<string, map<tuple<int, int, int>, expr_t>> AR;
  map<string, map<tuple<int, int>, expr_t>> A0;
  map<string, map<int, expr_t>> constants;
  for (const auto &[model_name, eqns] : var_model_table.getEqNums())
    {
      const vector<int> &lhs = var_model_table.getLhs(model_name);
      int max_lag = var_model_table.getMaxLag(model_name);
      for (auto lhs_symb_id : lhs)
        {
          // Fill autoregressive matrix (AR)
          for (int lag = 1; lag <= max_lag; lag++)
            {
              vector<int> deriv_ids = getVARDerivIDs(lhs_symb_id, -lag);;
              for (size_t i = 0; i < eqns.size(); i++)
                {
                  expr_t d = Zero;
                  for (int deriv_id : deriv_ids)
                    d = AddPlus(d, equations[eqns[i]]->getDerivative(deriv_id));
                  if (d != Zero)
                    {
                      if (!d->isConstant())
                        {
                          cerr << "ERROR: Equation " << equation_tags.getTagValueByEqnAndKey(eqns[i], "name").value_or(to_string(eqns[i]+1)) << " is not linear" << endl;
                          exit(EXIT_FAILURE);
                        }

                      AR[model_name][{ i, lag, lhs_symb_id }] = AddUMinus(d);
                    }
                }
            }

          // Fill A0 matrix (for contemporaneous variables)
          int lhs_deriv_id = getDerivID(lhs_symb_id, 0);
          for (size_t i = 0; i < eqns.size(); i++)
            {
              expr_t d = equations[eqns[i]]->getDerivative(lhs_deriv_id);
              if (d != Zero)
                {
                  if (!d->isConstant())
                    {
                      cerr << "ERROR: Equation " << equation_tags.getTagValueByEqnAndKey(eqns[i], "name").value_or(to_string(eqns[i]+1)) << " is not linear" << endl;
                      exit(EXIT_FAILURE);
                    }

                  A0[model_name][{ i, lhs_symb_id }] = d;
                }
            }

          // Fill constants vector
          // Constants are computed by replacing all (transformed) endos and exos by zero
          constants[model_name] = {}; // Ensure that the map exists, even if constants are all zero
          for (size_t i = 0; i < eqns.size(); i++)
            {
              auto rhs = equations[eqns[i]]->arg2;
              map<VariableNode *, NumConstNode *> subst_table;
              auto rhs_vars = var_model_table.getRhs(model_name)[i]; // All the (transformed) endogenous on RHS, as computed by updateVarAndTrendModel()
              rhs->collectDynamicVariables(SymbolType::exogenous, rhs_vars); // Add exos
              for (auto [symb_id, lag] : rhs_vars)
                subst_table[AddVariable(symb_id, lag)] = Zero;
              expr_t c = rhs->replaceVarsInEquation(subst_table);
              if (c != Zero)
                constants[model_name][i] = c;
            }
        }
    }
  var_model_table.setAR(move(AR));
  var_model_table.setA0(move(A0));
  var_model_table.setConstants(move(constants));
}

map<string, map<tuple<int, int, int>, expr_t>>
DynamicModel::computeAutoregressiveMatrices() const
{
  map<string, map<tuple<int, int, int>, expr_t>> ARr;
  for (const auto &[model_name, eqns] : trend_component_model_table.getNonTargetEqNums())
    {
      map<tuple<int, int, int>, expr_t> AR;
      const vector<int> &lhs = trend_component_model_table.getNonTargetLhs(model_name);
      for (int i{0};
           auto eqn : eqns)
        {
          auto bopn = dynamic_cast<BinaryOpNode *>(equations[eqn]->arg2);
          bopn->fillAutoregressiveRow(i++, lhs, AR);
        }
      ARr[model_name] = AR;
    }
  return ARr;
}

void
DynamicModel::fillTrendComponentModelTable() const
{
  map<string, vector<int>> eqnums, trend_eqnums, lhsr;
  map<string, vector<expr_t>> lhs_expr_tr;
  map<string, vector<set<pair<int, int>>>> rhsr;
  for (const auto &[model_name, eqtags] : trend_component_model_table.getTargetEqTags())
    {
      vector<int> trend_eqnumber;
      for (const auto &eqtag : eqtags)
        {
          optional<int> eqn { equation_tags.getEqnByTag("name", eqtag) };
          if (!eqn)
            {
              cerr << "ERROR: no equation is named '" << eqtag << "'" << endl;
              exit(EXIT_FAILURE);
            }
          trend_eqnumber.push_back(*eqn);
        }
      trend_eqnums[model_name] = trend_eqnumber;
    }

  for (const auto &[model_name, eqtags] : trend_component_model_table.getEqTags())
    {
      vector<int> eqnumber, lhs;
      vector<expr_t> lhs_expr_t;
      vector<set<pair<int, int>>> rhs;

      for (const auto &eqtag : eqtags)
        {
          set<pair<int, int>> lhs_set, lhs_tmp_set, rhs_set;
          optional<int> eqn { equation_tags.getEqnByTag("name", eqtag) };
          if (!eqn)
            {
              cerr << "ERROR: no equation is named '" << eqtag << "'" << endl;
              exit(EXIT_FAILURE);
            }

          equations[*eqn]->arg1->collectDynamicVariables(SymbolType::endogenous, lhs_set);
          equations[*eqn]->arg1->collectDynamicVariables(SymbolType::exogenous, lhs_tmp_set);
          equations[*eqn]->arg1->collectDynamicVariables(SymbolType::parameter, lhs_tmp_set);

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

          eqnumber.push_back(*eqn);
          lhs.push_back(itlhs->first);
          lhs_set.clear();
          set<expr_t> lhs_expr_t_set;
          equations[*eqn]->arg1->collectVARLHSVariable(lhs_expr_t_set);
          lhs_expr_t.push_back(*(lhs_expr_t_set.begin()));

          equations[*eqn]->arg2->collectDynamicVariables(SymbolType::endogenous, rhs_set);
          rhs.push_back(rhs_set);
        }
      eqnums[model_name] = eqnumber;
      lhsr[model_name] = lhs;
      lhs_expr_tr[model_name] = lhs_expr_t;
      rhsr[model_name] = rhs;
    }
  trend_component_model_table.setRhs(move(rhsr));
  trend_component_model_table.setVals(move(eqnums), move(trend_eqnums), move(lhsr), move(lhs_expr_tr));
}

pair<map<string, map<tuple<int, int>, expr_t>>, map<string, map<tuple<int, int>, expr_t>>>
DynamicModel::computeErrorComponentMatrices(const ExprNode::subst_table_t &diff_subst_table) const
{
  map<string, map<tuple<int, int>, expr_t>> A0r, A0starr;

  for (const auto &[model_name, eqns] : trend_component_model_table.getEqNums())
    {
      map<tuple<int, int>, expr_t> A0, A0star;
      vector<int> target_lhs = trend_component_model_table.getTargetLhs(model_name);
      vector<int> nontarget_eqnums = trend_component_model_table.getNonTargetEqNums(model_name);
      vector<int> undiff_nontarget_lhs = getUndiffLHSForPac(model_name, diff_subst_table);
      vector<int> parsed_undiff_nontarget_lhs;

      for (int i{0};
           auto eqn : eqns)
        {
          if (find(nontarget_eqnums.begin(), nontarget_eqnums.end(), eqn) != nontarget_eqnums.end())
            parsed_undiff_nontarget_lhs.push_back(undiff_nontarget_lhs.at(i));
          i++;
        }

      for (int i{0};
           auto eqn : eqns)
        if (find(nontarget_eqnums.begin(), nontarget_eqnums.end(), eqn) != nontarget_eqnums.end())
          equations[eqn]->arg2->fillErrorCorrectionRow(i++, parsed_undiff_nontarget_lhs, target_lhs, A0, A0star);
      A0r[model_name] = A0;
      A0starr[model_name] = A0star;
    }

  return { A0r, A0starr };
}

void
DynamicModel::fillTrendComponentModelTableFromOrigModel() const
{
  map<string, vector<int>> lags;
  map<string, vector<optional<int>>> orig_diff_var;
  map<string, vector<bool>> diff;
  for (const auto &[model_name, eqns] : trend_component_model_table.getEqNums())
    {
      set<expr_t> lhs;
      vector<optional<int>> orig_diff_var_vec;
      vector<bool> diff_vec;
      for (auto eqn : eqns)
        {
          // Perform some sanity checks on the RHS
          optional<string> eqtag { equation_tags.getTagValueByEqnAndKey(eqn, "name") };
          set<pair<int, int>> rhs_endo_set, rhs_exo_set;
          equations[eqn]->arg2->collectDynamicVariables(SymbolType::endogenous, rhs_endo_set);
          for (const auto &[symb_id, lag] : rhs_endo_set)
            if (lag >= 0)
              {
                cerr << "ERROR: in Equation " << eqtag.value_or(to_string(eqn+1))
                     << ". A trend component model may not have leaded or contemporaneous endogenous variables on the RHS. " << endl;
                exit(EXIT_FAILURE);
              }
          equations[eqn]->arg2->collectDynamicVariables(SymbolType::exogenous, rhs_exo_set);
          for (const auto &[symb_id, lag] : rhs_exo_set)
            if (lag != 0)
              {
                cerr << "ERROR: in Equation " << eqtag.value_or(to_string(eqn+1))
                     << ". A trend component model may not have lagged or leaded exogenous variables on the RHS. " << endl;
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
            orig_diff_var_vec.push_back(nullopt);
        }

      if (eqns.size() != lhs.size())
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
      for (auto eqn : eqns)
        max_lag.push_back(equations[eqn]->arg2->VarMaxLag(lhs_lag_equiv));
      lags[model_name] = max_lag;
      diff[model_name] = diff_vec;
      orig_diff_var[model_name] = orig_diff_var_vec;
    }
  trend_component_model_table.setDiff(move(diff));
  trend_component_model_table.setMaxLags(move(lags));
  trend_component_model_table.setOrigDiffVar(move(orig_diff_var));
}

void
DynamicModel::fillTrendComponentModelTableAREC(const ExprNode::subst_table_t &diff_subst_table) const
{
  auto ARr = computeAutoregressiveMatrices();
  trend_component_model_table.setAR(move(ARr));
  auto [A0r, A0starr] = computeErrorComponentMatrices(diff_subst_table);
  trend_component_model_table.setA0(move(A0r), move(A0starr));
}

vector<int>
DynamicModel::getUndiffLHSForPac(const string &aux_model_name,
                                 const ExprNode::subst_table_t &diff_subst_table) const
{
  vector<expr_t> lhs_expr_t = trend_component_model_table.getLhsExprT(aux_model_name);
  vector<int> lhs = trend_component_model_table.getLhs(aux_model_name);
  vector<bool> diff = trend_component_model_table.getDiff(aux_model_name);
  vector<int> eqnumber = trend_component_model_table.getEqNums(aux_model_name);
  vector<int> nontrend_eqnums = trend_component_model_table.getNonTargetEqNums(aux_model_name);

  for (auto eqn : nontrend_eqnums)
    {
      auto i = distance(eqnumber.begin(), find(eqnumber.begin(), eqnumber.end(), eqn));

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

void
DynamicModel::analyzePacEquationStructure(const string &name, map<string, string> &pac_eq_name, PacModelTable::equation_info_t &pac_equation_info)
{
  for (auto &equation : equations)
    if (equation->containsPacExpectation(name))
      {
        if (!pac_eq_name[name].empty())
          {
            cerr << "It is not possible to use 'pac_expectation(" << name << ")' in several equations." << endl;
            exit(EXIT_FAILURE);
          }
        optional<string> eqn { equation_tags.getTagValueByEqnAndKey(&equation - &equations[0], "name") };
        if (!eqn)
          {
            cerr << "Every equation with a 'pac_expectation' operator must have been assigned an equation tag name" << endl;
            exit(EXIT_FAILURE);
          }
        pac_eq_name[name] = *eqn;

        set<pair<int, int>> lhss;
        equation->arg1->collectDynamicVariables(SymbolType::endogenous, lhss);
        auto lhs = *lhss.begin();
        int lhs_symb_id = lhs.first;
        int lhs_orig_symb_id = lhs_symb_id;
        if (symbol_table.isDiffAuxiliaryVariable(lhs_orig_symb_id))
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
          = arg2->getPacOptimizingShareAndExprNodes(lhs_orig_symb_id);

        pair<int, vector<tuple<int, bool, int>>> ec_params_and_vars;
        vector<tuple<optional<int>, optional<int>, int>> ar_params_and_vars;
        vector<tuple<int, int, optional<int>, double>> non_optim_vars_params_and_constants, optim_additive_vars_params_and_constants, additive_vars_params_and_constants;
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

        if (lhs.first == -1)
          {
            cerr << "analyzePacEquationStructure: error obtaining LHS variable." << endl;
            exit(EXIT_FAILURE);
          }
        if (ec_params_and_vars.second.empty())
          {
            cerr << "analyzePacEquationStructure: error obtaining RHS parameters." << endl;
            exit(EXIT_FAILURE);
          }
        pac_equation_info[name] = { lhs, move(optim_share_index),
          move(ar_params_and_vars), move(ec_params_and_vars),
          move(non_optim_vars_params_and_constants),
          move(additive_vars_params_and_constants),
          move(optim_additive_vars_params_and_constants) };
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
          = barg2->getPacOptimizingShareAndExprNodes(undiff_lhs_symb_id);
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
DynamicModel::computePacModelConsistentExpectationSubstitution(const string &name,
                                                               int discount_symb_id,
                                                               int pac_eq_max_lag,
                                                               expr_t growth_correction_term,
                                                               string auxname,
                                                               ExprNode::subst_table_t &diff_subst_table,
                                                               map<string, int> &pac_aux_var_symb_ids,
                                                               map<string, vector<int>> &pac_aux_param_symb_ids,
                                                               map<string, expr_t> &pac_expectation_substitution)
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
  int neqs = 0;

  // Create the endogenous representing Z₁ (no orig_expr is given since its definition is recursive)
  if (auxname.empty())
    auxname = "mce_Z1_" + name;
  int mce_z1_symb_id = symbol_table.addPacExpectationAuxiliaryVar(auxname, nullptr);
  pac_aux_var_symb_ids[name] = mce_z1_symb_id;

  expr_t A = One;
  expr_t fp = Zero;
  expr_t beta = AddVariable(discount_symb_id);
  for (int i = 1; i <= pac_eq_max_lag+1; i++)
    {
      string param_name = "mce_alpha_" + name + "_" + to_string(i);
      try
        {
          int alpha_i_symb_id = symbol_table.addSymbol(param_name, SymbolType::parameter);
          pac_aux_param_symb_ids[name].push_back(alpha_i_symb_id);
          A = AddPlus(A, AddVariable(alpha_i_symb_id));
          fp = AddPlus(fp,
                       AddTimes(AddTimes(AddVariable(alpha_i_symb_id),
                                         AddPower(beta, AddPossiblyNegativeConstant(i))),
                                AddVariable(mce_z1_symb_id, i)));
        }
      catch (SymbolTable::AlreadyDeclaredException &e)
        {
          cerr << "The variable/parameter '" << param_name << "' conflicts with a parameter that will be generated for the '" << name << "' PAC model. Please rename it." << endl;
          exit(EXIT_FAILURE);
        }
    }

  // Add diff nodes and eqs for pac_target_symb_id
  const VariableNode *target_base_diff_node;
  auto create_target_lag = [&](int lag)
  {
    if (symbol_table.isAuxiliaryVariable(pac_target_symb_id))
      // We know it is a log, see ExprNode::matchParamTimesTargetMinusVariable()
      return AddLog(AddVariable(symbol_table.getOrigSymbIdForAuxVar(pac_target_symb_id), lag));
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
      addEquation(neweq, nullopt);
      addAuxEquation(neweq);
      neqs++;
    }

  map<int, VariableNode *> target_aux_var_to_add;
  const VariableNode *last_aux_var = target_base_diff_node;
  for (int i = 1; i <= pac_eq_max_lag; i++, neqs++)
    {
      expr_t this_diff_node = AddDiff(create_target_lag(i));
      int symb_id = symbol_table.addDiffLeadAuxiliaryVar(this_diff_node->idx, this_diff_node,
                                                         last_aux_var->symb_id, 1);
      VariableNode *current_aux_var = AddVariable(symb_id);
      auto neweq = AddEqual(current_aux_var, AddVariable(last_aux_var->symb_id, 1));
      addEquation(neweq, nullopt);
      addAuxEquation(neweq);
      last_aux_var = current_aux_var;
      target_aux_var_to_add[i] = current_aux_var;
    }

  expr_t fs = Zero;
  for (int k = 1; k <= pac_eq_max_lag; k++)
    {
      expr_t ssum = Zero;
      for (int j = k+1; j <= pac_eq_max_lag+1; j++)
        {
          int alpha_j_symb_id = -1;
          string param_name = "mce_alpha_" + name + "_" + to_string(j);
          try
            {
              alpha_j_symb_id = symbol_table.getID(param_name);
            }
          catch (SymbolTable::UnknownSymbolNameException &e)
            {
              alpha_j_symb_id = symbol_table.addSymbol(param_name, SymbolType::parameter);
            }
          ssum = AddPlus(ssum,
                         AddTimes(AddVariable(alpha_j_symb_id), AddPower(beta, AddPossiblyNegativeConstant(j))));
        }
      fs = AddPlus(fs, AddTimes(ssum, target_aux_var_to_add[k]));
    }
  auto neweq = AddEqual(AddVariable(mce_z1_symb_id),
                        AddMinus(AddTimes(A, AddMinus(const_cast<VariableNode *>(target_base_diff_node), fs)), fp));
  addEquation(neweq, nullopt);
  /* This equation is not added to the list of auxiliary equations, because it
     is recursive, and this would in particular break dynamic_set_auxiliary_series.m */
  neqs++;

  cout << "PAC Model Consistent Expectation: added " << neqs << " auxiliary variables and equations for model " << name << "." << endl;

  /* The growth correction term is not added to the definition of Z₁
     because the latter is recursive. Rather put it at the level of the
     substition of pac_expectation operator. */
  pac_expectation_substitution[name] = AddPlus(AddVariable(mce_z1_symb_id), growth_correction_term);
}

void
DynamicModel::computePacBackwardExpectationSubstitution(const string &name,
                                                        const vector<int> &lhs,
                                                        int max_lag,
                                                        const string &aux_model_type,
                                                        expr_t growth_correction_term,
                                                        string auxname,
                                                        map<string, int> &pac_aux_var_symb_ids,
                                                        map<string, vector<int>> &pac_aux_param_symb_ids,
                                                        map<string, expr_t> &pac_expectation_substitution)
{
  auto create_aux_param = [&](const string &param_name)
  {
    try
      {
        return symbol_table.addSymbol(param_name, SymbolType::parameter);
      }
    catch (SymbolTable::AlreadyDeclaredException)
      {
        cerr << "ERROR: the variable/parameter '" << param_name << "' conflicts with some auxiliary parameter that will be generated for the '" << name << "' PAC model. Please rename that parameter." << endl;
        exit(EXIT_FAILURE);
      }
  };

  expr_t subExpr = Zero;
  if (aux_model_type == "var")
    {
      /* If the auxiliary model is a VAR, add a parameter corresponding
         to the constant. */
      int new_param_symb_id = create_aux_param("h_" + name + "_constant");
      pac_aux_param_symb_ids[name].push_back(new_param_symb_id);
      subExpr = AddPlus(subExpr, AddVariable(new_param_symb_id));
    }
  for (int i = 1; i < max_lag + 1; i++)
    for (auto lhsit : lhs)
      {
        int new_param_symb_id = create_aux_param("h_" + name + "_var_"
                                                 + symbol_table.getName(lhsit)
                                                 + "_lag_" + to_string(i));
        pac_aux_param_symb_ids[name].push_back(new_param_symb_id);
        subExpr = AddPlus(subExpr,
                          AddTimes(AddVariable(new_param_symb_id),
                                   AddVariable(lhsit, -i)));
      }

  subExpr = AddPlus(subExpr, growth_correction_term);

  if (auxname.empty())
    auxname = "pac_expectation_" + name;
  int expect_var_id = symbol_table.addPacExpectationAuxiliaryVar(auxname, subExpr);
  expr_t neweq = AddEqual(AddVariable(expect_var_id), subExpr);
  addEquation(neweq, nullopt);
  addAuxEquation(neweq);
  pac_aux_var_symb_ids[name] = expect_var_id;
  pac_expectation_substitution[name] = AddVariable(expect_var_id);
}

void
DynamicModel::computePacBackwardExpectationSubstitutionWithComponents(const string &name,
                                                                      const vector<int> &lhs,
                                                                      int max_lag,
                                                                      const string &aux_model_type,
                                                                      vector<PacModelTable::target_component_t> &pac_target_components,
                                                                      map<string, expr_t> &pac_expectation_substitution)
{
  auto create_aux_param = [&](const string &param_name)
  {
    try
      {
        return symbol_table.addSymbol(param_name, SymbolType::parameter);
      }
    catch (SymbolTable::AlreadyDeclaredException)
      {
        cerr << "ERROR: the variable/parameter '" << param_name << "' conflicts with some auxiliary parameter that will be generated for the '" << name << "' PAC model. Please rename that parameter." << endl;
        exit(EXIT_FAILURE);
      }
  };

  expr_t substexpr = Zero;

  for (int component_idx{1};
       auto &[component, growth, auxname, kind, coeff, growth_neutrality_param, h_indices, original_growth, growth_info] : pac_target_components)
    {
      string name_component = name + "_component" + to_string(component_idx);

      // Create the linear combination of the variables from the auxiliary model
      expr_t auxdef = Zero;
      if (aux_model_type == "var")
        {
          /* If the auxiliary model is a VAR, add a parameter corresponding
             to the constant. */
          int new_param_symb_id = create_aux_param("h_" + name_component + "_constant");
          h_indices.push_back(new_param_symb_id);
          auxdef = AddPlus(auxdef, AddVariable(new_param_symb_id));
        }
      for (int i = 1; i < max_lag + 1; i++)
        for (auto lhsit : lhs)
          {
            int new_param_symb_id = create_aux_param("h_" + name_component + "_var_" + symbol_table.getName(lhsit) + "_lag_" + to_string(i));
            h_indices.push_back(new_param_symb_id);
            auxdef = AddPlus(auxdef, AddTimes(AddVariable(new_param_symb_id),
                                              AddVariable(lhsit, -i)));
          }

      // If needed, add the growth neutrality correction for this component
      if (growth)
        {
          growth_neutrality_param = create_aux_param(name_component + "_pac_growth_neutrality_correction");
          auxdef = AddPlus(auxdef, AddTimes(growth, AddVariable(growth_neutrality_param)));
        }
      else
        growth_neutrality_param = -1;

      // Create the auxiliary variable for this component
      int aux_id = symbol_table.addPacExpectationAuxiliaryVar(auxname, auxdef);
      expr_t auxvar = AddVariable(aux_id);

      // Add the equation defining the auxiliary variable for this component
      expr_t neweq = AddEqual(auxvar, auxdef);
      addEquation(neweq, nullopt);
      addAuxEquation(neweq);

      // Update the expression to be substituted for the pac_expectation operator
      substexpr = AddPlus(substexpr, AddTimes(coeff, auxvar));

      component_idx++;
    }

  pac_expectation_substitution[name] = substexpr;
}

void
DynamicModel::substitutePacExpectation(const map<string, expr_t> &pac_expectation_substitution,
                                       const map<string, string> &pac_eq_name)
{
  for (auto &[model_name, substexpr] : pac_expectation_substitution)
    {
      optional<int> eq { equation_tags.getEqnByTag("name", pac_eq_name.at(model_name)) };
      auto substeq = dynamic_cast<BinaryOpNode *>(equations[eq.value()]->substitutePacExpectation(model_name, substexpr));
      assert(substeq);
      equations[eq.value()] = substeq;
    }
}

void
DynamicModel::substitutePacTargetNonstationary(const string &pac_model_name, expr_t substexpr)
{
  for (auto &eq : equations)
    eq = dynamic_cast<BinaryOpNode *>(eq->substitutePacTargetNonstationary(pac_model_name, substexpr));
}

void
DynamicModel::computingPass(int derivsOrder, int paramsDerivsOrder, const eval_context_t &eval_context,
                            bool no_tmp_terms, bool block, bool use_dll)
{
  initializeVariablesAndEquations();

  // Prepare for derivation
  computeDerivIDs();

  // Computes dynamic jacobian columns, must be done after computeDerivIDs()
  computeDynJacobianCols();

  /* In both MATLAB and Julia, tensors for higher-order derivatives are stored
     in matrices whose columns correspond to variable multi-indices. Since we
     currently are limited to 32-bit signed integers (hence 31 bits) for matrix
     indices, check that we will not overflow (see #89). Note that such a check
     is not needed for parameter derivatives, since tensors for those are not
     stored as matrices. This check cannot be done before since
     getJacobianColsNbr() is not yet set.
     We only do the check for the legacy representation, since the sparse
     representation is not affected by this problem (TODO: thus the check can be
     removed once the legacy representation is dropped). */
  if (log2(getJacobianColsNbr(false))*derivsOrder >= numeric_limits<int>::digits)
    {
      cerr << "ERROR: The derivatives matrix of the " << modelClassName() << " is too large. Please decrease the approximation order." << endl;
      exit(EXIT_FAILURE);
    }

  // Compute derivatives w.r. to all endogenous, exogenous and exogenous deterministic
  set<int> vars;
  for (auto &it : deriv_id_table)
    {
      SymbolType type = symbol_table.getType(it.first.first);
      if (type == SymbolType::endogenous || type == SymbolType::exogenous
          || type == SymbolType::exogenousDet)
        vars.insert(it.second);
    }

  // Launch computations
  cout << "Computing " << modelClassName() << " derivatives (order " << derivsOrder << ")." << endl;

  computeDerivatives(derivsOrder, vars);

  if (derivsOrder > 1)
    for (const auto &[indices, d2] : derivatives[2])
      nonzero_hessian_eqs.insert(indices[0]);

  if (paramsDerivsOrder > 0)
    {
      cout << "Computing " << modelClassName() << " derivatives w.r.t. parameters (order " << paramsDerivsOrder << ")." << endl;
      computeParamsDerivatives(paramsDerivsOrder);
    }

  computeTemporaryTerms(!use_dll, no_tmp_terms);

  if (paramsDerivsOrder > 0 && !no_tmp_terms)
    computeParamsDerivativesTemporaryTerms();

  if (!block && (max_endo_lag == 0 || max_endo_lead == 0))
    {
      time_recursive_block_decomposition = true;
      mfs = 3; // FIXME: remove this line when mfs=3 becomes the global default
    }
  computingPassBlock(eval_context, no_tmp_terms);
  if (block_decomposed)
    computeBlockDynJacobianCols();
  if (!block_decomposed && block)
    {
      cerr << "ERROR: Block decomposition requested but failed." << endl;
      exit(EXIT_FAILURE);
    }
}

void
DynamicModel::computeXrefs()
{
  for (int i{0};
       auto &equation : equations)
    {
      ExprNode::EquationInfo ei;
      equation->computeXrefs(ei);
      xrefs[i++] = ei;
    }

  for (int i{0};
       const auto &[eq, eqinfo] : xrefs)
    {
      computeRevXref(xref_param, eqinfo.param, i);
      computeRevXref(xref_endo, eqinfo.endo, i);
      computeRevXref(xref_exo, eqinfo.exo, i);
      computeRevXref(xref_exo_det, eqinfo.exo_det, i);
      i++;
    }
}

void
DynamicModel::computeRevXref(map<pair<int, int>, set<int>> &xrefset, const set<pair<int, int>> &eiref, int eqn)
{
  for (const auto &it : eiref)
    {
      set<int> eq;
      if (xrefset.contains(it))
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
  for (int i{1};
       const auto &[eq, eqinfo] : xrefs)
    {
      output << "M_.xref1.param{" << i << "} = [ ";
      for (const auto &[id, lag] : eqinfo.param)
        output << symbol_table.getTypeSpecificID(id) + 1 << " ";
      output << "];" << endl;

      output << "M_.xref1.endo{" << i << "} = [ ";
      for (const auto &[id, lag] : eqinfo.endo)
        output << "struct('id', " << symbol_table.getTypeSpecificID(id) + 1 << ", 'shift', " << lag << ");";
      output << "];" << endl;

      output << "M_.xref1.exo{" << i << "} = [ ";
      for (const auto &[id, lag] : eqinfo.exo)
        output << "struct('id', " << symbol_table.getTypeSpecificID(id) + 1 << ", 'shift', " << lag << ");";
      output << "];" << endl;

      output << "M_.xref1.exo_det{" << i << "} = [ ";
      for (const auto &[id, lag] : eqinfo.exo_det)
        output << "struct('id', " << symbol_table.getTypeSpecificID(id) + 1 << ", 'shift', " << lag << ");";
      output << "];" << endl;

      i++;
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
  for (int last_tsid{-1};
       const auto &[key, eqs] : xrefmap)
    {
      auto &[id, lag] = key;
      int tsid = symbol_table.getTypeSpecificID(id) + 1;
      output << "M_.xref2." << type << "{" << tsid << "} = [ ";
      if (last_tsid == tsid)
        output << "M_.xref2." << type << "{" << tsid << "}; ";
      else
        last_tsid = tsid;

      for (int eq : eqs)
        if (type == "param")
          output << eq + 1 << " ";
        else
          output << "struct('shift', " << lag << ", 'eq', " << eq+1 << ");";
      output << "];" << endl;
    }
}


map<tuple<int, int, int>, DynamicModel::BlockDerivativeType>
DynamicModel::determineBlockDerivativesType(int blk)
{
  map<tuple<int, int, int>, BlockDerivativeType> derivType;
  int size = blocks[blk].size;
  int nb_recursive = blocks[blk].getRecursiveSize();
  for (int lag {time_recursive_block_decomposition ? 0 : -blocks[blk].max_endo_lag};
       lag <= (time_recursive_block_decomposition ? 0 : blocks[blk].max_endo_lead);
       lag++)
    for (int eq = 0; eq < size; eq++)
      {
        set<pair<int, int>> endos_and_lags;
        int eq_orig = getBlockEquationID(blk, eq);
        equations[eq_orig]->collectEndogenous(endos_and_lags);
        for (int var = 0; var < size; var++)
          if (int var_orig = getBlockVariableID(blk, var);
              endos_and_lags.contains({ var_orig, lag }))
            {
              if (getBlockEquationType(blk, eq) == EquationType::evaluateRenormalized
                  && eq < nb_recursive)
                /* It’s a normalized recursive equation, we have to recompute
                   the derivative using the chain rule */
                derivType[{ lag, eq, var }] = BlockDerivativeType::normalizedChainRule;
              else if (!derivType.contains({ lag, eq, var }))
                derivType[{ lag, eq, var }] = BlockDerivativeType::standard;

              if (var < nb_recursive)
                for (int feedback_var = nb_recursive; feedback_var < size; feedback_var++)
                  if (derivType.contains({ lag, var, feedback_var }))
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
  size_t nb_blocks { blocks.size() };

  blocks_derivatives.resize(nb_blocks);
  blocks_jacobian_sparse_column_major_order.resize(nb_blocks);
  blocks_jacobian_sparse_colptr.resize(nb_blocks);

  for (int blk {0}; blk < static_cast<int>(nb_blocks); blk++)
    {
      int nb_recursives = blocks[blk].getRecursiveSize();
      int mfs_size { blocks[blk].mfs_size };
      BlockSimulationType simulation_type { blocks[blk].simulation_type };

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
      map<pair<expr_t, int>, expr_t> chain_rule_deriv_cache;
      for (const auto &[indices, derivType] : determineBlockDerivativesType(blk))
        {
          auto [lag, eq, var] = indices;
          int eq_orig = getBlockEquationID(blk, eq), var_orig = getBlockVariableID(blk, var);
          int deriv_id = getDerivID(symbol_table.getID(SymbolType::endogenous, var_orig), lag);
          expr_t d{nullptr};
          switch (derivType)
            {
            case BlockDerivativeType::standard:
              if (auto it = derivatives[1].find({ eq_orig, deriv_id });
                  it != derivatives[1].end())
                d = it->second;
              else
                d = Zero;
              break;
            case BlockDerivativeType::chainRule:
              d = equations[eq_orig]->getChainRuleDerivative(deriv_id, recursive_vars, chain_rule_deriv_cache);
              break;
            case BlockDerivativeType::normalizedChainRule:
              d = equation_type_and_normalized_equation[eq_orig].second->getChainRuleDerivative(deriv_id, recursive_vars, chain_rule_deriv_cache);
              break;
            }

          if (d != Zero)
            blocks_derivatives[blk][{ eq, var, lag }] = d;
        }

      // Compute the sparse representation of the Jacobian
      if (simulation_type != BlockSimulationType::evaluateForward
          && simulation_type != BlockSimulationType::evaluateBackward)
        {
          const bool one_boundary {simulation_type == BlockSimulationType::solveBackwardSimple
                                   || simulation_type == BlockSimulationType::solveForwardSimple
                                   || simulation_type == BlockSimulationType::solveBackwardComplete
                                   || simulation_type == BlockSimulationType::solveForwardComplete};
          for (const auto &[indices, d1] : blocks_derivatives[blk])
            {
              auto &[eq, var, lag] { indices };
              assert(lag >= -1 && lag <= 1);
              if (eq >= nb_recursives && var >= nb_recursives
                  && !(one_boundary && lag != 0))
                blocks_jacobian_sparse_column_major_order[blk].emplace(pair{eq-nb_recursives, var-nb_recursives+static_cast<int>(!one_boundary)*(lag+1)*mfs_size}, d1);
            }
          blocks_jacobian_sparse_colptr[blk] = computeCSCColPtr(blocks_jacobian_sparse_column_major_order[blk], (one_boundary ? 1 : 3)*mfs_size);
        }
    }
}

void
DynamicModel::computeBlockDynJacobianCols()
{
  size_t nb_blocks { blocks.size() };
  // Structures used for lexicographic ordering over (lag, var ID)
  vector<set<pair<int, int>>> dynamic_endo(nb_blocks), dynamic_other_endo(nb_blocks),
    dynamic_exo(nb_blocks), dynamic_exo_det(nb_blocks);

  for (auto & [indices, d1] : derivatives[1])
    {
      auto [eq_orig, deriv_id] { vectorToTuple<2>(indices) };
      int block_eq { eq2block[eq_orig] };
      int var { getTypeSpecificIDByDerivID(deriv_id) };
      int lag { getLagByDerivID(deriv_id) };
      switch (getTypeByDerivID(deriv_id))
        {
        case SymbolType::endogenous:
          if (block_eq == endo2block[var])
            dynamic_endo[block_eq].emplace(lag, getBlockInitialVariableID(block_eq, var));
          else
            dynamic_other_endo[block_eq].emplace(lag, var);
          break;
        case SymbolType::exogenous:
          dynamic_exo[block_eq].emplace(lag, var);
          break;
        case SymbolType::exogenousDet:
          dynamic_exo_det[block_eq].emplace(lag, var);
          break;
        default:
          break;
        }
    }

  // Compute Jacobian column indices
  blocks_jacob_cols_endo.resize(nb_blocks);
  for (size_t blk {0}; blk < nb_blocks; blk++)
    for (int index{0};
         auto [lag, var] : dynamic_endo[blk])
      blocks_jacob_cols_endo[blk][{ var, lag }] = index++;
}

void
DynamicModel::writeDynamicFile(const string &basename, bool use_dll, const string &mexext, const filesystem::path &matlabroot, const filesystem::path &dynareroot, bool julia) const
{
  filesystem::path model_dir{basename};
  model_dir /= "model";
  if (use_dll)
    {
      create_directories(model_dir / "src" / "sparse");
      if (block_decomposed)
        create_directories(model_dir / "src" / "sparse" / "block");
    }
  if (julia)
    create_directories(model_dir / "julia");
  else
    {
      auto plusfolder {packageDir(basename)};
      /* The following is not a duplicate of the same call from
         ModFile::writeMOutput(), because of planner_objective which needs its
         +objective subdirectory */
      create_directories(plusfolder);

      auto sparsefolder {plusfolder / "+sparse"};
      create_directories(sparsefolder);
      if (block_decomposed)
        create_directories(sparsefolder / "+block");

      create_directories(plusfolder / "+debug");
    }
  create_directories(model_dir / "bytecode" / "block");

  // Legacy representation
  if (use_dll)
    writeModelCFile<true>(basename, mexext, matlabroot, dynareroot);
  else if (!julia) // M-files
    writeDynamicMFile(basename);
  // The legacy representation is no longer produced for Julia

  writeDynamicBytecode(basename);
  if (block_decomposed)
    writeDynamicBlockBytecode(basename);

  // Sparse representation
  if (use_dll)
    writeSparseModelCFiles<true>(basename, mexext, matlabroot, dynareroot);
  else if (julia)
    writeSparseModelJuliaFiles<true>(basename);
  else // MATLAB/Octave
    writeSparseModelMFiles<true>(basename);

  writeSetAuxiliaryVariables(basename, julia);

  // Support for model debugging
  if (!julia)
    writeDebugModelMFiles<true>(basename);
}

void
DynamicModel::writeSetAuxiliaryVariables(const string &basename, bool julia) const
{
  ostringstream output_func_body;
  writeAuxVarRecursiveDefinitions(output_func_body, julia ? ExprNodeOutputType::juliaTimeDataFrame
                                  : ExprNodeOutputType::matlabDseries);

  if (output_func_body.str().empty())
    return;

  string func_name = julia ? "dynamic_set_auxiliary_series!" : "dynamic_set_auxiliary_series";
  string comment = julia ? "#" : "%";

  stringstream output;
  output << "function ";
  if (!julia)
    output << "ds = ";
  output << func_name + "(ds, params)" << endl
         << comment << endl
         << comment << " Status : Computes Auxiliary variables of the " << modelClassName() << " and returns a dseries" << endl
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
    writeToFileIfModified(output, filesystem::path{basename} / "model" / "julia" / "DynamicSetAuxiliarySeries.jl");
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
  addEquation(static_model.equations[0]->clone(*this), nullopt);

  // Get max endo lead and max endo lag
  set<pair<int, int>> dynvars;
  int max_eq_lead = 0;
  int max_eq_lag = 0;
  for (auto &equation : equations)
    equation->collectDynamicVariables(SymbolType::endogenous, dynvars);

  for (const auto &[symb_id, lag] : dynvars)
    {
      max_eq_lead = max(lag, max_eq_lead);
      max_eq_lag = max(-lag, max_eq_lag);
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
  addEquation(AddEqual(lagrangian, Zero), nullopt);
  computeDerivIDs();

  /* Compute Lagrangian derivatives.
     Also restore line numbers and tags for FOCs w.r.t. a Lagrange multiplier
     (i.e. a FOC identical to an equation of the original model) */
  vector<expr_t> neweqs;
  vector<optional<int>> neweqs_lineno;
  map<int, map<string, string>> neweqs_tags;
  for (auto &[symb_id_and_lag, deriv_id] : deriv_id_table)
    {
      auto &[symb_id, lag] = symb_id_and_lag;
      if (symbol_table.getType(symb_id) == SymbolType::endogenous && lag == 0)
        {
          neweqs.push_back(AddEqual(equations[0]->getNonZeroPartofEquation()->getDerivative(deriv_id), Zero));
          if (optional<int> i = symbol_table.getEquationNumberForMultiplier(symb_id);
              i)
            {
              // This is a derivative w.r.t. a Lagrange multiplier
              neweqs_lineno.push_back(old_equations_lineno[*i]);
              neweqs_tags[neweqs.size()-1] = old_equation_tags.getTagsByEqn(*i);
            }
          else
            neweqs_lineno.push_back(nullopt);
        }
    }

  // Overwrite equations with the Lagrangian derivatives
  clearEquations();
  for (size_t i = 0; i < neweqs.size(); i++)
    addEquation(neweqs[i], neweqs_lineno[i], neweqs_tags[i]);
}

bool
DynamicModel::ParamUsedWithLeadLag() const
{
  return ParamUsedWithLeadLagInternal();
}

void
DynamicModel::createVariableMapping()
{
  for (size_t ii = 0; ii < equations.size(); ii++)
    {
      set<int> eqvars;
      equations[ii]->collectVariables(SymbolType::endogenous, eqvars);
      equations[ii]->collectVariables(SymbolType::exogenous, eqvars);
      for (auto eqvar : eqvars)
        variableMapping[symbol_table.getUltimateOrigSymbID(eqvar)].emplace(ii);
    }
}

void
DynamicModel::expandEqTags()
{
  set<int> existing_tags = equation_tags.getEqnsByKey("name");
  for (int eq = 0; eq < static_cast<int>(equations.size()); eq++)
    if (!existing_tags.contains(eq))
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
  for (auto &equation : static_only_equations)
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
  for (auto &equation : static_only_equations)
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

  for (const auto &[symb_id, lag] : dynvars)
    {
      SymbolType type = symbol_table.getType(symb_id);

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
    {
      equation->collectDynamicVariables(SymbolType::endogenous, dynvars);
      equation->collectDynamicVariables(SymbolType::exogenous, dynvars);
      equation->collectDynamicVariables(SymbolType::exogenousDet, dynvars);
      equation->collectDynamicVariables(SymbolType::parameter, dynvars);
      equation->collectDynamicVariables(SymbolType::trend, dynvars);
      equation->collectDynamicVariables(SymbolType::logTrend, dynvars);
    }

  for (const auto &[symb_id, lag] : dynvars)
    {
      SymbolType type = symbol_table.getType(symb_id);

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

      deriv_id_table[{symb_id, lag}] = deriv_id;
      inv_deriv_id_table.emplace_back(symb_id, lag);
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
DynamicModel::getTypeSpecificIDByDerivID(int deriv_id) const
{
  return symbol_table.getTypeSpecificID(getSymbIDByDerivID(deriv_id));
}

int
DynamicModel::getDerivID(int symb_id, int lag) const noexcept(false)
{
  if (auto it = deriv_id_table.find({ symb_id, lag });
      it == deriv_id_table.end())
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
DynamicModel::computeDynJacobianCols()
{
  // Sort the dynamic endogenous variables by lexicographic order over (lag, type_specific_symbol_id)
  map<pair<int, int>, int> ordered_dyn_endo;
  for (const auto &[symb_lag, deriv_id] : deriv_id_table)
    if (const auto &[symb_id, lag] = symb_lag;
        symbol_table.getType(symb_id) == SymbolType::endogenous)
      ordered_dyn_endo[{ lag, symbol_table.getTypeSpecificID(symb_id) }] = deriv_id;

  // Fill the dynamic jacobian columns for endogenous (legacy representation)
  for (int sorted_id{0};
       const auto &[ignore, deriv_id] : ordered_dyn_endo)
    dyn_jacobian_cols_table[deriv_id] = sorted_id++;

  /* Fill the dynamic columns for exogenous and exogenous deterministic (legacy
     representation) */
  for (const auto &[symb_lag, deriv_id] : deriv_id_table)
    {
      int symb_id{symb_lag.first};
      int tsid{symbol_table.getTypeSpecificID(symb_id)};
      if (SymbolType type{symbol_table.getType(symb_id)};
          type == SymbolType::exogenous)
        dyn_jacobian_cols_table[deriv_id] = ordered_dyn_endo.size() + tsid;
      else if (type == SymbolType::exogenousDet)
        dyn_jacobian_cols_table[deriv_id] = ordered_dyn_endo.size() + symbol_table.exo_nbr() + tsid;
    }

  /* NB: the following could differ from dyn_jacobian_cols_table.size() if
     there are unused exogenous (and “nostrict” option is given) */
  dyn_jacobian_ncols = ordered_dyn_endo.size() + symbol_table.exo_nbr() + symbol_table.exo_det_nbr();
}

void
DynamicModel::testTrendDerivativesEqualToZero(const eval_context_t &eval_context)
{
  for (auto &[symb_lag1, deriv_id1] : deriv_id_table)
    if (auto &[symb_id1, lag1] = symb_lag1;
        symbol_table.getType(symb_id1) == SymbolType::trend
        || symbol_table.getType(symb_id1) == SymbolType::logTrend)
      for (int eq = 0; eq < static_cast<int>(equations.size()); eq++)
        {
          expr_t homogeneq = AddMinus(equations[eq]->arg1,
                                      equations[eq]->arg2);

          // Do not run the test if the term inside the log is zero
          if (fabs(homogeneq->eval(eval_context)) > zero_band)
            {
              expr_t testeq = AddLog(homogeneq); // F = log(lhs-rhs)
              testeq = testeq->getDerivative(deriv_id1); // d F / d Trend
              for (auto &[symb_lag2, deriv_id2] : deriv_id_table)
                if (auto &[symb_id2, lag2] = symb_lag2;
                    symbol_table.getType(symb_id2) == SymbolType::endogenous)
                  {
                    double nearZero = testeq->getDerivative(deriv_id2)->eval(eval_context); // eval d F / d Trend d Endog
                    if (fabs(nearZero) > balanced_growth_test_tol)
                      {
                        cerr << "ERROR: trends not compatible with balanced growth path; the second-order cross partial of equation " << eq + 1;
                        if (equations_lineno[eq])
                          cerr << " (line " << *equations_lineno[eq] << ") ";
                        cerr << "w.r.t. trend variable "
                             << symbol_table.getName(symb_id1) << " and endogenous variable "
                             << symbol_table.getName(symb_id2) << " is not null (abs. value = "
                             << fabs(nearZero) << "). If you are confident that your trends are correctly specified, you can raise the value of option 'balanced_growth_test_tol' in the 'model' block." << endl;
                        exit(EXIT_FAILURE);
                      }
                  }
            }
        }
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
      const expr_t value = local_variables_table.at(used_local_var);
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
      addEquation(neweq, nullopt);
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
  /* Contrary to other substitution methods, we do the substitution in MLV
     definitions here, instead of doing it at the ExprNode method level,
     because otherwise this would substitute MLV in the original model (see
     #65). */
  for (auto &[id, definition] : local_variables_table)
    definition = definition->substituteAdl();

  for (auto &equation : equations)
    equation = dynamic_cast<BinaryOpNode *>(equation->substituteAdl());

  for (auto &equation : static_only_equations)
    equation = dynamic_cast<BinaryOpNode *>(equation->substituteAdl());
}

void
DynamicModel::substituteModelLocalVariables()
{
  for (auto &equation : equations)
    equation = dynamic_cast<BinaryOpNode *>(equation->substituteModelLocalVariables());

  for (auto &equation : static_only_equations)
    equation = dynamic_cast<BinaryOpNode *>(equation->substituteModelLocalVariables());

  /* We can’t clear local_variables_table at this point, because in case of
     ramsey_policy, the original model is saved via DynamicModel::operator=()
     before computing the FOC. But since DataTree::operator=() clones all
     nodes, it will try to clone nodes containing model-local variables, and
     this will fail at the point where DataTree methods try to evaluate those
     nodes to a numerical value. */
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

set<int>
DynamicModel::findPacExpectationEquationNumbers() const
{
  set<int> eqnumbers;
  for (int i{0};
       auto &equation : equations)
    {
      if (equation->containsPacExpectation())
        eqnumbers.insert(i);
      i++;
    }
  return eqnumbers;
}

pair<lag_equivalence_table_t, ExprNode::subst_table_t>
DynamicModel::substituteUnaryOps(VarExpectationModelTable &var_expectation_model_table, PacModelTable &pac_model_table)
{
  vector<int> eqnumbers(equations.size());
  iota(eqnumbers.begin(), eqnumbers.end(), 0);
  return substituteUnaryOps({ eqnumbers.begin(), eqnumbers.end() }, var_expectation_model_table, pac_model_table);
}

pair<lag_equivalence_table_t, ExprNode::subst_table_t>
DynamicModel::substituteUnaryOps(const set<int> &eqnumbers, VarExpectationModelTable &var_expectation_model_table, PacModelTable &pac_model_table)
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

  // Substitute in expressions of var_expectation_model
  var_expectation_model_table.substituteUnaryOpsInExpression(nodes, subst_table, neweqs);
  // Substitute in growth terms in pac_model and pac_target_info
  pac_model_table.substituteUnaryOpsInGrowth(nodes, subst_table, neweqs);

  // Add new equations
  for (auto &neweq : neweqs)
    {
      addEquation(neweq, nullopt);
      aux_equations.push_back(neweq);
    }

  if (subst_table.size() > 0)
    cout << "Substitution of Unary Ops: added " << neweqs.size() << " auxiliary variables and equations." << endl;

  return { nodes, subst_table };
}

pair<lag_equivalence_table_t, ExprNode::subst_table_t>
DynamicModel::substituteDiff(VarExpectationModelTable &var_expectation_model_table, PacModelTable &pac_model_table)
{
  /* Note: at this point, we know that there is no diff operator with a lead,
     because they have been expanded by DataTree::AddDiff().
     Hence we can go forward with the substitution without worrying about the
     expectation operator. */

  lag_equivalence_table_t diff_nodes;
  ExprNode::subst_table_t diff_subst_table;

  // Mark diff operators to be substituted in model local variables
  set<int> used_local_vars;
  for (auto equation : equations)
    equation->collectVariables(SymbolType::modelLocalVariable, used_local_vars);
  for (auto &[symb_id, expr] : local_variables_table)
    if (used_local_vars.contains(symb_id))
      expr->findDiffNodes(diff_nodes);

  // Mark diff operators to be substituted in equations
  for (auto equation : equations)
    equation->findDiffNodes(diff_nodes);

  pac_model_table.findDiffNodesInGrowth(diff_nodes);

  // Substitute in model local variables
  vector<BinaryOpNode *> neweqs;
  for (auto &[symb_id, expr] : local_variables_table)
    expr = expr->substituteDiff(diff_nodes, diff_subst_table, neweqs);

  // Substitute in equations
  for (auto &equation : equations)
    {
      auto substeq = dynamic_cast<BinaryOpNode *>(equation->
                                                  substituteDiff(diff_nodes, diff_subst_table, neweqs));
      assert(substeq);
      equation = substeq;
    }

  var_expectation_model_table.substituteDiffNodesInExpression(diff_nodes, diff_subst_table, neweqs);
  pac_model_table.substituteDiffNodesInGrowth(diff_nodes, diff_subst_table, neweqs);

  // Add new equations
  for (auto neweq : neweqs)
    {
      addEquation(neweq, nullopt);
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
  for (auto &[symb_id, expr] : local_variables_table)
    expr = expr->substituteExpectation(subst_table, neweqs, partial_information_model);

  // Substitute in equations
  for (auto &equation : equations)
    {
      equation = dynamic_cast<BinaryOpNode *>(equation->substituteExpectation(subst_table, neweqs, partial_information_model));
      assert(equation);
    }

  /* No need to substitute in static_only_equations, since expectation()
     operators in [static] equations are forbidden at the parsing level. */

  // Add new equations
  for (auto neweq : neweqs)
    {
      addEquation(neweq, nullopt);
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
  for (auto &[id, definition] : local_variables_table)
    definition = definition->decreaseLeadsLagsPredeterminedVariables();

  for (auto &equation : equations)
    {
      equation = dynamic_cast<BinaryOpNode *>(equation->decreaseLeadsLagsPredeterminedVariables());
      assert(equation);
    }

  // No need to handle static_only_equations, since there are no leads/lags there
}

void
DynamicModel::substituteLogTransform()
{
  for (int symb_id : symbol_table.getVariablesWithLogTransform())
    {
      expr_t aux_def = AddLog(AddVariable(symb_id));
      int aux_symb_id = symbol_table.addLogTransformAuxiliaryVar(symb_id, 0, aux_def);

      for (auto &[id, definition] : local_variables_table)
        definition = definition->substituteLogTransform(symb_id, aux_symb_id);

      for (auto &equation : equations)
        equation = dynamic_cast<BinaryOpNode *>(equation->substituteLogTransform(symb_id, aux_symb_id));

      for (auto &equation : static_only_equations)
        equation = dynamic_cast<BinaryOpNode *>(equation->substituteLogTransform(symb_id, aux_symb_id));

      /*
        We add the following new equations:
        + X=exp(log_X) to the model
        + log_X=log(X) to the list of auxiliary equations

        In this way:
        + statements like X=1 in initval/endval blocks will be correctly
          handled (i.e. log_X will be initialized to 0 in this case), through
          the set_auxiliary_variables.m and dynamic_set_auxiliary_series.m files
        + computation of X in perfect foresight simulations will be done by
          simple evaluation when using block decomposition (X will belong to an
          block of type “evaluate”, or maybe even the epilogue)
      */
      addAuxEquation(AddEqual(AddVariable(aux_symb_id), aux_def));
      addEquation(AddEqual(AddVariable(symb_id), AddExp(AddVariable(aux_symb_id))),
                  nullopt, {});
    }
}

void
DynamicModel::checkNoWithLogTransform(const set<int> &eqnumbers)
{
  set<int> endos;
  for (int eq : eqnumbers)
    equations[eq]->collectVariables(SymbolType::endogenous, endos);

  const set<int> &with_log_transform = symbol_table.getVariablesWithLogTransform();

  vector<int> intersect;
  set_intersection(endos.begin(), endos.end(),
                   with_log_transform.begin(), with_log_transform.end(),
                   back_inserter(intersect));
  if (!intersect.empty())
    {
      cerr << "ERROR: the following variables are declared with var(log) and therefore cannot appear in a VAR/TCM/PAC equation: ";
      for (int symb_id : intersect)
        cerr << symbol_table.getName(symb_id) << " ";
      cerr << endl;
      exit(EXIT_FAILURE);
    }
}

void
DynamicModel::detrendEquations()
{
  // We go backwards in the list of trend_vars, to deal correctly with I(2) processes
  for (auto it = nonstationary_symbols_map.crbegin();
       it != nonstationary_symbols_map.crend(); ++it)
    {
      for (auto &equation : equations)
        {
          equation = dynamic_cast<BinaryOpNode *>(equation->detrend(it->first, it->second.first, it->second.second));
          assert(equation);
        }
      for (auto &equation : static_only_equations)
        {
          equation = dynamic_cast<BinaryOpNode *>(equation->detrend(it->first, it->second.first, it->second.second));
          assert(equation);
        }
    }

  for (auto &equation : equations)
    {
      equation = dynamic_cast<BinaryOpNode *>(equation->removeTrendLeadLag(trend_symbols_map));
      assert(equation);
    }
  for (auto &equation : static_only_equations)
    {
      equation = dynamic_cast<BinaryOpNode *>(equation->removeTrendLeadLag(trend_symbols_map));
      assert(equation);
    }
}

void
DynamicModel::removeTrendVariableFromEquations()
{
  for (auto &equation : equations)
    {
      equation = dynamic_cast<BinaryOpNode *>(equation->replaceTrendVar());
      assert(equation);
    }
  for (auto &equation : static_only_equations)
    {
      equation = dynamic_cast<BinaryOpNode *>(equation->replaceTrendVar());
      assert(equation);
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
  for (auto &[symb_id, expression] : local_variables_table)
    {
      try
        {
          double val = expression->eval(eval_context);
          eval_context[symb_id] = val;
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

void
DynamicModel::addStaticOnlyEquation(expr_t eq, optional<int> lineno, const map<string, string> &eq_tags)
{
  auto beq = dynamic_cast<BinaryOpNode *>(eq);
  assert(beq && beq->op_code == BinaryOpcode::equal);

  static_only_equations_equation_tags.add(static_only_equations.size(), eq_tags);
  static_only_equations.push_back(beq);
  static_only_equations_lineno.push_back(move(lineno));
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

void
DynamicModel::addOccbinEquation(expr_t eq, optional<int> lineno, const map<string, string> &eq_tags, const vector<string> &regimes_bind, const vector<string> &regimes_relax)
{
  auto beq = dynamic_cast<BinaryOpNode *>(eq);
  assert(beq && beq->op_code == BinaryOpcode::equal);

  // Construct the term to be added to the corresponding equation
  expr_t basic_term = AddMinus(beq->arg1, beq->arg2);
  expr_t term = basic_term;
  for (auto &regime : regimes_bind)
    {
      int param_id = symbol_table.getID(ParsingDriver::buildOccbinBindParamName(regime));
      term = AddTimes(term, AddVariable(param_id));
    }
  for (auto &regime : regimes_relax)
    {
      int param_id = symbol_table.getID(ParsingDriver::buildOccbinBindParamName(regime));
      term = AddTimes(term, AddMinus(One, AddVariable(param_id)));
    }

  // Create or update the dynamic equation
  optional<int> eqn { equation_tags.getEqnByTag("name", eq_tags.at("name")) };
  if (eqn)
    {
      BinaryOpNode *orig_eq { equations[*eqn] };
      /* In the following, we could have kept only orig_eq->arg1, but the
         following adds a (somewhat bizarre) support for equation snippets
         without “bind” nor “relax” */
      equations[*eqn] = AddEqual(AddPlus(AddMinus(orig_eq->arg1, orig_eq->arg2), term), Zero);
      // It’s unclear how to update lineno and tags, so don’t do it
    }
  else
    {
      auto eq_tags_dynamic = eq_tags;
      eq_tags_dynamic["dynamic"] = "";
      addEquation(AddEqual(term, Zero), lineno, eq_tags_dynamic);
    }

  // Create or update the static equation (corresponding to the pure relax regime)
  if (regimes_bind.empty())
    {
      /* Similar remark as above. We could have entirely skipped this
         equation updating, since normally there is only one such clause,
         but the following adds a (somewhat bizarre) support for equation
         snippets without “bind” nor “relax” */
      optional<int> eqn_static { static_only_equations_equation_tags.getEqnByTag("name", eq_tags.at("name")) };
      if (eqn_static)
        {
          BinaryOpNode *orig_eq { static_only_equations[*eqn_static] };
          static_only_equations[*eqn_static] = AddEqual(AddPlus(AddMinus(orig_eq->arg1, orig_eq->arg2), basic_term), Zero);
          // It’s unclear how to update lineno and tags, so don’t do it
        }
      else
        {
          auto eq_tags_static = eq_tags;
          eq_tags_static["static"] = "";
          addStaticOnlyEquation(AddEqual(basic_term, Zero), lineno, eq_tags_static);
        }
    }
}

bool
DynamicModel::isChecksumMatching(const string &basename) const
{
  stringstream buffer;

  // Write equation tags
  equation_tags.writeCheckSumInfo(buffer);

  constexpr ExprNodeOutputType buffer_type{ExprNodeOutputType::CDynamicModel};

  deriv_node_temp_terms_t tef_terms;
  temporary_terms_t temp_term_union;

  writeTemporaryTerms<buffer_type>(temporary_terms_derivatives[0],
                                   temp_term_union, temporary_terms_idxs,
                                   buffer, tef_terms);

  writeModelEquations<buffer_type>(buffer, temp_term_union);

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
      cerr << "ERROR: Can't open file " << filename.string() << endl;
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
  output << R"(, "dynamic_tmp_nbr": [)";
  for (bool printed_something {false};
       const auto &tts : temporary_terms_derivatives)
    {
      if (exchange(printed_something, true))
        output << ", ";
      output << tts.size();
    }
  output << "], ";
  writeJsonSparseIndicesHelper<true>(output);
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

      output << R"({ "number":)" << eq;
      if (equations_lineno[eq])
        output << R"(, "line":)" << *equations_lineno[eq];

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
  for (bool printed_something{false};
       const auto &[var, eqs] : variableMapping)
    {
      if (exchange(printed_something, true))
        output << ", ";
      output << R"({"name": ")" << symbol_table.getName(var) << R"(", "equations":[)";
      for (bool printed_something2{false};
           int it2 : eqs)
        if (auto tmp = equation_tags.getTagValueByEqnAndKey(it2, "name");
            tmp)
          {
            if (exchange(printed_something2, true))
              output << ", ";
            output << '"' << *tmp << '"';
          }
      output << "]}" << endl;
    }
  output << "]";
}

void
DynamicModel::writeJsonXrefsHelper(ostream &output, const map<pair<int, int>, set<int>> &xrefmap) const
{
  for (bool printed_something{false};
       const auto &[symb_lag, eqs] : xrefmap)
    {
      if (exchange(printed_something, true))
        output << ", ";
      output << R"({"name": ")" << symbol_table.getName(symb_lag.first) << R"(")"
             << R"(, "shift": )" << symb_lag.second
             << R"(, "equations": [)";
      for (bool printed_something2{false};
           int eq : eqs)
        {
          if (exchange(printed_something2, true))
            output << ", ";
          output << eq + 1;
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
              output << " " << getJacobianCol(varID, false) + 1;
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
         << R"("nstatic": )" << nstatic << ", " << endl
         << R"("nfwrd": )" << nfwrd << ", " << endl
         << R"("npred": )" << npred << ", " << endl
         << R"("nboth": )" << nboth << ", " << endl
         << R"("nsfwrd": )" << nfwrd+nboth << ", " << endl
         << R"("nspred": )" << npred+nboth << ", " << endl
         << R"("ndynamic": )" << npred+nboth+nfwrd << ", " << endl
         << R"("maximum_endo_lag": )" << max_endo_lag << ", " << endl
         << R"("maximum_endo_lead": )" << max_endo_lead << ", " << endl
         << R"("maximum_exo_lag": )" << max_exo_lag << ", " << endl
         << R"("maximum_exo_lead": )" << max_exo_lead << ", " << endl
         << R"("maximum_exo_det_lag": )" << max_exo_det_lag << ", " << endl
         << R"("maximum_exo_det_lead": )" << max_exo_det_lead << ", " << endl
         << R"("maximum_lag": )" << max_lag << ", " << endl
         << R"("maximum_lead": )" << max_lead << ", " << endl
         << R"("orig_maximum_endo_lag": )" << max_endo_lag_orig << "," << endl
         << R"("orig_maximum_endo_lead": )" << max_endo_lead_orig << "," << endl
         << R"("orig_maximum_exo_lag": )" << max_exo_lag_orig << "," << endl
         << R"("orig_maximum_exo_lead": )" << max_exo_lead_orig << "," << endl
         << R"("orig_maximum_exo_det_lag": )" << max_exo_det_lag_orig << "," << endl
         << R"("orig_maximum_exo_det_lead": )" << max_exo_det_lead_orig << "," << endl
         << R"("orig_maximum_lag": )" << max_lag_orig << "," << endl
         << R"("orig_maximum_lead": )" << max_lead_orig << "," << endl
         << R"("orig_maximum_lag_with_diffs_expanded": )" << max_lag_with_diffs_expanded_orig
	 << "," <<endl
	 << R"("NNZDerivatives": [)"; 
  for (int i = 1; i < static_cast<int>(NNZDerivatives.size()); i++)
    {
      output << (i > computed_derivs_order ? -1 : NNZDerivatives[i]);
      if (i < static_cast<int>(NNZDerivatives.size()) - 1)
	output << ", ";
    }
  output << "]}" 
	 << endl;
}

void
DynamicModel::writeJsonComputingPassOutput(ostream &output, bool writeDetails) const
{
  auto [mlv_output, d_output] { writeJsonComputingPassOutputHelper<true>(writeDetails) };

  if (writeDetails)
    output << R"("dynamic_model": {)";
  else
    output << R"("dynamic_model_simple": {)";
  output << mlv_output.str();
  for (const auto &it : d_output)
    output << ", " << it.str();
  output << "}";
}

void
DynamicModel::writeJsonParamsDerivatives(ostream &output, bool writeDetails) const
{
  if (!params_derivatives.size())
    return;

  auto [mlv_output, tt_output, rp_output, gp_output, rpp_output, gpp_output, hp_output, g3p_output]
    { writeJsonParamsDerivativesHelper<true>(writeDetails) };

  if (writeDetails)
    output << R"("dynamic_model_params_derivative": {)";
  else
    output << R"("dynamic_model_params_derivatives_simple": {)";
  output << mlv_output.str()
         << ", " << tt_output.str()
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
        cerr << "ERROR: in equation " << equation_tags.getTagValueByEqnAndKey(eq, "name").value_or(to_string(eq+1))
             << ", the pac_expectation operator references an unknown pac_model" << endl;
        exit(EXIT_FAILURE);
      }
}

void
DynamicModel::simplifyEquations()
{
  size_t last_subst_table_size = 0;
  map<VariableNode *, NumConstNode *> subst_table;
  // Equations with “mcp” tag are excluded, see dynare#1697
  findConstantEquationsWithoutMcpTag(subst_table);
  while (subst_table.size() != last_subst_table_size)
    {
      last_subst_table_size = subst_table.size();
      for (auto &[id, definition] : local_variables_table)
        definition = definition->replaceVarsInEquation(subst_table);
      for (auto &equation : equations)
        equation = dynamic_cast<BinaryOpNode *>(equation->replaceVarsInEquation(subst_table));
      for (auto &equation : static_only_equations)
        equation = dynamic_cast<BinaryOpNode *>(equation->replaceVarsInEquation(subst_table));
      subst_table.clear();
      findConstantEquationsWithoutMcpTag(subst_table);
    }
}

void
DynamicModel::checkNoRemainingPacTargetNonstationary() const
{
  for (size_t eq = 0; eq < equations.size(); eq++)
    if (equations[eq]->containsPacTargetNonstationary())
      {
        cerr << "ERROR: in equation " << equation_tags.getTagValueByEqnAndKey(eq, "name").value_or(to_string(eq+1))
             << ", the pac_target_nonstationary operator does not match a corresponding 'pac_target_info' block" << endl;
        exit(EXIT_FAILURE);
      }
}

void
DynamicModel::checkIsLinear() const
{
  if (!nonzero_hessian_eqs.empty())
    {
      cerr << "ERROR: If the model is declared linear the second derivatives must be equal to zero." << endl
           << "       The following equations have non-zero second derivatives:" << endl;
      for (auto it : nonzero_hessian_eqs)
        {
          cerr << "       * Eq # " << it+1;
          if (optional<string> eqname { equation_tags.getTagValueByEqnAndKey(it, "name") };
              eqname)
            cerr << " [" << *eqname << "]";
          cerr << endl;
        }
      exit(EXIT_FAILURE);
    }
}
