/*
 * Copyright (C) 2003-2012 Dynare Team
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

using namespace std;

#include <cassert>
#include <cstdlib>
#include <iostream>

#include "Shocks.hh"

AbstractShocksStatement::AbstractShocksStatement(bool mshocks_arg,
                                                 const det_shocks_t &det_shocks_arg,
                                                 const SymbolTable &symbol_table_arg) :
  mshocks(mshocks_arg),
  det_shocks(det_shocks_arg),
  symbol_table(symbol_table_arg)
{
}

void
AbstractShocksStatement::writeDetShocks(ostream &output) const
{
  int exo_det_length = 0;

  for (det_shocks_t::const_iterator it = det_shocks.begin();
       it != det_shocks.end(); it++)
    {
      int id = symbol_table.getTypeSpecificID(it->first) + 1;
      bool exo_det = (symbol_table.getType(it->first) == eExogenousDet);
      int set_shocks_index = ((int) mshocks) + 2 * ((int) exo_det);

      for (size_t i = 0; i < it->second.size(); i++)
        {
          const int &period1 = it->second[i].period1;
          const int &period2 = it->second[i].period2;
          const expr_t value = it->second[i].value;

          if (period1 == period2)
            {
              output << "set_shocks(" << set_shocks_index << "," << period1
                     << ", " << id << ", ";
              value->writeOutput(output);
              output << ");" << endl;
            }
          else
            {
              output << "set_shocks(" << set_shocks_index << "," << period1
                     << ":" << period2 << ", " << id << ", ";
              value->writeOutput(output);
              output << ");" << endl;
            }

          if (exo_det && (period2 > exo_det_length))
            exo_det_length = period2;
        }
    }
  output << "M_.exo_det_length = " << exo_det_length << ";\n";
}

ShocksStatement::ShocksStatement(const det_shocks_t &det_shocks_arg,
                                 const var_and_std_shocks_t &var_shocks_arg,
                                 const var_and_std_shocks_t &std_shocks_arg,
                                 const covar_and_corr_shocks_t &covar_shocks_arg,
                                 const covar_and_corr_shocks_t &corr_shocks_arg,
                                 const SymbolTable &symbol_table_arg) :
  AbstractShocksStatement(false, det_shocks_arg, symbol_table_arg),
  var_shocks(var_shocks_arg),
  std_shocks(std_shocks_arg),
  covar_shocks(covar_shocks_arg),
  corr_shocks(corr_shocks_arg)
{
}

void
ShocksStatement::writeOutput(ostream &output, const string &basename) const
{
  output << "%" << endl
         << "% SHOCKS instructions" << endl
         << "%" << endl;

  // Write instruction that initializes a shock
  output << "make_ex_;" << endl;

  writeDetShocks(output);
  writeVarAndStdShocks(output);
  writeCovarAndCorrShocks(output);
  if (covar_shocks.size()+corr_shocks.size() > 0)
    output << "M_.sigma_e_is_diagonal = 0;" << endl;
  else
    output << "M_.sigma_e_is_diagonal = 1;" << endl;
}

void
ShocksStatement::writeVarOrStdShock(ostream &output, var_and_std_shocks_t::const_iterator &it,
                                    bool stddev) const
{
  SymbolType type = symbol_table.getType(it->first);
  assert(type == eExogenous || symbol_table.isObservedVariable(it->first));

  int id;
  if (type == eExogenous)
    {
      output << "M_.Sigma_e(";
      id = symbol_table.getTypeSpecificID(it->first) + 1;
    }
  else
    {
      output << "M_.H(";
      id = symbol_table.getObservedVariableIndex(it->first) + 1;
    }

  output << id << ", " << id << ") = ";
  if (stddev)
    output << "(";
  it->second->writeOutput(output);
  if (stddev)
    output << ")^2";
  output << ";" << endl;
}

void
ShocksStatement::writeVarAndStdShocks(ostream &output) const
{
  var_and_std_shocks_t::const_iterator it;

  for (it = var_shocks.begin(); it != var_shocks.end(); it++)
    writeVarOrStdShock(output, it, false);

  for (it = std_shocks.begin(); it != std_shocks.end(); it++)
    writeVarOrStdShock(output, it, true);
}

void
ShocksStatement::writeCovarOrCorrShock(ostream &output, covar_and_corr_shocks_t::const_iterator &it,
                                       bool corr) const
{
  SymbolType type1 = symbol_table.getType(it->first.first);
  SymbolType type2 = symbol_table.getType(it->first.second);
  assert((type1 == eExogenous && type2 == eExogenous)
         || (symbol_table.isObservedVariable(it->first.first) && symbol_table.isObservedVariable(it->first.second)));
  string matrix;
  int id1, id2;
  if (type1 == eExogenous)
    {
      matrix = "M_.Sigma_e";
      id1 = symbol_table.getTypeSpecificID(it->first.first) + 1;
      id2 = symbol_table.getTypeSpecificID(it->first.second) + 1;
    }
  else
    {
      matrix = "M_.H";
      id1 = symbol_table.getObservedVariableIndex(it->first.first) + 1;
      id2 = symbol_table.getObservedVariableIndex(it->first.second) + 1;
    }

  output << matrix << "(" << id1 << ", " << id2 << ") = ";
  it->second->writeOutput(output);
  if (corr)
    output << "*sqrt(" << matrix << "(" << id1 << ", " << id1 << ")*"
           << matrix << "(" << id2 << ", " << id2 << "))";
  output << ";" << endl
         << matrix << "(" << id2 << ", " << id1 << ") = "
         << matrix << "(" << id1 << ", " << id2 << ");" << endl;
}

void
ShocksStatement::writeCovarAndCorrShocks(ostream &output) const
{
  covar_and_corr_shocks_t::const_iterator it;

  for (it = covar_shocks.begin(); it != covar_shocks.end(); it++)
    writeCovarOrCorrShock(output, it, false);

  for (it = corr_shocks.begin(); it != corr_shocks.end(); it++)
    writeCovarOrCorrShock(output, it, true);
}

void
ShocksStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  // Workaround for trac ticket #35
  mod_file_struct.shocks_present_but_simul_not_yet = true;

  // Determine if there is a calibrated measurement error
  for (var_and_std_shocks_t::const_iterator it = var_shocks.begin();
       it != var_shocks.end(); it++)
    if (symbol_table.isObservedVariable(it->first))
      mod_file_struct.calibrated_measurement_errors = true;

  for (var_and_std_shocks_t::const_iterator it = std_shocks.begin();
       it != std_shocks.end(); it++)
    if (symbol_table.isObservedVariable(it->first))
      mod_file_struct.calibrated_measurement_errors = true;

  for (covar_and_corr_shocks_t::const_iterator it = covar_shocks.begin();
       it != covar_shocks.end(); it++)
    if (symbol_table.isObservedVariable(it->first.first)
        || symbol_table.isObservedVariable(it->first.second))
      mod_file_struct.calibrated_measurement_errors = true;

  for (covar_and_corr_shocks_t::const_iterator it = corr_shocks.begin();
       it != corr_shocks.end(); it++)
    if (symbol_table.isObservedVariable(it->first.first)
        || symbol_table.isObservedVariable(it->first.second))
      mod_file_struct.calibrated_measurement_errors = true;
}

MShocksStatement::MShocksStatement(const det_shocks_t &det_shocks_arg,
                                   const SymbolTable &symbol_table_arg) :
  AbstractShocksStatement(true, det_shocks_arg, symbol_table_arg)
{
}

void
MShocksStatement::writeOutput(ostream &output, const string &basename) const
{
  output << "%" << endl
         << "% MSHOCKS instructions" << endl
         << "%" << endl;

  // Write instruction that initializes a shock
  output << "make_ex_;" << endl;

  writeDetShocks(output);
}

void
MShocksStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  // Workaround for trac ticket #35
  mod_file_struct.shocks_present_but_simul_not_yet = true;
}

ConditionalForecastPathsStatement::ConditionalForecastPathsStatement(const AbstractShocksStatement::det_shocks_t &paths_arg, const SymbolTable &symbol_table_arg) :
  paths(paths_arg),
  symbol_table(symbol_table_arg),
  path_length(-1)
{
}

void
ConditionalForecastPathsStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  for (AbstractShocksStatement::det_shocks_t::const_iterator it = paths.begin();
       it != paths.end(); it++)
    {
      int this_path_length = 0;
      const vector<AbstractShocksStatement::DetShockElement> &elems = it->second;
      for (int i = 0; i < (int) elems.size(); i++)
        // Period1 < Period2, as enforced in ParsingDriver::add_period()
        this_path_length = max(this_path_length, elems[i].period2);
      if (path_length == -1)
        path_length = this_path_length;
      else if (path_length != this_path_length)
        {
          cerr << "conditional_forecast_paths: all constrained paths must have the same length!" << endl;
          exit(EXIT_FAILURE);
        }
    }
}

void
ConditionalForecastPathsStatement::writeOutput(ostream &output, const string &basename) const
{
  assert(path_length > 0);
  output << "constrained_vars_ = [];" << endl
         << "constrained_paths_ = zeros(" << paths.size() << ", " << path_length << ");" << endl;

  int k = 1;

  for (AbstractShocksStatement::det_shocks_t::const_iterator it = paths.begin();
       it != paths.end(); it++)
    {
      if (it == paths.begin())
        output << "constrained_vars_ = " << it->first +1 << ";" << endl;
      else
        output << "constrained_vars_ = [constrained_vars_; " << it->first +1 << "];" << endl;

      const vector<AbstractShocksStatement::DetShockElement> &elems = it->second;
      for (int i = 0; i < (int) elems.size(); i++)
        for (int j = elems[i].period1; j <= elems[i].period2; j++)
          {
            output << "constrained_paths_(" << k << "," << j << ")=";
            elems[i].value->writeOutput(output);
            output << ";" << endl;
          }
      k++;
    }
}
