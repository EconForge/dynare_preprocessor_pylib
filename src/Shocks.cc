/*
 * Copyright (C) 2003-2017 Dynare Team
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

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <utility>

#include "Shocks.hh"

AbstractShocksStatement::AbstractShocksStatement(bool mshocks_arg,
                                                 bool overwrite_arg,
                                                 det_shocks_t det_shocks_arg,
                                                 const SymbolTable &symbol_table_arg) :
  mshocks(mshocks_arg),
  overwrite(overwrite_arg),
  det_shocks(move(det_shocks_arg)),
  symbol_table(symbol_table_arg)
{
}

void
AbstractShocksStatement::writeDetShocks(ostream &output) const
{
  int exo_det_length = 0;

  for (const auto & det_shock : det_shocks)
    {
      int id = symbol_table.getTypeSpecificID(det_shock.first) + 1;
      bool exo_det = (symbol_table.getType(det_shock.first) == eExogenousDet);

      for (size_t i = 0; i < det_shock.second.size(); i++)
        {
          const int &period1 = det_shock.second[i].period1;
          const int &period2 = det_shock.second[i].period2;
          const expr_t value = det_shock.second[i].value;

          output << "M_.det_shocks = [ M_.det_shocks;" << endl
                 << "struct('exo_det'," << (int) exo_det
                 << ",'exo_id'," << id
                 << ",'multiplicative'," << (int) mshocks
                 << ",'periods'," << period1 << ":" << period2
                 << ",'value',";
          value->writeOutput(output);
          output << ") ];" << endl;

          if (exo_det && (period2 > exo_det_length))
            exo_det_length = period2;
        }
    }
  output << "M_.exo_det_length = " << exo_det_length << ";\n";
}

void
AbstractShocksStatement::writeJsonDetShocks(ostream &output) const
{
  output << "\"deterministic_shocks\": [";
  for (det_shocks_t::const_iterator it = det_shocks.begin();
       it != det_shocks.end(); it++)
    {
      if (it != det_shocks.begin())
        output << ", ";
      output << "{\"var\": \"" << symbol_table.getName(it->first) << "\", "
             << "\"values\": [";
      for (vector<DetShockElement>::const_iterator it1 = it->second.begin();
           it1 != it->second.end(); it1++)
        {
          if (it1 != it->second.begin())
            output << ", ";
          output << "{\"period1\": " << it1->period1 << ", "
                 << "\"period2\": " << it1->period2 << ", "
                 << "\"value\": \"";
          it1->value->writeJsonOutput(output, {}, {});
          output << "\"}";
        }
      output << "]}";
    }
  output << "]";
}

ShocksStatement::ShocksStatement(bool overwrite_arg,
                                 const det_shocks_t &det_shocks_arg,
                                 var_and_std_shocks_t var_shocks_arg,
                                 var_and_std_shocks_t std_shocks_arg,
                                 covar_and_corr_shocks_t covar_shocks_arg,
                                 covar_and_corr_shocks_t corr_shocks_arg,
                                 const SymbolTable &symbol_table_arg) :
  AbstractShocksStatement(false, overwrite_arg, det_shocks_arg, symbol_table_arg),
  var_shocks(move(var_shocks_arg)),
  std_shocks(move(std_shocks_arg)),
  covar_shocks(move(covar_shocks_arg)),
  corr_shocks(move(corr_shocks_arg))
{
}

void
ShocksStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "%" << endl
         << "% SHOCKS instructions" << endl
         << "%" << endl;

  if (overwrite)
    {
      output << "M_.det_shocks = [];" << endl;

      output << "M_.Sigma_e = zeros(" << symbol_table.exo_nbr() << ", "
             << symbol_table.exo_nbr() << ");" << endl
             << "M_.Correlation_matrix = eye(" << symbol_table.exo_nbr() << ", "
             << symbol_table.exo_nbr() << ");" << endl;

      if (has_calibrated_measurement_errors())
        output << "M_.H = zeros(" << symbol_table.observedVariablesNbr() << ", "
               << symbol_table.observedVariablesNbr() << ");" << endl
               << "M_.Correlation_matrix_ME = eye(" << symbol_table.observedVariablesNbr() << ", "
               << symbol_table.observedVariablesNbr() << ");" << endl;
      else
        output << "M_.H = 0;" << endl
               << "M_.Correlation_matrix_ME = 1;" << endl;

    }

  writeDetShocks(output);
  writeVarAndStdShocks(output);
  writeCovarAndCorrShocks(output);

  /* M_.sigma_e_is_diagonal is initialized to 1 by ModFile.cc.
     If there are no off-diagonal elements, and we are not in overwrite mode,
     then we don't reset it to 1, since there might be previous shocks blocks
     with off-diagonal elements. */
  if (covar_shocks.size()+corr_shocks.size() > 0)
    output << "M_.sigma_e_is_diagonal = 0;" << endl;
  else if (overwrite)
    output << "M_.sigma_e_is_diagonal = 1;" << endl;
}

void
ShocksStatement::writeJsonOutput(ostream &output) const
{
  output << "{\"statementName\": \"shocks\""
         << ", \"overwrite\": ";
  if (overwrite)
    output << "true";
  else
    output << "false";
  if (!det_shocks.empty())
    {
      output << ", ";
      writeJsonDetShocks(output);
    }
  output<< ", \"variance\": [";
  for (var_and_std_shocks_t::const_iterator it = var_shocks.begin(); it != var_shocks.end(); it++)
    {
      if (it != var_shocks.begin())
        output << ", ";
      output << "{\"name\": \"" << symbol_table.getName(it->first) << "\", "
             << "\"variance\": \"";
      it->second->writeJsonOutput(output, {}, {});
      output << "\"}";
    }
  output << "]"
         << ", \"stderr\": [";
  for (var_and_std_shocks_t::const_iterator it = std_shocks.begin(); it != std_shocks.end(); it++)
    {
      if (it != std_shocks.begin())
        output << ", ";
      output << "{\"name\": \"" << symbol_table.getName(it->first) << "\", "
             << "\"stderr\": \"";
      it->second->writeJsonOutput(output, {}, {});
      output << "\"}";
    }
  output << "]"
         << ", \"covariance\": [";
  for (covar_and_corr_shocks_t::const_iterator it = covar_shocks.begin(); it != covar_shocks.end(); it++)
    {
      if (it != covar_shocks.begin())
        output << ", ";
      output << "{"
             << "\"name\": \"" << symbol_table.getName(it->first.first) << "\", "
             << "\"name2\": \"" << symbol_table.getName(it->first.second) << "\", "
             << "\"covariance\": \"";
      it->second->writeJsonOutput(output, {}, {});
      output << "\"}";
    }
  output << "]"
         << ", \"correlation\": [";
  for (covar_and_corr_shocks_t::const_iterator it = corr_shocks.begin(); it != corr_shocks.end(); it++)
    {
      if (it != corr_shocks.begin())
        output << ", ";
      output << "{"
             << "\"name\": \"" << symbol_table.getName(it->first.first) << "\", "
             << "\"name2\": \"" << symbol_table.getName(it->first.second) << "\", "
             << "\"correlation\": \"";
      it->second->writeJsonOutput(output, {}, {});
      output << "\"}";
    }
  output << "]"
         << "}";
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
  string matrix, corr_matrix;
  int id1, id2;
  if (type1 == eExogenous)
    {
      matrix = "M_.Sigma_e";
      corr_matrix = "M_.Correlation_matrix";
      id1 = symbol_table.getTypeSpecificID(it->first.first) + 1;
      id2 = symbol_table.getTypeSpecificID(it->first.second) + 1;
    }
  else
    {
      matrix = "M_.H";
      corr_matrix = "M_.Correlation_matrix_ME";
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

  if (corr)
    {
      output << corr_matrix << "(" << id1 << ", " << id2 << ") = ";
      it->second->writeOutput(output);
      output << ";" << endl
             << corr_matrix << "(" << id2 << ", " << id1 << ") = "
             << corr_matrix << "(" << id1 << ", " << id2 << ");" << endl;
    }
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
  /* Error out if variables are not of the right type. This must be done here
     and not at parsing time (see #448).
     Also Determine if there is a calibrated measurement error */
  for (auto var_shock : var_shocks)
    {
      if (symbol_table.getType(var_shock.first) != eExogenous
          && !symbol_table.isObservedVariable(var_shock.first))
        {
          cerr << "shocks: setting a variance on '"
               << symbol_table.getName(var_shock.first) << "' is not allowed, because it is neither an exogenous variable nor an observed endogenous variable" << endl;
          exit(EXIT_FAILURE);
        }
    }

  for (auto std_shock : std_shocks)
    {
      if (symbol_table.getType(std_shock.first) != eExogenous
          && !symbol_table.isObservedVariable(std_shock.first))
        {
          cerr << "shocks: setting a standard error on '"
               << symbol_table.getName(std_shock.first) << "' is not allowed, because it is neither an exogenous variable nor an observed endogenous variable" << endl;
          exit(EXIT_FAILURE);
        }
    }

  for (const auto & covar_shock : covar_shocks)
    {
      int symb_id1 = covar_shock.first.first;
      int symb_id2 = covar_shock.first.second;

      if (!((symbol_table.getType(symb_id1) == eExogenous
             && symbol_table.getType(symb_id2) == eExogenous)
            || (symbol_table.isObservedVariable(symb_id1)
                && symbol_table.isObservedVariable(symb_id2))))
        {
          cerr << "shocks: setting a covariance between '"
               << symbol_table.getName(symb_id1) << "' and '"
               << symbol_table.getName(symb_id2) << "'is not allowed; covariances can only be specified for exogenous or observed endogenous variables of same type" << endl;
          exit(EXIT_FAILURE);
        }
    }

  for (const auto & corr_shock : corr_shocks)
    {
      int symb_id1 = corr_shock.first.first;
      int symb_id2 = corr_shock.first.second;

      if (!((symbol_table.getType(symb_id1) == eExogenous
             && symbol_table.getType(symb_id2) == eExogenous)
            || (symbol_table.isObservedVariable(symb_id1)
                && symbol_table.isObservedVariable(symb_id2))))
        {
          cerr << "shocks: setting a correlation between '"
               << symbol_table.getName(symb_id1) << "' and '"
               << symbol_table.getName(symb_id2) << "'is not allowed; correlations can only be specified for exogenous or observed endogenous variables of same type" << endl;
          exit(EXIT_FAILURE);
        }
    }

  // Determine if there is a calibrated measurement error
  mod_file_struct.calibrated_measurement_errors |= has_calibrated_measurement_errors();

  // Fill in mod_file_struct.parameters_with_shocks_values (related to #469)
  for (auto var_shock : var_shocks)
    var_shock.second->collectVariables(eParameter, mod_file_struct.parameters_within_shocks_values);
  for (auto std_shock : std_shocks)
    std_shock.second->collectVariables(eParameter, mod_file_struct.parameters_within_shocks_values);
  for (const auto & covar_shock : covar_shocks)
    covar_shock.second->collectVariables(eParameter, mod_file_struct.parameters_within_shocks_values);
  for (const auto & corr_shock : corr_shocks)
    corr_shock.second->collectVariables(eParameter, mod_file_struct.parameters_within_shocks_values);

}

bool
ShocksStatement::has_calibrated_measurement_errors() const
{
  for (auto var_shock : var_shocks)
    if (symbol_table.isObservedVariable(var_shock.first))
      return true;

  for (auto std_shock : std_shocks)
    if (symbol_table.isObservedVariable(std_shock.first))
      return true;

  for (const auto & covar_shock : covar_shocks)
    if (symbol_table.isObservedVariable(covar_shock.first.first)
        || symbol_table.isObservedVariable(covar_shock.first.second))
      return true;

  for (const auto & corr_shock : corr_shocks)
    if (symbol_table.isObservedVariable(corr_shock.first.first)
        || symbol_table.isObservedVariable(corr_shock.first.second))
      return true;

  return false;
}

MShocksStatement::MShocksStatement(bool overwrite_arg,
                                   const det_shocks_t &det_shocks_arg,
                                   const SymbolTable &symbol_table_arg) :
  AbstractShocksStatement(true, overwrite_arg, det_shocks_arg, symbol_table_arg)
{
}

void
MShocksStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "%" << endl
         << "% MSHOCKS instructions" << endl
         << "%" << endl;

  if (overwrite)
    output << "M_.det_shocks = [];" << endl;

  writeDetShocks(output);
}

ConditionalForecastPathsStatement::ConditionalForecastPathsStatement(AbstractShocksStatement::det_shocks_t paths_arg,
                                                                     const SymbolTable &symbol_table_arg) :
  paths(move(paths_arg)),
  symbol_table(symbol_table_arg),
  path_length(-1)
{
}

void
ConditionalForecastPathsStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  for (const auto & path : paths)
    {
      int this_path_length = 0;
      const vector<AbstractShocksStatement::DetShockElement> &elems = path.second;
      for (auto elem : elems)
        // Period1 < Period2, as enforced in ParsingDriver::add_period()
        this_path_length = max(this_path_length, elem.period2);
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
ConditionalForecastPathsStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  assert(path_length > 0);
  output << "constrained_vars_ = [];" << endl
         << "constrained_paths_ = zeros(" << paths.size() << ", " << path_length << ");" << endl;

  int k = 1;
  for (AbstractShocksStatement::det_shocks_t::const_iterator it = paths.begin();
       it != paths.end(); it++, k++)
    {
      if (it == paths.begin())
        output << "constrained_vars_ = " << symbol_table.getTypeSpecificID(it->first) + 1 << ";" << endl;
      else
        output << "constrained_vars_ = [constrained_vars_; " << symbol_table.getTypeSpecificID(it->first) + 1 << "];" << endl;

      const vector<AbstractShocksStatement::DetShockElement> &elems = it->second;
      for (auto elem : elems)
        for (int j = elem.period1; j <= elem.period2; j++)
          {
            output << "constrained_paths_(" << k << "," << j << ")=";
            elem.value->writeOutput(output);
            output << ";" << endl;
          }
    }
}

MomentCalibration::MomentCalibration(constraints_t constraints_arg,
                                     const SymbolTable &symbol_table_arg)
  : constraints(move(constraints_arg)), symbol_table(symbol_table_arg)
{
}

void
MomentCalibration::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "options_.endogenous_prior_restrictions.moment = {" << endl;
  for (const auto & c : constraints)
    {
      output << "'" << symbol_table.getName(c.endo1) << "', "
             << "'" << symbol_table.getName(c.endo2) << "', "
             << c.lags << ", "
             << "[ " << c.lower_bound << ", " << c.upper_bound << " ];"
             << endl;
    }
  output << "};" << endl;
}

void
MomentCalibration::writeJsonOutput(ostream &output) const
{
  output << "{\"statementName\": \"moment_calibration\""
         << ", \"moment_calibration_criteria\": [";
  for (constraints_t::const_iterator it = constraints.begin(); it != constraints.end(); it++)
    {
      if (it != constraints.begin())
        output << ", ";
      output << "{\"endogenous1\": \"" << symbol_table.getName(it->endo1) << "\""
             << ", \"endogenous2\": \"" << symbol_table.getName(it->endo2) << "\""
             << ", \"lags\": \"" << it->lags << "\""
             << ", \"lower_bound\": \"" << it->lower_bound << "\""
             << ", \"upper_bound\": \"" << it->upper_bound << "\""
             << "}";
    }
  output << "]"
         << "}";
}

IrfCalibration::IrfCalibration(constraints_t constraints_arg,
                               const SymbolTable &symbol_table_arg,
                               OptionsList options_list_arg)
  : constraints(move(constraints_arg)), symbol_table(symbol_table_arg), options_list(move(options_list_arg))
{
}

void
IrfCalibration::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  options_list.writeOutput(output);

  output << "options_.endogenous_prior_restrictions.irf = {" << endl;
  for (const auto & c : constraints)
    {
      output << "'" << symbol_table.getName(c.endo) << "', "
             << "'" << symbol_table.getName(c.exo) << "', "
             << c.periods << ", "
             << "[ " << c.lower_bound << ", " << c.upper_bound << " ];"
             << endl;
    }
  output << "};" << endl;
}

void
IrfCalibration::writeJsonOutput(ostream &output) const
{
  output << "{\"statementName\": \"irf_calibration\"";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }

  output << ", \"irf_restrictions\": [";
  for (constraints_t::const_iterator it = constraints.begin(); it != constraints.end(); it++)
    {
      if (it != constraints.begin())
        output << ", ";
      output << "{\"endogenous\": \"" << symbol_table.getName(it->endo) << "\""
             << ", \"exogenous\": \"" << symbol_table.getName(it->exo) << "\""
             << ", \"periods\": \"" << it->periods << "\""
             << ", \"lower_bound\": \"" << it->lower_bound << "\""
             << ", \"upper_bound\": \"" << it->upper_bound << "\""
             << "}";
    }
  output << "]"
         << "}";
}

ShockGroupsStatement::ShockGroupsStatement(group_t shock_groups_arg, string name_arg)
  : shock_groups(move(shock_groups_arg)), name(move(name_arg))
{
}

void
ShockGroupsStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  int i = 1;
  bool unique_label = true;
  for (vector<Group>::const_iterator it = shock_groups.begin(); it != shock_groups.end(); it++, unique_label = true)
    {
      for (vector<Group>::const_iterator it1 = it+1; it1 != shock_groups.end(); it1++)
        if (it->name == it1->name)
          {
            unique_label = false;
            cerr << "Warning: shock group label '" << it->name << "' has been reused. "
                 << "Only using the last definition." << endl;
            break;
          }

      if (unique_label)
        {
          output << "M_.shock_groups." << name
                 << ".group" << i << ".label = '" << it->name << "';" << endl
                 << "M_.shock_groups." << name
                 << ".group" << i << ".shocks = {";
          for (const auto & it1 : it->list)
            output << " '" << it1 << "'";
          output << "};" << endl;
          i++;
        }
    }
}
