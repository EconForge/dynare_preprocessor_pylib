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

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <utility>

#include "Shocks.hh"

AbstractShocksStatement::AbstractShocksStatement(bool mshocks_arg,
                                                 bool overwrite_arg,
                                                 det_shocks_t det_shocks_arg,
                                                 const SymbolTable &symbol_table_arg) :
  mshocks{mshocks_arg},
  overwrite{overwrite_arg},
  det_shocks{move(det_shocks_arg)},
  symbol_table{symbol_table_arg}
{
}

void
AbstractShocksStatement::writeDetShocks(ostream &output) const
{
  int exo_det_length = 0;

  for (const auto & [id, shock_vec] : det_shocks)
    {
      bool exo_det = (symbol_table.getType(id) == SymbolType::exogenousDet);

      for (const auto &[period1, period2, value] : shock_vec)
        {
          output << "M_.det_shocks = [ M_.det_shocks;" << endl
                 << boolalpha
                 << "struct('exo_det'," << exo_det
                 << ",'exo_id'," << symbol_table.getTypeSpecificID(id)+1
                 << ",'multiplicative'," << mshocks
                 << ",'periods'," << period1 << ":" << period2
                 << ",'value',";
          value->writeOutput(output);
          output << ") ];" << endl;

          if (exo_det && period2 > exo_det_length)
            exo_det_length = period2;
        }
    }
  output << "M_.exo_det_length = " << exo_det_length << ";\n";
}

void
AbstractShocksStatement::writeJsonDetShocks(ostream &output) const
{
  output << R"("deterministic_shocks": [)";
  for (auto it = det_shocks.begin(); it != det_shocks.end(); ++it)
    {
      if (it != det_shocks.begin())
        output << ", ";
      output << R"({"var": ")" << symbol_table.getName(it->first) << R"(", )"
             << R"("values": [)";
      for (auto it1 = it->second.begin(); it1 != it->second.end(); ++it1)
        {
          auto [period1, period2, value] = *it1;
          if (it1 != it->second.begin())
            output << ", ";
          output << R"({"period1": )" << period1 << ", "
                 << R"("period2": )" << period2 << ", "
                 << R"("value": ")";
          value->writeJsonOutput(output, {}, {});
          output << R"("})";
        }
      output << "]}";
    }
  output << "]";
}

ShocksStatement::ShocksStatement(bool overwrite_arg,
                                 det_shocks_t det_shocks_arg,
                                 var_and_std_shocks_t var_shocks_arg,
                                 var_and_std_shocks_t std_shocks_arg,
                                 covar_and_corr_shocks_t covar_shocks_arg,
                                 covar_and_corr_shocks_t corr_shocks_arg,
                                 const SymbolTable &symbol_table_arg) :
  AbstractShocksStatement{false, overwrite_arg, move(det_shocks_arg), symbol_table_arg},
  var_shocks{move(var_shocks_arg)},
  std_shocks{move(std_shocks_arg)},
  covar_shocks{move(covar_shocks_arg)},
  corr_shocks{move(corr_shocks_arg)}
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
  output << R"({"statementName": "shocks")"
         << R"(, "overwrite": )" << boolalpha << overwrite;
  if (!det_shocks.empty())
    {
      output << ", ";
      writeJsonDetShocks(output);
    }
  output<< R"(, "variance": [)";
  for (auto it = var_shocks.begin(); it != var_shocks.end(); ++it)
    {
      if (it != var_shocks.begin())
        output << ", ";
      output << R"({"name": ")" << symbol_table.getName(it->first) << R"(", )"
             << R"("variance": ")";
      it->second->writeJsonOutput(output, {}, {});
      output << R"("})";
    }
  output << "]"
         << R"(, "stderr": [)";
  for (auto it = std_shocks.begin(); it != std_shocks.end(); it++)
    {
      if (it != std_shocks.begin())
        output << ", ";
      output << R"({"name": ")" << symbol_table.getName(it->first) << R"(", )"
             << R"("stderr": ")";
      it->second->writeJsonOutput(output, {}, {});
      output << R"("})";
    }
  output << "]"
         << R"(, "covariance": [)";
  for (auto it = covar_shocks.begin(); it != covar_shocks.end(); ++it)
    {
      if (it != covar_shocks.begin())
        output << ", ";
      output << "{"
             << R"("name": ")" << symbol_table.getName(it->first.first) << R"(", )"
             << R"("name2": ")" << symbol_table.getName(it->first.second) << R"(", )"
             << R"("covariance": ")";
      it->second->writeJsonOutput(output, {}, {});
      output << R"("})";
    }
  output << "]"
         << R"(, "correlation": [)";
  for (auto it = corr_shocks.begin(); it != corr_shocks.end(); ++it)
    {
      if (it != corr_shocks.begin())
        output << ", ";
      output << "{"
             << R"("name": ")" << symbol_table.getName(it->first.first) << R"(", )"
             << R"("name2": ")" << symbol_table.getName(it->first.second) << R"(", )"
             << R"("correlation": ")";
      it->second->writeJsonOutput(output, {}, {});
      output << R"("})";
    }
  output << "]"
         << "}";
}

void
ShocksStatement::writeVarOrStdShock(ostream &output, var_and_std_shocks_t::const_iterator &it,
                                    bool stddev) const
{
  SymbolType type = symbol_table.getType(it->first);
  assert(type == SymbolType::exogenous || symbol_table.isObservedVariable(it->first));

  int id;
  if (type == SymbolType::exogenous)
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
  for (auto it = var_shocks.begin(); it != var_shocks.end(); ++it)
    writeVarOrStdShock(output, it, false);

  for (auto it = std_shocks.begin(); it != std_shocks.end(); ++it)
    writeVarOrStdShock(output, it, true);
}

void
ShocksStatement::writeCovarOrCorrShock(ostream &output, covar_and_corr_shocks_t::const_iterator &it,
                                       bool corr) const
{
  SymbolType type1 = symbol_table.getType(it->first.first);
  SymbolType type2 = symbol_table.getType(it->first.second);
  assert((type1 == SymbolType::exogenous && type2 == SymbolType::exogenous)
         || (symbol_table.isObservedVariable(it->first.first) && symbol_table.isObservedVariable(it->first.second)));
  string matrix, corr_matrix;
  int id1, id2;
  if (type1 == SymbolType::exogenous)
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
  for (auto it = covar_shocks.begin(); it != covar_shocks.end(); ++it)
    writeCovarOrCorrShock(output, it, false);

  for (auto it = corr_shocks.begin(); it != corr_shocks.end(); ++it)
    writeCovarOrCorrShock(output, it, true);
}

void
ShocksStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  /* Error out if variables are not of the right type. This must be done here
     and not at parsing time (see #448).
     Also Determine if there is a calibrated measurement error */
  for (auto [id, val] : var_shocks)
    {
      if (symbol_table.getType(id) != SymbolType::exogenous
          && !symbol_table.isObservedVariable(id))
        {
          cerr << "shocks: setting a variance on '"
               << symbol_table.getName(id) << "' is not allowed, because it is neither an exogenous variable nor an observed endogenous variable" << endl;
          exit(EXIT_FAILURE);
        }
    }

  for (auto [id, val] : std_shocks)
    {
      if (symbol_table.getType(id) != SymbolType::exogenous
          && !symbol_table.isObservedVariable(id))
        {
          cerr << "shocks: setting a standard error on '"
               << symbol_table.getName(id) << "' is not allowed, because it is neither an exogenous variable nor an observed endogenous variable" << endl;
          exit(EXIT_FAILURE);
        }
    }

  for (const auto & [ids, val] : covar_shocks)
    {
      int symb_id1 = ids.first, symb_id2 = ids.second;

      if (!((symbol_table.getType(symb_id1) == SymbolType::exogenous
             && symbol_table.getType(symb_id2) == SymbolType::exogenous)
            || (symbol_table.isObservedVariable(symb_id1)
                && symbol_table.isObservedVariable(symb_id2))))
        {
          cerr << "shocks: setting a covariance between '"
               << symbol_table.getName(symb_id1) << "' and '"
               << symbol_table.getName(symb_id2) << "'is not allowed; covariances can only be specified for exogenous or observed endogenous variables of same type" << endl;
          exit(EXIT_FAILURE);
        }
    }

  for (const auto & [ids, val] : corr_shocks)
    {
      int symb_id1 = ids.first, symb_id2 = ids.second;

      if (!((symbol_table.getType(symb_id1) == SymbolType::exogenous
             && symbol_table.getType(symb_id2) == SymbolType::exogenous)
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
  for (auto [id, val] : var_shocks)
    val->collectVariables(SymbolType::parameter, mod_file_struct.parameters_within_shocks_values);
  for (auto [id, val] : std_shocks)
    val->collectVariables(SymbolType::parameter, mod_file_struct.parameters_within_shocks_values);
  for (const auto &[ids, val] : covar_shocks)
    val->collectVariables(SymbolType::parameter, mod_file_struct.parameters_within_shocks_values);
  for (const auto &[ids, val] : corr_shocks)
    val->collectVariables(SymbolType::parameter, mod_file_struct.parameters_within_shocks_values);
}

bool
ShocksStatement::has_calibrated_measurement_errors() const
{
  for (auto [id, val] : var_shocks)
    if (symbol_table.isObservedVariable(id))
      return true;

  for (auto [id, val] : std_shocks)
    if (symbol_table.isObservedVariable(id))
      return true;

  for (const auto & [ids, val] : covar_shocks)
    if (symbol_table.isObservedVariable(ids.first)
        || symbol_table.isObservedVariable(ids.second))
      return true;

  for (const auto & [ids, val] : corr_shocks)
    if (symbol_table.isObservedVariable(ids.first)
        || symbol_table.isObservedVariable(ids.second))
      return true;

  return false;
}

MShocksStatement::MShocksStatement(bool overwrite_arg,
                                   det_shocks_t det_shocks_arg,
                                   const SymbolTable &symbol_table_arg) :
  AbstractShocksStatement{true, overwrite_arg, move(det_shocks_arg), symbol_table_arg}
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

void
MShocksStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "mshocks")"
         << R"(, "overwrite": )" << boolalpha << overwrite;
  if (!det_shocks.empty())
    {
      output << ", ";
      writeJsonDetShocks(output);
    }
  output << "}";
}

ShocksSurpriseStatement::ShocksSurpriseStatement(bool overwrite_arg,
                                                 AbstractShocksStatement::det_shocks_t surprise_shocks_arg,
                                                 const SymbolTable &symbol_table_arg) :
  overwrite{overwrite_arg}, surprise_shocks{move(surprise_shocks_arg)},
  symbol_table{symbol_table_arg}
{
}

void
ShocksSurpriseStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.shocks_surprise_present = true;
}

void
ShocksSurpriseStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  if (overwrite)
    output << "M_.surprise_shocks = [" << endl;
  else
    output << "M_.surprise_shocks = [ M_.surprise_shocks;" << endl;
  for (const auto &[id, shock_vec] : surprise_shocks)
    {
      for (const auto &[period1, period2, value] : shock_vec)
        {
          output << "struct('exo_id'," << symbol_table.getTypeSpecificID(id)+1
                 << ",'periods'," << period1 << ":" << period2
                 << ",'value',";
          value->writeOutput(output);
          output << ");" << endl;
        }
    }
  output << "];" << endl;
}

void
ShocksSurpriseStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "shocks")"
         << R"(, "surprise": true)"
         << R"(, "surprise_shocks": [)";
  for (auto it = surprise_shocks.begin(); it != surprise_shocks.end(); ++it)
    {
      if (it != surprise_shocks.begin())
        output << ", ";
      output << R"({"var": ")" << symbol_table.getName(it->first) << R"(", )"
             << R"("values": [)";
      for (auto it1 = it->second.begin(); it1 != it->second.end(); ++it1)
        {
          auto [period1, period2, value] = *it1;
          if (it1 != it->second.begin())
            output << ", ";
          output << R"({"period1": )" << period1 << ", "
                 << R"("period2": )" << period2 << ", "
                 << R"("value": ")";
          value->writeJsonOutput(output, {}, {});
          output << R"("})";
        }
      output << "]}";
    }
  output << "]}";
}

ShocksLearntInStatement::ShocksLearntInStatement(int learnt_in_period_arg, bool overwrite_arg,
                                                 learnt_shocks_t learnt_shocks_arg,
                                                 const SymbolTable &symbol_table_arg) :
  learnt_in_period{learnt_in_period_arg}, overwrite{overwrite_arg},
  learnt_shocks{move(learnt_shocks_arg)}, symbol_table{symbol_table_arg}
{
}

void
ShocksLearntInStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.shocks_learnt_in_present = true;
}

string
ShocksLearntInStatement::typeToString(LearntShockType type)
{
  switch (type)
    {
    case LearntShockType::level:
      return "level";
    case LearntShockType::add:
      return "add";
    case LearntShockType::multiply:
      return "multiply";
    }
  exit(EXIT_FAILURE); // Silence GCC warning
}

void
ShocksLearntInStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  if (overwrite)
    output << "if ~isempty(M_.learnt_shocks)" << endl
           << "  M_.learnt_shocks = M_.learnt_shocks([M_.learnt_shocks.learnt_in] ~= " << learnt_in_period << ");" << endl
           << "end" << endl;

  output << "M_.learnt_shocks = [ M_.learnt_shocks;" << endl;
  for (const auto &[id, shock_vec] : learnt_shocks)
    {
      for (const auto &[type, period1, period2, value] : shock_vec)
        {
          output << "struct('learnt_in'," << learnt_in_period
                 << ",'exo_id'," << symbol_table.getTypeSpecificID(id)+1
                 << ",'periods'," << period1 << ":" << period2
                 << ",'type','" << typeToString(type) << "'"
                 << ",'value',";
          value->writeOutput(output);
          output << ");" << endl;
        }
    }
  output << "];" << endl;
}

void
ShocksLearntInStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "shocks")"
         << R"(, "learnt_in": )" << learnt_in_period
         << R"(, "overwrite": )" << boolalpha << overwrite
         << R"(, "learnt_shocks": [)";
  for (auto it = learnt_shocks.begin(); it != learnt_shocks.end(); ++it)
    {
      if (it != learnt_shocks.begin())
        output << ", ";
      output << R"({"var": ")" << symbol_table.getName(it->first) << R"(", )"
             << R"("values": [)";
      for (auto it1 = it->second.begin(); it1 != it->second.end(); ++it1)
        {
          auto [type, period1, period2, value] = *it1;
          if (it1 != it->second.begin())
            output << ", ";
          output << R"({"period1": )" << period1 << ", "
                 << R"("period2": )" << period2 << ", "
                 << R"("type": ")" << typeToString(type) << R"(", )"
                 << R"("value": ")";
          value->writeJsonOutput(output, {}, {});
          output << R"("})";
        }
      output << "]}";
    }
  output << "]}";
}

ConditionalForecastPathsStatement::ConditionalForecastPathsStatement(AbstractShocksStatement::det_shocks_t paths_arg,
                                                                     const SymbolTable &symbol_table_arg) :
  paths{move(paths_arg)},
  symbol_table{symbol_table_arg}
{
}

void
ConditionalForecastPathsStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  for (const auto &[ignore, elems] : paths)
    {
      int this_path_length = 0;
      for (auto [period1, period2, value] : elems)
        // Period1 < Period2, as enforced in ParsingDriver::add_period()
        this_path_length = max(this_path_length, period2);
      path_length = max(this_path_length, path_length);
    }
}

void
ConditionalForecastPathsStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  assert(path_length > 0);
  output << "constrained_vars_ = [];" << endl
         << "constrained_paths_ = NaN(" << paths.size() << ", " << path_length << ");" << endl;

  int k = 1;
  for (auto it = paths.begin(); it != paths.end(); ++it, k++)
    {
      if (it == paths.begin())
        output << "constrained_vars_ = " << symbol_table.getTypeSpecificID(it->first) + 1 << ";" << endl;
      else
        output << "constrained_vars_ = [constrained_vars_; " << symbol_table.getTypeSpecificID(it->first) + 1 << "];" << endl;
      for (const auto &[period1, period2, value] : it->second)
        for (int j = period1; j <= period2; j++)
          {
            output << "constrained_paths_(" << k << "," << j << ")=";
            value->writeOutput(output);
            output << ";" << endl;
          }
    }
}

void
ConditionalForecastPathsStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "conditional_forecast_paths")"
         << R"(, "paths": [)";
  for (auto it = paths.begin(); it != paths.end(); ++it)
    {
      if (it != paths.begin())
        output << ", ";
      output << R"({"var": ")" << symbol_table.getName(it->first) << R"(", )"
             << R"("values": [)";
      for (auto it1 = it->second.begin(); it1 != it->second.end(); ++it1)
        {
          auto [period1, period2, value] = *it1;
          if (it1 != it->second.begin())
            output << ", ";
          output << R"({"period1": )" << period1 << ", "
                 << R"("period2": )" << period2 << ", "
                 << R"("value": ")";
          value->writeJsonOutput(output, {}, {});
          output << R"("})";
        }
      output << "]}";
    }
  output << "]}";
}

MomentCalibration::MomentCalibration(constraints_t constraints_arg,
                                     const SymbolTable &symbol_table_arg)
  : constraints{move(constraints_arg)}, symbol_table{symbol_table_arg}
{
}

void
MomentCalibration::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "options_.endogenous_prior_restrictions.moment = {" << endl;
  for (const auto &c : constraints)
    {
      output << "'" << symbol_table.getName(c.endo1) << "', "
             << "'" << symbol_table.getName(c.endo2) << "', "
             << c.lags << ", "
             << "[ ";
      c.lower_bound->writeOutput(output);
      output << ", ";
      c.upper_bound->writeOutput(output);
      output << " ];"
             << endl;
    }
  output << "};" << endl;
}

void
MomentCalibration::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "moment_calibration")"
         << R"(, "moment_calibration_criteria": [)";
  for (auto it = constraints.begin(); it != constraints.end(); ++it)
    {
      if (it != constraints.begin())
        output << ", ";
      output << R"({"endogenous1": ")" << symbol_table.getName(it->endo1) << R"(")"
             << R"(, "endogenous2": ")" << symbol_table.getName(it->endo2) << R"(")"
             << R"(, "lags": ")" << it->lags << R"(")"
             << R"(, "lower_bound": ")";
      it->lower_bound->writeJsonOutput(output, {}, {});
      output << R"(")"
             << R"(, "upper_bound": ")";
      it->upper_bound->writeJsonOutput(output, {}, {});
      output << R"(")"
             << "}";
    }
  output << "]"
         << "}";
}

IrfCalibration::IrfCalibration(constraints_t constraints_arg,
                               const SymbolTable &symbol_table_arg,
                               OptionsList options_list_arg)
  : constraints{move(constraints_arg)}, symbol_table{symbol_table_arg}, options_list{move(options_list_arg)}
{
}

void
IrfCalibration::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  options_list.writeOutput(output);

  output << "options_.endogenous_prior_restrictions.irf = {" << endl;
  for (const auto &c : constraints)
    {
      output << "'" << symbol_table.getName(c.endo) << "', "
             << "'" << symbol_table.getName(c.exo) << "', "
             << c.periods << ", "
             << "[ ";
      c.lower_bound->writeOutput(output);
      output << ", ";
      c.upper_bound->writeOutput(output);
      output << " ];"
             << endl;
    }
  output << "};" << endl;
}

void
IrfCalibration::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "irf_calibration")";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }

  output << R"(, "irf_restrictions": [)";
  for (auto it = constraints.begin(); it != constraints.end(); ++it)
    {
      if (it != constraints.begin())
        output << ", ";
      output << R"({"endogenous": ")" << symbol_table.getName(it->endo) << R"(")"
             << R"(, "exogenous": ")" << symbol_table.getName(it->exo) << R"(")"
             << R"(, "periods": ")" << it->periods << R"(")"
             << R"(, "lower_bound": ")";
      it->lower_bound->writeJsonOutput(output, {}, {});
      output << R"(")";
      output << R"(, "upper_bound": ")";
      it->upper_bound->writeJsonOutput(output, {}, {});
      output << R"(")"
             << "}";
    }
  output << "]"
         << "}";
}

ShockGroupsStatement::ShockGroupsStatement(group_t shock_groups_arg, string name_arg)
  : shock_groups{move(shock_groups_arg)}, name{move(name_arg)}
{
}

void
ShockGroupsStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  int i = 1;
  bool unique_label = true;
  for (auto it = shock_groups.begin(); it != shock_groups.end(); ++it, unique_label = true)
    {
      for (auto it1 = it+1; it1 != shock_groups.end(); ++it1)
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
          for (const auto &it1 : it->list)
            output << " '" << it1 << "'";
          output << "};" << endl;
          i++;
        }
    }
}

void
ShockGroupsStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "shock_groups", "name": ")" << name << R"(", "groups": [)";
  bool unique_label = true;
  bool printed_group = false;
  for (auto it = shock_groups.begin(); it != shock_groups.end(); ++it, unique_label = true)
    {
      for (auto it1 = it+1; it1 != shock_groups.end(); ++it1)
        if (it->name == it1->name)
          {
            unique_label = false;
            break;
          }

      if (unique_label)
        {
          if (printed_group)
            output << ", ";
          else
            printed_group = true;
          output << R"({"group_name": ")" << it->name << R"(",)"
                 << R"("shocks": [)";
          for (auto it1 = it->list.begin(); it1 != it->list.end(); ++it1)
            {
              if (it1 != it->list.begin())
                output << ", ";
              output << R"(")" << *it1 << R"(")";
            }
          output << "]}";
        }
    }
  output << "]}";
}

Init2shocksStatement::Init2shocksStatement(vector<pair<int, int>> init2shocks_arg, string name_arg,
                                           const SymbolTable &symbol_table_arg)
  : init2shocks{move(init2shocks_arg)}, name{move(name_arg)}, symbol_table{symbol_table_arg}
{
}

void
Init2shocksStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  for (size_t i = 0; i < init2shocks.size(); i++)
    for (size_t j = i + 1; j < init2shocks.size(); j++)
      if (init2shocks.at(i).first == init2shocks.at(j).first)
        {
          cerr << "Init2shocks(" << name << "): enogenous variable '"
               << symbol_table.getName(init2shocks.at(i).first)
               << "' appears more than once in the init2shocks statement" << endl;
          exit(EXIT_FAILURE);
        }
}

void
Init2shocksStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "M_.init2shocks." << name << " = {" << endl;
  for (auto &it : init2shocks)
    output << "{'" << symbol_table.getName(it.first) << "', '" << symbol_table.getName(it.second) << "'};" << endl;
  output << "};" << endl;
}

void
Init2shocksStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "init2shocks", "name": ")" << name << R"(", "groups": [)";
  for (auto &it : init2shocks)
    {
      if (it != *(init2shocks.begin()))
        output << ",";
      output << R"({"endogenous": ")" << symbol_table.getName(it.first) << R"(", )"
             << R"( "exogenous": ")" << symbol_table.getName(it.second) << R"("})";
    }
  output << "]}";
}

HeteroskedasticShocksStatement::HeteroskedasticShocksStatement(bool overwrite_arg,
                                                               const heteroskedastic_shocks_t &values_arg,
                                                               const heteroskedastic_shocks_t &scales_arg,
                                                               const SymbolTable &symbol_table_arg)
  : overwrite{overwrite_arg}, values{values_arg}, scales{scales_arg}, symbol_table{symbol_table_arg}
{
}

void
HeteroskedasticShocksStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  // NB: The first initialization of the fields is done in ModFile::writeMOutput()
  if (overwrite)
    output << "M_.heteroskedastic_shocks.Qvalue_orig = [];" << endl
           << "M_.heteroskedastic_shocks.Qscale_orig = [];" << endl;

  for (const auto &[symb_id, vec] : values)
    {
      int tsid = symbol_table.getTypeSpecificID(symb_id);
      for (const auto &[period1, period2, value] : vec)
        {
          output << "M_.heteroskedastic_shocks.Qvalue_orig = [M_.heteroskedastic_shocks.Qvalue_orig; struct('exo_id', "
                 << tsid+1 << ",'periods',"
                 << period1 << ":" << period2 << ",'value',";
          value->writeOutput(output);
          output << ")];" << endl;
        }
    }
  for (const auto &[symb_id, vec] : scales)
    {
      int tsid = symbol_table.getTypeSpecificID(symb_id);
      for (const auto &[period1, period2, scale] : vec)
        {
          output << "M_.heteroskedastic_shocks.Qscale_orig = [M_.heteroskedastic_shocks.Qscale_orig; struct('exo_id', "
                 << tsid+1 << ",'periods',"
                 << period1 << ":" << period2 << ",'scale',";
          scale->writeOutput(output);
          output << ")];" << endl;
        }
    }
}

void
HeteroskedasticShocksStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "heteroskedastic_shocks")"
         << R"(, "overwrite": )" << boolalpha << overwrite
         << R"(, "shocks_values": [)";
  for (auto it = values.begin(); it != values.end(); ++it)
    {
      if (it != values.begin())
        output << ", ";
      output << R"({"var": ")" << symbol_table.getName(it->first) << R"(", )"
             << R"("values": [)";
      for (auto it1 = it->second.begin(); it1 != it->second.end(); ++it1)
        {
          if (it1 != it->second.begin())
            output << ", ";
          auto [period1, period2, value] = *it1;
          output << R"({"period1": )" << period1 << ", "
                 << R"("period2": )" << period2 << ", "
                 << R"("value": ")";
          value->writeJsonOutput(output, {}, {});
          output << R"("})";
        }
      output << "]}";
    }
  output << R"(], "shocks_scales": [)";
  for (auto it = scales.begin(); it != scales.end(); ++it)
    {
      if (it != scales.begin())
        output << ", ";
      output << R"({"var": ")" << symbol_table.getName(it->first) << R"(", )"
             << R"("scales": [)";
      for (auto it1 = it->second.begin(); it1 != it->second.end(); ++it1)
        {
          if (it1 != it->second.begin())
            output << ", ";
          auto [period1, period2, value] = *it1;
          output << R"({"period1": )" << period1 << ", "
                 << R"("period2": )" << period2 << ", "
                 << R"("value": ")";
          value->writeJsonOutput(output, {}, {});
          output << R"("})";
        }
      output << "]}";
    }
  output << "]}";
}
