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
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <utility>

#include "NumericalInitialization.hh"

InitParamStatement::InitParamStatement(int symb_id_arg,
                                       const expr_t param_value_arg,
                                       const SymbolTable &symbol_table_arg) :
  symb_id{symb_id_arg},
  param_value{param_value_arg},
  symbol_table{symbol_table_arg}
{
}

void
InitParamStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  if (symbol_table.getName(symb_id) == "dsge_prior_weight")
    mod_file_struct.dsge_prior_weight_initialized = true;
}

void
InitParamStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  int id = symbol_table.getTypeSpecificID(symb_id) + 1;
  output << "M_.params(" << id << ") = ";
  param_value->writeOutput(output);
  output << ";" << endl;
  if (!minimal_workspace)
    output << symbol_table.getName(symb_id) << " = M_.params(" << id << ");" << endl;
}

void
InitParamStatement::writeJuliaOutput(ostream &output, const string &basename)
{
  int id = symbol_table.getTypeSpecificID(symb_id) + 1;
  output << "model_.params[ " << id << " ] = ";
  param_value->writeOutput(output);
  output << endl;
}

void
InitParamStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "param_init", "name": ")" << symbol_table.getName(symb_id) << R"(", )" << R"("value": ")";
  param_value->writeJsonOutput(output, {}, {});
  output << R"("})";
}

void
InitParamStatement::fillEvalContext(eval_context_t &eval_context) const
{
  try
    {
      eval_context[symb_id] = param_value->eval(eval_context);
    }
  catch (ExprNode::EvalException &e)
    {
      // Do nothing
    }
}

InitOrEndValStatement::InitOrEndValStatement(init_values_t init_values_arg,
                                             const SymbolTable &symbol_table_arg,
                                             bool all_values_required_arg) :
  init_values{move(init_values_arg)},
  symbol_table{symbol_table_arg},
  all_values_required{all_values_required_arg}
{
}

void
InitOrEndValStatement::fillEvalContext(eval_context_t &eval_context) const
{
  for (const auto & init_value : init_values)
    {
      try
        {
          eval_context[init_value.first] = (init_value.second)->eval(eval_context);
        }
      catch (ExprNode::EvalException &e)
        {
          // Do nothing
        }
    }
}

set<int>
InitOrEndValStatement::getUninitializedVariables(SymbolType type)
{
  set<int> unused;
  if (!all_values_required)
    return unused;

  if (type == SymbolType::endogenous)
    unused = symbol_table.getEndogenous();
  else if (type == SymbolType::exogenous)
    unused = symbol_table.getExogenous();
  else
    {
      cerr << "ERROR: Shouldn't arrive here." << endl;
      exit(EXIT_FAILURE);
    }

  for (const auto & init_value : init_values)
    if (auto sit = unused.find(init_value.first);
        sit != unused.end())
      unused.erase(sit);

  return unused;
}

void
InitOrEndValStatement::writeInitValues(ostream &output) const
{
  for (const auto & init_value : init_values)
    {
      const int symb_id = init_value.first;
      const expr_t expression = init_value.second;

      SymbolType type = symbol_table.getType(symb_id);
      int tsid = symbol_table.getTypeSpecificID(symb_id) + 1;

      switch (type)
        {
        case SymbolType::endogenous:
          output << "oo_.steady_state";
          break;
        case SymbolType::exogenous:
          output << "oo_.exo_steady_state";
          break;
        case SymbolType::exogenousDet:
          output << "oo_.exo_det_steady_state";
          break;
        case SymbolType::excludedVariable:
          cerr << "ERROR: Variable `" << symbol_table.getName(symb_id)
               << "` was excluded but found in an initval or endval statement" << endl;
          exit(EXIT_FAILURE);
        default:
          cerr << "Should not arrive here" << endl;
          exit(EXIT_FAILURE);
        }

      output << "(" << tsid << ") = ";
      expression->writeOutput(output);
      output << ";" << endl;
    }
}

void
InitOrEndValStatement::writeJsonInitValues(ostream &output) const
{
  for (auto it = init_values.begin();
       it != init_values.end(); it++)
    {
      if (it != init_values.begin())
        output << ", ";
      output << R"({"name": ")" << symbol_table.getName(it->first) << R"(", )" << R"("value": ")";
      it->second->writeJsonOutput(output, {}, {});
      output << R"("})";
    }
}

InitValStatement::InitValStatement(const init_values_t &init_values_arg,
                                   const SymbolTable &symbol_table_arg,
                                   bool all_values_required_arg) :
  InitOrEndValStatement{init_values_arg, symbol_table_arg, all_values_required_arg}
{
}

void
InitValStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  set<int> exogs = getUninitializedVariables(SymbolType::exogenous);
  set<int> endogs = getUninitializedVariables(SymbolType::endogenous);

  if (endogs.size() > 0)
    {
      cerr << "ERROR: You have not set the following endogenous variables in initval:";
      for (int endog : endogs)
        cerr << " " << symbol_table.getName(endog);
      cerr << endl;
    }

  if (exogs.size() > 0)
    {
      cerr << "ERROR: You have not set the following exogenous variables in initval:";
      for (int exog : exogs)
        cerr << " " << symbol_table.getName(exog);
      cerr << endl;
    }

  if (endogs.size() > 0 || exogs.size() > 0)
    exit(EXIT_FAILURE);
}

void
InitValStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "%" << endl
         << "% INITVAL instructions" << endl
         << "%" << endl;
  // Writing initval block to set initial values for variables
  output << "options_.initval_file = false;" << endl;

  writeInitValues(output);
}

void
InitValStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "initval", "vals": [)";
  writeJsonInitValues(output);
  output << "]}";
}

void
InitValStatement::writeOutputPostInit(ostream &output) const
{
  output << "if M_.exo_nbr > 0" << endl
         << "\too_.exo_simul = ones(M_.maximum_lag,1)*oo_.exo_steady_state';" << endl
         <<"end" << endl
         << "if M_.exo_det_nbr > 0" << endl
         << "\too_.exo_det_simul = ones(M_.maximum_lag,1)*oo_.exo_det_steady_state';" << endl
         <<"end" << endl;
}

EndValStatement::EndValStatement(const init_values_t &init_values_arg,
                                 const SymbolTable &symbol_table_arg,
                                 bool all_values_required_arg) :
  InitOrEndValStatement{init_values_arg, symbol_table_arg, all_values_required_arg}
{
}

void
EndValStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  set<int> exogs = getUninitializedVariables(SymbolType::exogenous);
  set<int> endogs = getUninitializedVariables(SymbolType::endogenous);

  if (endogs.size() > 0)
    {
      cerr << "ERROR: You have not set the following endogenous variables in endval:";
      for (int endog : endogs)
        cerr << " " << symbol_table.getName(endog);
      cerr << endl;
    }

  if (exogs.size() > 0)
    {
      cerr << "ERROR: You have not set the following exogenous variables in endval:";
      for (int exog : exogs)
        cerr << " " << symbol_table.getName(exog);
      cerr << endl;
    }

  if (endogs.size() > 0 || exogs.size() > 0)
    exit(EXIT_FAILURE);
}

void
EndValStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "%" << endl
         << "% ENDVAL instructions" << endl
         << "%" << endl;
  // Writing endval block to set terminal values for variables
  output << "ys0_= oo_.steady_state;" << endl
         << "ex0_ = oo_.exo_steady_state;" << endl;

  writeInitValues(output);
}

void
EndValStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "endval", "vals": [)";
  writeJsonInitValues(output);
  output << "]}";
}

HistValStatement::HistValStatement(hist_values_t hist_values_arg,
                                   const SymbolTable &symbol_table_arg,
                                   bool all_values_required_arg) :
  hist_values{move(hist_values_arg)},
  symbol_table{symbol_table_arg},
  all_values_required{all_values_required_arg}
{
}

void
HistValStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  if (all_values_required)
    {
      set<int> unused_endo = symbol_table.getEndogenous();
      set<int> unused_exo = symbol_table.getExogenous();

      for (const auto & hist_value : hist_values)
        {
          if (auto sit = unused_endo.find(hist_value.first.first);
              sit != unused_endo.end())
            unused_endo.erase(sit);

          if (auto sit = unused_exo.find(hist_value.first.first);
              sit != unused_exo.end())
            unused_exo.erase(sit);
        }

      if (unused_endo.size() > 0)
        {
          cerr << "ERROR: You have not set the following endogenous variables in histval:";
          for (int it : unused_endo)
            cerr << " " << symbol_table.getName(it);
          cerr << endl;
        }

      if (unused_exo.size() > 0)
        {
          cerr << "ERROR: You have not set the following exogenous variables in endval:";
          for (int it : unused_exo)
            cerr << " " << symbol_table.getName(it);
          cerr << endl;
        }

      if (unused_endo.size() > 0 || unused_exo.size() > 0)
        exit(EXIT_FAILURE);
    }
}

void
HistValStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "%" << endl
         << "% HISTVAL instructions" << endl
         << "%" << endl
         << "M_.histval_dseries = dseries(zeros(M_.orig_maximum_lag_with_diffs_expanded, M_.orig_endo_nbr"
         << (symbol_table.AuxVarsSize() > 0 ? "+sum([M_.aux_vars.type]==6)" : "")
         << (symbol_table.exo_nbr() > 0 ? "+M_.exo_nbr" : "")
         << (symbol_table.exo_det_nbr() > 0 ? "+M_.exo_det_nbr" : "")
         << "), dates(sprintf('%dY', -M_.orig_maximum_lag_with_diffs_expanded+1)), [ M_.endo_names(1:M_.orig_endo_nbr); "
         << (symbol_table.AuxVarsSize() > 0 ? "M_.endo_names([M_.aux_vars(find([M_.aux_vars.type]==6)).endo_index]); " : "")
         << (symbol_table.exo_nbr() > 0 ? "M_.exo_names; " : "")
         << (symbol_table.exo_det_nbr() > 0 ? "M_.exo_det_names; " : "")
         << "]);" << endl;

  for (const auto & hist_value : hist_values)
    {
      int symb_id = hist_value.first.first;
      int lag = hist_value.first.second;
      const expr_t expression = hist_value.second;

      output << "M_.histval_dseries{'" << symbol_table.getName(symb_id) << "'}(dates('" << lag << "Y'))=";
      expression->writeOutput(output);
      output << ";" << endl;
    }

  output << "if exist(['+' M_.fname '/dynamic_set_auxiliary_series.m'])" << endl
       << "  eval(['M_.histval_dseries = ' M_.fname '.dynamic_set_auxiliary_series(M_.histval_dseries, M_.params);']);" << endl
       << "end" << endl
         << "M_.endo_histval = M_.histval_dseries{M_.endo_names{:}}(dates(sprintf('%dY', 1-M_.maximum_lag)):dates('0Y')).data';" << endl
         << "M_.endo_histval(isnan(M_.endo_histval)) = 0;" << endl; // Ensure that lead aux variables do not have a NaN

  if (symbol_table.exo_nbr() > 0)
    output << "M_.exo_histval = M_.histval_dseries{M_.exo_names{:}}(dates(sprintf('%dY', 1-M_.maximum_lag)):dates('0Y')).data';" << endl;
  if (symbol_table.exo_det_nbr() > 0)
    output << "M_.exo_det_histval = M_.histval_dseries{M_.exo_det_names{:}}(dates(sprintf('%dY', 1-M_.maximum_lag)):dates('0Y')).data';" << endl;
}

void
HistValStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "histval", "vals": [)";
  for (auto it = hist_values.begin();
       it != hist_values.end(); ++it)
    {
      if (it != hist_values.begin())
        output << ", ";
      output << R"({ "name": ")" << symbol_table.getName(it->first.first) << R"(")"
             << R"(, "lag": )" << it->first.second
             << R"(, "value": ")";
      it->second->writeJsonOutput(output, {}, {});
      output << R"("})";
    }
  output << "]}";
}

InitvalFileStatement::InitvalFileStatement(string filename_arg) :
  filename{move(filename_arg)}
{
}

void
InitvalFileStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "%" << endl
         << "% INITVAL_FILE statement" << endl
         << "%" << endl
         << "options_.initval_file = true;" << endl
         << "initvalf('" << filename << "');" << endl;
}

void
InitvalFileStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "init_val_file")"
         << R"(, "filename": ")" << filename << R"(")"
         << "}";
}

HistvalFileStatement::HistvalFileStatement(string filename_arg) :
  filename{move(filename_arg)}
{
}

void
HistvalFileStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "histvalf('" << filename << "');" << endl;
}

void
HistvalFileStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "hist_val_file")"
         << R"(, "filename": ")" << filename << R"(")"
         << "}";
}

HomotopyStatement::HomotopyStatement(homotopy_values_t homotopy_values_arg,
                                     const SymbolTable &symbol_table_arg) :
  homotopy_values{move(homotopy_values_arg)},
  symbol_table{symbol_table_arg}
{
}

void
HomotopyStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "%" << endl
         << "% HOMOTOPY_SETUP instructions" << endl
         << "%" << endl
         << "options_.homotopy_values = [];" << endl;

  for (auto [symb_id, expression1, expression2] : homotopy_values)
    {
      const SymbolType type = symbol_table.getType(symb_id);
      const int tsid = symbol_table.getTypeSpecificID(symb_id) + 1;

      output << "options_.homotopy_values = vertcat(options_.homotopy_values, [ " << static_cast<int>(type) << ", " << tsid << ", ";
      if (expression1)
        expression1->writeOutput(output);
      else
        output << "NaN";
      output << ", ";
      expression2->writeOutput(output);
      output << "]);" << endl;
    }
}

void
HomotopyStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "homotopy", )"
         << R"("values": [)";
  for (auto it = homotopy_values.begin();
       it != homotopy_values.end(); ++it)
    {
      if (it != homotopy_values.begin())
        output << ", ";

      auto [symb_id, expression1, expression2] = *it;

      output << R"({"name": ")" << symbol_table.getName(symb_id) << R"(")"
             << R"(, "initial_value": ")";
      if (expression1)
        expression1->writeJsonOutput(output, {}, {});
      else
        output << "NaN";
      output << R"(", "final_value": ")";
      expression2->writeJsonOutput(output, {}, {});
      output << R"("})";
    }
  output << "]"
         << "}";
}

SaveParamsAndSteadyStateStatement::SaveParamsAndSteadyStateStatement(string filename_arg) :
  filename{move(filename_arg)}
{
}

void
SaveParamsAndSteadyStateStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "save_params_and_steady_state('" << filename << "');" << endl;
}

void
SaveParamsAndSteadyStateStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "save_params_and_steady_state")"
         << R"(, "filename": ")" << filename << R"(")"
         << "}";
}

LoadParamsAndSteadyStateStatement::LoadParamsAndSteadyStateStatement(const string &filename,
                                                                     const SymbolTable &symbol_table_arg,
                                                                     WarningConsolidation &warnings) :
  symbol_table{symbol_table_arg}
{
  cout << "Reading " << filename << "." << endl;

  ifstream f;
  f.open(filename, ios::in);
  if (f.fail())
    {
      cerr << "ERROR: Can't open " << filename << endl;
      exit(EXIT_FAILURE);
    }

  while (true)
    {
      string symb_name, value;
      f >> symb_name >> value;
      if (f.eof())
        break;

      try
        {
          int symb_id = symbol_table.getID(symb_name);
          content[symb_id] = value;
        }
      catch (SymbolTable::UnknownSymbolNameException &e)
        {
          warnings << "WARNING: Unknown symbol " << symb_name << " in " << filename << endl;
        }
    }
  f.close();
}

void
LoadParamsAndSteadyStateStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  for (const auto &[id, value] : content)
    {
      switch (symbol_table.getType(id))
        {
        case SymbolType::parameter:
          output << "M_.params";
          break;
        case SymbolType::endogenous:
          output << "oo_.steady_state";
          break;
        case SymbolType::exogenous:
          output << "oo_.exo_steady_state";
          break;
        case SymbolType::exogenousDet:
          output << "oo_.exo_det_steady_state";
          break;
        default:
          cerr << "ERROR: Unsupported variable type for " << symbol_table.getName(id) << " in load_params_and_steady_state" << endl;
          exit(EXIT_FAILURE);
        }

      int tsid = symbol_table.getTypeSpecificID(id) + 1;
      output << "(" << tsid << ") = " << value << ";" << endl;
    }
}

void
LoadParamsAndSteadyStateStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "load_params_and_steady_state")"
         << R"("values": [)";
  for (auto it = content.begin(); it != content.end(); ++it)
    {
      if (it != content.begin())
        output << ", ";
      output << R"({"name": ")" << symbol_table.getName(it->first) << R"(")"
             << R"(, "value": ")" << it->second << R"("})";
    }
  output << "]"
         << "}";
}

void
LoadParamsAndSteadyStateStatement::fillEvalContext(eval_context_t &eval_context) const
{
  for (const auto & [id, value] : content)
    /* We use strtod() instead of stod() because we want underflows and
       overflows to respectively yield 0 and ±Inf. See also the comment in
       NumericalConstants.cc */
    eval_context[id] = strtod(value.c_str(), nullptr);
}
