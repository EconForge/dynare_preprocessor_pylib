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
InitParamStatement::checkPass(ModFileStructure &mod_file_struct, [[maybe_unused]] WarningConsolidation &warnings)
{
  if (symbol_table.getName(symb_id) == "dsge_prior_weight")
    mod_file_struct.dsge_prior_weight_initialized = true;

  // Needed for the workaround discussed in dynare#1173
  if (symbol_table.getName(symb_id) == "optimal_policy_discount_factor")
    param_value->collectVariables(SymbolType::parameter, mod_file_struct.parameters_in_planner_discount);
}

void
InitParamStatement::writeOutput(ostream &output, [[maybe_unused]] const string &basename, bool minimal_workspace) const
{
  int id = symbol_table.getTypeSpecificID(symb_id) + 1;
  output << "M_.params(" << id << ") = ";
  param_value->writeOutput(output);
  output << ";" << endl;
  if (!minimal_workspace)
    output << symbol_table.getName(symb_id) << " = M_.params(" << id << ");" << endl;
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
  for (auto [symb_id, value] : init_values)
    try
      {
        eval_context[symb_id] = value->eval(eval_context);
      }
    catch (ExprNode::EvalException &e)
      {
        // Do nothing
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

  for (auto [symb_id, value] : init_values)
    unused.erase(symb_id);

  return unused;
}

void
InitOrEndValStatement::writeInitValues(ostream &output) const
{
  for (auto [symb_id, value] : init_values)
    {
      if (symbol_table.getType(symb_id) == SymbolType::unusedEndogenous) // See #82
        continue;

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
      value->writeOutput(output);
      output << ";" << endl;
    }
}

void
InitOrEndValStatement::writeJsonInitValues(ostream &output) const
{
  for (bool printed_something{false};
       auto &[symb_id, value] : init_values)
    {
      if (symbol_table.getType(symb_id) == SymbolType::unusedEndogenous) // See #82
        continue;
      if (exchange(printed_something, true))
        output << ", ";
      output << R"({"name": ")" << symbol_table.getName(symb_id) << R"(", )" << R"("value": ")";
      value->writeJsonOutput(output, {}, {});
      output << R"("})";
    }
}

InitValStatement::InitValStatement(init_values_t init_values_arg,
                                   const SymbolTable &symbol_table_arg,
                                   bool all_values_required_arg) :
  InitOrEndValStatement{move(init_values_arg), symbol_table_arg, all_values_required_arg}
{
}

void
InitValStatement::checkPass([[maybe_unused]] ModFileStructure &mod_file_struct,
                            [[maybe_unused]] WarningConsolidation &warnings)
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
InitValStatement::writeOutput(ostream &output, [[maybe_unused]] const string &basename,
                              [[maybe_unused]] bool minimal_workspace) const
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

EndValStatement::EndValStatement(init_values_t init_values_arg,
                                 const SymbolTable &symbol_table_arg,
                                 bool all_values_required_arg) :
  InitOrEndValStatement{move(init_values_arg), symbol_table_arg, all_values_required_arg}
{
}

void
EndValStatement::checkPass([[maybe_unused]] ModFileStructure &mod_file_struct,
                           [[maybe_unused]] WarningConsolidation &warnings)
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
EndValStatement::writeOutput(ostream &output, [[maybe_unused]] const string &basename,
                             [[maybe_unused]] bool minimal_workspace) const
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

EndValLearntInStatement::EndValLearntInStatement(int learnt_in_period_arg,
                                                 learnt_end_values_t learnt_end_values_arg,
                                                 const SymbolTable &symbol_table_arg) :
  learnt_in_period{learnt_in_period_arg},
  learnt_end_values{move(learnt_end_values_arg)},
  symbol_table{symbol_table_arg}
{
}

void
EndValLearntInStatement::checkPass(ModFileStructure &mod_file_struct,
                                   [[maybe_unused]] WarningConsolidation &warnings)
{
  mod_file_struct.endval_learnt_in_present = true;
}

string
EndValLearntInStatement::typeToString(LearntEndValType type)
{
  switch (type)
    {
    case LearntEndValType::level:
      return "level";
    case LearntEndValType::add:
      return "add";
    case LearntEndValType::multiply:
      return "multiply";
    }
  exit(EXIT_FAILURE); // Silence GCC warning
}

void
EndValLearntInStatement::writeOutput(ostream &output, [[maybe_unused]] const string &basename,
                                     [[maybe_unused]] bool minimal_workspace) const
{
  output << "M_.learnt_endval = [ M_.learnt_endval;" << endl;
  for (auto [type, symb_id, value] : learnt_end_values)
    {
      if (symbol_table.getType(symb_id) == SymbolType::unusedEndogenous) // See #82
        continue;
      output << "struct('learnt_in'," << learnt_in_period
             << ",'exo_id'," << symbol_table.getTypeSpecificID(symb_id)+1
             << ",'type','" << typeToString(type) << "'"
             << ",'value',";
      value->writeOutput(output);
      output << ");" << endl;
    }
  output << "];" << endl;
}

void
EndValLearntInStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "endval", "learnt_in": )"
         << learnt_in_period <<  R"(, "vals": [)";
  for (bool printed_something{false};
       auto &[type, symb_id, value] : learnt_end_values)
    {
      if (symbol_table.getType(symb_id) == SymbolType::unusedEndogenous) // See #82
        continue;
      if (exchange(printed_something, true))
        output << ", ";
      output << R"({"name": ")" << symbol_table.getName(symb_id) << R"(", )"
             << R"("type": ")" << typeToString(type) << R"(", )"
             << R"("value": ")";
      value->writeJsonOutput(output, {}, {});
      output << R"("})";
    }
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
HistValStatement::checkPass([[maybe_unused]] ModFileStructure &mod_file_struct,
                            [[maybe_unused]] WarningConsolidation &warnings)
{
  if (all_values_required)
    {
      set<int> unused_endo = symbol_table.getEndogenous();
      set<int> unused_exo = symbol_table.getExogenous();

      for (const auto &[key, value] : hist_values)
        {
          int symb_id = key.first;
          unused_endo.erase(symb_id);
          unused_exo.erase(symb_id);
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
HistValStatement::writeOutput(ostream &output, [[maybe_unused]] const string &basename,
                              [[maybe_unused]] bool minimal_workspace) const
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

  for (const auto &[key, value] : hist_values)
    {
      auto [symb_id, lag] = key;
      if (symbol_table.getType(symb_id) == SymbolType::unusedEndogenous) // See #82
        continue;

      output << "M_.histval_dseries{'" << symbol_table.getName(symb_id) << "'}(dates('" << lag << "Y'))=";
      value->writeOutput(output);
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
  for (bool printed_something{false};
       const auto &[key, value] : hist_values)
    {
      auto &[symb_id, lag] = key;
      if (symbol_table.getType(symb_id) == SymbolType::unusedEndogenous) // See #82
        continue;
      if (exchange(printed_something, true))
        output << ", ";
      output << R"({ "name": ")" << symbol_table.getName(symb_id) << R"(")"
             << R"(, "lag": )" << lag
             << R"(, "value": ")";
      value->writeJsonOutput(output, {}, {});
      output << R"("})";
    }
  output << "]}";
}

InitvalFileStatement::InitvalFileStatement(OptionsList options_list_arg) :
  options_list{move(options_list_arg)}
{
}

void
InitvalFileStatement::writeOutput(ostream &output, [[maybe_unused]] const string &basename,
                                  [[maybe_unused]] bool minimal_workspace) const
{
  output << "%" << endl
         << "% INITVAL_FILE statement" << endl
         << "%" << endl
         << "options_.initval_file = true;" << endl;
  options_list.writeOutput(output, "options_initvalf");
  output << "[oo_.initval_series, options_.periods] = initvalf(M_, options_initvalf);" << endl;
}

void
InitvalFileStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "initval_file")";
  if (!options_list.empty())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

HistvalFileStatement::HistvalFileStatement(OptionsList options_list_arg) :
  options_list{move(options_list_arg)}
{
}

void
HistvalFileStatement::writeOutput(ostream &output, [[maybe_unused]] const string &basename,
                                  [[maybe_unused]] bool minimal_workspace) const
{
  output << "%" << endl
         << "% HISTVAL_FILE statement" << endl
         << "%" << endl
         << "options_.histval_file = true;" << endl;
  options_list.writeOutput(output, "options_histvalf");
  output << "[M_.endo_histval, M_.exo_histval, M_.exo_det_histval] = histvalf(M_, options_histvalf);" << endl;
}

void
HistvalFileStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "histval_file")";
  if (!options_list.empty())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

HomotopyStatement::HomotopyStatement(homotopy_values_t homotopy_values_arg,
                                     const SymbolTable &symbol_table_arg) :
  homotopy_values{move(homotopy_values_arg)},
  symbol_table{symbol_table_arg}
{
}

void
HomotopyStatement::writeOutput(ostream &output, [[maybe_unused]] const string &basename,
                               [[maybe_unused]] bool minimal_workspace) const
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
  for (bool printed_something{false};
       const auto &[symb_id, expression1, expression2] : homotopy_values)
    {
      if (exchange(printed_something, true))
        output << ", ";

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
SaveParamsAndSteadyStateStatement::writeOutput(ostream &output, [[maybe_unused]] const string &basename,
                                               [[maybe_unused]] bool minimal_workspace) const
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

LoadParamsAndSteadyStateStatement::LoadParamsAndSteadyStateStatement(const filesystem::path &filename,
                                                                     const SymbolTable &symbol_table_arg,
                                                                     WarningConsolidation &warnings) :
  symbol_table{symbol_table_arg}
{
  cout << "Reading " << filename.string() << "." << endl;

  ifstream f;
  f.open(filename, ios::in);
  if (f.fail())
    {
      cerr << "ERROR: Can't open " << filename.string() << endl;
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
          warnings << "WARNING: Unknown symbol " << symb_name << " in " << filename.string() << endl;
        }
    }
  f.close();
}

void
LoadParamsAndSteadyStateStatement::writeOutput(ostream &output, [[maybe_unused]] const string &basename,
                                               [[maybe_unused]] bool minimal_workspace) const
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
  output << R"({"statementName": "load_params_and_steady_state",)"
         << R"("values": [)";
  for (bool printed_something{false};
       const auto &[id, value] : content)
    {
      if (exchange(printed_something, true))
        output << ", ";
      output << R"({"name": ")" << symbol_table.getName(id) << R"(")"
             << R"(, "value": ")" << value << R"("})";
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
