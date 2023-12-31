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

#include <fstream>
#include <iostream>
#include <cassert>
#include <sstream>
#include <cmath>
#include <numeric>

#include "ParsingDriver.hh"
#include "Statement.hh"
#include "ExprNode.hh"
#include "WarningConsolidation.hh"

bool
ParsingDriver::symbol_exists_and_is_not_modfile_local_or_external_function(const string &s)
{
  if (!mod_file->symbol_table.exists(s))
    return false;

  SymbolType type = mod_file->symbol_table.getType(s);

  return type != SymbolType::modFileLocalVariable && type != SymbolType::externalFunction;
}

void
ParsingDriver::check_symbol_existence_in_model_block(const string &name)
{
  if (!mod_file->symbol_table.exists(name) || undeclared_model_vars.contains(name))
    undeclared_model_variable_error("Unknown symbol: " + name, name);
}

void
ParsingDriver::check_symbol_existence(const string &name)
{
  if (!mod_file->symbol_table.exists(name))
    error("Unknown symbol: " + name + ".\nIf referenced from the 'initval', 'endval', 'histval', or 'shocks' block, you can pass the 'nostrict' option to dynare to have this line ignored.");
}

void
ParsingDriver::check_symbol_is_parameter(const string &name)
{
  check_symbol_existence(name);
  int symb_id = mod_file->symbol_table.getID(name);
  if (mod_file->symbol_table.getType(symb_id) != SymbolType::parameter)
    error(name + " is not a parameter");
}

void
ParsingDriver::set_current_data_tree(DataTree *data_tree_arg)
{
  data_tree = data_tree_arg;
  model_tree = dynamic_cast<ModelTree *>(data_tree_arg);
  dynamic_model = dynamic_cast<DynamicModel *>(data_tree_arg);
}

void
ParsingDriver::reset_data_tree()
{
  set_current_data_tree(&mod_file->expressions_tree);
}

void
ParsingDriver::reset_current_external_function_options()
{
  current_external_function_options.nargs = ExternalFunctionsTable::defaultNargs;
  current_external_function_options.firstDerivSymbID = ExternalFunctionsTable::IDNotSet;
  current_external_function_options.secondDerivSymbID = ExternalFunctionsTable::IDNotSet;
  current_external_function_id = ExternalFunctionsTable::IDNotSet;
}

unique_ptr<ModFile>
ParsingDriver::parse(istream &in, bool debug)
{
  mod_file = make_unique<ModFile>(warnings);

  reset_data_tree();
  estim_params.init(*data_tree);
  osr_params.init(*data_tree);
  reset_current_external_function_options();

  lexer = make_unique<DynareFlex>(&in);
  lexer->set_debug(debug);

  Dynare::parser parser(*this);
  parser.set_debug_level(debug);
  parser.parse();

  return move(mod_file);
}

void
ParsingDriver::error(const Dynare::parser::location_type &l, const string &m)
{
  create_error_string(l, m, cerr);
  exit(EXIT_FAILURE);
}

void
ParsingDriver::error(const string &m)
{
  error(location, m);
}

void
ParsingDriver::create_error_string(const Dynare::parser::location_type &l, const string &m, ostream &stream)
{
  stream << "ERROR: " << *l.begin.filename << ": line " << l.begin.line;
  if (l.begin.line == l.end.line)
    if (l.begin.column == l.end.column - 1)
      stream << ", col " << l.begin.column;
    else
      stream << ", cols " << l.begin.column << "-" << l.end.column - 1;
  else
    stream << ", col " << l.begin.column << " -"
           << " line " << l.end.line << ", col " << l.end.column - 1;
  stream << ": " << m << endl;
}

void
ParsingDriver::create_error_string(const Dynare::parser::location_type &l, const string &m, const string &var)
{
  ostringstream stream;
  create_error_string(l, m, stream);
  model_errors.emplace_back(var, stream.str());
}

void
ParsingDriver::model_error(const string &m, const string &var)
{
  create_error_string(location, m, var);
}

void
ParsingDriver::undeclared_model_variable_error(const string &m, const string &var)
{
  ostringstream stream;
  if (!nostrict)
    {
      stream << "ERROR: " << *location.begin.filename << ": line " << location.begin.line;
      if (location.begin.line == location.end.line)
        if (location.begin.column == location.end.column - 1)
          stream << ", col " << location.begin.column;
        else
          stream << ", cols " << location.begin.column << "-" << location.end.column - 1;
      else
        stream << ", col " << location.begin.column << " -"
               << " line " << location.end.line << ", col " << location.end.column - 1;
      stream << ": ";
    }
  stream << m;
  if (nostrict)
    stream << " automatically declared exogenous.";
  undeclared_model_variable_errors.emplace_back(var, stream.str());
}

void
ParsingDriver::warning(const string &m)
{
  warnings << "WARNING: " << location << ": " << m << endl;
}

int
ParsingDriver::declare_symbol(const string &name, SymbolType type, const string &tex_name, const vector<pair<string, string>> &partition_value)
{
  int symb_id;
  try
    {
      symb_id = mod_file->symbol_table.addSymbol(name, type, tex_name, partition_value);
    }
  catch (SymbolTable::AlreadyDeclaredException &e)
    {
      if (e.same_type)
        {
          warning("Symbol " + name + " declared twice.");
          symb_id = mod_file->symbol_table.getID(name);
        }
      else
        error("Symbol " + name + " declared twice with different types!");
    }
  return symb_id;
}

int
ParsingDriver::declare_endogenous(const string &name, const string &tex_name, const vector<pair<string, string>> &partition_value)
{
  return declare_symbol(name, SymbolType::endogenous, tex_name, partition_value);
}

void
ParsingDriver::var(const vector<tuple<string, string, vector<pair<string, string>>>> &symbol_list,
                   bool log_option)
{
  for (auto &[name, tex_name, partition] : symbol_list)
    {
      int symb_id = declare_endogenous(name, tex_name, partition);
      if (log_option)
        mod_file->symbol_table.markWithLogTransform(symb_id);
    }
}

int
ParsingDriver::declare_exogenous(const string &name, const string &tex_name, const vector<pair<string, string>> &partition_value)
{
  return declare_symbol(name, SymbolType::exogenous, tex_name, partition_value);
}

void
ParsingDriver::varexo(const vector<tuple<string, string, vector<pair<string, string>>>> &symbol_list)
{
  for (auto &[name, tex_name, partition] : symbol_list)
    declare_exogenous(name, tex_name, partition);
}

void
ParsingDriver::varexo_det(const vector<tuple<string, string, vector<pair<string, string>>>> &symbol_list)
{
  for (auto &[name, tex_name, partition] : symbol_list)
    declare_symbol(name, SymbolType::exogenousDet, tex_name, partition);
}

int
ParsingDriver::declare_parameter(const string &name, const string &tex_name, const vector<pair<string, string>> &partition_value)
{
  return declare_symbol(name, SymbolType::parameter, tex_name, partition_value);
}

void
ParsingDriver::parameters(const vector<tuple<string, string, vector<pair<string, string>>>> &symbol_list)
{
  for (auto &[name, tex_name, partition] : symbol_list)
    declare_parameter(name, tex_name, partition);
}

void
ParsingDriver::declare_statement_local_variable(const string &name)
{
  if (mod_file->symbol_table.exists(name))
    error("Symbol " + name + " cannot be assigned within a statement "
          +"while being assigned elsewhere in the modfile");
  declare_symbol(name, SymbolType::statementDeclaredVariable, "", {});
}

void
ParsingDriver::set_planner_discount(expr_t value)
{
  planner_discount = value;
}

void
ParsingDriver::set_planner_discount_latex_name(string tex_name)
{
  planner_discount_latex_name = move(tex_name);
}

void
ParsingDriver::end_trend_var(bool log_trend, expr_t growth_factor, const vector<pair<string, string>> &symbol_list)
{
  /* Run detrending engine if trend variables are present, even if unused in
     a var(deflator=…) statement (see #113). */
  mod_file->nonstationary_variables = true;

  vector<int> declared_trend_vars;
  for (auto &[name, tex_name] : symbol_list)
    {
      int symb_id = declare_symbol(name, log_trend ? SymbolType::logTrend : SymbolType::trend, tex_name, {});
      declared_trend_vars.push_back(symb_id);
    }

  try
    {
      dynamic_model->addTrendVariables(declared_trend_vars, growth_factor);
    }
  catch (DataTree::TrendException &e)
    {
      error("Trend variable " + e.name + " was declared twice.");
    }
  declared_trend_vars.clear();
  reset_data_tree();
}

void
ParsingDriver::predetermined_variables(const vector<string> &symbol_list)
{
  for (auto &name : symbol_list)
    {
      check_symbol_is_endogenous(name);
      int symb_id = mod_file->symbol_table.getID(name);
      mod_file->symbol_table.markPredetermined(symb_id);
    }
}

void
ParsingDriver::add_equation_tags(string key, string value)
{
  if (eq_tags.contains(key))
    error("Tag '" + key + "' cannot be declared twice for the same equation");

  eq_tags[key] = value;

  transform(key.begin(), key.end(), key.begin(), ::tolower);
  if (key == "endogenous")
    declare_or_change_type(SymbolType::endogenous, value);
}

expr_t
ParsingDriver::add_non_negative_constant(const string &constant)
{
  return data_tree->AddNonNegativeConstant(constant);
}

expr_t
ParsingDriver::add_nan_constant()
{
  return data_tree->NaN;
}

expr_t
ParsingDriver::add_inf_constant()
{
  return data_tree->Infinity;
}

expr_t
ParsingDriver::add_model_variable(const string &name)
{
  if (name.find(".") != string::npos)
    error(name + " treated as a variable, but it contains a '.'");

  check_symbol_existence_in_model_block(name);
  int symb_id;
  try
    {
      symb_id = mod_file->symbol_table.getID(name);
      if (mod_file->symbol_table.getType(symb_id) == SymbolType::excludedVariable)
        error("Variable '" + name + "' can no longer be used since it has been excluded by a previous 'model_remove' or 'var_remove' statement");
    }
  catch (SymbolTable::UnknownSymbolNameException &e)
    {
      /* Declare variable as exogenous to continue parsing. Processing will end
         at end of model block (or planner_objective statement) if nostrict
         option was not passed. */
      symb_id = declare_exogenous(name);
      undeclared_model_vars.insert(name);
    }
  return add_model_variable(symb_id, 0);
}

expr_t
ParsingDriver::declare_or_change_type(SymbolType new_type, const string &name)
{
  int symb_id;
  try
    {
      symb_id = mod_file->symbol_table.getID(name);
      mod_file->symbol_table.changeType(symb_id, new_type);

      // remove error messages
      undeclared_model_vars.erase(name);
      erase_if(undeclared_model_variable_errors,
               [&name](auto &v) { return v.first == name; });
    }
  catch (SymbolTable::UnknownSymbolNameException &e)
    {
      switch (new_type)
        {
        case SymbolType::endogenous:
          symb_id = declare_endogenous(name);
          break;
        case SymbolType::exogenous:
          symb_id = declare_exogenous(name);
          break;
        case SymbolType::parameter:
          symb_id = declare_parameter(name);
          break;
        default:
          error("Type not yet supported");
        }
    }
  return add_model_variable(symb_id, 0);

}

expr_t
ParsingDriver::add_model_variable(int symb_id, int lag)
{
  assert(symb_id >= 0);
  SymbolType type = mod_file->symbol_table.getType(symb_id);

  if (type == SymbolType::modFileLocalVariable)
    error("Variable " + mod_file->symbol_table.getName(symb_id)
          +" not allowed inside model declaration. Its scope is only outside model.");

  if (type == SymbolType::externalFunction)
    error("Symbol " + mod_file->symbol_table.getName(symb_id)
          +" is a function name external to Dynare. It cannot be used like a variable without input argument inside model.");

  // See dynare#1765
  if (type == SymbolType::exogenousDet && lag != 0)
    error("Exogenous deterministic variable " + mod_file->symbol_table.getName(symb_id) + " cannot be given a lead or a lag.");

  if (type == SymbolType::modelLocalVariable && lag != 0)
    error("Model local variable " + mod_file->symbol_table.getName(symb_id) + " cannot be given a lead or a lag.");

  if (data_tree == planner_objective.get())
    {
      if (lag != 0)
        error("Leads and lags on variables are forbidden in 'planner_objective'.");

      if (type == SymbolType::modelLocalVariable)
        error("Model local variable " + mod_file->symbol_table.getName(symb_id) + " cannot be used in 'planner_objective'.");
    }

  if (data_tree == occbin_constraints_tree.get())
    {
      if (lag != 0)
        error("Leads and lags on variables are forbidden in 'occbin_constraints'. Note that you can achieve the same effect by introducing an auxiliary variable in the model.");

      if (type == SymbolType::modelLocalVariable)
        error("Model local variable " + mod_file->symbol_table.getName(symb_id) + " cannot be used in 'occbin_constraints'.");

      if (type == SymbolType::exogenous || type == SymbolType::exogenousDet)
        error("Exogenous variable " + mod_file->symbol_table.getName(symb_id) + " cannot be used in 'occbin_constraints'.");
    }

  // It makes sense to allow a lead/lag on parameters: during steady state calibration, endogenous and parameters can be swapped
  // NB: we use data_tree here, to avoid a crash in the occbin_constraints case
  return data_tree->AddVariable(symb_id, lag);
}

expr_t
ParsingDriver::add_expression_variable(const string &name)
{
  if (name.find(".") != string::npos)
    error(name + " treated as a variable, but it contains a '.'");

  if (parsing_epilogue && !mod_file->symbol_table.exists(name))
    error("Variable " + name + " used in the epilogue block but was not declared.");

  // If symbol doesn't exist, then declare it as a mod file local variable
  if (!mod_file->symbol_table.exists(name))
    mod_file->symbol_table.addSymbol(name, SymbolType::modFileLocalVariable);

  // This check must come after the previous one!
  if (mod_file->symbol_table.getType(name) == SymbolType::modelLocalVariable)
    error("Variable " + name + " not allowed outside model declaration. Its scope is only inside model.");

  if (mod_file->symbol_table.getType(name) == SymbolType::trend
      || mod_file->symbol_table.getType(name) == SymbolType::logTrend)
    error("Variable " + name + " not allowed outside model declaration, because it is a trend variable.");

  if (mod_file->symbol_table.getType(name) == SymbolType::externalFunction)
    error("Symbol '" + name + "' is the name of a MATLAB/Octave function, and cannot be used as a variable.");

  int symb_id = mod_file->symbol_table.getID(name);
  return data_tree->AddVariable(symb_id);
}

void
ParsingDriver::end_nonstationary_var(bool log_deflator, expr_t deflator, const vector<tuple<string, string, vector<pair<string, string>>>> &symbol_list, bool log_option)
{
  mod_file->nonstationary_variables = true;

  vector<int> declared_nonstationary_vars;
  for (auto &[name, tex_name, partition] : symbol_list)
    {
      int symb_id = declare_endogenous(name, tex_name, partition);
      declared_nonstationary_vars.push_back(symb_id);
      if (log_option)
        mod_file->symbol_table.markWithLogTransform(symb_id);
    }

  try
    {
      dynamic_model->addNonstationaryVariables(declared_nonstationary_vars, log_deflator, deflator);
    }
  catch (DataTree::TrendException &e)
    {
      error("Variable " + e.name + " was listed more than once as following a trend.");
    }

  set<int> r;
  deflator->collectVariables(SymbolType::endogenous, r);
  for (int it : r)
    if (dynamic_model->isNonstationary(it))
      error("The deflator contains a non-stationary endogenous variable. This is not allowed. Please use only stationary endogenous and/or {log_}trend_vars.");

  declared_nonstationary_vars.clear();
  reset_data_tree();
}

void
ParsingDriver::dsample(const string &arg1)
{
  int arg1_val = stoi(arg1);
  mod_file->addStatement(make_unique<DsampleStatement>(arg1_val));
}

void
ParsingDriver::dsample(const string &arg1, const string &arg2)
{
  int arg1_val = stoi(arg1);
  int arg2_val = stoi(arg2);
  mod_file->addStatement(make_unique<DsampleStatement>(arg1_val, arg2_val));
}

void
ParsingDriver::init_param(const string &name, expr_t rhs)
{
  check_symbol_is_parameter(name);
  int symb_id = mod_file->symbol_table.getID(name);
  mod_file->addStatement(make_unique<InitParamStatement>(symb_id, rhs, mod_file->symbol_table));
}

void
ParsingDriver::init_val(const string &name, expr_t rhs)
{
  if (nostrict && !mod_file->symbol_table.exists(name))
    {
      warning("discarding '" + name + "' as it was not recognized in the initval statement");
      return;
    }

  check_symbol_is_endogenous_or_exogenous(name, true);
  int symb_id = mod_file->symbol_table.getID(name);
  init_values.emplace_back(symb_id, rhs);
}

void
ParsingDriver::initval_file()
{
  mod_file->addStatement(make_unique<InitvalFileStatement>(move(options_list)));
  options_list.clear(); 
}

void
ParsingDriver::end_val(EndValLearntInStatement::LearntEndValType type, const string &name, expr_t rhs)
{
  if (nostrict && !mod_file->symbol_table.exists(name))
    {
      warning("discarding '" + name + "' as it was not recognized in the endval statement");
      return;
    }

  check_symbol_is_endogenous_or_exogenous(name, false);
  int symb_id = mod_file->symbol_table.getID(name);
  end_values.emplace_back(type, symb_id, rhs);
}

void
ParsingDriver::hist_val(const string &name, const string &lag, expr_t rhs)
{
  if (nostrict && !mod_file->symbol_table.exists(name))
    {
      warning("discarding '" + name + "' as it was not recognized in the histval block");
      return;
    }

  check_symbol_is_endogenous_or_exogenous(name, true);
  int symb_id = mod_file->symbol_table.getID(name);

  int ilag = stoi(lag);
  if (ilag > 0)
    error("histval: the lag on " + name + " should be less than or equal to 0");

  pair key{symb_id, ilag};

  if (hist_values.contains(key))
    error("hist_val: (" + name + ", " + lag + ") declared twice");

  hist_values[move(key)] = rhs;
}

void
ParsingDriver::homotopy_val(const string &name, expr_t val1, expr_t val2)
{
  check_symbol_existence(name);
  int symb_id = mod_file->symbol_table.getID(name);
  SymbolType type = mod_file->symbol_table.getType(symb_id);

  if (type != SymbolType::parameter
      && type != SymbolType::exogenous
      && type != SymbolType::exogenousDet)
    error("homotopy_val: " + name + " should be a parameter or exogenous variable");

  homotopy_values.emplace_back(symb_id, val1, val2);
}

void
ParsingDriver::end_generate_irfs()
{
  mod_file->addStatement(make_unique<GenerateIRFsStatement>(move(options_list),
                                                            move(generate_irf_names),
                                                            move(generate_irf_elements)));

  generate_irf_elements.clear();
  generate_irf_names.clear();
  options_list.clear();
}

void
ParsingDriver::add_generate_irfs_element(string name)
{
  for (const auto &it : generate_irf_names)
    if (it == name)
      error("Names in the generate_irfs block must be unique but you entered '"
            + name + "' more than once.");

  generate_irf_names.push_back(move(name));
  generate_irf_elements.push_back(generate_irf_exos);

  generate_irf_exos.clear();
}

void
ParsingDriver::add_generate_irfs_exog_element(string exo, const string &value)
{
  check_symbol_is_exogenous(exo, false);
  if (generate_irf_exos.contains(exo))
    error("You have set the exogenous variable " + exo + " twice.");

  generate_irf_exos[move(exo)] = stod(value);
}

void
ParsingDriver::forecast(vector<string> symbol_list)
{
  mod_file->addStatement(make_unique<ForecastStatement>(move(symbol_list), move(options_list),
                                                        mod_file->symbol_table));
  options_list.clear();
}

void
ParsingDriver::use_dll()
{
  mod_file->use_dll = true;
}

void
ParsingDriver::block()
{
  mod_file->block = true;
}

void
ParsingDriver::no_static()
{
  mod_file->no_static = true;
}

void
ParsingDriver::bytecode()
{
  mod_file->bytecode = true;
}

void
ParsingDriver::differentiate_forward_vars_all()
{
  mod_file->differentiate_forward_vars = true;
}

void
ParsingDriver::differentiate_forward_vars_some(vector<string> symbol_list)
{
  mod_file->differentiate_forward_vars = true;
  mod_file->differentiate_forward_vars_subset = move(symbol_list);
  for (auto &it : mod_file->differentiate_forward_vars_subset)
    check_symbol_is_endogenous(it);
}

void
ParsingDriver::cutoff(const string &value)
{
  double val = stod(value);
  mod_file->dynamic_model.cutoff = val;
}

void
ParsingDriver::mfs(const string &value)
{
  int val = stoi(value);
  mod_file->dynamic_model.mfs = val;
}

void
ParsingDriver::compilation_setup_substitute_flags(const string &flags)
{
  mod_file->dynamic_model.user_set_subst_flags = flags;
}

void
ParsingDriver::compilation_setup_add_flags(const string &flags)
{
  mod_file->dynamic_model.user_set_add_flags = flags;
}

void
ParsingDriver::compilation_setup_substitute_libs(const string &libs)
{
  mod_file->dynamic_model.user_set_subst_libs = libs;
}

void
ParsingDriver::compilation_setup_add_libs(const string &libs)
{
  mod_file->dynamic_model.user_set_add_libs = libs;
}

void
ParsingDriver::compilation_setup_compiler(const string &compiler)
{
  mod_file->dynamic_model.user_set_compiler = compiler;
}

void
ParsingDriver::balanced_growth_test_tol(const string &value)
{
  mod_file->dynamic_model.balanced_growth_test_tol = stod(value);
}

void
ParsingDriver::end_initval(bool all_values_required)
{
  mod_file->addStatement(make_unique<InitValStatement>(move(init_values), mod_file->symbol_table,
                                                       all_values_required));
  init_values.clear();
}

void
ParsingDriver::end_endval(bool all_values_required)
{
  InitOrEndValStatement::init_values_t end_values_new;
  for (auto [type, symb_id, value] : end_values)
    switch (type)
      {
      case EndValLearntInStatement::LearntEndValType::level:
        end_values_new.emplace_back(symb_id, value);
        break;
      case EndValLearntInStatement::LearntEndValType::add:
        error("endval: '" + mod_file->symbol_table.getName(symb_id) + " += ...' line not allowed unless 'learnt_in' option with value >1 is passed");
      case EndValLearntInStatement::LearntEndValType::multiply:
        error("endval: '" + mod_file->symbol_table.getName(symb_id) + " *= ...' line not allowed unless 'learnt_in' option with value >1 is passed");
      }

  mod_file->addStatement(make_unique<EndValStatement>(move(end_values_new), mod_file->symbol_table,
                                                      all_values_required));
  end_values.clear();
}

void
ParsingDriver::end_endval_learnt_in(const string &learnt_in_period)
{
  int learnt_in_period_int = stoi(learnt_in_period);
  if (learnt_in_period_int < 1)
    error("endval: value '" + learnt_in_period + "' is not allowed for 'learnt_in' option");
  if (learnt_in_period_int == 1)
    {
      end_endval(false);
      return;
    }
  for (auto [type, symb_id, value] : end_values)
    if (mod_file->symbol_table.getType(symb_id) != SymbolType::exogenous)
      error("endval(learnt_in=...): " + mod_file->symbol_table.getName(symb_id) + " is not an exogenous variable");
  mod_file->addStatement(make_unique<EndValLearntInStatement>(learnt_in_period_int, move(end_values),
                                                              mod_file->symbol_table));
  end_values.clear();
}

void
ParsingDriver::end_histval(bool all_values_required)
{
  mod_file->addStatement(make_unique<HistValStatement>(move(hist_values), mod_file->symbol_table,
                                                       all_values_required));
  hist_values.clear();
}

void
ParsingDriver::end_homotopy()
{
  mod_file->addStatement(make_unique<HomotopyStatement>(move(homotopy_values), mod_file->symbol_table));
  homotopy_values.clear();
}

void
ParsingDriver::begin_epilogue()
{
  parsing_epilogue = true;
  set_current_data_tree(&mod_file->epilogue);
}

void
ParsingDriver::end_epilogue()
{
  parsing_epilogue = false;
  reset_data_tree();
}

void
ParsingDriver::add_epilogue_variable(const string &name)
{
  declare_symbol(name, SymbolType::epilogue, "", {});
}

void
ParsingDriver::add_epilogue_equal(const string &name, expr_t expr)
{
  mod_file->epilogue.addDefinition(mod_file->symbol_table.getID(name), expr);
}

void
ParsingDriver::begin_model()
{
  set_current_data_tree(&mod_file->dynamic_model);
}

void
ParsingDriver::end_model()
{
  bool exit_after_write = false;
  if (model_errors.size() > 0)
    for (auto &it : model_errors)
      {
        if (it.first.empty())
          exit_after_write = true;
        cerr << it.second;
      }

  if (undeclared_model_variable_errors.size() > 0)
    for (auto &it : undeclared_model_variable_errors)
      {
        if (nostrict)
          warning(it.second);
        else
          {
            exit_after_write = true;
            cerr << it.second << endl;
          }
      }
  undeclared_model_variable_errors.clear();

  if (exit_after_write)
    exit(EXIT_FAILURE);

  reset_data_tree();
}

void
ParsingDriver::end_shocks(bool overwrite)
{
  mod_file->addStatement(make_unique<ShocksStatement>(overwrite, move(det_shocks), move(var_shocks),
                                                      move(std_shocks), move(covar_shocks),
                                                      move(corr_shocks), mod_file->symbol_table));
  det_shocks.clear();
  if (!learnt_shocks_add.empty())
    error("shocks: 'add' keyword not allowed unless 'learnt_in' option with value >1 is passed");
  if (!learnt_shocks_multiply.empty())
    error("shocks: 'multiply' keyword not allowed unless 'learnt_in' option with value >1 is passed");
  var_shocks.clear();
  std_shocks.clear();
  covar_shocks.clear();
  corr_shocks.clear();
}

void
ParsingDriver::end_mshocks(bool overwrite)
{
  mod_file->addStatement(make_unique<MShocksStatement>(overwrite, move(det_shocks),
                                                       mod_file->symbol_table));
  det_shocks.clear();
  if (!learnt_shocks_add.empty())
    error("mshocks: 'add' keyword not allowed");
  if (!learnt_shocks_multiply.empty())
    error("mshocks: 'multiply' keyword not allowed");
}

void
ParsingDriver::end_shocks_surprise(bool overwrite)
{
  mod_file->addStatement(make_unique<ShocksSurpriseStatement>(overwrite, move(det_shocks),
                                                              mod_file->symbol_table));
  det_shocks.clear();
  if (!learnt_shocks_add.empty())
    error("shocks(surprise): 'add' keyword not allowed");
  if (!learnt_shocks_multiply.empty())
    error("shocks(surprise): 'multiply' keyword not allowed");
}

void
ParsingDriver::end_shocks_learnt_in(const string &learnt_in_period, bool overwrite)
{
  int learnt_in_period_int = stoi(learnt_in_period);
  if (learnt_in_period_int < 1)
    error("shocks: value '" + learnt_in_period + "' is not allowed for 'learnt_in' option");
  if (learnt_in_period_int == 1)
    {
      end_shocks(overwrite);
      return;
    }
  for (auto &[symb_id, vals] : det_shocks)
    for (auto [period1, period2, expr] : vals)
      if (period1 < learnt_in_period_int)
        error("shocks: for variable " + mod_file->symbol_table.getName(symb_id) + ", shock period (" + to_string(period1) + ") is earlier than the period in which the shock is learnt (" + learnt_in_period + ")");

  // Aggregate the three types of shocks
  ShocksLearntInStatement::learnt_shocks_t learnt_shocks;
  for (const auto &[id, v] : det_shocks)
    {
      vector<tuple<ShocksLearntInStatement::LearntShockType, int, int, expr_t>> v2;
      for (auto [period1, period2, value] : v)
        v2.emplace_back(ShocksLearntInStatement::LearntShockType::level, period1, period2, value);
      learnt_shocks[id] = v2;
    }
  for (const auto &[id, v] : learnt_shocks_add)
    {
      vector<tuple<ShocksLearntInStatement::LearntShockType, int, int, expr_t>> v2;
      for (auto [period1, period2, value] : v)
        v2.emplace_back(ShocksLearntInStatement::LearntShockType::add, period1, period2, value);
      learnt_shocks[id] = v2;
    }
  for (const auto &[id, v] : learnt_shocks_multiply)
    {
      vector<tuple<ShocksLearntInStatement::LearntShockType, int, int, expr_t>> v2;
      for (auto [period1, period2, value] : v)
        v2.emplace_back(ShocksLearntInStatement::LearntShockType::multiply, period1, period2, value);
      learnt_shocks[id] = v2;
    }

  mod_file->addStatement(make_unique<ShocksLearntInStatement>(learnt_in_period_int, overwrite,
                                                              move(learnt_shocks),
                                                              mod_file->symbol_table));
  det_shocks.clear();
  learnt_shocks_add.clear();
  learnt_shocks_multiply.clear();
}

void
ParsingDriver::end_heteroskedastic_shocks(bool overwrite)
{
  mod_file->addStatement(make_unique<HeteroskedasticShocksStatement>(overwrite,
                                                                     move(heteroskedastic_shocks_values),
                                                                     move(heteroskedastic_shocks_scales),
                                                                     mod_file->symbol_table));
  heteroskedastic_shocks_values.clear();
  heteroskedastic_shocks_scales.clear();
}

void
ParsingDriver::add_det_shock(const string &var, const vector<pair<int, int>> &periods, const vector<expr_t> &values, DetShockType type)
{
  switch (type)
    {
    case DetShockType::conditional_forecast:
      check_symbol_is_endogenous(var);
      break;
    case DetShockType::standard:
       // Allow exo_det, for stochastic context
      check_symbol_is_exogenous(var, true);
      break;
    case DetShockType::add:
    case DetShockType::multiply:
      check_symbol_is_exogenous(var, false);
      break;
    }

  int symb_id = mod_file->symbol_table.getID(var);

  if (det_shocks.contains(symb_id) || learnt_shocks_add.contains(symb_id)
      || learnt_shocks_multiply.contains(symb_id))
    error("shocks/conditional_forecast_paths: variable " + var + " declared twice");

  if (periods.size() != values.size())
    error("shocks/conditional_forecast_paths: variable " + var + ": number of periods is different from number of shock values");

  vector<tuple<int, int, expr_t>> v;

  for (size_t i = 0; i < periods.size(); i++)
    v.emplace_back(periods[i].first, periods[i].second, values[i]);

  switch (type)
    {
    case DetShockType::standard:
    case DetShockType::conditional_forecast:
      det_shocks[symb_id] = v;
      break;
    case DetShockType::add:
      learnt_shocks_add[symb_id] = v;
      break;
    case DetShockType::multiply:
      learnt_shocks_multiply[symb_id] = v;
      break;
    }
}

void
ParsingDriver::add_heteroskedastic_shock(const string &var, const vector<pair<int, int>> &periods, const vector<expr_t> &values, bool scales)
{
  check_symbol_is_exogenous(var, false);

  int symb_id = mod_file->symbol_table.getID(var);

  if ((!scales && heteroskedastic_shocks_values.contains(symb_id))
      || (scales && heteroskedastic_shocks_scales.contains(symb_id)))
    error("heteroskedastic_shocks: variable " + var + " declared twice");

  if (periods.size() != values.size())
    error("heteroskedastic_shocks: variable " + var + ": number of periods is different from number of shock values");

  vector<tuple<int, int, expr_t>> v;
  for (size_t i = 0; i < periods.size(); i++)
    v.emplace_back(periods[i].first, periods[i].second, values[i]);

  if (scales)
    heteroskedastic_shocks_scales[symb_id] = v;
  else
    heteroskedastic_shocks_values[symb_id] = v;
}

void
ParsingDriver::add_stderr_shock(const string &var, expr_t value)
{
  if (nostrict && !mod_file->symbol_table.exists(var))
    {
      warning("discarding shocks block declaration of the standard error of '" + var + "' as it was not declared");
      return;
    }

  check_symbol_existence(var);
  int symb_id = mod_file->symbol_table.getID(var);

  if (var_shocks.contains(symb_id) || std_shocks.contains(symb_id))
    error("shocks: variance or stderr of shock on " + var + " declared twice");

  std_shocks[symb_id] = value;
}

void
ParsingDriver::add_var_shock(const string &var, expr_t value)
{
  if (nostrict && !mod_file->symbol_table.exists(var))
    {
      warning("discarding shocks block declaration of the variance of '" + var + "' as it was not declared");
      return;
    }

  check_symbol_existence(var);
  int symb_id = mod_file->symbol_table.getID(var);

  if (var_shocks.contains(symb_id) || std_shocks.contains(symb_id))
    error("shocks: variance or stderr of shock on " + var + " declared twice");

  var_shocks[symb_id] = value;
}

void
ParsingDriver::add_covar_shock(const string &var1, const string &var2, expr_t value)
{
  if (nostrict &&
      (!mod_file->symbol_table.exists(var1) || !mod_file->symbol_table.exists(var2)))
    {
      warning("discarding shocks block declaration of the covariance of '" + var1 + "' and '" + var2 + "' as at least one was not declared");
      return;
    }

  check_symbol_existence(var1);
  check_symbol_existence(var2);
  int symb_id1 = mod_file->symbol_table.getID(var1);
  int symb_id2 = mod_file->symbol_table.getID(var2);

  pair key{symb_id1, symb_id2}, key_inv{symb_id2, symb_id1};

  if (covar_shocks.contains(key) || covar_shocks.contains(key_inv)
      || corr_shocks.contains(key) || corr_shocks.contains(key_inv))
    error("shocks: covariance or correlation shock on variable pair (" + var1 + ", "
          + var2 + ") declared twice");

  covar_shocks[key] = value;
}

void
ParsingDriver::add_correl_shock(const string &var1, const string &var2, expr_t value)
{
  if (nostrict &&
      (!mod_file->symbol_table.exists(var1) || !mod_file->symbol_table.exists(var2)))
    {
      warning("discarding shocks block declaration of the correlation of '" + var1 + "' and '" + var2 + "' as at least one was not declared");
      return;
    }

  check_symbol_existence(var1);
  check_symbol_existence(var2);
  int symb_id1 = mod_file->symbol_table.getID(var1);
  int symb_id2 = mod_file->symbol_table.getID(var2);

  pair key{symb_id1, symb_id2}, key_inv{symb_id2, symb_id1};

  if (covar_shocks.contains(key) || covar_shocks.contains(key_inv)
      || corr_shocks.contains(key) || corr_shocks.contains(key_inv))
    error("shocks: covariance or correlation shock on variable pair (" + var1 + ", "
          + var2 + ") declared twice");

  corr_shocks[key] = value;
}

void
ParsingDriver::begin_svar_identification()
{
  svar_upper_cholesky = false;
  svar_lower_cholesky = false;
  svar_constants_exclusion = false;
}

void
ParsingDriver::end_svar_identification()
{
  mod_file->addStatement(make_unique<SvarIdentificationStatement>(move(svar_ident_restrictions),
                                                                  svar_upper_cholesky,
                                                                  svar_lower_cholesky,
                                                                  svar_constants_exclusion,
                                                                  mod_file->symbol_table));
  svar_equation_restrictions.clear();
  svar_ident_restrictions.clear();
  svar_Qi_restriction_nbr.clear();
  svar_Ri_restriction_nbr.clear();
}

void
ParsingDriver::combine_lag_and_restriction(const string &lag)
{
  int current_lag = stoi(lag);

  for (const auto &it : svar_ident_restrictions)
    if (it.lag == current_lag)
      error("lag " + lag + " used more than once.");

  for (const auto &it : svar_equation_restrictions)
    for (auto it1 : it.second)
      {
        SvarIdentificationStatement::svar_identification_restriction new_restriction;
        new_restriction.equation = it.first;
        if (current_lag > 0)
          new_restriction.restriction_nbr = ++svar_Ri_restriction_nbr[it.first];
        else
          new_restriction.restriction_nbr = ++svar_Qi_restriction_nbr[it.first];
        new_restriction.lag = current_lag;
        new_restriction.variable = it1;
        new_restriction.value = data_tree->One;
        svar_ident_restrictions.push_back(new_restriction);
      }

  svar_upper_cholesky = false;
  svar_lower_cholesky = false;
  svar_equation_restrictions.clear();
}

void
ParsingDriver::add_restriction_in_equation(const string &equation, const vector<string> &symbol_list)
{
  int eqn = stoi(equation);
  if (eqn < 1)
    error("equation numbers must be greater than or equal to 1.");

  if (svar_equation_restrictions.contains(eqn))
    error("equation number " + equation + " referenced more than once under a single lag.");

  vector<int> svar_restriction_symbols;
  for (auto &name : symbol_list)
    {
      check_symbol_existence(name);
      int symb_id = mod_file->symbol_table.getID(name);

      for (const auto &viit : svar_restriction_symbols)
        if (symb_id == viit)
          error(name + " restriction added twice.");

      svar_restriction_symbols.push_back(symb_id);
    }
  svar_equation_restrictions[eqn] = svar_restriction_symbols;
}

void
ParsingDriver::add_restriction_equation_nbr(const string &eq_nbr)
{
  svar_equation_nbr = stoi(eq_nbr);
  svar_left_handside = true;
  // reinitialize restriction type that must be set from the first restriction element
  svar_restriction_type = SvarRestrictionType::NOT_SET;
}

void
ParsingDriver::add_restriction_equal()
{
  if (svar_left_handside)
    svar_left_handside = false;
  else
    error("svar_identification: there are more than one EQUAL sign in a restriction equation");
}

void
ParsingDriver::add_positive_restriction_element(expr_t value, const string &variable, const string &lag)
{
  // if the expression is not on the left handside, change its sign
  if (!svar_left_handside)
    value = add_uminus(value);

  add_restriction_element(value, variable, lag);
}

void
ParsingDriver::add_positive_restriction_element(const string &variable, const string &lag)
{
  expr_t value(data_tree->One);

  // if the expression is not on the left handside, change its sign
  if (!svar_left_handside)
    value = add_uminus(value);

  add_restriction_element(value, variable, lag);
}

void
ParsingDriver::add_negative_restriction_element(expr_t value, const string &variable, const string &lag)
{
  // if the expression is on the left handside, change its sign
  if (svar_left_handside)
    value = add_uminus(value);

  add_restriction_element(value, variable, lag);
}

void
ParsingDriver::add_negative_restriction_element(const string &variable, const string &lag)
{
  expr_t value(data_tree->One);

  // if the expression is on the left handside, change its sign
  if (svar_left_handside)
    value = add_uminus(value);

  add_restriction_element(value, variable, lag);
}

void
ParsingDriver::add_restriction_element(expr_t value, const string &variable, const string &lag)
{
  check_symbol_existence(variable);
  int symb_id = mod_file->symbol_table.getID(variable);

  int current_lag = stoi(lag);
  if (svar_restriction_type == SvarRestrictionType::NOT_SET)
    {
      if (current_lag == 0)
        {
          svar_restriction_type = SvarRestrictionType::Qi_TYPE;
          ++svar_Qi_restriction_nbr[svar_equation_nbr];
        }
      else
        {
          svar_restriction_type = SvarRestrictionType::Ri_TYPE;
          ++svar_Ri_restriction_nbr[svar_equation_nbr];
        }
    }
  else
    {
      if ((svar_restriction_type == SvarRestrictionType::Qi_TYPE && current_lag > 0)
          || (svar_restriction_type == SvarRestrictionType::Ri_TYPE && current_lag == 0))
        error("SVAR_IDENTIFICATION: a single restrictions must affect either Qi or Ri, but not both");
    }
  SvarIdentificationStatement::svar_identification_restriction new_restriction;
  new_restriction.equation = svar_equation_nbr;
  if (current_lag > 0)
    new_restriction.restriction_nbr = svar_Ri_restriction_nbr[svar_equation_nbr];
  else
    new_restriction.restriction_nbr = svar_Qi_restriction_nbr[svar_equation_nbr];
  new_restriction.lag = current_lag;
  new_restriction.variable = symb_id;
  new_restriction.value = value;

  svar_ident_restrictions.push_back(new_restriction);
}

void
ParsingDriver::check_restriction_expression_constant(expr_t value)
{
  if (value->eval({}) != 0)
    error("SVAR_INDENTIFICATION restrictions must be homogenous");
}

void
ParsingDriver::add_upper_cholesky()
{
  svar_upper_cholesky = true;
}

void
ParsingDriver::add_lower_cholesky()
{
  svar_lower_cholesky = true;
}

void
ParsingDriver::add_constants_exclusion()
{
  svar_constants_exclusion = true;
}

void
ParsingDriver::add_svar_global_identification_check()
{
  mod_file->addStatement(make_unique<SvarGlobalIdentificationCheckStatement>());
}

void
ParsingDriver::steady()
{
  mod_file->addStatement(make_unique<SteadyStatement>(move(options_list)));
  options_list.clear();
}

void
ParsingDriver::option_num(string name_option, string opt1, string opt2)
{
  if (options_list.contains(name_option))
    error("option " + name_option + " declared twice");

  options_list.set(move(name_option), pair{move(opt1), move(opt2)});
}

void
ParsingDriver::option_num(string name_option, string opt)
{
  if (options_list.contains(name_option))
    error("option " + name_option + " declared twice");

  options_list.set(move(name_option), OptionsList::NumVal{move(opt)});
}

void
ParsingDriver::option_str(string name_option, string opt)
{
  if (options_list.contains(name_option))
    error("option " + name_option + " declared twice");

  options_list.set(move(name_option), OptionsList::StringVal{move(opt)});
}

void
ParsingDriver::option_date(string name_option, string opt)
{
  if (options_list.contains(name_option))
    error("option " + name_option + " declared twice");

  options_list.set(move(name_option), OptionsList::DateVal{move(opt)});
}

void
ParsingDriver::option_symbol_list(string name_option, vector<string> symbol_list)
{
  if (options_list.contains(name_option))
    error("option " + name_option + " declared twice");

  if (name_option == "irf_shocks")
    for (auto &shock : symbol_list)
      {
        if (!mod_file->symbol_table.exists(shock))
          error("Unknown symbol: " + shock);
        if (mod_file->symbol_table.getType(shock) != SymbolType::exogenous)
          error("Variables passed to irf_shocks must be exogenous. Caused by: " + shock);
      }

  if (name_option == "ms.parameters")
    for (auto &it : symbol_list)
      if (mod_file->symbol_table.getType(it) != SymbolType::parameter)
        error("Variables passed to the parameters option of the markov_switching statement must be parameters. Caused by: " + it);

  options_list.set(move(name_option), OptionsList::SymbolListVal{move(symbol_list)});
}

void
ParsingDriver::option_vec_int(string name_option, vector<int> opt)
{
  if (options_list.contains(name_option))
    error("option " + name_option + " declared twice");

  if (opt.empty())
    error("option " + name_option + " was passed an empty vector.");

  options_list.set(move(name_option), move(opt));
}

void
ParsingDriver::option_vec_str(string name_option, vector<string> opt)
{
  if (options_list.contains(name_option))
    error("option " + name_option + " declared twice");

  if (opt.empty())
    error("option " + name_option + " was passed an empty vector.");

  options_list.set(move(name_option), OptionsList::VecStrVal{move(opt)});
}

void
ParsingDriver::option_vec_cellstr(string name_option, vector<string> opt)
{
  if (options_list.contains(name_option))
    error("option " + name_option + " declared twice");

  if (opt.empty())
    error("option " + name_option + " was passed an empty vector.");

  options_list.set(move(name_option), OptionsList::VecCellStrVal{move(opt)});
}

void
ParsingDriver::option_vec_value(string name_option, vector<string> opt)
{
  if (options_list.contains(name_option))
    error("option " + name_option + " declared twice");

  if (opt.empty())
    error("option " + name_option + " was passed an empty vector.");

  options_list.set(move(name_option), OptionsList::VecValueVal{move(opt)});
}

void
ParsingDriver::option_vec_of_vec_value(string name_option, vector<vector<string>> opt)
{
  if (options_list.contains(name_option))
    error("option " + name_option + " declared twice");

  if (opt.empty())
    error("option " + name_option + " was passed an empty vector.");

  options_list.set(move(name_option), move(opt));
}

void
ParsingDriver::linear()
{
  mod_file->linear = true;
}

void
ParsingDriver::rplot(vector<string> symbol_list)
{
  mod_file->addStatement(make_unique<RplotStatement>(move(symbol_list), mod_file->symbol_table));
}

void
ParsingDriver::stoch_simul(SymbolList symbol_list)
{
  //make sure default order is known to preprocessor, see #49
  if (!options_list.contains("order"))
    options_list.set("order", OptionsList::NumVal{"2"});

  symbol_list.removeDuplicates("stoch_simul", warnings);

  mod_file->addStatement(make_unique<StochSimulStatement>(move(symbol_list), move(options_list),
                                                          mod_file->symbol_table));
  options_list.clear();
}

void
ParsingDriver::trend_component_model()
{
  try
    {
      mod_file->trend_component_model_table.addTrendComponentModel(options_list.get<OptionsList::StringVal>("trend_component.name"),
                                                                   options_list.get<OptionsList::VecStrVal>("trend_component.eqtags"),
                                                                   options_list.get<OptionsList::VecStrVal>("trend_component.targets"));
    }
  catch (OptionsList::UnknownOptionException &e)
    {
      string name {e.name.substr(16)};
      if (name == "name")
        name = "model_name";
      error("You must pass the '" + name + "' option to the 'trend_component_model' statement.");
    }

  options_list.clear();
}

void
ParsingDriver::var_model()
{
  try
    {
      mod_file->var_model_table.addVarModel(options_list.get<OptionsList::StringVal>("var.model_name"),
                                            options_list.get_if<OptionsList::NumVal>("var.structural").value_or(OptionsList::NumVal{"false"}) == "true",
                                            options_list.get<OptionsList::VecStrVal>("var.eqtags"));
    }
  catch (OptionsList::UnknownOptionException &e)
    {
      error("You must pass the '" + e.name.substr(4) + "' option to the 'var_model' statement.");
    }
  options_list.clear();
}

void
ParsingDriver::simul()
{
  warning("The 'simul' statement is deprecated. Please use 'perfect_foresight_setup' and 'perfect_foresight_solver' instead.");
  mod_file->addStatement(make_unique<SimulStatement>(move(options_list)));
  options_list.clear();
}

void
ParsingDriver::model_info()
{
  mod_file->addStatement(make_unique<ModelInfoStatement>(move(options_list)));
  options_list.clear();
}

void
ParsingDriver::check()
{
  mod_file->addStatement(make_unique<CheckStatement>(move(options_list)));
  options_list.clear();
}

void
ParsingDriver::add_estimated_params_element()
{
  if (estim_params.name != "dsge_prior_weight")
    {
      check_symbol_existence(estim_params.name);
      SymbolType type = mod_file->symbol_table.getType(estim_params.name);
      switch (estim_params.type)
        {
        case 1:
          if (type != SymbolType::endogenous && type != SymbolType::exogenous)
            error(estim_params.name + " must be an endogenous or an exogenous variable");
          break;
        case 2:
          check_symbol_is_parameter(estim_params.name);
          break;
        case 3:
          check_symbol_existence(estim_params.name2);
          SymbolType type2 = mod_file->symbol_table.getType(estim_params.name2);
          if ((type != SymbolType::endogenous && type != SymbolType::exogenous) || type != type2)
            error(estim_params.name + " and " + estim_params.name2 + " must either be both endogenous variables or both exogenous");
          break;
        }
    }
  estim_params_list.push_back(estim_params);
  estim_params.init(*data_tree);
}

void
ParsingDriver::estimated_params(bool overwrite)
{
  mod_file->addStatement(make_unique<EstimatedParamsStatement>(move(estim_params_list),
                                                               mod_file->symbol_table, overwrite));
  estim_params_list.clear();
}

void
ParsingDriver::estimated_params_init(bool use_calibration)
{
  mod_file->addStatement(make_unique<EstimatedParamsInitStatement>(move(estim_params_list),
                                                                   mod_file->symbol_table,
                                                                   use_calibration));
  estim_params_list.clear();
}

void
ParsingDriver::estimated_params_bounds()
{
  mod_file->addStatement(make_unique<EstimatedParamsBoundsStatement>(move(estim_params_list),
                                                                     mod_file->symbol_table));
  estim_params_list.clear();
}

void
ParsingDriver::estimated_params_remove()
{
  mod_file->addStatement(make_unique<EstimatedParamsRemoveStatement>(move(estim_params_list),
                                                                     mod_file->symbol_table));
  estim_params_list.clear();
}

void
ParsingDriver::add_osr_params_element()
{
  check_symbol_existence(osr_params.name);
  SymbolType type = mod_file->symbol_table.getType(osr_params.name);
  if (type != SymbolType::parameter)
    error(osr_params.name + " must be a parameter to be used in the osr_bounds block");
  osr_params_list.push_back(osr_params);
  osr_params.init(*data_tree);
}

void
ParsingDriver::osr_params_bounds()
{
  mod_file->addStatement(make_unique<OsrParamsBoundsStatement>(move(osr_params_list)));
  osr_params_list.clear();
}

void
ParsingDriver::set_unit_root_vars()
{
  mod_file->addStatement(make_unique<UnitRootVarsStatement>());
  warning("''unit_root_vars'' is now obsolete; use the ''diffuse_filter'' option of ''estimation'' instead");
}

void
ParsingDriver::set_time(const string &arg)
{
  option_date("initial_period", arg);
  mod_file->addStatement(make_unique<SetTimeStatement>(move(options_list)));
  options_list.clear();
}

void
ParsingDriver::estimation_data()
{
  mod_file->addStatement(make_unique<EstimationDataStatement>(move(options_list)));
  options_list.clear();
}

void
ParsingDriver::set_subsamples(string name1, string name2)
{
  check_symbol_existence(name1);
  if (!name2.empty())
    check_symbol_existence(name2);

  mod_file->addStatement(make_unique<SubsamplesStatement>(name1, name2, subsample_declaration_map,
                                                          mod_file->symbol_table));
  subsample_declarations[{ move(name1), move(name2) }] = move(subsample_declaration_map);
  subsample_declaration_map.clear();
}

void
ParsingDriver::copy_subsamples(string to_name1, string to_name2, string from_name1, string from_name2)
{
  check_symbol_existence(to_name1);
  check_symbol_existence(from_name1);
  if (!to_name2.empty())
    check_symbol_existence(to_name2);
  if (!from_name2.empty())
    check_symbol_existence(from_name2);

  if (!subsample_declarations.contains({ from_name1, from_name2 }))
    {
      string err{from_name1};
      if (!from_name2.empty())
        err.append(",").append(from_name2);
      error(err + " does not have an associated subsample statement.");
    }

  mod_file->addStatement(make_unique<SubsamplesEqualStatement>(to_name1, to_name2, from_name1, from_name2,
                                                               mod_file->symbol_table));

  subsample_declarations[{ move(to_name1), move(to_name2) }]
    = subsample_declarations[{ move(from_name1), move(from_name2) }];
}

void
ParsingDriver::check_symbol_is_statement_variable(const string &name)
{
  check_symbol_existence(name);
  int symb_id = mod_file->symbol_table.getID(name);
  if (mod_file->symbol_table.getType(symb_id) != SymbolType::statementDeclaredVariable)
    error(name + " is not a variable assigned in a statement");
}

void
ParsingDriver::set_subsample_name_equal_to_date_range(string name, string date1, string date2)
{
  if (subsample_declaration_map.contains(name))
    error("Symbol " + name + " may only be assigned once in a SUBSAMPLE statement");
  subsample_declaration_map[move(name)] = { move(date1), move(date2) };
}

void
ParsingDriver::check_subsample_declaration_exists(const string &name1, const string &subsample_name)
{
  if (subsample_name.empty())
    return;

  check_subsample_declaration_exists(name1, "", subsample_name);
}

void
ParsingDriver::check_subsample_declaration_exists(const string &name1, const string &name2, const string &subsample_name)
{
  if (subsample_name.empty())
    return;

  check_symbol_existence(name1);
  if (!name2.empty())
    check_symbol_existence(name2);

  auto it = subsample_declarations.find({ name1, name2 });
  if (it == subsample_declarations.end())
    {
      it = subsample_declarations.find({ name2, name1 });
      if (it == subsample_declarations.end())
        {
          string err{name1};
          if (!name2.empty())
            err.append(",").append(name2);
          error("A subsample statement has not been issued for " + err);
        }
    }

  auto tmp_map = it->second;
  if (!tmp_map.contains(subsample_name))
    error("The subsample name " + subsample_name + " was not previously declared in a subsample statement.");
}

void
ParsingDriver::set_prior(string name, string subsample_name)
{
  check_symbol_is_parameter(name);
  check_subsample_declaration_exists(name, subsample_name);
  mod_file->addStatement(make_unique<PriorStatement>(move(name), move(subsample_name), prior_shape,
                                                     prior_variance, move(options_list)));
  options_list.clear();
  set_prior_variance();
  prior_shape = PriorDistributions::noShape;
}

void
ParsingDriver::set_joint_prior(const vector<string> &symbol_vec)
{
  for (auto &it : symbol_vec)
    add_joint_parameter(it);
  mod_file->addStatement(make_unique<JointPriorStatement>(move(joint_parameters), prior_shape,
                                                          move(options_list)));
  joint_parameters.clear();
  options_list.clear();
  prior_shape = PriorDistributions::noShape;
}

void
ParsingDriver::add_joint_parameter(string name)
{
  check_symbol_is_parameter(name);
  joint_parameters.push_back(move(name));
}

void
ParsingDriver::set_prior_variance(expr_t variance)
{
  prior_variance = variance;
}

void
ParsingDriver::copy_prior(string to_declaration_type, string to_name1,
                          string to_name2, string to_subsample_name,
                          string from_declaration_type, string from_name1,
                          string from_name2, string from_subsample_name)
{
  if (to_declaration_type == "par")
    check_symbol_is_parameter(to_name1);
  else
    {
      check_symbol_is_endogenous_or_exogenous(to_name1, false);
      if (!to_name2.empty())
        check_symbol_is_endogenous_or_exogenous(to_name2, false);
    }

  if (from_declaration_type == "par")
    check_symbol_is_parameter(from_name1);
  else
    {
      check_symbol_is_endogenous_or_exogenous(from_name1, false);
      if (!from_name2.empty())
        check_symbol_is_endogenous_or_exogenous(from_name2, false);
    }

  mod_file->addStatement(make_unique<PriorEqualStatement>(move(to_declaration_type), move(to_name1),
                                                          move(to_name2), move(to_subsample_name),
                                                          move(from_declaration_type), move(from_name1),
                                                          move(from_name2), move(from_subsample_name),
                                                          mod_file->symbol_table));
}

void
ParsingDriver::set_options(string name, string subsample_name)
{
  check_symbol_is_parameter(name);
  check_subsample_declaration_exists(name, subsample_name);
  mod_file->addStatement(make_unique<OptionsStatement>(move(name), move(subsample_name),
                                                       move(options_list)));
  options_list.clear();
}

void
ParsingDriver::copy_options(string to_declaration_type, string to_name1,
                            string to_name2, string to_subsample_name,
                            string from_declaration_type, string from_name1,
                            string from_name2, string from_subsample_name)
{
  if (to_declaration_type == "par")
    check_symbol_is_parameter(to_name1);
  else
    {
      check_symbol_is_endogenous_or_exogenous(to_name1, false);
      if (!to_name2.empty())
        check_symbol_is_endogenous_or_exogenous(to_name2, false);
    }

  if (from_declaration_type == "par")
    check_symbol_is_parameter(from_name1);
  else
    {
      check_symbol_is_endogenous_or_exogenous(from_name1, false);
      if (!from_name2.empty())
        check_symbol_is_endogenous_or_exogenous(from_name2, false);
    }

  mod_file->addStatement(make_unique<OptionsEqualStatement>(move(to_declaration_type), move(to_name1),
                                                            move(to_name2), move(to_subsample_name),
                                                            move(from_declaration_type), move(from_name1),
                                                            move(from_name2), move(from_subsample_name),
                                                            mod_file->symbol_table));
}

void
ParsingDriver::check_symbol_is_endogenous_or_exogenous(const string &name, bool allow_det)
{
  check_symbol_existence(name);
  switch (mod_file->symbol_table.getType(name))
    {
    case SymbolType::endogenous:
    case SymbolType::exogenous:
      break;
    case SymbolType::exogenousDet:
      if (!allow_det)
        error(name + " is an exogenous deterministic.");
      break;
    default:
      error(name + " is neither endogenous or exogenous.");
    }
}

void
ParsingDriver::check_symbol_is_endogenous(const string &name)
{
  check_symbol_existence(name);
  if (mod_file->symbol_table.getType(name) != SymbolType::endogenous)
    error(name + " is not endogenous.");
}

void
ParsingDriver::check_symbol_is_exogenous(const string &name, bool allow_exo_det)
{
  check_symbol_existence(name);
  switch (mod_file->symbol_table.getType(name))
    {
    case SymbolType::exogenous:
      break;
    case SymbolType::exogenousDet:
      if (!allow_exo_det)
        error(name + " is an exogenous deterministic.");
      break;
    default:
      error(name + " is not exogenous.");
    }
}

void
ParsingDriver::set_std_prior(string name, string subsample_name)
{
  check_symbol_is_endogenous_or_exogenous(name, false);
  check_subsample_declaration_exists(name, subsample_name);
  mod_file->addStatement(make_unique<StdPriorStatement>(move(name), move(subsample_name),
                                                        prior_shape, prior_variance,
                                                        move(options_list), mod_file->symbol_table));
  options_list.clear();
  set_prior_variance();
  prior_shape = PriorDistributions::noShape;
}

void
ParsingDriver::set_std_options(string name, string subsample_name)
{
  check_symbol_is_endogenous_or_exogenous(name, false);
  check_subsample_declaration_exists(name, subsample_name);
  mod_file->addStatement(make_unique<StdOptionsStatement>(move(name), move(subsample_name),
                                                          move(options_list), mod_file->symbol_table));
  options_list.clear();
}

void
ParsingDriver::set_corr_prior(string name1, string name2, string subsample_name)
{
  check_symbol_is_endogenous_or_exogenous(name1, false);
  check_symbol_is_endogenous_or_exogenous(name2, false);
  check_subsample_declaration_exists(name1, name2, subsample_name);
  mod_file->addStatement(make_unique<CorrPriorStatement>(move(name1), move(name2),
                                                         move(subsample_name), prior_shape, prior_variance,
                                                         move(options_list), mod_file->symbol_table));
  options_list.clear();
  set_prior_variance();
  prior_shape = PriorDistributions::noShape;
}

void
ParsingDriver::set_corr_options(string name1, string name2, string subsample_name)
{
  check_symbol_is_endogenous_or_exogenous(name1, false);
  check_symbol_is_endogenous_or_exogenous(name2, false);
  check_subsample_declaration_exists(name1, name2, subsample_name);
  mod_file->addStatement(make_unique<CorrOptionsStatement>(move(name1), move(name2),
                                                           move(subsample_name), move(options_list),
                                                           mod_file->symbol_table));
  options_list.clear();
}

void
ParsingDriver::run_estimation(vector<string> symbol_list)
{
  mod_file->addStatement(make_unique<EstimationStatement>(mod_file->symbol_table, move(symbol_list),
                                                          move(options_list)));
  options_list.clear();
}

void
ParsingDriver::dynare_sensitivity()
{
  mod_file->addStatement(make_unique<DynareSensitivityStatement>(move(options_list)));
  options_list.clear();
}

void
ParsingDriver::check_varobs()
{
  if (mod_file->symbol_table.observedVariablesNbr() > 0)
    error("varobs: you cannot have several 'varobs' statements in the same MOD file");
}

void
ParsingDriver::add_varobs(const string &name)
{
  check_symbol_is_endogenous(name);
  int symb_id = mod_file->symbol_table.getID(name);
  mod_file->symbol_table.addObservedVariable(symb_id);
}

void
ParsingDriver::check_varexobs()
{
  if (mod_file->symbol_table.observedExogenousVariablesNbr() > 0)
    error("varexobs: you cannot have several 'varexobs' statements in the same MOD file");
}

void
ParsingDriver::add_varexobs(const string &name)
{
  check_symbol_existence(name);
  int symb_id = mod_file->symbol_table.getID(name);
  if (mod_file->symbol_table.getType(symb_id) != SymbolType::exogenous)
    error("varexobs: " + name + " is not an exogenous variable");
  mod_file->symbol_table.addObservedExogenousVariable(symb_id);
}

void
ParsingDriver::set_trends()
{
  mod_file->addStatement(make_unique<ObservationTrendsStatement>(move(trend_elements),
                                                                 mod_file->symbol_table));
  trend_elements.clear();
}

void
ParsingDriver::set_deterministic_trends()
{
  mod_file->addStatement(make_unique<DeterministicTrendsStatement>(move(trend_elements),
                                                                   mod_file->symbol_table));
  trend_elements.clear();
}

void
ParsingDriver::set_trend_element(string arg1, expr_t arg2)
{
  check_symbol_existence(arg1);
  if (trend_elements.contains(arg1))
    error("observation_trends: " + arg1 + " declared twice");
  trend_elements[move(arg1)] = arg2;
}

void
ParsingDriver::set_filter_initial_state()
{
  mod_file->addStatement(make_unique<FilterInitialStateStatement>(move(filter_initial_state_elements),
                                                                  mod_file->symbol_table));
  filter_initial_state_elements.clear();
}

void
ParsingDriver::set_filter_initial_state_element(const string &name, const string &lag, expr_t rhs)
{
  check_symbol_existence(name);
  int symb_id = mod_file->symbol_table.getID(name);
  SymbolType type = mod_file->symbol_table.getType(symb_id);
  int ilag = stoi(lag);

  if (type != SymbolType::endogenous
      && type != SymbolType::exogenous
      && type != SymbolType::exogenousDet)
    error("filter_initial_state: " + name + " should be an endogenous or exogenous variable");

  if ((type == SymbolType::exogenous || type == SymbolType::exogenousDet) && ilag == 0)
    error("filter_initial_state: exogenous variable " + name + " must be provided with a lag");

  if (filter_initial_state_elements.contains({ symb_id, ilag }))
    error("filter_initial_state: (" + name + ", " + lag + ") declared twice");

  if (mod_file->dynamic_model.minLagForSymbol(symb_id) > ilag - 1)
    error("filter_initial_state: variable " + name + " does not appear in the model with the lag " + to_string(ilag-1) + " (see the reference manual for the timing convention in 'filter_initial_state')");

  filter_initial_state_elements[{ symb_id, ilag }] = rhs;
}

void
ParsingDriver::set_optim_weights(string name, expr_t value)
{
  check_symbol_is_endogenous(name);
  if (var_weights.contains(name))
    error("optim_weights: " + name + " declared twice");
  var_weights[move(name)] = move(value);
}

void
ParsingDriver::set_optim_weights(const string &name1, const string &name2, expr_t value)
{
  check_symbol_is_endogenous(name1);
  check_symbol_is_endogenous(name2);

  pair covar_key{name1, name2};

  if (covar_weights.contains(covar_key))
    error("optim_weights: pair of variables (" + name1 + ", " + name2
          + ") declared twice");

  covar_weights[covar_key] = value;
}

void
ParsingDriver::optim_weights()
{
  mod_file->addStatement(make_unique<OptimWeightsStatement>(move(var_weights), move(covar_weights),
                                                            mod_file->symbol_table));
  var_weights.clear();
  covar_weights.clear();
}

void
ParsingDriver::set_osr_params(vector<string> symbol_list)
{
  mod_file->addStatement(make_unique<OsrParamsStatement>(move(symbol_list), mod_file->symbol_table));
}

void
ParsingDriver::run_osr(vector<string> symbol_list)
{
  mod_file->addStatement(make_unique<OsrStatement>(move(symbol_list), move(options_list),
                                                   mod_file->symbol_table));
  options_list.clear();
}

void
ParsingDriver::run_dynatype(string filename, vector<string> symbol_list)
{
  mod_file->addStatement(make_unique<DynaTypeStatement>(move(symbol_list), move(filename),
                                                        mod_file->symbol_table));
}

void
ParsingDriver::run_dynasave(string filename, vector<string> symbol_list)
{
  mod_file->addStatement(make_unique<DynaSaveStatement>(move(symbol_list), move(filename),
                                                        mod_file->symbol_table));
}

void
ParsingDriver::run_load_params_and_steady_state(const string &filename)
{
  mod_file->addStatement(make_unique<LoadParamsAndSteadyStateStatement>(filename, mod_file->symbol_table, warnings));
}

void
ParsingDriver::run_save_params_and_steady_state(string filename)
{
  mod_file->addStatement(make_unique<SaveParamsAndSteadyStateStatement>(move(filename)));
}

void
ParsingDriver::run_identification()
{
  mod_file->addStatement(make_unique<IdentificationStatement>(move(options_list)));
  options_list.clear();
}

void
ParsingDriver::add_mc_filename(string filename, string prior)
{
  for (auto &it : filename_list)
    if (it.first == filename)
      error("model_comparison: filename " + filename + " declared twice");
  filename_list.emplace_back(move(filename), move(prior));
}

void
ParsingDriver::run_model_comparison()
{
  mod_file->addStatement(make_unique<ModelComparisonStatement>(move(filename_list),
                                                               move(options_list)));
  filename_list.clear();
  options_list.clear();
}

void
ParsingDriver::begin_planner_objective()
{
  planner_objective = make_unique<PlannerObjective>(mod_file->symbol_table,
                                                    mod_file->num_constants,
                                                    mod_file->external_functions_table);
  set_current_data_tree(planner_objective.get());
}

void
ParsingDriver::end_planner_objective(expr_t expr)
{
  // Add equation corresponding to expression
  expr_t eq = model_tree->AddEqual(expr, model_tree->Zero);
  model_tree->addEquation(eq, location.begin.line);

  mod_file->addStatement(make_unique<PlannerObjectiveStatement>(*planner_objective));

  // Handle undeclared variables (see #81)
  bool exit_after_write = false;
  if (undeclared_model_variable_errors.size() > 0)
    for (auto &it : undeclared_model_variable_errors)
      {
        if (nostrict)
          warning(it.second);
        else
          {
            exit_after_write = true;
            cerr << it.second << endl;
          }
      }
  undeclared_model_variable_errors.clear();
  if (exit_after_write)
    exit(EXIT_FAILURE);

  reset_data_tree();
}

void
ParsingDriver::ramsey_model()
{
  // Some checks to ensure correct error messages (see #90)
  if (ramsey_policy_seen)
    error("A 'ramsey_model' statement cannot follow a 'ramsey_policy' statement.");
  if (ramsey_model_seen)
    error("Several 'ramsey_model' statements cannot appear in a given .mod file.");
  ramsey_model_seen = true;

  if (!mod_file->symbol_table.exists("optimal_policy_discount_factor"))
    {
      if (!planner_discount)
        planner_discount = data_tree->One;
      declare_parameter("optimal_policy_discount_factor", planner_discount_latex_name);
      init_param("optimal_policy_discount_factor", planner_discount);
    }
  else if (planner_discount)
    error("ramsey_model: the 'planner_discount' option cannot be used when the 'optimal_policy_discount_factor' parameter is explicitly declared.");

  // Check that instruments are declared endogenous (#72)
  if (options_list.contains("instruments"))
    for (const auto &s : options_list.get<OptionsList::SymbolListVal>("instruments").getSymbols())
      check_symbol_is_endogenous(s);

  mod_file->addStatement(make_unique<RamseyModelStatement>(move(options_list)));
  options_list.clear();
  planner_discount = nullptr;
  planner_discount_latex_name.clear();
}

void
ParsingDriver::ramsey_policy(vector<string> symbol_list)
{
  warning("The 'ramsey_policy' statement is deprecated. Please use 'ramsey_model', 'stoch_simul', and 'evaluate_planner_objective' instead.");

  // Some checks to ensure correct error messages (see #90)
  if (ramsey_model_seen)
    error("A 'ramsey_policy' statement cannot follow a 'ramsey_model' statement.");
  if (ramsey_policy_seen)
    error("Several 'ramsey_policy' statements cannot appear in a given .mod file.");
  ramsey_policy_seen = true;

  if (!mod_file->symbol_table.exists("optimal_policy_discount_factor"))
    {
      if (!planner_discount)
        planner_discount = data_tree->One;
      declare_parameter("optimal_policy_discount_factor");
      init_param("optimal_policy_discount_factor", planner_discount);
    }
  else if (planner_discount)
    error("ramsey_policy: the 'planner_discount' option cannot be used when the 'optimal_policy_discount_factor' parameter is explicitly declared.");

  // Check that instruments are declared endogenous (#72)
  if (options_list.contains("instruments"))
    for (const auto &s : options_list.get<OptionsList::SymbolListVal>("instruments").getSymbols())
      check_symbol_is_endogenous(s);

  mod_file->addStatement(make_unique<RamseyPolicyStatement>(move(symbol_list), move(options_list),
                                                            mod_file->symbol_table));
  options_list.clear();
  planner_discount = nullptr;
}

void
ParsingDriver::evaluate_planner_objective()
{
  mod_file->addStatement(make_unique<EvaluatePlannerObjectiveStatement>(move(options_list)));
  options_list.clear();
}

void
ParsingDriver::occbin_setup()
{
  mod_file->addStatement(make_unique<OccbinSetupStatement>(move(options_list)));
  options_list.clear();
}

void
ParsingDriver::occbin_solver()
{
  mod_file->addStatement(make_unique<OccbinSolverStatement>(move(options_list)));
  options_list.clear();
}

void
ParsingDriver::occbin_write_regimes()
{
  mod_file->addStatement(make_unique<OccbinWriteRegimesStatement>(move(options_list)));
  options_list.clear();
}

void
ParsingDriver::occbin_graph(vector<string> symbol_list)
{
  mod_file->addStatement(make_unique<OccbinGraphStatement>(move(symbol_list), move(options_list)));
  options_list.clear();
}

void
ParsingDriver::discretionary_policy(vector<string> symbol_list)
{
  /* The logic here is different from “ramsey_policy” and “ramsey_model”,
     because we want to allow several instances of “discretionary_policy” in
     the same .mod file. */
  if (!mod_file->symbol_table.exists("optimal_policy_discount_factor"))
    declare_parameter("optimal_policy_discount_factor");

  if (!planner_discount)
    planner_discount = data_tree->One;
  init_param("optimal_policy_discount_factor", planner_discount);

  // Check that instruments are declared endogenous (#72)
  if (options_list.contains("instruments"))
    for (const auto &s : options_list.get<OptionsList::SymbolListVal>("instruments").getSymbols())
      check_symbol_is_endogenous(s);

  mod_file->addStatement(make_unique<DiscretionaryPolicyStatement>(move(symbol_list),
                                                                   move(options_list),
                                                                   mod_file->symbol_table));
  options_list.clear();
  planner_discount = nullptr;
}

void
ParsingDriver::write_latex_dynamic_model(bool write_equation_tags)
{
  mod_file->addStatement(make_unique<WriteLatexDynamicModelStatement>(mod_file->dynamic_model, write_equation_tags));
}

void
ParsingDriver::write_latex_static_model(bool write_equation_tags)
{
  mod_file->addStatement(make_unique<WriteLatexStaticModelStatement>(mod_file->static_model, write_equation_tags));
}

void
ParsingDriver::write_latex_original_model(bool write_equation_tags)
{
  mod_file->addStatement(make_unique<WriteLatexOriginalModelStatement>(mod_file->original_model, write_equation_tags));
}

void
ParsingDriver::write_latex_steady_state_model()
{
  mod_file->addStatement(make_unique<WriteLatexSteadyStateModelStatement>(mod_file->steady_state_model));
}

void
ParsingDriver::bvar_density(const string &maxnlags)
{
  mod_file->addStatement(make_unique<BVARDensityStatement>(stoi(maxnlags), move(options_list)));
  options_list.clear();
}

void
ParsingDriver::bvar_forecast(const string &nlags)
{
  mod_file->addStatement(make_unique<BVARForecastStatement>(stoi(nlags), move(options_list)));
  options_list.clear();
}

void
ParsingDriver::bvar_irf(const string &nirf, string identificationname)
{
  mod_file->addStatement(make_unique<BVARIRFStatement>(stoi(nirf), move(identificationname)));
}

void
ParsingDriver::sbvar()
{
  mod_file->addStatement(make_unique<SBVARStatement>(move(options_list)));
  options_list.clear();
}

void
ParsingDriver::ms_estimation()
{
  mod_file->addStatement(make_unique<MSSBVAREstimationStatement>(move(options_list)));
  options_list.clear();
}

void
ParsingDriver::ms_simulation()
{
  mod_file->addStatement(make_unique<MSSBVARSimulationStatement>(move(options_list)));
  options_list.clear();
}

void
ParsingDriver::ms_compute_mdd()
{
  mod_file->addStatement(make_unique<MSSBVARComputeMDDStatement>(move(options_list)));
  options_list.clear();
}

void
ParsingDriver::ms_compute_probabilities()
{
  mod_file->addStatement(make_unique<MSSBVARComputeProbabilitiesStatement>(move(options_list)));
  options_list.clear();
}

void
ParsingDriver::ms_irf(vector<string> symbol_list)
{
  mod_file->addStatement(make_unique<MSSBVARIrfStatement>(move(symbol_list), move(options_list),
                                                          mod_file->symbol_table));
  options_list.clear();
}

void
ParsingDriver::ms_forecast()
{
  mod_file->addStatement(make_unique<MSSBVARForecastStatement>(move(options_list)));
  options_list.clear();
}

void
ParsingDriver::ms_variance_decomposition()
{
  mod_file->addStatement(make_unique<MSSBVARVarianceDecompositionStatement>(move(options_list)));
  options_list.clear();
}

void
ParsingDriver::svar()
{
  bool has_coefficients = options_list.contains("ms.coefficients"),
    has_variances = options_list.contains("ms.variances"),
    has_constants = options_list.contains("ms.constants");
  if (!has_coefficients && !has_variances && !has_constants)
    error("You must pass one of 'coefficients', 'variances', or 'constants'.");

  if ((has_coefficients && has_variances) || (has_variances && has_constants)
      || (has_coefficients && has_constants))
    error("You may only pass one of 'coefficients', 'variances', or 'constants'.");

  try
    {
      if (stoi(options_list.get<OptionsList::NumVal>("ms.chain")) <= 0)
        error("The value passed to the 'chain' option must be greater than zero.");
    }
  catch (OptionsList::UnknownOptionException &)
    {
      error("A 'chain' option must be passed to the 'svar' statement.");
    }

  if (options_list.contains("ms.equations"))
    for (int viit : options_list.get<vector<int>>("ms.equations"))
      if (viit <= 0)
        error("The value(s) passed to the 'equations' option must be greater than zero.");

  mod_file->addStatement(make_unique<SvarStatement>(move(options_list)));
  options_list.clear();
}

void
ParsingDriver::markov_switching()
{
  try
    {
      if (stoi(options_list.get<OptionsList::NumVal>("ms.chain")) <= 0)
        error("The value passed to the chain option must be greater than zero.");
      if (stoi(options_list.get<OptionsList::NumVal>("ms.number_of_regimes")) <= 0)
        error("The value passed to the number_of_regimes option must be greater than zero.");
      options_list.get<OptionsList::NumVal>("ms.duration"); // Just check its presence
    }
  catch (OptionsList::UnknownOptionException &e)
    {
      error("A '" + e.name.substr(3) + "' option must be passed to the 'markov_switching' statement.");
    }

  mod_file->addStatement(make_unique<MarkovSwitchingStatement>(move(options_list)));
  options_list.clear();
}

void
ParsingDriver::shock_decomposition(vector<string> symbol_list)
{
  mod_file->addStatement(make_unique<ShockDecompositionStatement>(move(symbol_list),
                                                                  move(options_list),
                                                                  mod_file->symbol_table));
  options_list.clear();
}

void
ParsingDriver::realtime_shock_decomposition(vector<string> symbol_list)
{
  mod_file->addStatement(make_unique<RealtimeShockDecompositionStatement>(move(symbol_list),
                                                                          move(options_list),
                                                                          mod_file->symbol_table));
  options_list.clear();
}

void
ParsingDriver::plot_shock_decomposition(vector<string> symbol_list)
{
  mod_file->addStatement(make_unique<PlotShockDecompositionStatement>(move(symbol_list),
                                                                      move(options_list),
                                                                      mod_file->symbol_table));
  options_list.clear();
}

void
ParsingDriver::initial_condition_decomposition(vector<string> symbol_list)
{
  mod_file->addStatement(make_unique<InitialConditionDecompositionStatement>(move(symbol_list),
                                                                             move(options_list),
                                                                             mod_file->symbol_table));
  options_list.clear();
}

void
ParsingDriver::squeeze_shock_decomposition(vector<string> symbol_list)
{
  mod_file->addStatement(make_unique<SqueezeShockDecompositionStatement>(move(symbol_list),
                                                                         mod_file->symbol_table));
}

void
ParsingDriver::conditional_forecast()
{
  mod_file->addStatement(make_unique<ConditionalForecastStatement>(move(options_list)));
  options_list.clear();
}

void
ParsingDriver::plot_conditional_forecast(const optional<string> &periods, vector<string> symbol_list)
{
  optional<int> iperiods;
  if (periods)
    iperiods = stoi(*periods);
  mod_file->addStatement(make_unique<PlotConditionalForecastStatement>(move(iperiods), move(symbol_list),
                                                                       mod_file->symbol_table));
}

void
ParsingDriver::conditional_forecast_paths()
{
  mod_file->addStatement(make_unique<ConditionalForecastPathsStatement>(move(det_shocks),
                                                                        mod_file->symbol_table));
  det_shocks.clear();
  if (!learnt_shocks_add.empty())
    error("conditional_forecast_paths: 'add' keyword not allowed");
  if (!learnt_shocks_multiply.empty())
    error("conditional_forecast_paths: 'multiply' keyword not allowed");
}

void
ParsingDriver::calib_smoother(vector<string> symbol_list)
{
  mod_file->addStatement(make_unique<CalibSmootherStatement>(move(symbol_list), move(options_list),
                                                             mod_file->symbol_table));
  options_list.clear();
}

void
ParsingDriver::extended_path()
{
  mod_file->addStatement(make_unique<ExtendedPathStatement>(move(options_list)));
  options_list.clear();
}

expr_t
ParsingDriver::add_model_equal(expr_t arg1, expr_t arg2)
{
  expr_t id = model_tree->AddEqual(arg1, arg2);

  if (eq_tags.contains("static"))
    {
      // If the equation is tagged [static]
      if (!id->isInStaticForm())
        error("An equation tagged [static] cannot contain leads, lags, expectations or STEADY_STATE operators");

      dynamic_model->addStaticOnlyEquation(id, location.begin.line, eq_tags);
    }
  else if (eq_tags.contains("bind") || eq_tags.contains("relax"))
    {
      // If the equation has a “bind” or “relax” tag (occbin case)
      if (!eq_tags.contains("name"))
        error("An equation with a 'bind' or 'relax' tag must have a 'name' tag");
      auto regimes_bind = DataTree::strsplit(eq_tags["bind"], ',');
      auto regimes_relax = DataTree::strsplit(eq_tags["relax"], ',');
      auto regimes_all = regimes_bind;
      regimes_all.insert(regimes_all.end(), regimes_relax.begin(), regimes_relax.end()); // Concatenate the two vectors
      for (const auto &regime : regimes_all)
        {
          if (!isSymbolIdentifier(regime))
            error("The string '" + regime + "' is not a valid Occbin regime name (contains unauthorized characters)");
          string param_name = buildOccbinBindParamName(regime);
          try
            {
              if (mod_file->symbol_table.getType(param_name) != SymbolType::parameter)
                error("The name '" + param_name + "' is already used. Please use another name for Occbin regime '" + regime + "'");
            }
          catch (SymbolTable::UnknownSymbolNameException &e)
            {
              // Declare and initialize the new parameter
              int symb_id = mod_file->symbol_table.addSymbol(param_name, SymbolType::parameter);
              mod_file->addStatement(make_unique<InitParamStatement>(symb_id, dynamic_model->Zero, mod_file->symbol_table));
            }
        }
      eq_tags.erase("bind");
      eq_tags.erase("relax");
      dynamic_model->addOccbinEquation(id, location.begin.line, eq_tags, regimes_bind, regimes_relax);
    }
  else // General case
    model_tree->addEquation(id, location.begin.line, eq_tags);

  eq_tags.clear();
  return id;
}

expr_t
ParsingDriver::add_model_equal_with_zero_rhs(expr_t arg)
{
  return add_model_equal(arg, model_tree->Zero);
}

void
ParsingDriver::model_local_variable(const vector<pair<string, string>> &symbol_list)
{
  for (auto &[name, tex_name] : symbol_list)
    declare_symbol(name, SymbolType::modelLocalVariable, tex_name, {});
}

void
ParsingDriver::declare_and_init_model_local_variable(const string &name, expr_t rhs)
{
  int symb_id;
  try
    {
      symb_id = mod_file->symbol_table.addSymbol(name, SymbolType::modelLocalVariable);
    }
  catch (SymbolTable::AlreadyDeclaredException &e)
    {
      /* It can have already been declared in a steady_state_model block or
         model_local_variable statement, check that it is indeed a
         ModelLocalVariable */
      symb_id = mod_file->symbol_table.getID(name);
      if (mod_file->symbol_table.getType(symb_id) != SymbolType::modelLocalVariable)
        error(name + " has wrong type or was already used on the right-hand side. You cannot use it on the left-hand side of a pound ('#') expression");
    }

  try
    {
      model_tree->AddLocalVariable(symb_id, rhs);
    }
  catch (DataTree::LocalVariableException &e)
    {
      error("Local model variable " + name + " declared twice.");
    }
}

void
ParsingDriver::change_type(SymbolType new_type, const vector<string> &symbol_list)
{
  for (auto &it : symbol_list)
    {
      int id;
      try
        {
          id = mod_file->symbol_table.getID(it);
        }
      catch (SymbolTable::UnknownSymbolNameException &e)
        {
          error("Unknown variable " + it);
        }

      // Check if symbol already used in a VariableNode
      if (mod_file->expressions_tree.isSymbolUsed(id)
          || mod_file->dynamic_model.isSymbolUsed(id))
        error("You cannot modify the type of symbol " + it + " after having used it in an expression");

      mod_file->symbol_table.changeType(id, new_type);
    }
}

expr_t
ParsingDriver::add_plus(expr_t arg1, expr_t arg2)
{
  return data_tree->AddPlus(arg1, arg2);
}

expr_t
ParsingDriver::add_minus(expr_t arg1, expr_t arg2)
{
  return data_tree->AddMinus(arg1, arg2);
}

expr_t
ParsingDriver::add_uminus(expr_t arg1)
{
  return data_tree->AddUMinus(arg1);
}

expr_t
ParsingDriver::add_times(expr_t arg1, expr_t arg2)
{
  return data_tree->AddTimes(arg1, arg2);
}

expr_t
ParsingDriver::add_divide(expr_t arg1, expr_t arg2)
{
  try
    {
      return data_tree->AddDivide(arg1, arg2);
    }
  catch (DataTree::DivisionByZeroException)
    {
      error("Division by zero error encountered when reading model from .mod file");
    }
}

expr_t
ParsingDriver::add_less(expr_t arg1, expr_t arg2)
{
  return data_tree->AddLess(arg1, arg2);
}

expr_t
ParsingDriver::add_greater(expr_t arg1, expr_t arg2)
{
  return data_tree->AddGreater(arg1, arg2);
}

expr_t
ParsingDriver::add_less_equal(expr_t arg1, expr_t arg2)
{
  return data_tree->AddLessEqual(arg1, arg2);
}

expr_t
ParsingDriver::add_greater_equal(expr_t arg1, expr_t arg2)
{
  return data_tree->AddGreaterEqual(arg1, arg2);
}

expr_t
ParsingDriver::add_equal_equal(expr_t arg1, expr_t arg2)
{
  return data_tree->AddEqualEqual(arg1, arg2);
}

expr_t
ParsingDriver::add_different(expr_t arg1, expr_t arg2)
{
  return data_tree->AddDifferent(arg1, arg2);
}

expr_t
ParsingDriver::add_power(expr_t arg1, expr_t arg2)
{
  return data_tree->AddPower(arg1, arg2);
}

expr_t
ParsingDriver::add_expectation(const string &arg1, expr_t arg2)
{
  if (data_tree == occbin_constraints_tree.get())
    error("The 'expectation' operator is forbidden in 'occbin_constraints'.");

  return data_tree->AddExpectation(stoi(arg1), arg2);
}

expr_t
ParsingDriver::add_var_expectation(const string &model_name)
{
  if (data_tree == occbin_constraints_tree.get())
    error("The 'var_expectation' operator is forbidden in 'occbin_constraints'.");

  return data_tree->AddVarExpectation(model_name);
}

expr_t
ParsingDriver::add_pac_expectation(const string &model_name)
{
  if (data_tree == occbin_constraints_tree.get())
    error("The 'pac_expectation' operator is forbidden in 'occbin_constraints'.");

  return data_tree->AddPacExpectation(model_name);
}

expr_t
ParsingDriver::add_pac_target_nonstationary(const string &model_name)
{
  if (data_tree == occbin_constraints_tree.get())
    error("The 'pac_target_nonstationary' operator is forbidden in 'occbin_constraints'.");

  return data_tree->AddPacTargetNonstationary(model_name);
}

void
ParsingDriver::begin_pac_model()
{
  parsing_pac_model = true;
  pac_growth = nullptr;
  pac_auxname.clear();
  pac_kind = PacTargetKind::unspecified;
}

void
ParsingDriver::pac_model()
{
  try
    {
      auto discount {options_list.get<OptionsList::StringVal>("pac.discount")};
      check_symbol_is_parameter(discount);
      mod_file->pac_model_table.addPacModel(options_list.get<OptionsList::StringVal>("pac.model_name"),
                                            options_list.get_if<OptionsList::StringVal>("pac.aux_model_name").value_or(OptionsList::StringVal{}),
                                            move(discount), pac_growth,
                                            pac_auxname, pac_kind);
    }
  catch (OptionsList::UnknownOptionException &e)
    {
      error("You must pass the '" + e.name.substr(4) + "' option to the 'pac_model' statement.");
    }

  options_list.clear();
  parsing_pac_model = false;
}

void
ParsingDriver::set_pac_growth(expr_t pac_growth_arg)
{
  pac_growth = pac_growth_arg;
  reset_data_tree();
}

void
ParsingDriver::set_pac_auxname(string auxname)
{
  pac_auxname = move(auxname);
}

void
ParsingDriver::set_pac_kind(PacTargetKind kind)
{
  pac_kind = kind;
}

expr_t
ParsingDriver::add_exp(expr_t arg1)
{
  return data_tree->AddExp(arg1);
}

expr_t
ParsingDriver::add_diff(expr_t arg1)
{
  return data_tree->AddDiff(arg1);
}

expr_t
ParsingDriver::add_adl(expr_t arg1, const string &name, const string &lag)
{
  vector<int> lags(stoi(lag));
  iota(lags.begin(), lags.end(), 1);
  return add_adl(arg1, name, lags);
}

expr_t
ParsingDriver::add_adl(expr_t arg1, const string &name, const vector<int> &lags)
{
  expr_t id = data_tree->AddAdl(arg1, name, lags);

  // Declare parameters here so that parameters can be initialized after the model block
  for (auto i : lags)
    declare_parameter(name + "_lag_" + to_string(i));

  return id;
}

expr_t
ParsingDriver::add_log(expr_t arg1)
{
  return data_tree->AddLog(arg1);
}

expr_t
ParsingDriver::add_log10(expr_t arg1)
{
  return data_tree->AddLog10(arg1);
}

expr_t
ParsingDriver::add_cos(expr_t arg1)
{
  return data_tree->AddCos(arg1);
}

expr_t
ParsingDriver::add_sin(expr_t arg1)
{
  return data_tree->AddSin(arg1);
}

expr_t
ParsingDriver::add_tan(expr_t arg1)
{
  return data_tree->AddTan(arg1);
}

expr_t
ParsingDriver::add_acos(expr_t arg1)
{
  return data_tree->AddAcos(arg1);
}

expr_t
ParsingDriver::add_asin(expr_t arg1)
{
  return data_tree->AddAsin(arg1);
}

expr_t
ParsingDriver::add_atan(expr_t arg1)
{
  return data_tree->AddAtan(arg1);
}

expr_t
ParsingDriver::add_cosh(expr_t arg1)
{
  return data_tree->AddCosh(arg1);
}

expr_t
ParsingDriver::add_sinh(expr_t arg1)
{
  return data_tree->AddSinh(arg1);
}

expr_t
ParsingDriver::add_tanh(expr_t arg1)
{
  return data_tree->AddTanh(arg1);
}

expr_t
ParsingDriver::add_acosh(expr_t arg1)
{
  return data_tree->AddAcosh(arg1);
}

expr_t
ParsingDriver::add_asinh(expr_t arg1)
{
  return data_tree->AddAsinh(arg1);
}

expr_t
ParsingDriver::add_atanh(expr_t arg1)
{
  return data_tree->AddAtanh(arg1);
}

expr_t
ParsingDriver::add_sqrt(expr_t arg1)
{
  return data_tree->AddSqrt(arg1);
}

expr_t
ParsingDriver::add_cbrt(expr_t arg1)
{
  return data_tree->AddCbrt(arg1);
}

expr_t
ParsingDriver::add_abs(expr_t arg1)
{
  return data_tree->AddAbs(arg1);
}

expr_t
ParsingDriver::add_sign(expr_t arg1)
{
  return data_tree->AddSign(arg1);
}

expr_t
ParsingDriver::add_max(expr_t arg1, expr_t arg2)
{
  return data_tree->AddMax(arg1, arg2);
}

expr_t
ParsingDriver::add_min(expr_t arg1, expr_t arg2)
{
  return data_tree->AddMin(arg1, arg2);
}

expr_t
ParsingDriver::add_normcdf(expr_t arg1, expr_t arg2, expr_t arg3)
{
  return data_tree->AddNormcdf(arg1, arg2, arg3);
}

expr_t
ParsingDriver::add_normcdf(expr_t arg)
{
  return add_normcdf(arg, data_tree->Zero, data_tree->One);
}

expr_t
ParsingDriver::add_normpdf(expr_t arg1, expr_t arg2, expr_t arg3)
{
  return data_tree->AddNormpdf(arg1, arg2, arg3);
}

expr_t
ParsingDriver::add_normpdf(expr_t arg)
{
  return add_normpdf(arg, data_tree->Zero, data_tree->One);
}

expr_t
ParsingDriver::add_erf(expr_t arg1)
{
  return data_tree->AddErf(arg1);
}

expr_t
ParsingDriver::add_erfc(expr_t arg1)
{
  return data_tree->AddErfc(arg1);
}

expr_t
ParsingDriver::add_steady_state(expr_t arg1)
{
  // Forbid exogenous variables, see dynare#825
  if (arg1->hasExogenous())
    error("Exogenous variables are not allowed in the context of the STEADY_STATE() operator.");

  return data_tree->AddSteadyState(arg1);
}

void
ParsingDriver::external_function_option(const string &name_option, const string &opt)
{
  if (name_option == "name")
    {
      if (opt.empty())
        error("An argument must be passed to the 'name' option of the external_function() statement.");
      declare_symbol(opt, SymbolType::externalFunction, "", {});
      current_external_function_id = mod_file->symbol_table.getID(opt);
    }
  else if (name_option == "first_deriv_provided")
    {
      if (opt.empty())
        current_external_function_options.firstDerivSymbID = ExternalFunctionsTable::IDSetButNoNameProvided;
      else
        {
          int symb_id = declare_symbol(opt, SymbolType::externalFunction, "", {});
          current_external_function_options.firstDerivSymbID = symb_id;
        }
    }
  else if (name_option == "second_deriv_provided")
    {
      if (opt.empty())
        current_external_function_options.secondDerivSymbID = ExternalFunctionsTable::IDSetButNoNameProvided;
      else
        {
          int symb_id = declare_symbol(opt, SymbolType::externalFunction, "", {});
          current_external_function_options.secondDerivSymbID = symb_id;
        }
    }
  else if (name_option == "nargs")
    current_external_function_options.nargs = stoi(opt);
  else
    error("Unexpected error in ParsingDriver::external_function_option(): Please inform Dynare Team.");
}

void
ParsingDriver::external_function()
{
  if (current_external_function_id == ExternalFunctionsTable::IDNotSet)
    error("The 'name' option must be passed to external_function().");

  if (current_external_function_options.secondDerivSymbID >= 0
      && current_external_function_options.firstDerivSymbID == ExternalFunctionsTable::IDNotSet)
    error("If the second derivative is provided to the external_function command, the first derivative must also be provided.");

  if (current_external_function_options.secondDerivSymbID == ExternalFunctionsTable::IDSetButNoNameProvided
      && current_external_function_options.firstDerivSymbID != ExternalFunctionsTable::IDSetButNoNameProvided)
    error("If the second derivative is provided in the top-level function, the first derivative must also be provided in that function.");

  mod_file->external_functions_table.addExternalFunction(current_external_function_id, current_external_function_options, true);
  reset_current_external_function_options();
}

void
ParsingDriver::push_external_function_arg_vector_onto_stack()
{
  stack_external_function_args.push({});
}

void
ParsingDriver::add_external_function_arg(expr_t arg)
{
  stack_external_function_args.top().push_back(arg);
}

optional<int>
ParsingDriver::is_there_one_integer_argument() const
{
  if (stack_external_function_args.top().size() != 1)
    return nullopt;

  auto numNode = dynamic_cast<NumConstNode *>(stack_external_function_args.top().front());
  auto unaryNode = dynamic_cast<UnaryOpNode *>(stack_external_function_args.top().front());

  if (!numNode && !unaryNode)
    return nullopt;

  eval_context_t ectmp;
  double model_var_arg;
  if (!unaryNode)
    {
      try
        {
          model_var_arg = numNode->eval(ectmp);
        }
      catch (ExprNode::EvalException &e)
        {
          return nullopt;
        }
    }
  else
    if (unaryNode->op_code != UnaryOpcode::uminus)
      return nullopt;
    else
      {
        try
          {
            model_var_arg = unaryNode->eval(ectmp);
          }
        catch (ExprNode::EvalException &e)
          {
            return nullopt;
          }
      }

  if (model_var_arg != floor(model_var_arg))
    return nullopt;
  return static_cast<int>(model_var_arg);
}

expr_t
ParsingDriver::add_model_var_or_external_function(const string &function_name, bool in_model_block)
{
  expr_t nid;
  if (mod_file->symbol_table.exists(function_name))
    if (mod_file->symbol_table.getType(function_name) != SymbolType::externalFunction)
      if (!in_model_block && !parsing_epilogue && !parsing_pac_model)
        {
          if (stack_external_function_args.top().size() > 0)
            error("Symbol " + function_name + " cannot take arguments.");
          else
            return add_expression_variable(function_name);
        }
      else
        { // e.g. model_var(lag) => ADD MODEL VARIABLE WITH LEAD (NumConstNode)/LAG (UnaryOpNode)
          if (undeclared_model_vars.contains(function_name))
            undeclared_model_variable_error("Unknown symbol: " + function_name, function_name);

          optional<int> rv{is_there_one_integer_argument()};
          if (!rv)
            model_error("Symbol " + function_name
                        +" is being treated as if it were a function (i.e., takes an argument that is not an integer).", "");

          nid = add_model_variable(mod_file->symbol_table.getID(function_name), *rv);
          stack_external_function_args.pop();
          return nid;
        }
    else
      { // e.g. this function has already been referenced (either ad hoc or through the external_function() statement
        // => check that the information matches previously declared info
        int symb_id = mod_file->symbol_table.getID(function_name);
        if (!mod_file->external_functions_table.exists(symb_id))
          error("Using a derivative of an external function (" + function_name + ") in the model block is currently not allowed.");

        if (in_model_block || parsing_epilogue)
          {
            if (mod_file->external_functions_table.getNargs(symb_id) == ExternalFunctionsTable::IDNotSet)
              error("Before using " + function_name
                    +"() in the model block, you must first declare it via the external_function() statement");
            else if (static_cast<int>(stack_external_function_args.top().size()) != mod_file->external_functions_table.getNargs(symb_id))
              error("The number of arguments passed to " + function_name
                    +"() does not match those of a previous call or declaration of this function.");
          }
      }
  else
    { //First time encountering this external function i.e., not previously declared or encountered
      if (parsing_epilogue)
        error("Variable " + function_name + " used in the epilogue block but was not declared.");

      if (in_model_block)
        {
          // Continue processing, noting that it was not declared
          // Processing will end at the end of the model block if nostrict was not passed
          undeclared_model_vars.insert(function_name);
          undeclared_model_variable_error("Unknown symbol: " + function_name, function_name);

          optional<int>rv{is_there_one_integer_argument()};
          if (rv)
            {
              // assume it's a lead/lagged variable
              int symb_id = declare_exogenous(function_name);
              return add_model_variable(symb_id, *rv);
            }
          else
            error("To use an external function (" + function_name
                  +") within the model block, you must first declare it via the external_function() statement.");
        }
      int symb_id = declare_symbol(function_name, SymbolType::externalFunction, "", {});
      current_external_function_options.nargs = stack_external_function_args.top().size();
      mod_file->external_functions_table.addExternalFunction(symb_id,
                                                             current_external_function_options, in_model_block);
      reset_current_external_function_options();
    }

  //By this point, we're sure that this function exists in the External Functions Table and is not a mod var
  int symb_id = mod_file->symbol_table.getID(function_name);
  nid = data_tree->AddExternalFunction(symb_id, stack_external_function_args.top());
  stack_external_function_args.pop();
  return nid;
}

void
ParsingDriver::add_native(string s)
{
  mod_file->addStatement(make_unique<NativeStatement>(move(s)));
}

void
ParsingDriver::add_native_remove_charset(string_view str, string_view token)
{
  size_t found = str.find(token);
  assert(found != string_view::npos);
  add_native(string{str.substr(0, found)});
}

void
ParsingDriver::add_verbatim(string s)
{
  mod_file->addStatement(make_unique<VerbatimStatement>(move(s)));
}

void
ParsingDriver::add_verbatim_remove_charset(string_view str, string_view token)
{
  size_t found = str.find(token);
  assert(found != string_view::npos);
  add_verbatim(string{str.substr(0, found)});
}

void
ParsingDriver::begin_steady_state_model()
{
  set_current_data_tree(&mod_file->steady_state_model);
}

void
ParsingDriver::add_steady_state_model_equal(const string &varname, expr_t expr)
{
  int id;
  try
    {
      id = mod_file->symbol_table.getID(varname);
    }
  catch (SymbolTable::UnknownSymbolNameException &e)
    {
      // Unknown symbol, declare it as a ModFileLocalVariable
      id = mod_file->symbol_table.addSymbol(varname, SymbolType::modFileLocalVariable);
    }

  if (SymbolType type = mod_file->symbol_table.getType(id);
      type != SymbolType::endogenous && type != SymbolType::modFileLocalVariable && type != SymbolType::parameter)
    error(varname + " has incorrect type");

  mod_file->steady_state_model.addDefinition(id, expr);
}

void
ParsingDriver::add_steady_state_model_equal_multiple(const vector<string> &symbol_list, expr_t expr)
{
  vector<int> ids;

  for (const auto &symb : symbol_list)
    {
      int id;
      try
        {
          id = mod_file->symbol_table.getID(symb);
        }
      catch (SymbolTable::UnknownSymbolNameException &e)
        {
          // Unknown symbol, declare it as a ModFileLocalVariable
          id = mod_file->symbol_table.addSymbol(symb, SymbolType::modFileLocalVariable);
        }
      if (SymbolType type = mod_file->symbol_table.getType(id);
          type != SymbolType::endogenous && type != SymbolType::modFileLocalVariable && type != SymbolType::parameter)
        error(symb + " has incorrect type");
      ids.push_back(id);
    }

  mod_file->steady_state_model.addMultipleDefinitions(ids, expr);
}

void
ParsingDriver::add_graph_format(string name)
{
  graph_formats.emplace_back(move(name));
}

void
ParsingDriver::process_graph_format_option()
{
  options_list.set("graph_format", OptionsList::SymbolListVal{move(graph_formats)});
  graph_formats.clear();
}

void
ParsingDriver::initial_condition_decomp_process_graph_format_option()
{
  options_list.set("initial_condition_decomp.graph_format", OptionsList::SymbolListVal{move(graph_formats)});
  graph_formats.clear();
}

void
ParsingDriver::plot_shock_decomp_process_graph_format_option()
{
  options_list.set("plot_shock_decomp.graph_format", OptionsList::SymbolListVal{move(graph_formats)});
  graph_formats.clear();
}

void
ParsingDriver::model_diagnostics()
{
  mod_file->addStatement(make_unique<ModelDiagnosticsStatement>());
}

void
ParsingDriver::add_parallel_local_file(string filename)
{
  mod_file->parallel_local_files.push_back(move(filename));
}

void
ParsingDriver::add_moment_calibration_item(const string &endo1, const string &endo2, string lags, const pair<expr_t, expr_t> &range)
{
  MomentCalibration::Constraint c;

  check_symbol_is_endogenous(endo1);
  c.endo1 = mod_file->symbol_table.getID(endo1);

  check_symbol_is_endogenous(endo2);
  c.endo2 = mod_file->symbol_table.getID(endo2);

  c.lags = move(lags);

  c.lower_bound = range.first;
  c.upper_bound = range.second;

  moment_calibration_constraints.push_back(move(c));
}

void
ParsingDriver::end_moment_calibration()
{
  mod_file->addStatement(make_unique<MomentCalibration>(move(moment_calibration_constraints),
                                                        mod_file->symbol_table));
  moment_calibration_constraints.clear();
}

void
ParsingDriver::add_irf_calibration_item(const string &endo, string periods, const string &exo, const pair<expr_t, expr_t> &range)
{
  IrfCalibration::Constraint c;

  check_symbol_is_endogenous(endo);
  c.endo = mod_file->symbol_table.getID(endo);

  c.periods = move(periods);

  check_symbol_existence(exo);
  c.exo = mod_file->symbol_table.getID(exo);
  if (mod_file->symbol_table.getType(exo) != SymbolType::exogenous)
    error("Variable " + endo + " is not an exogenous.");

  c.lower_bound = range.first;
  c.upper_bound = range.second;

  irf_calibration_constraints.push_back(move(c));
}

void
ParsingDriver::end_irf_calibration()
{
  mod_file->addStatement(make_unique<IrfCalibration>(move(irf_calibration_constraints),
                                                     mod_file->symbol_table,
                                                     move(options_list)));
  irf_calibration_constraints.clear();
  options_list.clear();
}

void
ParsingDriver::smoother2histval()
{
  mod_file->addStatement(make_unique<Smoother2histvalStatement>(move(options_list)));
  options_list.clear();
}

void
ParsingDriver::histval_file()
{
  mod_file->addStatement(make_unique<HistvalFileStatement>(move(options_list)));
  options_list.clear();
}

void
ParsingDriver::perfect_foresight_setup()
{
  mod_file->addStatement(make_unique<PerfectForesightSetupStatement>(move(options_list)));
  options_list.clear();
}

void
ParsingDriver::perfect_foresight_solver()
{
  mod_file->addStatement(make_unique<PerfectForesightSolverStatement>(move(options_list)));
  options_list.clear();
}

void
ParsingDriver::perfect_foresight_with_expectation_errors_setup()
{
  mod_file->addStatement(make_unique<PerfectForesightWithExpectationErrorsSetupStatement>(move(options_list)));
  options_list.clear();
}

void
ParsingDriver::perfect_foresight_with_expectation_errors_solver()
{
  mod_file->addStatement(make_unique<PerfectForesightWithExpectationErrorsSolverStatement>(move(options_list)));
  options_list.clear();
}

void
ParsingDriver::method_of_moments()
{
  mod_file->addStatement(make_unique<MethodOfMomentsStatement>(move(options_list)));
  options_list.clear();
}

void
ParsingDriver::prior_posterior_function(bool prior_func)
{
  mod_file->addStatement(make_unique<PriorPosteriorFunctionStatement>(static_cast<bool>(prior_func),
                                                                      move(options_list)));
  options_list.clear();
}

void
ParsingDriver::add_ramsey_constraints_statement()
{
  mod_file->addStatement(make_unique<RamseyConstraintsStatement>(mod_file->symbol_table,
                                                                 move(ramsey_constraints)));
  ramsey_constraints.clear();
}

void
ParsingDriver::ramsey_constraint_add_less(const string &name, const expr_t rhs)
{
  add_ramsey_constraint(name, BinaryOpcode::less, rhs);
}

void
ParsingDriver::ramsey_constraint_add_greater(const string &name, const expr_t rhs)
{
  add_ramsey_constraint(name, BinaryOpcode::greater, rhs);
}

void
ParsingDriver::ramsey_constraint_add_less_equal(const string &name, const expr_t rhs)
{
  add_ramsey_constraint(name, BinaryOpcode::lessEqual, rhs);
}

void
ParsingDriver::ramsey_constraint_add_greater_equal(const string &name, const expr_t rhs)
{
  add_ramsey_constraint(name, BinaryOpcode::greaterEqual, rhs);
}

void
ParsingDriver::add_ramsey_constraint(const string &name, BinaryOpcode op_code, const expr_t rhs)
{
  check_symbol_is_endogenous(name);
  int symb_id = mod_file->symbol_table.getID(name);

  RamseyConstraintsStatement::Constraint C;
  C.endo = symb_id;
  C.code = op_code;
  C.expression = rhs;
  ramsey_constraints.push_back(C);
}

void
ParsingDriver::add_shock_group_element(string name)
{
  check_symbol_existence(name);
  int symb_id = mod_file->symbol_table.getID(name);

  if (mod_file->symbol_table.getType(symb_id) != SymbolType::exogenous)
    error("shock_groups: " + name + " should be an exogenous variable");

  shock_group.push_back(move(name));
}

void
ParsingDriver::add_shock_group(string name)
{
  ShockGroupsStatement::Group G;
  G.name = move(name);
  G.list = shock_group;
  shock_groups.push_back(move(G));

  shock_group.clear();
}

void
ParsingDriver::end_shock_groups(string name)
{
  mod_file->addStatement(make_unique<ShockGroupsStatement>(move(shock_groups), move(name)));
  shock_groups.clear();
}

void
ParsingDriver::add_init2shocks(const string &endo_name, const string &exo_name)
{
  check_symbol_existence(endo_name);
  check_symbol_existence(exo_name);
  int symb_id_endo = mod_file->symbol_table.getID(endo_name);
  if (mod_file->symbol_table.getType(symb_id_endo) != SymbolType::endogenous)
    error("init2shocks: " + endo_name + " should be an endogenous variable");

  int symb_id_exo = mod_file->symbol_table.getID(exo_name);
  if (mod_file->symbol_table.getType(symb_id_exo) != SymbolType::exogenous)
    error("init2shocks: " + exo_name + " should be an exogenous variable");

  init2shocks.emplace_back(symb_id_endo, symb_id_exo);
}

void
ParsingDriver::end_init2shocks(string name)
{
  mod_file->addStatement(make_unique<Init2shocksStatement>(move(init2shocks), move(name),
                                                           mod_file->symbol_table));
  init2shocks.clear();
}

void
ParsingDriver::var_expectation_model()
{
  try
    {
      string v {options_list.get<OptionsList::StringVal>("variable")};
      if (var_expectation_model_expression)
        error("You can't pass both the 'variable' or the 'expression' options to the var_expectation_model statement.");
      var_expectation_model_expression = data_tree->AddVariable(mod_file->symbol_table.getID(v));
    }
  catch (OptionsList::UnknownOptionException &)
    {
      if (!var_expectation_model_expression)
        error("You must pass either the 'variable' or the 'expression' option to the var_expectation_model statement.");
    }

  if (var_expectation_model_discount)
    {
      VariableNode *var;
      if (!dynamic_cast<NumConstNode *>(var_expectation_model_discount)
          && !((var = dynamic_cast<VariableNode *>(var_expectation_model_discount))
               && var->get_type() == SymbolType::parameter))
        error("The discount factor must be a constant expression or a parameter");
    }
  else
    var_expectation_model_discount = data_tree->One;

  int time_shift { stoi(options_list.get_if<OptionsList::NumVal>("time_shift").value_or(OptionsList::NumVal{"0"})) };
  if (time_shift > 0)
    error("The 'time_shift' option must be a non-positive integer");

  try
    {
      mod_file->var_expectation_model_table.addVarExpectationModel(options_list.get<OptionsList::StringVal>("model_name"),
                                                                   var_expectation_model_expression,
                                                                   options_list.get<OptionsList::StringVal>("auxiliary_model_name"),
                                                                   options_list.get<OptionsList::NumVal>("horizon"),
                                                                   var_expectation_model_discount,
                                                                   time_shift);
    }
  catch (OptionsList::UnknownOptionException &e)
    {
      error("You must pass the '" + e.name + "' option to the 'var_expectation_model' statement.");
    }

  options_list.clear();
  var_expectation_model_discount = nullptr;
  var_expectation_model_expression = nullptr;
}

void
ParsingDriver::begin_matched_moments()
{
  set_current_data_tree(&mod_file->dynamic_model);
}

void
ParsingDriver::end_matched_moments(const vector<expr_t> &moments)
{
  vector<tuple<vector<int>, vector<int>, vector<int>>> parsed_moments;
  for (auto m : moments)
    try
      {
        vector<int> symb_ids, lags, powers;
        m->matchMatchedMoment(symb_ids, lags, powers);
        parsed_moments.emplace_back(move(symb_ids), move(lags), move(powers));
      }
    catch (ExprNode::MatchFailureException &e)
      {
        error("Matched moment expression has incorrect format: " + e.message);
      }
  mod_file->addStatement(make_unique<MatchedMomentsStatement>(mod_file->symbol_table,
                                                              move(parsed_moments)));

  reset_data_tree();
}

void
ParsingDriver::begin_occbin_constraints()
{
  /* We use a separate data tree, because we will add inequality constraints,
     and those would trigger the non-linearity warning in a stochastic context
     if added to the main DynamicModel tree. It also simplifies the
     enforcement of various constraints at parsing time. */
  occbin_constraints_tree = make_unique<DataTree>(mod_file->symbol_table,
                                                  mod_file->num_constants,
                                                  mod_file->external_functions_table,
                                                  false);
  set_current_data_tree(occbin_constraints_tree.get());
}

void
ParsingDriver::end_occbin_constraints(vector<tuple<string, BinaryOpNode *, BinaryOpNode *, expr_t, expr_t>> constraints)
{
  // Perform a few checks
  for (const auto &[name, bind, relax, error_bind, error_relax] : constraints)
    {
      string param_name = buildOccbinBindParamName(name);
      if (!mod_file->symbol_table.exists(param_name))
        error("No equation has been declared for regime '" + name + "'");
      if (!bind)
        error("The 'bind' expression is missing in constraint '" + name + "'");
    }

  mod_file->addStatement(make_unique<OccbinConstraintsStatement>(*occbin_constraints_tree,
                                                                 move(constraints)));

  reset_data_tree();
}

void
ParsingDriver::begin_pac_target_info(string name)
{
  pac_target_info_name = move(name);
  set_current_data_tree(&mod_file->dynamic_model);
}

void
ParsingDriver::end_pac_target_info()
{
  reset_data_tree();
}

void
ParsingDriver::set_pac_target_info_target(expr_t target)
{
  mod_file->pac_model_table.setTargetExpr(pac_target_info_name, target);
}

void
ParsingDriver::set_pac_target_info_auxname_target_nonstationary(string auxname)
{
  mod_file->pac_model_table.setTargetAuxnameNonstationary(pac_target_info_name, move(auxname));
}

void
ParsingDriver::add_pac_target_info_component(expr_t component_expr)
{
  get<0>(pac_target_info_component) = component_expr;
  mod_file->pac_model_table.addTargetComponent(pac_target_info_name, exchange(pac_target_info_component, {}));
}

void
ParsingDriver::set_pac_target_info_component_growth(expr_t growth)
{
  get<1>(pac_target_info_component) = growth;
}

void
ParsingDriver::set_pac_target_info_component_auxname(string auxname)
{
  get<2>(pac_target_info_component) = move(auxname);
}

void
ParsingDriver::set_pac_target_info_component_kind(PacTargetKind kind)
{
  get<3>(pac_target_info_component) = kind;
}

bool
ParsingDriver::isSymbolIdentifier(const string &str)
{
  if (str.empty())
    return false;
  auto myisalpha = [](char ch)
  {
    // We cannot use std::isalpha(), because it is locale-dependent
    return (ch >= 'a' && ch <= 'z') || (ch >= 'A' && ch <= 'Z');
  };
  if (!(myisalpha(str[0]) || str[0] == '_'))
    return false;
  for (size_t i = 1; i < str.size(); i++)
    if (!(myisalpha(str[i]) || isdigit(static_cast<unsigned char>(str[i])) || str[i] == '_'))
      return false;
  return true;
}

void
ParsingDriver::model_remove(const vector<pair<string, string>> &listed_eqs_by_tags)
{
  mod_file->dynamic_model.removeEquations(listed_eqs_by_tags, true, true);
}

void
ParsingDriver::begin_model_replace(const vector<pair<string, string>> &listed_eqs_by_tags)
{
  mod_file->dynamic_model.removeEquations(listed_eqs_by_tags, true, false);
  set_current_data_tree(&mod_file->dynamic_model);
}

void
ParsingDriver::var_remove(const vector<string> &symbol_list)
{
  for (const auto &name : symbol_list)
    {
      check_symbol_existence(name);
      int symb_id = mod_file->symbol_table.getID(name);
      mod_file->symbol_table.changeType(symb_id, SymbolType::excludedVariable);
    }
}

void
ParsingDriver::resid()
{
  mod_file->addStatement(make_unique<ResidStatement>(move(options_list)));
  options_list.clear();
}
