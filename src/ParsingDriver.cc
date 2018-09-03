/*
 * Copyright (C) 2003-2018 Dynare Team
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

#include <fstream>
#include <iostream>
#include <cassert>
#include <sstream>
#include <cmath>

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

  return (type != SymbolType::modFileLocalVariable && type != SymbolType::externalFunction);
}

void
ParsingDriver::check_symbol_existence_in_model_block(const string &name)
{
  if (!mod_file->symbol_table.exists(name)
      || undeclared_model_vars.find(name) != undeclared_model_vars.end())
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
  current_external_function_options.nargs = eExtFunSetDefaultNargs;
  current_external_function_options.firstDerivSymbID = eExtFunNotSet;
  current_external_function_options.secondDerivSymbID = eExtFunNotSet;
  current_external_function_id = eExtFunNotSet;
}

unique_ptr<ModFile>
ParsingDriver::parse(istream &in, bool debug)
{
  mod_file = make_unique<ModFile>(warnings);

  symbol_list.clear();

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

void
ParsingDriver::declare_symbol(const string &name, SymbolType type, const string &tex_name, const vector<pair<string, string>> &partition_value)
{
  try
    {
      mod_file->symbol_table.addSymbol(name, type, tex_name, partition_value);
    }
  catch (SymbolTable::AlreadyDeclaredException &e)
    {
      if (e.same_type)
        warning("Symbol " + name + " declared twice.");
      else
        error("Symbol " + name + " declared twice with different types!");
    }
}

void
ParsingDriver::declare_endogenous(const string &name, const string &tex_name, const vector<pair<string, string>> &partition_value)
{
  declare_symbol(name, SymbolType::endogenous, tex_name, partition_value);
}

void
ParsingDriver::declare_epilogue_endogenous(const string &name,
                                           const string &tex_name,
                                           const vector<pair<string, string>> &partition_value)
{
  declare_symbol(name, SymbolType::endogenousEpilogue, tex_name, partition_value);
}

void
ParsingDriver::declare_var_endogenous(const string &name)
{
  if (mod_file->symbol_table.exists(name))
    {
      SymbolType type = mod_file->symbol_table.getType(name);
      if (type != SymbolType::endogenous && type != SymbolType::exogenous && type != SymbolType::exogenousDet)
        error("Symbol " + name + " used in a VAR must be either endogenous or "
              +"exogenous if it is also used elsewhere in the .mod file");
      add_in_symbol_list(name);
      return;
    }

  declare_symbol(name, SymbolType::endogenousVAR, "", {});
  add_in_symbol_list(name);
}

void
ParsingDriver::declare_exogenous(const string &name, const string &tex_name, const vector<pair<string, string>> &partition_value)
{
  declare_symbol(name, SymbolType::exogenous, tex_name, partition_value);
}

void
ParsingDriver::declare_exogenous_det(const string &name, const string &tex_name, const vector<pair<string, string>> &partition_value)
{
  declare_symbol(name, SymbolType::exogenousDet, tex_name, partition_value);
}

void
ParsingDriver::declare_epilogue_exogenous(const string &name,
                                          const string &tex_name,
                                          const vector<pair<string, string>> &partition_value)
{
  declare_symbol(name, SymbolType::exogenousEpilogue, tex_name, partition_value);
}

void
ParsingDriver::declare_parameter(const string &name, const string &tex_name, const vector<pair<string, string>> &partition_value)
{
  declare_symbol(name, SymbolType::parameter, tex_name, partition_value);
}

void
ParsingDriver::declare_epilogue_parameter(const string &name,
                                          const string &tex_name,
                                          const vector<pair<string, string>> &partition_value)
{
  declare_symbol(name, SymbolType::parameterEpilogue, tex_name, partition_value);
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
ParsingDriver::declare_optimal_policy_discount_factor_parameter(expr_t exprnode)
{
  declare_parameter("optimal_policy_discount_factor");
  init_param("optimal_policy_discount_factor", exprnode);
}

void
ParsingDriver::begin_trend()
{
  set_current_data_tree(&mod_file->dynamic_model);
}

void
ParsingDriver::declare_trend_var(bool log_trend, const string &name, const string &tex_name)
{
  declare_symbol(name, log_trend ? SymbolType::logTrend : SymbolType::trend, tex_name, {});
  declared_trend_vars.push_back(mod_file->symbol_table.getID(name));
}

void
ParsingDriver::end_trend_var(expr_t growth_factor)
{
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
ParsingDriver::add_predetermined_variable(const string &name)
{
  check_symbol_is_endogenous(name);
  int symb_id = mod_file->symbol_table.getID(name);
  mod_file->symbol_table.markPredetermined(symb_id);
}

void
ParsingDriver::add_equation_tags(string key, string value)
{
  eq_tags.emplace_back(key, value);

  transform(key.begin(), key.end(), key.begin(), ::tolower);
  if (key.compare("endogenous") == 0)
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
  check_symbol_existence_in_model_block(name);
  int symb_id;
  try
    {
      symb_id = mod_file->symbol_table.getID(name);
    }
  catch (SymbolTable::UnknownSymbolNameException &e)
    {
      // Declare variable as exogenous to continue parsing
      // processing will end at end of model block if nostrict option was not passed
      declare_exogenous(name);
      undeclared_model_vars.insert(name);
      symb_id = mod_file->symbol_table.getID(name);
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

      // change in equations in ModelTree
      auto dm = make_unique<DynamicModel>(mod_file->symbol_table,
                                          mod_file->num_constants,
                                          mod_file->external_functions_table,
                                          mod_file->trend_component_model_table,
                                          mod_file->var_model_table);
      mod_file->dynamic_model.updateAfterVariableChange(*dm);

      // remove error messages
      undeclared_model_vars.erase(name);
      for (auto it = undeclared_model_variable_errors.begin();
           it != undeclared_model_variable_errors.end();)
        if (it->first == name)
          it = undeclared_model_variable_errors.erase(it);
        else
          it++;
    }
  catch (SymbolTable::UnknownSymbolNameException &e)
    {
      switch (new_type)
        {
        case SymbolType::endogenous:
          declare_endogenous(name);
          break;
        case SymbolType::exogenous:
          declare_exogenous(name);
          break;
        case SymbolType::parameter:
          declare_parameter(name);
          break;
        default:
          error("Type not yet supported");
        }
      symb_id = mod_file->symbol_table.getID(name);
    }
  return add_model_variable(symb_id, 0);

}

expr_t
ParsingDriver::add_model_variable(int symb_id, int lag)
{
  assert(symb_id >= 0);
  SymbolType type = mod_file->symbol_table.getType(symb_id);

  if (type == SymbolType::modFileLocalVariable)
    error("Variable " + mod_file->symbol_table.getName(symb_id) +
          " not allowed inside model declaration. Its scope is only outside model.");

  if (type == SymbolType::externalFunction)
    error("Symbol " + mod_file->symbol_table.getName(symb_id) +
          " is a function name external to Dynare. It cannot be used like a variable without input argument inside model.");

  if (type == SymbolType::modelLocalVariable && lag != 0)
    error("Model local variable " + mod_file->symbol_table.getName(symb_id) + " cannot be given a lead or a lag.");

  if (dynamic_cast<StaticModel *>(model_tree) != nullptr && lag != 0)
    error("Leads and lags on variables are forbidden in 'planner_objective'.");

  if (dynamic_cast<StaticModel *>(model_tree) != nullptr && type == SymbolType::modelLocalVariable)
    error("Model local variable " + mod_file->symbol_table.getName(symb_id) + " cannot be used in 'planner_objective'.");

  // It makes sense to allow a lead/lag on parameters: during steady state calibration, endogenous and parameters can be swapped
  return model_tree->AddVariable(symb_id, lag);
}

expr_t
ParsingDriver::add_expression_variable(const string &name)
{
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
ParsingDriver::declare_nonstationary_var(const string &name, const string &tex_name, const vector<pair<string, string>> &partition_value)
{
  declare_endogenous(name, tex_name, partition_value);

  declared_nonstationary_vars.push_back(mod_file->symbol_table.getID(name));
  mod_file->nonstationary_variables = true;
}

void
ParsingDriver::end_nonstationary_var(bool log_deflator, expr_t deflator)
{
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
ParsingDriver::begin_VAR_restrictions()
{
  clear_VAR_storage();
}

void
ParsingDriver::end_VAR_restrictions(const string &var_model_name)
{
  mod_file->addStatement(new VarRestrictionsStatement(var_model_name,
                                                      var_map,
                                                      exclusion_restrictions,
                                                      equation_restrictions,
                                                      crossequation_restrictions,
                                                      covariance_number_restriction,
                                                      covariance_pair_restriction,
                                                      mod_file->symbol_table));
  clear_VAR_storage();
}

void
ParsingDriver::clear_VAR_storage()
{
  exclusion_restriction.clear();
  exclusion_restrictions.clear();
  symbol_list.clear();
  var_restriction_eq_or_crosseq.clear();
  equation_restrictions.clear();
  crossequation_restrictions.clear();
  covariance_number_restriction.clear();
  covariance_pair_restriction.clear();
}

void
ParsingDriver::add_VAR_exclusion_restriction(const string &lagstr)
{
  int lag = stoi(lagstr);
  auto it = exclusion_restrictions.find(lag);
  if (it == exclusion_restrictions.end())
    exclusion_restrictions[lag] = exclusion_restriction;
  else
    for (map<int, SymbolList>::const_iterator it1 = exclusion_restriction.begin();
         it1 != exclusion_restriction.end(); it1++)
      it->second[it1->first] = it1->second;

  exclusion_restriction.clear();
}

void
ParsingDriver::add_VAR_restriction_coeff(const string &name1, const string &name2, const string &lagstr)
{
  int symb_id1 = mod_file->symbol_table.getID(name1);
  int symb_id2 = name2.empty() ? -1 : mod_file->symbol_table.getID(name2);
  int lag = stoi(lagstr);

  var_restriction_coeff = { symb_id1, { symb_id2, lag } };
}

void
ParsingDriver::add_VAR_restriction_eq_or_crosseq(expr_t expr)
{
  var_restriction_eq_or_crosseq.emplace_back(var_restriction_coeff, expr);
}

void
ParsingDriver::add_VAR_restriction_equation_or_crossequation(const string &numberstr)
{
  assert(var_restriction_eq_or_crosseq.size() > 0 && var_restriction_eq_or_crosseq.size() < 3);
  double number = stod(numberstr);
  if (var_restriction_eq_or_crosseq.size() == 1)
    var_restriction_equation_or_crossequation = { { var_restriction_eq_or_crosseq[0], { { -1, { -1, -1 } }, nullptr } }, number };
  else
    var_restriction_equation_or_crossequation = { { var_restriction_eq_or_crosseq[0], var_restriction_eq_or_crosseq[1] }, number };
  var_restriction_eq_or_crosseq.clear();
}

void
ParsingDriver::multiply_arg2_by_neg_one()
{
  assert(var_restriction_eq_or_crosseq.size() == 2);
  expr_t exprtm1 = add_times(var_restriction_eq_or_crosseq[1].second,
                             add_uminus(add_non_negative_constant("-1")));
  var_restriction_eq_or_crosseq[1] = { var_restriction_eq_or_crosseq[1].first, exprtm1 };
}

void
ParsingDriver::add_VAR_restriction_equation_or_crossequation_final(const string &name)
{
  if (!name.empty())
    {
      int symb_id = mod_file->symbol_table.getID(name);
      equation_restrictions[symb_id] = var_restriction_equation_or_crossequation;
    }
  else
    crossequation_restrictions.push_back(var_restriction_equation_or_crossequation);
}

void
ParsingDriver::add_VAR_restriction_exclusion_equation(const string &name)
{
  int symb_id = mod_file->symbol_table.getID(name);
  exclusion_restriction[symb_id] = symbol_list;
  symbol_list.clear();
}

void
ParsingDriver::add_VAR_covariance_number_restriction(const string &name1, const string &name2, const string &valuestr)
{
  int symb_id1 = mod_file->symbol_table.getID(name1);
  int symb_id2 = mod_file->symbol_table.getID(name2);
  double value = stod(valuestr);
  covariance_number_restriction[{ symb_id1, symb_id2 }] = value;
}

void
ParsingDriver::add_VAR_covariance_pair_restriction(const string &name11, const string &name12, const string &name21, const string &name22)
{
  int symb_id11 = mod_file->symbol_table.getID(name11);
  int symb_id12 = mod_file->symbol_table.getID(name12);
  int symb_id21 = mod_file->symbol_table.getID(name21);
  int symb_id22 = mod_file->symbol_table.getID(name22);
  covariance_pair_restriction[{ symb_id11, symb_id12 }] = { symb_id21, symb_id22 };
}

void
ParsingDriver::run_var_estimation()
{
  mod_file->addStatement(new VarEstimationStatement(options_list));
}

void
ParsingDriver::periods(const string &periods)
{
  warning("periods: this command is now deprecated and may be removed in a future version of Dynare. Please use the ''periods'' option of the ''simul'' command instead.");

  int periods_val = stoi(periods);
  mod_file->addStatement(new PeriodsStatement(periods_val));
}

void
ParsingDriver::dsample(const string &arg1)
{
  int arg1_val = stoi(arg1);
  mod_file->addStatement(new DsampleStatement(arg1_val));
}

void
ParsingDriver::dsample(const string &arg1, const string &arg2)
{
  int arg1_val = stoi(arg1);
  int arg2_val = stoi(arg2);
  mod_file->addStatement(new DsampleStatement(arg1_val, arg2_val));
}

void
ParsingDriver::init_param(const string &name, expr_t rhs)
{
  check_symbol_is_parameter(name);
  int symb_id = mod_file->symbol_table.getID(name);
  mod_file->addStatement(new InitParamStatement(symb_id, rhs, mod_file->symbol_table));
}

void
ParsingDriver::init_val(const string &name, expr_t rhs)
{
  if (nostrict)
    if (!mod_file->symbol_table.exists(name))
      {
        warning("discarding '" + name + "' as it was not recognized in the initval or endval statement");
        return;
      }

  check_symbol_is_endogenous_or_exogenous(name);
  int symb_id = mod_file->symbol_table.getID(name);
  init_values.emplace_back(symb_id, rhs);
}

void
ParsingDriver::initval_file(const string &filename)
{
  mod_file->addStatement(new InitvalFileStatement(filename));
}

void
ParsingDriver::hist_val(const string &name, const string &lag, expr_t rhs)
{
  if (nostrict)
    if (!mod_file->symbol_table.exists(name))
      {
        warning("discarding '" + name + "' as it was not recognized in the histavl block");
        return;
      }

  check_symbol_is_endogenous_or_exogenous(name);
  int symb_id = mod_file->symbol_table.getID(name);

  int ilag = stoi(lag);
  if (ilag > 0)
    error("histval: the lag on " + name + " should be less than or equal to 0");

  pair<int, int> key(symb_id, ilag);

  if (mod_file->dynamic_model.minLagForSymbol(symb_id) > ilag - 1)
    hist_vals_wrong_lag[symb_id] = ilag;

  if (hist_values.find(key) != hist_values.end())
    error("hist_val: (" + name + ", " + lag + ") declared twice");

  hist_values[key] = rhs;
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

  homotopy_values.emplace_back(symb_id, make_pair(val1, val2));
}

void
ParsingDriver::end_generate_irfs()
{
  mod_file->addStatement(new GenerateIRFsStatement(options_list, generate_irf_names, generate_irf_elements));

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
  check_symbol_is_exogenous(exo);
  if (generate_irf_exos.find(exo) != generate_irf_exos.end())
    error("You have set the exogenous variable " + exo + " twice.");

  generate_irf_exos[move(exo)] = stod(value);
}

void
ParsingDriver::forecast()
{
  mod_file->addStatement(new ForecastStatement(symbol_list, options_list));
  symbol_list.clear();
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
ParsingDriver::byte_code()
{
  mod_file->byte_code = true;
}

void
ParsingDriver::differentiate_forward_vars_all()
{
  mod_file->differentiate_forward_vars = true;
}

void
ParsingDriver::differentiate_forward_vars_some()
{
  mod_file->differentiate_forward_vars = true;
  mod_file->differentiate_forward_vars_subset = symbol_list.get_symbols();
  for (vector<string>::const_iterator it = mod_file->differentiate_forward_vars_subset.begin();
       it != mod_file->differentiate_forward_vars_subset.end(); ++it)
    check_symbol_is_endogenous(*it);
  symbol_list.clear();
}

void
ParsingDriver::cutoff(const string &value)
{
  double val = stod(value);
  mod_file->dynamic_model.cutoff = val;
  mod_file->static_model.cutoff = val;
}

void
ParsingDriver::mfs(const string &value)
{
  int val = stoi(value);
  mod_file->dynamic_model.mfs = val;
  mod_file->static_model.mfs = val;
}

void
ParsingDriver::end_initval(bool all_values_required)
{
  mod_file->addStatement(new InitValStatement(init_values, mod_file->symbol_table, all_values_required));
  init_values.clear();
}

void
ParsingDriver::end_endval(bool all_values_required)
{
  mod_file->addStatement(new EndValStatement(init_values, mod_file->symbol_table, all_values_required));
  init_values.clear();
}

void
ParsingDriver::end_histval(bool all_values_required)
{
  mod_file->addStatement(new HistValStatement(hist_values, hist_vals_wrong_lag, mod_file->symbol_table, all_values_required));
  hist_values.clear();
}

void
ParsingDriver::end_homotopy()
{
  mod_file->addStatement(new HomotopyStatement(homotopy_values, mod_file->symbol_table));
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
ParsingDriver::add_epilogue_equal(const string &varname, expr_t expr)
{
  int id;
  try
    {
      id = mod_file->symbol_table.getID(varname);
    }
  catch (SymbolTable::UnknownSymbolNameException &e)
    {
      error("Variable " + varname + " used in the epilogue block but was not declared.");
    }
  mod_file->epilogue.addDefinition(id, expr);
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
    for (vector<pair<string, string>>::const_iterator it = model_errors.begin();
         it != model_errors.end(); it++)
      {
        if (it->first == "")
          exit_after_write = true;
        cerr << it->second;
      }

  if (undeclared_model_variable_errors.size() > 0)
    for (vector<pair<string, string>>::const_iterator it = undeclared_model_variable_errors.begin();
         it != undeclared_model_variable_errors.end(); it++)
      if (nostrict)
        warning(it->second);
      else
        {
          exit_after_write = true;
          cerr << it->second << endl;
        }

  if (exit_after_write)
    exit(EXIT_FAILURE);

  reset_data_tree();
}

void
ParsingDriver::end_shocks(bool overwrite)
{
  mod_file->addStatement(new ShocksStatement(overwrite, det_shocks, var_shocks, std_shocks,
                                             covar_shocks, corr_shocks, mod_file->symbol_table));
  det_shocks.clear();
  var_shocks.clear();
  std_shocks.clear();
  covar_shocks.clear();
  corr_shocks.clear();
}

void
ParsingDriver::end_mshocks(bool overwrite)
{
  mod_file->addStatement(new MShocksStatement(overwrite, det_shocks, mod_file->symbol_table));
  det_shocks.clear();
}

void
ParsingDriver::add_det_shock(const string &var, bool conditional_forecast)
{
  if (conditional_forecast)
    check_symbol_is_endogenous(var);
  else
    check_symbol_is_exogenous(var);

  int symb_id = mod_file->symbol_table.getID(var);

  if (det_shocks.find(symb_id) != det_shocks.end())
    error("shocks/conditional_forecast_paths: variable " + var + " declared twice");

  if (det_shocks_periods.size() != det_shocks_values.size())
    error("shocks/conditional_forecast_paths: variable " + var + ": number of periods is different from number of shock values");

  vector<ShocksStatement::DetShockElement> v;

  for (size_t i = 0; i < det_shocks_periods.size(); i++)
    {
      ShocksStatement::DetShockElement dse;
      dse.period1 = det_shocks_periods[i].first;
      dse.period2 = det_shocks_periods[i].second;
      dse.value = det_shocks_values[i];
      v.push_back(dse);
    }

  det_shocks[symb_id] = v;

  det_shocks_periods.clear();
  det_shocks_values.clear();
}

void
ParsingDriver::add_stderr_shock(const string &var, expr_t value)
{
  if (nostrict)
    if (!mod_file->symbol_table.exists(var))
      {
        warning("discarding shocks block declaration of the standard error of '" + var + "' as it was not declared");
        return;
      }

  check_symbol_existence(var);
  int symb_id = mod_file->symbol_table.getID(var);

  if (var_shocks.find(symb_id) != var_shocks.end()
      || std_shocks.find(symb_id) != std_shocks.end())
    error("shocks: variance or stderr of shock on " + var + " declared twice");

  std_shocks[symb_id] = value;
}

void
ParsingDriver::add_var_shock(const string &var, expr_t value)
{
  if (nostrict)
    if (!mod_file->symbol_table.exists(var))
      {
        warning("discarding shocks block declaration of the variance of '" + var + "' as it was not declared");
        return;
      }

  check_symbol_existence(var);
  int symb_id = mod_file->symbol_table.getID(var);

  if (var_shocks.find(symb_id) != var_shocks.end()
      || std_shocks.find(symb_id) != std_shocks.end())
    error("shocks: variance or stderr of shock on " + var + " declared twice");

  var_shocks[symb_id] = value;
}

void
ParsingDriver::add_covar_shock(const string &var1, const string &var2, expr_t value)
{
  if (nostrict)
    if (!mod_file->symbol_table.exists(var1) || !mod_file->symbol_table.exists(var2))
      {
        warning("discarding shocks block declaration of the covariance of '" + var1 + "' and '" + var2 + "' as at least one was not declared");
        return;
      }

  check_symbol_existence(var1);
  check_symbol_existence(var2);
  int symb_id1 = mod_file->symbol_table.getID(var1);
  int symb_id2 = mod_file->symbol_table.getID(var2);

  pair<int, int> key(symb_id1, symb_id2), key_inv(symb_id2, symb_id1);

  if (covar_shocks.find(key) != covar_shocks.end()
      || covar_shocks.find(key_inv) != covar_shocks.end()
      || corr_shocks.find(key) != corr_shocks.end()
      || corr_shocks.find(key_inv) != corr_shocks.end())
    error("shocks: covariance or correlation shock on variable pair (" + var1 + ", "
          + var2 + ") declared twice");

  covar_shocks[key] = value;
}

void
ParsingDriver::add_correl_shock(const string &var1, const string &var2, expr_t value)
{
  if (nostrict)
    if (!mod_file->symbol_table.exists(var1) || !mod_file->symbol_table.exists(var2))
      {
        warning("discarding shocks block declaration of the correlation of '" + var1 + "' and '" + var2 + "' as at least one was not declared");
        return;
      }

  check_symbol_existence(var1);
  check_symbol_existence(var2);
  int symb_id1 = mod_file->symbol_table.getID(var1);
  int symb_id2 = mod_file->symbol_table.getID(var2);

  pair<int, int> key(symb_id1, symb_id2), key_inv(symb_id2, symb_id1);

  if (covar_shocks.find(key) != covar_shocks.end()
      || covar_shocks.find(key_inv) != covar_shocks.end()
      || corr_shocks.find(key) != corr_shocks.end()
      || corr_shocks.find(key_inv) != corr_shocks.end())
    error("shocks: covariance or correlation shock on variable pair (" + var1 + ", "
          + var2 + ") declared twice");

  corr_shocks[key] = value;
}

void
ParsingDriver::add_period(const string &p1, const string &p2)
{
  int p1_val = stoi(p1);
  int p2_val = stoi(p2);
  if (p1_val > p2_val)
    error("shocks/conditional_forecast_paths: can't have first period index greater than second index in range specification");
  det_shocks_periods.emplace_back(p1_val, p2_val);
}

void
ParsingDriver::add_period(const string &p1)
{
  int p1_val = stoi(p1);
  det_shocks_periods.emplace_back(p1_val, p1_val);
}

void
ParsingDriver::add_value(expr_t value)
{
  det_shocks_values.push_back(value);
}

void
ParsingDriver::add_value(const string &v)
{
  expr_t id;

  if (v.at(0) == '-')
    id = data_tree->AddUMinus(data_tree->AddNonNegativeConstant(v.substr(1, string::npos)));
  else
    id = data_tree->AddNonNegativeConstant(v);

  det_shocks_values.push_back(id);
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
  mod_file->addStatement(new SvarIdentificationStatement(svar_ident_restrictions,
                                                         svar_upper_cholesky,
                                                         svar_lower_cholesky,
                                                         svar_constants_exclusion,
                                                         mod_file->symbol_table));
  svar_restriction_symbols.clear();
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
  //    svar_ident_exclusion_values[make_pair(current_lag, it->first)] = it->second;

  svar_upper_cholesky = false;
  svar_lower_cholesky = false;
  svar_equation_restrictions.clear();
}

void
ParsingDriver::add_restriction_in_equation(const string &equation)
{
  int eqn = stoi(equation);
  if (eqn < 1)
    error("equation numbers must be greater than or equal to 1.");

  if (svar_equation_restrictions.count(eqn) > 0)
    error("equation number " + equation + " referenced more than once under a single lag.");

  svar_equation_restrictions[eqn] = svar_restriction_symbols;

  svar_restriction_symbols.clear();
}

void
ParsingDriver::add_in_svar_restriction_symbols(const string &tmp_var)
{
  check_symbol_existence(tmp_var);
  int symb_id = mod_file->symbol_table.getID(tmp_var);

  for (const auto &viit : svar_restriction_symbols)
    if (symb_id == viit)
      error(tmp_var + " restriction added twice.");

  svar_restriction_symbols.push_back(symb_id);
}

void
ParsingDriver::add_restriction_equation_nbr(const string &eq_nbr)
{
  svar_equation_nbr = stoi(eq_nbr);
  svar_left_handside = true;
  // reinitialize restriction type that must be set from the first restriction element
  svar_restriction_type = ParsingDriver::NOT_SET;
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
  if (svar_restriction_type == ParsingDriver::NOT_SET)
    {
      if (current_lag == 0)
        {
          svar_restriction_type = ParsingDriver::Qi_TYPE;
          ++svar_Qi_restriction_nbr[svar_equation_nbr];
        }
      else
        {
          svar_restriction_type = ParsingDriver::Ri_TYPE;
          ++svar_Ri_restriction_nbr[svar_equation_nbr];
        }
    }
  else
    {
      if ((svar_restriction_type == Qi_TYPE && current_lag > 0)
          || (svar_restriction_type == Ri_TYPE && current_lag == 0))
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
  if (value->eval(eval_context_t()) != 0)
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
  mod_file->addStatement(new SvarGlobalIdentificationCheckStatement);
}

void
ParsingDriver::do_sigma_e()
{
  warning("Sigma_e: this command is now deprecated and may be removed in a future version of Dynare. Please use the ''shocks'' command instead.");

  try
    {
      mod_file->addStatement(new SigmaeStatement(sigmae_matrix));
    }
  catch (SigmaeStatement::MatrixFormException &e)
    {
      error("Sigma_e: matrix is neither upper triangular nor lower triangular");
    }
  sigmae_matrix.clear();
}

void
ParsingDriver::end_of_row()
{
  sigmae_matrix.push_back(sigmae_row);
  sigmae_row.clear();
}

void
ParsingDriver::add_to_row_const(const string &v)
{
  expr_t id;

  if (v.at(0) == '-')
    id = data_tree->AddUMinus(data_tree->AddNonNegativeConstant(v.substr(1, string::npos)));
  else
    id = data_tree->AddNonNegativeConstant(v);

  sigmae_row.push_back(id);
}

void
ParsingDriver::add_to_row(expr_t v)
{
  sigmae_row.push_back(v);
}

void
ParsingDriver::steady()
{
  mod_file->addStatement(new SteadyStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::option_num(string name_option, string opt1, string opt2)
{
  if (options_list.paired_num_options.find(name_option)
      != options_list.paired_num_options.end())
    error("option " + name_option + " declared twice");

  options_list.paired_num_options[move(name_option)] = { move(opt1), move(opt2) };
}

void
ParsingDriver::option_num(string name_option, string opt)
{
  if (options_list.num_options.find(name_option) != options_list.num_options.end())
    error("option " + name_option + " declared twice");

  options_list.num_options[move(name_option)] = move(opt);
}

void
ParsingDriver::option_str(string name_option, string opt)
{
  if (options_list.string_options.find(name_option)
      != options_list.string_options.end())
    error("option " + name_option + " declared twice");

  options_list.string_options[move(name_option)] = move(opt);
}

void
ParsingDriver::option_date(string name_option, string opt)
{
  if (options_list.date_options.find(name_option)
      != options_list.date_options.end())
    error("option " + name_option + " declared twice");

  options_list.date_options[move(name_option)] = move(opt);
}

void
ParsingDriver::option_symbol_list(string name_option)
{
  if (options_list.symbol_list_options.find(name_option)
      != options_list.symbol_list_options.end())
    error("option " + name_option + " declared twice");

  if (name_option.compare("irf_shocks") == 0)
    {
      vector<string> shocks = symbol_list.get_symbols();
      for (vector<string>::const_iterator it = shocks.begin();
           it != shocks.end(); it++)
        if (mod_file->symbol_table.getType(*it) != SymbolType::exogenous)
          error("Variables passed to irf_shocks must be exogenous. Caused by: " + *it);
    }

  if (name_option.compare("ms.parameters") == 0)
    {
      vector<string> parameters = symbol_list.get_symbols();
      for (vector<string>::const_iterator it = parameters.begin();
           it != parameters.end(); it++)
        if (mod_file->symbol_table.getType(*it) != SymbolType::parameter)
          error("Variables passed to the parameters option of the markov_switching statement must be parameters. Caused by: " + *it);
    }

  options_list.symbol_list_options[move(name_option)] = symbol_list;
  symbol_list.clear();
}

void
ParsingDriver::option_vec_int(string name_option, vector<int> opt)
{
  if (options_list.vector_int_options.find(name_option)
      != options_list.vector_int_options.end())
    error("option " + name_option + " declared twice");

  if (opt.empty())
    error("option " + name_option + " was passed an empty vector.");

  options_list.vector_int_options[move(name_option)] = move(opt);
}

void
ParsingDriver::option_vec_str(string name_option, vector<string> opt)
{
  if (options_list.vector_str_options.find(name_option)
      != options_list.vector_str_options.end())
    error("option " + name_option + " declared twice");

  if (opt.empty())
    error("option " + name_option + " was passed an empty vector.");

  options_list.vector_str_options[move(name_option)] = move(opt);
}

void
ParsingDriver::linear()
{
  mod_file->linear = true;
}

void
ParsingDriver::add_in_symbol_list(const string &tmp_var)
{
  if (tmp_var != ":")
    check_symbol_existence(tmp_var);
  symbol_list.addSymbol(tmp_var);
}

void
ParsingDriver::rplot()
{
  mod_file->addStatement(new RplotStatement(symbol_list));
  symbol_list.clear();
}

void
ParsingDriver::stoch_simul()
{
  mod_file->addStatement(new StochSimulStatement(symbol_list, options_list));
  symbol_list.clear();
  options_list.clear();
}

void
ParsingDriver::trend_component_model()
{
  const auto its = options_list.string_options.find("trend_component.name");
  if (its == options_list.string_options.end())
    error("You must pass the model_name option to the trend_component_model statement.");
  auto name = its->second;

  const auto itvs = options_list.vector_str_options.find("trend_component.eqtags");
  if (itvs == options_list.vector_str_options.end())
    error("You must pass the eqtags option to the trend_component_model statement.");
  auto eqtags = itvs->second;

  const auto itvs1 = options_list.vector_str_options.find("trend_component.trends");
  if (itvs1 == options_list.vector_str_options.end())
    error("You must pass the trends option to the trend_component_model statement.");
  auto trends = itvs1->second;

  mod_file->trend_component_model_table.addTrendComponentModel(name, eqtags, trends);
  options_list.clear();
}

void
ParsingDriver::var_model()
{
  const auto its = options_list.string_options.find("var.model_name");
  if (its == options_list.string_options.end())
    error("You must pass the model_name option to the var_model statement.");
  auto name = its->second;

  int order = 0;
  const auto itn = options_list.num_options.find("var.order");
  if (itn != options_list.num_options.end())
    order = stoi(itn->second);
  else
  if (!symbol_list.empty())
    error("You must pass the order option when passing a symbol list to the var_model statement");

  vector<string> eqtags;
  const auto itvs = options_list.vector_str_options.find("var.eqtags");
  if (itvs != options_list.vector_str_options.end())
    {
      eqtags = itvs->second;
      if (!symbol_list.empty())
        error("You cannot pass a symbol list when passing equation tags to the var_model statement");
      else if (itn != options_list.num_options.end())
        error("You cannot pass the order option when passing equation tags to the var_model statement");
    }

  mod_file->var_model_table.addVarModel(name, eqtags, make_pair(symbol_list, order));
  symbol_list.clear();
  options_list.clear();
  var_map[its->second] = symbol_list.getSymbols();
}

void
ParsingDriver::simul()
{
  mod_file->addStatement(new SimulStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::model_info()
{
  mod_file->addStatement(new ModelInfoStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::check()
{
  mod_file->addStatement(new CheckStatement(options_list));
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
ParsingDriver::estimated_params()
{
  mod_file->addStatement(new EstimatedParamsStatement(estim_params_list, mod_file->symbol_table));
  estim_params_list.clear();
}

void
ParsingDriver::estimated_params_init(bool use_calibration)
{
  mod_file->addStatement(new EstimatedParamsInitStatement(estim_params_list, mod_file->symbol_table, use_calibration));
  estim_params_list.clear();
}

void
ParsingDriver::estimated_params_bounds()
{
  mod_file->addStatement(new EstimatedParamsBoundsStatement(estim_params_list, mod_file->symbol_table));
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
  mod_file->addStatement(new OsrParamsBoundsStatement(osr_params_list));
  osr_params_list.clear();
}

void
ParsingDriver::set_unit_root_vars()
{
  mod_file->addStatement(new UnitRootVarsStatement());
  warning("''unit_root_vars'' is now obsolete; use the ''diffuse_filter'' option of ''estimation'' instead");
  symbol_list.clear();
}

void
ParsingDriver::set_time(const string &arg)
{
  option_date("initial_period", arg);
  mod_file->addStatement(new SetTimeStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::estimation_data()
{
  mod_file->addStatement(new EstimationDataStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::set_subsamples(string name1, string name2)
{
  check_symbol_existence(name1);
  if (!name2.empty())
    check_symbol_existence(name2);

  mod_file->addStatement(new SubsamplesStatement(name1, name2, subsample_declaration_map,
                                                 mod_file->symbol_table));
  subsample_declarations[{ move(name1), move(name2) }] = subsample_declaration_map;
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

  if (subsample_declarations.find({ from_name1, from_name2 }) == subsample_declarations.end())
    {
      string err{from_name1};
      if (!from_name2.empty())
        err.append(",").append(from_name2);
      error(err + " does not have an associated subsample statement.");
    }

  mod_file->addStatement(new SubsamplesEqualStatement(to_name1, to_name2, from_name1, from_name2,
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
  if (subsample_declaration_map.find(name) != subsample_declaration_map.end())
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
  if (tmp_map.find(subsample_name) == tmp_map.end())
    error("The subsample name " + subsample_name + " was not previously declared in a subsample statement.");
}

void
ParsingDriver::set_prior(const string &name, const string &subsample_name)
{
  check_symbol_is_parameter(name);
  check_subsample_declaration_exists(name, subsample_name);
  mod_file->addStatement(new PriorStatement(name, subsample_name, prior_shape, prior_variance, options_list));
  options_list.clear();
  set_prior_variance();
  prior_shape = PriorDistributions::noShape;
}

void
ParsingDriver::set_joint_prior(const vector<string> &symbol_vec)
{
  for (auto &it : symbol_vec)
    add_joint_parameter(it);
  mod_file->addStatement(new JointPriorStatement(joint_parameters, prior_shape, options_list));
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
ParsingDriver::copy_prior(const string &to_declaration_type, const string &to_name1,
                          const string &to_name2, const string &to_subsample_name,
                          const string &from_declaration_type, const string &from_name1,
                          const string &from_name2, const string &from_subsample_name)
{
  if (to_declaration_type == "par")
    check_symbol_is_parameter(to_name1);
  else
    {
      check_symbol_is_endogenous_or_exogenous(to_name1);
      if (!to_name2.empty())
        check_symbol_is_endogenous_or_exogenous(to_name2);
    }

  if (from_declaration_type == "par")
    check_symbol_is_parameter(from_name1);
  else
    {
      check_symbol_is_endogenous_or_exogenous(from_name1);
      if (!from_name2.empty())
        check_symbol_is_endogenous_or_exogenous(from_name2);
    }

  mod_file->addStatement(new PriorEqualStatement(to_declaration_type, to_name1,
                                                 to_name2, to_subsample_name,
                                                 from_declaration_type, from_name1,
                                                 from_name2, from_subsample_name,
                                                 mod_file->symbol_table));

}

void
ParsingDriver::set_options(const string &name, const string &subsample_name)
{
  check_symbol_is_parameter(name);
  check_subsample_declaration_exists(name, subsample_name);
  mod_file->addStatement(new OptionsStatement(name, subsample_name, options_list));
  options_list.clear();
}

void
ParsingDriver::copy_options(const string &to_declaration_type, const string &to_name1,
                            const string &to_name2, const string &to_subsample_name,
                            const string &from_declaration_type, const string &from_name1,
                            const string &from_name2, const string &from_subsample_name)
{
  if (to_declaration_type == "par")
    check_symbol_is_parameter(to_name1);
  else
    {
      check_symbol_is_endogenous_or_exogenous(to_name1);
      if (!to_name2.empty())
        check_symbol_is_endogenous_or_exogenous(to_name2);
    }

  if (from_declaration_type == "par")
    check_symbol_is_parameter(from_name1);
  else
    {
      check_symbol_is_endogenous_or_exogenous(from_name1);
      if (!from_name2.empty())
        check_symbol_is_endogenous_or_exogenous(from_name2);
    }

  mod_file->addStatement(new OptionsEqualStatement(to_declaration_type, to_name1,
                                                   to_name2, to_subsample_name,
                                                   from_declaration_type, from_name1,
                                                   from_name2, from_subsample_name,
                                                   mod_file->symbol_table));
}

void
ParsingDriver::check_symbol_is_endogenous_or_exogenous(const string &name)
{
  check_symbol_existence(name);
  switch (mod_file->symbol_table.getType(name))
    {
    case SymbolType::endogenous:
    case SymbolType::exogenous:
    case SymbolType::exogenousDet:
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
ParsingDriver::check_symbol_is_exogenous(const string &name)
{
  check_symbol_existence(name);
  switch (mod_file->symbol_table.getType(name))
    {
    case SymbolType::exogenous:
    case SymbolType::exogenousDet:
      break;
    default:
      error(name + " is not exogenous.");
    }
}

void
ParsingDriver::set_std_prior(const string &name, const string &subsample_name)
{
  check_symbol_is_endogenous_or_exogenous(name);
  check_subsample_declaration_exists(name, subsample_name);
  mod_file->addStatement(new StdPriorStatement(name, subsample_name, prior_shape, prior_variance,
                                               options_list, mod_file->symbol_table));
  options_list.clear();
  set_prior_variance();
  prior_shape = PriorDistributions::noShape;
}

void
ParsingDriver::set_std_options(const string &name, const string &subsample_name)
{
  check_symbol_is_endogenous_or_exogenous(name);
  check_subsample_declaration_exists(name, subsample_name);
  mod_file->addStatement(new StdOptionsStatement(name, subsample_name, options_list, mod_file->symbol_table));
  options_list.clear();
}

void
ParsingDriver::set_corr_prior(const string &name1, const string &name2, const string &subsample_name)
{
  check_symbol_is_endogenous_or_exogenous(name1);
  check_symbol_is_endogenous_or_exogenous(name2);
  check_subsample_declaration_exists(name1, name2, subsample_name);
  mod_file->addStatement(new CorrPriorStatement(name1, name2, subsample_name, prior_shape, prior_variance,
                                                options_list, mod_file->symbol_table));
  options_list.clear();
  set_prior_variance();
  prior_shape = PriorDistributions::noShape;
}

void
ParsingDriver::set_corr_options(const string &name1, const string &name2, const string &subsample_name)
{
  check_symbol_is_endogenous_or_exogenous(name1);
  check_symbol_is_endogenous_or_exogenous(name2);
  check_subsample_declaration_exists(name1, name2, subsample_name);
  mod_file->addStatement(new CorrOptionsStatement(name1, name2, subsample_name, options_list, mod_file->symbol_table));
  options_list.clear();
}

void
ParsingDriver::run_estimation()
{
  mod_file->addStatement(new EstimationStatement(symbol_list, options_list));
  symbol_list.clear();
  options_list.clear();
}

void
ParsingDriver::dynare_sensitivity()
{
  mod_file->addStatement(new DynareSensitivityStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::optim_options_helper(const string &name)
{
  if (options_list.string_options.find("optim_opt") == options_list.string_options.end())
    options_list.string_options["optim_opt"] = "";
  else
    options_list.string_options["optim_opt"] += ",";
  options_list.string_options["optim_opt"] += "''" + name + "'',";
}

void
ParsingDriver::optim_options_string(const string &name, const string &value)
{
  optim_options_helper(name);
  options_list.string_options["optim_opt"] += "''" + value + "''";
}

void
ParsingDriver::optim_options_num(const string &name, const string &value)
{
  optim_options_helper(name);
  options_list.string_options["optim_opt"] += value;
}

void
ParsingDriver::sampling_options_helper(const string &name)
{
  if (options_list.string_options.find("posterior_sampler_options.sampling_opt") ==
      options_list.string_options.end())
    options_list.string_options["posterior_sampler_options.sampling_opt"] = "";
  else
    options_list.string_options["posterior_sampler_options.sampling_opt"] += ",";
  options_list.string_options["posterior_sampler_options.sampling_opt"] += "''" + name + "'',";
}

void
ParsingDriver::sampling_options_string(const string &name, const string &value)
{
  sampling_options_helper(name);
  options_list.string_options["posterior_sampler_options.sampling_opt"] += "''" + value + "''";
}

void
ParsingDriver::sampling_options_num(const string &name, const string &value)
{
  sampling_options_helper(name);
  options_list.string_options["posterior_sampler_options.sampling_opt"] += value;
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
  mod_file->addStatement(new ObservationTrendsStatement(trend_elements, mod_file->symbol_table));
  trend_elements.clear();
}

void
ParsingDriver::set_trend_element(string arg1, expr_t arg2)
{
  check_symbol_existence(arg1);
  if (trend_elements.find(arg1) != trend_elements.end())
    error("observation_trends: " + arg1 + " declared twice");
  trend_elements[move(arg1)] = arg2;
}

void
ParsingDriver::set_optim_weights(string name, expr_t value)
{
  check_symbol_is_endogenous(name);
  if (var_weights.find(name) != var_weights.end())
    error("optim_weights: " + name + " declared twice");
  var_weights[move(name)] = move(value);
}

void
ParsingDriver::set_optim_weights(const string &name1, const string &name2, expr_t value)
{
  check_symbol_is_endogenous(name1);
  check_symbol_is_endogenous(name2);

  pair<string, string> covar_key{name1, name2};

  if (covar_weights.find(covar_key) != covar_weights.end())
    error("optim_weights: pair of variables (" + name1 + ", " + name2
          + ") declared twice");

  covar_weights[covar_key] = value;
}

void
ParsingDriver::optim_weights()
{
  mod_file->addStatement(new OptimWeightsStatement(var_weights, covar_weights, mod_file->symbol_table));
  var_weights.clear();
  covar_weights.clear();
}

void
ParsingDriver::set_osr_params()
{
  mod_file->addStatement(new OsrParamsStatement(symbol_list, mod_file->symbol_table));
  symbol_list.clear();
}

void
ParsingDriver::run_osr()
{
  mod_file->addStatement(new OsrStatement(symbol_list, options_list));
  symbol_list.clear();
  options_list.clear();
}

void
ParsingDriver::run_dynatype(const string &filename)
{
  mod_file->addStatement(new DynaTypeStatement(symbol_list, filename));
  symbol_list.clear();
}

void
ParsingDriver::run_dynasave(const string &filename)
{
  mod_file->addStatement(new DynaSaveStatement(symbol_list, filename));
  symbol_list.clear();
}

void
ParsingDriver::run_load_params_and_steady_state(const string &filename)
{
  mod_file->addStatement(new LoadParamsAndSteadyStateStatement(filename, mod_file->symbol_table, warnings));
}

void
ParsingDriver::run_save_params_and_steady_state(const string &filename)
{
  mod_file->addStatement(new SaveParamsAndSteadyStateStatement(filename));
}

void
ParsingDriver::run_identification()
{
  mod_file->addStatement(new IdentificationStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::add_mc_filename(string filename, string prior)
{
  for (auto & it : filename_list)
    if (it.first == filename)
      error("model_comparison: filename " + filename + " declared twice");
  filename_list.emplace_back(move(filename), move(prior));
}

void
ParsingDriver::run_model_comparison()
{
  mod_file->addStatement(new ModelComparisonStatement(filename_list, options_list));
  filename_list.clear();
  options_list.clear();
}

void
ParsingDriver::begin_planner_objective()
{
  planner_objective_statement = new PlannerObjectiveStatement(mod_file->symbol_table,
                                                              mod_file->num_constants,
                                                              mod_file->external_functions_table,
                                                              mod_file->trend_component_model_table,
                                                              mod_file->var_model_table);
  set_current_data_tree(&planner_objective_statement->getPlannerObjective());
}

void
ParsingDriver::end_planner_objective(expr_t expr)
{
  // Add equation corresponding to expression
  expr_t eq = model_tree->AddEqual(expr, model_tree->Zero);
  model_tree->addEquation(eq, location.begin.line);

  mod_file->addStatement(planner_objective_statement);
  planner_objective_statement = nullptr; // The object is now owned by mod_file

  reset_data_tree();
}

void
ParsingDriver::ramsey_model()
{
  if (!mod_file->symbol_table.exists("optimal_policy_discount_factor"))
    declare_optimal_policy_discount_factor_parameter(data_tree->One);
  mod_file->addStatement(new RamseyModelStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::ramsey_policy()
{
  if (!mod_file->symbol_table.exists("optimal_policy_discount_factor"))
    declare_optimal_policy_discount_factor_parameter(data_tree->One);
  mod_file->addStatement(new RamseyPolicyStatement(mod_file->symbol_table, ramsey_policy_list, options_list));
  options_list.clear();
  ramsey_policy_list.clear();
}

void
ParsingDriver::add_to_ramsey_policy_list(string name)
{
  ramsey_policy_list.push_back(move(name));
}

void
ParsingDriver::discretionary_policy()
{
  if (!mod_file->symbol_table.exists("optimal_policy_discount_factor"))
    declare_optimal_policy_discount_factor_parameter(data_tree->One);
  mod_file->addStatement(new DiscretionaryPolicyStatement(symbol_list, options_list));
  symbol_list.clear();
  options_list.clear();
}

void
ParsingDriver::write_latex_dynamic_model(bool write_equation_tags)
{
  mod_file->addStatement(new WriteLatexDynamicModelStatement(mod_file->dynamic_model, write_equation_tags));
}

void
ParsingDriver::write_latex_static_model(bool write_equation_tags)
{
  mod_file->addStatement(new WriteLatexStaticModelStatement(mod_file->static_model, write_equation_tags));
}

void
ParsingDriver::write_latex_original_model(bool write_equation_tags)
{
  mod_file->addStatement(new WriteLatexOriginalModelStatement(mod_file->original_model, write_equation_tags));
}

void
ParsingDriver::write_latex_steady_state_model()
{
  mod_file->addStatement(new WriteLatexSteadyStateModelStatement(mod_file->steady_state_model));
}

void
ParsingDriver::bvar_density(const string &maxnlags)
{
  mod_file->addStatement(new BVARDensityStatement(stoi(maxnlags), options_list));
  options_list.clear();
}

void
ParsingDriver::bvar_forecast(const string &nlags)
{
  mod_file->addStatement(new BVARForecastStatement(stoi(nlags), options_list));
  options_list.clear();
}

void
ParsingDriver::sbvar()
{
  mod_file->addStatement(new SBVARStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::ms_estimation()
{
  mod_file->addStatement(new MSSBVAREstimationStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::ms_simulation()
{
  mod_file->addStatement(new MSSBVARSimulationStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::ms_compute_mdd()
{
  mod_file->addStatement(new MSSBVARComputeMDDStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::ms_compute_probabilities()
{
  mod_file->addStatement(new MSSBVARComputeProbabilitiesStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::ms_irf()
{
  mod_file->addStatement(new MSSBVARIrfStatement(symbol_list, options_list));
  symbol_list.clear();
  options_list.clear();
}

void
ParsingDriver::ms_forecast()
{
  mod_file->addStatement(new MSSBVARForecastStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::ms_variance_decomposition()
{
  mod_file->addStatement(new MSSBVARVarianceDecompositionStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::svar()
{
  OptionsList::string_options_t::const_iterator it0, it1, it2;
  OptionsList::num_options_t::const_iterator itn;
  OptionsList::vec_int_options_t::const_iterator itv;

  it0 = options_list.string_options.find("ms.coefficients");
  it1 = options_list.string_options.find("ms.variances");
  it2 = options_list.string_options.find("ms.constants");
  if (it0 == options_list.string_options.end()
      && it1 == options_list.string_options.end()
      && it2 == options_list.string_options.end())
    error("You must pass one of 'coefficients', 'variances', or 'constants'.");

  if ((it0 != options_list.string_options.end()
       && it1 != options_list.string_options.end())
      || (it1 != options_list.string_options.end()
          && it2 != options_list.string_options.end())
      || (it0 != options_list.string_options.end()
          && it2 != options_list.string_options.end()))
    error("You may only pass one of 'coefficients', 'variances', or 'constants'.");

  itn = options_list.num_options.find("ms.chain");
  if (itn == options_list.num_options.end())
    error("A chain option must be passed to the svar statement.");
  else if (stoi(itn->second) <= 0)
    error("The value passed to the chain option must be greater than zero.");

  itv = options_list.vector_int_options.find("ms.equations");
  if (itv != options_list.vector_int_options.end())
    for (int viit : itv->second)
      if (viit <= 0)
        error("The value(s) passed to the equation option must be greater than zero.");

  mod_file->addStatement(new SvarStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::markov_switching()
{
  OptionsList::num_options_t::const_iterator it0;

  it0 = options_list.num_options.find("ms.chain");
  if (it0 == options_list.num_options.end())
    error("A chain option must be passed to the markov_switching statement.");
  else if (stoi(it0->second) <= 0)
    error("The value passed to the chain option must be greater than zero.");

  it0 = options_list.num_options.find("ms.number_of_regimes");
  if (it0 == options_list.num_options.end())
    error("A number_of_regimes option must be passed to the markov_switching statement.");
  else if (stoi(it0->second) <= 0)
    error("The value passed to the number_of_regimes option must be greater than zero.");

  it0 = options_list.num_options.find("ms.duration");
  if (it0 == options_list.num_options.end())
    error("A duration option must be passed to the markov_switching statement.");

  mod_file->addStatement(new MarkovSwitchingStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::shock_decomposition()
{
  mod_file->addStatement(new ShockDecompositionStatement(symbol_list, options_list));
  symbol_list.clear();
  options_list.clear();
}

void
ParsingDriver::realtime_shock_decomposition()
{
  mod_file->addStatement(new RealtimeShockDecompositionStatement(symbol_list, options_list));
  symbol_list.clear();
  options_list.clear();
}

void
ParsingDriver::plot_shock_decomposition()
{
  mod_file->addStatement(new PlotShockDecompositionStatement(symbol_list, options_list));
  symbol_list.clear();
  options_list.clear();
}

void
ParsingDriver::initial_condition_decomposition()
{
  mod_file->addStatement(new InitialConditionDecompositionStatement(symbol_list, options_list));
  symbol_list.clear();
  options_list.clear();
}

void
ParsingDriver::conditional_forecast()
{
  mod_file->addStatement(new ConditionalForecastStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::plot_conditional_forecast(const string &periods)
{
  int nperiods;
  if (periods.empty())
    nperiods = -1;
  else
    nperiods = stoi(periods);

  mod_file->addStatement(new PlotConditionalForecastStatement(nperiods, symbol_list));
  symbol_list.clear();
}

void
ParsingDriver::conditional_forecast_paths()
{
  mod_file->addStatement(new ConditionalForecastPathsStatement(det_shocks, mod_file->symbol_table));
  det_shocks.clear();
}

void
ParsingDriver::calib_smoother()
{
  mod_file->addStatement(new CalibSmootherStatement(symbol_list, options_list));
  symbol_list.clear();
  options_list.clear();
}

void
ParsingDriver::extended_path()
{
  mod_file->addStatement(new ExtendedPathStatement(options_list));
  options_list.clear();
}

expr_t
ParsingDriver::add_model_equal(expr_t arg1, expr_t arg2)
{
  expr_t id = model_tree->AddEqual(arg1, arg2);

  // Detect if the equation is tagged [static]
  bool is_static_only = false;
  for (vector<pair<string, string>>::const_iterator it = eq_tags.begin();
       it != eq_tags.end(); ++it)
    if (it->first == "static")
      {
        is_static_only = true;
        break;
      }

  if (is_static_only)
    {
      if (!id->isInStaticForm())
        error("An equation tagged [static] cannot contain leads, lags, expectations or STEADY_STATE operators");

      dynamic_model->addStaticOnlyEquation(id, location.begin.line, eq_tags);
    }
  else
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
ParsingDriver::declare_model_local_variable(const string &name, const string &tex_name)
{
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
      // It can have already been declared in a steady_state_model block, check that it is indeed a ModelLocalVariable
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
ParsingDriver::change_type(SymbolType new_type, const vector<string> &var_list)
{
  for (auto & it : var_list)
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
  return data_tree->AddExpectation(stoi(arg1), arg2);
}

expr_t
ParsingDriver::add_var_expectation(const string &model_name)
{
  return data_tree->AddVarExpectation(model_name);
}

expr_t
ParsingDriver::add_pac_expectation(const string &var_model_name)
{
  return data_tree->AddPacExpectation(var_model_name);
}

void
ParsingDriver::pac_model()
{
  OptionsList::string_options_t::const_iterator it = options_list.string_options.find("pac.model_name");
  if (it == options_list.string_options.end())
    error("You must pass the model_name option to the pac_model statement.");
  auto name = it->second;

  it = options_list.string_options.find("pac.aux_model_name");
  if (it == options_list.string_options.end())
    error("You must pass the auxiliary_model_name option to the pac_model statement.");
  auto aux_model_name = it->second;

  it = options_list.string_options.find("pac.discount");
  if (it == options_list.string_options.end())
    error("You must pass the discount option to the pac_model statement.");
  auto discount = it->second;

  string growth;
  it = options_list.string_options.find("pac.growth");
  if (it != options_list.string_options.end())
    growth = it->second;

  mod_file->addStatement(new PacModelStatement(name, aux_model_name, discount, growth,
                                               mod_file->symbol_table));

  symbol_list.clear();
  options_list.clear();
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
  vector<int> lags;
  for (int i = 1; i <= stoi(lag); i++)
    lags.push_back(i);

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
ParsingDriver::add_steady_state(expr_t arg1)
{
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
        current_external_function_options.firstDerivSymbID = eExtFunSetButNoNameProvided;
      else
        {
          declare_symbol(opt, SymbolType::externalFunction, "", {});
          current_external_function_options.firstDerivSymbID = mod_file->symbol_table.getID(opt);
        }
    }
  else if (name_option == "second_deriv_provided")
    {
      if (opt.empty())
        current_external_function_options.secondDerivSymbID = eExtFunSetButNoNameProvided;
      else
        {
          declare_symbol(opt, SymbolType::externalFunction, "", {});
          current_external_function_options.secondDerivSymbID = mod_file->symbol_table.getID(opt);
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
  if (current_external_function_id == eExtFunNotSet)
    error("The 'name' option must be passed to external_function().");

  if (current_external_function_options.secondDerivSymbID >= 0
      && current_external_function_options.firstDerivSymbID  == eExtFunNotSet)
    error("If the second derivative is provided to the external_function command, the first derivative must also be provided.");

  if (current_external_function_options.secondDerivSymbID == eExtFunSetButNoNameProvided
      && current_external_function_options.firstDerivSymbID  != eExtFunSetButNoNameProvided)
    error("If the second derivative is provided in the top-level function, the first derivative must also be provided in that function.");

  mod_file->external_functions_table.addExternalFunction(current_external_function_id, current_external_function_options, true);
  reset_current_external_function_options();
}

void
ParsingDriver::push_external_function_arg_vector_onto_stack()
{
  vector<expr_t> emptyvec;
  stack_external_function_args.push(emptyvec);
}

void
ParsingDriver::add_external_function_arg(expr_t arg)
{
  stack_external_function_args.top().push_back(arg);
}

pair<bool, double>
ParsingDriver::is_there_one_integer_argument() const
{
  if (stack_external_function_args.top().size() != 1)
    return { false, 0 };

  auto *numNode = dynamic_cast<NumConstNode *>(stack_external_function_args.top().front());
  auto *unaryNode = dynamic_cast<UnaryOpNode *>(stack_external_function_args.top().front());

  if (numNode == nullptr && unaryNode == nullptr)
    return { false, 0 };

  eval_context_t ectmp;
  double model_var_arg;
  if (unaryNode == nullptr)
    {
      try
        {
          model_var_arg = numNode->eval(ectmp);
        }
      catch (ExprNode::EvalException &e)
        {
          return { false, 0 };
        }
    }
  else
    if (unaryNode->get_op_code() != UnaryOpcode::uminus)
      return { false, 0 };
    else
      {
        try
          {
            model_var_arg = unaryNode->eval(ectmp);
          }
        catch (ExprNode::EvalException &e)
          {
            return { false, 0 };
          }
      }

  if (model_var_arg != floor(model_var_arg))
    return { false, 0 };
  return { true, model_var_arg };
}

expr_t
ParsingDriver::add_model_var_or_external_function(const string &function_name, bool in_model_block)
{
  expr_t nid;
  if (mod_file->symbol_table.exists(function_name))
    if (mod_file->symbol_table.getType(function_name) != SymbolType::externalFunction)
      if (!in_model_block && !parsing_epilogue)
        {
          if (stack_external_function_args.top().size() > 0)
            error(string("Symbol ") + function_name + string(" cannot take arguments."));
          else
            return add_expression_variable(function_name);
        }
      else
        { // e.g. model_var(lag) => ADD MODEL VARIABLE WITH LEAD (NumConstNode)/LAG (UnaryOpNode)
          if (undeclared_model_vars.find(function_name) != undeclared_model_vars.end())
            undeclared_model_variable_error("Unknown symbol: " + function_name, function_name);

          pair<bool, double> rv = is_there_one_integer_argument();
          if (!rv.first)
            model_error("Symbol " + function_name +
                        " is being treated as if it were a function (i.e., takes an argument that is not an integer).", "");

          nid = add_model_variable(mod_file->symbol_table.getID(function_name), (int) rv.second);
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
          if (mod_file->external_functions_table.getNargs(symb_id) == eExtFunNotSet)
            error("Before using " + function_name
                  +"() in the model block, you must first declare it via the external_function() statement");
          else if ((int) (stack_external_function_args.top().size()) != mod_file->external_functions_table.getNargs(symb_id))
            error("The number of arguments passed to " + function_name
                  +"() does not match those of a previous call or declaration of this function.");
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

          pair<bool, double> rv = is_there_one_integer_argument();
          if (rv.first)
            {
              // assume it's a lead/lagged variable
              declare_exogenous(function_name);
              return add_model_variable(mod_file->symbol_table.getID(function_name), (int) rv.second);
            }
          else
            error("To use an external function (" + function_name +
                  ") within the model block, you must first declare it via the external_function() statement.");
        }
      declare_symbol(function_name, SymbolType::externalFunction, "", {});
      current_external_function_options.nargs = stack_external_function_args.top().size();
      mod_file->external_functions_table.addExternalFunction(mod_file->symbol_table.getID(function_name),
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
ParsingDriver::add_native(const string &s)
{
  mod_file->addStatement(new NativeStatement(s));
}

void
ParsingDriver::add_native_remove_charset(string str, const string &token)
{
  size_t found = str.find(token);

  assert(found != string::npos);
  str.resize(found);
  add_native(str);
}

void
ParsingDriver::add_verbatim(const string &s)
{
  mod_file->addStatement(new VerbatimStatement(s));
}

void
ParsingDriver::add_verbatim_remove_charset(string str, const string &token)
{
  size_t found = str.find(token);

  assert(found != string::npos);
  str.resize(found);
  add_verbatim(str);
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

  SymbolType type = mod_file->symbol_table.getType(id);
  if (type != SymbolType::endogenous && type != SymbolType::modFileLocalVariable && type != SymbolType::parameter)
    error(varname + " has incorrect type");

  mod_file->steady_state_model.addDefinition(id, expr);
}

void
ParsingDriver::add_steady_state_model_equal_multiple(expr_t expr)
{
  const vector<string> &symbs = symbol_list.get_symbols();
  vector<int> ids;

  for (const auto & symb : symbs)
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
      SymbolType type = mod_file->symbol_table.getType(id);
      if (type != SymbolType::endogenous && type != SymbolType::modFileLocalVariable && type != SymbolType::parameter)
        error(symb + " has incorrect type");
      ids.push_back(id);
    }

  mod_file->steady_state_model.addMultipleDefinitions(ids, expr);

  symbol_list.clear();
}

void
ParsingDriver::add_graph_format(const string &name)
{
  graph_formats.addSymbol(name);
}

void
ParsingDriver::process_graph_format_option()
{
  options_list.symbol_list_options["graph_format"] = graph_formats;
  graph_formats.clear();
}

void
ParsingDriver::plot_shock_decomp_process_graph_format_option()
{
  options_list.symbol_list_options["plot_shock_decomp.graph_format"] = graph_formats;
  graph_formats.clear();
}

void
ParsingDriver::model_diagnostics()
{
  mod_file->addStatement(new ModelDiagnosticsStatement());
}

void
ParsingDriver::add_parallel_local_file(string filename)
{
  mod_file->parallel_local_files.push_back(move(filename));
}

void
ParsingDriver::add_moment_calibration_item(const string &endo1, const string &endo2, string lags, const pair<string, string> &range)
{
  MomentCalibration::Constraint c;

  check_symbol_is_endogenous(endo1);
  c.endo1 = mod_file->symbol_table.getID(endo1);

  check_symbol_is_endogenous(endo2);
  c.endo2 = mod_file->symbol_table.getID(endo2);

  c.lags = move(lags);

  c.lower_bound = range.first;
  c.upper_bound = range.second;

  moment_calibration_constraints.push_back(c);
}

void
ParsingDriver::end_moment_calibration()
{
  mod_file->addStatement(new MomentCalibration(moment_calibration_constraints,
                                               mod_file->symbol_table));
  moment_calibration_constraints.clear();
}

void
ParsingDriver::add_irf_calibration_item(const string &endo, string periods, const string &exo, const pair<string, string> &range)
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

  irf_calibration_constraints.push_back(c);
}

void
ParsingDriver::end_irf_calibration()
{
  mod_file->addStatement(new IrfCalibration(irf_calibration_constraints,
                                            mod_file->symbol_table,
                                            options_list));
  irf_calibration_constraints.clear();
}

void
ParsingDriver::smoother2histval()
{
  mod_file->addStatement(new Smoother2histvalStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::histval_file(const string &filename)
{
  mod_file->addStatement(new HistvalFileStatement(filename));
}

void
ParsingDriver::perfect_foresight_setup()
{
  mod_file->addStatement(new PerfectForesightSetupStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::perfect_foresight_solver()
{
  mod_file->addStatement(new PerfectForesightSolverStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::gmm_estimation()
{
  mod_file->addStatement(new GMMEstimationStatement(symbol_list, options_list));
  symbol_list.clear();
  options_list.clear();
}

void
ParsingDriver::smm_estimation()
{
  mod_file->addStatement(new SMMEstimationStatement(symbol_list, options_list));
  symbol_list.clear();
  options_list.clear();
}

void
ParsingDriver::prior_posterior_function(bool prior_func)
{
  mod_file->addStatement(new PriorPosteriorFunctionStatement((bool) prior_func, options_list));
  options_list.clear();
}

void
ParsingDriver::add_ramsey_constraints_statement()
{
  mod_file->addStatement(new RamseyConstraintsStatement(mod_file->symbol_table, ramsey_constraints));
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
  SymbolType type = mod_file->symbol_table.getType(symb_id);

  if (type != SymbolType::exogenous)
    error("shock_groups: " + name + " should be an exogenous variable");

  shock_group.push_back(move(name));
}

void
ParsingDriver::add_shock_group(string name)
{
  ShockGroupsStatement::Group G;
  G.name = move(name);
  G.list = shock_group;
  shock_groups.push_back(G);

  shock_group.clear();
}

void
ParsingDriver::end_shock_groups(const string &name)
{
  mod_file->addStatement(new ShockGroupsStatement(shock_groups, name));
  shock_groups.clear();
}

void
ParsingDriver::var_expectation_model()
{
  auto it = options_list.string_options.find("variable");
  if (it == options_list.string_options.end())
    error("You must pass the variable option to the var_expectation_model statement.");
  auto variable = it->second;
  check_symbol_is_endogenous(variable);

  it = options_list.string_options.find("auxiliary_model_name");
  if (it == options_list.string_options.end())
    error("You must pass the auxiliary_model_name option to the var_expectation_model statement.");
  auto var_model_name = it->second;

  it = options_list.string_options.find("model_name");
  if (it == options_list.string_options.end())
    error("You must pass the model_name option to the var_expectation_model statement.");
  auto model_name = it->second;

  it = options_list.num_options.find("horizon");
  if (it == options_list.num_options.end())
    error("You must pass the horizon option to the var_expectation_model statement.");
  auto horizon = it->second;

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

  mod_file->addStatement(new VarExpectationModelStatement(model_name, variable, var_model_name, horizon,
                                                          var_expectation_model_discount, mod_file->symbol_table));

  options_list.clear();
  var_expectation_model_discount = nullptr;
}
