/*
 * Copyright (C) 2003-2009 Dynare Team
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

#include <cstdlib>
#include <fstream>
#include <iostream>

#include "ParsingDriver.hh"
#include "Statement.hh"

bool
ParsingDriver::symbol_exists_and_is_not_modfile_local_or_unknown_function(const char *s)
{
  if (!mod_file->symbol_table.exists(s))
    return false;

  SymbolType type = mod_file->symbol_table.getType(s);

  return(type != eModFileLocalVariable && type != eUnknownFunction);
}

void
ParsingDriver::check_symbol_existence(const string &name)
{
  if (!mod_file->symbol_table.exists(name))
    error("Unknown symbol: " + name);
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

ModFile *
ParsingDriver::parse(istream &in, bool debug)
{
  mod_file = new ModFile();

  symbol_list.clear();

  reset_data_tree();
  estim_params.init(*data_tree);

  lexer = new DynareFlex(&in);
  lexer->set_debug(debug);

  Dynare::parser parser(*this);
  parser.set_debug_level(debug);
  parser.parse();

  delete lexer;

  return mod_file;
}

void
ParsingDriver::error(const Dynare::parser::location_type &l, const string &m)
{
  cerr << "ERROR: " << l << ": " << m << endl;
  exit(EXIT_FAILURE);
}

void
ParsingDriver::error(const string &m)
{
  error(location, m);
}

void
ParsingDriver::warning(const string &m)
{
  cerr << "WARNING: " << location << ": " << m << endl;
}

void
ParsingDriver::declare_symbol(string *name, SymbolType type, string *tex_name)
{
  try
    {
      if (tex_name == NULL)
        mod_file->symbol_table.addSymbol(*name, type);
      else
        mod_file->symbol_table.addSymbol(*name, type, *tex_name);
    }
  catch(SymbolTable::AlreadyDeclaredException &e)
    {
      if (e.same_type)
        warning("Symbol " + *name + " declared twice.");
      else
        error("Symbol " + *name + " declared twice with different types!");
    }

  delete name;
  if (tex_name != NULL)
    delete tex_name;
}

void
ParsingDriver::declare_endogenous(string *name, string *tex_name)
{
  declare_symbol(name, eEndogenous, tex_name);
}

void
ParsingDriver::declare_exogenous(string *name, string *tex_name)
{
  declare_symbol(name, eExogenous, tex_name);
}

void
ParsingDriver::declare_exogenous_det(string *name, string *tex_name)
{
  declare_symbol(name, eExogenousDet, tex_name);
}

void
ParsingDriver::declare_parameter(string *name, string *tex_name)
{
  declare_symbol(name, eParameter, tex_name);
}

void
ParsingDriver::add_predetermined_variable(string *name)
{
  try
    {
      int symb_id = mod_file->symbol_table.getID(*name);
      if (mod_file->symbol_table.getType(symb_id) != eEndogenous)
        error("Predetermined variables must be endogenous variables");

      mod_file->symbol_table.markPredetermined(symb_id);
    }
  catch(SymbolTable::UnknownSymbolNameException &e)
    {
      error("Undeclared symbol name: " + *name);
    }
  delete name;
}

void
ParsingDriver::add_equation_tags(string *key, string *value)
{
  int n = model_tree->equation_number();
  model_tree->addEquationTags(n, *key, *value);
  delete key;
  delete value;
}

NodeID
ParsingDriver::add_constant(string *constant)
{
  NodeID id = data_tree->AddNumConstant(*constant);
  delete constant;
  return id;
}

NodeID
ParsingDriver::add_nan_constant()
{
  return data_tree->NaN;
}

NodeID
ParsingDriver::add_inf_constant()
{
  return data_tree->Infinity;
}

NodeID
ParsingDriver::add_model_variable(string *name)
{
  return add_model_variable(name, new string("0"));
}

NodeID
ParsingDriver::add_model_variable(string *name, string *olag)
{
  check_symbol_existence(*name);
  SymbolType type = mod_file->symbol_table.getType(*name);
  int lag = atoi(olag->c_str());

  if (type == eModFileLocalVariable)
    error("Variable " + *name + " not allowed inside model declaration. Its scope is only outside model.");

  if (type == eUnknownFunction)
    error("Symbol " + *name + " is a function name unknown to Dynare. It cannot be used inside model.");

  if (type == eModelLocalVariable && lag != 0)
    error("Model local variable " + *name + " cannot be given a lead or a lag.");

  // It makes sense to allow a lead/lag on parameters: during steady state calibration, endogenous and parameters can be swapped

  NodeID id = model_tree->AddVariable(*name, lag);

  delete name;
  delete olag;
  return id;
}

NodeID
ParsingDriver::add_expression_variable(string *name)
{
  // If symbol doesn't exist, then declare it as a mod file local variable
  if (!mod_file->symbol_table.exists(*name))
    mod_file->symbol_table.addSymbol(*name, eModFileLocalVariable);

  // This check must come after the previous one!
  if (mod_file->symbol_table.getType(*name) == eModelLocalVariable)
    error("Variable " + *name + " not allowed outside model declaration. Its scope is only inside model.");

  NodeID id = data_tree->AddVariable(*name);

  delete name;
  return id;
}

void
ParsingDriver::periods(string *periods)
{
  warning("periods: this command is now deprecated and may be removed in a future version of Dynare. Please of the \"periods\" option of \"simul\" command instead.");

  int periods_val = atoi(periods->c_str());
  mod_file->addStatement(new PeriodsStatement(periods_val));
  delete periods;
}

void
ParsingDriver::dsample(string *arg1)
{
  int arg1_val = atoi(arg1->c_str());
  mod_file->addStatement(new DsampleStatement(arg1_val));
  delete arg1;
}

void
ParsingDriver::dsample(string *arg1, string *arg2)
{
  int arg1_val = atoi(arg1->c_str());
  int arg2_val = atoi(arg2->c_str());
  mod_file->addStatement(new DsampleStatement(arg1_val, arg2_val));
  delete arg1;
  delete arg2;
}


void
ParsingDriver::init_param(string *name, NodeID rhs)
{
  check_symbol_existence(*name);
  int symb_id = mod_file->symbol_table.getID(*name);
  if (mod_file->symbol_table.getType(symb_id) != eParameter)
    error(*name + " is not a parameter");

  mod_file->addStatement(new InitParamStatement(symb_id, rhs, mod_file->symbol_table));

  delete name;
}

void
ParsingDriver::init_val(string *name, NodeID rhs)
{
  check_symbol_existence(*name);
  int symb_id = mod_file->symbol_table.getID(*name);
  SymbolType type = mod_file->symbol_table.getType(symb_id);

  if (type != eEndogenous
      && type != eExogenous
      && type != eExogenousDet)
    error("initval/endval: " + *name + " should be an endogenous or exogenous variable");

  init_values.push_back(make_pair(symb_id, rhs));

  delete name;
}

void
ParsingDriver::initval_file(string *filename)
{
  mod_file->addStatement(new InitvalFileStatement(*filename));
  delete filename;
}

void
ParsingDriver::hist_val(string *name, string *lag, NodeID rhs)
{
  check_symbol_existence(*name);
  int symb_id = mod_file->symbol_table.getID(*name);
  SymbolType type = mod_file->symbol_table.getType(symb_id);

  if (type != eEndogenous
      && type != eExogenous
      && type != eExogenousDet)
    error("hist_val: " + *name + " should be an endogenous or exogenous variable");

  int ilag = atoi(lag->c_str());
  pair<int, int> key(symb_id, ilag);

  if (hist_values.find(key) != hist_values.end())
    error("hist_val: (" + *name + ", " + *lag + ") declared twice");

  hist_values[key] = rhs;

  delete name;
  delete lag;
}

void
ParsingDriver::homotopy_val(string *name, NodeID val1, NodeID val2)
{
  check_symbol_existence(*name);
  int symb_id = mod_file->symbol_table.getID(*name);
  SymbolType type = mod_file->symbol_table.getType(symb_id);

  if (type != eParameter
      && type != eExogenous
      && type != eExogenousDet)
    error("homotopy_val: " + *name + " should be a parameter or exogenous variable");

  homotopy_values.push_back(make_pair(symb_id, make_pair(val1, val2)));

  delete name;
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
ParsingDriver::byte_code()
{
  mod_file->byte_code = true;
}

void
ParsingDriver::cutoff(string *value)
{
  double val = atof(value->c_str());
  mod_file->dynamic_model.cutoff = val;
  mod_file->static_dll_model.cutoff = val;
  delete value;
}

void
ParsingDriver::mfs(string *value)
{
  int val = atoi(value->c_str());
  mod_file->dynamic_model.mfs = val;
  mod_file->static_dll_model.mfs = val;
  delete value;
}

void
ParsingDriver::end_initval()
{
  mod_file->addStatement(new InitValStatement(init_values, mod_file->symbol_table));
  init_values.clear();
}

void
ParsingDriver::end_endval()
{
  mod_file->addStatement(new EndValStatement(init_values, mod_file->symbol_table));
  init_values.clear();
}

void
ParsingDriver::end_histval()
{
  mod_file->addStatement(new HistValStatement(hist_values, mod_file->symbol_table));
  hist_values.clear();
}

void
ParsingDriver::end_homotopy()
{
  mod_file->addStatement(new HomotopyStatement(homotopy_values, mod_file->symbol_table));
  homotopy_values.clear();
}

void
ParsingDriver::begin_model()
{
  set_current_data_tree(&mod_file->dynamic_model);
}

void
ParsingDriver::end_shocks()
{
  mod_file->addStatement(new ShocksStatement(det_shocks, var_shocks, std_shocks,
                                             covar_shocks, corr_shocks, mod_file->symbol_table));
  det_shocks.clear();
  var_shocks.clear();
  std_shocks.clear();
  covar_shocks.clear();
  corr_shocks.clear();
}

void
ParsingDriver::end_mshocks()
{
  mod_file->addStatement(new MShocksStatement(det_shocks, mod_file->symbol_table));
  det_shocks.clear();
}

void
ParsingDriver::add_det_shock(string *var, bool conditional_forecast)
{
  check_symbol_existence(*var);
  SymbolType type = mod_file->symbol_table.getType(*var);

  if (conditional_forecast)
    {
      if (type != eEndogenous)
        error("conditional_forecast_paths: shocks can only be applied to exogenous variables");
    }
  else
    {
      if (type != eExogenous && type != eExogenousDet)
        error("shocks: shocks can only be applied to exogenous variables");
    }

  if (det_shocks.find(*var) != det_shocks.end())
    error("shocks/conditional_forecast_paths: variable " + *var + " declared twice");

  if (det_shocks_periods.size() != det_shocks_values.size())
    error("shocks/conditional_forecast_paths: variable " + *var + ": number of periods is different from number of shock values");

  vector<ShocksStatement::DetShockElement> v;

  for(unsigned int i = 0; i < det_shocks_periods.size(); i++)
    {
      ShocksStatement::DetShockElement dse;
      dse.period1 = det_shocks_periods[i].first;
      dse.period2 = det_shocks_periods[i].second;
      dse.value = det_shocks_values[i];
      v.push_back(dse);
    }

  det_shocks[*var] = v;

  det_shocks_periods.clear();
  det_shocks_values.clear();
  delete var;
}

void
ParsingDriver::add_stderr_shock(string *var, NodeID value)
{
  check_symbol_existence(*var);
  if (var_shocks.find(*var) != var_shocks.end()
      || std_shocks.find(*var) != std_shocks.end())
    error("shocks: variance or stderr of shock on " + *var + " declared twice");

  std_shocks[*var] = value;

  delete var;
}

void
ParsingDriver::add_var_shock(string *var, NodeID value)
{
  check_symbol_existence(*var);
  if (var_shocks.find(*var) != var_shocks.end()
      || std_shocks.find(*var) != std_shocks.end())
    error("shocks: variance or stderr of shock on " + *var + " declared twice");

  var_shocks[*var] = value;

  delete var;
}

void
ParsingDriver::add_covar_shock(string *var1, string *var2, NodeID value)
{
  check_symbol_existence(*var1);
  check_symbol_existence(*var2);

  pair<string, string> key(*var1, *var2), key_inv(*var2, *var1);

  if (covar_shocks.find(key) != covar_shocks.end()
      || covar_shocks.find(key_inv) != covar_shocks.end()
      || corr_shocks.find(key) != corr_shocks.end()
      || corr_shocks.find(key_inv) != corr_shocks.end())
    error("shocks: covariance or correlation shock on variable pair (" + *var1 + ", "
          + *var2 + ") declared twice");

  covar_shocks[key] = value;

  delete var1;
  delete var2;
}

void
ParsingDriver::add_correl_shock(string *var1, string *var2, NodeID value)
{
  check_symbol_existence(*var1);
  check_symbol_existence(*var2);

  pair<string, string> key(*var1, *var2), key_inv(*var2, *var1);

  if (covar_shocks.find(key) != covar_shocks.end()
      || covar_shocks.find(key_inv) != covar_shocks.end()
      || corr_shocks.find(key) != corr_shocks.end()
      || corr_shocks.find(key_inv) != corr_shocks.end())
    error("shocks: covariance or correlation shock on variable pair (" + *var1 + ", "
          + *var2 + ") declared twice");

  corr_shocks[key] = value;

  delete var1;
  delete var2;
}

void
ParsingDriver::add_period(string *p1, string *p2)
{
  int p1_val = atoi(p1->c_str());
  int p2_val = atoi(p2->c_str());
  if (p1_val > p2_val)
    error("shocks/conditional_forecast_paths: can't have first period index greater than second index in range specification");
  det_shocks_periods.push_back(make_pair(p1_val, p2_val));
  delete p1;
  delete p2;
}

void
ParsingDriver::add_period(string *p1)
{
  int p1_val = atoi(p1->c_str());
  det_shocks_periods.push_back(make_pair(p1_val, p1_val));
  delete p1;
}

void
ParsingDriver::add_value(NodeID value)
{
  det_shocks_values.push_back(value);
}

void
ParsingDriver::add_value(string *p1)
{
  det_shocks_values.push_back(add_constant(p1));
}

void
ParsingDriver::end_svar_identification()
{
  mod_file->addStatement(new SvarIdentificationStatement(svar_ident_exclusion_values,
                                                         svar_upper_cholesky,
                                                         svar_lower_cholesky,
                                                         mod_file->symbol_table));
  svar_upper_cholesky = false;
  svar_lower_cholesky = false;
  svar_restriction_symbols.clear();
  svar_equation_restrictions.clear();
  svar_ident_exclusion_values.clear();
}

void
ParsingDriver::combine_lag_and_restriction(string *lag)
{
  int current_lag = atoi(lag->c_str());

  for (SvarIdentificationStatement::svar_identification_exclusion_type::const_iterator it = svar_ident_exclusion_values.begin();
       it != svar_ident_exclusion_values.end(); it++)
    if (it->first.first == current_lag)
      error("lag " + *lag + " used more than once.");

  for(map<int, vector<string> >::const_iterator it = svar_equation_restrictions.begin();
      it != svar_equation_restrictions.end(); it++ )
    svar_ident_exclusion_values[make_pair(current_lag, it->first)] = it->second;

  svar_upper_cholesky = false;
  svar_lower_cholesky = false;
  svar_equation_restrictions.clear();
  delete lag;
}

void
ParsingDriver::add_restriction_in_equation(string *equation)
{
  int eqn = atoi(equation->c_str());
  if (eqn < 1)
    error("equation numbers must be greater than or equal to 1.");

  if (svar_equation_restrictions.count(eqn) > 0)
    error("equation number " + *equation + " referenced more than once under a single lag.");

  svar_equation_restrictions[eqn] = svar_restriction_symbols;

  svar_restriction_symbols.clear();
  delete equation;
}

void
ParsingDriver::add_in_svar_restriction_symbols(string *tmp_var)
{
  check_symbol_existence(*tmp_var);
  for (unsigned int i=0; i<svar_restriction_symbols.size(); i++)
    if (*tmp_var==svar_restriction_symbols.at(i))
      error(*tmp_var + " restriction added twice.");

  svar_restriction_symbols.push_back(*tmp_var);
  delete tmp_var;
}

void
ParsingDriver::add_upper_cholesky()
{
  svar_upper_cholesky = true;
  svar_lower_cholesky = false;
}

void
ParsingDriver::add_lower_cholesky()
{
  svar_upper_cholesky = false;
  svar_lower_cholesky = true;
}

void
ParsingDriver::do_sigma_e()
{
  warning("Sigma_e: this command is now deprecated and may be removed in a future version of Dynare. Please use the \"shocks\" command instead.");

  try
    {
      mod_file->addStatement(new SigmaeStatement(sigmae_matrix));
    }
  catch(SigmaeStatement::MatrixFormException &e)
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
ParsingDriver::add_to_row_const(string *s)
{
  sigmae_row.push_back(add_constant(s));
}

void
ParsingDriver::add_to_row(NodeID v)
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
ParsingDriver::option_num(const string &name_option, string *opt1, string *opt2)
{
  if (options_list.paired_num_options.find(name_option)
      != options_list.paired_num_options.end())
    error("option " + name_option + " declared twice");

  options_list.paired_num_options[name_option] = make_pair(*opt1, *opt2);
  delete opt1;
  delete opt2;
}

void
ParsingDriver::option_num(const string &name_option, string *opt)
{
  option_num(name_option, *opt);
  delete opt;
}

void
ParsingDriver::option_num(const string &name_option, const string &opt)
{
  if (options_list.num_options.find(name_option) != options_list.num_options.end())
    error("option " + name_option + " declared twice");

  if ((name_option == "periods") && mod_file->block)
    mod_file->dynamic_model.block_triangular.periods = atoi(opt.c_str());

  options_list.num_options[name_option] = opt;
}

void
ParsingDriver::option_str(const string &name_option, string *opt)
{
  option_str(name_option, *opt);
  delete opt;
}

void
ParsingDriver::option_str(const string &name_option, const string &opt)
{
  if (options_list.string_options.find(name_option)
      != options_list.string_options.end())
    error("option " + name_option + " declared twice");

  options_list.string_options[name_option] = opt;
}

void
ParsingDriver::option_symbol_list(const string &name_option)
{
  if (options_list.symbol_list_options.find(name_option)
      != options_list.symbol_list_options.end())
    error("option " + name_option + " declared twice");

  options_list.symbol_list_options[name_option] = symbol_list;
  symbol_list.clear();
}

void
ParsingDriver::linear()
{
  mod_file->linear = true;
}

void
ParsingDriver::add_in_symbol_list(string *tmp_var)
{
  if (*tmp_var != ":")
    check_symbol_existence(*tmp_var);
  symbol_list.addSymbol(*tmp_var);
  delete tmp_var;
}

void ParsingDriver::rplot()
{
  mod_file->addStatement(new RplotStatement(symbol_list, options_list));
  options_list.clear();
  symbol_list.clear();
}

void ParsingDriver::stoch_simul()
{
  mod_file->addStatement(new StochSimulStatement(symbol_list, options_list));
  symbol_list.clear();
  options_list.clear();
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
  check_symbol_existence(estim_params.name);
  if (estim_params.name2.size() > 0)
    check_symbol_existence(estim_params.name2);

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
ParsingDriver::estimated_params_init()
{
  mod_file->addStatement(new EstimatedParamsInitStatement(estim_params_list, mod_file->symbol_table));
  estim_params_list.clear();
}

void
ParsingDriver::estimated_params_bounds()
{
  mod_file->addStatement(new EstimatedParamsBoundsStatement(estim_params_list, mod_file->symbol_table));
  estim_params_list.clear();
}

void
ParsingDriver::set_unit_root_vars()
{
  mod_file->addStatement(new UnitRootVarsStatement(symbol_list));
  symbol_list.clear();
}

void
ParsingDriver::run_estimation()
{
  mod_file->addStatement(new EstimationStatement(symbol_list, options_list));
  symbol_list.clear();
  options_list.clear();
}

void
ParsingDriver::run_prior_analysis()
{
  mod_file->addStatement(new PriorAnalysisStatement(symbol_list, options_list));
  symbol_list.clear();
  options_list.clear();
}

void
ParsingDriver::run_posterior_analysis()
{
  mod_file->addStatement(new PosteriorAnalysisStatement(symbol_list, options_list));
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
ParsingDriver::optim_options_string(string *name, string *value)
{
  optim_options_helper(*name);
  options_list.string_options["optim_opt"] += "''" + *value + "''";
  delete name;
  delete value;
}

void
ParsingDriver::optim_options_num(string *name, string *value)
{
  optim_options_helper(*name);
  options_list.string_options["optim_opt"] += *value;
  delete name;
  delete value;
}

void
ParsingDriver::set_varobs()
{
  mod_file->addStatement(new VarobsStatement(symbol_list));
  symbol_list.clear();
}

void
ParsingDriver::set_trends()
{
  mod_file->addStatement(new ObservationTrendsStatement(trend_elements, mod_file->symbol_table));
  trend_elements.clear();
}

void
ParsingDriver::set_trend_element(string *arg1, NodeID arg2)
{
  check_symbol_existence(*arg1);
  if (trend_elements.find(*arg1) != trend_elements.end())
    error("observation_trends: " + *arg1 + " declared twice");
  trend_elements[*arg1] = arg2;
  delete arg1;
}

void
ParsingDriver::set_optim_weights(string *name, NodeID value)
{
  check_symbol_existence(*name);
  if (mod_file->symbol_table.getType(*name) != eEndogenous)
    error("optim_weights: " + *name + " isn't an endogenous variable");
  if (var_weights.find(*name) != var_weights.end())
    error("optim_weights: " + *name + " declared twice");
  var_weights[*name] = value;
  delete name;
}

void
ParsingDriver::set_optim_weights(string *name1, string *name2, NodeID value)
{
  check_symbol_existence(*name1);
  if (mod_file->symbol_table.getType(*name1) != eEndogenous)
    error("optim_weights: " + *name1 + " isn't an endogenous variable");

  check_symbol_existence(*name2);
  if (mod_file->symbol_table.getType(*name2) != eEndogenous)
    error("optim_weights: " + *name2 + " isn't an endogenous variable");

  pair<string, string> covar_key(*name1, *name2);

  if (covar_weights.find(covar_key) != covar_weights.end())
    error("optim_weights: pair of variables (" + *name1 + ", " + *name2
          + ") declared twice");

  covar_weights[covar_key] = value;
  delete name1;
  delete name2;
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
  mod_file->addStatement(new OsrParamsStatement(symbol_list));
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
ParsingDriver::set_calib_var(string *name, string *weight, NodeID expression)
{
  check_symbol_existence(*name);
  if (mod_file->symbol_table.getType(*name) != eEndogenous
      && mod_file->symbol_table.getType(*name) != eExogenous)
    error("calib_var: " + *name + " isn't an endogenous or exogenous variable");

  if (calib_var.find(*name) != calib_var.end())
    error("calib_var: " + *name + " declared twice");

  calib_var[*name] = make_pair(*weight, expression);

  delete name;
  delete weight;
}

void
ParsingDriver::set_calib_covar(string *name1, string *name2,
                               string *weight, NodeID expression)
{
  check_symbol_existence(*name1);
  check_symbol_existence(*name2);
  if (mod_file->symbol_table.getType(*name1) != mod_file->symbol_table.getType(*name2))
    error("calib_var: " + *name1 + " and " + *name2 + "dont't have the same type");
  if (mod_file->symbol_table.getType(*name1) != eEndogenous
      && mod_file->symbol_table.getType(*name1) != eExogenous)
    error("calib_var: " + *name1 + " and " + *name2 + "aren't endogenous or exogenous variables");

  pair<string, string> covar_key(*name1, *name2);

  if (calib_covar.find(covar_key) != calib_covar.end())
    error("calib_var: pair of variables (" + *name1 + ", " + *name2
          + ") declared twice");

  calib_covar[covar_key] = make_pair(*weight, expression);

  delete name1;
  delete name2;
  delete weight;
}

void
ParsingDriver::set_calib_ac(string *name, string *ar,
                            string *weight, NodeID expression)
{
  check_symbol_existence(*name);
  if (mod_file->symbol_table.getType(*name) != eEndogenous)
    error("calib_var: " + *name + "isn't an endogenous variable");

  int iar = atoi(ar->c_str());
  pair<string, int> ac_key(*name, iar);

  if (calib_ac.find(ac_key) != calib_ac.end())
    error("calib_var: autocorr " + *name + "(" + *ar + ") declared twice");

  calib_ac[ac_key] = make_pair(*weight, expression);

  delete name;
  delete ar;
  delete weight;
}

void
ParsingDriver::run_calib_var()
{
  mod_file->addStatement(new CalibVarStatement(calib_var, calib_covar, calib_ac,
                                               mod_file->symbol_table));
  calib_var.clear();
  calib_covar.clear();
  calib_ac.clear();
}

void
ParsingDriver::run_calib(int covar)
{
  mod_file->addStatement(new CalibStatement(covar));
}

void
ParsingDriver::run_dynatype(string *filename)
{
  mod_file->addStatement(new DynaTypeStatement(symbol_list, *filename));
  symbol_list.clear();
  delete filename;
}

void
ParsingDriver::run_dynasave(string *filename)
{
  mod_file->addStatement(new DynaSaveStatement(symbol_list, *filename));
  symbol_list.clear();
  delete filename;
}

void
ParsingDriver::run_load_params_and_steady_state(string *filename)
{
  mod_file->addStatement(new LoadParamsAndSteadyStateStatement(*filename, mod_file->symbol_table));
  delete filename;
}

void
ParsingDriver::run_save_params_and_steady_state(string *filename)
{
  mod_file->addStatement(new SaveParamsAndSteadyStateStatement(*filename));
  delete filename;
}

void
ParsingDriver::run_identification()
{
  mod_file->addStatement(new IdentificationStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::add_mc_filename(string *filename, string *prior)
{
  for(ModelComparisonStatement::filename_list_type::iterator it = filename_list.begin();
      it != filename_list.end(); it++)
    if ((*it).first == *filename)
      error("model_comparison: filename " + *filename + " declared twice");
  filename_list.push_back(make_pair(*filename, *prior));
  delete filename;
  delete prior;
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
  set_current_data_tree(new StaticModel(mod_file->symbol_table, mod_file->num_constants));
}

void
ParsingDriver::end_planner_objective(NodeID expr)
{
  // Add equation corresponding to expression
  NodeID eq = model_tree->AddEqual(expr, model_tree->Zero);
  model_tree->addEquation(eq);

  mod_file->addStatement(new PlannerObjectiveStatement(dynamic_cast<StaticModel *>(model_tree)));

  reset_data_tree();
}

void
ParsingDriver::ramsey_policy()
{
  mod_file->addStatement(new RamseyPolicyStatement(symbol_list, options_list));
  symbol_list.clear();
  options_list.clear();
}

void
ParsingDriver::write_latex_dynamic_model()
{
  mod_file->addStatement(new WriteLatexDynamicModelStatement(mod_file->dynamic_model));
}

void
ParsingDriver::write_latex_static_model()
{
  mod_file->addStatement(new WriteLatexStaticModelStatement(mod_file->static_model));
}

void
ParsingDriver::bvar_density(string *maxnlags)
{
  mod_file->addStatement(new BVARDensityStatement(atoi(maxnlags->c_str()), options_list));
  options_list.clear();
  delete maxnlags;
}

void
ParsingDriver::bvar_forecast(string *nlags)
{
  mod_file->addStatement(new BVARForecastStatement(atoi(nlags->c_str()), options_list));
  options_list.clear();
  delete nlags;
}

void
ParsingDriver::sbvar()
{
  mod_file->addStatement(new SBVARStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::ms_sbvar()
{
  mod_file->addStatement(new MS_SBVARStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::markov_switching()
{
  OptionsList::num_options_type::const_iterator it0, it1;

  it0 = options_list.num_options.find("ms.chain");
  if (it0 == options_list.num_options.end())
    error("A chain option must be passed to the markov_switching statement.");
  else if (atoi(it0->second.c_str()) <= 0)
    error("The value passed to the chain option must be greater than zero.");

  it0=options_list.num_options.find("ms.state");
  it1=options_list.num_options.find("ms.number_of_states");
  if ((it0 == options_list.num_options.end()) &&
      (it1 == options_list.num_options.end()))
    error("Either a state option or a number_of_states option must be passed to the markov_switching statement.");

  if ((it0 != options_list.num_options.end()) &&
      (it1 != options_list.num_options.end()))
    error("You cannot pass both a state option and a number_of_states option to the markov_switching statement.");

  if (it0 != options_list.num_options.end())
    if (atoi(it0->second.c_str()) <= 0)
      error("The value passed to the state option must be greater than zero.");

  if (it1 != options_list.num_options.end())
    if (atoi(it1->second.c_str()) <= 0)
      error("The value passed to the number_of_states option must be greater than zero.");

  string infStr ("Inf");
  it0 = options_list.num_options.find("ms.duration");
  if (it0 == options_list.num_options.end())
    error("A duration option must be passed to the markov_switching statement.");
  else if (infStr.compare(it0->second) != 0)
    if (atof(it0->second.c_str()) <= 0.0)
      error("The value passed to the duration option must be greater than zero.");

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
ParsingDriver::conditional_forecast()
{
  mod_file->addStatement(new ConditionalForecastStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::plot_conditional_forecast(string *periods)
{
  int nperiods;
  if (periods == NULL)
    nperiods = -1;
  else
    {
      nperiods = atoi(periods->c_str());
      delete periods;
    }
  mod_file->addStatement(new PlotConditionalForecastStatement(nperiods, symbol_list));
  symbol_list.clear();
}

void
ParsingDriver::conditional_forecast_paths()
{
  mod_file->addStatement(new ConditionalForecastPathsStatement(det_shocks, mod_file->symbol_table));
  det_shocks.clear();
}

NodeID
ParsingDriver::add_model_equal(NodeID arg1, NodeID arg2)
{
  NodeID id = model_tree->AddEqual(arg1, arg2);
  model_tree->addEquation(id);
  return id;
}

NodeID
ParsingDriver::add_model_equal_with_zero_rhs(NodeID arg)
{
  return add_model_equal(arg, model_tree->Zero);
}

void
ParsingDriver::declare_and_init_model_local_variable(string *name, NodeID rhs)
{
  try
    {
      mod_file->symbol_table.addSymbol(*name, eModelLocalVariable);
    }
  catch(SymbolTable::AlreadyDeclaredException &e)
    {
      error("Local model variable " + *name + " declared twice.");
    }

  model_tree->AddLocalVariable(*name, rhs);
  delete name;
}

void
ParsingDriver::change_type(SymbolType new_type, vector<string *> *var_list)
{
  for(vector<string *>::iterator it = var_list->begin();
      it != var_list->end(); it++)
    {
      int id;
      try
        {
          id = mod_file->symbol_table.getID(**it);
        }
      catch(SymbolTable::UnknownSymbolNameException &e)
        {
          error("Unknown variable " + **it);
        }

      // Check if symbol already used in a VariableNode
      if (mod_file->expressions_tree.isSymbolUsed(id)
          || mod_file->dynamic_model.isSymbolUsed(id))
        error("You cannot modify the type of symbol " + **it + " after having used it in an expression");

      mod_file->symbol_table.changeType(id, new_type);

      delete *it;
    }
  delete var_list;
}

NodeID
ParsingDriver::add_plus(NodeID arg1, NodeID arg2)
{
  return data_tree->AddPlus(arg1, arg2);
}

NodeID
ParsingDriver::add_minus(NodeID arg1, NodeID arg2)
{
  return data_tree->AddMinus(arg1, arg2);
}

NodeID
ParsingDriver::add_uminus(NodeID arg1)
{
  return data_tree->AddUMinus(arg1);
}

NodeID
ParsingDriver::add_times(NodeID arg1, NodeID arg2)
{
  return data_tree->AddTimes(arg1, arg2);
}

NodeID
ParsingDriver::add_divide(NodeID arg1, NodeID arg2)
{
  return data_tree->AddDivide(arg1, arg2);
}

NodeID
ParsingDriver::add_less(NodeID arg1, NodeID arg2)
{
  return data_tree->AddLess(arg1, arg2);
}

NodeID
ParsingDriver::add_greater(NodeID arg1, NodeID arg2)
{
  return data_tree->AddGreater(arg1, arg2);
}

NodeID
ParsingDriver::add_less_equal(NodeID arg1, NodeID arg2)
{
  return data_tree->AddLessEqual(arg1, arg2);
}

NodeID
ParsingDriver::add_greater_equal(NodeID arg1, NodeID arg2)
{
  return data_tree->AddGreaterEqual(arg1, arg2);
}

NodeID
ParsingDriver::add_equal_equal(NodeID arg1, NodeID arg2)
{
  return data_tree->AddEqualEqual(arg1, arg2);
}

NodeID
ParsingDriver::add_different(NodeID arg1, NodeID arg2)
{
  return data_tree->AddDifferent(arg1, arg2);
}

NodeID
ParsingDriver::add_power(NodeID arg1, NodeID arg2)
{
  return data_tree->AddPower(arg1, arg2);
}

NodeID
ParsingDriver::add_expectation(string *arg1, NodeID arg2)
{
  NodeID expectationNode = data_tree->AddExpectation(atoi(arg1->c_str()), arg2);
  delete arg1;
  return expectationNode;
}

NodeID
ParsingDriver::add_exp(NodeID arg1)
{
  return data_tree->AddExp(arg1);
}

NodeID
ParsingDriver::add_log(NodeID arg1)
{
  return data_tree->AddLog(arg1);
}

NodeID
ParsingDriver::add_log10(NodeID arg1)
{
  return data_tree->AddLog10(arg1);
}

NodeID
ParsingDriver::add_cos(NodeID arg1)
{
  return data_tree->AddCos(arg1);
}

NodeID
ParsingDriver::add_sin(NodeID arg1)
{
  return data_tree->AddSin(arg1);
}

NodeID
ParsingDriver::add_tan(NodeID arg1)
{
  return data_tree->AddTan(arg1);
}

NodeID
ParsingDriver::add_acos(NodeID arg1)
{
  return data_tree->AddAcos(arg1);
}

NodeID
ParsingDriver::add_asin(NodeID arg1)
{
  return data_tree->AddAsin(arg1);
}

NodeID
ParsingDriver::add_atan(NodeID arg1)
{
  return data_tree->AddAtan(arg1);
}

NodeID
ParsingDriver::add_cosh(NodeID arg1)
{
  return data_tree->AddCosh(arg1);
}

NodeID
ParsingDriver::add_sinh(NodeID arg1)
{
  return data_tree->AddSinh(arg1);
}

NodeID
ParsingDriver::add_tanh(NodeID arg1)
{
  return data_tree->AddTanh(arg1);
}

NodeID
ParsingDriver::add_acosh(NodeID arg1)
{
  return data_tree->AddAcosh(arg1);
}

NodeID
ParsingDriver::add_asinh(NodeID arg1)
{
  return data_tree->AddAsinh(arg1);
}

NodeID
ParsingDriver::add_atanh(NodeID arg1)
{
  return data_tree->AddAtanh(arg1);
}

NodeID
ParsingDriver::add_sqrt(NodeID arg1)
{
  return data_tree->AddSqrt(arg1);
}

NodeID
ParsingDriver::add_max(NodeID arg1, NodeID arg2)
{
  return data_tree->AddMax(arg1,arg2);
}

NodeID
ParsingDriver::add_min(NodeID arg1, NodeID arg2)
{
  return data_tree->AddMin(arg1,arg2);
}

NodeID
ParsingDriver::add_normcdf(NodeID arg1, NodeID arg2, NodeID arg3)
{
  return data_tree->AddNormcdf(arg1,arg2,arg3);
}

NodeID
ParsingDriver::add_normcdf(NodeID arg)
{
  return add_normcdf(arg, data_tree->Zero, data_tree->One);
}

NodeID
ParsingDriver::add_steady_state(NodeID arg1)
{
  return data_tree->AddSteadyState(arg1);
}

void
ParsingDriver::add_unknown_function_arg(NodeID arg)
{
  unknown_function_args.push_back(arg);
}

NodeID
ParsingDriver::add_unknown_function(string *function_name)
{
  if (mod_file->symbol_table.exists(*function_name))
    {
      if (mod_file->symbol_table.getType(*function_name) != eUnknownFunction)
        error("Symbol " + *function_name + " is not a function name.");
    }
  else
    mod_file->symbol_table.addSymbol(*function_name, eUnknownFunction);

  NodeID id = data_tree->AddUnknownFunction(*function_name, unknown_function_args);
  unknown_function_args.clear();
  return id;
}

void
ParsingDriver::add_native(const char *s)
{
  mod_file->addStatement(new NativeStatement(s));
}
