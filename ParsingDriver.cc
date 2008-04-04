/*
 * Copyright (C) 2003-2008 Dynare Team
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

#include "ParsingDriver.hh"
#include "Statement.hh"

ParsingDriver::ParsingDriver()
{
}

ParsingDriver::~ParsingDriver()
{
}

bool
ParsingDriver::symbol_exists_and_is_not_modfile_local_variable(const char *s)
{
  if (!mod_file->symbol_table.exists(s))
    return false;

  return(mod_file->symbol_table.getType(s) != eModFileLocalVariable);
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

  tmp_symbol_table = new TmpSymbolTable(mod_file->symbol_table);

  reset_data_tree();

  lexer = new DynareFlex(&in);
  lexer->set_debug(debug);

  Dynare::parser parser(*this);
  parser.set_debug_level(debug);
  parser.parse();

  delete lexer;
  delete tmp_symbol_table;

  return mod_file;
}

void
ParsingDriver::error(const Dynare::parser::location_type &l, const string &m)
{
  cerr << "ERROR: " << l << ": " << m << endl;
  exit(-1);
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
ParsingDriver::declare_symbol(string *name, Type type, string *tex_name)
{
  try
    {
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

NodeID
ParsingDriver::add_constant(string *constant)
{
  NodeID id = data_tree->AddNumConstant(*constant);
  delete constant;
  return id;
}

NodeID
ParsingDriver::add_model_variable(string *name)
{
  check_symbol_existence(*name);
  NodeID id = model_tree->AddVariable(*name);

  Type type = mod_file->symbol_table.getType(*name);

  if (type == eModFileLocalVariable)
    error("Variable " + *name + " not allowed inside model declaration. Its scope is only outside model.");

  if ((type == eEndogenous) && (model_tree->mode == eSparseDLLMode || model_tree->mode == eSparseMode))
    {
      int ID = mod_file->symbol_table.getID(*name);
      model_tree->block_triangular.fill_IM(model_tree->equation_number(), ID, 0);
    }

  delete name;
  return id;
}

NodeID
ParsingDriver::add_model_variable(string *name, string *olag)
{
  check_symbol_existence(*name);
  Type type = mod_file->symbol_table.getType(*name);
  int lag = atoi(olag->c_str());

  if (type == eModFileLocalVariable)
    error("Variable " + *name + " not allowed inside model declaration. Its scope is only outside model.");

  if ((type == eExogenous) && lag != 0)
    {
      ostringstream ost;
      ost << "Exogenous variable " << *name << " has lag " << lag;
      warning(ost.str());
    }

  NodeID id = model_tree->AddVariable(*name, lag);

  if ((type == eEndogenous) && (model_tree->mode == eSparseDLLMode || model_tree->mode == eSparseMode))
    model_tree->block_triangular.fill_IM(model_tree->equation_number(), mod_file->symbol_table.getID(*name), lag);

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

  NodeID id = data_tree->AddVariable(*name);

  delete name;
  return id;
}

void
ParsingDriver::periods(string *periods)
{
  int periods_val = atoi(periods->c_str());
  mod_file->addStatement(new PeriodsStatement(periods_val));
  delete periods;
}

void
ParsingDriver::cutoff(string *cutoff)
{
  double cutoff_val = atof(cutoff->c_str());
  mod_file->addStatement(new CutoffStatement(cutoff_val));
  delete cutoff;
}

void
ParsingDriver::markowitz(string *markowitz)
{
  double markowitz_val = atof(markowitz->c_str());
  mod_file->addStatement(new MarkowitzStatement(markowitz_val));
  delete markowitz;
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
  if (mod_file->symbol_table.getType(*name) != eParameter)
    error(*name + " is not a parameter");

  mod_file->addStatement(new InitParamStatement(*name, rhs, mod_file->symbol_table));

  // Update global eval context
  try
    {
      double val = rhs->eval(mod_file->global_eval_context);
      int symb_id = mod_file->symbol_table.getID(*name);
      mod_file->global_eval_context[make_pair(symb_id, eParameter)] = val;
     }
  catch(ExprNode::EvalException &e)
    {
    }

  delete name;
}

void
ParsingDriver::init_val(string *name, NodeID rhs)
{
  check_symbol_existence(*name);
  Type type = mod_file->symbol_table.getType(*name);

  if (type != eEndogenous
      && type != eExogenous
      && type != eExogenousDet)
    error("initval/endval: " + *name + " should be an endogenous or exogenous variable");

  init_values.push_back(make_pair(*name, rhs));

  // Update global evaluation context
  try
    {
      double val = rhs->eval(mod_file->global_eval_context);
      int symb_id = mod_file->symbol_table.getID(*name);
      mod_file->global_eval_context[make_pair(symb_id, type)] = val;
    }
  catch(ExprNode::EvalException &e)
    {
    }

  delete name;
}

void
ParsingDriver::init_val_filename(string *filename)
{
  options_list.num_options["INITVAL_FILE"] = 1;
  options_list.string_options["INITVAL_FILENAME"] = *filename;
  delete filename;
}

void
ParsingDriver::hist_val(string *name, string *lag, NodeID rhs)
{
  check_symbol_existence(*name);
  Type type = mod_file->symbol_table.getType(*name);

  if (type != eEndogenous
      && type != eExogenous
      && type != eExogenousDet)
    error("hist_val: " + *name + " should be an endogenous or exogenous variable");

  int ilag = atoi(lag->c_str());
  pair<string, int> key(*name, ilag);

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
  Type type = mod_file->symbol_table.getType(*name);

  if (type != eParameter
      && type != eExogenous
      && type != eExogenousDet)
    error("homotopy_val: " + *name + " should be a parameter or exogenous variable");

  homotopy_values.push_back(make_pair(*name, make_pair(val1, val2)));

  delete name;
}

void
ParsingDriver::use_dll()
{
  model_tree->mode = eDLLMode;
}

void
ParsingDriver::sparse_dll()
{
  model_tree->mode = eSparseDLLMode;
  model_tree->block_triangular.init_incidence_matrix(mod_file->symbol_table.endo_nbr);
}

void
ParsingDriver::sparse()
{
  model_tree->mode = eSparseMode;
  model_tree->block_triangular.init_incidence_matrix(mod_file->symbol_table.endo_nbr);
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
  set_current_data_tree(&mod_file->model_tree);
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
  mod_file->addStatement(new MShocksStatement(det_shocks, var_shocks, std_shocks,
                                              covar_shocks, corr_shocks, mod_file->symbol_table));
  det_shocks.clear();
  var_shocks.clear();
  std_shocks.clear();
  covar_shocks.clear();
  corr_shocks.clear();
}

void
ParsingDriver::add_det_shock(string *var)
{
  check_symbol_existence(*var);
  Type type = mod_file->symbol_table.getType(*var);
  if (type != eExogenous && type != eExogenousDet)
    error("shocks: shocks can only be applied to exogenous variables");

  if (det_shocks.find(*var) != det_shocks.end())
    error("shocks: variable " + *var + " declared twice");

  if (det_shocks_periods.size() != det_shocks_values.size())
    error("shocks: variable " + *var + ": number of periods is different from number of shock values");

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
ParsingDriver::do_sigma_e()
{
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
  if (mod_file->model_tree.mode == eSparseMode)
    mod_file->addStatement(new SteadySparseStatement(options_list));
  else
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
  if (options_list.num_options.find(name_option)
      != options_list.num_options.end())
    error("option " + name_option + " declared twice");

  if ((name_option == "periods") && (mod_file->model_tree.mode == eSparseDLLMode || mod_file->model_tree.mode == eSparseMode))
    mod_file->model_tree.block_triangular.periods = atoi(opt.c_str());
  else if (name_option == "cutoff")
    mod_file->model_tree.cutoff = atof(opt.c_str());

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
ParsingDriver::option_str_lst(const string &name_option)
{
  if (options_list.string_list_options.find(name_option)
      != options_list.string_list_options.end())
    error("option " + name_option + " declared twice");

  options_list.string_list_options[name_option] = new TmpSymbolTable::TmpSymbolTable(*tmp_symbol_table);
  tmp_symbol_table->clear();
}




void
ParsingDriver::linear()
{
  mod_file->linear = true;
}

void
ParsingDriver::add_tmp_var(string *tmp_var1, string *tmp_var2)
{
  check_symbol_existence(*tmp_var1);
  check_symbol_existence(*tmp_var2);
  tmp_symbol_table->AddTempSymbol(*tmp_var1, *tmp_var2);
  delete tmp_var1;
  delete tmp_var2;
}

void
ParsingDriver::add_tmp_var(string *tmp_var)
{
  check_symbol_existence(*tmp_var);
  tmp_symbol_table->AddTempSymbol(*tmp_var);
  delete tmp_var;
}

void ParsingDriver::rplot()
{
  mod_file->addStatement(new RplotStatement(*tmp_symbol_table, options_list));
  options_list.clear();
  tmp_symbol_table->clear();
}

void ParsingDriver::stoch_simul()
{
  mod_file->addStatement(new StochSimulStatement(*tmp_symbol_table, options_list));
  tmp_symbol_table->clear();
  options_list.clear();
}

void ParsingDriver::simulate()
{
  if (mod_file->model_tree.mode == eSparseDLLMode || mod_file->model_tree.mode == eSparseMode)
    simul_sparse();
  else
    simul();
}

void
ParsingDriver::simul_sparse()
{
  mod_file->addStatement(new SimulSparseStatement(options_list, mod_file->model_tree.compiler, mod_file->model_tree.mode));
  options_list.clear();
}

void
ParsingDriver::init_compiler(int compiler_type)
{
  mod_file->model_tree.compiler = compiler_type;
}

void
ParsingDriver::simul()
{
  mod_file->addStatement(new SimulStatement(options_list));
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
  estim_params.clear();
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
  mod_file->addStatement(new UnitRootVarsStatement(*tmp_symbol_table));
  tmp_symbol_table->clear();
}

void
ParsingDriver::run_estimation()
{
  mod_file->addStatement(new EstimationStatement(*tmp_symbol_table, options_list));
  tmp_symbol_table->clear();
  options_list.clear();
}

void
ParsingDriver::run_prior_analysis()
{
  mod_file->addStatement(new PriorAnalysisStatement(*tmp_symbol_table, options_list));
  tmp_symbol_table->clear();
  options_list.clear();
}

void
ParsingDriver::run_posterior_analysis()
{
  mod_file->addStatement(new PosteriorAnalysisStatement(*tmp_symbol_table, options_list));
  tmp_symbol_table->clear();
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
  mod_file->addStatement(new VarobsStatement(*tmp_symbol_table));
  tmp_symbol_table->clear();
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
  mod_file->addStatement(new OsrParamsStatement(*tmp_symbol_table));
  tmp_symbol_table->clear();
}

void
ParsingDriver::run_osr()
{
  mod_file->addStatement(new OsrStatement(*tmp_symbol_table, options_list));
  tmp_symbol_table->clear();
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
ParsingDriver::run_dynatype(string *filename, string *ext)
{
  mod_file->addStatement(new DynaTypeStatement(*tmp_symbol_table, *filename, *ext));
  tmp_symbol_table->clear();
  delete filename;
  delete ext;
}

void
ParsingDriver::run_dynasave(string *filename, string *ext)
{
  mod_file->addStatement(new DynaSaveStatement(*tmp_symbol_table, *filename, *ext));
  tmp_symbol_table->clear();
  delete filename;
  delete ext;
}

void
ParsingDriver::add_mc_filename(string *filename, string *prior)
{
  if (filename_list.find(*filename) != filename_list.end())
    error("model_comparison: filename " + *filename + " declared twice");
  filename_list[*filename] = *prior;
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
  set_current_data_tree(new ModelTree(mod_file->symbol_table, mod_file->num_constants));
}

void
ParsingDriver::end_planner_objective(NodeID expr)
{
  // Add equation corresponding to expression
  NodeID eq = model_tree->AddEqual(expr, model_tree->Zero);
  model_tree->addEquation(eq);

  mod_file->addStatement(new PlannerObjectiveStatement(model_tree));

  reset_data_tree();
}

void
ParsingDriver::ramsey_policy()
{
  mod_file->addStatement(new RamseyPolicyStatement(*tmp_symbol_table, options_list));
  tmp_symbol_table->clear();
  options_list.clear();
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

  model_tree->AddLocalParameter(*name, rhs);
  delete name;
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
  return data_tree->AddACos(arg1);
}

NodeID
ParsingDriver::add_asin(NodeID arg1)
{
  return data_tree->AddASin(arg1);
}

NodeID
ParsingDriver::add_atan(NodeID arg1)
{
  return data_tree->AddATan(arg1);
}

NodeID
ParsingDriver::add_cosh(NodeID arg1)
{
  return data_tree->AddCosH(arg1);
}

NodeID
ParsingDriver::add_sinh(NodeID arg1)
{
  return data_tree->AddSinH(arg1);
}

NodeID
ParsingDriver::add_tanh(NodeID arg1)
{
  return data_tree->AddTanH(arg1);
}

NodeID
ParsingDriver::add_acosh(NodeID arg1)
{
  return data_tree->AddACosH(arg1);
}

NodeID
ParsingDriver::add_asinh(NodeID arg1)
{
  return data_tree->AddASinH(arg1);
}

NodeID
ParsingDriver::add_atanh(NodeID arg1)
{
  return data_tree->AddATanH(arg1);
}

NodeID
ParsingDriver::add_sqrt(NodeID arg1)
{
  return data_tree->AddSqRt(arg1);
}

NodeID
ParsingDriver::add_max(NodeID arg1, NodeID arg2)
{
  return data_tree->AddMaX(arg1,arg2);
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
