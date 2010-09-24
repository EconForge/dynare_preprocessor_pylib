/*
 * Copyright (C) 2003-2010 Dynare Team
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
#include <cassert>
#include <iostream>
#include <sstream>

using namespace std;

#include "ComputingTasks.hh"
#include "Statement.hh"

SteadyStatement::SteadyStatement(const OptionsList &options_list_arg) :
  options_list(options_list_arg)
{
}

void
SteadyStatement::checkPass(ModFileStructure &mod_file_struct)
{
  mod_file_struct.steady_present = true;
}

void
SteadyStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output);
  output << "steady;\n";
}

CheckStatement::CheckStatement(const OptionsList &options_list_arg) :
  options_list(options_list_arg)
{
}

void
CheckStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output);
  output << "check;\n";
}

void
CheckStatement::checkPass(ModFileStructure &mod_file_struct)
{
  mod_file_struct.check_present = true;
}

ModelInfoStatement::ModelInfoStatement(const OptionsList &options_list_arg) :
  options_list(options_list_arg)
{
}

void
ModelInfoStatement::checkPass(ModFileStructure &mod_file_struct)
{
  //mod_file_struct.model_info_present = true;
}

void
ModelInfoStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output);
  output << "model_info();\n";
}

SimulStatement::SimulStatement(const OptionsList &options_list_arg) :
  options_list(options_list_arg)
{
}

void
SimulStatement::checkPass(ModFileStructure &mod_file_struct)
{
  mod_file_struct.simul_present = true;
}

void
SimulStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output);
  output << "simul();\n";
}

StochSimulStatement::StochSimulStatement(const SymbolList &symbol_list_arg,
                                         const OptionsList &options_list_arg) :
  symbol_list(symbol_list_arg),
  options_list(options_list_arg)
{
}

void
StochSimulStatement::checkPass(ModFileStructure &mod_file_struct)
{
  mod_file_struct.stoch_simul_present = true;

  // Fill in option_order of mod_file_struct
  OptionsList::num_options_t::const_iterator it = options_list.num_options.find("order");
  if (it != options_list.num_options.end())
    mod_file_struct.order_option = max(mod_file_struct.order_option, atoi(it->second.c_str()));

  // Fill in mod_file_struct.partial_information
  it = options_list.num_options.find("partial_information");
  if (it != options_list.num_options.end() && it->second == "1")
    mod_file_struct.partial_information = true;

  // Option k_order_solver (implicit when order >= 3)
  it = options_list.num_options.find("k_order_solver");
  if ((it != options_list.num_options.end() && it->second == "1")
      || mod_file_struct.order_option >= 3)
    mod_file_struct.k_order_solver = true;

  // Check that option "pruning" is not used with k-order
  it = options_list.num_options.find("pruning");
  if ((it != options_list.num_options.end() && it->second == "1")
      && mod_file_struct.k_order_solver)
    {
      cerr << "ERROR: in 'stoch_simul', you cannot use option 'pruning' with 'k_order_solver' option or with 3rd order approximation" << endl;
      exit(EXIT_FAILURE);
    }
}

void
StochSimulStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output);
  symbol_list.writeOutput("var_list_", output);
  output << "info = stoch_simul(var_list_);" << endl;
}

ForecastStatement::ForecastStatement(const SymbolList &symbol_list_arg,
                                     const OptionsList &options_list_arg) :
  symbol_list(symbol_list_arg),
  options_list(options_list_arg)
{
}

void
ForecastStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output);
  symbol_list.writeOutput("var_list_", output);
  output << "info = forecast(var_list_,'simul');\n";
}

RamseyPolicyStatement::RamseyPolicyStatement(const SymbolList &symbol_list_arg,
                                             const OptionsList &options_list_arg) :
  symbol_list(symbol_list_arg),
  options_list(options_list_arg)
{
}

void
RamseyPolicyStatement::checkPass(ModFileStructure &mod_file_struct)
{
  mod_file_struct.ramsey_policy_present = true;

  /* Fill in option_order of mod_file_struct
     Since ramsey policy needs one further order of derivation (for example, for 1st order
     approximation, it needs 2nd derivatives), we add 1 to the order declared by user */
  OptionsList::num_options_t::const_iterator it = options_list.num_options.find("order");
  if (it != options_list.num_options.end())
    {
      int order = atoi(it->second.c_str());
      if (order > 1)
        {
          cerr << "ERROR: ramsey_policy: order > 1 is not yet implemented" << endl;
          exit(EXIT_FAILURE);
        }
      mod_file_struct.order_option = max(mod_file_struct.order_option, order + 1);
    }

  // Fill in mod_file_struct.partial_information
  it = options_list.num_options.find("partial_information");
  if (it != options_list.num_options.end() && it->second == "1")
    mod_file_struct.partial_information = true;

  // Option k_order_solver (implicit when order >= 3)
  it = options_list.num_options.find("k_order_solver");
  if ((it != options_list.num_options.end() && it->second == "1")
      || mod_file_struct.order_option >= 3)
    mod_file_struct.k_order_solver = true;
}

void
RamseyPolicyStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output);
  symbol_list.writeOutput("var_list_", output);
  output << "ramsey_policy(var_list_);\n";
}

EstimationStatement::EstimationStatement(const SymbolList &symbol_list_arg,
                                         const OptionsList &options_list_arg,
                                         const SymbolTable &symbol_table_arg) :
  symbol_list(symbol_list_arg),
  options_list(options_list_arg),
  symbol_table(symbol_table_arg)
{
}

void
EstimationStatement::checkPass(ModFileStructure &mod_file_struct)
{
  mod_file_struct.estimation_present = true;

  // Fill in option_order of mod_file_struct
  OptionsList::num_options_t::const_iterator it = options_list.num_options.find("order");
  if (it != options_list.num_options.end())
    mod_file_struct.order_option = max(mod_file_struct.order_option, atoi(it->second.c_str()));

  // Fill in mod_file_struct.partial_information
  it = options_list.num_options.find("partial_information");
  if (it != options_list.num_options.end() && it->second == "1")
    mod_file_struct.partial_information = true;

  // Fill in mod_file_struct.dsge_var_calibrated
  it = options_list.num_options.find("dsge_var");
  if (it != options_list.num_options.end())
    mod_file_struct.dsge_var_calibrated = it->second;

  // Fill in mod_file_struct.dsge_var_estimated
  OptionsList::string_options_t::const_iterator it_str = options_list.string_options.find("dsge_var");
  if (it_str != options_list.string_options.end())
    mod_file_struct.dsge_var_estimated = true;

  // Fill in mod_file_struct.bayesian_irf_present
  it = options_list.num_options.find("bayesian_irf");
  if (it != options_list.num_options.end() && it->second == "1")
    mod_file_struct.bayesian_irf_present = true;

  it = options_list.num_options.find("dsge_varlag");
  if (it != options_list.num_options.end())
    if (mod_file_struct.dsge_var_calibrated.empty() &&
        !mod_file_struct.dsge_var_estimated)
      {
        cerr << "ERROR: The estimation statement requires a dsge_var option to be passed "
             << "if the dsge_varlag option is passed." << endl;
        exit(EXIT_FAILURE);
      }

  if (!mod_file_struct.dsge_var_calibrated.empty() &&
      mod_file_struct.dsge_var_estimated)
    {
      cerr << "ERROR: An estimation statement cannot take more than one dsge_var option." << endl;
      exit(EXIT_FAILURE);
    }
}

void
EstimationStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output);
  symbol_list.writeOutput("var_list_", output);
  output << "dynare_estimation(var_list_);\n";
}

DynareSensitivityStatement::DynareSensitivityStatement(const OptionsList &options_list_arg) :
  options_list(options_list_arg)
{
}

void
DynareSensitivityStatement::checkPass(ModFileStructure &mod_file_struct)
{
  OptionsList::num_options_t::const_iterator it = options_list.num_options.find("identification");
  if (it != options_list.num_options.end()
      && it->second == "1")
    mod_file_struct.identification_present = true;
}

void
DynareSensitivityStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output, "options_gsa");
  output << "dynare_sensitivity(options_gsa);" << endl;
}

RplotStatement::RplotStatement(const SymbolList &symbol_list_arg,
                               const OptionsList &options_list_arg) :
  symbol_list(symbol_list_arg),
  options_list(options_list_arg)
{
}

void
RplotStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output);
  symbol_list.writeOutput("var_list_", output);
  output << "rplot(var_list_);\n";
}

UnitRootVarsStatement::UnitRootVarsStatement(const SymbolList &symbol_list_arg) :
  symbol_list(symbol_list_arg)
{
}

void
UnitRootVarsStatement::writeOutput(ostream &output, const string &basename) const
{
  symbol_list.writeOutput("options_.unit_root_vars", output);
}

PeriodsStatement::PeriodsStatement(int periods_arg) : periods(periods_arg)
{
}

void
PeriodsStatement::writeOutput(ostream &output, const string &basename) const
{
  output << "options_.periods = " << periods << ";" << endl;
}

DsampleStatement::DsampleStatement(int val1_arg) : val1(val1_arg), val2(-1)
{
}

DsampleStatement::DsampleStatement(int val1_arg, int val2_arg) : val1(val1_arg), val2(val2_arg)
{
}

void
DsampleStatement::writeOutput(ostream &output, const string &basename) const
{
  if (val2 < 0)
    output << "dsample(" << val1 << ");" << endl;
  else
    output << "dsample(" << val1 << ", " << val2 << ");" << endl;
}

EstimatedParamsStatement::EstimatedParamsStatement(const vector<EstimationParams> &estim_params_list_arg,
                                                   const SymbolTable &symbol_table_arg) :
  estim_params_list(estim_params_list_arg),
  symbol_table(symbol_table_arg)
{
}

void
EstimatedParamsStatement::checkPass(ModFileStructure &mod_file_struct)
{
  for (vector<EstimationParams>::const_iterator it = estim_params_list.begin();
       it != estim_params_list.end(); it++)
    {
      if (it->name == "dsge_prior_weight")
        mod_file_struct.dsge_prior_weight_in_estimated_params = true;

      // Handle case of degenerate beta prior
      if (it->prior == "1") //BETA_PDF is associated with "1" in DynareBison.yy
        try
          {
            if (it->mean->eval(eval_context_t()) == 0.5 &&
                it->std->eval(eval_context_t()) == 0.5)
              {
                cerr << "ERROR: The prior density is not defined for the beta distribution when the mean = standard deviation = 0.5." << endl;
                exit(EXIT_FAILURE);
              }
          }
        catch (ExprNode::EvalException &e)
          {
            // We don't have enough information to compute the numerical value, skip the test
          }
    }
}

void
EstimatedParamsStatement::writeOutput(ostream &output, const string &basename) const
{
  output << "global estim_params_" << endl
         << "estim_params_.var_exo = [];" << endl
         << "estim_params_.var_endo = [];" << endl
         << "estim_params_.corrx = [];" << endl
         << "estim_params_.corrn = [];" << endl
         << "estim_params_.param_vals = [];" << endl;

  vector<EstimationParams>::const_iterator it;

  for (it = estim_params_list.begin(); it != estim_params_list.end(); it++)
    {
      int symb_id = symbol_table.getTypeSpecificID(it->name) + 1;
      SymbolType symb_type = symbol_table.getType(it->name);

      switch (it->type)
        {
        case 1:
          if (symb_type == eExogenous)
            output << "estim_params_.var_exo = [estim_params_.var_exo; ";
          else if (symb_type == eEndogenous)
            output << "estim_params_.var_endo = [estim_params_.var_endo; ";
          output << symb_id;
          break;
        case 2:
          output << "estim_params_.param_vals = [estim_params_.param_vals; "
                 << symb_id;
          break;
        case 3:
          if (symb_type == eExogenous)
            output << "estim_params_.corrx = [estim_params_.corrx; ";
          else if (symb_type == eEndogenous)
            output << "estim_params_.corrn = [estim_params_.corrn; ";
          output << symb_id << " " << symbol_table.getTypeSpecificID(it->name2)+1;
          break;
        }
      output << ", ";
      it->init_val->writeOutput(output);
      output << ", ";
      it->low_bound->writeOutput(output);
      output << ", ";
      it->up_bound->writeOutput(output);
      output << ", "
             << it->prior << ", ";
      it->mean->writeOutput(output);
      output << ", ";
      it->std->writeOutput(output);
      output << ", ";
      it->p3->writeOutput(output);
      output << ", ";
      it->p4->writeOutput(output);
      output << ", ";
      it->jscale->writeOutput(output);
      output << " ];" << endl;
    }
}

EstimatedParamsInitStatement::EstimatedParamsInitStatement(const vector<EstimationParams> &estim_params_list_arg,
                                                           const SymbolTable &symbol_table_arg) :
  estim_params_list(estim_params_list_arg),
  symbol_table(symbol_table_arg)
{
}

void
EstimatedParamsInitStatement::writeOutput(ostream &output, const string &basename) const
{
  vector<EstimationParams>::const_iterator it;

  for (it = estim_params_list.begin(); it != estim_params_list.end(); it++)
    {
      int symb_id = symbol_table.getTypeSpecificID(it->name) + 1;
      SymbolType symb_type = symbol_table.getType(it->name);

      if (it->type < 3)
        {
          if (symb_type == eExogenous)
            {
              output << "tmp1 = find(estim_params_.var_exo(:,1)==" << symb_id << ");" << endl;
              output << "estim_params_.var_exo(tmp1,2) = ";
              it->init_val->writeOutput(output);
              output << ";" << endl;
            }
          else if (symb_type == eEndogenous)
            {
              output << "tmp1 = find(estim_params_.var_endo(:,1)==" << symb_id << ");" << endl;
              output << "estim_params_.var_endo(tmp1,2) = ";
              it->init_val->writeOutput(output);
              output << ";" << endl;
            }
          else if (symb_type == eParameter)
            {
              output << "tmp1 = find(estim_params_.param_vals(:,1)==" << symb_id << ");" << endl;
              output << "estim_params_.param_vals(tmp1,2) = ";
              it->init_val->writeOutput(output);
              output << ";" << endl;
            }
        }
      else
        {
          if (symb_type == eExogenous)
            {
              output << "tmp1 = find((estim_params_.corrx(:,1)==" << symb_id << ")) & (estim_params_.corrx(:,2)==" << symbol_table.getTypeSpecificID(it->name2)+1 << ");" << endl;
              output << "estim_params_.corrx(tmp1,3) = ";
              it->init_val->writeOutput(output);
              output << ";" << endl;
            }
          else if (symb_type == eEndogenous)
            {
              output << "tmp1 = find((estim_params_.corrn(:,1)==" << symb_id << ")) & (estim_params_.corrn(:,2)==" << symbol_table.getTypeSpecificID(it->name2)+1 << ";" << endl;
              output << "estim_params_.corrn(tmp1,3) = ";
              it->init_val->writeOutput(output);
              output << ";" << endl;
            }
        }
    }
}

EstimatedParamsBoundsStatement::EstimatedParamsBoundsStatement(const vector<EstimationParams> &estim_params_list_arg,
                                                               const SymbolTable &symbol_table_arg) :
  estim_params_list(estim_params_list_arg),
  symbol_table(symbol_table_arg)
{
}

void
EstimatedParamsBoundsStatement::writeOutput(ostream &output, const string &basename) const
{
  vector<EstimationParams>::const_iterator it;

  for (it = estim_params_list.begin(); it != estim_params_list.end(); it++)
    {
      int symb_id = symbol_table.getTypeSpecificID(it->name) + 1;
      SymbolType symb_type = symbol_table.getType(it->name);

      if (it->type < 3)
        {
          if (symb_type == eExogenous)
            {
              output << "tmp1 = find(estim_params_.var_exo(:,1)==" << symb_id << ");" << endl;

              output << "estim_params_.var_exo(tmp1,3) = ";
              it->low_bound->writeOutput(output);
              output << ";" << endl;

              output << "estim_params_.var_exo(tmp1,4) = ";
              it->up_bound->writeOutput(output);
              output << ";" << endl;
            }
          else if (symb_type == eEndogenous)
            {
              output << "tmp1 = find(estim_params_.var_endo(:,1)==" << symb_id << ");" << endl;

              output << "estim_params_.var_endo(tmp1,3) = ";
              it->low_bound->writeOutput(output);
              output << ";" << endl;

              output << "estim_params_.var_endo(tmp1,4) = ";
              it->up_bound->writeOutput(output);
              output << ";" << endl;
            }
          else if (symb_type == eParameter)
            {
              output << "tmp1 = find(estim_params_.param_vals(:,1)==" << symb_id << ");" << endl;

              output << "estim_params_.param_vals(tmp1,3) = ";
              it->low_bound->writeOutput(output);
              output << ";" << endl;

              output << "estim_params_.param_vals(tmp1,4) = ";
              it->up_bound->writeOutput(output);
              output << ";" << endl;
            }
        }
      else
        {
          if (symb_type == eExogenous)
            {
              output << "tmp1 = find((estim_params_.corrx(:,1)==" << symb_id << ")) & (estim_params_.corrx(:,2)==" << symbol_table.getTypeSpecificID(it->name2)+1 << ");" << endl;

              output << "estim_params_.corrx(tmp1,4) = ";
              it->low_bound->writeOutput(output);
              output << ";" << endl;

              output << "estim_params_.corrx(tmp1,5) = ";
              it->up_bound->writeOutput(output);
              output << ";" << endl;
            }
          else if (symb_type == eEndogenous)
            {
              output << "tmp1 = find((estim_params_.corrn(:,1)==" << symb_id << ")) & (estim_params_.corrn(:,2)==" << symbol_table.getTypeSpecificID(it->name2)+1 << ";" << endl;

              output << "estim_params_.corrn(tmp1,4) = ";
              it->low_bound->writeOutput(output);
              output << ";" << endl;

              output << "estim_params_.corrn(tmp1,5) = ";
              it->up_bound->writeOutput(output);
              output << ";" << endl;
            }
        }
    }
}

ObservationTrendsStatement::ObservationTrendsStatement(const trend_elements_t &trend_elements_arg,
                                                       const SymbolTable &symbol_table_arg) :
  trend_elements(trend_elements_arg),
  symbol_table(symbol_table_arg)
{
}

void
ObservationTrendsStatement::writeOutput(ostream &output, const string &basename) const
{
  output << "options_.trend_coeff_ = {};" << endl;

  trend_elements_t::const_iterator it;

  for (it = trend_elements.begin(); it != trend_elements.end(); it++)
    {
      SymbolType type = symbol_table.getType(it->first);
      if (type == eEndogenous)
        {
          output << "tmp1 = strmatch('" << it->first << "',options_.varobs,'exact');\n";
          output << "options_.trend_coeffs{tmp1} = '";
          it->second->writeOutput(output);
          output << "';" << endl;
        }
      else
        cout << "Error : Non-variable symbol used in TREND_COEFF: " << it->first << endl;
    }
}

CalibVarStatement::CalibVarStatement(const calib_var_t &calib_var_arg,
                                     const calib_covar_t &calib_covar_arg,
                                     const calib_ac_t &calib_ac_arg,
                                     const SymbolTable &symbol_table_arg) :
  calib_var(calib_var_arg),
  calib_covar(calib_covar_arg),
  calib_ac(calib_ac_arg),
  symbol_table(symbol_table_arg)
{
}

void
CalibVarStatement::writeOutput(ostream &output, const string &basename) const
{

  output << "%" << endl
         << "% CALIB_VAR" << endl
         << "%" << endl;

  for (int i = 1; i < 4; i++)
    {
      output << "calib_var_index{" << i << "} = [];\n";
      output << "calib_targets{" << i << "} = [];\n";
      output << "calib_weights{" << i << "}=[];\n";
    }

  // Print calibration variances
  for (calib_var_t::const_iterator it = calib_var.begin();
       it != calib_var.end(); it++)
    {
      const string &name = it->first;
      const string &weight = it->second.first;
      const expr_t expression = it->second.second;

      int id = symbol_table.getTypeSpecificID(name) + 1;
      if (symbol_table.getType(name) == eEndogenous)
        {
          output << "calib_var_index{1} = [calib_var_index{1};" <<  id << "," << id << "];\n";
          output << "calib_weights{1} = [calib_weights{1}; " << weight << "];\n";
          output << "calib_targets{1} =[calib_targets{1}; ";
          expression->writeOutput(output);
          output << "];\n";
        }
      else if (symbol_table.getType(name) == eExogenous)
        {
          output << "calib_var_index{3} = [calib_var_index{3};" <<  id << "," << id << "];\n";
          output << "calib_weights{3} = [calib_weights{3}; " << weight << "];\n";
          output << "calib_targets{3} =[calib_targets{3}; ";
          expression->writeOutput(output);
          output << "];\n";
        }
    }

  // Print calibration covariances
  for (calib_covar_t::const_iterator it = calib_covar.begin();
       it != calib_covar.end(); it++)
    {
      const string &name1 = it->first.first;
      const string &name2 = it->first.second;
      const string &weight = it->second.first;
      const expr_t expression = it->second.second;

      int id1 = symbol_table.getTypeSpecificID(name1) + 1;
      int id2 = symbol_table.getTypeSpecificID(name2) + 1;
      if (symbol_table.getType(name1) == eEndogenous)
        {
          output << "calib_var_index{1} = [calib_var_index{1};" <<  id1 << "," << id2 << "];\n";
          output << "calib_weights{1} = [calib_weights{1}; " << weight << "];\n";
          output << "calib_targets{1} =[calib_targets{1}; ";
          expression->writeOutput(output);
          output << "];\n";
        }
      else if (symbol_table.getType(name1) == eExogenous)
        {
          output << "calib_var_index{3} = [calib_var_index{3};" <<  id1 << "," << id2 << "];\n";
          output << "calib_weights{3} = [calib_weights{3}; " << weight << "];\n";
          output << "calib_targets{3} =[calib_targets{3}; ";
          expression->writeOutput(output);
          output << "];\n";
        }
    }

  // Print calibration autocorrelations
  int max_iar = 3;

  for (calib_ac_t::const_iterator it = calib_ac.begin();
       it != calib_ac.end(); it++)
    {
      const string &name = it->first.first;
      int iar = it->first.second + 3;
      const string &weight = it->second.first;
      const expr_t expression = it->second.second;

      int id = symbol_table.getTypeSpecificID(name) + 1;

      if (iar > max_iar)
        {
          // Create new variables
          for (int i = max_iar + 1; i <= iar; i++)
            {
              output << "calib_var_index{" << i << "} = [];\n";
              output << "calib_targets{" << i << "} = [];\n";
              output << "calib_weights{" << i << "}=[];\n";
            }
          max_iar = iar;
        }

      output << "calib_var_index{" << iar << "} = [calib_var_index{" << iar << "};" <<  id << "];\n";
      output << "calib_weights{" << iar << "} = [calib_weights{" << iar << "}; " << weight << "];\n";
      output << "calib_targets{" << iar << "} =[calib_targets{" << iar << "}; ";
      expression->writeOutput(output);
      output << "];\n";
    }
}

CalibStatement::CalibStatement(int covar_arg) : covar(covar_arg)
{
}

void
CalibStatement::writeOutput(ostream &output, const string &basename) const
{
  output << "M_.Sigma_e=calib(calib_var_index,calib_targets,calib_weights," << covar << ",Sigma_e_);\n";
}

OsrParamsStatement::OsrParamsStatement(const SymbolList &symbol_list_arg) :
  symbol_list(symbol_list_arg)
{
}

void
OsrParamsStatement::writeOutput(ostream &output, const string &basename) const
{
  symbol_list.writeOutput("osr_params_", output);
}

OsrStatement::OsrStatement(const SymbolList &symbol_list_arg,
                           const OptionsList &options_list_arg) :
  symbol_list(symbol_list_arg),
  options_list(options_list_arg)
{
}

void
OsrStatement::checkPass(ModFileStructure &mod_file_struct)
{
  mod_file_struct.osr_present = true;

  // Fill in option_order of mod_file_struct
  OptionsList::num_options_t::const_iterator it = options_list.num_options.find("order");
  if (it != options_list.num_options.end())
    mod_file_struct.order_option = max(mod_file_struct.order_option, atoi(it->second.c_str()));

  // Fill in mod_file_struct.partial_information
  it = options_list.num_options.find("partial_information");
  if (it != options_list.num_options.end() && it->second == "1")
    mod_file_struct.partial_information = true;

  // Option k_order_solver (implicit when order >= 3)
  it = options_list.num_options.find("k_order_solver");
  if ((it != options_list.num_options.end() && it->second == "1")
      || mod_file_struct.order_option >= 3)
    mod_file_struct.k_order_solver = true;
}

void
OsrStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output);
  symbol_list.writeOutput("var_list_", output);
  output << "osr(var_list_,osr_params_,obj_var_,optim_weights_);\n";
}

OptimWeightsStatement::OptimWeightsStatement(const var_weights_t &var_weights_arg,
                                             const covar_weights_t &covar_weights_arg,
                                             const SymbolTable &symbol_table_arg) :
  var_weights(var_weights_arg),
  covar_weights(covar_weights_arg),
  symbol_table(symbol_table_arg)
{
}

void
OptimWeightsStatement::writeOutput(ostream &output, const string &basename) const
{
  output << "%" << endl
         << "% OPTIM_WEIGHTS" << endl
         << "%" << endl
         << "optim_weights_ = sparse(M_.endo_nbr,M_.endo_nbr);" << endl
         << "obj_var_ = [];" << endl << endl;

  for (var_weights_t::const_iterator it = var_weights.begin();
       it != var_weights.end(); it++)
    {
      const string &name = it->first;
      const expr_t value = it->second;
      int id = symbol_table.getTypeSpecificID(name) + 1;
      output <<  "optim_weights_(" << id << "," << id << ") = ";
      value->writeOutput(output);
      output << ";" << endl;
      output << "obj_var_ = [obj_var_; " << id << "];\n";
    }

  for (covar_weights_t::const_iterator it = covar_weights.begin();
       it != covar_weights.end(); it++)
    {
      const string &name1 = it->first.first;
      const string &name2 = it->first.second;
      const expr_t value = it->second;
      int id1 = symbol_table.getTypeSpecificID(name1) + 1;
      int id2 = symbol_table.getTypeSpecificID(name2) + 1;
      output <<  "optim_weights_(" << id1 << "," << id2 << ") = ";
      value->writeOutput(output);
      output << ";" << endl;
      output << "obj_var_ = [obj_var_; " << id1 << "; " << id2 << "];\n";
    }
}

DynaSaveStatement::DynaSaveStatement(const SymbolList &symbol_list_arg,
                                     const string &filename_arg) :
  symbol_list(symbol_list_arg),
  filename(filename_arg)
{
}

void
DynaSaveStatement::writeOutput(ostream &output, const string &basename) const
{
  symbol_list.writeOutput("var_list_", output);
  output << "dynasave('" << filename
         << "',var_list_);" << endl;
}

DynaTypeStatement::DynaTypeStatement(const SymbolList &symbol_list_arg,
                                     const string &filename_arg) :
  symbol_list(symbol_list_arg),
  filename(filename_arg)
{
}

void
DynaTypeStatement::writeOutput(ostream &output, const string &basename) const
{
  symbol_list.writeOutput("var_list_", output);
  output << "dynatype('" << filename
         << "',var_list_);" << endl;
}

ModelComparisonStatement::ModelComparisonStatement(const filename_list_t &filename_list_arg,
                                                   const OptionsList &options_list_arg) :
  filename_list(filename_list_arg),
  options_list(options_list_arg)
{
}

void
ModelComparisonStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output);

  output << "ModelNames_ = {};" << endl;
  output << "ModelPriors_ = [];" << endl;

  for (filename_list_t::const_iterator it = filename_list.begin();
       it != filename_list.end(); it++)
    {
      output << "ModelNames_ = { ModelNames_{:} '" << (*it).first << "'};" << endl;
      output << "ModelPriors_ = [ ModelPriors_ ; " << (*it).second << "];" << endl;
    }
  output << "model_comparison(ModelNames_,ModelPriors_,oo_,options_,M_.fname);" << endl;
}

PlannerObjectiveStatement::PlannerObjectiveStatement(StaticModel *model_tree_arg) :
  model_tree(model_tree_arg)
{
}

PlannerObjectiveStatement::~PlannerObjectiveStatement()
{
  delete model_tree;
}

void
PlannerObjectiveStatement::checkPass(ModFileStructure &mod_file_struct)
{
  assert(model_tree->equation_number() == 1);
}

void
PlannerObjectiveStatement::computingPass()
{
  model_tree->computingPass(eval_context_t(), false, true, false, false);
}

void
PlannerObjectiveStatement::writeOutput(ostream &output, const string &basename) const
{
  model_tree->writeStaticFile(basename + "_objective", false, false);
}

BVARDensityStatement::BVARDensityStatement(int maxnlags_arg, const OptionsList &options_list_arg) :
  maxnlags(maxnlags_arg),
  options_list(options_list_arg)
{
}

void
BVARDensityStatement::checkPass(ModFileStructure &mod_file_struct)
{
  mod_file_struct.bvar_present = true;
}

void
BVARDensityStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output);
  output << "bvar_density(" << maxnlags << ");" << endl;
}

BVARForecastStatement::BVARForecastStatement(int nlags_arg, const OptionsList &options_list_arg) :
  nlags(nlags_arg),
  options_list(options_list_arg)
{
}

void
BVARForecastStatement::checkPass(ModFileStructure &mod_file_struct)
{
  mod_file_struct.bvar_present = true;
}

void
BVARForecastStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output);
  output << "bvar_forecast(" << nlags << ");" << endl;
}

SBVARStatement::SBVARStatement(const OptionsList &options_list_arg) :
  options_list(options_list_arg)
{
}

void
SBVARStatement::checkPass(ModFileStructure &mod_file_struct)
{
  mod_file_struct.bvar_present = true;
}

void
SBVARStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output);
  output << "swz_sbvar(0,M_,options_);" << endl;
}

MS_SBVARStatement::MS_SBVARStatement(const OptionsList &options_list_arg) :
  options_list(options_list_arg)
{
}

void
MS_SBVARStatement::checkPass(ModFileStructure &mod_file_struct)
{
  mod_file_struct.bvar_present = true;
}

void
MS_SBVARStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output);
  output << "swz_sbvar(1,M_,options_);" << endl;
}

IdentificationStatement::IdentificationStatement(const OptionsList &options_list_arg) :
  options_list(options_list_arg)
{
}

void
IdentificationStatement::checkPass(ModFileStructure &mod_file_struct)
{
  mod_file_struct.identification_present = true;
}

void
IdentificationStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output, "options_ident");
  output << "dynare_identification(options_ident);" << endl;
}

WriteLatexDynamicModelStatement::WriteLatexDynamicModelStatement(const DynamicModel &dynamic_model_arg) :
  dynamic_model(dynamic_model_arg)
{
}

void
WriteLatexDynamicModelStatement::writeOutput(ostream &output, const string &basename) const
{
  dynamic_model.writeLatexFile(basename);
}

WriteLatexStaticModelStatement::WriteLatexStaticModelStatement(const StaticModel &static_model_arg) :
  static_model(static_model_arg)
{
}

void
WriteLatexStaticModelStatement::writeOutput(ostream &output, const string &basename) const
{
  static_model.writeLatexFile(basename);
}

ShockDecompositionStatement::ShockDecompositionStatement(const SymbolList &symbol_list_arg,
                                                         const OptionsList &options_list_arg) :
  symbol_list(symbol_list_arg),
  options_list(options_list_arg)
{
}

void
ShockDecompositionStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output);
  symbol_list.writeOutput("var_list_", output);
  output << "oo_ = shock_decomposition(M_,oo_,options_,var_list_);\n";
}

ConditionalForecastStatement::ConditionalForecastStatement(const OptionsList &options_list_arg) :
  options_list(options_list_arg)
{
}

void
ConditionalForecastStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output, "options_cond_fcst_");
  output << "imcforecast(constrained_paths_, constrained_vars_, options_cond_fcst_);" << endl;
}

PlotConditionalForecastStatement::PlotConditionalForecastStatement(int periods_arg, const SymbolList &symbol_list_arg) :
  periods(periods_arg),
  symbol_list(symbol_list_arg)
{
}

void
PlotConditionalForecastStatement::writeOutput(ostream &output, const string &basename) const
{
  symbol_list.writeOutput("var_list_", output);
  if (periods == -1)
    output << "plot_icforecast(var_list_);" << endl;
  else
    output << "plot_icforecast(var_list_, " << periods << ");" << endl;
}

SvarIdentificationStatement::SvarIdentificationStatement(const svar_identification_exclusion_t &exclusion_arg,
                                                         const bool &upper_cholesky_present_arg,
                                                         const bool &lower_cholesky_present_arg,
                                                         const SymbolTable &symbol_table_arg) :
  exclusion(exclusion_arg),
  upper_cholesky_present(upper_cholesky_present_arg),
  lower_cholesky_present(lower_cholesky_present_arg),
  symbol_table(symbol_table_arg)
{
}

int
SvarIdentificationStatement::getMaxLag() const
{
  int max_lag = 0;
  for (svar_identification_exclusion_t::const_iterator it = exclusion.begin(); it != exclusion.end(); it++)
    if (it->first.first > max_lag)
      max_lag = it->first.first;

  return max_lag;
}

void
SvarIdentificationStatement::checkPass(ModFileStructure &mod_file_struct)
{
  if (!mod_file_struct.svar_identification_present)
    mod_file_struct.svar_identification_present = true;
  else
    {
      cerr << "ERROR: You may only have one svar_identification block in your .mod file." << endl;
      exit(EXIT_FAILURE);
    }
}

void
SvarIdentificationStatement::writeOutput(ostream &output, const string &basename) const
{
  if (upper_cholesky_present && lower_cholesky_present)
    {
      cerr << "SvarIdentificationStatement::writeOutput() Should not arrive here (1). Please report this to the Dynare Team." << endl;
      exit(EXIT_FAILURE);
    }

  output << "%" << endl
         << "% SVAR IDENTIFICATION" << endl
         << "%" << endl;

  if (upper_cholesky_present)
    output << "options_.ms.upper_cholesky=1;" << endl;

  if (lower_cholesky_present)
    output << "options_.ms.lower_cholesky=1;" << endl;

  if (!upper_cholesky_present && !lower_cholesky_present)
    {
      int n = symbol_table.endo_nbr();
//       int m = symbol_table.exo_nbr();
      int m = 1; // this is the constant, not the shocks
      int r = getMaxLag();
      int k = r*n+m;

      if (k < 1)
	{
	  cerr << "ERROR: lag = " << r
	       << ", number of endogenous variables = " << n
	       << ", number of exogenous variables = " << m
	       << ". If this is not a logical error in the specification"
	       << " of the .mod file, please report it to the Dynare Team." << endl;
	  exit(EXIT_FAILURE);
	}
      if (n < 1)
        {
          cerr << "ERROR: Number of endogenous variables = " << n << "< 1. If this is not a logical "
               << "error in the specification of the .mod file, please report it to the Dynare Team." << endl;
          exit(EXIT_FAILURE);
        }
      output << "options_.ms.Qi = zeros(" << n << ", " << n << ", " << n << ");" << endl;
      output << "options_.ms.Ri = zeros(" << k << ", " << k << ", " << n << ");" << endl;

      for (svar_identification_exclusion_t::const_iterator it = exclusion.begin(); it != exclusion.end(); it++)
        {
          for (unsigned int h = 0; h < it->second.size(); h++)
            {
              int j = it->second.at(h) + 1;
              int i = it->first.second;
              if (j < 1 || j > n || (int) h+1 > n || i < 1)
                {
                  cerr << "SvarIdentificationStatement::writeOutput() Should not arrive here (2). Please report this to the Dynare Team." << endl;
                  exit(EXIT_FAILURE);
                }
              if (i > n)
                {
                  cerr << "ERROR: equation number " << i << " is greater than the number of endogenous variables, " << n << "." << endl;
                  exit(EXIT_FAILURE);
                }

              if (it->first.first == 0)
                output << "options_.ms.Qi(" << h+1 << ", " << j << ", "<< i << ") = 1;" << endl;
              else if (it->first.first > 0)
                {
                  if ((it->first.first-1)*n+j > k)
                    {
                      cerr << "ERROR: lag =" << it->first.first << ", num endog vars = " << n << "current endog var index = " << j << ". Index "
                           << "out of bounds. If the above does not represent a logical error, please report this to the Dyanre Team." << endl;
                    }
                  output << "options_.ms.Ri(" << h+1 << ", " << (it->first.first-1)*n+j << ", "<< i << ") = 1;" << endl;
                }
              else
                {
                  cerr << "SvarIdentificationStatement::writeOutput() Should not arrive here (3). Please report this to the Dynare Team." << endl;
                  exit(EXIT_FAILURE);
                }
            }
        }
    }
}

MarkovSwitchingStatement::MarkovSwitchingStatement(const OptionsList &options_list_arg) :
  options_list(options_list_arg)
{
}

void
MarkovSwitchingStatement::writeOutput(ostream &output, const string &basename) const
{
  OptionsList::num_options_t::const_iterator itChain, itState, itNOS, itDuration;

  itChain = options_list.num_options.find("ms.chain");
  if (itChain == options_list.num_options.end())
    {
      cerr << "MarkovSwitchingStatement::writeOutput() Should not arrive here (1). Please report this to the Dynare Team." << endl;
      exit(EXIT_FAILURE);
    }

  itDuration = options_list.num_options.find("ms.duration");
  if (itDuration == options_list.num_options.end())
    {
      cerr << "MarkovSwitchingStatement::writeOutput() Should not arrive here (2). Please report this to the Dynare Team." << endl;
      exit(EXIT_FAILURE);
    }

  itState = options_list.num_options.find("ms.state");
  itNOS = options_list.num_options.find("ms.number_of_states");
  if (itState != options_list.num_options.end()
      && itNOS == options_list.num_options.end())
    output << "options_.ms.ms_chain(" << itChain->second << ").state(" << itState->second << ").duration = " << itDuration->second << ";" << endl;
  else if (itState == options_list.num_options.end()
           && itNOS != options_list.num_options.end())
    for (int i = 0; i < atoi(itNOS->second.c_str()); i++)
      output << "options_.ms.ms_chain(" << itChain->second << ").state(" << i+1 << ").duration = " << itDuration->second << ";" << endl;
  else
    {
      cerr << "MarkovSwitchingStatement::writeOutput() Should not arrive here (3). Please report this to the Dynare Team." << endl;
      exit(EXIT_FAILURE);
    }
}

SvarStatement::SvarStatement(const OptionsList &options_list_arg) :
  options_list(options_list_arg)
{
}

void
SvarStatement::writeOutput(ostream &output, const string &basename) const
{
  OptionsList::num_options_t::const_iterator it0, it1, it2;
  OptionsList::vec_int_options_t::const_iterator itv;

  it0 = options_list.num_options.find("ms.chain");
  if (it0 != options_list.num_options.end())
    output << "options_.ms.ms_chain(" << it0->second << ")";
  else
    {
      cerr << "SvarStatement::writeOutput() Should not arrive here (1). Please report this to the Dynare Team." << endl;
      exit(EXIT_FAILURE);
    }

  it0 = options_list.string_options.find("ms.coefficients");
  it1 = options_list.string_options.find("ms.variances");
  it2 = options_list.string_options.find("ms.constants");
  if (it0 != options_list.string_options.end()
      && it1 == options_list.string_options.end()
      && it2 == options_list.string_options.end())
    output << "." << it0->second;
  else if (it0 == options_list.string_options.end()
           && it1 != options_list.string_options.end()
           && it2 == options_list.string_options.end())
    output << "." << it1->second;
  else if (it0 == options_list.string_options.end()
           && it1 == options_list.string_options.end()
           && it2 != options_list.string_options.end())
    output << "." << it2->second;
  else
    {
      cerr << "SvarStatement::writeOutput() Should not arrive here (2). Please report this to the Dynare Team." << endl;
      exit(EXIT_FAILURE);
    }

  itv = options_list.vector_int_options.find("ms.equations");
  output << ".equations = ";
  if (itv != options_list.vector_int_options.end())
    {
      if (itv->second.size() > 1)
        {
          output << "[";
          for (vector<int>::const_iterator viit = itv->second.begin();
               viit != itv->second.end(); viit++)
            output << *viit << ";";
          output << "];" << endl;
        }
      else if (itv->second.size() == 1)
        output << itv->second.front() << ";" << endl;
      else
        {
          cerr << "SvarStatement::writeOutput() Should not arrive here (3). Please report this to the Dynare Team." << endl;
          exit(EXIT_FAILURE);
        }
    }
  else
    output << "'ALL';" << endl;
}
