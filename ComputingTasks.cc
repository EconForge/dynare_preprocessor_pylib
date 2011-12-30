/*
 * Copyright (C) 2003-2011 Dynare Team
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

#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>

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
  output << "check(M_,options_,oo_);\n";
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
  output << "info = dyn_forecast(var_list_,'simul');" << endl;
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

DiscretionaryPolicyStatement::DiscretionaryPolicyStatement(const SymbolList &symbol_list_arg,
							   const OptionsList &options_list_arg) :
  symbol_list(symbol_list_arg),
  options_list(options_list_arg)
{
}

void
DiscretionaryPolicyStatement::checkPass(ModFileStructure &mod_file_struct)
{
  mod_file_struct.discretionary_policy_present = true;

  /* Fill in option_order of mod_file_struct
     Since discretionary policy needs one further order of derivation (for example, for 1st order
     approximation, it needs 2nd derivatives), we add 1 to the order declared by user */
  OptionsList::num_options_t::const_iterator it = options_list.num_options.find("order");
  if (it != options_list.num_options.end())
    {
      int order = atoi(it->second.c_str());
      if (order > 1)
        {
          cerr << "ERROR: discretionary_policy: order > 1 is not yet implemented" << endl;
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
DiscretionaryPolicyStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output);
  symbol_list.writeOutput("var_list_", output);
  output << "discretionary_policy(var_list_);\n";
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

  it = options_list.num_options.find("dsge_var");
  if (it != options_list.num_options.end())
    // Ensure that irf_shocks & dsge_var have not both been passed
    if (options_list.symbol_list_options.find("irf_shocks") != options_list.symbol_list_options.end())
      {
        cerr << "The irf_shocks and dsge_var options may not both be passed to estimation." << endl;
        exit(EXIT_FAILURE);
      }
    else
      // Fill in mod_file_struct.dsge_var_calibrated
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
    if (mod_file_struct.dsge_var_calibrated.empty()
        && !mod_file_struct.dsge_var_estimated)
      {
        cerr << "ERROR: The estimation statement requires a dsge_var option to be passed "
             << "if the dsge_varlag option is passed." << endl;
        exit(EXIT_FAILURE);
      }

  if (!mod_file_struct.dsge_var_calibrated.empty()
      && mod_file_struct.dsge_var_estimated)
    {
      cerr << "ERROR: An estimation statement cannot take more than one dsge_var option." << endl;
      exit(EXIT_FAILURE);
    }

  if (options_list.string_options.find("datafile") == options_list.string_options.end() &&
      !mod_file_struct.estimation_data_statement_present)
    {
      cerr << "ERROR: The estimation statement requires a data file to be supplied "
           << "either from the data statement or from the deprecated option datafile." << endl;
      exit(EXIT_FAILURE);
    }

  if (options_list.string_options.find("datafile") != options_list.string_options.end())
    cerr << "WARNING: The datafile option of estimation has been deprecated. "
         << "Use the data command instead." << endl;

  if (options_list.string_options.find("xls_sheet") != options_list.string_options.end())
    cerr << "WARNING: The xls_sheet option of estimation has been deprecated. "
         << "Use the data command instead." << endl;

  if (options_list.string_options.find("xls_range") != options_list.string_options.end())
    cerr << "WARNING: The xls_range option of estimation has been deprecated. "
         << "Use the data command instead." << endl;

  if (options_list.num_options.find("first_obs") != options_list.num_options.end())
    cerr << "WARNING: The first_obs option of estimation has been deprecated. "
         << "Use the data command instead." << endl;
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

UnitRootVarsStatement::UnitRootVarsStatement(void)
{
}

void
UnitRootVarsStatement::writeOutput(ostream &output, const string &basename) const
{
  output << "options_.diffuse_filter = 1;" << endl
	 << "options_.steadystate.nocheck = 1;" << endl;
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
            if (it->mean->eval(eval_context_t()) == 0.5
                && it->std->eval(eval_context_t()) == 0.5)
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

OsrParamsStatement::OsrParamsStatement(const SymbolList &symbol_list_arg) :
  symbol_list(symbol_list_arg)
{
}

void
OsrParamsStatement::checkPass(ModFileStructure &mod_file_struct)
{
  mod_file_struct.osr_params_present = true;
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
OptimWeightsStatement::checkPass(ModFileStructure &mod_file_struct)
{
  mod_file_struct.optim_weights_present = true;
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
  mod_file_struct.planner_objective_present = true;
}

StaticModel *
PlannerObjectiveStatement::getPlannerObjective() const
{
  return model_tree;
}

void
PlannerObjectiveStatement::computingPass()
{
  model_tree->computingPass(eval_context_t(), false, true, false, false);
}

void
PlannerObjectiveStatement::writeOutput(ostream &output, const string &basename) const
{
  model_tree->writeStaticFile(basename + "_objective", false, false, false);
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
  output << "sbvar(M_,options_);" << endl;
}

MSSBVAREstimationStatement::MSSBVAREstimationStatement(const OptionsList &options_list_arg) :
  options_list(options_list_arg)
{
}

void
MSSBVAREstimationStatement::checkPass(ModFileStructure &mod_file_struct)
{
  mod_file_struct.bvar_present = true;
}

void
MSSBVAREstimationStatement::writeOutput(ostream &output, const string &basename) const
{
  output << "options_ = initialize_ms_sbvar_options(M_, options_);" << endl;
  options_list.writeOutput(output);
  output << "[options_, oo_] = ms_estimation(M_, options_, oo_);" << endl;
}

MSSBVARSimulationStatement::MSSBVARSimulationStatement(const OptionsList &options_list_arg) :
  options_list(options_list_arg)
{
}

void
MSSBVARSimulationStatement::checkPass(ModFileStructure &mod_file_struct)
{
  mod_file_struct.bvar_present = true;
}

void
MSSBVARSimulationStatement::writeOutput(ostream &output, const string &basename) const
{
  output << "options_ = initialize_ms_sbvar_options(M_, options_);" << endl;
  options_list.writeOutput(output);

  // Redeclare drop option if necessary
  OptionsList::num_options_t::const_iterator mh_replic_it = options_list.num_options.find("ms.mh_replic");
  OptionsList::num_options_t::const_iterator thinning_factor_it = options_list.num_options.find("ms.thinning_factor");
  OptionsList::num_options_t::const_iterator drop_it = options_list.num_options.find("ms.drop");
  if (mh_replic_it != options_list.num_options.end() || thinning_factor_it != options_list.num_options.end())
    if (drop_it == options_list.num_options.end())
      output << "options_.ms.drop = 0.1*options_.ms.mh_replic*options_.ms.thinning_factor;" << endl;

  output << "[options_, oo_] = ms_simulation(M_, options_, oo_);" << endl;
}

MSSBVARComputeMDDStatement::MSSBVARComputeMDDStatement(const OptionsList &options_list_arg) :
  options_list(options_list_arg)
{
}

void
MSSBVARComputeMDDStatement::checkPass(ModFileStructure &mod_file_struct)
{
  mod_file_struct.bvar_present = true;
}

void
MSSBVARComputeMDDStatement::writeOutput(ostream &output, const string &basename) const
{
  output << "options_ = initialize_ms_sbvar_options(M_, options_);" << endl;
  options_list.writeOutput(output);
  output << "[options_, oo_] = ms_compute_mdd(M_, options_, oo_);" << endl;
}

MSSBVARComputeProbabilitiesStatement::MSSBVARComputeProbabilitiesStatement(const OptionsList &options_list_arg) :
  options_list(options_list_arg)
{
}

void
MSSBVARComputeProbabilitiesStatement::checkPass(ModFileStructure &mod_file_struct)
{
  mod_file_struct.bvar_present = true;

  if (options_list.num_options.find("ms.real_time_smoothed_probabilities") != options_list.num_options.end())
    if (options_list.num_options.find("ms.filtered_probabilities") != options_list.num_options.end())
      {
        cerr << "ERROR: You may only pass one of real_time_smoothed "
             << "and filtered_probabilities to ms_compute_probabilities." << endl;
        exit(EXIT_FAILURE);
      }
}

void
MSSBVARComputeProbabilitiesStatement::writeOutput(ostream &output, const string &basename) const
{
  output << "options_ = initialize_ms_sbvar_options(M_, options_);" << endl;
  options_list.writeOutput(output);
  output << "[options_, oo_] = ms_compute_probabilities(M_, options_, oo_);" << endl;
}

MSSBVARIrfStatement::MSSBVARIrfStatement(const SymbolList &symbol_list_arg,
					 const OptionsList &options_list_arg) :
  symbol_list(symbol_list_arg),
  options_list(options_list_arg)
{
}

void
MSSBVARIrfStatement::checkPass(ModFileStructure &mod_file_struct)
{
  mod_file_struct.bvar_present = true;
}

void
MSSBVARIrfStatement::writeOutput(ostream &output, const string &basename) const
{
  output << "options_ = initialize_ms_sbvar_options(M_, options_);" << endl;
  symbol_list.writeOutput("var_list_", output);
  options_list.writeOutput(output);
  output << "[options_, oo_] = ms_irf(var_list_,M_, options_, oo_);" << endl;
}

MSSBVARForecastStatement::MSSBVARForecastStatement(const OptionsList &options_list_arg) :
  options_list(options_list_arg)
{
}

void
MSSBVARForecastStatement::checkPass(ModFileStructure &mod_file_struct)
{
  mod_file_struct.bvar_present = true;
}

void
MSSBVARForecastStatement::writeOutput(ostream &output, const string &basename) const
{
  output << "options_ = initialize_ms_sbvar_options(M_, options_);" << endl;
  options_list.writeOutput(output);
  output << "[options_, oo_] = ms_forecast(M_, options_, oo_);" << endl;
}

MSSBVARVarianceDecompositionStatement::MSSBVARVarianceDecompositionStatement(const OptionsList &options_list_arg) :
  options_list(options_list_arg)
{
}

void
MSSBVARVarianceDecompositionStatement::checkPass(ModFileStructure &mod_file_struct)
{
  mod_file_struct.bvar_present = true;
}

void
MSSBVARVarianceDecompositionStatement::writeOutput(ostream &output, const string &basename) const
{
  output << "options_ = initialize_ms_sbvar_options(M_, options_);" << endl;
  options_list.writeOutput(output);
  output << "[options_, oo_] = ms_variance_decomposition(M_, options_, oo_);" << endl;
}

IdentificationStatement::IdentificationStatement(const OptionsList &options_list_arg)
{
  options_list = options_list_arg;
  if (options_list.num_options.find("max_dim_cova_group") != options_list.num_options.end())
    if (atoi(options_list.num_options["max_dim_cova_group"].c_str()) == 0)
      {
        cerr << "ERROR: The max_dim_cova_group option to identification only accepts integers > 0." << endl;
        exit(EXIT_FAILURE);
      }
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

SvarIdentificationStatement::SvarIdentificationStatement(const svar_identification_restrictions_t &restrictions_arg,
                                                         const bool &upper_cholesky_present_arg,
                                                         const bool &lower_cholesky_present_arg,
                                                         const bool &constants_exclusion_present_arg,
                                                         const SymbolTable &symbol_table_arg) :
  restrictions(restrictions_arg),
  upper_cholesky_present(upper_cholesky_present_arg),
  lower_cholesky_present(lower_cholesky_present_arg),
  constants_exclusion_present(constants_exclusion_present_arg),
  symbol_table(symbol_table_arg)
{
}

int
SvarIdentificationStatement::getMaxLag() const
{
  int max_lag = 0;
  for (svar_identification_restrictions_t::const_iterator it = restrictions.begin(); it != restrictions.end(); it++)
    if (it->lag > max_lag)
      max_lag = it->lag;

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

  if (upper_cholesky_present && lower_cholesky_present)
    {
      cerr << "ERROR: Within the svar_identification statement, you may only have one of "
           << "upper_cholesky and lower_cholesky." << endl;
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

  if (constants_exclusion_present)
    output << "options_.ms.constants_exclusion=1;" << endl;

  if (!upper_cholesky_present && !lower_cholesky_present)
    {
      int n = symbol_table.endo_nbr();
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
      output << "options_.ms.Qi = cell(" << n << ",1);" << endl;
      output << "options_.ms.Ri = cell(" << n << ",1);" << endl;

      // vector<int> rows(n);
      // fill(rows.begin(),rows.end(),1);

      for (svar_identification_restrictions_t::const_iterator it = restrictions.begin(); it != restrictions.end(); it++)
        {
	  if (it->lag == 0)
	    {
	      output << "options_.ms.Qi{" << it->equation << "}(" << it->restriction_nbr << ", " << it->variable + 1 << ") = ";
	      it->value->writeOutput(output);
	      output << ";" << endl;
	    }
	  else if (it->lag > 0)
	    {
	      int col = (it->lag-1)*n+it->variable+1;
	      if (col > k)
                     {
                       cerr << "ERROR: lag =" << it->lag << ", num endog vars = " << n << "current endog var index = " << it->variable << ". Index "
                            << "out of bounds. If the above does not represent a logical error, please report this to the Dyanre Team." << endl;
                     }
	      output << "options_.ms.Ri{" << it->equation << "}(" << it->restriction_nbr << ", " << col << ") = ";
	      it->value->writeOutput(output);
	      output << ";" << endl;
	    }
	  else
	    {
	      cerr << "SvarIdentificationStatement::writeOutput() Should not arrive here (3). Please report this to the Dynare Team." << endl;
	      exit(EXIT_FAILURE);
	    }
 
          // for (unsigned int h = 0; h < it->second.size(); h++)
          //   {
          //     int j = it->second.at(h) + 1;
          //     int i = it->first.second;
	  //     int lag = it->first.first;
          //     if (j < 1 || j > n || (int) h+1 > n || i < 1)
          //       {
          //         cerr << "SvarIdentificationStatement::writeOutput() Should not arrive here (2). Please report this to the Dynare Team." << endl;
          //         exit(EXIT_FAILURE);
          //       }
          //     if (i > n)
          //       {
          //         cerr << "ERROR: equation number " << i << " is greater than the number of endogenous variables, " << n << "." << endl;
          //         exit(EXIT_FAILURE);
          //       }

          //     if (lag == 0)
          //       output << "options_.ms.Qi{" << i << "}(" << h+1 << ", " << j << ") = 1;" << endl;
          //     else if (lag > 0)
          //       {
          //         if ((lag-1)*n+j > k)
          //           {
          //             cerr << "ERROR: lag =" << lag << ", num endog vars = " << n << "current endog var index = " << j << ". Index "
          //                  << "out of bounds. If the above does not represent a logical error, please report this to the Dyanre Team." << endl;
          //           }
          //         output << "options_.ms.Ri{" << i << "}(" << rows[i-1] << ", " << (lag-1)*n+j << ") = 1;" << endl;
	  // 	  rows[i-1]++;
          //       }
          //     else
          //       {
          //         cerr << "SvarIdentificationStatement::writeOutput() Should not arrive here (3). Please report this to the Dynare Team." << endl;
          //         exit(EXIT_FAILURE);
          //       }
          //   }
        }
    }
}

MarkovSwitchingStatement::MarkovSwitchingStatement(const OptionsList &options_list_arg) :
  options_list(options_list_arg)
{
  OptionsList::num_options_t::const_iterator it_num = options_list.num_options.find("ms.restrictions");
  if (it_num != options_list.num_options.end())
    {
      using namespace boost;
      OptionsList::num_options_t::const_iterator it_num_regimes =
        options_list.num_options.find("ms.number_of_regimes");
      if (it_num_regimes ==  options_list.num_options.end())
        {
          cerr << "ERROR: should not arrive here: MarkovSwitchingStatement constructor" << endl;
          exit(EXIT_FAILURE);
        }
      int num_regimes = lexical_cast< int >(it_num_regimes->second);

      vector<string> tokenizedRestrictions;
      split(tokenizedRestrictions, it_num->second, is_any_of("["), token_compress_on);
      for (vector<string>::iterator it = tokenizedRestrictions.begin();
            it != tokenizedRestrictions.end(); it++ )
        if (it->size() > 0)
          {
            vector<string> restriction;
            split(restriction, *it, is_any_of("], "));
            for (vector<string>::iterator it1 = restriction.begin();
                 it1 != restriction.end(); )
              if (it1->empty())
                restriction.erase(it1);
              else
                it1++;

            if (restriction.size() != 3)
              {
                cerr << "ERROR: restrictions in the subsample statement must be specified in the form "
                     << "[current_period_regime, next_period_regime, transition_probability]" << endl;
                exit(EXIT_FAILURE);
              }

            try
              {
                int from_regime = lexical_cast< int >(restriction[0]);
                int to_regime = lexical_cast< int >(restriction[1]);
                if (from_regime > num_regimes || to_regime > num_regimes)
                  {
                    cerr << "ERROR: the regimes specified in the restrictions option must be "
                         << "<= the number of regimes specified in the number_of_regimes option" << endl;
                    exit(EXIT_FAILURE);
                  }

                if (restriction_map.find(make_pair(from_regime, to_regime)) !=
                    restriction_map.end())
                  {
                    cerr << "ERROR: two restrictions were given for: " << from_regime << ", "
                         << to_regime << endl;
                    exit(EXIT_FAILURE);
                  }

                double transition_probability = lexical_cast< double >(restriction[2]);
                if (transition_probability > 1.0)
                  {
                    cerr << "ERROR: the transition probability, " << transition_probability
                         << " must be less than 1" << endl;
                    exit(EXIT_FAILURE);
                  }
                restriction_map[make_pair(from_regime, to_regime)] = transition_probability;
              }
            catch (const bad_lexical_cast &)
              {
                cerr << "ERROR: The first two arguments for a restriction must be integers "
                     << "specifying the regime and the last must be a double specifying the "
                     << "transition probability. You wrote [" << *it << endl;
                exit(EXIT_FAILURE);
              }
          }
    }
}

void
MarkovSwitchingStatement::checkPass(ModFileStructure &mod_file_struct)
{
  OptionsList::num_options_t::const_iterator it_num = options_list.num_options.find("ms.restrictions");
  if (it_num != options_list.num_options.end())
    {
      using namespace boost;
      OptionsList::num_options_t::const_iterator it_num_regimes =
        options_list.num_options.find("ms.number_of_regimes");
      int num_regimes = lexical_cast< int >(it_num_regimes->second);
      vector<double> col_trans_prob_sum (num_regimes, 0);
      vector<double> row_trans_prob_sum (num_regimes, 0);
      vector<bool> all_restrictions_in_row (num_regimes, true);
      vector<bool> all_restrictions_in_col (num_regimes, true);
      for (int row=0; row<num_regimes; row++)
        for (int col=0; col<num_regimes; col++)
          if (restriction_map.find(make_pair(row+1, col+1)) != restriction_map.end())
            {
              row_trans_prob_sum[row] += restriction_map[make_pair(row+1, col+1)];
              col_trans_prob_sum[col] += restriction_map[make_pair(row+1, col+1)];
            }
          else
            {
              all_restrictions_in_row[row] = false;
              all_restrictions_in_col[col] = false;
            }

      for (int i=0; i<num_regimes; i++)
        {
          if (all_restrictions_in_row[i])
          {
            if (row_trans_prob_sum[i] != 1.0)
              {
                cerr << "ERROR: When all transitions probabilities are specified for a certain "
                     << "regime, they must sum to 1" << endl;
                exit(EXIT_FAILURE);
              }
          }
        else
          if (row_trans_prob_sum[i] >= 1.0)
            {
              cerr << "ERROR: When transition probabilites are not specified for every regime, "
                   << "their sum must be < 1" << endl;
              exit(EXIT_FAILURE);
            }

        if (all_restrictions_in_col[i])
          {
            if (col_trans_prob_sum[i] != 1.0)
              {
                cerr << "ERROR: When all transitions probabilities are specified for a certain "
                     << "regime, they must sum to 1" << endl;
                exit(EXIT_FAILURE);
              }
          }
        else
          if (col_trans_prob_sum[i] >= 1.0)
            {
              cerr << "ERROR: When transition probabilites are not specified for every regime, "
                   << "their sum must be < 1" << endl;
              exit(EXIT_FAILURE);
            }
        }
    }
}

void
MarkovSwitchingStatement::writeOutput(ostream &output, const string &basename) const
{
  bool isDurationAVec = true;
  string infStr("Inf");
  OptionsList::num_options_t::const_iterator itChain, itNOR, itDuration;
  map<pair<int, int>, double >::const_iterator itR;

  itChain = options_list.num_options.find("ms.chain");
  if (itChain == options_list.num_options.end())
    {
      cerr << "MarkovSwitchingStatement::writeOutput() Should not arrive here (1). "
           << "Please report this to the Dynare Team." << endl;
      exit(EXIT_FAILURE);
    }

  itDuration = options_list.num_options.find("ms.duration");
  if (itDuration == options_list.num_options.end())
    {
      cerr << "MarkovSwitchingStatement::writeOutput() Should not arrive here (2). "
           << "Please report this to the Dynare Team." << endl;
      exit(EXIT_FAILURE);
    }
  else if (atof(itDuration->second.c_str()) || infStr.compare(itDuration->second) == 0)
    isDurationAVec = false;
  output << "options_.ms.duration = " << itDuration->second << ";" << endl;

  itNOR = options_list.num_options.find("ms.number_of_regimes");
  if (itNOR != options_list.num_options.end())
    for (int i = 0; i < atoi(itNOR->second.c_str()); i++)
      {
        output << "options_.ms.ms_chain(" << itChain->second << ").regime("
               << i+1 << ").duration = options_.ms.duration";
        if (isDurationAVec)
          output << "(" << i+1 << ")";
        output << ";" << endl;
      }
  else
    {
      cerr << "MarkovSwitchingStatement::writeOutput() Should not arrive here (3). "
           << "Please report this to the Dynare Team." << endl;
      exit(EXIT_FAILURE);
    }

  int restrictions_index = 0;
  for (itR=restriction_map.begin(); itR != restriction_map.end(); itR++)
    output << "options_.ms.ms_chain(" << itChain->second << ").restrictions("
           << ++restrictions_index << ") = {[" << itR->first.first << ", "
           << itR->first.second << ", " << itR->second << "]};" << endl;
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

SetTimeStatement::SetTimeStatement(const OptionsList &options_list_arg) :
  options_list(options_list_arg)
{
}

void
SetTimeStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output);
}

EstimationDataStatement::EstimationDataStatement(const OptionsList &options_list_arg) :
  options_list(options_list_arg)
{
}

void
EstimationDataStatement::checkPass(ModFileStructure &mod_file_struct)
{
  mod_file_struct.estimation_data_statement_present = true;

  OptionsList::num_options_t::const_iterator it = options_list.num_options.find("nobs");
  if (it != options_list.num_options.end())
    if (atoi(it->second.c_str()) <= 0)
      {
        cerr << "ERROR: The nobs option of the data statement only accepts positive integers." << endl;
        exit(EXIT_FAILURE);
      }

  if (options_list.string_options.find("file") == options_list.string_options.end())
    {
      cerr << "ERROR: The file option must be passed to the data statement." << endl;
      exit(EXIT_FAILURE);
    }
}

void
EstimationDataStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output, "options_.dataset");
  if (options_list.date_options.find("first_obs") == options_list.date_options.end())
    output << "options_.dataset.firstobs = options_.initial_period;" << endl;
}

BasicPriorStatement::~BasicPriorStatement()
{
}

BasicPriorStatement::BasicPriorStatement(const string &name_arg,
                                         const expr_t &variance_arg,
                                         const OptionsList &options_list_arg) :
  name(name_arg),
  variance(variance_arg),
  options_list(options_list_arg),
  first_statement_encountered(false)
{
}

void
BasicPriorStatement::checkPass(ModFileStructure &mod_file_struct)
{
  if (options_list.num_options.find("shape") == options_list.num_options.end())
    {
      cerr << "ERROR: You must pass the shape option to the prior statement." << endl;
      exit(EXIT_FAILURE);
    }
  if (options_list.num_options.find("date1") != options_list.num_options.end() ||
      options_list.num_options.find("date2") != options_list.num_options.end())
    if (options_list.num_options.find("date1") == options_list.num_options.end() ||
        options_list.num_options.find("date2") == options_list.num_options.end())
      {
        cerr << "ERROR: PriorStatement::checkPass(1). Should not arrive here. "
             << "Please inform Dynare Team." << endl;
        exit(EXIT_FAILURE);
      }
}

void
BasicPriorStatement::get_base_name(const SymbolType symb_type, string &lhs_field) const
{
  if (symb_type == eExogenous || symb_type == eExogenousDet)
    lhs_field = "structural_innovation";
  else
    lhs_field = "measurement_error";
}

void
BasicPriorStatement::writePriorIndex(ostream &output, const string &lhs_field) const
{
  if (first_statement_encountered)
    output << "prior_indx = 1;" << endl;
  else
    output << "prior_indx = size(estimation_info" << lhs_field << "_index, 2) + 1;" << endl;
}

void
BasicPriorStatement::writeVarianceOption(ostream &output, const string &lhs_field) const
{
  if (variance)
    {
      output << "estimation_info" << lhs_field << "(prior_indx).variance = ";
      variance->writeOutput(output);
      output << ";" << endl;
    }
}

void
BasicPriorStatement::writeOutputHelper(ostream &output, const string &field, const string &lhs_field) const
{

  OptionsList::num_options_t::const_iterator itn = options_list.num_options.find(field);
  if (itn != options_list.num_options.end())
    output << "estimation_info" << lhs_field << "(prior_indx)." << field
           << " = " << itn->second << ";" << endl;

  OptionsList::date_options_t::const_iterator itd = options_list.date_options.find(field);
  if (itd != options_list.date_options.end())
    output << "estimation_info" << lhs_field << "(prior_indx)." << field
           << " = '" << itd->second << "';" << endl;
}

PriorStatement::PriorStatement(const string &name_arg,
                               const expr_t &variance_arg,
                               const OptionsList &options_list_arg) :
  BasicPriorStatement(name_arg, variance_arg, options_list_arg)
{
}

void
PriorStatement::checkPass(ModFileStructure &mod_file_struct)
{
  BasicPriorStatement::checkPass(mod_file_struct);
  if (!mod_file_struct.prior_statement_present)
    first_statement_encountered = true;
  mod_file_struct.prior_statement_present = true;
}

void
PriorStatement::writeOutput(ostream &output, const string &basename) const
{
  string lhs_field = ".prior";

  writePriorIndex(output, lhs_field);
  output << "estimation_info" << lhs_field << "_index(prior_indx) = {'" << name << "'};" << endl
         << "estimation_info" << lhs_field <<"(prior_indx).name = '" << name << "';" << endl;

  writeOutputHelper(output, "mean", lhs_field);
  writeOutputHelper(output, "mode", lhs_field);
  writeOutputHelper(output, "stdev", lhs_field);
  writeOutputHelper(output, "shape", lhs_field);
  writeOutputHelper(output, "shift", lhs_field);
  writeOutputHelper(output, "date1", lhs_field);
  writeOutputHelper(output, "date2", lhs_field);
  writeOutputHelper(output, "domain", lhs_field);
  writeOutputHelper(output, "interval", lhs_field);
  writeVarianceOption(output, lhs_field);
}

StdPriorStatement::StdPriorStatement(const string &name_arg,
                                     const expr_t &variance_arg,
                                     const OptionsList &options_list_arg,
                                     const SymbolTable &symbol_table_arg ) :
  BasicPriorStatement(name_arg, variance_arg, options_list_arg),
  symbol_table(symbol_table_arg)
{
}

void
StdPriorStatement::checkPass(ModFileStructure &mod_file_struct)
{
  BasicPriorStatement::checkPass(mod_file_struct);
  if (!mod_file_struct.std_prior_statement_present)
    first_statement_encountered = true;
  mod_file_struct.std_prior_statement_present = true;
}

void
StdPriorStatement::writeOutput(ostream &output, const string &basename) const
{
  string lhs_field;
  get_base_name(symbol_table.getType(name), lhs_field);
  lhs_field = "." + lhs_field + ".prior";

  writePriorIndex(output, lhs_field);
  output << "estimation_info" << lhs_field << "_index(prior_indx) = {'" << name << "'};" << endl;
  output << "estimation_info" << lhs_field << "(prior_indx).name = '" << name << "';" << endl;

  writeOutputHelper(output, "mean", lhs_field);
  writeOutputHelper(output, "mode", lhs_field);
  writeOutputHelper(output, "stdev", lhs_field);
  writeOutputHelper(output, "shape", lhs_field);
  writeOutputHelper(output, "shift", lhs_field);
  writeOutputHelper(output, "domain", lhs_field);
  writeOutputHelper(output, "interval", lhs_field);
  writeVarianceOption(output, lhs_field);
}

CorrPriorStatement::CorrPriorStatement(const string &name_arg1, const string &name_arg2,
                                       const expr_t &variance_arg,
                                       const OptionsList &options_list_arg,
                                       const SymbolTable &symbol_table_arg ) :
  BasicPriorStatement(name_arg1, variance_arg, options_list_arg),
  name1(name_arg2),
  symbol_table(symbol_table_arg)
{
}

void
CorrPriorStatement::checkPass(ModFileStructure &mod_file_struct)
{
  BasicPriorStatement::checkPass(mod_file_struct);
  if (symbol_table.getType(name) != symbol_table.getType(name1))
    {
      cerr << "ERROR: In the corr(A,B).prior statement, A and B must be of the same type. "
           << "In your case, " << name << " and " << name1 << " are of different "
           << "types." << endl;
      exit(EXIT_FAILURE);
    }
  if (!mod_file_struct.corr_prior_statement_present)
    first_statement_encountered = true;
  mod_file_struct.corr_prior_statement_present = true;
}

void
CorrPriorStatement::writeOutput(ostream &output, const string &basename) const
{
  string lhs_field;
  get_base_name(symbol_table.getType(name), lhs_field);
  lhs_field = "." + lhs_field + "_corr.prior";

  writePriorIndex(output, lhs_field);
  output << "estimation_info" << lhs_field << "_index(prior_indx) = {'" << name << "_" << name1 << "'};" << endl;
  output << "estimation_info" << lhs_field << "(prior_indx).name1 = '" << name << "';" << endl;
  output << "estimation_info" << lhs_field << "(prior_indx).name2 = '" << name1 << "';" << endl;

  writeOutputHelper(output, "mean", lhs_field);
  writeOutputHelper(output, "mode", lhs_field);
  writeOutputHelper(output, "stdev", lhs_field);
  writeOutputHelper(output, "shape", lhs_field);
  writeOutputHelper(output, "shift", lhs_field);
  writeOutputHelper(output, "domain", lhs_field);
  writeOutputHelper(output, "interval", lhs_field);
  writeVarianceOption(output, lhs_field);
}

BasicOptionsStatement::~BasicOptionsStatement()
{
}

BasicOptionsStatement::BasicOptionsStatement(const string &name_arg,
                                         const OptionsList &options_list_arg) :
  name(name_arg),
  options_list(options_list_arg),
  first_statement_encountered(false)
{
}

void
BasicOptionsStatement::checkPass(ModFileStructure &mod_file_struct)
{
  if (options_list.num_options.find("date1") != options_list.num_options.end() ||
      options_list.num_options.find("date2") != options_list.num_options.end())
    if (options_list.num_options.find("date1") == options_list.num_options.end() ||
        options_list.num_options.find("date2") == options_list.num_options.end())
      {
        cerr << "ERROR: OptionsStatement::checkPass(1). Should not arrive here. "
             << "Please inform Dynare Team." << endl;
        exit(EXIT_FAILURE);
      }
}

void
BasicOptionsStatement::writeOptionsIndex(ostream &output, const string &lhs_field) const
{
  if (first_statement_encountered)
    output << "options_indx = 1;" << endl;
  else
    output << "options_indx = size(estimation_info" << lhs_field << "_index, 2) + 1;" << endl;
}

void
BasicOptionsStatement::get_base_name(const SymbolType symb_type, string &lhs_field) const
{
  if (symb_type == eExogenous || symb_type == eExogenousDet)
    lhs_field = "structural_innovation";
  else
    lhs_field = "measurement_error";
}

void
BasicOptionsStatement::writeOutputHelper(ostream &output, const string &field, const string &lhs_field) const
{
  OptionsList::num_options_t::const_iterator itn = options_list.num_options.find(field);
  if (itn != options_list.num_options.end())
    output << "estimation_info" << lhs_field << "(options_indx)." << field
           << " = " << itn->second << ";" << endl;

  OptionsList::date_options_t::const_iterator itd = options_list.date_options.find(field);
  if (itd != options_list.date_options.end())
    output << "estimation_info" << lhs_field << "(options_indx)." << field
           << " = '" << itd->second << "';" << endl;
}

OptionsStatement::OptionsStatement(const string &name_arg,
                                   const OptionsList &options_list_arg) :
  BasicOptionsStatement(name_arg, options_list_arg)
{
}

void
OptionsStatement::checkPass(ModFileStructure &mod_file_struct)
{
  BasicOptionsStatement::checkPass(mod_file_struct);
  if (!mod_file_struct.options_statement_present)
    first_statement_encountered = true;
  mod_file_struct.options_statement_present = true;
}

void
OptionsStatement::writeOutput(ostream &output, const string &basename) const
{
  string lhs_field = ".options";

  writeOptionsIndex(output, lhs_field);
  output << "estimation_info" << lhs_field <<"_index(options_indx) = {'" << name << "'};" << endl
         << "estimation_info" << lhs_field << "(options_indx).name = '" << name << "';" << endl;

  writeOutputHelper(output, "init", lhs_field);
  writeOutputHelper(output, "bounds", lhs_field);
  writeOutputHelper(output, "jscale", lhs_field);
  writeOutputHelper(output, "date1", lhs_field);
  writeOutputHelper(output, "date2", lhs_field);
}

StdOptionsStatement::StdOptionsStatement(const string &name_arg,
                                         const OptionsList &options_list_arg,
                                         const SymbolTable &symbol_table_arg ) :
  BasicOptionsStatement(name_arg, options_list_arg),
  symbol_table(symbol_table_arg)
{
}

void
StdOptionsStatement::checkPass(ModFileStructure &mod_file_struct)
{
  BasicOptionsStatement::checkPass(mod_file_struct);
  if (!mod_file_struct.std_options_statement_present)
    first_statement_encountered = true;
  mod_file_struct.std_options_statement_present = true;
}

void
StdOptionsStatement::writeOutput(ostream &output, const string &basename) const
{
  string lhs_field;
  get_base_name(symbol_table.getType(name), lhs_field);
  lhs_field = "." + lhs_field + ".options";

  writeOptionsIndex(output, lhs_field);
  output << "estimation_info" << lhs_field << "_index(options_indx) = {'" << name << "'};" << endl;
  output << "estimation_info" << lhs_field << "(options_indx).name = '" << name << "';" << endl;

  writeOutputHelper(output, "init", lhs_field);
  writeOutputHelper(output, "bounds", lhs_field);
  writeOutputHelper(output, "jscale", lhs_field);
  writeOutputHelper(output, "date1", lhs_field);
  writeOutputHelper(output, "date2", lhs_field);
}

CorrOptionsStatement::CorrOptionsStatement(const string &name_arg1, const string &name_arg2,
                                           const OptionsList &options_list_arg,
                                           const SymbolTable &symbol_table_arg ) :
  BasicOptionsStatement(name_arg1, options_list_arg),
  name1(name_arg2),
  symbol_table(symbol_table_arg)
{
}

void
CorrOptionsStatement::checkPass(ModFileStructure &mod_file_struct)
{
  if (symbol_table.getType(name) != symbol_table.getType(name1))
    {
      cerr << "ERROR: In the corr(A,B).options statement, A and B must be of the same type. "
           << "In your case, " << name << " and " << name1 << " are of different "
           << "types." << endl;
      exit(EXIT_FAILURE);
    }
  if (!mod_file_struct.corr_options_statement_present)
    first_statement_encountered = true;
  mod_file_struct.corr_prior_statement_present = true;
}

void
CorrOptionsStatement::writeOutput(ostream &output, const string &basename) const
{
  string lhs_field;
  get_base_name(symbol_table.getType(name), lhs_field);
  lhs_field = "." + lhs_field + "_corr.options";

  writeOptionsIndex(output, lhs_field);
  output << "estimation_info" << lhs_field << "_index(options_indx) = {'" << name << "_" << name1 << "'};" << endl;
  lhs_field += ".";
  output << "estimation_info" << lhs_field << "(options_indx).name1 = '" << name << "';" << endl;
  output << "estimation_info" << lhs_field << "(options_indx).name2 = '" << name1 << "';" << endl;

  writeOutputHelper(output, "init", lhs_field);
  writeOutputHelper(output, "bounds", lhs_field);
  writeOutputHelper(output, "jscale", lhs_field);
  writeOutputHelper(output, "date1", lhs_field);
  writeOutputHelper(output, "date2", lhs_field);
}
