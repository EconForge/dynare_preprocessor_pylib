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

#include <cstdlib>
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
SteadyStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output);
  output << "steady;\n";
}

SteadySparseStatement::SteadySparseStatement(const OptionsList &options_list_arg) :
  options_list(options_list_arg)
{
}

void
SteadySparseStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output);
  output << basename << "_static;\n";
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

Model_InfoStatement::Model_InfoStatement(const OptionsList &options_list_arg) :
  options_list(options_list_arg)
{
}

void Model_InfoStatement::checkPass(ModFileStructure &mod_file_struct)
{
  //mod_file_struct.model_info_present = true;
}

void Model_InfoStatement::writeOutput(ostream &output, const string &basename) const
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
  output << "simul(oo_.dr);\n";
}

SimulSparseStatement::SimulSparseStatement(const OptionsList &options_list_arg,
                                           int compiler_arg, int mode_arg) :
  options_list(options_list_arg),
  compiler(compiler_arg),
  mode(mode_arg)
{
}

void
SimulSparseStatement::checkPass(ModFileStructure &mod_file_struct)
{
  mod_file_struct.simul_present = true;
}

void
SimulSparseStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output);
  output << "if (~ options_.initval_file) & (size(oo_.endo_simul,2)<options_.periods)\n";
  output << "  if ~isfield(options_,'datafile')\n";
  output << "    make_y_;\n";
  output << "    make_ex_;\n";
  output << "  else\n";
  output << "    read_data_;\n";
  output << "  end\n";
  output << "end\n";
  if(mode==eSparseDLLMode)
    {
      if(compiler!=NO_COMPILE)
        {
          output << "disp('compiling...');\n";
          output << "t0=clock;\n";
          if (compiler == 0)
            output << "mex " << basename << "_dynamic.c;\n";
          else
            output << "mex " << basename << "_dynamic.cc;\n";
          output << "disp(['compiling time: ' num2str(etime(clock,t0))]);\n";
          output << "oo_.endo_simul=" << basename << "_dynamic;\n";
        }
      else
        {
          output << "oo_.endo_simul=simulate;\n";
          output << "clear simulate.dll;\n";
        }
    }
  else
    {
      //output << "oo_.endo_simul=" << basename << "_dynamic();\n";
      output << basename << "_dynamic;\n";
    }
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
  mod_file_struct.stoch_simul_or_similar_present = true;

  // Fill in option_order of mod_file_struct
  OptionsList::num_options_type::const_iterator it = options_list.num_options.find("order");
  if (it != options_list.num_options.end())
    mod_file_struct.order_option = max(mod_file_struct.order_option,atoi(it->second.c_str()));
}

void
StochSimulStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output);
  symbol_list.writeOutput("var_list_", output);
  output << "info = stoch_simul(var_list_);\n";
}

ForecastStatement::ForecastStatement(const SymbolList &symbol_list_arg,
                                         const OptionsList &options_list_arg) :
  symbol_list(symbol_list_arg),
  options_list(options_list_arg)
{
}

void
ForecastStatement::checkPass(ModFileStructure &mod_file_struct)
{
  mod_file_struct.stoch_simul_or_similar_present = true;

  // Fill in option_order of mod_file_struct
  OptionsList::num_options_type::const_iterator it = options_list.num_options.find("order");
  if (it != options_list.num_options.end())
    mod_file_struct.order_option = max(mod_file_struct.order_option,atoi(it->second.c_str()));
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
  mod_file_struct.stoch_simul_or_similar_present = true;

  /* Fill in option_order of mod_file_struct
     Since ramsey policy needs one further order of derivation (for example, for 1st order
     approximation, it needs 2nd derivatives), we add 1 to the order declared by user */
  OptionsList::num_options_type::const_iterator it = options_list.num_options.find("order");
  if (it != options_list.num_options.end())
    mod_file_struct.order_option = max(mod_file_struct.order_option,atoi(it->second.c_str()) + 1);
}

void
RamseyPolicyStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output);
  symbol_list.writeOutput("var_list_", output);
  output << "ramsey_policy(var_list_);\n";
}

EstimationStatement::EstimationStatement(const SymbolList &symbol_list_arg,
                                         const OptionsList &options_list_arg) :
  symbol_list(symbol_list_arg),
  options_list(options_list_arg)
{
}

void
EstimationStatement::checkPass(ModFileStructure &mod_file_struct)
{
  mod_file_struct.stoch_simul_or_similar_present = true;

  // Fill in option_order of mod_file_struct
  OptionsList::num_options_type::const_iterator it = options_list.num_options.find("order");
  if (it != options_list.num_options.end())
    mod_file_struct.order_option = max(mod_file_struct.order_option,atoi(it->second.c_str()));
}

void
EstimationStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output);
  symbol_list.writeOutput("var_list_", output);
  output << "dynare_estimation(var_list_);\n";
}

PriorAnalysisStatement::PriorAnalysisStatement(const SymbolList &symbol_list_arg,
                                               const OptionsList &options_list_arg) :
  symbol_list(symbol_list_arg),
  options_list(options_list_arg)
{
}

void
PriorAnalysisStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output);
  symbol_list.writeOutput("var_list_", output);
  output << "prior_analysis(var_list_);\n";
}

PosteriorAnalysisStatement::PosteriorAnalysisStatement(const SymbolList &symbol_list_arg,
                                                       const OptionsList &options_list_arg) :
  symbol_list(symbol_list_arg),
  options_list(options_list_arg)
{
}

void
PosteriorAnalysisStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output);
  symbol_list.writeOutput("var_list_", output);
  output << "posterior_analysis(var_list_);\n";
}

DynareSensitivityStatement::DynareSensitivityStatement(const OptionsList &options_list_arg) :
  options_list(options_list_arg)
{
}

void
DynareSensitivityStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output,"options_gsa");
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
  output << "options_.simul = 1;" << endl;
}

CutoffStatement::CutoffStatement(double cutoff_arg) : cutoff(cutoff_arg)
{
}

void
CutoffStatement::writeOutput(ostream &output, const string &basename) const
{
  output << "options_.cutoff = " << cutoff << ";" << endl;
}

MarkowitzStatement::MarkowitzStatement(double markowitz_arg) : markowitz(markowitz_arg)
{
}

void
MarkowitzStatement::writeOutput(ostream &output, const string &basename) const
{
  output << "options_.markowitz = " << markowitz << ";" << endl;
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
    output << "options_.dsample = " << val1 << ";" << endl;
  else
    output << "options_.dsample = [" << val1 << "; " << val2 << "];" << endl;
}

VarobsStatement::VarobsStatement(const SymbolList &symbol_list_arg) :
  symbol_list(symbol_list_arg)
{
}

void
VarobsStatement::writeOutput(ostream &output, const string &basename) const
{
  symbol_list.writeOutput("options_.varobs", output);
}

EstimatedParamsStatement::EstimatedParamsStatement(const vector<EstimationParams> &estim_params_list_arg,
                                                   const SymbolTable &symbol_table_arg) :
  estim_params_list(estim_params_list_arg),
  symbol_table(symbol_table_arg)
{
}

void
EstimatedParamsStatement::writeOutput(ostream &output, const string &basename) const
{
  output << "global estim_params_\n";
  output << "var_list_ = [];\n";
  output << "estim_params_.var_exo = [];\n";
  output << "estim_params_.var_endo = [];\n";
  output << "estim_params_.corrx = [];\n";
  output << "estim_params_.corrn = [];\n";
  output << "estim_params_.param_names = [];\n";
  output << "estim_params_.user_param_names = [];\n";
  output << "estim_params_.param_vals = [];\n";
  output << "M_.H = 0;\n";

  vector<EstimationParams>::const_iterator it;

  for(it = estim_params_list.begin(); it != estim_params_list.end(); it++)
    {
      switch(it->type)
        {
        case 1:
          if (symbol_table.getType(it->name) == eExogenous)
            output << "estim_params_.var_exo = [estim_params_.var_exo; ";
          else if (symbol_table.getType(it->name) == eEndogenous)
            output << "estim_params_.var_endo = [estim_params_.var_endo; ";
          output << symbol_table.getID(it->name)+1;
          break;
        case 2:
          output << "estim_params_.param_vals = [estim_params_.param_vals; ";
          output << symbol_table.getID(it->name)+1;
          break;
        case 3:
          if (symbol_table.getType(it->name) == eExogenous)
            output << "estim_params_.corrx = [estim_params_.corrx; ";
          else if (symbol_table.getType(it->name) == eEndogenous)
            output << "estim_params_.corrn = [estim_params_.corrn; ";
          output << symbol_table.getID(it->name)+1;
          output << " " << symbol_table.getID(it->name2)+1;
          break;
        }
      output << " " << it->init_val << " " <<  it->low_bound
             << " " << it->up_bound << " " <<  it->prior
             << " " << it->mean << " " <<  it->std
             << " " << it->p3 << " " <<  it->p4  << " " <<  it->jscale << "];\n";
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

  for(it = estim_params_list.begin(); it != estim_params_list.end(); it++)
    {
      if (it->type < 3)
        {
          if (symbol_table.getType(it->name) == eExogenous)
            {
              output << "tmp1 = find(estim_params_.var_exo(:,1)==" << symbol_table.getID(it->name)+1 << ");\n";
              output << "estim_params_.var_exo(tmp1,2) = " << it->init_val << ";\n";
            }
          else if (symbol_table.getType(it->name) == eEndogenous)
            {
              output << "tmp1 = find(estim_params_.var_endo(:,1)==" << symbol_table.getID(it->name)+1 << ");\n";
              output << "estim_params_.var_endo(tmp1,2) = " << it->init_val << ";\n";
            }
          else if (symbol_table.getType(it->name) == eParameter)
            {
              output << "tmp1 = find(estim_params_.param_vals(:,1)==" << symbol_table.getID(it->name)+1 << ");\n";
              output << "estim_params_.param_vals(tmp1,2) = " << it->init_val << ";\n";
            }
        }
      else
        {
          if (symbol_table.getType(it->name) == eExogenous)
            {
              output << "tmp1 = find((estim_params_.corrx(:,1)==" << symbol_table.getID(it->name)+1 << ")) & (estim_params_.corrx(:,2)==" << symbol_table.getID(it->name2)+1 << ");\n";
              output << "estim_params_.corrx(tmp1,3) = " << it->init_val << ";\n";
            }
          else if (symbol_table.getType(it->name) == eEndogenous)
            {
              output << "tmp1 = find((estim_params_.corrn(:,1)==" << symbol_table.getID(it->name)+1 << ")) & (estim_params_.corrn(:,2)==" << symbol_table.getID(it->name2)+1 << ";\n";
              output << "estim_params_.corrx(tmp1,3) = " << it->init_val << ";\n";
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

  for(it = estim_params_list.begin(); it != estim_params_list.end(); it++)
    {
      if (it->type < 3)
        {
          if (symbol_table.getType(it->name) == eExogenous)
            {
              output << "tmp1 = find(estim_params_.var_exo(:,1)==" << symbol_table.getID(it->name)+1 << ");\n";
              output << "estim_params_.var_exo(tmp1,3) = " << it->low_bound << ";\n";
              output << "estim_params_.var_exo(tmp1,4) = " << it->up_bound << ";\n";
            }
          else if (symbol_table.getType(it->name) == eEndogenous)
            {
              output << "tmp1 = find(estim_params_.var_endo(:,1)==" << symbol_table.getID(it->name)+1 << ");\n";
              output << "estim_params_.var_endo(tmp1,3) = " << it->low_bound << ";\n";
              output << "estim_params_.var_endo(tmp1,4) = " << it->up_bound << ";\n";
            }
          else if (symbol_table.getType(it->name) == eParameter)
            {
              output << "tmp1 = find(estim_params_.param_vals(:,1)==" << symbol_table.getID(it->name)+1 << ");\n";
              output << "estim_params_.param_vals(tmp1,3) = " << it->low_bound << ";\n";
              output << "estim_params_.param_vals(tmp1,4) = " << it->up_bound << ";\n";
            }
        }
      else
        {
          if (symbol_table.getType(it->name) == eExogenous)
            {
              output << "tmp1 = find((estim_params_.corrx(:,1)==" << symbol_table.getID(it->name)+1 << ")) & (estim_params_.corrx(:,2)==" << symbol_table.getID(it->name2)+1 << ");\n";
              output << "estim_params_.corrx(tmp1,4) = " << it->low_bound << ";\n";
              output << "estim_params_.corrx(tmp1,5) = " << it->up_bound << ";\n";
            }
          else if (symbol_table.getType(it->name) == eEndogenous)
            {
              output << "tmp1 = find((estim_params_.corrn(:,1)==" << symbol_table.getID(it->name)+1 << ")) & (estim_params_.corrn(:,2)==" << symbol_table.getID(it->name2)+1 << ";\n";
              output << "estim_params_.corrx(tmp1,4) = " << it->low_bound << ";\n";
              output << "estim_params_.corrx(tmp1,5) = " << it->up_bound << ";\n";
            }
        }
    }
}

ObservationTrendsStatement::ObservationTrendsStatement(const trend_elements_type &trend_elements_arg,
                                                       const SymbolTable &symbol_table_arg) :
  trend_elements(trend_elements_arg),
  symbol_table(symbol_table_arg)
{
}

void
ObservationTrendsStatement::writeOutput(ostream &output, const string &basename) const
{
  output << "options_.trend_coeff_ = {};" << endl;

  trend_elements_type::const_iterator it;

  for(it = trend_elements.begin(); it != trend_elements.end(); it++)
    {
      Type type = symbol_table.getType(it->first);
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

CalibVarStatement::CalibVarStatement(const calib_var_type &calib_var_arg,
                                     const calib_covar_type &calib_covar_arg,
                                     const calib_ac_type &calib_ac_arg,
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

  for(int i = 1; i < 4 ; i++)
    {
      output << "calib_var_index{" << i << "} = [];\n";
      output << "calib_targets{" << i << "} = [];\n";
      output << "calib_weights{" << i << "}=[];\n";
    }

  // Print calibration variances
  for(calib_var_type::const_iterator it = calib_var.begin();
      it != calib_var.end(); it++)
    {
      const string &name = it->first;
      const string &weight = it->second.first;
      const NodeID expression = it->second.second;

      int id = symbol_table.getID(name) + 1;
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
  for(calib_covar_type::const_iterator it = calib_covar.begin();
      it != calib_covar.end(); it++)
    {
      const string &name1 = it->first.first;
      const string &name2 = it->first.second;
      const string &weight = it->second.first;
      const NodeID expression = it->second.second;

      int id1 = symbol_table.getID(name1) + 1;
      int id2 = symbol_table.getID(name2) + 1;
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

  for(calib_ac_type::const_iterator it = calib_ac.begin();
      it != calib_ac.end(); it++)
    {
      const string &name = it->first.first;
      int iar = it->first.second + 3;
      const string &weight = it->second.first;
      const NodeID expression = it->second.second;

      int id = symbol_table.getID(name) + 1;

      if (iar > max_iar)
        {
          // Create new variables
          for(int i = max_iar + 1; i <= iar; i++)
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
  mod_file_struct.stoch_simul_or_similar_present = true;

  // Fill in option_order of mod_file_struct
  OptionsList::num_options_type::const_iterator it = options_list.num_options.find("order");
  if (it != options_list.num_options.end())
    mod_file_struct.order_option = max(mod_file_struct.order_option,atoi(it->second.c_str()));
}

void
OsrStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output);
  symbol_list.writeOutput("var_list_", output);
  output << "osr(var_list_,osr_params_,obj_var_,optim_weights_);\n";
}

OptimWeightsStatement::OptimWeightsStatement(const var_weights_type &var_weights_arg,
                                             const covar_weights_type &covar_weights_arg,
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

  for(var_weights_type::const_iterator it = var_weights.begin();
      it != var_weights.end(); it++)
    {
      const string &name = it->first;
      const NodeID value = it->second;
      int id = symbol_table.getID(name) + 1;
      output <<  "optim_weights_(" << id << "," << id << ") = ";
      value->writeOutput(output);
      output << ";" << endl;
      output << "obj_var_ = [obj_var_; " << id << "];\n";
    }

  for(covar_weights_type::const_iterator it = covar_weights.begin();
      it != covar_weights.end(); it++)
    {
      const string &name1 = it->first.first;
      const string &name2 = it->first.second;
      const NodeID value = it->second;
      int id1 = symbol_table.getID(name1) + 1;
      int id2 = symbol_table.getID(name2) + 1;
      output <<  "optim_weights_(" << id1 << "," << id2 << ") = ";
      value->writeOutput(output);
      output << ";" << endl;
      output << "obj_var_ = [obj_var_; " << id1 << " " << id2 << "];\n";
    }
}

DynaSaveStatement::DynaSaveStatement(const SymbolList &symbol_list_arg,
                                     const string &filename_arg, const string &ext_arg) :
  symbol_list(symbol_list_arg),
  filename(filename_arg),
  ext(ext_arg)
{
}

void
DynaSaveStatement::writeOutput(ostream &output, const string &basename) const
{
  symbol_list.writeOutput("var_list_", output);
  output << "dynasave('" << filename;
  if (ext.size() > 0)
    output << "," << ext;
  output << "',var_list_);\n";
}

DynaTypeStatement::DynaTypeStatement(const SymbolList &symbol_list_arg,
                                     const string &filename_arg, const string &ext_arg) :
  symbol_list(symbol_list_arg),
  filename(filename_arg),
  ext(ext_arg)
{
}

void
DynaTypeStatement::writeOutput(ostream &output, const string &basename) const
{
  symbol_list.writeOutput("var_list_", output);
  output << "dynatype(" << filename;
  if (ext.size() > 0)
    output << "," << ext;
  output << ",var_list_);\n";
}

ModelComparisonStatement::ModelComparisonStatement(const filename_list_type &filename_list_arg,
                                                   const OptionsList &options_list_arg) :
  filename_list(filename_list_arg),
  options_list(options_list_arg)
{
}

void
ModelComparisonStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output);

  output << "ModelNames_ = {};\n";
  output << "ModelPriors_ = {};\n";

  for(filename_list_type::const_iterator it = filename_list.begin();
      it != filename_list.end(); it++)
    {
      output << "ModelNames_ = { ModelNames_{:} '" << it->first << "};\n";
      output << "ModelPriors_ = { ModelPriors_{:} '" << it->second << "};\n";
    }
  output << "model_comparison(ModelNames_,ModelPriors_);\n";
}

PlannerObjectiveStatement::PlannerObjectiveStatement(ModelTree *model_tree_arg) :
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
  if (model_tree->equation_number() != 1)
    {
      cerr << "Error: planer_objective: should have only one equation!" << endl;
      exit(-1);
    }
}

void
PlannerObjectiveStatement::computingPass()
{
  model_tree->computeStaticHessian = true;
  model_tree->computingPass(eval_context_type());
}

void
PlannerObjectiveStatement::writeOutput(ostream &output, const string &basename) const
{
  model_tree->writeStaticFile(basename + "_objective");
}

BVARDensityStatement::BVARDensityStatement(int maxnlags_arg, const OptionsList &options_list_arg) :
  maxnlags(maxnlags_arg),
  options_list(options_list_arg)
{
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
BVARForecastStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output);
  output << "bvar_forecast(" << nlags << ");" << endl;
}
