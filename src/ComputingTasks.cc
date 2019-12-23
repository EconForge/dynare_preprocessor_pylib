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

#include <cassert>
#include <iostream>
#include <sstream>

using namespace std;

#include "ComputingTasks.hh"
#include "Statement.hh"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/tokenizer.hpp>
#pragma GCC diagnostic pop

#include <utility>

SteadyStatement::SteadyStatement(OptionsList options_list_arg) :
  options_list{move(options_list_arg)}
{
}

void
SteadyStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.steady_present = true;
}

void
SteadyStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  options_list.writeOutput(output);
  output << "steady;" << endl;
}

void
SteadyStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "steady")";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

CheckStatement::CheckStatement(OptionsList options_list_arg) :
  options_list{move(options_list_arg)}
{
}

void
CheckStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  options_list.writeOutput(output);
  output << "oo_.dr.eigval = check(M_,options_,oo_);" << endl;
}

void
CheckStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.check_present = true;
}

void
CheckStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "check")";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

ModelInfoStatement::ModelInfoStatement(OptionsList options_list_arg) :
  options_list{move(options_list_arg)}
{
}

void
ModelInfoStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  //mod_file_struct.model_info_present = true;
}

void
ModelInfoStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  options_list.writeOutput(output);
  output << "model_info();" << endl;
}

void
ModelInfoStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "model_info")";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

SimulStatement::SimulStatement(OptionsList options_list_arg) :
  options_list{move(options_list_arg)}
{
}

void
SimulStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.perfect_foresight_solver_present = true;
}

void
SimulStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  // Translate the “datafile” option into “initval_file” (see dynare#1663)
  auto options_list_new = options_list; // Need a copy, because of const
  if (auto it = options_list_new.string_options.find("datafile");
      it != options_list_new.string_options.end())
    {
      output << "options_.initval_file = true;" << endl
             << "initvalf('" << it->second << "');" << endl;
      options_list_new.string_options.erase(it);
    }
  options_list_new.writeOutput(output);
  output << "perfect_foresight_setup;" << endl
         << "perfect_foresight_solver;" << endl;
}

void
SimulStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "simul")";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

PerfectForesightSetupStatement::PerfectForesightSetupStatement(OptionsList options_list_arg) :
  options_list{move(options_list_arg)}
{
}

void
PerfectForesightSetupStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  auto options_list_new = options_list; // Need a copy, because of const
  if (auto it = options_list_new.string_options.find("datafile");
      it != options_list_new.string_options.end())
    {
      output << "options_.initval_file = true;" << endl
             << "initvalf('" << it->second << "');" << endl;
      options_list_new.string_options.erase(it);
    }
  options_list_new.writeOutput(output);
  output << "perfect_foresight_setup;" << endl;
}

void
PerfectForesightSetupStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "perfect_foresight_setup")";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

PerfectForesightSolverStatement::PerfectForesightSolverStatement(OptionsList options_list_arg) :
  options_list(move(options_list_arg))
{
}

void
PerfectForesightSolverStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.perfect_foresight_solver_present = true;
  // Fill in option_occbin of mod_file_struct
  if (options_list.num_options.find("occbin") != options_list.num_options.end())
    mod_file_struct.occbin_option = true;
}

void
PerfectForesightSolverStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  options_list.writeOutput(output);
  output << "perfect_foresight_solver;" << endl;
}

void
PerfectForesightSolverStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "perfect_foresight_solver")";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

PriorPosteriorFunctionStatement::PriorPosteriorFunctionStatement(const bool prior_func_arg,
                                                                 OptionsList options_list_arg) :
  prior_func{prior_func_arg},
  options_list{move(options_list_arg)}
{
}

void
PriorPosteriorFunctionStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  if (auto it2 = options_list.string_options.find("function");
      it2 == options_list.string_options.end() || it2->second.empty())
    {
      cerr << "ERROR: both the prior_function and posterior_function commands require the 'function' argument"
           << endl;
      exit(EXIT_FAILURE);
    }
}

void
PriorPosteriorFunctionStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  options_list.writeOutput(output);
  string type = prior_func ? "prior" : "posterior";

  output << "oo_ = execute_prior_posterior_function("
         << "'" << options_list.string_options.find("function")->second << "', "
         << "M_, options_, oo_, estim_params_, bayestopt_, dataset_, dataset_info, "
         << "'" << type << "');" << endl;
}

void
PriorPosteriorFunctionStatement::writeJsonOutput(ostream &output) const
{
  string type = prior_func ? "prior" : "posterior";
  output << R"({"statementName": "prior_posterior_function", "type": ")" << type << R"(")";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

PacModelStatement::PacModelStatement(string name_arg,
                                     string aux_model_name_arg,
                                     string discount_arg,
                                     expr_t growth_arg,
                                     double steady_state_growth_rate_number_arg,
                                     int steady_state_growth_rate_symb_id_arg,
                                     const SymbolTable &symbol_table_arg) :
  name{move(name_arg)},
  aux_model_name{move(aux_model_name_arg)},
  discount{move(discount_arg)},
  growth{growth_arg},
  original_growth{growth_arg},
  steady_state_growth_rate_number{steady_state_growth_rate_number_arg},
  steady_state_growth_rate_symb_id{steady_state_growth_rate_symb_id_arg},
  symbol_table{symbol_table_arg}
{
}

void
PacModelStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  if (growth)
    growth->collectVariables(SymbolType::exogenous, mod_file_struct.pac_params);
}

void
PacModelStatement::overwriteGrowth(expr_t new_growth)
{
  if (!new_growth || !growth)
    return;

  growth = new_growth;

  try
    {
      growth_info = growth->matchLinearCombinationOfVariables(false);
    }
  catch (ExprNode::MatchFailureException &e)
    {
      auto gv = dynamic_cast<const VariableNode *>(growth);
      if (gv)
        growth_info.emplace_back(gv->symb_id, gv->lag, -1, 1);
      else
        {
          cerr << "Pac growth must be a linear combination of variables" << endl;
          exit(EXIT_FAILURE);
        }
    }
}

void
PacModelStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "M_.pac." << name << ".auxiliary_model_name = '" << aux_model_name << "';" << endl
         << "M_.pac." << name << ".discount_index = " << symbol_table.getTypeSpecificID(discount) + 1 << ";" << endl;
  if (steady_state_growth_rate_symb_id < 0 && steady_state_growth_rate_number > 0)
    output << "M_.pac." << name << ".steady_state_growth_rate = "
           << steady_state_growth_rate_number << ";" << endl;
  else if (steady_state_growth_rate_symb_id >= 0)
    output << "M_.pac." << name << ".steady_state_growth_rate = "
           << symbol_table.getTypeSpecificID(steady_state_growth_rate_symb_id) + 1 << ";" << endl;

  if (growth)
    {
      output << "M_.pac." << name << ".growth_str = '";
      original_growth->writeJsonOutput(output, {}, {}, true);
      output << "';" << endl;
      int i = 0;
      for (auto [growth_symb_id, growth_lag, param_id, constant] : growth_info)
        {
          string structname = "M_.pac." + name + ".growth_linear_comb(" + to_string(++i) + ").";
          if (growth_symb_id >= 0)
            {
              string var_field = "endo_id";
              if (symbol_table.getType(growth_symb_id) == SymbolType::exogenous)
                {
                  var_field = "exo_id";
                  output << structname << "endo_id = 0;" << endl;
                }
              else
                output << structname << "exo_id = 0;" << endl;
              try
                {
                  // case when this is not the highest lag of the growth variable
                  int aux_symb_id = symbol_table.searchAuxiliaryVars(growth_symb_id, growth_lag);
                  output << structname << var_field << " = " << symbol_table.getTypeSpecificID(aux_symb_id) + 1 << ";" << endl
                         << structname << "lag = 0;" << endl;
                }
              catch (...)
                {
                  try
                    {
                      // case when this is the highest lag of the growth variable
                      int tmp_growth_lag = growth_lag + 1;
                      int aux_symb_id = symbol_table.searchAuxiliaryVars(growth_symb_id, tmp_growth_lag);
                      output << structname << var_field << " = " << symbol_table.getTypeSpecificID(aux_symb_id) + 1 << ";" << endl
                             << structname << "lag = -1;" << endl;
                    }
                  catch (...)
                    {
                      // case when there is no aux var for the variable
                      output << structname << var_field << " = "<< symbol_table.getTypeSpecificID(growth_symb_id) + 1 << ";" << endl
                             << structname << "lag = " << growth_lag << ";" << endl;
                    }
                }
            }
          else
            output << structname << "endo_id = 0;" << endl
                   << structname << "exo_id = 0;" << endl
                   << structname << "lag = 0;" << endl;
          output << structname << "param_id = "
                 << (param_id == -1 ? 0 : symbol_table.getTypeSpecificID(param_id) + 1) << ";" << endl
                 << structname << "constant = " << constant << ";" << endl;
        }
    }
}

void
PacModelStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "pac_model",)"
         << R"("model_name": ")" << name << R"(",)"
         << R"("auxiliary_model_name": ")" << aux_model_name << R"(",)"
         << R"("discount_index": )" << symbol_table.getTypeSpecificID(discount) + 1;
  if (growth)
    {
      output << R"(,"growth_str": ")";
      original_growth->writeJsonOutput(output, {}, {}, true);
      output << R"(")";
    }
  output << "}" << endl;
}

VarEstimationStatement::VarEstimationStatement(OptionsList options_list_arg) :
  options_list{move(options_list_arg)}
{
}

void
VarEstimationStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  if (auto it = options_list.string_options.find("var_estimation.model_name");
      it == options_list.string_options.end())
    {
      cerr << "ERROR: You must provide the model name to the var_estimation statement." << endl;
      exit(EXIT_FAILURE);
    }
}

void
VarEstimationStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  options_list.writeOutput(output);
  output << "oo_ = var_estimation(M_, options_, oo_);" << endl;
}

VarRestrictionsStatement::VarRestrictionsStatement(string var_model_name_arg,
                                                   map<string, vector<string>> var_map_arg,
                                                   map<int, map<int, SymbolList>> exclusion_restrictions_arg,
                                                   equation_restrictions_t equation_restrictions_arg,
                                                   crossequation_restrictions_t crossequation_restrictions_arg,
                                                   map<pair<int, int>, double> covariance_number_restriction_arg,
                                                   map<pair<int, int>, pair<int, int>> covariance_pair_restriction_arg,
                                                   const SymbolTable &symbol_table_arg) :
  var_model_name{move(var_model_name_arg)},
  var_map{move(var_map_arg)},
  exclusion_restrictions{move(exclusion_restrictions_arg)},
  equation_restrictions{move(equation_restrictions_arg)},
  crossequation_restrictions{move(crossequation_restrictions_arg)},
  covariance_number_restriction{move(covariance_number_restriction_arg)},
  covariance_pair_restriction{move(covariance_pair_restriction_arg)},
  symbol_table{symbol_table_arg}
{
}

int
VarRestrictionsStatement::findIdxInVector(const vector<string> &vecvars, const string &var) const
{
  int idx = 0;
  bool setflag = false;
  for (auto itvs = vecvars.begin(); itvs != vecvars.end(); ++itvs, idx++)
    if (*itvs == var)
      {
        setflag = true;
        break;
      }

  if (!setflag)
    {
      cerr << "ERROR: you are imposing an exclusion restriction on an equation or variable "
           << var << " that is not contained in VAR " << var_model_name;
      exit(EXIT_FAILURE);
    }
  return idx;
}

void
VarRestrictionsStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  auto itvs = var_map.find(var_model_name);
  if (itvs == var_map.end())
    {
      cerr << "ERROR: you are imposing restrictions on a VAR named " << var_model_name
           << " but this VAR has not been declared via thevar_model statement." << endl;
      exit(EXIT_FAILURE);
    }
  vector<string> vars = itvs->second;

  string Mstr("M_.var." + var_model_name + ".restrictions.");
  int nrestrictions = 0;

  // Exclusion Restrictions
  int idx = 1;
  for (auto it = exclusion_restrictions.begin();
       it != exclusion_restrictions.end(); ++it, idx++)
    {
      output << Mstr << "exclusion_restrictions{" << idx<< "}.lag = "
             << it->first << ";" << endl
             << Mstr << "exclusion_restrictions{" << idx << "}.restrictions = [";
      for (auto it1 = it->second.begin(); it1 != it->second.end(); ++it1)
        {
          if (it1 != it->second.begin())
            output << " ";

          output << "struct('eq', " << findIdxInVector(vars, symbol_table.getName(it1->first)) + 1
                 << ", 'vars', [";
          vector<string> excvars = it1->second.getSymbols();
          for (auto &excvar : excvars)
            output << findIdxInVector(vars, excvar) + 1 << " ";
          output << "])";
          nrestrictions += it1->second.getSize();
        }
      output << "];" << endl;
    }

  // Equation Restrictions
  idx = 1;
  for (auto it = equation_restrictions.begin();
       it != equation_restrictions.end(); ++it, idx++, nrestrictions++)
    {
      output << Mstr << "equation_restriction{" << idx << "}.eq = '"
             << symbol_table.getName(it->first) << "';" << endl
             << Mstr << "equation_restriction{" << idx << "}.val = "
             << it->second.second << ";" << endl;

      var_restriction_eq_crosseq_t ls = it->second.first.first;
      output << Mstr << "equation_restriction{" << idx << "}.ls = '"
             << symbol_table.getName(ls.first.first) << "';" << endl
             << Mstr << "equation_restriction{" << idx << "}.lslag = "
             << ls.first.second.second << ";" << endl
             << Mstr << "equation_restriction{" << idx << "}.lscoeff = ";
      ls.second->writeOutput(output);
      output << ";" << endl;

      var_restriction_eq_crosseq_t rs = it->second.first.second;
      if (rs.first.first >= 0)
        {
          output << Mstr << "equation_restriction{" << idx << "}.rs = '"
                 << symbol_table.getName(rs.first.first) << "';" << endl
                 << Mstr << "equation_restriction{" << idx << "}.rslag = "
                 << rs.first.second.second << ";" << endl
                 << Mstr << "equation_restriction{" << idx << "}.rscoeff = ";
          rs.second->writeOutput(output);
          output << ";" << endl;
        }
    }

  // Cross Equation Restrictions
  idx = 1;
  for (auto it = crossequation_restrictions.begin();
       it != crossequation_restrictions.end(); ++it, idx++, nrestrictions++)
    {
      output << Mstr << "crossequation_restriction{" << idx << "}.val = "
             << it->second << ";" << endl;

      var_restriction_eq_crosseq_t ls = it->first.first;
      output << Mstr << "crossequation_restriction{" << idx << "}.lseq = "
             << findIdxInVector(vars, symbol_table.getName(ls.first.first)) + 1 << ";" << endl
             << Mstr << "crossequation_restriction{" << idx << "}.lsvar = "
             << findIdxInVector(vars, symbol_table.getName(ls.first.second.first)) + 1 << ";" << endl
             << Mstr << "crossequation_restriction{" << idx << "}.lslag = "
             << ls.first.second.second << ";" << endl
             << Mstr << "crossequation_restriction{" << idx << "}.lscoeff = ";
      ls.second->writeOutput(output);
      output << ";" << endl;

      var_restriction_eq_crosseq_t rs = it->first.second;
      if (rs.first.first >= 0)
        {
          output << Mstr << "crossequation_restriction{" << idx << "}.rseq = "
                 << findIdxInVector(vars, symbol_table.getName(rs.first.first)) + 1 << ";" << endl
                 << Mstr << "crossequation_restriction{" << idx << "}.rsvar = "
                 << findIdxInVector(vars, symbol_table.getName(rs.first.second.first)) + 1 << ";" << endl
                 << Mstr << "crossequation_restriction{" << idx << "}.rslag = "
                 << rs.first.second.second << ";" << endl
                 << Mstr << "crossequation_restriction{" << idx << "}.rscoeff = ";
          rs.second->writeOutput(output);
          output << ";" << endl;
        }
    }

  // Covariance Const Restrictions
  idx = 1;
  for (auto it = covariance_number_restriction.begin();
       it != covariance_number_restriction.end(); ++it, idx++)
    output << Mstr << "covariance_const_restriction{" << idx << "}.var1 = '"
           << symbol_table.getName(it->first.first) << "';" << endl
           << Mstr << "covariance_const_restriction{" << idx << "}.var2 = '"
           << symbol_table.getName(it->first.second) << "';" << endl
           << Mstr << "covariance_const_restriction{" << idx << "}.val = "
           << it->second << ";" << endl;

  // Covariance Pair Restrictions
  idx = 1;
  for (auto it = covariance_pair_restriction.begin();
       it != covariance_pair_restriction.end(); ++it, idx++)
    output << Mstr << "covariance_pair_restriction{" << idx << "}.var11 = '"
           << symbol_table.getName(it->first.first) << "';" << endl
           << Mstr << "covariance_pair_restriction{" << idx << "}.var12 = '"
           << symbol_table.getName(it->first.second) << "';" << endl
           << Mstr << "covariance_pair_restriction{" << idx << "}.var21 = '"
           << symbol_table.getName(it->second.first) << "';" << endl
           << Mstr << "covariance_pair_restriction{" << idx << "}.var22 = '"
           << symbol_table.getName(it->second.second) << "';" << endl;

  output << Mstr << "N = " << nrestrictions << ";" << endl;
}

StochSimulStatement::StochSimulStatement(SymbolList symbol_list_arg,
                                         OptionsList options_list_arg) :
  symbol_list{move(symbol_list_arg)},
  options_list{move(options_list_arg)}
{
}

void
StochSimulStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.stoch_simul_present = true;

  // Fill in option_order of mod_file_struct
  if (auto it = options_list.num_options.find("order");
      it != options_list.num_options.end())
    mod_file_struct.order_option = max(mod_file_struct.order_option, stoi(it->second));

  // Fill in mod_file_struct.partial_information
  if (auto it = options_list.num_options.find("partial_information");
      it != options_list.num_options.end() && it->second == "true")
    mod_file_struct.partial_information = true;

  // Option k_order_solver (implicit when order >= 3)
  if (auto it = options_list.num_options.find("k_order_solver");
      (it != options_list.num_options.end() && it->second == "true")
      || mod_file_struct.order_option >= 3)
    mod_file_struct.k_order_solver = true;

  if (auto it = options_list.num_options.find("hp_filter"),
      it1 = options_list.num_options.find("bandpass.indicator"),
      it2 = options_list.num_options.find("one_sided_hp_filter");
      (it != options_list.num_options.end() && it1 != options_list.num_options.end())
      || (it != options_list.num_options.end() && it2 != options_list.num_options.end())
      || (it1 != options_list.num_options.end() && it2 != options_list.num_options.end()))
    {
      cerr << "ERROR: stoch_simul: can only use one of hp, one-sided hp, and bandpass filters"
           << endl;
      exit(EXIT_FAILURE);
    }

  symbol_list.removeDuplicates("stoch_simul", warnings);

  try
    {
      const vector<SymbolType> valid_symbol_list_types { SymbolType::endogenous };
      symbol_list.checkPass(warnings, valid_symbol_list_types);
    }
  catch (SymbolList::SymbolListException &e)
    {
      cerr << "ERROR: stoch_simul: " << e.message << endl;
      exit(EXIT_FAILURE);
    }
}

void
StochSimulStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  // Ensure that order 3 implies k_order (#844)
  if (auto it = options_list.num_options.find("order"),
      it1 = options_list.num_options.find("k_order_solver");
      (it1 != options_list.num_options.end() && it1->second == "true")
      || (it != options_list.num_options.end() && stoi(it->second) >= 3))
    output << "options_.k_order_solver = true;" << endl;

  options_list.writeOutput(output);
  symbol_list.writeOutput("var_list_", output);
  output << "[info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_);" << endl;
}

void
StochSimulStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "stoch_simul")";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  if (!symbol_list.empty())
    {
      output << ", ";
      symbol_list.writeJsonOutput(output);
    }
  output << "}";
}

ForecastStatement::ForecastStatement(SymbolList symbol_list_arg,
                                     OptionsList options_list_arg) :
  symbol_list{move(symbol_list_arg)},
  options_list{move(options_list_arg)}
{
}

void
ForecastStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  try
    {
      const vector<SymbolType> valid_symbol_list_types { SymbolType::endogenous };
      symbol_list.checkPass(warnings, valid_symbol_list_types);
    }
  catch (SymbolList::SymbolListException &e)
    {
      cerr << "ERROR: forecast: " << e.message << endl;
      exit(EXIT_FAILURE);
    }
}

void
ForecastStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  options_list.writeOutput(output);
  symbol_list.writeOutput("var_list_", output);
  output << "[oo_.forecast,info] = dyn_forecast(var_list_,M_,options_,oo_,'simul');" << endl;
}

void
ForecastStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "forecast")";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  if (!symbol_list.empty())
    {
      output << ", ";
      symbol_list.writeJsonOutput(output);
    }
  output << "}";
}

DetCondForecast::DetCondForecast(SymbolList symbol_list_arg,
                                 OptionsList options_list_arg,
                                 const bool linear_decomposition_arg) :
  options_list{move(options_list_arg)},
  symbol_list{move(symbol_list_arg)},
  linear_decomposition{linear_decomposition_arg}
{

}

void
DetCondForecast::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  options_list.writeOutput(output);
  if (linear_decomposition)
    {
      output << "first_order_solution_to_compute = 1;" << endl
             << "if isfield(oo_, 'dr')" << endl
             << "    if isfield(oo_.dr, 'ghx') && isfield(oo_.dr, 'ghu') && isfield(oo_.dr, 'state_var') && isfield(oo_.dr, 'order_var')" << endl
             << "        first_order_solution_to_compute = 0;" << endl
             << "    end" << endl
             << "end" << endl
             << "if first_order_solution_to_compute" << endl
             << "  fprintf('%s','Computing the first order solution ...');" << endl
             << "  options_.nograph = true;" << endl
             << "  options_.order = 1;" << endl
             << "  options_.noprint = true;" << endl
             << "  options_.nocorr = true;" << endl
             << "  options_.nomoments = true;" << endl
             << "  options_.nodecomposition = true;" << endl
             << "  options_.nofunctions = true;" << endl
             << "  options_.irf = 0;" << endl
             << "  tmp_periods = options_.periods;" << endl
             << "  options_.periods = 0;" << endl
             << "  var_list_ = char();" << endl
             << "  info = stoch_simul(var_list_);" << endl
             << R"(  fprintf('%s\n','done');)" << endl
             << "  options_.periods = tmp_periods;" << endl
             << "end" << endl;
    }
  vector<string> symbols = symbol_list.get_symbols();
  if (symbols.size() > 0)
    output << symbols[1] << " = det_cond_forecast(";
  for (unsigned int i = 0; i < symbols.size() - 1; i++)
    output << symbols[i] << ", ";
  if (symbols.size() > 0)
    output << symbols[symbols.size() - 1];
  output << ");" << endl;
}

RamseyModelStatement::RamseyModelStatement(OptionsList options_list_arg) :
  options_list{move(options_list_arg)}
{
}

void
RamseyModelStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.ramsey_model_present = true;

  /* Fill in option_order of mod_file_struct
     Since ramsey model needs one further order of derivation (for example, for 1st order
     approximation, it needs 2nd derivatives), we add 1 to the order declared by user */
  if (auto it = options_list.num_options.find("order");
      it != options_list.num_options.end())
    {
      int order = stoi(it->second);
      if (order > 2)
        {
          cerr << "ERROR: ramsey_model: order > 2 is not  implemented" << endl;
          exit(EXIT_FAILURE);
        }
      mod_file_struct.order_option = max(mod_file_struct.order_option, order + 1);
    }

  // Fill in mod_file_struct.partial_information
  if (auto it = options_list.num_options.find("partial_information");
      it != options_list.num_options.end() && it->second == "true")
    mod_file_struct.partial_information = true;

  // Option k_order_solver (implicit when order >= 3)
  if (auto it = options_list.num_options.find("k_order_solver");
      (it != options_list.num_options.end() && it->second == "true")
      || mod_file_struct.order_option >= 3)
    mod_file_struct.k_order_solver = true;

  // Fill list of instruments
  if (auto it = options_list.symbol_list_options.find("instruments");
      it != options_list.symbol_list_options.end())
    mod_file_struct.instruments = it->second;
}

void
RamseyModelStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  // options_.ramsey_policy indicates that a Ramsey model is present in the *.mod file
  // this affects the computation of the steady state that uses a special algorithm
  // It should probably rather be a M_ field, but we leave it in options_ for historical reason

  // Ensure that order 3 implies k_order (#844)
  if (auto it = options_list.num_options.find("order"),
      it1 = options_list.num_options.find("k_order_solver");
      (it1 != options_list.num_options.end() && it1->second == "true")
      || (it != options_list.num_options.end() && stoi(it->second) >= 3))
    output << "options_.k_order_solver = true;" << endl;

  output << "options_.ramsey_policy = true;" << endl;
  options_list.writeOutput(output);
}

void
RamseyModelStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "ramsey_model")";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

RamseyConstraintsStatement::RamseyConstraintsStatement(const SymbolTable &symbol_table_arg, constraints_t constraints_arg) :
  symbol_table{symbol_table_arg},
  constraints{move(constraints_arg)}
{
}

void
RamseyConstraintsStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  if (!mod_file_struct.ramsey_model_present || !mod_file_struct.ramsey_policy_present)
    cerr << "ramsey_constraints: can only be used with ramsey_model or ramsey_policy" << endl;
}

void
RamseyConstraintsStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "M_.ramsey_model_constraints = {" << endl;
  for (auto it = constraints.begin(); it != constraints.end(); ++it)
    {
      if (it != constraints.begin())
        output << ", ";
      output << "{" << it->endo + 1 << ", '";
      switch (it->code)
        {
        case BinaryOpcode::less:
          output << '<';
          break;
        case BinaryOpcode::greater:
          output << '>';
          break;
        case BinaryOpcode::lessEqual:
          output << "<=";
          break;
        case BinaryOpcode::greaterEqual:
          output << ">=";
          break;
        default:
          cerr << "Ramsey constraints: this shouldn't happen." << endl;
          exit(EXIT_FAILURE);
        }
      output << "', '";
      it->expression->writeOutput(output);
      output << "'}" << endl;
    }
  output << "};" << endl;
}

void
RamseyConstraintsStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "ramsey_constraints")"
         << R"(, "ramsey_model_constraints": [)" << endl;
  for (auto it = constraints.begin(); it != constraints.end(); ++it)
    {
      if (it != constraints.begin())
        output << ", ";
      output << R"({"constraint": ")" << symbol_table.getName(it->endo) << " ";
      switch (it->code)
        {
        case BinaryOpcode::less:
          output << '<';
          break;
        case BinaryOpcode::greater:
          output << '>';
          break;
        case BinaryOpcode::lessEqual:
          output << "<=";
          break;
        case BinaryOpcode::greaterEqual:
          output << ">=";
          break;
        default:
          cerr << "Ramsey constraints: this shouldn't happen." << endl;
          exit(EXIT_FAILURE);
        }
      output << " ";
      it->expression->writeJsonOutput(output, {}, {});
      output << R"("})" << endl;
    }
  output << "]" << endl;
  output << "}";
}

RamseyPolicyStatement::RamseyPolicyStatement(const SymbolTable &symbol_table_arg,
                                             SymbolList symbol_list_arg,
                                             OptionsList options_list_arg) :
  symbol_table{symbol_table_arg},
  symbol_list{move(symbol_list_arg)},
  options_list{move(options_list_arg)}
{
}

void
RamseyPolicyStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  // ramsey_model_present indicates that the model is augmented with the FOC of the planner problem
  mod_file_struct.ramsey_model_present = true;
  // ramsey_policy_present indicates that ramsey_policy instruction for computation of first order approximation
  // of  a stochastic Ramsey problem if present in the *.mod file
  mod_file_struct.ramsey_policy_present = true;

  /* Fill in option_order of mod_file_struct
     Since ramsey policy needs one further order of derivation (for example, for 1st order
     approximation, it needs 2nd derivatives), we add 1 to the order declared by user */
  if (auto it = options_list.num_options.find("order");
      it != options_list.num_options.end())
    {
      int order = stoi(it->second);
      if (order > 2)
        {
          cerr << "ERROR: ramsey_policy: order > 2 is not  implemented" << endl;
          exit(EXIT_FAILURE);
        }
      mod_file_struct.order_option = max(mod_file_struct.order_option, order + 1);
    }

  // Fill in mod_file_struct.partial_information
  if (auto it = options_list.num_options.find("partial_information");
      it != options_list.num_options.end() && it->second == "true")
    mod_file_struct.partial_information = true;

  // Option k_order_solver (implicit when order >= 3)
  if (auto it = options_list.num_options.find("k_order_solver");
      (it != options_list.num_options.end() && it->second == "true")
      || mod_file_struct.order_option >= 3)
    mod_file_struct.k_order_solver = true;

  // Fill list of instruments
  if (auto it = options_list.symbol_list_options.find("instruments");
      it != options_list.symbol_list_options.end())
    mod_file_struct.instruments = it->second;

  try
    {
      const vector<SymbolType> valid_symbol_list_types { SymbolType::endogenous };
      symbol_list.checkPass(warnings, valid_symbol_list_types);
    }
  catch (SymbolList::SymbolListException &e)
    {
      cerr << "ERROR: ramsey_policy: " << e.message << endl;
      exit(EXIT_FAILURE);
    }
}

void
RamseyPolicyStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  // Ensure that order 3 implies k_order (#844)
  if (auto it = options_list.num_options.find("order"),
      it1 = options_list.num_options.find("k_order_solver");
      (it1 != options_list.num_options.end() && it1->second == "true")
      || (it != options_list.num_options.end() && stoi(it->second) >= 3))
    output << "options_.k_order_solver = true;" << endl;

  options_list.writeOutput(output);
  symbol_list.writeOutput("var_list_", output);
  output << "ramsey_policy(var_list_);" << endl;
}

void
RamseyPolicyStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "ramsey_policy")";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  if (!symbol_list.empty())
    {
      output << ", ";
      symbol_list.writeJsonOutput(output);
    }
  output << "}";
}

void
EvaluatePlannerObjective::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "oo_.planner_objective_value = evaluate_planner_objective(M_, options_, oo_);" << endl;
}

void
EvaluatePlannerObjective::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "evaluate_planner_objective"})";
}

DiscretionaryPolicyStatement::DiscretionaryPolicyStatement(SymbolList symbol_list_arg,
                                                           OptionsList options_list_arg) :
  symbol_list{move(symbol_list_arg)},
  options_list{move(options_list_arg)}
{
}

void
DiscretionaryPolicyStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.discretionary_policy_present = true;

  if (options_list.symbol_list_options.find("instruments") == options_list.symbol_list_options.end())
    {
      cerr << "ERROR: discretionary_policy: the instruments option is required." << endl;
      exit(EXIT_FAILURE);
    }

  /* Fill in option_order of mod_file_struct
     Since discretionary policy needs one further order of derivation (for example, for 1st order
     approximation, it needs 2nd derivatives), we add 1 to the order declared by user */
  if (auto it = options_list.num_options.find("order");
      it != options_list.num_options.end())
    {
      int order = stoi(it->second);
      if (order > 1)
        {
          cerr << "ERROR: discretionary_policy: order > 1 is not yet implemented" << endl;
          exit(EXIT_FAILURE);
        }
      mod_file_struct.order_option = max(mod_file_struct.order_option, order + 1);
    }

  // Fill in mod_file_struct.partial_information
  if (auto it = options_list.num_options.find("partial_information");
      it != options_list.num_options.end() && it->second == "true")
    mod_file_struct.partial_information = true;

  // Option k_order_solver (implicit when order >= 3)
  if (auto it = options_list.num_options.find("k_order_solver");
      (it != options_list.num_options.end() && it->second == "true")
      || mod_file_struct.order_option >= 3)
    mod_file_struct.k_order_solver = true;

  // Fill list of instruments
  if (auto it = options_list.symbol_list_options.find("instruments");
      it != options_list.symbol_list_options.end())
    mod_file_struct.instruments = it->second;

  try
    {
      const vector<SymbolType> valid_symbol_list_types { SymbolType::endogenous };
      symbol_list.checkPass(warnings, valid_symbol_list_types);
    }
  catch (SymbolList::SymbolListException &e)
    {
      cerr << "ERROR: discretionary_policy: " << e.message << endl;
      exit(EXIT_FAILURE);
    }
}

void
DiscretionaryPolicyStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  // Ensure that order 3 implies k_order (#844)
  if (auto it = options_list.num_options.find("order"),
      it1 = options_list.num_options.find("k_order_solver");
      (it1 != options_list.num_options.end() && it1->second == "true")
      || (it != options_list.num_options.end() && stoi(it->second) >= 3))
    output << "options_.k_order_solver = true;" << endl;

  options_list.writeOutput(output);
  symbol_list.writeOutput("var_list_", output);
  output << "[info, oo_, options_] = discretionary_policy(M_, options_, oo_, var_list_);" << endl;
}

void
DiscretionaryPolicyStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "discretionary_policy")";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  if (!symbol_list.empty())
    {
      output << ", ";
      symbol_list.writeJsonOutput(output);
    }
  output << "}";
}

EstimationStatement::EstimationStatement(SymbolList symbol_list_arg,
                                         OptionsList options_list_arg) :
  symbol_list{move(symbol_list_arg)},
  options_list{move(options_list_arg)}
{
}

void
EstimationStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.estimation_present = true;

  // Fill in option_order of mod_file_struct
  if (auto it = options_list.num_options.find("order");
      it != options_list.num_options.end())
    {
      int order = stoi(it->second);

      if (order > 2)
        mod_file_struct.k_order_solver = true;

      mod_file_struct.order_option = max(mod_file_struct.order_option, order);
    }

  // Fill in mod_file_struct.partial_information
  if (auto it = options_list.num_options.find("partial_information");
      it != options_list.num_options.end() && it->second == "true")
    mod_file_struct.partial_information = true;

  // Fill in mod_file_struct.estimation_analytic_derivation
  if (auto it = options_list.num_options.find("analytic_derivation");
      it != options_list.num_options.end() && it->second == "1")
    mod_file_struct.estimation_analytic_derivation = true;

  if (auto it = options_list.num_options.find("dsge_var");
      it != options_list.num_options.end())
    // Fill in mod_file_struct.dsge_var_calibrated
    mod_file_struct.dsge_var_calibrated = it->second;

  // Fill in mod_file_struct.dsge_var_estimated
  if (options_list.string_options.find("dsge_var") != options_list.string_options.end())
    mod_file_struct.dsge_var_estimated = true;

  // Fill in mod_file_struct.bayesian_irf_present
  if (auto it = options_list.num_options.find("bayesian_irf");
      it != options_list.num_options.end() && it->second == "true")
    mod_file_struct.bayesian_irf_present = true;

  if (options_list.num_options.find("dsge_varlag") != options_list.num_options.end())
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

  if (options_list.string_options.find("datafile") == options_list.string_options.end()
      && !mod_file_struct.estimation_data_statement_present)
    {
      cerr << "ERROR: The estimation statement requires a data file to be supplied via the datafile option." << endl;
      exit(EXIT_FAILURE);
    }

  if (options_list.string_options.find("mode_file") != options_list.string_options.end()
      && mod_file_struct.estim_params_use_calib)
    {
      cerr << "ERROR: The mode_file option of the estimation statement is incompatible with the use_calibration option of the estimated_params_init block." << endl;
      exit(EXIT_FAILURE);
    }

  try
    {
      const vector<SymbolType> valid_symbol_list_types { SymbolType::endogenous };
      symbol_list.checkPass(warnings, valid_symbol_list_types);
    }
  catch (SymbolList::SymbolListException &e)
    {
      cerr << "ERROR: estimation: " << e.message << endl;
      exit(EXIT_FAILURE);
    }
}

void
EstimationStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  options_list.writeOutput(output);

  // Special treatment for order option and particle filter
  if (auto it = options_list.num_options.find("order");
      it == options_list.num_options.end())
    output << "options_.order = 1;" << endl;
  else if (stoi(it->second) >= 2)
    {
      output << "options_.particle.status = true;" << endl;
      if (stoi(it->second) > 2)
        output << "options_.k_order_solver = true;" << endl;
    }

  // Do not check for the steady state in diffuse filter mode (#400)
  if (auto it = options_list.num_options.find("diffuse_filter");
      it != options_list.num_options.end() && it->second == "true")
    output << "options_.steadystate.nocheck = true;" << endl;

  symbol_list.writeOutput("var_list_", output);
  output << "oo_recursive_=dynare_estimation(var_list_);" << endl;
}

void
EstimationStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "estimation")";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  if (!symbol_list.empty())
    {
      output << ", ";
      symbol_list.writeJsonOutput(output);
    }
  output << "}";
}

DynareSensitivityStatement::DynareSensitivityStatement(OptionsList options_list_arg) :
  options_list{move(options_list_arg)}
{
}

void
DynareSensitivityStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  if (auto it = options_list.num_options.find("identification");
      it != options_list.num_options.end() && it->second == "1")
    mod_file_struct.identification_present = true;
  mod_file_struct.sensitivity_present = true;
}

void
DynareSensitivityStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  options_list.writeOutput(output, "options_gsa");

  /* Ensure that nograph, nodisplay and graph_format are also set in top-level
     options_.
     \todo factorize this code between identification and dynare_sensitivity,
     and provide a generic mechanism for this situation (maybe using regexps) */
  if (auto it = options_list.num_options.find("nodisplay");
      it != options_list.num_options.end())
    output << "options_.nodisplay = " << it->second << ";" << endl;
  if (auto it = options_list.num_options.find("nograph");
      it != options_list.num_options.end())
    output << "options_.nograph = " << it->second << ";" << endl;
  if (auto it = options_list.string_options.find("graph_format");
      it != options_list.string_options.end())
    output << "options_.graph_format = '" << it->second << "';" << endl;

  output << "dynare_sensitivity(options_gsa);" << endl;
}

void
DynareSensitivityStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "dynare_sensitivity")";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

RplotStatement::RplotStatement(SymbolList symbol_list_arg) :
  symbol_list{move(symbol_list_arg)}
{
}

void
RplotStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  try
    {
      const vector<SymbolType> valid_symbol_list_types { SymbolType::endogenous,
                                                         SymbolType::exogenous};
      symbol_list.checkPass(warnings, valid_symbol_list_types);
    }
  catch (SymbolList::SymbolListException &e)
    {
      cerr << "ERROR: rplot: " << e.message << endl;
      exit(EXIT_FAILURE);
    }
}

void
RplotStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  symbol_list.writeOutput("var_list_", output);
  output << "rplot(var_list_);" << endl;
}

void
RplotStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "rplot")";
  if (!symbol_list.empty())
    {
      output << ", ";
      symbol_list.writeJsonOutput(output);
    }
  output << "}";
}

void
UnitRootVarsStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "options_.diffuse_filter = 1;" << endl
         << "options_.steadystate.nocheck = 1;" << endl;
}

void
UnitRootVarsStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "unit_root_vars", )"
         << R"("diffuse_filter": 1, )"
         << R"("steady_state.nocheck": 1})";
}

PeriodsStatement::PeriodsStatement(int periods_arg) : periods{periods_arg}
{
}

void
PeriodsStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "options_.periods = " << periods << ";" << endl;
}

void
PeriodsStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "periods", )"
         << R"("periods": )" << periods << "}";
}

DsampleStatement::DsampleStatement(int val1_arg) : val1{val1_arg}, val2{-1}
{
}

DsampleStatement::DsampleStatement(int val1_arg, int val2_arg) : val1{val1_arg}, val2{val2_arg}
{
}

void
DsampleStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  if (val2 < 0)
    output << "dsample(" << val1 << ");" << endl;
  else
    output << "dsample(" << val1 << ", " << val2 << ");" << endl;
}

void
DsampleStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "dsample", )"
         << R"("value1": )" << val1 << ", "
         << R"("value2": )" << val2 << "}";
}

EstimatedParamsStatement::EstimatedParamsStatement(vector<EstimationParams> estim_params_list_arg,
                                                   const SymbolTable &symbol_table_arg) :
  estim_params_list{move(estim_params_list_arg)},
  symbol_table{symbol_table_arg}
{
}

void
EstimatedParamsStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  for (const auto &it : estim_params_list)
    {
      if (it.name == "dsge_prior_weight")
        mod_file_struct.dsge_prior_weight_in_estimated_params = true;

      // Handle case of degenerate beta prior
      if (it.prior == PriorDistributions::beta)
        try
          {
            if (it.mean->eval(eval_context_t()) == 0.5
                && it.std->eval(eval_context_t()) == 0.5)
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

  // Check that no parameter/endogenous is declared twice in the block
  set<string> already_declared;
  set<pair<string, string>> already_declared_corr;
  for (const auto &it : estim_params_list)
    {
      if (it.type == 3) // Correlation
        {
          // Use lexical ordering for the pair of symbols
          auto x = it.name < it.name2 ? pair(it.name, it.name2) : pair(it.name2, it.name);

          if (already_declared_corr.find(x) == already_declared_corr.end())
            already_declared_corr.insert(x);
          else
            {
              cerr << "ERROR: in `estimated_params' block, the correlation between " << it.name << " and " << it.name2 << " is declared twice." << endl;
              exit(EXIT_FAILURE);
            }
        }
      else
        {
          if (already_declared.find(it.name) == already_declared.end())
            already_declared.insert(it.name);
          else
            {
              cerr << "ERROR: in `estimated_params' block, the symbol " << it.name << " is declared twice." << endl;
              exit(EXIT_FAILURE);
            }
        }
    }

  // Fill in mod_file_struct.estimated_parameters (related to #469)
  for (const auto &it : estim_params_list)
    if (it.type == 2 && it.name != "dsge_prior_weight")
      mod_file_struct.estimated_parameters.insert(symbol_table.getID(it.name));
}

void
EstimatedParamsStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "estim_params_.var_exo = zeros(0, 10);" << endl
         << "estim_params_.var_endo = zeros(0, 10);" << endl
         << "estim_params_.corrx = zeros(0, 11);" << endl
         << "estim_params_.corrn = zeros(0, 11);" << endl
         << "estim_params_.param_vals = zeros(0, 10);" << endl;

  for (const auto &it : estim_params_list)
    {
      int symb_id = symbol_table.getTypeSpecificID(it.name) + 1;
      SymbolType symb_type = symbol_table.getType(it.name);

      switch (it.type)
        {
        case 1:
          if (symb_type == SymbolType::exogenous)
            output << "estim_params_.var_exo = [estim_params_.var_exo; ";
          else if (symb_type == SymbolType::endogenous)
            output << "estim_params_.var_endo = [estim_params_.var_endo; ";
          output << symb_id;
          break;
        case 2:
          output << "estim_params_.param_vals = [estim_params_.param_vals; "
                 << symb_id;
          break;
        case 3:
          if (symb_type == SymbolType::exogenous)
            output << "estim_params_.corrx = [estim_params_.corrx; ";
          else if (symb_type == SymbolType::endogenous)
            output << "estim_params_.corrn = [estim_params_.corrn; ";
          output << symb_id << ", " << symbol_table.getTypeSpecificID(it.name2)+1;
          break;
        }
      output << ", ";
      it.init_val->writeOutput(output);
      output << ", ";
      it.low_bound->writeOutput(output);
      output << ", ";
      it.up_bound->writeOutput(output);
      output << ", "
             << static_cast<int>(it.prior) << ", ";
      it.mean->writeOutput(output);
      output << ", ";
      it.std->writeOutput(output);
      output << ", ";
      it.p3->writeOutput(output);
      output << ", ";
      it.p4->writeOutput(output);
      output << ", ";
      it.jscale->writeOutput(output);
      output << " ];" << endl;
    }
}

void
EstimatedParamsStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "estimated_params", )"
         << R"("params": [)";
  for (auto it = estim_params_list.begin(); it != estim_params_list.end(); ++it)
    {
      if (it != estim_params_list.begin())
        output << ", ";
      output << "{";
      switch (it->type)
        {
        case 1:
          output << R"("var": ")" << it->name << R"(")";
          break;
        case 2:
          output << R"("param": ")" << it->name << R"(")";
          break;
        case 3:
          output << R"("var1": ")" << it->name << R"(",)"
                 << R"("var2": ")" << it->name2 << R"(")";
          break;
        }

      output << R"(, "init_val": ")";
      it->init_val->writeJsonOutput(output, {}, {});
      output << R"(", "lower_bound": ")";
      it->low_bound->writeJsonOutput(output, {}, {});
      output << R"(", "upper_bound": ")";
      it->up_bound->writeJsonOutput(output, {}, {});
      output << R"(", "prior_distribution": )"
             << static_cast<int>(it->prior)
             << R"(, "mean": ")";
      it->mean->writeJsonOutput(output, {}, {});
      output << R"(", "std": ")";
      it->std->writeJsonOutput(output, {}, {});
      output << R"(", "p3": ")";
      it->p3->writeJsonOutput(output, {}, {});
      output << R"(", "p4": ")";
      it->p4->writeJsonOutput(output, {}, {});
      output << R"(", "jscale": ")";
      it->jscale->writeJsonOutput(output, {}, {});
      output << R"("})" << endl;
    }
  output << "]"
         << "}";
}

EstimatedParamsInitStatement::EstimatedParamsInitStatement(vector<EstimationParams> estim_params_list_arg,
                                                           const SymbolTable &symbol_table_arg,
                                                           const bool use_calibration_arg) :
  estim_params_list{move(estim_params_list_arg)},
  symbol_table{symbol_table_arg},
  use_calibration{use_calibration_arg}
{
}

void
EstimatedParamsInitStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  if (use_calibration)
    mod_file_struct.estim_params_use_calib = true;
}

void
EstimatedParamsInitStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  if (use_calibration)
    output << "options_.use_calibration_initialization = 1;" << endl;

  bool skipline = false;

  for (const auto &it : estim_params_list)
    {
      int symb_id = symbol_table.getTypeSpecificID(it.name) + 1;
      SymbolType symb_type = symbol_table.getType(it.name);

      if (it.type < 3)
        {
          if (symb_type == SymbolType::exogenous)
            {
              output << "tmp1 = find(estim_params_.var_exo(:,1)==" << symb_id << ");" << endl;
              output << "if isempty(tmp1)" << endl;
              output << "    disp(sprintf('The standard deviation of %s is not estimated (the value provided in estimated_params_init is not used).', M_.exo_names{" << symb_id << "}))" << endl;
              skipline = true;
              output << "else" << endl;
              output << "    estim_params_.var_exo(tmp1,2) = ";
              it.init_val->writeOutput(output);
              output << ";" << endl;
              output << "end" << endl;
            }
          else if (symb_type == SymbolType::endogenous)
            {
              output << "tmp1 = find(estim_params_.var_endo(:,1)==" << symb_id << ");" << endl;
              output << "if isempty(tmp1)" << endl;
              output << "    disp(sprintf('The standard deviation of the measurement error on %s is not estimated (the value provided in estimated_params_init is not used).', M_.endo_names{" << symb_id << "}))" << endl;
              skipline = true;
              output << "else" << endl;
              output << "    estim_params_.var_endo(tmp1,2) = ";
              it.init_val->writeOutput(output);
              output << ";" << endl;
              output << "end" << endl;
            }
          else if (symb_type == SymbolType::parameter)
            {
              output << "tmp1 = find(estim_params_.param_vals(:,1)==" << symb_id << ");" << endl;
              output << "if isempty(tmp1)" << endl;
              output << "    disp(sprintf('Parameter %s is not estimated (the value provided in estimated_params_init is not used).', M_.param_names{" << symb_id << "}))" << endl;
              skipline = true;
              output << "else" << endl;
              output << "    estim_params_.param_vals(tmp1,2) = ";
              it.init_val->writeOutput(output);
              output << ";" << endl;
              output << "end" << endl;
            }
        }
      else
        {
          if (symb_type == SymbolType::exogenous)
            {
              output << "tmp1 = find((estim_params_.corrx(:,1)==" << symb_id << " & estim_params_.corrx(:,2)==" << symbol_table.getTypeSpecificID(it.name2)+1 << ") | "
                     <<             "(estim_params_.corrx(:,2)==" << symb_id << " & estim_params_.corrx(:,1)==" << symbol_table.getTypeSpecificID(it.name2)+1 << "));" << endl;
              output << "if isempty(tmp1)" << endl;
              output << "    disp(sprintf('The correlation between %s and %s is not estimated (the value provided in estimated_params_init is not used).', M_.exo_names{"
                     << symb_id << "}, M_.exo_names{" << symbol_table.getTypeSpecificID(it.name2)+1 << "}))" << endl;
              skipline = true;
              output << "else" << endl;
              output << "    estim_params_.corrx(tmp1,3) = ";
              it.init_val->writeOutput(output);
              output << ";" << endl;
              output << "end" << endl;
            }
          else if (symb_type == SymbolType::endogenous)
            {
              output << "tmp1 = find((estim_params_.corrn(:,1)==" << symb_id << " & estim_params_.corrn(:,2)==" << symbol_table.getTypeSpecificID(it.name2)+1 << ") | "
                     <<             "(estim_params_.corrn(:,2)==" << symb_id << " & estim_params_.corrn(:,1)==" << symbol_table.getTypeSpecificID(it.name2)+1 << "));" << endl;
              output << "if isempty(tmp1)" << endl;
              output << "    disp(sprintf('The correlation between measurement errors on %s and %s is not estimated (the value provided in estimated_params_init is not used).', M_.endo_names{"
                     << symb_id << "}, M_.endo_names{" << symbol_table.getTypeSpecificID(it.name2)+1 << "}))" << endl;
              skipline = true;
              output << "else" << endl;
              output << "    estim_params_.corrn(tmp1,3) = ";
              it.init_val->writeOutput(output);
              output << ";" << endl;
              output << "end" << endl;
            }
        }
    }
  if (skipline == true)
    output << "skipline()" << endl;
}

void
EstimatedParamsInitStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "estimated_params_init")";

  if (use_calibration)
    output << R"(, "use_calibration_initialization": 1)";

  output << R"(, "params": [)";
  for (auto it = estim_params_list.begin(); it != estim_params_list.end(); ++it)
    {
      if (it != estim_params_list.begin())
        output << ", ";
      output << "{";
      switch (it->type)
        {
        case 1:
          output << R"("var": ")" << it->name << R"(")";
          break;
        case 2:
          output << R"("param": ")" << it->name << R"(")";
          break;
        case 3:
          output << R"("var1": ")" << it->name << R"(",)"
                 << R"("var2": ")" << it->name2 << R"(")";
          break;
        }
      output << R"(, "init_val": ")";
      it->init_val->writeJsonOutput(output, {}, {});
      output << R"("})";
    }
  output << "]"
         << "}";
}

EstimatedParamsBoundsStatement::EstimatedParamsBoundsStatement(vector<EstimationParams> estim_params_list_arg,
                                                               const SymbolTable &symbol_table_arg) :
  estim_params_list{move(estim_params_list_arg)},
  symbol_table{symbol_table_arg}
{
}

void
EstimatedParamsBoundsStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  for (const auto &it : estim_params_list)
    {
      int symb_id = symbol_table.getTypeSpecificID(it.name) + 1;
      SymbolType symb_type = symbol_table.getType(it.name);

      if (it.type < 3)
        {
          if (symb_type == SymbolType::exogenous)
            {
              output << "tmp1 = find(estim_params_.var_exo(:,1)==" << symb_id << ");" << endl;

              output << "estim_params_.var_exo(tmp1,3) = ";
              it.low_bound->writeOutput(output);
              output << ";" << endl;

              output << "estim_params_.var_exo(tmp1,4) = ";
              it.up_bound->writeOutput(output);
              output << ";" << endl;
            }
          else if (symb_type == SymbolType::endogenous)
            {
              output << "tmp1 = find(estim_params_.var_endo(:,1)==" << symb_id << ");" << endl;

              output << "estim_params_.var_endo(tmp1,3) = ";
              it.low_bound->writeOutput(output);
              output << ";" << endl;

              output << "estim_params_.var_endo(tmp1,4) = ";
              it.up_bound->writeOutput(output);
              output << ";" << endl;
            }
          else if (symb_type == SymbolType::parameter)
            {
              output << "tmp1 = find(estim_params_.param_vals(:,1)==" << symb_id << ");" << endl;

              output << "estim_params_.param_vals(tmp1,3) = ";
              it.low_bound->writeOutput(output);
              output << ";" << endl;

              output << "estim_params_.param_vals(tmp1,4) = ";
              it.up_bound->writeOutput(output);
              output << ";" << endl;
            }
        }
      else
        {
          if (symb_type == SymbolType::exogenous)
            {
              output << "tmp1 = find((estim_params_.corrx(:,1)==" << symb_id << " & estim_params_.corrx(:,2)==" << symbol_table.getTypeSpecificID(it.name2)+1 << ") | "
                     <<             "(estim_params_.corrx(:,2)==" << symb_id << " & estim_params_.corrx(:,1)==" << symbol_table.getTypeSpecificID(it.name2)+1 << "));" << endl;

              output << "estim_params_.corrx(tmp1,4) = ";
              it.low_bound->writeOutput(output);
              output << ";" << endl;

              output << "estim_params_.corrx(tmp1,5) = ";
              it.up_bound->writeOutput(output);
              output << ";" << endl;
            }
          else if (symb_type == SymbolType::endogenous)
            {
              output << "tmp1 = find((estim_params_.corrn(:,1)==" << symb_id << " & estim_params_.corrn(:,2)==" << symbol_table.getTypeSpecificID(it.name2)+1 << ") | "
                     <<             "(estim_params_.corrn(:,2)==" << symb_id << " & estim_params_.corrn(:,1)==" << symbol_table.getTypeSpecificID(it.name2)+1 << "));" << endl;

              output << "estim_params_.corrn(tmp1,4) = ";
              it.low_bound->writeOutput(output);
              output << ";" << endl;

              output << "estim_params_.corrn(tmp1,5) = ";
              it.up_bound->writeOutput(output);
              output << ";" << endl;
            }
        }
    }
}

void
EstimatedParamsBoundsStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "estimated_params_bounds", )"
         << R"("params": [)";

  for (auto it = estim_params_list.begin(); it != estim_params_list.end(); ++it)
    {
      if (it != estim_params_list.begin())
        output << ", ";
      output << "{";
      switch (it->type)
        {
        case 1:
          output << R"("var": ")" << it->name << R"(")";
        case 2:
          output << R"("param": ")" << it->name << R"(")";
          break;
        case 3:
          output << R"("var1": ")" << it->name << R"(",)"
                 << R"("var2": ")" << it->name2 << R"(")";
          break;
        }
      output << R"(, "lower_bound": )";
      it->low_bound->writeJsonOutput(output, {}, {});
      output << R"(, "upper_bound": )";
      it->up_bound->writeJsonOutput(output, {}, {});
      output << "}";
    }
  output << "]"
         << "}";
}

ObservationTrendsStatement::ObservationTrendsStatement(trend_elements_t trend_elements_arg,
                                                       const SymbolTable &symbol_table_arg) :
  trend_elements{move(trend_elements_arg)},
  symbol_table{symbol_table_arg}
{
}

void
ObservationTrendsStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "options_.trend_coeff = {};" << endl;
  for (const auto &trend_element : trend_elements)
    {
      SymbolType type = symbol_table.getType(trend_element.first);
      if (type == SymbolType::endogenous)
        {
          output << "tmp1 = strmatch('" << trend_element.first << "',options_.varobs,'exact');" << endl;
          output << "options_.trend_coeffs{tmp1} = '";
          trend_element.second->writeOutput(output);
          output << "';" << endl;
        }
      else
        cerr << "Warning : Non-variable symbol used in observation_trends: " << trend_element.first << endl;
    }
}

void
ObservationTrendsStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "observation_trends", )"
         << R"("trends" : {)";
  bool printed = false;
  for (const auto &trend_element : trend_elements)
    {
      if (symbol_table.getType(trend_element.first) == SymbolType::endogenous)
        {
          if (printed)
            output << ", ";
          output << R"(")" << trend_element.first << R"(": ")";
          trend_element.second->writeJsonOutput(output, {}, {});
          output << R"(")" << endl;
          printed = true;
        }
      else
        cerr << "Warning : Non-variable symbol used in observation_trends: " << trend_element.first << endl;
    }
  output << "}"
         << "}";
}

OsrParamsStatement::OsrParamsStatement(SymbolList symbol_list_arg, const SymbolTable &symbol_table_arg) :
  symbol_list{move(symbol_list_arg)},
  symbol_table{symbol_table_arg}
{
}

void
OsrParamsStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  if (mod_file_struct.osr_params_present)
    cerr << "WARNING: You have more than one osr_params statement in the .mod file." << endl;
  mod_file_struct.osr_params_present = true;

  try
    {
      const vector<SymbolType> valid_symbol_list_types { SymbolType::parameter };
      symbol_list.checkPass(warnings, valid_symbol_list_types);
    }
  catch (SymbolList::SymbolListException &e)
    {
      cerr << "ERROR: osr: " << e.message << endl;
      exit(EXIT_FAILURE);
    }
}

void
OsrParamsStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  symbol_list.writeOutput("M_.osr.param_names", output);
  output << "M_.osr.param_names = cellstr(M_.osr.param_names);" << endl
         << "M_.osr.param_indices = zeros(length(M_.osr.param_names), 1);" << endl;
  int i = 0;
  vector<string> symbols = symbol_list.get_symbols();
  for (auto &symbol : symbols)
    output << "M_.osr.param_indices(" << ++i <<") = " << symbol_table.getTypeSpecificID(symbol) + 1 << ";" << endl;
}

void
OsrParamsStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "osr_params")";
  if (!symbol_list.empty())
    {
      output << ", ";
      symbol_list.writeJsonOutput(output);
    }
  output << "}";
}

OsrParamsBoundsStatement::OsrParamsBoundsStatement(vector<OsrParams> osr_params_list_arg) :
  osr_params_list{move(osr_params_list_arg)}
{
}

void
OsrParamsBoundsStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  if (!mod_file_struct.osr_params_present)
    {
      cerr << "ERROR: you must have an osr_params statement before the osr_params_bounds block." << endl;
      exit(EXIT_FAILURE);
    }
}

void
OsrParamsBoundsStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{

  output << "M_.osr.param_bounds = [-inf(length(M_.osr.param_names), 1), inf(length(M_.osr.param_names), 1)];" << endl;

  for (const auto &it : osr_params_list)
    {
      output << "M_.osr.param_bounds(strcmp(M_.osr.param_names, '" << it.name << "'), :) = [";
      it.low_bound->writeOutput(output);
      output << ", ";
      it.up_bound->writeOutput(output);
      output << "];" << endl;
    }
}

void
OsrParamsBoundsStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "osr_params_bounds")"
         << R"(, "bounds": [)";
  for (auto it = osr_params_list.begin(); it != osr_params_list.end(); ++it)
    {
      if (it != osr_params_list.begin())
        output << ", ";
      output << R"({"parameter": ")" << it->name << R"(",)"
             << R"("bounds": [")";
      it->low_bound->writeJsonOutput(output, {}, {});
      output << R"(", ")";
      it->up_bound->writeJsonOutput(output, {}, {});
      output << R"("])"
             << "}";
    }
  output << "]"
         << "}";
}

OsrStatement::OsrStatement(SymbolList symbol_list_arg,
                           OptionsList options_list_arg) :
  symbol_list{move(symbol_list_arg)},
  options_list{move(options_list_arg)}
{
}

void
OsrStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.osr_present = true;

  // Fill in option_order of mod_file_struct
  if (auto it = options_list.num_options.find("order");
      it != options_list.num_options.end())
    mod_file_struct.order_option = max(mod_file_struct.order_option, stoi(it->second));

  // Fill in mod_file_struct.partial_information
  if (auto it = options_list.num_options.find("partial_information");
      it != options_list.num_options.end() && it->second == "true")
    mod_file_struct.partial_information = true;

  // Option k_order_solver (implicit when order >= 3)
  if (auto it = options_list.num_options.find("k_order_solver");
      (it != options_list.num_options.end() && it->second == "true")
      || mod_file_struct.order_option >= 3)
    mod_file_struct.k_order_solver = true;

  try
    {
      const vector<SymbolType> valid_symbol_list_types { SymbolType::endogenous };
      symbol_list.checkPass(warnings, valid_symbol_list_types);
    }
  catch (SymbolList::SymbolListException &e)
    {
      cerr << "ERROR: osr: " << e.message << endl;
      exit(EXIT_FAILURE);
    }
}

void
OsrStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  // Ensure that order 3 implies k_order (#844)
  if (auto it = options_list.num_options.find("order"),
      it1 = options_list.num_options.find("k_order_solver");
      (it1 != options_list.num_options.end() && it1->second == "true")
      || (it != options_list.num_options.end() && stoi(it->second) >= 3))
    output << "options_.k_order_solver = true;" << endl;

  options_list.writeOutput(output);
  symbol_list.writeOutput("var_list_", output);
  output << "oo_.osr = osr(var_list_,M_.osr.param_names,M_.osr.variable_indices,M_.osr.variable_weights);" << endl;
}

void
OsrStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "osr")";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  if (!symbol_list.empty())
    {
      output << ", ";
      symbol_list.writeJsonOutput(output);
    }
  output << "}";
}

OptimWeightsStatement::OptimWeightsStatement(var_weights_t var_weights_arg,
                                             covar_weights_t covar_weights_arg,
                                             const SymbolTable &symbol_table_arg) :
  var_weights{move(var_weights_arg)},
  covar_weights{move(covar_weights_arg)},
  symbol_table{symbol_table_arg}
{
}

void
OptimWeightsStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.optim_weights_present = true;
}

void
OptimWeightsStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "%" << endl
         << "% OPTIM_WEIGHTS" << endl
         << "%" << endl
         << "M_.osr.variable_weights = sparse(M_.endo_nbr,M_.endo_nbr);" << endl
         << "M_.osr.variable_indices = [];" << endl << endl;

  for (const auto & [name, value] : var_weights)
    {
      int id = symbol_table.getTypeSpecificID(name) + 1;
      output << "M_.osr.variable_weights(" << id << "," << id << ") = ";
      value->writeOutput(output);
      output << ";" << endl;
      output << "M_.osr.variable_indices = [M_.osr.variable_indices; " << id << "];" << endl;
    }

  for (const auto & [names, value] : covar_weights)
    {
      int id1 = symbol_table.getTypeSpecificID(names.first) + 1;
      int id2 = symbol_table.getTypeSpecificID(names.second) + 1;
      output << "M_.osr.variable_weights(" << id1 << "," << id2 << ") = ";
      value->writeOutput(output);
      output << ";" << endl;
      output << "M_.osr.variable_indices = [M_.osr.variable_indices; " << id1 << "; " << id2 << "];" << endl;
    }
}

void
OptimWeightsStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "optim_weights", )"
         << R"("weights": [)";
  for (auto it = var_weights.begin(); it != var_weights.end(); ++it)
    {
      if (it != var_weights.begin())
        output << ", ";
      output << R"({"name": ")" << it->first << R"(")"
             << R"(, "value": ")";
      it->second->writeJsonOutput(output, {}, {});
      output << R"("})";
    }

  for (auto it = covar_weights.begin(); it != covar_weights.end(); ++it)
    {
      if (it != covar_weights.begin() || !var_weights.empty())
        output << ", ";
      output << R"({"name1": ")" << it->first.first << R"(")"
             << R"(, "name2": ")" << it->first.second << R"(")"
             << R"(, "value": ")";
      it->second->writeJsonOutput(output, {}, {});
      output << R"("})";
    }
  output << "]"
         << "}";
}

DynaSaveStatement::DynaSaveStatement(SymbolList symbol_list_arg,
                                     string filename_arg) :
  symbol_list{move(symbol_list_arg)},
  filename{move(filename_arg)}
{
}

void
DynaSaveStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  try
    {
      const vector<SymbolType> valid_symbol_list_types { SymbolType::endogenous };
      symbol_list.checkPass(warnings, valid_symbol_list_types);
    }
  catch (SymbolList::SymbolListException &e)
    {
      cerr << "ERROR: dynasave: " << e.message << endl;
      exit(EXIT_FAILURE);
    }
}

void
DynaSaveStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  symbol_list.writeOutput("var_list_", output);
  output << "dynasave('" << filename
         << "',var_list_);" << endl;
}

void
DynaSaveStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "dynasave", )"
         << R"("filename": ")" << filename << R"(")";
  if (!symbol_list.empty())
    {
      output << ", ";
      symbol_list.writeJsonOutput(output);
    }
  output << "}";
}

DynaTypeStatement::DynaTypeStatement(SymbolList symbol_list_arg,
                                     string filename_arg) :
  symbol_list(move(symbol_list_arg)),
  filename(move(filename_arg))
{
}

void
DynaTypeStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  try
    {
      const vector<SymbolType> valid_symbol_list_types { SymbolType::endogenous };
      symbol_list.checkPass(warnings, valid_symbol_list_types);
    }
  catch (SymbolList::SymbolListException &e)
    {
      cerr << "ERROR: dynatype: " << e.message << endl;
      exit(EXIT_FAILURE);
    }
}

void
DynaTypeStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  symbol_list.writeOutput("var_list_", output);
  output << "dynatype('" << filename
         << "',var_list_);" << endl;
}

void
DynaTypeStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "dynatype", )"
         << R"("filename": ")" << filename << R"(")";
  if (!symbol_list.empty())
    {
      output << ", ";
      symbol_list.writeJsonOutput(output);
    }
  output << "}";
}

ModelComparisonStatement::ModelComparisonStatement(filename_list_t filename_list_arg,
                                                   OptionsList options_list_arg) :
  filename_list{move(filename_list_arg)},
  options_list{move(options_list_arg)}
{
}

void
ModelComparisonStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  options_list.writeOutput(output);

  output << "ModelNames_ = {};" << endl;
  output << "ModelPriors_ = [];" << endl;

  for (const auto &it : filename_list)
    {
      output << "ModelNames_ = { ModelNames_{:} '" << it.first << "'};" << endl;
      output << "ModelPriors_ = [ ModelPriors_ ; " << it.second << "];" << endl;
    }
  output << "oo_ = model_comparison(ModelNames_,ModelPriors_,oo_,options_,M_.fname);" << endl;
}

void
ModelComparisonStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "model_comparison")";
  if (!filename_list.empty())
    output << R"(, "filename_list": {)";

  for (auto it = filename_list.begin(); it != filename_list.end(); ++it)
    {
      if (it != filename_list.begin())
        output << ", ";
      output << R"("name": ")" << it->first << R"(")"
             << R"("prior": ")" << it->second << R"(")";
    }

  if (!filename_list.empty())
    output << "}";

  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }

  output << "}";
}

PlannerObjectiveStatement::PlannerObjectiveStatement(const StaticModel &model_tree_arg) :
  model_tree{model_tree_arg}
{
}

void
PlannerObjectiveStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  assert(model_tree.equation_number() == 1);
  if (model_tree.exoPresentInEqs())
    {
      cerr << "ERROR: You cannot include exogenous variables in the planner objective. Please "
           << "define an auxiliary endogenous variable like eps_aux=epsilon and use it instead "
           << "of the varexo." << endl;
      exit(EXIT_FAILURE);
    }
  mod_file_struct.planner_objective_present = true;
}

const StaticModel &
PlannerObjectiveStatement::getPlannerObjective() const
{
  return model_tree;
}

void
PlannerObjectiveStatement::computingPass()
{
  model_tree.computingPass(3, 0, {}, false, false, false);
  computing_pass_called = true;
}

void
PlannerObjectiveStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  model_tree.writeStaticFile(basename + ".objective", false, false, false, "", {}, {}, false);
}

void
PlannerObjectiveStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "planner_objective")"
         << ", ";
  if (computing_pass_called)
    model_tree.writeJsonComputingPassOutput(output, false);
  else
    model_tree.writeJsonOutput(output);

  output << "}";
}

BVARDensityStatement::BVARDensityStatement(int maxnlags_arg, OptionsList options_list_arg) :
  maxnlags{maxnlags_arg},
  options_list{move(options_list_arg)}
{
}

void
BVARDensityStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.bvar_present = true;
}

void
BVARDensityStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  options_list.writeOutput(output);
  output << "bvar_density(" << maxnlags << ");" << endl;
}

void
BVARDensityStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "bvar_density")";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

BVARForecastStatement::BVARForecastStatement(int nlags_arg, OptionsList options_list_arg) :
  nlags{nlags_arg},
  options_list{move(options_list_arg)}
{
}

void
BVARForecastStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.bvar_present = true;
}

void
BVARForecastStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  options_list.writeOutput(output);
  output << "bvar_forecast(" << nlags << ");" << endl;
}

void
BVARForecastStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "bvar_forecast")";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

SBVARStatement::SBVARStatement(OptionsList options_list_arg) :
  options_list{move(options_list_arg)}
{
}

void
SBVARStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.bvar_present = true;
}

void
SBVARStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  options_list.writeOutput(output);
  output << "sbvar(M_,options_);" << endl;
}

void
SBVARStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "sbvar")";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

MSSBVAREstimationStatement::MSSBVAREstimationStatement(OptionsList options_list_arg) :
  options_list{move(options_list_arg)}
{
}

void
MSSBVAREstimationStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.bvar_present = true;

  if (options_list.num_options.find("ms.create_init") == options_list.num_options.end())
    if (options_list.string_options.find("datafile") == options_list.string_options.end()
        || options_list.num_options.find("ms.initial_year") == options_list.num_options.end())
      {
        cerr << "ERROR: If you do not pass no_create_init to ms_estimation, "
             << "you must pass the datafile and initial_year options." << endl;
        exit(EXIT_FAILURE);
      }
}

void
MSSBVAREstimationStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "options_ = initialize_ms_sbvar_options(M_, options_);" << endl
         << "options_.datafile = '';" << endl;
  options_list.writeOutput(output);
  output << "[options_, oo_] = ms_estimation(M_, options_, oo_);" << endl;
}

void
MSSBVAREstimationStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "ms_sbvar_estimation")";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

MSSBVARSimulationStatement::MSSBVARSimulationStatement(OptionsList options_list_arg) :
  options_list(move(options_list_arg))
{
}

void
MSSBVARSimulationStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.bvar_present = true;
}

void
MSSBVARSimulationStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "options_ = initialize_ms_sbvar_options(M_, options_);" << endl;
  options_list.writeOutput(output);

  // Redeclare drop option if necessary
  if ((options_list.num_options.find("ms.mh_replic") != options_list.num_options.end()
       || options_list.num_options.find("ms.thinning_factor") != options_list.num_options.end())
      && (options_list.num_options.find("ms.drop") == options_list.num_options.end()))
    output << "options_.ms.drop = 0.1*options_.ms.mh_replic*options_.ms.thinning_factor;" << endl;

  output << "[options_, oo_] = ms_simulation(M_, options_, oo_);" << endl;
}

void
MSSBVARSimulationStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "ms_sbvar_simulation")";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

MSSBVARComputeMDDStatement::MSSBVARComputeMDDStatement(OptionsList options_list_arg) :
  options_list{move(options_list_arg)}
{
}

void
MSSBVARComputeMDDStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.bvar_present = true;
}

void
MSSBVARComputeMDDStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "options_ = initialize_ms_sbvar_options(M_, options_);" << endl;
  options_list.writeOutput(output);
  output << "[options_, oo_] = ms_compute_mdd(M_, options_, oo_);" << endl;
}

void
MSSBVARComputeMDDStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "ms_sbvar_compute_mdd")";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

MSSBVARComputeProbabilitiesStatement::MSSBVARComputeProbabilitiesStatement(OptionsList options_list_arg) :
  options_list{move(options_list_arg)}
{
}

void
MSSBVARComputeProbabilitiesStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.bvar_present = true;

  if (options_list.num_options.find("ms.real_time_smoothed_probabilities") != options_list.num_options.end()
      && options_list.num_options.find("ms.filtered_probabilities") != options_list.num_options.end())
    {
      cerr << "ERROR: You may only pass one of real_time_smoothed "
           << "and filtered_probabilities to ms_compute_probabilities." << endl;
      exit(EXIT_FAILURE);
    }
}

void
MSSBVARComputeProbabilitiesStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "options_ = initialize_ms_sbvar_options(M_, options_);" << endl;
  options_list.writeOutput(output);
  output << "[options_, oo_] = ms_compute_probabilities(M_, options_, oo_);" << endl;
}

void
MSSBVARComputeProbabilitiesStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "ms_sbvar_compute_probabilities")";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

MSSBVARIrfStatement::MSSBVARIrfStatement(SymbolList symbol_list_arg,
                                         OptionsList options_list_arg) :
  symbol_list{move(symbol_list_arg)},
  options_list{move(options_list_arg)}
{
}

void
MSSBVARIrfStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.bvar_present = true;

  bool regime_present = options_list.num_options.find("ms.regime") != options_list.num_options.end();
  bool regimes_present = options_list.num_options.find("ms.regimes") != options_list.num_options.end();
  bool filtered_probabilities_present = options_list.num_options.find("ms.filtered_probabilities") != options_list.num_options.end();

  if ((filtered_probabilities_present && regime_present)
      || (filtered_probabilities_present && regimes_present)
      || (regimes_present && regime_present))
    {
      cerr << "ERROR: You may only pass one of regime, regimes and "
           << "filtered_probabilities to ms_irf" << endl;
      exit(EXIT_FAILURE);
    }

  try
    {
      const vector<SymbolType> valid_symbol_list_types { SymbolType::endogenous };
      symbol_list.checkPass(warnings, valid_symbol_list_types);
    }
  catch (SymbolList::SymbolListException &e)
    {
      cerr << "ERROR: ms_irf: " << e.message << endl;
      exit(EXIT_FAILURE);
    }
}

void
MSSBVARIrfStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "options_ = initialize_ms_sbvar_options(M_, options_);" << endl;
  symbol_list.writeOutput("var_list_", output);
  options_list.writeOutput(output);
  output << "[options_, oo_] = ms_irf(var_list_,M_, options_, oo_);" << endl;
}

void
MSSBVARIrfStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "ms_sbvar_irf")";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  if (!symbol_list.empty())
    {
      output << ", ";
      symbol_list.writeJsonOutput(output);
    }
  output << "}";
}

MSSBVARForecastStatement::MSSBVARForecastStatement(OptionsList options_list_arg) :
  options_list{move(options_list_arg)}
{
}

void
MSSBVARForecastStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.bvar_present = true;

  if (options_list.num_options.find("ms.regimes") != options_list.num_options.end()
      && options_list.num_options.find("ms.regime") != options_list.num_options.end())
    {
      cerr << "ERROR: You may only pass one of regime and regimes to ms_forecast" << endl;
      exit(EXIT_FAILURE);
    }
}

void
MSSBVARForecastStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "options_ = initialize_ms_sbvar_options(M_, options_);" << endl;
  options_list.writeOutput(output);
  output << "[options_, oo_] = ms_forecast(M_, options_, oo_);" << endl;
}

void
MSSBVARForecastStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "ms_sbvar_forecast")";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

MSSBVARVarianceDecompositionStatement::MSSBVARVarianceDecompositionStatement(OptionsList options_list_arg) :
  options_list{move(options_list_arg)}
{
}

void
MSSBVARVarianceDecompositionStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.bvar_present = true;

  bool regime_present = options_list.num_options.find("ms.regime") != options_list.num_options.end();
  bool regimes_present = options_list.num_options.find("ms.regimes") != options_list.num_options.end();
  bool filtered_probabilities_present = options_list.num_options.find("ms.filtered_probabilities") != options_list.num_options.end();

  if ((filtered_probabilities_present && regime_present)
      || (filtered_probabilities_present && regimes_present)
      || (regimes_present && regime_present))
    {
      cerr << "ERROR: You may only pass one of regime, regimes and "
           << "filtered_probabilities to ms_variance_decomposition" << endl;
      exit(EXIT_FAILURE);
    }
}

void
MSSBVARVarianceDecompositionStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "options_ = initialize_ms_sbvar_options(M_, options_);" << endl;
  options_list.writeOutput(output);
  output << "[options_, oo_] = ms_variance_decomposition(M_, options_, oo_);" << endl;
}

void
MSSBVARVarianceDecompositionStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "ms_sbvar_variance_decomposition")";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

IdentificationStatement::IdentificationStatement(OptionsList options_list_arg)
  : options_list{move(options_list_arg)}
{
  if (auto it = options_list.num_options.find("max_dim_cova_group");
      it != options_list.num_options.end() && stoi(it->second) == 0)
    {
      cerr << "ERROR: The max_dim_cova_group option to identification only accepts integers > 0." << endl;
      exit(EXIT_FAILURE);
    }
}

void
IdentificationStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.identification_present = true;

  if (auto it = options_list.num_options.find("order");
      it != options_list.num_options.end())
    {
      int order = stoi(it->second);
      if (order < 1 || order > 3)
        {
          cerr << "ERROR: the order option of identification command must be between 1 and 3" << endl;

          exit(EXIT_FAILURE);
        }
      mod_file_struct.identification_order = max(mod_file_struct.identification_order, order);
    }
  else
    // The default value for order is 1 (which triggers 2nd order dynamic derivatives)
    mod_file_struct.identification_order = max(mod_file_struct.identification_order, 1);
}

void
IdentificationStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  options_list.writeOutput(output, "options_ident");

  /* Ensure that nograph, nodisplay and graph_format are also set in top-level
     options_.
     \todo factorize this code between identification and dynare_sensitivity,
     and provide a generic mechanism for this situation (maybe using regexps) */
  if (auto it = options_list.num_options.find("nodisplay");
      it != options_list.num_options.end())
    output << "options_.nodisplay = " << it->second << ";" << endl;
  if (auto it = options_list.num_options.find("nograph");
      it != options_list.num_options.end())
    output << "options_.nograph = " << it->second << ";" << endl;
  if (auto it = options_list.string_options.find("graph_format");
      it != options_list.string_options.end())
    output << "options_.graph_format = '" << it->second << "';" << endl;

  output << "dynare_identification(options_ident);" << endl;
}

void
IdentificationStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "identification")";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

WriteLatexDynamicModelStatement::WriteLatexDynamicModelStatement(const DynamicModel &dynamic_model_arg, bool write_equation_tags_arg) :
  dynamic_model{dynamic_model_arg},
  write_equation_tags{write_equation_tags_arg}
{
}

void
WriteLatexDynamicModelStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  dynamic_model.writeLatexFile(basename, write_equation_tags);
}

void
WriteLatexDynamicModelStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "write_latex_dynamic_model"})";
}

WriteLatexStaticModelStatement::WriteLatexStaticModelStatement(const StaticModel &static_model_arg, bool write_equation_tags_arg) :
  static_model(static_model_arg),
  write_equation_tags(write_equation_tags_arg)
{
}

void
WriteLatexStaticModelStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  static_model.writeLatexFile(basename, write_equation_tags);
}

void
WriteLatexStaticModelStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "write_latex_static_model"})";
}

WriteLatexOriginalModelStatement::WriteLatexOriginalModelStatement(const DynamicModel &original_model_arg, bool write_equation_tags_arg) :
  original_model{original_model_arg},
  write_equation_tags{write_equation_tags_arg}
{
}

void
WriteLatexOriginalModelStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  original_model.writeLatexOriginalFile(basename, write_equation_tags);
}

void
WriteLatexOriginalModelStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "write_latex_original_model"})";
}

WriteLatexSteadyStateModelStatement::WriteLatexSteadyStateModelStatement(const SteadyStateModel &steady_state_model_arg) :
  steady_state_model{steady_state_model_arg}
{
}

void
WriteLatexSteadyStateModelStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.write_latex_steady_state_model_present = true;
}

void
WriteLatexSteadyStateModelStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  steady_state_model.writeLatexSteadyStateFile(basename);
}

void
WriteLatexSteadyStateModelStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "write_latex_steady_state_model"})";
}

ShockDecompositionStatement::ShockDecompositionStatement(SymbolList symbol_list_arg,
                                                         OptionsList options_list_arg) :
  symbol_list{move(symbol_list_arg)},
  options_list{move(options_list_arg)}
{
}

void
ShockDecompositionStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  if (auto it = options_list.num_options.find("shock_decomp.with_epilogue");
      it != options_list.num_options.end() && it->second == "true")
    mod_file_struct.with_epilogue_option = true;

  try
    {
      const vector<SymbolType> valid_symbol_list_types { SymbolType::endogenous };
      symbol_list.checkPass(warnings, valid_symbol_list_types);
    }
  catch (SymbolList::SymbolListException &e)
    {
      cerr << "ERROR: shock_decomposition: " << e.message << endl;
      exit(EXIT_FAILURE);
    }
}

void
ShockDecompositionStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  options_list.writeOutput(output);
  symbol_list.writeOutput("var_list_", output);
  output << "oo_ = shock_decomposition(M_,oo_,options_,var_list_,bayestopt_,estim_params_);" << endl;
}

void
ShockDecompositionStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "shock_decomposition")";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  if (!symbol_list.empty())
    {
      output << ", ";
      symbol_list.writeJsonOutput(output);
    }
  output << "}";
}

RealtimeShockDecompositionStatement::RealtimeShockDecompositionStatement(SymbolList symbol_list_arg,
                                                                         OptionsList options_list_arg) :
  symbol_list{move(symbol_list_arg)},
  options_list{move(options_list_arg)}
{
}

void
RealtimeShockDecompositionStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  if (auto it = options_list.num_options.find("shock_decomp.with_epilogue");
      it != options_list.num_options.end() && it->second == "true")
    mod_file_struct.with_epilogue_option = true;

  try
    {
      const vector<SymbolType> valid_symbol_list_types { SymbolType::endogenous };
      symbol_list.checkPass(warnings, valid_symbol_list_types);
    }
  catch (SymbolList::SymbolListException &e)
    {
      cerr << "ERROR: realtime_shock_decomposition: " << e.message << endl;
      exit(EXIT_FAILURE);
    }
}

void
RealtimeShockDecompositionStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  options_list.writeOutput(output);
  symbol_list.writeOutput("var_list_", output);
  output << "oo_ = realtime_shock_decomposition(M_,oo_,options_,var_list_,bayestopt_,estim_params_);" << endl;
}

void
RealtimeShockDecompositionStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "realtime_shock_decomposition")";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  if (!symbol_list.empty())
    {
      output << ", ";
      symbol_list.writeJsonOutput(output);
    }
  output << "}";
}

PlotShockDecompositionStatement::PlotShockDecompositionStatement(SymbolList symbol_list_arg,
                                                                 OptionsList options_list_arg) :
  symbol_list{move(symbol_list_arg)},
  options_list{move(options_list_arg)}
{
}

void
PlotShockDecompositionStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  try
    {
      const vector<SymbolType> valid_symbol_list_types { SymbolType::endogenous,
                                                         SymbolType::epilogue };
      symbol_list.checkPass(warnings, valid_symbol_list_types);
    }
  catch (SymbolList::SymbolListException &e)
    {
      cerr << "ERROR: plot_shock_decomposition: " << e.message << endl;
      exit(EXIT_FAILURE);
    }
}

void
PlotShockDecompositionStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "options_ = set_default_plot_shock_decomposition_options(options_);" << endl;
  options_list.writeOutput(output);
  symbol_list.writeOutput("var_list_", output);
  output << "oo_ = plot_shock_decomposition(M_, oo_, options_, var_list_);" << endl;
}

void
PlotShockDecompositionStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "plot_shock_decomposition")";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  if (!symbol_list.empty())
    {
      output << ", ";
      symbol_list.writeJsonOutput(output);
    }
  output << "}";
}

InitialConditionDecompositionStatement::InitialConditionDecompositionStatement(SymbolList symbol_list_arg,
                                                                               OptionsList options_list_arg) :
  symbol_list{move(symbol_list_arg)},
  options_list{move(options_list_arg)}
{
}

void
InitialConditionDecompositionStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  if (auto it = options_list.num_options.find("initial_condition_decomp.with_epilogue");
      it != options_list.num_options.end() && it->second == "true")
    mod_file_struct.with_epilogue_option = true;

  try
    {
      const vector<SymbolType> valid_symbol_list_types { SymbolType::endogenous };
      symbol_list.checkPass(warnings, valid_symbol_list_types);
    }
  catch (SymbolList::SymbolListException &e)
    {
      cerr << "ERROR: initial_condition_decomposition: " << e.message << endl;
      exit(EXIT_FAILURE);
    }
}

void
InitialConditionDecompositionStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "options_ = set_default_initial_condition_decomposition_options(options_);" << endl;
  options_list.writeOutput(output);
  symbol_list.writeOutput("var_list_", output);
  output << "oo_ = initial_condition_decomposition(M_, oo_, options_, var_list_, bayestopt_, estim_params_);" << endl;
}

void
InitialConditionDecompositionStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "initial_condition_decomposition")";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  if (!symbol_list.empty())
    {
      output << ", ";
      symbol_list.writeJsonOutput(output);
    }
  output << "}";
}

SqueezeShockDecompositionStatement::SqueezeShockDecompositionStatement(SymbolList symbol_list_arg)
  : symbol_list{move(symbol_list_arg)}
{
}

void
SqueezeShockDecompositionStatement::checkPass(ModFileStructure &mod_file_struct,
                                              WarningConsolidation &warnings)
{
  try
    {
      const vector<SymbolType> valid_symbol_list_types { SymbolType::endogenous };
      symbol_list.checkPass(warnings, valid_symbol_list_types);
    }
  catch (SymbolList::SymbolListException &e)
    {
      cerr << "ERROR: squeeze_shock_decomposition: " << e.message << endl;
      exit(EXIT_FAILURE);
    }
}

void
SqueezeShockDecompositionStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  if (symbol_list.empty())
    output << "oo_ = squeeze_shock_decomposition(M_, oo_, options_);" << endl;
  else
    {
      symbol_list.writeOutput("var_list_", output);
      output << "oo_ = squeeze_shock_decomposition(M_, oo_, options_, var_list_);" << endl;
    }
}

void
SqueezeShockDecompositionStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "squeeze_shock_decomposition")";
  if (!symbol_list.empty())
    {
      output << ", ";
      symbol_list.writeJsonOutput(output);
    }
  output << "}";
}

ConditionalForecastStatement::ConditionalForecastStatement(OptionsList options_list_arg) :
  options_list{move(options_list_arg)}
{
}

void
ConditionalForecastStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  if (options_list.string_options.find("parameter_set") == options_list.string_options.end())
    {
      cerr << "ERROR: You must pass the `parameter_set` option to conditional_forecast" << endl;
      exit(EXIT_FAILURE);
    }
}

void
ConditionalForecastStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  options_list.writeOutput(output, "options_cond_fcst_");
  output << "imcforecast(constrained_paths_, constrained_vars_, options_cond_fcst_);" << endl;
}

void
ConditionalForecastStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "conditional_forecast")";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

PlotConditionalForecastStatement::PlotConditionalForecastStatement(int periods_arg, SymbolList symbol_list_arg) :
  periods{periods_arg},
  symbol_list{move(symbol_list_arg)}
{
}

void
PlotConditionalForecastStatement::checkPass(ModFileStructure &mod_file_struct,
                                            WarningConsolidation &warnings)
{
  try
    {
      const vector<SymbolType> valid_symbol_list_types { SymbolType::endogenous };
      symbol_list.checkPass(warnings, valid_symbol_list_types);
    }
  catch (SymbolList::SymbolListException &e)
    {
      cerr << "ERROR: plot_conditional_forecast: " << e.message << endl;
      exit(EXIT_FAILURE);
    }
}

void
PlotConditionalForecastStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  symbol_list.writeOutput("var_list_", output);
  if (periods == -1)
    output << "plot_icforecast(var_list_,[],options_,oo_);" << endl;
  else
    output << "plot_icforecast(var_list_, " << periods << ",options_,oo_);" << endl;
}

void
PlotConditionalForecastStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "plot_conditional_forecast", )"
         << R"("periods": )" << periods;
  if (!symbol_list.empty())
    {
      output << ", ";
      symbol_list.writeJsonOutput(output);
    }
  output << "}";
}

SvarIdentificationStatement::SvarIdentificationStatement(svar_identification_restrictions_t restrictions_arg,
                                                         bool upper_cholesky_present_arg,
                                                         bool lower_cholesky_present_arg,
                                                         bool constants_exclusion_present_arg,
                                                         const SymbolTable &symbol_table_arg) :
  restrictions{move(restrictions_arg)},
  upper_cholesky_present{upper_cholesky_present_arg},
  lower_cholesky_present{lower_cholesky_present_arg},
  constants_exclusion_present{constants_exclusion_present_arg},
  symbol_table{symbol_table_arg}
{
}

int
SvarIdentificationStatement::getMaxLag() const
{
  int max_lag = 0;
  for (const auto &restriction : restrictions)
    if (restriction.lag > max_lag)
      max_lag = restriction.lag;

  return max_lag;
}

void
SvarIdentificationStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  // no equations OK with Svar Identification
  mod_file_struct.bvar_present = true;
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
SvarIdentificationStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  assert(!(upper_cholesky_present && lower_cholesky_present));
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
      output << "options_.ms.Qi = cell(" << n << ",1);" << endl
             << "options_.ms.Ri = cell(" << n << ",1);" << endl;

      for (auto &it : restrictions)
        {
          assert(it.lag >= 0);
          if (it.lag == 0)
            output << "options_.ms.Qi{" << it.equation << "}(" << it.restriction_nbr << ", " << it.variable + 1 << ") = ";
          else
            {
              int col = (it.lag-1)*n+it.variable+1;
              if (col > k)
                {
                  cerr << "ERROR: lag =" << it.lag << ", num endog vars = " << n << "current endog var index = " << it.variable << ". Index "
                       << "out of bounds. If the above does not represent a logical error, please report this to the Dynare Team." << endl;
                  exit(EXIT_FAILURE);
                }
              output << "options_.ms.Ri{" << it.equation << "}(" << it.restriction_nbr << ", " << col << ") = ";
            }
          it.value->writeOutput(output);
          output << ";" << endl;
        }
      output << "options_.ms.nlags = " << r << ";" << endl;
    }
}

void
SvarIdentificationStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "svar_identification")";

  if (upper_cholesky_present)
    output << R"(, "upper_cholesky": 1)";

  if (lower_cholesky_present)
    output << R"(, "lower_cholesky": 1)";

  if (constants_exclusion_present)
    output << R"(, "constants_exclusion": 1)";

  if (!upper_cholesky_present && !lower_cholesky_present)
    {
      output << R"(, "nlags": )" << getMaxLag()
             << R"(, "restrictions": [)";

      for (auto it = restrictions.begin(); it != restrictions.end(); ++it)
        {
          if (it != restrictions.begin())
            output << ", ";
          output << "{"
                 << R"("equation_number": )" << it->equation << ", "
                 << R"("restriction_number": )" << it->restriction_nbr << ", "
                 << R"("variable": ")" << symbol_table.getName(it->variable) << R"(", )"
                 << R"("expression": ")";
          it->value->writeOutput(output);
          output << R"("})";
        }
      output << "]";
    }
  output << "}";
}

MarkovSwitchingStatement::MarkovSwitchingStatement(OptionsList options_list_arg) :
  options_list{move(options_list_arg)}
{
  if (auto it_num = options_list.num_options.find("ms.restrictions");
      it_num != options_list.num_options.end())
    {
      using namespace boost;
      auto it_num_regimes = options_list.num_options.find("ms.number_of_regimes");
      assert(it_num_regimes != options_list.num_options.end());
      auto num_regimes = stoi(it_num_regimes->second);

      vector<string> tokenizedRestrictions;
      split(tokenizedRestrictions, it_num->second, is_any_of("["), token_compress_on);
      for (auto &tokenizedRestriction : tokenizedRestrictions)
        if (tokenizedRestriction.size() > 0)
          {
            vector<string> restriction;
            split(restriction, tokenizedRestriction, is_any_of("], "));
            for (auto it1 = restriction.begin(); it1 != restriction.end();)
              if (it1->empty())
                restriction.erase(it1);
              else
                ++it1;

            if (restriction.size() != 3)
              {
                cerr << "ERROR: restrictions in the subsample statement must be specified in the form "
                     << "[current_period_regime, next_period_regime, transition_probability]" << endl;
                exit(EXIT_FAILURE);
              }

            try
              {
                auto from_regime = stoi(restriction[0]);
                auto to_regime = stoi(restriction[1]);
                if (from_regime > num_regimes || to_regime > num_regimes)
                  {
                    cerr << "ERROR: the regimes specified in the restrictions option must be "
                         << "<= the number of regimes specified in the number_of_regimes option" << endl;
                    exit(EXIT_FAILURE);
                  }

                if (restriction_map.find({ from_regime, to_regime }) !=
                    restriction_map.end())
                  {
                    cerr << "ERROR: two restrictions were given for: " << from_regime << ", "
                         << to_regime << endl;
                    exit(EXIT_FAILURE);
                  }

                auto transition_probability = stod(restriction[2]);
                if (transition_probability > 1.0)
                  {
                    cerr << "ERROR: the transition probability, " << transition_probability
                         << " must be less than 1" << endl;
                    exit(EXIT_FAILURE);
                  }
                restriction_map[{ from_regime, to_regime }] = transition_probability;
              }
            catch (const invalid_argument &)
              {
                cerr << "ERROR: The first two arguments for a restriction must be integers "
                     << "specifying the regime and the last must be a double specifying the "
                     << "transition probability. You wrote [" << tokenizedRestriction << endl;
                exit(EXIT_FAILURE);
              }
          }
    }
}

void
MarkovSwitchingStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  auto itChain = options_list.num_options.find("ms.chain");
  assert(itChain != options_list.num_options.end());
  int chainNumber = stoi(itChain->second);
  if (++mod_file_struct.last_markov_switching_chain != chainNumber)
    {
      cerr << "ERROR: The markov_switching chain option takes consecutive integers "
           << "beginning at 1." << endl;
      exit(EXIT_FAILURE);
    }

  if (auto it_num = options_list.num_options.find("ms.restrictions");
      it_num != options_list.num_options.end())
    {
      using namespace boost;
      auto it_num_regimes = options_list.num_options.find("ms.number_of_regimes");
      assert(it_num_regimes != options_list.num_options.end());
      auto num_regimes = stoi(it_num_regimes->second);
      vector<double> col_trans_prob_sum(num_regimes, 0);
      vector<double> row_trans_prob_sum(num_regimes, 0);
      vector<bool> all_restrictions_in_row(num_regimes, true);
      vector<bool> all_restrictions_in_col(num_regimes, true);
      for (int row = 0; row < num_regimes; row++)
        for (int col = 0; col < num_regimes; col++)
          if (restriction_map.find({ row+1, col+1 }) != restriction_map.end())
            {
              row_trans_prob_sum[row] += restriction_map[{ row+1, col+1 }];
              col_trans_prob_sum[col] += restriction_map[{ row+1, col+1 }];
            }
          else
            {
              all_restrictions_in_row[row] = false;
              all_restrictions_in_col[col] = false;
            }

      for (int i = 0; i < num_regimes; i++)
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

  if (options_list.symbol_list_options.find("ms.parameters") != options_list.symbol_list_options.end())
    mod_file_struct.ms_dsge_present = true;
}

void
MarkovSwitchingStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  bool isDurationAVec = true;
  string infStr("Inf");
  OptionsList::num_options_t::const_iterator itChain, itNOR, itDuration;
  map<pair<int, int>, double >::const_iterator itR;

  itChain = options_list.num_options.find("ms.chain");
  assert(itChain != options_list.num_options.end());

  itDuration = options_list.num_options.find("ms.duration");
  assert(itDuration != options_list.num_options.end());
  if (stod(itDuration->second) || infStr.compare(itDuration->second) == 0)
    isDurationAVec = false;
  output << "options_.ms.duration = " << itDuration->second << ";" << endl;

  itNOR = options_list.num_options.find("ms.number_of_regimes");
  assert(itNOR != options_list.num_options.end());
  for (int i = 0; i < stoi(itNOR->second); i++)
    {
      output << "options_.ms.ms_chain(" << itChain->second << ").regime("
             << i+1 << ").duration = options_.ms.duration";
      if (isDurationAVec)
        output << "(" << i+1 << ")";
      output << ";" << endl;
    }

  int restrictions_index = 0;
  for (itR = restriction_map.begin(); itR != restriction_map.end(); itR++)
    output << "options_.ms.ms_chain(" << itChain->second << ").restrictions("
           << ++restrictions_index << ") = {[" << itR->first.first << ", "
           << itR->first.second << ", " << itR->second << "]};" << endl;
}

void
MarkovSwitchingStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "markov_switching")";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }

  if (!restriction_map.empty())
    output << ", {";
  for (auto it = restriction_map.begin(); it != restriction_map.end(); ++it)
    {
      if (it != restriction_map.begin())
        output << ", ";
      output << R"({"current_period_regime": )" << it->first.first
             << R"(, "next_period_regime": )" << it->first.second
             << R"(, "transition_probability": )"<< it->second
             << "}";
    }
  if (!restriction_map.empty())
    output << "}";
  output << "}";
}

SvarStatement::SvarStatement(OptionsList options_list_arg) :
  options_list{move(options_list_arg)}
{
}

void
SvarStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  auto it0 = options_list.string_options.find("ms.coefficients"),
    it1 = options_list.string_options.find("ms.variances"),
    it2 = options_list.string_options.find("ms.constants");
  assert((it0 != options_list.string_options.end()
          && it1 == options_list.string_options.end()
          && it2 == options_list.string_options.end())
         || (it0 == options_list.string_options.end()
             && it1 != options_list.string_options.end()
             && it2 == options_list.string_options.end())
         || (it0 == options_list.string_options.end()
             && it1 == options_list.string_options.end()
             && it2 != options_list.string_options.end()));
}

void
SvarStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  auto it = options_list.num_options.find("ms.chain");
  assert(it != options_list.num_options.end());
  output << "options_.ms.ms_chain(" << it->second << ")";

  if (auto it0 = options_list.string_options.find("ms.coefficients");
      it0 != options_list.string_options.end())
    output << "." << it0->second;
  else if (auto it1 = options_list.string_options.find("ms.variances");
           it1 != options_list.string_options.end())
    output << "." << it1->second;
  else
    output << "." << options_list.string_options.find("ms.constants")->second;

  output << ".equations = ";
  if (auto itv = options_list.vector_int_options.find("ms.equations");
      itv != options_list.vector_int_options.end())
    {
      assert(itv->second.size() >= 1);
      if (itv->second.size() > 1)
        {
          output << "[";
          for (int viit : itv->second)
            output << viit << ";";
          output << "];" << endl;
        }
      else
        output << itv->second.front() << ";" << endl;
    }
  else
    output << "'ALL';" << endl;
}

void
SvarStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "svar")";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

void
SvarGlobalIdentificationCheckStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "svar_global_identification_check(options_);" << std::endl;
}

void
SvarGlobalIdentificationCheckStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "svar_global_identification"})";
}

SetTimeStatement::SetTimeStatement(OptionsList options_list_arg) :
  options_list{move(options_list_arg)}
{
}

void
SetTimeStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  options_list.writeOutput(output);
}

void
SetTimeStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "set_time")";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

EstimationDataStatement::EstimationDataStatement(OptionsList options_list_arg) :
  options_list{move(options_list_arg)}
{
}

void
EstimationDataStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.estimation_data_statement_present = true;

  if (auto it = options_list.num_options.find("nobs");
      it != options_list.num_options.end())
    if (stoi(it->second) <= 0)
      {
        cerr << "ERROR: The nobs option of the data statement only accepts positive integers." << endl;
        exit(EXIT_FAILURE);
      }

  if (options_list.string_options.find("file") == options_list.string_options.end()
      && options_list.string_options.find("series") == options_list.string_options.end())
    {
      cerr << "ERROR: The file or series option must be passed to the data statement." << endl;
      exit(EXIT_FAILURE);
    }

  if (options_list.string_options.find("file") != options_list.string_options.end()
      && options_list.string_options.find("series") != options_list.string_options.end())
    {
      cerr << "ERROR: The file and series options cannot be used simultaneously in the data statement." << endl;
      exit(EXIT_FAILURE);
    }
}

void
EstimationDataStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  options_list.writeOutput(output, "options_.dataset");
}

void
EstimationDataStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "estimation_data")";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

SubsamplesStatement::SubsamplesStatement(string name1_arg,
                                         string name2_arg,
                                         subsample_declaration_map_t subsample_declaration_map_arg,
                                         const SymbolTable &symbol_table_arg) :
  name1{move(name1_arg)},
  name2{move(name2_arg)},
  subsample_declaration_map{move(subsample_declaration_map_arg)},
  symbol_table{symbol_table_arg}
{
}

void
SubsamplesStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
}

void
SubsamplesStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "subsamples_indx = get_new_or_existing_ei_index('subsamples_index', '"
         << name1 << "','" << name2 << "');" << endl
         << "estimation_info.subsamples_index(subsamples_indx) = {'" << name1;
  if (!name2.empty())
    output << ":" << name2;
  output << "'};" << endl
         << "estimation_info.subsamples(subsamples_indx).range = {};" << endl
         << "estimation_info.subsamples(subsamples_indx).range_index = {};" << endl;

  int map_indx = 1;
  for (auto it = subsample_declaration_map.begin();
       it != subsample_declaration_map.end(); ++it, map_indx++)
    output << "estimation_info.subsamples(subsamples_indx).range_index(" << map_indx << ") = {'"
           << it->first << "'};" << endl
           << "estimation_info.subsamples(subsamples_indx).range(" << map_indx << ").date1 = "
           << it->second.first << ";" << endl
           << "estimation_info.subsamples(subsamples_indx).range(" << map_indx << ").date2 = "
           << it->second.second << ";" << endl;

  // Initialize associated subsample substructures in estimation_info
  const SymbolType symb_type = symbol_table.getType(name1);
  string lhs_field;
  if (symb_type == SymbolType::parameter)
    lhs_field = "parameter";
  else if (symb_type == SymbolType::exogenous || symb_type == SymbolType::exogenousDet)
    lhs_field = "structural_innovation";
  else
    lhs_field = "measurement_error";

  output << "eifind = get_new_or_existing_ei_index('" << lhs_field;

  if (!name2.empty())
    output << "_corr";
  output << "_prior_index', '"
         << name1 << "', '";
  if (!name2.empty())
    output << name2;
  output << "');" << endl;

  lhs_field = "estimation_info." + lhs_field;
  if (!name2.empty())
    lhs_field += "_corr";
  output << lhs_field << "_prior_index(eifind) = {'" << name1;
  if (!name2.empty())
    output << ":" << name2;
  output << "'};" << endl;

  output << lhs_field << "(eifind).subsample_prior = estimation_info.empty_prior;" << endl
         << lhs_field << "(eifind).subsample_prior(1:" << subsample_declaration_map.size()
         << ") = estimation_info.empty_prior;" << endl
         << lhs_field << "(eifind).range_index = estimation_info.subsamples(subsamples_indx).range_index;"
         << endl;
}

void
SubsamplesStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "subsamples")"
         << R"(, "name1": ")" << name1 << R"(")";
  if (!name2.empty())
    output << R"(, "name2": ")" << name2 << R"(")";

  output << R"(, "declarations": {)";
  for (auto it = subsample_declaration_map.begin();
       it != subsample_declaration_map.end(); ++it)
    {
      if (it != subsample_declaration_map.begin())
        output << ",";
      output << "{"
             << R"("range_index": ")" << it->first << R"(")"
             << R"(, "date1": ")" << it->second.first << R"(")"
             << R"(, "date2": ")" << it->second.second << R"(")"
             << "}";
    }
  output << "}"
         << "}";
}

SubsamplesEqualStatement::SubsamplesEqualStatement(string to_name1_arg,
                                                   string to_name2_arg,
                                                   string from_name1_arg,
                                                   string from_name2_arg,
                                                   const SymbolTable &symbol_table_arg) :
  to_name1{move(to_name1_arg)},
  to_name2{move(to_name2_arg)},
  from_name1{move(from_name1_arg)},
  from_name2{move(from_name2_arg)},
  symbol_table{symbol_table_arg}
{
}

void
SubsamplesEqualStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "subsamples_to_indx = get_new_or_existing_ei_index('subsamples_index', '"
         << to_name1 << "','" << to_name2 << "');" << endl
         << "estimation_info.subsamples_index(subsamples_to_indx) = {'" << to_name1;
  if (!to_name2.empty())
    output << ":" << to_name2;
  output << "'};" << endl
         << "subsamples_from_indx = get_existing_subsamples_indx('" << from_name1 << "','" << from_name2 << "');"
         << endl
         << "estimation_info.subsamples(subsamples_to_indx) = estimation_info.subsamples(subsamples_from_indx);"
         << endl;

  // Initialize associated subsample substructures in estimation_info
  const SymbolType symb_type = symbol_table.getType(to_name1);
  string lhs_field;
  if (symb_type == SymbolType::parameter)
    lhs_field = "parameter";
  else if (symb_type == SymbolType::exogenous || symb_type == SymbolType::exogenousDet)
    lhs_field = "structural_innovation";
  else
    lhs_field = "measurement_error";

  output << "eifind = get_new_or_existing_ei_index('" << lhs_field;

  if (!to_name2.empty())
    output << "_corr";
  output << "_prior_index', '"
         << to_name1 << "', '";
  if (!to_name2.empty())
    output << to_name2;
  output << "');" << endl;

  lhs_field = "estimation_info." + lhs_field;
  if (!to_name2.empty())
    lhs_field += "_corr";
  output << lhs_field << "_prior_index(eifind) = {'" << to_name1;
  if (!to_name2.empty())
    output << ":" << to_name2;
  output << "'};" << endl;

  output << lhs_field << "(eifind).subsample_prior = estimation_info.empty_prior;" << endl
         << lhs_field << "(eifind).subsample_prior(1:size(estimation_info.subsamples(subsamples_to_indx).range_index,2)) = estimation_info.empty_prior;"
         << endl
         << lhs_field << "(eifind).range_index = estimation_info.subsamples(subsamples_to_indx).range_index;"
         << endl;
}

void
SubsamplesEqualStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "subsamples_equal")"
         << R"(, "to_name1": ")" << to_name1 << R"(")";
  if (!to_name2.empty())
    output << R"(, "to_name2": ")" << to_name2 << R"(")";
  output << R"(, "from_name1": ")" << from_name1 << R"(")";
  if (!from_name2.empty())
    output << R"(, "from_name2": ")" << from_name2 << R"(")";
  output << "}";
}

JointPriorStatement::JointPriorStatement(vector<string> joint_parameters_arg,
                                         PriorDistributions prior_shape_arg,
                                         OptionsList options_list_arg) :
  joint_parameters{move(joint_parameters_arg)},
  prior_shape{prior_shape_arg},
  options_list{move(options_list_arg)}
{
}

void
JointPriorStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  if (joint_parameters.size() < 2)
    {
      cerr << "ERROR: you must pass at least two parameters to the joint prior statement" << endl;
      exit(EXIT_FAILURE);
    }

  if (prior_shape == PriorDistributions::noShape)
    {
      cerr << "ERROR: You must pass the shape option to the prior statement." << endl;
      exit(EXIT_FAILURE);
    }

  if (options_list.num_options.find("mean") == options_list.num_options.end()
      && options_list.num_options.find("mode") == options_list.num_options.end())
    {
      cerr << "ERROR: You must pass at least one of mean and mode to the prior statement." << endl;
      exit(EXIT_FAILURE);
    }

  if (auto it_num = options_list.num_options.find("domain");
      it_num != options_list.num_options.end())
    {
      using namespace boost;
      vector<string> tokenizedDomain;
      split(tokenizedDomain, it_num->second, is_any_of("[ ]"), token_compress_on);
      if (tokenizedDomain.size() != 4)
        {
          cerr << "ERROR: You must pass exactly two values to the domain option." << endl;
          exit(EXIT_FAILURE);
        }
    }
}

void
JointPriorStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  for (const auto &joint_parameter : joint_parameters)
    output << "eifind = get_new_or_existing_ei_index('joint_parameter_prior_index', '"
           << joint_parameter << "', '');" << endl
           << "estimation_info.joint_parameter_prior_index(eifind) = {'" << joint_parameter << "'};" << endl;

  output << "key = {[";
  for (const auto &joint_parameter : joint_parameters)
    output << "get_new_or_existing_ei_index('joint_parameter_prior_index', '" << joint_parameter << "', '') ..."
           << endl << "    ";
  output << "]};" << endl;

  string lhs_field("estimation_info.joint_parameter_tmp");

  writeOutputHelper(output, "domain", lhs_field);
  writeOutputHelper(output, "interval", lhs_field);
  writeOutputHelper(output, "mean", lhs_field);
  writeOutputHelper(output, "median", lhs_field);
  writeOutputHelper(output, "mode", lhs_field);

  assert(prior_shape != PriorDistributions::noShape);
  output << lhs_field << ".shape = " << static_cast<int>(prior_shape) << ";" << endl;

  writeOutputHelper(output, "shift", lhs_field);
  writeOutputHelper(output, "stdev", lhs_field);
  writeOutputHelper(output, "truncate", lhs_field);
  writeOutputHelper(output, "variance", lhs_field);

  output << "estimation_info.joint_parameter_tmp = [key, ..." << endl
         << "    " << lhs_field << ".domain , ..." << endl
         << "    " << lhs_field << ".interval , ..." << endl
         << "    " << lhs_field << ".mean , ..." << endl
         << "    " << lhs_field << ".median , ..." << endl
         << "    " << lhs_field << ".mode , ..." << endl
         << "    " << lhs_field << ".shape , ..." << endl
         << "    " << lhs_field << ".shift , ..." << endl
         << "    " << lhs_field << ".stdev , ..." << endl
         << "    " << lhs_field << ".truncate , ..." << endl
         << "    " << lhs_field << ".variance];" << endl
         << "estimation_info.joint_parameter = [estimation_info.joint_parameter; estimation_info.joint_parameter_tmp];" << endl
         << "estimation_info=rmfield(estimation_info, 'joint_parameter_tmp');" << endl;
}

void
JointPriorStatement::writeOutputHelper(ostream &output, const string &field, const string &lhs_field) const
{
  output << lhs_field << "." << field << " = {";
  if (field == "variance")
    output << "{";
  if (auto itn = options_list.num_options.find(field);
      itn != options_list.num_options.end())
    output << itn->second;
  else
    output << "{}";
  if (field == "variance")
    output << "}";
  output << "};" << endl;
}

void
JointPriorStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "joint_prior")"
         << R"(, "key": [)";
  for (auto it = joint_parameters.begin(); it != joint_parameters.end(); ++it)
    {
      if (it != joint_parameters.begin())
        output << ", ";
      output << R"(")" << *it << R"(")";
    }
  output << "]";

  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }

  output << R"(, "shape": )";
  switch (prior_shape)
    {
    case PriorDistributions::beta:
      output << R"("beta")";
      break;
    case PriorDistributions::gamma:
      output << R"("gamma")";
      break;
    case PriorDistributions::normal:
      output << R"("normal")";
      break;
    case PriorDistributions::invGamma:
      output << R"("inv_gamma")";
      break;
    case PriorDistributions::uniform:
      output << R"("uniform")";
      break;
    case PriorDistributions::invGamma2:
      output << R"("inv_gamma2")";
      break;
    case PriorDistributions::dirichlet:
      output << R"("dirichlet")";
      break;
    case PriorDistributions::weibull:
      output << R"("weibull")";
      break;
    case PriorDistributions::noShape:
      cerr << "Impossible case." << endl;
      exit(EXIT_FAILURE);
    }
  output << "}";
}

BasicPriorStatement::BasicPriorStatement(string name_arg,
                                         string subsample_name_arg,
                                         PriorDistributions prior_shape_arg,
                                         expr_t variance_arg,
                                         OptionsList options_list_arg) :
  name{move(name_arg)},
  subsample_name{move(subsample_name_arg)},
  prior_shape{prior_shape_arg},
  variance{variance_arg},
  options_list{move(options_list_arg)}
{
}

void
BasicPriorStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  if (prior_shape == PriorDistributions::noShape)
    {
      cerr << "ERROR: You must pass the shape option to the prior statement." << endl;
      exit(EXIT_FAILURE);
    }

  if (options_list.num_options.find("mean") == options_list.num_options.end()
      && options_list.num_options.find("mode") == options_list.num_options.end())
    {
      cerr << "ERROR: You must pass at least one of mean and mode to the prior statement." << endl;
      exit(EXIT_FAILURE);
    }

  if (auto it_stdev = options_list.num_options.find("stdev");
      (it_stdev == options_list.num_options.end() && !variance)
      || (it_stdev != options_list.num_options.end() && variance))
    {
      cerr << "ERROR: You must pass exactly one of stdev and variance to the prior statement." << endl;
      exit(EXIT_FAILURE);
    }

  if (auto it_num = options_list.num_options.find("domain");
      it_num != options_list.num_options.end())
    {
      using namespace boost;
      vector<string> tokenizedDomain;
      split(tokenizedDomain, it_num->second, is_any_of("[ ]"), token_compress_on);
      if (tokenizedDomain.size() != 4)
        {
          cerr << "ERROR: You must pass exactly two values to the domain option." << endl;
          exit(EXIT_FAILURE);
        }
    }
}

bool
BasicPriorStatement::is_structural_innovation(const SymbolType symb_type) const
{
  if (symb_type == SymbolType::exogenous || symb_type == SymbolType::exogenousDet)
    return true;
  return false;
}

void
BasicPriorStatement::get_base_name(const SymbolType symb_type, string &lhs_field) const
{
  if (symb_type == SymbolType::exogenous || symb_type == SymbolType::exogenousDet)
    lhs_field = "structural_innovation";
  else
    lhs_field = "measurement_error";
}

void
BasicPriorStatement::writeCommonOutput(ostream &output, const string &lhs_field) const
{
  output << lhs_field << " = estimation_info.empty_prior;" << endl;

  writeCommonOutputHelper(output, "domain", lhs_field);
  writeCommonOutputHelper(output, "interval", lhs_field);
  writeCommonOutputHelper(output, "mean", lhs_field);
  writeCommonOutputHelper(output, "median", lhs_field);
  writeCommonOutputHelper(output, "mode", lhs_field);

  assert(prior_shape != PriorDistributions::noShape);
  output << lhs_field << ".shape = " << static_cast<int>(prior_shape) << ";" << endl;

  writeCommonOutputHelper(output, "shift", lhs_field);
  writeCommonOutputHelper(output, "stdev", lhs_field);
  writeCommonOutputHelper(output, "truncate", lhs_field);

  if (variance)
    {
      output << lhs_field << ".variance = ";
      variance->writeOutput(output);
      output << ";" << endl;
    }
}

void
BasicPriorStatement::writeCommonOutputHelper(ostream &output, const string &field, const string &lhs_field) const
{
  if (auto itn = options_list.num_options.find(field);
      itn != options_list.num_options.end())
    output << lhs_field << "." << field << " = "<< itn->second << ";" << endl;
}

void
BasicPriorStatement::writePriorOutput(ostream &output, string &lhs_field, const string &name2) const
{
  if (subsample_name.empty())
    lhs_field += ".prior(1)";
  else
    {
      output << "subsamples_indx = get_existing_subsamples_indx('" << name << "','" << name2 << "');" << endl
             << "eisind = get_subsamples_range_indx(subsamples_indx, '" << subsample_name << "');" << endl;
      lhs_field += ".subsample_prior(eisind)";
    }
  writeCommonOutput(output, lhs_field);
}

void
BasicPriorStatement::writeJsonPriorOutput(ostream &output) const
{
  output << R"(, "name": ")" << name << R"(")"
         << R"(, "subsample": ")" << subsample_name << R"(")"
         << ", ";
  writeJsonShape(output);
  if (variance)
    {
      output << R"(, "variance": ")";
      variance->writeJsonOutput(output, {}, {});
      output << R"(")";
    }
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
}

void
BasicPriorStatement::writeJsonShape(ostream &output) const
{
  output << R"("shape": )";
  switch (prior_shape)
    {
    case PriorDistributions::beta:
      output << R"("beta")";
      break;
    case PriorDistributions::gamma:
      output << R"("gamma")";
      break;
    case PriorDistributions::normal:
      output << R"("normal")";
      break;
    case PriorDistributions::invGamma:
      output << R"("inv_gamma")";
      break;
    case PriorDistributions::uniform:
      output << R"("uniform")";
      break;
    case PriorDistributions::invGamma2:
      output << R"("inv_gamma2")";
      break;
    case PriorDistributions::dirichlet:
      output << R"("dirichlet")";
      break;
    case PriorDistributions::weibull:
      output << R"("weibull")";
      break;
    case PriorDistributions::noShape:
      assert(prior_shape != PriorDistributions::noShape);
    }
}

PriorStatement::PriorStatement(string name_arg,
                               string subsample_name_arg,
                               PriorDistributions prior_shape_arg,
                               expr_t variance_arg,
                               OptionsList options_list_arg) :
  BasicPriorStatement{move(name_arg), move(subsample_name_arg), prior_shape_arg, variance_arg, move(options_list_arg)}
{
}

void
PriorStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  string lhs_field = "estimation_info.parameter(eifind)";
  output << "eifind = get_new_or_existing_ei_index('parameter_prior_index', '"
         << name << "', '');" << endl
         << "estimation_info.parameter_prior_index(eifind) = {'" << name << "'};" << endl;
  writePriorOutput(output, lhs_field, "");
}

void
PriorStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "prior")";
  writeJsonPriorOutput(output);
  output << "}";
}

StdPriorStatement::StdPriorStatement(string name_arg,
                                     string subsample_name_arg,
                                     PriorDistributions prior_shape_arg,
                                     expr_t variance_arg,
                                     OptionsList options_list_arg,
                                     const SymbolTable &symbol_table_arg) :
  BasicPriorStatement{move(name_arg), move(subsample_name_arg), prior_shape_arg, variance_arg, move(options_list_arg)},
  symbol_table{symbol_table_arg}
{
}

void
StdPriorStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  string lhs_field;
  get_base_name(symbol_table.getType(name), lhs_field);
  output << "eifind = get_new_or_existing_ei_index('" << lhs_field << "_prior_index', '"
         << name << "', '');" << endl
         << "estimation_info." << lhs_field << "_prior_index(eifind) = {'" << name << "'};" << endl;

  lhs_field = "estimation_info." + lhs_field + "(eifind)";
  writePriorOutput(output, lhs_field, "");
}

void
StdPriorStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "std_prior")";
  writeJsonPriorOutput(output);
  output << "}";
}

CorrPriorStatement::CorrPriorStatement(string name_arg1, string name_arg2,
                                       string subsample_name_arg,
                                       PriorDistributions prior_shape_arg,
                                       expr_t variance_arg,
                                       OptionsList options_list_arg,
                                       const SymbolTable &symbol_table_arg) :
  BasicPriorStatement{move(name_arg1), move(subsample_name_arg), prior_shape_arg, variance_arg, move(options_list_arg)},
  name1{move(name_arg2)},
  symbol_table{symbol_table_arg}
{
}

void
CorrPriorStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  BasicPriorStatement::checkPass(mod_file_struct, warnings);
  if (symbol_table.getType(name) != symbol_table.getType(name1))
    {
      cerr << "ERROR: In the corr(A,B).prior statement, A and B must be of the same type. "
           << "In your case, " << name << " and " << name1 << " are of different "
           << "types." << endl;
      exit(EXIT_FAILURE);
    }
}

void
CorrPriorStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  string lhs_field;
  get_base_name(symbol_table.getType(name), lhs_field);

  output << "eifind = get_new_or_existing_ei_index('" << lhs_field << "_corr_prior_index', '"
         << name << "', '" << name1 << "');" << endl
         << "estimation_info." << lhs_field << "_corr_prior_index(eifind) = {'"
         << name << ":" << name1 << "'};" << endl;

  lhs_field = "estimation_info." + lhs_field + "_corr(eifind)";
  writePriorOutput(output, lhs_field, name1);
}

void
CorrPriorStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "corr_prior")"
         << R"(, "name2": ")" << name1 << R"(")";
  writeJsonPriorOutput(output);
  output << "}";
}

PriorEqualStatement::PriorEqualStatement(string to_declaration_type_arg,
                                         string to_name1_arg,
                                         string to_name2_arg,
                                         string to_subsample_name_arg,
                                         string from_declaration_type_arg,
                                         string from_name1_arg,
                                         string from_name2_arg,
                                         string from_subsample_name_arg,
                                         const SymbolTable &symbol_table_arg) :
  to_declaration_type{move(to_declaration_type_arg)},
  to_name1{move(to_name1_arg)},
  to_name2{move(to_name2_arg)},
  to_subsample_name{move(to_subsample_name_arg)},
  from_declaration_type{move(from_declaration_type_arg)},
  from_name1{move(from_name1_arg)},
  from_name2{move(from_name2_arg)},
  from_subsample_name{move(from_subsample_name_arg)},
  symbol_table{symbol_table_arg}
{
}

void
PriorEqualStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  if ((to_declaration_type != "par" && to_declaration_type != "std" && to_declaration_type != "corr")
      || (from_declaration_type != "par" && from_declaration_type != "std" && from_declaration_type != "corr"))
    {
      cerr << "Internal Dynare Error" << endl;
      exit(EXIT_FAILURE);
    }
}

void
PriorEqualStatement::get_base_name(const SymbolType symb_type, string &lhs_field) const
{
  if (symb_type == SymbolType::exogenous || symb_type == SymbolType::exogenousDet)
    lhs_field = "structural_innovation";
  else
    lhs_field = "measurement_error";
}

void
PriorEqualStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  string lhs_field, rhs_field;

  if (to_declaration_type == "par")
    lhs_field = "parameter";
  else
    get_base_name(symbol_table.getType(to_name1), lhs_field);

  if (from_declaration_type == "par")
    rhs_field = "parameter";
  else
    get_base_name(symbol_table.getType(from_name1), rhs_field);

  if (to_declaration_type == "corr")
    lhs_field += "_corr";

  if (from_declaration_type == "corr")
    rhs_field += "_corr";

  output << "ei_to_ind = get_new_or_existing_ei_index('" << lhs_field << "_prior_index', '"
         << to_name1 << "', '" << to_name2<< "');" << endl
         << "ei_from_ind = get_new_or_existing_ei_index('" << rhs_field << "_prior_index', '"
         << from_name1 << "', '" << from_name2<< "');" << endl
         << "estimation_info." << lhs_field << "_prior_index(ei_to_ind) = {'" << to_name1;

  if (to_declaration_type == "corr")
    output << ":" << to_name2;
  output << "'};" << endl;

  if (to_declaration_type == "par")
    lhs_field = "parameter";

  if (from_declaration_type == "par")
    rhs_field = "parameter";

  lhs_field = "estimation_info." + lhs_field + "(ei_to_ind)";
  rhs_field = "estimation_info." + rhs_field + "(ei_from_ind)";

  if (to_subsample_name.empty())
    lhs_field += ".prior";
  else
    {
      output << "subsamples_to_indx = get_existing_subsamples_indx('" << to_name1 << "','" << to_name2 << "');" << endl
             << "ei_to_ss_ind = get_subsamples_range_indx(subsamples_to_indx, '" << to_subsample_name << "');" << endl;
      lhs_field += ".subsample_prior(ei_to_ss_ind)";
    }

  if (from_subsample_name.empty())
    rhs_field += ".prior";
  else
    {
      output << "subsamples_from_indx = get_existing_subsamples_indx('" << from_name1 << "','" << from_name2 << "');" << endl
             << "ei_from_ss_ind = get_subsamples_range_indx(subsamples_from_indx, '" << from_subsample_name << "');" << endl;
      rhs_field += ".subsample_prior(ei_from_ss_ind)";
    }

  output << lhs_field << " = " << rhs_field << ";" << endl;
}

void
PriorEqualStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "prior_equal")"
         << R"(, "to_name1": ")" << to_name1 << R"(")";
  if (to_declaration_type == "corr")
    output << R"(, "to_name2": ")" << to_name2 << R"(")";
  output << R"(, "to_subsample": ")" << to_subsample_name << R"(")"
         << R"(, "from_name1": ")" << from_name1 << R"(")";
  if (to_declaration_type == "corr")
    output << R"(, "from_name2": ")" << from_name2 << R"(")";
  output << R"(, "from_subsample": ")" << from_subsample_name << R"(")"
         << "}";
}

BasicOptionsStatement::BasicOptionsStatement(string name_arg,
                                             string subsample_name_arg,
                                             OptionsList options_list_arg) :
  name{move(name_arg)},
  subsample_name{move(subsample_name_arg)},
  options_list{move(options_list_arg)}
{
}

void
BasicOptionsStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
}

bool
BasicOptionsStatement::is_structural_innovation(const SymbolType symb_type) const
{
  return symb_type == SymbolType::exogenous || symb_type == SymbolType::exogenousDet;
}

void
BasicOptionsStatement::get_base_name(const SymbolType symb_type, string &lhs_field) const
{
  if (symb_type == SymbolType::exogenous || symb_type == SymbolType::exogenousDet)
    lhs_field = "structural_innovation";
  else
    lhs_field = "measurement_error";
}

void
BasicOptionsStatement::writeCommonOutput(ostream &output, const string &lhs_field) const
{
  output << lhs_field << " = estimation_info.empty_options;" << endl;

  writeCommonOutputHelper(output, "bounds", lhs_field);
  writeCommonOutputHelper(output, "init", lhs_field);
  writeCommonOutputHelper(output, "jscale", lhs_field);
}

void
BasicOptionsStatement::writeCommonOutputHelper(ostream &output, const string &field, const string &lhs_field) const
{
  if (auto itn = options_list.num_options.find(field);
      itn != options_list.num_options.end())
    output << lhs_field << "." << field << " = " << itn->second << ";" << endl;
}

void
BasicOptionsStatement::writeOptionsOutput(ostream &output, string &lhs_field, const string &name2) const
{
  if (subsample_name.empty())
    lhs_field += ".options(1)";
  else
    {
      output << "subsamples_indx = get_existing_subsamples_indx('" << name << "','" << name2 << "');" << endl
             << "eisind = get_subsamples_range_indx(subsamples_indx, '" << subsample_name << "');" << endl;
      lhs_field += ".subsample_options(eisind)";
    }
  writeCommonOutput(output, lhs_field);
}

void
BasicOptionsStatement::writeJsonOptionsOutput(ostream &output) const
{
  output << R"(, "name": ")" << name << R"(")";
  if (!subsample_name.empty())
    output << R"(, "subsample_name": ")" << subsample_name << R"(")";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
}

OptionsStatement::OptionsStatement(string name_arg,
                                   string subsample_name_arg,
                                   OptionsList options_list_arg) :
  BasicOptionsStatement{move(name_arg), move(subsample_name_arg), move(options_list_arg)}
{
}

void
OptionsStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  string lhs_field = "estimation_info.parameter(eifind)";
  output << "eifind = get_new_or_existing_ei_index('parameter_options_index', '"
         << name << "', '');" << endl
         << "estimation_info.parameter_options_index(eifind) = {'" << name << "'};" << endl;
  writeOptionsOutput(output, lhs_field, "");
}

void
OptionsStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "options")";
  writeJsonOptionsOutput(output);
  output << "}";
}

StdOptionsStatement::StdOptionsStatement(string name_arg,
                                         string subsample_name_arg,
                                         OptionsList options_list_arg,
                                         const SymbolTable &symbol_table_arg) :
  BasicOptionsStatement{move(name_arg), move(subsample_name_arg), move(options_list_arg)},
  symbol_table{symbol_table_arg}
{
}

void
StdOptionsStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  string lhs_field;
  get_base_name(symbol_table.getType(name), lhs_field);
  output << "eifind = get_new_or_existing_ei_index('" << lhs_field << "_options_index', '"
         << name << "', '');" << endl
         << "estimation_info." << lhs_field << "_options_index(eifind) = {'" << name << "'};" << endl;

  lhs_field = "estimation_info." + lhs_field + "(eifind)";
  writeOptionsOutput(output, lhs_field, "");
}

void
StdOptionsStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "std_options")";
  writeJsonOptionsOutput(output);
  output << "}";
}

CorrOptionsStatement::CorrOptionsStatement(string name_arg1, string name_arg2,
                                           string subsample_name_arg,
                                           OptionsList options_list_arg,
                                           const SymbolTable &symbol_table_arg) :
  BasicOptionsStatement{move(name_arg1), move(subsample_name_arg), move(options_list_arg)},
  name1{move(name_arg2)},
  symbol_table{symbol_table_arg}
{
}

void
CorrOptionsStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  if (symbol_table.getType(name) != symbol_table.getType(name1))
    {
      cerr << "ERROR: In the corr(A,B).options statement, A and B must be of the same type. "
           << "In your case, " << name << " and " << name1 << " are of different "
           << "types." << endl;
      exit(EXIT_FAILURE);
    }
}

void
CorrOptionsStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  string lhs_field;
  get_base_name(symbol_table.getType(name), lhs_field);

  output << "eifind = get_new_or_existing_ei_index('" << lhs_field << "_corr_options_index', '"
         << name << "', '" << name1 << "');" << endl
         << "estimation_info." << lhs_field << "_corr_options_index(eifind) = {'"
         << name << ":" << name1 << "'};" << endl;

  lhs_field = "estimation_info." + lhs_field + "_corr(eifind)";
  writeOptionsOutput(output, lhs_field, name1);
}

void
CorrOptionsStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "corr_options")"
         << R"(, "name2": ")" << name1 << R"(")";
  writeJsonOptionsOutput(output);
  output << "}";
}

OptionsEqualStatement::OptionsEqualStatement(string to_declaration_type_arg,
                                             string to_name1_arg,
                                             string to_name2_arg,
                                             string to_subsample_name_arg,
                                             string from_declaration_type_arg,
                                             string from_name1_arg,
                                             string from_name2_arg,
                                             string from_subsample_name_arg,
                                             const SymbolTable &symbol_table_arg) :
  to_declaration_type{move(to_declaration_type_arg)},
  to_name1{move(to_name1_arg)},
  to_name2{move(to_name2_arg)},
  to_subsample_name{move(to_subsample_name_arg)},
  from_declaration_type{move(from_declaration_type_arg)},
  from_name1{move(from_name1_arg)},
  from_name2{move(from_name2_arg)},
  from_subsample_name{move(from_subsample_name_arg)},
  symbol_table{symbol_table_arg}
{
}

void
OptionsEqualStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  if ((to_declaration_type != "par" && to_declaration_type != "std" && to_declaration_type != "corr")
      || (from_declaration_type != "par" && from_declaration_type != "std" && from_declaration_type != "corr"))
    {
      cerr << "Internal Dynare Error" << endl;
      exit(EXIT_FAILURE);
    }
}

void
OptionsEqualStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "options_equal")"
         << R"(, "to_name1": ")" << to_name1 << R"(")";
  if (to_declaration_type == "corr")
    output << R"(, "to_name2": ")" << to_name2 << R"(")";
  output << R"(, "to_subsample": ")" << to_subsample_name << R"(")"
         << R"(, "from_name1": ")" << from_name1 << R"(")";
  if (to_declaration_type == "corr")
    output << R"(, "from_name2": ")" << from_name2 << R"(")";
  output << R"(, "from_subsample": ")" << from_subsample_name << R"(")"
         << "}";
}

void
OptionsEqualStatement::get_base_name(const SymbolType symb_type, string &lhs_field) const
{
  if (symb_type == SymbolType::exogenous || symb_type == SymbolType::exogenousDet)
    lhs_field = "structural_innovation";
  else
    lhs_field = "measurement_error";
}

void
OptionsEqualStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  string lhs_field, rhs_field;

  if (to_declaration_type == "par")
    lhs_field = "parameter";
  else
    get_base_name(symbol_table.getType(to_name1), lhs_field);

  if (from_declaration_type == "par")
    rhs_field = "parameter";
  else
    get_base_name(symbol_table.getType(from_name1), rhs_field);

  if (to_declaration_type == "corr")
    lhs_field += "_corr";

  if (from_declaration_type == "corr")
    rhs_field += "_corr";

  output << "ei_to_ind = get_new_or_existing_ei_index('" << lhs_field << "_options_index', '"
         << to_name1 << "', '" << to_name2<< "');" << endl
         << "ei_from_ind = get_new_or_existing_ei_index('" << rhs_field << "_options_index', '"
         << from_name1 << "', '" << from_name2<< "');" << endl
         << "estimation_info." << lhs_field << "_options_index(ei_to_ind) = {'" << to_name1;

  if (to_declaration_type == "corr")
    output << ":" << to_name2;
  output << "'};" << endl;

  if (to_declaration_type == "par")
    lhs_field = "parameter";

  if (from_declaration_type == "par")
    rhs_field = "parameter";

  lhs_field = "estimation_info." + lhs_field + "(ei_to_ind)";
  rhs_field = "estimation_info." + rhs_field + "(ei_from_ind)";

  if (to_subsample_name.empty())
    lhs_field += ".options";
  else
    {
      output << "subsamples_to_indx = get_existing_subsamples_indx('" << to_name1 << "','" << to_name2 << "');" << endl
             << "ei_to_ss_ind = get_subsamples_range_indx(subsamples_to_indx, '" << to_subsample_name << "');" << endl;
      lhs_field += ".subsample_options(ei_to_ss_ind)";
    }

  if (from_subsample_name.empty())
    rhs_field += ".options";
  else
    {
      output << "subsamples_from_indx = get_existing_subsamples_indx('" << from_name1 << "','" << from_name2 << "');" << endl
             << "ei_from_ss_ind = get_subsamples_range_indx(subsamples_from_indx, '" << from_subsample_name << "');" << endl;
      rhs_field += ".subsample_options(ei_from_ss_ind)";
    }

  output << lhs_field << " = " << rhs_field << ";" << endl;
}

CalibSmootherStatement::CalibSmootherStatement(SymbolList symbol_list_arg,
                                               OptionsList options_list_arg)
  : symbol_list{move(symbol_list_arg)}, options_list{move(options_list_arg)}
{
}

void
CalibSmootherStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.calib_smoother_present = true;
  try
    {
      const vector<SymbolType> valid_symbol_list_types { SymbolType::endogenous };
      symbol_list.checkPass(warnings, valid_symbol_list_types);
    }
  catch (SymbolList::SymbolListException &e)
    {
      cerr << "ERROR: calib_smoother: " << e.message << endl;
      exit(EXIT_FAILURE);
    }
}

void
CalibSmootherStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  options_list.writeOutput(output);
  if (options_list.string_options.find("parameter_set") == options_list.string_options.end())
    output << "options_.parameter_set = 'calibration';" << endl;
  symbol_list.writeOutput("var_list_", output);
  output << "options_.smoother = true;" << endl
         << "options_.order = 1;" << endl
         << "[oo_, M_, options_, bayestopt_] = evaluate_smoother(options_.parameter_set, var_list_, M_, oo_, options_, bayestopt_, estim_params_);" << endl;
}

void
CalibSmootherStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "calib_smoother")";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  if (!symbol_list.empty())
    {
      output << ", ";
      symbol_list.writeJsonOutput(output);
    }
  output << "}";
}

ExtendedPathStatement::ExtendedPathStatement(OptionsList options_list_arg)
  : options_list{move(options_list_arg)}
{
}

void
ExtendedPathStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.extended_path_present = true;

  if (options_list.num_options.find("periods") == options_list.num_options.end())
    {
      cerr << "ERROR: the 'periods' option of 'extended_path' is mandatory" << endl;
      exit(EXIT_FAILURE);
    }

  // Fill in option_occbin of mod_file_struct
  if (options_list.num_options.find("occbin") != options_list.string_options.end())
    mod_file_struct.occbin_option = true;
}

void
ExtendedPathStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  // Beware: options do not have the same name in the interface and in the M code...

  for (const auto &num_option : options_list.num_options)
    if (num_option.first != "periods")
      output << "options_." << num_option.first << " = " << num_option.second << ";" << endl;

  output << "extended_path([], " << options_list.num_options.find("periods")->second
         << ", [], options_, M_, oo_);" << endl;
}

void
ExtendedPathStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "extended_path")";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

void
ModelDiagnosticsStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "model_diagnostics(M_,options_,oo_);" << endl;
}

void
ModelDiagnosticsStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "model_diagnostics"})";
}

Smoother2histvalStatement::Smoother2histvalStatement(OptionsList options_list_arg) :
  options_list{move(options_list_arg)}
{
}

void
Smoother2histvalStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  options_list.writeOutput(output, "options_smoother2histval");
  output << "smoother2histval(options_smoother2histval);" << endl;
}

void
Smoother2histvalStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "smoother_2_histval")";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

GMMEstimationStatement::GMMEstimationStatement(SymbolList symbol_list_arg,
                                               OptionsList options_list_arg) :
  symbol_list{move(symbol_list_arg)},
  options_list{move(options_list_arg)}
{
}

void
GMMEstimationStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  try
    {
      const vector<SymbolType> valid_symbol_list_types { SymbolType::endogenous };
      symbol_list.checkPass(warnings, valid_symbol_list_types);
    }
  catch (SymbolList::SymbolListException &e)
    {
      cerr << "ERROR: gmm_estimation: " << e.message << endl;
      exit(EXIT_FAILURE);
    }
}

void
GMMEstimationStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  symbol_list.writeOutput("var_list_", output);
  options_list.writeOutput(output);
  output << "[M_, oo_, estim_params_, bayestopt_, dataset_, dataset_info] = "
         << "GMM_SMM_estimation_core(var_list_, M_, options_, oo_, estim_params_, bayestopt_, dataset_, dataset_info, 'GMM');" << endl;
}

void
GMMEstimationStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "gmm_estimation")";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  if (!symbol_list.empty())
    {
      output << ", ";
      symbol_list.writeJsonOutput(output);
    }
  output << "}";
}

SMMEstimationStatement::SMMEstimationStatement(SymbolList symbol_list_arg,
                                               OptionsList options_list_arg) :
  symbol_list{move(symbol_list_arg)},
  options_list{move(options_list_arg)}
{
}

void
SMMEstimationStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  try
    {
      const vector<SymbolType> valid_symbol_list_types { SymbolType::endogenous };
      symbol_list.checkPass(warnings, valid_symbol_list_types);
    }
  catch (SymbolList::SymbolListException &e)
    {
      cerr << "ERROR: smm_estimation: " << e.message << endl;
      exit(EXIT_FAILURE);
    }
}

void
SMMEstimationStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  symbol_list.writeOutput("var_list_", output);
  options_list.writeOutput(output);
  output << "[M_, oo_, estim_params_, bayestopt_, dataset_, dataset_info] = "
         << "GMM_SMM_estimation_core(var_list_, M_, options_, oo_, estim_params_, bayestopt_, dataset_, dataset_info, 'SMM');" << endl;
}

void
SMMEstimationStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "smm_estimation")";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  if (!symbol_list.empty())
    {
      output << ", ";
      symbol_list.writeJsonOutput(output);
    }
  output << "}";
}

GenerateIRFsStatement::GenerateIRFsStatement(OptionsList options_list_arg,
                                             vector<string> generate_irf_names_arg,
                                             vector<map<string, double>> generate_irf_elements_arg) :
  options_list{move(options_list_arg)},
  generate_irf_names{move(generate_irf_names_arg)},
  generate_irf_elements{move(generate_irf_elements_arg)}
{
}

void
GenerateIRFsStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  options_list.writeOutput(output);

  if (generate_irf_names.empty())
    return;

  output << "options_.irf_opt.irf_shock_graphtitles = { ";
  for (const auto &generate_irf_name : generate_irf_names)
    output << "'" << generate_irf_name << "'; ";
  output << "};" << endl;

  output << "options_.irf_opt.irf_shocks = zeros(M_.exo_nbr, "
         << generate_irf_names.size() << ");" << endl;

  for (size_t i = 0; i < generate_irf_names.size(); i++)
    {
      map<string, double> m = generate_irf_elements[i];
      for (auto &it : m)
        output << "options_.irf_opt.irf_shocks(M_.exo_names == '"
               << it.first << "', " << i + 1 << ") = "
               << it.second << ";" << endl;
    }
}

void
GenerateIRFsStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "generate_irfs")";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }

  if (!generate_irf_names.empty())
    {
      output << R"(, "irf_elements": [)";
      for (size_t i = 0; i < generate_irf_names.size(); i++)
        {
          output << R"({"name": ")" << generate_irf_names[i] << R"(", "shocks": [)";
          map<string, double> m = generate_irf_elements[i];
          size_t idx = 0;
          for (auto it = m.begin(); it != m.end(); ++it, idx++)
            {
              output << R"({"exogenous_variable": ")" << it->first << R"(", )"
                     << R"("exogenous_variable_value": ")" << it->second << R"("})";
              if (idx + 1 < m.size())
                output << ", ";
            }
          output << "]}";
          if (i + 1 < generate_irf_names.size())
            output << ", ";
        }
      output << "]";
    }
  output << "}";
}

VarExpectationModelStatement::VarExpectationModelStatement(string model_name_arg,
                                                           expr_t expression_arg,
                                                           string aux_model_name_arg,
                                                           string horizon_arg,
                                                           expr_t discount_arg,
                                                           const SymbolTable &symbol_table_arg) :
  model_name{move(model_name_arg)}, expression{expression_arg},
  aux_model_name{move(aux_model_name_arg)}, horizon{move(horizon_arg)},
  discount{discount_arg}, symbol_table{symbol_table_arg}
{
}

void
VarExpectationModelStatement::substituteUnaryOpNodes(const lag_equivalence_table_t &nodes, ExprNode::subst_table_t &subst_table)
{
  vector<BinaryOpNode *> neweqs;
  expression = expression->substituteUnaryOpNodes(nodes, subst_table, neweqs);
  if (neweqs.size() > 0)
    {
      cerr << "ERROR: the 'expression' option of var_expectation_model contains a variable with a unary operator that is not present in the VAR model" << endl;
      exit(EXIT_FAILURE);
    }
}

void
VarExpectationModelStatement::substituteDiff(const lag_equivalence_table_t &nodes, ExprNode::subst_table_t &subst_table)
{
  vector<BinaryOpNode *> neweqs;
  expression = expression->substituteDiff(nodes, subst_table, neweqs);
  if (neweqs.size() > 0)
    {
      cerr << "ERROR: the 'expression' option of var_expectation_model contains a diff'd variable that is not present in the VAR model" << endl;
      exit(EXIT_FAILURE);
    }
}

void
VarExpectationModelStatement::matchExpression()
{
  try
    {
      auto vpc = expression->matchLinearCombinationOfVariables();
      for (const auto &it : vpc)
        {
          if (get<1>(it) != 0)
            throw ExprNode::MatchFailureException{"lead/lags are not allowed"};
          if (symbol_table.getType(get<0>(it)) != SymbolType::endogenous)
            throw ExprNode::MatchFailureException{"Variable is not an endogenous"};
          vars_params_constants.emplace_back(get<0>(it), get<2>(it), get<3>(it));
        }
    }
  catch (ExprNode::MatchFailureException &e)
    {
      cerr << "ERROR: expression in var_expectation_model is not of the expected form: " << e.message << endl;
      exit(EXIT_FAILURE);
    }
}

void
VarExpectationModelStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  string mstruct = "M_.var_expectation." + model_name;
  output << mstruct << ".auxiliary_model_name = '" << aux_model_name << "';" << endl
         << mstruct << ".horizon = " << horizon << ';' << endl;

  if (!vars_params_constants.size())
    {
      cerr << "ERROR: VarExpectationModelStatement::writeOutput: matchExpression() has not been called" << endl;
      exit(EXIT_FAILURE);
    }

  ostringstream vars_list, params_list, constants_list;
  for (auto it = vars_params_constants.begin(); it != vars_params_constants.end(); ++it)
    {
      if (it != vars_params_constants.begin())
        {
          vars_list << ", ";
          params_list << ", ";
          constants_list << ", ";
        }
      vars_list << symbol_table.getTypeSpecificID(get<0>(*it))+1;
      if (get<1>(*it) == -1)
        params_list << "NaN";
      else
        params_list << symbol_table.getTypeSpecificID(get<1>(*it))+1;
      constants_list << get<2>(*it);
    }
  output << mstruct << ".expr.vars = [ " << vars_list.str() << " ];" << endl
         << mstruct << ".expr.params = [ " << params_list.str() << " ];" << endl
         << mstruct << ".expr.constants = [ " << constants_list.str() << " ];" << endl;

  if (auto disc_var = dynamic_cast<const VariableNode *>(discount);
      disc_var)
    output << mstruct << ".discount_index = " << symbol_table.getTypeSpecificID(disc_var->symb_id) + 1 << ';' << endl;
  else
    {
      output << mstruct << ".discount_value = ";
      discount->writeOutput(output);
      output << ';' << endl;
    }
  output << mstruct << ".param_indices = [ ";
  for (int param_id : aux_params_ids)
    output << symbol_table.getTypeSpecificID(param_id)+1 << ' ';
  output << "];" << endl;
}

void
VarExpectationModelStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "var_expectation_model",)"
         << R"("model_name": ")" << model_name << R"(", )"
         << R"("expression": ")";
  expression->writeOutput(output);
  output << R"(", )"
         << R"("auxiliary_model_name": ")" << aux_model_name << R"(", )"
         << R"("horizon": ")" << horizon << R"(", )"
         << R"("discount": ")";
  discount->writeOutput(output);
  output << R"("})";
}
