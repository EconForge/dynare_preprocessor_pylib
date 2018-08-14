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

#include <cassert>
#include <iostream>
#include <sstream>

using namespace std;

#include "ComputingTasks.hh"
#include "Statement.hh"

#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/tokenizer.hpp>
#include <utility>

SteadyStatement::SteadyStatement(OptionsList options_list_arg) :
  options_list(move(options_list_arg))
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
  output << "{\"statementName\": \"steady\"";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

CheckStatement::CheckStatement(OptionsList options_list_arg) :
  options_list(move(options_list_arg))
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
  output << "{\"statementName\": \"check\"";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

ModelInfoStatement::ModelInfoStatement(OptionsList options_list_arg) :
  options_list(move(options_list_arg))
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
  output << "{\"statementName\": \"model_info\"";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

SimulStatement::SimulStatement(OptionsList options_list_arg) :
  options_list(move(options_list_arg))
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
  options_list.writeOutput(output);
  output << "perfect_foresight_setup;" << endl
         << "perfect_foresight_solver;" << endl;
}

void
SimulStatement::writeJsonOutput(ostream &output) const
{
  output << "{\"statementName\": \"simul\"";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

PerfectForesightSetupStatement::PerfectForesightSetupStatement(OptionsList options_list_arg) :
  options_list(move(options_list_arg))
{
}

void
PerfectForesightSetupStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  options_list.writeOutput(output);
  output << "perfect_foresight_setup;" << endl;
}

void
PerfectForesightSetupStatement::writeJsonOutput(ostream &output) const
{
  output << "{\"statementName\": \"perfect_foresight_setup\"";
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
  output << "{\"statementName\": \"perfect_foresight_solver\"";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

PriorPosteriorFunctionStatement::PriorPosteriorFunctionStatement(const bool prior_func_arg,
                                                                 OptionsList options_list_arg) :
  prior_func(prior_func_arg),
  options_list(move(options_list_arg))
{
}

void
PriorPosteriorFunctionStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  auto it2 = options_list.string_options.find("function");
  if (it2 == options_list.string_options.end() || it2->second.empty())
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
  string type = "posterior";
  if (prior_func)
    type = "prior";

  output << "oo_ = execute_prior_posterior_function("
         << "'" << options_list.string_options.find("function")->second << "', "
         << "M_, options_, oo_, estim_params_, bayestopt_, dataset_, dataset_info, "
         << "'" << type << "');" << endl;
}

void
PriorPosteriorFunctionStatement::writeJsonOutput(ostream &output) const
{
  string type = "posterior";
  if (prior_func)
    type = "prior";
  output << "{\"statementName\": \"prior_posterior_function\", \"type\": \"" << type << "\"";
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
                                     string growth_arg,
                                     const SymbolTable &symbol_table_arg) :
  name{move(name_arg)},
  aux_model_name{move(aux_model_name_arg)},
  discount{move(discount_arg)},
  growth{move(growth_arg)},
  symbol_table{symbol_table_arg}
{
}

void
PacModelStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.pac_params.insert(symbol_table.getID(discount));
  if (!growth.empty())
    mod_file_struct.pac_params.insert(symbol_table.getID(growth));
}

void
PacModelStatement::fillUndiffedLHS(vector<int> &lhs_arg)
{
  lhs = lhs_arg;
}

void
PacModelStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "M_.pac." << name << ".auxiliary_model_name = '" << aux_model_name << "';" << endl
         << "M_.pac." << name << ".discount_index = " << symbol_table.getTypeSpecificID(discount) + 1 << ";" << endl;

  if (!growth.empty())
    {
      output << "M_.pac." << name << ".growth_index = " << symbol_table.getTypeSpecificID(growth) + 1 << ";" << endl
             << "M_.pac." << name << ".growth_type = ";
      switch(symbol_table.getType(growth))
        {
        case SymbolType::endogenous:
          output << "'endogenous';" << endl;
          break;
        case SymbolType::exogenous:
          output << "'exogenous';" << endl;
          break;
        case SymbolType::parameter:
          output << "'parameter';" << endl;
          break;
        default:
          cerr << "pac_model: error encountered in growth type" << endl;
          exit(EXIT_FAILURE);
        }
    }

  output << "M_.pac." << name << ".lhs = [";
  for (auto it = lhs.begin(); it !=lhs.end(); it++)
    {
      if (it != lhs.begin())
        output << " ";
      output << *it + 1;
    }
  output << "];" << endl;
  /*
         << "M_.pac." << name << ".undiff_eqtags = {";
  for (auto it = undiff.begin(); it != undiff.end(); it++)
    {
      if (it != undiff.begin())
        output << "; ";
      output << "'" << it->first << "'";
    }
  output << "};" << endl
         << "M_.pac." << name << ".undiff_num = [";
  for (auto it = undiff.begin(); it != undiff.end(); it++)
    {
      if (it != undiff.begin())
        output << " ";
      output << it->second;
    }
  output << "];" << endl;
  */
}

void
PacModelStatement::writeJsonOutput(ostream &output) const
{
  output << "{\"statementName\": \"pac_model\","
         << "\"model_name\": \"" << name << "\","
         << "\"auxiliary_model_name\": \"" << aux_model_name << "\","
         << "\"discount_index\": " << symbol_table.getTypeSpecificID(discount) + 1;

  if (!growth.empty())
    {
      output << ","
             << "\"growth_index\": " << symbol_table.getTypeSpecificID(growth) + 1 << ","
             << "\"growth_type\": ";
      switch(symbol_table.getType(growth))
        {
        case SymbolType::endogenous:
          output << "\"endogenous\"" << endl;
          break;
        case SymbolType::exogenous:
          output << "\"exogenous\"" << endl;
          break;
        case SymbolType::parameter:
          output << "\"parameter\"" << endl;
          break;
        default:
          cerr << "pac_model: error encountered in growth type" << endl;
          exit(EXIT_FAILURE);
        }
    }
  output << "}";
}

tuple<string, string, int>
PacModelStatement::getPacModelInfoForPacExpectation() const
{
  int growth_symb_id = -1;
  if (!growth.empty())
    growth_symb_id = symbol_table.getID(growth);
  return { name, aux_model_name, growth_symb_id };
}

VarModelStatement::VarModelStatement(SymbolList symbol_list_arg,
                                     OptionsList options_list_arg,
                                     string name_arg,
                                     const SymbolTable &symbol_table_arg) :
  symbol_list(move(symbol_list_arg)),
  options_list(move(options_list_arg)),
  name{move(name_arg)},
  symbol_table(symbol_table_arg)
{
}

void
VarModelStatement::getVarModelInfo(string &var_model_name,
                                   map<string, pair<SymbolList, int>> &var_model_info,
                                   map<string, vector<string>> &var_model_eqtags) const
{
  var_model_name = name;
  if (symbol_list.empty())
    {
      auto it = options_list.vector_str_options.find("var.eqtags");
      var_model_eqtags[name] = it->second;
    }
  else
    {
      auto it = options_list.num_options.find("var.order");
      var_model_info[name] = { symbol_list, stoi(it->second) };
    }
}

void
VarModelStatement::fillVarModelInfoFromEquations(vector<int> &eqnumber_arg, vector<int> &lhs_arg,
                                                 vector<set<pair<int, int>>> &rhs_arg,
                                                 vector<bool> &nonstationary_arg,
                                                 vector<bool> &diff_arg,
                                                 vector<int> &orig_diff_var_arg,
                                                 int max_lag_arg)
{
  eqnumber = eqnumber_arg;
  lhs = lhs_arg;
  rhs_by_eq = rhs_arg;
  nonstationary = nonstationary_arg;
  diff = diff_arg;
  orig_diff_var = orig_diff_var_arg;
  max_lag = max_lag_arg;
}

void
VarModelStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
}

void
VarModelStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  options_list.writeOutput(output);
  if (!symbol_list.empty())
    symbol_list.writeOutput("options_.var.var_list_", output);

  output << "M_.var." << name << " = options_.var;" << endl
         << "options_ = rmfield(options_, 'var');" << endl;

  output << "M_.var." << name << ".eqn = [";
  for (auto it = eqnumber.begin();
       it != eqnumber.end(); it++)
    {
      if (it != eqnumber.begin())
        output << " ";
      output << *it + 1;
    }
  output << "];" << endl
         << "M_.var." << name << ".lhs = [";
  for (auto it = lhs.begin();
       it != lhs.end(); it++)
    {
      if (it != lhs.begin())
        output << " ";
      output << symbol_table.getTypeSpecificID(*it) + 1;
    }
  output << "];" << endl
         << "M_.var." << name << ".max_lag = " << max_lag << ";" << endl
         << "M_.var." << name << ".nonstationary = [";
  for (auto it = nonstationary.begin();
       it != nonstationary.end(); it++)
    {
      if (it != nonstationary.begin())
        output << " ";
      if (*it)
        output << "true";
      else
        output << "false";
    }
  output << "];" << endl
         << "M_.var." << name << ".diff = [";
  for (auto it = diff.begin();
       it != diff.end(); it++)
    {
      if (it != diff.begin())
        output << " ";
      if (*it)
        output << "true";
      else
        output << "false";
    }
  output << "];" << endl
         << "M_.var." << name << ".orig_diff_var = [";
  for (auto it = orig_diff_var.begin();
       it != orig_diff_var.end(); it++)
    {
      if (it != orig_diff_var.begin())
        output << " ";
      if (*it == -1)
        output << -1;
      else
        output << symbol_table.getTypeSpecificID(*it) + 1;
    }
  output << "];" << endl;
  int i = 1;
  for (auto it = rhs_by_eq.begin();
       it != rhs_by_eq.end(); it++, i++)
    {
      output << "M_.var." << name << ".rhs.vars_at_eq{" << i << "}.var = [";
      for (auto it1 = it->begin();
           it1 != it->end(); it1++)
        {
          if (it1 != it->begin())
            output << " ";
          output << symbol_table.getTypeSpecificID(it1->first) + 1;
        }
      output << "];" << endl
             << "M_.var." << name << ".rhs.vars_at_eq{" << i << "}.lag = [";
      for (auto it1 = it->begin();
           it1 != it->end(); it1++)
        {
          if (it1 != it->begin())
            output << " ";
          output << it1->second;
        }
      output << "];" << endl;
    }
}

void
VarModelStatement::writeJsonOutput(ostream &output) const
{
  output << "{\"statementName\": \"var_model\","
         << "\"model_name\": \"" << name << "\"}";
}

VarEstimationStatement::VarEstimationStatement(OptionsList options_list_arg) :
  options_list(move(options_list_arg))
{
}

void
VarEstimationStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  auto it = options_list.string_options.find("var_estimation.model_name");
  if (it == options_list.string_options.end())
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
  exclusion_restrictions(move(exclusion_restrictions_arg)),
  equation_restrictions(move(equation_restrictions_arg)),
  crossequation_restrictions(move(crossequation_restrictions_arg)),
  covariance_number_restriction(move(covariance_number_restriction_arg)),
  covariance_pair_restriction(move(covariance_pair_restriction_arg)),
  symbol_table(symbol_table_arg)
{
}

int
VarRestrictionsStatement::findIdxInVector(const vector<string> &vecvars, const string &var) const
{
  int idx = 0;
  bool setflag = false;
  for (auto itvs = vecvars.begin();
       itvs != vecvars.end(); itvs++, idx++)
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
       it != exclusion_restrictions.end(); it++, idx++)
    {
      output << Mstr << "exclusion_restrictions{" << idx<< "}.lag = "
             << it->first << ";" << endl
             << Mstr << "exclusion_restrictions{" << idx << "}.restrictions = [";
      for (auto it1 = it->second.begin();
           it1 != it->second.end(); it1++)
        {
          if (it1 != it->second.begin())
            output << " ";

          output << "struct('eq', " << findIdxInVector(vars, symbol_table.getName(it1->first)) + 1
                 << ", 'vars', [";
          vector<string> excvars = it1->second.getSymbols();
          for (vector<string>::const_iterator itvs1 = excvars.begin();
               itvs1 != excvars.end(); itvs1++)
            output << findIdxInVector(vars, *itvs1) + 1 << " ";
          output << "])";
          nrestrictions += it1->second.getSize();
        }
      output << "];" << endl;
    }

  // Equation Restrictions
  idx = 1;
  for (auto it = equation_restrictions.begin();
       it != equation_restrictions.end(); it++, idx++, nrestrictions++)
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
       it != crossequation_restrictions.end(); it++, idx++, nrestrictions++)
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
       it != covariance_number_restriction.end(); it++, idx++)
    output << Mstr << "covariance_const_restriction{" << idx << "}.var1 = '"
           << symbol_table.getName(it->first.first) << "';" << endl
           << Mstr << "covariance_const_restriction{" << idx << "}.var2 = '"
           << symbol_table.getName(it->first.second) << "';" << endl
           << Mstr << "covariance_const_restriction{" << idx << "}.val = "
           << it->second << ";" << endl;

  // Covariance Pair Restrictions
  idx = 1;
  for (auto it = covariance_pair_restriction.begin();
       it != covariance_pair_restriction.end(); it++, idx++)
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
  symbol_list(move(symbol_list_arg)),
  options_list(move(options_list_arg))
{
}

void
StochSimulStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.stoch_simul_present = true;

  // Fill in option_order of mod_file_struct
  auto it = options_list.num_options.find("order");
  if (it != options_list.num_options.end())
    mod_file_struct.order_option = max(mod_file_struct.order_option, stoi(it->second));

  // Fill in mod_file_struct.partial_information
  it = options_list.num_options.find("partial_information");
  if (it != options_list.num_options.end() && it->second == "1")
    mod_file_struct.partial_information = true;

  // Option k_order_solver (implicit when order >= 3)
  it = options_list.num_options.find("k_order_solver");
  if ((it != options_list.num_options.end() && it->second == "1")
      || mod_file_struct.order_option >= 3)
    mod_file_struct.k_order_solver = true;

  it = options_list.num_options.find("hp_filter");
  auto it1 = options_list.num_options.find("bandpass.indicator");
  auto it2 = options_list.num_options.find("one_sided_hp_filter");
  if ((it != options_list.num_options.end() && it1 != options_list.num_options.end())
      || (it != options_list.num_options.end() && it2 != options_list.num_options.end())
      || (it1 != options_list.num_options.end() && it2 != options_list.num_options.end()))
    {
      cerr << "ERROR: stoch_simul: can only use one of hp, one-sided hp, and bandpass filters"
           << endl;
      exit(EXIT_FAILURE);
    }
}

void
StochSimulStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  // Ensure that order 3 implies k_order (#844)
  auto it = options_list.num_options.find("order");
  auto it1 = options_list.num_options.find("k_order_solver");
  if ((it1 != options_list.num_options.end() && it1->second == "1")
      || (it != options_list.num_options.end() && stoi(it->second) >= 3))
    output << "options_.k_order_solver = 1;" << endl;

  options_list.writeOutput(output);
  symbol_list.writeOutput("var_list_", output);
  output << "info = stoch_simul(var_list_);" << endl;
}

void
StochSimulStatement::writeJsonOutput(ostream &output) const
{
  output << "{\"statementName\": \"stoch_simul\"";
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
  symbol_list(move(symbol_list_arg)),
  options_list(move(options_list_arg))
{
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
  output << "{\"statementName\": \"forecast\"";
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

RamseyModelStatement::RamseyModelStatement(OptionsList options_list_arg) :
  options_list(move(options_list_arg))
{
}

void
RamseyModelStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.ramsey_model_present = true;

  /* Fill in option_order of mod_file_struct
     Since ramsey model needs one further order of derivation (for example, for 1st order
     approximation, it needs 2nd derivatives), we add 1 to the order declared by user */
  auto it = options_list.num_options.find("order");
  if (it != options_list.num_options.end())
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
RamseyModelStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  // options_.ramsey_policy indicates that a Ramsey model is present in the *.mod file
  // this affects the computation of the steady state that uses a special algorithm
  // It should probably rather be a M_ field, but we leave it in options_ for historical reason

  // Ensure that order 3 implies k_order (#844)
  auto it = options_list.num_options.find("order");
  auto it1 = options_list.num_options.find("k_order_solver");
  if ((it1 != options_list.num_options.end() && it1->second == "1")
      || (it != options_list.num_options.end() && stoi(it->second) >= 3))
    output << "options_.k_order_solver = 1;" << endl;

  output << "options_.ramsey_policy = 1;" << endl;
  options_list.writeOutput(output);
}

void
RamseyModelStatement::writeJsonOutput(ostream &output) const
{
  output << "{\"statementName\": \"ramsey_model\"";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

RamseyConstraintsStatement::RamseyConstraintsStatement(const SymbolTable &symbol_table_arg, constraints_t constraints_arg) :
  symbol_table(symbol_table_arg),
  constraints(move(constraints_arg))
{
}

void
RamseyConstraintsStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  if ((mod_file_struct.ramsey_model_present != true) || (mod_file_struct.ramsey_policy_present != true))
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
  output << "{\"statementName\": \"ramsey_constraints\""
         << ", \"ramsey_model_constraints\": [" << endl;
  for (auto it = constraints.begin(); it != constraints.end(); ++it)
    {
      if (it != constraints.begin())
        output << ", ";
      output << "{\"constraint\": \"" << symbol_table.getName(it->endo) << " ";
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
          exit(1);
        }
      output << " ";
      it->expression->writeJsonOutput(output, {}, {});
      output << "\"}" << endl;
    }
  output << "]" << endl;
  output << "}";
}

RamseyPolicyStatement::RamseyPolicyStatement(const SymbolTable &symbol_table_arg,
                                             vector<string> ramsey_policy_list_arg,
                                             OptionsList options_list_arg) :
  symbol_table(symbol_table_arg),
  ramsey_policy_list(move(ramsey_policy_list_arg)),
  options_list(move(options_list_arg))
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
  auto it = options_list.num_options.find("order");
  if (it != options_list.num_options.end())
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
RamseyPolicyStatement::checkRamseyPolicyList()
{
  for (const auto & it : ramsey_policy_list)
    {
      if (!symbol_table.exists(it))
        {
          cerr << "ERROR: ramsey_policy: " << it << " was not declared." << endl;
          exit(EXIT_FAILURE);
        }
      if (symbol_table.getType(it) != SymbolType::endogenous)
        {
          cerr << "ERROR: ramsey_policy: " << it << " is not endogenous." << endl;
          exit(EXIT_FAILURE);
        }
    }
}

void
RamseyPolicyStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  // Ensure that order 3 implies k_order (#844)
  auto it = options_list.num_options.find("order");
  auto it1 = options_list.num_options.find("k_order_solver");
  if ((it1 != options_list.num_options.end() && it1->second == "1")
      || (it != options_list.num_options.end() && stoi(it->second) >= 3))
    output << "options_.k_order_solver = 1;" << endl;

  options_list.writeOutput(output);
  output << "var_list_ = {";
  for (auto it = ramsey_policy_list.begin();
       it != ramsey_policy_list.end(); ++it)
    {
      if (it != ramsey_policy_list.begin())
        output << ";";
      output << "'" << *it << "'";
    }
  output << "};" << endl
         << "ramsey_policy(var_list_);" << endl;
}

void
RamseyPolicyStatement::writeJsonOutput(ostream &output) const
{
  output << "{\"statementName\": \"ramsey_policy\"";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << ", \"ramsey_policy_list\": [";
  for (auto it = ramsey_policy_list.begin();
       it != ramsey_policy_list.end(); ++it)
    {
      if (it != ramsey_policy_list.begin())
        output << ",";
      output << "\"" << *it << "\"";
    }
  output << "]"
         << "}";
}

DiscretionaryPolicyStatement::DiscretionaryPolicyStatement(SymbolList symbol_list_arg,
                                                           OptionsList options_list_arg) :
  symbol_list(move(symbol_list_arg)),
  options_list(move(options_list_arg))
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
  auto it = options_list.num_options.find("order");
  if (it != options_list.num_options.end())
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
DiscretionaryPolicyStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  // Ensure that order 3 implies k_order (#844)
  auto it = options_list.num_options.find("order");
  auto it1 = options_list.num_options.find("k_order_solver");
  if ((it1 != options_list.num_options.end() && it1->second == "1")
      || (it != options_list.num_options.end() && stoi(it->second) >= 3))
    output << "options_.k_order_solver = 1;" << endl;

  options_list.writeOutput(output);
  symbol_list.writeOutput("var_list_", output);
  output << "discretionary_policy(var_list_);" << endl;
}

void
DiscretionaryPolicyStatement::writeJsonOutput(ostream &output) const
{
  output << "{\"statementName\": \"discretionary_policy\"";
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
  symbol_list(move(symbol_list_arg)),
  options_list(move(options_list_arg))
{
}

void
EstimationStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.estimation_present = true;

  // Fill in option_order of mod_file_struct
  auto it = options_list.num_options.find("order");
  if (it != options_list.num_options.end())
    {
      int order = stoi(it->second);

      if (order > 2)
        {
          cerr << "ERROR: order > 2 is not supported in estimation" << endl;
          exit(EXIT_FAILURE);
        }

      mod_file_struct.order_option = max(mod_file_struct.order_option, order);
    }

  // Fill in mod_file_struct.partial_information
  it = options_list.num_options.find("partial_information");
  if (it != options_list.num_options.end() && it->second == "1")
    mod_file_struct.partial_information = true;

  // Fill in mod_file_struct.estimation_analytic_derivation
  it = options_list.num_options.find("analytic_derivation");
  if (it != options_list.num_options.end() && it->second == "1")
    mod_file_struct.estimation_analytic_derivation = true;

  it = options_list.num_options.find("dsge_var");
  if (it != options_list.num_options.end())
    // Fill in mod_file_struct.dsge_var_calibrated
    mod_file_struct.dsge_var_calibrated = it->second;

  // Fill in mod_file_struct.dsge_var_estimated
  auto it_str = options_list.string_options.find("dsge_var");
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
}

void
EstimationStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  options_list.writeOutput(output);

  // Special treatment for order option and particle filter
  auto it = options_list.num_options.find("order");
  if (it == options_list.num_options.end())
    output << "options_.order = 1;" << endl;
  else if (stoi(it->second) == 2)
    output << "options_.particle.status = 1;" << endl;

  // Do not check for the steady state in diffuse filter mode (#400)
  it = options_list.num_options.find("diffuse_filter");
  if (it != options_list.num_options.end() && it->second == "1")
    output << "options_.steadystate.nocheck = 1;" << endl;

  symbol_list.writeOutput("var_list_", output);
  output << "oo_recursive_=dynare_estimation(var_list_);" << endl;
}

void
EstimationStatement::writeJsonOutput(ostream &output) const
{
  output << "{\"statementName\": \"estimation\"";
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
  options_list(move(options_list_arg))
{
}

void
DynareSensitivityStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  auto it = options_list.num_options.find("identification");
  if (it != options_list.num_options.end()
      && it->second == "1")
    mod_file_struct.identification_present = true;
}

void
DynareSensitivityStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  options_list.writeOutput(output, "options_gsa");

  /* Ensure that nograph, nodisplay and graph_format are also set in top-level
     options_.
     \todo factorize this code between identification and dynare_sensitivity,
     and provide a generic mechanism for this situation (maybe using regexps) */
  auto it = options_list.num_options.find("nodisplay");
  if (it != options_list.num_options.end())
    output << "options_.nodisplay = " << it->second << ";" << endl;
  it = options_list.num_options.find("nograph");
  if (it != options_list.num_options.end())
    output << "options_.nograph = " << it->second << ";" << endl;
  auto it2 = options_list.string_options.find("graph_format");
  if (it2 != options_list.string_options.end())
    output << "options_.graph_format = '" << it2->second << "';" << endl;

  output << "dynare_sensitivity(options_gsa);" << endl;
}

void
DynareSensitivityStatement::writeJsonOutput(ostream &output) const
{
  output << "{\"statementName\": \"dynare_sensitivity\"";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

RplotStatement::RplotStatement(SymbolList symbol_list_arg) :
  symbol_list(move(symbol_list_arg))
{
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
  output << "{\"statementName\": \"rplot\"";
  if (!symbol_list.empty())
    {
      output << ", ";
      symbol_list.writeJsonOutput(output);
    }
  output << "}";
}

UnitRootVarsStatement::UnitRootVarsStatement()
= default;

void
UnitRootVarsStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "options_.diffuse_filter = 1;" << endl
         << "options_.steadystate.nocheck = 1;" << endl;
}

void
UnitRootVarsStatement::writeJsonOutput(ostream &output) const
{
  output << "{\"statementName\": \"unit_root_vars\", "
         << "\"diffuse_filter\": 1, "
         << "\"steady_state.nocheck\": 1}";
}

PeriodsStatement::PeriodsStatement(int periods_arg) : periods(periods_arg)
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
  output << "{\"statementName\": \"periods\", "
         << "\"periods\": " << periods << "}";
}

DsampleStatement::DsampleStatement(int val1_arg) : val1(val1_arg), val2(-1)
{
}

DsampleStatement::DsampleStatement(int val1_arg, int val2_arg) : val1(val1_arg), val2(val2_arg)
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
  output << "{\"statementName\": \"dsample\", "
         << "\"value1\": " << val1 << ", "
         << "\"value2\": " << val2 << "}";
}

EstimatedParamsStatement::EstimatedParamsStatement(vector<EstimationParams> estim_params_list_arg,
                                                   const SymbolTable &symbol_table_arg) :
  estim_params_list(move(estim_params_list_arg)),
  symbol_table(symbol_table_arg)
{
}

void
EstimatedParamsStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  for (const auto & it : estim_params_list)
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
  for (const auto & it : estim_params_list)
    {
      if (it.type == 3) // Correlation
        {
          // Use lexical ordering for the pair of symbols
          pair<string, string> x = it.name < it.name2 ? make_pair(it.name, it.name2) : make_pair(it.name2, it.name);

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
  for (const auto & it : estim_params_list)
    if (it.type == 2 && it.name != "dsge_prior_weight")
      mod_file_struct.estimated_parameters.insert(symbol_table.getID(it.name));
}

void
EstimatedParamsStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "estim_params_.var_exo = [];" << endl
         << "estim_params_.var_endo = [];" << endl
         << "estim_params_.corrx = [];" << endl
         << "estim_params_.corrn = [];" << endl
         << "estim_params_.param_vals = [];" << endl;

  for (const auto & it : estim_params_list)
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
          output << symb_id << " " << symbol_table.getTypeSpecificID(it.name2)+1;
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
  output << "{\"statementName\": \"estimated_params\", "
         << "\"params\": [";
  for (auto it = estim_params_list.begin(); it != estim_params_list.end(); it++)
    {
      if (it != estim_params_list.begin())
        output << ", ";
      output << "{";
      switch (it->type)
        {
        case 1:
          output << "\"var\": \"" << it->name << "\"";
          break;
        case 2:
          output << "\"param\": \"" << it->name << "\"";
          break;
        case 3:
          output << "\"var1\": \"" << it->name << "\","
                 << "\"var2\": \"" << it->name2 << "\"";
          break;
        }

      output << ", \"init_val\": \"";
      it->init_val->writeJsonOutput(output, {}, {});
      output << "\", \"lower_bound\": \"";
      it->low_bound->writeJsonOutput(output, {}, {});
      output << "\", \"upper_bound\": \"";
      it->up_bound->writeJsonOutput(output, {}, {});
      output << "\", \"prior_distribution\": "
             << static_cast<int>(it->prior)
             << ", \"mean\": \"";
      it->mean->writeJsonOutput(output, {}, {});
      output << "\", \"std\": \"";
      it->std->writeJsonOutput(output, {}, {});
      output << "\", \"p3\": \"";
      it->p3->writeJsonOutput(output, {}, {});
      output << "\", \"p4\": \"";
      it->p4->writeJsonOutput(output, {}, {});
      output << "\", \"jscale\": \"";
      it->jscale->writeJsonOutput(output, {}, {});
      output << "\"}" << endl;
    }
  output << "]"
         << "}";
}

EstimatedParamsInitStatement::EstimatedParamsInitStatement(vector<EstimationParams> estim_params_list_arg,
                                                           const SymbolTable &symbol_table_arg,
                                                           const bool use_calibration_arg) :
  estim_params_list(move(estim_params_list_arg)),
  symbol_table(symbol_table_arg),
  use_calibration(use_calibration_arg)
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

  for (const auto & it : estim_params_list)
    {
      int symb_id = symbol_table.getTypeSpecificID(it.name) + 1;
      SymbolType symb_type = symbol_table.getType(it.name);

      if (it.type < 3)
        {
          if (symb_type == SymbolType::exogenous)
            {
              output << "tmp1 = find(estim_params_.var_exo(:,1)==" << symb_id << ");" << endl;
              output << "estim_params_.var_exo(tmp1,2) = ";
              it.init_val->writeOutput(output);
              output << ";" << endl;
            }
          else if (symb_type == SymbolType::endogenous)
            {
              output << "tmp1 = find(estim_params_.var_endo(:,1)==" << symb_id << ");" << endl;
              output << "estim_params_.var_endo(tmp1,2) = ";
              it.init_val->writeOutput(output);
              output << ";" << endl;
            }
          else if (symb_type == SymbolType::parameter)
            {
              output << "tmp1 = find(estim_params_.param_vals(:,1)==" << symb_id << ");" << endl;
              output << "estim_params_.param_vals(tmp1,2) = ";
              it.init_val->writeOutput(output);
              output << ";" << endl;
            }
        }
      else
        {
          if (symb_type == SymbolType::exogenous)
            {
              output << "tmp1 = find((estim_params_.corrx(:,1)==" << symb_id << " & estim_params_.corrx(:,2)==" << symbol_table.getTypeSpecificID(it.name2)+1 << ") | "
                     <<             "(estim_params_.corrx(:,2)==" << symb_id << " & estim_params_.corrx(:,1)==" << symbol_table.getTypeSpecificID(it.name2)+1 << "));" << endl;
              output << "estim_params_.corrx(tmp1,3) = ";
              it.init_val->writeOutput(output);
              output << ";" << endl;
            }
          else if (symb_type == SymbolType::endogenous)
            {
              output << "tmp1 = find((estim_params_.corrn(:,1)==" << symb_id << " & estim_params_.corrn(:,2)==" << symbol_table.getTypeSpecificID(it.name2)+1 << ") | "
                     <<             "(estim_params_.corrn(:,2)==" << symb_id << " & estim_params_.corrn(:,1)==" << symbol_table.getTypeSpecificID(it.name2)+1 << "));" << endl;
              output << "estim_params_.corrn(tmp1,3) = ";
              it.init_val->writeOutput(output);
              output << ";" << endl;
            }
        }
    }
}

void
EstimatedParamsInitStatement::writeJsonOutput(ostream &output) const
{
  output << "{\"statementName\": \"estimated_params_init\"";

  if (use_calibration)
    output << ", \"use_calibration_initialization\": 1";

  output << ", \"params\": [";
  for (auto it = estim_params_list.begin(); it != estim_params_list.end(); it++)
    {
      if (it != estim_params_list.begin())
        output << ", ";
      output << "{";
      switch (it->type)
        {
        case 1:
          output << "\"var\": \"" << it->name << "\"";
          break;
        case 2:
          output << "\"param\": \"" << it->name << "\"";
          break;
        case 3:
          output << "\"var1\": \"" << it->name << "\","
                 << "\"var2\": \"" << it->name2 << "\"";
          break;
        }
      output << ", \"init_val\": \"";
      it->init_val->writeJsonOutput(output, {}, {});
      output << "\"}";
    }
  output << "]"
         << "}";
}

EstimatedParamsBoundsStatement::EstimatedParamsBoundsStatement(vector<EstimationParams> estim_params_list_arg,
                                                               const SymbolTable &symbol_table_arg) :
  estim_params_list(move(estim_params_list_arg)),
  symbol_table(symbol_table_arg)
{
}

void
EstimatedParamsBoundsStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  for (const auto & it : estim_params_list)
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
  output << "{\"statementName\": \"estimated_params_bounds\", "
         << "\"params\": [";

  for (auto it = estim_params_list.begin(); it != estim_params_list.end(); it++)
    {
      if (it != estim_params_list.begin())
        output << ", ";
      output << "{";
      switch (it->type)
        {
        case 1:
          output << "\"var\": \"" << it->name << "\"";
        case 2:
          output << "\"param\": \"" << it->name << "\"";
          break;
        case 3:
          output << "\"var1\": \"" << it->name << "\","
                 << "\"var2\": \"" << it->name2 << "\"";
          break;
        }
      output << ", \"lower_bound\": ";
      it->low_bound->writeJsonOutput(output, {}, {});
      output << ", \"upper_bound\": ";
      it->up_bound->writeJsonOutput(output, {}, {});
      output << "}";
    }
  output << "]"
         << "}";
}

ObservationTrendsStatement::ObservationTrendsStatement(trend_elements_t trend_elements_arg,
                                                       const SymbolTable &symbol_table_arg) :
  trend_elements(move(trend_elements_arg)),
  symbol_table(symbol_table_arg)
{
}

void
ObservationTrendsStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "options_.trend_coeff = {};" << endl;
  for (const auto & trend_element : trend_elements)
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
  output << "{\"statementName\": \"observation_trends\", "
         << "\"trends\" : {";
  bool printed = false;
  for (const auto & trend_element : trend_elements)
    {
      if (symbol_table.getType(trend_element.first) == SymbolType::endogenous)
        {
          if (printed)
            output << ", ";
          output << "\"" << trend_element.first << "\": \"";
          trend_element.second->writeJsonOutput(output, {}, {});
          output << "\"" << endl;
          printed = true;
        }
      else
        cerr << "Warning : Non-variable symbol used in observation_trends: " << trend_element.first << endl;
    }
  output << "}"
         << "}";
}

OsrParamsStatement::OsrParamsStatement(SymbolList symbol_list_arg, const SymbolTable &symbol_table_arg) :
  symbol_list(move(symbol_list_arg)),
  symbol_table(symbol_table_arg)
{
}

void
OsrParamsStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  if (mod_file_struct.osr_params_present)
    cerr << "WARNING: You have more than one osr_params statement in the .mod file." << endl;
  mod_file_struct.osr_params_present = true;
}

void
OsrParamsStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  symbol_list.writeOutput("M_.osr.param_names", output);
  output << "M_.osr.param_names = cellstr(M_.osr.param_names);" << endl
         << "M_.osr.param_indices = zeros(length(M_.osr.param_names), 1);" << endl;
  int i = 0;
  vector<string> symbols = symbol_list.get_symbols();
  for (vector<string>::const_iterator it = symbols.begin(); it != symbols.end(); it++)
    output << "M_.osr.param_indices(" << ++i <<") = " << symbol_table.getTypeSpecificID(*it) + 1 << ";" << endl;
}

void
OsrParamsStatement::writeJsonOutput(ostream &output) const
{
  output << "{\"statementName\": \"osr_params\"";
  if (!symbol_list.empty())
    {
      output << ", ";
      symbol_list.writeJsonOutput(output);
    }
  output << "}";
}

OsrStatement::OsrStatement(SymbolList symbol_list_arg,
                           OptionsList options_list_arg) :
  symbol_list(move(symbol_list_arg)),
  options_list(move(options_list_arg))
{
}

OsrParamsBoundsStatement::OsrParamsBoundsStatement(vector<OsrParams> osr_params_list_arg) :
  osr_params_list(move(osr_params_list_arg))
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

  for (const auto & it : osr_params_list)
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
  output << "{\"statementName\": \"osr_params_bounds\""
         << ", \"bounds\": [";
  for (auto it = osr_params_list.begin();
       it != osr_params_list.end(); it++)
    {
      if (it != osr_params_list.begin())
        output << ", ";
      output << "{\"parameter\": \"" << it->name << "\","
             << "\"bounds\": [\"";
      it->low_bound->writeJsonOutput(output, {}, {});
      output << "\", \"";
      it->up_bound->writeJsonOutput(output, {}, {});
      output << "\"]"
             << "}";
    }
  output << "]"
         << "}";
}

void
OsrStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.osr_present = true;

  // Fill in option_order of mod_file_struct
  auto it = options_list.num_options.find("order");
  if (it != options_list.num_options.end())
    mod_file_struct.order_option = max(mod_file_struct.order_option, stoi(it->second));

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
OsrStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  // Ensure that order 3 implies k_order (#844)
  auto it = options_list.num_options.find("order");
  auto it1 = options_list.num_options.find("k_order_solver");
  if ((it1 != options_list.num_options.end() && it1->second == "1")
      || (it != options_list.num_options.end() && stoi(it->second) >= 3))
    output << "options_.k_order_solver = 1;" << endl;

  options_list.writeOutput(output);
  symbol_list.writeOutput("var_list_", output);
  output << "oo_.osr = osr(var_list_,M_.osr.param_names,M_.osr.variable_indices,M_.osr.variable_weights);" << endl;
}

void
OsrStatement::writeJsonOutput(ostream &output) const
{
  output << "{\"statementName\": \"osr\"";
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
  var_weights(move(var_weights_arg)),
  covar_weights(move(covar_weights_arg)),
  symbol_table(symbol_table_arg)
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

  for (const auto & var_weight : var_weights)
    {
      const string &name = var_weight.first;
      const expr_t value = var_weight.second;
      int id = symbol_table.getTypeSpecificID(name) + 1;
      output << "M_.osr.variable_weights(" << id << "," << id << ") = ";
      value->writeOutput(output);
      output << ";" << endl;
      output << "M_.osr.variable_indices = [M_.osr.variable_indices; " << id << "];" << endl;
    }

  for (const auto & covar_weight : covar_weights)
    {
      const string &name1 = covar_weight.first.first;
      const string &name2 = covar_weight.first.second;
      const expr_t value = covar_weight.second;
      int id1 = symbol_table.getTypeSpecificID(name1) + 1;
      int id2 = symbol_table.getTypeSpecificID(name2) + 1;
      output << "M_.osr.variable_weights(" << id1 << "," << id2 << ") = ";
      value->writeOutput(output);
      output << ";" << endl;
      output << "M_.osr.variable_indices = [M_.osr.variable_indices; " << id1 << "; " << id2 << "];" << endl;
    }
}

void
OptimWeightsStatement::writeJsonOutput(ostream &output) const
{
  output << "{\"statementName\": \"optim_weights\", "
         << "\"weights\": [";
  for (auto it = var_weights.begin();
       it != var_weights.end(); it++)
    {
      if (it != var_weights.begin())
        output << ", ";
      output << "{\"name\": \"" << it->first << "\""
             << ", \"value\": \"";
      it->second->writeJsonOutput(output, {}, {});
      output << "\"}";
    }

  for (auto it = covar_weights.begin();
       it != covar_weights.end(); it++)
    {
      if (it != covar_weights.begin() || !var_weights.empty())
        output << ", ";
      output << "{\"name1\": \"" << it->first.first << "\""
             << ", \"name2\": \"" << it->first.second << "\""
             << ", \"value\": \"";
      it->second->writeJsonOutput(output, {}, {});
      output << "\"}";
    }
  output << "]"
         << "}";
}

DynaSaveStatement::DynaSaveStatement(SymbolList symbol_list_arg,
                                     string filename_arg) :
  symbol_list(move(symbol_list_arg)),
  filename(move(filename_arg))
{
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
  output << "{\"statementName\": \"dynasave\", "
         << "\"filename\": \"" << filename << "\"";
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
DynaTypeStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  symbol_list.writeOutput("var_list_", output);
  output << "dynatype('" << filename
         << "',var_list_);" << endl;
}

void
DynaTypeStatement::writeJsonOutput(ostream &output) const
{
  output << "{\"statementName\": \"dynatype\", "
         << "\"filename\": \"" << filename << "\"";
  if (!symbol_list.empty())
    {
      output << ", ";
      symbol_list.writeJsonOutput(output);
    }
  output << "}";
}

ModelComparisonStatement::ModelComparisonStatement(filename_list_t filename_list_arg,
                                                   OptionsList options_list_arg) :
  filename_list(move(filename_list_arg)),
  options_list(move(options_list_arg))
{
}

void
ModelComparisonStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  options_list.writeOutput(output);

  output << "ModelNames_ = {};" << endl;
  output << "ModelPriors_ = [];" << endl;

  for (const auto & it : filename_list)
    {
      output << "ModelNames_ = { ModelNames_{:} '" << it.first << "'};" << endl;
      output << "ModelPriors_ = [ ModelPriors_ ; " << it.second << "];" << endl;
    }
  output << "oo_ = model_comparison(ModelNames_,ModelPriors_,oo_,options_,M_.fname);" << endl;
}

void
ModelComparisonStatement::writeJsonOutput(ostream &output) const
{
  output << "{\"statementName\": \"model_comparison\"";
  if (!filename_list.empty())
    output << ", \"filename_list\": {";

  for (auto it = filename_list.begin();
       it != filename_list.end(); it++)
    {
      if (it != filename_list.begin())
        output << ", ";
      output << "\"name\": \"" << it->first << "\""
             << "\"prior\": \"" << it->second << "\"";
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

PlannerObjectiveStatement::PlannerObjectiveStatement(SymbolTable &symbol_table,
                                                     NumericalConstants &num_constants,
                                                     ExternalFunctionsTable &external_functions_table,
                                                     TrendComponentModelTable &trend_component_model_table_arg) :
  model_tree{symbol_table, num_constants,
             external_functions_table, trend_component_model_table_arg}
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

StaticModel &
PlannerObjectiveStatement::getPlannerObjective()
{
  return model_tree;
}

void
PlannerObjectiveStatement::computingPass()
{
  model_tree.computingPass({}, false, true, true, 0, false, false, false);
  computing_pass_called = true;
}

void
PlannerObjectiveStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  model_tree.writeStaticFile(basename + ".objective", false, false, false, false);
}

void
PlannerObjectiveStatement::writeJsonOutput(ostream &output) const
{
  output << "{\"statementName\": \"planner_objective\""
         << ", ";
  if (computing_pass_called)
    model_tree.writeJsonComputingPassOutput(output, false);
  else
    model_tree.writeJsonOutput(output);

  output << "}";
}

BVARDensityStatement::BVARDensityStatement(int maxnlags_arg, OptionsList options_list_arg) :
  maxnlags(maxnlags_arg),
  options_list(move(options_list_arg))
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
  output << "{\"statementName\": \"bvar_density\"";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

BVARForecastStatement::BVARForecastStatement(int nlags_arg, OptionsList options_list_arg) :
  nlags(nlags_arg),
  options_list(move(options_list_arg))
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
  output << "{\"statementName\": \"bvar_forecast\"";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

SBVARStatement::SBVARStatement(OptionsList options_list_arg) :
  options_list(move(options_list_arg))
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
  output << "{\"statementName\": \"sbvar\"";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

MSSBVAREstimationStatement::MSSBVAREstimationStatement(OptionsList options_list_arg) :
  options_list(move(options_list_arg))
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
  output << "{\"statementName\": \"ms_sbvar_estimation\"";
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
  auto mh_replic_it = options_list.num_options.find("ms.mh_replic");
  auto thinning_factor_it = options_list.num_options.find("ms.thinning_factor");
  auto drop_it = options_list.num_options.find("ms.drop");
  if (mh_replic_it != options_list.num_options.end() || thinning_factor_it != options_list.num_options.end())
    if (drop_it == options_list.num_options.end())
      output << "options_.ms.drop = 0.1*options_.ms.mh_replic*options_.ms.thinning_factor;" << endl;

  output << "[options_, oo_] = ms_simulation(M_, options_, oo_);" << endl;
}

void
MSSBVARSimulationStatement::writeJsonOutput(ostream &output) const
{
  output << "{\"statementName\": \"ms_sbvar_simulation\"";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

MSSBVARComputeMDDStatement::MSSBVARComputeMDDStatement(OptionsList options_list_arg) :
  options_list(move(options_list_arg))
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
  output << "{\"statementName\": \"ms_sbvar_compute_mdd\"";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

MSSBVARComputeProbabilitiesStatement::MSSBVARComputeProbabilitiesStatement(OptionsList options_list_arg) :
  options_list(move(options_list_arg))
{
}

void
MSSBVARComputeProbabilitiesStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
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
MSSBVARComputeProbabilitiesStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "options_ = initialize_ms_sbvar_options(M_, options_);" << endl;
  options_list.writeOutput(output);
  output << "[options_, oo_] = ms_compute_probabilities(M_, options_, oo_);" << endl;
}

void
MSSBVARComputeProbabilitiesStatement::writeJsonOutput(ostream &output) const
{
  output << "{\"statementName\": \"ms_sbvar_compute_probabilities\"";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

MSSBVARIrfStatement::MSSBVARIrfStatement(SymbolList symbol_list_arg,
                                         OptionsList options_list_arg) :
  symbol_list(move(symbol_list_arg)),
  options_list(move(options_list_arg))
{
}

void
MSSBVARIrfStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.bvar_present = true;

  bool regime_present = false;
  bool regimes_present = false;
  bool filtered_probabilities_present = false;

  auto it = options_list.num_options.find("ms.regimes");
  if (it != options_list.num_options.end())
    regimes_present = true;

  it = options_list.num_options.find("ms.regime");
  if (it != options_list.num_options.end())
    regime_present = true;

  it = options_list.num_options.find("ms.filtered_probabilities");
  if (it != options_list.num_options.end())
    filtered_probabilities_present = true;

  if ((filtered_probabilities_present && regime_present)
      || (filtered_probabilities_present && regimes_present)
      || (regimes_present && regime_present))
    {
      cerr << "ERROR: You may only pass one of regime, regimes and "
           << "filtered_probabilities to ms_irf" << endl;
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
  output << "{\"statementName\": \"ms_sbvar_irf\"";
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
  options_list(move(options_list_arg))
{
}

void
MSSBVARForecastStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.bvar_present = true;

  if (options_list.num_options.find("ms.regimes") != options_list.num_options.end())
    if (options_list.num_options.find("ms.regime") != options_list.num_options.end())
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
  output << "{\"statementName\": \"ms_sbvar_forecast\"";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

MSSBVARVarianceDecompositionStatement::MSSBVARVarianceDecompositionStatement(OptionsList options_list_arg) :
  options_list(move(options_list_arg))
{
}

void
MSSBVARVarianceDecompositionStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.bvar_present = true;

  bool regime_present = false;
  bool regimes_present = false;
  bool filtered_probabilities_present = false;

  auto it = options_list.num_options.find("ms.regimes");
  if (it != options_list.num_options.end())
    regimes_present = true;

  it = options_list.num_options.find("ms.regime");
  if (it != options_list.num_options.end())
    regime_present = true;

  it = options_list.num_options.find("ms.filtered_probabilities");
  if (it != options_list.num_options.end())
    filtered_probabilities_present = true;

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
  output << "{\"statementName\": \"ms_sbvar_variance_decomposition\"";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

IdentificationStatement::IdentificationStatement(const OptionsList &options_list_arg)
{
  options_list = options_list_arg;
  if (options_list.num_options.find("max_dim_cova_group") != options_list.num_options.end())
    if (stoi(options_list.num_options["max_dim_cova_group"]) == 0)
      {
        cerr << "ERROR: The max_dim_cova_group option to identification only accepts integers > 0." << endl;
        exit(EXIT_FAILURE);
      }
}

void
IdentificationStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.identification_present = true;
}

void
IdentificationStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  options_list.writeOutput(output, "options_ident");

  /* Ensure that nograph, nodisplay and graph_format are also set in top-level
     options_.
     \todo factorize this code between identification and dynare_sensitivity,
     and provide a generic mechanism for this situation (maybe using regexps) */
  auto it = options_list.num_options.find("nodisplay");
  if (it != options_list.num_options.end())
    output << "options_.nodisplay = " << it->second << ";" << endl;
  it = options_list.num_options.find("nograph");
  if (it != options_list.num_options.end())
    output << "options_.nograph = " << it->second << ";" << endl;
  auto it2 = options_list.string_options.find("graph_format");
  if (it2 != options_list.string_options.end())
    output << "options_.graph_format = '" << it2->second << "';" << endl;

  output << "dynare_identification(options_ident);" << endl;
}

void
IdentificationStatement::writeJsonOutput(ostream &output) const
{
  output << "{\"statementName\": \"identification\"";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

WriteLatexDynamicModelStatement::WriteLatexDynamicModelStatement(const DynamicModel &dynamic_model_arg, bool write_equation_tags_arg) :
  dynamic_model(dynamic_model_arg),
  write_equation_tags(write_equation_tags_arg)
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
  output << "{\"statementName\": \"write_latex_dynamic_model\"}";
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
  output << "{\"statementName\": \"write_latex_static_model\"}";
}

WriteLatexOriginalModelStatement::WriteLatexOriginalModelStatement(const DynamicModel &original_model_arg, bool write_equation_tags_arg) :
  original_model(original_model_arg),
  write_equation_tags(write_equation_tags_arg)
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
  output << "{\"statementName\": \"write_latex_original_model\"}";
}

WriteLatexSteadyStateModelStatement::WriteLatexSteadyStateModelStatement(const SteadyStateModel &steady_state_model_arg) :
  steady_state_model(steady_state_model_arg)
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
  output << "{\"statementName\": \"write_latex_steady_state_model\"}";
}

ShockDecompositionStatement::ShockDecompositionStatement(SymbolList symbol_list_arg,
                                                         OptionsList options_list_arg) :
  symbol_list(move(symbol_list_arg)),
  options_list(move(options_list_arg))
{
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
  output << "{\"statementName\": \"shock_decomposition\"";
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
  symbol_list(move(symbol_list_arg)),
  options_list(move(options_list_arg))
{
}

void
RealtimeShockDecompositionStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  options_list.writeOutput(output);
  symbol_list.writeOutput("var_list_", output);
  output << "oo_ = realtime_shock_decomposition(M_,oo_,options_,var_list_,bayestopt_,estim_params_);" << endl;
}

PlotShockDecompositionStatement::PlotShockDecompositionStatement(SymbolList symbol_list_arg,
                                                                 OptionsList options_list_arg) :
  symbol_list(move(symbol_list_arg)),
  options_list(move(options_list_arg))
{
}

void
PlotShockDecompositionStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "options_ = set_default_plot_shock_decomposition_options(options_);" << endl;
  options_list.writeOutput(output);
  symbol_list.writeOutput("var_list_", output);
  output << "plot_shock_decomposition(M_, oo_, options_, var_list_);" << endl;
}

InitialConditionDecompositionStatement::InitialConditionDecompositionStatement(SymbolList symbol_list_arg,
                                                                               OptionsList options_list_arg) :
  symbol_list(move(symbol_list_arg)),
  options_list(move(options_list_arg))
{
}

void
InitialConditionDecompositionStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "options_ = set_default_initial_condition_decomposition_options(options_);" << endl;
  options_list.writeOutput(output);
  symbol_list.writeOutput("var_list_", output);
  output << "oo_ = initial_condition_decomposition(M_, oo_, options_, var_list_, bayestopt_, estim_params_);" << endl;
}

ConditionalForecastStatement::ConditionalForecastStatement(OptionsList options_list_arg) :
  options_list(move(options_list_arg))
{
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
  output << "{\"statementName\": \"conditional_forecast\"";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

PlotConditionalForecastStatement::PlotConditionalForecastStatement(int periods_arg, SymbolList symbol_list_arg) :
  periods(periods_arg),
  symbol_list(move(symbol_list_arg))
{
}

void
PlotConditionalForecastStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  symbol_list.writeOutput("var_list_", output);
  if (periods == -1)
    output << "plot_icforecast(var_list_,[],options_);" << endl;
  else
    output << "plot_icforecast(var_list_, " << periods << ",options_);" << endl;
}

void
PlotConditionalForecastStatement::writeJsonOutput(ostream &output) const
{
  output << "{\"statementName\": \"plot_conditional_forecast\", "
         << "\"periods\": " << periods;
  if (!symbol_list.empty())
    {
      output << ", ";
      symbol_list.writeJsonOutput(output);
    }
  output << "}";
}

SvarIdentificationStatement::SvarIdentificationStatement(svar_identification_restrictions_t restrictions_arg,
                                                         const bool &upper_cholesky_present_arg,
                                                         const bool &lower_cholesky_present_arg,
                                                         const bool &constants_exclusion_present_arg,
                                                         const SymbolTable &symbol_table_arg) :
  restrictions(move(restrictions_arg)),
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
  for (const auto & restriction : restrictions)
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
      output << "options_.ms.Qi = cell(" << n << ",1);" << endl;
      output << "options_.ms.Ri = cell(" << n << ",1);" << endl;

      for (auto it = restrictions.begin(); it != restrictions.end(); it++)
        {
          assert(it->lag >= 0);
          if (it->lag == 0)
            output << "options_.ms.Qi{" << it->equation << "}(" << it->restriction_nbr << ", " << it->variable + 1 << ") = ";
          else
            {
              int col = (it->lag-1)*n+it->variable+1;
              if (col > k)
                {
                  cerr << "ERROR: lag =" << it->lag << ", num endog vars = " << n << "current endog var index = " << it->variable << ". Index "
                       << "out of bounds. If the above does not represent a logical error, please report this to the Dyanre Team." << endl;
                  exit(EXIT_FAILURE);
                }
              output << "options_.ms.Ri{" << it->equation << "}(" << it->restriction_nbr << ", " << col << ") = ";
            }
          it->value->writeOutput(output);
          output << ";" << endl;
        }
      output << "options_.ms.nlags = " << r << ";" << endl;
    }
}

void
SvarIdentificationStatement::writeJsonOutput(ostream &output) const
{
  output << "{\"statementName\": \"svar_identification\"";

  if (upper_cholesky_present)
    output << ", \"upper_cholesky\": 1";

  if (lower_cholesky_present)
    output << ", \"lower_cholesky\": 1";

  if (constants_exclusion_present)
    output << ", \"constants_exclusion\": 1";

  if (!upper_cholesky_present && !lower_cholesky_present)
    {
      output << ", \"nlags\": " << getMaxLag()
             << ", \"restrictions\": [";

      for (auto it = restrictions.begin(); it != restrictions.end(); it++)
        {
          if (it != restrictions.begin())
            output << ", ";
          output << "{"
                 << "\"equation_number\": " << it->equation << ", "
                 << "\"restriction_number\": " << it->restriction_nbr << ", "
                 << "\"variable\": \"" << symbol_table.getName(it->variable) << "\", "
                 << "\"expression\": \"";
          it->value->writeOutput(output);
          output << "\"}";
        }
      output << "]";
    }
  output << "}";
}

MarkovSwitchingStatement::MarkovSwitchingStatement(OptionsList options_list_arg) :
  options_list(move(options_list_arg))
{
  auto it_num = options_list.num_options.find("ms.restrictions");
  if (it_num != options_list.num_options.end())
    {
      using namespace boost;
      auto it_num_regimes
        = options_list.num_options.find("ms.number_of_regimes");
      assert(it_num_regimes !=  options_list.num_options.end());
      auto num_regimes = stoi(it_num_regimes->second);

      vector<string> tokenizedRestrictions;
      split(tokenizedRestrictions, it_num->second, is_any_of("["), token_compress_on);
      for (auto & tokenizedRestriction : tokenizedRestrictions)
        if (tokenizedRestriction.size() > 0)
          {
            vector<string> restriction;
            split(restriction, tokenizedRestriction, is_any_of("], "));
            for (auto it1 = restriction.begin();
                 it1 != restriction.end();)
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

  auto it_num = options_list.num_options.find("ms.restrictions");
  if (it_num != options_list.num_options.end())
    {
      using namespace boost;
      auto it_num_regimes
        = options_list.num_options.find("ms.number_of_regimes");
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
MarkovSwitchingStatement::writeCOutput(ostream &output, const string &basename)
{
  output << endl;

  auto it
    = options_list.num_options.find("ms.chain");
  assert(it !=  options_list.num_options.end());
  output << "chain = " << it->second << ";" << endl;

  it = options_list.num_options.find("ms.number_of_regimes");
  assert(it !=  options_list.num_options.end());
  output << "number_of_regimes = " << it->second << ";" << endl;

  it = options_list.num_options.find("ms.number_of_lags");
  if (it !=  options_list.num_options.end())
    output << "number_of_lags = " << it->second << ";" << endl
           << "number_of_lags_was_passed = true;" << endl;
  else
    output << "number_of_lags_was_passed = false;" << endl;

  it = options_list.num_options.find("ms.duration");
  assert(it != options_list.num_options.end());
  output << "duration.clear();" << endl;
  using namespace boost;
  vector<string> tokenizedDomain;
  split(tokenizedDomain, it->second, is_any_of("[ ]"), token_compress_on);
  for (auto & itvs : tokenizedDomain)
    if (!itvs.empty())
      output << "duration.push_back(" << itvs << ");" << endl;

  auto itsl
    = options_list.symbol_list_options.find("ms.parameters");
  assert(itsl != options_list.symbol_list_options.end());
  vector<string> parameters = itsl->second.get_symbols();
  output << "parameters.clear();" << endl;
  for (auto & parameter : parameters)
    output << "parameters.push_back(param_names[\"" << parameter << "\"]);" << endl;

  output << "restriction_map.clear();" << endl;
  for (auto & itrm : restriction_map)
    output << "restriction_map[make_pair(" << itrm.first.first << ","
           << itrm.first.second << ")] = " << itrm.second << ";" << endl;

  output << "msdsgeinfo->addMarkovSwitching(new MarkovSwitching(" << endl
         << "     chain, number_of_regimes, number_of_lags, number_of_lags_was_passed, parameters, duration, restriction_map));" << endl;
}

void
MarkovSwitchingStatement::writeJsonOutput(ostream &output) const
{
  output << "{\"statementName\": \"markov_switching\"";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }

  if (!restriction_map.empty())
    output << ", {";
  for (auto it = restriction_map.begin();
       it != restriction_map.end(); it++)
    {
      if (it != restriction_map.begin())
        output << ", ";
      output << "{\"current_period_regime\": " << it->first.first
             << ", \"next_period_regime\": " << it->first.second
             << ", \"transition_probability\": "<< it->second
             << "}";
    }
  if (!restriction_map.empty())
    output << "}";
  output << "}";
}

SvarStatement::SvarStatement(OptionsList options_list_arg) :
  options_list(move(options_list_arg))
{
}

void
SvarStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  OptionsList::num_options_t::const_iterator it0, it1, it2;
  it0 = options_list.string_options.find("ms.coefficients");
  it1 = options_list.string_options.find("ms.variances");
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
  OptionsList::num_options_t::const_iterator it0, it1, it2;
  OptionsList::vec_int_options_t::const_iterator itv;

  it0 = options_list.num_options.find("ms.chain");
  assert(it0 != options_list.num_options.end());
  output << "options_.ms.ms_chain(" << it0->second << ")";

  it0 = options_list.string_options.find("ms.coefficients");
  it1 = options_list.string_options.find("ms.variances");
  it2 = options_list.string_options.find("ms.constants");

  if (it0 != options_list.string_options.end())
    output << "." << it0->second;
  else if (it1 != options_list.string_options.end())
    output << "." << it1->second;
  else
    output << "." << it2->second;

  itv = options_list.vector_int_options.find("ms.equations");
  output << ".equations = ";
  if (itv != options_list.vector_int_options.end())
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
  output << "{\"statementName\": \"svar\"";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

SvarGlobalIdentificationCheckStatement::SvarGlobalIdentificationCheckStatement()
= default;

void
SvarGlobalIdentificationCheckStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "svar_global_identification_check(options_);" << std::endl;
}

void
SvarGlobalIdentificationCheckStatement::writeJsonOutput(ostream &output) const
{
  output << "{\"statementName\": \"svar_global_identification\"}";
}

SetTimeStatement::SetTimeStatement(OptionsList options_list_arg) :
  options_list(move(options_list_arg))
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
  output << "{\"statementName\": \"set_time\"";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

EstimationDataStatement::EstimationDataStatement(OptionsList options_list_arg) :
  options_list(move(options_list_arg))
{
}

void
EstimationDataStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.estimation_data_statement_present = true;

  auto it = options_list.num_options.find("nobs");
  if (it != options_list.num_options.end())
    if (stoi(it->second) <= 0)
      {
        cerr << "ERROR: The nobs option of the data statement only accepts positive integers." << endl;
        exit(EXIT_FAILURE);
      }

  if ((options_list.string_options.find("file") == options_list.string_options.end())
      && (options_list.string_options.find("series") == options_list.string_options.end()))
    {
      cerr << "ERROR: The file or series option must be passed to the data statement." << endl;
      exit(EXIT_FAILURE);
    }

  if ((options_list.string_options.find("file") != options_list.string_options.end())
      && (options_list.string_options.find("series") != options_list.string_options.end()))
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
  output << "{\"statementName\": \"estimation_data\"";
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
  name1(move(name1_arg)),
  name2(move(name2_arg)),
  subsample_declaration_map(subsample_declaration_map_arg),
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
       it != subsample_declaration_map.end(); it++, map_indx++)
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
  output << "{\"statementName\": \"subsamples\""
         << ", \"name1\": \"" << name1 << "\"";
  if (!name2.empty())
    output << ", \"name2\": \"" << name2 << "\"";

  output << ", \"declarations\": {";
  for (auto it = subsample_declaration_map.begin();
       it != subsample_declaration_map.end(); it++)
    {
      if (it != subsample_declaration_map.begin())
        output << ",";
      output << "{"
             << "\"range_index\": \"" << it->first << "\""
             << ", \"date1\": \"" << it->second.first << "\""
             << ", \"date2\": \"" << it->second.second << "\""
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
  to_name1(move(to_name1_arg)),
  to_name2(move(to_name2_arg)),
  from_name1(move(from_name1_arg)),
  from_name2(move(from_name2_arg)),
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
  output << "{\"statementName\": \"subsamples_equal\""
         << ", \"to_name1\": \"" << to_name1 << "\"";
  if (!to_name2.empty())
    output << ", \"to_name2\": \"" << to_name2 << "\"";
  output << ", \"from_name1\": \"" << from_name1 << "\"";
  if (!from_name2.empty())
    output << ", \"from_name2\": \"" << from_name2 << "\"";
  output << "}";
}

JointPriorStatement::JointPriorStatement(vector<string> joint_parameters_arg,
                                         PriorDistributions prior_shape_arg,
                                         OptionsList options_list_arg) :
  joint_parameters{move(joint_parameters_arg)},
  prior_shape(prior_shape_arg),
  options_list(move(options_list_arg))
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

  auto it_num = options_list.num_options.find("domain");
  if (it_num != options_list.num_options.end())
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
  for (const auto & joint_parameter : joint_parameters)
    output << "eifind = get_new_or_existing_ei_index('joint_parameter_prior_index', '"
           << joint_parameter << "', '');" << endl
           << "estimation_info.joint_parameter_prior_index(eifind) = {'" << joint_parameter << "'};" << endl;

  output << "key = {[";
  for (const auto & joint_parameter : joint_parameters)
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
  auto itn = options_list.num_options.find(field);
  output << lhs_field << "." << field << " = {";
  if (field == "variance")
    output << "{";
  if (itn != options_list.num_options.end())
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
  output << "{\"statementName\": \"joint_prior\""
         << ", \"key\": [";
  for (auto it = joint_parameters.begin(); it != joint_parameters.end(); it++)
    {
      if (it != joint_parameters.begin())
        output << ", ";
      output << "\"" << *it << "\"";
    }
  output << "]";

  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }

  output << ", \"shape\": ";
  switch (prior_shape)
    {
    case PriorDistributions::beta:
      output << "\"beta\"";
      break;
    case PriorDistributions::gamma:
      output << "\"gamma\"";
      break;
    case PriorDistributions::normal:
      output << "\"normal\"";
      break;
    case PriorDistributions::invGamma:
      output << "\"inv_gamma\"";
      break;
    case PriorDistributions::uniform:
      output << "\"uniform\"";
      break;
    case PriorDistributions::invGamma2:
      output << "\"inv_gamma2\"";
      break;
    case PriorDistributions::dirichlet:
      output << "\"dirichlet\"";
      break;
    case PriorDistributions::weibull:
      output << "\"weibull\"";
      break;
    case PriorDistributions::noShape:
      cerr << "Impossible case." << endl;
      exit(EXIT_FAILURE);
    }
  output << "}";
}

BasicPriorStatement::~BasicPriorStatement()
= default;

BasicPriorStatement::BasicPriorStatement(string name_arg,
                                         string subsample_name_arg,
                                         PriorDistributions prior_shape_arg,
                                         expr_t variance_arg,
                                         OptionsList options_list_arg) :
  name(move(name_arg)),
  subsample_name(move(subsample_name_arg)),
  prior_shape(prior_shape_arg),
  variance(variance_arg),
  options_list(move(options_list_arg))
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

  auto it_stdev = options_list.num_options.find("stdev");
  if ((it_stdev == options_list.num_options.end() && variance == nullptr)
      || (it_stdev != options_list.num_options.end() && variance != nullptr))
    {
      cerr << "ERROR: You must pass exactly one of stdev and variance to the prior statement." << endl;
      exit(EXIT_FAILURE);
    }

  auto it_num = options_list.num_options.find("domain");
  if (it_num != options_list.num_options.end())
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
  auto itn = options_list.num_options.find(field);
  if (itn != options_list.num_options.end())
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
  output << ", \"name\": \"" << name << "\""
         << ", \"subsample\": \"" << subsample_name << "\""
         << ", ";
  writeJsonShape(output);
  if (variance != nullptr)
    {
      output << ", \"variance\": \"";
      variance->writeJsonOutput(output, {}, {});
      output << "\"";
    }
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
}

void
BasicPriorStatement::writeCVarianceOption(ostream &output) const
{
  output << "variance = ";
  if (variance)
    variance->writeOutput(output);
  else
    output << "numeric_limits<double>::quiet_NaN()";
  output << ";" << endl;
}

void
BasicPriorStatement::writeCDomain(ostream &output) const
{
  output << "domain.clear();" << endl;
  auto it_num = options_list.num_options.find("domain");
  if (it_num != options_list.num_options.end())
    {
      using namespace boost;
      vector<string> tokenizedDomain;
      split(tokenizedDomain, it_num->second, is_any_of("[ ]"), token_compress_on);
      for (auto & it : tokenizedDomain)
        if (!it.empty())
          output << "domain.push_back(" << it << ");" << endl;
    }
}

void
BasicPriorStatement::writeCOutputHelper(ostream &output, const string &field) const
{
  auto itn = options_list.num_options.find(field);
  if (itn != options_list.num_options.end())
    output << field << " = " << itn->second << ";" << endl;
  else
    output << field << " = " << "numeric_limits<double>::quiet_NaN();" << endl;
}

void
BasicPriorStatement::writeCShape(ostream &output) const
{
  output << "shape = ";
  switch (prior_shape)
    {
    case PriorDistributions::beta:
      output << "\"beta\";" << endl;
      break;
    case PriorDistributions::gamma:
      output << "\"gamma\";" << endl;
      break;
    case PriorDistributions::normal:
      output << "\"normal\";" << endl;
      break;
    case PriorDistributions::invGamma:
      output << "\"inv_gamma\";" << endl;
      break;
    case PriorDistributions::uniform:
      output << "\"uniform\";" << endl;
      break;
    case PriorDistributions::invGamma2:
      output << "\"inv_gamma2\";" << endl;
      break;
    case PriorDistributions::dirichlet:
      output << "\"dirichlet\";" << endl;
      break;
    case PriorDistributions::weibull:
      output << "\"weibull\";" << endl;
      break;
    case PriorDistributions::noShape:
      assert(prior_shape != PriorDistributions::noShape);
    }
}

void
BasicPriorStatement::writeJsonShape(ostream &output) const
{
  output << "\"shape\": ";
  switch (prior_shape)
    {
    case PriorDistributions::beta:
      output << "\"beta\"";
      break;
    case PriorDistributions::gamma:
      output << "\"gamma\"";
      break;
    case PriorDistributions::normal:
      output << "\"normal\"";
      break;
    case PriorDistributions::invGamma:
      output << "\"inv_gamma\"";
      break;
    case PriorDistributions::uniform:
      output << "\"uniform\"";
      break;
    case PriorDistributions::invGamma2:
      output << "\"inv_gamma2\"";
      break;
    case PriorDistributions::dirichlet:
      output << "\"dirichlet\"";
      break;
    case PriorDistributions::weibull:
      output << "\"weibull\"";
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
  output << "{\"statementName\": \"prior\"";
  writeJsonPriorOutput(output);
  output << "}";
}

void
PriorStatement::writeCOutput(ostream &output, const string &basename)
{
  output << endl
         << "index = param_names[\""<< name << "\"];" << endl;
  writeCShape(output);
  writeCOutputHelper(output, "mean");
  writeCOutputHelper(output, "mode");
  writeCOutputHelper(output, "stdev");
  writeCVarianceOption(output);
  writeCDomain(output);

  output << "msdsgeinfo->addPrior(new ModFilePrior(" << endl
         << "     index, shape, mean, mode, stdev, variance, domain));" << endl;
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
  output << "{\"statementName\": \"std_prior\"";
  writeJsonPriorOutput(output);
  output << "}";
}

void
StdPriorStatement::writeCOutput(ostream &output, const string &basename)
{
  output << endl
         << "index = ";
  if (is_structural_innovation(symbol_table.getType(name)))
    output << "exo_names";
  else
    output << "endo_names";
  output << "[\""<< name << "\"];" << endl;

  writeCShape(output);
  writeCOutputHelper(output, "mean");
  writeCOutputHelper(output, "mode");
  writeCOutputHelper(output, "stdev");
  writeCVarianceOption(output);
  writeCDomain(output);

  if (is_structural_innovation(symbol_table.getType(name)))
    output << "msdsgeinfo->addStructuralInnovationPrior(new ModFileStructuralInnovationPrior(";
  else
    output << "msdsgeinfo->addMeasurementErrorPrior(new ModFileMeasurementErrorPrior(";
  output << endl << "     index, shape, mean, mode, stdev, variance, domain));" << endl;
}

CorrPriorStatement::CorrPriorStatement(string name_arg1, string name_arg2,
                                       string subsample_name_arg,
                                       PriorDistributions prior_shape_arg,
                                       expr_t variance_arg,
                                       OptionsList options_list_arg,
                                       const SymbolTable &symbol_table_arg) :
  BasicPriorStatement{move(name_arg1), move(subsample_name_arg), prior_shape_arg, variance_arg, move(options_list_arg)},
  name1(move(name_arg2)),
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
  output << "{\"statementName\": \"corr_prior\""
         << ", \"name2\": \"" << name1 << "\"";
  writeJsonPriorOutput(output);
  output << "}";
}

void
CorrPriorStatement::writeCOutput(ostream &output, const string &basename)
{
  output << endl
         << "index = ";
  if (is_structural_innovation(symbol_table.getType(name)))
    output << "exo_names";
  else
    output << "endo_names";
  output << "[\""<< name << "\"];" << endl;

  output << "index1 = ";
  if (is_structural_innovation(symbol_table.getType(name1)))
    output << "exo_names";
  else
    output << "endo_names";
  output << "[\""<< name1 << "\"];" << endl;

  writeCShape(output);
  writeCOutputHelper(output, "mean");
  writeCOutputHelper(output, "mode");
  writeCOutputHelper(output, "stdev");
  writeCVarianceOption(output);
  writeCDomain(output);

  if (is_structural_innovation(symbol_table.getType(name)))
    output << "msdsgeinfo->addStructuralInnovationCorrPrior(new ModFileStructuralInnovationCorrPrior(";
  else
    output << "msdsgeinfo->addMeasurementErrorCorrPrior(new ModFileMeasurementErrorCorrPrior(";
  output << endl <<"     index, index1, shape, mean, mode, stdev, variance, domain));" << endl;
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
  to_declaration_type(move(to_declaration_type_arg)),
  to_name1(move(to_name1_arg)),
  to_name2(move(to_name2_arg)),
  to_subsample_name(move(to_subsample_name_arg)),
  from_declaration_type(move(from_declaration_type_arg)),
  from_name1(move(from_name1_arg)),
  from_name2(move(from_name2_arg)),
  from_subsample_name(move(from_subsample_name_arg)),
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
  output << "{\"statementName\": \"prior_equal\""
         << ", \"to_name1\": \"" << to_name1 << "\"";
  if (to_declaration_type == "corr")
    output << ", \"to_name2\": \"" << to_name2 << "\"";
  output << ", \"to_subsample\": \"" << to_subsample_name << "\""
         << ", \"from_name1\": \"" << from_name1 << "\"";
  if (to_declaration_type == "corr")
    output << ", \"from_name2\": \"" << from_name2 << "\"";
  output << ", \"from_subsample\": \"" << from_subsample_name << "\""
         << "}";
}

BasicOptionsStatement::~BasicOptionsStatement()
= default;

BasicOptionsStatement::BasicOptionsStatement(string name_arg,
                                             string subsample_name_arg,
                                             OptionsList options_list_arg) :
  name(move(name_arg)),
  subsample_name(move(subsample_name_arg)),
  options_list(move(options_list_arg))
{
}

void
BasicOptionsStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
}

bool
BasicOptionsStatement::is_structural_innovation(const SymbolType symb_type) const
{
  if (symb_type == SymbolType::exogenous || symb_type == SymbolType::exogenousDet)
    return true;
  return false;
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
  auto itn = options_list.num_options.find(field);
  if (itn != options_list.num_options.end())
    output << lhs_field << "." << field << " = " << itn->second << ";" << endl;
}

void
BasicOptionsStatement::writeCOutputHelper(ostream &output, const string &field) const
{
  auto itn = options_list.num_options.find(field);
  if (itn != options_list.num_options.end())
    output << field << " = " << itn->second << ";" << endl;
  else
    output << field << " = " << "numeric_limits<double>::quiet_NaN();" << endl;
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
  output << ", \"name\": \"" << name << "\"";
  if (!subsample_name.empty())
    output << ", \"subsample_name\": \"" << subsample_name << "\"";
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
  output << "{\"statementName\": \"options\"";
  writeJsonOptionsOutput(output);
  output << "}";
}

void
OptionsStatement::writeCOutput(ostream &output, const string &basename)
{
  output << endl
         << "index = param_names[\""<< name << "\"];" << endl;
  writeCOutputHelper(output, "init");
  output << "msdsgeinfo->addOption(new ModFileOption(index, init));" << endl;
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
  output << "{\"statementName\": \"std_options\"";
  writeJsonOptionsOutput(output);
  output << "}";
}

void
StdOptionsStatement::writeCOutput(ostream &output, const string &basename)
{
  output << endl
         << "index = ";
  if (is_structural_innovation(symbol_table.getType(name)))
    output << "exo_names";
  else
    output << "endo_names";
  output << "[\""<< name << "\"];" << endl;

  writeCOutputHelper(output, "init");

  if (is_structural_innovation(symbol_table.getType(name)))
    output << "msdsgeinfo->addStructuralInnovationOption(new ModFileStructuralInnovationOption(";
  else
    output << "msdsgeinfo->addMeasurementErrorOption(new ModFileMeasurementErrorOption(";
  output << "index, init));" << endl;
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
  output << "{\"statementName\": \"corr_options\""
         << ", \"name2\": \"" << name1 << "\"";
  writeJsonOptionsOutput(output);
  output << "}";
}

void
CorrOptionsStatement::writeCOutput(ostream &output, const string &basename)
{
  output << endl
         << "index = ";
  if (is_structural_innovation(symbol_table.getType(name)))
    output << "exo_names";
  else
    output << "endo_names";
  output << "[\""<< name << "\"];" << endl;

  output << "index1 = ";
  if (is_structural_innovation(symbol_table.getType(name1)))
    output << "exo_names";
  else
    output << "endo_names";
  output << "[\""<< name1 << "\"];" << endl;

  writeCOutputHelper(output, "init");

  if (is_structural_innovation(symbol_table.getType(name)))
    output << "msdsgeinfo->addStructuralInnovationCorrOption(new ModFileStructuralInnovationCorrOption(";
  else
    output << "msdsgeinfo->addMeasurementErrorCorrOption(new ModFileMeasurementErrorCorrOption(";
  output << "index, index1, init));" << endl;
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
  to_declaration_type(move(to_declaration_type_arg)),
  to_name1(move(to_name1_arg)),
  to_name2(move(to_name2_arg)),
  to_subsample_name(move(to_subsample_name_arg)),
  from_declaration_type(move(from_declaration_type_arg)),
  from_name1(move(from_name1_arg)),
  from_name2(move(from_name2_arg)),
  from_subsample_name(move(from_subsample_name_arg)),
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
  output << "{\"statementName\": \"options_equal\""
         << ", \"to_name1\": \"" << to_name1 << "\"";
  if (to_declaration_type == "corr")
    output << ", \"to_name2\": \"" << to_name2 << "\"";
  output << ", \"to_subsample\": \"" << to_subsample_name << "\""
         << ", \"from_name1\": \"" << from_name1 << "\"";
  if (to_declaration_type == "corr")
    output << ", \"from_name2\": \"" << from_name2 << "\"";
  output << ", \"from_subsample\": \"" << from_subsample_name << "\""
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
  : symbol_list(move(symbol_list_arg)), options_list(move(options_list_arg))
{
}

void
CalibSmootherStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.calib_smoother_present = true;
}

void
CalibSmootherStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  options_list.writeOutput(output);
  auto it = options_list.string_options.find("parameter_set");
  if (it == options_list.string_options.end())
    output << "options_.parameter_set = 'calibration';" << endl;
  symbol_list.writeOutput("var_list_", output);
  output << "options_.smoother = 1;" << endl
         << "options_.order = 1;" << endl
         << "[oo_, M_, options_, bayestopt_] = evaluate_smoother(options_.parameter_set, var_list_, M_, oo_, options_, bayestopt_, estim_params_);" << endl;
}

void
CalibSmootherStatement::writeJsonOutput(ostream &output) const
{
  output << "{\"statementName\": \"calib_smoother\"";
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
  : options_list(move(options_list_arg))
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
  auto it = options_list.num_options.find("occbin");
  if (it != options_list.string_options.end())
    mod_file_struct.occbin_option = true;
}

void
ExtendedPathStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  // Beware: options do not have the same name in the interface and in the M code...

  for (const auto & num_option : options_list.num_options)
    if (num_option.first != string("periods"))
      output << "options_." << num_option.first << " = " << num_option.second << ";" << endl;

  output << "extended_path([], " << options_list.num_options.find("periods")->second
         << ", [], options_, M_, oo_);" << endl;
}

void
ExtendedPathStatement::writeJsonOutput(ostream &output) const
{
  output << "{\"statementName\": \"extended_path\"";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

ModelDiagnosticsStatement::ModelDiagnosticsStatement()
= default;

void
ModelDiagnosticsStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "model_diagnostics(M_,options_,oo_);" << endl;
}

void
ModelDiagnosticsStatement::writeJsonOutput(ostream &output) const
{
  output << "{\"statementName\": \"model_diagnostics\"}";
}

Smoother2histvalStatement::Smoother2histvalStatement(OptionsList options_list_arg) :
  options_list(move(options_list_arg))
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
  output << "{\"statementName\": \"smoother_2_histval\"";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }
  output << "}";
}

GMMEstimationStatement::GMMEstimationStatement(SymbolList symbol_list_arg,
                                               OptionsList options_list_arg) :
  symbol_list(move(symbol_list_arg)),
  options_list(move(options_list_arg))
{
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
  output << "{\"statementName\": \"gmm_estimation\"";
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
  symbol_list(move(symbol_list_arg)),
  options_list(move(options_list_arg))
{
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
  output << "{\"statementName\": \"smm_estimation\"";
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
  options_list(move(options_list_arg)),
  generate_irf_names(move(generate_irf_names_arg)),
  generate_irf_elements(move(generate_irf_elements_arg))
{
}

void
GenerateIRFsStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  options_list.writeOutput(output);

  if (generate_irf_names.empty())
    return;

  output << "options_.irf_opt.irf_shock_graphtitles = { ";
  for (const auto & generate_irf_name : generate_irf_names)
    output << "'" << generate_irf_name << "'; ";
  output << "};" << endl;

  output << "options_.irf_opt.irf_shocks = zeros(M_.exo_nbr, "
         << generate_irf_names.size() << ");" << endl;

  for (size_t i = 0; i < generate_irf_names.size(); i++)
    {
      map<string, double> m = generate_irf_elements[i];
      for (map<string, double>::const_iterator it = m.begin();
         it != m.end(); it++)
        output << "options_.irf_opt.irf_shocks(M_.exo_names == '"
               << it->first << "', " << i + 1 << ") = "
               << it->second << ";" << endl;
    }
}

void
GenerateIRFsStatement::writeJsonOutput(ostream &output) const
{
  output << "{\"statementName\": \"generate_irfs\"";
  if (options_list.getNumberOfOptions())
    {
      output << ", ";
      options_list.writeJsonOutput(output);
    }

  if (!generate_irf_names.empty())
    {
      output << ", \"irf_elements\": [";
      for (size_t i = 0; i < generate_irf_names.size(); i++)
        {
          output << "{\"name\": \"" << generate_irf_names[i] << "\", \"shocks\": [";
          map<string, double> m = generate_irf_elements[i];
          size_t idx = 0;
          for (map<string, double>::const_iterator it = m.begin();
               it != m.end(); it++, idx++)
            {
              output << "{\"exogenous_variable\": \"" << it->first << "\", "
                     << "\"exogenous_variable_value\": \"" << it->second << "\"}";
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

VarExpectationModelStatement::VarExpectationModelStatement(string model_name_arg, string variable_arg, string var_model_name_arg,
                                                           string horizon_arg, expr_t discount_arg, const SymbolTable &symbol_table_arg) :
  model_name{move(model_name_arg)}, variable{move(variable_arg)},
  var_model_name{move(var_model_name_arg)}, horizon{move(horizon_arg)},
  discount{discount_arg}, symbol_table{symbol_table_arg}
{
}

void
VarExpectationModelStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  string mstruct = "M_.var_expectation." + model_name;
  output << mstruct << ".var_model_name = '" << var_model_name << "';" << endl
         << mstruct << ".horizon = " << horizon << ';' << endl
         << mstruct << ".variable = '" << variable << "';" << endl
         << mstruct << ".variable_id = " << symbol_table.getTypeSpecificID(variable)+1 << ";" << endl;
  auto disc_var = dynamic_cast<const VariableNode *>(discount);
  if (disc_var)
    output << mstruct << ".discount_index = " << symbol_table.getTypeSpecificID(disc_var->get_symb_id()) + 1 << ';' << endl;
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
  output << "{\"statementName\": \"var_expectation_model\","
         << "\"model_name\": \"" << model_name << "\", "
         << "\"variable\": \"" << variable << "\", "
         << "\"var_model_name\": \"" << var_model_name << "\", "
         << "\"horizon\": \"" << horizon << "\", "
         << "\"discount\": \"";
  discount->writeOutput(output);
  output << "\"}";
}
