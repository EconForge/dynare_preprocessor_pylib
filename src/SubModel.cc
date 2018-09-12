/*
 * Copyright (C) 2018 Dynare Team
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

#include <algorithm>

#include "SubModel.hh"

TrendComponentModelTable::TrendComponentModelTable(SymbolTable &symbol_table_arg) :
  symbol_table(symbol_table_arg)
{
}

void
TrendComponentModelTable::addTrendComponentModel(string name_arg,
                                                 vector<string> eqtags_arg,
                                                 vector<string> trend_eqtags_arg)
{
  if (isExistingTrendComponentModelName(name_arg))
    {
      cerr << "Error: a trend component model already exists with the name " << name_arg << endl;
      exit(EXIT_FAILURE);
    }
  eqtags[name_arg] = move(eqtags_arg);
  trend_eqtags[name_arg] = move(trend_eqtags_arg);
  names.insert(move(name_arg));
}

void
TrendComponentModelTable::setVals(map<string, vector<int>> eqnums_arg, map<string, vector<int>> trend_eqnums_arg,
                                  map<string, vector<int>> lhs_arg,
                                  map<string, vector<expr_t>> lhs_expr_t_arg, map<string, vector<bool>> nonstationary_arg)
{
  eqnums = move(eqnums_arg);
  trend_eqnums = move(trend_eqnums_arg);
  lhs = move(lhs_arg);
  lhs_expr_t = move(lhs_expr_t_arg);
  nonstationary = move(nonstationary_arg);

  for (const auto &it : eqnums)
    {
      vector<int> nontrend_vec;
      for (auto eq : it.second)
        if (find(trend_eqnums[it.first].begin(), trend_eqnums[it.first].end(), eq) == trend_eqnums[it.first].end())
          nontrend_vec.push_back(eq);
      nontrend_eqnums[it.first] = nontrend_vec;
    }

  for (const auto &name : names)
    {
      vector<int> nontrend_lhs_vec, trend_lhs_vec;
      vector<int> lhsv = getLhs(name);
      vector<int> eqnumsv = getEqNums(name);
      for (int nontrend_it : getNonTrendEqNums(name))
        nontrend_lhs_vec.push_back(lhsv.at(distance(eqnumsv.begin(), find(eqnumsv.begin(), eqnumsv.end(), nontrend_it))));
      nontrend_lhs[name] = nontrend_lhs_vec;

      for (int trend_it : getTrendEqNums(name))
        trend_lhs_vec.push_back(lhsv.at(distance(eqnumsv.begin(), find(eqnumsv.begin(), eqnumsv.end(), trend_it))));
      trend_lhs[name] = trend_lhs_vec;
    }
}

void
TrendComponentModelTable::setRhs(map<string, vector<set<pair<int, int>>>> rhs_arg)
{
  rhs = move(rhs_arg);
}

void
TrendComponentModelTable::setTrendVar(map<string, vector<int>> trend_vars_arg)
{
  trend_vars = move(trend_vars_arg);
}

void
TrendComponentModelTable::setMaxLags(map<string, vector<int>> max_lags_arg)
{
  max_lags = move(max_lags_arg);
}

void
TrendComponentModelTable::setDiff(map<string, vector<bool>> diff_arg)
{
  diff = move(diff_arg);
}

void
TrendComponentModelTable::setOrigDiffVar(map<string, vector<int>> orig_diff_var_arg)
{
  orig_diff_var = move(orig_diff_var_arg);
}

void
TrendComponentModelTable::setAR(map<string, map<tuple<int, int, int>, expr_t>> AR_arg)
{
  AR = move(AR_arg);
}

void
TrendComponentModelTable::setEC(map<string, map<tuple<int, int, int>, expr_t>> EC_arg)
{
  for (const auto & itmap1 : EC_arg)
    for (const auto &itmap2 : itmap1.second)
      if (get<1>(itmap2.first) != 1)
        {
          cerr << "Error in EC component: lag must be equal to 1" << endl;
          exit(EXIT_FAILURE);
        }
  EC = move(EC_arg);
}

map<string, vector<string>>
TrendComponentModelTable::getEqTags() const
{
  return eqtags;
}

vector<string>
TrendComponentModelTable::getEqTags(const string &name_arg) const
{
  checkModelName(name_arg);
  return eqtags.find(name_arg)->second;
}

void
TrendComponentModelTable::checkModelName(const string &name_arg) const
{
  if (!isExistingTrendComponentModelName(name_arg))
    {
      cerr << name_arg
           << " is not a recognized equation tag of a trend component model equation" << endl;
      exit(EXIT_FAILURE);
    }
}

vector<bool>
TrendComponentModelTable::getNonstationary(const string &name_arg) const
{
  checkModelName(name_arg);
  return nonstationary.find(name_arg)->second;
}

vector<int>
TrendComponentModelTable::getNontrendLhs(const string &name_arg) const
{
  checkModelName(name_arg);
  return nontrend_lhs.find(name_arg)->second;
}

vector<int>
TrendComponentModelTable::getTrendLhs(const string &name_arg) const
{
  checkModelName(name_arg);
  return trend_lhs.find(name_arg)->second;
}

vector<int>
TrendComponentModelTable::getLhs(const string &name_arg) const
{
  checkModelName(name_arg);
  return lhs.find(name_arg)->second;
}

vector<expr_t>
TrendComponentModelTable::getLhsExprT(const string &name_arg) const
{
  checkModelName(name_arg);
  return lhs_expr_t.find(name_arg)->second;
}

map<string, vector<string>>
TrendComponentModelTable::getTrendEqTags() const
{
  return trend_eqtags;
}

map<string, vector<int>>
TrendComponentModelTable::getEqNums() const
{
  return eqnums;
}

map<string, vector<int>>
TrendComponentModelTable::getTrendEqNums() const
{
  return trend_eqnums;
}

vector<int>
TrendComponentModelTable::getTrendEqNums(const string &name_arg) const
{
  checkModelName(name_arg);
  return trend_eqnums.find(name_arg)->second;
}

map<string, vector<int>>
TrendComponentModelTable::getNonTrendEqNums() const
{
  return nontrend_eqnums;
}

vector<int>
TrendComponentModelTable::getNonTrendEqNums(const string &name_arg) const
{
  checkModelName(name_arg);
  return nontrend_eqnums.find(name_arg)->second;
}

vector<int>
TrendComponentModelTable::getEqNums(const string &name_arg) const
{
  checkModelName(name_arg);
  return eqnums.find(name_arg)->second;
}

vector<int>
TrendComponentModelTable::getMaxLags(const string &name_arg) const
{
  checkModelName(name_arg);
  return max_lags.find(name_arg)->second;
}

int
TrendComponentModelTable::getMaxLag(const string &name_arg) const
{
  int max_lag_int = 0;
  for (auto it : getMaxLags(name_arg))
    max_lag_int = max(max_lag_int, it);
  return max_lag_int;
}

vector<bool>
TrendComponentModelTable::getDiff(const string &name_arg) const
{
  checkModelName(name_arg);
  return diff.find(name_arg)->second;
}

vector<int>
TrendComponentModelTable::getOrigDiffVar(const string &name_arg) const
{
  checkModelName(name_arg);
  return orig_diff_var.find(name_arg)->second;
}

void
TrendComponentModelTable::writeOutput(const string &basename, ostream &output) const
{
  string filename = "+" + basename + "/trend_component_ar_ec.m";
  ofstream ar_ec_output;
  ar_ec_output.open(filename, ios::out | ios::binary);
  if (!ar_ec_output.is_open())
    {
      cerr << "Error: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }
  ar_ec_output << "function [ar, ec] = trend_component_ar_ec(model_name, params)" << endl
            << "%function [ar, ec] = trend_component_ar_ec(model_name, params)" << endl
            << "% File automatically generated by the Dynare preprocessor" << endl << endl;

  for (const auto &name : names)
    {
      output << "M_.trend_component." << name << ".model_name = '" << name << "';" << endl
             << "M_.trend_component." << name << ".eqtags = {";
      for (const auto &it : eqtags.at(name))
        output << "'" << it << "'; ";
      output << "};" << endl
             << "M_.trend_component." << name << ".eqn = [";
      for (auto it : eqnums.at(name))
        output << it + 1 << " ";
      output << "];" << endl
             << "M_.trend_component." << name << ".trend_eqn = [";
      for (auto it : trend_eqnums.at(name))
        output << it + 1 << " ";
      output << "];" << endl
             << "M_.trend_component." << name << ".trends = [";
      for (auto it : eqnums.at(name))
        if (find(trend_eqnums.at(name).begin(), trend_eqnums.at(name).end(), it)
            == trend_eqnums.at(name).end())
          output << "false ";
        else
          output << "true ";
      output << "];" << endl
             << "M_.trend_component." << name << ".lhs = [";
      for (auto it : lhs.at(name))
        output << symbol_table.getTypeSpecificID(it) + 1 << " ";
      output << "];" << endl
             << "M_.trend_component." << name << ".max_lag = [";
      for (auto it : max_lags.at(name))
        output << it << " ";
      output << "];" << endl
             << "M_.trend_component." << name << ".diff = [";
      for (const auto &it : diff.at(name))
        output << (it ? "true" : "false") << " ";
      output << "];" << endl
             << "M_.trend_component." << name << ".orig_diff_var = [";
      for (auto it : orig_diff_var.at(name))
        output << (it >= 0 ? symbol_table.getTypeSpecificID(it) + 1 : -1) << " ";
      output << "];" << endl
             << "M_.trend_component." << name << ".nonstationary = [";
      for (auto it : nonstationary.at(name))
        output << (it ? "true" : "false") << " ";
      output << "];" << endl;
      int i = 1;
      for (const auto &it : rhs.at(name))
        {
          output << "M_.trend_component." << name << ".rhs.vars_at_eq{" << i << "}.var = [";
          for (const auto &it1 : it)
            output << symbol_table.getTypeSpecificID(it1.first) + 1 << " ";
          output << "];" << endl
                 << "M_.trend_component." << name << ".rhs.vars_at_eq{" << i << "}.lag = [";
          for (const auto &it1 : it)
            output << it1.second << " ";
          output << "];" << endl;

          i++;
        }
      output << "M_.trend_component." << name << ".trend_vars = [";
      for (auto it : trend_vars.at(name))
        output << (it >= 0 ? symbol_table.getTypeSpecificID(it) + 1 : -1) << " ";
      output << "];" << endl;

      vector<int> trend_lhs_vec = getTrendLhs(name);
      vector<int> nontrend_lhs_vec = getNontrendLhs(name);

      ar_ec_output << "if strcmp(model_name, '" << name << "')" << endl
                << "    % AR" << endl
                << "    ar = zeros(" << nontrend_lhs_vec.size() << ", " << nontrend_lhs_vec.size() << ", " << getMaxLag(name) << ");" << endl;
      for (const auto & it : AR.at(name))
        {
          int eqn, lag, lhs_symb_id;
          tie (eqn, lag, lhs_symb_id) = it.first;
          int colidx = (int) distance(nontrend_lhs_vec.begin(), find(nontrend_lhs_vec.begin(), nontrend_lhs_vec.end(), lhs_symb_id));
          ar_ec_output << "    ar(" << eqn + 1 << ", " << colidx + 1 << ", " << lag << ") = ";
          it.second->writeOutput(ar_ec_output, ExprNodeOutputType::matlabDynamicModel);
          ar_ec_output << ";" << endl;
        }
      ar_ec_output << endl
                   << "    % EC" << endl
                   << "    ec = zeros(" << nontrend_lhs_vec.size() << ", " << nontrend_lhs_vec.size() << ", 1);" << endl;
      for (const auto & it : EC.at(name))
        {
          int eqn, lag, lhs_symb_id;
          tie (eqn, lag, lhs_symb_id) = it.first;
          int colidx = (int) distance(trend_lhs_vec.begin(), find(trend_lhs_vec.begin(), trend_lhs_vec.end(), lhs_symb_id));
          ar_ec_output << "    ec(" << eqn + 1 << ", " << colidx + 1 << ", 1) = ";
          it.second->writeOutput(ar_ec_output, ExprNodeOutputType::matlabDynamicModel);
          ar_ec_output << ";" << endl;
        }
      ar_ec_output << "    return" << endl
                << "end" << endl << endl;
    }
  ar_ec_output << "error([model_name ' is not a valid trend_component_model name'])" << endl
            << "end" << endl;
  ar_ec_output.close();
}

void
TrendComponentModelTable::writeJsonOutput(ostream &output) const
{
  for (const auto &name : names)
    {
      if (name != *(names.begin()))
        output << ", ";
      output << "{\"statementName\": \"trend_component_model\","
             << "\"model_name\": \"" << name << "\","
             << "\"eqtags\": [";
      for (const auto &it : eqtags.at(name))
        {
          output << "\"" << it << "\"";
          if (&it != &eqtags.at(name).back())
            output << ", ";
        }
      output << "], \"trend_eqtags\": [";
      for (const auto &it : trend_eqtags.at(name))
        {
          output << "\"" << it << "\"";
          if (&it != &trend_eqtags.at(name).back())
            output << ", ";
        }
      output << "]}";
    }
}

VarModelTable::VarModelTable(SymbolTable &symbol_table_arg) :
  symbol_table(symbol_table_arg)
{
}

void
VarModelTable::addVarModel(string name_arg, vector<string> eqtags_arg,
                           pair<SymbolList, int> symbol_list_and_order_arg)
{
  if (isExistingVarModelName(name_arg))
    {
      cerr << "Error: a VAR model already exists with the name " << name_arg << endl;
      exit(EXIT_FAILURE);
    }

  eqtags[name_arg] = move(eqtags_arg);
  symbol_list_and_order[name_arg] = move(symbol_list_and_order_arg);
  names.insert(move(name_arg));
}

map<string, pair<SymbolList, int>>
VarModelTable::getSymbolListAndOrder() const
{
  return symbol_list_and_order;
}

void
VarModelTable::writeOutput(const string &basename, ostream &output) const
{
  string filename = "+" + basename + "/var_ar.m";
  ofstream ar_output;
  ar_output.open(filename, ios::out | ios::binary);
  if (!ar_output.is_open())
    {
      cerr << "Error: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }
  ar_output << "function ar = var_ar(model_name, params)" << endl
            << "%function ar = var_ar(model_name, params)" << endl
            << "% File automatically generated by the Dynare preprocessor" << endl << endl;

  for (const auto &name : names)
    {
      output << "M_.var." << name << ".model_name = '" << name << "';" << endl;
      if (!symbol_list_and_order.at(name).first.empty())
        {
          symbol_list_and_order.at(name).first.writeOutput("M_.var." + name + ".var_list_", output);
          output  << "M_.var." << name << ".order = "
                  << symbol_list_and_order.at(name).second << ";" << endl;
        }
      output << "M_.var." << name << ".eqtags = {";
      for (const auto &it : eqtags.at(name))
        output << "'" << it << "'; ";
      output << "};" << endl
             << "M_.var." << name << ".eqn = [";
      for (auto it : eqnums.at(name))
        output << it + 1 << " ";
      output << "];" << endl
             << "M_.var." << name << ".lhs = [";
      for (auto it : lhs.at(name))
        output << symbol_table.getTypeSpecificID(it) + 1 << " ";
      output << "];" << endl
             << "M_.var." << name << ".max_lag = [";
      for (auto it : max_lags.at(name))
        output << it << " ";
      output << "];" << endl
             << "M_.var." << name << ".diff = [";
      for (const auto &it : diff.at(name))
        output << (it ? "true" : "false") << " ";
      output << "];" << endl
             << "M_.var." << name << ".orig_diff_var = [";
      for (auto it : orig_diff_var.at(name))
        output << (it >= 0 ? symbol_table.getTypeSpecificID(it) + 1 : -1) << " ";
      output << "];" << endl
             << "M_.var." << name << ".nonstationary = [";
      for (auto it : nonstationary.at(name))
        output << (it ? "true" : "false") << " ";
      output << "];" << endl;
      int i = 1;
      for (const auto &it : rhs.at(name))
        {
          output << "M_.var." << name << ".rhs.vars_at_eq{" << i << "}.var = [";
          for (const auto &it1 : it)
            output << symbol_table.getTypeSpecificID(it1.first) + 1 << " ";
          output << "];" << endl
                 << "M_.var." << name << ".rhs.vars_at_eq{" << i << "}.lag = [";
          for (const auto &it1 : it)
            output << it1.second << " ";
          output << "];" << endl;

          i++;
        }

      vector<int> lhs = getLhs(name);
      ar_output << "if strcmp(model_name, '" << name << "')" << endl
                << "    ar = zeros(" << lhs.size() << ", " << lhs.size() << ", " << getMaxLag(name) << ");" << endl;
      for (const auto & it : AR.at(name))
        {
          int eqn, lag, lhs_symb_id;
          tie (eqn, lag, lhs_symb_id) = it.first;
          int colidx = (int) distance(lhs.begin(), find(lhs.begin(), lhs.end(), lhs_symb_id));
          ar_output << "    ar(" << eqn + 1 << ", " << colidx + 1 << ", " << lag
                    << ") = ";
          it.second->writeOutput(ar_output, ExprNodeOutputType::matlabDynamicModel);
          ar_output << endl;
        }
      ar_output << "    return" << endl
                << "end" << endl << endl;
    }
  ar_output << "error([model_name ' is not a valid var_model name'])" << endl
            << "end" << endl;
  ar_output.close();
}

void
VarModelTable::writeJsonOutput(ostream &output) const
{
  for (const auto &name : names)
    {
      if (name != *(names.begin()))
        output << ", ";
      output << "{\"statementName\": \"var_model\","
             << "\"model_name\": \"" << name << "\",";
      if (symbol_list_and_order.empty())
        {
          output << "\"eqtags\": [";
          for (const auto &it : eqtags.at(name))
            {
              output << "\"" << it << "\"";
              if (&it != &eqtags.at(name).back())
                output << ", ";
            }
          output << "]";
        }
      else
        {
          output << "\"order\": \"" << symbol_list_and_order.at(name).second << "\",";
        }
      output << "}";
    }
}

map<string, vector<string>>
VarModelTable::getEqTags() const
{
  return eqtags;
}

vector<string>
VarModelTable::getEqTags(const string &name_arg) const
{
  checkModelName(name_arg);
  return eqtags.find(name_arg)->second;
}

void
VarModelTable::checkModelName(const string &name_arg) const
{
  if (!isExistingVarModelName(name_arg))
    {
      cerr << name_arg
           << " is not a recognized equation tag of a VAR model equation" << endl;
      exit(EXIT_FAILURE);
    }
}

void
VarModelTable::setEqNums(map<string, vector<int>> eqnums_arg)
{
  eqnums = move(eqnums_arg);
}

void
VarModelTable::setLhs(map<string, vector<int>> lhs_arg)
{
  lhs = move(lhs_arg);
}

void
VarModelTable::setRhs(map<string, vector<set<pair<int, int>>>> rhs_arg)
{
  rhs = move(rhs_arg);
}

void
VarModelTable::setLhsExprT(map<string, vector<expr_t>> lhs_expr_t_arg)
{
  lhs_expr_t = move(lhs_expr_t_arg);
}

void
VarModelTable::setNonstationary(map<string, vector<bool>> nonstationary_arg)
{
  nonstationary = move(nonstationary_arg);
}

map<string, vector<int>>
VarModelTable::getEqNums() const
{
  return eqnums;
}

vector<int>
VarModelTable::getEqNums(const string &name_arg) const
{
  checkModelName(name_arg);
  return eqnums.find(name_arg)->second;
}

void
VarModelTable::setMaxLags(map<string, vector<int>> max_lags_arg)
{
  max_lags = move(max_lags_arg);
}

void
VarModelTable::setDiff(map<string, vector<bool>> diff_arg)
{
  diff = move(diff_arg);
}

void
VarModelTable::setOrigDiffVar(map<string, vector<int>> orig_diff_var_arg)
{
  orig_diff_var = move(orig_diff_var_arg);
}

void
VarModelTable::setAR(map<string, map<tuple<int, int, int>, expr_t>> AR_arg)
{
  AR = move(AR_arg);
}

vector<int>
VarModelTable::getMaxLags(const string &name_arg) const
{
  checkModelName(name_arg);
  return max_lags.find(name_arg)->second;
}

int
VarModelTable::getMaxLag(const string &name_arg) const
{
  int max_lag_int = 0;
  for (auto it : getMaxLags(name_arg))
    max_lag_int = max(max_lag_int, it);
  return max_lag_int;
}

vector<int>
VarModelTable::getLhs(const string &name_arg) const
{
  checkModelName(name_arg);
  return lhs.find(name_arg)->second;
}

vector<bool>
VarModelTable::getNonstationary(const string &name_arg) const
{
  checkModelName(name_arg);
  return nonstationary.find(name_arg)->second;
}

vector<set<pair<int, int>>>
VarModelTable::getRhs(const string &name_arg) const
{
  checkModelName(name_arg);
  return rhs.find(name_arg)->second;
}


vector<expr_t>
VarModelTable::getLhsExprT(const string &name_arg) const
{
  checkModelName(name_arg);
  return lhs_expr_t.find(name_arg)->second;
}
