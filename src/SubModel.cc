/*
 * Copyright © 2018-2021 Dynare Team
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

#include <algorithm>

#include "SubModel.hh"

TrendComponentModelTable::TrendComponentModelTable(SymbolTable &symbol_table_arg) :
  symbol_table{symbol_table_arg}
{
}

void
TrendComponentModelTable::addTrendComponentModel(string name_arg,
                                                 vector<string> eqtags_arg,
                                                 vector<string> target_eqtags_arg)
{
  if (isExistingTrendComponentModelName(name_arg))
    {
      cerr << "Error: a trend component model already exists with the name " << name_arg << endl;
      exit(EXIT_FAILURE);
    }
  eqtags[name_arg] = move(eqtags_arg);
  target_eqtags[name_arg] = move(target_eqtags_arg);
  names.insert(move(name_arg));
}

void
TrendComponentModelTable::setVals(map<string, vector<int>> eqnums_arg, map<string, vector<int>> target_eqnums_arg,
                                  map<string, vector<int>> lhs_arg, map<string, vector<expr_t>> lhs_expr_t_arg)
{
  eqnums = move(eqnums_arg);
  target_eqnums = move(target_eqnums_arg);
  lhs = move(lhs_arg);
  lhs_expr_t = move(lhs_expr_t_arg);

  for (const auto &it : eqnums)
    {
      vector<int> nontrend_vec;
      for (auto eq : it.second)
        if (find(target_eqnums[it.first].begin(), target_eqnums[it.first].end(), eq) == target_eqnums[it.first].end())
          nontrend_vec.push_back(eq);
      nontarget_eqnums[it.first] = nontrend_vec;
    }

  for (const auto &name : names)
    {
      vector<int> nontarget_lhs_vec, target_lhs_vec;
      vector<int> lhsv = getLhs(name);
      vector<int> eqnumsv = getEqNums(name);
      for (int nontrend_it : getNonTargetEqNums(name))
        nontarget_lhs_vec.push_back(lhsv.at(distance(eqnumsv.begin(), find(eqnumsv.begin(), eqnumsv.end(), nontrend_it))));
      nontarget_lhs[name] = nontarget_lhs_vec;

      for (int trend_it : getTargetEqNums(name))
        target_lhs_vec.push_back(lhsv.at(distance(eqnumsv.begin(), find(eqnumsv.begin(), eqnumsv.end(), trend_it))));
      target_lhs[name] = target_lhs_vec;
    }
}

void
TrendComponentModelTable::setRhs(map<string, vector<set<pair<int, int>>>> rhs_arg)
{
  rhs = move(rhs_arg);
}

void
TrendComponentModelTable::setTargetVar(map<string, vector<int>> target_vars_arg)
{
  target_vars = move(target_vars_arg);
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
TrendComponentModelTable::setA0(map<string, map<tuple<int, int, int>, expr_t>> A0_arg,
                                map<string, map<tuple<int, int, int>, expr_t>> A0star_arg)
{
  A0 = move(A0_arg);
  A0star = move(A0star_arg);
}

const map<string, vector<string>> &
TrendComponentModelTable::getEqTags() const
{
  return eqtags;
}

const vector<string> &
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

const vector<int> &
TrendComponentModelTable::getNonTargetLhs(const string &name_arg) const
{
  checkModelName(name_arg);
  return nontarget_lhs.find(name_arg)->second;
}

const vector<int> &
TrendComponentModelTable::getTargetLhs(const string &name_arg) const
{
  checkModelName(name_arg);
  return target_lhs.find(name_arg)->second;
}

const vector<int> &
TrendComponentModelTable::getLhs(const string &name_arg) const
{
  checkModelName(name_arg);
  return lhs.find(name_arg)->second;
}

const vector<expr_t> &
TrendComponentModelTable::getLhsExprT(const string &name_arg) const
{
  checkModelName(name_arg);
  return lhs_expr_t.find(name_arg)->second;
}

const map<string, vector<string>> &
TrendComponentModelTable::getTargetEqTags() const
{
  return target_eqtags;
}

const map<string, vector<int>> &
TrendComponentModelTable::getEqNums() const
{
  return eqnums;
}

const map<string, vector<int>> &
TrendComponentModelTable::getTargetEqNums() const
{
  return target_eqnums;
}

const vector<int> &
TrendComponentModelTable::getTargetEqNums(const string &name_arg) const
{
  checkModelName(name_arg);
  return target_eqnums.find(name_arg)->second;
}

const map<string, vector<int>> &
TrendComponentModelTable::getNonTargetEqNums() const
{
  return nontarget_eqnums;
}

const vector<int> &
TrendComponentModelTable::getNonTargetEqNums(const string &name_arg) const
{
  checkModelName(name_arg);
  return nontarget_eqnums.find(name_arg)->second;
}

const vector<int> &
TrendComponentModelTable::getEqNums(const string &name_arg) const
{
  checkModelName(name_arg);
  return eqnums.find(name_arg)->second;
}

const vector<int> &
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

const vector<bool> &
TrendComponentModelTable::getDiff(const string &name_arg) const
{
  checkModelName(name_arg);
  return diff.find(name_arg)->second;
}

const vector<int> &
TrendComponentModelTable::getOrigDiffVar(const string &name_arg) const
{
  checkModelName(name_arg);
  return orig_diff_var.find(name_arg)->second;
}

void
TrendComponentModelTable::writeOutput(const string &basename, ostream &output) const
{
  if (names.empty())
    return;

  string filename = "+" + basename + "/trend_component_ar_a0.m";
  ofstream ar_ec_output;
  ar_ec_output.open(filename, ios::out | ios::binary);
  if (!ar_ec_output.is_open())
    {
      cerr << "Error: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }
  ar_ec_output << "function [AR, A0, A0star] = trend_component_ar_a0(model_name, params)" << endl
               << "%function [AR, A0, A0star] = trend_component_ar_a0(model_name, params)" << endl
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
             << "M_.trend_component." << name << ".targets = [";
      for (auto it : eqnums.at(name))
        if (find(target_eqnums.at(name).begin(), target_eqnums.at(name).end(), it)
            == target_eqnums.at(name).end())
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
      for (size_t i = 0; i < diff.at(name).size(); i++)
        output << "true ";
      output << "];" << endl;
      int i = 1;
      for (const auto &it : rhs.at(name))
        {
          output << "M_.trend_component." << name << ".rhs.vars_at_eq{" << i << "}.var = [";
          for (auto [var, lag] : it)
            output << symbol_table.getTypeSpecificID(var) + 1 << " ";
          output << "];" << endl
                 << "M_.trend_component." << name << ".rhs.vars_at_eq{" << i << "}.lag = [";
          for (auto [var, lag] : it)
            output << lag << " ";
          output << "];" << endl;

          i++;
        }
      output << "M_.trend_component." << name << ".target_vars = [";
      for (auto it : target_vars.at(name))
        output << (it >= 0 ? symbol_table.getTypeSpecificID(it) + 1 : -1) << " ";
      output << "];" << endl;

      vector<string> target_eqtags_vec = target_eqtags.at(name);
      output << "M_.trend_component." << name << ".target_eqtags = {";
      for (auto it : target_eqtags_vec)
        output << "'" << it << "';";
      output << "};" << endl;

      vector<string> eqtags_vec = eqtags.at(name);
      output << "M_.trend_component." << name << ".target_eqn = [";
      for (auto it : target_eqtags_vec)
        {
          int i = 0;
          for (auto it1 : eqtags_vec)
            {
              i++;
              if (it == it1)
                {
                  output << i << " ";
                  break;
                }
            }
        }
      output << "];" << endl;

      vector<int> target_lhs_vec = getTargetLhs(name);
      vector<int> nontarget_lhs_vec = getNonTargetLhs(name);

      ar_ec_output << "if strcmp(model_name, '" << name << "')" << endl
                   << "    % AR" << endl
                   << "    AR = zeros(" << nontarget_lhs_vec.size() << ", " << nontarget_lhs_vec.size() << ", " << getMaxLag(name) << ");" << endl;
      for (const auto &[key, expr] : AR.at(name))
        {
          auto [eqn, lag, lhs_symb_id] = key;
          int colidx = static_cast<int>(distance(nontarget_lhs_vec.begin(), find(nontarget_lhs_vec.begin(), nontarget_lhs_vec.end(), lhs_symb_id)));
          ar_ec_output << "    AR(" << eqn + 1 << ", " << colidx + 1 << ", " << lag << ") = ";
          expr->writeOutput(ar_ec_output, ExprNodeOutputType::matlabDynamicModel);
          ar_ec_output << ";" << endl;
        }

      int a0_lag = 0;
      for (const auto &[key, expr] : A0.at(name))
        a0_lag = max(a0_lag, get<1>(key));
      ar_ec_output << endl
                   << "    % A0" << endl
                   << "    A0 = zeros(" << nontarget_lhs_vec.size() << ", " << nontarget_lhs_vec.size() << ", " << a0_lag << ");" << endl;
      for (const auto &[key, expr] : A0.at(name))
        {
          auto [eqn, lag, colidx] = key;
          ar_ec_output << "    A0(" << eqn + 1 << ", " << colidx + 1 << ", " << lag << ") = ";
          expr->writeOutput(ar_ec_output, ExprNodeOutputType::matlabDynamicModel);
          ar_ec_output << ";" << endl;
        }

      int a0star_lag = 0;
      for (const auto &[key, expr] : A0star.at(name))
        a0star_lag = max(a0star_lag, get<1>(key));
      ar_ec_output << endl
                   << "    % A0star" << endl
                   << "    A0star = zeros(" << nontarget_lhs_vec.size() << ", " << target_lhs_vec.size() << ", " << a0star_lag << ");" << endl;
      for (const auto &[key, expr] : A0star.at(name))
        {
          auto [eqn, lag, colidx] = key;
          ar_ec_output << "    A0star(" << eqn + 1 << ", " << colidx + 1 << ", " << lag << ") = ";
          expr->writeOutput(ar_ec_output, ExprNodeOutputType::matlabDynamicModel);
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
      if (name != *names.begin())
        output << ", ";
      output << R"({"statementName": "trend_component_model",)"
             << R"("model_name": ")" << name << R"(",)"
             << R"("eqtags": [)";
      for (const auto &it : eqtags.at(name))
        {
          output << R"(")" << it << R"(")";
          if (&it != &eqtags.at(name).back())
            output << ", ";
        }
      output << R"(], "target_eqtags": [)";
      for (const auto &it : target_eqtags.at(name))
        {
          output << R"(")" << it << R"(")";
          if (&it != &target_eqtags.at(name).back())
            output << ", ";
        }
      output << "]}";
    }
}

VarModelTable::VarModelTable(SymbolTable &symbol_table_arg) :
  symbol_table{symbol_table_arg}
{
}

void
VarModelTable::addVarModel(string name_arg, bool structural_arg, vector<string> eqtags_arg)
{
  if (isExistingVarModelName(name_arg))
    {
      cerr << "Error: a VAR model already exists with the name " << name_arg << endl;
      exit(EXIT_FAILURE);
    }

  structural[name_arg] = structural_arg;
  eqtags[name_arg] = move(eqtags_arg);
  names.insert(move(name_arg));
}

void
VarModelTable::writeOutput(const string &basename, ostream &output) const
{
  if (names.empty())
    return;

  string filename = "+" + basename + "/varmatrices.m";
  ofstream ar_output;
  ar_output.open(filename, ios::out | ios::binary);
  if (!ar_output.is_open())
    {
      cerr << "Error: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }
  ar_output << "function [ar, a0, constants] = varmatrices(model_name, params, reducedform)" << endl
            << "% File automatically generated by the Dynare preprocessor" << endl << endl
            << "if nargin<3" << endl
            << "    reducedform = false;" << endl
            << "end" << endl << endl;

  for (const auto &name : names)
    {
      output << "M_.var." << name << ".model_name = '" << name << "';" << endl
             << "M_.var." << name << ".structural = " << (structural.at(name) ? "true" : "false") << ";" << endl
             << "M_.var." << name << ".eqtags = {";
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
             << "M_.var." << name << ".nonstationary = M_.var." << name << ".diff;" << endl
             << "M_.var." << name << ".orig_diff_var = [";
      for (auto it : orig_diff_var.at(name))
        output << (it >= 0 ? symbol_table.getTypeSpecificID(it) + 1 : -1) << " ";
      output << "];" << endl;
      int i = 1;
      for (const auto &it : rhs.at(name))
        {
          output << "M_.var." << name << ".rhs.vars_at_eq{" << i << "}.var = [";
          for (auto [var, lag] : it)
            output << symbol_table.getTypeSpecificID(var) + 1 << " ";
          output << "];" << endl
                 << "M_.var." << name << ".rhs.vars_at_eq{" << i << "}.lag = [";
          for (auto [var, lag] : it)
            output << lag << " ";
          output << "];" << endl;

          i++;
        }

      vector<int> lhs = getLhsOrigIds(name);
      ar_output << "if strcmp(model_name, '" << name << "')" << endl
                << "    ar = zeros(" << lhs.size() << ", " << lhs.size() << ", " << getMaxLag(name) << ");" << endl;
      for (const auto &[key, expr] : AR.at(name))
        {
          auto [eqn, lag, lhs_symb_id] = key;
          int colidx = static_cast<int>(distance(lhs.begin(), find(lhs.begin(), lhs.end(), lhs_symb_id)));
          ar_output << "    ar(" << eqn + 1 << "," << colidx + 1 << "," << lag << ") = ";
          expr->writeOutput(ar_output, ExprNodeOutputType::matlabDynamicModel);
          ar_output << ";" << endl;
        }
      ar_output << "    if nargout>1" << endl
                << "        a0 = eye(" << lhs.size() << ");" << endl;
      for (const auto &[key, expr] : A0.at(name))
        {
          auto [eqn, lhs_symb_id] = key;
          int colidx = static_cast<int>(distance(lhs.begin(), find(lhs.begin(), lhs.end(), lhs_symb_id)));
          if (eqn!=colidx)
            {
              ar_output << "        a0(" << eqn + 1 << "," << colidx + 1 << ") = ";
              expr->writeOutput(ar_output, ExprNodeOutputType::matlabDynamicModel);
              ar_output << ";" << endl;
            }
        }
      ar_output << "        if reducedform" << endl
                << "            for i=1:" << getMaxLag(name) << endl
                << "                ar(:,:,i) = a0\\ar(:,:,i);" << endl
                << "            end" << endl
                << "            a0 = eye(" << lhs.size() << ");" << endl
                << "        end" << endl
                << "        if nargout>2" << endl
                << "            constants = zeros(" << lhs.size() << ");" << endl;
      for (auto [eqn, expr] : constants.at(name))
        {
          ar_output << "            constants(" << eqn + 1 << ") = ";
          expr->writeOutput(ar_output, ExprNodeOutputType::matlabDynamicModel);
          ar_output << ";" << endl;
        }
      ar_output << "        end" << endl
                << "    end" << endl
                << "    return" << endl
                << "end" << endl << endl;
    }
  ar_output << "error('%s is not a valid var_model name', model_name)" << endl;
  ar_output.close();
}

void
VarModelTable::writeJsonOutput(ostream &output) const
{
  for (const auto &name : names)
    {
      if (name != *names.begin())
        output << ", ";
      output << R"({"statementName": "var_model",)"
             << R"("model_name": ")" << name << R"(",)"
             << R"("eqtags": [)";
      for (const auto &it : eqtags.at(name))
        {
          output << R"(")" << it << R"(")";
          if (&it != &eqtags.at(name).back())
            output << ", ";
        }
      output << "]}";
    }
}

const map<string, bool> &
VarModelTable::getStructural() const
{
  return structural;
}

const map<string, vector<string>> &
VarModelTable::getEqTags() const
{
  return eqtags;
}

const vector<string> &
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
  for (auto it : lhs)
    {
      vector<int> lhsvec;
      for (auto ids : it.second)
        {
          int lhs_last_orig_symb_id = ids;
          int lhs_orig_symb_id = ids;
          if (symbol_table.isAuxiliaryVariable(lhs_orig_symb_id))
            try
              {
                lhs_last_orig_symb_id = lhs_orig_symb_id;
                lhs_orig_symb_id = symbol_table.getOrigSymbIdForAuxVar(lhs_orig_symb_id);
              }
            catch (...)
              {
              }
          lhsvec.emplace_back(lhs_last_orig_symb_id);
        }
      lhs_orig_symb_ids[it.first] = lhsvec;
    }
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

const map<string, vector<int>> &
VarModelTable::getEqNums() const
{
  return eqnums;
}

const vector<bool> &
VarModelTable::getDiff(const string &name_arg) const
{
  checkModelName(name_arg);
  return diff.find(name_arg)->second;
}

const vector<int> &
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

void
VarModelTable::setA0(map<string, map<tuple<int, int>, expr_t>> A0_arg)
{
  A0 = move(A0_arg);
}

void
VarModelTable::setConstants(map<string, map<int, expr_t>> constants_arg)
{
  constants = move(constants_arg);
}

const vector<int> &
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

const vector<int> &
VarModelTable::getLhs(const string &name_arg) const
{
  checkModelName(name_arg);
  return lhs.find(name_arg)->second;
}

const vector<int> &
VarModelTable::getLhsOrigIds(const string &name_arg) const
{
  checkModelName(name_arg);
  return lhs_orig_symb_ids.find(name_arg)->second;
}

const vector<set<pair<int, int>>> &
VarModelTable::getRhs(const string &name_arg) const
{
  checkModelName(name_arg);
  return rhs.find(name_arg)->second;
}

const vector<expr_t> &
VarModelTable::getLhsExprT(const string &name_arg) const
{
  checkModelName(name_arg);
  return lhs_expr_t.find(name_arg)->second;
}
