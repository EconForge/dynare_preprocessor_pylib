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
#include "DynamicModel.hh"

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
                << "            if nargout<3" << endl
                << "                a0 = eye(" << lhs.size() << ");" << endl
                << "            end" << endl
                << "        end" << endl
                << "        if nargout>2" << endl
                << "            constants = zeros(" << lhs.size() << ",1);" << endl;
      for (auto [eqn, expr] : constants.at(name))
        {
          ar_output << "            constants(" << eqn + 1 << ") = ";
          expr->writeOutput(ar_output, ExprNodeOutputType::matlabDynamicModel);
          ar_output << ";" << endl;
        }
      ar_output << "        end" << endl
                << "        if reducedform" << endl
                << "            constants = a0\\constants;" << endl
                << "            a0 = eye(" << lhs.size() << ");" << endl
                << "        end" << endl
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

PacModelTable::PacModelTable(SymbolTable &symbol_table_arg) :
  symbol_table{symbol_table_arg}
{
}

void
PacModelTable::addPacModel(string name_arg, string aux_model_name_arg, string discount_arg, expr_t growth_arg)
{
  if (isExistingPacModelName(name_arg))
    {
      cerr << "Error: a PAC model already exists with the name " << name_arg << endl;
      exit(EXIT_FAILURE);
    }

  aux_model_name[name_arg] = move(aux_model_name_arg);
  discount[name_arg] = move(discount_arg);
  growth[name_arg] = growth_arg;
  original_growth[name_arg] = growth_arg;
  names.insert(move(name_arg));
}

bool
PacModelTable::isExistingPacModelName(const string &name_arg) const
{
  return names.find(name_arg) != names.end();
}

bool
PacModelTable::empty() const
{
  return names.empty();
}

void
PacModelTable::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  for (auto &[name, gv] : growth)
    if (gv)
      gv->collectVariables(SymbolType::exogenous, mod_file_struct.pac_params);
}

void
PacModelTable::findDiffNodesInGrowth(lag_equivalence_table_t &diff_nodes) const
{
  for (auto &[name, gv] : growth)
    if (gv)
      gv->findDiffNodes(diff_nodes);
}

void
PacModelTable::substituteDiffNodesInGrowth(const lag_equivalence_table_t &diff_nodes, ExprNode::subst_table_t &diff_subst_table, vector<BinaryOpNode *> &neweqs)
{
  for (auto &[name, gv] : growth)
    if (gv)
      gv = gv->substituteDiff(diff_nodes, diff_subst_table, neweqs);
}

void
PacModelTable::transformPass(ExprNode::subst_table_t &diff_subst_table,
                             DynamicModel &dynamic_model, const VarModelTable &var_model_table,
                             const TrendComponentModelTable &trend_component_model_table)
{
  // model name → expression for pac_expectation
  map<string, expr_t> pac_expectation_substitution;

  for (const auto &name : names)
    {
      /* Fill the growth_info structure.
         Cannot be done in an earlier pass since growth terms can be
         transformed by DynamicModel::substituteDiff(). */
      if (growth[name])
        try
          {
            growth_info[name] = growth[name]->matchLinearCombinationOfVariables(false);
          }
        catch (ExprNode::MatchFailureException &e)
          {
            auto gv = dynamic_cast<const VariableNode *>(growth[name]);
            if (gv)
              growth_info[name].emplace_back(gv->symb_id, gv->lag, -1, 1);
            else
              {
                cerr << "Pac growth must be a linear combination of variables" << endl;
                exit(EXIT_FAILURE);
              }
          }

      // Collect some information about PAC models
      int max_lag;
      vector<bool> nonstationary;
      if (trend_component_model_table.isExistingTrendComponentModelName(aux_model_name[name]))
        {
          aux_model_type[name] = "trend_component";
          max_lag = trend_component_model_table.getMaxLag(aux_model_name[name]) + 1;
          lhs[name] = dynamic_model.getUndiffLHSForPac(aux_model_name[name], diff_subst_table);
          // All lhs variables in a trend component model are nonstationary
          nonstationary.insert(nonstationary.end(), trend_component_model_table.getDiff(aux_model_name[name]).size(), true);
        }
      else if (var_model_table.isExistingVarModelName(aux_model_name[name]))
        {
          aux_model_type[name] = "var";
          max_lag = var_model_table.getMaxLag(aux_model_name[name]);
          lhs[name] = var_model_table.getLhs(aux_model_name[name]);
          // nonstationary variables in a VAR are those that are in diff
          nonstationary = var_model_table.getDiff(aux_model_name[name]);
        }
      else if (aux_model_name[name].empty())
        max_lag = 0;
      else
        {
          cerr << "Error: aux_model_name not recognized as VAR model or Trend Component model" << endl;
          exit(EXIT_FAILURE);
        }
      dynamic_model.analyzePacEquationStructure(name, eq_name, equation_info);

      // Declare parameter for growth neutrality correction
      if (growth[name])
        {
          string param_name = name + "_pac_growth_neutrality_correction";
          try
            {
              int param_idx = symbol_table.addSymbol(param_name, SymbolType::parameter);
              growth_neutrality_params[name] = param_idx;
            }
          catch (SymbolTable::AlreadyDeclaredException)
            {
              cerr << "The variable/parameter '" << param_name << "' conflicts with the auxiliary parameter that will be generated for the growth neutrality correction of the '" << name << "' PAC model. Please rename that parameter." << endl;
              exit(EXIT_FAILURE);
            }
        }

      // In the MCE case, add the variable and the equation defining Z₁
      if (aux_model_name[name].empty())
        dynamic_model.addPacModelConsistentExpectationEquation(name, symbol_table.getID(discount[name]),
                                                               pacEquationMaxLag(name),
                                                               diff_subst_table,
                                                               mce_z1_symb_ids, mce_alpha_symb_ids);

      // Compute the expressions that will be substituted for the pac_expectation operators
      expr_t growth_correction_term = dynamic_model.Zero;
      if (growth[name])
        growth_correction_term = dynamic_model.AddTimes(growth[name], dynamic_model.AddVariable(growth_neutrality_params[name]));
      if (aux_model_name[name].empty())
        dynamic_model.computePacModelConsistentExpectationSubstitution(name,
                                                                       growth_correction_term,
                                                                       mce_z1_symb_ids[name],
                                                                       pac_expectation_substitution);
      else
        dynamic_model.computePacBackwardExpectationSubstitution(name, lhs[name], max_lag,
                                                                aux_model_type[name],
                                                                nonstationary,
                                                                growth_correction_term,
                                                                h0_indices, h1_indices,
                                                                pac_expectation_substitution);
    }

  // Actually perform the substitution of pac_expectation
  dynamic_model.substitutePacExpectation(pac_expectation_substitution, eq_name);
  dynamic_model.checkNoRemainingPacExpectation();
}

void
PacModelTable::writeOutput(const string &basename, ostream &output) const
{
  for (const auto &name : names)
    {
      output << "M_.pac." << name << ".auxiliary_model_name = '" << aux_model_name.at(name) << "';" << endl
             << "M_.pac." << name << ".discount_index = " << symbol_table.getTypeSpecificID(discount.at(name)) + 1 << ";" << endl;

      if (growth.at(name))
        {
          output << "M_.pac." << name << ".growth_str = '";
          original_growth.at(name)->writeJsonOutput(output, {}, {}, true);
          output << "';" << endl;
          int i = 0;
          for (auto [growth_symb_id, growth_lag, param_id, constant] : growth_info.at(name))
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

  // Write PAC Model Consistent Expectation parameter info
  for (auto &[name, ids] : mce_alpha_symb_ids)
    {
      output << "M_.pac." << name << ".mce.alpha = [";
      for (auto id : ids)
        output << symbol_table.getTypeSpecificID(id) + 1 << " ";
      output << "];" << endl;
    }

  // Write PAC Model Consistent Expectation Z1 info
  for (auto &[name, id] : mce_z1_symb_ids)
    output << "M_.pac." << name << ".mce.z1 = "
           << symbol_table.getTypeSpecificID(id) + 1 << ";" << endl;

  // Write PAC equation name info
  for (auto &[name, eq] : eq_name)
    output << "M_.pac." << name << ".eq_name = '" << eq << "';" << endl;

  for (auto &[model, growth_neutrality_param_index] : growth_neutrality_params)
    output << "M_.pac." << model << ".growth_neutrality_param_index = "
           << symbol_table.getTypeSpecificID(growth_neutrality_param_index) + 1 << ";" << endl;

  for (auto &[model, lhs] : lhs)
    {
      output << "M_.pac." << model << ".lhs = [";
      for (auto id : lhs)
        output << id + 1 << " ";
      output << "];" << endl;
    }

  for (auto &[model, type] : aux_model_type)
      output << "M_.pac." << model << ".auxiliary_model_type = '" << type << "';" << endl;

  for (auto &[name, val] : equation_info)
    {
      auto [lhs_pac_var, optim_share_index, ar_params_and_vars, ec_params_and_vars, non_optim_vars_params_and_constants, additive_vars_params_and_constants, optim_additive_vars_params_and_constants] = val;
      output << "M_.pac." << name << ".lhs_var = "
             << symbol_table.getTypeSpecificID(lhs_pac_var.first) + 1 << ";" << endl;

      if (optim_share_index >= 0)
        output << "M_.pac." << name << ".share_of_optimizing_agents_index = "
               << symbol_table.getTypeSpecificID(optim_share_index) + 1 << ";" << endl;

      output << "M_.pac." << name << ".ec.params = "
             << symbol_table.getTypeSpecificID(ec_params_and_vars.first) + 1 << ";" << endl
             << "M_.pac." << name << ".ec.vars = [";
      for (auto it : ec_params_and_vars.second)
        output << symbol_table.getTypeSpecificID(get<0>(it)) + 1 << " ";
      output << "];" << endl
             << "M_.pac." << name << ".ec.istarget = [";
      for (auto it : ec_params_and_vars.second)
        output << (get<1>(it) ? "true " : "false ");
      output << "];" << endl
             << "M_.pac." << name << ".ec.scale = [";
      for (auto it : ec_params_and_vars.second)
        output << get<2>(it) << " ";
      output << "];" << endl
             << "M_.pac." << name << ".ec.isendo = [";
      for (auto it : ec_params_and_vars.second)
        switch (symbol_table.getType(get<0>(it)))
          {
          case SymbolType::endogenous:
            output << "true ";
            break;
          case SymbolType::exogenous:
            output << "false ";
            break;
          default:
            cerr << "expecting endogenous or exogenous" << endl;
            exit(EXIT_FAILURE);
          }
      output << "];" << endl
             << "M_.pac." << name << ".ar.params = [";
      for (auto &[pid, vid, vlag] : ar_params_and_vars)
        output << (pid != -1 ? symbol_table.getTypeSpecificID(pid) + 1 : -1) << " ";
      output << "];" << endl
             << "M_.pac." << name << ".ar.vars = [";
      for (auto &[pid, vid, vlag] : ar_params_and_vars)
        output << (vid != -1 ? symbol_table.getTypeSpecificID(vid) + 1 : -1) << " ";
      output << "];" << endl
             << "M_.pac." << name << ".ar.lags = [";
      for (auto &[pid, vid, vlag] : ar_params_and_vars)
        output << vlag << " ";
      output << "];" << endl
             << "M_.pac." << name << ".max_lag = " << pacEquationMaxLag(name) << ";" << endl;
      if (!non_optim_vars_params_and_constants.empty())
        {
          output << "M_.pac." << name << ".non_optimizing_behaviour.params = [";
          for (auto &it : non_optim_vars_params_and_constants)
            if (get<2>(it) >= 0)
              output << symbol_table.getTypeSpecificID(get<2>(it)) + 1 << " ";
            else
              output << "NaN ";
          output << "];" << endl
                 << "M_.pac." << name << ".non_optimizing_behaviour.vars = [";
          for (auto &it : non_optim_vars_params_and_constants)
            output << symbol_table.getTypeSpecificID(get<0>(it)) + 1 << " ";
          output << "];" << endl
                 << "M_.pac." << name << ".non_optimizing_behaviour.isendo = [";
          for (auto &it : non_optim_vars_params_and_constants)
            switch (symbol_table.getType(get<0>(it)))
              {
              case SymbolType::endogenous:
                output << "true ";
                break;
              case SymbolType::exogenous:
                output << "false ";
                break;
              default:
                cerr << "expecting endogenous or exogenous" << endl;
                exit(EXIT_FAILURE);
              }
          output << "];" << endl
                 << "M_.pac." << name << ".non_optimizing_behaviour.lags = [";
          for (auto &it : non_optim_vars_params_and_constants)
            output << get<1>(it) << " ";
          output << "];" << endl
                 << "M_.pac." << name << ".non_optimizing_behaviour.scaling_factor = [";
          for (auto &it : non_optim_vars_params_and_constants)
            output << get<3>(it) << " ";
          output << "];" << endl;
        }
      if (!additive_vars_params_and_constants.empty())
        {
          output << "M_.pac." << name << ".additive.params = [";
          for (auto &it : additive_vars_params_and_constants)
            if (get<2>(it) >= 0)
              output << symbol_table.getTypeSpecificID(get<2>(it)) + 1 << " ";
            else
              output << "NaN ";
          output << "];" << endl
                 << "M_.pac." << name << ".additive.vars = [";
          for (auto &it : additive_vars_params_and_constants)
            output << symbol_table.getTypeSpecificID(get<0>(it)) + 1 << " ";
          output << "];" << endl
                 << "M_.pac." << name << ".additive.isendo = [";
          for (auto &it : additive_vars_params_and_constants)
            switch (symbol_table.getType(get<0>(it)))
              {
              case SymbolType::endogenous:
                output << "true ";
                break;
              case SymbolType::exogenous:
                output << "false ";
                break;
              default:
                cerr << "expecting endogenous or exogenous" << endl;
                exit(EXIT_FAILURE);
              }
          output << "];" << endl
                 << "M_.pac." << name << ".additive.lags = [";
          for (auto &it : additive_vars_params_and_constants)
            output << get<1>(it) << " ";
          output << "];" << endl
                 << "M_.pac." << name << ".additive.scaling_factor = [";
          for (auto &it : additive_vars_params_and_constants)
            output << get<3>(it) << " ";
          output << "];" << endl;
        }
      if (!optim_additive_vars_params_and_constants.empty())
        {
          output << "M_.pac." << name << ".optim_additive.params = [";
          for (auto &it : optim_additive_vars_params_and_constants)
            if (get<2>(it) >= 0)
              output << symbol_table.getTypeSpecificID(get<2>(it)) + 1 << " ";
            else
              output << "NaN ";
          output << "];" << endl
                 << "M_.pac." << name << ".optim_additive.vars = [";
          for (auto &it : optim_additive_vars_params_and_constants)
            output << symbol_table.getTypeSpecificID(get<0>(it)) + 1 << " ";
          output << "];" << endl
                 << "M_.pac." << name << ".optim_additive.isendo = [";
          for (auto &it : optim_additive_vars_params_and_constants)
            switch (symbol_table.getType(get<0>(it)))
              {
              case SymbolType::endogenous:
                output << "true ";
                break;
              case SymbolType::exogenous:
                output << "false ";
                break;
              default:
                cerr << "expecting endogenous or exogenous" << endl;
                exit(EXIT_FAILURE);
              }
          output << "];" << endl
                 << "M_.pac." << name << ".optim_additive.lags = [";
          for (auto &it : optim_additive_vars_params_and_constants)
            output << get<1>(it) << " ";
          output << "];" << endl
                 << "M_.pac." << name << ".optim_additive.scaling_factor = [";
          for (auto &it : optim_additive_vars_params_and_constants)
            output << get<3>(it) << " ";
          output << "];" << endl;
        }
      // Create empty h0 and h1 substructures that will be overwritten later if not empty
      output << "M_.pac." << name << ".h0_param_indices = [];" << endl
             << "M_.pac." << name << ".h1_param_indices = [];" << endl;
    }

  for (auto &[name, symb_ids] : h0_indices)
    {
      output << "M_.pac." << name << ".h0_param_indices = [";
      for (auto it : symb_ids)
        output << symbol_table.getTypeSpecificID(it) + 1 << " ";
      output << "];" << endl;
    }

  for (auto &[name, symb_ids] : h1_indices)
    {
      output << "M_.pac." << name << ".h1_param_indices = [";
      for (auto it : symb_ids)
        output << symbol_table.getTypeSpecificID(it) + 1 << " ";
      output << "];" << endl;
    }
}

void
PacModelTable::writeJsonOutput(ostream &output) const
{
  for (const auto &name : names)
    {
      if (name != *names.begin())
        output << ", ";
      output << R"({"statementName": "pac_model",)"
             << R"("model_name": ")" << name << R"(",)"
             << R"("auxiliary_model_name": ")" << aux_model_name.at(name) << R"(",)"
             << R"("discount_index": )" << symbol_table.getTypeSpecificID(discount.at(name)) + 1;
      if (growth.at(name))
        {
          output << R"(,"growth_str": ")";
          original_growth.at(name)->writeJsonOutput(output, {}, {}, true);
          output << R"(")";
        }
      output << "}" << endl;
    }
}

int
PacModelTable::pacEquationMaxLag(const string &name_arg) const
{
  return get<2>(equation_info.at(name_arg)).size();
}
