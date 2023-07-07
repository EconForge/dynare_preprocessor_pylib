/*
 * Copyright © 2018-2023 Dynare Team
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
#include <cassert>
#include <numeric>

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
TrendComponentModelTable::setTargetVar(map<string, vector<optional<int>>> target_vars_arg)
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
TrendComponentModelTable::setOrigDiffVar(map<string, vector<optional<int>>> orig_diff_var_arg)
{
  orig_diff_var = move(orig_diff_var_arg);
}

void
TrendComponentModelTable::setAR(map<string, map<tuple<int, int, int>, expr_t>> AR_arg)
{
  AR = move(AR_arg);
}

void
TrendComponentModelTable::setA0(map<string, map<tuple<int, int>, expr_t>> A0_arg,
                                map<string, map<tuple<int, int>, expr_t>> A0star_arg)
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
  return eqtags.at(name_arg);
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
  return nontarget_lhs.at(name_arg);
}

const vector<int> &
TrendComponentModelTable::getTargetLhs(const string &name_arg) const
{
  checkModelName(name_arg);
  return target_lhs.at(name_arg);
}

const vector<int> &
TrendComponentModelTable::getLhs(const string &name_arg) const
{
  checkModelName(name_arg);
  return lhs.at(name_arg);
}

const vector<expr_t> &
TrendComponentModelTable::getLhsExprT(const string &name_arg) const
{
  checkModelName(name_arg);
  return lhs_expr_t.at(name_arg);
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
  return target_eqnums.at(name_arg);
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
  return nontarget_eqnums.at(name_arg);
}

const vector<int> &
TrendComponentModelTable::getEqNums(const string &name_arg) const
{
  checkModelName(name_arg);
  return eqnums.at(name_arg);
}

const vector<int> &
TrendComponentModelTable::getMaxLags(const string &name_arg) const
{
  checkModelName(name_arg);
  return max_lags.at(name_arg);
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
  return diff.at(name_arg);
}

void
TrendComponentModelTable::writeOutput(const string &basename, ostream &output) const
{
  if (names.empty())
    return;

  const filesystem::path filename {DataTree::packageDir(basename) / "trend_component_ar_a0.m"};
  ofstream ar_ec_output{filename, ios::out | ios::binary};
  if (!ar_ec_output.is_open())
    {
      cerr << "ERROR: Can't open file " << filename.string() << " for writing" << endl;
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
      for (bool it : diff.at(name))
        output << boolalpha << it << " ";
      output << "];" << endl
             << "M_.trend_component." << name << ".orig_diff_var = [";
      for (const auto &it : orig_diff_var.at(name))
        output << (it ? symbol_table.getTypeSpecificID(*it) + 1 : -1) << " ";
      output << "];" << endl
             << "M_.trend_component." << name << ".nonstationary = [";
      for (size_t i = 0; i < diff.at(name).size(); i++)
        output << "true ";
      output << "];" << endl;
      for (int i{1};
           const auto &it : rhs.at(name))
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
      for (const optional<int> &it : target_vars.at(name))
        output << (it ? symbol_table.getTypeSpecificID(*it) + 1 : -1) << " ";
      output << "];" << endl;

      vector<string> target_eqtags_vec = target_eqtags.at(name);
      output << "M_.trend_component." << name << ".target_eqtags = {";
      for (auto it : target_eqtags_vec)
        output << "'" << it << "';";
      output << "};" << endl;

      vector<string> eqtags_vec = eqtags.at(name);
      output << "M_.trend_component." << name << ".target_eqn = [";
      for (auto it : target_eqtags_vec)
        output << distance(eqtags_vec.begin(), find(eqtags_vec.begin(), eqtags_vec.end(), it)) + 1 << " ";
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

      ar_ec_output << endl
                   << "    % A0" << endl
                   << "    A0 = zeros(" << nontarget_lhs_vec.size() << ", " << nontarget_lhs_vec.size() << ");" << endl;
      for (const auto &[key, expr] : A0.at(name))
        {
          auto [eqn, colidx] = key;
          ar_ec_output << "    A0(" << eqn + 1 << ", " << colidx + 1 << ") = ";
          expr->writeOutput(ar_ec_output, ExprNodeOutputType::matlabDynamicModel);
          ar_ec_output << ";" << endl;
        }

      ar_ec_output << endl
                   << "    % A0star" << endl
                   << "    A0star = zeros(" << nontarget_lhs_vec.size() << ", " << target_lhs_vec.size() << ");" << endl;
      for (const auto &[key, expr] : A0star.at(name))
        {
          auto [eqn, colidx] = key;
          ar_ec_output << "    A0star(" << eqn + 1 << ", " << colidx + 1 << ") = ";
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
  for (bool printed_something{false};
       const auto &name : names)
    {
      if (exchange(printed_something, true))
        output << ", ";
      output << R"({"statementName": "trend_component_model",)"
             << R"("model_name": ")" << name << R"(",)"
             << R"("eqtags": [)";
      for (bool printed_something2{false};
           const auto &it : eqtags.at(name))
        {
          if (exchange(printed_something2, true))
            output << ", ";
          output << R"(")" << it << R"(")";
        }
      output << R"(], "target_eqtags": [)";
      for (bool printed_something2{false};
           const auto &it : target_eqtags.at(name))
        {
          if (exchange(printed_something2, true))
            output << ", ";
          output << R"(")" << it << R"(")";
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

  const filesystem::path filename {DataTree::packageDir(basename) / "varmatrices.m"};
  ofstream ar_output{filename, ios::out | ios::binary};
  if (!ar_output.is_open())
    {
      cerr << "ERROR: Can't open file " << filename.string() << " for writing" << endl;
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
             << "M_.var." << name << ".structural = " << boolalpha << structural.at(name) << ";" << endl
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
      for (bool it : diff.at(name))
        output << boolalpha << it << " ";
      output << "];" << endl
             << "M_.var." << name << ".nonstationary = M_.var." << name << ".diff;" << endl
             << "M_.var." << name << ".orig_diff_var = [";
      for (const auto &it : orig_diff_var.at(name))
        output << (it ? symbol_table.getTypeSpecificID(*it) + 1 : -1) << " ";
      output << "];" << endl;
      for (int i{1};
           const auto &it : rhs.at(name))
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
  for (bool printed_something{false};
       const auto &name : names)
    {
      if (exchange(printed_something, true))
        output << ", ";
      output << R"({"statementName": "var_model",)"
             << R"("model_name": ")" << name << R"(",)"
             << R"("eqtags": [)";
      for (bool printed_something2{false};
           const auto &it : eqtags.at(name))
        {
          if (exchange(printed_something2, true))
            output << ", ";
          output << R"(")" << it << R"(")";
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
  return eqtags.at(name_arg);
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
          if (symbol_table.isDiffAuxiliaryVariable(lhs_orig_symb_id))
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
  return diff.at(name_arg);
}

const vector<int> &
VarModelTable::getEqNums(const string &name_arg) const
{
  checkModelName(name_arg);
  return eqnums.at(name_arg);
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
VarModelTable::setOrigDiffVar(map<string, vector<optional<int>>> orig_diff_var_arg)
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
  return max_lags.at(name_arg);
}

int
VarModelTable::getMaxLag(const string &name_arg) const
{
  vector<int> maxlags { getMaxLags(name_arg) };
  return reduce(maxlags.begin(), maxlags.end(), 0, [](int a, int b) { return max(a, b); });
}

const vector<int> &
VarModelTable::getLhs(const string &name_arg) const
{
  checkModelName(name_arg);
  return lhs.at(name_arg);
}

const vector<int> &
VarModelTable::getLhsOrigIds(const string &name_arg) const
{
  checkModelName(name_arg);
  return lhs_orig_symb_ids.at(name_arg);
}

const vector<set<pair<int, int>>> &
VarModelTable::getRhs(const string &name_arg) const
{
  checkModelName(name_arg);
  return rhs.at(name_arg);
}

const vector<expr_t> &
VarModelTable::getLhsExprT(const string &name_arg) const
{
  checkModelName(name_arg);
  return lhs_expr_t.at(name_arg);
}


VarExpectationModelTable::VarExpectationModelTable(SymbolTable &symbol_table_arg) :
  symbol_table{symbol_table_arg}
{
}

void
VarExpectationModelTable::addVarExpectationModel(string name_arg, expr_t expression_arg, string aux_model_name_arg, string horizon_arg, expr_t discount_arg, int time_shift_arg)
{
  if (isExistingVarExpectationModelName(name_arg))
    {
      cerr << "Error: a var_expectation_model already exists with the name " << name_arg << endl;
      exit(EXIT_FAILURE);
    }

  expression[name_arg] = expression_arg;
  aux_model_name[name_arg] = move(aux_model_name_arg);
  horizon[name_arg] = move(horizon_arg);
  discount[name_arg] = discount_arg;
  time_shift[name_arg] = time_shift_arg;
  names.insert(move(name_arg));
}

bool
VarExpectationModelTable::isExistingVarExpectationModelName(const string &name_arg) const
{
  return names.contains(name_arg);
}

bool
VarExpectationModelTable::empty() const
{
  return names.empty();
}

void
VarExpectationModelTable::writeOutput(ostream &output) const
{
  for (const auto &name : names)
    {
      string mstruct = "M_.var_expectation." + name;
      output << mstruct << ".auxiliary_model_name = '" << aux_model_name.at(name) << "';" << endl
             << mstruct << ".horizon = " << horizon.at(name) << ';' << endl
             << mstruct << ".time_shift = " << time_shift.at(name) << ';' << endl;

      auto &vpc = vars_params_constants.at(name);
      if (!vpc.size())
        {
          cerr << "ERROR: VarExpectationModelStatement::writeOutput: matchExpression() has not been called" << endl;
          exit(EXIT_FAILURE);
        }

      ostringstream vars_list, params_list, constants_list;
      for (bool printed_something{false};
           const auto &[variable_id, param_id, constant] : vpc)
        {
          if (exchange(printed_something, true))
            {
              vars_list << ", ";
              params_list << ", ";
              constants_list << ", ";
            }
          vars_list << symbol_table.getTypeSpecificID(variable_id)+1;
          if (param_id)
            params_list << symbol_table.getTypeSpecificID(*param_id)+1;
          else
            params_list << "NaN";
          constants_list << constant;
        }
      output << mstruct << ".expr.vars = [ " << vars_list.str() << " ];" << endl
             << mstruct << ".expr.params = [ " << params_list.str() << " ];" << endl
             << mstruct << ".expr.constants = [ " << constants_list.str() << " ];" << endl;

      if (auto disc_var = dynamic_cast<const VariableNode *>(discount.at(name));
          disc_var)
        output << mstruct << ".discount_index = " << symbol_table.getTypeSpecificID(disc_var->symb_id) + 1 << ';' << endl;
      else
        {
          output << mstruct << ".discount_value = ";
          discount.at(name)->writeOutput(output);
          output << ';' << endl;
        }
      output << mstruct << ".param_indices = [ ";
      for (int param_id : aux_param_symb_ids.at(name))
        output << symbol_table.getTypeSpecificID(param_id)+1 << ' ';
      output << "];" << endl;
    }
}

void
VarExpectationModelTable::substituteUnaryOpsInExpression(const lag_equivalence_table_t &nodes, ExprNode::subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs)
{
  for (const auto &name : names)
    expression[name] = expression[name]->substituteUnaryOpNodes(nodes, subst_table, neweqs);
}

void
VarExpectationModelTable::substituteDiffNodesInExpression(const lag_equivalence_table_t &nodes, ExprNode::subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs)
{
  for (const auto &name : names)
    expression[name] = expression[name]->substituteDiff(nodes, subst_table, neweqs);
}

void
VarExpectationModelTable::transformPass(ExprNode::subst_table_t &diff_subst_table,
                                        DynamicModel &dynamic_model, const VarModelTable &var_model_table,
                                        const TrendComponentModelTable &trend_component_model_table)
{
  map<string, expr_t> var_expectation_subst_table;

  for (const auto &name : names)
    {
      // Collect information about the auxiliary model

      int max_lag;
      vector<int> lhs;
      if (var_model_table.isExistingVarModelName(aux_model_name[name]))
        {
          max_lag = var_model_table.getMaxLag(aux_model_name[name]);
          lhs = var_model_table.getLhs(aux_model_name[name]);
        }
      else if (trend_component_model_table.isExistingTrendComponentModelName(aux_model_name[name]))
        {
          max_lag = trend_component_model_table.getMaxLag(aux_model_name[name]) + 1;
          lhs = dynamic_model.getUndiffLHSForPac(aux_model_name[name], diff_subst_table);
        }
      else
        {
          cerr << "ERROR: var_expectation_model " << name
               << " refers to nonexistent auxiliary model " << aux_model_name[name] << endl;
          exit(EXIT_FAILURE);
        }

      // Match the linear combination in the expression option
      try
        {
          auto vpc = expression[name]->matchLinearCombinationOfVariables();
          for (const auto &[variable_id, lag, param_id, constant] : vpc)
            {
              if (lag != 0)
                throw ExprNode::MatchFailureException{"lead/lags are not allowed"};
              if (symbol_table.getType(variable_id) != SymbolType::endogenous)
                throw ExprNode::MatchFailureException{"Variable is not an endogenous"};
              vars_params_constants[name].emplace_back(variable_id, param_id, constant);
            }
        }
      catch (ExprNode::MatchFailureException &e)
        {
          cerr << "ERROR: expression in var_expectation_model " << name << " is not of the expected form: " << e.message << endl;
          exit(EXIT_FAILURE);
        }

      /* Create auxiliary parameters and the expression to be substituted into
         the var_expectations statement */
      expr_t subst_expr = dynamic_model.Zero;
      if (var_model_table.isExistingVarModelName(aux_model_name[name]))
        {
          /* If the auxiliary model is a VAR, add a parameter corresponding to
             the constant. */
          string constant_param_name = "var_expectation_model_" + name + "_constant";
          int constant_param_id = symbol_table.addSymbol(constant_param_name, SymbolType::parameter);
          aux_param_symb_ids[name].push_back(constant_param_id);
          subst_expr = dynamic_model.AddPlus(subst_expr, dynamic_model.AddVariable(constant_param_id));
        }
      for (int lag = 0; lag < max_lag; lag++)
        for (auto variable : lhs)
          {
            string param_name = "var_expectation_model_" + name + '_' + symbol_table.getName(variable) + '_' + to_string(lag);
            int new_param_id = symbol_table.addSymbol(param_name, SymbolType::parameter);
            aux_param_symb_ids[name].push_back(new_param_id);

            subst_expr = dynamic_model.AddPlus(subst_expr,
                                               dynamic_model.AddTimes(dynamic_model.AddVariable(new_param_id),
                                                                      dynamic_model.AddVariable(variable, -lag + time_shift[name])));
          }

      if (var_expectation_subst_table.contains(name))
        {
          cerr << "ERROR: model name '" << name << "' is used by several var_expectation_model statements" << endl;
          exit(EXIT_FAILURE);
        }
      var_expectation_subst_table[name] = subst_expr;
    }

  // Actually substitute var_expectation statements
  dynamic_model.substituteVarExpectation(var_expectation_subst_table);
  /* At this point, we know that all var_expectation operators have been
     substituted, because of the error check performed in
     VarExpectationNode::substituteVarExpectation(). */
}

void
VarExpectationModelTable::writeJsonOutput(ostream &output) const
{
  for (bool printed_something{false};
       const auto &name : names)
    {
      if (exchange(printed_something, true))
        output << ", ";
      output << R"({"statementName": "var_expectation_model",)"
             << R"("model_name": ")" << name << R"(", )"
             << R"("expression": ")";
      expression.at(name)->writeOutput(output);
      output << R"(", )"
             << R"("auxiliary_model_name": ")" << aux_model_name.at(name) << R"(", )"
             << R"("horizon": ")" << horizon.at(name) << R"(", )"
             << R"("discount": ")";
      discount.at(name)->writeOutput(output);
      output << R"(", )"
             << R"("time_shift": )" << time_shift.at(name)
             << R"(})";
    }
}


PacModelTable::PacModelTable(SymbolTable &symbol_table_arg) :
  symbol_table{symbol_table_arg}
{
}

void
PacModelTable::addPacModel(string name_arg, string aux_model_name_arg, string discount_arg, expr_t growth_arg, string auxname_arg, PacTargetKind kind_arg)
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
  auxname[name_arg] = move(auxname_arg);
  kind[name_arg] = kind_arg;
  names.insert(move(name_arg));
}

bool
PacModelTable::isExistingPacModelName(const string &name_arg) const
{
  return names.contains(name_arg);
}

bool
PacModelTable::empty() const
{
  return names.empty();
}

void
PacModelTable::checkPass(ModFileStructure &mod_file_struct)
{
  for (auto &[name, gv] : growth)
    if (gv)
      {
        if (target_info.contains(name))
          {
            cerr << "ERROR: for PAC model '" << name << "', it is not possible to declare a 'growth' option in the 'pac_model' command when there is also a 'pac_target_info' block" << endl;
            exit(EXIT_FAILURE);
          }
        gv->collectVariables(SymbolType::exogenous, mod_file_struct.pac_params);
      }

  for (auto &[name, auxn] : auxname)
    if (!auxn.empty() && target_info.contains(name))
      {
        cerr << "ERROR: for PAC model '" << name << "', it is not possible to declare an 'auxname' option in the 'pac_model' command when there is also a 'pac_target_info' block" << endl;
        exit(EXIT_FAILURE);
      }

  for (auto &[name, k] : kind)
    if (k != PacTargetKind::unspecified)
      {
        if (target_info.contains(name))
          {
            cerr << "ERROR: for PAC model '" << name << "', it is not possible to declare a 'kind' option in the 'pac_model' command when there is also a 'pac_target_info' block" << endl;
            exit(EXIT_FAILURE);
          }
        if (aux_model_name[name].empty())
          {
            cerr << "ERROR: for PAC model '" << name << "', it is not possible to declare a 'kind' option in the 'pac_model' command since this is a MCE model" << endl;
            exit(EXIT_FAILURE);
          }
      }

  for (const auto &[name, ti] : target_info)
    for (auto &[expr, gv, auxname, kind, coeff, growth_neutrality_param, h_indices, original_gv, gv_info] : get<2>(ti))
      if (gv)
        gv->collectVariables(SymbolType::exogenous, mod_file_struct.pac_params);

  for (const auto &[name, ti] : target_info)
    {
      auto &[target, auxname_target_nonstationary, components] = ti;
      if (!target)
        {
          cerr << "ERROR: the block 'pac_target_info(" << name << ")' is missing the 'target' statement" << endl;
          exit(EXIT_FAILURE);
        }
      if (auxname_target_nonstationary.empty())
        {
          cerr << "ERROR: the block 'pac_target_info(" << name << ")' is missing the 'auxname_target_nonstationary' statement" << endl;
          exit(EXIT_FAILURE);
        }
      int nonstationary_nb = 0;
      for (auto &[component, growth_component, auxname, kind, coeff, growth_neutrality_param, h_indices, original_growth_component, growth_component_info] : components)
        {
          if (auxname.empty())
            {
              cerr << "ERROR: the block 'pac_target_info(" << name << ")' is missing the 'auxname' statement in some 'component'" << endl;
              exit(EXIT_FAILURE);
            }
          if (kind == PacTargetKind::unspecified)
            {
              cerr << "ERROR: the block 'pac_target_info(" << name << ")' is missing the 'kind' statement in some 'component'" << endl;
              exit(EXIT_FAILURE);
            }
          if (kind == PacTargetKind::ll && growth_component)
            {
              cerr << "ERROR: in the block 'pac_target_info(" << name << ")', a component of 'kind ll' (i.e. stationary) has a 'growth' option. This is not permitted." << endl;
              exit(EXIT_FAILURE);
            }
          if (kind == PacTargetKind::dd || kind == PacTargetKind::dl)
            nonstationary_nb++;
        }
      if (!nonstationary_nb)
        {
          cerr << "ERROR: the block 'pac_target_info(" << name << ")' must contain at least one nonstationary component (i.e. of 'kind' equal to either 'dd' or 'dl')." << endl;
          exit(EXIT_FAILURE);
        }
    }
}

void
PacModelTable::substituteUnaryOpsInGrowth(const lag_equivalence_table_t &nodes, ExprNode::subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs)
{
  for (auto &[name, gv] : growth)
    if (gv)
      gv = gv->substituteUnaryOpNodes(nodes, subst_table, neweqs);

  for (auto &[name, ti] : target_info)
    for (auto &[expr, gv, auxname, kind, coeff, growth_neutrality_param, h_indices, original_gv, gv_info] : get<2>(ti))
      if (gv)
        gv = gv->substituteUnaryOpNodes(nodes, subst_table, neweqs);
}

void
PacModelTable::findDiffNodesInGrowth(lag_equivalence_table_t &diff_nodes) const
{
  for (auto &[name, gv] : growth)
    if (gv)
      gv->findDiffNodes(diff_nodes);

  for (const auto &[name, ti] : target_info)
    for (auto &[expr, gv, auxname, kind, coeff, growth_neutrality_param, h_indices, original_gv, gv_info] : get<2>(ti))
      if (gv)
        gv->findDiffNodes(diff_nodes);
}

void
PacModelTable::substituteDiffNodesInGrowth(const lag_equivalence_table_t &diff_nodes, ExprNode::subst_table_t &diff_subst_table, vector<BinaryOpNode *> &neweqs)
{
  for (auto &[name, gv] : growth)
    if (gv)
      gv = gv->substituteDiff(diff_nodes, diff_subst_table, neweqs);

  for (auto &[name, ti] : target_info)
    for (auto &[expr, gv, auxname, kind, coeff, growth_neutrality_param, h_indices, original_gv, gv_info] : get<2>(ti))
      if (gv)
        gv = gv->substituteDiff(diff_nodes, diff_subst_table, neweqs);
}

void
PacModelTable::transformPass(const lag_equivalence_table_t &unary_ops_nodes,
                             ExprNode::subst_table_t &unary_ops_subst_table,
                             const lag_equivalence_table_t &diff_nodes,
                             ExprNode::subst_table_t &diff_subst_table,
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
            growth_info[name] = growth[name]->matchLinearCombinationOfVariablesPlusConstant();
          }
        catch (ExprNode::MatchFailureException &e)
          {
            cerr << "ERROR: PAC growth must be a linear combination of variables" << endl;
            exit(EXIT_FAILURE);
          }

      // Perform transformations for the pac_target_info block (if any)
      if (target_info.contains(name))
        {
          // Substitute unary ops and diffs in the target…
          expr_t &target = get<0>(target_info[name]);
          vector<BinaryOpNode *> neweqs;
          target = target->substituteUnaryOpNodes(unary_ops_nodes, unary_ops_subst_table, neweqs);
          if (neweqs.size() > 0)
            {
              cerr << "ERROR: the 'target' expression of 'pac_target_info(" << name << ")' contains a variable with a unary operator that is not present in the model" << endl;
              exit(EXIT_FAILURE);
            }
          target = target->substituteDiff(diff_nodes, diff_subst_table, neweqs);
          if (neweqs.size() > 0)
            {
              cerr << "ERROR: the 'target' expression of 'pac_target_info(" << name << ")' contains a diff'd variable that is not present in the model" << endl;
              exit(EXIT_FAILURE);
            }

          // …and in component expressions
          auto &components = get<2>(target_info[name]);
          for (auto &[component, growth_component, auxname, kind, coeff, growth_neutrality_param, h_indices, original_growth_component, growth_component_info] : components)
            {
              component = component->substituteUnaryOpNodes(unary_ops_nodes, unary_ops_subst_table, neweqs);
              if (neweqs.size() > 0)
                {
                  cerr << "ERROR: a 'component' expression of 'pac_target_info(" << name << ")' contains a variable with a unary operator that is not present in the model" << endl;
                  exit(EXIT_FAILURE);
                }
              component = component->substituteDiff(diff_nodes, diff_subst_table, neweqs);
              if (neweqs.size() > 0)
                {
                  cerr << "ERROR: a 'component' expression of 'pac_target_info(" << name << ")' contains a diff'd variable that is not present in the model" << endl;
                  exit(EXIT_FAILURE);
                }
            }

          /* Fill the growth_info structure.
             Cannot be done in an earlier pass since growth terms can be
             transformed by DynamicModel::substituteDiff(). */
          for (auto &[component, growth_component, auxname, kind, coeff, growth_neutrality_param, h_indices, original_growth_component, growth_component_info] : components)
            {
              if (growth_component)
                try
                  {
                    growth_component_info = growth_component->matchLinearCombinationOfVariablesPlusConstant();
                  }
                catch (ExprNode::MatchFailureException &e)
                  {
                    cerr << "ERROR: PAC growth must be a linear combination of variables" << endl;
                    exit(EXIT_FAILURE);
                  }
            }

          // Identify the model equation defining the target
          expr_t target_expr;
          try
            {
              target_expr = dynamic_model.getRHSFromLHS(target);
            }
          catch (ExprNode::MatchFailureException)
            {
              cerr << "ERROR: there is no equation whose LHS is equal to the 'target' of 'pac_target_info(" << name << ")'" << endl;
              exit(EXIT_FAILURE);
            }

          // Substitute unary ops and diffs in that equation, before parsing (see dynare#1837)
          target_expr = target_expr->substituteUnaryOpNodes(unary_ops_nodes, unary_ops_subst_table, neweqs);
          if (neweqs.size() > 0)
            {
              cerr << "ERROR: the equation defining the target of 'pac_target_info(" << name << ")' contains a variable with a unary operator that is not present in the model" << endl;
              exit(EXIT_FAILURE);
            }
          target_expr = target_expr->substituteDiff(diff_nodes, diff_subst_table, neweqs);
          if (neweqs.size() > 0)
            {
              cerr << "ERROR: the equation defining the target of 'pac_target_info(" << name << ")' contains a diff'd variable that is not present in the model" << endl;
              exit(EXIT_FAILURE);
            }

          // Parse that model equation
          vector<pair<int, expr_t>> terms;
          expr_t constant;
          try
            {
              tie(terms, constant) = target_expr->matchLinearCombinationOfEndogenousWithConstant();
            }
          catch (ExprNode::MatchFailureException)
            {
              cerr << "ERROR: the model equation defining the 'target' of 'pac_target_info(" << name << ")' is not of the right form (should be a linear combination of endogenous variables)" << endl;
              exit(EXIT_FAILURE);
            }

          // Associate the coefficients of the linear combination with the right components
          for (auto [var, coeff] : terms)
            if (auto it = find_if(components.begin(), components.end(),
                                  [&, &var = var](const auto &v) { return get<0>(v) == dynamic_model.AddVariable(var); });
                it != components.end())
              get<4>(*it) = coeff;
            else
              {
                cerr << "ERROR: the model equation defining the 'target' of 'pac_target_info(" << name << ")' contains a variable (" << symbol_table.getName(var) << ") that is not declared as a 'component'" << endl;
                exit(EXIT_FAILURE);
              }

          // Verify that all declared components appear in that equation
          for (const auto &[component, growth_component, auxname, kind, coeff, growth_neutrality_param, h_indices, original_growth_component, growth_component_info] : components)
            if (!coeff)
              {
                cerr << "ERROR: a 'component' of 'pac_target_info(" << name << ")' does not appear in the model equation defining the 'target'" << endl;
                exit(EXIT_FAILURE);
              }

          /* Add the variable and equation defining the stationary part of the
             target. Note that it includes the constant. */
          expr_t yns = constant;
          for (const auto &[component, growth_component, auxname, kind, coeff, growth_neutrality_param, h_indices, original_growth_component, growth_component_info] : components)
            if (kind != PacTargetKind::ll)
              yns = dynamic_model.AddPlus(yns, dynamic_model.AddTimes(coeff, component));
          int target_nonstationary_id = symbol_table.addPacTargetNonstationaryAuxiliaryVar(get<1>(target_info[name]), yns);
          expr_t neweq = dynamic_model.AddEqual(dynamic_model.AddVariable(target_nonstationary_id), yns);
          dynamic_model.addEquation(neweq, nullopt);
          dynamic_model.addAuxEquation(neweq);

          /* Perform the substitution of the pac_target_nonstationary operator.
             This needs to be done here, otherwise
             DynamicModel::analyzePacEquationStructure() will not be able to
             identify the error-correction part */
          dynamic_model.substitutePacTargetNonstationary(name, dynamic_model.AddVariable(target_nonstationary_id, -1));
        }

      // Collect some information about PAC models
      int max_lag;
      if (trend_component_model_table.isExistingTrendComponentModelName(aux_model_name[name]))
        {
          aux_model_type[name] = "trend_component";
          max_lag = trend_component_model_table.getMaxLag(aux_model_name[name]) + 1;
          lhs[name] = dynamic_model.getUndiffLHSForPac(aux_model_name[name], diff_subst_table);
        }
      else if (var_model_table.isExistingVarModelName(aux_model_name[name]))
        {
          aux_model_type[name] = "var";
          max_lag = var_model_table.getMaxLag(aux_model_name[name]);
          lhs[name] = var_model_table.getLhs(aux_model_name[name]);
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

      // Compute the expressions that will be substituted for the pac_expectation operators
      expr_t growth_correction_term = dynamic_model.Zero;
      if (growth[name])
        growth_correction_term = dynamic_model.AddTimes(growth[name], dynamic_model.AddVariable(growth_neutrality_params[name]));
      if (aux_model_name[name].empty())
        {
          if (target_info.contains(name))
            {
              cerr << "ERROR: the block 'pac_target_info(" << name << ")' is not supported in the context of a PAC model with model-consistent expectations (MCE)." << endl;
              exit(EXIT_FAILURE);
            }
          else
            dynamic_model.computePacModelConsistentExpectationSubstitution(name,
                                                                           symbol_table.getID(discount[name]),
                                                                           pacEquationMaxLag(name),
                                                                           growth_correction_term,
                                                                           auxname[name],
                                                                           diff_subst_table,
                                                                           aux_var_symb_ids,
                                                                           aux_param_symb_ids,
                                                                           pac_expectation_substitution);
        }
      else
        {
          if (target_info.contains(name))
            {
              assert(growth_correction_term == dynamic_model.Zero);
              dynamic_model.computePacBackwardExpectationSubstitutionWithComponents(name, lhs[name],
                                                                                    max_lag,
                                                                                    aux_model_type[name],
                                                                                    get<2>(target_info[name]),
                                                                                    pac_expectation_substitution);
            }
          else
            dynamic_model.computePacBackwardExpectationSubstitution(name, lhs[name], max_lag,
                                                                    aux_model_type[name],
                                                                    growth_correction_term,
                                                                    auxname[name],
                                                                    aux_var_symb_ids,
                                                                    aux_param_symb_ids,
                                                                    pac_expectation_substitution);
        }
    }

  // Actually perform the substitution of pac_expectation
  dynamic_model.substitutePacExpectation(pac_expectation_substitution, eq_name);
  dynamic_model.checkNoRemainingPacExpectation();

  // Check that there is no remaining pac_target_nonstationary operator
  dynamic_model.checkNoRemainingPacTargetNonstationary();
}

void
PacModelTable::writeOutput(ostream &output) const
{
  // Helper to print the “growth_info” structure (linear decomposition of growth)
  auto growth_info_helper = [&](const string &fieldname, const growth_info_t &gi)
  {
    for (int i{1};
         auto [growth_symb_id, growth_lag, param_id, constant] : gi)
      {
        string structname = fieldname + "(" + to_string(i++) + ").";
        if (growth_symb_id)
          {
            string var_field = "endo_id";
            if (symbol_table.getType(*growth_symb_id) == SymbolType::exogenous)
              {
                var_field = "exo_id";
                output << structname << "endo_id = 0;" << endl;
              }
            else
              output << structname << "exo_id = 0;" << endl;
            try
              {
                // case when this is not the highest lag of the growth variable
                int aux_symb_id = symbol_table.searchAuxiliaryVars(*growth_symb_id, growth_lag);
                output << structname << var_field << " = " << symbol_table.getTypeSpecificID(aux_symb_id) + 1 << ";" << endl
                       << structname << "lag = 0;" << endl;
              }
            catch (...)
              {
                try
                  {
                    // case when this is the highest lag of the growth variable
                    int tmp_growth_lag = growth_lag + 1;
                    int aux_symb_id = symbol_table.searchAuxiliaryVars(*growth_symb_id, tmp_growth_lag);
                    output << structname << var_field << " = " << symbol_table.getTypeSpecificID(aux_symb_id) + 1 << ";" << endl
                           << structname << "lag = -1;" << endl;
                  }
                catch (...)
                  {
                    // case when there is no aux var for the variable
                    output << structname << var_field << " = "<< symbol_table.getTypeSpecificID(*growth_symb_id) + 1 << ";" << endl
                           << structname << "lag = " << growth_lag << ";" << endl;
                  }
              }
          }
        else
          output << structname << "endo_id = 0;" << endl
                 << structname << "exo_id = 0;" << endl
                 << structname << "lag = 0;" << endl;
        output << structname << "param_id = "
               << (param_id ? symbol_table.getTypeSpecificID(*param_id) + 1 : 0) << ";" << endl
               << structname << "constant = " << constant << ";" << endl;
      }
  };

  for (const auto &name : names)
    {
      output << "M_.pac." << name << ".auxiliary_model_name = '" << aux_model_name.at(name) << "';" << endl
             << "M_.pac." << name << ".discount_index = " << symbol_table.getTypeSpecificID(discount.at(name)) + 1 << ";" << endl;

      if (growth.at(name))
        {
          output << "M_.pac." << name << ".growth_str = '";
          original_growth.at(name)->writeJsonOutput(output, {}, {}, true);
          output << "';" << endl;
          growth_info_helper("M_.pac." + name + ".growth_linear_comb", growth_info.at(name));
        }
    }

  // Write the auxiliary parameter IDs created for the pac_expectation operator
  for (auto &[name, ids] : aux_param_symb_ids)
    {
      output << "M_.pac." << name << "." << (aux_model_name.at(name).empty() ? "mce.alpha" : "h_param_indices") << " = [";
      for (auto id : ids)
        output << symbol_table.getTypeSpecificID(id) + 1 << " ";
      output << "];" << endl;
    }

  // Write the auxiliary variable IDs created for the pac_expectation operator
  for (auto &[name, id] : aux_var_symb_ids)
    output << "M_.pac." << name << "." << (aux_model_name.at(name).empty() ? "mce.z1" : "aux_id")
           << " = " << symbol_table.getTypeSpecificID(id) + 1 << ";" << endl;

  // Write PAC equation name info
  for (auto &[name, eq] : eq_name)
    output << "M_.pac." << name << ".eq_name = '" << eq << "';" << endl;

  for (auto &[model, growth_neutrality_param_index] : growth_neutrality_params)
    output << "M_.pac." << model << ".growth_neutrality_param_index = "
           << symbol_table.getTypeSpecificID(growth_neutrality_param_index) + 1 << ";" << endl;

  for (auto &[model, type] : aux_model_type)
      output << "M_.pac." << model << ".auxiliary_model_type = '" << type << "';" << endl;

  for (auto &[name, k] : kind)
    if (!aux_model_name.empty())
      output << "M_.pac." << name << ".kind = '"
             << (k == PacTargetKind::unspecified ? "" : kindToString(k)) << "';" << endl;

  for (auto &[name, val] : equation_info)
    {
      auto &[lhs_pac_var, optim_share_index, ar_params_and_vars, ec_params_and_vars, non_optim_vars_params_and_constants, additive_vars_params_and_constants, optim_additive_vars_params_and_constants] = val;
      output << "M_.pac." << name << ".lhs_var = "
             << symbol_table.getTypeSpecificID(lhs_pac_var.first) + 1 << ";" << endl;

      if (optim_share_index)
        output << "M_.pac." << name << ".share_of_optimizing_agents_index = "
               << symbol_table.getTypeSpecificID(*optim_share_index) + 1 << ";" << endl;

      output << "M_.pac." << name << ".ec.params = "
             << symbol_table.getTypeSpecificID(ec_params_and_vars.first) + 1 << ";" << endl
             << "M_.pac." << name << ".ec.vars = [";
      for (auto &it : ec_params_and_vars.second)
        output << symbol_table.getTypeSpecificID(get<0>(it)) + 1 << " ";
      output << "];" << endl
             << "M_.pac." << name << ".ec.istarget = [";
      for (auto &it : ec_params_and_vars.second)
        output << boolalpha << get<1>(it) << " ";
      output << "];" << endl
             << "M_.pac." << name << ".ec.scale = [";
      for (auto &it : ec_params_and_vars.second)
        output << get<2>(it) << " ";
      output << "];" << endl
             << "M_.pac." << name << ".ec.isendo = [";
      for (auto &it : ec_params_and_vars.second)
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
        output << (pid ? symbol_table.getTypeSpecificID(*pid) + 1 : -1) << " ";
      output << "];" << endl
             << "M_.pac." << name << ".ar.vars = [";
      for (auto &[pid, vid, vlag] : ar_params_and_vars)
        output << (vid ? symbol_table.getTypeSpecificID(*vid) + 1 : -1) << " ";
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
            if (get<2>(it))
              output << symbol_table.getTypeSpecificID(*get<2>(it)) + 1 << " ";
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
            if (get<2>(it))
              output << symbol_table.getTypeSpecificID(*get<2>(it)) + 1 << " ";
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
            if (get<2>(it))
              output << symbol_table.getTypeSpecificID(*get<2>(it)) + 1 << " ";
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
    }

  for (auto &[name, val] : target_info)
    for (int component_idx{1};
         auto &[component, growth_component, auxname, kind, coeff, growth_neutrality_param, h_indices, original_growth_component, growth_component_info] : get<2>(val))
      {
        string fieldname = "M_.pac." + name + ".components(" + to_string(component_idx) + ")";
        output << fieldname << ".aux_id = " << symbol_table.getTypeSpecificID(auxname) + 1 << ";" << endl
               << fieldname << ".endo_var = " << symbol_table.getTypeSpecificID(dynamic_cast<VariableNode *>(component)->symb_id) + 1 << ";" << endl
               << fieldname << ".kind = '" << kindToString(kind) << "';" << endl
               << fieldname << ".h_param_indices = [";
        for (int id : h_indices)
          output << symbol_table.getTypeSpecificID(id) + 1 << " ";
        output << "];" << endl
               << fieldname << ".coeff_str = '";
        coeff->writeJsonOutput(output, {}, {}, true);
        output << "';" << endl;
        if (growth_component)
          {
            output << fieldname << ".growth_neutrality_param_index = " << symbol_table.getTypeSpecificID(growth_neutrality_param) + 1 << ";" << endl
                   << fieldname << ".growth_str = '";
            original_growth_component->writeJsonOutput(output, {}, {}, true);
            output << "';" << endl;
            growth_info_helper(fieldname + ".growth_linear_comb", growth_component_info);
          }
        component_idx++;
      }
}

void
PacModelTable::writeJsonOutput(ostream &output) const
{
  for (bool printed_something{false};
       const auto &name : names)
    {
      /* The calling method has already added a comma, so don’t output one for
         the first statement */
      if (exchange(printed_something, true))
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

  for (auto &[name, val] : target_info)
    {
      output << R"(, {"statementName": "pac_target_info", "model_name": ")" << name
             << R"(", "target": ")";
      get<0>(val)->writeJsonOutput(output, {}, {}, true);
      output << R"(", "auxname_target_nonstationary": ")" << get<1>(val)
             << R"(", "components": [)";
      for (auto &[component, growth_component, auxname, kind, coeff, growth_neutrality_param, h_indices, original_growth_component, growth_component_info] : get<2>(val))
        {
          if (component != get<0>(get<2>(val).front()))
            output << ", ";
          output << R"({"component": ")";
          component->writeJsonOutput(output, {}, {}, true);
          output << R"(", "auxname": ")" << auxname
                 << R"(", "kind": ")" << kindToString(kind);
          if (growth_component)
            {
              output << R"(", "growth_str": ")";
              original_growth_component->writeJsonOutput(output, {}, {}, true);
            }
          output << R"("})";
        }
      output << "]}" << endl;
    }
}

int
PacModelTable::pacEquationMaxLag(const string &name_arg) const
{
  return get<2>(equation_info.at(name_arg)).size();
}

string
PacModelTable::kindToString(PacTargetKind kind)
{
  switch (kind)
    {
    case PacTargetKind::unspecified:
      cerr << "Internal error: kind should not be unspecified" << endl;
      exit(EXIT_FAILURE);
    case PacTargetKind::ll:
      return "ll";
    case PacTargetKind::dl:
      return "dl";
    case PacTargetKind::dd:
      return "dd";
    }
  // Silent GCC warning
  assert(false);
}

void
PacModelTable::setTargetExpr(const string &name_arg, expr_t target)
{
  get<0>(target_info[name_arg]) = target;
}

void
PacModelTable::setTargetAuxnameNonstationary(const string &name_arg, string auxname)
{
  get<1>(target_info[name_arg]) = move(auxname);
}

void
PacModelTable::addTargetComponent(const string &name_arg, target_component_t component)
{
  get<7>(component) = get<1>(component); // original_growth = growth
  get<2>(target_info[name_arg]).emplace_back(move(component));
}

void
PacModelTable::writeTargetCoefficientsFile(const string &basename) const
{
  if (target_info.empty())
    return;

  filesystem::path filename {DataTree::packageDir(basename) / "pac_target_coefficients.m"};
  ofstream output{filename, ios::out | ios::binary};
  if (!output.is_open())
    {
      cerr << "ERROR: Can't open file " << filename.string() << " for writing" << endl;
      exit(EXIT_FAILURE);
    }
  output << "function coeffs = pac_target_coefficients(model_name, params)" << endl;
  for (auto &[model_name, val] : target_info)
    {
      output << "  if strcmp(model_name, '" << model_name << "')" << endl
             << "    coeffs = NaN(" << get<2>(val).size() << ",1);" << endl;
      for (int i{1};
           auto &[component, growth_component, auxname, kind, coeff, growth_neutrality_param, h_indices, original_growth_component, growth_component_info] : get<2>(val))
        {
          output << "    coeffs(" << i++ << ") = ";
          coeff->writeOutput(output, ExprNodeOutputType::matlabDynamicModel);
          output << ";" << endl;
        }
      output << "    return" << endl
             << "  end" << endl;
    }
  output << "  error([ 'Unknown PAC model: ' model_name ])" << endl
         << "end" << endl;
  output.close();
}
