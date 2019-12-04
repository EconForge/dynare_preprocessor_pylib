/*
 * Copyright © 2010-2019 Dynare Team
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
#include <algorithm>

#include "ModelEquationBlock.hh"

SteadyStateModel::SteadyStateModel(SymbolTable &symbol_table_arg,
                                   NumericalConstants &num_constants_arg,
                                   ExternalFunctionsTable &external_functions_table_arg,
                                   const StaticModel &static_model_arg) :
  DataTree{symbol_table_arg, num_constants_arg, external_functions_table_arg},
  static_model{static_model_arg}
{
}

SteadyStateModel::SteadyStateModel(const SteadyStateModel &m) :
  DataTree {m},
  static_model {m.static_model}
{
  for (const auto &it : m.def_table)
    def_table.emplace_back(it.first, it.second->clone(*this));
}

SteadyStateModel &
SteadyStateModel::operator=(const SteadyStateModel &m)
{
  DataTree::operator=(m);

  assert(&static_model == &m.static_model);

  def_table.clear();
  for (const auto &it : m.def_table)
    def_table.emplace_back(it.first, it.second->clone(*this));

  return *this;
}

void
SteadyStateModel::addDefinition(int symb_id, expr_t expr)
{
  AddVariable(symb_id); // Create the variable node to be used in write method

  assert(symbol_table.getType(symb_id) == SymbolType::endogenous
         || symbol_table.getType(symb_id) == SymbolType::modFileLocalVariable
         || symbol_table.getType(symb_id) == SymbolType::parameter);

  // Add the variable
  vector<int> v;
  v.push_back(symb_id);
  def_table.emplace_back(v, expr);
}

void
SteadyStateModel::addMultipleDefinitions(const vector<int> &symb_ids, expr_t expr)
{
  for (int symb_id : symb_ids)
    {
      AddVariable(symb_id); // Create the variable nodes to be used in write method
      assert(symbol_table.getType(symb_id) == SymbolType::endogenous
             || symbol_table.getType(symb_id) == SymbolType::modFileLocalVariable
             || symbol_table.getType(symb_id) == SymbolType::parameter);
    }
  def_table.emplace_back(symb_ids, expr);
}

void
SteadyStateModel::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) const
{
  if (def_table.size() == 0)
    return;

  mod_file_struct.steady_state_model_present = true;
  vector<int> so_far_defined;

  for (const auto & i : def_table)
    {
      const vector<int> &symb_ids = i.first;

      // Check that symbols are not already defined
      for (int symb_id : symb_ids)
        if (find(so_far_defined.begin(), so_far_defined.end(), symb_id)
            != so_far_defined.end())
          warnings << "WARNING: in the 'steady_state_model' block, variable '" << symbol_table.getName(symb_id) << "' is declared twice" << endl;

      // Check that expression has no undefined symbol
      if (!mod_file_struct.ramsey_model_present)
        {
          set<int> used_symbols;
          const expr_t &expr = i.second;
          expr->collectVariables(SymbolType::endogenous, used_symbols);
          expr->collectVariables(SymbolType::modFileLocalVariable, used_symbols);
          for (int used_symbol : used_symbols)
            if (find(so_far_defined.begin(), so_far_defined.end(), used_symbol)
                == so_far_defined.end())
              {
                cerr << "ERROR: in the 'steady_state_model' block, variable '" << symbol_table.getName(used_symbol)
                     << "' is undefined in the declaration of variable '" << symbol_table.getName(symb_ids[0]) << "'" << endl;
                exit(EXIT_FAILURE);
              }
        }

      copy(symb_ids.begin(), symb_ids.end(), back_inserter(so_far_defined));
    }

  set<int> orig_endogs = symbol_table.getOrigEndogenous();
  for (int orig_endog : orig_endogs)
    {
      if (find(so_far_defined.begin(), so_far_defined.end(), orig_endog)
          == so_far_defined.end())
        warnings << "WARNING: in the 'steady_state_model' block, variable '" << symbol_table.getName(orig_endog) << "' is not assigned a value" << endl;
    }
}

void
SteadyStateModel::writeLatexSteadyStateFile(const string &basename) const
{
  filesystem::create_directories(basename + "/latex");

  ofstream output, content_output;
  string filename = basename + "/latex/steady_state.tex";
  string content_filename = basename + "/latex/steady_state_content.tex";

  output.open(filename, ios::out | ios::binary);
  if (!output.is_open())
    {
      cerr << "ERROR: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }

  content_output.open(content_filename, ios::out | ios::binary);
  if (!content_output.is_open())
    {
      cerr << "ERROR: Can't open file " << content_filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }

  output << "\\documentclass[10pt,a4paper]{article}" << endl
         << "\\usepackage[landscape]{geometry}" << endl
         << "\\usepackage{fullpage}" << endl
         << "\\usepackage{amsfonts}" << endl
         << "\\usepackage{breqn}" << endl
         << "\\begin{document}" << endl
         << "\\footnotesize" << endl;

  for (const auto & it : def_table)
    for (auto it1 = it.first.begin(); it1 != it.first.end(); it1++)
      {
        int id = *it1;
        expr_t value = it.second;
        content_output << "\\begin{dmath}" << endl
                       << symbol_table.getTeXName(id) << " = ";
        value->writeOutput(content_output, ExprNodeOutputType::latexStaticModel);
        content_output << endl << "\\end{dmath}" << endl;
      }

  static_model.writeLatexAuxVarRecursiveDefinitions(content_output);

  output << "\\include{steady_state_content.tex}" << endl
         << "\\end{document}" << endl;

  output.close();
  content_output.close();
}

void
SteadyStateModel::writeSteadyStateFile(const string &basename, bool ramsey_model, bool julia) const
{
  if (def_table.size() == 0)
    return;

  string filename = julia ? basename + "SteadyState2.jl" : packageDir(basename) + "/steadystate.m";
  ofstream output;
  output.open(filename, ios::out | ios::binary);
  if (!output.is_open())
    {
      cerr << "ERROR: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }

  ExprNodeOutputType output_type = (julia ? ExprNodeOutputType::juliaSteadyStateFile : ExprNodeOutputType::steadyStateFile);

  if (!julia)
    output << "function [ys_, params, info] = steadystate("
           << "ys_, exo_, params)" << endl
           << "% Steady state generated by Dynare preprocessor" << endl
           << "    info = 0;" << endl;
  else
    output << "module " << basename << "SteadyState2" << endl
           << "#" << endl
           << "# NB: this file was automatically generated by Dynare" << endl
           << "#     from " << basename << ".mod" << endl
           << "#" << endl
           << "export steady_state!" << endl << endl
           << "function steady_state!(ys_::Vector{Float64}, exo_::Vector{Float64}, "
           << "params::Vector{Float64})" << endl;

  for (const auto & i : def_table)
    {
      const vector<int> &symb_ids = i.first;
      output << "    ";
      if (symb_ids.size() > 1)
        output << "[";
      for (size_t j = 0; j < symb_ids.size(); j++)
        {
          getVariable(symb_ids[j])->ExprNode::writeOutput(output, output_type);
          if (j < symb_ids.size()-1)
            output << ",";
        }
      if (symb_ids.size() > 1)
        output << "]";

      output << "=";
      i.second->writeOutput(output, output_type);
      output << ";" << endl;
    }
  if (!julia)
    output << "    % Auxiliary equations" << endl;
  else
    output << "    # Auxiliary equations" << endl;
  static_model.writeAuxVarRecursiveDefinitions(output, output_type);

  output << "end" << endl;
  if (julia)
    output << "end" << endl;
  output.close();
}

void
SteadyStateModel::writeJsonSteadyStateFile(ostream &output, bool transformComputingPass) const
{
  if (def_table.size() == 0)
    return;

  vector<pair<string, string>> eqtags;

  output << "{\"steady_state_model\": [";

  for (size_t i = 0; i < def_table.size(); i++)
    {
      const vector<int> &symb_ids = def_table[i].first;
      if (i != 0)
        output << ",";
      output << "{\"lhs\": ";
      if (symb_ids.size() > 1)
        output << "[";
      for (size_t j = 0; j < symb_ids.size(); j++)
        {
          if (j != 0)
            output << ",";
          output << "\"";
          getVariable(symb_ids[j])->writeJsonOutput(output, {}, {}, false);
          output << "\"";
        }
      if (symb_ids.size() > 1)
        output << "]";
      output << ", \"rhs\":\"";
      def_table[i].second->writeJsonOutput(output, {}, {}, false);
      output << "\"}" << endl;
    }

  if (transformComputingPass)
    static_model.writeJsonAuxVarRecursiveDefinitions(output);

  output << "]}";
}

Epilogue::Epilogue(SymbolTable &symbol_table_arg,
                   NumericalConstants &num_constants_arg,
                   ExternalFunctionsTable &external_functions_table_arg,
                   TrendComponentModelTable &trend_component_model_table_arg,
                   VarModelTable &var_model_table_arg) :
  DynamicModel{symbol_table_arg, num_constants_arg, external_functions_table_arg,
               trend_component_model_table_arg, var_model_table_arg}
{
}

Epilogue::Epilogue(const Epilogue &m) :
  DynamicModel {m}
{
  for (const auto &it : m.def_table)
    def_table.emplace_back(it.first, it.second->clone(*this));
}

Epilogue &
Epilogue::operator=(const Epilogue &m)
{
  DynamicModel::operator=(m);

  def_table.clear();
  for (const auto &it : m.def_table)
    def_table.emplace_back(it.first, it.second->clone(*this));

  return *this;
}

void
Epilogue::addDefinition(int symb_id, expr_t expr)
{
  def_table.emplace_back(symb_id, expr);
}

void
Epilogue::checkPass(WarningConsolidation &warnings) const
{
  if (def_table.size() == 0)
    return;

  vector<int> so_far_defined;
  for (const auto & it : def_table)
    if (find(so_far_defined.begin(), so_far_defined.end(), it.first) != so_far_defined.end())
      {
        cerr << "WARNING: in the 'epilogue' block, variable '" << it.first
             << "' is declared twice" << endl;
        exit(EXIT_FAILURE);
      }
    else
      so_far_defined.push_back(it.first);
}

void
Epilogue::detrend(const map<int, expr_t> & trend_symbols_map,
                  const nonstationary_symbols_map_t & nonstationary_symbols_map)
{
  for (auto it = nonstationary_symbols_map.crbegin();
       it != nonstationary_symbols_map.crend(); it++)
    for (auto & [symb_id, expr] : def_table)
      {
        expr = expr->detrend(it->first, it->second.first, it->second.second);
        assert(expr != nullptr);
      }

  for (auto & [symb_id, expr] : def_table)
    {
      expr = expr->removeTrendLeadLag(trend_symbols_map);
      assert(expr != nullptr);
    }

  for (auto & [symb_id, expr] : def_table)
    {
      expr = expr->replaceTrendVar();
      assert(expr != nullptr);
    }
}

void
Epilogue::writeEpilogueFile(const string &basename) const
{
  if (def_table.size() == 0)
    return;

  string filename = packageDir(basename) + "/epilogue.m";
  ofstream output;
  output.open(filename, ios::out | ios::binary);
  if (!output.is_open())
    {
      cerr << "ERROR: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }

  ExprNodeOutputType output_type = ExprNodeOutputType::epilogueFile;
  output << "function dseries__ = epilogue(params, dseries__)" << endl
         << "% function dseries__ = epilogue(params, dseries__)" << endl
         << "% Epilogue file generated by Dynare preprocessor" << endl << endl
         << "simul_end_date = lastdate(dseries__);" << endl;

  deriv_node_temp_terms_t tef_terms;
  temporary_terms_t temporary_terms;
  temporary_terms_idxs_t temporary_terms_idxs;
  for (const auto & it : def_table)
    if (it.second->containsExternalFunction())
      it.second->writeExternalFunctionOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
  output << endl;
  for (const auto & it : def_table)
    {
      int max_lag = it.second->maxLagWithDiffsExpanded();
      set<int> used_symbols;
      it.second->collectVariables(SymbolType::endogenous, used_symbols);
      it.second->collectVariables(SymbolType::exogenous, used_symbols);
      it.second->collectVariables(SymbolType::epilogue, used_symbols);

      output << "simul_begin_date = firstobservedperiod(dseries__{";
      for (auto it1 = used_symbols.begin(); it1 != used_symbols.end(); it1++)
        {
          if (it1 != used_symbols.begin())
            output << ", ";
          output << "'" << symbol_table.getName(*it1) << "'";
        }
      output << "}) + " << max_lag << " + 1;" << endl
             << "if ~dseries__.exist('" << symbol_table.getName(it.first) << "')" << endl
             << "    dseries__ = [dseries__ dseries(NaN(dseries__.nobs,1), dseries__.firstdate, '" << symbol_table.getName(it.first)<< "')];" << endl
             << "end" << endl
             << "from simul_begin_date to simul_end_date do "
             << "dseries__." << symbol_table.getName(it.first) << "(t) = ";
      it.second->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
      output << ";" << endl << endl;
    }
  output << "end" << endl;
  output.close();
}
