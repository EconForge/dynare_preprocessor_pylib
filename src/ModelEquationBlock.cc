/*
 * Copyright © 2010-2021 Dynare Team
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

#include <cassert>
#include <algorithm>
#include <sstream>

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
  DataTree{m},
  static_model{m.static_model}
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
  vector v{symb_id};
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
  set<int> so_far_defined;

  for (const auto &[symb_ids, expr] : def_table)
    {
      // Check that symbols are not already defined
      for (int symb_id : symb_ids)
        if (so_far_defined.contains(symb_id))
          warnings << "WARNING: in the 'steady_state_model' block, variable '" << symbol_table.getName(symb_id) << "' is declared twice" << endl;

      // Check that expression has no undefined symbol
      if (!mod_file_struct.ramsey_model_present)
        {
          set<int> used_symbols;
          expr->collectVariables(SymbolType::endogenous, used_symbols);
          expr->collectVariables(SymbolType::modFileLocalVariable, used_symbols);
          for (int used_symbol : used_symbols)
            if (!so_far_defined.contains(used_symbol))
              {
                cerr << "ERROR: in the 'steady_state_model' block, variable '" << symbol_table.getName(used_symbol)
                     << "' is undefined in the declaration of variable '" << symbol_table.getName(symb_ids[0]) << "'" << endl;
                exit(EXIT_FAILURE);
              }
        }

      so_far_defined.insert(symb_ids.begin(), symb_ids.end());
    }

  /* Check that all original endogous are defined (except the instruments of a
     Ramsey model, since the steady_state_block should give the steady state
     *conditional* to those instruments) */
  set<int> should_be_defined = symbol_table.getOrigEndogenous();
  if (mod_file_struct.ramsey_policy_present || mod_file_struct.ramsey_model_present)
    for (const auto &s : mod_file_struct.instruments.getSymbols())
      should_be_defined.erase(symbol_table.getID(s));
  for (int v : should_be_defined)
    if (!so_far_defined.contains(v))
      warnings << "WARNING: in the 'steady_state_model' block, variable '" << symbol_table.getName(v) << "' is not assigned a value" << endl;
}

void
SteadyStateModel::writeLatexSteadyStateFile(const string &basename) const
{
  filesystem::create_directories(basename + "/latex");

  string filename = basename + "/latex/steady_state.tex";
  string content_filename = basename + "/latex/steady_state_content.tex";

  ofstream output{filename, ios::out | ios::binary};
  if (!output.is_open())
    {
      cerr << "ERROR: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }

  ofstream content_output{content_filename, ios::out | ios::binary};
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

  for (const auto & [ids, value] : def_table)
    for (int id : ids)
      {
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
SteadyStateModel::writeSteadyStateFile(const string &basename, bool julia) const
{
  if (def_table.size() == 0)
    return;

  ExprNodeOutputType output_type = (julia ? ExprNodeOutputType::juliaSteadyStateFile : ExprNodeOutputType::steadyStateFile);

  stringstream output;
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
           << "function steady_state!(ys_::Vector{<: Real}, exo_::Vector{<: Real}, "
           << "params::Vector{<: Real})" << endl;

  for (const auto & [symb_ids, value] : def_table)
    {
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
      value->writeOutput(output, output_type);
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

  writeToFileIfModified(output, julia ? basename + "SteadyState2.jl" : packageDir(basename) + "/steadystate.m");
}

void
SteadyStateModel::writeJsonSteadyStateFile(ostream &output, bool transformComputingPass) const
{
  if (def_table.size() == 0)
    return;

  vector<pair<string, string>> eqtags;

  output << "{\"steady_state_model\": [";

  for (bool printed_something{false};
       const auto &[symb_ids, value] : def_table)
    {
      if (exchange(printed_something, true))
        output << ",";
      output << "{\"lhs\": ";
      if (symb_ids.size() > 1)
        output << "[";
      for (bool printed_something2{false};
           int symb_id : symb_ids)
        {
          if (exchange(printed_something2, true))
            output << ",";
          output << "\"";
          getVariable(symb_id)->writeJsonOutput(output, {}, {}, false);
          output << "\"";
        }
      if (symb_ids.size() > 1)
        output << "]";
      output << R"(, "rhs":")";
      value->writeJsonOutput(output, {}, {}, false);
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
  DynamicModel{m}
{
  for (const auto &it : m.dynamic_def_table)
    dynamic_def_table.emplace_back(it.first, it.second->clone(*this));
}

Epilogue &
Epilogue::operator=(const Epilogue &m)
{
  DynamicModel::operator=(m);

  dynamic_def_table.clear();
  for (const auto &it : m.dynamic_def_table)
    dynamic_def_table.emplace_back(it.first, it.second->clone(*this));

  return *this;
}

void
Epilogue::addDefinition(int symb_id, expr_t expr)
{
  dynamic_def_table.emplace_back(symb_id, expr);
}

void
Epilogue::checkPass(ModFileStructure &mod_file_struct) const
{
  if (dynamic_def_table.size() == 0)
    {
      if (mod_file_struct.with_epilogue_option)
        {
          cerr << "ERROR: the 'with_epilogue' option cannot be specified when there is no 'epilogue' block" << endl;
          exit(EXIT_FAILURE);
        }
      return;
    }

  set<int> so_far_defined;
  for (const auto &[symb_id, expr] : dynamic_def_table)
    if (so_far_defined.contains(symb_id))
      {
        cerr << "WARNING: in the 'epilogue' block, variable '" << symbol_table.getName(symb_id)
             << "' is declared twice" << endl;
        exit(EXIT_FAILURE);
      }
    else
      so_far_defined.insert(symb_id);
}

void
Epilogue::toStatic()
{
  for (const auto & [symb_id, expr] : dynamic_def_table)
    static_def_table.emplace_back(symb_id, expr->toStatic(*this));
}

void
Epilogue::detrend(const map<int, expr_t> &trend_symbols_map,
                  const nonstationary_symbols_map_t &nonstationary_symbols_map)
{
  for (auto it = nonstationary_symbols_map.crbegin();
       it != nonstationary_symbols_map.crend(); ++it)
    for (auto & [symb_id, expr] : dynamic_def_table)
      {
        expr = expr->detrend(it->first, it->second.first, it->second.second);
        assert(expr);
      }

  for (auto & [symb_id, expr] : dynamic_def_table)
    {
      expr = expr->removeTrendLeadLag(trend_symbols_map);
      assert(expr);
    }

  for (auto & [symb_id, expr] : dynamic_def_table)
    {
      expr = expr->replaceTrendVar();
      assert(expr);
    }
}

void
Epilogue::writeEpilogueFile(const string &basename) const
{
  if (dynamic_def_table.empty())
    return;

  writeDynamicEpilogueFile(basename);
  writeStaticEpilogueFile(basename);
}

void
Epilogue::writeStaticEpilogueFile(const string &basename) const
{
  string filename = packageDir(basename) + "/epilogue_static.m";
  ofstream output{filename, ios::out | ios::binary};
  if (!output.is_open())
    {
      cerr << "ERROR: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }

  output << "function ds = epilogue_static(params, ds)" << endl
         << "% function ds = epilogue_static(params, ds)" << endl
         << "% Epilogue file generated by Dynare preprocessor" << endl;

  for (const auto & [symb_id, expr] : static_def_table)
    {
      // Rewrite external function TEF term for every equation as argument values could have been changed
      // in between two calls to the same function;
      deriv_node_temp_terms_t tef_terms;
      temporary_terms_t temporary_terms;
      temporary_terms_idxs_t temporary_terms_idxs;
      output << endl;
      if (expr->containsExternalFunction())
        expr->writeExternalFunctionOutput(output, ExprNodeOutputType::matlabDseries, temporary_terms, temporary_terms_idxs, tef_terms);
      output << "epilogue_static_tmp_term = ";
      expr->writeOutput(output, ExprNodeOutputType::matlabDseries, temporary_terms, temporary_terms_idxs, tef_terms);
      output << ";" << endl
             << "if isdseries(epilogue_static_tmp_term)" << endl
             << "    ds." << symbol_table.getName(symb_id) << " = epilogue_static_tmp_term;" << endl
             << "else" << endl
             << "    ds." << symbol_table.getName(symb_id) << " = dseries(ones(ds.nobs,1)*epilogue_static_tmp_term, ds.firstdate, '" << symbol_table.getName(symb_id) << "');" << endl
             << "end" << endl;
    }
  output << "end" << endl;
  output.close();
}

void
Epilogue::writeDynamicEpilogueFile(const string &basename) const
{
  string filename = packageDir(basename) + "/epilogue_dynamic.m";
  ofstream output{filename, ios::out | ios::binary};
  if (!output.is_open())
    {
      cerr << "ERROR: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }

  output << "function ds = epilogue_dynamic(params, ds)" << endl
         << "% function ds = epilogue_dynamic(params, ds)" << endl
         << "% Epilogue file generated by Dynare preprocessor" << endl << endl
         << "simul_end_date = lastdate(ds);" << endl;

  deriv_node_temp_terms_t tef_terms;
  temporary_terms_t temporary_terms;
  temporary_terms_idxs_t temporary_terms_idxs;
  for (const auto & [symb_id, expr] : dynamic_def_table)
    {
      int max_lag = expr->maxLagWithDiffsExpanded();
      set<int> used_symbols;
      expr->collectVariables(SymbolType::endogenous, used_symbols);
      expr->collectVariables(SymbolType::exogenous, used_symbols);
      expr->collectVariables(SymbolType::epilogue, used_symbols);

      output << endl
             << "if ~ds.exist('" << symbol_table.getName(symb_id) << "')" << endl
             << "    ds = [ds dseries(NaN(ds.nobs,1), ds.firstdate, '" << symbol_table.getName(symb_id)<< "')];" << endl
             << "end" << endl
             << "try" << endl
             << "    simul_begin_date = firstobservedperiod(ds{";
      for (bool printed_something{false};
           int symb_id : used_symbols)
        {
          if (exchange(printed_something, true))
            output << ", ";
          output << "'" << symbol_table.getName(symb_id) << "'";
        }
      output << "}) + " << max_lag << ";" << endl
             << "    from simul_begin_date to simul_end_date do "
             << "ds." << symbol_table.getName(symb_id) << "(t) = ";
      expr->writeOutput(output, ExprNodeOutputType::epilogueFile, temporary_terms, temporary_terms_idxs, tef_terms);
      output << ";" << endl
             << "catch" << endl
             << "end" << endl;
    }
  output << "end" << endl;
  output.close();
}

void
Epilogue::writeOutput(ostream &output) const
{
  if (dynamic_def_table.empty())
    {
      output << "M_.epilogue_names = {};" << endl
             << "M_.epilogue_var_list_ = {};" << endl;
      return;
    }

  output << "M_.epilogue_names = cell(" << dynamic_def_table.size() << ",1);" << endl;
  for (int idx{1};
       const auto &[symb_id, expr] : dynamic_def_table)
    output << "M_.epilogue_names{" << idx++ << "} = '"
           << symbol_table.getName(symb_id) << "';" << endl;

  set<int> endogs;
  for (const auto & [symb_id, expr] : dynamic_def_table)
    expr->collectVariables(SymbolType::endogenous, endogs);

  vector<string> symbol_list;
  for (auto symb_id : endogs)
    symbol_list.push_back(symbol_table.getName(symb_id));
  SymbolList{move(symbol_list)}.writeOutput("M_.epilogue_var_list_", output);
}
