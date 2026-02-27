/*
 * Copyright © 2010-2026 Dynare Team
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
#include <ranges>
#include <sstream>

#include "ModelEquationBlock.hh"
#include "Utils.hh"

PlannerObjective::PlannerObjective(SymbolTable& symbol_table_arg,
                                   NumericalConstants& num_constants_arg,
                                   ExternalFunctionsTable& external_functions_table_arg,
                                   HeterogeneityTable& heterogeneity_table_arg,
                                   DatabaseTable& database_table_arg) :
    StaticModel {symbol_table_arg, num_constants_arg, external_functions_table_arg,
                 heterogeneity_table_arg, database_table_arg}
{
}

void
PlannerObjective::writeDriverOutput(ostream& output) const
{
  output << "M_.objective_tmp_nbr = [";
  for (const auto& it : temporary_terms_derivatives)
    output << it.size() << "; ";
  output << "];" << '\n';
  writeDriverSparseIndicesHelper("objective", output);
}

void
PlannerObjective::computingPassBlock([[maybe_unused]] const eval_context_t& eval_context,
                                     [[maybe_unused]] bool no_tmp_terms)
{
  // Disable block decomposition on planner objective
}

OrigRamseyDynamicModel::OrigRamseyDynamicModel(
    SymbolTable& symbol_table_arg, NumericalConstants& num_constants_arg,
    ExternalFunctionsTable& external_functions_table_arg,
    HeterogeneityTable& heterogeneity_table_arg, DatabaseTable& database_table_arg,
    TrendComponentModelTable& trend_component_model_table_arg, VarModelTable& var_model_table_arg) :
    DynamicModel {symbol_table_arg,        num_constants_arg,  external_functions_table_arg,
                  heterogeneity_table_arg, database_table_arg, trend_component_model_table_arg,
                  var_model_table_arg}
{
}

OrigRamseyDynamicModel&
OrigRamseyDynamicModel::operator=(const DynamicModel& m)
{
  DynamicModel::operator=(m);
  return *this;
}

SteadyStateModel::SteadyStateModel(SymbolTable& symbol_table_arg,
                                   NumericalConstants& num_constants_arg,
                                   ExternalFunctionsTable& external_functions_table_arg,
                                   HeterogeneityTable& heterogeneity_table_arg,
                                   DatabaseTable& database_table_arg,
                                   const StaticModel& static_model_arg) :
    DataTree {symbol_table_arg, num_constants_arg, external_functions_table_arg,
              heterogeneity_table_arg, database_table_arg},
    static_model {static_model_arg}
{
}

SteadyStateModel::SteadyStateModel(const SteadyStateModel& m) :
    DataTree {m}, static_model {m.static_model}
{
  assert(m.def_table.size() == m.def_table_lineno.size());
  for (size_t i {0}; i < m.def_table.size(); i++)
    {
      const auto& [ids, expr] = m.def_table[i];
      def_table.emplace_back(ids, expr->clone(*this));
      def_table_lineno.push_back(m.def_table_lineno[i]);
    }
}

SteadyStateModel&
SteadyStateModel::operator=(const SteadyStateModel& m)
{
  DataTree::operator=(m);

  assert(&static_model == &m.static_model);

  def_table.clear();
  def_table_lineno.clear();

  assert(m.def_table.size() == m.def_table_lineno.size());
  for (size_t i {0}; i < m.def_table.size(); i++)
    {
      const auto& [ids, expr] = m.def_table[i];
      def_table.emplace_back(ids, expr->clone(*this));
      def_table_lineno.push_back(m.def_table_lineno[i]);
    }

  return *this;
}

void
SteadyStateModel::addDefinition(int symb_id, expr_t expr, optional<int> lineno)
{
  AddVariable(symb_id); // Create the variable node to be used in write method

  assert(symbol_table.getType(symb_id) == SymbolType::endogenous
         || symbol_table.getType(symb_id) == SymbolType::modFileLocalVariable
         || symbol_table.getType(symb_id) == SymbolType::parameter);

  // Add the variable
  vector v {symb_id};
  def_table.emplace_back(v, expr);
  def_table_lineno.push_back(lineno);
}

void
SteadyStateModel::addMultipleDefinitions(const vector<int>& symb_ids, expr_t expr,
                                         optional<int> lineno)
{
  for (int symb_id : symb_ids)
    {
      AddVariable(symb_id); // Create the variable nodes to be used in write method
      assert(symbol_table.getType(symb_id) == SymbolType::endogenous
             || symbol_table.getType(symb_id) == SymbolType::modFileLocalVariable
             || symbol_table.getType(symb_id) == SymbolType::parameter);
    }
  def_table.emplace_back(symb_ids, expr);
  def_table_lineno.push_back(lineno);
}

void
SteadyStateModel::checkPass(ModFileStructure& mod_file_struct, WarningConsolidation& warnings) const
{
  if (def_table.size() == 0)
    return;

  assert(def_table.size() == def_table_lineno.size());

  mod_file_struct.steady_state_model_present = true;
  set<int> so_far_defined;

  for (size_t def_idx {0}; def_idx < def_table.size(); def_idx++)
    {
      const auto& [symb_ids, expr] = def_table[def_idx];
      const auto& lineno = def_table_lineno[def_idx];

      // Check that symbols are not already defined
      for (int symb_id : symb_ids)
        if (so_far_defined.contains(symb_id))
          warnings << "WARNING: in the 'steady_state_model' block, variable '"
                   << symbol_table.getName(symb_id) << "' is declared twice" << '\n';

      // Check that expression has no undefined symbol
      if (!mod_file_struct.ramsey_model_present)
        {
          set<int> used_symbols;
          expr->collectVariables(SymbolType::endogenous, used_symbols);
          expr->collectVariables(SymbolType::modFileLocalVariable, used_symbols);
          for (int used_symbol : used_symbols)
            if (!so_far_defined.contains(used_symbol))
              {
                cerr << "ERROR: in the 'steady_state_model' block";
                if (lineno)
                  cerr << ", line " << *lineno;
                cerr << ", variable '" << symbol_table.getName(used_symbol)
                     << "' is undefined in the declaration of variable '"
                     << symbol_table.getName(symb_ids[0]) << "'" << '\n';
                exit(EXIT_FAILURE);
              }
        }

      so_far_defined.insert(symb_ids.begin(), symb_ids.end());
    }

  /* Check that all original endogous are defined (except the instruments of a
     Ramsey model, since the steady_state_block should give the steady state
     *conditional* to those instruments) */
  set<int> should_be_defined = symbol_table.getOrigEndogenous();
  if (mod_file_struct.ramsey_model_present)
    for (const auto& s : mod_file_struct.instruments.getSymbols())
      should_be_defined.erase(symbol_table.getID(s));
  for (int v : should_be_defined)
    if (!so_far_defined.contains(v))
      warnings << "WARNING: in the 'steady_state_model' block, variable '"
               << symbol_table.getName(v) << "' is not assigned a value" << '\n';
}

void
SteadyStateModel::writeLatexSteadyStateFile(const string& basename) const
{
  filesystem::create_directories(basename + "/latex");

  const filesystem::path filename {basename + "/latex/steady_state.tex"},
      content_filename {basename + "/latex/steady_state_content.tex"};

  ofstream output {filename, ios::out | ios::binary};
  if (!output.is_open())
    {
      cerr << "ERROR: Can't open file " << filename.string() << " for writing" << '\n';
      exit(EXIT_FAILURE);
    }

  ofstream content_output {content_filename, ios::out | ios::binary};
  if (!content_output.is_open())
    {
      cerr << "ERROR: Can't open file " << content_filename.string() << " for writing" << '\n';
      exit(EXIT_FAILURE);
    }

  output << "\\documentclass[10pt,a4paper]{article}" << '\n'
         << "\\usepackage[landscape]{geometry}" << '\n'
         << "\\usepackage{fullpage}" << '\n'
         << "\\usepackage{amsfonts}" << '\n'
         << "\\usepackage{breqn}" << '\n'
         << "\\begin{document}" << '\n'
         << "\\footnotesize" << '\n';

  for (const auto& [ids, value] : def_table)
    for (int id : ids)
      {
        content_output << "\\begin{dmath}" << '\n' << symbol_table.getTeXName(id) << " = ";
        value->writeOutput(content_output, ExprNodeOutputType::latexStaticModel);
        content_output << '\n' << "\\end{dmath}" << '\n';
      }

  static_model.writeLatexAuxVarRecursiveDefinitions(content_output);

  output << "\\include{steady_state_content.tex}" << '\n' << "\\end{document}" << '\n';

  output.close();
  content_output.close();
}

void
SteadyStateModel::writeSteadyStateFile(const string& basename, bool julia) const
{
  if (def_table.size() == 0)
    return;

  ExprNodeOutputType output_type
      = (julia ? ExprNodeOutputType::juliaSteadyStateFile : ExprNodeOutputType::steadyStateFile);

  stringstream output;
  if (!julia)
    output << "function [ys_, params, info] = steadystate("
           << "ys_, exo_, params)" << '\n'
           << "% Steady state generated by Dynare preprocessor" << '\n'
           << "    info = 0;" << '\n';
  else
    output << "# NB: this file was automatically generated by Dynare" << '\n'
           << "#     from " << basename << ".mod" << '\n'
           << "#" << '\n'
           << "function steady_state!(ys_::Vector{<: Real}, exo_::Vector{<: Real}, "
           << "params::Vector{<: Real})" << '\n'
           << "@inbounds begin" << '\n';

  for (const auto& [symb_ids, value] : def_table)
    {
      output << "    ";
      if (symb_ids.size() > 1)
        output << "[";
      for (size_t j = 0; j < symb_ids.size(); j++)
        {
          getVariable(symb_ids[j])->ExprNode::writeOutput(output, output_type);
          if (j < symb_ids.size() - 1)
            output << ",";
        }
      if (symb_ids.size() > 1)
        output << "]";

      output << "=";
      value->writeOutput(output, output_type);
      output << ";" << '\n';
    }
  if (!julia)
    output << "    % Auxiliary equations" << '\n';
  else
    output << "    # Auxiliary equations" << '\n';
  static_model.writeAuxVarRecursiveDefinitions(output, output_type);

  output << "end" << '\n';
  if (julia)
    output << "end" << '\n';

  if (julia)
    writeToFileIfModified(output,
                          filesystem::path {basename} / "model" / "julia" / "SteadyState2.jl");
  else
    {
      /* Calling writeToFileIfModified() is useless here since we write inside
         a subdirectory deleted at each preprocessor run. */
      filesystem::path filename {packageDir(basename) / "steadystate.m"};
      ofstream output_file {filename, ios::out | ios::binary};
      if (!output_file.is_open())
        {
          cerr << "ERROR: Can't open file " << filename.string() << " for writing" << '\n';
          exit(EXIT_FAILURE);
        }
      output_file << output.str();
      output_file.close();
    }
}

void
SteadyStateModel::writeJsonSteadyStateFile(ostream& output, bool transformComputingPass) const
{
  if (def_table.size() == 0)
    return;

  vector<pair<string, string>> eqtags;

  output << "{\"steady_state_model\": [";

  for (bool printed_something {false}; const auto& [symb_ids, value] : def_table)
    {
      if (exchange(printed_something, true))
        output << ",";
      output << "{\"lhs\": ";
      if (symb_ids.size() > 1)
        output << "[";
      for (bool printed_something2 {false}; int symb_id : symb_ids)
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
      output << "\"}" << '\n';
    }

  if (transformComputingPass)
    static_model.writeJsonAuxVarRecursiveDefinitions(output);

  output << "]}";
}

set<int>
SteadyStateModel::getUsedParameters() const
{
  set<int> used;

  for (const auto& [symb_ids, expr] : def_table)
    {
      for (int symb_id : symb_ids)
        if (auto type = symbol_table.getType(symb_id);
            type == SymbolType::parameter || type == SymbolType::heterogeneousParameter)
          used.insert(symb_id);

      expr->collectVariables(SymbolType::parameter, used);
      expr->collectVariables(SymbolType::heterogeneousParameter, used);
    }

  return used;
}

Epilogue::Epilogue(SymbolTable& symbol_table_arg, NumericalConstants& num_constants_arg,
                   ExternalFunctionsTable& external_functions_table_arg,
                   HeterogeneityTable& heterogeneity_table_arg, DatabaseTable& database_table_arg,
                   TrendComponentModelTable& trend_component_model_table_arg,
                   VarModelTable& var_model_table_arg) :
    DynamicModel {symbol_table_arg,        num_constants_arg,  external_functions_table_arg,
                  heterogeneity_table_arg, database_table_arg, trend_component_model_table_arg,
                  var_model_table_arg}
{
}

Epilogue::Epilogue(const Epilogue& m) : DynamicModel {m}
{
  for (const auto& it : m.dynamic_def_table)
    dynamic_def_table.emplace_back(it.first, it.second->clone(*this));
}

Epilogue&
Epilogue::operator=(const Epilogue& m)
{
  DynamicModel::operator=(m);

  dynamic_def_table.clear();
  for (const auto& it : m.dynamic_def_table)
    dynamic_def_table.emplace_back(it.first, it.second->clone(*this));

  return *this;
}

void
Epilogue::addDefinition(int symb_id, expr_t expr)
{
  dynamic_def_table.emplace_back(symb_id, expr);
}

void
Epilogue::checkPass(ModFileStructure& mod_file_struct) const
{
  if (dynamic_def_table.size() == 0)
    {
      if (mod_file_struct.with_epilogue_option)
        {
          cerr << "ERROR: the 'with_epilogue' option cannot be specified when there is no "
                  "'epilogue' block"
               << '\n';
          exit(EXIT_FAILURE);
        }
      return;
    }

  set<int> so_far_defined;
  for (const auto& [symb_id, expr] : dynamic_def_table)
    if (so_far_defined.contains(symb_id))
      {
        cerr << "ERROR: in the 'epilogue' block, variable '" << symbol_table.getName(symb_id)
             << "' is declared twice" << '\n';
        exit(EXIT_FAILURE);
      }
    else
      so_far_defined.insert(symb_id);
}

void
Epilogue::toStatic()
{
  for (const auto& [symb_id, expr] : dynamic_def_table)
    static_def_table.emplace_back(symb_id, expr->toStatic(*this));
}

void
Epilogue::detrend(const map<int, expr_t>& trend_symbols_map,
                  const nonstationary_symbols_map_t& nonstationary_symbols_map)
{
  for (const auto& [symb_id, deflator] : ranges::reverse_view(nonstationary_symbols_map))
    for (auto& [symb_id, expr] : dynamic_def_table)
      {
        expr = expr->detrend(symb_id, deflator.first, deflator.second);
        assert(expr);
      }

  for (auto& [symb_id, expr] : dynamic_def_table)
    {
      expr = expr->removeTrendLeadLag(trend_symbols_map);
      assert(expr);
    }

  for (auto& [symb_id, expr] : dynamic_def_table)
    {
      expr = expr->replaceTrendVar();
      assert(expr);
    }
}

void
Epilogue::writeEpilogueFile(const string& basename) const
{
  if (dynamic_def_table.empty())
    return;

  writeDynamicEpilogueFile(basename);
  writeStaticEpilogueFile(basename);
}

void
Epilogue::writeStaticEpilogueFile(const string& basename) const
{
  filesystem::path filename {packageDir(basename) / "epilogue_static.m"};
  ofstream output {filename, ios::out | ios::binary};
  if (!output.is_open())
    {
      cerr << "ERROR: Can't open file " << filename.string() << " for writing" << '\n';
      exit(EXIT_FAILURE);
    }

  output << "function ds = epilogue_static(params, ds)" << '\n'
         << "% function ds = epilogue_static(params, ds)" << '\n'
         << "% Epilogue file generated by Dynare preprocessor" << '\n';

  for (const auto& [symb_id, expr] : static_def_table)
    {
      // Rewrite external function TEF term for every equation as argument values could have been
      // changed in between two calls to the same function;
      deriv_node_temp_terms_t tef_terms;
      temporary_terms_t temporary_terms;
      temporary_terms_idxs_t temporary_terms_idxs;
      output << '\n';
      if (expr->containsExternalFunction())
        expr->writeExternalFunctionOutput(output, ExprNodeOutputType::matlabDseries,
                                          temporary_terms, temporary_terms_idxs, tef_terms);
      output << "epilogue_static_tmp_term = ";
      expr->writeOutput(output, ExprNodeOutputType::matlabDseries, temporary_terms,
                        temporary_terms_idxs, tef_terms);
      output << ";" << '\n'
             << "if isdseries(epilogue_static_tmp_term)" << '\n'
             << "    ds." << symbol_table.getName(symb_id) << " = epilogue_static_tmp_term;" << '\n'
             << "else" << '\n'
             << "    ds." << symbol_table.getName(symb_id)
             << " = dseries(ones(ds.nobs,1)*epilogue_static_tmp_term, ds.firstdate, '"
             << symbol_table.getName(symb_id) << "');" << '\n'
             << "end" << '\n';
    }
  output << "end" << '\n';
  output.close();
}

void
Epilogue::writeDynamicEpilogueFile(const string& basename) const
{
  filesystem::path filename {packageDir(basename) / "epilogue_dynamic.m"};
  ofstream output {filename, ios::out | ios::binary};
  if (!output.is_open())
    {
      cerr << "ERROR: Can't open file " << filename.string() << " for writing" << '\n';
      exit(EXIT_FAILURE);
    }

  output << "function ds = epilogue_dynamic(params, ds)" << '\n'
         << "% function ds = epilogue_dynamic(params, ds)" << '\n'
         << "% Epilogue file generated by Dynare preprocessor" << '\n'
         << '\n'
         << "simul_end_date = lastdate(ds);" << '\n';

  deriv_node_temp_terms_t tef_terms;
  temporary_terms_t temporary_terms;
  temporary_terms_idxs_t temporary_terms_idxs;
  for (const auto& [symb_id, expr] : dynamic_def_table)
    {
      int max_lag = expr->maxLagWithDiffsExpanded();
      set<int> used_symbols;
      expr->collectVariables(SymbolType::endogenous, used_symbols);
      expr->collectVariables(SymbolType::exogenous, used_symbols);
      expr->collectVariables(SymbolType::epilogue, used_symbols);

      output << '\n'
             << "if ~ds.exist('" << symbol_table.getName(symb_id) << "')" << '\n'
             << "    ds = [ds dseries(NaN(ds.nobs,1), ds.firstdate, '"
             << symbol_table.getName(symb_id) << "')];" << '\n'
             << "end" << '\n'
             << "try" << '\n'
             << "    simul_begin_date = firstobservedperiod(ds{";
      for (bool printed_something {false}; int symb_id : used_symbols)
        {
          if (exchange(printed_something, true))
            output << ", ";
          output << "'" << symbol_table.getName(symb_id) << "'";
        }
      output << "}) + " << max_lag << ";" << '\n'
             << "    from simul_begin_date to simul_end_date do "
             << "ds." << symbol_table.getName(symb_id) << "(t) = ";
      expr->writeOutput(output, ExprNodeOutputType::epilogueFile, temporary_terms,
                        temporary_terms_idxs, tef_terms);
      output << ";" << '\n' << "catch" << '\n' << "end" << '\n';
    }
  output << "end" << '\n';
  output.close();
}

void
Epilogue::writeOutput(ostream& output) const
{
  if (dynamic_def_table.empty())
    {
      output << "M_.epilogue_names = {};" << '\n' << "M_.epilogue_var_list_ = {};" << '\n';
      return;
    }

  output << "M_.epilogue_names = cell(" << dynamic_def_table.size() << ",1);" << '\n';
  for (int idx {1}; const auto& [symb_id, expr] : dynamic_def_table)
    output << "M_.epilogue_names{" << idx++ << "} = '" << symbol_table.getName(symb_id) << "';"
           << '\n';

  set<int> endogs;
  for (const auto& [symb_id, expr] : dynamic_def_table)
    expr->collectVariables(SymbolType::endogenous, endogs);

  vector<string> symbol_list;
  symbol_list.reserve(endogs.size());
  for (auto symb_id : endogs)
    symbol_list.push_back(symbol_table.getName(symb_id));
  SymbolList {move(symbol_list)}.writeOutput("M_.epilogue_var_list_", output);
}

void
Epilogue::computingPassBlock([[maybe_unused]] const eval_context_t& eval_context,
                             [[maybe_unused]] bool no_tmp_terms)
{
  // Disable block decomposition on epilogue blocks
}

void
SteadyStateModel::writePythonSteadyStateFile(const string& basename, bool use_jax,
                                             bool use_numba) const
{
  filesystem::path model_dir {basename};
  model_dir /= "python";
  filesystem::create_directories(model_dir);

  string filename = model_dir.string() + "/steadystate.py";
  ofstream output {filename, ios::out | ios::binary};
  if (!output.is_open())
    {
      cerr << "ERROR: Can't open file " << filename << " for writing" << '\n';
      exit(EXIT_FAILURE);
    }

  output << writePythonHelpersImports(use_jax, use_numba) << '\n'
         << writePythonHelpersImportFromInit() << '\n'
         << '\n';

  output << "def steadystate(y, exo, params):" << '\n';
  if (def_table.size() == 0)
    {
      output << "    return y, params, 0" << '\n';
      output.close();
      return;
    }

  temporary_terms_t temp_terms_union;
  deriv_node_temp_terms_t tef_terms;
  temporary_terms_idxs_t temporary_terms_idxs;

  for (const auto& [symb_ids, value] : def_table)
    {
      output << "    ";
      if (symb_ids.size() > 1)
        {
          for (size_t j = 0; j < symb_ids.size(); j++)
            {
              getVariable(symb_ids[j])
                  ->writeOutput(output, ExprNodeOutputType::pythonStaticModel, temp_terms_union,
                                temporary_terms_idxs, tef_terms);
              if (j < symb_ids.size() - 1)
                output << ", ";
            }
          output << " = ";
        }
      else
        {
          getVariable(symb_ids[0])
              ->writeOutput(output, ExprNodeOutputType::pythonStaticModel, temp_terms_union,
                            temporary_terms_idxs, tef_terms);
          output << " = ";
        }

      value->writeOutput(output, ExprNodeOutputType::pythonStaticModel, temp_terms_union,
                         temporary_terms_idxs, tef_terms);
      output << '\n';
    }

  output << "    # Auxiliary equations" << '\n';
  static_model.writeAuxVarRecursiveDefinitions(output, ExprNodeOutputType::pythonStaticModel);

  output << "    return y, params, 0" << '\n';
  output.close();
}
