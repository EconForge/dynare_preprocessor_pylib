/*
 * Copyright (C) 2010-2018 Dynare Team
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
    def_table.push_back(make_pair(it.first, it.second->cloneDynamic(*this)));
}

SteadyStateModel &
SteadyStateModel::operator=(const SteadyStateModel &m)
{
  DataTree::operator=(m);

  assert(&static_model == &m.static_model);

  def_table.clear();
  for (const auto &it : m.def_table)
    def_table.push_back(make_pair(it.first, it.second->cloneDynamic(*this)));

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
  ofstream output, content_output;
  string filename = basename + "_steady_state.tex";
  string content_basename = basename + "_steady_state_content";
  string content_filename = content_basename + ".tex";

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

  output << "\\include{" << content_basename << "}" << endl
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
          auto it = variable_node_map.find({ symb_ids[j], 0 });
          assert(it != variable_node_map.end());
          dynamic_cast<ExprNode *>(it->second)->writeOutput(output, output_type);
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

  if (!julia)
    output << "    check_=0;" << endl;

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
          auto it =
            variable_node_map.find({ symb_ids[j], 0 });
          assert(it != variable_node_map.end());
          output << "\"";
          dynamic_cast<ExprNode *>(it->second)->writeJsonOutput(output, {}, {}, false);
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
  DynamicModel {m},
  endogs {m.endogs},
  exogs {m.exogs}
{
  for (const auto &it : m.def_table)
    def_table.push_back(make_pair(it.first, it.second->cloneDynamic(*this)));
}

Epilogue &
Epilogue::operator=(const Epilogue &m)
{
  DynamicModel::operator=(m);

  endogs = m.endogs;
  exogs = m.exogs;

  def_table.clear();
  for (const auto &it : m.def_table)
    def_table.push_back(make_pair(it.first, it.second->cloneDynamic(*this)));

  return *this;
}

void
Epilogue::addDefinition(int symb_id, expr_t expr)
{
  AddVariable(symb_id); // Create the variable node to be used in write method
  def_table.emplace_back(symb_id, expr);
  endogs.emplace(symb_id);
  expr->collectVariables(SymbolType::endogenous, exogs);
  expr->collectVariables(SymbolType::exogenous, exogs);
  expr->collectVariables(SymbolType::endogenousEpilogue, exogs);
  expr->collectVariables(SymbolType::exogenousEpilogue, exogs);
  for (auto it : endogs)
    exogs.erase(it);
}

void
Epilogue::checkPass(WarningConsolidation &warnings) const
{
  if (def_table.size() == 0)
    return;

  vector<int> so_far_defined;
  for (const auto & it : def_table)
    {
      if (find(so_far_defined.begin(), so_far_defined.end(), it.first) != so_far_defined.end())
        {
          cerr << "WARNING: in the 'epilogue' block, variable '" << symbol_table.getName(it.first)
               << "' is declared twice" << endl;
          exit(EXIT_FAILURE);
        }
      so_far_defined.push_back(it.first);
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
  output << "function ds = epilogue(params, ds)" << endl
         << "% function ds = epilogue(params, ds)" << endl
         << "% Epilogue file generated by Dynare preprocessor" << endl << endl
         << "epilogue_ds_first_date__ = ds.firstdate;" << endl
         << "epilogue_loop_begin_idx__ = lastdate(ds) - ds.lastobservedperiod;" << endl
         << "epilogue_loop_end_idx__ = lastdate(ds) - firstdate(ds) + 1;" << endl << endl;

  output << "% endogenous" << endl;
  for (auto symb_id : endogs)
    output << symbol_table.getName(symb_id) << " = ds." << symbol_table.getName(symb_id) << ".data;" << endl;
  output << endl
         << "% exogenous" << endl;
  for (auto symb_id : exogs)
    output << symbol_table.getName(symb_id) << " = ds." << symbol_table.getName(symb_id) << ".data;" << endl;
  output << endl
         << "for epilogue_it__ = epilogue_loop_begin_idx__::epilogue_loop_end_idx__" << endl;

  deriv_node_temp_terms_t tef_terms;
  temporary_terms_t temporary_terms;
  temporary_terms_idxs_t temporary_terms_idxs;
  for (const auto & it : def_table)
    if (it.second->containsExternalFunction())
      {
        output << "    ";
        it.second->writeExternalFunctionOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
      }
  output << endl;
  for (const auto & it : def_table)
    {
      auto node = variable_node_map.find({ it.first, 0 });
      assert(node != variable_node_map.end());

      output << "    ";
      dynamic_cast<ExprNode *>(node->second)->writeOutput(output, output_type);
      output << " = ";
      it.second->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
      output << ";" << endl;
    }
  output << "end" << endl << endl;
  for (auto symb_id : endogs)
    output << "ds." << symbol_table.getName(symb_id) << " = dseries(" << symbol_table.getName(symb_id)
           << ", epilogue_ds_first_date__);" << endl;
  output << endl
         << "end" << endl;
  output.close();
}
