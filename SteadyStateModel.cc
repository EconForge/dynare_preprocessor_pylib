/*
 * Copyright (C) 2010-2014 Dynare Team
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

#include "SteadyStateModel.hh"

SteadyStateModel::SteadyStateModel(SymbolTable &symbol_table_arg, NumericalConstants &num_constants, ExternalFunctionsTable &external_functions_table_arg, const StaticModel &static_model_arg) :
  DataTree(symbol_table_arg, num_constants, external_functions_table), static_model(static_model_arg)
{
}

void
SteadyStateModel::addDefinition(int symb_id, expr_t expr)
{
  AddVariable(symb_id); // Create the variable node to be used in write method

  assert(symbol_table.getType(symb_id) == eEndogenous
         || symbol_table.getType(symb_id) == eModFileLocalVariable
         || symbol_table.getType(symb_id) == eParameter);

  // Add the variable
  vector<int> v;
  v.push_back(symb_id);
  def_table.push_back(make_pair(v, expr));
}

void
SteadyStateModel::addMultipleDefinitions(const vector<int> &symb_ids, expr_t expr)
{
  for (size_t i = 0; i < symb_ids.size(); i++)
    {
      AddVariable(symb_ids[i]); // Create the variable nodes to be used in write method
      assert(symbol_table.getType(symb_ids[i]) == eEndogenous
             || symbol_table.getType(symb_ids[i]) == eModFileLocalVariable
             || symbol_table.getType(symb_ids[i]) == eParameter);
    }
  def_table.push_back(make_pair(symb_ids, expr));
}

void
SteadyStateModel::checkPass(bool ramsey_model, WarningConsolidation &warnings) const
{
  if (def_table.size() == 0)
    return;

  vector<int> so_far_defined;

  for (size_t i = 0; i < def_table.size(); i++)
    {
      const vector<int> &symb_ids = def_table[i].first;

      // Check that symbols are not already defined
      for (size_t j = 0; j < symb_ids.size(); j++)
        if (find(so_far_defined.begin(), so_far_defined.end(), symb_ids[j])
            != so_far_defined.end())
          warnings << "WARNING: in the 'steady_state_model' block, variable '" << symbol_table.getName(symb_ids[j]) << "' is declared twice" << endl;
      
      // Check that expression has no undefined symbol
      if (!ramsey_model)
        {
          set<int> used_symbols;
          const expr_t &expr = def_table[i].second;
          expr->collectVariables(eEndogenous, used_symbols);
          expr->collectVariables(eModFileLocalVariable, used_symbols);
          for (set<int>::const_iterator it = used_symbols.begin();
               it != used_symbols.end(); ++it)
            if (find(so_far_defined.begin(), so_far_defined.end(), *it)
                == so_far_defined.end())
              {
                cerr << "ERROR: in the 'steady_state_model' block, variable '" << symbol_table.getName(*it)
                     << "' is undefined in the declaration of variable '" << symbol_table.getName(symb_ids[0]) << "'" << endl;
                exit(EXIT_FAILURE);
              }
        }

      copy(symb_ids.begin(), symb_ids.end(), back_inserter(so_far_defined));
    }

  set<int> orig_endogs = symbol_table.getOrigEndogenous();
  for (set<int>::const_iterator it = orig_endogs.begin();
       it != orig_endogs.end(); ++it)
    {
      if (find(so_far_defined.begin(), so_far_defined.end(), *it)
          == so_far_defined.end())
        warnings << "WARNING: in the 'steady_state_model' block, variable '" << symbol_table.getName(*it) << "' is not assigned a value" << endl;
    }
}

void
SteadyStateModel::writeSteadyStateFile(const string &basename, bool ramsey_model) const
{
  if (def_table.size() == 0)
    return;

  string filename = basename + "_steadystate2.m";

  ofstream output;
  output.open(filename.c_str(), ios::out | ios::binary);
  if (!output.is_open())
    {
      cerr << "ERROR: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }

  output << "function [ys_, params, info] = " << basename << "_steadystate2("
	 << "ys_, exo_, params)" << endl
         << "% Steady state generated by Dynare preprocessor" << endl
         << "    info = 0;" << endl;

  for (size_t i = 0; i < def_table.size(); i++)
    {
      const vector<int> &symb_ids = def_table[i].first;
      output << "    ";
      if (symb_ids.size() > 1)
        output << "[";
      for (size_t j = 0; j < symb_ids.size(); j++)
        {
          variable_node_map_t::const_iterator it = variable_node_map.find(make_pair(symb_ids[j], 0));
          assert(it != variable_node_map.end());
          dynamic_cast<ExprNode *>(it->second)->writeOutput(output, oSteadyStateFile);
          if (j < symb_ids.size()-1)
            output << ",";
        }
      if (symb_ids.size() > 1)
        output << "]";

      output << "=";
      def_table[i].second->writeOutput(output, oSteadyStateFile);
      output << ";" << endl;
    }
  output << "    % Auxiliary equations" << endl;
  static_model.writeAuxVarInitval(output, oSteadyStateFile);
  output << "    check_=0;" << endl
         << "end" << endl;
}

void
SteadyStateModel::writeSteadyStateFileC(const string &basename, bool ramsey_model, bool cuda) const
{
  string filename = basename + "_steadystate.c";

  ofstream output;
  output.open(filename.c_str(), ios::out | ios::binary);
  if (!output.is_open())
    {
      cerr << "ERROR: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }

  output << "#include <math.h>" << endl;

  if (cuda)
    output << "__global__ ";

  output << "void steadystate("
	 << "const double *exo_, const double *params, double *ys_, int *info)" << endl
         << "// Steady state file generated by Dynare preprocessor" << endl
	 << "{" << endl
         << "    *info = 0;" << endl;

  if (def_table.size() == 0)
    {
      output << "    return;" << endl
	     << "}" << endl;
      return;
    }

  for (size_t i = 0; i < def_table.size(); i++)
    {
      const vector<int> &symb_ids = def_table[i].first;
      output << "    ";
      if (symb_ids.size() > 1)
	std::cout << "Error: in C, multiple returns are not permitted in steady_state_model" << std::endl;
      variable_node_map_t::const_iterator it = variable_node_map.find(make_pair(symb_ids[0], 0));
      assert(it != variable_node_map.end());
      if (it->second->get_type() == eModFileLocalVariable)
	output << "double ";
      dynamic_cast<ExprNode *>(it->second)->writeOutput(output, oCSteadyStateFile);
      output << "=";
      def_table[i].second->writeOutput(output, oCSteadyStateFile);
      output << ";" << endl;
    }
  output << "    // Auxiliary equations" << endl;
  static_model.writeAuxVarInitval(output, oCSteadyStateFile);
  output << "}" << endl;
}

