/*
 * Copyright © 2010-2022 Dynare Team
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

#ifndef _MODEL_EQUATION_BLOCK_HH
#define _MODEL_EQUATION_BLOCK_HH

#include "DataTree.hh"
#include "Statement.hh"
#include "StaticModel.hh"
#include "DynamicModel.hh"
#include "WarningConsolidation.hh"

class PlannerObjective : public StaticModel
{
public:
  PlannerObjective(SymbolTable &symbol_table_arg,
                   NumericalConstants &num_constants_arg,
                   ExternalFunctionsTable &external_functions_table_arg);
protected:
  string
  modelClassName() const override
  {
    return "planner objective";
  }

private:
  void computingPassBlock(const eval_context_t &eval_context, bool no_tmp_terms) override;
};

class OrigRamseyDynamicModel : public DynamicModel
{
public:
  OrigRamseyDynamicModel(SymbolTable &symbol_table_arg,
                         NumericalConstants &num_constants_arg,
                         ExternalFunctionsTable &external_functions_table_arg,
                         TrendComponentModelTable &trend_component_model_table_arg,
                         VarModelTable &var_model_table_arg);
  OrigRamseyDynamicModel &operator=(const DynamicModel &m);

protected:
  string
  modelClassName() const override
  {
    return "original Ramsey model";
  }
};

class SteadyStateModel : public DataTree
{
private:
  //! Associates a set of symbol IDs (the variable(s) assigned in a given statement) to an expression (their assigned value)
  vector<pair<vector<int>, expr_t>> def_table;

  //! Reference to static model (for writing auxiliary equations)
  const StaticModel &static_model;

public:
  SteadyStateModel(SymbolTable &symbol_table_arg,
                   NumericalConstants &num_constants_arg,
                   ExternalFunctionsTable &external_functions_table_arg,
                   const StaticModel &static_model_arg);

  SteadyStateModel(const SteadyStateModel &m);
  SteadyStateModel &operator=(const SteadyStateModel &m);

  //! Add an expression of the form "var = expr;"
  void addDefinition(int symb_id, expr_t expr);
  //! Add an expression of the form "[ var1, var2, ... ] = expr;"
  void addMultipleDefinitions(const vector<int> &symb_ids, expr_t expr);
  //! Checks that definitions are in a recursive order, and that no variable is declared twice
  /*!
    \param[in] ramsey_model Is there a Ramsey model in the MOD file? If yes, then disable the check on the recursivity of the declarations
  */
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) const;
  //! Write the steady state file
  void writeSteadyStateFile(const string &basename, bool julia) const;
  //! Writes LaTeX file with the equations of the dynamic model (for the steady state model)
  void writeLatexSteadyStateFile(const string &basename) const;
  //! Writes JSON output
  void writeJsonSteadyStateFile(ostream &output, bool transformComputingPass) const;
};

class Epilogue : public DynamicModel
{
private:
  //! Associates a symbol ID (the variable assigned in a given statement) to an expression (its assigned value)
  vector<pair<int, expr_t>> dynamic_def_table, static_def_table;
public:
  Epilogue(SymbolTable &symbol_table_arg,
           NumericalConstants &num_constants_arg,
           ExternalFunctionsTable &external_functions_table_arg,
           TrendComponentModelTable &trend_component_model_table_arg,
           VarModelTable &var_model_table_arg);

  Epilogue(const Epilogue &m);
  Epilogue &operator=(const Epilogue &m);

  //! Add an expression of the form "var = expr;"
  void addDefinition(int symb_id, expr_t expr);

  //! Checks that no variable is declared twice, and that “with_epilogue” is not misused
  void checkPass(ModFileStructure &mod_file_struct) const;

  //! Creates static epilogue equations
  void toStatic();

  //! Deal with trend variables in the epilogue block
  void detrend(const map<int, expr_t> &trend_symbols_map,
               const nonstationary_symbols_map_t &nonstationary_symbols_map);

  //! Write the steady state file
  void writeEpilogueFile(const string &basename) const;

  //! Write Output
  void writeOutput(ostream &output) const;

protected:
  string
  modelClassName() const override
  {
    return "epilogue";
  }

private:
  //! Helper for public writeEpilogueFile
  void writeStaticEpilogueFile(const string &basename) const;
  void writeDynamicEpilogueFile(const string &basename) const;
  void computingPassBlock(const eval_context_t &eval_context, bool no_tmp_terms) override;
};

#endif
