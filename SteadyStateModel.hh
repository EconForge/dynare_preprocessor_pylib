/*
 * Copyright (C) 2010-2015 Dynare Team
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

#ifndef _STEADY_STATE_MODEL_HH
#define _STEADY_STATE_MODEL_HH

#include "DataTree.hh"
#include "StaticModel.hh"
#include "WarningConsolidation.hh"

class SteadyStateModel : public DataTree
{
private:
  //! Associates a set of symbol IDs (the variable(s) assigned in a given statement) to an expression (their assigned value)
  vector<pair<vector<int>, expr_t> > def_table;

  //! Reference to static model (for writing auxiliary equations)
  const StaticModel &static_model;

public:
  SteadyStateModel(SymbolTable &symbol_table_arg, NumericalConstants &num_constants, ExternalFunctionsTable &external_functions_table_arg, const StaticModel &static_model_arg);
  //! Add an expression of the form "var = expr;"
  void addDefinition(int symb_id, expr_t expr);
  //! Add an expression of the form "[ var1, var2, ... ] = expr;"
  void addMultipleDefinitions(const vector<int> &symb_ids, expr_t expr);
  //! Checks that definitions are in a recursive order, and that no variable is declared twice
  /*!
    \param[in] ramsey_model Is there a Ramsey model in the MOD file? If yes, then disable the check on the recursivity of the declarations
  */
  void checkPass(bool ramsey_model, WarningConsolidation &warnings) const;
  //! Write the steady state file
  /*!
    \param[in] ramsey_model Is there a Ramsey model in the MOD file? If yes, then use the "ys" in argument of the steady state file as initial values
  */
  void writeSteadyStateFile(const string &basename, bool ramsey_model) const;
  // in ExternalFiles.cc
  void writeSteadyStateFileC(const string &basename, bool ramsey_model) const;
  //! Used to reindex def_table after removal of extra exogenous
  void reindex(SymbolTable orig_symbol_table);
};

#endif
