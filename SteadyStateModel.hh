/*
 * Copyright (C) 2010 Dynare Team
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

class SteadyStateModel : public DataTree
{
private:
  //! Associates a symbol ID to an expression of the form "var = expr"
  map<int, NodeID> def_table;
  vector<int> recursive_order;

  //! Reference to static model (for writing auxiliary equations)
  const StaticModel &static_model;

public:
  class AlreadyDefinedException
  {
  public:
    const string &varname;
    AlreadyDefinedException(const string &varname_arg) : varname(varname_arg) {}
  };
  class UndefinedVariableException
  {
  public:
    const string &varname;
    UndefinedVariableException(const string &varname_arg) : varname(varname_arg) {}
  };

  SteadyStateModel(SymbolTable &symbol_table_arg, NumericalConstants &num_constants, ExternalFunctionsTable &external_functions_table_arg, const StaticModel &static_model_arg);
  //! Add an expression of the form "var = expr;"
  void addDefinition(int symb_id, NodeID expr) throw (UndefinedVariableException, AlreadyDefinedException);
  //! Write the steady state file
  void writeSteadyStateFile(const string &basename) const;
};

#endif
