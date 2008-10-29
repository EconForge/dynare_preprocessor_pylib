/*
 * Copyright (C) 2003-2008 Dynare Team
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

#include <cstdlib>

#include "VariableTable.hh"

VariableTable::VariableTable(const SymbolTable &symbol_table_arg) :
  symbol_table(symbol_table_arg),
  var_endo_nbr(0), var_exo_nbr(0), var_exo_det_nbr(0),
  max_lag(0), max_lead(0),
  max_endo_lag(0), max_endo_lead(0),
  max_exo_lag(0), max_exo_lead(0),
  max_exo_det_lag(0), max_exo_det_lead(0)
{
}

int
VariableTable::addVariable(SymbolType type, int symb_id, int lag) throw (DynJacobianColsAlreadyComputedException)
{
  if (dyn_jacobian_cols_table.size() != 0)
    throw DynJacobianColsAlreadyComputedException();

  var_key_type key = make_pair(make_pair(type, lag), symb_id);

  // Testing if variable already exists
  variable_table_type::const_iterator it = variable_table.find(key);
  if (it != variable_table.end())
    return it->second;

  int var_id = size();

  variable_table[key] = var_id;
  inv_variable_table[var_id] = key;

  // Setting maximum and minimum lags
  if (max_lead < lag)
    max_lead = lag;
  else if (-max_lag > lag)
    max_lag = -lag;

  switch(type)
    {
    case eEndogenous:
      var_endo_nbr++;
      if (max_endo_lead < lag)
        max_endo_lead = lag;
      else if (-max_endo_lag > lag)
        max_endo_lag = -lag;
      break;
    case eExogenous:
      var_exo_nbr++;
      if (max_exo_lead < lag)
        max_exo_lead = lag;
      else if (-max_exo_lag > lag)
        max_exo_lag = -lag;
      break;
    case eExogenousDet:
      var_exo_det_nbr++;
      if (max_exo_det_lead < lag)
        max_exo_det_lead = lag;
      else if (-max_exo_det_lag > lag)
        max_exo_det_lag = -lag;
      break;
    default:
      cerr << "VariableTable::addVariable(): forbidden variable type" << endl;
      exit(EXIT_FAILURE);
    }
  return var_id;
}

void
VariableTable::computeDynJacobianCols() throw (DynJacobianColsAlreadyComputedException)
{
  if (dyn_jacobian_cols_table.size() != 0)
    throw DynJacobianColsAlreadyComputedException();

  dyn_jacobian_cols_table.resize(size());

  variable_table_type::const_iterator it = variable_table.begin();

  // Assign the first columns to endogenous, using the lexicographic order over (lag, symbol_id) implemented in variable_table map
  int sorted_id = 0;
  while(it->first.first.first == eEndogenous && it != variable_table.end())
    {
      dyn_jacobian_cols_table[it->second] = sorted_id++;
      it++;
    }

  // Assign subsequent columns to exogenous and then exogenous deterministic, using an offset + symbol_id
  while(it->first.first.first == eExogenous && it != variable_table.end())
    {
      dyn_jacobian_cols_table[it->second] = var_endo_nbr + it->first.second;
      it++;
    }
  while(it->first.first.first == eExogenousDet && it != variable_table.end())
    {
      dyn_jacobian_cols_table[it->second] = var_endo_nbr + symbol_table.exo_nbr + it->first.second;
      it++;
    }
}
