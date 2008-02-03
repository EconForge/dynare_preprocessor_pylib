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

#include "VariableTable.hh"

VariableTable::VariableTable(const SymbolTable &symbol_table_arg) :
  symbol_table(symbol_table_arg),
  var_endo_nbr(0), var_exo_nbr(0), var_exo_det_nbr(0),
  max_lag(0), max_lead(0),
  max_endo_lag(0), max_endo_lead(0),
  max_exo_lag(0), max_exo_lead(0),
  max_exo_det_lag(0), max_exo_det_lead(0),
  max_recur_lag(0), max_recur_lead(0)
{
}

int
VariableTable::addVariable(Type type, int symb_id, int lag) throw (AlreadySortedException)
{
  if (sorted_ids_table.size() != 0)
    throw AlreadySortedException();

  var_key_type key = make_pair(make_pair(type, lag), symb_id);

  // Testing if variable already exists
  variable_table_type::const_iterator it = variable_table.find(key);
  if (it != variable_table.end())
    return it->second;

  int var_id = size();

  variable_table[key] = var_id;
  inv_variable_table[var_id] = key;

  if (type == eEndogenous)
    var_endo_nbr++;
  if (type == eExogenous)
    var_exo_nbr++;
  if (type == eExogenousDet)
    var_exo_det_nbr++;

  // Setting maximum and minimum lags
  if (max_lead < lag)
    max_lead = lag;
  else if (-max_lag > lag)
    max_lag = -lag;

  switch(type)
    {
    case eEndogenous:
      if (max_endo_lead < lag)
        max_endo_lead = lag;
      else if (-max_endo_lag > lag)
        max_endo_lag = -lag;
      break;
    case eExogenous:
      if (max_exo_lead < lag)
        max_exo_lead = lag;
      else if (-max_exo_lag > lag)
        max_exo_lag = -lag;
      break;
    case eExogenousDet:
      if (max_exo_det_lead < lag)
        max_exo_det_lead = lag;
      else if (-max_exo_det_lag > lag)
        max_exo_det_lag = -lag;
      break;
    case eRecursiveVariable:
      if (max_recur_lead < lag)
        max_recur_lead = lag;
      else if (-max_recur_lag > lag)
        max_recur_lag = -lag;
      break;
    default:
      ;
    }
  return var_id;
}

void
VariableTable::sort() throw (AlreadySortedException)
{
  if (sorted_ids_table.size() != 0)
    throw AlreadySortedException();

  int sorted_id = 0;
  sorted_ids_table.resize(size());

  for(variable_table_type::const_iterator it = variable_table.begin();
      it != variable_table.end(); it++)
    sorted_ids_table[it->second] = sorted_id++;
}
