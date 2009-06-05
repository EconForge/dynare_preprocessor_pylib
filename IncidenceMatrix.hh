/*
 * Copyright (C) 2007-2008 Dynare Team
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

#ifndef _INCIDENCEMATRIX_HH
#define _INCIDENCEMATRIX_HH


#include <map>
#include "ExprNode.hh"
#include "SymbolTable.hh"





//! List of incidence matrix (one matrix per lead/lag)
typedef bool* pbool;
typedef map<int,pbool> IncidenceList;

//! create and manage the incidence matrix
class IncidenceMatrix
{
public:
  const SymbolTable &symbol_table;
  IncidenceMatrix(const SymbolTable &symbol_table_arg);
  bool* Build_IM(int lead_lag, SymbolType type);
  bool* Get_IM(int lead_lag, SymbolType type) const;
  void fill_IM(int equation, int variable_endo, int lead_lag, SymbolType type);
  void unfill_IM(int equation, int variable_endo, int lead_lag, SymbolType type);
  void Free_IM() const;
  void Print_IM(SymbolType type) const;
  void Print_SIM(bool* IM, SymbolType type) const;
  void swap_IM_c(bool *SIM, int pos1, int pos2, int pos3, vector<int> &Index_Var_IM, vector<int> &Index_Equ_IM, int n) const;
  int Model_Max_Lead, Model_Max_Lag;
  int Model_Max_Lead_Endo, Model_Max_Lag_Endo, Model_Max_Lead_Exo, Model_Max_Lag_Exo;
private:
  IncidenceList  List_IM, List_IM_X;
};


#endif
