/*
 * Copyright (C) 2007-2009 Dynare Team
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

#include <iostream>
#include <cstdlib>
#include <cstring>

#include "IncidenceMatrix.hh"


IncidenceMatrix::IncidenceMatrix(const SymbolTable &symbol_table_arg) :
  symbol_table(symbol_table_arg)
{
  Model_Max_Lead = Model_Max_Lead_Endo = Model_Max_Lead_Exo = 0;
  Model_Max_Lag = Model_Max_Lag_Endo = Model_Max_Lag_Exo = 0;
}
//------------------------------------------------------------------------------
//For a lead or a lag build the Incidence Matrix structures
bool*
IncidenceMatrix::Build_IM(int lead_lag, SymbolType type)
{
  int size;
  bool *IM;
  if(type==eEndogenous)
    {
      size = symbol_table.endo_nbr() * symbol_table.endo_nbr() * sizeof(IM[0]);
      List_IM[lead_lag] = IM = (bool*)malloc(size);
      for(int i = 0; i< symbol_table.endo_nbr() * symbol_table.endo_nbr(); i++) IM[i] = 0;
      if(lead_lag > 0)
        {
          if(lead_lag > Model_Max_Lead_Endo)
            {
              Model_Max_Lead_Endo = lead_lag;
              if(lead_lag > Model_Max_Lead)
                Model_Max_Lead = lead_lag;
            }
        }
      else
        {
          if( -lead_lag > Model_Max_Lag_Endo)
            {
              Model_Max_Lag_Endo = -lead_lag;
              if(-lead_lag > Model_Max_Lag)
                Model_Max_Lag = -lead_lag;
            }
        }
    }
  else
    {  //eExogenous
      size = symbol_table.endo_nbr() * symbol_table.exo_nbr() * sizeof(IM[0]);
      List_IM_X[lead_lag] = IM = (bool*)malloc(size);
      for(int i = 0; i< symbol_table.endo_nbr() * symbol_table.exo_nbr(); i++) IM[i] = 0;
      if(lead_lag > 0)
        {
          if(lead_lag > Model_Max_Lead_Exo)
            {
              Model_Max_Lead_Exo = lead_lag;
              if(lead_lag > Model_Max_Lead)
                Model_Max_Lead = lead_lag;
            }
        }
      else
        {
          if( -lead_lag > Model_Max_Lag_Exo)
            {
              Model_Max_Lag_Exo = -lead_lag;
              if(-lead_lag > Model_Max_Lag)
                Model_Max_Lag = -lead_lag;
            }
        }
    }
  return (IM);
}


void
IncidenceMatrix::Free_IM() const
{
  IncidenceList::const_iterator it = List_IM.begin();
  for(it = List_IM.begin(); it != List_IM.end(); it++)
    free(it->second);
  for(it = List_IM_X.begin(); it != List_IM_X.end(); it++)
    free(it->second);
}

//------------------------------------------------------------------------------
// Return the incidence matrix related to a lead or a lag
bool*
IncidenceMatrix::Get_IM(int lead_lag, SymbolType type) const
{
  IncidenceList::const_iterator it;
  if(type==eEndogenous)
    {
      it = List_IM.find(lead_lag);
      if(it!=List_IM.end())
        return(it->second);
      else
        return(NULL);
    }
  else  //eExogenous
    {
      it = List_IM_X.find(lead_lag);
      if(it!=List_IM_X.end())
        return(it->second);
      else
        return(NULL);
    }
}


//------------------------------------------------------------------------------
// Fill the incidence matrix related to a lead or a lag
void
IncidenceMatrix::fill_IM(int equation, int variable, int lead_lag, SymbolType type)
{
  bool* Cur_IM;
  Cur_IM = Get_IM(lead_lag, type);
  if(equation >= symbol_table.endo_nbr())
    {
      cout << "Error : The model has more equations (at least " << equation + 1 << ") than declared endogenous variables (" << symbol_table.endo_nbr() << ")\n";
      exit(EXIT_FAILURE);
    }
  if (!Cur_IM)
    Cur_IM = Build_IM(lead_lag, type);
  if(type==eEndogenous)
    Cur_IM[equation*symbol_table.endo_nbr() + variable] = 1;
  else
    Cur_IM[equation*symbol_table.exo_nbr() + variable] = 1;
}

//------------------------------------------------------------------------------
// unFill the incidence matrix related to a lead or a lag
void
IncidenceMatrix::unfill_IM(int equation, int variable, int lead_lag, SymbolType type)
{
  bool* Cur_IM;
  Cur_IM = Get_IM(lead_lag, type);
  if (!Cur_IM)
    Cur_IM = Build_IM(lead_lag, type);
  if(type==eEndogenous)
    Cur_IM[equation*symbol_table.endo_nbr() + variable] = 0;
  else
    Cur_IM[equation*symbol_table.exo_nbr() + variable] = 0;
}


//------------------------------------------------------------------------------
//Print azn incidence matrix
void
IncidenceMatrix::Print_SIM(bool* IM, SymbolType type) const
{
  int i, j, n;
  if(type == eEndogenous)
    n = symbol_table.endo_nbr();
  else
    n = symbol_table.exo_nbr();
  for(i = 0;i < symbol_table.endo_nbr();i++)
    {
      cout << " ";
      for(j = 0;j < n;j++)
        cout << IM[i*n + j] << " ";
      cout << "\n";
    }
}

//------------------------------------------------------------------------------
//Print all incidence matrix
void
IncidenceMatrix::Print_IM(SymbolType type) const
{
  IncidenceList::const_iterator it;
  cout << "-------------------------------------------------------------------\n";
  if(type == eEndogenous)
    for(int k=-Model_Max_Lag_Endo; k <= Model_Max_Lead_Endo; k++)
      {
        it = List_IM.find(k);
        if(it!=List_IM.end())
          {
            cout << "Incidence matrix for lead_lag = " << k << "\n";
            Print_SIM(it->second, type);
          }
      }
  else // eExogenous
    for(int k=-Model_Max_Lag_Exo; k <= Model_Max_Lead_Exo; k++)
      {
        it = List_IM_X.find(k);
        if(it!=List_IM_X.end())
          {
            cout << "Incidence matrix for lead_lag = " << k << "\n";
            Print_SIM(it->second, type);
          }
      }
}


//------------------------------------------------------------------------------
// Swap rows and columns of the incidence matrix
void
IncidenceMatrix::swap_IM_c(bool *SIM, int pos1, int pos2, int pos3, vector<int> &Index_Var_IM, vector<int> &Index_Equ_IM, int n) const
{
  int tmp_i, j;
  bool tmp_b;
  /* We exchange equation (row)...*/
  if(pos1 != pos2)
    {
      tmp_i = Index_Equ_IM[pos1];
      Index_Equ_IM[pos1] = Index_Equ_IM[pos2];
      Index_Equ_IM[pos2] = tmp_i;
      for(j = 0;j < n;j++)
        {
          tmp_b = SIM[pos1 * n + j];
          SIM[pos1*n + j] = SIM[pos2 * n + j];
          SIM[pos2*n + j] = tmp_b;
        }
    }
  /* ...and variables (column)*/
  if(pos1 != pos3)
    {
      tmp_i = Index_Var_IM[pos1];
      Index_Var_IM[pos1] = Index_Var_IM[pos3];
      Index_Var_IM[pos3] = tmp_i;
      for(j = 0;j < n;j++)
        {
          tmp_b = SIM[j * n + pos1];
          SIM[j*n + pos1] = SIM[j * n + pos3];
          SIM[j*n + pos3] = tmp_b;
        }
    }
}
