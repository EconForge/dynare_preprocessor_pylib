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

// TODO Apply Block Decomposition to the static model

#include <iostream>
#include <sstream>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <cstring>
#include <cmath>
using namespace std;
//------------------------------------------------------------------------------
#include "BlockTriangular.hh"
//------------------------------------------------------------------------------

BlockTriangular::BlockTriangular(const SymbolTable &symbol_table_arg) :
  symbol_table(symbol_table_arg),
  normalization(symbol_table_arg)
{
  bt_verbose = 0;
  ModelBlock = NULL;
  Model_Max_Lead = 0;
  Model_Max_Lag = 0;
  periods = 0;
}



//------------------------------------------------------------------------------
BlockTriangular::~BlockTriangular()
{
  // Empty
}

//------------------------------------------------------------------------------
//For a lead or a lag build the Incidence Matrix structures
List_IM*
BlockTriangular::Build_IM(int lead_lag)
{
  List_IM* pIM = new List_IM;
  int i;
  Last_IM->pNext = pIM;
  pIM->IM = (bool*)malloc(endo_nbr * endo_nbr * sizeof(pIM->IM[0]));
  for(i = 0;i < endo_nbr*endo_nbr;i++)
    pIM->IM[i] = 0;
  pIM->lead_lag = lead_lag;
  if(lead_lag > 0)
    {
      if(lead_lag > Model_Max_Lead)
        Model_Max_Lead = lead_lag;
    }
  else
    {
      if( -lead_lag > Model_Max_Lag)
        Model_Max_Lag = -lead_lag;
    }
  pIM->pNext = NULL;
  Last_IM = pIM;
  return (pIM);
}

//------------------------------------------------------------------------------
// initialize all the incidence matrix structures
void
BlockTriangular::init_incidence_matrix(int nb_endo)
{
  int i;
  endo_nbr = nb_endo;
  First_IM = new List_IM;
  First_IM->IM = (bool*)malloc(nb_endo * nb_endo * sizeof(First_IM->IM[0]));
  for(i = 0;i < nb_endo*nb_endo;i++)
    First_IM->IM[i] = 0;
  First_IM->lead_lag = 0;
  First_IM->pNext = NULL;
  Last_IM = First_IM;
  //cout << "init_incidence_matrix done \n";
}


void
BlockTriangular::Free_IM(List_IM* First_IM)
{
  List_IM *Cur_IM, *SFirst_IM;
  Cur_IM = SFirst_IM = First_IM;
  while(Cur_IM)
    {
      First_IM = Cur_IM->pNext;
      free(Cur_IM->IM);
      Cur_IM = First_IM;
    }
  free(SFirst_IM);
}

//------------------------------------------------------------------------------
// Return the inceidence matrix related to a lead or a lag
List_IM*
BlockTriangular::Get_IM(int lead_lag)
{
  List_IM* Cur_IM;
  Cur_IM = First_IM;
  while ((Cur_IM != NULL) && (Cur_IM->lead_lag != lead_lag))
    Cur_IM = Cur_IM->pNext;
  return (Cur_IM);
}

bool*
BlockTriangular::bGet_IM(int lead_lag)
{
  List_IM* Cur_IM;
  Cur_IM = First_IM;
  while ((Cur_IM != NULL) && (Cur_IM->lead_lag != lead_lag))
    {
      Cur_IM = Cur_IM->pNext;
    }
  if((Cur_IM->lead_lag != lead_lag) || (Cur_IM==NULL))
    {
      cout << "the incidence matrix with lag " << lead_lag << " does not exist !!";
      exit(-1);
    }
  return (Cur_IM->IM);
}

bool*
BlockTriangular::bGet_IM(int lead_lag) const
{
  List_IM* Cur_IM;
  Cur_IM = First_IM;
  while ((Cur_IM != NULL) && (Cur_IM->lead_lag != lead_lag))
    {
      Cur_IM = Cur_IM->pNext;
    }
  if((Cur_IM->lead_lag != lead_lag) || (Cur_IM==NULL))
    {
      cout << "the incidence matrix with lag " << lead_lag << " does not exist !!";
      exit(-1);
    }
  return (Cur_IM->IM);
}


//------------------------------------------------------------------------------
// Fill the incidence matrix related to a lead or a lag
void
BlockTriangular::fill_IM(int equation, int variable_endo, int lead_lag)
{
  List_IM* Cur_IM;
  //cout << "equation=" << equation << " variable_endo=" << variable_endo << " lead_lag=" << lead_lag << "\n";
  Cur_IM = Get_IM(lead_lag);
  if(equation >= endo_nbr)
    {
      cout << "Error : The model has more equations (at least " << equation + 1 << ") than declared endogenous variables (" << endo_nbr << ")\n";
      system("PAUSE");
      exit( -1);
    }
  if (!Cur_IM)
    Cur_IM = Build_IM(lead_lag);
  Cur_IM->IM[equation*endo_nbr + variable_endo] = 1;
}

//------------------------------------------------------------------------------
// unFill the incidence matrix related to a lead or a lag
void
BlockTriangular::unfill_IM(int equation, int variable_endo, int lead_lag)
{
  List_IM* Cur_IM;
  //cout << "lead_lag=" << lead_lag << "\n";
  Cur_IM = Get_IM(lead_lag);
  /*if(equation >= endo_nbr)
    {
    cout << "Error : The model has more equations (at least " << equation + 1 << ") than declared endogenous variables (" << endo_nbr << ")\n";
    system("PAUSE");
    exit( -1);
    }*/
  if (!Cur_IM)
    Cur_IM = Build_IM(lead_lag);
  Cur_IM->IM[equation*endo_nbr + variable_endo] = 0;
  /*system("pause");*/
}


//------------------------------------------------------------------------------
//Print azn incidence matrix
void
BlockTriangular::Print_SIM(bool* IM, int n) const
{
  int i, j;
  for(i = 0;i < n;i++)
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
BlockTriangular::Print_IM(int n) const
{
  List_IM* Cur_IM;
  Cur_IM = First_IM;
  cout << "-------------------------------------------------------------------\n";
  while(Cur_IM)
    {
      cout << "Incidence matrix for lead_lag = " << Cur_IM->lead_lag << "\n";
      Print_SIM(Cur_IM->IM, n);
      Cur_IM = Cur_IM->pNext;
    }
}


//------------------------------------------------------------------------------
// Swap rows and columns of the incidence matrix
void
BlockTriangular::swap_IM_c(bool *SIM, int pos1, int pos2, int pos3, simple* Index_Var_IM, simple* Index_Equ_IM, int n)
{
  int tmp_i, j;
  bool tmp_b;
  /* We exchange equation (row)...*/
  if(pos1 != pos2)
    {
      tmp_i = Index_Equ_IM[pos1].index;
      Index_Equ_IM[pos1].index = Index_Equ_IM[pos2].index;
      Index_Equ_IM[pos2].index = tmp_i;
      for(j = 0;j < n;j++)
        {
          tmp_b = SIM[pos1 * n + j];
          SIM[pos1*n + j] = SIM[pos2 * n + j];
          SIM[pos2*n + j] = tmp_b;
        }
    }
  /* ...and variables (colomn)*/
  if(pos1 != pos3)
    {
      tmp_i = Index_Var_IM[pos1].index;
      Index_Var_IM[pos1].index = Index_Var_IM[pos3].index;
      Index_Var_IM[pos3].index = tmp_i;
      for(j = 0;j < n;j++)
        {
          tmp_b = SIM[j * n + pos1];
          SIM[j*n + pos1] = SIM[j * n + pos3];
          SIM[j*n + pos3] = tmp_b;
        }
    }
}

//------------------------------------------------------------------------------
// Find the prologue and the epilogue of the model
void
BlockTriangular::Prologue_Epilogue(bool* IM, int* prologue, int* epilogue, int n, simple* Index_Var_IM, simple* Index_Equ_IM)
{
  bool modifie = 1;
  int i, j, k, l = 0;
  /*Looking for a prologue */
  *prologue = 0;
  while(modifie)
    {
      modifie = 0;
      for(i = *prologue;i < n;i++)
        {
          k = 0;
          for(j = *prologue;j < n;j++)
            {
              if(IM[i*n + j])
                {
                  k++;
                  l = j;
                }
            }
          if ((k == 1) /* && (l==i)*/)
            {
              modifie = 1;
              swap_IM_c(IM, *prologue, i, l, Index_Var_IM, Index_Equ_IM, n);
              Index_Equ_IM[*prologue].available = 0;
              Index_Var_IM[*prologue].available = 0;
              (*prologue)++;
            }
        }
    }
  *epilogue = 0;
  modifie = 1;
  while(modifie)
    {
      modifie = 0;
      for(i = *prologue;i < n - *epilogue;i++)
        {
          k = 0;
          for(j = *prologue;j < n - *epilogue;j++)
            {
              if(IM[j*n + i])
                {
                  k++;
                  l = j;
                }
            }
          if ((k == 1) /* && (l==i)*/)
            {
              modifie = 1;
              swap_IM_c(IM, n - (1 + *epilogue), l, i, Index_Var_IM, Index_Equ_IM, n);
              Index_Equ_IM[n - (1 + *epilogue)].available = 0;
              Index_Var_IM[n - (1 + *epilogue)].available = 0;
              (*epilogue)++;
            }
        }
    }
}


void
BlockTriangular::getMax_Lead_Lag(int var, int equ, int *lead, int *lag)
{
  List_IM* Cur_IM;
  Cur_IM = First_IM->pNext;
  (*lead) = (*lag) = 0;
  while(Cur_IM)
    {
      if(Cur_IM->IM[equ*endo_nbr + var])
        {
          if(Cur_IM->lead_lag > 0)
            {
              if ((*lead) < Cur_IM->lead_lag)
                *lead = Cur_IM->lead_lag;
            }
          else
            {
              if ((*lag) < abs(Cur_IM->lead_lag))
                *lag = abs(Cur_IM->lead_lag);
            }
        }
      Cur_IM = Cur_IM->pNext;
    }
}

void
BlockTriangular::getMax_Lead_Lag_B(int size, int* Equation, int *Variable, int *lead, int *lag)
{
  List_IM* Cur_IM;
  int i, j;
  Cur_IM = First_IM->pNext;
  (*lead) = (*lag) = 0;
  while(Cur_IM)
    {
      for(i = 0;i < size;i++)
        {
          for(j = 0;j < size;j++)
            {
              if(Cur_IM->IM[Equation[i]*endo_nbr + Variable[j]])
                {
                  if(Cur_IM->lead_lag > 0)
                    {
                      if ((*lead) < Cur_IM->lead_lag)
                        *lead = Cur_IM->lead_lag;
                    }
                  else
                    {
                      if ((*lag) < abs(Cur_IM->lead_lag))
                        *lag = abs(Cur_IM->lead_lag);
                    }
                }
            }
        }
      Cur_IM = Cur_IM->pNext;
    }
}


void
BlockTriangular::Allocate_Block(int size, int *count_Equ, int *count_Block, int type, Model_Block * ModelBlock)
{
  int i, j, k, l, ls, m, i_1, Lead, Lag, size_list_lead_var, first_count_equ, i1;
  int *list_lead_var, *tmp_size, *tmp_var, *tmp_endo, nb_lead_lag_endo;
  List_IM *Cur_IM;
  bool *IM, OK;
  ModelBlock->Periods = periods;
  if ((type == PROLOGUE) || (type == EPILOGUE))
    {
      for(i = 0;i < size;i++)
        {
          ModelBlock->Block_List[*count_Block].is_linear=true;
          ModelBlock->Block_List[*count_Block].Size = 1;
          ModelBlock->Block_List[*count_Block].Type = type;
          ModelBlock->Block_List[*count_Block].Simulation_Type = UNKNOWN;
          ModelBlock->Block_List[*count_Block].Temporary_terms=new temporary_terms_type ();
          ModelBlock->Block_List[*count_Block].Temporary_terms->clear();
          list_lead_var = (int*)malloc(Model_Max_Lead * endo_nbr * sizeof(int));
          size_list_lead_var = 0;
          tmp_endo = (int*)malloc((Model_Max_Lead + Model_Max_Lag + 1) * sizeof(int));
          tmp_size = (int*)malloc((Model_Max_Lead + Model_Max_Lag + 1) * sizeof(int));
          tmp_var = (int*)malloc(sizeof(int));
          memset(tmp_size, 0, (Model_Max_Lead + Model_Max_Lag + 1)*sizeof(int));
          memset(tmp_endo, 0, (Model_Max_Lead + Model_Max_Lag + 1)*sizeof(int));
          nb_lead_lag_endo = Lead = Lag = 0;
          Cur_IM = First_IM;
          while(Cur_IM)
            {
              k = Cur_IM->lead_lag;
              i_1 = Index_Var_IM[*count_Equ].index * endo_nbr;
              if(k > 0)
                {
                  if(Cur_IM->IM[i_1 + Index_Equ_IM[ /*j*/*count_Equ].index])
                    {
                      nb_lead_lag_endo++;
                      tmp_size[Model_Max_Lag + k]++;
                      if(k > Lead)
                        {
                          Lead = k;
                          list_lead_var[size_list_lead_var] = Index_Var_IM[*count_Equ].index + size * (k - 1);
                          size_list_lead_var++;
                        }
                    }
                }
              else
                {
                  k = -k;
                  if(Cur_IM->IM[i_1 + Index_Equ_IM[ /*j*/*count_Equ].index])
                    {
                      tmp_size[Model_Max_Lag - k]++;
                      nb_lead_lag_endo++;
                      if(k > Lag)
                        {
                          Lag = k;
                        }
                    }
                }
              Cur_IM = Cur_IM->pNext;
            }
          ModelBlock->Block_List[*count_Block].Max_Lag = Lag;
          ModelBlock->Block_List[*count_Block].Max_Lead = Lead;
          free(list_lead_var);
          ModelBlock->Block_List[*count_Block].Equation = (int*)malloc(sizeof(int));
          ModelBlock->Block_List[*count_Block].Variable = (int*)malloc(sizeof(int));
          ModelBlock->Block_List[*count_Block].Variable_Sorted = (int*)malloc(sizeof(int));
          ModelBlock->Block_List[*count_Block].Own_Derivative = (int*)malloc(sizeof(int));
          ModelBlock->Block_List[*count_Block].Equation[0] = Index_Equ_IM[*count_Equ].index;
          ModelBlock->Block_List[*count_Block].Variable[0] = Index_Var_IM[*count_Equ].index;
          ModelBlock->Block_List[*count_Block].Variable_Sorted[0] = -1;
          ModelBlock->in_Block_Equ[Index_Equ_IM[*count_Equ].index] = *count_Block;
          ModelBlock->in_Block_Var[Index_Var_IM[*count_Equ].index] = *count_Block;
          ModelBlock->in_Equ_of_Block[Index_Equ_IM[*count_Equ].index] = ModelBlock->in_Var_of_Block[Index_Var_IM[*count_Equ].index] = 0;
          Index_Equ_IM[*count_Equ].block = *count_Block;
          cout << "Lead=" << Lead << " Lag=" << Lag << "\n";
          if ((Lead > 0) && (Lag > 0))
            {
              ModelBlock->Block_List[*count_Block].Simulation_Type = SOLVE_TWO_BOUNDARIES_SIMPLE;
              cout << "alloc ModelBlock->Block_List[" << *count_Block << "].IM_lead_lag = (" << (Lead + Lag + 1) * sizeof(IM_compact) << ")\n";
              ModelBlock->Block_List[*count_Block].IM_lead_lag = (IM_compact*)malloc((Lead + Lag + 1) * sizeof(IM_compact));
              ModelBlock->Block_List[*count_Block].Nb_Lead_Lag_Endo = nb_lead_lag_endo;
              ModelBlock->Block_List[*count_Block].variable_dyn_index = (int*)malloc(nb_lead_lag_endo * sizeof(int));
              ModelBlock->Block_List[*count_Block].variable_dyn_leadlag = (int*)malloc(nb_lead_lag_endo * sizeof(int));
              ls = l = 1;
              i1 = 0;
              for(i = 0;i < Lead + Lag + 1;i++)
                {
                  ModelBlock->Block_List[*count_Block].IM_lead_lag[i].size = tmp_size[Model_Max_Lag - Lag + i];
                  ModelBlock->Block_List[*count_Block].IM_lead_lag[i].nb_endo = tmp_size[Model_Max_Lag - Lag + i];
                  ModelBlock->Block_List[*count_Block].IM_lead_lag[i].u = (int*)malloc(tmp_size[Model_Max_Lag - Lag + i] * sizeof(int));
                  ModelBlock->Block_List[*count_Block].IM_lead_lag[i].us = (int*)malloc(tmp_size[Model_Max_Lag - Lag + i] * sizeof(int));
                  ModelBlock->Block_List[*count_Block].IM_lead_lag[i].Var = (int*)malloc(tmp_size[Model_Max_Lag - Lag + i] * sizeof(int));
                  ModelBlock->Block_List[*count_Block].IM_lead_lag[i].Equ = (int*)malloc(tmp_size[Model_Max_Lag - Lag + i] * sizeof(int));
                  ModelBlock->Block_List[*count_Block].IM_lead_lag[i].Var_Index = (int*)malloc(tmp_size[Model_Max_Lag - Lag + i] * sizeof(int));
                  ModelBlock->Block_List[*count_Block].IM_lead_lag[i].Equ_Index = (int*)malloc(tmp_size[Model_Max_Lag - Lag + i] * sizeof(int));
                  ModelBlock->Block_List[*count_Block].IM_lead_lag[i].Var_dyn_Index = (int*)malloc(tmp_size[Model_Max_Lag - Lag + i] * sizeof(int));
                  ModelBlock->Block_List[*count_Block].IM_lead_lag[i].u_init = l;
                  IM = bGet_IM(i - Lag);
                  if(IM == NULL)
                    {
                      cout << "Error IM(" << i - Lag << ") doesn't exist\n";
                      exit( -1);
                    }
                  if(IM[Index_Var_IM[*count_Equ].index + Index_Equ_IM[*count_Equ].index*endo_nbr])
                    {
                      ModelBlock->Block_List[*count_Block].variable_dyn_index[i1] = Index_Var_IM[*count_Equ].index;
                      ModelBlock->Block_List[*count_Block].variable_dyn_leadlag[i1] = i - Lag;
                      tmp_var[0] = i1;
                      i1++;
                    }
                  m = 0;
                  i_1 = Index_Equ_IM[*count_Equ].index * endo_nbr;
                  if(IM[Index_Var_IM[*count_Equ].index + i_1])
                    {
                      if(i == Lag)
                        {
                          ModelBlock->Block_List[*count_Block].IM_lead_lag[i].us[m] = ls;
                          ls++;
                        }
                      ModelBlock->Block_List[*count_Block].IM_lead_lag[i].u[m] = l;
                      ModelBlock->Block_List[*count_Block].IM_lead_lag[i].Equ[m] = 0;
                      ModelBlock->Block_List[*count_Block].IM_lead_lag[i].Var[m] = 0;
                      ModelBlock->Block_List[*count_Block].IM_lead_lag[i].Equ_Index[m] = Index_Equ_IM[*count_Equ].index;
                      ModelBlock->Block_List[*count_Block].IM_lead_lag[i].Var_Index[m] = Index_Var_IM[*count_Equ].index;
                      ModelBlock->Block_List[*count_Block].IM_lead_lag[i].Var_dyn_Index[m] = ModelBlock->Block_List[*count_Block].variable_dyn_index[tmp_var[0]];
                      l++;
                      m++;
                    }
                  ModelBlock->Block_List[*count_Block].IM_lead_lag[i].u_finish = l - 1;
                }
            }
          else if((Lead > 0) && (Lag == 0))
            ModelBlock->Block_List[*count_Block].Simulation_Type = SOLVE_BACKWARD_SIMPLE;
          else
            ModelBlock->Block_List[*count_Block].Simulation_Type = SOLVE_FOREWARD_SIMPLE;
          (*count_Equ)++;
          (*count_Block)++;
          free(tmp_size);
          free(tmp_endo);
        }
    }
  else
    {
      ModelBlock->Block_List[*count_Block].is_linear=true;
      ModelBlock->Block_List[*count_Block].Size = size;
      ModelBlock->Block_List[*count_Block].Type = type;
      ModelBlock->Block_List[*count_Block].Temporary_terms=new temporary_terms_type ();
      ModelBlock->Block_List[*count_Block].Temporary_terms->clear();
      ModelBlock->Block_List[*count_Block].Simulation_Type = UNKNOWN;
      ModelBlock->Block_List[*count_Block].Equation = (int*)malloc(ModelBlock->Block_List[*count_Block].Size * sizeof(int));
      ModelBlock->Block_List[*count_Block].Variable = (int*)malloc(ModelBlock->Block_List[*count_Block].Size * sizeof(int));
      ModelBlock->Block_List[*count_Block].Variable_Sorted = (int*)malloc(ModelBlock->Block_List[*count_Block].Size * sizeof(int));
      ModelBlock->Block_List[*count_Block].Own_Derivative = (int*)malloc(ModelBlock->Block_List[*count_Block].Size * sizeof(int));
      Lead = Lag = 0;
      first_count_equ = *count_Equ;
      tmp_var = (int*)malloc(size * sizeof(int));
      tmp_endo = (int*)malloc((Model_Max_Lead + Model_Max_Lag + 1) * sizeof(int));
      tmp_size = (int*)malloc((Model_Max_Lead + Model_Max_Lag + 1) * sizeof(int));
      memset(tmp_size, 0, (Model_Max_Lead + Model_Max_Lag + 1)*sizeof(int));
      memset(tmp_endo, 0, (Model_Max_Lead + Model_Max_Lag + 1)*sizeof(int));
      nb_lead_lag_endo = 0;
      for(i = 0;i < size;i++)
        {
          Index_Equ_IM[*count_Equ].block = *count_Block;
          ModelBlock->Block_List[*count_Block].Equation[i] = Index_Equ_IM[*count_Equ].index;
          ModelBlock->Block_List[*count_Block].Variable[i] = Index_Var_IM[*count_Equ].index;
          ModelBlock->Block_List[*count_Block].Variable_Sorted[i] = -1;
          ModelBlock->in_Block_Equ[Index_Equ_IM[*count_Equ].index] = *count_Block;
          ModelBlock->in_Block_Var[Index_Var_IM[*count_Equ].index] = *count_Block;
          ModelBlock->in_Equ_of_Block[Index_Equ_IM[*count_Equ].index] = ModelBlock->in_Var_of_Block[Index_Var_IM[*count_Equ].index] = i;
          Cur_IM = First_IM;
          i_1 = Index_Var_IM[*count_Equ].index;
          while(Cur_IM)
            {
              k = Cur_IM->lead_lag;
              OK = false;
              if(k >= 0)
                {
                  for(j = 0;j < size;j++)
                    {
                      if(Cur_IM->IM[i_1 + Index_Equ_IM[first_count_equ + j].index*endo_nbr])
                        {
                          tmp_size[Model_Max_Lag + k]++;
                          if (!OK)
                            {
                              tmp_endo[Model_Max_Lag + k]++;
                              nb_lead_lag_endo++;
                              OK = true;
                            }
                          if(k > Lead)
                            Lead = k;
                        }
                    }
                }
              else
                {
                  k = -k;
                  for(j = 0;j < size;j++)
                    {
                      if(Cur_IM->IM[i_1 + Index_Equ_IM[first_count_equ + j].index*endo_nbr])
                        {
                          tmp_size[Model_Max_Lag - k]++;
                          if (!OK)
                            {
                              tmp_endo[Model_Max_Lag - k]++;
                              nb_lead_lag_endo++;
                              OK = true;
                            }
                          if(k > Lag)
                            Lag = k;
                        }
                    }
                }
              Cur_IM = Cur_IM->pNext;
            }
          (*count_Equ)++;
        }
      if ((Lag > 0) && (Lead > 0))
        ModelBlock->Block_List[*count_Block].Simulation_Type = SOLVE_TWO_BOUNDARIES_COMPLETE;
      else if(size > 1)
        {
          if(Lead > 0)
            ModelBlock->Block_List[*count_Block].Simulation_Type = SOLVE_BACKWARD_COMPLETE;
          else
            ModelBlock->Block_List[*count_Block].Simulation_Type = SOLVE_FOREWARD_COMPLETE;
        }
      else
        {
          if(Lead > 0)
            ModelBlock->Block_List[*count_Block].Simulation_Type = SOLVE_BACKWARD_SIMPLE;
          else
            ModelBlock->Block_List[*count_Block].Simulation_Type = SOLVE_FOREWARD_SIMPLE;
        }
      ModelBlock->Block_List[*count_Block].Max_Lag = Lag;
      ModelBlock->Block_List[*count_Block].Max_Lead = Lead;
      cout << "alloc ModelBlock->Block_List[" << *count_Block << "].IM_lead_lag = (" << (Lead + Lag + 1) * sizeof(IM_compact) << ")\n";
      ModelBlock->Block_List[*count_Block].IM_lead_lag = (IM_compact*)malloc((Lead + Lag + 1) * sizeof(IM_compact));
      ls = l = size;
      i1 = 0;
      ModelBlock->Block_List[*count_Block].Nb_Lead_Lag_Endo = nb_lead_lag_endo;
      ModelBlock->Block_List[*count_Block].variable_dyn_index = (int*)malloc(nb_lead_lag_endo * sizeof(int));
      ModelBlock->Block_List[*count_Block].variable_dyn_leadlag = (int*)malloc(nb_lead_lag_endo * sizeof(int));
      for(i = 0;i < Lead + Lag + 1;i++)
        {
          ModelBlock->Block_List[*count_Block].IM_lead_lag[i].size = tmp_size[Model_Max_Lag - Lag + i];
          ModelBlock->Block_List[*count_Block].IM_lead_lag[i].nb_endo = tmp_endo[Model_Max_Lag - Lag + i];
          ModelBlock->Block_List[*count_Block].IM_lead_lag[i].u = (int*)malloc(tmp_size[Model_Max_Lag - Lag + i] * sizeof(int));
          ModelBlock->Block_List[*count_Block].IM_lead_lag[i].us = (int*)malloc(tmp_size[Model_Max_Lag - Lag + i] * sizeof(int));
          ModelBlock->Block_List[*count_Block].IM_lead_lag[i].Var = (int*)malloc(tmp_size[Model_Max_Lag - Lag + i] * sizeof(int));
          ModelBlock->Block_List[*count_Block].IM_lead_lag[i].Equ = (int*)malloc(tmp_size[Model_Max_Lag - Lag + i] * sizeof(int));
          ModelBlock->Block_List[*count_Block].IM_lead_lag[i].Var_Index = (int*)malloc(tmp_size[Model_Max_Lag - Lag + i] * sizeof(int));
          ModelBlock->Block_List[*count_Block].IM_lead_lag[i].Equ_Index = (int*)malloc(tmp_size[Model_Max_Lag - Lag + i] * sizeof(int));
          ModelBlock->Block_List[*count_Block].IM_lead_lag[i].Var_dyn_Index = (int*)malloc(tmp_size[Model_Max_Lag - Lag + i] * sizeof(int));
          ModelBlock->Block_List[*count_Block].IM_lead_lag[i].u_init = l;
          IM = bGet_IM(i - Lag);
          if(IM != NULL)
            {
            }
          else
            {
              cout << "Error IM(" << i - Lag << ") doesn't exist\n";
              exit( -1);
            }
          for(j = first_count_equ;j < size + first_count_equ;j++)
            {
              i_1 = Index_Var_IM[j].index;
              m = 0;
              for(k = first_count_equ;k < size + first_count_equ;k++)
                if(IM[i_1 + Index_Equ_IM[k].index*endo_nbr])
                  m++;
              if(m > 0)
                {
                  ModelBlock->Block_List[*count_Block].variable_dyn_index[i1] = i_1;
                  ModelBlock->Block_List[*count_Block].variable_dyn_leadlag[i1] = i - Lag;
                  tmp_var[j - first_count_equ] = i1;
                  i1++;
                }
            }
          m = 0;
          for(j = first_count_equ;j < size + first_count_equ;j++)
            {
              i_1 = Index_Equ_IM[j].index * endo_nbr;
              for(k = first_count_equ;k < size + first_count_equ;k++)
                if(IM[Index_Var_IM[k].index + i_1])
                  {
                    if(i == Lag)
                      {
                        ModelBlock->Block_List[*count_Block].IM_lead_lag[i].us[m] = ls;
                        ls++;
#ifdef DEBUG
                        printf("j=%d ModelBlock->Block_List[%d].Variable[%d]=%d Index_Var_IM[%d].index=%d", j, *count_Block, j - first_count_equ, ModelBlock->Block_List[*count_Block].Variable[j - first_count_equ], k, Index_Var_IM[k].index);
                        if(ModelBlock->Block_List[*count_Block].Variable[j - first_count_equ] == Index_Var_IM[k].index)
                          {
                            ModelBlock->Block_List[*count_Block].Own_Derivative[j - first_count_equ]=l;
                            printf(" l=%d OK\n",l);
                          }
                        else
                          printf("\n");
#endif
                      }
                    ModelBlock->Block_List[*count_Block].IM_lead_lag[i].u[m] = l;
                    ModelBlock->Block_List[*count_Block].IM_lead_lag[i].Equ[m] = j - first_count_equ;
                    ModelBlock->Block_List[*count_Block].IM_lead_lag[i].Var[m] = k - first_count_equ;
                    ModelBlock->Block_List[*count_Block].IM_lead_lag[i].Equ_Index[m] = Index_Equ_IM[j].index;
                    ModelBlock->Block_List[*count_Block].IM_lead_lag[i].Var_Index[m] = Index_Var_IM[k].index;
                    ModelBlock->Block_List[*count_Block].IM_lead_lag[i].Var_dyn_Index[m] = ModelBlock->Block_List[*count_Block].variable_dyn_index[tmp_var[k - first_count_equ]];
                    l++;
                    m++;
                  }
            }
          ModelBlock->Block_List[*count_Block].IM_lead_lag[i].u_finish = l - 1;
        }
      (*count_Block)++;
      free(tmp_size);
      free(tmp_endo);
      free(tmp_var);
    }
}


void
BlockTriangular::Free_Block(Model_Block* ModelBlock)
{
  int blk, i;
#ifdef DEBUG

  cout << "Free_Block\n";
#endif

  for(blk = 0;blk < ModelBlock->Size;blk++)
    {

      if ((ModelBlock->Block_List[blk].Type == PROLOGUE) || (ModelBlock->Block_List[blk].Type == EPILOGUE))
        {
          free(ModelBlock->Block_List[blk].Equation);
          free(ModelBlock->Block_List[blk].Variable);
          free(ModelBlock->Block_List[blk].Variable_Sorted);
        }
      else
        {
          free(ModelBlock->Block_List[blk].Equation);
          free(ModelBlock->Block_List[blk].Variable);
          free(ModelBlock->Block_List[blk].Variable_Sorted);
          for(i = 0;i < ModelBlock->Block_List[blk].Max_Lag + ModelBlock->Block_List[blk].Max_Lead + 1;i++)
            {
              free(ModelBlock->Block_List[blk].IM_lead_lag[i].u);
              free(ModelBlock->Block_List[blk].IM_lead_lag[i].us);
              free(ModelBlock->Block_List[blk].IM_lead_lag[i].Equ);
              free(ModelBlock->Block_List[blk].IM_lead_lag[i].Var);
              free(ModelBlock->Block_List[blk].IM_lead_lag[i].Equ_Index);
              free(ModelBlock->Block_List[blk].IM_lead_lag[i].Var_Index);
            }
          free(ModelBlock->Block_List[blk].IM_lead_lag);
        }
    }
  free(ModelBlock->in_Block_Equ);
  free(ModelBlock->in_Block_Var);
  free(ModelBlock->in_Equ_of_Block);
  free(ModelBlock->in_Var_of_Block);
  free(ModelBlock->Block_List);
  free(ModelBlock);
}

string
BlockTriangular::getnamebyID(Type type, int id)
{
  return symbol_table.getNameByID(type,id);
}

//------------------------------------------------------------------------------
// Normalize each equation of the model (endgenous_i = f_i(endogenous_1, ..., endogenous_n) - in order to apply strong connex components search algorithm -
// and find the optimal blocks triangular decomposition
bool
BlockTriangular::Normalize_and_BlockDecompose(bool* IM, Model_Block* ModelBlock, int n, int* prologue, int* epilogue, simple* Index_Var_IM, simple* Index_Equ_IM, bool Do_Normalization, bool mixing, bool* IM_0, jacob_map j_m )
{
  int i, j, Nb_TotalBlocks, Nb_RecursBlocks;
  int count_Block, count_Equ;
  block_result_t* res;
  List_IM * p_First_IM, *p_Cur_IM, *Cur_IM;
  Equation_set* Equation_gr = (Equation_set*) malloc(sizeof(Equation_set));
  bool* SIM0, *SIM00;
  p_First_IM = (List_IM*)malloc(sizeof(*p_First_IM));
  p_Cur_IM = p_First_IM;
  Cur_IM = First_IM;
  i = endo_nbr * endo_nbr;
  while(Cur_IM)
    {
      p_Cur_IM->lead_lag = Cur_IM->lead_lag;
      p_Cur_IM->IM = (bool*)malloc(i * sizeof(bool));
      memcpy ( p_Cur_IM->IM, Cur_IM->IM, i);
      Cur_IM = Cur_IM->pNext;
      if(Cur_IM)
        {
          p_Cur_IM->pNext = (List_IM*)malloc(sizeof(*p_Cur_IM));
          p_Cur_IM = p_Cur_IM->pNext;
        }
      else
        p_Cur_IM->pNext = NULL;
    }
  Prologue_Epilogue(IM, prologue, epilogue, n, Index_Var_IM, Index_Equ_IM);
  if(bt_verbose)
    {
      cout << "prologue : " << *prologue << " epilogue : " << *epilogue << "\n";
      Print_SIM(IM, n);
      for(i = 0;i < n;i++)
        cout << "Index_Var_IM[" << i << "]=" << Index_Var_IM[i].index << " Index_Equ_IM[" << i << "]=" << Index_Equ_IM[i].index << "\n";
    }
  if(Do_Normalization)
    {
      cout << "Normalizing the model ...\n";
      if(mixing)
        {
          double* max_val=(double*)malloc(n*sizeof(double));
          memset(max_val,0,n*sizeof(double));
          for( map< pair< int, int >, double >::iterator iter = j_m.begin(); iter != j_m.end(); iter++ )
            {
              if(fabs(iter->second)>max_val[iter->first.first])
                max_val[iter->first.first]=fabs(iter->second);
            }
          for( map< pair< int, int >, double >::iterator iter = j_m.begin(); iter != j_m.end(); iter++ )
              iter->second/=max_val[iter->first.first];
          free(max_val);
          bool OK=false;
          double bi=0.99999999;
          int suppressed=0;
          while(!OK && bi>1e-14)
            {
              int suppress=0;
              SIM0 = (bool*)malloc(n * n * sizeof(bool));
              memset(SIM0,0,n*n*sizeof(bool));
              SIM00 = (bool*)malloc(n * n * sizeof(bool));
              memset(SIM00,0,n*n*sizeof(bool));
              for( map< pair< int, int >, double >::iterator iter = j_m.begin(); iter != j_m.end(); iter++ )
               {
                if(fabs(iter->second)>bi)
                  {
                    SIM0[iter->first.first*n+iter->first.second]=1;
                    if(!IM_0[iter->first.first*n+iter->first.second])
                      {
                        cout << "Error nothing at IM_0[" << iter->first.first << ", " << iter->first.second << "]=" << IM_0[iter->first.first*n+iter->first.second] << "\n";
                      }
                  }
                else
                  suppress++;
               }
              for(i = 0;i < n;i++)
                for(j = 0;j < n;j++)
                  SIM00[i*n + j] = SIM0[Index_Equ_IM[i].index * n + Index_Var_IM[j].index];
              free(SIM0);
              if(suppress!=suppressed)
                OK=normalization.Normalize(n, *prologue, *epilogue, SIM00, Index_Equ_IM, Equation_gr, 1, IM);
              suppressed=suppress;
              if(!OK)
                bi/=1.07;
              if(bi>1e-14)
                free(SIM00);
            }
          if(!OK)
            {
              normalization.Set_fp_verbose(true);
              OK=normalization.Normalize(n, *prologue, *epilogue, SIM00, Index_Equ_IM, Equation_gr, 1, IM);
              cout << "Error\n";
              exit(-1);
            }
        }
      else
        normalization.Normalize(n, *prologue, *epilogue, IM, Index_Equ_IM, Equation_gr, 0, 0);
    }
  else
    normalization.Gr_to_IM_basic(n, *prologue, *epilogue, IM, Equation_gr, false);
  cout << "Finding the optimal block decomposition of the model ...\n";
  if(bt_verbose)
    blocks.Print_Equation_gr(Equation_gr);
  res = blocks.sc(Equation_gr);
  if(bt_verbose)
    blocks.block_result_print(res);
  if ((*prologue) || (*epilogue))
    j = 1;
  else
    j = 0;
  for(i = 0;i < res->n_sets;i++)
    if ((res->sets_f[i] - res->sets_s[i] + 1) > j)
      j = res->sets_f[i] - res->sets_s[i] + 1;
  Nb_RecursBlocks = *prologue + *epilogue;
  Nb_TotalBlocks = res->n_sets + Nb_RecursBlocks;
  cout << Nb_TotalBlocks << " block(s) found:\n";
  cout << "  " << Nb_RecursBlocks << " recursive block(s) and " << res->n_sets << " simultaneous block(s). \n";
  cout << "  the largest simultaneous block has " << j << " equation(s). \n";
  ModelBlock->Size = Nb_TotalBlocks;
  ModelBlock->Periods = periods;
  ModelBlock->in_Block_Equ = (int*)malloc(n * sizeof(int));
  ModelBlock->in_Block_Var = (int*)malloc(n * sizeof(int));
  ModelBlock->in_Equ_of_Block = (int*)malloc(n * sizeof(int));
  ModelBlock->in_Var_of_Block = (int*)malloc(n * sizeof(int));
  ModelBlock->Block_List = (Block*)malloc(sizeof(ModelBlock->Block_List[0]) * Nb_TotalBlocks);
  blocks.block_result_to_IM(res, IM, *prologue, endo_nbr, Index_Equ_IM, Index_Var_IM);
  Free_IM(p_First_IM);
  count_Equ = count_Block = 0;
  if (*prologue)
    Allocate_Block(*prologue, &count_Equ, &count_Block, PROLOGUE, ModelBlock);
  for(j = 0;j < res->n_sets;j++)
    {
      if(res->sets_f[res->ordered[j]] == res->sets_s[res->ordered[j]])
        Allocate_Block(res->sets_f[res->ordered[j]] - res->sets_s[res->ordered[j]] + 1, &count_Equ, &count_Block, PROLOGUE, ModelBlock);
      else
        Allocate_Block(res->sets_f[res->ordered[j]] - res->sets_s[res->ordered[j]] + 1, &count_Equ, &count_Block, SIMULTANS, ModelBlock);
    }
  if (*epilogue)
    Allocate_Block(*epilogue, &count_Equ, &count_Block, EPILOGUE, ModelBlock);
  return 0;
}



//------------------------------------------------------------------------------
// Normalize each equation of the model (endgenous_i = f_i(endogenous_1, ..., endogenous_n) - in order to apply strong connex components search algorithm -
// and find the optimal blocks triangular decomposition
bool
BlockTriangular::Normalize_and_BlockDecompose(bool* IM, Model_Block* ModelBlock, int n, int* prologue, int* epilogue, simple* Index_Var_IM, simple* Index_Equ_IM, bool Do_Normalization, bool mixing, bool* IM_0)
{
  int i, j, Nb_TotalBlocks, Nb_RecursBlocks;
  int count_Block, count_Equ;
  block_result_t* res;
  List_IM * p_First_IM, *p_Cur_IM, *Cur_IM;
  Equation_set* Equation_gr = (Equation_set*) malloc(sizeof(Equation_set));
  p_First_IM = (List_IM*)malloc(sizeof(*p_First_IM));
  p_Cur_IM = p_First_IM;
  Cur_IM = First_IM;
  i = endo_nbr * endo_nbr;
  while(Cur_IM)
    {
      p_Cur_IM->lead_lag = Cur_IM->lead_lag;
      p_Cur_IM->IM = (bool*)malloc(i * sizeof(bool));
      memcpy ( p_Cur_IM->IM, Cur_IM->IM, i);
      Cur_IM = Cur_IM->pNext;
      if(Cur_IM)
        {
          p_Cur_IM->pNext = (List_IM*)malloc(sizeof(*p_Cur_IM));
          p_Cur_IM = p_Cur_IM->pNext;
        }
      else
        p_Cur_IM->pNext = NULL;
    }
  Prologue_Epilogue(IM, prologue, epilogue, n, Index_Var_IM, Index_Equ_IM);
  if(bt_verbose)
    {
      cout << "prologue : " << *prologue << " epilogue : " << *epilogue << "\n";
      Print_SIM(IM, n);
      for(i = 0;i < n;i++)
        cout << "Index_Var_IM[" << i << "]=" << Index_Var_IM[i].index << " Index_Equ_IM[" << i << "]=" << Index_Equ_IM[i].index << "\n";
    }
  if(Do_Normalization)
    {
      cout << "Normalizing the model ...\n";
      if(mixing)
        {
          bool* SIM0;
          SIM0 = (bool*)malloc(n * n * sizeof(bool));

          for(i = 0;i < n*n;i++)
            SIM0[i] = IM_0[i];
          for(i = 0;i < n;i++)
            for(j = 0;j < n;j++)
              IM_0[i*n + j] = SIM0[Index_Equ_IM[i].index * n + Index_Var_IM[j].index];
          free(SIM0);
          normalization.Normalize(n, *prologue, *epilogue, IM_0, Index_Equ_IM, Equation_gr, 1, IM);
        }
      else
        normalization.Normalize(n, *prologue, *epilogue, IM, Index_Equ_IM, Equation_gr, 0, 0);
    }
  else
    normalization.Gr_to_IM_basic(n, *prologue, *epilogue, IM, Equation_gr, false);
  cout << "Finding the optimal block decomposition of the model ...\n";
  if(bt_verbose)
    blocks.Print_Equation_gr(Equation_gr);
  res = blocks.sc(Equation_gr);
  if(bt_verbose)
    blocks.block_result_print(res);
  if ((*prologue) || (*epilogue))
    j = 1;
  else
    j = 0;
  for(i = 0;i < res->n_sets;i++)
    if ((res->sets_f[i] - res->sets_s[i] + 1) > j)
      j = res->sets_f[i] - res->sets_s[i] + 1;
  Nb_RecursBlocks = *prologue + *epilogue;
  Nb_TotalBlocks = res->n_sets + Nb_RecursBlocks;
  cout << Nb_TotalBlocks << " block(s) found:\n";
  cout << "  " << Nb_RecursBlocks << " recursive block(s) and " << res->n_sets << " simultaneous block(s). \n";
  cout << "  the largest simultaneous block has " << j << " equation(s). \n";
  ModelBlock->Size = Nb_TotalBlocks;
  ModelBlock->Periods = periods;
  ModelBlock->in_Block_Equ = (int*)malloc(n * sizeof(int));
  ModelBlock->in_Block_Var = (int*)malloc(n * sizeof(int));
  ModelBlock->in_Equ_of_Block = (int*)malloc(n * sizeof(int));
  ModelBlock->in_Var_of_Block = (int*)malloc(n * sizeof(int));
  ModelBlock->Block_List = (Block*)malloc(sizeof(ModelBlock->Block_List[0]) * Nb_TotalBlocks);
  blocks.block_result_to_IM(res, IM, *prologue, endo_nbr, Index_Equ_IM, Index_Var_IM);
  Free_IM(p_First_IM);
  count_Equ = count_Block = 0;
  if (*prologue)
    Allocate_Block(*prologue, &count_Equ, &count_Block, PROLOGUE, ModelBlock);
  for(j = 0;j < res->n_sets;j++)
    {
      if(res->sets_f[res->ordered[j]] == res->sets_s[res->ordered[j]])
        Allocate_Block(res->sets_f[res->ordered[j]] - res->sets_s[res->ordered[j]] + 1, &count_Equ, &count_Block, PROLOGUE, ModelBlock);
      else
        Allocate_Block(res->sets_f[res->ordered[j]] - res->sets_s[res->ordered[j]] + 1, &count_Equ, &count_Block, SIMULTANS, ModelBlock);
    }
  if (*epilogue)
    Allocate_Block(*epilogue, &count_Equ, &count_Block, EPILOGUE, ModelBlock);
  return 0;
}



//------------------------------------------------------------------------------
// For the contemparenous simultaneities
// normalize each equation of the model
// and find the optimal block triangular decomposition
void
BlockTriangular::Normalize_and_BlockDecompose_0()
{
  int i;
  List_IM* Cur_IM;
#ifdef DEBUG

  cout << "Normalize_and_BlockDecompose_0 \n";
#endif

  Index_Equ_IM = (simple*)malloc(endo_nbr * sizeof(*Index_Equ_IM));
  for(i = 0;i < endo_nbr;i++)
    {
      Index_Equ_IM[i].index = i;
      Index_Equ_IM[i].available = 1;
    }
  Index_Var_IM = (simple*)malloc(endo_nbr * sizeof(*Index_Var_IM));
  for(i = 0;i < endo_nbr;i++)
    {
      Index_Var_IM[i].index = i;
      Index_Var_IM[i].available = 1;
    }
  Cur_IM = Get_IM(0);
  Normalize_and_BlockDecompose(Cur_IM->IM, ModelBlock, endo_nbr, &prologue, &epilogue, Index_Var_IM, Index_Equ_IM, 1, 0, NULL);
}


//------------------------------------------------------------------------------
// normalize each equation of the statitic model
// and find the optimal block triangular decomposition
void
BlockTriangular::Normalize_and_BlockDecompose_Inside_Earth()
{
  bool *SIM, *SIM1;
  int i, j, blk, Size, Nb_blk, k, k1, k2, n0, n, k1_block;
  int prologue, epilogue;
  bool OK;
  simple *Index_Equ_IM, *Index_Var_IM;
  Model_Block *ModelBlock_Earth;
  //First create a the incidence matrix of each stack block
  n0 = endo_nbr;
  for(blk = 0;blk < ModelBlock->Size;blk++)
    {
      if(ModelBlock->Block_List[blk].Type == SIMULTANS)
        {
          Size = (n = ModelBlock->Block_List[blk].Size) * (Nb_blk = (max(ModelBlock->Block_List[blk].Max_Lag, ModelBlock->Block_List[blk].Max_Lead)) + 1);
          SIM = (bool*)malloc(Size * Size * sizeof(*SIM));
          memset (SIM, 0, Size*Size*sizeof(*SIM));
          k1_block = n * n * Nb_blk;
          for(k1 = 0;k1 < Nb_blk;k1++)
            {
              for(k2 = 0;k2 < Nb_blk;k2++)
                {
                  k = k2 - k1;
                  OK = 1;
                  if(k < 0)
                    if( -k > ModelBlock->Block_List[blk].Max_Lag)
                      OK = 0;
                  if(k > 0)
                    if(k > ModelBlock->Block_List[blk].Max_Lead)
                      OK = 0;
                  if(OK)
                    {
                      SIM1 = Get_IM(k)->IM;
                      for(i = 0;i < n;i++)
                        for(j = 0;j < n;j++)
                          SIM[i*Size + k1*k1_block + j + k2*n] = SIM1[ModelBlock->Block_List[blk].Equation[i] * n0 + ModelBlock->Block_List[blk].Variable[j]];
                    }
                }
            }
          if(bt_verbose)
            cout << "incidence matrix for the static model (unsorted) \n";
          Index_Equ_IM = (simple*)malloc(Size * sizeof(*Index_Equ_IM));
          for(i = 0;i < Size;i++)
            {
              Index_Equ_IM[i].index = i;
              Index_Equ_IM[i].available = 1;
            }
          Index_Var_IM = (simple*)malloc(Size * sizeof(*Index_Var_IM));
          for(i = 0;i < Size;i++)
            {
              Index_Var_IM[i].index = i;
              Index_Var_IM[i].available = 1;
            }
          ModelBlock_Earth = (Model_Block*)malloc(sizeof(Model_Block));
          Normalize_and_BlockDecompose(SIM, ModelBlock_Earth, Size, &prologue, &epilogue, Index_Var_IM, Index_Equ_IM, 0, 0, NULL);
        }
    }
}

//------------------------------------------------------------------------------
// normalize each equation of the statitic model
// and find the optimal block triangular decomposition
void
BlockTriangular::Normalize_and_BlockDecompose_Static_Model()
{
  bool* SIM;
  List_IM* Cur_IM;
  int i;
  //First create a static model incidence matrix
  SIM = (bool*)malloc(endo_nbr * endo_nbr * sizeof(*SIM));
  for(i = 0;i < endo_nbr*endo_nbr;i++)
    {
      SIM[i] = 0;
      Cur_IM = First_IM;
      while(Cur_IM)
        {
          SIM[i] = (SIM[i]) || (Cur_IM->IM[i]);
          Cur_IM = Cur_IM->pNext;
        }
    }
  if(bt_verbose)
    {
      cout << "incidence matrix for the static model (unsorted) \n";
      Print_SIM(SIM, endo_nbr);
    }
  Index_Equ_IM = (simple*)malloc(endo_nbr * sizeof(*Index_Equ_IM));
  for(i = 0;i < endo_nbr;i++)
    {
      Index_Equ_IM[i].index = i;
      Index_Equ_IM[i].available = 1;
    }
  Index_Var_IM = (simple*)malloc(endo_nbr * sizeof(*Index_Var_IM));
  for(i = 0;i < endo_nbr;i++)
    {
      Index_Var_IM[i].index = i;
      Index_Var_IM[i].available = 1;
    }
  if(ModelBlock != NULL)
    Free_Block(ModelBlock);
  ModelBlock = (Model_Block*)malloc(sizeof(*ModelBlock));
  Normalize_and_BlockDecompose(SIM, ModelBlock, endo_nbr, &prologue, &epilogue, Index_Var_IM, Index_Equ_IM, 1, 0, NULL);
}


//------------------------------------------------------------------------------
// normalize each equation of the dynamic model
// and find the optimal block triangular decomposition of the static model
void
BlockTriangular::Normalize_and_BlockDecompose_Static_0_Model(const jacob_map &j_m)
{
  bool* SIM, *SIM_0;
  List_IM* Cur_IM;
  int i;
  //First create a static model incidence matrix
  SIM = (bool*)malloc(endo_nbr * endo_nbr * sizeof(*SIM));
  for(i = 0;i < endo_nbr*endo_nbr;i++)
    {
      SIM[i] = 0;
      Cur_IM = First_IM;
      while(Cur_IM)
        {
          SIM[i] = (SIM[i]) || (Cur_IM->IM[i]);
          Cur_IM = Cur_IM->pNext;
        }
    }
  if(bt_verbose)
    {
      cout << "incidence matrix for the static model (unsorted) \n";
      Print_SIM(SIM, endo_nbr);
    }
  Index_Equ_IM = (simple*)malloc(endo_nbr * sizeof(*Index_Equ_IM));
  for(i = 0;i < endo_nbr;i++)
    {
      Index_Equ_IM[i].index = i;
      Index_Equ_IM[i].available = 1;
    }
  Index_Var_IM = (simple*)malloc(endo_nbr * sizeof(*Index_Var_IM));
  for(i = 0;i < endo_nbr;i++)
    {
      Index_Var_IM[i].index = i;
      Index_Var_IM[i].available = 1;
    }
  if(ModelBlock != NULL)
    Free_Block(ModelBlock);
  ModelBlock = (Model_Block*)malloc(sizeof(*ModelBlock));
  Cur_IM = Get_IM(0);
  SIM_0 = (bool*)malloc(endo_nbr * endo_nbr * sizeof(*SIM_0));
  for(i = 0;i < endo_nbr*endo_nbr;i++)
    SIM_0[i] = Cur_IM->IM[i];
  Normalize_and_BlockDecompose(SIM, ModelBlock, endo_nbr, &prologue, &epilogue, Index_Var_IM, Index_Equ_IM, 1, 1, SIM_0, j_m);
  if(bt_verbose)
    for(i = 0;i < endo_nbr;i++)
      cout << "Block=" << Index_Equ_IM[i].block << " Equ=" << Index_Equ_IM[i].index << " Var= " << Index_Var_IM[i].index << "  " << symbol_table.getNameByID(eEndogenous, Index_Var_IM[i].index) << "\n";
}
