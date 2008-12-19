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
  normalization(symbol_table_arg),
  incidencematrix(symbol_table_arg)
{
  bt_verbose = 0;
  ModelBlock = NULL;
  periods = 0;
}



//------------------------------------------------------------------------------
// Find the prologue and the epilogue of the model
void
BlockTriangular::Prologue_Epilogue(bool* IM, int* prologue, int* epilogue, int n, simple* Index_Var_IM, simple* Index_Equ_IM, bool* IM0)
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
          if ((k == 1) && IM0[Index_Equ_IM[i].index*n + Index_Var_IM[l].index])
            {
              modifie = 1;
              incidencematrix.swap_IM_c(IM, *prologue, i, l, Index_Var_IM, Index_Equ_IM, n);
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
          if ((k == 1) && IM0[Index_Equ_IM[l].index*n + Index_Var_IM[i].index])
            {
              modifie = 1;
              incidencematrix.swap_IM_c(IM, n - (1 + *epilogue), l, i, Index_Var_IM, Index_Equ_IM, n);
              (*epilogue)++;
            }
        }
    }
}


void
BlockTriangular::Allocate_Block(int size, int *count_Equ, int *count_Block, BlockType type, Model_Block * ModelBlock)
{
  int i, j, k, l, ls, m, i_1, Lead, Lag, first_count_equ, i1, li;
  int *tmp_size, *tmp_size_other_endo, *tmp_size_exo, *tmp_var, *tmp_endo, *tmp_other_endo, *tmp_exo, tmp_nb_other_endo, tmp_nb_exo, nb_lead_lag_endo;
  bool *tmp_variable_evaluated;
  bool *Cur_IM;
  bool *IM, OK;
  ModelBlock->Periods = periods;
  int Lag_Endo, Lead_Endo, Lag_Exo, Lead_Exo, Lag_Other_Endo, Lead_Other_Endo;
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
          tmp_endo = (int*)malloc((incidencematrix.Model_Max_Lead + incidencematrix.Model_Max_Lag + 1) * sizeof(int));
          tmp_other_endo = (int*)malloc(symbol_table.endo_nbr * sizeof(int));
          tmp_size_other_endo = (int*)malloc((incidencematrix.Model_Max_Lead + incidencematrix.Model_Max_Lag + 1) * sizeof(int));
          tmp_size = (int*)malloc((incidencematrix.Model_Max_Lead + incidencematrix.Model_Max_Lag + 1) * sizeof(int));
          tmp_var = (int*)malloc(sizeof(int));
          tmp_size_exo = (int*)malloc((incidencematrix.Model_Max_Lead + incidencematrix.Model_Max_Lag + 1) * sizeof(int));

          memset(tmp_size_exo, 0, (incidencematrix.Model_Max_Lead + incidencematrix.Model_Max_Lag + 1)*sizeof(int));
          memset(tmp_size, 0, (incidencematrix.Model_Max_Lead + incidencematrix.Model_Max_Lag + 1)*sizeof(int));
          memset(tmp_endo, 0, (incidencematrix.Model_Max_Lead + incidencematrix.Model_Max_Lag + 1)*sizeof(int));
          memset(tmp_size_other_endo, 0, (incidencematrix.Model_Max_Lead + incidencematrix.Model_Max_Lag + 1)*sizeof(int));
          memset(tmp_other_endo, 0, symbol_table.endo_nbr*sizeof(int));


          nb_lead_lag_endo = Lead = Lag = 0;
          Lag_Endo = Lead_Endo = Lag_Exo = Lead_Exo = Lag_Other_Endo = Lead_Other_Endo = 0;

          tmp_variable_evaluated = (bool*)malloc(symbol_table.endo_nbr*sizeof(bool));
          memset(tmp_variable_evaluated, 0, symbol_table.endo_nbr*sizeof(bool));
          for(k = -incidencematrix.Model_Max_Lag_Endo; k<=incidencematrix.Model_Max_Lead_Endo; k++)
            {
              Cur_IM = incidencematrix.Get_IM(k, eEndogenous);
              if(Cur_IM)
                {
                  i_1 = Index_Equ_IM[*count_Equ].index * symbol_table.endo_nbr;
                  if(k > 0)
                    {
                      if(Cur_IM[i_1 + Index_Var_IM[*count_Equ].index])
                        {
                          tmp_variable_evaluated[Index_Var_IM[*count_Equ].index] = true;
                          nb_lead_lag_endo++;
                          tmp_size[incidencematrix.Model_Max_Lag_Endo + k]++;
                          if(k > Lead)
                            Lead = k;
                        }
                    }
                  else
                    {
                      if(Cur_IM[i_1 + Index_Var_IM[*count_Equ].index])
                        {
                          tmp_variable_evaluated[Index_Var_IM[*count_Equ].index] = true;
                          tmp_size[incidencematrix.Model_Max_Lag_Endo + k]++;
                          nb_lead_lag_endo++;
                          if(-k > Lag)
                            Lag = -k;
                        }
                    }
                }
            }

          Lag_Endo = Lag;
          Lead_Endo = Lead;

          tmp_nb_other_endo = 0;
          for(k = -incidencematrix.Model_Max_Lag_Endo; k<=incidencematrix.Model_Max_Lead_Endo; k++)
            {
              Cur_IM = incidencematrix.Get_IM(k, eEndogenous);
              if(Cur_IM)
                {
                  i_1 = Index_Equ_IM[*count_Equ].index * symbol_table.endo_nbr;
                  for(j = 0;j < symbol_table.endo_nbr;j++)
                    {
                      int ij = Index_Var_IM[j].index;
                      if(Cur_IM[i_1 + ij])
                        {
                          if(!tmp_variable_evaluated[ij])
                            {
                              if(!tmp_other_endo[ij])
                                {
                                  tmp_other_endo[ij] = 1;
                                  tmp_nb_other_endo++;
                                }
                              if(k>0 && k>Lead_Other_Endo)
                                Lead_Other_Endo = k;
                              else if(k<0 && (-k)>Lag_Other_Endo)
                                Lag_Other_Endo = -k;
                              if(k>0 && k>Lead)
                                Lead = k;
                              else if(k<0 && (-k)>Lag)
                               Lag = -k;
                              tmp_size_other_endo[k+incidencematrix.Model_Max_Lag_Endo]++;
                            }
                        }
                    }
                }
            }
          ModelBlock->Block_List[*count_Block].nb_other_endo = tmp_nb_other_endo;
          ModelBlock->Block_List[*count_Block].Other_Endogenous = (int*)malloc(tmp_nb_other_endo * sizeof(int));



          tmp_exo = (int*)malloc(symbol_table.exo_nbr * sizeof(int));
          memset(tmp_exo, 0, symbol_table.exo_nbr * sizeof(int));
          tmp_nb_exo = 0;

          for(k = -incidencematrix.Model_Max_Lag_Exo; k<=incidencematrix.Model_Max_Lead_Exo; k++)
            {
              Cur_IM = incidencematrix.Get_IM(k, eExogenous);
              if(Cur_IM)
                {
                  i_1 = Index_Equ_IM[*count_Equ].index * symbol_table.exo_nbr;
                  for(j=0;j<symbol_table.exo_nbr;j++)
                    if(Cur_IM[i_1 + j])
                      {
                        if(!tmp_exo[j])
                          {
                            tmp_exo[j] = 1;
                            tmp_nb_exo++;
                          }
                        if(k>0 && k>Lead_Exo)
                          Lead_Exo = k;
                        else if(k<0 && (-k)>Lag_Exo)
                          Lag_Exo = -k;
                        if(k>0 && k>Lead)
                          Lead = k;
                        else if(k<0 && (-k)>Lag)
                          Lag = -k;
                        tmp_size_exo[k+incidencematrix.Model_Max_Lag_Exo]++;
                      }
                }
            }

          ModelBlock->Block_List[*count_Block].nb_exo = tmp_nb_exo;
          ModelBlock->Block_List[*count_Block].Exogenous = (int*)malloc(tmp_nb_exo * sizeof(int));
          k = 0;
          for(j=0;j<symbol_table.exo_nbr;j++)
            if(tmp_exo[j])
              {
                ModelBlock->Block_List[*count_Block].Exogenous[k] = j;
                k++;
              }

          ModelBlock->Block_List[*count_Block].nb_exo_det = 0;

          ModelBlock->Block_List[*count_Block].Max_Lag = Lag;
          ModelBlock->Block_List[*count_Block].Max_Lead = Lead;
          ModelBlock->Block_List[*count_Block].Max_Lag_Endo = Lag_Endo;
          ModelBlock->Block_List[*count_Block].Max_Lead_Endo = Lead_Endo;
          ModelBlock->Block_List[*count_Block].Max_Lag_Exo = Lag_Exo;
          ModelBlock->Block_List[*count_Block].Max_Lead_Exo = Lead_Exo;
          ModelBlock->Block_List[*count_Block].Equation = (int*)malloc(sizeof(int));
          ModelBlock->Block_List[*count_Block].Variable = (int*)malloc(sizeof(int));
          ModelBlock->Block_List[*count_Block].Own_Derivative = (int*)malloc(sizeof(int));
          ModelBlock->Block_List[*count_Block].Equation[0] = Index_Equ_IM[*count_Equ].index;
          ModelBlock->Block_List[*count_Block].Variable[0] = Index_Var_IM[*count_Equ].index;



          if ((Lead > 0) && (Lag > 0))
            ModelBlock->Block_List[*count_Block].Simulation_Type = SOLVE_TWO_BOUNDARIES_SIMPLE;
          else if((Lead > 0) && (Lag == 0))
            ModelBlock->Block_List[*count_Block].Simulation_Type = SOLVE_BACKWARD_SIMPLE;
          else
            ModelBlock->Block_List[*count_Block].Simulation_Type = SOLVE_FORWARD_SIMPLE;

          tmp_exo = (int*)malloc(symbol_table.exo_nbr * sizeof(int));
          memset(tmp_exo, 0, symbol_table.exo_nbr * sizeof(int));
          tmp_nb_exo = 0;
          for(k = -incidencematrix.Model_Max_Lag_Exo; k <=incidencematrix.Model_Max_Lead_Exo; k++)
            {
              Cur_IM = incidencematrix.Get_IM(k, eExogenous);
              if(Cur_IM)
                {
                  i_1 = Index_Equ_IM[*count_Equ].index * symbol_table.exo_nbr;
                  for(j=0;j<symbol_table.exo_nbr;j++)
                    if(Cur_IM[i_1 + j] && (!tmp_exo[j]))
                      {
                        tmp_exo[j] = 1;
                        tmp_nb_exo++;
                      }
                }
            }
          ModelBlock->Block_List[*count_Block].nb_exo = tmp_nb_exo;
          ModelBlock->Block_List[*count_Block].Exogenous = (int*)malloc(tmp_nb_exo * sizeof(int));

          ModelBlock->Block_List[*count_Block].IM_lead_lag = (IM_compact*)malloc((Lead + Lag + 1) * sizeof(IM_compact));
          ModelBlock->Block_List[*count_Block].Nb_Lead_Lag_Endo = nb_lead_lag_endo;


          k = 0;
          for(j=0;j<symbol_table.exo_nbr;j++)
            if(tmp_exo[j])
              {
                ModelBlock->Block_List[*count_Block].Exogenous[k] = j;
                k++;
              }
          ls = l = 1;
          i1 = 0;
          for(int li = 0;li < Lead + Lag + 1;li++)
            {
              if(incidencematrix.Model_Max_Lag_Endo - Lag + li>=0)
                {
                  ModelBlock->Block_List[*count_Block].IM_lead_lag[li].size = tmp_size[incidencematrix.Model_Max_Lag_Endo - Lag + li];
                  ModelBlock->Block_List[*count_Block].IM_lead_lag[li].nb_endo = tmp_size[incidencematrix.Model_Max_Lag_Endo - Lag + li];
                  ModelBlock->Block_List[*count_Block].IM_lead_lag[li].u = (int*)malloc(tmp_size[incidencematrix.Model_Max_Lag_Endo - Lag + li] * sizeof(int));
                  ModelBlock->Block_List[*count_Block].IM_lead_lag[li].us = (int*)malloc(tmp_size[incidencematrix.Model_Max_Lag_Endo - Lag + li] * sizeof(int));
                  ModelBlock->Block_List[*count_Block].IM_lead_lag[li].Var = (int*)malloc(tmp_size[incidencematrix.Model_Max_Lag_Endo - Lag + li] * sizeof(int));
                  ModelBlock->Block_List[*count_Block].IM_lead_lag[li].Equ = (int*)malloc(tmp_size[incidencematrix.Model_Max_Lag_Endo - Lag + li] * sizeof(int));
                  ModelBlock->Block_List[*count_Block].IM_lead_lag[li].Var_Index = (int*)malloc(tmp_size[incidencematrix.Model_Max_Lag_Endo - Lag + li] * sizeof(int));
                  ModelBlock->Block_List[*count_Block].IM_lead_lag[li].Equ_Index = (int*)malloc(tmp_size[incidencematrix.Model_Max_Lag_Endo - Lag + li] * sizeof(int));

                  ModelBlock->Block_List[*count_Block].IM_lead_lag[li].size_other_endo = tmp_size_other_endo[incidencematrix.Model_Max_Lag_Endo - Lag + li];
                  ModelBlock->Block_List[*count_Block].IM_lead_lag[li].nb_other_endo = tmp_other_endo[incidencematrix.Model_Max_Lag_Endo - Lag + li];
                  ModelBlock->Block_List[*count_Block].IM_lead_lag[li].u_other_endo = (int*)malloc(tmp_size_other_endo[incidencematrix.Model_Max_Lag_Endo - Lag + li] * sizeof(int));
                  ModelBlock->Block_List[*count_Block].IM_lead_lag[li].Var_other_endo = (int*)malloc(tmp_size_other_endo[incidencematrix.Model_Max_Lag_Endo - Lag + li] * sizeof(int));
                  ModelBlock->Block_List[*count_Block].IM_lead_lag[li].Equ_other_endo = (int*)malloc(tmp_size_other_endo[incidencematrix.Model_Max_Lag_Endo - Lag + li] * sizeof(int));
                  ModelBlock->Block_List[*count_Block].IM_lead_lag[li].Var_Index_other_endo = (int*)malloc(tmp_size_other_endo[incidencematrix.Model_Max_Lag_Endo - Lag + li] * sizeof(int));
                  ModelBlock->Block_List[*count_Block].IM_lead_lag[li].Equ_Index_other_endo = (int*)malloc(tmp_size_other_endo[incidencematrix.Model_Max_Lag_Endo - Lag + li] * sizeof(int));



                  ModelBlock->Block_List[*count_Block].IM_lead_lag[li].u_init = l;
                  IM = incidencematrix.Get_IM(li - Lag, eEndogenous);
                  if(IM)
                    {
                      if(IM[Index_Var_IM[*count_Equ].index + Index_Equ_IM[*count_Equ].index*symbol_table.endo_nbr] && nb_lead_lag_endo)
                        {
                          tmp_var[0] = i1;
                          i1++;
                        }
                      m = 0;
                      i_1 = Index_Equ_IM[*count_Equ].index * symbol_table.endo_nbr;
                      if(IM[Index_Var_IM[*count_Equ].index + i_1])
                        {
                          if(li == Lag)
                            {
                              ModelBlock->Block_List[*count_Block].IM_lead_lag[li].us[m] = ls;
                              ls++;
                            }
                          ModelBlock->Block_List[*count_Block].IM_lead_lag[li].u[m] = l;
                          ModelBlock->Block_List[*count_Block].IM_lead_lag[li].Equ[m] = 0;
                          ModelBlock->Block_List[*count_Block].IM_lead_lag[li].Var[m] = 0;
                          ModelBlock->Block_List[*count_Block].IM_lead_lag[li].Equ_Index[m] = Index_Equ_IM[*count_Equ].index;
                          ModelBlock->Block_List[*count_Block].IM_lead_lag[li].Var_Index[m] = Index_Var_IM[*count_Equ].index;
                          tmp_variable_evaluated[Index_Var_IM[*count_Equ].index] = true;
                          l++;
                          m++;
                        }
                      ModelBlock->Block_List[*count_Block].IM_lead_lag[li].u_finish = l - 1;

                      m = 0;
                      for(k = 0;k < symbol_table.endo_nbr;k++)
                        if((!tmp_variable_evaluated[Index_Var_IM[k].index]) && IM[Index_Var_IM[k].index + i_1])
                          {
                            ModelBlock->Block_List[*count_Block].IM_lead_lag[li].u_other_endo[m] = l;
                            ModelBlock->Block_List[*count_Block].IM_lead_lag[li].Equ_other_endo[m] = 0; //j - first_count_equ;
                            ModelBlock->Block_List[*count_Block].IM_lead_lag[li].Var_other_endo[m] = k ;
                            ModelBlock->Block_List[*count_Block].IM_lead_lag[li].Equ_Index_other_endo[m] = Index_Equ_IM[*count_Equ].index;
                            ModelBlock->Block_List[*count_Block].IM_lead_lag[li].Var_Index_other_endo[m] = Index_Var_IM[k].index;
                            l++;
                            m++;
                          }
                    }
                }
              else
                ModelBlock->Block_List[*count_Block].IM_lead_lag[li].size = 0;
              if(incidencematrix.Model_Max_Lag_Exo - Lag + li>=0)
                {
                  ModelBlock->Block_List[*count_Block].IM_lead_lag[li].size_exo = tmp_size_exo[incidencematrix.Model_Max_Lag_Exo - Lag + li];
                  ModelBlock->Block_List[*count_Block].IM_lead_lag[li].Exogenous = (int*)malloc(tmp_size_exo[incidencematrix.Model_Max_Lag_Exo - Lag + li] * sizeof(int));
                  ModelBlock->Block_List[*count_Block].IM_lead_lag[li].Exogenous_Index = (int*)malloc(tmp_size_exo[incidencematrix.Model_Max_Lag_Exo - Lag + li] * sizeof(int));
                  ModelBlock->Block_List[*count_Block].IM_lead_lag[li].Equ_X = (int*)malloc(tmp_size_exo[incidencematrix.Model_Max_Lag_Exo - Lag + li] * sizeof(int));
                  ModelBlock->Block_List[*count_Block].IM_lead_lag[li].Equ_X_Index = (int*)malloc(tmp_size_exo[incidencematrix.Model_Max_Lag_Exo - Lag + li] * sizeof(int));
                  IM = incidencematrix.Get_IM(li - Lag, eExogenous);
                  if(IM)
                    {
                      m = 0;
                      i_1 = Index_Equ_IM[*count_Equ].index * symbol_table.exo_nbr;
                      for(k = 0; k<tmp_nb_exo; k++)
                        {
                          if(IM[ModelBlock->Block_List[*count_Block].Exogenous[k]+i_1])
                            {
                              ModelBlock->Block_List[*count_Block].IM_lead_lag[li].Exogenous[m] = k;
                              ModelBlock->Block_List[*count_Block].IM_lead_lag[li].Exogenous_Index[m] = ModelBlock->Block_List[*count_Block].Exogenous[k];
                              ModelBlock->Block_List[*count_Block].IM_lead_lag[li].Equ_X[m] = 0;
                              ModelBlock->Block_List[*count_Block].IM_lead_lag[li].Equ_X_Index[m] = Index_Equ_IM[*count_Equ].index;
                              m++;
                            }
                        }
                    }
                }
              else
                ModelBlock->Block_List[*count_Block].IM_lead_lag[li].size_exo = 0;
            }
          (*count_Equ)++;
          (*count_Block)++;
          free(tmp_size);
          free(tmp_size_exo);
          free(tmp_endo);
          free(tmp_exo);
          free(tmp_var);
          free(tmp_size_other_endo);
          free(tmp_other_endo);
          free(tmp_variable_evaluated);
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
      ModelBlock->Block_List[*count_Block].Own_Derivative = (int*)malloc(ModelBlock->Block_List[*count_Block].Size * sizeof(int));
      Lead = Lag = 0;
      first_count_equ = *count_Equ;
      tmp_var = (int*)malloc(size * sizeof(int));
      tmp_endo = (int*)malloc((incidencematrix.Model_Max_Lead + incidencematrix.Model_Max_Lag + 1) * sizeof(int));
      tmp_other_endo = (int*)malloc(symbol_table.endo_nbr * sizeof(int));
      tmp_size = (int*)malloc((incidencematrix.Model_Max_Lead + incidencematrix.Model_Max_Lag + 1) * sizeof(int));
      tmp_size_other_endo = (int*)malloc((incidencematrix.Model_Max_Lead + incidencematrix.Model_Max_Lag + 1) * sizeof(int));
      tmp_size_exo = (int*)malloc((incidencematrix.Model_Max_Lead + incidencematrix.Model_Max_Lag + 1) * sizeof(int));
      memset(tmp_size_exo, 0, (incidencematrix.Model_Max_Lead + incidencematrix.Model_Max_Lag + 1)*sizeof(int));
      memset(tmp_size_other_endo, 0, (incidencematrix.Model_Max_Lead + incidencematrix.Model_Max_Lag + 1)*sizeof(int));
      memset(tmp_size, 0, (incidencematrix.Model_Max_Lead + incidencematrix.Model_Max_Lag + 1)*sizeof(int));
      memset(tmp_endo, 0, (incidencematrix.Model_Max_Lead + incidencematrix.Model_Max_Lag + 1)*sizeof(int));
      memset(tmp_other_endo, 0, symbol_table.endo_nbr*sizeof(int));
      nb_lead_lag_endo = 0;
      Lag_Endo = Lead_Endo = Lag_Other_Endo = Lead_Other_Endo = Lag_Exo = Lead_Exo = 0;

      //Variable by variable looking for all leads and lags its occurence in each equation of the block
      tmp_variable_evaluated = (bool*)malloc(symbol_table.endo_nbr*sizeof(bool));
      memset(tmp_variable_evaluated, 0, symbol_table.endo_nbr*sizeof(bool));
      for(i = 0;i < size;i++)
        {
          ModelBlock->Block_List[*count_Block].Equation[i] = Index_Equ_IM[*count_Equ].index;
          ModelBlock->Block_List[*count_Block].Variable[i] = Index_Var_IM[*count_Equ].index;
          i_1 = Index_Var_IM[*count_Equ].index;
          for(k = -incidencematrix.Model_Max_Lag_Endo; k<=incidencematrix.Model_Max_Lead_Endo; k++)
            {
              Cur_IM = incidencematrix.Get_IM(k, eEndogenous);
              if(Cur_IM)
                {
                  OK = false;
                  if(k >= 0)
                    {
                      for(j = 0;j < size;j++)
                        {
                          if(Cur_IM[i_1 + Index_Equ_IM[first_count_equ + j].index*symbol_table.endo_nbr])
                            {
                              tmp_variable_evaluated[i_1] = true;
                              tmp_size[incidencematrix.Model_Max_Lag_Endo + k]++;
                              if (!OK)
                                {
                                  tmp_endo[incidencematrix.Model_Max_Lag + k]++;
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
                      for(j = 0;j < size;j++)
                        {
                          if(Cur_IM[i_1 + Index_Equ_IM[first_count_equ + j].index*symbol_table.endo_nbr])
                            {
                              tmp_size[incidencematrix.Model_Max_Lag_Endo + k]++;
                              if (!OK)
                                {
                                  tmp_variable_evaluated[i_1] = true;
                                  tmp_endo[incidencematrix.Model_Max_Lag + k]++;
                                  nb_lead_lag_endo++;
                                  OK = true;
                                }
                              if(-k > Lag)
                                Lag = -k;
                            }
                        }
                    }
               }
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
            ModelBlock->Block_List[*count_Block].Simulation_Type = SOLVE_FORWARD_COMPLETE;
        }
      else
        {
          if(Lead > 0)
            ModelBlock->Block_List[*count_Block].Simulation_Type = SOLVE_BACKWARD_SIMPLE;
          else
            ModelBlock->Block_List[*count_Block].Simulation_Type = SOLVE_FORWARD_SIMPLE;
        }
      Lag_Endo = Lag;
      Lead_Endo = Lead;

      tmp_nb_other_endo = 0;
      for(i = 0;i < size;i++)
        {
          for(k = -incidencematrix.Model_Max_Lag_Endo; k<=incidencematrix.Model_Max_Lead_Endo; k++)
            {
              Cur_IM = incidencematrix.Get_IM(k, eEndogenous);
              if(Cur_IM)
                {
                  i_1 = Index_Equ_IM[first_count_equ+i].index * symbol_table.endo_nbr;
                  for(j = 0;j < symbol_table.endo_nbr;j++)
                    if(Cur_IM[i_1 + j])
                      {
                        if(!tmp_variable_evaluated[j])
                          {
                            tmp_other_endo[j] = 1;
                            tmp_nb_other_endo++;
                          }
                        if(k>0 && k>Lead_Other_Endo)
                          Lead_Other_Endo = k;
                        else if(k<0 && (-k)>Lag_Other_Endo)
                          Lag_Other_Endo = -k;
                        if(k>0 && k>Lead)
                          Lead = k;
                        else if(k<0 && (-k)>Lag)
                          Lag = -k;
                        tmp_size_other_endo[k+incidencematrix.Model_Max_Lag_Endo]++;
                      }
                }
            }
        }
      ModelBlock->Block_List[*count_Block].nb_other_endo = tmp_nb_other_endo;
      ModelBlock->Block_List[*count_Block].Other_Endogenous = (int*)malloc(tmp_nb_other_endo * sizeof(int));


      tmp_exo = (int*)malloc(symbol_table.exo_nbr * sizeof(int));
      memset(tmp_exo, 0, symbol_table.exo_nbr *     sizeof(int));
      tmp_nb_exo = 0;
      for(i = 0;i < size;i++)
        {
          for(k = -incidencematrix.Model_Max_Lag_Exo; k<=incidencematrix.Model_Max_Lead_Exo; k++)
            {
              Cur_IM = incidencematrix.Get_IM(k, eExogenous);
              if(Cur_IM)
                {
                  i_1 = Index_Equ_IM[first_count_equ+i].index * symbol_table.exo_nbr;
                  for(j=0;j<symbol_table.exo_nbr;j++)
                    if(Cur_IM[i_1 + j])
                      {
                        if(!tmp_exo[j])
                          {
                            tmp_exo[j] = 1;
                            tmp_nb_exo++;
                          }
                        if(k>0 && k>Lead_Exo)
                          Lead_Exo = k;
                        else if(k<0 && (-k)>Lag_Exo)
                          Lag_Exo = -k;
                        if(k>0 && k>Lead)
                          Lead = k;
                        else if(k<0 && (-k)>Lag)
                          Lag = -k;
                        tmp_size_exo[k+incidencematrix.Model_Max_Lag_Exo]++;
                      }
                }
            }
        }


      ModelBlock->Block_List[*count_Block].nb_exo = tmp_nb_exo;
      ModelBlock->Block_List[*count_Block].Exogenous = (int*)malloc(tmp_nb_exo * sizeof(int));
      k = 0;
      for(j=0;j<symbol_table.exo_nbr;j++)
        if(tmp_exo[j])
          {
            ModelBlock->Block_List[*count_Block].Exogenous[k] = j;
            k++;
          }

      ModelBlock->Block_List[*count_Block].nb_exo_det = 0;

      ModelBlock->Block_List[*count_Block].Max_Lag = Lag;
      ModelBlock->Block_List[*count_Block].Max_Lead = Lead;
      ModelBlock->Block_List[*count_Block].Max_Lag_Endo = Lag_Endo;
      ModelBlock->Block_List[*count_Block].Max_Lead_Endo = Lead_Endo;
      ModelBlock->Block_List[*count_Block].Max_Lag_Other_Endo = Lag_Other_Endo;
      ModelBlock->Block_List[*count_Block].Max_Lead_Other_Endo = Lead_Other_Endo;
      ModelBlock->Block_List[*count_Block].Max_Lag_Exo = Lag_Exo;
      ModelBlock->Block_List[*count_Block].Max_Lead_Exo = Lead_Exo;
      ModelBlock->Block_List[*count_Block].IM_lead_lag = (IM_compact*)malloc((Lead + Lag + 1) * sizeof(IM_compact));
      ls = l = li = size;
      i1 = 0;
      ModelBlock->Block_List[*count_Block].Nb_Lead_Lag_Endo = nb_lead_lag_endo;
      for(i = 0;i < Lead + Lag + 1;i++)
        {
          if(incidencematrix.Model_Max_Lag_Endo-Lag+i>=0)
            {
              ModelBlock->Block_List[*count_Block].IM_lead_lag[i].size = tmp_size[incidencematrix.Model_Max_Lag_Endo - Lag + i];
              ModelBlock->Block_List[*count_Block].IM_lead_lag[i].nb_endo = tmp_endo[incidencematrix.Model_Max_Lag_Endo - Lag + i];
              ModelBlock->Block_List[*count_Block].IM_lead_lag[i].u = (int*)malloc(tmp_size[incidencematrix.Model_Max_Lag_Endo - Lag + i] * sizeof(int));
              ModelBlock->Block_List[*count_Block].IM_lead_lag[i].us = (int*)malloc(tmp_size[incidencematrix.Model_Max_Lag_Endo - Lag + i] * sizeof(int));
              ModelBlock->Block_List[*count_Block].IM_lead_lag[i].Var = (int*)malloc(tmp_size[incidencematrix.Model_Max_Lag_Endo - Lag + i] * sizeof(int));
              ModelBlock->Block_List[*count_Block].IM_lead_lag[i].Equ = (int*)malloc(tmp_size[incidencematrix.Model_Max_Lag_Endo - Lag + i] * sizeof(int));
              ModelBlock->Block_List[*count_Block].IM_lead_lag[i].Var_Index = (int*)malloc(tmp_size[incidencematrix.Model_Max_Lag_Endo - Lag + i] * sizeof(int));
              ModelBlock->Block_List[*count_Block].IM_lead_lag[i].Equ_Index = (int*)malloc(tmp_size[incidencematrix.Model_Max_Lag_Endo - Lag + i] * sizeof(int));

              ModelBlock->Block_List[*count_Block].IM_lead_lag[i].size_other_endo = tmp_size_other_endo[incidencematrix.Model_Max_Lag_Endo - Lag + i];
              ModelBlock->Block_List[*count_Block].IM_lead_lag[i].nb_other_endo = tmp_other_endo[incidencematrix.Model_Max_Lag_Endo - Lag + i];
              ModelBlock->Block_List[*count_Block].IM_lead_lag[i].u_other_endo = (int*)malloc(tmp_size_other_endo[incidencematrix.Model_Max_Lag_Endo - Lag + i] * sizeof(int));
              ModelBlock->Block_List[*count_Block].IM_lead_lag[i].Var_other_endo = (int*)malloc(tmp_size_other_endo[incidencematrix.Model_Max_Lag_Endo - Lag + i] * sizeof(int));
              ModelBlock->Block_List[*count_Block].IM_lead_lag[i].Equ_other_endo = (int*)malloc(tmp_size_other_endo[incidencematrix.Model_Max_Lag_Endo - Lag + i] * sizeof(int));
              ModelBlock->Block_List[*count_Block].IM_lead_lag[i].Var_Index_other_endo = (int*)malloc(tmp_size_other_endo[incidencematrix.Model_Max_Lag_Endo - Lag + i] * sizeof(int));
              ModelBlock->Block_List[*count_Block].IM_lead_lag[i].Equ_Index_other_endo = (int*)malloc(tmp_size_other_endo[incidencematrix.Model_Max_Lag_Endo - Lag + i] * sizeof(int));
            }
          else
            ModelBlock->Block_List[*count_Block].IM_lead_lag[i].size = 0;
          if(incidencematrix.Model_Max_Lag_Exo-Lag+i>=0)
            {
              ModelBlock->Block_List[*count_Block].IM_lead_lag[i].size_exo = tmp_size_exo[incidencematrix.Model_Max_Lag_Exo - Lag + i];
              ModelBlock->Block_List[*count_Block].IM_lead_lag[i].Exogenous = (int*)malloc(tmp_size_exo[incidencematrix.Model_Max_Lag_Exo - Lag + i] * sizeof(int));
              ModelBlock->Block_List[*count_Block].IM_lead_lag[i].Exogenous_Index = (int*)malloc(tmp_size_exo[incidencematrix.Model_Max_Lag_Exo - Lag + i] * sizeof(int));
              ModelBlock->Block_List[*count_Block].IM_lead_lag[i].Equ_X = (int*)malloc(tmp_size_exo[incidencematrix.Model_Max_Lag_Exo - Lag + i] * sizeof(int));
              ModelBlock->Block_List[*count_Block].IM_lead_lag[i].Equ_X_Index = (int*)malloc(tmp_size_exo[incidencematrix.Model_Max_Lag_Exo - Lag + i] * sizeof(int));
            }
          else
            ModelBlock->Block_List[*count_Block].IM_lead_lag[i].size_exo = 0;
          ModelBlock->Block_List[*count_Block].IM_lead_lag[i].u_init = l;
          memset(tmp_variable_evaluated, 0, symbol_table.endo_nbr*sizeof(bool));
          IM = incidencematrix.Get_IM(i - Lag, eEndogenous);
          if(IM)
            {
              for(j = first_count_equ;j < size + first_count_equ;j++)
                {
                  i_1 = Index_Var_IM[j].index;
                  m = 0;
                  for(k = first_count_equ;k < size + first_count_equ;k++)
                    if(IM[i_1 + Index_Equ_IM[k].index*symbol_table.endo_nbr])
                      m++;
                  if(m > 0)
                    {
                      tmp_var[j - first_count_equ] = i1;
                      i1++;
                    }
                }
              m = 0;
              for(j = first_count_equ;j < size + first_count_equ;j++)
                {
                  i_1 = Index_Equ_IM[j].index * symbol_table.endo_nbr;
                  for(k = first_count_equ;k < size + first_count_equ;k++)
                    if(IM[Index_Var_IM[k].index + i_1])
                      {
                        if(i == Lag)
                          {
                            ModelBlock->Block_List[*count_Block].IM_lead_lag[i].us[m] = ls;
                            ls++;
                          }
                        ModelBlock->Block_List[*count_Block].IM_lead_lag[i].u[m] = li;
                        ModelBlock->Block_List[*count_Block].IM_lead_lag[i].Equ[m] = j - first_count_equ;
                        ModelBlock->Block_List[*count_Block].IM_lead_lag[i].Var[m] = k - first_count_equ;
                        ModelBlock->Block_List[*count_Block].IM_lead_lag[i].Equ_Index[m] = Index_Equ_IM[j].index;
                        ModelBlock->Block_List[*count_Block].IM_lead_lag[i].Var_Index[m] = Index_Var_IM[k].index;
                        tmp_variable_evaluated[Index_Var_IM[k].index] = true;
                        l++;
                        m++;
                        li++;
                      }
                }
              ModelBlock->Block_List[*count_Block].IM_lead_lag[i].u_finish = li - 1;
              m = 0;
              for(j = first_count_equ;j < size + first_count_equ;j++)
                {
                  i_1 = Index_Equ_IM[j].index * symbol_table.endo_nbr;
                  for(k = 0;k < symbol_table.endo_nbr;k++)
                    if((!tmp_variable_evaluated[Index_Var_IM[k].index]) && IM[Index_Var_IM[k].index + i_1])
                      {
                        ModelBlock->Block_List[*count_Block].IM_lead_lag[i].u_other_endo[m] = l;
                        ModelBlock->Block_List[*count_Block].IM_lead_lag[i].Equ_other_endo[m] = j - first_count_equ;
                        ModelBlock->Block_List[*count_Block].IM_lead_lag[i].Var_other_endo[m] = k - first_count_equ;
                        ModelBlock->Block_List[*count_Block].IM_lead_lag[i].Equ_Index_other_endo[m] = Index_Equ_IM[j].index;
                        ModelBlock->Block_List[*count_Block].IM_lead_lag[i].Var_Index_other_endo[m] = Index_Var_IM[k].index;
                        l++;
                        m++;
                      }
                }
            }
          IM = incidencematrix.Get_IM(i - Lag, eExogenous);
          if(IM)
            {
              m = 0;
              for(j = first_count_equ;j < size + first_count_equ;j++)
                {
                  i_1 = Index_Equ_IM[j].index * symbol_table.exo_nbr;
                  for(k = 0; k<tmp_nb_exo; k++)
                    {
                      if(IM[ModelBlock->Block_List[*count_Block].Exogenous[k]+i_1])
                        {
                          ModelBlock->Block_List[*count_Block].IM_lead_lag[i].Exogenous[m] = k;
                          ModelBlock->Block_List[*count_Block].IM_lead_lag[i].Exogenous_Index[m] = ModelBlock->Block_List[*count_Block].Exogenous[k];
                          ModelBlock->Block_List[*count_Block].IM_lead_lag[i].Equ_X[m] = j - first_count_equ;
                          ModelBlock->Block_List[*count_Block].IM_lead_lag[i].Equ_X_Index[m] = Index_Equ_IM[j].index;
                          m++;
                        }
                    }
                }
            }
        }
      (*count_Block)++;
      free(tmp_size);
      free(tmp_size_other_endo);
      free(tmp_size_exo);
      free(tmp_endo);
      free(tmp_other_endo);
      free(tmp_exo);
      free(tmp_var);
      free(tmp_variable_evaluated);
    }
}


void
BlockTriangular::Free_Block(Model_Block* ModelBlock) const
{
  int blk, i;
  for(blk = 0;blk < ModelBlock->Size;blk++)
    {

      if ((ModelBlock->Block_List[blk].Type == PROLOGUE) || (ModelBlock->Block_List[blk].Type == EPILOGUE))
        {
          free(ModelBlock->Block_List[blk].Equation);
          free(ModelBlock->Block_List[blk].Variable);
          free(ModelBlock->Block_List[blk].Exogenous);
          free(ModelBlock->Block_List[blk].Own_Derivative);
          for(i = 0;i < ModelBlock->Block_List[blk].Max_Lag + ModelBlock->Block_List[blk].Max_Lead + 1;i++)
            {
              if(ModelBlock->Block_List[blk].IM_lead_lag[i].size)
                {
                  free(ModelBlock->Block_List[blk].IM_lead_lag[i].u);
                  free(ModelBlock->Block_List[blk].IM_lead_lag[i].us);
                  free(ModelBlock->Block_List[blk].IM_lead_lag[i].Var);
                  free(ModelBlock->Block_List[blk].IM_lead_lag[i].Equ);
                  free(ModelBlock->Block_List[blk].IM_lead_lag[i].Var_Index);
                  free(ModelBlock->Block_List[blk].IM_lead_lag[i].Equ_Index);
                }
              if(ModelBlock->Block_List[blk].IM_lead_lag[i].size_exo)
                {
                  free(ModelBlock->Block_List[blk].IM_lead_lag[i].Exogenous);
                  free(ModelBlock->Block_List[blk].IM_lead_lag[i].Exogenous_Index);
                  free(ModelBlock->Block_List[blk].IM_lead_lag[i].Equ_X_Index);
                  free(ModelBlock->Block_List[blk].IM_lead_lag[i].Equ_X);
                }
            }
          free(ModelBlock->Block_List[blk].IM_lead_lag);
          delete(ModelBlock->Block_List[blk].Temporary_terms);
        }
      else
        {
          free(ModelBlock->Block_List[blk].Equation);
          free(ModelBlock->Block_List[blk].Variable);
          free(ModelBlock->Block_List[blk].Exogenous);
          free(ModelBlock->Block_List[blk].Own_Derivative);
          for(i = 0;i < ModelBlock->Block_List[blk].Max_Lag + ModelBlock->Block_List[blk].Max_Lead + 1;i++)
            {
              if(incidencematrix.Model_Max_Lag_Endo-ModelBlock->Block_List[blk].Max_Lag+i>=0 && ModelBlock->Block_List[blk].IM_lead_lag[i].size)
                {
                  free(ModelBlock->Block_List[blk].IM_lead_lag[i].u);
                  free(ModelBlock->Block_List[blk].IM_lead_lag[i].us);
                  free(ModelBlock->Block_List[blk].IM_lead_lag[i].Equ);
                  free(ModelBlock->Block_List[blk].IM_lead_lag[i].Var);
                  free(ModelBlock->Block_List[blk].IM_lead_lag[i].Equ_Index);
                  free(ModelBlock->Block_List[blk].IM_lead_lag[i].Var_Index);
                }
              if(incidencematrix.Model_Max_Lag_Exo-ModelBlock->Block_List[blk].Max_Lag+i>=0 && ModelBlock->Block_List[blk].IM_lead_lag[i].size_exo)
                {
                  free(ModelBlock->Block_List[blk].IM_lead_lag[i].Exogenous);
                  free(ModelBlock->Block_List[blk].IM_lead_lag[i].Exogenous_Index);
                  free(ModelBlock->Block_List[blk].IM_lead_lag[i].Equ_X_Index);
                  free(ModelBlock->Block_List[blk].IM_lead_lag[i].Equ_X);
                }
            }
          free(ModelBlock->Block_List[blk].IM_lead_lag);
          delete(ModelBlock->Block_List[blk].Temporary_terms);
        }
    }
  free(ModelBlock->Block_List);
  free(ModelBlock);
  free(Index_Equ_IM);
  free(Index_Var_IM);
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
  Equation_set* Equation_gr = (Equation_set*) malloc(sizeof(Equation_set));
  bool* SIM0, *SIM00;
  SIM0 = (bool*)malloc(n * n * sizeof(bool));
  memcpy(SIM0,IM_0,n*n*sizeof(bool));
  Prologue_Epilogue(IM, prologue, epilogue, n, Index_Var_IM, Index_Equ_IM, SIM0);
  free(SIM0);
  if(bt_verbose)
    {
      cout << "prologue : " << *prologue << " epilogue : " << *epilogue << "\n";
      cout << "IM_0\n";
      incidencematrix.Print_SIM(IM_0, eEndogenous);
      cout << "IM\n";
      incidencematrix.Print_SIM(IM, eEndogenous);
      for(i = 0;i < n;i++)
        cout << "Index_Var_IM[" << i << "]=" << Index_Var_IM[i].index << " Index_Equ_IM[" << i << "]=" << Index_Equ_IM[i].index << "\n";
    }
  if(*prologue+*epilogue<n)
    {
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
                              cout << "Error nothing at IM_0[" << iter->first.first << ", " << iter->first.second << "]=" << IM_0[iter->first.first*n+iter->first.second] << "  " << iter->second << "\n";
                            }
                        }
                      else
                        suppress++;
                    }
                  for(i = 0;i < n;i++)
                    for(j = 0;j < n;j++)
                      {
                        SIM00[i*n + j] = SIM0[Index_Equ_IM[i].index * n + Index_Var_IM[j].index];
                      }
                  free(SIM0);
                  if(suppress!=suppressed)
                    {
                      OK=normalization.Normalize(n, *prologue, *epilogue, SIM00, Index_Equ_IM, Equation_gr, 1, IM);
                      if(!OK)
                        normalization.Free_Equation(n-*prologue-*epilogue,Equation_gr);
                    }
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
                  exit(EXIT_FAILURE);
                }
            }
          else
            normalization.Normalize(n, *prologue, *epilogue, IM, Index_Equ_IM, Equation_gr, 0, 0);
        }
      else
        normalization.Gr_to_IM_basic(n, *prologue, *epilogue, IM, Equation_gr, false);
    }
  cout << "Finding the optimal block decomposition of the model ...\n";
  if(*prologue+*epilogue<n)
    {
      if(bt_verbose)
        blocks.Print_Equation_gr(Equation_gr);
      res = blocks.sc(Equation_gr);
      normalization.Free_Equation(n-*prologue-*epilogue,Equation_gr);
      if(bt_verbose)
        blocks.block_result_print(res);
    }
  else
    {
      res = (block_result_t*)malloc(sizeof(*res));
      res->n_sets=0;
    }
  free(Equation_gr);
  if ((*prologue) || (*epilogue))
    j = 1;
  else
    j = 0;
  for(i = 0;i < res->n_sets;i++)
    {
      if ((res->sets_f[i] - res->sets_s[i] + 1) > j)
        j = res->sets_f[i] - res->sets_s[i] + 1;
    }
  Nb_RecursBlocks = *prologue + *epilogue;
  Nb_TotalBlocks = res->n_sets + Nb_RecursBlocks;
  cout << Nb_TotalBlocks << " block(s) found:\n";
  cout << "  " << Nb_RecursBlocks << " recursive block(s) and " << res->n_sets << " simultaneous block(s). \n";
  cout << "  the largest simultaneous block has " << j << " equation(s). \n";
  ModelBlock->Size = Nb_TotalBlocks;
  ModelBlock->Periods = periods;
  ModelBlock->Block_List = (Block*)malloc(sizeof(ModelBlock->Block_List[0]) * Nb_TotalBlocks);
  blocks.block_result_to_IM(res, IM, *prologue, symbol_table.endo_nbr, Index_Equ_IM, Index_Var_IM);
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
  if(res->n_sets)
    blocks.block_result_free(res);
  else
    free(res);
  return 0;
}

//------------------------------------------------------------------------------
// normalize each equation of the dynamic model
// and find the optimal block triangular decomposition of the static model
void
BlockTriangular::Normalize_and_BlockDecompose_Static_0_Model(const jacob_map &j_m)
{
  bool* SIM, *SIM_0;
  bool* Cur_IM;
  int i, k, size;
  //First create a static model incidence matrix
  size = symbol_table.endo_nbr * symbol_table.endo_nbr * sizeof(*SIM);
  SIM = (bool*)malloc(size);
  for(i = 0; i< symbol_table.endo_nbr * symbol_table.endo_nbr; i++) SIM[i] = 0;
  for(k = -incidencematrix.Model_Max_Lag_Endo; k<=incidencematrix.Model_Max_Lead_Endo; k++)
    {
      Cur_IM = incidencematrix.Get_IM(k, eEndogenous);
      if(Cur_IM)
        {
          for(i = 0;i < symbol_table.endo_nbr*symbol_table.endo_nbr;i++)
            {
              SIM[i] = (SIM[i]) || (Cur_IM[i]);
            }
        }
    }
  if(bt_verbose)
    {
      cout << "incidence matrix for the static model (unsorted) \n";
      incidencematrix.Print_SIM(SIM, eEndogenous);
    }
  Index_Equ_IM = (simple*)malloc(symbol_table.endo_nbr * sizeof(*Index_Equ_IM));
  for(i = 0;i < symbol_table.endo_nbr;i++)
    {
      Index_Equ_IM[i].index = i;
    }
  Index_Var_IM = (simple*)malloc(symbol_table.endo_nbr * sizeof(*Index_Var_IM));
  for(i = 0;i < symbol_table.endo_nbr;i++)
    {
      Index_Var_IM[i].index = i;
    }
  if(ModelBlock != NULL)
    Free_Block(ModelBlock);
  ModelBlock = (Model_Block*)malloc(sizeof(*ModelBlock));
  Cur_IM = incidencematrix.Get_IM(0, eEndogenous);
  SIM_0 = (bool*)malloc(symbol_table.endo_nbr * symbol_table.endo_nbr * sizeof(*SIM_0));
  for(i = 0;i < symbol_table.endo_nbr*symbol_table.endo_nbr;i++)
    SIM_0[i] = Cur_IM[i];
  Normalize_and_BlockDecompose(SIM, ModelBlock, symbol_table.endo_nbr, &prologue, &epilogue, Index_Var_IM, Index_Equ_IM, 1, 1, SIM_0, j_m);
  free(SIM_0);
  free(SIM);
}
