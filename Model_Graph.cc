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

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <string>
#include <ctime>
#include <stack>
#include <cmath>
#include "ModelTree.hh"
#include "Model_Graph.hh"
#include "BlockTriangular.hh"

using namespace std;

void
free_model_graph(t_model_graph* model_graph)
{
  int i;
  for(i = 0;i < model_graph->nb_vertices;i++)
    {
      free(model_graph->vertex[i].in_degree_edge);
      free(model_graph->vertex[i].out_degree_edge);
    }
  free(model_graph->vertex);
  free(model_graph);
}

void
print_Graph(t_model_graph* model_graph)
{
  int i, j;
  for(i = 0;i < model_graph->nb_vertices;i++)
    {
      cout << "vertex " << model_graph->vertex[i].index << "(" << i << " ," << model_graph->vertex[i].nb_out_degree_edges << ")\n";
      cout << "      -> ";
      for(j = 0;j < model_graph->vertex[i].nb_out_degree_edges;j++)
        cout << model_graph->vertex[model_graph->vertex[i].out_degree_edge[j].index].index <<  /*" -" << model_graph->vertex[i].out_degree_edge[j].index <<  "-*/" (" << model_graph->vertex[i].out_degree_edge[j].u_count << "), ";
      cout << "\n";
      cout << "      <- ";
      for(j = 0;j < model_graph->vertex[i].nb_in_degree_edges;j++)
        cout << model_graph->vertex[model_graph->vertex[i].in_degree_edge[j].index].index <<  /*" -" << model_graph->vertex[i].in_degree_edge[j].index <<  "-*/" (" << model_graph->vertex[i].in_degree_edge[j].u_count << "), ";
      cout << "\n";
    }
}



void Check_Graph(t_model_graph* model_graph)
{
  int i, j, k, i1, i2;
  bool OK, OK_u_count;
  for(i = 0;i < model_graph->nb_vertices;i++)
    {
      for(j = 0;j < model_graph->vertex[i].nb_in_degree_edges;j++)
        {
          i1 = model_graph->vertex[i].in_degree_edge[j].index;
          i2 = model_graph->vertex[i].in_degree_edge[j].u_count;
          OK = 0;
          OK_u_count = 0;
          for(k = 0;(k < model_graph->vertex[i1].nb_out_degree_edges) && (!OK);k++)
            {
              if(model_graph->vertex[i1].out_degree_edge[k].index == i)
                {
                  OK = 1;
                  if(model_graph->vertex[i1].out_degree_edge[k].u_count == i2)
                    OK_u_count = 1;
                }
            }
          if(!OK)
            {
              cout << "not symetric for edge between vertices " << model_graph->vertex[i1].index << " and " << model_graph->vertex[i].index << " (in_degree)\n";
              print_Graph(model_graph);
              system("pause");
              exit( -1);
            }
          if(!OK_u_count)
            {
              cout << "valeur de u_count non symétrique sur l'arc entre " << model_graph->vertex[i1].index << " et " << model_graph->vertex[i].index << " (in_degree)\n";
              print_Graph(model_graph);
              system("pause");
              exit( -1);
            }
        }
      for(j = 0;j < model_graph->vertex[i].nb_out_degree_edges;j++)
        {
          i1 = model_graph->vertex[i].out_degree_edge[j].index;
          i2 = model_graph->vertex[i].out_degree_edge[j].u_count;
          OK = 0;
          OK_u_count = 0;
          for(k = 0;(k < model_graph->vertex[i1].nb_in_degree_edges) && (!OK);k++)
            {
              if(model_graph->vertex[i1].in_degree_edge[k].index == i)
                {
                  OK = 1;
                  if(model_graph->vertex[i1].in_degree_edge[k].u_count == i2)
                    OK_u_count = 1;
                }
            }
          if(!OK)
            {
              cout << "pas symétrique sur l'arc entre " << model_graph->vertex[i1].index << " et " << model_graph->vertex[i].index << " (out_degree)\n";
              print_Graph(model_graph);
              system("pause");
              exit( -1);
            }
          if(!OK_u_count)
            {
              cout << "valeur de u_count non symétrique sur l'arc entre " << model_graph->vertex[i1].index << " et " << model_graph->vertex[i].index << " (out_degree)\n";
              print_Graph(model_graph);
              system("pause");
              exit( -1);
            }
        }
    }
}

int
ModelBlock_Graph(Model_Block *ModelBlock, int Blck_num, bool dynamic, t_model_graph* model_graph, int nb_endo, int* block_u_count, int *starting_vertex, int *periods, int *nb_table_y, int *mean_var_in_equ)
{
  int i, j, k, l, m, lag, per, lag1, k2, complete_size = 0, u_count;
  int max_lead, max_lag, size, Lead, Lag;
  int *Used, *todo_lag, *todo_lead, *vertex_ref, *vertex_index, *todo_lag1, *todo_lead1 ;
  max_lag = ModelBlock->Block_List[Blck_num].Max_Lag;
  max_lead = ModelBlock->Block_List[Blck_num].Max_Lead;
  if(!dynamic)
    {
      /*It's a static model that have to be solved at each period*/
      /*size=ModelBlock->Block_List[Blck_num].IM_lead_lag[max_lag].size;*/
      size = ModelBlock->Block_List[Blck_num].Size;
      /*We add an extra vertex to take into account of the f(x0) constant term in f(x)=0 approximated by f(x0) + (x-x0) f'(x0) = 0*/
      //cout << "Static, Blck_num= " << Blck_num << "size= " << size << "\n";
      model_graph->nb_vertices = size + 1;
      *starting_vertex = 0;
      model_graph->vertex = (t_vertex*)malloc(model_graph->nb_vertices * sizeof(*model_graph->vertex));
      for(i = 0;i < size;i++)
        {
          /*It's not f(x0) vertex*/
          model_graph->vertex[i].in_degree_edge = (t_edge*)malloc((size + 1) * sizeof(t_edge));
          model_graph->vertex[i].out_degree_edge = (t_edge*)malloc((size + 1) * sizeof(t_edge));
          model_graph->vertex[i].nb_in_degree_edges = 0;
          model_graph->vertex[i].nb_out_degree_edges = 0;
          model_graph->vertex[i].index = ModelBlock->Block_List[Blck_num].Variable[i];
          model_graph->vertex[i].lag_lead = 0;
        }
      /*It's f(x0) vertex*/
      model_graph->vertex[size].in_degree_edge = (t_edge*)malloc(0 * sizeof(t_edge));
      model_graph->vertex[size].out_degree_edge = (t_edge*)malloc((size) * sizeof(t_edge));
      model_graph->vertex[size].nb_in_degree_edges = 0;
      model_graph->vertex[size].index = -1;
      model_graph->vertex[size].lag_lead = 0;
      for(i = 0;i < ModelBlock->Block_List[Blck_num].IM_lead_lag[max_lag].size;i++)
        {
          k = ModelBlock->Block_List[Blck_num].IM_lead_lag[max_lag].Equ[i];
          m = ModelBlock->Block_List[Blck_num].IM_lead_lag[max_lag].Var[i];
          j = model_graph->vertex[k].nb_in_degree_edges++;
          l = model_graph->vertex[m].nb_out_degree_edges++;
          model_graph->vertex[k].in_degree_edge[j].index = m;
          model_graph->vertex[m].out_degree_edge[l].index = k;
          model_graph->vertex[k].in_degree_edge[j].u_count = ModelBlock->Block_List[Blck_num].IM_lead_lag[max_lag].us[i];
          model_graph->vertex[m].out_degree_edge[l].u_count = ModelBlock->Block_List[Blck_num].IM_lead_lag[max_lag].us[i];
        }
      model_graph->vertex[size].nb_out_degree_edges = size;
      for(i = 0;i < size;i++)
        {
          j = model_graph->vertex[i].nb_in_degree_edges++;
          model_graph->vertex[i].in_degree_edge[j].index = size;
          model_graph->vertex[i].in_degree_edge[j].u_count = i;
          model_graph->vertex[size].out_degree_edge[i].index = i;
          model_graph->vertex[size].out_degree_edge[i].u_count = i;
        }
      u_count = ModelBlock->Block_List[Blck_num].IM_lead_lag[max_lag].u_finish - ModelBlock->Block_List[Blck_num].IM_lead_lag[max_lag].u_init + 1
        + ModelBlock->Block_List[Blck_num].Size;
      *block_u_count = u_count;
      *nb_table_y = size;
      return (u_count);
    }
  else
    {
      int sup;
      Lead = ModelBlock->Block_List[Blck_num].Max_Lead;
      Lag = ModelBlock->Block_List[Blck_num].Max_Lag;
      cout << "---> *periods=" << *periods << "\n";
      if(*periods>3)
        {
          sup = Lead + Lag +3;
          *periods = Lead + Lag  + sup;
        }
#ifdef PRINT_OUT
      cout << "Lag=" << Lag << " Lead=" << Lead << "\n";
      cout << "periods=Lead+2*Lag+2= " << *periods << "\n";
#endif
      size = ModelBlock->Block_List[Blck_num].Size;
      /*It's a dynamic model that have to be solved for all periods.
        So we consider the incidence matrice for all lead and lags plus the current value*/
      model_graph->nb_vertices = 0;
      vertex_ref = (int*)malloc(size * (Lag + Lead + *periods) * sizeof(int));
      memset(vertex_ref, -1, size*(Lag + Lead + *periods)*sizeof(int));
      vertex_index = (int*)malloc(size * (Lag + Lead + *periods) * sizeof(int));
      complete_size = ModelBlock->Block_List[Blck_num].IM_lead_lag[Lag].size * (*periods);
      if(Lag > 0)
        {
          todo_lag = (int*)malloc(size * Lag * sizeof(int));
          todo_lag1 = (int*)malloc(size * Lag * sizeof(int));
          memset(todo_lag, -1, size*Lag*sizeof(int));
          memset(todo_lag1, -1, size*Lag*sizeof(int));
          Used = (int*)malloc(size * Lag * sizeof(int));
          for(lag = 0;lag < Lag;lag++)
            {
              memset(Used, -1, size*Lag*sizeof(int));
              complete_size += ModelBlock->Block_List[Blck_num].IM_lead_lag[lag].size;
              for(i = 0;i < ModelBlock->Block_List[Blck_num].IM_lead_lag[lag].size;i++)
                {
                  if(Used[ModelBlock->Block_List[Blck_num].IM_lead_lag[lag].Var[i]] < 0)
                    {
                      k = ModelBlock->Block_List[Blck_num].IM_lead_lag[lag].Var[i];
                      todo_lag[lag*size + k] = k;
                      vertex_ref[lag*size + k] = model_graph->nb_vertices;
                      vertex_index[model_graph->nb_vertices] = lag * nb_endo + ModelBlock->Block_List[Blck_num].IM_lead_lag[lag].Var_Index[i];
                      todo_lag1[lag*size + k] = i;
                      model_graph->nb_vertices++;
                      Used[k] = i;
                    }
                }
              if(lag > 0)
                {
                  for(lag1 = 0;lag1 < lag;lag1++)
                    for(i = 0;i < size;i++)
                      if(todo_lag[(lag1)*size + i] >= 0)
                        {
                          if(Used[i] < 0)
                            {
                              todo_lag[lag*size + i] = i;
                              k = todo_lag[(lag1) * size + i];
                              vertex_ref[lag*size + k] = model_graph->nb_vertices;
                              j = todo_lag1[(lag1) * size + i];
                              vertex_index[model_graph->nb_vertices] = lag * nb_endo + ModelBlock->Block_List[Blck_num].IM_lead_lag[lag1].Var_Index[k];
                              model_graph->nb_vertices++;
                            }
                        }
                }
            }
          *starting_vertex = model_graph->nb_vertices;
          free(Used);
          free(todo_lag);
          free(todo_lag1);
        }
      int nb_vertices_1=model_graph->nb_vertices;
#ifdef PRINT_OUT
      cout << "nb_vertices in the first part: " << nb_vertices_1 << "\n";
#endif
      for(per = Lag;per < Lag + *periods;per++)
        for(i = 0;i < size;i++)
          {
            vertex_ref[per*size + i] = model_graph->nb_vertices;
            vertex_index[model_graph->nb_vertices] = (per) * nb_endo + ModelBlock->Block_List[Blck_num].Variable[i];
            model_graph->nb_vertices++;
          }
      int nb_vertices_2=model_graph->nb_vertices-nb_vertices_1;
#ifdef PRINT_OUT
      cout << "nb_vertices in the second part: " << nb_vertices_2 << "\n";
#endif
      if(Lead > 0)
        {
          todo_lead = (int*)malloc(size * Lead * sizeof(int));
          todo_lead1 = (int*)malloc(size * Lead * sizeof(int));
          memset(todo_lead, -1, size*Lead*sizeof(int));
          memset(todo_lead1, -1, size*Lead*sizeof(int));
          Used = (int*)malloc(size * Lead * sizeof(int));
          k2 = model_graph->nb_vertices;
          for(lag = Lag + Lead;lag > Lag;lag--)
            {
              complete_size += ModelBlock->Block_List[Blck_num].IM_lead_lag[lag].size;
              memset(Used, -1, size*Lead*sizeof(int));
              for(i = 0;i < ModelBlock->Block_List[Blck_num].IM_lead_lag[lag].size;i++)
                {
                  if(Used[ModelBlock->Block_List[Blck_num].IM_lead_lag[lag].Var[i]] < 0)
                    {
                      k = ModelBlock->Block_List[Blck_num].IM_lead_lag[lag].Var[i];
                      todo_lead[(lag - Lag - 1)*size + k] = k;
                      todo_lead1[(lag - Lag - 1)*size + k] = i;
                      Used[k] = i;
                      model_graph->nb_vertices++;
                    }
                }
              if(lag < Lag + Lead)
                {
                  for(lag1 = Lag + Lead;lag1 > lag;lag1--)
                    for(i = 0;i < size;i++)
                      {
                        if(todo_lead[(lag1 - Lag - 1)*size + i] >= 0)
                          {
                            if(Used[i] < 0)
                              {
                                k = todo_lead[(lag1 - Lag - 1) * size + i];
                                model_graph->nb_vertices++;
                              }
                          }
                      }
                }
            }
          k2 = model_graph->nb_vertices;
          memset(todo_lead, -1, size*Lead*sizeof(int));
          for(lag = Lag + Lead;lag > Lag;lag--)
            {
              complete_size += ModelBlock->Block_List[Blck_num].IM_lead_lag[lag].size;
              memset(Used, -1, size*Lead*sizeof(int));
              for(i = ModelBlock->Block_List[Blck_num].IM_lead_lag[lag].size - 1;i >= 0;i--)
                {
                  if(Used[ModelBlock->Block_List[Blck_num].IM_lead_lag[lag].Var[i]] < 0)
                    {
                      k2--;
                      k = ModelBlock->Block_List[Blck_num].IM_lead_lag[lag].Var[i];
                      todo_lead[(lag - Lag - 1)*size + k] = k;
                      todo_lead1[(lag - Lag - 1)*size + k] = i;
                      vertex_ref[(lag + *periods - 1)*size + k] = k2;
                      vertex_index[k2] = (lag + *periods - 1) * nb_endo + ModelBlock->Block_List[Blck_num].IM_lead_lag[lag].Var_Index[i];
                      Used[k] = i;
                    }
                }
              if(lag < Lag + Lead)
                {
                  for(lag1 = Lag + Lead;lag1 > lag;lag1--)
                    {
                      for(i = size - 1;i >= 0;i--)
                        {
                          if(todo_lead[(lag1 - Lag - 1)*size + i] >= 0)
                            {
                              if(Used[i] < 0)
                                {
                                  k2--;
                                  todo_lead[(lag - Lag - 1)*size + i] = i;
                                  todo_lead1[(lag - Lag - 1)*size + i] = todo_lead1[(lag1 - Lag - 1)*size + i];
                                  k = todo_lead[(lag1 - Lag - 1) * size + i];
                                  vertex_ref[(lag + *periods - 1)*size + k] = k2;
                                  //#ifdef PRINT_OUT

                                  //#endif

                                  j = todo_lead1[(lag1 - Lag - 1) * size + i];
                                  //#ifdef PRINT_OUT
                                  if(j>ModelBlock->Block_List[Blck_num].IM_lead_lag[lag1].size||j==-1)
                                    {
                                      cout << "Error in model graph construction (lead part): j (" << j << ")>size (" << ModelBlock->Block_List[Blck_num].IM_lead_lag[lag1].size << ")\n";
                                      system("pause");
                                      exit(-1);
                                    }
                                  //#endif
                                  vertex_index[k2] = (lag + *periods - 1) * nb_endo + ModelBlock->Block_List[Blck_num].IM_lead_lag[lag1].Var_Index[j];
                                }
                            }
                        }
                    }
                }
            }
          free(Used);
          free(todo_lead);
          free(todo_lead1);
        }
      int nb_vertices_3=model_graph->nb_vertices-nb_vertices_2-nb_vertices_1;
#ifdef PRINT_OUT
      cout << "nb_vertices in the last part: " << nb_vertices_3 << "\n";
#endif
      /*We add an extra vertex to take into account of the f(x0) constant term in f(x)=0 approx f(x0) + (x-x0) f'(x0) = 0*/
      model_graph->nb_vertices++;
      model_graph->vertex = (t_vertex*)malloc(model_graph->nb_vertices * sizeof(*model_graph->vertex));
      vertex_index[model_graph->nb_vertices - 1] = -1;
#ifdef PRINT_OUT
      cout << "ok0\n";
      cout << "model_graph->nb_vertices=" << model_graph->nb_vertices << " Lag=" << Lag << " Lead=" << Lead << "\n";
      cout << "*periods=" << *periods << " size=" << size << "\n";
      cout << "allocated / vertex = " << (size + nb_vertices_1 + nb_vertices_3+ 1) << "\n";
#endif
      int nb_table_u= size+nb_vertices_1+nb_vertices_3+2;
      for(k = 0;k < model_graph->nb_vertices-1;k++)
        {
          model_graph->vertex[k].index = vertex_index[k];
          model_graph->vertex[k].in_degree_edge = (t_edge*)malloc(nb_table_u * sizeof(t_edge));
          model_graph->vertex[k].out_degree_edge = (t_edge*)malloc(nb_table_u * sizeof(t_edge));
          model_graph->vertex[k].nb_in_degree_edges = 0;
          model_graph->vertex[k].nb_out_degree_edges = 0;
          model_graph->vertex[k].max_nb_in_degree_edges = nb_table_u;
          model_graph->vertex[k].max_nb_out_degree_edges = nb_table_u;
#ifdef PRINT_OUT
          //if(k==8)
          {
            cout << " model_graph->vertex[" << k << "].in_degree_edge=" << model_graph->vertex[k].in_degree_edge << "\n";
          }
#endif
        }
      model_graph->vertex[model_graph->nb_vertices-1].index = vertex_index[model_graph->nb_vertices-1];
      model_graph->vertex[model_graph->nb_vertices-1].in_degree_edge = (t_edge*)malloc(/*model_graph->nb_vertices **/ sizeof(t_edge));
      model_graph->vertex[model_graph->nb_vertices-1].out_degree_edge = (t_edge*)malloc(model_graph->nb_vertices * sizeof(t_edge));
      model_graph->vertex[model_graph->nb_vertices-1].nb_in_degree_edges = 0;
      model_graph->vertex[model_graph->nb_vertices-1].nb_out_degree_edges = 0;
      model_graph->vertex[model_graph->nb_vertices-1].max_nb_in_degree_edges = 0;
      model_graph->vertex[model_graph->nb_vertices-1].max_nb_out_degree_edges = model_graph->nb_vertices;
#ifdef PRINT_OUT
      cout << "ok1\n";
      system("pause");
#endif
      u_count = 0;
      *mean_var_in_equ = 0;
      for(per = 0;per < *periods;per++)
        {
          j = model_graph->nb_vertices - 1;
          for(i = 0;i < size;i++)
            {
              k = vertex_ref[(Lag + per) * size + i];
              model_graph->vertex[k].in_degree_edge[model_graph->vertex[k].nb_in_degree_edges].index = j;
              model_graph->vertex[j].out_degree_edge[model_graph->vertex[j].nb_out_degree_edges].index = k;
              model_graph->vertex[k].in_degree_edge[model_graph->vertex[k].nb_in_degree_edges].u_count = u_count;
              model_graph->vertex[j].out_degree_edge[model_graph->vertex[j].nb_out_degree_edges].u_count = u_count;
              model_graph->vertex[k].nb_in_degree_edges++;
              model_graph->vertex[j].nb_out_degree_edges++;
              u_count++;
            }
          for(lag = 0;lag < Lag + Lead + 1;lag++)
            {
#ifdef PRINT_OUT
              cout << "ModelBlock->Block_List[" << Blck_num << "].IM_lead_lag[" << lag << "].size = " << ModelBlock->Block_List[Blck_num].IM_lead_lag[lag].size << "\n";
#endif

              for(i = 0;i < ModelBlock->Block_List[Blck_num].IM_lead_lag[lag].size;i++)
                {
                  j = vertex_ref[(lag + per) * size + ModelBlock->Block_List[Blck_num].IM_lead_lag[lag].Var[i]];
                  k = vertex_ref[(Lag + per) * size + ModelBlock->Block_List[Blck_num].IM_lead_lag[lag].Equ[i]];
#ifdef PRINT_OUT

                  cout << "per=" << per << " lag=" << lag << " i=" << i << " j=" << j << " k=" << k << "\n";
#endif

                  model_graph->vertex[k].in_degree_edge[model_graph->vertex[k].nb_in_degree_edges].index = j;
                  model_graph->vertex[j].out_degree_edge[model_graph->vertex[j].nb_out_degree_edges].index = k;
                  model_graph->vertex[k].in_degree_edge[model_graph->vertex[k].nb_in_degree_edges].u_count = u_count;
                  model_graph->vertex[j].out_degree_edge[model_graph->vertex[j].nb_out_degree_edges].u_count = u_count;
                  if(per==(Lag+2))/*&&(lag==(Lag+1))*/
                    (*mean_var_in_equ)++;
                  model_graph->vertex[k].nb_in_degree_edges++;
                  model_graph->vertex[j].nb_out_degree_edges++;
                  u_count++;
                }
            }
        }
      (*mean_var_in_equ) += size;
      //cout << "Total variables=" << *mean_var_in_equ << " nb_endo=" << size << "\n";
      i=*mean_var_in_equ ;
      i =int(ceil(double(i)/size));
      *mean_var_in_equ = i;

      //cout << "Mean var in equation=" << *mean_var_in_equ  << "\n";
      *block_u_count = u_count / (*periods);
      free(vertex_index);
      free(vertex_ref);
      if(nb_vertices_1+nb_vertices_3+1>size)
        *nb_table_y = nb_vertices_1+nb_vertices_3+1;
      else
        *nb_table_y = size;
      return (u_count);
    }
}
