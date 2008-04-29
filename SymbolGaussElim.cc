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
#include <iostream>
#include <string>
#include <ctime>
#include <stack>
#include <cmath>
#include <fstream>

#include "ModelTree.hh"
//#include "pctimer_h.hh"
#include "Model_Graph.hh"
#include "SymbolGaussElim.hh"
//#define SIMPLIFY

using namespace std;
int max_nb_table_y, max_nb_in_degree_edges=0, max_nb_out_degree_edges=0;

SymbolicGaussElimination::SymbolicGaussElimination()
{
#ifdef SAVE
  file_open = false;
#endif
}

#ifdef SIMPLIFY

void
SymbolicGaussElimination::set_free_u_list(int index)
{
#ifndef UNSORTED
  int i = 0, j;
  if(nb_free_u_list > 0)
    {
      while((free_u_list[i] > index) && (i < nb_free_u_list))
        i++;
      for(j = nb_free_u_list;j > i;j--)
        {
          free_u_list[j] = free_u_list[j - 1];
          if(j >= MAX_FREE_U_LIST)
            {
              cout << "Error to much j=" << j << " Maximum allowed=" << MAX_FREE_U_LIST << "\n";
              system("pause");
              exit( -1);
            }
        }
    }
  if(i >= MAX_FREE_U_LIST)
    {
      cout << "Error to much i=" << i << " Maximum allowed=" << MAX_FREE_U_LIST << "\n";
      system("pause");
      exit( -1);
    }
  free_u_list[i] = index;
  nb_free_u_list++;
  if(nb_free_u_list >= MAX_FREE_U_LIST)
    {
      cout << "Error to much free_u_list=" << nb_free_u_list << " Maximum allowed=" << MAX_FREE_U_LIST << "\n";
      system("pause");
      exit( -1);
    }
  if(nb_free_u_list > max_nb_free_u_list)
    max_nb_free_u_list = nb_free_u_list;
#else
  free_u_list[nb_free_u_list++] = index;
#endif
}

int
SymbolicGaussElimination::get_free_u_list(bool dynamic)
{
  int i;
  if((nb_free_u_list > 0) && (simplification_allowed))
    {
      i = free_u_list[--nb_free_u_list];
      return (i);
    }
  else
    if(dynamic)
      return (u_count++);
    else
      return (u_count++);
}
#endif /**SIMPLIFY**/

//==================================================================================
void
SymbolicGaussElimination::init_glb()
{
  tol = TOL;
  vertex_count = 0;
  nb_table_u = 0;
  nb_prologue_save_table_y = 0;
  nb_first_save_table_y = 0;
  nb_middle_save_table_y = 0;
  nb_last_save_table_y = 0;
  nb_prologue_save_table_u = 0;
  nb_first_save_table_u = 0;
  nb_middle_save_table_u = 0;
  nb_last_save_table_u = 0;
  middle_count_loop = 0;
  u_count = 0;
  y_kmin = y_kmax = 0;
  prologue_save_table_u = NULL;
  first_save_table_u = NULL;
  first_save_i_table_u = NULL;
  middle_save_table_u = NULL;
  middle_save_i_table_u = NULL;
  last_save_table_u = NULL;
  prologue_save_table_y = NULL;
  first_save_table_y = NULL;
  first_save_i_table_y = NULL;
  middle_save_table_y = NULL;
  middle_save_i_table_y = NULL;
  last_save_table_y = NULL;
  pos_nb_first_save_table_u = 0;
#ifdef SIMPLIFY
  nb_free_u_list = 0;
  max_nb_free_u_list1 = 0;
  max_nb_free_u_list = 0;
#endif
  save_direct=true;
}


void
SymbolicGaussElimination::write_to_file_table_y( t_table_y *save_table_y, t_table_y *save_i_table_y, int nb_save_table_y, int nb_save_i_table_y)
{
  int i, k;
#ifdef PRINT_OUT
  cout << "nb_save_table_y=" << nb_save_table_y << "\n";
#endif
  SaveCode.write(reinterpret_cast<char *> (&nb_save_table_y), sizeof(nb_save_table_y));
  for(i = 0;i < nb_save_table_y;i++)
    {
      SaveCode.write(reinterpret_cast<char *>(&save_table_y[i].index), sizeof(save_table_y[i].index));
#ifdef PRINT_OUT
      cout << "+y[" << save_table_y[i].index << "]=";
#endif
      SaveCode.write(reinterpret_cast<char *>(&save_table_y[i].nb), sizeof(save_table_y[i].nb));
      for(k = 0;k < save_table_y[i].nb;k++)
        {
          SaveCode.write(reinterpret_cast<char *>(&save_table_y[i].u_index[k]), sizeof(save_table_y[i].u_index[k]));
          SaveCode.write(reinterpret_cast<char *>(&save_table_y[i].y_index[k]), sizeof(save_table_y[i].y_index[k]));
#ifdef PRINT_OUT
          cout << "u[" << save_table_y[i].u_index[k] << "]*y[" << save_table_y[i].y_index[k] << "]";
          if(k + 1 < save_table_y[i].nb)
            cout << "+";
          else
            cout << "\n";
#endif
        }
    }
  for(i = 0;i < nb_save_i_table_y;i++)
    {
      SaveCode.write(reinterpret_cast<char *>(&save_i_table_y[i].index), sizeof(save_i_table_y[i].index));
#ifdef PRINT_OUT
      cout << "+i y[" << save_i_table_y[i].index << "]=";
#endif
      SaveCode.write(reinterpret_cast<char *>(&save_i_table_y[i].nb), sizeof(save_i_table_y[i].nb));
      for(k = 0;k < save_i_table_y[i].nb;k++)
        {
          SaveCode.write(reinterpret_cast<char *>(&save_i_table_y[i].u_index[k]), sizeof(save_i_table_y[i].u_index[k]));
          SaveCode.write(reinterpret_cast<char *>(&save_i_table_y[i].y_index[k]), sizeof(save_i_table_y[i].y_index[k]));
#ifdef PRINT_OUT
          cout << "u[" << save_i_table_y[i].u_index[k] << "]*y[" << save_i_table_y[i].y_index[k] << "]";
          if(k + 1 < save_i_table_y[i].nb)
            cout << "+";
          else
            cout << "\n";
#endif
        }
    }
}

void
SymbolicGaussElimination::write_to_file_table_u_b(t_table_u *save_table_u, t_table_u *last_table_u, int *nb_save_table_u, bool chk)
{
  t_table_u *table_u;
  bool OK=true;
  while(/*save_table_u!=last_table_u*/OK)
    {
#ifdef PRINT_OUT
      cout << "**save_table_u->type=" << int(save_table_u->type) << "\n";
#endif
      switch (save_table_u->type)
        {
        case 3:
        case 7:
          (*nb_save_table_u)++;
          //SaveCode << save_table_u->type << save_table_u->index << save_table_u->op1;
          SaveCode.write(reinterpret_cast<char *>(&save_table_u->type), sizeof(save_table_u->type));
          SaveCode.write(reinterpret_cast<char *>(&save_table_u->index), sizeof(save_table_u->index));
          SaveCode.write(reinterpret_cast<char *>(&save_table_u->op1), sizeof(save_table_u->op1));
#ifdef PRINT_OUT
          if(save_table_u->type == 3)
            cout << "+u[" << save_table_u->index << "]=1/(1-u[" << save_table_u->op1 << "])\n";
          else
            cout << "+u[" << save_table_u->index << "]*=u[" << save_table_u->op1 << "]\n";
#endif /**PRINT_OUT**/
          break;
        case 1:
        case 2:
        case 6:
          (*nb_save_table_u)++;
          //SaveCode << save_table_u->type << save_table_u->index << save_table_u->op1 << save_table_u->op2;
          SaveCode.write(reinterpret_cast<char *>(&save_table_u->type), sizeof(save_table_u->type));
          SaveCode.write(reinterpret_cast<char *>(&save_table_u->index), sizeof(save_table_u->index));
          SaveCode.write(reinterpret_cast<char *>(&save_table_u->op1), sizeof(save_table_u->op1));
          SaveCode.write(reinterpret_cast<char *>(&save_table_u->op2), sizeof(save_table_u->op2));
#ifdef PRINT_OUT
          if(save_table_u->type == 1)
            cout << "+u[" << save_table_u->index << "]=" << "u[" << save_table_u->op1 << "]*u[" << save_table_u->op2 << "]\n";
          else if(save_table_u->type == 2)
            cout << "+u[" << save_table_u->index << "]+=u[" << save_table_u->op1 << "]*u[" << save_table_u->op2 << "]\n";
          else
            cout << "+u[" << save_table_u->index << "]=1/(1-u[" << save_table_u->op1 << "]*u[" << save_table_u->op2 << "])\n";
#endif /**PRINT_OUT**/
          break;
        case 5:
          (*nb_save_table_u)++;
          //SaveCode << save_table_u->type << save_table_u->index;
          SaveCode.write(reinterpret_cast<char *>(&save_table_u->type), sizeof(save_table_u->type));
          SaveCode.write(reinterpret_cast<char *>(&save_table_u->index), sizeof(save_table_u->index));
#ifdef PRINT_OUT
          cout << "+push(u[" << save_table_u->index << "])\n";
#endif /**PRINT_OUT**/
          break;
        }
      if(chk)
        {
          OK=(save_table_u!=last_table_u);
          if(OK)
            {
              table_u = save_table_u->pNext;
              free(save_table_u);
              save_table_u = table_u;
            }
        }
      else
        {
          table_u = save_table_u->pNext;
          free(save_table_u);
          save_table_u = table_u;
          OK=(save_table_u!=last_table_u);
        }
    }
}


void
store_code(t_table_u **save_table_u, t_table_u *First_table_u, t_table_u *stop_table_u, t_table_y *table_y, t_table_y **s_table_y, int *nb_save_table_u, int *nb_save_table_y, int vertex_count)
{
  int i, k;
  t_table_u *s_table_u, *table_u;
  table_u = First_table_u->pNext;
  s_table_u = (t_table_u*)malloc(sizeof(t_table_u) - 3 * sizeof(int));
  (*save_table_u) = s_table_u;
  (*nb_save_table_u) = 0;
  while(table_u != stop_table_u)
    {
      if((table_u->type > 7) || (table_u->type < 1))
        {
          cout << "Error : table_u->type=" << int(table_u->type) << " *nb_save_table_u=" << *nb_save_table_u << "\n";
          exit( -1);
        }
      else
        {
          (*nb_save_table_u)++;
#ifdef PRINT_OUT
          cout << "--table_u->type=" << int(table_u->type) << "\n";
#endif
          switch (table_u->type)
            {
            case 3:
            case 7:
              s_table_u->pNext = (t_table_u*)malloc(sizeof(t_table_u) - sizeof(int));
              s_table_u = s_table_u->pNext;
              s_table_u->type = table_u->type;
              s_table_u->index = table_u->index;
              s_table_u->op1 = table_u->op1;
#ifdef PRINT_OUT
              if(s_table_u->type == 3)
                cout << "st u[" << s_table_u->index << "]=-1/u[" << s_table_u->op1 << "]\n";
              else
                cout << "st u[" << s_table_u->index << "]*=u[" << s_table_u->op1 << "]\n";
#endif
              break;
            case 1:
            case 2:
            case 6:
              s_table_u->pNext = (t_table_u*)malloc(sizeof(t_table_u) );
              s_table_u = s_table_u->pNext;
              s_table_u->type = table_u->type;
              s_table_u->index = table_u->index;
              s_table_u->op1 = table_u->op1;
              s_table_u->op2 = table_u->op2;
#ifdef PRINT_OUT
              if(s_table_u->type == 1)
                cout << "st u[" << s_table_u->index << "]=" << "u[" << s_table_u->op1 << "]*u[" << s_table_u->op2 << "]\n";
              else if(s_table_u->type == 2)
                cout << "st u[" << s_table_u->index << "]+=u[" << s_table_u->op1 << "]*u[" << s_table_u->op2 << "]\n";
              else
                cout << "st u[" << s_table_u->index << "]=1/(1-u[" << s_table_u->op1 << "]*u[" << s_table_u->op2 << "])\n";
#endif /**PRINT_OUT**/
              break;
            case 5:
              s_table_u->pNext = (t_table_u*)malloc(sizeof(t_table_u) - 2 * sizeof(int));
              s_table_u = s_table_u->pNext;
              s_table_u->type = table_u->type;
              s_table_u->index = table_u->index;
#ifdef PRINT_OUT
              cout << "st push(u[" << s_table_u->index << "])\n";
#endif /**PRINT_OUT**/
              break;
            }
        }
      table_u = table_u->pNext;
    }
  s_table_u->pNext = NULL;
#ifdef PRINT_OUT
  cout << "*nb_save_table_u=" << *nb_save_table_u << "\n";
  cout << "vertex_count=" << vertex_count << " nb_save_table_y=" << nb_save_table_y << "\n";
#endif
  (*s_table_y) = (t_table_y*)malloc((vertex_count - *nb_save_table_y) * sizeof(t_table_y));
  for(i = 0;i < vertex_count - *nb_save_table_y;i++)
    {
      (*s_table_y)[i].u_index = (int*)malloc(table_y[i + *nb_save_table_y].nb * sizeof(int));
      (*s_table_y)[i].y_index = (int*)malloc(table_y[i + *nb_save_table_y].nb * sizeof(int));
      /*(*s_table_y)[i].y_lead_lag = (int*)malloc(table_y[i + *nb_save_table_y].nb * sizeof(int));*/
      (*s_table_y)[i].index = table_y[i + *nb_save_table_y].index;
      (*s_table_y)[i].nb = table_y[i + *nb_save_table_y].nb;
      for(k = 0;k < table_y[i + *nb_save_table_y].nb;k++)
        {
          (*s_table_y)[i].u_index[k] = table_y[i + *nb_save_table_y].u_index[k];
          (*s_table_y)[i].y_index[k] = table_y[i + *nb_save_table_y].y_index[k];
        }
    }
  *nb_save_table_y = vertex_count - *nb_save_table_y;
}



bool
SymbolicGaussElimination::Check_Regularity(t_table_u *first_u_blck, t_table_u *second_u_blck, t_table_u *third_u_blck)
{
  t_table_u *c_first_table_u, *c_second_table_u, *c_third_table_u;
  bool OK, Regular;
  int i=0;
  c_first_table_u = first_u_blck->pNext;
  c_second_table_u = second_u_blck->pNext;
  c_third_table_u = third_u_blck->pNext;
  OK = true;
  Regular = true;
  while(OK)
    {
      i++;
#ifdef PRINT_OUT
      cout << "c_first_table_u=" << c_first_table_u << " c_second_table_u=" << c_second_table_u << " c_third_table_u=" << c_third_table_u << "\n";
#endif /**PRINT_OUT**/
      if((c_first_table_u->type != c_second_table_u->type) ||
         (c_second_table_u->index + (c_second_table_u->index - c_first_table_u->index) != c_third_table_u->index))
        {
          Regular = false;
#ifdef PRINT_OUT
          cout << "c_first_table_u->type=" << (int)c_first_table_u->type << "?=c_second_table_u->type=" << (int)c_second_table_u->type << "\n";
          cout << "c_second_table_u->index+(c_second_table_u->index-c_first_table_u->index)=" << c_second_table_u->index + (c_second_table_u->index - c_first_table_u->index)
               << "?=c_third_table_u->index=" << c_third_table_u->index << "\n";
#endif
          break;
        }
      if(c_first_table_u->type != 5)
        {
          if(c_second_table_u->op1 + (c_second_table_u->op1 - c_first_table_u->op1) != c_third_table_u->op1)
            {
#ifdef PRINT_OUT
              cout << "c_second_table_u->op1+(c_second_table_u->op1-c_first_table_u->op1)=" << c_second_table_u->op1 + (c_second_table_u->op1 - c_first_table_u->op1)
                   << "?=c_third_table_u->op1=" << c_third_table_u->op1 << "\n";
#endif
              Regular = false;
              break;
            }
          if((c_first_table_u->type != 3) && (c_first_table_u->type != 7))
            {
              if(c_second_table_u->op2 + (c_second_table_u->op2 - c_first_table_u->op2) != c_third_table_u->op2)
                {
#ifdef PRINT_OUT
                  cout << "c_second_table_u->op2+(c_second_table_u->op2-c_first_table_u->op2)=" << c_second_table_u->op2 + (c_second_table_u->op2 - c_first_table_u->op2)
                       << "?=c_third_table_u->op2=" << c_third_table_u->op2 << "\n";
#endif
                  Regular = false;
                  break;
                }
            }
        }
      if(c_first_table_u->pNext != second_u_blck->pNext)
        {
          c_second_table_u = c_second_table_u->pNext;
          c_first_table_u = c_first_table_u->pNext;
          c_third_table_u = c_third_table_u->pNext;
        }
      else
        OK = 0;
    }
#ifdef PRINT_OUT
  cout << "Regular=" << Regular << " OK=" << OK << " i=" << i << "\n";
  cout << "Check_Regularity Regular=" << Regular << "\n";
#endif
  return (Regular);
}


//==================================================================================
t_table_u*
SymbolicGaussElimination::interpolation(t_model_graph* model_graph, t_table_y* table_y, int to_add, bool middle, int per)
{
  t_table_u *tmp_table_u=NULL, *old_table_u=NULL;
  t_table_u *c_first_table_u, *c_second_table_u, *c_third_table_u;
  int c_first_y_blck, c_second_y_blck;
  int i, j, k, up_to=0, op_count, k1, k2;
  bool OK, modify_u_count;
  int cur_pos, nb_table_u=0;
#ifdef PRINT_OUT
  cout << "in interpolation\n";
  cout << "--- u_count=" << u_count << "\n";
  print_Graph(model_graph);
#endif /**PRINT_OUT**/
  cur_pos=SaveCode.tellp();
  SaveCode.write(reinterpret_cast<char *>(&up_to), sizeof(up_to));
  if(middle)
    {
      up_to = 2;
      middle_count_loop = 1 + per;
      //middle_count_loop = 4;
    }
  else
    up_to = 2;
#ifdef SAVE
  if(middle)
    {
      middle_save_table_u = NULL;
#ifdef PRINT_OUT
      cout << "up_to=" << up_to << " middle_count_loop=" << middle_count_loop << "\n";
#endif
    }
#endif /**SAVE**/
  /*Adding the Gaussian Elimination coefficients for the uncomputed periods*/
  /*table_u = tmp_table_u = (t_table_u*)malloc(sizeof(t_table_u) - 3 * sizeof(int));*/
  c_third_table_u = third_u_blck;
  c_first_y_blck = first_y_blck;
  c_second_y_blck = second_y_blck;
  modify_u_count = true;
  for(i = 1;i < up_to;i++)
    {
      c_first_table_u = first_u_blck->pNext;
      c_second_table_u = second_u_blck->pNext;
#ifdef PRINT_OUT
      cout << "c_first_table_u=" << c_first_table_u << " c_second_table_u=" << c_second_table_u << "\n";
#endif
      old_table_u = tmp_table_u;
      op_count = 0;
      OK = true;
      while(OK)
        {
          if(c_first_table_u->type != c_second_table_u->type)
            {
              cout << "Error : Interpolation the two lists are not synchronized (middle=" << middle << ", op_count=" << op_count << ")\n";
              cout << "c_first_table_u=" << c_first_table_u << " c_second_table_u=" << c_second_table_u << "\n";
              cout << "y_kmin=" << y_kmin << " y_kmax=" << y_kmax << "\n";
              cout << "first_u_blck=" << first_u_blck << " second_u_blck=" << second_u_blck << "\n";
              cout << "c_first_table_u->type=" << int(c_first_table_u->type) << " c_second_table_u->type=" << int(c_second_table_u->type) << "\n";
              system("pause");
              exit( -1);
            }
          switch (c_first_table_u->type)
            {
            case 3:
            case 7:
              nb_table_u++;
#ifdef SAVE
              if(i == 1)
                {
                  if(middle)
                    {
                      nb_middle_save_table_u++;
                    }
                  else
                    {
#ifdef PRINT_OUT
                      cout << "!!!!!nb_first_save_table_u=" << nb_first_save_table_u << "\n";
#endif
                      nb_first_save_table_u++;
                    }
                  SaveCode.write(reinterpret_cast<char *>(&c_first_table_u->type), sizeof(c_first_table_u->type));
                  SaveCode.write(reinterpret_cast<char *>(&c_first_table_u->index), sizeof(c_first_table_u->index));
                  SaveCode.write(reinterpret_cast<char *>(&c_first_table_u->op1), sizeof(c_first_table_u->op1));
#ifdef PRINT_OUT
                  if(c_first_table_u->type == 3)
                    cout << ">u[" << c_first_table_u->index << "]=1/(1-u[" << c_first_table_u->op1 << "])\n"; //"]) (" << tmp_table_u << ", " << tmp_table_u->pNext << ")\n";
                  else
                    cout << ">u[" << c_first_table_u->index << "]*=u[" << c_first_table_u->op1 << "]\n"; //"] (" << tmp_table_u << ", " << tmp_table_u->pNext << ")\n";
#endif
                  k=c_second_table_u->index - c_first_table_u->index;
                  SaveCode.write(reinterpret_cast<char *>(&k), sizeof(k));
                  k1=c_second_table_u->op1- c_first_table_u->op1;
                  SaveCode.write(reinterpret_cast<char *>(&k1), sizeof(k1));
#ifdef PRINT_OUT
                  if(c_first_table_u->type == 3)
                    cout << ">i u[" << k << "]=1/(1-u[" << k1 << "])\n"; //"]) (" << tmp_table_u << ", " << tmp_table_u->pNext << ")\n";
                  else
                    cout << ">i u[" << k << "]*=u[" << k1 << "]\n"; //"] (" << tmp_table_u << ", " << tmp_table_u->pNext << ")\n";
#endif
                }
#endif /**SAVE**/
              break;
            case 1:
            case 2:
            case 6:
              nb_table_u++;
#ifdef SAVE
              if((i == 1) && (middle))
                {
                  nb_middle_save_table_u++;
                  SaveCode.write(reinterpret_cast<char *>(&c_first_table_u->type), sizeof(c_first_table_u->type));
                  SaveCode.write(reinterpret_cast<char *>(&c_first_table_u->index), sizeof(c_first_table_u->index));
                  SaveCode.write(reinterpret_cast<char *>(&c_first_table_u->op1), sizeof(c_first_table_u->op1));
                  SaveCode.write(reinterpret_cast<char *>(&c_first_table_u->op2), sizeof(c_first_table_u->op2));
#ifdef PRINT_OUT
                  if((c_first_table_u->type == 1) && (middle))
                    cout << ">u[" << c_first_table_u->index << "]=u[" << c_first_table_u->op1 << "]*u[" << c_first_table_u->op2 << "]\n"; //"] (" << tmp_table_u << ", " << tmp_table_u->pNext << ")\n";
                  else if(c_first_table_u->type == 2)
                    cout << ">u[" << c_first_table_u->index << "]+=u[" << c_first_table_u->op1 << "]*u[" << c_first_table_u->op2 << "]\n"; //"] (" << tmp_table_u << ", " << tmp_table_u->pNext << ")\n";
                  else
                    cout << ">u[" << c_first_table_u->index << "]=1/(1-u[" << c_first_table_u->op1 << "]*u[" << c_first_table_u->op2 << "])\n"; //"]) (" << tmp_table_u << ", " << tmp_table_u->pNext << ")\n";
#endif //PRINT_OUT

                  k=c_second_table_u->index - c_first_table_u->index;
                  SaveCode.write(reinterpret_cast<char *>(&k), sizeof(c_first_table_u->index));
                  k1=c_second_table_u->op1- c_first_table_u->op1;
                  SaveCode.write(reinterpret_cast<char *>(&k1), sizeof(c_first_table_u->op1));
                  k2=c_second_table_u->op2- c_first_table_u->op2;
                  SaveCode.write(reinterpret_cast<char *>(&k2), sizeof(c_first_table_u->op2));
#ifdef PRINT_OUT
                  if((c_first_table_u->type == 1) && (middle))
                    cout << ">u[" << k << "]=u[" << k1 << "]*u[" << k2 << "]\n"; //"] (" << tmp_table_u << ", " << tmp_table_u->pNext << ")\n";
                  else if(c_first_table_u->type == 2)
                    cout << ">u[" << k << "]+=u[" << k1 << "]*u[" << k2 << "]\n"; //"] (" << tmp_table_u << ", " << tmp_table_u->pNext << ")\n";
                  else
                    cout << ">u[" << k << "]=1/(1-u[" << k1 << "]*u[" << k2 << "])\n"; //"]) (" << tmp_table_u << ", " << tmp_table_u->pNext << ")\n";
#endif //PRINT_OUT
                }
#endif /**SAVE**/
              break;
            case 5:
              nb_table_u++;
#ifdef SAVE
              if((i == 1) && (middle))
                {
                  nb_middle_save_table_u++;
                  SaveCode.write(reinterpret_cast<char *>(&c_first_table_u->type), sizeof(c_first_table_u->type));
                  SaveCode.write(reinterpret_cast<char *>(&c_first_table_u->index), sizeof(c_first_table_u->index));
#ifdef PRINT_OUT
                  cout << ">push(u[" << c_first_table_u->index << "])\n"; // "]) (" << tmp_table_u << ", " << tmp_table_u->pNext << ")\n";
#endif //PRINT_OUT
                  k=c_second_table_u->index - c_first_table_u->index;
                  SaveCode.write(reinterpret_cast<char *>(&k), sizeof(c_first_table_u->index));
#ifdef PRINT_OUT
                  cout << ">push(u[" << k << "])\n"; // "]) (" << tmp_table_u << ", " << tmp_table_u->pNext << ")\n";
#endif //PRINT_OUT

                }
#endif /**SAVE**/
              break;
            }
          if((c_first_table_u->pNext != second_u_blck->pNext) && (c_first_table_u->pNext != NULL) && (c_second_table_u->pNext != NULL))
            {
              tmp_table_u=c_first_table_u ;
              c_first_table_u = c_first_table_u->pNext;
              free(tmp_table_u);
              tmp_table_u=c_second_table_u ;
              c_second_table_u = c_second_table_u->pNext;
              free(tmp_table_u);
            }
          else
            {
              if(c_first_table_u->pNext != second_u_blck->pNext) /*||(second_u_blck->pNext!=NULL))*/
                {
                  cout << "c_first_table_u->pNext=" << c_first_table_u->pNext << " second_u_blck->pNext=" << second_u_blck->pNext << "\n";
                  cout << "Error not synchronize graph interpolation\n";
                  exit( -1);
                }
              OK = 0;
            }
          op_count++;
        }
      SaveCode.flush();
      k=SaveCode.tellp();
      SaveCode.seekp(cur_pos);
      SaveCode.write(reinterpret_cast<char *>(&nb_table_u), sizeof(nb_table_u));
      SaveCode.seekp(k);
      if(middle)
        {
          first_y_blck = c_first_y_blck;
          second_y_blck = c_second_y_blck;
#ifdef SAVE
          if((i == 1) && (middle))
            {
              middle_save_table_y = (t_table_y*)malloc((nb_endo) * sizeof(t_table_y));
              middle_save_i_table_y = (t_table_y*)malloc((nb_endo) * sizeof(t_table_y));
              nb_middle_save_table_y = nb_endo;
            }
#endif /**SAVE**/
          for(j = 0;j < nb_endo;j++)
            {
#ifdef SAVE
              if(i == 1)
                {
                  if(middle)
                    {
                      middle_save_i_table_y[j].u_index = (int*)malloc(table_y[first_y_blck].nb * sizeof(int));
                      middle_save_i_table_y[j].y_index = (int*)malloc(table_y[first_y_blck].nb * sizeof(int));
                      middle_save_i_table_y[j].index = table_y[second_y_blck].index - table_y[first_y_blck].index;
                      middle_save_i_table_y[j].nb = table_y[first_y_blck].nb;
                      for(k = 0;k < table_y[first_y_blck].nb;k++)
                        {
                          middle_save_i_table_y[j].u_index[k] = table_y[second_y_blck].u_index[k] - table_y[first_y_blck].u_index[k];
                          middle_save_i_table_y[j].y_index[k] = table_y[second_y_blck].y_index[k] - table_y[first_y_blck].y_index[k];
                        }
                      middle_save_table_y[j].u_index = (int*)malloc(table_y[first_y_blck].nb * sizeof(int));
                      middle_save_table_y[j].y_index = (int*)malloc(table_y[first_y_blck].nb * sizeof(int));
                      /*middle_save_table_y[j].y_lead_lag = (int*)malloc(table_y[first_y_blck].nb * sizeof(int));*/
                      middle_save_table_y[j].index = table_y[first_y_blck].index;
                      middle_save_table_y[j].nb = table_y[first_y_blck].nb;
                      for(k = 0;k < table_y[first_y_blck].nb;k++)
                        {
                          middle_save_table_y[j].u_index[k] = table_y[first_y_blck].u_index[k];
                          middle_save_table_y[j].y_index[k] = table_y[first_y_blck].y_index[k];
                        }
                    }
                  else
                    {
                      first_save_i_table_y[j].u_index = (int*)malloc(table_y[first_y_blck].nb * sizeof(int));
                      first_save_i_table_y[j].y_index = (int*)malloc(table_y[first_y_blck].nb * sizeof(int));
                      first_save_i_table_y[j].index = table_y[second_y_blck].index - table_y[first_y_blck].index;
                      first_save_i_table_y[j].nb = table_y[first_y_blck].nb;
                      for(k = 0;k < table_y[first_y_blck].nb;k++)
                        {
                          first_save_i_table_y[j].u_index[k] = table_y[second_y_blck].u_index[k] - table_y[first_y_blck].u_index[k];
                          first_save_i_table_y[j].y_index[k] = table_y[second_y_blck].y_index[k] - table_y[first_y_blck].y_index[k];
                        }
                      first_save_table_y[j].u_index = (int*)malloc(table_y[first_y_blck].nb * sizeof(int));
                      first_save_table_y[j].y_index = (int*)malloc(table_y[first_y_blck].nb * sizeof(int));
                      first_save_table_y[j].index = table_y[first_y_blck].index;
                      first_save_table_y[j].nb = table_y[first_y_blck].nb;
                      for(k = 0;k < table_y[first_y_blck].nb;k++)
                        {
                          first_save_table_y[j].u_index[k] = table_y[first_y_blck].u_index[k];
                          first_save_table_y[j].y_index[k] = table_y[first_y_blck].y_index[k];
                        }
                    }
                }
#endif /**SAVE**/
              table_y[vertex_count].nb = table_y[first_y_blck].nb;
              table_y[vertex_count].index = table_y[first_y_blck].index + (i - 1) * (table_y[second_y_blck].index - table_y[first_y_blck].index);
#ifdef PRINT_OUT
              cout << ">y[" << table_y[vertex_count].index << "]=";
#endif /**PRINT_OUT**/
              for(k = 0;k < table_y[vertex_count].nb;k++)
                {
                  table_y[vertex_count].u_index[k] = table_y[first_y_blck].u_index[k] + (i - 1) * (table_y[second_y_blck].u_index[k] - table_y[first_y_blck].u_index[k]);
                  table_y[vertex_count].y_index[k] = table_y[first_y_blck].y_index[k] + (i - 1) * (table_y[second_y_blck].y_index[k] - table_y[first_y_blck].y_index[k]);
#ifdef PRINT_OUT
                  cout << "u[" << table_y[vertex_count].u_index[k] << "]*y[" << table_y[vertex_count].y_index[k] << "]";
                  if(k + 1 == table_y[vertex_count].nb)
                    cout << "\n";
                  else
                    cout << "+";
#endif /**PRINT_OUT**/
                }
              vertex_count++;
              first_y_blck++;
              second_y_blck++;
            }
        }
    }
#ifdef PRINT_OUT
  cout << "+++ u_count=" << u_count << "\n";
  cout << "out interpolation\n";
#endif

  return (old_table_u);
}


//==================================================================================
bool
SymbolicGaussElimination::Loop_Elimination(t_model_graph* model_graph)
{
  int i, j, pos, i1, i_per;
  bool there_is_a_loop, go_on = true;
  /* destroy the loop*/
  int per = 0;
  there_is_a_loop = 0;
  i_per = 0;
  for(i = 0;i < model_graph->nb_vertices;i++)
    {
      if(nstacked)
        {
          if(model_graph->vertex[i].nb_in_degree_edges > 0)
            {
              if(!(i_per % nb_endo))
                {
                  per++;
                  if(per == 1)
                    {
                      first_u_blck = table_u;
#ifdef PRINT_OUT
                      cout << "first_u_blck=" << first_u_blck << "\n";
#endif
                    }
                  if(per == 2)
                    {
                      pos_nb_first_save_table_u = nb_table_u;
                      second_u_blck = table_u;
#ifdef PRINT_OUT
                      cout << "second_u_blck=" << second_u_blck << "\n";
#endif
                    }
                  if(per == 3)
                    {
                      go_on = false;
                      third_u_blck = table_u;
#ifdef PRINT_OUT
                      cout << "third_u_blck\n";
#endif
                    }
                }
              i_per++;
            }
        }
      pos = -1;
      for(j = 0;((j < model_graph->vertex[i].nb_in_degree_edges) && (pos < 0));j++)
        if(model_graph->vertex[i].in_degree_edge[j].index == i)
          {
            // at least one loop
            pos = j;
          }
      if(pos >= 0)
        {
          there_is_a_loop = 1;
#ifdef DEBUGR
          cout << "--------------------------------------------------------------- \n";
          cout << "Loop elimination on vertex " << model_graph->vertex[i].index << "\n";
#endif /**DEBUGR**/
          if(nstacked)
            {
              /*Store the u_count associated to the loop to re-use it in case of
                loop creation during vertex supression step*/
              loop_table_u_count[i] = model_graph->vertex[i].in_degree_edge[pos].u_count;
              loop_table_vertex_index[i] = model_graph->vertex[i].index;
            }
          if(go_on)
            {
              nb_table_u++;
              table_u->pNext = (t_table_u*)malloc(sizeof(t_table_u) - sizeof(int));
              table_u = table_u->pNext;
              table_u->pNext = NULL;
              table_u->type = 3;
              table_u->index = model_graph->vertex[i].in_degree_edge[pos].u_count;
              table_u->op1 = model_graph->vertex[i].in_degree_edge[pos].u_count;
#ifdef PRINT_OUT
              cout << "u[" << table_u->op1 << "]=1/(-u[" << table_u->op1 << "]) " << table_u << "\n"; // "])  (nb_free_u_list : " << nb_free_u_list  << " u_count=" << u_count << ")\n";
#endif /**PRINT_OUT**/
              for(j = 0;j < model_graph->vertex[i].nb_in_degree_edges;j++)
                {
                  if(j != pos)
                    {
                      i1 = model_graph->vertex[i].in_degree_edge[j].index;
#ifdef DIRECT_COMPUTE
                      nb_table_u++;
                      table_u->pNext = (t_table_u*)malloc(sizeof(t_table_u) - sizeof(int));
                      table_u = table_u->pNext;
                      table_u->pNext = NULL;
                      table_u->type = 7;
                      table_u->index = model_graph->vertex[i].in_degree_edge[j].u_count;
                      table_u->op1 = model_graph->vertex[i].in_degree_edge[pos].u_count;
#ifdef PRINT_OUT
                      cout << "u[" << model_graph->vertex[i].in_degree_edge[j].u_count << "]*=u[" << table_u->op1 << "] " << table_u << "\n"; // "] (nb_free_u_list : " << nb_free_u_list  << " u_count=" << u_count << ")\n";
#endif /**PRINT_OUT**/
#endif /**DIRECT_COMPUTE**/
                    }
                }
            }
          /*elimination of the loop*/
          /*in the out_degree*/
          pos = -1;
          for(j = 0;((j < model_graph->vertex[i].nb_out_degree_edges) && (pos < 0));j++)
            if(model_graph->vertex[i].out_degree_edge[j].index == i)
              pos = j;
          if(pos < 0)
            {
              cout << "Error: not symetric on a loop on vertex " << model_graph->vertex[i].index << "\n";
              print_Graph(model_graph);
              system("pause");
              exit( -1);
            }
          for(j = pos + 1;j < model_graph->vertex[i].nb_out_degree_edges;j++)
            {
              model_graph->vertex[i].out_degree_edge[j - 1].index = model_graph->vertex[i].out_degree_edge[j].index;
              model_graph->vertex[i].out_degree_edge[j - 1].u_count = model_graph->vertex[i].out_degree_edge[j].u_count;
            }
          model_graph->vertex[i].nb_out_degree_edges--;
          /*in the in_degree*/
          pos = -1;
          for(j = 0;((j < model_graph->vertex[i].nb_in_degree_edges) && (pos < 0));j++)
            if(model_graph->vertex[i].in_degree_edge[j].index == i)
              pos = j;
          if(pos < 0)
            {
              cout << "Error: not symetric on a loop on vertex " << model_graph->vertex[i].index << "\n";
              print_Graph(model_graph);
              system("pause");
              exit( -1);
            }
          for(j = pos + 1;j < model_graph->vertex[i].nb_in_degree_edges;j++)
            {
              model_graph->vertex[i].in_degree_edge[j - 1].index = model_graph->vertex[i].in_degree_edge[j].index;
              model_graph->vertex[i].in_degree_edge[j - 1].u_count = model_graph->vertex[i].in_degree_edge[j].u_count;
            }
          model_graph->vertex[i].nb_in_degree_edges--;
        }
    }
  return (there_is_a_loop);
}

//==================================================================================
bool
SymbolicGaussElimination::Vertex_Elimination(t_model_graph* model_graph, int pos, bool* interpolate, int length_markowitz, bool dynamic)
{
  int i, j, k, min_edge, vertex_to_eliminate = -1, to_add, to_add_index, curr_u_count;
  int i1, i2, j1, j2, i3, i4, i5, size;
  int i_max, j_max;
  bool OK;
  t_vertex* lvertex = model_graph->vertex;
  t_vertex* ilvertex;
  int nb_new = 0, nb_free = 0;
  int a_loop = 0;
  char usi;
#ifdef DIRECT_COMPUTE
  int in_y = 0;
#endif /**DIRECT_COMPUTE**/
  size = model_graph->nb_vertices;
  min_edge = size * size + 1;
  vertex_to_eliminate = -1;
  // minimize the product of in and out degree of the remaining vertices (Markowitz criteria)
  for(i = starting_vertex;i < starting_vertex + length_markowitz;i++)
    if(((lvertex[i].nb_out_degree_edges*lvertex[i].nb_in_degree_edges) < min_edge) && (lvertex[i].nb_in_degree_edges > 0))
      {
        min_edge = lvertex[i].nb_out_degree_edges * lvertex[i].nb_in_degree_edges;
        vertex_to_eliminate = i;
      }
  OK = 0;
  //cout << "save_direct=" << save_direct << "\n";
  if(vertex_to_eliminate >= 0)
    {
#ifdef PRINT_OUT_1
      cout << "--------------------------------------------------------------- \n";
      cout << "Elimination of vertex " << lvertex[vertex_to_eliminate].index << " interpolate=" << *interpolate << "\n";
      cout << "min_edge=" << min_edge << " length_markowitz=" << length_markowitz << "\n";
      print_Graph(model_graph);
#endif /**PRINT_OUT**/
#ifdef DIRECT_COMPUTE
      table_y[vertex_count].index = lvertex[vertex_to_eliminate].index;
#ifdef PRINT_OUT
      cout << "vertex_count=" << vertex_count << " size=" << size << " vertex_to_eliminate=" << vertex_to_eliminate << "\n";
#endif
#endif /**DIRECT_COMPUTE**/
      ilvertex = &(lvertex[vertex_to_eliminate]);
      i_max = ilvertex->nb_in_degree_edges;
      j_max = ilvertex->nb_out_degree_edges;
      for(i = 0;i < i_max;i++)
        {
          i1 = ilvertex->in_degree_edge[i].index;
          i2 = ilvertex->in_degree_edge[i].u_count;
#ifdef DIRECT_COMPUTE
          table_y[vertex_count].u_index[in_y] = i2;
          table_y[vertex_count].y_index[in_y] = lvertex[i1].index;
          in_y++;
          nb_table_u++;
          if(save_direct)
            {
              usi=5;
              //SaveCode << (&usi) << i2;
              SaveCode.write(reinterpret_cast<char *>(&usi), sizeof(usi));
              SaveCode.write(reinterpret_cast<char *>(&i2), sizeof(i2));
#ifdef PRINT_OUT_1
              cout << "push(u[" << i2 << "]) \n";
#endif
            }
          else
            {
              table_u->pNext = (t_table_u*)malloc(sizeof(t_table_u) - 2 * sizeof(int));
              table_u = table_u->pNext;
              table_u->pNext = NULL;
              table_u->type = 5;
              table_u->index = i2;
            }
          //cout << "ok0\n";
#ifdef PRINT_OUT_1
          cout << "push(u[" << i2 << "]) (" << table_u << ")\n";
#endif /**PRINT_OUT**/
#endif /**DIRECT_COMPUTE**/
          for(j = 0;j < j_max;j++)
            {
              to_add = -9999999;
              j1 = ilvertex->out_degree_edge[j].index;
              j2 = ilvertex->out_degree_edge[j].u_count;
              //cout << "ok1\n";
              /* Is there already an edge from vertex i1 to vertex j1? */
              for(k = 0;((k < lvertex[i1].nb_out_degree_edges) && (to_add < 0));k++)
                if(lvertex[i1].out_degree_edge[k].index == j1)
                  {
                    /*yes*/
                    to_add = lvertex[i1].out_degree_edge[k].u_count;
                    to_add_index = k;
                  }
              if(to_add != -9999999)
                {
#ifdef DEBUGR
                  cout << " modification of edge between vertices " << lvertex[i1].index << " and " << lvertex[j1].index << "\n";
#endif /**DEBUGR**/
#ifdef DIRECT_COMPUTE
                  if(save_direct)
                    {
                      nb_table_u++;
                      usi=2;
                      SaveCode.write(reinterpret_cast<char *>(&usi), sizeof(usi));
                      SaveCode.write(reinterpret_cast<char *>(&to_add), sizeof(to_add));
                      SaveCode.write(reinterpret_cast<char *>(&j2), sizeof(j2));
                      SaveCode.write(reinterpret_cast<char *>(&i2), sizeof(i2));
                      //SaveCode << (unsigned short int)(2) << to_add  << j2 << i2;
#ifdef PRINT_OUT_1
                      cout << "u[" << to_add << "]+=u[" << j2 << "]*u[" << i2 << "]; \n";
#endif
                    }
                  else
                    {
                      nb_table_u++;
                      table_u->pNext = (t_table_u*)malloc(sizeof(t_table_u));
                      table_u = table_u->pNext;
                      table_u->pNext = NULL;
                      table_u->type = 2;
                      table_u->index = to_add;
                      table_u->op1 = j2;
                      table_u->op2 = i2;
                    }
#ifdef PRINT_OUT_1
                  cout << "u[" << to_add << "]+=u[" << j2 << "]*u[" << i2 << "]; (" << table_u << ")\n";
#endif /**PRINT_OUT**/
#endif /**DIRECT_COMPUTE**/
                }
              else
                {
                  /*now modifie the in_degree edge from vertex_to_elminate to j1 */
                  if(j1 == i1)
                    {
                      /* its a loop  */
                      /*Store it to operate after*/
                      if(a_loop >= size)
                        {
                          cout << "Error : a_loop (" << a_loop << ") >= " << size << "\n";
                          exit( -1);
                        }
                      s_j2[a_loop] = j2;
                      s_i2[a_loop] = i2;
                      s_i1[a_loop] = i1;
                      a_loop++;
                    }
                  else
                    {
#ifdef DEBUGR
                      cout << " creation of edge between vertices " << lvertex[i1].index << " and " << lvertex[j1].index << "\n";
#endif /**DEBUGR**/
#ifdef SYMPLIFY  /*no reuse of free numbers, to be sure that there is no cyclcial behaviour in numbering*/
                      curr_u_count = get_free_u_list(dynamic);
#else
                      curr_u_count = u_count;
                      u_count++;
#endif
#ifdef DIRECT_COMPUTE
                      nb_table_u++;
                      nb_new++;
                      if(save_direct)
                        {
                          usi=1;
                          SaveCode.write(reinterpret_cast<char *>(&usi), sizeof(usi));
                          SaveCode.write(reinterpret_cast<char *>(&curr_u_count), sizeof(curr_u_count));
                          SaveCode.write(reinterpret_cast<char *>(&j2), sizeof(j2));
                          SaveCode.write(reinterpret_cast<char *>(&i2), sizeof(i2));
#ifdef PRINT_OUT_1
                          cout << "u[" << curr_u_count << "]=u[" << j2 << "]*u[" << i2 << "]; \n";
#endif
                        }
                      else
                        {
                          table_u->pNext = (t_table_u*)malloc(sizeof(t_table_u));
                          table_u = table_u->pNext;
                          table_u->pNext = NULL;
                          table_u->type = 1;
                          table_u->index = curr_u_count;
                          table_u->op1 = j2;
                          table_u->op2 = i2;
                        }
#ifdef PRINT_OUT_1
                      cout << "u[" << curr_u_count << "]=u[" << j2 << "]*u[" << i2 << "]; (" << table_u << ")\n";
#endif /**PRINT_OUT**/
#endif /**DIRECT_COMPUTE**/
                      /*its not a loop so adding a new adge*/
                      /*modify the in_degree edge from vertex_to_elminate to j1*/
#ifdef SORTED
                      k = lvertex[j1].nb_in_degree_edges;
                      while((i1 < lvertex[j1].in_degree_edge[k - 1].index) && (k > 0))
                        {
                          lvertex[j1].in_degree_edge[k].index = lvertex[j1].in_degree_edge[k - 1].index;
                          lvertex[j1].in_degree_edge[k].u_count = lvertex[j1].in_degree_edge[k - 1].u_count;
                          k--;
                        }
                      lvertex[j1].in_degree_edge[k].index = i1;
                      lvertex[j1].in_degree_edge[k].u_count = curr_u_count;
                      (lvertex[j1].nb_in_degree_edges)++;
                      /*Tested*/
                      if(lvertex[j1].max_nb_in_degree_edges<=lvertex[j1].nb_in_degree_edges)
                        {
                          lvertex[j1].max_nb_in_degree_edges=lvertex[j1].nb_in_degree_edges;
                          t_edge *tmp_edge;
                          tmp_edge=lvertex[j1].in_degree_edge;
                          int taille=lvertex[j1].nb_in_degree_edges*sizeof(t_edge);
                          tmp_edge = (t_edge*)realloc(tmp_edge, taille);
                          lvertex[j1].in_degree_edge=tmp_edge;
                        }
                      /*EndTested*/
#else
                      k = lvertex[j1].nb_in_degree_edges;
                      lvertex[j1].in_degree_edge[k].index = i1;
                      lvertex[j1].in_degree_edge[k].u_count = curr_u_count;
                      (lvertex[j1].nb_in_degree_edges)++;
#endif
                      /*now modify the out_degree edge from i1 to vertex_to_elminate */
#ifdef SORTED
                      k = lvertex[i1].nb_out_degree_edges;
                      while((j1 < lvertex[i1].out_degree_edge[k - 1].index) && (k > 0))
                        {
                          lvertex[i1].out_degree_edge[k].index = lvertex[i1].out_degree_edge[k - 1].index;
                          lvertex[i1].out_degree_edge[k].u_count = lvertex[i1].out_degree_edge[k - 1].u_count;
                          k--;
                        }
                      lvertex[i1].out_degree_edge[k].index = j1;
                      lvertex[i1].out_degree_edge[k].u_count = curr_u_count;
                      (lvertex[i1].nb_out_degree_edges)++;
                      /*Tested*/
                      if(lvertex[j1].max_nb_out_degree_edges<lvertex[j1].nb_out_degree_edges)
                        {
                          lvertex[j1].max_nb_out_degree_edges=lvertex[j1].nb_out_degree_edges;
                          lvertex[j1].out_degree_edge=(t_edge*)realloc(lvertex[j1].out_degree_edge,lvertex[j1].nb_out_degree_edges*sizeof(t_edge));
                        }
                      /*EndTested*/
#else
                      k = lvertex[i1].nb_out_degree_edges;
                      lvertex[i1].out_degree_edge[k].index = j1;
                      lvertex[i1].out_degree_edge[k].u_count = curr_u_count;
                      (lvertex[i1].nb_out_degree_edges)++;
#endif
                    }
                }

            }
        }
#ifdef SIMPLIFY
      nb_free_u_list1 = 0;
#endif /**SIMPLIFY**/
      for(i = 0;i < i_max;i++)
        {
          /*now supress the out_degree edge from i1 to vertex_to_elminate */
          i1 = ilvertex->in_degree_edge[i].index;
#ifdef DEBUGR
          cout << " supress the out_degree edge from " << lvertex[i1].index << " to " << ilvertex->index << "\n";
#endif /**DEBUGR**/
          to_add = -1;
          for(k = 0;((k < lvertex[i1].nb_out_degree_edges) && (to_add < 0));k++)
            if(lvertex[i1].out_degree_edge[k].index == vertex_to_eliminate)
              to_add = k;
          if(to_add < 0)
            {
              cout << "error: Model_Graph not correctly filled in out_degree (edge from " << lvertex[i1].index << " to " << lvertex[vertex_to_eliminate].index << ")\n";
              print_Graph(model_graph);
              system("pause");
              exit( -1);
            }
#ifdef SIMPLIFYS
          nb_free++;
          set_free_u_list(lvertex[i1].out_degree_edge[to_add].u_count);
#endif /**SIMPLIFY**/
#ifdef SIMPLIFYT
          if(simplification_allowed)
            {
              free_u_list1[nb_free_u_list1] = lvertex[i1].out_degree_edge[to_add].u_count;
              if(nb_free_u_list1 > max_nb_free_u_list1)
                max_nb_free_u_list1 = nb_free_u_list1;
              nb_free_u_list1++;
            }
#endif /**SIMPLIFY**/
          for(k = to_add + 1;k < lvertex[i1].nb_out_degree_edges;k++)
            {
              lvertex[i1].out_degree_edge[k - 1].index = lvertex[i1].out_degree_edge[k].index;
              lvertex[i1].out_degree_edge[k - 1].u_count = lvertex[i1].out_degree_edge[k].u_count;
            }
          (lvertex[i1].nb_out_degree_edges)--;
        }
      for(j = 0;j < j_max;j++)
        {
          /*now supress the in_degree edge from vertex_to_elminate to j1 */
          j1 = ilvertex->out_degree_edge[j].index;
#ifdef DEBUGR
          cout << " supress the in_degree edge from " << ilvertex->index << " to " << lvertex[j1].index << "\n";
#endif /**DEBUGR**/
          to_add = -1;
          for(k = 0;((k < lvertex[j1].nb_in_degree_edges) && (to_add < 0));k++)
            if(lvertex[j1].in_degree_edge[k].index == vertex_to_eliminate)
              to_add = k;
          if(to_add < 0)
            {
              cout << "error: Model_Graph not correctly filled in in_degree (edge from " << lvertex[vertex_to_eliminate].index << " to " << lvertex[j1].index << ")\n";
              print_Graph(model_graph);
              system("pause");
              exit( -1);
            }
#ifdef SIMPLIFYS
          set_free_u_list(lvertex[j1].in_degree_edge[to_add].u_count);
          nb_free++;
#endif /**SIMPLIFY**/
#ifdef SIMPLIFYT
          if(simplification_allowed)
            {
              free_u_list1[nb_free_u_list1] = lvertex[j1].in_degree_edge[to_add].u_count;
              if(nb_free_u_list1 > max_nb_free_u_list1)
                max_nb_free_u_list1 = nb_free_u_list1;
              nb_free_u_list1++;
            }
#endif /**SIMPLIFY**/
          for(k = to_add + 1;k < lvertex[j1].nb_in_degree_edges;k++)
            {
              lvertex[j1].in_degree_edge[k - 1].index = lvertex[j1].in_degree_edge[k].index;
              lvertex[j1].in_degree_edge[k - 1].u_count = lvertex[j1].in_degree_edge[k].u_count;
            }
          (lvertex[j1].nb_in_degree_edges)--;
        }
      if(a_loop)
        {
          for(i = 0;i < a_loop;i++)
            {
              i1 = s_i1[i];
              if(nstacked)
                {
                  /*re-use the u_count index eliminated during the loop elmination*/
                  k = 0;
                  while((k <= nb_loop_table) && (loop_table_vertex_index[k] != model_graph->vertex[i1].index))
                    k++;
                  if(loop_table_vertex_index[k] != model_graph->vertex[i1].index)
                    {
                      cout << "Error can't find the loop associated to vertex " << model_graph->vertex[i1].index << "\n";
                      k = 0;
                      while((k <= nb_loop_table) && (loop_table_vertex_index[k] != model_graph->vertex[i1].index))
                        {
                          cout << "loop_table_vertex_index[" << k << "]=" << loop_table_vertex_index[k] << " =? model_graph->vertex[" << i1 << "].index=" << model_graph->vertex[i1].index << "\n";
                          k++;
                        }
                      exit( -1);
                    }
                  curr_u_count = loop_table_u_count[k] ;
                }
              else
                {
#ifdef SIMPLIFY
                  curr_u_count = get_free_u_list(dynamic);
#else /**SIMPLIFY**/
                  curr_u_count = u_count;
                  u_count++;
#endif /**SIMPLIFY**/
                }
              nb_table_u++;
              if(save_direct)
                {
                  usi=6;
                  SaveCode.write(reinterpret_cast<char *>(&usi), sizeof(usi));
                  SaveCode.write(reinterpret_cast<char *>(&curr_u_count), sizeof(curr_u_count));
                  SaveCode.write(reinterpret_cast<char *>(&s_i2[i]), sizeof(s_i2[i]));
                  SaveCode.write(reinterpret_cast<char *>(&s_j2[i]), sizeof(s_j2[i]));
#ifdef PRINT_OUT_1
                  cout << "u[" << curr_u_count << "]=1/(1-u[" << s_i2[i] << "]*u[" << s_j2[i] << "]) ;\n";
#endif
                }
              else
                {
                  table_u->pNext = (t_table_u*)malloc(sizeof(t_table_u) );
                  table_u = table_u->pNext;
                  table_u->pNext = NULL;
                  table_u->type = 6;
                  table_u->index = curr_u_count;
                  table_u->op1 = s_i2[i];
                  table_u->op2 = s_j2[i];
                }
              nb_new++;
#ifdef PRINT_OUT_1
              cout << "u[" << curr_u_count << "]=1/(1-u[" << s_i2[i] << "]*u[" << s_j2[i] << "]) (" << table_u << ");\n";
#endif /**PRINT_OUT**/
              i5 = curr_u_count;
              OK = 0;
              for(k = 0;k < lvertex[i1].nb_in_degree_edges;k++)
                {
                  i3 = lvertex[i1].in_degree_edge[k].index;
                  if(i3 != i1)
                    {
                      i4 = lvertex[i1].in_degree_edge[k].u_count;
#ifdef DIRECT_COMPUTE
                      nb_table_u++;
                      if(save_direct)
                        {
                          usi=7;
                          SaveCode.write(reinterpret_cast<char *>(&usi), sizeof(usi));
                          SaveCode.write(reinterpret_cast<char *>(&i4), sizeof(i4));
                          SaveCode.write(reinterpret_cast<char *>(&i5), sizeof(i5));
#ifdef PRINT_OUT_1
                          cout << "u[" << i4 << "]*=u[" << i5 << "]; \n";
#endif
                        }
                      else
                        {
                          table_u->pNext = (t_table_u*)malloc(sizeof(t_table_u) - sizeof(int));
                          table_u = table_u->pNext;
                          table_u->pNext = NULL;
                          table_u->type = 7;
                          table_u->index = i4;
                          table_u->op1 = i5;
                        }
#ifdef PRINT_OUT_1
                      cout << "u[" << i4 << "]*=u[" << i5 << "]; (" << table_u << ")\n";
#endif /**PRINT_OUT**/
#endif /**DIRECT_COMPUTE**/
                    }
                }
#ifdef SIMPLIFY
              if(simplification_allowed)
                {
                  set_free_u_list(i5);
                  nb_free++;
                }
#endif /**SIMPLIFY**/
            }
        }
#ifdef SIMPLIFYS
      if(simplification_allowed)
        for(i = 0;i < nb_free_u_list1;i++)
          {
            set_free_u_list(free_u_list1[i]);
            nb_free++;
          }
#endif /**SIMPLIFY**/
      /*Change all indexes*/
      free(ilvertex->in_degree_edge);
      free(ilvertex->out_degree_edge);
      for(i = vertex_to_eliminate + 1;i < size;i++)
        lvertex[i - 1] = lvertex[i];
      (model_graph->nb_vertices)--;
      OK = 0;
      for(i = 0;i < model_graph->nb_vertices;i++)
        {
          for(j = 0;j < lvertex[i].nb_in_degree_edges;j++)
            {
              OK = 1;
              if(lvertex[i].in_degree_edge[j].index > vertex_to_eliminate)
                (lvertex[i].in_degree_edge[j].index)--;
            }
          for(j = 0;j < lvertex[i].nb_out_degree_edges;j++)
            {
              if(lvertex[i].out_degree_edge[j].index > vertex_to_eliminate)
                (lvertex[i].out_degree_edge[j].index)--;
            }
        }
      if(in_y>max_nb_table_y)
        max_nb_table_y=in_y;
#ifdef DIRECT_COMPUTE
      table_y[vertex_count].nb = in_y;
#ifdef PRINT_OUT_1
      cout << "y[" << table_y[vertex_count].index << "]=";
      for(i = 0;i < table_y[vertex_count].nb;i++)
        {
          cout << "u[" << table_y[vertex_count].u_index[i] << "]*y[" << table_y[vertex_count].y_index[i] << "]";
          if(i + 1 < table_y[vertex_count].nb)
            cout << "+";
          else
            cout << "\n";
        }
#endif /**PRINT_OUT**/
      vertex_count++;
#endif /**DIRECT_COMPUTE**/
    }
#ifdef PRINT_OUT
  cout << "nb_new=" << nb_new << " nb_free=" << nb_free << "\n";
  cout <<"end of vertex elimination\n";
#endif
  return (OK);
}


//==================================================================================
void
SymbolicGaussElimination::Gaussian_Elimination(t_model_graph* model_graph
#ifdef SAVE
                                               , string file_name
#endif
                                               , bool dynamic)
{
  int i, size, j, k, per;
  bool OK, interpolate = false;
  t_table_u *First_prologue_table_u;
  int length_markowitz=0, per_next = 3;
  bool try_to_interpolate=false;
  int cur_pos=0;
  int prologue_nb_table_u=0, first_nb_prologue_save_table_y=0;
  int nb_first_u_blck, nb_second_u_blck=0, nb_third_u_blck=0;
  int nb_first_y_blck, nb_second_y_blck=0, nb_third_y_blck=0;

  size = model_graph->nb_vertices;
#ifdef PRINT_OUT
  cout << "going to open file file_open=" << file_open << " file_name={" << file_name << "}\n";
#endif
  if(file_open)
    SaveCode.open((file_name + ".bin").c_str(), ios::out | ios::in | ios::binary | ios ::ate );
  else
    SaveCode.open((file_name + ".bin").c_str(), ios::out | ios::binary);
  file_open = true;
  if(!SaveCode.is_open())
    {
      cout << "Error : Can't open file \"" << file_name << ".bin\" for writing\n";
      exit( -1);
    }
#ifdef PRINT_OUT
  print_Graph(model_graph);
#endif /**PRINT_OUT**/
  s_i1 = (int*)malloc(size * sizeof(int));
  s_i2 = (int*)malloc(size * sizeof(int));
  s_j2 = (int*)malloc(size * sizeof(int));
#ifdef SIMPLIFY
  simplification_allowed = 1;
  MAX_FREE_U_LIST = size * /*size *//*12*/50;
  free_u_list = (int*)malloc(MAX_FREE_U_LIST* sizeof(int));
#ifdef PRINT_OUT
  cout << "free_u_list=" << free_u_list << " length=" << ceil((double)4 / (double)2*(double)u_count)*sizeof(int) << "\n";
  cout << "free_u_list length0=" << (*((short int*)&(*(((char*)free_u_list) - 4))) & ~(3)) << "\n";
  cout << "free_u_list1=(int*)malloc(" << u_count << ");\n";
#endif
  free_u_list1 = (int*)malloc(u_count * sizeof(int));
#ifdef PRINT_OUT
  cout << "after free alloc u_count=" << u_count << "\n";
#endif
#endif /**SIMPLIFY**/
#ifdef PRINT_OUT
  cout << "periods=" << periods << "\n";
#endif
  if(dynamic)
    u_count = (periods + y_kmin + y_kmax) * u_count / (y_kmin + y_kmax + 2);
  if(nstacked)
    {
      nb_loop_table = model_graph->nb_vertices;
      loop_table_u_count = (int*)malloc(nb_loop_table * sizeof(int));
      loop_table_vertex_index = (int*)malloc(nb_loop_table * sizeof(int));
    }
#ifdef PRINT_OUT
  cout << "after nb_loop_table=" << nb_loop_table  << "\n";
  system("pause");
#endif

  SaveCode.write(reinterpret_cast<char *>(&nb_endo), sizeof(nb_endo));
  SaveCode.write(reinterpret_cast<char *>(&u_count), sizeof(u_count));
  SaveCode.write(reinterpret_cast<char *>(&u_count_init), sizeof(u_count_init));


#ifdef PRINT_OUT
  cout << "u_count=" << u_count << "\n";
  //We start with the elimination of all loops in the graph
  cout << "going to loop elminate\n";
#endif
  save_direct=false;
  if(nstacked)
    save_direct=false;
  else
    {
      i=0;
      SaveCode.write(reinterpret_cast<char *>(&i), sizeof(i));
      SaveCode.write(reinterpret_cast<char *>(&i), sizeof(i));
      SaveCode.write(reinterpret_cast<char *>(&i), sizeof(i));
      SaveCode.write(reinterpret_cast<char *>(&i), sizeof(i));
    }
  if(Loop_Elimination(model_graph))
    {
#ifdef PRINT_OUT
      cout << "->table_u=" << table_u << "\n";
#endif
      if(nstacked)
        {
          //if there are some loops eliminated, we have to update it for the next blocks
          // So we call the interpolation procedure
          //cout << "nb_first_table_u=";
          interpolation(model_graph, table_y,     1, false, 0);
        }
#ifdef PRINT_OUT
      else
        {
          cout << "nb_first_save_table_u=" << nb_first_save_table_u << " nb_table_u=" << nb_table_u << "\n";
        }
#endif
#ifdef PRINT_OUT
      cout << "->table_u=" << table_u << "\n";
#endif
    }
  if(nstacked)
    {
      first_nb_prologue_save_table_y = vertex_count;
      nb_prologue_save_table_y=0;
      First_prologue_table_u = table_u;
      cur_pos=SaveCode.tellp();
      prologue_nb_table_u=nb_table_u;
      i=0;
      SaveCode.write(reinterpret_cast<char *>(&i), sizeof(i));
#ifdef PRINT_OUT
      cout << "--> fixe prologue\n";
#endif
    }
#ifdef PRINT_OUT
#ifdef SIMPLIFY
  cout << "1 free_u_list length0=" << (*((short int*)&(*(((char*)free_u_list) - 4))) & ~(3)) << "\n";
#endif /**SIMPLIFY**/
#endif
#ifdef DEBUGR
  //Check the graph
#ifdef PRINT_OUT
  cout << "going to Check Graph\n";
#endif
  Check_Graph(model_graph);
#endif
  per = 0;
  OK = true;
#ifdef PRINT_OUT
  cout << "nb_endo=" << nb_endo << " y_kmin=" << y_kmin << " y_kmax=" << y_kmax << " size=" << size << "\n";
#endif
  int pos = 0;
  j = 0;
  if(nstacked)
    try_to_interpolate = true;
#ifdef PRINT_OUT
  cout << "try_to_interpolate=" << try_to_interpolate << "\n";
  cout << "y_kmin=" << y_kmin << " y_kmax=" << y_kmax << "\n";
  cout << "(0) nb_table_u=" << nb_table_u << "\n";
#endif
  if(nstacked)
    save_direct=true;
  while(OK)
    {
#ifdef DEBUGR
      Check_Graph(model_graph);
#endif /**DEBUGR**/
      if(!(j % nb_endo))
        length_markowitz = nb_endo;
      else
        length_markowitz--;
      if(nstacked)
        {
#ifdef PRINT_OUT
          cout << "try_to_interpolate=" << try_to_interpolate << "\n";
#endif
          if(try_to_interpolate)
            {
#ifdef PRINT_OUT
              cout << j << " % " << nb_endo << " = " << (j % nb_endo) << "\n";
#endif
              if(!(j % nb_endo))
                {
                  per++;
                  if(per == y_kmin + 1)
                    {
                      nb_prologue_save_table_y=vertex_count-first_nb_prologue_save_table_y;
                      prologue_save_table_y=(t_table_y*)malloc(nb_prologue_save_table_y*sizeof(*prologue_save_table_y));
                      for(i=first_nb_prologue_save_table_y;i<vertex_count;i++)
                        prologue_save_table_y[i-first_nb_prologue_save_table_y]=table_y[i];
                      k=SaveCode.tellp();
                      SaveCode.seekp(cur_pos);
                      i=nb_table_u-prologue_nb_table_u;
                      SaveCode.write(reinterpret_cast<char *>(&i), sizeof(i));
                      SaveCode.seekp(k);
                      if(nstacked)
                        save_direct=false;
                      first_u_blck = table_u;
                      nb_first_u_blck= nb_table_u;
                      nb_first_y_blck= vertex_count;
                      first_y_blck = vertex_count;
#ifdef PRINT_OUT
                      cout << "(1) nb_table_u=" << nb_table_u << "\n";
                      cout << "first_u_blck=" << first_u_blck << "\n";
                      system("pause");
#endif
                    }
                  else if(per == y_kmin + 2)
                    {
                      second_u_blck = table_u;
                      nb_second_u_blck= nb_table_u;
                      second_y_blck = vertex_count;
                      nb_second_y_blck= vertex_count;
#ifdef PRINT_OUT
                      cout << "(2) nb_table_u=" << nb_table_u << "\n";
                      cout << "second_u_blck=" << second_u_blck << "\n";
                      system("pause");
#endif
                    }
                  else if(per >= y_kmin + 3)
                    {
                      third_u_blck = table_u;
                      nb_third_u_blck= nb_table_u;
                      third_y_blck = vertex_count;
                      nb_third_y_blck= vertex_count;
#ifdef PRINT_OUT
                      cout << "(3) nb_table_u=" << nb_table_u << "\n";
                      cout << "third_u_blck=" << third_u_blck << "\n";
                      system("pause");
#endif
                    }
                }
            }
          OK = Vertex_Elimination(model_graph, pos, &interpolate, length_markowitz, dynamic);
        }
      else
        {
#ifdef PRINT_OUT
          cout << "j=" << j << "\n";
#endif
          OK = Vertex_Elimination(model_graph, pos, &interpolate, length_markowitz, dynamic);
        }
#ifdef SIMPLIFY
#endif /**SIMPLIFY**/
      j++;
#ifdef SIMPLIFY
      if(!(j % nb_endo))
        nb_free_u_list = 0;
#endif /**SIMPLIFY**/
      if(nstacked)
        {
#ifdef PRINT_OUT
          cout << "try_to_interpolate=" << try_to_interpolate << "\n";
#endif
          if(try_to_interpolate)
            {
#ifdef PRINT_OUT
              cout << j << " % " << nb_endo << " = " << (j % nb_endo) << "\n";
#endif
              if(!((j) % nb_endo))
                {
#ifdef PRINT_OUT
                  cout << "per (" << per << ") >= y_kmin (" << y_kmin << ")+3\n";
#endif
                  if( per >= y_kmin + 3)
                    {
                      // pctimer_t t1, t2;
#ifdef PRINT_OUT
                      cout << "bef check regularity\n";
                      system("pause");
#endif
                      if(Check_Regularity(first_u_blck, second_u_blck, third_u_blck))
                        {
#ifdef PRINT_OUT
                          cout << "af check regularity OK \n";
                          system("pause");
                          cout << "nb_first_save_table_u=" << nb_first_save_table_u << "\n";
#endif
#ifdef PRINT_OUT
                          cout << "table_u=" << table_u << "  table_u->pNext=" << table_u->pNext << "\n";
                          cout << "(3) nb_table_u=" << nb_table_u << "\n";
#endif
                          // t1 = pctimer();

                          k=SaveCode.tellp();
#ifdef PRINT_OUT
                          cout << "middle_count_loop= " << middle_count_loop << "\n";
#endif
                          i=0;
                          SaveCode.write(reinterpret_cast<char *>(&i), sizeof(i));
                          table_u = interpolation(model_graph, table_y, per_next /*+1*/ + 2, true, per - y_kmin - 3);
#ifdef PRINT_OUT
                          cout << "middle_count_loop= " << middle_count_loop << "\n";
#endif
                          i=SaveCode.tellp();
                          SaveCode.seekp(k);
                          SaveCode.write(reinterpret_cast<char *>(&middle_count_loop), sizeof(middle_count_loop));
                          SaveCode.seekp(i);
                          i=0;
                          SaveCode.write(reinterpret_cast<char *>(&i), sizeof(i));  //last_u
                          OK = false;
                          vertex_count -= nb_endo;
                          // t2 = pctimer();
#ifdef PRINT_OUT
                          // cout << "interpolate=" << 1000*(t2 - t1) << "\n";
                          cout << "table_u=" << table_u << "  table_u->pNext=" << table_u->pNext << "\n";
#endif
#ifdef SAVE
                          last_save_table_u = table_u;
                          nb_last_save_table_y = vertex_count;
#endif /**SAVE**/
                        }
                      else
                        {
                          nb_prologue_save_table_y=nb_second_y_blck-first_nb_prologue_save_table_y;
                          prologue_save_table_y=(t_table_y*)malloc(nb_prologue_save_table_y*sizeof(*prologue_save_table_y));
                          for(i=first_nb_prologue_save_table_y;i<nb_second_y_blck;i++)
                            prologue_save_table_y[i-first_nb_prologue_save_table_y]=table_y[i];
                          k=SaveCode.tellp();
                          SaveCode.seekp(cur_pos);
                          i=nb_second_u_blck-prologue_nb_table_u+1;
                          //cout << "nb_prologue_save_table_u=" << i << "\n";
                          SaveCode.write(reinterpret_cast<char *>(&i), sizeof(i));
                          SaveCode.seekp(k);
                          write_to_file_table_u_b(first_u_blck,second_u_blck, &nb_prologue_save_table_u,true);
                          //cout << "nb_prologue_save_table_u(1)=" << nb_prologue_save_table_u << "\n";
                          nb_first_u_blck=nb_second_u_blck;
                          nb_second_u_blck=nb_third_u_blck+1;
                          nb_first_y_blck=nb_second_y_blck;
                          nb_second_y_blck=nb_third_y_blck;
                          first_u_blck = second_u_blck;
                          second_u_blck = third_u_blck;
                          first_y_blck = second_y_blck;
                          second_y_blck = third_y_blck;
                        }
                    }
                }
#ifdef DEBUGR
              cout << 100*j / size << "%\n";
#endif /**DEBUGR**/
            }
        }
    }
#ifdef PRINT_OUT
  cout << "nstacked=" << nstacked << "\n";
  cout << "nb_first_save_table_y=" << nb_first_save_table_y << "\n";
  cout << "nb_prologue_save_table_y=" << nb_prologue_save_table_y << "\n";
  cout << "nb_middle_save_table_y=" << nb_middle_save_table_y << "\n";
  cout << "nb_last_save_table_y=" << nb_last_save_table_y << "\n";
#endif
  if((nstacked)&&(!nb_last_save_table_y))
    {
      cout << "not synchronized per=" << per << "\n";
      exit(-1);
    }

#ifdef PRINT_OUT
  cout << "end vertex supress\n";
#endif
#ifdef SAVE
  if(!nstacked)
    table_u->pNext = NULL;
  if(nstacked)
    {
      table_u = NULL;
      last_save_table_u = NULL;
      nb_last_save_table_y = 0;
      OK = true;
    }
  else
    {
      SaveCode.write(reinterpret_cast<char *>(&nb_table_u), sizeof(nb_table_u));
      write_to_file_table_u_b(First_table_u->pNext, table_u->pNext, &nb_last_save_table_u, false );
      nb_last_save_table_y = vertex_count;
      last_save_table_y=(t_table_y*)malloc(nb_last_save_table_y*sizeof(*last_save_table_y));
      for(i=0;i<nb_last_save_table_y;i++)
        last_save_table_y[i]=table_y[i];
    }
#ifdef PRINT_OUT
  cout << "goint to write nb_endo=" << nb_endo << "\n";
#endif
#ifdef PRINT_OUT
  cout << "-->write prologue\n";
#endif
#ifdef PRINT_OUT
  cout << "-->write first\n";
#endif
#ifdef PRINT_OUT
  cout << "middle_count_loop=" << middle_count_loop << "\n";
  cout << "-->write middle\n";
#endif
#ifdef PRINT_OUT
  cout << "-->write last\n";
#endif
  nb_last_save_table_u = i;
#ifdef PRINT_OUT
  cout << "nb_prologue_save_table_y=" << nb_prologue_save_table_y << "\n";
#endif
  write_to_file_table_y( prologue_save_table_y, NULL, nb_prologue_save_table_y, 0);
#ifdef PRINT_OUT
  cout << "nb_first_save_table_y=" << nb_first_save_table_y << "\n";
#endif
  write_to_file_table_y( first_save_table_y, NULL, nb_first_save_table_y, 0);
#ifdef PRINT_OUT
  cout << "nb_middle_save_table_y=" << nb_first_save_table_y << "\n";
#endif
  write_to_file_table_y( middle_save_table_y, middle_save_i_table_y, nb_middle_save_table_y, nb_middle_save_table_y);
#ifdef PRINT_OUT
  cout << "//// last_save_table_y ===\n";
  cout << "nb_last_save_table_y=" << nb_last_save_table_y << "\n";
#endif
  write_to_file_table_y( last_save_table_y, NULL, nb_last_save_table_y, 0);
  SaveCode.close();
#endif /**SAVE**/
#ifdef SIMPLIFY
#ifdef PRINT_OUT
  cout << "try_to_interpolate=" << try_to_interpolate << "\n";
#endif
  free(free_u_list);
  free(free_u_list1);
#endif /**SIMPLIFY**/
  free(s_i1);
  free(s_i2);
  free(s_j2);
}

void
SymbolicGaussElimination::file_is_open()
{
  file_open=true;
  //file_is_open1();
}

void
SymbolicGaussElimination::SGE_compute(Model_Block *ModelBlock, int blck, bool dynamic, string file_name, int endo_nbr)
{
  t_model_graph *model_graph;
  int block_u_count, nb_table_y;
  // pctimer_t t1, t2;
  int i;
  int mean_var_in_equation;

  init_glb();
  model_graph = (t_model_graph*)malloc(sizeof(*model_graph));
  nstacked = dynamic;
#ifdef PRINT_OUT
  periods = ModelBlock->Periods;
  cout << "nstacked=" << nstacked << "\n";
  // t1 = pctimer();
  cout << "periods=" << periods << "\n";
#endif

  u_count = ModelBlock_Graph(ModelBlock, blck, dynamic, model_graph, endo_nbr, &block_u_count, &starting_vertex, &periods, &nb_table_y, &mean_var_in_equation);
  if(dynamic)
    {
      cout << "Mean endogenous variable per equation: " << mean_var_in_equation << ", density indicator=" << double(mean_var_in_equation)/endo_nbr*100 << "%\n";
      cout << "Coding the model...";
    }
  Time = periods;
#ifdef PRINT_OUT
  // t2 = pctimer();
  // cout << "/* ModelBlock_Graph : " << 1000*(t2 - t1) << " milliseconds u_count : " << u_count << "*/\n";
#endif
  int size = model_graph->nb_vertices;
  nb_endo = ModelBlock->Block_List[blck].Size;
  y_kmin = ModelBlock->Block_List[blck].Max_Lag;
  y_kmax = ModelBlock->Block_List[blck].Max_Lead;
  periods = ModelBlock->Periods;
  if(periods<=3)
    nstacked=false;
  //cout << "periods=" << periods << " y_kmin=" << y_kmin << " y_kmax=" << y_kmax << "\n";
  u_count_init = u_count;
#ifdef PRINT_OUT
  cout << "size : " << size << "\n";
#endif
  table_y = (t_table_y*)malloc((size + 1) * sizeof(t_table_y));
  for(i = 0;i <= size;i++)
    {
      table_y[i].u_index = (int*)malloc(nb_table_y * sizeof(int));
      table_y[i].y_index = (int*)malloc(nb_table_y * sizeof(int));
      table_y[i].index = i;
      table_y[i].nb = 0;
    }
  table_u = (t_table_u*)malloc(sizeof(*table_u)-3*sizeof(int));
  table_u->pNext = NULL;
  First_table_u = table_u;
#ifdef PRINT_OUT
  cout << "où en est-on ? nb_table_y=" << nb_table_y << "\n";
  system("pause");
  cout << "u_count_init=" << u_count_init << "\n";
  cout << "table_u=" << table_u << "\n";
  t1 = pctimer();
#endif
  Gaussian_Elimination(model_graph, file_name, dynamic);
  free_model_graph(model_graph);
#ifdef PRINT_OUT
  cout << "max_nb_table_y=" << max_nb_table_y << "\n";
  cout << "max_nb_in_degree_edges=" << max_nb_in_degree_edges << "\n";
  cout << "max_nb_out_degree_edges=" << max_nb_out_degree_edges << "\n";
  t2 = pctimer();
  cout << "/* Gaussian Elimination : " << 1000*(t2 - t1) << " milliseconds u_count : " << u_count << "*/\n";
#endif
}

