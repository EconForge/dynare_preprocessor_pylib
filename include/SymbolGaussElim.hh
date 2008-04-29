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

#ifndef SYMBGAUSSELIM
#define SYMBGAUSSELIM

#include <iostream>
#include <fstream>
#include "Model_Graph.hh"

using namespace std;
#define TOL 1e-9

typedef struct t_table_y
{
  int index,nb;
  int *u_index, *y_index;
};

typedef struct t_table_u
{
  t_table_u* pNext;
  unsigned char type;
  int index;
  int op1,op2;

};

class SymbolicGaussElimination
{
public:
  int y_kmin, y_kmax, nb_endo, Time, u_count, periods;
  int *s_i1, *s_i2, *s_j2;
  bool nstacked, save_direct;
  int nb_loop_table;
  int *loop_table_u_count, *loop_table_vertex_index;
#ifdef DIRECT_COMPUTE
  int u_count_init;
  double tol;
  t_table_y *table_y;
  int vertex_count;
  t_table_u* table_u;
  int nb_table_u;
  t_table_u* First_table_u;
#endif /**DIRECT_COMPUTE**/
  int starting_vertex;
  t_table_u *first_u_blck, *second_u_blck, *third_u_blck;
  int first_y_blck, second_y_blck, third_y_blck;
#ifdef SAVE
  t_table_u *prologue_save_table_u, *first_save_table_u, *first_save_i_table_u, *middle_save_table_u, *middle_save_i_table_u, *last_save_table_u, *save_table_u;
  t_table_y *prologue_save_table_y, *first_save_table_y, *first_save_i_table_y, *middle_save_table_y, *middle_save_i_table_y, *last_save_table_y;
  int nb_prologue_save_table_y, nb_first_save_table_y, nb_middle_save_table_y, nb_last_save_table_y;
  int nb_prologue_save_table_u, nb_first_save_table_u, nb_middle_save_table_u, nb_last_save_table_u, middle_count_loop;
  int pos_nb_first_save_table_u;
  std::ofstream SaveCode;
  bool file_open;
#endif /**SAVE**/
#ifdef SIMPLIFY
  int* free_u_list;
  int* free_u_list1;
  int MAX_FREE_U_LIST;
  int nb_free_u_list,nb_free_u_list1,max_nb_free_u_list1,max_nb_free_u_list;
  bool simplification_allowed;
  void set_free_u_list(int index);
  int get_free_u_list(bool dynamic);
#endif /**SIMPLIFY**/
  SymbolicGaussElimination();
  void init_glb();
  void write_to_file_table_y( t_table_y *save_table_y, t_table_y *save_i_table_y, int nb_save_table_y, int nb_save_i_table_y);
  void write_to_file_table_u_b(t_table_u *save_table_u, t_table_u *last_table_u, int *nb_save_table_u, bool chk);
  bool Check_Regularity(t_table_u *first_u_blck, t_table_u *second_u_blck, t_table_u *third_u_blck);
  t_table_u* interpolation(t_model_graph* model_graph,t_table_y* table_y, int to_add, bool middle,int per);
  bool Loop_Elimination(t_model_graph* model_graph);
  bool Vertex_Elimination(t_model_graph* model_graph, int pos,bool* interpolate, int length_markowitz, bool dynamic);
  void Gaussian_Elimination(t_model_graph* model_graph
#ifdef SAVE
                            , string file_name
#endif
                            , bool dynamic);
  void SGE_compute(Model_Block *ModelBlock, int blck, bool dynamic, string file_name, int endo_nbr);
  void file_is_open();
};

#endif
