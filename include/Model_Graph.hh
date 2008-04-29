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

#ifndef MODEL_GRAPH
#define MODEL_GRAPH
#define DIRECT_COMPUTE
#define SORTED
#define SIMPLIFY
#define SIMPLIFYS
#define SAVE
#define COMPUTE
//#define PRINT_OUT_OUT
//#define PRINT_OUT_1
#define DIRECT_SAVE
#include "ModelTree.hh"
#include "BlockTriangular.hh"

typedef struct t_edge
{
  int index, u_count;
};

typedef struct t_vertex
{
  t_edge *out_degree_edge, *in_degree_edge;
  int nb_out_degree_edges, nb_in_degree_edges;
  int max_nb_out_degree_edges, max_nb_in_degree_edges;
  int index, lag_lead;
};

typedef struct t_model_graph
{
  int nb_vertices;
  t_vertex* vertex;
};

void free_model_graph(t_model_graph* model_graph);
void print_Graph(t_model_graph* model_graph);
void Check_Graph(t_model_graph* model_graph);
int ModelBlock_Graph(Model_Block *ModelBlock, int Blck_num,bool dynamic, t_model_graph* model_graph, int nb_endo, int *block_u_count, int *starting_vertex, int* periods, int *nb_table_y, int *mean_var_in_equ);
#endif
