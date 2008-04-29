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

#ifndef MODELBLOCKS
#define MODELBLOCKS
#include "ModelNormalization.hh"


typedef struct block_result
{
  int size, n_sets;
  int *vertices;
  int *sets_s, *sets_f;
  int *order, *ordered;
}
  block_result_t;



class Blocks
{
public:
  Blocks();
  void block_depth_search(int v);
  block_result_t* sc(Equation_set *g);
  void block_result_free(block_result_t *r);
  void block_result_print(block_result_t *r);
  void Print_Equation_gr(Equation_set* Equation);
  void block_result_to_IM(block_result_t *r,bool* IM,int prologue, int n,simple* Index_Equ_IM,simple* Index_Var_IM);
  Equation_vertex *vertices;
  int *block_vertices, *sets_s, *sets_f;
  int *block_stack, *sp, tos;
  int *visit_nos, *low_link_nos;
  int n_visited, n_written, n_sets;
  int *pos_sc;
};
#endif

