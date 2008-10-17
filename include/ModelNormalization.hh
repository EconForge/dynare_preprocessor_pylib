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

#ifndef MODELNORMALIZATION
#define MODELNORMALIZATION
#include "SymbolTable.hh"
#include "CodeInterpreter.hh"


struct Edge
{
  Edge *next;
  int Vertex_Index;
};

struct Equation_vertex
{
  Edge *First_Edge;
  Edge *Next_Edge;
  int matched;
};

struct Equation_set
{
  Equation_vertex *Number;
  int size;
  int edges;
};

//! Stores result of block decomposition for ONE equation or ONE variable
struct simple
{
  //! New {variable, equation} index after reordering
  int index;
};

class Normalization
{
private:
  struct Variable_vertex
  {
    int  matched;
  };
  struct Variable_set
  {
    Variable_vertex *Number;
    int size;
  };
  struct t_Heap
  {
    int u;        /* vertex */
    int i_parent; /* index in t_Heap of parent vertex in tree of u */
    int v;        /* current matched of u */
  };
public:
  Normalization(const SymbolTable &symbol_table_arg);
  bool Normalize(int n, int prologue, int epilogue, bool* IM, simple* Index_Var_IM, Equation_set* Equation,bool mixing, bool* IM_s);
  void Gr_to_IM_basic(int n0, int prologue, int epilogue, bool* IM, Equation_set *Equation,bool transpose);
  const SymbolTable &symbol_table;
  void Set_fp_verbose(bool ok);
  void Free_Equation(int n, Equation_set* Equation);
private:
  void IM_to_Gr(int n0, int prologue, int epilogue, bool* IM, Equation_set *Equation, Variable_set *Variable );
  void Inits(Equation_set *Equation);
  void UpdatePath(Equation_set *Equation, Variable_set *Variable, int i1, int i2);
  void FindAugmentingPaths(Equation_set *Equation, Variable_set *Variable);
  void CheapMatching(Equation_set *Equation, Variable_set *Variable);
  void MaximumMatching(Equation_set *Equation, Variable_set *Variable);
  int MeasureMatching(Equation_set *Equation);
  void OutputMatching(Equation_set* Equation);
  void Gr_to_IM(int n0, int prologue, int epilogue, bool* IM, simple* Index_Var_IM, Equation_set *Equation,bool mixing, bool* IM_s);
  void Free_Other(Variable_set* Variable);
  void Free_All(int n, Equation_set* Equation, Variable_set* Variable);
  int eq, eex;
  int IndexUnmatched;
  bool fp_verbose;
  bool* visited;
  t_Heap* Local_Heap;
};
#endif
