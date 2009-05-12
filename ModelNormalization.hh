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

#ifndef _MODELNORMALIZATION_HH
#define _MODELNORMALIZATION_HH
#include "SymbolTable.hh"
#include "CodeInterpreter.hh"


//! One edge in the bi-partite graph (equation side), stored in a chained-list
struct Edge
{
  Edge *next;
  int Vertex_Index; //!< Variable linked to the equation
};

//! Set of the edges going to a given equation
struct Equation_vertex
{
  Edge *First_Edge;
  Edge *Next_Edge;
  int matched;
};

//! Bi-partite graph, seen from the equation side
struct Equation_set
{
  Equation_vertex *Number;
  int size;
  int edges;
};


//! Computes the model normalization
class Normalization
{
private:
  //! Indicates if a given variable vertex is matched
  struct Variable_vertex
  {
    int  matched;
  };
  //! Data extracted from the bi-partite graph, seen from the variable side
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
  bool Normalize(int n, int prologue, int epilogue, bool* IM, vector<int> &Index_Var_IM, Equation_set* Equation,bool mixing, bool* IM_s);
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
  void Gr_to_IM(int n0, int prologue, int epilogue, bool* IM, vector<int> &Index_Var_IM, Equation_set *Equation,bool mixing, bool* IM_s);
  void Free_Other(Variable_set* Variable);
  void Free_All(int n, Equation_set* Equation, Variable_set* Variable);
  int eq, eex;
  int IndexUnmatched;
  bool fp_verbose;
  bool* visited;
  t_Heap* Local_Heap;
};
#endif
