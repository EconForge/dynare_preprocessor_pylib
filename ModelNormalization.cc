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

//#define DEBUG
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include "ModelNormalization.hh"


using namespace std;


Normalization::Normalization(const SymbolTable &symbol_table_arg) :
  symbol_table(symbol_table_arg), fp_verbose(false)
{
}

void
Normalization::IM_to_Gr(int n0, int prologue, int epilogue, bool* IM, Equation_set *Equation, Variable_set *Variable )
// Create a non-oriented graph of the model from the incidence matrix
{
  int i, j, edges, n;
  Edge *e1;
#ifdef DEBUG
  cout << "in IM_to_Gr\n";
#endif
  //Normalize only the earth block (the prologue and the epilogue are still normalized)
  n = n0 - prologue - epilogue;
  Equation->size = n;
  Variable->size = n;
  Equation->Number = (Equation_vertex*)malloc(n * sizeof(Equation_vertex));
  Variable->Number = (Variable_vertex*)malloc(n * sizeof(Variable_vertex));
  edges = 0;
  for(i = 0;i < n;i++)
    {
      Equation->Number[i].First_Edge = NULL;
      Equation->Number[i].matched = -1;
      Variable->Number[i].matched = -1;
      for(j = 0;j < n;j++)
        {
          if(IM[(j + prologue)*n0 + (i + prologue)])
            {
              edges++;
              e1 = (Edge *) malloc(sizeof(Edge));
              e1->next = Equation->Number[i].First_Edge;
              Equation->Number[i].First_Edge = e1;
              e1->Vertex_Index = j;
            }
        }
    }
  //The maximum number of vertex in each equation is set to the total amount of edges in the model
  Equation->edges = edges;
#ifdef DEBUG
  cout << "end of IM_to_Gr\n";
#endif
}

void
Normalization::Inits(Equation_set *Equation)
{
  int i;
#ifdef DEBUG
  cout << "in Inits\n";
#endif
  eq = eex = 0;
  IndexUnmatched = Equation->edges * 2;
  Local_Heap = (t_Heap*)malloc(IndexUnmatched * sizeof(t_Heap));
  for(i = 0; i < Equation->size; i++)
    {
      Equation->Number[i].Next_Edge = Equation->Number[i].First_Edge;
      visited[i] = 0;
      // we put all unmatched vertices from Equation at the other end of the Local_Heap
      if(Equation->Number[i].matched == -1)
        {
          Local_Heap[--IndexUnmatched].u = i;
#ifdef DEBUG
          cout << i << " is unmatched\n";
#endif
        }
    }
#ifdef DEBUG
  cout << "end of Inits\n";
#endif
}

void
Normalization::UpdatePath(Equation_set *Equation, Variable_set *Variable, int i1, int i2)
{
  int i, j;
#ifdef DEBUG
  cout << "in UpdatePath \n";
#endif
  while(i2 >= 0)
    {
      i = Local_Heap[i2].u;
      j = Local_Heap[i1].v;
      Variable->Number[j].matched = i;
      Equation->Number[i].matched = j;
      i1 = i2;
      i2 = Local_Heap[i2].i_parent;
      eex++;
    }
#ifdef DEBUG
  cout << "end of UpdatePath \n";
#endif
}

void
Normalization::FindAugmentingPaths(Equation_set *Equation, Variable_set *Variable)
{
  // augmenting paths using breadth-first search.
  int Bottom;
  int Top;
  int u, i;
  Edge *e, *e2;
#ifdef DEBUG
  cout << "in FindAugmentingPaths\n";
#endif
  // external loop gets unmatched u vertices from far end of array Local_Heap
  while(IndexUnmatched < Equation->edges*2)
    {
      Top = Bottom = 0;
      Local_Heap[Top].u = Local_Heap[IndexUnmatched++].u;
      Local_Heap[Top].i_parent = -1; /* root of BFS tree */
#ifdef DEBUG
      cout << "unmatched u" << Local_Heap[Top].u << " will be processed\n";
#endif
      // Local_Heap processing
      while(Bottom >= Top)
        {
          u = Local_Heap[Top++].u;
          e = Equation->Number[u].First_Edge;
          eq++;
          // adjacency list scanning
          while(e != NULL)
            {
              if (!visited[Variable->Number[e->Vertex_Index].matched])
                {
                  // extend tree
                  Local_Heap[++Bottom].u = u = Variable->Number[e->Vertex_Index].matched;
                  Local_Heap[Bottom].i_parent = Top - 1;
                  Local_Heap[Bottom].v = e->Vertex_Index;
                  visited[u] = 1;
                  e2 = Equation->Number[u].Next_Edge;
                  eq++;
                  while ((e2 != NULL) && (Variable->Number[e2->Vertex_Index].matched != -1))
                    {
                      e2 = e2->next;
                      eq++;
                    }
                  Equation->Number[u].Next_Edge = e2;
                  if(e2 != NULL)
                    {
#ifdef DEBUG
                      cout << "augmenting path found\n";
#endif
                      // u in the Local_Heap but not the edge to v
                      Variable->Number[e2->Vertex_Index].matched = u;
                      Equation->Number[u].matched = e2->Vertex_Index;
                      // now for the rest of the path
                      UpdatePath(Equation, Variable, Bottom, Top - 1);
                      // temporary cut is emptied
                      for(i = 0; i <= Bottom; i++)
                        visited[Local_Heap[i].u] = 0;
                      Bottom = Top - 1;
                      // to get off from Local_Heap loop
                      // to get off from adj list scan loop
                      break;
                    }
                }
              e = e->next;
              eq++;
            }
        }
    }
#ifdef DEBUG
  cout << "end of FindAugmentingPaths\n";
#endif
}


void
Normalization::CheapMatching(Equation_set *Equation, Variable_set *Variable)
{
  int i;
  Edge *e;
  int count = 0;
#ifdef DEBUG
  cout << "in CheapMatching Equation->size : " << Equation->size << "\n";
#endif
  for(i = 0; i < Equation->size; i++)
    {
      e = Equation->Number[i].First_Edge;
      while(e != (Edge *) NULL)
        {
          if(Variable->Number[e->Vertex_Index].matched == -1)
            {
              Variable->Number[e->Vertex_Index].matched = i;
              Equation->Number[i].matched = e->Vertex_Index;
#ifdef DEBUG
              cout << i << " matched to " << e->Vertex_Index << "\n";
#endif
              count++;
              break;
            }
          e = e->next;
        }
    }
  if(fp_verbose)
    cout << count << " vertices in Equation were initially matched (" << (float) 100*count / Equation->size << "%)\n";
#ifdef DEBUG
  cout << "end of CheapMatching\n";
#endif
}


void
Normalization::MaximumMatching(Equation_set *Equation, Variable_set *Variable)
{
#ifdef DEBUG
  cout << "in MaximumMatching\n";
#endif
  CheapMatching(Equation, Variable);
  Inits(Equation);
  FindAugmentingPaths(Equation, Variable);
#ifdef DEBUG
  cout << "end of MaximumMatching\n";
#endif
}

int
Normalization::MeasureMatching(Equation_set *Equation)
{
  int size = 0, i;
  for(i = 0; i < Equation->size; i++)
    if(Equation->Number[i].matched != -1)
      size++;
  return size;
}

void
Normalization::OutputMatching(Equation_set* Equation)
{
  int i;
  Edge* e1;
  cout << "Maximum Matching Results for |Equation|=" << Equation->size << " |Edges|=" << Equation->edges << "\n";
  for(i = 0; i < Equation->size; i++)
    {
      if(Equation->Number[i].matched != -1)
        cout << "equation " << i << " matched to variable " << Equation->Number[i].matched;
      else
        cout << "equation " << i << " not matched \n";
      e1 = Equation->Number[i].First_Edge;
      while(e1 != NULL)
        {
          cout << " " << e1->Vertex_Index;
          e1 = e1->next;
        }
      cout << "\n";
    }
}

void
Normalization::Gr_to_IM_basic(int n0, int prologue, int epilogue, bool* IM, Equation_set *Equation, bool transpose)
{
  int i, j, edges, n;
  Edge *e1, *e2;
  n = n0 - prologue - epilogue;
  Equation->size = n;
  if(Equation->Number)
    {
      for(i = 0;i < n;i++)
        {
          e1 = Equation->Number[i].First_Edge;
          while(e1 != NULL)
            {
              e2 = e1->next;
              free(e1);
              e1 = e2;
            }
        }
      free(Equation->Number);
    }
  Equation->Number = (Equation_vertex*)malloc(n * sizeof(Equation_vertex));
  edges = 0;
  if(transpose)
    {
      for(i = 0;i < n;i++)
        {
          Equation->Number[i].First_Edge = NULL;
          Equation->Number[i].matched = -1;
          for(j = 0;j < n;j++)
            {
              if ((IM[(j + prologue)*n0 + (i + prologue)]) && (i != j))
                {
                  edges++;
                  e1 = (Edge *) malloc(sizeof(Edge));
                  e1->next = Equation->Number[i].First_Edge;
                  Equation->Number[i].First_Edge = e1;
                  e1->Vertex_Index = j;
                }
            }
        }
    }
  else
    {
      for(i = 0;i < n;i++)
        {
          Equation->Number[i].First_Edge = NULL;
          Equation->Number[i].matched = -1;
          for(j = 0;j < n;j++)
            {
              if ((IM[(i + prologue)*n0 + (j + prologue)]) && (i != j))
                {
                  edges++;
                  e1 = (Edge *) malloc(sizeof(Edge));
                  e1->next = Equation->Number[i].First_Edge;
                  Equation->Number[i].First_Edge = e1;
                  e1->Vertex_Index = j;
                }
            }
        }
    }
  //The maximum number of vertex in each equation is set to the total amount of edges in the model
  Equation->edges = edges;
#ifdef DEBUG
  cout << "end of IM_to_Gr\n";
#endif
}

void
Normalization::Gr_to_IM(int n0, int prologue, int epilogue, bool* IM, simple* Index_Equ_IM, Equation_set *Equation, bool mixing, bool* IM_s)
{
  int i, j, n, l;
  Edge *e1, *e2;
  Equation_set* Equation_p;
  simple* Index_Equ_IM_tmp = (simple*)malloc(n0 * sizeof(*Index_Equ_IM_tmp));
  bool* SIM = (bool*)malloc(n0 * n0 * sizeof(bool));
#ifdef DEBUG
  cout << "in Gr_to_IM\n";
#endif
  n = n0 - prologue - epilogue;
  if(mixing)
    {
      for(i = 0;i < n0*n0;i++)
        SIM[i] = IM_s[i];
      for(i = 0;i < n0;i++)
        Index_Equ_IM_tmp[i].index = Index_Equ_IM[i].index;
      for(i = 0;i < n;i++)
        {
          /*Index_Var_IM[j+prologue].index=Index_Var_IM_tmp[Equation->Number[j].matched+prologue].index;*/
          if(fp_verbose)
            cout << "Equation->Number[" << i << "].matched=" << Equation->Number[i].matched << "\n";
          Index_Equ_IM[i + prologue].index = Index_Equ_IM_tmp[Equation->Number[i].matched + prologue].index;
          for(j = 0;j < n0;j++)
            SIM[(i + prologue)*n0 + j] = IM_s[(Equation->Number[i].matched + prologue) * n0 + j];
        }
      for(i = 0;i < n0*n0;i++)
        IM[i] = SIM[i];
    }
  else
    {
      for(i = 0;i < n0*n0;i++)
        SIM[i] = IM[i];
      for(i = 0;i < n0;i++)
        Index_Equ_IM_tmp[i].index = Index_Equ_IM[i].index;
      for(j = 0;j < n;j++)
        {
          if(fp_verbose)
            cout << "Equation->Number[" << j << "].matched=" << Equation->Number[j].matched << "\n";
          Index_Equ_IM[j + prologue].index = Index_Equ_IM_tmp[Equation->Number[j].matched + prologue].index;
          for(i = 0;i < n0;i++)
            SIM[(i)*n0 + j + prologue] = IM[(i) * n0 + Equation->Number[j].matched + prologue];
        }
      for(i = 0;i < n0*n0;i++)
        IM[i] = SIM[i];
    }
  free(SIM);
  free(Index_Equ_IM_tmp);
  //cout << "mixing=" << mixing << "\n";
  if(mixing)
    {
      //Free_Equation(n,Equation);
      Gr_to_IM_basic(n0, prologue, epilogue, IM, Equation, true);
    }
  else
    {
      //  In this step we :
      //  1) get ride of the edge from the equation to its explain variable
      //  2) resort the equation in the order of the matched variable
      //  3) transpose the graph
      //  in order to get the oriented graph needed to find strong connex components
      Equation_p = (Equation_set*)malloc(sizeof(Equation_set));
      Equation_p->size = Equation->size;
      Equation_p->edges = Equation->edges;
      Equation_p->Number = (Equation_vertex*)malloc(n * sizeof(Equation_vertex));
      for(i = 0;i < n;i++)
        {
          Equation_p->Number[i].First_Edge = NULL;
          Equation_p->Number[i].Next_Edge = NULL;
        }
      for(i = 0;i < n;i++)
        {
          l = Equation->Number[i].matched;
          e1 = Equation->Number[l].First_Edge;
          while(e1 != NULL)
            {
              if(e1->Vertex_Index != i)
                {
                  j = e1->Vertex_Index;
                  if(Equation_p->Number[j].First_Edge != NULL)
                    {
                      Equation_p->Number[j].Next_Edge->next = (Edge*)malloc(sizeof(Edge*));
                      Equation_p->Number[j].Next_Edge = Equation_p->Number[j].Next_Edge->next;
                    }
                  else
                    {
                      Equation_p->Number[j].First_Edge = (Edge*)malloc(sizeof(Edge*));
                      Equation_p->Number[j].Next_Edge = Equation_p->Number[j].First_Edge;
                    }
                  Equation_p->Number[j].Next_Edge->next = NULL;
                  Equation_p->Number[j].Next_Edge->Vertex_Index = i;
                }
              e2 = e1->next;
              free(e1);
              e1 = e2;
            }
        }
      for(i = 0;i < n;i++)
        {
          Equation->Number[i].matched = Equation_p->Number[i].matched;
          Equation->Number[i].First_Edge = Equation_p->Number[i].First_Edge;
          Equation->Number[i].Next_Edge = Equation_p->Number[i].Next_Edge;
        }
      free(Equation_p->Number);
      free(Equation_p);
    }
#ifdef DEBUG
  cout << "end of Gr_to_IM\n";
#endif
}

void
Normalization::Free_Equation(int n, Equation_set* Equation)
{
  //free unused space
  Edge *e1, *e2;
  int i;
  for(i = 0;i < n;i++)
    {
      e1 = Equation->Number[i].First_Edge;
      while(e1 != NULL)
        {
          e2 = e1->next;
          free(e1);
          e1 = e2;
        }
    }
  free(Equation->Number);
  //free(Equation);
}

void
Normalization::Free_Other(Variable_set* Variable)
{
  //free unused space
#ifdef DEBUG
  cout << "Free_Other\n";
#endif
  free(Local_Heap);
  free(Variable->Number);
  free(Variable);
  free(visited);
}

void
Normalization::Free_All(int n, Equation_set* Equation, Variable_set* Variable)
{
  Free_Equation(n, Equation);
  Free_Other(Variable);
}

void
Normalization::Set_fp_verbose(bool ok)
{
  fp_verbose=ok;
}

bool
Normalization::Normalize(int n, int prologue, int epilogue, bool* IM, simple* Index_Equ_IM, Equation_set* Equation, bool mixing, bool* IM_s)
{
  int matchingSize, effective_n;
  int save_fp_verbose=fp_verbose;
  fp_verbose = 0;
  Variable_set* Variable = (Variable_set*) malloc(sizeof(Variable_set));
#ifdef DEBUG
  cout << "in Normalize\n";
#endif
  visited = (bool*)malloc(n * sizeof(*visited));
  IM_to_Gr(n, prologue, epilogue, IM, Equation, Variable);
  MaximumMatching(Equation, Variable);
  matchingSize = MeasureMatching(Equation);
  effective_n = n - prologue - epilogue;
  fp_verbose=save_fp_verbose;
  if(matchingSize < effective_n && fp_verbose)
    {
      cout << "Error: dynare could not normalize the model.\n The following equations:\n - ";
      int i;
      for(i = 0; i < Equation->size; i++)
        if(Equation->Number[i].matched == -1)
          cout << i << " ";
      cout << "\n and the following variables:\n - ";
      for(i = 0; i < Variable->size; i++)
        if(Variable->Number[i].matched == -1)
          cout << symbol_table.getName(Index_Equ_IM[i].index) << " ";
      cout << "\n could not be normalized\n";
      //ErrorHandling(n, IM, Index_Equ_IM);
      //system("PAUSE");
      exit(EXIT_FAILURE);
    }
  if(matchingSize >= effective_n )
    {
      Gr_to_IM(n, prologue, epilogue, IM, Index_Equ_IM, Equation, mixing, IM_s);
      if(fp_verbose)
        {
          OutputMatching(Equation);
          for(int i = 0;i < n;i++)
            cout << "Index_Equ_IM[" << i << "]=" << Index_Equ_IM[i].index /*<< " == " "Index_Var_IM[" << i << "]=" << Index_Var_IM[i].index*/ << "\n";
        }
    }
  Free_Other(Variable);
  //Free_All(n,Equation,Variable);
#ifdef DEBUG
  cout << "end of Normalize\n";
#endif
  if(matchingSize < effective_n )
    return(0);
  else
    return(1);
}

