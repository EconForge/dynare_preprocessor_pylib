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
#include <sstream>
#include <fstream>
#include "ModelBlocks.hh"

using namespace std;

#define UNDEFINED -1

Blocks::Blocks()
{
  //Empty
}

Blocks::~Blocks()
{
  //Empty
}

int n_sc_set=0;

void
Blocks::block_depth_search(int v)
// block_depth_search()
// find the strong components of the graph using a recursive depth first search
// The results are stored in the global variables block_vertices, sets_s, sets_f.
{
  Edge *edge_ptr;
  int w;
  // Add the vertex v to the visited vertex and store it in the result (low_link_nos)
  // and increase the number of visited vertex
  low_link_nos[v] = visit_nos[v] = n_visited;
  n_visited++;
  // Put v in the stack.
  block_stack[tos] = v;
  sp[v] = tos;
  tos++;
  // Going to visite the edges from vertex v starting
  //   from the First edge of vetexe v
  edge_ptr = vertices[v].First_Edge;
  // While there is edge
  while(edge_ptr)
    {
      w = edge_ptr->Vertex_Index;
      // if the vertex w hasen't been visited
      if(visit_nos[w] == UNDEFINED)
        {
          // visits the vertex w
          block_depth_search(w);
          // Update low_link no.
          if(low_link_nos[w] < low_link_nos[v])
            low_link_nos[v] = low_link_nos[w];
        }
      else if(visit_nos[w] < visit_nos[v] && sp[w] != UNDEFINED)
        {
          // Update low_link no. */
          if(visit_nos[w] < low_link_nos[v])
            if(visit_nos[w]>=0)
              low_link_nos[v] = visit_nos[w];
            else
              {
                // Check for hierarchic structure accross strong connex components
                if(pos_sc[-(visit_nos[w]+2)]<pos_sc[n_sets])
                  {
                    int j=pos_sc[-(visit_nos[w]+2)];
                    pos_sc[-(visit_nos[w]+2)]=pos_sc[n_sets];
                    pos_sc[n_sets]=j;
                  }
              }
        }
      edge_ptr = edge_ptr->next;
    }
  // If all vertices in v's SC component have been found.
  if(low_link_nos[v] == visit_nos[v])
    {
      int vpos = sp[v];
      int i;
      sets_s[n_sets] = n_written;
      // The SC component vertices are stored from the top of the stack, down
      // to v.  Write these to the result structure.
      for(i = vpos; i < tos; i++)
        {
          block_vertices[n_written] = block_stack[i];
          n_written++;
        }
      if(n_sc_set>0)
        for(i=0;i<n_sc_set;i++)
          if(pos_sc[i]<pos_sc[n_sc_set])
            {
              int j=pos_sc[i];
              pos_sc[i]=pos_sc[n_sc_set];
              pos_sc[n_sc_set]=j;
              for(j=sets_s[i];j<=sets_f[i];j++)
                visit_nos[block_vertices[j]] = -(2+pos_sc[i]);
            }
      n_sc_set++;
      for(i = vpos; i < tos; i++)
        {
          visit_nos[block_stack[i]] = -(2+pos_sc[n_sets]);
        }
      // Now remove these vertices from the stack.
      for(i = vpos; i < tos; i++)
        block_stack[i] = UNDEFINED;
      tos = vpos;
      sets_f[n_sets] = n_written - 1;
      n_sets++;
    }
  //  stsz.erase(stsz.length()-1,1);
}



block_result_t*
Blocks::sc(Equation_set *g)
// Generates SC components using Tarjan's algorithm.
// The result is returned as a pointer to a block_result_t structure.  The SC
// components are stored as two arrays:
//   - sets_s[i] gives the starting position in the vertices[] array of SC
//     component i.
//   - sets_f[i] gives the finishing position in the vertices[] array of SC
//     component i.
//   - vertices[] is used for storing the vertex numbers of vertices in the
//     SC components.
// For example if there are three SC components the vertices in each are stored
// as follows:
//   - SC0:  vertices[sets_s[0]] ... vertices[sets_f[0]].
//   - SC1:  vertices[sets_s[1]] ... vertices[sets_f[1]].
//   - SC2:  vertices[sets_s[2]] ... vertices[sets_f[2]].
// Note that array entries sets[3] onwards are set to UNDEFINED.
{
  int i, v, n;
  block_result_t *result;
  n = g->size;
  // accessed by block_depth_search()
  vertices = g->Number;
  // Allocate space for arrays to represent the search result.
  result = (block_result*)malloc(sizeof(block_result_t));
  block_vertices = result->vertices = (int*)malloc(n * sizeof(int));
  sets_s = result->sets_s = (int*)malloc(n *sizeof(int));
  sets_f = result->sets_f = (int*)malloc(n *sizeof(int));
  pos_sc = result->order = (int*)malloc(n * sizeof(int));
  result->ordered = (int*)malloc(n * sizeof(int));
  // Allocate space for arrays used while generating the result.
  block_stack = (int*)malloc(n * sizeof(int));
  sp = (int*)malloc(n * sizeof(int));
  visit_nos = (int*)malloc(n * sizeof(int));
  low_link_nos = (int*)malloc(n * sizeof(int));
  // Initialise necessary array entries to UNDEFINED.
  //  - sets_s[] and sets_f[] array entries are UNDEFINED, until data is
  //    written into them.
  //  - visit_nos[] array entries are UNDEFINED, until a vertex has been
  //    visited,
  //  - sp[v] is UNDEFINED unless v is in the stack.
  for(i = 0; i < n; i++)
    {
      sets_s[i] = sets_f[i] = visit_nos[i] = sp[i] = UNDEFINED;
      pos_sc[i] = i;
    }

  // Array sizes in the result structure.
  result->size = n;
  // Tarjan's algorithm proceeds as a recursive depth first search.  Note
  // that the block_depth_search() function accesses the current graph through the
  // global variable `vertices'.  If parts of the graph were not reached
  // block_depth_search() will be called again, until all vertices have been
  // reached.
  tos = 0;
  n_written = n_visited = 0;
  n_sets = 0;
  for(v = 0; v < n; v++)
    {
      n_sc_set=0;
      if(visit_nos[v] == UNDEFINED)
        block_depth_search(v);
    }
  result->n_sets = n_sets;
  for(i = 0; i < n_sets; i++)
    result->ordered[result->order[i]]=i;
  // free space taken up by arrays used while generating the result.
  free(block_stack);
  free(sp);
  free(visit_nos);
  free(low_link_nos);
  return result;
}



void
Blocks::block_result_free(block_result_t *r)
{
  free(r->vertices);
  free(r->sets_s);
  free(r->sets_f);
  free(r->order);
  free(r->ordered);
  free(r);
}


void
Blocks::block_result_print(block_result_t *r)
{
  int i, j, n_sets;

  n_sets = r->n_sets;

  cout << n_sets << " SC components:\n\n";
  for(i = 0; i < n_sets; i++)
    {
      cout << "SC" << r->order[i] << " = ";
      for(j = r->sets_s[i]; j <= r->sets_f[i]; j++)
        {
          cout << r->vertices[j] << " ";
        }
      cout << "\n";
    }
  for(i = 0; i < n_sets; i++)
    {
      cout << "SC" << i << " = ";
      for(j = r->sets_s[r->ordered[i]]; j <= r->sets_f[r->ordered[i]]; j++)
        {
          cout << r->vertices[j] << " ";
        }
      cout << "\n";
    }
}


void
Blocks::block_result_to_IM(block_result_t *r,bool* IM,int prologue, int n,simple* Index_Equ_IM,simple* Index_Var_IM)
{
  int i, j, k, l;
  bool* SIM=(bool*)malloc(n*n*sizeof(*SIM));
  simple* Index_Equ_IM_tmp=(simple*)malloc(n*sizeof(*Index_Equ_IM_tmp));
  simple* Index_Var_IM_tmp=(simple*)malloc(n*sizeof(*Index_Var_IM_tmp));
  for(i=0;i<n*n;i++)
    SIM[i]=IM[i];
  for(i=0;i<n;i++)
    {
      Index_Equ_IM_tmp[i].index=Index_Equ_IM[i].index;
      Index_Var_IM_tmp[i].index=Index_Var_IM[i].index;
    }
  l=prologue;
  for(i = 0; i < r->n_sets; i++)
    {
      for(j = r->sets_s[r->ordered[i]]; j <= r->sets_f[r->ordered[i]]; j++)
        {
          Index_Equ_IM[l].index=Index_Equ_IM_tmp[r->vertices[j]+prologue].index;
          for(k=0;k<n;k++)
            SIM[l*n+k]=IM[(r->vertices[j]+prologue)*n+k];
          l++;
        }
    }
  for(i=0;i<n*n;i++)
    IM[i]=SIM[i];
  l=prologue;
  for(i = 0; i < r->n_sets; i++)
    {
      for(j = r->sets_s[r->ordered[i]]; j <= r->sets_f[r->ordered[i]]; j++)
        {
          Index_Var_IM[l].index=Index_Var_IM_tmp[r->vertices[j]+prologue].index;
          for(k=0;k<n;k++)
            IM[k*n+l]=SIM[(k*n+r->vertices[j]+prologue)];
          l++;
        }
    }
  free(Index_Equ_IM_tmp);
  free(Index_Var_IM_tmp);
  free(SIM);
}


Equation_set*
Blocks::Equation_gr_IM( int n , bool* IM)
{
  Equation_set *g;
  Equation_vertex *vertices;
  Edge *edge_ptr;
  int i,j;
  g = (Equation_set*)malloc(sizeof(Equation_set));
  vertices = g->Number = (Equation_vertex*)malloc(n*sizeof(Equation_vertex));
  g->size = n;
  for(i = 0; i < n; i++)
    {
      vertices[i].First_Edge = NULL;
      for(j=0; j<n;j++)
        {
          if (IM[j*n+i])
            {
              if (vertices[i].First_Edge==NULL)
                {
                  vertices[i].First_Edge=(Edge*)malloc(sizeof(Edge));
                  edge_ptr=vertices[i].First_Edge;
                  edge_ptr->Vertex_Index=j;
                  edge_ptr->next= NULL;
                }
              else
                {
                  edge_ptr=(Edge*)malloc(sizeof(Edge));
                  edge_ptr->Vertex_Index=j;
                  edge_ptr->next= NULL;
                }

            }
        }
    }
  return g;
}

void
Blocks::Print_Equation_gr(Equation_set* Equation)
{
  int i;
  Edge *e1, *e2;
  cout << "The oriented graph of the model (earth blocks only) \n";
  cout << "equation | links\n";
  for(i=0;i<Equation->size;i++)
    {
      cout << "  " << i << "        ";
      e1=Equation->Number[i].First_Edge;
      while(e1!=NULL)
        {
          e2=e1->next;
          cout << e1->Vertex_Index << " ";
          e1=e2;
        }
      cout << "\n";
    }
}
