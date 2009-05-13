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


#include <iostream>
#include <sstream>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <functional>
#include "MinimumFeedbackSet.hh"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/max_cardinality_matching.hpp>
#include <boost/graph/strong_components.hpp>
//------------------------------------------------------------------------------
#include "BlockTriangular.hh"
//------------------------------------------------------------------------------
using namespace std;
using namespace boost;
using namespace MFS;




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
BlockTriangular::Prologue_Epilogue(bool* IM, int &prologue, int &epilogue, int n, vector<int> &Index_Var_IM, vector<int> &Index_Equ_IM, bool* IM0)
{
  bool modifie = 1;
  int i, j, k, l = 0;
  /*Looking for a prologue */
  prologue = 0;
  while (modifie)
    {
      modifie = 0;
      for (i = prologue;i < n;i++)
        {
          k = 0;
          for (j = prologue;j < n;j++)
            {
              if (IM[i*n + j])
                {
                  k++;
                  l = j;
                }
            }
          if ((k == 1) && IM0[Index_Equ_IM[i]*n + Index_Var_IM[l]])
            {
              modifie = 1;
              incidencematrix.swap_IM_c(IM, prologue, i, l, Index_Var_IM, Index_Equ_IM, n);
              prologue++;
            }
        }
    }
  epilogue = 0;
  modifie = 1;
  while (modifie)
    {
      modifie = 0;
      for (i = prologue;i < n - epilogue;i++)
        {
          k = 0;
          for (j = prologue;j < n - epilogue;j++)
            {
              if (IM[j*n + i])
                {
                  k++;
                  l = j;
                }
            }
          if ((k == 1) && IM0[Index_Equ_IM[l]*n + Index_Var_IM[i]])
            {
              modifie = 1;
              incidencematrix.swap_IM_c(IM, n - (1 + epilogue), l, i, Index_Var_IM, Index_Equ_IM, n);
              epilogue++;
            }
        }
    }
}


//------------------------------------------------------------------------------
// Find a matching between equations and endogenous variables
bool
BlockTriangular::Compute_Normalization(bool *IM, int equation_number, int prologue, int epilogue, bool verbose, bool *IM0, vector<int> &Index_Equ_IM) const
  {
    int n = equation_number - prologue - epilogue;


    typedef adjacency_list<vecS, vecS, undirectedS> BipartiteGraph;

    /*
      Vertices 0 to n-1 are for endogenous (using type specific ID)
      Vertices n to 2*n-1 are for equations (using equation no.)
    */
    BipartiteGraph g(2 * n);

    // Fill in the graph
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        if (IM0[(i+prologue) * equation_number+j+prologue])
          add_edge(i + n, j, g);


    // Compute maximum cardinality matching
    typedef vector<graph_traits<BipartiteGraph>::vertex_descriptor> mate_map_t;
    mate_map_t mate_map(2*n);

    bool check = checked_edmonds_maximum_cardinality_matching(g, &mate_map[0]);
    //cout << "check = " << check << "\n";
    if (check)
      {
        // Check if all variables are normalized
        mate_map_t::const_iterator it = find(mate_map.begin(), mate_map.begin() + n, graph_traits<BipartiteGraph>::null_vertex());
        if (it != mate_map.begin() + n)
          {
            if (verbose)
              {
                cerr << "ERROR: Could not normalize static model. Variable "
                << symbol_table.getName(symbol_table.getID(eEndogenous, it - mate_map.begin()))
                << " is not in the maximum cardinality matching." << endl;
                exit(EXIT_FAILURE);
              }
            return false;
          }
        for (int i = 0; i < n; i++)
          {
            int j = Index_Equ_IM[mate_map[i]-n+prologue];
            Index_Equ_IM[mate_map[i]-n+prologue] = Index_Equ_IM[i+prologue];
            Index_Equ_IM[i+prologue] = j;
            for (int k=0; k<n; k++)
              {
                j = IM[(mate_map[i]-n+prologue)*equation_number+k +prologue];
                IM[(mate_map[i]-n+prologue)*equation_number+k +prologue] = IM[(i+prologue)*equation_number+k +prologue];
                IM[(i+prologue)*equation_number+k +prologue] = j;
              }
          }
      }
    return check;
  }

void
BlockTriangular::Compute_Block_Decomposition_and_Feedback_Variables_For_Each_Block(bool *IM, int nb_var, int prologue, int epilogue, vector<int> &Index_Equ_IM, vector<int> &Index_Var_IM, vector<pair<int, int> > &blocks, bool verbose_) const
  {
    int n = nb_var - prologue - epilogue;
    bool *AMp;
    AMp = (bool*) malloc(n*n*sizeof(bool));
    //transforms the incidence matrix of the complet model into an adjancent matrix of the non-recursive part of the model
    for (int i=prologue; i<nb_var - epilogue; i++)
      for (int j=prologue; j<nb_var - epilogue; j++)
        if (Index_Equ_IM[i]!=Index_Var_IM[j])
          AMp[(i-prologue)*n+j-prologue] = IM[i*nb_var + j];
        else
          AMp[(i-prologue)*n+j-prologue] = 0;
    //In a first step we compute the strong components of the graph representation of the static model.
    // This insures that block are dynamically recursives.
    GraphvizDigraph G2 = AM_2_GraphvizDigraph(AMp, n);
    vector<int> component(num_vertices(G2)), discover_time(num_vertices(G2));

    //write_graphviz(cout, G2);

    int num = strong_components(G2, &component[0]);

    blocks = vector<pair<int, int> >(num, make_pair(0,0));

    //This vector contains for each block:
    //   - first set = equations belonging to the block,
    //   - second set = the feeback variables,
    //   - third vector = the reordered non-feedback variables.
    vector<pair<set<int>, pair<set<int>, vector<int> > > > components_set(num);

    for (unsigned int i=0; i<component.size(); i++)
      {
        blocks[component[i]].first++;
        components_set[component[i]].first.insert(i);
      }

    //For each block, the minimum set of feedback variable is computed
    // and the non-feedback variables are reordered to get
    // a sub-recursive block without feedback variables
    vector<int> tmp_Index_Equ_IM(Index_Equ_IM), tmp_Index_Var_IM(Index_Var_IM);
    int order = prologue;
    bool *SIM;
    SIM = (bool*)malloc(nb_var*nb_var*sizeof(bool));
    memcpy(SIM, IM, nb_var*nb_var*sizeof(bool));
    for (int i=0; i<num; i++)
      {
        AdjacencyList_type G = GraphvizDigraph_2_AdjacencyList(G2, components_set[i].first);
        set<int> feed_back_vertices;
        AdjacencyList_type G1 = Minimal_set_of_feedback_vertex(feed_back_vertices, G);
        property_map<AdjacencyList_type, vertex_index_t>::type v_index = get(vertex_index, G);
        components_set[i].second.first = feed_back_vertices;
        blocks[i].second = feed_back_vertices.size();
        vector<int> Reordered_Vertice;
        Reordered_Vertice = Reorder_the_recursive_variables(G, feed_back_vertices);
        for (vector<int>::iterator its=Reordered_Vertice.begin();its != Reordered_Vertice.end(); its++)
          {
            Index_Equ_IM[order] = tmp_Index_Equ_IM[*its+prologue];
            Index_Var_IM[order] = tmp_Index_Var_IM[*its+prologue];
            order++;
          }
        components_set[i].second.second = Reordered_Vertice;

        for (set<int>::iterator its=feed_back_vertices.begin();its != feed_back_vertices.end(); its++)
          {
            Index_Equ_IM[order] = tmp_Index_Equ_IM[v_index[vertex(*its, G)]+prologue];
            Index_Var_IM[order] = tmp_Index_Var_IM[v_index[vertex(*its, G)]+prologue];
            order++;
          }
      }
  }

void
BlockTriangular::Allocate_Block(int size, int *count_Equ, int count_Block, BlockType type, BlockSimulationType SimType, Model_Block * ModelBlock)
{
  int i, j, k, l, ls, m, i_1, Lead, Lag, first_count_equ, i1, li;
  int *tmp_size, *tmp_size_other_endo, *tmp_size_exo, *tmp_var, *tmp_endo, *tmp_other_endo, *tmp_exo, tmp_nb_other_endo, tmp_nb_exo, nb_lead_lag_endo;
  bool *tmp_variable_evaluated;
  bool *Cur_IM;
  bool *IM, OK;
  int Lag_Endo, Lead_Endo, Lag_Exo, Lead_Exo, Lag_Other_Endo, Lead_Other_Endo;

  ModelBlock->Periods = periods;
  ModelBlock->Block_List[count_Block].is_linear=true;
  ModelBlock->Block_List[count_Block].Size = size;
  ModelBlock->Block_List[count_Block].Type = type;
  ModelBlock->Block_List[count_Block].Temporary_InUse=new temporary_terms_inuse_type ();
  ModelBlock->Block_List[count_Block].Temporary_InUse->clear();
  ModelBlock->Block_List[count_Block].Simulation_Type = SimType;
  ModelBlock->Block_List[count_Block].Equation = (int*)malloc(ModelBlock->Block_List[count_Block].Size * sizeof(int));
  ModelBlock->Block_List[count_Block].Variable = (int*)malloc(ModelBlock->Block_List[count_Block].Size * sizeof(int));
  ModelBlock->Block_List[count_Block].Temporary_Terms_in_Equation = (temporary_terms_type**)malloc(ModelBlock->Block_List[count_Block].Size * sizeof(temporary_terms_type));
  ModelBlock->Block_List[count_Block].Own_Derivative = (int*)malloc(ModelBlock->Block_List[count_Block].Size * sizeof(int));
  Lead = Lag = 0;
  first_count_equ = *count_Equ;
  tmp_var = (int*)malloc(size * sizeof(int));
  tmp_endo = (int*)malloc((incidencematrix.Model_Max_Lead + incidencematrix.Model_Max_Lag + 1) * sizeof(int));
  tmp_other_endo = (int*)malloc(symbol_table.endo_nbr() * sizeof(int));
  tmp_size = (int*)malloc((incidencematrix.Model_Max_Lead + incidencematrix.Model_Max_Lag + 1) * sizeof(int));
  //cout << "tmp_size = (int*)malloc((incidencematrix.Model_Max_Lead + incidencematrix.Model_Max_Lag + 1= " << incidencematrix.Model_Max_Lead + incidencematrix.Model_Max_Lag + 1 << ") * sizeof(int))\n";
  tmp_size_other_endo = (int*)malloc((incidencematrix.Model_Max_Lead + incidencematrix.Model_Max_Lag + 1) * sizeof(int));
  tmp_size_exo = (int*)malloc((incidencematrix.Model_Max_Lead + incidencematrix.Model_Max_Lag + 1) * sizeof(int));
  memset(tmp_size_exo, 0, (incidencematrix.Model_Max_Lead + incidencematrix.Model_Max_Lag + 1)*sizeof(int));
  memset(tmp_size_other_endo, 0, (incidencematrix.Model_Max_Lead + incidencematrix.Model_Max_Lag + 1)*sizeof(int));
  memset(tmp_size, 0, (incidencematrix.Model_Max_Lead + incidencematrix.Model_Max_Lag + 1)*sizeof(int));
  memset(tmp_endo, 0, (incidencematrix.Model_Max_Lead + incidencematrix.Model_Max_Lag + 1)*sizeof(int));
  memset(tmp_other_endo, 0, symbol_table.endo_nbr()*sizeof(int));
  nb_lead_lag_endo = 0;
  Lag_Endo = Lead_Endo = Lag_Other_Endo = Lead_Other_Endo = Lag_Exo = Lead_Exo = 0;

  //Variable by variable looking for all leads and lags its occurence in each equation of the block
  tmp_variable_evaluated = (bool*)malloc(symbol_table.endo_nbr()*sizeof(bool));
  memset(tmp_variable_evaluated, 0, symbol_table.endo_nbr()*sizeof(bool));
  for (i = 0;i < size;i++)
    {
      ModelBlock->Block_List[count_Block].Temporary_Terms_in_Equation[i]=new temporary_terms_type ();
      ModelBlock->Block_List[count_Block].Temporary_Terms_in_Equation[i]->clear();
      ModelBlock->Block_List[count_Block].Equation[i] = Index_Equ_IM[*count_Equ];
      ModelBlock->Block_List[count_Block].Variable[i] = Index_Var_IM[*count_Equ];
      i_1 = Index_Var_IM[*count_Equ];
      for (k = -incidencematrix.Model_Max_Lag_Endo; k<=incidencematrix.Model_Max_Lead_Endo; k++)
        {
          Cur_IM = incidencematrix.Get_IM(k, eEndogenous);
          if (Cur_IM)
            {
              OK = false;
              if (k >= 0)
                {
                  for (j = 0;j < size;j++)
                    {
                      if (Cur_IM[i_1 + Index_Equ_IM[first_count_equ + j]*symbol_table.endo_nbr()])
                        {
                          tmp_variable_evaluated[i_1] = true;
                          tmp_size[incidencematrix.Model_Max_Lag_Endo + k]++;
                          if (!OK)
                            {
                              tmp_endo[incidencematrix.Model_Max_Lag + k]++;
                              nb_lead_lag_endo++;
                              OK = true;
                            }
                          if (k > Lead)
                            Lead = k;
                        }
                    }
                }
              else
                {
                  for (j = 0;j < size;j++)
                    {
                      if (Cur_IM[i_1 + Index_Equ_IM[first_count_equ + j]*symbol_table.endo_nbr()])
                        {
                          tmp_variable_evaluated[i_1] = true;
                          tmp_size[incidencematrix.Model_Max_Lag_Endo + k]++;
                          if (!OK)
                            {
                              tmp_variable_evaluated[i_1] = true;
                              tmp_endo[incidencematrix.Model_Max_Lag + k]++;
                              nb_lead_lag_endo++;
                              OK = true;
                            }
                          if (-k > Lag)
                            Lag = -k;
                        }
                    }
                }
            }
        }
      (*count_Equ)++;
    }
  Lag_Endo = Lag;
  Lead_Endo = Lead;

  tmp_nb_other_endo = 0;
  for (i = 0;i < size;i++)
    {
      for (k = -incidencematrix.Model_Max_Lag_Endo; k<=incidencematrix.Model_Max_Lead_Endo; k++)
        {
          Cur_IM = incidencematrix.Get_IM(k, eEndogenous);
          if (Cur_IM)
            {
              i_1 = Index_Equ_IM[first_count_equ+i] * symbol_table.endo_nbr();
              for (j = 0;j < symbol_table.endo_nbr();j++)
                if (Cur_IM[i_1 + j])
                  {
                    if (!tmp_variable_evaluated[j])
                      {
                        tmp_other_endo[j] = 1;
                        tmp_nb_other_endo++;
                      }
                    if (k>0 && k>Lead_Other_Endo)
                      Lead_Other_Endo = k;
                    else if (k<0 && (-k)>Lag_Other_Endo)
                      Lag_Other_Endo = -k;
                    if (k>0 && k>Lead)
                      Lead = k;
                    else if (k<0 && (-k)>Lag)
                      Lag = -k;
                    tmp_size_other_endo[k+incidencematrix.Model_Max_Lag_Endo]++;
                  }
            }
        }
    }
  ModelBlock->Block_List[count_Block].nb_other_endo = tmp_nb_other_endo;
  ModelBlock->Block_List[count_Block].Other_Endogenous = (int*)malloc(tmp_nb_other_endo * sizeof(int));


  tmp_exo = (int*)malloc(symbol_table.exo_nbr() * sizeof(int));
  memset(tmp_exo, 0, symbol_table.exo_nbr() *     sizeof(int));
  tmp_nb_exo = 0;
  for (i = 0;i < size;i++)
    {
      for (k = -incidencematrix.Model_Max_Lag_Exo; k<=incidencematrix.Model_Max_Lead_Exo; k++)
        {
          Cur_IM = incidencematrix.Get_IM(k, eExogenous);
          if (Cur_IM)
            {
              i_1 = Index_Equ_IM[first_count_equ+i] * symbol_table.exo_nbr();
              for (j=0;j<symbol_table.exo_nbr();j++)
                if (Cur_IM[i_1 + j])
                  {
                    if (!tmp_exo[j])
                      {
                        tmp_exo[j] = 1;
                        tmp_nb_exo++;
                      }
                    if (k>0 && k>Lead_Exo)
                      Lead_Exo = k;
                    else if (k<0 && (-k)>Lag_Exo)
                      Lag_Exo = -k;
                    if (k>0 && k>Lead)
                      Lead = k;
                    else if (k<0 && (-k)>Lag)
                      Lag = -k;
                    tmp_size_exo[k+incidencematrix.Model_Max_Lag_Exo]++;
                  }
            }
        }
    }


  ModelBlock->Block_List[count_Block].nb_exo = tmp_nb_exo;
  ModelBlock->Block_List[count_Block].Exogenous = (int*)malloc(tmp_nb_exo * sizeof(int));
  k = 0;
  for (j=0;j<symbol_table.exo_nbr();j++)
    if (tmp_exo[j])
      {
        ModelBlock->Block_List[count_Block].Exogenous[k] = j;
        k++;
      }

  ModelBlock->Block_List[count_Block].nb_exo_det = 0;

  ModelBlock->Block_List[count_Block].Max_Lag = Lag;
  ModelBlock->Block_List[count_Block].Max_Lead = Lead;
  ModelBlock->Block_List[count_Block].Max_Lag_Endo = Lag_Endo;
  ModelBlock->Block_List[count_Block].Max_Lead_Endo = Lead_Endo;
  ModelBlock->Block_List[count_Block].Max_Lag_Other_Endo = Lag_Other_Endo;
  ModelBlock->Block_List[count_Block].Max_Lead_Other_Endo = Lead_Other_Endo;
  ModelBlock->Block_List[count_Block].Max_Lag_Exo = Lag_Exo;
  ModelBlock->Block_List[count_Block].Max_Lead_Exo = Lead_Exo;
  ModelBlock->Block_List[count_Block].IM_lead_lag = (IM_compact*)malloc((Lead + Lag + 1) * sizeof(IM_compact));
  ls = l = li = size;
  i1 = 0;
  ModelBlock->Block_List[count_Block].Nb_Lead_Lag_Endo = nb_lead_lag_endo;
  for (i = 0;i < Lead + Lag + 1;i++)
    {
      if (incidencematrix.Model_Max_Lag_Endo-Lag+i>=0)
        {
          ModelBlock->Block_List[count_Block].IM_lead_lag[i].size = tmp_size[incidencematrix.Model_Max_Lag_Endo - Lag + i];
          ModelBlock->Block_List[count_Block].IM_lead_lag[i].nb_endo = tmp_endo[incidencematrix.Model_Max_Lag_Endo - Lag + i];
          ModelBlock->Block_List[count_Block].IM_lead_lag[i].u = (int*)malloc(tmp_size[incidencematrix.Model_Max_Lag_Endo - Lag + i] * sizeof(int));
          ModelBlock->Block_List[count_Block].IM_lead_lag[i].us = (int*)malloc(tmp_size[incidencematrix.Model_Max_Lag_Endo - Lag + i] * sizeof(int));
          ModelBlock->Block_List[count_Block].IM_lead_lag[i].Var = (int*)malloc(tmp_size[incidencematrix.Model_Max_Lag_Endo - Lag + i] * sizeof(int));
          ModelBlock->Block_List[count_Block].IM_lead_lag[i].Equ = (int*)malloc(tmp_size[incidencematrix.Model_Max_Lag_Endo - Lag + i] * sizeof(int));
          ModelBlock->Block_List[count_Block].IM_lead_lag[i].Var_Index = (int*)malloc(tmp_size[incidencematrix.Model_Max_Lag_Endo - Lag + i] * sizeof(int));
          ModelBlock->Block_List[count_Block].IM_lead_lag[i].Equ_Index = (int*)malloc(tmp_size[incidencematrix.Model_Max_Lag_Endo - Lag + i] * sizeof(int));
          ModelBlock->Block_List[count_Block].IM_lead_lag[i].size_other_endo = tmp_size_other_endo[incidencematrix.Model_Max_Lag_Endo - Lag + i];
          ModelBlock->Block_List[count_Block].IM_lead_lag[i].nb_other_endo = tmp_other_endo[incidencematrix.Model_Max_Lag_Endo - Lag + i];
          ModelBlock->Block_List[count_Block].IM_lead_lag[i].u_other_endo = (int*)malloc(tmp_size_other_endo[incidencematrix.Model_Max_Lag_Endo - Lag + i] * sizeof(int));
          ModelBlock->Block_List[count_Block].IM_lead_lag[i].Var_other_endo = (int*)malloc(tmp_size_other_endo[incidencematrix.Model_Max_Lag_Endo - Lag + i] * sizeof(int));
          ModelBlock->Block_List[count_Block].IM_lead_lag[i].Equ_other_endo = (int*)malloc(tmp_size_other_endo[incidencematrix.Model_Max_Lag_Endo - Lag + i] * sizeof(int));
          ModelBlock->Block_List[count_Block].IM_lead_lag[i].Var_Index_other_endo = (int*)malloc(tmp_size_other_endo[incidencematrix.Model_Max_Lag_Endo - Lag + i] * sizeof(int));
          ModelBlock->Block_List[count_Block].IM_lead_lag[i].Equ_Index_other_endo = (int*)malloc(tmp_size_other_endo[incidencematrix.Model_Max_Lag_Endo - Lag + i] * sizeof(int));
        }
      else
        ModelBlock->Block_List[count_Block].IM_lead_lag[i].size = 0;
      /*if (incidencematrix.Model_Max_Lag_Exo-Lag+i>=0)
        {
          ModelBlock->Block_List[count_Block].IM_lead_lag[i].size_exo = tmp_size_exo[incidencematrix.Model_Max_Lag_Exo - Lag + i];
          ModelBlock->Block_List[count_Block].IM_lead_lag[i].Exogenous = (int*)malloc(tmp_size_exo[incidencematrix.Model_Max_Lag_Exo - Lag + i] * sizeof(int));
          ModelBlock->Block_List[count_Block].IM_lead_lag[i].Exogenous_Index = (int*)malloc(tmp_size_exo[incidencematrix.Model_Max_Lag_Exo - Lag + i] * sizeof(int));
          ModelBlock->Block_List[count_Block].IM_lead_lag[i].Equ_X = (int*)malloc(tmp_size_exo[incidencematrix.Model_Max_Lag_Exo - Lag + i] * sizeof(int));
          ModelBlock->Block_List[count_Block].IM_lead_lag[i].Equ_X_Index = (int*)malloc(tmp_size_exo[incidencematrix.Model_Max_Lag_Exo - Lag + i] * sizeof(int));
        }
      else
        ModelBlock->Block_List[count_Block].IM_lead_lag[i].size_exo = 0;*/
      ModelBlock->Block_List[count_Block].IM_lead_lag[i].u_init = l;
      memset(tmp_variable_evaluated, 0, symbol_table.endo_nbr()*sizeof(bool));
      IM = incidencematrix.Get_IM(i - Lag, eEndogenous);
      if (IM)
        {
          for (j = first_count_equ;j < size + first_count_equ;j++)
            {
              i_1 = Index_Var_IM[j];
              m = 0;
              for (k = first_count_equ;k < size + first_count_equ;k++)
                if (IM[i_1 + Index_Equ_IM[k] * symbol_table.endo_nbr()])
                  m++;
              if (m > 0)
                {
                  tmp_var[j - first_count_equ] = i1;
                  i1++;
                }
            }
          m = 0;
          for (j = first_count_equ;j < size + first_count_equ;j++)
            {
              i_1 = Index_Equ_IM[j] * symbol_table.endo_nbr();
              for (k = first_count_equ;k < size + first_count_equ;k++)
                if (IM[Index_Var_IM[k] + i_1])
                  {
                    if (i == Lag)
                      {
                        ModelBlock->Block_List[count_Block].IM_lead_lag[i].us[m] = ls;
                        ls++;
                      }
                    ModelBlock->Block_List[count_Block].IM_lead_lag[i].u[m] = li;
                    ModelBlock->Block_List[count_Block].IM_lead_lag[i].Equ[m] = j - first_count_equ;
                    ModelBlock->Block_List[count_Block].IM_lead_lag[i].Var[m] = k - first_count_equ;
                    ModelBlock->Block_List[count_Block].IM_lead_lag[i].Equ_Index[m] = Index_Equ_IM[j];
                    ModelBlock->Block_List[count_Block].IM_lead_lag[i].Var_Index[m] = Index_Var_IM[k];
                    tmp_variable_evaluated[Index_Var_IM[k]] = true;
                    l++;
                    m++;
                    li++;
                  }
            }
          ModelBlock->Block_List[count_Block].IM_lead_lag[i].u_finish = li - 1;
          m = 0;
          for (j = first_count_equ;j < size + first_count_equ;j++)
            {
              i_1 = Index_Equ_IM[j] * symbol_table.endo_nbr();
              for (k = 0;k < symbol_table.endo_nbr();k++)
                if ((!tmp_variable_evaluated[Index_Var_IM[k]]) && IM[Index_Var_IM[k] + i_1])
                  {
                    ModelBlock->Block_List[count_Block].IM_lead_lag[i].u_other_endo[m] = l;
                    ModelBlock->Block_List[count_Block].IM_lead_lag[i].Equ_other_endo[m] = j - first_count_equ;
                    ModelBlock->Block_List[count_Block].IM_lead_lag[i].Var_other_endo[m] = k;
                    ModelBlock->Block_List[count_Block].IM_lead_lag[i].Equ_Index_other_endo[m] = Index_Equ_IM[j];
                    ModelBlock->Block_List[count_Block].IM_lead_lag[i].Var_Index_other_endo[m] = Index_Var_IM[k];
                    l++;
                    m++;
                  }
            }
          ModelBlock->Block_List[count_Block].IM_lead_lag[i].size_other_endo = m;
        }
      /*IM = incidencematrix.Get_IM(i - Lag, eExogenous);
      if (IM)
        {
          m = 0;
          for (j = first_count_equ;j < size + first_count_equ;j++)
            {
              i_1 = Index_Equ_IM[j] * symbol_table.exo_nbr();
              for (k = 0; k<tmp_nb_exo; k++)
                {
                  if (IM[ModelBlock->Block_List[count_Block].Exogenous[k]+i_1])
                    {
                      ModelBlock->Block_List[count_Block].IM_lead_lag[i].Exogenous[m] = k;
                      ModelBlock->Block_List[count_Block].IM_lead_lag[i].Exogenous_Index[m] = ModelBlock->Block_List[count_Block].Exogenous[k];
                      ModelBlock->Block_List[count_Block].IM_lead_lag[i].Equ_X[m] = j - first_count_equ;
                      ModelBlock->Block_List[count_Block].IM_lead_lag[i].Equ_X_Index[m] = Index_Equ_IM[j];
                      m++;
                    }
                }
            }
        }*/
    }
  free(tmp_size);
  free(tmp_size_other_endo);
  free(tmp_size_exo);
  free(tmp_endo);
  free(tmp_other_endo);
  free(tmp_exo);
  free(tmp_var);
  free(tmp_variable_evaluated);
}


void
BlockTriangular::Free_Block(Model_Block* ModelBlock) const
  {
    int blk, i;
    for (blk = 0;blk < ModelBlock->Size;blk++)
      {


        free(ModelBlock->Block_List[blk].Equation);
        free(ModelBlock->Block_List[blk].Variable);
        free(ModelBlock->Block_List[blk].Exogenous);
        free(ModelBlock->Block_List[blk].Own_Derivative);
        free(ModelBlock->Block_List[blk].Other_Endogenous);
        for (i = 0;i < ModelBlock->Block_List[blk].Max_Lag + ModelBlock->Block_List[blk].Max_Lead + 1;i++)
          {
            if (incidencematrix.Model_Max_Lag_Endo-ModelBlock->Block_List[blk].Max_Lag+i>=0 /*&& ModelBlock->Block_List[blk].IM_lead_lag[i].size*/)
              {
                free(ModelBlock->Block_List[blk].IM_lead_lag[i].u);
                free(ModelBlock->Block_List[blk].IM_lead_lag[i].us);
                free(ModelBlock->Block_List[blk].IM_lead_lag[i].Equ);
                free(ModelBlock->Block_List[blk].IM_lead_lag[i].Var);
                free(ModelBlock->Block_List[blk].IM_lead_lag[i].Equ_Index);
                free(ModelBlock->Block_List[blk].IM_lead_lag[i].Var_Index);
                free(ModelBlock->Block_List[blk].IM_lead_lag[i].u_other_endo);
                free(ModelBlock->Block_List[blk].IM_lead_lag[i].Var_other_endo);
                free(ModelBlock->Block_List[blk].IM_lead_lag[i].Equ_other_endo);
                free(ModelBlock->Block_List[blk].IM_lead_lag[i].Var_Index_other_endo);
                free(ModelBlock->Block_List[blk].IM_lead_lag[i].Equ_Index_other_endo);
              }
            /*if (incidencematrix.Model_Max_Lag_Exo-ModelBlock->Block_List[blk].Max_Lag+i>=0 )
              {
                free(ModelBlock->Block_List[blk].IM_lead_lag[i].Exogenous);
                free(ModelBlock->Block_List[blk].IM_lead_lag[i].Exogenous_Index);
                free(ModelBlock->Block_List[blk].IM_lead_lag[i].Equ_X_Index);
                free(ModelBlock->Block_List[blk].IM_lead_lag[i].Equ_X);
              }*/
          }
        free(ModelBlock->Block_List[blk].IM_lead_lag);
        for (i=0; i<ModelBlock->Block_List[blk].Size; i++)
          delete ModelBlock->Block_List[blk].Temporary_Terms_in_Equation[i];
        free(ModelBlock->Block_List[blk].Temporary_Terms_in_Equation);
        delete(ModelBlock->Block_List[blk].Temporary_InUse);
      }
    free(ModelBlock->Block_List);
    free(ModelBlock);
  }


t_type
BlockTriangular::Reduce_Blocks_and_type_determination(int prologue, int epilogue, vector<pair<int, int> > &blocks, vector<BinaryOpNode *> equations )
{
  int i=0;
  NodeID lhs, rhs;
  ostringstream tmp_output;
  BinaryOpNode *eq_node;
  ostringstream tmp_s;
  temporary_terms_type temporary_terms;
  int count_equ = 0, blck_count_simult =0;
  int Blck_Size;
  int Lead, Lag;
  t_type Type;
  bool *Cur_IM;
  BlockSimulationType Simulation_Type , prev_Type=UNKNOWN;
  for ( i=0; i<prologue+(int)blocks.size()+epilogue; i++)
    {
      int first_count_equ = count_equ;
      if (i < prologue)
        Blck_Size = 1;
      else if (i < prologue+(int)blocks.size())
        {
          Blck_Size = blocks[blck_count_simult].first;
          blck_count_simult++;
        }
      else if (i < prologue+(int)blocks.size()+epilogue)
        Blck_Size = 1;

      Lag = Lead = 0;
      for (count_equ = first_count_equ; count_equ < Blck_Size+first_count_equ; count_equ++)
        {
          int i_1 = Index_Var_IM[count_equ];
          for (int k = -incidencematrix.Model_Max_Lag_Endo; k<=incidencematrix.Model_Max_Lead_Endo; k++)
            {
              Cur_IM = incidencematrix.Get_IM(k, eEndogenous);
              if (Cur_IM)
                {
                  for (int j = 0;j < Blck_Size;j++)
                    {
                      if (Cur_IM[i_1 + Index_Equ_IM[first_count_equ + j] * symbol_table.endo_nbr()])
                        {
                          if (k > Lead)
                            Lead = k;
                          else if (-k > Lag)
                            Lag = -k;
                        }
                    }
                }
            }
        }
      if ((Lag > 0) && (Lead > 0))
        {
          if (Blck_Size == 1)
            Simulation_Type = SOLVE_TWO_BOUNDARIES_SIMPLE;
          else
            Simulation_Type = SOLVE_TWO_BOUNDARIES_COMPLETE;
        }
      else if (Blck_Size > 1)
        {
          if (Lead > 0)
            Simulation_Type = SOLVE_BACKWARD_COMPLETE;
          else
            Simulation_Type = SOLVE_FORWARD_COMPLETE;
        }
      else
        {
          if (Lead > 0)
            Simulation_Type = SOLVE_BACKWARD_SIMPLE;
          else
            Simulation_Type = SOLVE_FORWARD_SIMPLE;
        }
      if (Blck_Size == 1)
        {
          temporary_terms.clear();
          eq_node = equations[Index_Equ_IM[first_count_equ]];
          lhs = eq_node->get_arg1();
          rhs = eq_node->get_arg2();
          tmp_s.str("");
          tmp_output.str("");
          lhs->writeOutput(tmp_output, oMatlabDynamicModelSparse, temporary_terms);
          tmp_s << "y(it_, " << Index_Var_IM[first_count_equ]+1 << ")";
          //Determine whether the equation could be evaluated rather than to be solved
          if (tmp_output.str()==tmp_s.str())
            {
              if (Simulation_Type==SOLVE_BACKWARD_SIMPLE)
                Simulation_Type=EVALUATE_BACKWARD;
              else if (Simulation_Type==SOLVE_FORWARD_SIMPLE)
                Simulation_Type=EVALUATE_FORWARD;
            }
          else
            {
              tmp_output.str("");
              rhs->writeOutput(tmp_output, oMatlabDynamicModelSparse, temporary_terms);
              if (tmp_output.str()==tmp_s.str())
                {
                  if (Simulation_Type==SOLVE_BACKWARD_SIMPLE)
                    Simulation_Type=EVALUATE_BACKWARD_R;
                  else if (Simulation_Type==SOLVE_FORWARD_SIMPLE)
                    Simulation_Type=EVALUATE_FORWARD_R;
                }
            }
          if (i > 0)
            {
              if ( ((prev_Type ==  EVALUATE_FORWARD_R || prev_Type ==  EVALUATE_FORWARD) && (Simulation_Type == EVALUATE_FORWARD_R || Simulation_Type == EVALUATE_FORWARD))
                   || ((prev_Type ==  EVALUATE_BACKWARD_R || prev_Type ==  EVALUATE_BACKWARD) && (Simulation_Type == EVALUATE_BACKWARD_R || Simulation_Type == EVALUATE_BACKWARD))
                 )
                {
                  BlockSimulationType c_Type = (Type[Type.size()-1]).first;
                  int c_Size = (Type[Type.size()-1]).second;
                  Type[Type.size()-1]=make_pair(c_Type, ++c_Size);
                }
              else
                Type.push_back(make_pair(Simulation_Type, Blck_Size));
            }
          else
            Type.push_back(make_pair(Simulation_Type, Blck_Size));
        }
      else
        {
          Type.push_back(make_pair(Simulation_Type, Blck_Size));
        }
      prev_Type = Simulation_Type;
    }
  return(Type);
}


void
BlockTriangular::Normalize_and_BlockDecompose(bool* IM, Model_Block* ModelBlock, int n, int &prologue, int &epilogue, vector<int> &Index_Var_IM, vector<int> &Index_Equ_IM, bool* IM_0, jacob_map j_m, vector<BinaryOpNode *> equations )
{
  int i, j, Nb_TotalBlocks, Nb_RecursBlocks, Nb_SimulBlocks;
  BlockType Btype;
  int count_Block, count_Equ;
  //block_result_t* res;
  //Equation_set* Equation_gr = (Equation_set*) malloc(sizeof(Equation_set));
  bool* SIM0, *SIM00;
  SIM0 = (bool*)malloc(n * n * sizeof(bool));
  memcpy(SIM0,IM_0,n*n*sizeof(bool));
  Prologue_Epilogue(IM, prologue, epilogue, n, Index_Var_IM, Index_Equ_IM, SIM0);
  free(SIM0);
  int counted=0;
  if (prologue+epilogue<n)
    {
      cout << "Normalizing the model ...\n";
      double* max_val=(double*)malloc(n*sizeof(double));
      memset(max_val,0,n*sizeof(double));
      for ( map< pair< int, int >, double >::iterator iter = j_m.begin(); iter != j_m.end(); iter++ )
        {
          if (fabs(iter->second)>max_val[iter->first.first])
            max_val[iter->first.first]=fabs(iter->second);
        }
      for ( map< pair< int, int >, double >::iterator iter = j_m.begin(); iter != j_m.end(); iter++ )
        iter->second/=max_val[iter->first.first];
      free(max_val);
      bool OK=false;
      double bi=0.99999999;
      int suppressed=0;
      while (!OK && bi>1e-14)
        {
          int suppress=0;
          SIM0 = (bool*)malloc(n * n * sizeof(bool));
          memset(SIM0,0,n*n*sizeof(bool));
          SIM00 = (bool*)malloc(n * n * sizeof(bool));
          memset(SIM00,0,n*n*sizeof(bool));
          for ( map< pair< int, int >, double >::iterator iter = j_m.begin(); iter != j_m.end(); iter++ )
            {
              if (fabs(iter->second)>bi)
                {
                  SIM0[iter->first.first*n+iter->first.second]=1;
                  if (!IM_0[iter->first.first*n+iter->first.second])
                    {
                      cout << "Error nothing at IM_0[" << iter->first.first << ", " << iter->first.second << "]=" << IM_0[iter->first.first*n+iter->first.second] << "  " << iter->second << "\n";
                    }
                }
              else
                suppress++;
            }
          for (i = 0;i < n;i++)
            for (j = 0;j < n;j++)
              {
                SIM00[i*n + j] = SIM0[Index_Equ_IM[i] * n + Index_Var_IM[j]];
              }
          free(SIM0);
          if (suppress!=suppressed)
            OK = Compute_Normalization(IM, n, prologue, epilogue, false, SIM00, Index_Equ_IM);
          suppressed=suppress;
          if (!OK)
            //bi/=1.07;
            bi/=3;
          counted++;
          if (bi>1e-14)
            free(SIM00);
        }
      if (!OK)
        Compute_Normalization(IM, n, prologue, epilogue, true, SIM00, Index_Equ_IM);

    }

  cout << "Finding the optimal block decomposition of the model ...\n";
  vector<pair<int, int> > blocks;
  if (prologue+epilogue<n)
    Compute_Block_Decomposition_and_Feedback_Variables_For_Each_Block(IM, n, prologue, epilogue, Index_Equ_IM, Index_Var_IM, blocks, false);

  t_type  Type = Reduce_Blocks_and_type_determination(prologue, epilogue, blocks, equations);

  i = 0;
  j = 0;
  Nb_SimulBlocks = 0;
  int Nb_fv = 0;
  for (t_type::const_iterator it = Type.begin(); it!=Type.end(); it++)
    {
      if (it->first==SOLVE_FORWARD_COMPLETE || it->first==SOLVE_BACKWARD_COMPLETE || it->first==SOLVE_TWO_BOUNDARIES_COMPLETE)
        {
          Nb_SimulBlocks++;
          if (it->second>j)
            {
              j=it->second;
              Nb_fv = blocks[Nb_SimulBlocks-1].second;
            }
        }
    }

  Nb_TotalBlocks = Type.size();
  Nb_RecursBlocks = Nb_TotalBlocks - Nb_SimulBlocks;
  cout << Nb_TotalBlocks << " block(s) found:\n";
  cout << "  " << Nb_RecursBlocks << " recursive block(s) and " << blocks.size() << " simultaneous block(s). \n";
  cout << "  the largest simultaneous block has " << j     << " equation(s)\n" <<
  "                                 and " << Nb_fv << " feedback variable(s).\n";

  ModelBlock->Size = Nb_TotalBlocks;
  ModelBlock->Periods = periods;
  ModelBlock->Block_List = (Block*)malloc(sizeof(ModelBlock->Block_List[0]) * Nb_TotalBlocks);

  count_Equ = count_Block = 0;
  for (t_type::const_iterator it = Type.begin(); it!=Type.end(); it++)
    {
      if (count_Equ<prologue)
        Btype = PROLOGUE;
      else if (count_Equ<n-epilogue)
        if (it->second==1)
          Btype = PROLOGUE;
        else
          Btype = SIMULTANS;
      else
        Btype = EPILOGUE;
      Allocate_Block(it->second, &count_Equ, count_Block++, Btype, it->first, ModelBlock);
    }
}


//------------------------------------------------------------------------------
// normalize each equation of the dynamic model
// and find the optimal block triangular decomposition of the static model
void
BlockTriangular::Normalize_and_BlockDecompose_Static_0_Model(const jacob_map &j_m, vector<BinaryOpNode *> equations)
{
  bool* SIM, *SIM_0;
  bool* Cur_IM;
  int i, k, size;
  //First create a static model incidence matrix
  size = symbol_table.endo_nbr() * symbol_table.endo_nbr() * sizeof(*SIM);
  SIM = (bool*)malloc(size);
  for (i = 0; i< symbol_table.endo_nbr() * symbol_table.endo_nbr(); i++) SIM[i] = 0;
  for (k = -incidencematrix.Model_Max_Lag_Endo; k<=incidencematrix.Model_Max_Lead_Endo; k++)
    {
      Cur_IM = incidencematrix.Get_IM(k, eEndogenous);
      if (Cur_IM)
        {
          for (i = 0;i < symbol_table.endo_nbr()*symbol_table.endo_nbr();i++)
            {
              SIM[i] = (SIM[i]) || (Cur_IM[i]);
            }
        }
    }
  if (bt_verbose)
    {
      cout << "incidence matrix for the static model (unsorted) \n";
      incidencematrix.Print_SIM(SIM, eEndogenous);
    }
  //Index_Equ_IM = (simple*)malloc(symbol_table.endo_nbr() * sizeof(*Index_Equ_IM));
  Index_Equ_IM = vector<int>(symbol_table.endo_nbr());
  for (i = 0;i < symbol_table.endo_nbr();i++)
    {
      Index_Equ_IM[i] = i;
    }
  Index_Var_IM = vector<int>(symbol_table.endo_nbr());
  for (i = 0;i < symbol_table.endo_nbr();i++)
    {
      Index_Var_IM[i] = i;
    }
  if (ModelBlock != NULL)
    Free_Block(ModelBlock);
  ModelBlock = (Model_Block*)malloc(sizeof(*ModelBlock));
  Cur_IM = incidencematrix.Get_IM(0, eEndogenous);
  SIM_0 = (bool*)malloc(symbol_table.endo_nbr() * symbol_table.endo_nbr() * sizeof(*SIM_0));
  for (i = 0;i < symbol_table.endo_nbr()*symbol_table.endo_nbr();i++)
    SIM_0[i] = Cur_IM[i];
  Normalize_and_BlockDecompose(SIM, ModelBlock, symbol_table.endo_nbr(), prologue, epilogue, Index_Var_IM, Index_Equ_IM, SIM_0, j_m, equations);
  free(SIM_0);
  free(SIM);
}
