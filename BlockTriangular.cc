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
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <functional>
#include "MinimumFeedbackSet.hh"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/max_cardinality_matching.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/topological_sort.hpp>

//------------------------------------------------------------------------------
#include "BlockTriangular.hh"
//------------------------------------------------------------------------------
using namespace std;
using namespace boost;
using namespace MFS;

BlockTriangular::BlockTriangular(SymbolTable &symbol_table_arg, NumericalConstants &num_const_arg)
  : symbol_table(symbol_table_arg),
  //normalization(symbol_table_arg),
  incidencematrix(symbol_table_arg),
  num_const(num_const_arg)
{
  bt_verbose = 0;
  ModelBlock = NULL;
  periods = 0;
  prologue = epilogue = 0;
  Normalized_Equation = new DataTree(symbol_table, num_const);
}

BlockTriangular::~BlockTriangular()
{
  delete Normalized_Equation;
}

//------------------------------------------------------------------------------
// Find the prologue and the epilogue of the model
void
BlockTriangular::Prologue_Epilogue(bool *IM, int &prologue, int &epilogue, int n, vector<int> &Index_Var_IM, vector<int> &Index_Equ_IM, bool *IM0)
{
  bool modifie = 1;
  int i, j, k, l = 0;
  /*Looking for a prologue */
  prologue = 0;
  while (modifie)
    {
      modifie = 0;
      for (i = prologue; i < n; i++)
        {
          k = 0;
          for (j = prologue; j < n; j++)
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
      for (i = prologue; i < n - epilogue; i++)
        {
          k = 0;
          for (j = prologue; j < n - epilogue; j++)
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
BlockTriangular::Compute_Normalization(bool *IM, int equation_number, int prologue, int epilogue, int verbose, bool *IM0, vector<int> &Index_Equ_IM) const
{
  int n = equation_number - prologue - epilogue;

  typedef adjacency_list<vecS, vecS, undirectedS> BipartiteGraph;


  if(verbose == 2)
    cout << "trying to normlized even in singular case\n";
  /*
     Vertices 0 to n-1 are for endogenous (using type specific ID)
     Vertices n to 2*n-1 are for equations (using equation no.)
   */
  BipartiteGraph g(2 * n);
  // Fill in the graph
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      if (IM0[(i+prologue) * equation_number+j+prologue])
        {
        	//printf("equation=%3d variable=%3d\n",i,j);
          add_edge(i + n, j, g);
        }

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
          if (verbose == 1)
            {
              cout << "WARNING: Could not normalize dynamic model. Variable "
                   << symbol_table.getName(symbol_table.getID(eEndogenous, it - mate_map.begin()))
                   << " is not in the maximum cardinality matching. Trying to find a singular normalization." << endl;
              //exit(EXIT_FAILURE);
              return false;
            }
          else if (verbose == 2)
            {
              cerr << "ERROR: Could not normalize dynamic model (even with a singularity). Variable "
                   << symbol_table.getName(symbol_table.getID(eEndogenous, it - mate_map.begin()))
                   << " is not in the maximum cardinality matching." << endl;
              exit(EXIT_FAILURE);
            }
          return false;
        }
      vector<int> Index_Equ_IM_tmp(Index_Equ_IM);
      bool *SIM;
      SIM = (bool *) malloc(equation_number*equation_number*sizeof(bool));
      memcpy(SIM, IM, equation_number*equation_number*sizeof(bool));
      for (int i = 0; i < n; i++)
        {
          //printf("match equation %4d with variable %4d \n", mate_map[i] - n, i);
          Index_Equ_IM[i + prologue] = Index_Equ_IM_tmp[mate_map[i] - n + prologue];
          for (int k = 0; k < n; k++)
            IM[(i+prologue)*equation_number+k +prologue] = SIM[(mate_map[i]-n+prologue)*equation_number+k +prologue];
        }
      free(SIM);
    }
  return check;
}




t_vtype
BlockTriangular::Get_Variable_LeadLag_By_Block(vector<int > &components_set, int nb_blck_sim, int prologue, int epilogue, t_vtype &equation_lead_lag) const
{
  int nb_endo = symbol_table.endo_nbr();
  vector<int> variable_blck(nb_endo), equation_blck(nb_endo);
  t_vtype Variable_Type(nb_endo);
  for (int i = 0; i < nb_endo; i++)
    {
      if (i < prologue)
        {
          variable_blck[Index_Var_IM[i]] = i;
          equation_blck[Index_Equ_IM[i]] = i;
        }
      else if (i < (int)components_set.size() + prologue)
        {
          variable_blck[Index_Var_IM[i]] = components_set[i-prologue] + prologue;
          equation_blck[Index_Equ_IM[i]] = components_set[i-prologue] + prologue;
        }
      else
        {
          variable_blck[Index_Var_IM[i]] = i- (nb_endo - nb_blck_sim - prologue - epilogue);
          equation_blck[Index_Equ_IM[i]] = i- (nb_endo - nb_blck_sim - prologue - epilogue);
        }
      Variable_Type[i] = make_pair(0, 0);
    }
  equation_lead_lag = Variable_Type;
  for (int k = -incidencematrix.Model_Max_Lag_Endo; k <= incidencematrix.Model_Max_Lead_Endo; k++)
    {
      bool *Cur_IM = incidencematrix.Get_IM(k, eEndogenous);
      if (Cur_IM)
        {
          for (int i = 0; i < nb_endo; i++)
            {
              int i_1 = Index_Var_IM[i];
              for (int j = 0; j < nb_endo; j++)
                {
                  int j_l = Index_Equ_IM[ j];
                  if (Cur_IM[i_1 + Index_Equ_IM[ j] * nb_endo] and variable_blck[i_1] == equation_blck[j_l])
                    {
                      if (k > Variable_Type[i_1].second)
                        Variable_Type[i_1] = make_pair(Variable_Type[i_1].first, k);
                      if (k < -Variable_Type[i_1].first)
                        Variable_Type[i_1] = make_pair(-k, Variable_Type[i_1].second);
                      if (k > equation_lead_lag[j_l].second)
                        equation_lead_lag[j_l] = make_pair(equation_lead_lag[j_l].first, k);
                      if (k < -equation_lead_lag[j_l].first)
                        equation_lead_lag[j_l] = make_pair(-k, equation_lead_lag[j_l].second);
                    }
                }
            }
        }
    }
  return (Variable_Type);
}

void
BlockTriangular::Compute_Block_Decomposition_and_Feedback_Variables_For_Each_Block(bool *IM, int nb_var, int prologue, int epilogue, vector<int> &Index_Equ_IM, vector<int> &Index_Var_IM, vector<pair<int, int> > &blocks, t_etype &Equation_Type, bool verbose_, bool select_feedback_variable, int mfs) const
{
  t_vtype V_Variable_Type;
  int n = nb_var - prologue - epilogue;
  bool *AMp;
  AMp = (bool *) malloc(n*n*sizeof(bool));
  //transforms the incidence matrix of the complet model into an adjancent matrix of the non-recursive part of the model
  for (int i = prologue; i < nb_var - epilogue; i++)
    for (int j = prologue; j < nb_var - epilogue; j++)
      if (j != i)
        AMp[(i-prologue)*n+j-prologue] = IM[i*nb_var + j];
      else
        AMp[(i-prologue)*n+j-prologue] = 0;

  //In a first step we compute the strong components of the graph representation of the static model.
  // This insures that block are dynamically recursives.
  GraphvizDigraph G2 = AM_2_GraphvizDigraph(AMp, n);
  vector<int> endo2block(num_vertices(G2)), discover_time(num_vertices(G2));

  int num = strong_components(G2, &endo2block[0]);

  blocks = vector<pair<int, int> >(num, make_pair(0, 0));



/*New*/
 // Compute strongly connected components
  // Create directed acyclic graph associated to the strongly connected components
  typedef adjacency_list<vecS, vecS, directedS> DirectedGraph;
  DirectedGraph dag(num);
  /*graph_traits<DirectedGraph>::edge_iterator ei, ei_end;
  for(tie(ei, ei_end) = edges(G2); ei != ei_end; ++ei)
    {
      int s = endo2block[source(*ei, G2)];
      int t = endo2block[target(*ei, G2)];
      if (s != t)
        add_edge(s, t, dag);
    }*/
  for (int i = 0;i < num_vertices(G2);i++)
    {
      GraphvizDigraph::out_edge_iterator it_out, out_end;
      GraphvizDigraph::vertex_descriptor vi = vertex(i, G2);
      for (tie(it_out, out_end) = out_edges(vi, G2); it_out != out_end; ++it_out)
        {
          int t_b = endo2block[target(*it_out, G2)];
          int s_b = endo2block[source(*it_out, G2)];
          if (s_b != t_b)
            add_edge(s_b, t_b, dag);
        }
    }
  // Compute topological sort of DAG (ordered list of unordered SCC)
  deque<int> ordered2unordered;
  topological_sort(dag, front_inserter(ordered2unordered)); // We use a front inserter because topological_sort returns the inverse order
  // Construct mapping from unordered SCC to ordered SCC
  vector<int> unordered2ordered(num);
  for(int i = 0; i < num; i++)
    unordered2ordered[ordered2unordered[i]] = i;
/*EndNew*/


  //This vector contains for each block:
  //   - first set = equations belonging to the block,
  //   - second set = the feeback variables,
  //   - third vector = the reordered non-feedback variables.
  vector<pair<set<int>, pair<set<int>, vector<int> > > > components_set(num);

  for (unsigned int i = 0; i < endo2block.size(); i++)
    {
      endo2block[i] = unordered2ordered[endo2block[i]];
      blocks[endo2block[i]].first++;
      components_set[endo2block[i]].first.insert(i);
    }


  t_vtype equation_lead_lag;
  V_Variable_Type = Get_Variable_LeadLag_By_Block(endo2block, num, prologue, epilogue, equation_lead_lag);

  vector<int> tmp_Index_Equ_IM(Index_Equ_IM), tmp_Index_Var_IM(Index_Var_IM);
  int order = prologue;
  bool *SIM;
  SIM = (bool *) malloc(nb_var*nb_var*sizeof(bool));
  memcpy(SIM, IM, nb_var*nb_var*sizeof(bool));

  //Add a loop on vertices which could not be normalized or vertices related to lead variables => force those vertices to belong to the feedback set
  if(select_feedback_variable)
    for (int i = 0; i < n; i++)
      if (Equation_Type[Index_Equ_IM[i+prologue]].first == E_SOLVE or V_Variable_Type[Index_Var_IM[i+prologue]].second > 0 or V_Variable_Type[Index_Var_IM[i+prologue]].first > 0
                                                               or equation_lead_lag[Index_Equ_IM[i+prologue]].second > 0 or equation_lead_lag[Index_Equ_IM[i+prologue]].first > 0
                                                               or mfs == 0)
        add_edge(i, i, G2);
  else
    for (int i = 0; i < n; i++)
      if (Equation_Type[Index_Equ_IM[i+prologue]].first == E_SOLVE or mfs == 0)
        add_edge(i, i, G2);
  //For each block, the minimum set of feedback variable is computed
  // and the non-feedback variables are reordered to get
  // a sub-recursive block without feedback variables
  for (int i = 0; i < num; i++)
    {
      AdjacencyList_type G = GraphvizDigraph_2_AdjacencyList(G2, components_set[i].first);
      set<int> feed_back_vertices;
      //Print(G);
      AdjacencyList_type G1 = Minimal_set_of_feedback_vertex(feed_back_vertices, G);
      property_map<AdjacencyList_type, vertex_index_t>::type v_index = get(vertex_index, G);
      components_set[i].second.first = feed_back_vertices;
      blocks[i].second = feed_back_vertices.size();
      vector<int> Reordered_Vertice;
      Reorder_the_recursive_variables(G, feed_back_vertices, Reordered_Vertice);
      //First we have the recursive equations conditional on feedback variables
      for (vector<int>::iterator its = Reordered_Vertice.begin(); its != Reordered_Vertice.end(); its++)
        {
          Index_Equ_IM[order] = tmp_Index_Equ_IM[*its+prologue];
          Index_Var_IM[order] = tmp_Index_Var_IM[*its+prologue];
          order++;
        }
      components_set[i].second.second = Reordered_Vertice;
      //Second we have the equations related to the feedback variables
      for (set<int>::iterator its = feed_back_vertices.begin(); its != feed_back_vertices.end(); its++)
        {
          Index_Equ_IM[order] = tmp_Index_Equ_IM[v_index[vertex(*its, G)]+prologue];
          Index_Var_IM[order] = tmp_Index_Var_IM[v_index[vertex(*its, G)]+prologue];
          order++;
        }

    }
  free(AMp);
  free(SIM);
}

void
BlockTriangular::Allocate_Block(int size, int *count_Equ, int count_Block, BlockType type, BlockSimulationType SimType, Model_Block *ModelBlock, t_etype &Equation_Type, int recurs_Size, vector<int> &Index_Var_IM, vector<int> &Index_Equ_IM)
{
  int i, j, k, l, ls, m, i_1, Lead, Lag, first_count_equ, i1, li;
  int *tmp_size, *tmp_size_other_endo, *tmp_size_exo, *tmp_var, *tmp_endo, *tmp_other_endo, *tmp_exo, tmp_nb_other_endo, tmp_nb_exo, nb_lead_lag_endo;
  bool *tmp_variable_evaluated;
  bool *Cur_IM;
  bool *IM, OK;
  int Lag_Endo, Lead_Endo, Lag_Exo, Lead_Exo, Lag_Other_Endo, Lead_Other_Endo;
  //cout << "block " << count_Block << " size " << size << " SimType=" << BlockSim(SimType) << "\n";

  ModelBlock->Periods = periods;
  ModelBlock->Block_List[count_Block].is_linear = true;
  ModelBlock->Block_List[count_Block].Size = size;
  ModelBlock->Block_List[count_Block].Type = type;
  ModelBlock->Block_List[count_Block].Nb_Recursives = recurs_Size;
  ModelBlock->Block_List[count_Block].Temporary_InUse = new temporary_terms_inuse_type();
  ModelBlock->Block_List[count_Block].Chain_Rule_Derivatives = new chain_rule_derivatives_type();
  ModelBlock->Block_List[count_Block].Temporary_InUse->clear();
  ModelBlock->Block_List[count_Block].Simulation_Type = SimType;
  ModelBlock->Block_List[count_Block].Equation = (int *) malloc(ModelBlock->Block_List[count_Block].Size * sizeof(int));
  ModelBlock->Block_List[count_Block].Equation_Type = (EquationType *) malloc(ModelBlock->Block_List[count_Block].Size * sizeof(EquationType));
  ModelBlock->Block_List[count_Block].Equation_Normalized = (NodeID*)malloc(ModelBlock->Block_List[count_Block].Size * sizeof(NodeID));
  ModelBlock->Block_List[count_Block].Variable = (int *) malloc(ModelBlock->Block_List[count_Block].Size * sizeof(int));
  ModelBlock->Block_List[count_Block].Temporary_Terms_in_Equation = (temporary_terms_type **) malloc(ModelBlock->Block_List[count_Block].Size * sizeof(temporary_terms_type));
  ModelBlock->Block_List[count_Block].Own_Derivative = (int *) malloc(ModelBlock->Block_List[count_Block].Size * sizeof(int));
  Lead = Lag = 0;
  first_count_equ = *count_Equ;
  tmp_var = (int *) malloc(size * sizeof(int));
  tmp_endo = (int *) malloc((incidencematrix.Model_Max_Lead + incidencematrix.Model_Max_Lag + 1) * sizeof(int));
  tmp_other_endo = (int *) malloc((incidencematrix.Model_Max_Lead + incidencematrix.Model_Max_Lag + 1) * sizeof(int));
  tmp_size = (int *) malloc((incidencematrix.Model_Max_Lead + incidencematrix.Model_Max_Lag + 1) * sizeof(int));
  tmp_size_other_endo = (int *) malloc((incidencematrix.Model_Max_Lead + incidencematrix.Model_Max_Lag + 1) * sizeof(int));
  tmp_size_exo = (int *) malloc((incidencematrix.Model_Max_Lead + incidencematrix.Model_Max_Lag + 1) * sizeof(int));
  memset(tmp_size_exo, 0, (incidencematrix.Model_Max_Lead + incidencematrix.Model_Max_Lag + 1)*sizeof(int));
  memset(tmp_size_other_endo, 0, (incidencematrix.Model_Max_Lead + incidencematrix.Model_Max_Lag + 1)*sizeof(int));
  memset(tmp_size, 0, (incidencematrix.Model_Max_Lead + incidencematrix.Model_Max_Lag + 1)*sizeof(int));
  memset(tmp_endo, 0, (incidencematrix.Model_Max_Lead + incidencematrix.Model_Max_Lag + 1)*sizeof(int));
  memset(tmp_other_endo, 0, (incidencematrix.Model_Max_Lead + incidencematrix.Model_Max_Lag + 1)*sizeof(int));
  nb_lead_lag_endo = 0;
  Lag_Endo = Lead_Endo = Lag_Other_Endo = Lead_Other_Endo = Lag_Exo = Lead_Exo = 0;

  //Variable by variable looking for all leads and lags its occurence in each equation of the block
  tmp_variable_evaluated = (bool *) malloc(symbol_table.endo_nbr()*sizeof(bool));
  memset(tmp_variable_evaluated, 0, symbol_table.endo_nbr()*sizeof(bool));
  for (i = 0; i < size; i++)
    {
      ModelBlock->Block_List[count_Block].Temporary_Terms_in_Equation[i] = new temporary_terms_type();
      ModelBlock->Block_List[count_Block].Temporary_Terms_in_Equation[i]->clear();
      ModelBlock->Block_List[count_Block].Equation[i] = Index_Equ_IM[*count_Equ];
      ModelBlock->Block_List[count_Block].Variable[i] = Index_Var_IM[*count_Equ];
      ModelBlock->Block_List[count_Block].Equation_Type[i] = Equation_Type[Index_Equ_IM[*count_Equ]].first;
      ModelBlock->Block_List[count_Block].Equation_Normalized[i] = Equation_Type[Index_Equ_IM[*count_Equ]].second;
      /*if(Equation_Type[Index_Equ_IM[*count_Equ]].second)
         {
          temporary_terms_type temporary_terms;
          cout << "Equation_Type[Index_Equ_IM[*count_Equ]].second->get_op_code()=" << Equation_Type[Index_Equ_IM[*count_Equ]].second->get_op_code() << "\n";
          cout << "ModelBlock->Block_List[" << count_Block << "].Equation_Normalized[" << i << "]->get_op_code()=" << ModelBlock->Block_List[count_Block].Equation_Normalized[i]->get_op_code() << "\n";
          ModelBlock->Block_List[count_Block].Equation_Normalized[i]->writeOutput(cout, oMatlabDynamicModelSparse, temporary_terms);
         }*/
      i_1 = Index_Var_IM[*count_Equ];
      for (k = -incidencematrix.Model_Max_Lag_Endo; k <= incidencematrix.Model_Max_Lead_Endo; k++)
        {
          Cur_IM = incidencematrix.Get_IM(k, eEndogenous);
          if (Cur_IM)
            {
              OK = false;
              if (k >= 0)
                {
                  for (j = 0; j < size; j++)
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
                  for (j = 0; j < size; j++)
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
  for (i = 0; i < size; i++)
    {
      for (k = -incidencematrix.Model_Max_Lag_Endo; k <= incidencematrix.Model_Max_Lead_Endo; k++)
        {
          Cur_IM = incidencematrix.Get_IM(k, eEndogenous);
          if (Cur_IM)
            {
              i_1 = Index_Equ_IM[first_count_equ+i] * symbol_table.endo_nbr();
              for (j = 0; j < symbol_table.endo_nbr(); j++)
                if (Cur_IM[i_1 + j])
                  {
                    if (!tmp_variable_evaluated[j])
                      {
                      	tmp_other_endo[incidencematrix.Model_Max_Lag + k]++;
                        tmp_nb_other_endo++;
                      }
                    if (k > 0 && k > Lead_Other_Endo)
                      Lead_Other_Endo = k;
                    else if (k < 0 && (-k) > Lag_Other_Endo)
                      Lag_Other_Endo = -k;
                    if (k > 0 && k > Lead)
                      Lead = k;
                    else if (k < 0 && (-k) > Lag)
                      Lag = -k;
                    tmp_size_other_endo[k+incidencematrix.Model_Max_Lag_Endo]++;
                  }
            }
        }
    }
  ModelBlock->Block_List[count_Block].nb_other_endo = tmp_nb_other_endo;
  ModelBlock->Block_List[count_Block].Other_Endogenous = (int *) malloc(tmp_nb_other_endo * sizeof(int));

  tmp_exo = (int *) malloc(symbol_table.exo_nbr() * sizeof(int));
  memset(tmp_exo, 0, symbol_table.exo_nbr() *     sizeof(int));
  tmp_nb_exo = 0;
  for (i = 0; i < size; i++)
    {
      for (k = -incidencematrix.Model_Max_Lag_Exo; k <= incidencematrix.Model_Max_Lead_Exo; k++)
        {
          Cur_IM = incidencematrix.Get_IM(k, eExogenous);
          if (Cur_IM)
            {
              i_1 = Index_Equ_IM[first_count_equ+i] * symbol_table.exo_nbr();
              for (j = 0; j < symbol_table.exo_nbr(); j++)
                if (Cur_IM[i_1 + j])
                  {
                    if (!tmp_exo[j])
                      {
                        tmp_exo[j] = 1;
                        tmp_nb_exo++;
                      }
                    if (k > 0 && k > Lead_Exo)
                      Lead_Exo = k;
                    else if (k < 0 && (-k) > Lag_Exo)
                      Lag_Exo = -k;
                    if (k > 0 && k > Lead)
                      Lead = k;
                    else if (k < 0 && (-k) > Lag)
                      Lag = -k;
                    tmp_size_exo[k+incidencematrix.Model_Max_Lag_Exo]++;
                  }
            }
        }
    }

  ModelBlock->Block_List[count_Block].nb_exo = tmp_nb_exo;
  ModelBlock->Block_List[count_Block].Exogenous = (int *) malloc(tmp_nb_exo * sizeof(int));
  k = 0;
  for (j = 0; j < symbol_table.exo_nbr(); j++)
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
  ModelBlock->Block_List[count_Block].IM_lead_lag = (IM_compact *) malloc((Lead + Lag + 1) * sizeof(IM_compact));
  ls = l = li = size;
  i1 = 0;
  ModelBlock->Block_List[count_Block].Nb_Lead_Lag_Endo = nb_lead_lag_endo;
  for (i = 0; i < Lead + Lag + 1; i++)
    {
      if (incidencematrix.Model_Max_Lag_Endo-Lag+i >= 0)
        {
          ModelBlock->Block_List[count_Block].IM_lead_lag[i].size = tmp_size[incidencematrix.Model_Max_Lag_Endo - Lag + i];
          ModelBlock->Block_List[count_Block].IM_lead_lag[i].nb_endo = tmp_endo[incidencematrix.Model_Max_Lag_Endo - Lag + i];
          ModelBlock->Block_List[count_Block].IM_lead_lag[i].u = (int *) malloc(tmp_size[incidencematrix.Model_Max_Lag_Endo - Lag + i] * sizeof(int));
          ModelBlock->Block_List[count_Block].IM_lead_lag[i].us = (int *) malloc(tmp_size[incidencematrix.Model_Max_Lag_Endo - Lag + i] * sizeof(int));
          ModelBlock->Block_List[count_Block].IM_lead_lag[i].Var = (int *) malloc(tmp_size[incidencematrix.Model_Max_Lag_Endo - Lag + i] * sizeof(int));
          ModelBlock->Block_List[count_Block].IM_lead_lag[i].Equ = (int *) malloc(tmp_size[incidencematrix.Model_Max_Lag_Endo - Lag + i] * sizeof(int));
          ModelBlock->Block_List[count_Block].IM_lead_lag[i].Var_Index = (int *) malloc(tmp_size[incidencematrix.Model_Max_Lag_Endo - Lag + i] * sizeof(int));
          ModelBlock->Block_List[count_Block].IM_lead_lag[i].Equ_Index = (int *) malloc(tmp_size[incidencematrix.Model_Max_Lag_Endo - Lag + i] * sizeof(int));
          ModelBlock->Block_List[count_Block].IM_lead_lag[i].size_other_endo = tmp_size_other_endo[incidencematrix.Model_Max_Lag_Endo - Lag + i];
          ModelBlock->Block_List[count_Block].IM_lead_lag[i].nb_other_endo = tmp_other_endo[incidencematrix.Model_Max_Lag_Endo - Lag + i];
          ModelBlock->Block_List[count_Block].IM_lead_lag[i].u_other_endo = (int *) malloc(tmp_size_other_endo[incidencematrix.Model_Max_Lag_Endo - Lag + i] * sizeof(int));
          ModelBlock->Block_List[count_Block].IM_lead_lag[i].Var_other_endo = (int *) malloc(tmp_size_other_endo[incidencematrix.Model_Max_Lag_Endo - Lag + i] * sizeof(int));
          ModelBlock->Block_List[count_Block].IM_lead_lag[i].Equ_other_endo = (int *) malloc(tmp_size_other_endo[incidencematrix.Model_Max_Lag_Endo - Lag + i] * sizeof(int));
          ModelBlock->Block_List[count_Block].IM_lead_lag[i].Var_Index_other_endo = (int *) malloc(tmp_size_other_endo[incidencematrix.Model_Max_Lag_Endo - Lag + i] * sizeof(int));
          ModelBlock->Block_List[count_Block].IM_lead_lag[i].Equ_Index_other_endo = (int *) malloc(tmp_size_other_endo[incidencematrix.Model_Max_Lag_Endo - Lag + i] * sizeof(int));
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
          for (j = first_count_equ; j < size + first_count_equ; j++)
            {
              i_1 = Index_Var_IM[j];
              m = 0;
              for (k = first_count_equ; k < size + first_count_equ; k++)
                if (IM[i_1 + Index_Equ_IM[k] * symbol_table.endo_nbr()])
                  m++;
              if (m > 0)
                {
                  tmp_var[j - first_count_equ] = i1;
                  i1++;
                }
            }
          m = 0;
          for (j = first_count_equ; j < size + first_count_equ; j++)
            {
              i_1 = Index_Equ_IM[j] * symbol_table.endo_nbr();
              for (k = first_count_equ; k < size + first_count_equ; k++)
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
          for (j = first_count_equ; j < size + first_count_equ; j++)
            {
              i_1 = Index_Equ_IM[j] * symbol_table.endo_nbr();
              for (k = 0; k < symbol_table.endo_nbr(); k++)
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
BlockTriangular::Free_Block(Model_Block *ModelBlock) const
{
  int blk, i;
  for (blk = 0; blk < ModelBlock->Size; blk++)
    {
      free(ModelBlock->Block_List[blk].Equation);
      free(ModelBlock->Block_List[blk].Variable);
      free(ModelBlock->Block_List[blk].Exogenous);
      free(ModelBlock->Block_List[blk].Own_Derivative);
      free(ModelBlock->Block_List[blk].Other_Endogenous);
      free(ModelBlock->Block_List[blk].Equation_Type);
      free(ModelBlock->Block_List[blk].Equation_Normalized);
      for (i = 0; i < ModelBlock->Block_List[blk].Max_Lag + ModelBlock->Block_List[blk].Max_Lead + 1; i++)
        {
          if (incidencematrix.Model_Max_Lag_Endo-ModelBlock->Block_List[blk].Max_Lag+i >= 0 /*&& ModelBlock->Block_List[blk].IM_lead_lag[i].size*/)
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
      for (i = 0; i < ModelBlock->Block_List[blk].Size; i++)
        delete ModelBlock->Block_List[blk].Temporary_Terms_in_Equation[i];
      free(ModelBlock->Block_List[blk].Temporary_Terms_in_Equation);
      delete (ModelBlock->Block_List[blk].Temporary_InUse);
      delete ModelBlock->Block_List[blk].Chain_Rule_Derivatives;
    }
  free(ModelBlock->Block_List);
  free(ModelBlock);
}

t_etype
BlockTriangular::Equation_Type_determination(vector<BinaryOpNode *> &equations, map<pair<int, pair<int, int> >, NodeID> &first_order_endo_derivatives, vector<int> &Index_Var_IM, vector<int> &Index_Equ_IM, int mfs)
{
  NodeID lhs, rhs;
  ostringstream tmp_output;
  BinaryOpNode *eq_node;
  ostringstream tmp_s;
  temporary_terms_type temporary_terms;
  EquationType Equation_Simulation_Type;
  t_etype V_Equation_Simulation_Type(equations.size());
  for (unsigned int i = 0; i < equations.size(); i++)
    {
      temporary_terms.clear();
      int eq = Index_Equ_IM[i];
      int var = Index_Var_IM[i];
      eq_node = equations[eq];
      lhs = eq_node->get_arg1();
      rhs = eq_node->get_arg2();
      Equation_Simulation_Type = E_SOLVE;
      tmp_s.str("");
      tmp_output.str("");
      lhs->writeOutput(tmp_output, oMatlabDynamicModelSparse, temporary_terms);
      tmp_s << "y(it_, " << Index_Var_IM[i]+1 << ")";
      map<pair<int, pair<int, int> >, NodeID>::iterator derivative = first_order_endo_derivatives.find(make_pair(eq, make_pair(var, 0)));
      pair<bool, NodeID> res;
      if(derivative != first_order_endo_derivatives.end())
        {
          set<pair<int, int> > result;
          derivative->second->collectEndogenous(result);
          set<pair<int, int> >::const_iterator d_endo_variable = result.find(make_pair(var, 0));
          //Determine whether the equation could be evaluated rather than to be solved
          ostringstream tt("");
          derivative->second->writeOutput(tt, oMatlabDynamicModelSparse, temporary_terms);
          if (tmp_output.str() == tmp_s.str() and tt.str()=="1")
            {
              Equation_Simulation_Type = E_EVALUATE;
            }
          else
            {
        	    vector<pair<int, pair<NodeID, NodeID> > > List_of_Op_RHS;
              res =  equations[eq]->normalizeEquation(var, List_of_Op_RHS);
              if(mfs==2)
                {
                  if(d_endo_variable == result.end() && res.second)
                    Equation_Simulation_Type = E_EVALUATE_S;
                }
              else if(mfs==3)
                {
                  if(res.second) // The equation could be solved analytically
                    Equation_Simulation_Type = E_EVALUATE_S;
                }
            }
        }
      V_Equation_Simulation_Type[eq] = make_pair(Equation_Simulation_Type, dynamic_cast<BinaryOpNode *>(res.second));
    }
  return (V_Equation_Simulation_Type);
}

t_type
BlockTriangular::Reduce_Blocks_and_type_determination(int prologue, int epilogue, vector<pair<int, int> > &blocks, vector<BinaryOpNode *> &equations, t_etype &Equation_Type, vector<int> &Index_Var_IM, vector<int> &Index_Equ_IM)
{
  int i = 0;
  int count_equ = 0, blck_count_simult = 0;
  int Blck_Size, Recurs_Size;
  int Lead, Lag;
  t_type Type;
  bool *Cur_IM;
  BlockSimulationType Simulation_Type, prev_Type = UNKNOWN;
  int eq = 0;
  for (i = 0; i < prologue+(int) blocks.size()+epilogue; i++)
    {
      int first_count_equ = count_equ;
      if (i < prologue)
        {
          Blck_Size = 1;
          Recurs_Size = 0;
        }
      else if (i < prologue+(int) blocks.size())
        {
          Blck_Size = blocks[blck_count_simult].first;
          Recurs_Size = Blck_Size - blocks[blck_count_simult].second;
          blck_count_simult++;
        }
      else if (i < prologue+(int) blocks.size()+epilogue)
        {
          Blck_Size = 1;
          Recurs_Size = 0;
        }

      Lag = Lead = 0;
      for (count_equ = first_count_equ; count_equ < Blck_Size+first_count_equ; count_equ++)
        {
          int i_1 = Index_Var_IM[count_equ];
          for (int k = -incidencematrix.Model_Max_Lag_Endo; k <= incidencematrix.Model_Max_Lead_Endo; k++)
            {
              Cur_IM = incidencematrix.Get_IM(k, eEndogenous);
              if (Cur_IM)
                {
                  for (int j = 0; j < Blck_Size; j++)
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
          if (Equation_Type[Index_Equ_IM[eq]].first == E_EVALUATE /*or Equation_Type[Index_Equ_IM[eq]].first==E_EVALUATE_R*/ or Equation_Type[Index_Equ_IM[eq]].first == E_EVALUATE_S)
            {
              if (Simulation_Type == SOLVE_BACKWARD_SIMPLE)
                Simulation_Type = EVALUATE_BACKWARD;
              else if (Simulation_Type == SOLVE_FORWARD_SIMPLE)
                Simulation_Type = EVALUATE_FORWARD;
            }
          if (i > 0)
            {
              if ((prev_Type ==  EVALUATE_FORWARD and Simulation_Type == EVALUATE_FORWARD)
                  or (prev_Type ==  EVALUATE_BACKWARD and Simulation_Type == EVALUATE_BACKWARD))
                {
                  BlockSimulationType c_Type = (Type[Type.size()-1]).first;
                  int c_Size = (Type[Type.size()-1]).second.first;
                  Type[Type.size()-1] = make_pair(c_Type, make_pair(++c_Size, Type[Type.size()-1].second.second));
                }
              else
                Type.push_back(make_pair(Simulation_Type, make_pair(Blck_Size, Recurs_Size)));
            }
          else
            Type.push_back(make_pair(Simulation_Type, make_pair(Blck_Size, Recurs_Size)));
        }
      else
        {
          Type.push_back(make_pair(Simulation_Type, make_pair(Blck_Size, Recurs_Size)));
        }
      prev_Type = Simulation_Type;
      eq += Blck_Size;
    }
  return (Type);
}


map<pair<pair<int, pair<int, int> >, pair<int, int> >, int>
BlockTriangular::get_Derivatives(Model_Block *ModelBlock, int blck)
{
  map<pair<pair<int, pair<int, int> >, pair<int, int> >, int> Derivatives;
  Derivatives.clear();
  int nb_endo = symbol_table.endo_nbr();
  /*ModelBlock.Block_List[Blck].first_order_determinstic_simulation_derivatives = new*/
  for(int lag = -ModelBlock->Block_List[blck].Max_Lag; lag <= ModelBlock->Block_List[blck].Max_Lead; lag++)
    {
      bool *IM=incidencematrix.Get_IM(lag, eEndogenous);
      if(IM)
        {
          for(int eq = 0; eq < ModelBlock->Block_List[blck].Size; eq++)
            {
              int eqr = ModelBlock->Block_List[blck].Equation[eq];
              for(int var = 0; var < ModelBlock->Block_List[blck].Size; var++)
                {
                  int varr = ModelBlock->Block_List[blck].Variable[var];
                  /*cout << "IM=" << IM << "\n";
                  cout << "varr=" << varr << " eqr=" << eqr << " lag=" << lag << "\n";*/
                  if(IM[varr+eqr*nb_endo])
                    {
                      bool OK = true;
                      map<pair<pair<int, pair<int, int> >, pair<int, int> >, int>::const_iterator its = Derivatives.find(make_pair(make_pair(lag, make_pair(eq, var)), make_pair(eqr, varr)));
                      if(its!=Derivatives.end())
                        {
                        	if(its->second == 2)
                        	  OK=false;
                        }

                      if(OK)
                        {
                          if (ModelBlock->Block_List[blck].Equation_Type[eq] == E_EVALUATE_S and eq<ModelBlock->Block_List[blck].Nb_Recursives)
                            //It's a normalized equation, we have to recompute the derivative using chain rule derivative function*/
                            Derivatives[make_pair(make_pair(lag, make_pair(eq, var)), make_pair(eqr, varr))] = 1;
                          else
                            //It's a feedback equation we can use the derivatives
                            Derivatives[make_pair(make_pair(lag, make_pair(eq, var)), make_pair(eqr, varr))] = 0;
                        }
                      if(var<ModelBlock->Block_List[blck].Nb_Recursives)
                        {
                          int eqs = ModelBlock->Block_List[blck].Equation[var];
                          for(int vars=ModelBlock->Block_List[blck].Nb_Recursives; vars<ModelBlock->Block_List[blck].Size; vars++)
                            {
                              int varrs = ModelBlock->Block_List[blck].Variable[vars];
                              //A new derivative need to be computed using the chain rule derivative function (a feedback variable appear in a recursive equation)
                              if(Derivatives.find(make_pair(make_pair(lag, make_pair(var, vars)), make_pair(eqs, varrs)))!=Derivatives.end())
                                Derivatives[make_pair(make_pair(lag, make_pair(eq, vars)), make_pair(eqr, varrs))] = 2;
                            }
                        }
                    }
                }
            }
        }
    }
	return(Derivatives);
}

void
BlockTriangular::Normalize_and_BlockDecompose(bool *IM, Model_Block *ModelBlock, int n, int &prologue, int &epilogue, vector<int> &Index_Var_IM, vector<int> &Index_Equ_IM, bool *IM_0, jacob_map &j_m, vector<BinaryOpNode *> &equations, t_etype &V_Equation_Type, map<pair<int, pair<int, int> >, NodeID> &first_order_endo_derivatives, bool dynamic, int mfs, double cutoff)
{
  int i, j, Nb_TotalBlocks, Nb_RecursBlocks, Nb_SimulBlocks;
  BlockType Btype;
  int count_Block, count_Equ;
  bool *SIM0, *SIM00;



  int counted = 0;
  if (prologue+epilogue < n)
    {
      cout << "Normalizing the model ...\n";
      double *max_val = (double *) malloc(n*sizeof(double));
      memset(max_val, 0, n*sizeof(double));
      for (map< pair< int, int >, double >::iterator iter = j_m.begin(); iter != j_m.end(); iter++)
        {
          if (fabs(iter->second) > max_val[iter->first.first])
            max_val[iter->first.first] = fabs(iter->second);
        }
      for (map< pair< int, int >, double >::iterator iter = j_m.begin(); iter != j_m.end(); iter++)
        iter->second /= max_val[iter->first.first];
      free(max_val);
      bool OK = false;
      double bi = 0.99999999;
      //double bi=1e-13;
      int suppressed = 0;
      vector<int> Index_Equ_IM_save(Index_Equ_IM);
      while (!OK && bi > 1e-19)
        {
          int suppress = 0;
          Index_Equ_IM = Index_Equ_IM_save;
          SIM0 = (bool *) malloc(n * n * sizeof(bool));
          memset(SIM0, 0, n*n*sizeof(bool));
          SIM00 = (bool *) malloc(n * n * sizeof(bool));
          memset(SIM00, 0, n*n*sizeof(bool));
          //cout << "---------------------------------\n";
          for (map< pair< int, int >, double >::iterator iter = j_m.begin(); iter != j_m.end(); iter++)
            {
            	//printf("iter->second=% 1.10f iter->first.first=%3d iter->first.second=%3d  bi=%f\n", iter->second, iter->first.first, iter->first.second, bi);
              if (fabs(iter->second) > max(bi, cutoff))
                {
                  SIM0[iter->first.first*n+iter->first.second] = 1;
                  if (!IM_0[iter->first.first*n+iter->first.second])
                    {
                      cout << "Error nothing at IM_0[" << iter->first.first << ", " << iter->first.second << "]=" << IM_0[iter->first.first*n+iter->first.second] << "  " << iter->second << "\n";
                    }
                }
              else
                suppress++;
            }
          for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
              {
                SIM00[i*n + j] = SIM0[Index_Equ_IM[i] * n + Index_Var_IM[j]];
              }
          free(SIM0);
          if (suppress != suppressed)
            OK = Compute_Normalization(IM, n, prologue, epilogue, 0, SIM00, Index_Equ_IM);
          suppressed = suppress;
          if (!OK)
            //bi/=1.07;
            bi /= 2;
          counted++;
          if (bi > 1e-19)
            free(SIM00);
        }
      if (!OK)
        {
          Compute_Normalization(IM, n, prologue, epilogue, 1, SIM00, Index_Equ_IM);
          Compute_Normalization(IM, n, prologue, epilogue, 2, IM_0, Index_Equ_IM);
        }
    }
  SIM0 = (bool *) malloc(n * n * sizeof(bool));
  memcpy(SIM0, IM_0, n*n*sizeof(bool));
  Prologue_Epilogue(IM, prologue, epilogue, n, Index_Var_IM, Index_Equ_IM, SIM0);

  free(SIM0);

  V_Equation_Type = Equation_Type_determination(equations, first_order_endo_derivatives, Index_Var_IM, Index_Equ_IM, mfs);

  cout << "Finding the optimal block decomposition of the model ...\n";
  vector<pair<int, int> > blocks;
  if (prologue+epilogue < n)
    {
      if(dynamic)
        Compute_Block_Decomposition_and_Feedback_Variables_For_Each_Block(IM, n, prologue, epilogue, Index_Equ_IM, Index_Var_IM, blocks, V_Equation_Type, false, true, mfs);
      else
        Compute_Block_Decomposition_and_Feedback_Variables_For_Each_Block(IM, n, prologue, epilogue, Index_Equ_IM, Index_Var_IM, blocks, V_Equation_Type, false, false, mfs);
    }

  t_type  Type = Reduce_Blocks_and_type_determination(prologue, epilogue, blocks, equations, V_Equation_Type, Index_Var_IM, Index_Equ_IM);

  i = 0;
  j = 0;
  Nb_SimulBlocks = 0;
  int Nb_feedback_variable = 0;
  for (t_type::const_iterator it = Type.begin(); it != Type.end(); it++)
    {
      if (it->first == SOLVE_FORWARD_COMPLETE || it->first == SOLVE_BACKWARD_COMPLETE || it->first == SOLVE_TWO_BOUNDARIES_COMPLETE)
        {
          Nb_SimulBlocks++;
          if (it->second.first > j)
            {
              j = it->second.first;
              Nb_feedback_variable = blocks[Nb_SimulBlocks-1].second;
            }
        }
    }

  Nb_TotalBlocks = Type.size();
  Nb_RecursBlocks = Nb_TotalBlocks - Nb_SimulBlocks;
  cout << Nb_TotalBlocks << " block(s) found:\n";
  cout << "  " << Nb_RecursBlocks << " recursive block(s) and " << blocks.size() << " simultaneous block(s). \n";
  cout << "  the largest simultaneous block has " << j     << " equation(s)\n"
       <<"                                 and " << Nb_feedback_variable << " feedback variable(s).\n";

  ModelBlock->Size = Nb_TotalBlocks;
  ModelBlock->Periods = periods;
  ModelBlock->Block_List = (Block *) malloc(sizeof(ModelBlock->Block_List[0]) * Nb_TotalBlocks);

  count_Equ = count_Block = 0;



  for (t_type::const_iterator it = Type.begin(); it != Type.end(); it++)
    {
      if (count_Equ < prologue)
        Btype = PROLOGUE;
      else if (count_Equ < n-epilogue)
        if (it->second.first == 1)
          Btype = PROLOGUE;
        else
          {
            Btype = SIMULTANS;
          }
      else
        Btype = EPILOGUE;
      Allocate_Block(it->second.first, &count_Equ, count_Block++, Btype, it->first, ModelBlock, V_Equation_Type, it->second.second, Index_Var_IM, Index_Equ_IM);
    }
}

//------------------------------------------------------------------------------
// normalize each equation of the dynamic model
// and find the optimal block triangular decomposition of the static model
void
BlockTriangular::Normalize_and_BlockDecompose_Static_0_Model(jacob_map &j_m, vector<BinaryOpNode *> &equations, t_etype &equation_simulation_type, map<pair<int, pair<int, int> >, NodeID> &first_order_endo_derivatives, int mfs, double cutoff)
{
  bool *SIM, *SIM_0;
  bool *Cur_IM;
  int i, k, size;
  //First create a static model incidence matrix
  size = symbol_table.endo_nbr() * symbol_table.endo_nbr() * sizeof(*SIM);
  SIM = (bool *) malloc(size);
  for (i = 0; i < symbol_table.endo_nbr() * symbol_table.endo_nbr(); i++) SIM[i] = 0;
  for (k = -incidencematrix.Model_Max_Lag_Endo; k <= incidencematrix.Model_Max_Lead_Endo; k++)
    {
      Cur_IM = incidencematrix.Get_IM(k, eEndogenous);
      if (Cur_IM)
        {
          for (i = 0; i < symbol_table.endo_nbr()*symbol_table.endo_nbr(); i++)
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
  Index_Equ_IM = vector<int>(symbol_table.endo_nbr());
  for (i = 0; i < symbol_table.endo_nbr(); i++)
    {
      Index_Equ_IM[i] = i;
    }
  Index_Var_IM = vector<int>(symbol_table.endo_nbr());
  for (i = 0; i < symbol_table.endo_nbr(); i++)
    {
      Index_Var_IM[i] = i;
    }
  if (ModelBlock != NULL)
    Free_Block(ModelBlock);
  ModelBlock = (Model_Block *) malloc(sizeof(*ModelBlock));
  Cur_IM = incidencematrix.Get_IM(0, eEndogenous);
  SIM_0 = (bool *) malloc(symbol_table.endo_nbr() * symbol_table.endo_nbr() * sizeof(*SIM_0));
  for (i = 0; i < symbol_table.endo_nbr()*symbol_table.endo_nbr(); i++)
    SIM_0[i] = Cur_IM[i];
  Normalize_and_BlockDecompose(SIM, ModelBlock, symbol_table.endo_nbr(), prologue, epilogue, Index_Var_IM, Index_Equ_IM, SIM_0, j_m, equations, equation_simulation_type, first_order_endo_derivatives, true, mfs, cutoff);
  free(SIM_0);
  free(SIM);
}
