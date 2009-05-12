#ifndef MFE_BOOST
#define MFE_BOOST

#include <iostream>
#include <map>
#include <vector>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <stdio.h>
#include <stdlib.h>

using namespace std;
using namespace boost;

//#define verbose

  typedef property<vertex_index_t, int,
          property<vertex_index1_t, int,
          property<vertex_degree_t, int,
          property<vertex_in_degree_t, int,
          property<vertex_out_degree_t, int > > > > > VertexProperty;
  typedef adjacency_list<listS, listS, bidirectionalS, VertexProperty> AdjacencyList_type;
  typedef map<graph_traits<AdjacencyList_type>::vertex_descriptor,default_color_type> color_type;
  typedef vector<AdjacencyList_type::vertex_descriptor> vector_vertex_descriptor;

namespace MFS
{
  void Eliminate(AdjacencyList_type::vertex_descriptor vertex_to_eliminate, AdjacencyList_type& G);
  vector_vertex_descriptor Collect_Doublet(AdjacencyList_type::vertex_descriptor vertex, AdjacencyList_type& G);
  bool Vertex_Belong_to_a_Clique(AdjacencyList_type::vertex_descriptor vertex, AdjacencyList_type& G);
  bool Elimination_of_Vertex_With_One_or_Less_Indegree_or_Outdegree_Step(AdjacencyList_type& G);  //Graph reduction: eliminating purely intermediate variables or variables outside of any circuit
  bool Elimination_of_Vertex_belonging_to_a_clique_Step(AdjacencyList_type& G);                   //Graphe reduction: eliminaion of Cliques
  bool Suppression_of_Edge_i_j_if_not_a_loop_and_if_for_all_i_k_edge_we_have_a_k_j_edge_Step(AdjacencyList_type& G);  //Suppression
  bool Suppression_of_all_in_Edge_in_i_if_not_a_loop_and_if_all_doublet_i_eq_Min_inDegree_outDegree_Step(AdjacencyList_type& G);
  bool Suppression_of_Vertex_X_if_it_loops_store_in_set_of_feedback_vertex_Step(vector<pair<int, AdjacencyList_type::vertex_descriptor> > &looping_variable, AdjacencyList_type& G);
  void Print(AdjacencyList_type& G);
  AdjacencyList_type AM_2_AdjacencyList(bool* AMp,unsigned int n);
  void Print(GraphvizDigraph& G);
  GraphvizDigraph AM_2_GraphvizDigraph(bool* AM, unsigned int n);
  AdjacencyList_type GraphvizDigraph_2_AdjacencyList(GraphvizDigraph& G1, set<int> select_index);
  bool has_cycle_dfs(AdjacencyList_type& g, AdjacencyList_type::vertex_descriptor u, color_type& color, vector<int> &circuit_stack);
  bool has_cylce(AdjacencyList_type& g, vector<int> &circuit_stack, int size);
  bool has_cycle(vector<int> &circuit_stack, AdjacencyList_type& G);
  AdjacencyList_type Minimal_set_of_feedback_vertex(set<int> &feed_back_vertices, const AdjacencyList_type& G);
  void Suppress(AdjacencyList_type::vertex_descriptor vertex_to_eliminate, AdjacencyList_type& G);
  void Suppress(int vertex_num, AdjacencyList_type& G);
  vector<int> Reorder_the_recursive_variables(const AdjacencyList_type& G1, set<int> &feed_back_vertices);
};

#endif
