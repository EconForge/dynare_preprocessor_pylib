/*
 * Copyright (C) 2009 Dynare Team
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

#ifndef _MINIMUMFEEDBACKSET_HH
#define _MINIMUMFEEDBACKSET_HH

#include <map>
#include <vector>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/adjacency_list.hpp>

using namespace std;
using namespace boost;

namespace MFS
{
  typedef property<vertex_index_t, int,
          property<vertex_index1_t, int,
          property<vertex_degree_t, int,
          property<vertex_in_degree_t, int,
          property<vertex_out_degree_t, int > > > > > VertexProperty;
  typedef adjacency_list<listS, listS, bidirectionalS, VertexProperty> AdjacencyList_type;
  typedef map<graph_traits<AdjacencyList_type>::vertex_descriptor,default_color_type> color_type;
  typedef vector<AdjacencyList_type::vertex_descriptor> vector_vertex_descriptor;

  //! Eliminate a vertex i
  /*! For a vertex i replace all edges e_k_i and e_i_j by a shorcut e_k_j and then Suppress the vertex i*/
  void Eliminate(AdjacencyList_type::vertex_descriptor vertex_to_eliminate, AdjacencyList_type& G);
  //! Collect all doublets (edges e_i_k such that there is an edge e_k_i with k!=i in the graph)
  /*! Returns the vector of doublets */
  vector_vertex_descriptor Collect_Doublet(AdjacencyList_type::vertex_descriptor vertex, AdjacencyList_type& G);
  //! Detect all the clique (all vertex in a clique are related to each other) in the graph
  bool Vertex_Belong_to_a_Clique(AdjacencyList_type::vertex_descriptor vertex, AdjacencyList_type& G);
  //! Graph reduction: eliminating purely intermediate variables or variables outside of any circuit
  bool Elimination_of_Vertex_With_One_or_Less_Indegree_or_Outdegree_Step(AdjacencyList_type& G);
  //! Graph reduction: elimination of a vertex inside a clique
  bool Elimination_of_Vertex_belonging_to_a_clique_Step(AdjacencyList_type& G);
  //! A vertex belong to the feedback vertex set if the vertex loops on itself.
  /*! We have to suppress this vertex and store it into the feedback set.*/
  bool Suppression_of_Vertex_X_if_it_loops_store_in_set_of_feedback_vertex_Step(set<int> &feed_back_vertices, AdjacencyList_type& G1);
  //! Print the Graph
  void Print(GraphvizDigraph& G);
  void Print(AdjacencyList_type& G);
  //! Create a GraphvizDigraph from a Adjacency Matrix (an incidence Matrix without the diagonal terms)
  GraphvizDigraph AM_2_GraphvizDigraph(bool* AM, unsigned int n);
  //! Create an adjacency graph from a Adjacency Matrix (an incidence Matrix without the diagonal terms)
  AdjacencyList_type AM_2_AdjacencyList(bool* AMp,unsigned int n);
  //! Create an adjacency graph from a GraphvizDigraph
  AdjacencyList_type GraphvizDigraph_2_AdjacencyList(GraphvizDigraph& G1, set<int> select_index);
  //! Check if the graph contains any cycle (true if the model contains at least one cycle, false otherwise)
  bool has_cycle(vector<int> &circuit_stack, AdjacencyList_type& g);
  bool has_cycle_dfs(AdjacencyList_type& g, AdjacencyList_type::vertex_descriptor u, color_type& color, vector<int> &circuit_stack);
  //! Return the feedback set
  AdjacencyList_type Minimal_set_of_feedback_vertex(set<int> &feed_back_vertices, const AdjacencyList_type& G);
  //! Clear all in and out edges of vertex_to_eliminate and remove vertex_to_eliminate from the graph
  void Suppress(AdjacencyList_type::vertex_descriptor vertex_to_eliminate, AdjacencyList_type& G);
  void Suppress(int vertex_num, AdjacencyList_type& G);
  //! Reorder the recursive variables
  /*! They appear first in a quasi triangular form and they are followed by the feedback variables */
	void Reorder_the_recursive_variables(const AdjacencyList_type& G1, set<int> &feedback_vertices, vector< int> &Reordered_Vertices);
};

#endif // _MINIMUMFEEDBACKSET_HH
