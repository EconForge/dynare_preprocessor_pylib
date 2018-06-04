/*
 * Copyright (C) 2009-2010 Dynare Team
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
#include <boost/graph/adjacency_list.hpp>

using namespace std;
using namespace boost;

namespace MFS
{
  typedef property<vertex_index_t, int,
                   property<vertex_index1_t, int,
                            property<vertex_degree_t, int,
                                     property<vertex_in_degree_t, int,
                                              property<vertex_out_degree_t, int >>>>> VertexProperty_t;
  typedef adjacency_list<listS, listS, bidirectionalS, VertexProperty_t> AdjacencyList_t;
  typedef map<graph_traits<AdjacencyList_t>::vertex_descriptor, default_color_type> color_t;
  using vector_vertex_descriptor_t = vector<AdjacencyList_t::vertex_descriptor>;

  //! Eliminate a vertex i
  /*! For a vertex i replace all edges e_k_i and e_i_j by a shorcut e_k_j and then Suppress the vertex i*/
  void Eliminate(AdjacencyList_t::vertex_descriptor vertex_to_eliminate, AdjacencyList_t &G);
  //! Collect all doublets (edges e_i_k such that there is an edge e_k_i with k!=i in the graph)
  /*! Returns the vector of doublets */
  vector_vertex_descriptor_t Collect_Doublet(AdjacencyList_t::vertex_descriptor vertex, AdjacencyList_t &G);
  //! Detect all the clique (all vertex in a clique are related to each other) in the graph
  bool Vertex_Belong_to_a_Clique(AdjacencyList_t::vertex_descriptor vertex, AdjacencyList_t &G);
  //! Graph reduction: eliminating purely intermediate variables or variables outside of any circuit
  bool Elimination_of_Vertex_With_One_or_Less_Indegree_or_Outdegree_Step(AdjacencyList_t &G);
  //! Graph reduction: elimination of a vertex inside a clique
  bool Elimination_of_Vertex_belonging_to_a_clique_Step(AdjacencyList_t &G);
  //! A vertex belong to the feedback vertex set if the vertex loops on itself.
  /*! We have to suppress this vertex and store it into the feedback set.*/
  bool Suppression_of_Vertex_X_if_it_loops_store_in_set_of_feedback_vertex_Step(set<int> &feed_back_vertices, AdjacencyList_t &G1);
  //! Print the Graph
  void Print(AdjacencyList_t &G);
  //! Create an adjacency graph from a Adjacency Matrix (an incidence Matrix without the diagonal terms)
  AdjacencyList_t AM_2_AdjacencyList(bool *AMp, unsigned int n);
  //! Extracts a subgraph
  /*!
    \param[in] G1 The original graph
    \param[in] select_index The vertex indices to select
    \return The subgraph

    The property vertex_index of the subgraph contains indices of the original
    graph, the property vertex_index1 contains new contiguous indices specific
    to the subgraph.
  */
  AdjacencyList_t extract_subgraph(AdjacencyList_t &G1, set<int> select_index);
  //! Check if the graph contains any cycle (true if the model contains at least one cycle, false otherwise)
  bool has_cycle(vector<int> &circuit_stack, AdjacencyList_t &g);
  bool has_cycle_dfs(AdjacencyList_t &g, AdjacencyList_t::vertex_descriptor u, color_t &color, vector<int> &circuit_stack);
  //! Return the feedback set
  AdjacencyList_t Minimal_set_of_feedback_vertex(set<int> &feed_back_vertices, const AdjacencyList_t &G);
  //! Clear all in and out edges of vertex_to_eliminate and remove vertex_to_eliminate from the graph
  void Suppress(AdjacencyList_t::vertex_descriptor vertex_to_eliminate, AdjacencyList_t &G);
  void Suppress(int vertex_num, AdjacencyList_t &G);
  //! Reorder the recursive variables
  /*! They appear first in a quasi triangular form and they are followed by the feedback variables */
  void Reorder_the_recursive_variables(const AdjacencyList_t &G1, set<int> &feedback_vertices, vector< int> &Reordered_Vertices);
};

#endif // _MINIMUMFEEDBACKSET_HH
