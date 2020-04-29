/*
 * Copyright © 2009-2020 Dynare Team
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

#ifndef _VARIABLEDEPENDENCYGRAPH_HH
#define _VARIABLEDEPENDENCYGRAPH_HH

#include <map>
#include <vector>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#include <boost/graph/adjacency_list.hpp>
#pragma GCC diagnostic pop

using namespace std;

using VertexProperty_t = boost::property<boost::vertex_index_t, int,
                                           boost::property<boost::vertex_index1_t, int,
                                                           boost::property<boost::vertex_degree_t, int,
                                                                           boost::property<boost::vertex_in_degree_t, int,
                                                                                           boost::property<boost::vertex_out_degree_t, int>>>>>;

/* Class used to store a graph representing dependencies between variables.
   Used in the block decomposition. */
class VariableDependencyGraph
  : public boost::adjacency_list<boost::listS, boost::listS, boost::bidirectionalS, VertexProperty_t>
{
public:
  using color_t = map<boost::graph_traits<VariableDependencyGraph>::vertex_descriptor, boost::default_color_type>;
  using base = boost::adjacency_list<boost::listS, boost::listS, boost::bidirectionalS, VertexProperty_t>;

  VariableDependencyGraph(int n);
  //! Extracts a subgraph
  /*!
    \param[in] select_index The vertex indices to select
    \return The subgraph

    The property vertex_index1 of the subgraph contains indices of the original
    graph.
  */
  VariableDependencyGraph extractSubgraph(const vector<int> &select_index) const;
  //! Return the feedback set
  set<int> minimalSetOfFeedbackVertices() const;
  //! Reorder the recursive variables
  /*! They appear first in a quasi triangular form and they are followed by the feedback variables */
  vector<int> reorderRecursiveVariables(const set<int> &feedback_vertices) const;
  /* Computes the strongly connected components (SCCs) of the graph, and sort them
     topologically.
     Returns the number of SCCs, and a mapping of vertex indices to sorted SCC
     indices. */
  pair<int, vector<int>> sortedStronglyConnectedComponents() const;
  // Print on stdout a description of the graph
  void print() const;
private:
  // Remove a vertex (including all edges to and from it); takes a vertex descriptor
  void suppress(vertex_descriptor vertex_to_eliminate);
  // Remove a vertex (including all edges to and from it); takes a vertex index
  void suppress(int vertex_num);
  /* Remove a vertex, but keeping the paths that go through it (i.e. by adding
     edges that directly connect vertices that would otherwise be connected
     through the vertex to be removed) */
  void eliminate(vertex_descriptor vertex_to_eliminate);
  // Internal helper for hasCycle()
  bool hasCycleDFS(vertex_descriptor u, color_t &color, vector<int> &circuit_stack) const;
  // Determine whether the graph has a cycle
  bool hasCycle() const;
  bool vertexBelongsToAClique(vertex_descriptor vertex) const;
  bool eliminationOfVerticesWithOneOrLessIndegreeOrOutdegree();
  bool eliminationOfVerticesBelongingToAClique();
  // The suppressed vertices are stored in feedback set
  bool suppressionOfVerticesWithLoop(set<int> &feed_back_vertices);
};

#endif // _VARIABLEDEPENDENCYGRAPH_HH
