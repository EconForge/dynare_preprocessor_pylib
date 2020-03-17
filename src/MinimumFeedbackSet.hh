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

#ifndef _MINIMUMFEEDBACKSET_HH
#define _MINIMUMFEEDBACKSET_HH

#include <map>
#include <vector>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#include <boost/graph/adjacency_list.hpp>
#pragma GCC diagnostic pop

using namespace std;

namespace MFS
{
  using VertexProperty_t = boost::property<boost::vertex_index_t, int,
                                           boost::property<boost::vertex_index1_t, int,
                                                           boost::property<boost::vertex_degree_t, int,
                                                                           boost::property<boost::vertex_in_degree_t, int,
                                                                                           boost::property<boost::vertex_out_degree_t, int >>>>>;
  using AdjacencyList_t = boost::adjacency_list<boost::listS, boost::listS, boost::bidirectionalS, VertexProperty_t>;

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
  //! Return the feedback set
  AdjacencyList_t Minimal_set_of_feedback_vertex(set<int> &feed_back_vertices, const AdjacencyList_t &G);
  //! Reorder the recursive variables
  /*! They appear first in a quasi triangular form and they are followed by the feedback variables */
  void Reorder_the_recursive_variables(const AdjacencyList_t &G1, set<int> &feedback_vertices, vector<int> &Reordered_Vertices);
};

#endif // _MINIMUMFEEDBACKSET_HH
