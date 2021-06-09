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
 * along with Dynare.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <algorithm>

#include "VariableDependencyGraph.hh"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#include <boost/graph/strong_components.hpp>
#include <boost/graph/topological_sort.hpp>
#pragma GCC diagnostic pop

using namespace boost;

VariableDependencyGraph::VariableDependencyGraph(int n) : base(n)
{
  /* It is necessary to manually initialize the vertex_index property since
     this graph uses listS and not vecS as underlying vertex container */
  auto v_index = get(vertex_index, *this);
  for (int i = 0; i < n; i++)
    put(v_index, vertex(i, *this), i);
}

void
VariableDependencyGraph::suppress(vertex_descriptor vertex_to_eliminate)
{
  clear_vertex(vertex_to_eliminate, *this);
  remove_vertex(vertex_to_eliminate, *this);
}

void
VariableDependencyGraph::suppress(int vertex_num)
{
  suppress(vertex(vertex_num, *this));
}

void
VariableDependencyGraph::eliminate(vertex_descriptor vertex_to_eliminate)
{
  if (in_degree(vertex_to_eliminate, *this) > 0 && out_degree(vertex_to_eliminate, *this) > 0)
    for (auto [it_in, in_end] = in_edges(vertex_to_eliminate, *this); it_in != in_end; ++it_in)
      for (auto [it_out, out_end] = out_edges(vertex_to_eliminate, *this); it_out != out_end; ++it_out)
        if (auto [ed, exist] = edge(source(*it_in, *this), target(*it_out, *this), *this);
            !exist)
          add_edge(source(*it_in, *this), target(*it_out, *this), *this);

  suppress(vertex_to_eliminate);
}

bool
VariableDependencyGraph::hasCycleDFS(vertex_descriptor u, color_t &color, vector<int> &circuit_stack) const
{
  auto v_index = get(vertex_index, *this);
  color[u] = gray_color;
  for (auto [vi, vi_end] = out_edges(u, *this); vi != vi_end; ++vi)
    if (color[target(*vi, *this)] == white_color
        && hasCycleDFS(target(*vi, *this), color, circuit_stack))
      {
        // cycle detected, return immediately
        circuit_stack.push_back(v_index[target(*vi, *this)]);
        return true;
      }
    else if (color[target(*vi, *this)] == gray_color)
      {
        // *vi is an ancestor!
        circuit_stack.push_back(v_index[target(*vi, *this)]);
        return true;
      }
  color[u] = black_color;
  return false;
}

bool
VariableDependencyGraph::hasCycle() const
{
  // Initialize color map to white
  color_t color;
  vector<int> circuit_stack;
  for (auto [vi, vi_end] = vertices(*this); vi != vi_end; ++vi)
    color[*vi] = white_color;

  // Perform depth-first search
  for (auto [vi, vi_end] = vertices(*this); vi != vi_end; ++vi)
    if (color[*vi] == white_color && hasCycleDFS(*vi, color, circuit_stack))
      return true;

  return false;
}

void
VariableDependencyGraph::print() const
{
  auto v_index = get(vertex_index, *this);
  cout << "Graph\n"
       << "-----\n";
  for (auto [it, it_end] = vertices(*this); it != it_end; ++it)
    {
      cout << "vertex[" << v_index[*it] + 1 << "] <-";
      for (auto [it_in, in_end] = in_edges(*it, *this); it_in != in_end; ++it_in)
        cout << v_index[source(*it_in, *this)] + 1 << " ";
      cout << "\n       ->";
      for (auto [it_out, out_end] = out_edges(*it, *this); it_out != out_end; ++it_out)
        cout << v_index[target(*it_out, *this)] + 1 << " ";
      cout << "\n";
    }
}

VariableDependencyGraph
VariableDependencyGraph::extractSubgraph(const vector<int> &select_index) const
{
  int n = select_index.size();
  VariableDependencyGraph G(n);
  auto v_index = get(vertex_index, *this);
  auto v_index1_G = get(vertex_index1, G); // Maps new vertices to original indices
  map<int, int> reverse_index; // Maps orig indices to new ones
  for (int i = 0; i < n; i++)
    {
      reverse_index[select_index[i]] = i;
      v_index1_G[vertex(i, G)] = select_index[i];
    }
  for (int i = 0; i < n; i++)
    {
      auto vi = vertex(select_index[i], *this);
      for (auto [it_out, out_end] = out_edges(vi, *this); it_out != out_end; ++it_out)
        if (auto it = reverse_index.find(v_index[target(*it_out, *this)]);
            it != reverse_index.end())
          add_edge(vertex(i, G), vertex(it->second, G), G);
    }
  return G;
}

bool
VariableDependencyGraph::vertexBelongsToAClique(vertex_descriptor vertex) const
{
  vector<vertex_descriptor> liste;
  bool agree = true;
  auto [it_in, in_end] = in_edges(vertex, *this);
  auto [it_out, out_end] = out_edges(vertex, *this);
  while (it_in != in_end && it_out != out_end && agree)
    {
      agree = (source(*it_in, *this) == target(*it_out, *this)
               && source(*it_in, *this) != target(*it_in, *this)); //not a loop
      liste.push_back(source(*it_in, *this));
      ++it_in;
      ++it_out;
    }
  if (agree)
    {
      if (it_in != in_end || it_out != out_end)
        agree = false;
      int i = 1;
      while (i < static_cast<int>(liste.size()) && agree)
        {
          int j = i + 1;
          while (j < static_cast<int>(liste.size()) && agree)
            {
              auto [ed1, exist1] = edge(liste[i], liste[j], *this);
              auto [ed2, exist2] = edge(liste[j], liste[i], *this);
              agree = exist1 && exist2;
              j++;
            }
          i++;
        }
    }
  return agree;
}

bool
VariableDependencyGraph::eliminationOfVerticesWithOneOrLessIndegreeOrOutdegree()
{
  bool something_has_been_done = false;
  bool not_a_loop;
  int i;
  vertex_iterator it, ita, it_end;
  for (tie(it, it_end) = vertices(*this), i = 0; it != it_end; ++it, i++)
    {
      int in_degree_n = in_degree(*it, *this);
      int out_degree_n = out_degree(*it, *this);
      if (in_degree_n <= 1 || out_degree_n <= 1)
        {
          not_a_loop = true;
          if (in_degree_n >= 1 && out_degree_n >= 1) // Do not eliminate a vertex if it loops on itself!
            for (auto [it_in, in_end] = in_edges(*it, *this); it_in != in_end; ++it_in)
              if (source(*it_in, *this) == target(*it_in, *this))
                not_a_loop = false;
          if (not_a_loop)
            {
              eliminate(*it);
              something_has_been_done = true;
              if (i > 0)
                it = ita;
              else
                {
                  tie(it, it_end) = vertices(*this);
                  i--;
                }
            }
        }
      ita = it;
    }
  return something_has_been_done;
}

bool
VariableDependencyGraph::eliminationOfVerticesBelongingToAClique()
{
  vertex_iterator it, ita, it_end;
  bool something_has_been_done = false;
  int i;
  for (tie(it, it_end) = vertices(*this), i = 0; it != it_end; ++it, i++)
    {
      if (vertexBelongsToAClique(*it))
        {
          eliminate(*it);
          something_has_been_done = true;
          if (i > 0)
            it = ita;
          else
            {
              tie(it, it_end) = vertices(*this);
              i--;
            }
        }
      ita = it;
    }
  return something_has_been_done;
}

bool
VariableDependencyGraph::suppressionOfVerticesWithLoop(set<int> &feed_back_vertices)
{
  bool something_has_been_done = false;
  vertex_iterator ita;
  int i = 0;
  for (auto [it, it_end] = vertices(*this); it != it_end; ++it, i++)
    {
      auto [ed, exist] = edge(*it, *it, *this);
      if (exist)
        {
          auto v_index = get(vertex_index, *this);
          feed_back_vertices.insert(v_index[*it]);
          suppress(*it);
          something_has_been_done = true;
          if (i > 0)
            it = ita;
          else
            {
              tie(it, it_end) = vertices(*this);
              i--;
            }
        }
      ita = it;
    }
  return something_has_been_done;
}

set<int>
VariableDependencyGraph::minimalSetOfFeedbackVertices() const
{
  set<int> feed_back_vertices;
  VariableDependencyGraph G(*this);
  while (num_vertices(G) > 0)
    {
      bool something_has_been_done = true;
      while (something_has_been_done && num_vertices(G) > 0)
        {
          something_has_been_done = G.eliminationOfVerticesWithOneOrLessIndegreeOrOutdegree();
          something_has_been_done = G.eliminationOfVerticesBelongingToAClique() || something_has_been_done;
          something_has_been_done = G.suppressionOfVerticesWithLoop(feed_back_vertices) || something_has_been_done;
        }

      if (!G.hasCycle())
        return feed_back_vertices;

      if (num_vertices(G) > 0)
        {
          /* If nothing has been done in the five previous rule then cut the
             vertex with the maximum in_degree+out_degree */
          int max_degree = 0, num = 0;
          vertex_iterator max_degree_index;
          for (auto [it, it_end] = vertices(G); it != it_end; ++it, num++)
            if (static_cast<int>(in_degree(*it, G) + out_degree(*it, G)) > max_degree)
              {
                max_degree = in_degree(*it, G) + out_degree(*it, G);
                max_degree_index = it;
              }

          auto v_index = get(vertex_index, G);
          feed_back_vertices.insert(v_index[*max_degree_index]);
          G.suppress(*max_degree_index);
        }
    }
  return feed_back_vertices;
}

vector<int>
VariableDependencyGraph::reorderRecursiveVariables(const set<int> &feedback_vertices) const
{
  vector<int> reordered_vertices;
  VariableDependencyGraph G(*this);
  auto v_index = get(vertex_index, G);

  // Suppress feedback vertices, in decreasing order
  for (auto it = feedback_vertices.rbegin(); it != feedback_vertices.rend(); ++it)
    G.suppress(*it);

  bool something_has_been_done = true;
  while (something_has_been_done)
    {
      something_has_been_done = false;
      vertex_iterator it, it_end, ita;
      int i;
      for (tie(it, it_end) = vertices(G), i = 0; it != it_end; ++it, i++)
        {
          if (in_degree(*it, G) == 0)
            {
              reordered_vertices.push_back(v_index[*it]);
              G.suppress(*it);
              something_has_been_done = true;
              if (i > 0)
                it = ita;
              else
                {
                  tie(it, it_end) = vertices(G);
                  i--;
                }
            }
          ita = it;
        }
    }

  if (num_vertices(G))
    cout << "Error in the computation of feedback vertex set\n";

  return reordered_vertices;
}

pair<int, vector<int>>
VariableDependencyGraph::sortedStronglyConnectedComponents() const
{
  vector<int> vertex2scc(num_vertices(*this));
  auto v_index = get(vertex_index, *this);

  // Compute SCCs and create mapping from vertices to unordered SCCs
  int num_scc = strong_components(static_cast<base>(*this), make_iterator_property_map(vertex2scc.begin(), v_index));

  // Create directed acyclic graph (DAG) associated to the SCCs
  adjacency_list<vecS, vecS, directedS> dag(num_scc);
  for (int i = 0; i < static_cast<int>(num_vertices(*this)); i++)
    {
      auto vi = vertex(i, *this);
      for (auto [it_out, out_end] = out_edges(vi, *this); it_out != out_end; ++it_out)
        if (int t_b = vertex2scc[v_index[target(*it_out, *this)]],
            s_b = vertex2scc[v_index[source(*it_out, *this)]];
            s_b != t_b)
          add_edge(s_b, t_b, dag);
    }

  /* Compute topological sort of DAG (ordered list of unordered SCC)
     Note: the order is reversed. */
  vector<int> reverseOrdered2unordered;
  topological_sort(dag, back_inserter(reverseOrdered2unordered));

  // Construct mapping from unordered SCC to ordered SCC
  vector<int> unordered2ordered(num_scc);
  for (int j = 0; j < num_scc; j++)
    unordered2ordered[reverseOrdered2unordered[num_scc-j-1]] = j;

  // Update the mapping of vertices to (now sorted) SCCs
  for (int i = 0; i < static_cast<int>(num_vertices(*this)); i++)
    vertex2scc[i] = unordered2ordered[vertex2scc[i]];

  return { num_scc, vertex2scc };
}
