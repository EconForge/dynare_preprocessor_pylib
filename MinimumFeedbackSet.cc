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

#include "MinimumFeedbackSet.hh"

namespace MFS
  {
    void
    Suppress(AdjacencyList_type::vertex_descriptor vertex_to_eliminate, AdjacencyList_type& G)
    {
      /*clear all in and out edges of vertex_to_eliminate
      and remove vertex_to_eliminate from the graph*/
      clear_vertex(vertex_to_eliminate, G);
      remove_vertex(vertex_to_eliminate, G);
    }


    void
    Suppress(int vertex_num, AdjacencyList_type& G)
    {
      Suppress(vertex(vertex_num, G), G);
    }


    void
    Eliminate(AdjacencyList_type::vertex_descriptor vertex_to_eliminate, AdjacencyList_type& G)
    {
      /*before the vertex i suppression replace all edges e_k_i and e_i_j by e_k_j*/
      if (in_degree (vertex_to_eliminate, G) > 0 and out_degree (vertex_to_eliminate, G) > 0)
        {
          AdjacencyList_type::in_edge_iterator it_in, in_end;
          AdjacencyList_type::out_edge_iterator it_out, out_end;
          for (tie(it_in, in_end) = in_edges(vertex_to_eliminate, G); it_in != in_end; ++it_in)
            for (tie(it_out, out_end) = out_edges(vertex_to_eliminate, G); it_out != out_end; ++it_out)
              {
                AdjacencyList_type::edge_descriptor ed;
                bool exist;
                tie(ed, exist) = edge(source(*it_in, G) , target(*it_out, G), G);
                if (!exist)
                  add_edge(source(*it_in, G) , target(*it_out, G), G);
              }
        }
      Suppress(vertex_to_eliminate, G);
    }


    bool
    has_cycle_dfs(AdjacencyList_type& g, AdjacencyList_type::vertex_descriptor u, color_type& color, vector<int> &circuit_stack)
    {
      property_map<AdjacencyList_type, vertex_index_t>::type v_index = get(vertex_index, g);
      color[u] = gray_color;
      graph_traits<AdjacencyList_type>::out_edge_iterator vi, vi_end;
      for (tie(vi, vi_end) = out_edges(u, g); vi != vi_end; ++vi)
        if (color[target(*vi, g)] == white_color)
          {
            if (has_cycle_dfs(g, target(*vi, g), color, circuit_stack))
              {
                circuit_stack.push_back(v_index[target(*vi, g)]);
                return true; // cycle detected, return immediately
              }
          }
        else if (color[target(*vi, g)] == gray_color) // *vi is an ancestor!
          {
            circuit_stack.push_back(v_index[target(*vi, g)]);
            return true;
          }
      color[u] = black_color;
      return false;
    }

    bool
    has_cylce(AdjacencyList_type& g, vector<int> &circuit_stack)
    {
      color_type color;
      graph_traits<AdjacencyList_type>::vertex_iterator vi, vi_end;
      for (tie(vi, vi_end) = vertices(g); vi != vi_end; vi++)
        color[*vi] = white_color;
      property_map<AdjacencyList_type, vertex_index_t>::type v_index = get(vertex_index, g);
      for (tie(vi, vi_end) = vertices(g); vi != vi_end; vi++)
        if (color[*vi] == white_color)
          if (has_cycle_dfs(g, *vi, color, circuit_stack))
            return true;
      return false;
    }

    bool
    has_cycle(vector<int> &circuit_stack, AdjacencyList_type& G)
    {
      return has_cylce(G, circuit_stack);
    }


    void
    Print(AdjacencyList_type& G)
    {
      AdjacencyList_type::vertex_iterator  it, it_end, it_begin;
      property_map<AdjacencyList_type, vertex_index_t>::type v_index = get(vertex_index, G);
      cout << "Graph\n";
      cout << "-----\n";
      for (tie(it, it_end) = vertices(G);it != it_end; ++it)
        {
          cout << "vertex[" << v_index[*it] + 1 << "] <-";
          AdjacencyList_type::in_edge_iterator it_in, in_end;
          for (tie(it_in, in_end) = in_edges(*it, G); it_in != in_end; ++it_in)
            cout << v_index[source(*it_in, G)] + 1 << " ";
          cout << "\n       ->";
          AdjacencyList_type::out_edge_iterator it_out, out_end;
          for (tie(it_out, out_end) = out_edges(*it, G); it_out != out_end; ++it_out)
            cout << v_index[target(*it_out, G)] + 1 << " ";
          cout << "\n";
        }
    }


    AdjacencyList_type
    AM_2_AdjacencyList(bool* AM, unsigned int n)
    {
      AdjacencyList_type G(n);
      property_map<AdjacencyList_type, vertex_index_t>::type v_index = get(vertex_index, G);
      property_map<AdjacencyList_type, vertex_index1_t>::type v_index1 = get(vertex_index1, G);
      for (unsigned int i = 0;i < n;i++)
        {
          put(v_index, vertex(i, G), i);
          put(v_index1, vertex(i, G), i);
        }
      for (unsigned int i = 0;i < n;i++)
        for (unsigned int j = 0;j < n;j++)
          if (AM[i*n+j])
            add_edge(vertex(j, G), vertex(i, G), G);
      return G;
    }


    void
    Print(GraphvizDigraph& G)
    {
      GraphvizDigraph::vertex_iterator  it, it_end, it_begin;
      property_map<GraphvizDigraph, vertex_index_t>::type v_index = get(vertex_index, G);
      cout << "Graph\n";
      cout << "-----\n";
      for (tie(it, it_end) = vertices(G);it != it_end; ++it)
        {
          cout << "vertex[" << v_index[*it] + 1 << "] ->";
          GraphvizDigraph::out_edge_iterator it_out, out_end;
          for (tie(it_out, out_end) = out_edges(*it, G); it_out != out_end; ++it_out)
            cout << v_index[target(*it_out, G)] + 1 << " ";
          cout << "\n";
        }
    }


    GraphvizDigraph
    AM_2_GraphvizDigraph(bool* AM, unsigned int n)
    {
      GraphvizDigraph G(n);
      property_map<GraphvizDigraph, vertex_index_t>::type v_index = get(vertex_index, G);
      /*for (unsigned int i = 0;i < n;i++)
        cout << "v_index[" << i << "] = " << v_index[i] << "\n";*/
      //put(v_index, vertex(i, G), i);
      //v_index[/*vertex(i,G)*/i]["v_index"]=i;
      for (unsigned int i = 0;i < n;i++)
        for (unsigned int j = 0;j < n;j++)
          if (AM[i*n+j])
            add_edge(vertex(j, G), vertex(i, G), G);
      return G;
    }


    AdjacencyList_type
    GraphvizDigraph_2_AdjacencyList(GraphvizDigraph& G1, set<int> select_index)
    {
      unsigned int n = select_index.size();
      AdjacencyList_type G(n);
      property_map<AdjacencyList_type, vertex_index_t>::type v_index = get(vertex_index, G);
      property_map<AdjacencyList_type, vertex_index1_t>::type v_index1 = get(vertex_index1, G);
      property_map<GraphvizDigraph, vertex_index_t>::type v1_index = get(vertex_index, G1);
      set<int>::iterator it = select_index.begin();
      map<int,int> reverse_index;
      for (unsigned int i = 0;i < n;i++, it++)
        {
          reverse_index[v1_index[*it]]=i;
          put(v_index, vertex(i, G), v1_index[*it]);
          put(v_index1, vertex(i, G), i);
        }
      unsigned int i;
      for (it = select_index.begin(), i = 0;i < n;i++, it++)
        {
          GraphvizDigraph::out_edge_iterator it_out, out_end;
          GraphvizDigraph::vertex_descriptor vi = vertex(*it, G1);
          for (tie(it_out, out_end) = out_edges(vi, G1); it_out != out_end; ++it_out)
            {
              int ii = target(*it_out, G1);
              if (select_index.find(ii) != select_index.end())
                add_edge( vertex(reverse_index[source(*it_out, G1)],G), vertex(reverse_index[target(*it_out, G1)], G), G);
            }
        }
      return G;
    }

    vector_vertex_descriptor
    Collect_Doublet(AdjacencyList_type::vertex_descriptor vertex, AdjacencyList_type& G)
    {
      /*collect all doublet (for each edge e_i_k there is an edge e_k_i with k!=i) in the graph
        and return the vector of doublet*/
      AdjacencyList_type::in_edge_iterator it_in, in_end;
      AdjacencyList_type::out_edge_iterator it_out, out_end;
      vector<AdjacencyList_type::vertex_descriptor> Doublet;
      if (in_degree(vertex, G) > 0 and out_degree(vertex, G) > 0)
        for (tie(it_in, in_end) = in_edges(vertex, G); it_in != in_end; ++it_in)
          for (tie(it_out, out_end) = out_edges(vertex, G); it_out != out_end; ++it_out)
            if (source(*it_in, G) == target(*it_out, G) and source(*it_in, G) != target(*it_in, G))  // not a loop
              Doublet.push_back(source(*it_in, G));
      return Doublet;
    }

    bool
    Vertex_Belong_to_a_Clique(AdjacencyList_type::vertex_descriptor vertex, AdjacencyList_type& G)
    {
      /*Detect all the clique (all vertex in a clique are related to each other) in the graph*/
      vector<AdjacencyList_type::vertex_descriptor> liste;
      bool agree = true;
      AdjacencyList_type::in_edge_iterator it_in, in_end;
      AdjacencyList_type::out_edge_iterator it_out, out_end;
      tie(it_in, in_end) = in_edges(vertex, G);
      tie(it_out, out_end) = out_edges(vertex, G);
      while (it_in != in_end and it_out != out_end and agree)
        {
          agree = (source(*it_in, G) == target(*it_out, G) and source(*it_in, G) != target(*it_in, G));  //not a loop
          liste.push_back(source(*it_in, G));
          it_in++;
          it_out++;
        }
      if (agree)
        {
          if (it_in != in_end or it_out != out_end)
            agree = false;
          unsigned int i = 1;
          while (i < liste.size() and agree)
            {
              unsigned int j = i + 1;
              while (j < liste.size() and agree)
                {
                  AdjacencyList_type::edge_descriptor ed;
                  bool exist1, exist2;
                  tie(ed, exist1) = edge(liste[i], liste[j] , G);
                  tie(ed, exist2) = edge(liste[j], liste[i] , G);
                  agree = (exist1 and exist2);
                  j++;
                }
              i++;
            }
        }
      return agree;
    }




    bool
    Elimination_of_Vertex_With_One_or_Less_Indegree_or_Outdegree_Step(AdjacencyList_type& G)
    {
      /*Graph reduction: eliminating purely intermediate variables or variables outside of any circuit*/
      bool something_has_been_done = false;
      bool not_a_loop;
      int i;
      AdjacencyList_type::vertex_iterator  it, it1, ita, it_end, it_begin;
      tie(it, it_end) = vertices(G);
      it_begin = it;
      property_map<AdjacencyList_type, vertex_index_t>::type v_index = get(vertex_index, G);
      for ( i = 0; it != it_end; ++it, i++)
        {
          int in_degree_n = in_degree(*it, G);
          int out_degree_n = out_degree(*it, G);
          if (in_degree_n <= 1 or out_degree_n <= 1)
            {
              not_a_loop = true;
              if (in_degree_n >= 1 and out_degree_n >= 1) //do not eliminate a vertex if it loops on its self!
                {
                  AdjacencyList_type::in_edge_iterator it_in, in_end;
                  for (tie(it_in, in_end) = in_edges(*it, G); it_in != in_end; ++it_in)
                    if (source(*it_in, G) == target(*it_in, G))
                      {
#ifdef verbose
                        cout << v_index[source(*it_in, G)] << " == " << v_index[target(*it_in, G)] << "\n";
#endif
                        not_a_loop = false;
                      }
                }
              if (not_a_loop)
                {
#ifdef verbose
                  property_map<AdjacencyList_type, vertex_index_t>::type v_index = get(vertex_index, G);
                  cout << "->eliminate vertex[" << v_index[*it] + 1 << "]\n";
#endif
                  Eliminate(*it, G);
#ifdef verbose
                  Print(G);
#endif
                  something_has_been_done = true;
                  if (i > 0)
                    it = ita;
                  else
                    {
                      tie(it, it_end) = vertices(G);
                      i--;
                    }
                }
            }
          ita = it;
        }
      return something_has_been_done;
    }


    bool
    Elimination_of_Vertex_belonging_to_a_clique_Step(AdjacencyList_type& G)
    {
      /*Graphe reduction: elimination of a vertex inside a clique*/
      AdjacencyList_type::vertex_iterator  it, it1, ita, it_end, it_begin;
      bool something_has_been_done = false;
      int i;
      tie(it, it_end) = vertices(G);
      it_begin = it;
      for (i = 0;it != it_end; ++it, i++)
        {
          if (Vertex_Belong_to_a_Clique(*it, G))
            {
#ifdef verbose
              property_map<AdjacencyList_type, vertex_index_t>::type v_index = get(vertex_index, G);
              cout << "eliminate vertex[" << v_index[*it] + 1 << "]\n";
#endif
              Eliminate(*it, G);
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
      return something_has_been_done;
    }




    bool
    Suppression_of_Vertex_X_if_it_loops_store_in_set_of_feedback_vertex_Step(set<int> &feed_back_vertices, AdjacencyList_type& G)
    {
      /*If a vertex loop on itself it's a feedback variable
        we eliminate it from the graph and store the vertex
        in the minimum feedback set*/
      bool something_has_been_done = false;
      AdjacencyList_type::vertex_iterator  it, it_end, it_begin, ita;
      int i = 0;
      tie(it, it_end) = vertices(G);
      it_begin = it;
      for (;it != it_end; ++it, i++)
        {
          AdjacencyList_type::edge_descriptor ed;
          bool exist;
          tie(ed, exist) = edge(*it, *it , G);
          if (exist)
            {
#ifdef verbose
              property_map<AdjacencyList_type, vertex_index_t>::type v_index = get(vertex_index, G);
              cout << "store v[*it] = " << v_index[*it]+1 << "\n";
#endif
              property_map<AdjacencyList_type, vertex_index1_t>::type v_index1 = get(vertex_index1, G);
              feed_back_vertices.insert(v_index1[*it] );
              /*property_map<AdjacencyList_type, vertex_index_t>::type v_index = get(vertex_index, G);
              feed_back_vertices.insert(v_index[*it] );*/
              Suppress(*it, G);
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
      return something_has_been_done;
    }


    AdjacencyList_type
    Minimal_set_of_feedback_vertex(set<int> &feed_back_vertices,const AdjacencyList_type& G1)
    {
      bool something_has_been_done = true;
      int cut_ = 0;
      AdjacencyList_type G(G1);
      while (num_vertices(G) > 0)
        {
          while (something_has_been_done and num_vertices(G) > 0)
            {
              //Rule 1
              something_has_been_done = (Elimination_of_Vertex_With_One_or_Less_Indegree_or_Outdegree_Step(G) /*or something_has_been_done*/);
#ifdef verbose
              cout << "1 something_has_been_done=" << something_has_been_done << "\n";
#endif

              //Rule 2
              something_has_been_done = (Elimination_of_Vertex_belonging_to_a_clique_Step(G) or something_has_been_done);
#ifdef verbose
              cout << "2 something_has_been_done=" << something_has_been_done << "\n";
#endif

              //Rule 3
              something_has_been_done = (Suppression_of_Vertex_X_if_it_loops_store_in_set_of_feedback_vertex_Step(feed_back_vertices, G) or something_has_been_done);
#ifdef verbose
              cout << "3 something_has_been_done=" << something_has_been_done << "\n";
#endif
            }
          vector<int> circuit;
          if (!has_cycle(circuit, G))
            {
#ifdef verobse
              cout << "has_cycle=false\n";
#endif
              //sort(feed_back_vertices.begin(), feed_back_vertices.end());
              return G;
            }
          if (num_vertices(G) > 0)
            {
              /*if nothing has been done in the five previous rule then cut the vertex with the maximum in_degree+out_degree*/
              unsigned int max_degree = 0, num = 0;
              AdjacencyList_type::vertex_iterator  it, it_end, max_degree_index;
              for (tie(it, it_end) = vertices(G);it != it_end; ++it, num++)
                {
                  if (in_degree(*it, G) + out_degree(*it, G) > max_degree)
                    {
                      max_degree = in_degree(*it, G) + out_degree(*it, G);
                      max_degree_index = it;
                    }
                }
              property_map<AdjacencyList_type, vertex_index1_t>::type v_index1 = get(vertex_index1, G);
              feed_back_vertices.insert(v_index1[*max_degree_index]);
              /*property_map<AdjacencyList_type, vertex_index_t>::type v_index = get(vertex_index, G);
              feed_back_vertices.insert(v_index[*max_degree_index]);*/
              //cout << "v_index1[*max_degree_index] = " << v_index1[*max_degree_index] << "\n";
              cut_++;
#ifdef verbose
              property_map<AdjacencyList_type, vertex_index_t>::type v_index = get(vertex_index, G);
              cout << "--> cut vertex " << v_index[*max_degree_index] + 1 << "\n";
#endif
              Suppress(*max_degree_index, G);
              something_has_been_done = true;
            }
        }
#ifdef verbose
      cout << "cut_=" << cut_ << "\n";
#endif
      //sort(feed_back_vertices.begin(), feed_back_vertices.end());
      return G;
    }

    struct rev
      {
        bool operator()(const int a, const int b) const
          {
            return (a>b);
          }
      };

    vector<int>
    Reorder_the_recursive_variables(const AdjacencyList_type& G1, set<int> &feedback_vertices)
    {
      AdjacencyList_type G(G1);
      property_map<AdjacencyList_type, vertex_index_t>::type v_index = get(vertex_index, G);
      set<int>::iterator its, ita;
      set<int, rev> fv;
      for (its = feedback_vertices.begin(); its != feedback_vertices.end(); its++)
        fv.insert(*its);
      int i=0;
      for (its = fv.begin(); its != fv.end(); its++, i++)
        {
          //cout << "supress " << v_index[vertex(*its, G)]+1 << "  " << *its << "\n";
          Suppress(*its, G);
        }
      vector< int> Reordered_Vertices;
      bool something_has_been_done = true;
      while (something_has_been_done)
        {
          something_has_been_done = false;
          AdjacencyList_type::vertex_iterator  it, it_end, it_begin, ita;
          tie(it, it_end) = vertices(G);
          int i = 0;
          for (it_begin = it;it != it_end; ++it, i++)
            {
              if (in_degree(*it, G) == 0)
                {
                  Reordered_Vertices.push_back(v_index[*it]);
                  Suppress(*it, G);
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
      return Reordered_Vertices;
    }

  }


