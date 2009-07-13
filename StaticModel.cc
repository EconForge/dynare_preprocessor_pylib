/*
 * Copyright (C) 2003-2009 Dynare Team
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
#include <cassert>

#include <deque>
#include <algorithm>
#include <iterator>
#include <functional>

#ifdef DEBUG
// For select2nd()
# ifdef __GNUC__
#  include <ext/functional>
using namespace __gnu_cxx;
# endif
#endif

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/max_cardinality_matching.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/topological_sort.hpp>

#include "StaticModel.hh"
#include "MinimumFeedbackSet.hh"

using namespace boost;

StaticModel::StaticModel(SymbolTable &symbol_table_arg,
                         NumericalConstants &num_constants_arg) :
  ModelTree(symbol_table_arg, num_constants_arg),
  block_mfs(false)
{
}

void
StaticModel::writeStaticMFile(ostream &output, const string &func_name) const
{
  // Writing comments and function definition command
  output << "function [residual, g1, g2] = " << func_name << "(y, x, params)" << endl
         << "%" << endl
         << "% Status : Computes static model for Dynare" << endl
         << "%" << endl
         << "% Warning : this file is generated automatically by Dynare" << endl
         << "%           from model file (.mod)" << endl
         << endl
         << "residual = zeros( " << equations.size() << ", 1);" << endl
         << endl
         << "%" << endl
         << "% Model equations" << endl
         << "%" << endl
         << endl;

  writeModelLocalVariables(output, oMatlabStaticModel);

  writeTemporaryTerms(temporary_terms, output, oMatlabStaticModel);

  writeModelEquations(output, oMatlabStaticModel);

  output << "if ~isreal(residual)" << endl
         << "  residual = real(residual)+imag(residual).^2;" << endl
         << "end" << endl
         << "if nargout >= 2," << endl
         << "  g1 = zeros(" << equations.size() << ", " << symbol_table.endo_nbr() << ");" << endl
         << endl
         << "%" << endl
         << "% Jacobian matrix" << endl
         << "%" << endl
         << endl;

  // Write Jacobian w.r. to endogenous only
  for (first_derivatives_type::const_iterator it = first_derivatives.begin();
       it != first_derivatives.end(); it++)
    {
      int eq = it->first.first;
      int symb_id = it->first.second;
      NodeID d1 = it->second;

      ostringstream g1;
      g1 << "  g1(" << eq+1 << "," << symbol_table.getTypeSpecificID(symb_id)+1 << ")";
      output << g1.str() << "=" << g1.str() << "+";
      d1->writeOutput(output, oMatlabStaticModel, temporary_terms);
      output << ";" << endl;
    }

  output << "  if ~isreal(g1)" << endl
         << "    g1 = real(g1)+2*imag(g1);" << endl
         << "  end" << endl
         << "end" << endl
         << "if nargout >= 3," << endl
         << "%" << endl
         << "% Hessian matrix" << endl
         << "%" << endl
         << endl;

  int g2ncols = symbol_table.endo_nbr() * symbol_table.endo_nbr();
  if (second_derivatives.size())
    {
      output << "  v2 = zeros(" << NNZDerivatives[1] << ",3);" << endl;

      // Write Hessian w.r. to endogenous only (only if 2nd order derivatives have been computed)
      int k = 0; // Keep the line of a 2nd derivative in v2
      for (second_derivatives_type::const_iterator it = second_derivatives.begin();
           it != second_derivatives.end(); it++)
        {
          int eq = it->first.first;
          int symb_id1 = it->first.second.first;
          int symb_id2 = it->first.second.second;
          NodeID d2 = it->second;

          int tsid1 = symbol_table.getTypeSpecificID(symb_id1);
          int tsid2 = symbol_table.getTypeSpecificID(symb_id2);

          int col_nb = tsid1*symbol_table.endo_nbr()+tsid2;
          int col_nb_sym = tsid2*symbol_table.endo_nbr()+tsid1;

          output << "v2(" << k+1 << ",1)=" << eq + 1 << ";" << endl
                 << "v2(" << k+1 << ",2)=" << col_nb + 1 << ";" << endl
                 << "v2(" << k+1 << ",3)=";
          d2->writeOutput(output, oMatlabStaticModel, temporary_terms);
          output << ";" << endl;

          k++;

          // Treating symetric elements
          if (symb_id1 != symb_id2)
            {
              output << "v2(" << k+1 << ",1)=" << eq + 1 << ";" << endl
                     << "v2(" << k+1 << ",2)=" << col_nb_sym + 1 << ";" << endl
                     << "v2(" << k+1 << ",3)=v2(" << k << ",3);" << endl;
              k++;
            }
        }

      output << "  g2 = sparse(v2(:,1),v2(:,2),v2(:,3)," << equations.size() << "," << g2ncols << ");" << endl;
    }
  else // Either hessian is all zero, or we didn't compute it
    output << "  g2 = sparse([],[],[]," << equations.size() << "," << g2ncols << ");" << endl;

  output << "end;" << endl; // Close the if nargout >= 3 statement
}

void
StaticModel::writeStaticFile(const string &basename) const
{
  string filename = basename + "_static.m";

  ofstream output;
  output.open(filename.c_str(), ios::out | ios::binary);
  if (!output.is_open())
    {
      cerr << "ERROR: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }

  if (block_mfs)
    writeStaticBlockMFSFile(output, basename + "_static");
  else
    writeStaticMFile(output, basename + "_static");

  output.close();
}

void
StaticModel::computingPass(bool block_mfs_arg, bool hessian, bool no_tmp_terms)
{
  block_mfs = block_mfs_arg;

  // Compute derivatives w.r. to all endogenous
  set<int> vars;
  for(int i = 0; i < symbol_table.endo_nbr(); i++)
    vars.insert(symbol_table.getID(eEndogenous, i));

  // Launch computations
  cout << "Computing static model derivatives:" << endl
       << " - order 1" << endl;
  computeJacobian(vars);

  if (hessian)
    {
      cout << " - order 2" << endl;
      computeHessian(vars);
    }

  if (block_mfs)
    {
      computeNormalization();
      computeSortedBlockDecomposition();
      computeMFS();
      computeSortedRecursive();
      computeBlockMFSJacobian();
    }
  else if (!no_tmp_terms)
    computeTemporaryTerms(true);
}

int
StaticModel::computeDerivID(int symb_id, int lag)
{
  if (symbol_table.getType(symb_id) == eEndogenous)
    return symb_id;
  else
    return -1;
}

int
StaticModel::getDerivID(int symb_id, int lag) const throw (UnknownDerivIDException)
{
  if (symbol_table.getType(symb_id) == eEndogenous)
    return symb_id;
  else
    throw UnknownDerivIDException();
}

void
StaticModel::computeNormalization()
{
  const int n = equation_number();

  assert(n == symbol_table.endo_nbr());

  typedef adjacency_list<vecS, vecS, undirectedS> BipartiteGraph;

  /*
    Vertices 0 to n-1 are for endogenous (using type specific ID)
    Vertices n to 2*n-1 are for equations (using equation no.)
  */
  BipartiteGraph g(2 * n);

  // Fill in the graph
  set<pair<int, int> > endo;
  for(int i = 0; i < n; i++)
    {
      endo.clear();
      equations[i]->collectEndogenous(endo);
      for(set<pair<int, int> >::const_iterator it = endo.begin();
          it != endo.end(); it++)
        add_edge(i + n, it->first, g);
    }

  // Compute maximum cardinality matching
  vector<int> mate_map(2*n);

#if 1
  bool check = checked_edmonds_maximum_cardinality_matching(g, &mate_map[0]);
#else // Alternative way to compute normalization, by giving an initial matching using natural normalizations
  fill(mate_map.begin(), mate_map.end(), graph_traits<BipartiteGraph>::null_vertex());

  multimap<int, int> natural_endo2eqs;
  computeNormalizedEquations(natural_endo2eqs);

  for(int i = 0; i < symbol_table.endo_nbr(); i++)
    {
      if (natural_endo2eqs.count(i) == 0)
        continue;

      int j = natural_endo2eqs.find(i)->second;

      put(&mate_map[0], i, n+j);
      put(&mate_map[0], n+j, i);
    }

  edmonds_augmenting_path_finder<BipartiteGraph, size_t *, property_map<BipartiteGraph, vertex_index_t>::type> augmentor(g, &mate_map[0], get(vertex_index, g));
  bool not_maximum_yet = true;
  while(not_maximum_yet)
    {
      not_maximum_yet = augmentor.augment_matching();
    }
  augmentor.get_current_matching(&mate_map[0]);

  bool check = maximum_cardinality_matching_verifier<BipartiteGraph, size_t *, property_map<BipartiteGraph, vertex_index_t>::type>::verify_matching(g, &mate_map[0], get(vertex_index, g));
#endif

  assert(check);

  // Check if all variables are normalized
  vector<int>::const_iterator it = find(mate_map.begin(), mate_map.begin() + n, graph_traits<BipartiteGraph>::null_vertex());
  if (it != mate_map.begin() + n)
    {
      cerr << "ERROR: Could not normalize static model. Variable "
           << symbol_table.getName(symbol_table.getID(eEndogenous, it - mate_map.begin()))
           << " is not in the maximum cardinality matching." << endl;
      exit(EXIT_FAILURE);
    }

#ifdef DEBUG
  for(int i = 0; i < n; i++)
    cout << "Endogenous " << symbol_table.getName(symbol_table.getID(eEndogenous, i))
         << " matched with equation " << (mate_map[i]-n+1) << endl;
#endif

  // Create the resulting map, by copying the n first elements of mate_map, and substracting n to them
  endo2eq.resize(equation_number());
  transform(mate_map.begin(), mate_map.begin() + n, endo2eq.begin(), bind2nd(minus<int>(), n));

#ifdef DEBUG
  multimap<int, int> natural_endo2eqs;
  computeNormalizedEquations(natural_endo2eqs);

  int n1 = 0, n2 = 0;

  for(int i = 0; i < symbol_table.endo_nbr(); i++)
    {
      if (natural_endo2eqs.count(i) == 0)
        continue;

      n1++;

      pair<multimap<int, int>::const_iterator, multimap<int, int>::const_iterator> x = natural_endo2eqs.equal_range(i);
      if (find_if(x.first, x.second, compose1(bind2nd(equal_to<int>(), endo2eq[i]), select2nd<multimap<int, int>::value_type>())) == x.second)
        cout << "Natural normalization of variable " << symbol_table.getName(symbol_table.getID(eEndogenous, i))
             << " not used." << endl;
      else
        n2++;
    }

  cout << "Used " << n2 << " natural normalizations out of " << n1 << ", for a total of " << n << " equations." << endl;
#endif
}

void
StaticModel::computeNormalizedEquations(multimap<int, int> &endo2eqs) const
{
  for(int i = 0; i < equation_number(); i++)
    {
      VariableNode *lhs = dynamic_cast<VariableNode *>(equations[i]->get_arg1());
      if (lhs == NULL)
        continue;

      int symb_id = lhs->get_symb_id();
      if (symbol_table.getType(symb_id) != eEndogenous)
        continue;

      set<pair<int, int> > endo;
      equations[i]->get_arg2()->collectEndogenous(endo);
      if (endo.find(make_pair(symbol_table.getTypeSpecificID(symb_id), 0)) != endo.end())
        continue;

      endo2eqs.insert(make_pair(symbol_table.getTypeSpecificID(symb_id), i));
      cout << "Endogenous " << symbol_table.getName(symb_id) << " normalized in equation " << (i+1) << endl;
    }
}

void
StaticModel::writeLatexFile(const string &basename) const
{
  writeLatexModelFile(basename + "_static.tex", oLatexStaticModel);
}

void
StaticModel::computeSortedBlockDecomposition()
{
  const int n = equation_number();

  assert((int) endo2eq.size() == n);

  // Compute graph representation of static model
  typedef adjacency_list<vecS, vecS, directedS> DirectedGraph;
  DirectedGraph g(n);

  set<pair<int, int> > endo;
  for(int i = 0; i < n; i++)
    {
      endo.clear();
      equations[endo2eq[i]]->collectEndogenous(endo);
      for(set<pair<int, int> >::const_iterator it = endo.begin();
          it != endo.end(); it++)
        add_edge(it->first, i, g);
    }

  // Compute strongly connected components
  vector<int> endo2block(n);
  int m = strong_components(g, &endo2block[0]);

  // Create directed acyclic graph associated to the strongly connected components
  DirectedGraph dag(m);
  graph_traits<DirectedGraph>::edge_iterator ei, ei_end;
  for(tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
    {
      int s = endo2block[source(*ei, g)];
      int t = endo2block[target(*ei, g)];
      if (s != t)
        add_edge(s, t, dag);
    }

  // Compute topological sort of DAG (ordered list of unordered SCC)
  deque<int> ordered2unordered;
  topological_sort(dag, front_inserter(ordered2unordered)); // We use a front inserter because topological_sort returns the inverse order
  // Construct mapping from unordered SCC to ordered SCC
  vector<int> unordered2ordered(m);
  for(int i = 0; i < m; i++)
    unordered2ordered[ordered2unordered[i]] = i;

  // Fill in data structure representing blocks
  blocks.clear();
  blocks.resize(m);
  for(int i = 0; i < n; i++)
    blocks[unordered2ordered[endo2block[i]]].insert(i);

#ifdef DEBUG
  cout << "Found " << m << " blocks" << endl;
  for(int i = 0; i < m; i++)
    cout << " Block " << i << " of size " << blocks[i].size() << endl;
#endif
}

void
StaticModel::computeMFS()
{
  const int n = equation_number();
  assert((int) endo2eq.size() == n);

  const int nblocks = blocks.size();
  blocksMFS.clear();
  blocksMFS.resize(nblocks);

  // Iterate over blocks
  for(int b = 0; b < nblocks; b++)
    {
      // Construct subgraph for MFS computation, where vertex number is position in the block
      int p = blocks[b].size();
      MFS::AdjacencyList_type g(p);

      // Construct v_index and v_index1 properties, and a mapping between type specific IDs and vertex descriptors
      property_map<MFS::AdjacencyList_type, vertex_index_t>::type v_index = get(vertex_index, g);
      property_map<MFS::AdjacencyList_type, vertex_index1_t>::type v_index1 = get(vertex_index1, g);
      map<int, graph_traits<MFS::AdjacencyList_type>::vertex_descriptor> tsid2vertex;
      int j = 0;
      for(set<int>::const_iterator it = blocks[b].begin(); it != blocks[b].end(); ++it)
        {
          tsid2vertex[*it] = vertex(j, g);
          put(v_index, vertex(j, g), *it);
          put(v_index1, vertex(j, g), *it);
          j++;
        }

      // Add edges, loop over endogenous in the block
      set<pair<int, int> > endo;
      for(set<int>::const_iterator it = blocks[b].begin(); it != blocks[b].end(); ++it)
        {
          endo.clear();

          // Test if associated equation is in normalized form, and compute set of endogenous appearing in it
          ExprNode *lhs = equations[endo2eq[*it]]->get_arg1();
          VariableNode *lhs_var = dynamic_cast<VariableNode *>(lhs);
          if (lhs_var == NULL || lhs_var->get_symb_id() != symbol_table.getID(eEndogenous, *it))
            lhs->collectEndogenous(endo); // Only collect endogenous of LHS if not normalized form
          ExprNode *rhs = equations[endo2eq[*it]]->get_arg2();
          rhs->collectEndogenous(endo);

          for(set<pair<int, int> >::const_iterator it2 = endo.begin();
              it2 != endo.end(); ++it2)
            if (blocks[b].find(it2->first) != blocks[b].end()) // Add edge only if vertex member of this block
              add_edge(tsid2vertex[it2->first], tsid2vertex[*it], g);
        }

      // Compute minimum feedback set
      MFS::Minimal_set_of_feedback_vertex(blocksMFS[b], g);

#ifdef DEBUG
      cout << "Block " << b << ": " << blocksMFS[b].size() << "/" << blocks[b].size() << " in MFS" << endl;
#endif
    }
}

void
StaticModel::computeSortedRecursive()
{
  const int nblocks = blocks.size();
  blocksRecursive.clear();
  blocksRecursive.resize(nblocks);

  for(int b = 0; b < nblocks; b++)
    {
      // Construct the set of recursive vars
      // The index in this vector will be the vertex descriptor in the graph
      vector<int> recurs_vars;
      set_difference(blocks[b].begin(), blocks[b].end(),
                     blocksMFS[b].begin(), blocksMFS[b].end(),
                     back_inserter(recurs_vars));

      // Construct graph representation of recursive vars
      typedef adjacency_list<vecS, vecS, directedS> DirectedGraph;
      DirectedGraph dag(recurs_vars.size());
      set<pair<int, int> > endo;
      for(int i = 0; i < (int) recurs_vars.size(); i++)
        {
          endo.clear();
          equations[endo2eq[recurs_vars[i]]]->get_arg2()->collectEndogenous(endo);
          for(set<pair<int, int> >::const_iterator it = endo.begin();
              it != endo.end(); it++)
            {
              vector<int>::const_iterator it2 = find(recurs_vars.begin(), recurs_vars.end(), it->first);
              if (it2 != recurs_vars.end())
                {
                  int source_vertex = it2 - recurs_vars.begin();
                  add_edge(source_vertex, i, dag);
                }
            }
        }
      // Compute topological sort
      deque<int> ordered_recurs_vertices;
      topological_sort(dag, front_inserter(ordered_recurs_vertices)); // We use a front inserter because topological_sort returns the inverse order

      // Construct the final order
      for(deque<int>::const_iterator it = ordered_recurs_vertices.begin();
          it != ordered_recurs_vertices.end(); it++)
        blocksRecursive[b].push_back(recurs_vars[*it]);
    }
}

void
StaticModel::computeBlockMFSJacobian()
{
  blocksMFSJacobian.clear();
  for(int b = 0; b < (int) blocks.size(); b++)
    {
      // Create the map of recursive vars to their normalized equation
      map<int, NodeID> recursive2eq;
      for(vector<int>::const_iterator it = blocksRecursive[b].begin();
          it != blocksRecursive[b].end(); it++)
        recursive2eq[symbol_table.getID(eEndogenous, *it)] = equations[endo2eq[*it]];

      for(set<int>::const_iterator it = blocksMFS[b].begin();
          it != blocksMFS[b].end(); it++)
        {
          int eq_no = endo2eq[*it];
          for(set<int>::const_iterator it2 = blocksMFS[b].begin();
              it2 != blocksMFS[b].end(); it2++)
            {
              int deriv_id = symbol_table.getID(eEndogenous, *it2);
              NodeID d = equations[eq_no]->getChainRuleDerivative(deriv_id, recursive2eq);
              if (d != Zero)
                blocksMFSJacobian[make_pair(eq_no, *it2)] = d;
            }
        }
    }
}

void
StaticModel::writeOutput(ostream &output) const
{
  if (!block_mfs)
    return;

  output << "M_.blocksMFS = cell(" << blocksMFS.size() << ", 1);" << endl;
  for(int b = 0; b < (int) blocks.size(); b++)
    {
      output << "M_.blocksMFS{" << b+1 << "} = [ ";
      transform(blocksMFS[b].begin(), blocksMFS[b].end(), ostream_iterator<int>(output, "; "), bind2nd(plus<int>(), 1));
      output << "];" << endl;
    }
}

void
StaticModel::writeStaticBlockMFSFile(ostream &output, const string &func_name) const
{
  output << "function [residual, g1, y] = " << func_name << "(nblock, y, x, params)" << endl
         << "  switch nblock" << endl;

  for(int b = 0; b < (int) blocks.size(); b++)
    {
      output << "    case " << b+1 << endl
             << "      % Variables not in minimum feedback set" << endl;
      for(vector<int>::const_iterator it = blocksRecursive[b].begin();
          it != blocksRecursive[b].end(); it++)
        {
          equations[endo2eq[*it]]->writeOutput(output, oMatlabStaticModel, temporary_terms_type());
          output << ";" << endl;
        }

      output << "      % Model residuals" << endl
             << "residual = zeros(" << blocksMFS[b].size() << ",1);" << endl;

      int i = 1;
      for (set<int>::const_iterator it = blocksMFS[b].begin();
           it != blocksMFS[b].end(); it++)
        {
          output << "residual(" << i << ")=(";

          BinaryOpNode *eq_node = equations[endo2eq[*it]];

          NodeID lhs = eq_node->get_arg1();
          lhs->writeOutput(output, oMatlabStaticModel, temporary_terms_type());
          output << ")-(";

          NodeID rhs = eq_node->get_arg2();
          rhs->writeOutput(output, oMatlabStaticModel, temporary_terms_type());
          output << ");" << endl;

          i++;
        }

      output << "      % Jacobian matrix" << endl
             << "g1 = zeros(" << blocksMFS[b].size() << ", " << blocksMFS[b].size() << ");" << endl;
      i = 1;
      for (set<int>::const_iterator it = blocksMFS[b].begin();
           it != blocksMFS[b].end(); it++)
        {
          int eq_no = endo2eq[*it];
          int j = 1;
          for (set<int>::const_iterator it2 = blocksMFS[b].begin();
               it2 != blocksMFS[b].end(); it2++)
            {
              map<pair<int, int>, NodeID>::const_iterator itd = blocksMFSJacobian.find(make_pair(eq_no, *it2));
              if (itd != blocksMFSJacobian.end())
                {
                  output << "g1(" << i << "," << j << ")=";
                  itd->second->writeOutput(output, oMatlabStaticModel, temporary_terms_type());
                  output << ";" << endl;
                }
              j++;
            }
          i++;
        }
    }

  output << "  end" << endl
         << "end" << endl;
}
