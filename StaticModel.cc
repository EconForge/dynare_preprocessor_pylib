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
#include <functional>

#ifdef DEBUG
# include <ext/functional>
using namespace __gnu_cxx;
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
  ModelTree(symbol_table_arg, num_constants_arg)
{
}

void
StaticModel::writeStaticMFile(const string &static_basename) const
{
  string filename = static_basename + ".m";

  ofstream mStaticModelFile;
  mStaticModelFile.open(filename.c_str(), ios::out | ios::binary);
  if (!mStaticModelFile.is_open())
    {
      cerr << "Error: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }
  // Writing comments and function definition command
  mStaticModelFile << "function [residual, g1, g2] = " << static_basename << "(y, x, params)" << endl
                   << "%" << endl
                   << "% Status : Computes static model for Dynare" << endl
                   << "%" << endl
                   << "% Warning : this file is generated automatically by Dynare" << endl
                   << "%           from model file (.mod)" << endl << endl;

  writeStaticModel(mStaticModelFile);

  mStaticModelFile.close();
}

void
StaticModel::writeStaticCFile(const string &static_basename) const
{
  string filename = static_basename + ".c";

  ofstream mStaticModelFile;
  mStaticModelFile.open(filename.c_str(), ios::out | ios::binary);
  if (!mStaticModelFile.is_open())
    {
      cerr << "Error: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }
  mStaticModelFile << "/*" << endl
                   << " * " << filename << " : Computes static model for Dynare" << endl
                   << " * Warning : this file is generated automatically by Dynare" << endl
                   << " *           from model file (.mod)" << endl
                   << endl
                   << " */" << endl
                   << "#include <math.h>" << endl
                   << "#include \"mex.h\"" << endl;

  // Writing the function Static
  writeStaticModel(mStaticModelFile);

  // Writing the gateway routine
  mStaticModelFile << "/* The gateway routine */" << endl
                   << "void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])" << endl
                   << "{" << endl
                   << "  double *y, *x, *params;" << endl
                   << "  double *residual, *g1;" << endl
                   << endl
                   << "  /* Create a pointer to the input matrix y. */" << endl
                   << "  y = mxGetPr(prhs[0]);" << endl
                   << endl
                   << "  /* Create a pointer to the input matrix x. */" << endl
                   << "  x = mxGetPr(prhs[1]);" << endl
                   << endl
                   << "  /* Create a pointer to the input matrix params. */" << endl
                   << "  params = mxGetPr(prhs[2]);" << endl
                   << endl
                   << "  residual = NULL;" << endl
                   << "  if (nlhs >= 1)" << endl
                   << "  {" << endl
                   << "      /* Set the output pointer to the output matrix residual. */" << endl
                   << "      plhs[0] = mxCreateDoubleMatrix(" << equations.size() << ",1, mxREAL);" << endl
                   << "     /* Create a C pointer to a copy of the output matrix residual. */" << endl
                   << "     residual = mxGetPr(plhs[0]);" << endl
                   << "  }" << endl
                   << endl
                   << "  g1 = NULL;" << endl
                   << "  if (nlhs >= 2)" << endl
                   << "  {" << endl
                   << "      /* Set the output pointer to the output matrix g1. */" << endl
                   << "      plhs[1] = mxCreateDoubleMatrix(" << equations.size() << ", " << symbol_table.endo_nbr() << ", mxREAL);" << endl
                   << "      /* Create a C pointer to a copy of the output matrix g1. */" << endl
                   << "      g1 = mxGetPr(plhs[1]);" << endl
                   << "  }" << endl
                   << endl
                   << "  /* Call the C Static. */" << endl
                   << "  Static(y, x, params, residual, g1);" << endl
                   << "}" << endl;

  mStaticModelFile.close();
}

void
StaticModel::writeStaticModel(ostream &StaticOutput) const
{
  ostringstream model_output;    // Used for storing model equations
  ostringstream jacobian_output; // Used for storing jacobian equations
  ostringstream hessian_output;

  ExprNodeOutputType output_type = (mode == eDLLMode ? oCStaticModel : oMatlabStaticModel);

  writeModelLocalVariables(model_output, output_type);

  writeTemporaryTerms(temporary_terms, model_output, output_type);

  writeModelEquations(model_output, output_type);

  // Write Jacobian w.r. to endogenous only
  for (first_derivatives_type::const_iterator it = first_derivatives.begin();
       it != first_derivatives.end(); it++)
    {
      int eq = it->first.first;
      int symb_id = it->first.second;
      NodeID d1 = it->second;

      ostringstream g1;
      g1 << "  g1";
      jacobianHelper(g1, eq, symbol_table.getTypeSpecificID(symb_id), output_type);

      jacobian_output << g1.str() << "=" << g1.str() << "+";
      d1->writeOutput(jacobian_output, output_type, temporary_terms);
      jacobian_output << ";" << endl;
    }

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

      hessian_output << "v2";
      hessianHelper(hessian_output, k, 0, output_type);
      hessian_output << "=" << eq + 1 << ";" << endl;

      hessian_output << "v2";
      hessianHelper(hessian_output, k, 1, output_type);
      hessian_output << "=" << col_nb + 1 << ";" << endl;

      hessian_output << "v2";
      hessianHelper(hessian_output, k, 2, output_type);
      hessian_output << "=";
      d2->writeOutput(hessian_output, output_type, temporary_terms);
      hessian_output << ";" << endl;

      k++;

      // Treating symetric elements
      if (symb_id1 != symb_id2)
        {
          hessian_output << "v2";
          hessianHelper(hessian_output, k, 0, output_type);
          hessian_output << "=" << eq + 1 << ";" << endl;

          hessian_output << "v2";
          hessianHelper(hessian_output, k, 1, output_type);
          hessian_output << "=" << col_nb_sym + 1 << ";" << endl;

          hessian_output << "v2";
          hessianHelper(hessian_output, k, 2, output_type);
          hessian_output << "=v2"; 
          hessianHelper(hessian_output, k-1, 2, output_type);
          hessian_output << ";" << endl;

          k++;
        }
    }

  // Writing ouputs
  if (mode != eDLLMode)
    {
      StaticOutput << "residual = zeros( " << equations.size() << ", 1);" << endl << endl
                   << "%" << endl
                   << "% Model equations" << endl
                   << "%" << endl
                   << endl
                   << model_output.str()
                   << "if ~isreal(residual)" << endl
                   << "  residual = real(residual)+imag(residual).^2;" << endl
                   << "end" << endl
                   << "if nargout >= 2," << endl
                   << "  g1 = zeros(" << equations.size() << ", " << symbol_table.endo_nbr() << ");" << endl
                   << endl
                   << "%" << endl
                   << "% Jacobian matrix" << endl
                   << "%" << endl
                   << endl
                   << jacobian_output.str()
                   << "  if ~isreal(g1)" << endl
                   << "    g1 = real(g1)+2*imag(g1);" << endl
                   << "  end" << endl
                   << "end" << endl;

      // Initialize g2 matrix
      StaticOutput << "if nargout >= 3," << endl
                   << "%" << endl
                   << "% Hessian matrix" << endl
                   << "%" << endl
                   << endl;
      int ncols = symbol_table.endo_nbr() * symbol_table.endo_nbr();
      if (second_derivatives.size())
        StaticOutput << "  v2 = zeros(" << NNZDerivatives[1] << ",3);" << endl
                     << hessian_output.str()
                     << "  g2 = sparse(v2(:,1),v2(:,2),v2(:,3)," << equations.size() << "," << ncols << ");" << endl;
      else // Either hessian is all zero, or we didn't compute it
        StaticOutput << "  g2 = sparse([],[],[]," << equations.size() << "," << ncols << ");" << endl;
      StaticOutput << "end;" << endl;
    }
  else
    {
      StaticOutput << "void Static(double *y, double *x, double *params, double *residual, double *g1)" << endl
                   << "{" << endl
                   << "  double lhs, rhs;" << endl
        // Writing residual equations
                   << "  /* Residual equations */" << endl
                   << "  if (residual == NULL)" << endl
                   << "    return;" << endl
                   << "  else" << endl
                   << "    {" << endl
                   << model_output.str()
        // Writing Jacobian
                   << "     /* Jacobian for endogenous variables without lag */" << endl
                   << "     if (g1 == NULL)" << endl
                   << "       return;" << endl
                   << "     else" << endl
                   << "       {" << endl
                   << jacobian_output.str()
                   << "       }" << endl
                   << "    }" << endl
                   << "}" << endl << endl;
    }
}

void
StaticModel::writeStaticFile(const string &basename) const
{
  switch (mode)
    {
    case eStandardMode:
    case eSparseDLLMode:
    case eSparseMode:
      writeStaticMFile(basename + "_static");
      break;
    case eDLLMode:
      writeStaticCFile(basename + "_static");
      break;
    }
}

void
StaticModel::computingPass(bool block_mfs, bool hessian, bool no_tmp_terms)
{
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

  if (!no_tmp_terms)
    computeTemporaryTerms();

  if (block_mfs)
    {
      vector<int> endo2eq(equation_number());
      computeNormalization(endo2eq);

      vector<set<int> > blocks;
      computeSortedBlockDecomposition(blocks, endo2eq);

      vector<set<int> > blocksMFS;
      computeMFS(blocksMFS, blocks, endo2eq);
    }
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
StaticModel::computeNormalization(vector<int> &endo2eq) const
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
        add_edge(i + n, symbol_table.getTypeSpecificID(it->first), g);
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

  assert((int) endo2eq.size() == n);

  // Create the resulting map, by copying the n first elements of mate_map, and substracting n to them
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
      if (endo.find(make_pair(symb_id, 0)) != endo.end())
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
StaticModel::computeSortedBlockDecomposition(vector<set<int> > &blocks, const vector<int> &endo2eq) const
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
        add_edge(symbol_table.getTypeSpecificID(it->first), i, g);
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
StaticModel::computeMFS(vector<set<int> > &blocksMFS, const vector<set<int> > &blocks, const vector<int> &endo2eq) const
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
            {
              const int tsid = symbol_table.getTypeSpecificID(it2->first);
              if (blocks[b].find(tsid) != blocks[b].end()) // Add edge only if vertex member of this block
                add_edge(tsid2vertex[tsid], tsid2vertex[*it], g);
            }
        }

      // Compute minimum feedback set
      MFS::Minimal_set_of_feedback_vertex(blocksMFS[b], g);

      cout << "Block " << b << ": " << blocksMFS[b].size() << "/" << blocks[b].size() << " in MFS" << endl;
    }
}
