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
#include <iostream>
#include <fstream>

#include "ModelTree.hh"
#include "MinimumFeedbackSet.hh"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/max_cardinality_matching.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/topological_sort.hpp>

using namespace boost;
using namespace MFS;

bool
ModelTree::computeNormalization(const jacob_map &contemporaneous_jacobian, bool verbose)
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

  for (jacob_map::const_iterator it = contemporaneous_jacobian.begin(); it != contemporaneous_jacobian.end(); it++)
    add_edge(it->first.first + n, it->first.second, g);

  // Compute maximum cardinality matching
  vector<int> mate_map(2*n);

#if 1
  bool check = checked_edmonds_maximum_cardinality_matching(g, &mate_map[0]);
#else // Alternative way to compute normalization, by giving an initial matching using natural normalizations
  fill(mate_map.begin(), mate_map.end(), graph_traits<BipartiteGraph>::null_vertex());

  multimap<int, int> natural_endo2eqs;
  computeNormalizedEquations(natural_endo2eqs);

  for (int i = 0; i < symbol_table.endo_nbr(); i++)
    {
      if (natural_endo2eqs.count(i) == 0)
        continue;

      int j = natural_endo2eqs.find(i)->second;

      put(&mate_map[0], i, n+j);
      put(&mate_map[0], n+j, i);
    }

  edmonds_augmenting_path_finder<BipartiteGraph, size_t *, property_map<BipartiteGraph, vertex_index_t>::type> augmentor(g, &mate_map[0], get(vertex_index, g));
  bool not_maximum_yet = true;
  while (not_maximum_yet)
    {
      not_maximum_yet = augmentor.augment_matching();
    }
  augmentor.get_current_matching(&mate_map[0]);

  bool check = maximum_cardinality_matching_verifier<BipartiteGraph, size_t *, property_map<BipartiteGraph, vertex_index_t>::type>::verify_matching(g, &mate_map[0], get(vertex_index, g));
#endif

  assert(check);

#ifdef DEBUG
  for (int i = 0; i < n; i++)
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

  for (int i = 0; i < symbol_table.endo_nbr(); i++)
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

  // Check if all variables are normalized
  vector<int>::const_iterator it = find(mate_map.begin(), mate_map.begin() + n, graph_traits<BipartiteGraph>::null_vertex());
  if (it != mate_map.begin() + n)
    {
      if (verbose)
        cerr << "ERROR: Could not normalize the model. Variable "
             << symbol_table.getName(symbol_table.getID(eEndogenous, it - mate_map.begin()))
             << " is not in the maximum cardinality matching." << endl;
      check = false;
    }
  return check;
}

void
ModelTree::computeNonSingularNormalization(jacob_map &contemporaneous_jacobian, double cutoff, jacob_map &static_jacobian, dynamic_jacob_map &dynamic_jacobian)
{
  bool check = false;

  cout << "Normalizing the model..." << endl;

  int n = equation_number();

  // compute the maximum value of each row of the contemporaneous Jacobian matrix
  //jacob_map normalized_contemporaneous_jacobian;
  jacob_map normalized_contemporaneous_jacobian(contemporaneous_jacobian);
  vector<double> max_val(n, 0.0);
  for (jacob_map::const_iterator iter = contemporaneous_jacobian.begin(); iter != contemporaneous_jacobian.end(); iter++)
    if (fabs(iter->second) > max_val[iter->first.first])
      max_val[iter->first.first] = fabs(iter->second);

  for (jacob_map::iterator iter = normalized_contemporaneous_jacobian.begin(); iter != normalized_contemporaneous_jacobian.end(); iter++)
    iter->second /= max_val[iter->first.first];

  //We start with the highest value of the cutoff and try to normalize the model
  double current_cutoff = 0.99999999;

  int suppressed = 0;
  while (!check && current_cutoff > 1e-19)
    {
      jacob_map tmp_normalized_contemporaneous_jacobian;
      int suppress = 0;
      for (jacob_map::iterator iter = normalized_contemporaneous_jacobian.begin(); iter != normalized_contemporaneous_jacobian.end(); iter++)
        if (fabs(iter->second) > max(current_cutoff, cutoff))
          tmp_normalized_contemporaneous_jacobian[make_pair(iter->first.first, iter->first.second)] = iter->second;
        else
          suppress++;

      if (suppress != suppressed)
        check = computeNormalization(tmp_normalized_contemporaneous_jacobian, false);
      suppressed = suppress;
      if (!check)
        {
          current_cutoff /= 2;
          // In this last case try to normalize with the complete jacobian
          if (current_cutoff <= 1e-19)
            check = computeNormalization(normalized_contemporaneous_jacobian, false);
        }
    }

  if (!check)
    {
      cout << "Normalization failed with cutoff, trying symbolic normalization..." << endl;
      //if no non-singular normalization can be found, try to find a normalization even with a potential singularity
      jacob_map tmp_normalized_contemporaneous_jacobian;
      set<pair<int, int> > endo;
      for (int i = 0; i < n; i++)
        {
          endo.clear();
          equations[i]->collectEndogenous(endo);
          for (set<pair<int, int> >::const_iterator it = endo.begin(); it != endo.end(); it++)
            tmp_normalized_contemporaneous_jacobian[make_pair(i, it->first)] = 1;
        }
      check = computeNormalization(tmp_normalized_contemporaneous_jacobian, true);
      if (check)
        {
          // Update the jacobian matrix
          for (jacob_map::const_iterator it = tmp_normalized_contemporaneous_jacobian.begin(); it != tmp_normalized_contemporaneous_jacobian.end(); it++)
            {
              if (static_jacobian.find(make_pair(it->first.first, it->first.second)) == static_jacobian.end())
                static_jacobian[make_pair(it->first.first, it->first.second)] = 0;
              if (dynamic_jacobian.find(make_pair(0, make_pair(it->first.first, it->first.second))) == dynamic_jacobian.end())
                dynamic_jacobian[make_pair(0, make_pair(it->first.first, it->first.second))] = 0;
              if (contemporaneous_jacobian.find(make_pair(it->first.first, it->first.second)) == contemporaneous_jacobian.end())
                contemporaneous_jacobian[make_pair(it->first.first, it->first.second)] = 0;
            }
        }
    }

  if (!check)
    {
      cerr << "No normalization could be computed. Aborting." << endl;
      exit(EXIT_FAILURE);
    }
}

void
ModelTree::computeNormalizedEquations(multimap<int, int> &endo2eqs) const
{
  for (int i = 0; i < equation_number(); i++)
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
ModelTree::evaluateAndReduceJacobian(const eval_context_type &eval_context, jacob_map &contemporaneous_jacobian, jacob_map &static_jacobian, dynamic_jacob_map &dynamic_jacobian, double cutoff, bool verbose)
{
  int nb_elements_contemparenous_Jacobian = 0;
  set<pair<int, int> > jacobian_elements_to_delete;
  for (first_derivatives_type::iterator it = first_derivatives.begin();
       it != first_derivatives.end(); it++)
    {
      int deriv_id = it->first.second;
      if (getTypeByDerivID(deriv_id) == eEndogenous)
        {
          NodeID Id = it->second;
          int eq = it->first.first;
          int symb = getSymbIDByDerivID(deriv_id);
          int var = symbol_table.getTypeSpecificID(symb);
          int lag = getLagByDerivID(deriv_id);
          double val = 0;
          try
            {
              val = Id->eval(eval_context);
            }
          catch (ExprNode::EvalException &e)
            {
              cerr << "ERROR: evaluation of Jacobian failed for equation " << eq+1 << " and variable " << symbol_table.getName(symb) << "(" << lag << ") [" << symb << "] !" << endl;
              Id->writeOutput(cerr, oMatlabDynamicModelSparse, temporary_terms);
              cerr << endl;
              exit(EXIT_FAILURE);
            }

          if (fabs(val) < cutoff)
            {
              if (verbose)
                cout << "the coefficient related to variable " << var << " with lag " << lag << " in equation " << eq << " is equal to " << val << " and is set to 0 in the incidence matrix (size=" << symbol_table.endo_nbr() << ")" << endl;
              jacobian_elements_to_delete.insert(make_pair(eq, deriv_id));
            }
          else
            {
              if (lag == 0)
                {
                  nb_elements_contemparenous_Jacobian++;
                  contemporaneous_jacobian[make_pair(eq, var)] = val;
                }
              if (static_jacobian.find(make_pair(eq, var)) != static_jacobian.end())
                static_jacobian[make_pair(eq, var)] += val;
              else
                static_jacobian[make_pair(eq, var)] = val;
              dynamic_jacobian[make_pair(lag, make_pair(eq, var))] = Id;
            }
        }
    }

  // Get rid of the elements of the Jacobian matrix below the cutoff
  for (set<pair<int, int> >::const_iterator it = jacobian_elements_to_delete.begin(); it != jacobian_elements_to_delete.end(); it++)
    first_derivatives.erase(*it);

  if (jacobian_elements_to_delete.size() > 0)
    {
      cout << jacobian_elements_to_delete.size() << " elements among " << first_derivatives.size() << " in the incidence matrices are below the cutoff (" << cutoff << ") and are discarded" << endl
           << "The contemporaneous incidence matrix has " << nb_elements_contemparenous_Jacobian << " elements" << endl;
    }
}

void
ModelTree::computePrologueAndEpilogue(jacob_map &static_jacobian_arg, vector<int> &equation_reordered, vector<int> &variable_reordered, unsigned int &prologue, unsigned int &epilogue)
{
  vector<int> eq2endo(equation_number(), 0);
  equation_reordered.resize(equation_number());
  variable_reordered.resize(equation_number());
  bool *IM;
  int n = equation_number();
  IM = (bool *) calloc(n*n, sizeof(bool));
  int i = 0;
  for (vector<int>::const_iterator it = endo2eq.begin(); it != endo2eq.end(); it++, i++)
    {
      eq2endo[*it] = i;
      equation_reordered[i] = i;
      variable_reordered[*it] = i;
    }
  for (jacob_map::const_iterator it = static_jacobian_arg.begin(); it != static_jacobian_arg.end(); it++)
    IM[it->first.first * n + endo2eq[it->first.second]] = true;
  bool something_has_been_done = true;
  prologue = 0;
  int k = 0;
  // Find the prologue equations and place first the AR(1) shock equations first
  while (something_has_been_done)
    {
      int tmp_prologue = prologue;
      something_has_been_done = false;
      for (int i = prologue; i < n; i++)
        {
          int nze = 0;
          for (int j = tmp_prologue; j < n; j++)
            if (IM[i * n + j])
              {
                nze++;
                k = j;
              }
          if (nze == 1)
            {
              for (int j = 0; j < n; j++)
                {
                  bool tmp_bool = IM[tmp_prologue * n + j];
                  IM[tmp_prologue * n + j] = IM[i * n + j];
                  IM[i * n + j] = tmp_bool;
                }
              int tmp = equation_reordered[tmp_prologue];
              equation_reordered[tmp_prologue] = equation_reordered[i];
              equation_reordered[i] = tmp;
              for (int j = 0; j < n; j++)
                {
                  bool tmp_bool = IM[j * n + tmp_prologue];
                  IM[j * n + tmp_prologue] = IM[j * n + k];
                  IM[j * n + k] = tmp_bool;
                }
              tmp = variable_reordered[tmp_prologue];
              variable_reordered[tmp_prologue] = variable_reordered[k];
              variable_reordered[k] = tmp;
              tmp_prologue++;
              something_has_been_done = true;
            }
        }
      prologue = tmp_prologue;
    }

  something_has_been_done = true;
  epilogue = 0;
  // Find the epilogue equations
  while (something_has_been_done)
    {
      int tmp_epilogue = epilogue;
      something_has_been_done = false;
      for (int i = prologue; i < n - (int) epilogue; i++)
        {
          int nze = 0;
          for (int j = prologue; j < n - tmp_epilogue; j++)
            if (IM[j * n + i])
              {
                nze++;
                k = j;
              }
          if (nze == 1)
            {
              for (int j = 0; j < n; j++)
                {
                  bool tmp_bool = IM[(n - 1 - tmp_epilogue) * n + j];
                  IM[(n - 1 - tmp_epilogue) * n + j] = IM[k * n + j];
                  IM[k * n + j] = tmp_bool;
                }
              int tmp = equation_reordered[n - 1 - tmp_epilogue];
              equation_reordered[n - 1 - tmp_epilogue] = equation_reordered[k];
              equation_reordered[k] = tmp;
              for (int j = 0; j < n; j++)
                {
                  bool tmp_bool = IM[j * n + n - 1 - tmp_epilogue];
                  IM[j * n + n - 1 - tmp_epilogue] = IM[j * n + i];
                  IM[j * n + i] = tmp_bool;
                }
              tmp = variable_reordered[n - 1 - tmp_epilogue];
              variable_reordered[n - 1 - tmp_epilogue] = variable_reordered[i];
              variable_reordered[i] = tmp;
              tmp_epilogue++;
              something_has_been_done = true;
            }
        }
      epilogue = tmp_epilogue;
    }
  free(IM);
}

t_equation_type_and_normalized_equation
ModelTree::equationTypeDetermination(vector<BinaryOpNode *> &equations, map<pair<int, pair<int, int> >, NodeID> &first_order_endo_derivatives, vector<int> &Index_Var_IM, vector<int> &Index_Equ_IM, int mfs)
{
  NodeID lhs, rhs;
  ostringstream tmp_output;
  BinaryOpNode *eq_node;
  ostringstream tmp_s;
  temporary_terms_type temporary_terms;
  EquationType Equation_Simulation_Type;
  t_equation_type_and_normalized_equation V_Equation_Simulation_Type(equations.size());
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
      if (derivative != first_order_endo_derivatives.end())
        {
          set<pair<int, int> > result;
          derivative->second->collectEndogenous(result);
          set<pair<int, int> >::const_iterator d_endo_variable = result.find(make_pair(var, 0));
          //Determine whether the equation could be evaluated rather than to be solved
          ostringstream tt("");
          derivative->second->writeOutput(tt, oMatlabDynamicModelSparse, temporary_terms);
          if (tmp_output.str() == tmp_s.str() and tt.str() == "1")
            {
              Equation_Simulation_Type = E_EVALUATE;
            }
          else
            {
              vector<pair<int, pair<NodeID, NodeID> > > List_of_Op_RHS;
              res =  equations[eq]->normalizeEquation(var, List_of_Op_RHS);
              if (mfs == 2)
                {
                  if (d_endo_variable == result.end() && res.second)
                    Equation_Simulation_Type = E_EVALUATE_S;
                }
              else if (mfs == 3)
                {
                  if (res.second) // The equation could be solved analytically
                    Equation_Simulation_Type = E_EVALUATE_S;
                }
            }
        }
      V_Equation_Simulation_Type[eq] = make_pair(Equation_Simulation_Type, dynamic_cast<BinaryOpNode *>(res.second));
    }
  return (V_Equation_Simulation_Type);
}

void
ModelTree::getVariableLeadLagByBlock(dynamic_jacob_map &dynamic_jacobian, vector<int > &components_set, int nb_blck_sim, int prologue, int epilogue, t_lag_lead_vector &equation_lead_lag, t_lag_lead_vector &variable_lead_lag, vector<int> equation_reordered, vector<int> variable_reordered) const
{
  int nb_endo = symbol_table.endo_nbr();
  variable_lead_lag = t_lag_lead_vector(nb_endo, make_pair(0, 0));
  equation_lead_lag = t_lag_lead_vector(nb_endo, make_pair(0, 0));
  vector<int> variable_blck(nb_endo), equation_blck(nb_endo);
  for (int i = 0; i < nb_endo; i++)
    {
      if (i < prologue)
        {
          variable_blck[variable_reordered[i]] = i;
          equation_blck[equation_reordered[i]] = i;
        }
      else if (i < (int) components_set.size() + prologue)
        {
          variable_blck[variable_reordered[i]] = components_set[i-prologue] + prologue;
          equation_blck[equation_reordered[i]] = components_set[i-prologue] + prologue;
        }
      else
        {
          variable_blck[variable_reordered[i]] = i- (nb_endo - nb_blck_sim - prologue - epilogue);
          equation_blck[equation_reordered[i]] = i- (nb_endo - nb_blck_sim - prologue - epilogue);
        }
    }
  for (dynamic_jacob_map::const_iterator it = dynamic_jacobian.begin(); it != dynamic_jacobian.end(); it++)
    {
      int lag = it->first.first;
      int j_1 = it->first.second.second;
      int i_1 = it->first.second.second;
      if (variable_blck[i_1] == equation_blck[j_1])
        {
          if (lag > variable_lead_lag[i_1].second)
            variable_lead_lag[i_1] = make_pair(variable_lead_lag[i_1].first, lag);
          if (lag < -variable_lead_lag[i_1].first)
            variable_lead_lag[i_1] = make_pair(-lag, variable_lead_lag[i_1].second);
          if (lag > equation_lead_lag[j_1].second)
            equation_lead_lag[j_1] = make_pair(equation_lead_lag[j_1].first, lag);
          if (lag < -equation_lead_lag[j_1].first)
            equation_lead_lag[j_1] = make_pair(-lag, equation_lead_lag[j_1].second);
        }
    }
}

void
ModelTree::computeBlockDecompositionAndFeedbackVariablesForEachBlock(jacob_map &static_jacobian, dynamic_jacob_map &dynamic_jacobian, int prologue, int epilogue, vector<int> &equation_reordered, vector<int> &variable_reordered, vector<pair<int, int> > &blocks, t_equation_type_and_normalized_equation &Equation_Type, bool verbose_, bool select_feedback_variable, int mfs, vector<int> &inv_equation_reordered, vector<int> &inv_variable_reordered) const
{
  int nb_var = variable_reordered.size();
  int n = nb_var - prologue - epilogue;
  typedef adjacency_list<vecS, vecS, directedS> DirectedGraph;

  GraphvizDigraph G2(n);

  vector<int> reverse_equation_reordered(nb_var), reverse_variable_reordered(nb_var);

  for (int i = 0; i < nb_var; i++)
    {
      reverse_equation_reordered[equation_reordered[i]] = i;
      reverse_variable_reordered[variable_reordered[i]] = i;
    }

  for (jacob_map::const_iterator it = static_jacobian.begin(); it != static_jacobian.end(); it++)
    if (reverse_equation_reordered[it->first.first] >= prologue && reverse_equation_reordered[it->first.first] < nb_var - epilogue
        && reverse_variable_reordered[it->first.second] >= prologue && reverse_variable_reordered[it->first.second] < nb_var - epilogue
        && it->first.first != endo2eq[it->first.second])
      add_edge(reverse_equation_reordered[it->first.first]-prologue, reverse_equation_reordered[endo2eq[it->first.second]]-prologue, G2);

  vector<int> endo2block(num_vertices(G2)), discover_time(num_vertices(G2));

  // Compute strongly connected components
  int num = strong_components(G2, &endo2block[0]);

  blocks = vector<pair<int, int> >(num, make_pair(0, 0));

  // Create directed acyclic graph associated to the strongly connected components
  typedef adjacency_list<vecS, vecS, directedS> DirectedGraph;
  DirectedGraph dag(num);

  for (unsigned int i = 0; i < num_vertices(G2); i++)
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
  for (int i = 0; i < num; i++)
    unordered2ordered[ordered2unordered[i]] = i;

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

  t_lag_lead_vector equation_lag_lead, variable_lag_lead;

  getVariableLeadLagByBlock(dynamic_jacobian, endo2block, num, prologue, epilogue, equation_lag_lead, variable_lag_lead, equation_reordered, variable_reordered);

  vector<int> tmp_equation_reordered(equation_reordered), tmp_variable_reordered(variable_reordered);
  int order = prologue;
  //Add a loop on vertices which could not be normalized or vertices related to lead variables => force those vertices to belong to the feedback set
  if (select_feedback_variable)
    {
      for (int i = 0; i < n; i++)
        if (Equation_Type[equation_reordered[i+prologue]].first == E_SOLVE
            or variable_lag_lead[variable_reordered[i+prologue]].second > 0 or variable_lag_lead[variable_reordered[i+prologue]].first > 0
            or equation_lag_lead[equation_reordered[i+prologue]].second > 0 or equation_lag_lead[equation_reordered[i+prologue]].first > 0
            or mfs == 0)
          add_edge(i, i, G2);
    }
  else
    {
      for (int i = 0; i < n; i++)
        if (Equation_Type[equation_reordered[i+prologue]].first == E_SOLVE || mfs == 0)
          add_edge(i, i, G2);
    }
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
          equation_reordered[order] = tmp_equation_reordered[*its+prologue];
          variable_reordered[order] = tmp_variable_reordered[*its+prologue];
          order++;
        }
      components_set[i].second.second = Reordered_Vertice;
      //Second we have the equations related to the feedback variables
      for (set<int>::iterator its = feed_back_vertices.begin(); its != feed_back_vertices.end(); its++)
        {
          equation_reordered[order] = tmp_equation_reordered[v_index[vertex(*its, G)]+prologue];
          variable_reordered[order] = tmp_variable_reordered[v_index[vertex(*its, G)]+prologue];
          order++;
        }
    }
  inv_equation_reordered = vector<int>(nb_var);
  inv_variable_reordered = vector<int>(nb_var);
  for (int i = 0; i < nb_var; i++)
    {
      inv_variable_reordered[variable_reordered[i]] = i;
      inv_equation_reordered[equation_reordered[i]] = i;
    }
}

void
ModelTree::printBlockDecomposition(vector<pair<int, int> > blocks)
{
  int largest_block = 0;
  int Nb_SimulBlocks = 0;
  int Nb_feedback_variable = 0;
  unsigned int Nb_TotalBlocks = getNbBlocks();
  for (unsigned int block = 0; block < Nb_TotalBlocks; block++)
    {
      BlockSimulationType simulation_type = getBlockSimulationType(block);
      if (simulation_type == SOLVE_FORWARD_COMPLETE || simulation_type == SOLVE_BACKWARD_COMPLETE || simulation_type == SOLVE_TWO_BOUNDARIES_COMPLETE)
        {
          Nb_SimulBlocks++;
          int size = getBlockSize(block);
          if (size > largest_block)
            {
              largest_block = size;
              Nb_feedback_variable = blocks[Nb_SimulBlocks-1].second;
            }
        }
    }

  int Nb_RecursBlocks = Nb_TotalBlocks - Nb_SimulBlocks;
  cout << Nb_TotalBlocks << " block(s) found:" << endl
       << "  " << Nb_RecursBlocks << " recursive block(s) and " << Nb_SimulBlocks << " simultaneous block(s)." << endl
       << "  the largest simultaneous block has " << largest_block << " equation(s)" << endl
       << "                                 and " << Nb_feedback_variable << " feedback variable(s)." << endl;
}

t_block_type_firstequation_size_mfs
ModelTree::reduceBlocksAndTypeDetermination(dynamic_jacob_map &dynamic_jacobian, int prologue, int epilogue, vector<pair<int, int> > &blocks, vector<BinaryOpNode *> &equations, t_equation_type_and_normalized_equation &Equation_Type, vector<int> &variable_reordered, vector<int> &equation_reordered)
{
  int i = 0;
  int count_equ = 0, blck_count_simult = 0;
  int Blck_Size, MFS_Size;
  int Lead, Lag;
  t_block_type_firstequation_size_mfs block_type_size_mfs;
  BlockSimulationType Simulation_Type, prev_Type = UNKNOWN;
  int eq = 0;
  for (i = 0; i < prologue+(int) blocks.size()+epilogue; i++)
    {
      int first_count_equ = count_equ;
      if (i < prologue)
        {
          Blck_Size = 1;
          MFS_Size = 1;
        }
      else if (i < prologue+(int) blocks.size())
        {
          Blck_Size = blocks[blck_count_simult].first;
          MFS_Size = blocks[blck_count_simult].second;
          blck_count_simult++;
        }
      else if (i < prologue+(int) blocks.size()+epilogue)
        {
          Blck_Size = 1;
          MFS_Size = 1;
        }

      Lag = Lead = 0;
      set<pair<int, int> > endo;
      for (count_equ  = first_count_equ; count_equ  < Blck_Size+first_count_equ; count_equ++)
        {
          endo.clear();
          equations[equation_reordered[count_equ]]->collectEndogenous(endo);
          for (set<pair<int, int> >::const_iterator it = endo.begin(); it != endo.end(); it++)
            {
              int curr_variable = it->first;
              int curr_lag = it->second;
              vector<int>::const_iterator it = find(variable_reordered.begin()+first_count_equ, variable_reordered.begin()+(first_count_equ+Blck_Size), curr_variable);
              if (it != variable_reordered.begin()+(first_count_equ+Blck_Size))
                if (dynamic_jacobian.find(make_pair(curr_lag, make_pair(equation_reordered[count_equ], curr_variable))) != dynamic_jacobian.end())
                  {
                    if (curr_lag > Lead)
                      Lead = curr_lag;
                    else if (-curr_lag > Lag)
                      Lag = -curr_lag;
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
          if (Equation_Type[equation_reordered[eq]].first == E_EVALUATE or Equation_Type[equation_reordered[eq]].first == E_EVALUATE_S)
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
                  //merge the current block with the previous one
                  BlockSimulationType c_Type = (block_type_size_mfs[block_type_size_mfs.size()-1]).first.first;
                  int c_Size = (block_type_size_mfs[block_type_size_mfs.size()-1]).second.first;
                  int first_equation = (block_type_size_mfs[block_type_size_mfs.size()-1]).first.second;
                  block_type_size_mfs[block_type_size_mfs.size()-1] = make_pair(make_pair(c_Type, first_equation), make_pair(++c_Size, block_type_size_mfs[block_type_size_mfs.size()-1].second.second));
                  if (block_lag_lead[block_type_size_mfs.size()-1].first > Lag)
                    Lag = block_lag_lead[block_type_size_mfs.size()-1].first;
                  if (block_lag_lead[block_type_size_mfs.size()-1].second > Lead)
                    Lead = block_lag_lead[block_type_size_mfs.size()-1].second;
                  block_lag_lead[block_type_size_mfs.size()-1] = make_pair(Lag, Lead);
                }
              else
                {
                  block_type_size_mfs.push_back(make_pair(make_pair(Simulation_Type, eq), make_pair(Blck_Size, MFS_Size)));
                  block_lag_lead.push_back(make_pair(Lag, Lead));
                }
            }
          else
            {
              block_type_size_mfs.push_back(make_pair(make_pair(Simulation_Type, eq), make_pair(Blck_Size, MFS_Size)));
              block_lag_lead.push_back(make_pair(Lag, Lead));
            }
        }
      else
        {
          block_type_size_mfs.push_back(make_pair(make_pair(Simulation_Type, eq), make_pair(Blck_Size, MFS_Size)));
          block_lag_lead.push_back(make_pair(Lag, Lead));
        }
      prev_Type = Simulation_Type;
      eq += Blck_Size;
    }
  return (block_type_size_mfs);
}

vector<bool>
ModelTree::BlockLinear(t_blocks_derivatives &blocks_derivatives, vector<int> &variable_reordered)
{
  unsigned int nb_blocks = getNbBlocks();
  vector<bool> blocks_linear(nb_blocks, true);
  for (unsigned int block = 0; block < nb_blocks; block++)
    {
      BlockSimulationType simulation_type = getBlockSimulationType(block);
      int block_size = getBlockSize(block);
      t_block_derivatives_equation_variable_laglead_nodeid derivatives_block = blocks_derivatives[block];
      int first_variable_position = getBlockFirstEquation(block);
      if (simulation_type == SOLVE_BACKWARD_COMPLETE || simulation_type == SOLVE_FORWARD_COMPLETE)
        {
          for (t_block_derivatives_equation_variable_laglead_nodeid::const_iterator it = derivatives_block.begin(); it != derivatives_block.end(); it++)
            {
              int lag = it->second.first;
              if (lag == 0)
                {
                  NodeID Id = it->second.second;
                  set<pair<int, int> > endogenous;
                  Id->collectEndogenous(endogenous);
                  if (endogenous.size() > 0)
                    {
                      for (int l = 0; l < block_size; l++)
                        {
                          if (endogenous.find(make_pair(variable_reordered[first_variable_position+l], 0)) != endogenous.end())
                            {
                              blocks_linear[block] = false;
                              goto the_end;
                            }
                        }
                    }
                }
            }
        }
      else if (simulation_type == SOLVE_TWO_BOUNDARIES_COMPLETE || simulation_type == SOLVE_TWO_BOUNDARIES_SIMPLE)
        {
          for (t_block_derivatives_equation_variable_laglead_nodeid::const_iterator it = derivatives_block.begin(); it != derivatives_block.end(); it++)
            {
              int lag = it->second.first;
              NodeID Id = it->second.second; //
              set<pair<int, int> > endogenous;
              Id->collectEndogenous(endogenous);
              if (endogenous.size() > 0)
                {
                  for (int l = 0; l < block_size; l++)
                    {
                      if (endogenous.find(make_pair(variable_reordered[first_variable_position+l], lag)) != endogenous.end())
                        {
                          blocks_linear[block] = false;
                          goto the_end;
                        }
                    }
                }
            }
        }
    the_end:
      ;
    }
  return (blocks_linear);
}

ModelTree::ModelTree(SymbolTable &symbol_table_arg,
                     NumericalConstants &num_constants_arg) :
  DataTree(symbol_table_arg, num_constants_arg)
{
  for (int i = 0; i < 3; i++)
    NNZDerivatives[i] = 0;
}

int
ModelTree::equation_number() const
{
  return (equations.size());
}

void
ModelTree::writeDerivative(ostream &output, int eq, int symb_id, int lag,
                           ExprNodeOutputType output_type,
                           const temporary_terms_type &temporary_terms) const
{
  first_derivatives_type::const_iterator it = first_derivatives.find(make_pair(eq, getDerivID(symb_id, lag)));
  if (it != first_derivatives.end())
    (it->second)->writeOutput(output, output_type, temporary_terms);
  else
    output << 0;
}

void
ModelTree::computeJacobian(const set<int> &vars)
{
  for (set<int>::const_iterator it = vars.begin();
       it != vars.end(); it++)
    for (int eq = 0; eq < (int) equations.size(); eq++)
      {
        NodeID d1 = equations[eq]->getDerivative(*it);
        if (d1 == Zero)
          continue;
        first_derivatives[make_pair(eq, *it)] = d1;
        ++NNZDerivatives[0];
      }
}

void
ModelTree::computeHessian(const set<int> &vars)
{
  for (first_derivatives_type::const_iterator it = first_derivatives.begin();
       it != first_derivatives.end(); it++)
    {
      int eq = it->first.first;
      int var1 = it->first.second;
      NodeID d1 = it->second;

      // Store only second derivatives with var2 <= var1
      for (set<int>::const_iterator it2 = vars.begin();
           it2 != vars.end(); it2++)
        {
          int var2 = *it2;
          if (var2 > var1)
            continue;

          NodeID d2 = d1->getDerivative(var2);
          if (d2 == Zero)
            continue;
          second_derivatives[make_pair(eq, make_pair(var1, var2))] = d2;
          if (var2 == var1)
            ++NNZDerivatives[1];
          else
            NNZDerivatives[1] += 2;
        }
    }
}

void
ModelTree::computeThirdDerivatives(const set<int> &vars)
{
  for (second_derivatives_type::const_iterator it = second_derivatives.begin();
       it != second_derivatives.end(); it++)
    {
      int eq = it->first.first;

      int var1 = it->first.second.first;
      int var2 = it->first.second.second;
      // By construction, var2 <= var1

      NodeID d2 = it->second;

      // Store only third derivatives such that var3 <= var2 <= var1
      for (set<int>::const_iterator it2 = vars.begin();
           it2 != vars.end(); it2++)
        {
          int var3 = *it2;
          if (var3 > var2)
            continue;

          NodeID d3 = d2->getDerivative(var3);
          if (d3 == Zero)
            continue;
          third_derivatives[make_pair(eq, make_pair(var1, make_pair(var2, var3)))] = d3;
          if (var3 == var2 && var2 == var1)
            ++NNZDerivatives[2];
          else if (var3 == var2 || var2 == var1)
            NNZDerivatives[2] += 3;
          else
            NNZDerivatives[2] += 6;
        }
    }
}

void
ModelTree::computeTemporaryTerms(bool is_matlab)
{
  map<NodeID, int> reference_count;
  temporary_terms.clear();

  for (vector<BinaryOpNode *>::iterator it = equations.begin();
       it != equations.end(); it++)
    (*it)->computeTemporaryTerms(reference_count, temporary_terms, is_matlab);

  for (first_derivatives_type::iterator it = first_derivatives.begin();
       it != first_derivatives.end(); it++)
    it->second->computeTemporaryTerms(reference_count, temporary_terms, is_matlab);

  for (second_derivatives_type::iterator it = second_derivatives.begin();
       it != second_derivatives.end(); it++)
    it->second->computeTemporaryTerms(reference_count, temporary_terms, is_matlab);

  for (third_derivatives_type::iterator it = third_derivatives.begin();
       it != third_derivatives.end(); it++)
    it->second->computeTemporaryTerms(reference_count, temporary_terms, is_matlab);
}

void
ModelTree::writeTemporaryTerms(const temporary_terms_type &tt, ostream &output,
                               ExprNodeOutputType output_type) const
{
  // Local var used to keep track of temp nodes already written
  temporary_terms_type tt2;

  if (tt.size() > 0 && (IS_C(output_type)))
    output << "double" << endl;

  for (temporary_terms_type::const_iterator it = tt.begin();
       it != tt.end(); it++)
    {
      if (IS_C(output_type) && it != tt.begin())
        output << "," << endl;

      (*it)->writeOutput(output, output_type, tt);
      output << " = ";

      (*it)->writeOutput(output, output_type, tt2);

      // Insert current node into tt2
      tt2.insert(*it);

      if (IS_MATLAB(output_type))
        output << ";" << endl;
    }
  if (IS_C(output_type))
    output << ";" << endl;
}

void
ModelTree::compileTemporaryTerms(ostream &code_file, const temporary_terms_type &tt, map_idx_type map_idx, bool dynamic, bool steady_dynamic) const
{
  // Local var used to keep track of temp nodes already written
  temporary_terms_type tt2;
  for (temporary_terms_type::const_iterator it = tt.begin();
       it != tt.end(); it++)
    {
      (*it)->compile(code_file, false, tt2, map_idx, dynamic, steady_dynamic);
      if (dynamic)
        {
          FSTPT_ fstpt((int)(map_idx.find((*it)->idx)->second));
          fstpt.write(code_file);
        }
      else
        {
          FSTPST_ fstpst((int)(map_idx.find((*it)->idx)->second));
          fstpst.write(code_file);
        }
      // Insert current node into tt2
      tt2.insert(*it);
    }
}



void
ModelTree::writeModelLocalVariables(ostream &output, ExprNodeOutputType output_type) const
{
  for (map<int, NodeID>::const_iterator it = local_variables_table.begin();
       it != local_variables_table.end(); it++)
    {
      int id = it->first;
      NodeID value = it->second;

      if (IS_C(output_type))
        output << "double ";

      output << symbol_table.getName(id) << " = ";
      // Use an empty set for the temporary terms
      value->writeOutput(output, output_type, temporary_terms_type());
      output << ";" << endl;
    }
}

void
ModelTree::writeModelEquations(ostream &output, ExprNodeOutputType output_type) const
{
  for (int eq = 0; eq < (int) equations.size(); eq++)
    {
      BinaryOpNode *eq_node = equations[eq];
      NodeID lhs = eq_node->get_arg1();
      NodeID rhs = eq_node->get_arg2();

      // Test if the right hand side of the equation is empty.
      double vrhs = 1.0;
      try
        {
          vrhs = rhs->eval(eval_context_type());
        }
      catch (ExprNode::EvalException &e)
        {
        }

      if (vrhs != 0) // The right hand side of the equation is not empty ==> residual=lhs-rhs;
        {
          output << "lhs =";
          lhs->writeOutput(output, output_type, temporary_terms);
          output << ";" << endl;

          output << "rhs =";
          rhs->writeOutput(output, output_type, temporary_terms);
          output << ";" << endl;

          output << "residual" << LEFT_ARRAY_SUBSCRIPT(output_type)
                 << eq + ARRAY_SUBSCRIPT_OFFSET(output_type)
                 << RIGHT_ARRAY_SUBSCRIPT(output_type)
                 << "= lhs-rhs;" << endl;
        }
      else // The right hand side of the equation is empty ==> residual=lhs;
        {
          output << "residual" << LEFT_ARRAY_SUBSCRIPT(output_type)
                 << eq + ARRAY_SUBSCRIPT_OFFSET(output_type)
                 << RIGHT_ARRAY_SUBSCRIPT(output_type)
                 << " = ";
          lhs->writeOutput(output, output_type, temporary_terms);
          output << ";" << endl;
        }
    }
}

void
ModelTree::compileModelEquations(ostream &code_file, const temporary_terms_type &tt, map_idx_type &map_idx, bool dynamic, bool steady_dynamic) const
{
  for (int eq = 0; eq < (int) equations.size(); eq++)
    {
      BinaryOpNode *eq_node = equations[eq];
      NodeID lhs = eq_node->get_arg1();
      NodeID rhs = eq_node->get_arg2();

      // Test if the right hand side of the equation is empty.
      double vrhs = 1.0;
      try
        {
          vrhs = rhs->eval(eval_context_type());
        }
      catch (ExprNode::EvalException &e)
        {
        }

      if (vrhs != 0) // The right hand side of the equation is not empty ==> residual=lhs-rhs;
        {
          lhs->compile(code_file, false, temporary_terms, map_idx, dynamic, steady_dynamic);
          rhs->compile(code_file, false, temporary_terms, map_idx, dynamic, steady_dynamic);

          FBINARY_ fbinary(oMinus);
          fbinary.write(code_file);

          FSTPR_ fstpr(eq);
          fstpr.write(code_file);
        }
      else // The right hand side of the equation is empty ==> residual=lhs;
        {
          lhs->compile(code_file, false, temporary_terms, map_idx, dynamic, steady_dynamic);
          FSTPR_ fstpr(eq);
          fstpr.write(code_file);
        }
    }
}

void
ModelTree::Write_Inf_To_Bin_File(const string &basename,
                                   int &u_count_int, bool &file_open, bool is_two_boundaries, int block_mfs) const
{
  int j;
  std::ofstream SaveCode;
  const string bin_basename = basename + ".bin";
  if (file_open)
    SaveCode.open(bin_basename.c_str(), ios::out | ios::in | ios::binary | ios::ate);
  else
    SaveCode.open(bin_basename.c_str(), ios::out | ios::binary);
  if (!SaveCode.is_open())
    {
      cout << "Error : Can't open file \"" << bin_basename << "\" for writing\n";
      exit(EXIT_FAILURE);
    }
  u_count_int = 0;
  for (first_derivatives_type::const_iterator it = first_derivatives.begin();it != first_derivatives.end(); it++)
    {
      int deriv_id = it->first.second;
      if (getTypeByDerivID(deriv_id) == eEndogenous)
        {
          int eq = it->first.first;
          int symb = getSymbIDByDerivID(deriv_id);
          int var = symbol_table.getTypeSpecificID(symb);
          int lag = getLagByDerivID(deriv_id);
          SaveCode.write(reinterpret_cast<char *>(&eq), sizeof(eq));
          int varr = var + lag * block_mfs;
          SaveCode.write(reinterpret_cast<char *>(&varr), sizeof(varr));
          SaveCode.write(reinterpret_cast<char *>(&lag), sizeof(lag));
          int u = u_count_int + block_mfs;
          SaveCode.write(reinterpret_cast<char *>(&u), sizeof(u));
          u_count_int++;
        }
    }
  if (is_two_boundaries)
    u_count_int +=  symbol_table.endo_nbr();
  for (j = 0; j < (int) symbol_table.endo_nbr(); j++)
    SaveCode.write(reinterpret_cast<char *>(&j), sizeof(j));
  for (j = 0; j < (int) symbol_table.endo_nbr(); j++)
    SaveCode.write(reinterpret_cast<char *>(&j), sizeof(j));
  SaveCode.close();
}

void
ModelTree::writeLatexModelFile(const string &filename, ExprNodeOutputType output_type) const
{
  ofstream output;
  output.open(filename.c_str(), ios::out | ios::binary);
  if (!output.is_open())
    {
      cerr << "ERROR: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }

  output << "\\documentclass[10pt,a4paper]{article}" << endl
         << "\\usepackage[landscape]{geometry}" << endl
         << "\\usepackage{fullpage}" << endl
         << "\\begin{document}" << endl
         << "\\footnotesize" << endl;

  // Write model local variables
  for (map<int, NodeID>::const_iterator it = local_variables_table.begin();
       it != local_variables_table.end(); it++)
    {
      int id = it->first;
      NodeID value = it->second;

      output << "\\begin{equation*}" << endl
             << symbol_table.getName(id) << " = ";
      // Use an empty set for the temporary terms
      value->writeOutput(output, output_type, temporary_terms_type());
      output << endl << "\\end{equation*}" << endl;
    }

  for (int eq = 0; eq < (int) equations.size(); eq++)
    {
      output << "\\begin{equation}" << endl
             << "% Equation " << eq+1 << endl;
      equations[eq]->writeOutput(output, output_type, temporary_terms_type());
      output << endl << "\\end{equation}" << endl;
    }

  output << "\\end{document}" << endl;

  output.close();
}

void
ModelTree::addEquation(NodeID eq)
{
  BinaryOpNode *beq = dynamic_cast<BinaryOpNode *>(eq);
  assert(beq != NULL && beq->get_op_code() == oEqual);

  equations.push_back(beq);
}

void
ModelTree::addEquationTags(int i, const string &key, const string &value)
{
  equation_tags.push_back(make_pair(i, make_pair(key, value)));
}

void
ModelTree::addAuxEquation(NodeID eq)
{
  BinaryOpNode *beq = dynamic_cast<BinaryOpNode *>(eq);
  assert(beq != NULL && beq->get_op_code() == oEqual);

  aux_equations.push_back(beq);
}
