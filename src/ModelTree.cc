/*
 * Copyright (C) 2003-2018 Dynare Team
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
#include <cmath>
#include <iostream>
#include <fstream>

#include "ModelTree.hh"
#include "MinimumFeedbackSet.hh"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/max_cardinality_matching.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/topological_sort.hpp>

using namespace MFS;

void
ModelTree::copyHelper(const ModelTree &m)
{
  auto f = [this](expr_t e) { return e->clone(*this); };

  // Equations
  for (const auto & it : m.equations)
    equations.push_back(dynamic_cast<BinaryOpNode *>(f(it)));
  for (const auto & it : m.aux_equations)
    aux_equations.push_back(dynamic_cast<BinaryOpNode *>(f(it)));

  auto convert_deriv_map = [f](map<vector<int>, expr_t> dm)
    {
      map<vector<int>, expr_t> dm2;
      for (const auto &it : dm)
        dm2.emplace(it.first, f(it.second));
      return dm2;
    };

  // Derivatives
  for (const auto &it : m.derivatives)
    derivatives.push_back(convert_deriv_map(it));
  for (const auto &it : m.params_derivatives)
    params_derivatives[it.first] = convert_deriv_map(it.second);

  auto convert_temporary_terms_t = [f](temporary_terms_t tt)
    {
      temporary_terms_t tt2;
      for (const auto &it : tt)
        tt2.insert(f(it));
      return tt2;
    };

  // Temporary terms
  for (const auto & it : m.temporary_terms)
    temporary_terms.insert(f(it));
  for (const auto & it : m.temporary_terms_mlv)
    temporary_terms_mlv[f(it.first)] = f(it.second);
  for (const auto &it : m.temporary_terms_derivatives)
    temporary_terms_derivatives.push_back(convert_temporary_terms_t(it));
  for (const auto & it : m.temporary_terms_idxs)
    temporary_terms_idxs[f(it.first)] = it.second;
  for (const auto & it : m.params_derivs_temporary_terms)
    params_derivs_temporary_terms[it.first] = convert_temporary_terms_t(it.second);
  for (const auto & it : m.params_derivs_temporary_terms_idxs)
    params_derivs_temporary_terms_idxs[f(it.first)] = it.second;

  // Other stuff
  for (const auto & it : m.trend_symbols_map)
    trend_symbols_map[it.first] = f(it.second);
  for (const auto & it : m.nonstationary_symbols_map)
    nonstationary_symbols_map[it.first] = make_pair(it.second.first, f(it.second.second));
}

ModelTree::ModelTree(SymbolTable &symbol_table_arg,
                     NumericalConstants &num_constants_arg,
                     ExternalFunctionsTable &external_functions_table_arg,
                     bool is_dynamic_arg) :
  DataTree {symbol_table_arg, num_constants_arg, external_functions_table_arg, is_dynamic_arg},
  derivatives(4),
  NNZDerivatives(4, 0),
  temporary_terms_derivatives(4)
{
}

ModelTree::ModelTree(const ModelTree &m) :
  DataTree {m},
  equations_lineno {m.equations_lineno},
  equation_tags {m.equation_tags},
  NNZDerivatives {m.NNZDerivatives},
  equation_reordered {m.equation_reordered},
  variable_reordered {m.variable_reordered},
  inv_equation_reordered {m.inv_equation_reordered},
  inv_variable_reordered {m.inv_variable_reordered},
  is_equation_linear {m.is_equation_linear},
  endo2eq {m.endo2eq},
  epilogue {m.epilogue},
  prologue {m.prologue},
  block_lag_lead {m.block_lag_lead},
  cutoff {m.cutoff},
  mfs {m.mfs}
{
  copyHelper(m);
}

ModelTree &
ModelTree::operator=(const ModelTree &m)
{
  DataTree::operator=(m);

  equations.clear();
  equations_lineno = m.equations_lineno;
  aux_equations.clear();
  equation_tags = m.equation_tags;
  NNZDerivatives = m.NNZDerivatives;

  derivatives.clear();
  params_derivatives.clear();

  temporary_terms.clear();
  temporary_terms_mlv.clear();
  temporary_terms_derivatives.clear();
  params_derivs_temporary_terms.clear();
  params_derivs_temporary_terms_idxs.clear();

  trend_symbols_map.clear();
  nonstationary_symbols_map.clear();

  equation_reordered = m.equation_reordered;
  variable_reordered = m.variable_reordered;
  inv_equation_reordered = m.inv_equation_reordered;
  inv_variable_reordered = m.inv_variable_reordered;
  is_equation_linear = m.is_equation_linear;
  endo2eq = m.endo2eq;
  epilogue = m.epilogue;
  prologue = m.prologue;
  block_lag_lead = m.block_lag_lead;
  cutoff = m.cutoff;
  mfs = m.mfs;

  copyHelper(m);

  return *this;
}


bool
ModelTree::computeNormalization(const jacob_map_t &contemporaneous_jacobian, bool verbose)
{
  const int n = equations.size();

  assert(n == symbol_table.endo_nbr());

  using BipartiteGraph = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS>;

  /*
    Vertices 0 to n-1 are for endogenous (using type specific ID)
    Vertices n to 2*n-1 are for equations (using equation no.)
  */
  BipartiteGraph g(2 * n);

  // Fill in the graph
  set<pair<int, int>> endo;

  for (const auto & it : contemporaneous_jacobian)
    add_edge(it.first.first + n, it.first.second, g);

  // Compute maximum cardinality matching
  vector<int> mate_map(2*n);

#if 1
  bool check = checked_edmonds_maximum_cardinality_matching(g, &mate_map[0]);
#else // Alternative way to compute normalization, by giving an initial matching using natural normalizations
  fill(mate_map.begin(), mate_map.end(), boost::graph_traits<BipartiteGraph>::null_vertex());

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

  boost::edmonds_augmenting_path_finder<BipartiteGraph, int *, boost::property_map<BipartiteGraph, boost::vertex_index_t>::type> augmentor(g, &mate_map[0], get(boost::vertex_index, g));
  while (augmentor.augment_matching())
    {
    };

  augmentor.get_current_matching(&mate_map[0]);

  bool check = boost::maximum_cardinality_matching_verifier<BipartiteGraph, int *, boost::property_map<BipartiteGraph, boost::vertex_index_t>::type>::verify_matching(g, &mate_map[0], get(boost::vertex_index, g));
#endif

  assert(check);

#ifdef DEBUG
  for (int i = 0; i < n; i++)
    cout << "Endogenous " << symbol_table.getName(symbol_table.getID(eEndogenous, i))
         << " matched with equation " << (mate_map[i]-n+1) << endl;
#endif

  // Create the resulting map, by copying the n first elements of mate_map, and substracting n to them
  endo2eq.resize(equations.size());
  transform(mate_map.begin(), mate_map.begin() + n, endo2eq.begin(), [=](int i) { return i-n; });

#ifdef DEBUG
  multimap<int, int> natural_endo2eqs;
  computeNormalizedEquations(natural_endo2eqs);

  int n1 = 0, n2 = 0;

  for (int i = 0; i < symbol_table.endo_nbr(); i++)
    {
      if (natural_endo2eqs.count(i) == 0)
        continue;

      n1++;

      auto x = natural_endo2eqs.equal_range(i);
      if (find_if(x.first, x.second, [=](auto y) { return y.second == endo2eq[i]; }) == x.second)
        cout << "Natural normalization of variable " << symbol_table.getName(symbol_table.getID(SymbolType::endogenous, i))
             << " not used." << endl;
      else
        n2++;
    }

  cout << "Used " << n2 << " natural normalizations out of " << n1 << ", for a total of " << n << " equations." << endl;
#endif

  // Check if all variables are normalized
  auto it = find(mate_map.begin(), mate_map.begin() + n, boost::graph_traits<BipartiteGraph>::null_vertex());
  if (it != mate_map.begin() + n)
    {
      if (verbose)
        cerr << "ERROR: Could not normalize the model. Variable "
             << symbol_table.getName(symbol_table.getID(SymbolType::endogenous, it - mate_map.begin()))
             << " is not in the maximum cardinality matching." << endl;
      check = false;
    }
  return check;
}

void
ModelTree::computeNonSingularNormalization(jacob_map_t &contemporaneous_jacobian, double cutoff, jacob_map_t &static_jacobian, dynamic_jacob_map_t &dynamic_jacobian)
{
  bool check = false;

  cout << "Normalizing the model..." << endl;

  int n = equations.size();

  // compute the maximum value of each row of the contemporaneous Jacobian matrix
  //jacob_map normalized_contemporaneous_jacobian;
  jacob_map_t normalized_contemporaneous_jacobian(contemporaneous_jacobian);
  vector<double> max_val(n, 0.0);
  for (const auto &it : contemporaneous_jacobian)
    if (fabs(it.second) > max_val[it.first.first])
      max_val[it.first.first] = fabs(it.second);

  for (auto & iter : normalized_contemporaneous_jacobian)
    iter.second /= max_val[iter.first.first];

  //We start with the highest value of the cutoff and try to normalize the model
  double current_cutoff = 0.99999999;

  int suppressed = 0;
  while (!check && current_cutoff > 1e-19)
    {
      jacob_map_t tmp_normalized_contemporaneous_jacobian;
      int suppress = 0;
      for (auto & iter : normalized_contemporaneous_jacobian)
        if (fabs(iter.second) > max(current_cutoff, cutoff))
          tmp_normalized_contemporaneous_jacobian[{ iter.first.first, iter.first.second }] = iter.second;
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
      jacob_map_t tmp_normalized_contemporaneous_jacobian;
      set<pair<int, int>> endo;
      for (int i = 0; i < n; i++)
        {
          endo.clear();
          equations[i]->collectEndogenous(endo);
          for (const auto & it : endo)
            tmp_normalized_contemporaneous_jacobian[{ i, it.first }] = 1;
        }
      check = computeNormalization(tmp_normalized_contemporaneous_jacobian, true);
      if (check)
        {
          // Update the jacobian matrix
          for (const auto &it : tmp_normalized_contemporaneous_jacobian)
            {
              if (static_jacobian.find({ it.first.first, it.first.second }) == static_jacobian.end())
                static_jacobian[{ it.first.first, it.first.second }] = 0;
              if (dynamic_jacobian.find({ 0, it.first.first, it.first.second }) == dynamic_jacobian.end())
                dynamic_jacobian[{ 0, it.first.first, it.first.second }] = nullptr;
              if (contemporaneous_jacobian.find({ it.first.first, it.first.second }) == contemporaneous_jacobian.end())
                contemporaneous_jacobian[{ it.first.first, it.first.second }] = 0;
              try
                {
                  if (derivatives[1].find({ it.first.first, getDerivID(symbol_table.getID(SymbolType::endogenous, it.first.second), 0) }) == derivatives[1].end())
                    derivatives[1][{ it.first.first, getDerivID(symbol_table.getID(SymbolType::endogenous, it.first.second), 0) }] = Zero;
                }
              catch (DataTree::UnknownDerivIDException &e)
                {
                  cerr << "The variable " << symbol_table.getName(symbol_table.getID(SymbolType::endogenous, it.first.second))
                       << " does not appear at the current period (i.e. with no lead and no lag); this case is not handled by the 'block' option of the 'model' block." << endl;
                  exit(EXIT_FAILURE);
                }
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
  for (size_t i = 0; i < equations.size(); i++)
    {
      auto *lhs = dynamic_cast<VariableNode *>(equations[i]->arg1);
      if (lhs == nullptr)
        continue;

      int symb_id = lhs->symb_id;
      if (symbol_table.getType(symb_id) != SymbolType::endogenous)
        continue;

      set<pair<int, int>> endo;
      equations[i]->arg2->collectEndogenous(endo);
      if (endo.find({ symbol_table.getTypeSpecificID(symb_id), 0 }) != endo.end())
        continue;

      endo2eqs.emplace(symbol_table.getTypeSpecificID(symb_id), (int) i);
      cout << "Endogenous " << symbol_table.getName(symb_id) << " normalized in equation " << (i+1) << endl;
    }
}

void
ModelTree::evaluateAndReduceJacobian(const eval_context_t &eval_context, jacob_map_t &contemporaneous_jacobian, jacob_map_t &static_jacobian, dynamic_jacob_map_t &dynamic_jacobian, double cutoff, bool verbose)
{
  int nb_elements_contemparenous_Jacobian = 0;
  set<vector<int>> jacobian_elements_to_delete;
  for (const auto &it : derivatives[1])
    {
      int deriv_id = it.first[1];
      if (getTypeByDerivID(deriv_id) == SymbolType::endogenous)
        {
          expr_t Id = it.second;
          int eq = it.first[0];
          int symb = getSymbIDByDerivID(deriv_id);
          int var = symbol_table.getTypeSpecificID(symb);
          int lag = getLagByDerivID(deriv_id);
          double val = 0;
          try
            {
              val = Id->eval(eval_context);
            }
          catch (ExprNode::EvalExternalFunctionException &e)
            {
              val = 1;
            }
          catch (ExprNode::EvalException &e)
            {
              cerr << "ERROR: evaluation of Jacobian failed for equation " << eq+1 << " (line " << equations_lineno[eq] << ") and variable " << symbol_table.getName(symb) << "(" << lag << ") [" << symb << "] !" << endl;
              Id->writeOutput(cerr, ExprNodeOutputType::matlabDynamicModelSparse, temporary_terms, {});
              cerr << endl;
              exit(EXIT_FAILURE);
            }
          if (fabs(val) < cutoff)
            {
              if (verbose)
                cout << "the coefficient related to variable " << var << " with lag " << lag << " in equation " << eq << " is equal to " << val << " and is set to 0 in the incidence matrix (size=" << symbol_table.endo_nbr() << ")" << endl;
              jacobian_elements_to_delete.insert({ eq, deriv_id });
            }
          else
            {
              if (lag == 0)
                {
                  nb_elements_contemparenous_Jacobian++;
                  contemporaneous_jacobian[{ eq, var }] = val;
                }
              if (static_jacobian.find({ eq, var }) != static_jacobian.end())
                static_jacobian[{ eq, var }] += val;
              else
                static_jacobian[{ eq, var }] = val;
              dynamic_jacobian[{ lag, eq, var }] = Id;
            }
        }
    }

  // Get rid of the elements of the Jacobian matrix below the cutoff
  for (const auto & it : jacobian_elements_to_delete)
    derivatives[1].erase(it);

  if (jacobian_elements_to_delete.size() > 0)
    {
      cout << jacobian_elements_to_delete.size() << " elements among " << derivatives[1].size() << " in the incidence matrices are below the cutoff (" << cutoff << ") and are discarded" << endl
           << "The contemporaneous incidence matrix has " << nb_elements_contemparenous_Jacobian << " elements" << endl;
    }
}

vector<pair<int, int> >
ModelTree::select_non_linear_equations_and_variables(vector<bool> is_equation_linear, const dynamic_jacob_map_t &dynamic_jacobian, vector<int> &equation_reordered, vector<int> &variable_reordered,
                                                     vector<int> &inv_equation_reordered, vector<int> &inv_variable_reordered,
                                                     lag_lead_vector_t &equation_lag_lead, lag_lead_vector_t &variable_lag_lead,
                                                     vector<unsigned int> &n_static, vector<unsigned int> &n_forward, vector<unsigned int> &n_backward, vector<unsigned int> &n_mixed)
{
  vector<int> eq2endo(equations.size(), 0);
  /*equation_reordered.resize(equations.size());
  variable_reordered.resize(equations.size());*/
  unsigned int num = 0;
  for (auto it : endo2eq)
    if (!is_equation_linear[it])
      num++;
  vector<int> endo2block = vector<int>(endo2eq.size(), 1);
  vector<pair<set<int>, pair<set<int>, vector<int> > > > components_set(num);
  int i = 0, j = 0;
  for (auto it : endo2eq)
    {
      if (!is_equation_linear[it])
        {
          equation_reordered[i] = it;
          variable_reordered[i] = j;
          endo2block[j] = 0;
          components_set[endo2block[j]].first.insert(i);
          /*cout << " -----------------------------------------" << endl;
          cout.flush();
          cout << "equation_reordered[" << *it << "] = " << i << endl;
          cout.flush();
          cout << "variable_reordered[" << j << "] = " << i << endl;
          cout.flush();
          cout << "components_set[" << endo2block[j] << "].first[" << i << "] = " << i << endl;
          cout << "endo2block[" << j << "] = " << 0 << endl;
          cout.flush();
          */
          i++;
          j++;
        }
    }
/*  for (unsigned int j = 0; j < is_equation_linear.size() ; j++)
     cout << "endo2block[" << j << "] = " << endo2block[j] << endl;*/
  /*cout << "before getVariableLeadLAgByBlock\n";
  cout.flush();*/
  getVariableLeadLagByBlock(dynamic_jacobian, endo2block, endo2block.size(), equation_lag_lead, variable_lag_lead, equation_reordered, variable_reordered);
  n_static = vector<unsigned int>(endo2eq.size(), 0);
  n_forward = vector<unsigned int>(endo2eq.size(), 0);
  n_backward = vector<unsigned int>(endo2eq.size(), 0);
  n_mixed = vector<unsigned int>(endo2eq.size(), 0);
  for (unsigned int i = 0; i < endo2eq.size(); i++)
    {
      if      (variable_lag_lead[variable_reordered[i]].first != 0 && variable_lag_lead[variable_reordered[i]].second != 0)
        n_mixed[i]++;
      else if (variable_lag_lead[variable_reordered[i]].first == 0 && variable_lag_lead[variable_reordered[i]].second != 0)
        n_forward[i]++;
      else if (variable_lag_lead[variable_reordered[i]].first != 0 && variable_lag_lead[variable_reordered[i]].second == 0)
        n_backward[i]++;
      else if (variable_lag_lead[variable_reordered[i]].first == 0 && variable_lag_lead[variable_reordered[i]].second == 0)
        n_static[i]++;
    }
  cout.flush();
  int nb_endo = is_equation_linear.size();
  vector<pair<int, int>> blocks(1, make_pair(i, i));
  inv_equation_reordered.resize(nb_endo);
  inv_variable_reordered.resize(nb_endo);
  for (int i = 0; i < nb_endo; i++)
    {
      inv_variable_reordered[variable_reordered[i]] = i;
      inv_equation_reordered[equation_reordered[i]] = i;
    }
  return blocks;
}

bool
ModelTree::computeNaturalNormalization()
{
  bool bool_result = true;
  set<pair<int, int> > result;
  endo2eq.resize(equations.size());
  for (int eq = 0; eq < (int) equations.size(); eq++)
    if (!is_equation_linear[eq])
      {
        BinaryOpNode *eq_node = equations[eq];
        expr_t lhs = eq_node->arg1;
        result.clear();
        lhs->collectDynamicVariables(SymbolType::endogenous, result);
        if (result.size() == 1 && result.begin()->second == 0)
          {
            //check if the endogenous variable has not been already used in an other match !
             auto it = find(endo2eq.begin(), endo2eq.end(), result.begin()->first);
             if (it == endo2eq.end())
               endo2eq[result.begin()->first] = eq;
             else
              {
                bool_result = false;
                break;
              }
          }
      }
  return bool_result;
}

void
ModelTree::computePrologueAndEpilogue(const jacob_map_t &static_jacobian_arg, vector<int> &equation_reordered, vector<int> &variable_reordered)
{
  vector<int> eq2endo(equations.size(), 0);
  equation_reordered.resize(equations.size());
  variable_reordered.resize(equations.size());
  int n = equations.size();
  vector<bool> IM(n*n);
  int i = 0;
  for (auto it : endo2eq)
    {
      eq2endo[it] = i;
      equation_reordered[i] = i;
      variable_reordered[it] = i;
      i++;
    }
  if (cutoff == 0)
    {
      set<pair<int, int>> endo;
      for (int i = 0; i < n; i++)
        {
          endo.clear();
          equations[i]->collectEndogenous(endo);
          for (const auto & it : endo)
            IM[i * n + endo2eq[it.first]] = true;
        }
    }
  else
    for (const auto & it : static_jacobian_arg)
      IM[it.first.first * n + endo2eq[it.first.second]] = true;
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
}

equation_type_and_normalized_equation_t
ModelTree::equationTypeDetermination(const map<tuple<int, int, int>, expr_t> &first_order_endo_derivatives, const vector<int> &Index_Var_IM, const vector<int> &Index_Equ_IM, int mfs) const
{
  expr_t lhs;
  BinaryOpNode *eq_node;
  EquationType Equation_Simulation_Type;
  equation_type_and_normalized_equation_t V_Equation_Simulation_Type(equations.size());
  for (unsigned int i = 0; i < equations.size(); i++)
    {
      int eq = Index_Equ_IM[i];
      int var = Index_Var_IM[i];
      eq_node = equations[eq];
      lhs = eq_node->arg1;
      Equation_Simulation_Type = E_SOLVE;
      auto derivative = first_order_endo_derivatives.find({ eq, var, 0 });
      pair<bool, expr_t> res;
      if (derivative != first_order_endo_derivatives.end())
        {
          set<pair<int, int>> result;
          derivative->second->collectEndogenous(result);
          auto d_endo_variable = result.find({ var, 0 });
          //Determine whether the equation could be evaluated rather than to be solved
          if (lhs->isVariableNodeEqualTo(SymbolType::endogenous, Index_Var_IM[i], 0) && derivative->second->isNumConstNodeEqualTo(1))
            Equation_Simulation_Type = E_EVALUATE;
          else
            {
              vector<tuple<int, expr_t, expr_t>> List_of_Op_RHS;
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
      V_Equation_Simulation_Type[eq] = { Equation_Simulation_Type, dynamic_cast<BinaryOpNode *>(res.second) };
    }
  return V_Equation_Simulation_Type;
}

void
ModelTree::getVariableLeadLagByBlock(const dynamic_jacob_map_t &dynamic_jacobian, const vector<int> &components_set, int nb_blck_sim, lag_lead_vector_t &equation_lead_lag, lag_lead_vector_t &variable_lead_lag, const vector<int> &equation_reordered, const vector<int> &variable_reordered) const
{
  int nb_endo = symbol_table.endo_nbr();
  variable_lead_lag = lag_lead_vector_t(nb_endo, { 0, 0 });
  equation_lead_lag = lag_lead_vector_t(nb_endo, { 0, 0 });
  vector<int> variable_blck(nb_endo), equation_blck(nb_endo);
  for (int i = 0; i < nb_endo; i++)
    {
      if (i < (int) prologue)
        {
          variable_blck[variable_reordered[i]] = i;
          equation_blck[equation_reordered[i]] = i;
        }
      else if (i < (int) (components_set.size() + prologue))
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
  for (const auto & it : dynamic_jacobian)
    {
      int lag, j_1, i_1;
      tie(lag, j_1, i_1) = it.first;
      if (variable_blck[i_1] == equation_blck[j_1])
        {
          if (lag > variable_lead_lag[i_1].second)
            variable_lead_lag[i_1] = { variable_lead_lag[i_1].first, lag };
          if (lag < -variable_lead_lag[i_1].first)
            variable_lead_lag[i_1] = { -lag, variable_lead_lag[i_1].second };
          if (lag > equation_lead_lag[j_1].second)
            equation_lead_lag[j_1] = { equation_lead_lag[j_1].first, lag };
          if (lag < -equation_lead_lag[j_1].first)
            equation_lead_lag[j_1] = { -lag, equation_lead_lag[j_1].second };
        }
    }
}

void
ModelTree::computeBlockDecompositionAndFeedbackVariablesForEachBlock(const jacob_map_t &static_jacobian, const dynamic_jacob_map_t &dynamic_jacobian, vector<int> &equation_reordered, vector<int> &variable_reordered, vector<pair<int, int>> &blocks, const equation_type_and_normalized_equation_t &Equation_Type, bool verbose_, bool select_feedback_variable, int mfs, vector<int> &inv_equation_reordered, vector<int> &inv_variable_reordered, lag_lead_vector_t &equation_lag_lead, lag_lead_vector_t &variable_lag_lead, vector<unsigned int> &n_static, vector<unsigned int> &n_forward, vector<unsigned int> &n_backward, vector<unsigned int> &n_mixed) const
{
  int nb_var = variable_reordered.size();
  int n = nb_var - prologue - epilogue;

  AdjacencyList_t G2(n);

  // It is necessary to manually initialize vertex_index property since this graph uses listS and not vecS as underlying vertex container
  auto v_index = get(boost::vertex_index, G2);
  for (int i = 0; i < n; i++)
    put(v_index, vertex(i, G2), i);

  vector<int> reverse_equation_reordered(nb_var), reverse_variable_reordered(nb_var);

  for (int i = 0; i < nb_var; i++)
    {
      reverse_equation_reordered[equation_reordered[i]] = i;
      reverse_variable_reordered[variable_reordered[i]] = i;
    }
  jacob_map_t tmp_normalized_contemporaneous_jacobian;
  if (cutoff == 0)
    {
      set<pair<int, int>> endo;
      for (int i = 0; i < nb_var; i++)
        {
          endo.clear();
          equations[i]->collectEndogenous(endo);
          for (const auto & it : endo)
            tmp_normalized_contemporaneous_jacobian[{ i, it.first }] = 1;

        }
    }
  else
    tmp_normalized_contemporaneous_jacobian = static_jacobian;
  for (const auto &it : tmp_normalized_contemporaneous_jacobian)
    if (reverse_equation_reordered[it.first.first] >= (int) prologue && reverse_equation_reordered[it.first.first] < (int) (nb_var - epilogue)
        && reverse_variable_reordered[it.first.second] >= (int) prologue && reverse_variable_reordered[it.first.second] < (int) (nb_var - epilogue)
        && it.first.first != endo2eq[it.first.second])
      add_edge(vertex(reverse_equation_reordered[endo2eq[it.first.second]]-prologue, G2),
               vertex(reverse_equation_reordered[it.first.first]-prologue, G2),
               G2);

  vector<int> endo2block(num_vertices(G2)), discover_time(num_vertices(G2));
  boost::iterator_property_map<int *, boost::property_map<AdjacencyList_t, boost::vertex_index_t>::type, int, int &> endo2block_map(&endo2block[0], get(boost::vertex_index, G2));

  // Compute strongly connected components
  int num = strong_components(G2, endo2block_map);

  blocks = vector<pair<int, int>>(num, { 0, 0 });

  // Create directed acyclic graph associated to the strongly connected components
  using DirectedGraph = boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS>;
  DirectedGraph dag(num);

  for (unsigned int i = 0; i < num_vertices(G2); i++)
    {
      AdjacencyList_t::out_edge_iterator it_out, out_end;
      AdjacencyList_t::vertex_descriptor vi = vertex(i, G2);
      for (tie(it_out, out_end) = out_edges(vi, G2); it_out != out_end; ++it_out)
        {
          int t_b = endo2block_map[target(*it_out, G2)];
          int s_b = endo2block_map[source(*it_out, G2)];
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
  vector<tuple<set<int>, set<int>, vector<int>>> components_set(num);
  for (unsigned int i = 0; i < endo2block.size(); i++)
    {
      endo2block[i] = unordered2ordered[endo2block[i]];
      blocks[endo2block[i]].first++;
      get<0>(components_set[endo2block[i]]).insert(i);
    }

  getVariableLeadLagByBlock(dynamic_jacobian, endo2block, num, equation_lag_lead, variable_lag_lead, equation_reordered, variable_reordered);

  vector<int> tmp_equation_reordered(equation_reordered), tmp_variable_reordered(variable_reordered);
  int order = prologue;
  //Add a loop on vertices which could not be normalized or vertices related to lead variables => force those vertices to belong to the feedback set
  if (select_feedback_variable)
    {
      for (int i = 0; i < n; i++)
        if (Equation_Type[equation_reordered[i+prologue]].first == E_SOLVE
            || variable_lag_lead[variable_reordered[i+prologue]].second > 0
            || variable_lag_lead[variable_reordered[i+prologue]].first > 0
            || equation_lag_lead[equation_reordered[i+prologue]].second > 0
            || equation_lag_lead[equation_reordered[i+prologue]].first > 0
            || mfs == 0)
          add_edge(vertex(i, G2), vertex(i, G2), G2);
    }
  else
    {
      for (int i = 0; i < n; i++)
        if (Equation_Type[equation_reordered[i+prologue]].first == E_SOLVE || mfs == 0)
          add_edge(vertex(i, G2), vertex(i, G2), G2);
    }
  //Determines the dynamic structure of each equation
  n_static = vector<unsigned int>(prologue+num+epilogue, 0);
  n_forward = vector<unsigned int>(prologue+num+epilogue, 0);
  n_backward = vector<unsigned int>(prologue+num+epilogue, 0);
  n_mixed = vector<unsigned int>(prologue+num+epilogue, 0);

  for (int i = 0; i < (int) prologue; i++)
    {
      if      (variable_lag_lead[tmp_variable_reordered[i]].first != 0 && variable_lag_lead[tmp_variable_reordered[i]].second != 0)
        n_mixed[i]++;
      else if (variable_lag_lead[tmp_variable_reordered[i]].first == 0 && variable_lag_lead[tmp_variable_reordered[i]].second != 0)
        n_forward[i]++;
      else if (variable_lag_lead[tmp_variable_reordered[i]].first != 0 && variable_lag_lead[tmp_variable_reordered[i]].second == 0)
        n_backward[i]++;
      else if (variable_lag_lead[tmp_variable_reordered[i]].first == 0 && variable_lag_lead[tmp_variable_reordered[i]].second == 0)
        n_static[i]++;
    }
  //For each block, the minimum set of feedback variable is computed
  // and the non-feedback variables are reordered to get
  // a sub-recursive block without feedback variables

  for (int i = 0; i < num; i++)
    {
      AdjacencyList_t G = extract_subgraph(G2, get<0>(components_set[i]));
      set<int> feed_back_vertices;
      //Print(G);
      AdjacencyList_t G1 = Minimal_set_of_feedback_vertex(feed_back_vertices, G);
      auto v_index = get(boost::vertex_index, G);
      get<1>(components_set[i]) = feed_back_vertices;
      blocks[i].second = feed_back_vertices.size();
      vector<int> Reordered_Vertice;
      Reorder_the_recursive_variables(G, feed_back_vertices, Reordered_Vertice);

      //First we have the recursive equations conditional on feedback variables
      for (int j = 0; j < 4; j++)
        {
          for (int its : Reordered_Vertice)
            {
              bool something_done = false;
              if      (j == 2 && variable_lag_lead[tmp_variable_reordered[its +prologue]].first != 0 && variable_lag_lead[tmp_variable_reordered[its +prologue]].second != 0)
                {
                  n_mixed[prologue+i]++;
                  something_done = true;
                }
              else if (j == 3 && variable_lag_lead[tmp_variable_reordered[its +prologue]].first == 0 && variable_lag_lead[tmp_variable_reordered[its +prologue]].second != 0)
                {
                  n_forward[prologue+i]++;
                  something_done = true;
                }
              else if (j == 1 && variable_lag_lead[tmp_variable_reordered[its +prologue]].first != 0 && variable_lag_lead[tmp_variable_reordered[its +prologue]].second == 0)
                {
                  n_backward[prologue+i]++;
                  something_done = true;
                }
              else if (j == 0 && variable_lag_lead[tmp_variable_reordered[its +prologue]].first == 0 && variable_lag_lead[tmp_variable_reordered[its +prologue]].second == 0)
                {
                  n_static[prologue+i]++;
                  something_done = true;
                }
              if (something_done)
                {
                  equation_reordered[order] = tmp_equation_reordered[its+prologue];
                  variable_reordered[order] = tmp_variable_reordered[its+prologue];
                  order++;
                }
            }
        }
      get<2>(components_set[i]) = Reordered_Vertice;
      //Second we have the equations related to the feedback variables
      for (int j = 0; j < 4; j++)
        {
          for (int feed_back_vertice : feed_back_vertices)
            {
              bool something_done = false;
              if      (j == 2 && variable_lag_lead[tmp_variable_reordered[v_index[vertex(feed_back_vertice, G)]+prologue]].first != 0 && variable_lag_lead[tmp_variable_reordered[v_index[vertex(feed_back_vertice, G)]+prologue]].second != 0)
                {
                  n_mixed[prologue+i]++;
                  something_done = true;
                }
              else if (j == 3 && variable_lag_lead[tmp_variable_reordered[v_index[vertex(feed_back_vertice, G)]+prologue]].first == 0 && variable_lag_lead[tmp_variable_reordered[v_index[vertex(feed_back_vertice, G)]+prologue]].second != 0)
                {
                  n_forward[prologue+i]++;
                  something_done = true;
                }
              else if (j == 1 && variable_lag_lead[tmp_variable_reordered[v_index[vertex(feed_back_vertice, G)]+prologue]].first != 0 && variable_lag_lead[tmp_variable_reordered[v_index[vertex(feed_back_vertice, G)]+prologue]].second == 0)
                {
                  n_backward[prologue+i]++;
                  something_done = true;
                }
              else if (j == 0 && variable_lag_lead[tmp_variable_reordered[v_index[vertex(feed_back_vertice, G)]+prologue]].first == 0 && variable_lag_lead[tmp_variable_reordered[v_index[vertex(feed_back_vertice, G)]+prologue]].second == 0)
                {
                  n_static[prologue+i]++;
                  something_done = true;
                }
              if (something_done)
                {
                  equation_reordered[order] = tmp_equation_reordered[v_index[vertex(feed_back_vertice, G)]+prologue];
                  variable_reordered[order] = tmp_variable_reordered[v_index[vertex(feed_back_vertice, G)]+prologue];
                  order++;
                }
            }
        }
    }

  for (int i = 0; i < (int) epilogue; i++)
    {
      if      (variable_lag_lead[tmp_variable_reordered[prologue+n+i]].first != 0 && variable_lag_lead[tmp_variable_reordered[prologue+n+i]].second != 0)
        n_mixed[prologue+num+i]++;
      else if (variable_lag_lead[tmp_variable_reordered[prologue+n+i]].first == 0 && variable_lag_lead[tmp_variable_reordered[prologue+n+i]].second != 0)
        n_forward[prologue+num+i]++;
      else if (variable_lag_lead[tmp_variable_reordered[prologue+n+i]].first != 0 && variable_lag_lead[tmp_variable_reordered[prologue+n+i]].second == 0)
        n_backward[prologue+num+i]++;
      else if (variable_lag_lead[tmp_variable_reordered[prologue+n+i]].first == 0 && variable_lag_lead[tmp_variable_reordered[prologue+n+i]].second == 0)
        n_static[prologue+num+i]++;
    }

  inv_equation_reordered.resize(nb_var);
  inv_variable_reordered.resize(nb_var);
  for (int i = 0; i < nb_var; i++)
    {
      inv_variable_reordered[variable_reordered[i]] = i;
      inv_equation_reordered[equation_reordered[i]] = i;
    }
}

void
ModelTree::printBlockDecomposition(const vector<pair<int, int>> &blocks) const
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
              Nb_feedback_variable = getBlockMfs(block);
            }
        }
    }

  int Nb_RecursBlocks = Nb_TotalBlocks - Nb_SimulBlocks;
  cout << Nb_TotalBlocks << " block(s) found:" << endl
       << "  " << Nb_RecursBlocks << " recursive block(s) and " << Nb_SimulBlocks << " simultaneous block(s)." << endl
       << "  the largest simultaneous block has " << largest_block << " equation(s)" << endl
       << "                                 and " << Nb_feedback_variable << " feedback variable(s)." << endl;
}

block_type_firstequation_size_mfs_t
ModelTree::reduceBlocksAndTypeDetermination(const dynamic_jacob_map_t &dynamic_jacobian, vector<pair<int, int>> &blocks, const equation_type_and_normalized_equation_t &Equation_Type, const vector<int> &variable_reordered, const vector<int> &equation_reordered, vector<unsigned int> &n_static, vector<unsigned int> &n_forward, vector<unsigned int> &n_backward, vector<unsigned int> &n_mixed, vector<tuple<int, int, int, int>> &block_col_type, bool linear_decomposition)
{
  int i = 0;
  int count_equ = 0, blck_count_simult = 0;
  int Blck_Size, MFS_Size;
  int Lead, Lag;
  block_type_firstequation_size_mfs_t block_type_size_mfs;
  BlockSimulationType Simulation_Type, prev_Type = UNKNOWN;
  int eq = 0;
  unsigned int l_n_static = 0;
  unsigned int l_n_forward = 0;
  unsigned int l_n_backward = 0;
  unsigned int l_n_mixed = 0;
  for (i = 0; i < (int) (prologue+blocks.size()+epilogue); i++)
    {
      int first_count_equ = count_equ;
      if (i < (int) prologue)
        {
          Blck_Size = 1;
          MFS_Size = 1;
        }
      else if (i < (int) (prologue+blocks.size()))
        {
          Blck_Size = blocks[blck_count_simult].first;
          MFS_Size = blocks[blck_count_simult].second;
          blck_count_simult++;
        }
      else if (i < (int) (prologue+blocks.size()+epilogue))
        {
          Blck_Size = 1;
          MFS_Size = 1;
        }

      Lag = Lead = 0;
      set<pair<int, int>> endo;
      for (count_equ  = first_count_equ; count_equ  < Blck_Size+first_count_equ; count_equ++)
        {
          endo.clear();
          equations[equation_reordered[count_equ]]->collectEndogenous(endo);
          for (const auto & it : endo)
            {
              int curr_variable = it.first;
              int curr_lag = it.second;
              auto it1 = find(variable_reordered.begin()+first_count_equ, variable_reordered.begin()+(first_count_equ+Blck_Size), curr_variable);
              if (linear_decomposition)
                {
                  if (dynamic_jacobian.find({ curr_lag, equation_reordered[count_equ], curr_variable }) != dynamic_jacobian.end())
                    {
                      if (curr_lag > Lead)
                        Lead = curr_lag;
                      else if (-curr_lag > Lag)
                        Lag = -curr_lag;
                    }
                }
              else
                {
                  if (it1 != variable_reordered.begin()+(first_count_equ+Blck_Size))
                    if (dynamic_jacobian.find({ curr_lag, equation_reordered[count_equ], curr_variable }) != dynamic_jacobian.end())
                      {
                        if (curr_lag > Lead)
                          Lead = curr_lag;
                        else if (-curr_lag > Lag)
                          Lag = -curr_lag;
                      }
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
      l_n_static = n_static[i];
      l_n_forward = n_forward[i];
      l_n_backward = n_backward[i];
      l_n_mixed = n_mixed[i];
      if (Blck_Size == 1)
        {
          if (Equation_Type[equation_reordered[eq]].first == E_EVALUATE || Equation_Type[equation_reordered[eq]].first == E_EVALUATE_S)
            {
              if (Simulation_Type == SOLVE_BACKWARD_SIMPLE)
                Simulation_Type = EVALUATE_BACKWARD;
              else if (Simulation_Type == SOLVE_FORWARD_SIMPLE)
                Simulation_Type = EVALUATE_FORWARD;
            }
          if (i > 0)
            {
              bool is_lead = false, is_lag = false;
              int c_Size = get<2>(block_type_size_mfs[block_type_size_mfs.size()-1]);
              int first_equation = get<1>(block_type_size_mfs[block_type_size_mfs.size()-1]);
              if (c_Size > 0 && ((prev_Type ==  EVALUATE_FORWARD && Simulation_Type == EVALUATE_FORWARD && !is_lead)
                                 || (prev_Type ==  EVALUATE_BACKWARD && Simulation_Type == EVALUATE_BACKWARD && !is_lag)))
                {
                  for (int j = first_equation; j < first_equation+c_Size; j++)
                    {
                      auto it = dynamic_jacobian.find({ -1, equation_reordered[eq], variable_reordered[j] });
                      if (it != dynamic_jacobian.end())
                        is_lag = true;
                      it = dynamic_jacobian.find({ +1, equation_reordered[eq], variable_reordered[j] });
                      if (it != dynamic_jacobian.end())
                        is_lead = true;
                    }
                }
              if ((prev_Type ==  EVALUATE_FORWARD && Simulation_Type == EVALUATE_FORWARD && !is_lead)
                  || (prev_Type ==  EVALUATE_BACKWARD && Simulation_Type == EVALUATE_BACKWARD && !is_lag))
                {
                  //merge the current block with the previous one
                  BlockSimulationType c_Type = get<0>(block_type_size_mfs[block_type_size_mfs.size()-1]);
                  c_Size++;
                  block_type_size_mfs[block_type_size_mfs.size()-1] = { c_Type, first_equation, c_Size, c_Size };
                  if (block_lag_lead[block_type_size_mfs.size()-1].first > Lag)
                    Lag = block_lag_lead[block_type_size_mfs.size()-1].first;
                  if (block_lag_lead[block_type_size_mfs.size()-1].second > Lead)
                    Lead = block_lag_lead[block_type_size_mfs.size()-1].second;
                  block_lag_lead[block_type_size_mfs.size()-1] = { Lag, Lead };
                  auto tmp = block_col_type[block_col_type.size()-1];
                  block_col_type[block_col_type.size()-1] = { get<0>(tmp)+l_n_static, get<1>(tmp)+l_n_forward, get<2>(tmp)+l_n_backward, get<3>(tmp)+l_n_mixed };
                }
              else
                {
                  block_type_size_mfs.emplace_back(Simulation_Type, eq, Blck_Size, MFS_Size);
                  block_lag_lead.emplace_back(Lag, Lead);
                  block_col_type.emplace_back(l_n_static, l_n_forward, l_n_backward, l_n_mixed);
                }
            }
          else
            {
              block_type_size_mfs.emplace_back(Simulation_Type, eq, Blck_Size, MFS_Size);
              block_lag_lead.emplace_back(Lag, Lead);
              block_col_type.emplace_back(l_n_static, l_n_forward, l_n_backward, l_n_mixed);
            }
        }
      else
        {
          block_type_size_mfs.emplace_back(Simulation_Type, eq, Blck_Size, MFS_Size);
          block_lag_lead.emplace_back(Lag, Lead);
          block_col_type.emplace_back(l_n_static, l_n_forward, l_n_backward, l_n_mixed);
        }
      prev_Type = Simulation_Type;
      eq += Blck_Size;
    }
  return block_type_size_mfs;
}

vector<bool>
ModelTree::equationLinear(map<tuple<int, int, int>, expr_t> first_order_endo_derivatives) const
{
  vector<bool> is_linear(symbol_table.endo_nbr(), true);
  for (const auto &it : first_order_endo_derivatives)
    {
       expr_t Id = it.second;
       set<pair<int, int>> endogenous;
       Id->collectEndogenous(endogenous);
       if (endogenous.size() > 0)
         {
           int eq = get<0>(it.first);
           is_linear[eq] = false;
         }
    }
  return is_linear;
}

vector<bool>
ModelTree::BlockLinear(const blocks_derivatives_t &blocks_derivatives, const vector<int> &variable_reordered) const
{
  unsigned int nb_blocks = getNbBlocks();
  vector<bool> blocks_linear(nb_blocks, true);
  for (unsigned int block = 0; block < nb_blocks; block++)
    {
      BlockSimulationType simulation_type = getBlockSimulationType(block);
      int block_size = getBlockSize(block);
      block_derivatives_equation_variable_laglead_nodeid_t derivatives_block = blocks_derivatives[block];
      int first_variable_position = getBlockFirstEquation(block);
      if (simulation_type == SOLVE_BACKWARD_COMPLETE || simulation_type == SOLVE_FORWARD_COMPLETE)
        {
          for (const auto &it : derivatives_block)
            {
              int lag = get<2>(it);
              if (lag == 0)
                {
                  expr_t Id = get<3>(it);
                  set<pair<int, int>> endogenous;
                  Id->collectEndogenous(endogenous);
                  if (endogenous.size() > 0)
                    {
                      for (int l = 0; l < block_size; l++)
                        {
                          if (endogenous.find({ variable_reordered[first_variable_position+l], 0 }) != endogenous.end())
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
          for (const auto &it : derivatives_block)
            {
              int lag = get<2>(it);
              expr_t Id = get<3>(it); //
              set<pair<int, int>> endogenous;
              Id->collectEndogenous(endogenous);
              if (endogenous.size() > 0)
                {
                  for (int l = 0; l < block_size; l++)
                    {
                      if (endogenous.find({ variable_reordered[first_variable_position+l], lag }) != endogenous.end())
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
  return blocks_linear;
}

int
ModelTree::equation_number() const
{
  return (equations.size());
}

void
ModelTree::writeDerivative(ostream &output, int eq, int symb_id, int lag,
                           ExprNodeOutputType output_type,
                           const temporary_terms_t &temporary_terms) const
{
  auto it = derivatives[1].find({ eq, getDerivID(symb_id, lag) });
  if (it != derivatives[1].end())
    (it->second)->writeOutput(output, output_type, temporary_terms, {});
  else
    output << 0;
}

void
ModelTree::computeDerivatives(int order, const set<int> &vars)
{
  assert (order >= 1);

  // Do not shrink the vectors, since they have a minimal size of 4 (see constructor)
  derivatives.resize(max(static_cast<size_t>(order), derivatives.size()));
  NNZDerivatives.resize(max(static_cast<size_t>(order), NNZDerivatives.size()), 0);

  // First-order derivatives
  for (int var : vars)
    for (int eq = 0; eq < (int) equations.size(); eq++)
      {
        expr_t d1 = equations[eq]->getDerivative(var);
        if (d1 == Zero)
          continue;
        derivatives[1][{ eq, var }] = d1;
        ++NNZDerivatives[1];
      }

  // Higher-order derivatives
  for (int o = 2; o <= order; o++)
    for (const auto &it : derivatives[o-1])
      for (int var : vars)
        {
          if (it.first.back() > var)
            continue;

          expr_t d = it.second->getDerivative(var);
          if (d == Zero)
            continue;

          vector<int> indices{it.first};
          indices.push_back(var);
          // At this point, indices of endogenous variables are sorted in non-decreasing order
          derivatives[o][indices] = d;
          do
            NNZDerivatives[o]++;
          while (next_permutation(next(indices.begin()), indices.end()));
        }
}

void
ModelTree::computeTemporaryTerms(bool is_matlab, bool no_tmp_terms)
{
  /* Collect all model local variables appearing in equations (and only those,
     because printing unused model local variables can lead to a crash,
     see Dynare/dynare#101).
     Then store them in a dedicated structure (temporary_terms_mlv), that will
     be treated as the rest of temporary terms. */
  temporary_terms_mlv.clear();
  set<int> used_local_vars;
  for (auto & equation : equations)
    equation->collectVariables(SymbolType::modelLocalVariable, used_local_vars);
  for (int used_local_var : used_local_vars)
    {
      VariableNode *v = AddVariable(used_local_var);
      temporary_terms_mlv[v] = local_variables_table.find(used_local_var)->second;
    }

  // Compute the temporary terms in equations and derivatives
  map<pair<int, int>, temporary_terms_t> temp_terms_map;
  if (!no_tmp_terms)
    {
      map<expr_t, pair<int, pair<int, int>>> reference_count;

      for (auto & equation : equations)
        equation->computeTemporaryTerms({ 0, 0 },
                                        temp_terms_map,
                                        reference_count,
                                        is_matlab);

      for (int order = 1; order < (int) derivatives.size(); order++)
        for (const auto &it : derivatives[order])
          it.second->computeTemporaryTerms({ 0, order },
                                           temp_terms_map,
                                           reference_count,
                                           is_matlab);
    }

  // Fill the (now obsolete) temporary_terms structure
  temporary_terms.clear();
  for (const auto &it : temp_terms_map)
    temporary_terms.insert(it.second.begin(), it.second.end());

  // Fill the new structure
  temporary_terms_derivatives.clear();
  temporary_terms_derivatives.resize(derivatives.size());
  for (int order = 0; order < (int) derivatives.size(); order++)
    temporary_terms_derivatives[order] = move(temp_terms_map[{ 0, order }]);

  // Compute indices in MATLAB/Julia vector
  int idx = 0;
  for (auto &it : temporary_terms_mlv)
    temporary_terms_idxs[it.first] = idx++;
  for (int order = 0; order < (int) derivatives.size(); order++)
    for (const auto &it : temporary_terms_derivatives[order])
      temporary_terms_idxs[it] = idx++;
}

void
ModelTree::writeModelLocalVariableTemporaryTerms(temporary_terms_t &temp_term_union,
                                                 const temporary_terms_idxs_t &tt_idxs,
                                                 ostream &output, ExprNodeOutputType output_type,
                                                 deriv_node_temp_terms_t &tef_terms) const
{
  temporary_terms_t tto;
  for (auto it : temporary_terms_mlv)
    tto.insert(it.first);

  for (auto &it : temporary_terms_mlv)
    {
      if (isJuliaOutput(output_type))
        output << "    @inbounds const ";

      it.first->writeOutput(output, output_type, tto, tt_idxs, tef_terms);
      output << " = ";
      it.second->writeOutput(output, output_type, temp_term_union, tt_idxs, tef_terms);

      if (isCOutput(output_type) || isMatlabOutput(output_type))
        output << ";";
      output << endl;

      /* We put in temp_term_union the VariableNode corresponding to the MLV,
         not its definition, so that when equations use the MLV,
         T(XXX) is printed instead of the MLV name */
      temp_term_union.insert(it.first);
    }
}

void
ModelTree::writeTemporaryTerms(const temporary_terms_t &tt,
                               temporary_terms_t &temp_term_union,
                               const temporary_terms_idxs_t &tt_idxs,
                               ostream &output, ExprNodeOutputType output_type, deriv_node_temp_terms_t &tef_terms) const
{
  for (auto it : tt)
    {
      if (dynamic_cast<AbstractExternalFunctionNode *>(it) != nullptr)
        it->writeExternalFunctionOutput(output, output_type, temp_term_union, tt_idxs, tef_terms);

      if (isJuliaOutput(output_type))
        output << "    @inbounds ";

      it->writeOutput(output, output_type, tt, tt_idxs, tef_terms);
      output << " = ";
      it->writeOutput(output, output_type, temp_term_union, tt_idxs, tef_terms);

      if (isCOutput(output_type) || isMatlabOutput(output_type))
        output << ";";
      output << endl;

      temp_term_union.insert(it);
    }
}

void
ModelTree::writeJsonTemporaryTerms(const temporary_terms_t &tt, const temporary_terms_t &ttm1, ostream &output,
                                   deriv_node_temp_terms_t &tef_terms, string &concat) const
{
  // Local var used to keep track of temp nodes already written
  bool wrote_term = false;
  temporary_terms_t tt2 = ttm1;

  output << "\"external_functions_temporary_terms_" << concat << "\": [";
  for (auto it : tt)
    if (ttm1.find(it) == ttm1.end())
      {
        if (dynamic_cast<AbstractExternalFunctionNode *>(it) != nullptr)
          {
            if (wrote_term)
              output << ", ";
            vector<string> efout;
            it->writeJsonExternalFunctionOutput(efout, tt2, tef_terms);
            for (vector<string>::const_iterator it1 = efout.begin(); it1 != efout.end(); it1++)
              {
                if (it1 != efout.begin())
                  output << ", ";
                output << *it1;
              }
            wrote_term = true;
          }
        tt2.insert(it);
      }

  tt2 = ttm1;
  wrote_term = false;
  output << "]"
         << ", \"temporary_terms_" << concat << "\": [";
  for (auto it = tt.begin();
       it != tt.end(); it++)
    if (ttm1.find(*it) == ttm1.end())
      {
        if (wrote_term)
          output << ", ";
        output << "{\"temporary_term\": \"";
        (*it)->writeJsonOutput(output, tt, tef_terms);
        output << "\""
               << ", \"value\": \"";
        (*it)->writeJsonOutput(output, tt2, tef_terms);
        output << "\"}" << endl;
        wrote_term = true;

        // Insert current node into tt2
        tt2.insert(*it);
      }
  output << "]";
}

void
ModelTree::fixNestedParenthesis(ostringstream &output, map<string, string> &tmp_paren_vars, bool &message_printed) const
{
  string str = output.str();
  if (!testNestedParenthesis(str))
    return;
  int open = 0;
  int first_open_paren = 0;
  int matching_paren = 0;
  bool hit_limit = false;
  int i1 = 0;
  map<string, string>::iterator it;
  for (size_t i = 0; i < str.length(); i++)
    {
      if (str.at(i) == '(')
        {
          if (open == 0)
            first_open_paren = i;
          open++;
        }
      else if (str.at(i) == ')')
        {
          open--;
          if (open == 0)
            matching_paren = i;
        }
      if (open > 32)
        hit_limit = true;

      if (hit_limit && open == 0)
        {
          if (!message_printed)
            {
              cerr << "Warning: A .m file created by Dynare will have more than 32 nested parenthesis. Matlab cannot support this. " << endl
                   << "         We are going to modify, albeit inefficiently, this output to have fewer than 32 nested parenthesis. " << endl
                   << "         It would hence behoove you to use the use_dll option of the model block to circumnavigate this problem." << endl
                   << "         If you have not yet set up a compiler on your system, see the Matlab documentation for doing so." << endl
                   << "         For Windows, see: https://www.mathworks.com/help/matlab/matlab_external/install-mingw-support-package.html" << endl << endl;
              message_printed = true;
            }
          string str1 = str.substr(first_open_paren, matching_paren - first_open_paren + 1);
          string repstr = "";
          string varname;
          while (testNestedParenthesis(str1))
            {
              size_t open_paren_idx  = string::npos;
              size_t match_paren_idx = string::npos;
              size_t last_open_paren = string::npos;
              for (size_t j = 0; j < str1.length(); j++)
                {
                  if (str1.at(j) == '(')
                    {
                      // don't match, e.g. y(1)
                      size_t idx = str1.find_last_of("*/-+", j - 1);
                      if (j == 0 || (idx != string::npos && idx == j - 1))
                        open_paren_idx = j;
                      last_open_paren = j;
                    }
                  else if (str1.at(j) == ')')
                    {
                      // don't match, e.g. y(1)
                      size_t idx = str1.find_last_not_of("0123456789", j - 1);
                      if (idx != string::npos && idx != last_open_paren)
                        match_paren_idx = j;
                    }

                  if (open_paren_idx != string::npos && match_paren_idx != string::npos)
                    {
                      string val = str1.substr(open_paren_idx, match_paren_idx - open_paren_idx + 1);
                      it = tmp_paren_vars.find(val);
                      if (it == tmp_paren_vars.end())
                        {
                          ostringstream ptvstr;
                          ptvstr << i1++;
                          varname = "paren32_tmp_var_" + ptvstr.str();
                          repstr = repstr + varname + " = " + val + ";\n";
                          tmp_paren_vars[val] = varname;
                        }
                      else
                        varname = it->second;
                      str1.replace(open_paren_idx, match_paren_idx - open_paren_idx + 1, varname);
                      break;
                    }
                }
            }
          it = tmp_paren_vars.find(str1);
          if (it == tmp_paren_vars.end())
            {
              ostringstream ptvstr;
              ptvstr << i1++;
              varname = "paren32_tmp_var_" + ptvstr.str();
              repstr = repstr + varname + " = " + str1 + ";\n";
            }
          else
            varname = it->second;
          str.replace(first_open_paren, matching_paren - first_open_paren + 1, varname);
          size_t insertLoc = str.find_last_of("\n", first_open_paren);
          str.insert(insertLoc + 1, repstr);
          hit_limit = false;
          i = -1;
          first_open_paren = matching_paren = open = 0;
        }
    }
  output.str(str);
}

bool
ModelTree::testNestedParenthesis(const string &str) const
{
  int open = 0;
  for (char i : str)
    {
      if (i == '(')
        open++;
      else if (i == ')')
        open--;
      if (open > 32)
        return true;
    }
  return false;
}

void
ModelTree::compileTemporaryTerms(ostream &code_file, unsigned int &instruction_number, const temporary_terms_t &tt, map_idx_t map_idx, bool dynamic, bool steady_dynamic) const
{
  // Local var used to keep track of temp nodes already written
  temporary_terms_t tt2;
  // To store the functions that have already been written in the form TEF* = ext_fun();
  deriv_node_temp_terms_t tef_terms;
  for (auto it : tt)
    {
      if (dynamic_cast<AbstractExternalFunctionNode *>(it) != nullptr)
        {
          it->compileExternalFunctionOutput(code_file, instruction_number, false, tt2, map_idx, dynamic, steady_dynamic, tef_terms);
        }

      FNUMEXPR_ fnumexpr(TemporaryTerm, (int)(map_idx.find(it->idx)->second));
      fnumexpr.write(code_file, instruction_number);
      it->compile(code_file, instruction_number, false, tt2, map_idx, dynamic, steady_dynamic, tef_terms);
      if (dynamic)
        {
          FSTPT_ fstpt((int)(map_idx.find(it->idx)->second));
          fstpt.write(code_file, instruction_number);
        }
      else
        {
          FSTPST_ fstpst((int)(map_idx.find(it->idx)->second));
          fstpst.write(code_file, instruction_number);
        }
      // Insert current node into tt2
      tt2.insert(it);
    }
}

void
ModelTree::writeJsonModelLocalVariables(ostream &output, deriv_node_temp_terms_t &tef_terms) const
{
  /* Collect all model local variables appearing in equations, and print only
     them. Printing unused model local variables can lead to a crash (see
     ticket #101). */
  set<int> used_local_vars;

  // Use an empty set for the temporary terms
  const temporary_terms_t tt;

  for (auto equation : equations)
    equation->collectVariables(SymbolType::modelLocalVariable, used_local_vars);

  output << "\"model_local_variables\": [";
  bool printed = false;
  for (int it : local_variables_vector)
    if (used_local_vars.find(it) != used_local_vars.end())
      {
        if (printed)
          output << ", ";
        else
          printed = true;

        int id = it;
        vector<string> efout;
        expr_t value = local_variables_table.find(id)->second;
        value->writeJsonExternalFunctionOutput(efout, tt, tef_terms);
        for (vector<string>::const_iterator it1 = efout.begin(); it1 != efout.end(); it1++)
          {
            if (it1 != efout.begin())
              output << ", ";
            output << *it1;
          }

        if (!efout.empty())
          output << ", ";

        /* We append underscores to avoid name clashes with "g1" or "oo_" (see
           also VariableNode::writeOutput) */
        output << "{\"variable\": \"" << symbol_table.getName(id) << "__\""
               << ", \"value\": \"";
        value->writeJsonOutput(output, tt, tef_terms);
        output << "\"}" << endl;
      }
  output << "]";
}

void
ModelTree::writeModelEquations(ostream &output, ExprNodeOutputType output_type) const
{
  temporary_terms_t tt;
  temporary_terms_idxs_t ttidxs;
  writeModelEquations(output, output_type, tt);
}

void
ModelTree::writeModelEquations(ostream &output, ExprNodeOutputType output_type,
                               const temporary_terms_t &temporary_terms) const
{
  for (int eq = 0; eq < (int) equations.size(); eq++)
    {
      BinaryOpNode *eq_node = equations[eq];
      expr_t lhs = eq_node->arg1;
      expr_t rhs = eq_node->arg2;

      // Test if the right hand side of the equation is empty.
      double vrhs = 1.0;
      try
        {
          vrhs = rhs->eval(eval_context_t());
        }
      catch (ExprNode::EvalException &e)
        {
        }

      if (vrhs != 0) // The right hand side of the equation is not empty ==> residual=lhs-rhs;
        if (isJuliaOutput(output_type))
          {
            output << "    @inbounds residual" << LEFT_ARRAY_SUBSCRIPT(output_type)
                   << eq + ARRAY_SUBSCRIPT_OFFSET(output_type)
                   << RIGHT_ARRAY_SUBSCRIPT(output_type)
                   << " = (";
            lhs->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs);
            output << ") - (";
            rhs->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs);
            output << ")" << endl;
          }
        else
          {
            output << "lhs = ";
            lhs->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs);
            output << ";" << endl
                   << "rhs = ";
            rhs->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs);
            output << ";" << endl
                   << "residual" << LEFT_ARRAY_SUBSCRIPT(output_type)
                   << eq + ARRAY_SUBSCRIPT_OFFSET(output_type)
                   << RIGHT_ARRAY_SUBSCRIPT(output_type)
                   << " = lhs - rhs;" << endl;
          }
      else // The right hand side of the equation is empty ==> residual=lhs;
        {
          if (isJuliaOutput(output_type))
            output << "    @inbounds ";
          output << "residual" << LEFT_ARRAY_SUBSCRIPT(output_type)
                 << eq + ARRAY_SUBSCRIPT_OFFSET(output_type)
                 << RIGHT_ARRAY_SUBSCRIPT(output_type)
                 << " = ";
          lhs->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs);
          output << ";" << endl;
        }
    }
}

void
ModelTree::compileModelEquations(ostream &code_file, unsigned int &instruction_number, const temporary_terms_t &tt, const map_idx_t &map_idx, bool dynamic, bool steady_dynamic) const
{
  for (int eq = 0; eq < (int) equations.size(); eq++)
    {
      BinaryOpNode *eq_node = equations[eq];
      expr_t lhs = eq_node->arg1;
      expr_t rhs = eq_node->arg2;
      FNUMEXPR_ fnumexpr(ModelEquation, eq);
      fnumexpr.write(code_file, instruction_number);
      // Test if the right hand side of the equation is empty.
      double vrhs = 1.0;
      try
        {
          vrhs = rhs->eval(eval_context_t());
        }
      catch (ExprNode::EvalException &e)
        {
        }

      if (vrhs != 0) // The right hand side of the equation is not empty ==> residual=lhs-rhs;
        {
          lhs->compile(code_file, instruction_number, false, temporary_terms, map_idx, dynamic, steady_dynamic);
          rhs->compile(code_file, instruction_number, false, temporary_terms, map_idx, dynamic, steady_dynamic);

          FBINARY_ fbinary{static_cast<int>(BinaryOpcode::minus)};
          fbinary.write(code_file, instruction_number);

          FSTPR_ fstpr(eq);
          fstpr.write(code_file, instruction_number);
        }
      else // The right hand side of the equation is empty ==> residual=lhs;
        {
          lhs->compile(code_file, instruction_number, false, temporary_terms, map_idx, dynamic, steady_dynamic);
          FSTPR_ fstpr(eq);
          fstpr.write(code_file, instruction_number);
        }
    }
}

void
ModelTree::Write_Inf_To_Bin_File(const string &filename,
                                 int &u_count_int, bool &file_open, bool is_two_boundaries, int block_mfs) const
{
  int j;
  std::ofstream SaveCode;
  if (file_open)
    SaveCode.open(filename, ios::out | ios::in | ios::binary | ios::ate);
  else
    SaveCode.open(filename, ios::out | ios::binary);
  if (!SaveCode.is_open())
    {
      cerr << "Error : Can't open file \"" << filename << "\" for writing" << endl;
      exit(EXIT_FAILURE);
    }
  u_count_int = 0;
  for (const auto & first_derivative : derivatives[1])
    {
      int deriv_id = first_derivative.first[1];
      if (getTypeByDerivID(deriv_id) == SymbolType::endogenous)
        {
          int eq = first_derivative.first[0];
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
ModelTree::writeLatexModelFile(const string &basename, ExprNodeOutputType output_type, const bool write_equation_tags) const
{
  ofstream output, content_output;
  string filename = basename + ".tex";
  string content_basename = basename + "_content";
  string content_filename = content_basename + ".tex";
  output.open(filename, ios::out | ios::binary);
  if (!output.is_open())
    {
      cerr << "ERROR: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }

  content_output.open(content_filename, ios::out | ios::binary);
  if (!content_output.is_open())
    {
      cerr << "ERROR: Can't open file " << content_filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }

  output << "\\documentclass[10pt,a4paper]{article}" << endl
         << "\\usepackage[landscape]{geometry}" << endl
         << "\\usepackage{fullpage}" << endl
         << "\\usepackage{amsfonts}" << endl
         << "\\usepackage{breqn}" << endl
         << "\\begin{document}" << endl
         << "\\footnotesize" << endl;

  // Write model local variables
  for (int id : local_variables_vector)
    {
      expr_t value = local_variables_table.find(id)->second;

      content_output << "\\begin{dmath*}" << endl
                     << symbol_table.getTeXName(id) << " = ";
      // Use an empty set for the temporary terms
      value->writeOutput(content_output, output_type);
      content_output << endl << "\\end{dmath*}" << endl;
    }

  for (int eq = 0; eq < (int) equations.size(); eq++)
    {
      content_output << "% Equation " << eq + 1 << endl;
      if (write_equation_tags)
        {
          bool wrote_eq_tag = false;
          for (const auto & equation_tag : equation_tags)
            if (equation_tag.first == eq)
              {
                if (!wrote_eq_tag)
                  content_output << "\\noindent[";
                else
                  content_output << ", ";

                content_output << equation_tag.second.first;

                if (!(equation_tag.second.second.empty()))
                  content_output << "= `" << equation_tag.second.second << "'";

                wrote_eq_tag = true;
              }

          if (wrote_eq_tag)
            content_output << "]";
        }

      content_output << "\\begin{dmath}" << endl;
      // Here it is necessary to cast to superclass ExprNode, otherwise the overloaded writeOutput() method is not found
      dynamic_cast<ExprNode *>(equations[eq])->writeOutput(content_output, output_type);
      content_output << endl << "\\end{dmath}" << endl;
    }

  output << "\\include{" << content_basename << "}" << endl
         << "\\end{document}" << endl;

  output.close();
  content_output.close();
}

void
ModelTree::addEquation(expr_t eq, int lineno)
{
  auto *beq = dynamic_cast<BinaryOpNode *>(eq);
  assert(beq != nullptr && beq->op_code == BinaryOpcode::equal);

  equations.push_back(beq);
  equations_lineno.push_back(lineno);
}

void
ModelTree::addEquation(expr_t eq, int lineno, const vector<pair<string, string>> &eq_tags)
{
  int n = equations.size();
  for (const auto & eq_tag : eq_tags)
    equation_tags.emplace_back(n, eq_tag);
  addEquation(eq, lineno);
}

void
ModelTree::addAuxEquation(expr_t eq)
{
  auto *beq = dynamic_cast<BinaryOpNode *>(eq);
  assert(beq != nullptr && beq->op_code == BinaryOpcode::equal);

  aux_equations.push_back(beq);
}

void
ModelTree::addTrendVariables(vector<int> trend_vars, expr_t growth_factor) noexcept(false)
{
  while (!trend_vars.empty())
    if (trend_symbols_map.find(trend_vars.back()) != trend_symbols_map.end())
      throw TrendException(symbol_table.getName(trend_vars.back()));
    else
      {
        trend_symbols_map[trend_vars.back()] = growth_factor;
        trend_vars.pop_back();
      }
}

void
ModelTree::addNonstationaryVariables(vector<int> nonstationary_vars, bool log_deflator, expr_t deflator) noexcept(false)
{
  while (!nonstationary_vars.empty())
    if (nonstationary_symbols_map.find(nonstationary_vars.back()) != nonstationary_symbols_map.end())
      throw TrendException(symbol_table.getName(nonstationary_vars.back()));
    else
      {
        nonstationary_symbols_map[nonstationary_vars.back()] = { log_deflator, deflator };
        nonstationary_vars.pop_back();
      }
}

void
ModelTree::initializeVariablesAndEquations()
{
  for (size_t j = 0; j < equations.size(); j++)
    equation_reordered.push_back(j);

  for (int j = 0; j < symbol_table.endo_nbr(); j++)
    variable_reordered.push_back(j);
}

void
ModelTree::set_cutoff_to_zero()
{
  cutoff = 0;
}

void
ModelTree::jacobianHelper(ostream &output, int eq_nb, int col_nb, ExprNodeOutputType output_type) const
{
  if (isJuliaOutput(output_type))
    output << "    @inbounds ";
  output << "g1" << LEFT_ARRAY_SUBSCRIPT(output_type);
  if (isMatlabOutput(output_type) || isJuliaOutput(output_type))
    output << eq_nb + 1 << "," << col_nb + 1;
  else
    output << eq_nb + col_nb *equations.size();
  output << RIGHT_ARRAY_SUBSCRIPT(output_type);
}

void
ModelTree::sparseHelper(int order, ostream &output, int row_nb, int col_nb, ExprNodeOutputType output_type) const
{
  output << "v" << order << LEFT_ARRAY_SUBSCRIPT(output_type);
  if (isMatlabOutput(output_type) || isJuliaOutput(output_type))
    output << row_nb + 1 << "," << col_nb + 1;
  else
    output << row_nb + col_nb * NNZDerivatives[order];
  output << RIGHT_ARRAY_SUBSCRIPT(output_type);
}

void
ModelTree::computeParamsDerivatives(int paramsDerivsOrder)
{
  assert(paramsDerivsOrder >= 1);

  set<int> deriv_id_set;
  addAllParamDerivId(deriv_id_set);

  // First-order derivatives w.r.t. params
  for (int param : deriv_id_set)
    {
      for (int eq = 0; eq < (int) equations.size(); eq++)
        {
          expr_t d = equations[eq]->getDerivative(param);
          if (d == Zero)
            continue;
          params_derivatives[{ 0, 1 }][{ eq, param }] = d;
        }

      for (int endoOrd = 1; endoOrd < (int) derivatives.size(); endoOrd++)
        for (const auto &it : derivatives[endoOrd])
          {
            expr_t d = it.second->getDerivative(param);
            if (d == Zero)
              continue;
            vector<int> indices{it.first};
            indices.push_back(param);
            params_derivatives[{ endoOrd, 1 }][indices] = d;
          }
    }

  // Higher-order derivatives w.r.t. parameters
  for (int endoOrd = 0; endoOrd < (int) derivatives.size(); endoOrd++)
    for (int paramOrd = 2; paramOrd <= paramsDerivsOrder; paramOrd++)
      for (const auto &it : params_derivatives[{ endoOrd, paramOrd-1 }])
        for (int param : deriv_id_set)
          {
            if (it.first.back() > param)
              continue;

            expr_t d = it.second->getDerivative(param);
            if (d == Zero)
              continue;
            vector<int> indices{it.first};
            indices.push_back(param);
            // At this point, indices of both endogenous and parameters are sorted in non-decreasing order
            params_derivatives[{ endoOrd, paramOrd }][indices] = d;
          }
}

void
ModelTree::computeParamsDerivativesTemporaryTerms()
{
  map<expr_t, pair<int, pair<int, int>>> reference_count;

  /* The temp terms should be constructed in the same order as the for loops in
     {Static,Dynamic}Model::write{Json,}ParamsDerivativesFile() */
  params_derivs_temporary_terms.clear();
  for (const auto &it : params_derivatives)
    for (const auto &it2 : it.second)
      it2.second->computeTemporaryTerms(it.first,
                                        params_derivs_temporary_terms,
                                        reference_count,
                                        true);

  int idx = 0;
  for (auto &it : temporary_terms_mlv)
    params_derivs_temporary_terms_idxs[it.first] = idx++;
  for (const auto &it : params_derivs_temporary_terms)
    for (const auto &tt : it.second)
      params_derivs_temporary_terms_idxs[tt] = idx++;
}

bool
ModelTree::isNonstationary(int symb_id) const
{
  return (nonstationary_symbols_map.find(symb_id)
          != nonstationary_symbols_map.end());
}

void
ModelTree::writeJsonModelEquations(ostream &output, bool residuals) const
{
  vector<pair<string, string>> eqtags;
  temporary_terms_t tt_empty;
  if (residuals)
    output << endl << "\"residuals\":[" << endl;
  else
    output << endl << "\"model\":[" << endl;
  for (int eq = 0; eq < (int) equations.size(); eq++)
    {
      if (eq > 0)
        output << ", ";

      BinaryOpNode *eq_node = equations[eq];
      expr_t lhs = eq_node->arg1;
      expr_t rhs = eq_node->arg2;

      if (residuals)
        {
          output << "{\"residual\": {"
                 << "\"lhs\": \"";
          lhs->writeJsonOutput(output, temporary_terms, {});
          output << "\"";

          output << ", \"rhs\": \"";
          rhs->writeJsonOutput(output, temporary_terms, {});
          output << "\"";
          try
            {
              // Test if the right hand side of the equation is empty.
              if (rhs->eval(eval_context_t()) != 0)
                {
                  output << ", \"rhs\": \"";
                  rhs->writeJsonOutput(output, temporary_terms, {});
                  output << "\"";
                }
            }
          catch (ExprNode::EvalException &e)
            {
            }
          output << "}";
        }
      else
        {
          output << "{\"lhs\": \"";
          lhs->writeJsonOutput(output, tt_empty, {});
          output << "\", \"rhs\": \"";
          rhs->writeJsonOutput(output, tt_empty, {});
          output << "\""
                 << ", \"line\": " << equations_lineno[eq];

          for (const auto & equation_tag : equation_tags)
            if (equation_tag.first == eq)
              eqtags.push_back(equation_tag.second);

          if (!eqtags.empty())
            {
              output << ", \"tags\": {";
              int i = 0;
              for (vector<pair<string, string>>::const_iterator it = eqtags.begin(); it != eqtags.end(); it++, i++)
                {
                  if (i != 0)
                    output << ", ";
                  output << "\"" << it->first << "\": \"" << it->second << "\"";
                }
              output << "}";
              eqtags.clear();
            }
        }
      output << "}" << endl;
    }
  output << endl << "]" << endl;
}

string
ModelTree::matlab_arch(const string &mexext)
{
  if (mexext == "mexglx")
    return "glnx86";
  else if (mexext == "mexa64")
    return "glnxa64";
  if (mexext == "mexw32")
    return "win32";
  else if (mexext == "mexw64")
    return "win64";
  else if (mexext == "mexmaci")
    return "maci";
  else if (mexext == "mexmaci64")
    return "maci64";
  else
    {
      cerr << "ERROR: 'mexext' option to preprocessor incorrectly set, needed with 'use_dll'" << endl;
      exit(EXIT_FAILURE);
    }
}

void
ModelTree::compileDll(const string &basename, const string &static_or_dynamic, const string &mexext, const boost::filesystem::path &matlabroot, const boost::filesystem::path &dynareroot)
{
  const string opt_flags = "-O3 -g0 --param ira-max-conflict-table-size=1 -fno-forward-propagate -fno-gcse -fno-dce -fno-dse -fno-tree-fre -fno-tree-pre -fno-tree-cselim -fno-tree-dse -fno-tree-dce -fno-tree-pta -fno-gcse-after-reload";

  boost::filesystem::path compiler;
  ostringstream flags;
  string libs;

  if (mexext == "mex")
    {
      // Octave
      compiler = matlabroot / "bin" / "mkoctfile";
      flags << "--mex";
    }
  else
    {
      // MATLAB
      compiler = "gcc";
      string arch = matlab_arch(mexext);
      auto include_dir = matlabroot / "extern" / "include";
      flags << "-I " << include_dir;
      auto bin_dir = matlabroot / "bin" / arch;
      flags << " -L " << bin_dir;
      flags << " -fexceptions -DNDEBUG";
      libs = "-lmex -lmx -lmat -lmwlapack -lmwblas";
      if (mexext == "mexglx" || mexext == "mexa64")
        {
          // GNU/Linux
          flags << " -D_GNU_SOURCE -fPIC -pthread"
                << " -shared -Wl,--no-undefined -Wl,-rpath-link," << bin_dir;
          libs += " -lm -lstdc++";

          if (mexext == "mexglx")
            flags << " -D_FILE_OFFSET_BITS=64 -m32";
          else
            flags << " -fno-omit-frame-pointer";
        }
      else if (mexext == "mexw32" || mexext == "mexw64")
        {
          // Windows
          flags << " -static-libgcc -static-libstdc++ -shared";
          // Put the MinGW environment shipped with Dynare in the path
          boost::filesystem::path mingwpath = dynareroot / (string{"mingw"} + (mexext == "mexw32" ? "32" : "64")) / "bin";
          string newpath = "PATH=" + mingwpath.string() + ';' + string{getenv("PATH")};
          if (putenv(const_cast<char *>(newpath.c_str())) != 0)
            {
              cerr << "Can't set PATH" << endl;
              exit(EXIT_FAILURE);
            }
          }
      else
        {
          // macOS
          string archs = (mexext == "maci" ? "i386" : "x86_64");
          flags << " -fno-common -arch " << archs << " -mmacosx-version-min=10.7 -Wl,-twolevel_namespace -undefined error -bundle";
          libs += " -lm -lstdc++";
        }
    }

  auto model_dir = boost::filesystem::path{basename} / "model" / "src";
  boost::filesystem::path main_src{model_dir / (static_or_dynamic + ".c")},
    mex_src{model_dir / (static_or_dynamic + "_mex.c")};

  boost::filesystem::path mex_dir{"+" + basename};
  boost::filesystem::path binary{mex_dir / (static_or_dynamic + "." + mexext)};

  ostringstream cmd;

#ifdef _WIN32
  /* On Windows, system() hands the command over to "cmd.exe /C". We need to
     enclose the whole command line within double quotes if we want the inner
     quotes to be correctly handled. See "cmd /?" for more details. */
  cmd << '"';
#endif
  cmd << compiler << " " << opt_flags << " " << flags.str() << " " << main_src << " " << mex_src << " -o " << binary << " " << libs;

#ifdef _WIN32
  cmd << '"';
#endif

  cout << "Compiling " << static_or_dynamic << " MEX..." << endl << cmd.str() << endl;

  if (system(cmd.str().c_str()))
    {
      cerr << "Compilation failed" << endl;
      exit(EXIT_FAILURE);
    }
}
