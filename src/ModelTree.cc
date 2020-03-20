/*
 * Copyright © 2003-2020 Dynare Team
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

#include "ModelTree.hh"
#include "MinimumFeedbackSet.hh"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/max_cardinality_matching.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/topological_sort.hpp>
#pragma GCC diagnostic pop

#ifdef __APPLE__
# include <mach-o/dyld.h>
#endif

#include <regex>

void
ModelTree::copyHelper(const ModelTree &m)
{
  auto f = [this](expr_t e) { return e->clone(*this); };
  auto convert_vector_tt = [f](vector<temporary_terms_t> vtt)
                           {
                             vector<temporary_terms_t> vtt2;
                             for (const auto &tt : vtt)
                               {
                                 temporary_terms_t tt2;
                                 for (const auto &it : tt)
                                   tt2.insert(f(it));
                                 vtt2.push_back(tt2);
                               }
                             return vtt2;
                           };

  // Equations
  for (const auto &it : m.equations)
    equations.push_back(dynamic_cast<BinaryOpNode *>(f(it)));
  for (const auto &it : m.aux_equations)
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
  for (const auto &it : m.temporary_terms)
    temporary_terms.insert(f(it));
  for (const auto &it : m.temporary_terms_mlv)
    temporary_terms_mlv[f(it.first)] = f(it.second);
  for (const auto &it : m.temporary_terms_derivatives)
    temporary_terms_derivatives.push_back(convert_temporary_terms_t(it));
  for (const auto &it : m.temporary_terms_idxs)
    temporary_terms_idxs[f(it.first)] = it.second;
  for (const auto &it : m.v_temporary_terms)
    v_temporary_terms.push_back(convert_vector_tt(it));
  for (const auto &it : m.params_derivs_temporary_terms)
    params_derivs_temporary_terms[it.first] = convert_temporary_terms_t(it.second);
  for (const auto &it : m.params_derivs_temporary_terms_idxs)
    params_derivs_temporary_terms_idxs[f(it.first)] = it.second;

  // Other stuff
  for (const auto &it : m.trend_symbols_map)
    trend_symbols_map[it.first] = f(it.second);
  for (const auto &it : m.nonstationary_symbols_map)
    nonstationary_symbols_map[it.first] = {it.second.first, f(it.second.second)};
  for (const auto &it : m.dynamic_jacobian)
    dynamic_jacobian[it.first] = f(it.second);
  for (const auto &it : m.first_chain_rule_derivatives)
    first_chain_rule_derivatives[it.first] = f(it.second);

  for (const auto &it : m.equation_type_and_normalized_equation)
    equation_type_and_normalized_equation.emplace_back(it.first, f(it.second));

  for (const auto &it : m.blocks_derivatives)
    {
      block_derivatives_equation_variable_laglead_nodeid_t v;
      for (const auto &it2 : it)
        v.emplace_back(get<0>(it2), get<1>(it2), get<2>(it2), f(get<3>(it2)));
      blocks_derivatives.push_back(v);
    }

  auto convert_derivative_t = [f](derivative_t dt)
                              {
                                derivative_t dt2;
                                for (const auto &it : dt)
                                  dt2[it.first] = f(it.second);
                                return dt2;
                              };
  for (const auto &it : m.derivative_endo)
    derivative_endo.push_back(convert_derivative_t(it));
  for (const auto &it : m.derivative_other_endo)
    derivative_other_endo.push_back(convert_derivative_t(it));
  for (const auto &it : m.derivative_exo)
    derivative_exo.push_back(convert_derivative_t(it));
  for (const auto &it : m.derivative_exo_det)
    derivative_exo_det.push_back(convert_derivative_t(it));
}

ModelTree::ModelTree(SymbolTable &symbol_table_arg,
                     NumericalConstants &num_constants_arg,
                     ExternalFunctionsTable &external_functions_table_arg,
                     bool is_dynamic_arg) :
  DataTree{symbol_table_arg, num_constants_arg, external_functions_table_arg, is_dynamic_arg},
  derivatives(4),
  NNZDerivatives(4, 0),
  temporary_terms_derivatives(4)
{
}

ModelTree::ModelTree(const ModelTree &m) :
  DataTree{m},
  user_set_add_flags{m.user_set_add_flags},
  user_set_subst_flags{m.user_set_subst_flags},
  user_set_add_libs{m.user_set_add_libs},
  user_set_subst_libs{m.user_set_subst_libs},
  user_set_compiler{m.user_set_compiler},
  equations_lineno{m.equations_lineno},
  equation_tags{m.equation_tags},
  computed_derivs_order{m.computed_derivs_order},
  NNZDerivatives{m.NNZDerivatives},
  v_temporary_terms_inuse{m.v_temporary_terms_inuse},
  equation_reordered{m.equation_reordered},
  variable_reordered{m.variable_reordered},
  inv_equation_reordered{m.inv_equation_reordered},
  inv_variable_reordered{m.inv_variable_reordered},
  map_idx{m.map_idx},
  block_type_firstequation_size_mfs{m.block_type_firstequation_size_mfs},
  blocks_linear{m.blocks_linear},
  block_col_type{m.block_col_type},
  endo_max_leadlag_block{m.endo_max_leadlag_block},
  other_endo_max_leadlag_block{m.other_endo_max_leadlag_block},
  exo_max_leadlag_block{m.exo_max_leadlag_block},
  exo_det_max_leadlag_block{m.exo_det_max_leadlag_block},
  max_leadlag_block{m.max_leadlag_block},
  is_equation_linear{m.is_equation_linear},
  endo2eq{m.endo2eq},
  epilogue{m.epilogue},
  prologue{m.prologue},
  block_lag_lead{m.block_lag_lead},
  cutoff{m.cutoff},
  mfs{m.mfs}
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
  computed_derivs_order = m.computed_derivs_order;
  NNZDerivatives = m.NNZDerivatives;

  derivatives.clear();
  params_derivatives.clear();

  temporary_terms.clear();
  temporary_terms_mlv.clear();
  temporary_terms_derivatives.clear();
  v_temporary_terms.clear();
  v_temporary_terms_inuse = m.v_temporary_terms_inuse;
  params_derivs_temporary_terms.clear();
  params_derivs_temporary_terms_idxs.clear();

  trend_symbols_map.clear();
  nonstationary_symbols_map.clear();

  dynamic_jacobian.clear();
  equation_reordered = m.equation_reordered;
  variable_reordered = m.variable_reordered;
  inv_equation_reordered = m.inv_equation_reordered;
  inv_variable_reordered = m.inv_variable_reordered;
  first_chain_rule_derivatives.clear();
  map_idx = m.map_idx;
  equation_type_and_normalized_equation.clear();
  block_type_firstequation_size_mfs = m.block_type_firstequation_size_mfs;
  blocks_derivatives.clear();
  blocks_linear = m.blocks_linear;
  derivative_endo.clear();
  derivative_other_endo.clear();
  derivative_exo.clear();
  derivative_exo_det.clear();
  block_col_type = m.block_col_type;
  endo_max_leadlag_block = m.endo_max_leadlag_block;
  other_endo_max_leadlag_block = m.other_endo_max_leadlag_block;
  exo_max_leadlag_block = m.exo_max_leadlag_block;
  exo_det_max_leadlag_block = m.exo_det_max_leadlag_block;
  max_leadlag_block = m.max_leadlag_block;

  is_equation_linear = m.is_equation_linear;
  endo2eq = m.endo2eq;
  epilogue = m.epilogue;
  prologue = m.prologue;
  block_lag_lead = m.block_lag_lead;
  cutoff = m.cutoff;
  mfs = m.mfs;

  user_set_add_flags = m.user_set_add_flags;
  user_set_subst_flags = m.user_set_subst_flags;
  user_set_add_libs = m.user_set_add_libs;
  user_set_subst_libs = m.user_set_subst_libs;
  user_set_compiler = m.user_set_compiler;

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

  for (const auto &it : contemporaneous_jacobian)
    add_edge(it.first.first + n, it.first.second, g);

  // Compute maximum cardinality matching
  vector<int> mate_map(2*n);

#if 1
  bool check = checked_edmonds_maximum_cardinality_matching(g, &mate_map[0]);
#else // Alternative way to compute normalization, by giving an initial matching using natural normalizations
  fill(mate_map.begin(), mate_map.end(), boost::graph_traits<BipartiteGraph>::null_vertex());

  auto natural_endo2eqs = computeNormalizedEquations();

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
  auto natural_endo2eqs = computeNormalizedEquations(natural_endo2eqs);

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
  if (auto it = find(mate_map.begin(), mate_map.begin() + n, boost::graph_traits<BipartiteGraph>::null_vertex());
      it != mate_map.begin() + n)
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
ModelTree::computeNonSingularNormalization(jacob_map_t &contemporaneous_jacobian, double cutoff, jacob_map_t &static_jacobian)
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

  for (auto &iter : normalized_contemporaneous_jacobian)
    iter.second /= max_val[iter.first.first];

  //We start with the highest value of the cutoff and try to normalize the model
  double current_cutoff = 0.99999999;

  int suppressed = 0;
  while (!check && current_cutoff > 1e-19)
    {
      jacob_map_t tmp_normalized_contemporaneous_jacobian;
      int suppress = 0;
      for (auto &iter : normalized_contemporaneous_jacobian)
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
          for (const auto &it : endo)
            tmp_normalized_contemporaneous_jacobian[{ i, it.first }] = 1;
        }
      check = computeNormalization(tmp_normalized_contemporaneous_jacobian, true);
      if (check)
        {
          // Update the jacobian matrix
          for (const auto &[key, ignore] : tmp_normalized_contemporaneous_jacobian)
            {
              if (static_jacobian.find({ key.first, key.second }) == static_jacobian.end())
                static_jacobian[{ key.first, key.second }] = 0;
              if (dynamic_jacobian.find({ 0, key.first, key.second }) == dynamic_jacobian.end())
                dynamic_jacobian[{ 0, key.first, key.second }] = nullptr;
              if (contemporaneous_jacobian.find({ key.first, key.second }) == contemporaneous_jacobian.end())
                contemporaneous_jacobian[{ key.first, key.second }] = 0;
              try
                {
                  if (derivatives[1].find({ key.first, getDerivID(symbol_table.getID(SymbolType::endogenous, key.second), 0) }) == derivatives[1].end())
                    derivatives[1][{ key.first, getDerivID(symbol_table.getID(SymbolType::endogenous, key.second), 0) }] = Zero;
                }
              catch (DataTree::UnknownDerivIDException &e)
                {
                  cerr << "The variable " << symbol_table.getName(symbol_table.getID(SymbolType::endogenous, key.second))
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

multimap<int, int>
ModelTree::computeNormalizedEquations() const
{
  multimap<int, int> endo2eqs;
  for (size_t i = 0; i < equations.size(); i++)
    {
      auto lhs = dynamic_cast<VariableNode *>(equations[i]->arg1);
      if (!lhs)
        continue;

      int symb_id = lhs->symb_id;
      if (symbol_table.getType(symb_id) != SymbolType::endogenous)
        continue;

      set<pair<int, int>> endo;
      equations[i]->arg2->collectEndogenous(endo);
      if (endo.find({ symbol_table.getTypeSpecificID(symb_id), 0 }) != endo.end())
        continue;

      endo2eqs.emplace(symbol_table.getTypeSpecificID(symb_id), static_cast<int>(i));
      cout << "Endogenous " << symbol_table.getName(symb_id) << " normalized in equation " << i+1 << endl;
    }
  return endo2eqs;
}

pair<ModelTree::jacob_map_t, ModelTree::jacob_map_t>
ModelTree::evaluateAndReduceJacobian(const eval_context_t &eval_context, double cutoff, bool verbose)
{
  jacob_map_t contemporaneous_jacobian, static_jacobian;
  int nb_elements_contemparenous_Jacobian = 0;
  set<vector<int>> jacobian_elements_to_delete;
  for (const auto &[indices, d1] : derivatives[1])
    {
      int deriv_id = indices[1];
      if (getTypeByDerivID(deriv_id) == SymbolType::endogenous)
        {
          int eq = indices[0];
          int symb = getSymbIDByDerivID(deriv_id);
          int var = symbol_table.getTypeSpecificID(symb);
          int lag = getLagByDerivID(deriv_id);
          double val = 0;
          try
            {
              val = d1->eval(eval_context);
            }
          catch (ExprNode::EvalExternalFunctionException &e)
            {
              val = 1;
            }
          catch (ExprNode::EvalException &e)
            {
              cerr << "ERROR: evaluation of Jacobian failed for equation " << eq+1 << " (line " << equations_lineno[eq] << ") and variable " << symbol_table.getName(symb) << "(" << lag << ") [" << symb << "] !" << endl;
              d1->writeOutput(cerr, ExprNodeOutputType::matlabDynamicModelSparse, temporary_terms, {});
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
              dynamic_jacobian[{ lag, eq, var }] = d1;
            }
        }
    }

  // Get rid of the elements of the Jacobian matrix below the cutoff
  for (const auto &it : jacobian_elements_to_delete)
    derivatives[1].erase(it);

  if (jacobian_elements_to_delete.size() > 0)
    {
      cout << jacobian_elements_to_delete.size() << " elements among " << derivatives[1].size() << " in the incidence matrices are below the cutoff (" << cutoff << ") and are discarded" << endl
           << "The contemporaneous incidence matrix has " << nb_elements_contemparenous_Jacobian << " elements" << endl;
    }

  return { contemporaneous_jacobian, static_jacobian };
}

tuple<vector<pair<int, int>>, lag_lead_vector_t, lag_lead_vector_t,
      vector<unsigned int>, vector<unsigned int>, vector<unsigned int>, vector<unsigned int>>
ModelTree::select_non_linear_equations_and_variables(const vector<bool> &is_equation_linear)
{
  vector<int> eq2endo(equations.size(), 0);
  unsigned int num = 0;
  for (auto it : endo2eq)
    if (!is_equation_linear[it])
      num++;
  vector<int> endo2block(endo2eq.size(), 1);
  vector<pair<set<int>, pair<set<int>, vector<int>>>> components_set(num);
  int i = 0, j = 0;
  for (auto it : endo2eq)
    if (!is_equation_linear[it])
      {
        equation_reordered[i] = it;
        variable_reordered[i] = j;
        endo2block[j] = 0;
        components_set[endo2block[j]].first.insert(i);
        i++;
        j++;
      }
  auto [equation_lag_lead, variable_lag_lead] = getVariableLeadLagByBlock(endo2block, endo2block.size());
  vector<unsigned int> n_static(endo2eq.size(), 0), n_forward(endo2eq.size(), 0),
    n_backward(endo2eq.size(), 0), n_mixed(endo2eq.size(), 0);
  for (unsigned int i = 0; i < endo2eq.size(); i++)
    {
      if (variable_lag_lead[variable_reordered[i]].first != 0 && variable_lag_lead[variable_reordered[i]].second != 0)
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
  vector<pair<int, int>> blocks(1, {i, i});
  inv_equation_reordered.resize(nb_endo);
  inv_variable_reordered.resize(nb_endo);
  for (int i = 0; i < nb_endo; i++)
    {
      inv_variable_reordered[variable_reordered[i]] = i;
      inv_equation_reordered[equation_reordered[i]] = i;
    }
  return { blocks, equation_lag_lead, variable_lag_lead,
           n_static, n_forward, n_backward, n_mixed };
}

bool
ModelTree::computeNaturalNormalization()
{
  bool bool_result = true;
  set<pair<int, int>> result;
  endo2eq.resize(equations.size());
  for (int eq = 0; eq < static_cast<int>(equations.size()); eq++)
    if (!is_equation_linear[eq])
      {
        BinaryOpNode *eq_node = equations[eq];
        expr_t lhs = eq_node->arg1;
        result.clear();
        lhs->collectDynamicVariables(SymbolType::endogenous, result);
        if (result.size() == 1 && result.begin()->second == 0)
          {
            //check if the endogenous variable has not been already used in an other match !
            if (find(endo2eq.begin(), endo2eq.end(), result.begin()->first) == endo2eq.end())
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
ModelTree::computePrologueAndEpilogue(const jacob_map_t &static_jacobian_arg)
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
          for (const auto &it : endo)
            IM[i * n + endo2eq[it.first]] = true;
        }
    }
  else
    for (const auto &it : static_jacobian_arg)
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
      for (int i = prologue; i < n - static_cast<int>(epilogue); i++)
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

void
ModelTree::equationTypeDetermination(const map<tuple<int, int, int>, expr_t> &first_order_endo_derivatives, int mfs)
{
  expr_t lhs;
  BinaryOpNode *eq_node;
  EquationType Equation_Simulation_Type;
  equation_type_and_normalized_equation.clear();
  equation_type_and_normalized_equation.resize(equations.size());
  for (unsigned int i = 0; i < equations.size(); i++)
    {
      int eq = equation_reordered[i];
      int var = variable_reordered[i];
      eq_node = equations[eq];
      lhs = eq_node->arg1;
      Equation_Simulation_Type = EquationType::solve;
      pair<bool, expr_t> res;
      if (auto derivative = first_order_endo_derivatives.find({ eq, var, 0 });
          derivative != first_order_endo_derivatives.end())
        {
          set<pair<int, int>> result;
          derivative->second->collectEndogenous(result);
          auto d_endo_variable = result.find({ var, 0 });
          //Determine whether the equation could be evaluated rather than to be solved
          if (lhs->isVariableNodeEqualTo(SymbolType::endogenous, variable_reordered[i], 0) && derivative->second->isNumConstNodeEqualTo(1))
            Equation_Simulation_Type = EquationType::evaluate;
          else
            {
              vector<tuple<int, expr_t, expr_t>> List_of_Op_RHS;
              res = equations[eq]->normalizeEquation(var, List_of_Op_RHS);
              if (mfs == 2)
                {
                  if (d_endo_variable == result.end() && res.second)
                    Equation_Simulation_Type = EquationType::evaluate_s;
                }
              else if (mfs == 3)
                {
                  if (res.second) // The equation could be solved analytically
                    Equation_Simulation_Type = EquationType::evaluate_s;
                }
            }
        }
      equation_type_and_normalized_equation[eq] = { Equation_Simulation_Type, dynamic_cast<BinaryOpNode *>(res.second) };
    }
}

pair<lag_lead_vector_t, lag_lead_vector_t>
ModelTree::getVariableLeadLagByBlock(const vector<int> &components_set, int nb_blck_sim) const
{
  int nb_endo = symbol_table.endo_nbr();
  lag_lead_vector_t variable_lead_lag(nb_endo, { 0, 0 }), equation_lead_lag(nb_endo, { 0, 0 });
  vector<int> variable_blck(nb_endo), equation_blck(nb_endo);
  for (int i = 0; i < nb_endo; i++)
    {
      if (i < static_cast<int>(prologue))
        {
          variable_blck[variable_reordered[i]] = i;
          equation_blck[equation_reordered[i]] = i;
        }
      else if (i < static_cast<int>(components_set.size() + prologue))
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
  for (const auto &it : dynamic_jacobian)
    {
      auto [lag, j_1, i_1] = it.first;
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
  return { equation_lead_lag, variable_lead_lag };
}

tuple<vector<pair<int, int>>, lag_lead_vector_t, lag_lead_vector_t,
      vector<unsigned int>, vector<unsigned int>, vector<unsigned int>, vector<unsigned int>>
ModelTree::computeBlockDecompositionAndFeedbackVariablesForEachBlock(const jacob_map_t &static_jacobian, const equation_type_and_normalized_equation_t &Equation_Type, bool verbose_, bool select_feedback_variable)
{
  int nb_var = variable_reordered.size();
  int n = nb_var - prologue - epilogue;

  MFS::AdjacencyList_t G2(n);

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
          for (const auto &it : endo)
            tmp_normalized_contemporaneous_jacobian[{ i, it.first }] = 1;
        }
    }
  else
    tmp_normalized_contemporaneous_jacobian = static_jacobian;
  for (const auto &[key, value] : tmp_normalized_contemporaneous_jacobian)
    if (reverse_equation_reordered[key.first] >= static_cast<int>(prologue) && reverse_equation_reordered[key.first] < static_cast<int>(nb_var - epilogue)
        && reverse_variable_reordered[key.second] >= static_cast<int>(prologue) && reverse_variable_reordered[key.second] < static_cast<int>(nb_var - epilogue)
        && key.first != endo2eq[key.second])
      add_edge(vertex(reverse_equation_reordered[endo2eq[key.second]]-prologue, G2),
               vertex(reverse_equation_reordered[key.first]-prologue, G2),
               G2);

  vector<int> endo2block(num_vertices(G2)), discover_time(num_vertices(G2));
  boost::iterator_property_map<int *, boost::property_map<MFS::AdjacencyList_t, boost::vertex_index_t>::type, int, int &> endo2block_map(&endo2block[0], get(boost::vertex_index, G2));

  // Compute strongly connected components
  int num = strong_components(G2, endo2block_map);

  vector<pair<int, int>> blocks(num, { 0, 0 });

  // Create directed acyclic graph associated to the strongly connected components
  using DirectedGraph = boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS>;
  DirectedGraph dag(num);

  for (unsigned int i = 0; i < num_vertices(G2); i++)
    {
      MFS::AdjacencyList_t::out_edge_iterator it_out, out_end;
      MFS::AdjacencyList_t::vertex_descriptor vi = vertex(i, G2);
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

  auto [equation_lag_lead, variable_lag_lead] = getVariableLeadLagByBlock(endo2block, num);

  vector<int> tmp_equation_reordered(equation_reordered), tmp_variable_reordered(variable_reordered);
  int order = prologue;
  //Add a loop on vertices which could not be normalized or vertices related to lead variables => force those vertices to belong to the feedback set
  if (select_feedback_variable)
    {
      for (int i = 0; i < n; i++)
        if (Equation_Type[equation_reordered[i+prologue]].first == EquationType::solve
            || variable_lag_lead[variable_reordered[i+prologue]].second > 0
            || variable_lag_lead[variable_reordered[i+prologue]].first > 0
            || equation_lag_lead[equation_reordered[i+prologue]].second > 0
            || equation_lag_lead[equation_reordered[i+prologue]].first > 0
            || mfs == 0)
          add_edge(vertex(i, G2), vertex(i, G2), G2);
    }
  else
    for (int i = 0; i < n; i++)
      if (Equation_Type[equation_reordered[i+prologue]].first == EquationType::solve || mfs == 0)
        add_edge(vertex(i, G2), vertex(i, G2), G2);

  //Determines the dynamic structure of each equation
  vector<unsigned int> n_static(prologue+num+epilogue, 0), n_forward(prologue+num+epilogue, 0),
    n_backward(prologue+num+epilogue, 0), n_mixed(prologue+num+epilogue, 0);

  for (int i = 0; i < static_cast<int>(prologue); i++)
    if (variable_lag_lead[tmp_variable_reordered[i]].first != 0 && variable_lag_lead[tmp_variable_reordered[i]].second != 0)
      n_mixed[i]++;
    else if (variable_lag_lead[tmp_variable_reordered[i]].first == 0 && variable_lag_lead[tmp_variable_reordered[i]].second != 0)
      n_forward[i]++;
    else if (variable_lag_lead[tmp_variable_reordered[i]].first != 0 && variable_lag_lead[tmp_variable_reordered[i]].second == 0)
      n_backward[i]++;
    else if (variable_lag_lead[tmp_variable_reordered[i]].first == 0 && variable_lag_lead[tmp_variable_reordered[i]].second == 0)
      n_static[i]++;

  //For each block, the minimum set of feedback variable is computed
  // and the non-feedback variables are reordered to get
  // a sub-recursive block without feedback variables

  for (int i = 0; i < num; i++)
    {
      MFS::AdjacencyList_t G = MFS::extract_subgraph(G2, get<0>(components_set[i]));
      set<int> feed_back_vertices;
      MFS::AdjacencyList_t G1 = MFS::Minimal_set_of_feedback_vertex(feed_back_vertices, G);
      auto v_index = get(boost::vertex_index, G);
      get<1>(components_set[i]) = feed_back_vertices;
      blocks[i].second = feed_back_vertices.size();
      vector<int> Reordered_Vertice;
      MFS::Reorder_the_recursive_variables(G, feed_back_vertices, Reordered_Vertice);

      //First we have the recursive equations conditional on feedback variables
      for (int j = 0; j < 4; j++)
        for (int its : Reordered_Vertice)
          {
            bool something_done = false;
            if (j == 2 && variable_lag_lead[tmp_variable_reordered[its+prologue]].first != 0 && variable_lag_lead[tmp_variable_reordered[its+prologue]].second != 0)
              {
                n_mixed[prologue+i]++;
                something_done = true;
              }
            else if (j == 3 && variable_lag_lead[tmp_variable_reordered[its+prologue]].first == 0 && variable_lag_lead[tmp_variable_reordered[its+prologue]].second != 0)
              {
                n_forward[prologue+i]++;
                something_done = true;
              }
            else if (j == 1 && variable_lag_lead[tmp_variable_reordered[its+prologue]].first != 0 && variable_lag_lead[tmp_variable_reordered[its+prologue]].second == 0)
              {
                n_backward[prologue+i]++;
                something_done = true;
              }
            else if (j == 0 && variable_lag_lead[tmp_variable_reordered[its+prologue]].first == 0 && variable_lag_lead[tmp_variable_reordered[its+prologue]].second == 0)
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

      get<2>(components_set[i]) = Reordered_Vertice;
      //Second we have the equations related to the feedback variables
      for (int j = 0; j < 4; j++)
        for (int feed_back_vertice : feed_back_vertices)
          {
            bool something_done = false;
            if (j == 2 && variable_lag_lead[tmp_variable_reordered[v_index[vertex(feed_back_vertice, G)]+prologue]].first != 0 && variable_lag_lead[tmp_variable_reordered[v_index[vertex(feed_back_vertice, G)]+prologue]].second != 0)
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

  for (int i = 0; i < static_cast<int>(epilogue); i++)
    if (variable_lag_lead[tmp_variable_reordered[prologue+n+i]].first != 0 && variable_lag_lead[tmp_variable_reordered[prologue+n+i]].second != 0)
      n_mixed[prologue+num+i]++;
    else if (variable_lag_lead[tmp_variable_reordered[prologue+n+i]].first == 0 && variable_lag_lead[tmp_variable_reordered[prologue+n+i]].second != 0)
      n_forward[prologue+num+i]++;
    else if (variable_lag_lead[tmp_variable_reordered[prologue+n+i]].first != 0 && variable_lag_lead[tmp_variable_reordered[prologue+n+i]].second == 0)
      n_backward[prologue+num+i]++;
    else if (variable_lag_lead[tmp_variable_reordered[prologue+n+i]].first == 0 && variable_lag_lead[tmp_variable_reordered[prologue+n+i]].second == 0)
      n_static[prologue+num+i]++;

  inv_equation_reordered.resize(nb_var);
  inv_variable_reordered.resize(nb_var);
  for (int i = 0; i < nb_var; i++)
    {
      inv_variable_reordered[variable_reordered[i]] = i;
      inv_equation_reordered[equation_reordered[i]] = i;
    }

  return { blocks, equation_lag_lead, variable_lag_lead, n_static, n_forward, n_backward, n_mixed };
}

void
ModelTree::printBlockDecomposition(const vector<pair<int, int>> &blocks) const
{
  int largest_block = 0,
    Nb_SimulBlocks = 0,
    Nb_feedback_variable = 0;
  unsigned int Nb_TotalBlocks = getNbBlocks();
  for (unsigned int block = 0; block < Nb_TotalBlocks; block++)
    {
      BlockSimulationType simulation_type = getBlockSimulationType(block);
      if (simulation_type == BlockSimulationType::solveForwardComplete
          || simulation_type == BlockSimulationType::solveBackwardComplete
          || simulation_type == BlockSimulationType::solveTwoBoundariesComplete)
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

void
ModelTree::reduceBlocksAndTypeDetermination(const vector<pair<int, int>> &blocks, const equation_type_and_normalized_equation_t &Equation_Type, const vector<unsigned int> &n_static, const vector<unsigned int> &n_forward, const vector<unsigned int> &n_backward, const vector<unsigned int> &n_mixed, bool linear_decomposition)
{
  int i = 0;
  int count_equ = 0, blck_count_simult = 0;
  int Blck_Size, MFS_Size;
  int Lead, Lag;
  block_type_firstequation_size_mfs.clear();
  BlockSimulationType Simulation_Type, prev_Type = BlockSimulationType::unknown;
  int eq = 0;
  unsigned int l_n_static = 0, l_n_forward = 0, l_n_backward = 0, l_n_mixed = 0;
  for (i = 0; i < static_cast<int>(prologue+blocks.size()+epilogue); i++)
    {
      int first_count_equ = count_equ;
      if (i < static_cast<int>(prologue))
        {
          Blck_Size = 1;
          MFS_Size = 1;
        }
      else if (i < static_cast<int>(prologue+blocks.size()))
        {
          Blck_Size = blocks[blck_count_simult].first;
          MFS_Size = blocks[blck_count_simult].second;
          blck_count_simult++;
        }
      else if (i < static_cast<int>(prologue+blocks.size()+epilogue))
        {
          Blck_Size = 1;
          MFS_Size = 1;
        }

      Lag = Lead = 0;
      set<pair<int, int>> endo;
      for (count_equ = first_count_equ; count_equ < Blck_Size+first_count_equ; count_equ++)
        {
          endo.clear();
          equations[equation_reordered[count_equ]]->collectEndogenous(endo);
          for (const auto &it : endo)
            {
              int curr_variable = it.first;
              int curr_lag = it.second;
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
                  if (find(variable_reordered.begin()+first_count_equ, variable_reordered.begin()+(first_count_equ+Blck_Size), curr_variable)
                      != variable_reordered.begin()+(first_count_equ+Blck_Size)
                      && dynamic_jacobian.find({ curr_lag, equation_reordered[count_equ], curr_variable }) != dynamic_jacobian.end())
                    {
                      if (curr_lag > Lead)
                        Lead = curr_lag;
                      else if (-curr_lag > Lag)
                        Lag = -curr_lag;
                    }
                }
            }
        }
      if (Lag > 0 && Lead > 0)
        {
          if (Blck_Size == 1)
            Simulation_Type = BlockSimulationType::solveTwoBoundariesSimple;
          else
            Simulation_Type = BlockSimulationType::solveTwoBoundariesComplete;
        }
      else if (Blck_Size > 1)
        {
          if (Lead > 0)
            Simulation_Type = BlockSimulationType::solveBackwardComplete;
          else
            Simulation_Type = BlockSimulationType::solveForwardComplete;
        }
      else
        {
          if (Lead > 0)
            Simulation_Type = BlockSimulationType::solveBackwardSimple;
          else
            Simulation_Type = BlockSimulationType::solveForwardSimple;
        }
      l_n_static = n_static[i];
      l_n_forward = n_forward[i];
      l_n_backward = n_backward[i];
      l_n_mixed = n_mixed[i];
      if (Blck_Size == 1)
        {
          if (Equation_Type[equation_reordered[eq]].first == EquationType::evaluate
              || Equation_Type[equation_reordered[eq]].first == EquationType::evaluate_s)
            {
              if (Simulation_Type == BlockSimulationType::solveBackwardSimple)
                Simulation_Type = BlockSimulationType::evaluateBackward;
              else if (Simulation_Type == BlockSimulationType::solveForwardSimple)
                Simulation_Type = BlockSimulationType::evaluateForward;
            }
          if (i > 0)
            {
              bool is_lead = false, is_lag = false;
              int c_Size = get<2>(block_type_firstequation_size_mfs[block_type_firstequation_size_mfs.size()-1]);
              int first_equation = get<1>(block_type_firstequation_size_mfs[block_type_firstequation_size_mfs.size()-1]);
              if (c_Size > 0
                  && ((prev_Type == BlockSimulationType::evaluateForward && Simulation_Type == BlockSimulationType::evaluateForward && !is_lead)
                      || (prev_Type == BlockSimulationType::evaluateBackward && Simulation_Type == BlockSimulationType::evaluateBackward && !is_lag)))
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
              if ((prev_Type == BlockSimulationType::evaluateForward && Simulation_Type == BlockSimulationType::evaluateForward && !is_lead)
                  || (prev_Type == BlockSimulationType::evaluateBackward && Simulation_Type == BlockSimulationType::evaluateBackward && !is_lag))
                {
                  //merge the current block with the previous one
                  BlockSimulationType c_Type = get<0>(block_type_firstequation_size_mfs[block_type_firstequation_size_mfs.size()-1]);
                  c_Size++;
                  block_type_firstequation_size_mfs[block_type_firstequation_size_mfs.size()-1] = { c_Type, first_equation, c_Size, c_Size };
                  if (block_lag_lead[block_type_firstequation_size_mfs.size()-1].first > Lag)
                    Lag = block_lag_lead[block_type_firstequation_size_mfs.size()-1].first;
                  if (block_lag_lead[block_type_firstequation_size_mfs.size()-1].second > Lead)
                    Lead = block_lag_lead[block_type_firstequation_size_mfs.size()-1].second;
                  block_lag_lead[block_type_firstequation_size_mfs.size()-1] = { Lag, Lead };
                  auto tmp = block_col_type[block_col_type.size()-1];
                  block_col_type[block_col_type.size()-1] = { get<0>(tmp)+l_n_static, get<1>(tmp)+l_n_forward, get<2>(tmp)+l_n_backward, get<3>(tmp)+l_n_mixed };
                }
              else
                {
                  block_type_firstequation_size_mfs.emplace_back(Simulation_Type, eq, Blck_Size, MFS_Size);
                  block_lag_lead.emplace_back(Lag, Lead);
                  block_col_type.emplace_back(l_n_static, l_n_forward, l_n_backward, l_n_mixed);
                }
            }
          else
            {
              block_type_firstequation_size_mfs.emplace_back(Simulation_Type, eq, Blck_Size, MFS_Size);
              block_lag_lead.emplace_back(Lag, Lead);
              block_col_type.emplace_back(l_n_static, l_n_forward, l_n_backward, l_n_mixed);
            }
        }
      else
        {
          block_type_firstequation_size_mfs.emplace_back(Simulation_Type, eq, Blck_Size, MFS_Size);
          block_lag_lead.emplace_back(Lag, Lead);
          block_col_type.emplace_back(l_n_static, l_n_forward, l_n_backward, l_n_mixed);
        }
      prev_Type = Simulation_Type;
      eq += Blck_Size;
    }
}

vector<bool>
ModelTree::equationLinear(const map<tuple<int, int, int>, expr_t> &first_order_endo_derivatives) const
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

void
ModelTree::determineLinearBlocks()
{
  unsigned int nb_blocks = getNbBlocks();
  blocks_linear.clear();
  blocks_linear.resize(nb_blocks, true);
  for (unsigned int block = 0; block < nb_blocks; block++)
    {
      BlockSimulationType simulation_type = getBlockSimulationType(block);
      int block_size = getBlockSize(block);
      block_derivatives_equation_variable_laglead_nodeid_t derivatives_block = blocks_derivatives[block];
      int first_variable_position = getBlockFirstEquation(block);
      if (simulation_type == BlockSimulationType::solveBackwardComplete
          || simulation_type == BlockSimulationType::solveForwardComplete)
        for (const auto &[ignore, ignore2, lag, d1] : derivatives_block)
          {
            if (lag == 0)
              {
                set<pair<int, int>> endogenous;
                d1->collectEndogenous(endogenous);
                if (endogenous.size() > 0)
                  for (int l = 0; l < block_size; l++)
                    if (endogenous.find({ variable_reordered[first_variable_position+l], 0 }) != endogenous.end())
                      {
                        blocks_linear[block] = false;
                        goto the_end;
                      }
              }
          }
      else if (simulation_type == BlockSimulationType::solveTwoBoundariesComplete
               || simulation_type == BlockSimulationType::solveTwoBoundariesSimple)
        for (const auto &[ignore, ignore2, lag, d1] : derivatives_block)
          {
            set<pair<int, int>> endogenous;
            d1->collectEndogenous(endogenous);
            if (endogenous.size() > 0)
              for (int l = 0; l < block_size; l++)
                if (endogenous.find({ variable_reordered[first_variable_position+l], lag }) != endogenous.end())
                  {
                    blocks_linear[block] = false;
                    goto the_end;
                  }
          }
    the_end:
      ;
    }
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
  if (auto it = derivatives[1].find({ eq, getDerivID(symb_id, lag) });
      it != derivatives[1].end())
    it->second->writeOutput(output, output_type, temporary_terms, {});
  else
    output << 0;
}

void
ModelTree::computeDerivatives(int order, const set<int> &vars)
{
  assert(order >= 1);

  computed_derivs_order = order;

  // Do not shrink the vectors, since they have a minimal size of 4 (see constructor)
  derivatives.resize(max(static_cast<size_t>(order+1), derivatives.size()));
  NNZDerivatives.resize(max(static_cast<size_t>(order+1), NNZDerivatives.size()), 0);

  // First-order derivatives
  for (int var : vars)
    for (int eq = 0; eq < static_cast<int>(equations.size()); eq++)
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
          // We output symmetric elements at order = 2
          if (o == 2 && indices[1] != indices[2])
            NNZDerivatives[o] += 2;
          else
            NNZDerivatives[o]++;
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
  for (auto &equation : equations)
    equation->collectVariables(SymbolType::modelLocalVariable, used_local_vars);
  for (int used_local_var : used_local_vars)
    {
      VariableNode *v = AddVariable(used_local_var);
      temporary_terms_mlv[v] = local_variables_table.find(used_local_var)->second;
    }

  // Compute the temporary terms in equations and derivatives
  map<pair<int, int>, temporary_terms_t> temp_terms_map;
  map<expr_t, pair<int, pair<int, int>>> reference_count;

  for (auto &equation : equations)
    equation->computeTemporaryTerms({ 0, 0 },
                                    temp_terms_map,
                                    reference_count,
                                    is_matlab);

  for (int order = 1; order < static_cast<int>(derivatives.size()); order++)
    for (const auto &it : derivatives[order])
      it.second->computeTemporaryTerms({ 0, order },
                                       temp_terms_map,
                                       reference_count,
                                       is_matlab);

  /* If the user has specified the notmpterms option, clear all temporary
     terms, except those that correspond to external functions (since they are
     not optional) */
  if (no_tmp_terms)
    for (auto &it : temp_terms_map)
      // The following loop can be simplified with std::erase_if() in C++20
      for (auto it2 = it.second.begin(); it2 != it.second.end();)
        if (!dynamic_cast<AbstractExternalFunctionNode *>(*it2))
          it2 = it.second.erase(it2);
        else
          ++it2;

  // Fill the (now obsolete) temporary_terms structure
  temporary_terms.clear();
  for (const auto &it : temp_terms_map)
    temporary_terms.insert(it.second.begin(), it.second.end());

  // Fill the new structure
  temporary_terms_derivatives.clear();
  temporary_terms_derivatives.resize(derivatives.size());
  for (int order = 0; order < static_cast<int>(derivatives.size()); order++)
    temporary_terms_derivatives[order] = move(temp_terms_map[{ 0, order }]);

  // Compute indices in MATLAB/Julia vector
  int idx = 0;
  for (auto &it : temporary_terms_mlv)
    temporary_terms_idxs[it.first] = idx++;
  for (int order = 0; order < static_cast<int>(derivatives.size()); order++)
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
      if (dynamic_cast<AbstractExternalFunctionNode *>(it))
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
ModelTree::writeJsonTemporaryTerms(const temporary_terms_t &tt,
                                   temporary_terms_t &temp_term_union,
                                   ostream &output,
                                   deriv_node_temp_terms_t &tef_terms, const string &concat) const
{
  // Local var used to keep track of temp nodes already written
  bool wrote_term = false;
  temporary_terms_t tt2 = temp_term_union;

  output << R"("external_functions_temporary_terms_)" << concat << R"(": [)";
  for (auto it : tt)
    {
      if (dynamic_cast<AbstractExternalFunctionNode *>(it))
        {
          if (wrote_term)
            output << ", ";
          vector<string> efout;
          it->writeJsonExternalFunctionOutput(efout, tt2, tef_terms);
          for (auto it1 = efout.begin(); it1 != efout.end(); ++it1)
            {
              if (it1 != efout.begin())
                output << ", ";
              output << *it1;
            }
          wrote_term = true;
        }
      tt2.insert(it);
    }

  wrote_term = false;
  output << "]"
         << R"(, "temporary_terms_)" << concat << R"(": [)";
  for (const auto &it : tt)
    {
      if (wrote_term)
        output << ", ";
      output << R"({"temporary_term": ")";
      it->writeJsonOutput(output, tt, tef_terms);
      output << R"(")"
             << R"(, "value": ")";
      it->writeJsonOutput(output, temp_term_union, tef_terms);
      output << R"("})" << endl;
      wrote_term = true;

      temp_term_union.insert(it);
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
              cerr << "Warning: A .m file created by Dynare will have more than 32 nested parenthesis. MATLAB cannot support this. " << endl
                   << "         We are going to modify, albeit inefficiently, this output to have fewer than 32 nested parenthesis. " << endl
                   << "         It would hence behoove you to use the use_dll option of the model block to circumnavigate this problem." << endl
                   << "         If you have not yet set up a compiler on your system, see the MATLAB documentation for doing so." << endl
                   << "         For Windows, see: https://www.mathworks.com/help/matlab/matlab_external/install-mingw-support-package.html" << endl << endl;
              message_printed = true;
            }
          string str1 = str.substr(first_open_paren, matching_paren - first_open_paren + 1);
          string repstr, varname;
          while (testNestedParenthesis(str1))
            {
              size_t open_paren_idx = string::npos;
              size_t match_paren_idx = string::npos;
              size_t last_open_paren = string::npos;
              for (size_t j = 0; j < str1.length(); j++)
                {
                  if (str1.at(j) == '(')
                    {
                      // don't match, e.g. y(1)
                      if (size_t idx = str1.find_last_of("*/-+", j - 1);
                          j == 0 || (idx != string::npos && idx == j - 1))
                        open_paren_idx = j;
                      last_open_paren = j;
                    }
                  else if (str1.at(j) == ')')
                    {
                      // don't match, e.g. y(1)
                      if (size_t idx = str1.find_last_not_of("0123456789", j - 1);
                          idx != string::npos && idx != last_open_paren)
                        match_paren_idx = j;
                    }

                  if (open_paren_idx != string::npos && match_paren_idx != string::npos)
                    {
                      string val = str1.substr(open_paren_idx, match_paren_idx - open_paren_idx + 1);
                      if (auto it = tmp_paren_vars.find(val);
                          it == tmp_paren_vars.end())
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
          if (auto it = tmp_paren_vars.find(str1);
              it == tmp_paren_vars.end())
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
      if (dynamic_cast<AbstractExternalFunctionNode *>(it))
        {
          it->compileExternalFunctionOutput(code_file, instruction_number, false, tt2, map_idx, dynamic, steady_dynamic, tef_terms);
        }

      FNUMEXPR_ fnumexpr(TemporaryTerm, static_cast<int>(map_idx.find(it->idx)->second));
      fnumexpr.write(code_file, instruction_number);
      it->compile(code_file, instruction_number, false, tt2, map_idx, dynamic, steady_dynamic, tef_terms);
      if (dynamic)
        {
          FSTPT_ fstpt(static_cast<int>(map_idx.find(it->idx)->second));
          fstpt.write(code_file, instruction_number);
        }
      else
        {
          FSTPST_ fstpst(static_cast<int>(map_idx.find(it->idx)->second));
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

  output << R"("model_local_variables": [)";
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
        for (auto it1 = efout.begin(); it1 != efout.end(); ++it1)
          {
            if (it1 != efout.begin())
              output << ", ";
            output << *it1;
          }

        if (!efout.empty())
          output << ", ";

        /* We append underscores to avoid name clashes with "g1" or "oo_" (see
           also VariableNode::writeOutput) */
        output << R"({"variable": ")" << symbol_table.getName(id) << R"(__")"
               << R"(, "value": ")";
        value->writeJsonOutput(output, tt, tef_terms);
        output << R"("})" << endl;
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
  for (int eq = 0; eq < static_cast<int>(equations.size()); eq++)
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
  for (int eq = 0; eq < static_cast<int>(equations.size()); eq++)
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
      cerr << R"(Error : Can't open file ")" << filename << R"(" for writing)" << endl;
      exit(EXIT_FAILURE);
    }
  u_count_int = 0;
  for (const auto & [indices, d1] : derivatives[1])
    {
      int deriv_id = indices[1];
      if (getTypeByDerivID(deriv_id) == SymbolType::endogenous)
        {
          int eq = indices[0];
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
    u_count_int += symbol_table.endo_nbr();
  for (j = 0; j < static_cast<int>(symbol_table.endo_nbr()); j++)
    SaveCode.write(reinterpret_cast<char *>(&j), sizeof(j));
  for (j = 0; j < static_cast<int>(symbol_table.endo_nbr()); j++)
    SaveCode.write(reinterpret_cast<char *>(&j), sizeof(j));
  SaveCode.close();
}

void
ModelTree::writeLatexModelFile(const string &mod_basename, const string &latex_basename, ExprNodeOutputType output_type, bool write_equation_tags) const
{
  filesystem::create_directories(mod_basename + "/latex");

  ofstream output, content_output;
  string filename = mod_basename + "/latex/" + latex_basename + ".tex";
  string content_filename = mod_basename + "/latex/" + latex_basename + "_content" + ".tex";
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

  output << R"(\documentclass[10pt,a4paper]{article})" << endl
         << R"(\usepackage[landscape]{geometry})" << endl
         << R"(\usepackage{fullpage})" << endl
         << R"(\usepackage{amsfonts})" << endl
         << R"(\usepackage{breqn})" << endl
         << R"(\begin{document})" << endl
         << R"(\footnotesize)" << endl;

  // Write model local variables
  for (int id : local_variables_vector)
    {
      expr_t value = local_variables_table.find(id)->second;

      content_output << R"(\begin{dmath*})" << endl
                     << symbol_table.getTeXName(id) << " = ";
      // Use an empty set for the temporary terms
      value->writeOutput(content_output, output_type);
      content_output << endl << R"(\end{dmath*})" << endl;
    }

  for (int eq = 0; eq < static_cast<int>(equations.size()); eq++)
    {
      content_output << "% Equation " << eq + 1 << endl;
      if (write_equation_tags)
        equation_tags.writeLatexOutput(content_output, eq);

      content_output << R"(\begin{dmath})" << endl;
      // Here it is necessary to cast to superclass ExprNode, otherwise the overloaded writeOutput() method is not found
      dynamic_cast<ExprNode *>(equations[eq])->writeOutput(content_output, output_type);
      content_output << endl << R"(\end{dmath})" << endl;
    }

  output << R"(\include{)" << latex_basename + "_content" << "}" << endl
         << R"(\end{document})" << endl;

  output.close();
  content_output.close();
}

void
ModelTree::addEquation(expr_t eq, int lineno)
{
  auto beq = dynamic_cast<BinaryOpNode *>(eq);
  assert(beq && beq->op_code == BinaryOpcode::equal);

  equations.push_back(beq);
  equations_lineno.push_back(lineno);
}

vector<int>
ModelTree::includeExcludeEquations(set<pair<string, string>> &eqs, bool exclude_eqs,
                                   vector<BinaryOpNode *> &equations, vector<int> &equations_lineno,
                                   EquationTags &equation_tags, bool static_equations) const
{
  vector<int> excluded_vars;
  if (equations.empty())
    return excluded_vars;

  // Get equation numbers of tags
  set<int> tag_eqns;
  for (auto &it : eqs)
    if (auto tmp = equation_tags.getEqnsByTag(it.first, it.second); !tmp.empty())
      {
        tag_eqns.insert(tmp.begin(), tmp.end());
        eqs.erase(it);
      }

  if (tag_eqns.empty())
    return excluded_vars;

  set<int> eqns;
  if (exclude_eqs)
    eqns = tag_eqns;
  else
    for (size_t i = 0; i < equations.size(); i++)
      if (tag_eqns.find(i) == tag_eqns.end())
        eqns.insert(i);

  // remove from equations, equations_lineno, equation_tags
  vector<BinaryOpNode *> new_eqns;
  vector<int> new_equations_lineno;
  map<int, int> old_eqn_num_2_new;
  for (size_t i = 0; i < equations.size(); i++)
    if (eqns.find(i) != eqns.end())
      {
        if (auto tmp = equation_tags.getTagValueByEqnAndKey(i, "endogenous"); !tmp.empty())
          {
            excluded_vars.push_back(symbol_table.getID(tmp));
            set<pair<int, int>> result;
            equations[i]->arg1->collectDynamicVariables(SymbolType::endogenous, result);
            if (result.size() == 1)
              excluded_vars.push_back(result.begin()->first);
            else
              {
                cerr << "ERROR: Equation " << i
                     << " has been excluded but does not have a single variable on LHS or `endogenous` tag" << endl;
                exit(EXIT_FAILURE);
              }
          }
      }
    else
      {
        new_eqns.emplace_back(equations[i]);
        old_eqn_num_2_new[i] = new_eqns.size() - 1;
        new_equations_lineno.emplace_back(equations_lineno[i]);
      }
  int n_excl = equations.size() - new_eqns.size();

  equations = new_eqns;
  equations_lineno = new_equations_lineno;

  equation_tags.erase(eqns, old_eqn_num_2_new);

  if (!static_equations)
    for (size_t i = 0; i < excluded_vars.size(); i++)
      for (size_t j = i+1; j < excluded_vars.size(); j++)
        if (excluded_vars[i] == excluded_vars[j])
          {
            cerr << "Error: Variable " << symbol_table.getName(i) << " was excluded twice"
                 << " via in/exclude_eqs option" << endl;
            exit(EXIT_FAILURE);
          }

  cout << "Excluded " << n_excl << (static_equations ? " static " : " dynamic ")
       << "equation" << (n_excl > 1 ? "s" : "") << " via in/exclude_eqs option" << endl;

  return excluded_vars;
}

void
ModelTree::simplifyEquations()
{
  size_t last_subst_table_size = 0;
  map<VariableNode *, NumConstNode *> subst_table;
  // Equations with tags are excluded, in particular because of MCPs, see dynare#1697
  findConstantEquationsWithoutTags(subst_table);
  while (subst_table.size() != last_subst_table_size)
    {
      last_subst_table_size = subst_table.size();
      for (auto &equation : equations)
        equation = dynamic_cast<BinaryOpNode *>(equation->replaceVarsInEquation(subst_table));
      subst_table.clear();
      findConstantEquationsWithoutTags(subst_table);
    }
}

void
ModelTree::findConstantEquationsWithoutTags(map<VariableNode *, NumConstNode *> &subst_table) const
{
  for (size_t i = 0; i < equations.size(); i++)
    if (getEquationTags(i).empty())
      equations[i]->findConstantEquations(subst_table);
}

void
ModelTree::addEquation(expr_t eq, int lineno, const map<string, string> &eq_tags)
{
  equation_tags.add(equations.size(), eq_tags);
  addEquation(eq, lineno);
}

void
ModelTree::addAuxEquation(expr_t eq)
{
  auto beq = dynamic_cast<BinaryOpNode *>(eq);
  assert(beq && beq->op_code == BinaryOpcode::equal);

  aux_equations.push_back(beq);
}

void
ModelTree::addTrendVariables(const vector<int> &trend_vars, expr_t growth_factor) noexcept(false)
{
  for (int id : trend_vars)
    if (trend_symbols_map.find(id) != trend_symbols_map.end())
      throw TrendException(symbol_table.getName(id));
    else
      trend_symbols_map[id] = growth_factor;
}

void
ModelTree::addNonstationaryVariables(const vector<int> &nonstationary_vars, bool log_deflator, expr_t deflator) noexcept(false)
{
  for (int id : nonstationary_vars)
    if (nonstationary_symbols_map.find(id) != nonstationary_symbols_map.end())
      throw TrendException(symbol_table.getName(id));
    else
      nonstationary_symbols_map[id] = { log_deflator, deflator };
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
      for (int eq = 0; eq < static_cast<int>(equations.size()); eq++)
        {
          expr_t d = equations[eq]->getDerivative(param);
          if (d == Zero)
            continue;
          params_derivatives[{ 0, 1 }][{ eq, param }] = d;
        }

      for (int endoOrd = 1; endoOrd < static_cast<int>(derivatives.size()); endoOrd++)
        for (const auto &[indices, dprev] : derivatives[endoOrd])
          {
            expr_t d = dprev->getDerivative(param);
            if (d == Zero)
              continue;
            vector<int> new_indices = indices;
            new_indices.push_back(param);
            params_derivatives[{ endoOrd, 1 }][new_indices] = d;
          }
    }

  // Higher-order derivatives w.r.t. parameters
  for (int endoOrd = 0; endoOrd < static_cast<int>(derivatives.size()); endoOrd++)
    for (int paramOrd = 2; paramOrd <= paramsDerivsOrder; paramOrd++)
      for (const auto &[indices, dprev] : params_derivatives[{ endoOrd, paramOrd-1 }])
        for (int param : deriv_id_set)
          {
            if (indices.back() > param)
              continue;

            expr_t d = dprev->getDerivative(param);
            if (d == Zero)
              continue;
            vector<int> new_indices = indices;
            new_indices.push_back(param);
            // At this point, indices of both endogenous and parameters are sorted in non-decreasing order
            params_derivatives[{ endoOrd, paramOrd }][new_indices] = d;
          }
}

void
ModelTree::computeParamsDerivativesTemporaryTerms()
{
  map<expr_t, pair<int, pair<int, int>>> reference_count;

  /* The temp terms should be constructed in the same order as the for loops in
     {Static,Dynamic}Model::write{Json,}ParamsDerivativesFile() */
  params_derivs_temporary_terms.clear();
  for (const auto &[order, derivs] : params_derivatives)
    for (const auto &[indices, d] : derivs)
      d->computeTemporaryTerms(order, params_derivs_temporary_terms,
                               reference_count, true);

  int idx = 0;
  for (auto &[mlv, value] : temporary_terms_mlv)
    params_derivs_temporary_terms_idxs[mlv] = idx++;
  for (const auto &[order, tts] : params_derivs_temporary_terms)
    for (const auto &tt : tts)
      params_derivs_temporary_terms_idxs[tt] = idx++;
}

bool
ModelTree::isNonstationary(int symb_id) const
{
  return nonstationary_symbols_map.find(symb_id) != nonstationary_symbols_map.end();
}

void
ModelTree::writeJsonModelEquations(ostream &output, bool residuals) const
{
  if (residuals)
    output << endl << R"("residuals":[)" << endl;
  else
    output << endl << R"("model":[)" << endl;
  for (int eq = 0; eq < static_cast<int>(equations.size()); eq++)
    {
      if (eq > 0)
        output << ", ";

      BinaryOpNode *eq_node = equations[eq];
      expr_t lhs = eq_node->arg1;
      expr_t rhs = eq_node->arg2;

      if (residuals)
        {
          output << R"({"residual": {)"
                 << R"("lhs": ")";
          lhs->writeJsonOutput(output, temporary_terms, {});
          output << R"(")";

          output << R"(, "rhs": ")";
          rhs->writeJsonOutput(output, temporary_terms, {});
          output << R"(")";
          try
            {
              // Test if the right hand side of the equation is empty.
              if (rhs->eval(eval_context_t()) != 0)
                {
                  output << R"(, "rhs": ")";
                  rhs->writeJsonOutput(output, temporary_terms, {});
                  output << R"(")";
                }
            }
          catch (ExprNode::EvalException &e)
            {
            }
          output << "}";
        }
      else
        {
          output << R"({"lhs": ")";
          lhs->writeJsonOutput(output, {}, {});
          output << R"(", "rhs": ")";
          rhs->writeJsonOutput(output, {}, {});
          output << R"(")"
                 << R"(, "line": )" << equations_lineno[eq];

          if (auto eqtags = getEquationTags(eq);
              !eqtags.empty())
            {
              output << R"(, "tags": {)";
              int i = 0;
              for (const auto &[name, value] : eqtags)
                {
                  if (i != 0)
                    output << ", ";
                  output << R"(")" << name << R"(": ")" << value << R"(")";
                  i++;
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
    {
      cerr << "32-bit MATLAB not supported on macOS" << endl;
      exit(EXIT_FAILURE);
    }
  else if (mexext == "mexmaci64")
    return "maci64";
  else
    {
      cerr << "ERROR: 'mexext' option to preprocessor incorrectly set, needed with 'use_dll'" << endl;
      exit(EXIT_FAILURE);
    }
}

void
ModelTree::compileDll(const string &basename, const string &static_or_dynamic, const string &mexext, const filesystem::path &matlabroot, const filesystem::path &dynareroot) const
{
  const string opt_flags = "-O3 -g0 --param ira-max-conflict-table-size=1 -fno-forward-propagate -fno-gcse -fno-dce -fno-dse -fno-tree-fre -fno-tree-pre -fno-tree-cselim -fno-tree-dse -fno-tree-dce -fno-tree-pta -fno-gcse-after-reload";

  filesystem::path compiler;
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
      libs = "-lmex -lmx";
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
          auto mingwpath = dynareroot / (string{"mingw"} + (mexext == "mexw32" ? "32" : "64")) / "bin";
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
#ifdef __APPLE__
          char dynare_m_path[PATH_MAX];
          uint32_t size = PATH_MAX;
          string gcc_relative_path = "";
          if (_NSGetExecutablePath(dynare_m_path, &size) == 0)
            {
              string str = dynare_m_path;
              gcc_relative_path = str.substr(0, str.find_last_of("/")) + "/../../.brew/bin/gcc-9";
            }

          if (filesystem::exists(gcc_relative_path))
            compiler = gcc_relative_path;
          else if (filesystem::exists("/usr/local/bin/gcc-9"))
            compiler = "/usr/local/bin/gcc-9";
          else
            {
              cerr << "ERROR: You must install gcc-9 on your system before using the `use_dll` option of Dynare. "
                   << "You can do this via the Dynare installation package." << endl;
              exit(EXIT_FAILURE);
            }
#endif
          flags << " -fno-common -arch x86_64 -mmacosx-version-min=10.9 -Wl,-twolevel_namespace -undefined error -bundle";
          libs += " -lm -lstdc++";
        }
    }

  auto model_dir = filesystem::path{basename} / "model" / "src";
  filesystem::path main_src{model_dir / (static_or_dynamic + ".c")},
    mex_src{model_dir / (static_or_dynamic + "_mex.c")};

  filesystem::path mex_dir{"+" + basename};
  filesystem::path binary{mex_dir / (static_or_dynamic + "." + mexext)};

  ostringstream cmd;

#ifdef _WIN32
  /* On Windows, system() hands the command over to "cmd.exe /C". We need to
     enclose the whole command line within double quotes if we want the inner
     quotes to be correctly handled. See "cmd /?" for more details. */
  cmd << '"';
#endif

  if (user_set_compiler.empty())
    cmd << compiler << " ";
  else
    if (!filesystem::exists(user_set_compiler))
      {
        cerr << "Error: The specified compiler '" << user_set_compiler << "' cannot be found on your system" << endl;
        exit(EXIT_FAILURE);
      }
    else
      cmd << user_set_compiler << " ";

  if (user_set_subst_flags.empty())
    cmd << opt_flags << " " << flags.str() << " ";
  else
    cmd << user_set_subst_flags << " ";

  if (!user_set_add_flags.empty())
    cmd << user_set_add_flags << " ";

  cmd << main_src << " " << mex_src << " -o " << binary << " ";

  if (user_set_subst_libs.empty())
    cmd << libs;
  else
    cmd << user_set_subst_libs;

  if (!user_set_add_libs.empty())
    cmd << " " << user_set_add_libs;

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

void
ModelTree::reorderAuxiliaryEquations()
{
  using namespace boost;

  // Create the mapping between auxiliary variables and auxiliary equations
  int n = static_cast<int>(aux_equations.size());
  map<int, int> auxEndoToEq;
  for (int i = 0; i < n; i++)
    {
      auto varexpr = dynamic_cast<VariableNode *>(aux_equations[i]->arg1);
      assert(varexpr && symbol_table.getType(varexpr->symb_id) == SymbolType::endogenous);
      auxEndoToEq[varexpr->symb_id] = i;
    }
  assert(static_cast<int>(auxEndoToEq.size()) == n);

  /* Construct the directed acyclic graph where auxiliary equations are
     vertices and edges represent dependency relationships. */
  using Graph = adjacency_list<vecS, vecS, directedS>;
  Graph g(n);
  for (int i = 0; i < n; i++)
    {
      set<int> endos;
      aux_equations[i]->collectVariables(SymbolType::endogenous, endos);
      for (int endo : endos)
        if (auto it = auxEndoToEq.find(endo);
            it != auxEndoToEq.end() && it->second != i)
          add_edge(i, it->second, g);
    }

  // Topological sort of the graph
  using Vertex = graph_traits<Graph>::vertex_descriptor;
  vector<Vertex> ordered;
  topological_sort(g, back_inserter(ordered));

  // Reorder auxiliary equations accordingly
  auto aux_equations_old = aux_equations;
  auto index = get(vertex_index, g); // Maps vertex descriptors to their index
  for (int i = 0; i < n; i++)
    aux_equations[i] = aux_equations_old[index[ordered[i]]];
}
