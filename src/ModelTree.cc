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
#include "VariableDependencyGraph.hh"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/max_cardinality_matching.hpp>
#include <boost/graph/topological_sort.hpp>
#pragma GCC diagnostic pop

#ifdef __APPLE__
# include <mach-o/dyld.h>
#endif

#include <regex>
#include <utility>

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
  eq_idx_block2orig{m.eq_idx_block2orig},
  endo_idx_block2orig{m.endo_idx_block2orig},
  eq_idx_orig2block{m.eq_idx_orig2block},
  endo_idx_orig2block{m.endo_idx_orig2block},
  map_idx{m.map_idx},
  endo_max_leadlag_block{m.endo_max_leadlag_block},
  other_endo_max_leadlag_block{m.other_endo_max_leadlag_block},
  exo_max_leadlag_block{m.exo_max_leadlag_block},
  exo_det_max_leadlag_block{m.exo_det_max_leadlag_block},
  max_leadlag_block{m.max_leadlag_block},
  blocks{m.blocks},
  is_equation_linear{m.is_equation_linear},
  endo2eq{m.endo2eq},
  epilogue{m.epilogue},
  prologue{m.prologue},
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
  eq_idx_block2orig = m.eq_idx_block2orig;
  endo_idx_block2orig = m.endo_idx_block2orig;
  eq_idx_orig2block = m.eq_idx_orig2block;
  endo_idx_orig2block = m.endo_idx_orig2block;
  first_chain_rule_derivatives.clear();
  map_idx = m.map_idx;
  equation_type_and_normalized_equation.clear();
  blocks_derivatives.clear();
  derivative_endo.clear();
  derivative_other_endo.clear();
  derivative_exo.clear();
  derivative_exo_det.clear();
  endo_max_leadlag_block = m.endo_max_leadlag_block;
  other_endo_max_leadlag_block = m.other_endo_max_leadlag_block;
  exo_max_leadlag_block = m.exo_max_leadlag_block;
  exo_det_max_leadlag_block = m.exo_det_max_leadlag_block;
  max_leadlag_block = m.max_leadlag_block;
  blocks = m.blocks;
  is_equation_linear = m.is_equation_linear;
  endo2eq = m.endo2eq;
  epilogue = m.epilogue;
  prologue = m.prologue;
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
  for (const auto &[eq_and_endo, val] : contemporaneous_jacobian)
    add_edge(eq_and_endo.first + n, eq_and_endo.second, g);

  // Compute maximum cardinality matching
  vector<int> mate_map(2*n);

  bool check = checked_edmonds_maximum_cardinality_matching(g, &mate_map[0]);

  assert(check);

#ifdef DEBUG
  for (int i = 0; i < n; i++)
    cout << "Endogenous " << symbol_table.getName(symbol_table.getID(eEndogenous, i))
         << " matched with equation " << (mate_map[i]-n+1) << endl;
#endif

  // Create the resulting map, by copying the n first elements of mate_map, and substracting n to them
  endo2eq.resize(equations.size());
  transform(mate_map.begin(), mate_map.begin() + n, endo2eq.begin(), [=](int i) { return i-n; });

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
  cout << "Normalizing the model..." << endl;

  int n = equations.size();

  // Compute the maximum value of each row of the contemporaneous Jacobian matrix
  vector<double> max_val(n, 0.0);
  for (const auto &[eq_and_endo, val] : contemporaneous_jacobian)
    max_val[eq_and_endo.first] = max(max_val[eq_and_endo.first], fabs(val));

  // Compute normalized contemporaneous Jacobian
  jacob_map_t normalized_contemporaneous_jacobian(contemporaneous_jacobian);
  for (auto &[eq_and_endo, val] : normalized_contemporaneous_jacobian)
    val /= max_val[eq_and_endo.first];

  // We start with the highest value of the cutoff and try to normalize the model
  double current_cutoff = 0.99999999;
  const double cutoff_lower_limit = 1e-19;

  bool found_normalization = false;
  int last_suppressed = 0;
  while (!found_normalization && current_cutoff > cutoff_lower_limit)
    {
      // Drop elements below cutoff from normalized contemporaneous Jacobian
      jacob_map_t normalized_contemporaneous_jacobian_above_cutoff;
      int suppressed = 0;
      for (const auto &[eq_and_endo, val] : normalized_contemporaneous_jacobian)
        if (fabs(val) > max(current_cutoff, cutoff))
          normalized_contemporaneous_jacobian_above_cutoff[eq_and_endo] = val;
        else
          suppressed++;

      if (suppressed != last_suppressed)
        found_normalization = computeNormalization(normalized_contemporaneous_jacobian_above_cutoff, false);
      last_suppressed = suppressed;
      if (!found_normalization)
        {
          current_cutoff /= 2;
          // In this last case try to normalize with the complete jacobian
          if (current_cutoff <= cutoff_lower_limit)
            found_normalization = computeNormalization(normalized_contemporaneous_jacobian, false);
        }
    }

  if (!found_normalization)
    {
      cout << "Normalization failed with cutoff, trying symbolic normalization..." << endl;
      /* If no non-singular normalization can be found, try to find a
         normalization even with a potential singularity.
         TODO: Explain why symbolic_jacobian is not contemporaneous. */
      auto symbolic_jacobian = computeSymbolicJacobian();
      found_normalization = computeNormalization(symbolic_jacobian, true);
      if (found_normalization)
        {
          /* Update the Jacobian matrices by ensuring that they have a
             numerical value associated to each symbolic occurrence.
             TODO: Explain why this is needed. Maybe in order to have the
             right incidence matrix in computePrologueAndEpilogue? */
          for (const auto &[eq_and_endo, ignore] : symbolic_jacobian)
            {
              if (static_jacobian.find(eq_and_endo) == static_jacobian.end())
                static_jacobian[eq_and_endo] = 0;
              if (dynamic_jacobian.find({ 0, eq_and_endo.first, eq_and_endo.second }) == dynamic_jacobian.end())
                dynamic_jacobian[{ 0, eq_and_endo.first, eq_and_endo.second }] = Zero;
              if (contemporaneous_jacobian.find(eq_and_endo) == contemporaneous_jacobian.end())
                contemporaneous_jacobian[eq_and_endo] = 0;
              try
                {
                  if (derivatives[1].find({ eq_and_endo.first, getDerivID(symbol_table.getID(SymbolType::endogenous, eq_and_endo.second), 0) }) == derivatives[1].end())
                    derivatives[1][{ eq_and_endo.first, getDerivID(symbol_table.getID(SymbolType::endogenous, eq_and_endo.second), 0) }] = Zero;
                }
              catch (DataTree::UnknownDerivIDException &e)
                {
                  cerr << "The variable " << symbol_table.getName(symbol_table.getID(SymbolType::endogenous, eq_and_endo.second))
                       << " does not appear at the current period (i.e. with no lead and no lag); this case is not handled by the 'block' option of the 'model' block." << endl;
                  exit(EXIT_FAILURE);
                }
            }
        }
    }

  if (!found_normalization)
    {
      cerr << "No normalization could be computed. Aborting." << endl;
      exit(EXIT_FAILURE);
    }
}

pair<ModelTree::jacob_map_t, ModelTree::jacob_map_t>
ModelTree::evaluateAndReduceJacobian(const eval_context_t &eval_context, double cutoff, bool verbose)
{
  jacob_map_t contemporaneous_jacobian, static_jacobian;
  int nb_elements_contemporaneous_jacobian = 0;
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
                cout << "The coefficient related to variable " << var << " with lag " << lag << " in equation " << eq << " is equal to " << val << " and is set to 0 in the incidence matrix (size=" << symbol_table.endo_nbr() << ")." << endl;
              jacobian_elements_to_delete.insert({ eq, deriv_id });
            }
          else
            {
              if (lag == 0)
                {
                  nb_elements_contemporaneous_jacobian++;
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
  // TODO: Why?
  for (const auto &it : jacobian_elements_to_delete)
    derivatives[1].erase(it);

  if (jacobian_elements_to_delete.size() > 0)
    {
      cout << jacobian_elements_to_delete.size() << " elements among " << derivatives[1].size() << " in the incidence matrices are below the cutoff (" << cutoff << ") and are discarded." << endl
           << "The contemporaneous incidence matrix has " << nb_elements_contemporaneous_jacobian << " elements." << endl;
    }

  return { contemporaneous_jacobian, static_jacobian };
}

void
ModelTree::select_non_linear_equations_and_variables()
{
  prologue = 0;
  epilogue = 0;

  vector<int> endo2block(endo2eq.size(), 1); // The 1 is a dummy value, distinct from 0
  int i = 0;
  for (int endo = 0; endo < static_cast<int>(endo2eq.size()); endo++)
    {
      int eq = endo2eq[endo];
      if (!is_equation_linear[eq])
        {
          eq_idx_block2orig[i] = eq;
          endo_idx_block2orig[i] = endo;
          endo2block[i] = 0;
          i++;
        }
    }
  updateReverseVariableEquationOrderings();

  blocks.clear();
  blocks.resize(1);
  blocks[0].size = i;
  blocks[0].mfs_size = i;

  auto [equation_lag_lead, variable_lag_lead] = getVariableLeadLagByBlock(endo2block);

  for (int i = 0; i < blocks[0].size; i++)
    {
      auto [max_lag, max_lead] = variable_lag_lead[endo_idx_block2orig[i]];
      if (max_lag != 0 && max_lead != 0)
        blocks[0].n_mixed++;
      else if (max_lag == 0 && max_lead != 0)
        blocks[0].n_forward++;
      else if (max_lag != 0 && max_lead == 0)
        blocks[0].n_backward++;
      else
        blocks[0].n_static++;
    }
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
ModelTree::computePrologueAndEpilogue(const jacob_map_t &static_jacobian)
{
  const int n = equations.size();

  /* Initialize “eq_idx_block2orig” and “endo_idx_block2orig” to the identity
     permutation. */
  eq_idx_block2orig.resize(n);
  endo_idx_block2orig.resize(n);
  for (int i = 0; i < n; i++)
    {
      eq_idx_block2orig[i] = i;
      endo_idx_block2orig[endo2eq[i]] = i;
    }

  /* Compute incidence matrix, equations in rows, variables in columns. Row
     (resp. column) indices are to be interpreted according to
     “eq_idx_block2orig” (resp. “endo_idx_block2orig”). Stored in row-major
     order. */
  vector<bool> IM(n*n, false);
  if (cutoff == 0)
    for (int i = 0; i < n; i++)
      {
        set<pair<int, int>> endos_and_lags;
        equations[i]->collectEndogenous(endos_and_lags);
        for (const auto &[endo, lag] : endos_and_lags)
          IM[i * n + endo2eq[endo]] = true;
      }
  else
    for (const auto &[eq_and_endo, val] : static_jacobian)
      IM[eq_and_endo.first * n + endo2eq[eq_and_endo.second]] = true;

  bool something_has_been_done;
  // Find the prologue equations and place first the AR(1) shock equations first
  prologue = 0;
  do
    {
      something_has_been_done = false;
      int new_prologue = prologue;
      for (int i = prologue; i < n; i++)
        {
          int nze = 0;
          int k = 0;
          for (int j = new_prologue; j < n; j++)
            if (IM[i * n + j])
              {
                nze++;
                k = j;
              }
          if (nze == 1)
            {
              // Swap equations indexed by “new_prologue” and i
              for (int j = 0; j < n; j++)
                swap(IM[new_prologue * n + j], IM[i * n + j]);
              swap(eq_idx_block2orig[new_prologue], eq_idx_block2orig[i]);

              // Swap variables indexed by “new_prologue” and k (in the matching)
              for (int j = 0; j < n; j++)
                swap(IM[j * n + new_prologue], IM[j * n + k]);
              swap(endo_idx_block2orig[new_prologue], endo_idx_block2orig[k]);

              new_prologue++;
              something_has_been_done = true;
            }
        }
      prologue = new_prologue;
    }
  while (something_has_been_done);
  
  // Find the epilogue equations
  epilogue = 0;
  do
    {
      something_has_been_done = false;
      int new_epilogue = epilogue;
      for (int i = prologue; i < n - epilogue; i++)
        {
          int nze = 0;
          int k = 0;
          for (int j = prologue; j < n - new_epilogue; j++)
            if (IM[j * n + i])
              {
                nze++;
                k = j;
              }
          if (nze == 1)
            {
              for (int j = 0; j < n; j++)
                swap(IM[(n - 1 - new_epilogue) * n + j], IM[k * n + j]);
              swap(eq_idx_block2orig[n - 1 - new_epilogue], eq_idx_block2orig[k]);

              for (int j = 0; j < n; j++)
                swap(IM[j * n + n - 1 - new_epilogue], IM[j * n + i]);
              swap(endo_idx_block2orig[n - 1 - new_epilogue], endo_idx_block2orig[i]);

              new_epilogue++;
              something_has_been_done = true;
            }
        }
      epilogue = new_epilogue;
    }
  while (something_has_been_done);

  updateReverseVariableEquationOrderings();
}

void
ModelTree::equationTypeDetermination(const map<tuple<int, int, int>, expr_t> &first_order_endo_derivatives, int mfs)
{
  equation_type_and_normalized_equation.clear();
  equation_type_and_normalized_equation.resize(equations.size());
  for (int i = 0; i < static_cast<int>(equations.size()); i++)
    {
      int eq = eq_idx_block2orig[i];
      int var = endo_idx_block2orig[i];
      expr_t lhs = equations[eq]->arg1;
      EquationType Equation_Simulation_Type = EquationType::solve;
      BinaryOpNode *normalized_eq = nullptr;
      if (auto it = first_order_endo_derivatives.find({ eq, var, 0 });
          it != first_order_endo_derivatives.end())
        {
          expr_t derivative = it->second;
          // Determine whether the equation can be evaluated rather than solved
          if (lhs->isVariableNodeEqualTo(SymbolType::endogenous, endo_idx_block2orig[i], 0)
              && derivative->isNumConstNodeEqualTo(1))
            Equation_Simulation_Type = EquationType::evaluate;
          else
            {
              set<pair<int, int>> result;
              derivative->collectEndogenous(result);
              bool variable_not_in_derivative = result.find({ var, 0 }) == result.end();

              try
                {
                  normalized_eq = equations[eq]->normalizeEquation(symbol_table.getID(SymbolType::endogenous, var), 0);
                  if ((mfs == 2 && variable_not_in_derivative) || mfs == 3)
                    Equation_Simulation_Type = EquationType::evaluate_s;
                }
              catch (ExprNode::NormalizationFailed &e)
                {
                }
            }
        }
      equation_type_and_normalized_equation[eq] = { Equation_Simulation_Type, normalized_eq };
    }
}

pair<lag_lead_vector_t, lag_lead_vector_t>
ModelTree::getVariableLeadLagByBlock(const vector<int> &endo2simblock) const
{
  int nb_endo = symbol_table.endo_nbr();

  auto belong_to_same_block = [&](int endo, int eq)
                              {
                                int endo2 = endo_idx_orig2block[endo];
                                int eq2 = eq_idx_orig2block[eq];
                                if (endo2 < prologue || endo2 >= nb_endo-epilogue
                                    || eq2 < prologue || eq2 >= nb_endo-epilogue)
                                  return endo2 == eq2;
                                else
                                  return endo2simblock[endo2-prologue] == endo2simblock[eq2-prologue];
                              };

  lag_lead_vector_t variable_lag_lead(nb_endo, { 0, 0 }), equation_lag_lead(nb_endo, { 0, 0 });
  for (const auto &[key, value] : dynamic_jacobian)
    {
      auto [lag, eq, endo] = key;
      if (belong_to_same_block(endo, eq))
        {
          variable_lag_lead[endo].first = max(variable_lag_lead[endo].first, -lag);
          variable_lag_lead[endo].second = max(variable_lag_lead[endo].second, lag);
          equation_lag_lead[eq].first = max(equation_lag_lead[eq].first, -lag);
          equation_lag_lead[eq].second = max(equation_lag_lead[eq].second, lag);
        }
    }
  return { equation_lag_lead, variable_lag_lead };
}

lag_lead_vector_t
ModelTree::computeBlockDecompositionAndFeedbackVariablesForEachBlock(const jacob_map_t &static_jacobian, bool verbose_)
{
  int nb_var = symbol_table.endo_nbr();
  int nb_simvars = nb_var - prologue - epilogue;

  /* Construct the graph representing the dependencies between all
     variables that do not belong to the prologue or the epilogue.

     For detecting dependencies between variables, use the static jacobian,
     except when the cutoff is zero, in which case use the symbolic adjacency
     matrix */
  VariableDependencyGraph G(nb_simvars);
  for (const auto &[key, value] : cutoff == 0 ? computeSymbolicJacobian() : static_jacobian)
    {
      auto [eq, endo] = key;
      if (eq_idx_orig2block[eq] >= prologue
          && eq_idx_orig2block[eq] < nb_var - epilogue
          && endo_idx_orig2block[endo] >= prologue
          && endo_idx_orig2block[endo] < nb_var - epilogue
          && eq != endo2eq[endo])
        add_edge(vertex(eq_idx_orig2block[endo2eq[endo]]-prologue, G),
                 vertex(eq_idx_orig2block[eq]-prologue, G), G);
    }

  /* Identify the simultaneous blocks. Each simultaneous block is given an
     index, starting from 0, in recursive order */
  auto [num_simblocks, endo2simblock] = G.sortedStronglyConnectedComponents();

  int num_blocks = prologue+num_simblocks+epilogue;

  blocks.clear();
  blocks.resize(num_blocks);

  // Initialize size and mfs_size for prologue and epilogue
  for (int i = 0; i < prologue; i++)
    {
      blocks[i].size = 1;
      blocks[i].mfs_size = 1;
    }
  for (int i = 0; i < epilogue; i++)
    {
      blocks[prologue+num_simblocks+i].size = 1;
      blocks[prologue+num_simblocks+i].mfs_size = 1;
    }

  // Compute size and list of equations for simultaneous blocks
  vector<set<int>> eqs_in_simblock(num_simblocks);
  for (int i = 0; i < static_cast<int>(endo2simblock.size()); i++)
    {
      blocks[prologue+endo2simblock[i]].size++;
      eqs_in_simblock[endo2simblock[i]].insert(i);
    }

  // Determine the dynamic structure of each block
  auto [equation_lag_lead, variable_lag_lead] = getVariableLeadLagByBlock(endo2simblock);

  for (int var = 0; var < nb_var; var++)
    {
      auto [max_lag, max_lead] = variable_lag_lead[endo_idx_block2orig[var]];
      int blk;
      if (var < prologue)
        blk = var;
      else if (var >= nb_var - epilogue)
        blk = var - nb_simvars + num_simblocks;
      else
        blk = endo2simblock[var-prologue] + prologue;

      if (max_lag != 0 && max_lead != 0)
        blocks[blk].n_mixed++;
      else if (max_lag == 0 && max_lead != 0)
        blocks[blk].n_forward++;
      else if (max_lag != 0 && max_lead == 0)
        blocks[blk].n_backward++;
      else
        blocks[blk].n_static++;
    }

  /* For each simultaneous block, the minimum set of feedback variable is computed.
     Then, the variables within the blocks are reordered so that recursive
     (non-feedback) appear first, to get a sub-recursive block without feedback variables.
     Within each of the two sub-blocks, variables are reordered depending
     on their dynamic status: static first, then backward, mixed and forward. */

  /* First, add a loop on vertices which could not be normalized or vertices
     related to lead/lag variables. This forces those vertices to belong to the
     feedback set */
  for (int i = 0; i < nb_simvars; i++)
    if (equation_type_and_normalized_equation[eq_idx_block2orig[i+prologue]].first == EquationType::solve
        || variable_lag_lead[endo_idx_block2orig[i+prologue]].first > 0
        || variable_lag_lead[endo_idx_block2orig[i+prologue]].second > 0
        || equation_lag_lead[eq_idx_block2orig[i+prologue]].first > 0
        || equation_lag_lead[eq_idx_block2orig[i+prologue]].second > 0
        || mfs == 0)
      add_edge(vertex(i, G), vertex(i, G), G);

  const vector<int> old_eq_idx_block2orig(eq_idx_block2orig), old_endo_idx_block2orig(endo_idx_block2orig);
  int ordidx = prologue;
  for (int i = 0; i < num_simblocks; i++)
    {
      auto subG = G.extractSubgraph(eqs_in_simblock[i]);
      auto feed_back_vertices = subG.minimalSetOfFeedbackVertices();
      blocks[prologue+i].mfs_size = feed_back_vertices.size();
      auto reordered_vertices = subG.reorderRecursiveVariables(feed_back_vertices);

      const vector<pair<int, int>> dynamic_order{ make_pair(0, 0), make_pair(1, 0),
                                                  make_pair(1, 1), make_pair(0, 1) };

      // First the recursive equations conditional on feedback variables
      for (auto max_lag_lead : dynamic_order)
        for (int its : reordered_vertices)
            if (variable_lag_lead[old_endo_idx_block2orig[its+prologue]] == max_lag_lead)
              {
                eq_idx_block2orig[ordidx] = old_eq_idx_block2orig[its+prologue];
                endo_idx_block2orig[ordidx] = old_endo_idx_block2orig[its+prologue];
                ordidx++;
              }

      // Then the equations related to the feedback variables
      auto v_index1 = get(boost::vertex_index1, subG);
      for (auto max_lag_lead : dynamic_order)
        for (int fbvertex : feed_back_vertices)
          if (int idx = v_index1[vertex(fbvertex, subG)];
              variable_lag_lead[old_endo_idx_block2orig[idx+prologue]] == max_lag_lead)
            {
              eq_idx_block2orig[ordidx] = old_eq_idx_block2orig[idx+prologue];
              endo_idx_block2orig[ordidx] = old_endo_idx_block2orig[idx+prologue];
              ordidx++;
            }
    }

  updateReverseVariableEquationOrderings();

  return variable_lag_lead;
}

void
ModelTree::printBlockDecomposition() const
{
  int largest_block = 0, Nb_SimulBlocks = 0, Nb_feedback_variable = 0;
  int Nb_TotalBlocks = blocks.size();
  for (int block = 0; block < Nb_TotalBlocks; block++)
    if (BlockSimulationType simulation_type = blocks[block].simulation_type;
        simulation_type == BlockSimulationType::solveForwardComplete
        || simulation_type == BlockSimulationType::solveBackwardComplete
        || simulation_type == BlockSimulationType::solveTwoBoundariesComplete)
      {
        Nb_SimulBlocks++;
        if (int size = blocks[block].size;
            size > largest_block)
          {
            largest_block = size;
            Nb_feedback_variable = blocks[block].mfs_size;
          }
      }

  int Nb_RecursBlocks = Nb_TotalBlocks - Nb_SimulBlocks;
  cout << Nb_TotalBlocks << " block(s) found:" << endl
       << "  " << Nb_RecursBlocks << " recursive block(s) and " << Nb_SimulBlocks << " simultaneous block(s)." << endl
       << "  the largest simultaneous block has " << largest_block << " equation(s)" << endl
       << "                                 and " << Nb_feedback_variable << " feedback variable(s)." << endl;
}

void
ModelTree::reduceBlocksAndTypeDetermination(bool linear_decomposition)
{
  for (int i = 0; i < static_cast<int>(blocks.size()); i++)
    {
      int first_eq = (i == 0) ? 0 : blocks[i-1].first_equation+blocks[i-1].size;

      /* Compute the maximum lead and lag across all endogenous that appear in
         this block and that belong to it */
      int max_lag = 0, max_lead = 0;
      for (int eq = first_eq; eq < first_eq+blocks[i].size; eq++)
        {
          set<pair<int, int>> endos_and_lags;
          equations[eq_idx_block2orig[eq]]->collectEndogenous(endos_and_lags);
          for (const auto &[endo, lag] : endos_and_lags)
            {
              if (linear_decomposition)
                {
                  if (dynamic_jacobian.find({ lag, eq_idx_block2orig[eq], endo })
                      != dynamic_jacobian.end())
                    {
                      max_lead = max(lag, max_lead);
                      max_lag = max(-lag, max_lag);
                    }
                }
              else
                {
                  if (find(endo_idx_block2orig.begin()+first_eq,
                           endo_idx_block2orig.begin()+first_eq+blocks[i].size, endo)
                      != endo_idx_block2orig.begin()+first_eq+blocks[i].size
                      && dynamic_jacobian.find({ lag, eq_idx_block2orig[eq], endo })
                      != dynamic_jacobian.end())
                    {
                      max_lead = max(lag, max_lead);
                      max_lag = max(-lag, max_lag);
                    }
                }
            }
        }

      // Determine the block type
      BlockSimulationType Simulation_Type;
      if (max_lag > 0 && max_lead > 0)
        {
          if (blocks[i].size == 1)
            Simulation_Type = BlockSimulationType::solveTwoBoundariesSimple;
          else
            Simulation_Type = BlockSimulationType::solveTwoBoundariesComplete;
        }
      else if (blocks[i].size > 1)
        {
          if (max_lead > 0)
            Simulation_Type = BlockSimulationType::solveBackwardComplete;
          else
            Simulation_Type = BlockSimulationType::solveForwardComplete;
        }
      else
        {
          if (max_lead > 0)
            Simulation_Type = BlockSimulationType::solveBackwardSimple;
          else
            Simulation_Type = BlockSimulationType::solveForwardSimple;
        }

      if (blocks[i].size == 1)
        {
          // Determine if the block can simply be evaluated
          if (equation_type_and_normalized_equation[eq_idx_block2orig[first_eq]].first == EquationType::evaluate
              || equation_type_and_normalized_equation[eq_idx_block2orig[first_eq]].first == EquationType::evaluate_s)
            {
              if (Simulation_Type == BlockSimulationType::solveBackwardSimple)
                Simulation_Type = BlockSimulationType::evaluateBackward;
              else if (Simulation_Type == BlockSimulationType::solveForwardSimple)
                Simulation_Type = BlockSimulationType::evaluateForward;
            }

          /* Try to merge this block with the previous one.
             This is only possible if the two blocks can simply be evaluated
             (in the same direction), and if the merge does not break the
             restrictions on leads/lags. */
          if (i > 0)
            {
              bool is_lead = false, is_lag = false;
              for (int j = blocks[i-1].first_equation;
                   j < blocks[i-1].first_equation+blocks[i-1].size; j++)
                {
                  if (dynamic_jacobian.find({ -1, eq_idx_block2orig[first_eq], endo_idx_block2orig[j] })
                      != dynamic_jacobian.end())
                    is_lag = true;
                  if (dynamic_jacobian.find({ +1, eq_idx_block2orig[first_eq], endo_idx_block2orig[j] })
                      != dynamic_jacobian.end())
                    is_lead = true;
                }

              BlockSimulationType prev_Type = blocks[i-1].simulation_type;
              if ((prev_Type == BlockSimulationType::evaluateForward
                   && Simulation_Type == BlockSimulationType::evaluateForward
                   && !is_lead)
                  || (prev_Type == BlockSimulationType::evaluateBackward
                      && Simulation_Type == BlockSimulationType::evaluateBackward
                      && !is_lag))
                {
                  // Merge the current block into the previous one
                  blocks[i-1].size++;
                  blocks[i-1].mfs_size = blocks[i-1].size;
                  blocks[i-1].max_lag = max(blocks[i-1].max_lag, max_lag);
                  blocks[i-1].max_lead = max(blocks[i-1].max_lead, max_lead);
                  blocks[i-1].n_static += blocks[i].n_static;
                  blocks[i-1].n_forward += blocks[i].n_forward;
                  blocks[i-1].n_backward += blocks[i].n_backward;
                  blocks[i-1].n_mixed += blocks[i].n_mixed;
                  blocks.erase(blocks.begin()+i);
                  i--;
                  continue;
                }
            }
        }

      blocks[i].simulation_type = Simulation_Type;
      blocks[i].first_equation = first_eq;
      blocks[i].max_lag = max_lag;
      blocks[i].max_lead = max_lead;
    }
}

void
ModelTree::equationLinear(const map<tuple<int, int, int>, expr_t> &first_order_endo_derivatives)
{
  is_equation_linear.clear();
  is_equation_linear.resize(symbol_table.endo_nbr(), true);
  for (const auto &[indices, expr] : first_order_endo_derivatives)
    {
      set<pair<int, int>> endogenous;
      expr->collectEndogenous(endogenous);
      if (endogenous.size() > 0)
        {
          int eq = get<0>(indices);
          is_equation_linear[eq] = false;
        }
    }
}

void
ModelTree::determineLinearBlocks()
{
  // Note that field “linear” in class BlockInfo defaults to true
  for (int block = 0; block < static_cast<int>(blocks.size()); block++)
    {
      BlockSimulationType simulation_type = blocks[block].simulation_type;
      int block_size = blocks[block].size;
      block_derivatives_equation_variable_laglead_nodeid_t derivatives_block = blocks_derivatives[block];
      int first_variable_position = blocks[block].first_equation;
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
                    if (endogenous.find({ endo_idx_block2orig[first_variable_position+l], 0 }) != endogenous.end())
                      {
                        blocks[block].linear = false;
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
                if (endogenous.find({ endo_idx_block2orig[first_variable_position+l], lag }) != endogenous.end())
                  {
                    blocks[block].linear = false;
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
  for (j = 0; j < symbol_table.endo_nbr(); j++)
    SaveCode.write(reinterpret_cast<char *>(&j), sizeof(j));
  for (j = 0; j < symbol_table.endo_nbr(); j++)
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
  // Equations with “mcp” tag are excluded, see dynare#1697
  findConstantEquationsWithoutMcpTag(subst_table);
  while (subst_table.size() != last_subst_table_size)
    {
      last_subst_table_size = subst_table.size();
      for (auto &equation : equations)
        equation = dynamic_cast<BinaryOpNode *>(equation->replaceVarsInEquation(subst_table));
      subst_table.clear();
      findConstantEquationsWithoutMcpTag(subst_table);
    }
}

void
ModelTree::findConstantEquationsWithoutMcpTag(map<VariableNode *, NumConstNode *> &subst_table) const
{
  for (size_t i = 0; i < equations.size(); i++)
    if (auto tags = getEquationTags(i);
        tags.find("mcp") == tags.end())
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
    eq_idx_block2orig.push_back(j);

  for (int j = 0; j < symbol_table.endo_nbr(); j++)
    endo_idx_block2orig.push_back(j);
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

map<tuple<int, int, int>, expr_t>
ModelTree::collectFirstOrderDerivativesEndogenous()
{
  map<tuple<int, int, int>, expr_t> endo_derivatives;
  for (auto &[indices, d1] : derivatives[1])
    if (getTypeByDerivID(indices[1]) == SymbolType::endogenous)
      {
        int eq = indices[0];
        int var = symbol_table.getTypeSpecificID(getSymbIDByDerivID(indices[1]));
        int lag = getLagByDerivID(indices[1]);
        endo_derivatives[{ eq, var, lag }] = d1;
      }
  return endo_derivatives;
}

ModelTree::jacob_map_t
ModelTree::computeSymbolicJacobian() const
{
  jacob_map_t symbolic_jacobian;
  for (int i = 0; i < static_cast<int>(equations.size()); i++)
    {
      set<pair<int, int>> endos_and_lags;
      equations[i]->collectEndogenous(endos_and_lags);
      for (const auto &[endo, lag] : endos_and_lags)
        symbolic_jacobian[{ i, endo }] = 1;
    }
  return symbolic_jacobian;
}

void
ModelTree::updateReverseVariableEquationOrderings()
{
  int n = equations.size();
  eq_idx_orig2block.resize(n);
  endo_idx_orig2block.resize(n);
  for (int i = 0; i < n; i++)
    {
      endo_idx_orig2block[endo_idx_block2orig[i]] = i;
      eq_idx_orig2block[eq_idx_block2orig[i]] = i;
    }
}
