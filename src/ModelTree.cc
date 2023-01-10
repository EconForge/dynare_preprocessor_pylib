/*
 * Copyright © 2003-2023 Dynare Team
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
#include <algorithm>

/* NB: The workers must be listed *after* all the other static variables
   related to MEX compilation, so that when the preprocessor exits, the workers
   are destroyed *before* those variables (since the former rely on the latter
   for their functioning). */
condition_variable_any ModelTree::mex_compilation_cv;
mutex ModelTree::mex_compilation_mut;
vector<tuple<filesystem::path, set<filesystem::path>, string>> ModelTree::mex_compilation_queue;
set<filesystem::path> ModelTree::mex_compilation_ongoing, ModelTree::mex_compilation_done, ModelTree::mex_compilation_failed;
vector<jthread> ModelTree::mex_compilation_workers;

void
ModelTree::copyHelper(const ModelTree &m)
{
  auto f = [this](expr_t e) { return e->clone(*this); };

  // Equations
  for (const auto &it : m.equations)
    equations.push_back(dynamic_cast<BinaryOpNode *>(f(it)));
  for (const auto &it : m.aux_equations)
    aux_equations.push_back(dynamic_cast<BinaryOpNode *>(f(it)));

  auto convert_deriv_map = [f](const map<vector<int>, expr_t> &dm)
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
    params_derivatives.emplace(it.first, convert_deriv_map(it.second));
  for (const auto &it : m.jacobian_sparse_column_major_order)
    jacobian_sparse_column_major_order.emplace(it.first, f(it.second));

  auto convert_temporary_terms_t = [f](const temporary_terms_t &tt)
                                   {
                                     temporary_terms_t tt2;
                                     for (const auto &it : tt)
                                       tt2.insert(f(it));
                                     return tt2;
                                   };

  // Temporary terms
  for (const auto &it : m.temporary_terms_derivatives)
    temporary_terms_derivatives.push_back(convert_temporary_terms_t(it));
  for (const auto &it : m.temporary_terms_idxs)
    temporary_terms_idxs.emplace(f(it.first), it.second);
  for (const auto &it : m.params_derivs_temporary_terms)
    params_derivs_temporary_terms.emplace(it.first, convert_temporary_terms_t(it.second));
  for (const auto &it : m.params_derivs_temporary_terms_idxs)
    params_derivs_temporary_terms_idxs.emplace(f(it.first), it.second);

  // Other stuff
  for (const auto &it : m.trend_symbols_map)
    trend_symbols_map.emplace(it.first, f(it.second));
  for (const auto &it : m.nonstationary_symbols_map)
    nonstationary_symbols_map.emplace(it.first, pair{it.second.first, f(it.second.second)});

  for (const auto &it : m.equation_type_and_normalized_equation)
    equation_type_and_normalized_equation.emplace_back(it.first, dynamic_cast<BinaryOpNode *>(f(it.second)));

  for (const auto &it : m.blocks_derivatives)
    {
      map<tuple<int, int, int>, expr_t> v;
      for (const auto &it2 : it)
        v.emplace(it2.first, f(it2.second));
      blocks_derivatives.push_back(v);
    }

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
  for (const auto &it : m.blocks_temporary_terms)
    blocks_temporary_terms.push_back(convert_vector_tt(it));
  for (const auto &it : m.blocks_temporary_terms_idxs)
    blocks_temporary_terms_idxs.emplace(f(it.first), it.second);

  for (const auto &it : m.blocks_jacobian_sparse_column_major_order)
    {
      map<pair<int, int>, expr_t, columnMajorOrderLess> v;
      for (const auto &it2 : it)
        v.emplace(it2.first, f(it2.second));
      blocks_jacobian_sparse_column_major_order.push_back(v);
    }
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
  // Ensure that elements accessed by writeParamsDerivativesFileHelper() exist
  for (const auto &ord : {pair{0, 1}, pair{1, 1}, pair{0, 2}, pair{1, 2}, pair{2, 1}, pair{3, 1}})
    params_derivatives.emplace(ord, decltype(params_derivatives)::mapped_type{});
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
  jacobian_sparse_colptr{m.jacobian_sparse_colptr},
  eq_idx_block2orig{m.eq_idx_block2orig},
  endo_idx_block2orig{m.endo_idx_block2orig},
  eq_idx_orig2block{m.eq_idx_orig2block},
  endo_idx_orig2block{m.endo_idx_orig2block},
  block_decomposed{m.block_decomposed},
  time_recursive_block_decomposition{m.time_recursive_block_decomposition},
  blocks{m.blocks},
  endo2block{m.endo2block},
  eq2block{m.eq2block},
  blocks_jacobian_sparse_colptr{m.blocks_jacobian_sparse_colptr},
  endo2eq{m.endo2eq},
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

  jacobian_sparse_column_major_order.clear();
  jacobian_sparse_colptr = m.jacobian_sparse_colptr;

  params_derivatives.clear();

  temporary_terms_derivatives.clear();
  params_derivs_temporary_terms.clear();
  params_derivs_temporary_terms_idxs.clear();

  trend_symbols_map.clear();
  nonstationary_symbols_map.clear();

  eq_idx_block2orig = m.eq_idx_block2orig;
  endo_idx_block2orig = m.endo_idx_block2orig;
  eq_idx_orig2block = m.eq_idx_orig2block;
  endo_idx_orig2block = m.endo_idx_orig2block;
  equation_type_and_normalized_equation.clear();
  blocks_derivatives.clear();
  block_decomposed = m.block_decomposed;
  time_recursive_block_decomposition = m.time_recursive_block_decomposition;
  blocks = m.blocks;
  endo2block = m.endo2block;
  eq2block = m.eq2block;
  blocks_temporary_terms.clear();
  blocks_temporary_terms_idxs.clear();
  blocks_jacobian_sparse_column_major_order.clear();
  blocks_jacobian_sparse_colptr = m.blocks_jacobian_sparse_colptr;
  endo2eq = m.endo2eq;
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

  // Create the resulting map, by copying the n first elements of mate_map, and substracting n to them
  endo2eq.resize(equations.size());
  transform(mate_map.begin(), mate_map.begin() + n, endo2eq.begin(), [=](int i) { return i-n; });

  // Check if all variables are normalized
  if (auto it = find(mate_map.begin(), mate_map.begin() + n, boost::graph_traits<BipartiteGraph>::null_vertex());
      it != mate_map.begin() + n)
    {
      if (verbose)
        cerr << "Could not normalize the " << modelClassName() << ". Variable "
             << symbol_table.getName(symbol_table.getID(SymbolType::endogenous, it - mate_map.begin()))
             << " is not in the maximum cardinality matching." << endl;
      check = false;
    }
  return check;
}

bool
ModelTree::computeNonSingularNormalization(const jacob_map_t &contemporaneous_jacobian)
{
  int n = equations.size();

  /* Optimal policy models (discretionary, or Ramsey before computing FOCs) do
     not have as many equations as variables. */
  if (n != symbol_table.endo_nbr())
    {
      cout << "The " << modelClassName() << " cannot be normalized, since it does not have as many equations as variables." << endl;
      return false;
    }

  cout << "Normalizing the " << modelClassName() << "..." << endl;

  // Compute the maximum value of each row of the contemporaneous Jacobian matrix
  vector max_val(n, 0.0);
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
         normalization even with a potential singularity. */
      auto symbolic_jacobian = computeSymbolicJacobian(true);
      found_normalization = computeNormalization(symbolic_jacobian, true);
    }

#ifdef DEBUG
  for (size_t i {0}; i < equations.size(); i++)
    cout << "Variable " << symbol_table.getName(symbol_table.getID(SymbolType::endogenous, i))
         << " associated to equation " << endo2eq[i] << " (" << equation_tags.getTagValueByEqnAndKey(endo2eq[i], "name") << ")" << endl;
#endif

  /* NB: If normalization failed, an explanatory message has been printed by the last call
     to computeNormalization(), which has verbose=true */
  return found_normalization;
}

ModelTree::jacob_map_t
ModelTree::evaluateAndReduceJacobian(const eval_context_t &eval_context) const
{
  jacob_map_t contemporaneous_jacobian;
  for (const auto &[indices, d1] : derivatives[1])
    {
      int deriv_id = indices[1];
      if (getTypeByDerivID(deriv_id) == SymbolType::endogenous)
        {
          int eq = indices[0];
          int var { getTypeSpecificIDByDerivID(deriv_id) };
          int lag = getLagByDerivID(deriv_id);
          double val { [&]
          {
            try
              {
                return d1->eval(eval_context);
              }
            catch (ExprNode::EvalExternalFunctionException &e)
              {
                return 1.0;
              }
            /* Other types of EvalException should not happen (all symbols should
               have a value; we don’t evaluate an equal sign) */
          }() };

          if ((isnan(val) || fabs(val) >= cutoff) && lag == 0)
            contemporaneous_jacobian[{ eq, var }] = val;
        }
    }

  return contemporaneous_jacobian;
}

pair<int, int>
ModelTree::computePrologueAndEpilogue()
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
  vector IM(n*n, false);
  for (int i = 0; i < n; i++)
    {
      set<pair<int, int>> endos_and_lags;
      equations[i]->collectEndogenous(endos_and_lags);
      for (auto [endo, lag] : endos_and_lags)
        if (!time_recursive_block_decomposition || lag == 0)
          IM[i * n + endo2eq[endo]] = true;
    }

  bool something_has_been_done;
  // Find the prologue equations and place first the AR(1) shock equations first
  int prologue = 0;
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
  int epilogue = 0;
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

  return { prologue, epilogue };
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
              bool variable_not_in_derivative = !result.contains({ var, 0 });

              try
                {
                  normalized_eq = equations[eq]->normalizeEquation(symbol_table.getID(SymbolType::endogenous, var), 0);
                  if ((mfs == 2 && variable_not_in_derivative) || mfs == 3)
                    Equation_Simulation_Type = EquationType::evaluateRenormalized;
                }
              catch (ExprNode::NormalizationFailed &e)
                {
                }
            }
        }
      equation_type_and_normalized_equation[eq] = { Equation_Simulation_Type, normalized_eq };
    }
}

void
ModelTree::computeDynamicStructureOfBlock(int blk)
{
  vector max_endo_lag_lead(blocks[blk].size, pair{0, 0});
  blocks[blk].max_endo_lag = blocks[blk].max_endo_lead = 0;
  blocks[blk].max_other_endo_lag = blocks[blk].max_other_endo_lead = 0;
  blocks[blk].max_exo_lag = blocks[blk].max_exo_lead = 0;
  blocks[blk].max_exo_det_lag = blocks[blk].max_exo_det_lead = 0;
  for (int eq = 0; eq < blocks[blk].size; eq++)
    {
      set<pair<int, int>> endos_and_lags;
      expr_t e = getBlockEquationExpr(blk, eq);

      /* Compute max lags/leads for endogenous. Also fill per-variable structure
         for endos belonging to this block */
      e->collectEndogenous(endos_and_lags);
      for (auto [endo, lag] : endos_and_lags)
        if (endo2block[endo] == blk)
          {
            blocks[blk].max_endo_lag = max(blocks[blk].max_endo_lag, -lag);
            blocks[blk].max_endo_lead = max(blocks[blk].max_endo_lead, lag);
            auto &[max_endo_lag, max_endo_lead] = max_endo_lag_lead[getBlockInitialVariableID(blk, endo)];
            max_endo_lag = max(max_endo_lag, -lag);
            max_endo_lead = max(max_endo_lead, lag);
          }
        else
          {
            blocks[blk].max_other_endo_lag = max(blocks[blk].max_other_endo_lag, -lag);
            blocks[blk].max_other_endo_lead = max(blocks[blk].max_other_endo_lead, lag);
          }

      // Compute max lags/leads for exogenous
      blocks[blk].max_exo_lag = max(e->maxExoLag(), blocks[blk].max_exo_lag);
      blocks[blk].max_exo_lead = max(e->maxExoLead(), blocks[blk].max_exo_lead);

      // Compute max lags/leads for deterministic exogenous
      set<pair<int, int>> dynvars;
      e->collectDynamicVariables(SymbolType::exogenousDet, dynvars);
      for (auto [symb_id, lag] : dynvars)
        {
          blocks[blk].max_exo_det_lag = max(-lag, blocks[blk].max_exo_det_lag);
          blocks[blk].max_exo_det_lead = max(lag, blocks[blk].max_exo_det_lead);
        }
    }

  // Compute max lags/leads over all variables
  blocks[blk].max_lag = max(blocks[blk].max_endo_lag, max(blocks[blk].max_other_endo_lag,
                                                          max(blocks[blk].max_exo_lag,
                                                              blocks[blk].max_exo_det_lag)));
  blocks[blk].max_lead = max(blocks[blk].max_endo_lead, max(blocks[blk].max_other_endo_lead,
                                                            max(blocks[blk].max_exo_lead,
                                                                blocks[blk].max_exo_det_lead)));

  // Categorize endos that belong to the block
  blocks[blk].n_mixed = blocks[blk].n_forward = blocks[blk].n_backward = blocks[blk].n_static = 0;
  for (int var = 0; var < blocks[blk].size; var++)
    {
      auto [max_lag, max_lead] = max_endo_lag_lead[var];
      if (max_lag != 0 && max_lead != 0)
        blocks[blk].n_mixed++;
      else if (max_lag == 0 && max_lead != 0)
        blocks[blk].n_forward++;
      else if (max_lag != 0 && max_lead == 0)
        blocks[blk].n_backward++;
      else
        blocks[blk].n_static++;
    }
}

void
ModelTree::computeSimulationTypeOfBlock(int blk)
{
  auto &type = blocks[blk].simulation_type;
  if (blocks[blk].max_endo_lag > 0 && blocks[blk].max_endo_lead > 0)
    {
      if (blocks[blk].size == 1)
        type = BlockSimulationType::solveTwoBoundariesSimple;
      else
        type = BlockSimulationType::solveTwoBoundariesComplete;
    }
  else if (blocks[blk].size > 1)
    {
      if (blocks[blk].max_endo_lead > 0)
        type = BlockSimulationType::solveBackwardComplete;
      else
        type = BlockSimulationType::solveForwardComplete;
    }
  else
    {
      bool can_eval = (getBlockEquationType(blk, 0) == EquationType::evaluate
                       || getBlockEquationType(blk, 0) == EquationType::evaluateRenormalized);
      if (blocks[blk].max_endo_lead > 0)
        type = can_eval ? BlockSimulationType::evaluateBackward :
          BlockSimulationType::solveBackwardSimple;
      else
        type = can_eval ? BlockSimulationType::evaluateForward :
          BlockSimulationType::solveForwardSimple;
    }
}

pair<lag_lead_vector_t, lag_lead_vector_t>
ModelTree::getVariableLeadLagByBlock() const
{
  int nb_endo = symbol_table.endo_nbr();

  lag_lead_vector_t variable_lag_lead(nb_endo, { 0, 0 }), equation_lag_lead(nb_endo, { 0, 0 });
  for (int eq = 0; eq < nb_endo; eq++)
    {
      set<pair<int, int>> endos_and_lags;
      equations[eq]->collectEndogenous(endos_and_lags);
      for (auto [endo, lag] : endos_and_lags)
        if (endo2block[endo] == eq2block[eq])
          {
            variable_lag_lead[endo].first = max(variable_lag_lead[endo].first, -lag);
            variable_lag_lead[endo].second = max(variable_lag_lead[endo].second, lag);
            equation_lag_lead[eq].first = max(equation_lag_lead[eq].first, -lag);
            equation_lag_lead[eq].second = max(equation_lag_lead[eq].second, lag);
          }
    }
  return { equation_lag_lead, variable_lag_lead };
}

void
ModelTree::computeBlockDecomposition(int prologue, int epilogue)
{
  int nb_var = symbol_table.endo_nbr();
  int nb_simvars = nb_var - prologue - epilogue;

  /* Construct the graph representing the dependencies between all
     variables that do not belong to the prologue or the epilogue.

     For detecting dependencies between variables, use the symbolic adjacency
     matrix */
  VariableDependencyGraph G(nb_simvars);
  for (const auto &[key, value] : computeSymbolicJacobian(time_recursive_block_decomposition))
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
  auto [num_simblocks, simvar2simblock] = G.sortedStronglyConnectedComponents();

  int num_blocks = prologue+num_simblocks+epilogue;

  blocks.clear();
  blocks.resize(num_blocks);
  endo2block.resize(nb_var);
  eq2block.resize(nb_var);

  // Initialize size and mfs_size for prologue and epilogue, plus eq/endo→block mappings
  for (int blk = 0; blk < num_blocks; blk++)
    if (blk < prologue || blk >= num_blocks-epilogue)
      {
        int var_eq = (blk < prologue ? blk : blk-num_simblocks+nb_simvars);
        blocks[blk].size = 1;
        blocks[blk].mfs_size = 1;
        blocks[blk].first_equation = var_eq;
        endo2block[endo_idx_block2orig[var_eq]] = blk;
        eq2block[eq_idx_block2orig[var_eq]] = blk;
      }

  // Initialize size for simultaneous blocks, plus eq/endo→block mappings
  vector<vector<int>> simblock2simvars(num_simblocks);
  for (int i = 0; i < static_cast<int>(simvar2simblock.size()); i++)
    {
      simblock2simvars[simvar2simblock[i]].push_back(i);
      int blk = prologue+simvar2simblock[i];
      blocks[blk].size++;
      endo2block[endo_idx_block2orig[prologue+i]] = blk;
      eq2block[eq_idx_block2orig[prologue+i]] = blk;
    }

  // Determine the dynamic structure of each block
  auto [equation_lag_lead, variable_lag_lead] = getVariableLeadLagByBlock();

  /* For each simultaneous block, the minimum set of feedback variables is
     computed. Then, the variables within the blocks are reordered so that
     recursive (non-feedback) appear first, in recursive order. They are
     followed by feedback variables, which are reordered according to their
     dynamic status (static first, then backward, mixed and forward). */

  /* First, add a loop on vertices which could not be normalized or vertices
     related to lead/lag variables. This forces those vertices to belong to the
     feedback set */
  for (int i = 0; i < nb_simvars; i++)
    if (equation_type_and_normalized_equation[eq_idx_block2orig[i+prologue]].first == EquationType::solve
        || (!time_recursive_block_decomposition &&
            (variable_lag_lead[endo_idx_block2orig[i+prologue]].first > 0
             || variable_lag_lead[endo_idx_block2orig[i+prologue]].second > 0
             || equation_lag_lead[eq_idx_block2orig[i+prologue]].first > 0
             || equation_lag_lead[eq_idx_block2orig[i+prologue]].second > 0))
        || mfs == 0)
      add_edge(vertex(i, G), vertex(i, G), G);

  const vector<int> old_eq_idx_block2orig(eq_idx_block2orig), old_endo_idx_block2orig(endo_idx_block2orig);
  int ordidx = prologue;
  for (int blk = prologue; blk < prologue+num_simblocks; blk++)
    {
      blocks[blk].first_equation = (blk == 0 ? 0 : blocks[blk-1].first_equation + blocks[blk-1].size);
      auto subG = G.extractSubgraph(simblock2simvars[blk-prologue]);
      auto feed_back_vertices = subG.minimalSetOfFeedbackVertices();
      blocks[blk].mfs_size = feed_back_vertices.size();
      auto recursive_vertices = subG.reorderRecursiveVariables(feed_back_vertices);
      auto v_index1 = get(boost::vertex_index1, subG);

      /* First the recursive variables conditional on feedback variables, in
         recursive order */
      for (int vtx : recursive_vertices)
        {
          int simvar { v_index1[vertex(vtx, subG)] };
          eq_idx_block2orig[ordidx] = old_eq_idx_block2orig[simvar+prologue];
          endo_idx_block2orig[ordidx] = old_endo_idx_block2orig[simvar+prologue];
          ordidx++;
        }

      // Then the feedback variables, reordered by dynamic status
      for (auto max_lag_lead : { pair{0, 0}, pair{1, 0}, pair{1, 1}, pair{0, 1} })
        for (int vtx : feed_back_vertices)
          if (int simvar = v_index1[vertex(vtx, subG)];
              variable_lag_lead[old_endo_idx_block2orig[simvar+prologue]] == max_lag_lead)
            {
              eq_idx_block2orig[ordidx] = old_eq_idx_block2orig[simvar+prologue];
              endo_idx_block2orig[ordidx] = old_endo_idx_block2orig[simvar+prologue];
              ordidx++;
            }
    }

  updateReverseVariableEquationOrderings();

  for (int blk = 0; blk < static_cast<int>(blocks.size()); blk++)
    {
      computeDynamicStructureOfBlock(blk);
      computeSimulationTypeOfBlock(blk);
    }
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
ModelTree::reduceBlockDecomposition()
{
  for (int blk = 1; blk < static_cast<int>(blocks.size()); blk++)
    if (blocks[blk].size == 1)
      {
        /* Try to merge this block with the previous one.
           This is only possible if the two blocks can simply be evaluated
           (in the same direction), and if the merge does not break the
           restrictions on leads/lags. */
        set<pair<int, int>> endos_and_lags;
        getBlockEquationExpr(blk, 0)->collectEndogenous(endos_and_lags);
        bool is_lead = false, is_lag = false;
        for (int var = 0; var < blocks[blk-1].size; var++)
          {
            is_lag = is_lag || endos_and_lags.contains({ getBlockVariableID(blk-1, var), -1 });
            is_lead = is_lead || endos_and_lags.contains({ getBlockVariableID(blk-1, var), 1 });
          }

        if ((blocks[blk-1].simulation_type == BlockSimulationType::evaluateForward
             && blocks[blk].simulation_type == BlockSimulationType::evaluateForward
             && !is_lead)
            || (blocks[blk-1].simulation_type == BlockSimulationType::evaluateBackward
                && blocks[blk].simulation_type == BlockSimulationType::evaluateBackward
                && !is_lag))
          {
            // Merge the current block into the previous one
            blocks[blk-1].size++;
            blocks[blk-1].mfs_size = blocks[blk-1].size;
            computeDynamicStructureOfBlock(blk-1);
            blocks.erase(blocks.begin()+blk);
            for (auto &b : endo2block)
              if (b >= blk)
                b--;
            for (auto &b : eq2block)
              if (b >= blk)
                b--;
            blk--;
            continue;
          }
      }
}

void
ModelTree::determineLinearBlocks()
{
  // Note that field “linear” in class BlockInfo defaults to true
  for (int blk = 0; blk < static_cast<int>(blocks.size()); blk++)
    switch (blocks[blk].simulation_type)
      {
      case BlockSimulationType::solveBackwardSimple:
      case BlockSimulationType::solveBackwardComplete:
      case BlockSimulationType::solveForwardSimple:
      case BlockSimulationType::solveForwardComplete:
        for (const auto &[indices, d1] : blocks_derivatives[blk])
          {
            int lag = get<2>(indices);
            if (lag == 0)
              {
                set<pair<int, int>> endogenous;
                d1->collectEndogenous(endogenous);
                for (int l = 0; l < blocks[blk].size; l++)
                  if (endogenous.contains({ endo_idx_block2orig[blocks[blk].first_equation+l], 0 }))
                    {
                      blocks[blk].linear = false;
                      goto the_end;
                    }
              }
          }
      the_end:
        break;
      case BlockSimulationType::solveTwoBoundariesComplete:
      case BlockSimulationType::solveTwoBoundariesSimple:
        for (const auto &[indices, d1] : blocks_derivatives[blk])
          {
            int lag = get<2>(indices);
            set<pair<int, int>> endogenous;
            d1->collectEndogenous(endogenous);
            for (int l = 0; l < blocks[blk].size; l++)
              if (endogenous.contains({ endo_idx_block2orig[blocks[blk].first_equation+l], lag }))
                {
                  blocks[blk].linear = false;
                  goto the_end2;
                }
          }
      the_end2:
        break;
      default:
        break;
      }
}

int
ModelTree::equation_number() const
{
  return (equations.size());
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

  // Compute the sparse representation of the Jacobian
  for (const auto &[indices, d1] : derivatives[1])
    jacobian_sparse_column_major_order.emplace(pair{indices[0], getJacobianCol(indices[1], true)}, d1);
  jacobian_sparse_colptr = computeCSCColPtr(jacobian_sparse_column_major_order, getJacobianColsNbr(true));

  // Higher-order derivatives
  for (int o = 2; o <= order; o++)
    for (const auto &[lower_indices, lower_d] : derivatives[o-1])
      for (int var : vars)
        {
          if (lower_indices.back() > var)
            continue;

          expr_t d = lower_d->getDerivative(var);
          if (d == Zero)
            continue;

          vector<int> indices{lower_indices};
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
  /* Ensure that we don’t have any model-local variable in the model at this
     point (we used to treat them as temporary terms) */
  assert([&]
  {
    set<int> used_local_vars;
    for (auto &equation : equations)
      equation->collectVariables(SymbolType::modelLocalVariable, used_local_vars);
    return used_local_vars.empty();
  }());

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
      erase_if(it.second,
               [](expr_t e) { return !dynamic_cast<AbstractExternalFunctionNode *>(e); });

  // Fill the structures
  temporary_terms_derivatives.clear();
  temporary_terms_derivatives.resize(derivatives.size());
  for (int order = 0; order < static_cast<int>(derivatives.size()); order++)
    temporary_terms_derivatives[order] = move(temp_terms_map[{ 0, order }]);

  // Compute indices in MATLAB/Julia vector
  for (int order {0}, idx {0}; order < static_cast<int>(derivatives.size()); order++)
    for (auto it : temporary_terms_derivatives[order])
      temporary_terms_idxs[it] = idx++;
}

void
ModelTree::computeBlockTemporaryTerms(bool no_tmp_terms)
{
  int nb_blocks = blocks.size();
  blocks_temporary_terms.resize(nb_blocks);

  map<expr_t, tuple<int, int, int>> reference_count;
  for (int blk = 0; blk < nb_blocks; blk++)
    {
      blocks_temporary_terms[blk].resize(blocks[blk].size + 1);
      for (int eq = 0; eq < blocks[blk].size; eq++)
        {
          /* It is important to compute temporary terms of the renormalized
             equation if the latter is going to be used in the output files.
             Otherwise, for an equation of the form log(x) = RHS, a temporary
             term could be associated to log(x), and since it would be
             associated to this equation, it would be printed and thus computed
             *before* x is actually evaluated, and thus would be incorrect. */
          if ((blocks[blk].simulation_type == BlockSimulationType::evaluateBackward
               || blocks[blk].simulation_type == BlockSimulationType::evaluateForward
               || eq < blocks[blk].getRecursiveSize())
              && isBlockEquationRenormalized(blk, eq))
            getBlockEquationRenormalizedExpr(blk, eq)->computeBlockTemporaryTerms(blk, eq, blocks_temporary_terms, reference_count);
          else
            getBlockEquationExpr(blk, eq)->computeBlockTemporaryTerms(blk, eq, blocks_temporary_terms, reference_count);
        }
      for (const auto &[ignore, d] : blocks_derivatives[blk])
        d->computeBlockTemporaryTerms(blk, blocks[blk].size, blocks_temporary_terms, reference_count);

      additionalBlockTemporaryTerms(blk, blocks_temporary_terms, reference_count);
    }

  /* If the user has specified the notmpterms option, clear all temporary
     terms, except those that correspond to external functions (since they are
     not optional) */
  if (no_tmp_terms)
    for (auto &it : blocks_temporary_terms)
      for (auto &it2 : it)
        erase_if(it2, [](expr_t e) { return !dynamic_cast<AbstractExternalFunctionNode *>(e); });

  // Compute indices in the temporary terms vector
  blocks_temporary_terms_idxs.clear();
  for (int idx{0};
       auto &blk_tt : blocks_temporary_terms)
    for (auto &eq_tt : blk_tt)
      for (auto tt : eq_tt)
        blocks_temporary_terms_idxs[tt] = idx++;
}

void
ModelTree::additionalBlockTemporaryTerms([[maybe_unused]] int blk,
                                         [[maybe_unused]] vector<vector<temporary_terms_t>> &blocks_temporary_terms,
                                         [[maybe_unused]] map<expr_t, tuple<int, int, int>> &reference_count) const
{
}

void
ModelTree::writeJsonTemporaryTerms(const temporary_terms_t &tt,
                                   temporary_terms_t &temp_term_union,
                                   ostream &output,
                                   deriv_node_temp_terms_t &tef_terms, const string &concat) const
{
  // Local var used to keep track of temp nodes already written
  temporary_terms_t tt2 = temp_term_union;

  output << R"("external_functions_temporary_terms_)" << concat << R"(": [)";
  for (bool printed_term{false};
       auto it : tt)
    {
      if (dynamic_cast<AbstractExternalFunctionNode *>(it))
        {
          if (exchange(printed_term, true))
            output << ", ";
          vector<string> efout;
          it->writeJsonExternalFunctionOutput(efout, tt2, tef_terms);
          for (bool printed_efout{false};
               auto &it : efout)
            {
              if (exchange(printed_efout, true))
                output << ", ";
              output << it;
            }
        }
      tt2.insert(it);
    }

  output << "]"
         << R"(, "temporary_terms_)" << concat << R"(": [)";
  for (bool printed_term{false};
       const auto &it : tt)
    {
      if (exchange(printed_term, true))
        output << ", ";
      output << R"({"temporary_term": ")";
      it->writeJsonOutput(output, tt, tef_terms);
      output << R"(")"
             << R"(, "value": ")";
      it->writeJsonOutput(output, temp_term_union, tef_terms);
      output << R"("})" << endl;

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
                          varname = "paren32_tmp_var_" + to_string(i1++);
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
              varname = "paren32_tmp_var_" + to_string(i1++);
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
  for (int open{0};
       char i : str)
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
ModelTree::writeJsonModelLocalVariables(ostream &output, bool write_tef_terms, deriv_node_temp_terms_t &tef_terms) const
{
  /* Collect all model local variables appearing in equations, and print only
     them. Printing unused model local variables can lead to a crash (see
     ticket #101). */
  set<int> used_local_vars;

  for (auto equation : equations)
    equation->collectVariables(SymbolType::modelLocalVariable, used_local_vars);

  output << R"("model_local_variables": [)";
  for (bool printed_something{false};
       int id : local_variables_vector)
    if (used_local_vars.contains(id))
      {
        if (exchange(printed_something, true))
          output << ", ";

        expr_t value = local_variables_table.at(id);
        if (write_tef_terms)
          {
            vector<string> efout;
            value->writeJsonExternalFunctionOutput(efout, {}, tef_terms);
            for (bool printed_efout{false};
                 auto &it : efout)
              {
                if (exchange(printed_efout, true))
                  output << ", ";
                output << it;
              }

            if (!efout.empty())
              output << ", ";
          }

        output << R"({"variable": ")" << symbol_table.getName(id)
               << R"(", "value": ")";
        value->writeJsonOutput(output, {}, tef_terms);
        output << R"("})" << endl;
      }
  output << "]";
}

int
ModelTree::writeBytecodeBinFile(const filesystem::path &filename, bool is_two_boundaries) const
{
  ofstream SaveCode { filename, ios::out | ios::binary };
  if (!SaveCode.is_open())
    {
      cerr << R"(Error : Can't open file ")" << filename.string() << R"(" for writing)" << endl;
      exit(EXIT_FAILURE);
    }
  int u_count {0};
  for (const auto &[indices, d1] : derivatives[1])
    if (int deriv_id {indices[1]};
        getTypeByDerivID(deriv_id) == SymbolType::endogenous)
      {
        int eq {indices[0]};
        SaveCode.write(reinterpret_cast<char *>(&eq), sizeof eq);
        int tsid {getTypeSpecificIDByDerivID(deriv_id)};
        int lag {getLagByDerivID(deriv_id)};
        int varr {tsid + lag * symbol_table.endo_nbr()};
        SaveCode.write(reinterpret_cast<char *>(&varr), sizeof varr);
        SaveCode.write(reinterpret_cast<char *>(&lag), sizeof lag);
        int u {u_count + symbol_table.endo_nbr()};
        SaveCode.write(reinterpret_cast<char *>(&u), sizeof u);
        u_count++;
      }
  if (is_two_boundaries)
    u_count += symbol_table.endo_nbr();
  for (int j {0}; j < symbol_table.endo_nbr(); j++)
    SaveCode.write(reinterpret_cast<char *>(&j), sizeof j);
  for (int j {0}; j < symbol_table.endo_nbr(); j++)
    SaveCode.write(reinterpret_cast<char *>(&j), sizeof j);
  SaveCode.close();
  return u_count;
}

int
ModelTree::writeBlockBytecodeBinFile(ofstream &bin_file, int block) const
{
  int u_count {0};
  int block_size {blocks[block].size};
  int block_mfs {blocks[block].mfs_size};
  int block_recursive {blocks[block].getRecursiveSize()};
  BlockSimulationType simulation_type {blocks[block].simulation_type};
  bool is_two_boundaries {simulation_type == BlockSimulationType::solveTwoBoundariesComplete
                          || simulation_type == BlockSimulationType::solveTwoBoundariesSimple};
  for (const auto &[indices, ignore] : blocks_derivatives[block])
    {
      const auto &[eq, var, lag] {indices};
      if (lag != 0 && !is_two_boundaries)
        continue;
      if (eq >= block_recursive && var >= block_recursive)
        {
          int v {eq - block_recursive};
          bin_file.write(reinterpret_cast<char *>(&v), sizeof v);
          int varr {var - block_recursive + lag * block_mfs};
          bin_file.write(reinterpret_cast<char *>(&varr), sizeof varr);
          bin_file.write(reinterpret_cast<const char *>(&lag), sizeof lag);
          int u {u_count + block_mfs};
          bin_file.write(reinterpret_cast<char *>(&u), sizeof u);
          u_count++;
        }
    }

  if (is_two_boundaries)
    u_count += block_mfs;
  for (int j {block_recursive}; j < block_size; j++)
    {
      int varr {getBlockVariableID(block, j)};
      bin_file.write(reinterpret_cast<char *>(&varr), sizeof varr);
    }
  for (int j {block_recursive}; j < block_size; j++)
    {
      int eqr {getBlockEquationID(block, j)};
      bin_file.write(reinterpret_cast<char *>(&eqr), sizeof eqr);
    }
  return u_count;
}

void
ModelTree::writeLatexModelFile(const string &mod_basename, const string &latex_basename, ExprNodeOutputType output_type, bool write_equation_tags) const
{
  filesystem::create_directories(mod_basename + "/latex");

  const filesystem::path filename {mod_basename + "/latex/" + latex_basename + ".tex"},
    content_filename {mod_basename + "/latex/" + latex_basename + "_content" + ".tex"};
  ofstream output{filename, ios::out | ios::binary};
  if (!output.is_open())
    {
      cerr << "ERROR: Can't open file " << filename.string() << " for writing" << endl;
      exit(EXIT_FAILURE);
    }

  ofstream content_output{content_filename, ios::out | ios::binary};
  if (!content_output.is_open())
    {
      cerr << "ERROR: Can't open file " << content_filename.string() << " for writing" << endl;
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
      expr_t value = local_variables_table.at(id);

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
ModelTree::addEquation(expr_t eq, optional<int> lineno)
{
  auto beq = dynamic_cast<BinaryOpNode *>(eq);
  assert(beq && beq->op_code == BinaryOpcode::equal);

  equations.push_back(beq);
  equations_lineno.push_back(move(lineno));
}

void
ModelTree::findConstantEquationsWithoutMcpTag(map<VariableNode *, NumConstNode *> &subst_table) const
{
  for (size_t i = 0; i < equations.size(); i++)
    if (!equation_tags.exists(i, "mcp"))
      equations[i]->findConstantEquations(subst_table);
}

void
ModelTree::addEquation(expr_t eq, optional<int> lineno, const map<string, string> &eq_tags)
{
  equation_tags.add(equations.size(), eq_tags);
  addEquation(eq, move(lineno));
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
    if (trend_symbols_map.contains(id))
      throw TrendException(symbol_table.getName(id));
    else
      trend_symbols_map[id] = growth_factor;
}

void
ModelTree::addNonstationaryVariables(const vector<int> &nonstationary_vars, bool log_deflator, expr_t deflator) noexcept(false)
{
  for (int id : nonstationary_vars)
    if (nonstationary_symbols_map.contains(id))
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
        for (const auto &[lower_indices, lower_d] : derivatives[endoOrd])
          {
            expr_t d = lower_d->getDerivative(param);
            if (d == Zero)
              continue;
            vector<int> indices{lower_indices};
            indices.push_back(param);
            params_derivatives[{ endoOrd, 1 }][indices] = d;
          }
    }

  // Higher-order derivatives w.r.t. parameters
  for (int endoOrd = 0; endoOrd < static_cast<int>(derivatives.size()); endoOrd++)
    for (int paramOrd = 2; paramOrd <= paramsDerivsOrder; paramOrd++)
      for (const auto &[lower_indices, lower_d] : params_derivatives[{ endoOrd, paramOrd-1 }])
        for (int param : deriv_id_set)
          {
            if (lower_indices.back() > param)
              continue;

            expr_t d = lower_d->getDerivative(param);
            if (d == Zero)
              continue;
            vector<int> indices{lower_indices};
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
  for (const auto &[order, derivs] : params_derivatives)
    for (const auto &[indices, d] : derivs)
      d->computeTemporaryTerms(order, params_derivs_temporary_terms,
                               reference_count, true);

  for (int idx {0};
       const auto &[order, tts] : params_derivs_temporary_terms)
    for (const auto &tt : tts)
      params_derivs_temporary_terms_idxs[tt] = idx++;
}

bool
ModelTree::isNonstationary(int symb_id) const
{
  return nonstationary_symbols_map.contains(symb_id);
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
          lhs->writeJsonOutput(output, {}, {});
          output << R"(")";

          output << R"(, "rhs": ")";
          rhs->writeJsonOutput(output, {}, {});
          output << R"("})";
        }
      else
        {
          output << R"({"lhs": ")";
          lhs->writeJsonOutput(output, {}, {});
          output << R"(", "rhs": ")";
          rhs->writeJsonOutput(output, {}, {});
          output << R"(")";
          if (equations_lineno[eq])
            output << R"(, "line": )" << *equations_lineno[eq];

          if (auto eqtags = equation_tags.getTagsByEqn(eq);
              !eqtags.empty())
            {
              output << R"(, "tags": {)";
              for (bool printed_something{false};
                   const auto &[name, value] : eqtags)
                {
                  if (exchange(printed_something, true))
                    output << ", ";
                  output << R"(")" << name << R"(": ")" << value << R"(")";
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
  else if (mexext == "mexmaca64")
    return "maca64";
  else
    {
      cerr << "ERROR: 'mexext' option to preprocessor incorrectly set, needed with 'use_dll'" << endl;
      exit(EXIT_FAILURE);
    }
}

#ifdef __APPLE__
filesystem::path
ModelTree::findGccOnMacos(const string &mexext)
{
  const string macos_gcc_version {"12"}; // doc/manual/source/installation-and-configuration.rst
                                         // should be updated when this is changed
  char dynare_preprocessor_path[PATH_MAX];
  uint32_t size = PATH_MAX;
  filesystem::path local_gcc_path;
  if (_NSGetExecutablePath(dynare_preprocessor_path, &size) == 0)
    {
      string s = dynare_preprocessor_path;
      local_gcc_path = s.substr(0, s.find_last_of("/")) + "/../.brew/bin/gcc-" + macos_gcc_version;
    }

  // if user did not choose to install gcc locally via the pkg-installer then we need to find GNU gcc
  // homebrew binaries are located in /usr/local/bin/ on x86_64 systems and in /opt/homebrew/bin/ on arm64 systems
  if (exists(local_gcc_path))
    return local_gcc_path;
  else if (filesystem::path global_gcc_path {"/usr/local/bin/gcc-" + macos_gcc_version};
           exists(global_gcc_path) && mexext == "mexmaci64")
    return global_gcc_path;
  else if (filesystem::path global_gcc_path {"/opt/homebrew/bin/gcc-" + macos_gcc_version};
           exists(global_gcc_path) && mexext == "mexmaca64")
    return global_gcc_path;
  else
    {
      cerr << "ERROR: You must install gcc-" << macos_gcc_version
           << " on your system before using the `use_dll` option of Dynare. "
           << "If using MATLAB, you can do this via the Dynare installation package. If using Octave, you should run `brew install gcc-" << macos_gcc_version << "` in a terminal." << endl;
      exit(EXIT_FAILURE);
    }
}
#endif

filesystem::path
ModelTree::compileMEX(const filesystem::path &output_dir, const string &output_basename, const string &mexext, const vector<filesystem::path> &input_files, const filesystem::path &matlabroot, const filesystem::path &dynareroot, bool link) const
{
  assert(!mex_compilation_workers.empty());

  const string opt_flags = "-O3 -g0 --param ira-max-conflict-table-size=1 -fno-forward-propagate -fno-gcse -fno-dce -fno-dse -fno-tree-fre -fno-tree-pre -fno-tree-cselim -fno-tree-dse -fno-tree-dce -fno-tree-pta -fno-gcse-after-reload";

  filesystem::path compiler;
  ostringstream flags;
  string libs;

  if (matlabroot.empty())
    {
      cerr << "ERROR: 'matlabroot' option to preprocessor is not set, needed with 'use_dll'" << endl;
      exit(EXIT_FAILURE);
    }

  if (mexext == "mex")
    {
      // Octave
      compiler = matlabroot / "bin" / "mkoctfile";
      flags << "--mex";
#ifdef __APPLE__
      /* On macOS, enforce GCC, otherwise Clang will be used, and it does not
         accept our custom optimization flags (see dynare#1797) */
      filesystem::path gcc_path {findGccOnMacos(mexext)};
      if (setenv("CC", gcc_path.c_str(), 1) != 0)
        {
          cerr << "Can't set CC environment variable" << endl;
          exit(EXIT_FAILURE);
        }
      // We also define CXX, because that is used for linking
      if (setenv("CXX", gcc_path.c_str(), 1) != 0)
        {
          cerr << "Can't set CXX environment variable" << endl;
          exit(EXIT_FAILURE);
        }
#endif
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
      if (mexext == "mexa64")
        {
          // GNU/Linux
          flags << " -D_GNU_SOURCE -fPIC -pthread"
                << " -shared -Wl,--no-undefined -Wl,-rpath-link," << bin_dir;
          libs += " -lm";
        }
      else if (mexext == "mexw64")
        {
          // Windows
          flags << " -static-libgcc -shared";
          // Put the MinGW environment shipped with Dynare in the path
          auto mingwpath = dynareroot / "mingw64" / "bin";
          string newpath = "PATH=" + mingwpath.string() + ';' + getenv("PATH");
          /* We can’t use setenv() since it is not available on MinGW. Note
            that putenv() seems to make a copy of the string on MinGW, contrary
            to what is done on GNU/Linux and macOS. */
          if (putenv(const_cast<char *>(newpath.c_str())) != 0)
            {
              cerr << "Can't set PATH" << endl;
              exit(EXIT_FAILURE);
            }
        }
#ifdef __APPLE__
      else if (mexext == "mexmaci64" || mexext == "mexmaca64")
        {
          compiler = findGccOnMacos(mexext);
          flags << " -fno-common -Wl,-twolevel_namespace -undefined error -bundle";
          libs += " -lm";
        }
#endif
      else
        {
          cerr << "ERROR: unsupported value '" << mexext << "' for 'mexext' option" << endl;
          exit(EXIT_FAILURE);
        }
    }

  filesystem::path output_filename {output_dir / (output_basename + "." + (link ? mexext : "o"))};

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

  for (auto &f : input_files)
    cmd << f << " ";
  cmd << "-o " << output_filename << " ";

  if (link)
    {
      if (user_set_subst_libs.empty())
        cmd << libs;
      else
        cmd << user_set_subst_libs;
      if (!user_set_add_libs.empty())
        cmd << " " << user_set_add_libs;
    }
  else
    cmd << " -c";

#ifdef _WIN32
  cmd << '"';
#endif

  cout << "Compiling " << output_filename.string() << endl;

  // The prerequisites are the object files among the input files
  set<filesystem::path> prerequisites;
  copy_if(input_files.begin(), input_files.end(),
          inserter(prerequisites, prerequisites.end()), [](const auto &p)
          {
            return p.extension() == ".o";
          });

  unique_lock<mutex> lk {mex_compilation_mut};
  mex_compilation_queue.emplace_back(output_filename, prerequisites, cmd.str());
  lk.unlock();
  mex_compilation_cv.notify_one();

  return output_filename;
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
        int var { getTypeSpecificIDByDerivID(indices[1]) };
        int lag = getLagByDerivID(indices[1]);
        endo_derivatives[{ eq, var, lag }] = d1;
      }
  return endo_derivatives;
}

ModelTree::jacob_map_t
ModelTree::computeSymbolicJacobian(bool contemporaneous_only) const
{
  jacob_map_t symbolic_jacobian;
  for (int i = 0; i < static_cast<int>(equations.size()); i++)
    {
      set<pair<int, int>> endos_and_lags;
      equations[i]->collectEndogenous(endos_and_lags);
      for (const auto &[endo, lag] : endos_and_lags)
        if (!contemporaneous_only || lag == 0)
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

expr_t
ModelTree::getRHSFromLHS(expr_t lhs) const
{
  for (auto eq : equations)
    if (eq->arg1 == lhs)
      return eq->arg2;
  throw ExprNode::MatchFailureException{"Cannot find an equation with the requested LHS"};
}

void
ModelTree::writeBlockBytecodeAdditionalDerivatives([[maybe_unused]] BytecodeWriter &code_file,
                                                   [[maybe_unused]] int block,
                                                   [[maybe_unused]] const temporary_terms_t &temporary_terms_union,
                                                   [[maybe_unused]] const deriv_node_temp_terms_t &tef_terms) const
{
}

void
ModelTree::initializeMEXCompilationWorkers(int numworkers)
{
  assert(numworkers > 0);
  assert(mex_compilation_workers.empty());

  cout << "Spawning " << numworkers << " threads for compiling MEX files." << endl;

  for (int i {0}; i < numworkers; i++)
    mex_compilation_workers.emplace_back([](stop_token stoken)
    {
      unique_lock<mutex> lk {mex_compilation_mut};
      filesystem::path output;
      string cmd;

      /* Look for an object to compile, whose prerequisites are already
         compiled. If found, remove it from the queue, save the output path and
         the compilation command, and return true. Must be run under the lock. */
      auto pick_job = [&cmd, &output]
      {
        for (auto it {mex_compilation_queue.begin()}; it != mex_compilation_queue.end(); ++it)
          if (const auto &prerequisites {get<1>(*it)}; // Will become dangling after erase
              includes(mex_compilation_done.begin(), mex_compilation_done.end(),
                       prerequisites.begin(), prerequisites.end()))
            {
              output = get<0>(*it);
              cmd = get<2>(*it);
              mex_compilation_queue.erase(it);
              mex_compilation_ongoing.insert(output);
              return true;
            }
        return false;
      };

      while (!stoken.stop_requested())
        if (mex_compilation_cv.wait(lk, stoken, pick_job))
          {
            lk.unlock();
            int r { system(cmd.c_str()) };
            lk.lock();
            mex_compilation_ongoing.erase(output);
            if (r)
              mex_compilation_failed.insert(output);
            else
              mex_compilation_done.insert(output);
            /* The object just compiled may be a prerequisite for several
               other objects, so notify all waiting workers. Also needed to
               notify the main thread when in
               ModelTree::waitForMEXCompilationWorkers().*/
            mex_compilation_cv.notify_all();
          }
    });
}

void
ModelTree::waitForMEXCompilationWorkers()
{
  unique_lock<mutex> lk {mex_compilation_mut};
  mex_compilation_cv.wait(lk, [] {
    return (mex_compilation_queue.empty() && mex_compilation_ongoing.empty())
      || !mex_compilation_failed.empty(); });
  if (!mex_compilation_failed.empty())
    {
      cerr << "Compilation failed for: ";
      for (const auto &p : mex_compilation_failed)
        cerr << p.string() << " ";
      cerr << endl;
      lk.unlock(); // So that threads can process their stoken
      exit(EXIT_FAILURE);
    }
}

void
ModelTree::computingPassBlock(const eval_context_t &eval_context, bool no_tmp_terms)
{
  auto contemporaneous_jacobian = evaluateAndReduceJacobian(eval_context);
  if (!computeNonSingularNormalization(contemporaneous_jacobian))
    return;
  auto [prologue, epilogue] = computePrologueAndEpilogue();
  auto first_order_endo_derivatives = collectFirstOrderDerivativesEndogenous();
  equationTypeDetermination(first_order_endo_derivatives, mfs);
  cout << "Finding the optimal block decomposition of the " << modelClassName() << "..." << endl;
  computeBlockDecomposition(prologue, epilogue);
  reduceBlockDecomposition();
  printBlockDecomposition();
  computeChainRuleJacobian();
  determineLinearBlocks();
  computeBlockTemporaryTerms(no_tmp_terms);
  block_decomposed = true;
}

vector<int>
ModelTree::computeCSCColPtr(const SparseColumnMajorOrderMatrix &matrix, int ncols)
{
  vector<int> colptr(ncols+1, matrix.size());
  for (int k {0}, current_col {0};
       const auto &[indices, d1] : matrix)
    {
      while (indices.second >= current_col)
        colptr[current_col++] = k;
      k++;
    }
  return colptr;
}
