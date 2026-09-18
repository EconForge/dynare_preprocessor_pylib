/*
 * Copyright © 2024-2026 Dynare Team
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

#include <cstdlib>
#include <format>
#include <fstream>
#include <iostream>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/topological_sort.hpp>

#include "DynamicModel.hh"
#include "Exceptions.hh"
#include "HeterogeneousModel.hh"

HeterogeneousModel::HeterogeneousModel(SymbolTable& symbol_table_arg,
                                       NumericalConstants& num_constants_arg,
                                       ExternalFunctionsTable& external_functions_table_arg,
                                       HeterogeneityTable& heterogeneity_table_arg,
                                       DatabaseTable& database_table_arg,
                                       int heterogeneity_dimension_arg) :
    ModelTree {symbol_table_arg,        num_constants_arg,  external_functions_table_arg,
               heterogeneity_table_arg, database_table_arg, true},
    heterogeneity_dimension {heterogeneity_dimension_arg}
{
}

HeterogeneousModel&
HeterogeneousModel::operator=(const HeterogeneousModel& m)
{
  ModelTree::operator=(m);

  assert(heterogeneity_dimension == m.heterogeneity_dimension);

  deriv_id_table = m.deriv_id_table;
  inv_deriv_id_table = m.inv_deriv_id_table;

  return *this;
}

set<int>
HeterogeneousModel::getUsedParameters() const
{
  set<int> used;

  for (auto* equation : equations)
    {
      equation->collectVariables(SymbolType::parameter, used);
      equation->collectVariables(SymbolType::heterogeneousParameter, used);
    }

  return used;
}

set<int>
HeterogeneousModel::getUsedAggregateEndogenous() const
{
  set<int> used;
  for (auto* equation : equations)
    equation->collectVariables(SymbolType::endogenous, used);
  return used;
}

set<int>
HeterogeneousModel::getUsedAggregateExogenous() const
{
  set<int> used;
  for (auto* equation : equations)
    equation->collectVariables(SymbolType::exogenous, used);
  return used;
}

void
HeterogeneousModel::computeChainRuleJacobian()
{
  throw InternalCompilerException("HeterogeneousModel::computeChainRuleJacobian(): unimplemented");
}

int
HeterogeneousModel::getLegacyBlockJacobianEndoCol([[maybe_unused]] int blk,
                                                  [[maybe_unused]] int var,
                                                  [[maybe_unused]] int lead_lag) const
{
  throw InternalCompilerException(
      "HeterogeneousModel::getLegacyBlockJacobianEndoCol(): unimplemented");
}

int
HeterogeneousModel::getMFS() const
{
  throw InternalCompilerException("HeterogeneousModel::getMFS(): unimplemented");
}

void
HeterogeneousModel::computeDerivIDs()
{
  set<pair<int, int>> dynvars;

  for (auto& equation : equations)
    {
      equation->collectDynamicVariables(SymbolType::heterogeneousEndogenous, dynvars);
      equation->collectDynamicVariables(SymbolType::heterogeneousExogenous, dynvars);
      equation->collectDynamicVariables(SymbolType::endogenous, dynvars);
      equation->collectDynamicVariables(SymbolType::exogenous, dynvars);
      equation->collectDynamicVariables(SymbolType::parameter, dynvars);
    }

  for (const auto& [symb_id, lead_lag] : dynvars)
    {
      auto type {symbol_table.getType(symb_id)};
      if (isHeterogeneous(type))
        assert(symbol_table.getHeterogeneityDimension(symb_id) == heterogeneity_dimension);
      if (type == SymbolType::heterogeneousEndogenous || type == SymbolType::endogenous)
        assert(abs(lead_lag) <= 1);
      if (type == SymbolType::heterogeneousExogenous || type == SymbolType::exogenous)
        assert(lead_lag == 0);
      int deriv_id {static_cast<int>(deriv_id_table.size())};
      deriv_id_table.emplace(pair {symb_id, lead_lag}, deriv_id);
      inv_deriv_id_table.emplace_back(symb_id, lead_lag);
    }
}

void
HeterogeneousModel::checkPass() const
{
  const string dim_name {heterogeneity_table.getName(heterogeneity_dimension)};

  for (int i = 0; i < static_cast<int>(equations.size()); ++i)
    {
      set<pair<int, int>> het_exo_vars, het_endo_vars;
      equations[i]->collectDynamicVariables(SymbolType::heterogeneousExogenous, het_exo_vars);
      equations[i]->collectDynamicVariables(SymbolType::heterogeneousEndogenous, het_endo_vars);

      for (const auto& [symb_id, lead_lag] : het_exo_vars)
        if (symbol_table.getHeterogeneityDimension(symb_id) == heterogeneity_dimension
            && lead_lag < 0)
          throw ModelSemanticException(
              format("In model(heterogeneity={}), equation {}: lagged heterogeneous exogenous "
                     "variable '{}' is not supported.",
                     dim_name, i + 1, symbol_table.getName(symb_id)));

      // Check for het exo leads (not supported)
      for (const auto& [symb_id, lead_lag] : het_exo_vars)
        if (symbol_table.getHeterogeneityDimension(symb_id) == heterogeneity_dimension
            && lead_lag > 0)
          throw ModelSemanticException(
              format("In model(heterogeneity={}), equation {}: lead on heterogeneous exogenous "
                     "variable '{}({})' is not supported.",
                     dim_name, i + 1, symbol_table.getName(symb_id),
                     (lead_lag > 0 ? "+" : "") + to_string(lead_lag)));

      for (const auto& [symb_id, lead_lag] : het_endo_vars)
        if (symbol_table.getHeterogeneityDimension(symb_id) == heterogeneity_dimension
            && lead_lag < -1)
          throw ModelSemanticException(
              format("In model(heterogeneity={}), equation {}: heterogeneous endogenous variable "
                     "'{}' with lag {} is not supported (maximum lag is -1).",
                     dim_name, i + 1, symbol_table.getName(symb_id), lead_lag));

      // Check for het endo leads > 1 (not supported)
      for (const auto& [symb_id, lead_lag] : het_endo_vars)
        if (symbol_table.getHeterogeneityDimension(symb_id) == heterogeneity_dimension
            && lead_lag > 1)
          throw ModelSemanticException(
              format("In model(heterogeneity={}), equation {}: heterogeneous endogenous variable "
                     "'{}' with lead {} is not supported (maximum lead is +1).",
                     dim_name, i + 1, symbol_table.getName(symb_id), lead_lag));

      // Check for non-separable het lead/lag expressions: non-separable combinations of leads with
      // lagged states. A non-separable het lead/lag pattern occurs when an expression:
      // 1. Contains a heterogeneous endogenous variable with lead >= 1
      // 2. Contains a heterogeneous endogenous variable with lag == -1 (lagged state)
      // 3. These appear together in a non-separable subexpression
      //
      // Example violations: log(k(-1) + c(+1)), exp(a(-1) * c(+1))
      // Valid separable forms: k(-1) + c(+1), beta * c(+1) + a(-1)
      if (auto subexpr = equations[i]->findNonSeparableHetLeadLagSubexpr(heterogeneity_dimension);
          subexpr)
        {
          ostringstream oss;
          oss << "In model(heterogeneity=" << dim_name << "):\n"
              << "  Non-separable expression '" << subexpr->toString() << "'"
              << "  combines forward-looking variables with lagged states and is not supported.";
          throw ModelSemanticException(oss.str(), static_cast<int>(i) + 1, equations_lineno[i]);
        }
    }
}

/*
 * Unfold complementarity conditions: (i) declare the multipliers associated
 * with each bound constraint μ_l and μ_u ; (ii) add or substract the
 * multiplier into the associated condition; (iii) add the the complementarity
 * slackness conditions into the set of equations. For example,
 * households choose {cₜ, aₜ₊₁} to maximize expected lifetime utility:
 *     max 𝐸ₜ [∑ₛ₌₀^∞ βˢ · u(cₜ₊ₛ)]
 *
 * Subject to:
 *   1. Budget constraint:      cₜ + aₜ₊₁ = yₜ + (1 + rₜ) · aₜ
 *   2. Borrowing constraint:   aₜ₊₁ ≥ aₘᵢₙ
 *
 * Let u'(cₜ) denote the marginal utility of consumption.
 * Let μₜ ≥ 0 be the Lagrange multiplier on the borrowing constraint.
 *
 * Then, the Euler equation becomes:
 *     u′(cₜ) = β · (1 + rₜ₊₁) · u′(cₜ₊₁) − μₜ
 *
 * Together with:
 *     aₜ₊₁ ≥ aₘᵢₙ                 [primal feasibility]
 *     μₜ ≥ 0                      [dual feasibility]
 *     μₜ · (aₜ₊₁ − aₘᵢₙ) = 0      [complementarity slackness]
 * Note that the primal feasibility and dual feasibility constraints are not
 * introduced here, but Bhandari et al. (2023) show in Appendix B.1 that they
 * are redundant.
 */
void
HeterogeneousModel::transformPass()
{
  // Unfold complementarity conditions (MCP)
  for (int i = 0; i < static_cast<int>(equations.size()); ++i)
    {
      if (!complementarity_conditions[i])
        continue;

      /*
       * `const auto& [symb_id, lb, ub] = *complementarity_conditions[i];` was not used here because
       * the call to `addEquation` may eventually lead to a resize of the
       * `complementarity_conditions` vector, which may invalidate the reference to its element. We
       * take a copy instead for safety.
       */
      auto [symb_id, lb, ub] = *complementarity_conditions[i];

      VariableNode* var = getVariable(symb_id);
      if (lb)
        {
          int mu_id = symbol_table.addHeterogeneousMultiplierAuxiliaryVar(
              heterogeneity_dimension, i, "MULT_L_" + symbol_table.getName(symb_id));
          expr_t mu_L = AddVariable(mu_id);
          // For lower bound: F >= 0 ⟂ var >= lb, KKT gives F - μ = 0 with μ >= 0
          expr_t orig_resid = AddMinus(equations[i]->arg1, equations[i]->arg2);
          auto substeq = AddEqual(AddMinus(equations[i]->arg1, mu_L), equations[i]->arg2);
          assert(substeq);
          equations[i] = substeq;
          addEquation(AddEqual(AddTimes(mu_L, AddMinus(var, lb)), Zero), nullopt);
          addAuxEquation(AddEqual(mu_L, AddTimes(orig_resid, AddLessEqual(var, lb))));
        }
      if (ub)
        {
          int mu_id = symbol_table.addHeterogeneousMultiplierAuxiliaryVar(
              heterogeneity_dimension, i, "MULT_U_" + symbol_table.getName(symb_id));
          auto mu_U = AddVariable(mu_id);
          // For upper bound: F <= 0 ⟂ var <= ub, KKT gives F + μ = 0 with μ >= 0
          expr_t orig_resid = AddMinus(equations[i]->arg2, equations[i]->arg1);
          auto substeq = AddEqual(AddPlus(equations[i]->arg1, mu_U), equations[i]->arg2);
          assert(substeq);
          equations[i] = substeq;
          addEquation(AddEqual(AddTimes(mu_U, AddMinus(ub, var)), Zero), nullopt);
          addAuxEquation(AddEqual(mu_U, AddTimes(orig_resid, AddGreaterEqual(var, ub))));
        }
    }
  for (size_t idx : het_aux_equations_indices)
    aux_equations.push_back(equations[idx]);

  // Reorder auxiliary equations in dependency order
  reorderHetAuxiliaryEquations();
}

void
HeterogeneousModel::reorderHetAuxiliaryEquations()
{
  using namespace boost;

  if (aux_equations.empty())
    return;

  // Create the mapping between auxiliary variables and auxiliary equations
  int n = static_cast<int>(aux_equations.size());
  map<int, int> auxHetEndoToEq;
  for (int i = 0; i < n; i++)
    {
      auto varexpr = dynamic_cast<VariableNode*>(aux_equations[i]->arg1);
      assert(varexpr
             && symbol_table.getType(varexpr->symb_id) == SymbolType::heterogeneousEndogenous);
      auxHetEndoToEq[varexpr->symb_id] = i;
    }
  assert(static_cast<int>(auxHetEndoToEq.size()) == n);

  /* Construct the directed acyclic graph where auxiliary equations are
     vertices and edges represent dependency relationships. */
  using Graph = adjacency_list<vecS, vecS, directedS>;
  Graph g(n);
  for (int i = 0; i < n; i++)
    {
      set<int> het_endos;
      aux_equations[i]->collectVariables(SymbolType::heterogeneousEndogenous, het_endos);
      for (int het_endo : het_endos)
        if (auto it = auxHetEndoToEq.find(het_endo); it != auxHetEndoToEq.end() && it->second != i)
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

void
HeterogeneousModel::computeHetAuxTopologicalLevels()
{
  het_aux_levels.clear();

  if (aux_equations.empty())
    return;

  int n = static_cast<int>(aux_equations.size());

  // Build mapping from aux var symbol ID to equation index
  map<int, int> auxSymbToEqIdx;
  for (int i = 0; i < n; i++)
    {
      auto varexpr = dynamic_cast<VariableNode*>(aux_equations[i]->arg1);
      assert(varexpr);
      auxSymbToEqIdx[varexpr->symb_id] = i;
    }

  // For each aux equation, find which aux vars it depends on at (+1) or (-1) timing
  // A variable at (+1) requires E[.] computation; (-1) requires time-shifted values
  vector<set<int>> deps(n); // deps[i] = set of equation indices that i depends on at (+1) or (-1)

  for (int i = 0; i < n; i++)
    {
      // Collect all heterogeneous endogenous in the RHS with their (symb_id, lag)
      set<pair<int, int>> het_endos;
      aux_equations[i]->arg2->collectDynamicVariables(SymbolType::heterogeneousEndogenous,
                                                      het_endos);

      for (const auto& [symb_id, lag] : het_endos)
        {
          if (lag != 1 && lag != -1)
            continue; // Only interested in (+1) or (-1) timing

          // Check if this symbol is an auxiliary variable
          if (auto it = auxSymbToEqIdx.find(symb_id); it != auxSymbToEqIdx.end() && it->second != i)
            deps[i].insert(it->second);
        }
    }

  // Compute levels iteratively
  // Level 0: aux vars with no (+1) or (-1) aux dependencies
  // Level k: aux vars whose (+1) or (-1) dependencies are all at level < k
  vector<int> levels(n, -1); // -1 means unassigned
  int current_level = 0;
  int assigned = 0;

  while (assigned < n)
    {
      vector<int> level_vars;

      for (int i = 0; i < n; i++)
        {
          if (levels[i] >= 0)
            continue; // Already assigned

          // Check if all dependencies are at level < current_level
          if (ranges::all_of(deps[i], [&](int dep) {
                return levels[dep] >= 0 && levels[dep] < current_level;
              }))
            {
              levels[i] = current_level;
              // Get the het endo type-specific ID for this aux var (1-based for MATLAB)
              auto varexpr = dynamic_cast<VariableNode*>(aux_equations[i]->arg1);
              int tsid = symbol_table.getTypeSpecificID(varexpr->symb_id) + 1;
              level_vars.push_back(tsid);
              assigned++;
            }
        }

      if (!level_vars.empty())
        het_aux_levels.push_back(level_vars);

      current_level++;

      // Safety check: if no progress made, there's a cycle
      if (level_vars.empty() && assigned < n)
        throw ModelSemanticException(
            "Circular dependency detected in heterogeneous auxiliary variables");
    }
}

void
HeterogeneousModel::substituteEndoLead(DynamicModel& dynamic_model)
{
  substituteLeadLagInternal(dynamic_model, AuxVarType::endoLead);
  substituteLeadLagInternal(dynamic_model, AuxVarType::heterogeneousEndoLead);
}

void
HeterogeneousModel::substituteEndoLagGreaterThanTwo(DynamicModel& dynamic_model)
{
  substituteLeadLagInternal(dynamic_model, AuxVarType::endoLag);
}

void
HeterogeneousModel::substituteExoLead(DynamicModel& dynamic_model)
{
  substituteLeadLagInternal(dynamic_model, AuxVarType::exoLead);
}

void
HeterogeneousModel::substituteExoLag(DynamicModel& dynamic_model)
{
  substituteLeadLagInternal(dynamic_model, AuxVarType::exoLag);
}

void
HeterogeneousModel::substituteLeadLagInternal(DynamicModel& dynamic_model, AuxVarType type)
{
  ExprNode::subst_table_t subst_table;
  vector<BinaryOpNode*> neweqs;

  // Substitute in used model local variables
  set<int> used_local_vars;
  for (auto& equation : equations)
    equation->collectVariables(SymbolType::modelLocalVariable, used_local_vars);

  for (int used_local_var : used_local_vars)
    {
      const expr_t value = local_variables_table.at(used_local_var);
      expr_t subst;
      switch (type)
        {
        case AuxVarType::heterogeneousEndoLead:
          subst
              = value->substituteHetEndoNonLinearLead(heterogeneity_dimension, subst_table, neweqs);
          break;
        case AuxVarType::endoLead:
          subst = value->substituteEndoLeadGreaterThanTwo(subst_table, neweqs, true);
          break;
        case AuxVarType::endoLag:
          subst = value->substituteEndoLagGreaterThanTwo(subst_table, neweqs);
          break;
        case AuxVarType::exoLead:
          subst = value->substituteExoLead(subst_table, neweqs, true);
          break;
        case AuxVarType::exoLag:
          subst = value->substituteExoLag(subst_table, neweqs);
          break;
        default:
          throw InternalCompilerException(
              "DynamicModel::substituteLeadLagInternal: impossible case");
        }
      local_variables_table[used_local_var] = subst;
    }

  // Substitute in equations
  for (auto& equation : equations)
    {
      expr_t subst;
      switch (type)
        {
        case AuxVarType::heterogeneousEndoLead:
          subst = equation->substituteHetEndoNonLinearLead(heterogeneity_dimension, subst_table,
                                                           neweqs);
          break;
        case AuxVarType::endoLead:
          subst = equation->substituteEndoLeadGreaterThanTwo(subst_table, neweqs, true);
          break;
        case AuxVarType::endoLag:
          subst = equation->substituteEndoLagGreaterThanTwo(subst_table, neweqs);
          break;
        case AuxVarType::exoLead:
          subst = equation->substituteExoLead(subst_table, neweqs, true);
          break;
        case AuxVarType::exoLag:
          subst = equation->substituteExoLag(subst_table, neweqs);
          break;
        default:
          throw InternalCompilerException(
              "HeterogeneousModel::substituteLeadLagInternal: impossible case");
        }
      auto substeq = dynamic_cast<BinaryOpNode*>(subst);
      assert(substeq);
      equation = substeq;
    }

  // Add new equations
  for (auto& neweq : neweqs)
    {
      if (isHeterogeneousAux(type))
        {
          addEquation(neweq, nullopt);
          het_aux_equations_indices.push_back(equations.size() - 1);
        }
      else
        {
          dynamic_model.addEquation(dynamic_cast<BinaryOpNode*>(neweq->clone(dynamic_model)),
                                    nullopt);
          dynamic_model.addAuxEquation(dynamic_cast<BinaryOpNode*>(neweq->clone(dynamic_model)));
        }
    }

  if (neweqs.size() > 0)
    {
      cout << "Substitution of ";
      switch (type)
        {
        case AuxVarType::heterogeneousEndoLead:
          cout << "het endo leads (non-separable >= 1, vars >= 2)";
          break;
        case AuxVarType::endoLead:
          cout << "endo leads >= 2";
          break;
        case AuxVarType::endoLag:
          cout << "endo lags >= 2";
          break;
        case AuxVarType::exoLead:
          cout << "exo leads";
          break;
        case AuxVarType::exoLag:
          cout << "exo lags";
          break;
        default:
          throw InternalCompilerException(
              "DynamicModel::substituteLeadLagInternal: impossible case");
        }
      cout << ": added " << neweqs.size() << " auxiliary variables and equations." << '\n';
    }
}

void
HeterogeneousModel::computingPass(int derivsOrder, bool no_tmp_terms, bool use_dll)
{
  assert(!use_dll); // Not yet implemented

  // Compute topological levels for auxiliary variables (requires frozen symbol table)
  computeHetAuxTopologicalLevels();

  computeDerivIDs();

  set<int> vars;
  for (auto& [symb_lag, deriv_id] : deriv_id_table)
    if (symbol_table.getType(symb_lag.first) != SymbolType::parameter)
      vars.insert(deriv_id);

  cout << "Computing " << modelClassName() << " derivatives (order " << derivsOrder << ")." << '\n'
       << flush;

  computeDerivatives(derivsOrder, vars);

  computeTemporaryTerms(!use_dll, no_tmp_terms);

  computeMCPEquationsReordering(heterogeneity_dimension);
}

void
HeterogeneousModel::writeModelFiles(const string& basename, bool julia) const
{
  assert(!julia); // Not yet implemented
  writeModelMFiles<true>(basename, heterogeneity_dimension);
  writeComplementarityConditionsFile<true>(basename, heterogeneity_dimension);
  writeSetHetAuxiliaryVariablesFile(
      basename); // Always call; early return inside if nothing to write
}

void
HeterogeneousModel::writeSetHetAuxiliaryVariablesFile(const string& basename) const
{
  if (aux_equations.empty())
    return;

  string filename {
      (packageDir(basename)
       / format("dynamic_het{}_set_auxiliary_variables.m", heterogeneity_dimension + 1))
          .string()};
  ofstream output {filename, ios::out | ios::binary};
  if (!output.is_open())
    throw FileIOException(filename, "writing");

  output << "function yh = dynamic_het" << heterogeneity_dimension + 1
         << "_set_auxiliary_variables(y, x, params, steady_state, yh, xh, paramsh, step)" << '\n'
         << "%" << '\n'
         << "% Sets auxiliary variables for heterogeneous model dimension "
         << heterogeneity_dimension + 1 << '\n'
         << "% step: level of auxiliary variables to compute (0-based)" << '\n'
         << "%" << '\n'
         << '\n';

  // Build mapping from aux var type-specific ID to equation index
  map<int, int> tsidToEqIdx;
  for (size_t i = 0; i < aux_equations.size(); i++)
    {
      auto varexpr = dynamic_cast<VariableNode*>(aux_equations[i]->arg1);
      assert(varexpr);
      int tsid = symbol_table.getTypeSpecificID(varexpr->symb_id) + 1; // 1-based
      tsidToEqIdx[tsid] = static_cast<int>(i);
    }

  // Write auxiliary equations grouped by level
  deriv_node_temp_terms_t tef_terms;

  // First pass: write any external function outputs
  for (auto aux_equation : aux_equations)
    if (aux_equation->containsExternalFunction())
      aux_equation->writeExternalFunctionOutput(output, ExprNodeOutputType::matlabDynamicModel, {},
                                                {}, tef_terms);

  // Second pass: write equations grouped by level with conditionals
  for (size_t level = 0; level < het_aux_levels.size(); level++)
    {
      output << "if step == " << level << '\n';

      for (int tsid : het_aux_levels[level])
        {
          auto it = tsidToEqIdx.find(tsid);
          if (it != tsidToEqIdx.end())
            {
              int eq_idx = it->second;
              output << "    ";
              aux_equations[eq_idx]->writeOutput(output, ExprNodeOutputType::matlabDynamicModel, {},
                                                 {}, tef_terms);
              output << ";" << '\n';
            }
        }

      output << "end" << '\n' << '\n';
    }

  output << "end" << '\n';
  output.close();
}

int
HeterogeneousModel::getJacobianCol(int deriv_id) const
{
  SymbolType type {getTypeByDerivID(deriv_id)};
  int tsid {getTypeSpecificIDByDerivID(deriv_id)};
  int lag {getLagByDerivID(deriv_id)};

  if (type == SymbolType::heterogeneousEndogenous)
    return tsid + (lag + 1) * symbol_table.het_endo_nbr(heterogeneity_dimension);

  int shift {3 * symbol_table.het_endo_nbr(heterogeneity_dimension)};

  if (type == SymbolType::heterogeneousExogenous)
    return shift + tsid;

  shift += symbol_table.het_exo_nbr(heterogeneity_dimension);

  if (type == SymbolType::endogenous)
    return shift + tsid + (lag + 1) * symbol_table.endo_nbr();

  shift += 3 * symbol_table.endo_nbr();

  if (type == SymbolType::exogenous)
    return shift + tsid;

  throw UnknownDerivIDException();
}

int
HeterogeneousModel::getJacobianColsNbr() const
{
  return 3 * (symbol_table.het_endo_nbr(heterogeneity_dimension) + symbol_table.endo_nbr())
         + symbol_table.het_exo_nbr(heterogeneity_dimension) + symbol_table.exo_nbr();
}

int
HeterogeneousModel::getLegacyJacobianCol([[maybe_unused]] int deriv_id) const
{
  throw InternalCompilerException("HeterogeneousModel::getLegacyJacobianCol(): unimplemented");
}

SymbolType
HeterogeneousModel::getTypeByDerivID(int deriv_id) const noexcept(false)
{
  return symbol_table.getType(getSymbIDByDerivID(deriv_id));
}

int
HeterogeneousModel::getLagByDerivID(int deriv_id) const noexcept(false)
{
  if (deriv_id < 0 || deriv_id >= static_cast<int>(inv_deriv_id_table.size()))
    throw UnknownDerivIDException();

  return inv_deriv_id_table[deriv_id].second;
}

int
HeterogeneousModel::getSymbIDByDerivID(int deriv_id) const noexcept(false)
{
  if (deriv_id < 0 || deriv_id >= static_cast<int>(inv_deriv_id_table.size()))
    throw UnknownDerivIDException();

  return inv_deriv_id_table[deriv_id].first;
}

int
HeterogeneousModel::getTypeSpecificIDByDerivID(int deriv_id) const
{
  return symbol_table.getTypeSpecificID(getSymbIDByDerivID(deriv_id));
}

int
HeterogeneousModel::getDerivID(int symb_id, int lead_lag) const noexcept(false)
{
  if (auto it = deriv_id_table.find({symb_id, lead_lag}); it == deriv_id_table.end())
    throw UnknownDerivIDException();
  else
    return it->second;
}

void
HeterogeneousModel::writeDriverOutput(ostream& output) const
{
  std::vector<int> state_var;
  for (int endoID = 0; endoID < symbol_table.het_orig_endo_nbr(heterogeneity_dimension); endoID++)
    try
      {
        // The return value of the following call is ignored, since we care about the exception
        // NOLINTNEXTLINE(clang-diagnostic-unused-result)
        getDerivID(symbol_table.getID(SymbolType::heterogeneousEndogenous, endoID,
                                      heterogeneity_dimension),
                   -1);
        if (ranges::find(state_var, endoID) == state_var.end())
          state_var.push_back(endoID);
      }
    catch (UnknownDerivIDException& e)
      {
      }

  output << "M_.heterogeneity(" << heterogeneity_dimension + 1 << ").state_var = [";
  for (int it : state_var)
    output << it + 1 << " ";
  output << "];" << '\n';

  output << "M_.heterogeneity(" << heterogeneity_dimension + 1 << ").dynamic_tmp_nbr = [";
  for (const auto& it : temporary_terms_derivatives)
    output << it.size() << "; ";
  output << "];" << '\n';
  writeDriverSparseIndicesHelper(format("heterogeneity({}).dynamic", heterogeneity_dimension + 1),
                                 output);
  output << "M_.heterogeneity(" << heterogeneity_dimension + 1
         << ").dynamic_mcp_equations_reordering = [";
  for (auto i : mcp_equations_reordering)
    output << i + 1 << "; ";
  output << "];" << '\n';

  output << "M_.heterogeneity(" << heterogeneity_dimension + 1
         << ").set_auxiliary_variables = exist(['./+' M_.fname '/dynamic_het"
         << heterogeneity_dimension + 1 << "_set_auxiliary_variables.m'], 'file') == 2;" << '\n';

  // Output auxiliary variable level structure
  output << "M_.heterogeneity(" << heterogeneity_dimension + 1
         << ").n_aux_levels = " << het_aux_levels.size() << ";" << '\n';

  output << "M_.heterogeneity(" << heterogeneity_dimension + 1 << ").het_aux_levels = {";
  for (size_t level = 0; level < het_aux_levels.size(); level++)
    {
      output << "[";
      for (size_t i = 0; i < het_aux_levels[level].size(); i++)
        {
          if (i > 0)
            output << " ";
          output << het_aux_levels[level][i];
        }
      output << "]";
      if (level < het_aux_levels.size() - 1)
        output << ", ";
    }
  output << "};" << '\n';

  equation_tags.writeOutput(
      output, format("M_.heterogeneity({}).equations_tags", heterogeneity_dimension + 1));
}
