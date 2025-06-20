/*
 * Copyright © 2024-2025 Dynare Team
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
#include <iostream>

#include "HeterogeneousModel.hh"

HeterogeneousModel::HeterogeneousModel(SymbolTable& symbol_table_arg,
                                       NumericalConstants& num_constants_arg,
                                       ExternalFunctionsTable& external_functions_table_arg,
                                       HeterogeneityTable& heterogeneity_table_arg,
                                       int heterogeneity_dimension_arg) :
    ModelTree {symbol_table_arg, num_constants_arg, external_functions_table_arg,
               heterogeneity_table_arg, true},
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

void
HeterogeneousModel::computeChainRuleJacobian()
{
  cerr << "Heterogeneous::computeChainRuleJacobian(): unimplemented" << endl;
  exit(EXIT_FAILURE);
}

int
HeterogeneousModel::getBlockJacobianEndoCol([[maybe_unused]] int blk, [[maybe_unused]] int var,
                                            [[maybe_unused]] int lead_lag) const
{
  cerr << "Heterogeneous::getBlockJacobianEndoCol(): unimplemented" << endl;
  exit(EXIT_FAILURE);
}

int
HeterogeneousModel::getMFS() const
{
  cerr << "Heterogeneous::getMFS(): unimplemented" << endl;
  exit(EXIT_FAILURE);
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
          auto substeq = AddEqual(AddPlus(equations[i]->arg1, mu_L), equations[i]->arg2);
          assert(substeq);
          equations[i] = substeq;
          addEquation(AddEqual(AddTimes(mu_L, AddMinus(var, lb)), Zero), nullopt);
        }
      if (ub)
        {
          int mu_id = symbol_table.addHeterogeneousMultiplierAuxiliaryVar(
              heterogeneity_dimension, i, "MULT_U_" + symbol_table.getName(symb_id));
          auto mu_U = AddVariable(mu_id);
          auto substeq = AddEqual(AddMinus(equations[i]->arg1, mu_U), equations[i]->arg2);
          assert(substeq);
          equations[i] = substeq;
          addEquation(AddEqual(AddTimes(mu_U, AddMinus(ub, var)), Zero), nullopt);
        }
    }
}

void
HeterogeneousModel::computingPass(int derivsOrder, bool no_tmp_terms, bool use_dll)
{
  assert(!use_dll); // Not yet implemented

  computeDerivIDs();

  set<int> vars;
  for (auto& [symb_lag, deriv_id] : deriv_id_table)
    if (symbol_table.getType(symb_lag.first) != SymbolType::parameter)
      vars.insert(deriv_id);

  cout << "Computing " << modelClassName() << " derivatives (order " << derivsOrder << ")." << endl;

  computeDerivatives(derivsOrder, vars);

  computeTemporaryTerms(!use_dll, no_tmp_terms);

  computeMCPEquationsReordering(heterogeneity_dimension);
}

void
HeterogeneousModel::writeModelFiles(const string& basename, bool julia) const
{
  assert(!julia); // Not yet implemented
  writeSparseModelMFiles<true>(basename, heterogeneity_dimension);
  writeComplementarityConditionsFile<true>(basename, heterogeneity_dimension);
}

int
HeterogeneousModel::getJacobianCol(int deriv_id, bool sparse) const
{
  assert(sparse);

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

  shift += symbol_table.endo_nbr();

  if (type == SymbolType::exogenous)
    return shift + tsid;

  throw UnknownDerivIDException();
}

int
HeterogeneousModel::getJacobianColsNbr(bool sparse) const
{
  assert(sparse);
  return 3 * (symbol_table.het_endo_nbr(heterogeneity_dimension) + symbol_table.endo_nbr())
         + symbol_table.het_exo_nbr(heterogeneity_dimension) + symbol_table.exo_nbr();
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
  for (int endoID = 0; endoID < symbol_table.het_endo_nbr(heterogeneity_dimension); endoID++)
    try
      {
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
  output << "];" << endl;

  output << "M_.heterogeneity(" << heterogeneity_dimension + 1 << ").dynamic_tmp_nbr = [";
  for (const auto& it : temporary_terms_derivatives)
    output << it.size() << "; ";
  output << "];" << endl;
  writeDriverSparseIndicesHelper(
      "heterogeneity("s + to_string(heterogeneity_dimension + 1) + ").dynamic", output);
  output << "M_.heterogeneity(" << heterogeneity_dimension + 1
         << ").dynamic_mcp_equations_reordering = [";
  for (auto i : mcp_equations_reordering)
    output << i + 1 << "; ";
  output << "];" << endl;
}
