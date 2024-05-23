/*
 * Copyright © 2024 Dynare Team
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

HeterogeneousModel::HeterogeneousModel(const HeterogeneousModel& m) :
    ModelTree {m},
    heterogeneity_dimension {m.heterogeneity_dimension},
    deriv_id_table {m.deriv_id_table},
    inv_deriv_id_table {m.inv_deriv_id_table}
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

  if (ranges::any_of(complementarity_conditions, [](const auto& x) { return x.has_value(); }))
    {
      // Implementing it requires modifications in ModelTree::computeMCPEquationsReordering()
      cerr << "ERROR: Complementarity conditions are not yet implemented in "
              "model(heterogeneity=...) blocks"
           << endl;
      exit(EXIT_FAILURE);
    }
}

void
HeterogeneousModel::writeModelFiles(const string& basename, bool julia) const
{
  assert(!julia); // Not yet implemented
  writeSparseModelMFiles<true>(basename, heterogeneity_dimension);
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
  output << "M_.heterogeneity(" << heterogeneity_dimension + 1 << ").dynamic_tmp_nbr = [";
  for (const auto& it : temporary_terms_derivatives)
    output << it.size() << "; ";
  output << "];" << endl;
  writeDriverSparseIndicesHelper(
      "heterogeneity("s + to_string(heterogeneity_dimension + 1) + ").dynamic", output);
}
