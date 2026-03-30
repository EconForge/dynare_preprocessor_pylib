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

#ifndef HETEROGENEOUS_MODEL_HH
#define HETEROGENEOUS_MODEL_HH

#include <set>
#include <string>

#include "ModelTree.hh"

using namespace std;

class DynamicModel;

class HeterogeneousModel : public ModelTree
{
private:
  vector<size_t> het_aux_equations_indices;

  // Auxiliary variables grouped by topological level
  // Level 0: aux vars with no (+1) dependencies on other aux vars
  // Level k: aux vars whose (+1) dependencies are at level < k
  vector<vector<int>> het_aux_levels;

  //! Factorized code for substitutions of leads/lags
  void substituteLeadLagInternal(DynamicModel& dynamic_model, AuxVarType type);

public:
  const int heterogeneity_dimension;

  HeterogeneousModel(SymbolTable& symbol_table_arg, NumericalConstants& num_constants_arg,
                     ExternalFunctionsTable& external_functions_table_arg,
                     HeterogeneityTable& heterogeneity_table_arg, DatabaseTable& database_table_arg,
                     int heterogeneity_dimension_arg);

  HeterogeneousModel(const HeterogeneousModel& m) = default;
  HeterogeneousModel& operator=(const HeterogeneousModel& m);

  // Checks for unsupported lead/lag patterns in user-written het equations
  void checkPass() const;

  void transformPass();

  void computingPass(int derivsOrder, bool no_tmp_terms, bool use_dll);

  void writeModelFiles(const string& basename, bool julia) const;
  void writeDriverOutput(ostream& output) const;

  [[nodiscard]] set<int> getUsedParameters() const;
  [[nodiscard]] set<int> getUsedAggregateEndogenous() const;
  [[nodiscard]] set<int> getUsedAggregateExogenous() const;

  [[nodiscard]] int getJacobianCol(int deriv_id) const override;
  [[nodiscard]] int getJacobianColsNbr() const override;
  [[nodiscard]] int getLegacyJacobianCol(int deriv_id) const override;

  /* These methods substitute variables with leads/lags in heterogeneous model
     equations, replacing them with auxiliary variables.
     For aggregate endogenous and exogenous leads, substitution uses deterministic mode.
     For heterogeneous endogenous leads, a specialized substitution handles both
     individual variables (lead >= 2) and non-separable expressions (lead >= 1). */

  /*! Substitutes:
     - aggregate endogenous variables with lead >= 2 (deterministic mode).
     - heterogeneous endogenous variables with lead >= 2 and
       heterogeneous endogenous non-separable expressions with lead >= 1. */
  void substituteEndoLead(DynamicModel& dynamic_model);

  //! Substitutes aggregate endogenous with lag >= 2
  void substituteEndoLagGreaterThanTwo(DynamicModel& dynamic_model);

  //! Substitutes aggregate exogenous with any lead (deterministic mode)
  void substituteExoLead(DynamicModel& dynamic_model);

  //! Substitutes aggregate exogenous with any lag
  void substituteExoLag(DynamicModel& dynamic_model);

  // FIXME: the following 5 functions are identical to those in DynamicModel. Factorization?
  [[nodiscard]] int getDerivID(int symb_id, int lead_lag) const noexcept(false) override;
  [[nodiscard]] SymbolType getTypeByDerivID(int deriv_id) const noexcept(false) override;
  [[nodiscard]] int getLagByDerivID(int deriv_id) const noexcept(false) override;
  [[nodiscard]] int getSymbIDByDerivID(int deriv_id) const noexcept(false) override;
  [[nodiscard]] int getTypeSpecificIDByDerivID(int deriv_id) const override;

protected:
  void computeChainRuleJacobian() override;
  int getLegacyBlockJacobianEndoCol(int blk, int var, int lead_lag) const override;
  string
  modelClassName() const override
  {
    return "dynamic model for heterogeneity dimension '"
           + heterogeneity_table.getName(heterogeneity_dimension) + "'";
  }
  int getMFS() const override;

private:
  // Maps a pair (symbol ID, lead/lag) to a deriv ID
  map<pair<int, int>, int> deriv_id_table;
  // Maps a deriv ID to a pair (symbol ID, lead/lag)
  vector<pair<int, int>> inv_deriv_id_table;

  // Reorders auxiliary equations based on dependencies (adapted for heterogeneousEndogenous)
  void reorderHetAuxiliaryEquations();

  // Computes topological levels for auxiliary variables based on (+1) dependencies
  void computeHetAuxTopologicalLevels();

  // Allocates the derivation IDs for all endogenous variables for this heterogeneity dimension
  void computeDerivIDs();

  // Writes the file for setting heterogeneous auxiliary variables
  void writeSetHetAuxiliaryVariablesFile(const string& basename) const;
};

#endif
