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

#ifndef HETEROGENEOUS_MODEL_HH
#define HETEROGENEOUS_MODEL_HH

#include <string>

#include "ModelTree.hh"

using namespace std;

class HeterogeneousModel : public ModelTree
{
public:
  const int heterogeneity_dimension;

  HeterogeneousModel(SymbolTable& symbol_table_arg, NumericalConstants& num_constants_arg,
                     ExternalFunctionsTable& external_functions_table_arg,
                     HeterogeneityTable& heterogeneity_table_arg, int heterogeneity_dimension_arg);

  HeterogeneousModel(const HeterogeneousModel& m) = default;
  HeterogeneousModel& operator=(const HeterogeneousModel& m);

  void transformPass();

  void computingPass(int derivsOrder, bool no_tmp_terms, bool use_dll);

  void writeModelFiles(const string& basename, bool julia) const;
  void writeDriverOutput(ostream& output) const;

  [[nodiscard]] int getJacobianCol(int deriv_id, bool sparse) const override;
  [[nodiscard]] int getJacobianColsNbr(bool sparse) const override;

#if 0
  void substituteEndoLeadGreaterThanTwo();

  //! Transforms the model by removing all lags greater or equal than 2 on endos
  void substituteEndoLagGreaterThanTwo();

  //! Transforms the model by removing all leads on exos
  /*! Note that this can create new lags on endos and exos */
  void substituteExoLead();

  //! Transforms the model by removing all lags on exos
  void substituteExoLag();

  //! Transforms the model by removing all UnaryOpcode::expectation
  void substituteExpectation(bool partial_information_model);

  //! Transforms the model by decreasing the lead/lag of predetermined variables in model equations
  //! by one
  void transformPredeterminedVariables();

  //! Substitutes out all model-local variables
  void substituteModelLocalVariables();
#endif

  // FIXME: the following 5 functions are identical to those in DynamicModel. Factorization?
  [[nodiscard]] int getDerivID(int symb_id, int lead_lag) const noexcept(false) override;
  [[nodiscard]] SymbolType getTypeByDerivID(int deriv_id) const noexcept(false) override;
  [[nodiscard]] int getLagByDerivID(int deriv_id) const noexcept(false) override;
  [[nodiscard]] int getSymbIDByDerivID(int deriv_id) const noexcept(false) override;
  [[nodiscard]] int getTypeSpecificIDByDerivID(int deriv_id) const override;

protected:
  void computeChainRuleJacobian() override;
  int getBlockJacobianEndoCol(int blk, int var, int lead_lag) const override;
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

  // Allocates the derivation IDs for all endogenous variables for this heterogeneity dimension
  void computeDerivIDs();
};

#endif
