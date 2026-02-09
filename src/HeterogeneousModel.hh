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
public:
  const int heterogeneity_dimension;

  HeterogeneousModel(SymbolTable& symbol_table_arg, NumericalConstants& num_constants_arg,
                     ExternalFunctionsTable& external_functions_table_arg,
                     HeterogeneityTable& heterogeneity_table_arg, DatabaseTable& database_table_arg,
                     int heterogeneity_dimension_arg);

  HeterogeneousModel(const HeterogeneousModel& m) = default;
  HeterogeneousModel& operator=(const HeterogeneousModel& m);

  void transformPass();

  void computingPass(int derivsOrder, bool no_tmp_terms, bool use_dll);

  void writeModelFiles(const string& basename, bool julia) const;
  void writeDriverOutput(ostream& output) const;

  [[nodiscard]] set<int> getUsedParameters() const;

  [[nodiscard]] int getJacobianCol(int deriv_id) const override;
  [[nodiscard]] int getJacobianColsNbr() const override;
  [[nodiscard]] int getLegacyJacobianCol(int deriv_id) const override;

  /* These methods substitute both aggregate and heterogeneous variables with leads/lags
   beyond the standard range (leads >= 2, lags >= 2 for endos; any lead/lag for exos).
   For lead substitution, aggregate variables are treated as known and the
   substitution operates in deterministic mode. As for heterogeneous variables,
   a substitution in stochastic mode is necessary to properly handle the
   expectation operator */
  //! Transforms the model by removing all leads on aggregate and het endos >= 2
  void substituteEndoLeadGreaterThanTwo(DynamicModel& dynamic_model);

  //! Transforms the model by removing all lags >= 2 on aggregate and het endos
  void substituteEndoLagGreaterThanTwo(DynamicModel& dynamic_model);

  //! Transforms the model by removing all leads on aggregate and het exos
  /*! Note that this can create new lags on endos and exos */
  void substituteExoLead(DynamicModel& dynamic_model);

  //! Transforms the model by removing all lags on aggregate and het exos
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

  // Auxiliary equations created by nonlinear expectation substitution
  vector<BinaryOpNode*> het_nonlinear_expectation_aux_equations;

  // Information about MCP multipliers for output generation
  struct MCPMultiplierInfo
  {
    const int multiplier_symb_id; // Symbol ID of the multiplier (MULT_L_* or MULT_U_*)
    const int bound_var_symb_id;  // Symbol ID of the bound variable (e.g., 'a' in a >= 0)
    const expr_t bound_expr;      // The bound expression (lb or ub)
    const bool is_lower_bound;    // true for lower bound (>=), false for upper bound (<=)
    const expr_t
        original_residual; // The original equation's residual (LHS - RHS) before MCP transformation
  };
  vector<MCPMultiplierInfo> mcp_multiplier_info;

  // Allocates the derivation IDs for all endogenous variables for this heterogeneity dimension
  void computeDerivIDs();

  // Writes the file for setting heterogeneous auxiliary variables
  void writeSetHetAuxiliaryVariablesFile(const string& basename) const;
};

#endif
