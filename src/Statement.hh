/*
 * Copyright © 2006-2023 Dynare Team
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

#ifndef _STATEMENT_HH
#define _STATEMENT_HH

#include <ostream>
#include <string>
#include <map>
#include <set>
#include <optional>
#include <variant>

#include "SymbolList.hh"
#include "WarningConsolidation.hh"

struct ModFileStructure
{
  //! Whether check is present
  bool check_present{false};
  //! Whether steady is present
  bool steady_present{false};
  //! Whether a perfect_foresight_solver/simul statement is present
  bool perfect_foresight_solver_present{false};
  //! Whether a perfect_foresight_with_expectation_errors_solver statement is present
  bool perfect_foresight_with_expectation_errors_solver_present{false};
  //! Whether a stoch_simul statement is present
  bool stoch_simul_present{false};
  //! Whether an estimation statement is present
  bool estimation_present{false};
  //! Whether an osr statement is present
  bool osr_present{false};
  //! Whether an osr params statement is present
  bool osr_params_present{false};
  //! Whether an optim weight statement is present
  bool optim_weights_present{false};
  //! Whether a ramsey_model statement is present
  bool ramsey_model_present{false};
  //! Whether a ramsey_policy statement is present
  bool ramsey_policy_present{false};
  //! Whether a discretionary_objective statement is present
  bool discretionary_policy_present{false};
  //! Whether a planner_objective statement is present
  bool planner_objective_present{false};
  //! Whether an extended_path statement is present
  bool extended_path_present{false};
  //! The value of the "order" option of stoch_simul, estimation, osr, ramsey_policy
  //! Derivation order
  /*! First initialized to zero. If user sets order option somewhere in the MOD file, it will be equal to the maximum of order options. Otherwise will default to 2 */
  int order_option{0};
  //! Whether a bvar_density, bvar_forecast, sbvar, ms_sbvar statement is present
  bool bvar_present{false};
  //! Whether an svar_identification statement is present
  bool svar_identification_present{false};
  //! Whether an identification statement is present or the identification option of dynare_sensitivity statement is equal to one
  bool identification_present{false};
  //! The maximum of the “order” option in identification statements
  int identification_order{0};
  //! Whether a sensitivity statement is present
  bool sensitivity_present{false};
  //! Whether the option analytic_derivation is given to estimation
  bool estimation_analytic_derivation{false};
  //! Whether the option partial_information is given to stoch_simul/estimation/osr/ramsey_policy
  bool partial_information{false};
  //! Whether the "k_order_solver" option is used (explictly, or implicitly if order >= 3)
  bool k_order_solver{false};
  //! Whether an method_of_moments statement is present
  bool mom_estimation_present{false};
  //! Whether an GMM-option is present
  bool GMM_present{false};
  //! Whether an analytic_standard_errors-option is present
  bool analytic_standard_errors_present{false};
  //! Whether an analytic_jacobian-option is present
  bool analytic_jacobian_present{false};
  //! The maximum of the “order” option in method_of_moments statements
  int mom_order{0};
  //! Whether there is a calibrated measurement error
  bool calibrated_measurement_errors{false};
  //! Whether dsge_prior_weight was initialized as a parameter
  bool dsge_prior_weight_initialized;
  //! Whether dsge_prior_weight is in the estimated_params block
  bool dsge_prior_weight_in_estimated_params{false};
  //! Whether there is a dsge_var, with calibrated prior weight
  string dsge_var_calibrated;
  //! Whether there is a dsge_var, with prior weight that must be estimated
  bool dsge_var_estimated{false};
  //! Whether there is a bayesian_irf option passed to the estimation statement
  bool bayesian_irf_present{false};
  //! Whether there is a data statement present
  bool estimation_data_statement_present{false};
  //! Last chain number for Markov Switching statement2
  int last_markov_switching_chain{0};
  //! Whether a calib_smoother statement is present
  bool calib_smoother_present{false};
  //! Whether there is an estimated_params_init with use_calibration
  bool estim_params_use_calib{false};
  //! Set of parameters used within shocks blocks, inside the expressions
  //! defining the values of covariances (stored as symbol ids)
  set<int> parameters_within_shocks_values;
  //! Set of estimated parameters (stored as symbol ids)
  set<int> estimated_parameters;
  //! Whether there is a prior statement present
  bool prior_statement_present{false};
  //! Whether there is a std prior statement present
  bool std_prior_statement_present{false};
  //! Whether there is a corr prior statement present
  bool corr_prior_statement_present{false};
  //! Whether there is a options statement present
  bool options_statement_present{false};
  //! Whether there is a std options statement present
  bool std_options_statement_present{false};
  //! Whether there is a corr options statement present
  bool corr_options_statement_present{false};
  //! Whether a Markov Switching DSGE is present
  bool ms_dsge_present{false};
  //! Stores the original number of equations in the model_block
  int orig_eq_nbr{0};
  //! Stores the number of equations added to the Ramsey model
  int ramsey_eq_nbr{0};
  /* The number of endogenous variables in the model present just before adding
     the Lagrange multipliers and computing the Ramsey FOC; it is by
     construction equal to the number of equations that will be added by the
     process of computing the FOCs */
  int ramsey_orig_endo_nbr {0};
  //! Whether there was a steady_state_model block
  bool steady_state_model_present{false};
  //! Whether there is a write_latex_steady_state_model statement present
  bool write_latex_steady_state_model_present{false};
  //! Pac growth and discount
  set<int> pac_params;
  //! Instruments if ramsey_model, ramsey_policy or discretionary_policy is present
  SymbolList instruments;
  /* Whether any of shock_decomposition, realtime_shock_decomposition and
     initial_condition_decomposition has the “with_epilogue” option */
  bool with_epilogue_option{false};
  /* Lists symbol IDs of parameters that appear in a “planner_discount” option.
     See dynare#1173 for more details. */
  set<int> parameters_in_planner_discount;
  // Whether a shocks(surprise) block appears
  bool shocks_surprise_present{false};
  // Whether a shocks(learnt_in=…) block appears
  bool shocks_learnt_in_present{false};
  // Whether an endval(learnt_in=…) block appears
  bool endval_learnt_in_present{false};
  // Whether an occbin_constraints block appears
  bool occbin_constraints_present{false};
  // Whether a ramsey_constraints block appears
  bool ramsey_constraints_present{false};
};

class Statement
{
public:
  Statement() = default;
  virtual ~Statement() = default;

  Statement(const Statement &) = delete;
  Statement &operator=(const Statement &) = delete;

  //! Do some internal check, and fill the ModFileStructure class
  /*! Don't forget to update ComputingTasks.hh, Shocks.hh and
    NumericalInitialization.hh if you modify the signature of this
    method. Otherwise the default implementation (i.e. a no-op) will apply and
    some checks won't be run. */
  virtual void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings);
  virtual void computingPass(const ModFileStructure &mod_file_struct);
  //! Write Matlab output code
  /*!
    \param output is the output stream of the main matlab file
    \param basename is the name of the modfile (without extension) which can be used to build auxiliary files
  */
  virtual void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const = 0;
  virtual void writeJsonOutput(ostream &output) const = 0;
};

class NativeStatement : public Statement
{
private:
  const string native_statement;
public:
  explicit NativeStatement(string native_statement_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class VerbatimStatement : public Statement
{
private:
  const string verbatim_statement;
public:
  explicit VerbatimStatement(string verbatim_statement_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

/* Stores a list of named options with their values.
   The values are stored using an std::variant; see the “options” data member
   for the list of available types. */
class OptionsList
{
public:
  // Some types to lift ambiguities
  struct NumVal : string
  {
  };
  struct StringVal : string
  {
  };
  struct DateVal : string
  {
  };
  struct SymbolListVal : SymbolList
  {
    /* This one is needed because vector<string> implicitly converts to
       SymbolList. Otherwise adding a vector<string> to the variant would add a
       SymbolList, which is probably not the intended meaning. */
  };
  struct VecStrVal : vector<string>
  {
  };
  struct VecCellStrVal : vector<string>
  {
  };
  struct VecValueVal : vector<string>
  {
  };

  bool empty() const;
  void clear();
  // Whether there is an option with that name that has been given a value
  bool contains(const string &name) const;
  // Erase the option with that name
  void erase(const string &name);

  /* Declares an option with a name and value. Overwrite any previous value for
     that name. */
  template<class T>
  void
  set(string name, T &&val)
  {
    options.insert_or_assign(move(name), forward<T>(val));
  }

  struct UnknownOptionException
  {
    const string name;
    UnknownOptionException(string name_arg) : name{move(name_arg)}
    {
    }
  };

  /* Retrieves the value of the option with that name.
     Throws UnknownOptionException if there is no option with that name.
     Throws bad_variant_access if the option has a value of a different type. */
  template<class T>
  T
  get(const string &name) const
  {
    auto it = options.find(name);
    if (it != options.end())
      return std::get<T>(it->second);
    else
      throw UnknownOptionException{name};
  }

  /* Retrieves the value of the option with that name.
     Returns nullopt if there is no option with that name.
     Throws bad_variant_access if the option has a value of a different type. */
  template<class T>
  optional<T>
  get_if(const string &name) const
  {
    auto it = options.find(name);
    if (it != options.end())
      return std::get<T>(it->second);
    else
      return nullopt;
  }

  /* Applies a variant visitor to the value of the option with that name.
     Throws UnknownOptionException if there is no option with that name. */
  template<class Visitor>
  decltype(auto)
  visit(const string &name, Visitor &&vis) const
  {
    auto it = options.find(name);
    if (it != options.end())
      return std::visit(forward<Visitor>(vis), it->second);
    else
      throw UnknownOptionException{name};
  }

  void writeOutput(ostream &output) const;
  void writeOutput(ostream &output, const string &option_group) const;
  void writeJsonOutput(ostream &output) const;

private:
  // pair<string, string> corresponds to a pair of numerical values
  // vector<vector<string>> corresponds to a vector of vectors of numerical values
  map<string, variant<NumVal, pair<string, string>, StringVal, DateVal, SymbolListVal, vector<int>,
                      VecStrVal, VecCellStrVal, VecValueVal, vector<vector<string>>>> options;
  void writeOutputCommon(ostream &output, const string &option_group) const;
  // Helper constant for visitors
  template<class> static constexpr bool always_false_v {false};
};

#endif // ! _STATEMENT_HH
