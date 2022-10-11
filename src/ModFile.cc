/*
 * Copyright © 2006-2022 Dynare Team
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
#include <fstream>
#include <typeinfo>
#include <cassert>
#include <random>

#include <filesystem>

#include "ModFile.hh"
#include "ConfigFile.hh"
#include "ComputingTasks.hh"
#include "Shocks.hh"

ModFile::ModFile(WarningConsolidation &warnings_arg)
  : var_model_table{symbol_table},
    trend_component_model_table{symbol_table},
    var_expectation_model_table{symbol_table},
    pac_model_table{symbol_table},
    expressions_tree{symbol_table, num_constants, external_functions_table},
    original_model{symbol_table, num_constants, external_functions_table,
                   trend_component_model_table, var_model_table},
    dynamic_model{symbol_table, num_constants, external_functions_table,
                  trend_component_model_table, var_model_table},
    trend_dynamic_model{symbol_table, num_constants, external_functions_table,
                        trend_component_model_table, var_model_table},
    orig_ramsey_dynamic_model{symbol_table, num_constants, external_functions_table,
                              trend_component_model_table, var_model_table},
    epilogue{symbol_table, num_constants, external_functions_table,
             trend_component_model_table, var_model_table},
    static_model{symbol_table, num_constants, external_functions_table},
    steady_state_model{symbol_table, num_constants, external_functions_table, static_model},
    warnings{warnings_arg}
{
}

void
ModFile::evalAllExpressions(bool warn_uninit)
{
  cout << "Evaluating expressions..." << endl;

  // Loop over all statements, and fill global eval context if relevant
  for (auto &st : statements)
    {
      if (auto ips = dynamic_cast<InitParamStatement *>(st.get()); ips)
        ips->fillEvalContext(global_eval_context);

      if (auto ies = dynamic_cast<InitOrEndValStatement *>(st.get()); ies)
        ies->fillEvalContext(global_eval_context);

      if (auto lpass = dynamic_cast<LoadParamsAndSteadyStateStatement *>(st.get()); lpass)
        lpass->fillEvalContext(global_eval_context);
    }

  // Evaluate model local variables
  dynamic_model.fillEvalContext(global_eval_context);

  // Check if some symbols are not initialized, and give them a zero value then
  for (int id = 0; id <= symbol_table.maxID(); id++)
    if (auto type = symbol_table.getType(id);
        (type == SymbolType::endogenous || type == SymbolType::exogenous || type == SymbolType::exogenousDet
         || type == SymbolType::parameter || type == SymbolType::modelLocalVariable)
        && !global_eval_context.contains(id))
      {
        if (warn_uninit)
          warnings << "WARNING: Can't find a numeric initial value for "
                   << symbol_table.getName(id) << ", using zero" << endl;
        global_eval_context[id] = 0;
      }
}

void
ModFile::addStatement(unique_ptr<Statement> st)
{
  statements.push_back(move(st));
}

void
ModFile::addStatementAtFront(unique_ptr<Statement> st)
{
  statements.insert(statements.begin(), move(st));
}

void
ModFile::checkPass(bool nostrict, bool stochastic)
{
  for (auto &statement : statements)
    statement->checkPass(mod_file_struct, warnings);

  // Check the steady state block
  steady_state_model.checkPass(mod_file_struct, warnings);

  // Check epilogue block
  epilogue.checkPass(mod_file_struct);

  pac_model_table.checkPass(mod_file_struct);

  if (mod_file_struct.write_latex_steady_state_model_present
      && !mod_file_struct.steady_state_model_present)
    {
      cerr << "ERROR: You cannot have a write_latex_steady_state_model statement without a steady_state_model block." << endl;
      exit(EXIT_FAILURE);
    }

  // If order option has not been set, default to 2
  if (!mod_file_struct.order_option)
    mod_file_struct.order_option = 2;

  param_used_with_lead_lag = dynamic_model.ParamUsedWithLeadLag();
  if (param_used_with_lead_lag)
    warnings << "WARNING: A parameter was used with a lead or a lag in the model block" << endl;

  bool stochastic_statement_present = mod_file_struct.stoch_simul_present
    || mod_file_struct.estimation_present
    || mod_file_struct.osr_present
    || mod_file_struct.ramsey_policy_present
    || mod_file_struct.discretionary_policy_present
    || mod_file_struct.calib_smoother_present
    || mod_file_struct.identification_present
    || mod_file_struct.mom_estimation_present
    || mod_file_struct.sensitivity_present
    || stochastic;

  // Allow empty model only when doing a standalone BVAR estimation
  if (dynamic_model.equation_number() == 0
      && (mod_file_struct.check_present
          || mod_file_struct.perfect_foresight_solver_present
          || mod_file_struct.perfect_foresight_with_expectation_errors_solver_present
          || stochastic_statement_present))
    {
      cerr << "ERROR: At least one model equation must be declared!" << endl;
      exit(EXIT_FAILURE);
    }

  if ((mod_file_struct.ramsey_model_present || mod_file_struct.ramsey_policy_present)
      && mod_file_struct.discretionary_policy_present)
    {
      cerr << "ERROR: You cannot use the discretionary_policy command when you use either ramsey_model or ramsey_policy and vice versa" << endl;
      exit(EXIT_FAILURE);
    }

  if (((mod_file_struct.ramsey_model_present || mod_file_struct.discretionary_policy_present)
       && !mod_file_struct.planner_objective_present)
      || (!(mod_file_struct.ramsey_model_present || mod_file_struct.discretionary_policy_present)
          && mod_file_struct.planner_objective_present))
    {
      cerr << "ERROR: A planner_objective statement must be used with a ramsey_model, a ramsey_policy or a discretionary_policy statement and vice versa." << endl;
      exit(EXIT_FAILURE);
    }

  if (mod_file_struct.ramsey_constraints_present && !mod_file_struct.ramsey_model_present && !mod_file_struct.ramsey_policy_present)
    {
      cerr << "ERROR: A ramsey_constraints block requires the presence of a ramsey_model or ramsey_policy statement" << endl;
      exit(EXIT_FAILURE);
    }

  if ((mod_file_struct.osr_present && (!mod_file_struct.osr_params_present || !mod_file_struct.optim_weights_present))
      || ((!mod_file_struct.osr_present || !mod_file_struct.osr_params_present) && mod_file_struct.optim_weights_present)
      || ((!mod_file_struct.osr_present || !mod_file_struct.optim_weights_present) && mod_file_struct.osr_params_present))
    {
      cerr << "ERROR: The osr statement must be used with osr_params and optim_weights." << endl;
      exit(EXIT_FAILURE);
    }

  if ((mod_file_struct.perfect_foresight_solver_present || mod_file_struct.perfect_foresight_with_expectation_errors_solver_present)
      && stochastic_statement_present)
    {
      cerr << "ERROR: A .mod file cannot contain both one of {perfect_foresight_solver, simul, perfect_foresight_with_expectation_errors_solver} and one of {stoch_simul, estimation, osr, ramsey_policy, discretionary_policy}. This is not possible: one cannot mix perfect foresight context with stochastic context in the same file." << endl;
      exit(EXIT_FAILURE);
    }

  if (mod_file_struct.k_order_solver && bytecode)
    {
      cerr << "ERROR: 'k_order_solver' (which is implicit if order >= 3), is not yet compatible with 'bytecode'." << endl;
      exit(EXIT_FAILURE);
    }

  if (use_dll && bytecode)
    {
      cerr << "ERROR: In 'model' block, 'use_dll' option is not compatible with 'bytecode'" << endl;
      exit(EXIT_FAILURE);
    }

  if ((stochastic_statement_present || mod_file_struct.check_present || mod_file_struct.steady_present) && no_static)
    {
      cerr << "ERROR: no_static option is incompatible with stoch_simul, estimation, osr, ramsey_policy, discretionary_policy, steady and check commands" << endl;
      exit(EXIT_FAILURE);
    }

  if (mod_file_struct.dsge_var_estimated && !mod_file_struct.dsge_prior_weight_in_estimated_params)
    {
      cerr << "ERROR: When estimating a DSGE-VAR model and estimating the weight of the prior, dsge_prior_weight must "
           << "be referenced in the estimated_params block." << endl;
      exit(EXIT_FAILURE);
    }

  if (symbol_table.exists("dsge_prior_weight"))
    {
      if (symbol_table.getType("dsge_prior_weight") != SymbolType::parameter)
        {
          cerr << "ERROR: dsge_prior_weight may only be used as a parameter." << endl;
          exit(EXIT_FAILURE);
        }
      else
        warnings << "WARNING: When estimating a DSGE-Var, declaring dsge_prior_weight as a "
                 << "parameter is deprecated. The preferred method is to do this via "
                 << "the dsge_var option in the estimation statement." << endl;

      if (mod_file_struct.dsge_var_estimated || !mod_file_struct.dsge_var_calibrated.empty())
        {
          cerr << "ERROR: dsge_prior_weight can either be declared as a parameter (deprecated) or via the dsge_var option "
               << "to the estimation statement (preferred), but not both." << endl;
          exit(EXIT_FAILURE);
        }

      if (!mod_file_struct.dsge_prior_weight_initialized && !mod_file_struct.dsge_prior_weight_in_estimated_params)
        {
          cerr << "ERROR: If dsge_prior_weight is declared as a parameter, it must either be initialized or placed in the "
               << "estimated_params block." << endl;
          exit(EXIT_FAILURE);
        }

      if (mod_file_struct.dsge_prior_weight_initialized && mod_file_struct.dsge_prior_weight_in_estimated_params)
        {
          cerr << "ERROR: dsge_prior_weight cannot be both initialized and estimated." << endl;
          exit(EXIT_FAILURE);
        }
    }

  if (mod_file_struct.dsge_prior_weight_in_estimated_params)
    if (!mod_file_struct.dsge_var_estimated && !mod_file_struct.dsge_var_calibrated.empty())
      {
        cerr << "ERROR: If dsge_prior_weight is in the estimated_params block, the prior weight cannot be calibrated "
             << "via the dsge_var option in the estimation statement." << endl;
        exit(EXIT_FAILURE);
      }
    else if (!mod_file_struct.dsge_var_estimated && !symbol_table.exists("dsge_prior_weight"))
      {
        cerr << "ERROR: If dsge_prior_weight is in the estimated_params block, it must either be declared as a parameter "
             << "(deprecated) or the dsge_var option must be passed to the estimation statement (preferred)." << endl;
        exit(EXIT_FAILURE);
      }

  if (dynamic_model.staticOnlyEquationsNbr() != dynamic_model.dynamicOnlyEquationsNbr())
    {
      cerr << "ERROR: the number of equations marked [static] must be equal to the number of equations marked [dynamic]" << endl;
      exit(EXIT_FAILURE);
    }

  if (dynamic_model.staticOnlyEquationsNbr() > 0
      && (mod_file_struct.ramsey_model_present || mod_file_struct.discretionary_policy_present))
    {
      cerr << "ERROR: marking equations as [static] or [dynamic] is not possible with ramsey_model, ramsey_policy or discretionary_policy" << endl;
      exit(EXIT_FAILURE);
    }

  if (stochastic_statement_present
      && (dynamic_model.isUnaryOpUsed(UnaryOpcode::sign)
          || dynamic_model.isUnaryOpUsed(UnaryOpcode::abs)
          || dynamic_model.isBinaryOpUsed(BinaryOpcode::max)
          || dynamic_model.isBinaryOpUsed(BinaryOpcode::min)
          || dynamic_model.isBinaryOpUsed(BinaryOpcode::greater)
          || dynamic_model.isBinaryOpUsed(BinaryOpcode::less)
          || dynamic_model.isBinaryOpUsed(BinaryOpcode::greaterEqual)
          || dynamic_model.isBinaryOpUsed(BinaryOpcode::lessEqual)
          || dynamic_model.isBinaryOpUsed(BinaryOpcode::equalEqual)
          || dynamic_model.isBinaryOpUsed(BinaryOpcode::different)))
    warnings << R"(WARNING: you are using a function (max, min, abs, sign) or an operator (<, >, <=, >=, ==, !=) which is unsuitable for a stochastic context; see the reference manual, section about "Expressions", for more details.)" << endl;

  if (linear
      && (dynamic_model.isUnaryOpUsedOnType(SymbolType::endogenous, UnaryOpcode::sign)
          || dynamic_model.isUnaryOpUsedOnType(SymbolType::endogenous, UnaryOpcode::abs)
          || dynamic_model.isBinaryOpUsedOnType(SymbolType::endogenous, BinaryOpcode::max)
          || dynamic_model.isBinaryOpUsedOnType(SymbolType::endogenous, BinaryOpcode::min)
          || dynamic_model.isBinaryOpUsedOnType(SymbolType::endogenous, BinaryOpcode::greater)
          || dynamic_model.isBinaryOpUsedOnType(SymbolType::endogenous, BinaryOpcode::less)
          || dynamic_model.isBinaryOpUsedOnType(SymbolType::endogenous, BinaryOpcode::greaterEqual)
          || dynamic_model.isBinaryOpUsedOnType(SymbolType::endogenous, BinaryOpcode::lessEqual)
          || dynamic_model.isBinaryOpUsedOnType(SymbolType::endogenous, BinaryOpcode::equalEqual)
          || dynamic_model.isBinaryOpUsedOnType(SymbolType::endogenous, BinaryOpcode::different)))
    {
      cerr << "ERROR: you have declared your model 'linear' but you are using a function "
           << "(max, min, abs, sign) or an operator (<, >, <=, >=, ==, !=) on an "
           << "endogenous variable." << endl;
      exit(EXIT_FAILURE);
    }

  if (linear
      && !mod_file_struct.perfect_foresight_solver_present
      && !mod_file_struct.perfect_foresight_with_expectation_errors_solver_present
      && (dynamic_model.isUnaryOpUsedOnType(SymbolType::exogenous, UnaryOpcode::sign)
          || dynamic_model.isUnaryOpUsedOnType(SymbolType::exogenous, UnaryOpcode::abs)
          || dynamic_model.isBinaryOpUsedOnType(SymbolType::exogenous, BinaryOpcode::max)
          || dynamic_model.isBinaryOpUsedOnType(SymbolType::exogenous, BinaryOpcode::min)
          || dynamic_model.isBinaryOpUsedOnType(SymbolType::exogenous, BinaryOpcode::greater)
          || dynamic_model.isBinaryOpUsedOnType(SymbolType::exogenous, BinaryOpcode::less)
          || dynamic_model.isBinaryOpUsedOnType(SymbolType::exogenous, BinaryOpcode::greaterEqual)
          || dynamic_model.isBinaryOpUsedOnType(SymbolType::exogenous, BinaryOpcode::lessEqual)
          || dynamic_model.isBinaryOpUsedOnType(SymbolType::exogenous, BinaryOpcode::equalEqual)
          || dynamic_model.isBinaryOpUsedOnType(SymbolType::exogenous, BinaryOpcode::different)))
    {
      cerr << "ERROR: you have declared your model 'linear' but you are using a function "
           << "(max, min, abs, sign) or an operator (<, >, <=, >=, ==, !=) on an "
           << "exogenous variable in a non-perfect-foresight context." << endl;
      exit(EXIT_FAILURE);
    }

  // Test if some estimated parameters are used within the values of shocks
  // statements (see issue #469)
  set<int> parameters_intersect;
  set_intersection(mod_file_struct.parameters_within_shocks_values.begin(),
                   mod_file_struct.parameters_within_shocks_values.end(),
                   mod_file_struct.estimated_parameters.begin(),
                   mod_file_struct.estimated_parameters.end(),
                   inserter(parameters_intersect, parameters_intersect.begin()));
  if (parameters_intersect.size() > 0)
    {
      cerr << "ERROR: some estimated parameters (";
      for (bool printed_something{false};
           int symb_id : parameters_intersect)
        {
          if (exchange(printed_something, true))
            cerr << ", ";
          cerr << symbol_table.getName(symb_id);
        }
      cerr << ") also appear in the expressions defining the variance/covariance matrix of shocks; this is not allowed." << endl;
      exit(EXIT_FAILURE);
    }

  // Check if some exogenous is not used in the model block, Issue #841
  set<int> unusedExo0 = dynamic_model.findUnusedExogenous();
  set<int> unusedExo;
  set_difference(unusedExo0.begin(), unusedExo0.end(),
                 mod_file_struct.pac_params.begin(), mod_file_struct.pac_params.end(),
                 inserter(unusedExo, unusedExo.begin()));
  if (unusedExo.size() > 0)
    {
      ostringstream unused_exos;
      for (int it : unusedExo)
        unused_exos << symbol_table.getName(it) << " ";

      if (nostrict)
        warnings << "WARNING: " << unused_exos.str()
                 << "not used in model block, removed by nostrict command-line option" << endl;
      else
        {
          cerr << "ERROR: " << unused_exos.str() << "not used in model block. To bypass this error, use the `nostrict` option. This may lead to crashes or unexpected behavior." << endl;
          exit(EXIT_FAILURE);
        }
    }

  // See dynare#1726
  if ((stochastic_statement_present || mod_file_struct.check_present) && dynamic_model.mfs > 0)
    {
      cerr << "ERROR: mfs > 0 is incompatible with check, stoch_simul, estimation, osr, ramsey_policy, discretionary_policy, calib_smoother, identification, methods_of_moments and sensitivity commands" << endl;
      exit(EXIT_FAILURE);
    }
}

void
ModFile::transformPass(bool nostrict, bool stochastic, bool compute_xrefs, bool transform_unary_ops,
                       const string &exclude_eqs, const string &include_eqs)
{
  /* Save the original model (must be done before any model transformations by preprocessor)
     — except predetermined variables (which must be handled before the call to
       setLeadsLagsOrig(), see #47, and also before equation simplification,
       since the latter uses leads/lags, see #83)
     — except substituting out variables which we know are constant (they
       appear in an equation of the form: X = constant)
     — except adl operators which we always want expanded
     — except diff operators with a lead which have been expanded by
       DataTree:AddDiff()
  */
  dynamic_model.includeExcludeEquations(exclude_eqs, true);
  dynamic_model.includeExcludeEquations(include_eqs, false);
  if (symbol_table.predeterminedNbr() > 0)
    dynamic_model.transformPredeterminedVariables();
  dynamic_model.simplifyEquations();
  dynamic_model.substituteAdl();
  dynamic_model.setLeadsLagsOrig();
  original_model = dynamic_model;
  dynamic_model.expandEqTags();

  // Replace all model-local variables by their expression
  dynamic_model.substituteModelLocalVariables();

  // Check that all declared endogenous are used in equations
  set<int> unusedEndogs = dynamic_model.findUnusedEndogenous();
  bool unusedEndogsIsErr = !nostrict && !mod_file_struct.bvar_present && unusedEndogs.size();
  for (int unusedEndog : unusedEndogs)
    if (nostrict)
      {
        symbol_table.changeType(unusedEndog, SymbolType::unusedEndogenous);
        warnings << "WARNING: '" << symbol_table.getName(unusedEndog)
                 << "' not used in model block, removed by nostrict command-line option" << endl;
      }
    else if (unusedEndogsIsErr)
      cerr << "Error: " << symbol_table.getName(unusedEndog) << " not used in the model block"<< endl;

  if (unusedEndogsIsErr)
    exit(EXIT_FAILURE);

  /* Get the list of equations in which to scan for and substitute unary ops:
     – equations which are part of VARs and Trend Component Models
     – PAC equations (those with a pac_expectation operator) */
  set<string> var_tcm_eqtags;
  for (const auto &[name, tags] : trend_component_model_table.getEqTags())
    for (auto &tag : tags)
      var_tcm_eqtags.insert(tag);
  for (const auto &[name, tags] : var_model_table.getEqTags())
    for (auto &tag : tags)
      var_tcm_eqtags.insert(tag);

  set<int> unary_ops_eqs = dynamic_model.getEquationNumbersFromTags(var_tcm_eqtags);
  unary_ops_eqs.merge(dynamic_model.findPacExpectationEquationNumbers());

  // Check that no variable in VAR/TCM/PAC equations was declared with “var(log)”
  dynamic_model.checkNoWithLogTransform(unary_ops_eqs);

  // Create auxiliary variables and equations for unary ops
  lag_equivalence_table_t unary_ops_nodes;
  ExprNode::subst_table_t unary_ops_subst_table;
  if (transform_unary_ops)
    tie(unary_ops_nodes, unary_ops_subst_table) = dynamic_model.substituteUnaryOps(var_expectation_model_table, pac_model_table);
  else
    // substitute only those unary ops that appear in VAR, TCM and PAC model equations
    tie(unary_ops_nodes, unary_ops_subst_table) = dynamic_model.substituteUnaryOps(unary_ops_eqs, var_expectation_model_table, pac_model_table);

  // Create auxiliary variable and equations for Diff operators
  auto [diff_nodes, diff_subst_table] = dynamic_model.substituteDiff(var_expectation_model_table, pac_model_table);

  // Fill trend component and VAR model tables
  dynamic_model.fillTrendComponentModelTable();
  original_model.fillTrendComponentModelTableFromOrigModel();
  dynamic_model.fillTrendComponentModelTableAREC(diff_subst_table);
  dynamic_model.fillVarModelTable();
  original_model.fillVarModelTableFromOrigModel();

  // VAR expectation models
  var_expectation_model_table.transformPass(diff_subst_table, dynamic_model, var_model_table,
                                            trend_component_model_table);

  // PAC model
  pac_model_table.transformPass(unary_ops_nodes, unary_ops_subst_table,
                                diff_nodes, diff_subst_table,
                                dynamic_model, var_model_table,
                                trend_component_model_table);

  // Create auxiliary vars for Expectation operator
  dynamic_model.substituteExpectation(mod_file_struct.partial_information);

  if (nonstationary_variables)
    {
      dynamic_model.detrendEquations();
      trend_dynamic_model = dynamic_model;
      dynamic_model.removeTrendVariableFromEquations();
      const auto &trend_symbols = dynamic_model.getTrendSymbolsMap();
      const auto &nonstationary_symbols = dynamic_model.getNonstationarySymbolsMap();
      epilogue.detrend(trend_symbols, nonstationary_symbols);
    }

  epilogue.toStatic();

  mod_file_struct.orig_eq_nbr = dynamic_model.equation_number();
  if (mod_file_struct.ramsey_model_present)
    {
      PlannerObjectiveStatement *pos = nullptr;
      for (auto &statement : statements)
        if (auto pos2 = dynamic_cast<PlannerObjectiveStatement *>(statement.get()); pos2)
          if (pos)
            {
              cerr << "ERROR: there can only be one planner_objective statement" << endl;
              exit(EXIT_FAILURE);
            }
          else
            pos = pos2;
      assert(pos);
      const PlannerObjective &planner_objective = pos->getPlannerObjective();

      /*
        clone the model then clone the new equations back to the original because
        we have to call computeDerivIDs (in computeRamseyPolicyFOCs and computingPass)
      */
      if (linear)
        orig_ramsey_dynamic_model = dynamic_model;
      DynamicModel ramsey_FOC_equations_dynamic_model {symbol_table, num_constants, external_functions_table, trend_component_model_table, var_model_table};
      ramsey_FOC_equations_dynamic_model = dynamic_model;
      ramsey_FOC_equations_dynamic_model.computeRamseyPolicyFOCs(planner_objective);
      ramsey_FOC_equations_dynamic_model.replaceMyEquations(dynamic_model);
      mod_file_struct.ramsey_eq_nbr = dynamic_model.equation_number() - mod_file_struct.orig_eq_nbr;
    }

  dynamic_model.createVariableMapping();

  // Must come after detrending of variables and Ramsey policy transformation
  dynamic_model.substituteLogTransform();

  /* Create auxiliary vars for leads and lags greater than 2, on both endos and
     exos. The transformation is not exactly the same on stochastic and
     deterministic models, because there is no need to take into account the
     Jensen inequality on the latter. */
  bool deterministic_model = !(mod_file_struct.stoch_simul_present
                               || mod_file_struct.estimation_present
                               || mod_file_struct.osr_present
                               || mod_file_struct.ramsey_policy_present
                               || mod_file_struct.discretionary_policy_present
                               || mod_file_struct.calib_smoother_present
                               || mod_file_struct.identification_present
                               || mod_file_struct.mom_estimation_present
                               || mod_file_struct.sensitivity_present
                               || stochastic);
  dynamic_model.substituteEndoLeadGreaterThanTwo(deterministic_model);
  dynamic_model.substituteExoLead(deterministic_model);
  dynamic_model.substituteEndoLagGreaterThanTwo(deterministic_model);
  dynamic_model.substituteExoLag(deterministic_model);

  dynamic_model.updateVarAndTrendModel();

  if (differentiate_forward_vars)
    dynamic_model.differentiateForwardVars(differentiate_forward_vars_subset);

  if (mod_file_struct.dsge_var_estimated || !mod_file_struct.dsge_var_calibrated.empty())
    try
      {
        int sid = symbol_table.addSymbol("dsge_prior_weight", SymbolType::parameter);
        if (!mod_file_struct.dsge_var_calibrated.empty())
          addStatementAtFront(make_unique<InitParamStatement>(sid,
                                                              expressions_tree.AddNonNegativeConstant(mod_file_struct.dsge_var_calibrated),
                                                              symbol_table));
      }
    catch (SymbolTable::AlreadyDeclaredException &e)
      {
        cerr << "ERROR: dsge_prior_weight should not be declared as a model variable / parameter "
             << "when the dsge_var option is passed to the estimation statement." << endl;
        exit(EXIT_FAILURE);
      }

  dynamic_model.reorderAuxiliaryEquations();

  // Freeze the symbol table
  symbol_table.freeze();

  if (compute_xrefs)
    dynamic_model.computeXrefs();

  /*
    Enforce the same number of equations and endogenous, except in three cases:
    - ramsey_model, ramsey_policy or discretionary_policy is used
    - a BVAR command is used and there is no equation (standalone BVAR estimation)
    - nostrict option is passed and there are more endogs than equations (dealt with before freeze)
  */
  if (!(mod_file_struct.ramsey_model_present || mod_file_struct.discretionary_policy_present)
      && !(mod_file_struct.bvar_present && dynamic_model.equation_number() == 0)
      && (dynamic_model.equation_number() != symbol_table.endo_nbr()))
    {
      cerr << "ERROR: There are " << dynamic_model.equation_number() << " equations but " << symbol_table.endo_nbr() << " endogenous variables!" << endl;
      exit(EXIT_FAILURE);
    }

  if (symbol_table.exo_det_nbr() > 0
      && (mod_file_struct.perfect_foresight_solver_present || mod_file_struct.perfect_foresight_with_expectation_errors_solver_present))
    {
      cerr << "ERROR: A .mod file cannot contain both one of {perfect_foresight_solver, simul, perfect_foresight_with_expectation_errors_solver} and varexo_det declaration (all exogenous variables are deterministic in this case)" << endl;
      exit(EXIT_FAILURE);
    }

  if (mod_file_struct.ramsey_policy_present && symbol_table.exo_det_nbr() > 0)
    {
      cerr << "ERROR: ramsey_policy is incompatible with deterministic exogenous variables" << endl;
      exit(EXIT_FAILURE);
    }

  if (mod_file_struct.identification_present && symbol_table.exo_det_nbr() > 0)
    {
      cerr << "ERROR: identification is incompatible with deterministic exogenous variables" << endl;
      exit(EXIT_FAILURE);
    }

  if (mod_file_struct.occbin_constraints_present
      && (mod_file_struct.osr_present || mod_file_struct.mom_estimation_present
          || mod_file_struct.ramsey_model_present || mod_file_struct.ramsey_policy_present
          || mod_file_struct.discretionary_policy_present || mod_file_struct.extended_path_present
          || mod_file_struct.identification_present || mod_file_struct.sensitivity_present))
    {
      cerr << "ERROR: the 'occbin_constraints' block is not compatible with commands other than 'estimation', 'stoch_simul', and 'calib_smoother'." << endl;
      exit(EXIT_FAILURE);
    }

  if (mod_file_struct.shocks_surprise_present && !mod_file_struct.occbin_constraints_present)
    {
      cerr << "ERROR: the 'shocks(surprise)' block can only be used in conjunction with the 'occbin_constraints' block." << endl;
      exit(EXIT_FAILURE);
    }

  if (mod_file_struct.shocks_learnt_in_present && !mod_file_struct.perfect_foresight_with_expectation_errors_solver_present)
    {
      cerr << "ERROR: the 'shocks(learnt_in=…)' block can only be used in conjunction with the 'perfect_foresight_with_expectation_errors_solver' command." << endl;
      exit(EXIT_FAILURE);
    }

  if (mod_file_struct.endval_learnt_in_present && !mod_file_struct.perfect_foresight_with_expectation_errors_solver_present)
    {
      cerr << "ERROR: the 'endval(learnt_in=…)' block can only be used in conjunction with the 'perfect_foresight_with_expectation_errors_solver' command." << endl;
      exit(EXIT_FAILURE);
    }

  if (!mod_file_struct.ramsey_model_present)
    cout << "Found " << dynamic_model.equation_number() << " equation(s)." << endl;
  else
    {
      cout << "Found " << mod_file_struct.orig_eq_nbr  << " equation(s)." << endl;
      cout << "Found " << dynamic_model.equation_number() << " FOC equation(s) for Ramsey Problem." << endl;
    }

  if (symbol_table.exists("dsge_prior_weight"))
    if (mod_file_struct.bayesian_irf_present)
      {
        if (symbol_table.exo_nbr() != symbol_table.observedVariablesNbr())
          {
            cerr << "ERROR: When estimating a DSGE-Var and the bayesian_irf option is passed to the estimation "
                 << "statement, the number of shocks must equal the number of observed variables." << endl;
            exit(EXIT_FAILURE);
          }
      }
    else
      if (symbol_table.exo_nbr() < symbol_table.observedVariablesNbr())
        {
          cerr << "ERROR: When estimating a DSGE-Var, the number of shocks must be "
               << "greater than or equal to the number of observed variables." << endl;
          exit(EXIT_FAILURE);
        }
}

void
ModFile::computingPass(bool no_tmp_terms, OutputType output, int params_derivs_order)
{
  // Mod file may have no equation (for example in a standalone BVAR estimation)
  if (dynamic_model.equation_number() > 0)
    {
      if (nonstationary_variables)
        trend_dynamic_model.runTrendTest(global_eval_context);

      // Compute static model and its derivatives
      static_model = static_cast<StaticModel>(dynamic_model);
      if (!no_static)
        {
          if (mod_file_struct.stoch_simul_present
              || mod_file_struct.estimation_present || mod_file_struct.osr_present
              || mod_file_struct.ramsey_model_present || mod_file_struct.identification_present
              || mod_file_struct.calib_smoother_present || mod_file_struct.mom_estimation_present)
            static_model.set_cutoff_to_zero();

          int derivsOrder = 1;
          int paramsDerivsOrder = 0;
          if (mod_file_struct.identification_present || mod_file_struct.estimation_analytic_derivation)
            derivsOrder = 2;

          if (mod_file_struct.identification_present 
              || mod_file_struct.estimation_analytic_derivation
              || (mod_file_struct.GMM_present && (mod_file_struct.analytic_standard_errors_present || mod_file_struct.analytic_jacobian_present)))
            paramsDerivsOrder = params_derivs_order;

          static_model.computingPass(derivsOrder, paramsDerivsOrder, global_eval_context, no_tmp_terms, block);
        }
      // Set things to compute for dynamic model
      if (mod_file_struct.perfect_foresight_solver_present
          || mod_file_struct.perfect_foresight_with_expectation_errors_solver_present
          || mod_file_struct.check_present
          || mod_file_struct.stoch_simul_present
          || mod_file_struct.estimation_present || mod_file_struct.osr_present
          || mod_file_struct.ramsey_model_present || mod_file_struct.identification_present
          || mod_file_struct.calib_smoother_present || mod_file_struct.mom_estimation_present)
        {
          if (mod_file_struct.perfect_foresight_solver_present
              || mod_file_struct.perfect_foresight_with_expectation_errors_solver_present)
            {
              int derivsOrder = 1;
              if (output == OutputType::second)
                derivsOrder = 2;
              else if  (output == OutputType::third)
                derivsOrder = 3;
              dynamic_model.computingPass(derivsOrder, 0, global_eval_context, no_tmp_terms, block, use_dll);
            }
          else
            {
              if (mod_file_struct.stoch_simul_present
                  || mod_file_struct.estimation_present || mod_file_struct.osr_present
                  || mod_file_struct.ramsey_model_present || mod_file_struct.identification_present
                  || mod_file_struct.calib_smoother_present || mod_file_struct.mom_estimation_present)
                dynamic_model.set_cutoff_to_zero();
              if (mod_file_struct.order_option < 1)
                {
                  cerr << "ERROR: Incorrect order option..." << endl;
                  exit(EXIT_FAILURE);
                }
              int derivsOrder = max(mod_file_struct.order_option,mod_file_struct.identification_order + 1); // See preprocessor#40
              if (mod_file_struct.GMM_present 
                  && (mod_file_struct.analytic_standard_errors_present || mod_file_struct.analytic_jacobian_present)) //analytic_standard_errors or analytic_jacobian require one order more
                derivsOrder = max(mod_file_struct.order_option,
                                    max(mod_file_struct.identification_order,mod_file_struct.mom_order) + 1); // See preprocessor#40

              if (mod_file_struct.sensitivity_present || linear || output == OutputType::second)
                derivsOrder = max(derivsOrder, 2);
              if (mod_file_struct.estimation_analytic_derivation || output == OutputType::third)
                derivsOrder = max(derivsOrder, 3);
              int paramsDerivsOrder = 0;
              if (mod_file_struct.identification_present 
                  || mod_file_struct.estimation_analytic_derivation
                  || (mod_file_struct.GMM_present && (mod_file_struct.analytic_standard_errors_present || mod_file_struct.analytic_jacobian_present)))
                paramsDerivsOrder = params_derivs_order;
              dynamic_model.computingPass(derivsOrder, paramsDerivsOrder, global_eval_context, no_tmp_terms, block, use_dll);
              if (linear && mod_file_struct.ramsey_model_present)
                orig_ramsey_dynamic_model.computingPass(2, paramsDerivsOrder, global_eval_context, no_tmp_terms, block, use_dll);
            }
        }
      else // No computing task requested, compute derivatives up to 2nd order by default
        dynamic_model.computingPass(2, 0, global_eval_context, no_tmp_terms, block, use_dll);

      /* Check that the model is linear.
         FIXME: this check always passes if derivsOrder = 1, i.e. for a perfect
         foresight model, because the Hessian is not computed in that case. */
      if (linear)
        {
          set<int> eqs = mod_file_struct.ramsey_model_present ?
            orig_ramsey_dynamic_model.getNonZeroHessianEquations() :
            dynamic_model.getNonZeroHessianEquations();

          if (!eqs.empty())
            {
              cerr << "ERROR: If the model is declared linear the second derivatives must be equal to zero." << endl
                   << "       The following equations have non-zero second derivatives:" << endl;
              for (const auto &it : eqs)
                {
                  cerr << "       * Eq # " << it+1;
                  auto tags = dynamic_model.getEquationTags(it);
                  if (auto it2 = tags.find("name"); it2 != tags.end())
                    cerr << " [" << it2->second << "]";
                  cerr << endl;
                }
              exit(EXIT_FAILURE);
            }
        }
    }

  // Those matrices can only be filled here, because we use derivatives
  dynamic_model.fillVarModelTableMatrices();

  for (auto &statement : statements)
    statement->computingPass(mod_file_struct);

  // Compute epilogue derivatives (but silence standard output)
  streambuf *oldcout = cout.rdbuf();
  cout.rdbuf(nullptr);
  epilogue.computingPass(2, 0, global_eval_context, true, false, false);
  cout.rdbuf(oldcout);
}

void
ModFile::writeMOutput(const string &basename, bool clear_all, bool clear_global, bool no_warn,
                      bool console, bool nograph, bool nointeractive, const ConfigFile &config_file,
                      bool check_model_changes, bool minimal_workspace, bool compute_xrefs,
                      const string &mexext,
                      const filesystem::path &matlabroot,
                      const filesystem::path &dynareroot, bool onlymodel, bool gui, bool notime) const
{
  if (basename.empty())
    {
      cerr << "ERROR: Missing file name" << endl;
      exit(EXIT_FAILURE);
    }

  auto plusfolder {DataTree::packageDir(basename)};

  bool hasModelChanged = !dynamic_model.isChecksumMatching(basename) || !check_model_changes;
  if (hasModelChanged)
    {
      // Erase possible remnants of previous runs
      /* Under MATLAB+Windows (but not under Octave nor under GNU/Linux or
         macOS), if we directly remove the "+" subdirectory, then the
         preprocessor is not able to recreate it afterwards (presumably because
         MATLAB maintains some sort of lock on it). The workaround is to rename
         it before deleting it (the renaming must occur in the same directory,
         otherwise it may file if the destination is not on the same
         filesystem). */
      if (exists(plusfolder))
        {
          if (exists(plusfolder / "+objective"))
            {
              // Do it recursively for the +objective folder, created by ramsey_policy
              auto tmp2 = unique_path();
              rename(plusfolder / "+objective", tmp2);
              remove_all(tmp2);
            }

          auto tmp = unique_path();
          rename(plusfolder, tmp);
          remove_all(tmp);
        }
      filesystem::remove_all(basename + "/model/src");
      filesystem::remove_all(basename + "/model/bytecode");
    }

  create_directory(plusfolder);
  filesystem::path fname {plusfolder / "driver.m"};
  ofstream mOutputFile{fname, ios::out | ios::binary};
  if (!mOutputFile.is_open())
    {
      cerr << "ERROR: Can't open file " << fname.string() << " for writing" << endl;
      exit(EXIT_FAILURE);
    }

  mOutputFile << "%" << endl
              << "% Status : main Dynare file" << endl
              << "%" << endl
              << "% Warning : this file is generated automatically by Dynare" << endl
              << "%           from model file (.mod)" << endl << endl;

  if (no_warn)
    mOutputFile << "warning off" << endl; // This will be executed *after* function warning_config()

  if (clear_all)
    mOutputFile << "if isoctave || matlab_ver_less_than('8.6')" << endl
                << "    clear all" << endl
                << "else" << endl
                << "    clearvars -global" << endl
                << "    clear_persistent_variables(fileparts(which('dynare')), false)" << endl
                << "end" << endl;
  else if (clear_global)
    mOutputFile << "clear M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_;" << endl;

  if (!notime)
    mOutputFile << "tic0 = tic;" << endl;

  mOutputFile << "% Define global variables." << endl
              << "global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_" << endl
              << "options_ = [];" << endl
              << "M_.fname = '" << basename << "';" << endl
              << "M_.dynare_version = '" << PACKAGE_VERSION << "';" << endl
              << "oo_.dynare_version = '" << PACKAGE_VERSION << "';" << endl
              << "options_.dynare_version = '" << PACKAGE_VERSION << "';" << endl
              << "%" << endl
              << "% Some global variables initialization" << endl
              << "%" << endl;
  if (!onlymodel)
    config_file.writeHooks(mOutputFile);
  mOutputFile << "global_initialization;" << endl;

  if (minimal_workspace)
    mOutputFile << "options_.minimal_workspace = true;" << endl;

  if (console)
    mOutputFile << "options_.console_mode = true;" << endl
                << "options_.nodisplay = true;" << endl;
  if (nograph)
    mOutputFile << "options_.nograph = true;" << endl;

  if (nointeractive)
    mOutputFile << "options_.nointeractive = true;" << endl;

  if (param_used_with_lead_lag)
    mOutputFile << "M_.parameter_used_with_lead_lag = true;" << endl;

  symbol_table.writeOutput(mOutputFile);

  // Fill the fields related to solve_algo={12,14} if possible.
  mOutputFile << "M_ = setup_solvers(M_);" << endl;

  var_model_table.writeOutput(basename, mOutputFile);
  trend_component_model_table.writeOutput(basename, mOutputFile);
  var_expectation_model_table.writeOutput(mOutputFile);
  pac_model_table.writeOutput(mOutputFile);

  // Initialize M_.Sigma_e, M_.Correlation_matrix, M_.H, and M_.Correlation_matrix_ME
  mOutputFile << "M_.Sigma_e = zeros(" << symbol_table.exo_nbr() << ", "
              << symbol_table.exo_nbr() << ");" << endl
              << "M_.Correlation_matrix = eye(" << symbol_table.exo_nbr() << ", "
              << symbol_table.exo_nbr() << ");" << endl;

  if (mod_file_struct.calibrated_measurement_errors)
    mOutputFile << "M_.H = zeros(" << symbol_table.observedVariablesNbr() << ", "
                << symbol_table.observedVariablesNbr() << ");" << endl
                << "M_.Correlation_matrix_ME = eye(" << symbol_table.observedVariablesNbr() << ", "
                << symbol_table.observedVariablesNbr() << ");" << endl;
  else
    mOutputFile << "M_.H = 0;" << endl
                << "M_.Correlation_matrix_ME = 1;" << endl;

  // May be later modified by a shocks block
  mOutputFile << "M_.sigma_e_is_diagonal = true;" << endl;

  // Initialize M_.det_shocks, M_.surprise_shocks, M_.learnt_shocks, M_.learnt_endval and M_.heteroskedastic_shocks
  mOutputFile << "M_.det_shocks = [];" << endl
              << "M_.surprise_shocks = [];" << endl
              << "M_.learnt_shocks = [];" << endl
              << "M_.learnt_endval = [];" << endl
              << "M_.heteroskedastic_shocks.Qvalue_orig = [];" << endl
              << "M_.heteroskedastic_shocks.Qscale_orig = [];" << endl;

  mOutputFile << boolalpha
              << "options_.linear = " << linear << ";" << endl
              << "options_.block = " << block << ";" << endl
              << "options_.bytecode = " << bytecode << ";" << endl
              << "options_.use_dll = " << use_dll << ";" << endl;

  if (parallel_local_files.size() > 0)
    {
      mOutputFile << "options_.parallel_info.local_files = {" << endl;
      for (const auto &parallel_local_file : parallel_local_files)
        {
          size_t j = parallel_local_file.find_last_of(R"(/\)");
          if (j == string::npos)
            mOutputFile << "'', '" << parallel_local_file << "';" << endl;
          else
            mOutputFile << "'" << parallel_local_file.substr(0, j+1) << "', '"
                        << parallel_local_file.substr(j+1) << "';" << endl;
        }
      mOutputFile << "};" << endl;
    }

  if (dynamic_model.isHessianComputed())
    {
      mOutputFile << "M_.nonzero_hessian_eqs = ";
      dynamic_model.printNonZeroHessianEquations(mOutputFile);
      mOutputFile << ";" << endl
                  << "M_.hessian_eq_zero = isempty(M_.nonzero_hessian_eqs);" << endl;
    }

  if (!onlymodel)
    config_file.writeCluster(mOutputFile);

  if (bytecode)
    mOutputFile << "if exist('bytecode') ~= 3" << endl
                << "  error('DYNARE: Can''t find bytecode DLL. Please compile it or remove the ''bytecode'' option.')" << endl
                << "end" << endl;

  mOutputFile << "M_.orig_eq_nbr = " << mod_file_struct.orig_eq_nbr << ";" << endl
              << "M_.eq_nbr = " << dynamic_model.equation_number() << ";" << endl
              << "M_.ramsey_eq_nbr = " << mod_file_struct.ramsey_eq_nbr << ";" << endl
              << "M_.set_auxiliary_variables = exist(['./+' M_.fname '/set_auxiliary_variables.m'], 'file') == 2;" << endl;

  epilogue.writeOutput(mOutputFile);

  if (dynamic_model.equation_number() > 0)
    {
      dynamic_model.writeDriverOutput(mOutputFile, basename, block, mod_file_struct.estimation_present, compute_xrefs);
      if (!no_static)
        static_model.writeDriverOutput(mOutputFile, block);
    }

  if (onlymodel || gui)
    for (const auto &statement : statements)
      {
        /* Special treatment for initval block: insert initial values for the
           auxiliary variables and initialize exo det */
        if (auto ivs = dynamic_cast<InitValStatement *>(statement.get()); ivs)
          {
            ivs->writeOutput(mOutputFile, basename, minimal_workspace);
            static_model.writeAuxVarInitval(mOutputFile, ExprNodeOutputType::matlabOutsideModel);
            ivs->writeOutputPostInit(mOutputFile);
          }

        // Special treatment for endval block: insert initial values for the auxiliary variables
        if (auto evs = dynamic_cast<EndValStatement *>(statement.get()); evs)
          {
            evs->writeOutput(mOutputFile, basename, minimal_workspace);
            static_model.writeAuxVarInitval(mOutputFile, ExprNodeOutputType::matlabOutsideModel);
          }

        if (auto ips = dynamic_cast<InitParamStatement *>(statement.get()); ips)
          ips->writeOutput(mOutputFile, basename, minimal_workspace);

        if (auto ss = dynamic_cast<ShocksStatement *>(statement.get()); ss)
          ss->writeOutput(mOutputFile, basename, minimal_workspace);

        if (auto eps = dynamic_cast<EstimatedParamsStatement *>(statement.get()); eps)
          eps->writeOutput(mOutputFile, basename, minimal_workspace);

        if (auto sgs = dynamic_cast<ShockGroupsStatement *>(statement.get()); sgs)
          sgs->writeOutput(mOutputFile, basename, minimal_workspace);

        if (gui)
          if (auto it = dynamic_cast<NativeStatement *>(statement.get()); it)
            it->writeOutput(mOutputFile, basename, minimal_workspace);
      }
  else
    {
      for (const auto &statement : statements)
        {
          statement->writeOutput(mOutputFile, basename, minimal_workspace);

          /* Special treatment for initval block: insert initial values for the
             auxiliary variables and initialize exo det */
          if (auto ivs = dynamic_cast<InitValStatement *>(statement.get()); ivs)
            {
              static_model.writeAuxVarInitval(mOutputFile, ExprNodeOutputType::matlabOutsideModel);
              ivs->writeOutputPostInit(mOutputFile);
            }

          // Special treatment for endval block: insert initial values for the auxiliary variables
          if (auto evs = dynamic_cast<EndValStatement *>(statement.get()); evs)
            static_model.writeAuxVarInitval(mOutputFile, ExprNodeOutputType::matlabOutsideModel);

          // Special treatment for load params and steady state statement: insert initial values for the auxiliary variables
          if (auto lpass = dynamic_cast<LoadParamsAndSteadyStateStatement *>(statement.get()); lpass && !no_static)
            static_model.writeAuxVarInitval(mOutputFile, ExprNodeOutputType::matlabOutsideModel);
        }

      if (!notime)
          mOutputFile << endl << endl
                      << "oo_.time = toc(tic0);" << endl
                      << "disp(['Total computing time : ' dynsec2hms(oo_.time) ]);" << endl;

      mOutputFile << "if ~exist([M_.dname filesep 'Output'],'dir')" << endl
                  << "    mkdir(M_.dname,'Output');" << endl
                  << "end" << endl
                  << "save([M_.dname filesep 'Output' filesep '" << basename << "_results.mat'], 'oo_', 'M_', 'options_');" << endl
                  << "if exist('estim_params_', 'var') == 1" << endl
                  << "  save([M_.dname filesep 'Output' filesep '" << basename << "_results.mat'], 'estim_params_', '-append');" << endl << "end" << endl
                  << "if exist('bayestopt_', 'var') == 1" << endl
                  << "  save([M_.dname filesep 'Output' filesep '" << basename << "_results.mat'], 'bayestopt_', '-append');" << endl << "end" << endl
                  << "if exist('dataset_', 'var') == 1" << endl
                  << "  save([M_.dname filesep 'Output' filesep '" << basename << "_results.mat'], 'dataset_', '-append');" << endl << "end" << endl
                  << "if exist('estimation_info', 'var') == 1" << endl
                  << "  save([M_.dname filesep 'Output' filesep '" << basename << "_results.mat'], 'estimation_info', '-append');" << endl << "end" << endl
                  << "if exist('dataset_info', 'var') == 1" << endl
                  << "  save([M_.dname filesep 'Output' filesep '" << basename << "_results.mat'], 'dataset_info', '-append');" << endl << "end" << endl
                  << "if exist('oo_recursive_', 'var') == 1" << endl
                  << "  save([M_.dname filesep 'Output' filesep '" << basename << "_results.mat'], 'oo_recursive_', '-append');" << endl << "end" << endl;

      config_file.writeEndParallel(mOutputFile);

      if (!no_warn)
        {
          if (warnings.countWarnings() > 0)
            mOutputFile << "disp('Note: " << warnings.countWarnings() << " warning(s) encountered in the preprocessor')" << endl;

          mOutputFile << "if ~isempty(lastwarn)" << endl
                      << "  disp('Note: warning(s) encountered in MATLAB/Octave code')" << endl
                      << "end" << endl;
        }
    }

  mOutputFile.close();

  if (hasModelChanged)
    {
      // Create static and dynamic files
      if (dynamic_model.equation_number() > 0)
        {
          if (!no_static)
            {
              static_model.writeStaticFile(basename, block, use_dll, mexext, matlabroot, dynareroot, false);
              static_model.writeParamsDerivativesFile<false>(basename);
            }

          dynamic_model.writeDynamicFile(basename, block, use_dll, mexext, matlabroot, dynareroot, false);

          dynamic_model.writeParamsDerivativesFile<false>(basename);

          dynamic_model.writeDynamicJacobianNonZeroEltsFile(basename);
        }

      // Create steady state file
      steady_state_model.writeSteadyStateFile(basename, false);

      // Create epilogue file
      epilogue.writeEpilogueFile(basename);

      pac_model_table.writeTargetCoefficientsFile(basename);
    }
}

void
ModFile::writeJuliaOutput(const string &basename) const
{
  if (dynamic_model.equation_number() > 0)
    {
      if (!no_static)
        {
          static_model.writeStaticFile(basename, false, false, "", {}, {}, true);
          static_model.writeParamsDerivativesFile<true>(basename);
        }
      dynamic_model.writeDynamicFile(basename, block, use_dll, "", {}, {}, true);
      dynamic_model.writeParamsDerivativesFile<true>(basename);
    }
  steady_state_model.writeSteadyStateFile(basename, true);
}

void
ModFile::writeJsonOutput(const string &basename, JsonOutputPointType json, JsonFileOutputType json_output_mode, bool onlyjson, bool jsonderivsimple)
{
  if (json == JsonOutputPointType::nojson)
    return;

  if (json == JsonOutputPointType::parsing || json == JsonOutputPointType::checkpass)
    symbol_table.freeze();

  if (json_output_mode == JsonFileOutputType::standardout)
    cout << "//-- BEGIN JSON --// " << endl
         << "{" << endl;

  writeJsonOutputParsingCheck(basename, json_output_mode, json == JsonOutputPointType::transformpass, json == JsonOutputPointType::computingpass);

  if (json == JsonOutputPointType::parsing || json == JsonOutputPointType::checkpass)
    symbol_table.unfreeze();

  if (json == JsonOutputPointType::computingpass)
    writeJsonComputingPassOutput(basename, json_output_mode, jsonderivsimple);

  if (json_output_mode == JsonFileOutputType::standardout)
    cout << "}" << endl
         << "//-- END JSON --// " << endl;

  switch (json)
    {
    case JsonOutputPointType::parsing:
      cout << "JSON written after Parsing step." << endl;
      break;
    case JsonOutputPointType::checkpass:
      cout << "JSON written after Check step." << endl;
      break;
    case JsonOutputPointType::transformpass:
      cout << "JSON written after Transform step." << endl;
      break;
    case JsonOutputPointType::computingpass:
      cout << "JSON written after Computing step." << endl;
      break;
    case JsonOutputPointType::nojson:
      cerr << "ModFile::writeJsonOutput: should not arrive here." << endl;
      exit(EXIT_FAILURE);
    }

  if (onlyjson)
    exit(EXIT_SUCCESS);
}

void
ModFile::writeJsonOutputParsingCheck(const string &basename, JsonFileOutputType json_output_mode, bool transformpass, bool computingpass) const
{
  ostringstream output;
  output << "{" << endl;

  symbol_table.writeJsonOutput(output);
  output << ", ";
  dynamic_model.writeJsonOutput(output);

  if (!statements.empty()
      || !var_model_table.empty()
      || !trend_component_model_table.empty())
    {
      output << R"(, "statements": [)";
      if (!var_model_table.empty())
        {
          var_model_table.writeJsonOutput(output);
          output << ", ";
        }

      if (!trend_component_model_table.empty())
        {
          trend_component_model_table.writeJsonOutput(output);
          output << ", ";
        }

      if (!var_expectation_model_table.empty())
        {
          var_expectation_model_table.writeJsonOutput(output);
          output << ", ";
        }

      if (!pac_model_table.empty())
        {
          pac_model_table.writeJsonOutput(output);
          output << ", ";
        }

      for (bool printed_something{false};
           auto &it : statements)
        {
          if (exchange(printed_something, true))
            output << ", " << endl;
          it->writeJsonOutput(output);
        }
      output << "]" << endl;
    }

  if (computingpass)
    {
      output << ",";
      dynamic_model.writeJsonDynamicModelInfo(output);
    }
  output << "}" << endl;

  ostringstream original_model_output;
  original_model_output << "";
  if (transformpass || computingpass)
    {
      original_model_output << "{";
      original_model.writeJsonOriginalModelOutput(original_model_output);
      if (!statements.empty() || !var_model_table.empty() || !trend_component_model_table.empty())
        {
          original_model_output << endl << R"(, "statements": [)";
          if (!var_model_table.empty())
            {
              var_model_table.writeJsonOutput(original_model_output);
              original_model_output << ", ";
            }
          if (!trend_component_model_table.empty())
            {
              trend_component_model_table.writeJsonOutput(original_model_output);
              original_model_output << ", ";
            }
          if (!pac_model_table.empty())
            {
              pac_model_table.writeJsonOutput(original_model_output);
              original_model_output << ", ";
            }
          for (bool printed_something{false};
               const auto &it : statements)
            {
              original_model_output << (exchange(printed_something, true) ? "," : "") << endl;
              it->writeJsonOutput(original_model_output);
            }
          original_model_output << "]" << endl;
        }
      original_model_output << "}" << endl;
    }

  ostringstream steady_state_model_output;
  steady_state_model_output << "";
  if (dynamic_model.equation_number() > 0)
    steady_state_model.writeJsonSteadyStateFile(steady_state_model_output,
                                                transformpass || computingpass);

  if (json_output_mode == JsonFileOutputType::standardout)
    {
      if (transformpass || computingpass)
        cout << R"("transformed_modfile": )";
      else
        cout << R"("modfile": )";
      cout << output.str();
      if (!original_model_output.str().empty())
        cout << R"(, "original_model": )" << original_model_output.str();
      if (!steady_state_model_output.str().empty())
        cout << R"(, "steady_state_model": )" << steady_state_model_output.str();
    }
  else
    {
      if (!basename.size())
        {
          cerr << "ERROR: Missing file name" << endl;
          exit(EXIT_FAILURE);
        }

      filesystem::create_directories(basename + "/model/json");
      string fname{basename + "/model/json/modfile.json"};
      ofstream jsonOutputFile{fname, ios::out | ios::binary};
      if (!jsonOutputFile.is_open())
        {
          cerr << "ERROR: Can't open file " << fname << " for writing" << endl;
          exit(EXIT_FAILURE);
        }

      jsonOutputFile << output.str();
      jsonOutputFile.close();

      if (!original_model_output.str().empty())
        {
          if (basename.size())
            {
              string fname{basename + "/model/json/modfile-original.json"};
              jsonOutputFile.open(fname, ios::out | ios::binary);
              if (!jsonOutputFile.is_open())
                {
                  cerr << "ERROR: Can't open file " << fname << " for writing" << endl;
                  exit(EXIT_FAILURE);
                }
            }
          else
            {
              cerr << "ERROR: Missing file name" << endl;
              exit(EXIT_FAILURE);
            }

          jsonOutputFile << original_model_output.str();
          jsonOutputFile.close();
        }
      if (!steady_state_model_output.str().empty())
        {
          if (basename.size())
            {
              string fname{basename + "/model/json/steady_state_model.json"};
              jsonOutputFile.open(fname, ios::out | ios::binary);
              if (!jsonOutputFile.is_open())
                {
                  cerr << "ERROR: Can't open file " << fname << " for writing" << endl;
                  exit(EXIT_FAILURE);
                }
            }
          else
            {
              cerr << "ERROR: Missing file name" << endl;
              exit(EXIT_FAILURE);
            }

          jsonOutputFile << steady_state_model_output.str();
          jsonOutputFile.close();
        }
    }
}

void
ModFile::writeJsonComputingPassOutput(const string &basename, JsonFileOutputType json_output_mode, bool jsonderivsimple) const
{
  if (basename.empty() && json_output_mode != JsonFileOutputType::standardout)
    {
      cerr << "ERROR: Missing file name" << endl;
      exit(EXIT_FAILURE);
    }

  ostringstream tmp_out, static_output, dynamic_output, static_paramsd_output, dynamic_paramsd_output;

  static_output << "{";
  static_model.writeJsonComputingPassOutput(static_output, !jsonderivsimple);
  static_output << "}";

  dynamic_output << "{";
  dynamic_model.writeJsonComputingPassOutput(dynamic_output, !jsonderivsimple);
  dynamic_output << "}";

  static_model.writeJsonParamsDerivatives(tmp_out, !jsonderivsimple);
  if (!tmp_out.str().empty())
    static_paramsd_output << "{" << tmp_out.str() << "}" << endl;

  tmp_out.str("");
  dynamic_model.writeJsonParamsDerivatives(tmp_out, !jsonderivsimple);
  if (!tmp_out.str().empty())
    dynamic_paramsd_output << "{" << tmp_out.str() << "}" << endl;

  if (json_output_mode == JsonFileOutputType::standardout)
    {
      cout << R"(, "static_model": )" << static_output.str() << endl
           << R"(, "dynamic_model": )" << dynamic_output.str() << endl;

      if (!static_paramsd_output.str().empty())
        cout << R"(, "static_params_deriv": )" << static_paramsd_output.str() << endl;

      if (!dynamic_paramsd_output.str().empty())
        cout << R"(, "dynamic_params_deriv": )" << dynamic_paramsd_output.str() << endl;
    }
  else
    {
      filesystem::create_directories(basename + "/model/json");

      writeJsonFileHelper(basename + "/model/json/static.json", static_output);
      writeJsonFileHelper(basename + "/model/json/dynamic.json", dynamic_output);

      if (!static_paramsd_output.str().empty())
        writeJsonFileHelper(basename + "/model/json/static_params_derivs.json", static_paramsd_output);

      if (!dynamic_paramsd_output.str().empty())
        writeJsonFileHelper(basename + "/model/json/params_derivs.json", dynamic_paramsd_output);
    }
}

void
ModFile::writeJsonFileHelper(const string &fname, ostringstream &output) const
{
  ofstream jsonOutput{fname, ios::out | ios::binary};
  if (!jsonOutput.is_open())
    {
      cerr << "ERROR: Can't open file " << fname << " for writing" << endl;
      exit(EXIT_FAILURE);
    }
  jsonOutput << output.str();
  jsonOutput.close();
}

filesystem::path
ModFile::unique_path()
{
  filesystem::path path;
  string possible_characters = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
  random_device rd;
  mt19937 generator(rd());
  uniform_int_distribution distribution{0, static_cast<int>(possible_characters.size())-1};
  do
    {
      constexpr int rand_length = 10;
      string rand_str(rand_length, '\0');
      for (auto &dis : rand_str)
        dis = possible_characters[distribution(generator)];
      path = rand_str;
    }
  while (exists(path));

  return path;
}
