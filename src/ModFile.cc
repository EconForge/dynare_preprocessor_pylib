/*
 * Copyright © 2006-2019 Dynare Team
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
    expressions_tree{symbol_table, num_constants, external_functions_table},
    original_model{symbol_table, num_constants, external_functions_table,
                   trend_component_model_table, var_model_table},
    dynamic_model{symbol_table, num_constants, external_functions_table,
                  trend_component_model_table, var_model_table},
    trend_dynamic_model{symbol_table, num_constants, external_functions_table,
                        trend_component_model_table, var_model_table},
    ramsey_FOC_equations_dynamic_model{symbol_table, num_constants, external_functions_table,
                                       trend_component_model_table, var_model_table},
    orig_ramsey_dynamic_model{symbol_table, num_constants, external_functions_table,
                              trend_component_model_table, var_model_table},
    non_linear_equations_dynamic_model{symbol_table, num_constants, external_functions_table,
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
  cout << "Evaluating expressions...";

  // Loop over all statements, and fill global eval context if relevant
  for (auto &st : statements)
    {
      auto ips = dynamic_cast<InitParamStatement *>(st.get());
      if (ips)
        ips->fillEvalContext(global_eval_context);

      auto ies = dynamic_cast<InitOrEndValStatement *>(st.get());
      if (ies)
        ies->fillEvalContext(global_eval_context);

      auto lpass = dynamic_cast<LoadParamsAndSteadyStateStatement *>(st.get());
      if (lpass)
        lpass->fillEvalContext(global_eval_context);
    }

  // Evaluate model local variables
  dynamic_model.fillEvalContext(global_eval_context);

  cout << "done" << endl;

  // Check if some symbols are not initialized, and give them a zero value then
  for (int id = 0; id <= symbol_table.maxID(); id++)
    {
      SymbolType type = symbol_table.getType(id);
      if ((type == SymbolType::endogenous || type == SymbolType::exogenous || type == SymbolType::exogenousDet
           || type == SymbolType::parameter || type == SymbolType::modelLocalVariable)
          && global_eval_context.find(id) == global_eval_context.end())
        {
          if (warn_uninit)
            warnings << "WARNING: Can't find a numeric initial value for "
                     << symbol_table.getName(id) << ", using zero" << endl;
          global_eval_context[id] = 0;
        }
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
  for (auto & statement : statements)
    statement->checkPass(mod_file_struct, warnings);

  // Check the steady state block
  steady_state_model.checkPass(mod_file_struct, warnings);

  // Check epilogue block
  epilogue.checkPass(warnings);

  if (mod_file_struct.write_latex_steady_state_model_present &&
      !mod_file_struct.steady_state_model_present)
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
    || mod_file_struct.sensitivity_present
    || stochastic;

  // Allow empty model only when doing a standalone BVAR estimation
  if (dynamic_model.equation_number() == 0
      && (mod_file_struct.check_present
          || mod_file_struct.perfect_foresight_solver_present
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

  if ((mod_file_struct.osr_present && (!mod_file_struct.osr_params_present || !mod_file_struct.optim_weights_present))
      || ((!mod_file_struct.osr_present || !mod_file_struct.osr_params_present) && mod_file_struct.optim_weights_present)
      || ((!mod_file_struct.osr_present || !mod_file_struct.optim_weights_present) && mod_file_struct.osr_params_present))
    {
      cerr << "ERROR: The osr statement must be used with osr_params and optim_weights." << endl;
      exit(EXIT_FAILURE);
    }

  if (mod_file_struct.perfect_foresight_solver_present && stochastic_statement_present)
    {
      cerr << "ERROR: A .mod file cannot contain both one of {perfect_foresight_solver,simul} and one of {stoch_simul, estimation, osr, ramsey_policy, discretionary_policy}. This is not possible: one cannot mix perfect foresight context with stochastic context in the same file." << endl;
      exit(EXIT_FAILURE);
    }

  if (mod_file_struct.k_order_solver && byte_code)
    {
      cerr << "ERROR: 'k_order_solver' (which is implicit if order >= 3), is not yet compatible with 'bytecode'." << endl;
      exit(EXIT_FAILURE);
    }

  if (use_dll && (block || byte_code || linear_decomposition))
    {
      cerr << "ERROR: In 'model' block, 'use_dll' option is not compatible with 'block', 'bytecode' or 'linear_decomposition'" << endl;
      exit(EXIT_FAILURE);
    }
  if (block && linear_decomposition)
    {
      cerr << "ERROR: In 'model' block, 'block' option is not compatible with 'linear_decomposition'" << endl;
      exit(EXIT_FAILURE);
    }

  if (!byte_code && linear_decomposition)
    {
      cerr << "ERROR: For the moment in 'model' block, 'linear_decomposition' option is compatible only with 'bytecode' option" << endl;
      exit(EXIT_FAILURE);
    }

  if (block || byte_code)
    if (dynamic_model.isModelLocalVariableUsed())
      {
        cerr << "ERROR: In 'model' block, 'block' or 'bytecode' options are not yet compatible with pound expressions" << endl;
        exit(EXIT_FAILURE);
      }

  if ((stochastic_statement_present || mod_file_struct.check_present || mod_file_struct.steady_present) && no_static)
    {
      cerr << "ERROR: no_static option is incompatible with stoch_simul, estimation, osr, ramsey_policy, discretionary_policy, steady and check commands" << endl;
      exit(EXIT_FAILURE);
    }

  if (mod_file_struct.dsge_var_estimated)
    if (!mod_file_struct.dsge_prior_weight_in_estimated_params)
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
      for (auto it = parameters_intersect.begin();
           it != parameters_intersect.end();)
        {
          cerr << symbol_table.getName(*it);
          if (++it != parameters_intersect.end())
            cerr << ", ";
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
}

void
ModFile::transformPass(bool nostrict, bool stochastic, bool compute_xrefs, const bool transform_unary_ops)
{
  /* Save the original model (must be done before any model transformations by preprocessor)
     — except substituting out variables which we know are constant (they
       appear in an equation of the form: X = constant)
     — except adl operators which we always want expanded
     — except diff operators with a lead which have been expanded by
       DataTree:AddDiff()
  */
  dynamic_model.simplifyEquations();
  dynamic_model.substituteAdl();
  dynamic_model.setLeadsLagsOrig();
  original_model = dynamic_model;

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

  // Declare endogenous used for PAC model-consistent expectations
  for (auto & statement : statements)
    {
      auto pms = dynamic_cast<PacModelStatement *>(statement.get());
      if (pms != nullptr)
        {
          if (pms->growth != nullptr)
            pac_growth.push_back(pms->growth);
          if (pms->aux_model_name == "")
            dynamic_model.declarePacModelConsistentExpectationEndogs(pms->name);
        }
    }

  // Get all equation tags associated with VARs and Trend Component Models
  set<string> eqtags;
  for (auto const & it : trend_component_model_table.getEqTags())
    for (auto & it1 : it.second)
      eqtags.insert(it1);

  for (auto const & it : var_model_table.getEqTags())
    for (auto & it1 : it.second)
      eqtags.insert(it1);

  // Create auxiliary variables and equations for unary ops
  lag_equivalence_table_t unary_ops_nodes;
  ExprNode::subst_table_t unary_ops_subst_table;
  if (transform_unary_ops)
    tie(unary_ops_nodes, unary_ops_subst_table) = dynamic_model.substituteUnaryOps();
  else
    // substitute only those unary ops that appear in auxiliary model equations
    tie(unary_ops_nodes, unary_ops_subst_table) = dynamic_model.substituteUnaryOps(eqtags);

  // Create auxiliary variable and equations for Diff operators
  auto [diff_nodes, diff_subst_table] = dynamic_model.substituteDiff(pac_growth);

  // Fill Trend Component Model Table
  dynamic_model.fillTrendComponentModelTable();
  original_model.fillTrendComponentModelTableFromOrigModel();
  dynamic_model.fillTrendComponentmodelTableAREC(diff_subst_table);
  dynamic_model.fillVarModelTable();
  original_model.fillVarModelTableFromOrigModel();

  // Pac Model
  int i = 0;
  for (auto & statement : statements)
    {
      auto pms = dynamic_cast<PacModelStatement *>(statement.get());
      if (pms != nullptr)
         {
           if (pms->growth != nullptr)
             pms->overwriteGrowth(pac_growth.at(i++));

           int max_lag;
           vector<int> lhs;
           vector<bool> nonstationary;
           string aux_model_type;
           if (trend_component_model_table.isExistingTrendComponentModelName(pms->aux_model_name))
             {
               aux_model_type = "trend_component";
               max_lag = trend_component_model_table.getMaxLag(pms->aux_model_name) + 1;
               lhs = dynamic_model.getUndiffLHSForPac(pms->aux_model_name, diff_subst_table);
               // All lhs variables in a trend component model are nonstationary
               nonstationary.insert(nonstationary.end(), trend_component_model_table.getDiff(pms->aux_model_name).size(), true);
             }
           else if (var_model_table.isExistingVarModelName(pms->aux_model_name))
             {
               aux_model_type = "var";
               max_lag = var_model_table.getMaxLag(pms->aux_model_name);
               lhs = var_model_table.getLhs(pms->aux_model_name);
               // nonstationary variables in a VAR are those that are in diff
               nonstationary = var_model_table.getDiff(pms->aux_model_name);
             }
           else if (pms->aux_model_name == "")
             max_lag = 0;
           else
             {
               cerr << "Error: aux_model_name not recognized as VAR model or Trend Component model" << endl;
               exit(EXIT_FAILURE);
             }
           auto eqtag_and_lag = dynamic_model.walkPacParameters(pms->name);
           original_model.getPacMaxLag(pms->name, eqtag_and_lag);
           if (pms->aux_model_name == "")
             dynamic_model.addPacModelConsistentExpectationEquation(pms->name, symbol_table.getID(pms->discount),
                                                                    eqtag_and_lag, diff_subst_table);
           else
             dynamic_model.fillPacModelInfo(pms->name, lhs, max_lag, aux_model_type,
                                            eqtag_and_lag, nonstationary, pms->growth);
           dynamic_model.substitutePacExpectation(pms->name);
         }
     }

  dynamic_model.addEquationsForVar();

  if (symbol_table.predeterminedNbr() > 0)
    dynamic_model.transformPredeterminedVariables();

  // Create auxiliary vars for Expectation operator
  dynamic_model.substituteExpectation(mod_file_struct.partial_information);

  if (nonstationary_variables)
    {
      dynamic_model.detrendEquations();
      trend_dynamic_model = dynamic_model;
      dynamic_model.removeTrendVariableFromEquations();
    }

  mod_file_struct.orig_eq_nbr = dynamic_model.equation_number();
  if (mod_file_struct.ramsey_model_present)
    {
      PlannerObjectiveStatement *pos = nullptr;
      for (auto & statement : statements)
        {
          auto pos2 = dynamic_cast<PlannerObjectiveStatement *>(statement.get());
          if (pos2 != nullptr)
            if (pos != nullptr)
              {
                cerr << "ERROR: there can only be one planner_objective statement" << endl;
                exit(EXIT_FAILURE);
              }
            else
              pos = pos2;
        }
      assert(pos != nullptr);
      const StaticModel &planner_objective = pos->getPlannerObjective();

      /*
        clone the model then clone the new equations back to the original because
        we have to call computeDerivIDs (in computeRamseyPolicyFOCs and computingPass)
      */
      if (linear)
        orig_ramsey_dynamic_model = dynamic_model;
      ramsey_FOC_equations_dynamic_model = dynamic_model;
      ramsey_FOC_equations_dynamic_model.computeRamseyPolicyFOCs(planner_objective);
      ramsey_FOC_equations_dynamic_model.replaceMyEquations(dynamic_model);
      mod_file_struct.ramsey_eq_nbr = dynamic_model.equation_number() - mod_file_struct.orig_eq_nbr;
    }

  /* Handle var_expectation_model statements: collect information about them,
     create the new corresponding parameters, and the expressions to replace
     the var_expectation statements.
     TODO: move information collection to checkPass(), within a new
     VarModelTable class */
  map<string, expr_t> var_expectation_subst_table;
  for (auto & statement : statements)
    {
      auto vems = dynamic_cast<VarExpectationModelStatement *>(statement.get());
      if (!vems)
        continue;

      int max_lag;
      vector<int> lhs;
      auto &model_name = vems->model_name;
      if (var_model_table.isExistingVarModelName(vems->aux_model_name))
        {
          max_lag = var_model_table.getMaxLag(vems->aux_model_name);
          lhs = var_model_table.getLhs(vems->aux_model_name);
        }
      else if (trend_component_model_table.isExistingTrendComponentModelName(vems->aux_model_name))
        {
          max_lag = trend_component_model_table.getMaxLag(vems->aux_model_name) + 1;
          lhs = dynamic_model.getUndiffLHSForPac(vems->aux_model_name, diff_subst_table);
        }
      else
        {
          cerr << "ERROR: var_expectation_model " << model_name
               << " refers to nonexistent auxiliary model " << vems->aux_model_name << endl;
          exit(EXIT_FAILURE);
        }

      /* Substitute unary and diff operators in the 'expression' option, then
         match the linear combination in the expression option */
      vems->substituteUnaryOpNodes(unary_ops_nodes, unary_ops_subst_table);
      vems->substituteDiff(diff_nodes, diff_subst_table);
      vems->matchExpression();

      /* Create auxiliary parameters and the expression to be substituted into
         the var_expectations statement */
      auto subst_expr = dynamic_model.Zero;
      for (int lag = 0; lag < max_lag; lag++)
        for (auto variable : lhs)
          {
            string param_name = "var_expectation_model_" + model_name + '_' + symbol_table.getName(variable) + '_' + to_string(lag);
            int new_param_id = symbol_table.addSymbol(param_name, SymbolType::parameter);
            vems->aux_params_ids.push_back(new_param_id);

            subst_expr = dynamic_model.AddPlus(subst_expr,
                                               dynamic_model.AddTimes(dynamic_model.AddVariable(new_param_id),
                                                                      dynamic_model.AddVariable(variable, -lag)));
          }

      if (var_expectation_subst_table.find(model_name) != var_expectation_subst_table.end())
        {
          cerr << "ERROR: model name '" << model_name << "' is used by several var_expectation_model statements" << endl;
          exit(EXIT_FAILURE);
        }
      var_expectation_subst_table[model_name] = subst_expr;
    }
  // And finally perform the substitutions
  dynamic_model.substituteVarExpectation(var_expectation_subst_table);

  if (mod_file_struct.stoch_simul_present
      || mod_file_struct.estimation_present
      || mod_file_struct.osr_present
      || mod_file_struct.ramsey_policy_present
      || mod_file_struct.discretionary_policy_present
      || mod_file_struct.calib_smoother_present
      || mod_file_struct.identification_present
      || mod_file_struct.sensitivity_present
      || stochastic )
    {
      // In stochastic models, create auxiliary vars for leads and lags greater than 2, on both endos and exos
      dynamic_model.substituteEndoLeadGreaterThanTwo(false);
      dynamic_model.substituteExoLead(false);
      dynamic_model.substituteEndoLagGreaterThanTwo(false);
      dynamic_model.substituteExoLag(false);
    }
  else
    {
      // In deterministic models, create auxiliary vars for leads and lags endogenous greater than 2, only on endos (useless on exos)
      dynamic_model.substituteEndoLeadGreaterThanTwo(true);
      dynamic_model.substituteEndoLagGreaterThanTwo(true);
    }

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
      && !(mod_file_struct.occbin_option)
      && (dynamic_model.equation_number() != symbol_table.endo_nbr()))
    {
      cerr << "ERROR: There are " << dynamic_model.equation_number() << " equations but " << symbol_table.endo_nbr() << " endogenous variables!" << endl;
      exit(EXIT_FAILURE);
    }

  if (symbol_table.exo_det_nbr() > 0 && mod_file_struct.perfect_foresight_solver_present)
    {
      cerr << "ERROR: A .mod file cannot contain both one of {perfect_foresight_solver, simul} and varexo_det declaration (all exogenous variables are deterministic in this case)" << endl;
      exit(EXIT_FAILURE);
    }

  if (mod_file_struct.ramsey_policy_present && symbol_table.exo_det_nbr() > 0)
    {
      cerr << "ERROR: ramsey_policy is incompatible with deterministic exogenous variables" << endl;
      exit(EXIT_FAILURE);
    }

  if (mod_file_struct.ramsey_policy_present)
    for (auto & statement : statements)
      {
        auto *rps = dynamic_cast<RamseyPolicyStatement *>(statement.get());
        if (rps != nullptr)
          rps->checkRamseyPolicyList();
      }

  if (mod_file_struct.identification_present && symbol_table.exo_det_nbr() > 0)
    {
      cerr << "ERROR: identification is incompatible with deterministic exogenous variables" << endl;
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
ModFile::computingPass(bool no_tmp_terms, FileOutputType output, int params_derivs_order)
{
  // Mod file may have no equation (for example in a standalone BVAR estimation)
  if (dynamic_model.equation_number() > 0)
    {
      if (nonstationary_variables)
        trend_dynamic_model.runTrendTest(global_eval_context);

      // Compute static model and its derivatives
      static_model = static_cast<StaticModel>(dynamic_model);
	  if (linear_decomposition)
        {
          non_linear_equations_dynamic_model = dynamic_model;
          non_linear_equations_dynamic_model.set_cutoff_to_zero();
          non_linear_equations_dynamic_model.computingPass(true, 1, 0, global_eval_context, no_tmp_terms, block, use_dll, byte_code, linear_decomposition);
        }
      if (!no_static)
        {
          if (mod_file_struct.stoch_simul_present
              || mod_file_struct.estimation_present || mod_file_struct.osr_present
              || mod_file_struct.ramsey_model_present || mod_file_struct.identification_present
              || mod_file_struct.calib_smoother_present)
            static_model.set_cutoff_to_zero();

          int derivsOrder = 1;
          int paramsDerivsOrder = 0;
          if (mod_file_struct.identification_present || mod_file_struct.estimation_analytic_derivation)
            {
              derivsOrder = 2;
              paramsDerivsOrder = params_derivs_order;
            }
          static_model.computingPass(derivsOrder, paramsDerivsOrder, global_eval_context, no_tmp_terms, block, byte_code);
        }
      // Set things to compute for dynamic model
      if (mod_file_struct.perfect_foresight_solver_present || mod_file_struct.check_present
          || mod_file_struct.stoch_simul_present
          || mod_file_struct.estimation_present || mod_file_struct.osr_present
          || mod_file_struct.ramsey_model_present || mod_file_struct.identification_present
          || mod_file_struct.calib_smoother_present)
        {
          if (mod_file_struct.perfect_foresight_solver_present)
            dynamic_model.computingPass(true, 1, 0, global_eval_context, no_tmp_terms, block, use_dll, byte_code, linear_decomposition);
          else
            {
              if (mod_file_struct.stoch_simul_present
                  || mod_file_struct.estimation_present || mod_file_struct.osr_present
                  || mod_file_struct.ramsey_model_present || mod_file_struct.identification_present
                  || mod_file_struct.calib_smoother_present)
                dynamic_model.set_cutoff_to_zero();
              if (mod_file_struct.order_option < 1)
                {
                  cerr << "ERROR: Incorrect order option..." << endl;
                  exit(EXIT_FAILURE);
                }
              int derivsOrder = max(mod_file_struct.order_option,
                                    mod_file_struct.identification_order + 1);
              if (mod_file_struct.sensitivity_present || linear || output == FileOutputType::second)
                derivsOrder = max(derivsOrder, 2);
              if (mod_file_struct.estimation_analytic_derivation || output == FileOutputType::third)
                derivsOrder = max(derivsOrder, 3);
              int paramsDerivsOrder = 0;
              if (mod_file_struct.identification_present || mod_file_struct.estimation_analytic_derivation)
                paramsDerivsOrder = params_derivs_order;
              dynamic_model.computingPass(true, derivsOrder, paramsDerivsOrder, global_eval_context, no_tmp_terms, block, use_dll, byte_code, linear_decomposition);
              if (linear && mod_file_struct.ramsey_model_present)
                orig_ramsey_dynamic_model.computingPass(true, 2, paramsDerivsOrder, global_eval_context, no_tmp_terms, block, use_dll, byte_code, linear_decomposition);
            }
        }
      else // No computing task requested, compute derivatives up to 2nd order by default
        dynamic_model.computingPass(true, 2, 0, global_eval_context, no_tmp_terms, block, use_dll, byte_code, linear_decomposition);

      map<int, string> eqs;
      if (mod_file_struct.ramsey_model_present)
        orig_ramsey_dynamic_model.setNonZeroHessianEquations(eqs);
      else
        dynamic_model.setNonZeroHessianEquations(eqs);

      if (linear && !eqs.empty())
        {
          cerr << "ERROR: If the model is declared linear the second derivatives must be equal to zero." << endl
               << "       The following equations had non-zero second derivatives:" << endl;
          for (map<int, string >::const_iterator it = eqs.begin(); it != eqs.end(); it++)
            {
              cerr << "       * Eq # " << it->first+1;
              if (!it->second.empty())
                cerr << " [" << it->second << "]";
              cerr << endl;
            }
          exit(EXIT_FAILURE);
        }
    }

  for (auto & statement : statements)
    statement->computingPass();

  // Compute epilogue derivatives (but silence standard output)
  streambuf *oldcout = cout.rdbuf();
  cout.rdbuf(nullptr);
  epilogue.computingPass(true, 2, 0, global_eval_context, true, false, false, false, false);
  cout.rdbuf(oldcout);
}

void
ModFile::writeOutputFiles(const string &basename, bool clear_all, bool clear_global, bool no_log, bool no_warn,
                          bool console, bool nograph, bool nointeractive, const ConfigFile &config_file,
                          bool check_model_changes, bool minimal_workspace, bool compute_xrefs,
                          const string &mexext,
                          const filesystem::path &matlabroot,
                          const filesystem::path &dynareroot, bool onlymodel) const
{
  bool hasModelChanged = !dynamic_model.isChecksumMatching(basename, block) || !check_model_changes;
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
      if (filesystem::exists("+" + basename))
	{
	  auto tmp = unique_path();
	  filesystem::rename("+" + basename, tmp);
	  filesystem::remove_all(tmp);
	}
      filesystem::remove_all(basename + "/model/src");
      filesystem::remove_all(basename + "/model/bytecode");
    }

  ofstream mOutputFile;

  if (basename.size())
    {
      filesystem::create_directory("+" + basename);
      string fname = "+" + basename + "/driver.m";
      mOutputFile.open(fname, ios::out | ios::binary);
      if (!mOutputFile.is_open())
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

  mOutputFile << "tic0 = tic;" << endl
              << "% Define global variables." << endl
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
  mOutputFile << "global_initialization;" << endl
              << "diary off;" << endl;
  if (!no_log)
    mOutputFile << "diary('" << basename << ".log');" << endl;

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

  cout << "Processing outputs ..." << endl;

  symbol_table.writeOutput(mOutputFile);

  var_model_table.writeOutput(basename, mOutputFile);
  trend_component_model_table.writeOutput(basename, mOutputFile);

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

  // Initialize M_.det_shocks
  mOutputFile << "M_.det_shocks = [];" << endl;

  auto to_matlab_logical = [](bool m) { return m ? "true" : "false"; };

  mOutputFile << "options_.linear = " << to_matlab_logical(linear) << ";" << endl
              << "options_.block = " << to_matlab_logical(block) << ";" << endl
              << "options_.bytecode = " << to_matlab_logical(byte_code) << ";" << endl
              << "options_.use_dll = " << to_matlab_logical(use_dll) << ";" << endl
              << "options_.linear_decomposition = " << to_matlab_logical(linear_decomposition) << ";" << endl;

  if (parallel_local_files.size() > 0)
    {
      mOutputFile << "options_.parallel_info.local_files = {" << endl;
      for (const auto & parallel_local_file : parallel_local_files)
        {
          size_t j = parallel_local_file.find_last_of(R"(/\)");
          if (j == string::npos)
            mOutputFile << "'', '" << parallel_local_file << "';" << endl;
          else
            mOutputFile << "'" << parallel_local_file.substr(0, j+1) << "', '"
                        << parallel_local_file.substr(j+1, string::npos) << "';" << endl;
        }
      mOutputFile << "};" << endl;
    }

  mOutputFile << "M_.nonzero_hessian_eqs = ";
  if (mod_file_struct.ramsey_model_present)
    orig_ramsey_dynamic_model.printNonZeroHessianEquations(mOutputFile);
  else
    dynamic_model.printNonZeroHessianEquations(mOutputFile);
  mOutputFile << ";" << endl
              << "M_.hessian_eq_zero = isempty(M_.nonzero_hessian_eqs);" << endl;

  if (!onlymodel)
    config_file.writeCluster(mOutputFile);

  if (byte_code)
    mOutputFile << "if exist('bytecode') ~= 3" << endl
                << "  error('DYNARE: Can''t find bytecode DLL. Please compile it or remove the ''bytecode'' option.')" << endl
                << "end" << endl;

  mOutputFile << "M_.orig_eq_nbr = " << mod_file_struct.orig_eq_nbr << ";" << endl
              << "M_.eq_nbr = " << dynamic_model.equation_number() << ";" << endl
              << "M_.ramsey_eq_nbr = " << mod_file_struct.ramsey_eq_nbr << ";" << endl
              << "M_.set_auxiliary_variables = exist(['./+' M_.fname '/set_auxiliary_variables.m'], 'file') == 2;" << endl;

  if (dynamic_model.equation_number() > 0)
    {
      if (linear_decomposition)
        non_linear_equations_dynamic_model.writeOutput(mOutputFile, basename, block, true, byte_code, use_dll, mod_file_struct.order_option, mod_file_struct.estimation_present, compute_xrefs, false);
      dynamic_model.writeOutput(mOutputFile, basename, block, false, byte_code, use_dll, mod_file_struct.order_option, mod_file_struct.estimation_present, compute_xrefs, false);
      if (!no_static)
        static_model.writeOutput(mOutputFile, block);
    }

  if (onlymodel)
    for (auto &statement : statements)
      {
        /* Special treatment for initval block: insert initial values for the
           auxiliary variables and initialize exo det */
        auto *ivs = dynamic_cast<InitValStatement *>(statement.get());
        if (ivs != nullptr)
          {
            ivs->writeOutput(mOutputFile, basename, minimal_workspace);
            static_model.writeAuxVarInitval(mOutputFile, ExprNodeOutputType::matlabOutsideModel);
            ivs->writeOutputPostInit(mOutputFile);
          }

        // Special treatment for endval block: insert initial values for the auxiliary variables
        auto *evs = dynamic_cast<EndValStatement *>(statement.get());
        if (evs != nullptr)
          {
            evs->writeOutput(mOutputFile, basename, minimal_workspace);
            static_model.writeAuxVarInitval(mOutputFile, ExprNodeOutputType::matlabOutsideModel);
          }

        auto *ips = dynamic_cast<InitParamStatement *>(statement.get());
        if (ips != nullptr)
          ips->writeOutput(mOutputFile, basename, minimal_workspace);

        auto *ss = dynamic_cast<ShocksStatement *>(statement.get());
        if (ss != nullptr)
          ss->writeOutput(mOutputFile, basename, minimal_workspace);

        auto *eps = dynamic_cast<EstimatedParamsStatement *>(statement.get());
        if (eps != nullptr)
          eps->writeOutput(mOutputFile, basename, minimal_workspace);

        auto *sgs = dynamic_cast<ShockGroupsStatement *>(statement.get());
        if (sgs != nullptr)
          sgs->writeOutput(mOutputFile, basename, minimal_workspace);
      }
  else
    {
      for (auto &statement : statements)
        {
          statement->writeOutput(mOutputFile, basename, minimal_workspace);

          /* Special treatment for initval block: insert initial values for the
             auxiliary variables and initialize exo det */
          auto ivs = dynamic_cast<InitValStatement *>(statement.get());
          if (ivs != nullptr)
            {
              static_model.writeAuxVarInitval(mOutputFile, ExprNodeOutputType::matlabOutsideModel);
              ivs->writeOutputPostInit(mOutputFile);
            }

          // Special treatment for endval block: insert initial values for the auxiliary variables
          auto evs = dynamic_cast<EndValStatement *>(statement.get());
          if (evs != nullptr)
            static_model.writeAuxVarInitval(mOutputFile, ExprNodeOutputType::matlabOutsideModel);

          // Special treatment for load params and steady state statement: insert initial values for the auxiliary variables
          auto lpass = dynamic_cast<LoadParamsAndSteadyStateStatement *>(statement.get());
          if (lpass && !no_static)
            static_model.writeAuxVarInitval(mOutputFile, ExprNodeOutputType::matlabOutsideModel);
        }

      mOutputFile << "save('" << basename << "_results.mat', 'oo_', 'M_', 'options_');" << endl
                  << "if exist('estim_params_', 'var') == 1" << endl
                  << "  save('" << basename << "_results.mat', 'estim_params_', '-append');" << endl << "end" << endl
                  << "if exist('bayestopt_', 'var') == 1" << endl
                  << "  save('" << basename << "_results.mat', 'bayestopt_', '-append');" << endl << "end" << endl
                  << "if exist('dataset_', 'var') == 1" << endl
                  << "  save('" << basename << "_results.mat', 'dataset_', '-append');" << endl << "end" << endl
                  << "if exist('estimation_info', 'var') == 1" << endl
                  << "  save('" << basename << "_results.mat', 'estimation_info', '-append');" << endl << "end" << endl
                  << "if exist('dataset_info', 'var') == 1" << endl
                  << "  save('" << basename << "_results.mat', 'dataset_info', '-append');" << endl << "end" << endl
                  << "if exist('oo_recursive_', 'var') == 1" << endl
                  << "  save('" << basename << "_results.mat', 'oo_recursive_', '-append');" << endl << "end" << endl;

      config_file.writeEndParallel(mOutputFile);

      mOutputFile << endl << endl
                  << "disp(['Total computing time : ' dynsec2hms(toc(tic0)) ]);" << endl;

      if (!no_warn)
        {
          if (warnings.countWarnings() > 0)
            mOutputFile << "disp('Note: " << warnings.countWarnings() << " warning(s) encountered in the preprocessor')" << endl;

          mOutputFile << "if ~isempty(lastwarn)" << endl
                      << "  disp('Note: warning(s) encountered in MATLAB/Octave code')" << endl
                      << "end" << endl;
        }
    }

  if (!no_log)
    mOutputFile << "diary off" << endl;

  mOutputFile.close();

  if (hasModelChanged)
    {
      // Create static and dynamic files
      if (dynamic_model.equation_number() > 0)
        {
          if (!no_static)
            {
              static_model.writeStaticFile(basename, block, byte_code, use_dll, mexext, matlabroot, dynareroot, false);
              static_model.writeParamsDerivativesFile(basename, false);
            }

          if (linear_decomposition)
            {
              non_linear_equations_dynamic_model.writeDynamicFile(basename, block, linear_decomposition, byte_code, use_dll, mexext, matlabroot, dynareroot, mod_file_struct.order_option, false);
              non_linear_equations_dynamic_model.writeParamsDerivativesFile(basename, false);
            }

          dynamic_model.writeDynamicFile(basename, block, false, byte_code, use_dll, mexext, matlabroot, dynareroot, mod_file_struct.order_option, false);

          dynamic_model.writeParamsDerivativesFile(basename, false);

          dynamic_model.writeDynamicJacobianNonZeroElts(basename);
        }

      // Create steady state file
      steady_state_model.writeSteadyStateFile(basename, mod_file_struct.ramsey_model_present, false);

      // Create epilogue file
      epilogue.writeEpilogueFile(basename);
    }

  cout << "done" << endl;
}

void
ModFile::writeExternalFiles(const string &basename, LanguageOutputType language) const
{
  switch (language)
    {
    case LanguageOutputType::julia:
      writeExternalFilesJulia(basename);
      break;
    case LanguageOutputType::matlab:
      cerr << "The 'output' option cannot be used when language=matlab" << endl;
      exit(EXIT_FAILURE);
    }
}

void
ModFile::writeExternalFilesJulia(const string &basename) const
{
  ofstream jlOutputFile;
  if (basename.size())
    {
      string fname(basename);
      fname += ".jl";
      jlOutputFile.open(fname, ios::out | ios::binary);
      if (!jlOutputFile.is_open())
        {
          cerr << "ERROR: Can't open file " << fname
               << " for writing" << endl;
          exit(EXIT_FAILURE);
        }
    }
  else
    {
      cerr << "ERROR: Missing file name" << endl;
      exit(EXIT_FAILURE);
    }

  jlOutputFile << "module " << basename << endl
               << "#" << endl
               << "# NB: this file was automatically generated by Dynare" << endl
               << "#     from " << basename << ".mod" << endl
               << "#" << endl << endl
               << "using DynareModel" << endl
               << "using DynareOptions" << endl
               << "using DynareOutput" << endl << endl
               << "using Utils" << endl
               << "using SteadyState" << endl << endl
               << "using " << basename << "Static" << endl
               << "using " << basename << "Dynamic" << endl
               << R"(if isfile(")" << basename << R"(SteadyState.jl"))" << endl
               << "    using " << basename << "SteadyState" << endl
               << "end" << endl
               << R"(if isfile(")" << basename << R"(SteadyState2.jl"))" << endl
               << "    using " << basename << "SteadyState2" << endl
               << "end" << endl << endl
               << "export model_, options_, oo_" << endl;

  // Write Output
  jlOutputFile << endl
               << "oo_ = dynare_output()" << endl
               << R"(oo_.dynare_version = ")" << PACKAGE_VERSION << R"(")" << endl;

  // Write Options
  jlOutputFile << endl
               << "options_ = dynare_options()" << endl
               << R"(options_.dynare_version = ")" << PACKAGE_VERSION << R"(")" << endl;
  if (linear)
    jlOutputFile << "options_.linear = true" << endl;

  // Write Model
  jlOutputFile << endl
               << "model_ = dynare_model()" << endl
               << R"(model_.fname = ")" << basename << R"(")" << endl
               << R"(model_.dynare_version = ")" << PACKAGE_VERSION << R"(")" << endl
               << "model_.sigma_e = zeros(Float64, " << symbol_table.exo_nbr() << ", "
               << symbol_table.exo_nbr() << ")" << endl
               << "model_.correlation_matrix = ones(Float64, " << symbol_table.exo_nbr() << ", "
               << symbol_table.exo_nbr() << ")" << endl
               << "model_.orig_eq_nbr = " << mod_file_struct.orig_eq_nbr << endl
               << "model_.eq_nbr = " << dynamic_model.equation_number() << endl
               << "model_.ramsey_eq_nbr = " << mod_file_struct.ramsey_eq_nbr << endl;

  if (mod_file_struct.calibrated_measurement_errors)
    jlOutputFile << "model_.h = zeros(Float64,"
                 << symbol_table.observedVariablesNbr() << ", "
                 << symbol_table.observedVariablesNbr() << ");" << endl
                 << "model_.correlation_matrix_me = ones(Float64, "
                 << symbol_table.observedVariablesNbr() << ", "
                 << symbol_table.observedVariablesNbr() << ");" << endl;
  else
    jlOutputFile << "model_.h = zeros(Float64, 1, 1)" << endl
                 << "model_.correlation_matrix_me = ones(Float64, 1, 1)" << endl;

  cout << "Processing outputs ..." << endl;
  symbol_table.writeJuliaOutput(jlOutputFile);

  if (dynamic_model.equation_number() > 0)
    {
      dynamic_model.writeOutput(jlOutputFile, basename, false, false, false, false,
                                mod_file_struct.order_option,
                                mod_file_struct.estimation_present, false, true);
      if (!no_static)
        {
          static_model.writeStaticFile(basename, false, false, false, "", {}, {}, true);
          static_model.writeParamsDerivativesFile(basename, true);
        }
      dynamic_model.writeDynamicFile(basename, block, linear_decomposition, byte_code, use_dll,
                                     "", {}, {}, mod_file_struct.order_option, true);
      dynamic_model.writeParamsDerivativesFile(basename, true);
    }
  steady_state_model.writeSteadyStateFile(basename, mod_file_struct.ramsey_model_present, true);

  // Print statements (includes parameter values)
  for (auto &statement : statements)
    statement->writeJuliaOutput(jlOutputFile, basename);

  jlOutputFile << "model_.static = " << basename << "Static.static!" << endl
               << "model_.dynamic = " << basename << "Dynamic.dynamic!" << endl
               << "model_.temporaries.static = " << basename << "Static.tmp_nbr" << endl
               << "model_.temporaries.dynamic = " << basename << "Dynamic.tmp_nbr" << endl
               << R"(if isfile(")" << basename << R"(SteadyState.jl"))" << endl
               << "    model_.user_written_analytical_steady_state = true" << endl
               << "    model_.steady_state = " << basename << "SteadyState.steady_state!" << endl
               << "end" << endl
               << R"(if isfile(")" << basename << R"(SteadyState2.jl"))" << endl
               << "    model_.analytical_steady_state = true" << endl
               << "    model_.steady_state = " << basename << "SteadyState2.steady_state!" << endl
               << "end" << endl
               << R"(if isfile(")" << basename << R"(StaticParamsDerivs.jl"))" << endl
               << "    using " << basename << "StaticParamsDerivs" << endl
               << "    model_.static_params_derivs = " << basename << "StaticParamsDerivs.params_derivs" << endl
               << "end" << endl
               << R"(if isfile(")" << basename << R"(DynamicParamsDerivs.jl"))" << endl
               << "    using " << basename << "DynamicParamsDerivs" << endl
               << "    model_.dynamic_params_derivs = " << basename << "DynamicParamsDerivs.params_derivs" << endl
               << "end" << endl
               << "end" << endl;
  jlOutputFile.close();
  cout << "done" << endl;
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

      for (auto it = statements.begin();
           it != statements.end(); it++)
        {
          if (it != statements.begin())
            output << ", " << endl;
          (*it)->writeJsonOutput(output);
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
          int i = 0;
          for (const auto & it : statements)
            {
              original_model_output << (i++ > 0 ? "," : "") << endl;
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
      ofstream jsonOutputFile;

      if (basename.size())
        {
          filesystem::create_directories(basename + "/model/json");
          string fname{basename + "/model/json/modfile.json"};
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

  tmp_out << "";
  static_paramsd_output << "";
  static_model.writeJsonParamsDerivativesFile(tmp_out, !jsonderivsimple);
  if (!tmp_out.str().empty())
    static_paramsd_output << "{" << tmp_out.str() << "}" << endl;

  tmp_out.str("");
  dynamic_paramsd_output << "";
  dynamic_model.writeJsonParamsDerivativesFile(tmp_out, !jsonderivsimple);
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
  ofstream jsonOutput;
  jsonOutput.open(fname, ios::out | ios::binary);
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
  uniform_int_distribution<int> distribution{0, static_cast<int>(possible_characters.size())-1};
  do
    {
      constexpr int rand_length = 10;
      string rand_str(rand_length, '\0');
      for (auto &dis : rand_str)
        dis = possible_characters[distribution(generator)];
      path = rand_str;
    }
  while (filesystem::exists(path));

  return path;
}
