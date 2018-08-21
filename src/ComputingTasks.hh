/*
 * Copyright (C) 2003-2018 Dynare Team
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

#ifndef _COMPUTINGTASKS_HH
#define _COMPUTINGTASKS_HH

#include <ostream>

#include "SymbolList.hh"
#include "SymbolTable.hh"
#include "Statement.hh"
#include "StaticModel.hh"
#include "DynamicModel.hh"
#include "SteadyStateModel.hh"

class SteadyStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  SteadyStatement(OptionsList options_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class CheckStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  CheckStatement(OptionsList options_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class SimulStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  SimulStatement(OptionsList options_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class PerfectForesightSetupStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  PerfectForesightSetupStatement(OptionsList options_list_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class PerfectForesightSolverStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  PerfectForesightSolverStatement(OptionsList options_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class PriorPosteriorFunctionStatement : public Statement
{
private:
  const bool prior_func;
  const OptionsList options_list;
public:
  PriorPosteriorFunctionStatement(const bool prior_func_arg, OptionsList options_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class ModelInfoStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  ModelInfoStatement(OptionsList options_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class StochSimulStatement : public Statement
{
private:
  const SymbolList symbol_list;
  const OptionsList options_list;
public:
  StochSimulStatement(SymbolList symbol_list_arg,
                      OptionsList options_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class PacModelStatement : public Statement
{
private:
  const string name;
  const string aux_model_name;
  const string discount;
  const string growth;
  const SymbolTable &symbol_table;
  vector<int> lhs;
public:
  PacModelStatement(string name_arg,
                    string aux_model_name_arg,
                    string discount_arg,
                    string growth_arg,
                    const SymbolTable &symbol_table_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
  void fillUndiffedLHS(vector<int> &lhs);
  tuple<string, string, int> getPacModelInfoForPacExpectation() const;
};

class VarRestrictionsStatement : public Statement
{
private:
  using var_restriction_eq_crosseq_t = pair<pair<int, pair<int, int>>, expr_t>;
  const string var_model_name;
  const map<string, vector<string>> var_map;
  const map<int, map<int, SymbolList>> exclusion_restrictions;
  using equation_restrictions_t = map<int, pair<pair<var_restriction_eq_crosseq_t, var_restriction_eq_crosseq_t>, double>>;
  const equation_restrictions_t equation_restrictions;
  using crossequation_restrictions_t = vector<pair<pair<var_restriction_eq_crosseq_t, var_restriction_eq_crosseq_t>, double>>;
  const crossequation_restrictions_t crossequation_restrictions;
  const map<pair<int, int>, double> covariance_number_restriction;
  const map<pair<int, int>, pair<int, int>> covariance_pair_restriction;
  const SymbolTable &symbol_table;
  int findIdxInVector(const vector<string> &vecvars, const string &var) const;
public:
  VarRestrictionsStatement(string var_model_name_arg,
                           map<string, vector<string>> var_map_arg,
                           map<int, map<int, SymbolList>> exclusion_restrictions_arg,
                           equation_restrictions_t equation_restrictions_arg,
                           crossequation_restrictions_t crossequation_restrictions_arg,
                           map<pair<int, int>, double> covariance_number_restriction_arg,
                           map<pair<int, int>, pair<int, int>> covariance_pair_restriction_arg,
                           const SymbolTable &symbol_table_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
};

class VarEstimationStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  VarEstimationStatement(OptionsList options_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
};

class ForecastStatement : public Statement
{
private:
  const SymbolList symbol_list;
  const OptionsList options_list;
public:
  ForecastStatement(SymbolList symbol_list_arg,
                    OptionsList options_list_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class RamseyModelStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  RamseyModelStatement(OptionsList options_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class RamseyConstraintsStatement : public Statement
{
public:
  struct Constraint
  {
    int endo;
    BinaryOpcode code;
    expr_t expression;
  };
  using constraints_t = vector<Constraint>;
private:
  const SymbolTable &symbol_table;
  const constraints_t constraints;
public:
  RamseyConstraintsStatement(const SymbolTable &symbol_table_arg, constraints_t constraints_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class RamseyPolicyStatement : public Statement
{
private:
  const SymbolTable &symbol_table;
  const vector<string> ramsey_policy_list;
  const OptionsList options_list;
public:
  RamseyPolicyStatement(const SymbolTable &symbol_table_arg,
                        vector<string> ramsey_policy_list_arg,
                        OptionsList options_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void checkRamseyPolicyList();
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class DiscretionaryPolicyStatement : public Statement
{
private:
  const SymbolList symbol_list;
  const OptionsList options_list;
public:
  DiscretionaryPolicyStatement(SymbolList symbol_list_arg,
                               OptionsList options_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class RplotStatement : public Statement
{
private:
  const SymbolList symbol_list;
public:
  RplotStatement(SymbolList symbol_list_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class UnitRootVarsStatement : public Statement
{
public:
  UnitRootVarsStatement();
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class PeriodsStatement : public Statement
{
private:
  const int periods;
public:
  PeriodsStatement(int periods_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class DsampleStatement : public Statement
{
private:
  const int val1, val2;
public:
  DsampleStatement(int val1_arg);
  DsampleStatement(int val1_arg, int val2_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class EstimationStatement : public Statement
{
private:
  const SymbolList symbol_list;
  const OptionsList options_list;
public:
  EstimationStatement(SymbolList symbol_list_arg,
                      OptionsList options_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class DynareSensitivityStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  DynareSensitivityStatement(OptionsList options_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class ObservationTrendsStatement : public Statement
{
public:
  using trend_elements_t = map<string, expr_t>;
private:
  const trend_elements_t trend_elements;
  const SymbolTable &symbol_table;
public:
  ObservationTrendsStatement(trend_elements_t trend_elements_arg,
                             const SymbolTable &symbol_table_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class OsrParamsStatement : public Statement
{
private:
  const SymbolList symbol_list;
  const SymbolTable &symbol_table;
public:
  OsrParamsStatement(SymbolList symbol_list_arg, const SymbolTable &symbol_table_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class OsrStatement : public Statement
{
private:
  const SymbolList symbol_list;
  const OptionsList options_list;
public:
  OsrStatement(SymbolList symbol_list_arg,
               OptionsList options_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

//! Temporary structure used when parsing estimation_params* statements
class OsrParams
{
public:
  string name;
  expr_t low_bound, up_bound;

  void
  init(const DataTree &datatree)
  {
    name = "";
    low_bound = datatree.MinusInfinity;
    up_bound = datatree.Infinity;
  }
};

class OsrParamsBoundsStatement : public Statement
{
private:
  const vector<OsrParams> osr_params_list;
public:
  OsrParamsBoundsStatement(vector<OsrParams> osr_params_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class DynaTypeStatement : public Statement
{
private:
  const SymbolList symbol_list;
  const string filename;
public:
  DynaTypeStatement(SymbolList symbol_list_arg,
                    string filename_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class DynaSaveStatement : public Statement
{
private:
  const SymbolList symbol_list;
  const string filename;
public:
  DynaSaveStatement(SymbolList symbol_list_arg,
                    string filename_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class ModelComparisonStatement : public Statement
{
public:
  using filename_list_t = vector<pair<string, string>>;
private:
  filename_list_t filename_list;
  OptionsList options_list;
public:
  ModelComparisonStatement(filename_list_t filename_list_arg,
                           OptionsList options_list_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

//! Temporary structure used when parsing estimation_params* statements
class EstimationParams
{
public:
  int type;
  string name, name2;
  PriorDistributions prior;
  expr_t init_val, low_bound, up_bound, mean, std, p3, p4, jscale;

  void
  init(const DataTree &datatree)
  {
    type = 0;
    name = "";
    name2 = "";
    prior = PriorDistributions::noShape;
    init_val = datatree.NaN;
    low_bound = datatree.MinusInfinity;
    up_bound = datatree.Infinity;
    mean = datatree.NaN;
    std = datatree.NaN;
    p3 = datatree.NaN;
    p4 = datatree.NaN;
    jscale = datatree.NaN;
  }
};

class EstimatedParamsStatement : public Statement
{
private:
  const vector<EstimationParams> estim_params_list;
  const SymbolTable &symbol_table;
public:
  EstimatedParamsStatement(vector<EstimationParams> estim_params_list_arg,
                           const SymbolTable &symbol_table_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class EstimatedParamsInitStatement : public Statement
{
private:
  const vector<EstimationParams> estim_params_list;
  const SymbolTable &symbol_table;
  const bool use_calibration;
public:
  EstimatedParamsInitStatement(vector<EstimationParams> estim_params_list_arg,
                               const SymbolTable &symbol_table_arg,
                               const bool use_calibration_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class EstimatedParamsBoundsStatement : public Statement
{
private:
  const vector<EstimationParams> estim_params_list;
  const SymbolTable &symbol_table;
public:
  EstimatedParamsBoundsStatement(vector<EstimationParams> estim_params_list_arg,
                                 const SymbolTable &symbol_table_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class OptimWeightsStatement : public Statement
{
public:
  using var_weights_t = map<string, expr_t>;
  using covar_weights_t = map<pair<string, string>, expr_t>;
private:
  const var_weights_t var_weights;
  const covar_weights_t covar_weights;
  const SymbolTable &symbol_table;
public:
  OptimWeightsStatement(var_weights_t var_weights_arg,
                        covar_weights_t covar_weights_arg,
                        const SymbolTable &symbol_table_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class PlannerObjectiveStatement : public Statement
{
private:
  StaticModel model_tree;
  bool computing_pass_called{false};
public:
  PlannerObjectiveStatement(SymbolTable &symbol_table,
                            NumericalConstants &num_constants,
                            ExternalFunctionsTable &external_functions_table,
                            TrendComponentModelTable &trend_component_model_table_arg,
                            VarModelTable &var_model_table_arg);
  /*! \todo check there are only endogenous variables at the current period in the objective
    (no exogenous, no lead/lag) */
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  /*! \todo allow for the possibility of disabling temporary terms */
  void computingPass() override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
  //! Return a reference the Planner Objective model tree
  StaticModel &getPlannerObjective();
};

class BVARDensityStatement : public Statement
{
private:
  const int maxnlags;
  const OptionsList options_list;
public:
  BVARDensityStatement(int maxnlags_arg, OptionsList options_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class BVARForecastStatement : public Statement
{
private:
  const int nlags;
  const OptionsList options_list;
public:
  BVARForecastStatement(int nlags_arg, OptionsList options_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class SBVARStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  SBVARStatement(OptionsList options_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class MSSBVAREstimationStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  MSSBVAREstimationStatement(OptionsList options_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class MSSBVARSimulationStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  MSSBVARSimulationStatement(OptionsList options_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class MSSBVARComputeMDDStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  MSSBVARComputeMDDStatement(OptionsList options_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class MSSBVARComputeProbabilitiesStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  MSSBVARComputeProbabilitiesStatement(OptionsList options_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class MSSBVARIrfStatement : public Statement
{
private:
  const SymbolList symbol_list;
  const OptionsList options_list;
public:
  MSSBVARIrfStatement(SymbolList symbol_list_arg,
                      OptionsList options_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class MSSBVARForecastStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  MSSBVARForecastStatement(OptionsList options_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class MSSBVARVarianceDecompositionStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  MSSBVARVarianceDecompositionStatement(OptionsList options_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class IdentificationStatement : public Statement
{
private:
  OptionsList options_list;
public:
  IdentificationStatement(const OptionsList &options_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class WriteLatexDynamicModelStatement : public Statement
{
private:
  const DynamicModel &dynamic_model;
  const bool write_equation_tags;
public:
  WriteLatexDynamicModelStatement(const DynamicModel &dynamic_model_arg, bool write_equation_tags_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class WriteLatexStaticModelStatement : public Statement
{
private:
  const StaticModel &static_model;
  const bool write_equation_tags;
public:
  WriteLatexStaticModelStatement(const StaticModel &static_model_arg, bool write_equation_tags_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class WriteLatexOriginalModelStatement : public Statement
{
private:
  const DynamicModel &original_model;
  const bool write_equation_tags;
public:
  WriteLatexOriginalModelStatement(const DynamicModel &original_model_arg, bool write_equation_tags_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class WriteLatexSteadyStateModelStatement : public Statement
{
private:
  const SteadyStateModel &steady_state_model;
public:
  WriteLatexSteadyStateModelStatement(const SteadyStateModel &steady_state_model_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class ShockDecompositionStatement : public Statement
{
private:
  const SymbolList symbol_list;
  const OptionsList options_list;
public:
  ShockDecompositionStatement(SymbolList symbol_list_arg,
                              OptionsList options_list_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class RealtimeShockDecompositionStatement : public Statement
{
private:
  const SymbolList symbol_list;
  const OptionsList options_list;
public:
  RealtimeShockDecompositionStatement(SymbolList symbol_list_arg,
                                      OptionsList options_list_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
};

class PlotShockDecompositionStatement : public Statement
{
private:
  const SymbolList symbol_list;
  const OptionsList options_list;
public:
  PlotShockDecompositionStatement(SymbolList symbol_list_arg,
                                  OptionsList options_list_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
};

class InitialConditionDecompositionStatement : public Statement
{
private:
  const SymbolList symbol_list;
  const OptionsList options_list;
public:
  InitialConditionDecompositionStatement(SymbolList symbol_list_arg,
                                         OptionsList options_list_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
};

class ConditionalForecastStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  ConditionalForecastStatement(OptionsList options_list_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class PlotConditionalForecastStatement : public Statement
{
private:
  //! A value of -1 indicates that the user didn't specify a value
  const int periods;
  const SymbolList symbol_list;
public:
  PlotConditionalForecastStatement(int periods_arg, SymbolList symbol_list_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class CalibSmootherStatement : public Statement
{
private:
  const SymbolList symbol_list;
  const OptionsList options_list;
public:
  CalibSmootherStatement(SymbolList symbol_list_arg,
                         OptionsList options_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class ExtendedPathStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  ExtendedPathStatement(OptionsList options_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class SvarIdentificationStatement : public Statement
{
public:
  //  using svar_identification_exclusion_t = map<pair<int, int>, vector<int>>;
  struct svar_identification_restriction
  {
    int equation;
    int restriction_nbr;
    int lag;
    int variable;
    expr_t value;
  };

  using svar_identification_restrictions_t = vector<svar_identification_restriction>;
private:
  const svar_identification_restrictions_t restrictions;
  const bool upper_cholesky_present;
  const bool lower_cholesky_present;
  const bool constants_exclusion_present;
  const SymbolTable &symbol_table;
  int getMaxLag() const;
public:
  SvarIdentificationStatement(svar_identification_restrictions_t restrictions_arg,
                              const bool &upper_cholesky_present_arg,
                              const bool &lower_cholesky_present_arg,
                              const bool &constants_exclusion_present_arg,
                              const SymbolTable &symbol_table_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class MarkovSwitchingStatement : public Statement
{
private:
  const OptionsList options_list;
  map <pair<int, int >, double > restriction_map;
public:
  MarkovSwitchingStatement(OptionsList options_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeCOutput(ostream &output, const string &basename) override;
  void writeJsonOutput(ostream &output) const override;
};

class SvarStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  SvarStatement(OptionsList options_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class SvarGlobalIdentificationCheckStatement : public Statement
{
public:
  SvarGlobalIdentificationCheckStatement();
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class SetTimeStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  SetTimeStatement(OptionsList options_list_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class EstimationDataStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  EstimationDataStatement(OptionsList options_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class SubsamplesStatement : public Statement
{
public:
  //! Storage for declaring subsamples: map<subsample_name, <date1, date2 >
  using subsample_declaration_map_t = map<string, pair<string, string>>;
private:
  const string name1;
  const string name2;
  const subsample_declaration_map_t subsample_declaration_map;
  const SymbolTable &symbol_table;
public:
  SubsamplesStatement(string name1_arg,
                      string name2_arg,
                      subsample_declaration_map_t subsample_declaration_map_arg,
                      const SymbolTable &symbol_table_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class SubsamplesEqualStatement : public Statement
{
private:
  const string to_name1;
  const string to_name2;
  const string from_name1;
  const string from_name2;
  const SymbolTable &symbol_table;
public:
  SubsamplesEqualStatement(string to_name1_arg,
                           string to_name2_arg,
                           string from_name1_arg,
                           string from_name2_arg,
                           const SymbolTable &symbol_table_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class JointPriorStatement : public Statement
{
private:
  const vector<string> joint_parameters;
  const PriorDistributions prior_shape;
  const OptionsList options_list;
public:
  JointPriorStatement(vector<string> joint_parameters_arg,
                      PriorDistributions prior_shape_arg,
                      OptionsList options_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeOutputHelper(ostream &output, const string &field, const string &lhs_field) const;
  void writeJsonOutput(ostream &output) const override;
};

class BasicPriorStatement : public Statement
{
public:
  
  ~BasicPriorStatement() override;
protected:
  const string name;
  const string subsample_name;
  const PriorDistributions prior_shape;
  const expr_t variance;
  const OptionsList options_list;
  BasicPriorStatement(string name_arg,
                      string subsample_name_arg,
                      PriorDistributions prior_shape_arg,
                      expr_t variance_arg,
                      OptionsList options_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void get_base_name(const SymbolType symb_type, string &lhs_field) const;
  void writeCommonOutput(ostream &output, const string &lhs_field) const;
  void writeCommonOutputHelper(ostream &output, const string &field, const string &lhs_field) const;
  void writePriorOutput(ostream &output, string &lhs_field, const string &name2) const;
  bool is_structural_innovation(const SymbolType symb_type) const;
  void writePriorIndex(ostream &output, const string &lhs_field) const;
  void writeVarianceOption(ostream &output, const string &lhs_field) const;
  void writeOutputHelper(ostream &output, const string &field, const string &lhs_field) const;
  void writeShape(ostream &output, const string &lhs_field) const;
  void writeCOutputHelper(ostream &output, const string &field) const;
  void writeCShape(ostream &output) const;
  void writeCVarianceOption(ostream &output) const;
  void writeCDomain(ostream &output) const;
  void writeJsonShape(ostream &output) const;
  void writeJsonPriorOutput(ostream &output) const;
};

class PriorStatement : public BasicPriorStatement
{
public:
  PriorStatement(string name_arg,
                 string subsample_name_arg,
                 PriorDistributions prior_shape_arg,
                 expr_t variance_arg,
                 OptionsList options_list_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeCOutput(ostream &output, const string &basename) override;
  void writeJsonOutput(ostream &output) const override;
};

class StdPriorStatement : public BasicPriorStatement
{
private:
  const SymbolTable &symbol_table;
public:
  StdPriorStatement(string name_arg,
                    string subsample_name_arg,
                    PriorDistributions prior_shape_arg,
                    expr_t variance_arg,
                    OptionsList options_list_arg,
                    const SymbolTable &symbol_table_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeCOutput(ostream &output, const string &basename) override;
  void writeJsonOutput(ostream &output) const override;
};

class CorrPriorStatement : public BasicPriorStatement
{
private:
  const string name1;
  const SymbolTable &symbol_table;
public:
  CorrPriorStatement(string name_arg1,
                     string name_arg2,
                     string subsample_name_arg,
                     PriorDistributions prior_shape_arg,
                     expr_t variance_arg,
                     OptionsList options_list_arg,
                     const SymbolTable &symbol_table_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeCOutput(ostream &output, const string &basename) override;
  void writeJsonOutput(ostream &output) const override;
};

class PriorEqualStatement : public Statement
{
private:
  const string to_declaration_type;
  const string to_name1;
  const string to_name2;
  const string to_subsample_name;
  const string from_declaration_type;
  const string from_name1;
  const string from_name2;
  const string from_subsample_name;
  const SymbolTable &symbol_table;
public:
  PriorEqualStatement(string to_declaration_type_arg,
                      string to_name1_arg,
                      string to_name2_arg,
                      string to_subsample_name_arg,
                      string from_declaration_type_arg,
                      string from_name1_arg,
                      string from_name2_arg,
                      string from_subsample_name_arg,
                      const SymbolTable &symbol_table_arg);
  void get_base_name(const SymbolType symb_type, string &lhs_field) const;
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class BasicOptionsStatement : public Statement
{
public:
  
  ~BasicOptionsStatement() override;
protected:
  const string name;
  const string subsample_name;
  const OptionsList options_list;
  BasicOptionsStatement(string name_arg,
                        string subsample_name_arg,
                        OptionsList options_list_arg);
  void get_base_name(const SymbolType symb_type, string &lhs_field) const;
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOptionsOutput(ostream &output, string &lhs_field, const string &name2) const;
  void writeCommonOutput(ostream &output, const string &lhs_field) const;
  void writeCommonOutputHelper(ostream &output, const string &field, const string &lhs_field) const;
  bool is_structural_innovation(const SymbolType symb_type) const;
  void writeOptionsIndex(ostream &output, const string &lhs_field) const;
  void writeOutputHelper(ostream &output, const string &field, const string &lhs_field) const;
  void writeCOutputHelper(ostream &output, const string &field) const;
  void writeJsonOptionsOutput(ostream &output) const;
};

class OptionsStatement : public BasicOptionsStatement
{
public:
  OptionsStatement(string name_arg, string subsample_name_arg, OptionsList options_list_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeCOutput(ostream &output, const string &basename) override;
  void writeJsonOutput(ostream &output) const override;
};

class StdOptionsStatement : public BasicOptionsStatement
{
private:
  const SymbolTable &symbol_table;
public:
  StdOptionsStatement(string name_arg,
                      string subsample_name_arg,
                      OptionsList options_list_arg,
                      const SymbolTable &symbol_table_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeCOutput(ostream &output, const string &basename) override;
  void writeJsonOutput(ostream &output) const override;
};

class CorrOptionsStatement : public BasicOptionsStatement
{
private:
  const string name1;
  const SymbolTable &symbol_table;
public:
  CorrOptionsStatement(string name_arg1, string name_arg2,
                       string subsample_name_arg,
                       OptionsList options_list_arg,
                       const SymbolTable &symbol_table_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeCOutput(ostream &output, const string &basename) override;
  void writeJsonOutput(ostream &output) const override;
};

class OptionsEqualStatement : public Statement
{
private:
  const string to_declaration_type;
  const string to_name1;
  const string to_name2;
  const string to_subsample_name;
  const string from_declaration_type;
  const string from_name1;
  const string from_name2;
  const string from_subsample_name;
  const SymbolTable &symbol_table;
public:
  OptionsEqualStatement(string to_declaration_type_arg,
                        string to_name1_arg,
                        string to_name2_arg,
                        string to_subsample_name_arg,
                        string from_declaration_type_arg,
                        string from_name1_arg,
                        string from_name2_arg,
                        string from_subsample_name_arg,
                        const SymbolTable &symbol_table_arg);
  void get_base_name(const SymbolType symb_type, string &lhs_field) const;
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class ModelDiagnosticsStatement : public Statement
{
public:
  ModelDiagnosticsStatement();
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class Smoother2histvalStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  Smoother2histvalStatement(OptionsList options_list_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class GMMEstimationStatement : public Statement
{
private:
  const SymbolList symbol_list;
  const OptionsList options_list;
public:
  GMMEstimationStatement(SymbolList symbol_list_arg, OptionsList options_list_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class SMMEstimationStatement : public Statement
{
private:
  const SymbolList symbol_list;
  const OptionsList options_list;
public:
  SMMEstimationStatement(SymbolList symbol_list_arg, OptionsList options_list_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class GenerateIRFsStatement : public Statement
{
public:
private:
  const OptionsList options_list;
  const vector<string> generate_irf_names;
  const vector<map<string, double>> generate_irf_elements;
public:
  GenerateIRFsStatement(OptionsList options_list_arg,
                        vector<string> generate_irf_names_arg,
                        vector<map<string, double>> generate_irf_elements_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class VarExpectationModelStatement : public Statement
{
public:
  const string model_name, variable, var_model_name, horizon;
  const expr_t discount;
  const SymbolTable &symbol_table;
  // List of generated auxiliary param ids, in variable-major order
  vector<int> aux_params_ids; // TODO: move this to some new VarModelTable object

public:
  VarExpectationModelStatement(string model_name_arg, string variable_arg, string var_model_name_arg,
                               string horizon_arg, expr_t discount_arg, const SymbolTable &symbol_table_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

#endif
