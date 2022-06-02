/*
 * Copyright © 2003-2022 Dynare Team
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

#ifndef _COMPUTINGTASKS_HH
#define _COMPUTINGTASKS_HH

#include <ostream>
#include <optional>

#include "SymbolList.hh"
#include "SymbolTable.hh"
#include "Statement.hh"
#include "StaticModel.hh"
#include "DynamicModel.hh"
#include "ModelEquationBlock.hh"

class SteadyStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  explicit SteadyStatement(OptionsList options_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class CheckStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  explicit CheckStatement(OptionsList options_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class SimulStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  explicit SimulStatement(OptionsList options_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class PerfectForesightSetupStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  explicit PerfectForesightSetupStatement(OptionsList options_list_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class PerfectForesightSolverStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  explicit PerfectForesightSolverStatement(OptionsList options_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class PerfectForesightWithExpectationErrorsSetupStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  explicit PerfectForesightWithExpectationErrorsSetupStatement(OptionsList options_list_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class PerfectForesightWithExpectationErrorsSolverStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  explicit PerfectForesightWithExpectationErrorsSolverStatement(OptionsList options_list_arg);
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
  explicit ModelInfoStatement(OptionsList options_list_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class StochSimulStatement : public Statement
{
private:
  const SymbolList symbol_list;
  const OptionsList options_list;
  const SymbolTable &symbol_table;
public:
  StochSimulStatement(SymbolList symbol_list_arg, OptionsList options_list_arg,
                      const SymbolTable &symbol_table_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class ForecastStatement : public Statement
{
private:
  const SymbolList symbol_list;
  const OptionsList options_list;
  const SymbolTable &symbol_table;
public:
  ForecastStatement(SymbolList symbol_list_arg, OptionsList options_list_arg,
                    const SymbolTable &symbol_table_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class RamseyModelStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  explicit RamseyModelStatement(OptionsList options_list_arg);
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
  const SymbolList symbol_list;
  const OptionsList options_list;
  const SymbolTable &symbol_table;
public:
  RamseyPolicyStatement(SymbolList symbol_list_arg, OptionsList options_list_arg,
                        const SymbolTable &symbol_table_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class EvaluatePlannerObjectiveStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  explicit EvaluatePlannerObjectiveStatement(OptionsList options_list_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class OccbinSetupStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  explicit OccbinSetupStatement(OptionsList options_list_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class OccbinSolverStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  explicit OccbinSolverStatement(OptionsList options_list_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class OccbinWriteRegimesStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  explicit OccbinWriteRegimesStatement(OptionsList options_list_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class OccbinGraphStatement : public Statement
{
private:
  const SymbolList symbol_list;
  const OptionsList options_list;
public:
  OccbinGraphStatement(SymbolList symbol_list_arg,
                       OptionsList options_list_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class DiscretionaryPolicyStatement : public Statement
{
private:
  const SymbolList symbol_list;
  const OptionsList options_list;
  const SymbolTable &symbol_table;
public:
  DiscretionaryPolicyStatement(SymbolList symbol_list_arg, OptionsList options_list_arg,
                               const SymbolTable &symbol_table_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class RplotStatement : public Statement
{
private:
  const SymbolList symbol_list;
  const SymbolTable &symbol_table;
public:
  RplotStatement(SymbolList symbol_list_arg, const SymbolTable &symbol_table_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class UnitRootVarsStatement : public Statement
{
public:
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class PeriodsStatement : public Statement
{
private:
  const int periods;
public:
  explicit PeriodsStatement(int periods_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class DsampleStatement : public Statement
{
private:
  const int val1, val2;
public:
  explicit DsampleStatement(int val1_arg);
  DsampleStatement(int val1_arg, int val2_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class EstimationStatement : public Statement
{
private:
  const SymbolTable &symbol_table;
  const SymbolList symbol_list;
  const OptionsList options_list;
public:
  EstimationStatement(const SymbolTable &symbol_table_arg,
                      SymbolList symbol_list_arg,
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
  explicit DynareSensitivityStatement(OptionsList options_list_arg);
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

class DeterministicTrendsStatement : public Statement
{
public:
  using trend_elements_t = map<string, expr_t>;
private:
  const trend_elements_t trend_elements;
  const SymbolTable &symbol_table;
public:
  DeterministicTrendsStatement(trend_elements_t trend_elements_arg,
                             const SymbolTable &symbol_table_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class FilterInitialStateStatement : public Statement
{
public:
  using filter_initial_state_elements_t = map<pair<int, int>, expr_t>;
private:
  const filter_initial_state_elements_t filter_initial_state_elements;
  const SymbolTable &symbol_table;
public:
  FilterInitialStateStatement(filter_initial_state_elements_t filter_initial_state_elements_arg,
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
  const SymbolTable &symbol_table;
public:
  OsrStatement(SymbolList symbol_list_arg, OptionsList options_list_arg,
               const SymbolTable &symbol_table_arg);
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
  explicit OsrParamsBoundsStatement(vector<OsrParams> osr_params_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class DynaTypeStatement : public Statement
{
private:
  const SymbolList symbol_list;
  const string filename;
  const SymbolTable &symbol_table;
public:
  DynaTypeStatement(SymbolList symbol_list_arg, string filename_arg,
                    const SymbolTable &symbol_table_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class DynaSaveStatement : public Statement
{
private:
  const SymbolList symbol_list;
  const string filename;
  const SymbolTable &symbol_table;
public:
  DynaSaveStatement(SymbolList symbol_list_arg, string filename_arg,
                    const SymbolTable &symbol_table_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
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

class AbstractEstimatedParamsStatement : public Statement
{
protected:
  const vector<EstimationParams> estim_params_list;
  const SymbolTable &symbol_table;
  AbstractEstimatedParamsStatement(vector<EstimationParams> estim_params_list_arg,
                                   const SymbolTable &symbol_table_arg);
  virtual string blockName() const = 0;
  // Part of the check pass that is common to the three estimated_params{,_init,bounds} blocks
  void commonCheckPass() const;
};

class EstimatedParamsStatement : public AbstractEstimatedParamsStatement
{
private:
  const bool overwrite;
public:
  EstimatedParamsStatement(vector<EstimationParams> estim_params_list_arg,
                           const SymbolTable &symbol_table_arg,
                           bool overwrite_arg);
  string blockName() const override { return "estimated_params"; };
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class EstimatedParamsInitStatement : public AbstractEstimatedParamsStatement
{
private:
  const bool use_calibration;
public:
  EstimatedParamsInitStatement(vector<EstimationParams> estim_params_list_arg,
                               const SymbolTable &symbol_table_arg,
                               const bool use_calibration_arg);
  string blockName() const override { return "estimated_params_init"; };
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class EstimatedParamsBoundsStatement : public AbstractEstimatedParamsStatement
{
public:
  EstimatedParamsBoundsStatement(vector<EstimationParams> estim_params_list_arg,
                                 const SymbolTable &symbol_table_arg);
  string blockName() const override { return "estimated_params_bounds"; };
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class EstimatedParamsRemoveStatement : public Statement
{
public:
  // Only the type, name and name2 fields of EstimationParams are used
  const vector<EstimationParams> estim_params_list;
  const SymbolTable &symbol_table;
  EstimatedParamsRemoveStatement(vector<EstimationParams> estim_params_list_arg,
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
  explicit PlannerObjectiveStatement(const StaticModel &model_tree_arg);
  /*! \todo check there are only endogenous variables at the current period in the objective
    (no exogenous, no lead/lag) */
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  /*! \todo allow for the possibility of disabling temporary terms */
  void computingPass(const ModFileStructure &mod_file_struct) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
  //! Return a reference the Planner Objective model tree
  const StaticModel &getPlannerObjective() const;
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
  explicit SBVARStatement(OptionsList options_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class MSSBVAREstimationStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  explicit MSSBVAREstimationStatement(OptionsList options_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class MSSBVARSimulationStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  explicit MSSBVARSimulationStatement(OptionsList options_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class MSSBVARComputeMDDStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  explicit MSSBVARComputeMDDStatement(OptionsList options_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class MSSBVARComputeProbabilitiesStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  explicit MSSBVARComputeProbabilitiesStatement(OptionsList options_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class MSSBVARIrfStatement : public Statement
{
private:
  const SymbolList symbol_list;
  const OptionsList options_list;
  const SymbolTable &symbol_table;
public:
  MSSBVARIrfStatement(SymbolList symbol_list_arg, OptionsList options_list_arg,
                      const SymbolTable &symbol_table_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class MSSBVARForecastStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  explicit MSSBVARForecastStatement(OptionsList options_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class MSSBVARVarianceDecompositionStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  explicit MSSBVARVarianceDecompositionStatement(OptionsList options_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class IdentificationStatement : public Statement
{
private:
  OptionsList options_list;
public:
  explicit IdentificationStatement(OptionsList options_list_arg);
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
  explicit WriteLatexSteadyStateModelStatement(const SteadyStateModel &steady_state_model_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class ShockDecompositionStatement : public Statement
{
private:
  const SymbolList symbol_list;
  const OptionsList options_list;
  const SymbolTable &symbol_table;
public:
  ShockDecompositionStatement(SymbolList symbol_list_arg, OptionsList options_list_arg,
                              const SymbolTable &symbol_table_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class RealtimeShockDecompositionStatement : public Statement
{
private:
  const SymbolList symbol_list;
  const OptionsList options_list;
  const SymbolTable &symbol_table;
public:
  RealtimeShockDecompositionStatement(SymbolList symbol_list_arg, OptionsList options_list_arg,
                                      const SymbolTable &symbol_table_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class PlotShockDecompositionStatement : public Statement
{
private:
  const SymbolList symbol_list;
  const OptionsList options_list;
  const SymbolTable &symbol_table;
public:
  PlotShockDecompositionStatement(SymbolList symbol_list_arg, OptionsList options_list_arg,
                                  const SymbolTable &symbol_table_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class InitialConditionDecompositionStatement : public Statement
{
private:
  const SymbolList symbol_list;
  const OptionsList options_list;
  const SymbolTable &symbol_table;
public:
  InitialConditionDecompositionStatement(SymbolList symbol_list_arg, OptionsList options_list_arg,
                                         const SymbolTable &symbol_table_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class SqueezeShockDecompositionStatement : public Statement
{
private:
  const SymbolList symbol_list;
  const SymbolTable &symbol_table;
public:
  SqueezeShockDecompositionStatement(SymbolList symbol_list_arg,
                                     const SymbolTable &symbol_table_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class ConditionalForecastStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  explicit ConditionalForecastStatement(OptionsList options_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class PlotConditionalForecastStatement : public Statement
{
private:
  const optional<int> periods; // The user is allowed not to declare periods
  const SymbolList symbol_list;
  const SymbolTable &symbol_table;
public:
  PlotConditionalForecastStatement(optional<int> periods_arg, SymbolList symbol_list_arg,
                                   const SymbolTable &symbol_table_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class CalibSmootherStatement : public Statement
{
private:
  const SymbolList symbol_list;
  const OptionsList options_list;
  const SymbolTable &symbol_table;
public:
  CalibSmootherStatement(SymbolList symbol_list_arg, OptionsList options_list_arg,
                         const SymbolTable &symbol_table_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class ExtendedPathStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  explicit ExtendedPathStatement(OptionsList options_list_arg);
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
  const bool upper_cholesky_present, lower_cholesky_present, constants_exclusion_present;
  const SymbolTable &symbol_table;
  int getMaxLag() const;
public:
  SvarIdentificationStatement(svar_identification_restrictions_t restrictions_arg,
                              bool upper_cholesky_present_arg,
                              bool lower_cholesky_present_arg,
                              bool constants_exclusion_present_arg,
                              const SymbolTable &symbol_table_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class MarkovSwitchingStatement : public Statement
{
private:
  const OptionsList options_list;
  map<pair<int, int>, double> restriction_map;
public:
  explicit MarkovSwitchingStatement(OptionsList options_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class SvarStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  explicit SvarStatement(OptionsList options_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class SvarGlobalIdentificationCheckStatement : public Statement
{
public:
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class SetTimeStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  explicit SetTimeStatement(OptionsList options_list_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class EstimationDataStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  explicit EstimationDataStatement(OptionsList options_list_arg);
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
  const string name1, name2;
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
  const string to_name1, to_name2, from_name1, from_name2;
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
protected:
  const string name, subsample_name;
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
  void writeJsonOutput(ostream &output) const override;
};

class PriorEqualStatement : public Statement
{
private:
  const string to_declaration_type, to_name1, to_name2, to_subsample_name;
  const string from_declaration_type, from_name1, from_name2, from_subsample_name;
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
protected:
  const string name, subsample_name;
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
  void writeJsonOptionsOutput(ostream &output) const;
};

class OptionsStatement : public BasicOptionsStatement
{
public:
  OptionsStatement(string name_arg, string subsample_name_arg, OptionsList options_list_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
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
  void writeJsonOutput(ostream &output) const override;
};

class OptionsEqualStatement : public Statement
{
private:
  const string to_declaration_type, to_name1, to_name2, to_subsample_name;
  const string from_declaration_type, from_name1, from_name2, from_subsample_name;
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
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class Smoother2histvalStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  explicit Smoother2histvalStatement(OptionsList options_list_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class MethodOfMomentsStatement : public Statement
{
private:  
  const OptionsList options_list;
public:
  explicit MethodOfMomentsStatement(OptionsList options_list_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
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

class MatchedMomentsStatement : public Statement
{
private:
  const SymbolTable &symbol_table;
public:
  /* Each moment is represented by a three vectors: symb_ids, lags, powers.
     See the definition of ExprNode::matchMatchedMoment() for more details */
  const vector<tuple<vector<int>, vector<int>, vector<int>>> moments;
  MatchedMomentsStatement(const SymbolTable &symbol_table_arg,
                          vector<tuple<vector<int>, vector<int>, vector<int>>> moments_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class OccbinConstraintsStatement : public Statement
{
private:
  DataTree data_tree;
public:
  // The tuple is (name, bind, relax, error_bind, error_relax) (where relax and error_{bind,relax} can be nullptr)
  const vector<tuple<string, BinaryOpNode *, BinaryOpNode *, expr_t, expr_t>> constraints;
  OccbinConstraintsStatement(const DataTree &data_tree_arg,
                             const vector<tuple<string, BinaryOpNode *, BinaryOpNode *, expr_t, expr_t>> constraints_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class ResidStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  explicit ResidStatement(OptionsList options_list_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

#endif
