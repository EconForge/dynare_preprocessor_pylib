/*
 * Copyright © 2003-2024 Dynare Team
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

#ifndef SHOCKS_HH
#define SHOCKS_HH

#include <map>
#include <string>
#include <vector>

#include "ExprNode.hh"
#include "HeterogeneityTable.hh"
#include "Statement.hh"
#include "SymbolTable.hh"

using namespace std;

class AbstractShocksStatement : public Statement
{
public:
  // The tuple is (period1, period2, value)
  using det_shocks_t = map<int, vector<tuple<int, int, expr_t>>>;
  enum class ShockType
  {
    level, // The value is the level of the exogenous (“values” statement in “shocks”)
    multiplySteadyState, // The value is the ratio of the exogenous over its (terminal) steady state
                         // (“values” statement in “mshocks”)
    multiplyInitialSteadyState // The value is the ratio of the exogenous over its initial steady
                               // state (“values” statement in “mshocks(relative_to_initval)”)
  };

protected:
  //! Does this "shocks" statement replace the previous ones?
  const bool overwrite;
  const ShockType type; // Type of shocks represented by this block
  const det_shocks_t det_shocks;
  const SymbolTable& symbol_table;
  void writeDetShocks(ostream& output) const;
  void writeJsonDetShocks(ostream& output) const;
  static string typeToString(ShockType type);

  AbstractShocksStatement(bool overwrite_arg, ShockType type_arg, det_shocks_t det_shocks_arg,
                          const SymbolTable& symbol_table_arg);
};

class ShocksStatement : public AbstractShocksStatement
{
public:
  using var_and_std_shocks_t = map<int, expr_t>;
  using covar_and_corr_shocks_t = map<pair<int, int>, expr_t>;

private:
  const var_and_std_shocks_t var_shocks, std_shocks;
  const covar_and_corr_shocks_t covar_shocks, corr_shocks;
  void writeVarOrStdShock(ostream& output, const pair<int, expr_t>& it, bool stddev) const;
  void writeVarAndStdShocks(ostream& output) const;
  void writeCovarOrCorrShock(ostream& output, const pair<pair<int, int>, expr_t>& it,
                             bool corr) const;
  void writeCovarAndCorrShocks(ostream& output) const;
  [[nodiscard]] bool has_calibrated_measurement_errors() const;

public:
  ShocksStatement(bool overwrite_arg, det_shocks_t det_shocks_arg,
                  var_and_std_shocks_t var_shocks_arg, var_and_std_shocks_t std_shocks_arg,
                  covar_and_corr_shocks_t covar_shocks_arg, covar_and_corr_shocks_t corr_shocks_arg,
                  const SymbolTable& symbol_table_arg);
  void checkPass(ModFileStructure& mod_file_struct, WarningConsolidation& warnings) override;
  void writeOutput(ostream& output, const string& basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream& output) const override;
};

class MShocksStatement : public AbstractShocksStatement
{
public:
  const bool relative_to_initval;
  MShocksStatement(bool overwrite_arg, bool relative_to_initval_arg, det_shocks_t det_shocks_arg,
                   const SymbolTable& symbol_table_arg);
  void writeOutput(ostream& output, const string& basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream& output) const override;
};

/* Represents a shocks(surprise) block.
   Given the differences with the plain “shocks” block, it was easier to make
   it a separate class. */
class ShocksSurpriseStatement : public Statement
{
public:
  //! Does this "shocks(surprise)" statement replace the previous ones?
  const bool overwrite;
  const AbstractShocksStatement::det_shocks_t surprise_shocks;

private:
  const SymbolTable& symbol_table;

public:
  ShocksSurpriseStatement(bool overwrite_arg,
                          AbstractShocksStatement::det_shocks_t surprise_shocks_arg,
                          const SymbolTable& symbol_table_arg);
  void checkPass(ModFileStructure& mod_file_struct, WarningConsolidation& warnings) override;
  void writeOutput(ostream& output, const string& basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream& output) const override;
};

/* Represents a shocks(learnt_in=…) or mshocks(learnt_in=…) block.
   Given the differences with the plain “shocks” and “mshocks” blocks,
   it was easier to make it a separate class. */
class ShocksLearntInStatement : public Statement
{
public:
  const int learnt_in_period;
  //! Does this “shocks(learnt_in=…)” or “mshocks(learnt_in=…)” block replace the previous ones?
  const bool overwrite;
  enum class LearntShockType
  {
    level, // The value is the level of the exogenous (“values” statement in “shocks(learnt_in=…)”)
    add,      // The value is the additive change of the exogenous compared to previous information
              // period (“add” statement in “shocks(learnt_in=…)”)
    multiply, // The value is the multiplicative change of the exogenous compared to previous
              // information period (“multiply” statement in “shocks(learnt_in=…)”)
    multiplySteadyState, // The value is the ratio of the exogenous over its (terminal) steady state
                         // as anticipated in the same informational period (“values” statement in
                         // “mshocks(learnt_in=…)”)
    multiplyInitialSteadyState // The value is the ratio of the exogenous over its initial steady
                               // state as anticipated in the same informational period (“values”
                               // statement in “mshocks(learnt_in=…, relative_to_initval)”)
  };
  // The tuple is (type, period1, period2, value)
  using learnt_shocks_t = map<int, vector<tuple<LearntShockType, int, int, expr_t>>>;
  const learnt_shocks_t learnt_shocks;

private:
  const SymbolTable& symbol_table;
  static string typeToString(LearntShockType type);

public:
  ShocksLearntInStatement(int learnt_in_period_arg, bool overwrite_arg,
                          learnt_shocks_t learnt_shocks_arg, const SymbolTable& symbol_table_arg);
  void checkPass(ModFileStructure& mod_file_struct, WarningConsolidation& warnings) override;
  void writeOutput(ostream& output, const string& basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream& output) const override;
};

class HeterogeneousShocksStatement : public Statement
{
public:
  const int heterogeneity_dimension;
  const bool overwrite;

  using var_and_std_shocks_t = map<int, expr_t>;
  using covar_and_corr_shocks_t = map<pair<int, int>, expr_t>;

  const var_and_std_shocks_t var_shocks, std_shocks;
  const covar_and_corr_shocks_t covar_shocks, corr_shocks;

private:
  const SymbolTable& symbol_table;
  const HeterogeneityTable& heterogeneity_table;

  void writeVarOrStdShock(ostream& output, const pair<int, expr_t>& it, bool stddev) const;
  void writeVarAndStdShocks(ostream& output) const;
  void writeCovarOrCorrShock(ostream& output, const pair<pair<int, int>, expr_t>& it,
                             bool corr) const;
  void writeCovarAndCorrShocks(ostream& output) const;

  string
  sigmaeName() const
  {
    return "M_.heterogeneity("s + to_string(heterogeneity_dimension + 1) + ").Sigma_e"s;
  }

public:
  HeterogeneousShocksStatement(int heterogeneity_dimension_arg, bool overwrite_arg,
                               var_and_std_shocks_t var_shocks_arg,
                               var_and_std_shocks_t std_shocks_arg,
                               covar_and_corr_shocks_t covar_shocks_arg,
                               covar_and_corr_shocks_t corr_shocks_arg,
                               const SymbolTable& symbol_table_arg,
                               const HeterogeneityTable& heterogeneity_table_arg);
  void checkPass(ModFileStructure& mod_file_struct, WarningConsolidation& warnings) override;
  void writeOutput(ostream& output, const string& basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream& output) const override;
};

class ConditionalForecastPathsStatement : public Statement
{
private:
  const AbstractShocksStatement::det_shocks_t paths;
  const SymbolTable& symbol_table;
  const int path_length;

public:
  ConditionalForecastPathsStatement(AbstractShocksStatement::det_shocks_t paths_arg,
                                    const SymbolTable& symbol_table_arg);
  void writeOutput(ostream& output, const string& basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream& output) const override;
  static int computePathLength(const AbstractShocksStatement::det_shocks_t& paths);
};

class MomentCalibration : public Statement
{
public:
  struct Constraint
  {
    int endo1, endo2;
    string lags;
    expr_t lower_bound, upper_bound;
  };
  using constraints_t = vector<Constraint>;

private:
  constraints_t constraints;
  const SymbolTable& symbol_table;

public:
  MomentCalibration(constraints_t constraints_arg, const SymbolTable& symbol_table_arg);
  void writeOutput(ostream& output, const string& basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream& output) const override;
};

class IrfCalibration : public Statement
{
public:
  struct Constraint
  {
    int endo;
    int exo;
    string periods;
    expr_t lower_bound, upper_bound;
  };
  using constraints_t = vector<Constraint>;

private:
  constraints_t constraints;
  const SymbolTable& symbol_table;
  const OptionsList options_list;

public:
  IrfCalibration(constraints_t constraints_arg, const SymbolTable& symbol_table_arg,
                 OptionsList options_list_arg);
  void writeOutput(ostream& output, const string& basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream& output) const override;
};

class ShockGroupsStatement : public Statement
{
public:
  struct Group
  {
    string name;
    vector<string> list;
  };
  using group_t = vector<Group>;

private:
  group_t shock_groups;
  string name;

public:
  ShockGroupsStatement(group_t shock_groups_arg, string name_arg);
  void writeOutput(ostream& output, const string& basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream& output) const override;
};

class Init2shocksStatement : public Statement
{
private:
  const vector<pair<int, int>> init2shocks;
  const string name;
  const SymbolTable& symbol_table;

public:
  Init2shocksStatement(vector<pair<int, int>> init2shocks_arg, string name_arg,
                       const SymbolTable& symbol_table_arg);
  void checkPass(ModFileStructure& mod_file_struct, WarningConsolidation& warnings) override;
  void writeOutput(ostream& output, const string& basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream& output) const override;
};

class HeteroskedasticShocksStatement : public Statement
{
public:
  // Maps exo symb_id to list of tuples (period1, period2, value/scale)
  using heteroskedastic_shocks_t = map<int, vector<tuple<int, int, expr_t>>>;

private:
  const bool overwrite;
  const heteroskedastic_shocks_t values, scales;
  const SymbolTable& symbol_table;

public:
  HeteroskedasticShocksStatement(bool overwrite_arg, heteroskedastic_shocks_t values_arg,
                                 heteroskedastic_shocks_t scales_arg,
                                 const SymbolTable& symbol_table_arg);
  void writeOutput(ostream& output, const string& basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream& output) const override;
};

#endif
