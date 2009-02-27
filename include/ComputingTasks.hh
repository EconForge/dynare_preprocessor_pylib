/*
 * Copyright (C) 2003-2009 Dynare Team
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
#include "ModelTree.hh"

class SteadyStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  SteadyStatement(const OptionsList &options_list_arg);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

class SteadySparseStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  SteadySparseStatement(const OptionsList &options_list_arg);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

class CheckStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  CheckStatement(const OptionsList &options_list_arg);
  virtual void checkPass(ModFileStructure &mod_file_struct);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

class SimulStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  SimulStatement(const OptionsList &options_list_arg);
  virtual void checkPass(ModFileStructure &mod_file_struct);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

class SimulSparseStatement : public Statement
{
private:
  const OptionsList options_list;
  const int mode;
public:
  SimulSparseStatement(const OptionsList &options_list_arg, int mode_arg);
  virtual void checkPass(ModFileStructure &mod_file_struct);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

class ModelInfoStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  ModelInfoStatement(const OptionsList &options_list_arg);
  virtual void checkPass(ModFileStructure &mod_file_struct);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

class StochSimulStatement : public Statement
{
private:
  const SymbolList symbol_list;
  const OptionsList options_list;
public:
  StochSimulStatement(const SymbolList &symbol_list_arg,
                      const OptionsList &options_list_arg);
  virtual void checkPass(ModFileStructure &mod_file_struct);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

class ForecastStatement : public Statement
{
private:
  const SymbolList symbol_list;
  const OptionsList options_list;
public:
  ForecastStatement(const SymbolList &symbol_list_arg,
                    const OptionsList &options_list_arg);
  virtual void checkPass(ModFileStructure &mod_file_struct);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

class RamseyPolicyStatement : public Statement
{
private:
  const SymbolList symbol_list;
  const OptionsList options_list;
public:
  RamseyPolicyStatement(const SymbolList &symbol_list_arg,
                        const OptionsList &options_list_arg);
  virtual void checkPass(ModFileStructure &mod_file_struct);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

class RplotStatement : public Statement
{
private:
  const SymbolList symbol_list;
  const OptionsList options_list;
public:
  RplotStatement(const SymbolList &symbol_list_arg,
                 const OptionsList &options_list_arg);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

class UnitRootVarsStatement : public Statement
{
private:
  const SymbolList symbol_list;
public:
  UnitRootVarsStatement(const SymbolList &symbol_list_arg);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

class PeriodsStatement : public Statement
{
private:
  const int periods;
public:
  PeriodsStatement(int periods_arg);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

class CutoffStatement : public Statement
{
private:
  const double cutoff;
public:
  CutoffStatement(double cutoff_arg);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

class MarkowitzStatement : public Statement
{
private:
  const double markowitz;
public:
  MarkowitzStatement(double markowitz_arg);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

class DsampleStatement : public Statement
{
private:
  const int val1, val2;
public:
  DsampleStatement(int val1_arg);
  DsampleStatement(int val1_arg, int val2_arg);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

class EstimationStatement : public Statement
{
private:
  const SymbolList symbol_list;
  const OptionsList options_list;
public:
  EstimationStatement(const SymbolList &symbol_list_arg,
                      const OptionsList &options_list_arg);
  virtual void checkPass(ModFileStructure &mod_file_struct);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

class PriorAnalysisStatement : public Statement
{
private:
  const SymbolList symbol_list;
  const OptionsList options_list;
public:
  PriorAnalysisStatement(const SymbolList &symbol_list_arg,
                         const OptionsList &options_list_arg);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

class PosteriorAnalysisStatement : public Statement
{
private:
  const SymbolList symbol_list;
  const OptionsList options_list;
public:
  PosteriorAnalysisStatement(const SymbolList &symbol_list_arg,
                             const OptionsList &options_list_arg);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

class DynareSensitivityStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  DynareSensitivityStatement(const OptionsList &options_list_arg);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

class VarobsStatement : public Statement
{
private:
  const SymbolList symbol_list;
public:
  VarobsStatement(const SymbolList &symbol_list_arg);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

class ObservationTrendsStatement : public Statement
{
public:
  typedef map<string, NodeID> trend_elements_type;
private:
  const trend_elements_type trend_elements;
  const SymbolTable &symbol_table;
public:
  ObservationTrendsStatement(const trend_elements_type &trend_elements_arg,
                             const SymbolTable &symbol_table_arg);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

class OsrParamsStatement : public Statement
{
private:
  const SymbolList symbol_list;
public:
  OsrParamsStatement(const SymbolList &symbol_list_arg);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

class OsrStatement : public Statement
{
private:
  const SymbolList symbol_list;
  const OptionsList options_list;
public:
  OsrStatement(const SymbolList &symbol_list_arg,
               const OptionsList &options_list_arg);
  virtual void checkPass(ModFileStructure &mod_file_struct);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

class DynaTypeStatement : public Statement
{
private:
  const SymbolList symbol_list;
  const string filename;
public:
  DynaTypeStatement(const SymbolList &symbol_list_arg,
                    const string &filename_arg);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

class DynaSaveStatement : public Statement
{
private:
  const SymbolList symbol_list;
  const string filename;
public:
  DynaSaveStatement(const SymbolList &symbol_list_arg,
                    const string &filename_arg);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

class ModelComparisonStatement : public Statement
{
public:
  typedef vector<pair<string, string> > filename_list_type;
private:
  filename_list_type filename_list;
  OptionsList options_list;
public:
  ModelComparisonStatement(const filename_list_type &filename_list_arg,
                           const OptionsList &options_list_arg);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

//! Temporary structure used when parsing estimation_params* statements
class EstimationParams
{
public:
  int type;
  string name, name2, prior;
  NodeID init_val, low_bound, up_bound, mean, std, p3, p4, jscale;

  void init(const DataTree &datatree)
  {
    type = 0;
    name = "";
    name2 = "";
    prior = "NaN";
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
  EstimatedParamsStatement(const vector<EstimationParams> &estim_params_list_arg,
                           const SymbolTable &symbol_table_arg);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

class EstimatedParamsInitStatement : public Statement
{
private:
  const vector<EstimationParams> estim_params_list;
  const SymbolTable &symbol_table;
public:
  EstimatedParamsInitStatement(const vector<EstimationParams> &estim_params_list_arg,
                               const SymbolTable &symbol_table_arg);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

class EstimatedParamsBoundsStatement : public Statement
{
private:
  const vector<EstimationParams> estim_params_list;
  const SymbolTable &symbol_table;
public:
  EstimatedParamsBoundsStatement(const vector<EstimationParams> &estim_params_list_arg,
                                 const SymbolTable &symbol_table_arg);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

class OptimWeightsStatement : public Statement
{
public:
  typedef map<string, NodeID> var_weights_type;
  typedef map<pair<string, string>, NodeID> covar_weights_type;
private:
  const var_weights_type var_weights;
  const covar_weights_type covar_weights;
  const SymbolTable &symbol_table;
public:
  OptimWeightsStatement(const var_weights_type &var_weights_arg,
                        const covar_weights_type &covar_weights_arg,
                        const SymbolTable &symbol_table_arg);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

class CalibStatement : public Statement
{
private:
  const int covar;
public:
  CalibStatement(int covar_arg);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

class CalibVarStatement : public Statement
{
public:
  //! Maps a variable to a pair (weight, expression)
  typedef map<string, pair<string, NodeID> > calib_var_type;
  //! Maps a pair of variables to a pair (weight, expression)
  typedef map<pair<string, string>, pair<string, NodeID> > calib_covar_type;
  //! Maps a pair (variable, autocorr) to a pair (weight, expression)
  typedef map<pair<string, int>, pair<string, NodeID> > calib_ac_type;
private:
  const calib_var_type calib_var;
  const calib_covar_type calib_covar;
  const calib_ac_type calib_ac;
  const SymbolTable &symbol_table;
public:
  CalibVarStatement(const calib_var_type &calib_var_arg,
                    const calib_covar_type &calib_covar_arg,
                    const calib_ac_type &calib_ac_arg,
                    const SymbolTable &symbol_table_arg);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

/*! \todo Make model_tree a member instead of a pointer */
class PlannerObjectiveStatement : public Statement
{
private:
  ModelTree *model_tree;
public:
  //! Constructor
  /*! \param model_tree_arg the model tree used to store the objective function.
    It is owned by the PlannerObjectiveStatement, and will be deleted by its destructor */
  PlannerObjectiveStatement(ModelTree *model_tree_arg);
  virtual ~PlannerObjectiveStatement();
  /*! \todo check there are only endogenous variables at the current period in the objective
    (no exogenous, no lead/lag) */
  virtual void checkPass(ModFileStructure &mod_file_struct);
  /*! \todo allow for the possibility of disabling temporary terms */
  virtual void computingPass();
  virtual void writeOutput(ostream &output, const string &basename) const;
};

class BVARDensityStatement : public Statement
{
private:
  const int maxnlags;
  const OptionsList options_list;
public:
  BVARDensityStatement(int maxnlags_arg, const OptionsList &options_list_arg);
  virtual void checkPass(ModFileStructure &mod_file_struct);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

class BVARForecastStatement : public Statement
{
private:
  const int nlags;
  const OptionsList options_list;
public:
  BVARForecastStatement(int nlags_arg, const OptionsList &options_list_arg);
  virtual void checkPass(ModFileStructure &mod_file_struct);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

#endif
