/*
 * Copyright © 2003-2023 Dynare Team
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

#ifndef _DYNAMICMODEL_HH
#define _DYNAMICMODEL_HH

#include <fstream>
#include <filesystem>

#include "StaticModel.hh"
#include "Bytecode.hh"

using namespace std;

//! Stores a dynamic model
class DynamicModel : public ModelTree
{
public:
  //! A reference to the trend component model table
  TrendComponentModelTable &trend_component_model_table;
  //! A reference to the VAR model table
  VarModelTable &var_model_table;
  /* Used in the balanced growth test, for determining whether the
     cross-derivative of a given equation, w.r.t. an endogenous and a trend
     variable is zero. Controlled by option “balanced_growth_test_tol” of the
     “model” block. The default should not be too small (see dynare#1389). */
  double balanced_growth_test_tol{1e-6};
private:
  /* Used in the balanced growth test, for skipping equations where the test
     cannot be performed (i.e. when LHS=RHS at the initial values). Should not
     be too large, otherwise the test becomes less powerful. */
  constexpr static double zero_band{1e-8};

  //! Stores equations declared as [static]
  /*! They will be used in the conversion to StaticModel to replace equations marked as [dynamic] */
  vector<BinaryOpNode *> static_only_equations;

  //! Stores line numbers of equations declared as [static]
  vector<optional<int>> static_only_equations_lineno;

  //! Stores the equation tags of equations declared as [static]
  EquationTags static_only_equations_equation_tags;

  using deriv_id_table_t = map<pair<int, int>, int>;
  //! Maps a pair (symbol_id, lag) to a deriv ID
  deriv_id_table_t deriv_id_table;
  //! Maps a deriv ID to a pair (symbol_id, lag)
  vector<pair<int, int>> inv_deriv_id_table;

  /* Maps a deriv_id to the column index of the dynamic Jacobian, in the legacy
     representation.
     Contains only endogenous, exogenous and exogenous deterministic */
  map<int, int> dyn_jacobian_cols_table;
  // Number of columns of the dynamic Jacobian (legacy representation)
  int dyn_jacobian_ncols;

  //! Maximum lag and lead over all types of variables (positive values)
  /*! Set by computeDerivIDs() */
  int max_lag{0}, max_lead{0};
  //! Maximum lag and lead over endogenous variables (positive values)
  /*! Set by computeDerivIDs() */
  int max_endo_lag{0}, max_endo_lead{0};
  //! Maximum lag and lead over exogenous variables (positive values)
  /*! Set by computeDerivIDs() */
  int max_exo_lag{0}, max_exo_lead{0};
  //! Maximum lag and lead over deterministic exogenous variables (positive values)
  /*! Set by computeDerivIDs() */
  int max_exo_det_lag{0}, max_exo_det_lead{0};
  //! Maximum lag and lead over all types of variables (positive values) of original model
  int max_lag_orig{0}, max_lead_orig{0}, max_lag_with_diffs_expanded_orig{0};
  //! Maximum lag and lead over endogenous variables (positive values) of original model
  int max_endo_lag_orig{0}, max_endo_lead_orig{0};
  //! Maximum lag and lead over exogenous variables (positive values) of original model
  int max_exo_lag_orig{0}, max_exo_lead_orig{0};
  //! Maximum lag and lead over deterministic exogenous variables (positive values) of original model
  int max_exo_det_lag_orig{0}, max_exo_det_lead_orig{0};

  // Cross reference information: eq → set of (symb_id, lag) for each symbol type
  map<int, ExprNode::EquationInfo> xrefs;
  // Reverse cross reference information: (symb_id, lag) → eqs
  map<pair<int, int>, set<int>> xref_param, xref_endo, xref_exo, xref_exo_det;

  //! Nonzero equations in the Hessian
  set<int> nonzero_hessian_eqs;

  //! Creates mapping for variables and equations they are present in
  map<int, set<int>> variableMapping;

  /* For each block, and for each variable type, maps (variable ID, lag) to
     Jacobian column. The variable ID is the index within the block. */
  vector<map<pair<int, int>, int>> blocks_jacob_cols_endo;

  //! Used for var_expectation and var_model
  map<string, set<int>> var_expectation_functions_to_write;

  // Writes dynamic model file (MATLAB/Octave version, legacy representation)
  void writeDynamicMFile(const string &basename) const;
  //! Writes the code of the block-decomposed model in virtual machine bytecode
  void writeDynamicBlockBytecode(const string &basename) const;
  //! Writes the code of the model in virtual machine bytecode
  void writeDynamicBytecode(const string &basename) const;

  // Write the block structure of the model in the driver file
  void writeBlockDriverOutput(ostream &output) const;

  // Used by determineBlockDerivativesType()
  enum class BlockDerivativeType
    {
     standard,
     chainRule,
     normalizedChainRule
    };

  /* For each tuple (lag, eq, var) within the given block, determine the type
     of the derivative to be computed. Indices are within the block (i.e.
     between 0 and blocks[blk].size-1). */
  map<tuple<int, int, int>, BlockDerivativeType> determineBlockDerivativesType(int blk);

  void computeChainRuleJacobian() override;

  string reform(const string &name) const;

  SymbolType getTypeByDerivID(int deriv_id) const noexcept(false) override;
  int getLagByDerivID(int deriv_id) const noexcept(false) override;
  int getSymbIDByDerivID(int deriv_id) const noexcept(false) override;
  int getTypeSpecificIDByDerivID(int deriv_id) const override;

  //! Compute the column indices of the dynamic Jacobian
  void computeDynJacobianCols();
  //! Computes derivatives of the Jacobian w.r. to trend vars and tests that they are equal to zero
  void testTrendDerivativesEqualToZero(const eval_context_t &eval_context);

  //! Allocates the derivation IDs for all dynamic variables of the model
  /*! Also computes max_{endo,exo}_{lead_lag}, and initializes dynJacobianColsNbr to the number of dynamic endos */
  void computeDerivIDs();

  /* Compute the Jacobian column indices in the block decomposition case
     (stored in blocks_jacob_cols_*).
     Also fills auxiliary structures related to “other” endogenous and
     exogenous: blocks{,_derivatives}_{other_endo,exo_exo_det} */
  void computeBlockDynJacobianCols();

  //! Factorized code for substitutions of leads/lags
  /*! \param[in] type determines which type of variables is concerned
    \param[in] deterministic_model whether we are in a deterministic model (only for exogenous leads/lags)
    \param[in] subset variables to which to apply the transformation (only for diff of forward vars)
  */
  void substituteLeadLagInternal(AuxVarType type, bool deterministic_model, const vector<string> &subset);

  //! Help computeXrefs to compute the reverse references (i.e. param->eqs, endo->eqs, etc)
  void computeRevXref(map<pair<int, int>, set<int>> &xrefset, const set<pair<int, int>> &eiref, int eqn);

  //! Write reverse cross references
  void writeRevXrefs(ostream &output, const map<pair<int, int>, set<int>> &xrefmap, const string &type) const;

  /* Writes MATLAB/Octave wrapper function for computing residuals and
     derivatives at the same time (legacy representation) */
  void writeDynamicMWrapperFunction(const string &name, const string &ending) const;
  /* Helper for writing MATLAB/Octave functions for residuals/derivatives and
     their temporary terms (legacy representation) */
  void writeDynamicMFileHelper(const string &basename,
                               const string &name, const string &retvalname,
                               const string &name_tt, size_t ttlen,
                               const string &previous_tt_name,
                               const ostringstream &init_s, const ostringstream &end_s,
                               const ostringstream &s, const ostringstream &s_tt) const;

  /* Create the compatibility dynamic.m file for MATLAB/Octave not yet using
     the temporary terms array interface (legacy representation) */
  void writeDynamicMCompatFile(const string &basename) const;

  //! Internal helper for the copy constructor and assignment operator
  /*! Copies all the structures that contain ExprNode*, by the converting the
      pointers into their equivalent in the new tree */
  void copyHelper(const DynamicModel &m);

  /* Handles parsing of argument passed to exclude_eqs/include_eqs.

    The argument inc_exc_option_value should be of one of the following forms:
      * filename.txt
      * eq1
      * ['eq 1', 'eq 2']
      * [tagname='eq 1']
      * [tagname=('eq 1', 'eq 2')]
    If argument is a filename, the file should be formatted as:
        eq 1
        eq 2
    OR
        tagname=
        X
        Y

    The boolean exclude_eqs should be true if we are in the exclude_eqs case,
    false in the include_eqs case (this only affects error messages).

    Returns a set of pairs (tag name, tag value) corresponding to the set of
    equations to be included or excluded.
   */
  static vector<pair<string, string>> parseIncludeExcludeEquations(const string &inc_exc_option_value, bool exclude_eqs);

  /* Helper for the removeEquations() method.
     listed_eqs_by_tag is the list of (tag name, tag value) pairs corresponding
     to the option value, exclude_eqs is a boolean indicating whether we’re
     excluding or including, and excluded_vars_change_type is a boolean
     indicating whether to compute variables to be excluded.

     The all_equations* arguments will be modified by the routine by excluding
     equations. They are either the main structures for storing equations in
     ModelTree, or their counterpart for static-only equations. The
     static_equations boolean indicates when we are in the latter case.
     The listed_eqs_by_tag structure will be updated by removing those tag
     pairs that have been matched with equations in the all_equations*
     argument*.

     Returns a list of excluded variables (empty if
     excluded_vars_change_type=false) */
  vector<int> removeEquationsHelper(set<pair<string, string>> &listed_eqs_by_tag,
                                    bool exclude_eqs, bool excluded_vars_change_type,
                                    vector<BinaryOpNode *> &all_equations,
                                    vector<optional<int>> &all_equations_lineno,
                                    EquationTags &all_equation_tags,
                                    bool static_equations) const;

  //! Compute autoregressive matrices of trend component models
  /* The algorithm uses matching rules over expression trees. It cannot handle
     arbitrarily-written expressions. */
  map<string, map<tuple<int, int, int>, expr_t>> computeAutoregressiveMatrices() const;

  //! Compute error component matrices of trend component_models
  /*! Returns a pair (A0r, A0starr) */
  pair<map<string, map<tuple<int, int>, expr_t>>, map<string, map<tuple<int, int>, expr_t>>> computeErrorComponentMatrices(const ExprNode::subst_table_t &diff_subst_table) const;

  /* For a VAR model, given the symbol ID of a LHS variable, and a (negative)
     lag, returns all the corresponding deriv_ids (by properly dealing with two
     types of auxiliary variables: endo lags and diff lags). It returns a
     vector because in some cases there may be sereval corresponding deriv_ids
     (for example, in the deriv_id table, AUX_DIFF_nn(-1) may appear as itself
     (with a lag), and also as a contemporaneous diff lag auxvar). */
  vector<int> getVARDerivIDs(int lhs_symb_id, int lead_lag) const;

  int
  getBlockJacobianEndoCol(int blk, int var, int lag) const override
  {
    return blocks_jacob_cols_endo[blk].at({ var, lag });
  }

protected:
  string
  modelClassName() const override
  {
    return "dynamic model";
  }

public:
  DynamicModel(SymbolTable &symbol_table_arg,
               NumericalConstants &num_constants_arg,
               ExternalFunctionsTable &external_functions_table_arg,
               TrendComponentModelTable &trend_component_model_table_arg,
               VarModelTable &var_model_table_arg);

  DynamicModel(const DynamicModel &m);
  DynamicModel &operator=(const DynamicModel &m);

  //! Compute cross references
  void computeXrefs();

  //! Write cross references
  void writeXrefs(ostream &output) const;

  //! Execute computations (variable sorting + derivation + block decomposition)
  /*!
    \param derivsOrder order of derivatives w.r. to exo, exo_det and endo should be computed (implies jacobianExo = true when order >= 2)
    \param paramsDerivsOrder order of derivatives w.r. to a pair (endo/exo/exo_det, parameter) to be computed (>0 implies jacobianExo = true)
    \param eval_context evaluation context for normalization
    \param no_tmp_terms if true, no temporary terms will be computed in the dynamic files
  */
  void computingPass(int derivsOrder, int paramsDerivsOrder, const eval_context_t &eval_context,
                     bool no_tmp_terms, bool block, bool use_dll);
  //! Writes information about the dynamic model to the driver file
  void writeDriverOutput(ostream &output, bool compute_xrefs) const;

  //! Write JSON AST
  void writeJsonAST(ostream &output) const;

  //! Write JSON variable mapping
  void writeJsonVariableMapping(ostream &output) const;

  //! Write JSON Output
  void writeJsonOutput(ostream &output) const;

  //! Write JSON Output representation of original dynamic model
  void writeJsonOriginalModelOutput(ostream &output) const;

  //! Write JSON Output representation of model info (useful stuff from M_)
  void writeJsonDynamicModelInfo(ostream &output) const;

  //! Write JSON Output representation of dynamic model after computing pass
  void writeJsonComputingPassOutput(ostream &output, bool writeDetails) const;

  //! Write JSON params derivatives
  void writeJsonParamsDerivatives(ostream &output, bool writeDetails) const;

  //! Write cross reference output if the xref maps have been filed
  void writeJsonXrefs(ostream &output) const;
  void writeJsonXrefsHelper(ostream &output, const map<pair<int, int>, set<int>> &xrefmap) const;

  //! Print equations that have non-zero second derivatives
  void printNonZeroHessianEquations(ostream &output) const;

  //! Tells whether Hessian has been computed
  /*! This is needed to know whether no non-zero equation in Hessian means a
    zero Hessian or Hessian not computed */
  bool
  isHessianComputed() const
  {
    return computed_derivs_order >= 2;
  }

  /* Check whether the model is linear, by verifying that the hessian is zero,
     and error out otherwise.
     Must be called after computingPass().
     FIXME: this check always passes if derivsOrder = 1, i.e. for a perfect
     foresight model, because the Hessian is not computed in that case. */
  void checkIsLinear() const;

  //! Fill the trend component model table with information available from the transformed model
  void fillTrendComponentModelTable() const;
  //! Fill the trend component model table with information available from the original model
  void fillTrendComponentModelTableFromOrigModel() const;
  /* Fill the trend component model table with information about AR/EC
     components, available from the transformed model. Needs to be called after
     fillTrendComponentModelTableFromOrigModel() has been called on the
     original model */
  void fillTrendComponentModelTableAREC(const ExprNode::subst_table_t &diff_subst_table) const;

  //! Fill the VAR model table with information available from the transformed model
  // NB: Does not fill the AR and A0 matrices
  void fillVarModelTable() const;
  //! Fill the VAR model table with information available from the original model
  void fillVarModelTableFromOrigModel() const;
  //! Fill the AR and A0 matrices of the VAR model table
  // Uses derivatives, hence must be called after computingPass()
  void fillVarModelTableMatrices();

  //! Update the rhs references in the var model and trend component tables
  //! after substitution of auxiliary variables and find the trend variables
  //! in the trend_component model
  void updateVarAndTrendModel() const;

  //! Writes dynamic model file (+ bytecode)
  void writeDynamicFile(const string &basename, bool use_dll, const string &mexext, const filesystem::path &matlabroot, bool julia) const;

  //! Writes file containing parameters derivatives
  template<bool julia>
  void writeParamsDerivativesFile(const string &basename) const;

  //! Creates mapping for variables and equations they are present in
  void createVariableMapping();

  //! Expands equation tags with default equation names (available "name" tag or LHS variable or equation ID)
  void expandEqTags();

  //! Find endogenous variables not used in model
  set<int> findUnusedEndogenous();
  //! Find exogenous variables not used in model
  set<int> findUnusedExogenous();

  //! Set the max leads/lags of the original model
  void setLeadsLagsOrig();

  //! Implements the include_eqs/exclude_eqs options
  void includeExcludeEquations(const string &inc_exc_option_value, bool exclude_eqs);

  /* Removes equations from the model (identified by their name tags).
     Used for include_eqs/exclude_eqs options and for model_remove and
     model_replace blocks */
  void removeEquations(const vector<pair<string, string>> &listed_eqs_by_tag, bool exclude_eqs,
                       bool excluded_vars_change_type);

  /* Replaces model equations with derivatives of Lagrangian w.r.t. endogenous.
     The optimality FOCs (derivatives w.r.t. ordinary endogenous) will appear
     first, followed by the constraints (derivatives w.r.t. multipliers).
     Returns the number of optimality FOCs, which is by construction equal to
     the number of endogenous before adding the Lagrange multipliers
     (internally called ramsey_endo_nbr). */
  int computeRamseyPolicyFOCs(const StaticModel &static_model);

  //! Clears all equations
  void clearEquations();

  //! Replaces the model equations in dynamic_model with those in this model
  void replaceMyEquations(DynamicModel &dynamic_model) const;

  //! Adds an equation marked as [static]
  void addStaticOnlyEquation(expr_t eq, optional<int> lineno, const map<string, string> &eq_tags);

  //! Returns number of static only equations
  size_t staticOnlyEquationsNbr() const;

  //! Returns number of dynamic only equations
  size_t dynamicOnlyEquationsNbr() const;

  // Adds an occbin equation (with “bind” and/or “relax” tag)
  /* This function assumes that there is a “name” tag, and that the relevant
     auxiliary parameters have already been added to the symbol table.
     It also assumes that the “bind” and “relax” tags have been cleared from
     eq_tags. */
  void addOccbinEquation(expr_t eq, optional<int> lineno, const map<string, string> &eq_tags, const vector<string> &regimes_bind, const vector<string> &regimes_relax);

  //! Writes LaTeX file with the equations of the dynamic model
  void writeLatexFile(const string &basename, bool write_equation_tags) const;

  //! Writes LaTeX file with the equations of the dynamic model (for the original model)
  void writeLatexOriginalFile(const string &basename, bool write_equation_tags) const;

  int getDerivID(int symb_id, int lag) const noexcept(false) override;

  int
  getJacobianCol(int deriv_id, bool sparse) const override
  {
    if (sparse)
      {
        SymbolType type {getTypeByDerivID(deriv_id)};
        int tsid {getTypeSpecificIDByDerivID(deriv_id)};
        int lag {getLagByDerivID(deriv_id)};
        if (type == SymbolType::endogenous)
          {
            assert(lag >= -1 && lag <= 1);
            return tsid+(lag+1)*symbol_table.endo_nbr();
          }
        else if (type == SymbolType::exogenous)
          {
            assert(lag == 0);
            return tsid+3*symbol_table.endo_nbr();
          }
        else if (type == SymbolType::exogenousDet)
          {
            assert(lag == 0);
            return tsid+3*symbol_table.endo_nbr()+symbol_table.exo_nbr();
          }
        else
          throw UnknownDerivIDException();
      }
    else
      {
        if (auto it = dyn_jacobian_cols_table.find(deriv_id);
            it == dyn_jacobian_cols_table.end())
          throw UnknownDerivIDException();
        else
          return it->second;
      }
  }
  int
  getJacobianColsNbr(bool sparse) const override
  {
    return sparse ?
      3*symbol_table.endo_nbr() + symbol_table.exo_nbr() + symbol_table.exo_det_nbr() :
      dyn_jacobian_ncols;
  }

  void addAllParamDerivId(set<int> &deriv_id_set) override;

  //! Returns true indicating that this is a dynamic model
  bool
  isDynamic() const override
  {
    return true;
  };

  //! Drive test of detrended equations
  void runTrendTest(const eval_context_t &eval_context);

  //! Transforms the model by removing all leads greater or equal than 2 on endos
  /*! Note that this can create new lags on endos and exos */
  void substituteEndoLeadGreaterThanTwo(bool deterministic_model);

  //! Transforms the model by removing all lags greater or equal than 2 on endos
  void substituteEndoLagGreaterThanTwo(bool deterministic_model);

  //! Transforms the model by removing all leads on exos
  /*! Note that this can create new lags on endos and exos */
  void substituteExoLead(bool deterministic_model);

  //! Transforms the model by removing all lags on exos
  void substituteExoLag(bool deterministic_model);

  //! Transforms the model by removing all UnaryOpcode::expectation
  void substituteExpectation(bool partial_information_model);

  //! Transforms the model by decreasing the lead/lag of predetermined variables in model equations by one
  void transformPredeterminedVariables();

  // Performs the transformations associated to variables declared with “var(log)”
  void substituteLogTransform();

  // Check that no variable was declared with “var(log)” in the given equations
  void checkNoWithLogTransform(const set<int> &eqnumbers);

  //! Transforms the model by removing trends specified by the user
  void detrendEquations();

  const nonstationary_symbols_map_t &
  getNonstationarySymbolsMap() const
  {
    return nonstationary_symbols_map;
  }

  const map<int, expr_t> &
  getTrendSymbolsMap() const
  {
    return trend_symbols_map;
  }

  //! Substitutes adl operator
  void substituteAdl();

  //! Substitutes out all model-local variables
  void substituteModelLocalVariables();

  /* Creates aux vars for all unary operators in all equations. Also makes the
     substitution in growth terms of pac_model/pac_target_info and in
     expressions of var_expectation_model. */
  pair<lag_equivalence_table_t, ExprNode::subst_table_t> substituteUnaryOps(VarExpectationModelTable &var_expectation_model_table, PacModelTable &pac_model_table);

  /* Creates aux vars for all unary operators in specified equations. Also makes the
     substitution in growth terms of pac_model/pac_target_info and in
     expressions of var_expectation_model. */
  pair<lag_equivalence_table_t, ExprNode::subst_table_t> substituteUnaryOps(const set<int> &eqnumbers, VarExpectationModelTable &var_expectation_model_table, PacModelTable &pac_model_table);

  //! Substitutes diff operator
  pair<lag_equivalence_table_t, ExprNode::subst_table_t> substituteDiff(VarExpectationModelTable &var_expectation_model_table, PacModelTable &pac_model_table);

  //! Substitute VarExpectation operators
  void substituteVarExpectation(const map<string, expr_t> &subst_table);

  void analyzePacEquationStructure(const string &name, map<string, string> &pac_eq_name, PacModelTable::equation_info_t &pac_equation_info);

  // Exception thrown by getPacTargetSymbId()
  struct PacTargetNotIdentifiedException
  {
    const string model_name, message;
  };

  //! Return target of the pac equation
  int getPacTargetSymbId(const string &pac_model_name) const;

  /* For a PAC MCE model, fill pac_expectation_substitution with the
     expression that will be substituted for the pac_expectation operator.
     In the process, add the variable and the equation defining Z₁.
     The symbol IDs of the new endogenous are added to pac_aux_var_symb_ids,
     and the new auxiliary parameters to pac_mce_alpha_symb_ids.
  */
  void computePacModelConsistentExpectationSubstitution(const string &name,
                                                        int discount_symb_id, int pac_eq_max_lag,
                                                        expr_t growth_correction_term,
                                                        string auxname,
                                                        ExprNode::subst_table_t &diff_subst_table,
                                                        map<string, int> &pac_aux_var_symb_ids,
                                                        map<string, vector<int>> &pac_aux_param_symb_ids,
                                                        map<string, expr_t> &pac_expectation_substitution);


  /* For a PAC backward model, fill pac_expectation_substitution with the
     expression that will be substituted for the pac_expectation operator.
     The symbol IDs of the new parameters are also added to pac_aux_param_symb_ids.
     The symbol ID of the new auxiliary variable is added to pac_aux_var_symb_ids. */
  void computePacBackwardExpectationSubstitution(const string &name,
                                                 const vector<int> &lhs,
                                                 int max_lag,
                                                 const string &aux_model_type,
                                                 expr_t growth_correction_term,
                                                 string auxname,
                                                 map<string, int> &pac_aux_var_symb_ids,
                                                 map<string, vector<int>> &pac_aux_param_symb_ids,
                                                 map<string, expr_t> &pac_expectation_substitution);

  /* Same as above, but for PAC models which have an associated
     pac_target_info.
     Contrary to the above routine, this one will create the growth correction
     parameters as needed.
     Those parameter IDs, as well as the IDs for the h parameters, are stored
     in target_components.
     The routine also creates the auxiliary variables for the components, and
     adds the corresponding equations. */
  void computePacBackwardExpectationSubstitutionWithComponents(const string &name,
                                                               const vector<int> &lhs,
                                                               int max_lag,
                                                               const string &aux_model_type,
                                                               vector<PacModelTable::target_component_t> &pac_target_components,
                                                               map<string, expr_t> &pac_expectation_substitution);

  //! Substitutes pac_expectation operator with expectation based on auxiliary model
  void substitutePacExpectation(const map<string, expr_t> &pac_expectation_substitution,
                                const map<string, string> &pac_eq_name);

  //! Substitutes the pac_target_nonstationary operator of a given pac_model
  void substitutePacTargetNonstationary(const string &pac_model_name, expr_t substexpr);

  //! Table to undiff LHS variables for pac vector z
  vector<int> getUndiffLHSForPac(const string &aux_model_name,
                                 const ExprNode::subst_table_t &diff_subst_table) const;

  //! Transforms the model by replacing trend variables with a 1
  void removeTrendVariableFromEquations();

  //! Transforms the model by creating aux vars for the diff of forward vars
  /*! If subset is empty, does the transformation for all fwrd vars; otherwise
    restrict it to the vars in subset */
  void differentiateForwardVars(const vector<string> &subset);

  //! Fills eval context with values of model local variables and auxiliary variables
  void fillEvalContext(eval_context_t &eval_context) const;

  /*! Checks that all pac_expectation operators have been substituted, error
    out otherwise */
  void checkNoRemainingPacExpectation() const;

  /*! Checks that all pac_target_nonstationary operators have been substituted, error
    out otherwise */
  void checkNoRemainingPacTargetNonstationary() const;

  auto
  getStaticOnlyEquationsInfo() const
  {
    return tuple{static_only_equations, static_only_equations_lineno, static_only_equations_equation_tags};
  };

  //! Returns true if a parameter was used in the model block with a lead or lag
  bool ParamUsedWithLeadLag() const;

  bool isChecksumMatching(const string &basename) const;

  //! Simplify model equations: if a variable is equal to a constant, replace that variable elsewhere in the model
  /*! Equations with MCP tags are excluded, see dynare#1697 */
  void simplifyEquations();

  // Converts a set of equation tags into the corresponding set of equation numbers
  set<int> getEquationNumbersFromTags(const set<string> &eqtags) const;

  // Returns the set of equations (as numbers) which have a pac_expectation operator
  set<int> findPacExpectationEquationNumbers() const;
};

template<bool julia>
void
DynamicModel::writeParamsDerivativesFile(const string &basename) const
{
  if (!params_derivatives.size())
    return;

  constexpr ExprNodeOutputType output_type { julia ? ExprNodeOutputType::juliaDynamicModel : ExprNodeOutputType::matlabDynamicModel };

  auto [tt_output, rp_output, gp_output, rpp_output, gpp_output, hp_output, g3p_output]
    { writeParamsDerivativesFileHelper<output_type>() };

  if constexpr(!julia)
    {
      filesystem::path filename {packageDir(basename) / "dynamic_params_derivs.m"};
      ofstream paramsDerivsFile {filename, ios::out | ios::binary};
      if (!paramsDerivsFile.is_open())
        {
          cerr << "ERROR: Can't open file " << filename.string() << " for writing" << endl;
          exit(EXIT_FAILURE);
        }
      paramsDerivsFile << "function [rp, gp, rpp, gpp, hp, g3p] = dynamic_params_derivs(y, x, params, steady_state, it_, ss_param_deriv, ss_param_2nd_deriv)" << endl
                       << "%" << endl
                       << "% Compute the derivatives of the dynamic model with respect to the parameters" << endl
                       << "% Inputs :" << endl
                       << "%   y         [#dynamic variables by 1] double    vector of endogenous variables in the order stored" << endl
                       << "%                                                 in M_.lead_lag_incidence; see the Manual" << endl
                       << "%   x         [nperiods by M_.exo_nbr] double     matrix of exogenous variables (in declaration order)" << endl
                       << "%                                                 for all simulation periods" << endl
                       << "%   params    [M_.param_nbr by 1] double          vector of parameter values in declaration order" << endl
                       << "%   steady_state  [M_.endo_nbr by 1] double       vector of steady state values" << endl
                       << "%   it_       scalar double                       time period for exogenous variables for which to evaluate the model" << endl
                       << "%   ss_param_deriv     [M_.eq_nbr by #params]     Jacobian matrix of the steady states values with respect to the parameters" << endl
                       << "%   ss_param_2nd_deriv [M_.eq_nbr by #params by #params] Hessian matrix of the steady states values with respect to the parameters" << endl
                       << "%" << endl
                       << "% Outputs:" << endl
                       << "%   rp        [M_.eq_nbr by #params] double    Jacobian matrix of dynamic model equations with respect to parameters " << endl
                       << "%                                              Dynare may prepend or append auxiliary equations, see M_.aux_vars" << endl
                       << "%   gp        [M_.endo_nbr by #dynamic variables by #params] double    Derivative of the Jacobian matrix of the dynamic model equations with respect to the parameters" << endl
                       << "%                                                           rows: equations in order of declaration" << endl
                       << "%                                                           columns: variables in order stored in M_.lead_lag_incidence" << endl
                       << "%   rpp       [#second_order_residual_terms by 4] double   Hessian matrix of second derivatives of residuals with respect to parameters;" << endl
                       << "%                                                              rows: respective derivative term" << endl
                       << "%                                                              1st column: equation number of the term appearing" << endl
                       << "%                                                              2nd column: number of the first parameter in derivative" << endl
                       << "%                                                              3rd column: number of the second parameter in derivative" << endl
                       << "%                                                              4th column: value of the Hessian term" << endl
                       << "%   gpp      [#second_order_Jacobian_terms by 5] double   Hessian matrix of second derivatives of the Jacobian with respect to the parameters;" << endl
                       << "%                                                              rows: respective derivative term" << endl
                       << "%                                                              1st column: equation number of the term appearing" << endl
                       << "%                                                              2nd column: column number of variable in Jacobian of the dynamic model" << endl
                       << "%                                                              3rd column: number of the first parameter in derivative" << endl
                       << "%                                                              4th column: number of the second parameter in derivative" << endl
                       << "%                                                              5th column: value of the Hessian term" << endl
                       << "%   hp      [#first_order_Hessian_terms by 5] double   Jacobian matrix of derivatives of the dynamic Hessian with respect to the parameters;" << endl
                       << "%                                                              rows: respective derivative term" << endl
                       << "%                                                              1st column: equation number of the term appearing" << endl
                       << "%                                                              2nd column: column number of first variable in Hessian of the dynamic model" << endl
                       << "%                                                              3rd column: column number of second variable in Hessian of the dynamic model" << endl
                       << "%                                                              4th column: number of the parameter in derivative" << endl
                       << "%                                                              5th column: value of the Hessian term" << endl
                       << "%   g3p      [#first_order_g3_terms by 6] double   Jacobian matrix of derivatives of g3 (dynamic 3rd derivs) with respect to the parameters;" << endl
                       << "%                                                              rows: respective derivative term" << endl
                       << "%                                                              1st column: equation number of the term appearing" << endl
                       << "%                                                              2nd column: column number of first variable in g3 of the dynamic model" << endl
                       << "%                                                              3rd column: column number of second variable in g3 of the dynamic model" << endl
                       << "%                                                              4th column: column number of third variable in g3 of the dynamic model" << endl
                       << "%                                                              5th column: number of the parameter in derivative" << endl
                       << "%                                                              6th column: value of the Hessian term" << endl
                       << "%" << endl
                       << "%" << endl
                       << "% Warning : this file is generated automatically by Dynare" << endl
                       << "%           from model file (.mod)" << endl << endl
                       << "T = NaN(" << params_derivs_temporary_terms_idxs.size() << ",1);" << endl
                       << tt_output.str()
                       << "rp = zeros(" << equations.size() << ", "
                       << symbol_table.param_nbr() << ");" << endl
                       << rp_output.str()
                       << "gp = zeros(" << equations.size() << ", " << getJacobianColsNbr(false) << ", " << symbol_table.param_nbr() << ");" << endl
                       << gp_output.str()
                       << "if nargout >= 3" << endl
                       << "rpp = zeros(" << params_derivatives.at({ 0, 2 }).size() << ",4);" << endl
                       << rpp_output.str()
                       << "gpp = zeros(" << params_derivatives.at({ 1, 2 }).size() << ",5);" << endl
                       << gpp_output.str()
                       << "end" << endl
                       << "if nargout >= 5" << endl
                       << "hp = zeros(" << params_derivatives.at({ 2, 1 }).size() << ",5);" << endl
                       << hp_output.str()
                       << "end" << endl
                       << "if nargout >= 6" << endl
                       << "g3p = zeros(" << params_derivatives.at({ 3, 1 }).size() << ",6);" << endl
                       << g3p_output.str()
                       << "end" << endl
                       << "end" << endl;
      paramsDerivsFile.close();
    }
  else
    {
      stringstream output;
      output << "# NB: this file was automatically generated by Dynare" << endl
             << "#     from " << basename << ".mod" << endl
             << "#" << endl
             << "function dynamic_params_derivs(y, x, params, steady_state, it_,"
             << "ss_param_deriv, ss_param_2nd_deriv)" << endl
             << "@inbounds begin" << endl
             << tt_output.str()
             << "rp = zeros(" << equations.size() << ", "
             << symbol_table.param_nbr() << ");" << endl
             << rp_output.str()
             << "gp = zeros(" << equations.size() << ", " << getJacobianColsNbr(false) << ", " << symbol_table.param_nbr() << ");" << endl
             << gp_output.str()
             << "rpp = zeros(" << params_derivatives.at({ 0, 2 }).size() << ",4);" << endl
             << rpp_output.str()
             << "gpp = zeros(" << params_derivatives.at({ 1, 2 }).size() << ",5);" << endl
             << gpp_output.str()
             << "hp = zeros(" << params_derivatives.at({ 2, 1 }).size() << ",5);" << endl
             << hp_output.str()
             << "g3p = zeros(" << params_derivatives.at({ 3, 1 }).size() << ",6);" << endl
             << g3p_output.str()
             << "end" << endl
             << "return (rp, gp, rpp, gpp, hp, g3p)" << endl
             << "end" << endl;

      writeToFileIfModified(output, filesystem::path{basename} / "model" / "julia" / "DynamicParamsDerivs.jl");
    }
}

#endif
