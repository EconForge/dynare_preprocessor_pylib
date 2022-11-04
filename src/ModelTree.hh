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

#ifndef _MODELTREE_HH
#define _MODELTREE_HH

#include <string>
#include <vector>
#include <deque>
#include <map>
#include <ostream>
#include <array>
#include <filesystem>
#include <optional>
#include <cassert>
#include <thread>
#include <mutex>
#include <condition_variable>

#include "DataTree.hh"
#include "EquationTags.hh"
#include "ExtendedPreprocessorTypes.hh"
#include "Bytecode.hh"

using namespace std;

// Helper to convert a vector into a tuple
template<typename T, size_t... Indices>
auto
vectorToTupleHelper(const vector<T> &v, index_sequence<Indices...>)
{
  return tuple{v[Indices]...};
}
template<size_t N, typename T>
auto
vectorToTuple(const vector<T> &v)
{
  assert(v.size() >= N);
  return vectorToTupleHelper(v, make_index_sequence<N>());
}

//! Vector describing equations: BlockSimulationType, if BlockSimulationType == EVALUATE_s then a expr_t on the new normalized equation
using equation_type_and_normalized_equation_t = vector<pair<EquationType, BinaryOpNode *>>;

//! Vector describing variables: max_lag in the block, max_lead in the block
using lag_lead_vector_t = vector<pair<int, int>>;

//! Shared code for static and dynamic models
class ModelTree : public DataTree
{
  friend class DynamicModel;
  friend class StaticModel;
public:
  // Set via the `compiler` command
  string user_set_add_flags, user_set_subst_flags, user_set_add_libs, user_set_subst_libs, user_set_compiler;
protected:
  /*
   * ************** BEGIN **************
   * The following structures keep track of the model equations and must all be updated
   * when adding or removing an equation. Hence, if a new parallel structure is added
   * in the future, it must be maintained whereever these structures are updated
   * See in particular methods clearEquations(), replaceMyEquations() and
   * computeRamseyPolicyFOCs() of DynamicModel class.
   * NB: This message added with the introduction of the `exclude_eqs` option, hence
   *     that's a place to update future structures.
   */
  //! Stores declared and generated auxiliary equations
  vector<BinaryOpNode *> equations;
  /* Stores line numbers of declared equations; undefined in some cases (e.g.
     auxiliary equations) */
  vector<optional<int>> equations_lineno;
  //! Stores equation tags
  EquationTags equation_tags;
  /*
   * ************** END **************
   */

  //! Only stores generated auxiliary equations, in an order meaningful for evaluation
  /*! These equations only contain the definition of auxiliary variables, and
      may diverge from those in the main model (equations), if other model
      transformations applied subsequently. This is not a problem, since
      aux_equations is only used for regenerating the values of auxiliaries
      given the others.

      For example, such a divergence appears when there is an expectation
      operator in a ramsey model, see
      tests/optimal_policy/nk_ramsey_expectation.mod */
  vector<BinaryOpNode *> aux_equations;

  //! Maximum order at which (endogenous) derivatives have been computed
  int computed_derivs_order{0};

  //! Stores derivatives
  /*! Index 0 is not used, index 1 contains first derivatives, ...
     For each derivation order, stores a map whose key is a vector of integer: the
     first integer is the equation index, the remaining ones are the derivation
     IDs of variables (in non-decreasing order, to avoid storing symmetric
     elements several times). Only non-zero derivatives are stored. */
  vector<map<vector<int>, expr_t>> derivatives;

  //! Number of non-zero derivatives
  /*! Index 0 is not used, index 1 contains number of non-zero first
    derivatives, ... */
  vector<int> NNZDerivatives;

  // Used to order pairs of indices (row, col) according to column-major order
  struct columnMajorOrderLess
  {
    bool
    operator()(const pair<int, int> &p1, const pair<int, int> &p2) const
    {
      return p1.second < p2.second || (p1.second == p2.second && p1.first < p2.first);
    }
  };
  using SparseColumnMajorOrderMatrix = map<pair<int, int>, expr_t, columnMajorOrderLess>;
  /* The nonzero values of the sparse Jacobian in column-major order (which is
     the natural order for Compressed Sparse Column (CSC) storage).
     The pair of indices is (row, column). */
  SparseColumnMajorOrderMatrix jacobian_sparse_column_major_order;
  /* Column indices for the sparse Jacobian in Compressed Sparse Column (CSC)
     storage (corresponds to the “jc” vector in MATLAB terminology) */
  vector<int> jacobian_sparse_colptr;

  //! Derivatives with respect to parameters
  /*! The key of the outer map is a pair (derivation order w.r.t. endogenous,
  derivation order w.r.t. parameters). For e.g., { 1, 2 } corresponds to the jacobian
  differentiated twice w.r.t. to parameters.
  In inner maps, the vector of integers consists of: the equation index, then
  the derivation IDs of endogenous (in non-decreasing order),
  then the IDs of parameters (in non-decreasing order)*/
  map<pair<int,int>, map<vector<int>, expr_t>> params_derivatives;

  //! Temporary terms for residuals and derivatives
  /*! Index 0 is temp. terms of residuals, index 1 for first derivatives, ... */
  vector<temporary_terms_t> temporary_terms_derivatives;

  //! Stores, for each temporary term, its index in the MATLAB/Julia vector
  temporary_terms_idxs_t temporary_terms_idxs;

  //! Temporary terms for parameter derivatives, under a disaggregated form
  /*! The pair of integers is to be interpreted as in param_derivatives */
  map<pair<int, int>, temporary_terms_t> params_derivs_temporary_terms;

  //! Stores, for each temporary term in param. derivs, its index in the MATLAB/Julia vector
  temporary_terms_idxs_t params_derivs_temporary_terms_idxs;

  //! Trend variables and their growth factors
  map<int, expr_t> trend_symbols_map;

  //! for all trends; the boolean is true if this is a log-trend, false otherwise
  using nonstationary_symbols_map_t = map<int, pair<bool, expr_t>>;

  //! Nonstationary variables and their deflators
  nonstationary_symbols_map_t nonstationary_symbols_map;

  /* Maps indices of equations in the block-decomposition order into original
     equation IDs */
  vector<int> eq_idx_block2orig;
  /* Maps indices of (endogenous) variables in the block-decomposition order into original
     type-specific endogenous IDs */
  vector<int> endo_idx_block2orig;
  /* Maps original variable and equation indices into the block-decomposition order.
     Set by updateReverseVariableEquationOrderings() */
  vector<int> eq_idx_orig2block, endo_idx_orig2block;

  //! Vector describing equations: BlockSimulationType, if BlockSimulationType == EVALUATE_s then a expr_t on the new normalized equation
  equation_type_and_normalized_equation_t equation_type_and_normalized_equation;

  /* Stores derivatives of each block w.r.t. endogenous that belong to it.
     The tuple is: equation number (inside the block), variable number (inside
     the block), lead/lag */
  vector<map<tuple<int, int, int>, expr_t>> blocks_derivatives;

  class BlockInfo
  {
  public:
    BlockSimulationType simulation_type;
    int first_equation; // Stores a block-ordered equation ID
    int size{0};
    int mfs_size{0}; // Size of the minimal feedback set
    bool linear{true}; // Whether the block is linear in endogenous variable
    int n_static{0}, n_forward{0}, n_backward{0}, n_mixed{0};
    int max_endo_lag{0}, max_endo_lead{0}; // Maximum lag/lead on endos that appear in and *that belong to* the block
    int max_other_endo_lag{0}, max_other_endo_lead{0}; // Maximum lag/lead on endos that appear in but do not belong to the block
    int max_exo_lag{0}, max_exo_lead{0};
    int max_exo_det_lag{0}, max_exo_det_lead{0};
    int max_lag{0}, max_lead{0}; // The max over all endo/exo variables

    int
    getRecursiveSize() const
    {
      return size - mfs_size;
    };
  };

  // Whether block decomposition has been successfully computed
  bool block_decomposed {false};

  // Stores various informations on the blocks
  vector<BlockInfo> blocks;

  // Maps endogenous type-specific IDs to the block number to which it belongs
  vector<int> endo2block;
  /* Maps (original) equation number to the block number to which it belongs.
     It verifies: ∀i, eq2block[endo2eq[i]] = endo2block[i] */
  vector<int> eq2block;

  /* Temporary terms for block decomposed models.
     - the outer vector has as many elements as there are blocks in the model
     - the inner vector has as many elements as there are equations in the
       block, plus a last one which contains the temporary terms for
       derivatives

     It’s necessary to track temporary terms per equation, because some
     equations are evaluated instead of solved, and an equation E1 may depend
     on the value of an endogenous Y computed by a previously evaluated equation
     E2; in this case, if some temporary term TT of equation E2 contains Y,
     then TT needs to be computed after E1, but before E2. */
  vector<vector<temporary_terms_t>> blocks_temporary_terms;

  /* Stores, for each temporary term in block decomposed models, its index in
     the vector of all temporary terms */
  temporary_terms_idxs_t blocks_temporary_terms_idxs;

  //! Computes derivatives
  /*! \param order the derivation order
      \param vars the derivation IDs w.r.t. which compute the derivatives */
  void computeDerivatives(int order, const set<int> &vars);
  //! Computes derivatives of the Jacobian and Hessian w.r. to parameters
  void computeParamsDerivatives(int paramsDerivsOrder);
  //! Computes temporary terms (for all equations and derivatives)
  void computeTemporaryTerms(bool is_matlab, bool no_tmp_terms);
  //! Computes temporary terms per block
  void computeBlockTemporaryTerms(bool no_tmp_terms);

private:
  /* Add additional temporary terms for a given block. This method is called by
     computeBlockTemporaryTerms(). It does nothing by default, but is meant to
     be overriden by subclasses (actually by DynamicModel, who needs extra
     temporary terms for derivatives w.r.t. exogenous and other endogenous) */
  virtual void additionalBlockTemporaryTerms(int blk,
                                             vector<vector<temporary_terms_t>> &blocks_temporary_terms,
                                             map<expr_t, tuple<int, int, int>> &reference_count) const;

protected:
  //! Computes temporary terms for the file containing parameters derivatives
  void computeParamsDerivativesTemporaryTerms();
  //! Writes temporary terms
  template<ExprNodeOutputType output_type>
  void writeTemporaryTerms(const temporary_terms_t &tt, temporary_terms_t &temp_term_union, const temporary_terms_idxs_t &tt_idxs, ostream &output, deriv_node_temp_terms_t &tef_terms) const;
  void writeJsonTemporaryTerms(const temporary_terms_t &tt, temporary_terms_t &temp_term_union, ostream &output, deriv_node_temp_terms_t &tef_terms, const string &concat) const;
  //! Writes temporary terms in bytecode
  template<ExprNodeBytecodeOutputType output_type>
  void writeBytecodeTemporaryTerms(const temporary_terms_t &tt,
                                   temporary_terms_t &temporary_terms_union,
                                   BytecodeWriter &code_file,
                                   deriv_node_temp_terms_t &tef_terms) const;
  /* Adds information for (non-block) bytecode simulation in a separate .bin
     file.
     Returns the number of first derivatives w.r.t. endogenous variables */
  int writeBytecodeBinFile(const string &filename, bool is_two_boundaries) const;
  //! Adds per-block information for bytecode simulation in a separate .bin file
  int writeBlockBytecodeBinFile(ofstream &bin_file, int block) const;

  //! Fixes output when there are more than 32 nested parens, Issue #1201
  void fixNestedParenthesis(ostringstream &output, map<string, string> &tmp_paren_vars, bool &message_printed) const;
  //! Tests if string contains more than 32 nested parens, Issue #1201
  bool testNestedParenthesis(const string &str) const;

  //! Writes model equations
  template<ExprNodeOutputType output_type>
  void writeModelEquations(ostream &output, const temporary_terms_t &temporary_terms) const;

  // Returns outputs for derivatives and temporary terms at each derivation order
  template<ExprNodeOutputType output_type>
  pair<vector<ostringstream>, vector<ostringstream>> writeModelFileHelper() const;

  // Writes and compiles dynamic/static file (C version)
  template<bool dynamic>
  void writeModelCFile(const string &basename, const string &mexext, const filesystem::path &matlabroot, const filesystem::path &dynareroot) const;

  // Writes per-block residuals and temporary terms (incl. for derivatives)
  template<ExprNodeOutputType output_type>
  void writePerBlockHelper(int blk, ostream &output, temporary_terms_t &temporary_terms) const;

  /* Helper for writing derivatives w.r.t. parameters.
     Returns { tt, rp, gp, rpp, gpp, hp, g3p }.
     g3p is empty if requesting a static output type. */
  template<ExprNodeOutputType output_type>
  tuple<ostringstream, ostringstream, ostringstream, ostringstream,
        ostringstream, ostringstream, ostringstream> writeParamsDerivativesFileHelper() const;

  // Helper for writing bytecode (without block decomposition)
  template<bool dynamic>
  void writeBytecodeHelper(BytecodeWriter &code_file) const;

  // Helper for writing blocks in bytecode
  template<bool dynamic>
  void writeBlockBytecodeHelper(BytecodeWriter &code_file, int block) const;

  /* Write additional derivatives w.r.t. to exogenous, exogenous det and other endo
     in block+bytecode mode. Does nothing by default, but overriden by
     DynamicModel which needs those. */
  virtual void writeBlockBytecodeAdditionalDerivatives(BytecodeWriter &code_file, int block,
                                                       const temporary_terms_t &temporary_terms_union,
                                                       const deriv_node_temp_terms_t &tef_terms) const;

  // Helper for writing sparse derivatives indices in MATLAB/Octave driver file
  template<bool dynamic>
  void writeDriverSparseIndicesHelper(ostream &output) const;

  // Helper for writing sparse derivatives indices in JSON
  template<bool dynamic>
  void writeJsonSparseIndicesHelper(ostream &output) const;

  /* Helper for writing JSON output for residuals and derivatives.
     Returns mlv and derivatives output at each derivation order. */
  template<bool dynamic>
  pair<ostringstream, vector<ostringstream>> writeJsonComputingPassOutputHelper(bool writeDetails) const;

  /* Helper for writing JSON derivatives w.r.t. parameters.
     Returns { mlv, tt, rp, gp, rpp, gpp, hp, g3p }.
     g3p is empty if requesting a static output type. */
  template<bool dynamic>
  tuple<ostringstream, ostringstream, ostringstream, ostringstream, ostringstream,
        ostringstream, ostringstream, ostringstream> writeJsonParamsDerivativesHelper(bool writeDetails) const;

  //! Writes JSON model equations
  //! if residuals = true, we are writing the dynamic/static model.
  //! Otherwise, just the model equations (with line numbers, no tmp terms)
  void writeJsonModelEquations(ostream &output, bool residuals) const;
  /* Writes JSON model local variables.
     Optionally put the external function variable calls into TEF terms */
  void writeJsonModelLocalVariables(ostream &output, bool write_tef_terms, deriv_node_temp_terms_t &tef_terms) const;

  //! Writes model equations in bytecode
  template<ExprNodeBytecodeOutputType output_type>
  void writeBytecodeModelEquations(BytecodeWriter &code_file, const temporary_terms_t &temporary_terms, const deriv_node_temp_terms_t &tef_terms) const;

  // Writes the sparse representation of the model in Julia
  // Assumes that the directory <MODFILE>/model/julia/ already exists
  template<bool dynamic>
  void writeSparseModelJuliaFiles(const string &basename) const;

  //! Writes LaTeX model file
  void writeLatexModelFile(const string &mod_basename, const string &latex_basename, ExprNodeOutputType output_type, bool write_equation_tags) const;

private:
  //! Sparse matrix of double to store the values of the static Jacobian
  /*! First index is equation number, second index is endogenous type specific ID */
  using jacob_map_t = map<pair<int, int>, double>;

  //! Normalization of equations, as computed by computeNonSingularNormalization()
  /*! Maps endogenous type specific IDs to equation numbers */
  vector<int> endo2eq;

  // Stores workers used for compiling MEX files in parallel
  static vector<jthread> mex_compilation_workers;

  /* The following variables implement the thread synchronization mechanism for
     limiting the number of concurrent GCC processes and tracking dependencies
     between object files. */
  static condition_variable_any mex_compilation_cv;
  static mutex mex_compilation_mut;
  /* Object/MEX files waiting to be compiled (with their prerequisites as 2nd
     element and compilation command as the 3rd element) */
  static vector<tuple<filesystem::path, set<filesystem::path>, string>> mex_compilation_queue;
  // Object/MEX files already compiled
  static set<filesystem::path> mex_compilation_done;

  /* Compute a pseudo-Jacobian whose all elements are either zero or one,
     depending on whether the variable symbolically appears in the equation */
  jacob_map_t computeSymbolicJacobian() const;

  // Compute {var,eq}_idx_orig2block from {var,eq}_idx_block2orig
  void updateReverseVariableEquationOrderings();

  //! Compute the matching between endogenous and variable using the jacobian contemporaneous_jacobian
  /*!
    \param contemporaneous_jacobian Jacobian used as an incidence matrix: all elements declared in the map (even if they are zero), are used as vertices of the incidence matrix
    \return True if a complete normalization has been achieved
  */
  bool computeNormalization(const jacob_map_t &contemporaneous_jacobian, bool verbose);

  //! Try to compute the matching between endogenous and variable using a decreasing cutoff
  /*!
    Applied to the jacobian contemporaneous_jacobian and stop when a matching is found.
    If no matching is found using a strictly positive cutoff, then a zero cutoff is applied (i.e. use a symbolic normalization); in that case, the method adds zeros in the jacobian matrices to reflect all the edges in the symbolic incidence matrix.
    If no matching is found with a zero cutoff, an error message is printed.
    The resulting normalization is stored in endo2eq.
    Returns a boolean indicating success.
  */
  bool computeNonSingularNormalization(const jacob_map_t &contemporaneous_jacobian);
  //! Evaluate the jacobian (w.r.t. endogenous) and suppress all the elements below the cutoff
  /*! Returns the contemporaneous_jacobian.
      Elements below the cutoff are discarded. External functions are evaluated to 1. */
  jacob_map_t evaluateAndReduceJacobian(const eval_context_t &eval_context) const;
  /* Search the equations and variables belonging to the prologue and the
     epilogue of the model.
     Initializes “eq_idx_block2orig” and “endo_idx_block2orig”.
     Returns the sizes of the prologue and epilogue. */
  pair<int, int> computePrologueAndEpilogue();
  //! Determine the type of each equation of model and try to normalize the unnormalized equation
  void equationTypeDetermination(const map<tuple<int, int, int>, expr_t> &first_order_endo_derivatives, int mfs);
  /* Fills the max lags/leads and n_{static,mixed,forward,backward} fields of a
     given block.
     Needs the fields size and first_equation. */
  void computeDynamicStructureOfBlock(int blk);
  /* Fills the simulation_type field of a given block.
     Needs the fields size, max_endo_lag and max_endo_lead. */
  void computeSimulationTypeOfBlock(int blk);
  /* Compute the block decomposition and for a non-recusive block find the minimum feedback set

     Initializes the “blocks”, “endo2block” and “eq2block” structures. */
  void computeBlockDecomposition(int prologue, int epilogue);
  /* Reduce the number of block by merging the same type of equations in the
     prologue and the epilogue */
  void reduceBlockDecomposition();
  /* The 1st output gives, for each equation (in original order) the (max_lag,
     max_lead) across all endogenous that appear in the equation and that
     belong to the same block (i.e. those endogenous are solved in the same
     block).

     The 2nd output gives, for each type-specific endo IDs, its (max_lag,
     max_lead) across all its occurences inside the equations of the block to
     which it belongs. */
  pair<lag_lead_vector_t, lag_lead_vector_t> getVariableLeadLagByBlock() const;
  //! Print an abstract of the block structure of the model
  void printBlockDecomposition() const;
  //! Determine for each block if it is linear or not
  void determineLinearBlocks();

protected:
  //! Return the type of equation belonging to the block
  EquationType
  getBlockEquationType(int blk, int eq) const
  {
    return equation_type_and_normalized_equation[eq_idx_block2orig[blocks[blk].first_equation+eq]].first;
  };
  //! Return true if the equation has been normalized
  bool
  isBlockEquationRenormalized(int blk, int eq) const
  {
    return equation_type_and_normalized_equation[eq_idx_block2orig[blocks[blk].first_equation + eq]].first == EquationType::evaluateRenormalized;
  };
  //! Return the expr_t of equation belonging to the block
  BinaryOpNode *
  getBlockEquationExpr(int blk, int eq) const
  {
    return equations[eq_idx_block2orig[blocks[blk].first_equation + eq]];
  };
  //! Return the expr_t of renormalized equation belonging to the block
  BinaryOpNode *
  getBlockEquationRenormalizedExpr(int blk, int eq) const
  {
    return equation_type_and_normalized_equation[eq_idx_block2orig[blocks[blk].first_equation + eq]].second;
  };
  //! Return the original number of equation belonging to the block
  int
  getBlockEquationID(int blk, int eq) const
  {
    return eq_idx_block2orig[blocks[blk].first_equation + eq];
  };
  //! Return the original number of variable belonging to the block
  int
  getBlockVariableID(int blk, int var) const
  {
    return endo_idx_block2orig[blocks[blk].first_equation + var];
  };
  //! Return the position of an equation (given by its original index) inside its block
  int
  getBlockInitialEquationID(int blk, int eq) const
  {
    return eq_idx_orig2block[eq] - blocks[blk].first_equation;
  };
  //! Return the position of a variable (given by its original index) inside its block
  int
  getBlockInitialVariableID(int blk, int var) const
  {
    return endo_idx_orig2block[var] - blocks[blk].first_equation;
  };
  //! Initialize equation_reordered & variable_reordered
  void initializeVariablesAndEquations();

private:
  //! Returns the 1st derivatives w.r.t. endogenous in a different format
  /*! Returns a map (equation, type-specific ID, lag) → derivative.
      Assumes that derivatives have already been computed. */
  map<tuple<int, int, int>, expr_t> collectFirstOrderDerivativesEndogenous();

protected:
  //! Computes chain rule derivatives of the Jacobian w.r. to endogenous variables
  virtual void computeChainRuleJacobian() = 0;

  /* Compute block decomposition, its derivatives and temporary terms. Meant to
     be overriden in derived classes which don’t support block decomposition
     (currently Epilogue and PlannerObjective). Sets “block_decomposed” to true
     in case of success. */
  virtual void computingPassBlock(const eval_context_t &eval_context, bool no_tmp_terms);

  /* Get column number within Jacobian of a given block.
     “var” is the block-specific endogenous variable index. */
  virtual int getBlockJacobianEndoCol(int blk, int var, int lag) const = 0;

  // Returns a human-readable string describing the model class (e.g. “dynamic model”…)
  virtual string modelClassName() const = 0;

  /* Given a sparse matrix in column major order, returns the colptr pointer for
     the CSC storage */
  static vector<int> computeCSCColPtr(const SparseColumnMajorOrderMatrix &matrix, int ncols);

private:
  //! Internal helper for the copy constructor and assignment operator
  /*! Copies all the structures that contain ExprNode*, by the converting the
      pointers into their equivalent in the new tree */
  void copyHelper(const ModelTree &m);
  //! Returns the name of the MATLAB architecture given the extension used for MEX files
  static string matlab_arch(const string &mexext);
#ifdef __APPLE__
  //! Finds a suitable GCC compiler on macOS
  static string findGccOnMacos(const string &mexext);
#endif
  /* Compiles a MEX file (if link=true) or an object file to be linked later
     into a MEX file (if link=false). The compilation is done in separate
     worker threads working in parallel, so the call to this function is not
     blocking. The dependency of a linked MEX file upon intermediary objects is
     nicely handled. Returns the name of the output file (to be reused later as
     input file if link=false). */
  filesystem::path compileMEX(const filesystem::path &output_dir, const string &output_basename, const string &mexext, const vector<filesystem::path> &input_files, const filesystem::path &matlabroot, const filesystem::path &dynareroot, bool link = true) const;

public:
  ModelTree(SymbolTable &symbol_table_arg,
            NumericalConstants &num_constants_arg,
            ExternalFunctionsTable &external_functions_table_arg,
            bool is_dynamic_arg = false);

protected:
  ModelTree(const ModelTree &m);
  ModelTree &operator=(const ModelTree &m);

public:
  //! Absolute value under which a number is considered to be zero
  double cutoff{1e-15};
  //! Compute the minimum feedback set
  /*!   0 : all endogenous variables are considered as feedback variables
    1 : the variables belonging to non normalized equation are considered as feedback variables
    2 : the variables belonging to a non linear equation are considered as feedback variables
    3 : the variables belonging to a non normalizable non linear equation are considered as feedback variables
    default value = 0 */
  int mfs{0};
  //! Declare a node as an equation of the model; also give its line number
  void addEquation(expr_t eq, optional<int> lineno);
  //! Declare a node as an equation of the model, also giving its tags
  void addEquation(expr_t eq, optional<int> lineno, const map<string, string> &eq_tags);
  //! Declare a node as an auxiliary equation of the model, adding it at the end of the list of auxiliary equations
  void addAuxEquation(expr_t eq);
  //! Returns the number of equations in the model
  int equation_number() const;
  //! Adds a trend variable with its growth factor
  void addTrendVariables(const vector<int> &trend_vars, expr_t growth_factor) noexcept(false);
  //! Adds a nonstationary variables with their (common) deflator
  void addNonstationaryVariables(const vector<int> &nonstationary_vars, bool log_deflator, expr_t deflator) noexcept(false);
  //! Is a given variable non-stationary?
  bool isNonstationary(int symb_id) const;
  void set_cutoff_to_zero();
  /*! Reorder auxiliary variables so that they appear in recursive order in
      set_auxiliary_variables.m and dynamic_set_auxiliary_series.m */
  void reorderAuxiliaryEquations();
  //! Find equations of the form “variable=constant”, excluding equations with “mcp” tag (see dynare#1697)
  void findConstantEquationsWithoutMcpTag(map<VariableNode *, NumConstNode *> &subst_table) const;
  /* Given an expression, searches for the first equation that has exactly this
     expression on the LHS, and returns the RHS of that equation.
     If no such equation can be found, throws an ExprNode::MatchFailureExpression */
  expr_t getRHSFromLHS(expr_t lhs) const;

  // Initialize the MEX compilation workers
  static void initializeMEXCompilationWorkers(int numworkers);

  // Waits until the MEX compilation queue is empty
  static void waitForMEXCompilationWorkers();

  //! Returns all the equation tags associated to an equation
  map<string, string>
  getEquationTags(int eq) const
  {
    return equation_tags.getTagsByEqn(eq);
  }

  //! Returns the vector of non-zero derivative counts
  const vector<int> &
  getNNZDerivatives() const
  {
    return NNZDerivatives;
  }

  //! Returns the vector of temporary terms derivatives
  const vector<temporary_terms_t> &
  getTemporaryTermsDerivatives() const
  {
    return temporary_terms_derivatives;
  }

  //!Returns the maximum order of computed derivatives
  int
  getComputedDerivsOrder() const
  {
    return computed_derivs_order;
  }

  static string
  BlockSim(BlockSimulationType type)
  {
    switch (type)
      {
      case BlockSimulationType::evaluateForward:
        return "EVALUATE FORWARD             ";
      case BlockSimulationType::evaluateBackward:
        return "EVALUATE BACKWARD            ";
      case BlockSimulationType::solveForwardSimple:
        return "SOLVE FORWARD SIMPLE         ";
      case BlockSimulationType::solveBackwardSimple:
        return "SOLVE BACKWARD SIMPLE        ";
      case BlockSimulationType::solveTwoBoundariesSimple:
        return "SOLVE TWO BOUNDARIES SIMPLE  ";
      case BlockSimulationType::solveForwardComplete:
        return "SOLVE FORWARD COMPLETE       ";
      case BlockSimulationType::solveBackwardComplete:
        return "SOLVE BACKWARD COMPLETE      ";
      case BlockSimulationType::solveTwoBoundariesComplete:
        return "SOLVE TWO BOUNDARIES COMPLETE";
      default:
        return "UNKNOWN                      ";
      }
  }
};

template<ExprNodeOutputType output_type>
void
ModelTree::writeTemporaryTerms(const temporary_terms_t &tt,
                               temporary_terms_t &temp_term_union,
                               const temporary_terms_idxs_t &tt_idxs,
                               ostream &output, deriv_node_temp_terms_t &tef_terms) const
{
  for (auto it : tt)
    {
      if (dynamic_cast<AbstractExternalFunctionNode *>(it))
        it->writeExternalFunctionOutput(output, output_type, temp_term_union, tt_idxs, tef_terms);

      it->writeOutput(output, output_type, tt, tt_idxs, tef_terms);
      output << " = ";
      it->writeOutput(output, output_type, temp_term_union, tt_idxs, tef_terms);

      if constexpr(isCOutput(output_type) || isMatlabOutput(output_type))
        output << ";";
      output << endl;

      temp_term_union.insert(it);
    }
}

template<ExprNodeOutputType output_type>
void
ModelTree::writeModelEquations(ostream &output, const temporary_terms_t &temporary_terms) const
{
  for (int eq {0}; eq < static_cast<int>(equations.size()); eq++)
    {
      BinaryOpNode *eq_node { equations[eq] };
      expr_t lhs { eq_node->arg1 }, rhs { eq_node->arg2 };

      // Test if the right hand side of the equation is empty.
      double vrhs {1.0};
      try
        {
          vrhs = rhs->eval({});
        }
      catch (ExprNode::EvalException &e)
        {
        }

      if (vrhs != 0) // The right hand side of the equation is not empty ==> residual=lhs-rhs;
        if constexpr(isJuliaOutput(output_type))
          {
            output << "    residual" << LEFT_ARRAY_SUBSCRIPT(output_type)
                   << eq + ARRAY_SUBSCRIPT_OFFSET(output_type)
                   << RIGHT_ARRAY_SUBSCRIPT(output_type)
                   << " = (";
            lhs->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs);
            output << ") - (";
            rhs->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs);
            output << ")" << endl;
          }
        else
          {
            output << "lhs = ";
            lhs->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs);
            output << ";" << endl
                   << "rhs = ";
            rhs->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs);
            output << ";" << endl
                   << "residual" << LEFT_ARRAY_SUBSCRIPT(output_type)
                   << eq + ARRAY_SUBSCRIPT_OFFSET(output_type)
                   << RIGHT_ARRAY_SUBSCRIPT(output_type)
                   << " = lhs - rhs;" << endl;
          }
      else // The right hand side of the equation is empty ==> residual=lhs;
        {
          output << "residual" << LEFT_ARRAY_SUBSCRIPT(output_type)
                 << eq + ARRAY_SUBSCRIPT_OFFSET(output_type)
                 << RIGHT_ARRAY_SUBSCRIPT(output_type)
                 << " = ";
          lhs->writeOutput(output, output_type, temporary_terms, temporary_terms_idxs);
          output << ";" << endl;
        }
    }
}

template<ExprNodeOutputType output_type>
pair<vector<ostringstream>, vector<ostringstream>>
ModelTree::writeModelFileHelper() const
{
  constexpr bool sparse {isSparseModelOutput(output_type)};

  vector<ostringstream> d_output(derivatives.size()); // Derivatives output (at all orders, including 0=residual)
  vector<ostringstream> tt_output(derivatives.size()); // Temp terms output (at all orders)

  deriv_node_temp_terms_t tef_terms;
  temporary_terms_t temp_term_union;

  writeTemporaryTerms<output_type>(temporary_terms_derivatives[0], temp_term_union,
                                   temporary_terms_idxs, tt_output[0], tef_terms);

  writeModelEquations<output_type>(d_output[0], temp_term_union);

  // Writing Jacobian
  if (!derivatives[1].empty())
    {
      writeTemporaryTerms<output_type>(temporary_terms_derivatives[1], temp_term_union,
                                       temporary_terms_idxs, tt_output[1], tef_terms);

      if constexpr(sparse)
        {
          // NB: we iterate over the Jacobian reordered in column-major order
          // Indices of rows and columns are output in M_ and the JSON file (since they are constant)
          for (int k {0};
               const auto &[row_col, d1] : jacobian_sparse_column_major_order)
            {
              d_output[1] << "g1_v" << LEFT_ARRAY_SUBSCRIPT(output_type)
                          << k + ARRAY_SUBSCRIPT_OFFSET(output_type)
                          << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=";
              d1->writeOutput(d_output[1], output_type, temp_term_union, temporary_terms_idxs, tef_terms);
              d_output[1] << ";" << endl;
              k++;
            }
        }
      else // Legacy representation (dense matrix)
        {
          for (const auto &[indices, d1] : derivatives[1])
            {
              auto [eq, var] = vectorToTuple<2>(indices);

              d_output[1] << "g1" << LEFT_ARRAY_SUBSCRIPT(output_type);
              if constexpr(isMatlabOutput(output_type) || isJuliaOutput(output_type))
                d_output[1] << eq + 1 << "," << getJacobianCol(var, sparse) + 1;
              else
                d_output[1] << eq + getJacobianCol(var, sparse)*equations.size();
              d_output[1] << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=";
              d1->writeOutput(d_output[1], output_type,
                              temp_term_union, temporary_terms_idxs, tef_terms);
              d_output[1] << ";" << endl;
            }
        }
    }

  // Write derivatives for order ≥ 2
  for (size_t i = 2; i < derivatives.size(); i++)
    if (!derivatives[i].empty())
      {
        writeTemporaryTerms<output_type>(temporary_terms_derivatives[i], temp_term_union,
                                         temporary_terms_idxs, tt_output[i], tef_terms);

        if constexpr(sparse)
          {
            /* List non-zero elements of the tensor in row-major order (this is
               suitable for the k-order solver according to Normann). */
            // Tensor indices are output in M_ and the JSON file (since they are constant)
            for (int k {0};
                 const auto &[vidx, d] : derivatives[i])
              {
                d_output[i] << "g" << i << "_v" << LEFT_ARRAY_SUBSCRIPT(output_type)
                            << k + ARRAY_SUBSCRIPT_OFFSET(output_type)
                            << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=";
                d->writeOutput(d_output[i], output_type, temp_term_union, temporary_terms_idxs, tef_terms);
                d_output[i] << ";" << endl;
                k++;
              }
          }
        else // Legacy representation
          {
            /* When creating the sparse matrix (in MATLAB or C mode), since storage
               is in column-major order, output the first column, then the second,
               then the third. This gives a significant performance boost in use_dll
               mode (at both compilation and runtime), because it facilitates memory
               accesses and expression reusage. */
            ostringstream i_output, j_output, v_output;

            for (int k{0}; // Current line index in the 3-column matrix
                 const auto &[vidx, d] : derivatives[i])
              {
                int eq{vidx[0]};

                int col_idx{0};
                for (size_t j = 1; j < vidx.size(); j++)
                  {
                    col_idx *= getJacobianColsNbr(sparse);
                    col_idx += getJacobianCol(vidx[j], sparse);
                  }

                if constexpr(isJuliaOutput(output_type))
                  {
                    d_output[i] << "    g" << i << "[" << eq + 1 << "," << col_idx + 1 << "] = ";
                    d->writeOutput(d_output[i], output_type, temp_term_union, temporary_terms_idxs, tef_terms);
                    d_output[i] << endl;
                  }
                else
                  {
                    i_output << "g" << i << "_i" << LEFT_ARRAY_SUBSCRIPT(output_type)
                             << k + ARRAY_SUBSCRIPT_OFFSET(output_type)
                             << RIGHT_ARRAY_SUBSCRIPT(output_type)
                             << "=" << eq + 1 << ";" << endl;
                    j_output << "g" << i << "_j" << LEFT_ARRAY_SUBSCRIPT(output_type)
                             << k + ARRAY_SUBSCRIPT_OFFSET(output_type)
                             << RIGHT_ARRAY_SUBSCRIPT(output_type)
                             << "=" << col_idx + 1 << ";" << endl;
                    v_output << "g" << i << "_v" << LEFT_ARRAY_SUBSCRIPT(output_type)
                             << k + ARRAY_SUBSCRIPT_OFFSET(output_type)
                             << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=";
                    d->writeOutput(v_output, output_type, temp_term_union, temporary_terms_idxs, tef_terms);
                    v_output << ";" << endl;

                    k++;
                  }

                // Output symetric elements at order 2
                if (i == 2 && vidx[1] != vidx[2])
                  {
                    int col_idx_sym{getJacobianCol(vidx[2], sparse) * getJacobianColsNbr(sparse) + getJacobianCol(vidx[1], sparse)};

                    if constexpr(isJuliaOutput(output_type))
                      d_output[2] << "    g2[" << eq + 1 << "," << col_idx_sym + 1 << "] = "
                                  << "g2[" << eq + 1 << "," << col_idx + 1 << "]" << endl;
                    else
                      {
                        i_output << "g" << i << "_i" << LEFT_ARRAY_SUBSCRIPT(output_type)
                                 << k + ARRAY_SUBSCRIPT_OFFSET(output_type)
                                 << RIGHT_ARRAY_SUBSCRIPT(output_type)
                                 << "=" << eq + 1 << ";" << endl;
                        j_output << "g" << i << "_j" << LEFT_ARRAY_SUBSCRIPT(output_type)
                                 << k + ARRAY_SUBSCRIPT_OFFSET(output_type)
                                 << RIGHT_ARRAY_SUBSCRIPT(output_type)
                                 << "=" << col_idx_sym + 1 << ";" << endl;
                        v_output << "g" << i << "_v" << LEFT_ARRAY_SUBSCRIPT(output_type)
                                 << k + ARRAY_SUBSCRIPT_OFFSET(output_type)
                                 << RIGHT_ARRAY_SUBSCRIPT(output_type) << "="
                                 << "g" << i << "_v" << LEFT_ARRAY_SUBSCRIPT(output_type)
                                 << k-1 + ARRAY_SUBSCRIPT_OFFSET(output_type)
                                 << RIGHT_ARRAY_SUBSCRIPT(output_type) << ";" << endl;

                        k++;
                      }
                  }
              }
            if constexpr(!isJuliaOutput(output_type))
              d_output[i] << i_output.str() << j_output.str() << v_output.str();
          }
      }

  if constexpr(isMatlabOutput(output_type))
    {
      // Check that we don't have more than 32 nested parenthesis because MATLAB does not suppor this. See Issue #1201
      map<string, string> tmp_paren_vars;
      bool message_printed {false};
      for (auto &it : tt_output)
        fixNestedParenthesis(it, tmp_paren_vars, message_printed);
      for (auto &it : d_output)
        fixNestedParenthesis(it, tmp_paren_vars, message_printed);
    }

  return { move(d_output), move(tt_output) };
}

template<bool dynamic>
void
ModelTree::writeModelCFile(const string &basename, const string &mexext,
                           const filesystem::path &matlabroot,
                           const filesystem::path &dynareroot) const
{
  ofstream output;
  auto open_file = [&output](const filesystem::path &p)
  {
    output.open(p, ios::out | ios::binary);
    if (!output.is_open())
      {
        cerr << "ERROR: Can't open file " << p.string() << " for writing" << endl;
        exit(EXIT_FAILURE);
      }
  };

  const filesystem::path model_src_dir { filesystem::path{basename} / "model" / "src" };

  auto [d_output, tt_output] = writeModelFileHelper<dynamic ? ExprNodeOutputType::CDynamicModel : ExprNodeOutputType::CStaticModel>();
  vector<filesystem::path> header_files, object_files;

  // TODO: when C++20 support is complete, mark the following strings constexpr
  const string prefix { dynamic ? "dynamic_" : "static_" };
  const string ss_it_argin { dynamic ? ", const double *restrict steady_state, int it_" : "" };
  const string ss_it_argout { dynamic ? ", steady_state, it_" : "" };
  const string nb_row_x_argin { dynamic ? ", int nb_row_x" : "" };
  const string nb_row_x_argout { dynamic ? ", nb_row_x" : "" };

  for (size_t i {0}; i < d_output.size(); i++)
    {
      const string funcname { prefix + (i == 0 ? "resid" : "g" + to_string(i))};

      const string prototype_tt { "void " + funcname + "_tt(const double *restrict y, const double *restrict x" + nb_row_x_argin + ", const double *restrict params" + ss_it_argin + ", double *restrict T)" };

      const filesystem::path header_tt { model_src_dir / (funcname + "_tt.h") };
      open_file(header_tt);
      output << prototype_tt << ";" << endl;
      output.close();
      header_files.push_back(header_tt);

      const filesystem::path source_tt { model_src_dir / (funcname + "_tt.c") };
      open_file(source_tt);
      output << "#include <math.h>" << endl
             << R"(#include "mex.h")" << endl // Needed for calls to external functions
             << endl;
      writePowerDerivHeader(output);
      output << endl
             << prototype_tt << endl
             << "{" << endl
             << tt_output[i].str()
             << "}" << endl
             << endl;
      output.close();
      object_files.push_back(compileMEX(model_src_dir, funcname + "_tt" , mexext, { source_tt },
                                        matlabroot, dynareroot, false));

      const string prototype_main
        {
          [&funcname, &ss_it_argin, &nb_row_x_argin, i]
          {
            string p = "void " + funcname + "(const double *restrict y, const double *restrict x" + nb_row_x_argin + ", const double *restrict params" + ss_it_argin + ", const double *restrict T, ";
            if (i == 0)
              p += "double *restrict residual";
            else if (i == 1)
              p += "double *restrict g1";
            else
              p += "double *restrict g" + to_string(i) + "_i, double *restrict g" +
                to_string(i) + "_j, double *restrict g" + to_string(i) + "_v";
            p += ")";
            return p;
          }()
        };

      const filesystem::path header_main { model_src_dir / (funcname + ".h") };
      open_file(header_main);
      output << prototype_main << ";" << endl;
      output.close();
      header_files.push_back(header_main);

      const filesystem::path source_main { model_src_dir / (funcname + ".c") };
      open_file(source_main);
      output << "#include <math.h>" << endl
             << R"(#include "mex.h")" << endl // Needed for calls to external functions
             << endl;
      writePowerDerivHeader(output);
      output << endl
             << prototype_main << endl
             << "{" << endl;
      if (i == 0)
        output << "  double lhs, rhs;" << endl;
      output << d_output[i].str()
             << "}" << endl
             << endl;
      output.close();
      object_files.push_back(compileMEX(model_src_dir, funcname, mexext, { source_main },
                                        matlabroot, dynareroot, false));
    }

  const filesystem::path filename { model_src_dir / (dynamic ? "dynamic.c" : "static.c") };

  const int ntt { static_cast<int>(temporary_terms_derivatives[0].size() + temporary_terms_derivatives[1].size() + temporary_terms_derivatives[2].size() + temporary_terms_derivatives[3].size()) };

  open_file(filename);
  output << "/*" << endl
         << " * " << filename << " : Computes " << modelClassName() << " for Dynare" << endl
         << " *" << endl
         << " * Warning : this file is generated automatically by Dynare" << endl
         << " *           from model file (.mod)" << endl
         << " */" << endl
         << endl
         << "#include <math.h>" << endl // Needed for getPowerDeriv()
         << "#include <stdlib.h>" << endl // Needed for malloc() and free()
         << R"(#include "mex.h")" << endl;
  for (const auto &it : header_files)
    output << "#include " << it.filename() << endl;
  output << endl;

  // Write function definition if BinaryOpcode::powerDeriv is used
  writePowerDeriv(output);

  output << endl
         << "void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])" << endl
         << "{" << endl;
  if constexpr(dynamic)
    output << "  if (nlhs > " << min(computed_derivs_order + 1, 4) << ")" << endl
           << R"(    mexErrMsgTxt("Derivatives of higher order than computed have been requested");)" << endl
           << "  if (nrhs != 5)" << endl
           << R"(    mexErrMsgTxt("Requires exactly 5 input arguments");)" << endl;
  else
    output << "  if (nrhs > 3)" << endl
           << R"(    mexErrMsgTxt("Accepts at most 3 output arguments");)" << endl
           << "  if (nrhs != 3)" << endl
           << R"(    mexErrMsgTxt("Requires exactly 3 input arguments");)" << endl;
  output << endl
         << "  double *y = mxGetPr(prhs[0]);" << endl
         << "  double *x = mxGetPr(prhs[1]);" << endl
         << "  double *params = mxGetPr(prhs[2]);" << endl;
  if constexpr(dynamic)
    output << "  double *steady_state = mxGetPr(prhs[3]);" << endl
           << "  int it_ = (int) mxGetScalar(prhs[4]) - 1;" << endl
           << "  int nb_row_x = mxGetM(prhs[1]);" << endl;
  output << endl
         << "  double *T = (double *) malloc(sizeof(double)*" << ntt << ");" << endl
         << endl
         << "  if (nlhs >= 1)" << endl
         << "    {" << endl
         << "       plhs[0] = mxCreateDoubleMatrix(" << equations.size() << ",1, mxREAL);" << endl
         << "       double *residual = mxGetPr(plhs[0]);" << endl
         << "       " << prefix << "resid_tt(y, x" << nb_row_x_argout << ", params" << ss_it_argout << ", T);" << endl
         << "       " << prefix << "resid(y, x" << nb_row_x_argout << ", params" << ss_it_argout << ", T, residual);" << endl
         << "    }" << endl
         << endl
         << "  if (nlhs >= 2)" << endl
         << "    {" << endl
         << "       plhs[1] = mxCreateDoubleMatrix(" << equations.size() << ", " << getJacobianColsNbr(false) << ", mxREAL);" << endl
         << "       double *g1 = mxGetPr(plhs[1]);" << endl
         << "       " << prefix << "g1_tt(y, x" << nb_row_x_argout << ", params" << ss_it_argout << ", T);" << endl
         << "       " << prefix << "g1(y, x" << nb_row_x_argout << ", params" << ss_it_argout << ", T, g1);" << endl
         << "    }" << endl
         << endl
         << "  if (nlhs >= 3)" << endl
         << "    {" << endl
         << "      mxArray *g2_i = mxCreateDoubleMatrix(" << NNZDerivatives[2] << ", " << 1 << ", mxREAL);" << endl
         << "      mxArray *g2_j = mxCreateDoubleMatrix(" << NNZDerivatives[2] << ", " << 1 << ", mxREAL);" << endl
         << "      mxArray *g2_v = mxCreateDoubleMatrix(" << NNZDerivatives[2] << ", " << 1 << ", mxREAL);" << endl
         << "      " << prefix << "g2_tt(y, x" << nb_row_x_argout << ", params" << ss_it_argout << ", T);" << endl
         << "      " << prefix << "g2(y, x" << nb_row_x_argout << ", params" << ss_it_argout << ", T, mxGetPr(g2_i), mxGetPr(g2_j), mxGetPr(g2_v));" << endl
         << "      mxArray *m = mxCreateDoubleScalar(" << equations.size() << ");" << endl
         << "      mxArray *n = mxCreateDoubleScalar(" << getJacobianColsNbr(false)*getJacobianColsNbr(false) << ");" << endl
         << "      mxArray *plhs_sparse[1], *prhs_sparse[5] = { g2_i, g2_j, g2_v, m, n };" << endl
         << R"(      mexCallMATLAB(1, plhs_sparse, 5, prhs_sparse, "sparse");)" << endl
         << "      plhs[2] = plhs_sparse[0];" << endl
         << "      mxDestroyArray(g2_i);" << endl
         << "      mxDestroyArray(g2_j);" << endl
         << "      mxDestroyArray(g2_v);" << endl
         << "      mxDestroyArray(m);" << endl
         << "      mxDestroyArray(n);" << endl
         << "    }" << endl
         << endl;
  if constexpr(dynamic)
    output << "  if (nlhs >= 4)" << endl
           << "    {" << endl
           << "      mxArray *g3_i = mxCreateDoubleMatrix(" << NNZDerivatives[3] << ", " << 1 << ", mxREAL);" << endl
           << "      mxArray *g3_j = mxCreateDoubleMatrix(" << NNZDerivatives[3] << ", " << 1 << ", mxREAL);" << endl
           << "      mxArray *g3_v = mxCreateDoubleMatrix(" << NNZDerivatives[3] << ", " << 1 << ", mxREAL);" << endl
           << "      " << prefix << "g3_tt(y, x" << nb_row_x_argout << ", params" << ss_it_argout << ", T);" << endl
           << "      " << prefix << "g3(y, x" << nb_row_x_argout << ", params" << ss_it_argout << ", T, mxGetPr(g3_i), mxGetPr(g3_j), mxGetPr(g3_v));" << endl
           << "      mxArray *m = mxCreateDoubleScalar(" << equations.size() << ");" << endl
           << "      mxArray *n = mxCreateDoubleScalar(" << getJacobianColsNbr(false)*getJacobianColsNbr(false)*getJacobianColsNbr(false) << ");" << endl
           << "      mxArray *plhs_sparse[1], *prhs_sparse[5] = { g3_i, g3_j, g3_v, m, n };" << endl
           << R"(      mexCallMATLAB(1, plhs_sparse, 5, prhs_sparse, "sparse");)" << endl
           << "      plhs[3] = plhs_sparse[0];" << endl
           << "      mxDestroyArray(g3_i);" << endl
           << "      mxDestroyArray(g3_j);" << endl
           << "      mxDestroyArray(g3_v);" << endl
           << "      mxDestroyArray(m);" << endl
           << "      mxDestroyArray(n);" << endl
           << "    }" << endl
           << endl;

  output << "  free(T);" << endl
         << "}" << endl;
  output.close();

  object_files.push_back(filename);
  compileMEX(packageDir(basename), dynamic ? "dynamic" : "static", mexext, object_files, matlabroot,
             dynareroot);
}

template<ExprNodeOutputType output_type>
void
ModelTree::writePerBlockHelper(int blk, ostream &output, temporary_terms_t &temporary_terms) const
{
  int block_recursive_size { blocks[blk].getRecursiveSize() };

  // The equations
  deriv_node_temp_terms_t tef_terms;

  auto write_eq_tt = [&](int eq)
                     {
                       for (auto it : blocks_temporary_terms[blk][eq])
                         {
                           if (dynamic_cast<AbstractExternalFunctionNode *>(it))
                             it->writeExternalFunctionOutput(output, output_type, temporary_terms, blocks_temporary_terms_idxs, tef_terms);

                           output << "  ";
                           it->writeOutput(output, output_type, blocks_temporary_terms[blk][eq], blocks_temporary_terms_idxs, tef_terms);
                           output << '=';
                           it->writeOutput(output, output_type, temporary_terms, blocks_temporary_terms_idxs, tef_terms);
                           temporary_terms.insert(it);
                           output << ';' << endl;
                         }
                     };

  for (int eq {0}; eq < blocks[blk].size; eq++)
    {
      write_eq_tt(eq);

      EquationType equ_type { getBlockEquationType(blk, eq) };
      BinaryOpNode *e { getBlockEquationExpr(blk, eq) };
      expr_t lhs { e->arg1 }, rhs { e->arg2 };
      switch (blocks[blk].simulation_type)
        {
        case BlockSimulationType::evaluateBackward:
        case BlockSimulationType::evaluateForward:
          evaluation:
          if (equ_type == EquationType::evaluateRenormalized)
            {
              e = getBlockEquationRenormalizedExpr(blk, eq);
              lhs = e->arg1;
              rhs = e->arg2;
            }
          else if (equ_type != EquationType::evaluate)
            {
              cerr << "Type mismatch for equation " << getBlockEquationID(blk, eq)+1  << endl;
              exit(EXIT_FAILURE);
            }
          output << "  ";
          lhs->writeOutput(output, output_type, temporary_terms, blocks_temporary_terms_idxs);
          output << '=';
          rhs->writeOutput(output, output_type, temporary_terms, blocks_temporary_terms_idxs);
          output << ';' << endl;
          break;
        case BlockSimulationType::solveBackwardSimple:
        case BlockSimulationType::solveForwardSimple:
        case BlockSimulationType::solveBackwardComplete:
        case BlockSimulationType::solveForwardComplete:
        case BlockSimulationType::solveTwoBoundariesComplete:
        case BlockSimulationType::solveTwoBoundariesSimple:
          if (eq < block_recursive_size)
            goto evaluation;
          output << "  residual" << LEFT_ARRAY_SUBSCRIPT(output_type)
                 << eq-block_recursive_size+ARRAY_SUBSCRIPT_OFFSET(output_type)
                 << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=(";
          lhs->writeOutput(output, output_type, temporary_terms, blocks_temporary_terms_idxs);
          output << ")-(";
          rhs->writeOutput(output, output_type, temporary_terms, blocks_temporary_terms_idxs);
          output << ");" << endl;
          break;
        default:
          cerr << "Incorrect type for block " << blk+1 << endl;
          exit(EXIT_FAILURE);
        }
    }

  /* Write temporary terms for derivatives.
     This is done even for “evaluate” blocks, whose derivatives are not
     always computed at runtime; still those temporary terms may be needed by
     subsequent blocks. */
  write_eq_tt(blocks[blk].size);
}

template<ExprNodeOutputType output_type>
tuple<ostringstream, ostringstream, ostringstream, ostringstream,
      ostringstream, ostringstream, ostringstream>
ModelTree::writeParamsDerivativesFileHelper() const
{
  static_assert(!isCOutput(output_type), "C output is not implemented");

  constexpr bool sparse {isSparseModelOutput(output_type)};

  ostringstream tt_output; // Used for storing model temp vars and equations
  ostringstream rp_output; // 1st deriv. of residuals w.r.t. parameters
  ostringstream gp_output; // 1st deriv. of Jacobian w.r.t. parameters
  ostringstream rpp_output; // 2nd deriv of residuals w.r.t. parameters
  ostringstream gpp_output; // 2nd deriv of Jacobian w.r.t. parameters
  ostringstream hp_output; // 1st deriv. of Hessian w.r.t. parameters
  ostringstream g3p_output; // 1st deriv. of 3rd deriv. matrix w.r.t. parameters (only in dynamic case)

  temporary_terms_t temp_term_union;
  deriv_node_temp_terms_t tef_terms;

  for (const auto &[order, tts] : params_derivs_temporary_terms)
    writeTemporaryTerms<output_type>(tts, temp_term_union, params_derivs_temporary_terms_idxs,
                                     tt_output, tef_terms);

  for (const auto &[indices, d1] : params_derivatives.at({ 0, 1 }))
    {
      auto [eq, param] { vectorToTuple<2>(indices) };

      int param_col { getTypeSpecificIDByDerivID(param) + 1 };

      rp_output << "rp" << LEFT_ARRAY_SUBSCRIPT(output_type) << eq+1 << ", " << param_col
                << RIGHT_ARRAY_SUBSCRIPT(output_type) << " = ";
      d1->writeOutput(rp_output, output_type, temp_term_union, params_derivs_temporary_terms_idxs, tef_terms);
      rp_output << ";" << endl;
    }

  for (const auto &[indices, d2] : params_derivatives.at({ 1, 1 }))
    {
      auto [eq, var, param] { vectorToTuple<3>(indices) };

      int var_col { getJacobianCol(var, sparse) + 1 };
      int param_col { getTypeSpecificIDByDerivID(param) + 1 };

      gp_output << "gp" << LEFT_ARRAY_SUBSCRIPT(output_type) << eq+1 << ", " << var_col
                << ", " << param_col << RIGHT_ARRAY_SUBSCRIPT(output_type) << " = ";
      d2->writeOutput(gp_output, output_type, temp_term_union, params_derivs_temporary_terms_idxs, tef_terms);
      gp_output << ";" << endl;
    }

  for (int i {1};
       const auto &[indices, d2] : params_derivatives.at({ 0, 2 }))
    {
      auto [eq, param1, param2] { vectorToTuple<3>(indices) };

      int param1_col { getTypeSpecificIDByDerivID(param1) + 1 };
      int param2_col { getTypeSpecificIDByDerivID(param2) + 1 };

      rpp_output << "rpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",1"
                 << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << eq+1 << ";" << endl
                 << "rpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",2"
                 << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << param1_col << ";" << endl
                 << "rpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",3"
                 << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << param2_col << ";" << endl
                 << "rpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",4"
                 << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=";
      d2->writeOutput(rpp_output, output_type, temp_term_union, params_derivs_temporary_terms_idxs, tef_terms);
      rpp_output << ";" << endl;

      i++;

      if (param1 != param2)
        {
          // Treat symmetric elements
          rpp_output << "rpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",1"
                     << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << eq+1 << ";" << endl
                     << "rpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",2"
                     << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << param2_col << ";" << endl
                     << "rpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",3"
                     << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << param1_col << ";" << endl
                     << "rpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",4"
                     << RIGHT_ARRAY_SUBSCRIPT(output_type)
                     << "=rpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i-1 << ",4"
                     << RIGHT_ARRAY_SUBSCRIPT(output_type) << ";" << endl;
          i++;
        }
    }

  for (int i {1};
       const auto &[indices, d2] : params_derivatives.at({ 1, 2 }))
    {
      auto [eq, var, param1, param2] { vectorToTuple<4>(indices) };

      int var_col { getJacobianCol(var, sparse) + 1 };
      int param1_col { getTypeSpecificIDByDerivID(param1) + 1 };
      int param2_col { getTypeSpecificIDByDerivID(param2) + 1 };

      gpp_output << "gpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",1"
                 << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << eq+1 << ";" << endl
                 << "gpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",2"
                 << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << var_col << ";" << endl
                 << "gpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",3"
                 << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << param1_col << ";" << endl
                 << "gpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",4"
                 << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << param2_col << ";" << endl
                 << "gpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",5"
                 << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=";
      d2->writeOutput(gpp_output, output_type, temp_term_union, params_derivs_temporary_terms_idxs, tef_terms);
      gpp_output << ";" << endl;

      i++;

      if (param1 != param2)
        {
          // Treat symmetric elements
          gpp_output << "gpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",1"
                     << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << eq+1 << ";" << endl
                     << "gpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",2"
                     << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << var_col << ";" << endl
                     << "gpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",3"
                     << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << param2_col << ";" << endl
                     << "gpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",4"
                     << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << param1_col << ";" << endl
                     << "gpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",5"
                     << RIGHT_ARRAY_SUBSCRIPT(output_type)
                     << "=gpp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i-1 << ",5"
                     << RIGHT_ARRAY_SUBSCRIPT(output_type) << ";" << endl;
          i++;
        }
    }

  for (int i {1};
       const auto &[indices, d2] : params_derivatives.at({ 2, 1 }))
    {
      auto [eq, var1, var2, param] { vectorToTuple<4>(indices) };

      int var1_col { getJacobianCol(var1, sparse) + 1 };
      int var2_col { getJacobianCol(var2, sparse) + 1 };
      int param_col { getTypeSpecificIDByDerivID(param) + 1 };

      hp_output << "hp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",1"
                << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << eq+1 << ";" << endl
                << "hp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",2"
                << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << var1_col << ";" << endl
                << "hp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",3"
                << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << var2_col << ";" << endl
                << "hp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",4"
                << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << param_col << ";" << endl
                << "hp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",5"
                << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=";
      d2->writeOutput(hp_output, output_type, temp_term_union, params_derivs_temporary_terms_idxs, tef_terms);
      hp_output << ";" << endl;

      i++;

      if (var1 != var2)
        {
          // Treat symmetric elements
          hp_output << "hp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",1"
                    << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << eq+1 << ";" << endl
                    << "hp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",2"
                    << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << var2_col << ";" << endl
                    << "hp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",3"
                    << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << var1_col << ";" << endl
                    << "hp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",4"
                    << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << param_col << ";" << endl
                    << "hp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",5"
                    << RIGHT_ARRAY_SUBSCRIPT(output_type)
                    << "=hp" << LEFT_ARRAY_SUBSCRIPT(output_type) << i-1 << ",5"
                    << RIGHT_ARRAY_SUBSCRIPT(output_type) << ";" << endl;
          i++;
        }
    }

  if constexpr(output_type == ExprNodeOutputType::matlabDynamicModel
               || output_type == ExprNodeOutputType::juliaDynamicModel)
    for (int i {1};
         const auto &[indices, d2] : params_derivatives.at({ 3, 1 }))
      {
        auto [eq, var1, var2, var3, param] { vectorToTuple<5>(indices) };

        int var1_col { getJacobianCol(var1, sparse) + 1 };
        int var2_col { getJacobianCol(var2, sparse) + 1 };
        int var3_col { getJacobianCol(var3, sparse) + 1 };
        int param_col { getTypeSpecificIDByDerivID(param) + 1 };

        g3p_output << "g3p" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",1"
                   << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << eq+1 << ";" << endl
                   << "g3p" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",2"
                   << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << var1_col << ";" << endl
                   << "g3p" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",3"
                   << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << var2_col << ";" << endl
                   << "g3p" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",4"
                   << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << var3_col << ";" << endl
                   << "g3p" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",5"
                   << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=" << param_col << ";" << endl
                   << "g3p" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << ",6"
                   << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=";
        d2->writeOutput(g3p_output, output_type, temp_term_union, params_derivs_temporary_terms_idxs, tef_terms);
        g3p_output << ";" << endl;

        i++;
      }

  if constexpr(isMatlabOutput(output_type))
    {
      // Check that we don't have more than 32 nested parenthesis because MATLAB does not support this. See Issue #1201
      map<string, string> tmp_paren_vars;
      bool message_printed {false};
      fixNestedParenthesis(tt_output, tmp_paren_vars, message_printed);
      fixNestedParenthesis(rp_output, tmp_paren_vars, message_printed);
      fixNestedParenthesis(gp_output, tmp_paren_vars, message_printed);
      fixNestedParenthesis(rpp_output, tmp_paren_vars, message_printed);
      fixNestedParenthesis(gpp_output, tmp_paren_vars, message_printed);
      fixNestedParenthesis(hp_output, tmp_paren_vars, message_printed);
      fixNestedParenthesis(g3p_output, tmp_paren_vars, message_printed);
    }

  return { move(tt_output), move(rp_output), move(gp_output),
    move(rpp_output), move(gpp_output), move(hp_output), move(g3p_output) };
}

template<ExprNodeBytecodeOutputType output_type>
void
ModelTree::writeBytecodeTemporaryTerms(const temporary_terms_t &tt,
                                       temporary_terms_t &temporary_terms_union,
                                       BytecodeWriter &code_file,
                                       deriv_node_temp_terms_t &tef_terms) const
{
  for (auto it : tt)
    {
      if (dynamic_cast<AbstractExternalFunctionNode *>(it))
        it->writeBytecodeExternalFunctionOutput(code_file, output_type, temporary_terms_union, temporary_terms_idxs, tef_terms);

      int idx {temporary_terms_idxs.at(it)};
      code_file << FNUMEXPR_{ExpressionType::TemporaryTerm, idx};
      it->writeBytecodeOutput(code_file, output_type, temporary_terms_union, temporary_terms_idxs, tef_terms);

      static_assert(output_type == ExprNodeBytecodeOutputType::dynamicModel
                    || output_type == ExprNodeBytecodeOutputType::staticModel);
      if constexpr(output_type == ExprNodeBytecodeOutputType::dynamicModel)
        code_file << FSTPT_{idx};
      else
        code_file << FSTPST_{idx};

      temporary_terms_union.insert(it);
    }
}

template<ExprNodeBytecodeOutputType output_type>
void
ModelTree::writeBytecodeModelEquations(BytecodeWriter &code_file, const temporary_terms_t &temporary_terms, const deriv_node_temp_terms_t &tef_terms) const
{
  for (int eq {0}; eq < static_cast<int>(equations.size()); eq++)
    {
      BinaryOpNode *eq_node {equations[eq]};
      expr_t lhs {eq_node->arg1}, rhs {eq_node->arg2};
      code_file << FNUMEXPR_{ExpressionType::ModelEquation, eq};
      // Test if the right hand side of the equation is empty.
      double vrhs {1.0};
      try
        {
          vrhs = rhs->eval({});
        }
      catch (ExprNode::EvalException &e)
        {
        }

      if (vrhs != 0) // The right hand side of the equation is not empty ⇒ residual=lhs-rhs
        {
          lhs->writeBytecodeOutput(code_file, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
          rhs->writeBytecodeOutput(code_file, output_type, temporary_terms, temporary_terms_idxs, tef_terms);

          code_file << FBINARY_{BinaryOpcode::minus} << FSTPR_{eq};
        }
      else // The right hand side of the equation is empty ⇒ residual=lhs
        {
          lhs->writeBytecodeOutput(code_file, output_type, temporary_terms, temporary_terms_idxs, tef_terms);
          code_file << FSTPR_{eq};
        }
    }
}

template<bool dynamic>
void
ModelTree::writeBytecodeHelper(BytecodeWriter &code_file) const
{
  constexpr ExprNodeBytecodeOutputType output_type { dynamic ? ExprNodeBytecodeOutputType::dynamicModel : ExprNodeBytecodeOutputType::staticModel };

  temporary_terms_t temporary_terms_union;
  deriv_node_temp_terms_t tef_terms;

  writeBytecodeTemporaryTerms<output_type>(temporary_terms_derivatives[0], temporary_terms_union, code_file, tef_terms);
  writeBytecodeModelEquations<output_type>(code_file, temporary_terms_union, tef_terms);

  code_file << FENDEQU_{};

  // Temporary terms for the Jacobian
  writeBytecodeTemporaryTerms<output_type>(temporary_terms_derivatives[1], temporary_terms_union, code_file, tef_terms);

  // Get the current code_file position and jump if “evaluate” mode
  int pos_jmpifeval {code_file.getInstructionCounter()};
  code_file << FJMPIFEVAL_{0}; // Use 0 as jump offset for the time being

  // The Jacobian in “simulate” mode
  vector<vector<tuple<int, int, int>>> my_derivatives(symbol_table.endo_nbr());;
  int count_u {symbol_table.endo_nbr()};
  for (const auto &[indices, d1] : derivatives[1])
    {
      auto [eq, deriv_id] {vectorToTuple<2>(indices)};
      if (getTypeByDerivID(deriv_id) == SymbolType::endogenous)
        {
          int tsid {getTypeSpecificIDByDerivID(deriv_id)};
          int lag {getLagByDerivID(deriv_id)};
          if constexpr(dynamic)
            code_file << FNUMEXPR_{ExpressionType::FirstEndoDerivative, eq, tsid, lag};
          else
            code_file << FNUMEXPR_{ExpressionType::FirstEndoDerivative, eq, tsid};
          if (!my_derivatives[eq].size())
            my_derivatives[eq].clear();
          my_derivatives[eq].emplace_back(tsid, lag, count_u);
          d1->writeBytecodeOutput(code_file, output_type, temporary_terms_union, temporary_terms_idxs, tef_terms);
          if constexpr(dynamic)
            code_file << FSTPU_{count_u};
          else
            code_file << FSTPSU_{count_u};
          count_u++;
        }
    }
  for (int i {0}; i < symbol_table.endo_nbr(); i++)
    {
      code_file << FLDR_{i};
      if (my_derivatives[i].size())
        {
          for (bool first_term {true};
               const auto &[tsid, lag, uidx] : my_derivatives[i])
            {
              if constexpr(dynamic)
                code_file << FLDU_{uidx} << FLDV_{SymbolType::endogenous, tsid, lag};
              else
                code_file << FLDSU_{uidx} << FLDSV_{SymbolType::endogenous, tsid};
              code_file << FBINARY_{BinaryOpcode::times};
              if (!exchange(first_term, false))
                code_file << FBINARY_{BinaryOpcode::plus};
            }
          code_file << FBINARY_{BinaryOpcode::minus};
        }
      if constexpr(dynamic)
        code_file << FSTPU_{i};
      else
        code_file << FSTPSU_{i};
    }

  // Jump unconditionally after the block
  int pos_jmp {code_file.getInstructionCounter()};
  code_file << FJMP_{0}; // Use 0 as jump offset for the time being
  // Update jump offset for previous JMPIFEVAL
  code_file.overwriteInstruction(pos_jmpifeval, FJMPIFEVAL_{pos_jmp-pos_jmpifeval});

  // The Jacobian in “evaluate” mode
  for (const auto &[indices, d1] : derivatives[1])
    {
      auto [eq, deriv_id] {vectorToTuple<2>(indices)};
      int tsid {getTypeSpecificIDByDerivID(deriv_id)};
      int lag {getLagByDerivID(deriv_id)};
      SymbolType type {getTypeByDerivID(deriv_id)};

      if constexpr(dynamic)
        {
          ExpressionType expr_type;
          switch (type)
            {
            case SymbolType::endogenous:
              expr_type = ExpressionType::FirstEndoDerivative;
              break;
            case SymbolType::exogenous:
              expr_type = ExpressionType::FirstExoDerivative;
              break;
            case SymbolType::exogenousDet:
              expr_type = ExpressionType::FirstExodetDerivative;
              break;
            default:
              assert(false);
              break;
            }
          code_file << FNUMEXPR_{expr_type, eq, tsid, lag};
        }
      else
        {
          assert(type == SymbolType::endogenous);
          code_file << FNUMEXPR_{ExpressionType::FirstEndoDerivative, eq, tsid};
        }

      d1->writeBytecodeOutput(code_file, output_type, temporary_terms_union, temporary_terms_idxs, tef_terms);
      if constexpr(dynamic)
         {
           // Bytecode MEX uses a separate matrix for exogenous and exodet Jacobians
           int jacob_col { type == SymbolType::endogenous ? getJacobianCol(deriv_id, false) : tsid };
           code_file << FSTPG3_{eq, tsid, lag, jacob_col};
         }
      else
        code_file << FSTPG2_{eq, tsid};
    }

  // Update jump offset for previous JMP
  int pos_end_block {code_file.getInstructionCounter()};
  code_file.overwriteInstruction(pos_jmp, FJMP_{pos_end_block-pos_jmp-1});

  code_file << FENDBLOCK_{} << FEND_{};
}

template<bool dynamic>
void
ModelTree::writeBlockBytecodeHelper(BytecodeWriter &code_file, int block) const
{
  constexpr ExprNodeBytecodeOutputType output_type
    { dynamic ? ExprNodeBytecodeOutputType::dynamicModel : ExprNodeBytecodeOutputType::staticModel };
  constexpr ExprNodeBytecodeOutputType assignment_lhs_output_type
    { dynamic ? ExprNodeBytecodeOutputType::dynamicAssignmentLHS : ExprNodeBytecodeOutputType::staticAssignmentLHS };

  const BlockSimulationType simulation_type {blocks[block].simulation_type};
  const int block_size {blocks[block].size};
  const int block_mfs {blocks[block].mfs_size};
  const int block_recursive {blocks[block].getRecursiveSize()};

  temporary_terms_t temporary_terms_union;
  deriv_node_temp_terms_t tef_terms;

  auto write_eq_tt = [&](int eq)
  {
    for (auto it : blocks_temporary_terms[block][eq])
      {
        if (dynamic_cast<AbstractExternalFunctionNode *>(it))
          it->writeBytecodeExternalFunctionOutput(code_file, output_type, temporary_terms_union, blocks_temporary_terms_idxs, tef_terms);

        code_file << FNUMEXPR_{ExpressionType::TemporaryTerm, blocks_temporary_terms_idxs.at(it)};
        it->writeBytecodeOutput(code_file, output_type, temporary_terms_union, blocks_temporary_terms_idxs, tef_terms);
        if constexpr(dynamic)
          code_file << FSTPT_{blocks_temporary_terms_idxs.at(it)};
        else
          code_file << FSTPST_{blocks_temporary_terms_idxs.at(it)};
        temporary_terms_union.insert(it);
      }
  };

  // The equations
  for (int i {0}; i < block_size; i++)
    {
      write_eq_tt(i);

      switch (simulation_type)
        {
        evaluation:
        case BlockSimulationType::evaluateBackward:
        case BlockSimulationType::evaluateForward:
          code_file << FNUMEXPR_{ExpressionType::ModelEquation, getBlockEquationID(block, i)};
          if (EquationType equ_type {getBlockEquationType(block, i)};
              equ_type == EquationType::evaluate)
            {
              BinaryOpNode *eq_node {getBlockEquationExpr(block, i)};
              expr_t lhs {eq_node->arg1};
              expr_t rhs {eq_node->arg2};
              rhs->writeBytecodeOutput(code_file, output_type, temporary_terms_union, blocks_temporary_terms_idxs, tef_terms);
              lhs->writeBytecodeOutput(code_file, assignment_lhs_output_type, temporary_terms_union, blocks_temporary_terms_idxs, tef_terms);
            }
          else if (equ_type == EquationType::evaluateRenormalized)
            {
              BinaryOpNode *eq_node {getBlockEquationRenormalizedExpr(block, i)};
              expr_t lhs {eq_node->arg1};
              expr_t rhs {eq_node->arg2};
              rhs->writeBytecodeOutput(code_file, output_type, temporary_terms_union, blocks_temporary_terms_idxs, tef_terms);
              lhs->writeBytecodeOutput(code_file, assignment_lhs_output_type, temporary_terms_union, blocks_temporary_terms_idxs, tef_terms);
            }
          break;
        case BlockSimulationType::solveBackwardComplete:
        case BlockSimulationType::solveForwardComplete:
        case BlockSimulationType::solveTwoBoundariesComplete:
        case BlockSimulationType::solveTwoBoundariesSimple:
          if (i < block_recursive)
            goto evaluation;
          [[fallthrough]];
        default:
          code_file << FNUMEXPR_{ExpressionType::ModelEquation, getBlockEquationID(block, i)};
          BinaryOpNode *eq_node {getBlockEquationExpr(block, i)};
          expr_t lhs {eq_node->arg1};
          expr_t rhs {eq_node->arg2};
          lhs->writeBytecodeOutput(code_file, output_type, temporary_terms_union, blocks_temporary_terms_idxs, tef_terms);
          rhs->writeBytecodeOutput(code_file, output_type, temporary_terms_union, blocks_temporary_terms_idxs, tef_terms);
          code_file << FBINARY_{BinaryOpcode::minus} << FSTPR_{i - block_recursive};
        }
    }

  /* Write temporary terms for derivatives. This is done before FENDEQU,
     because residuals of a subsequent block may depend on temporary terms for
     the derivatives of the present block.

     Also note that in the case of “evaluate” blocks, derivatives are not
     computed in the “evaluate” mode; still their temporary terms must be
     computed even in that mode, because for the same reason as above they may
     be needed in subsequent blocks. */
  write_eq_tt(blocks[block].size);

  code_file << FENDEQU_{};

  // Get the current code_file position and jump if evaluating
  int pos_jmpifeval {code_file.getInstructionCounter()};
  code_file << FJMPIFEVAL_{0}; // Use 0 as jump offset for the time being

  /* Write the derivatives for the “simulate” mode (not needed if the block
     is of type “evaluate backward/forward”) */
  if (simulation_type != BlockSimulationType::evaluateBackward
      && simulation_type != BlockSimulationType::evaluateForward)
    {
      switch (simulation_type)
        {
        case BlockSimulationType::solveBackwardSimple:
        case BlockSimulationType::solveForwardSimple:
          {
            int eqr {getBlockEquationID(block, 0)};
            int varr {getBlockVariableID(block, 0)};
            code_file << FNUMEXPR_{ExpressionType::FirstEndoDerivative, eqr, varr, 0};
            // Get contemporaneous derivative of the single variable in the block
            if (auto it { blocks_derivatives[block].find({ 0, 0, 0 }) };
                it != blocks_derivatives[block].end())
              it->second->writeBytecodeOutput(code_file, output_type, temporary_terms_union, blocks_temporary_terms_idxs, tef_terms);
            else
              code_file << FLDZ_{};
            code_file << FSTPG_{0};
          }
          break;

        case BlockSimulationType::solveBackwardComplete:
        case BlockSimulationType::solveForwardComplete:
        case BlockSimulationType::solveTwoBoundariesComplete:
        case BlockSimulationType::solveTwoBoundariesSimple:
          {
            // For each equation, stores a list of tuples (index_u, var, lag)
            vector<vector<tuple<int, int, int>>> Uf(symbol_table.endo_nbr());

            for (int count_u {block_mfs};
                 const auto &[indices, ignore] : blocks_derivatives[block])
              {
                const auto &[eq, var, lag] {indices};
                int eqr {getBlockEquationID(block, eq)};
                int varr {getBlockVariableID(block, var)};
                if (eq >= block_recursive && var >= block_recursive)
                  {
                    if constexpr(dynamic)
                      if (lag != 0
                          && (simulation_type == BlockSimulationType::solveForwardComplete
                              || simulation_type == BlockSimulationType::solveBackwardComplete))
                        continue;
                    code_file << FNUMEXPR_{ExpressionType::FirstEndoDerivative, eqr, varr, lag};
                    if (auto it { blocks_derivatives[block].find({ eq, var, lag }) };
                        it != blocks_derivatives[block].end())
                      it->second->writeBytecodeOutput(code_file, output_type, temporary_terms_union, blocks_temporary_terms_idxs, tef_terms);
                    else
                      code_file << FLDZ_{};
                    if constexpr(dynamic)
                      code_file << FSTPU_{count_u};
                    else
                      code_file << FSTPSU_{count_u};
                    Uf[eqr].emplace_back(count_u, varr, lag);
                    count_u++;
                  }
              }
            for (int i {0}; i < block_size; i++)
              if (i >= block_recursive)
                {
                  code_file << FLDR_{i-block_recursive} << FLDZ_{};

                  int eqr {getBlockEquationID(block, i)};
                  for (const auto &[index_u, var, lag] : Uf[eqr])
                    {
                      if constexpr(dynamic)
                        code_file << FLDU_{index_u}
                                  << FLDV_{SymbolType::endogenous, var, lag};
                      else
                        code_file << FLDSU_{index_u}
                                  << FLDSV_{SymbolType::endogenous, var};
                      code_file << FBINARY_{BinaryOpcode::times}
                                << FBINARY_{BinaryOpcode::plus};
                    }
                  code_file << FBINARY_{BinaryOpcode::minus};
                  if constexpr(dynamic)
                    code_file << FSTPU_{i - block_recursive};
                  else
                    code_file << FSTPSU_{i - block_recursive};
                }
          }
          break;
        default:
          break;
        }
    }

  // Jump unconditionally after the block
  int pos_jmp {code_file.getInstructionCounter()};
  code_file << FJMP_{0}; // Use 0 as jump offset for the time being
  // Update jump offset for previous JMPIFEVAL
  code_file.overwriteInstruction(pos_jmpifeval, FJMPIFEVAL_{pos_jmp-pos_jmpifeval});

  // Write the derivatives for the “evaluate” mode
  for (const auto &[indices, d] : blocks_derivatives[block])
    {
      const auto &[eq, var, lag] {indices};
      int eqr {getBlockEquationID(block, eq)};
      int varr {getBlockVariableID(block, var)};
      code_file << FNUMEXPR_{ExpressionType::FirstEndoDerivative, eqr, varr, lag};
      d->writeBytecodeOutput(code_file, output_type, temporary_terms_union, blocks_temporary_terms_idxs, tef_terms);
      if constexpr(dynamic)
        code_file << FSTPG3_{eq, var, lag, getBlockJacobianEndoCol(block, var, lag)};
      else
        code_file << FSTPG2_{eq, getBlockJacobianEndoCol(block, var, lag)};
    }

  /* Write derivatives w.r.t. exo, exo det and other endogenous, but only in
     dynamic mode */
  writeBlockBytecodeAdditionalDerivatives(code_file, block, temporary_terms_union, tef_terms);

  // Update jump offset for previous JMP
  int pos_end_block {code_file.getInstructionCounter()};
  code_file.overwriteInstruction(pos_jmp, FJMP_{pos_end_block-pos_jmp-1});

  code_file << FENDBLOCK_{};
}

template<bool dynamic>
pair<ostringstream, vector<ostringstream>>
ModelTree::writeJsonComputingPassOutputHelper(bool writeDetails) const
{
  ostringstream mlv_output; // Used for storing model local vars
  vector<ostringstream> d_output(derivatives.size()); // Derivatives output (at all orders, including 0=residual)

  deriv_node_temp_terms_t tef_terms;
  temporary_terms_t temp_term_union;

  writeJsonModelLocalVariables(mlv_output, true, tef_terms);

  writeJsonTemporaryTerms(temporary_terms_derivatives[0], temp_term_union, d_output[0], tef_terms, "");
  d_output[0] << ", ";
  writeJsonModelEquations(d_output[0], true);

  int ncols { getJacobianColsNbr(false) };
  for (size_t i {1}; i < derivatives.size(); i++)
    {
      string matrix_name { i == 1 ? "jacobian" : i == 2 ? "hessian" : i == 3 ? "third_derivative" : to_string(i) + "th_derivative"};
      writeJsonTemporaryTerms(temporary_terms_derivatives[i], temp_term_union, d_output[i], tef_terms, matrix_name);
      temp_term_union.insert(temporary_terms_derivatives[i].begin(), temporary_terms_derivatives[i].end());
      d_output[i] << R"(, ")" << matrix_name  << R"(": {)"
                  << R"(  "nrows": )" << equations.size()
                  << R"(, "ncols": )" << ncols
                  << R"(, "entries": [)";

      for (bool printed_something {false};
           const auto &[vidx, d] : derivatives[i])
        {
          if (exchange(printed_something, true))
            d_output[i] << ", ";

          int eq { vidx[0] };

          int col_idx {0};
          for (size_t j {1}; j < vidx.size(); j++)
            {
              col_idx *= getJacobianColsNbr(false);
              col_idx += getJacobianCol(vidx[j], false);
            }

          if (writeDetails)
            d_output[i] << R"({"eq": )" << eq + 1;
          else
            d_output[i] << R"({"row": )" << eq + 1;

          d_output[i] << R"(, "col": )" << (i > 1 ? "[" : "") << col_idx + 1;

          if (i == 2 && vidx[1] != vidx[2]) // Symmetric elements in hessian
            {
              int col_idx_sym { getJacobianCol(vidx[2], false) * getJacobianColsNbr(false) + getJacobianCol(vidx[1], false)};
              d_output[i] << ", " << col_idx_sym + 1;
            }
          if (i > 1)
            d_output[i] << "]";

          if (writeDetails)
            for (size_t j = 1; j < vidx.size(); j++)
              {
                d_output[i] << R"(, "var)" << (i > 1 ? to_string(j) : "") << R"(": ")" << getNameByDerivID(vidx[j]) << R"(")";
                if constexpr(dynamic)
                   d_output[i] << R"(, "shift)" << (i > 1 ? to_string(j) : "") << R"(": )" << getLagByDerivID(vidx[j]);
              }

          d_output[i] << R"(, "val": ")";
          d->writeJsonOutput(d_output[i], temp_term_union, tef_terms);
          d_output[i] << R"("})" << endl;
        }
      d_output[i] << "]}";

      ncols *= getJacobianColsNbr(false);
    }

  return { move(mlv_output), move(d_output) };
}

template<bool dynamic>
tuple<ostringstream, ostringstream, ostringstream, ostringstream,
      ostringstream, ostringstream, ostringstream, ostringstream>
ModelTree::writeJsonParamsDerivativesHelper(bool writeDetails) const
{
  ostringstream mlv_output; // Used for storing model local vars
  ostringstream tt_output; // Used for storing model temp vars and equations
  ostringstream rp_output; // 1st deriv. of residuals w.r.t. parameters
  ostringstream gp_output; // 1st deriv. of Jacobian w.r.t. parameters
  ostringstream rpp_output; // 2nd deriv of residuals w.r.t. parameters
  ostringstream gpp_output; // 2nd deriv of Jacobian w.r.t. parameters
  ostringstream hp_output; // 1st deriv. of Hessian w.r.t. parameters
  ostringstream g3p_output; // 1st deriv. of 3rd deriv. matrix w.r.t. parameters

  deriv_node_temp_terms_t tef_terms;
  writeJsonModelLocalVariables(mlv_output, true, tef_terms);

  temporary_terms_t temp_term_union;
  for (const auto &[order, tts] : params_derivs_temporary_terms)
    writeJsonTemporaryTerms(tts, temp_term_union, tt_output, tef_terms, "all");

  rp_output << R"("deriv_wrt_params": {)"
            << R"(  "neqs": )" << equations.size()
            << R"(, "nparamcols": )" << symbol_table.param_nbr()
            << R"(, "entries": [)";
  for (bool printed_something {false};
       const auto &[vidx, d] : params_derivatives.at({ 0, 1 }))
    {
      if (exchange(printed_something, true))
        rp_output << ", ";

      auto [eq, param] { vectorToTuple<2>(vidx) };

      int param_col { getTypeSpecificIDByDerivID(param) + 1 };

      if (writeDetails)
        rp_output << R"({"eq": )" << eq + 1;
      else
        rp_output << R"({"row": )" << eq + 1;

      rp_output << R"(, "param_col": )" << param_col;

      if (writeDetails)
        rp_output << R"(, "param": ")" << getNameByDerivID(param) << R"(")";

      rp_output << R"(, "val": ")";
      d->writeJsonOutput(rp_output, temp_term_union, tef_terms);
      rp_output << R"("})" << endl;
    }
  rp_output << "]}";

  gp_output << R"("deriv_jacobian_wrt_params": {)"
            << R"(  "neqs": )" << equations.size()
            << R"(, "nvarcols": )" << getJacobianColsNbr(false)
            << R"(, "nparamcols": )" << symbol_table.param_nbr()
            << R"(, "entries": [)";
  for (bool printed_something {false};
       const auto &[vidx, d] : params_derivatives.at({ 1, 1 }))
    {
      if (exchange(printed_something, true))
        gp_output << ", ";

      auto [eq, var, param] { vectorToTuple<3>(vidx) };

      int var_col { getJacobianCol(var, false) + 1 };
      int param_col { getTypeSpecificIDByDerivID(param) + 1 };

      if (writeDetails)
        gp_output << R"({"eq": )" << eq + 1;
      else
        gp_output << R"({"row": )" << eq + 1;

      gp_output << R"(, "var_col": )" << var_col
                << R"(, "param_col": )" << param_col;

      if (writeDetails)
        {
          gp_output << R"(, "var": ")" << getNameByDerivID(var) << R"(")";
          if constexpr(dynamic)
            gp_output << R"(, "lag": )" << getLagByDerivID(var);
          gp_output << R"(, "param": ")" << getNameByDerivID(param) << R"(")";
        }

      gp_output << R"(, "val": ")";
      d->writeJsonOutput(gp_output, temp_term_union, tef_terms);
      gp_output << R"("})" << endl;
    }
  gp_output << "]}";

  rpp_output << R"("second_deriv_residuals_wrt_params": {)"
             << R"(  "nrows": )" << equations.size()
             << R"(, "nparam1cols": )" << symbol_table.param_nbr()
             << R"(, "nparam2cols": )" << symbol_table.param_nbr()
             << R"(, "entries": [)";
  for (bool printed_something {false};
       const auto &[vidx, d] : params_derivatives.at({ 0, 2 }))
    {
      if (exchange(printed_something, true))
        rpp_output << ", ";

      auto [eq, param1, param2] { vectorToTuple<3>(vidx) };

      int param1_col { getTypeSpecificIDByDerivID(param1) + 1 };
      int param2_col { getTypeSpecificIDByDerivID(param2) + 1 };

      if (writeDetails)
        rpp_output << R"({"eq": )" << eq + 1;
      else
        rpp_output << R"({"row": )" << eq + 1;
      rpp_output << R"(, "param1_col": )" << param1_col
                 << R"(, "param2_col": )" << param2_col;

      if (writeDetails)
        rpp_output << R"(, "param1": ")" << getNameByDerivID(param1) << R"(")"
                   << R"(, "param2": ")" << getNameByDerivID(param2) << R"(")";

      rpp_output << R"(, "val": ")";
      d->writeJsonOutput(rpp_output, temp_term_union, tef_terms);
      rpp_output << R"("})" << endl;
    }
  rpp_output << "]}";

  gpp_output << R"("second_deriv_jacobian_wrt_params": {)"
             << R"(  "neqs": )" << equations.size()
             << R"(, "nvarcols": )" << getJacobianColsNbr(false)
             << R"(, "nparam1cols": )" << symbol_table.param_nbr()
             << R"(, "nparam2cols": )" << symbol_table.param_nbr()
             << R"(, "entries": [)";
  for (bool printed_something {false};
       const auto &[vidx, d] : params_derivatives.at({ 1, 2 }))
    {
      if (exchange(printed_something, true))
        gpp_output << ", ";

      auto [eq, var, param1, param2] { vectorToTuple<4>(vidx) };

      int var_col { getJacobianCol(var, false) + 1 };
      int param1_col { getTypeSpecificIDByDerivID(param1) + 1 };
      int param2_col { getTypeSpecificIDByDerivID(param2) + 1 };

      if (writeDetails)
        gpp_output << R"({"eq": )" << eq + 1;
      else
        gpp_output << R"({"row": )" << eq + 1;

      gpp_output << R"(, "var_col": )" << var_col
                 << R"(, "param1_col": )" << param1_col
                 << R"(, "param2_col": )" << param2_col;

      if (writeDetails)
        {
          gpp_output << R"(, "var": ")" << getNameByDerivID(var) << R"(")";
          if constexpr(dynamic)
            gpp_output << R"(, "lag": )" << getLagByDerivID(var);
          gpp_output << R"(, "param1": ")" << getNameByDerivID(param1) << R"(")"
                     << R"(, "param2": ")" << getNameByDerivID(param2) << R"(")";
        }

      gpp_output << R"(, "val": ")";
      d->writeJsonOutput(gpp_output, temp_term_union, tef_terms);
      gpp_output << R"("})" << endl;
    }
  gpp_output << "]}" << endl;

  hp_output << R"("derivative_hessian_wrt_params": {)"
            << R"(  "neqs": )" << equations.size()
            << R"(, "nvar1cols": )" << getJacobianColsNbr(false)
            << R"(, "nvar2cols": )" << getJacobianColsNbr(false)
            << R"(, "nparamcols": )" << symbol_table.param_nbr()
            << R"(, "entries": [)";
  for (bool printed_something {false};
       const auto &[vidx, d] : params_derivatives.at({ 2, 1 }))
    {
      if (exchange(printed_something, true))
        hp_output << ", ";

      auto [eq, var1, var2, param] { vectorToTuple<4>(vidx) };

      int var1_col { getJacobianCol(var1, false) + 1 };
      int var2_col { getJacobianCol(var2, false) + 1 };
      int param_col { getTypeSpecificIDByDerivID(param) + 1 };

      if (writeDetails)
        hp_output << R"({"eq": )" << eq + 1;
      else
        hp_output << R"({"row": )" << eq + 1;

      hp_output << R"(, "var1_col": )" << var1_col
                << R"(, "var2_col": )" << var2_col
                << R"(, "param_col": )" << param_col;

      if (writeDetails)
        {
          hp_output << R"(, "var1": ")" << getNameByDerivID(var1) << R"(")";
          if constexpr(dynamic)
            hp_output << R"(, "lag1": )" << getLagByDerivID(var1);
          hp_output << R"(, "var2": ")" << getNameByDerivID(var2) << R"(")";
          if constexpr(dynamic)
            hp_output << R"(, "lag2": )" << getLagByDerivID(var2);
          hp_output << R"(, "param": ")" << getNameByDerivID(param) << R"(")";
        }

      hp_output << R"(, "val": ")";
      d->writeJsonOutput(hp_output, temp_term_union, tef_terms);
      hp_output << R"("})" << endl;
    }
  hp_output << "]}" << endl;

  if constexpr(dynamic)
    {
      g3p_output << R"("derivative_g3_wrt_params": {)"
                 << R"(  "neqs": )" << equations.size()
                 << R"(, "nvar1cols": )" << getJacobianColsNbr(false)
                 << R"(, "nvar2cols": )" << getJacobianColsNbr(false)
                 << R"(, "nvar3cols": )" << getJacobianColsNbr(false)
                 << R"(, "nparamcols": )" << symbol_table.param_nbr()
                 << R"(, "entries": [)";
      for (bool printed_something {false};
           const auto &[vidx, d] : params_derivatives.at({ 3, 1 }))
        {
          if (exchange(printed_something, true))
            g3p_output << ", ";

          auto [eq, var1, var2, var3, param] { vectorToTuple<5>(vidx) };

          int var1_col { getJacobianCol(var1, false) + 1 };
          int var2_col { getJacobianCol(var2, false) + 1 };
          int var3_col { getJacobianCol(var3, false) + 1 };
          int param_col { getTypeSpecificIDByDerivID(param) + 1 };

          if (writeDetails)
            g3p_output << R"({"eq": )" << eq + 1;
          else
            g3p_output << R"({"row": )" << eq + 1;

          g3p_output << R"(, "var1_col": )" << var1_col + 1
                     << R"(, "var2_col": )" << var2_col + 1
                     << R"(, "var3_col": )" << var3_col + 1
                     << R"(, "param_col": )" << param_col + 1;

          if (writeDetails)
            g3p_output << R"(, "var1": ")" << getNameByDerivID(var1) << R"(")"
                       << R"(, "lag1": )" << getLagByDerivID(var1)
                       << R"(, "var2": ")" << getNameByDerivID(var2) << R"(")"
                       << R"(, "lag2": )" << getLagByDerivID(var2)
                       << R"(, "var3": ")" << getNameByDerivID(var3) << R"(")"
                       << R"(, "lag3": )" << getLagByDerivID(var3)
                       << R"(, "param": ")" << getNameByDerivID(param) << R"(")";

          g3p_output << R"(, "val": ")";
          d->writeJsonOutput(g3p_output, temp_term_union, tef_terms);
          g3p_output << R"("})" << endl;
        }
      g3p_output << "]}" << endl;
    }

  return { move(mlv_output), move(tt_output), move(rp_output), move(gp_output),
    move(rpp_output), move(gpp_output), move(hp_output), move(g3p_output) };
}

template<bool dynamic>
void
ModelTree::writeDriverSparseIndicesHelper(ostream &output) const
{
  // TODO: when C++20 support is complete, mark this constexpr
  const string model_name {dynamic ? "dynamic" : "static"};

  // Write indices for the sparse Jacobian (both naive and CSC storage)
  output << "M_." << model_name << "_g1_sparse_rowval = int32([";
  for (const auto &[indices, d1] : jacobian_sparse_column_major_order)
    output << indices.first+1 << ' ';
  output << "]);" << endl
         << "M_." << model_name << "_g1_sparse_colval = int32([";
  for (const auto &[indices, d1] : jacobian_sparse_column_major_order)
    output << indices.second+1 << ' ';
  output << "]);" << endl
         << "M_." << model_name << "_g1_sparse_colptr = int32([";
  for (int it : jacobian_sparse_colptr)
    output << it+1 << ' ';
  output << "]);" << endl;

  // Write indices for the sparse higher-order derivatives
  for (int i {2}; i < computed_derivs_order; i++)
    {
      output << "M_." << model_name << "_g" << i << "_sparse_indices = int32([";
      for (const auto &[vidx, d] : derivatives[i])
        {
          for (int it : vidx)
            output << it+1 << ' ';
          output << ';' << endl;
        }
      output << "]);" << endl;
    }
}

template<bool dynamic>
void
ModelTree::writeJsonSparseIndicesHelper(ostream &output) const
{
  // TODO: when C++20 support is complete, mark this constexpr
  const string model_name {dynamic ? "dynamic" : "static"};

  // Write indices for the sparse Jacobian (both naive and CSC storage)
  output << '"' << model_name << R"(_g1_sparse_rowval": [)";
  for (bool printed_something {false};
       const auto &[indices, d1] : jacobian_sparse_column_major_order)
    {
      if (exchange(printed_something, true))
        output << ", ";
      output << indices.first+1;
    }
  output << "], " << endl
         << '"' << model_name << R"(_g1_sparse_colval": [)";
  for (bool printed_something {false};
       const auto &[indices, d1] : jacobian_sparse_column_major_order)
    {
      if (exchange(printed_something, true))
        output << ", ";
      output << indices.second+1;
    }
  output << "], " << endl
         << '"' << model_name << R"(_g1_sparse_colptr": [)";
  for (bool printed_something {false};
       int it : jacobian_sparse_colptr)
    {
      if (exchange(printed_something, true))
        output << ", ";
      output << it+1;
    }
  output << ']' << endl;

  // Write indices for the sparse higher-order derivatives
  for (int i {2}; i < computed_derivs_order; i++)
    {
      output << R"(, ")" << model_name << "_g" << i << R"(_sparse_indices": [)";
      for (bool printed_something {false};
           const auto &[vidx, d] : derivatives[i])
        {
          if (exchange(printed_something, true))
            output << ", ";
          output << '[';
          for (bool printed_something2 {false};
               int it : vidx)
            {
              if (exchange(printed_something2, true))
                output << ", ";
              output << it+1;
            }
          output << ']' << endl;
        }
      output << ']' << endl;
    }
}

template<bool dynamic>
void
ModelTree::writeSparseModelJuliaFiles(const string &basename) const
{
  auto [d_sparse_output, tt_sparse_output] = writeModelFileHelper<dynamic ? ExprNodeOutputType::juliaSparseDynamicModel : ExprNodeOutputType::juliaSparseStaticModel>();

  filesystem::path julia_dir {filesystem::path{basename} / "model" / "julia"};
  // TODO: when C++20 support is complete, mark the following strings constexpr
  const string prefix { dynamic ? "SparseDynamic" : "SparseStatic" };
  const string ss_argin { dynamic ? ", steady_state::Vector{<: Real}" : "" };
  const string ss_argout { dynamic ? ", steady_state" : "" };
  const int ylen {(dynamic ? 3 : 1)*symbol_table.endo_nbr()};
  const int xlen {symbol_table.exo_nbr()+symbol_table.exo_det_nbr()};

  size_t ttlen {0};

  stringstream output;

  // ResidTT!
  output << "function " << prefix << "ResidTT!(T::Vector{<: Real}, "
         << "y::Vector{<: Real}, x::Vector{<: Real}, params::Vector{<: Real}"
         << ss_argin << ")" << endl
         << "@inbounds begin" << endl
         << tt_sparse_output[0].str()
	 << "end" << endl
         << "    return nothing" << endl
         << "end" << endl << endl;
  writeToFileIfModified(output, julia_dir / (prefix + "ResidTT!.jl"));
  ttlen += temporary_terms_derivatives[0].size();

  // Resid!
  output.str("");
  output << "function " << prefix << "Resid!(T::Vector{<: Real}, residual::AbstractVector{<: Real}, "
         << "y::Vector{<: Real}, x::Vector{<: Real}, params::Vector{<: Real}"
         << ss_argin << ")" << endl
         << "    @assert length(T) >= " << ttlen << endl
         << "    @assert length(residual) == " << equations.size() << endl
         << "    @assert length(y) == " << ylen << endl
         << "    @assert length(x) == " << xlen << endl
         << "    @assert length(params) == " << symbol_table.param_nbr() << endl
         << "@inbounds begin" << endl
         << d_sparse_output[0].str()
	 << "end" << endl;
  if constexpr(!dynamic)
    output << "    if ~isreal(residual)" << endl
           << "        residual = real(residual)+imag(residual).^2;" << endl
           << "    end" << endl;
  output << "    return nothing" << endl
         << "end" << endl << endl;
  writeToFileIfModified(output, julia_dir / (prefix + "Resid!.jl"));

  // G1TT!
  output.str("");
  output << "function " << prefix << "G1TT!(T::Vector{<: Real}, y::Vector{<: Real}, "
         << "x::Vector{<: Real}, params::Vector{<: Real}" << ss_argin << ")" << endl
         << "    " << prefix << "ResidTT!(T, y, x, params" << ss_argout << ")" << endl
         << "@inbounds begin" << endl
         << tt_sparse_output[1].str()
	 << "end" << endl
         << "    return nothing" << endl
         << "end" << endl << endl;
  writeToFileIfModified(output, julia_dir / (prefix + "G1TT!.jl"));
  ttlen += temporary_terms_derivatives[1].size();

  // G1!
  output.str("");
  output << "function " << prefix << "G1!(T::Vector{<: Real}, g1_v::Vector{<: Real}, "
         << "y::Vector{<: Real}, x::Vector{<: Real}, params::Vector{<: Real}"
         << ss_argin << ")" << endl
         << "    @assert length(T) >= " << ttlen << endl
         << "    @assert length(g1_v) == " << derivatives[1].size() << endl
         << "    @assert length(y) == " << ylen << endl
         << "    @assert length(x) == " << xlen << endl
         << "    @assert length(params) == " << symbol_table.param_nbr() << endl
         << "@inbounds begin" << endl
         << d_sparse_output[1].str()
	 << "end" << endl;
  if constexpr(!dynamic)
    output << "    if ~isreal(g1_v)" << endl
           << "        g1_v = real(g1_v)+2*imag(g1_v);" << endl
           << "    end" << endl;
  output << "    return nothing" << endl
         << "end" << endl << endl;
  writeToFileIfModified(output, julia_dir / (prefix + "G1!.jl"));

  for (int i {2}; i <= computed_derivs_order; i++)
    {
      // G<i>TT!
      output.str("");
      output << "function " << prefix << "G" << i << "TT!(T::Vector{<: Real}, y::Vector{<: Real}, "
             << "x::Vector{<: Real}, params::Vector{<: Real}" << ss_argin << ")" << endl
             << "    " << prefix << "G" << to_string(i-1) << "TT!(T, y, x, params" << ss_argout << ")" << endl
             << "@inbounds begin" << endl
             << tt_sparse_output[i].str()
             << "end" << endl
             << "    return nothing" << endl
             << "end" << endl << endl;
      writeToFileIfModified(output, julia_dir / (prefix + "G" + to_string(i) + "TT!.jl"));
      ttlen += temporary_terms_derivatives[i].size();

      // G<i>!
      output.str("");
      output << "function " << prefix << "G" << i << "!(T::Vector{<: Real}, g" << i << "_v::Vector{<: Real}, "
             << "y::Vector{<: Real}, x::Vector{<: Real}, params::Vector{<: Real}"
             << ss_argin << ")" << endl
             << "    @assert length(T) >= " << ttlen << endl
             << "    @assert length(g" << i << "_v) == " << derivatives[i].size() << endl
             << "    @assert length(y) == " << ylen << endl
             << "    @assert length(x) == " << xlen << endl
             << "    @assert length(params) == " << symbol_table.param_nbr() << endl
             << "@inbounds begin" << endl
             << d_sparse_output[i].str()
             << "end" << endl
             << "    return nothing" << endl
             << "end" << endl << endl;
      writeToFileIfModified(output, julia_dir / (prefix + "G" + to_string(i) + "!.jl"));
    }
}
#endif
