/*
 * Copyright © 2003-2021 Dynare Team
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

#ifndef _MODELTREE_HH
#define _MODELTREE_HH

using namespace std;

#include <string>
#include <vector>
#include <deque>
#include <map>
#include <ostream>
#include <array>
#include <filesystem>

#include "DataTree.hh"
#include "EquationTags.hh"
#include "ExtendedPreprocessorTypes.hh"

// Helper to convert a vector into a tuple
template<typename T, size_t... Indices>
auto
vectorToTupleHelper(const vector<T> &v, index_sequence<Indices...>)
{
  return tuple(v[Indices] ...);
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
  //! Stores line numbers of declared equations; -1 means undefined
  vector<int> equations_lineno;
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
     elements several times) */
  vector<map<vector<int>, expr_t>> derivatives;

  //! Number of non-zero derivatives
  /*! Index 0 is not used, index 1 contains number of non-zero first
    derivatives, ... */
  vector<int> NNZDerivatives;

  //! Derivatives with respect to parameters
  /*! The key of the outer map is a pair (derivation order w.r.t. endogenous,
  derivation order w.r.t. parameters). For e.g., { 1, 2 } corresponds to the jacobian
  differentiated twice w.r.t. to parameters.
  In inner maps, the vector of integers consists of: the equation index, then
  the derivation IDs of endogenous (in non-decreasing order),
  then the IDs of parameters (in non-decreasing order)*/
  map<pair<int,int>, map<vector<int>, expr_t>> params_derivatives;

  //! Used model local variables, that will be treated as temporary terms
  /*! See the comments in ModelTree::computeTemporaryTerms() */
  map<VariableNode *, expr_t, ExprNodeLess> temporary_terms_mlv;

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

    inline int getRecursiveSize() const { return size - mfs_size; };
  };

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
  //! Write derivative of an equation w.r. to a variable
  void writeDerivative(ostream &output, int eq, int symb_id, int lag, ExprNodeOutputType output_type, const temporary_terms_t &temporary_terms) const;
  //! Computes temporary terms (for all equations and derivatives)
  void computeTemporaryTerms(bool is_matlab, bool no_tmp_terms);
  //! Computes temporary terms per block
  void computeBlockTemporaryTerms();
  /* Add additional temporary terms for a given block. This method is called by
     computeBlockTemporaryTerms(). It does nothing by default, but is meant to
     be overriden by subclasses (actually by DynamicModel, who needs extra
     temporary terms for derivatives w.r.t. exogenous and other endogenous) */
  virtual void additionalBlockTemporaryTerms(int blk,
                                             vector<vector<temporary_terms_t>> &blocks_temporary_terms,
                                             map<expr_t, tuple<int, int, int>> &reference_count) const;
  //! Computes temporary terms for the file containing parameters derivatives
  void computeParamsDerivativesTemporaryTerms();
  //! Writes temporary terms
  void writeTemporaryTerms(const temporary_terms_t &tt, temporary_terms_t &temp_term_union, const temporary_terms_idxs_t &tt_idxs, ostream &output, ExprNodeOutputType output_type, deriv_node_temp_terms_t &tef_terms) const;
  void writeJsonTemporaryTerms(const temporary_terms_t &tt, temporary_terms_t &temp_term_union, ostream &output, deriv_node_temp_terms_t &tef_terms, const string &concat) const;
  //! Compiles temporary terms
  void compileTemporaryTerms(ostream &code_file, unsigned int &instruction_number, bool dynamic, bool steady_dynamic, temporary_terms_t &temporary_terms_union, const temporary_terms_idxs_t &temporary_terms_idxs) const;
  //! Adds information for (non-block) bytecode simulation in a separate .bin file
  void writeBytecodeBinFile(const string &filename, int &u_count_int, bool &file_open, bool is_two_boundaries) const;
  //! Fixes output when there are more than 32 nested parens, Issue #1201
  void fixNestedParenthesis(ostringstream &output, map<string, string> &tmp_paren_vars, bool &message_printed) const;
  //! Tests if string contains more than 32 nested parens, Issue #1201
  bool testNestedParenthesis(const string &str) const;
  void writeModelLocalVariableTemporaryTerms(temporary_terms_t &temp_term_union,
                                             const temporary_terms_idxs_t &tt_idxs,
                                             ostream &output, ExprNodeOutputType output_type,
                                             deriv_node_temp_terms_t &tef_terms) const;
  //! Writes model equations
  void writeModelEquations(ostream &output, ExprNodeOutputType output_type) const;
  void writeModelEquations(ostream &output, ExprNodeOutputType output_type,
                           const temporary_terms_t &temporary_terms) const;
  //! Writes JSON model equations
  //! if residuals = true, we are writing the dynamic/static model.
  //! Otherwise, just the model equations (with line numbers, no tmp terms)
  void writeJsonModelEquations(ostream &output, bool residuals) const;
  /* Writes JSON model local variables.
     Optionally put the external function variable calls into TEF terms */
  void writeJsonModelLocalVariables(ostream &output, bool write_tef_terms, deriv_node_temp_terms_t &tef_terms) const;
  //! Compiles model equations
  void compileModelEquations(ostream &code_file, unsigned int &instruction_number, bool dynamic, bool steady_dynamic, const temporary_terms_t &temporary_terms_union, const temporary_terms_idxs_t &temporary_terms_idxs) const;

  //! Writes LaTeX model file
  void writeLatexModelFile(const string &mod_basename, const string &latex_basename, ExprNodeOutputType output_type, bool write_equation_tags) const;

  //! Sparse matrix of double to store the values of the static Jacobian
  /*! First index is equation number, second index is endogenous type specific ID */
  using jacob_map_t = map<pair<int, int>, double>;

  //! Normalization of equations, as computed by computeNonSingularNormalization()
  /*! Maps endogenous type specific IDs to equation numbers */
  vector<int> endo2eq;

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
  */
  void computeNonSingularNormalization(const jacob_map_t &contemporaneous_jacobian);
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
  //! Remove equations specified by exclude_eqs
  vector<int> includeExcludeEquations(set<pair<string, string>> &eqs, bool exclude_eqs,
                                      vector<BinaryOpNode *> &equations, vector<int> &equations_lineno,
                                      EquationTags &equation_tags, bool static_equations) const;

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
  //! Returns the 1st derivatives w.r.t. endogenous in a different format
  /*! Returns a map (equation, type-specific ID, lag) → derivative.
      Assumes that derivatives have already been computed. */
  map<tuple<int, int, int>, expr_t> collectFirstOrderDerivativesEndogenous();

private:
  //! Internal helper for the copy constructor and assignment operator
  /*! Copies all the structures that contain ExprNode*, by the converting the
      pointers into their equivalent in the new tree */
  void copyHelper(const ModelTree &m);
  //! Returns the name of the MATLAB architecture given the extension used for MEX files
  static string matlab_arch(const string &mexext);
  //! Compiles a MEX file
  void compileMEX(const string &basename, const string &funcname, const string &mexext, const vector<filesystem::path> &src_files, const filesystem::path &matlabroot, const filesystem::path &dynareroot) const;

public:
  ModelTree(SymbolTable &symbol_table_arg,
            NumericalConstants &num_constants_arg,
            ExternalFunctionsTable &external_functions_table_arg,
            bool is_dynamic_arg = false);

  ModelTree(const ModelTree &m);
  ModelTree(ModelTree &&) = delete;
  ModelTree &operator=(const ModelTree &m);
  ModelTree &operator=(ModelTree &&) = delete;

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
  void addEquation(expr_t eq, int lineno);
  //! Declare a node as an equation of the model, also giving its tags
  void addEquation(expr_t eq, int lineno, const map<string, string> &eq_tags);
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
  //! Simplify model equations: if a variable is equal to a constant, replace that variable elsewhere in the model
  /*! Equations with MCP tags are excluded, see dynare#1697 */
  void simplifyEquations();
  /*! Reorder auxiliary variables so that they appear in recursive order in
      set_auxiliary_variables.m and dynamic_set_auxiliary_series.m */
  void reorderAuxiliaryEquations();
  //! Find equations of the form “variable=constant”, excluding equations with “mcp” tag (see dynare#1697)
  void findConstantEquationsWithoutMcpTag(map<VariableNode *, NumConstNode *> &subst_table) const;
  //! Helper for writing the Jacobian elements in MATLAB and C
  /*! Writes either (i+1,j+1) or [i+j*no_eq] */
  void jacobianHelper(ostream &output, int eq_nb, int col_nb, ExprNodeOutputType output_type) const;

  //! Returns all the equation tags associated to an equation
  inline map<string, string>
  getEquationTags(int eq) const
  {
    return equation_tags.getTagsByEqn(eq);
  }

  //! Returns the vector of non-zero derivative counts
  inline const vector<int> &
  getNNZDerivatives() const
  {
    return NNZDerivatives;
  }

  //!Returns the maximum order of computed derivatives
  inline int
  getComputedDerivsOrder() const
  {
    return computed_derivs_order;
  }

  inline static string
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

  /* Returns true if at least one equation has a tag associated to occbin
     (bind/relax/pswitch/pcrit) */
  bool hasOccbinTags() const;
};

#endif
