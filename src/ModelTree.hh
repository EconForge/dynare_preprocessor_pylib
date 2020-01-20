/*
 * Copyright © 2003-2020 Dynare Team
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
using equation_type_and_normalized_equation_t = vector<pair<EquationType, expr_t>>;

//! Vector describing variables: max_lag in the block, max_lead in the block
using lag_lead_vector_t = vector<pair< int, int>>;

//! for each block contains tuple<Simulation_Type, first_equation, Block_Size, Recursive_part_Size>
using block_type_firstequation_size_mfs_t = vector<tuple<BlockSimulationType, int, int, int>>;

//! for a block contains derivatives tuple<block_equation_number, block_variable_number, lead_lag, expr_t>
using block_derivatives_equation_variable_laglead_nodeid_t = vector<tuple<int, int, int, expr_t>>;

//! for all blocks derivatives description
using blocks_derivatives_t = vector<block_derivatives_equation_variable_laglead_nodeid_t>;

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
  vector<pair<int, pair<string, string>>> equation_tags;
  //! Stores mapping from equation tags to equation number
  multimap<pair<string, string>, int> equation_tags_xref;
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

  //! Storage for temporary terms in block/bytecode mode
  temporary_terms_t temporary_terms;

  //! Used model local variables, that will be treated as temporary terms
  /*! See the comments in ModelTree::computeTemporaryTerms() */
  map<expr_t, expr_t, ExprNodeLess> temporary_terms_mlv;

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

  //! vector of block reordered variables and equations
  vector<int> equation_reordered, variable_reordered, inv_equation_reordered, inv_variable_reordered;

  //! the file containing the model and the derivatives code
  ofstream code_file;

  //! Vector indicating if the equation is linear in endogenous variable (true) or not (false)
  vector<bool> is_equation_linear;

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
  //! Computes temporary terms for the file containing parameters derivatives
  void computeParamsDerivativesTemporaryTerms();
  //! Writes temporary terms
  void writeTemporaryTerms(const temporary_terms_t &tt, temporary_terms_t &temp_term_union, const temporary_terms_idxs_t &tt_idxs, ostream &output, ExprNodeOutputType output_type, deriv_node_temp_terms_t &tef_terms) const;
  void writeJsonTemporaryTerms(const temporary_terms_t &tt, temporary_terms_t &temp_term_union, ostream &output, deriv_node_temp_terms_t &tef_terms, const string &concat) const;
  //! Compiles temporary terms
  void compileTemporaryTerms(ostream &code_file, unsigned int &instruction_number, const temporary_terms_t &tt, map_idx_t map_idx, bool dynamic, bool steady_dynamic) const;
  //! Adds informations for simulation in a binary file
  void Write_Inf_To_Bin_File(const string &filename, int &u_count_int, bool &file_open, bool is_two_boundaries, int block_mfs) const;
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
  void writeJsonModelLocalVariables(ostream &output, deriv_node_temp_terms_t &tef_terms) const;
  //! Compiles model equations
  void compileModelEquations(ostream &code_file, unsigned int &instruction_number, const temporary_terms_t &tt, const map_idx_t &map_idx, bool dynamic, bool steady_dynamic) const;

  //! Writes LaTeX model file
  void writeLatexModelFile(const string &mod_basename, const string &latex_basename, ExprNodeOutputType output_type, bool write_equation_tags) const;

  //! Sparse matrix of double to store the values of the Jacobian
  /*! First index is equation number, second index is endogenous type specific ID */
  using jacob_map_t = map<pair<int, int>, double>;

  //! Sparse matrix of double to store the values of the Jacobian
  /*! First index is lag, second index is equation number, third index is endogenous type specific ID */
  using dynamic_jacob_map_t = map<tuple<int, int, int>, expr_t>;

  //! Normalization of equations
  /*! Maps endogenous type specific IDs to equation numbers */
  vector<int> endo2eq;

  //! number of equation in the prologue and in the epilogue
  unsigned int epilogue, prologue;

  //! for each block contains pair< max_lag, max_lead>
  lag_lead_vector_t block_lag_lead;

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
    If no matching is found with a zero cutoff close to zero an error message is printout.
  */
  void computeNonSingularNormalization(jacob_map_t &contemporaneous_jacobian, double cutoff, jacob_map_t &static_jacobian, dynamic_jacob_map_t &dynamic_jacobian);
  //! Try to find a natural normalization if all equations are matched to an endogenous variable on the LHS
  bool computeNaturalNormalization();
  //! Try to normalized each unnormalized equation (matched endogenous variable only on the LHS)
  multimap<int, int> computeNormalizedEquations() const;
  //! Evaluate the jacobian and suppress all the elements below the cutoff
  void evaluateAndReduceJacobian(const eval_context_t &eval_context, jacob_map_t &contemporaneous_jacobian, jacob_map_t &static_jacobian, dynamic_jacob_map_t &dynamic_jacobian, double cutoff, bool verbose);
  //! Select and reorder the non linear equations of the model
  vector<pair<int, int>> select_non_linear_equations_and_variables(vector<bool> is_equation_linear, const dynamic_jacob_map_t &dynamic_jacobian, vector<int> &equation_reordered, vector<int> &variable_reordered,
                                                                   vector<int> &inv_equation_reordered, vector<int> &inv_variable_reordered,
                                                                   lag_lead_vector_t &equation_lag_lead, lag_lead_vector_t &variable_lag_lead,
                                                                   vector<unsigned int> &n_static, vector<unsigned int> &n_forward, vector<unsigned int> &n_backward, vector<unsigned int> &n_mixed);
  //! Search the equations and variables belonging to the prologue and the epilogue of the model
  void computePrologueAndEpilogue(const jacob_map_t &static_jacobian, vector<int> &equation_reordered, vector<int> &variable_reordered);
  //! Determine the type of each equation of model and try to normalized the unnormalized equation using computeNormalizedEquations
  equation_type_and_normalized_equation_t equationTypeDetermination(const map<tuple<int, int, int>, expr_t> &first_order_endo_derivatives, const vector<int> &Index_Var_IM, const vector<int> &Index_Equ_IM, int mfs) const;
  //! Compute the block decomposition and for a non-recusive block find the minimum feedback set
  void computeBlockDecompositionAndFeedbackVariablesForEachBlock(const jacob_map_t &static_jacobian, const dynamic_jacob_map_t &dynamic_jacobian, vector<int> &equation_reordered, vector<int> &variable_reordered, vector<pair<int, int>> &blocks, const equation_type_and_normalized_equation_t &Equation_Type, bool verbose_, bool select_feedback_variable, int mfs, vector<int> &inv_equation_reordered, vector<int> &inv_variable_reordered, lag_lead_vector_t &equation_lag_lead, lag_lead_vector_t &variable_lag_lead_t, vector<unsigned int> &n_static, vector<unsigned int> &n_forward, vector<unsigned int> &n_backward, vector<unsigned int> &n_mixed) const;
  //! Reduce the number of block merging the same type equation in the prologue and the epilogue and determine the type of each block
  block_type_firstequation_size_mfs_t reduceBlocksAndTypeDetermination(const dynamic_jacob_map_t &dynamic_jacobian, vector<pair<int, int>> &blocks, const equation_type_and_normalized_equation_t &Equation_Type, const vector<int> &variable_reordered, const vector<int> &equation_reordered, vector<unsigned int> &n_static, vector<unsigned int> &n_forward, vector<unsigned int> &n_backward, vector<unsigned int> &n_mixed, vector<tuple<int, int, int, int>> &block_col_type, bool linear_decomposition);
  //! Determine the maximum number of lead and lag for the endogenous variable in a bloc
  void getVariableLeadLagByBlock(const dynamic_jacob_map_t &dynamic_jacobian, const vector<int> &components_set, int nb_blck_sim, lag_lead_vector_t &equation_lead_lag, lag_lead_vector_t &variable_lead_lag, const vector<int> &equation_reordered, const vector<int> &variable_reordered) const;
  //! For each equation determine if it is linear or not
  vector<bool> equationLinear(map<tuple<int, int, int>, expr_t> first_order_endo_derivatives) const;
  //! Print an abstract of the block structure of the model
  void printBlockDecomposition(const vector<pair<int, int>> &blocks) const;
  //! Determine for each block if it is linear or not
  vector<bool> BlockLinear(const blocks_derivatives_t &blocks_derivatives, const vector<int> &variable_reordered) const;
  //! Remove equations specified by exclude_eqs
  vector<int> includeExcludeEquations(set<pair<string, string>> &eqs, bool exclude_eqs,
                                      vector<BinaryOpNode *> &equations, vector<int> &equations_lineno,
                                      vector<pair<int, pair<string, string>>> &equation_tags,
                                      multimap<pair<string, string>, int> &equation_tags_xref, bool static_equations) const;

  //! Determine the simulation type of each block
  virtual BlockSimulationType getBlockSimulationType(int block_number) const = 0;
  //! Return the number of blocks
  virtual unsigned int getNbBlocks() const = 0;
  //! Return the first equation number of a block
  virtual unsigned int getBlockFirstEquation(int block_number) const = 0;
  //! Return the size of the block block_number
  virtual unsigned int getBlockSize(int block_number) const = 0;
  //! Return the number of exogenous variable in the block block_number
  virtual unsigned int getBlockExoSize(int block_number) const = 0;
  //! Return the number of colums in the jacobian matrix for exogenous variable in the block block_number
  virtual unsigned int getBlockExoColSize(int block_number) const = 0;
  //! Return the number of feedback variable of the block block_number
  virtual unsigned int getBlockMfs(int block_number) const = 0;
  //! Return the maximum lag in a block
  virtual unsigned int getBlockMaxLag(int block_number) const = 0;
  //! Return the maximum lead in a block
  virtual unsigned int getBlockMaxLead(int block_number) const = 0;
  inline void
  setBlockLeadLag(int block, int max_lag, int max_lead)
  {
    block_lag_lead[block] = { max_lag, max_lead };
  };

  //! Return the type of equation (equation_number) belonging to the block block_number
  virtual EquationType getBlockEquationType(int block_number, int equation_number) const = 0;
  //! Return true if the equation has been normalized
  virtual bool isBlockEquationRenormalized(int block_number, int equation_number) const = 0;
  //! Return the expr_t of the equation equation_number belonging to the block block_number
  virtual expr_t getBlockEquationExpr(int block_number, int equation_number) const = 0;
  //! Return the expr_t of the renormalized equation equation_number belonging to the block block_number
  virtual expr_t getBlockEquationRenormalizedExpr(int block_number, int equation_number) const = 0;
  //! Return the original number of equation equation_number belonging to the block block_number
  virtual int getBlockEquationID(int block_number, int equation_number) const = 0;
  //! Return the original number of variable variable_number belonging to the block block_number
  virtual int getBlockVariableID(int block_number, int variable_number) const = 0;
  //! Return the original number of the exogenous variable varexo_number belonging to the block block_number
  virtual int getBlockVariableExoID(int block_number, int variable_number) const = 0;
  //! Return the position of equation_number in the block number belonging to the block block_number
  virtual int getBlockInitialEquationID(int block_number, int equation_number) const = 0;
  //! Return the position of variable_number in the block number belonging to the block block_number
  virtual int getBlockInitialVariableID(int block_number, int variable_number) const = 0;
  //! Return the position of variable_number in the block number belonging to the block block_number
  virtual int getBlockInitialExogenousID(int block_number, int variable_number) const = 0;
  //! Return the position of the deterministic exogenous variable_number in the block number belonging to the block block_number
  virtual int getBlockInitialDetExogenousID(int block_number, int variable_number) const = 0;
  //! Return the position of the other endogenous variable_number in the block number belonging to the block block_number
  virtual int getBlockInitialOtherEndogenousID(int block_number, int variable_number) const = 0;
  //! Initialize equation_reordered & variable_reordered
  void initializeVariablesAndEquations();

private:
  //! Internal helper for the copy constructor and assignment operator
  /*! Copies all the structures that contain ExprNode*, by the converting the
      pointers into their equivalent in the new tree */
  void copyHelper(const ModelTree &m);
  //! Returns the name of the MATLAB architecture given the extension used for MEX files
  static string matlab_arch(const string &mexext);
  //! Compiles the MEX file
  void compileDll(const string &basename, const string &static_or_dynamic, const string &mexext, const filesystem::path &matlabroot, const filesystem::path &dynareroot) const;

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
  void addEquation(expr_t eq, int lineno, const vector<pair<string, string>> &eq_tags);
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
  void simplifyEquations();
  /*! Reorder auxiliary variables so that they appear in recursive order in
      set_auxiliary_variables.m and dynamic_set_auxiliary_series.m */
  void reorderAuxiliaryEquations();
  //! Find equations where variable is equal to a constant
  void findConstantEquations(map<VariableNode *, NumConstNode *> &subst_table) const;
  //! Helper for writing the Jacobian elements in MATLAB and C
  /*! Writes either (i+1,j+1) or [i+j*no_eq] */
  void jacobianHelper(ostream &output, int eq_nb, int col_nb, ExprNodeOutputType output_type) const;
  //! Helper for writing the sparse Hessian or third derivatives in MATLAB and C
  /*! If order=2, writes either v2(i+1,j+1) or v2[i+j*NNZDerivatives[2]]
    If order=3, writes either v3(i+1,j+1) or v3[i+j*NNZDerivatives[3]] */
  void sparseHelper(int order, ostream &output, int row_nb, int col_nb, ExprNodeOutputType output_type) const;

  //! Returns all the equation tags associated to an equation
  inline map<string, string>
  getEquationTags(int eq) const
  {
    map<string, string> r;
    for (auto &[eq2, tagpair] : equation_tags)
      if (eq2 == eq)
        r[tagpair.first] = tagpair.second;
    return r;
  }

  inline static string
  c_Equation_Type(int type)
  {
    vector<string> c_Equation_Type =
      {
       "E_UNKNOWN   ",
       "E_EVALUATE  ",
       "E_EVALUATE_S",
       "E_SOLVE     "
      };
    return c_Equation_Type[type];
  };

  inline static string
  BlockType0(BlockType type)
  {
    switch (type)
      {
      case SIMULTANS:
        return "SIMULTANEOUS TIME SEPARABLE  ";
      case PROLOGUE:
        return "PROLOGUE                     ";
      case EPILOGUE:
        return "EPILOGUE                     ";
      case SIMULTAN:
        return "SIMULTANEOUS TIME UNSEPARABLE";
      default:
        return "UNKNOWN                      ";
      }
  };

  inline static string
  BlockSim(int type)
  {
    switch (type)
      {
      case EVALUATE_FORWARD:
        return "EVALUATE FORWARD             ";
      case EVALUATE_BACKWARD:
        return "EVALUATE BACKWARD            ";
      case SOLVE_FORWARD_SIMPLE:
        return "SOLVE FORWARD SIMPLE         ";
      case SOLVE_BACKWARD_SIMPLE:
        return "SOLVE BACKWARD SIMPLE        ";
      case SOLVE_TWO_BOUNDARIES_SIMPLE:
        return "SOLVE TWO BOUNDARIES SIMPLE  ";
      case SOLVE_FORWARD_COMPLETE:
        return "SOLVE FORWARD COMPLETE       ";
      case SOLVE_BACKWARD_COMPLETE:
        return "SOLVE BACKWARD COMPLETE      ";
      case SOLVE_TWO_BOUNDARIES_COMPLETE:
        return "SOLVE TWO BOUNDARIES COMPLETE";
      default:
        return "UNKNOWN                      ";
      }
  };
};

#endif
