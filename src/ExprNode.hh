/*
 * Copyright (C) 2007-2019 Dynare Team
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

#ifndef _EXPR_NODE_HH
#define _EXPR_NODE_HH

#include <set>
#include <map>
#include <vector>
#include <ostream>
#include <functional>

using namespace std;

#include "SymbolTable.hh"
#include "CodeInterpreter.hh"
#include "ExternalFunctionsTable.hh"
#include "SymbolList.hh"

class DataTree;
class NumConstNode;
class VariableNode;
class UnaryOpNode;
class BinaryOpNode;
class PacExpectationNode;

using expr_t = class ExprNode *;

struct ExprNodeLess;

//! Type for set of temporary terms
/*! The ExprNodeLess ordering is important for the temporary terms algorithm,
    see the definition of ExprNodeLess */
using temporary_terms_t = set<expr_t, ExprNodeLess>;
/*! Keeps track of array indices of temporary_terms for writing */
using temporary_terms_idxs_t = map<expr_t, int>;

//! set of temporary terms used in a block
using temporary_terms_inuse_t = set<int>;

using map_idx_t = map<int, int>;

//! Type for evaluation contexts
/*! The key is a symbol id. Lags are assumed to be null */
using eval_context_t = map<int, double>;

//! Type for tracking first/second derivative functions that have already been written as temporary terms
using deriv_node_temp_terms_t = map<pair<int, vector<expr_t>>, int>;

//! Type for the substitution map used in the process of substitutitng diff expressions
//! diff_table[static_expr_t][lag] -> [dynamic_expr_t]
using diff_table_t = map<expr_t, map<int, expr_t>>;

//! Possible types of output when writing ExprNode(s)
enum class ExprNodeOutputType
  {
    matlabStaticModel,                           //!< Matlab code, static model
    matlabDynamicModel,                          //!< Matlab code, dynamic model
    matlabStaticModelSparse,                     //!< Matlab code, static block decomposed model
    matlabDynamicModelSparse,                    //!< Matlab code, dynamic block decomposed model
    CDynamicModel,                               //!< C code, dynamic model
    CStaticModel,                                //!< C code, static model
    juliaStaticModel,                            //!< Julia code, static model
    juliaDynamicModel,                           //!< Julia code, dynamic model
    matlabOutsideModel,                          //!< Matlab code, outside model block (for example in initval)
    latexStaticModel,                            //!< LaTeX code, static model
    latexDynamicModel,                           //!< LaTeX code, dynamic model
    latexDynamicSteadyStateOperator,             //!< LaTeX code, dynamic model, inside a steady state operator
    matlabDynamicSteadyStateOperator,            //!< Matlab code, dynamic model, inside a steady state operator
    matlabDynamicSparseSteadyStateOperator,      //!< Matlab code, dynamic block decomposed model, inside a steady state operator
    CDynamicSteadyStateOperator,                 //!< C code, dynamic model, inside a steady state operator
    juliaDynamicSteadyStateOperator,             //!< Julia code, dynamic model, inside a steady state operator
    steadyStateFile,                             //!< Matlab code, in the generated steady state file
    juliaSteadyStateFile,                        //!< Julia code, in the generated steady state file
    matlabDseries,                               //!< Matlab code for dseries
    epilogueFile                                 //!< Matlab code, in the generated epilogue file
  };

inline bool
isMatlabOutput(ExprNodeOutputType output_type)
{
  return output_type == ExprNodeOutputType::matlabStaticModel
    || output_type == ExprNodeOutputType::matlabDynamicModel
    || output_type == ExprNodeOutputType::matlabOutsideModel
    || output_type == ExprNodeOutputType::matlabStaticModelSparse
    || output_type == ExprNodeOutputType::matlabDynamicModelSparse
    || output_type == ExprNodeOutputType::matlabDynamicSteadyStateOperator
    || output_type == ExprNodeOutputType::matlabDynamicSparseSteadyStateOperator
    || output_type == ExprNodeOutputType::steadyStateFile
    || output_type == ExprNodeOutputType::matlabDseries
    || output_type == ExprNodeOutputType::epilogueFile;
}

inline bool
isJuliaOutput(ExprNodeOutputType output_type)
{
  return output_type == ExprNodeOutputType::juliaStaticModel
    || output_type == ExprNodeOutputType::juliaDynamicModel
    || output_type == ExprNodeOutputType::juliaDynamicSteadyStateOperator
    || output_type == ExprNodeOutputType::juliaSteadyStateFile;
}

inline bool
isCOutput(ExprNodeOutputType output_type)
{
  return output_type == ExprNodeOutputType::CDynamicModel
    || output_type == ExprNodeOutputType::CStaticModel
    || output_type == ExprNodeOutputType::CDynamicSteadyStateOperator;
}

inline bool
isLatexOutput(ExprNodeOutputType output_type)
{
  return output_type == ExprNodeOutputType::latexStaticModel
    || output_type == ExprNodeOutputType::latexDynamicModel
    || output_type == ExprNodeOutputType::latexDynamicSteadyStateOperator;
}

/* Equal to 1 for Matlab langage or Julia, or to 0 for C language. Not defined for LaTeX.
   In Matlab and Julia, array indexes begin at 1, while they begin at 0 in C */
#define ARRAY_SUBSCRIPT_OFFSET(output_type) ((int) (isMatlabOutput(output_type) || isJuliaOutput(output_type)))

// Left and right array subscript delimiters: '(' and ')' for Matlab, '[' and ']' for C
#define LEFT_ARRAY_SUBSCRIPT(output_type) (isMatlabOutput(output_type) ? '(' : '[')
#define RIGHT_ARRAY_SUBSCRIPT(output_type) (isMatlabOutput(output_type) ? ')' : ']')

// Left and right parentheses
#define LEFT_PAR(output_type) (isLatexOutput(output_type) ? "\\left(" : "(")
#define RIGHT_PAR(output_type) (isLatexOutput(output_type) ? "\\right)" : ")")

//! Base class for expression nodes
class ExprNode
    {
      friend class DataTree;
      friend class DynamicModel;
      friend class StaticModel;
      friend class ModelTree;
      friend struct ExprNodeLess;
      friend class NumConstNode;
      friend class VariableNode;
      friend class UnaryOpNode;
      friend class BinaryOpNode;
      friend class TrinaryOpNode;
      friend class AbstractExternalFunctionNode;
      friend class VarExpectationNode;
      friend class PacExpectationNode;
    private:
      //! Computes derivative w.r. to a derivation ID (but doesn't store it in derivatives map)
      /*! You shoud use getDerivative() to get the benefit of symbolic a priori and of caching */
      virtual expr_t computeDerivative(int deriv_id) = 0;

    protected:
      //! Reference to the enclosing DataTree
      DataTree &datatree;

      //! Index number
      const int idx;

      //! Is the data member non_null_derivatives initialized ?
      bool preparedForDerivation{false};

      //! Set of derivation IDs with respect to which the derivative is potentially non-null
      set<int> non_null_derivatives;

      //! Used for caching of first order derivatives (when non-null)
      map<int, expr_t> derivatives;

      constexpr static int min_cost_matlab{40*90};
      constexpr static int min_cost_c{40*4};
      inline static int min_cost(bool is_matlab) { return(is_matlab ? min_cost_matlab : min_cost_c); };

      //! Cost of computing current node
      /*! Nodes included in temporary_terms are considered having a null cost */
      virtual int cost(int cost, bool is_matlab) const;
      virtual int cost(const temporary_terms_t &temporary_terms, bool is_matlab) const;
      virtual int cost(const map<pair<int, int>, temporary_terms_t> &temp_terms_map, bool is_matlab) const;

      //! For creating equation cross references
      struct EquationInfo
      {
        set<pair<int, int>> param;
        set<pair<int, int>> endo;
        set<pair<int, int>> exo;
        set<pair<int, int>> exo_det;
      };

      //! If this node is a temporary term, writes its temporary term representation
      /*! Returns true if node is a temporary term and has therefore been
          written to output*/
      bool checkIfTemporaryTermThenWrite(ostream &output, ExprNodeOutputType output_type,
                                         const temporary_terms_t &temporary_terms,
                                         const temporary_terms_idxs_t &temporary_terms_idxs) const;

      // Internal helper for matchVariableTimesConstantTimesParam()
      virtual void matchVTCTPHelper(int &var_id, int &lag, int &param_id, double &constant, bool at_denominator) const;

    public:
      ExprNode(DataTree &datatree_arg, int idx_arg);
      virtual ~ExprNode() = default;

      ExprNode(const ExprNode &) = delete;
      ExprNode(ExprNode &&) = delete;
      ExprNode & operator=(const ExprNode &) = delete;
      ExprNode & operator=(ExprNode &&) = delete;

      //! Initializes data member non_null_derivatives
      virtual void prepareForDerivation() = 0;

      //! Returns derivative w.r. to derivation ID
      /*! Uses a symbolic a priori to pre-detect null derivatives, and caches the result for other derivatives (to avoid computing it several times)
        For an equal node, returns the derivative of lhs minus rhs */
      expr_t getDerivative(int deriv_id);

      //! Computes derivatives by applying the chain rule for some variables
      /*!
        \param deriv_id The derivation ID with respect to which we are derivating
        \param recursive_variables Contains the derivation ID for which chain rules must be applied. Keys are derivation IDs, values are equations of the form x=f(y) where x is the key variable and x doesn't appear in y
      */
      virtual expr_t getChainRuleDerivative(int deriv_id, const map<int, expr_t> &recursive_variables) = 0;

      //! Returns precedence of node
      /*! Equals 100 for constants, variables, unary ops, and temporary terms */
      virtual int precedence(ExprNodeOutputType output_t, const temporary_terms_t &temporary_terms) const;

      //! Compute temporary terms in this expression
      /*!
        \param[in] derivOrder the derivation order (first w.r.t. endo/exo,
                   second w.r.t. params)
        \param[out] temp_terms_map the computed temporary terms, associated
                    with their derivation order
        \param[out] reference_count a temporary structure, used to count
                    references to each node (integer in outer pair is the
                    reference count, the inner pair is the derivation order)
        \param[in] is_matlab whether we are in a MATLAB context, since that
                    affects the cost of each operator

        A node will be marked as a temporary term if it is referenced at least
        two times (i.e. has at least two parents), and has a computing cost
        (multiplied by reference count) greater to datatree.min_cost
      */
      virtual void computeTemporaryTerms(const pair<int, int> &derivOrder,
                                         map<pair<int, int>, temporary_terms_t> &temp_terms_map,
                                         map<expr_t, pair<int, pair<int, int>>> &reference_count,
                                         bool is_matlab) const;

      //! Writes output of node, using a Txxx notation for nodes in temporary_terms, and specifiying the set of already written external functions
      /*!
        \param[in] output the output stream
        \param[in] output_type the type of output (MATLAB, C, LaTeX...)
        \param[in] temporary_terms the nodes that are marked as temporary terms
        \param[in] a map from temporary_terms to integers indexes (in the
                   MATLAB, C or Julia vector of temporary terms); can be empty
                   when writing MATLAB with block decomposition)
        \param[in] tef_terms the set of already written external function nodes
      */
      virtual void writeOutput(ostream &output, ExprNodeOutputType output_type, const temporary_terms_t &temporary_terms, const temporary_terms_idxs_t &temporary_terms_idxs, const deriv_node_temp_terms_t &tef_terms) const = 0;

      //! returns true if the expr node contains an external function
      virtual bool containsExternalFunction() const = 0;

      //! Writes output of node (with no temporary terms and with "outside model" output type)
      void writeOutput(ostream &output) const;

      //! Writes output of node (with no temporary terms)
      void writeOutput(ostream &output, ExprNodeOutputType output_type) const;

      //! Writes output of node, using a Txxx notation for nodes in temporary_terms
      void writeOutput(ostream &output, ExprNodeOutputType output_type, const temporary_terms_t &temporary_terms, const temporary_terms_idxs_t &temporary_terms_idxs) const;

      //! Writes output of node in JSON syntax
      virtual void writeJsonOutput(ostream &output, const temporary_terms_t &temporary_terms, const deriv_node_temp_terms_t &tef_terms, const bool isdynamic = true) const = 0;

      //! Writes the Abstract Syntax Tree in JSON
      virtual void writeJsonAST(ostream &output) const = 0;

      virtual int precedenceJson(const temporary_terms_t &temporary_terms) const;

      //! Writes the output for an external function, ensuring that the external function is called as few times as possible using temporary terms
      virtual void writeExternalFunctionOutput(ostream &output, ExprNodeOutputType output_type,
                                               const temporary_terms_t &temporary_terms,
                                               const temporary_terms_idxs_t &temporary_terms_idxs,
                                               deriv_node_temp_terms_t &tef_terms) const;

      //! Write the JSON output of an external function in a string vector
      //! Allows the insertion of commas if necessary
      virtual void writeJsonExternalFunctionOutput(vector<string> &efout,
                                                   const temporary_terms_t &temporary_terms,
                                                   deriv_node_temp_terms_t &tef_terms,
                                                   const bool isdynamic = true) const;

      virtual void compileExternalFunctionOutput(ostream &CompileCode, unsigned int &instruction_number,
                                                 bool lhs_rhs, const temporary_terms_t &temporary_terms,
                                                 const map_idx_t &map_idx, bool dynamic, bool steady_dynamic,
                                                 deriv_node_temp_terms_t &tef_terms) const;

      //! Computes the set of all variables of a given symbol type in the expression (with information on lags)
      /*!
        Variables are stored as integer pairs of the form (symb_id, lag).
        They are added to the set given in argument.
        Note that model local variables are substituted by their expression in the computation
        (and added if type_arg = ModelLocalVariable).
      */
      virtual void collectDynamicVariables(SymbolType type_arg, set<pair<int, int>> &result) const = 0;

      //! Find lowest lag for VAR
      virtual int VarMinLag() const = 0;

      //! Find the maximum lag in a VAR: handles case where LHS is diff
      virtual int VarMaxLag(DataTree &static_datatree, set<expr_t> &static_lhs) const = 0;

      //! Finds LHS variable in a VAR equation
      virtual void collectVARLHSVariable(set<expr_t> &result) const = 0;

      //! Computes the set of all variables of a given symbol type in the expression (without information on lags)
      /*!
        Variables are stored as symb_id.
        They are added to the set given in argument.
        Note that model local variables are substituted by their expression in the computation
        (and added if type_arg = ModelLocalVariable).
      */
      void collectVariables(SymbolType type_arg, set<int> &result) const;

      //! Computes the set of endogenous variables in the expression
      /*!
        Endogenous are stored as integer pairs of the form (type_specific_id, lag).
        They are added to the set given in argument.
        Note that model local variables are substituted by their expression in the computation.
      */
      virtual void collectEndogenous(set<pair<int, int>> &result) const;

      //! Computes the set of exogenous variables in the expression
      /*!
        Exogenous are stored as integer pairs of the form (type_specific_id, lag).
        They are added to the set given in argument.
        Note that model local variables are substituted by their expression in the computation.
      */
      virtual void collectExogenous(set<pair<int, int>> &result) const;

      virtual void collectTemporary_terms(const temporary_terms_t &temporary_terms, temporary_terms_inuse_t &temporary_terms_inuse, int Curr_Block) const = 0;

      virtual void computeTemporaryTerms(map<expr_t, int> &reference_count,
                                         temporary_terms_t &temporary_terms,
                                         map<expr_t, pair<int, int>> &first_occurence,
                                         int Curr_block,
                                         vector< vector<temporary_terms_t>> &v_temporary_terms,
                                         int equation) const;

      class EvalException

      {
      };

      class EvalExternalFunctionException : public EvalException

      {
      };

      virtual double eval(const eval_context_t &eval_context) const noexcept(false) = 0;
      virtual void compile(ostream &CompileCode, unsigned int &instruction_number, bool lhs_rhs, const temporary_terms_t &temporary_terms, const map_idx_t &map_idx, bool dynamic, bool steady_dynamic, const deriv_node_temp_terms_t &tef_terms) const = 0;
      void compile(ostream &CompileCode, unsigned int &instruction_number, bool lhs_rhs, const temporary_terms_t &temporary_terms, const map_idx_t &map_idx, bool dynamic, bool steady_dynamic) const;
      //! Creates a static version of this node
      /*!
        This method duplicates the current node by creating a similar node from which all leads/lags have been stripped,
        adds the result in the static_datatree argument (and not in the original datatree), and returns it.
      */
      virtual expr_t toStatic(DataTree &static_datatree) const = 0;

      /*!
        Compute cross references for equations
      */
      //  virtual void computeXrefs(set<int> &param, set<int> &endo, set<int> &exo, set<int> &exo_det) const = 0;
      virtual void computeXrefs(EquationInfo &ei) const = 0;
      //! Try to normalize an equation linear in its endogenous variable
      virtual pair<int, expr_t> normalizeEquation(int symb_id_endo, vector<tuple<int, expr_t, expr_t>> &List_of_Op_RHS) const = 0;

      //! Returns the maximum lead of endogenous in this expression
      /*! Always returns a non-negative value */
      virtual int maxEndoLead() const = 0;

      //! Returns the maximum lead of exogenous in this expression
      /*! Always returns a non-negative value */
      virtual int maxExoLead() const = 0;

      //! Returns the maximum lag of endogenous in this expression
      /*! Always returns a non-negative value */
      virtual int maxEndoLag() const = 0;

      //! Returns the maximum lag of exogenous in this expression
      /*! Always returns a non-negative value */
      virtual int maxExoLag() const = 0;

      //! Returns the maximum lead of endo/exo/exodet in this expression
      /*! A negative value means that the expression contains only lagged
          variables. */
      virtual int maxLead() const = 0;

      //! Returns the maximum lag of endo/exo/exodet in this expression
      /*! A negative value means that the expression contains only leaded
          variables. */
      virtual int maxLag() const = 0;

      //! Returns the maximum lag of endo/exo/exodet, as if diffs were expanded
      /*! This function behaves as maxLag(), except that it treats diff()
          differently. For e.g., on diff(diff(x(-1))), maxLag() returns 1 while
          maxLagWithDiffsExpanded() returns 3. */
      virtual int maxLagWithDiffsExpanded() const = 0;

      //! Get Max lag of var associated with Pac model
      //! Takes account of undiffed LHS variables in calculating the max lag
      virtual int PacMaxLag(int lhs_symb_id) const = 0;

      //! Get the target variable of the PAC model
      virtual int getPacTargetSymbId(int lhs_symb_id, int undiff_lhs_symb_id) const = 0;

      virtual expr_t undiff() const = 0;

      //! Returns a new expression where all the leads/lags have been shifted backwards by the same amount
      /*!
        Only acts on endogenous, exogenous, exogenous det
        \param[in] n The number of lags by which to shift
        \return The same expression except that leads/lags have been shifted backwards
      */
      virtual expr_t decreaseLeadsLags(int n) const = 0;

      //! Type for the substitution map used in the process of creating auxiliary vars
      using subst_table_t = map<const ExprNode *, const VariableNode *>;

      //! Type for the substitution map used in the process of substituting adl expressions
      using subst_table_adl_t = map<const ExprNode *, const expr_t>;

      //! Creates auxiliary endo lead variables corresponding to this expression
      /*!
        If maximum endogenous lead >= 3, this method will also create intermediary auxiliary var, and will add the equations of the form aux1 = aux2(+1) to the substitution table.
        \pre This expression is assumed to have maximum endogenous lead >= 2
        \param[in,out] subst_table The table to which new auxiliary variables and their correspondance will be added
        \param[out] neweqs Equations to be added to the model to match the creation of auxiliary variables.
        \return The new variable node corresponding to the current expression
      */
      VariableNode *createEndoLeadAuxiliaryVarForMyself(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const;

      //! Creates auxiliary exo lead variables corresponding to this expression
      /*!
        If maximum exogenous lead >= 2, this method will also create intermediary auxiliary var, and will add the equations of the form aux1 = aux2(+1) to the substitution table.
        \pre This expression is assumed to have maximum exogenous lead >= 1
        \param[in,out] subst_table The table to which new auxiliary variables and their correspondance will be added
        \param[out] neweqs Equations to be added to the model to match the creation of auxiliary variables.
        \return The new variable node corresponding to the current expression
      */
      VariableNode *createExoLeadAuxiliaryVarForMyself(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const;

      //! Constructs a new expression where sub-expressions with max endo lead >= 2 have been replaced by auxiliary variables
      /*!
        \param[in,out] subst_table Map used to store expressions that have already be substituted and their corresponding variable, in order to avoid creating two auxiliary variables for the same sub-expr.
        \param[out] neweqs Equations to be added to the model to match the creation of auxiliary variables.

        If the method detects a sub-expr which needs to be substituted, two cases are possible:
        - if this expr is in the table, then it will use the corresponding variable and return the substituted expression
        - if this expr is not in the table, then it will create an auxiliary endogenous variable, add the substitution in the table and return the substituted expression

        \return A new equivalent expression where sub-expressions with max endo lead >= 2 have been replaced by auxiliary variables
      */
      virtual expr_t substituteEndoLeadGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const = 0;

      //! Constructs a new expression where endo variables with max endo lag >= 2 have been replaced by auxiliary variables
      /*!
        \param[in,out] subst_table Map used to store expressions that have already be substituted and their corresponding variable, in order to avoid creating two auxiliary variables for the same sub-expr.
        \param[out] neweqs Equations to be added to the model to match the creation of auxiliary variables.
      */
      virtual expr_t substituteEndoLagGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const = 0;

      //! Constructs a new expression where exogenous variables with a lead have been replaced by auxiliary variables
      /*!
        \param[in,out] subst_table Map used to store expressions that have already be substituted and their corresponding variable, in order to avoid creating two auxiliary variables for the same sub-expr.
        \param[out] neweqs Equations to be added to the model to match the creation of auxiliary variables.
      */
      virtual expr_t substituteExoLead(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const = 0;
      //! Constructs a new expression where exogenous variables with a lag have been replaced by auxiliary variables
      /*!
        \param[in,out] subst_table Map used to store expressions that have already be substituted and their corresponding variable, in order to avoid creating two auxiliary variables for the same sub-expr.
        \param[out] neweqs Equations to be added to the model to match the creation of auxiliary variables.
      */
      virtual expr_t substituteExoLag(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const = 0;

      //! Constructs a new expression where the expectation operator has been replaced by auxiliary variables
      /*!
        \param[in,out] subst_table Map used to store expressions that have already be substituted and their corresponding variable, in order to avoid creating two auxiliary variables for the same sub-expr.
        \param[out] neweqs Equations to be added to the model to match the creation of auxiliary variables.
        \param[in] partial_information_model Are we substituting in a partial information model?
      */
      virtual expr_t substituteExpectation(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool partial_information_model) const = 0;

      virtual expr_t decreaseLeadsLagsPredeterminedVariables() const = 0;

      //! Constructs a new expression where forward variables (supposed to be at most in t+1) have been replaced by themselves at t, plus a new aux var representing their (time) differentiate
      /*!
        \param[in] subset variables to which to limit the transformation; transform
        all fwrd vars if empty
        \param[in,out] subst_table Map used to store mapping between a given
        forward variable and the aux var that contains its differentiate
        \param[out] neweqs Equations to be added to the model to match the creation of auxiliary variables.
      */
      virtual expr_t differentiateForwardVars(const vector<string> &subset, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const = 0;

      //! Return true if the nodeID is a numerical constant equal to value and false otherwise
      /*!
        \param[in] value of the numerical constante
        \param[out] the boolean equal to true if NodeId is a constant equal to value
      */
      virtual bool isNumConstNodeEqualTo(double value) const = 0;

      //! Returns true if the expression contains one or several endogenous variable
      virtual bool containsEndogenous() const = 0;

      //! Returns true if the expression contains one or several exogenous variable
      virtual bool containsExogenous() const = 0;

      //! Returns the maximum number of nested diffs in the expression
      virtual int countDiffs() const = 0;

      //! Return true if the nodeID is a variable withe a type equal to type_arg, a specific variable id aqual to varfiable_id and a lag equal to lag_arg and false otherwise
      /*!
        \param[in] the type (type_arg), specifique variable id (variable_id and the lag (lag_arg)
        \param[out] the boolean equal to true if NodeId is the variable
      */
      virtual bool isVariableNodeEqualTo(SymbolType type_arg, int variable_id, int lag_arg) const = 0;

      //! Replaces the Trend var with datatree.One
      virtual expr_t replaceTrendVar() const = 0;

      //! Constructs a new expression where the variable indicated by symb_id has been detrended
      /*!
        \param[in] symb_id indicating the variable to be detrended
        \param[in] log_trend indicates if the trend is in log
        \param[in] trend indicating the trend
        \return the new binary op pointing to a detrended variable
      */
      virtual expr_t detrend(int symb_id, bool log_trend, expr_t trend) const = 0;

      //! Substitute adl operator
      virtual expr_t substituteAdl() const = 0;

      //! Substitute VarExpectation nodes
      virtual expr_t substituteVarExpectation(const map<string, expr_t> &subst_table) const = 0;

      //! Substitute diff operator
      virtual void findDiffNodes(DataTree &static_datatree, diff_table_t &diff_table) const = 0;
      virtual void findUnaryOpNodesForAuxVarCreation(DataTree &static_datatree, diff_table_t &nodes) const = 0;
      virtual int findTargetVariable(int lhs_symb_id) const = 0;
      virtual expr_t substituteDiff(DataTree &static_datatree, diff_table_t &diff_table, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const = 0;
      virtual expr_t substituteUnaryOpNodes(DataTree &static_datatree, diff_table_t &nodes, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const = 0;

      //! Substitute pac_expectation operator
      virtual expr_t substitutePacExpectation(const string & name, expr_t subexpr) = 0;

      //! Add ExprNodes to the provided datatree
      virtual expr_t clone(DataTree &datatree) const = 0;

      //! Move a trend variable with lag/lead to time t by dividing/multiplying by its growth factor
      virtual expr_t removeTrendLeadLag(map<int, expr_t> trend_symbols_map) const = 0;

      //! Returns true if the expression is in static form (no lead, no lag, no expectation, no STEADY_STATE)
      virtual bool isInStaticForm() const = 0;

      //! Substitute auxiliary variables by their expression in static model
      virtual expr_t substituteStaticAuxiliaryVariable() const = 0;

      //! Returns true if model_info_name is referenced by a VarExpectationNode
      virtual bool isVarModelReferenced(const string &model_info_name) const = 0;

      //! Matches a linear combination of variables, where scalars can be constant*parameter
      /*! Returns a list of (variable_id, lag, param_id, constant)
          corresponding to the terms in the expression. When there is no
          parameter in a term, param_id == -1.
          Can throw a MatchFailureException. */
      vector<tuple<int, int, int, double>> matchLinearCombinationOfVariables() const;

      pair<int, vector<tuple<int, int, int, double>>> matchParamTimesLinearCombinationOfVariables() const;

      //! Returns true if expression is of the form:
      //! param * (endog op endog op ...) + param * (endog op endog op ...) + ...
      virtual bool isParamTimesEndogExpr() const = 0;

      //! Fills the EC matrix structure
      void fillErrorCorrectionRow(int eqn, const vector<int> &nontarget_lhs, const vector<int> &target_lhs,
                                  map<tuple<int, int, int>, expr_t> &A0,
                                  map<tuple<int, int, int>, expr_t> &A0star) const;

      //! Finds equations where a variable is equal to a constant
      virtual void findConstantEquations(map<VariableNode *, NumConstNode *> &table) const = 0;

      //! Replaces variables found in findConstantEquations() with their constant values
      virtual expr_t replaceVarsInEquation(map<VariableNode *, NumConstNode *> &table) const = 0;

      //! Returns true if PacExpectationNode encountered
      virtual bool containsPacExpectation(const string &pac_model_name = "") const = 0;

      //! Fills map
      virtual void getEndosAndMaxLags(map<string, int> &model_endos_and_lags) const = 0;

      //! Decompose an expression into its additive terms
      /*! Returns a list of terms, with their sign (either 1 or -1, depending
        on whether the terms appears with a plus or a minus).
        The current_sign argument should normally be left to 1.
        If current_sign == -1, then all signs are inverted */
      virtual void decomposeAdditiveTerms(vector<pair<expr_t, int>> &terms, int current_sign = 1) const;

      // Matches an expression of the form variable*constant*parameter
      /* Returns a tuple (variable_id, lag, param_id, constant).
         The variable must be an exogenous or an endogenous.
         The constant is optional (in which case 1 is returned); there can be
         several multiplicative constants; constants can also appear at the
         denominator (i.e. after a divide sign).
         The parameter is optional (in which case param_id == -1).
         If the expression is not of the expected form, throws a
         MatchFailureException */
      tuple<int, int, int, double> matchVariableTimesConstantTimesParam() const;

      //! Exception thrown by matchVariableTimesConstantTimesParam when matching fails
      class MatchFailureException
      {
      public:
        const string message;
        MatchFailureException(string message_arg) : message{move(message_arg)} {};
      };
};

//! Object used to compare two nodes (using their indexes)
/*! Note that in this ordering, a subexpression is always less than the
    expression from which it is extracted. This property is used extensively in
    the temporary terms computations. */
struct ExprNodeLess
{
  bool
  operator()(expr_t arg1, expr_t arg2) const
  {
    return arg1->idx < arg2->idx;
  }
};

//! Numerical constant node
/*! The constant is necessarily non-negative (this is enforced at the NumericalConstants class level) */
class NumConstNode : public ExprNode
{
public:
  //! Id from numerical constants table
  const int id;
private:
  expr_t computeDerivative(int deriv_id) override;
protected:
  void matchVTCTPHelper(int &var_id, int &lag, int &param_id, double &constant, bool at_denominator) const override;
public:
  NumConstNode(DataTree &datatree_arg, int idx_arg, int id_arg);
  void prepareForDerivation() override;
  void writeOutput(ostream &output, ExprNodeOutputType output_type, const temporary_terms_t &temporary_terms, const temporary_terms_idxs_t &temporary_terms_idxs, const deriv_node_temp_terms_t &tef_terms) const override;
  void writeJsonAST(ostream &output) const override;
  void writeJsonOutput(ostream &output, const temporary_terms_t &temporary_terms, const deriv_node_temp_terms_t &tef_terms, const bool isdynamic) const override;
  bool containsExternalFunction() const override;
  void collectVARLHSVariable(set<expr_t> &result) const override;
  void collectDynamicVariables(SymbolType type_arg, set<pair<int, int>> &result) const override;
  void collectTemporary_terms(const temporary_terms_t &temporary_terms, temporary_terms_inuse_t &temporary_terms_inuse, int Curr_Block) const override;
  double eval(const eval_context_t &eval_context) const noexcept(false) override;
  void compile(ostream &CompileCode, unsigned int &instruction_number, bool lhs_rhs, const temporary_terms_t &temporary_terms, const map_idx_t &map_idx, bool dynamic, bool steady_dynamic, const deriv_node_temp_terms_t &tef_terms) const override;
  expr_t toStatic(DataTree &static_datatree) const override;
  void computeXrefs(EquationInfo &ei) const override;
  pair<int, expr_t> normalizeEquation(int symb_id_endo, vector<tuple<int, expr_t, expr_t>>  &List_of_Op_RHS) const override;
  expr_t getChainRuleDerivative(int deriv_id, const map<int, expr_t> &recursive_variables) override;
  int maxEndoLead() const override;
  int maxExoLead() const override;
  int maxEndoLag() const override;
  int maxExoLag() const override;
  int maxLead() const override;
  int maxLag() const override;
  int maxLagWithDiffsExpanded() const override;
  int VarMinLag() const override;
  int VarMaxLag(DataTree &static_datatree, set<expr_t> &static_lhs) const override;
  int PacMaxLag(int lhs_symb_id) const override;
  int getPacTargetSymbId(int lhs_symb_id, int undiff_lhs_symb_id) const override;
  expr_t undiff() const override;
  expr_t decreaseLeadsLags(int n) const override;
  expr_t substituteEndoLeadGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const override;
  expr_t substituteEndoLagGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const override;
  expr_t substituteExoLead(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const override;
  expr_t substituteExoLag(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const override;
  expr_t substituteExpectation(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool partial_information_model) const override;
  expr_t substituteAdl() const override;
  expr_t substituteVarExpectation(const map<string, expr_t> &subst_table) const override;
  void findDiffNodes(DataTree &static_datatree, diff_table_t &diff_table) const override;
  void findUnaryOpNodesForAuxVarCreation(DataTree &static_datatree, diff_table_t &nodes) const override;
  int findTargetVariable(int lhs_symb_id) const override;
  expr_t substituteDiff(DataTree &static_datatree, diff_table_t &diff_table, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const override;
  expr_t substituteUnaryOpNodes(DataTree &static_datatree, diff_table_t &nodes, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const override;
  expr_t substitutePacExpectation(const string & name, expr_t subexpr) override;
  expr_t decreaseLeadsLagsPredeterminedVariables() const override;
  expr_t differentiateForwardVars(const vector<string> &subset, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const override;
  bool isNumConstNodeEqualTo(double value) const override;
  bool containsEndogenous() const override;
  bool containsExogenous() const override;
  int countDiffs() const override;
  bool isVariableNodeEqualTo(SymbolType type_arg, int variable_id, int lag_arg) const override;
  expr_t replaceTrendVar() const override;
  expr_t detrend(int symb_id, bool log_trend, expr_t trend) const override;
  expr_t clone(DataTree &datatree) const override;
  expr_t removeTrendLeadLag(map<int, expr_t> trend_symbols_map) const override;
  bool isInStaticForm() const override;
  void findConstantEquations(map<VariableNode *, NumConstNode *> &table) const override;
  expr_t replaceVarsInEquation(map<VariableNode *, NumConstNode *> &table) const override;
  bool containsPacExpectation(const string &pac_model_name = "") const override;
  bool isParamTimesEndogExpr() const override;
  bool isVarModelReferenced(const string &model_info_name) const override;
  void getEndosAndMaxLags(map<string, int> &model_endos_and_lags) const override;
  expr_t substituteStaticAuxiliaryVariable() const override;
};

//! Symbol or variable node
class VariableNode : public ExprNode
{
  friend class UnaryOpNode;
public:
  //! Id from the symbol table
  const int symb_id;
  //! A positive value is a lead, a negative is a lag
  const int lag;
private:
  expr_t computeDerivative(int deriv_id) override;
protected:
  void matchVTCTPHelper(int &var_id, int &lag, int &param_id, double &constant, bool at_denominator) const override;
public:
  VariableNode(DataTree &datatree_arg, int idx_arg, int symb_id_arg, int lag_arg);
  void prepareForDerivation() override;
  void writeOutput(ostream &output, ExprNodeOutputType output_type, const temporary_terms_t &temporary_terms, const temporary_terms_idxs_t &temporary_terms_idxs, const deriv_node_temp_terms_t &tef_terms) const override;
  void writeJsonAST(ostream &output) const override;
  void writeJsonOutput(ostream &output, const temporary_terms_t &temporary_terms, const deriv_node_temp_terms_t &tef_terms, const bool isdynamic) const override;
  bool containsExternalFunction() const override;
  void collectVARLHSVariable(set<expr_t> &result) const override;
  void collectDynamicVariables(SymbolType type_arg, set<pair<int, int>> &result) const override;
  void computeTemporaryTerms(map<expr_t, int > &reference_count,
                                     temporary_terms_t &temporary_terms,
                                     map<expr_t, pair<int, int>> &first_occurence,
                                     int Curr_block,
                                     vector< vector<temporary_terms_t>> &v_temporary_terms,
                                     int equation) const override;
  void collectTemporary_terms(const temporary_terms_t &temporary_terms, temporary_terms_inuse_t &temporary_terms_inuse, int Curr_Block) const override;
  double eval(const eval_context_t &eval_context) const noexcept(false) override;
  void compile(ostream &CompileCode, unsigned int &instruction_number, bool lhs_rhs, const temporary_terms_t &temporary_terms, const map_idx_t &map_idx, bool dynamic, bool steady_dynamic, const deriv_node_temp_terms_t &tef_terms) const override;
  expr_t toStatic(DataTree &static_datatree) const override;
  void computeXrefs(EquationInfo &ei) const override;
  SymbolType get_type() const;
  pair<int, expr_t> normalizeEquation(int symb_id_endo, vector<tuple<int, expr_t, expr_t>>  &List_of_Op_RHS) const override;
  expr_t getChainRuleDerivative(int deriv_id, const map<int, expr_t> &recursive_variables) override;
  int maxEndoLead() const override;
  int maxExoLead() const override;
  int maxEndoLag() const override;
  int maxExoLag() const override;
  int maxLead() const override;
  int maxLag() const override;
  int maxLagWithDiffsExpanded() const override;
  int VarMinLag() const override;
  int VarMaxLag(DataTree &static_datatree, set<expr_t> &static_lhs) const override;
  int PacMaxLag(int lhs_symb_id) const override;
  int getPacTargetSymbId(int lhs_symb_id, int undiff_lhs_symb_id) const override;
  expr_t undiff() const override;
  expr_t decreaseLeadsLags(int n) const override;
  expr_t substituteEndoLeadGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const override;
  expr_t substituteEndoLagGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const override;
  expr_t substituteExoLead(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const override;
  expr_t substituteExoLag(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const override;
  expr_t substituteExpectation(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool partial_information_model) const override;
  expr_t substituteAdl() const override;
  expr_t substituteVarExpectation(const map<string, expr_t> &subst_table) const override;
  void findDiffNodes(DataTree &static_datatree, diff_table_t &diff_table) const override;
  void findUnaryOpNodesForAuxVarCreation(DataTree &static_datatree, diff_table_t &nodes) const override;
  int findTargetVariable(int lhs_symb_id) const override;
  expr_t substituteDiff(DataTree &static_datatree, diff_table_t &diff_table, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const override;
  expr_t substituteUnaryOpNodes(DataTree &static_datatree, diff_table_t &nodes, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const override;
  expr_t substitutePacExpectation(const string & name, expr_t subexpr) override;
  expr_t decreaseLeadsLagsPredeterminedVariables() const override;
  expr_t differentiateForwardVars(const vector<string> &subset, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const override;
  bool isNumConstNodeEqualTo(double value) const override;
  bool containsEndogenous() const override;
  bool containsExogenous() const override;
  int countDiffs() const override;
  bool isVariableNodeEqualTo(SymbolType type_arg, int variable_id, int lag_arg) const override;
  expr_t replaceTrendVar() const override;
  expr_t detrend(int symb_id, bool log_trend, expr_t trend) const override;
  expr_t clone(DataTree &datatree) const override;
  expr_t removeTrendLeadLag(map<int, expr_t> trend_symbols_map) const override;
  bool isInStaticForm() const override;
  void findConstantEquations(map<VariableNode *, NumConstNode *> &table) const override;
  expr_t replaceVarsInEquation(map<VariableNode *, NumConstNode *> &table) const override;
  bool containsPacExpectation(const string &pac_model_name = "") const override;
  bool isParamTimesEndogExpr() const override;
  bool isVarModelReferenced(const string &model_info_name) const override;
  void getEndosAndMaxLags(map<string, int> &model_endos_and_lags) const override;
  //! Substitute auxiliary variables by their expression in static model
  expr_t substituteStaticAuxiliaryVariable() const override;
};

//! Unary operator node
class UnaryOpNode : public ExprNode
{
protected:
  void matchVTCTPHelper(int &var_id, int &lag, int &param_id, double &constant, bool at_denominator) const override;
public:
  const expr_t arg;
  //! Stores the information set. Only used for expectation operator
  const int expectation_information_set;
  //! Only used for UnaryOpcode::steadyStateParamDeriv and UnaryOpcode::steadyStateParam2ndDeriv
  const int param1_symb_id, param2_symb_id;
  const UnaryOpcode op_code;
  const string adl_param_name;
  const vector<int> adl_lags;
private:
  expr_t computeDerivative(int deriv_id) override;
  int cost(int cost, bool is_matlab) const override;
  int cost(const temporary_terms_t &temporary_terms, bool is_matlab) const override;
  int cost(const map<pair<int, int>, temporary_terms_t> &temp_terms_map, bool is_matlab) const override;
  //! Returns the derivative of this node if darg is the derivative of the argument
  expr_t composeDerivatives(expr_t darg, int deriv_id);
public:
  UnaryOpNode(DataTree &datatree_arg, int idx_arg, UnaryOpcode op_code_arg, const expr_t arg_arg, int expectation_information_set_arg, int param1_symb_id_arg, int param2_symb_id_arg, string adl_param_name_arg, vector<int> adl_lags_arg);
  void prepareForDerivation() override;
  void computeTemporaryTerms(const pair<int, int> &derivOrder,
                             map<pair<int, int>, temporary_terms_t> &temp_terms_map,
                             map<expr_t, pair<int, pair<int, int>>> &reference_count,
                             bool is_matlab) const override;
  void writeOutput(ostream &output, ExprNodeOutputType output_type, const temporary_terms_t &temporary_terms, const temporary_terms_idxs_t &temporary_terms_idxs, const deriv_node_temp_terms_t &tef_terms) const override;
  void writeJsonAST(ostream &output) const override;
  void writeJsonOutput(ostream &output, const temporary_terms_t &temporary_terms, const deriv_node_temp_terms_t &tef_terms, const bool isdynamic) const override;
  bool containsExternalFunction() const override;
  void writeExternalFunctionOutput(ostream &output, ExprNodeOutputType output_type,
                                           const temporary_terms_t &temporary_terms,
                                           const temporary_terms_idxs_t &temporary_terms_idxs,
                                           deriv_node_temp_terms_t &tef_terms) const override;
  void writeJsonExternalFunctionOutput(vector<string> &efout,
                                               const temporary_terms_t &temporary_terms,
                                               deriv_node_temp_terms_t &tef_terms,
                                               const bool isdynamic) const override;
  void compileExternalFunctionOutput(ostream &CompileCode, unsigned int &instruction_number,
                                             bool lhs_rhs, const temporary_terms_t &temporary_terms,
                                             const map_idx_t &map_idx, bool dynamic, bool steady_dynamic,
                                             deriv_node_temp_terms_t &tef_terms) const override;
  void computeTemporaryTerms(map<expr_t, int> &reference_count,
                                     temporary_terms_t &temporary_terms,
                                     map<expr_t, pair<int, int>> &first_occurence,
                                     int Curr_block,
                                     vector< vector<temporary_terms_t>> &v_temporary_terms,
                                     int equation) const override;
  void collectVARLHSVariable(set<expr_t> &result) const override;
  void collectDynamicVariables(SymbolType type_arg, set<pair<int, int>> &result) const override;
  void collectTemporary_terms(const temporary_terms_t &temporary_terms, temporary_terms_inuse_t &temporary_terms_inuse, int Curr_Block) const override;
  static double eval_opcode(UnaryOpcode op_code, double v) noexcept(false);
  double eval(const eval_context_t &eval_context) const noexcept(false) override;
  void compile(ostream &CompileCode, unsigned int &instruction_number, bool lhs_rhs, const temporary_terms_t &temporary_terms, const map_idx_t &map_idx, bool dynamic, bool steady_dynamic, const deriv_node_temp_terms_t &tef_terms) const override;
  expr_t toStatic(DataTree &static_datatree) const override;
  void computeXrefs(EquationInfo &ei) const override;
  pair<int, expr_t> normalizeEquation(int symb_id_endo, vector<tuple<int, expr_t, expr_t>>  &List_of_Op_RHS) const override;
  expr_t getChainRuleDerivative(int deriv_id, const map<int, expr_t> &recursive_variables) override;
  int maxEndoLead() const override;
  int maxExoLead() const override;
  int maxEndoLag() const override;
  int maxExoLag() const override;
  int maxLead() const override;
  int maxLag() const override;
  int maxLagWithDiffsExpanded() const override;
  int VarMinLag() const override;
  int VarMaxLag(DataTree &static_datatree, set<expr_t> &static_lhs) const override;
  int PacMaxLag(int lhs_symb_id) const override;
  int getPacTargetSymbId(int lhs_symb_id, int undiff_lhs_symb_id) const override;
  expr_t undiff() const override;
  expr_t decreaseLeadsLags(int n) const override;
  expr_t substituteEndoLeadGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const override;
  //! Creates another UnaryOpNode with the same opcode, but with a possibly different datatree and argument
  expr_t buildSimilarUnaryOpNode(expr_t alt_arg, DataTree &alt_datatree) const;
  expr_t substituteEndoLagGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const override;
  expr_t substituteExoLead(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const override;
  expr_t substituteExoLag(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const override;
  expr_t substituteExpectation(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool partial_information_model) const override;
  expr_t substituteAdl() const override;
  expr_t substituteVarExpectation(const map<string, expr_t> &subst_table) const override;
  void findDiffNodes(DataTree &static_datatree, diff_table_t &diff_table) const override;
  bool createAuxVarForUnaryOpNode() const;
  void findUnaryOpNodesForAuxVarCreation(DataTree &static_datatree, diff_table_t &nodes) const override;
  int findTargetVariable(int lhs_symb_id) const override;
  expr_t substituteDiff(DataTree &static_datatree, diff_table_t &diff_table, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const override;
  expr_t substituteUnaryOpNodes(DataTree &static_datatree, diff_table_t &nodes, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const override;
  expr_t substitutePacExpectation(const string & name, expr_t subexpr) override;
  expr_t decreaseLeadsLagsPredeterminedVariables() const override;
  expr_t differentiateForwardVars(const vector<string> &subset, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const override;
  bool isNumConstNodeEqualTo(double value) const override;
  bool containsEndogenous() const override;
  bool containsExogenous() const override;
  int countDiffs() const override;
  bool isVariableNodeEqualTo(SymbolType type_arg, int variable_id, int lag_arg) const override;
  expr_t replaceTrendVar() const override;
  expr_t detrend(int symb_id, bool log_trend, expr_t trend) const override;
  expr_t clone(DataTree &datatree) const override;
  expr_t removeTrendLeadLag(map<int, expr_t> trend_symbols_map) const override;
  bool isInStaticForm() const override;
  void findConstantEquations(map<VariableNode *, NumConstNode *> &table) const override;
  expr_t replaceVarsInEquation(map<VariableNode *, NumConstNode *> &table) const override;
  bool containsPacExpectation(const string &pac_model_name = "") const override;
  bool isParamTimesEndogExpr() const override;
  bool isVarModelReferenced(const string &model_info_name) const override;
  void getEndosAndMaxLags(map<string, int> &model_endos_and_lags) const override;
  //! Substitute auxiliary variables by their expression in static model
  expr_t substituteStaticAuxiliaryVariable() const override;
  void decomposeAdditiveTerms(vector<pair<expr_t, int>> &terms, int current_sign) const override;
};

//! Binary operator node
class BinaryOpNode : public ExprNode
{
protected:
  void matchVTCTPHelper(int &var_id, int &lag, int &param_id, double &constant, bool at_denominator) const override;
public:
  const expr_t arg1, arg2;
  const BinaryOpcode op_code;
  const int powerDerivOrder;
  const string adlparam;
private:
  expr_t computeDerivative(int deriv_id) override;
  int cost(int cost, bool is_matlab) const override;
  int cost(const temporary_terms_t &temporary_terms, bool is_matlab) const override;
  int cost(const map<pair<int, int>, temporary_terms_t> &temp_terms_map, bool is_matlab) const override;
  //! Returns the derivative of this node if darg1 and darg2 are the derivatives of the arguments
  expr_t composeDerivatives(expr_t darg1, expr_t darg2);
public:
  BinaryOpNode(DataTree &datatree_arg, int idx_arg, const expr_t arg1_arg,
               BinaryOpcode op_code_arg, const expr_t arg2_arg, int powerDerivOrder);
  void prepareForDerivation() override;
  int precedenceJson(const temporary_terms_t &temporary_terms) const override;
  int precedence(ExprNodeOutputType output_type, const temporary_terms_t &temporary_terms) const override;
  void computeTemporaryTerms(const pair<int, int> &derivOrder,
                             map<pair<int, int>, temporary_terms_t> &temp_terms_map,
                             map<expr_t, pair<int, pair<int, int>>> &reference_count,
                             bool is_matlab) const override;
  void writeOutput(ostream &output, ExprNodeOutputType output_type, const temporary_terms_t &temporary_terms, const temporary_terms_idxs_t &temporary_terms_idxs, const deriv_node_temp_terms_t &tef_terms) const override;
  void writeJsonAST(ostream &output) const override;
  void writeJsonOutput(ostream &output, const temporary_terms_t &temporary_terms, const deriv_node_temp_terms_t &tef_terms, const bool isdynamic) const override;
  bool containsExternalFunction() const override;
  void writeExternalFunctionOutput(ostream &output, ExprNodeOutputType output_type,
                                           const temporary_terms_t &temporary_terms,
                                           const temporary_terms_idxs_t &temporary_terms_idxs,
                                           deriv_node_temp_terms_t &tef_terms) const override;
  void writeJsonExternalFunctionOutput(vector<string> &efout,
                                               const temporary_terms_t &temporary_terms,
                                               deriv_node_temp_terms_t &tef_terms,
                                               const bool isdynamic) const override;
  void compileExternalFunctionOutput(ostream &CompileCode, unsigned int &instruction_number,
                                             bool lhs_rhs, const temporary_terms_t &temporary_terms,
                                             const map_idx_t &map_idx, bool dynamic, bool steady_dynamic,
                                             deriv_node_temp_terms_t &tef_terms) const override;
  void computeTemporaryTerms(map<expr_t, int> &reference_count,
                                     temporary_terms_t &temporary_terms,
                                     map<expr_t, pair<int, int>> &first_occurence,
                                     int Curr_block,
                                     vector< vector<temporary_terms_t>> &v_temporary_terms,
                                     int equation) const override;
  void collectVARLHSVariable(set<expr_t> &result) const override;
  void collectDynamicVariables(SymbolType type_arg, set<pair<int, int>> &result) const override;
  void collectTemporary_terms(const temporary_terms_t &temporary_terms, temporary_terms_inuse_t &temporary_terms_inuse, int Curr_Block) const override;
  static double eval_opcode(double v1, BinaryOpcode op_code, double v2, int derivOrder) noexcept(false);
  double eval(const eval_context_t &eval_context) const noexcept(false) override;
  void compile(ostream &CompileCode, unsigned int &instruction_number, bool lhs_rhs, const temporary_terms_t &temporary_terms, const map_idx_t &map_idx, bool dynamic, bool steady_dynamic, const deriv_node_temp_terms_t &tef_terms) const override;
  expr_t Compute_RHS(expr_t arg1, expr_t arg2, int op, int op_type) const;
  expr_t toStatic(DataTree &static_datatree) const override;
  void computeXrefs(EquationInfo &ei) const override;
  pair<int, expr_t> normalizeEquation(int symb_id_endo, vector<tuple<int, expr_t, expr_t>>  &List_of_Op_RHS) const override;
  expr_t getChainRuleDerivative(int deriv_id, const map<int, expr_t> &recursive_variables) override;
  int maxEndoLead() const override;
  int maxExoLead() const override;
  int maxEndoLag() const override;
  int maxExoLag() const override;
  int maxLead() const override;
  int maxLag() const override;
  int maxLagWithDiffsExpanded() const override;
  int VarMinLag() const override;
  int VarMaxLag(DataTree &static_datatree, set<expr_t> &static_lhs) const override;
  int PacMaxLag(int lhs_symb_id) const override;
  int getPacTargetSymbIdHelper(int lhs_symb_id, int undiff_lhs_symb_id, const set<pair<int, int>> & endogs) const;
  int getPacTargetSymbId(int lhs_symb_id, int undiff_lhs_symb_id) const override;
  expr_t undiff() const override;
  expr_t decreaseLeadsLags(int n) const override;
  expr_t substituteEndoLeadGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const override;
  //! Creates another BinaryOpNode with the same opcode, but with a possibly different datatree and arguments
  expr_t buildSimilarBinaryOpNode(expr_t alt_arg1, expr_t alt_arg2, DataTree &alt_datatree) const;
  expr_t substituteEndoLagGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const override;
  expr_t substituteExoLead(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const override;
  expr_t substituteExoLag(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const override;
  expr_t substituteExpectation(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool partial_information_model) const override;
  expr_t substituteAdl() const override;
  expr_t substituteVarExpectation(const map<string, expr_t> &subst_table) const override;
  void findDiffNodes(DataTree &static_datatree, diff_table_t &diff_table) const override;
  void findUnaryOpNodesForAuxVarCreation(DataTree &static_datatree, diff_table_t &nodes) const override;
  bool findTargetVariableHelper1(int lhs_symb_id, int rhs_symb_id) const;
  int findTargetVariableHelper(const expr_t arg1, const expr_t arg2, int lhs_symb_id) const;
  int findTargetVariable(int lhs_symb_id) const override;
  expr_t substituteDiff(DataTree &static_datatree, diff_table_t &diff_table, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const override;
  expr_t substituteUnaryOpNodes(DataTree &static_datatree, diff_table_t &nodes, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const override;
  expr_t substitutePacExpectation(const string & name, expr_t subexpr) override;
  expr_t decreaseLeadsLagsPredeterminedVariables() const override;
  expr_t differentiateForwardVars(const vector<string> &subset, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const override;
  bool isNumConstNodeEqualTo(double value) const override;
  bool containsEndogenous() const override;
  bool containsExogenous() const override;
  int countDiffs() const override;
  bool isVariableNodeEqualTo(SymbolType type_arg, int variable_id, int lag_arg) const override;
  expr_t replaceTrendVar() const override;
  expr_t detrend(int symb_id, bool log_trend, expr_t trend) const override;
  expr_t clone(DataTree &datatree) const override;
  expr_t removeTrendLeadLag(map<int, expr_t> trend_symbols_map) const override;
  //! Function to write out the oPowerNode in expr_t terms as opposed to writing out the function itself
  expr_t unpackPowerDeriv() const;
  //! Returns MULT_i*(lhs-rhs) = 0, creating multiplier MULT_i
  expr_t addMultipliersToConstraints(int i);
  //! Returns the non-zero hand-side of an equation (that must have a hand side equal to zero)
  expr_t getNonZeroPartofEquation() const;
  bool isInStaticForm() const override;
  void fillAutoregressiveRow(int eqn, const vector<int> &lhs, map<tuple<int, int, int>, expr_t> &AR) const;
  void findConstantEquations(map<VariableNode *, NumConstNode *> &table) const override;
  expr_t replaceVarsInEquation(map<VariableNode *, NumConstNode *> &table) const override;
  bool containsPacExpectation(const string &pac_model_name = "") const override;
  pair<int, vector<tuple<int, bool, int>>> getPacEC(BinaryOpNode *bopn, int lhs_symb_id, int lhs_orig_symb_id) const;
  void getPacAREC(int lhs_symb_id, int lhs_orig_symb_id,
                  pair<int, vector<tuple<int, bool, int>>> &ec_params_and_vars,
                  set<pair<int, pair<int, int>>> &ar_params_and_vars,
                  vector<tuple<int, int, int, double>> &additive_vars_params_and_constants) const;

  //! Finds the share of optimizing agents in the PAC equation,
  //! the expr node associated with it,
  //! and the expr node associated with the non-optimizing part
  tuple<int, expr_t, expr_t, expr_t> getPacOptimizingShareAndExprNodes(int lhs_symb_id, int lhs_orig_symb_id) const;
  pair<int, expr_t> getPacOptimizingShareAndExprNodesHelper(BinaryOpNode *bopn, int lhs_symb_id, int lhs_orig_symb_id) const;
  expr_t getPacNonOptimizingPart(BinaryOpNode *bopn, int optim_share) const;
  bool getPacNonOptimizingPartHelper(BinaryOpNode *bopn, int optim_share) const;
  bool isParamTimesEndogExpr() const override;
  bool isVarModelReferenced(const string &model_info_name) const override;
  void getEndosAndMaxLags(map<string, int> &model_endos_and_lags) const override;
  //! Substitute auxiliary variables by their expression in static model
  expr_t substituteStaticAuxiliaryVariable() const override;
  //! Substitute auxiliary variables by their expression in static model auxiliary variable definition
  expr_t substituteStaticAuxiliaryDefinition() const;
  void decomposeAdditiveTerms(vector<pair<expr_t, int>> &terms, int current_sign) const override;
};

//! Trinary operator node
class TrinaryOpNode : public ExprNode
{
  friend class ModelTree;
public:
  const expr_t arg1, arg2, arg3;
  const TrinaryOpcode op_code;
private:
  expr_t computeDerivative(int deriv_id) override;
  int cost(int cost, bool is_matlab) const override;
  int cost(const temporary_terms_t &temporary_terms, bool is_matlab) const override;
  int cost(const map<pair<int, int>, temporary_terms_t> &temp_terms_map, bool is_matlab) const override;
  //! Returns the derivative of this node if darg1, darg2 and darg3 are the derivatives of the arguments
  expr_t composeDerivatives(expr_t darg1, expr_t darg2, expr_t darg3);
public:
  TrinaryOpNode(DataTree &datatree_arg, int idx_arg, const expr_t arg1_arg,
                TrinaryOpcode op_code_arg, const expr_t arg2_arg, const expr_t arg3_arg);
  void prepareForDerivation() override;
  int precedence(ExprNodeOutputType output_type, const temporary_terms_t &temporary_terms) const override;
  void computeTemporaryTerms(const pair<int, int> &derivOrder,
                             map<pair<int, int>, temporary_terms_t> &temp_terms_map,
                             map<expr_t, pair<int, pair<int, int>>> &reference_count,
                             bool is_matlab) const override;
  void writeOutput(ostream &output, ExprNodeOutputType output_type, const temporary_terms_t &temporary_terms, const temporary_terms_idxs_t &temporary_terms_idxs, const deriv_node_temp_terms_t &tef_terms) const override;
  void writeJsonAST(ostream &output) const override;
  void writeJsonOutput(ostream &output, const temporary_terms_t &temporary_terms, const deriv_node_temp_terms_t &tef_terms, const bool isdynamic) const override;
  bool containsExternalFunction() const override;
  void writeExternalFunctionOutput(ostream &output, ExprNodeOutputType output_type,
                                           const temporary_terms_t &temporary_terms,
                                           const temporary_terms_idxs_t &temporary_terms_idxs,
                                           deriv_node_temp_terms_t &tef_terms) const override;
  void writeJsonExternalFunctionOutput(vector<string> &efout,
                                               const temporary_terms_t &temporary_terms,
                                               deriv_node_temp_terms_t &tef_terms,
                                               const bool isdynamic) const override;
  void compileExternalFunctionOutput(ostream &CompileCode, unsigned int &instruction_number,
                                             bool lhs_rhs, const temporary_terms_t &temporary_terms,
                                             const map_idx_t &map_idx, bool dynamic, bool steady_dynamic,
                                             deriv_node_temp_terms_t &tef_terms) const override;
  void computeTemporaryTerms(map<expr_t, int> &reference_count,
                                     temporary_terms_t &temporary_terms,
                                     map<expr_t, pair<int, int>> &first_occurence,
                                     int Curr_block,
                                     vector< vector<temporary_terms_t>> &v_temporary_terms,
                                     int equation) const override;
  void collectVARLHSVariable(set<expr_t> &result) const override;
  void collectDynamicVariables(SymbolType type_arg, set<pair<int, int>> &result) const override;
  void collectTemporary_terms(const temporary_terms_t &temporary_terms, temporary_terms_inuse_t &temporary_terms_inuse, int Curr_Block) const override;
  static double eval_opcode(double v1, TrinaryOpcode op_code, double v2, double v3) noexcept(false);
  double eval(const eval_context_t &eval_context) const noexcept(false) override;
  void compile(ostream &CompileCode, unsigned int &instruction_number, bool lhs_rhs, const temporary_terms_t &temporary_terms, const map_idx_t &map_idx, bool dynamic, bool steady_dynamic, const deriv_node_temp_terms_t &tef_terms) const override;
  expr_t toStatic(DataTree &static_datatree) const override;
  void computeXrefs(EquationInfo &ei) const override;
  pair<int, expr_t> normalizeEquation(int symb_id_endo, vector<tuple<int, expr_t, expr_t>>  &List_of_Op_RHS) const override;
  expr_t getChainRuleDerivative(int deriv_id, const map<int, expr_t> &recursive_variables) override;
  int maxEndoLead() const override;
  int maxExoLead() const override;
  int maxEndoLag() const override;
  int maxExoLag() const override;
  int maxLead() const override;
  int maxLag() const override;
  int maxLagWithDiffsExpanded() const override;
  int VarMinLag() const override;
  int VarMaxLag(DataTree &static_datatree, set<expr_t> &static_lhs) const override;
  int PacMaxLag(int lhs_symb_id) const override;
  int getPacTargetSymbId(int lhs_symb_id, int undiff_lhs_symb_id) const override;
  expr_t undiff() const override;
  expr_t decreaseLeadsLags(int n) const override;
  expr_t substituteEndoLeadGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const override;
  //! Creates another TrinaryOpNode with the same opcode, but with a possibly different datatree and arguments
  expr_t buildSimilarTrinaryOpNode(expr_t alt_arg1, expr_t alt_arg2, expr_t alt_arg3, DataTree &alt_datatree) const;
  expr_t substituteEndoLagGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const override;
  expr_t substituteExoLead(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const override;
  expr_t substituteExoLag(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const override;
  expr_t substituteExpectation(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool partial_information_model) const override;
  expr_t substituteAdl() const override;
  expr_t substituteVarExpectation(const map<string, expr_t> &subst_table) const override;
  void findDiffNodes(DataTree &static_datatree, diff_table_t &diff_table) const override;
  void findUnaryOpNodesForAuxVarCreation(DataTree &static_datatree, diff_table_t &nodes) const override;
  int findTargetVariable(int lhs_symb_id) const override;
  expr_t substituteDiff(DataTree &static_datatree, diff_table_t &diff_table, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const override;
  expr_t substituteUnaryOpNodes(DataTree &static_datatree, diff_table_t &nodes, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const override;
  expr_t substitutePacExpectation(const string & name, expr_t subexpr) override;
  expr_t decreaseLeadsLagsPredeterminedVariables() const override;
  expr_t differentiateForwardVars(const vector<string> &subset, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const override;
  bool isNumConstNodeEqualTo(double value) const override;
  bool containsEndogenous() const override;
  bool containsExogenous() const override;
  int countDiffs() const override;
  bool isVariableNodeEqualTo(SymbolType type_arg, int variable_id, int lag_arg) const override;
  expr_t replaceTrendVar() const override;
  expr_t detrend(int symb_id, bool log_trend, expr_t trend) const override;
  expr_t clone(DataTree &datatree) const override;
  expr_t removeTrendLeadLag(map<int, expr_t> trend_symbols_map) const override;
  bool isInStaticForm() const override;
  void findConstantEquations(map<VariableNode *, NumConstNode *> &table) const override;
  expr_t replaceVarsInEquation(map<VariableNode *, NumConstNode *> &table) const override;
  bool containsPacExpectation(const string &pac_model_name = "") const override;
  bool isParamTimesEndogExpr() const override;
  bool isVarModelReferenced(const string &model_info_name) const override;
  void getEndosAndMaxLags(map<string, int> &model_endos_and_lags) const override;
  //! Substitute auxiliary variables by their expression in static model
  expr_t substituteStaticAuxiliaryVariable() const override;
};

//! External function node
class AbstractExternalFunctionNode : public ExprNode
{
public:
  const int symb_id;
  const vector<expr_t> arguments;
private:
  expr_t computeDerivative(int deriv_id) override;
  virtual expr_t composeDerivatives(const vector<expr_t> &dargs) = 0;
protected:
  //! Thrown when trying to access an unknown entry in external_function_node_map
  class UnknownFunctionNameAndArgs
  {
  };
  //! Returns true if the given external function has been written as a temporary term
  bool alreadyWrittenAsTefTerm(int the_symb_id, const deriv_node_temp_terms_t &tef_terms) const;
  //! Returns the index in the tef_terms map of this external function
  int getIndxInTefTerms(int the_symb_id, const deriv_node_temp_terms_t &tef_terms) const noexcept(false);
  //! Helper function to write output arguments of any given external function
  void writeExternalFunctionArguments(ostream &output, ExprNodeOutputType output_type, const temporary_terms_t &temporary_terms, const temporary_terms_idxs_t &temporary_terms_idxs, const deriv_node_temp_terms_t &tef_terms) const;
  void writeJsonASTExternalFunctionArguments(ostream &output) const;
  void writeJsonExternalFunctionArguments(ostream &output, const temporary_terms_t &temporary_terms, const deriv_node_temp_terms_t &tef_terms, const bool isdynamic) const;
  /*! Returns a predicate that tests whether an other ExprNode is an external
    function which is computed by the same external function call (i.e. it has
    the same so-called "Tef" index) */
  virtual function<bool (expr_t)> sameTefTermPredicate() const = 0;
public:
  AbstractExternalFunctionNode(DataTree &datatree_arg, int idx_arg, int symb_id_arg,
                               vector<expr_t> arguments_arg);
  void prepareForDerivation() override;
  void computeTemporaryTerms(const pair<int, int> &derivOrder,
                             map<pair<int, int>, temporary_terms_t> &temp_terms_map,
                             map<expr_t, pair<int, pair<int, int>>> &reference_count,
                             bool is_matlab) const override;
  void writeOutput(ostream &output, ExprNodeOutputType output_type, const temporary_terms_t &temporary_terms, const temporary_terms_idxs_t &temporary_terms_idxs, const deriv_node_temp_terms_t &tef_terms) const override = 0;
  void writeJsonAST(ostream &output) const override = 0;
  void writeJsonOutput(ostream &output, const temporary_terms_t &temporary_terms, const deriv_node_temp_terms_t &tef_terms, const bool isdynamic = true) const override = 0;
  bool containsExternalFunction() const override;
  void writeExternalFunctionOutput(ostream &output, ExprNodeOutputType output_type,
                                           const temporary_terms_t &temporary_terms,
                                           const temporary_terms_idxs_t &temporary_terms_idxs,
                                           deriv_node_temp_terms_t &tef_terms) const override = 0;
  void writeJsonExternalFunctionOutput(vector<string> &efout,
                                               const temporary_terms_t &temporary_terms,
                                               deriv_node_temp_terms_t &tef_terms,
                                               const bool isdynamic = true) const override = 0;
  void compileExternalFunctionOutput(ostream &CompileCode, unsigned int &instruction_number,
                                             bool lhs_rhs, const temporary_terms_t &temporary_terms,
                                             const map_idx_t &map_idx, bool dynamic, bool steady_dynamic,
                                             deriv_node_temp_terms_t &tef_terms) const override = 0;
  void computeTemporaryTerms(map<expr_t, int> &reference_count,
                                     temporary_terms_t &temporary_terms,
                                     map<expr_t, pair<int, int>> &first_occurence,
                                     int Curr_block,
                                     vector< vector<temporary_terms_t>> &v_temporary_terms,
                                     int equation) const override = 0;
  void collectVARLHSVariable(set<expr_t> &result) const override;
  void collectDynamicVariables(SymbolType type_arg, set<pair<int, int>> &result) const override;
  void collectTemporary_terms(const temporary_terms_t &temporary_terms, temporary_terms_inuse_t &temporary_terms_inuse, int Curr_Block) const override;
  double eval(const eval_context_t &eval_context) const noexcept(false) override;
  unsigned int compileExternalFunctionArguments(ostream &CompileCode, unsigned int &instruction_number,
                                                bool lhs_rhs, const temporary_terms_t &temporary_terms,
                                                const map_idx_t &map_idx, bool dynamic, bool steady_dynamic,
                                                const deriv_node_temp_terms_t &tef_terms) const;

  void compile(ostream &CompileCode, unsigned int &instruction_number, bool lhs_rhs, const temporary_terms_t &temporary_terms, const map_idx_t &map_idx, bool dynamic, bool steady_dynamic, const deriv_node_temp_terms_t &tef_terms) const override = 0;
  expr_t toStatic(DataTree &static_datatree) const override = 0;
  void computeXrefs(EquationInfo &ei) const override = 0;
  pair<int, expr_t> normalizeEquation(int symb_id_endo, vector<tuple<int, expr_t, expr_t>>  &List_of_Op_RHS) const override;
  expr_t getChainRuleDerivative(int deriv_id, const map<int, expr_t> &recursive_variables) override;
  int maxEndoLead() const override;
  int maxExoLead() const override;
  int maxEndoLag() const override;
  int maxExoLag() const override;
  int maxLead() const override;
  int maxLag() const override;
  int maxLagWithDiffsExpanded() const override;
  int VarMinLag() const override;
  int VarMaxLag(DataTree &static_datatree, set<expr_t> &static_lhs) const override;
  int PacMaxLag(int lhs_symb_id) const override;
  int getPacTargetSymbId(int lhs_symb_id, int undiff_lhs_symb_id) const override;
  expr_t undiff() const override;
  expr_t decreaseLeadsLags(int n) const override;
  expr_t substituteEndoLeadGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const override;
  expr_t substituteEndoLagGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const override;
  expr_t substituteExoLead(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const override;
  expr_t substituteExoLag(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const override;
  expr_t substituteExpectation(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool partial_information_model) const override;
  expr_t substituteAdl() const override;
  expr_t substituteVarExpectation(const map<string, expr_t> &subst_table) const override;
  void findDiffNodes(DataTree &static_datatree, diff_table_t &diff_table) const override;
  void findUnaryOpNodesForAuxVarCreation(DataTree &static_datatree, diff_table_t &nodes) const override;
  int findTargetVariable(int lhs_symb_id) const override;
  expr_t substituteDiff(DataTree &static_datatree, diff_table_t &diff_table, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const override;
  expr_t substituteUnaryOpNodes(DataTree &static_datatree, diff_table_t &nodes, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const override;
  expr_t substitutePacExpectation(const string & name, expr_t subexpr) override;
  virtual expr_t buildSimilarExternalFunctionNode(vector<expr_t> &alt_args, DataTree &alt_datatree) const = 0;
  expr_t decreaseLeadsLagsPredeterminedVariables() const override;
  expr_t differentiateForwardVars(const vector<string> &subset, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const override;
  bool isNumConstNodeEqualTo(double value) const override;
  bool containsEndogenous() const override;
  bool containsExogenous() const override;
  int countDiffs() const override;
  bool isVariableNodeEqualTo(SymbolType type_arg, int variable_id, int lag_arg) const override;
  void writePrhs(ostream &output, ExprNodeOutputType output_type, const temporary_terms_t &temporary_terms, const temporary_terms_idxs_t &temporary_terms_idxs, const deriv_node_temp_terms_t &tef_terms, const string &ending) const;
  expr_t replaceTrendVar() const override;
  expr_t detrend(int symb_id, bool log_trend, expr_t trend) const override;
  expr_t clone(DataTree &datatree) const override = 0;
  expr_t removeTrendLeadLag(map<int, expr_t> trend_symbols_map) const override;
  bool isInStaticForm() const override;
  void findConstantEquations(map<VariableNode *, NumConstNode *> &table) const override;
  expr_t replaceVarsInEquation(map<VariableNode *, NumConstNode *> &table) const override;
  bool containsPacExpectation(const string &pac_model_name = "") const override;
  bool isParamTimesEndogExpr() const override;
  bool isVarModelReferenced(const string &model_info_name) const override;
  void getEndosAndMaxLags(map<string, int> &model_endos_and_lags) const override;
  //! Substitute auxiliary variables by their expression in static model
  expr_t substituteStaticAuxiliaryVariable() const override;
};

class ExternalFunctionNode : public AbstractExternalFunctionNode
{
  friend class FirstDerivExternalFunctionNode;
  friend class SecondDerivExternalFunctionNode;
private:
  expr_t composeDerivatives(const vector<expr_t> &dargs) override;
protected:
  function<bool (expr_t)> sameTefTermPredicate() const override;
public:
  ExternalFunctionNode(DataTree &datatree_arg, int idx_arg, int symb_id_arg,
                       const vector<expr_t> &arguments_arg);
  void writeOutput(ostream &output, ExprNodeOutputType output_type, const temporary_terms_t &temporary_terms, const temporary_terms_idxs_t &temporary_terms_idxs, const deriv_node_temp_terms_t &tef_terms) const override;
  void writeJsonAST(ostream &output) const override;
  void writeJsonOutput(ostream &output, const temporary_terms_t &temporary_terms, const deriv_node_temp_terms_t &tef_terms, const bool isdynamic) const override;
  void writeExternalFunctionOutput(ostream &output, ExprNodeOutputType output_type,
                                           const temporary_terms_t &temporary_terms,
                                           const temporary_terms_idxs_t &temporary_terms_idxs,
                                           deriv_node_temp_terms_t &tef_terms) const override;
  void writeJsonExternalFunctionOutput(vector<string> &efout,
                                               const temporary_terms_t &temporary_terms,
                                               deriv_node_temp_terms_t &tef_terms,
                                               const bool isdynamic) const override;
  void compileExternalFunctionOutput(ostream &CompileCode, unsigned int &instruction_number,
                                             bool lhs_rhs, const temporary_terms_t &temporary_terms,
                                             const map_idx_t &map_idx, bool dynamic, bool steady_dynamic,
                                             deriv_node_temp_terms_t &tef_terms) const override;
  void computeTemporaryTerms(map<expr_t, int> &reference_count,
                                     temporary_terms_t &temporary_terms,
                                     map<expr_t, pair<int, int>> &first_occurence,
                                     int Curr_block,
                                     vector< vector<temporary_terms_t>> &v_temporary_terms,
                                     int equation) const override;
  void compile(ostream &CompileCode, unsigned int &instruction_number, bool lhs_rhs, const temporary_terms_t &temporary_terms, const map_idx_t &map_idx, bool dynamic, bool steady_dynamic, const deriv_node_temp_terms_t &tef_terms) const override;
  expr_t toStatic(DataTree &static_datatree) const override;
  void computeXrefs(EquationInfo &ei) const override;
  expr_t buildSimilarExternalFunctionNode(vector<expr_t> &alt_args, DataTree &alt_datatree) const override;
  expr_t clone(DataTree &datatree) const override;
};

class FirstDerivExternalFunctionNode : public AbstractExternalFunctionNode
{
public:
  const int inputIndex;
private:
  expr_t composeDerivatives(const vector<expr_t> &dargs) override;
protected:
  function<bool (expr_t)> sameTefTermPredicate() const override;
public:
  FirstDerivExternalFunctionNode(DataTree &datatree_arg, int idx_arg,
                                 int top_level_symb_id_arg,
                                 const vector<expr_t> &arguments_arg,
                                 int inputIndex_arg);
  void computeTemporaryTerms(map<expr_t, int> &reference_count,
                                     temporary_terms_t &temporary_terms,
                                     map<expr_t, pair<int, int>> &first_occurence,
                                     int Curr_block,
                                     vector< vector<temporary_terms_t>> &v_temporary_terms,
                                     int equation) const override;
  void writeOutput(ostream &output, ExprNodeOutputType output_type, const temporary_terms_t &temporary_terms, const temporary_terms_idxs_t &temporary_terms_idxs, const deriv_node_temp_terms_t &tef_terms) const override;
  void writeJsonAST(ostream &output) const override;
  void writeJsonOutput(ostream &output, const temporary_terms_t &temporary_terms, const deriv_node_temp_terms_t &tef_terms, const bool isdynamic) const override;
  void compile(ostream &CompileCode, unsigned int &instruction_number,
                       bool lhs_rhs, const temporary_terms_t &temporary_terms,
                       const map_idx_t &map_idx, bool dynamic, bool steady_dynamic,
                       const deriv_node_temp_terms_t &tef_terms) const override;
  void writeExternalFunctionOutput(ostream &output, ExprNodeOutputType output_type,
                                           const temporary_terms_t &temporary_terms,
                                           const temporary_terms_idxs_t &temporary_terms_idxs,
                                           deriv_node_temp_terms_t &tef_terms) const override;
  void writeJsonExternalFunctionOutput(vector<string> &efout,
                                               const temporary_terms_t &temporary_terms,
                                               deriv_node_temp_terms_t &tef_terms,
                                               const bool isdynamic) const override;
  void compileExternalFunctionOutput(ostream &CompileCode, unsigned int &instruction_number,
                                             bool lhs_rhs, const temporary_terms_t &temporary_terms,
                                             const map_idx_t &map_idx, bool dynamic, bool steady_dynamic,
                                             deriv_node_temp_terms_t &tef_terms) const override;
  expr_t toStatic(DataTree &static_datatree) const override;
  void computeXrefs(EquationInfo &ei) const override;
  expr_t buildSimilarExternalFunctionNode(vector<expr_t> &alt_args, DataTree &alt_datatree) const override;
  expr_t clone(DataTree &datatree) const override;
};

class SecondDerivExternalFunctionNode : public AbstractExternalFunctionNode
{
public:
  const int inputIndex1;
  const int inputIndex2;
private:
  expr_t composeDerivatives(const vector<expr_t> &dargs) override;
protected:
  function<bool (expr_t)> sameTefTermPredicate() const override;
public:
  SecondDerivExternalFunctionNode(DataTree &datatree_arg, int idx_arg,
                                  int top_level_symb_id_arg,
                                  const vector<expr_t> &arguments_arg,
                                  int inputIndex1_arg,
                                  int inputIndex2_arg);
  void computeTemporaryTerms(map<expr_t, int> &reference_count,
                                     temporary_terms_t &temporary_terms,
                                     map<expr_t, pair<int, int>> &first_occurence,
                                     int Curr_block,
                                     vector< vector<temporary_terms_t>> &v_temporary_terms,
                                     int equation) const override;
  void writeOutput(ostream &output, ExprNodeOutputType output_type, const temporary_terms_t &temporary_terms, const temporary_terms_idxs_t &temporary_terms_idxs, const deriv_node_temp_terms_t &tef_terms) const override;
  void writeJsonAST(ostream &output) const override;
  void writeJsonOutput(ostream &output, const temporary_terms_t &temporary_terms, const deriv_node_temp_terms_t &tef_terms, const bool isdynamic) const override;
  void compile(ostream &CompileCode, unsigned int &instruction_number,
                       bool lhs_rhs, const temporary_terms_t &temporary_terms,
                       const map_idx_t &map_idx, bool dynamic, bool steady_dynamic,
                       const deriv_node_temp_terms_t &tef_terms) const override;
  void writeExternalFunctionOutput(ostream &output, ExprNodeOutputType output_type,
                                           const temporary_terms_t &temporary_terms,
                                           const temporary_terms_idxs_t &temporary_terms_idxs,
                                           deriv_node_temp_terms_t &tef_terms) const override;
  void writeJsonExternalFunctionOutput(vector<string> &efout,
                                               const temporary_terms_t &temporary_terms,
                                               deriv_node_temp_terms_t &tef_terms,
                                               const bool isdynamic) const override;
  void compileExternalFunctionOutput(ostream &CompileCode, unsigned int &instruction_number,
                                             bool lhs_rhs, const temporary_terms_t &temporary_terms,
                                             const map_idx_t &map_idx, bool dynamic, bool steady_dynamic,
                                             deriv_node_temp_terms_t &tef_terms) const override;
  expr_t toStatic(DataTree &static_datatree) const override;
  void computeXrefs(EquationInfo &ei) const override;
  expr_t buildSimilarExternalFunctionNode(vector<expr_t> &alt_args, DataTree &alt_datatree) const override;
  expr_t clone(DataTree &datatree) const override;
};

class VarExpectationNode : public ExprNode
{
public:
  const string model_name;
  VarExpectationNode(DataTree &datatree_arg, int idx_arg, string model_name_arg);
  void computeTemporaryTerms(const pair<int, int> &derivOrder,
                             map<pair<int, int>, temporary_terms_t> &temp_terms_map,
                             map<expr_t, pair<int, pair<int, int>>> &reference_count,
                             bool is_matlab) const override;
  void writeOutput(ostream &output, ExprNodeOutputType output_type, const temporary_terms_t &temporary_terms, const temporary_terms_idxs_t &temporary_terms_idxs, const deriv_node_temp_terms_t &tef_terms) const override;
  void computeTemporaryTerms(map<expr_t, int> &reference_count,
                                     temporary_terms_t &temporary_terms,
                                     map<expr_t, pair<int, int>> &first_occurence,
                                     int Curr_block,
                                     vector< vector<temporary_terms_t>> &v_temporary_terms,
                                     int equation) const override;
  expr_t toStatic(DataTree &static_datatree) const override;
  expr_t clone(DataTree &datatree) const override;
  int maxEndoLead() const override;
  int maxExoLead() const override;
  int maxEndoLag() const override;
  int maxExoLag() const override;
  int maxLead() const override;
  int maxLag() const override;
  int maxLagWithDiffsExpanded() const override;
  int VarMinLag() const override;
  int VarMaxLag(DataTree &static_datatree, set<expr_t> &static_lhs) const override;
  int PacMaxLag(int lhs_symb_id) const override;
  int getPacTargetSymbId(int lhs_symb_id, int undiff_lhs_symb_id) const override;
  expr_t undiff() const override;
  expr_t decreaseLeadsLags(int n) const override;
  void prepareForDerivation() override;
  expr_t computeDerivative(int deriv_id) override;
  expr_t getChainRuleDerivative(int deriv_id, const map<int, expr_t> &recursive_variables) override;
  bool containsExternalFunction() const override;
  double eval(const eval_context_t &eval_context) const noexcept(false) override;
  void computeXrefs(EquationInfo &ei) const override;
  expr_t substituteEndoLeadGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const override;
  expr_t substituteEndoLagGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const override;
  expr_t substituteExoLead(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const override;
  expr_t substituteExoLag(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const override;
  expr_t substituteExpectation(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool partial_information_model) const override;
  expr_t substituteAdl() const override;
  expr_t substituteVarExpectation(const map<string, expr_t> &subst_table) const override;
  void findDiffNodes(DataTree &static_datatree, diff_table_t &diff_table) const override;
  void findUnaryOpNodesForAuxVarCreation(DataTree &static_datatree, diff_table_t &nodes) const override;
  int findTargetVariable(int lhs_symb_id) const override;
  expr_t substituteDiff(DataTree &static_datatree, diff_table_t &diff_table, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const override;
  expr_t substituteUnaryOpNodes(DataTree &static_datatree, diff_table_t &nodes, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const override;
  expr_t substitutePacExpectation(const string & name, expr_t subexpr) override;
  pair<int, expr_t> normalizeEquation(int symb_id_endo, vector<tuple<int, expr_t, expr_t>>  &List_of_Op_RHS) const override;
  void compile(ostream &CompileCode, unsigned int &instruction_number,
                       bool lhs_rhs, const temporary_terms_t &temporary_terms,
                       const map_idx_t &map_idx, bool dynamic, bool steady_dynamic,
                       const deriv_node_temp_terms_t &tef_terms) const override;
  void collectTemporary_terms(const temporary_terms_t &temporary_terms, temporary_terms_inuse_t &temporary_terms_inuse, int Curr_Block) const override;
  void collectVARLHSVariable(set<expr_t> &result) const override;
  void collectDynamicVariables(SymbolType type_arg, set<pair<int, int>> &result) const override;
  bool containsEndogenous() const override;
  bool containsExogenous() const override;
  int countDiffs() const override;
  bool isNumConstNodeEqualTo(double value) const override;
  expr_t differentiateForwardVars(const vector<string> &subset, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const override;
  expr_t decreaseLeadsLagsPredeterminedVariables() const override;
  bool isVariableNodeEqualTo(SymbolType type_arg, int variable_id, int lag_arg) const override;
  expr_t replaceTrendVar() const override;
  expr_t detrend(int symb_id, bool log_trend, expr_t trend) const override;
  expr_t removeTrendLeadLag(map<int, expr_t> trend_symbols_map) const override;
  bool isInStaticForm() const override;
  void findConstantEquations(map<VariableNode *, NumConstNode *> &table) const override;
  expr_t replaceVarsInEquation(map<VariableNode *, NumConstNode *> &table) const override;
  bool containsPacExpectation(const string &pac_model_name = "") const override;
  bool isParamTimesEndogExpr() const override;
  bool isVarModelReferenced(const string &model_info_name) const override;
  void getEndosAndMaxLags(map<string, int> &model_endos_and_lags) const override;
  expr_t substituteStaticAuxiliaryVariable() const override;
  void writeJsonAST(ostream &output) const override;
  void writeJsonOutput(ostream &output, const temporary_terms_t &temporary_terms, const deriv_node_temp_terms_t &tef_terms, const bool isdynamic) const override;
};

class PacExpectationNode : public ExprNode
{
public:
  const string model_name;
public:
  PacExpectationNode(DataTree &datatree_arg, int idx_arg, string model_name);
  void computeTemporaryTerms(const pair<int, int> &derivOrder,
                             map<pair<int, int>, temporary_terms_t> &temp_terms_map,
                             map<expr_t, pair<int, pair<int, int>>> &reference_count,
                             bool is_matlab) const override;
  void writeOutput(ostream &output, ExprNodeOutputType output_type, const temporary_terms_t &temporary_terms, const temporary_terms_idxs_t &temporary_terms_idxs, const deriv_node_temp_terms_t &tef_terms) const override;
  void computeTemporaryTerms(map<expr_t, int> &reference_count,
                                     temporary_terms_t &temporary_terms,
                                     map<expr_t, pair<int, int>> &first_occurence,
                                     int Curr_block,
                                     vector< vector<temporary_terms_t>> &v_temporary_terms,
                                     int equation) const override;
  expr_t toStatic(DataTree &static_datatree) const override;
  expr_t clone(DataTree &datatree) const override;
  int maxEndoLead() const override;
  int maxExoLead() const override;
  int maxEndoLag() const override;
  int maxExoLag() const override;
  int maxLead() const override;
  int maxLag() const override;
  int maxLagWithDiffsExpanded() const override;
  int VarMinLag() const override;
  int VarMaxLag(DataTree &static_datatree, set<expr_t> &static_lhs) const override;
  int PacMaxLag(int lhs_symb_id) const override;
  int getPacTargetSymbId(int lhs_symb_id, int undiff_lhs_symb_id) const override;
  expr_t undiff() const override;
  expr_t decreaseLeadsLags(int n) const override;
  void prepareForDerivation() override;
  expr_t computeDerivative(int deriv_id) override;
  expr_t getChainRuleDerivative(int deriv_id, const map<int, expr_t> &recursive_variables) override;
  bool containsExternalFunction() const override;
  double eval(const eval_context_t &eval_context) const noexcept(false) override;
  void computeXrefs(EquationInfo &ei) const override;
  expr_t substituteEndoLeadGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const override;
  expr_t substituteEndoLagGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const override;
  expr_t substituteExoLead(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const override;
  expr_t substituteExoLag(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const override;
  expr_t substituteExpectation(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool partial_information_model) const override;
  expr_t substituteAdl() const override;
  expr_t substituteVarExpectation(const map<string, expr_t> &subst_table) const override;
  void findDiffNodes(DataTree &static_datatree, diff_table_t &diff_table) const override;
  void findUnaryOpNodesForAuxVarCreation(DataTree &static_datatree, diff_table_t &nodes) const override;
  int findTargetVariable(int lhs_symb_id) const override;
  expr_t substituteDiff(DataTree &static_datatree, diff_table_t &diff_table, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const override;
  expr_t substituteUnaryOpNodes(DataTree &static_datatree, diff_table_t &nodes, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const override;
  expr_t substitutePacExpectation(const string & name, expr_t subexpr) override;
  pair<int, expr_t> normalizeEquation(int symb_id_endo, vector<tuple<int, expr_t, expr_t>>  &List_of_Op_RHS) const override;
  void compile(ostream &CompileCode, unsigned int &instruction_number,
                       bool lhs_rhs, const temporary_terms_t &temporary_terms,
                       const map_idx_t &map_idx, bool dynamic, bool steady_dynamic,
                       const deriv_node_temp_terms_t &tef_terms) const override;
  void collectTemporary_terms(const temporary_terms_t &temporary_terms, temporary_terms_inuse_t &temporary_terms_inuse, int Curr_Block) const override;
  void collectVARLHSVariable(set<expr_t> &result) const override;
  void collectDynamicVariables(SymbolType type_arg, set<pair<int, int>> &result) const override;
  bool containsEndogenous() const override;
  bool containsExogenous() const override;
  int countDiffs() const override;
  bool isNumConstNodeEqualTo(double value) const override;
  expr_t differentiateForwardVars(const vector<string> &subset, subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const override;
  expr_t decreaseLeadsLagsPredeterminedVariables() const override;
  bool isVariableNodeEqualTo(SymbolType type_arg, int variable_id, int lag_arg) const override;
  expr_t replaceTrendVar() const override;
  expr_t detrend(int symb_id, bool log_trend, expr_t trend) const override;
  expr_t removeTrendLeadLag(map<int, expr_t> trend_symbols_map) const override;
  bool isInStaticForm() const override;
  void findConstantEquations(map<VariableNode *, NumConstNode *> &table) const override;
  expr_t replaceVarsInEquation(map<VariableNode *, NumConstNode *> &table) const override;
  bool containsPacExpectation(const string &pac_model_name = "") const override;
  bool isParamTimesEndogExpr() const override;
  bool isVarModelReferenced(const string &model_info_name) const override;
  void getEndosAndMaxLags(map<string, int> &model_endos_and_lags) const override;
  expr_t substituteStaticAuxiliaryVariable() const override;
  void writeJsonAST(ostream &output) const override;
  void writeJsonOutput(ostream &output, const temporary_terms_t &temporary_terms, const deriv_node_temp_terms_t &tef_terms, const bool isdynamic) const override;
};

#endif
