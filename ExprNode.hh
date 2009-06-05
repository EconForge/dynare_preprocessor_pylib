/*
 * Copyright (C) 2007-2009 Dynare Team
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

using namespace std;

#include <set>
#include <map>
#include <vector>
#include <iostream>
#include <fstream>

#include "SymbolTable.hh"
#include "CodeInterpreter.hh"

class DataTree;

typedef class ExprNode *NodeID;

struct Model_Block;

struct ExprNodeLess;

//! Type for set of temporary terms
/*! They are ordered by index number thanks to ExprNodeLess */
typedef set<NodeID, ExprNodeLess> temporary_terms_type;

typedef map<int,int> map_idx_type;
typedef set<int> temporary_terms_inuse_type;

//! Type for evaluation contexts
/*! The key is a symbol id. Lags are assumed to be null */
typedef map<int, double> eval_context_type;

//! Possible types of output when writing ExprNode(s)
enum ExprNodeOutputType
  {
    oMatlabStaticModel,       //!< Matlab code, static model declarations
    oMatlabDynamicModel,      //!< Matlab code, dynamic model declarations
    oMatlabStaticModelSparse, //!< Matlab code, static block decomposed mode declaration
    oMatlabDynamicModelSparse,//!< Matlab code, dynamic block decomposed mode declaration
    oCStaticModel,            //!< C code, static model declarations
    oCDynamicModel,           //!< C code, dynamic model declarations
    oMatlabOutsideModel,      //!< Matlab code, outside model block (for example in initval)
    oLatexStaticModel,        //!< LaTeX code, static model declarations
    oLatexDynamicModel        //!< LaTeX code, dynamic model declarations
  };

#define IS_MATLAB(output_type) ((output_type) == oMatlabStaticModel     \
                                || (output_type) == oMatlabDynamicModel \
                                || (output_type) == oMatlabOutsideModel \
                                || (output_type) == oMatlabStaticModelSparse \
                                || (output_type) == oMatlabDynamicModelSparse)

#define IS_C(output_type) ((output_type) == oCStaticModel       \
                           || (output_type) == oCDynamicModel)

#define IS_LATEX(output_type) ((output_type) == oLatexStaticModel       \
                               || (output_type) == oLatexDynamicModel)

/* Equal to 1 for Matlab langage, or to 0 for C language. Not defined for LaTeX.
   In Matlab, array indexes begin at 1, while they begin at 0 in C */
#define ARRAY_SUBSCRIPT_OFFSET(output_type) ((int) IS_MATLAB(output_type))

// Left and right array subscript delimiters: '(' and ')' for Matlab, '[' and ']' for C
#define LEFT_ARRAY_SUBSCRIPT(output_type) (IS_MATLAB(output_type) ? '(' : '[')
#define RIGHT_ARRAY_SUBSCRIPT(output_type) (IS_MATLAB(output_type) ? ')' : ']')

// Left and right parentheses
#define LEFT_PAR(output_type) (IS_LATEX(output_type) ? "\\left(" : "(")
#define RIGHT_PAR(output_type) (IS_LATEX(output_type) ? "\\right)" : ")")

// Computing cost above which a node can be declared a temporary term
#define MIN_COST_MATLAB (40*90)
#define MIN_COST_C (40*4)
#define MIN_COST(is_matlab) ((is_matlab) ? MIN_COST_MATLAB : MIN_COST_C)

//! Base class for expression nodes
class ExprNode
{
  friend class DataTree;
  friend class DynamicModel;
  friend class ExprNodeLess;
  friend class NumConstNode;
  friend class VariableNode;
  friend class UnaryOpNode;
  friend class BinaryOpNode;
  friend class TrinaryOpNode;

private:
  //! Computes derivative w.r. to a derivation ID (but doesn't store it in derivatives map)
  /*! You shoud use getDerivative() to get the benefit of symbolic a priori and of caching */
  virtual NodeID computeDerivative(int deriv_id) = 0;
  //! Computes derivative w.r. to a derivation ID and use chaine rule derivatives (but doesn't store it in derivatives map)
  /*! You shoud use getDerivative() to get the benefit of symbolic a priori and of caching */
  virtual NodeID computeChaineRuleDerivative(int deriv_id, map<int, NodeID> &recursive_variables, int var, int lag_) = 0;


protected:
  //! Reference to the enclosing DataTree
  DataTree &datatree;

  //! Index number
  int idx;

  //! Set of derivation IDs with respect to which the derivative is potentially non-null
  set<int> non_null_derivatives;

  //! Used for caching of first order derivatives (when non-null)
  map<int, NodeID> derivatives;

  //! Cost of computing current node
  /*! Nodes included in temporary_terms are considered having a null cost */
  virtual int cost(const temporary_terms_type &temporary_terms, bool is_matlab) const;

public:
  ExprNode(DataTree &datatree_arg);
  virtual ~ExprNode();

  //! Returns derivative w.r. to derivation ID
  /*! Uses a symbolic a priori to pre-detect null derivatives, and caches the result for other derivatives (to avoid computing it several times)
    For an equal node, returns the derivative of lhs minus rhs */
  NodeID getDerivative(int deriv_id);

  //! Returns derivative w.r. to derivation ID and use if it possible chaine rule derivatives
  NodeID getChaineRuleDerivative(int deriv_id, map<int, NodeID> &recursive_variables, int var, int lag_);


  //! Returns precedence of node
  /*! Equals 100 for constants, variables, unary ops, and temporary terms */
  virtual int precedence(ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms) const;

  //! Fills temporary_terms set, using reference counts
  /*! A node will be marked as a temporary term if it is referenced at least two times (i.e. has at least two parents), and has a computing cost (multiplied by reference count) greater to datatree.min_cost */
  virtual void computeTemporaryTerms(map<NodeID, int> &reference_count, temporary_terms_type &temporary_terms, bool is_matlab) const;

  //! Writes output of node, using a Txxx notation for nodes in temporary_terms
  virtual void writeOutput(ostream &output, ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms) const /*= 0*/;

  //! Writes output of node (with no temporary terms and with "outside model" output type)
  void writeOutput(ostream &output);

  //! Computes the set of endogenous variables in the expression
  /*!
    Endogenous are stored as integer pairs of the form (symb_id, lag).
    They are added to the set given in argument.
  */
  virtual void collectEndogenous(set<pair<int, int> > &result) const = 0;

  //! Computes the set of exogenous variables in the expression
  /*!
    Exogenous are stored as integer pairs of the form (symb_id, lag).
    They are added to the set given in argument.
  */
  virtual void collectExogenous(set<pair<int, int> > &result) const = 0;

  virtual void collectTemporary_terms(const temporary_terms_type &temporary_terms, Model_Block *ModelBlock, int Curr_Block) const = 0;
  virtual void computeTemporaryTerms(map<NodeID, int> &reference_count,
                                     temporary_terms_type &temporary_terms,
                                     map<NodeID, pair<int, int> > &first_occurence,
                                     int Curr_block,
                                     Model_Block *ModelBlock,
                                     int equation,
                                     map_idx_type &map_idx) const;

  class EvalException

  {
  };

  virtual double eval(const eval_context_type &eval_context) const throw (EvalException) = 0;
  virtual void compile(ofstream &CompileCode, bool lhs_rhs, const temporary_terms_type &temporary_terms, map_idx_type &map_idx) const = 0;
  //! Creates a static version of this node
  /*!
    This method duplicates the current node by creating a similar node from which all leads/lags have been stripped,
    adds the result in the static_datatree argument (and not in the original datatree), and returns it.
  */
  virtual NodeID toStatic(DataTree &static_datatree) const = 0;
  //! Try to normalize an equation linear in its endogenous variable
  virtual pair<bool, NodeID> normalizeLinearInEndoEquation(int symb_id_endo, NodeID Derivative) const;
};

//! Object used to compare two nodes (using their indexes)
struct ExprNodeLess
{
  bool operator()(NodeID arg1, NodeID arg2) const
  {
    return arg1->idx < arg2->idx;
  }
};

//! Numerical constant node
/*! The constant is necessarily non-negative (this is enforced at the NumericalConstants class level) */
class NumConstNode : public ExprNode
{
private:
  //! Id from numerical constants table
  const int id;
  virtual NodeID computeDerivative(int deriv_id);
  virtual NodeID computeChaineRuleDerivative(int deriv_id, map<int, NodeID> &recursive_variables, int var, int lag_);
public:
  NumConstNode(DataTree &datatree_arg, int id_arg);
  virtual void writeOutput(ostream &output, ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms) const;
  virtual void collectEndogenous(set<pair<int, int> > &result) const;
  virtual void collectExogenous(set<pair<int, int> > &result) const;
  virtual void collectTemporary_terms(const temporary_terms_type &temporary_terms, Model_Block *ModelBlock, int Curr_Block) const;
  virtual double eval(const eval_context_type &eval_context) const throw (EvalException);
  virtual void compile(ofstream &CompileCode, bool lhs_rhs, const temporary_terms_type &temporary_terms, map_idx_type &map_idx) const;
  virtual NodeID toStatic(DataTree &static_datatree) const;
  virtual pair<bool, NodeID> normalizeLinearInEndoEquation(int symb_id_endo, NodeID Derivative) const;
};

//! Symbol or variable node
class VariableNode : public ExprNode
{
private:
  //! Id from the symbol table
  const int symb_id;
  const SymbolType type;
  const int lag;
  //! Derivation ID
  const int deriv_id;
  virtual NodeID computeDerivative(int deriv_id_arg);
  virtual NodeID computeChaineRuleDerivative(int deriv_id, map<int, NodeID> &recursive_variables, int var, int lag_);
public:
  VariableNode(DataTree &datatree_arg, int symb_id_arg, int lag_arg, int deriv_id_arg);
  virtual void writeOutput(ostream &output, ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms = temporary_terms_type()) const;
  virtual void collectEndogenous(set<pair<int, int> > &result) const;
  virtual void collectExogenous(set<pair<int, int> > &result) const;
  virtual void computeTemporaryTerms(map<NodeID, int> &reference_count,
                                     temporary_terms_type &temporary_terms,
                                     map<NodeID, pair<int, int> > &first_occurence,
                                     int Curr_block,
                                     Model_Block *ModelBlock,
                                     int equation,
                                     map_idx_type &map_idx) const;
  virtual void collectTemporary_terms(const temporary_terms_type &temporary_terms, Model_Block *ModelBlock, int Curr_Block) const;
  virtual double eval(const eval_context_type &eval_context) const throw (EvalException);
  virtual void compile(ofstream &CompileCode, bool lhs_rhs, const temporary_terms_type &temporary_terms, map_idx_type &map_idx) const;
  virtual NodeID toStatic(DataTree &static_datatree) const;
  int get_symb_id() const { return symb_id; };
  virtual pair<bool, NodeID> normalizeLinearInEndoEquation(int symb_id_endo, NodeID Derivative) const;
};

//! Unary operator node
class UnaryOpNode : public ExprNode
{
private:
  const NodeID arg;
  const UnaryOpcode op_code;
  virtual NodeID computeDerivative(int deriv_id);
  virtual NodeID computeChaineRuleDerivative(int deriv_id, map<int, NodeID> &recursive_variables, int var, int lag_);

  virtual int cost(const temporary_terms_type &temporary_terms, bool is_matlab) const;
public:
  UnaryOpNode(DataTree &datatree_arg, UnaryOpcode op_code_arg, const NodeID arg_arg);
  virtual void computeTemporaryTerms(map<NodeID, int> &reference_count, temporary_terms_type &temporary_terms, bool is_matlab) const;
  virtual void writeOutput(ostream &output, ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms) const;
  virtual void computeTemporaryTerms(map<NodeID, int> &reference_count,
                                     temporary_terms_type &temporary_terms,
                                     map<NodeID, pair<int, int> > &first_occurence,
                                     int Curr_block,
                                     Model_Block *ModelBlock,
                                     int equation,
                                     map_idx_type &map_idx) const;
  virtual void collectEndogenous(set<pair<int, int> > &result) const;
  virtual void collectExogenous(set<pair<int, int> > &result) const;
  virtual void collectTemporary_terms(const temporary_terms_type &temporary_terms, Model_Block *ModelBlock, int Curr_Block) const;
  static double eval_opcode(UnaryOpcode op_code, double v) throw (EvalException);
  virtual double eval(const eval_context_type &eval_context) const throw (EvalException);
  virtual void compile(ofstream &CompileCode, bool lhs_rhs, const temporary_terms_type &temporary_terms, map_idx_type &map_idx) const;
  //! Returns operand
  NodeID get_arg() const { return(arg); };
  //! Returns op code
  UnaryOpcode get_op_code() const { return(op_code); };
  virtual NodeID toStatic(DataTree &static_datatree) const;
  virtual pair<bool, NodeID> normalizeLinearInEndoEquation(int symb_id_endo, NodeID Derivative) const;
};

//! Binary operator node
class BinaryOpNode : public ExprNode
{
private:
  const NodeID arg1, arg2;
  const BinaryOpcode op_code;
  virtual NodeID computeDerivative(int deriv_id);
  virtual NodeID computeChaineRuleDerivative(int deriv_id, map<int, NodeID> &recursive_variables, int var, int lag_);

  virtual int cost(const temporary_terms_type &temporary_terms, bool is_matlab) const;
public:
  BinaryOpNode(DataTree &datatree_arg, const NodeID arg1_arg,
               BinaryOpcode op_code_arg, const NodeID arg2_arg);
  virtual int precedence(ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms) const;
  virtual void computeTemporaryTerms(map<NodeID, int> &reference_count, temporary_terms_type &temporary_terms, bool is_matlab) const;
  virtual void writeOutput(ostream &output, ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms) const;
  virtual void computeTemporaryTerms(map<NodeID, int> &reference_count,
                                     temporary_terms_type &temporary_terms,
                                     map<NodeID, pair<int, int> > &first_occurence,
                                     int Curr_block,
                                     Model_Block *ModelBlock,
                                     int equation,
                                     map_idx_type &map_idx) const;
  virtual void collectEndogenous(set<pair<int, int> > &result) const;
  virtual void collectExogenous(set<pair<int, int> > &result) const;
  virtual void collectTemporary_terms(const temporary_terms_type &temporary_terms, Model_Block *ModelBlock, int Curr_Block) const;
  static double eval_opcode(double v1, BinaryOpcode op_code, double v2) throw (EvalException);
  virtual double eval(const eval_context_type &eval_context) const throw (EvalException);
  virtual void compile(ofstream &CompileCode, bool lhs_rhs, const temporary_terms_type &temporary_terms, map_idx_type &map_idx) const;
  //! Returns first operand
  NodeID get_arg1() const { return(arg1); };
  //! Returns second operand
  NodeID get_arg2() const { return(arg2); };
  //! Returns op code
  BinaryOpcode get_op_code() const { return(op_code); };
  virtual NodeID toStatic(DataTree &static_datatree) const;
  pair<bool, NodeID> normalizeLinearInEndoEquation(int symb_id_endo, NodeID Derivative) const;
};

//! Trinary operator node
class TrinaryOpNode : public ExprNode
{
  friend class ModelTree;
private:
  const NodeID arg1, arg2, arg3;
  const TrinaryOpcode op_code;
  virtual NodeID computeDerivative(int deriv_id);
  virtual NodeID computeChaineRuleDerivative(int deriv_id, map<int, NodeID> &recursive_variables, int var, int lag_);

  virtual int cost(const temporary_terms_type &temporary_terms, bool is_matlab) const;
public:
  TrinaryOpNode(DataTree &datatree_arg, const NodeID arg1_arg,
		TrinaryOpcode op_code_arg, const NodeID arg2_arg, const NodeID arg3_arg);
  virtual int precedence(ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms) const;
  virtual void computeTemporaryTerms(map<NodeID, int> &reference_count, temporary_terms_type &temporary_terms, bool is_matlab) const;
  virtual void writeOutput(ostream &output, ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms) const;
  virtual void computeTemporaryTerms(map<NodeID, int> &reference_count,
                                     temporary_terms_type &temporary_terms,
                                     map<NodeID, pair<int, int> > &first_occurence,
                                     int Curr_block,
                                     Model_Block *ModelBlock,
                                     int equation,
                                     map_idx_type &map_idx) const;
  virtual void collectEndogenous(set<pair<int, int> > &result) const;
  virtual void collectExogenous(set<pair<int, int> > &result) const;
  virtual void collectTemporary_terms(const temporary_terms_type &temporary_terms, Model_Block *ModelBlock, int Curr_Block) const;
  static double eval_opcode(double v1, TrinaryOpcode op_code, double v2, double v3) throw (EvalException);
  virtual double eval(const eval_context_type &eval_context) const throw (EvalException);
  virtual void compile(ofstream &CompileCode, bool lhs_rhs, const temporary_terms_type &temporary_terms, map_idx_type &map_idx) const;
  virtual NodeID toStatic(DataTree &static_datatree) const;
  virtual pair<bool, NodeID> normalizeLinearInEndoEquation(int symb_id_endo, NodeID Derivative) const;
};

//! Unknown function node
class UnknownFunctionNode : public ExprNode
{
private:
  const int symb_id;
  const vector<NodeID> arguments;
  virtual NodeID computeDerivative(int deriv_id);
  virtual NodeID computeChaineRuleDerivative(int deriv_id, map<int, NodeID> &recursive_variables, int var, int lag_);
public:
  UnknownFunctionNode(DataTree &datatree_arg, int symb_id_arg,
                      const vector<NodeID> &arguments_arg);
  virtual void computeTemporaryTerms(map<NodeID, int> &reference_count, temporary_terms_type &temporary_terms, bool is_matlab) const;
  virtual void writeOutput(ostream &output, ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms) const;
  virtual void computeTemporaryTerms(map<NodeID, int> &reference_count,
                                     temporary_terms_type &temporary_terms,
                                     map<NodeID, pair<int, int> > &first_occurence,
                                     int Curr_block,
                                     Model_Block *ModelBlock,
                                     int equation,
                                     map_idx_type &map_idx) const;
  virtual void collectEndogenous(set<pair<int, int> > &result) const;
  virtual void collectExogenous(set<pair<int, int> > &result) const;
  virtual void collectTemporary_terms(const temporary_terms_type &temporary_terms, Model_Block *ModelBlock, int Curr_Block) const;
  virtual double eval(const eval_context_type &eval_context) const throw (EvalException);
  virtual void compile(ofstream &CompileCode, bool lhs_rhs, const temporary_terms_type &temporary_terms, map_idx_type &map_idx) const;
  virtual NodeID toStatic(DataTree &static_datatree) const;
  virtual pair<bool, NodeID> normalizeLinearInEndoEquation(int symb_id_endo, NodeID Derivative) const;
};

//! For one lead/lag of one block, stores mapping of information between original model and block-decomposed model
struct IM_compact
{
  int size, u_init, u_finish, nb_endo, nb_other_endo, size_exo, size_other_endo;
  int *u, *us, *Var, *Equ, *Var_Index, *Equ_Index, *Exogenous, *Exogenous_Index, *Equ_X, *Equ_X_Index;
  int *u_other_endo, *Var_other_endo, *Equ_other_endo, *Var_Index_other_endo, *Equ_Index_other_endo;
};

//! One block of the model
struct Block
{
  int Size, Sized, nb_exo, nb_exo_det, nb_other_endo, Nb_Recursives;
  BlockType Type;
  BlockSimulationType Simulation_Type;
  int Max_Lead, Max_Lag, Nb_Lead_Lag_Endo;
  int Max_Lag_Endo, Max_Lead_Endo;
  int Max_Lag_Other_Endo, Max_Lead_Other_Endo;
  int Max_Lag_Exo, Max_Lead_Exo;
  bool is_linear;
  int *Equation, *Own_Derivative;
  EquationType *Equation_Type;
  NodeID *Equation_Normalized;
  int *Variable, *Other_Endogenous, *Exogenous;
  temporary_terms_type **Temporary_Terms_in_Equation;
  //temporary_terms_type *Temporary_terms;
  temporary_terms_inuse_type *Temporary_InUse;
  IM_compact *IM_lead_lag;
  int Code_Start, Code_Length;
};



//! The set of all blocks of the model
struct Model_Block
{
  int Size, Periods;
  Block* Block_List;
  //int *in_Block_Equ, *in_Block_Var, *in_Equ_of_Block, *in_Var_of_Block;
};

#endif
