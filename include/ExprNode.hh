/*
 * Copyright (C) 2007-2008 Dynare Team
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

typedef struct Model_Block;

struct ExprNodeLess;

//! Type for set of temporary terms
/*! They are ordered by index number thanks to ExprNodeLess */
typedef set<NodeID, ExprNodeLess> temporary_terms_type;
typedef map<int,int> map_idx_type;

//! Possible types of output when writing ExprNode(s)
enum ExprNodeOutputType
  {
    oMatlabStaticModel,       //!< Matlab code, static model declarations
    oMatlabDynamicModel,      //!< Matlab code, dynamic model declarations
    oMatlabStaticModelSparse, //!< Matlab code, static block decomposed mode declaration
    oMatlabDynamicModelSparse, //!< Matlab code, dynamic block decomposed mode declaration
    oCStaticModel,            //!< C code, static model declarations
    oCDynamicModel,           //!< C code, dynamic model declarations
    oCDynamicModelSparseDLL,  //!< C code, dynamic model declarations in SparseDLL module
    oMatlabOutsideModel       //!< Matlab code, outside model block (for example in initval)
  };

//! Type for evaluation contexts
/*! The key is a pair (symbol id, symbol type)
  Lags are assumed to be null */
typedef map<pair<int, Type>, double> eval_context_type;

/* Equal to 1 for Matlab langage, or to 0 for C language
   In Matlab, array indexes begin at 1, while they begin at 0 in C */
#define OFFSET(output_type) ((output_type == oMatlabStaticModel)      \
                             || (output_type == oMatlabDynamicModel)  \
                             || (output_type == oMatlabOutsideModel)  \
                             || (output_type == oMatlabStaticModelSparse)  \
                             || (output_type == oMatlabDynamicModelSparse))

// Left parenthesis: '(' for Matlab, '[' for C
#define LPAR(output_type) (OFFSET(output_type) ? '(' : '[')

// Right parenthesis: ')' for Matlab, ']' for C
#define RPAR(output_type) (OFFSET(output_type) ? ')' : ']')

// Computing cost above which a node can be declared a temporary term
#define MIN_COST_MATLAB (40*90)
#define MIN_COST_C (40*4)
#define MIN_COST(is_matlab) (is_matlab ? MIN_COST_MATLAB : MIN_COST_C)

//! Base class for expression nodes
class ExprNode
{
  friend class DataTree;
  friend class ModelTree;
  friend class ExprNodeLess;
  friend class NumConstNode;
  friend class VariableNode;
  friend class UnaryOpNode;
  friend class BinaryOpNode;
  friend class TrinaryOpNode;

private:
  //! Computes derivative w.r. to variable varID (but doesn't store it in derivatives map)
  /*! You shoud use getDerivative() to get the benefit of symbolic a priori and of caching */
  virtual NodeID computeDerivative(int varID) = 0;

protected:
  //! Reference to the enclosing DataTree
  DataTree &datatree;

  //! Index number
  int idx;

  //! Set of variable IDs with respect to which the derivative is potentially non-null
  set<int> non_null_derivatives;

  //! Used for caching of first order derivatives (when non-null)
  map<int, NodeID> derivatives;

  //! Cost of computing current node
  /*! Nodes included in temporary_terms are considered having a null cost */
  virtual int cost(const temporary_terms_type &temporary_terms, bool is_matlab) const;

public:
  ExprNode(DataTree &datatree_arg);
  virtual ~ExprNode();

  //! Returns derivative w.r. to variable varID
  /*! Uses a symbolic a priori to pre-detect null derivatives, and caches the result for other derivatives (to avoid computing it several times)
   For an equal node, returns the derivative of lhs minus rhs */
  NodeID getDerivative(int varID);

  //! Returns precedence of node
  /*! Equals 100 for constants, variables, unary ops, and temporary terms */
  virtual int precedence(ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms) const;

  //! Fills temporary_terms set, using reference counts
  /*! A node will be marked as a temporary term if it is referenced at least two times (i.e. has at least two parents), and has a computing cost (multiplied by reference count) greater to datatree.min_cost */
  virtual void computeTemporaryTerms(map<NodeID, int> &reference_count, temporary_terms_type &temporary_terms, bool is_matlab) const;

  //! Writes output of node, using a Txxx notation for nodes in temporary_terms
  virtual void writeOutput(ostream &output, ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms) const = 0;

  //! Writes output of node (with no temporary terms and with "outside model" output type)
  void writeOutput(ostream &output);

  //! Computes the set of endogenous variables in the expression
  /*! Endogenous are stored as integer pairs of the form (symb_id, lag)
      They are added to the set given in argument */
  virtual void collectEndogenous(set<pair<int, int> > &result) const = 0;
  virtual void computeTemporaryTerms(map<NodeID, int> &reference_count,
                                     temporary_terms_type &temporary_terms,
                                     map<NodeID, int> &first_occurence,
                                     int Curr_block,
                                     Model_Block *ModelBlock,
                                     map_idx_type &map_idx) const;

  class EvalException
  {
  };

  virtual double eval(const eval_context_type &eval_context) const throw (EvalException) = 0;
  virtual void compile(ofstream &CompileCode, bool lhs_rhs, ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms, map_idx_type map_idx) const = 0;
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
  virtual NodeID computeDerivative(int varID);
public:
  NumConstNode(DataTree &datatree_arg, int id_arg);
  virtual void writeOutput(ostream &output, ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms) const;
  virtual void collectEndogenous(set<pair<int, int> > &result) const;
  virtual double eval(const eval_context_type &eval_context) const throw (EvalException);
  virtual void compile(ofstream &CompileCode, bool lhs_rhs, ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms, map_idx_type map_idx) const;
};

//! Symbol or variable node
class VariableNode : public ExprNode
{
private:
  //! Id from the symbol table
  const int symb_id;
  const Type type;
  const int lag;
  //! Id from the variable table (-1 if not a endogenous/exogenous/recursive)
  int var_id;
  virtual NodeID computeDerivative(int varID);
public:
  VariableNode(DataTree &datatree_arg, int symb_id_arg, Type type_arg, int lag_arg);
  virtual void writeOutput(ostream &output, ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms = temporary_terms_type()) const;
  virtual void collectEndogenous(set<pair<int, int> > &result) const;
  virtual double eval(const eval_context_type &eval_context) const throw (EvalException);
  virtual void compile(ofstream &CompileCode, bool lhs_rhs, ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms, map_idx_type map_idx) const;
};

//! Unary operator node
class UnaryOpNode : public ExprNode
{
  friend class DataTree;
private:
  const NodeID arg;
  const UnaryOpcode op_code;
  virtual NodeID computeDerivative(int varID);

  virtual int cost(const temporary_terms_type &temporary_terms, bool is_matlab) const;
public:
  UnaryOpNode(DataTree &datatree_arg, UnaryOpcode op_code_arg, const NodeID arg_arg);
  virtual void computeTemporaryTerms(map<NodeID, int> &reference_count, temporary_terms_type &temporary_terms, bool is_matlab) const;
  virtual void writeOutput(ostream &output, ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms) const;
  virtual void computeTemporaryTerms(map<NodeID, int> &reference_count,
                                     temporary_terms_type &temporary_terms,
                                     map<NodeID, int> &first_occurence,
                                     int Curr_block,
                                     Model_Block *ModelBlock,
                                     map_idx_type &map_idx) const;
  virtual void collectEndogenous(set<pair<int, int> > &result) const;
  static double eval_opcode(UnaryOpcode op_code, double v) throw (EvalException);
  virtual double eval(const eval_context_type &eval_context) const throw (EvalException);
  virtual void compile(ofstream &CompileCode, bool lhs_rhs, ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms, map_idx_type map_idx) const;
};

//! Binary operator node
class BinaryOpNode : public ExprNode
{
  friend class ModelTree;
private:
  const NodeID arg1, arg2;
  const BinaryOpcode op_code;
  virtual NodeID computeDerivative(int varID);
  virtual int cost(const temporary_terms_type &temporary_terms, bool is_matlab) const;
public:
  BinaryOpNode(DataTree &datatree_arg, const NodeID arg1_arg,
               BinaryOpcode op_code_arg, const NodeID arg2_arg);
  virtual int precedence(ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms) const;
  virtual void computeTemporaryTerms(map<NodeID, int> &reference_count, temporary_terms_type &temporary_terms, bool is_matlab) const;
  virtual void writeOutput(ostream &output, ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms) const;
  virtual void computeTemporaryTerms(map<NodeID, int> &reference_count,
                                     temporary_terms_type &temporary_terms,
                                     map<NodeID, int> &first_occurence,
                                     int Curr_block,
                                     Model_Block *ModelBlock,
                                     map_idx_type &map_idx) const;
  virtual void collectEndogenous(set<pair<int, int> > &result) const;
  static double eval_opcode(double v1, BinaryOpcode op_code, double v2) throw (EvalException);
  virtual double eval(const eval_context_type &eval_context) const throw (EvalException);
  virtual void compile(ofstream &CompileCode, bool lhs_rhs, ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms, map_idx_type map_idx) const;
};

enum TrinaryOpcode
  {
    oNormcdf
  };

//! Trinary operator node
class TrinaryOpNode : public ExprNode
{
  friend class ModelTree;
private:
  const NodeID arg1, arg2, arg3;
  const TrinaryOpcode op_code;
  virtual NodeID computeDerivative(int varID);
  virtual int cost(const temporary_terms_type &temporary_terms, bool is_matlab) const;
public:
  TrinaryOpNode(DataTree &datatree_arg, const NodeID arg1_arg,
		TrinaryOpcode op_code_arg, const NodeID arg2_arg, const NodeID arg3_arg);
  virtual int precedence(ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms) const;
  virtual void computeTemporaryTerms(map<NodeID, int> &reference_count, temporary_terms_type &temporary_terms, bool is_matlab) const;
  virtual void writeOutput(ostream &output, ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms) const;
  virtual void computeTemporaryTerms(map<NodeID, int> &reference_count,
                                     temporary_terms_type &temporary_terms,
                                     map<NodeID, int> &first_occurence,
                                     int Curr_block,
                                     Model_Block *ModelBlock,
                                     map_idx_type &map_idx) const;
  virtual void collectEndogenous(set<pair<int, int> > &result) const;
  static double eval_opcode(double v1, TrinaryOpcode op_code, double v2, double v3) throw (EvalException);
  virtual double eval(const eval_context_type &eval_context) const throw (EvalException);
  virtual void compile(ofstream &CompileCode, bool lhs_rhs, ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms, map_idx_type map_idx) const;
};

//! Unknown function node
class UnknownFunctionNode : public ExprNode
{
private:
  //! Symbol ID (no need to store type: it is necessary eUnknownFunction)
  const int symb_id;
  const vector<NodeID> arguments;
  virtual NodeID computeDerivative(int varID);
public:
  UnknownFunctionNode(DataTree &datatree_arg, int symb_id_arg,
                      const vector<NodeID> &arguments_arg);
  virtual void computeTemporaryTerms(map<NodeID, int> &reference_count, temporary_terms_type &temporary_terms, bool is_matlab) const;
  virtual void writeOutput(ostream &output, ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms) const;
  virtual void computeTemporaryTerms(map<NodeID, int> &reference_count,
                                     temporary_terms_type &temporary_terms,
                                     map<NodeID, int> &first_occurence,
                                     int Curr_block,
                                     Model_Block *ModelBlock,
                                     map_idx_type &map_idx) const;
  virtual void collectEndogenous(set<pair<int, int> > &result) const;
  virtual double eval(const eval_context_type &eval_context) const throw (EvalException);
  virtual void compile(ofstream &CompileCode, bool lhs_rhs, ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms, map_idx_type map_idx) const;
};

typedef struct IM_compact
{
  int size, u_init, u_finish, nb_endo;
  int *u, *us, *Var, *Equ, *Var_Index, *Equ_Index, *Var_dyn_Index;
};

typedef struct Block
{
  int Size, Sized, Type, Simulation_Type, Max_Lead, Max_Lag, Nb_Lead_Lag_Endo;
  bool is_linear;
  int *Equation, *Own_Derivative;
  int *Variable, *Variable_Sorted, *dVariable;
  int *variable_dyn_index, *variable_dyn_leadlag;
  temporary_terms_type *Temporary_terms;
  IM_compact *IM_lead_lag;
  int Code_Start, Code_Length;
};

typedef struct Model_Block
{
  int Size, Periods;
  Block* Block_List;
  int *in_Block_Equ, *in_Block_Var, *in_Equ_of_Block, *in_Var_of_Block;
};

#endif
