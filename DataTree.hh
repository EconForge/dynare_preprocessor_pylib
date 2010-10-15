/*
 * Copyright (C) 2003-2010 Dynare Team
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

#ifndef _DATATREE_HH
#define _DATATREE_HH

using namespace std;

#include <string>
#include <map>
#include <list>
#include <sstream>
#include <iomanip>

#include "SymbolTable.hh"
#include "NumericalConstants.hh"
#include "ExternalFunctionsTable.hh"
#include "ExprNode.hh"

#define CONSTANTS_PRECISION 16

class DataTree
{
  friend class ExprNode;
  friend class NumConstNode;
  friend class VariableNode;
  friend class UnaryOpNode;
  friend class BinaryOpNode;
  friend class TrinaryOpNode;
  friend class ExternalFunctionNode;
  friend class FirstDerivExternalFunctionNode;
  friend class SecondDerivExternalFunctionNode;
protected:
  //! A reference to the symbol table
  SymbolTable &symbol_table;
  //! Reference to numerical constants table
  NumericalConstants &num_constants;
  //! A reference to the external functions table
  ExternalFunctionsTable &external_functions_table;

  typedef map<int, NumConstNode *> num_const_node_map_t;
  num_const_node_map_t num_const_node_map;
  //! Pair (symbol_id, lag) used as key
  typedef map<pair<int, int>, VariableNode *> variable_node_map_t;
  variable_node_map_t variable_node_map;
  typedef map<pair<expr_t, UnaryOpcode>, UnaryOpNode *> unary_op_node_map_t;
  unary_op_node_map_t unary_op_node_map;
  typedef map<pair<pair<expr_t, expr_t>, BinaryOpcode>, BinaryOpNode *> binary_op_node_map_t;
  binary_op_node_map_t binary_op_node_map;
  typedef map<pair<pair<pair<expr_t, expr_t>, expr_t>, TrinaryOpcode>, TrinaryOpNode *> trinary_op_node_map_t;
  trinary_op_node_map_t trinary_op_node_map;
  typedef map<pair<vector<expr_t>, int>, ExternalFunctionNode *> external_function_node_map_t;
  external_function_node_map_t external_function_node_map;
  typedef map<pair<pair<vector<expr_t>, int>, int>, FirstDerivExternalFunctionNode *> first_deriv_external_function_node_map_t;
  first_deriv_external_function_node_map_t first_deriv_external_function_node_map;
  typedef map<pair<pair<vector<expr_t>, pair<int, int> >, int>, SecondDerivExternalFunctionNode *> second_deriv_external_function_node_map_t;
  second_deriv_external_function_node_map_t second_deriv_external_function_node_map;

  //! Stores local variables value (maps symbol ID to corresponding node)
  map<int, expr_t> local_variables_table;

  //! Internal implementation of AddVariable(), without the check on the lag
  VariableNode *AddVariableInternal(int symb_id, int lag);

private:
  typedef list<expr_t> node_list_t;
  //! The list of nodes
  node_list_t node_list;
  //! A counter for filling ExprNode's idx field
  int node_counter;

  inline expr_t AddPossiblyNegativeConstant(double val);
  inline expr_t AddUnaryOp(UnaryOpcode op_code, expr_t arg, int arg_exp_info_set = 0, const string &arg_exp_info_set_name="");
  inline expr_t AddBinaryOp(expr_t arg1, BinaryOpcode op_code, expr_t arg2);
  inline expr_t AddTrinaryOp(expr_t arg1, TrinaryOpcode op_code, expr_t arg2, expr_t arg3);

public:
  DataTree(SymbolTable &symbol_table_arg, NumericalConstants &num_constants_arg, ExternalFunctionsTable &external_functions_table_arg);
  virtual ~DataTree();

  //! Some predefined constants
  expr_t Zero, One, Two, MinusOne, NaN, Infinity, MinusInfinity, Pi;

  //! Raised when a local parameter is declared twice
  class LocalVariableException
  {
  public:
    string name;
    LocalVariableException(const string &name_arg) : name(name_arg)
    {
    }
  };

  //! Adds a numerical constant
  expr_t AddNumConstant(const string &value);
  //! Adds a variable
  /*! The default implementation of the method refuses any lag != 0 */
  virtual VariableNode *AddVariable(int symb_id, int lag = 0);
  //! Adds "arg1+arg2" to model tree
  expr_t AddPlus(expr_t iArg1, expr_t iArg2);
  //! Adds "arg1-arg2" to model tree
  expr_t AddMinus(expr_t iArg1, expr_t iArg2);
  //! Adds "-arg" to model tree
  expr_t AddUMinus(expr_t iArg1);
  //! Adds "arg1*arg2" to model tree
  expr_t AddTimes(expr_t iArg1, expr_t iArg2);
  //! Adds "arg1/arg2" to model tree
  expr_t AddDivide(expr_t iArg1, expr_t iArg2);
  //! Adds "arg1<arg2" to model tree
  expr_t AddLess(expr_t iArg1, expr_t iArg2);
  //! Adds "arg1>arg2" to model tree
  expr_t AddGreater(expr_t iArg1, expr_t iArg2);
  //! Adds "arg1<=arg2" to model tree
  expr_t AddLessEqual(expr_t iArg1, expr_t iArg2);
  //! Adds "arg1>=arg2" to model tree
  expr_t AddGreaterEqual(expr_t iArg1, expr_t iArg2);
  //! Adds "arg1==arg2" to model tree
  expr_t AddEqualEqual(expr_t iArg1, expr_t iArg2);
  //! Adds "arg1!=arg2" to model tree
  expr_t AddDifferent(expr_t iArg1, expr_t iArg2);
  //! Adds "arg1^arg2" to model tree
  expr_t AddPower(expr_t iArg1, expr_t iArg2);
  //! Adds "E(arg1)(arg2)" to model tree
  expr_t AddExpectation(int iArg1, expr_t iArg2);
  //! Adds "E(arg1)(arg2)" to model tree
  expr_t AddExpectation(string *iArg1, expr_t iArg2);
  //! Adds "exp(arg)" to model tree
  expr_t AddExp(expr_t iArg1);
  //! Adds "log(arg)" to model tree
  expr_t AddLog(expr_t iArg1);
  //! Adds "log10(arg)" to model tree
  expr_t AddLog10(expr_t iArg1);
  //! Adds "cos(arg)" to model tree
  expr_t AddCos(expr_t iArg1);
  //! Adds "sin(arg)" to model tree
  expr_t AddSin(expr_t iArg1);
  //! Adds "tan(arg)" to model tree
  expr_t AddTan(expr_t iArg1);
  //! Adds "acos(arg)" to model tree
  expr_t AddAcos(expr_t iArg1);
  //! Adds "asin(arg)" to model tree
  expr_t AddAsin(expr_t iArg1);
  //! Adds "atan(arg)" to model tree
  expr_t AddAtan(expr_t iArg1);
  //! Adds "cosh(arg)" to model tree
  expr_t AddCosh(expr_t iArg1);
  //! Adds "sinh(arg)" to model tree
  expr_t AddSinh(expr_t iArg1);
  //! Adds "tanh(arg)" to model tree
  expr_t AddTanh(expr_t iArg1);
  //! Adds "acosh(arg)" to model tree
  expr_t AddAcosh(expr_t iArg1);
  //! Adds "asinh(arg)" to model tree
  expr_t AddAsinh(expr_t iArg1);
  //! Adds "atanh(args)" to model tree
  expr_t AddAtanh(expr_t iArg1);
  //! Adds "sqrt(arg)" to model tree
  expr_t AddSqrt(expr_t iArg1);
  //! Adds "erf(arg)" to model tree
  expr_t AddErf(expr_t iArg1);
  //! Adds "max(arg1,arg2)" to model tree
  expr_t AddMax(expr_t iArg1, expr_t iArg2);
  //! Adds "min(arg1,arg2)" to model tree
  expr_t AddMin(expr_t iArg1, expr_t iArg2);
  //! Adds "normcdf(arg1,arg2,arg3)" to model tree
  expr_t AddNormcdf(expr_t iArg1, expr_t iArg2, expr_t iArg3);
  //! Adds "normpdf(arg1,arg2,arg3)" to model tree
  expr_t AddNormpdf(expr_t iArg1, expr_t iArg2, expr_t iArg3);
  //! Adds "steadyState(arg)" to model tree
  expr_t AddSteadyState(expr_t iArg1);
  //! Adds "arg1=arg2" to model tree
  expr_t AddEqual(expr_t iArg1, expr_t iArg2);
  //! Adds a model local variable with its value
  void AddLocalVariable(int symb_id, expr_t value) throw (LocalVariableException);
  //! Adds an external function node
  expr_t AddExternalFunction(int symb_id, const vector<expr_t> &arguments);
  //! Adds an external function node for the first derivative of an external function
  expr_t AddFirstDerivExternalFunctionNode(int top_level_symb_id, const vector<expr_t> &arguments, int input_index);
  //! Adds an external function node for the second derivative of an external function
  expr_t AddSecondDerivExternalFunctionNode(int top_level_symb_id, const vector<expr_t> &arguments, int input_index1, int input_index2);
  //! Checks if a given symbol is used somewhere in the data tree
  bool isSymbolUsed(int symb_id) const;
  //! Checks if a given unary op is used somewhere in the data tree
  bool isUnaryOpUsed(UnaryOpcode opcode) const;
  //! Checks if a given binary op is used somewhere in the data tree
  bool isBinaryOpUsed(BinaryOpcode opcode) const;
  //! Checks if a given trinary op is used somewhere in the data tree
  bool isTrinaryOpUsed(TrinaryOpcode opcode) const;
  //! Checks if a given external function is used somewhere in the data tree
  bool isExternalFunctionUsed(int symb_id) const;
  //! Checks if a given first derivative external function is used somewhere in the data tree
  bool isFirstDerivExternalFunctionUsed(int symb_id) const;
  //! Checks if a given second derivative external function is used somewhere in the data tree
  bool isSecondDerivExternalFunctionUsed(int symb_id) const;
  //! Thrown when trying to access an unknown variable by deriv_id
  class UnknownDerivIDException
  {
  };

  //! Raised when a trend is declared twice
  class TrendException
  {
  public:
    string name;
    TrendException(const string &name_arg) : name(name_arg)
    {
    }
  };

  //! Returns the derivation ID, or throws an exception if the derivation ID does not exist
  virtual int getDerivID(int symb_id, int lag) const throw (UnknownDerivIDException);
  //! Returns the column of the dynamic Jacobian associated to a derivation ID
  virtual int getDynJacobianCol(int deriv_id) const throw (UnknownDerivIDException);

  //! Returns bool indicating whether DataTree represents a Dynamic Model (returns true in DynamicModel.hh)
  virtual bool
  isDynamic() const
  {
    return false;
  };
};

inline expr_t
DataTree::AddPossiblyNegativeConstant(double v)
{
  bool neg = false;
  if (v < 0)
    {
      v = -v;
      neg = true;
    }
  ostringstream ost;
  ost << setprecision(CONSTANTS_PRECISION) << v;

  expr_t cnode = AddNumConstant(ost.str());

  if (neg)
    return AddUMinus(cnode);
  else
    return cnode;
}

inline expr_t
DataTree::AddUnaryOp(UnaryOpcode op_code, expr_t arg, int arg_exp_info_set, const string &arg_exp_info_set_name)
{
  // If the node already exists in tree, share it
  unary_op_node_map_t::iterator it = unary_op_node_map.find(make_pair(arg, op_code));
  if (it != unary_op_node_map.end())
    return it->second;

  // Try to reduce to a constant
  // Case where arg is a constant and op_code == oUminus (i.e. we're adding a negative constant) is skipped
  NumConstNode *carg = dynamic_cast<NumConstNode *>(arg);
  if (op_code != oUminus || carg == NULL)
    {
      try
        {
          double argval = arg->eval(eval_context_t());
          double val = UnaryOpNode::eval_opcode(op_code, argval);
          return AddPossiblyNegativeConstant(val);
        }
      catch (ExprNode::EvalException &e)
        {
        }
    }
  return new UnaryOpNode(*this, op_code, arg, arg_exp_info_set, arg_exp_info_set_name);
}

inline expr_t
DataTree::AddBinaryOp(expr_t arg1, BinaryOpcode op_code, expr_t arg2)
{
  binary_op_node_map_t::iterator it = binary_op_node_map.find(make_pair(make_pair(arg1, arg2), op_code));
  if (it != binary_op_node_map.end())
    return it->second;

  // Try to reduce to a constant
  try
    {
      double argval1 = arg1->eval(eval_context_t());
      double argval2 = arg2->eval(eval_context_t());
      double val = BinaryOpNode::eval_opcode(argval1, op_code, argval2);
      return AddPossiblyNegativeConstant(val);
    }
  catch (ExprNode::EvalException &e)
    {
    }
  return new BinaryOpNode(*this, arg1, op_code, arg2);
}

inline expr_t
DataTree::AddTrinaryOp(expr_t arg1, TrinaryOpcode op_code, expr_t arg2, expr_t arg3)
{
  trinary_op_node_map_t::iterator it = trinary_op_node_map.find(make_pair(make_pair(make_pair(arg1, arg2), arg3), op_code));
  if (it != trinary_op_node_map.end())
    return it->second;

  // Try to reduce to a constant
  try
    {
      double argval1 = arg1->eval(eval_context_t());
      double argval2 = arg2->eval(eval_context_t());
      double argval3 = arg3->eval(eval_context_t());
      double val = TrinaryOpNode::eval_opcode(argval1, op_code, argval2, argval3);
      return AddPossiblyNegativeConstant(val);
    }
  catch (ExprNode::EvalException &e)
    {
    }
  return new TrinaryOpNode(*this, arg1, op_code, arg2, arg3);
}

#endif
