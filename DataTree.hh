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
#include "VariableTable.hh"
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
  friend class UnknownFunctionNode;
protected:
  //! A reference to the symbol table
  SymbolTable &symbol_table;
  //! Reference to numerical constants table
  NumericalConstants &num_constants;

  typedef map<int, NodeID> num_const_node_map_type;
  num_const_node_map_type num_const_node_map;
  //! Pair (symbol_id, lag) used as key
  typedef map<pair<int, int>, NodeID> variable_node_map_type;
  variable_node_map_type variable_node_map;
  typedef map<pair<NodeID, int>, NodeID> unary_op_node_map_type;
  unary_op_node_map_type unary_op_node_map;
  typedef map<pair<pair<NodeID, NodeID>, int>, NodeID> binary_op_node_map_type;
  binary_op_node_map_type binary_op_node_map;
  typedef map<pair<pair<pair<NodeID, NodeID>,NodeID>, int>, NodeID> trinary_op_node_map_type;
  trinary_op_node_map_type trinary_op_node_map;

  //! Stores local variables value (maps symbol ID to corresponding node)
  map<int, NodeID> local_variables_table;

  //! Internal implementation of AddVariable(), without the check on the lag
  NodeID AddVariableInternal(const string &name, int lag);

private:
  typedef list<NodeID> node_list_type;
  //! The list of nodes
  node_list_type node_list;
  //! A counter for filling ExprNode's idx field
  int node_counter;

  inline NodeID AddPossiblyNegativeConstant(double val);
  inline NodeID AddUnaryOp(UnaryOpcode op_code, NodeID arg);
  inline NodeID AddBinaryOp(NodeID arg1, BinaryOpcode op_code, NodeID arg2);
  inline NodeID AddTrinaryOp(NodeID arg1, TrinaryOpcode op_code, NodeID arg2, NodeID arg3);

public:
  DataTree(SymbolTable &symbol_table_arg, NumericalConstants &num_constants_arg);
  virtual ~DataTree();
  //! The variable table
  VariableTable variable_table;
  //! Some predefined constants
  NodeID Zero, One, Two, MinusOne, NaN, Infinity, MinusInfinity, Pi;

  //! Raised when a local parameter is declared twice
  class LocalVariableException
  {
  public:
    string name;
    LocalVariableException(const string &name_arg) : name(name_arg) {}
  };

  //! Adds a numerical constant
  NodeID AddNumConstant(const string &value);
  //! Adds a variable
  /*! The default implementation of the method refuses any lag != 0 */
  virtual NodeID AddVariable(const string &name, int lag = 0);
  //! Adds "arg1+arg2" to model tree
  NodeID AddPlus(NodeID iArg1, NodeID iArg2);
  //! Adds "arg1-arg2" to model tree
  NodeID AddMinus(NodeID iArg1, NodeID iArg2);
  //! Adds "-arg" to model tree
  NodeID AddUMinus(NodeID iArg1);
  //! Adds "arg1*arg2" to model tree
  NodeID AddTimes(NodeID iArg1, NodeID iArg2);
  //! Adds "arg1/arg2" to model tree
  NodeID AddDivide(NodeID iArg1, NodeID iArg2);
  //! Adds "arg1<arg2" to model tree
  NodeID AddLess(NodeID iArg1, NodeID iArg2);
  //! Adds "arg1>arg2" to model tree
  NodeID AddGreater(NodeID iArg1, NodeID iArg2);
  //! Adds "arg1<=arg2" to model tree
  NodeID AddLessEqual(NodeID iArg1, NodeID iArg2);
  //! Adds "arg1>=arg2" to model tree
  NodeID AddGreaterEqual(NodeID iArg1, NodeID iArg2);
  //! Adds "arg1==arg2" to model tree
  NodeID AddEqualEqual(NodeID iArg1, NodeID iArg2);
  //! Adds "arg1!=arg2" to model tree
  NodeID AddDifferent(NodeID iArg1, NodeID iArg2);
  //! Adds "arg1^arg2" to model tree
  NodeID AddPower(NodeID iArg1, NodeID iArg2);
  //! Adds "exp(arg)" to model tree
  NodeID AddExp(NodeID iArg1);
  //! Adds "log(arg)" to model tree
  NodeID AddLog(NodeID iArg1);
  //! Adds "log10(arg)" to model tree
  NodeID AddLog10(NodeID iArg1);
  //! Adds "cos(arg)" to model tree
  NodeID AddCos(NodeID iArg1);
  //! Adds "sin(arg)" to model tree
  NodeID AddSin(NodeID iArg1);
  //! Adds "tan(arg)" to model tree
  NodeID AddTan(NodeID iArg1);
  //! Adds "acos(arg)" to model tree
  NodeID AddAcos(NodeID iArg1);
  //! Adds "asin(arg)" to model tree
  NodeID AddAsin(NodeID iArg1);
  //! Adds "atan(arg)" to model tree
  NodeID AddAtan(NodeID iArg1);
  //! Adds "cosh(arg)" to model tree
  NodeID AddCosh(NodeID iArg1);
  //! Adds "sinh(arg)" to model tree
  NodeID AddSinh(NodeID iArg1);
  //! Adds "tanh(arg)" to model tree
  NodeID AddTanh(NodeID iArg1);
  //! Adds "acosh(arg)" to model tree
  NodeID AddAcosh(NodeID iArg1);
  //! Adds "asinh(arg)" to model tree
  NodeID AddAsinh(NodeID iArg1);
  //! Adds "atanh(args)" to model tree
  NodeID AddAtanh(NodeID iArg1);
  //! Adds "sqrt(arg)" to model tree
  NodeID AddSqrt(NodeID iArg1);
  //! Adds "max(arg1,arg2)" to model tree
  NodeID AddMax(NodeID iArg1, NodeID iArg2);
  //! Adds "min(arg1,arg2)" to model tree
  NodeID AddMin(NodeID iArg1, NodeID iArg2);
  //! Adds "normcdf(arg1,arg2,arg3)" to model tree
  NodeID AddNormcdf(NodeID iArg1, NodeID iArg2, NodeID iArg3);
  //! Adds "arg1=arg2" to model tree
  NodeID AddEqual(NodeID iArg1, NodeID iArg2);
  //! Adds a model local variable with its value
  void AddLocalVariable(const string &name, NodeID value) throw (LocalVariableException);
  //! Adds an unknown function node
  /*! \todo Use a map to share identical nodes */
  NodeID AddUnknownFunction(const string &function_name, const vector<NodeID> &arguments);
  //! Fill eval context with values of local variables
  void fillEvalContext(eval_context_type &eval_context) const;
  //! Checks if a given symbol is used somewhere in the data tree
  bool isSymbolUsed(int symb_id) const;
};

inline NodeID
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

  NodeID cnode = AddNumConstant(ost.str());

  if (neg)
    return AddUMinus(cnode);
  else
    return cnode;
}

inline NodeID
DataTree::AddUnaryOp(UnaryOpcode op_code, NodeID arg)
{
  // If the node already exists in tree, share it
  unary_op_node_map_type::iterator it = unary_op_node_map.find(make_pair(arg, op_code));
  if (it != unary_op_node_map.end())
    return it->second;

  // Try to reduce to a constant
  // Case where arg is a constant and op_code == oUminus (i.e. we're adding a negative constant) is skipped
  NumConstNode *carg = dynamic_cast<NumConstNode *>(arg);
  if (op_code != oUminus || carg == NULL)
    {
      try
        {
          double argval = arg->eval(eval_context_type());
          double val = UnaryOpNode::eval_opcode(op_code, argval);
          return AddPossiblyNegativeConstant(val);
        }
      catch(ExprNode::EvalException &e)
        {
        }
    }
  return new UnaryOpNode(*this, op_code, arg);
}

inline NodeID
DataTree::AddBinaryOp(NodeID arg1, BinaryOpcode op_code, NodeID arg2)
{
  binary_op_node_map_type::iterator it = binary_op_node_map.find(make_pair(make_pair(arg1, arg2), op_code));
  if (it != binary_op_node_map.end())
    return it->second;

  // Try to reduce to a constant
  try
    {
      double argval1 = arg1->eval(eval_context_type());
      double argval2 = arg2->eval(eval_context_type());
      double val = BinaryOpNode::eval_opcode(argval1, op_code, argval2);
      return AddPossiblyNegativeConstant(val);
    }
  catch(ExprNode::EvalException &e)
    {
    }
  return new BinaryOpNode(*this, arg1, op_code, arg2);
}

inline NodeID
DataTree::AddTrinaryOp(NodeID arg1, TrinaryOpcode op_code, NodeID arg2, NodeID arg3)
{
  trinary_op_node_map_type::iterator it = trinary_op_node_map.find(make_pair(make_pair(make_pair(arg1, arg2), arg3), op_code));
  if (it != trinary_op_node_map.end())
    return it->second;

  // Try to reduce to a constant
  try
    {
      double argval1 = arg1->eval(eval_context_type());
      double argval2 = arg2->eval(eval_context_type());
      double argval3 = arg3->eval(eval_context_type());
      double val = TrinaryOpNode::eval_opcode(argval1, op_code, argval2, argval3);
      return AddPossiblyNegativeConstant(val);
    }
  catch(ExprNode::EvalException &e)
    {
    }
  return new TrinaryOpNode(*this, arg1, op_code, arg2, arg3);
}

#endif
