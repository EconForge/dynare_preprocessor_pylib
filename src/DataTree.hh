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

#ifndef _DATATREE_HH
#define _DATATREE_HH

#include <string>
#include <map>
#include <vector>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <utility>
#include <memory>
#include <filesystem>
#include <string_view>

#include "SymbolTable.hh"
#include "NumericalConstants.hh"
#include "ExternalFunctionsTable.hh"
#include "ExprNode.hh"
#include "SubModel.hh"

using namespace std;

class DataTree
{
public:
  //! A reference to the symbol table
  SymbolTable &symbol_table;
  //! Reference to numerical constants table
  NumericalConstants &num_constants;
  //! A reference to the external functions table
  ExternalFunctionsTable &external_functions_table;
  //! Is it possible to use leads/lags on variable nodes?
  const bool is_dynamic;

private:
  //! num_constant_id -> NumConstNode
  using num_const_node_map_t = map<int, NumConstNode *>;
  num_const_node_map_t num_const_node_map;

  //! (symbol_id, lag) -> VariableNode
  using variable_node_map_t = map<pair<int, int>, VariableNode *>;
  variable_node_map_t variable_node_map;

  //! (arg, op_code, arg_exp_info_set, param1_symb_id, param2_symb_id, adl_param_name, adl_lags) -> UnaryOpNode
  using unary_op_node_map_t = map<tuple<expr_t, UnaryOpcode, int, int, int, string, vector<int>>, UnaryOpNode *>;
  unary_op_node_map_t unary_op_node_map;

  //! ( arg1, arg2, opCode, order of Power Derivative) -> BinaryOpNode
  using binary_op_node_map_t = map<tuple<expr_t, expr_t, BinaryOpcode, int>, BinaryOpNode *>;
  binary_op_node_map_t binary_op_node_map;

  //! ( arg1, arg2, arg3, opCode) -> TrinaryOpNode
  using trinary_op_node_map_t = map<tuple<expr_t, expr_t, expr_t, TrinaryOpcode>, TrinaryOpNode *>;
  trinary_op_node_map_t trinary_op_node_map;

  // (arguments, symb_id) -> ExternalFunctionNode
  using external_function_node_map_t = map<pair<vector<expr_t>, int>, ExternalFunctionNode *>;
  external_function_node_map_t external_function_node_map;

  // (model_name, symb_id, forecast_horizon) -> VarExpectationNode
  using var_expectation_node_map_t = map<string, VarExpectationNode *>;
  var_expectation_node_map_t var_expectation_node_map;

  // model_name -> PacExpectationNode
  using pac_expectation_node_map_t = map<string, PacExpectationNode *>;
  pac_expectation_node_map_t pac_expectation_node_map;

  // model_name -> PacTargetNonstationaryNode
  using pac_target_nonstationary_node_map_t = map<string, PacTargetNonstationaryNode *>;
  pac_target_nonstationary_node_map_t pac_target_nonstationary_node_map;

  // (arguments, deriv_idx, symb_id) -> FirstDerivExternalFunctionNode
  using first_deriv_external_function_node_map_t = map<tuple<vector<expr_t>, int, int>, FirstDerivExternalFunctionNode *>;
  first_deriv_external_function_node_map_t first_deriv_external_function_node_map;

  // (arguments, deriv_idx1, deriv_idx2, symb_id) -> SecondDerivExternalFunctionNode
  using second_deriv_external_function_node_map_t = map<tuple<vector<expr_t>, int, int, int>, SecondDerivExternalFunctionNode *>;
  second_deriv_external_function_node_map_t second_deriv_external_function_node_map;

  // Flag to disable simplifications related to commutativity of addition and multiplication
  static bool no_commutativity;

protected:
  //! Stores local variables value (maps symbol ID to corresponding node)
  map<int, expr_t> local_variables_table;
  //! Stores the order of appearance of local variables in the model block. Needed following change in #563
  vector<int> local_variables_vector;

  //! Internal implementation of ParamUsedWithLeadLag()
  bool ParamUsedWithLeadLagInternal() const;

  /* Writes the contents of “new_contents” to the file “filename”. However, if
     the file already exists and would not be modified by this operation, then do
     nothing. */
  static void writeToFileIfModified(stringstream &new_contents, const filesystem::path &filename);

private:
  constexpr static int constants_precision{16};

  //! The list of nodes
  vector<unique_ptr<ExprNode>> node_list;

  inline expr_t AddUnaryOp(UnaryOpcode op_code, expr_t arg, int arg_exp_info_set = 0, int param1_symb_id = 0, int param2_symb_id = 0, const string &adl_param_name = "", const vector<int> &adl_lags = vector<int>());
  inline expr_t AddBinaryOp(expr_t arg1, BinaryOpcode op_code, expr_t arg2, int powerDerivOrder = 0);
  inline expr_t AddTrinaryOp(expr_t arg1, TrinaryOpcode op_code, expr_t arg2, expr_t arg3);

  //! Initializes the predefined constants, used only from the constructors
  void initConstants();

public:
  DataTree(SymbolTable &symbol_table_arg,
           NumericalConstants &num_constants_arg,
           ExternalFunctionsTable &external_functions_table_arg,
           bool is_static_args = false);

  virtual ~DataTree() = default;

  DataTree(const DataTree &d);
  DataTree &operator=(const DataTree &d);

  //! Some predefined constants
  NumConstNode *Zero, *One, *Two, *Three, *NaN, *Infinity, *Pi;
  expr_t MinusOne, MinusInfinity;

  //! Raised when a local parameter is declared twice
  struct LocalVariableException
  {
    string name;
  };

  class DivisionByZeroException
  {
  };

  inline expr_t AddPossiblyNegativeConstant(double val);
  //! Adds a non-negative numerical constant (possibly Inf or NaN)
  NumConstNode *AddNonNegativeConstant(const string &value);
  //! Adds a variable
  VariableNode *AddVariable(int symb_id, int lag = 0);
  //! Gets a variable
  /*! Same as AddVariable, except that it fails if the variable node has not
    already been created */
  VariableNode *getVariable(int symb_id, int lag = 0) const;
  //! Adds "arg1+arg2" to model tree
  expr_t AddPlus(expr_t iArg1, expr_t iArg2);
  //! Adds "arg1-arg2" to model tree
  expr_t AddMinus(expr_t iArg1, expr_t iArg2);
  //! Adds "-arg" to model tree
  expr_t AddUMinus(expr_t iArg1);
  //! Adds "arg1*arg2" to model tree
  expr_t AddTimes(expr_t iArg1, expr_t iArg2);
  //! Adds "arg1/arg2" to model tree
  expr_t AddDivide(expr_t iArg1, expr_t iArg2) noexcept(false);
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
  //! Adds "getPowerDeriv(arg1, arg2, powerDerivOrder)" to model tree
  expr_t AddPowerDeriv(expr_t iArg1, expr_t iArg2, int powerDerivOrder);
  //! Adds "E(arg1)(arg2)" to model tree
  expr_t AddExpectation(int iArg1, expr_t iArg2);
  //! Adds "diff(arg)" to model tree
  expr_t AddDiff(expr_t iArg1);
  //! Adds "adl(arg1, name, lag/lags)" to model tree
  expr_t AddAdl(expr_t iArg1, const string &name, const vector<int> &lags);
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
  //! Adds "cbrt(arg)" to model tree
  expr_t AddCbrt(expr_t iArg1);
  //! Adds "abs(arg)" to model tree
  expr_t AddAbs(expr_t iArg1);
  //! Adds "sign(arg)" to model tree
  expr_t AddSign(expr_t iArg1);
  //! Adds "erf(arg)" to model tree
  expr_t AddErf(expr_t iArg1);
  //! Adds "erfc(arg)" to model tree
  expr_t AddErfc(expr_t iArg1);
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
  //! Add derivative of steady state w.r.t. parameter to model tree
  expr_t AddSteadyStateParamDeriv(expr_t iArg1, int param_symb_id);
  //! Add 2nd derivative of steady state w.r.t. parameter to model tree
  expr_t AddSteadyStateParam2ndDeriv(expr_t iArg1, int param1_symb_id, int param2_symb_id);
  //! Adds "arg1=arg2" to model tree
  BinaryOpNode *AddEqual(expr_t iArg1, expr_t iArg2);
  //! Adds "var_expectation(model_name)" to model tree
  expr_t AddVarExpectation(const string &model_name);
  //! Adds pac_expectation command to model tree
  expr_t AddPacExpectation(const string &model_name);
  //! Adds a pac_target_nonstationary node to model tree
  expr_t AddPacTargetNonstationary(const string &model_name);
  //! Adds a model local variable with its value
  void AddLocalVariable(int symb_id, expr_t value) noexcept(false);
  //! Adds an external function node
  expr_t AddExternalFunction(int symb_id, const vector<expr_t> &arguments);
  //! Adds an external function node for the first derivative of an external function
  expr_t AddFirstDerivExternalFunction(int top_level_symb_id, const vector<expr_t> &arguments, int input_index);
  //! Adds an external function node for the second derivative of an external function
  expr_t AddSecondDerivExternalFunction(int top_level_symb_id, const vector<expr_t> &arguments, int input_index1, int input_index2);
  //! Checks if a given symbol is used somewhere in the data tree
  bool isSymbolUsed(int symb_id) const;
  //! Checks if a given unary op is used somewhere in the data tree
  bool isUnaryOpUsed(UnaryOpcode opcode) const;
  //! Checks if a given unary op is used somewhere in the data tree on an endogenous variable
  bool isUnaryOpUsedOnType(SymbolType type, UnaryOpcode opcode) const;
  //! Checks if a given binary op is used somewhere in the data tree
  bool isBinaryOpUsed(BinaryOpcode opcode) const;
  //! Checks if a given binary op is used somewhere in the data tree on an endogenous variable
  bool isBinaryOpUsedOnType(SymbolType type, BinaryOpcode opcode) const;
  //! Returns the minimum lag (as a negative number) of the given symbol in the whole data tree (and not only in the equations !!)
  /*! Returns 0 if the symbol is not used */
  int minLagForSymbol(int symb_id) const;
  /* Writes definitions of C function helpers (getPowerDeriv(), sign()) as
     inline functions */
  void writeCHelpersDefinition(ostream &output) const;
  /* Writes declarations of C function helpers (getPowerDeriv(), sign()) as
     extern inline (external definition). Those need to be included in exactly
     one translation unit. That external definition will be used or not,
     depending on the optimization decision by the compiler.
     See https://en.cppreference.com/w/c/language/inline */
  void writeCHelpersDeclaration(ostream &output) const;
  //! Thrown when trying to access an unknown variable by deriv_id
  class UnknownDerivIDException
  {
  };

  //! Raised when a trend is declared twice
  struct TrendException
  {
    string name;
  };

  // Returns the derivation ID, or throws an exception if the derivation ID does not exist
  virtual int getDerivID(int symb_id, int lag) const noexcept(false);
  // Get the type corresponding to a derivation ID
  virtual SymbolType getTypeByDerivID(int deriv_id) const noexcept(false);
  // Get the lag corresponding to a derivation ID
  virtual int getLagByDerivID(int deriv_id) const noexcept(false);
  // Get the symbol ID corresponding to a derivation ID
  virtual int getSymbIDByDerivID(int deriv_id) const noexcept(false);
  // Get the type-specific ID corresponding to a derivation ID
  virtual int getTypeSpecificIDByDerivID(int deriv_id) const;
  // Get the symbol name corresponding to a derivation ID
  string
  getNameByDerivID(int deriv_id) const
  {
    return symbol_table.getName(getSymbIDByDerivID(deriv_id));
  }

  /* Returns the column of the Jacobian associated to a derivation ID.
     The “sparse” argument selects between the legacy representation and the
     sparse representation. */
  virtual int
  getJacobianCol([[maybe_unused]] int deriv_id, [[maybe_unused]] bool sparse) const
  {
    throw UnknownDerivIDException();
  }

  /* Returns the number of columns of the Jacobian
     The “sparse” argument selects between the legacy representation and the
     sparse representation. */
  virtual int
  getJacobianColsNbr([[maybe_unused]] bool sparse) const
  {
    throw UnknownDerivIDException();
  }

  //! Adds to the set all the deriv IDs corresponding to parameters
  virtual void addAllParamDerivId(set<int> &deriv_id_set);

  //! Returns bool indicating whether DataTree represents a Dynamic Model (returns true in DynamicModel.hh)
  virtual bool
  isDynamic() const
  {
    return false;
  };

  struct UnknownLocalVariableException
  {
    //! Symbol ID
    int id;
  };

  expr_t
  getLocalVariable(int symb_id) const
  {
    auto it = local_variables_table.find(symb_id);
    if (it == local_variables_table.end())
      throw UnknownLocalVariableException{symb_id};

    return it->second;
  }

  static void
  setNoCommutativity()
  {
    no_commutativity = true;
  }

  /* Equivalent of MATLAB/Octave’s strsplit, except that it ignores empty
     substring components (MATLAB/Octave adds them to the output); in
     particular, returns an empty vector given an empty string. */
  static vector<string> strsplit(string_view str, char delim);

  /*! Takes a MATLAB/Octave package name (possibly with several levels nested using dots),
    and returns the path to the corresponding filesystem directory.
    In practice the package nesting is used for the planner_objective (stored
    inside +objective subdir). */
  static filesystem::path packageDir(string_view package);
};

inline expr_t
DataTree::AddPossiblyNegativeConstant(double v)
{
  /* Treat NaN and Inf separately. In particular, under Windows, converting
     them to a string does not work as expected */
  if (isnan(v))
    return NaN;
  if (isinf(v))
    return (v < 0 ? MinusInfinity : Infinity);

  bool neg = false;
  if (v < 0)
    {
      v = -v;
      neg = true;
    }
  ostringstream ost;
  ost << setprecision(constants_precision) << v;

  expr_t cnode = AddNonNegativeConstant(ost.str());

  if (neg)
    return AddUMinus(cnode);
  else
    return cnode;
}

inline expr_t
DataTree::AddUnaryOp(UnaryOpcode op_code, expr_t arg, int arg_exp_info_set, int param1_symb_id, int param2_symb_id, const string &adl_param_name, const vector<int> &adl_lags)
{
  // If the node already exists in tree, share it
  if (auto it = unary_op_node_map.find({ arg, op_code, arg_exp_info_set, param1_symb_id, param2_symb_id, adl_param_name, adl_lags });
      it != unary_op_node_map.end())
    return it->second;

  // Try to reduce to a constant
  // Case where arg is a constant and op_code == UnaryOpcode::uminus (i.e. we're adding a negative constant) is skipped
  if (auto carg = dynamic_cast<NumConstNode *>(arg);
      op_code != UnaryOpcode::uminus || !carg)
    {
      try
        {
          double argval = arg->eval({});
          double val = UnaryOpNode::eval_opcode(op_code, argval);
          return AddPossiblyNegativeConstant(val);
        }
      catch (ExprNode::EvalException &e)
        {
        }
    }

  auto sp = make_unique<UnaryOpNode>(*this, node_list.size(), op_code, arg, arg_exp_info_set, param1_symb_id, param2_symb_id, adl_param_name, adl_lags);
  auto p = sp.get();
  node_list.push_back(move(sp));
  unary_op_node_map.try_emplace({ arg, op_code, arg_exp_info_set, param1_symb_id, param2_symb_id, adl_param_name, adl_lags }, p);
  return p;
}

inline expr_t
DataTree::AddBinaryOp(expr_t arg1, BinaryOpcode op_code, expr_t arg2, int powerDerivOrder)
{
  if (auto it = binary_op_node_map.find({ arg1, arg2, op_code, powerDerivOrder });
      it != binary_op_node_map.end())
    return it->second;

  // Try to reduce to a constant
  try
    {
      double argval1 = arg1->eval({});
      double argval2 = arg2->eval({});
      double val = BinaryOpNode::eval_opcode(argval1, op_code, argval2, powerDerivOrder);
      return AddPossiblyNegativeConstant(val);
    }
  catch (ExprNode::EvalException &e)
    {
    }

  auto sp = make_unique<BinaryOpNode>(*this, node_list.size(), arg1, op_code, arg2, powerDerivOrder);
  auto p = sp.get();
  node_list.push_back(move(sp));
  binary_op_node_map.try_emplace({ arg1, arg2, op_code, powerDerivOrder }, p);
  return p;
}

inline expr_t
DataTree::AddTrinaryOp(expr_t arg1, TrinaryOpcode op_code, expr_t arg2, expr_t arg3)
{
  if (auto it = trinary_op_node_map.find({ arg1, arg2, arg3, op_code });
      it != trinary_op_node_map.end())
    return it->second;

  // Try to reduce to a constant
  try
    {
      double argval1 = arg1->eval({});
      double argval2 = arg2->eval({});
      double argval3 = arg3->eval({});
      double val = TrinaryOpNode::eval_opcode(argval1, op_code, argval2, argval3);
      return AddPossiblyNegativeConstant(val);
    }
  catch (ExprNode::EvalException &e)
    {
    }

  auto sp = make_unique<TrinaryOpNode>(*this, node_list.size(), arg1, op_code, arg2, arg3);
  auto p = sp.get();
  node_list.push_back(move(sp));
  trinary_op_node_map.try_emplace({ arg1, arg2, arg3, op_code }, p);
  return p;
}

#endif
