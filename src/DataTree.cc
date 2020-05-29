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

#include <cstdlib>
#include <cassert>
#include <iostream>
#include <regex>

#include <filesystem>

#include "DataTree.hh"

void
DataTree::initConstants()
{
  Zero = AddNonNegativeConstant("0");
  One = AddNonNegativeConstant("1");
  Two = AddNonNegativeConstant("2");
  Three = AddNonNegativeConstant("3");

  MinusOne = AddUMinus(One);

  NaN = AddNonNegativeConstant("NaN");
  Infinity = AddNonNegativeConstant("Inf");
  MinusInfinity = AddUMinus(Infinity);

  Pi = AddNonNegativeConstant("3.141592653589793");
}

DataTree::DataTree(SymbolTable &symbol_table_arg,
                   NumericalConstants &num_constants_arg,
                   ExternalFunctionsTable &external_functions_table_arg,
                   bool is_dynamic_arg) :
  symbol_table{symbol_table_arg},
  num_constants{num_constants_arg},
  external_functions_table{external_functions_table_arg},
  is_dynamic{is_dynamic_arg}
{
  initConstants();
}

DataTree::DataTree(const DataTree &d) :
  symbol_table{d.symbol_table},
  num_constants{d.num_constants},
  external_functions_table{d.external_functions_table},
  is_dynamic{d.is_dynamic},
  local_variables_vector{d.local_variables_vector}
{
  // Constants must be initialized first because they are used in some Add* methods
  initConstants();

  for (const auto &it : d.node_list)
    it->clone(*this);

  assert(node_list.size() == d.node_list.size());

  for (const auto &it : d.local_variables_table)
    local_variables_table[it.first] = it.second->clone(*this);
}

DataTree &
DataTree::operator=(const DataTree &d)
{
  assert(&symbol_table == &d.symbol_table);
  assert(&num_constants == &d.num_constants);
  assert(&external_functions_table == &d.external_functions_table);
  assert(is_dynamic == d.is_dynamic);

  num_const_node_map.clear();
  variable_node_map.clear();
  unary_op_node_map.clear();
  binary_op_node_map.clear();
  trinary_op_node_map.clear();
  external_function_node_map.clear();
  var_expectation_node_map.clear();
  pac_expectation_node_map.clear();
  first_deriv_external_function_node_map.clear();
  second_deriv_external_function_node_map.clear();

  node_list.clear();

  // Constants must be initialized first because they are used in some Add* methods
  initConstants();

  for (const auto &it : d.node_list)
    it->clone(*this);

  assert(node_list.size() == d.node_list.size());

  local_variables_vector = d.local_variables_vector;

  for (const auto &it : d.local_variables_table)
    local_variables_table[it.first] = it.second->clone(*this);

  return *this;
}

expr_t
DataTree::AddNonNegativeConstant(const string &value)
{
  int id = num_constants.AddNonNegativeConstant(value);

  if (auto it = num_const_node_map.find(id);
      it != num_const_node_map.end())
    return it->second;

  auto sp = make_unique<NumConstNode>(*this, node_list.size(), id);
  auto p = sp.get();
  node_list.push_back(move(sp));
  num_const_node_map[id] = p;
  return p;
}

VariableNode *
DataTree::AddVariable(int symb_id, int lag)
{
  if (lag != 0 && !is_dynamic)
    {
      cerr << "Leads/lags not authorized in this DataTree" << endl;
      exit(EXIT_FAILURE);
    }

  if (auto it = variable_node_map.find({ symb_id, lag });
      it != variable_node_map.end())
    return it->second;

  auto sp = make_unique<VariableNode>(*this, node_list.size(), symb_id, lag);
  auto p = sp.get();
  node_list.push_back(move(sp));
  variable_node_map[{ symb_id, lag }] = p;
  return p;
}

VariableNode *
DataTree::getVariable(int symb_id, int lag) const
{
  auto it = variable_node_map.find({ symb_id, lag });
  if (it == variable_node_map.end())
    {
      cerr << "DataTree::getVariable: unknown variable node for symb_id=" << symb_id << " and lag=" << lag << endl;
      exit(EXIT_FAILURE);
    }
  return it->second;
}

bool
DataTree::ParamUsedWithLeadLagInternal() const
{
  for (const auto &it : variable_node_map)
    if (symbol_table.getType(it.first.first) == SymbolType::parameter && it.first.second != 0)
      return true;
  return false;
}

expr_t
DataTree::AddPlus(expr_t iArg1, expr_t iArg2)
{
  if (iArg2 == Zero)
    return iArg1;

  if (iArg1 == Zero)
    return iArg2;

  // Simplify x+(-y) in x-y
  if (auto uarg2 = dynamic_cast<UnaryOpNode *>(iArg2);
      uarg2 && uarg2->op_code == UnaryOpcode::uminus)
    return AddMinus(iArg1, uarg2->arg);

  // Simplify (-x)+y in y-x
  if (auto uarg1 = dynamic_cast<UnaryOpNode *>(iArg1);
      uarg1 && uarg1->op_code == UnaryOpcode::uminus)
    return AddMinus(iArg2, uarg1->arg);

  // Simplify (x-y)+y in x
  if (auto barg1 = dynamic_cast<BinaryOpNode *>(iArg1);
      barg1 && barg1->op_code == BinaryOpcode::minus && barg1->arg2 == iArg2)
    return barg1->arg1;

  // Simplify y+(x-y) in x
  if (auto barg2 = dynamic_cast<BinaryOpNode *>(iArg2);
      barg2 && barg2->op_code == BinaryOpcode::minus && barg2->arg2 == iArg1)
    return barg2->arg1;

  // To treat commutativity of "+"
  // Nodes iArg1 and iArg2 are sorted by index
  if (iArg1->idx > iArg2->idx)
    swap(iArg1, iArg2);
  return AddBinaryOp(iArg1, BinaryOpcode::plus, iArg2);
}

expr_t
DataTree::AddMinus(expr_t iArg1, expr_t iArg2)
{
  if (iArg2 == Zero)
    return iArg1;

  if (iArg1 == Zero)
    return AddUMinus(iArg2);

  if (iArg1 == iArg2)
    return Zero;

  // Simplify x-(-y) in x+y
  if (auto uarg2 = dynamic_cast<UnaryOpNode *>(iArg2);
      uarg2 && uarg2->op_code == UnaryOpcode::uminus)
    return AddPlus(iArg1, uarg2->arg);

  // Simplify (x+y)-y and (y+x)-y in x
  if (auto barg1 = dynamic_cast<BinaryOpNode *>(iArg1);
      barg1 && barg1->op_code == BinaryOpcode::plus)
    {
      if (barg1->arg2 == iArg2)
        return barg1->arg1;
      if (barg1->arg1 == iArg2)
        return barg1->arg2;
    }

  return AddBinaryOp(iArg1, BinaryOpcode::minus, iArg2);
}

expr_t
DataTree::AddUMinus(expr_t iArg1)
{
  if (iArg1 == Zero)
    return Zero;

  // Simplify -(-x) in x
  if (auto uarg = dynamic_cast<UnaryOpNode *>(iArg1);
      uarg && uarg->op_code == UnaryOpcode::uminus)
    return uarg->arg;

  return AddUnaryOp(UnaryOpcode::uminus, iArg1);
}

expr_t
DataTree::AddTimes(expr_t iArg1, expr_t iArg2)
{
  if (iArg1 == Zero || iArg2 == Zero)
    return Zero;

  if (iArg1 == One)
    return iArg2;

  if (iArg2 == One)
    return iArg1;

  if (iArg1 == MinusOne)
    return AddUMinus(iArg2);

  if (iArg2 == MinusOne)
    return AddUMinus(iArg1);

  // Simplify (x/y)*y in x
  if (auto barg1 = dynamic_cast<BinaryOpNode *>(iArg1);
      barg1 && barg1->op_code == BinaryOpcode::divide && barg1->arg2 == iArg2)
    return barg1->arg1;

  // Simplify y*(x/y) in x
  if (auto barg2 = dynamic_cast<BinaryOpNode *>(iArg2);
      barg2 && barg2->op_code == BinaryOpcode::divide && barg2->arg2 == iArg1)
    return barg2->arg1;

  // To treat commutativity of "*"
  // Nodes iArg1 and iArg2 are sorted by index
  if (iArg1->idx > iArg2->idx)
    swap(iArg1, iArg2);
  return AddBinaryOp(iArg1, BinaryOpcode::times, iArg2);
}

expr_t
DataTree::AddDivide(expr_t iArg1, expr_t iArg2) noexcept(false)
{
  if (iArg2 == One)
    return iArg1;

  // This test should be before the next two, otherwise 0/0 won't be rejected
  if (iArg2 == Zero)
    {
      cerr << "ERROR: Division by zero!" << endl;
      throw DivisionByZeroException();
    }

  if (iArg1 == Zero)
    return Zero;

  if (iArg1 == iArg2)
    return One;

  // Simplify x/(1/y) in x*y
  if (auto barg2 = dynamic_cast<BinaryOpNode *>(iArg2);
      barg2 && barg2->op_code == BinaryOpcode::divide && barg2->arg1 == One)
    return AddTimes(iArg1, barg2->arg2);

  // Simplify (x*y)/y and (y*x)/y in x
  if (auto barg1 = dynamic_cast<BinaryOpNode *>(iArg1);
      barg1 && barg1->op_code == BinaryOpcode::times)
    {
      if (barg1->arg2 == iArg2)
        return barg1->arg1;
      if (barg1->arg1 == iArg2)
        return barg1->arg2;
    }

  return AddBinaryOp(iArg1, BinaryOpcode::divide, iArg2);
}

expr_t
DataTree::AddLess(expr_t iArg1, expr_t iArg2)
{
  return AddBinaryOp(iArg1, BinaryOpcode::less, iArg2);
}

expr_t
DataTree::AddGreater(expr_t iArg1, expr_t iArg2)
{
  return AddBinaryOp(iArg1, BinaryOpcode::greater, iArg2);
}

expr_t
DataTree::AddLessEqual(expr_t iArg1, expr_t iArg2)
{
  return AddBinaryOp(iArg1, BinaryOpcode::lessEqual, iArg2);
}

expr_t
DataTree::AddGreaterEqual(expr_t iArg1, expr_t iArg2)
{
  return AddBinaryOp(iArg1, BinaryOpcode::greaterEqual, iArg2);
}

expr_t
DataTree::AddEqualEqual(expr_t iArg1, expr_t iArg2)
{
  return AddBinaryOp(iArg1, BinaryOpcode::equalEqual, iArg2);
}

expr_t
DataTree::AddDifferent(expr_t iArg1, expr_t iArg2)
{
  return AddBinaryOp(iArg1, BinaryOpcode::different, iArg2);
}

expr_t
DataTree::AddPower(expr_t iArg1, expr_t iArg2)
{
  // This one comes first, because 0⁰=1
  if (iArg2 == Zero)
    return One;

  if (iArg1 == Zero)
    return Zero;

  if (iArg1 == One)
    return One;

  if (iArg2 == One)
    return iArg1;

  return AddBinaryOp(iArg1, BinaryOpcode::power, iArg2);
}

expr_t
DataTree::AddPowerDeriv(expr_t iArg1, expr_t iArg2, int powerDerivOrder)
{
  assert(powerDerivOrder > 0);
  return AddBinaryOp(iArg1, BinaryOpcode::powerDeriv, iArg2, powerDerivOrder);
}

expr_t
DataTree::AddDiff(expr_t iArg1)
{
  if (iArg1->maxLead() > 0)
    // Issue preprocessor#21: always expand diffs with leads
    return AddMinus(iArg1, iArg1->decreaseLeadsLags(1));
  return AddUnaryOp(UnaryOpcode::diff, iArg1);
}

expr_t
DataTree::AddAdl(expr_t iArg1, const string &name, const vector<int> &lags)
{
  return AddUnaryOp(UnaryOpcode::adl, iArg1, 0, 0, 0, string(name), lags);
}

expr_t
DataTree::AddExp(expr_t iArg1)
{
  if (iArg1 == Zero)
    return One;

  return AddUnaryOp(UnaryOpcode::exp, iArg1);
}

expr_t
DataTree::AddLog(expr_t iArg1)
{
  if (iArg1 == One)
    return Zero;

  if (iArg1 == Zero)
    {
      cerr << "ERROR: log(0) not defined!" << endl;
      exit(EXIT_FAILURE);
    }

  // Simplify log(1/x) in −log(x)
  if (auto barg1 = dynamic_cast<BinaryOpNode *>(iArg1);
      barg1 && barg1->op_code == BinaryOpcode::divide && barg1->arg1 == One)
    return AddUMinus(AddLog(barg1->arg2));

  return AddUnaryOp(UnaryOpcode::log, iArg1);
}

expr_t
DataTree::AddLog10(expr_t iArg1)
{
  if (iArg1 == One)
    return Zero;

  if (iArg1 == Zero)
    {
      cerr << "ERROR: log10(0) not defined!" << endl;
      exit(EXIT_FAILURE);
    }

  // Simplify log₁₀(1/x) in −log₁₀(x)
  if (auto barg1 = dynamic_cast<BinaryOpNode *>(iArg1);
      barg1 && barg1->op_code == BinaryOpcode::divide && barg1->arg1 == One)
    return AddUMinus(AddLog10(barg1->arg2));

  return AddUnaryOp(UnaryOpcode::log10, iArg1);
}

expr_t
DataTree::AddCos(expr_t iArg1)
{
  if (iArg1 == Zero)
    return One;

  return AddUnaryOp(UnaryOpcode::cos, iArg1);
}

expr_t
DataTree::AddSin(expr_t iArg1)
{
  if (iArg1 == Zero)
    return Zero;

  return AddUnaryOp(UnaryOpcode::sin, iArg1);
}

expr_t
DataTree::AddTan(expr_t iArg1)
{
  if (iArg1 == Zero)
    return Zero;

  return AddUnaryOp(UnaryOpcode::tan, iArg1);
}

expr_t
DataTree::AddAcos(expr_t iArg1)
{
  if (iArg1 == One)
    return Zero;

  return AddUnaryOp(UnaryOpcode::acos, iArg1);
}

expr_t
DataTree::AddAsin(expr_t iArg1)
{
  if (iArg1 == Zero)
    return Zero;

  return AddUnaryOp(UnaryOpcode::asin, iArg1);
}

expr_t
DataTree::AddAtan(expr_t iArg1)
{
  if (iArg1 == Zero)
    return Zero;

  return AddUnaryOp(UnaryOpcode::atan, iArg1);
}

expr_t
DataTree::AddCosh(expr_t iArg1)
{
  if (iArg1 == Zero)
    return One;

  return AddUnaryOp(UnaryOpcode::cosh, iArg1);
}

expr_t
DataTree::AddSinh(expr_t iArg1)
{
  if (iArg1 == Zero)
    return Zero;

  return AddUnaryOp(UnaryOpcode::sinh, iArg1);
}

expr_t
DataTree::AddTanh(expr_t iArg1)
{
  if (iArg1 == Zero)
    return Zero;

  return AddUnaryOp(UnaryOpcode::tanh, iArg1);
}

expr_t
DataTree::AddAcosh(expr_t iArg1)
{
  if (iArg1 == One)
    return Zero;

  return AddUnaryOp(UnaryOpcode::acosh, iArg1);
}

expr_t
DataTree::AddAsinh(expr_t iArg1)
{
  if (iArg1 == Zero)
    return Zero;

  return AddUnaryOp(UnaryOpcode::asinh, iArg1);
}

expr_t
DataTree::AddAtanh(expr_t iArg1)
{
  if (iArg1 == Zero)
    return Zero;

  return AddUnaryOp(UnaryOpcode::atanh, iArg1);
}

expr_t
DataTree::AddSqrt(expr_t iArg1)
{
  if (iArg1 == Zero)
    return Zero;

  if (iArg1 == One)
    return One;

  return AddUnaryOp(UnaryOpcode::sqrt, iArg1);
}

expr_t
DataTree::AddCbrt(expr_t iArg1)
{
  if (iArg1 == Zero)
    return Zero;

  if (iArg1 == One)
    return One;

  return AddUnaryOp(UnaryOpcode::cbrt, iArg1);
}

expr_t
DataTree::AddAbs(expr_t iArg1)
{
  if (iArg1 == Zero)
    return Zero;

  if (iArg1 == One)
    return One;

  return AddUnaryOp(UnaryOpcode::abs, iArg1);
}

expr_t
DataTree::AddSign(expr_t iArg1)
{
  if (iArg1 == Zero)
    return Zero;

  if (iArg1 == One)
    return One;

  return AddUnaryOp(UnaryOpcode::sign, iArg1);
}

expr_t
DataTree::AddErf(expr_t iArg1)
{
  if (iArg1 == Zero)
    return Zero;

  return AddUnaryOp(UnaryOpcode::erf, iArg1);
}

expr_t
DataTree::AddMax(expr_t iArg1, expr_t iArg2)
{
  return AddBinaryOp(iArg1, BinaryOpcode::max, iArg2);
}

expr_t
DataTree::AddMin(expr_t iArg1, expr_t iArg2)
{
  return AddBinaryOp(iArg1, BinaryOpcode::min, iArg2);
}

expr_t
DataTree::AddNormcdf(expr_t iArg1, expr_t iArg2, expr_t iArg3)
{
  return AddTrinaryOp(iArg1, TrinaryOpcode::normcdf, iArg2, iArg3);
}

expr_t
DataTree::AddNormpdf(expr_t iArg1, expr_t iArg2, expr_t iArg3)
{
  return AddTrinaryOp(iArg1, TrinaryOpcode::normpdf, iArg2, iArg3);
}

expr_t
DataTree::AddSteadyState(expr_t iArg1)
{
  return AddUnaryOp(UnaryOpcode::steadyState, iArg1);
}

expr_t
DataTree::AddSteadyStateParamDeriv(expr_t iArg1, int param_symb_id)
{
  return AddUnaryOp(UnaryOpcode::steadyStateParamDeriv, iArg1, 0, param_symb_id);
}

expr_t
DataTree::AddSteadyStateParam2ndDeriv(expr_t iArg1, int param1_symb_id, int param2_symb_id)
{
  return AddUnaryOp(UnaryOpcode::steadyStateParam2ndDeriv, iArg1, 0, param1_symb_id, param2_symb_id);
}

expr_t
DataTree::AddExpectation(int iArg1, expr_t iArg2)
{
  return AddUnaryOp(UnaryOpcode::expectation, iArg2, iArg1);
}

expr_t
DataTree::AddVarExpectation(const string &model_name)
{
  if (auto it = var_expectation_node_map.find(model_name);
      it != var_expectation_node_map.end())
    return it->second;

  auto sp = make_unique<VarExpectationNode>(*this, node_list.size(), model_name);
  auto p = sp.get();
  node_list.push_back(move(sp));
  var_expectation_node_map[model_name] = p;
  return p;
}

expr_t
DataTree::AddPacExpectation(const string &model_name)
{
  if (auto it = pac_expectation_node_map.find(model_name);
      it != pac_expectation_node_map.end())
    return it->second;

  auto sp = make_unique<PacExpectationNode>(*this, node_list.size(), model_name);
  auto p = sp.get();
  node_list.push_back(move(sp));
  pac_expectation_node_map[model_name] = p;
  return p;
}

BinaryOpNode *
DataTree::AddEqual(expr_t iArg1, expr_t iArg2)
{
  /* We know that we can safely cast to BinaryOpNode because
     BinaryOpCode::equal can never be reduced to a constant. */
  return dynamic_cast<BinaryOpNode *>(AddBinaryOp(iArg1, BinaryOpcode::equal, iArg2));
}

void
DataTree::AddLocalVariable(int symb_id, expr_t value) noexcept(false)
{
  assert(symbol_table.getType(symb_id) == SymbolType::modelLocalVariable);

  // Throw an exception if symbol already declared
  if (auto it = local_variables_table.find(symb_id);
      it != local_variables_table.end())
    throw LocalVariableException(symbol_table.getName(symb_id));

  local_variables_table[symb_id] = value;
  local_variables_vector.push_back(symb_id);
}

expr_t
DataTree::AddExternalFunction(int symb_id, const vector<expr_t> &arguments)
{
  assert(symbol_table.getType(symb_id) == SymbolType::externalFunction);

  if (auto it = external_function_node_map.find({ arguments, symb_id });
      it != external_function_node_map.end())
    return it->second;

  auto sp = make_unique<ExternalFunctionNode>(*this, node_list.size(), symb_id, arguments);
  auto p = sp.get();
  node_list.push_back(move(sp));
  external_function_node_map[{ arguments, symb_id }] = p;
  return p;
}

expr_t
DataTree::AddFirstDerivExternalFunction(int top_level_symb_id, const vector<expr_t> &arguments, int input_index)
{
  assert(symbol_table.getType(top_level_symb_id) == SymbolType::externalFunction);

  if (auto it = first_deriv_external_function_node_map.find({ arguments, input_index, top_level_symb_id });
      it != first_deriv_external_function_node_map.end())
    return it->second;

  auto sp = make_unique<FirstDerivExternalFunctionNode>(*this, node_list.size(), top_level_symb_id, arguments, input_index);
  auto p = sp.get();
  node_list.push_back(move(sp));
  first_deriv_external_function_node_map[{ arguments, input_index, top_level_symb_id }] = p;
  return p;
}

expr_t
DataTree::AddSecondDerivExternalFunction(int top_level_symb_id, const vector<expr_t> &arguments, int input_index1, int input_index2)
{
  assert(symbol_table.getType(top_level_symb_id) == SymbolType::externalFunction);

  if (auto it = second_deriv_external_function_node_map.find({ arguments, input_index1, input_index2,
                                                               top_level_symb_id });
    it != second_deriv_external_function_node_map.end())
    return it->second;

  auto sp = make_unique<SecondDerivExternalFunctionNode>(*this, node_list.size(), top_level_symb_id, arguments, input_index1, input_index2);
  auto p = sp.get();
  node_list.push_back(move(sp));
  second_deriv_external_function_node_map[{ arguments, input_index1, input_index2, top_level_symb_id }] = p;
  return p;
}

bool
DataTree::isSymbolUsed(int symb_id) const
{
  for (const auto &it : variable_node_map)
    if (it.first.first == symb_id)
      return true;

  if (local_variables_table.find(symb_id) != local_variables_table.end())
    return true;

  return false;
}

int
DataTree::getDerivID(int symb_id, int lag) const noexcept(false)
{
  throw UnknownDerivIDException();
}

SymbolType
DataTree::getTypeByDerivID(int deriv_id) const noexcept(false)
{
  throw UnknownDerivIDException();
}

int
DataTree::getLagByDerivID(int deriv_id) const noexcept(false)
{
  throw UnknownDerivIDException();
}

int
DataTree::getSymbIDByDerivID(int deriv_id) const noexcept(false)
{
  throw UnknownDerivIDException();
}

void
DataTree::addAllParamDerivId(set<int> &deriv_id_set)
{
}

int
DataTree::getDynJacobianCol(int deriv_id) const noexcept(false)
{
  throw UnknownDerivIDException();
}

bool
DataTree::isUnaryOpUsed(UnaryOpcode opcode) const
{
  for (const auto &it : unary_op_node_map)
    if (get<1>(it.first) == opcode)
      return true;

  return false;
}

bool
DataTree::isUnaryOpUsedOnType(SymbolType type, UnaryOpcode opcode) const
{
  set<int> var;
  for (const auto &it : unary_op_node_map)
    if (get<1>(it.first) == opcode)
      {
        it.second->collectVariables(type, var);
        if (!var.empty())
          return true;
      }
  return false;
}

bool
DataTree::isBinaryOpUsed(BinaryOpcode opcode) const
{
  for (const auto &it : binary_op_node_map)
    if (get<2>(it.first) == opcode)
      return true;

  return false;
}

bool
DataTree::isBinaryOpUsedOnType(SymbolType type, BinaryOpcode opcode) const
{
  set<int> var;
  for (const auto &it : binary_op_node_map)
    if (get<2>(it.first) == opcode)
      {
        it.second->collectVariables(type, var);
        if (!var.empty())
          return true;
      }
  return false;
}

bool
DataTree::isTrinaryOpUsed(TrinaryOpcode opcode) const
{
  for (const auto &it : trinary_op_node_map)
    if (get<3>(it.first) == opcode)
      return true;

  return false;
}

bool
DataTree::isExternalFunctionUsed(int symb_id) const
{
  for (const auto &it : external_function_node_map)
    if (it.first.second == symb_id)
      return true;

  return false;
}

bool
DataTree::isFirstDerivExternalFunctionUsed(int symb_id) const
{
  for (const auto &it : first_deriv_external_function_node_map)
    if (get<2>(it.first) == symb_id)
      return true;

  return false;
}

bool
DataTree::isSecondDerivExternalFunctionUsed(int symb_id) const
{
  for (const auto &it : second_deriv_external_function_node_map)
    if (get<3>(it.first) == symb_id)
      return true;

  return false;
}

int
DataTree::minLagForSymbol(int symb_id) const
{
  int r = 0;
  for (const auto &it : variable_node_map)
    if (it.first.first == symb_id && it.first.second < r)
      r = it.first.second;
  return r;
}

void
DataTree::writePowerDeriv(ostream &output) const
{
  if (isBinaryOpUsed(BinaryOpcode::powerDeriv))
    output << "/*" << endl
           << " * The k-th derivative of x^p" << endl
           << " */" << endl
           << "double getPowerDeriv(double x, double p, int k)" << endl
           << "{" << endl
           << "#ifdef _MSC_VER" << endl
           << "# define nearbyint(x) (fabs((x)-floor(x)) < fabs((x)-ceil(x)) ? floor(x) : ceil(x))" << endl
           << "#endif" << endl
           << "  if ( fabs(x) < " << near_zero << " && p > 0 && k > p && fabs(p-nearbyint(p)) < " << near_zero << " )" << endl
           << "    return 0.0;" << endl
           << "  else" << endl
           << "    {" << endl
           << "      int i = 0;" << endl
           << "      double dxp = pow(x, p-k);" << endl
           << "      for (; i<k; i++)" << endl
           << "        dxp *= p--;" << endl
           << "      return dxp;" << endl
           << "    }" << endl
           << "}" << endl;
}

string
DataTree::packageDir(const string &package)
{
  regex pat{R"(\.)"};
  string dirname = "+" + regex_replace(package, pat, "/+");
  filesystem::create_directories(dirname);
  return dirname;
}
