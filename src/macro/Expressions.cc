/*
 * Copyright © 2019-2023 Dynare Team
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

#include <utility>
#include <numbers>

#include "Expressions.hh"

using namespace macro;

BoolPtr
BaseType::is_different(const BaseTypePtr &btp) const
{
  if (*(this->is_equal(btp)))
    return make_shared<Bool>(false);
  return make_shared<Bool>(true);
}

BoolPtr
Bool::is_equal(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<Bool>(btp);
  if (!btp2)
    return make_shared<Bool>(false);
  return make_shared<Bool>(value == btp2->value);
}

BoolPtr
Bool::logical_and(const ExpressionPtr &ep, Environment &env) const
{
  if (!value)
    return make_shared<Bool>(false);

  auto btp = ep->eval(env);
  if (auto btp2 = dynamic_pointer_cast<Bool>(btp); btp2)
    return make_shared<Bool>(*btp2);

  if (auto btp2 = dynamic_pointer_cast<Real>(btp); btp2)
    return make_shared<Bool>(*btp2);

  throw StackTrace("Type mismatch for operands of && operator");
}

BoolPtr
Bool::logical_or(const ExpressionPtr &ep, Environment &env) const
{
  if (value)
    return make_shared<Bool>(true);

  auto btp = ep->eval(env);
  if (auto btp2 = dynamic_pointer_cast<Bool>(btp); btp2)
    return make_shared<Bool>(*btp2);

  if (auto btp2 = dynamic_pointer_cast<Real>(btp); btp2)
    return make_shared<Bool>(*btp2);

  throw StackTrace("Type mismatch for operands of || operator");
}

BoolPtr
Bool::logical_not() const
{
  return make_shared<Bool>(!value);
}

BaseTypePtr
Real::plus(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<Real>(btp);
  if (!btp2)
    throw StackTrace("Type mismatch for operands of + operator");
  return make_shared<Real>(value + btp2->value);
}

BaseTypePtr
Real::minus(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<Real>(btp);
  if (!btp2)
    throw StackTrace("Type mismatch for operands of - operator");
  return make_shared<Real>(value - btp2->value);
}

BaseTypePtr
Real::times(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<Real>(btp);
  if (!btp2)
    throw StackTrace("Type mismatch for operands of * operator");
  return make_shared<Real>(value * btp2->value);
}

BaseTypePtr
Real::divide(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<Real>(btp);
  if (!btp2)
    throw StackTrace("Type mismatch for operands of / operator");
  return make_shared<Real>(value / btp2->value);
}

BaseTypePtr
Real::power(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<Real>(btp);
  if (!btp2)
    throw StackTrace("Type mismatch for operands of ^ operator");
  return make_shared<Real>(pow(value, btp2->value));
}

BoolPtr
Real::is_less(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<Real>(btp);
  if (!btp2)
    throw StackTrace("Type mismatch for operands of < operator");
  return make_shared<Bool>(isless(value, btp2->value));
}

BoolPtr
Real::is_greater(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<Real>(btp);
  if (!btp2)
    throw StackTrace("Type mismatch for operands of > operator");
  return make_shared<Bool>(isgreater(value, btp2->value));
}

BoolPtr
Real::is_less_equal(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<Real>(btp);
  if (!btp2)
    throw StackTrace("Type mismatch for operands of <= operator");
  return make_shared<Bool>(islessequal(value, btp2->value));
}

BoolPtr
Real::is_greater_equal(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<Real>(btp);
  if (!btp2)
    throw StackTrace("Type mismatch for operands of >= operator");
  return make_shared<Bool>(isgreaterequal(value, btp2->value));
}

BoolPtr
Real::is_equal(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<Real>(btp);
  if (!btp2)
    return make_shared<Bool>(false);
  return make_shared<Bool>(value == btp2->value);
}

BoolPtr
Real::logical_and(const ExpressionPtr &ep, Environment &env) const
{
  if (!value)
    return make_shared<Bool>(false);

  auto btp = ep->eval(env);
  if (auto btp2 = dynamic_pointer_cast<Real>(btp); btp2)
    return make_shared<Bool>(*btp2);

  if (auto btp2 = dynamic_pointer_cast<Bool>(btp); btp2)
    return make_shared<Bool>(*btp2);

  throw StackTrace("Type mismatch for operands of && operator");
}

BoolPtr
Real::logical_or(const ExpressionPtr &ep, Environment &env) const
{
  if (value)
    return make_shared<Bool>(true);

  auto btp = ep->eval(env);
  if (auto btp2 = dynamic_pointer_cast<Real>(btp); btp2)
    return make_shared<Bool>(*btp2);

  if (auto btp2 = dynamic_pointer_cast<Bool>(btp); btp2)
    return make_shared<Bool>(*btp2);

  throw StackTrace("Type mismatch for operands of || operator");
}

BoolPtr
Real::logical_not() const
{
  return make_shared<Bool>(!value);
}

RealPtr
Real::max(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<Real>(btp);
  if (!btp2)
    throw StackTrace("Type mismatch for operands of `max` operator");
  return make_shared<Real>(std::max(value, btp2->value));
}

RealPtr
Real::min(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<Real>(btp);
  if (!btp2)
    throw StackTrace("Type mismatch for operands of `min` operator");
  return make_shared<Real>(std::min(value, btp2->value));
}

RealPtr
Real::mod(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<Real>(btp);
  if (!btp2)
    throw StackTrace("Type mismatch for operands of `mod` operator");
  return make_shared<Real>(std::fmod(value, btp2->value));
}

RealPtr
Real::normpdf(const BaseTypePtr &btp1, const BaseTypePtr &btp2) const
{
  auto btp12 = dynamic_pointer_cast<Real>(btp1);
  auto btp22 = dynamic_pointer_cast<Real>(btp2);
  if (!btp12 || !btp22)
    throw StackTrace("Type mismatch for operands of `normpdf` operator");
  return make_shared<Real>((1/(btp22->value*std::sqrt(2*numbers::pi)*std::exp(pow((value-btp12->value)/btp22->value, 2)/2))));
}

RealPtr
Real::normcdf(const BaseTypePtr &btp1, const BaseTypePtr &btp2) const
{
  auto btp12 = dynamic_pointer_cast<Real>(btp1);
  auto btp22 = dynamic_pointer_cast<Real>(btp2);
  if (!btp12 || !btp22)
    throw StackTrace("Type mismatch for operands of `normpdf` operator");
  return make_shared<Real>((0.5*(1+std::erf((value-btp12->value)/btp22->value/numbers::sqrt2))));
}

BaseTypePtr
String::plus(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<String>(btp);
  if (!btp2)
    throw StackTrace("Type mismatch for operands of + operator");
  return make_shared<String>(value + btp2->value);
}

BoolPtr
String::is_less(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<String>(btp);
  if (!btp2)
    throw StackTrace("Type mismatch for operands of < operator");
  return make_shared<Bool>(value < btp2->value);
}

BoolPtr
String::is_greater(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<String>(btp);
  if (!btp2)
    throw StackTrace("Type mismatch for operands of > operator");
  return make_shared<Bool>(value > btp2->value);
}

BoolPtr
String::is_less_equal(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<String>(btp);
  if (!btp2)
    throw StackTrace("Type mismatch for operands of <= operator");
  return make_shared<Bool>(value <= btp2->value);
}

BoolPtr
String::is_greater_equal(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<String>(btp);
  if (!btp2)
    throw StackTrace("Type mismatch for operands of >= operator");
  return make_shared<Bool>(value >= btp2->value);
}

BoolPtr
String::is_equal(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<String>(btp);
  if (!btp2)
    return make_shared<Bool>(false);
  return make_shared<Bool>(value == btp2->value);
}

BoolPtr
String::cast_bool([[maybe_unused]] Environment &env) const
{
  auto f = [](const char &a, const char &b) { return (tolower(a) == tolower(b)); };

  if (string tf = "true"; equal(value.begin(), value.end(), tf.begin(), tf.end(), f))
    return make_shared<Bool>(true);

  if (string tf = "false"; equal(value.begin(), value.end(), tf.begin(), tf.end(), f))
    return make_shared<Bool>(false);

  try
    {
      size_t pos = 0;
      double value_d = stod(value, &pos);
      if (pos != value.length())
        throw StackTrace("Entire string not converted");
      return make_shared<Bool>(static_cast<bool>(value_d));
    }
  catch (...)
    {
      throw StackTrace(R"(")" + value + R"(" cannot be converted to a boolean)");
    }
}

RealPtr
String::cast_real([[maybe_unused]] Environment &env) const
{
  try
    {
      size_t pos = 0;
      double value_d = stod(value, &pos);
      if (pos != value.length())
        throw StackTrace("Entire string not converted");
      return make_shared<Real>(value_d);
    }
  catch (...)
    {
      throw StackTrace(R"(")" + value + R"(" cannot be converted to a real)");
    }
}

BaseTypePtr
Array::plus(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<Array>(btp);
  if (!btp2)
    throw StackTrace("Type mismatch for operands of + operator");

  vector<ExpressionPtr> arr_copy{arr};
  arr_copy.insert(arr_copy.end(), btp2->arr.begin(), btp2->arr.end());
  return make_shared<Array>(arr_copy);
}

BaseTypePtr
Array::minus(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<Array>(btp);
  if (!btp2)
    throw StackTrace("Type mismatch for operands of - operator");

  /* Highly inefficient algorithm for computing set difference
     (but vector<T> is not suited for that...) */
  vector<ExpressionPtr> arr_copy;
  for (const auto &it : arr)
    {
      auto itbtp = dynamic_pointer_cast<BaseType>(it);
      auto it2 = btp2->arr.cbegin();
      for (; it2 != btp2->arr.cend(); ++it2)
        if (*(itbtp->is_equal(dynamic_pointer_cast<BaseType>(*it2))))
          break;
      if (it2 == btp2->arr.cend())
        arr_copy.emplace_back(itbtp);
    }
  return make_shared<Array>(arr_copy);
}

BaseTypePtr
Array::times(const BaseTypePtr &btp) const
{
  vector<ExpressionPtr> values;
  auto btp2 = dynamic_pointer_cast<Array>(btp);
  if (!btp2)
    throw StackTrace("Type mismatch for operands of * operator");

  for (const auto &itl : arr)
    for (const auto &itr : btp2->getValue())
      {
        vector<ExpressionPtr> new_tuple;
        if (dynamic_pointer_cast<Real>(itl) || dynamic_pointer_cast<String>(itl))
          new_tuple.push_back(itl);
        else if (dynamic_pointer_cast<Tuple>(itl))
          new_tuple = dynamic_pointer_cast<Tuple>(itl)->getValue();
        else
          throw StackTrace("Array::times: unsupported type on lhs");

        if (dynamic_pointer_cast<Real>(itr) || dynamic_pointer_cast<String>(itr))
          new_tuple.push_back(itr);
        else if (dynamic_pointer_cast<Tuple>(itr))
          for (const auto &tit : dynamic_pointer_cast<Tuple>(itr)->getValue())
            new_tuple.push_back(tit);
        else
          throw StackTrace("Array::times: unsupported type on rhs");

        values.emplace_back(make_shared<Tuple>(new_tuple));
      }

  return make_shared<Array>(values);
}

BaseTypePtr
Array::power(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<Real>(btp);
  if (!btp2 || !*(btp2->isinteger()))
    throw StackTrace("The second argument of the power operator (^) must be an integer");

  auto retval = make_shared<Array>(arr);
  for (int i = 1; i < *btp2; i++)
    {
      auto btpv = retval->times(make_shared<Array>(arr));
      retval = make_shared<Array>(dynamic_pointer_cast<Array>(btpv)->getValue());
    }
  return retval;
}

BoolPtr
Array::is_equal(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<Array>(btp);
  if (!btp2)
    return make_shared<Bool>(false);

  if (arr.size() != btp2->arr.size())
    return make_shared<Bool>(false);

  for (size_t i = 0; i < arr.size(); i++)
    {
      auto bt = dynamic_pointer_cast<BaseType>(arr[i]);
      auto bt2 = dynamic_pointer_cast<BaseType>(btp2->arr[i]);
      if (!*(bt->is_equal(bt2)))
        return make_shared<Bool>(false);
    }
  return make_shared<Bool>(true);
}

ArrayPtr
Array::set_union(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<Array>(btp);
  if (!btp2)
    throw StackTrace("Arguments of the union operator (|) must be sets");

  vector<ExpressionPtr> new_values = arr;
  for (const auto &it : btp2->arr)
    {
      bool found = false;
      auto it2 = dynamic_pointer_cast<BaseType>(it);
      if (!it2)
        throw StackTrace("Type mismatch for operands of in operator");
      for (const auto &nvit : new_values)
        {
          auto v2 = dynamic_pointer_cast<BaseType>(nvit);
          if (!v2)
            throw StackTrace("Type mismatch for operands of in operator");
          if (*(v2->is_equal(it2)))
            {
              found = true;
              break;
            }
        }
      if (!found)
        new_values.push_back(it);
    }
  return make_shared<Array>(new_values);
}

ArrayPtr
Array::set_intersection(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<Array>(btp);
  if (!btp2)
    throw StackTrace("Arguments of the intersection operator (|) must be sets");

  vector<ExpressionPtr> new_values;
  for (const auto &it : btp2->arr)
    {
      auto it2 = dynamic_pointer_cast<BaseType>(it);
      if (!it2)
        throw StackTrace("Type mismatch for operands of in operator");
      for (const auto &nvit : arr)
        {
          auto v2 = dynamic_pointer_cast<BaseType>(nvit);
          if (!v2)
            throw StackTrace("Type mismatch for operands of in operator");
          if (*(v2->is_equal(it2)))
            {
              new_values.push_back(it);
              break;
            }
        }
    }
  return make_shared<Array>(new_values);
}

BoolPtr
Array::contains(const BaseTypePtr &btp) const
{
  for (const auto &v : arr)
    {
      auto v2 = dynamic_pointer_cast<BaseType>(v);
      if (!v2)
        throw StackTrace("Type mismatch for operands of in operator");
      if (*(v2->is_equal(btp)))
        return make_shared<Bool>(true);
    }
  return make_shared<Bool>(false);
}

RealPtr
Array::sum() const
{
  double retval = 0;
  for (const auto &v : arr)
    {
      auto v2 = dynamic_pointer_cast<Real>(v);
      if (!v2)
        throw StackTrace("Type mismatch for operands of in operator");
      retval += *v2;
    }
  return make_shared<Real>(retval);
}

BoolPtr
Array::cast_bool(Environment &env) const
{
  if (arr.size() != 1)
    throw StackTrace("Array must be of size 1 to be cast to a boolean");
  return arr.at(0)->eval(env)->cast_bool(env);
}

RealPtr
Array::cast_real(Environment &env) const
{
  if (arr.size() != 1)
    throw StackTrace("Array must be of size 1 to be cast to a real");
  return arr.at(0)->eval(env)->cast_real(env);
}

BoolPtr
Tuple::is_equal(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<Tuple>(btp);
  if (!btp2)
    return make_shared<Bool>(false);

  if (tup.size() != btp2->tup.size())
    return make_shared<Bool>(false);

  for (size_t i = 0; i < tup.size(); i++)
    {
      auto bt = dynamic_pointer_cast<BaseType>(tup[i]);
      auto bt2 = dynamic_pointer_cast<BaseType>(btp2->tup[i]);
      if (!*(bt->is_equal(bt2)))
        return make_shared<Bool>(false);
    }
  return make_shared<Bool>(true);
}

BoolPtr
Tuple::contains(const BaseTypePtr &btp) const
{
  for (const auto &v : tup)
    {
      auto v2 = dynamic_pointer_cast<BaseType>(v);
      if (!v2)
        throw StackTrace("Type mismatch for operands of in operator");
      if (*(v2->is_equal(btp)))
        return make_shared<Bool>(true);
    }
  return make_shared<Bool>(false);
}

BoolPtr
Tuple::cast_bool(Environment &env) const
{
  if (tup.size() != 1)
    throw StackTrace("Tuple must be of size 1 to be cast to a boolean");
  return tup.at(0)->eval(env)->cast_bool(env);
}

RealPtr
Tuple::cast_real(Environment &env) const
{
  if (tup.size() != 1)
    throw StackTrace("Tuple must be of size 1 to be cast to a real");
  return tup.at(0)->eval(env)->cast_real(env);
}

BaseTypePtr
Range::eval(Environment &env) const
{
  RealPtr incdbl = make_shared<Real>(1);
  if (inc)
    incdbl = dynamic_pointer_cast<Real>(inc->eval(env));
  RealPtr startdbl = dynamic_pointer_cast<Real>(start->eval(env));
  RealPtr enddbl = dynamic_pointer_cast<Real>(end->eval(env));
  if (!startdbl || !enddbl || !incdbl)
    throw StackTrace("To create an array from a range using the colon operator, "
                     "the arguments must evaluate to reals");

  vector<ExpressionPtr> arr;
  if (*incdbl > 0 && *startdbl <= *enddbl)
    for (double i = *startdbl; i <= *enddbl; i += *incdbl)
      arr.emplace_back(make_shared<Real>(i));
  else if (*startdbl >= *enddbl && *incdbl < 0)
    for (double i = *startdbl; i >= *enddbl; i += *incdbl)
      arr.emplace_back(make_shared<Real>(i));

  return make_shared<Array>(arr, location);
}

BaseTypePtr
Array::eval(Environment &env) const
{
  vector<ExpressionPtr> retval;
  for (const auto &it : arr)
    retval.emplace_back(it->eval(env));
  return make_shared<Array>(retval);
}

BaseTypePtr
Tuple::eval(Environment &env) const
{
  vector<ExpressionPtr> retval;
  for (const auto &it : tup)
    retval.emplace_back(it->eval(env));
  return make_shared<Tuple>(retval);
}

BaseTypePtr
Variable::eval(Environment &env) const
{
  if (indices && !indices->empty())
    {
      ArrayPtr map = dynamic_pointer_cast<Array>(indices->eval(env));
      vector<int> ind;
      for (const auto &it : map->getValue())
        // Necessary to handle indexes like: y[1:2,2]
        // In general this evaluates to [[1:2],2] but when subscripting we want to expand it to [1,2,2]
        if (auto db = dynamic_pointer_cast<Real>(it); db)
          {
            if (!*(db->isinteger()))
              throw StackTrace("variable", "When indexing a variable you must pass "
                               "an int or an int array", location);
            ind.emplace_back(*db);
          }
        else if (dynamic_pointer_cast<Array>(it))
          for (const auto &it1 : dynamic_pointer_cast<Array>(it)->getValue())
            if (db = dynamic_pointer_cast<Real>(it1); db)
              {
                if (!*(db->isinteger()))
                  throw StackTrace("variable", "When indexing a variable you must pass "
                                   "an int or an int array", location);
                ind.emplace_back(*db);
              }
            else
              throw StackTrace("variable", "You cannot index a variable with a "
                               "nested array", location);
        else
          throw StackTrace("variable", "You can only index a variable with an int or "
                           "an int array", location);

      switch (env.getType(name))
        {
        case codes::BaseType::Bool:
          throw StackTrace("variable", "You cannot index a boolean", location);
        case codes::BaseType::Real:
          throw StackTrace("variable", "You cannot index a real", location);
        case codes::BaseType::Tuple:
          throw StackTrace("variable", "You cannot index a tuple", location);
        case codes::BaseType::Range:
          throw StackTrace("variable", "Internal Error: Range: should not arrive here", location);
        case codes::BaseType::String:
          {
            string orig_string
              = dynamic_pointer_cast<String>(env.getVariable(name))->to_string();
            string retvals;
            for (auto it : ind)
              try
                {
                  retvals += orig_string.substr(it - 1, 1);
                }
              catch (const out_of_range &ex)
                {
                  throw StackTrace("variable", "Index out of range", location);
                }
            return make_shared<String>(retvals);
          }
        case codes::BaseType::Array:
          {
            ArrayPtr ap = dynamic_pointer_cast<Array>(env.getVariable(name));
            vector<BaseTypePtr> retval;
            for (auto it : ind)
              try
                {
                  retval.emplace_back(ap->at(it - 1)->eval(env));
                }
              catch (const out_of_range &ex)
                {
                  throw StackTrace("variable", "Index out of range", location);
                }

            if (retval.size() == 1)
              return retval.at(0);
            vector<ExpressionPtr> retvala(retval.begin(), retval.end());
            return make_shared<Array>(retvala);
          }
        }
    }
  return env.getVariable(name)->eval(env);
}

BaseTypePtr
Function::eval(Environment &env) const
{
  FunctionPtr func;
  ExpressionPtr body;
  Environment env_orig = env;
  env = new Environment(env);
  try
    {
      tie(func, body) = env.getFunction(name);
    }
  catch (StackTrace &ex)
    {
      ex.push("Function", location);
      throw;
    }

  if (func->args.size() != args.size())
    throw StackTrace("Function", "The number of arguments used to call " + name
                     +" does not match the number used in its definition", location);

  try
    {
      for (size_t i = 0; i < func->args.size(); i++)
        {
          VariablePtr mvp = dynamic_pointer_cast<Variable>(func->args.at(i));
          env.define(mvp, args.at(i)->eval(env));
        }
      auto retval = body->eval(env);
      env = env_orig;
      return retval;
    }
  catch (StackTrace &ex)
    {
      ex.push("Function", location);
      throw;
    }
}

BaseTypePtr
UnaryOp::eval(Environment &env) const
{
  try
    {
      switch (op_code)
        {
        case codes::UnaryOp::cast_bool:
          return arg->eval(env)->cast_bool(env);
        case codes::UnaryOp::cast_real:
          return arg->eval(env)->cast_real(env);
        case codes::UnaryOp::cast_string:
          return arg->eval(env)->cast_string();
        case codes::UnaryOp::cast_tuple:
          return arg->eval(env)->cast_tuple();
        case codes::UnaryOp::cast_array:
          return arg->eval(env)->cast_array();
        case codes::UnaryOp::logical_not:
          return arg->eval(env)->logical_not();
        case codes::UnaryOp::unary_minus:
          return arg->eval(env)->unary_minus();
        case codes::UnaryOp::unary_plus:
          return arg->eval(env)->unary_plus();
        case codes::UnaryOp::length:
          return arg->eval(env)->length();
        case codes::UnaryOp::isempty:
          return arg->eval(env)->isempty();
        case codes::UnaryOp::isboolean:
          return arg->eval(env)->isboolean();
        case codes::UnaryOp::isreal:
          return arg->eval(env)->isreal();
        case codes::UnaryOp::isstring:
          return arg->eval(env)->isstring();
        case codes::UnaryOp::istuple:
          return arg->eval(env)->istuple();
        case codes::UnaryOp::isarray:
          return arg->eval(env)->isarray();
        case codes::UnaryOp::exp:
          return arg->eval(env)->exp();
        case codes::UnaryOp::ln:
          return arg->eval(env)->ln();
        case codes::UnaryOp::log10:
          return arg->eval(env)->log10();
        case codes::UnaryOp::sin:
          return arg->eval(env)->sin();
        case codes::UnaryOp::cos:
          return arg->eval(env)->cos();
        case codes::UnaryOp::tan:
          return arg->eval(env)->tan();
        case codes::UnaryOp::asin:
          return arg->eval(env)->asin();
        case codes::UnaryOp::acos:
          return arg->eval(env)->acos();
        case codes::UnaryOp::atan:
          return arg->eval(env)->atan();
        case codes::UnaryOp::sqrt:
          return arg->eval(env)->sqrt();
        case codes::UnaryOp::cbrt:
          return arg->eval(env)->cbrt();
        case codes::UnaryOp::sign:
          return arg->eval(env)->sign();
        case codes::UnaryOp::floor:
          return arg->eval(env)->floor();
        case codes::UnaryOp::ceil:
          return arg->eval(env)->ceil();
        case codes::UnaryOp::trunc:
          return arg->eval(env)->trunc();
        case codes::UnaryOp::sum:
          return arg->eval(env)->sum();
        case codes::UnaryOp::erf:
          return arg->eval(env)->erf();
        case codes::UnaryOp::erfc:
          return arg->eval(env)->erfc();
        case codes::UnaryOp::gamma:
          return arg->eval(env)->gamma();
        case codes::UnaryOp::lgamma:
          return arg->eval(env)->lgamma();
        case codes::UnaryOp::round:
          return arg->eval(env)->round();
        case codes::UnaryOp::normpdf:
          return arg->eval(env)->normpdf();
        case codes::UnaryOp::normcdf:
          return arg->eval(env)->normcdf();
        case codes::UnaryOp::defined:
          return arg->eval(env)->defined(env);
        }
    }
  catch (StackTrace &ex)
    {
      ex.push("unary operation", location);
      throw;
    }
  catch (exception &e)
    {
      throw StackTrace("unary operation", e.what(), location);
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

BaseTypePtr
BinaryOp::eval(Environment &env) const
{
  try
    {
      switch (op_code)
        {
        case codes::BinaryOp::plus:
          return arg1->eval(env)->plus(arg2->eval(env));
        case codes::BinaryOp::minus:
          return arg1->eval(env)->minus(arg2->eval(env));
        case codes::BinaryOp::times:
          return arg1->eval(env)->times(arg2->eval(env));
        case codes::BinaryOp::divide:
          return arg1->eval(env)->divide(arg2->eval(env));
        case codes::BinaryOp::power:
          return arg1->eval(env)->power(arg2->eval(env));
        case codes::BinaryOp::equal_equal:
          return arg1->eval(env)->is_equal(arg2->eval(env));
        case codes::BinaryOp::not_equal:
          return arg1->eval(env)->is_different(arg2->eval(env));
        case codes::BinaryOp::less:
          return arg1->eval(env)->is_less(arg2->eval(env));
        case codes::BinaryOp::greater:
          return arg1->eval(env)->is_greater(arg2->eval(env));
        case codes::BinaryOp::less_equal:
          return arg1->eval(env)->is_less_equal(arg2->eval(env));
        case codes::BinaryOp::greater_equal:
          return arg1->eval(env)->is_greater_equal(arg2->eval(env));
        case codes::BinaryOp::logical_and:
          return arg1->eval(env)->logical_and(arg2, env);
        case codes::BinaryOp::logical_or:
          return arg1->eval(env)->logical_or(arg2, env);
        case codes::BinaryOp::in:
          return arg2->eval(env)->contains(arg1->eval(env));
        case codes::BinaryOp::set_union:
          return arg1->eval(env)->set_union(arg2->eval(env));
        case codes::BinaryOp::set_intersection:
          return arg1->eval(env)->set_intersection(arg2->eval(env));
        case codes::BinaryOp::max:
          return arg1->eval(env)->max(arg2->eval(env));
        case codes::BinaryOp::min:
          return arg1->eval(env)->min(arg2->eval(env));
        case codes::BinaryOp::mod:
          return arg1->eval(env)->mod(arg2->eval(env));
        }
    }
  catch (StackTrace &ex)
    {
      ex.push("binary operation", location);
      throw;
    }
  catch (exception &e)
    {
      throw StackTrace("binary operation", e.what(), location);
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

BaseTypePtr
TrinaryOp::eval(Environment &env) const
{
  try
    {
      switch (op_code)
        {
        case codes::TrinaryOp::normpdf:
          return arg1->eval(env)->normpdf(arg2->eval(env), arg3->eval(env));
        case codes::TrinaryOp::normcdf:
          return arg1->eval(env)->normcdf(arg2->eval(env), arg3->eval(env));
        }
    }
  catch (StackTrace &ex)
    {
      ex.push("trinary operation", location);
      throw;
    }
  catch (exception &e)
    {
      throw StackTrace("trinary operation", e.what(), location);
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

BaseTypePtr
Comprehension::eval(Environment &env) const
{
  ArrayPtr input_set;
  VariablePtr vp;
  TuplePtr mt;
  try
    {
      input_set = dynamic_pointer_cast<Array>(c_set->eval(env));
      if (!input_set)
        throw StackTrace("Comprehension", "The input set must evaluate to an array", location);
      vp = dynamic_pointer_cast<Variable>(c_vars);
      mt = dynamic_pointer_cast<Tuple>(c_vars);
      if ((!vp && !mt) || (vp && mt))
        throw StackTrace("Comprehension", "the loop variables must be either "
                         "a tuple or a variable", location);
    }
  catch (StackTrace &ex)
    {
      ex.push("Comprehension: ", location);
      throw;
    }

  vector<ExpressionPtr> values;
  for (size_t i = 0; i < input_set->size(); i++)
    {
      auto btp = dynamic_pointer_cast<BaseType>(input_set->at(i));
      if (vp)
        env.define(vp, btp);
      else
        if (btp->getType() == codes::BaseType::Tuple)
          {
            auto mt2 = dynamic_pointer_cast<Tuple>(btp);
            if (mt->size() != mt2->size())
              throw StackTrace("Comprehension", "The number of elements in the input "
                               " set tuple are not the same as the number of elements in "
                               "the output expression tuple", location);

            for (size_t j = 0; j < mt->size(); j++)
              {
                auto vp2 = dynamic_pointer_cast<Variable>(mt->at(j));
                if (!vp2)
                  throw StackTrace("Comprehension", "Output expression tuple must be "
                                   "comprised of variable names", location);
                env.define(vp2, mt2->at(j));
              }
          }
        else
          throw StackTrace("Comprehension", "assigning to tuple in output expression "
                           "but input expression does not contain tuples", location);

      if (!c_when)
        if (!c_expr)
          throw StackTrace("Comprehension", "Internal Error: Impossible case", location);
        else
          values.emplace_back(c_expr->eval(env));
      else
        {
          RealPtr dp;
          BoolPtr bp;
          try
            {
              auto tmp = c_when->eval(env);
              dp = dynamic_pointer_cast<Real>(tmp);
              bp = dynamic_pointer_cast<Bool>(tmp);
              if (!bp && !dp)
                throw StackTrace("The condition must evaluate to a boolean or a real");
            }
          catch (StackTrace &ex)
            {
              ex.push("Comprehension", location);
              throw;
            }
          if ((bp && *bp) || (dp && *dp))
            {
              if (c_expr)
                values.emplace_back(c_expr->eval(env));
              else
                values.emplace_back(btp);
            }
        }
    }
  return make_shared<Array>(values);
}

string
Array::to_string() const noexcept
{
  if (arr.empty())
    return "[]";
  string retval = "[";
  for (const auto &it : arr)
    retval += it->to_string() + ", ";
  return retval.substr(0, retval.size()-2) + "]";
}

string
Tuple::to_string() const noexcept
{
  string retval = "(";
  for (const auto &it : tup)
    retval += it->to_string() + ", ";
  return retval.substr(0, retval.size()-2) + ")";
}

string
Function::to_string() const noexcept
{
  string retval = name + "(";
  for (const auto &it : args)
    retval += it->to_string() + ", ";
  return retval.substr(0, retval.size()-2) + ")";
}

string
UnaryOp::to_string() const noexcept
{
  string retval = arg->to_string();
  switch (op_code)
    {
    case codes::UnaryOp::cast_bool:
      return "(bool)" + retval;
    case codes::UnaryOp::cast_real:
      return "(real)" + retval;
    case codes::UnaryOp::cast_string:
      return "(string)" + retval;
    case codes::UnaryOp::cast_tuple:
      return "(tuple)" + retval;
    case codes::UnaryOp::cast_array:
      return "(array)" + retval;
    case codes::UnaryOp::logical_not:
      return "!" + retval;
    case codes::UnaryOp::unary_minus:
      return "-" + retval;
    case codes::UnaryOp::unary_plus:
      return "+" + retval;
    case codes::UnaryOp::length:
      return "length(" + retval + ")";
    case codes::UnaryOp::isempty:
      return "isempty(" + retval + ")";
    case codes::UnaryOp::isboolean:
      return "isboolean(" + retval + ")";
    case codes::UnaryOp::isreal:
      return "isreal(" + retval + ")";
    case codes::UnaryOp::isstring:
      return "isstring(" + retval + ")";
    case codes::UnaryOp::istuple:
      return "istuple(" + retval + ")";
    case codes::UnaryOp::isarray:
      return "isarray(" + retval + ")";
    case codes::UnaryOp::exp:
      return "exp(" + retval + ")";
    case codes::UnaryOp::ln:
      return "ln(" + retval + ")";
    case codes::UnaryOp::log10:
      return "log10(" + retval + ")";
    case codes::UnaryOp::sin:
      return "sin(" + retval + ")";
    case codes::UnaryOp::cos:
      return "cos(" + retval + ")";
    case codes::UnaryOp::tan:
      return "tan(" + retval + ")";
    case codes::UnaryOp::asin:
      return "asin(" + retval + ")";
    case codes::UnaryOp::acos:
      return "acos(" + retval + ")";
    case codes::UnaryOp::atan:
      return "atan(" + retval + ")";
    case codes::UnaryOp::sqrt:
      return "sqrt(" + retval + ")";
    case codes::UnaryOp::cbrt:
      return "cbrt(" + retval + ")";
    case codes::UnaryOp::sign:
      return "sign(" + retval + ")";
    case codes::UnaryOp::floor:
      return "floor(" + retval + ")";
    case codes::UnaryOp::ceil:
      return "ceil(" + retval + ")";
    case codes::UnaryOp::trunc:
      return "trunc(" + retval + ")";
    case codes::UnaryOp::sum:
      return "sum(" + retval + ")";
    case codes::UnaryOp::erf:
      return "erf(" + retval + ")";
    case codes::UnaryOp::erfc:
      return "erfc(" + retval + ")";
    case codes::UnaryOp::gamma:
      return "gamma(" + retval + ")";
    case codes::UnaryOp::lgamma:
      return "lgamma(" + retval + ")";
    case codes::UnaryOp::round:
      return "round(" + retval + ")";
    case codes::UnaryOp::normpdf:
      return "normpdf(" + retval + ")";
    case codes::UnaryOp::normcdf:
      return "normcdf(" + retval + ")";
    case codes::UnaryOp::defined:
      return "defined(" + retval + ")";
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

string
BinaryOp::to_string() const noexcept
{
  string retval = "(" + arg1->to_string();
  switch (op_code)
    {
    case codes::BinaryOp::plus:
      retval += " + ";
      break;
    case codes::BinaryOp::minus:
      retval += " - ";
      break;
    case codes::BinaryOp::times:
      retval += " * ";
      break;
    case codes::BinaryOp::divide:
      retval += " / ";
      break;
    case codes::BinaryOp::power:
      retval += " ^ ";
      break;
    case codes::BinaryOp::equal_equal:
      retval += " == ";
      break;
    case codes::BinaryOp::not_equal:
      retval += " != ";
      break;
    case codes::BinaryOp::less:
      retval += " < ";
      break;
    case codes::BinaryOp::greater:
      retval += " > ";
      break;
    case codes::BinaryOp::less_equal:
      retval += " <= ";
      break;
    case codes::BinaryOp::greater_equal:
      retval += " >= ";
      break;
    case codes::BinaryOp::logical_and:
      retval += " && ";
      break;
    case codes::BinaryOp::logical_or:
      retval += " || ";
      break;
    case codes::BinaryOp::in:
      retval += " in ";
      break;
    case codes::BinaryOp::set_union:
      return "union(" + retval + ", " + arg2->to_string() + ")";
    case codes::BinaryOp::set_intersection:
      return "intersection(" + retval + ", " + arg2->to_string() + ")";
    case codes::BinaryOp::max:
      return "max(" + retval + ", " + arg2->to_string() + ")";
    case codes::BinaryOp::min:
      return "min(" + retval + ", " + arg2->to_string() + ")";
    case codes::BinaryOp::mod:
      return "mod(" + retval + ", " + arg2->to_string() + ")";
    }
  return retval + arg2->to_string() + ")";
}

string
TrinaryOp::to_string() const noexcept
{
  switch (op_code)
    {
    case codes::TrinaryOp::normpdf:
      return "normpdf(" + arg1->to_string() + ", " + arg2->to_string() + ", " + arg3->to_string() + ")";
    case codes::TrinaryOp::normcdf:
      return "normcdf(" + arg1->to_string() + ", " + arg2->to_string() + ", " + arg3->to_string() + ")";
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

string
Comprehension::to_string() const noexcept
{
  string retval = "[";
  if (c_expr)
    retval += c_expr->to_string() + " for ";
  retval += c_vars->to_string() + " in " + c_set->to_string();
  if (c_when)
    retval += " when " + c_when->to_string();
  return retval + "]";
}

void
String::print(ostream &output, bool matlab_output) const noexcept
{
  output << (matlab_output ? "'" : R"(")")
         << value
         << (matlab_output ? "'" : R"(")");
}

void
Array::print(ostream &output, bool matlab_output) const noexcept
{
  output << (matlab_output ? "{" : "[");
  for (bool printed_something{false};
       auto e : arr)
    {
      if (exchange(printed_something, true))
        output << ", ";
      e->print(output, matlab_output);
    }
  output << (matlab_output ? "}" : "]");
}

void
Tuple::print(ostream &output, bool matlab_output) const noexcept
{
  output << (matlab_output ? "{" : "(");
  for (bool printed_something{false};
       auto e : tup)
    {
      if (exchange(printed_something, true))
        output << ", ";
      e->print(output, matlab_output);
    }
  output << (matlab_output ? "}" : ")");
}

void
Function::printArgs(ostream &output) const noexcept
{
  output << "(";
  for (bool printed_something{false};
       auto e : args)
    {
      if (exchange(printed_something, true))
        output << ", ";
      e->print(output);
    }
  output << ")";
}

void
UnaryOp::print(ostream &output, bool matlab_output) const noexcept
{
  switch (op_code)
    {
    case codes::UnaryOp::cast_bool:
      output << "(bool)";
      break;
    case codes::UnaryOp::cast_real:
      output << "(real)";
      break;
    case codes::UnaryOp::cast_string:
      output << "(string)";
      break;
    case codes::UnaryOp::cast_tuple:
      output << "(tuple)";
      break;
    case codes::UnaryOp::cast_array:
      output << "(array)";
      break;
    case codes::UnaryOp::logical_not:
      output << "!";
      break;
    case codes::UnaryOp::unary_minus:
      output << "-";
      break;
    case codes::UnaryOp::unary_plus:
      output << "+";
      break;
    case codes::UnaryOp::length:
      output << "length(";
      break;
    case codes::UnaryOp::isempty:
      output << "isempty(";
      break;
    case codes::UnaryOp::isboolean:
      output << "isboolean(";
      break;
    case codes::UnaryOp::isreal:
      output << "isreal(";
      break;
    case codes::UnaryOp::isstring:
      output << "isstring(";
      break;
    case codes::UnaryOp::istuple:
      output << "istuple(";
      break;
    case codes::UnaryOp::isarray:
      output << "isarray(";
      break;
    case codes::UnaryOp::exp:
      output << "exp(";
      break;
    case codes::UnaryOp::ln:
      output << "ln(";
      break;
    case codes::UnaryOp::log10:
      output << "log10(";
      break;
    case codes::UnaryOp::sin:
      output << "sin(";
      break;
    case codes::UnaryOp::cos:
      output << "cos(";
      break;
    case codes::UnaryOp::tan:
      output << "tan(";
      break;
    case codes::UnaryOp::asin:
      output << "asin(";
      break;
    case codes::UnaryOp::acos:
      output << "acos(";
      break;
    case codes::UnaryOp::atan:
      output << "atan(";
      break;
    case codes::UnaryOp::sqrt:
      output << "sqrt(";
      break;
    case codes::UnaryOp::cbrt:
      output << "cbrt(";
      break;
    case codes::UnaryOp::sign:
      output << "sign(";
      break;
    case codes::UnaryOp::floor:
      output << "floor(";
      break;
    case codes::UnaryOp::ceil:
      output << "ceil(";
      break;
    case codes::UnaryOp::trunc:
      output << "trunc(";
      break;
    case codes::UnaryOp::sum:
      output << "sum(";
      break;
    case codes::UnaryOp::erf:
      output << "erf(";
      break;
    case codes::UnaryOp::erfc:
      output << "erfc(";
      break;
    case codes::UnaryOp::gamma:
      output << "gamma(";
      break;
    case codes::UnaryOp::lgamma:
      output << "lgamma(";
      break;
    case codes::UnaryOp::round:
      output << "round(";
      break;
    case codes::UnaryOp::normpdf:
      output << "normpdf(";
      break;
    case codes::UnaryOp::normcdf:
      output << "normcdf(";
      break;
    case codes::UnaryOp::defined:
      output << "defined(";
      break;
    }

  arg->print(output, matlab_output);

  if (op_code != codes::UnaryOp::cast_bool
      && op_code != codes::UnaryOp::cast_real
      && op_code != codes::UnaryOp::cast_string
      && op_code != codes::UnaryOp::cast_tuple
      && op_code != codes::UnaryOp::cast_array
      && op_code != codes::UnaryOp::logical_not
      && op_code != codes::UnaryOp::unary_plus
      && op_code != codes::UnaryOp::unary_minus)
    output << ")";
}

void
BinaryOp::print(ostream &output, bool matlab_output) const noexcept
{
  if (op_code == codes::BinaryOp::set_union
      || op_code == codes::BinaryOp::set_intersection
      || op_code == codes::BinaryOp::max
      || op_code == codes::BinaryOp::min
      || op_code == codes::BinaryOp::mod)
    {
      switch (op_code)
        {
        case codes::BinaryOp::set_union:
          output << "union(";
          break;
        case codes::BinaryOp::set_intersection:
          output << "intersection(";
          break;
        case codes::BinaryOp::max:
          output << "max(";
          break;
        case codes::BinaryOp::min:
          output << "min(";
          break;
        case codes::BinaryOp::mod:
          output << "mod(";
          break;
        default:
          break;
        }
      arg1->print(output, matlab_output);
      output << ", ";
      arg2->print(output, matlab_output);
      output << ")";
      return;
    }

  output << "(";
  arg1->print(output, matlab_output);
  switch (op_code)
    {
    case codes::BinaryOp::plus:
      output << " + ";
      break;
    case codes::BinaryOp::minus:
      output << " - ";
      break;
    case codes::BinaryOp::times:
      output << " * ";
      break;
    case codes::BinaryOp::divide:
      output << " / ";
      break;
    case codes::BinaryOp::power:
      output << " ^ ";
      break;
    case codes::BinaryOp::equal_equal:
      output << " == ";
      break;
    case codes::BinaryOp::not_equal:
      output << " != ";
      break;
    case codes::BinaryOp::less:
      output << " < ";
      break;
    case codes::BinaryOp::greater:
      output << " > ";
      break;
    case codes::BinaryOp::less_equal:
      output << " <= ";
      break;
    case codes::BinaryOp::greater_equal:
      output << " >= ";
      break;
    case codes::BinaryOp::logical_and:
      output << " && ";
      break;
    case codes::BinaryOp::logical_or:
      output << " || ";
      break;
    case codes::BinaryOp::in:
      output << " in ";
      break;
    case codes::BinaryOp::set_union:
    case codes::BinaryOp::set_intersection:
    case codes::BinaryOp::max:
    case codes::BinaryOp::min:
    case codes::BinaryOp::mod:
      cerr << "macro::BinaryOp::print: Should not arrive here" << endl;
      exit(EXIT_FAILURE);
    }
  arg2->print(output, matlab_output);
  output << ")";
}

void
TrinaryOp::print(ostream &output, bool matlab_output) const noexcept
{
  switch (op_code)
    {
    case codes::TrinaryOp::normpdf:
      output << "normpdf(";
      break;
    case codes::TrinaryOp::normcdf:
      output << "normcdf(";
      break;
    }
  arg1->print(output, matlab_output);
  output << ", ";
  arg2->print(output, matlab_output);
  output << ", ";
  arg3->print(output, matlab_output);
  output << ")";
}

void
Comprehension::print(ostream &output, bool matlab_output) const noexcept
{
  output << "[";
  if (c_expr)
    {
      c_expr->print(output, matlab_output);
      output << " for ";
    }
  c_vars->print(output, matlab_output);
  output << " in ";
  c_set->print(output, matlab_output);
  if (c_when)
    {
      output << " when ";
      c_when->print(output, matlab_output);
    }
  output << "]";
}
