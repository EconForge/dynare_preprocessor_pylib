/*
 * Copyright © 2019 Dynare Team
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

#include "Expressions.hh"

using namespace macro;

BoolPtr
BaseType::is_different(const BaseTypePtr &btp) const
{
  if (*(this->is_equal(btp)))
    return make_shared<Bool>(false, env);
  return make_shared<Bool>(true, env);
}

BoolPtr
Bool::is_equal(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<Bool>(btp);
  if (!btp2)
    return make_shared<Bool>(false, env);
  return make_shared<Bool>(value == btp2->value, env);
}

BoolPtr
Bool::logical_and(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<Bool>(btp);
  if (btp2)
    return make_shared<Bool>(value && btp2->value, env);

  auto btp3 = dynamic_pointer_cast<Real>(btp);
  if (btp3)
    return make_shared<Bool>(value && *btp3, env);

  throw StackTrace("Type mismatch for operands of && operator");
}

BoolPtr
Bool::logical_or(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<Bool>(btp);
  if (btp2)
    return make_shared<Bool>(value || btp2->value, env);

  auto btp3 = dynamic_pointer_cast<Real>(btp);
  if (btp3)
    return make_shared<Bool>(value || *btp3, env);

  throw StackTrace("Type mismatch for operands of || operator");
}

BoolPtr
Bool::logical_not() const
{
  return make_shared<Bool>(!value, env);
}

BaseTypePtr
Real::plus(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<Real>(btp);
  if (!btp2)
    throw StackTrace("Type mismatch for operands of + operator");
  return make_shared<Real>(value + btp2->value, env);
}

BaseTypePtr
Real::minus(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<Real>(btp);
  if (!btp2)
    throw StackTrace("Type mismatch for operands of - operator");
  return make_shared<Real>(value - btp2->value, env);
}

BaseTypePtr
Real::times(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<Real>(btp);
  if (!btp2)
    throw StackTrace("Type mismatch for operands of * operator");
  return make_shared<Real>(value * btp2->value, env);
}

BaseTypePtr
Real::divide(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<Real>(btp);
  if (!btp2)
    throw StackTrace("Type mismatch for operands of / operator");
  return make_shared<Real>(value / btp2->value, env);
}

BaseTypePtr
Real::power(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<Real>(btp);
  if (!btp2)
    throw StackTrace("Type mismatch for operands of ^ operator");
  return make_shared<Real>(pow(value, btp2->value), env);
}

BoolPtr
Real::is_less(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<Real>(btp);
  if (!btp2)
    throw StackTrace("Type mismatch for operands of < operator");
  return make_shared<Bool>(isless(value, btp2->value), env);
}

BoolPtr
Real::is_greater(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<Real>(btp);
  if (!btp2)
    throw StackTrace("Type mismatch for operands of > operator");
  return make_shared<Bool>(isgreater(value, btp2->value), env);
}

BoolPtr
Real::is_less_equal(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<Real>(btp);
  if (!btp2)
    throw StackTrace("Type mismatch for operands of <= operator");
  return make_shared<Bool>(islessequal(value, btp2->value), env);
}

BoolPtr
Real::is_greater_equal(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<Real>(btp);
  if (!btp2)
    throw StackTrace("Type mismatch for operands of >= operator");
  return make_shared<Bool>(isgreaterequal(value, btp2->value), env);
}

BoolPtr
Real::is_equal(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<Real>(btp);
  if (!btp2)
    return make_shared<Bool>(false, env);
  return make_shared<Bool>(value == btp2->value, env);
}

BoolPtr
Real::logical_and(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<Real>(btp);
  if (btp2)
    return make_shared<Bool>(value && btp2->value, env);

  auto btp3 = dynamic_pointer_cast<Bool>(btp);
  if (btp3)
    return make_shared<Bool>(value && *btp3, env);

  throw StackTrace("Type mismatch for operands of && operator");
}

BoolPtr
Real::logical_or(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<Real>(btp);
  if (!btp2)
    return make_shared<Bool>(value || btp2->value, env);

  auto btp3 = dynamic_pointer_cast<Bool>(btp);
  if (btp3)
    return make_shared<Bool>(value || *btp3, env);

  throw StackTrace("Type mismatch for operands of || operator");
}

BoolPtr
Real::logical_not() const
{
  return make_shared<Bool>(!value, env);
}

RealPtr
Real::max(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<Real>(btp);
  if (!btp2)
    throw StackTrace("Type mismatch for operands of `max` operator");
  return make_shared<Real>(std::max(value, btp2->value), env);
}

RealPtr
Real::min(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<Real>(btp);
  if (!btp2)
    throw StackTrace("Type mismatch for operands of `min` operator");
  return make_shared<Real>(std::min(value, btp2->value), env);
}

RealPtr
Real::mod(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<Real>(btp);
  if (!btp2)
    throw StackTrace("Type mismatch for operands of `mod` operator");
  return make_shared<Real>(std::fmod(value, btp2->value), env);
}

RealPtr
Real::normpdf(const BaseTypePtr &btp1, const BaseTypePtr &btp2) const
{
  auto btp12 = dynamic_pointer_cast<Real>(btp1);
  auto btp22 = dynamic_pointer_cast<Real>(btp2);
  if (!btp12 || !btp22)
    throw StackTrace("Type mismatch for operands of `normpdf` operator");
  return make_shared<Real>((1/(btp22->value*std::sqrt(2*M_PI)*std::exp(pow((value-btp12->value)/btp22->value, 2)/2))), env);
}

RealPtr
Real::normcdf(const BaseTypePtr &btp1, const BaseTypePtr &btp2) const
{
  auto btp12 = dynamic_pointer_cast<Real>(btp1);
  auto btp22 = dynamic_pointer_cast<Real>(btp2);
  if (!btp12 || !btp22)
    throw StackTrace("Type mismatch for operands of `normpdf` operator");
  return make_shared<Real>((0.5*(1+std::erf((value-btp12->value)/btp22->value/M_SQRT2))), env);
}

BaseTypePtr
String::plus(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<String>(btp);
  if (!btp2)
    throw StackTrace("Type mismatch for operands of + operator");
  return make_shared<String>(value + btp2->value, env);
}

BoolPtr
String::is_less(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<String>(btp);
  if (!btp2)
    throw StackTrace("Type mismatch for operands of < operator");
  return make_shared<Bool>(value < btp2->value, env);
}

BoolPtr
String::is_greater(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<String>(btp);
  if (!btp2)
    throw StackTrace("Type mismatch for operands of > operator");
  return make_shared<Bool>(value > btp2->value, env);
}

BoolPtr
String::is_less_equal(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<String>(btp);
  if (!btp2)
    throw StackTrace("Type mismatch for operands of <= operator");
  return make_shared<Bool>(value <= btp2->value, env);
}

BoolPtr
String::is_greater_equal(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<String>(btp);
  if (!btp2)
    throw StackTrace("Type mismatch for operands of >= operator");
  return make_shared<Bool>(value >= btp2->value, env);
}

BoolPtr
String::is_equal(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<String>(btp);
  if (!btp2)
    return make_shared<Bool>(false, env);
  return make_shared<Bool>(value == btp2->value, env);
}

BoolPtr
String::cast_bool() const
{
  string tf = "true";
  if (equal(value.begin(), value.end(), tf.begin(),
            [] (const char& a, const char& b) { return (tolower(a) == tolower(b)); }))
    return make_shared<Bool>(true, env);

  tf = "false";
  if (equal(value.begin(), value.end(), tf.begin(),
            [] (const char& a, const char& b) { return (tolower(a) == tolower(b)); }))
    return make_shared<Bool>(false, env);

  try
    {
      size_t pos = 0;
      double value_d = stod(value, &pos);
      if (pos != value.length())
        throw StackTrace("Entire string not converted");
      return make_shared<Bool>(static_cast<bool>(value_d), env);
    }
  catch (...)
    {
      throw StackTrace(R"(")" + value + R"(" cannot be converted to a boolean)");
    }
}

RealPtr
String::cast_real() const
{
  try
    {
      size_t pos = 0;
      double value_d = stod(value, &pos);
      if (pos != value.length())
        throw StackTrace("Entire string not converted");
      return make_shared<Real>(value_d, env);
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
  return make_shared<Array>(arr_copy, env);
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
  for (auto & it : arr)
    {
      auto itbtp = dynamic_pointer_cast<BaseType>(it);
      auto it2 = btp2->arr.cbegin();
      for (; it2 != btp2->arr.cend(); ++it2)
        {
          auto it2btp = dynamic_pointer_cast<BaseType>(*it2);
          if (*(itbtp->is_equal(it2btp)))
            break;
        }
      if (it2 == btp2->arr.cend())
        arr_copy.emplace_back(itbtp);
    }
  return make_shared<Array>(arr_copy, env);
}

BaseTypePtr
Array::times(const BaseTypePtr &btp) const
{
  vector<ExpressionPtr> values;
  auto btp2 = dynamic_pointer_cast<Array>(btp);
  if (!btp2)
    throw StackTrace("Type mismatch for operands of * operator");

  for (auto & itl : arr)
    for (auto & itr : btp2->getValue())
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
          for (auto & tit : dynamic_pointer_cast<Tuple>(itr)->getValue())
            new_tuple.push_back(tit);
        else
          throw StackTrace("Array::times: unsupported type on rhs");

        values.emplace_back(make_shared<Tuple>(new_tuple, env));
      }

  return make_shared<Array>(values, env);
}

BaseTypePtr
Array::power(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<Real>(btp);
  if (!btp2 || !*(btp2->isinteger()))
    throw StackTrace("The second argument of the power operator (^) must be an integer");

  auto retval = make_shared<Array>(arr, env);
  for (int i = 1; i < *btp2; i++)
    {
      auto btpv = retval->times(make_shared<Array>(arr, env));
      retval = make_shared<Array>(dynamic_pointer_cast<Array>(btpv)->getValue(), env);
    }
  return retval;
}

BoolPtr
Array::is_equal(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<Array>(btp);
  if (!btp2)
    return make_shared<Bool>(false, env);

  if (arr.size() != btp2->arr.size())
    return make_shared<Bool>(false, env);

  for (size_t i = 0; i < arr.size(); i++)
    {
      auto bt = dynamic_pointer_cast<BaseType>(arr[i]);
      auto bt2 = dynamic_pointer_cast<BaseType>(btp2->arr[i]);
      if (!*(bt->is_equal(bt2)))
        return make_shared<Bool>(false, env);
    }
  return make_shared<Bool>(true, env);
}

ArrayPtr
Array::set_union(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<Array>(btp);
  if (!btp2)
    throw StackTrace("Arguments of the union operator (|) must be sets");

  vector<ExpressionPtr> new_values = arr;
  for (auto & it : btp2->arr)
    {
      bool found = false;
      auto it2 = dynamic_pointer_cast<BaseType>(it);
      if (!it2)
        throw StackTrace("Type mismatch for operands of in operator");
      for (auto & nvit : new_values)
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
  return make_shared<Array>(new_values, env);
}

ArrayPtr
Array::set_intersection(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<Array>(btp);
  if (!btp2)
    throw StackTrace("Arguments of the intersection operator (|) must be sets");

  vector<ExpressionPtr> new_values;
  for (auto & it : btp2->arr)
    {
      auto it2 = dynamic_pointer_cast<BaseType>(it);
      if (!it2)
        throw StackTrace("Type mismatch for operands of in operator");
      for (auto & nvit : arr)
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
  return make_shared<Array>(new_values, env);
}

BoolPtr
Array::contains(const BaseTypePtr &btp) const
{
  for (auto & v : arr)
    {
      auto v2 = dynamic_pointer_cast<BaseType>(v);
      if (!v2)
        throw StackTrace("Type mismatch for operands of in operator");
      if (*(v2->is_equal(btp)))
        return make_shared<Bool>(true, env);
    }
  return make_shared<Bool>(false, env);
}

RealPtr
Array::sum() const
{
  double retval = 0;
  for (auto & v : arr)
    {
      auto v2 = dynamic_pointer_cast<Real>(v);
      if (!v2)
        throw StackTrace("Type mismatch for operands of in operator");
      retval += *v2;
    }
  return make_shared<Real>(retval, env);
}

BoolPtr
Array::cast_bool() const
{
  if (arr.size() != 1)
    throw StackTrace("Array must be of size 1 to be cast to a boolean");
  return arr.at(0)->eval()->cast_bool();
}

RealPtr
Array::cast_real() const
{
  if (arr.size() != 1)
    throw StackTrace("Array must be of size 1 to be cast to a real");
  return arr.at(0)->eval()->cast_real();
}

BoolPtr
Tuple::is_equal(const BaseTypePtr &btp) const
{
  auto btp2 = dynamic_pointer_cast<Tuple>(btp);
  if (!btp2)
    return make_shared<Bool>(false, env);

  if (tup.size() != btp2->tup.size())
    return make_shared<Bool>(false, env);

  for (size_t i = 0; i < tup.size(); i++)
    {
      auto bt = dynamic_pointer_cast<BaseType>(tup[i]);
      auto bt2 = dynamic_pointer_cast<BaseType>(btp2->tup[i]);
      if (!*(bt->is_equal(bt2)))
        return make_shared<Bool>(false, env);
    }
  return make_shared<Bool>(true, env);
}

BoolPtr
Tuple::contains(const BaseTypePtr &btp) const
{
  for (auto & v : tup)
    {
      auto v2 = dynamic_pointer_cast<BaseType>(v);
      if (!v2)
        throw StackTrace("Type mismatch for operands of in operator");
      if (*(v2->is_equal(btp)))
        return make_shared<Bool>(true, env);
    }
  return make_shared<Bool>(false, env);
}

BoolPtr
Tuple::cast_bool() const
{
  if (tup.size() != 1)
    throw StackTrace("Tuple must be of size 1 to be cast to a boolean");
  return tup.at(0)->eval()->cast_bool();
}

RealPtr
Tuple::cast_real() const
{
  if (tup.size() != 1)
    throw StackTrace("Tuple must be of size 1 to be cast to a real");
  return tup.at(0)->eval()->cast_real();
}

BaseTypePtr
Range::eval()
{
  RealPtr incdbl = make_shared<Real>(1, env);
  if (inc)
    incdbl = dynamic_pointer_cast<Real>(inc->eval());
  RealPtr startdbl = dynamic_pointer_cast<Real>(start->eval());
  RealPtr enddbl = dynamic_pointer_cast<Real>(end->eval());
  if (!startdbl || !enddbl || !incdbl)
    throw StackTrace("To create an array from a range using the colon operator, "
                     "the arguments must evaluate to reals");

  vector<ExpressionPtr> arr;
  if (*incdbl > 0 && *startdbl <= *enddbl)
    for (double i = *startdbl; i <= *enddbl; i += *incdbl)
      arr.emplace_back(make_shared<Real>(i, env));
  else if (*startdbl >= *enddbl && *incdbl < 0)
    for (double i = *startdbl; i >= *enddbl; i += *incdbl)
      arr.emplace_back(make_shared<Real>(i, env));

  return make_shared<Array>(arr, env, location);
}

BaseTypePtr
Array::eval()
{
  vector<ExpressionPtr> retval;
  for (const auto & it : arr)
    retval.emplace_back(it->eval());
  return make_shared<Array>(retval, env);
}

BaseTypePtr
Tuple::eval()
{
  vector<ExpressionPtr> retval;
  for (const auto & it : tup)
    retval.emplace_back(it->eval());
  return make_shared<Tuple>(retval, env);
}

BaseTypePtr
Variable::eval()
{
  if (indices && !indices->empty())
    {
      ArrayPtr map = dynamic_pointer_cast<Array>(indices->eval());
      vector<ExpressionPtr> index = map->getValue();
      vector<int> ind;
      for (auto it : index)
        {
          // Necessary to handle indexes like: y[1:2,2]
          // In general this evaluates to [[1:2],2] but when subscripting we want to expand it to [1,2,2]
          auto db = dynamic_pointer_cast<Real>(it);
          if (db)
            {
              if (!*(db->isinteger()))
                throw StackTrace("variable", "When indexing a variable you must pass "
                                 "an int or an int array", location);
              ind.emplace_back(*db);
            }
          else if (dynamic_pointer_cast<Array>(it))
            for (auto it1 : dynamic_pointer_cast<Array>(it)->getValue())
              {
                db = dynamic_pointer_cast<Real>(it1);
                if (db)
                  {
                    if (!*(db->isinteger()))
                      throw StackTrace("variable", "When indexing a variable you must pass "
                                       "an int or an int array", location);
                    ind.emplace_back(*db);
                  }
                else
                  throw StackTrace("variable", "You cannot index a variable with a "
                                   "nested array", location);
              }
          else
            throw StackTrace("variable", "You can only index a variable with an int or "
                             "an int array", location);
        }

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
            string orig_string =
              dynamic_pointer_cast<String>(env.getVariable(name))->to_string();
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
            return make_shared<String>(retvals, env);
          }
        case codes::BaseType::Array:
          {
            ArrayPtr ap = dynamic_pointer_cast<Array>(env.getVariable(name));
            vector<BaseTypePtr> retval;
            for (auto it : ind)
              try
                {
                  retval.emplace_back(ap->at(it - 1)->eval());
                }
              catch (const out_of_range &ex)
                {
                  throw StackTrace("variable", "Index out of range", location);
                }

            if (retval.size() == 1)
              return retval.at(0);
            vector<ExpressionPtr> retvala(retval.begin(), retval.end());
            return make_shared<Array>(retvala, env);
          }
        }
    }
  return env.getVariable(name)->eval();
}

BaseTypePtr
Function::eval()
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
    throw StackTrace("Function", "The number of arguments used to call " + name +
                     " does not match the number used in its definition", location);

  try
    {
      for (size_t i = 0; i < func->args.size(); i++)
        {
          VariablePtr mvp = dynamic_pointer_cast<Variable>(func->args.at(i));
          env.define(mvp, args.at(i)->eval());
        }
      auto retval = body->eval();
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
UnaryOp::eval()
{
  try
    {
      auto argbt = arg->eval();
      switch (op_code)
        {
        case codes::UnaryOp::cast_bool:
          return argbt->cast_bool();
        case codes::UnaryOp::cast_real:
          return argbt->cast_real();
        case codes::UnaryOp::cast_string:
          return argbt->cast_string();
        case codes::UnaryOp::cast_tuple:
          return argbt->cast_tuple();
        case codes::UnaryOp::cast_array:
          return argbt->cast_array();
        case codes::UnaryOp::logical_not:
          return argbt->logical_not();
        case codes::UnaryOp::unary_minus:
          return argbt->unary_minus();
        case codes::UnaryOp::unary_plus:
          return argbt->unary_plus();
        case codes::UnaryOp::length:
          return argbt->length();
        case codes::UnaryOp::isempty:
          return argbt->isempty();
        case codes::UnaryOp::isboolean:
          return argbt->isboolean();
        case codes::UnaryOp::isreal:
          return argbt->isreal();
        case codes::UnaryOp::isstring:
          return argbt->isstring();
        case codes::UnaryOp::istuple:
          return argbt->istuple();
        case codes::UnaryOp::isarray:
          return argbt->isarray();
        case codes::UnaryOp::exp:
          return argbt->exp();
        case codes::UnaryOp::ln:
          return argbt->ln();
        case codes::UnaryOp::log10:
          return argbt->log10();
        case codes::UnaryOp::sin:
          return argbt->sin();
        case codes::UnaryOp::cos:
          return argbt->cos();
        case codes::UnaryOp::tan:
          return argbt->tan();
        case codes::UnaryOp::asin:
          return argbt->asin();
        case codes::UnaryOp::acos:
          return argbt->acos();
        case codes::UnaryOp::atan:
          return argbt->atan();
        case codes::UnaryOp::sqrt:
          return argbt->sqrt();
        case codes::UnaryOp::cbrt:
          return argbt->cbrt();
        case codes::UnaryOp::sign:
          return argbt->sign();
        case codes::UnaryOp::floor:
          return argbt->floor();
        case codes::UnaryOp::ceil:
          return argbt->ceil();
        case codes::UnaryOp::trunc:
          return argbt->trunc();
        case codes::UnaryOp::sum:
          return argbt->sum();
        case codes::UnaryOp::erf:
          return argbt->erf();
        case codes::UnaryOp::erfc:
          return argbt->erfc();
        case codes::UnaryOp::gamma:
          return argbt->gamma();
        case codes::UnaryOp::lgamma:
          return argbt->lgamma();
        case codes::UnaryOp::round:
          return argbt->round();
        case codes::UnaryOp::normpdf:
          return argbt->normpdf();
        case codes::UnaryOp::normcdf:
          return argbt->normcdf();
        case codes::UnaryOp::defined:
          return argbt->defined();
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
BinaryOp::eval()
{
  try
    {
      auto arg1bt = arg1->eval();
      auto arg2bt = arg2->eval();
      switch(op_code)
        {
        case codes::BinaryOp::plus:
          return arg1bt->plus(arg2bt);
        case codes::BinaryOp::minus:
          return arg1bt->minus(arg2bt);
        case codes::BinaryOp::times:
          return arg1bt->times(arg2bt);
        case codes::BinaryOp::divide:
          return arg1bt->divide(arg2bt);
        case codes::BinaryOp::power:
          return arg1bt->power( arg2bt);
        case codes::BinaryOp::equal_equal:
          return arg1bt->is_equal(arg2bt);
        case codes::BinaryOp::not_equal:
          return arg1bt->is_different(arg2bt);
        case codes::BinaryOp::less:
          return arg1bt->is_less(arg2bt);
        case codes::BinaryOp::greater:
          return arg1bt->is_greater(arg2bt);
        case codes::BinaryOp::less_equal:
          return arg1bt->is_less_equal(arg2bt);
        case codes::BinaryOp::greater_equal:
          return arg1bt->is_greater_equal(arg2bt);
        case codes::BinaryOp::logical_and:
          return arg1bt->logical_and(arg2bt);
        case codes::BinaryOp::logical_or:
          return arg1bt->logical_or(arg2bt);
        case codes::BinaryOp::in:
          return arg2bt->contains(arg1bt);
        case codes::BinaryOp::set_union:
          return arg1bt->set_union(arg2bt);
        case codes::BinaryOp::set_intersection:
          return arg1bt->set_intersection(arg2bt);
        case codes::BinaryOp::max:
          return arg1bt->max(arg2bt);
        case codes::BinaryOp::min:
          return arg1bt->min(arg2bt);
        case codes::BinaryOp::mod:
          return arg1bt->mod(arg2bt);
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
TrinaryOp::eval()
{
  try
    {
      auto arg1bt = arg1->eval();
      auto arg2bt = arg2->eval();
      auto arg3bt = arg3->eval();
      switch(op_code)
        {
        case codes::TrinaryOp::normpdf:
          return arg1bt->normpdf(arg2bt, arg3bt);
        case codes::TrinaryOp::normcdf:
          return arg1bt->normcdf(arg2bt, arg3bt);
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
Comprehension::eval()
{
  ArrayPtr input_set;
  VariablePtr vp;
  TuplePtr mt;
  try
    {
      input_set = dynamic_pointer_cast<Array>(c_set->eval());
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
          values.emplace_back(c_expr->clone()->eval());
      else
        {
          RealPtr dp;
          BoolPtr bp;
          try
            {
              auto tmp = c_when->eval();
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
            if (c_expr)
              values.emplace_back(c_expr->clone()->eval());
            else
              values.emplace_back(btp);
        }
    }
  return make_shared<Array>(values, env);
}

ExpressionPtr
Tuple::clone() const noexcept
{
  vector<ExpressionPtr> tup_copy;
  for (auto & it : tup)
    tup_copy.emplace_back(it->clone());
  return make_shared<Tuple>(tup_copy, env, location);
}

ExpressionPtr
Array::clone() const noexcept
{
  vector<ExpressionPtr> arr_copy;
  for (auto & it : arr)
    arr_copy.emplace_back(it->clone());
  return make_shared<Array>(arr_copy, env, location);
}

ExpressionPtr
Function::clone() const noexcept
{
  vector<ExpressionPtr> args_copy;
  for (auto & it : args)
    args_copy.emplace_back(it->clone());
  return make_shared<Function>(name, args_copy, env, location);
}

ExpressionPtr
Comprehension::clone() const noexcept
{
  if (c_expr && c_when)
    return make_shared<Comprehension>(c_expr->clone(), c_vars->clone(), c_set->clone(), c_when->clone(), env, location);
  else if (c_expr)
    return make_shared<Comprehension>(c_expr->clone(), c_vars->clone(), c_set->clone(), env, location);
  else
    return make_shared<Comprehension>(true, c_vars->clone(), c_set->clone(), c_when->clone(), env, location);
}

string
Array::to_string() const noexcept
{
  if (arr.empty())
    return "[]";
  string retval = "[";
  for (const auto & it : arr)
    retval += it->to_string() + ", ";
  return retval.substr(0, retval.size()-2) + "]";
}

string
Tuple::to_string() const noexcept
{
  string retval = "(";
  for (const auto & it : tup)
    retval += it->to_string() + ", ";
  return retval.substr(0, retval.size()-2) + ")";
}

string
Function::to_string() const noexcept
{
  string retval = name + "(";
  for (auto & it : args)
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
  switch(op_code)
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
  switch(op_code)
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
  output << (matlab_output ? "'" : R"(")");
  output << value;
  output << (matlab_output ? "'" : R"(")");
}

void
Array::print(ostream &output, bool matlab_output) const noexcept
{
  output << (matlab_output ? "{" : "[");
  for (auto it = arr.begin(); it != arr.end(); it++)
    {
      if (it != arr.begin())
        output << ", ";
      (*it)->print(output, matlab_output);
    }
  output << (matlab_output ? "}" : "]");
}

void
Tuple::print(ostream &output, bool matlab_output) const noexcept
{
  output << (matlab_output ? "{" : "(");
  for (auto it = tup.begin(); it != tup.end(); it++)
    {
      if (it != tup.begin())
        output << ", ";
      (*it)->print(output, matlab_output);
    }
  output << (matlab_output ? "}" : ")");
}

void
Function::printArgs(ostream &output) const noexcept
{
  output << "(";
  for (auto it = args.begin(); it != args.end(); it++)
    {
      if (it != args.begin())
        output << ", ";
      (*it)->print(output);
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
      switch(op_code)
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
  switch(op_code)
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
  switch(op_code)
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

