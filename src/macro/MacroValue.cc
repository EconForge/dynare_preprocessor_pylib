/*
 * Copyright (C) 2008-2018 Dynare Team
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

#include <utility>

#include "MacroDriver.hh"

MacroValuePtr
MacroValue::plus(const MacroValuePtr &mv) noexcept(false)
{
  throw TypeError("Operator + does not exist for this type");
}

MacroValuePtr
MacroValue::unary_plus() noexcept(false)
{
  throw TypeError("Unary operator + does not exist for this type");
}

MacroValuePtr
MacroValue::minus(const MacroValuePtr &mv) noexcept(false)
{
  throw TypeError("Operator - does not exist for this type");
}

MacroValuePtr
MacroValue::unary_minus() noexcept(false)
{
  throw TypeError("Unary operator - does not exist for this type");
}

MacroValuePtr
MacroValue::times(const MacroValuePtr &mv) noexcept(false)
{
  throw TypeError("Operator * does not exist for this type");
}

MacroValuePtr
MacroValue::divide(const MacroValuePtr &mv) noexcept(false)
{
  throw TypeError("Operator / does not exist for this type");
}

shared_ptr<IntMV>
MacroValue::is_less(const MacroValuePtr &mv) noexcept(false)
{
  throw TypeError("Operator < does not exist for this type");
}

shared_ptr<IntMV>
MacroValue::is_greater(const MacroValuePtr &mv) noexcept(false)
{
  throw TypeError("Operator > does not exist for this type");
}

shared_ptr<IntMV>
MacroValue::is_less_equal(const MacroValuePtr &mv) noexcept(false)
{
  throw TypeError("Operator <= does not exist for this type");
}

shared_ptr<IntMV>
MacroValue::is_greater_equal(const MacroValuePtr &mv) noexcept(false)
{
  throw TypeError("Operator >= does not exist for this type");
}

shared_ptr<IntMV>
MacroValue::is_different(const MacroValuePtr &mv)
{
  if (is_equal(mv)->value)
    return make_shared<IntMV>(0);
  else
    return make_shared<IntMV>(1);
}

shared_ptr<IntMV>
MacroValue::logical_and(const MacroValuePtr &mv) noexcept(false)
{
  throw TypeError("Operator && does not exist for this type");
}

shared_ptr<IntMV>
MacroValue::logical_or(const MacroValuePtr &mv) noexcept(false)
{
  throw TypeError("Operator || does not exist for this type");
}

shared_ptr<IntMV>
MacroValue::logical_not() noexcept(false)
{
  throw TypeError("Operator ! does not exist for this type");
}

MacroValuePtr
MacroValue::subscript(const MacroValuePtr &mv) noexcept(false)
{
  throw TypeError("Operator [] does not exist for this type");
}

shared_ptr<IntMV>
MacroValue::length() noexcept(false)
{
  throw TypeError("Length not supported for this type");
}

shared_ptr<IntMV>
MacroValue::in(const MacroValuePtr &mv) noexcept(false)
{
  throw TypeError("Second argument of 'in' operator must be an array");
}

IntMV::IntMV(int value_arg) : value{value_arg}
{
}

MacroValuePtr
IntMV::plus(const MacroValuePtr &mv) noexcept(false)
{
  auto mv2 = dynamic_pointer_cast<IntMV>(mv);
  if (!mv2)
    throw TypeError("Type mismatch for operands of + operator");
  return make_shared<IntMV>(value + mv2->value);
}

MacroValuePtr
IntMV::unary_plus() noexcept(false)
{
  return make_shared<IntMV>(value);
}

MacroValuePtr
IntMV::minus(const MacroValuePtr &mv) noexcept(false)
{
  auto mv2 = dynamic_pointer_cast<IntMV>(mv);
  if (!mv2)
    throw TypeError("Type mismatch for operands of - operator");
  return make_shared<IntMV>(value - mv2->value);
}

MacroValuePtr
IntMV::unary_minus() noexcept(false)
{
  return make_shared<IntMV>(-value);
}

MacroValuePtr
IntMV::times(const MacroValuePtr &mv) noexcept(false)
{
  auto mv2 = dynamic_pointer_cast<IntMV>(mv);
  if (!mv2)
    throw TypeError("Type mismatch for operands of * operator");
  return make_shared<IntMV>(value * mv2->value);
}

MacroValuePtr
IntMV::divide(const MacroValuePtr &mv) noexcept(false)
{
  auto mv2 = dynamic_pointer_cast<IntMV>(mv);
  if (!mv2)
    throw TypeError("Type mismatch for operands of / operator");
  if (!mv2->value)
    throw DivisionByZeroError();
  return make_shared<IntMV>(value / mv2->value);
}

shared_ptr<IntMV>
IntMV::is_less(const MacroValuePtr &mv) noexcept(false)
{
  auto mv2 = dynamic_pointer_cast<IntMV>(mv);
  if (!mv2)
    throw TypeError("Type mismatch for operands of < operator");
  return make_shared<IntMV>(value < mv2->value);
}

shared_ptr<IntMV>
IntMV::is_greater(const MacroValuePtr &mv) noexcept(false)
{
  auto mv2 = dynamic_pointer_cast<IntMV>(mv);
  if (!mv2)
    throw TypeError("Type mismatch for operands of > operator");
  return make_shared<IntMV>(value > mv2->value);
}

shared_ptr<IntMV>
IntMV::is_less_equal(const MacroValuePtr &mv) noexcept(false)
{
  auto mv2 = dynamic_pointer_cast<IntMV>(mv);
  if (!mv2)
    throw TypeError("Type mismatch for operands of <= operator");
  return make_shared<IntMV>(value <= mv2->value);
}

shared_ptr<IntMV>
IntMV::is_greater_equal(const MacroValuePtr &mv) noexcept(false)
{
  auto mv2 = dynamic_pointer_cast<IntMV>(mv);
  if (!mv2)
    throw TypeError("Type mismatch for operands of >= operator");
  return make_shared<IntMV>(value >= mv2->value);
}

shared_ptr<IntMV>
IntMV::is_equal(const MacroValuePtr &mv)
{
  auto mv2 = dynamic_pointer_cast<IntMV>(mv);
  if (!mv2)
    return make_shared<IntMV>(0);
  else
    return make_shared<IntMV>(value == mv2->value);
}

shared_ptr<IntMV>
IntMV::logical_and(const MacroValuePtr &mv) noexcept(false)
{
  auto mv2 = dynamic_pointer_cast<IntMV>(mv);
  if (!mv2)
    throw TypeError("Type mismatch for operands of && operator");
  return make_shared<IntMV>(value && mv2->value);
}

shared_ptr<IntMV>
IntMV::logical_or(const MacroValuePtr &mv) noexcept(false)
{
  auto mv2 = dynamic_pointer_cast<IntMV>(mv);
  if (!mv2)
    throw TypeError("Type mismatch for operands of || operator");
  return make_shared<IntMV>(value || mv2->value);
}

shared_ptr<IntMV>
IntMV::logical_not() noexcept(false)
{
  return make_shared<IntMV>(!value);
}

string
IntMV::toString()
{
  return to_string(value);
}

string
IntMV::print()
{
  return toString();
}

StringMV::StringMV(string value_arg) :
  value{move(value_arg)}
{
}

MacroValuePtr
StringMV::plus(const MacroValuePtr &mv) noexcept(false)
{
  auto mv2 = dynamic_pointer_cast<StringMV>(mv);
  if (mv2)
    return make_shared<StringMV>(value + mv2->value);

  throw TypeError("Type mismatch for operands of + operator");
}

shared_ptr<IntMV>
StringMV::is_equal(const MacroValuePtr &mv)
{
  auto mv2 = dynamic_pointer_cast<StringMV>(mv);
  if (mv2)
    return make_shared<IntMV>(value == mv2->value);
  else
    return make_shared<IntMV>(0);
}

MacroValuePtr
StringMV::subscript(const MacroValuePtr &mv) noexcept(false)
{
  string result;

  auto copy_element = [&](int i) {
    if (i < 1 || i > static_cast<int>(value.length()))
      throw OutOfBoundsError();
    result.append(1, value.at(i - 1));
  };

  auto mv2 = dynamic_pointer_cast<IntMV>(mv);
  auto mv3 = dynamic_pointer_cast<ArrayMV>(mv);

  if (mv2)
    copy_element(mv2->value);
  else if (mv3)
    for (auto &v : mv3->values)
      {
        auto v2 = dynamic_pointer_cast<IntMV>(v);
        if (!v2)
          throw TypeError("Expression inside [] must be an integer or an integer array");
        copy_element(v2->value);
      }
  else
    throw TypeError("Expression inside [] must be an integer or an integer array");

  return make_shared<StringMV>(result);
}

string
StringMV::toString()
{
  return value;
}

string
StringMV::print()
{
  return "'" + value + "'";
}

shared_ptr<IntMV>
StringMV::length() noexcept(false)
{
  return make_shared<IntMV>(value.length());
}

FuncMV::FuncMV(vector<string> args_arg, string body_arg) :
  args{move(args_arg)}, body{move(body_arg)}
{
}

shared_ptr<IntMV>
FuncMV::is_equal(const MacroValuePtr &mv)
{
  auto mv2 = dynamic_pointer_cast<FuncMV>(mv);
  if (!mv2 || body != mv2->body)
    return make_shared<IntMV>(0);

  if (args.size() == mv2->args.size())
    for (size_t i = 0; i < args.size(); i++)
      if (args[i] != mv2->args[i])
        return make_shared<IntMV>(0);

  return make_shared<IntMV>(1);
}

string
FuncMV::toString()
{
  return body;
}

string
FuncMV::print()
{
  bool comma_flag = false;
  string retval = "(";
  for (auto it : args)
    {
      if (comma_flag)
          retval += ", ";
      retval += it;
      comma_flag = true;
    }
  retval += ")";
  return retval + " = '" + body + "'";
}

ArrayMV::ArrayMV(vector<MacroValuePtr> values_arg) : values{move(values_arg)}
{
}

MacroValuePtr
ArrayMV::plus(const MacroValuePtr &mv) noexcept(false)
{
  auto mv2 = dynamic_pointer_cast<ArrayMV>(mv);
  if (!mv2)
    throw TypeError("Type mismatch for operands of + operator");

  vector<MacroValuePtr> values_copy{values};
  values_copy.insert(values_copy.end(), mv2->values.begin(), mv2->values.end());
  return make_shared<ArrayMV>(values_copy);
}

MacroValuePtr
ArrayMV::minus(const MacroValuePtr &mv) noexcept(false)
{
  auto mv2 = dynamic_pointer_cast<ArrayMV>(mv);
  if (!mv2)
    throw TypeError("Type mismatch for operands of - operator");

  /* Highly inefficient algorithm for computing set difference
     (but vector<T> is not suited for that...) */
  vector<MacroValuePtr> new_values;
  for (auto &it : values)
    {
      auto it2 = mv2->values.cbegin();
      for (; it2 != mv2->values.cend(); ++it2)
        if (it->is_different(*it2)->value)
          break;
      if (it2 == mv2->values.cend())
        new_values.push_back(it);
    }

  return make_shared<ArrayMV>(new_values);
}

shared_ptr<IntMV>
ArrayMV::is_equal(const MacroValuePtr &mv)
{
  auto mv2 = dynamic_pointer_cast<ArrayMV>(mv);
  if (!mv2 || values.size() != mv2->values.size())
    return make_shared<IntMV>(0);

  auto it = values.cbegin();
  auto it2 = mv2->values.cbegin();
  while (it != values.cend())
    {
      if ((*it)->is_different(*it2)->value)
        return make_shared<IntMV>(0);
      ++it;
      ++it2;
    }
  return make_shared<IntMV>(1);
}

MacroValuePtr
ArrayMV::subscript(const MacroValuePtr &mv) noexcept(false)
{
  vector<MacroValuePtr> result;

  auto copy_element = [&](int i) {
    if (i < 1 || i > static_cast<int>(values.size()))
      throw OutOfBoundsError();
    result.push_back(values[i - 1]);
  };

  auto mv2 = dynamic_pointer_cast<IntMV>(mv);
  auto mv3 = dynamic_pointer_cast<ArrayMV>(mv);

  if (mv2)
    copy_element(mv2->value);
  else if (mv3)
    for (auto &v : mv3->values)
      {
        auto v2 = dynamic_pointer_cast<IntMV>(v);
        if (!v2)
          throw TypeError("Expression inside [] must be an integer or an integer array");
        copy_element(v2->value);
      }
  else
    throw TypeError("Expression inside [] must be an integer or an integer array");

  if (result.size() > 1 || result.size() == 0)
    return make_shared<ArrayMV>(result);
  else
    return result[0];
}

string
ArrayMV::toString()
{
  ostringstream ss;
  for (auto &v : values)
    ss << v->toString();
  return ss.str();
}

shared_ptr<IntMV>
ArrayMV::length() noexcept(false)
{
  return make_shared<IntMV>(values.size());
}

string
ArrayMV::print()
{
  ostringstream ss;
  ss << "[";
  for (auto it = values.begin();
       it != values.end(); it++)
    {
      if (it != values.begin())
        ss << ", ";

      ss << (*it)->print();
    }
  ss << "]";
  return ss.str();
}

shared_ptr<ArrayMV>
ArrayMV::append(MacroValuePtr mv) noexcept(false)
{
  vector<MacroValuePtr> v{values};
  v.push_back(move(mv));
  return make_shared<ArrayMV>(v);
}

shared_ptr<IntMV>
ArrayMV::in(const MacroValuePtr &mv) noexcept(false)
{
  for (auto &v : values)
    if (v->is_equal(mv)->value)
      return make_shared<IntMV>(1);

  return make_shared<IntMV>(0);
}

shared_ptr<ArrayMV>
ArrayMV::range(const MacroValuePtr &mv1, const MacroValuePtr &mv2) noexcept(false)
{
  auto mv1i = dynamic_pointer_cast<IntMV>(mv1);
  auto mv2i = dynamic_pointer_cast<IntMV>(mv2);
  if (!mv1i || !mv2i)
    throw TypeError("Arguments of range operator (:) must be integers");

  vector<MacroValuePtr> result;
  for (int v1 = mv1i->value, v2 = mv2i->value; v1 <= v2; v1++)
    result.push_back(make_shared<IntMV>(v1));
  return make_shared<ArrayMV>(result);
}
