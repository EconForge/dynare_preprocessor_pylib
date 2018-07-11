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

MacroValue::MacroValue(MacroDriver &driver_arg) : driver(driver_arg)
{
  driver.values.insert(this);
}

MacroValue::~MacroValue()
= default;

const MacroValue *
MacroValue::operator+() const noexcept(false)
{
  throw TypeError("Unary operator + does not exist for this type");
}

const MacroValue *
MacroValue::operator-(const MacroValue &mv) const noexcept(false)
{
  throw TypeError("Operator - does not exist for this type");
}

const MacroValue *
MacroValue::operator-() const noexcept(false)
{
  throw TypeError("Unary operator - does not exist for this type");
}

const MacroValue *
MacroValue::operator*(const MacroValue &mv) const noexcept(false)
{
  throw TypeError("Operator * does not exist for this type");
}

const MacroValue *
MacroValue::operator/(const MacroValue &mv) const noexcept(false)
{
  throw TypeError("Operator / does not exist for this type");
}

const MacroValue *
MacroValue::operator<(const MacroValue &mv) const noexcept(false)
{
  throw TypeError("Operator < does not exist for this type");
}

const MacroValue *
MacroValue::operator>(const MacroValue &mv) const noexcept(false)
{
  throw TypeError("Operator > does not exist for this type");
}

const MacroValue *
MacroValue::operator<=(const MacroValue &mv) const noexcept(false)
{
  throw TypeError("Operator <= does not exist for this type");
}

const MacroValue *
MacroValue::operator>=(const MacroValue &mv) const noexcept(false)
{
  throw TypeError("Operator >= does not exist for this type");
}

const MacroValue *
MacroValue::operator&&(const MacroValue &mv) const noexcept(false)
{
  throw TypeError("Operator && does not exist for this type");
}

const MacroValue *
MacroValue::operator||(const MacroValue &mv) const noexcept(false)
{
  throw TypeError("Operator || does not exist for this type");
}

const MacroValue *
MacroValue::operator!() const noexcept(false)
{
  throw TypeError("Operator ! does not exist for this type");
}

const MacroValue *
MacroValue::operator[](const MacroValue &mv) const noexcept(false)
{
  throw TypeError("Operator [] does not exist for this type");
}

const MacroValue *
MacroValue::length() const noexcept(false)
{
  throw TypeError("Length not supported for this type");
}

const MacroValue *
MacroValue::at(int i) const noexcept(false)
{
  throw TypeError("Length not supported for this type");
}

const MacroValue *
MacroValue::append(const MacroValue *mv) const noexcept(false)
{
  throw TypeError("Cannot append an array at the end of another one. Should use concatenation.");
}

const MacroValue *
MacroValue::in(const MacroValue *array) const noexcept(false)
{
  throw TypeError("First argument of 'in' operator cannot be an array");
}

const MacroValue *
MacroValue::new_base_value(MacroDriver &driver, int i)
{
  return new IntMV(driver, i);
}

const MacroValue *
MacroValue::new_base_value(MacroDriver &driver, const string &s)
{
  return new StringMV(driver, s);
}

IntMV::IntMV(MacroDriver &driver, int value_arg) : MacroValue(driver), value(value_arg)
{
}

IntMV::~IntMV()
= default;

const MacroValue *
IntMV::operator+(const MacroValue &mv) const noexcept(false)
{
  const auto *mv2 = dynamic_cast<const IntMV *>(&mv);
  if (mv2 == nullptr)
    throw TypeError("Type mismatch for operands of + operator");
  return new IntMV(driver, value + mv2->value);
}

const MacroValue *
IntMV::operator+() const noexcept(false)
{
  return this;
}

const MacroValue *
IntMV::operator-(const MacroValue &mv) const noexcept(false)
{
  const auto *mv2 = dynamic_cast<const IntMV *>(&mv);
  if (mv2 == nullptr)
    throw TypeError("Type mismatch for operands of - operator");
  return new IntMV(driver, value - mv2->value);
}

const MacroValue *
IntMV::operator-() const noexcept(false)
{
  return new IntMV(driver, -value);
}

const MacroValue *
IntMV::operator*(const MacroValue &mv) const noexcept(false)
{
  const auto *mv2 = dynamic_cast<const IntMV *>(&mv);
  if (mv2 == nullptr)
    throw TypeError("Type mismatch for operands of * operator");
  return new IntMV(driver, value * mv2->value);
}

const MacroValue *
IntMV::operator/(const MacroValue &mv) const noexcept(false)
{
  const auto *mv2 = dynamic_cast<const IntMV *>(&mv);
  if (mv2 == nullptr)
    throw TypeError("Type mismatch for operands of / operator");
  return new IntMV(driver, value / mv2->value);
}

const MacroValue *
IntMV::operator<(const MacroValue &mv) const noexcept(false)
{
  const auto *mv2 = dynamic_cast<const IntMV *>(&mv);
  if (mv2 == nullptr)
    throw TypeError("Type mismatch for operands of < operator");
  return new IntMV(driver, value < mv2->value);
}

const MacroValue *
IntMV::operator>(const MacroValue &mv) const noexcept(false)
{
  const auto *mv2 = dynamic_cast<const IntMV *>(&mv);
  if (mv2 == nullptr)
    throw TypeError("Type mismatch for operands of > operator");
  return new IntMV(driver, value > mv2->value);
}

const MacroValue *
IntMV::operator<=(const MacroValue &mv) const noexcept(false)
{
  const auto *mv2 = dynamic_cast<const IntMV *>(&mv);
  if (mv2 == nullptr)
    throw TypeError("Type mismatch for operands of <= operator");
  return new IntMV(driver, value <= mv2->value);
}

const MacroValue *
IntMV::operator>=(const MacroValue &mv) const noexcept(false)
{
  const auto *mv2 = dynamic_cast<const IntMV *>(&mv);
  if (mv2 == nullptr)
    throw TypeError("Type mismatch for operands of >= operator");
  return new IntMV(driver, value >= mv2->value);
}

const MacroValue *
IntMV::operator==(const MacroValue &mv) const noexcept(false)
{
  const auto *mv2 = dynamic_cast<const IntMV *>(&mv);
  if (mv2 == nullptr)
    return new IntMV(driver, 0);
  else
    return new IntMV(driver, value == mv2->value);
}

const MacroValue *
IntMV::operator!=(const MacroValue &mv) const noexcept(false)
{
  const auto *mv2 = dynamic_cast<const IntMV *>(&mv);
  if (mv2 == nullptr)
    return new IntMV(driver, 1);
  else
    return new IntMV(driver, value != mv2->value);
}

const MacroValue *
IntMV::operator&&(const MacroValue &mv) const noexcept(false)
{
  const auto *mv2 = dynamic_cast<const IntMV *>(&mv);
  if (mv2 == nullptr)
    throw TypeError("Type mismatch for operands of && operator");
  return new IntMV(driver, value && mv2->value);
}

const MacroValue *
IntMV::operator||(const MacroValue &mv) const noexcept(false)
{
  const auto *mv2 = dynamic_cast<const IntMV *>(&mv);
  if (mv2 == nullptr)
    throw TypeError("Type mismatch for operands of || operator");
  return new IntMV(driver, value || mv2->value);
}

const MacroValue *
IntMV::operator!() const noexcept(false)
{
  return new IntMV(driver, !value);
}

string
IntMV::toString() const
{
  return to_string(value);
}

string
IntMV::print() const
{
  return toString();
}

const MacroValue *
IntMV::toArray() const
{
  vector<int> v;
  v.push_back(value);
  return new ArrayMV<int>(driver, v);
}

const MacroValue *
IntMV::append(const MacroValue *array) const noexcept(false)
{
  const auto *array2 = dynamic_cast<const ArrayMV<int> *>(array);
  if (array2 == nullptr)
    throw TypeError("Type mismatch for append operation");

  vector<int> v(array2->values);
  v.push_back(value);
  return new ArrayMV<int>(driver, v);
}

const MacroValue *
IntMV::in(const MacroValue *array) const noexcept(false)
{
  const auto *array2 = dynamic_cast<const ArrayMV<int> *>(array);
  if (array2 == nullptr)
    throw TypeError("Type mismatch for 'in' operator");

  int result = 0;
  for (int v : array2->values)
    if (v == value)
      {
        result = 1;
        break;
      }

  return new IntMV(driver, result);
}

const MacroValue *
IntMV::new_range(MacroDriver &driver, const MacroValue *mv1, const MacroValue *mv2) noexcept(false)
{
  const auto *mv1i = dynamic_cast<const IntMV *>(mv1);
  const auto *mv2i = dynamic_cast<const IntMV *>(mv2);
  if (mv1i == nullptr || mv2i == nullptr)
    throw TypeError("Arguments of range operator (:) must be integers");

  int v1 = mv1i->value;
  int v2 = mv2i->value;

  vector<int> result;
  for (; v1 <= v2; v1++)
    result.push_back(v1);
  return new ArrayMV<int>(driver, result);
}

StringMV::StringMV(MacroDriver &driver, string value_arg) :
  MacroValue(driver), value(move(value_arg))
{
}

StringMV::~StringMV()
= default;

const MacroValue *
StringMV::operator+(const MacroValue &mv) const noexcept(false)
{
  const auto *mv2 = dynamic_cast<const StringMV *>(&mv);
  if (mv2 != nullptr)
    return new StringMV(driver, value + mv2->value);

  const auto *mv3 = dynamic_cast<const FuncMV *>(&mv);
  if (mv3 != nullptr)
    return new StringMV(driver, value + mv3->toString());

  throw TypeError("Type mismatch for operands of + operator");
}

const MacroValue *
StringMV::operator==(const MacroValue &mv) const noexcept(false)
{
  const auto *mv2 = dynamic_cast<const StringMV *>(&mv);
  if (mv2 == nullptr)
    return new IntMV(driver, 0);
  else
    return new IntMV(driver, value == mv2->value);
}

const MacroValue *
StringMV::operator!=(const MacroValue &mv) const noexcept(false)
{
  const auto *mv2 = dynamic_cast<const StringMV *>(&mv);
  if (mv2 == nullptr)
    return new IntMV(driver, 1);
  else
    return new IntMV(driver, value != mv2->value);
}

const MacroValue *
StringMV::operator[](const MacroValue &mv) const noexcept(false)
{
  const auto *mv2 = dynamic_cast<const ArrayMV<int> *>(&mv);
  if (mv2 == nullptr)
    throw TypeError("Expression inside [] must be an integer array");
  string result;
  for (int v : mv2->values)
    {
      if (v < 1 || v > (int) value.length())
        throw OutOfBoundsError();
      char c = value.at(v - 1);
      result.append(1, c);
    }
  return new StringMV(driver, result);
}

string
StringMV::toString() const
{
  return value;
}

string
StringMV::print() const
{
  return "'" + value + "'";
}

const MacroValue *
StringMV::toArray() const
{
  vector<string> v;
  v.push_back(value);
  return new ArrayMV<string>(driver, v);
}

const MacroValue *
StringMV::append(const MacroValue *array) const noexcept(false)
{
  const auto *array2 = dynamic_cast<const ArrayMV<string> *>(array);
  if (array2 == nullptr)
    throw TypeError("Type mismatch for append operation");

  vector<string> v(array2->values);
  v.push_back(value);
  return new ArrayMV<string>(driver, v);
}

const MacroValue *
StringMV::in(const MacroValue *array) const noexcept(false)
{
  const auto *array2 = dynamic_cast<const ArrayMV<string> *>(array);
  if (array2 == nullptr)
    throw TypeError("Type mismatch for 'in' operator");

  int result = 0;
  for (const auto &v : array2->values)
    if (v == value)
      {
        result = 1;
        break;
      }

  return new IntMV(driver, result);
}

template<>
string
ArrayMV<int>::print() const
{
  ostringstream ss;
  ss << "[";
  for (auto it = values.begin();
       it != values.end(); it++)
    {
      if (it != values.begin())
        ss << ", ";

        ss << *it;
    }
  ss << "]";
  return ss.str();
}

template<>
string
ArrayMV<string>::print() const
{
  ostringstream ss;
  ss << "{";
  for (auto it = values.begin();
       it != values.end(); it++)
    {
      if (it != values.begin())
        ss << ", ";

      ss << "'" << *it << "'";
    }
  ss << "}";
  return ss.str();
}

FuncMV::FuncMV(MacroDriver &driver, vector<string> &args_arg, StringMV &value_arg) :
  MacroValue(driver), args(args_arg), value(value_arg)
{
}

FuncMV::~FuncMV()
= default;

const MacroValue *
FuncMV::operator+(const MacroValue &mv) const noexcept(false)
{
  const auto *mv2 = dynamic_cast<const FuncMV *>(&mv);
  if (mv2 != nullptr)
    return value + mv2->value;

  const auto *mv3 = dynamic_cast<const StringMV *>(&mv);
  if (mv3 != nullptr)
    return value + *mv3;

  throw TypeError("Type mismatch for operands of + operator");
}

const MacroValue *
FuncMV::operator==(const MacroValue &mv) const noexcept(false)
{
  const auto *mv2 = dynamic_cast<const FuncMV *>(&mv);
  if (mv2 == nullptr)
    return new IntMV(driver, 0);

  if (value != mv2->value)
    return new IntMV(driver, 0);

  if (args.size() == mv2->args.size())
    for (size_t i = 0; i < args.size(); i++)
      if (args[i] != mv2->args[i])
        return new IntMV(driver, 0);

  return new IntMV(driver, 1);
}

const MacroValue *
FuncMV::operator!=(const MacroValue &mv) const noexcept(false)
{
  if (dynamic_cast<const IntMV *>(*this == mv)->value == 1)
    return new IntMV(driver, 0);
  return new IntMV(driver, 1);
}

string
FuncMV::toString() const
{
  return value.toString();
}

string
FuncMV::print() const
{
  bool comma_flag = false;
  string retval = "(";
  for (const auto it : args)
    {
      if (comma_flag)
          retval += ", ";
      retval += it;
      comma_flag = true;
    }
  retval += ")";
  return retval + " = '" + value.toString() + "'";
}

const MacroValue *
FuncMV::toArray() const
{
  vector<string> v;
  v.push_back(value.toString());
  return new ArrayMV<string>(driver, v);
}
