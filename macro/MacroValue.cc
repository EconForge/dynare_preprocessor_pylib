/*
 * Copyright (C) 2008-2013 Dynare Team
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

#include "MacroDriver.hh"

MacroValue::MacroValue(MacroDriver &driver_arg) : driver(driver_arg)
{
  driver.values.insert(this);
}

MacroValue::~MacroValue()
{
}

const MacroValue *
MacroValue::operator+() const throw (TypeError)
{
  throw TypeError("Unary operator + does not exist for this type");
}

const MacroValue *
MacroValue::operator-(const MacroValue &mv) const throw (TypeError)
{
  throw TypeError("Operator - does not exist for this type");
}

const MacroValue *
MacroValue::operator-() const throw (TypeError)
{
  throw TypeError("Unary operator - does not exist for this type");
}

const MacroValue *
MacroValue::operator*(const MacroValue &mv) const throw (TypeError)
{
  throw TypeError("Operator * does not exist for this type");
}

const MacroValue *
MacroValue::operator/(const MacroValue &mv) const throw (TypeError)
{
  throw TypeError("Operator / does not exist for this type");
}

const MacroValue *
MacroValue::operator<(const MacroValue &mv) const throw (TypeError)
{
  throw TypeError("Operator < does not exist for this type");
}

const MacroValue *
MacroValue::operator>(const MacroValue &mv) const throw (TypeError)
{
  throw TypeError("Operator > does not exist for this type");
}

const MacroValue *
MacroValue::operator<=(const MacroValue &mv) const throw (TypeError)
{
  throw TypeError("Operator <= does not exist for this type");
}

const MacroValue *
MacroValue::operator>=(const MacroValue &mv) const throw (TypeError)
{
  throw TypeError("Operator >= does not exist for this type");
}

const MacroValue *
MacroValue::operator&&(const MacroValue &mv) const throw (TypeError)
{
  throw TypeError("Operator && does not exist for this type");
}

const MacroValue *
MacroValue::operator||(const MacroValue &mv) const throw (TypeError)
{
  throw TypeError("Operator || does not exist for this type");
}

const MacroValue *
MacroValue::operator!() const throw (TypeError)
{
  throw TypeError("Operator ! does not exist for this type");
}

const MacroValue *
MacroValue::operator[](const MacroValue &mv) const throw (TypeError, OutOfBoundsError)
{
  throw TypeError("Operator [] does not exist for this type");
}

const MacroValue *
MacroValue::length() const throw (TypeError)
{
  throw TypeError("Length not supported for this type");
}

const MacroValue *
MacroValue::append(const MacroValue *mv) const throw (TypeError)
{
  throw TypeError("Cannot append an array at the end of another one. Should use concatenation.");
}

const MacroValue *
MacroValue::in(const MacroValue *array) const throw (TypeError)
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
{
}

const MacroValue *
IntMV::operator+(const MacroValue &mv) const throw (TypeError)
{
  const IntMV *mv2 = dynamic_cast<const IntMV *>(&mv);
  if (mv2 == NULL)
    throw TypeError("Type mismatch for operands of + operator");
  return new IntMV(driver, value + mv2->value);
}

const MacroValue *
IntMV::operator+() const throw (TypeError)
{
  return this;
}

const MacroValue *
IntMV::operator-(const MacroValue &mv) const throw (TypeError)
{
  const IntMV *mv2 = dynamic_cast<const IntMV *>(&mv);
  if (mv2 == NULL)
    throw TypeError("Type mismatch for operands of - operator");
  return new IntMV(driver, value - mv2->value);
}

const MacroValue *
IntMV::operator-() const throw (TypeError)
{
  return new IntMV(driver, -value);
}

const MacroValue *
IntMV::operator*(const MacroValue &mv) const throw (TypeError)
{
  const IntMV *mv2 = dynamic_cast<const IntMV *>(&mv);
  if (mv2 == NULL)
    throw TypeError("Type mismatch for operands of * operator");
  return new IntMV(driver, value * mv2->value);
}

const MacroValue *
IntMV::operator/(const MacroValue &mv) const throw (TypeError)
{
  const IntMV *mv2 = dynamic_cast<const IntMV *>(&mv);
  if (mv2 == NULL)
    throw TypeError("Type mismatch for operands of / operator");
  return new IntMV(driver, value / mv2->value);
}

const MacroValue *
IntMV::operator<(const MacroValue &mv) const throw (TypeError)
{
  const IntMV *mv2 = dynamic_cast<const IntMV *>(&mv);
  if (mv2 == NULL)
    throw TypeError("Type mismatch for operands of < operator");
  return new IntMV(driver, value < mv2->value);
}

const MacroValue *
IntMV::operator>(const MacroValue &mv) const throw (TypeError)
{
  const IntMV *mv2 = dynamic_cast<const IntMV *>(&mv);
  if (mv2 == NULL)
    throw TypeError("Type mismatch for operands of > operator");
  return new IntMV(driver, value > mv2->value);
}

const MacroValue *
IntMV::operator<=(const MacroValue &mv) const throw (TypeError)
{
  const IntMV *mv2 = dynamic_cast<const IntMV *>(&mv);
  if (mv2 == NULL)
    throw TypeError("Type mismatch for operands of <= operator");
  return new IntMV(driver, value <= mv2->value);
}

const MacroValue *
IntMV::operator>=(const MacroValue &mv) const throw (TypeError)
{
  const IntMV *mv2 = dynamic_cast<const IntMV *>(&mv);
  if (mv2 == NULL)
    throw TypeError("Type mismatch for operands of >= operator");
  return new IntMV(driver, value >= mv2->value);
}

const MacroValue *
IntMV::operator==(const MacroValue &mv) const throw (TypeError)
{
  const IntMV *mv2 = dynamic_cast<const IntMV *>(&mv);
  if (mv2 == NULL)
    return new IntMV(driver, 0);
  else
    return new IntMV(driver, value == mv2->value);
}

const MacroValue *
IntMV::operator!=(const MacroValue &mv) const throw (TypeError)
{
  const IntMV *mv2 = dynamic_cast<const IntMV *>(&mv);
  if (mv2 == NULL)
    return new IntMV(driver, 1);
  else
    return new IntMV(driver, value != mv2->value);
}

const MacroValue *
IntMV::operator&&(const MacroValue &mv) const throw (TypeError)
{
  const IntMV *mv2 = dynamic_cast<const IntMV *>(&mv);
  if (mv2 == NULL)
    throw TypeError("Type mismatch for operands of && operator");
  return new IntMV(driver, value && mv2->value);
}

const MacroValue *
IntMV::operator||(const MacroValue &mv) const throw (TypeError)
{
  const IntMV *mv2 = dynamic_cast<const IntMV *>(&mv);
  if (mv2 == NULL)
    throw TypeError("Type mismatch for operands of || operator");
  return new IntMV(driver, value || mv2->value);
}

const MacroValue *
IntMV::operator!() const throw (TypeError)
{
  return new IntMV(driver, !value);
}

string
IntMV::toString() const
{
  ostringstream ss;
  ss << value;
  return ss.str();
}

const MacroValue *
IntMV::toArray() const
{
  vector<int> v;
  v.push_back(value);
  return new ArrayMV<int>(driver, v);
}

const MacroValue *
IntMV::append(const MacroValue *array) const throw (TypeError)
{
  const ArrayMV<int> *array2 = dynamic_cast<const ArrayMV<int> *>(array);
  if (array2 == NULL)
    throw TypeError("Type mismatch for append operation");

  vector<int> v(array2->values);
  v.push_back(value);
  return new ArrayMV<int>(driver, v);
}

const MacroValue *
IntMV::in(const MacroValue *array) const throw (TypeError)
{
  const ArrayMV<int> *array2 = dynamic_cast<const ArrayMV<int> *>(array);
  if (array2 == NULL)
    throw TypeError("Type mismatch for 'in' operator");

  int result = 0;
  for (vector<int>::const_iterator it = array2->values.begin();
       it != array2->values.end(); it++)
    if (*it == value)
      {
        result = 1;
        break;
      }

  return new IntMV(driver, result);
}

const MacroValue *
IntMV::new_range(MacroDriver &driver, const MacroValue *mv1, const MacroValue *mv2) throw (TypeError)
{
  const IntMV *mv1i = dynamic_cast<const IntMV *>(mv1);
  const IntMV *mv2i = dynamic_cast<const IntMV *>(mv2);
  if (mv1i == NULL || mv2i == NULL)
    throw TypeError("Arguments of range operator (:) must be integers");

  int v1 = mv1i->value;
  int v2 = mv2i->value;

  vector<int> result;
  for (; v1 <= v2; v1++)
    result.push_back(v1);
  return new ArrayMV<int>(driver, result);
}

StringMV::StringMV(MacroDriver &driver, const string &value_arg) :
  MacroValue(driver), value(value_arg)
{
}

StringMV::~StringMV()
{
}

const MacroValue *
StringMV::operator+(const MacroValue &mv) const throw (TypeError)
{
  const StringMV *mv2 = dynamic_cast<const StringMV *>(&mv);
  if (mv2 == NULL)
    throw TypeError("Type mismatch for operands of + operator");
  return new StringMV(driver, value + mv2->value);
}

const MacroValue *
StringMV::operator==(const MacroValue &mv) const throw (TypeError)
{
  const StringMV *mv2 = dynamic_cast<const StringMV *>(&mv);
  if (mv2 == NULL)
    return new IntMV(driver, 0);
  else
    return new IntMV(driver, value == mv2->value);
}

const MacroValue *
StringMV::operator!=(const MacroValue &mv) const throw (TypeError)
{
  const StringMV *mv2 = dynamic_cast<const StringMV *>(&mv);
  if (mv2 == NULL)
    return new IntMV(driver, 1);
  else
    return new IntMV(driver, value != mv2->value);
}

const MacroValue *
StringMV::operator[](const MacroValue &mv) const throw (TypeError, OutOfBoundsError)
{
  const ArrayMV<int> *mv2 = dynamic_cast<const ArrayMV<int> *>(&mv);
  if (mv2 == NULL)
    throw TypeError("Expression inside [] must be an integer array");
  string result;
  for (vector<int>::const_iterator it = mv2->values.begin();
       it != mv2->values.end(); it++)
    {
      if (*it < 1 || *it > (int) value.length())
        throw OutOfBoundsError();
      char c = value.at(*it - 1);
      result.append(1, c);
    }
  return new StringMV(driver, result);
}

string
StringMV::toString() const
{
  return value;
}

const MacroValue *
StringMV::toArray() const
{
  vector<string> v;
  v.push_back(value);
  return new ArrayMV<string>(driver, v);
}

const MacroValue *
StringMV::append(const MacroValue *array) const throw (TypeError)
{
  const ArrayMV<string> *array2 = dynamic_cast<const ArrayMV<string> *>(array);
  if (array2 == NULL)
    throw TypeError("Type mismatch for append operation");

  vector<string> v(array2->values);
  v.push_back(value);
  return new ArrayMV<string>(driver, v);
}

const MacroValue *
StringMV::in(const MacroValue *array) const throw (TypeError)
{
  const ArrayMV<string> *array2 = dynamic_cast<const ArrayMV<string> *>(array);
  if (array2 == NULL)
    throw TypeError("Type mismatch for 'in' operator");

  int result = 0;
  for (vector<string>::const_iterator it = array2->values.begin();
       it != array2->values.end(); it++)
    if (*it == value)
      {
        result = 1;
        break;
      }

  return new IntMV(driver, result);
}
