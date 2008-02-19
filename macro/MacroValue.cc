/*
 * Copyright (C) 2008 Dynare Team
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

#include "MacroValue.hh"

MacroValue::~MacroValue()
{
}

MacroValue *
MacroValue::operator-(const MacroValue &mv) const throw (TypeError)
{
  throw TypeError("Operator - does not exist for this type");
}

MacroValue *
MacroValue::operator-() const throw (TypeError)
{
  throw TypeError("Operator - does not exist for this type");
}

MacroValue *
MacroValue::operator*(const MacroValue &mv) const throw (TypeError)
{
  throw TypeError("Operator * does not exist for this type");
}

MacroValue *
MacroValue::operator/(const MacroValue &mv) const throw (TypeError)
{
  throw TypeError("Operator / does not exist for this type");
}

MacroValue *
MacroValue::operator<(const MacroValue &mv) const throw (TypeError)
{
  throw TypeError("Operator < does not exist for this type");
}

MacroValue *
MacroValue::operator>(const MacroValue &mv) const throw (TypeError)
{
  throw TypeError("Operator > does not exist for this type");
}

MacroValue *
MacroValue::operator<=(const MacroValue &mv) const throw (TypeError)
{
  throw TypeError("Operator <= does not exist for this type");
}

MacroValue *
MacroValue::operator>=(const MacroValue &mv) const throw (TypeError)
{
  throw TypeError("Operator >= does not exist for this type");
}

MacroValue *
MacroValue::operator&&(const MacroValue &mv) const throw (TypeError)
{
  throw TypeError("Operator && does not exist for this type");
}

MacroValue *
MacroValue::operator||(const MacroValue &mv) const throw (TypeError)
{
  throw TypeError("Operator || does not exist for this type");
}

MacroValue *
MacroValue::operator!() const throw (TypeError)
{
  throw TypeError("Operator ! does not exist for this type");
}

MacroValue *
MacroValue::operator[](const MacroValue &mv) const throw (TypeError)
{
  throw TypeError("Operator [] does not exist for this type");
}

IntMV::IntMV(int value_arg) : value(value_arg)
{
}

IntMV::~IntMV()
{
}

MacroValue *
IntMV::operator+(const MacroValue &mv) const throw (TypeError)
{
  const IntMV *mv2 = dynamic_cast<const IntMV *>(&mv);
  if (mv2 == NULL)
    throw TypeError("Type mismatch for operands of + operator");
  else
    return new IntMV(value + mv2->value);
}

MacroValue *
IntMV::operator-(const MacroValue &mv) const throw (TypeError)
{
  const IntMV *mv2 = dynamic_cast<const IntMV *>(&mv);
  if (mv2 == NULL)
    throw TypeError("Type mismatch for operands of - operator");
  else
    return new IntMV(value - mv2->value);
}

MacroValue *
IntMV::operator-() const throw (TypeError)
{
  return new IntMV(-value);
}

MacroValue *
IntMV::operator*(const MacroValue &mv) const throw (TypeError)
{
  const IntMV *mv2 = dynamic_cast<const IntMV *>(&mv);
  if (mv2 == NULL)
    throw TypeError("Type mismatch for operands of * operator");
  else
    return new IntMV(value * mv2->value);
}

MacroValue *
IntMV::operator/(const MacroValue &mv) const throw (TypeError)
{
  const IntMV *mv2 = dynamic_cast<const IntMV *>(&mv);
  if (mv2 == NULL)
    throw TypeError("Type mismatch for operands of / operator");
  else
    return new IntMV(value / mv2->value);
}

MacroValue *
IntMV::operator<(const MacroValue &mv) const throw (TypeError)
{
  const IntMV *mv2 = dynamic_cast<const IntMV *>(&mv);
  if (mv2 == NULL)
    throw TypeError("Type mismatch for operands of < operator");
  else
    return new IntMV(value < mv2->value);
}

MacroValue *
IntMV::operator>(const MacroValue &mv) const throw (TypeError)
{
  const IntMV *mv2 = dynamic_cast<const IntMV *>(&mv);
  if (mv2 == NULL)
    throw TypeError("Type mismatch for operands of > operator");
  else
    return new IntMV(value > mv2->value);
}

MacroValue *
IntMV::operator<=(const MacroValue &mv) const throw (TypeError)
{
  const IntMV *mv2 = dynamic_cast<const IntMV *>(&mv);
  if (mv2 == NULL)
    throw TypeError("Type mismatch for operands of <= operator");
  else
    return new IntMV(value <= mv2->value);
}

MacroValue *
IntMV::operator>=(const MacroValue &mv) const throw (TypeError)
{
  const IntMV *mv2 = dynamic_cast<const IntMV *>(&mv);
  if (mv2 == NULL)
    throw TypeError("Type mismatch for operands of >= operator");
  else
    return new IntMV(value >= mv2->value);
}

MacroValue *
IntMV::operator==(const MacroValue &mv) const throw (TypeError)
{
  const IntMV *mv2 = dynamic_cast<const IntMV *>(&mv);
  if (mv2 == NULL)
    throw new IntMV(0);
  else
    return new IntMV(value == mv2->value);
}

MacroValue *
IntMV::operator!=(const MacroValue &mv) const throw (TypeError)
{
  const IntMV *mv2 = dynamic_cast<const IntMV *>(&mv);
  if (mv2 == NULL)
    throw new IntMV(1);
  else
    return new IntMV(value != mv2->value);
}

MacroValue *
IntMV::operator&&(const MacroValue &mv) const throw (TypeError)
{
  const IntMV *mv2 = dynamic_cast<const IntMV *>(&mv);
  if (mv2 == NULL)
    throw TypeError("Type mismatch for operands of && operator");
  else
    return new IntMV(value && mv2->value);
}

MacroValue *
IntMV::operator||(const MacroValue &mv) const throw (TypeError)
{
  const IntMV *mv2 = dynamic_cast<const IntMV *>(&mv);
  if (mv2 == NULL)
    throw TypeError("Type mismatch for operands of || operator");
  else
    return new IntMV(value || mv2->value);
}

MacroValue *
IntMV::operator!() const throw (TypeError)
{
  return new IntMV(!value);
}

string
IntMV::toString() const
{
  ostringstream ss;
  ss << value;
  return ss.str();
}

MacroValue *
IntMV::clone() const
{
  return new IntMV(value);
}

MacroValue *
IntMV::toArray() const
{
  vector<int> v;
  v.push_back(value);
  return new ArrayMV<int>(v);
}

MacroValue *
IntMV::append(const MacroValue &array) const
{
  const ArrayMV<int> *array2 = dynamic_cast<const ArrayMV<int> *>(&array);
  if (array2 == NULL)
    throw TypeError("Type mismatch for append operation");
  else
    {
      vector<int> v(array2->values);
      v.push_back(value);
      return new ArrayMV<int>(v);
    }
}

StringMV::StringMV(const string &value_arg) : value(value_arg)
{
}

StringMV::~StringMV()
{
}

MacroValue *
StringMV::operator+(const MacroValue &mv) const throw (TypeError)
{
  const StringMV *mv2 = dynamic_cast<const StringMV *>(&mv);
  if (mv2 == NULL)
    throw TypeError("Type mismatch for operands of + operator");
  else
    return new StringMV(value + mv2->value);
}

MacroValue *
StringMV::operator==(const MacroValue &mv) const throw (TypeError)
{
  const StringMV *mv2 = dynamic_cast<const StringMV *>(&mv);
  if (mv2 == NULL)
    throw new IntMV(0);
  else
    return new IntMV(value == mv2->value);
}

MacroValue *
StringMV::operator!=(const MacroValue &mv) const throw (TypeError)
{
  const StringMV *mv2 = dynamic_cast<const StringMV *>(&mv);
  if (mv2 == NULL)
    throw new IntMV(1);
  else
    return new IntMV(value != mv2->value);
}

MacroValue *
StringMV::operator[](const MacroValue &mv) const throw (TypeError)
{
  const ArrayMV<int> *mv2 = dynamic_cast<const ArrayMV<int> *>(&mv);
  if (mv2 == NULL)
    throw TypeError("Expression inside [] must be an integer array");
  string result;
  for(vector<int>::const_iterator it = mv2->values.begin();
      it != mv2->values.end(); it++)
    {
      char c = value.at(*it - 1);
      result.append(&c);
    }
  return new StringMV(result);
}

string
StringMV::toString() const
{
  return value;
}

MacroValue *
StringMV::clone() const
{
  return new StringMV(value);
}

MacroValue *
StringMV::toArray() const
{
  vector<string> v;
  v.push_back(value);
  return new ArrayMV<string>(v);
}

MacroValue *
StringMV::append(const MacroValue &array) const
{
  const ArrayMV<string> *array2 = dynamic_cast<const ArrayMV<string> *>(&array);
  if (array2 == NULL)
    throw TypeError("Type mismatch for append operation");
  else
    {
      vector<string> v(array2->values);
      v.push_back(value);
      return new ArrayMV<string>(v);
    }
}
