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

#include <sstream>

#include "MacroValue.hh"

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

StringMV::StringMV(const string &value_arg) : value(value_arg)
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
