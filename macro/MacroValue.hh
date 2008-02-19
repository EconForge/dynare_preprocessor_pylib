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

#ifndef _MACRO_VALUE_HH
#define _MACRO_VALUE_HH

using namespace std;

#include <string>
#include <vector>
#include <sstream>

class MacroValue
{
public:
  class TypeError
  {
  public:
    const string message;
    TypeError(const string &message_arg) : message(message_arg) {};
  };
  virtual ~MacroValue();
  virtual MacroValue *operator+(const MacroValue &mv) const throw (TypeError) = 0;
  virtual MacroValue *operator-(const MacroValue &mv) const throw (TypeError);
  virtual MacroValue *operator-() const throw (TypeError);
  virtual MacroValue *operator*(const MacroValue &mv) const throw (TypeError);
  virtual MacroValue *operator/(const MacroValue &mv) const throw (TypeError);
  virtual MacroValue *operator<(const MacroValue &mv) const throw (TypeError);
  virtual MacroValue *operator>(const MacroValue &mv) const throw (TypeError);
  virtual MacroValue *operator<=(const MacroValue &mv) const throw (TypeError);
  virtual MacroValue *operator>=(const MacroValue &mv) const throw (TypeError);
  virtual MacroValue *operator==(const MacroValue &mv) const throw (TypeError) = 0;
  virtual MacroValue *operator!=(const MacroValue &mv) const throw (TypeError) = 0;
  virtual MacroValue *operator&&(const MacroValue &mv) const throw (TypeError);
  virtual MacroValue *operator||(const MacroValue &mv) const throw (TypeError);
  virtual MacroValue *operator!() const throw (TypeError);
  virtual MacroValue *operator[](const MacroValue &mv) const throw (TypeError);
  virtual string toString() const = 0;
  virtual MacroValue *clone() const = 0;
  virtual MacroValue *toArray() const = 0;
  virtual MacroValue *append(const MacroValue &array) const = 0;
};

class IntMV : public MacroValue
{
  friend class StringMV;
private:
  int value;
public:
  IntMV(int value_arg);
  virtual ~IntMV();
  virtual MacroValue *operator+(const MacroValue &mv) const throw (TypeError);
  virtual MacroValue *operator-(const MacroValue &mv) const throw (TypeError);
  virtual MacroValue *operator-() const throw (TypeError);
  virtual MacroValue *operator*(const MacroValue &mv) const throw (TypeError);
  virtual MacroValue *operator/(const MacroValue &mv) const throw (TypeError);
  virtual MacroValue *operator<(const MacroValue &mv) const throw (TypeError);
  virtual MacroValue *operator>(const MacroValue &mv) const throw (TypeError);
  virtual MacroValue *operator<=(const MacroValue &mv) const throw (TypeError);
  virtual MacroValue *operator>=(const MacroValue &mv) const throw (TypeError);
  virtual MacroValue *operator==(const MacroValue &mv) const throw (TypeError);
  virtual MacroValue *operator!=(const MacroValue &mv) const throw (TypeError);
  virtual MacroValue *operator&&(const MacroValue &mv) const throw (TypeError);
  virtual MacroValue *operator||(const MacroValue &mv) const throw (TypeError);
  virtual MacroValue *operator!() const throw (TypeError);
  virtual string toString() const;
  virtual MacroValue *clone() const;
  virtual MacroValue *toArray() const;
  virtual MacroValue *append(const MacroValue &array) const;
};

class StringMV : public MacroValue
{
private:
  string value;
public:
  StringMV(const string &value_arg);
  virtual ~StringMV();
  virtual MacroValue *operator+(const MacroValue &mv) const throw (TypeError);
//   virtual MacroValue *operator<(const MacroValue &mv) const throw (TypeError);
//   virtual MacroValue *operator>(const MacroValue &mv) const throw (TypeError);
//   virtual MacroValue *operator<=(const MacroValue &mv) const throw (TypeError);
//   virtual MacroValue *operator>=(const MacroValue &mv) const throw (TypeError);
  virtual MacroValue *operator==(const MacroValue &mv) const throw (TypeError);
  virtual MacroValue *operator!=(const MacroValue &mv) const throw (TypeError);
  //! Subscripting operator
  /*!
    Argument must be an ArrayMV<int>. Indexes begin at 1.
    \todo Add bound error checking
  */
  virtual MacroValue *operator[](const MacroValue &mv) const throw (TypeError);
  virtual string toString() const;
  virtual MacroValue *clone() const;
  virtual MacroValue *toArray() const;
  virtual MacroValue *append(const MacroValue &array) const;
};

template<typename T>
class ArrayMV : public MacroValue
{
  friend class IntMV;
  friend class StringMV;
  friend class ArrayMV<string>; // Necessary for operator[] to access values of integer array when subscripting a string array
private:
  vector<T> values;
public:
  ArrayMV(const vector<T> &values_arg);
  virtual ~ArrayMV();
  virtual MacroValue *operator+(const MacroValue &mv) const throw (TypeError);
 virtual MacroValue *operator==(const MacroValue &mv) const throw (TypeError);
  virtual MacroValue *operator!=(const MacroValue &mv) const throw (TypeError);
  //! Subscripting operator
  /*!
    Argument must be an ArrayMV<int>. Indexes begin at 1.
    \todo Add bound error checking
  */
  virtual MacroValue *operator[](const MacroValue &mv) const throw (TypeError);
  virtual string toString() const;
  virtual MacroValue *clone() const;
  virtual MacroValue *toArray() const;
  virtual MacroValue *append(const MacroValue &array) const;
};

template<typename T>
ArrayMV<T>::ArrayMV(const vector<T> &values_arg) : values(values_arg)
{
}

template<typename T>
ArrayMV<T>::~ArrayMV()
{
}

template<typename T>
MacroValue *
ArrayMV<T>::operator+(const MacroValue &mv) const throw (TypeError)
{
  const ArrayMV<T> *mv2 = dynamic_cast<const ArrayMV<T> *>(&mv);
  if (mv2 == NULL)
    throw TypeError("Type mismatch for operands of + operator");
  else
    {
      vector<T> values_copy(values);
      values_copy.insert(values_copy.end(),
                         mv2->values.begin(),
                         mv2->values.end());
      return new ArrayMV<T>(values_copy);
    }
}

template<typename T>
MacroValue *
ArrayMV<T>::operator==(const MacroValue &mv) const throw (TypeError)
{
  const ArrayMV<T> *mv2 = dynamic_cast<const ArrayMV<T> *>(&mv);
  if (mv2 == NULL)
    throw new IntMV(0);
  else
    return new IntMV(values == mv2->values);
}

template<typename T>
MacroValue *
ArrayMV<T>::operator!=(const MacroValue &mv) const throw (TypeError)
{
  const ArrayMV<T> *mv2 = dynamic_cast<const ArrayMV<T> *>(&mv);
  if (mv2 == NULL)
    throw new IntMV(1);
  else
    return new IntMV(values != mv2->values);
}

template<typename T>
MacroValue *
ArrayMV<T>::operator[](const MacroValue &mv) const throw (TypeError)
{
  const ArrayMV<int> *mv2 = dynamic_cast<const ArrayMV<int> *>(&mv);
  if (mv2 == NULL)
    throw TypeError("Expression inside [] must be an integer array");
  vector<T> result;
  for(vector<int>::const_iterator it = mv2->values.begin();
      it != mv2->values.end(); it++)
    result.push_back(values[*it - 1]);
  return new ArrayMV<T>(result);
}

template<typename T>
string
ArrayMV<T>::toString() const
{
  ostringstream ss;
  for(typename vector<T>::const_iterator it = values.begin();
      it != values.end(); it++)
    ss << *it;
  return ss.str();
}

template<typename T>
MacroValue *
ArrayMV<T>::clone() const
{
  return new ArrayMV<T>(values);
}

template<typename T>
MacroValue *
ArrayMV<T>::toArray() const
{
  return clone();
}

template<typename T>
MacroValue *
ArrayMV<T>::append(const MacroValue &mv) const
{
  throw TypeError("Cannot append an array at the end of another one. Should use concatenation.");
}

#endif
