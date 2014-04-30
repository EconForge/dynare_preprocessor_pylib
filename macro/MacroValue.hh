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

#ifndef _MACRO_VALUE_HH
#define _MACRO_VALUE_HH

#include <string>
#include <vector>
#include <sstream>

using namespace std;

class MacroDriver;

//! Base class for representing values in macro language
class MacroValue
{
protected:
  //! Reference to enclosing MacroDriver
  MacroDriver &driver;
public:
  //! Exception thrown when type error occurs in macro language
  class TypeError
  {
  public:
    const string message;
    TypeError(const string &message_arg) : message(message_arg)
    {
    };
  };
  //! Exception thrown when doing an out-of-bounds access through [] operator
  class OutOfBoundsError
  {
  };
  MacroValue(MacroDriver &driver_arg);
  virtual ~MacroValue();
  //! Applies + operator
  virtual const MacroValue *operator+(const MacroValue &mv) const throw (TypeError) = 0;
  //! Applies unary + operator
  virtual const MacroValue *operator+() const throw (TypeError);
  //! Applies - operator
  virtual const MacroValue *operator-(const MacroValue &mv) const throw (TypeError);
  //! Applies unary - operator
  virtual const MacroValue *operator-() const throw (TypeError);
  //! Applies * operator
  virtual const MacroValue *operator*(const MacroValue &mv) const throw (TypeError);
  //! Applies / operator
  virtual const MacroValue *operator/(const MacroValue &mv) const throw (TypeError);
  //! Less comparison
  /*! Returns an IntMV, equal to 0 or 1 */
  virtual const MacroValue *operator<(const MacroValue &mv) const throw (TypeError);
  //! Greater comparision
  /*! Returns an IntMV, equal to 0 or 1 */
  virtual const MacroValue *operator>(const MacroValue &mv) const throw (TypeError);
  //! Less or equal comparison
  /*! Returns an IntMV, equal to 0 or 1 */
  virtual const MacroValue *operator<=(const MacroValue &mv) const throw (TypeError);
  //! Greater or equal comparison
  /*! Returns an IntMV, equal to 0 or 1 */
  virtual const MacroValue *operator>=(const MacroValue &mv) const throw (TypeError);
  //! Equal comparison
  /*! Returns an IntMV, equal to 0 or 1 */
  virtual const MacroValue *operator==(const MacroValue &mv) const throw (TypeError) = 0;
  //! Not equal comparison
  /*! Returns an IntMV, equal to 0 or 1 */
  virtual const MacroValue *operator!=(const MacroValue &mv) const throw (TypeError) = 0;
  //! Applies && operator
  virtual const MacroValue *operator&&(const MacroValue &mv) const throw (TypeError);
  //! Applies || operator
  virtual const MacroValue *operator||(const MacroValue &mv) const throw (TypeError);
  //! Applies unary ! operator
  virtual const MacroValue *operator!() const throw (TypeError);
  //! Applies [] operator
  virtual const MacroValue *operator[](const MacroValue &mv) const throw (TypeError, OutOfBoundsError);
  //! Converts value to string
  virtual string toString() const = 0;
  //! Converts value to array form
  virtual const MacroValue *toArray() const = 0;
  //! Gets length
  virtual const MacroValue *length() const throw (TypeError);
  //! Appends value at the end of an array
  /*! The argument must be an array. */
  virtual const MacroValue *append(const MacroValue *array) const throw (TypeError);
  //! Applies "in" operator
  /*! The argument must be an array. Returns an IntMV, equal to 0 or 1 */
  virtual const MacroValue *in(const MacroValue *array) const throw (TypeError);
  //! Returns a new IntMV
  /*! Necessary for ArrayMV::operator[] (template issue) */
  static const MacroValue *new_base_value(MacroDriver &driver, int i);
  //! Returns a new StringMV
  /*! Necessary for ArrayMV::operator[] (template issue) */
  static const MacroValue *new_base_value(MacroDriver &driver, const string &s);
};

//! Represents an integer value in macro language
class IntMV : public MacroValue
{
  friend class StringMV;
  friend class MacroDriver;
private:
  //! Underlying integer value
  const int value;
public:
  IntMV(MacroDriver &driver, int value_arg);
  virtual ~IntMV();
  //! Computes arithmetic addition
  virtual const MacroValue *operator+(const MacroValue &mv) const throw (TypeError);
  //! Unary plus
  /*! Returns itself */
  virtual const MacroValue *operator+() const throw (TypeError);
  //! Computes arithmetic substraction
  virtual const MacroValue *operator-(const MacroValue &mv) const throw (TypeError);
  //! Computes opposite
  virtual const MacroValue *operator-() const throw (TypeError);
  //! Computes arithmetic multiplication
  virtual const MacroValue *operator*(const MacroValue &mv) const throw (TypeError);
  //! Computes arithmetic division
  virtual const MacroValue *operator/(const MacroValue &mv) const throw (TypeError);
  virtual const MacroValue *operator<(const MacroValue &mv) const throw (TypeError);
  virtual const MacroValue *operator>(const MacroValue &mv) const throw (TypeError);
  virtual const MacroValue *operator<=(const MacroValue &mv) const throw (TypeError);
  virtual const MacroValue *operator>=(const MacroValue &mv) const throw (TypeError);
  virtual const MacroValue *operator==(const MacroValue &mv) const throw (TypeError);
  virtual const MacroValue *operator!=(const MacroValue &mv) const throw (TypeError);
  //! Computes logical and
  virtual const MacroValue *operator&&(const MacroValue &mv) const throw (TypeError);
  //! Computes logical or
  virtual const MacroValue *operator||(const MacroValue &mv) const throw (TypeError);
  //! Computes logical negation
  virtual const MacroValue *operator!() const throw (TypeError);
  virtual string toString() const;
  //! Converts value to array form
  /*! Returns an integer array containing a single value */
  virtual const MacroValue *toArray() const;
  //! Appends value at the end of an array
  /*! The first argument must be an integer array. */
  virtual const MacroValue *append(const MacroValue *array) const throw (TypeError);
  virtual const MacroValue *in(const MacroValue *array) const throw (TypeError);
  //! Creates a integer range
  /*! Arguments must be of type IntMV.
    Returns an integer array containing all integers between mv1 and mv2.
    If mv2 < mv1, returns an empty range (for consistency with MATLAB).
  */
  static const MacroValue *new_range(MacroDriver &driver, const MacroValue *mv1, const MacroValue *mv2) throw (TypeError);
};

//! Represents a string value in macro language
class StringMV : public MacroValue
{
  friend class MacroDriver;
private:
  //! Underlying string value
  const string value;
public:
  StringMV(MacroDriver &driver, const string &value_arg);
  virtual ~StringMV();
  //! Computes string concatenation
  virtual const MacroValue *operator+(const MacroValue &mv) const throw (TypeError);
  virtual const MacroValue *operator==(const MacroValue &mv) const throw (TypeError);
  virtual const MacroValue *operator!=(const MacroValue &mv) const throw (TypeError);
  //! Subscripting operator
  /*! Argument must be an ArrayMV<int>. Indexes begin at 1. Returns a StringMV. */
  virtual const MacroValue *operator[](const MacroValue &mv) const throw (TypeError, OutOfBoundsError);
  //! Returns underlying string value
  virtual string toString() const;
  //! Converts value to array form
  /*! Returns a string array containing a single value */
  virtual const MacroValue *toArray() const;
  //! Appends value at the end of an array
  /*! The first argument must be a string array. Returns a string array. */
  virtual const MacroValue *append(const MacroValue *array) const throw (TypeError);
  virtual const MacroValue *in(const MacroValue *array) const throw (TypeError);
};

//! Represents an array in macro language
template<typename T>
class ArrayMV : public MacroValue
{
  friend class IntMV;
  friend class StringMV;
  friend class ArrayMV<string>; // Necessary for operator[] to access values of integer array when subscripting a string array
  friend class MacroDriver;
private:
  //! Underlying vector
  const vector<T> values;
public:
  ArrayMV(MacroDriver &driver, const vector<T> &values_arg);
  virtual ~ArrayMV();
  //! Computes array concatenation
  /*! Both array must be of same type */
  virtual const MacroValue *operator+(const MacroValue &mv) const throw (TypeError);
  //! Returns an array in which the elements of the second array have been removed from the first
  /*! It is close to a set difference operation, except that if an element appears two times in the first array, it will also be in the returned value (provided it is not in the second array) */
  virtual const MacroValue *operator-(const MacroValue &mv) const throw (TypeError);
  virtual const MacroValue *operator==(const MacroValue &mv) const throw (TypeError);
  virtual const MacroValue *operator!=(const MacroValue &mv) const throw (TypeError);
  //! Subscripting operator
  /*! Argument must be an ArrayMV<int>. Indexes begin at 1.
    If argument is a one-element array, returns an IntMV or StringMV.
    Otherwise returns an array. */
  virtual const MacroValue *operator[](const MacroValue &mv) const throw (TypeError, OutOfBoundsError);
  //! Returns a string containing the concatenation of string representations of elements
  virtual string toString() const;
  //! Returns itself
  virtual const MacroValue *toArray() const;
  //! Gets length
  virtual const MacroValue *length() const throw (TypeError);
};

template<typename T>
ArrayMV<T>::ArrayMV(MacroDriver &driver, const vector<T> &values_arg) : MacroValue(driver), values(values_arg)
{
}

template<typename T>
ArrayMV<T>::~ArrayMV()
{
}

template<typename T>
const MacroValue *
ArrayMV<T>::operator+(const MacroValue &mv) const throw (TypeError)
{
  const ArrayMV<T> *mv2 = dynamic_cast<const ArrayMV<T> *>(&mv);
  if (mv2 == NULL)
    throw TypeError("Type mismatch for operands of + operator");

  vector<T> values_copy(values);
  values_copy.insert(values_copy.end(), mv2->values.begin(), mv2->values.end());
  return new ArrayMV<T>(driver, values_copy);
}

template<typename T>
const MacroValue *
ArrayMV<T>::operator-(const MacroValue &mv) const throw (TypeError)
{
  const ArrayMV<T> *mv2 = dynamic_cast<const ArrayMV<T> *>(&mv);
  if (mv2 == NULL)
    throw TypeError("Type mismatch for operands of - operator");

  /* Highly inefficient algorithm for computing set difference
     (but vector<T> is not suited for that...) */
  vector<T> new_values;
  for (typename vector<T>::const_iterator it = values.begin();
       it != values.end(); it++)
    {
      typename vector<T>::const_iterator it2;
      for (it2 = mv2->values.begin(); it2 != mv2->values.end(); it2++)
        if (*it == *it2)
          break;
      if (it2 == mv2->values.end())
        new_values.push_back(*it);
    }

  return new ArrayMV<T>(driver, new_values);
}

template<typename T>
const MacroValue *
ArrayMV<T>::operator==(const MacroValue &mv) const throw (TypeError)
{
  const ArrayMV<T> *mv2 = dynamic_cast<const ArrayMV<T> *>(&mv);
  if (mv2 == NULL)
    return new IntMV(driver, 0);
  else
    return new IntMV(driver, values == mv2->values);
}

template<typename T>
const MacroValue *
ArrayMV<T>::operator!=(const MacroValue &mv) const throw (TypeError)
{
  const ArrayMV<T> *mv2 = dynamic_cast<const ArrayMV<T> *>(&mv);
  if (mv2 == NULL)
    return new IntMV(driver, 1);
  else
    return new IntMV(driver, values != mv2->values);
}

template<typename T>
const MacroValue *
ArrayMV<T>::operator[](const MacroValue &mv) const throw (TypeError, OutOfBoundsError)
{
  const ArrayMV<int> *mv2 = dynamic_cast<const ArrayMV<int> *>(&mv);
  if (mv2 == NULL)
    throw TypeError("Expression inside [] must be an integer array");
  vector<T> result;
  for (vector<int>::const_iterator it = mv2->values.begin();
       it != mv2->values.end(); it++)
    {
      if (*it < 1 || *it > (int) values.size())
        throw OutOfBoundsError();
      result.push_back(values[*it - 1]);
    }

  if (result.size() > 1 || result.size() == 0)
    return new ArrayMV<T>(driver, result);
  else
    return MacroValue::new_base_value(driver, result[0]);
}

template<typename T>
string
ArrayMV<T>::toString() const
{
  ostringstream ss;
  for (typename vector<T>::const_iterator it = values.begin();
       it != values.end(); it++)
    ss << *it;
  return ss.str();
}

template<typename T>
const MacroValue *
ArrayMV<T>::toArray() const
{
  return this;
}

template<typename T>
const MacroValue *
ArrayMV<T>::length() const throw (TypeError)
{
  return new IntMV(driver, values.size());
}

#endif
