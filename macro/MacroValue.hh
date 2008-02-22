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

//! Base class for representing values in macro language
class MacroValue
{
public:
  //! Exception thrown when type error occurs in macro language
  class TypeError
  {
  public:
    const string message;
    TypeError(const string &message_arg) : message(message_arg) {};
  };
  //! Exception thrown when doing an out-of-bounds access through [] operator
  class OutOfBoundsError
  {
  };
  virtual ~MacroValue();
  //! Applies + operator
  /*! Returns a newly allocated value */
  virtual MacroValue *operator+(const MacroValue &mv) const throw (TypeError) = 0;
  //! Applies unary + operator
  /*! Returns a newly allocated value */
  virtual MacroValue *operator+() const throw (TypeError);
  //! Applies - operator
  /*! Returns a newly allocated value */
  virtual MacroValue *operator-(const MacroValue &mv) const throw (TypeError);
  //! Applies unary - operator
  /*! Returns a newly allocated value */
  virtual MacroValue *operator-() const throw (TypeError);
  //! Applies * operator
  /*! Returns a newly allocated value */
  virtual MacroValue *operator*(const MacroValue &mv) const throw (TypeError);
  //! Applies / operator
  /*! Returns a newly allocated value */
  virtual MacroValue *operator/(const MacroValue &mv) const throw (TypeError);
  //! Less comparison
  /*! Returns a newly allocated IntMV, equal to 0 or 1 */
  virtual MacroValue *operator<(const MacroValue &mv) const throw (TypeError);
  //! Greater comparision
  /*! Returns a newly allocated IntMV, equal to 0 or 1 */
  virtual MacroValue *operator>(const MacroValue &mv) const throw (TypeError);
  //! Less or equal comparison
  /*! Returns a newly allocated IntMV, equal to 0 or 1 */
  virtual MacroValue *operator<=(const MacroValue &mv) const throw (TypeError);
  //! Greater or equal comparison
  /*! Returns a newly allocated IntMV, equal to 0 or 1 */
  virtual MacroValue *operator>=(const MacroValue &mv) const throw (TypeError);
  //! Equal comparison
  /*! Returns a newly allocated IntMV, equal to 0 or 1 */
  virtual MacroValue *operator==(const MacroValue &mv) const throw (TypeError) = 0;
  //! Not equal comparison
  /*! Returns a newly allocated IntMV, equal to 0 or 1 */
  virtual MacroValue *operator!=(const MacroValue &mv) const throw (TypeError) = 0;
  //! Applies && operator
  /*! Returns a newly allocated value */
  virtual MacroValue *operator&&(const MacroValue &mv) const throw (TypeError);
  //! Applies || operator
  /*! Returns a newly allocated value */
  virtual MacroValue *operator||(const MacroValue &mv) const throw (TypeError);
  //! Applies unary ! operator
  /*! Returns a newly allocated value */
  virtual MacroValue *operator!() const throw (TypeError);
  //! Applies [] operator
  /*! Returns a newly allocated value */
  virtual MacroValue *operator[](const MacroValue &mv) const throw (TypeError, OutOfBoundsError);
  //! Converts value to string
  virtual string toString() const = 0;
  //! Clones value
  /*! Returns a newly allocated value */
  virtual MacroValue *clone() const = 0;
  //! Converts value to array form
  /*! Returns a newly allocated array value */
  virtual MacroValue *toArray() const = 0;
  //! Appends value at the end of an array
  /*! The first argument must be an array. Returns a newly allocated array. */
  virtual MacroValue *append(const MacroValue &array) const throw (TypeError);
  //! Returns a new IntMV
  /*! Necessary for ArrayMV::operator[] (template issue) */
  static MacroValue *new_base_value(int i);
  //! Returns a new StringMV
  /*! Necessary for ArrayMV::operator[] (template issue) */
  static MacroValue *new_base_value(const string &s);
};

//! Represents an integer value in macro language
class IntMV : public MacroValue
{
  friend class StringMV;
private:
  //! Underlying integer value
  int value;
public:
  IntMV(int value_arg);
  virtual ~IntMV();
  //! Computes arithmetic addition
  /*! Returns a newly allocated value */
  virtual MacroValue *operator+(const MacroValue &mv) const throw (TypeError);
  //! Unary plus
  /*! Returns a clone of itself */
  virtual MacroValue *operator+() const throw (TypeError);
  //! Computes arithmetic substraction
  /*! Returns a newly allocated value */
  virtual MacroValue *operator-(const MacroValue &mv) const throw (TypeError);
  //! Computes opposite
  /*! Returns a newly allocated value */
  virtual MacroValue *operator-() const throw (TypeError);
  //! Computes arithmetic multiplication
  /*! Returns a newly allocated value */
  virtual MacroValue *operator*(const MacroValue &mv) const throw (TypeError);
  //! Computes arithmetic division
  /*! Returns a newly allocated value */
  virtual MacroValue *operator/(const MacroValue &mv) const throw (TypeError);
  virtual MacroValue *operator<(const MacroValue &mv) const throw (TypeError);
  virtual MacroValue *operator>(const MacroValue &mv) const throw (TypeError);
  virtual MacroValue *operator<=(const MacroValue &mv) const throw (TypeError);
  virtual MacroValue *operator>=(const MacroValue &mv) const throw (TypeError);
  virtual MacroValue *operator==(const MacroValue &mv) const throw (TypeError);
  virtual MacroValue *operator!=(const MacroValue &mv) const throw (TypeError);
  //! Computes logical and
  /*! Returns a newly allocated value */
  virtual MacroValue *operator&&(const MacroValue &mv) const throw (TypeError);
  //! Computes logical or
  /*! Returns a newly allocated value */
  virtual MacroValue *operator||(const MacroValue &mv) const throw (TypeError);
  //! Computes logical negation
  /*! Returns a newly allocated value */
  virtual MacroValue *operator!() const throw (TypeError);
  virtual string toString() const;
  virtual MacroValue *clone() const;
  //! Converts value to array form
  /*! Returns a newly allocated integer array containing a single value */
  virtual MacroValue *toArray() const;
  //! Appends value at the end of an array
  /*! The first argument must be an integer array. Returns a newly allocated integer array. */
  virtual MacroValue *append(const MacroValue &array) const throw (TypeError);
  //! Creates a integer range
  /*! Arguments must be of type IntMV.
      Returns a newly allocated integer array containing all integers between mv1 and mv2.
      If mv2 < mv1, constructs the range in decreasing order.
  */
  static MacroValue *new_range(const MacroValue &mv1, const MacroValue &mv2) throw (TypeError);
};

//! Represents a string value in macro language
class StringMV : public MacroValue
{
private:
  //! Underlying string value
  string value;
public:
  StringMV(const string &value_arg);
  virtual ~StringMV();
  //! Computes string concatenation
  /*! Returns a newly allocated value */
  virtual MacroValue *operator+(const MacroValue &mv) const throw (TypeError);
  virtual MacroValue *operator==(const MacroValue &mv) const throw (TypeError);
  virtual MacroValue *operator!=(const MacroValue &mv) const throw (TypeError);
  //! Subscripting operator
  /*! Argument must be an ArrayMV<int>. Indexes begin at 1.
      Returns a newly allocated string. */
  virtual MacroValue *operator[](const MacroValue &mv) const throw (TypeError, OutOfBoundsError);
  //! Returns underlying string value
  virtual string toString() const;
  virtual MacroValue *clone() const;
  //! Converts value to array form
  /*! Returns a newly allocated string array containing a single value */
  virtual MacroValue *toArray() const;
  //! Appends value at the end of an array
  /*! The first argument must be a string array. Returns a newly allocated string array. */
  virtual MacroValue *append(const MacroValue &array) const throw (TypeError);
};

//! Represents an array in macro language
/*! Empty arrays are forbidden */
template<typename T>
class ArrayMV : public MacroValue
{
  friend class IntMV;
  friend class StringMV;
  friend class ArrayMV<string>; // Necessary for operator[] to access values of integer array when subscripting a string array
private:
  //! Underlying vector
  vector<T> values;
public:
  ArrayMV(const vector<T> &values_arg);
  virtual ~ArrayMV();
  //! Computes array concatenation
  /*! Both array must be of same type
      Returns a newly allocated array */
  virtual MacroValue *operator+(const MacroValue &mv) const throw (TypeError);
  virtual MacroValue *operator==(const MacroValue &mv) const throw (TypeError);
  virtual MacroValue *operator!=(const MacroValue &mv) const throw (TypeError);
  //! Subscripting operator
  /*! Argument must be an ArrayMV<int>. Indexes begin at 1.
      If argument is a one-element array, returns a newly-allocated IntMV or String.
      Otherwise returns a newly allocated array. */
  virtual MacroValue *operator[](const MacroValue &mv) const throw (TypeError, OutOfBoundsError);
  //! Returns a string containing the concatenation of string representations of elements
  virtual string toString() const;
  virtual MacroValue *clone() const;
  //! Returns a clone of itself
  virtual MacroValue *toArray() const;
};

template<typename T>
ArrayMV<T>::ArrayMV(const vector<T> &values_arg) : values(values_arg)
{
  if (values.size() == 0)
    throw "Empty arrays forbidden";
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

  vector<T> values_copy(values);
  values_copy.insert(values_copy.end(), mv2->values.begin(), mv2->values.end());
  return new ArrayMV<T>(values_copy);
}

template<typename T>
MacroValue *
ArrayMV<T>::operator==(const MacroValue &mv) const throw (TypeError)
{
  const ArrayMV<T> *mv2 = dynamic_cast<const ArrayMV<T> *>(&mv);
  if (mv2 == NULL)
    return new IntMV(0);
  else
    return new IntMV(values == mv2->values);
}

template<typename T>
MacroValue *
ArrayMV<T>::operator!=(const MacroValue &mv) const throw (TypeError)
{
  const ArrayMV<T> *mv2 = dynamic_cast<const ArrayMV<T> *>(&mv);
  if (mv2 == NULL)
    return new IntMV(1);
  else
    return new IntMV(values != mv2->values);
}

template<typename T>
MacroValue *
ArrayMV<T>::operator[](const MacroValue &mv) const throw (TypeError, OutOfBoundsError)
{
  const ArrayMV<int> *mv2 = dynamic_cast<const ArrayMV<int> *>(&mv);
  if (mv2 == NULL)
    throw TypeError("Expression inside [] must be an integer array");
  vector<T> result;
  for(vector<int>::const_iterator it = mv2->values.begin();
      it != mv2->values.end(); it++)
    {
      if (*it < 1 || *it > (int) values.size())
        throw OutOfBoundsError();
      result.push_back(values[*it - 1]);
    }

  if (result.size() > 1)
    return new ArrayMV<T>(result);
  else
    return MacroValue::new_base_value(result[0]);
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

#endif
