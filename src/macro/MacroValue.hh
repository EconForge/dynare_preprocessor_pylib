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

#ifndef _MACRO_VALUE_HH
#define _MACRO_VALUE_HH

#include <string>
#include <utility>
#include <vector>
#include <sstream>
#include <memory>

using namespace std;

class MacroDriver;
class MacroValue;
class IntMV;
class ArrayMV;

using MacroValuePtr = shared_ptr<MacroValue>;

//! Base class for representing values in macro language
/*! All its derived types are immutable, so that we don't need to carry const
    qualifiers everywhere in the code */
class MacroValue
{
public:
  //! Exception thrown when type error occurs in macro language
  class TypeError
  {
  public:
    const string message;
    TypeError(string message_arg) : message{move(message_arg)}
    {
    };
  };
  //! Exception thrown when doing an out-of-bounds access through [] operator
  class OutOfBoundsError
  {
  };
  //! Exception thrown when dividing by zero
  class DivisionByZeroError
  {
  };
  //! Applies + operator
  virtual MacroValuePtr plus(const MacroValuePtr &mv) noexcept(false);
  //! Applies unary + operator
  virtual MacroValuePtr unary_plus() noexcept(false);
  //! Applies - operator
  virtual MacroValuePtr minus(const MacroValuePtr &mv) noexcept(false);
  //! Applies unary - operator
  virtual MacroValuePtr unary_minus() noexcept(false);
  //! Applies * operator
  virtual MacroValuePtr times(const MacroValuePtr &mv) noexcept(false);
  //! Applies / operator
  virtual MacroValuePtr divide(const MacroValuePtr &mv) noexcept(false);
  //! Less comparison
  /*! Returns an IntMV, equal to 0 or 1 */
  virtual shared_ptr<IntMV> is_less(const MacroValuePtr &mv) noexcept(false);
  //! Greater comparision
  /*! Returns an IntMV, equal to 0 or 1 */
  virtual shared_ptr<IntMV> is_greater(const MacroValuePtr &mv) noexcept(false);
  //! Less or equal comparison
  /*! Returns an IntMV, equal to 0 or 1 */
  virtual shared_ptr<IntMV> is_less_equal(const MacroValuePtr &mv) noexcept(false);
  //! Greater or equal comparison
  /*! Returns an IntMV, equal to 0 or 1 */
  virtual shared_ptr<IntMV> is_greater_equal(const MacroValuePtr &mv) noexcept(false);
  //! Equal comparison
  /*! Returns an IntMV, equal to 0 or 1 */
  virtual shared_ptr<IntMV> is_equal(const MacroValuePtr &mv) = 0;
  //! Not equal comparison
  /*! Returns an IntMV, equal to 0 or 1 */
  shared_ptr<IntMV> is_different(const MacroValuePtr &mv);
  //! Applies && operator
  virtual shared_ptr<IntMV> logical_and(const MacroValuePtr &mv) noexcept(false);
  //! Applies || operator
  virtual shared_ptr<IntMV> logical_or(const MacroValuePtr &mv) noexcept(false);
  //! Applies unary ! operator
  virtual shared_ptr<IntMV> logical_not() noexcept(false);
  //! Applies [] operator
  virtual MacroValuePtr subscript(const MacroValuePtr &mv) noexcept(false);
  //! Converts value to string
  virtual string toString() = 0;
  //! Converts value to be printed
  virtual string print() = 0;
  //! Gets length
  virtual shared_ptr<IntMV> length() noexcept(false);
  //! Applies "in" operator
  /*! The argument is the element to be tested for inclusion. Returns an IntMV, equal to 0 or 1 */
  virtual shared_ptr<IntMV> in(const MacroValuePtr &mv) noexcept(false);
  //! Creates the union of two sets
  virtual shared_ptr<ArrayMV> set_union(const MacroValuePtr &mv) noexcept(false);
  //! Creates the intersection of two sets
  virtual shared_ptr<ArrayMV> set_intersection(const MacroValuePtr &mv) noexcept(false);
  //! Power as shortcut for Cartesian product
  virtual MacroValuePtr power(const MacroValuePtr &mv) noexcept(false);
};

//! Represents an integer value in macro language
class IntMV : public MacroValue
{
public:
  IntMV(int value_arg);
  
  //! Underlying integer value
  const int value;

  //! Computes arithmetic addition
  MacroValuePtr plus(const MacroValuePtr &mv) noexcept(false) override;
  //! Unary plus
  /*! Returns itself */
  MacroValuePtr unary_plus() noexcept(false) override;
  //! Computes arithmetic substraction
  MacroValuePtr minus(const MacroValuePtr &mv) noexcept(false) override;
  //! Computes opposite
  MacroValuePtr unary_minus() noexcept(false) override;
  //! Computes arithmetic multiplication
  MacroValuePtr times(const MacroValuePtr &mv) noexcept(false) override;
  //! Computes arithmetic division
  MacroValuePtr divide(const MacroValuePtr &mv) noexcept(false) override;
  MacroValuePtr power(const MacroValuePtr &mv) noexcept(false) override;
  shared_ptr<IntMV> is_less(const MacroValuePtr &mv) noexcept(false) override;
  shared_ptr<IntMV> is_greater(const MacroValuePtr &mv) noexcept(false) override;
  shared_ptr<IntMV> is_less_equal(const MacroValuePtr &mv) noexcept(false) override;
  shared_ptr<IntMV> is_greater_equal(const MacroValuePtr &mv) noexcept(false) override;
  shared_ptr<IntMV> is_equal(const MacroValuePtr &mv) override;
  //! Computes logical and
  shared_ptr<IntMV> logical_and(const MacroValuePtr &mv) noexcept(false) override;
  //! Computes logical or
  shared_ptr<IntMV> logical_or(const MacroValuePtr &mv) noexcept(false) override;
  //! Computes logical negation
  shared_ptr<IntMV> logical_not() noexcept(false) override;
  string toString() override;
  string print() override;
};

//! Represents a string value in macro language
class StringMV : public MacroValue
{
public:
  StringMV(string value_arg);
  //! Underlying string value
  const string value;
  
  //! Computes string concatenation
  MacroValuePtr plus(const MacroValuePtr &mv) noexcept(false) override;
  shared_ptr<IntMV> is_equal(const MacroValuePtr &mv) override;
  //! Subscripting operator
  /*! Argument must be an ArrayMV<int>. Indexes begin at 1. Returns a StringMV. */
  MacroValuePtr subscript(const MacroValuePtr &mv) noexcept(false) override;
  //! Returns underlying string value
  string toString() override;
  string print() override;
  shared_ptr<IntMV> length() noexcept(false) override;
};

class FuncMV : public MacroValue
{
public:
  FuncMV(vector<string> args, string body_arg);
  //! Function args & body
  const vector<string> args;
  const string body;

  shared_ptr<IntMV> is_equal(const MacroValuePtr &mv) override;
  string toString() override;
  string print() override;
};

//! Represents an array in macro language
class ArrayMV : public MacroValue
{
public:
  ArrayMV(vector<MacroValuePtr> values_arg);

  //! Underlying vector
  const vector<MacroValuePtr> values;

  //! Computes array concatenation
  /*! Both array must be of same type */
  MacroValuePtr plus(const MacroValuePtr &mv) noexcept(false) override;
  //! Returns an array in which the elements of the second array have been removed from the first
  /*! It is close to a set difference operation, except that if an element appears two times in the first array, it will also be in the returned value (provided it is not in the second array) */
  MacroValuePtr minus(const MacroValuePtr &mv) noexcept(false) override;
  shared_ptr<IntMV> is_equal(const MacroValuePtr &mv) override;
  //! Subscripting operator
  /*! Argument must be an ArrayMV<int>. Indexes begin at 1.
    If argument is a one-element array, returns an IntMV or StringMV.
    Otherwise returns an array. */
  MacroValuePtr subscript(const MacroValuePtr &mv) noexcept(false) override;
  //! Returns a string containing the concatenation of string representations of elements
  string toString() override;
  string print() override;
  //! Gets length
  shared_ptr<IntMV> length() noexcept(false) override;
  shared_ptr<ArrayMV> append(MacroValuePtr mv) noexcept(false);
  shared_ptr<IntMV> in(const MacroValuePtr &mv) noexcept(false) override;
  /*! Arguments must be of type IntMV.
    Returns an integer array containing all integers between mv1 and mv2.
    If mv2 < mv1, returns an empty range (for consistency with MATLAB).
  */
  static shared_ptr<ArrayMV> range(const MacroValuePtr &mv1, const MacroValuePtr &mv2) noexcept(false);
  shared_ptr<ArrayMV> set_union(const MacroValuePtr &mvp) noexcept(false) override;
  shared_ptr<ArrayMV> set_intersection(const MacroValuePtr &mvp) noexcept(false) override;
  // Computes the Cartesian product of two sets
  MacroValuePtr times(const MacroValuePtr &mv) noexcept(false) override;
  // Shortcut for Cartesian product of two sets
  MacroValuePtr power(const MacroValuePtr &mv) noexcept(false) override;
};

//! Represents a tuple value in macro language
class TupleMV : public MacroValue
{
public:
  TupleMV(vector<MacroValuePtr> values_arg);

  //! Underlying vector
  const vector<MacroValuePtr> values;

  shared_ptr<IntMV> is_equal(const MacroValuePtr &mv) override;
  //! Subscripting operator
  /*! Argument must be an ArrayMV<int>. Indexes begin at 1. Returns a StringMV. */
  MacroValuePtr subscript(const MacroValuePtr &mv) noexcept(false) override;
  //! Returns underlying string value
  string toString() override;
  string print() override;
  shared_ptr<IntMV> length() noexcept(false) override;
  shared_ptr<IntMV> in(const MacroValuePtr &mv) noexcept(false) override;
};

#endif
