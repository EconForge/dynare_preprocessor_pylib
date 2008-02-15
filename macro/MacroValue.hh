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

class MacroValue
{
public:
  class TypeError
  {
  public:
    const string message;
    TypeError(const string &message_arg) : message(message_arg) {};
  };
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
};

class IntMV : public MacroValue
{
private:
  int value;
public:
  IntMV(int value_arg);
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
};

class StringMV : public MacroValue
{
private:
  string value;
public:
  StringMV(const string &value_arg);
  virtual MacroValue *operator+(const MacroValue &mv) const throw (TypeError);
//   virtual MacroValue *operator<(const MacroValue &mv) const throw (TypeError);
//   virtual MacroValue *operator>(const MacroValue &mv) const throw (TypeError);
//   virtual MacroValue *operator<=(const MacroValue &mv) const throw (TypeError);
//   virtual MacroValue *operator>=(const MacroValue &mv) const throw (TypeError);
  virtual MacroValue *operator==(const MacroValue &mv) const throw (TypeError);
  virtual MacroValue *operator!=(const MacroValue &mv) const throw (TypeError);
// virtual MacroValue *operator[](const MacroValue &mv) const throw (TypeError);
  virtual string toString() const;
  virtual MacroValue *clone() const;
};

/*
class IntArrayMV : public MacroValue
{
private:
  vector<int> values;
public:
  IntArrayMV(const vector<int> &values_arg);
  virtual MacroValue *operator+(const MacroValue &mv) const throw (TypeError);
  virtual MacroValue *operator==(const MacroValue &mv) const throw (TypeError);
  virtual MacroValue *operator!=(const MacroValue &mv) const throw (TypeError);
  virtual MacroValue *operator[](const MacroValue &mv) const throw (TypeError);
  virtual string toString() const;
  virtual MacroValue *clone() const;
};

class StringArrayMV : public MacroValue
{
private:
  vector<string> values;
public:
  StringArrayMV(const vector<string> &values_arg);
  virtual MacroValue *operator+(const MacroValue &mv) const throw (TypeError);
  virtual MacroValue *operator==(const MacroValue &mv) const throw (TypeError);
  virtual MacroValue *operator!=(const MacroValue &mv) const throw (TypeError);
  virtual MacroValue *operator[](const MacroValue &mv) const throw (TypeError);
  virtual string toString() const;
  virtual MacroValue *clone() const;
};
*/

#endif
