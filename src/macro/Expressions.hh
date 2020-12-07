/*
 * Copyright © 2019-2020 Dynare Team
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

#ifndef _EXPRESSIONS_HH
#define _EXPRESSIONS_HH

#include "ForwardDeclarationsAndEnums.hh"
#include "Environment.hh"
#include "location.hh"

#include <cmath>
#include <vector>
#include <sstream>
#include <iomanip>

namespace macro
{
  class StackTrace final : public exception
  {
  private:
    vector<string> message;
  public:
    StackTrace (string message_arg) : message{move(message_arg)} { }
    StackTrace (const string &prefix, const char *standard_exception_message, const Tokenizer::location &location)
    {
      stringstream ss;
      ss << prefix << ": " << location << " " << standard_exception_message;
      message = {ss.str()};
    }
    StackTrace (const string &prefix, const string &msg, const Tokenizer::location &location)
    {
      stringstream ss;
      ss << prefix << ": " << location << " " << msg;
      message = {ss.str()};
    }
    void push(const string &prefix, const Tokenizer::location &location)
    {
      stringstream ss;
      auto end_col = 0 < location.end.column ? location.end.column - 1 : 0;

      ss << prefix << ": "
         << R"(")" << *location.begin.filename << R"(" line )" << location.begin.line
         << ", col " << location.begin.column;
      if (location.end.filename
          && (!location.begin.filename
              || *location.begin.filename != *location.end.filename))
        ss << R"( to ")" << location.end.filename << R"(")" << " line " << location.end.line
           << ", col " << end_col;
      else if (location.begin.line < location.end.line)
        ss << " to line " << location.end.line << ", col " << end_col;
      else if (location.begin.column < end_col)
        ss << "-" << end_col;
      message.emplace_back(ss.str());
    }
    string trace() const
    {
      stringstream ss;
      for (auto &msg : message)
        ss << "- " << msg << endl;
      return ss.str();
    }
  };


  class Node
  {
  protected:
    const Tokenizer::location location;
  public:
    explicit Node(Tokenizer::location location_arg) :
      location{move(location_arg)} { }
    virtual ~Node() = default;
  public:
    inline Tokenizer::location getLocation() const noexcept { return location; }
    inline void error(const StackTrace &e) const noexcept
    {
      cerr << endl << "Macro-processing error: backtrace..." << endl << e.trace();
      exit(EXIT_FAILURE);
    }
    inline void warning(const StackTrace &e) const noexcept
    {
      cerr << endl << "Macro-processing warning: backtrace..." << endl << e.trace();
    }
    inline void printLineInfo(ostream &output) const noexcept
    {
      output << R"(@#line ")" << *(location.begin.filename) << R"(" )" << location.begin.line << endl;
    }
    inline void printEndLineInfo(ostream &output) const noexcept
    {
      // Add one to end line because we want to print the line number of the line *following* the end statement
      output << R"(@#line ")" << *(location.begin.filename) << R"(" )" << location.end.line + 1 << endl;
    }
  };


  class Expression : public Node
  {
  public:
    explicit Expression(Tokenizer::location location_arg) :
      Node(move(location_arg)) { }
    virtual string to_string() const noexcept = 0;
    virtual void print(ostream &output, bool matlab_output = false) const noexcept = 0;
    virtual BaseTypePtr eval(Environment &env) = 0;
    virtual ExpressionPtr clone() const noexcept = 0;
  };


  class BaseType : public Expression, public enable_shared_from_this<BaseType>
  {
  public:
    explicit BaseType(Tokenizer::location location_arg = Tokenizer::location()) :
      Expression(move(location_arg)) { }
    virtual codes::BaseType getType() const noexcept = 0;
    inline BaseTypePtr eval(Environment &env) override { return shared_from_this(); }
  public:
    virtual BaseTypePtr plus(const BaseTypePtr &bt) const { throw StackTrace("Operator + does not exist for this type"); }
    virtual BaseTypePtr unary_plus() const { throw StackTrace("Unary operator + does not exist for this type"); }
    virtual BaseTypePtr minus(const BaseTypePtr &bt) const { throw StackTrace("Operator - does not exist for this type"); }
    virtual BaseTypePtr unary_minus() const { throw StackTrace("Unary operator - does not exist for this type"); }
    virtual BaseTypePtr times(const BaseTypePtr &bt) const { throw StackTrace("Operator * does not exist for this type"); }
    virtual BaseTypePtr divide(const BaseTypePtr &bt) const { throw StackTrace("Operator / does not exist for this type"); }
    virtual BaseTypePtr power(const BaseTypePtr &btp) const { throw StackTrace("Operator ^ does not exist for this type"); }
    virtual BoolPtr is_less(const BaseTypePtr &btp) const { throw StackTrace("Operator < does not exist for this type"); }
    virtual BoolPtr is_greater(const BaseTypePtr &btp) const { throw StackTrace("Operator > does not exist for this type"); }
    virtual BoolPtr is_less_equal(const BaseTypePtr &btp) const { throw StackTrace("Operator <= does not exist for this type"); }
    virtual BoolPtr is_greater_equal(const BaseTypePtr &btp) const { throw StackTrace("Operator >= does not exist for this type"); }
    virtual BoolPtr is_equal(const BaseTypePtr &btp) const = 0;
    virtual BoolPtr is_different(const BaseTypePtr &btp) const final;
    virtual BoolPtr logical_and(const ExpressionPtr &ep, Environment &env) const { throw StackTrace("Operator && does not exist for this type"); }
    virtual BoolPtr logical_or(const ExpressionPtr &ep, Environment &env) const { throw StackTrace("Operator || does not exist for this type"); }
    virtual BoolPtr logical_not() const { throw StackTrace("Operator ! does not exist for this type"); }
    virtual ArrayPtr set_union(const BaseTypePtr &btp) const { throw StackTrace("Operator | does not exist for this type"); }
    virtual ArrayPtr set_intersection(const BaseTypePtr &btp) const { throw StackTrace("Operator & does not exist for this type"); }
    virtual BoolPtr contains(const BaseTypePtr &btp) const { throw StackTrace("Second argument of `in` operator must be an array"); }
    virtual RealPtr length() const { throw StackTrace("Operator `length` does not exist for this type"); }
    virtual BoolPtr isempty() const { throw StackTrace("Operator `isempty` does not exist for this type"); }
    virtual BoolPtr isboolean() const noexcept { return make_shared<Bool>(false, location); }
    virtual BoolPtr isreal() const noexcept { return make_shared<Bool>(false, location); }
    virtual BoolPtr isinteger() const noexcept { return make_shared<Bool>(false, location); }
    virtual BoolPtr isstring() const noexcept { return make_shared<Bool>(false, location); }
    virtual BoolPtr istuple() const noexcept { return make_shared<Bool>(false, location); }
    virtual BoolPtr isarray() const noexcept { return make_shared<Bool>(false, location); }
    virtual RealPtr max(const BaseTypePtr &btp) const { throw StackTrace("Operator `max` does not exist for this type"); }
    virtual RealPtr min(const BaseTypePtr &btp) const { throw StackTrace("Operator `min` does not exist for this type"); }
    virtual RealPtr mod(const BaseTypePtr &btp) const { throw StackTrace("Operator `mod` does not exist for this type"); }
    virtual RealPtr exp() const { throw StackTrace("Operator `exp` does not exist for this type"); }
    virtual RealPtr ln() const { throw StackTrace("Operator `ln` does not exist for this type"); }
    virtual RealPtr log10() const { throw StackTrace("Operator `log10` does not exist for this type"); }
    virtual BoolPtr isinf() const { throw StackTrace("Operator `isinf` does not exist for this type"); }
    virtual BoolPtr isnan() const { throw StackTrace("Operator `isnan` does not exist for this type"); }
    virtual BoolPtr isfinite() const { throw StackTrace("Operator `isfinite` does not exist for this type"); }
    virtual BoolPtr isnormal() const { throw StackTrace("Operator `isnormal` does not exist for this type"); }
    virtual RealPtr sin() const { throw StackTrace("Operator `sin` does not exist for this type"); }
    virtual RealPtr cos() const { throw StackTrace("Operator `cos` does not exist for this type"); }
    virtual RealPtr tan() const { throw StackTrace("Operator `tan` does not exist for this type"); }
    virtual RealPtr asin() const { throw StackTrace("Operator `asin` does not exist for this type"); }
    virtual RealPtr acos() const { throw StackTrace("Operator `acos` does not exist for this type"); }
    virtual RealPtr atan() const { throw StackTrace("Operator `atan` does not exist for this type"); }
    virtual RealPtr sqrt() const { throw StackTrace("Operator `sqrt` does not exist for this type"); }
    virtual RealPtr cbrt() const { throw StackTrace("Operator `cbrt` does not exist for this type"); }
    virtual RealPtr sign() const { throw StackTrace("Operator `sign` does not exist for this type"); }
    virtual RealPtr floor() const { throw StackTrace("Operator `floor` does not exist for this type"); }
    virtual RealPtr ceil() const { throw StackTrace("Operator `ceil` does not exist for this type"); }
    virtual RealPtr trunc() const { throw StackTrace("Operator `trunc` does not exist for this type"); }
    virtual RealPtr sum() const { throw StackTrace("Operator `sum` does not exist for this type"); }
    virtual RealPtr erf() const { throw StackTrace("Operator `erf` does not exist for this type"); }
    virtual RealPtr erfc() const { throw StackTrace("Operator `erfc` does not exist for this type"); }
    virtual RealPtr gamma() const { throw StackTrace("Operator `gamma` does not exist for this type"); }
    virtual RealPtr lgamma() const { throw StackTrace("Operator `lgamma` does not exist for this type"); }
    virtual RealPtr round() const { throw StackTrace("Operator `round` does not exist for this type"); }
    virtual RealPtr normpdf() const { throw StackTrace("Operator `normpdf` does not exist for this type"); }
    virtual RealPtr normpdf(const BaseTypePtr &btp1, const BaseTypePtr &btp2) const { throw StackTrace("Operator `normpdf` does not exist for this type"); }
    virtual RealPtr normcdf() const { throw StackTrace("Operator `normcdf` does not exist for this type"); }
    virtual RealPtr normcdf(const BaseTypePtr &btp1, const BaseTypePtr &btp2) const { throw StackTrace("Operator `normcdf` does not exist for this type"); }
    virtual BoolPtr cast_bool(Environment &env) const { throw StackTrace("This type cannot be cast to a boolean"); }
    virtual RealPtr cast_real(Environment &env) const { throw StackTrace("This type cannot be cast to a real"); }
    virtual StringPtr cast_string() const { throw StackTrace("This type cannot be cast to a string"); }
    virtual TuplePtr cast_tuple() const { throw StackTrace("This type cannot be cast to a tuple"); }
    virtual ArrayPtr cast_array() const { throw StackTrace("This type cannot be cast to an array"); }
    virtual BoolPtr defined(const Environment &env) const { throw StackTrace("Operator `defined` does not exist for this type"); }
  };


  class Bool final : public BaseType
  {
  private:
    const bool value;
  public:
    Bool(bool value_arg,
         Tokenizer::location location_arg = Tokenizer::location()) :
      BaseType(move(location_arg)),
      value{value_arg} { }
    inline codes::BaseType getType() const noexcept override { return codes::BaseType::Bool; }
    inline string to_string() const noexcept override { return value ? "true" : "false"; }
    inline void print(ostream &output, bool matlab_output = false) const noexcept override { output << to_string(); }
    inline ExpressionPtr clone() const noexcept override { return make_shared<Bool>(value, location); }
  public:
    operator bool() const { return value; }
    BoolPtr is_equal(const BaseTypePtr &btp) const override;
    BoolPtr logical_and(const ExpressionPtr &ep, Environment &env) const override;
    BoolPtr logical_or(const ExpressionPtr &ep, Environment &env) const override;
    BoolPtr logical_not() const override;
    inline BoolPtr isboolean() const noexcept override { return make_shared<Bool>(true, location); }
    inline BoolPtr cast_bool(Environment &env) const override { return make_shared<Bool>(value); }
    inline RealPtr cast_real(Environment &env) const override { return value ? make_shared<Real>(1) : make_shared<Real>(0); }
    inline StringPtr cast_string() const override { return make_shared<String>(this->to_string()); }
    inline TuplePtr cast_tuple() const override
    {
      return make_shared<Tuple>(vector<ExpressionPtr>{make_shared<Bool>(value)});
    }
    inline ArrayPtr cast_array() const override
    {
      return make_shared<Array>(vector<ExpressionPtr>{make_shared<Bool>(value)});
    }
  };


  class Real final : public BaseType
  {
  private:
    const double value;
  public:
    // Use strtod to handle extreme cases (e.g. 1e500, 1e-500), nan, inf
    // See Note in NumericalConstants::AddNonNegativeConstant
    Real(const string &value_arg,
         Tokenizer::location location_arg = Tokenizer::location()) :
      BaseType(move(location_arg)),
      value{strtod(value_arg.c_str(), nullptr)} { }
    Real(double value_arg,
         Tokenizer::location location_arg = Tokenizer::location()) :
      BaseType(move(location_arg)),
      value{value_arg} { }
    inline codes::BaseType getType() const noexcept override { return codes::BaseType::Real; }
    inline string to_string() const noexcept override
    {
      ostringstream strs;
      strs << setprecision(15) << value;
      return strs.str();
    }
    inline void print(ostream &output, bool matlab_output = false) const noexcept override { output << to_string(); }
    inline ExpressionPtr clone() const noexcept override { return make_shared<Real>(value, location); }
  public:
    operator double() const { return value; }
    BaseTypePtr plus(const BaseTypePtr &bt) const override;
    inline BaseTypePtr unary_plus() const override { return make_shared<Real>(value); }
    BaseTypePtr minus(const BaseTypePtr &bt) const override;
    inline BaseTypePtr unary_minus() const override { return make_shared<Real>(-value); }
    BaseTypePtr times(const BaseTypePtr &bt) const override;
    BaseTypePtr divide(const BaseTypePtr &bt) const override;
    BaseTypePtr power(const BaseTypePtr &btp) const override;
    BoolPtr is_less(const BaseTypePtr &btp) const override;
    BoolPtr is_greater(const BaseTypePtr &btp) const override;
    BoolPtr is_less_equal(const BaseTypePtr &btp) const override;
    BoolPtr is_greater_equal(const BaseTypePtr &btp) const override;
    BoolPtr is_equal(const BaseTypePtr &btp) const override;
    inline BoolPtr isreal() const noexcept override { return make_shared<Bool>(true, location); }
    inline BoolPtr isinteger() const noexcept override
    {
      double intpart;
      return make_shared<Bool>(modf(value, &intpart) == 0.0, location);
    }
    BoolPtr logical_and(const ExpressionPtr &ep, Environment &env) const override;
    BoolPtr logical_or(const ExpressionPtr &ep, Environment &env) const override;
    BoolPtr logical_not() const override;
    RealPtr max(const BaseTypePtr &btp) const override;
    RealPtr min(const BaseTypePtr &btp) const override;
    RealPtr mod(const BaseTypePtr &btp) const override;
    inline RealPtr exp() const override { return make_shared<Real>(std::exp(value)); }
    inline RealPtr ln() const override { return make_shared<Real>(std::log(value)); }
    inline RealPtr log10() const override { return make_shared<Real>(std::log10(value)); }
    inline BoolPtr isinf() const override { return make_shared<Bool>(std::isinf(value)); }
    inline BoolPtr isnan() const override { return make_shared<Bool>(std::isnan(value)); }
    inline BoolPtr isfinite() const override { return make_shared<Bool>(std::isfinite(value)); }
    inline BoolPtr isnormal() const override { return make_shared<Bool>(std::isnormal(value)); }
    inline RealPtr sin() const override { return make_shared<Real>(std::sin(value)); }
    inline RealPtr cos() const override { return make_shared<Real>(std::cos(value)); }
    inline RealPtr tan() const override { return make_shared<Real>(std::tan(value)); }
    inline RealPtr asin() const override { return make_shared<Real>(std::asin(value)); }
    inline RealPtr acos() const override { return make_shared<Real>(std::acos(value)); }
    inline RealPtr atan() const override { return make_shared<Real>(std::atan(value)); }
    inline RealPtr sqrt() const override { return make_shared<Real>(std::sqrt(value)); }
    inline RealPtr cbrt() const override { return make_shared<Real>(std::cbrt(value)); }
    inline RealPtr sign() const override
    {
      return make_shared<Real>((value > 0) ? 1. : ((value < 0) ? -1. : 0.));
    }
    inline RealPtr floor() const override { return make_shared<Real>(std::floor(value)); }
    inline RealPtr ceil() const override { return make_shared<Real>(std::ceil(value)); }
    inline RealPtr trunc() const override { return make_shared<Real>(std::trunc(value)); }
    inline RealPtr erf() const override { return make_shared<Real>(std::erf(value)); }
    inline RealPtr erfc() const override { return make_shared<Real>(std::erfc(value)); }
    inline RealPtr gamma() const override { return make_shared<Real>(std::tgamma(value)); }
    inline RealPtr lgamma() const override { return make_shared<Real>(std::lgamma(value)); }
    inline RealPtr round() const override { return make_shared<Real>(std::round(value)); }
    inline RealPtr normpdf() const override
    {
      return normpdf(make_shared<Real>(0), make_shared<Real>(1));
    }
    RealPtr normpdf(const BaseTypePtr &btp1, const BaseTypePtr &btp2) const override;
    inline RealPtr normcdf() const override
    {
      return normcdf(make_shared<Real>(0), make_shared<Real>(1));
    }
    RealPtr normcdf(const BaseTypePtr &btp1, const BaseTypePtr &btp2) const override;
    inline BoolPtr cast_bool(Environment &env) const override { return make_shared<Bool>(static_cast<bool>(value)); }
    inline RealPtr cast_real(Environment &env) const override { return make_shared<Real>(value); }
    inline StringPtr cast_string() const override { return make_shared<String>(this->to_string()); }
    inline TuplePtr cast_tuple() const override
    {
      return make_shared<Tuple>(vector<ExpressionPtr>{make_shared<Real>(value)});
    }
    inline ArrayPtr cast_array() const override
    {
      return make_shared<Array>(vector<ExpressionPtr>{make_shared<Real>(value)});
    }
  };

  class String final : public BaseType
  {
  private:
    const string value;
  public:
    String(string value_arg,
           Tokenizer::location location_arg = Tokenizer::location()) :
      BaseType(move(location_arg)),
      value{move(value_arg)} { }
    inline codes::BaseType getType() const noexcept override { return codes::BaseType::String; }
    inline string to_string() const noexcept override { return value; }
    void print(ostream &output, bool matlab_output = false) const noexcept override;
    inline ExpressionPtr clone() const noexcept override { return make_shared<String>(value, location); }
  public:
    operator string() const { return value; }
    BaseTypePtr plus(const BaseTypePtr &bt) const override;
    BoolPtr is_less(const BaseTypePtr &btp) const override;
    BoolPtr is_greater(const BaseTypePtr &btp) const override;
    BoolPtr is_less_equal(const BaseTypePtr &btp) const override;
    BoolPtr is_greater_equal(const BaseTypePtr &btp) const override;
    BoolPtr is_equal(const BaseTypePtr &btp) const override;
    inline BoolPtr isstring() const noexcept override { return make_shared<Bool>(true, location); }
    inline RealPtr length() const override { return make_shared<Real>(value.size()); }
    inline BoolPtr isempty() const override { return make_shared<Bool>(value.empty()); }
    BoolPtr cast_bool(Environment &env) const override;
    RealPtr cast_real(Environment &env) const override;
    inline StringPtr cast_string() const override { return make_shared<String>(value); }
    inline TuplePtr cast_tuple() const override
    {
      return make_shared<Tuple>(vector<ExpressionPtr>{make_shared<String>(value)});
    }
    inline ArrayPtr cast_array() const override
    {
      return make_shared<Array>(vector<ExpressionPtr>{make_shared<String>(value)});
    }
    inline BoolPtr defined(const Environment &env) const override
    {
      return make_shared<Bool>(env.isSymbolDefined(value));
    }
  };


  class Tuple final : public BaseType
  {
  private:
    const vector<ExpressionPtr> tup;
  public:
    Tuple(vector<ExpressionPtr> tup_arg,
          Tokenizer::location location_arg = Tokenizer::location()) :
      BaseType(move(location_arg)),
      tup{move(tup_arg)} { }
    inline codes::BaseType getType() const noexcept override { return codes::BaseType::Tuple; }
    string to_string() const noexcept override;
    void print(ostream &output, bool matlab_output = false) const noexcept override;
    BaseTypePtr eval(Environment &env) override;
    ExpressionPtr clone() const noexcept override;
  public:
    inline size_t size() const { return tup.size(); }
    inline bool empty() const { return tup.empty(); }
    inline const vector<ExpressionPtr> &getValue() const { return tup; }
    inline const ExpressionPtr &at(int i) const { return tup.at(i); }
    BoolPtr is_equal(const BaseTypePtr &btp) const override;
    inline BoolPtr istuple() const noexcept override { return make_shared<Bool>(true, location); }
    BoolPtr contains(const BaseTypePtr &btp) const override;
    inline RealPtr length() const override { return make_shared<Real>(tup.size()); }
    inline BoolPtr isempty() const override { return make_shared<Bool>(empty()); }
    BoolPtr cast_bool(Environment &env) const override;
    RealPtr cast_real(Environment &env) const override;
    inline StringPtr cast_string() const override { return make_shared<String>(this->to_string()); }
    inline TuplePtr cast_tuple() const override { return make_shared<Tuple>(tup); }
    inline ArrayPtr cast_array() const override { return make_shared<Array>(tup); }
  };


  class Array final : public BaseType
  {
  private:
    const vector<ExpressionPtr> arr;
  public:
    Array(vector<ExpressionPtr> arr_arg,
          Tokenizer::location location_arg = Tokenizer::location()) :
      BaseType(move(location_arg)), arr{move(arr_arg)} { }
    inline codes::BaseType getType() const noexcept override { return codes::BaseType::Array; }
    string to_string() const noexcept override;
    void print(ostream &output, bool matlab_output = false) const noexcept override;
    BaseTypePtr eval(Environment &env) override;
    ExpressionPtr clone() const noexcept override;
  public:
    inline size_t size() const { return arr.size(); }
    inline const vector<ExpressionPtr> &getValue() const { return arr; }
    inline const ExpressionPtr &at(int i) const { return arr.at(i); }
    inline bool empty() const { return arr.empty(); }
    BaseTypePtr plus(const BaseTypePtr &bt) const override;
    BaseTypePtr minus(const BaseTypePtr &bt) const override;
    BaseTypePtr times(const BaseTypePtr &bt) const override;
    BaseTypePtr power(const BaseTypePtr &btp) const override;
    BoolPtr is_equal(const BaseTypePtr &btp) const override;
    inline BoolPtr isarray() const noexcept override { return make_shared<Bool>(true, location); }
    ArrayPtr set_union(const BaseTypePtr &btp) const override;
    ArrayPtr set_intersection(const BaseTypePtr &btp) const override;
    BoolPtr contains(const BaseTypePtr &btp) const override;
    inline RealPtr length() const override { return make_shared<Real>(arr.size()); }
    inline BoolPtr isempty() const override { return make_shared<Bool>(empty()); }
    RealPtr sum() const override;
    BoolPtr cast_bool(Environment &env) const override;
    RealPtr cast_real(Environment &env) const override;
    inline StringPtr cast_string() const override { return make_shared<String>(this->to_string()); }
    inline TuplePtr cast_tuple() const override { return make_shared<Tuple>(arr); }
    inline ArrayPtr cast_array() const override { return make_shared<Array>(arr); }
  };


  class Range final : public BaseType
  {
  private:
    const ExpressionPtr start, inc, end;
  public:
    Range(ExpressionPtr start_arg, ExpressionPtr end_arg,
          Tokenizer::location location_arg) :
      BaseType(move(location_arg)), start{move(start_arg)}, end{move(end_arg)} { }
    Range(ExpressionPtr start_arg, ExpressionPtr inc_arg, ExpressionPtr end_arg,
          Tokenizer::location location_arg) :
      BaseType(move(location_arg)),
      start{move(start_arg)}, inc{move(inc_arg)}, end{move(end_arg)} { }
    inline codes::BaseType getType() const noexcept override { return codes::BaseType::Range; }
    inline string to_string() const noexcept override
    {
      string retval = "[" + start->to_string() + ":";
      if (inc)
        retval += inc->to_string() + ":";
      return retval + end->to_string() + "]";
    }
    inline void print(ostream &output, bool matlab_output = false) const noexcept override { output << to_string(); }
    BaseTypePtr eval(Environment &env) override;
    inline ExpressionPtr clone() const noexcept override
    {
      return inc ?
        make_shared<Range>(start, inc, end, location)
        : make_shared<Range>(start, end, location);
    }
  public:
    inline BoolPtr is_equal(const BaseTypePtr &btp) const override
    {
      throw StackTrace("Internal error: Range: Should not arrive here: is_equal");
    }
  };


  class Variable final : public Expression
  {
  private:
    const string name;
    const ArrayPtr indices; // for indexing strings/arrays
  public:
    Variable(string name_arg,
             Tokenizer::location location_arg) :
      Expression(move(location_arg)), name{move(name_arg)} { }
    Variable(string name_arg, ArrayPtr indices_arg,
             Tokenizer::location location_arg) :
      Expression(move(location_arg)), name{move(name_arg)}, indices{move(indices_arg)} { }
    inline string to_string() const noexcept override { return name; }
    inline void print(ostream &output, bool matlab_output = false) const noexcept override { output << name; }
    BaseTypePtr eval(Environment &env) override;
    inline ExpressionPtr clone() const noexcept override
    {
      return indices ? make_shared<Variable>(name, indices, location) :
        make_shared<Variable>(name, location);
    }
  public:
    inline const string &getName() const noexcept { return name; }
    inline codes::BaseType getType(const Environment &env) const { return env.getType(name); }
  };


  class Function final : public Expression
  {
  private:
    const string name;
    const vector<ExpressionPtr> args;
  public:
    Function(string name_arg,
             vector<ExpressionPtr> args_arg,
             Tokenizer::location location_arg) :
      Expression(move(location_arg)), name{move(name_arg)}, args{move(args_arg)} { }
    string to_string() const noexcept override;
    inline void print(ostream &output, bool matlab_output = false) const noexcept override
    {
      printName(output); printArgs(output);
    }
    BaseTypePtr eval(Environment &env) override;
    ExpressionPtr clone() const noexcept override;
  public:
    inline void printName(ostream &output) const noexcept { output << name; }
    void printArgs(ostream &output) const noexcept;
    inline const string &getName() const { return name; }
    inline const vector<ExpressionPtr> &getArgs() const { return args; }
  };


  class UnaryOp final : public Expression
  {
  private:
    const codes::UnaryOp op_code;
    const ExpressionPtr arg;
  public:
    UnaryOp(codes::UnaryOp op_code_arg,
            ExpressionPtr arg_arg,
            Tokenizer::location location_arg) :
      Expression(move(location_arg)), op_code{move(op_code_arg)}, arg{move(arg_arg)} { }
    string to_string() const noexcept override;
    void print(ostream &output, bool matlab_output = false) const noexcept override;
    BaseTypePtr eval(Environment &env) override;
    inline ExpressionPtr clone() const noexcept override
    {
      return make_shared<UnaryOp>(op_code, arg->clone(), location);
    }
  };


  class BinaryOp final : public Expression
  {
  private:
    const codes::BinaryOp op_code;
    const ExpressionPtr arg1, arg2;
  public:
    BinaryOp(codes::BinaryOp op_code_arg,
             ExpressionPtr arg1_arg, ExpressionPtr arg2_arg,
             Tokenizer::location location_arg) :
      Expression(move(location_arg)), op_code{op_code_arg},
      arg1{move(arg1_arg)}, arg2{move(arg2_arg)} { }
  public:
    string to_string() const noexcept override;
    void print(ostream &output, bool matlab_output = false) const noexcept override;
    BaseTypePtr eval(Environment &env) override;
    inline ExpressionPtr clone() const noexcept override
    {
      return make_shared<BinaryOp>(op_code, arg1->clone(), arg2->clone(), location);
    }
  };


  class TrinaryOp final : public Expression
  {
  private:
    const codes::TrinaryOp op_code;
    const ExpressionPtr arg1, arg2, arg3;
  public:
    TrinaryOp(codes::TrinaryOp op_code_arg,
              ExpressionPtr arg1_arg, ExpressionPtr arg2_arg, ExpressionPtr arg3_arg,
              Tokenizer::location location_arg) :
      Expression(move(location_arg)), op_code{op_code_arg},
      arg1{move(arg1_arg)}, arg2{move(arg2_arg)}, arg3{move(arg3_arg)} { }
    string to_string() const noexcept override;
    void print(ostream &output, bool matlab_output = false) const noexcept override;
    BaseTypePtr eval(Environment &env) override;
    inline ExpressionPtr clone() const noexcept override
    {
      return make_shared<TrinaryOp>(op_code, arg1->clone(), arg2->clone(), arg3->clone(), location);
    }
  };


  class Comprehension final : public Expression
  {
    /*
     * Filter:       [c_vars IN c_set WHEN c_when]             => c_expr == nullptr
     * Map:          [c_expr FOR c_vars IN c_set]              => c_when == nullptr
     * Filter + Map: [c_expr FOR c_vars IN c_set WHEN c_when]  => all members assigned
     */
  private:
    const ExpressionPtr c_expr, c_vars, c_set, c_when;
  public:
    Comprehension(ExpressionPtr c_expr_arg,
                  ExpressionPtr c_vars_arg,
                  ExpressionPtr c_set_arg,
                  ExpressionPtr c_when_arg,
                  Tokenizer::location location_arg) :
      Expression(move(location_arg)),
      c_expr{move(c_expr_arg)}, c_vars{move(c_vars_arg)},
      c_set{move(c_set_arg)}, c_when{move(c_when_arg)} { }
    Comprehension(ExpressionPtr c_expr_arg,
                  ExpressionPtr c_vars_arg,
                  ExpressionPtr c_set_arg,
                  Tokenizer::location location_arg) :
      Expression(move(location_arg)),
      c_expr{move(c_expr_arg)}, c_vars{move(c_vars_arg)}, c_set{move(c_set_arg)} { }
    Comprehension(bool filter_only_arg,
                  ExpressionPtr c_vars_arg,
                  ExpressionPtr c_set_arg,
                  ExpressionPtr c_when_arg,
                  Tokenizer::location location_arg) :
      Expression(move(location_arg)),
      c_vars{move(c_vars_arg)}, c_set{move(c_set_arg)}, c_when{move(c_when_arg)} { }
    string to_string() const noexcept override;
    void print(ostream &output, bool matlab_output = false) const noexcept override;
    BaseTypePtr eval(Environment &env) override;
    ExpressionPtr clone() const noexcept override;
  };
}
#endif
