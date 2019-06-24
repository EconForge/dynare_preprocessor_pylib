/*
 * Copyright © 2019 Dynare Team
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

using namespace std;

namespace macro
{
  class StackTrace final : public exception
  {
  private:
    vector<string> message;
  public:
    StackTrace (string message_arg) : message{{move(message_arg)}} { }
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
      unsigned end_col = 0 < location.end.column ? location.end.column - 1 : 0;

      ss << prefix << ": "
         << R"(")" << *location.begin.filename << R"(" line )" << location.begin.line << ", col " << location.begin.column;
      if (location.end.filename
          && (!location.begin.filename
              || *location.begin.filename != *location.end.filename))
        ss << R"( to ")" << location.end.filename << R"(")" << " line " << location.end.line << ", col " << end_col;
      else if (location.begin.line < location.end.line)
        ss << " to line " << location.end.line << ", col " << end_col;
      else if (location.begin.column < end_col)
        ss << "-" << end_col;
      message.emplace_back(ss.str());
    }
    string trace() const
    {
      stringstream ss;
      for (auto & msg : message)
        ss << "- " << msg << endl;
      return ss.str();
    }
  };


  class Node
  {
  protected:
    Environment &env;
    const Tokenizer::location location;
  public:
    Node(Environment &env_arg, Tokenizer::location location_arg) :
      env{env_arg}, location{move(location_arg)} { }
    virtual ~Node() = default;
  public:
    inline Tokenizer::location getLocation() const noexcept { return location; }
    inline void error(const StackTrace &e) const noexcept
    {
      cerr << endl << "Macro-processing error: backtrace..." << endl << e.trace();
      exit(EXIT_FAILURE);
    }
    inline void printLineInfo(ostream &output, bool no_line_macro) const noexcept
    {
      if (!no_line_macro)
        output << R"(@#line ")" << *(location.begin.filename) << R"(" )" << location.begin.line << endl;
    }
    inline void printEndLineInfo(ostream &output, bool no_line_macro) const noexcept
    {
      if (!no_line_macro)
        // Add one to end line because we want to print the line number of the line *following* the end statement
        output << R"(@#line ")" << *(location.begin.filename) << R"(" )" << location.end.line + 1 << endl;
    }
  };


  class Expression : public Node
  {
  public:
    Expression(Environment &env_arg, const Tokenizer::location location_arg) :
      Node(env_arg, move(location_arg)) { }
    virtual string to_string() const noexcept = 0;
    virtual void print(ostream &output, bool matlab_output = false) const noexcept = 0;
    virtual BaseTypePtr eval() = 0;
  };


  class BaseType : public Expression, public enable_shared_from_this<BaseType>
  {
  public:
    BaseType(Environment &env_arg, Tokenizer::location location_arg = Tokenizer::location()) :
      Expression(env_arg, move(location_arg)) { }
    virtual codes::BaseType getType() const noexcept = 0;
    inline BaseTypePtr eval() override { return shared_from_this(); }
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
    virtual BoolPtr logical_and(const BaseTypePtr &btp) const { throw StackTrace("Operator && does not exist for this type"); }
    virtual BoolPtr logical_or(const BaseTypePtr &btp) const { throw StackTrace("Operator || does not exist for this type"); }
    virtual BoolPtr logical_not() const { throw StackTrace("Operator ! does not exist for this type"); }
    virtual ArrayPtr set_union(const BaseTypePtr &btp) const { throw StackTrace("Operator | does not exist for this type"); }
    virtual ArrayPtr set_intersection(const BaseTypePtr &btp) const { throw StackTrace("Operator & does not exist for this type"); }
    virtual BoolPtr contains(const BaseTypePtr &btp) const { throw StackTrace("Second argument of `in` operator must be an array"); }
    virtual DoublePtr length() const { throw StackTrace("Operator `length` does not exist for this type"); }
    virtual DoublePtr max(const BaseTypePtr &btp) const { throw StackTrace("Operator `max` does not exist for this type"); }
    virtual DoublePtr min(const BaseTypePtr &btp) const { throw StackTrace("Operator `min` does not exist for this type"); }
    virtual DoublePtr mod(const BaseTypePtr &btp) const { throw StackTrace("Operator `mod` does not exist for this type"); }
    virtual DoublePtr exp() const { throw StackTrace("Operator `exp` does not exist for this type"); }
    virtual DoublePtr ln() const { throw StackTrace("Operator `ln` does not exist for this type"); }
    virtual DoublePtr log10() const { throw StackTrace("Operator `log10` does not exist for this type"); }
    virtual BoolPtr isinf() const { throw StackTrace("Operator `isinf` does not exist for this type"); }
    virtual BoolPtr isnan() const { throw StackTrace("Operator `isnan` does not exist for this type"); }
    virtual BoolPtr isfinite() const { throw StackTrace("Operator `isfinite` does not exist for this type"); }
    virtual BoolPtr isnormal() const { throw StackTrace("Operator `isnormal` does not exist for this type"); }
    virtual DoublePtr sin() const { throw StackTrace("Operator `sin` does not exist for this type"); }
    virtual DoublePtr cos() const { throw StackTrace("Operator `cos` does not exist for this type"); }
    virtual DoublePtr tan() const { throw StackTrace("Operator `tan` does not exist for this type"); }
    virtual DoublePtr asin() const { throw StackTrace("Operator `asin` does not exist for this type"); }
    virtual DoublePtr acos() const { throw StackTrace("Operator `acos` does not exist for this type"); }
    virtual DoublePtr atan() const { throw StackTrace("Operator `atan` does not exist for this type"); }
    virtual DoublePtr sqrt() const { throw StackTrace("Operator `sqrt` does not exist for this type"); }
    virtual DoublePtr cbrt() const { throw StackTrace("Operator `cbrt` does not exist for this type"); }
    virtual DoublePtr sign() const { throw StackTrace("Operator `sign` does not exist for this type"); }
    virtual DoublePtr floor() const { throw StackTrace("Operator `floor` does not exist for this type"); }
    virtual DoublePtr ceil() const { throw StackTrace("Operator `ceil` does not exist for this type"); }
    virtual DoublePtr trunc() const { throw StackTrace("Operator `trunc` does not exist for this type"); }
    virtual DoublePtr sum() const { throw StackTrace("Operator `sum` does not exist for this type"); }
    virtual DoublePtr erf() const { throw StackTrace("Operator `erf` does not exist for this type"); }
    virtual DoublePtr erfc() const { throw StackTrace("Operator `erfc` does not exist for this type"); }
    virtual DoublePtr gamma() const { throw StackTrace("Operator `gamma` does not exist for this type"); }
    virtual DoublePtr lgamma() const { throw StackTrace("Operator `lgamma` does not exist for this type"); }
    virtual DoublePtr round() const { throw StackTrace("Operator `round` does not exist for this type"); }
    virtual DoublePtr normpdf() const { throw StackTrace("Operator `normpdf` does not exist for this type"); }
    virtual DoublePtr normpdf(const BaseTypePtr &btp1, const BaseTypePtr &btp2) const { throw StackTrace("Operator `normpdf` does not exist for this type"); }
    virtual DoublePtr normcdf() const { throw StackTrace("Operator `normcdf` does not exist for this type"); }
    virtual DoublePtr normcdf(const BaseTypePtr &btp1, const BaseTypePtr &btp2) const { throw StackTrace("Operator `normcdf` does not exist for this type"); }
  };


  class Bool final : public BaseType
  {
  private:
    bool value;
  public:
    Bool(const bool value_arg,
         Environment &env_arg, Tokenizer::location location_arg = Tokenizer::location()) :
      BaseType(env_arg, move(location_arg)),
      value{value_arg} { }
    inline codes::BaseType getType() const noexcept override { return codes::BaseType::Bool; }
    inline string to_string() const noexcept override { return value ? "true" : "false"; }
    inline void print(ostream &output, bool matlab_output = false) const noexcept override { output << to_string(); }
  public:
    operator bool() const { return value; }
    BoolPtr is_equal(const BaseTypePtr &btp) const override;
    BoolPtr logical_and(const BaseTypePtr &btp) const override;
    BoolPtr logical_or(const BaseTypePtr &btp) const override;
    BoolPtr logical_not() const override;
  };


  class Double final : public BaseType
  {
  private:
    double value;
  public:
    // Use strtod to handle extreme cases (e.g. 1e500, 1e-500), nan, inf
    // See Note in NumericalConstants::AddNonNegativeConstant
    Double(const string value_arg,
           Environment &env_arg, Tokenizer::location location_arg = Tokenizer::location()) :
      BaseType(env_arg, move(location_arg)),
      value{strtod(value_arg.c_str(), nullptr)} { }
    Double(double value_arg,
           Environment &env_arg) :
      BaseType(env_arg),
      value{value_arg} { }
    inline codes::BaseType getType() const noexcept override { return codes::BaseType::Double; }
    inline string to_string() const noexcept override
    {
      ostringstream strs;
      strs << setprecision(15) << value;
      return strs.str();
    }
    inline void print(ostream &output, bool matlab_output = false) const noexcept override { output << to_string(); }
  public:
    operator double() const { return value; }
    BaseTypePtr plus(const BaseTypePtr &bt) const override;
    inline BaseTypePtr unary_plus() const override { return make_shared<Double>(value, env); }
    BaseTypePtr minus(const BaseTypePtr &bt) const override;
    inline BaseTypePtr unary_minus() const override { return make_shared<Double>(-value, env); }
    BaseTypePtr times(const BaseTypePtr &bt) const override;
    BaseTypePtr divide(const BaseTypePtr &bt) const override;
    BaseTypePtr power(const BaseTypePtr &btp) const override;
    BoolPtr is_less(const BaseTypePtr &btp) const override;
    BoolPtr is_greater(const BaseTypePtr &btp) const override;
    BoolPtr is_less_equal(const BaseTypePtr &btp) const override;
    BoolPtr is_greater_equal(const BaseTypePtr &btp) const override;
    BoolPtr is_equal(const BaseTypePtr &btp) const override;
    BoolPtr logical_and(const BaseTypePtr &btp) const override;
    BoolPtr logical_or(const BaseTypePtr &btp) const override;
    BoolPtr logical_not() const override;
    DoublePtr max(const BaseTypePtr &btp) const override;
    DoublePtr min(const BaseTypePtr &btp) const override;
    DoublePtr mod(const BaseTypePtr &btp) const override;
    inline DoublePtr exp() const override { return make_shared<Double>(std::exp(value), env); }
    inline DoublePtr ln() const override { return make_shared<Double>(std::log(value), env); }
    inline DoublePtr log10() const override { return make_shared<Double>(std::log10(value), env); }
    inline BoolPtr isinf() const override { return make_shared<Bool>(std::isinf(value), env); }
    inline BoolPtr isnan() const override { return make_shared<Bool>(std::isnan(value), env); }
    inline BoolPtr isfinite() const override { return make_shared<Bool>(std::isfinite(value), env); }
    inline BoolPtr isnormal() const override { return make_shared<Bool>(std::isnormal(value), env); }
    inline DoublePtr sin() const override { return make_shared<Double>(std::sin(value), env); }
    inline DoublePtr cos() const override { return make_shared<Double>(std::cos(value), env); }
    inline DoublePtr tan() const override { return make_shared<Double>(std::tan(value), env); }
    inline DoublePtr asin() const override { return make_shared<Double>(std::asin(value), env); }
    inline DoublePtr acos() const override { return make_shared<Double>(std::acos(value), env); }
    inline DoublePtr atan() const override { return make_shared<Double>(std::atan(value), env); }
    inline DoublePtr sqrt() const override { return make_shared<Double>(std::sqrt(value), env); }
    inline DoublePtr cbrt() const override { return make_shared<Double>(std::cbrt(value), env); }
    inline DoublePtr sign() const override { return make_shared<Double>((value > 0) ? 1. : ((value < 0) ? -1. : 0.), env); }
    inline DoublePtr floor() const override { return make_shared<Double>(std::floor(value), env); }
    inline DoublePtr ceil() const override { return make_shared<Double>(std::ceil(value), env); }
    inline DoublePtr trunc() const override { return make_shared<Double>(std::trunc(value), env); }
    inline DoublePtr erf() const override { return make_shared<Double>(std::erf(value), env); }
    inline DoublePtr erfc() const override { return make_shared<Double>(std::erfc(value), env); }
    inline DoublePtr gamma() const override { return make_shared<Double>(std::tgamma(value), env); }
    inline DoublePtr lgamma() const override { return make_shared<Double>(std::lgamma(value), env); }
    inline DoublePtr round() const override { return make_shared<Double>(std::round(value), env); }
    inline DoublePtr normpdf() const override { return normpdf(make_shared<Double>(0, env), make_shared<Double>(1, env)); }
    DoublePtr normpdf(const BaseTypePtr &btp1, const BaseTypePtr &btp2) const override;
    inline DoublePtr normcdf() const override { return normcdf(make_shared<Double>(0, env), make_shared<Double>(1, env)); }
    DoublePtr normcdf(const BaseTypePtr &btp1, const BaseTypePtr &btp2) const override;
  };

  class String final : public BaseType
  {
  private:
    string value;
  public:
    String(const string value_arg,
           Environment &env_arg, Tokenizer::location location_arg = Tokenizer::location()) :
      BaseType(env_arg, move(location_arg)),
      value{move(value_arg)} { }
    inline codes::BaseType getType() const noexcept override { return codes::BaseType::String; }
    inline string to_string() const noexcept override { return value; }
    void print(ostream &output, bool matlab_output = false) const noexcept override;
  public:
    operator string() const { return value; }
    BaseTypePtr plus(const BaseTypePtr &bt) const override;
    BoolPtr is_less(const BaseTypePtr &btp) const override;
    BoolPtr is_greater(const BaseTypePtr &btp) const override;
    BoolPtr is_less_equal(const BaseTypePtr &btp) const override;
    BoolPtr is_greater_equal(const BaseTypePtr &btp) const override;
    BoolPtr is_equal(const BaseTypePtr &btp) const override;
    inline DoublePtr length() const override { return make_shared<Double>(value.size(), env); }
  };


  class Tuple final : public BaseType
  {
  private:
    vector<ExpressionPtr> tup;
  public:
    Tuple(const vector<ExpressionPtr> tup_arg,
          Environment &env_arg, Tokenizer::location location_arg = Tokenizer::location()) :
      BaseType(env_arg, move(location_arg)),
      tup{move(tup_arg)} { }
    inline codes::BaseType getType() const noexcept override { return codes::BaseType::Tuple; }
    string to_string() const noexcept override;
    void print(ostream &output, bool matlab_output = false) const noexcept override;
    BaseTypePtr eval() override;
  public:
    inline size_t size() const { return tup.size(); }
    inline bool empty() const { return tup.empty(); }
    inline vector<ExpressionPtr> getValue() const { return tup; }
    inline ExpressionPtr at(int i) const { return tup.at(i); }
    BoolPtr is_equal(const BaseTypePtr &btp) const override;
    BoolPtr contains(const BaseTypePtr &btp) const override;
    inline DoublePtr length() const override { return make_shared<Double>(tup.size(), env); }
  };


  class Array final : public BaseType
  {
  private:
    vector<ExpressionPtr> arr;
    ExpressionPtr range1, increment, range2;
  public:
    Array(const vector<ExpressionPtr> arr_arg,
          Environment &env_arg, Tokenizer::location location_arg = Tokenizer::location()) :
      BaseType(env_arg, move(location_arg)),
      arr{move(arr_arg)} { }
    Array(const ExpressionPtr range1_arg, const ExpressionPtr range2_arg,
          Environment &env_arg, Tokenizer::location location_arg) :
      BaseType(env_arg, move(location_arg)),
      range1{move(range1_arg)}, range2{move(range2_arg)} { }
    Array(const ExpressionPtr range1_arg, const ExpressionPtr increment_arg, const ExpressionPtr range2_arg,
          Environment &env_arg, Tokenizer::location location_arg) :
      BaseType(env_arg, move(location_arg)),
      range1{move(range1_arg)}, increment{move(increment_arg)}, range2{move(range2_arg)} { }
    inline codes::BaseType getType() const noexcept override { return codes::BaseType::Array; }
    string to_string() const noexcept override;
    void print(ostream &output, bool matlab_output = false) const noexcept override;
    BaseTypePtr eval() override;
  private:
    BaseTypePtr evalArray();
  public:
    inline size_t size() const { return arr.size(); }
    inline vector<ExpressionPtr> getValue() const { return arr; }
    inline ExpressionPtr at(int i) const { return arr.at(i); }
    inline bool empty() const { return arr.empty() && !range1 && !range2; }
    BaseTypePtr plus(const BaseTypePtr &bt) const override;
    BaseTypePtr minus(const BaseTypePtr &bt) const override;
    BaseTypePtr times(const BaseTypePtr &bt) const override;
    BaseTypePtr power(const BaseTypePtr &btp) const override;
    BoolPtr is_equal(const BaseTypePtr &btp) const override;
    ArrayPtr set_union(const BaseTypePtr &btp) const override;
    ArrayPtr set_intersection(const BaseTypePtr &btp) const override;
    BoolPtr contains(const BaseTypePtr &btp) const override;
    inline DoublePtr length() const override { return make_shared<Double>(arr.size(), env); }
    DoublePtr sum() const override;
  };


  class Variable final : public Expression
  {
  private:
    const string name;
    ArrayPtr indices; // for strings/arrays
  public:
    Variable(const string name_arg,
             Environment &env_arg, const Tokenizer::location location_arg) :
      Expression(env_arg, move(location_arg)), name{move(name_arg)} { }
    inline void addIndexing(const vector<ExpressionPtr> indices_arg) { indices = make_shared<Array>(indices_arg, env); }
    inline string to_string() const noexcept override { return name; }
    inline void print(ostream &output, bool matlab_output = false) const noexcept override { output << name; }
    BaseTypePtr eval() override;
  public:
    inline string getName() const noexcept { return name; }
    inline codes::BaseType getType() const { return env.getType(name); }
  };


  class Function final : public Expression
  {
  private:
    const string name;
    const vector<ExpressionPtr> args;
  public:
    Function(const string &name_arg,
             const vector<ExpressionPtr> &args_arg,
             Environment &env_arg, const Tokenizer::location location_arg) :
      Expression(env_arg, move(location_arg)), name{move(name_arg)}, args{move(args_arg)} { }
    string to_string() const noexcept override;
    inline void print(ostream &output, bool matlab_output = false) const noexcept override { printName(output); printArgs(output); }
    BaseTypePtr eval() override;
  public:
    inline void printName(ostream &output) const noexcept { output << name; }
    void printArgs(ostream &output) const noexcept;
    inline string getName() const { return name; }
    inline vector<ExpressionPtr> getArgs() const { return args; }
  };


  class UnaryOp final : public Expression
  {
  private:
    const codes::UnaryOp op_code;
    const ExpressionPtr arg;
  public:
    UnaryOp(codes::UnaryOp op_code_arg,
            const ExpressionPtr arg_arg,
            Environment &env_arg, const Tokenizer::location location_arg) :
      Expression(env_arg, move(location_arg)), op_code{move(op_code_arg)}, arg{move(arg_arg)} { }
    string to_string() const noexcept override;
    void print(ostream &output, bool matlab_output = false) const noexcept override;
    BaseTypePtr eval() override;
  };


  class BinaryOp final : public Expression
  {
  private:
    const codes::BinaryOp op_code;
    const ExpressionPtr arg1, arg2;
  public:
    BinaryOp(codes::BinaryOp op_code_arg,
             const ExpressionPtr arg1_arg, const ExpressionPtr arg2_arg,
             Environment &env_arg, const Tokenizer::location location_arg) :
      Expression(env_arg, move(location_arg)), op_code{move(op_code_arg)},
      arg1{move(arg1_arg)}, arg2{move(arg2_arg)} { }
  public:
    string to_string() const noexcept override;
    void print(ostream &output, bool matlab_output = false) const noexcept override;
    BaseTypePtr eval() override;
  };


  class TrinaryOp final : public Expression
  {
  private:
    const codes::TrinaryOp op_code;
    const ExpressionPtr arg1, arg2, arg3;
  public:
    TrinaryOp(codes::TrinaryOp op_code_arg,
              const ExpressionPtr arg1_arg, const ExpressionPtr arg2_arg, const ExpressionPtr arg3_arg,
              Environment &env_arg, const Tokenizer::location location_arg) :
      Expression(env_arg, move(location_arg)), op_code{move(op_code_arg)},
      arg1{move(arg1_arg)}, arg2{move(arg2_arg)}, arg3{move(arg3_arg)} { }
    string to_string() const noexcept override;
    void print(ostream &output, bool matlab_output = false) const noexcept override;
    BaseTypePtr eval() override;
  };


  class Comprehension final : public Expression
  {
  private:
    const ExpressionPtr c_vars, c_set, c_when;
  public:
    Comprehension(const ExpressionPtr c_vars_arg,
                  const ExpressionPtr c_set_arg,
                  const ExpressionPtr c_when_arg,
                  Environment &env_arg, const Tokenizer::location location_arg) :
      Expression(env_arg, move(location_arg)),
      c_vars{move(c_vars_arg)}, c_set{move(c_set_arg)}, c_when{move(c_when_arg)} { }
    string to_string() const noexcept override;
    void print(ostream &output, bool matlab_output = false) const noexcept override;
    BaseTypePtr eval() override;
  };
}
#endif
