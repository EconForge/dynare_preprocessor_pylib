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

#ifndef _FORWARDDECLARATIONSANDENUMS_HH
#define _FORWARDDECLARATIONSANDENUMS_HH

#include <memory>

using namespace std;

namespace macro
{
  // For Expressions.hh
  class BaseType;
  using BaseTypePtr = shared_ptr<BaseType>;
  class Bool;
  using BoolPtr = shared_ptr<Bool>;
  class Real;
  using RealPtr = shared_ptr<Real>;
  class String;
  using StringPtr = shared_ptr<String>;
  class Tuple;
  using TuplePtr = shared_ptr<Tuple>;
  class Array;
  using ArrayPtr = shared_ptr<Array>;

  // For Environment.hh
  class Expression;
  using ExpressionPtr = shared_ptr<Expression>;
  class Variable;
  using VariablePtr = shared_ptr<Variable>;
  class Function;
  using FunctionPtr = shared_ptr<Function>;

  // For Parser.yy
  class Directive;
  using DirectivePtr = shared_ptr<Directive>;
  class Eval;
  using EvalPtr = shared_ptr<Eval>;

  namespace codes
  {
    enum class BaseType
      {
       Bool,
       Real,
       String,
       Array,
       Tuple
      };

    enum class UnaryOp
      {
       cast_bool,
       cast_real,
       cast_string,
       cast_tuple,
       cast_array,
       logical_not,
       unary_minus,
       unary_plus,
       length,
       isempty,
       isboolean,
       isreal,
       isstring,
       istuple,
       isarray,
       exp,
       ln,
       log10,
       sin,
       cos,
       tan,
       asin,
       acos,
       atan,
       sqrt,
       cbrt,
       sign,
       floor,
       ceil,
       trunc,
       sum,
       erf,
       erfc,
       gamma,
       lgamma,
       round,
       normpdf,
       normcdf
      };

    enum class BinaryOp
      {
       plus,
       minus,
       times,
       divide,
       power,
       equal_equal,
       not_equal,
       less,
       greater,
       less_equal,
       greater_equal,
       logical_and,
       logical_or,
       in,
       set_union,
       set_intersection,
       max,
       min,
       mod
      };

    enum class TrinaryOp
      {
       normpdf,
       normcdf
      };
  }
}

#endif
