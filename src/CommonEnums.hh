/*
 * Copyright © 2007-2022 Dynare Team
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
 * along with Dynare.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef _COMMON_ENUMS_HH
#define _COMMON_ENUMS_HH

//! Enumeration of possible symbol types
/*! Warning: do not to change existing values for 0 to 4: the values matter for homotopy_setup command */
enum class SymbolType
  {
   endogenous = 0, //!< Endogenous
   exogenous = 1, //!< Exogenous
   exogenousDet = 2, //!< Exogenous deterministic
   parameter = 4, //!< Parameter
   modelLocalVariable = 10, //!< Local variable whose scope is model (pound expression)
   modFileLocalVariable = 11, //!< Local variable whose scope is mod file (model excluded)
   externalFunction = 12, //!< External (user-defined) function
   trend = 13, //!< Trend variable
   statementDeclaredVariable = 14, //!< Local variable assigned within a Statement (see subsample statement for example)
   logTrend = 15, //!< Log-trend variable
   unusedEndogenous = 16, //!< Type to mark unused endogenous variables when `nostrict` option is passed

   // Value 17 is unused for the time being (but could be reused)

   epilogue = 18, //!< Variables created in epilogue block
   excludedVariable = 19 //!< Variable excluded via model_remove/var_remove/include_eqs/exclude_eqs
  };

enum class UnaryOpcode
  {
   uminus,
   exp,
   log,
   log10,
   cos,
   sin,
   tan,
   acos,
   asin,
   atan,
   cosh,
   sinh,
   tanh,
   acosh,
   asinh,
   atanh,
   sqrt,
   cbrt,
   abs,
   sign,
   steadyState,
   steadyStateParamDeriv, // for the derivative of the STEADY_STATE operator w.r.t. to a parameter
   steadyStateParam2ndDeriv, // for the 2nd derivative of the STEADY_STATE operator w.r.t. to a parameter
   expectation,
   erf,
   erfc,
   diff,
   adl
  };

enum class BinaryOpcode
  {
   plus,
   minus,
   times,
   divide,
   power,
   powerDeriv, // for the derivative of the power function (see trac ticket #78)
   equal,
   max,
   min,
   less,
   greater,
   lessEqual,
   greaterEqual,
   equalEqual,
   different
  };

// Small number value used when evaluating powerDeriv opcodes.
// Put here instead of inside BinaryOpNode class, because needed by bytecode MEX.
constexpr double power_deriv_near_zero{1e-12};

enum class TrinaryOpcode
  {
   normcdf,
   normpdf
  };

enum class ExternalFunctionType
  {
   withoutDerivative,
   withFirstDerivative,
   withFirstAndSecondDerivative,
   numericalFirstDerivative,
   firstDerivative,
   numericalSecondDerivative,
   secondDerivative
  };

enum class PriorDistributions
  {
   noShape = 0,
   beta = 1,
   gamma = 2,
   normal = 3,
   invGamma = 4,
   invGamma1 = 4,
   uniform = 5,
   invGamma2 = 6,
   dirichlet = 7,
   weibull = 8
  };

enum class EquationType
  {
   unknown, //!< Unknown equation type
   evaluate, //!< Simple evaluation, normalized variable on left-hand side (written as such by the user)
   evaluateRenormalized, //!< Simple evaluation, normalized variable on left-hand side (normalization computed by the preprocessor)
   solve //!< No simple evaluation of the equation, it has to be solved
  };

enum class BlockSimulationType
  {
   unknown, //!< Unknown simulation type
   evaluateForward, //!< Simple evaluation, normalized variable on left-hand side, forward
   evaluateBackward, //!< Simple evaluation, normalized variable on left-hand side, backward
   solveForwardSimple, //!< Block of one equation, newton solver needed, forward
   solveBackwardSimple, //!< Block of one equation, newton solver needed, backward
   solveTwoBoundariesSimple, //!< Block of one equation, Newton solver needed, forward and backward
   solveForwardComplete, //!< Block of several equations, Newton solver needed, forward
   solveBackwardComplete, //!< Block of several equations, Newton solver needed, backward
   solveTwoBoundariesComplete //!< Block of several equations, Newton solver needed, forward and backwar
  };

enum class PacTargetKind
  {
    unspecified, // Must be the first one, because it’s the default initializer
    ll,
    dl,
    dd
  };

#endif // _COMMON_ENUMS_HH
