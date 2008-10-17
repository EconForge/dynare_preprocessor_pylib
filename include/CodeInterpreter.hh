/*
 * Copyright (C) 2007-2008 Dynare Team
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

#ifndef _CODEINTERPRETER_HH
#define _CODEINTERPRETER_HH

const char FLDZ=0;
const char FLDT=1;
const char FLDU=2;
const char FLDV=3;
const char FLDR=4;
const char FLDC=5;
const char FSTPT=6;
const char FSTPV=7;
const char FSTPR=8;
const char FSTPU=9;
const char FSTPG=10;
const char FUNARY=11;
const char FBINARY=12;
const char FCUML=13;
const char FBEGINBLOCK=14;
const char FENDBLOCK=15;
const char FDIMT=16;
const char FEND=17;
const char FOK=18;
const char FENDEQU=19;

enum BlockType
  {
    SIMULTANS = 0, //<! Simultaneous time separable block
    PROLOGUE = 1,  //<! Prologue block (one equation at the beginning, later merged)
    EPILOGUE = 2,  //<! Epilogue block (one equation at the beginning, later merged)
    SIMULTAN = 3   //<! Simultaneous time unseparable block
  };

enum BlockSimulationType
  {
    UNKNOWN = -1,                      //!< Unknown simulation type
    EVALUATE_FOREWARD = 0,             //!< Simple evaluation, normalized variable on left-hand side, forward
    EVALUATE_BACKWARD = 1,             //!< Simple evaluation, normalized variable on left-hand side, backward
    SOLVE_FOREWARD_SIMPLE = 2,         //!< Block of one equation, newton solver needed, forward
    SOLVE_BACKWARD_SIMPLE = 3,         //!< Block of one equation, newton solver needed, backward
    SOLVE_TWO_BOUNDARIES_SIMPLE = 4,   //!< Block of one equation, newton solver needed, forward & ackward
    SOLVE_FOREWARD_COMPLETE = 5,       //!< Block of several equations, newton solver needed, forward
    SOLVE_BACKWARD_COMPLETE = 6,       //!< Block of several equations, newton solver needed, backward
    SOLVE_TWO_BOUNDARIES_COMPLETE = 7, //!< Block of several equations, newton solver needed, forward and backwar
    EVALUATE_FOREWARD_R = 8,           //!< Simple evaluation, normalized variable on right-hand side, forward
    EVALUATE_BACKWARD_R = 9            //!< Simple evaluation, normalized variable on right-hand side, backward
  };

//! Enumeration of possible symbol types
/*! Warning: do not to change existing values: the order matters for VariableTable (at least for endogenous and exogenous types), and the values matter for homotopy_setup command */
enum SymbolType
  {
    eEndogenous = 0,               //!< Endogenous
    eExogenous = 1,                //!< Exogenous
    eExogenousDet = 2,             //!< Exogenous deterministic
    eParameter = 4,                //!< Parameter
    eModelLocalVariable = 10,      //!< Local variable whose scope is model (pound expression)
    eModFileLocalVariable = 11,    //!< Local variable whose scope is mod file (model excluded)
    eUnknownFunction = 12          //!< Function unknown to the preprocessor
  };

enum UnaryOpcode
  {
    oUminus,
    oExp,
    oLog,
    oLog10,
    oCos,
    oSin,
    oTan,
    oAcos,
    oAsin,
    oAtan,
    oCosh,
    oSinh,
    oTanh,
    oAcosh,
    oAsinh,
    oAtanh,
    oSqrt
  };

enum BinaryOpcode
  {
    oPlus,
    oMinus,
    oTimes,
    oDivide,
    oPower,
    oEqual,
    oMax,
    oMin,
    oLess,
    oGreater,
    oLessEqual,
    oGreaterEqual,
    oEqualEqual,
    oDifferent
  };

#endif
