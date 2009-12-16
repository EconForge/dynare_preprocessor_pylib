/*
 * Copyright (C) 2006-2009 Dynare Team
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

#ifndef _MOD_FILE_HH
#define _MOD_FILE_HH

using namespace std;

#include <ostream>
#include <ctime>

#include "SymbolTable.hh"
#include "NumericalConstants.hh"
#include "NumericalInitialization.hh"
#include "StaticModel.hh"
#include "DynamicModel.hh"
#include "Statement.hh"

//! The abstract representation of a "mod" file
class ModFile
{
public:
  ModFile();
  ~ModFile();
  //! Symbol table
  SymbolTable symbol_table;
  //! Numerical constants table
  NumericalConstants num_constants;
  //! Expressions outside model block
  DataTree expressions_tree;
  //! Static Dll model
  StaticModel static_model;
  //! Dynamic model
  DynamicModel dynamic_model;
  //! Option linear
  bool linear;

  //! Is the model block decomposed?
  bool block;

  //! Is the model stored in bytecode format (byte_code=true) or in a M-file (byte_code=false)
  bool byte_code;

  //! Is the model stored in a MEX file ? (option "use_dll" of "model")
  bool use_dll;

  //! Global evaluation context
  /*! Filled using initval blocks and parameters initializations */
  eval_context_type global_eval_context;

private:
  //! List of statements
  vector<Statement *> statements;
  //! Structure of the mod file
  ModFileStructure mod_file_struct;

public:
  //! Add a statement
  void addStatement(Statement *st);
  //! Evaluate all the statements
  /*! \param warn_uninit Should a warning be displayed for uninitialized endogenous/exogenous/parameters ? */
  void evalAllExpressions(bool warn_uninit);
  //! Do some checking and fills mod_file_struct
  /*! \todo add check for number of equations and endogenous if ramsey_policy is present */
  void checkPass();
  //! Perform some transformations on the model (creation of auxiliary vars and equations)
  void transformPass();
  //! Execute computations
  /*! \param no_tmp_terms if true, no temporary terms will be computed in the static and dynamic files */
  void computingPass(bool no_tmp_terms);
  //! Writes Matlab/Octave output files
  /*!
    \param basename The base name used for writing output files. Should be the name of the mod file without its extension
    \param clear_all Should a "clear all" instruction be written to output ?
    \param msvc Should the MEX command of use_dll be adapted for MSVC?
  */
  void writeOutputFiles(const string &basename, bool clear_all
#if defined(_WIN32) || defined(__CYGWIN32__)
                        , bool cygwin, bool msvc
#endif
                        ) const;
};

#endif // ! MOD_FILE_HH
