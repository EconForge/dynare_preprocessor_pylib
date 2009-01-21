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
#include "ModelTree.hh"
#include "VariableTable.hh"
#include "Statement.hh"
#include "MatlabFile.hh"

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
  //! Model equations and their derivatives
  ModelTree model_tree;
  //! MatFile reading
  MatlabFile matlab_file;
  //! Option linear
  bool linear;
  //! Global evaluation context
  /*! Filled using initval blocks and parameters initializations */
  eval_context_type global_eval_context;
  //! Temporary storage for initval/endval blocks
  InitOrEndValStatement::init_values_type init_values;

private:
  //! List of statements
  vector<Statement *> statements;
  //! Structure of the mod file
  ModFileStructure mod_file_struct;

public:
  //! Add a statement
  void addStatement(Statement *st);
  //! Evaluate all the statements
  void evalAllExpressions();
  //! Do some checking and fills mod_file_struct
  /*! \todo add check for number of equations and endogenous if ramsey_policy is present */
  void checkPass();
  //! Execute computations
  /*! \param no_tmp_terms if true, no temporary terms will be computed in the static and dynamic files */
  void computingPass(bool no_tmp_terms);
  //! Writes Matlab/Scilab output files
  /*!
    \param basename The base name used for writing output files. Should be the name of the mod file without its extension
    \param clear_all Should a "clear all" instruction be written to output ?
  */
  void writeOutputFiles(const string &basename, bool clear_all) const;
};

#endif // ! MOD_FILE_HH
