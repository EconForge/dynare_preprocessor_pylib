/*
 * Copyright (C) 2003-2009 Dynare Team
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

#ifndef _STATICMODEL_HH
#define _STATICMODEL_HH

#include "ModelTree.hh"

//! Stores a static model
/*! Derivation IDs are allocated only for endogenous, and are equal to symbol ID in that case */
class StaticModel : public ModelTree
{
private:
  //! Are we in block decomposition + min. feedback set mode ?
  bool block_mfs;

  //! Normalization of equations
  /*! Maps endogenous type specific IDs to equation numbers */
  vector<int> endo2eq;

  //! Block decomposition of the model
  /*! List of blocks in topological order. Lists the set of endogenous type specific IDs for each block. */
  vector<set<int> > blocks;

  //! Minimum feedback set for each block
  /*! Elements of blocksMFS are subset of elements of blocks */
  vector<set<int> > blocksMFS;

  //! Variables not in minimum feedback set for each block, sorted in topological order
  /*! This is the set difference blocks - blocksMFS. The variables are sorted in topological order. */
  vector<vector<int> > blocksRecursive;

  //! Jacobian for matrix restricted to MFS
  /*! Maps a pair (equation ID, endogenous type specific ID) to the derivative expression. Stores only non-null derivatives. */
  map<pair<int, int>, NodeID> blocksMFSJacobian;

  //! Writes static model file (standard Matlab version)
  void writeStaticMFile(ostream &output, const string &func_name) const;

  //! Writes static model file (block+MFS version)
  void writeStaticBlockMFSFile(ostream &output, const string &func_name) const;

  virtual int computeDerivID(int symb_id, int lag);

  //! Computes normalization of the static model
  void computeNormalization();

  //! Computes blocks of the static model, sorted in topological order
  /*! Must be called after computeNormalization() */
  void computeSortedBlockDecomposition();

  //! For each block of the static model, computes minimum feedback set (MFS)
  /*! Must be called after computeSortedBlockDecomposition() */
  void computeMFS();

  //! For each block of the static model, computes resursive variables (those not in minimum feedback set), and sort them in topological order
  /*! Must be called after computeMFS() */
  void computeSortedRecursive();

  //! Computes derivatives of each MFS
  /*! Must be called after computeSortedRecursive() */
  void computeBlockMFSJacobian();

  //! Computes the list of equations which are already in normalized form
  /*! Returns a multimap mapping endogenous which are normalized (represented by their type specific ID) to the equation(s) which define it */
  void computeNormalizedEquations(multimap<int, int> &endo2eqs) const;

  //! Helper for writing model local variables in block+MFS mode
  /*!
    Write the definition of model local variables which are used in expr, except those in local_var_written.
    Add these variables to local_var_written at the end.
  */
  void writeLocalVars(ostream &output, NodeID expr, set<int> &local_var_written) const;

public:
  StaticModel(SymbolTable &symbol_table_arg, NumericalConstants &num_constants);
  //! Execute computations (derivation)
  /*!
    \param block_mfs whether block decomposition and minimum feedback set should be computed
    \param hessian whether Hessian (w.r. to endogenous only) should be computed
    \param no_tmp_terms if true, no temporary terms will be computed in the static and dynamic files */
  void computingPass(bool block_mfs_arg, bool hessian, bool no_tmp_terms);

  //! Writes information on block decomposition when relevant
  void writeOutput(ostream &output) const;

  //! Writes static model file
  void writeStaticFile(const string &basename) const;

  //! Writes LaTeX file with the equations of the static model
  void writeLatexFile(const string &basename) const;

  virtual int getDerivID(int symb_id, int lag) const throw (UnknownDerivIDException);
};

#endif
