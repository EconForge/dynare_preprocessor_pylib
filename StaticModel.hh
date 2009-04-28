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
  //! Writes the static model equations and its derivatives
  /*! \todo handle hessian in C output */
  void writeStaticModel(ostream &StaticOutput) const;
  //! Writes static model file (Matlab version)
  void writeStaticMFile(const string &static_basename) const;
  //! Writes static model file (C version)
  void writeStaticCFile(const string &static_basename) const;

  virtual int computeDerivID(int symb_id, int lag);

  //! Computes normalization of the static model
  /*! Maps each endogenous type specific ID to the equation which defines it */
  void computeNormalization(vector<int> &endo_to_eq) const;

  //! Computes the list of equations which are already in normalized form
  /*! Returns a multimap mapping endogenous which are normalized (represented by their type specific ID) to the equation(s) which define it */
  void computeNormalizedEquations(multimap<int, int> &endo_to_eqs) const;

public:
  StaticModel(SymbolTable &symbol_table_arg, NumericalConstants &num_constants);
  //! Execute computations (derivation)
  /*! You must set computeStaticHessian before calling this function
    \param hessian whether Hessian (w.r. to endogenous only) should be computed
    \param no_tmp_terms if true, no temporary terms will be computed in the static and dynamic files */
  void computingPass(bool hessian, bool no_tmp_terms);
  //! Writes static model file
  void writeStaticFile(const string &basename) const;

  virtual int getDerivID(int symb_id, int lag) const throw (UnknownDerivIDException);
};

#endif
