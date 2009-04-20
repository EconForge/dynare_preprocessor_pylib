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

#ifndef _DYNAMICMODEL_HH
#define _DYNAMICMODEL_HH

using namespace std;

#include "StaticModel.hh"
#include "BlockTriangular.hh"

//! Stores a dynamic model
class DynamicModel : public ModelTree
{
private:
  typedef map<pair<int, int>, int> deriv_id_table_t;
  //! Maps a pair (symbol_id, lag) to a deriv ID
  deriv_id_table_t deriv_id_table;
  //! Maps a deriv ID to a pair (symbol_id, lag)
  vector<pair<int, int> > inv_deriv_id_table;

  //! Maps a deriv_id to the column index of the dynamic Jacobian
  /*! Contains only endogenous, exogenous and exogenous deterministic */
  map<int, int> dyn_jacobian_cols_table;

  //! Maximum lag and lead over all types of variables (positive values)
  /*! Set by computeDerivID() */
  int max_lag, max_lead;
  //! Maximum lag and lead over endogenous variables (positive values)
  /*! Set by computeDerivID() */
  int max_endo_lag, max_endo_lead;
  //! Maximum lag and lead over exogenous variables (positive values)
  /*! Set by computeDerivID() */
  int max_exo_lag, max_exo_lead;
  //! Maximum lag and lead over deterministic exogenous variables (positive values)
  /*! Set by computeDerivID() */
  int max_exo_det_lag, max_exo_det_lead;

  //! Number of columns of dynamic jacobian
  /*! Set by computeDerivID() and computeDynJacobianCols() */
  int dynJacobianColsNbr;

  //! Derivatives of the jacobian w.r. to parameters
  /*! First index is equation number, second is endo/exo/exo_det variable, and third is parameter.
    Only non-null derivatives are stored in the map.
    Variable and parameter indices are those of the getDerivID() method.
  */
  second_derivatives_type params_derivatives;

  //! Temporary terms for the file containing parameters dervicatives
  temporary_terms_type params_derivs_temporary_terms;

  //! Writes dynamic model file (Matlab version)
  void writeDynamicMFile(const string &dynamic_basename) const;
  //! Writes dynamic model file (C version)
  /*! \todo add third derivatives handling */
  void writeDynamicCFile(const string &dynamic_basename) const;
  //! Writes dynamic model file when SparseDLL option is on
  void writeSparseDynamicMFile(const string &dynamic_basename, const string &basename, const int mode) const;
  //! Writes the dynamic model equations and its derivatives
  /*! \todo add third derivatives handling in C output */
  void writeDynamicModel(ostream &DynamicOutput) const;
  //! Writes the Block reordred structure of the model in M output
  void writeModelEquationsOrdered_M(Model_Block *ModelBlock, const string &dynamic_basename) const;
  //! Writes the code of the Block reordred structure of the model in virtual machine bytecode
  void writeModelEquationsCodeOrdered(const string file_name, const Model_Block *ModelBlock, const string bin_basename, ExprNodeOutputType output_type, map_idx_type map_idx) const;
  //! Computes jacobian and prepares for equation normalization
  /*! Using values from initval/endval blocks and parameter initializations:
    - computes the jacobian for the model w.r. to contemporaneous variables
    - removes edges of the incidence matrix when derivative w.r. to the corresponding variable is too close to zero (below the cutoff)
  */
  void evaluateJacobian(const eval_context_type &eval_context, jacob_map *j_m);
  void BlockLinear(Model_Block *ModelBlock);
  string reform(string name) const;
  map_idx_type map_idx;
  //! Build The incidence matrix form the modeltree
  void BuildIncidenceMatrix();

  void computeTemporaryTermsOrdered(Model_Block *ModelBlock);
  //! Write derivative code of an equation w.r. to a variable
  void compileDerivative(ofstream &code_file, int eq, int symb_id, int lag, ExprNodeOutputType output_type, map_idx_type &map_idx) const;

  virtual int computeDerivID(int symb_id, int lag);
  //! Get the type corresponding to a derivation ID
  SymbolType getTypeByDerivID(int deriv_id) const throw (UnknownDerivIDException);
  //! Get the lag corresponding to a derivation ID
  int getLagByDerivID(int deriv_id) const throw (UnknownDerivIDException);
  //! Get the symbol ID corresponding to a derivation ID
  int getSymbIDByDerivID(int deriv_id) const throw (UnknownDerivIDException);
  //! Compute the column indices of the dynamic Jacobian
  void computeDynJacobianCols(bool jacobianExo);
  //! Computes derivatives of the Jacobian w.r. to parameters
  void computeParamsDerivatives();
  //! Computes temporary terms for the file containing parameters derivatives
  void computeParamsDerivativesTemporaryTerms();

public:
  DynamicModel(SymbolTable &symbol_table_arg, NumericalConstants &num_constants);
  //! Adds a variable node
  /*! This implementation allows for non-zero lag */
  virtual NodeID AddVariable(const string &name, int lag = 0);
  //! Absolute value under which a number is considered to be zero
  double cutoff;
  //! The weight of the Markowitz criteria to determine the pivot in the linear solver (simul_NG1 from simulate.cc)
  double markowitz;
  //! the file containing the model and the derivatives code
  ofstream code_file;
  //! Execute computations (variable sorting + derivation)
  /*!
    \param jacobianExo whether derivatives w.r. to exo and exo_det should be in the Jacobian (derivatives w.r. to endo are always computed)
    \param hessian whether 2nd derivatives w.r. to exo, exo_det and endo should be computed (implies jacobianExo = true)
    \param thirdDerivatives whether 3rd derivatives w.r. to endo/exo/exo_det should be computed (implies jacobianExo = true)
    \param paramsDerivatives whether 2nd derivatives w.r. to a pair (endo/exo/exo_det, parameter) should be computed (implies jacobianExo = true)
    \param eval_context evaluation context for normalization
    \param no_tmp_terms if true, no temporary terms will be computed in the dynamic files
  */
  void computingPass(bool jacobianExo, bool hessian, bool thirdDerivatives, bool paramsDerivatives,
                     const eval_context_type &eval_context, bool no_tmp_terms);
  //! Writes model initialization and lead/lag incidence matrix to output
  void writeOutput(ostream &output) const;
  //! Complete set to block decompose the model
  BlockTriangular block_triangular;
  //! Adds informations for simulation in a binary file
  void Write_Inf_To_Bin_File(const string &dynamic_basename, const string &bin_basename,
                             const int &num, int &u_count_int, bool &file_open, bool is_two_boundaries) const;
  //! Writes dynamic model file
  void writeDynamicFile(const string &basename) const;
  //! Writes file containing parameters derivatives
  void writeParamsDerivativesFile(const string &basename) const;
  //! Converts to static model (only the equations)
  /*! It assumes that the static model given in argument has just been allocated */
  void toStatic(StaticModel &static_model) const;
  virtual int getDerivID(int symb_id, int lag) const throw (UnknownDerivIDException);
  virtual int getDynJacobianCol(int deriv_id) const throw (UnknownDerivIDException);
};

#endif
