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

#ifndef _STATICDLLMODEL_HH
#define _STATICDLLMODEL_HH

using namespace std;

#include <fstream>

#include "StaticModel.hh"
#include "BlockTriangular.hh"

//! Stores a static model
class StaticDllModel : public ModelTree
{
private:
  typedef map<pair<int, int>, int> deriv_id_table_t;
  //! Maps a pair (symbol_id, lag) to a deriv ID
  deriv_id_table_t deriv_id_table;
  //! Maps a deriv ID to a pair (symbol_id, lag)
  vector<pair<int, int> > inv_deriv_id_table;

  //! Temporary terms for the file containing parameters dervicatives
  temporary_terms_type params_derivs_temporary_terms;

  typedef map< pair< int, pair< int, int> >, NodeID> first_chain_rule_derivatives_type;
  first_chain_rule_derivatives_type first_chain_rule_derivatives;


  //! Writes static model file (Matlab version)
  void writeStaticMFile(const string &static_basename) const;
  //! Writes static model file (C version)
  /*! \todo add third derivatives handling */
  void writeStaticCFile(const string &static_basename) const;
  //! Writes the Block reordred structure of the model in M output
  void writeModelEquationsOrdered_M(Model_Block *ModelBlock, const string &static_basename) const;
  //! Writes the code of the Block reordred structure of the model in virtual machine bytecode
  void writeModelEquationsCodeOrdered(const string file_name, const Model_Block *ModelBlock, const string bin_basename, map_idx_type map_idx) const;
  //! Computes jacobian and prepares for equation normalization
  /*! Using values from initval/endval blocks and parameter initializations:
    - computes the jacobian for the model w.r. to contemporaneous variables
    - removes edges of the incidence matrix when derivative w.r. to the corresponding variable is too close to zero (below the cutoff)
  */
  void evaluateJacobian(const eval_context_type &eval_context, jacob_map *j_m, bool dynamic);
  void BlockLinear(Model_Block *ModelBlock);
  map_idx_type map_idx;
  //! Build The incidence matrix form the modeltree
  void BuildIncidenceMatrix();

  void computeTemporaryTermsOrdered(Model_Block *ModelBlock);
  //! Write derivative code of an equation w.r. to a variable
  void compileDerivative(ofstream &code_file, int eq, int symb_id, int lag, map_idx_type &map_idx) const;
  //! Write chain rule derivative code of an equation w.r. to a variable
  void compileChainRuleDerivative(ofstream &code_file, int eq, int var, int lag, map_idx_type &map_idx) const;

  //! Get the type corresponding to a derivation ID
  SymbolType getTypeByDerivID(int deriv_id) const throw (UnknownDerivIDException);
  //! Get the lag corresponding to a derivation ID
  int getLagByDerivID(int deriv_id) const throw (UnknownDerivIDException);
  //! Get the symbol ID corresponding to a derivation ID
  int getSymbIDByDerivID(int deriv_id) const throw (UnknownDerivIDException);
  //! Compute the column indices of the static Jacobian
  void computeStatJacobianCols();
  //! Computes chain rule derivatives of the Jacobian w.r. to endogenous variables
  void computeChainRuleJacobian(Model_Block *ModelBlock);
  //! Collect only the first derivatives
  map<pair<int, pair<int, int> >, NodeID> collect_first_order_derivatives_endogenous();

  //! Helper for writing the Jacobian elements in MATLAB and C
  /*! Writes either (i+1,j+1) or [i+j*no_eq] */
  void jacobianHelper(ostream &output, int eq_nb, int col_nb, ExprNodeOutputType output_type) const;

  //! Helper for writing the sparse Hessian elements in MATLAB and C
  /*! Writes either (i+1,j+1) or [i+j*NNZDerivatives[1]] */
  void hessianHelper(ostream &output, int row_nb, int col_nb, ExprNodeOutputType output_type) const;

  //! Write chain rule derivative of a recursive equation w.r. to a variable
  void writeChainRuleDerivative(ostream &output, int eq, int var, int lag, ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms) const;


public:
  StaticDllModel(SymbolTable &symbol_table_arg, NumericalConstants &num_constants);
  //! Absolute value under which a number is considered to be zero
  double cutoff;
  //! Compute the minimum feedback set in the static model:
  /*!   0 : all endogenous variables are considered as feedback variables
				1 : the variables belonging to a non linear equation are considered as feedback variables
        2 : the variables belonging to a non normalizable non linear equation are considered as feedback variables
        default value = 0 */
  int mfs;
  //! the file containing the model and the derivatives code
  ofstream code_file;
  //! Execute computations (variable sorting + derivation)
  /*!
    \param jacobianExo whether derivatives w.r. to exo and exo_det should be in the Jacobian (derivatives w.r. to endo are always computed)
    \param hessian whether 2nd derivatives w.r. to exo, exo_det and endo should be computed (implies jacobianExo = true)
    \param thirdDerivatives whether 3rd derivatives w.r. to endo/exo/exo_det should be computed (implies jacobianExo = true)
    \param paramsDerivatives whether 2nd derivatives w.r. to a pair (endo/exo/exo_det, parameter) should be computed (implies jacobianExo = true)
    \param eval_context evaluation context for normalization
    \param no_tmp_terms if true, no temporary terms will be computed in the static files
  */
  void computingPass(const eval_context_type &eval_context, bool no_tmp_terms, bool block);
  //! Complete set to block decompose the model
  BlockTriangular block_triangular;
  //! Adds informations for simulation in a binary file
  void Write_Inf_To_Bin_File(const string &static_basename, const string &bin_basename,
                             const int &num, int &u_count_int, bool &file_open) const;
  //! Writes static model file
  void writeStaticFile(const string &basename, bool block) const;

  //! Writes LaTeX file with the equations of the static model
  void writeLatexFile(const string &basename) const;

  //! Writes initializations in oo_.steady_state for the auxiliary variables
  void writeAuxVarInitval(ostream &output) const;

  virtual int getDerivID(int symb_id, int lag) const throw (UnknownDerivIDException);
};

#endif
