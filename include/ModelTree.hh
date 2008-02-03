/*
 * Copyright (C) 2003-2008 Dynare Team
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

#ifndef _MODELTREE_HH
#define _MODELTREE_HH

using namespace std;

#include <string>
#include <vector>
#include <map>
#include <ostream>

#include "SymbolTable.hh"
#include "NumericalConstants.hh"
#include "DataTree.hh"
#include "BlockTriangular.hh"
#include "SymbolGaussElim.hh"

#define LCC_COMPILE 0
#define GCC_COMPILE 1
#define NO_COMPILE 2
//#define CONDITION

//! The three in which ModelTree can work
enum ModelTreeMode
  {
    eStandardMode, //!< Standard mode (static and dynamic files in Matlab)
    eSparseMode,  //!< Sparse mode (static file in Matlab, dynamic file in Matlab with block decomposition)
    eDLLMode,      //!< DLL mode (static and dynamic files in C)
    eSparseDLLMode //!< Sparse DLL mode (static file in Matlab, dynamic file in C with block decomposition plus a binary file)
  };

//! Stores a model's equations and derivatives
class ModelTree : public DataTree
{
private:
  //! Stores declared equations
  vector<BinaryOpNode *> equations;

  typedef map<pair<int, int>, NodeID> first_derivatives_type;
  //! First order derivatives
  /*! First index is equation number, second is variable w.r. to which is computed the derivative.
    Only non-null derivatives are stored in the map.
    Variable indexes used are those of the variable_table, before sorting.
  */
  first_derivatives_type first_derivatives;

  typedef map<pair<int, pair<int, int> >, NodeID> second_derivatives_type;
  //! Second order derivatives
  /*! First index is equation number, second and third are variables w.r. to which is computed the derivative.
    Only non-null derivatives are stored in the map.
    Contains only second order derivatives where var1 >= var2 (for obvious symmetry reasons).
    Variable indexes used are those of the variable_table, before sorting.
  */
  second_derivatives_type second_derivatives;

  typedef map<pair<int, pair<int, pair<int, int> > >, NodeID> third_derivatives_type;
  //! Third order derivatives
  /*! First index is equation number, second, third and fourth are variables w.r. to which is computed the derivative.
    Only non-null derivatives are stored in the map.
    Contains only third order derivatives where var1 >= var2 >= var3 (for obvious symmetry reasons).
    Variable indexes used are those of the variable_table, before sorting.
  */
  third_derivatives_type third_derivatives;

  //! Temporary terms (those which will be noted Txxxx)
  temporary_terms_type temporary_terms;
  map_idx_type map_idx;

  //! Computes derivatives of ModelTree
  void derive(int order);
  //! Write derivative of an equation w.r. to a variable
  void writeDerivative(ostream &output, int eq, int symb_id, int lag, ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms) const;
  //! Write derivative code of an equation w.r. to a variable
  void compileDerivative(ofstream &code_file, int eq, int symb_id, int lag, ExprNodeOutputType output_type, map_idx_type map_idx) const;
  //! Computes temporary terms
  void computeTemporaryTerms(int order);
  void computeTemporaryTermsOrdered(int order, Model_Block *ModelBlock);
  //! Writes temporary terms
  void writeTemporaryTerms(ostream &output, ExprNodeOutputType output_type) const;
  //! Writes model local variables
  /*! No temporary term is used in the output, so that local parameters declarations can be safely put before temporary terms declaration in the output files */
  void writeModelLocalVariables(ostream &output, ExprNodeOutputType output_type) const;
  //! Writes model equations
  void writeModelEquations(ostream &output, ExprNodeOutputType output_type) const;
  //! Writes the static model equations and its derivatives
  /*! \todo handle hessian in C output */
  void writeStaticModel(ostream &StaticOutput) const;
  //! Writes the dynamic model equations and its derivatives
  /*! \todo add third derivatives handling in C output */
  void writeDynamicModel(ostream &DynamicOutput) const;
  //! Writes the Block reordred structure of the model in C output
  void writeModelEquationsOrdered_C(ostream &output, Model_Block *ModelBlock) const;
  //! Writes the Block reordred structure of the model in M output
  void writeModelEquationsOrdered_M(ostream &output, Model_Block *ModelBlock, const string &dynamic_basename) const;
  //! Writes the Block reordred structure of the static model in M output
  void writeModelStaticEquationsOrdered_M(ostream &output, Model_Block *ModelBlock, const string &static_basename) const;
  //! Writes the code of the Block reordred structure of the model in C output
  void writeModelEquationsCodeOrdered(const string file_name, const Model_Block *ModelBlock, const string bin_basename, ExprNodeOutputType output_type) const;
  //! Writes static model file (Matlab version)
  void writeStaticMFile(const string &static_basename) const;
  //! Writes static model file (C version)
  void writeStaticCFile(const string &static_basename) const;
  //! Writes static model file when Sparse option is on (Matlab version)
  void writeSparseStaticMFile(const string &static_basename, const string &bin_basename, const int mode) const;
  //! Writes dynamic model file (Matlab version)
  void writeDynamicMFile(const string &dynamic_basename) const;
  //! Writes dynamic model file (C version)
  /*! \todo add third derivatives handling */
  void writeDynamicCFile(const string &dynamic_basename) const;
  //! Writes dynamic model header file when SparseDLL option is on
  void writeSparseDLLDynamicHFile(const string &dynamic_basename) const;
  //! Writes dynamic model file when SparseDLL option is on
  void writeSparseDynamicFileAndBinFile(const string &dynamic_basename, const string &bin_basename, ExprNodeOutputType output_type, const int mode) const;
  //! Computes jacobian and prepares for equation normalization
  /*! Using values from initval/endval blocks and parameter initializations:
    - computes the jacobian for the model w.r. to contemporaneous variables
    - removes edges of the incidence matrix when derivative w.r. to the corresponding variable is too close to zero (below the cutoff)
  */
  void evaluateJacobian(const eval_context_type &eval_context, jacob_map *j_m);
  void BlockLinear(Model_Block *ModelBlock);
  string reform(string name) const;

  //! Writes either (i+1,j+1) or [i+j*n_i] whether we are in Matlab or C mode
  void matrixHelper(ostream &output, int eq_nb, int col_nb, ExprNodeOutputType output_type) const;

public:
  ModelTree(SymbolTable &symbol_table_arg, NumericalConstants &num_constants);
  //! Mode in which the ModelTree is supposed to work (Matlab, DLL or SparseDLL)
  ModelTreeMode mode;
  //! Type of compiler used in matlab for SPARSE_DLL option: 0 = LCC or 1 = GCC or 2 = NO
  int compiler;
  //! Absolute value under which a number is considered to be zero
  double cutoff;
  //! The weight of the Markowitz criteria to determine the pivot in the linear solver (simul_NG1 from simulate.cc)
  double markowitz;
  //! Use a graphical and symbolic version of the symbolic gaussian elimination (new_SGE = false) or use direct gaussian elimination (new_SGE = true)
  bool new_SGE;
  //! the file containing the model and the derivatives code
  ofstream code_file;
  //! Declare a node as an equation of the model
  void addEquation(NodeID eq);
  //! Whether dynamic Jacobian (w.r. to endogenous) should be written
  bool computeJacobian;
  //! Whether dynamic Jacobian (w.r. to endogenous and exogenous) should be written
  bool computeJacobianExo;
  //! Whether dynamic Hessian (w.r. to endogenous and exogenous) should be written
  bool computeHessian;
  //! Whether static Hessian (w.r. to endogenous only) should be written
  bool computeStaticHessian;
  //! Whether dynamic third order derivatives (w.r. to endogenous and exogenous) should be written
  bool computeThirdDerivatives;
  //! Execute computations (variable sorting + derivation)
  /*! You must set computeJacobian, computeJacobianExo, computeHessian, computeStaticHessian and computeThirdDerivatives to correct values before calling this function */
  void computingPass(const eval_context_type &eval_context);
  //! Writes model initialization and lead/lag incidence matrix to output
  void writeOutput(ostream &output) const;
  //! Writes static model file
  void writeStaticFile(const string &basename) const;
  //! Writes dynamic model file
  void writeDynamicFile(const string &basename) const;
  //! Complete set to block decompose the model
  BlockTriangular block_triangular;
  //! Adds informations for simulation in a binary file
  void Write_Inf_To_Bin_File(const string &dynamic_basename, const string &bin_basename,
                             const int &num, int &u_count_int, bool &file_open) const;

  int equation_number() const;
};

#endif
