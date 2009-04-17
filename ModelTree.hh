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

#ifndef _MODELTREE_HH
#define _MODELTREE_HH

using namespace std;

#include <string>
#include <vector>
#include <map>
#include <ostream>

#include "DataTree.hh"

//! The three in which ModelTree can work
enum ModelTreeMode
  {
    eStandardMode, //!< Standard mode (static and dynamic files in Matlab)
    eSparseMode,  //!< Sparse mode (static file in Matlab, dynamic file in Matlab with block decomposition)
    eDLLMode,      //!< DLL mode (static and dynamic files in C)
    eSparseDLLMode //!< Sparse DLL mode (static file in Matlab, dynamic file in C with block decomposition plus a binary file)
  };

//! Shared code for static and dynamic models
class ModelTree : public DataTree
{
protected:
  //! Stores declared equations
  vector<BinaryOpNode *> equations;

  typedef map<pair<int, int>, NodeID> first_derivatives_type;
  //! First order derivatives
  /*! First index is equation number, second is variable w.r. to which is computed the derivative.
    Only non-null derivatives are stored in the map.
    Variable indices are those of the getDerivID() method.
  */
  first_derivatives_type first_derivatives;

  typedef map<pair<int, pair<int, int> >, NodeID> second_derivatives_type;
  //! Second order derivatives
  /*! First index is equation number, second and third are variables w.r. to which is computed the derivative.
    Only non-null derivatives are stored in the map.
    Contains only second order derivatives where var1 >= var2 (for obvious symmetry reasons).
    Variable indices are those of the getDerivID() method.
  */
  second_derivatives_type second_derivatives;

  typedef map<pair<int, pair<int, pair<int, int> > >, NodeID> third_derivatives_type;
  //! Third order derivatives
  /*! First index is equation number, second, third and fourth are variables w.r. to which is computed the derivative.
    Only non-null derivatives are stored in the map.
    Contains only third order derivatives where var1 >= var2 >= var3 (for obvious symmetry reasons).
    Variable indices are those of the getDerivID() method.
  */
  third_derivatives_type third_derivatives;

  //! Temporary terms (those which will be noted Txxxx)
  temporary_terms_type temporary_terms;

  //! Computes derivatives of ModelTree
  void derive(int order);
  //! Write derivative of an equation w.r. to a variable
  void writeDerivative(ostream &output, int eq, int symb_id, int lag, ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms) const;
  //! Computes temporary terms
  void computeTemporaryTerms(int order);
  //! Writes temporary terms
  void writeTemporaryTerms(ostream &output, ExprNodeOutputType output_type) const;
  //! Writes model local variables
  /*! No temporary term is used in the output, so that local parameters declarations can be safely put before temporary terms declaration in the output files */
  void writeModelLocalVariables(ostream &output, ExprNodeOutputType output_type) const;
  //! Writes model equations
  void writeModelEquations(ostream &output, ExprNodeOutputType output_type) const;

  //! Writes either (i+1,j+1) or [i+j*n_i] whether we are in Matlab or C mode
  void matrixHelper(ostream &output, int eq_nb, int col_nb, ExprNodeOutputType output_type) const;

public:
  ModelTree(SymbolTable &symbol_table_arg, NumericalConstants &num_constants);
  //! Mode in which the ModelTree is supposed to work (Matlab, DLL or SparseDLL)
  ModelTreeMode mode;
  //! Declare a node as an equation of the model
  void addEquation(NodeID eq);
  //! Returns the number of equations in the model
  int equation_number() const;
};

#endif
