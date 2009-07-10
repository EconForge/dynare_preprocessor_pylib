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

//! Shared code for static and dynamic models
class ModelTree : public DataTree
{
protected:
  //! Stores declared equations
  vector<BinaryOpNode *> equations;

  //! Number of non-zero derivatives
  int NNZDerivatives[3];

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

  //! Computes 1st derivatives
  /*! \param vars the derivation IDs w.r. to which compute the derivatives */
  void computeJacobian(const set<int> &vars);
  //! Computes 2nd derivatives
  /*! \param vars the derivation IDs w.r. to which derive the 1st derivatives */
  void computeHessian(const set<int> &vars);
  //! Computes 3rd derivatives
  /*! \param vars the derivation IDs w.r. to which derive the 2nd derivatives */
  void computeThirdDerivatives(const set<int> &vars);

  //! Write derivative of an equation w.r. to a variable
  void writeDerivative(ostream &output, int eq, int symb_id, int lag, ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms) const;
  //! Computes temporary terms (for all equations and derivatives)
  void computeTemporaryTerms(bool is_matlab);
  //! Writes temporary terms
  void writeTemporaryTerms(const temporary_terms_type &tt, ostream &output, ExprNodeOutputType output_type) const;
  //! Writes model local variables
  /*! No temporary term is used in the output, so that local parameters declarations can be safely put before temporary terms declaration in the output files */
  void writeModelLocalVariables(ostream &output, ExprNodeOutputType output_type) const;
  //! Writes model equations
  void writeModelEquations(ostream &output, ExprNodeOutputType output_type) const;

  //! Writes LaTeX model file
  void writeLatexModelFile(const string &filename, ExprNodeOutputType output_type) const;

public:
  ModelTree(SymbolTable &symbol_table_arg, NumericalConstants &num_constants);
  //! Declare a node as an equation of the model
  void addEquation(NodeID eq);
  //! Returns the number of equations in the model
  int equation_number() const;
};

#endif
