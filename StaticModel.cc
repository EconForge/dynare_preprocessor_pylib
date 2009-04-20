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

#include <cstdlib>

#include "StaticModel.hh"

StaticModel::StaticModel(SymbolTable &symbol_table_arg,
                         NumericalConstants &num_constants_arg) :
  ModelTree(symbol_table_arg, num_constants_arg)
{
}

void
StaticModel::writeStaticMFile(const string &static_basename) const
{
  string filename = static_basename + ".m";

  ofstream mStaticModelFile;
  mStaticModelFile.open(filename.c_str(), ios::out | ios::binary);
  if (!mStaticModelFile.is_open())
    {
      cerr << "Error: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }
  // Writing comments and function definition command
  mStaticModelFile << "function [residual, g1, g2] = " << static_basename << "(y, x, params)" << endl
                   << "%" << endl
                   << "% Status : Computes static model for Dynare" << endl
                   << "%" << endl
                   << "% Warning : this file is generated automatically by Dynare" << endl
                   << "%           from model file (.mod)" << endl << endl;

  writeStaticModel(mStaticModelFile);

  mStaticModelFile.close();
}

void
StaticModel::writeStaticCFile(const string &static_basename) const
{
  string filename = static_basename + ".c";

  ofstream mStaticModelFile;
  mStaticModelFile.open(filename.c_str(), ios::out | ios::binary);
  if (!mStaticModelFile.is_open())
    {
      cerr << "Error: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }
  mStaticModelFile << "/*" << endl
                   << " * " << filename << " : Computes static model for Dynare" << endl
                   << " * Warning : this file is generated automatically by Dynare" << endl
                   << " *           from model file (.mod)" << endl
                   << endl
                   << " */" << endl
                   << "#include <math.h>" << endl
                   << "#include \"mex.h\"" << endl;

  // Writing the function Static
  writeStaticModel(mStaticModelFile);

  // Writing the gateway routine
  mStaticModelFile << "/* The gateway routine */" << endl
                   << "void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])" << endl
                   << "{" << endl
                   << "  double *y, *x, *params;" << endl
                   << "  double *residual, *g1;" << endl
                   << endl
                   << "  /* Create a pointer to the input matrix y. */" << endl
                   << "  y = mxGetPr(prhs[0]);" << endl
                   << endl
                   << "  /* Create a pointer to the input matrix x. */" << endl
                   << "  x = mxGetPr(prhs[1]);" << endl
                   << endl
                   << "  /* Create a pointer to the input matrix params. */" << endl
                   << "  params = mxGetPr(prhs[2]);" << endl
                   << endl
                   << "  residual = NULL;" << endl
                   << "  if (nlhs >= 1)" << endl
                   << "  {" << endl
                   << "      /* Set the output pointer to the output matrix residual. */" << endl
                   << "      plhs[0] = mxCreateDoubleMatrix(" << equations.size() << ",1, mxREAL);" << endl
                   << "     /* Create a C pointer to a copy of the output matrix residual. */" << endl
                   << "     residual = mxGetPr(plhs[0]);" << endl
                   << "  }" << endl
                   << endl
                   << "  g1 = NULL;" << endl
                   << "  if (nlhs >= 2)" << endl
                   << "  {" << endl
                   << "      /* Set the output pointer to the output matrix g1. */" << endl
                   << "      plhs[1] = mxCreateDoubleMatrix(" << equations.size() << ", " << symbol_table.endo_nbr() << ", mxREAL);" << endl
                   << "      /* Create a C pointer to a copy of the output matrix g1. */" << endl
                   << "      g1 = mxGetPr(plhs[1]);" << endl
                   << "  }" << endl
                   << endl
                   << "  /* Call the C Static. */" << endl
                   << "  Static(y, x, params, residual, g1);" << endl
                   << "}" << endl;

  mStaticModelFile.close();
}

void
StaticModel::writeStaticModel(ostream &StaticOutput) const
{
  ostringstream model_output;    // Used for storing model equations
  ostringstream jacobian_output; // Used for storing jacobian equations
  ostringstream hessian_output;
  ostringstream lsymetric;       // For symmetric elements in hessian

  ExprNodeOutputType output_type = (mode == eDLLMode ? oCStaticModel : oMatlabStaticModel);

  writeModelLocalVariables(model_output, output_type);

  writeTemporaryTerms(temporary_terms, model_output, output_type);

  writeModelEquations(model_output, output_type);

  // Write Jacobian w.r. to endogenous only
  for (first_derivatives_type::const_iterator it = first_derivatives.begin();
       it != first_derivatives.end(); it++)
    {
      int eq = it->first.first;
      int symb_id = inv_deriv_id_table[it->first.second];
      NodeID d1 = it->second;

      ostringstream g1;
      g1 << "  g1";
      matrixHelper(g1, eq, symbol_table.getTypeSpecificID(symb_id), output_type);

      jacobian_output << g1.str() << "=" << g1.str() << "+";
      d1->writeOutput(jacobian_output, output_type, temporary_terms);
      jacobian_output << ";" << endl;
    }

  // Write Hessian w.r. to endogenous only (only if 2nd order derivatives have been computed)
  for (second_derivatives_type::const_iterator it = second_derivatives.begin();
       it != second_derivatives.end(); it++)
    {
      int eq = it->first.first;
      int symb_id1 = inv_deriv_id_table[it->first.second.first];
      int symb_id2 = inv_deriv_id_table[it->first.second.second];
      NodeID d2 = it->second;

      int tsid1 = symbol_table.getTypeSpecificID(symb_id1);
      int tsid2 = symbol_table.getTypeSpecificID(symb_id2);

      int col_nb = tsid1*symbol_table.endo_nbr()+tsid2;
      int col_nb_sym = tsid2*symbol_table.endo_nbr()+tsid1;

      hessian_output << "  g2";
      matrixHelper(hessian_output, eq, col_nb, output_type);
      hessian_output << " = ";
      d2->writeOutput(hessian_output, output_type, temporary_terms);
      hessian_output << ";" << endl;

      // Treating symetric elements
      if (symb_id1 != symb_id2)
        {
          lsymetric <<  "  g2";
          matrixHelper(lsymetric, eq, col_nb_sym, output_type);
          lsymetric << " = " <<  "g2";
          matrixHelper(lsymetric, eq, col_nb, output_type);
          lsymetric << ";" << endl;
        }
    }

  // Writing ouputs
  if (mode != eDLLMode)
    {
      StaticOutput << "residual = zeros( " << equations.size() << ", 1);" << endl << endl
                   << "%" << endl
                   << "% Model equations" << endl
                   << "%" << endl
                   << endl
                   << model_output.str()
                   << "if ~isreal(residual)" << endl
                   << "  residual = real(residual)+imag(residual).^2;" << endl
                   << "end" << endl
                   << "if nargout >= 2," << endl
                   << "  g1 = zeros(" << equations.size() << ", " << symbol_table.endo_nbr() << ");" << endl
                   << endl
                   << "%" << endl
                   << "% Jacobian matrix" << endl
                   << "%" << endl
                   << endl
                   << jacobian_output.str()
                   << "  if ~isreal(g1)" << endl
                   << "    g1 = real(g1)+2*imag(g1);" << endl
                   << "  end" << endl
                   << "end" << endl;

      // If 2nd order derivatives have been computed
      if (second_derivatives.size())
        {
          StaticOutput << "if nargout >= 3," << endl;
          // Writing initialization instruction for matrix g2
          int ncols = symbol_table.endo_nbr() * symbol_table.endo_nbr();
          StaticOutput << "  g2 = sparse([],[],[], " << equations.size() << ", " << ncols << ", " << 5*ncols << ");" << endl
                       << endl
                       << "%" << endl
                       << "% Hessian matrix" << endl
                       << "%" << endl
                       << endl
                       << hessian_output.str()
                       << lsymetric.str()
                       << "end;" << endl;
        }
    }
  else
    {
      StaticOutput << "void Static(double *y, double *x, double *params, double *residual, double *g1)" << endl
                   << "{" << endl
                   << "  double lhs, rhs;" << endl
        // Writing residual equations
                   << "  /* Residual equations */" << endl
                   << "  if (residual == NULL)" << endl
                   << "    return;" << endl
                   << "  else" << endl
                   << "    {" << endl
                   << model_output.str()
        // Writing Jacobian
                   << "     /* Jacobian for endogenous variables without lag */" << endl
                   << "     if (g1 == NULL)" << endl
                   << "       return;" << endl
                   << "     else" << endl
                   << "       {" << endl
                   << jacobian_output.str()
                   << "       }" << endl
                   << "    }" << endl
                   << "}" << endl << endl;
    }
}

void
StaticModel::writeStaticFile(const string &basename) const
{
  switch (mode)
    {
    case eStandardMode:
    case eSparseDLLMode:
    case eSparseMode:
      writeStaticMFile(basename + "_static");
      break;
    case eDLLMode:
      writeStaticCFile(basename + "_static");
      break;
    }
}

void
StaticModel::computingPass(bool hessian, bool no_tmp_terms)
{
  // Compute derivatives w.r. to all derivation IDs (i.e. all endogenous)
  set<int> vars;
  for(int i = 0; i < getDerivIDNbr(); i++)
    vars.insert(i);

  // Launch computations
  cout << "Computing static model derivatives:" << endl
       << " - order 1" << endl;
  computeJacobian(vars);

  if (hessian)
    {
      cout << " - order 2" << endl;
      computeHessian(vars);
    }

  if (!no_tmp_terms)
    computeTemporaryTerms();
}

int
StaticModel::computeDerivID(int symb_id, int lag)
{
  // Only create derivation ID for endogenous
  if (symbol_table.getType(symb_id) != eEndogenous)
    return -1;

  deriv_id_table_t::const_iterator it = deriv_id_table.find(symb_id);
  if (it != deriv_id_table.end())
    return it->second;

  // Create a new deriv_id
  int deriv_id = deriv_id_table.size();

  deriv_id_table[symb_id] = deriv_id;
  inv_deriv_id_table.push_back(symb_id);

  return deriv_id;
}

int
StaticModel::getDerivID(int symb_id, int lag) const throw (UnknownDerivIDException)
{
  deriv_id_table_t::const_iterator it = deriv_id_table.find(symb_id);
  if (it == deriv_id_table.end())
    throw UnknownDerivIDException();
  else
    return it->second;
}

int
StaticModel::getDerivIDNbr() const
{
  return deriv_id_table.size();
}
