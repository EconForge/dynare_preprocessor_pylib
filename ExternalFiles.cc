/*
 * Copyright (C) 2006-2013 Dynare Team
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
#include <cassert>
#include "ModFile.hh"
#include "DynamicModel.hh"
#include "StaticModel.hh"
#include "SteadyStateModel.hh"

void
ModFile::writeExternalFiles(const string &basename, FileOutputType output, LanguageOutputType language) const
{
  switch(language)
    {
    case c:
      writeExternalFilesC(basename, output);
      break;
    case cpp:
      writeExternalFilesCC(basename, output);
      break;
    default:
      cerr << "This case shouldn't happen. Contact the authors of Dynare" << endl;
      exit(EXIT_FAILURE);
    }
}

// C interface

void
ModFile::writeExternalFilesC(const string &basename, FileOutputType output) const
{
  writeModelC(basename);
  steady_state_model.writeSteadyStateFileC(basename, mod_file_struct.ramsey_model_present);

  dynamic_model.writeDynamicFile(basename, block, byte_code, use_dll, mod_file_struct.order_option);

  if (!no_static)
    static_model.writeStaticFile(basename, false, false, true);


  //  static_model.writeStaticCFile(basename, block, byte_code, use_dll);
  //  static_model.writeParamsDerivativesFileC(basename, cuda);
  //  static_model.writeAuxVarInitvalC(mOutputFile, oMatlabOutsideModel, cuda);

  // dynamic_model.writeResidualsC(basename, cuda);
  // dynamic_model.writeParamsDerivativesFileC(basename, cuda);
  dynamic_model.writeFirstDerivativesC(basename, cuda);
  
  if (output == second)
    dynamic_model.writeSecondDerivativesC_csr(basename, cuda);
  else if (output == third)
    {
        dynamic_model.writeSecondDerivativesC_csr(basename, cuda);
  	dynamic_model.writeThirdDerivativesC_csr(basename, cuda);
    }
}

void
ModFile::writeModelC(const string &basename) const
{
  string filename = basename + ".c";

  ofstream mDriverCFile;
  mDriverCFile.open(filename.c_str(), ios::out | ios::binary);
  if (!mDriverCFile.is_open())
    {
      cerr << "Error: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }

  mDriverCFile << "/*" << endl
               << " * " << filename << " : Driver file for Dynare C code" << endl
               << " *" << endl
               << " * Warning : this file is generated automatically by Dynare" << endl
               << " *           from model file (.mod)" << endl
               << " */" << endl
               << endl
               << "#include \"dynare_driver.h\"" << endl
               << endl
               << "struct" << endl
               << "{" << endl;

  // Write basic info
  symbol_table.writeCOutput(mDriverCFile);

  mDriverCFile << endl << "params.resize(param_nbr);" << endl;

  if (dynamic_model.equation_number() > 0)
    {
      dynamic_model.writeCOutput(mDriverCFile, basename, block, byte_code, use_dll, mod_file_struct.order_option, mod_file_struct.estimation_present);
      //      if (!no_static)
      //        static_model.writeCOutput(mOutputFile, block);
    }

  // Print statements
  for (vector<Statement *>::const_iterator it = statements.begin();
       it != statements.end(); it++)
      (*it)->writeCOutput(mDriverCFile, basename);

  mDriverCFile << "} DynareInfo;" << endl;
  mDriverCFile.close();

  // Write informational m file
  ofstream mOutputFile;

  if (basename.size())
    {
      string fname(basename);
      fname += ".m";
      mOutputFile.open(fname.c_str(), ios::out | ios::binary);
      if (!mOutputFile.is_open())
        {
          cerr << "ERROR: Can't open file " << fname
               << " for writing" << endl;
          exit(EXIT_FAILURE);
        }
    }
  else
    {
      cerr << "ERROR: Missing file name" << endl;
      exit(EXIT_FAILURE);
    }

  mOutputFile << "%" << endl
              << "% Status : informational m file" << endl
              << "%" << endl
              << "% Warning : this file is generated automatically by Dynare" << endl
              << "%           from model file (.mod)" << endl << endl
              << "disp('The following C file was successfully created:');" << endl
              << "ls preprocessorOutput.c" << endl << endl;
  mOutputFile.close();
}


void
DynamicModel::writeCOutput(ostream &output, const string &basename, bool block_decomposition, bool byte_code, bool use_dll, int order, bool estimation_present) const
{
  int lag_presence[3];
  // Loop on endogenous variables
  vector<int> zeta_back, zeta_mixed, zeta_fwrd, zeta_static;
  for (int endoID = 0; endoID < symbol_table.endo_nbr(); endoID++)
    {
      // Loop on periods
      for (int lag = 0; lag <= 2; lag++)
	{
	  lag_presence[lag] = 1;
          try
            {
              getDerivID(symbol_table.getID(eEndogenous, endoID), lag-1);
            }
          catch (UnknownDerivIDException &e)
            {
	      lag_presence[lag] = 0;
            }
        }
      if (lag_presence[0] == 1)
	if (lag_presence[2] == 1)
	  zeta_mixed.push_back(endoID);
	else
	  zeta_back.push_back(endoID);
      else if (lag_presence[2] == 1)
	zeta_fwrd.push_back(endoID);
      else
	zeta_static.push_back(endoID);
      
    }
  output << "size_t nstatic = " << zeta_static.size() << ";" << endl
         << "size_t nfwrd   = " << zeta_fwrd.size() << ";" << endl
         << "size_t nback   = " << zeta_back.size() << ";" << endl
         << "size_t nmixed  = " << zeta_mixed.size() << ";" << endl;
  output << "size_t zeta_static[" << zeta_static.size() << "] = {";
  for (vector<int>::iterator i = zeta_static.begin(); i != zeta_static.end(); ++i)
    {
      if ( i != zeta_static.begin() )
	output << ",";
      output << *i;
    }
  output << "};" << endl;

  output << "size_t zeta_back[" << zeta_back.size() << "] = {";
  for (vector<int>::iterator i = zeta_back.begin(); i != zeta_back.end(); ++i)
    {
      if ( i != zeta_back.begin() )
	output << ",";
      output << *i;
    }
  output << "};" << endl;

  output << "size_t zeta_fwrd[" << zeta_fwrd.size() << "] = {";
  for (vector<int>::iterator i = zeta_fwrd.begin(); i != zeta_fwrd.end(); ++i)
    {
      if ( i != zeta_fwrd.begin() )
	output << ",";
      output << *i;
    }
  output << "};" << endl;

  output << "size_t zeta_mixed[" << zeta_mixed.size() << "] = {";
  for (vector<int>::iterator i = zeta_mixed.begin(); i != zeta_mixed.end(); ++i)
    {
      if ( i != zeta_mixed.begin() )
	output << ",";
      output << *i;
    }
  output << "};" << endl;

  // Write number of non-zero derivatives
  // Use -1 if the derivatives have not been computed
  output << "int *NNZDerivatives[3] = {";
  switch (order)
    {
    case 0:
      output << NNZDerivatives[0] << ",-1,-1};" << endl;
      break;
    case 1:
      output << NNZDerivatives[0] << "," << NNZDerivatives[1] << ",-1};" << endl;
      break;
    case 2:
      output << NNZDerivatives[0] << "," << NNZDerivatives[1] << "," << NNZDerivatives[2] << "};" << endl;
      break;
    default:
	cerr << "Order larger than 3 not implemented" << endl;
	exit(EXIT_FAILURE);
    }
}

void
DynamicModel::writeFirstDerivativesC(const string &basename, bool cuda) const
{
  string filename = basename + "_first_derivatives.c";
  ofstream mDynamicModelFile, mDynamicMexFile;

  mDynamicModelFile.open(filename.c_str(), ios::out | ios::binary);
  if (!mDynamicModelFile.is_open())
    {
      cerr << "Error: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }
  mDynamicModelFile << "/*" << endl
                    << " * " << filename << " : Computes first order derivatives of the model for Dynare" << endl
                    << " *" << endl
                    << " * Warning : this file is generated automatically by Dynare" << endl
                    << " *           from model " << basename << "(.mod)" << endl
                    << " */" << endl
                    << "#include <math.h>" << endl;

  mDynamicModelFile << "#include <stdlib.h>" << endl;

  mDynamicModelFile << "#define max(a, b) (((a) > (b)) ? (a) : (b))" << endl
                    << "#define min(a, b) (((a) > (b)) ? (b) : (a))" << endl;

  // Write function definition if oPowerDeriv is used
  writePowerDerivCHeader(mDynamicModelFile);

  mDynamicModelFile << "void FirstDerivatives(const double *y, double *x, int nb_row_x, double *params, double *steady_state, int it_, double *residual, double *g1, double *v2, double *v3)" << endl
                    << "{" << endl;

  // this is always empty here, but needed by d1->writeOutput
  deriv_node_temp_terms_t tef_terms;

  // Writing Jacobian
  for (first_derivatives_t::const_iterator it = first_derivatives.begin();
       it != first_derivatives.end(); it++)
    {
      int eq = it->first.first;
      int var = it->first.second;
      expr_t d1 = it->second;

      jacobianHelper(mDynamicModelFile, eq, getDynJacobianCol(var), oCDynamicModel);
      mDynamicModelFile << "=";
      // oCstaticModel makes reference to the static variables
      d1->writeOutput(mDynamicModelFile, oCStaticModel, temporary_terms, tef_terms);
      mDynamicModelFile << ";" << endl;
    }
  
  mDynamicModelFile << "}" << endl;

  writePowerDeriv(mDynamicModelFile, true);
  mDynamicModelFile.close();

}

// using compressed sparse row format (CSR)
void
DynamicModel::writeSecondDerivativesC_csr(const string &basename, bool cuda) const
{

  string filename = basename + "_second_derivatives.c";
  ofstream mDynamicModelFile, mDynamicMexFile;

  mDynamicModelFile.open(filename.c_str(), ios::out | ios::binary);
  if (!mDynamicModelFile.is_open())
    {
      cerr << "Error: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }
  mDynamicModelFile << "/*" << endl
                    << " * " << filename << " : Computes second order derivatives of the model for Dynare" << endl
                    << " *" << endl
                    << " * Warning : this file is generated automatically by Dynare" << endl
                    << " *           from model " << basename << "(.mod)" << endl
                    << " */" << endl
                    << "#include <math.h>" << endl;

  mDynamicModelFile << "#include <stdlib.h>" << endl;

  mDynamicModelFile << "#define max(a, b) (((a) > (b)) ? (a) : (b))" << endl
                    << "#define min(a, b) (((a) > (b)) ? (b) : (a))" << endl;

  // write function definition if oPowerDeriv is used
  writePowerDerivCHeader(mDynamicModelFile);

  mDynamicModelFile << "void SecondDerivatives(const double *y, double *x, int nb_row_x, double *params, double *steady_state, int it_, double *residual, int *row_ptr, int *col_ptr, double *value)" << endl
                    << "{" << endl;

  // this is always empty here, but needed by d1->writeOutput
  deriv_node_temp_terms_t tef_terms;

  // Indexing derivatives in column order
  vector<derivative> D;
  int hessianColsNbr = dynJacobianColsNbr*dynJacobianColsNbr;
  for (second_derivatives_t::const_iterator it = second_derivatives.begin();
       it != second_derivatives.end(); it++)
    {
      int eq = it->first.first;
      int var1 = it->first.second.first;
      int var2 = it->first.second.second;

      int id1 = getDynJacobianCol(var1);
      int id2 = getDynJacobianCol(var2);

      int col_nb = id1 * dynJacobianColsNbr + id2;

      derivative deriv(col_nb + eq*hessianColsNbr,col_nb,eq,it->second);
      D.push_back(deriv);
      if (id1 != id2)
	{
	  col_nb = id2 * dynJacobianColsNbr + id1;
	  derivative deriv(col_nb + eq*hessianColsNbr,col_nb,eq,it->second);
	  D.push_back(deriv);
	}
    }
  sort(D.begin(), D.end(), derivative_less_than() );

  // Writing Hessian
  vector<int> row_ptr(equations.size());
  fill(row_ptr.begin(),row_ptr.end(),0.0);
  int k = 0;
  for(vector<derivative>::const_iterator it = D.begin(); it != D.end(); ++it)
    {
      row_ptr[it->row_nbr]++;
      mDynamicModelFile << "col_ptr[" << k << "] "
			<< "=" << it->col_nbr << ";" << endl;
      mDynamicModelFile << "value[" << k << "] = ";
      // oCstaticModel makes reference to the static variables
      it->value->writeOutput(mDynamicModelFile, oCStaticModel, temporary_terms, tef_terms);
      mDynamicModelFile << ";" << endl;
      k++;
    }
  
  // row_ptr must point to the relative address of the first element of the row
  int cumsum = 0;
  mDynamicModelFile << "row_ptr = [ 0";
  for (vector<int>::iterator it=row_ptr.begin(); it != row_ptr.end(); ++it)
    {
      cumsum += *it;
      mDynamicModelFile << ", " << cumsum;
    }
  mDynamicModelFile << "];" << endl;   

  mDynamicModelFile << "}" << endl;

  writePowerDeriv(mDynamicModelFile, true);
  mDynamicModelFile.close();

}

void
DynamicModel::writeThirdDerivativesC_csr(const string &basename, bool cuda) const
{
  string filename = basename + "_third_derivatives.c";
  ofstream mDynamicModelFile, mDynamicMexFile;

  mDynamicModelFile.open(filename.c_str(), ios::out | ios::binary);
  if (!mDynamicModelFile.is_open())
    {
      cerr << "Error: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }
  mDynamicModelFile << "/*" << endl
                    << " * " << filename << " : Computes third order derivatives of the model for Dynare" << endl
                    << " *" << endl
                    << " * Warning : this file is generated automatically by Dynare" << endl
                    << " *           from model " << basename << "(.mod)" << endl
                    << " */" << endl
                    << "#include <math.h>" << endl;

  mDynamicModelFile << "#include <stdlib.h>" << endl;

  mDynamicModelFile << "#define max(a, b) (((a) > (b)) ? (a) : (b))" << endl
                    << "#define min(a, b) (((a) > (b)) ? (b) : (a))" << endl;

  // Write function definition if oPowerDeriv is used
  writePowerDerivCHeader(mDynamicModelFile);

  mDynamicModelFile << "void ThirdDerivatives(const double *y, double *x, int nb_row_x, double *params, double *steady_state, int it_, double *residual, double *g1, double *v2, double *v3)" << endl
                    << "{" << endl;

  // this is always empty here, but needed by d1->writeOutput
  deriv_node_temp_terms_t tef_terms;

  vector<derivative> D;
  int hessianColsNbr = dynJacobianColsNbr*dynJacobianColsNbr;
  int thirdDerivativesColsNbr = hessianColsNbr*dynJacobianColsNbr;
  for (third_derivatives_t::const_iterator it = third_derivatives.begin();
       it != third_derivatives.end(); it++)
    {
      int eq = it->first.first;
      int var1 = it->first.second.first;
      int var2 = it->first.second.second.first;
      int var3 = it->first.second.second.second;

      int id1 = getDynJacobianCol(var1);
      int id2 = getDynJacobianCol(var2);
      int id3 = getDynJacobianCol(var3);

      // Reference column number for the g3 matrix (with symmetrical derivatives)
      vector<long unsigned int>  cols;
      long unsigned int col_nb = id1 * hessianColsNbr + id2 * dynJacobianColsNbr + id3;
      int thirdDColsNbr = hessianColsNbr*dynJacobianColsNbr;
      derivative deriv(col_nb + eq*thirdDColsNbr,col_nb,eq,it->second);
      D.push_back(deriv);
      cols.push_back(col_nb);
      col_nb = id1 * hessianColsNbr + id3 * dynJacobianColsNbr + id2;
      if (find(cols.begin(),cols.end(),col_nb) == cols.end())
	{
	  derivative deriv(col_nb + eq*thirdDerivativesColsNbr,col_nb,eq,it->second);
	  D.push_back(deriv);
	  cols.push_back(col_nb);
	}
      col_nb = id2 * hessianColsNbr + id1 * dynJacobianColsNbr + id3;
      if (find(cols.begin(),cols.end(),col_nb) == cols.end())
	{
	  derivative deriv(col_nb + eq*thirdDerivativesColsNbr,col_nb,eq,it->second);
	  D.push_back(deriv);
	  cols.push_back(col_nb);
	}
      col_nb = id2 * hessianColsNbr + id3 * dynJacobianColsNbr + id1;
      if (find(cols.begin(),cols.end(),col_nb) == cols.end())
	{
	  derivative deriv(col_nb + eq*thirdDerivativesColsNbr,col_nb,eq,it->second);
	  D.push_back(deriv);
	  cols.push_back(col_nb);
	}
      col_nb = id3 * hessianColsNbr + id1 * dynJacobianColsNbr + id2;
      if (find(cols.begin(),cols.end(),col_nb) == cols.end())
	{
	  derivative deriv(col_nb + eq*thirdDerivativesColsNbr,col_nb,eq,it->second);
	  D.push_back(deriv);
	  cols.push_back(col_nb);
	}
      col_nb = id3 * hessianColsNbr + id2 * dynJacobianColsNbr + id1;
      if (find(cols.begin(),cols.end(),col_nb) == cols.end())
	{
	  derivative deriv(col_nb + eq*thirdDerivativesColsNbr,col_nb,eq,it->second);
	  D.push_back(deriv);
	}
    }

  sort(D.begin(), D.end(), derivative_less_than() );

  vector<int> row_ptr(equations.size());
  fill(row_ptr.begin(),row_ptr.end(),0.0);
  int k = 0;
  for(vector<derivative>::const_iterator it = D.begin(); it != D.end(); ++it)
    {
      row_ptr[it->row_nbr]++;
      mDynamicModelFile << "col_ptr[" << k << "] "
			<< "=" << it->col_nbr << ";" << endl;
      mDynamicModelFile << "value[" << k << "] = ";
      // oCstaticModel makes reference to the static variables
      it->value->writeOutput(mDynamicModelFile, oCStaticModel, temporary_terms, tef_terms);
      mDynamicModelFile << ";" << endl;
      k++;
    }

  // row_ptr must point to the relative address of the first element of the row
  int cumsum = 0;
  mDynamicModelFile << "row_ptr = [ 0";
  for (vector<int>::iterator it=row_ptr.begin(); it != row_ptr.end(); ++it)
    {
      cumsum += *it;
      mDynamicModelFile << ", " << cumsum;
    }
  mDynamicModelFile << "];" << endl;   

  mDynamicModelFile << "}" << endl;

  writePowerDeriv(mDynamicModelFile, true);
  mDynamicModelFile.close();

}

void
SteadyStateModel::writeSteadyStateFileC(const string &basename, bool ramsey_model) const
{
  string filename = basename + "_steadystate.c";

  ofstream output;
  output.open(filename.c_str(), ios::out | ios::binary);
  if (!output.is_open())
    {
      cerr << "ERROR: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }

  output << "#include <math.h>" << endl;

  output << "void steadystate("
	 << "const double *exo_, const double *params, double *ys_, int *info)" << endl
         << "// Steady state file generated by Dynare preprocessor" << endl
	 << "{" << endl
         << "    *info = 0;" << endl;

  if (def_table.size() == 0)
    {
      output << "    return;" << endl
	     << "}" << endl;
      return;
    }

  for (size_t i = 0; i < def_table.size(); i++)
    {
      const vector<int> &symb_ids = def_table[i].first;
      output << "    ";
      if (symb_ids.size() > 1)
	std::cout << "Error: in C, multiple returns are not permitted in steady_state_model" << std::endl;
      variable_node_map_t::const_iterator it = variable_node_map.find(make_pair(symb_ids[0], 0));
      assert(it != variable_node_map.end());
      if (it->second->get_type() == eModFileLocalVariable)
	output << "double ";
      dynamic_cast<ExprNode *>(it->second)->writeOutput(output, oCSteadyStateFile);
      output << "=";
      def_table[i].second->writeOutput(output, oCSteadyStateFile);
      output << ";" << endl;
    }
  output << "    // Auxiliary equations" << endl;
  static_model.writeAuxVarInitval(output, oCSteadyStateFile);
  output << "}" << endl;
}


//
// C++ interface
//
void
ModFile::writeExternalFilesCC(const string &basename, FileOutputType output) const
{
  writeModelCC(basename);
  steady_state_model.writeSteadyStateFileC(basename, mod_file_struct.ramsey_model_present);

  dynamic_model.writeDynamicFile(basename, block, byte_code, use_dll, mod_file_struct.order_option);

  if (!no_static)
    static_model.writeStaticFile(basename, false, false, true);


  //  static_model.writeStaticCFile(basename, block, byte_code, use_dll);
  //  static_model.writeParamsDerivativesFileC(basename, cuda);
  //  static_model.writeAuxVarInitvalC(mOutputFile, oMatlabOutsideModel, cuda);

  // dynamic_model.writeResidualsC(basename, cuda);
  // dynamic_model.writeParamsDerivativesFileC(basename, cuda);
  dynamic_model.writeFirstDerivativesC(basename, cuda);
  
  if (output == second)
    dynamic_model.writeSecondDerivativesC_csr(basename, cuda);
  else if (output == third)
    {
        dynamic_model.writeSecondDerivativesC_csr(basename, cuda);
  	dynamic_model.writeThirdDerivativesC_csr(basename, cuda);
    }
}

void
ModFile::writeModelCC(const string &basename) const
{
  string filename = basename + ".cc";

  ofstream mDriverCFile;
  mDriverCFile.open(filename.c_str(), ios::out | ios::binary);
  if (!mDriverCFile.is_open())
    {
      cerr << "Error: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }

  mDriverCFile << "/*" << endl
               << " * " << filename << " : Driver file for Dynare C++ code" << endl
               << " *" << endl
               << " * Warning : this file is generated automatically by Dynare" << endl
               << " *           from model file (.mod)" << endl
               << " */" << endl
               << endl
               << "#include \"dynare_cpp_driver.hh\"" << endl
               << endl
               << "DynareInfo::DynareInfo(void)" << endl
               << "{" << endl;

  // Write basic info
  symbol_table.writeCCOutput(mDriverCFile);

  mDriverCFile << endl << "params.resize(param_nbr);" << endl;

  if (dynamic_model.equation_number() > 0)
    {
      dynamic_model.writeCCOutput(mDriverCFile, basename, block, byte_code, use_dll, mod_file_struct.order_option, mod_file_struct.estimation_present);
      //      if (!no_static)
      //        static_model.writeCOutput(mOutputFile, block);
    }

  // Print statements
  for (vector<Statement *>::const_iterator it = statements.begin();
       it != statements.end(); it++)
      (*it)->writeCOutput(mDriverCFile, basename);

  mDriverCFile << "};" << endl;
  mDriverCFile.close();

  // Write informational m file
  ofstream mOutputFile;

  if (basename.size())
    {
      string fname(basename);
      fname += ".m";
      mOutputFile.open(fname.c_str(), ios::out | ios::binary);
      if (!mOutputFile.is_open())
        {
          cerr << "ERROR: Can't open file " << fname
               << " for writing" << endl;
          exit(EXIT_FAILURE);
        }
    }
  else
    {
      cerr << "ERROR: Missing file name" << endl;
      exit(EXIT_FAILURE);
    }

  mOutputFile << "%" << endl
              << "% Status : informational m file" << endl
              << "%" << endl
              << "% Warning : this file is generated automatically by Dynare" << endl
              << "%           from model file (.mod)" << endl << endl
              << "disp('The following C++ file was successfully created:');" << endl
              << "ls preprocessorOutput.cc" << endl << endl;
  mOutputFile.close();
}


void
DynamicModel::writeCCOutput(ostream &output, const string &basename, bool block_decomposition, bool byte_code, bool use_dll, int order, bool estimation_present) const
{
  int lag_presence[3];
  // Loop on endogenous variables
  for (int endoID = 0; endoID < symbol_table.endo_nbr(); endoID++)
    {
      // Loop on periods
      for (int lag = 0; lag <= 2; lag++)
	{
	  lag_presence[lag] = 1;
          try
            {
              getDerivID(symbol_table.getID(eEndogenous, endoID), lag-1);
            }
          catch (UnknownDerivIDException &e)
            {
	      lag_presence[lag] = 0;
            }
        }
      if (lag_presence[0] == 1)
	if (lag_presence[2] == 1)
	  output << "zeta_mixed.push_back(" << endoID << ");" << endl;
	else
	  output << "zeta_back.push_back(" << endoID << ");" << endl;
      else if (lag_presence[2] == 1)
	output << "zeta_fwrd.push_back(" << endoID << ");" << endl;
      else
	output << "zeta_static.push_back(" << endoID << ");" << endl;
      
    }
  output << "nstatic = zeta_static.size();" << endl
         << "nfwrd   = zeta_fwrd.size();" << endl
         << "nback   = zeta_back.size();" << endl
         << "nmixed  = zeta_mixed.size();" << endl;

  // Write number of non-zero derivatives
  // Use -1 if the derivatives have not been computed
  output << endl
         << "NNZDerivatives.push_back(" << NNZDerivatives[0] << ");" << endl;
  if (order > 1)
    {
      output << "NNZDerivatives.push_back(" << NNZDerivatives[1] << ");" << endl;
      if (order > 2)
        output << "NNZDerivatives.push_back(" << NNZDerivatives[2] << ");" << endl;
      else
        output << "NNZDerivatives.push_back(-1);" << endl;
    }
  else
    output << "NNZDerivatives.push_back(-1);" << endl
           << "NNZDerivatives.push_back(-1);" << endl;
}


