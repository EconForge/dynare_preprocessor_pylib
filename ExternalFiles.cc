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
      int varID;
      // Loop on periods
      for (int lag = 0; lag <= 2; lag++)
	{
	  lag_presence[lag] = 1;
          try
            {
              varID = getDerivID(symbol_table.getID(eEndogenous, endoID), lag-1);
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

// C++ interface
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
      int varID;
      // Loop on periods
      for (int lag = 0; lag <= 2; lag++)
	{
	  lag_presence[lag] = 1;
          try
            {
              varID = getDerivID(symbol_table.getID(eEndogenous, endoID), lag-1);
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


