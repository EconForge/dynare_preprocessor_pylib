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

#input "ExternalFiles.hh"

ExternalFiles::writeModelCCconst string &basename, bool cuda)
{
  ofstream mOutputFile;

  if (basename.size())
    {
      string fname(basename);
      fname += ".cc";
      cOutputFile.open(fname.c_str(), ios::out | ios::binary);
      if (!cOutputFile.is_open())
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

  cOutputFile << "%" << endl
              << "% Status : dynare_model class initializing function " << endl
              << "%" << endl
              << "% Warning : this file is generated automatically by Dynare" << endl
              << "%           from model file " << baseline << "(.mod)" << endl << endl;

  cOutputFile << "#include \"dynare_model.hh\"" << endl << endl;

  cOutputFile << "dynare_model::dynare_model(void)" << endl
              << "{" << endl
              << "    model_name = \""<< basename << "\";" << endl;

  symbol_table.writeOutputCC(hOutputFile);

  // Initialize M_.Sigma_e and M_.H
  mOutputFile << "M_.Sigma_e = zeros(" << symbol_table.exo_nbr() << ", "
              << symbol_table.exo_nbr() << ");" << endl;

  if (mod_file_struct.calibrated_measurement_errors)
    mOutputFile << "M_.H = zeros(" << symbol_table.observedVariablesNbr() << ", "
                << symbol_table.observedVariablesNbr() << ");" << endl;
  else
    mOutputFile << "M_.H = 0;" << endl;

  if (linear == 1)
    mOutputFile << "options_.linear = 1;" << endl;

  mOutputFile << "options_.block=" << block << ";" << endl
              << "options_.bytecode=" << byte_code << ";" << endl
              << "options_.use_dll=" << use_dll << ";" << endl;

  if (parallel_local_files.size() > 0)
    {
      mOutputFile << "options_.parallel_info.local_files = {" << endl;
      for (size_t i = 0; i < parallel_local_files.size(); i++)
        {
          size_t j = parallel_local_files[i].find_last_of("/\\");
          if (j == string::npos)
            mOutputFile << "'', '" << parallel_local_files[i] << "';" << endl;
          else
            mOutputFile << "'" << parallel_local_files[i].substr(0, j+1) << "', '"
                        << parallel_local_files[i].substr(j+1, string::npos) << "';" << endl;
        }
      mOutputFile << "};" << endl;
    }

  config_file.writeCluster(mOutputFile);

  if (byte_code)
    mOutputFile << "if exist('bytecode') ~= 3" << endl
                << "  error('DYNARE: Can''t find bytecode DLL. Please compile it or remove the ''bytecode'' option.')" << endl
                << "end" << endl;

  // Erase possible remnants of previous runs
  string dynfile = basename + "_dynamic.m";
  unlink(dynfile.c_str());

  string statfile = basename + "_static.m";
  unlink(statfile.c_str());

  string steadystatefile = basename + "_steadystate2.m";
  unlink(steadystatefile.c_str());

  if (!use_dll)
    {
      mOutputFile << "erase_compiled_function('" + basename + "_static');" << endl;
      mOutputFile << "erase_compiled_function('" + basename + "_dynamic');" << endl;
    }

#if defined(_WIN32) || defined(__CYGWIN32__)
  // If using USE_DLL with MSVC, check that the user didn't use a function not supported by MSVC (because MSVC doesn't comply with C99 standard)
  if (use_dll && msvc)
    {
      if (dynamic_model.isUnaryOpUsed(oAcosh))
        {
          cerr << "ERROR: acosh() function is not supported with USE_DLL option and MSVC compiler; use Cygwin compiler instead." << endl;
          exit(EXIT_FAILURE);
        }
      if (dynamic_model.isUnaryOpUsed(oAsinh))
        {
          cerr << "ERROR: asinh() function is not supported with USE_DLL option and MSVC compiler; use Cygwin compiler instead." << endl;
          exit(EXIT_FAILURE);
        }
      if (dynamic_model.isUnaryOpUsed(oAtanh))
        {
          cerr << "ERROR: atanh() function is not supported with USE_DLL option and MSVC compiler; use Cygwin compiler instead." << endl;
          exit(EXIT_FAILURE);
        }
      if (dynamic_model.isTrinaryOpUsed(oNormcdf))
        {
          cerr << "ERROR: normcdf() function is not supported with USE_DLL option and MSVC compiler; use Cygwin compiler instead." << endl;
          exit(EXIT_FAILURE);
        }
    }
#endif

  // Compile the dynamic MEX file for use_dll option
  if (use_dll)
    {
      mOutputFile << "if ~exist('OCTAVE_VERSION')" << endl;
      // Some mex commands are enclosed in an eval(), because otherwise it will make Octave fail
#if defined(_WIN32) || defined(__CYGWIN32__)
      if (msvc)
        // MATLAB/Windows + Microsoft Visual C++
        mOutputFile << "    eval('mex -O LINKFLAGS=\"$LINKFLAGS /export:Dynamic\" " << basename << "_dynamic.c " << basename << "_dynamic_mex.c')" << endl
                    << "    eval('mex -O LINKFLAGS=\"$LINKFLAGS /export:Static\" " << basename << "_static.c "<< basename << "_static_mex.c')" << endl;
      else if (cygwin)
        // MATLAB/Windows + Cygwin g++
        mOutputFile << "    eval('mex -O PRELINK_CMDS1=\"echo EXPORTS > mex.def & echo mexFunction >> mex.def & echo Dynamic >> mex.def\" " << basename << "_dynamic.c " << basename << "_dynamic_mex.c')" << endl
                    << "    eval('mex -O PRELINK_CMDS1=\"echo EXPORTS > mex.def & echo mexFunction >> mex.def & echo Static >> mex.def\" " << basename << "_static.c "<< basename << "_static_mex.c')" << endl;
      else
        mOutputFile << "    error('When using the USE_DLL option, you must give either ''cygwin'' or ''msvc'' option to the ''dynare'' command')" << endl;
#else
# ifdef __linux__
      // MATLAB/Linux
      mOutputFile << "    eval('mex -O LDFLAGS=''-pthread -shared -Wl,--no-undefined'' " << basename << "_dynamic.c " << basename << "_dynamic_mex.c')" << endl
                  << "    eval('mex -O LDFLAGS=''-pthread -shared -Wl,--no-undefined'' " << basename << "_static.c "<< basename << "_static_mex.c')" << endl;
# else // MacOS
      // MATLAB/MacOS
      mOutputFile << "    eval('mex -O LDFLAGS=''-Wl,-twolevel_namespace -undefined error -arch \\$ARCHS -Wl,-syslibroot,\\$SDKROOT -mmacosx-version-min=\\$MACOSX_DEPLOYMENT_TARGET -bundle'' "
                  << basename << "_dynamic.c " << basename << "_dynamic_mex.c')" << endl
                  << "    eval('mex -O LDFLAGS=''-Wl,-twolevel_namespace -undefined error -arch \\$ARCHS -Wl,-syslibroot,\\$SDKROOT -mmacosx-version-min=\\$MACOSX_DEPLOYMENT_TARGET -bundle'' "
                  << basename << "_static.c " << basename << "_static_mex.c')" << endl;
# endif
#endif
      mOutputFile << "else" << endl // Octave
                  << "    mex " << basename << "_dynamic.c " << basename << "_dynamic_mex.c" << endl
                  << "    mex " << basename << "_static.c " << basename << "_static_mex.c" << endl
                  << "end" << endl;
    }

  // Add path for block option with M-files
  if (block && !byte_code)
    mOutputFile << "addpath " << basename << ";" << endl;

  if (mod_file_struct.ramsey_policy_present)
    mOutputFile << "M_.orig_eq_nbr = " << ramsey_policy_orig_eqn_nbr << ";" << endl;

  if (dynamic_model.equation_number() > 0)
    {
      dynamic_model.writeOutput(mOutputFile, basename, block, byte_code, use_dll, mod_file_struct.order_option, mod_file_struct.estimation_present);
      if (!no_static)
        static_model.writeOutput(mOutputFile, block);
    }

  // Print statements
  for (vector<Statement *>::const_iterator it = statements.begin();
       it != statements.end(); it++)
    {
      (*it)->writeOutput(mOutputFile, basename);

      /* Special treatment for initval block: insert initial values for the
         auxiliary variables and initialize exo det */
      InitValStatement *ivs = dynamic_cast<InitValStatement *>(*it);
      if (ivs != NULL)
        {
          static_model.writeAuxVarInitval(mOutputFile, oMatlabOutsideModel);
          ivs->writeOutputPostInit(mOutputFile);
        }

      // Special treatment for endval block: insert initial values for the auxiliary variables
      EndValStatement *evs = dynamic_cast<EndValStatement *>(*it);
      if (evs != NULL)
        static_model.writeAuxVarInitval(mOutputFile, oMatlabOutsideModel);

      // Special treatment for load params and steady state statement: insert initial values for the auxiliary variables
      LoadParamsAndSteadyStateStatement *lpass = dynamic_cast<LoadParamsAndSteadyStateStatement *>(*it);
      if (lpass && !no_static)
        static_model.writeAuxVarInitval(mOutputFile, oMatlabOutsideModel);
    }

  // Remove path for block option with M-files
  if (block && !byte_code)
    mOutputFile << "rmpath " << basename << ";" << endl;

  mOutputFile << "save('" << basename << "_results.mat', 'oo_', 'M_', 'options_');" << endl;

  config_file.writeEndParallel(mOutputFile);

  mOutputFile << endl << endl
	      << "disp(['Total computing time : ' dynsec2hms(toc) ]);" << endl;

  if (!no_warn)
    {
      if (warnings.countWarnings() > 0)
        mOutputFile << "disp('Note: " << warnings.countWarnings() << " warning(s) encountered in the preprocessor')" << endl;

      mOutputFile << "if ~isempty(lastwarn)" << endl
                  << "  disp('Note: warning(s) encountered in MATLAB/Octave code')" << endl
                  << "end" << endl;
    }
  

  if (!no_log)
    mOutputFile << "diary off" << endl;

  mOutputFile.close();

  // Create static and dynamic files
  if (dynamic_model.equation_number() > 0)
    {
      if (!no_static)
        {
          static_model.writeStaticFile(basename, block, byte_code, use_dll);
          static_model.writeParamsDerivativesFile(basename);
        }

      dynamic_model.writeDynamicFile(basename, block, byte_code, use_dll, mod_file_struct.order_option);
      dynamic_model.writeParamsDerivativesFile(basename);
    }

  // Create steady state file
  steady_state_model.writeSteadyStateFile(basename, mod_file_struct.ramsey_policy_present);

  cout << "done" << endl;
}
