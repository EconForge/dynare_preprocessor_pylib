/*
 * Copyright (C) 2006-2009 Dynare Team
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
#include <iostream>
#include <fstream>
#include <typeinfo>
#include "ModFile.hh"

ModFile::ModFile() : expressions_tree(symbol_table, num_constants),
                     static_model(symbol_table, num_constants),
                     dynamic_model(symbol_table, num_constants),
                     linear(false), block(false), byte_code(false),
                     use_dll(false), no_static(false)
{
}

ModFile::~ModFile()
{
  for (vector<Statement *>::iterator it = statements.begin();
       it != statements.end(); it++)
    delete (*it);
}

void
ModFile::evalAllExpressions(bool warn_uninit)
{
  cout << "Evaluating expressions...";

  // Loop over all statements, and fill global eval context if relevant
  for (vector<Statement *>::const_iterator it = statements.begin(); it != statements.end(); it++)
    {
      InitParamStatement *ips = dynamic_cast<InitParamStatement *>(*it);
      if (ips)
        ips->fillEvalContext(global_eval_context);

      InitOrEndValStatement *ies = dynamic_cast<InitOrEndValStatement *>(*it);
      if (ies)
        ies->fillEvalContext(global_eval_context);

      LoadParamsAndSteadyStateStatement *lpass = dynamic_cast<LoadParamsAndSteadyStateStatement *>(*it);
      if (lpass)
        lpass->fillEvalContext(global_eval_context);
    }

  // Evaluate model local variables
  dynamic_model.fillEvalContext(global_eval_context);

  cout << "done" << endl;

  // Check if some symbols are not initialized, and give them a zero value then
  for (int id = 0; id <= symbol_table.maxID(); id++)
    {
      SymbolType type = symbol_table.getType(id);
      if ((type == eEndogenous || type == eExogenous || type == eExogenousDet
           || type == eParameter || type == eModelLocalVariable)
          && global_eval_context.find(id) == global_eval_context.end())
        {
          if (warn_uninit)
            cerr << "WARNING: can't find a numeric initial value for " << symbol_table.getName(id) << ", using zero" << endl;
          global_eval_context[id] = 0;
        }
    }
}

void
ModFile::addStatement(Statement *st)
{
  statements.push_back(st);
}

void
ModFile::checkPass()
{
  for (vector<Statement *>::iterator it = statements.begin();
       it != statements.end(); it++)
    (*it)->checkPass(mod_file_struct);

  // If order option has not been set, default to 2
  if (!mod_file_struct.order_option)
    mod_file_struct.order_option = 2;

  bool stochastic_statement_present = mod_file_struct.stoch_simul_present
    || mod_file_struct.estimation_present
    || mod_file_struct.osr_present
    || mod_file_struct.ramsey_policy_present;

  // Allow empty model only when doing a standalone BVAR estimation
  if (dynamic_model.equation_number() == 0
      && (mod_file_struct.check_present
          || mod_file_struct.simul_present
          || stochastic_statement_present))
    {
      cerr << "ERROR: At least one model equation must be declared!" << endl;
      exit(EXIT_FAILURE);
    }

  if (mod_file_struct.simul_present && stochastic_statement_present)
    {
      cerr << "ERROR: A .mod file cannot contain both a simul command and one of {stoch_simul, estimation, osr, ramsey_policy}" << endl;
      exit(EXIT_FAILURE);
    }

  if (mod_file_struct.k_order_solver && !use_dll)
    {
      cerr << "ERROR: When using option 'k_order_solver' (which is implicit if order >= 3), you must specify option 'use_dll' on the 'model' block" << endl;
      exit(EXIT_FAILURE);
    }

  if (use_dll && (block || byte_code))
    {
      cerr << "ERROR: In 'model' block, 'use_dll' option is not compatible with 'block' or 'bytecode'" << endl;
      exit(EXIT_FAILURE);
    }

  if ( (stochastic_statement_present || mod_file_struct.check_present || mod_file_struct.steady_present) && no_static)
    {
      cerr << "no_static option is incompatible with stochastic simulation, estimation, optimal policy, steady or check command" << endl;
      exit(EXIT_FAILURE);
    }
}

void
ModFile::transformPass()
{
  if (symbol_table.predeterminedNbr() > 0)
    dynamic_model.transformPredeterminedVariables();

  // Create auxiliary vars for Expectation operator
  dynamic_model.substituteExpectation(mod_file_struct.partial_information);

  if (mod_file_struct.stoch_simul_present
      || mod_file_struct.estimation_present
      || mod_file_struct.osr_present
      || mod_file_struct.ramsey_policy_present)
    {
      // In stochastic models, create auxiliary vars for leads and lags greater than 2
      dynamic_model.substituteEndoLeadGreaterThanTwo();
      dynamic_model.substituteExoLead();
      dynamic_model.substituteEndoLagGreaterThanTwo();
      dynamic_model.substituteExoLag();
    }

  // Freeze the symbol table
  symbol_table.freeze();

  /*
    Enforce the same number of equations and endogenous, except in two cases:
    - ramsey_policy is used
    - a BVAR command is used and there is no equation (standalone BVAR estimation)
  */
  if (!mod_file_struct.ramsey_policy_present
      && !(mod_file_struct.bvar_present && dynamic_model.equation_number() == 0)
      && (dynamic_model.equation_number() != symbol_table.endo_nbr()))
    {
      cerr << "ERROR: There are " << dynamic_model.equation_number() << " equations but " << symbol_table.endo_nbr() << " endogenous variables!" << endl;
      exit(EXIT_FAILURE);
    }

  cout << "Found " << dynamic_model.equation_number() << " equation(s)." << endl;
}

void
ModFile::computingPass(bool no_tmp_terms)
{
  // Mod file may have no equation (for example in a standalone BVAR estimation)
  bool dynamic_model_needed = mod_file_struct.simul_present || mod_file_struct.check_present || mod_file_struct.stoch_simul_present
    || mod_file_struct.estimation_present || mod_file_struct.osr_present
    || mod_file_struct.ramsey_policy_present || mod_file_struct.identification_present;
  if (dynamic_model.equation_number() > 0)
    {
      // Compute static model and its derivatives
      dynamic_model.toStatic(static_model);
      if(!no_static)
        {
          static_model.initializeVariablesAndEquations();
          static_model.computingPass(global_eval_context, no_tmp_terms, false, block, byte_code);
        }
      // Set things to compute for dynamic model
      if (dynamic_model_needed)
        {
          if (mod_file_struct.simul_present)
            {
              dynamic_model.initializeVariablesAndEquations();
              dynamic_model.computingPass(false, false, false, false, global_eval_context, no_tmp_terms, block, use_dll, byte_code);
            }
          else
            {
              if (mod_file_struct.order_option < 1 || mod_file_struct.order_option > 3)
                {
                  cerr << "ERROR: Incorrect order option..." << endl;
                  exit(EXIT_FAILURE);
                }
              bool hessian = mod_file_struct.order_option >= 2 || mod_file_struct.identification_present;
              bool thirdDerivatives = mod_file_struct.order_option == 3;
              bool paramsDerivatives = mod_file_struct.identification_present;
              dynamic_model.computingPass(true, hessian, thirdDerivatives, paramsDerivatives, global_eval_context, no_tmp_terms, false, use_dll, byte_code);
            }
        }
      else
        dynamic_model.computingPass(true, true, false, false, global_eval_context, no_tmp_terms, false, false, byte_code);
    }

  for (vector<Statement *>::iterator it = statements.begin();
       it != statements.end(); it++)
    (*it)->computingPass();
}

void
ModFile::writeOutputFiles(const string &basename, bool clear_all
#if defined(_WIN32) || defined(__CYGWIN32__)
                          , bool cygwin, bool msvc
#endif
                          ) const
{
  ofstream mOutputFile;
  bool dynamic_model_needed = mod_file_struct.simul_present || mod_file_struct.check_present || mod_file_struct.stoch_simul_present
    || mod_file_struct.estimation_present || mod_file_struct.osr_present
    || mod_file_struct.ramsey_policy_present || mod_file_struct.identification_present;

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
              << "% Status : main Dynare file " << endl
              << "%" << endl
              << "% Warning : this file is generated automatically by Dynare" << endl
              << "%           from model file (.mod)" << endl << endl;

  if (clear_all)
    mOutputFile << "clear all" << endl;

  mOutputFile << "tic;" << endl
              << "global M_ oo_ options_ ys0_ ex0_" << endl
              << "options_ = [];" << endl
              << "M_.fname = '" << basename << "';" << endl
              << "%" << endl
              << "% Some global variables initialization" << endl
              << "%" << endl
              << "global_initialization;" << endl
              << "diary off;" << endl
              << "warning_old_state = warning;" << endl
              << "warning off;" << endl
              << "delete " << basename << ".log;" << endl
              << "warning warning_old_state" << endl
              << "logname_ = '" << basename << ".log';" << endl
              << "diary " << basename << ".log" << endl;

  cout << "Processing outputs ...";

  symbol_table.writeOutput(mOutputFile);

  if (linear == 1)
    mOutputFile << "options_.linear = 1;" << endl;

  mOutputFile << "options_.block=" << block << ";" << endl
              << "options_.bytecode=" << byte_code << ";" << endl
              << "options_.use_dll=" << use_dll << ";" << endl;

  if (byte_code)
    mOutputFile << "if exist('bytecode') ~= 3" << endl
                << "  error('DYNARE: Can''t find bytecode DLL. Please compile it or remove the ''bytecode'' option.')" << endl
                << "end" << endl;

  // Erase possible remnants of previous runs
  if (block || byte_code)
    mOutputFile << "delete('" << basename << "_dynamic.m');" << endl;

  if (byte_code)
    mOutputFile << "delete('" << basename << "_static.m');" << endl;

  if (!use_dll)
    mOutputFile << "erase_compiled_function('" + basename + "_dynamic');" << endl;

  // Compile the dynamic MEX file for use_dll option
  if (use_dll)
    {
      mOutputFile << "if ~exist('OCTAVE_VERSION')" << endl;
      // Some mex commands are enclosed in an eval(), because otherwise it will make Octave fail
#if defined(_WIN32) || defined(__CYGWIN32__)
      if (msvc)
        mOutputFile << "    eval('mex -O LINKFLAGS=\"$LINKFLAGS /export:Dynamic\" " << basename << "_dynamic.c')" << endl;                                                                                                                                                                                                                                                       // MATLAB/Windows + Microsoft Visual C++
      else if (cygwin)
        mOutputFile << "    eval('mex -O PRELINK_CMDS1=\"echo EXPORTS > mex.def & echo mexFunction >> mex.def & echo Dynamic >> mex.def\" " << basename << "_dynamic.c')" << endl;                                                                                                                                                                                                                                                                                                                                                                        // MATLAB/Windows + Cygwin g++
      else
        {
          cerr << "ERROR: When using the USE_DLL option, you must give either 'cygwin' or 'msvc' option to the 'dynare' command" << endl;
          exit(EXIT_FAILURE);
        }
#else
# ifdef __linux__
      mOutputFile << "    eval('mex -O LDFLAGS=''-pthread -shared -Wl,--no-undefined'' " << basename << "_dynamic.c')" << endl; // MATLAB/Linux
# else // MacOS
      mOutputFile << "    eval('mex -O LDFLAGS=''-Wl,-twolevel_namespace -undefined error -arch \\$ARCHS -Wl,-syslibroot,\\$SDKROOT -mmacosx-version-min=\\$MACOSX_DEPLOYMENT_TARGET -bundle'' " << basename << "_dynamic.c')" << endl; // MATLAB/MacOS
# endif
#endif
      mOutputFile << "else" << endl // Octave
                  << "    if ~octave_ver_less_than('3.2.0')" << endl // Workaround for bug in Octave >= 3.2, see http://bugs.debian.org/cgi-bin/bugreport.cgi?bug=550823
                  << "        sleep(2)" << endl
                  << "    end" << endl
                  << "    mex " << basename << "_dynamic.c" << endl
                  << "end" << endl;
    }

  // Add path for block option with M-files
  if (block && !byte_code)
    mOutputFile << "addpath " << basename << ";" << endl;

  if (dynamic_model.equation_number() > 0)
    {
      if (dynamic_model_needed)
        dynamic_model.writeOutput(mOutputFile, basename, block, byte_code, use_dll);
      else
        dynamic_model.writeOutput(mOutputFile, basename, false, false, false);
      if (!byte_code && !no_static)
        static_model.writeOutput(mOutputFile, block);
    }

  // Print statements
  for (vector<Statement *>::const_iterator it = statements.begin();
       it != statements.end(); it++)
    {
      (*it)->writeOutput(mOutputFile, basename);

      // Special treatment for initval block: insert initial values for the auxiliary variables
      InitValStatement *ivs = dynamic_cast<InitValStatement *>(*it);
      if (ivs != NULL)
        {
          static_model.writeAuxVarInitval(mOutputFile);
          ivs->writeOutputPostInit(mOutputFile);
        }

      // Special treatment for load params and steady state statement: insert initial values for the auxiliary variables
      LoadParamsAndSteadyStateStatement *lpass = dynamic_cast<LoadParamsAndSteadyStateStatement *>(*it);
      if (lpass && !no_static)
        static_model.writeAuxVarInitval(mOutputFile);
    }

  // Remove path for block option with M-files
  if (block && !byte_code)
    mOutputFile << "rmpath " << basename << ";" << endl;

  mOutputFile << "save('" << basename << "_results.mat', 'oo_', 'M_', 'options_');" << endl
              << "diary off" << endl
              << endl << "disp(['Total computing time : ' dynsec2hms(toc) ]);" << endl;

  mOutputFile.close();

  // Create static and dynamic files
  if (dynamic_model.equation_number() > 0)
    {
      if (!no_static)
        static_model.writeStaticFile(basename, block, byte_code);

      if (dynamic_model_needed)
        {
          dynamic_model.writeDynamicFile(basename, block, byte_code, use_dll);
          dynamic_model.writeParamsDerivativesFile(basename);
        }
      else
        {
          dynamic_model.writeDynamicFile(basename, false, false, false);
          dynamic_model.writeParamsDerivativesFile(basename);
        }
    }

  cout << "done" << endl;
}
