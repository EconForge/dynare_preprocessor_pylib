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
                     linear(false)
{
}

ModFile::~ModFile()
{
  for(vector<Statement *>::iterator it = statements.begin();
      it != statements.end(); it++)
    delete (*it);
}

void
ModFile::evalAllExpressions()
{
  cout << "Evaluating expressions...";

  // Loop over all statements, and fill global eval context if relevant
  for(vector<Statement *>::const_iterator it = statements.begin(); it != statements.end(); it++)
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
  for(int id = 0; id <= symbol_table.maxID(); id++)
    {
      SymbolType type = symbol_table.getType(id);
      if ((type == eEndogenous || type == eExogenous || type == eExogenousDet
          || type == eParameter || type == eModelLocalVariable)
          && global_eval_context.find(id) == global_eval_context.end())
        {
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
  for(vector<Statement *>::iterator it = statements.begin();
      it != statements.end(); it++)
    (*it)->checkPass(mod_file_struct);

  // If order option has not been set, default to 2
  if (!mod_file_struct.order_option)
    mod_file_struct.order_option = 2;

  bool stochastic_statement_present = mod_file_struct.stoch_simul_present
    || mod_file_struct.estimation_present
    || mod_file_struct.forecast_present
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
      cerr << "ERROR: A .mod file cannot contain both a simul command and one of {stoch_simul, estimation, forecast, osr, ramsey_policy}" << endl;
      exit(EXIT_FAILURE);
    }

  // Freeze the symbol table
  symbol_table.freeze();

  /*
    Enforce the same number of equations and endogenous, except in two cases:
    - ramsey_policy is used
    - a BVAR command is used and there is no equation (standalone BVAR estimation)
  */
  if (!mod_file_struct.ramsey_policy_present
      && !((mod_file_struct.bvar_density_present || mod_file_struct.bvar_forecast_present)
           && dynamic_model.equation_number() == 0)
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
  if (dynamic_model.equation_number() > 0)
    {
      // Compute static model and its derivatives
      dynamic_model.toStatic(static_model);
      static_model.computingPass(false, no_tmp_terms);

      // Set things to compute for dynamic model
      
      if (mod_file_struct.simul_present)
        dynamic_model.computingPass(false, false, false, false, global_eval_context, no_tmp_terms);
      else
        {
          if (mod_file_struct.order_option < 1 || mod_file_struct.order_option > 3)
            {
              cerr << "ERROR: Incorrect order option..." << endl;
              exit(EXIT_FAILURE);
            }
          bool hessian = mod_file_struct.order_option >= 2;
          bool thirdDerivatives = mod_file_struct.order_option == 3;
          bool paramsDerivatives = mod_file_struct.identification_present;
          dynamic_model.computingPass(true, hessian, thirdDerivatives, paramsDerivatives, global_eval_context, no_tmp_terms);
        }
    }

  for(vector<Statement *>::iterator it = statements.begin();
      it != statements.end(); it++)
    (*it)->computingPass();
}

void
ModFile::writeOutputFiles(const string &basename, bool clear_all) const
{
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
              << "% Status : main Dynare file " << endl
              << "%" << endl
              << "% Warning : this file is generated automatically by Dynare" << endl
              << "%           from model file (.mod)" << endl << endl;

  if (clear_all)
    mOutputFile << "clear all" << endl;
  mOutputFile << "tic;" << endl;
  mOutputFile << "global M_ oo_ options_" << endl;
  mOutputFile << "global ys0_ ex0_ ct_" << endl;
  mOutputFile << "options_ = [];" << endl;
  mOutputFile << "M_.fname = '" << basename << "';" << endl;
  mOutputFile << "%" << endl;
  mOutputFile << "% Some global variables initialization" << endl;
  mOutputFile << "%" << endl;
  mOutputFile << "global_initialization;" << endl;
  mOutputFile << "diary off;" << endl
              << "warning_old_state = warning;" << endl
              << "warning off;" << endl
              << "delete " << basename << ".log;" << endl
              << "warning warning_old_state" << endl;
  mOutputFile << "logname_ = '" << basename << ".log';" << endl;
  mOutputFile << "diary " << basename << ".log" << endl;
  mOutputFile << "options_.model_mode = " << dynamic_model.mode << ";\n";
  if (dynamic_model.mode == eSparseMode || dynamic_model.mode == eSparseDLLMode)
    {
      mOutputFile << "addpath " << basename << ";\n";
      mOutputFile << "delete('" << basename << "_dynamic.m');\n";
    }


  if (dynamic_model.equation_number() > 0)
    {
      if (dynamic_model.mode == eDLLMode || dynamic_model.mode == eSparseDLLMode)
        {
          mOutputFile << "if exist('" << basename << "_static.c')" << endl;
          mOutputFile << "   clear " << basename << "_static" << endl;
          mOutputFile << "   mex -O " << basename << "_static.c" << endl;
          mOutputFile << "end" << endl;
          mOutputFile << "if exist('" << basename << "_dynamic.c')" << endl;
          mOutputFile << "   clear " << basename << "_dynamic" << endl;
          mOutputFile << "   mex -O  LDFLAGS='-pthread -shared -Wl,--no-undefined' " << basename << "_dynamic.c" << endl;
          mOutputFile << "end" << endl;
        }
      else
        {
          mOutputFile << "erase_compiled_function('" + basename +"_static');" << endl;
          mOutputFile << "erase_compiled_function('" + basename +"_dynamic');" << endl;
        }
    }
  cout << "Processing outputs ...";

  symbol_table.writeOutput(mOutputFile);

  if (linear == 1)
    mOutputFile << "options_.linear = 1;" << endl;

  if (dynamic_model.equation_number() > 0)
    {
      dynamic_model.writeOutput(mOutputFile);
      static_model.writeStaticFile(basename);
      dynamic_model.writeDynamicFile(basename);
      dynamic_model.writeParamsDerivativesFile(basename);
    }

  // Print statements
  for(vector<Statement *>::const_iterator it = statements.begin();
      it != statements.end(); it++)
    (*it)->writeOutput(mOutputFile, basename);

  mOutputFile << "save('" << basename << "_results.mat', 'oo_', 'M_', 'options_');" << endl;
  mOutputFile << "diary off" << endl;

  if (dynamic_model.mode == eSparseMode || dynamic_model.mode == eSparseDLLMode)
    mOutputFile << "rmpath " << basename << ";\n";

  mOutputFile << endl << "disp(['Total computing time : ' dynsec2hms(toc) ]);" << endl;
  mOutputFile.close();
  cout << "done\n";
}
