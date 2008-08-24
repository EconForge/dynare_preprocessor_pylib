/*
 * Copyright (C) 2006-2008 Dynare Team
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

#include "ModFile.hh"

ModFile::ModFile() : expressions_tree(symbol_table, num_constants),
                     model_tree(symbol_table, num_constants),
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

  // Allow empty model only when doing a standalone BVAR estimation
  if (model_tree.equation_number() == 0
      && (mod_file_struct.check_present
          || mod_file_struct.simul_present
          || mod_file_struct.stoch_simul_or_similar_present))
    {
      cerr << "Error: you must declare at least one model equation!" << endl;
      exit(-1);
    }

  if (mod_file_struct.simul_present
      && mod_file_struct.stoch_simul_or_similar_present)
    {
      cerr << "Error: a mod file cannot contain both a simul command and one of {stoch_simul, estimation, osr, ramsey_policy}" << endl;
      exit(-1);
    }
}

void
ModFile::computingPass()
{
  // Mod file may have no equation (for example in a standalone BVAR estimation)
  if (model_tree.equation_number() > 0)
    {
      // Set things to compute
      if (mod_file_struct.simul_present)
        model_tree.computeJacobian = true;
      else
        {
          if (mod_file_struct.order_option < 1 || mod_file_struct.order_option > 3)
            {
              cerr << "Incorrect order option..." << endl;
              exit(-1);
            }
          model_tree.computeJacobianExo = true;
          if (mod_file_struct.order_option >= 2)
            model_tree.computeHessian = true;
          if (mod_file_struct.order_option == 3)
            model_tree.computeThirdDerivatives = true;
        }

      model_tree.computingPass(global_eval_context);
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
          cerr << "Error: Can't open file " << fname
               << " for writing" << endl;
          exit(-1);
        }
    }
  else
    {
      cerr << "Error: Missing file name" << endl;
      exit(-1);
    }

  mOutputFile << "%" << endl
              << "% Status : main Dynare file " << endl
              << "%" << endl
              << "% Warning : this file is generated automatically by Dynare" << endl
              << "%           from model file (.mod)" << endl << endl;

  if (clear_all)
    mOutputFile << "clear all" << endl;
  mOutputFile << "tic;" << endl;
  mOutputFile << "global M_ oo_ exedet_ exdet_ recur_ recurs_ " << endl;
  mOutputFile << "global options_ endval_" << endl;
  mOutputFile << "global ys0_ recurs0_ ex0_ ct_" << endl;
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


  if (model_tree.equation_number() > 0)
    {
      if (model_tree.mode == eDLLMode)
        {
          mOutputFile << "if exist('" << basename << "_static.c')" << endl;
          mOutputFile << "   clear " << basename << "_static" << endl;
          mOutputFile << "   mex -O " << basename << "_static.c" << endl;
          mOutputFile << "end" << endl;
          mOutputFile << "if exist('" << basename << "_dynamic.c')" << endl;
          mOutputFile << "   clear " << basename << "_dynamic" << endl;
          mOutputFile << "   mex -O " << basename << "_dynamic.c" << endl;
          mOutputFile << "end" << endl;
        }
      else
        {
          mOutputFile << "erase_compiled_function('" + basename +"_static');" << endl;
          mOutputFile << "erase_compiled_function('" + basename +"_dynamic');" << endl;
        }
    }

  cout << "Processing outputs ..." << endl;

  symbol_table.writeOutput(mOutputFile);

  if (linear == 1)
    mOutputFile << "options_.linear = 1;" << endl;

  if (model_tree.equation_number() > 0)
    {
      model_tree.writeOutput(mOutputFile);
      model_tree.writeStaticFile(basename);
      model_tree.writeDynamicFile(basename);
    }

  // Print statements
  for(vector<Statement *>::const_iterator it = statements.begin();
      it != statements.end(); it++)
    (*it)->writeOutput(mOutputFile, basename);

  mOutputFile << "save('" << basename << "_results.mat', 'oo_');" << endl;
  mOutputFile << "diary off" << endl;

  mOutputFile << endl << "disp(['Total computing time : ' sec2hms(toc) ]);" << endl;
  mOutputFile.close();
}
