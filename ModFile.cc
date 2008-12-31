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
#include <typeinfo>
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
ModFile::evalAllExpressions()
{
  //Evaluate Parameters

  InitParamStatement *it;
  vector< vector<double> >::iterator it2;
  ostringstream constant;
  NodeID tmp_id;
  CollectStruct collect_struct;
  int j=0, k;
  if(mod_file_struct.load_params_and_steady_state_present)
    {
      cout << "Reading " << mod_file_struct.load_params_and_steady_state_filename << " ...";
      matlab_file.MatFileRead(mod_file_struct.load_params_and_steady_state_filename);
      string sname="stored_values";
      bool tmp_b=matlab_file.Collect(sname, collect_struct);
      matlab_file.Delete();
      if(!tmp_b)
        {
          cout << "The structure " << sname << " is not found in " << mod_file_struct.load_params_and_steady_state_filename << "\n";
        }
      cout << "done\n";
    }
  cout << "Evaluating expressions ...";
  for(vector<Statement *>::const_iterator it1=statements.begin();it1!=statements.end(); it1++)
    {
      it=dynamic_cast<InitParamStatement *>(*it1);
      if(it)
        {
          try
            {
              const NodeID expression = it->get_expression();
              double val = expression->eval(global_eval_context);
              int symb_id = symbol_table.getID(it->get_name());
              global_eval_context[make_pair(symb_id, eParameter)] = val;
              j++;
            }
          catch(ExprNode::EvalException &e)
           {
             cout << "error in evaluation of param\n";
           }
        }
    }
  if(mod_file_struct.load_params_and_steady_state_present && j!=symbol_table.parameter_nbr)
    {
      //Reading a Mat-File
      for(k=0;k <symbol_table.parameter_nbr; k++)
        {
          if(global_eval_context.find(make_pair(k, eParameter))==global_eval_context.end())
            {
              map<string,vector<double> >::iterator it2=collect_struct.variable_double_name.find(symbol_table.getNameByID(eParameter, k));
              if(it2!=collect_struct.variable_double_name.end())
                {
                  j++;
                  vector<double>::iterator it=it2->second.begin();
                  global_eval_context[make_pair(k, eParameter)]=*it;
                }
            }
        }
    }
  if (j!=symbol_table.parameter_nbr)
    {
      cout << "Warning: Uninitialized parameters: \n";
      for(j=0;j <symbol_table.parameter_nbr; j++)
        {
          if(global_eval_context.find(make_pair(j, eParameter))==global_eval_context.end())
            cout << " " << symbol_table.getNameByID(eParameter, j) << "\n";
        }
    }
  //Evaluate variables
  for(InitOrEndValStatement::init_values_type::const_iterator it=init_values.begin(); it!=init_values.end(); it++)
    {
      try
        {
          const string &name = it->first;
          const NodeID expression = it->second;
          SymbolType type = symbol_table.getType(name);
          double val = expression->eval(global_eval_context);
          int symb_id = symbol_table.getID(name);
          global_eval_context[make_pair(symb_id, type)] = val;
        }
      catch(ExprNode::EvalException &e)
        {
          cout << "error in evaluation of variable\n";
        }
    }
  if(mod_file_struct.load_params_and_steady_state_present && init_values.size()<symbol_table.endo_nbr+symbol_table.exo_nbr+symbol_table.exo_det_nbr)
    {
      for(j=0;j <symbol_table.endo_nbr; j++)
        {
          if(global_eval_context.find(make_pair(j, eEndogenous))==global_eval_context.end())
            {
              //it2=mat_file.variable.find(symbol_table.getNameByID(eEndogenous, j));
              map<string,vector<double> >::iterator it2=collect_struct.variable_double_name.find(symbol_table.getNameByID(eEndogenous, j));
              if(it2!=collect_struct.variable_double_name.end())
                {
                  vector<double>::iterator it=it2->second.begin();
                  global_eval_context[make_pair(j, eEndogenous)]=*it;
                  constant.str("");
                  if(*it>=0)
                    {
                      constant << *it;
                      tmp_id=expressions_tree.AddNumConstant(constant.str());
                    }
                  else
                    {
                      constant << -*it;
                      tmp_id=expressions_tree.AddUMinus(expressions_tree.AddNumConstant(constant.str()));
                    }
                  init_values.push_back(make_pair(it2->first, tmp_id));
                }
            }
        }
      for(j=0;j <symbol_table.exo_nbr; j++)
        {
          if(global_eval_context.find(make_pair(j, eExogenous))==global_eval_context.end())
            {
              map<string,vector<double> >::iterator it2=collect_struct.variable_double_name.find(symbol_table.getNameByID(eExogenous, j));
              if(it2!=collect_struct.variable_double_name.end())
                {
                  vector<double>::iterator it=it2->second.begin();
                  global_eval_context[make_pair(j, eExogenous)]=*it;
                  constant.str("");
                  if(*it>=0)
                    {
                      constant << *it;
                      tmp_id=expressions_tree.AddNumConstant(constant.str());
                    }
                  else
                    {
                      constant << -*it;
                      tmp_id=expressions_tree.AddUMinus(expressions_tree.AddNumConstant(constant.str()));
                    }
                  init_values.push_back(make_pair(it2->first, tmp_id));
                }
            }
        }
      for(j=0;j <symbol_table.exo_det_nbr; j++)
        {
          if(global_eval_context.find(make_pair(j, eExogenous))==global_eval_context.end())
            {
              map<string,vector<double> >::iterator it2=collect_struct.variable_double_name.find(symbol_table.getNameByID(eExogenous, j));
              if(it2!=collect_struct.variable_double_name.end())
                {
                  vector<double>::iterator it=it2->second.begin();
                  global_eval_context[make_pair(j, eExogenous)]=*it;
                  constant.str("");
                  if(*it>=0)
                    {
                      constant << *it;
                      tmp_id=expressions_tree.AddNumConstant(constant.str());
                    }
                  else
                    {
                      constant << -*it;
                      tmp_id=expressions_tree.AddUMinus(expressions_tree.AddNumConstant(constant.str()));
                    }
                  init_values.push_back(make_pair(it2->first, tmp_id));
                }
            }
        }
    }
  if(init_values.size()<symbol_table.endo_nbr+symbol_table.exo_nbr+symbol_table.exo_det_nbr)
    {
      cout << "\nWarning: Uninitialized variable: \n";
      cout << "Endogenous\n";
      for(j=0;j <symbol_table.endo_nbr; j++)
        {
          if(global_eval_context.find(make_pair(j, eEndogenous))==global_eval_context.end())
            {
              cout << " " << symbol_table.getNameByID(eEndogenous, j) << "\n";
              global_eval_context[make_pair(j, eEndogenous)] = 0;
            }
        }
      cout << "Exogenous\n";
      for(j=0;j <symbol_table.exo_nbr; j++)
        {
          if(global_eval_context.find(make_pair(j, eExogenous))==global_eval_context.end())
            {
              cout << " " << symbol_table.getNameByID(eExogenous, j) << "\n";
              global_eval_context[make_pair(j, eExogenous)]=0;
            }
        }
      cout << "Deterministic exogenous\n";
      for(j=0;j <symbol_table.exo_det_nbr; j++)
        {
          if(global_eval_context.find(make_pair(j, eExogenousDet))==global_eval_context.end())
            {
              cout << " " << symbol_table.getNameByID(eExogenousDet, j) << "\n";
              global_eval_context[make_pair(j, eExogenousDet)]=0;
            }
        }
    }
  //Evaluate Local variables
  for(map<int, NodeID>::const_iterator it = model_tree.local_variables_table.begin(); it !=model_tree.local_variables_table.end(); it++)
    {
      try
        {
          const NodeID expression = it->second;
          double val = expression->eval(global_eval_context);
          //cout << it->first << "  " << symbol_table.getNameByID(eModelLocalVariable, it->first) << " = " << val << "\n";
          global_eval_context[make_pair(it->first, eModelLocalVariable)] = val;
        }
      catch(ExprNode::EvalException &e)
        {
          cout << "error in evaluation of pound\n";
        }
    }
  if(model_tree.local_variables_table.size()!=symbol_table.model_local_variable_nbr+symbol_table.modfile_local_variable_nbr)
    {
      cout << "Warning: Unitilialized pound: \n";
      cout << "Local variable in a model\n";
      for(j=0;j <symbol_table.model_local_variable_nbr; j++)
        {
          if(global_eval_context.find(make_pair(j, eModelLocalVariable))==global_eval_context.end())
            cout << " " << symbol_table.getNameByID(eModelLocalVariable, j) << "\n";
        }
      cout << "Local variable in a model file\n";
      for(j=0;j <symbol_table.modfile_local_variable_nbr; j++)
        {
          if(global_eval_context.find(make_pair(j, eModFileLocalVariable))==global_eval_context.end())
            cout << " " << symbol_table.getNameByID(eModFileLocalVariable, j) << "\n";
        }
    }
  cout << "done\n";
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
  if (model_tree.equation_number() == 0
      && (mod_file_struct.check_present
          || mod_file_struct.simul_present
          || stochastic_statement_present))
    {
      cerr << "ERROR: At least one model equation must be declared!" << endl;
      exit(EXIT_FAILURE);
    }

  if (mod_file_struct.simul_present && stochastic_statement_present && model_tree.mode==0)
    {
      cerr << "ERROR: A .mod file cannot contain both a simul command and one of {stoch_simul, estimation, forecast, osr, ramsey_policy}" << endl;
      exit(EXIT_FAILURE);
    }

  /*
    Enforce the same number of equations and endogenous, except in two cases:
    - ramsey_policy is used
    - a BVAR command is used and there is no equation (standalone BVAR estimation)
  */
  if (!mod_file_struct.ramsey_policy_present
      && !((mod_file_struct.bvar_density_present || mod_file_struct.bvar_forecast_present)
           && model_tree.equation_number() == 0)
      && (model_tree.equation_number() != symbol_table.endo_nbr))
    {
      cerr << "ERROR: There are " << model_tree.equation_number() << " equations but " << symbol_table.endo_nbr << " endogenous variables!" << endl;
      exit(EXIT_FAILURE);
    }
}

void
ModFile::computingPass(bool no_tmp_terms)
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
              exit(EXIT_FAILURE);
            }
          model_tree.computeJacobianExo = true;
          if (mod_file_struct.order_option >= 2)
            model_tree.computeHessian = true;
          if (mod_file_struct.order_option == 3)
            model_tree.computeThirdDerivatives = true;
        }
      //evalAllExpressions();
      model_tree.computingPass(global_eval_context, no_tmp_terms);
    }
  for(vector<Statement *>::iterator it = statements.begin();
      it != statements.end(); it++)
    (*it)->computingPass();
  //evalAllExpressions();
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
  mOutputFile << "options_.model_mode = " << model_tree.mode << ";\n";
  if (model_tree.mode == eSparseMode)
    {
      mOutputFile << "addpath " << basename << ";\n";
      mOutputFile << "delete('" << basename << "_static.m');\n";
      mOutputFile << "delete('" << basename << "_dynamic.m');\n";
    }


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

  cout << "Processing outputs ...";

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

  mOutputFile << "save('" << basename << "_results.mat', 'oo_', 'M_', 'options_');" << endl;
  mOutputFile << "diary off" << endl;

  if (model_tree.mode == eSparseMode)
    mOutputFile << "rmpath " << basename << ";\n";

  mOutputFile << endl << "disp(['Total computing time : ' dynsec2hms(toc) ]);" << endl;
  mOutputFile.close();
}
