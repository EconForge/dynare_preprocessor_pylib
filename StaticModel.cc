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

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <cstdio>
#include <cerrno>
#include <algorithm>
#include "StaticModel.hh"

// For mkdir() and chdir()
#ifdef _WIN32
# include <direct.h>
#else
# include <unistd.h>
# include <sys/stat.h>
# include <sys/types.h>
#endif

StaticModel::StaticModel(SymbolTable &symbol_table_arg,
                         NumericalConstants &num_constants_arg) :
  ModelTree(symbol_table_arg, num_constants_arg),
  global_temporary_terms(true),
  cutoff(1e-15),
  mfs(0)
{
}

void
StaticModel::compileDerivative(ofstream &code_file, int eq, int symb_id, map_idx_type &map_idx) const
{
  first_derivatives_type::const_iterator it = first_derivatives.find(make_pair(eq, symbol_table.getID(eEndogenous, symb_id)));
  if (it != first_derivatives.end())
    (it->second)->compile(code_file, false, temporary_terms, map_idx, false, false);
  else
    {
      FLDZ_ fldz;
      fldz.write(code_file);
    }
}

void
StaticModel::compileChainRuleDerivative(ofstream &code_file, int eqr, int varr, int lag, map_idx_type &map_idx) const
{
  map<pair<int, pair<int, int> >, NodeID>::const_iterator it = first_chain_rule_derivatives.find(make_pair(eqr, make_pair(varr, lag)));
  if (it != first_chain_rule_derivatives.end())
    (it->second)->compile(code_file, false, temporary_terms, map_idx, false, false);
  else
    {
      FLDZ_ fldz;
      fldz.write(code_file);
    }
}

void
StaticModel::initializeVariablesAndEquations()
{
  for(int j = 0; j < equation_number(); j++)
    {
      equation_reordered.push_back(j);
      variable_reordered.push_back(j);
    }
}

void
StaticModel::computeTemporaryTermsOrdered()
{
  map<NodeID, pair<int, int> > first_occurence;
  map<NodeID, int> reference_count;
  BinaryOpNode *eq_node;
  first_derivatives_type::const_iterator it;
  first_chain_rule_derivatives_type::const_iterator it_chr;
  ostringstream tmp_s;
  v_temporary_terms.clear();
  map_idx.clear();

  unsigned int nb_blocks = getNbBlocks();
  v_temporary_terms = vector< vector<temporary_terms_type> >(nb_blocks);

  v_temporary_terms_inuse = vector<temporary_terms_inuse_type>(nb_blocks);

  temporary_terms.clear();
  if (!global_temporary_terms)
    {
      for (unsigned int block = 0; block < nb_blocks; block++)
        {

          reference_count.clear();
          temporary_terms.clear();
          unsigned int block_size = getBlockSize(block);
          unsigned int block_nb_mfs = getBlockMfs(block);
          unsigned int block_nb_recursives = block_size - block_nb_mfs;
          v_temporary_terms[block] = vector<temporary_terms_type>(block_size);
          for (unsigned int i = 0; i < block_size; i++)
            {
              if (i < block_nb_recursives && isBlockEquationRenormalized(block, i))
                getBlockEquationRenormalizedNodeID(block, i)->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, block, v_temporary_terms,  i);
              else
                {
                  eq_node = (BinaryOpNode *) getBlockEquationNodeID(block, i);
                  eq_node->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, block, v_temporary_terms,  i);
                }
            }
          for (t_block_derivatives_equation_variable_laglead_nodeid::const_iterator it = blocks_derivatives[block].begin(); it != (blocks_derivatives[block]).end(); it++)
            {
              NodeID id = it->second.second;
              id->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, block, v_temporary_terms,  block_size-1);
            }
          set<int> temporary_terms_in_use;
          temporary_terms_in_use.clear();
          v_temporary_terms_inuse[block] = temporary_terms_in_use;
        }
    }
  else
    {
      for (unsigned int block = 0; block < nb_blocks; block++)
        {
          // Compute the temporary terms reordered
          unsigned int block_size = getBlockSize(block);
          unsigned int block_nb_mfs = getBlockMfs(block);
          unsigned int block_nb_recursives = block_size - block_nb_mfs;
          v_temporary_terms[block] = vector<temporary_terms_type>(block_size);
          for (unsigned int i = 0; i < block_size; i++)
            {
              if (i < block_nb_recursives && isBlockEquationRenormalized(block, i))
                getBlockEquationRenormalizedNodeID(block, i)->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, block, v_temporary_terms,  i);
              else
                {
                  eq_node = (BinaryOpNode *) getBlockEquationNodeID(block, i);
                  eq_node->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, block, v_temporary_terms, i);
                }
            }
          for (t_block_derivatives_equation_variable_laglead_nodeid::const_iterator it = blocks_derivatives[block].begin(); it != (blocks_derivatives[block]).end(); it++)
            {
              NodeID id = it->second.second;
              id->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, block, v_temporary_terms, block_size-1);
            }

        }
      for (unsigned int block = 0; block < nb_blocks; block++)
        {
          // Collecte the temporary terms reordered
          unsigned int block_size = getBlockSize(block);
          unsigned int block_nb_mfs = getBlockMfs(block);
          unsigned int block_nb_recursives = block_size - block_nb_mfs;
          set<int> temporary_terms_in_use;
          for (unsigned int i = 0; i < block_size; i++)
            {
              if (i < block_nb_recursives && isBlockEquationRenormalized(block, i))
                getBlockEquationRenormalizedNodeID(block, i)->collectTemporary_terms(temporary_terms, temporary_terms_in_use, block);
              else
                {
                  eq_node = (BinaryOpNode *) getBlockEquationNodeID(block, i);
                  eq_node->collectTemporary_terms(temporary_terms, temporary_terms_in_use, block);
                }
            }
          for (t_block_derivatives_equation_variable_laglead_nodeid::const_iterator it = blocks_derivatives[block].begin(); it != (blocks_derivatives[block]).end(); it++)
            {
              NodeID id = it->second.second;
              id->collectTemporary_terms(temporary_terms, temporary_terms_in_use, block);
            }
          for (int i = 0; i < (int) getBlockSize(block); i++)
            for (temporary_terms_type::const_iterator it = v_temporary_terms[block][i].begin();
                 it != v_temporary_terms[block][i].end(); it++)
              (*it)->collectTemporary_terms(temporary_terms, temporary_terms_in_use, block);
          v_temporary_terms_inuse[block] = temporary_terms_in_use;
        }
      computeTemporaryTermsMapping();
    }
}

void
StaticModel::computeTemporaryTermsMapping()
{
  // Add a mapping form node ID to temporary terms order
  int j = 0;
  for (temporary_terms_type::const_iterator it = temporary_terms.begin();
      it != temporary_terms.end(); it++)
    map_idx[(*it)->idx] = j++;
}

void
StaticModel::writeModelEquationsOrdered_M(const string &static_basename) const
{
  string tmp_s, sps;
  ostringstream tmp_output, tmp1_output, global_output;
  NodeID lhs = NULL, rhs = NULL;
  BinaryOpNode *eq_node;
  map<NodeID, int> reference_count;
  temporary_terms_type local_temporary_terms;
  ofstream  output;
  int nze;
  vector<int> feedback_variables;
  ExprNodeOutputType local_output_type;

  if (global_temporary_terms)
    {
      local_output_type = oMatlabStaticModelSparse;
      local_temporary_terms = temporary_terms;
    }
  else
    local_output_type = oMatlabDynamicModelSparseLocalTemporaryTerms;

  //----------------------------------------------------------------------
  //For each block
  for (unsigned int block = 0; block < getNbBlocks(); block++)
    {
      //recursive_variables.clear();
      feedback_variables.clear();
      //For a block composed of a single equation determines wether we have to evaluate or to solve the equation
      nze = derivative_endo[block].size();
      BlockSimulationType simulation_type = getBlockSimulationType(block);
      unsigned int block_size = getBlockSize(block);
      unsigned int block_mfs = getBlockMfs(block);
      unsigned int block_recursive = block_size - block_mfs;

      tmp1_output.str("");
      tmp1_output << static_basename << "_" << block+1 << ".m";
      output.open(tmp1_output.str().c_str(), ios::out | ios::binary);
      output << "%\n";
      output << "% " << tmp1_output.str() << " : Computes static model for Dynare\n";
      output << "%\n";
      output << "% Warning : this file is generated automatically by Dynare\n";
      output << "%           from model file (.mod)\n\n";
      output << "%/\n";
      if (simulation_type == EVALUATE_BACKWARD || simulation_type == EVALUATE_FORWARD)
        output << "function y = " << static_basename << "_" << block+1 << "(y, x, params)\n";
      else
        output << "function [residual, y, g1] = " << static_basename << "_" << block+1 << "(y, x, params)\n";

      BlockType block_type;
      if (simulation_type == SOLVE_FORWARD_COMPLETE || simulation_type == SOLVE_BACKWARD_COMPLETE)
        block_type = SIMULTANS;
      else if ((simulation_type == SOLVE_FORWARD_SIMPLE || simulation_type == SOLVE_BACKWARD_SIMPLE
                || simulation_type == EVALUATE_BACKWARD    || simulation_type == EVALUATE_FORWARD)
               && getBlockFirstEquation(block) < prologue)
        block_type = PROLOGUE;
      else if ((simulation_type == SOLVE_FORWARD_SIMPLE || simulation_type == SOLVE_BACKWARD_SIMPLE
                || simulation_type == EVALUATE_BACKWARD    || simulation_type == EVALUATE_FORWARD)
               && getBlockFirstEquation(block) >= equations.size() - epilogue)
        block_type = EPILOGUE;
      else
        block_type = SIMULTANS;
      output << "  % ////////////////////////////////////////////////////////////////////////" << endl
             << "  % //" << string("                     Block ").substr(int (log10(block + 1))) << block + 1 << " " << BlockType0(block_type)
             << "          //" << endl
             << "  % //                     Simulation type "
             << BlockSim(simulation_type) << "  //" << endl
             << "  % ////////////////////////////////////////////////////////////////////////" << endl;
      output << "  global options_;" << endl;
      //The Temporary terms
      if (simulation_type != EVALUATE_BACKWARD  && simulation_type != EVALUATE_FORWARD)
        output << "  g1 = zeros(" << block_mfs << ", " << block_mfs << ");" << endl;

      if (v_temporary_terms_inuse[block].size())
        {
          tmp_output.str("");
          for (temporary_terms_inuse_type::const_iterator it = v_temporary_terms_inuse[block].begin();
               it != v_temporary_terms_inuse[block].end(); it++)
            tmp_output << " T" << *it;
          output << "  global" << tmp_output.str() << ";\n";
        }

      if (simulation_type != EVALUATE_BACKWARD && simulation_type != EVALUATE_FORWARD)
        output << "  residual=zeros(" << block_mfs << ",1);\n";

      // The equations
      for (unsigned int i = 0; i < block_size; i++)
        {
          if (!global_temporary_terms)
            local_temporary_terms = v_temporary_terms[block][i];
          temporary_terms_type tt2;
          tt2.clear();
          if (v_temporary_terms[block].size())
            {
              output << "  " << "% //Temporary variables" << endl;
              for (temporary_terms_type::const_iterator it = v_temporary_terms[block][i].begin();
                   it != v_temporary_terms[block][i].end(); it++)
                {
                  output << "  " <<  sps;
                  (*it)->writeOutput(output, local_output_type, local_temporary_terms);
                  output << " = ";
                  (*it)->writeOutput(output, local_output_type, tt2);
                  // Insert current node into tt2
                  tt2.insert(*it);
                  output << ";" << endl;
                }
            }

          int variable_ID = getBlockVariableID(block, i);
          int equation_ID = getBlockEquationID(block, i);
          EquationType equ_type = getBlockEquationType(block, i);
          string sModel = symbol_table.getName(symbol_table.getID(eEndogenous, variable_ID));
          eq_node = (BinaryOpNode *) getBlockEquationNodeID(block, i);
          lhs = eq_node->get_arg1();
          rhs = eq_node->get_arg2();
          tmp_output.str("");
          lhs->writeOutput(tmp_output, local_output_type, local_temporary_terms);
          switch (simulation_type)
            {
            case EVALUATE_BACKWARD:
            case EVALUATE_FORWARD:
            evaluation:
              output << "  % equation " << getBlockEquationID(block, i)+1 << " variable : " << sModel
                     << " (" << variable_ID+1 << ") " << c_Equation_Type(equ_type) << endl;
              output << "  ";
              if (equ_type == E_EVALUATE)
                {
                  output << tmp_output.str();
                  output << " = ";
                  rhs->writeOutput(output, local_output_type, local_temporary_terms);
                }
              else if (equ_type == E_EVALUATE_S)
                {
                  output << "%" << tmp_output.str();
                  output << " = ";
                  if (isBlockEquationRenormalized(block, i))
                    {
                      rhs->writeOutput(output, local_output_type, local_temporary_terms);
                      output << "\n  ";
                      tmp_output.str("");
                      eq_node = (BinaryOpNode *) getBlockEquationRenormalizedNodeID(block, i);
                      lhs = eq_node->get_arg1();
                      rhs = eq_node->get_arg2();
                      lhs->writeOutput(output, local_output_type, local_temporary_terms);
                      output << " = ";
                      rhs->writeOutput(output, local_output_type, local_temporary_terms);
                    }
                }
              else
                {
                  cerr << "Type missmatch for equation " << equation_ID+1  << "\n";
                  exit(EXIT_FAILURE);
                }
              output << ";\n";
              break;
            case SOLVE_BACKWARD_SIMPLE:
            case SOLVE_FORWARD_SIMPLE:
            case SOLVE_BACKWARD_COMPLETE:
            case SOLVE_FORWARD_COMPLETE:
              if (i < block_recursive)
                goto evaluation;
              feedback_variables.push_back(variable_ID);
              output << "  % equation " << equation_ID+1 << " variable : " << sModel
                     << " (" << variable_ID+1 << ") " << c_Equation_Type(equ_type) << endl;
              output << "  " << "residual(" << i+1-block_recursive << ") = (";
              goto end;
            default:
            end:
              output << tmp_output.str();
              output << ") - (";
              rhs->writeOutput(output, local_output_type, local_temporary_terms);
              output << ");\n";
            }
        }
      // The Jacobian if we have to solve the block
      if (simulation_type == SOLVE_BACKWARD_SIMPLE   || simulation_type == SOLVE_FORWARD_SIMPLE
          || simulation_type == SOLVE_BACKWARD_COMPLETE || simulation_type == SOLVE_FORWARD_COMPLETE)
        output << "  " << sps << "% Jacobian  " << endl;
      switch (simulation_type)
        {
        case SOLVE_BACKWARD_SIMPLE:
        case SOLVE_FORWARD_SIMPLE:
        case SOLVE_BACKWARD_COMPLETE:
        case SOLVE_FORWARD_COMPLETE:
          for (t_block_derivatives_equation_variable_laglead_nodeid::const_iterator it = blocks_derivatives[block].begin(); it != (blocks_derivatives[block]).end(); it++)
            {
              unsigned int eq = it->first.first;
              unsigned int var = it->first.second;
              unsigned int eqr = getBlockEquationID(block, eq);
              unsigned int varr = getBlockVariableID(block, var);
              NodeID id = it->second.second;
              output << "    g1(" << eq+1-block_recursive << ", " << var+1-block_recursive << ") = ";
              id->writeOutput(output, local_output_type, local_temporary_terms);
              output << "; % variable=" << symbol_table.getName(symbol_table.getID(eEndogenous, varr))
                     << "(" << 0
                     << ") " << varr+1
                     << ", equation=" << eqr+1 << endl;
            }
          break;
        default:
          break;
        }
      output.close();
    }
}

void
StaticModel::writeModelEquationsCode(const string file_name, const string bin_basename, map_idx_type map_idx) const
{

  ostringstream tmp_output;
  ofstream code_file;
  bool file_open = false;

  string main_name = file_name;
  main_name += ".cod";
  code_file.open(main_name.c_str(), ios::out | ios::binary | ios::ate);
  if (!code_file.is_open())
    {
      cout << "Error : Can't open file \"" << main_name << "\" for writing\n";
      exit(EXIT_FAILURE);
    }
  int count_u;
  int u_count_int = 0;

  Write_Inf_To_Bin_File(file_name, u_count_int, file_open, false, symbol_table.endo_nbr());
  file_open = true;

  //Temporary variables declaration
  FDIMT_ fdimt(temporary_terms.size());
  fdimt.write(code_file);

  FBEGINBLOCK_ fbeginblock(symbol_table.endo_nbr(),
                           SOLVE_FORWARD_COMPLETE,
                           0,
                           symbol_table.endo_nbr(),
                           variable_reordered,
                           equation_reordered,
                           false,
                           symbol_table.endo_nbr(),
                           0,
                           0,
                           u_count_int
                           );
  fbeginblock.write(code_file);


  // Add a mapping form node ID to temporary terms order
  int j = 0;
  for (temporary_terms_type::const_iterator it = temporary_terms.begin();
       it != temporary_terms.end(); it++)
    map_idx[(*it)->idx] = j++;
  compileTemporaryTerms(code_file, temporary_terms, map_idx, false, false);

  compileModelEquations(code_file, temporary_terms, map_idx, false, false);

  FENDEQU_ fendequ;
  fendequ.write(code_file);

  vector<vector<pair<int, int> > > derivatives;
  derivatives.resize(symbol_table.endo_nbr());
  count_u = symbol_table.endo_nbr();
  for (first_derivatives_type::const_iterator it = first_derivatives.begin();
       it != first_derivatives.end(); it++)
    {
      int deriv_id = it->first.second;
      if (getTypeByDerivID(deriv_id) == eEndogenous)
        {
          NodeID d1 = it->second;
          unsigned int eq = it->first.first;
          int symb = getSymbIDByDerivID(deriv_id);
          unsigned int var = symbol_table.getTypeSpecificID(symb);
          if (!derivatives[eq].size())
            derivatives[eq].clear();
          derivatives[eq].push_back(make_pair(var, count_u));

          d1->compile(code_file, false, temporary_terms, map_idx, false, false);

          FSTPSU_ fstpsu(count_u);
          fstpsu.write(code_file);
          count_u++;
        }
    }
  for (int i = 0; i < symbol_table.endo_nbr(); i++)
    {
      FLDR_ fldr(i);
      fldr.write(code_file);
      for(vector<pair<int, int> >::const_iterator it = derivatives[i].begin();
          it != derivatives[i].end(); it++)
        {
          FLDSU_ fldsu(it->second);
          fldsu.write(code_file);
          FLDSV_ fldsv(eEndogenous, it->first);
          fldsv.write(code_file);
          FBINARY_ fbinary(oTimes);
          fbinary.write(code_file);
          if (it != derivatives[i].begin())
            {
              FBINARY_ fbinary(oPlus);
              fbinary.write(code_file);
            }
        }
      FBINARY_ fbinary(oMinus);
      fbinary.write(code_file);
      FSTPSU_ fstpsu(i);
      fstpsu.write(code_file);
    }
  FENDBLOCK_ fendblock;
  fendblock.write(code_file);
  FEND_ fend;
  fend.write(code_file);
  code_file.close();
}

void
StaticModel::writeModelEquationsCode_Block(const string file_name, const string bin_basename, map_idx_type map_idx) const
{
  struct Uff_l
  {
    int u, var, lag;
    Uff_l *pNext;
  };

  struct Uff
  {
    Uff_l *Ufl, *Ufl_First;
  };

  int i, v;
  string tmp_s;
  ostringstream tmp_output;
  ofstream code_file;
  NodeID lhs = NULL, rhs = NULL;
  BinaryOpNode *eq_node;
  Uff Uf[symbol_table.endo_nbr()];
  map<NodeID, int> reference_count;
  vector<int> feedback_variables;
  bool file_open = false;

  string main_name = file_name;
  main_name += ".cod";
  code_file.open(main_name.c_str(), ios::out | ios::binary | ios::ate);
  if (!code_file.is_open())
    {
      cout << "Error : Can't open file \"" << main_name << "\" for writing\n";
      exit(EXIT_FAILURE);
    }
  //Temporary variables declaration

  FDIMT_ fdimt(temporary_terms.size());
  fdimt.write(code_file);

  for (unsigned int block = 0; block < getNbBlocks(); block++)
    {
      feedback_variables.clear();
      if (block > 0)
        {
          FENDBLOCK_ fendblock;
          fendblock.write(code_file);
        }
      int count_u;
      int u_count_int = 0;
      BlockSimulationType simulation_type = getBlockSimulationType(block);
      unsigned int block_size = getBlockSize(block);
      unsigned int block_mfs = getBlockMfs(block);
      unsigned int block_recursive = block_size - block_mfs;

      if (simulation_type == SOLVE_TWO_BOUNDARIES_SIMPLE || simulation_type == SOLVE_TWO_BOUNDARIES_COMPLETE
          || simulation_type == SOLVE_BACKWARD_COMPLETE || simulation_type == SOLVE_FORWARD_COMPLETE)
        {
          Write_Inf_To_Bin_File_Block(file_name, bin_basename, block, u_count_int, file_open);
          file_open = true;
        }

      FBEGINBLOCK_ fbeginblock(block_mfs,
                               simulation_type,
                               getBlockFirstEquation(block),
                               block_size,
                               variable_reordered,
                               equation_reordered,
                               blocks_linear[block],
                               symbol_table.endo_nbr(),
                               0,
                               0,
                               u_count_int
                               );
      fbeginblock.write(code_file);

      // The equations
      for (i = 0; i < (int) block_size; i++)
        {
          //The Temporary terms
          temporary_terms_type tt2;
          tt2.clear();
          if (v_temporary_terms[block].size())
            {
              for (temporary_terms_type::const_iterator it = v_temporary_terms[block][i].begin();
                   it != v_temporary_terms[block][i].end(); it++)
                {
                  (*it)->compile(code_file, false, tt2, map_idx, false, false);
                  FSTPST_ fstpst((int)(map_idx.find((*it)->idx)->second));
                  fstpst.write(code_file);
                  // Insert current node into tt2
                  tt2.insert(*it);
                }
            }

          int variable_ID, equation_ID;
          EquationType equ_type;
          switch (simulation_type)
            {
            evaluation:
            case EVALUATE_BACKWARD:
            case EVALUATE_FORWARD:
              equ_type = getBlockEquationType(block, i);
              if (equ_type == E_EVALUATE)
                {
                  eq_node = (BinaryOpNode *) getBlockEquationNodeID(block, i);
                  lhs = eq_node->get_arg1();
                  rhs = eq_node->get_arg2();
                  rhs->compile(code_file, false, temporary_terms, map_idx, false, false);
                  lhs->compile(code_file, true, temporary_terms, map_idx, false, false);
                }
              else if (equ_type == E_EVALUATE_S)
                {
                  eq_node = (BinaryOpNode *) getBlockEquationRenormalizedNodeID(block, i);
                  lhs = eq_node->get_arg1();
                  rhs = eq_node->get_arg2();
                  rhs->compile(code_file, false, temporary_terms, map_idx, false, false);
                  lhs->compile(code_file, true, temporary_terms, map_idx, false, false);
                }
              break;
            case SOLVE_BACKWARD_COMPLETE:
            case SOLVE_FORWARD_COMPLETE:
              if (i < (int) block_recursive)
                goto evaluation;
              variable_ID = getBlockVariableID(block, i);
              equation_ID = getBlockEquationID(block, i);
              feedback_variables.push_back(variable_ID);
              Uf[equation_ID].Ufl = NULL;
              goto end;
            default:
            end:
              eq_node = (BinaryOpNode *) getBlockEquationNodeID(block, i);
              lhs = eq_node->get_arg1();
              rhs = eq_node->get_arg2();
              lhs->compile(code_file, false, temporary_terms, map_idx, false, false);
              rhs->compile(code_file, false, temporary_terms, map_idx, false, false);

              FBINARY_ fbinary(oMinus);
              fbinary.write(code_file);

              FSTPR_ fstpr(i - block_recursive);
              fstpr.write(code_file);
            }
        }
      FENDEQU_ fendequ;
      fendequ.write(code_file);
      // The Jacobian if we have to solve the block
      if    (simulation_type != EVALUATE_BACKWARD
             && simulation_type != EVALUATE_FORWARD)
        {
          switch (simulation_type)
            {
            case SOLVE_BACKWARD_SIMPLE:
            case SOLVE_FORWARD_SIMPLE:
              compileDerivative(code_file, getBlockEquationID(block, 0), getBlockVariableID(block, 0), map_idx);
              {
                FSTPG_ fstpg(0);
                fstpg.write(code_file);
              }
              break;

            case SOLVE_BACKWARD_COMPLETE:
            case SOLVE_FORWARD_COMPLETE:
              count_u = feedback_variables.size();
              for (t_block_derivatives_equation_variable_laglead_nodeid::const_iterator it = blocks_derivatives[block].begin(); it != (blocks_derivatives[block]).end(); it++)
                {
                  unsigned int eq = it->first.first;
                  unsigned int var = it->first.second;
                  unsigned int eqr = getBlockEquationID(block, eq);
                  unsigned int varr = getBlockVariableID(block, var);
                  if (eq >= block_recursive and var >= block_recursive)
                    {
                      if (!Uf[eqr].Ufl)
                        {
                          Uf[eqr].Ufl = (Uff_l *) malloc(sizeof(Uff_l));
                          Uf[eqr].Ufl_First = Uf[eqr].Ufl;
                        }
                      else
                        {
                          Uf[eqr].Ufl->pNext = (Uff_l *) malloc(sizeof(Uff_l));
                          Uf[eqr].Ufl = Uf[eqr].Ufl->pNext;
                        }
                      Uf[eqr].Ufl->pNext = NULL;
                      Uf[eqr].Ufl->u = count_u;
                      Uf[eqr].Ufl->var = varr;
                      compileChainRuleDerivative(code_file, eqr, varr, 0, map_idx);
                      FSTPSU_ fstpsu(count_u);
                      fstpsu.write(code_file);
                      count_u++;
                    }
                }
              for (i = 0; i < (int) block_size; i++)
                {
                  if (i >= (int) block_recursive)
                    {
                      FLDR_ fldr(i-block_recursive);
                      fldr.write(code_file);

                      FLDZ_ fldz;
                      fldz.write(code_file);

                      v = getBlockEquationID(block, i);
                      for (Uf[v].Ufl = Uf[v].Ufl_First; Uf[v].Ufl; Uf[v].Ufl = Uf[v].Ufl->pNext)
                        {
                          FLDSU_ fldsu(Uf[v].Ufl->u);
                          fldsu.write(code_file);
                          FLDSV_ fldsv(eEndogenous, Uf[v].Ufl->var);
                          fldsv.write(code_file);

                          FBINARY_ fbinary(oTimes);
                          fbinary.write(code_file);

                          FCUML_ fcuml;
                          fcuml.write(code_file);
                        }
                      Uf[v].Ufl = Uf[v].Ufl_First;
                      while (Uf[v].Ufl)
                        {
                          Uf[v].Ufl_First = Uf[v].Ufl->pNext;
                          free(Uf[v].Ufl);
                          Uf[v].Ufl = Uf[v].Ufl_First;
                        }
                      FBINARY_ fbinary(oMinus);
                      fbinary.write(code_file);

                      FSTPSU_ fstpsu(i - block_recursive);
                      fstpsu.write(code_file);

                    }
                }
              break;
            default:
              break;
            }
        }
    }
  FENDBLOCK_ fendblock;
  fendblock.write(code_file);
  FEND_ fend;
  fend.write(code_file);
  code_file.close();
}

void
StaticModel::Write_Inf_To_Bin_File_Block(const string &static_basename, const string &bin_basename, const int &num,
                                   int &u_count_int, bool &file_open) const
{
  int j;
  std::ofstream SaveCode;
  if (file_open)
    SaveCode.open((bin_basename + "_static.bin").c_str(), ios::out | ios::in | ios::binary | ios::ate);
  else
    SaveCode.open((bin_basename + "_static.bin").c_str(), ios::out | ios::binary);
  if (!SaveCode.is_open())
    {
      cout << "Error : Can't open file \"" << bin_basename << "_static.bin\" for writing\n";
      exit(EXIT_FAILURE);
    }
  u_count_int = 0;
  unsigned int block_size = getBlockSize(num);
  unsigned int block_mfs = getBlockMfs(num);
  unsigned int block_recursive = block_size - block_mfs;
  for (t_block_derivatives_equation_variable_laglead_nodeid::const_iterator it = blocks_derivatives[num].begin(); it != (blocks_derivatives[num]).end(); it++)
    {
      unsigned int eq = it->first.first;
      unsigned int var = it->first.second;
      int lag = 0;
      if (eq >= block_recursive and var >= block_recursive)
        {
          int v = eq - block_recursive;
          SaveCode.write(reinterpret_cast<char *>(&v), sizeof(v));
          int varr = var - block_recursive;
          SaveCode.write(reinterpret_cast<char *>(&varr), sizeof(varr));
          SaveCode.write(reinterpret_cast<char *>(&lag), sizeof(lag));
          int u = u_count_int + block_mfs;
          SaveCode.write(reinterpret_cast<char *>(&u), sizeof(u));
          u_count_int++;
        }
    }

  for (j = block_recursive; j < (int) block_size; j++)
    {
      unsigned int varr = getBlockVariableID(num, j);
      SaveCode.write(reinterpret_cast<char *>(&varr), sizeof(varr));
    }
  for (j = block_recursive; j < (int) block_size; j++)
    {
      unsigned int eqr = getBlockEquationID(num, j);
      SaveCode.write(reinterpret_cast<char *>(&eqr), sizeof(eqr));
    }
  SaveCode.close();
}

map<pair<int, pair<int, int > >, NodeID>
StaticModel::collect_first_order_derivatives_endogenous()
{
  map<pair<int, pair<int, int > >, NodeID> endo_derivatives;
  for (first_derivatives_type::iterator it2 = first_derivatives.begin();
       it2 != first_derivatives.end(); it2++)
    {
      if (getTypeByDerivID(it2->first.second) == eEndogenous)
        {
          int eq = it2->first.first;
          int var = symbol_table.getTypeSpecificID(it2->first.second);
          int lag = 0;
          endo_derivatives[make_pair(eq, make_pair(var, lag))] = it2->second;
        }
    }
  return endo_derivatives;
}

void
StaticModel::computingPass(const eval_context_type &eval_context, bool no_tmp_terms, bool hessian, bool block, bool bytecode)
{
  // Compute derivatives w.r. to all endogenous, and possibly exogenous and exogenous deterministic
  set<int> vars;

  for (int i = 0; i < symbol_table.endo_nbr(); i++)
    vars.insert(symbol_table.getID(eEndogenous, i));

  // Launch computations
  cout << "Computing static model derivatives:" << endl
       << " - order 1" << endl;
  first_derivatives.clear();

  computeJacobian(vars);

  if (hessian)
    {
      cout << " - order 2" << endl;
      computeHessian(vars);
    }

  if (block)
    {
      jacob_map contemporaneous_jacobian, static_jacobian;

      // for each block contains pair<Size, Feddback_variable>
      vector<pair<int, int> > blocks;

      evaluateAndReduceJacobian(eval_context, contemporaneous_jacobian, static_jacobian, dynamic_jacobian, cutoff, false);

      computeNonSingularNormalization(contemporaneous_jacobian, cutoff, static_jacobian, dynamic_jacobian);

      computePrologueAndEpilogue(static_jacobian, equation_reordered, variable_reordered, prologue, epilogue);

      map<pair<int, pair<int, int> >, NodeID> first_order_endo_derivatives = collect_first_order_derivatives_endogenous();

      equation_type_and_normalized_equation = equationTypeDetermination(equations, first_order_endo_derivatives, variable_reordered, equation_reordered, mfs);

      cout << "Finding the optimal block decomposition of the model ...\n";

      if (prologue+epilogue < (unsigned int) equation_number())
        computeBlockDecompositionAndFeedbackVariablesForEachBlock(static_jacobian, dynamic_jacobian, prologue, epilogue, equation_reordered, variable_reordered, blocks, equation_type_and_normalized_equation, false, false, mfs, inv_equation_reordered, inv_variable_reordered);

      block_type_firstequation_size_mfs = reduceBlocksAndTypeDetermination(dynamic_jacobian, prologue, epilogue, blocks, equations, equation_type_and_normalized_equation, variable_reordered, equation_reordered);

      printBlockDecomposition(blocks);

      computeChainRuleJacobian(blocks_derivatives);

      blocks_linear = BlockLinear(blocks_derivatives, variable_reordered);

      collect_block_first_order_derivatives();

      global_temporary_terms = true;
      if (!no_tmp_terms)
        computeTemporaryTermsOrdered();
    }
  else
    {
      if (!no_tmp_terms)
        {
          computeTemporaryTerms(true);
          if (bytecode)
            computeTemporaryTermsMapping();
        }
    }
}

void
StaticModel::writeStaticMFile(const string &func_name) const
{
  // Writing comments and function definition command
  string filename = func_name + "_static.m";

  ofstream output;
  output.open(filename.c_str(), ios::out | ios::binary);
  if (!output.is_open())
    {
      cerr << "ERROR: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }

  output << "function [residual, g1, g2] = " << func_name + "_static(y, x, params)" << endl
         << "%" << endl
         << "% Status : Computes static model for Dynare" << endl
         << "%" << endl
         << "% Warning : this file is generated automatically by Dynare" << endl
         << "%           from model file (.mod)" << endl
         << endl
         << "residual = zeros( " << equations.size() << ", 1);" << endl
         << endl
         << "%" << endl
         << "% Model equations" << endl
         << "%" << endl
         << endl;

  writeModelLocalVariables(output, oMatlabStaticModel);

  writeTemporaryTerms(temporary_terms, output, oMatlabStaticModel);

  writeModelEquations(output, oMatlabStaticModel);

  output << "if ~isreal(residual)" << endl
         << "  residual = real(residual)+imag(residual).^2;" << endl
         << "end" << endl
         << "if nargout >= 2," << endl
         << "  g1 = zeros(" << equations.size() << ", " << symbol_table.endo_nbr() << ");" << endl
         << endl
         << "%" << endl
         << "% Jacobian matrix" << endl
         << "%" << endl
         << endl;

  // Write Jacobian w.r. to endogenous only
  for (first_derivatives_type::const_iterator it = first_derivatives.begin();
       it != first_derivatives.end(); it++)
    {
      int eq = it->first.first;
      int symb_id = it->first.second;
      NodeID d1 = it->second;

      output << "  g1(" << eq+1 << "," << symbol_table.getTypeSpecificID(symb_id)+1 << ")=";
      d1->writeOutput(output, oMatlabStaticModel, temporary_terms);
      output << ";" << endl;
    }

  output << "  if ~isreal(g1)" << endl
         << "    g1 = real(g1)+2*imag(g1);" << endl
         << "  end" << endl
         << "end" << endl
         << "if nargout >= 3," << endl
         << "%" << endl
         << "% Hessian matrix" << endl
         << "%" << endl
         << endl;

  int g2ncols = symbol_table.endo_nbr() * symbol_table.endo_nbr();
  if (second_derivatives.size())
    {
      output << "  v2 = zeros(" << NNZDerivatives[1] << ",3);" << endl;

      // Write Hessian w.r. to endogenous only (only if 2nd order derivatives have been computed)
      int k = 0; // Keep the line of a 2nd derivative in v2
      for (second_derivatives_type::const_iterator it = second_derivatives.begin();
           it != second_derivatives.end(); it++)
        {
          int eq = it->first.first;
          int symb_id1 = it->first.second.first;
          int symb_id2 = it->first.second.second;
          NodeID d2 = it->second;

          int tsid1 = symbol_table.getTypeSpecificID(symb_id1);
          int tsid2 = symbol_table.getTypeSpecificID(symb_id2);

          int col_nb = tsid1*symbol_table.endo_nbr()+tsid2;
          int col_nb_sym = tsid2*symbol_table.endo_nbr()+tsid1;

          output << "v2(" << k+1 << ",1)=" << eq + 1 << ";" << endl
                 << "v2(" << k+1 << ",2)=" << col_nb + 1 << ";" << endl
                 << "v2(" << k+1 << ",3)=";
          d2->writeOutput(output, oMatlabStaticModel, temporary_terms);
          output << ";" << endl;

          k++;

          // Treating symetric elements
          if (symb_id1 != symb_id2)
            {
              output << "v2(" << k+1 << ",1)=" << eq + 1 << ";" << endl
                     << "v2(" << k+1 << ",2)=" << col_nb_sym + 1 << ";" << endl
                     << "v2(" << k+1 << ",3)=v2(" << k << ",3);" << endl;
              k++;
            }
        }

      output << "  g2 = sparse(v2(:,1),v2(:,2),v2(:,3)," << equations.size() << "," << g2ncols << ");" << endl;
    }
  else // Either hessian is all zero, or we didn't compute it
    output << "  g2 = sparse([],[],[]," << equations.size() << "," << g2ncols << ");" << endl;

  output << "end;" << endl; // Close the if nargout >= 3 statement
  output.close();
}

void
StaticModel::writeStaticFile(const string &basename, bool block, bool bytecode) const
{
  int r;

  //assert(block);

#ifdef _WIN32
  r = mkdir(basename.c_str());
#else
  r = mkdir(basename.c_str(), 0777);
#endif
  if (r < 0 && errno != EEXIST)
    {
      perror("ERROR");
      exit(EXIT_FAILURE);
    }
  if (block && bytecode)
    writeModelEquationsCode_Block(basename + "_static", basename, map_idx);
  else if (!block && bytecode)
    writeModelEquationsCode(basename + "_static", basename, map_idx);
  else if (block && !bytecode)
    {
      chdir(basename.c_str());
      writeModelEquationsOrdered_M(basename + "_static");
      chdir("..");
      writeStaticBlockMFSFile(basename);
    }
  else
    writeStaticMFile(basename);
}

void
StaticModel::writeStaticBlockMFSFile(const string &basename) const
{
  string filename = basename + "_static.m";

  ofstream output;
  output.open(filename.c_str(), ios::out | ios::binary);
  if (!output.is_open())
    {
      cerr << "ERROR: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }

  string func_name = basename + "_static";

  output << "function [residual, g1, y] = " << func_name << "(nblock, y, x, params)" << endl
         << "  residual = [];" << endl
         << "  g1 = [];" << endl
         << "  switch nblock" << endl;

  unsigned int nb_blocks = getNbBlocks();

  for (int b = 0; b < (int) nb_blocks; b++)
    {

      set<int> local_var;

      output << "    case " << b+1 << endl;

      BlockSimulationType simulation_type = getBlockSimulationType(b);

      if (simulation_type == EVALUATE_BACKWARD || simulation_type == EVALUATE_FORWARD)
        output << "      y = " << func_name << "_" << b+1 << "(y, x, params);\n";
      else
        output << "      [residual, y, g1] = " << func_name << "_" << b+1 << "(y, x, params);\n";
    }
  output << "  end" << endl
         << "end" << endl;
  output.close();

}

void
StaticModel::writeOutput(ostream &output, bool block) const
{
  if (!block)
    return;

  unsigned int nb_blocks = getNbBlocks();
  output << "M_.blocksMFS = cell(" << nb_blocks << ", 1);" << endl;
  for (int b = 0; b < (int) nb_blocks; b++)
    {
      output << "M_.blocksMFS{" << b+1 << "} = [ ";
      unsigned int block_size = getBlockSize(b);
      unsigned int block_mfs = getBlockMfs(b);
      unsigned int block_recursive = block_size - block_mfs;
      BlockSimulationType simulation_type = getBlockSimulationType(b);

      if (simulation_type != EVALUATE_BACKWARD && simulation_type != EVALUATE_FORWARD)
        for (int i = block_recursive; i < (int) block_size; i++)
          output << getBlockVariableID(b, i)+1 << "; ";

      output << "];" << endl;
    }
}

SymbolType
StaticModel::getTypeByDerivID(int deriv_id) const throw (UnknownDerivIDException)
{
  return symbol_table.getType(getSymbIDByDerivID(deriv_id));
}

int
StaticModel::getLagByDerivID(int deriv_id) const throw (UnknownDerivIDException)
{
  return 0;
}

int
StaticModel::getSymbIDByDerivID(int deriv_id) const throw (UnknownDerivIDException)
{
  return deriv_id;
}

int
StaticModel::getDerivID(int symb_id, int lag) const throw (UnknownDerivIDException)
{
  if (symbol_table.getType(symb_id) == eEndogenous)
    return symb_id;
  else
    return -1;
}

map<pair<pair<int, pair<int, int> >, pair<int, int> >, int>
StaticModel::get_Derivatives(int block)
{
  map<pair<pair<int, pair<int, int> >, pair<int, int> >, int> Derivatives;
  Derivatives.clear();
  int block_size = getBlockSize(block);
  int block_nb_recursive = block_size - getBlockMfs(block);
  int lag = 0;
  for (int eq = 0; eq < block_size; eq++)
    {
      int eqr = getBlockEquationID(block, eq);
      for (int var = 0; var < block_size; var++)
        {
          int varr = getBlockVariableID(block, var);
          if (dynamic_jacobian.find(make_pair(lag, make_pair(eqr, varr))) != dynamic_jacobian.end())
            {
              bool OK = true;
              map<pair<pair<int, pair<int, int> >, pair<int, int> >, int>::const_iterator its = Derivatives.find(make_pair(make_pair(lag, make_pair(eq, var)), make_pair(eqr, varr)));
              if (its != Derivatives.end())
                {
                  if (its->second == 2)
                    OK = false;
                }

              if (OK)
                {
                  if (getBlockEquationType(block, eq) == E_EVALUATE_S and eq < block_nb_recursive)
                    //It's a normalized equation, we have to recompute the derivative using chain rule derivative function
                    Derivatives[make_pair(make_pair(lag, make_pair(eq, var)), make_pair(eqr, varr))] = 1;
                  else
                    //It's a feedback equation we can use the derivatives
                    Derivatives[make_pair(make_pair(lag, make_pair(eq, var)), make_pair(eqr, varr))] = 0;
                }
              if (var < block_nb_recursive)
                {
                  int eqs = getBlockEquationID(block, var);
                  for (int vars = block_nb_recursive; vars < block_size; vars++)
                    {
                      int varrs = getBlockVariableID(block, vars);
                      //A new derivative needs to be computed using the chain rule derivative function (a feedback variable appears in a recursive equation)
                      if (Derivatives.find(make_pair(make_pair(lag, make_pair(var, vars)), make_pair(eqs, varrs))) != Derivatives.end())
                        Derivatives[make_pair(make_pair(lag, make_pair(eq, vars)), make_pair(eqr, varrs))] = 2;
                    }
                }
            }
        }
    }

  return (Derivatives);
}

void
StaticModel::computeChainRuleJacobian(t_blocks_derivatives &blocks_derivatives)
{
  map<int, NodeID> recursive_variables;
  unsigned int nb_blocks = getNbBlocks();
  blocks_derivatives = t_blocks_derivatives(nb_blocks);
  for (unsigned int block = 0; block < nb_blocks; block++)
    {
      t_block_derivatives_equation_variable_laglead_nodeid tmp_derivatives;
      recursive_variables.clear();
      BlockSimulationType simulation_type = getBlockSimulationType(block);
      int block_size = getBlockSize(block);
      int block_nb_mfs = getBlockMfs(block);
      int block_nb_recursives = block_size - block_nb_mfs;
      if (simulation_type == SOLVE_TWO_BOUNDARIES_COMPLETE or simulation_type == SOLVE_TWO_BOUNDARIES_SIMPLE)
        {
          blocks_derivatives.push_back(t_block_derivatives_equation_variable_laglead_nodeid(0));
          for (int i = 0; i < block_nb_recursives; i++)
            {
              if (getBlockEquationType(block, i) == E_EVALUATE_S)
                recursive_variables[getDerivID(symbol_table.getID(eEndogenous, getBlockVariableID(block, i)), 0)] = getBlockEquationRenormalizedNodeID(block, i);
              else
                recursive_variables[getDerivID(symbol_table.getID(eEndogenous, getBlockVariableID(block, i)), 0)] = getBlockEquationNodeID(block, i);
            }
          map<pair<pair<int, pair<int, int> >, pair<int, int> >, int> Derivatives = get_Derivatives(block);
          map<pair<pair<int, pair<int, int> >, pair<int, int> >, int>::const_iterator it = Derivatives.begin();
          for (int i = 0; i < (int) Derivatives.size(); i++)
            {
              int Deriv_type = it->second;
              pair<pair<int, pair<int, int> >, pair<int, int> > it_l(it->first);
              it++;
              int lag = it_l.first.first;
              int eq = it_l.first.second.first;
              int var = it_l.first.second.second;
              int eqr = it_l.second.first;
              int varr = it_l.second.second;
              if (Deriv_type == 0)
                first_chain_rule_derivatives[make_pair(eqr, make_pair(varr, lag))] = first_derivatives[make_pair(eqr, getDerivID(symbol_table.getID(eEndogenous, varr), lag))];
              else if (Deriv_type == 1)
                first_chain_rule_derivatives[make_pair(eqr, make_pair(varr, lag))] = (equation_type_and_normalized_equation[eqr].second)->getChainRuleDerivative(getDerivID(symbol_table.getID(eEndogenous, varr), lag), recursive_variables);
              else if (Deriv_type == 2)
                {
                  if (getBlockEquationType(block, eq) == E_EVALUATE_S && eq < block_nb_recursives)
                    first_chain_rule_derivatives[make_pair(eqr, make_pair(varr, lag))] = (equation_type_and_normalized_equation[eqr].second)->getChainRuleDerivative(getDerivID(symbol_table.getID(eEndogenous, varr), lag), recursive_variables);
                  else
                    first_chain_rule_derivatives[make_pair(eqr, make_pair(varr, lag))] = equations[eqr]->getChainRuleDerivative(getDerivID(symbol_table.getID(eEndogenous, varr), lag), recursive_variables);
                }
              tmp_derivatives.push_back(make_pair(make_pair(eq, var), make_pair(lag, first_chain_rule_derivatives[make_pair(eqr, make_pair(varr, lag))])));
            }
        }
      else if (simulation_type == SOLVE_BACKWARD_SIMPLE or simulation_type == SOLVE_FORWARD_SIMPLE
               or simulation_type == SOLVE_BACKWARD_COMPLETE or simulation_type == SOLVE_FORWARD_COMPLETE)
        {
          blocks_derivatives.push_back(t_block_derivatives_equation_variable_laglead_nodeid(0));
          for (int i = 0; i < block_nb_recursives; i++)
            {
              if (getBlockEquationType(block, i) == E_EVALUATE_S)
                recursive_variables[getDerivID(symbol_table.getID(eEndogenous, getBlockVariableID(block, i)), 0)] = getBlockEquationRenormalizedNodeID(block, i);
              else
                recursive_variables[getDerivID(symbol_table.getID(eEndogenous, getBlockVariableID(block, i)), 0)] = getBlockEquationNodeID(block, i);
            }
          for (int eq = block_nb_recursives; eq < block_size; eq++)
            {
              int eqr = getBlockEquationID(block, eq);
              for (int var = block_nb_recursives; var < block_size; var++)
                {
                  int varr = getBlockVariableID(block, var);
                  NodeID d1 = equations[eqr]->getChainRuleDerivative(getDerivID(symbol_table.getID(eEndogenous, varr), 0), recursive_variables);
                  if (d1 == Zero)
                    continue;
                  first_chain_rule_derivatives[make_pair(eqr, make_pair(varr, 0))] = d1;
                  tmp_derivatives.push_back(
                                            make_pair(make_pair(eq, var), make_pair(0, first_chain_rule_derivatives[make_pair(eqr, make_pair(varr, 0))])));
                }
            }
        }
      blocks_derivatives[block] = tmp_derivatives;
    }
}

void
StaticModel::collect_block_first_order_derivatives()
{
  //! vector for an equation or a variable indicates the block number
  vector<int> equation_2_block, variable_2_block;
  unsigned int nb_blocks = getNbBlocks();
  equation_2_block = vector<int>(equation_reordered.size());
  variable_2_block = vector<int>(variable_reordered.size());
  for (unsigned int block = 0; block < nb_blocks; block++)
    {
      unsigned int block_size = getBlockSize(block);
      for (unsigned int i = 0; i < block_size; i++)
        {
          equation_2_block[getBlockEquationID(block, i)] = block;
          variable_2_block[getBlockVariableID(block, i)] = block;
        }
    }
  derivative_endo = vector<t_derivative>(nb_blocks);
  endo_max_leadlag_block = vector<pair<int, int> >(nb_blocks, make_pair(0, 0));
  max_leadlag_block = vector<pair<int, int> >(nb_blocks, make_pair(0, 0));
  for (first_derivatives_type::iterator it2 = first_derivatives.begin();
       it2 != first_derivatives.end(); it2++)
    {
      int eq = it2->first.first;
      int var = symbol_table.getTypeSpecificID(getSymbIDByDerivID(it2->first.second));
      int lag = 0;
      int block_eq = equation_2_block[eq];
      int block_var = variable_2_block[var];
      max_leadlag_block[block_eq] = make_pair(0, 0);
      max_leadlag_block[block_eq] = make_pair(0, 0);
      endo_max_leadlag_block[block_eq] = make_pair(0, 0);
      endo_max_leadlag_block[block_eq] = make_pair(0, 0);
      t_derivative tmp_derivative;
      t_lag_var lag_var;
      if (getTypeByDerivID(it2->first.second) == eEndogenous && block_eq == block_var)
        {
          tmp_derivative = derivative_endo[block_eq];
          tmp_derivative[make_pair(lag, make_pair(eq, var))] = first_derivatives[make_pair(eq, getDerivID(symbol_table.getID(eEndogenous, var), lag))];
          derivative_endo[block_eq] = tmp_derivative;
        }
    }
}

void
StaticModel::writeChainRuleDerivative(ostream &output, int eqr, int varr, int lag,
                                      ExprNodeOutputType output_type,
                                      const temporary_terms_type &temporary_terms) const
{
  map<pair<int, pair<int, int> >, NodeID>::const_iterator it = first_chain_rule_derivatives.find(make_pair(eqr, make_pair(varr, lag)));
  if (it != first_chain_rule_derivatives.end())
    (it->second)->writeOutput(output, output_type, temporary_terms);
  else
    output << 0;
}

void
StaticModel::writeLatexFile(const string &basename) const
{
  writeLatexModelFile(basename + "_static.tex", oLatexStaticModel);
}

void
StaticModel::jacobianHelper(ostream &output, int eq_nb, int col_nb, ExprNodeOutputType output_type) const
{
  output << LEFT_ARRAY_SUBSCRIPT(output_type);
  if (IS_MATLAB(output_type))
    output << eq_nb + 1 << ", " << col_nb + 1;
  else
    output << eq_nb + col_nb * equations.size();
  output << RIGHT_ARRAY_SUBSCRIPT(output_type);
}

void
StaticModel::hessianHelper(ostream &output, int row_nb, int col_nb, ExprNodeOutputType output_type) const
{
  output << LEFT_ARRAY_SUBSCRIPT(output_type);
  if (IS_MATLAB(output_type))
    output << row_nb + 1 << ", " << col_nb + 1;
  else
    output << row_nb + col_nb * NNZDerivatives[1];
  output << RIGHT_ARRAY_SUBSCRIPT(output_type);
}

void
StaticModel::writeAuxVarInitval(ostream &output) const
{
  for (int i = 0; i < (int) aux_equations.size(); i++)
    {
      dynamic_cast<ExprNode *>(aux_equations[i])->writeOutput(output);
      output << ";" << endl;
    }
}
