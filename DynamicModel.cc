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
#include "DynamicModel.hh"

// For mkdir() and chdir()
#ifdef _WIN32
# include <direct.h>
#else
# include <unistd.h>
# include <sys/stat.h>
# include <sys/types.h>
#endif

DynamicModel::DynamicModel(SymbolTable &symbol_table_arg,
                           NumericalConstants &num_constants_arg) :
    ModelTree(symbol_table_arg, num_constants_arg),
    max_lag(0), max_lead(0),
    max_endo_lag(0), max_endo_lead(0),
    max_exo_lag(0), max_exo_lead(0),
    max_exo_det_lag(0), max_exo_det_lead(0),
    dynJacobianColsNbr(0),
    cutoff(1e-15),
    mfs(0),
    block_triangular(symbol_table_arg, num_constants_arg)
{
}

VariableNode *
DynamicModel::AddVariable(int symb_id, int lag)
{
  return AddVariableInternal(symb_id, lag);
}

void
DynamicModel::compileDerivative(ofstream &code_file, int eq, int symb_id, int lag, map_idx_type &map_idx) const
  {
    //first_derivatives_type::const_iterator it = first_derivatives.find(make_pair(eq, getDerivID(symb_id, lag)));
    first_derivatives_type::const_iterator it = first_derivatives.find(make_pair(eq, getDerivID(symbol_table.getID(eEndogenous, symb_id), lag)));
    if (it != first_derivatives.end())
      (it->second)->compile(code_file, false, temporary_terms, map_idx, true, false);
    else
      /*code_file.write(&FLDZ, sizeof(FLDZ));*/
      {
        FLDZ_ fldz;
        fldz.write(code_file);
      }
  }


void
DynamicModel::compileChainRuleDerivative(ofstream &code_file, int eqr, int varr, int lag, map_idx_type &map_idx) const
{
  map<pair<int, pair<int, int> >, NodeID>::const_iterator it = first_chain_rule_derivatives.find(make_pair(eqr, make_pair(varr, lag)));
  if (it != first_chain_rule_derivatives.end())
    (it->second)->compile(code_file, false, temporary_terms, map_idx, true, false);
  else
    {
      FLDZ_ fldz;
      fldz.write(code_file);
    }
    //code_file.write(&FLDZ, sizeof(FLDZ));
}


void
DynamicModel::BuildIncidenceMatrix()
{
  set<pair<int, int> > endogenous, exogenous;
  for (int eq = 0; eq < (int) equations.size(); eq++)
    {
      BinaryOpNode *eq_node = equations[eq];
      endogenous.clear();
      NodeID Id = eq_node->get_arg1();
      Id->collectEndogenous(endogenous);
      Id = eq_node->get_arg2();
      Id->collectEndogenous(endogenous);
      for (set<pair<int, int> >::iterator it_endogenous=endogenous.begin();it_endogenous!=endogenous.end();it_endogenous++)
        {
          block_triangular.incidencematrix.fill_IM(eq, it_endogenous->first, it_endogenous->second, eEndogenous);
        }
      exogenous.clear();
      Id = eq_node->get_arg1();
      Id->collectExogenous(exogenous);
      Id = eq_node->get_arg2();
      Id->collectExogenous(exogenous);
      for (set<pair<int, int> >::iterator it_exogenous=exogenous.begin();it_exogenous!=exogenous.end();it_exogenous++)
        {
          block_triangular.incidencematrix.fill_IM(eq, it_exogenous->first, it_exogenous->second, eExogenous);
        }
    }
}

void
DynamicModel::computeTemporaryTermsOrdered(Model_Block *ModelBlock)
{
  map<NodeID, pair<int, int> > first_occurence;
  map<NodeID, int> reference_count;
  int i, j, m, eq, var, eqr, varr, lag;
  temporary_terms_type vect;
  ostringstream tmp_output;
  BinaryOpNode *eq_node;
  first_derivatives_type::const_iterator it;
  first_chain_rule_derivatives_type::const_iterator it_chr;
  ostringstream tmp_s;

  temporary_terms.clear();
  map_idx.clear();
  for (j = 0;j < ModelBlock->Size;j++)
    {
      // Compute the temporary terms reordered
      for (i = 0;i < ModelBlock->Block_List[j].Size;i++)
        {
          if (ModelBlock->Block_List[j].Equation_Type[i] == E_EVALUATE_S && i<ModelBlock->Block_List[j].Nb_Recursives && ModelBlock->Block_List[j].Equation_Normalized[i])
              ModelBlock->Block_List[j].Equation_Normalized[i]->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, j, ModelBlock, i, map_idx);
          else
            {
              eq_node = equations[ModelBlock->Block_List[j].Equation[i]];
              eq_node->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, j, ModelBlock, i, map_idx);
            }
        }
      for(i=0; i<(int)ModelBlock->Block_List[j].Chain_Rule_Derivatives->size();i++)
        {
          pair< pair<int, pair<int, int> >, pair<int, int> > it = ModelBlock->Block_List[j].Chain_Rule_Derivatives->at(i);
          lag=it.first.first;
          int eqr=it.second.first;
          int varr=it.second.second;
          it_chr=first_chain_rule_derivatives.find(make_pair(eqr, make_pair( varr, lag)));
          it_chr->second->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, j, ModelBlock, ModelBlock->Block_List[j].Size-1, map_idx);
        }

      for (m=0;m<=ModelBlock->Block_List[j].Max_Lead+ModelBlock->Block_List[j].Max_Lag;m++)
        {
          lag=m-ModelBlock->Block_List[j].Max_Lag;
          for (i=0;i<ModelBlock->Block_List[j].IM_lead_lag[m].size;i++)
            {
              eq=ModelBlock->Block_List[j].IM_lead_lag[m].Equ_Index[i];
              var=ModelBlock->Block_List[j].IM_lead_lag[m].Var_Index[i];
              it=first_derivatives.find(make_pair(eq,getDerivID(symbol_table.getID(eEndogenous, var), lag)));
              it->second->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, j, ModelBlock, ModelBlock->Block_List[j].Size-1, map_idx);
            }
        }
      /*for(i=0; i<(int)ModelBlock->Block_List[j].Chain_Rule_Derivatives->size();i++)
        {
          pair< pair<int, pair<int, int> >, pair<int, int> > it = ModelBlock->Block_List[j].Chain_Rule_Derivatives->at(i);
          lag=it.first.first;
          eqr=it.second.first;
          varr=it.second.second;
          it_chr=first_chain_rule_derivatives.find(make_pair(eqr, make_pair( varr, lag)));
          it_chr->second->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, j, ModelBlock, ModelBlock->Block_List[j].Size-1, map_idx);
        }*/
      /*for (m=0;m<=ModelBlock->Block_List[j].Max_Lead+ModelBlock->Block_List[j].Max_Lag;m++)
        {
          lag=m-ModelBlock->Block_List[j].Max_Lag;
          for (i=0;i<ModelBlock->Block_List[j].IM_lead_lag[m].size_exo;i++)
            {
              eq=ModelBlock->Block_List[j].IM_lead_lag[m].Equ_X_Index[i];
              var=ModelBlock->Block_List[j].IM_lead_lag[m].Exogenous_Index[i];
              it=first_derivatives.find(make_pair(eq,getDerivID(symbol_table.getID(eExogenous, var), lag)));
              it->second->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, j, ModelBlock, ModelBlock->Block_List[j].Size-1, map_idx);
            }
        }*/
      //jacobian_max_exo_col=(variable_table.max_exo_lag+variable_table.max_exo_lead+1)*symbol_table.exo_nbr;
      for (m=0;m<=ModelBlock->Block_List[j].Max_Lead+ModelBlock->Block_List[j].Max_Lag;m++)
        {
          lag=m-ModelBlock->Block_List[j].Max_Lag;
          if (block_triangular.incidencematrix.Model_Max_Lag_Endo - ModelBlock->Block_List[j].Max_Lag +m >=0)
            {
              for (i=0;i<ModelBlock->Block_List[j].IM_lead_lag[m].size_other_endo;i++)
                {
                  eq=ModelBlock->Block_List[j].IM_lead_lag[m].Equ_Index_other_endo[i];
                  var=ModelBlock->Block_List[j].IM_lead_lag[m].Var_Index_other_endo[i];
                  it=first_derivatives.find(make_pair(eq,getDerivID(symbol_table.getID(eEndogenous, var), lag)));
                  it->second->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, j, ModelBlock, ModelBlock->Block_List[j].Size-1, map_idx);
                }
            }
        }
    }
  for (j = 0;j < ModelBlock->Size;j++)
    {
      // Collecte the temporary terms reordered
      for (i = 0;i < ModelBlock->Block_List[j].Size;i++)
        {
          if (ModelBlock->Block_List[j].Equation_Type[i] == E_EVALUATE_S && i<ModelBlock->Block_List[j].Nb_Recursives && ModelBlock->Block_List[j].Equation_Normalized[i])
              ModelBlock->Block_List[j].Equation_Normalized[i]->collectTemporary_terms(temporary_terms, ModelBlock, j);
          else
            {
              eq_node = equations[ModelBlock->Block_List[j].Equation[i]];
              eq_node->collectTemporary_terms(temporary_terms, ModelBlock, j);
            }

          /*eq_node = equations[ModelBlock->Block_List[j].Equation[i]];
          eq_node->collectTemporary_terms(temporary_terms, ModelBlock, j);
          if (ModelBlock->Block_List[j].Equation_Type[i] == E_EVALUATE_S)
            if(ModelBlock->Block_List[j].Equation_Normalized[i])
              ModelBlock->Block_List[j].Equation_Normalized[i]->collectTemporary_terms(temporary_terms, ModelBlock, j);
          for(temporary_terms_type::const_iterator it = ModelBlock->Block_List[j].Temporary_Terms_in_Equation[i]->begin(); it!= ModelBlock->Block_List[j].Temporary_Terms_in_Equation[i]->end(); it++)
            (*it)->collectTemporary_terms(temporary_terms, ModelBlock, j);*/
        }
      for (m=0;m<=ModelBlock->Block_List[j].Max_Lead+ModelBlock->Block_List[j].Max_Lag;m++)
        {
          lag=m-ModelBlock->Block_List[j].Max_Lag;
          for (i=0;i<ModelBlock->Block_List[j].IM_lead_lag[m].size;i++)
            {
              eq=ModelBlock->Block_List[j].IM_lead_lag[m].Equ_Index[i];
              var=ModelBlock->Block_List[j].IM_lead_lag[m].Var_Index[i];
              it=first_derivatives.find(make_pair(eq,getDerivID(symbol_table.getID(eEndogenous, var), lag)));
              //it=first_derivatives.find(make_pair(eq,variable_table.getID(var, lag)));
              //if(it!=first_derivatives.end())
              it->second->collectTemporary_terms(temporary_terms, ModelBlock, j);
            }
        }
      for(i=0; i<(int)ModelBlock->Block_List[j].Chain_Rule_Derivatives->size();i++)
        {
          pair< pair<int, pair<int, int> >, pair<int, int> > it = ModelBlock->Block_List[j].Chain_Rule_Derivatives->at(i);
          lag=it.first.first;
          eqr=it.second.first;
          varr=it.second.second;
          it_chr=first_chain_rule_derivatives.find(make_pair(eqr, make_pair( varr, lag)));
          it_chr->second->collectTemporary_terms(temporary_terms, ModelBlock, j);
        }
      /*for (m=0;m<=ModelBlock->Block_List[j].Max_Lead+ModelBlock->Block_List[j].Max_Lag;m++)
        {
          lag=m-ModelBlock->Block_List[j].Max_Lag;
          for (i=0;i<ModelBlock->Block_List[j].IM_lead_lag[m].size_exo;i++)
            {
              eq=ModelBlock->Block_List[j].IM_lead_lag[m].Equ_X_Index[i];
              var=ModelBlock->Block_List[j].IM_lead_lag[m].Exogenous_Index[i];
              it=first_derivatives.find(make_pair(eq,getDerivID(symbol_table.getID(eExogenous, var), lag)));
              //it=first_derivatives.find(make_pair(eq,variable_table.getID(var, lag)));
              it->second->collectTemporary_terms(temporary_terms, ModelBlock, j);
            }
        }*/
      //jacobian_max_exo_col=(variable_table.max_exo_lag+variable_table.max_exo_lead+1)*symbol_table.exo_nbr;
      for (m=0;m<=ModelBlock->Block_List[j].Max_Lead+ModelBlock->Block_List[j].Max_Lag;m++)
        {
          lag=m-ModelBlock->Block_List[j].Max_Lag;
          if (block_triangular.incidencematrix.Model_Max_Lag_Endo - ModelBlock->Block_List[j].Max_Lag +m >=0)
            {
              for (i=0;i<ModelBlock->Block_List[j].IM_lead_lag[m].size_other_endo;i++)
                {
                  eq=ModelBlock->Block_List[j].IM_lead_lag[m].Equ_Index_other_endo[i];
                  var=ModelBlock->Block_List[j].IM_lead_lag[m].Var_Index_other_endo[i];
                  it=first_derivatives.find(make_pair(eq,getDerivID(symbol_table.getID(eEndogenous, var), lag)));
                  //it=first_derivatives.find(make_pair(eq,variable_table.getID(var, lag)));
                  //if(it!=first_derivatives.end())
                  it->second->collectTemporary_terms(temporary_terms, ModelBlock, j);
                }
            }
        }
    }
  // Add a mapping form node ID to temporary terms order
  j=0;
  for (temporary_terms_type::const_iterator it = temporary_terms.begin();
       it != temporary_terms.end(); it++)
    map_idx[(*it)->idx]=j++;
}

void
DynamicModel::writeModelEquationsOrdered_M( Model_Block *ModelBlock, const string &dynamic_basename) const
  {
    int i,j,k,m;
    string tmp_s, sps;
    ostringstream tmp_output, tmp1_output, global_output;
    NodeID lhs=NULL, rhs=NULL;
    BinaryOpNode *eq_node;
    ostringstream Uf[symbol_table.endo_nbr()];
    map<NodeID, int> reference_count;
    //int prev_Simulation_Type=-1, count_derivates=0;
    int jacobian_max_endo_col;
    ofstream  output;
    //temporary_terms_type::const_iterator it_temp=temporary_terms.begin();
    int nze, nze_exo, nze_other_endo;
    //map<int, NodeID> recursive_variables;
    vector<int> feedback_variables;
    //----------------------------------------------------------------------
    //For each block
    for (j = 0;j < ModelBlock->Size;j++)
      {
        //recursive_variables.clear();
        feedback_variables.clear();
        //For a block composed of a single equation determines wether we have to evaluate or to solve the equation
        nze = nze_exo = nze_other_endo = 0;
        for (m=0;m<=ModelBlock->Block_List[j].Max_Lead+ModelBlock->Block_List[j].Max_Lag;m++)
          nze+=ModelBlock->Block_List[j].IM_lead_lag[m].size;
        /*for (m=0;m<=ModelBlock->Block_List[j].Max_Lead_Exo+ModelBlock->Block_List[j].Max_Lag_Exo;m++)
          nze_exo+=ModelBlock->Block_List[j].IM_lead_lag[m].size_exo;*/
        for (m=0;m<=ModelBlock->Block_List[j].Max_Lead+ModelBlock->Block_List[j].Max_Lag;m++)
          {
            k=m-ModelBlock->Block_List[j].Max_Lag;
            if (block_triangular.incidencematrix.Model_Max_Lag_Endo - ModelBlock->Block_List[j].Max_Lag +m >=0)
              nze_other_endo+=ModelBlock->Block_List[j].IM_lead_lag[m].size_other_endo;
          }
        tmp1_output.str("");
        tmp1_output << dynamic_basename << "_" << j+1 << ".m";
        output.open(tmp1_output.str().c_str(), ios::out | ios::binary);
        output << "%\n";
        output << "% " << tmp1_output.str() << " : Computes dynamic model for Dynare\n";
        output << "%\n";
        output << "% Warning : this file is generated automatically by Dynare\n";
        output << "%           from model file (.mod)\n\n";
        output << "%/\n";
        if (ModelBlock->Block_List[j].Simulation_Type==EVALUATE_BACKWARD
            ||ModelBlock->Block_List[j].Simulation_Type==EVALUATE_FORWARD
            /*||ModelBlock->Block_List[j].Simulation_Type==EVALUATE_BACKWARD_R
            ||ModelBlock->Block_List[j].Simulation_Type==EVALUATE_FORWARD_R*/)
          {
            output << "function [y, g1, g2, g3, varargout] = " << dynamic_basename << "_" << j+1 << "(y, x, params, jacobian_eval, y_kmin, periods)\n";
          }
        else if (ModelBlock->Block_List[j].Simulation_Type==SOLVE_FORWARD_COMPLETE
                 ||   ModelBlock->Block_List[j].Simulation_Type==SOLVE_BACKWARD_COMPLETE)
          output << "function [residual, y, g1, g2, g3, varargout] = " << dynamic_basename << "_" << j+1 << "(y, x, params, it_, jacobian_eval)\n";
        else if (ModelBlock->Block_List[j].Simulation_Type==SOLVE_BACKWARD_SIMPLE
                 ||   ModelBlock->Block_List[j].Simulation_Type==SOLVE_FORWARD_SIMPLE)
          output << "function [residual, y, g1, g2, g3, varargout] = " << dynamic_basename << "_" << j+1 << "(y, x, params, it_, jacobian_eval)\n";
        else
          output << "function [residual, y, g1, g2, g3, b, varargout] = " << dynamic_basename << "_" << j+1 << "(y, x, params, periods, jacobian_eval, y_kmin, y_size)\n";
        output << "  % ////////////////////////////////////////////////////////////////////////" << endl
        << "  % //" << string("                     Block ").substr(int(log10(j + 1))) << j + 1 << " " << BlockTriangular::BlockType0(ModelBlock->Block_List[j].Type)
        << "          //" << endl
        << "  % //                     Simulation type "
        << BlockTriangular::BlockSim(ModelBlock->Block_List[j].Simulation_Type) << "  //" << endl
        << "  % ////////////////////////////////////////////////////////////////////////" << endl;
        output << "  global options_;" << endl;
        //The Temporary terms
        //output << "  relax = 1;\n";
        if (ModelBlock->Block_List[j].Simulation_Type==EVALUATE_BACKWARD
            ||ModelBlock->Block_List[j].Simulation_Type==EVALUATE_FORWARD)
          {
            output << "  if(jacobian_eval)\n";
            output << "    g1 = spalloc(" << ModelBlock->Block_List[j].Size-ModelBlock->Block_List[j].Nb_Recursives
            << ", " << (ModelBlock->Block_List[j].Size-ModelBlock->Block_List[j].Nb_Recursives)*(1+ModelBlock->Block_List[j].Max_Lag_Endo+ModelBlock->Block_List[j].Max_Lead_Endo)
            << ", " << nze << ");\n";
            output << "    g1_x=spalloc(" << ModelBlock->Block_List[j].Size << ", " << (ModelBlock->Block_List[j].nb_exo + ModelBlock->Block_List[j].nb_exo_det)*(1+ModelBlock->Block_List[j].Max_Lag_Exo+ModelBlock->Block_List[j].Max_Lead_Exo) << ", " << nze_exo << ");\n";
            output << "    g1_o=spalloc(" << ModelBlock->Block_List[j].Size << ", " << ModelBlock->Block_List[j].nb_other_endo*(1+ModelBlock->Block_List[j].Max_Lag_Other_Endo+ModelBlock->Block_List[j].Max_Lead_Other_Endo) << ", " << nze_other_endo << ");\n";
            output << "  end;\n";
          }
        else
          {
            output << "  if(jacobian_eval)\n";
            output << "    g1 = spalloc(" << ModelBlock->Block_List[j].Size << ", " << ModelBlock->Block_List[j].Size*(1+ModelBlock->Block_List[j].Max_Lag_Endo+ModelBlock->Block_List[j].Max_Lead_Endo) << ", " << nze << ");\n";
            output << "    g1_x=spalloc(" << ModelBlock->Block_List[j].Size << ", " << (ModelBlock->Block_List[j].nb_exo + ModelBlock->Block_List[j].nb_exo_det)*(1+ModelBlock->Block_List[j].Max_Lag_Exo+ModelBlock->Block_List[j].Max_Lead_Exo) << ", " << nze_exo << ");\n";
            output << "    g1_o=spalloc(" << ModelBlock->Block_List[j].Size << ", " << ModelBlock->Block_List[j].nb_other_endo*(1+ModelBlock->Block_List[j].Max_Lag_Other_Endo+ModelBlock->Block_List[j].Max_Lead_Other_Endo) << ", " << nze_other_endo << ");\n";
            output << "  else\n";
            if (ModelBlock->Block_List[j].Simulation_Type==SOLVE_TWO_BOUNDARIES_COMPLETE || ModelBlock->Block_List[j].Simulation_Type==SOLVE_TWO_BOUNDARIES_SIMPLE)
              {
                output << "    g1 = spalloc(" << (ModelBlock->Block_List[j].Size-ModelBlock->Block_List[j].Nb_Recursives) << "*options_.periods, "
                << (ModelBlock->Block_List[j].Size-ModelBlock->Block_List[j].Nb_Recursives) << "*(options_.periods+" << ModelBlock->Block_List[j].Max_Lag+ModelBlock->Block_List[j].Max_Lead+1 << ")"
                << ", " << nze << "*options_.periods);\n";
                /*output << "    g1_tmp_r = spalloc(" << (ModelBlock->Block_List[j].Nb_Recursives)
                << ", " << (ModelBlock->Block_List[j].Size)*(ModelBlock->Block_List[j].Max_Lag+ModelBlock->Block_List[j].Max_Lead+1)
                << ", " << nze << ");\n";
                output << "    g1_tmp_b = spalloc(" << (ModelBlock->Block_List[j].Size-ModelBlock->Block_List[j].Nb_Recursives)
                << ", " << (ModelBlock->Block_List[j].Size)*(ModelBlock->Block_List[j].Max_Lag+ModelBlock->Block_List[j].Max_Lead+1)
                << ", " << nze << ");\n";*/
              }
            else
              {
                output << "    g1 = spalloc(" << ModelBlock->Block_List[j].Size-ModelBlock->Block_List[j].Nb_Recursives
                << ", " << ModelBlock->Block_List[j].Size-ModelBlock->Block_List[j].Nb_Recursives << ", " << nze << ");\n";
                output << "    g1_tmp_r = spalloc(" << ModelBlock->Block_List[j].Nb_Recursives
                << ", " << ModelBlock->Block_List[j].Size << ", " << nze << ");\n";
                output << "    g1_tmp_b = spalloc(" << ModelBlock->Block_List[j].Size-ModelBlock->Block_List[j].Nb_Recursives
                << ", " << ModelBlock->Block_List[j].Size << ", " << nze << ");\n";
              }
            output << "  end;\n";
          }

        output << "  g2=0;g3=0;\n";
        if (ModelBlock->Block_List[j].Temporary_InUse->size())
          {
            tmp_output.str("");
            for (temporary_terms_inuse_type::const_iterator it = ModelBlock->Block_List[j].Temporary_InUse->begin();
                 it != ModelBlock->Block_List[j].Temporary_InUse->end(); it++)
              tmp_output << " T" << *it;
            output << "  global" << tmp_output.str() << ";\n";
          }
        if (ModelBlock->Block_List[j].Simulation_Type!=EVALUATE_BACKWARD && ModelBlock->Block_List[j].Simulation_Type!=EVALUATE_FORWARD)
          output << "  residual=zeros(" << ModelBlock->Block_List[j].Size-ModelBlock->Block_List[j].Nb_Recursives << ",1);\n";
        if (ModelBlock->Block_List[j].Simulation_Type==EVALUATE_BACKWARD)
          output << "  for it_ = (y_kmin+periods):y_kmin+1\n";
        if (ModelBlock->Block_List[j].Simulation_Type==EVALUATE_FORWARD)
          output << "  for it_ = y_kmin+1:(y_kmin+periods)\n";

        if (ModelBlock->Block_List[j].Simulation_Type==SOLVE_TWO_BOUNDARIES_COMPLETE || ModelBlock->Block_List[j].Simulation_Type==SOLVE_TWO_BOUNDARIES_SIMPLE)
          {
            output << "  b = zeros(periods*y_size,1);\n";
            output << "  for it_ = y_kmin+1:(periods+y_kmin)\n";
            output << "    Per_y_=it_*y_size;\n";
            output << "    Per_J_=(it_-y_kmin-1)*y_size;\n";
            output << "    Per_K_=(it_-1)*y_size;\n";
            sps="  ";
          }
        else
          if (ModelBlock->Block_List[j].Simulation_Type==EVALUATE_BACKWARD || ModelBlock->Block_List[j].Simulation_Type==EVALUATE_FORWARD )
            sps = "  ";
          else
            sps="";
        // The equations
        for (i = 0;i < ModelBlock->Block_List[j].Size;i++)
          {
            temporary_terms_type tt2;
            tt2.clear();
            if (ModelBlock->Block_List[j].Temporary_Terms_in_Equation[i]->size())
              output << "  " << sps << "% //Temporary variables" << endl;
            for (temporary_terms_type::const_iterator it = ModelBlock->Block_List[j].Temporary_Terms_in_Equation[i]->begin();
                 it != ModelBlock->Block_List[j].Temporary_Terms_in_Equation[i]->end(); it++)
              {
                output << "  " <<  sps;
                (*it)->writeOutput(output, oMatlabDynamicModelSparse, temporary_terms);
                output << " = ";
                (*it)->writeOutput(output, oMatlabDynamicModelSparse, tt2);
                // Insert current node into tt2
                tt2.insert(*it);
                output << ";" << endl;
              }
            string sModel = symbol_table.getName(symbol_table.getID(eEndogenous, ModelBlock->Block_List[j].Variable[i])) ;
            eq_node = equations[ModelBlock->Block_List[j].Equation[i]];
            lhs = eq_node->get_arg1();
            rhs = eq_node->get_arg2();
            tmp_output.str("");
            /*if((ModelBlock->Block_List[j].Simulation_Type!=EVALUATE_BACKWARD or ModelBlock->Block_List[j].Simulation_Type!=EVALUATE_FORWARD) and (i<ModelBlock->Block_List[j].Nb_Recursives))
              lhs->writeOutput(tmp_output, oMatlabDynamicModelSparse, temporary_terms);
            else*/
						lhs->writeOutput(tmp_output, oMatlabDynamicModelSparse, temporary_terms);
            switch (ModelBlock->Block_List[j].Simulation_Type)
              {
              case EVALUATE_BACKWARD:
              case EVALUATE_FORWARD:
evaluation:     if (ModelBlock->Block_List[j].Simulation_Type==SOLVE_TWO_BOUNDARIES_COMPLETE || ModelBlock->Block_List[j].Simulation_Type==SOLVE_TWO_BOUNDARIES_SIMPLE)
                  output << "    % equation " << ModelBlock->Block_List[j].Equation[i]+1 << " variable : " << sModel
                  << " (" << ModelBlock->Block_List[j].Variable[i]+1 << ") " << block_triangular.c_Equation_Type(ModelBlock->Block_List[j].Equation_Type[i]) << endl;
                output << "    ";
                if (ModelBlock->Block_List[j].Equation_Type[i] == E_EVALUATE)
                  {
                    output << tmp_output.str();
                    output << " = ";
                    /*if(!(ModelBlock->Block_List[j].Simulation_Type==EVALUATE_BACKWARD or ModelBlock->Block_List[j].Simulation_Type==EVALUATE_FORWARD))
                      {
                        lhs->writeOutput(output, oMatlabDynamicModelSparse, temporary_terms);
                        output << "-relax*(";
                        lhs->writeOutput(output, oMatlabDynamicModelSparse, temporary_terms);
                        output << "-(";
                        rhs->writeOutput(output, oMatlabDynamicModelSparse, temporary_terms);
                        output << "))";
                      }
                    else*/
                    rhs->writeOutput(output, oMatlabDynamicModelSparse, temporary_terms);
                  }
                else if (ModelBlock->Block_List[j].Equation_Type[i] == E_EVALUATE_S)
                  {
                    output << "%" << tmp_output.str();
                    output << " = ";
                    if (ModelBlock->Block_List[j].Equation_Normalized[i])
                      {
                        rhs->writeOutput(output, oMatlabDynamicModelSparse, temporary_terms);
                        output << "\n    ";
                        tmp_output.str("");
                        eq_node = (BinaryOpNode *)ModelBlock->Block_List[j].Equation_Normalized[i];
                        lhs = eq_node->get_arg1();
                        rhs = eq_node->get_arg2();
                        lhs->writeOutput(output, oMatlabDynamicModelSparse, temporary_terms);
                        output << " = ";
                        /*if(!(ModelBlock->Block_List[j].Simulation_Type==EVALUATE_BACKWARD or ModelBlock->Block_List[j].Simulation_Type==EVALUATE_FORWARD))
                          {
                            lhs->writeOutput(output, oMatlabDynamicModelSparse, temporary_terms);
                            output << "-relax*(";
                            lhs->writeOutput(output, oMatlabDynamicModelSparse, temporary_terms);
                            output << "-(";
                            rhs->writeOutput(output, oMatlabDynamicModelSparse, temporary_terms);
                            output << "))";
                          }
                        else*/
                          rhs->writeOutput(output, oMatlabDynamicModelSparse, temporary_terms);
                      }
                  }
                else
                  {
                    cerr << "Type missmatch for equation " << ModelBlock->Block_List[j].Equation[i]+1  << "\n";
                    exit(EXIT_FAILURE);
                  }
                output << ";\n";
                break;
              case SOLVE_BACKWARD_SIMPLE:
              case SOLVE_FORWARD_SIMPLE:
              case SOLVE_BACKWARD_COMPLETE:
              case SOLVE_FORWARD_COMPLETE:
                if (i<ModelBlock->Block_List[j].Nb_Recursives)
                  {
                    /*if (ModelBlock->Block_List[j].Equation_Type[i] == E_EVALUATE_S)
                      recursive_variables[getDerivID(symbol_table.getID(eEndogenous, ModelBlock->Block_List[j].Variable[i]), 0)] = ModelBlock->Block_List[j].Equation_Normalized[i];
                    else
                      recursive_variables[getDerivID(symbol_table.getID(eEndogenous, ModelBlock->Block_List[j].Variable[i]), 0)] = equations[ModelBlock->Block_List[j].Equation[i]];*/
                    goto evaluation;
                  }
                feedback_variables.push_back(ModelBlock->Block_List[j].Variable[i]);
                output << "  % equation " << ModelBlock->Block_List[j].Equation[i]+1 << " variable : " << sModel
                << " (" << ModelBlock->Block_List[j].Variable[i]+1 << ") " << block_triangular.c_Equation_Type(ModelBlock->Block_List[j].Equation_Type[i]) << endl;
                output << "  " << "residual(" << i+1-ModelBlock->Block_List[j].Nb_Recursives << ") = (";
                goto end;
              case SOLVE_TWO_BOUNDARIES_COMPLETE:
              case SOLVE_TWO_BOUNDARIES_SIMPLE:
                if (i<ModelBlock->Block_List[j].Nb_Recursives)
                  {
                    /*if (ModelBlock->Block_List[j].Equation_Type[i] == E_EVALUATE_S)
                      recursive_variables[getDerivID(symbol_table.getID(eEndogenous, ModelBlock->Block_List[j].Variable[i]), 0)] = ModelBlock->Block_List[j].Equation_Normalized[i];
                    else
                      recursive_variables[getDerivID(symbol_table.getID(eEndogenous, ModelBlock->Block_List[j].Variable[i]), 0)] = equations[ModelBlock->Block_List[j].Equation[i]];*/
                    goto evaluation;
                  }
                feedback_variables.push_back(ModelBlock->Block_List[j].Variable[i]);
                output << "    % equation " << ModelBlock->Block_List[j].Equation[i]+1 << " variable : " << sModel
                << " (" << ModelBlock->Block_List[j].Variable[i]+1 << ") " << block_triangular.c_Equation_Type(ModelBlock->Block_List[j].Equation_Type[i]) << endl;
                Uf[ModelBlock->Block_List[j].Equation[i]] << "    b(" << i+1-ModelBlock->Block_List[j].Nb_Recursives << "+Per_J_) = -residual(" << i+1-ModelBlock->Block_List[j].Nb_Recursives << ", it_)";
                output << "    residual(" << i+1-ModelBlock->Block_List[j].Nb_Recursives << ", it_) = (";
                goto end;
              default:
end:
                output << tmp_output.str();
                output << ") - (";
                rhs->writeOutput(output, oMatlabDynamicModelSparse, temporary_terms);
                output << ");\n";
#ifdef CONDITION
                if (ModelBlock->Block_List[j].Simulation_Type==SOLVE_TWO_BOUNDARIES_COMPLETE || ModelBlock->Block_List[j].Simulation_Type==SOLVE_TWO_BOUNDARIES_SIMPLE)
                  output << "  condition(" << i+1 << ")=0;\n";
#endif
              }
          }
        // The Jacobian if we have to solve the block
        if (ModelBlock->Block_List[j].Simulation_Type==SOLVE_TWO_BOUNDARIES_SIMPLE
            ||  ModelBlock->Block_List[j].Simulation_Type==SOLVE_TWO_BOUNDARIES_COMPLETE)
          output << "  " << sps << "% Jacobian  " << endl;
        else
          if (ModelBlock->Block_List[j].Simulation_Type==SOLVE_BACKWARD_SIMPLE || ModelBlock->Block_List[j].Simulation_Type==SOLVE_FORWARD_SIMPLE ||
              ModelBlock->Block_List[j].Simulation_Type==SOLVE_BACKWARD_COMPLETE || ModelBlock->Block_List[j].Simulation_Type==SOLVE_FORWARD_COMPLETE)
            output << "  % Jacobian  " << endl << "  if jacobian_eval" << endl;
          else
            output << "    % Jacobian  " << endl << "    if jacobian_eval" << endl;
        switch (ModelBlock->Block_List[j].Simulation_Type)
          {
          case EVALUATE_BACKWARD:
          case EVALUATE_FORWARD:
            for (m=0;m<ModelBlock->Block_List[j].Max_Lead+ModelBlock->Block_List[j].Max_Lag+1;m++)
              {
                k=m-ModelBlock->Block_List[j].Max_Lag;
                for (i=0;i<ModelBlock->Block_List[j].IM_lead_lag[m].size;i++)
                  {
                    int eq=ModelBlock->Block_List[j].IM_lead_lag[m].Equ_Index[i];
                    int var=ModelBlock->Block_List[j].IM_lead_lag[m].Var_Index[i];
                    int eqr=ModelBlock->Block_List[j].IM_lead_lag[m].Equ[i];
                    int varr=ModelBlock->Block_List[j].IM_lead_lag[m].Var[i];
                    output << "      g1(" << eqr+1 << ", " << /*varr+1+(m+variable_table.max_lag-ModelBlock->Block_List[j].Max_Lag)*symbol_table.endo_nbr*/
                    varr+1+m*ModelBlock->Block_List[j].Size << ") = ";
                    writeDerivative(output, eq, symbol_table.getID(eEndogenous, var), k, oMatlabDynamicModelSparse, temporary_terms);
                    output << "; % variable=" << symbol_table.getName(symbol_table.getID(eEndogenous, var))
                    << "(" << k//variable_table.getLag(variable_table.getSymbolID(ModelBlock->Block_List[j].Variable[0]))
                    << ") " << var+1
                    << ", equation=" << eq+1 << endl;
                  }
              }
            //jacobian_max_endo_col=(variable_table.max_endo_lag+variable_table.max_endo_lead+1)*symbol_table.endo_nbr;
            /*for (m=0;m<=ModelBlock->Block_List[j].Max_Lead+ModelBlock->Block_List[j].Max_Lag;m++)
              {
                k=m-ModelBlock->Block_List[j].Max_Lag;
                for (i=0;i<ModelBlock->Block_List[j].IM_lead_lag[m].size_exo;i++)
                  {
                    int eq=ModelBlock->Block_List[j].IM_lead_lag[m].Equ_X_Index[i];
                    int var=ModelBlock->Block_List[j].IM_lead_lag[m].Exogenous_Index[i];
                    int eqr=ModelBlock->Block_List[j].IM_lead_lag[m].Equ_X[i];
                    int varr=ModelBlock->Block_List[j].IM_lead_lag[m].Exogenous[i];
                    output << "      g1_x(" << eqr+1 << ", "
                           << varr+1+(m+max_exo_lag-ModelBlock->Block_List[j].Max_Lag)*symbol_table.exo_nbr() << ") = ";
                    writeDerivative(output, eq, symbol_table.getID(eExogenous, var), k, oMatlabDynamicModelSparse, temporary_terms);
                    output << "; % variable=" << symbol_table.getName(var)
                           << "(" << k << ") " << var+1
                           << ", equation=" << eq+1 << endl;
                  }
              }*/
            for (m=0;m<=ModelBlock->Block_List[j].Max_Lead+ModelBlock->Block_List[j].Max_Lag;m++)
              {
                k=m-ModelBlock->Block_List[j].Max_Lag;
                if (block_triangular.incidencematrix.Model_Max_Lag_Endo - ModelBlock->Block_List[j].Max_Lag +m >=0)
                  {
                    for (i=0;i<ModelBlock->Block_List[j].IM_lead_lag[m].size_other_endo;i++)
                      {
                        int eq=ModelBlock->Block_List[j].IM_lead_lag[m].Equ_Index_other_endo[i];
                        int var=ModelBlock->Block_List[j].IM_lead_lag[m].Var_Index_other_endo[i];
                        int eqr=ModelBlock->Block_List[j].IM_lead_lag[m].Equ_other_endo[i];
                        int varr=ModelBlock->Block_List[j].IM_lead_lag[m].Var_other_endo[i];
                        output << "      g1_o(" << eqr+1 << ", "
                        << varr+1+(m+max_endo_lag-ModelBlock->Block_List[j].Max_Lag)*symbol_table.endo_nbr() << ") = ";
                        writeDerivative(output, eq, symbol_table.getID(eEndogenous, var), k, oMatlabDynamicModelSparse, temporary_terms);
                        output << "; % variable=" << symbol_table.getName(symbol_table.getID(eEndogenous, var))
                        << "(" << k << ") " << var+1
                        << ", equation=" << eq+1 << endl;
                      }
                  }
              }
            output << "      varargout{1}=g1_x;\n";
            output << "      varargout{2}=g1_o;\n";
            output << "    end;" << endl;
            //output << "    ya = y;\n";
            output << "  end;" << endl;
            break;
          case SOLVE_BACKWARD_SIMPLE:
          case SOLVE_FORWARD_SIMPLE:
          case SOLVE_BACKWARD_COMPLETE:
          case SOLVE_FORWARD_COMPLETE:
            for (m=0;m<ModelBlock->Block_List[j].Max_Lead+ModelBlock->Block_List[j].Max_Lag+1;m++)
              {
                k=m-ModelBlock->Block_List[j].Max_Lag;
                for (i=0;i<ModelBlock->Block_List[j].IM_lead_lag[m].size;i++)
                  {
                    int eqr=ModelBlock->Block_List[j].IM_lead_lag[m].Equ_Index[i];
                    int varr=ModelBlock->Block_List[j].IM_lead_lag[m].Var_Index[i];
                    int eq=ModelBlock->Block_List[j].IM_lead_lag[m].Equ[i];
                    int var=ModelBlock->Block_List[j].IM_lead_lag[m].Var[i];
                    output << "    g1(" << eq+1 << ", "
                    << var+1 + m*(ModelBlock->Block_List[j].Size) << ") = ";
                    writeDerivative(output, eqr, symbol_table.getID(eEndogenous, varr), k, oMatlabDynamicModelSparse, temporary_terms);
                    output << "; % variable=" << symbol_table.getName(symbol_table.getID(eEndogenous, varr))
                    << "(" << k << ") " << varr+1
                    << ", equation=" << eqr+1 << endl;
                  }
              }
            /*for (m=0;m<=ModelBlock->Block_List[j].Max_Lead+ModelBlock->Block_List[j].Max_Lag;m++)
              {
                k=m-ModelBlock->Block_List[j].Max_Lag;
                for (i=0;i<ModelBlock->Block_List[j].IM_lead_lag[m].size_exo;i++)
                  {
                    int eq=ModelBlock->Block_List[j].IM_lead_lag[m].Equ_X_Index[i];
                    int var=ModelBlock->Block_List[j].IM_lead_lag[m].Exogenous_Index[i];
                    int eqr=ModelBlock->Block_List[j].IM_lead_lag[m].Equ_X[i];
                    int varr=ModelBlock->Block_List[j].IM_lead_lag[m].Exogenous[i];
                    output << "    g1_x(" << eqr+1 << ", " << varr+1+(m+max_exo_lag-ModelBlock->Block_List[j].Max_Lag)*ModelBlock->Block_List[j].nb_exo << ") = ";
                    writeDerivative(output, eq, symbol_table.getID(eExogenous, var), k, oMatlabDynamicModelSparse, temporary_terms);
                    output << "; % variable=" << symbol_table.getName(var)
                           << "(" << k << ") " << var+1
                           << ", equation=" << eq+1 << endl;
                  }
              }*/
            for (m=0;m<=ModelBlock->Block_List[j].Max_Lead+ModelBlock->Block_List[j].Max_Lag;m++)
              {
                k=m-ModelBlock->Block_List[j].Max_Lag;
                if (block_triangular.incidencematrix.Model_Max_Lag_Endo - ModelBlock->Block_List[j].Max_Lag +m >=0)
                  {
                    for (i=0;i<ModelBlock->Block_List[j].IM_lead_lag[m].size_other_endo;i++)
                      {
                        int eq=ModelBlock->Block_List[j].IM_lead_lag[m].Equ_Index_other_endo[i];
                        int var=ModelBlock->Block_List[j].IM_lead_lag[m].Var_Index_other_endo[i];
                        int eqr=ModelBlock->Block_List[j].IM_lead_lag[m].Equ_other_endo[i];
                        int varr=ModelBlock->Block_List[j].IM_lead_lag[m].Var_other_endo[i];
                        output << "    g1_o(" << eqr+1/*-ModelBlock->Block_List[j].Nb_Recursives*/ << ", "
                        << varr+1+(m+max_endo_lag-ModelBlock->Block_List[j].Max_Lag)*symbol_table.endo_nbr() << ") = ";
                        writeDerivative(output, eq, symbol_table.getID(eEndogenous, var), k, oMatlabDynamicModelSparse, temporary_terms);
                        output << "; % variable=" << symbol_table.getName(symbol_table.getID(eEndogenous, var))
                        << "(" << k << ") " << var+1
                        << ", equation=" << eq+1 << endl;
                      }
                  }
              }
            output << "    varargout{1}=g1_x;\n";
            output << "    varargout{2}=g1_o;\n";
            output << "  else" << endl;

            for(i=0; i<(int)ModelBlock->Block_List[j].Chain_Rule_Derivatives->size();i++)
              {
                pair< pair<int, pair<int, int> >, pair<int, int> > it = ModelBlock->Block_List[j].Chain_Rule_Derivatives->at(i);
                k=it.first.first;
                int eq=it.first.second.first;
                int var=it.first.second.second;
                int eqr=it.second.first;
                int varr=it.second.second;
                output << "    g1(" << eq+1-ModelBlock->Block_List[j].Nb_Recursives << ", "
                       << var+1-ModelBlock->Block_List[j].Nb_Recursives  << ") = ";
                writeChainRuleDerivative(output, eqr, varr, k, oMatlabDynamicModelSparse, temporary_terms);
                output << "; %2 variable=" << symbol_table.getName(symbol_table.getID(eEndogenous, varr))
                       << "(" << k << ") " << varr+1 << ", equation=" << eqr+1 << endl;
              }
            output << "  end;\n";
            break;
          case SOLVE_TWO_BOUNDARIES_SIMPLE:
          case SOLVE_TWO_BOUNDARIES_COMPLETE:
            output << "    if ~jacobian_eval" << endl;
            for(i=0; i<(int)ModelBlock->Block_List[j].Chain_Rule_Derivatives->size();i++)
              {
                pair< pair<int, pair<int, int> >, pair<int, int> > it = ModelBlock->Block_List[j].Chain_Rule_Derivatives->at(i);
                k=it.first.first;
                int eq=it.first.second.first;
                int var=it.first.second.second;
                int eqr=it.second.first;
                int varr=it.second.second;
                ostringstream tmp_output;
                if(eq>=ModelBlock->Block_List[j].Nb_Recursives and var>=ModelBlock->Block_List[j].Nb_Recursives)
                  {
                    if (k==0)
                      Uf[ModelBlock->Block_List[j].Equation[eq]] << "+g1(" << eq+1-ModelBlock->Block_List[j].Nb_Recursives
                        << "+Per_J_, " << var+1-ModelBlock->Block_List[j].Nb_Recursives
                        << "+Per_K_)*y(it_, " << varr+1 << ")";
                    else if (k==1)
                      Uf[ModelBlock->Block_List[j].Equation[eq]] << "+g1(" << eq+1-ModelBlock->Block_List[j].Nb_Recursives
                        << "+Per_J_, " << var+1-ModelBlock->Block_List[j].Nb_Recursives
                        << "+Per_y_)*y(it_+1, " << varr+1 << ")";
                    else if (k>0)
                      Uf[ModelBlock->Block_List[j].Equation[eq]] << "+g1(" << eq+1-ModelBlock->Block_List[j].Nb_Recursives
                        << "+Per_J_, " << var+1-ModelBlock->Block_List[j].Nb_Recursives
                        << "+y_size*(it_+" << k-1 << "))*y(it_+" << k << ", " << varr+1 << ")";
                    else if (k<0)
                      Uf[ModelBlock->Block_List[j].Equation[eq]] << "+g1(" << eq+1-ModelBlock->Block_List[j].Nb_Recursives
                        << "+Per_J_, " << var+1-ModelBlock->Block_List[j].Nb_Recursives
                        << "+y_size*(it_" << k-1 << "))*y(it_" << k << ", " << varr+1 << ")";
                    if (k==0)
                      tmp_output << "     g1(" << eq+1-ModelBlock->Block_List[j].Nb_Recursives << "+Per_J_, "
                        << var+1-ModelBlock->Block_List[j].Nb_Recursives << "+Per_K_) = ";
                    else if (k==1)
                      tmp_output << "     g1(" << eq+1-ModelBlock->Block_List[j].Nb_Recursives << "+Per_J_, "
                        << var+1-ModelBlock->Block_List[j].Nb_Recursives << "+Per_y_) = ";
                    else if (k>0)
                      tmp_output << "     g1(" << eq+1-ModelBlock->Block_List[j].Nb_Recursives << "+Per_J_, "
                        << var+1-ModelBlock->Block_List[j].Nb_Recursives << "+y_size*(it_+" << k-1 << ")) = ";
                    else if (k<0)
                      tmp_output << "     g1(" << eq+1-ModelBlock->Block_List[j].Nb_Recursives << "+Per_J_, "
                        << var+1-ModelBlock->Block_List[j].Nb_Recursives << "+y_size*(it_" << k-1 << ")) = ";
                    output << " " << tmp_output.str();

                    writeChainRuleDerivative(output, eqr, varr, k, oMatlabDynamicModelSparse, temporary_terms);

                    output << ";";
                    output << " %2 variable=" << symbol_table.getName(symbol_table.getID(eEndogenous, varr))
                      << "(" << k << ") " << varr+1
                      << ", equation=" << eqr+1 << " (" << eq+1 << ")" << endl;
                  }
#ifdef CONDITION
                output << "  if (fabs(condition[" << eqr << "])<fabs(u[" << u << "+Per_u_]))\n";
                output << "    condition(" << eqr << ")=u(" << u << "+Per_u_);\n";
#endif
                  //}
              }
            for (i = 0;i < ModelBlock->Block_List[j].Size;i++)
              {
                if (i>=ModelBlock->Block_List[j].Nb_Recursives)
                  output << "  " << Uf[ModelBlock->Block_List[j].Equation[i]].str() << ";\n";
#ifdef CONDITION
                output << "  if (fabs(condition(" << i+1 << "))<fabs(u(" << i << "+Per_u_)))\n";
                output << "    condition(" << i+1 << ")=u(" << i+1 << "+Per_u_);\n";
#endif
              }
#ifdef CONDITION
            for (m=0;m<=ModelBlock->Block_List[j].Max_Lead+ModelBlock->Block_List[j].Max_Lag;m++)
              {
                k=m-ModelBlock->Block_List[j].Max_Lag;
                for (i=0;i<ModelBlock->Block_List[j].IM_lead_lag[m].size;i++)
                  {
                    int eq=ModelBlock->Block_List[j].IM_lead_lag[m].Equ_Index[i];
                    int var=ModelBlock->Block_List[j].IM_lead_lag[m].Var_Index[i];
                    int u=ModelBlock->Block_List[j].IM_lead_lag[m].u[i];
                    int eqr=ModelBlock->Block_List[j].IM_lead_lag[m].Equ[i];
                    output << "  u(" << u+1 << "+Per_u_) = u(" << u+1 << "+Per_u_) / condition(" << eqr+1 << ");\n";
                  }
              }
            for (i = 0;i < ModelBlock->Block_List[j].Size;i++)
              output << "  u(" << i+1 << "+Per_u_) = u(" << i+1 << "+Per_u_) / condition(" << i+1 << ");\n";
#endif

            output << "    else" << endl;
            for (m=0;m<=ModelBlock->Block_List[j].Max_Lead+ModelBlock->Block_List[j].Max_Lag;m++)
              {
                k=m-ModelBlock->Block_List[j].Max_Lag;
                for (i=0;i<ModelBlock->Block_List[j].IM_lead_lag[m].size;i++)
                  {
                    int eq=ModelBlock->Block_List[j].IM_lead_lag[m].Equ_Index[i];
                    int var=ModelBlock->Block_List[j].IM_lead_lag[m].Var_Index[i];
                    int eqr=ModelBlock->Block_List[j].IM_lead_lag[m].Equ[i];
                    int varr=ModelBlock->Block_List[j].IM_lead_lag[m].Var[i];
                    output << "      g1(" << eqr+1 << ", " << varr+1+(m-ModelBlock->Block_List[j].Max_Lag+ModelBlock->Block_List[j].Max_Lag_Endo)*ModelBlock->Block_List[j].Size << ") = ";
                    writeDerivative(output, eq, symbol_table.getID(eEndogenous, var), k, oMatlabDynamicModelSparse, temporary_terms);
                    output << "; % variable=" << symbol_table.getName(symbol_table.getID(eEndogenous, var))
                    << "(" << k << ") " << var+1
                    << ", equation=" << eq+1 << endl;
                  }
              }
            jacobian_max_endo_col=(ModelBlock->Block_List[j].Max_Lead_Endo+ModelBlock->Block_List[j].Max_Lag_Endo+1)*ModelBlock->Block_List[j].Size;
            /*for (m=0;m<=ModelBlock->Block_List[j].Max_Lead+ModelBlock->Block_List[j].Max_Lag;m++)
              {
                k=m-ModelBlock->Block_List[j].Max_Lag;
                for (i=0;i<ModelBlock->Block_List[j].IM_lead_lag[m].size_exo;i++)
                  {
                    int eq=ModelBlock->Block_List[j].IM_lead_lag[m].Equ_X_Index[i];
                    int eqr=ModelBlock->Block_List[j].IM_lead_lag[m].Equ_X[i];
                    int varr=ModelBlock->Block_List[j].IM_lead_lag[m].Exogenous[i];
                    int var=ModelBlock->Block_List[j].IM_lead_lag[m].Exogenous_Index[i];
                    output << "      g1_x(" << eqr+1 << ", "
                           << jacobian_max_endo_col+(m-(ModelBlock->Block_List[j].Max_Lag-ModelBlock->Block_List[j].Max_Lag_Exo))*ModelBlock->Block_List[j].nb_exo+varr+1 << ") = ";
                    writeDerivative(output, eq, symbol_table.getID(eExogenous, var), k, oMatlabDynamicModelSparse, temporary_terms);
                    output << "; % variable (exogenous)=" << symbol_table.getName(var)
                           << "(" << k << ") " << var+1 << " " << varr+1
                           << ", equation=" << eq+1 << endl;
                  }
              }*/
            for (m=0;m<=ModelBlock->Block_List[j].Max_Lead+ModelBlock->Block_List[j].Max_Lag;m++)
              {
                k=m-ModelBlock->Block_List[j].Max_Lag;
                if (block_triangular.incidencematrix.Model_Max_Lag_Endo - ModelBlock->Block_List[j].Max_Lag +m >=0)
                  {
                    for (i=0;i<ModelBlock->Block_List[j].IM_lead_lag[m].size_other_endo;i++)
                      {
                        int eq=ModelBlock->Block_List[j].IM_lead_lag[m].Equ_Index_other_endo[i];
                        int var=ModelBlock->Block_List[j].IM_lead_lag[m].Var_Index_other_endo[i];
                        int eqr=ModelBlock->Block_List[j].IM_lead_lag[m].Equ_other_endo[i];
                        int varr=ModelBlock->Block_List[j].IM_lead_lag[m].Var_other_endo[i];
                        output << "      g1_o(" << eqr+1 << ", "
                        << varr+1+(m+max_endo_lag-ModelBlock->Block_List[j].Max_Lag)*symbol_table.endo_nbr() << ") = ";
                        writeDerivative(output, eq, symbol_table.getID(eEndogenous, var), k, oMatlabDynamicModelSparse, temporary_terms);
                        output << "; % variable=" << symbol_table.getName(symbol_table.getID(eEndogenous, var))
                        << "(" << k << ") " << var+1
                        << ", equation=" << eq+1 << endl;
                      }
                  }
              }
            output << "      varargout{1}=g1_x;\n";
            output << "      varargout{2}=g1_o;\n";
            output << "    end;\n";
            //output << "    ya = y;\n";
            output << "  end;\n";
            break;
          default:
            break;
          }
        output.close();
      }
  }

void
DynamicModel::writeModelEquationsCodeOrdered(const string file_name, const Model_Block *ModelBlock, const string bin_basename, map_idx_type map_idx) const
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

    int i,j,k,v;
    string tmp_s;
    ostringstream tmp_output;
    ofstream code_file;
    NodeID lhs=NULL, rhs=NULL;
    BinaryOpNode *eq_node;
    Uff Uf[symbol_table.endo_nbr()];
    map<NodeID, int> reference_count;
    vector<int> feedback_variables;
    bool file_open=false;
    string main_name=file_name;
    main_name+=".cod";
    code_file.open(main_name.c_str(), ios::out | ios::binary | ios::ate );
    if (!code_file.is_open())
      {
        cout << "Error : Can't open file \"" << main_name << "\" for writing\n";
        exit(EXIT_FAILURE);
      }
    //Temporary variables declaration
    /*code_file.write(&FDIMT, sizeof(FDIMT));
    k=temporary_terms.size();
    code_file.write(reinterpret_cast<char *>(&k),sizeof(k));*/
    FDIMT_ fdimt(temporary_terms.size());
    fdimt.write(code_file);

    for (j = 0; j < ModelBlock->Size ;j++)
      {
        feedback_variables.clear();
        if (j>0)
          {
            FENDBLOCK_ fendblock;
            fendblock.write(code_file);
            //code_file.write(&FENDBLOCK, sizeof(FENDBLOCK));
          }
        int count_u;
        int u_count_int=0;
        if (ModelBlock->Block_List[j].Simulation_Type==SOLVE_TWO_BOUNDARIES_SIMPLE || ModelBlock->Block_List[j].Simulation_Type==SOLVE_TWO_BOUNDARIES_COMPLETE ||
            ModelBlock->Block_List[j].Simulation_Type==SOLVE_BACKWARD_COMPLETE || ModelBlock->Block_List[j].Simulation_Type==SOLVE_FORWARD_COMPLETE)
          {
            //cout << "ModelBlock->Block_List[j].Nb_Recursives = " << ModelBlock->Block_List[j].Nb_Recursives << "\n";
            Write_Inf_To_Bin_File(file_name, bin_basename, j, u_count_int,file_open,
                                  ModelBlock->Block_List[j].Simulation_Type==SOLVE_TWO_BOUNDARIES_COMPLETE || ModelBlock->Block_List[j].Simulation_Type==SOLVE_TWO_BOUNDARIES_SIMPLE);
            //cout << "u_count_int=" << u_count_int << "\n";


            /*code_file.write(reinterpret_cast<char *>(&ModelBlock->Block_List[j].is_linear),sizeof(ModelBlock->Block_List[j].is_linear));
            //v=block_triangular.ModelBlock->Block_List[j].IM_lead_lag[block_triangular.ModelBlock->Block_List[j].Max_Lag + block_triangular.ModelBlock->Block_List[j].Max_Lead].u_finish + 1;
            v=symbol_table.endo_nbr();
            code_file.write(reinterpret_cast<char *>(&v),sizeof(v));
            v=block_triangular.ModelBlock->Block_List[j].Max_Lag;
            code_file.write(reinterpret_cast<char *>(&v),sizeof(v));
            v=block_triangular.ModelBlock->Block_List[j].Max_Lead;
            code_file.write(reinterpret_cast<char *>(&v),sizeof(v));

            v=u_count_int;
            code_file.write(reinterpret_cast<char *>(&v),sizeof(v));*/
            file_open=true;
          }
        FBEGINBLOCK_ fbeginblock(ModelBlock->Block_List[j].Size - ModelBlock->Block_List[j].Nb_Recursives,
                                 ModelBlock->Block_List[j].Simulation_Type,
                                 ModelBlock->Block_List[j].Variable,
                                 ModelBlock->Block_List[j].Equation,
                                 ModelBlock->Block_List[j].Own_Derivative,
                                 ModelBlock->Block_List[j].is_linear,
                                 symbol_table.endo_nbr(),
                                 ModelBlock->Block_List[j].Max_Lag,
                                 ModelBlock->Block_List[j].Max_Lead,
                                 u_count_int
                                 );
        fbeginblock.write(code_file);
        /*code_file.write(&FBEGINBLOCK, sizeof(FBEGINBLOCK));
        v=ModelBlock->Block_List[j].Size - ModelBlock->Block_List[j].Nb_Recursives;
        //cout << "v (Size) = " << v  << "\n";
        code_file.write(reinterpret_cast<char *>(&v),sizeof(v));
        v=ModelBlock->Block_List[j].Simulation_Type;
        code_file.write(reinterpret_cast<char *>(&v),sizeof(v));

        for (i=ModelBlock->Block_List[j].Nb_Recursives; i < ModelBlock->Block_List[j].Size;i++)
          {
            code_file.write(reinterpret_cast<char *>(&ModelBlock->Block_List[j].Variable[i]),sizeof(ModelBlock->Block_List[j].Variable[i]));
            code_file.write(reinterpret_cast<char *>(&ModelBlock->Block_List[j].Equation[i]),sizeof(ModelBlock->Block_List[j].Equation[i]));
            code_file.write(reinterpret_cast<char *>(&ModelBlock->Block_List[j].Own_Derivative[i]),sizeof(ModelBlock->Block_List[j].Own_Derivative[i]));
          }*/

            // The equations
            for (i = 0;i < ModelBlock->Block_List[j].Size;i++)
              {
                //The Temporary terms
                temporary_terms_type tt2;
                tt2.clear();
#ifdef DEBUGC
                k=0;
#endif
                for (temporary_terms_type::const_iterator it = ModelBlock->Block_List[j].Temporary_Terms_in_Equation[i]->begin();
                     it != ModelBlock->Block_List[j].Temporary_Terms_in_Equation[i]->end(); it++)
                  {
                    (*it)->compile(code_file, false, tt2, map_idx, true, false);

                    FSTPT_ fstpt((int)(map_idx.find((*it)->idx)->second));
                    fstpt.write(code_file);

                    /*code_file.write(&FSTPT, sizeof(FSTPT));
                    map_idx_type::const_iterator ii=map_idx.find((*it)->idx);
                    v=(int)ii->second;
                    code_file.write(reinterpret_cast<char *>(&v), sizeof(v));*/

                    // Insert current node into tt2
                    tt2.insert(*it);
#ifdef DEBUGC
                    cout << "FSTPT " << v << "\n";
                    code_file.write(&FOK, sizeof(FOK));
                    code_file.write(reinterpret_cast<char *>(&k), sizeof(k));
                    ki++;
#endif

                  }
#ifdef DEBUGC
                for (temporary_terms_type::const_iterator it = ModelBlock->Block_List[j].Temporary_terms->begin();
                     it != ModelBlock->Block_List[j].Temporary_terms->end(); it++)
                  {
                    map_idx_type::const_iterator ii=map_idx.find((*it)->idx);
                    cout << "map_idx[" << (*it)->idx <<"]=" << ii->second << "\n";
                  }
#endif
                switch (ModelBlock->Block_List[j].Simulation_Type)
                  {
evaluation:
                  case EVALUATE_BACKWARD:
                  case EVALUATE_FORWARD:
                    if (ModelBlock->Block_List[j].Equation_Type[i] == E_EVALUATE)
                      {
                        eq_node = equations[ModelBlock->Block_List[j].Equation[i]];
                        lhs = eq_node->get_arg1();
                        rhs = eq_node->get_arg2();
                        rhs->compile(code_file, false, temporary_terms, map_idx, true, false);
                        lhs->compile(code_file, true, temporary_terms, map_idx, true, false);
                      }
                    else if (ModelBlock->Block_List[j].Equation_Type[i] == E_EVALUATE_S)
                      {
                        eq_node = (BinaryOpNode*)ModelBlock->Block_List[j].Equation_Normalized[i];
                        lhs = eq_node->get_arg1();
                        rhs = eq_node->get_arg2();
                        rhs->compile(code_file, false, temporary_terms, map_idx, true, false);
                        lhs->compile(code_file, true, temporary_terms, map_idx, true, false);
                      }
                    break;
                  case SOLVE_BACKWARD_COMPLETE:
                  case SOLVE_FORWARD_COMPLETE:
                  case SOLVE_TWO_BOUNDARIES_COMPLETE:
                  case SOLVE_TWO_BOUNDARIES_SIMPLE:
                    if (i<ModelBlock->Block_List[j].Nb_Recursives)
                      goto evaluation;
                    feedback_variables.push_back(ModelBlock->Block_List[j].Variable[i]);
                    v=ModelBlock->Block_List[j].Equation[i];
                    Uf[v].Ufl=NULL;
                    goto end;
                  default:
end:
                    eq_node = equations[ModelBlock->Block_List[j].Equation[i]];
                    lhs = eq_node->get_arg1();
                    rhs = eq_node->get_arg2();
                    lhs->compile(code_file, false, temporary_terms, map_idx, true, false);
                    rhs->compile(code_file, false, temporary_terms, map_idx, true, false);

                    FBINARY_ fbinary(oMinus);
                    fbinary.write(code_file);
                    /*code_file.write(&FBINARY, sizeof(FBINARY));
                    int v=oMinus;
                    code_file.write(reinterpret_cast<char *>(&v),sizeof(v));*/
                    FSTPR_ fstpr(i - ModelBlock->Block_List[j].Nb_Recursives);
                    fstpr.write(code_file);
                    /*code_file.write(&FSTPR, sizeof(FSTPR));
                    v = i - ModelBlock->Block_List[j].Nb_Recursives;
                    code_file.write(reinterpret_cast<char *>(&v), sizeof(v));*/
                  }
              }
            FENDEQU_ fendequ;
            fendequ.write(code_file);
            //code_file.write(&FENDEQU, sizeof(FENDEQU));
            // The Jacobian if we have to solve the block
            if (ModelBlock->Block_List[j].Simulation_Type!=EVALUATE_BACKWARD
                && ModelBlock->Block_List[j].Simulation_Type!=EVALUATE_FORWARD)
              {
                switch (ModelBlock->Block_List[j].Simulation_Type)
                  {
                  case SOLVE_BACKWARD_SIMPLE:
                  case SOLVE_FORWARD_SIMPLE:
                    compileDerivative(code_file, ModelBlock->Block_List[j].Equation[0], ModelBlock->Block_List[j].Variable[0], 0, map_idx);
                      {
                        FSTPG_ fstpg(0);
                        fstpg.write(code_file);
                      }
                    /*code_file.write(&FSTPG, sizeof(FSTPG));
                    v=0;
                    code_file.write(reinterpret_cast<char *>(&v), sizeof(v));*/
                    break;

                  case SOLVE_BACKWARD_COMPLETE:
                  case SOLVE_FORWARD_COMPLETE:
                  case SOLVE_TWO_BOUNDARIES_COMPLETE:
                  case SOLVE_TWO_BOUNDARIES_SIMPLE:
                    //count_u=ModelBlock->Block_List[j].Size - ModelBlock->Block_List[j].Nb_Recursives;
                    count_u = feedback_variables.size();
                    for(i=0; i<(int)ModelBlock->Block_List[j].Chain_Rule_Derivatives->size();i++)
                      {
                        pair< pair<int, pair<int, int> >, pair<int, int> > it = ModelBlock->Block_List[j].Chain_Rule_Derivatives->at(i);
                        k=it.first.first;
                        int eq=it.first.second.first;
                        int var=it.first.second.second;
                        int eqr=it.second.first;
                        int varr=it.second.second;
                        //cout << "k=" << k << " eq=" << eq << " (" << eq-ModelBlock->Block_List[j].Nb_Recursives << ") var=" << var << " (" << var-ModelBlock->Block_List[j].Nb_Recursives << ") eqr=" << eqr << " varr=" << varr << " count_u=" << count_u << "\n";
                        int v=ModelBlock->Block_List[j].Equation[eq];
                        /*m = ModelBlock->Block_List[j].Max_Lag + k;
                        int u=ModelBlock->Block_List[j].IM_lead_lag[m].u[i];*/
                        if(eq>=ModelBlock->Block_List[j].Nb_Recursives and var>=ModelBlock->Block_List[j].Nb_Recursives)
                          {
                            if (!Uf[v].Ufl)
                              {
                                Uf[v].Ufl=(Uff_l*)malloc(sizeof(Uff_l));
                                Uf[v].Ufl_First=Uf[v].Ufl;
                              }
                            else
                              {
                                Uf[v].Ufl->pNext=(Uff_l*)malloc(sizeof(Uff_l));
                                Uf[v].Ufl=Uf[v].Ufl->pNext;
                              }
                            Uf[v].Ufl->pNext=NULL;
                            Uf[v].Ufl->u=count_u;
                            Uf[v].Ufl->var=varr;
                            Uf[v].Ufl->lag=k;
                            compileChainRuleDerivative(code_file, eqr, varr, k, map_idx);

                            FSTPU_ fstpu(count_u);
                            fstpu.write(code_file);

                            /*code_file.write(&FSTPU, sizeof(FSTPU));
                            code_file.write(reinterpret_cast<char *>(&count_u), sizeof(count_u));*/
                            count_u++;
												  }
											}
                    for (i = 0;i < ModelBlock->Block_List[j].Size;i++)
                      {
                        if(i>=ModelBlock->Block_List[j].Nb_Recursives)
                          {
                            FLDR_ fldr(i-ModelBlock->Block_List[j].Nb_Recursives);
                            fldr.write(code_file);
                            /*code_file.write(&FLDR, sizeof(FLDR));
                            v = i-ModelBlock->Block_List[j].Nb_Recursives;
                            code_file.write(reinterpret_cast<char *>(&v), sizeof(v));*/

                            FLDZ_ fldz;
                            fldz.write(code_file);
                            //code_file.write(&FLDZ, sizeof(FLDZ));

                            v=ModelBlock->Block_List[j].Equation[i];
                            for (Uf[v].Ufl=Uf[v].Ufl_First; Uf[v].Ufl; Uf[v].Ufl=Uf[v].Ufl->pNext)
                              {
                                FLDU_ fldu(Uf[v].Ufl->u);
                                fldu.write(code_file);
                                /*code_file.write(&FLDU, sizeof(FLDU));
                                code_file.write(reinterpret_cast<char *>(&Uf[v].Ufl->u), sizeof(Uf[v].Ufl->u));*/
                                FLDV_ fldv(eEndogenous, Uf[v].Ufl->var, Uf[v].Ufl->lag);
                                fldv.write(code_file);

                                /*code_file.write(&FLDV, sizeof(FLDV));
                                char vc=eEndogenous;
                                code_file.write(reinterpret_cast<char *>(&vc), sizeof(vc));
                                int v1=Uf[v].Ufl->var;
                                code_file.write(reinterpret_cast<char *>(&v1), sizeof(v1));
                                v1=Uf[v].Ufl->lag;
                                code_file.write(reinterpret_cast<char *>(&v1), sizeof(v1));*/
                                FBINARY_ fbinary(oTimes);
                                fbinary.write(code_file);
                                /*code_file.write(&FBINARY, sizeof(FBINARY));
                                v1=oTimes;
                                code_file.write(reinterpret_cast<char *>(&v1), sizeof(v1));*/

                                FCUML_ fcuml;
                                fcuml.write(code_file);
                                //code_file.write(&FCUML, sizeof(FCUML));
                              }
                            Uf[v].Ufl=Uf[v].Ufl_First;
                            while (Uf[v].Ufl)
                              {
                                Uf[v].Ufl_First=Uf[v].Ufl->pNext;
                                free(Uf[v].Ufl);
                                Uf[v].Ufl=Uf[v].Ufl_First;
                              }
                            FBINARY_ fbinary(oMinus);
                            fbinary.write(code_file);
                            /*code_file.write(&FBINARY, sizeof(FBINARY));
                            v=oMinus;
                            code_file.write(reinterpret_cast<char *>(&v), sizeof(v));*/

                            FSTPU_ fstpu(i - ModelBlock->Block_List[j].Nb_Recursives);
                            fstpu.write(code_file);
                            /*code_file.write(&FSTPU, sizeof(FSTPU));
                            v = i - ModelBlock->Block_List[j].Nb_Recursives;
                            code_file.write(reinterpret_cast<char *>(&v), sizeof(v));*/
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
    //code_file.write(&FENDBLOCK, sizeof(FENDBLOCK));
    FEND_ fend;
    fend.write(code_file);
    //code_file.write(&FEND, sizeof(FEND));
    code_file.close();
  }

void
DynamicModel::writeDynamicMFile(const string &dynamic_basename) const
  {
    string filename = dynamic_basename + ".m";

    ofstream mDynamicModelFile;
    mDynamicModelFile.open(filename.c_str(), ios::out | ios::binary);
    if (!mDynamicModelFile.is_open())
      {
        cerr << "Error: Can't open file " << filename << " for writing" << endl;
        exit(EXIT_FAILURE);
      }
    mDynamicModelFile << "function [residual, g1, g2, g3] = " << dynamic_basename << "(y, x, params, it_)" << endl
    << "%" << endl
    << "% Status : Computes dynamic model for Dynare" << endl
    << "%" << endl
    << "% Warning : this file is generated automatically by Dynare" << endl
    << "%           from model file (.mod)" << endl << endl;
    
    if (containsSteadyStateOperator())
      mDynamicModelFile << "global oo_;" << endl << endl;

    writeDynamicModel(mDynamicModelFile, false);

    mDynamicModelFile.close();
  }

void
DynamicModel::writeDynamicCFile(const string &dynamic_basename) const
  {
    string filename = dynamic_basename + ".c";
    ofstream mDynamicModelFile;
    
    mDynamicModelFile.open(filename.c_str(), ios::out | ios::binary);
    if (!mDynamicModelFile.is_open())
      {
        cerr << "Error: Can't open file " << filename << " for writing" << endl;
        exit(EXIT_FAILURE);
      }
    mDynamicModelFile << "/*" << endl
                      << " * " << filename << " : Computes dynamic model for Dynare" << endl
                      << " *" << endl
                      << " * Warning : this file is generated automatically by Dynare" << endl
                      << " *           from model file (.mod)" << endl
                      << endl
                      << " */" << endl
                      << "#include <math.h>" << endl
                      << "#include \"mex.h\"" << endl
                      << endl
                      << "#define max(a, b) (((a) > (b)) ? (a) : (b))" << endl
                      << "#define min(a, b) (((a) > (b)) ? (b) : (a))" << endl;

    // Writing the function body
    writeDynamicModel(mDynamicModelFile, true);

    // Writing the gateway routine
    mDynamicModelFile << "/* The gateway routine */" << endl
    << "void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])" << endl
    << "{" << endl
    << "  double *y, *x, *params;" << endl
    << "  double *residual, *g1, *v2, *v3;" << endl
    << "  int nb_row_x, it_;" << endl
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
    << "  /* Fetch time index */" << endl
    << "  it_ = (int) mxGetScalar(prhs[3]) - 1;" << endl
    << endl
    << "  /* Gets number of rows of matrix x. */" << endl
    << "  nb_row_x = mxGetM(prhs[1]);" << endl
    << endl
    << "  residual = NULL;" << endl
    << "  if (nlhs >= 1)" << endl
    << "  {" << endl
    << "     /* Set the output pointer to the output matrix residual. */" << endl
    << "     plhs[0] = mxCreateDoubleMatrix(" << equations.size() << ",1, mxREAL);" << endl
    << "     /* Create a C pointer to a copy of the output matrix residual. */" << endl
    << "     residual = mxGetPr(plhs[0]);" << endl
    << "  }" << endl
    << endl
    << "  g1 = NULL;" << endl
    << "  if (nlhs >= 2)" << endl
    << "  {" << endl
    << "     /* Set the output pointer to the output matrix g1. */" << endl

    << "     plhs[1] = mxCreateDoubleMatrix(" << equations.size() << ", " << dynJacobianColsNbr << ", mxREAL);" << endl
    << "     /* Create a C pointer to a copy of the output matrix g1. */" << endl
    << "     g1 = mxGetPr(plhs[1]);" << endl
    << "  }" << endl
    << endl
    << "  v2 = NULL;" << endl
    << " if (nlhs >= 3)" << endl
    << "  {" << endl
    << "     /* Set the output pointer to the output matrix v2. */" << endl
    << "     plhs[2] = mxCreateDoubleMatrix(" << NNZDerivatives[1] << ", " << 3
    << ", mxREAL);" << endl
    << "     v2 = mxGetPr(plhs[2]);" << endl
    << "  }" << endl
    << endl
    << " if (nlhs >= 4)" << endl
    << "  {" << endl
    << "     /* Set the output pointer to the output matrix v3. */" << endl
    << "     plhs[3] = mxCreateDoubleMatrix(" << NNZDerivatives[2] << ", " << 3 << ", mxREAL);" << endl
    << "     v3 = mxGetPr(plhs[3]);" << endl
    << "  }" << endl
    << endl
    << "  /* Call the C subroutines. */" << endl
    << "  Dynamic(y, x, nb_row_x, params, it_, residual, g1, v2, v3);" << endl
    << "}" << endl;
    mDynamicModelFile.close();
  }

string
DynamicModel::reform(const string name1) const
  {
    string name=name1;
    int pos = name.find("\\", 0);
    while (pos >= 0)
      {
        if (name.substr(pos + 1, 1) != "\\")
          {
            name = name.insert(pos, "\\");
            pos++;
          }
        pos++;
        pos = name.find("\\", pos);
      }
    return (name);
  }

void
DynamicModel::Write_Inf_To_Bin_File(const string &dynamic_basename, const string &bin_basename, const int &num,
                                    int &u_count_int, bool &file_open, bool is_two_boundaries) const
  {
    int j;
    std::ofstream SaveCode;
    if (file_open)
      SaveCode.open((bin_basename + "_dynamic.bin").c_str(), ios::out | ios::in | ios::binary | ios ::ate );
    else
      SaveCode.open((bin_basename + "_dynamic.bin").c_str(), ios::out | ios::binary);
    if (!SaveCode.is_open())
      {
        cout << "Error : Can't open file \"" << bin_basename << "_dynamic.bin\" for writing\n";
        exit(EXIT_FAILURE);
      }
    u_count_int=0;
    int Size = block_triangular.ModelBlock->Block_List[num].Size - block_triangular.ModelBlock->Block_List[num].Nb_Recursives;
    for(int i=0; i<(int)block_triangular.ModelBlock->Block_List[num].Chain_Rule_Derivatives->size();i++)
			{
        //Chain_Rule_Derivatives.insert(make_pair( make_pair(eq, eqr), make_pair(var, make_pair(varr, lag))));
        pair< pair<int, pair<int, int> >, pair<int, int> > it = block_triangular.ModelBlock->Block_List[num].Chain_Rule_Derivatives->at(i);
        int k=it.first.first;
        int eq=it.first.second.first;

        int var_init=it.first.second.second;
        /*int eqr=it.second.first;
        int varr=it.second.second;*/
        if(eq>=block_triangular.ModelBlock->Block_List[num].Nb_Recursives and var_init>=block_triangular.ModelBlock->Block_List[num].Nb_Recursives)
					{
            int v=eq-block_triangular.ModelBlock->Block_List[num].Nb_Recursives;
            SaveCode.write(reinterpret_cast<char *>(&v), sizeof(v));
						int var=it.first.second.second-block_triangular.ModelBlock->Block_List[num].Nb_Recursives + k * Size;
				    SaveCode.write(reinterpret_cast<char *>(&var), sizeof(var));
            SaveCode.write(reinterpret_cast<char *>(&k), sizeof(k));
            int u = u_count_int + Size;
            SaveCode.write(reinterpret_cast<char *>(&u), sizeof(u));
            //cout << "eq=" << eq << " var=" << var << " k=" << k << " u=" << u << "\n";
            u_count_int++;
					}
			}


    /*for (int m=0;m<=block_triangular.ModelBlock->Block_List[num].Max_Lead+block_triangular.ModelBlock->Block_List[num].Max_Lag;m++)
      {
        int k1=m-block_triangular.ModelBlock->Block_List[num].Max_Lag;
        for (j=0;j<block_triangular.ModelBlock->Block_List[num].IM_lead_lag[m].size;j++)
          {
            int varr=block_triangular.ModelBlock->Block_List[num].IM_lead_lag[m].Var[j]+k1*block_triangular.ModelBlock->Block_List[num].Size;
            int u=block_triangular.ModelBlock->Block_List[num].IM_lead_lag[m].u[j];
            int eqr1=block_triangular.ModelBlock->Block_List[num].IM_lead_lag[m].Equ[j];
            SaveCode.write(reinterpret_cast<char *>(&eqr1), sizeof(eqr1));
            SaveCode.write(reinterpret_cast<char *>(&varr), sizeof(varr));
            SaveCode.write(reinterpret_cast<char *>(&k1), sizeof(k1));
            SaveCode.write(reinterpret_cast<char *>(&u), sizeof(u));
            u_count_int++;
          }
      }*/
    if (is_two_boundaries)
      {
        /*for (j=0;j<Size;j++)
          {
            int eqr1=j;
            int varr=Size*(block_triangular.periods
                     +block_triangular.incidencematrix.Model_Max_Lead_Endo);
            int k1=0;
            SaveCode.write(reinterpret_cast<char *>(&eqr1), sizeof(eqr1));
            SaveCode.write(reinterpret_cast<char *>(&varr), sizeof(varr));
            SaveCode.write(reinterpret_cast<char *>(&k1), sizeof(k1));
            SaveCode.write(reinterpret_cast<char *>(&eqr1), sizeof(eqr1));
            u_count_int++;
          }*/
				u_count_int+=Size;
      }
    //cout << "u_count_int=" << u_count_int << "\n";
    for (j=block_triangular.ModelBlock->Block_List[num].Nb_Recursives;j<block_triangular.ModelBlock->Block_List[num].Size;j++)
      {
        int varr=block_triangular.ModelBlock->Block_List[num].Variable[j];
        SaveCode.write(reinterpret_cast<char *>(&varr), sizeof(varr));
      }
    for (j=block_triangular.ModelBlock->Block_List[num].Nb_Recursives;j<block_triangular.ModelBlock->Block_List[num].Size;j++)
      {
        int eqr1=block_triangular.ModelBlock->Block_List[num].Equation[j];
        SaveCode.write(reinterpret_cast<char *>(&eqr1), sizeof(eqr1));
      }
    SaveCode.close();
  }

void
DynamicModel::writeSparseDynamicMFile(const string &dynamic_basename, const string &basename) const
  {
    string sp;
    ofstream mDynamicModelFile;
    ostringstream tmp, tmp1, tmp_eq;
    int prev_Simulation_Type, tmp_i;
    //SymbolicGaussElimination SGE;
    bool OK;
    chdir(basename.c_str());
    string filename = dynamic_basename + ".m";
    mDynamicModelFile.open(filename.c_str(), ios::out | ios::binary);
    if (!mDynamicModelFile.is_open())
      {
        cerr << "Error: Can't open file " << filename << " for writing" << endl;
        exit(EXIT_FAILURE);
      }
    mDynamicModelFile << "%\n";
    mDynamicModelFile << "% " << filename << " : Computes dynamic model for Dynare\n";
    mDynamicModelFile << "%\n";
    mDynamicModelFile << "% Warning : this file is generated automatically by Dynare\n";
    mDynamicModelFile << "%           from model file (.mod)\n\n";
    mDynamicModelFile << "%/\n";

    int i, k, Nb_SGE=0;
    bool skip_head, open_par=false;

    mDynamicModelFile << "function [varargout] = " << dynamic_basename << "(varargin)\n";
    mDynamicModelFile << "  global oo_ options_ M_ ;\n";
    mDynamicModelFile << "  g2=[];g3=[];\n";
    //Temporary variables declaration
    OK=true;
    ostringstream tmp_output;
    for (temporary_terms_type::const_iterator it = temporary_terms.begin();
         it != temporary_terms.end(); it++)
      {
        if (OK)
          OK=false;
        else
          tmp_output << " ";
        (*it)->writeOutput(tmp_output, oMatlabStaticModelSparse, temporary_terms);
      }
    if (tmp_output.str().length()>0)
      mDynamicModelFile << "  global " << tmp_output.str() << " M_ ;\n";

    mDynamicModelFile << "  T_init=zeros(1,options_.periods+M_.maximum_lag+M_.maximum_lead);\n";
    tmp_output.str("");
    for (temporary_terms_type::const_iterator it = temporary_terms.begin();
         it != temporary_terms.end(); it++)
      {
        tmp_output << "  ";
        (*it)->writeOutput(tmp_output, oMatlabDynamicModel, temporary_terms);
        tmp_output << "=T_init;\n";
      }
    if (tmp_output.str().length()>0)
      mDynamicModelFile << tmp_output.str();

    mDynamicModelFile << "  y_kmin=M_.maximum_lag;\n";
    mDynamicModelFile << "  y_kmax=M_.maximum_lead;\n";
    mDynamicModelFile << "  y_size=M_.endo_nbr;\n";
    mDynamicModelFile << "  if(length(varargin)>0)\n";
    mDynamicModelFile << "    %it is a simple evaluation of the dynamic model for time _it\n";
    mDynamicModelFile << "    params=varargin{3};\n";
    mDynamicModelFile << "    it_=varargin{4};\n";
    /*i = symbol_table.endo_nbr*(variable_table.max_endo_lag+variable_table.max_endo_lead+1)+
      symbol_table.exo_nbr*(variable_table.max_exo_lag+variable_table.max_exo_lead+1);
      mDynamicModelFile << "    g1=spalloc(" << symbol_table.endo_nbr << ", " << i << ", " << i*symbol_table.endo_nbr << ");\n";*/
    mDynamicModelFile << "    Per_u_=0;\n";
    mDynamicModelFile << "    Per_y_=it_*y_size;\n";
    mDynamicModelFile << "    y=varargin{1};\n";
    mDynamicModelFile << "    ys=y(it_,:);\n";
    mDynamicModelFile << "    x=varargin{2};\n";
    prev_Simulation_Type=-1;
    tmp.str("");
    tmp_eq.str("");
    for (int count_call=1, i = 0;i < block_triangular.ModelBlock->Size;i++, count_call++)
      {
        mDynamicModelFile << "    %block_triangular.ModelBlock->Block_List[i].Nb_Recursives=" << block_triangular.ModelBlock->Block_List[i].Nb_Recursives << " block_triangular.ModelBlock->Block_List[i].Size=" <<  block_triangular.ModelBlock->Block_List[i].Size << "\n";
        k=block_triangular.ModelBlock->Block_List[i].Simulation_Type;
        if(k==EVALUATE_FORWARD || k==EVALUATE_BACKWARD)
          {
            for (int ik=0 ;ik<block_triangular.ModelBlock->Block_List[i].Size;ik++)
              {
                tmp << " " << block_triangular.ModelBlock->Block_List[i].Variable[ik]+1;
                tmp_eq << " " << block_triangular.ModelBlock->Block_List[i].Equation[ik]+1;
              }
          }
        else
          {
            for (int ik=block_triangular.ModelBlock->Block_List[i].Nb_Recursives ;ik<block_triangular.ModelBlock->Block_List[i].Size;ik++)
              {
                tmp << " " << block_triangular.ModelBlock->Block_List[i].Variable[ik]+1;
                tmp_eq << " " << block_triangular.ModelBlock->Block_List[i].Equation[ik]+1;
              }
          }
        mDynamicModelFile << "    y_index_eq=[" << tmp_eq.str() << "];\n";
        mDynamicModelFile << "    y_index=[" << tmp.str() << "];\n";

        switch (k)
          {
          case EVALUATE_FORWARD:
          case EVALUATE_BACKWARD:
            mDynamicModelFile << "    [y, dr(" << count_call << ").g1, dr(" << count_call << ").g2, dr(" << count_call << ").g3, dr(" << count_call << ").g1_x, dr(" << count_call << ").g1_o]=" << dynamic_basename << "_" << i + 1 << "(y, x, params, 1, it_-1, 1);\n";
            mDynamicModelFile << "    residual(y_index_eq)=ys(y_index)-y(it_, y_index);\n";
           break;
          case SOLVE_FORWARD_SIMPLE:
          case SOLVE_BACKWARD_SIMPLE:
            //mDynamicModelFile << "    y_index_eq = " << block_triangular.ModelBlock->Block_List[i].Equation[0]+1 << ";\n";
            mDynamicModelFile << "    [r, dr(" << count_call << ").g1, dr(" << count_call << ").g2, dr(" << count_call << ").g3, dr(" << count_call << ").g1_x, dr(" << count_call << ").g1_o]=" << dynamic_basename << "_" << i + 1 << "(y, x, params, it_, 1);\n";
            mDynamicModelFile << "    residual(y_index_eq)=r;\n";
            break;
          case SOLVE_FORWARD_COMPLETE:
          case SOLVE_BACKWARD_COMPLETE:
            //mDynamicModelFile << "    y_index_eq = [" << tmp_eq.str() << "];\n";
            mDynamicModelFile << "    [r, dr(" << count_call << ").g1, dr(" << count_call << ").g2, dr(" << count_call << ").g3, dr(" << count_call << ").g1_x, dr(" << count_call << ").g1_o]=" << dynamic_basename << "_" << i + 1 << "(y, x, params, it_, 1);\n";
            mDynamicModelFile << "    residual(y_index_eq)=r;\n";
            break;
          case SOLVE_TWO_BOUNDARIES_COMPLETE:
          case SOLVE_TWO_BOUNDARIES_SIMPLE:
            int j;
            /*mDynamicModelFile << "    y_index_eq = [" << tmp_eq.str() << "];\n";
            tmp_i=block_triangular.ModelBlock->Block_List[i].Max_Lag_Endo+block_triangular.ModelBlock->Block_List[i].Max_Lead_Endo+1;
            mDynamicModelFile << "    y_index = [";
            for (j=0;j<tmp_i;j++)
              for (int ik=0;ik<block_triangular.ModelBlock->Block_List[i].Size;ik++)
                {
                  mDynamicModelFile << " " << block_triangular.ModelBlock->Block_List[i].Variable[ik]+1+j*symbol_table.endo_nbr();
                }
            int tmp_ix=block_triangular.ModelBlock->Block_List[i].Max_Lag_Exo+block_triangular.ModelBlock->Block_List[i].Max_Lead_Exo+1;
            for (j=0;j<tmp_ix;j++)
              for (int ik=0;ik<block_triangular.ModelBlock->Block_List[i].nb_exo;ik++)
                mDynamicModelFile << " " << block_triangular.ModelBlock->Block_List[i].Exogenous[ik]+1+j*symbol_table.exo_nbr()+symbol_table.endo_nbr()*tmp_i;
            mDynamicModelFile << " ];\n";*/
            //mDynamicModelFile << "    ga = [];\n";
            j = block_triangular.ModelBlock->Block_List[i].Size*(block_triangular.ModelBlock->Block_List[i].Max_Lag_Endo+block_triangular.ModelBlock->Block_List[i].Max_Lead_Endo+1)
                + block_triangular.ModelBlock->Block_List[i].nb_exo*(block_triangular.ModelBlock->Block_List[i].Max_Lag_Exo+block_triangular.ModelBlock->Block_List[i].Max_Lead_Exo+1);
            /*mDynamicModelFile << "    ga=spalloc(" << block_triangular.ModelBlock->Block_List[i].Size << ", " << j << ", " <<
              block_triangular.ModelBlock->Block_List[i].Size*j << ");\n";*/
            tmp_i=block_triangular.ModelBlock->Block_List[i].Max_Lag_Endo+block_triangular.ModelBlock->Block_List[i].Max_Lead_Endo+1;
            mDynamicModelFile << "    [r, dr(" << count_call << ").g1, dr(" << count_call << ").g2, dr(" << count_call << ").g3, b, dr(" << count_call << ").g1_x, dr(" << count_call << ").g1_o]=" << dynamic_basename << "_" <<  i + 1 << "(y, x, params, it_-" << max_lag << ", 1, " << max_lag << ", " << block_triangular.ModelBlock->Block_List[i].Size-block_triangular.ModelBlock->Block_List[i].Nb_Recursives << ");\n";
            /*if(block_triangular.ModelBlock->Block_List[i].Max_Lag==variable_table.max_lag && block_triangular.ModelBlock->Block_List[i].Max_Lead==variable_table.max_lead)
              mDynamicModelFile << "    g1(y_index_eq,y_index) = ga;\n";
              else
              mDynamicModelFile << "    g1(y_index_eq,y_index) = ga(:," << 1+(variable_table.max_lag-block_triangular.ModelBlock->Block_List[i].Max_Lag)*block_triangular.ModelBlock->Block_List[i].Size << ":" << (variable_table.max_lag+1+block_triangular.ModelBlock->Block_List[i].Max_Lead)*block_triangular.ModelBlock->Block_List[i].Size << ");\n";*/
            mDynamicModelFile << "    residual(y_index_eq)=r(:,M_.maximum_lag+1);\n";
            break;
          }
        tmp_eq.str("");
        tmp.str("");
      }
    if (tmp1.str().length())
      {
        mDynamicModelFile << tmp1.str();
        tmp1.str("");
      }
    mDynamicModelFile << "    varargout{1}=residual;\n";
    mDynamicModelFile << "    varargout{2}=dr;\n";
    mDynamicModelFile << "    return;\n";
    mDynamicModelFile << "  end;\n";
    mDynamicModelFile << "  %it is the deterministic simulation of the block decomposed dynamic model\n";
    mDynamicModelFile << "  if(options_.stack_solve_algo==1)\n";
    mDynamicModelFile << "    mthd='Sparse LU';\n";
    mDynamicModelFile << "  elseif(options_.stack_solve_algo==2)\n";
    mDynamicModelFile << "    mthd='GMRES';\n";
    mDynamicModelFile << "  elseif(options_.stack_solve_algo==3)\n";
    mDynamicModelFile << "    mthd='BICGSTAB';\n";
    mDynamicModelFile << "  elseif(options_.stack_solve_algo==4)\n";
    mDynamicModelFile << "    mthd='OPTIMPATH';\n";
    mDynamicModelFile << "  else\n";
    mDynamicModelFile << "    mthd='UNKNOWN';\n";
    mDynamicModelFile << "  end;\n";
    mDynamicModelFile << "  disp (['-----------------------------------------------------']) ;\n";
    mDynamicModelFile << "  disp (['MODEL SIMULATION: (method=' mthd ')']) ;\n";
    mDynamicModelFile << "  fprintf('\\n') ;\n";
    mDynamicModelFile << "  periods=options_.periods;\n";
    mDynamicModelFile << "  maxit_=options_.maxit_;\n";
    mDynamicModelFile << "  solve_tolf=options_.solve_tolf;\n";
    mDynamicModelFile << "  y=oo_.endo_simul';\n";
    mDynamicModelFile << "  x=oo_.exo_simul;\n";

    prev_Simulation_Type=-1;
    mDynamicModelFile << "  params=M_.params;\n";
    mDynamicModelFile << "  oo_.deterministic_simulation.status = 0;\n";
    for (i = 0;i < block_triangular.ModelBlock->Size;i++)
      {
        k = block_triangular.ModelBlock->Block_List[i].Simulation_Type;
        if (BlockTriangular::BlockSim(prev_Simulation_Type)==BlockTriangular::BlockSim(k) &&
            (k==EVALUATE_FORWARD || k==EVALUATE_BACKWARD /*|| k==EVALUATE_FORWARD_R || k==EVALUATE_BACKWARD_R*/))
          skip_head=true;
        else
          skip_head=false;
        if ((k == EVALUATE_FORWARD /*|| k == EVALUATE_FORWARD_R*/) && (block_triangular.ModelBlock->Block_List[i].Size))
          {
            if (!skip_head)
              {
                if (open_par)
                  {
                    mDynamicModelFile << "  end\n";
                  }
                mDynamicModelFile << "  oo_.deterministic_simulation.status = 1;\n";
                mDynamicModelFile << "  oo_.deterministic_simulation.error = 0;\n";
                mDynamicModelFile << "  oo_.deterministic_simulation.iterations = 0;\n";
                mDynamicModelFile << "  if(isfield(oo_.deterministic_simulation,'block'))\n";
                mDynamicModelFile << "    blck_num = length(oo_.deterministic_simulation.block)+1;\n";
                mDynamicModelFile << "  else\n";
                mDynamicModelFile << "    blck_num = 1;\n";
                mDynamicModelFile << "  end;\n";
                mDynamicModelFile << "  oo_.deterministic_simulation.block(blck_num).status = 1;\n";
                mDynamicModelFile << "  oo_.deterministic_simulation.block(blck_num).error = 0;\n";
                mDynamicModelFile << "  oo_.deterministic_simulation.block(blck_num).iterations = 0;\n";
                mDynamicModelFile << "  g1=[];g2=[];g3=[];\n";
                //mDynamicModelFile << "  for it_ = y_kmin+1:(periods+y_kmin)\n";
                mDynamicModelFile << "  y=" << dynamic_basename << "_" << i + 1 << "(y, x, params, 0, y_kmin, periods);\n";
                mDynamicModelFile << "  tmp = y(:,M_.block_structure.block(" << i + 1 << ").variable);\n";
                mDynamicModelFile << "  if(isnan(tmp) | isinf(tmp))\n";
                mDynamicModelFile << "    disp(['Inf or Nan value during the evaluation of block " << i <<"']);\n";
                mDynamicModelFile << "    return;\n";
                mDynamicModelFile << "  end;\n";
              }
            //open_par=true;
          }
        else if ((k == EVALUATE_BACKWARD /*|| k == EVALUATE_BACKWARD_R*/) && (block_triangular.ModelBlock->Block_List[i].Size))
          {
            if (!skip_head)
              {
                if (open_par)
                  {
                    mDynamicModelFile << "  end\n";
                  }
                mDynamicModelFile << "  oo_.deterministic_simulation.status = 1;\n";
                mDynamicModelFile << "  oo_.deterministic_simulation.error = 0;\n";
                mDynamicModelFile << "  oo_.deterministic_simulation.iterations = 0;\n";
                mDynamicModelFile << "  if(isfield(oo_.deterministic_simulation,'block'))\n";
                mDynamicModelFile << "    blck_num = length(oo_.deterministic_simulation.block)+1;\n";
                mDynamicModelFile << "  else\n";
                mDynamicModelFile << "    blck_num = 1;\n";
                mDynamicModelFile << "  end;\n";
                mDynamicModelFile << "  oo_.deterministic_simulation.block(blck_num).status = 1;\n";
                mDynamicModelFile << "  oo_.deterministic_simulation.block(blck_num).error = 0;\n";
                mDynamicModelFile << "  oo_.deterministic_simulation.block(blck_num).iterations = 0;\n";
                mDynamicModelFile << "  g1=[];g2=[];g3=[];\n";
                mDynamicModelFile << "  " << dynamic_basename << "_" << i + 1 << "(y, x, params, 0, y_kmin, periods);\n";
                mDynamicModelFile << "  tmp = y(:,M_.block_structure.block(" << i + 1 << ").variable);\n";
                mDynamicModelFile << "  if(isnan(tmp) | isinf(tmp))\n";
                mDynamicModelFile << "    disp(['Inf or Nan value during the evaluation of block " << i <<"']);\n";
                mDynamicModelFile << "    return;\n";
                mDynamicModelFile << "  end;\n";
              }
          }
        else if ((k == SOLVE_FORWARD_COMPLETE || k == SOLVE_FORWARD_SIMPLE) && (block_triangular.ModelBlock->Block_List[i].Size))
          {
            if (open_par)
              mDynamicModelFile << "  end\n";
            open_par=false;
            mDynamicModelFile << "  g1=0;\n";
            mDynamicModelFile << "  r=0;\n";
            tmp.str("");
            for (int ik=block_triangular.ModelBlock->Block_List[i].Nb_Recursives ;ik<block_triangular.ModelBlock->Block_List[i].Size;ik++)
              {
                tmp << " " << block_triangular.ModelBlock->Block_List[i].Variable[ik]+1;
              }
            mDynamicModelFile << "  y_index = [" << tmp.str() << "];\n";
            int nze, m;
            for (nze=0,m=0;m<=block_triangular.ModelBlock->Block_List[i].Max_Lead+block_triangular.ModelBlock->Block_List[i].Max_Lag;m++)
              nze+=block_triangular.ModelBlock->Block_List[i].IM_lead_lag[m].size;
            mDynamicModelFile << "  if(isfield(oo_.deterministic_simulation,'block'))\n";
            mDynamicModelFile << "    blck_num = length(oo_.deterministic_simulation.block)+1;\n";
            mDynamicModelFile << "  else\n";
            mDynamicModelFile << "    blck_num = 1;\n";
            mDynamicModelFile << "  end;\n";
            mDynamicModelFile << "  y = solve_one_boundary('"  << dynamic_basename << "_" <<  i + 1 << "'" <<
            ", y, x, params, y_index, " << nze <<
            ", options_.periods, " << block_triangular.ModelBlock->Block_List[i].is_linear <<
              ", blck_num, y_kmin, options_.maxit_, options_.solve_tolf, options_.slowc, " << cutoff << ", options_.stack_solve_algo, 1, 1, 0);\n";
            mDynamicModelFile << "  tmp = y(:,M_.block_structure.block(" << i + 1 << ").variable);\n";
            mDynamicModelFile << "  if(isnan(tmp) | isinf(tmp))\n";
            mDynamicModelFile << "    disp(['Inf or Nan value during the resolution of block " << i <<"']);\n";
            mDynamicModelFile << "    return;\n";
            mDynamicModelFile << "  end;\n";
          }
        else if ((k == SOLVE_BACKWARD_COMPLETE || k == SOLVE_BACKWARD_SIMPLE) && (block_triangular.ModelBlock->Block_List[i].Size))
          {
            if (open_par)
              mDynamicModelFile << "  end\n";
            open_par=false;
            mDynamicModelFile << "  g1=0;\n";
            mDynamicModelFile << "  r=0;\n";
            tmp.str("");
            for (int ik=block_triangular.ModelBlock->Block_List[i].Nb_Recursives ;ik<block_triangular.ModelBlock->Block_List[i].Size;ik++)
              {
                 tmp << " " << block_triangular.ModelBlock->Block_List[i].Variable[ik]+1;
              }
            mDynamicModelFile << "  y_index = [" << tmp.str() << "];\n";
            int nze, m;
            for (nze=0,m=0;m<=block_triangular.ModelBlock->Block_List[i].Max_Lead+block_triangular.ModelBlock->Block_List[i].Max_Lag;m++)
              nze+=block_triangular.ModelBlock->Block_List[i].IM_lead_lag[m].size;
            mDynamicModelFile << "  if(isfield(oo_.deterministic_simulation,'block'))\n";
            mDynamicModelFile << "    blck_num = length(oo_.deterministic_simulation.block)+1;\n";
            mDynamicModelFile << "  else\n";
            mDynamicModelFile << "    blck_num = 1;\n";
            mDynamicModelFile << "  end;\n";
            mDynamicModelFile << "  y = solve_one_boundary('"  << dynamic_basename << "_" <<  i + 1 << "'" <<
            ", y, x, params, y_index, " << nze <<
            ", options_.periods, " << block_triangular.ModelBlock->Block_List[i].is_linear <<
            ", blck_num, y_kmin, options_.maxit_, options_.solve_tolf, options_.slowc, " << cutoff << ", options_.stack_solve_algo, 1, 1, 0);\n";
            mDynamicModelFile << "  tmp = y(:,M_.block_structure.block(" << i + 1 << ").variable);\n";
            mDynamicModelFile << "  if(isnan(tmp) | isinf(tmp))\n";
            mDynamicModelFile << "    disp(['Inf or Nan value during the resolution of block " << i <<"']);\n";
            mDynamicModelFile << "    return;\n";
            mDynamicModelFile << "  end;\n";
          }
        else if ((k == SOLVE_TWO_BOUNDARIES_COMPLETE || k == SOLVE_TWO_BOUNDARIES_SIMPLE) && (block_triangular.ModelBlock->Block_List[i].Size))
          {
            if (open_par)
              mDynamicModelFile << "  end\n";
            open_par=false;
            Nb_SGE++;
            int nze, m;
            for (nze=0,m=0;m<=block_triangular.ModelBlock->Block_List[i].Max_Lead+block_triangular.ModelBlock->Block_List[i].Max_Lag;m++)
              nze+=block_triangular.ModelBlock->Block_List[i].IM_lead_lag[m].size;
            mDynamicModelFile << "  y_index=[";
            for (int ik=block_triangular.ModelBlock->Block_List[i].Nb_Recursives ;ik<block_triangular.ModelBlock->Block_List[i].Size;ik++)
              {
                mDynamicModelFile << " " << block_triangular.ModelBlock->Block_List[i].Variable[ik]+1;
              }
            mDynamicModelFile << "  ];\n";
            mDynamicModelFile << "  if(isfield(oo_.deterministic_simulation,'block'))\n";
            mDynamicModelFile << "    blck_num = length(oo_.deterministic_simulation.block)+1;\n";
            mDynamicModelFile << "  else\n";
            mDynamicModelFile << "    blck_num = 1;\n";
            mDynamicModelFile << "  end;\n";
            mDynamicModelFile << "  y = solve_two_boundaries('" << dynamic_basename << "_" <<  i + 1 << "'" <<
            ", y, x, params, y_index, " << nze <<
            ", options_.periods, " << block_triangular.ModelBlock->Block_List[i].Max_Lag <<
            ", " << block_triangular.ModelBlock->Block_List[i].Max_Lead <<
            ", " << block_triangular.ModelBlock->Block_List[i].is_linear <<
            ", blck_num, y_kmin, options_.maxit_, options_.solve_tolf, options_.slowc, " << cutoff << ", options_.stack_solve_algo);\n";
            mDynamicModelFile << "  tmp = y(:,M_.block_structure.block(" << i + 1 << ").variable);\n";
            mDynamicModelFile << "  if(isnan(tmp) | isinf(tmp))\n";
            mDynamicModelFile << "    disp(['Inf or Nan value during the resolution of block " << i <<"']);\n";
            mDynamicModelFile << "    return;\n";
            mDynamicModelFile << "  end;\n";
          }
        prev_Simulation_Type=k;
      }
    if (open_par)
      mDynamicModelFile << "  end;\n";
    open_par=false;
    mDynamicModelFile << "  oo_.endo_simul = y';\n";
    mDynamicModelFile << "return;\n";

    mDynamicModelFile.close();

    writeModelEquationsOrdered_M( block_triangular.ModelBlock, dynamic_basename);

    chdir("..");
  }

void
DynamicModel::writeDynamicModel(ostream &DynamicOutput, bool use_dll) const
  {
    ostringstream model_output;    // Used for storing model equations
    ostringstream jacobian_output; // Used for storing jacobian equations
    ostringstream hessian_output;  // Used for storing Hessian equations
    ostringstream third_derivatives_output;

    ExprNodeOutputType output_type = (use_dll ? oCDynamicModel : oMatlabDynamicModel);

    writeModelLocalVariables(model_output, output_type);

    writeTemporaryTerms(temporary_terms, model_output, output_type);

    writeModelEquations(model_output, output_type);

    int nrows = equations.size();
    int hessianColsNbr = dynJacobianColsNbr * dynJacobianColsNbr;

    // Writing Jacobian
    for (first_derivatives_type::const_iterator it = first_derivatives.begin();
         it != first_derivatives.end(); it++)
      {
        int eq = it->first.first;
        int var = it->first.second;
        NodeID d1 = it->second;

        jacobian_output << "g1";
        jacobianHelper(jacobian_output, eq, getDynJacobianCol(var), output_type);
        jacobian_output << "=";
        d1->writeOutput(jacobian_output, output_type, temporary_terms);
        jacobian_output << ";" << endl;
      }

    // Writing Hessian
    int k = 0; // Keep the line of a 2nd derivative in v2
    for (second_derivatives_type::const_iterator it = second_derivatives.begin();
         it != second_derivatives.end(); it++)
      {
        int eq = it->first.first;
        int var1 = it->first.second.first;
        int var2 = it->first.second.second;
        NodeID d2 = it->second;

        int id1 = getDynJacobianCol(var1);
        int id2 = getDynJacobianCol(var2);

        int col_nb = id1 * dynJacobianColsNbr + id2;
        int col_nb_sym = id2 * dynJacobianColsNbr + id1;

        sparseHelper(2, hessian_output, k, 0, output_type);
        hessian_output << "=" << eq + 1 << ";" << endl;

        sparseHelper(2, hessian_output, k, 1, output_type);
        hessian_output << "=" << col_nb + 1 << ";" << endl;

        sparseHelper(2, hessian_output, k, 2, output_type);
        hessian_output << "=";
        d2->writeOutput(hessian_output, output_type, temporary_terms);
        hessian_output << ";" << endl;

        k++;

        // Treating symetric elements
        if (id1 != id2)
          {
            sparseHelper(2, hessian_output, k, 0, output_type);
            hessian_output << "=" << eq + 1 << ";" << endl;

            sparseHelper(2, hessian_output, k, 1, output_type);
            hessian_output << "=" << col_nb_sym + 1 << ";" << endl;

            sparseHelper(2, hessian_output, k, 2, output_type);
            hessian_output << "=";
            sparseHelper(2, hessian_output, k-1, 2, output_type);
            hessian_output << ";" << endl;

            k++;
          }
      }

    // Writing third derivatives
    k = 0; // Keep the line of a 3rd derivative in v3
    for (third_derivatives_type::const_iterator it = third_derivatives.begin();
         it != third_derivatives.end(); it++)
      {
        int eq = it->first.first;
        int var1 = it->first.second.first;
        int var2 = it->first.second.second.first;
        int var3 = it->first.second.second.second;
        NodeID d3 = it->second;

        int id1 = getDynJacobianCol(var1);
        int id2 = getDynJacobianCol(var2);
        int id3 = getDynJacobianCol(var3);

        // Reference column number for the g3 matrix
        int ref_col = id1 * hessianColsNbr + id2 * dynJacobianColsNbr + id3;

        sparseHelper(3, third_derivatives_output, k, 0, output_type);
        third_derivatives_output << "=" << eq + 1 << ";" << endl;

        sparseHelper(3, third_derivatives_output, k, 1, output_type);
        third_derivatives_output << "=" << ref_col + 1 << ";" << endl;

        sparseHelper(3, third_derivatives_output, k, 2, output_type);
        third_derivatives_output << "=";
        d3->writeOutput(third_derivatives_output, output_type, temporary_terms);
        third_derivatives_output << ";" << endl;

        k++;

        // Compute the column numbers for the 5 other permutations of (id1,id2,id3) and store them in a set (to avoid duplicates if two indexes are equal)
        set<int> cols;
        cols.insert(id1 * hessianColsNbr + id3 * dynJacobianColsNbr + id2);
        cols.insert(id2 * hessianColsNbr + id1 * dynJacobianColsNbr + id3);
        cols.insert(id2 * hessianColsNbr + id3 * dynJacobianColsNbr + id1);
        cols.insert(id3 * hessianColsNbr + id1 * dynJacobianColsNbr + id2);
        cols.insert(id3 * hessianColsNbr + id2 * dynJacobianColsNbr + id1);

        int k2 = 0; // Keeps the offset of the permutation relative to k
        for (set<int>::iterator it2 = cols.begin(); it2 != cols.end(); it2++)
          if (*it2 != ref_col)
            {
              sparseHelper(3, third_derivatives_output, k+k2, 0, output_type);
              third_derivatives_output << "=" << eq + 1 << ";" << endl;

              sparseHelper(3, third_derivatives_output, k+k2, 1, output_type);
              third_derivatives_output << "=" << *it2 + 1 << ";" << endl;

              sparseHelper(3, third_derivatives_output, k+k2, 2, output_type);
              third_derivatives_output << "=";
              sparseHelper(3, third_derivatives_output, k, 2, output_type);
              third_derivatives_output << ";" << endl;

              k2++;
            }
        k += k2;
      }

    if (!use_dll)
      {
        DynamicOutput << "%" << endl
        << "% Model equations" << endl
        << "%" << endl
        << endl
        << "residual = zeros(" << nrows << ", 1);" << endl
        << model_output.str()
        // Writing initialization instruction for matrix g1
        << "if nargout >= 2," << endl
        << "  g1 = zeros(" << nrows << ", " << dynJacobianColsNbr << ");" << endl
        << endl
        << "%" << endl
        << "% Jacobian matrix" << endl
        << "%" << endl
        << endl
        << jacobian_output.str()
        << "end" << endl;

        // Initialize g2 matrix
        DynamicOutput << "if nargout >= 3," << endl
                      << "%" << endl
                      << "% Hessian matrix" << endl
                      << "%" << endl
                      << endl;
        if (second_derivatives.size())
          DynamicOutput << "  v2 = zeros(" << NNZDerivatives[1] << ",3);" << endl
                        << hessian_output.str()
                        << "  g2 = sparse(v2(:,1),v2(:,2),v2(:,3)," << nrows << "," << hessianColsNbr << ");" << endl;
        else // Either hessian is all zero, or we didn't compute it
          DynamicOutput << "  g2 = sparse([],[],[]," << nrows << "," << hessianColsNbr << ");" << endl;
        DynamicOutput << "end;" << endl;

        // Initialize g3 matrix
        DynamicOutput << "if nargout >= 4," << endl
                      << "%" << endl
                      << "% Third order derivatives" << endl
                      << "%" << endl
                      << endl;
        int ncols = hessianColsNbr * dynJacobianColsNbr;
        if (third_derivatives.size())
          DynamicOutput << "  v3 = zeros(" << NNZDerivatives[2] << ",3);" << endl
                        << third_derivatives_output.str()
                        << "  g3 = sparse(v3(:,1),v3(:,2),v3(:,3)," << nrows << "," << ncols << ");" << endl;
        else // Either 3rd derivatives is all zero, or we didn't compute it
          DynamicOutput << "  g3 = sparse([],[],[]," << nrows << "," << ncols << ");" << endl;

        DynamicOutput << "end;" << endl;
      }
    else
      {
        DynamicOutput << "void Dynamic(double *y, double *x, int nb_row_x, double *params, int it_, double *residual, double *g1, double *v2, double *v3)" << endl
        << "{" << endl
        << "  double lhs, rhs;" << endl
        << endl
        << "  /* Residual equations */" << endl
        << model_output.str()
        << "  /* Jacobian  */" << endl
        << "  if (g1 == NULL)" << endl
        << "    return;" << endl
        << "  else" << endl
        << "    {" << endl
        << jacobian_output.str()
        << "    }" << endl;

        if (second_derivatives.size())
            DynamicOutput << "  /* Hessian for endogenous and exogenous variables */" << endl
            << "  if (v2 == NULL)" << endl
            << "    return;" << endl
            << "  else" << endl
            << "    {" << endl
            << hessian_output.str()
            << "    }" << endl;

        if (third_derivatives.size())
          DynamicOutput << "  /* Third derivatives for endogenous and exogenous variables */" << endl
                        << "  if (v3 == NULL)" << endl
                        << "    return;" << endl
                        << "  else" << endl
                        << "    {" << endl
                        << third_derivatives_output.str()
                        << "    }" << endl;

        DynamicOutput << "}" << endl << endl;
      }
  }

void
DynamicModel::writeOutput(ostream &output, const string &basename, bool block, bool byte_code, bool use_dll) const
  {
    /* Writing initialisation for M_.lead_lag_incidence matrix
       M_.lead_lag_incidence is a matrix with as many columns as there are
       endogenous variables and as many rows as there are periods in the
       models (nbr of rows = M_.max_lag+M_.max_lead+1)

       The matrix elements are equal to zero if a variable isn't present in the
       model at a given period.
    */

    output << "M_.lead_lag_incidence = [";
    // Loop on endogenous variables
    for (int endoID = 0; endoID < symbol_table.endo_nbr(); endoID++)
      {
        output << endl;
        // Loop on periods
        for (int lag = -max_endo_lag; lag <= max_endo_lead; lag++)
          {
            // Print variableID if exists with current period, otherwise print 0
            try
              {
                int varID = getDerivID(symbol_table.getID(eEndogenous, endoID), lag);
                output << " " << getDynJacobianCol(varID) + 1;
              }
            catch (UnknownDerivIDException &e)
              {
                output << " 0";
              }
          }
        output << ";";
      }
    output << "]';" << endl;

    // Write equation tags
    output << "M_.equations_tags = {" << endl;
    for (unsigned int i = 0; i < equation_tags.size(); i++)
      output << "  " << equation_tags[i].first + 1 << " , '"
             << equation_tags[i].second.first << "' , '"
             << equation_tags[i].second.second << "' ;" << endl;
    output << "};" << endl;

    //In case of sparse model, writes the block structure of the model
    if (block)
      {
        //int prev_Simulation_Type=-1;
        //bool skip_the_head;
        int k=0;
        int count_lead_lag_incidence = 0;
        int max_lead, max_lag, max_lag_endo, max_lead_endo, max_lag_exo, max_lead_exo;
        for (int j = 0;j < block_triangular.ModelBlock->Size;j++)
          {
            //For a block composed of a single equation determines wether we have to evaluate or to solve the equation
            //skip_the_head=false;
            k++;
            count_lead_lag_incidence = 0;
            int Block_size=block_triangular.ModelBlock->Block_List[j].Size;
            max_lag =block_triangular.ModelBlock->Block_List[j].Max_Lag ;
            max_lead=block_triangular.ModelBlock->Block_List[j].Max_Lead;
            max_lag_endo =block_triangular.ModelBlock->Block_List[j].Max_Lag_Endo ;
            max_lead_endo=block_triangular.ModelBlock->Block_List[j].Max_Lead_Endo;
            max_lag_exo =block_triangular.ModelBlock->Block_List[j].Max_Lag_Exo ;
            max_lead_exo=block_triangular.ModelBlock->Block_List[j].Max_Lead_Exo;
            bool evaluate=false;
            vector<int> exogenous;
            vector<int>::iterator it_exogenous;
            exogenous.clear();
            ostringstream tmp_s, tmp_s_eq;
            tmp_s.str("");
            tmp_s_eq.str("");
            for (int i=0;i<block_triangular.ModelBlock->Block_List[j].Size;i++)
              {
                tmp_s << " " << block_triangular.ModelBlock->Block_List[j].Variable[i]+1;
                tmp_s_eq << " " << block_triangular.ModelBlock->Block_List[j].Equation[i]+1;
              }
            for (int i=0;i<block_triangular.ModelBlock->Block_List[j].nb_exo;i++)
              {
                int ii=block_triangular.ModelBlock->Block_List[j].Exogenous[i];
                for (it_exogenous=exogenous.begin();it_exogenous!=exogenous.end() && *it_exogenous!=ii;it_exogenous++) /*cout << "*it_exogenous=" << *it_exogenous << "\n"*/;
                if (it_exogenous==exogenous.end() || exogenous.begin()==exogenous.end())
                  exogenous.push_back(ii);
              }
            output << "M_.block_structure.block(" << k << ").num = " << j+1 << ";\n";
            output << "M_.block_structure.block(" << k << ").Simulation_Type = " << block_triangular.ModelBlock->Block_List[j].Simulation_Type << ";\n";
            output << "M_.block_structure.block(" << k << ").maximum_lag = " << max_lag << ";\n";
            output << "M_.block_structure.block(" << k << ").maximum_lead = " << max_lead << ";\n";
            output << "M_.block_structure.block(" << k << ").maximum_endo_lag = " << max_lag_endo << ";\n";
            output << "M_.block_structure.block(" << k << ").maximum_endo_lead = " << max_lead_endo << ";\n";
            output << "M_.block_structure.block(" << k << ").maximum_exo_lag = " << max_lag_exo << ";\n";
            output << "M_.block_structure.block(" << k << ").maximum_exo_lead = " << max_lead_exo << ";\n";
            output << "M_.block_structure.block(" << k << ").endo_nbr = " << Block_size << ";\n";
            output << "M_.block_structure.block(" << k << ").equation = [" << tmp_s_eq.str() << "];\n";
            output << "M_.block_structure.block(" << k << ").variable = [" << tmp_s.str() << "];\n";
            output << "M_.block_structure.block(" << k << ").exogenous = [";
            int i=0;
            for (it_exogenous=exogenous.begin();it_exogenous!=exogenous.end();it_exogenous++)
              if (*it_exogenous>=0)
                {
                  output << " " << *it_exogenous+1;
                  i++;
                }
            output << "];\n";
            output << "M_.block_structure.block(" << k << ").exo_nbr = " << i << ";\n";

            output << "M_.block_structure.block(" << k << ").exo_det_nbr = " << block_triangular.ModelBlock->Block_List[j].nb_exo_det << ";\n";

            tmp_s.str("");

            bool done_IM=false;
            if (!evaluate)
              {
                output << "M_.block_structure.block(" << k << ").lead_lag_incidence = [];\n";
                for (int l=-max_lag_endo;l<max_lead_endo+1;l++)
                  {
                    bool *tmp_IM;
                    tmp_IM=block_triangular.incidencematrix.Get_IM(l, eEndogenous);
                    if (tmp_IM)
                      {
                        for (int l_var=0;l_var<block_triangular.ModelBlock->Block_List[j].Size;l_var++)
                          {
                            for (int l_equ=0;l_equ<block_triangular.ModelBlock->Block_List[j].Size;l_equ++)
                              if (tmp_IM[block_triangular.ModelBlock->Block_List[j].Equation[l_equ]*symbol_table.endo_nbr()+block_triangular.ModelBlock->Block_List[j].Variable[l_var]])
                                {
                                  count_lead_lag_incidence++;
                                  if (tmp_s.str().length())
                                    tmp_s << " ";
                                  tmp_s << count_lead_lag_incidence;
                                  done_IM=true;
                                  break;
                                }
                            if (!done_IM)
                              tmp_s << " 0";
                            done_IM=false;
                          }
                        output << "M_.block_structure.block(" << k << ").lead_lag_incidence = [ M_.block_structure.block(" << k << ").lead_lag_incidence; " << tmp_s.str() << "];\n";
                        tmp_s.str("");
                      }
                  }
              }
            else
              {
                bool done_some_where;
                output << "M_.block_structure.block(" << k << ").lead_lag_incidence = [\n";
                for (int l=-max_lag_endo;l<max_lead_endo+1;l++)
                  {
                    bool not_increm=true;
                    bool *tmp_IM;
                    tmp_IM=block_triangular.incidencematrix.Get_IM(l, eEndogenous);
                    int ii=j;
                    if (tmp_IM)
                      {
                        done_some_where = false;
                        while (ii-j<Block_size)
                          {
                            for (int l_var=0;l_var<block_triangular.ModelBlock->Block_List[ii].Size;l_var++)
                              {
                                for (int l_equ=0;l_equ<block_triangular.ModelBlock->Block_List[ii].Size;l_equ++)
                                  if (tmp_IM[block_triangular.ModelBlock->Block_List[ii].Equation[l_equ]*symbol_table.endo_nbr()+block_triangular.ModelBlock->Block_List[ii].Variable[l_var]])
                                    {
                                      //if(not_increm && l==-max_lag)
                                      count_lead_lag_incidence++;
                                      not_increm=false;
                                      if (tmp_s.str().length())
                                        tmp_s << " ";
                                      //tmp_s << count_lead_lag_incidence+(l+max_lag)*Block_size;
                                      tmp_s << count_lead_lag_incidence;
                                      done_IM=true;
                                      break;
                                    }
                                if (!done_IM)
                                  tmp_s << " 0";
                                else
                                  done_some_where = true;
                                done_IM=false;
                              }
                            ii++;
                          }
                        output << tmp_s.str() << "\n";
                        tmp_s.str("");
                      }
                  }
                output << "];\n";
              }

          }
        for (int j=-block_triangular.incidencematrix.Model_Max_Lag_Endo;j<=block_triangular.incidencematrix.Model_Max_Lead_Endo;j++)
          {
            bool* IM = block_triangular.incidencematrix.Get_IM(j, eEndogenous);
            if (IM)
              {
                bool new_entry=true;
                output << "M_.block_structure.incidence(" << block_triangular.incidencematrix.Model_Max_Lag_Endo+j+1 << ").lead_lag = " << j << ";\n";
                output << "M_.block_structure.incidence(" << block_triangular.incidencematrix.Model_Max_Lag_Endo+j+1 << ").sparse_IM = [";
                for (int i=0;i<symbol_table.endo_nbr()*symbol_table.endo_nbr();i++)
                  {
                    if (IM[i])
                      {
                        if (!new_entry)
                          output << " ; ";
                        else
                          output << " ";
                        output << i/symbol_table.endo_nbr()+1 << " " << i % symbol_table.endo_nbr()+1;
                        new_entry=false;
                      }
                  }
                output << "];\n";
              }
          }
      }
    // Writing initialization for some other variables
    output << "M_.exo_names_orig_ord = [1:" << symbol_table.exo_nbr() << "];" << endl
    << "M_.maximum_lag = " << max_lag << ";" << endl
    << "M_.maximum_lead = " << max_lead << ";" << endl;
    if (symbol_table.endo_nbr())
      {
        output << "M_.maximum_endo_lag = " << max_endo_lag << ";" << endl
        << "M_.maximum_endo_lead = " << max_endo_lead << ";" << endl
        << "oo_.steady_state = zeros(" << symbol_table.endo_nbr() << ", 1);" << endl;
      }
    if (symbol_table.exo_nbr())
      {
        output << "M_.maximum_exo_lag = " << max_exo_lag << ";" << endl
        << "M_.maximum_exo_lead = " << max_exo_lead << ";" << endl
        << "oo_.exo_steady_state = zeros(" << symbol_table.exo_nbr() << ", 1);" << endl;
      }
    if (symbol_table.exo_det_nbr())
      {
        output << "M_.maximum_exo_det_lag = " << max_exo_det_lag << ";" << endl
        << "M_.maximum_exo_det_lead = " << max_exo_det_lead << ";" << endl
        << "oo_.exo_det_steady_state = zeros(" << symbol_table.exo_det_nbr() << ", 1);" << endl;
      }
    if (symbol_table.param_nbr())
      output << "M_.params = repmat(NaN," << symbol_table.param_nbr() << ", 1);" << endl;

    // Write number of non-zero derivatives
    output << "M_.NNZDerivatives = zeros(3, 1);" << endl
           << "M_.NNZDerivatives(1) = " << NNZDerivatives[0] << ";" << endl
           << "M_.NNZDerivatives(2) = " << NNZDerivatives[1] << ";" << endl
           << "M_.NNZDerivatives(3) = " << NNZDerivatives[2] << ";" << endl;
  }

void
DynamicModel::evaluateJacobian(const eval_context_type &eval_context, jacob_map *j_m, bool dynamic)
{
  int i=0;
  int j=0;
  bool *IM=NULL;
  int a_variable_lag=-9999;
  for (first_derivatives_type::iterator it = first_derivatives.begin();
       it != first_derivatives.end(); it++)
    {
      //cout << "it->first.second=" << it->first.second << " variable_table.getSymbolID(it->first.second)=" << variable_table.getSymbolID(it->first.second) << " Type=" << variable_table.getType(it->first.second) << " eEndogenous=" << eEndogenous << " eExogenous=" << eExogenous << " variable_table.getLag(it->first.second)=" << variable_table.getLag(it->first.second) << "\n";
      if (getTypeByDerivID(it->first.second) == eEndogenous)
        {
          NodeID Id = it->second;
          double val = 0;
          try
            {
              val = Id->eval(eval_context);
            }
          catch (ExprNode::EvalException &e)
            {
              cout << "evaluation of Jacobian failed for equation " << it->first.first+1 << " and variable " << symbol_table.getName(getSymbIDByDerivID(it->first.second)) << "(" << getLagByDerivID(it->first.second) << ") [" << getSymbIDByDerivID(it->first.second) << "] !" << endl;
              Id->writeOutput(cout, oMatlabDynamicModelSparse, temporary_terms);
              cout << "\n";
              cerr << "DynamicModel::evaluateJacobian: evaluation of Jacobian failed for equation " << it->first.first+1 << " and variable " << symbol_table.getName(getSymbIDByDerivID(it->first.second)) << "(" << getLagByDerivID(it->first.second) << ")!" << endl;
            }
          int eq=it->first.first;
          int var = symbol_table.getTypeSpecificID(getSymbIDByDerivID(it->first.second));///symbol_table.getID(eEndogenous,it->first.second);//variable_table.getSymbolID(it->first.second);
          int k1 = getLagByDerivID(it->first.second);
          if (a_variable_lag!=k1)
            {
              IM=block_triangular.incidencematrix.Get_IM(k1, eEndogenous);
              a_variable_lag=k1;
            }
          if (k1==0 or !dynamic)
            {
              j++;
              (*j_m)[make_pair(eq,var)]+=val;
            }
          if (IM[eq*symbol_table.endo_nbr()+var] && (fabs(val) < cutoff))
            {
              if (block_triangular.bt_verbose)
                cout << "the coefficient related to variable " << var << " with lag " << k1 << " in equation " << eq << " is equal to " << val << " and is set to 0 in the incidence matrix (size=" << symbol_table.endo_nbr() << ")\n";
              block_triangular.incidencematrix.unfill_IM(eq, var, k1, eEndogenous);
              i++;
            }
        }
    }
  //Get ride of the elements of the incidence matrix equal to Zero
  IM=block_triangular.incidencematrix.Get_IM(0, eEndogenous);
  for (int i=0;i<symbol_table.endo_nbr();i++)
    for (int j=0;j<symbol_table.endo_nbr();j++)
      if (IM[i*symbol_table.endo_nbr()+j])
        if (first_derivatives.find(make_pair(i,getDerivID(symbol_table.getID(eEndogenous, j), 0)))==first_derivatives.end())
          block_triangular.incidencematrix.unfill_IM(i, j, 0, eEndogenous);
  if (i>0)
    {
      cout << i << " elements among " << first_derivatives.size() << " in the incidence matrices are below the cutoff (" << cutoff << ") and are discarded\n";
      cout << "the contemporaneous incidence matrix has " << j << " elements\n";
    }
}

void
DynamicModel::BlockLinear(Model_Block *ModelBlock)
{
  int i,j,l,m,ll;
  for (j = 0;j < ModelBlock->Size;j++)
    {
      if (ModelBlock->Block_List[j].Simulation_Type==SOLVE_BACKWARD_COMPLETE ||
          ModelBlock->Block_List[j].Simulation_Type==SOLVE_FORWARD_COMPLETE)
        {
          ll=ModelBlock->Block_List[j].Max_Lag;
          for (i=0;i<ModelBlock->Block_List[j].IM_lead_lag[ll].size;i++)
            {
              int eq=ModelBlock->Block_List[j].IM_lead_lag[ll].Equ_Index[i];
              int var=ModelBlock->Block_List[j].IM_lead_lag[ll].Var_Index[i];
              //first_derivatives_type::const_iterator it=first_derivatives.find(make_pair(eq,variable_table.getID(var,0)));
              first_derivatives_type::const_iterator it=first_derivatives.find(make_pair(eq,getDerivID(symbol_table.getID(eEndogenous, var),0)));
              if (it!= first_derivatives.end())
                {
                  NodeID Id = it->second;
                  set<pair<int, int> > endogenous;
                  Id->collectEndogenous(endogenous);
                  if (endogenous.size() > 0)
                    {
                      for (l=0;l<ModelBlock->Block_List[j].Size;l++)
                        {
                          if (endogenous.find(make_pair(ModelBlock->Block_List[j].Variable[l], 0)) != endogenous.end())
                            {
                              ModelBlock->Block_List[j].is_linear=false;
                              goto follow;
                            }
                        }
                    }
                }
            }
        }
      else if (ModelBlock->Block_List[j].Simulation_Type==SOLVE_TWO_BOUNDARIES_COMPLETE || ModelBlock->Block_List[j].Simulation_Type==SOLVE_TWO_BOUNDARIES_SIMPLE)
        {
          for (m=0;m<=ModelBlock->Block_List[j].Max_Lead+ModelBlock->Block_List[j].Max_Lag;m++)
            {
              int k1=m-ModelBlock->Block_List[j].Max_Lag;
              for (i=0;i<ModelBlock->Block_List[j].IM_lead_lag[m].size;i++)
                {
                  int eq=ModelBlock->Block_List[j].IM_lead_lag[m].Equ_Index[i];
                  int var=ModelBlock->Block_List[j].IM_lead_lag[m].Var_Index[i];
                  //first_derivatives_type::const_iterator it=first_derivatives.find(make_pair(eq,variable_table.getID(var,k1)));
                  first_derivatives_type::const_iterator it=first_derivatives.find(make_pair(eq,getDerivID(symbol_table.getID(eEndogenous, var),k1)));
                  NodeID Id = it->second;
                  if (it!= first_derivatives.end())
                    {
                      set<pair<int, int> > endogenous;
                      Id->collectEndogenous(endogenous);
                      if (endogenous.size() > 0)
                        {
                          for (l=0;l<ModelBlock->Block_List[j].Size;l++)
                            {
                              if (endogenous.find(make_pair(ModelBlock->Block_List[j].Variable[l], k1)) != endogenous.end())
                                {
                                  ModelBlock->Block_List[j].is_linear=false;
                                  goto follow;
                                }
                            }
                        }
                    }
                }
            }
        }
follow:
      i=0;
    }
}


map<pair<int, pair<int, int > >, NodeID>
DynamicModel::collect_first_order_derivatives_endogenous()
{
  map<pair<int, pair<int, int > >, NodeID> endo_derivatives;
  for (first_derivatives_type::iterator it2 = first_derivatives.begin();
       it2 != first_derivatives.end(); it2++)
    {
      if (getTypeByDerivID(it2->first.second)==eEndogenous)
        {
          int eq = it2->first.first;
          int var=symbol_table.getTypeSpecificID(getSymbIDByDerivID(it2->first.second));
          int lag=getLagByDerivID(it2->first.second);
          //if (lag==0)
          endo_derivatives[make_pair(eq, make_pair(var, lag))] = it2->second;
        }
    }
  return  endo_derivatives;
}



void
DynamicModel::computingPass(bool jacobianExo, bool hessian, bool thirdDerivatives, bool paramsDerivatives,
                            const eval_context_type &eval_context, bool no_tmp_terms, bool block, bool use_dll)
{
  assert(jacobianExo || !(hessian || thirdDerivatives || paramsDerivatives));

  // Prepare for derivation
  computeDerivIDs();

  // Computes dynamic jacobian columns, must be done after computeDerivIDs()
  computeDynJacobianCols(jacobianExo);

  // Compute derivatives w.r. to all endogenous, and possibly exogenous and exogenous deterministic
  set<int> vars;
  for (deriv_id_table_t::const_iterator it = deriv_id_table.begin();
       it != deriv_id_table.end(); it++)
    {
      SymbolType type = symbol_table.getType(it->first.first);
      if (type == eEndogenous || (jacobianExo && (type == eExogenous || type == eExogenousDet)))
        vars.insert(it->second);
    }

  // Launch computations
  cout << "Computing dynamic model derivatives:" << endl
  << " - order 1" << endl;
  computeJacobian(vars);

  if (hessian)
    {
      cout << " - order 2" << endl;
      computeHessian(vars);
    }

  if (paramsDerivatives)
    {
      cout << " - order 2 (derivatives of Jacobian w.r. to parameters)" << endl;
      computeParamsDerivatives();

      if (!no_tmp_terms)
        computeParamsDerivativesTemporaryTerms();
    }

  if (thirdDerivatives)
    {
      cout << " - order 3" << endl;
      computeThirdDerivatives(vars);
    }

  if (block)
    {
      BuildIncidenceMatrix();

      jacob_map j_m;
      evaluateJacobian(eval_context, &j_m, true);


      if (block_triangular.bt_verbose)
        {
          cout << "The gross incidence matrix \n";
          block_triangular.incidencematrix.Print_IM(eEndogenous);
        }
      t_etype equation_simulation_type;
      map<pair<int, pair<int, int> >, NodeID> first_order_endo_derivatives = collect_first_order_derivatives_endogenous();

      block_triangular.Normalize_and_BlockDecompose_Static_0_Model(j_m, equations, equation_simulation_type, first_order_endo_derivatives, mfs, cutoff);

      BlockLinear(block_triangular.ModelBlock);

      computeChainRuleJacobian(block_triangular.ModelBlock);

      if (!no_tmp_terms)
        computeTemporaryTermsOrdered(block_triangular.ModelBlock);

    }
  else
    if (!no_tmp_terms)
      computeTemporaryTerms(!use_dll);
}

void
DynamicModel::writeDynamicFile(const string &basename, bool block, bool bytecode, bool use_dll) const
  {
    int r;
    if(block && bytecode)
      {
        writeModelEquationsCodeOrdered(basename + "_dynamic", block_triangular.ModelBlock, basename, map_idx);
        block_triangular.Free_Block(block_triangular.ModelBlock);
        block_triangular.incidencematrix.Free_IM();
        //block_triangular.Free_IM_X(block_triangular.First_IM_X);
      }
		else if(block && !bytecode)
		  {
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
        writeSparseDynamicMFile(basename + "_dynamic", basename);
        block_triangular.Free_Block(block_triangular.ModelBlock);
        block_triangular.incidencematrix.Free_IM();
        //block_triangular.Free_IM_X(block_triangular.First_IM_X);
		  }
		else if (use_dll)
     	writeDynamicCFile(basename + "_dynamic");
		else
		  writeDynamicMFile(basename + "_dynamic");
  }

void
DynamicModel::toStatic(StaticModel &static_model) const
  {
    // Convert model local variables (need to be done first)
    for (map<int, NodeID>::const_iterator it = local_variables_table.begin();
         it != local_variables_table.end(); it++)
      static_model.AddLocalVariable(symbol_table.getName(it->first), it->second->toStatic(static_model));

    // Convert equations
    for (vector<BinaryOpNode *>::const_iterator it = equations.begin();
         it != equations.end(); it++)
      static_model.addEquation((*it)->toStatic(static_model));

    // Convert auxiliary equations
    for (deque<BinaryOpNode *>::const_iterator it = aux_equations.begin();
         it != aux_equations.end(); it++)
      static_model.addAuxEquation((*it)->toStatic(static_model));
  }

void
DynamicModel::toStaticDll(StaticDllModel &static_model) const
  {
    // Convert model local variables (need to be done first)
    for (map<int, NodeID>::const_iterator it = local_variables_table.begin();
         it != local_variables_table.end(); it++)
      static_model.AddLocalVariable(symbol_table.getName(it->first), it->second->toStatic(static_model));

    // Convert equations
    for (vector<BinaryOpNode *>::const_iterator it = equations.begin();
         it != equations.end(); it++)
      static_model.addEquation((*it)->toStatic(static_model));
  }

void
DynamicModel::computeDerivIDs()
{
  set<pair<int, int> > dynvars;

  for(int i = 0; i < (int) equations.size(); i++)
    equations[i]->collectVariables(eEndogenous, dynvars);

  dynJacobianColsNbr = dynvars.size();

  for(int i = 0; i < (int) equations.size(); i++)
    {
      equations[i]->collectVariables(eExogenous, dynvars);
      equations[i]->collectVariables(eExogenousDet, dynvars);
      equations[i]->collectVariables(eParameter, dynvars);
    }

  for(set<pair<int, int> >::const_iterator it = dynvars.begin();
      it != dynvars.end(); it++)
    {
      int lag = it->second;
      SymbolType type = symbol_table.getType(it->first);

      // Setting maximum and minimum lags
      if (max_lead < lag)
        max_lead = lag;
      else if (-max_lag > lag)
        max_lag = -lag;

      switch (type)
        {
        case eEndogenous:
          if (max_endo_lead < lag)
            max_endo_lead = lag;
          else if (-max_endo_lag > lag)
            max_endo_lag = -lag;
          break;
        case eExogenous:
          if (max_exo_lead < lag)
            max_exo_lead = lag;
          else if (-max_exo_lag > lag)
            max_exo_lag = -lag;
          break;
        case eExogenousDet:
          if (max_exo_det_lead < lag)
            max_exo_det_lead = lag;
          else if (-max_exo_det_lag > lag)
            max_exo_det_lag = -lag;
          break;
        default:
          break;
        }

      // Create a new deriv_id
      int deriv_id = deriv_id_table.size();

      deriv_id_table[*it] = deriv_id;
      inv_deriv_id_table.push_back(*it);
    }
}

SymbolType
DynamicModel::getTypeByDerivID(int deriv_id) const throw (UnknownDerivIDException)
{
  return symbol_table.getType(getSymbIDByDerivID(deriv_id));
}

int
DynamicModel::getLagByDerivID(int deriv_id) const throw (UnknownDerivIDException)
{
  if (deriv_id < 0 || deriv_id >= (int) inv_deriv_id_table.size())
    throw UnknownDerivIDException();

  return inv_deriv_id_table[deriv_id].second;
}

int
DynamicModel::getSymbIDByDerivID(int deriv_id) const throw (UnknownDerivIDException)
{
  if (deriv_id < 0 || deriv_id >= (int) inv_deriv_id_table.size())
    throw UnknownDerivIDException();

  return inv_deriv_id_table[deriv_id].first;
}

int
DynamicModel::getDerivID(int symb_id, int lag) const throw (UnknownDerivIDException)
{
  deriv_id_table_t::const_iterator it = deriv_id_table.find(make_pair(symb_id, lag));
  if (it == deriv_id_table.end())
    throw UnknownDerivIDException();
  else
    return it->second;
}

void
DynamicModel::computeDynJacobianCols(bool jacobianExo)
{
  /* Sort the dynamic endogenous variables by lexicographic order over (lag, type_specific_symbol_id)
     and fill the dynamic columns for exogenous and exogenous deterministic */
  map<pair<int, int>, int> ordered_dyn_endo;

  for (deriv_id_table_t::const_iterator it = deriv_id_table.begin();
       it != deriv_id_table.end(); it++)
    {
      const int &symb_id = it->first.first;
      const int &lag = it->first.second;
      const int &deriv_id = it->second;
      SymbolType type = symbol_table.getType(symb_id);
      int tsid = symbol_table.getTypeSpecificID(symb_id);

      switch (type)
        {
        case eEndogenous:
          ordered_dyn_endo[make_pair(lag, tsid)] = deriv_id;
          break;
        case eExogenous:
          // At this point, dynJacobianColsNbr contains the number of dynamic endogenous
          if (jacobianExo)
            dyn_jacobian_cols_table[deriv_id] = dynJacobianColsNbr + tsid;
          break;
        case eExogenousDet:
          // At this point, dynJacobianColsNbr contains the number of dynamic endogenous
          if (jacobianExo)
            dyn_jacobian_cols_table[deriv_id] = dynJacobianColsNbr + symbol_table.exo_nbr() + tsid;
          break;
        case eParameter:
          // We don't assign a dynamic jacobian column to parameters
          break;
        default:
          // Shut up GCC
          cerr << "DynamicModel::computeDynJacobianCols: impossible case" << endl;
          exit(EXIT_FAILURE);
        }
    }

  // Fill in dynamic jacobian columns for endogenous
  int sorted_id = 0;
  for (map<pair<int, int>, int>::const_iterator it = ordered_dyn_endo.begin();
       it != ordered_dyn_endo.end(); it++)
    dyn_jacobian_cols_table[it->second] = sorted_id++;

  // Set final value for dynJacobianColsNbr
  if (jacobianExo)
    dynJacobianColsNbr += symbol_table.exo_nbr() + symbol_table.exo_det_nbr();
}

int
DynamicModel::getDynJacobianCol(int deriv_id) const throw (UnknownDerivIDException)
{
  map<int, int>::const_iterator it = dyn_jacobian_cols_table.find(deriv_id);
  if (it == dyn_jacobian_cols_table.end())
    throw UnknownDerivIDException();
  else
    return it->second;
}


void
DynamicModel::computeChainRuleJacobian(Model_Block *ModelBlock)
{
  map<int, NodeID> recursive_variables;
  first_chain_rule_derivatives.clear();
  for(int blck = 0; blck<ModelBlock->Size; blck++)
    {
      recursive_variables.clear();
      if (ModelBlock->Block_List[blck].Simulation_Type==SOLVE_TWO_BOUNDARIES_COMPLETE or ModelBlock->Block_List[blck].Simulation_Type==SOLVE_TWO_BOUNDARIES_SIMPLE)
        {
          ModelBlock->Block_List[blck].Chain_Rule_Derivatives->clear();
          for(int i = 0; i < ModelBlock->Block_List[blck].Nb_Recursives; i++)
            {
              if (ModelBlock->Block_List[blck].Equation_Type[i] == E_EVALUATE_S)
                recursive_variables[getDerivID(symbol_table.getID(eEndogenous, ModelBlock->Block_List[blck].Variable[i]), 0)] = ModelBlock->Block_List[blck].Equation_Normalized[i];
              else
                recursive_variables[getDerivID(symbol_table.getID(eEndogenous, ModelBlock->Block_List[blck].Variable[i]), 0)] = equations[ModelBlock->Block_List[blck].Equation[i]];
            }
          map<pair<pair<int, pair<int, int> >, pair<int, int> >, int> Derivatives = block_triangular.get_Derivatives(ModelBlock, blck);

          map<pair<pair<int, pair<int, int> >, pair<int, int> >, int>::const_iterator it = Derivatives.begin();
          //#pragma omp parallel for shared(it, blck)
          for(int i=0; i<(int)Derivatives.size(); i++)
            {
              int Deriv_type = it->second;
              pair<pair<int, pair<int, int> >, pair<int, int> > it_l(it->first);
              it++;
              int lag = it_l.first.first;
              int eq = it_l.first.second.first;
              int var = it_l.first.second.second;
              int eqr = it_l.second.first;
              int varr = it_l.second.second;
              if(Deriv_type == 0)
                {
                  first_chain_rule_derivatives[make_pair(eqr, make_pair(varr, lag))] = first_derivatives[make_pair(eqr, getDerivID(symbol_table.getID(eEndogenous, varr), lag))];
                }
              else if (Deriv_type == 1)
                {
                  first_chain_rule_derivatives[make_pair(eqr, make_pair(varr, lag))] = ModelBlock->Block_List[blck].Equation_Normalized[eq]->getChainRuleDerivative(getDerivID(symbol_table.getID(eEndogenous, varr), lag), recursive_variables);
                }
              else if (Deriv_type == 2)
                {
                  if(ModelBlock->Block_List[blck].Equation_Type[eq] == E_EVALUATE_S && eq<ModelBlock->Block_List[blck].Nb_Recursives)
                    first_chain_rule_derivatives[make_pair(eqr, make_pair(varr, lag))] = ModelBlock->Block_List[blck].Equation_Normalized[eq]->getChainRuleDerivative(getDerivID(symbol_table.getID(eEndogenous, varr), lag), recursive_variables);
                  else
                    first_chain_rule_derivatives[make_pair(eqr, make_pair(varr, lag))] = equations[eqr]->getChainRuleDerivative(getDerivID(symbol_table.getID(eEndogenous, varr), lag), recursive_variables);
                }
              ModelBlock->Block_List[blck].Chain_Rule_Derivatives->push_back(make_pair( make_pair(lag, make_pair(eq, var)), make_pair(eqr, varr)));
            }
        }
      else if(   ModelBlock->Block_List[blck].Simulation_Type==SOLVE_BACKWARD_SIMPLE or ModelBlock->Block_List[blck].Simulation_Type==SOLVE_FORWARD_SIMPLE
              or ModelBlock->Block_List[blck].Simulation_Type==SOLVE_BACKWARD_COMPLETE or ModelBlock->Block_List[blck].Simulation_Type==SOLVE_FORWARD_COMPLETE)
        {
          ModelBlock->Block_List[blck].Chain_Rule_Derivatives->clear();
          for(int i = 0; i < ModelBlock->Block_List[blck].Nb_Recursives; i++)
            {
              if (ModelBlock->Block_List[blck].Equation_Type[i] == E_EVALUATE_S)
                recursive_variables[getDerivID(symbol_table.getID(eEndogenous, ModelBlock->Block_List[blck].Variable[i]), 0)] = ModelBlock->Block_List[blck].Equation_Normalized[i];
              else
                recursive_variables[getDerivID(symbol_table.getID(eEndogenous, ModelBlock->Block_List[blck].Variable[i]), 0)] = equations[ModelBlock->Block_List[blck].Equation[i]];
            }
          for(int eq = ModelBlock->Block_List[blck].Nb_Recursives; eq < ModelBlock->Block_List[blck].Size; eq++)
            {
              int eqr = ModelBlock->Block_List[blck].Equation[eq];
              for(int var = ModelBlock->Block_List[blck].Nb_Recursives; var < ModelBlock->Block_List[blck].Size; var++)
                {
                  int varr = ModelBlock->Block_List[blck].Variable[var];
                  NodeID d1 = equations[eqr]->getChainRuleDerivative(getDerivID(symbol_table.getID(eEndogenous, varr), 0), recursive_variables);
                  if (d1 == Zero)
                    continue;
                  first_chain_rule_derivatives[make_pair(eqr, make_pair(varr, 0))] = d1;
                  ModelBlock->Block_List[blck].Chain_Rule_Derivatives->push_back(make_pair( make_pair(0, make_pair(eq, var)), make_pair(eqr, varr)));
                }
            }
        }
    }
}



void
DynamicModel::computeParamsDerivatives()
{
  for (deriv_id_table_t::const_iterator it = deriv_id_table.begin();
       it != deriv_id_table.end(); it++)
    {
      if (symbol_table.getType(it->first.first) != eParameter)
        continue;

      int param = it->second;

      for (int eq = 0; eq < (int) equations.size(); eq++)
        {
          NodeID d1 = equations[eq]->getDerivative(param);
          if (d1 == Zero)
            continue;
          residuals_params_derivatives[make_pair(eq, param)] = d1;
        }

      for (first_derivatives_type::const_iterator it2 = first_derivatives.begin();
           it2 != first_derivatives.end(); it2++)
        {
          int eq = it2->first.first;
          int var = it2->first.second;
          NodeID d1 = it2->second;

          NodeID d2 = d1->getDerivative(param);
          if (d2 == Zero)
            continue;
          jacobian_params_derivatives[make_pair(eq, make_pair(var, param))] = d2;
        }
    }
}

void
DynamicModel::computeParamsDerivativesTemporaryTerms()
{
  map<NodeID, int> reference_count;
  params_derivs_temporary_terms.clear();

  for (first_derivatives_type::iterator it = residuals_params_derivatives.begin();
       it != residuals_params_derivatives.end(); it++)
    it->second->computeTemporaryTerms(reference_count, params_derivs_temporary_terms, true);

  for (second_derivatives_type::iterator it = jacobian_params_derivatives.begin();
       it != jacobian_params_derivatives.end(); it++)
    it->second->computeTemporaryTerms(reference_count, params_derivs_temporary_terms, true);
}

void
DynamicModel::writeParamsDerivativesFile(const string &basename) const
  {
    if (!residuals_params_derivatives.size()
        && !jacobian_params_derivatives.size())
      return;

    string filename = basename + "_params_derivs.m";

    ofstream paramsDerivsFile;
    paramsDerivsFile.open(filename.c_str(), ios::out | ios::binary);
    if (!paramsDerivsFile.is_open())
      {
        cerr << "ERROR: Can't open file " << filename << " for writing" << endl;
        exit(EXIT_FAILURE);
      }
    paramsDerivsFile << "function [rp, gp] = " << basename << "_params_derivs(y, x, params, it_)" << endl
    << "%" << endl
    << "% Warning : this file is generated automatically by Dynare" << endl
    << "%           from model file (.mod)" << endl << endl;


    writeTemporaryTerms(params_derivs_temporary_terms, paramsDerivsFile, oMatlabDynamicModel);

    // Write parameter derivative
    paramsDerivsFile << "rp = zeros(" << equation_number() << ", "
                     << symbol_table.param_nbr() << ");" << endl;

    for (first_derivatives_type::const_iterator it = residuals_params_derivatives.begin();
         it != residuals_params_derivatives.end(); it++)
      {
        int eq = it->first.first;
        int param = it->first.second;
        NodeID d1 = it->second;

        int param_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param)) + 1;

        paramsDerivsFile << "rp(" << eq+1 << ", " << param_col << ") = ";
        d1->writeOutput(paramsDerivsFile, oMatlabDynamicModel, params_derivs_temporary_terms);
        paramsDerivsFile << ";" << endl;
      }

    // Write jacobian derivatives
    paramsDerivsFile << "gp = zeros(" << equation_number() << ", " << dynJacobianColsNbr << ", "
    << symbol_table.param_nbr() << ");" << endl;

    for (second_derivatives_type::const_iterator it = jacobian_params_derivatives.begin();
         it != jacobian_params_derivatives.end(); it++)
      {
        int eq = it->first.first;
        int var = it->first.second.first;
        int param = it->first.second.second;
        NodeID d2 = it->second;

        int var_col = getDynJacobianCol(var) + 1;
        int param_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param)) + 1;

        paramsDerivsFile << "gp(" << eq+1 << ", " << var_col << ", " << param_col << ") = ";
        d2->writeOutput(paramsDerivsFile, oMatlabDynamicModel, params_derivs_temporary_terms);
        paramsDerivsFile << ";" << endl;
      }

    paramsDerivsFile.close();
  }



void
DynamicModel::writeChainRuleDerivative(ostream &output, int eqr, int varr, int lag,
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
DynamicModel::writeLatexFile(const string &basename) const
  {
    writeLatexModelFile(basename + "_dynamic.tex", oLatexDynamicModel);
  }

void
DynamicModel::jacobianHelper(ostream &output, int eq_nb, int col_nb, ExprNodeOutputType output_type) const
{
  output << LEFT_ARRAY_SUBSCRIPT(output_type);
  if (IS_MATLAB(output_type))
    output << eq_nb + 1 << "," << col_nb + 1;
  else
    output << eq_nb + col_nb * equations.size();
  output << RIGHT_ARRAY_SUBSCRIPT(output_type);
}

void
DynamicModel::sparseHelper(int order, ostream &output, int row_nb, int col_nb, ExprNodeOutputType output_type) const
{
  output << "v" << order << LEFT_ARRAY_SUBSCRIPT(output_type);
  if (IS_MATLAB(output_type))
    output << row_nb + 1 << "," << col_nb + 1;
  else
    output << row_nb + col_nb * NNZDerivatives[order-1];
  output << RIGHT_ARRAY_SUBSCRIPT(output_type);
}

void
DynamicModel::substituteEndoLeadGreaterThanTwo()
{
  substituteLeadLagInternal(avEndoLead);
}

void
DynamicModel::substituteEndoLagGreaterThanTwo()
{
  substituteLeadLagInternal(avEndoLag);
}

void
DynamicModel::substituteExoLead()
{
  substituteLeadLagInternal(avExoLead);
}

void
DynamicModel::substituteExoLag()
{
  substituteLeadLagInternal(avExoLag);
}

void
DynamicModel::substituteLeadLagInternal(aux_var_t type)
{
  ExprNode::subst_table_t subst_table;
  vector<BinaryOpNode *> neweqs;

  // Substitute in model local variables
  for(map<int, NodeID>::iterator it = local_variables_table.begin();
      it != local_variables_table.end(); it++)
    {
      NodeID subst;
      switch(type)
        {
        case avEndoLead:
          subst = it->second->substituteEndoLeadGreaterThanTwo(subst_table, neweqs);
          break;
        case avEndoLag:
          subst = it->second->substituteEndoLagGreaterThanTwo(subst_table, neweqs);
          break;
        case avExoLead:
          subst = it->second->substituteExoLead(subst_table, neweqs);
          break;
        case avExoLag:
          subst = it->second->substituteExoLag(subst_table, neweqs);
          break;
        default:
          cerr << "DynamicModel::substituteLeadLagInternal: impossible case" << endl;
          exit(EXIT_FAILURE);
        }
      it->second = subst;
    }

  // Substitute in equations
  for(int i = 0; i < (int) equations.size(); i++)
    {
      NodeID subst;
      switch(type)
        {
        case avEndoLead:
          subst = equations[i]->substituteEndoLeadGreaterThanTwo(subst_table, neweqs);
          break;
        case avEndoLag:
          subst = equations[i]->substituteEndoLagGreaterThanTwo(subst_table, neweqs);
          break;
        case avExoLead:
          subst = equations[i]->substituteExoLead(subst_table, neweqs);
          break;
        case avExoLag:
          subst = equations[i]->substituteExoLag(subst_table, neweqs);
          break;
        default:
          cerr << "DynamicModel::substituteLeadLagInternal: impossible case" << endl;
          exit(EXIT_FAILURE);
        }
      BinaryOpNode *substeq = dynamic_cast<BinaryOpNode *>(subst);
      assert(substeq != NULL);
      equations[i] = substeq;
    }

  // Add new equations
  for(int i = 0; i < (int) neweqs.size(); i++)
    addEquation(neweqs[i]);

  // Add the new set of equations at the *beginning* of aux_equations
  copy(neweqs.rbegin(), neweqs.rend(), front_inserter(aux_equations));

  if (neweqs.size() > 0)
    {
      cout << "Substitution of ";
      switch(type)
        {
        case avEndoLead:
          cout << "endo leads >= 2";
          break;
        case avEndoLag:
          cout << "endo lags >= 2";
          break;
        case avExoLead:
          cout << "exo leads";
          break;
        case avExoLag:
          cout << "exo lags";
          break;
        case avExpectation:
          cout << "expectation";
          break;
        }
      cout << ": added " << neweqs.size() << " auxiliary variables and equations." << endl;
    }
}

void
DynamicModel::substituteExpectation(bool partial_information_model)
{
  ExprNode::subst_table_t subst_table;
  vector<BinaryOpNode *> neweqs;

  // Substitute in model binary op node map
  for(unary_op_node_map_type::reverse_iterator it = unary_op_node_map.rbegin();
      it != unary_op_node_map.rend(); it++)
    it->second->substituteExpectation(subst_table, neweqs, partial_information_model);

  // Substitute in equations
  for(int i = 0; i < (int) equations.size(); i++)
    {
      BinaryOpNode *substeq = dynamic_cast<BinaryOpNode *>(equations[i]->substituteExpectation(subst_table, neweqs, partial_information_model));
      assert(substeq != NULL);
      equations[i] = substeq;
    }

  // Add new equations
  for(int i = 0; i < (int) neweqs.size(); i++)
    addEquation(neweqs[i]);

  // Add the new set of equations at the *beginning* of aux_equations
  copy(neweqs.rbegin(), neweqs.rend(), front_inserter(aux_equations));

  if (neweqs.size() > 0)
    if (partial_information_model)
      cout << "Substitution of Expectation operator: added " << subst_table.size() << " auxiliary variables and " << neweqs.size() << " auxiliary equations." << endl;
    else
      cout << "Substitution of Expectation operator: added " << neweqs.size() << " auxiliary variables and equations." << endl;
}

void
DynamicModel::transformPredeterminedVariables()
{
  for(int i = 0; i < (int) equations.size(); i++)
    {
      BinaryOpNode *substeq = dynamic_cast<BinaryOpNode *>(equations[i]->decreaseLeadsLagsPredeterminedVariables());
      assert(substeq != NULL);
      equations[i] = substeq;
    }
}

void
DynamicModel::fillEvalContext(eval_context_type &eval_context) const
{
  // First, auxiliary variables
  for(deque<BinaryOpNode *>::const_iterator it = aux_equations.begin();
      it != aux_equations.end(); it++)
    {
      assert((*it)->get_op_code() == oEqual);
      VariableNode *auxvar = dynamic_cast<VariableNode *>((*it)->get_arg1());
      assert(auxvar != NULL);
      try
        {
          double val = (*it)->get_arg2()->eval(eval_context);
          eval_context[auxvar->get_symb_id()] = val;
        }
      catch(ExprNode::EvalException &e)
        {
          // Do nothing
        }
    }

  // Second, model local variables
  for(map<int, NodeID>::const_iterator it = local_variables_table.begin();
      it != local_variables_table.end(); it++)
    {
      try
        {
          const NodeID expression = it->second;
          double val = expression->eval(eval_context);
          eval_context[it->first] = val;
        }
      catch(ExprNode::EvalException &e)
        {
          // Do nothing
        }
    }
}
