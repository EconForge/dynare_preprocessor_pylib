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
#include "StaticDllModel.hh"

// For mkdir() and chdir()
#ifdef _WIN32
# include <direct.h>
#else
# include <unistd.h>
# include <sys/stat.h>
# include <sys/types.h>
#endif

StaticDllModel::StaticDllModel(SymbolTable &symbol_table_arg,
                           NumericalConstants &num_constants_arg) :
    ModelTree(symbol_table_arg, num_constants_arg),
    cutoff(1e-15),
    mfs(0),
    block_triangular(symbol_table_arg, num_constants_arg)
{
}

void
StaticDllModel::compileDerivative(ofstream &code_file, int eq, int symb_id, int lag, map_idx_type &map_idx) const
{
  //first_derivatives_type::const_iterator it = first_derivatives.find(make_pair(eq, getDerivID(symb_id, lag)));
  //first_derivatives_type::const_iterator it = first_derivatives.find(make_pair(eq, getDerivID(symbol_table.getID(eEndogenous, symb_id), lag)));
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
StaticDllModel::compileChainRuleDerivative(ofstream &code_file, int eqr, int varr, int lag, map_idx_type &map_idx) const
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
StaticDllModel::BuildIncidenceMatrix()
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
          block_triangular.incidencematrix.fill_IM(eq, it_endogenous->first, 0, eEndogenous);
        }
      exogenous.clear();
      Id = eq_node->get_arg1();
      Id->collectExogenous(exogenous);
      Id = eq_node->get_arg2();
      Id->collectExogenous(exogenous);
      for (set<pair<int, int> >::iterator it_exogenous=exogenous.begin();it_exogenous!=exogenous.end();it_exogenous++)
        {
          block_triangular.incidencematrix.fill_IM(eq, it_exogenous->first, 0, eExogenous);
        }
    }
}

void
StaticDllModel::computeTemporaryTermsOrdered(Model_Block *ModelBlock)
{
  map<NodeID, pair<int, int> > first_occurence;
  map<NodeID, int> reference_count;
  int i, j, eqr, varr, lag;
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
    }
  // Add a mapping form node ID to temporary terms order
  j=0;
  for (temporary_terms_type::const_iterator it = temporary_terms.begin();
       it != temporary_terms.end(); it++)
    map_idx[(*it)->idx]=j++;
}

void
StaticDllModel::writeModelEquationsOrdered_M( Model_Block *ModelBlock, const string &static_basename) const
  {
    int i,j,k,m;
    string tmp_s, sps;
    ostringstream tmp_output, tmp1_output, global_output;
    NodeID lhs=NULL, rhs=NULL;
    BinaryOpNode *eq_node;
    map<NodeID, int> reference_count;
    ofstream  output;
    int nze, nze_exo, nze_other_endo;
    vector<int> feedback_variables;
    //For each block
    for (j = 0;j < ModelBlock->Size;j++)
      {
        //recursive_variables.clear();
        feedback_variables.clear();
        //For a block composed of a single equation determines wether we have to evaluate or to solve the equation
        nze = nze_exo = nze_other_endo = 0;
        for (m=0;m<=ModelBlock->Block_List[j].Max_Lead+ModelBlock->Block_List[j].Max_Lag;m++)
          nze+=ModelBlock->Block_List[j].IM_lead_lag[m].size;
        tmp1_output.str("");
        tmp1_output << static_basename << "_" << j+1 << ".m";
        output.open(tmp1_output.str().c_str(), ios::out | ios::binary);
        output << "%\n";
        output << "% " << tmp1_output.str() << " : Computes static model for Dynare\n";
        output << "%\n";
        output << "% Warning : this file is generated automatically by Dynare\n";
        output << "%           from model file (.mod)\n\n";
        output << "%/\n";
        if (ModelBlock->Block_List[j].Simulation_Type==EVALUATE_BACKWARD
            ||ModelBlock->Block_List[j].Simulation_Type==EVALUATE_FORWARD)
          output << "function y = " << static_basename << "_" << j+1 << "(y, x, params)\n";
        else if (ModelBlock->Block_List[j].Simulation_Type==SOLVE_FORWARD_COMPLETE
                 ||   ModelBlock->Block_List[j].Simulation_Type==SOLVE_BACKWARD_COMPLETE
                 ||   ModelBlock->Block_List[j].Simulation_Type==SOLVE_BACKWARD_SIMPLE
                 ||   ModelBlock->Block_List[j].Simulation_Type==SOLVE_FORWARD_SIMPLE)
          output << "function [residual, y, g1] = " << static_basename << "_" << j+1 << "(y, x, params)\n";
        output << "  % ////////////////////////////////////////////////////////////////////////" << endl
        << "  % //" << string("                     Block ").substr(int(log10(j + 1))) << j + 1 << " " << BlockTriangular::BlockType0(ModelBlock->Block_List[j].Type)
        << "          //" << endl
        << "  % //                     Simulation type "
        << BlockTriangular::BlockSim(ModelBlock->Block_List[j].Simulation_Type) << "  //" << endl
        << "  % ////////////////////////////////////////////////////////////////////////" << endl;
        //The Temporary terms
        if (ModelBlock->Block_List[j].Simulation_Type!=EVALUATE_BACKWARD
            && ModelBlock->Block_List[j].Simulation_Type!=EVALUATE_FORWARD)
          output << "    g1 = spalloc(" << ModelBlock->Block_List[j].Size-ModelBlock->Block_List[j].Nb_Recursives
                 << ", " << ModelBlock->Block_List[j].Size-ModelBlock->Block_List[j].Nb_Recursives << ", " << nze << ");\n";

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
                (*it)->writeOutput(output, oMatlabStaticModelSparse, temporary_terms);
                output << " = ";
                (*it)->writeOutput(output, oMatlabStaticModelSparse, tt2);
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
              lhs->writeOutput(tmp_output, oMatlabStaticModelSparse, temporary_terms);
            else*/
						lhs->writeOutput(tmp_output, oMatlabStaticModelSparse, temporary_terms);
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
  									rhs->writeOutput(output, oMatlabStaticModelSparse, temporary_terms);
                  }
                else if (ModelBlock->Block_List[j].Equation_Type[i] == E_EVALUATE_S)
                  {
                    output << "%" << tmp_output.str();
                    output << " = ";
                    if (ModelBlock->Block_List[j].Equation_Normalized[i])
                      {
                        rhs->writeOutput(output, oMatlabStaticModelSparse, temporary_terms);
                        output << "\n    ";
                        tmp_output.str("");
                        eq_node = (BinaryOpNode *)ModelBlock->Block_List[j].Equation_Normalized[i];
                        lhs = eq_node->get_arg1();
                        rhs = eq_node->get_arg2();
                        lhs->writeOutput(output, oMatlabStaticModelSparse, temporary_terms);
                        output << " = ";
											  rhs->writeOutput(output, oMatlabStaticModelSparse, temporary_terms);
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
                  goto evaluation;
                feedback_variables.push_back(ModelBlock->Block_List[j].Variable[i]);
                output << "  % equation " << ModelBlock->Block_List[j].Equation[i]+1 << " variable : " << sModel
                << " (" << ModelBlock->Block_List[j].Variable[i]+1 << ") " << block_triangular.c_Equation_Type(ModelBlock->Block_List[j].Equation_Type[i]) << endl;
                output << "  " << "residual(" << i+1-ModelBlock->Block_List[j].Nb_Recursives << ") = (";
                goto end;
              default:
end:
                output << tmp_output.str();
                output << ") - (";
                rhs->writeOutput(output, oMatlabStaticModelSparse, temporary_terms);
                output << ");\n";
              }
          }
        // The Jacobian if we have to solve the block
        output << "  " << sps << "% Jacobian  " << endl;
        switch (ModelBlock->Block_List[j].Simulation_Type)
          {
          case EVALUATE_BACKWARD:
          case EVALUATE_FORWARD:
            break;
          case SOLVE_BACKWARD_SIMPLE:
          case SOLVE_FORWARD_SIMPLE:
          case SOLVE_BACKWARD_COMPLETE:
          case SOLVE_FORWARD_COMPLETE:
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
                writeChainRuleDerivative(output, eqr, varr, k, oMatlabStaticModelSparse, temporary_terms);
                output << "; % variable=" << symbol_table.getName(symbol_table.getID(eEndogenous, varr))
                       << " " << varr+1 << ", equation=" << eqr+1 << endl;
              }
            break;
          default:
            break;
          }
        output.close();
      }
  }

void
StaticDllModel::writeModelEquationsCodeOrdered(const string file_name, const Model_Block *ModelBlock, const string bin_basename, map_idx_type map_idx) const
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
    FDIMT_ fdimt(temporary_terms.size());
    fdimt.write(code_file);

    for (j = 0; j < ModelBlock->Size ;j++)
      {
        feedback_variables.clear();
        if (j>0)
          {
            FENDBLOCK_ fendblock;
            fendblock.write(code_file);
          }
        int count_u;

        int u_count_int=0;
        if (ModelBlock->Block_List[j].Simulation_Type==SOLVE_BACKWARD_COMPLETE || ModelBlock->Block_List[j].Simulation_Type==SOLVE_FORWARD_COMPLETE)
          {
            Write_Inf_To_Bin_File(file_name, bin_basename, j, u_count_int,file_open);
            file_open=true;
          }

        FBEGINBLOCK_ fbeginblock(ModelBlock->Block_List[j].Size - ModelBlock->Block_List[j].Nb_Recursives,
                                 ModelBlock->Block_List[j].Simulation_Type,
                                 ModelBlock->Block_List[j].Variable,
                                 ModelBlock->Block_List[j].Equation,
                                 ModelBlock->Block_List[j].Own_Derivative,
                                 ModelBlock->Block_List[j].is_linear,
                                 symbol_table.endo_nbr(),
                                 0,
                                 0,
                                 u_count_int
                                 );
        fbeginblock.write(code_file);
            // The equations
            //cout << block_triangular.BlockSim(ModelBlock->Block_List[j].Simulation_Type) << "  j=" << j << endl;
            for (i = 0;i < ModelBlock->Block_List[j].Size;i++)
              {
                //The Temporary terms
                temporary_terms_type tt2;
                tt2.clear();
                for (temporary_terms_type::const_iterator it = ModelBlock->Block_List[j].Temporary_Terms_in_Equation[i]->begin();
                     it != ModelBlock->Block_List[j].Temporary_Terms_in_Equation[i]->end(); it++)
                  {
                    (*it)->compile(code_file, false, tt2, map_idx, false, false);
                    FSTPST_ fstpst((int)(map_idx.find((*it)->idx))->second);
                    fstpst.write(code_file);
                    // Insert current node into tt2
                    tt2.insert(*it);
                  }
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
                        rhs->compile(code_file, false, temporary_terms, map_idx, false, false);
                        lhs->compile(code_file, true, temporary_terms, map_idx, false, false);
                      }
                    else if (ModelBlock->Block_List[j].Equation_Type[i] == E_EVALUATE_S)
                      {
                        eq_node = (BinaryOpNode*)ModelBlock->Block_List[j].Equation_Normalized[i];
                        lhs = eq_node->get_arg1();
                        rhs = eq_node->get_arg2();
                        rhs->compile(code_file, false, temporary_terms, map_idx, false, false);
                        lhs->compile(code_file, true, temporary_terms, map_idx, false, false);
                      }
                    break;
                  case SOLVE_BACKWARD_COMPLETE:
                  case SOLVE_FORWARD_COMPLETE:
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
                    lhs->compile(code_file, false, temporary_terms, map_idx, false, false);
                    rhs->compile(code_file, false, temporary_terms, map_idx, false, false);
                    FBINARY_ fbinary(oMinus);
                    fbinary.write(code_file);

                    FSTPR_ fstpr(i - ModelBlock->Block_List[j].Nb_Recursives);
                    fstpr.write(code_file);
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
                        FSTPG_ fstpg;
                        fstpg.write(code_file);
                      }
                    break;

                  case SOLVE_BACKWARD_COMPLETE:
                  case SOLVE_FORWARD_COMPLETE:
										count_u = feedback_variables.size();
                    for(i=0; i<(int)ModelBlock->Block_List[j].Chain_Rule_Derivatives->size();i++)
											{
                        pair< pair<int, pair<int, int> >, pair<int, int> > it = ModelBlock->Block_List[j].Chain_Rule_Derivatives->at(i);
                        k=it.first.first;
                        int eq=it.first.second.first;
                        int var=it.first.second.second;
                        int eqr=it.second.first;
                        int varr=it.second.second;
                        int v=ModelBlock->Block_List[j].Equation[eq];
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
                            FSTPSU_ fstpsu(count_u);
                            fstpsu.write(code_file);
                            count_u++;
												  }
											}
                    for (i = 0;i < ModelBlock->Block_List[j].Size;i++)
                      {
                      	if(i>=ModelBlock->Block_List[j].Nb_Recursives)
                      	  {
                            FLDR_ fldr(i-ModelBlock->Block_List[j].Nb_Recursives);
                            fldr.write(code_file);

                            FLDZ_ fldz;
                            fldz.write(code_file);

                            v=ModelBlock->Block_List[j].Equation[i];
                            for (Uf[v].Ufl=Uf[v].Ufl_First; Uf[v].Ufl; Uf[v].Ufl=Uf[v].Ufl->pNext)
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
                            Uf[v].Ufl=Uf[v].Ufl_First;
                            while (Uf[v].Ufl)
                              {
                                Uf[v].Ufl_First=Uf[v].Ufl->pNext;
                                free(Uf[v].Ufl);
                                Uf[v].Ufl=Uf[v].Ufl_First;
                              }
                            FBINARY_ fbinary(oMinus);
                            fbinary.write(code_file);
                            FSTPSU_ fstpsu(i - ModelBlock->Block_List[j].Nb_Recursives);
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
StaticDllModel::Write_Inf_To_Bin_File(const string &static_basename, const string &bin_basename, const int &num,
                                    int &u_count_int, bool &file_open) const
  {
    int j;
    std::ofstream SaveCode;
    if (file_open)
      SaveCode.open((bin_basename + "_static.bin").c_str(), ios::out | ios::in | ios::binary | ios ::ate );
    else
      SaveCode.open((bin_basename + "_static.bin").c_str(), ios::out | ios::binary);
    if (!SaveCode.is_open())
      {
        cout << "Error : Can't open file \"" << bin_basename << "_static.bin\" for writing\n";
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
            //cout << "eq=" << v << ", var=" << var << ", lag=" << k << " u=" << u << "\n";
            u_count_int++;
					}
			}
		/*cout << "u_count_int=" << u_count_int << endl;
		cout << "block_triangular.ModelBlock->Block_List[" << num << "].Nb_Recursives=" << block_triangular.ModelBlock->Block_List[num].Nb_Recursives << " block_triangular.ModelBlock->Block_List[" << num << "].Size=" << block_triangular.ModelBlock->Block_List[num].Size << endl;*/
    for (j=block_triangular.ModelBlock->Block_List[num].Nb_Recursives;j<block_triangular.ModelBlock->Block_List[num].Size;j++)
      {
        int varr=block_triangular.ModelBlock->Block_List[num].Variable[j];
        //cout << "j=" << j << " varr=" << varr << "\n";
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
StaticDllModel::evaluateJacobian(const eval_context_type &eval_context, jacob_map *j_m, bool dynamic)
{
  int i=0;
  int j=0;
  bool *IM=NULL;
  int a_variable_lag=-9999;
  for (first_derivatives_type::iterator it = first_derivatives.begin();
       it != first_derivatives.end(); it++)
    {
      //cout << "it->first.second=" << it->first.second << " variable_table.getSymbolID(it->first.second)=" << variable_table.getSymbolID(it->first.second) << " Type=" << variable_table.getType(it->first.second) << " eEndogenous=" << eEndogenous << " eExogenous=" << eExogenous << " variable_table.getLag(it->first.second)=" << variable_table.getLag(it->first.second) << "\n";
      /*if (getTypeByDerivID(it->first.second) == eEndogenous)
        {*/
          NodeID Id = it->second;
          double val = 0;
          try
            {
              val = Id->eval(eval_context);
            }
          catch (ExprNode::EvalException &e)
            {
              //cout << "evaluation of Jacobian failed for equation " << it->first.first+1 << " and variable " << symbol_table.getName(getSymbIDByDerivID(it->first.second)) << "(" << getLagByDerivID(it->first.second) << ") [" << getSymbIDByDerivID(it->first.second) << "] !" << endl;
              cout << "evaluation of Jacobian failed for equation " << it->first.first+1 << " and variable " << symbol_table.getName(it->first.second) << "(" << 0 << ") [" << it->first.second << "] !" << endl;
              Id->writeOutput(cout, oMatlabStaticModelSparse, temporary_terms);
              cout << "\n";
              //cerr << "StaticDllModel::evaluateJacobian: evaluation of Jacobian failed for equation " << it->first.first+1 << " and variable " << symbol_table.getName(getSymbIDByDerivID(it->first.second)) << "(" << getLagByDerivID(it->first.second) << ")!" << endl;
              cerr << "StaticDllModel::evaluateJacobian: evaluation of Jacobian failed for equation " << it->first.first+1 << " and variable " << symbol_table.getName(it->first.second) << "(" << 0 << ")!" << endl;
            }
          int eq=it->first.first;
          //int var = symbol_table.getTypeSpecificID(getSymbIDByDerivID(it->first.second));///symbol_table.getID(eEndogenous,it->first.second);//variable_table.getSymbolID(it->first.second);
          int var = symbol_table.getTypeSpecificID(it->first.second);///symbol_table.getID(eEndogenous,it->first.second);//variable_table.getSymbolID(it->first.second);
          int k1 = 0;//getLagByDerivID(it->first.second);
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
          /*if (IM[eq*symbol_table.endo_nbr()+var] && (fabs(val) < cutoff))
            {
              //if (block_triangular.bt_verbose)
                cout << "the coefficient related to variable " << var << " with lag " << k1 << " in equation " << eq << " is equal to " << val << " and is set to 0 in the incidence matrix (size=" << symbol_table.endo_nbr() << ")\n";
              block_triangular.incidencematrix.unfill_IM(eq, var, k1, eEndogenous);
              i++;
            }*/
        /*}*/
    }
  //Get ride of the elements of the incidence matrix equal to Zero
  IM=block_triangular.incidencematrix.Get_IM(0, eEndogenous);
  /*for (int i=0;i<symbol_table.endo_nbr();i++)
    for (int j=0;j<symbol_table.endo_nbr();j++)
      if (IM[i*symbol_table.endo_nbr()+j])
        //if (first_derivatives.find(make_pair(i,getDerivID(symbol_table.getID(eEndogenous, j), 0)))==first_derivatives.end())
        if (first_derivatives.find(make_pair(i,symbol_table.getID(eEndogenous, j)))==first_derivatives.end())
          {
            block_triangular.incidencematrix.unfill_IM(i, j, 0, eEndogenous);
            cout << "eliminating equation " << i << " and variable " << j << " in the incidence matrix\n";
          }*/
  if (i>0)
    {
      cout << i << " elements among " << first_derivatives.size() << " in the incidence matrices are below the cutoff (" << cutoff << ") and are discarded\n";
      cout << "the contemporaneous incidence matrix has " << j << " elements\n";
    }
}

void
StaticDllModel::BlockLinear(Model_Block *ModelBlock)
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
              //first_derivatives_type::const_iterator it=first_derivatives.find(make_pair(eq,getDerivID(symbol_table.getID(eEndogenous, var),0)));
              first_derivatives_type::const_iterator it=first_derivatives.find(make_pair(eq,symbol_table.getID(eEndogenous, var)));
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
                  //first_derivatives_type::const_iterator it=first_derivatives.find(make_pair(eq,getDerivID(symbol_table.getID(eEndogenous, var),k1)));
                  first_derivatives_type::const_iterator it=first_derivatives.find(make_pair(eq,symbol_table.getID(eEndogenous, var)));
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
StaticDllModel::collect_first_order_derivatives_endogenous()
{
  map<pair<int, pair<int, int > >, NodeID> endo_derivatives;
  for (first_derivatives_type::iterator it2 = first_derivatives.begin();
       it2 != first_derivatives.end(); it2++)
    {
      //if (getTypeByDerivID(it2->first.second)==eEndogenous)
      /*if (it2->first.second)==eEndogenous)
        {*/
          int eq = it2->first.first;
          //int var=symbol_table.getTypeSpecificID(getSymbIDByDerivID(it2->first.second));
          int var=symbol_table.getTypeSpecificID(it2->first.second);
          //int lag=getLagByDerivID(it2->first.second);
          int lag = 0;
          //if (lag==0)
          endo_derivatives[make_pair(eq, make_pair(var, lag))] = it2->second;
        //}
    }
  return  endo_derivatives;
}



void
StaticDllModel::computingPass(const eval_context_type &eval_context, bool no_tmp_terms, bool block)
{
  assert(block);

  // Compute derivatives w.r. to all endogenous, and possibly exogenous and exogenous deterministic
  set<int> vars;
  /*for (deriv_id_table_t::const_iterator it = deriv_id_table.begin();
       it != deriv_id_table.end(); it++)
    {
      SymbolType type = symbol_table.getType(it->first.first);
      if (type == eEndogenous)
        vars.insert(it->second);
    }*/
  for(int i = 0; i < symbol_table.endo_nbr(); i++)
    vars.insert(symbol_table.getID(eEndogenous, i));

  // Launch computations
  cout << "Computing static model derivatives:" << endl
  << " - order 1" << endl;
  computeJacobian(vars);
  //cout << "mode=" << mode << " eSparseDLLMode=" << eSparseDLLMode << " eSparseMode=" << eSparseMode << "\n";
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
      /*for (int j = 0;j < block_triangular.ModelBlock->Size;j++)
        {
          for (int i = 0;i < block_triangular.ModelBlock->Block_List[j].Size;i++)
            {
        	    if (i<block_triangular.ModelBlock->Block_List[j].Nb_Recursives )
        	    	cout << "block=" << j << " R i=" << i << " equation=" << block_triangular.ModelBlock->Block_List[j].Equation[i]+1 << " variable=" << block_triangular.ModelBlock->Block_List[j].Variable[i]+1 << endl;
							else
							  cout << "block=" << j << " S i=" << i << " equation=" << block_triangular.ModelBlock->Block_List[j].Equation[i]+1 << " variable=" << block_triangular.ModelBlock->Block_List[j].Variable[i]+1 << endl;
            }
        }*/

      BlockLinear(block_triangular.ModelBlock);

      computeChainRuleJacobian(block_triangular.ModelBlock);

      if (!no_tmp_terms)
        computeTemporaryTermsOrdered(block_triangular.ModelBlock);

}

void
StaticDllModel::writeStaticFile(const string &basename, bool block) const
  {
    int r;

		assert(block);

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
    writeModelEquationsCodeOrdered(basename + "_static", block_triangular.ModelBlock, basename, map_idx);
    block_triangular.Free_Block(block_triangular.ModelBlock);
    block_triangular.incidencematrix.Free_IM();
  }

SymbolType
StaticDllModel::getTypeByDerivID(int deriv_id) const throw (UnknownDerivIDException)
{
  return symbol_table.getType(getSymbIDByDerivID(deriv_id));
}

int
StaticDllModel::getLagByDerivID(int deriv_id) const throw (UnknownDerivIDException)
{
  if (deriv_id < 0 || deriv_id >= (int) inv_deriv_id_table.size())
    throw UnknownDerivIDException();

  return inv_deriv_id_table[deriv_id].second;
}

int
StaticDllModel::getSymbIDByDerivID(int deriv_id) const throw (UnknownDerivIDException)
{
  if (deriv_id < 0 || deriv_id >= (int) inv_deriv_id_table.size())
    throw UnknownDerivIDException();

  return inv_deriv_id_table[deriv_id].first;
}

int
StaticDllModel::getDerivID(int symb_id, int lag) const throw (UnknownDerivIDException)
{
  if (symbol_table.getType(symb_id) == eEndogenous)
    return symb_id;
  else
    return -1;
}

void
StaticDllModel::computeChainRuleJacobian(Model_Block *ModelBlock)
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
                //recursive_variables[getDerivID(symbol_table.getID(eEndogenous, ModelBlock->Block_List[blck].Variable[i]), 0)] = ModelBlock->Block_List[blck].Equation_Normalized[i];
                recursive_variables[symbol_table.getID(eEndogenous, ModelBlock->Block_List[blck].Variable[i])] = ModelBlock->Block_List[blck].Equation_Normalized[i];
              else
                //recursive_variables[getDerivID(symbol_table.getID(eEndogenous, ModelBlock->Block_List[blck].Variable[i]), 0)] = equations[ModelBlock->Block_List[blck].Equation[i]];
                recursive_variables[symbol_table.getID(eEndogenous, ModelBlock->Block_List[blck].Variable[i])] = equations[ModelBlock->Block_List[blck].Equation[i]];
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
            	    //first_chain_rule_derivatives[make_pair(eqr, make_pair(varr, lag))] = first_derivatives[make_pair(eqr, getDerivID(symbol_table.getID(eEndogenous, varr), lag))];
            	    first_chain_rule_derivatives[make_pair(eqr, make_pair(varr, lag))] = first_derivatives[make_pair(eqr, symbol_table.getID(eEndogenous, varr))];
            	  }
							else if (Deriv_type == 1)
							  {
							    first_chain_rule_derivatives[make_pair(eqr, make_pair(varr, lag))] = ModelBlock->Block_List[blck].Equation_Normalized[eq]->getChainRuleDerivative(symbol_table.getID(eEndogenous, varr), recursive_variables);
							  }
							else if (Deriv_type == 2)
							  {
							  	if(ModelBlock->Block_List[blck].Equation_Type[eq] == E_EVALUATE_S && eq<ModelBlock->Block_List[blck].Nb_Recursives)
   						      first_chain_rule_derivatives[make_pair(eqr, make_pair(varr, lag))] = ModelBlock->Block_List[blck].Equation_Normalized[eq]->getChainRuleDerivative(symbol_table.getID(eEndogenous, varr), recursive_variables);
									else
									  //first_chain_rule_derivatives[make_pair(eqr, make_pair(varr, lag))] = equations[eqr]->getChainRuleDerivative(getDerivID(symbol_table.getID(eEndogenous, varr), lag), recursive_variables);
									  first_chain_rule_derivatives[make_pair(eqr, make_pair(varr, lag))] = equations[eqr]->getChainRuleDerivative(symbol_table.getID(eEndogenous, varr), recursive_variables);
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
                //recursive_variables[getDerivID(symbol_table.getID(eEndogenous, ModelBlock->Block_List[blck].Variable[i]), 0)] = ModelBlock->Block_List[blck].Equation_Normalized[i];
                recursive_variables[symbol_table.getID(eEndogenous, ModelBlock->Block_List[blck].Variable[i])] = ModelBlock->Block_List[blck].Equation_Normalized[i];
              else
                //recursive_variables[getDerivID(symbol_table.getID(eEndogenous, ModelBlock->Block_List[blck].Variable[i]), 0)] = equations[ModelBlock->Block_List[blck].Equation[i]];
                recursive_variables[symbol_table.getID(eEndogenous, ModelBlock->Block_List[blck].Variable[i])] = equations[ModelBlock->Block_List[blck].Equation[i]];
            }
          for(int eq = ModelBlock->Block_List[blck].Nb_Recursives; eq < ModelBlock->Block_List[blck].Size; eq++)
            {
              int eqr = ModelBlock->Block_List[blck].Equation[eq];
              for(int var = ModelBlock->Block_List[blck].Nb_Recursives; var < ModelBlock->Block_List[blck].Size; var++)
                {
                  int varr = ModelBlock->Block_List[blck].Variable[var];
                  //NodeID d1 = equations[eqr]->getChainRuleDerivative(getDerivID(symbol_table.getID(eEndogenous, varr), 0), recursive_variables);
                  NodeID d1 = equations[eqr]->getChainRuleDerivative(symbol_table.getID(eEndogenous, varr), recursive_variables);
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
StaticDllModel::writeChainRuleDerivative(ostream &output, int eqr, int varr, int lag,
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
StaticDllModel::writeLatexFile(const string &basename) const
  {
    writeLatexModelFile(basename + "_static.tex", oLatexStaticModel);
  }

void
StaticDllModel::jacobianHelper(ostream &output, int eq_nb, int col_nb, ExprNodeOutputType output_type) const
{
  output << LEFT_ARRAY_SUBSCRIPT(output_type);
  if (IS_MATLAB(output_type))
    output << eq_nb + 1 << ", " << col_nb + 1;
  else
    output << eq_nb + col_nb * equations.size();
  output << RIGHT_ARRAY_SUBSCRIPT(output_type);
}

void
StaticDllModel::hessianHelper(ostream &output, int row_nb, int col_nb, ExprNodeOutputType output_type) const
{
  output << LEFT_ARRAY_SUBSCRIPT(output_type);
  if (IS_MATLAB(output_type))
    output << row_nb + 1 << ", " << col_nb + 1;
  else
    output << row_nb + col_nb * NNZDerivatives[1];
  output << RIGHT_ARRAY_SUBSCRIPT(output_type);
}

void
StaticDllModel::writeAuxVarInitval(ostream &output) const
{
  for(int i = 0; i < (int) aux_equations.size(); i++)
    {
      dynamic_cast<ExprNode *>(aux_equations[i])->writeOutput(output);
      output << ";" << endl;
    }
}
