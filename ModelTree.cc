/*
 * Copyright (C) 2003-2008 Dynare Team
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
#include <sstream>
#include <cstring>
#include <cmath>

#include "ModelTree.hh"

#include "Model_Graph.hh"

ModelTree::ModelTree(SymbolTable &symbol_table_arg,
                     NumericalConstants &num_constants_arg) :
  DataTree(symbol_table_arg, num_constants_arg),
  mode(eStandardMode),
  cutoff(1e-15),
  markowitz(0.7),
  new_SGE(true),
  computeJacobian(false),
  computeJacobianExo(false),
  computeHessian(false),
  computeStaticHessian(false),
  computeThirdDerivatives(false),
  block_triangular(symbol_table_arg)
{
}

int
ModelTree::equation_number() const
{
  return(equations.size());
}

void
ModelTree::writeDerivative(ostream &output, int eq, int symb_id, int lag,
                           ExprNodeOutputType output_type,
                           const temporary_terms_type &temporary_terms) const
{
  first_derivatives_type::const_iterator it = first_derivatives.find(make_pair(eq, variable_table.getID(eEndogenous, symb_id, lag)));
  if (it != first_derivatives.end())
    (it->second)->writeOutput(output, output_type, temporary_terms);
  else
    output << 0;
}

void
ModelTree::compileDerivative(ofstream &code_file, int eq, int symb_id, int lag, ExprNodeOutputType output_type, map_idx_type map_idx) const
{
  first_derivatives_type::const_iterator it = first_derivatives.find(make_pair(eq, variable_table.getID(eEndogenous, symb_id, lag)));
  if (it != first_derivatives.end())
    {
      /*NodeID Id = it->second;*/
      (it->second)->compile(code_file,false, output_type, temporary_terms, map_idx);
    }
  else
    {
      code_file.write(&FLDZ, sizeof(FLDZ));
    }
}


void
ModelTree::derive(int order)
{
  cout << "Processing derivation ..." << endl;

  cout << "  Processing Order 1... ";
  for(int var = 0; var < variable_table.size(); var++)
    for(int eq = 0; eq < (int) equations.size(); eq++)
      {
        NodeID d1 = equations[eq]->getDerivative(var);
        if (d1 == Zero)
          continue;
        first_derivatives[make_pair(eq, var)] = d1;
      }
  cout << "done" << endl;

  if (order >= 2)
    {
      cout << "  Processing Order 2... ";
      for(first_derivatives_type::const_iterator it = first_derivatives.begin();
          it != first_derivatives.end(); it++)
        {
          int eq = it->first.first;
          int var1 = it->first.second;
          NodeID d1 = it->second;

          // Store only second derivatives with var2 <= var1
          for(int var2 = 0; var2 <= var1; var2++)
            {
              NodeID d2 = d1->getDerivative(var2);
              if (d2 == Zero)
                continue;
              second_derivatives[make_pair(eq, make_pair(var1, var2))] = d2;
            }
        }
      cout << "done" << endl;
    }

  if (order >= 3)
    {
      cout << "  Processing Order 3... ";
      for(second_derivatives_type::const_iterator it = second_derivatives.begin();
          it != second_derivatives.end(); it++)
        {
          int eq = it->first.first;

          int var1 = it->first.second.first;
          int var2 = it->first.second.second;
          // By construction, var2 <= var1

          NodeID d2 = it->second;

          // Store only third derivatives such that var3 <= var2 <= var1
          for(int var3 = 0; var3 <= var2; var3++)
            {
              NodeID d3 = d2->getDerivative(var3);
              if (d3 == Zero)
                continue;
              third_derivatives[make_pair(eq, make_pair(var1, make_pair(var2, var3)))] = d3;
            }
        }
      cout << "done" << endl;
    }
}

void
ModelTree::computeTemporaryTerms(int order)
{
  map<NodeID, int> reference_count;
  temporary_terms.clear();

  bool is_matlab = (mode != eDLLMode);

  for(vector<BinaryOpNode *>::iterator it = equations.begin();
      it != equations.end(); it++)
    (*it)->computeTemporaryTerms(reference_count, temporary_terms, is_matlab);

  for(first_derivatives_type::iterator it = first_derivatives.begin();
      it != first_derivatives.end(); it++)
    it->second->computeTemporaryTerms(reference_count, temporary_terms, is_matlab);

  if (order >= 2)
    for(second_derivatives_type::iterator it = second_derivatives.begin();
        it != second_derivatives.end(); it++)
      it->second->computeTemporaryTerms(reference_count, temporary_terms, is_matlab);

  if (order >= 3)
    for(third_derivatives_type::iterator it = third_derivatives.begin();
        it != third_derivatives.end(); it++)
      it->second->computeTemporaryTerms(reference_count, temporary_terms, is_matlab);
}

void
ModelTree::writeTemporaryTerms(ostream &output, ExprNodeOutputType output_type) const
{
  // A copy of temporary terms
  temporary_terms_type tt2;

  if (temporary_terms.size() > 0 && (!OFFSET(output_type)))
    output << "double\n";

  for(temporary_terms_type::const_iterator it = temporary_terms.begin();
      it != temporary_terms.end(); it++)
    {
      if (!OFFSET(output_type) && it != temporary_terms.begin())
        output << "," << endl;

      (*it)->writeOutput(output, output_type, temporary_terms);
      output << " = ";

      (*it)->writeOutput(output, output_type, tt2);

      // Insert current node into tt2
      tt2.insert(*it);

      if (OFFSET(output_type))
        output << ";" << endl;
    }
  if (!OFFSET(output_type))
    output << ";" << endl;
}

void
ModelTree::writeModelLocalVariables(ostream &output, ExprNodeOutputType output_type) const
{
  for(map<int, NodeID>::const_iterator it = local_variables_table.begin();
      it != local_variables_table.end(); it++)
    {
      int id = it->first;
      NodeID value = it->second;

      if (!OFFSET(output_type))
        output << "double ";

      output << symbol_table.getNameByID(eModelLocalVariable, id) << " = ";
      // Use an empty set for the temporary terms
      value->writeOutput(output, output_type, temporary_terms_type());
      output << ";" << endl;
    }
}

void
ModelTree::writeModelEquations(ostream &output, ExprNodeOutputType output_type) const
{
  for(int eq = 0; eq < (int) equations.size(); eq++)
    {
      BinaryOpNode *eq_node = equations[eq];

      NodeID lhs = eq_node->arg1;
      output << "lhs =";
      lhs->writeOutput(output, output_type, temporary_terms);
      output << ";" << endl;

      NodeID rhs = eq_node->arg2;
      output << "rhs =";
      rhs->writeOutput(output, output_type, temporary_terms);
      output << ";" << endl;

      output << "residual" << LPAR(output_type) << eq + OFFSET(output_type) << RPAR(output_type) << "= lhs-rhs;" << endl;
    }
}

void
ModelTree::computeTemporaryTermsOrdered(int order, Model_Block *ModelBlock)
{
  map<NodeID, int> reference_count, first_occurence;
  int i, j, m, eq, var, lag/*, prev_size=0*/;
  temporary_terms_type vect;
  ostringstream tmp_output;
  BinaryOpNode *eq_node;
  NodeID lhs, rhs;
  first_derivatives_type::const_iterator it;
  ostringstream tmp_s;

  temporary_terms.clear();
  map_idx.clear();
  for(j = 0;j < ModelBlock->Size;j++)
    {
      if (ModelBlock->Block_List[j].Size==1)
        {
          eq_node = equations[ModelBlock->Block_List[j].Equation[0]];
          lhs = eq_node->arg1;
          rhs = eq_node->arg2;
          tmp_s.str("");
          tmp_output.str("");
          lhs->writeOutput(tmp_output, oCDynamicModelSparseDLL, temporary_terms);
          tmp_s << "y[Per_y_+" << ModelBlock->Block_List[j].Variable[0] << "]";
          if (tmp_output.str()==tmp_s.str())
            {
              if (ModelBlock->Block_List[j].Simulation_Type==SOLVE_BACKWARD_SIMPLE)
                ModelBlock->Block_List[j].Simulation_Type=EVALUATE_BACKWARD;
              else if (ModelBlock->Block_List[j].Simulation_Type==SOLVE_FORWARD_SIMPLE)
                ModelBlock->Block_List[j].Simulation_Type=EVALUATE_FORWARD;
            }
          else
            {
              tmp_output.str("");
              rhs->writeOutput(tmp_output, oCDynamicModelSparseDLL, temporary_terms);
              if (tmp_output.str()==tmp_s.str())
                {
                  if (ModelBlock->Block_List[j].Simulation_Type==SOLVE_BACKWARD_SIMPLE)
                    ModelBlock->Block_List[j].Simulation_Type=EVALUATE_BACKWARD_R;
                  else if (ModelBlock->Block_List[j].Simulation_Type==SOLVE_FORWARD_SIMPLE)
                    ModelBlock->Block_List[j].Simulation_Type=EVALUATE_FORWARD_R;
                }
            }
        }
      for(i = 0;i < ModelBlock->Block_List[j].Size;i++)
        {
          eq_node = equations[ModelBlock->Block_List[j].Equation[i]];
          eq_node->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, j, ModelBlock, map_idx);
        }
      if (ModelBlock->Block_List[j].Simulation_Type!=EVALUATE_BACKWARD
          && ModelBlock->Block_List[j].Simulation_Type!=EVALUATE_FORWARD
          &&ModelBlock->Block_List[j].Simulation_Type!=EVALUATE_BACKWARD_R
          && ModelBlock->Block_List[j].Simulation_Type!=EVALUATE_FORWARD_R)
        {
          if (ModelBlock->Block_List[j].Simulation_Type==SOLVE_TWO_BOUNDARIES_COMPLETE ||
              ModelBlock->Block_List[j].Simulation_Type==SOLVE_TWO_BOUNDARIES_SIMPLE)
            {
              for(m=0;m<=ModelBlock->Block_List[j].Max_Lead+ModelBlock->Block_List[j].Max_Lag;m++)
                {
                  lag=m-ModelBlock->Block_List[j].Max_Lag;
                  for(i=0;i<ModelBlock->Block_List[j].IM_lead_lag[m].size;i++)
                    {
                      eq=ModelBlock->Block_List[j].IM_lead_lag[m].Equ_Index[i];
                      var=ModelBlock->Block_List[j].IM_lead_lag[m].Var_Index[i];
                      it=first_derivatives.find(make_pair(eq,variable_table.getID(eEndogenous, var,lag)));
                      it->second->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, j, ModelBlock, map_idx);
                    }
                }
            }
          else if (ModelBlock->Block_List[j].Simulation_Type!=SOLVE_BACKWARD_SIMPLE
                   && ModelBlock->Block_List[j].Simulation_Type!=SOLVE_FORWARD_SIMPLE)
            {
              m=ModelBlock->Block_List[j].Max_Lag;
              for(i=0;i<ModelBlock->Block_List[j].IM_lead_lag[m].size;i++)
                {
                  eq=ModelBlock->Block_List[j].IM_lead_lag[m].Equ_Index[i];
                  var=ModelBlock->Block_List[j].IM_lead_lag[m].Var_Index[i];
                  it=first_derivatives.find(make_pair(eq,variable_table.getID(eEndogenous,var,0)));
                  it->second->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, j, ModelBlock, map_idx);
                }
            }
          else
            {
              eq=ModelBlock->Block_List[j].Equation[0];
              var=ModelBlock->Block_List[j].Variable[0];
              it=first_derivatives.find(make_pair(eq,variable_table.getID(eEndogenous,var,0)));
              it->second->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, j, ModelBlock, map_idx);
            }
        }
    }
  if (order == 2)
    for(second_derivatives_type::iterator it = second_derivatives.begin();
        it != second_derivatives.end(); it++)
      it->second->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, 0, ModelBlock, map_idx);
  /*New*/
  j=0;
  for(temporary_terms_type::const_iterator it = temporary_terms.begin();
       it != temporary_terms.end(); it++)
    map_idx[(*it)->idx]=j++;
  /*EndNew*/
}

void
ModelTree::writeModelEquationsOrdered_M(ostream &output, Model_Block *ModelBlock, const string &dynamic_basename) const
{
  int i,j,k,m;
  string tmp_s, sps;
  ostringstream tmp_output, global_output;
  NodeID lhs=NULL, rhs=NULL;
  BinaryOpNode *eq_node;
  bool OK, lhs_rhs_done, skip_the_head;
  ostringstream Uf[symbol_table.endo_nbr];
  map<NodeID, int> reference_count;
  int prev_Simulation_Type=-1, count_derivates=0;
  temporary_terms_type::const_iterator it_temp=temporary_terms.begin();
  //----------------------------------------------------------------------
  //Temporary variables declaration
  OK=true;
  for(temporary_terms_type::const_iterator it = temporary_terms.begin();
      it != temporary_terms.end(); it++)
    {
      if (OK)
        OK=false;
      else
        tmp_output << " ";

      (*it)->writeOutput(tmp_output, oMatlabDynamicModel, temporary_terms);

      /*tmp_output << "[" << block_triangular.periods + variable_table.max_lag+variable_table.max_lead << "]";*/
    }
  global_output << "  global " << tmp_output.str() << " M_ ;\n";
  //For each block
  int gen_blocks=0;
  for(j = 0;j < ModelBlock->Size;j++)
    {
      //For a block composed of a single equation determines wether we have to evaluate or to solve the equation
      if (ModelBlock->Block_List[j].Size==1)
        {
          lhs_rhs_done=true;
          eq_node = equations[ModelBlock->Block_List[j].Equation[0]];
          lhs = eq_node->arg1;
          rhs = eq_node->arg2;
          tmp_output.str("");
          lhs->writeOutput(tmp_output, oMatlabDynamicModelSparse, temporary_terms);
        }
      else
        lhs_rhs_done=false;
      if (BlockTriangular::BlockSim(prev_Simulation_Type)==BlockTriangular::BlockSim(ModelBlock->Block_List[j].Simulation_Type)
          && (ModelBlock->Block_List[j].Simulation_Type==EVALUATE_BACKWARD
              ||ModelBlock->Block_List[j].Simulation_Type==EVALUATE_FORWARD
              ||ModelBlock->Block_List[j].Simulation_Type==EVALUATE_BACKWARD_R
              ||ModelBlock->Block_List[j].Simulation_Type==EVALUATE_FORWARD_R ))
        skip_the_head=true;
      else
        skip_the_head=false;
      if (!skip_the_head)
        {
          count_derivates=0;
          gen_blocks++;
          if (j>0)
            {
              output << "return;\n\n\n";
            }
          else
            output << "\n\n";
          if (ModelBlock->Block_List[j].Simulation_Type==EVALUATE_BACKWARD
              ||ModelBlock->Block_List[j].Simulation_Type==EVALUATE_FORWARD
              ||ModelBlock->Block_List[j].Simulation_Type==EVALUATE_BACKWARD_R
              ||ModelBlock->Block_List[j].Simulation_Type==EVALUATE_FORWARD_R)
            output << "function [y, g1, g2, g3] = " << dynamic_basename << "_" << j+1 << "(y, x, it_, jacobian_eval, g1, g2, g3)\n";
          else if (ModelBlock->Block_List[j].Simulation_Type==SOLVE_FORWARD_COMPLETE
              ||   ModelBlock->Block_List[j].Simulation_Type==SOLVE_BACKWARD_COMPLETE)
            output << "function [residual, g1, g2, g3, b] = " << dynamic_basename << "_" << j+1 << "(y, x, it_, jacobian_eval, g1, g2, g3)\n";
          else if (ModelBlock->Block_List[j].Simulation_Type==SOLVE_BACKWARD_SIMPLE
              ||   ModelBlock->Block_List[j].Simulation_Type==SOLVE_FORWARD_SIMPLE)
            output << "function [residual, g1, g2, g3, b] = " << dynamic_basename << "_" << j+1 << "(y, x, it_, g1, g2, g3, y_index, jacobian_eval)\n";
          else
            output << "function [residual, g1, g2, g3, b] = " << dynamic_basename << "_" << j+1 << "(y, x, y_kmin, y_size, periods, g1, g2, g3)\n";
          output << "  % ////////////////////////////////////////////////////////////////////////" << endl
                 << "  % //" << string("                     Block ").substr(int(log10(j + 1))) << j + 1 << " " << BlockTriangular::BlockType0(ModelBlock->Block_List[j].Type)
                 << "          //" << endl
                 << "  % //                     Simulation type "
                 << BlockTriangular::BlockSim(ModelBlock->Block_List[j].Simulation_Type) << "  //" << endl
                 << "  % ////////////////////////////////////////////////////////////////////////" << endl;
          //The Temporary terms
          output << global_output.str();
          output << "  if M_.param_nbr > 0\n";
          output << "    params =  M_.params;\n";
          output << "  end\n";
        }


      temporary_terms_type tt2;
      if(ModelBlock->Block_List[j].Simulation_Type==SOLVE_TWO_BOUNDARIES_COMPLETE)
        {
          int nze;
          for(nze=0,m=0;m<=ModelBlock->Block_List[j].Max_Lead+ModelBlock->Block_List[j].Max_Lag;m++)
            nze+=ModelBlock->Block_List[j].IM_lead_lag[m].size;
          //output << "  Jacobian_Size=" << ModelBlock->Block_List[j].Size << "*(y_kmin+" << ModelBlock->Block_List[j].Max_Lead << " +periods);\n";
          //output << "  g1=spalloc( y_size*periods, Jacobian_Size, " << nze << "*periods" << ");\n";
          output << "  for it_ = y_kmin+1:(periods+y_kmin)\n";
          output << "    Per_y_=it_*y_size;\n";
          output << "    Per_J_=(it_-y_kmin-1)*y_size;\n";
          output << "    Per_K_=(it_-1)*y_size;\n";
          sps="  ";
        }
      else
        sps="";
      if (ModelBlock->Block_List[j].Temporary_terms->size())
        output << "  " << sps << "% //Temporary variables" << endl;
      i=0;
      for(temporary_terms_type::const_iterator it = ModelBlock->Block_List[j].Temporary_terms->begin();
          it != ModelBlock->Block_List[j].Temporary_terms->end(); it++)
        {
          output << "  " <<  sps;
          (*it)->writeOutput(output, oMatlabDynamicModelSparse, temporary_terms);
          output << " = ";
          (*it)->writeOutput(output, oMatlabDynamicModelSparse, tt2);
          // Insert current node into tt2
          tt2.insert(*it);
          output << ";" << endl;
          i++;
        }
      // The equations
      for(i = 0;i < ModelBlock->Block_List[j].Size;i++)
        {
          ModelBlock->Block_List[j].Variable_Sorted[i] = variable_table.getID(eEndogenous, ModelBlock->Block_List[j].Variable[i], 0);
          string sModel = symbol_table.getNameByID(eEndogenous, ModelBlock->Block_List[j].Variable[i]) ;
          output << sps << "  % equation " << ModelBlock->Block_List[j].Equation[i]+1 << " variable : " << sModel
                 << " (" << ModelBlock->Block_List[j].Variable[i]+1 << ")" << endl;
          if (!lhs_rhs_done)
            {
              eq_node = equations[ModelBlock->Block_List[j].Equation[i]];
              lhs = eq_node->arg1;
              rhs = eq_node->arg2;
              tmp_output.str("");
              lhs->writeOutput(tmp_output, oMatlabDynamicModelSparse, temporary_terms);
            }
          output << "  ";
          switch(ModelBlock->Block_List[j].Simulation_Type)
            {
            case EVALUATE_BACKWARD:
            case EVALUATE_FORWARD:
              output << tmp_output.str();
              output << " = ";
              rhs->writeOutput(output, oMatlabDynamicModelSparse, temporary_terms);
              output << ";\n";
              break;
            case EVALUATE_BACKWARD_R:
            case EVALUATE_FORWARD_R:
              rhs->writeOutput(output, oMatlabDynamicModelSparse, temporary_terms);
              output << " = ";
              lhs->writeOutput(output, oMatlabDynamicModelSparse, temporary_terms);
              output << ";\n";
              break;
            case SOLVE_BACKWARD_SIMPLE:
            case SOLVE_FORWARD_SIMPLE:
              output << sps << "residual(" << i+1 << ") = (";
              goto end;
            case SOLVE_BACKWARD_COMPLETE:
            case SOLVE_FORWARD_COMPLETE:
              Uf[ModelBlock->Block_List[j].Equation[i]] << "  b(" << i+1 << ") = residual(" << i+1 << ")";
              output << sps << "residual(" << i+1 << ") = (";
              goto end;
            case SOLVE_TWO_BOUNDARIES_COMPLETE:
              Uf[ModelBlock->Block_List[j].Equation[i]] << "    b(" << i+1 << "+Per_J_) = -residual(" << i+1 << ", it_)";
              output << sps << "residual(" << i+1 << ", it_) = (";
              goto end;
            default:
            end:
              output << tmp_output.str();
              output << ") - (";
              rhs->writeOutput(output, oMatlabDynamicModelSparse, temporary_terms);
              output << ");\n";
#ifdef CONDITION
              if (ModelBlock->Block_List[j].Simulation_Type==SOLVE_TWO_BOUNDARIES_COMPLETE)
                output << "  condition(" << i+1 << ")=0;\n";
#endif
            }
        }
      // The Jacobian if we have to solve the block
      if (ModelBlock->Block_List[j].Simulation_Type==SOLVE_TWO_BOUNDARIES_SIMPLE
          ||  ModelBlock->Block_List[j].Simulation_Type==SOLVE_TWO_BOUNDARIES_COMPLETE)
          output << "  " << sps << "% Jacobian  " << endl;
      else
          output << "  " << sps << "% Jacobian  " << endl << "  if jacobian_eval" << endl;
      switch(ModelBlock->Block_List[j].Simulation_Type)
        {
        case SOLVE_BACKWARD_SIMPLE:
        case SOLVE_FORWARD_SIMPLE:
        case EVALUATE_BACKWARD:
        case EVALUATE_FORWARD:
        case EVALUATE_BACKWARD_R:
        case EVALUATE_FORWARD_R:
          count_derivates++;
          for(m=0;m<ModelBlock->Block_List[j].Max_Lead+ModelBlock->Block_List[j].Max_Lag+1;m++)
            {
              k=m-ModelBlock->Block_List[j].Max_Lag;
              for(i=0;i<ModelBlock->Block_List[j].IM_lead_lag[m].size;i++)
                {
                  if(ModelBlock->Block_List[j].IM_lead_lag[m].Var_Index[i]==ModelBlock->Block_List[j].Variable[0])
                    {
                      //output << "    g1(M_.block_structure.block(" << gen_blocks << ").equation(" << count_derivates << "), M_.block_structure.block(" << gen_blocks << ").variable(" << count_derivates << ")+" << (m+variable_table.max_endo_lag-ModelBlock->Block_List[j].Max_Lag)*symbol_table.endo_nbr << ")=";
                      //output << "    g1(M_.block_structure.block(" << gen_blocks << ").equation(" << count_derivates << "), M_.block_structure.block(" << gen_blocks << ").variable(" << count_derivates << ")+" << (m+variable_table.max_endo_lag-ModelBlock->Block_List[j].Max_Lag)*symbol_table.endo_nbr << ")=";
                      output << "    g1(" << ModelBlock->Block_List[j].Equation[0]+1 << ", " << ModelBlock->Block_List[j].Variable[0]+1 + (m+variable_table.max_endo_lag-ModelBlock->Block_List[j].Max_Lag)*symbol_table.endo_nbr << ")=";
                      writeDerivative(output, ModelBlock->Block_List[j].Equation[0], ModelBlock->Block_List[j].Variable[0], k, oMatlabDynamicModelSparse, temporary_terms);
                      output << "; % variable=" << symbol_table.getNameByID(eEndogenous, ModelBlock->Block_List[j].Variable[0])
                             << "(" << variable_table.getLag(variable_table.getSymbolID(ModelBlock->Block_List[j].Variable[0]))
                             << ") " << ModelBlock->Block_List[j].Variable[0]+1
                             << ", equation=" << ModelBlock->Block_List[j].Equation[0]+1 << endl;
                    }
                }
            }
          if (ModelBlock->Block_List[j].Simulation_Type==SOLVE_BACKWARD_SIMPLE
          || ModelBlock->Block_List[j].Simulation_Type==SOLVE_FORWARD_SIMPLE)
            {
              output << "  else\n";
              output << "    g1=";
              writeDerivative(output, ModelBlock->Block_List[j].Equation[0], ModelBlock->Block_List[j].Variable[0], 0, oMatlabDynamicModelSparse, temporary_terms);
              output << "; % variable=" << symbol_table.getNameByID(eEndogenous, ModelBlock->Block_List[j].Variable[0])
                     << "(" << variable_table.getLag(variable_table.getSymbolID(ModelBlock->Block_List[j].Variable[0]))
                     << ") " << ModelBlock->Block_List[j].Variable[0]+1
                     << ", equation=" << ModelBlock->Block_List[j].Equation[0]+1 << endl;
            }
          if (ModelBlock->Block_List[j].Simulation_Type==EVALUATE_BACKWARD
          || ModelBlock->Block_List[j].Simulation_Type==EVALUATE_FORWARD
          || ModelBlock->Block_List[j].Simulation_Type==EVALUATE_BACKWARD_R
          || ModelBlock->Block_List[j].Simulation_Type==EVALUATE_FORWARD_R
          || ModelBlock->Block_List[j].Simulation_Type==SOLVE_BACKWARD_SIMPLE
          || ModelBlock->Block_List[j].Simulation_Type==SOLVE_FORWARD_SIMPLE)
            output << "  end;" << endl;
          break;
        case SOLVE_BACKWARD_COMPLETE:
        case SOLVE_FORWARD_COMPLETE:
          count_derivates++;
          for(m=0;m<ModelBlock->Block_List[j].Max_Lead+ModelBlock->Block_List[j].Max_Lag+1;m++)
            {
              k=m-ModelBlock->Block_List[j].Max_Lag;
              for(i=0;i<ModelBlock->Block_List[j].IM_lead_lag[m].size;i++)
                {
                  int eq=ModelBlock->Block_List[j].IM_lead_lag[m].Equ_Index[i];
                  int var=ModelBlock->Block_List[j].IM_lead_lag[m].Var_Index[i];
                  int varr=ModelBlock->Block_List[j].IM_lead_lag[m].Var[i];
                  int eqr=ModelBlock->Block_List[j].IM_lead_lag[m].Equ[i];
                  output << "    g1(" << eqr+1 << ", " << varr+1+m*ModelBlock->Block_List[j].Size << ")=";
                  writeDerivative(output, eq, var, k, oMatlabDynamicModelSparse, temporary_terms);
                  output << "; % variable=" << symbol_table.getNameByID(eEndogenous, var)
                         << "(" << /*variable_table.getLag(variable_table.getSymbolID(var))*/k
                         << ") " << var+1
                         << ", equation=" << eq+1 << endl;
                }
            }
          output << "  else\n";
          m=ModelBlock->Block_List[j].Max_Lag;
          for(i=0;i<ModelBlock->Block_List[j].IM_lead_lag[m].size;i++)
            {
              int eq=ModelBlock->Block_List[j].IM_lead_lag[m].Equ_Index[i];
              int var=ModelBlock->Block_List[j].IM_lead_lag[m].Var_Index[i];
              int eqr=ModelBlock->Block_List[j].IM_lead_lag[m].Equ[i];
              int varr=ModelBlock->Block_List[j].IM_lead_lag[m].Var[i];
              //Uf[ModelBlock->Block_List[j].Equation[eqr]] << "-u(" << u << ")*y(Per_y_+" << var << ")";
              Uf[ModelBlock->Block_List[j].Equation[eqr]] << "+g1(" << eqr+1 << ", " << varr+1 << ")*y(it_, " << var+1 << ")";
              //output << "  u(" << u+1 << ") = ";
              output << "    g1(" << eqr+1 << ", " << varr+1 << ") = ";
              writeDerivative(output, eq, var, 0, oMatlabDynamicModelSparse, temporary_terms);
              output << "; % variable=" << symbol_table.getNameByID(eEndogenous, var)
                     << "(" << variable_table.getLag(variable_table.getSymbolID(var)) << ") " << var+1
                     << ", equation=" << eq+1 << endl;
            }
          output << "  end;\n";
          for(i = 0;i < ModelBlock->Block_List[j].Size;i++)
            output << Uf[ModelBlock->Block_List[j].Equation[i]].str() << ";\n";
          break;
        case SOLVE_TWO_BOUNDARIES_SIMPLE:
        case SOLVE_TWO_BOUNDARIES_COMPLETE:
          output << "    g2=0;g3=0;\n";
          for(m=0;m<=ModelBlock->Block_List[j].Max_Lead+ModelBlock->Block_List[j].Max_Lag;m++)
            {
              k=m-ModelBlock->Block_List[j].Max_Lag;
              for(i=0;i<ModelBlock->Block_List[j].IM_lead_lag[m].size;i++)
                {
                  int eq=ModelBlock->Block_List[j].IM_lead_lag[m].Equ_Index[i];
                  int var=ModelBlock->Block_List[j].IM_lead_lag[m].Var_Index[i];
                  //int u=ModelBlock->Block_List[j].IM_lead_lag[m].u[i];
                  int eqr=ModelBlock->Block_List[j].IM_lead_lag[m].Equ[i];
                  int varr=ModelBlock->Block_List[j].IM_lead_lag[m].Var[i];
                  if (k==0)
                    Uf[ModelBlock->Block_List[j].Equation[eqr]] << "+g1(" << eqr+1 << "+Per_J_, " << varr+1 << "+Per_K_)*y(it_, " << var+1 << ")";
                  else if (k==1)
                    Uf[ModelBlock->Block_List[j].Equation[eqr]] << "+g1(" << eqr+1 << "+Per_J_, " << varr+1 << "+Per_y_)*y(it_+1, " << var+1 << ")";
                  else if (k>0)
                    Uf[ModelBlock->Block_List[j].Equation[eqr]] << "+g1(" << eqr+1 << "+Per_J_, " << varr+1 << "+y_size*(it_+" << k-1 << "))*y(it_+" << k << ", " << var+1 << ")";
                  else if (k<0)
                    Uf[ModelBlock->Block_List[j].Equation[eqr]] << "+g1(" << eqr+1 << "+Per_J_, " << varr+1 << "+y_size*(it_" << k-1 << "))*y(it_" << k << ", " << var+1 << ")";
                  //output << "  u(" << u+1 << "+Per_u_) = ";
                  if(k==0)
                    output << "    g1(" << eqr+1 << "+Per_J_, " << varr+1 << "+Per_K_) = ";
                  else if(k==1)
                    output << "    g1(" << eqr+1 << "+Per_J_, " << varr+1 << "+Per_y_) = ";
                  else if(k>0)
                    output << "    g1(" << eqr+1 << "+Per_J_, " << varr+1 << "+y_size*(it_+" << k-1 << ")) = ";
                  else if(k<0)
                    output << "    g1(" << eqr+1 << "+Per_J_, " << varr+1 << "+y_size*(it_" << k-1 << ")) = ";
                  writeDerivative(output, eq, var, k, oMatlabDynamicModelSparse, temporary_terms);
                  output << "; % variable=" << symbol_table.getNameByID(eEndogenous, var)
                         << "(" << k << ") " << var+1
                         << ", equation=" << eq+1 << endl;
#ifdef CONDITION
                  output << "  if (fabs(condition[" << eqr << "])<fabs(u[" << u << "+Per_u_]))\n";
                  output << "    condition(" << eqr << ")=u(" << u << "+Per_u_);\n";
#endif
                }
            }
          for(i = 0;i < ModelBlock->Block_List[j].Size;i++)
            {
              output << Uf[ModelBlock->Block_List[j].Equation[i]].str() << ";\n";
#ifdef CONDITION
              output << "  if (fabs(condition(" << i+1 << "))<fabs(u(" << i << "+Per_u_)))\n";
              output << "    condition(" << i+1 << ")=u(" << i+1 << "+Per_u_);\n";
#endif
            }
#ifdef CONDITION
          for(m=0;m<=ModelBlock->Block_List[j].Max_Lead+ModelBlock->Block_List[j].Max_Lag;m++)
            {
              k=m-ModelBlock->Block_List[j].Max_Lag;
              for(i=0;i<ModelBlock->Block_List[j].IM_lead_lag[m].size;i++)
                {
                  int eq=ModelBlock->Block_List[j].IM_lead_lag[m].Equ_Index[i];
                  int var=ModelBlock->Block_List[j].IM_lead_lag[m].Var_Index[i];
                  int u=ModelBlock->Block_List[j].IM_lead_lag[m].u[i];
                  int eqr=ModelBlock->Block_List[j].IM_lead_lag[m].Equ[i];
                  output << "  u(" << u+1 << "+Per_u_) = u(" << u+1 << "+Per_u_) / condition(" << eqr+1 << ");\n";
                }
            }
          for(i = 0;i < ModelBlock->Block_List[j].Size;i++)
            output << "  u(" << i+1 << "+Per_u_) = u(" << i+1 << "+Per_u_) / condition(" << i+1 << ");\n";
#endif
          output << "  end;\n";
          break;
        }
      prev_Simulation_Type=ModelBlock->Block_List[j].Simulation_Type;
    }
  output << "return;\n\n\n";
}

void
ModelTree::writeModelStaticEquationsOrdered_M(ostream &output, Model_Block *ModelBlock, const string &static_basename) const
{
  int i,j,k,m, var, eq;
  string tmp_s, sps;
  ostringstream tmp_output, global_output;
  NodeID lhs=NULL, rhs=NULL;
  BinaryOpNode *eq_node;
  bool OK, lhs_rhs_done, skip_the_head;
  ostringstream Uf[symbol_table.endo_nbr];
  map<NodeID, int> reference_count;
  int prev_Simulation_Type=-1;
  int nze=0;
  bool *IM, *IMl;
  temporary_terms_type::const_iterator it_temp=temporary_terms.begin();
  //----------------------------------------------------------------------
  //Temporary variables declaration
  OK=true;
  for(temporary_terms_type::const_iterator it = temporary_terms.begin();
      it != temporary_terms.end(); it++)
    {
      if (OK)
        OK=false;
      else
        tmp_output << " ";
      (*it)->writeOutput(tmp_output, oMatlabStaticModelSparse, temporary_terms);
    }
  global_output << "  global " << tmp_output.str() << " M_ ;\n";
  //For each block
  for(j = 0;j < ModelBlock->Size;j++)
    {
      //For a block composed of a single equation determines wether we have to evaluate or to solve the equation
      if (ModelBlock->Block_List[j].Size==1)
        {
          lhs_rhs_done=true;
          eq_node = equations[ModelBlock->Block_List[j].Equation[0]];
          lhs = eq_node->arg1;
          rhs = eq_node->arg2;
          tmp_output.str("");
          lhs->writeOutput(tmp_output, oMatlabStaticModelSparse, temporary_terms);
        }
      else
        lhs_rhs_done=false;
      if (BlockTriangular::BlockSim(prev_Simulation_Type)==BlockTriangular::BlockSim(ModelBlock->Block_List[j].Simulation_Type)
          && (ModelBlock->Block_List[j].Simulation_Type==EVALUATE_BACKWARD
              ||ModelBlock->Block_List[j].Simulation_Type==EVALUATE_FORWARD
              ||ModelBlock->Block_List[j].Simulation_Type==EVALUATE_BACKWARD_R
              ||ModelBlock->Block_List[j].Simulation_Type==EVALUATE_FORWARD_R ))
        skip_the_head=true;
      else
        skip_the_head=false;
      if (!skip_the_head)
        {
          if (j>0)
            {
              output << "return;\n\n\n";
            }
          else
            output << "\n\n";
          if(ModelBlock->Block_List[j].Simulation_Type==EVALUATE_BACKWARD
           ||ModelBlock->Block_List[j].Simulation_Type==EVALUATE_FORWARD
           ||ModelBlock->Block_List[j].Simulation_Type==EVALUATE_BACKWARD_R
           ||ModelBlock->Block_List[j].Simulation_Type==EVALUATE_FORWARD_R )
            output << "function [y] = " << static_basename << "_" << j+1 << "(y, x)\n";
          else
            output << "function [residual, g1, g2, g3, b] = " << static_basename << "_" << j+1 << "(y, x)\n";
          output << "  % ////////////////////////////////////////////////////////////////////////" << endl
                 << "  % //" << string("                     Block ").substr(int(log10(j + 1))) << j + 1 << " "
                 << BlockTriangular::BlockType0(ModelBlock->Block_List[j].Type) << "          //" << endl
                 << "  % //                     Simulation type ";
          output << BlockTriangular::BlockSim(ModelBlock->Block_List[j].Simulation_Type) << "  //" << endl
                 << "  % ////////////////////////////////////////////////////////////////////////" << endl;
          //The Temporary terms
          output << global_output.str();
          output << "  if M_.param_nbr > 0\n";
          output << "    params =  M_.params;\n";
          output << "  end\n";
        }

      temporary_terms_type tt2;

      int n=ModelBlock->Block_List[j].Size;
      int n1=symbol_table.endo_nbr;
      IM=(bool*)malloc(n*n*sizeof(bool));
      memset(IM, 0, n*n*sizeof(bool));
      for(m=-ModelBlock->Block_List[j].Max_Lag;m<=ModelBlock->Block_List[j].Max_Lead;m++)
        {
          IMl=block_triangular.bGet_IM(m);
          for(i=0;i<n;i++)
            {
              eq=ModelBlock->Block_List[j].Equation[i];
              for(k=0;k<n;k++)
                {
                  var=ModelBlock->Block_List[j].Variable[k];
                  IM[i*n+k]=IM[i*n+k] || IMl[eq*n1+var];
                }
            }
        }
      for(nze=0, i=0;i<n*n;i++)
        {
          nze+=IM[i];
        }
      memset(IM, 0, n*n*sizeof(bool));
      if( ModelBlock->Block_List[j].Simulation_Type!=EVALUATE_BACKWARD && ModelBlock->Block_List[j].Simulation_Type!=EVALUATE_FORWARD
       && ModelBlock->Block_List[j].Simulation_Type!=EVALUATE_BACKWARD_R && ModelBlock->Block_List[j].Simulation_Type!=EVALUATE_FORWARD_R)
         {
          output << "  g1=spalloc(" << ModelBlock->Block_List[j].Size << ", " << ModelBlock->Block_List[j].Size << ", " << nze << ");\n";
          output << "  residual=zeros(" << ModelBlock->Block_List[j].Size << ",1);\n";
         }
      sps="";
      if (ModelBlock->Block_List[j].Temporary_terms->size())
        output << "  " << sps << "% //Temporary variables" << endl;
      i=0;
      for(temporary_terms_type::const_iterator it = ModelBlock->Block_List[j].Temporary_terms->begin();
          it != ModelBlock->Block_List[j].Temporary_terms->end(); it++)
        {
          output << "  " <<  sps;
          (*it)->writeOutput(output, oMatlabStaticModelSparse, temporary_terms);
          output << " = ";
          (*it)->writeOutput(output, oMatlabStaticModelSparse, tt2);
          // Insert current node into tt2
          tt2.insert(*it);
          output << ";" << endl;
          i++;
        }
      // The equations
      for(i = 0;i < ModelBlock->Block_List[j].Size;i++)
        {
          ModelBlock->Block_List[j].Variable_Sorted[i] = variable_table.getID(eEndogenous, ModelBlock->Block_List[j].Variable[i], 0);
          string sModel = symbol_table.getNameByID(eEndogenous, ModelBlock->Block_List[j].Variable[i]) ;
          output << sps << "  % equation " << ModelBlock->Block_List[j].Equation[i]+1 << " variable : "
                 << sModel << " (" << ModelBlock->Block_List[j].Variable[i]+1 << ")" << endl;
          if (!lhs_rhs_done)
            {
              eq_node = equations[ModelBlock->Block_List[j].Equation[i]];
              lhs = eq_node->arg1;
              rhs = eq_node->arg2;
              tmp_output.str("");
              lhs->writeOutput(tmp_output, oMatlabStaticModelSparse, temporary_terms);
            }
          output << "  ";
          switch(ModelBlock->Block_List[j].Simulation_Type)
            {
            case EVALUATE_BACKWARD:
            case EVALUATE_FORWARD:
              output << tmp_output.str();
              output << " = ";
              rhs->writeOutput(output, oMatlabStaticModelSparse, temporary_terms);
              output << ";\n";
              break;
            case EVALUATE_BACKWARD_R:
            case EVALUATE_FORWARD_R:
              rhs->writeOutput(output, oMatlabStaticModelSparse, temporary_terms);
              output << " = ";
              lhs->writeOutput(output, oMatlabStaticModelSparse, temporary_terms);
              output << ";\n";
              break;
            case SOLVE_BACKWARD_COMPLETE:
            case SOLVE_FORWARD_COMPLETE:
            case SOLVE_TWO_BOUNDARIES_COMPLETE:
              Uf[ModelBlock->Block_List[j].Equation[i]] << "  b(" << i+1 << ") = - residual(" << i+1 << ")";
              goto end;
            default:
            end:
              output << sps << "residual(" << i+1 << ") = (";
              output << tmp_output.str();
              output << ") - (";
              rhs->writeOutput(output, oMatlabStaticModelSparse, temporary_terms);
              output << ");\n";
#ifdef CONDITION
              if (ModelBlock->Block_List[j].Simulation_Type==SOLVE_TWO_BOUNDARIES_COMPLETE)
                output << "  condition(" << i+1 << ")=0;\n";
#endif
            }
        }
      // The Jacobian if we have to solve the block
      if (ModelBlock->Block_List[j].Simulation_Type!=EVALUATE_BACKWARD
          && ModelBlock->Block_List[j].Simulation_Type!=EVALUATE_FORWARD
          && ModelBlock->Block_List[j].Simulation_Type!=EVALUATE_BACKWARD_R
          && ModelBlock->Block_List[j].Simulation_Type!=EVALUATE_FORWARD_R)
        {
          output << "  " << sps << "% Jacobian  " << endl;
          switch(ModelBlock->Block_List[j].Simulation_Type)
            {
            case SOLVE_BACKWARD_SIMPLE:
            case SOLVE_FORWARD_SIMPLE:
              output << "  g1(1)=";
              writeDerivative(output, ModelBlock->Block_List[j].Equation[0], ModelBlock->Block_List[j].Variable[0], 0, oMatlabStaticModelSparse, temporary_terms);
              output << "; % variable=" << symbol_table.getNameByID(eEndogenous, ModelBlock->Block_List[j].Variable[0])
                     << "(" << variable_table.getLag(variable_table.getSymbolID(ModelBlock->Block_List[j].Variable[0]))
                     << ") " << ModelBlock->Block_List[j].Variable[0]+1
                     << ", equation=" << ModelBlock->Block_List[j].Equation[0]+1 << endl;
              break;
            case SOLVE_BACKWARD_COMPLETE:
            case SOLVE_FORWARD_COMPLETE:
              output << "  g2=0;g3=0;\n";
              m=ModelBlock->Block_List[j].Max_Lag;
              for(i=0;i<ModelBlock->Block_List[j].IM_lead_lag[m].size;i++)
                {
                  int eq=ModelBlock->Block_List[j].IM_lead_lag[m].Equ_Index[i];
                  int var=ModelBlock->Block_List[j].IM_lead_lag[m].Var_Index[i];
                  //int u=ModelBlock->Block_List[j].IM_lead_lag[m].us[i];
                  int eqr=ModelBlock->Block_List[j].IM_lead_lag[m].Equ[i];
                  int varr=ModelBlock->Block_List[j].IM_lead_lag[m].Var[i];
                  //Uf[ModelBlock->Block_List[j].Equation[eqr]] << "-u(" << u << ")*y(Per_y_+" << var << ")";
                  Uf[ModelBlock->Block_List[j].Equation[eqr]] << "-g1(" << eqr+1 << ", " << varr+1 << ")*y(" << var+1 << ")";
                  //output << "  u(" << u+1 << ") = ";
                  output << "  g1(" << eqr+1 << ", " << varr+1 << ") = ";
                  writeDerivative(output, eq, var, 0, oMatlabStaticModelSparse, temporary_terms);
                  output << "; % variable=" << symbol_table.getNameByID(eEndogenous, var)
                         << "(" << variable_table.getLag(variable_table.getSymbolID(var)) << ") " << var+1
                         << ", equation=" << eq+1 << endl;
                }
              for(i = 0;i < ModelBlock->Block_List[j].Size;i++)
                output << Uf[ModelBlock->Block_List[j].Equation[i]].str() << ";\n";
              break;
            case SOLVE_TWO_BOUNDARIES_COMPLETE:
              output << "  g2=0;g3=0;\n";
              for(m=0;m<=ModelBlock->Block_List[j].Max_Lead+ModelBlock->Block_List[j].Max_Lag;m++)
                {
                  k=m-ModelBlock->Block_List[j].Max_Lag;
                  for(i=0;i<ModelBlock->Block_List[j].IM_lead_lag[m].size;i++)
                    {
                      int eq=ModelBlock->Block_List[j].IM_lead_lag[m].Equ_Index[i];
                      int var=ModelBlock->Block_List[j].IM_lead_lag[m].Var_Index[i];
                      //int u=ModelBlock->Block_List[j].IM_lead_lag[m].u[i];
                      int eqr=ModelBlock->Block_List[j].IM_lead_lag[m].Equ[i];
                      int varr=ModelBlock->Block_List[j].IM_lead_lag[m].Var[i];
                      //output << "% i=" << i << " eq=" << eq << " var=" << var << " eqr=" << eqr << " varr=" << varr << "\n";
                      if(!IM[eqr*ModelBlock->Block_List[j].Size+varr])
                        {
                          Uf[ModelBlock->Block_List[j].Equation[eqr]] << "+g1(" << eqr+1
                                                                        << ", " << varr+1 << ")*y( " << var+1 << ")";
                          IM[eqr*ModelBlock->Block_List[j].Size+varr]=1;
                        }
                      output << "  g1(" << eqr+1 << ", " << varr+1 << ") = g1(" << eqr+1 << ", " << varr+1 << ") + ";
                      writeDerivative(output, eq, var, k, oMatlabStaticModelSparse, temporary_terms);
                      output << "; % variable=" << symbol_table.getNameByID(eEndogenous, var)
                             << "(" << k << ") " << var+1
                             << ", equation=" << eq+1 << endl;
#ifdef CONDITION
                      output << "  if (fabs(condition[" << eqr << "])<fabs(u[" << u << "+Per_u_]))\n";
                      output << "    condition(" << eqr << ")=u(" << u << "+Per_u_);\n";
#endif
                    }
                }
              for(i = 0;i < ModelBlock->Block_List[j].Size;i++)
                {
                  output << Uf[ModelBlock->Block_List[j].Equation[i]].str() << ";\n";
#ifdef CONDITION
                  output << "  if (fabs(condition(" << i+1 << "))<fabs(u(" << i << "+Per_u_)))\n";
                  output << "    condition(" << i+1 << ")=u(" << i+1 << "+Per_u_);\n";
#endif
                }
#ifdef CONDITION
              for(m=0;m<=ModelBlock->Block_List[j].Max_Lead+ModelBlock->Block_List[j].Max_Lag;m++)
                {
                  k=m-ModelBlock->Block_List[j].Max_Lag;
                  for(i=0;i<ModelBlock->Block_List[j].IM_lead_lag[m].size;i++)
                    {
                      int eq=ModelBlock->Block_List[j].IM_lead_lag[m].Equ_Index[i];
                      int var=ModelBlock->Block_List[j].IM_lead_lag[m].Var_Index[i];
                      int u=ModelBlock->Block_List[j].IM_lead_lag[m].u[i];
                      int eqr=ModelBlock->Block_List[j].IM_lead_lag[m].Equ[i];
                      output << "  u(" << u+1 << "+Per_u_) = u(" << u+1 << "+Per_u_) / condition(" << eqr+1 << ");\n";
                    }
                }
              for(i = 0;i < ModelBlock->Block_List[j].Size;i++)
                output << "  u(" << i+1 << "+Per_u_) = u(" << i+1 << "+Per_u_) / condition(" << i+1 << ");\n";
#endif
              break;
            }
        }
      prev_Simulation_Type=ModelBlock->Block_List[j].Simulation_Type;
      free(IM);
    }
  output << "return;\n\n\n";
  //free(IM);
}


void
ModelTree::writeModelEquationsCodeOrdered(const string file_name, const Model_Block *ModelBlock, const string bin_basename, ExprNodeOutputType output_type) const
  {
    struct Uff_l
      {
        int u, var, lag;
        Uff_l *pNext;
      };

    struct Uff
      {
        Uff_l *Ufl, *Ufl_First;
        int eqr;
      };

    int i,j,k,m, v, ModelBlock_Aggregated_Count, k0, k1;
    string tmp_s;
    ostringstream tmp_output;
    ofstream code_file;
    NodeID lhs=NULL, rhs=NULL;
    BinaryOpNode *eq_node;
    bool lhs_rhs_done;
    Uff Uf[symbol_table.endo_nbr];
    map<NodeID, int> reference_count;
    map<int,int> ModelBlock_Aggregated_Size, ModelBlock_Aggregated_Number;
    int prev_Simulation_Type=-1;
    SymbolicGaussElimination SGE;
    temporary_terms_type::const_iterator it_temp=temporary_terms.begin();
    //----------------------------------------------------------------------
    string main_name=file_name;
    main_name+=".cod";
    code_file.open(main_name.c_str(), ios::out | ios::binary | ios::ate );
    if (!code_file.is_open())
      {
        cout << "Error : Can't open file \"" << main_name << "\" for writing\n";
        exit( -1);
      }
    //Temporary variables declaration
    code_file.write(&FDIMT, sizeof(FDIMT));
    k=temporary_terms.size();
    code_file.write(reinterpret_cast<char *>(&k),sizeof(k));
    //search for successive and identical blocks
    i=k=k0=0;
    ModelBlock_Aggregated_Count=-1;
    for(j = 0;j < ModelBlock->Size;j++)
      {
        if (BlockTriangular::BlockSim(prev_Simulation_Type)==BlockTriangular::BlockSim(ModelBlock->Block_List[j].Simulation_Type)
              && (ModelBlock->Block_List[j].Simulation_Type==EVALUATE_BACKWARD
                  ||ModelBlock->Block_List[j].Simulation_Type==EVALUATE_FORWARD
                  ||ModelBlock->Block_List[j].Simulation_Type==EVALUATE_BACKWARD_R
                  ||ModelBlock->Block_List[j].Simulation_Type==EVALUATE_FORWARD_R ))
          {
          }
        else
          {
            k=k0=0;
            ModelBlock_Aggregated_Count++;
          }
        k0+=ModelBlock->Block_List[j].Size;
        ModelBlock_Aggregated_Number[ModelBlock_Aggregated_Count]=k0;
        ModelBlock_Aggregated_Size[ModelBlock_Aggregated_Count]=++k;
        prev_Simulation_Type=ModelBlock->Block_List[j].Simulation_Type;
      }
    ModelBlock_Aggregated_Count++;
    //cout << "ModelBlock_Aggregated_Count=" << ModelBlock_Aggregated_Count << "\n";
    //For each block
    j=0;
    for(k0 = 0;k0 < ModelBlock_Aggregated_Count;k0++)
      {
        k1=j;
        if (k0>0)
          code_file.write(&FENDBLOCK, sizeof(FENDBLOCK));
        code_file.write(&FBEGINBLOCK, sizeof(FBEGINBLOCK));
        v=ModelBlock_Aggregated_Number[k0];
        code_file.write(reinterpret_cast<char *>(&v),sizeof(v));
        v=ModelBlock->Block_List[j].Simulation_Type;
        code_file.write(reinterpret_cast<char *>(&v),sizeof(v));
        //cout << "FBEGINBLOCK j=" << j << " size=" << ModelBlock_Aggregated_Number[k0] << " type=" << v << "\n";
        for(k=0; k<ModelBlock_Aggregated_Size[k0]; k++)
          {
            for(i=0; i < ModelBlock->Block_List[j].Size;i++)
              {
                code_file.write(reinterpret_cast<char *>(&ModelBlock->Block_List[j].Variable[i]),sizeof(ModelBlock->Block_List[j].Variable[i]));
                code_file.write(reinterpret_cast<char *>(&ModelBlock->Block_List[j].Equation[i]),sizeof(ModelBlock->Block_List[j].Equation[i]));
                code_file.write(reinterpret_cast<char *>(&ModelBlock->Block_List[j].Own_Derivative[i]),sizeof(ModelBlock->Block_List[j].Own_Derivative[i]));
              }
            j++;
          }
        j=k1;
        if (ModelBlock->Block_List[j].Simulation_Type==SOLVE_TWO_BOUNDARIES_SIMPLE || ModelBlock->Block_List[j].Simulation_Type==SOLVE_TWO_BOUNDARIES_COMPLETE ||
            ModelBlock->Block_List[j].Simulation_Type==SOLVE_BACKWARD_COMPLETE || ModelBlock->Block_List[j].Simulation_Type==SOLVE_FORWARD_COMPLETE)
          {
            code_file.write(reinterpret_cast<char *>(&ModelBlock->Block_List[j].is_linear),sizeof(ModelBlock->Block_List[j].is_linear));
            v=block_triangular.ModelBlock->Block_List[j].IM_lead_lag[block_triangular.ModelBlock->Block_List[j].Max_Lag + block_triangular.ModelBlock->Block_List[j].Max_Lead].u_finish + 1;
            code_file.write(reinterpret_cast<char *>(&v),sizeof(v));
            v=symbol_table.endo_nbr;
            code_file.write(reinterpret_cast<char *>(&v),sizeof(v));
            v=block_triangular.ModelBlock->Block_List[j].Max_Lag;
            code_file.write(reinterpret_cast<char *>(&v),sizeof(v));
            v=block_triangular.ModelBlock->Block_List[j].Max_Lead;
            code_file.write(reinterpret_cast<char *>(&v),sizeof(v));
            if (ModelBlock->Block_List[j].Simulation_Type==SOLVE_TWO_BOUNDARIES_COMPLETE)
              {
                int u_count_int=0;
                Write_Inf_To_Bin_File(file_name, bin_basename, j, u_count_int,SGE.file_open);
                v=u_count_int;
                code_file.write(reinterpret_cast<char *>(&v),sizeof(v));
                SGE.file_is_open();
              }
          }
        for(k1 = 0; k1 < ModelBlock_Aggregated_Size[k0]; k1++)
          {
            //For a block composed of a single equation determines wether we have to evaluate or to solve the equation
            if (ModelBlock->Block_List[j].Size==1)
              {
                lhs_rhs_done=true;
                eq_node = equations[ModelBlock->Block_List[j].Equation[0]];
                lhs = eq_node->arg1;
                rhs = eq_node->arg2;
              }
            else
              lhs_rhs_done=false;
            /*if (ModelBlock->Block_List[j].Size==1)
              lhs_rhs_done=true;
            else
              lhs_rhs_done=false;*/
            //The Temporary terms
            temporary_terms_type tt2;
            i=0;
            for(temporary_terms_type::const_iterator it = ModelBlock->Block_List[j].Temporary_terms->begin();
                 it != ModelBlock->Block_List[j].Temporary_terms->end(); it++)
              {
                (*it)->compile(code_file,false, output_type, tt2, map_idx);
                code_file.write(&FSTPT, sizeof(FSTPT));
                map_idx_type::const_iterator ii=map_idx.find((*it)->idx);
                v=(int)ii->second;
                code_file.write(reinterpret_cast<char *>(&v), sizeof(v));
                // Insert current node into tt2
                tt2.insert(*it);
#ifdef DEBUGC
                cout << "FSTPT " << v << "\n";
                code_file.write(&FOK, sizeof(FOK));
                code_file.write(reinterpret_cast<char *>(&i), sizeof(i));
#endif
                i++;
              }
            for(temporary_terms_type::const_iterator it = ModelBlock->Block_List[j].Temporary_terms->begin();
                 it != ModelBlock->Block_List[j].Temporary_terms->end(); it++)
              {
                map_idx_type::const_iterator ii=map_idx.find((*it)->idx);
#ifdef DEBUGC
                cout << "map_idx[" << (*it)->idx <<"]=" << ii->second << "\n";
#endif
              }
            // The equations
            for(i = 0;i < ModelBlock->Block_List[j].Size;i++)
              {
                ModelBlock->Block_List[j].Variable_Sorted[i] = variable_table.getID(eEndogenous, ModelBlock->Block_List[j].Variable[i], 0);
                if (!lhs_rhs_done)
                  {
                    eq_node = equations[ModelBlock->Block_List[j].Equation[i]];
                    lhs = eq_node->arg1;
                    rhs = eq_node->arg2;
                  }
                switch (ModelBlock->Block_List[j].Simulation_Type)
                  {
                    case EVALUATE_BACKWARD:
                    case EVALUATE_FORWARD:
                      rhs->compile(code_file,false, output_type, temporary_terms, map_idx);
                      lhs->compile(code_file,true, output_type, temporary_terms, map_idx);
                      break;
                    case EVALUATE_BACKWARD_R:
                    case EVALUATE_FORWARD_R:
                      lhs->compile(code_file,false, output_type, temporary_terms, map_idx);
                      rhs->compile(code_file,true, output_type, temporary_terms, map_idx);
                      break;
                    case SOLVE_TWO_BOUNDARIES_SIMPLE:
                      v=ModelBlock->Block_List[j].Equation[i];
                      Uf[v].eqr=i;
                      Uf[v].Ufl=NULL;
                      goto end;
                    case SOLVE_BACKWARD_COMPLETE:
                    case SOLVE_FORWARD_COMPLETE:
                      v=ModelBlock->Block_List[j].Equation[i];
                      Uf[v].eqr=i;
                      Uf[v].Ufl=NULL;
                      goto end;
                    case SOLVE_TWO_BOUNDARIES_COMPLETE:
                      v=ModelBlock->Block_List[j].Equation[i];
                      Uf[v].eqr=i;
                      Uf[v].Ufl=NULL;
                      goto end;
                    default:
                      end:
                      lhs->compile(code_file,false, output_type, temporary_terms, map_idx);
                      rhs->compile(code_file,false, output_type, temporary_terms, map_idx);
                      code_file.write(&FBINARY, sizeof(FBINARY));
                      int v=oMinus;
                      code_file.write(reinterpret_cast<char *>(&v),sizeof(v));
                      code_file.write(&FSTPR, sizeof(FSTPR));
                      code_file.write(reinterpret_cast<char *>(&i), sizeof(i));
#ifdef CONDITION
                      if (ModelBlock->Block_List[j].Simulation_Type==SOLVE_TWO_BOUNDARIES_COMPLETE)
                        output << "  condition[" << i << "]=0;\n";
#endif
                  }
              }
            code_file.write(&FENDEQU, sizeof(FENDEQU));
            // The Jacobian if we have to solve the block
            if (ModelBlock->Block_List[j].Simulation_Type!=EVALUATE_BACKWARD
                && ModelBlock->Block_List[j].Simulation_Type!=EVALUATE_FORWARD
                && ModelBlock->Block_List[j].Simulation_Type!=EVALUATE_BACKWARD_R
                && ModelBlock->Block_List[j].Simulation_Type!=EVALUATE_FORWARD_R)
              {
                switch (ModelBlock->Block_List[j].Simulation_Type)
                  {
                    case SOLVE_BACKWARD_SIMPLE:
                    case SOLVE_FORWARD_SIMPLE:
                      compileDerivative(code_file, ModelBlock->Block_List[j].Equation[0], ModelBlock->Block_List[j].Variable[0], 0, output_type, map_idx);
                      code_file.write(&FSTPG, sizeof(FSTPG));
                      v=0;
                      code_file.write(reinterpret_cast<char *>(&v), sizeof(v));
                      break;
                    case SOLVE_BACKWARD_COMPLETE:
                    case SOLVE_FORWARD_COMPLETE:
                      m=ModelBlock->Block_List[j].Max_Lag;
                      for(i=0;i<ModelBlock->Block_List[j].IM_lead_lag[m].size;i++)
                        {
                          int eq=ModelBlock->Block_List[j].IM_lead_lag[m].Equ_Index[i];
                          int var=ModelBlock->Block_List[j].IM_lead_lag[m].Var_Index[i];
                          int u=ModelBlock->Block_List[j].IM_lead_lag[m].us[i];
                          int eqr=ModelBlock->Block_List[j].IM_lead_lag[m].Equ[i];
                          int v=ModelBlock->Block_List[j].Equation[eqr];
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
                          Uf[v].Ufl->u=u;
                          Uf[v].Ufl->var=var;
                          compileDerivative(code_file, eq, var, 0, output_type, map_idx);
                          code_file.write(&FSTPU, sizeof(FSTPU));
                          code_file.write(reinterpret_cast<char *>(&u), sizeof(u));
                        }
                      for(i = 0;i < ModelBlock->Block_List[j].Size;i++)
                        {
                          code_file.write(&FLDR, sizeof(FLDR));
                          code_file.write(reinterpret_cast<char *>(&i), sizeof(i));
                          code_file.write(&FLDZ, sizeof(FLDZ));
                          int v=ModelBlock->Block_List[j].Equation[i];
                          for(Uf[v].Ufl=Uf[v].Ufl_First;Uf[v].Ufl;Uf[v].Ufl=Uf[v].Ufl->pNext)
                            {
                              code_file.write(&FLDU, sizeof(FLDU));
                              code_file.write(reinterpret_cast<char *>(&Uf[v].Ufl->u), sizeof(Uf[v].Ufl->u));
                              code_file.write(&FLDV, sizeof(FLDV));
                              char vc=eEndogenous;
                              code_file.write(reinterpret_cast<char *>(&vc), sizeof(vc));
                              code_file.write(reinterpret_cast<char *>(&Uf[v].Ufl->var), sizeof(Uf[v].Ufl->var));
                              int v1=0;
                              code_file.write(reinterpret_cast<char *>(&v1), sizeof(v1));
                              code_file.write(&FBINARY, sizeof(FBINARY));
                              v1=oTimes;
                              code_file.write(reinterpret_cast<char *>(&v1), sizeof(v1));
                              code_file.write(&FCUML, sizeof(FCUML));
                            }
                          Uf[v].Ufl=Uf[v].Ufl_First;
                          while(Uf[v].Ufl)
                            {
                              Uf[v].Ufl_First=Uf[v].Ufl->pNext;
                              free(Uf[v].Ufl);
                              Uf[v].Ufl=Uf[v].Ufl_First;
                            }
                          code_file.write(&FBINARY, sizeof(FBINARY));
                          v=oMinus;
                          code_file.write(reinterpret_cast<char *>(&v), sizeof(v));
                          code_file.write(&FSTPU, sizeof(FSTPU));
                          code_file.write(reinterpret_cast<char *>(&i), sizeof(i));
                        }
                      break;
                    case SOLVE_TWO_BOUNDARIES_COMPLETE:
                    case SOLVE_TWO_BOUNDARIES_SIMPLE:
                      for(m=0;m<=ModelBlock->Block_List[j].Max_Lead+ModelBlock->Block_List[j].Max_Lag;m++)
                        {
                          k=m-ModelBlock->Block_List[j].Max_Lag;
                          for(i=0;i<ModelBlock->Block_List[j].IM_lead_lag[m].size;i++)
                            {
                              int eq=ModelBlock->Block_List[j].IM_lead_lag[m].Equ_Index[i];
                              int var=ModelBlock->Block_List[j].IM_lead_lag[m].Var_Index[i];
                              int u=ModelBlock->Block_List[j].IM_lead_lag[m].u[i];
                              int eqr=ModelBlock->Block_List[j].IM_lead_lag[m].Equ[i];
                              int v=ModelBlock->Block_List[j].Equation[eqr];
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
                              Uf[v].Ufl->u=u;
                              Uf[v].Ufl->var=var;
                              Uf[v].Ufl->lag=k;
                              compileDerivative(code_file, eq, var, k, output_type, map_idx);
                              code_file.write(&FSTPU, sizeof(FSTPU));
                              code_file.write(reinterpret_cast<char *>(&u), sizeof(u));
#ifdef CONDITION
                              output << "  if (fabs(condition[" << eqr << "])<fabs(u[" << u << "+Per_u_]))\n";
                              output << "    condition[" << eqr << "]=u[" << u << "+Per_u_];\n";
#endif
                            }
                        }
                      for(i = 0;i < ModelBlock->Block_List[j].Size;i++)
                        {
                          code_file.write(&FLDR, sizeof(FLDR));
                          code_file.write(reinterpret_cast<char *>(&i), sizeof(i));
                          code_file.write(&FLDZ, sizeof(FLDZ));
                          int v=ModelBlock->Block_List[j].Equation[i];
                          for(Uf[v].Ufl=Uf[v].Ufl_First;Uf[v].Ufl;Uf[v].Ufl=Uf[v].Ufl->pNext)
                            {
                              code_file.write(&FLDU, sizeof(FLDU));
                              code_file.write(reinterpret_cast<char *>(&Uf[v].Ufl->u), sizeof(Uf[v].Ufl->u));
                              code_file.write(&FLDV, sizeof(FLDV));
                              char vc=eEndogenous;
                              code_file.write(reinterpret_cast<char *>(&vc), sizeof(vc));
                              int v1=Uf[v].Ufl->var;
                              code_file.write(reinterpret_cast<char *>(&v1), sizeof(v1));
                              v1=Uf[v].Ufl->lag;
                              code_file.write(reinterpret_cast<char *>(&v1), sizeof(v1));
                              code_file.write(&FBINARY, sizeof(FBINARY));
                              v1=oTimes;
                              code_file.write(reinterpret_cast<char *>(&v1), sizeof(v1));
                              code_file.write(&FCUML, sizeof(FCUML));
                            }
                          Uf[v].Ufl=Uf[v].Ufl_First;
                          while(Uf[v].Ufl)
                            {
                              Uf[v].Ufl_First=Uf[v].Ufl->pNext;
                              free(Uf[v].Ufl);
                              Uf[v].Ufl=Uf[v].Ufl_First;
                            }
                          code_file.write(&FBINARY, sizeof(FBINARY));
                          v=oMinus;
                          code_file.write(reinterpret_cast<char *>(&v), sizeof(v));
                          code_file.write(&FSTPU, sizeof(FSTPU));
                          code_file.write(reinterpret_cast<char *>(&i), sizeof(i));
#ifdef CONDITION
                          output << "  if (fabs(condition[" << i << "])<fabs(u[" << i << "+Per_u_]))\n";
                          output << "    condition[" << i << "]=u[" << i << "+Per_u_];\n";
#endif
                        }
#ifdef CONDITION
                      for(m=0;m<=ModelBlock->Block_List[j].Max_Lead+ModelBlock->Block_List[j].Max_Lag;m++)
                        {
                          k=m-ModelBlock->Block_List[j].Max_Lag;
                          for(i=0;i<ModelBlock->Block_List[j].IM_lead_lag[m].size;i++)
                            {
                              int eq=ModelBlock->Block_List[j].IM_lead_lag[m].Equ_Index[i];
                              int var=ModelBlock->Block_List[j].IM_lead_lag[m].Var_Index[i];
                              int u=ModelBlock->Block_List[j].IM_lead_lag[m].u[i];
                              int eqr=ModelBlock->Block_List[j].IM_lead_lag[m].Equ[i];
                              output << "  u[" << u << "+Per_u_] /= condition[" << eqr << "];\n";
                            }
                        }
                      for(i = 0;i < ModelBlock->Block_List[j].Size;i++)
                        output << "  u[" << i << "+Per_u_] /= condition[" << i << "];\n";
#endif
                      break;
                  }

                prev_Simulation_Type=ModelBlock->Block_List[j].Simulation_Type;
              }
            j++;
          }
      }
    code_file.write(&FENDBLOCK, sizeof(FENDBLOCK));
    code_file.write(&FEND, sizeof(FEND));
    code_file.close();
  }


void
ModelTree::writeStaticMFile(const string &static_basename) const
{
  string filename = static_basename + ".m";

  ofstream mStaticModelFile;
  mStaticModelFile.open(filename.c_str(), ios::out | ios::binary);
  if (!mStaticModelFile.is_open())
    {
      cerr << "Error: Can't open file " << filename << " for writing" << endl;
      exit(-1);
    }
  // Writing comments and function definition command
  mStaticModelFile << "function [residual, g1, g2] = " << static_basename << "(y, x, params)" << endl
                   << "%" << endl
                   << "% Status : Computes static model for Dynare" << endl
                   << "%" << endl
                   << "% Warning : this file is generated automatically by Dynare" << endl
                   << "%           from model file (.mod)" << endl << endl;

  writeStaticModel(mStaticModelFile);

  mStaticModelFile.close();
}


void
ModelTree::writeDynamicMFile(const string &dynamic_basename) const
{
  string filename = dynamic_basename + ".m";

  ofstream mDynamicModelFile;
  mDynamicModelFile.open(filename.c_str(), ios::out | ios::binary);
  if (!mDynamicModelFile.is_open())
    {
      cerr << "Error: Can't open file " << filename << " for writing" << endl;
      exit(-1);
    }
  mDynamicModelFile << "function [residual, g1, g2, g3] = " << dynamic_basename << "(y, x, params, it_)" << endl
                    << "%" << endl
                    << "% Status : Computes dynamic model for Dynare" << endl
                    << "%" << endl
                    << "% Warning : this file is generated automatically by Dynare" << endl
                    << "%           from model file (.mod)" << endl << endl;

  writeDynamicModel(mDynamicModelFile);

  mDynamicModelFile.close();
}

void
ModelTree::writeStaticCFile(const string &static_basename) const
{
  string filename = static_basename + ".c";

  ofstream mStaticModelFile;
  mStaticModelFile.open(filename.c_str(), ios::out | ios::binary);
  if (!mStaticModelFile.is_open())
    {
      cerr << "Error: Can't open file " << filename << " for writing" << endl;
      exit(-1);
    }
  mStaticModelFile << "/*" << endl
                   << " * " << filename << " : Computes static model for Dynare" << endl
                   << " * Warning : this file is generated automatically by Dynare" << endl
                   << " *           from model file (.mod)" << endl
                   << endl
                   << " */" << endl
                   << "#include <math.h>" << endl
                   << "#include \"mex.h\"" << endl;

  // Writing the function Static
  writeStaticModel(mStaticModelFile);

  // Writing the gateway routine
  mStaticModelFile << "/* The gateway routine */" << endl
                   << "void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])" << endl
                   << "{" << endl
                   << "  double *y, *x, *params;" << endl
                   << "  double *residual, *g1;" << endl
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
                   << "  residual = NULL;" << endl
                   << "  if (nlhs >= 1)" << endl
                   << "  {" << endl
                   << "      /* Set the output pointer to the output matrix residual. */" << endl
                   << "      plhs[0] = mxCreateDoubleMatrix(" << equations.size() << ",1, mxREAL);" << endl
                   << "     /* Create a C pointer to a copy of the output matrix residual. */" << endl
                   << "     residual = mxGetPr(plhs[0]);" << endl
                   << "  }" << endl
                   << endl
                   << "  g1 = NULL;" << endl
                   << "  if (nlhs >= 2)" << endl
                   << "  {" << endl
                   << "      /* Set the output pointer to the output matrix g1. */" << endl
                   << "      plhs[1] = mxCreateDoubleMatrix(" << equations.size() << ", " << symbol_table.endo_nbr << ", mxREAL);" << endl
                   << "      /* Create a C pointer to a copy of the output matrix g1. */" << endl
                   << "      g1 = mxGetPr(plhs[1]);" << endl
                   << "  }" << endl
                   << endl
                   << "  /* Call the C Static. */" << endl
                   << "  Static(y, x, params, residual, g1);" << endl
                   << "}" << endl;

  mStaticModelFile.close();
}

void
ModelTree::writeDynamicCFile(const string &dynamic_basename) const
{
  string filename = dynamic_basename + ".c";
  ofstream mDynamicModelFile;

  mDynamicModelFile.open(filename.c_str(), ios::out | ios::binary);
  if (!mDynamicModelFile.is_open())
    {
      cerr << "Error: Can't open file " << filename << " for writing" << endl;
      exit(-1);
    }
  mDynamicModelFile << "/*" << endl
                    << " * " << filename << " : Computes dynamic model for Dynare" << endl
                    << " *" << endl
                    << " * Warning : this file is generated automatically by Dynare" << endl
                    << " *           from model file (.mod)" << endl
                    << endl
                    << " */" << endl
                    << "#include <math.h>" << endl
                    << "#include \"mex.h\"" << endl;

  // Writing the function body
  writeDynamicModel(mDynamicModelFile);

  // Writing the gateway routine
  mDynamicModelFile << "/* The gateway routine */" << endl
                    << "void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])" << endl
                    << "{" << endl
                    << "  double *y, *x, *params;" << endl
                    << "  double *residual, *g1, *g2;" << endl
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

                    << "     plhs[1] = mxCreateDoubleMatrix(" << equations.size() << ", " << variable_table.getDynJacobianColsNbr(computeJacobianExo) << ", mxREAL);" << endl
                    << "     /* Create a C pointer to a copy of the output matrix g1. */" << endl
                    << "     g1 = mxGetPr(plhs[1]);" << endl
                    << "  }" << endl
                    << endl
                    << "  g2 = NULL;" << endl
                    << " if (nlhs >= 3)" << endl
                    << "  {" << endl
                    << "     /* Set the output pointer to the output matrix g2. */" << endl
                    << "     plhs[2] = mxCreateDoubleMatrix(" << equations.size() << ", " << variable_table.getDynJacobianColsNbr(computeJacobianExo)*variable_table.getDynJacobianColsNbr(computeJacobianExo) << ", mxREAL);" << endl
                    << "     /* Create a C pointer to a copy of the output matrix g1. */" << endl
                    << "     g2 = mxGetPr(plhs[2]);" << endl
                    << "  }" << endl
                    << endl
                    << "  /* Call the C subroutines. */" << endl
                    << "  Dynamic(y, x, nb_row_x, params, it_, residual, g1, g2);" << endl
                    << "}" << endl;
  mDynamicModelFile.close();
}

void
ModelTree::writeStaticModel(ostream &StaticOutput) const
{
  ostringstream model_output;    // Used for storing model equations
  ostringstream jacobian_output; // Used for storing jacobian equations
  ostringstream hessian_output;
  ostringstream lsymetric;       // For symmetric elements in hessian

  ExprNodeOutputType output_type = (mode == eDLLMode ? oCStaticModel : oMatlabStaticModel);

  writeModelLocalVariables(model_output, output_type);

  writeTemporaryTerms(model_output, output_type);

  writeModelEquations(model_output, output_type);

  // Write Jacobian w.r. to endogenous only
  for(first_derivatives_type::const_iterator it = first_derivatives.begin();
      it != first_derivatives.end(); it++)
    {
      int eq = it->first.first;
      int var = it->first.second;
      NodeID d1 = it->second;

      if (variable_table.getType(var) == eEndogenous)
        {
          ostringstream g1;
          g1 << "  g1";
          matrixHelper(g1, eq, variable_table.getSymbolID(var), output_type);

          jacobian_output << g1.str() << "=" << g1.str() << "+";
          d1->writeOutput(jacobian_output, output_type, temporary_terms);
          jacobian_output << ";" << endl;
        }
    }

  // Write Hessian w.r. to endogenous only
  if (computeStaticHessian)
    for(second_derivatives_type::const_iterator it = second_derivatives.begin();
        it != second_derivatives.end(); it++)
      {
        int eq = it->first.first;
        int var1 = it->first.second.first;
        int var2 = it->first.second.second;
        NodeID d2 = it->second;

        // Keep only derivatives w.r. to endogenous variables
        if (variable_table.getType(var1) == eEndogenous
            && variable_table.getType(var2) == eEndogenous)
          {
            int id1 = variable_table.getSymbolID(var1);
            int id2 = variable_table.getSymbolID(var2);

            int col_nb = id1*symbol_table.endo_nbr+id2;
            int col_nb_sym = id2*symbol_table.endo_nbr+id1;

            hessian_output << "  g2";
            matrixHelper(hessian_output, eq, col_nb, output_type);
            hessian_output << " = ";
            d2->writeOutput(hessian_output, output_type, temporary_terms);
            hessian_output << ";" << endl;

            // Treating symetric elements
            if (var1 != var2)
              {
                lsymetric <<  "  g2";
                matrixHelper(lsymetric, eq, col_nb_sym, output_type);
                lsymetric << " = " <<  "g2";
                matrixHelper(lsymetric, eq, col_nb, output_type);
                lsymetric << ";" << endl;
              }
          }
      }

  // Writing ouputs
  if (mode != eDLLMode)
    {
      StaticOutput << "residual = zeros( " << equations.size() << ", 1);" << endl << endl
                   << "%" << endl
                   << "% Model equations" << endl
                   << "%" << endl
                   << endl
                   << model_output.str()
                   << "if ~isreal(residual)" << endl
                   << "  residual = real(residual)+imag(residual).^2;" << endl
                   << "end" << endl
                   << "if nargout >= 2," << endl
                   << "  g1 = zeros(" << equations.size() << ", " << symbol_table.endo_nbr << ");" << endl
                   << endl
                   << "%" << endl
                   << "% Jacobian matrix" << endl
                   << "%" << endl
                   << endl
                   << jacobian_output.str()
                   << "  if ~isreal(g1)" << endl
                   << "    g1 = real(g1)+2*imag(g1);" << endl
                   << "  end" << endl
                   << "end" << endl;
      if (computeStaticHessian)
        {
          StaticOutput << "if nargout >= 3,\n";
          // Writing initialization instruction for matrix g2
          int ncols = symbol_table.endo_nbr * symbol_table.endo_nbr;
          StaticOutput << "  g2 = sparse([],[],[], " << equations.size() << ", " << ncols << ", " << 5*ncols << ");" << endl
                       << endl
                       << "%" << endl
                       << "% Hessian matrix" << endl
                       << "%" << endl
                       << endl
                       << hessian_output.str()
                       << lsymetric.str()
                       << "end;" << endl;
        }
    }
  else
    {
      StaticOutput << "void Static(double *y, double *x, double *params, double *residual, double *g1)" << endl
                   << "{" << endl
                   << "  double lhs, rhs;" << endl
        // Writing residual equations
                   << "  /* Residual equations */" << endl
                   << "  if (residual == NULL)" << endl
                   << "    return;" << endl
                   << "  else" << endl
                   << "    {" << endl
                   << model_output.str()
        // Writing Jacobian
                   << "     /* Jacobian for endogenous variables without lag */" << endl
                   << "     if (g1 == NULL)" << endl
                   << "       return;" << endl
                   << "     else" << endl
                   << "       {" << endl
                   << jacobian_output.str()
                   << "       }" << endl
                   << "    }" << endl
                   << "}" << endl << endl;
    }
}

string
ModelTree::reform(const string name1) const
{
  string name=name1;
  int pos = name.find("\\", 0);
  while(pos >= 0)
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
ModelTree::Write_Inf_To_Bin_File(const string &dynamic_basename, const string &bin_basename, const int &num,
                                 int &u_count_int, bool &file_open) const
{
  int j;
  std::ofstream SaveCode;
  if (file_open)
    SaveCode.open((bin_basename + ".bin").c_str(), ios::out | ios::in | ios::binary | ios ::ate );
  else
    SaveCode.open((bin_basename + ".bin").c_str(), ios::out | ios::binary);
  if (!SaveCode.is_open())
    {
      cout << "Error : Can't open file \"" << bin_basename << ".bin\" for writing\n";
      exit( -1);
    }
  u_count_int=0;
  for(int m=0;m<=block_triangular.ModelBlock->Block_List[num].Max_Lead+block_triangular.ModelBlock->Block_List[num].Max_Lag;m++)
    {
      int k1=m-block_triangular.ModelBlock->Block_List[num].Max_Lag;
      for(j=0;j<block_triangular.ModelBlock->Block_List[num].IM_lead_lag[m].size;j++)
        {
          int varr=block_triangular.ModelBlock->Block_List[num].IM_lead_lag[m].Var[j]+k1*block_triangular.ModelBlock->Block_List[num].Size;
          int u=block_triangular.ModelBlock->Block_List[num].IM_lead_lag[m].u[j];
          int eqr1=block_triangular.ModelBlock->Block_List[num].IM_lead_lag[m].Equ[j];
          /*cout << "   ! IM_i[std::make_pair(std::make_pair(" << eqr1 << ", " << varr+k1*block_triangular.ModelBlock->Block_List[num].Size << "), " << k1 << ")] = " << u << ";\n";
          cout << "   ? IM_i[std::make_pair(std::make_pair(" << eqr1 << ", " << varr << "), " << k1 << ")] = " << u << ";\n";*/
          SaveCode.write(reinterpret_cast<char *>(&eqr1), sizeof(eqr1));
          SaveCode.write(reinterpret_cast<char *>(&varr), sizeof(varr));
          SaveCode.write(reinterpret_cast<char *>(&k1), sizeof(k1));
          SaveCode.write(reinterpret_cast<char *>(&u), sizeof(u));
          //cout << "eqr1=" << eqr1 << " varr=" << varr << " k1=" << k1 << " u=" << u << "\n";
          u_count_int++;
        }
    }
  for(j=0;j<block_triangular.ModelBlock->Block_List[num].Size;j++)
    {
       int eqr1=j;
       int varr=block_triangular.ModelBlock->Block_List[num].Size*(block_triangular.periods
                                                                   +block_triangular.Model_Max_Lead);
       int k1=0;
       SaveCode.write(reinterpret_cast<char *>(&eqr1), sizeof(eqr1));
       SaveCode.write(reinterpret_cast<char *>(&varr), sizeof(varr));
       SaveCode.write(reinterpret_cast<char *>(&k1), sizeof(k1));
       SaveCode.write(reinterpret_cast<char *>(&eqr1), sizeof(eqr1));
       //cout << "eqr1=" << eqr1 << " varr=" << varr << " k1=" << k1 << " eqr1=" << eqr1 << "\n";
       u_count_int++;
    }
  for(j=0;j<block_triangular.ModelBlock->Block_List[num].Size;j++)
    {
      int varr=block_triangular.ModelBlock->Block_List[num].Variable[j];
      SaveCode.write(reinterpret_cast<char *>(&varr), sizeof(varr));
    }
  for(j=0;j<block_triangular.ModelBlock->Block_List[num].Size;j++)
    {
      int eqr1=block_triangular.ModelBlock->Block_List[num].Equation[j];
      SaveCode.write(reinterpret_cast<char *>(&eqr1), sizeof(eqr1));
    }
  SaveCode.close();
}

void
ModelTree::writeSparseStaticMFile(const string &static_basename, const string &bin_basename, const int mode) const
{
  string filename;
  ofstream mStaticModelFile;
  int i, k, prev_Simulation_Type;
  bool /*printed = false,*/ skip_head, open_par=false;
  filename = static_basename + ".m";
  mStaticModelFile.open(filename.c_str(), ios::out | ios::binary);
  if (!mStaticModelFile.is_open())
    {
      cerr << "Error: Can't open file " << filename << " for writing" << endl;
      exit(-1);
    }
  mStaticModelFile << "%\n";
  mStaticModelFile << "% " << filename << " : Computes static model for Dynare\n";
  mStaticModelFile << "%\n";
  mStaticModelFile << "% Warning : this file is generated automatically by Dynare\n";
  mStaticModelFile << "%           from model file (.mod)\n\n";
  mStaticModelFile << "%/\n";
  mStaticModelFile << "function [varargout] = " << static_basename << "(varargin)\n";
  mStaticModelFile << "  global oo_ options_ M_ ys0_ ;\n";
  mStaticModelFile << "  y_kmin=M_.maximum_lag;\n";
  mStaticModelFile << "  y_kmax=M_.maximum_lead;\n";
  mStaticModelFile << "  y_size=M_.endo_nbr;\n";
  mStaticModelFile << "  if(length(varargin)>0)\n";
  mStaticModelFile << "    %it is a simple evaluation of the dynamic model for time _it\n";
  mStaticModelFile << "    global it_;\n";
  mStaticModelFile << "    y=varargin{1}(:);\n";
  mStaticModelFile << "    ys=y;\n";
  mStaticModelFile << "    g1=[];\n";
  mStaticModelFile << "    x=varargin{2}(:);\n";
  mStaticModelFile << "    residual=zeros(1, " << symbol_table.endo_nbr << ");\n";
  prev_Simulation_Type=-1;
  for(i=0;i<block_triangular.ModelBlock->Size;i++)
    {
      mStaticModelFile << "    y_index=[";
      for(int ik=0;ik<block_triangular.ModelBlock->Block_List[i].Size;ik++)
        {
          mStaticModelFile << " " << block_triangular.ModelBlock->Block_List[i].Variable[ik]+1;
        }
      mStaticModelFile << " ];\n";
      k=block_triangular.ModelBlock->Block_List[i].Simulation_Type;
      if (BlockTriangular::BlockSim(prev_Simulation_Type)==BlockTriangular::BlockSim(k) &&
         (k==EVALUATE_FORWARD || k==EVALUATE_BACKWARD || k==EVALUATE_FORWARD_R || k==EVALUATE_BACKWARD_R))
        skip_head=true;
      else
        skip_head=false;
      switch(k)
        {
           case EVALUATE_FORWARD:
           case EVALUATE_BACKWARD:
           case EVALUATE_FORWARD_R:
           case EVALUATE_BACKWARD_R:
             if(!skip_head)
               mStaticModelFile << "    y=" << static_basename << "_" << i + 1 << "(y, x);\n";
             mStaticModelFile << "    residual(y_index)=ys(y_index)-y(y_index);\n";
             break;
           case SOLVE_FORWARD_COMPLETE:
           case SOLVE_BACKWARD_COMPLETE:
           case SOLVE_FORWARD_SIMPLE:
           case SOLVE_BACKWARD_SIMPLE:
           case SOLVE_TWO_BOUNDARIES_COMPLETE:
             mStaticModelFile << "    [r, g1]=" << static_basename << "_" <<  i + 1 << "(y, x);\n";
             mStaticModelFile << "    residual(y_index)=r;\n";
             break;
        }
      prev_Simulation_Type=k;
    }
  mStaticModelFile << "    varargout{1}=residual;\n";
  mStaticModelFile << "    varargout{2}=g1;\n";
  mStaticModelFile << "    return;\n";
  mStaticModelFile << "  end;\n";
  mStaticModelFile << "  %it is the deterministic simulation of the block decomposed static model\n";
  mStaticModelFile << "  periods=options_.periods;\n";
  mStaticModelFile << "  maxit_=options_.maxit_;\n";
  mStaticModelFile << "  solve_tolf=options_.solve_tolf;\n";
  mStaticModelFile << "  y=oo_.steady_state;\n";
  mStaticModelFile << "  x=oo_.exo_steady_state;\n";
  mStaticModelFile << "  varargout{2}=1;\n";
  prev_Simulation_Type=-1;
  for(i = 0;i < block_triangular.ModelBlock->Size;i++)
    {
      k = block_triangular.ModelBlock->Block_List[i].Simulation_Type;
      if (BlockTriangular::BlockSim(prev_Simulation_Type)==BlockTriangular::BlockSim(k) &&
         (k==EVALUATE_FORWARD || k==EVALUATE_BACKWARD || k==EVALUATE_FORWARD_R || k==EVALUATE_BACKWARD_R))
        skip_head=true;
      else
        skip_head=false;
      if ((k == EVALUATE_FORWARD || k == EVALUATE_FORWARD_R || k == EVALUATE_BACKWARD || k == EVALUATE_BACKWARD_R) && (block_triangular.ModelBlock->Block_List[i].Size))
        {
          if (!skip_head)
            {
              if (open_par)
                {
                  mStaticModelFile << "  end\n";
                }
             mStaticModelFile << "  y=" << static_basename << "_" << i + 1 << "(y, x);\n";
            }
          open_par=false;
        }
       else if ((k == SOLVE_FORWARD_SIMPLE || k == SOLVE_BACKWARD_SIMPLE) && (block_triangular.ModelBlock->Block_List[i].Size))
        {
          if (open_par)
            {
               mStaticModelFile << "  end\n";
            }
          open_par=false;
          mStaticModelFile << "  g1=0;\n";
          mStaticModelFile << "  r=0;\n";
          /*mStaticModelFile << "    for it_=y_kmin+1:periods+y_kmin\n";
          mStaticModelFile << "        cvg=0;\n";
          mStaticModelFile << "        iter=0;\n";
          mStaticModelFile << "        Per_y_=it_*y_size;\n";
          mStaticModelFile << "        while ~(cvg==1 | iter>maxit_),\n";
          mStaticModelFile << "            [r, g1] = " << static_basename << "_" << i + 1 << "(y, x, it_, Per_y_, y_size);\n";
          mStaticModelFile << "            y[it_, " << block_triangular.ModelBlock->Block_List[i].Variable[0] << "] = y[it_, " << block_triangular.ModelBlock->Block_List[i].Variable[0] << "]-r/g1;\n";
          mStaticModelFile << "            cvg=((r[0]*r[0])<solve_tolf);\n";
          mStaticModelFile << "            iter=iter+1;\n";
          mStaticModelFile << "        end\n";
          mStaticModelFile << "        if cvg==0\n";
          mStaticModelFile << "            fprintf('Convergence not achieved in block " << i << ", at time %d after %d iterations\\n',it_,iter);\n";
          mStaticModelFile << "            return;\n";
          mStaticModelFile << "        end\n";
          mStaticModelFile << "    end\n";*/
          mStaticModelFile << "  cvg=0;\n";
          mStaticModelFile << "  iter=0;\n";
          /*mStaticModelFile << "    Per_y_=it_*y_size;\n";*/
          mStaticModelFile << "  while ~(cvg==1 | iter>maxit_),\n";
          mStaticModelFile << "    [r, g1] = " << static_basename << "_" << i + 1 << "(y, x);\n";
          mStaticModelFile << "    y(" << block_triangular.ModelBlock->Block_List[i].Variable[0]+1 << ") = y(" << block_triangular.ModelBlock->Block_List[i].Variable[0]+1 << ")-r/g1;\n";
          mStaticModelFile << "    cvg=((r*r)<solve_tolf);\n";
          mStaticModelFile << "    iter=iter+1;\n";
          mStaticModelFile << "  end\n";
          mStaticModelFile << "  if cvg==0\n";
          mStaticModelFile << "     fprintf('Convergence not achieved in block " << i << ", after %d iterations\\n',iter);\n";
          mStaticModelFile << "     return;\n";
          mStaticModelFile << "  end\n";
        }
      else if ((k == SOLVE_FORWARD_COMPLETE || k == SOLVE_BACKWARD_COMPLETE || k == SOLVE_TWO_BOUNDARIES_COMPLETE) && (block_triangular.ModelBlock->Block_List[i].Size))
        {
          if (open_par)
            {
              mStaticModelFile << "end\n";
            }
          open_par=false;
          mStaticModelFile << "  y_index=[";
          for(int ik=0;ik<block_triangular.ModelBlock->Block_List[i].Size;ik++)
            {
              mStaticModelFile << " " << block_triangular.ModelBlock->Block_List[i].Variable[ik]+1;
            }
          mStaticModelFile << " ];\n";
          mStaticModelFile << "  g1=0;g2=0;g3=0;\n";
          mStaticModelFile << "  r=0;\n";
          /*mStaticModelFile << "  for it_=y_kmin+1:periods+y_kmin\n";
          mStaticModelFile << "      cvg=0;\n";
          mStaticModelFile << "      iter=0;\n";
          mStaticModelFile << "      Per_y_=it_*y_size;\n";
          mStaticModelFile << "      while ~(cvg==1 | iter>maxit_),\n";
          mStaticModelFile << "          [r, g1, g2, g3, b] = " << static_basename << "_" << i + 1 << "(y, x, it_, Per_y_, y_size);\n";
          mStaticModelFile << "          [L, U] = LU(g1);\n";
          mStaticModelFile << "          y(it_, y_index) = U\\(L\\b);\n";
          mStaticModelFile << "          cvg=((r'*r)<solve_tolf);\n";
          mStaticModelFile << "          iter=iter+1;\n";
          mStaticModelFile << "      end\n";
          mStaticModelFile << "      if cvg==0\n";
          mStaticModelFile << "          fprintf('Convergence not achieved in block " << i << ", at time %d after %d iterations\\n',it_,iter);\n";
          mStaticModelFile << "          return;\n";
          mStaticModelFile << "      end\n";
          mStaticModelFile << "  end\n";*/
          mStaticModelFile << "  cvg=0;\n";
          mStaticModelFile << "  iter=0;\n";
          /*mStaticModelFile << "  Per_y_=it_*y_size;\n";*/
          mStaticModelFile << "  lambda=1;\n";
          mStaticModelFile << "  stpmx = 100 ;\n";
          mStaticModelFile << "  stpmax = stpmx*max([sqrt(y'*y);size(y_index,2)]);\n";
          mStaticModelFile << "  nn=1:size(y_index,2);\n";
          mStaticModelFile << "  while ~(cvg==1 | iter>maxit_),\n";
          mStaticModelFile << "    [r, g1, g2, g3, b] = " << static_basename << "_" << i + 1 << "(y, x);\n";
          mStaticModelFile << "    max_res=max(abs(r));\n";
          mStaticModelFile << "    cvg=(max_res<solve_tolf);\n";
          mStaticModelFile << "    if (cvg==0),\n";
          mStaticModelFile << "      g = (r'*g1)';\n";
          mStaticModelFile << "      f = 0.5*r'*r;\n";
          mStaticModelFile << "      p = -g1\\r ;\n";
          mStaticModelFile << "      [y,f,r,check]=lnsrch1(y,f,g,p,stpmax,@" << static_basename << "_" << i + 1 << ",nn,y_index,x);\n";
          mStaticModelFile << "    end;\n";
          mStaticModelFile << "    iter=iter+1;\n";
          mStaticModelFile << "    disp(['iter=' num2str(iter,'%d') ' err=' num2str(max_res,'%f')]);\n";
          mStaticModelFile << "  end\n";
          mStaticModelFile << "  if cvg==0\n";
          mStaticModelFile << "    fprintf('Error in steady: Convergence not achieved in block " << i << ", after %d iterations\\n',iter);\n";
          mStaticModelFile << "    return;\n";
          mStaticModelFile << "  else\n";
          mStaticModelFile << "    fprintf('convergence achieved after %d iterations\\n',iter);\n";
          mStaticModelFile << "  end\n";
        }
      prev_Simulation_Type=k;
    }
  if(open_par)
    mStaticModelFile << "  end;\n";
  mStaticModelFile << "  oo_.steady_state = y;\n";
  /*mStaticModelFile << "  if isempty(ys0_)\n";
  mStaticModelFile << "    oo_.endo_simul(:,1:M_.maximum_lag) = oo_.steady_state * ones(1,M_.maximum_lag);\n";
  mStaticModelFile << "  else\n";
  mStaticModelFile << "    options_ =set_default_option(options_,'periods',1);\n";
  mStaticModelFile << "    oo_.endo_simul(:,M_.maximum_lag+1:M_.maximum_lag+options_.periods+M_.maximum_lead) = oo_.steady_state * ones(1,options_.periods+M_.maximum_lead);\n";
  mStaticModelFile << "  end;\n";*/
  mStaticModelFile << "  disp('Steady State value');\n";
  mStaticModelFile << "  disp([strcat(M_.endo_names,' : ') num2str(oo_.steady_state,'%f')]);\n";
  mStaticModelFile << "  varargout{2}=0;\n";
  mStaticModelFile << "  varargout{1}=oo_.steady_state;\n";
  mStaticModelFile << "return;\n";
  writeModelStaticEquationsOrdered_M(mStaticModelFile, block_triangular.ModelBlock, static_basename);
  mStaticModelFile.close();
}

void
ModelTree::writeSparseDynamicMFile(const string &dynamic_basename, const string &basename, const int mode) const
{
  string sp;
  ofstream mDynamicModelFile;
  ostringstream tmp, tmp1, tmp_eq;
  int prev_Simulation_Type, tmp_i;
  SymbolicGaussElimination SGE;
  bool OK;
  string filename = dynamic_basename + ".m";
  mDynamicModelFile.open(filename.c_str(), ios::out | ios::binary);
  if (!mDynamicModelFile.is_open())
    {
      cerr << "Error: Can't open file " << filename << " for writing" << endl;
      exit(-1);
    }
  mDynamicModelFile << "%\n";
  mDynamicModelFile << "% " << filename << " : Computes dynamic model for Dynare\n";
  mDynamicModelFile << "%\n";
  mDynamicModelFile << "% Warning : this file is generated automatically by Dynare\n";
  mDynamicModelFile << "%           from model file (.mod)\n\n";
  mDynamicModelFile << "%/\n";

  int i, k, Nb_SGE=0;
  bool printed = false, skip_head, open_par=false;
  if (computeJacobian || computeJacobianExo || computeHessian)
    {
      mDynamicModelFile << "function [varargout] = " << dynamic_basename << "(varargin)\n";
      mDynamicModelFile << "  global oo_ options_ M_ ;\n";
      mDynamicModelFile << "  g2=[];g3=[];\n";
      //Temporary variables declaration
      OK=true;
      ostringstream tmp_output;
      for(temporary_terms_type::const_iterator it = temporary_terms.begin();
          it != temporary_terms.end(); it++)
        {
          if (OK)
            OK=false;
          else
            tmp_output << " ";
          (*it)->writeOutput(tmp_output, oMatlabDynamicModel, temporary_terms);
        }
      if (tmp_output.str().length()>0)
        mDynamicModelFile << "  global " << tmp_output.str() << " M_ ;\n";

      mDynamicModelFile << "  T_init=zeros(1,options_.periods+M_.maximum_lag+M_.maximum_lead);\n";
      tmp_output.str("");
      //OK=true;
      for(temporary_terms_type::const_iterator it = temporary_terms.begin();
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
      //mDynamicModelFile << "    global it_;\n";
      mDynamicModelFile << "    it_=varargin{3};\n";
      //mDynamicModelFile << "    g=zeros(y_size,y_size*(M_.maximum_endo_lag+M_.maximum_endo_lead+1));\n";
      mDynamicModelFile << "    g1=spalloc(y_size,y_size*(M_.maximum_endo_lag+M_.maximum_endo_lead+1),y_size*y_size*(M_.maximum_endo_lag+M_.maximum_endo_lead+1));\n";
      mDynamicModelFile << "    Per_u_=0;\n";
      mDynamicModelFile << "    Per_y_=it_*y_size;\n";
      mDynamicModelFile << "    y=varargin{1};\n";
      mDynamicModelFile << "    ys=y(it_,:);\n";
      /*mDynamicModelFile << "    y1=varargin{1};\n";
        mDynamicModelFile << "    cnb_nz_elem=1;\n";
        mDynamicModelFile << "    for i = -y_kmin:y_kmax\n";
        mDynamicModelFile << "      nz_elem=find(M_.lead_lag_incidence(:,1+i+y_kmin));\n";
        mDynamicModelFile << "      nb_nz_elem=length(nz_elem);\n";
        mDynamicModelFile << "      y(it_+i, nz_elem)=y1(cnb_nz_elem:(cnb_nz_elem+nb_nz_elem));\n";
        mDynamicModelFile << "      if(i==0)\n";
        mDynamicModelFile << "        ys(nz_elem)=y(it_, nz_elem);\n";
        mDynamicModelFile << "        nz_elem_s=nz_elem;\n";
        mDynamicModelFile << "      end;\n";
        mDynamicModelFile << "      cnb_nz_elem=cnb_nz_elem+nb_nz_elem;\n";
        mDynamicModelFile << "    end;\n";*/
      mDynamicModelFile << "    x=varargin{2};\n";
      prev_Simulation_Type=-1;
      tmp.str("");
      tmp_eq.str("");
      for(i = 0;i < block_triangular.ModelBlock->Size;i++)
        {
          k=block_triangular.ModelBlock->Block_List[i].Simulation_Type;
          if ((BlockTriangular::BlockSim(prev_Simulation_Type)!=BlockTriangular::BlockSim(k))  &&
              ((prev_Simulation_Type==EVALUATE_FORWARD || prev_Simulation_Type==EVALUATE_BACKWARD || prev_Simulation_Type==EVALUATE_FORWARD_R || prev_Simulation_Type==EVALUATE_BACKWARD_R)
               || (k==EVALUATE_FORWARD || k==EVALUATE_BACKWARD || k==EVALUATE_FORWARD_R || k==EVALUATE_BACKWARD_R)))
            {
              mDynamicModelFile << "    y_index_eq=[" << tmp_eq.str() << "];\n";
              tmp_eq.str("");
              mDynamicModelFile << "    y_index=[" << tmp.str() << "];\n";
              tmp.str("");
              mDynamicModelFile << tmp1.str();
              tmp1.str("");
            }
          //mDynamicModelFile << "    y_index=[";
          for(int ik=0;ik<block_triangular.ModelBlock->Block_List[i].Size;ik++)
            {
              tmp << " " << block_triangular.ModelBlock->Block_List[i].Variable[ik]+1;
              tmp_eq << " " << block_triangular.ModelBlock->Block_List[i].Equation[ik]+1;
            }
          //mDynamicModelFile << " ];\n";
          if(k==EVALUATE_FORWARD || k==EVALUATE_BACKWARD || k==EVALUATE_FORWARD_R || k==EVALUATE_BACKWARD_R)
            {
              if(i==block_triangular.ModelBlock->Size-1)
                {
                  mDynamicModelFile << "    y_index_eq=[" << tmp_eq.str() << "];\n";
                  tmp_eq.str("");
                  mDynamicModelFile << "    y_index=[" << tmp.str() << "];\n";
                  tmp.str("");
                  mDynamicModelFile << tmp1.str();
                  tmp1.str("");
                }
            }
          if (BlockTriangular::BlockSim(prev_Simulation_Type)==BlockTriangular::BlockSim(k) &&
              (k==EVALUATE_FORWARD || k==EVALUATE_BACKWARD || k==EVALUATE_FORWARD_R || k==EVALUATE_BACKWARD_R))
            skip_head=true;
          else
            skip_head=false;
          switch(k)
            {
            case EVALUATE_FORWARD:
            case EVALUATE_BACKWARD:
            case EVALUATE_FORWARD_R:
            case EVALUATE_BACKWARD_R:
              if(!skip_head)
                {
                  tmp1 << "    [y, g1, g2, g3]=" << dynamic_basename << "_" << i + 1 << "(y, x, it_, 1, g1, g2, g3);\n";
                  tmp1 << "    residual(y_index_eq)=ys(y_index)-y(it_, y_index);\n";
                }
              break;
            case SOLVE_FORWARD_SIMPLE:
            case SOLVE_BACKWARD_SIMPLE:
              mDynamicModelFile << "    y_index_eq = " << block_triangular.ModelBlock->Block_List[i].Equation[0]+1 << ";\n";
              mDynamicModelFile << "    y_index = " << block_triangular.ModelBlock->Block_List[i].Variable[0]+1 << ";\n";
              mDynamicModelFile << "    [r, g1, g2, g3]=" << dynamic_basename << "_" << i + 1 << "(y, x, it_, g1, g2, g3, y_index, 1);\n";
              mDynamicModelFile << "    residual(y_index_eq)=r;\n";
              tmp_eq.str("");
              tmp.str("");
              break;
            case SOLVE_FORWARD_COMPLETE:
            case SOLVE_BACKWARD_COMPLETE:
              mDynamicModelFile << "    y_index_eq = [" << tmp_eq.str() << "];\n";
              mDynamicModelFile << "    y_index = [" << tmp.str() << "];\n";
              tmp_i=block_triangular.ModelBlock->Block_List[i].Max_Lag+block_triangular.ModelBlock->Block_List[i].Max_Lead+1;
              mDynamicModelFile << "    y_index_c = y_index;\n";
              mDynamicModelFile << "    for i=1:" << tmp_i-1 << ",\n";
              mDynamicModelFile << "      y_index_c = [y_index_c (y_index+i*y_size)];\n";
              mDynamicModelFile << "    end;\n";
              //tmp_i=variable_table.max_lag+variable_table.max_lead+1;
              mDynamicModelFile << "    ga = [];\n";
              mDynamicModelFile << "    ga=spalloc(" << block_triangular.ModelBlock->Block_List[i].Size << ", " << block_triangular.ModelBlock->Block_List[i].Size*tmp_i << ", " << block_triangular.ModelBlock->Block_List[i].Size*block_triangular.ModelBlock->Block_List[i].Size*tmp_i << ");\n";
              //mDynamicModelFile << "    [r, g1, g2, g3, b]=" << dynamic_basename << "_" << i + 1 << "(y, x, it_, 1, ga, g2, g3);\n";
              mDynamicModelFile << "    [r, ga, g2, g3, b]=" << dynamic_basename << "_" << i + 1 << "(y, x, it_, 1, ga, g2, g3);\n";
              mDynamicModelFile << "    g1(y_index_eq,y_index_c) = ga;\n";
              mDynamicModelFile << "    residual(y_index_eq)=r;\n";
              break;
            case SOLVE_TWO_BOUNDARIES_COMPLETE:
            case SOLVE_TWO_BOUNDARIES_SIMPLE:
              //mDynamicModelFile << "    [r, g1, g2, g3, b]=" << dynamic_basename << "_" <<  i + 1 << "(y, x, it_, y_size, 1);\n";
              mDynamicModelFile << "    y_index_eq = [" << tmp_eq.str() << "];\n";
              mDynamicModelFile << "    y_index = [" << tmp.str() << "];\n";
              tmp.str("");
              tmp_eq.str("");
              tmp_i=variable_table.max_lag+variable_table.max_lead+1;
              mDynamicModelFile << "    ga = [];\n";
              mDynamicModelFile << "    ga=spalloc(" << block_triangular.ModelBlock->Block_List[i].Size << ", " << block_triangular.ModelBlock->Block_List[i].Size*tmp_i << ", " << block_triangular.ModelBlock->Block_List[i].Size*block_triangular.ModelBlock->Block_List[i].Size*tmp_i << ");\n";
              mDynamicModelFile << "    y_index_c = y_index;\n";
              tmp_i=block_triangular.ModelBlock->Block_List[i].Max_Lag+block_triangular.ModelBlock->Block_List[i].Max_Lead+1;
              mDynamicModelFile << "    for i=1:" << tmp_i-1 << ",\n";
              mDynamicModelFile << "      y_index_c = [y_index_c (y_index+i*y_size)];\n";
              mDynamicModelFile << "    end;\n";
              //mDynamicModelFile << "    [r, g1, g2, g3, b]=" << dynamic_basename << "_" <<  i + 1 << "(y, x, it_-1, " << block_triangular.ModelBlock->Block_List[i].Size << ", 1, ga, g2, g3);\n";
              mDynamicModelFile << "    [r, ga, g2, g3, b]=" << dynamic_basename << "_" <<  i + 1 << "(y, x, it_-1, " << block_triangular.ModelBlock->Block_List[i].Size << ", 1, ga, g2, g3);\n";
              if(block_triangular.ModelBlock->Block_List[i].Max_Lag==variable_table.max_lag && block_triangular.ModelBlock->Block_List[i].Max_Lead==variable_table.max_lead)
                mDynamicModelFile << "    g1(y_index_eq,y_index_c) = ga;\n";
              else
                mDynamicModelFile << "    g1(y_index_eq,y_index_c) = ga(:," << 1+(variable_table.max_lag-block_triangular.ModelBlock->Block_List[i].Max_Lag)*block_triangular.ModelBlock->Block_List[i].Size << ":" << (variable_table.max_lag+1+block_triangular.ModelBlock->Block_List[i].Max_Lead)*block_triangular.ModelBlock->Block_List[i].Size << ");\n";
              mDynamicModelFile << "    residual(y_index_eq)=r(:,M_.maximum_lag+1);\n";
              break;
            }
          prev_Simulation_Type=k;
        }
      if(tmp1.str().length())
        {
          mDynamicModelFile << tmp1.str();
          tmp1.str("");
        }
      mDynamicModelFile << "    varargout{1}=residual;\n";
      mDynamicModelFile << "    varargout{2}=g1;\n";
      mDynamicModelFile << "    return;\n";
      mDynamicModelFile << "  end;\n";
      mDynamicModelFile << "  %it is the deterministic simulation of the block decomposed dynamic model\n";
      mDynamicModelFile << "  if(options_.simulation_method==0)\n";
      mDynamicModelFile << "    mthd='Sparse LU';\n";
      mDynamicModelFile << "  elseif(options_.simulation_method==2)\n";
      mDynamicModelFile << "    mthd='GMRES';\n";
      mDynamicModelFile << "  elseif(options_.simulation_method==3)\n";
      mDynamicModelFile << "    mthd='BICGSTAB';\n";
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
    }
  prev_Simulation_Type=-1;
  for(i = 0;i < block_triangular.ModelBlock->Size;i++)
    {
      k = block_triangular.ModelBlock->Block_List[i].Simulation_Type;
      if (BlockTriangular::BlockSim(prev_Simulation_Type)==BlockTriangular::BlockSim(k) &&
          (k==EVALUATE_FORWARD || k==EVALUATE_BACKWARD || k==EVALUATE_FORWARD_R || k==EVALUATE_BACKWARD_R))
        skip_head=true;
      else
        skip_head=false;
      if ((k == EVALUATE_FORWARD || k == EVALUATE_FORWARD_R) && (block_triangular.ModelBlock->Block_List[i].Size))
        {
          if (!skip_head)
            {
              if (open_par)
                {
                  mDynamicModelFile << "  end\n";
                }
              mDynamicModelFile << "  Per_u_=0;\n";
              mDynamicModelFile << "  for it_ = y_kmin+1:(periods+y_kmin)\n";
              mDynamicModelFile << "    Per_y_=it_*y_size;\n";
              mDynamicModelFile << "    g1=[];g2=[];g3=[];\n";
              mDynamicModelFile << "    y=" << dynamic_basename << "_" << i + 1 << "(y, x, it_, 0, g1, g2, g3);\n";
            }
          open_par=true;
        }
      else if ((k == EVALUATE_BACKWARD || k == EVALUATE_BACKWARD_R) && (block_triangular.ModelBlock->Block_List[i].Size))
        {
          if (!skip_head)
            {
              if (open_par)
                {
                  mDynamicModelFile << "  end\n";
                }
              mDynamicModelFile << "  Per_u_=0;\n";
              mDynamicModelFile << "  for it_ = y_kmin+1:(periods+y_kmin)\n";
              mDynamicModelFile << "    Per_y_=it_*y_size;\n";
              mDynamicModelFile << "    " << dynamic_basename << "_" << i + 1 << "(y, x, it_, 0, g1, g2, g3);\n";
            }
          open_par=true;
        }
      else if ((k == SOLVE_FORWARD_SIMPLE) && (block_triangular.ModelBlock->Block_List[i].Size))
        {
          if (open_par)
            mDynamicModelFile << "  end\n";
          open_par=false;
          mDynamicModelFile << "  g1=0;\n";
          mDynamicModelFile << "  r=0;\n";
          mDynamicModelFile << "  for it_=y_kmin+1:periods+y_kmin\n";
          mDynamicModelFile << "    cvg=0;\n";
          mDynamicModelFile << "    iter=0;\n";
          mDynamicModelFile << "    Per_y_=it_*y_size;\n";
          mDynamicModelFile << "    while ~(cvg==1 | iter>maxit_),\n";
          mDynamicModelFile << "      [r, g1] = " << dynamic_basename << "_" << i + 1 << "(y, x, it_, g1, g2, g3, 1, 0);\n";
          mDynamicModelFile << "      y(it_, " << block_triangular.ModelBlock->Block_List[i].Variable[0]+1 << ") = y(it_, " << block_triangular.ModelBlock->Block_List[i].Variable[0]+1 << ")-r/g1;\n";
          mDynamicModelFile << "      cvg=((r*r)<solve_tolf);\n";
          mDynamicModelFile << "      iter=iter+1;\n";
          mDynamicModelFile << "    end\n";
          mDynamicModelFile << "    if cvg==0\n";
          mDynamicModelFile << "      fprintf('Convergence not achieved in block " << i << ", at time %d after %d iterations\\n',it_,iter);\n";
          mDynamicModelFile << "      return;\n";
          mDynamicModelFile << "    end\n";
          mDynamicModelFile << "  end\n";
        }
      else if ((k == SOLVE_BACKWARD_SIMPLE) && (block_triangular.ModelBlock->Block_List[i].Size))
        {
          if (open_par)
            mDynamicModelFile << "  end\n";
          open_par=false;
          mDynamicModelFile << "  g1=0;\n";
          mDynamicModelFile << "  r=0;\n";
          mDynamicModelFile << "  for it_=periods+y_kmin:-1:y_kmin+1\n";
          mDynamicModelFile << "    cvg=0;\n";
          mDynamicModelFile << "    iter=0;\n";
          mDynamicModelFile << "    Per_y_=it_*y_size;\n";
          mDynamicModelFile << "    while ~(cvg==1 | iter>maxit_),\n";
          mDynamicModelFile << "      [r, g1] = " << dynamic_basename << "_" << i + 1 << "(y, x, it_, g1, g2, g3, 1, 0);\n";
          mDynamicModelFile << "      y(it_, " << block_triangular.ModelBlock->Block_List[i].Variable[0]+1 << ") = y(it_, " << block_triangular.ModelBlock->Block_List[i].Variable[0]+1 << ")-r/g1;\n";
          mDynamicModelFile << "      cvg=((r*r)<solve_tolf);\n";
          mDynamicModelFile << "      iter=iter+1;\n";
          mDynamicModelFile << "    end\n";
          mDynamicModelFile << "    if cvg==0\n";
          mDynamicModelFile << "      fprintf('Convergence not achieved in block " << i << ", at time %d after %d iterations\\n',it_,iter);\n";
          mDynamicModelFile << "      return;\n";
          mDynamicModelFile << "    end\n";
          mDynamicModelFile << "  end\n";
        }
      /*else if ((k == SOLVE_TWO_BOUNDARIES_SIMPLE) && (block_triangular.ModelBlock->Block_List[i].Size))
        {
          if (open_par)
            mDynamicModelFile << "  end\n";
          open_par=false;
          if (!printed)
            {
              printed = true;
            }
          SGE.SGE_compute(block_triangular.ModelBlock, i, true, basename, symbol_table.endo_nbr);
          Nb_SGE++;
#ifdef PRINT_OUT
          cout << "end of Gaussian elimination\n";
#endif
          mDynamicModelFile << "    Read_file(\"" << reform(basename) << "\",periods," <<
            block_triangular.ModelBlock->Block_List[i].IM_lead_lag[block_triangular.ModelBlock->Block_List[i].Max_Lag + block_triangular.ModelBlock->Block_List[i].Max_Lead].u_finish + 1 << ", " << symbol_table.endo_nbr <<
            ", " << block_triangular.ModelBlock->Block_List[i].Max_Lag << ", " << block_triangular.ModelBlock->Block_List[i].Max_Lead << ");\n";
          cerr << "Not implemented block SOLVE_TWO_BOUNDARIES_COMPLETE" << endl;
          exit(-1);
        }*/
      else if ((k == SOLVE_FORWARD_COMPLETE) && (block_triangular.ModelBlock->Block_List[i].Size))
        {
          /*if (open_par)
            mDynamicModelFile << "  end\n";
          open_par=false;
          if (!printed)
            {
              printed = true;
            }
          SGE.SGE_compute(block_triangular.ModelBlock, i, false, basename, symbol_table.endo_nbr);
          Nb_SGE++;
          mDynamicModelFile << "    Read_file(\"" << reform(basename) << "\", periods, 0, " << symbol_table.endo_nbr << ", " << block_triangular.ModelBlock->Block_List[i].Max_Lag << ", " << block_triangular.ModelBlock->Block_List[i].Max_Lead << " );\n";*/
          if (open_par)
            mDynamicModelFile << "  end\n";
          open_par=false;
          mDynamicModelFile << "  g1=0;\n";
          mDynamicModelFile << "  r=0;\n";
          tmp_eq.str("");
          for(int ik=0;ik<block_triangular.ModelBlock->Block_List[i].Size;ik++)
            {
              tmp_eq << " " << block_triangular.ModelBlock->Block_List[i].Equation[ik]+1;
            }
          mDynamicModelFile << "  y_index_eq = [" << tmp_eq.str() << "];\n";
          mDynamicModelFile << "  for it_=periods+y_kmin:-1:y_kmin+1\n";
          mDynamicModelFile << "    cvg=0;\n";
          mDynamicModelFile << "    iter=0;\n";
          mDynamicModelFile << "    Per_y_=it_*y_size;\n";
          mDynamicModelFile << "    while ~(cvg==1 | iter>maxit_),\n";
          mDynamicModelFile << "      [r, g1] = " << dynamic_basename << "_" << i + 1 << "(y, x, it_, 0, g1, g2, g3);\n";
          mDynamicModelFile << "      y(it_, y_index_eq) = y(it_, y_index_eq)-r/g1;\n";
          mDynamicModelFile << "      cvg=((r'*r)<solve_tolf);\n";
          mDynamicModelFile << "      iter=iter+1;\n";
          mDynamicModelFile << "    end\n";
          mDynamicModelFile << "    if cvg==0\n";
          mDynamicModelFile << "      fprintf('Convergence not achieved in block " << i << ", at time %d after %d iterations\\n',it_,iter);\n";
          mDynamicModelFile << "      return;\n";
          mDynamicModelFile << "    end\n";
          mDynamicModelFile << "  end\n";
          /*cerr << "Not implemented block SOLVE_FORWARD_COMPLETE" << endl;
          exit(-1);*/
        }
      else if ((k == SOLVE_BACKWARD_COMPLETE) && (block_triangular.ModelBlock->Block_List[i].Size))
        {
          /*if (open_par)
            mDynamicModelFile << "  end\n";
          open_par=false;
          SGE.SGE_compute(block_triangular.ModelBlock, i, false, basename, symbol_table.endo_nbr);
          Nb_SGE++;
          mDynamicModelFile << "    Read_file(\"" << reform(basename) << "\", periods, 0, " << symbol_table.endo_nbr << ", " << block_triangular.ModelBlock->Block_List[i].Max_Lag << ", " << block_triangular.ModelBlock->Block_List[i].Max_Lead << " );\n";*/
          if (open_par)
            mDynamicModelFile << "  end\n";
          open_par=false;
          mDynamicModelFile << "  g1=0;\n";
          mDynamicModelFile << "  r=0;\n";
          tmp_eq.str("");
          for(int ik=0;ik<block_triangular.ModelBlock->Block_List[i].Size;ik++)
            {
              tmp_eq << " " << block_triangular.ModelBlock->Block_List[i].Equation[ik]+1;
            }
          mDynamicModelFile << "  y_index_eq = [" << tmp_eq.str() << "];\n";
          mDynamicModelFile << "  for it_=periods+y_kmin:-1:y_kmin+1\n";
          mDynamicModelFile << "    cvg=0;\n";
          mDynamicModelFile << "    iter=0;\n";
          mDynamicModelFile << "    Per_y_=it_*y_size;\n";
          mDynamicModelFile << "    while ~(cvg==1 | iter>maxit_),\n";
          mDynamicModelFile << "      [r, g1] = " << dynamic_basename << "_" << i + 1 << "(y, x, it_, 0, g1, g2, g3);\n";
          mDynamicModelFile << "      y(it_, y_index_eq) = y(it_, y_index_eq)-r/g1;\n";
          mDynamicModelFile << "      cvg=((r'*r)<solve_tolf);\n";
          mDynamicModelFile << "      iter=iter+1;\n";
          mDynamicModelFile << "    end\n";
          mDynamicModelFile << "    if cvg==0\n";
          mDynamicModelFile << "      fprintf('Convergence not achieved in block " << i << ", at time %d after %d iterations\\n',it_,iter);\n";
          mDynamicModelFile << "      return;\n";
          mDynamicModelFile << "    end\n";
          mDynamicModelFile << "  end\n";
          /*cerr << "Not implemented block SOLVE_BACKWARD_COMPLETE" << endl;
          exit(-1);*/
        }
      else if ((k == SOLVE_TWO_BOUNDARIES_COMPLETE || k == SOLVE_TWO_BOUNDARIES_SIMPLE) && (block_triangular.ModelBlock->Block_List[i].Size))
        {
          if (open_par)
            mDynamicModelFile << "  end\n";
          open_par=false;
          if (!printed)
            {
              printed = true;
            }
          Nb_SGE++;
          //cout << "new_SGE=" << new_SGE << "\n";

          mDynamicModelFile << "  cvg=0;\n";
          mDynamicModelFile << "  iter=0;\n";
          mDynamicModelFile << "  Per_u_=0;\n";
          mDynamicModelFile << "  y_index=[";
          for(int ik=0;ik<block_triangular.ModelBlock->Block_List[i].Size;ik++)
            {
              mDynamicModelFile << " " << block_triangular.ModelBlock->Block_List[i].Variable[ik]+1;
            }
          mDynamicModelFile << "  ];\n";
          mDynamicModelFile << "  Blck_size=" << block_triangular.ModelBlock->Block_List[i].Size << ";\n";
          /*mDynamicModelFile << "  if(options_.simulation_method==2 | options_.simulation_method==3),\n";
            mDynamicModelFile << "    [r, g1]= " << basename << "_static(y, x);\n";
            mDynamicModelFile << "    [L1,U1] = lu(g1,1e-5);\n";
            mDynamicModelFile << "    I = speye(periods);\n";
            mDynamicModelFile << "    L1=kron(I,L1);\n";
            mDynamicModelFile << "    U1=kron(I,U1);\n";
            mDynamicModelFile << "  end;\n";*/
          mDynamicModelFile << "  y_kmin_l=" << block_triangular.ModelBlock->Block_List[i].Max_Lag << ";\n";
          mDynamicModelFile << "  y_kmax_l=" << block_triangular.ModelBlock->Block_List[i].Max_Lead << ";\n";
          mDynamicModelFile << "  lambda=options_.slowc;\n";
          mDynamicModelFile << "  correcting_factor=0.01;\n";
          mDynamicModelFile << "  luinc_tol=1e-10;\n";
          mDynamicModelFile << "  max_resa=1e100;\n";
          int nze, m;
          for(nze=0,m=0;m<=block_triangular.ModelBlock->Block_List[i].Max_Lead+block_triangular.ModelBlock->Block_List[i].Max_Lag;m++)
            nze+=block_triangular.ModelBlock->Block_List[i].IM_lead_lag[m].size;
          mDynamicModelFile << "  Jacobian_Size=" << block_triangular.ModelBlock->Block_List[i].Size << "*(y_kmin+" << block_triangular.ModelBlock->Block_List[i].Max_Lead << " +periods);\n";
          mDynamicModelFile << "  g1=spalloc( length(y_index)*periods, Jacobian_Size, " << nze << "*periods" << ");\n";
          /*mDynamicModelFile << "  cpath=path;\n";
            mDynamicModelFile << "  addpath(fullfile(matlabroot,'toolbox','matlab','sparfun'));\n";*/
          mDynamicModelFile << "  bicgstabh=@bicgstab;\n";
          //mDynamicModelFile << "  path(cpath);\n";
          mDynamicModelFile << sp << "  reduced = 0;\n";
          //mDynamicModelFile << "  functions(bicgstabh)\n";
          if (!block_triangular.ModelBlock->Block_List[i].is_linear)
            {
              sp="  ";
              mDynamicModelFile << "  while ~(cvg==1 | iter>maxit_),\n";
            }
          else
            {
              sp="";
            }
          mDynamicModelFile << sp << "  [r, g1, g2, g3, b]=" << dynamic_basename << "_" <<  i + 1 << "(y, x, y_kmin, Blck_size, periods, g1, g2, g3);\n";
          mDynamicModelFile << sp << "  g1a=g1(:, y_kmin*Blck_size+1:(periods+y_kmin)*Blck_size);\n";
          mDynamicModelFile << sp << "  b = b' -g1(:, 1+(y_kmin-y_kmin_l)*Blck_size:y_kmin*Blck_size)*reshape(y(1+y_kmin-y_kmin_l:y_kmin,y_index)',1,y_kmin_l*Blck_size)'-g1(:, (periods+y_kmin)*Blck_size+1:(periods+y_kmin+y_kmax_l)*Blck_size)*reshape(y(periods+y_kmin+1:periods+y_kmin+y_kmax_l,y_index)',1,y_kmax_l*Blck_size)';\n";
          mDynamicModelFile << sp << "  if(~isreal(r))\n";
          mDynamicModelFile << sp << "    max_res=(-(max(max(abs(r))))^2)^0.5;\n";
          mDynamicModelFile << sp << "  else\n";
          mDynamicModelFile << sp << "    max_res=max(max(abs(r)));\n";
          mDynamicModelFile << sp << "  end;\n";
          mDynamicModelFile << sp << "  if(iter>0)\n";
          mDynamicModelFile << sp << "    if(~isreal(max_res) | isnan(max_res) | (max_resa<max_res && iter>1))\n";
          mDynamicModelFile << sp << "      if(isnan(max_res))\n";
          mDynamicModelFile << sp << "        detJ=det(g1aa);\n";
          mDynamicModelFile << sp << "        if(abs(detJ)<1e-7)\n";
          mDynamicModelFile << sp << "          max_factor=max(max(abs(g1aa)));\n";
          mDynamicModelFile << sp << "          ze_elem=sum(diag(g1aa)<options_.cutoff);\n";
          mDynamicModelFile << sp << "          disp([num2str(full(ze_elem),'%d') ' elements on the Jacobian diagonal are below the cutoff (' num2str(options_.cutoff,'%f') ')']);\n";
          mDynamicModelFile << sp << "          if(correcting_factor<max_factor)\n";
          mDynamicModelFile << sp << "            correcting_factor=correcting_factor*4;\n";
          mDynamicModelFile << sp << "            disp(['The Jacobain matrix is singular, det(Jacobian)=' num2str(detJ,'%f') '.']);\n";
          mDynamicModelFile << sp << "            disp(['    trying to correct the Jacobian matrix:']);\n";
          mDynamicModelFile << sp << "            disp(['    correcting_factor=' num2str(correcting_factor,'%f') ' max(Jacobian)=' num2str(full(max_factor),'%f')]);\n";
          mDynamicModelFile << sp << "            dx = (g1aa+correcting_factor*speye(periods*Blck_size))\\ba- ya;\n";
          mDynamicModelFile << sp << "            y(1+y_kmin:periods+y_kmin,y_index)=reshape((ya_save+lambda*dx)',length(y_index),periods)';\n";
          if (!block_triangular.ModelBlock->Block_List[i].is_linear)
            mDynamicModelFile << sp << "            continue;\n";
          mDynamicModelFile << sp << "          else\n";
          mDynamicModelFile << sp << "            disp('The singularity of the jacobian matrix could not be corrected');\n";
          mDynamicModelFile << sp << "            return;\n";
          mDynamicModelFile << sp << "          end;\n";
          mDynamicModelFile << sp << "        end;\n";
          mDynamicModelFile << sp << "      elseif(lambda>1e-8)\n";
          mDynamicModelFile << sp << "        lambda=lambda/2;\n";
          mDynamicModelFile << sp << "        reduced = 1;\n";
          mDynamicModelFile << sp << "        disp(['reducing the path length: lambda=' num2str(lambda,'%f')]);\n";
          mDynamicModelFile << sp << "        y(1+y_kmin:periods+y_kmin,y_index)=reshape((ya_save+lambda*dx)',length(y_index),periods)';\n";
          if (!block_triangular.ModelBlock->Block_List[i].is_linear)
            mDynamicModelFile << sp << "        continue;\n";
          mDynamicModelFile << sp << "      else\n";
          mDynamicModelFile << sp << "        disp(['No convergence after ' num2str(iter,'%d') ' iterations']);\n";
          mDynamicModelFile << sp << "        return;\n";
          mDynamicModelFile << sp << "      end;\n";
          mDynamicModelFile << sp << "    else\n";
          mDynamicModelFile << sp << "      if(lambda<1)\n";
          mDynamicModelFile << sp << "        lambda=max(lambda*2, 1);\n";
          mDynamicModelFile << sp << "      end;\n";
          mDynamicModelFile << sp << "    end;\n";
          mDynamicModelFile << sp << "  end;\n";

          mDynamicModelFile << sp << "  ya = reshape(y(y_kmin+1:y_kmin+periods,y_index)',1,periods*Blck_size)';\n";
          mDynamicModelFile << sp << "  ya_save=ya;\n";
          mDynamicModelFile << sp << "  g1aa=g1a;\n";
          mDynamicModelFile << sp << "  ba=b;\n";
          mDynamicModelFile << sp << "  max_resa=max_res;\n";
          mDynamicModelFile << sp << "  if(options_.simulation_method==0),\n";
          mDynamicModelFile << sp << "    dx = g1a\\b- ya;\n";
          mDynamicModelFile << sp << "    ya = ya + lambda*dx;\n";
          mDynamicModelFile << sp << "    y(1+y_kmin:periods+y_kmin,y_index)=reshape(ya',length(y_index),periods)';\n";
          mDynamicModelFile << sp << "  elseif(options_.simulation_method==2),\n";
          mDynamicModelFile << sp << "    flag1=1;\n";
          mDynamicModelFile << sp << "    while(flag1>0)\n";
          mDynamicModelFile << sp << "      [L1, U1]=luinc(g1a,luinc_tol);\n";
          mDynamicModelFile << sp << "      [za,flag1] = gmres(g1a,b," << block_triangular.ModelBlock->Block_List[i].Size << ",1e-6," << block_triangular.ModelBlock->Block_List[i].Size << "*periods,L1,U1);\n";
          mDynamicModelFile << sp << "      if (flag1>0 | reduced)\n";
          mDynamicModelFile << sp << "        if(flag1==1)\n";
          mDynamicModelFile << sp << "          disp(['No convergence inside GMRES after ' num2str(periods*" <<  block_triangular.ModelBlock->Block_List[i].Size << ",'%6d') ' iterations']);\n";
          mDynamicModelFile << sp << "        elseif(flag1==2)\n";
          mDynamicModelFile << sp << "          disp(['Preconditioner is ill-conditioned ']);\n";
          mDynamicModelFile << sp << "        elseif(flag1==3)\n";
          mDynamicModelFile << sp << "          disp(['GMRES stagnated. (Two consecutive iterates were the same.)']);\n";
          mDynamicModelFile << sp << "        end;\n";
          mDynamicModelFile << sp << "        luinc_tol = luinc_tol/10;\n";
          mDynamicModelFile << sp << "        reduced = 0;\n";
          mDynamicModelFile << sp << "      else\n";
          mDynamicModelFile << sp << "        dx = za - ya;\n";
          mDynamicModelFile << sp << "        ya = ya + lambda*dx;\n";
          mDynamicModelFile << sp << "        y(1+y_kmin:periods+y_kmin,y_index)=reshape(ya',length(y_index),periods)';\n";
          mDynamicModelFile << sp << "      end;\n";
          mDynamicModelFile << sp << "    end;\n";
          mDynamicModelFile << sp << "  elseif(options_.simulation_method==3),\n";
          mDynamicModelFile << sp << "    flag1=1;\n";
          mDynamicModelFile << sp << "    while(flag1>0)\n";
          mDynamicModelFile << sp << "      [L1, U1]=luinc(g1a,luinc_tol);\n";
          mDynamicModelFile << sp << "      [za,flag1] = bicgstabh(g1a,b,1e-7," << block_triangular.ModelBlock->Block_List[i].Size << "*periods,L1,U1);\n";
          mDynamicModelFile << sp << "      if (flag1>0 | reduced)\n";
          mDynamicModelFile << sp << "        if(flag1==1)\n";
          mDynamicModelFile << sp << "          disp(['No convergence inside BICGSTAB after ' num2str(periods*" <<  block_triangular.ModelBlock->Block_List[i].Size << ",'%6d') ' iterations']);\n";
          mDynamicModelFile << sp << "        elseif(flag1==2)\n";
          mDynamicModelFile << sp << "          disp(['Preconditioner is ill-conditioned ']);\n";
          mDynamicModelFile << sp << "        elseif(flag1==3)\n";
          mDynamicModelFile << sp << "          disp(['BICGSTAB stagnated. (Two consecutive iterates were the same.)']);\n";
          mDynamicModelFile << sp << "        end;\n";
          mDynamicModelFile << sp << "        luinc_tol = luinc_tol/10;\n";
          mDynamicModelFile << sp << "        reduced = 0;\n";
          mDynamicModelFile << sp << "      else\n";
          mDynamicModelFile << sp << "        dx = za - ya;\n";
          mDynamicModelFile << sp << "        ya = ya + lambda*dx;\n";
          mDynamicModelFile << sp << "        y(1+y_kmin:periods+y_kmin,y_index)=reshape(ya',length(y_index),periods)';\n";
          mDynamicModelFile << sp << "      end;\n";
          mDynamicModelFile << sp << "    end;\n";
          mDynamicModelFile << sp << "  end;\n";
          if(!block_triangular.ModelBlock->Block_List[i].is_linear)
            {
              mDynamicModelFile << "    cvg=(max_res<solve_tolf);\n";
              mDynamicModelFile << "    iter=iter+1;\n";
            }
          mDynamicModelFile << "    disp(['iteration: ' num2str(iter,'%d') ' error: ' num2str(max_res,'%e')]);\n";
          if(!block_triangular.ModelBlock->Block_List[i].is_linear)
            {
              mDynamicModelFile << "  end\n";
              mDynamicModelFile << "  if (iter>maxit_)\n";
              mDynamicModelFile << "    disp(['No convergence after ' num2str(iter,'%4d') ' iterations']);\n";
              mDynamicModelFile << "    return;\n";
              mDynamicModelFile << "  end;\n";
            }
        }

      prev_Simulation_Type=k;
    }
  if(open_par)
    mDynamicModelFile << "  end;\n";
  open_par=false;
  mDynamicModelFile << "  oo_.endo_simul = y';\n";
  mDynamicModelFile << "return;\n";

  writeModelEquationsOrdered_M(mDynamicModelFile, block_triangular.ModelBlock, dynamic_basename);

  mDynamicModelFile.close();

  if (printed)
    cout << "done\n";
}

void
ModelTree::writeDynamicModel(ostream &DynamicOutput) const
{
  ostringstream lsymetric;       // Used when writing symetric elements in Hessian
  ostringstream model_output;    // Used for storing model equations
  ostringstream jacobian_output; // Used for storing jacobian equations
  ostringstream hessian_output;  // Used for storing Hessian equations
  ostringstream third_derivatives_output;

  ExprNodeOutputType output_type = (mode == eStandardMode || mode==eSparseMode ? oMatlabDynamicModel : oCDynamicModel);

  writeModelLocalVariables(model_output, output_type);

  writeTemporaryTerms(model_output, output_type);

  writeModelEquations(model_output, output_type);

  int nrows = equations.size();
  int nvars = variable_table.getDynJacobianColsNbr(computeJacobianExo);
  int nvars_sq = nvars * nvars;

  // Writing Jacobian
  if (computeJacobian || computeJacobianExo)
    for(first_derivatives_type::const_iterator it = first_derivatives.begin();
        it != first_derivatives.end(); it++)
      {
        int eq = it->first.first;
        int var = it->first.second;
        NodeID d1 = it->second;

        if (computeJacobianExo || variable_table.getType(var) == eEndogenous)
          {
            ostringstream g1;
            g1 << "  g1";
            matrixHelper(g1, eq, variable_table.getDynJacobianCol(var), output_type);

            jacobian_output << g1.str() << "=" << g1.str() << "+";
            d1->writeOutput(jacobian_output, output_type, temporary_terms);
            jacobian_output << ";" << endl;
          }
      }

  // Writing Hessian
  if (computeHessian)
    for(second_derivatives_type::const_iterator it = second_derivatives.begin();
        it != second_derivatives.end(); it++)
      {
        int eq = it->first.first;
        int var1 = it->first.second.first;
        int var2 = it->first.second.second;
        NodeID d2 = it->second;

        int id1 = variable_table.getDynJacobianCol(var1);
        int id2 = variable_table.getDynJacobianCol(var2);

        int col_nb = id1*nvars+id2;
        int col_nb_sym = id2*nvars+id1;

        hessian_output << "  g2";
        matrixHelper(hessian_output, eq, col_nb, output_type);
        hessian_output << " = ";
        d2->writeOutput(hessian_output, output_type, temporary_terms);
        hessian_output << ";" << endl;

        // Treating symetric elements
        if (id1 != id2)
          {
            lsymetric <<  "  g2";
            matrixHelper(lsymetric, eq, col_nb_sym, output_type);
            lsymetric << " = " <<  "g2";
            matrixHelper(lsymetric, eq, col_nb, output_type);
            lsymetric << ";" << endl;
          }
      }

  // Writing third derivatives
  if (computeThirdDerivatives)
    for(third_derivatives_type::const_iterator it = third_derivatives.begin();
        it != third_derivatives.end(); it++)
      {
        int eq = it->first.first;
        int var1 = it->first.second.first;
        int var2 = it->first.second.second.first;
        int var3 = it->first.second.second.second;
        NodeID d3 = it->second;

        int id1 = variable_table.getDynJacobianCol(var1);
        int id2 = variable_table.getDynJacobianCol(var2);
        int id3 = variable_table.getDynJacobianCol(var3);

        // Reference column number for the g3 matrix
        int ref_col = id1 * nvars_sq + id2 * nvars + id3;

        third_derivatives_output << "  g3";
        matrixHelper(third_derivatives_output, eq, ref_col, output_type);
        third_derivatives_output << " = ";
        d3->writeOutput(third_derivatives_output, output_type, temporary_terms);
        third_derivatives_output << ";" << endl;

        // Compute the column numbers for the 5 other permutations of (id1,id2,id3) and store them in a set (to avoid duplicates if two indexes are equal)
        set<int> cols;
        cols.insert(id1 * nvars_sq + id3 * nvars + id2);
        cols.insert(id2 * nvars_sq + id1 * nvars + id3);
        cols.insert(id2 * nvars_sq + id3 * nvars + id1);
        cols.insert(id3 * nvars_sq + id1 * nvars + id2);
        cols.insert(id3 * nvars_sq + id2 * nvars + id1);

        for(set<int>::iterator it2 = cols.begin(); it2 != cols.end(); it2++)
          if (*it2 != ref_col)
            {
              third_derivatives_output << "  g3";
              matrixHelper(third_derivatives_output, eq, *it2, output_type);
              third_derivatives_output << " = " << "g3";
              matrixHelper(third_derivatives_output, eq, ref_col, output_type);
              third_derivatives_output << ";" << endl;
            }
      }

  if (mode == eStandardMode)
    {
      DynamicOutput << "%" << endl
                    << "% Model equations" << endl
                    << "%" << endl
                    << endl
                    << "residual = zeros(" << nrows << ", 1);" << endl
                    << model_output.str();

      if (computeJacobian || computeJacobianExo)
        {
          // Writing initialization instruction for matrix g1
          DynamicOutput << "if nargout >= 2," << endl
                        << "  g1 = zeros(" << nrows << ", " << nvars << ");" << endl
                        << endl
                        << "%" << endl
                        << "% Jacobian matrix" << endl
                        << "%" << endl
                        << endl
                        << jacobian_output.str()
                        << "end" << endl;
        }
      if (computeHessian)
        {
          // Writing initialization instruction for matrix g2
          int ncols = nvars_sq;
          DynamicOutput << "if nargout >= 3," << endl
                        << "  g2 = sparse([],[],[], " << nrows << ", " << ncols << ", " << 5*ncols << ");" << endl
                        << endl
                        << "%" << endl
                        << "% Hessian matrix" << endl
                        << "%" << endl
                        << endl
                        << hessian_output.str()
                        << lsymetric.str()
                        << "end;" << endl;
        }
      if (computeThirdDerivatives)
        {
          int ncols = nvars_sq * nvars;
          DynamicOutput << "if nargout >= 4," << endl
                        << "  g3 = sparse([],[],[], " << nrows << ", " << ncols << ", " << 5*ncols << ");" << endl
                        << endl
                        << "%" << endl
                        << "% Third order derivatives" << endl
                        << "%" << endl
                        << endl
                        << third_derivatives_output.str()
                        << "end;" << endl;
        }
    }
  else
    {
      DynamicOutput << "void Dynamic(double *y, double *x, int nb_row_x, double *params, int it_, double *residual, double *g1, double *g2)" << endl
                    << "{" << endl
                    << "  double lhs, rhs;" << endl
                    << endl
                    << "  /* Residual equations */" << endl
                    << model_output.str();

      if (computeJacobian || computeJacobianExo)
        {
          DynamicOutput << "  /* Jacobian  */" << endl
                        << "  if (g1 == NULL)" << endl
                        << "    return;" << endl
                        << "  else" << endl
                        << "    {" << endl
                        << jacobian_output.str()
                        << "    }" << endl;
        }
      if (computeHessian)
        {
          DynamicOutput << "  /* Hessian for endogenous and exogenous variables */" << endl
                        << "  if (g2 == NULL)" << endl
                        << "    return;" << endl
                        << "  else" << endl
                        << "    {" << endl
                        << hessian_output.str()
                        << lsymetric.str()
                        << "    }" << endl;
        }
      DynamicOutput << "}" << endl << endl;
    }
}

void
ModelTree::writeOutput(ostream &output) const
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
  for(int endoID = 0; endoID < symbol_table.endo_nbr; endoID++)
    {
      output << "\n\t";
      // Loop on periods
      for(int lag = -variable_table.max_endo_lag; lag <= variable_table.max_endo_lead; lag++)
        {
          // Print variableID if exists with current period, otherwise print 0
          try
            {
              int varID = variable_table.getID(eEndogenous, endoID, lag);
              output << " " << variable_table.getDynJacobianCol(varID) + 1;
            }
          catch(VariableTable::UnknownVariableKeyException &e)
            {
              output << " 0";
            }
        }
      output << ";";
    }
  output << "]';\n";
  //In case of sparse model, writes the block structure of the model

  if(mode==eSparseMode || mode==eSparseDLLMode)
    {
      int prev_Simulation_Type=-1;
      bool skip_the_head;
      int k=0;
      int count_lead_lag_incidence = 0;
      int max_lead, max_lag;
      for(int j = 0;j < block_triangular.ModelBlock->Size;j++)
        {
          //For a block composed of a single equation determines wether we have to evaluate or to solve the equation
          if (BlockTriangular::BlockSim(prev_Simulation_Type)==BlockTriangular::BlockSim(block_triangular.ModelBlock->Block_List[j].Simulation_Type)
              && (block_triangular.ModelBlock->Block_List[j].Simulation_Type==EVALUATE_BACKWARD
              ||block_triangular.ModelBlock->Block_List[j].Simulation_Type==EVALUATE_FORWARD
              ||block_triangular.ModelBlock->Block_List[j].Simulation_Type==EVALUATE_BACKWARD_R
              ||block_triangular.ModelBlock->Block_List[j].Simulation_Type==EVALUATE_FORWARD_R ))
             skip_the_head=true;
          else
            {
              skip_the_head=false;
              k++;
              count_lead_lag_incidence = 0;
              int Block_size=block_triangular.ModelBlock->Block_List[j].Size;
              max_lag =block_triangular.ModelBlock->Block_List[j].Max_Lag ;
              max_lead=block_triangular.ModelBlock->Block_List[j].Max_Lead;
              bool evaluate=false;
              ostringstream tmp_s, tmp_s_eq;
              tmp_s.str("");
              tmp_s_eq.str("");
              for(int i=0;i<block_triangular.ModelBlock->Block_List[j].Size;i++)
                {
                  tmp_s << " " << block_triangular.ModelBlock->Block_List[j].Variable[i]+1;
                  tmp_s_eq << " " << block_triangular.ModelBlock->Block_List[j].Equation[i]+1;
                }
              if (block_triangular.ModelBlock->Block_List[j].Simulation_Type==EVALUATE_BACKWARD
                ||block_triangular.ModelBlock->Block_List[j].Simulation_Type==EVALUATE_FORWARD
                ||block_triangular.ModelBlock->Block_List[j].Simulation_Type==EVALUATE_BACKWARD_R
                ||block_triangular.ModelBlock->Block_List[j].Simulation_Type==EVALUATE_FORWARD_R
                && j+Block_size<(block_triangular.ModelBlock->Size))
                {
                  bool OK=true;
                  evaluate=true;
                  while(j+Block_size<(block_triangular.ModelBlock->Size) && OK)
                    {
                      if(BlockTriangular::BlockSim(block_triangular.ModelBlock->Block_List[j].Simulation_Type)!=BlockTriangular::BlockSim(block_triangular.ModelBlock->Block_List[j+Block_size].Simulation_Type))
                        OK=false;
                      else
                        {
                          if(max_lag <block_triangular.ModelBlock->Block_List[j+Block_size].Max_Lag )
                            max_lag =block_triangular.ModelBlock->Block_List[j+Block_size].Max_Lag ;
                          if(max_lead<block_triangular.ModelBlock->Block_List[j+Block_size].Max_Lead)
                            max_lead=block_triangular.ModelBlock->Block_List[j+Block_size].Max_Lead;
                          //cout << "block_triangular.ModelBlock->Block_List[" << j+Block_size << "].Size=" << block_triangular.ModelBlock->Block_List[j+Block_size].Size << "\n";
                          for(int i=0;i<block_triangular.ModelBlock->Block_List[j+Block_size].Size;i++)
                            {
                              tmp_s << " " << block_triangular.ModelBlock->Block_List[j+Block_size].Variable[i]+1;
                              tmp_s_eq << " " << block_triangular.ModelBlock->Block_List[j+Block_size].Equation[i]+1;
                            }
                          Block_size+=block_triangular.ModelBlock->Block_List[j+Block_size].Size;
                        }
                      //cout << "i=" << i << " max_lag=" << max_lag << " max_lead=" << max_lead << "\n";
                    }
                }
              output << "M_.block_structure.block(" << k << ").num = " << j+1 << ";\n";
              //output << "M_.block_structure.block(" << k << ").size = " << block_triangular.ModelBlock->Block_List[j].Size << ";\n";
              output << "M_.block_structure.block(" << k << ").Simulation_Type = " << block_triangular.ModelBlock->Block_List[j].Simulation_Type << ";\n";
              output << "M_.block_structure.block(" << k << ").maximum_endo_lag = " << max_lag << ";\n";
              output << "M_.block_structure.block(" << k << ").maximum_endo_lead = " << max_lead << ";\n";
              output << "M_.block_structure.block(" << k << ").endo_nbr = " << Block_size << ";\n";
              output << "M_.block_structure.block(" << k << ").equation = [" << tmp_s_eq.str() << "];\n";
              output << "M_.block_structure.block(" << k << ").variable = [" << tmp_s.str() << "];\n";
              tmp_s.str("");
              cout << "begining of lead_lag_incidence\n";

              bool done_IM=false;
              if(!evaluate)
                {
                  output << "M_.block_structure.block(" << k << ").lead_lag_incidence = [];\n";
                  for(int l=-max_lag;l<max_lead+1;l++)
                    {
                      bool *tmp_IM;
                      tmp_IM=block_triangular.bGet_IM(l);
                      for(int l_var=0;l_var<block_triangular.ModelBlock->Block_List[j].Size;l_var++)
                        {
                          for(int l_equ=0;l_equ<block_triangular.ModelBlock->Block_List[j].Size;l_equ++)
                            if(tmp_IM[block_triangular.ModelBlock->Block_List[j].Equation[l_equ]*symbol_table.endo_nbr+block_triangular.ModelBlock->Block_List[j].Variable[l_var]])
                              {
                                count_lead_lag_incidence++;
                                if(tmp_s.str().length())
                                  tmp_s << " ";
                                tmp_s << count_lead_lag_incidence;
                                done_IM=true;
                                break;
                              }
                          if(!done_IM)
                            tmp_s << " 0";
                          done_IM=false;
                        }
                       output << "M_.block_structure.block(" << k << ").lead_lag_incidence = [ M_.block_structure.block(" << k << ").lead_lag_incidence; " << tmp_s.str() << "];\n";
                       tmp_s.str("");
                    }
                }
              else
                {
                  output << "M_.block_structure.block(" << k << ").lead_lag_incidence = [\n";
                  for(int l=-max_lag;l<max_lead+1;l++)
                    {
                      bool not_increm=true;
                      bool *tmp_IM;
                      tmp_IM=block_triangular.bGet_IM(l);
                      int ii=j;
                      while(ii-j<Block_size)
                        {
                          for(int l_var=0;l_var<block_triangular.ModelBlock->Block_List[ii].Size;l_var++)
                            {
                              for(int l_equ=0;l_equ<block_triangular.ModelBlock->Block_List[ii].Size;l_equ++)
                                if(tmp_IM[block_triangular.ModelBlock->Block_List[ii].Equation[l_equ]*symbol_table.endo_nbr+block_triangular.ModelBlock->Block_List[ii].Variable[l_var]])
                                  {
                                    //if(not_increm && l==-max_lag)
                                      count_lead_lag_incidence++;
                                    not_increm=false;
                                    if(tmp_s.str().length())
                                      tmp_s << " ";
                                    //tmp_s << count_lead_lag_incidence+(l+max_lag)*Block_size;
                                    tmp_s << count_lead_lag_incidence;
                                    done_IM=true;
                                    break;
                                  }
                              if(!done_IM)
                                tmp_s << " 0";
                              done_IM=false;
                            }
                          ii++;
                        }
                      output << tmp_s.str() << "\n";
                      tmp_s.str("");
                    }
                  output << "];\n";
                }
            }
          prev_Simulation_Type=block_triangular.ModelBlock->Block_List[j].Simulation_Type;

        }
      for(int j=-block_triangular.Model_Max_Lag;j<block_triangular.Model_Max_Lead+1;j++)
        {
          bool* IM = block_triangular.bGet_IM(j);
          if(IM)
            {
              bool new_entry=true;
              output << "M_.block_structure.incidence(" << block_triangular.Model_Max_Lag+j+1 << ").sparse_IM = [";
              for(int i=0;i<symbol_table.endo_nbr*symbol_table.endo_nbr;i++)
                {
                  if(IM[i])
                    {
                      if(!new_entry)
                        output << " ; ";
                      else
                        output << " ";
                      output << i/symbol_table.endo_nbr+1 << " " << i % symbol_table.endo_nbr+1;
                      new_entry=false;
                    }
                }
              output << "];\n";
            }
        }
    }
  // Writing initialization for some other variables
  output << "M_.exo_names_orig_ord = [1:" << symbol_table.exo_nbr << "];\n";
  output << "M_.maximum_lag = " << variable_table.max_lag << ";\n";
  output << "M_.maximum_lead = " << variable_table.max_lead << ";\n";
  if (symbol_table.endo_nbr)
    {
      output << "M_.maximum_endo_lag = " << variable_table.max_endo_lag << ";\n";
      output << "M_.maximum_endo_lead = " << variable_table.max_endo_lead << ";\n";
      output << "oo_.steady_state = zeros(" << symbol_table.endo_nbr << ", 1);\n";
    }
  if (symbol_table.exo_nbr)
    {
      output << "M_.maximum_exo_lag = " << variable_table.max_exo_lag << ";\n";
      output << "M_.maximum_exo_lead = " << variable_table.max_exo_lead << ";\n";
      output << "oo_.exo_steady_state = zeros(" << symbol_table.exo_nbr << ", 1);\n";
    }
  if (symbol_table.exo_det_nbr)
    {
      output << "M_.maximum_exo_det_lag = " << variable_table.max_exo_det_lag << ";\n";
      output << "M_.maximum_exo_det_lead = " << variable_table.max_exo_det_lead << ";\n";
      output << "oo_.exo_det_steady_state = zeros(" << symbol_table.exo_det_nbr << ", 1);\n";
    }
  if (symbol_table.parameter_nbr)
    output << "M_.params = repmat(NaN," << symbol_table.parameter_nbr << ", 1);\n";
}

void
ModelTree::addEquation(NodeID eq)
{
  BinaryOpNode *beq = dynamic_cast<BinaryOpNode *>(eq);

  if (beq == NULL || beq->op_code != oEqual)
    {
      cerr << "ModelTree::addEquation: you didn't provide an equal node!" << endl;
      exit(-1);
    }

  equations.push_back(beq);
}

void
ModelTree::evaluateJacobian(const eval_context_type &eval_context, jacob_map *j_m)
{
  int i=0;
  int j=0;
  bool *IM=NULL;
  int a_variable_lag=-9999;
  //block_triangular.Print_IM(2);
  for(first_derivatives_type::iterator it = first_derivatives.begin();
      it != first_derivatives.end(); it++)
    {
      if (variable_table.getType(it->first.second) == eEndogenous)
        {
          NodeID Id = it->second;
          double val;
          try
            {
              val = Id->eval(eval_context);
            }
          catch(ExprNode::EvalException &e)
            {
              cerr << "ModelTree::evaluateJacobian: evaluation of Jacobian failed!" << endl;
            }
          int eq=it->first.first;
          int var=variable_table.getSymbolID(it->first.second);
          int k1=variable_table.getLag(it->first.second);
          if (a_variable_lag!=k1)
            {
              IM=block_triangular.bGet_IM(k1);
              a_variable_lag=k1;
            }
          if (k1==0)
            {
              j++;
              (*j_m)[make_pair(eq,var)]=val;
            }
          if (IM[eq*symbol_table.endo_nbr+var] && (fabs(val) < cutoff))
            {
              if(block_triangular.bt_verbose)
                cout << "the coefficient related to variable " << var << " with lag " << k1 << " in equation " << eq << " is equal to " << val << " and is set to 0 in the incidence matrix (size=" << symbol_table.endo_nbr << ")\n";
              block_triangular.unfill_IM(eq, var, k1);
              i++;
            }
        }
    }
  if (i>0)
    {
      cout << i << " elements among " << first_derivatives.size() << " in the incidence matrices are below the cutoff (" << cutoff << ") and are discarded\n";
      cout << "the contemporaneous incidence matrix has " << j << " elements\n";
    }
}

void
ModelTree::BlockLinear(Model_Block *ModelBlock)
{
  int i,j,l,m,ll;
  for(j = 0;j < ModelBlock->Size;j++)
    {
      if (ModelBlock->Block_List[j].Simulation_Type==SOLVE_BACKWARD_COMPLETE ||
          ModelBlock->Block_List[j].Simulation_Type==SOLVE_FORWARD_COMPLETE)
        {
          ll=ModelBlock->Block_List[j].Max_Lag;
          for(i=0;i<ModelBlock->Block_List[j].IM_lead_lag[ll].size;i++)
            {
              int eq=ModelBlock->Block_List[j].IM_lead_lag[ll].Equ_Index[i];
              int var=ModelBlock->Block_List[j].IM_lead_lag[ll].Var_Index[i];
              first_derivatives_type::const_iterator it=first_derivatives.find(make_pair(eq,variable_table.getID(eEndogenous,var,0)));
              if (it!= first_derivatives.end())
                {
                  NodeID Id = it->second;
                  set<pair<int, int> > endogenous;
                  Id->collectEndogenous(endogenous);
                  if (endogenous.size() > 0)
                    {
                      for(l=0;l<ModelBlock->Block_List[j].Size;l++)
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
      else if (ModelBlock->Block_List[j].Simulation_Type==SOLVE_TWO_BOUNDARIES_COMPLETE)
        {
          for(m=0;m<=ModelBlock->Block_List[j].Max_Lead+ModelBlock->Block_List[j].Max_Lag;m++)
            {
              int k1=m-ModelBlock->Block_List[j].Max_Lag;
              for(i=0;i<ModelBlock->Block_List[j].IM_lead_lag[m].size;i++)
                {
                  int eq=ModelBlock->Block_List[j].IM_lead_lag[m].Equ_Index[i];
                  int var=ModelBlock->Block_List[j].IM_lead_lag[m].Var_Index[i];
                  first_derivatives_type::const_iterator it=first_derivatives.find(make_pair(eq,variable_table.getID(eEndogenous,var,k1)));
                  NodeID Id = it->second;
                  if (it!= first_derivatives.end())
                    {
                      set<pair<int, int> > endogenous;
                      Id->collectEndogenous(endogenous);
                      if (endogenous.size() > 0)
                        {
                          for(l=0;l<ModelBlock->Block_List[j].Size;l++)
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

void
ModelTree::computingPass(const eval_context_type &eval_context)
{
  cout << equations.size() << " equation(s) found" << endl;

  // Computes dynamic jacobian columns
  variable_table.computeDynJacobianCols();

  // Determine derivation order
  int order = 1;
  if (computeThirdDerivatives)
    order = 3;
  else if (computeHessian || computeStaticHessian)
    order = 2;

  // Launch computations
  derive(order);

  if (mode == eSparseDLLMode || mode == eSparseMode)
    {
      jacob_map j_m;

      evaluateJacobian(eval_context, &j_m);

      if (block_triangular.bt_verbose)
        {
          cout << "The gross incidence matrix \n";
          block_triangular.Print_IM( symbol_table.endo_nbr);
        }
      block_triangular.Normalize_and_BlockDecompose_Static_0_Model(j_m);
      BlockLinear(block_triangular.ModelBlock);

      computeTemporaryTermsOrdered(order, block_triangular.ModelBlock);
    }
  else
    computeTemporaryTerms(order);
}

void
ModelTree::writeStaticFile(const string &basename) const
{
  switch(mode)
    {
    case eStandardMode:
    case eSparseDLLMode:
      writeStaticMFile(basename + "_static");
      break;
    case eSparseMode:
      writeSparseStaticMFile(basename + "_static", basename, mode);
      break;
    case eDLLMode:
      writeStaticCFile(basename + "_static");
      break;
    }
}

void
ModelTree::writeDynamicFile(const string &basename) const
{
  switch(mode)
    {
    case eStandardMode:
      writeDynamicMFile(basename + "_dynamic");
      break;
    case eSparseMode:
      writeSparseDynamicMFile(basename + "_dynamic", basename, mode);
      block_triangular.Free_Block(block_triangular.ModelBlock);
      block_triangular.Free_IM(block_triangular.First_IM);
      break;
    case eDLLMode:
      writeDynamicCFile(basename + "_dynamic");
      break;
    case eSparseDLLMode:
      writeModelEquationsCodeOrdered(basename + "_dynamic", block_triangular.ModelBlock, basename, oCDynamicModelSparseDLL);
      block_triangular.Free_Block(block_triangular.ModelBlock);
      block_triangular.Free_IM(block_triangular.First_IM);
      break;
    }
}

void
ModelTree::matrixHelper(ostream &output, int eq_nb, int col_nb, ExprNodeOutputType output_type) const
{
  output << LPAR(output_type);
  if (OFFSET(output_type))
    output << eq_nb + 1 << ", " << col_nb + 1;
  else
    output << eq_nb + col_nb * equations.size();
  output << RPAR(output_type);
}
