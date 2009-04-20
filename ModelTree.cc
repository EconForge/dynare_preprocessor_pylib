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

#include <cstdlib>
#include <iostream>

#include "ModelTree.hh"

ModelTree::ModelTree(SymbolTable &symbol_table_arg,
                     NumericalConstants &num_constants_arg) :
  DataTree(symbol_table_arg, num_constants_arg),
  mode(eStandardMode)
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
  first_derivatives_type::const_iterator it = first_derivatives.find(make_pair(eq, getDerivID(symb_id, lag)));
  if (it != first_derivatives.end())
    (it->second)->writeOutput(output, output_type, temporary_terms);
  else
    output << 0;
}

void
ModelTree::computeJacobian(const set<int> &vars)
{
  for(set<int>::const_iterator it = vars.begin();
      it != vars.end(); it++)
    for (int eq = 0; eq < (int) equations.size(); eq++)
      {
        NodeID d1 = equations[eq]->getDerivative(*it);
        if (d1 == Zero)
          continue;
        first_derivatives[make_pair(eq, *it)] = d1;
      }
}

void
ModelTree::computeHessian(const set<int> &vars)
{
  for (first_derivatives_type::const_iterator it = first_derivatives.begin();
       it != first_derivatives.end(); it++)
    {
      int eq = it->first.first;
      int var1 = it->first.second;
      NodeID d1 = it->second;

      // Store only second derivatives with var2 <= var1
      for(set<int>::const_iterator it2 = vars.begin();
          it2 != vars.end(); it2++)
        {
          int var2 = *it2;
          if (var2 > var1)
            continue;

          NodeID d2 = d1->getDerivative(var2);
          if (d2 == Zero)
            continue;
          second_derivatives[make_pair(eq, make_pair(var1, var2))] = d2;
        }
    }
}

void
ModelTree::computeThirdDerivatives(const set<int> &vars)
{
  for (second_derivatives_type::const_iterator it = second_derivatives.begin();
       it != second_derivatives.end(); it++)
    {
      int eq = it->first.first;

      int var1 = it->first.second.first;
      int var2 = it->first.second.second;
      // By construction, var2 <= var1

      NodeID d2 = it->second;

      // Store only third derivatives such that var3 <= var2 <= var1
      for(set<int>::const_iterator it2 = vars.begin();
          it2 != vars.end(); it2++)
        {
          int var3 = *it2;
          if (var3 > var2)
            continue;

          NodeID d3 = d2->getDerivative(var3);
          if (d3 == Zero)
            continue;
          third_derivatives[make_pair(eq, make_pair(var1, make_pair(var2, var3)))] = d3;
        }
    }
}

void
ModelTree::computeTemporaryTerms()
{
  map<NodeID, int> reference_count;
  temporary_terms.clear();

  bool is_matlab = (mode != eDLLMode);

  for (vector<BinaryOpNode *>::iterator it = equations.begin();
       it != equations.end(); it++)
    (*it)->computeTemporaryTerms(reference_count, temporary_terms, is_matlab);

  for (first_derivatives_type::iterator it = first_derivatives.begin();
       it != first_derivatives.end(); it++)
    it->second->computeTemporaryTerms(reference_count, temporary_terms, is_matlab);

  for (second_derivatives_type::iterator it = second_derivatives.begin();
       it != second_derivatives.end(); it++)
    it->second->computeTemporaryTerms(reference_count, temporary_terms, is_matlab);

  for (third_derivatives_type::iterator it = third_derivatives.begin();
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

  for (temporary_terms_type::const_iterator it = temporary_terms.begin();
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
  for (map<int, NodeID>::const_iterator it = local_variables_table.begin();
       it != local_variables_table.end(); it++)
    {
      int id = it->first;
      NodeID value = it->second;

      if (!OFFSET(output_type))
        output << "double ";

      output << symbol_table.getName(id) << " = ";
      // Use an empty set for the temporary terms
      value->writeOutput(output, output_type, temporary_terms_type());
      output << ";" << endl;
    }
}

void
ModelTree::writeModelEquations(ostream &output, ExprNodeOutputType output_type) const
{
  for (int eq = 0; eq < (int) equations.size(); eq++)
    {
      BinaryOpNode *eq_node = equations[eq];

      NodeID lhs = eq_node->get_arg1();
      output << "lhs =";
      lhs->writeOutput(output, output_type, temporary_terms);
      output << ";" << endl;

      NodeID rhs = eq_node->get_arg2();
      output << "rhs =";
      rhs->writeOutput(output, output_type, temporary_terms);
      output << ";" << endl;

      output << "residual" << LPAR(output_type) << eq + OFFSET(output_type) << RPAR(output_type) << "= lhs-rhs;" << endl;
    }
}

void
ModelTree::addEquation(NodeID eq)
{
  BinaryOpNode *beq = dynamic_cast<BinaryOpNode *>(eq);

  if (beq == NULL || beq->get_op_code() != oEqual)
    {
      cerr << "ModelTree::addEquation: you didn't provide an equal node!" << endl;
      exit(EXIT_FAILURE);
    }

  equations.push_back(beq);
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
