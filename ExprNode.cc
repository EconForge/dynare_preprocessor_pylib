/*
 * Copyright (C) 2007-2009 Dynare Team
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
#include <iterator>
#include <algorithm>

// For select1st()
#ifdef __GNUC__
# include <ext/functional>
using namespace __gnu_cxx;
#endif

#include <cassert>
#include <cmath>

#include "ExprNode.hh"
#include "DataTree.hh"
#include "BlockTriangular.hh"

ExprNode::ExprNode(DataTree &datatree_arg) : datatree(datatree_arg)
{
  // Add myself to datatree
  datatree.node_list.push_back(this);

  // Set my index and increment counter
  idx = datatree.node_counter++;
}

ExprNode::~ExprNode()
{
}

NodeID
ExprNode::getDerivative(int deriv_id)
{
  // Return zero if derivative is necessarily null (using symbolic a priori)
  set<int>::const_iterator it = non_null_derivatives.find(deriv_id);
  if (it == non_null_derivatives.end())
    return datatree.Zero;

  // If derivative is stored in cache, use the cached value, otherwise compute it (and cache it)
  map<int, NodeID>::const_iterator it2 = derivatives.find(deriv_id);
  if (it2 != derivatives.end())
    return it2->second;
  else
    {
      NodeID d = computeDerivative(deriv_id);
      derivatives[deriv_id] = d;
      return d;
    }
}

int
ExprNode::precedence(ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms) const
  {
    // For a constant, a variable, or a unary op, the precedence is maximal
    return 100;
  }

int
ExprNode::cost(const temporary_terms_type &temporary_terms, bool is_matlab) const
  {
    // For a terminal node, the cost is null
    return 0;
  }

void
ExprNode::collectEndogenous(set<pair<int, int> > &result) const
{
  set<pair<int, int> > symb_ids;
  collectVariables(eEndogenous, symb_ids);
  for(set<pair<int, int> >::const_iterator it = symb_ids.begin();
      it != symb_ids.end(); it++)
    result.insert(make_pair(datatree.symbol_table.getTypeSpecificID(it->first), it->second));
}

void
ExprNode::collectExogenous(set<pair<int, int> > &result) const
{
  set<pair<int, int> > symb_ids;
  collectVariables(eExogenous, symb_ids);
  for(set<pair<int, int> >::const_iterator it = symb_ids.begin();
      it != symb_ids.end(); it++)
    result.insert(make_pair(datatree.symbol_table.getTypeSpecificID(it->first), it->second));
}

void
ExprNode::collectModelLocalVariables(set<int> &result) const
{
  set<pair<int, int> > symb_ids;
  collectVariables(eModelLocalVariable, symb_ids);
  transform(symb_ids.begin(), symb_ids.end(), inserter(result, result.begin()),
            select1st<pair<int, int> >());
}

void
ExprNode::computeTemporaryTerms(map<NodeID, int> &reference_count,
                                temporary_terms_type &temporary_terms,
                                bool is_matlab) const
  {
    // Nothing to do for a terminal node
  }

void
ExprNode::computeTemporaryTerms(map<NodeID, int> &reference_count,
                                temporary_terms_type &temporary_terms,
                                map<NodeID, pair<int, int> > &first_occurence,
                                int Curr_block,
                                Model_Block *ModelBlock,
                                int equation,
                                map_idx_type &map_idx) const
  {
    // Nothing to do for a terminal node
  }


pair<int, NodeID >
ExprNode::normalizeEquation(int var_endo, vector<pair<int, pair<NodeID, NodeID> > > &List_of_Op_RHS) const
  {
    return(make_pair(0, (NodeID)NULL));
  }


void
ExprNode::writeOutput(ostream &output)
{
  writeOutput(output, oMatlabOutsideModel, temporary_terms_type());
}


NumConstNode::NumConstNode(DataTree &datatree_arg, int id_arg) :
    ExprNode(datatree_arg),
    id(id_arg)
{
  // Add myself to the num const map
  datatree.num_const_node_map[id] = this;

  // All derivatives are null, so non_null_derivatives is left empty
}

NodeID
NumConstNode::computeDerivative(int deriv_id)
{
  return datatree.Zero;
}

void
NumConstNode::collectTemporary_terms(const temporary_terms_type &temporary_terms, Model_Block *ModelBlock, int Curr_Block) const
  {
    temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<NumConstNode *>(this));
    if (it != temporary_terms.end())
      ModelBlock->Block_List[Curr_Block].Temporary_InUse->insert(idx);
  }

void
NumConstNode::writeOutput(ostream &output, ExprNodeOutputType output_type,
                          const temporary_terms_type &temporary_terms) const
  {
    //cout << "writeOutput constante\n";
    temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<NumConstNode *>(this));
    if (it != temporary_terms.end())
      if (output_type == oMatlabDynamicModelSparse)
        output << "T" << idx << "(it_)";
      else
        output << "T" << idx;
    else
      output << datatree.num_constants.get(id);
  }

double
NumConstNode::eval(const eval_context_type &eval_context) const throw (EvalException)
{
  return(datatree.num_constants.getDouble(id));
}

void
NumConstNode::compile(ostream &CompileCode, bool lhs_rhs, const temporary_terms_type &temporary_terms, map_idx_type &map_idx, bool dynamic, bool steady_dynamic) const
  {
    CompileCode.write(&FLDC, sizeof(FLDC));
    double vard = datatree.num_constants.getDouble(id);
    CompileCode.write(reinterpret_cast<char *>(&vard),sizeof(vard));
  }

void
NumConstNode::collectVariables(SymbolType type_arg, set<pair<int, int> > &result) const
{
}

pair<int, NodeID >
NumConstNode::normalizeEquation(int var_endo, vector<pair<int, pair<NodeID, NodeID> > > &List_of_Op_RHS) const
  {
    return(make_pair(0, datatree.AddNumConstant(datatree.num_constants.get(id))));
  }

NodeID
NumConstNode::getChainRuleDerivative(int deriv_id, const map<int, NodeID> &recursive_variables)
{
  return datatree.Zero;
}

NodeID
NumConstNode::toStatic(DataTree &static_datatree) const
  {
    return static_datatree.AddNumConstant(datatree.num_constants.get(id));
  }


VariableNode::VariableNode(DataTree &datatree_arg, int symb_id_arg, int lag_arg, int deriv_id_arg) :
    ExprNode(datatree_arg),
    symb_id(symb_id_arg),
    type(datatree.symbol_table.getType(symb_id_arg)),
    lag(lag_arg),
    deriv_id(deriv_id_arg)
{
  // Add myself to the variable map
  datatree.variable_node_map[make_pair(symb_id, lag)] = this;

  // It makes sense to allow a lead/lag on parameters: during steady state calibration, endogenous and parameters can be swapped
  assert(lag == 0 || (type != eModelLocalVariable && type != eModFileLocalVariable && type != eUnknownFunction));

  // Fill in non_null_derivatives
  switch (type)
    {
    case eEndogenous:
    case eExogenous:
    case eExogenousDet:
    case eParameter:
      // For a variable or a parameter, the only non-null derivative is with respect to itself
      non_null_derivatives.insert(deriv_id);
      break;
    case eModelLocalVariable:
      // Non null derivatives are those of the value of the local parameter
      non_null_derivatives = datatree.local_variables_table[symb_id]->non_null_derivatives;
      break;
    case eModFileLocalVariable:
      // Such a variable is never derived
      break;
    case eUnknownFunction:
      cerr << "Attempt to construct a VariableNode with an unknown function name" << endl;
      exit(EXIT_FAILURE);
    }
}

NodeID
VariableNode::computeDerivative(int deriv_id_arg)
{
  switch (type)
    {
    case eEndogenous:
    case eExogenous:
    case eExogenousDet:
    case eParameter:
      if (deriv_id == deriv_id_arg)
        return datatree.One;
      else
        return datatree.Zero;
    case eModelLocalVariable:
      return datatree.local_variables_table[symb_id]->getDerivative(deriv_id_arg);
    case eModFileLocalVariable:
      cerr << "ModFileLocalVariable is not derivable" << endl;
      exit(EXIT_FAILURE);
    case eUnknownFunction:
      cerr << "Impossible case!" << endl;
      exit(EXIT_FAILURE);
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

void
VariableNode::collectTemporary_terms(const temporary_terms_type &temporary_terms, Model_Block *ModelBlock, int Curr_Block) const
  {
    temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<VariableNode *>(this));
    if (it != temporary_terms.end())
      ModelBlock->Block_List[Curr_Block].Temporary_InUse->insert(idx);
    if (type== eModelLocalVariable)
      datatree.local_variables_table[symb_id]->collectTemporary_terms(temporary_terms, ModelBlock, Curr_Block);
  }

void
VariableNode::writeOutput(ostream &output, ExprNodeOutputType output_type,
                          const temporary_terms_type &temporary_terms) const
  {
    // If node is a temporary term
    temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<VariableNode *>(this));
    if (it != temporary_terms.end())
      {
        if (output_type == oMatlabDynamicModelSparse)
          output << "T" << idx << "(it_)";
        else
          output << "T" << idx;
        return;
      }

    if (IS_LATEX(output_type))
      {
        if (output_type == oLatexDynamicSteadyStateOperator)
          output << "\\bar{";
        output << datatree.symbol_table.getTeXName(symb_id);
        if (output_type == oLatexDynamicModel
            && (type == eEndogenous || type == eExogenous || type == eExogenousDet || type == eModelLocalVariable))
          {
            output << "_{t";
            if (lag != 0)
              {
                if (lag > 0)
                  output << "+";
                output << lag;
              }
            output << "}";
          }
        else if (output_type == oLatexDynamicSteadyStateOperator)
          output << "}";
        return;
      }

    int i;
    int tsid = datatree.symbol_table.getTypeSpecificID(symb_id);
    switch (type)
      {
      case eParameter:
        if (output_type == oMatlabOutsideModel)
          output << "M_.params" << "(" << tsid + 1 << ")";
        else
          output << "params" << LEFT_ARRAY_SUBSCRIPT(output_type) << tsid + ARRAY_SUBSCRIPT_OFFSET(output_type) << RIGHT_ARRAY_SUBSCRIPT(output_type);
        break;

      case eModelLocalVariable:
      case eModFileLocalVariable:
        if (output_type==oMatlabDynamicModelSparse || output_type==oMatlabStaticModelSparse)
          {
            output << "(";
            datatree.local_variables_table[symb_id]->writeOutput(output, output_type,temporary_terms);
            output << ")";
          }
        else
          output << datatree.symbol_table.getName(symb_id);
        break;

      case eEndogenous:
        switch (output_type)
          {
          case oMatlabDynamicModel:
          case oCDynamicModel:
            i = datatree.getDynJacobianCol(deriv_id) + ARRAY_SUBSCRIPT_OFFSET(output_type);
            output <<  "y" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << RIGHT_ARRAY_SUBSCRIPT(output_type);
            break;
          case oMatlabStaticModel:
          case oMatlabStaticModelSparse:
            i = tsid + ARRAY_SUBSCRIPT_OFFSET(output_type);
            output <<  "y" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << RIGHT_ARRAY_SUBSCRIPT(output_type);
            break;
          case oMatlabDynamicModelSparse:
            i = tsid + ARRAY_SUBSCRIPT_OFFSET(output_type);
            if (lag > 0)
              output << "y" << LEFT_ARRAY_SUBSCRIPT(output_type) << "it_+" << lag << ", " << i << RIGHT_ARRAY_SUBSCRIPT(output_type);
            else if (lag < 0)
              output << "y" << LEFT_ARRAY_SUBSCRIPT(output_type) << "it_" << lag << ", " << i << RIGHT_ARRAY_SUBSCRIPT(output_type);
            else
              output << "y" << LEFT_ARRAY_SUBSCRIPT(output_type) << "it_, " << i << RIGHT_ARRAY_SUBSCRIPT(output_type);
            break;
          case oMatlabOutsideModel:
            output << "oo_.steady_state(" << tsid + 1 << ")";
            break;
          case oMatlabDynamicSteadyStateOperator:
            output << "oo_.steady_state(" << tsid + 1 << ")";
            break;
          default:
            assert(false);
          }
        break;

      case eExogenous:
        i = tsid + ARRAY_SUBSCRIPT_OFFSET(output_type);
        switch (output_type)
          {
          case oMatlabDynamicModel:
          case oMatlabDynamicModelSparse:
            if (lag > 0)
              output <<  "x(it_+" << lag << ", " << i << ")";
            else if (lag < 0)
              output <<  "x(it_" << lag << ", " << i << ")";
            else
              output <<  "x(it_, " << i << ")";
            break;
          case oCDynamicModel:
            if (lag == 0)
              output <<  "x[it_+" << i << "*nb_row_x]";
            else if (lag > 0)
              output <<  "x[it_+" << lag << "+" << i << "*nb_row_x]";
            else
              output <<  "x[it_" << lag << "+" << i << "*nb_row_x]";
            break;
          case oMatlabStaticModel:
          case oMatlabStaticModelSparse:
            output << "x" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << RIGHT_ARRAY_SUBSCRIPT(output_type);
            break;
          case oMatlabOutsideModel:
            assert(lag == 0);
            output <<  "oo_.exo_steady_state(" << i << ")";
            break;
          case oMatlabDynamicSteadyStateOperator:
            output <<  "oo_.exo_steady_state(" << i << ")";
            break;
          default:
            assert(false);
          }
        break;

      case eExogenousDet:
        i = tsid + datatree.symbol_table.exo_nbr() + ARRAY_SUBSCRIPT_OFFSET(output_type);
        switch (output_type)
          {
          case oMatlabDynamicModel:
          case oMatlabDynamicModelSparse:
            if (lag > 0)
              output <<  "x(it_+" << lag << ", " << i << ")";
            else if (lag < 0)
              output <<  "x(it_" << lag << ", " << i << ")";
            else
              output <<  "x(it_, " << i << ")";
            break;
          case oCDynamicModel:
            if (lag == 0)
              output <<  "x[it_+" << i << "*nb_row_xd]";
            else if (lag > 0)
              output <<  "x[it_+" << lag << "+" << i << "*nb_row_xd]";
            else
              output <<  "x[it_" << lag << "+" << i << "*nb_row_xd]";
            break;
          case oMatlabStaticModel:
          case oMatlabStaticModelSparse:
            output << "x" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << RIGHT_ARRAY_SUBSCRIPT(output_type);
            break;
          case oMatlabOutsideModel:
            assert(lag == 0);
            output <<  "oo_.exo_det_steady_state(" << tsid + 1 << ")";
            break;
		  case oMatlabDynamicSteadyStateOperator:
			output <<  "oo_.exo_det_steady_state(" << tsid + 1 << ")";
			break;
          default:
            assert(false);
          }
        break;

      case eUnknownFunction:
        cerr << "Impossible case" << endl;
        exit(EXIT_FAILURE);
      }
  }

double
VariableNode::eval(const eval_context_type &eval_context) const throw (EvalException)
{
  eval_context_type::const_iterator it = eval_context.find(symb_id);
  if (it == eval_context.end())
    throw EvalException();

  return it->second;
}

void
VariableNode::compile(ostream &CompileCode, bool lhs_rhs, const temporary_terms_type &temporary_terms, map_idx_type &map_idx, bool dynamic, bool steady_dynamic) const
  {
    int i, lagl;
    if (!lhs_rhs)
      {
        if(dynamic)
          {
            if(steady_dynamic)  // steady state values in a dynamic model
              CompileCode.write(&FLDVS, sizeof(FLDVS));
            else
              CompileCode.write(&FLDV, sizeof(FLDV));
          }
        else
          CompileCode.write(&FLDSV, sizeof(FLDSV));
      }
    else
      {
        if(dynamic)
          {
            if(steady_dynamic)  // steady state values in a dynamic model
              {
                /*CompileCode.write(&FLDVS, sizeof(FLDVS));*/
                cerr << "Impossible case: steady_state in rhs of equation" << endl;
                exit(EXIT_FAILURE);
              }
            else
              CompileCode.write(&FSTPV, sizeof(FSTPV));
          }
        else
          CompileCode.write(&FSTPSV, sizeof(FSTPSV));
      }
    char typel=(char)type;
    CompileCode.write(&typel, sizeof(typel));
    int tsid = datatree.symbol_table.getTypeSpecificID(symb_id);
    switch (type)
      {
      case eParameter:
        //cout << "Parameter=" << tsid << "\n";
        i = tsid;
        CompileCode.write(reinterpret_cast<char *>(&i), sizeof(i));
        break;
      case eEndogenous :
        //cout << "Endogenous=" << symb_id << "\n";
        i = tsid;//symb_id;
        CompileCode.write(reinterpret_cast<char *>(&i), sizeof(i));
        if(dynamic && !steady_dynamic)
          {
            lagl=lag;
            CompileCode.write(reinterpret_cast<char *>(&lagl), sizeof(lagl));
          }
        break;
      case eExogenous :
        //cout << "Exogenous=" << tsid << "\n";
        i = tsid;
        CompileCode.write(reinterpret_cast<char *>(&i), sizeof(i));
        if(dynamic && !steady_dynamic)
          {
            lagl=lag;
            CompileCode.write(reinterpret_cast<char *>(&lagl), sizeof(lagl));
          }
        break;
      case eExogenousDet:
        i = tsid + datatree.symbol_table.exo_nbr();
        //cout << "ExogenousDet=" << i << "\n";
        CompileCode.write(reinterpret_cast<char *>(&i), sizeof(i));
        if(dynamic && !steady_dynamic)
          {
            lagl=lag;
            CompileCode.write(reinterpret_cast<char *>(&lagl), sizeof(lagl));
          }
        break;
      case eModelLocalVariable:
      case eModFileLocalVariable:
        //cout << "eModelLocalVariable=" << symb_id << "\n";
        datatree.local_variables_table[symb_id]->compile(CompileCode, lhs_rhs, temporary_terms, map_idx, dynamic, steady_dynamic);
        break;
      case eUnknownFunction:
        cerr << "Impossible case: eUnknownFuncion" << endl;
        exit(EXIT_FAILURE);
      }
  }

void
VariableNode::computeTemporaryTerms(map<NodeID, int> &reference_count,
                                    temporary_terms_type &temporary_terms,
                                    map<NodeID, pair<int, int> > &first_occurence,
                                    int Curr_block,
                                    Model_Block *ModelBlock,
                                    int equation,
                                    map_idx_type &map_idx) const
  {
    if (type== eModelLocalVariable)
      datatree.local_variables_table[symb_id]->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, Curr_block, ModelBlock, equation, map_idx);
  }

void
VariableNode::collectVariables(SymbolType type_arg, set<pair<int, int> > &result) const
{
  if (type == type_arg)
    result.insert(make_pair(symb_id, lag));
  if (type == eModelLocalVariable)
    datatree.local_variables_table[symb_id]->collectVariables(type_arg, result);
}

pair<int, NodeID>
VariableNode::normalizeEquation(int var_endo, vector<pair<int, pair<NodeID, NodeID> > > &List_of_Op_RHS) const
  {
    if (type ==eEndogenous)
      {
        if (datatree.symbol_table.getTypeSpecificID(symb_id)==var_endo && lag==0)
          return(make_pair(1, (NodeID)NULL ));
        else
          return(make_pair(0, datatree.AddVariableInternal(datatree.symbol_table.getName(symb_id), lag) ));
      }
    else
      {
        if (type == eParameter)
          return(make_pair(0, datatree.AddVariableInternal(datatree.symbol_table.getName(symb_id), 0) ));
        else
          return(make_pair(0, datatree.AddVariableInternal(datatree.symbol_table.getName(symb_id), lag) ));
      }
  }

NodeID
VariableNode::getChainRuleDerivative(int deriv_id_arg, const map<int, NodeID> &recursive_variables)
{
  switch (type)
    {
    case eEndogenous:
    case eExogenous:
    case eExogenousDet:
    case eParameter:
      if (deriv_id == deriv_id_arg)
        return datatree.One;
      else
        {
          //if there is in the equation a recursive variable we could use a chaine rule derivation
          map<int, NodeID>::const_iterator it = recursive_variables.find(deriv_id);
          if (it != recursive_variables.end())
            {
              map<int, NodeID>::const_iterator it2 = derivatives.find(deriv_id_arg);
              if (it2 != derivatives.end())
                return it2->second;
              else
                {
                  map<int, NodeID> recursive_vars2(recursive_variables);
                  recursive_vars2.erase(it->first);
                  //NodeID c = datatree.AddNumConstant("1");
                  NodeID d = datatree.AddUMinus(it->second->getChainRuleDerivative(deriv_id_arg, recursive_vars2));
                  //d = datatree.AddTimes(c, d);
                  derivatives[deriv_id_arg] = d;
                  return d;
                }
            }
          else
            return datatree.Zero;
        }
    case eModelLocalVariable:
      return datatree.local_variables_table[symb_id]->getChainRuleDerivative(deriv_id_arg, recursive_variables);
    case eModFileLocalVariable:
      cerr << "ModFileLocalVariable is not derivable" << endl;
      exit(EXIT_FAILURE);
    case eUnknownFunction:
      cerr << "Impossible case!" << endl;
      exit(EXIT_FAILURE);
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}



NodeID
VariableNode::toStatic(DataTree &static_datatree) const
  {
    return static_datatree.AddVariable(datatree.symbol_table.getName(symb_id));
  }


UnaryOpNode::UnaryOpNode(DataTree &datatree_arg, UnaryOpcode op_code_arg, const NodeID arg_arg) :
    ExprNode(datatree_arg),
    arg(arg_arg),
    op_code(op_code_arg)
{
  // Add myself to the unary op map
  datatree.unary_op_node_map[make_pair(arg, op_code)] = this;

  // Non-null derivatives are those of the argument
  non_null_derivatives = arg->non_null_derivatives;
}

NodeID
UnaryOpNode::composeDerivatives(NodeID darg)
{
  NodeID t11, t12, t13;

  switch (op_code)
    {
    case oUminus:
      return datatree.AddUMinus(darg);
    case oExp:
      return datatree.AddTimes(darg, this);
    case oLog:
      return datatree.AddDivide(darg, arg);
    case oLog10:
      t11 = datatree.AddExp(datatree.One);
      t12 = datatree.AddLog10(t11);
      t13 = datatree.AddDivide(darg, arg);
      return datatree.AddTimes(t12, t13);
    case oCos:
      t11 = datatree.AddSin(arg);
      t12 = datatree.AddUMinus(t11);
      return datatree.AddTimes(darg, t12);
    case oSin:
      t11 = datatree.AddCos(arg);
      return datatree.AddTimes(darg, t11);
    case oTan:
      t11 = datatree.AddTimes(this, this);
      t12 = datatree.AddPlus(t11, datatree.One);
      return datatree.AddTimes(darg, t12);
    case oAcos:
      t11 = datatree.AddSin(this);
      t12 = datatree.AddDivide(darg, t11);
      return datatree.AddUMinus(t12);
    case oAsin:
      t11 = datatree.AddCos(this);
      return datatree.AddDivide(darg, t11);
    case oAtan:
      t11 = datatree.AddTimes(arg, arg);
      t12 = datatree.AddPlus(datatree.One, t11);
      return datatree.AddDivide(darg, t12);
    case oCosh:
      t11 = datatree.AddSinh(arg);
      return datatree.AddTimes(darg, t11);
    case oSinh:
      t11 = datatree.AddCosh(arg);
      return datatree.AddTimes(darg, t11);
    case oTanh:
      t11 = datatree.AddTimes(this, this);
      t12 = datatree.AddMinus(datatree.One, t11);
      return datatree.AddTimes(darg, t12);
    case oAcosh:
      t11 = datatree.AddSinh(this);
      return datatree.AddDivide(darg, t11);
    case oAsinh:
      t11 = datatree.AddCosh(this);
      return datatree.AddDivide(darg, t11);
    case oAtanh:
      t11 = datatree.AddTimes(arg, arg);
      t12 = datatree.AddMinus(datatree.One, t11);
      return datatree.AddTimes(darg, t12);
    case oSqrt:
      t11 = datatree.AddPlus(this, this);
      return datatree.AddDivide(darg, t11);
    case oSteadyState:
	  if (datatree.isDynamic())
		return datatree.Zero;
	  else
		return darg;
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

NodeID
UnaryOpNode::computeDerivative(int deriv_id)
{
  NodeID darg = arg->getDerivative(deriv_id);
  return composeDerivatives(darg);
}

int
UnaryOpNode::cost(const temporary_terms_type &temporary_terms, bool is_matlab) const
  {
    // For a temporary term, the cost is null
    temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<UnaryOpNode *>(this));
    if (it != temporary_terms.end())
      return 0;

    int cost = arg->cost(temporary_terms, is_matlab);

    if (is_matlab)
      // Cost for Matlab files
      switch (op_code)
        {
        case oUminus:
          return cost + 70;
        case oExp:
          return cost + 160;
        case oLog:
          return cost + 300;
        case oLog10:
          return cost + 16000;
        case oCos:
        case oSin:
        case oCosh:
          return cost + 210;
        case oTan:
          return cost + 230;
        case oAcos:
          return cost + 300;
        case oAsin:
          return cost + 310;
        case oAtan:
          return cost + 140;
        case oSinh:
          return cost + 240;
        case oTanh:
          return cost + 190;
        case oAcosh:
          return cost + 770;
        case oAsinh:
          return cost + 460;
        case oAtanh:
          return cost + 350;
        case oSqrt:
          return cost + 570;
		case oSteadyState:
          return cost;
        }
    else
      // Cost for C files
      switch (op_code)
        {
        case oUminus:
          return cost + 3;
        case oExp:
        case oAcosh:
          return cost + 210;
        case oLog:
          return cost + 137;
        case oLog10:
          return cost + 139;
        case oCos:
        case oSin:
          return cost + 160;
        case oTan:
          return cost + 170;
        case oAcos:
        case oAtan:
          return cost + 190;
        case oAsin:
          return cost + 180;
        case oCosh:
        case oSinh:
        case oTanh:
          return cost + 240;
        case oAsinh:
          return cost + 220;
        case oAtanh:
          return cost + 150;
        case oSqrt:
          return cost + 90;
		case oSteadyState:
          return cost;
        }
    // Suppress GCC warning
    exit(EXIT_FAILURE);
  }

void
UnaryOpNode::computeTemporaryTerms(map<NodeID, int> &reference_count,
                                   temporary_terms_type &temporary_terms,
                                   bool is_matlab) const
  {
    NodeID this2 = const_cast<UnaryOpNode *>(this);

    map<NodeID, int>::iterator it = reference_count.find(this2);
    if (it == reference_count.end())
      {
        reference_count[this2] = 1;
        arg->computeTemporaryTerms(reference_count, temporary_terms, is_matlab);
      }
    else
      {
        reference_count[this2]++;
        if (reference_count[this2] * cost(temporary_terms, is_matlab) > MIN_COST(is_matlab))
          temporary_terms.insert(this2);
      }
  }

void
UnaryOpNode::computeTemporaryTerms(map<NodeID, int> &reference_count,
                                   temporary_terms_type &temporary_terms,
                                   map<NodeID, pair<int, int> > &first_occurence,
                                   int Curr_block,
                                   Model_Block *ModelBlock,
                                   int equation,
                                   map_idx_type &map_idx) const
  {
    NodeID this2 = const_cast<UnaryOpNode *>(this);
    map<NodeID, int>::iterator it = reference_count.find(this2);
    if (it == reference_count.end())
      {
        reference_count[this2] = 1;
        first_occurence[this2] = make_pair(Curr_block,equation);
        arg->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, Curr_block, ModelBlock, equation, map_idx);
      }
    else
      {
        reference_count[this2]++;
        if (reference_count[this2] * cost(temporary_terms, false) > MIN_COST_C)
          {
            temporary_terms.insert(this2);
            ModelBlock->Block_List[first_occurence[this2].first].Temporary_Terms_in_Equation[first_occurence[this2].second]->insert(this2);
          }
      }
  }

void
UnaryOpNode::collectTemporary_terms(const temporary_terms_type &temporary_terms, Model_Block *ModelBlock, int Curr_Block) const
  {
    temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<UnaryOpNode*>(this));
    if (it != temporary_terms.end())
      ModelBlock->Block_List[Curr_Block].Temporary_InUse->insert(idx);
    else
      arg->collectTemporary_terms(temporary_terms, ModelBlock, Curr_Block);
  }

void
UnaryOpNode::writeOutput(ostream &output, ExprNodeOutputType output_type,
                         const temporary_terms_type &temporary_terms) const
  {
    // If node is a temporary term
    temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<UnaryOpNode *>(this));
    if (it != temporary_terms.end())
      {
        if (output_type == oMatlabDynamicModelSparse)
          output << "T" << idx << "(it_)";
        else
          output << "T" << idx;
        return;
      }

    // Always put parenthesis around uminus nodes
    if (op_code == oUminus)
      output << LEFT_PAR(output_type);

    switch (op_code)
      {
      case oUminus:
        output << "-";
        break;
      case oExp:
        output << "exp";
        break;
      case oLog:
        output << "log";
        break;
      case oLog10:
        if (IS_LATEX(output_type))
          output << "log_{10}";
        else
          output << "log10";
        break;
      case oCos:
        output << "cos";
        break;
      case oSin:
        output << "sin";
        break;
      case oTan:
        output << "tan";
        break;
      case oAcos:
        output << "acos";
        break;
      case oAsin:
        output << "asin";
        break;
      case oAtan:
        output << "atan";
        break;
      case oCosh:
        output << "cosh";
        break;
      case oSinh:
        output << "sinh";
        break;
      case oTanh:
        output << "tanh";
        break;
      case oAcosh:
        output << "acosh";
        break;
      case oAsinh:
        output << "asinh";
        break;
      case oAtanh:
        output << "atanh";
        break;
      case oSqrt:
        output << "sqrt";
        break;
      case oSteadyState:
        ExprNodeOutputType new_output_type;
        switch(output_type)
          {
          case oMatlabDynamicModel:
            new_output_type = oMatlabDynamicSteadyStateOperator;
            break;
          case oLatexDynamicModel:
            new_output_type = oLatexDynamicSteadyStateOperator;
            break;
          case oCDynamicModel:
            cerr << "Steady State Operator not implemented for oCDynamicModel." << endl;
            exit(EXIT_FAILURE);
          case oMatlabDynamicModelSparse:
            cerr << "Steady State Operator not implemented for oMatlabDynamicModelSparse." << endl;
            exit(EXIT_FAILURE);
          default:
            new_output_type = output_type;
            break;
          }
        arg->writeOutput(output, new_output_type, temporary_terms);
        return;
      }

    bool close_parenthesis = false;

    /* Enclose argument with parentheses if:
       - current opcode is not uminus, or
       - current opcode is uminus and argument has lowest precedence
    */
    if (op_code != oUminus
        || (op_code == oUminus
            && arg->precedence(output_type, temporary_terms) < precedence(output_type, temporary_terms)))
      {
        output << LEFT_PAR(output_type);
        close_parenthesis = true;
      }

    // Write argument
    arg->writeOutput(output, output_type, temporary_terms);

    if (close_parenthesis)
      output << RIGHT_PAR(output_type);

    // Close parenthesis for uminus
    if (op_code == oUminus)
      output << RIGHT_PAR(output_type);
  }

double
UnaryOpNode::eval_opcode(UnaryOpcode op_code, double v) throw (EvalException)
{
  switch (op_code)
    {
    case oUminus:
      return(-v);
    case oExp:
      return(exp(v));
    case oLog:
      return(log(v));
    case oLog10:
      return(log10(v));
    case oCos:
      return(cos(v));
    case oSin:
      return(sin(v));
    case oTan:
      return(tan(v));
    case oAcos:
      return(acos(v));
    case oAsin:
      return(asin(v));
    case oAtan:
      return(atan(v));
    case oCosh:
      return(cosh(v));
    case oSinh:
      return(sinh(v));
    case oTanh:
      return(tanh(v));
    case oAcosh:
      return(acosh(v));
    case oAsinh:
      return(asinh(v));
    case oAtanh:
      return(atanh(v));
    case oSqrt:
      return(sqrt(v));
	case oSteadyState:
      return(v);
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

double
UnaryOpNode::eval(const eval_context_type &eval_context) const throw (EvalException)
{
  double v = arg->eval(eval_context);

  return eval_opcode(op_code, v);
}

void
UnaryOpNode::compile(ostream &CompileCode, bool lhs_rhs, const temporary_terms_type &temporary_terms, map_idx_type &map_idx, bool dynamic, bool steady_dynamic) const
  {
    temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<UnaryOpNode *>(this));
    if (it != temporary_terms.end())
      {
        if(dynamic)
          CompileCode.write(&FLDT, sizeof(FLDT));
        else
          CompileCode.write(&FLDST, sizeof(FLDST));
        int var=map_idx[idx];
        CompileCode.write(reinterpret_cast<char *>(&var), sizeof(var));
        return;
      }
    if (op_code == oSteadyState)
      arg->compile(CompileCode, lhs_rhs, temporary_terms, map_idx, dynamic, true);
    else
      {
        arg->compile(CompileCode, lhs_rhs, temporary_terms, map_idx, dynamic, steady_dynamic);
        CompileCode.write(&FUNARY, sizeof(FUNARY));
        UnaryOpcode op_codel=op_code;
        CompileCode.write(reinterpret_cast<char *>(&op_codel), sizeof(op_codel));
      }
  }

void
UnaryOpNode::collectVariables(SymbolType type_arg, set<pair<int, int> > &result) const
{
  arg->collectVariables(type_arg, result);
}

pair<int, NodeID>
UnaryOpNode::normalizeEquation(int var_endo, vector<pair<int, pair<NodeID, NodeID> > > &List_of_Op_RHS) const
  {
    pair<bool, NodeID > res = arg->normalizeEquation(var_endo, List_of_Op_RHS);
    int is_endogenous_present = res.first;
    NodeID New_NodeID = res.second;
    /*if(res.second.second)*/
    if(is_endogenous_present==2)
      return(make_pair(2, (NodeID)NULL));
    else if (is_endogenous_present)
      {
        switch (op_code)
          {
          case oUminus:
            List_of_Op_RHS.push_back(make_pair(oUminus, make_pair((NodeID)NULL, (NodeID)NULL)));
            return(make_pair(1, (NodeID)NULL));
          case oExp:
            List_of_Op_RHS.push_back(make_pair(oLog, make_pair((NodeID)NULL, (NodeID)NULL)));
            return(make_pair(1, (NodeID)NULL));
          case oLog:
            List_of_Op_RHS.push_back(make_pair(oExp, make_pair((NodeID)NULL, (NodeID)NULL)));
            return(make_pair(1, (NodeID)NULL));
          case oLog10:
            List_of_Op_RHS.push_back(make_pair(oPower, make_pair((NodeID)NULL, datatree.AddNumConstant("10"))));
            return(make_pair(1, (NodeID)NULL));
          case oCos:
            return(make_pair(1, (NodeID)NULL));
          case oSin:
            return(make_pair(1, (NodeID)NULL));
          case oTan:
            return(make_pair(1, (NodeID)NULL));
          case oAcos:
            return(make_pair(1, (NodeID)NULL));
          case oAsin:
            return(make_pair(1, (NodeID)NULL));
          case oAtan:
            return(make_pair(1, (NodeID)NULL));
          case oCosh:
            return(make_pair(1, (NodeID)NULL));
          case oSinh:
            return(make_pair(1, (NodeID)NULL));
          case oTanh:
            return(make_pair(1, (NodeID)NULL));
          case oAcosh:
            return(make_pair(1, (NodeID)NULL));
          case oAsinh:
            return(make_pair(1, (NodeID)NULL));
          case oAtanh:
            return(make_pair(1, (NodeID)NULL));
          case oSqrt:
            List_of_Op_RHS.push_back(make_pair(oPower, make_pair((NodeID)NULL, datatree.AddNumConstant("2"))));
            return(make_pair(1, (NodeID)NULL));
          case oSteadyState:
            return(make_pair(1, (NodeID)NULL));
          }
      }
    else
      {
        switch (op_code)
          {
          case oUminus:
            return(make_pair(0, datatree.AddUMinus(New_NodeID)));
          case oExp:
            return(make_pair(0, datatree.AddExp(New_NodeID)));
          case oLog:
            return(make_pair(0, datatree.AddLog(New_NodeID)));
          case oLog10:
            return(make_pair(0, datatree.AddLog10(New_NodeID)));
          case oCos:
            return(make_pair(0, datatree.AddCos(New_NodeID)));
          case oSin:
            return(make_pair(0, datatree.AddSin(New_NodeID)));
          case oTan:
            return(make_pair(0, datatree.AddTan(New_NodeID)));
          case oAcos:
            return(make_pair(0, datatree.AddAcos(New_NodeID)));
          case oAsin:
            return(make_pair(0, datatree.AddAsin(New_NodeID)));
          case oAtan:
            return(make_pair(0, datatree.AddAtan(New_NodeID)));
          case oCosh:
            return(make_pair(0, datatree.AddCosh(New_NodeID)));
          case oSinh:
            return(make_pair(0, datatree.AddSinh(New_NodeID)));
          case oTanh:
            return(make_pair(0, datatree.AddTanh(New_NodeID)));
          case oAcosh:
            return(make_pair(0, datatree.AddAcosh(New_NodeID)));
          case oAsinh:
            return(make_pair(0, datatree.AddAsinh(New_NodeID)));
          case oAtanh:
            return(make_pair(0, datatree.AddAtanh(New_NodeID)));
          case oSqrt:
            return(make_pair(0, datatree.AddSqrt(New_NodeID)));
          case oSteadyState:
            return(make_pair(0, datatree.AddSteadyState(New_NodeID)));
          }
      }
    return(make_pair(1, (NodeID)NULL));
  }


NodeID
UnaryOpNode::getChainRuleDerivative(int deriv_id, const map<int, NodeID> &recursive_variables)
{
  NodeID darg = arg->getChainRuleDerivative(deriv_id, recursive_variables);
  return composeDerivatives(darg);
}

NodeID
UnaryOpNode::toStatic(DataTree &static_datatree) const
  {
    NodeID sarg = arg->toStatic(static_datatree);
    switch (op_code)
      {
      case oUminus:
        return static_datatree.AddUMinus(sarg);
      case oExp:
        return static_datatree.AddExp(sarg);
      case oLog:
        return static_datatree.AddLog(sarg);
      case oLog10:
        return static_datatree.AddLog10(sarg);
      case oCos:
        return static_datatree.AddCos(sarg);
      case oSin:
        return static_datatree.AddSin(sarg);
      case oTan:
        return static_datatree.AddTan(sarg);
      case oAcos:
        return static_datatree.AddAcos(sarg);
      case oAsin:
        return static_datatree.AddAsin(sarg);
      case oAtan:
        return static_datatree.AddAtan(sarg);
      case oCosh:
        return static_datatree.AddCosh(sarg);
      case oSinh:
        return static_datatree.AddSinh(sarg);
      case oTanh:
        return static_datatree.AddTanh(sarg);
      case oAcosh:
        return static_datatree.AddAcosh(sarg);
      case oAsinh:
        return static_datatree.AddAsinh(sarg);
      case oAtanh:
        return static_datatree.AddAtanh(sarg);
      case oSqrt:
        return static_datatree.AddSqrt(sarg);
      case oSteadyState:
        return static_datatree.AddSteadyState(sarg);
      }
    // Suppress GCC warning
    exit(EXIT_FAILURE);
  }


BinaryOpNode::BinaryOpNode(DataTree &datatree_arg, const NodeID arg1_arg,
                           BinaryOpcode op_code_arg, const NodeID arg2_arg) :
    ExprNode(datatree_arg),
    arg1(arg1_arg),
    arg2(arg2_arg),
    op_code(op_code_arg)
{
  datatree.binary_op_node_map[make_pair(make_pair(arg1, arg2), op_code)] = this;

  // Non-null derivatives are the union of those of the arguments
  // Compute set union of arg1->non_null_derivatives and arg2->non_null_derivatives
  set_union(arg1->non_null_derivatives.begin(),
            arg1->non_null_derivatives.end(),
            arg2->non_null_derivatives.begin(),
            arg2->non_null_derivatives.end(),
            inserter(non_null_derivatives, non_null_derivatives.begin()));
}

NodeID
BinaryOpNode::composeDerivatives(NodeID darg1, NodeID darg2)
{
  NodeID t11, t12, t13, t14, t15;

  switch (op_code)
    {
    case oPlus:
      return datatree.AddPlus(darg1, darg2);
    case oMinus:
      return datatree.AddMinus(darg1, darg2);
    case oTimes:
      t11 = datatree.AddTimes(darg1, arg2);
      t12 = datatree.AddTimes(darg2, arg1);
      return datatree.AddPlus(t11, t12);
    case oDivide:
      if (darg2!=datatree.Zero)
        {
          t11 = datatree.AddTimes(darg1, arg2);
          t12 = datatree.AddTimes(darg2, arg1);
          t13 = datatree.AddMinus(t11, t12);
          t14 = datatree.AddTimes(arg2, arg2);
          return datatree.AddDivide(t13, t14);
        }
      else
        return datatree.AddDivide(darg1, arg2);
    case oLess:
    case oGreater:
    case oLessEqual:
    case oGreaterEqual:
    case oEqualEqual:
    case oDifferent:
      return datatree.Zero;
    case oPower:
      if (darg2 == datatree.Zero)
        {
          if (darg1 == datatree.Zero)
            return datatree.Zero;
          else
            {
              t11 = datatree.AddMinus(arg2, datatree.One);
              t12 = datatree.AddPower(arg1, t11);
              t13 = datatree.AddTimes(arg2, t12);
              return datatree.AddTimes(darg1, t13);
            }
        }
      else
        {
          t11 = datatree.AddLog(arg1);
          t12 = datatree.AddTimes(darg2, t11);
          t13 = datatree.AddTimes(darg1, arg2);
          t14 = datatree.AddDivide(t13, arg1);
          t15 = datatree.AddPlus(t12, t14);
          return datatree.AddTimes(t15, this);
        }
    case oMax:
      t11 = datatree.AddGreater(arg1,arg2);
      t12 = datatree.AddTimes(t11,darg1);
      t13 = datatree.AddMinus(datatree.One,t11);
      t14 = datatree.AddTimes(t13,darg2);
      return datatree.AddPlus(t14,t12);
    case oMin:
      t11 = datatree.AddGreater(arg2,arg1);
      t12 = datatree.AddTimes(t11,darg1);
      t13 = datatree.AddMinus(datatree.One,t11);
      t14 = datatree.AddTimes(t13,darg2);
      return datatree.AddPlus(t14,t12);
    case oEqual:
      return datatree.AddMinus(darg1, darg2);
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

NodeID
BinaryOpNode::computeDerivative(int deriv_id)
{
  NodeID darg1 = arg1->getDerivative(deriv_id);
  NodeID darg2 = arg2->getDerivative(deriv_id);
  return composeDerivatives(darg1, darg2);
}

int
BinaryOpNode::precedence(ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms) const
  {
    temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<BinaryOpNode *>(this));
    // A temporary term behaves as a variable
    if (it != temporary_terms.end())
      return 100;

    switch (op_code)
      {
      case oEqual:
        return 0;
      case oEqualEqual:
      case oDifferent:
        return 1;
      case oLessEqual:
      case oGreaterEqual:
      case oLess:
      case oGreater:
        return 2;
      case oPlus:
      case oMinus:
        return 3;
      case oTimes:
      case oDivide:
        return 4;
      case oPower:
        if (IS_C(output_type))
          // In C, power operator is of the form pow(a, b)
          return 100;
        else
          return 5;
      case oMin:
      case oMax:
        return 100;
      }
    // Suppress GCC warning
    exit(EXIT_FAILURE);
  }

int
BinaryOpNode::cost(const temporary_terms_type &temporary_terms, bool is_matlab) const
  {
    temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<BinaryOpNode *>(this));
    // For a temporary term, the cost is null
    if (it != temporary_terms.end())
      return 0;

    int cost = arg1->cost(temporary_terms, is_matlab);
    cost += arg2->cost(temporary_terms, is_matlab);

    if (is_matlab)
      // Cost for Matlab files
      switch (op_code)
        {
        case oLess:
        case oGreater:
        case oLessEqual:
        case oGreaterEqual:
        case oEqualEqual:
        case oDifferent:
          return cost + 60;
        case oPlus:
        case oMinus:
        case oTimes:
          return cost + 90;
        case oMax:
        case oMin:
          return cost + 110;
        case oDivide:
          return cost + 990;
        case oPower:
          return cost + 1160;
        case oEqual:
          return cost;
        }
    else
      // Cost for C files
      switch (op_code)
        {
        case oLess:
        case oGreater:
        case oLessEqual:
        case oGreaterEqual:
        case oEqualEqual:
        case oDifferent:
          return cost + 2;
        case oPlus:
        case oMinus:
        case oTimes:
          return cost + 4;
        case oMax:
        case oMin:
          return cost + 5;
        case oDivide:
          return cost + 15;
        case oPower:
          return cost + 520;
        case oEqual:
          return cost;
        }
    // Suppress GCC warning
    exit(EXIT_FAILURE);
  }

void
BinaryOpNode::computeTemporaryTerms(map<NodeID, int> &reference_count,
                                    temporary_terms_type &temporary_terms,
                                    bool is_matlab) const
  {
    NodeID this2 = const_cast<BinaryOpNode *>(this);
    map<NodeID, int>::iterator it = reference_count.find(this2);
    if (it == reference_count.end())
      {
        // If this node has never been encountered, set its ref count to one,
        //  and travel through its children
        reference_count[this2] = 1;
        arg1->computeTemporaryTerms(reference_count, temporary_terms, is_matlab);
        arg2->computeTemporaryTerms(reference_count, temporary_terms, is_matlab);
      }
    else
      {
        // If the node has already been encountered, increment its ref count
        //  and declare it as a temporary term if it is too costly
        reference_count[this2]++;
        if (reference_count[this2] * cost(temporary_terms, is_matlab) > MIN_COST(is_matlab))
          temporary_terms.insert(this2);
      }
  }

void
BinaryOpNode::computeTemporaryTerms(map<NodeID, int> &reference_count,
                                    temporary_terms_type &temporary_terms,
                                    map<NodeID, pair<int, int> > &first_occurence,
                                    int Curr_block,
                                    Model_Block *ModelBlock,
                                    int equation,
                                    map_idx_type &map_idx) const
  {
    NodeID this2 = const_cast<BinaryOpNode *>(this);
    map<NodeID, int>::iterator it = reference_count.find(this2);
    if (it == reference_count.end())
      {
        reference_count[this2] = 1;
        first_occurence[this2] = make_pair(Curr_block, equation);
        arg1->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, Curr_block, ModelBlock, equation, map_idx);
        arg2->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, Curr_block, ModelBlock, equation, map_idx);
      }
    else
      {
        reference_count[this2]++;
        if (reference_count[this2] * cost(temporary_terms, false) > MIN_COST_C)
          {
            temporary_terms.insert(this2);
            ModelBlock->Block_List[first_occurence[this2].first].Temporary_Terms_in_Equation[first_occurence[this2].second]->insert(this2);
          }
      }
  }

double
BinaryOpNode::eval_opcode(double v1, BinaryOpcode op_code, double v2) throw (EvalException)
{
  switch (op_code)
    {
    case oPlus:
      return(v1 + v2);
    case oMinus:
      return(v1 - v2);
    case oTimes:
      return(v1 * v2);
    case oDivide:
      return(v1 / v2);
    case oPower:
      return(pow(v1, v2));
    case oMax:
      if (v1 < v2)
        return v2;
      else
        return v1;
    case oMin:
      if (v1 > v2)
        return v2;
      else
        return v1;
    case oLess:
      return (v1 < v2);
    case oGreater:
      return (v1 > v2);
    case oLessEqual:
      return (v1 <= v2);
    case oGreaterEqual:
      return (v1 >= v2);
    case oEqualEqual:
      return (v1 == v2);
    case oDifferent:
      return (v1 != v2);
    case oEqual:
      throw EvalException();
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

double
BinaryOpNode::eval(const eval_context_type &eval_context) const throw (EvalException)
{
  double v1 = arg1->eval(eval_context);
  double v2 = arg2->eval(eval_context);

  return eval_opcode(v1, op_code, v2);
}

void
BinaryOpNode::compile(ostream &CompileCode, bool lhs_rhs, const temporary_terms_type &temporary_terms, map_idx_type &map_idx, bool dynamic, bool steady_dynamic) const
  {
    // If current node is a temporary term
    temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<BinaryOpNode *>(this));
    if (it != temporary_terms.end())
      {
        if(dynamic)
          CompileCode.write(&FLDT, sizeof(FLDT));
        else
          CompileCode.write(&FLDST, sizeof(FLDST));
        int var=map_idx[idx];
        CompileCode.write(reinterpret_cast<char *>(&var), sizeof(var));
        return;
      }
    arg1->compile(CompileCode, lhs_rhs, temporary_terms, map_idx, dynamic, steady_dynamic);
    arg2->compile(CompileCode, lhs_rhs, temporary_terms, map_idx, dynamic, steady_dynamic);
    CompileCode.write(&FBINARY, sizeof(FBINARY));
    BinaryOpcode op_codel=op_code;
    CompileCode.write(reinterpret_cast<char *>(&op_codel),sizeof(op_codel));
  }

void
BinaryOpNode::collectTemporary_terms(const temporary_terms_type &temporary_terms, Model_Block *ModelBlock, int Curr_Block) const
  {
    temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<BinaryOpNode *>(this));
    if (it != temporary_terms.end())
      ModelBlock->Block_List[Curr_Block].Temporary_InUse->insert(idx);
    else
      {
        arg1->collectTemporary_terms(temporary_terms, ModelBlock, Curr_Block);
        arg2->collectTemporary_terms(temporary_terms, ModelBlock, Curr_Block);
      }
  }


void
BinaryOpNode::writeOutput(ostream &output, ExprNodeOutputType output_type,
                          const temporary_terms_type &temporary_terms) const
  {
    //cout << "writeOutput binary\n";
    // If current node is a temporary term
    temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<BinaryOpNode *>(this));
    if (it != temporary_terms.end())
      {
        if (output_type == oMatlabDynamicModelSparse)
          output << "T" << idx << "(it_)";
        else
          output << "T" << idx;
        return;
      }

    // Treat special case of power operator in C, and case of max and min operators
    if ((op_code == oPower && IS_C(output_type)) || op_code == oMax || op_code == oMin )
      {
        switch (op_code)
          {
          case oPower:
            output << "pow(";
            break;
          case oMax:
            output << "max(";
            break;
          case oMin:
            output << "min(";
            break;
          default:
            ;
          }
        arg1->writeOutput(output, output_type, temporary_terms);
        output << ",";
        arg2->writeOutput(output, output_type, temporary_terms);
        output << ")";
        return;
      }

    int prec = precedence(output_type, temporary_terms);

    bool close_parenthesis = false;

    if (IS_LATEX(output_type) && op_code == oDivide)
      output << "\\frac{";
    else
      {
        // If left argument has a lower precedence, or if current and left argument are both power operators, add parenthesis around left argument
        BinaryOpNode *barg1 = dynamic_cast<BinaryOpNode *>(arg1);
        if (arg1->precedence(output_type, temporary_terms) < prec
            || (op_code == oPower && barg1 != NULL && barg1->op_code == oPower))
          {
            output << LEFT_PAR(output_type);
            close_parenthesis = true;
          }
      }

    // Write left argument
    arg1->writeOutput(output, output_type, temporary_terms);

    if (close_parenthesis)
      output << RIGHT_PAR(output_type);

    if (IS_LATEX(output_type) && op_code == oDivide)
      output << "}";


    // Write current operator symbol
    switch (op_code)
      {
      case oPlus:
        output << "+";
        break;
      case oMinus:
        output << "-";
        break;
      case oTimes:
        if (IS_LATEX(output_type))
          output << "\\, ";
        else
          output << "*";
        break;
      case oDivide:
        if (!IS_LATEX(output_type))
          output << "/";
        break;
      case oPower:
        output << "^";
        break;
      case oLess:
        output << "<";
        break;
      case oGreater:
        output << ">";
        break;
      case oLessEqual:
        if (IS_LATEX(output_type))
          output << "\\leq ";
        else
          output << "<=";
        break;
      case oGreaterEqual:
        if (IS_LATEX(output_type))
          output << "\\geq ";
        else
          output << ">=";
        break;
      case oEqualEqual:
        output << "==";
        break;
      case oDifferent:
        if (IS_MATLAB(output_type))
          output << "~=";
        else
          {
            if (IS_C(output_type))
              output << "!=";
            else
              output << "\\neq ";
          }
        break;
      case oEqual:
        output << "=";
        break;
      default:
        ;
      }

    close_parenthesis = false;

    if (IS_LATEX(output_type) && (op_code == oPower || op_code == oDivide))
      output << "{";
    else
      {
        /* Add parenthesis around right argument if:
           - its precedence is lower than those of the current node
           - it is a power operator and current operator is also a power operator
           - it is a minus operator with same precedence than current operator
           - it is a divide operator with same precedence than current operator */
        BinaryOpNode *barg2 = dynamic_cast<BinaryOpNode *>(arg2);
        int arg2_prec = arg2->precedence(output_type, temporary_terms);
        if (arg2_prec < prec
            || (op_code == oPower && barg2 != NULL && barg2->op_code == oPower && !IS_LATEX(output_type))
            || (op_code == oMinus && arg2_prec == prec)
            || (op_code == oDivide && arg2_prec == prec && !IS_LATEX(output_type)))
          {
            output << LEFT_PAR(output_type);
            close_parenthesis = true;
          }
      }

    // Write right argument
    arg2->writeOutput(output, output_type, temporary_terms);

    if (IS_LATEX(output_type) && (op_code == oPower || op_code == oDivide))
      output << "}";

    if (close_parenthesis)
      output << RIGHT_PAR(output_type);
  }

void
BinaryOpNode::collectVariables(SymbolType type_arg, set<pair<int, int> > &result) const
{
  arg1->collectVariables(type_arg, result);
  arg2->collectVariables(type_arg, result);
}

NodeID
BinaryOpNode::Compute_RHS(NodeID arg1, NodeID arg2, int op, int op_type) const
{
  temporary_terms_type temp;
  switch(op_type)
    {
    case 0: /*Unary Operator*/
      switch(op)
        {
        case oUminus:
          return(datatree.AddUMinus(arg1));
          break;
        case oExp:
          return(datatree.AddExp(arg1));
          break;
        case oLog:
          return(datatree.AddLog(arg1));
          break;
        case oLog10:
          return(datatree.AddLog10(arg1));
          break;
        }
      break;
    case 1: /*Binary Operator*/
      switch(op)
        {
        case oPlus:
          return(datatree.AddPlus(arg1, arg2));
          break;
        case oMinus:
          return(datatree.AddMinus(arg1, arg2));
          break;
        case oTimes:
          return(datatree.AddTimes(arg1, arg2));
          break;
        case oDivide:
          return(datatree.AddDivide(arg1, arg2));
          break;
        case oPower:
          return(datatree.AddPower(arg1, arg2));
          break;
        }
      break;
    }
  return((NodeID)NULL);
}

pair<int, NodeID>
BinaryOpNode::normalizeEquation(int var_endo, vector<pair<int, pair<NodeID, NodeID> > > &List_of_Op_RHS) const
  {
    vector<pair<int, pair<NodeID, NodeID> > > List_of_Op_RHS1, List_of_Op_RHS2;
    int is_endogenous_present_1, is_endogenous_present_2;
    pair<int, NodeID> res;
    NodeID NodeID_1, NodeID_2;
    res = arg1->normalizeEquation(var_endo, List_of_Op_RHS1);
    is_endogenous_present_1 = res.first;
    NodeID_1 = res.second;

    res = arg2->normalizeEquation(var_endo, List_of_Op_RHS2);
    is_endogenous_present_2 = res.first;
    NodeID_2 = res.second;
    if(is_endogenous_present_1==2 || is_endogenous_present_2==2)
      return(make_pair(2,(NodeID)NULL));
    else if(is_endogenous_present_1 && is_endogenous_present_2)
      return(make_pair(2,(NodeID)NULL));
    else if(is_endogenous_present_1)
      {
        if(op_code==oEqual)
          {
            pair<int, pair<NodeID, NodeID> > it;
            int oo=List_of_Op_RHS1.size();
            for(int i=0;i<oo;i++)
              {
                it = List_of_Op_RHS1.back();
                List_of_Op_RHS1.pop_back();
                if(it.second.first && !it.second.second) /*Binary operator*/
                  NodeID_2 = Compute_RHS(NodeID_2, (BinaryOpNode*)it.second.first, it.first, 1);
                else if(it.second.second && !it.second.first) /*Binary operator*/
                  NodeID_2 = Compute_RHS(it.second.second, NodeID_2, it.first, 1);
                else if(it.second.second && it.second.first) /*Binary operator*/
                  NodeID_2 = Compute_RHS(it.second.first, it.second.second, it.first, 1);
                else  /*Unary operator*/
                  NodeID_2 = Compute_RHS((UnaryOpNode*)NodeID_2, (UnaryOpNode*)it.second.first, it.first, 0);
              }
          }
        else
          List_of_Op_RHS = List_of_Op_RHS1;
      }
    else if(is_endogenous_present_2)
      {
        if(op_code==oEqual)
          {
            int oo=List_of_Op_RHS2.size();
            for(int i=0;i<oo;i++)
              {
                pair<int, pair<NodeID, NodeID> > it;
                it = List_of_Op_RHS2.back();
                List_of_Op_RHS2.pop_back();
                if(it.second.first && !it.second.second) /*Binary operator*/
                  NodeID_1 = Compute_RHS((BinaryOpNode*)NodeID_1, (BinaryOpNode*)it.second.first, it.first, 1);
                else if(it.second.second && !it.second.first) /*Binary operator*/
                  NodeID_1 = Compute_RHS((BinaryOpNode*)it.second.second, (BinaryOpNode*)NodeID_1, it.first, 1);
                else if(it.second.second && it.second.first) /*Binary operator*/
                  NodeID_1 = Compute_RHS(it.second.first, it.second.second, it.first, 1);
                else
                  NodeID_1 = Compute_RHS((UnaryOpNode*)NodeID_1, (UnaryOpNode*)it.second.first, it.first, 0);
              }
          }
        else
          List_of_Op_RHS =List_of_Op_RHS2;
      }
    switch (op_code)
      {
      case oPlus:
        if (!is_endogenous_present_1 && !is_endogenous_present_2)
          {
            List_of_Op_RHS.push_back(make_pair(oMinus, make_pair(datatree.AddPlus(NodeID_1, NodeID_2), (NodeID)NULL)));
            return(make_pair(0, datatree.AddPlus(NodeID_1, NodeID_2)));
          }
        else if (is_endogenous_present_1 && is_endogenous_present_2)
          return(make_pair(1, (NodeID)NULL));
        else if (!is_endogenous_present_1 && is_endogenous_present_2)
          {
            List_of_Op_RHS.push_back(make_pair(oMinus, make_pair(NodeID_1, (NodeID)NULL)));
            return(make_pair(1, NodeID_1));
          }
        else if (is_endogenous_present_1 && !is_endogenous_present_2)
          {
            List_of_Op_RHS.push_back(make_pair(oMinus, make_pair(NodeID_2, (NodeID)NULL) ));
            return(make_pair(1, NodeID_2));
          }
        break;
      case oMinus:
        if (!is_endogenous_present_1 && !is_endogenous_present_2)
          {
            List_of_Op_RHS.push_back(make_pair(oMinus, make_pair(datatree.AddMinus(NodeID_1, NodeID_2), (NodeID)NULL) ));
            return(make_pair(0, datatree.AddMinus(NodeID_1, NodeID_2)));
          }
        else if (is_endogenous_present_1 && is_endogenous_present_2)
          return(make_pair(1, (NodeID)NULL));
        else if (!is_endogenous_present_1 && is_endogenous_present_2)
          {
            List_of_Op_RHS.push_back(make_pair(oUminus, make_pair((NodeID)NULL, (NodeID)NULL)));
            List_of_Op_RHS.push_back(make_pair(oMinus, make_pair(NodeID_1, (NodeID)NULL) ));
            return(make_pair(1, NodeID_1));
          }
        else if (is_endogenous_present_1 && !is_endogenous_present_2)
          {
            List_of_Op_RHS.push_back(make_pair(oPlus, make_pair(NodeID_2, (NodeID) NULL) ));
            return(make_pair(1, datatree.AddUMinus(NodeID_2)));
          }
        break;
      case oTimes:
        if (!is_endogenous_present_1 && !is_endogenous_present_2)
          return(make_pair(0, datatree.AddTimes(NodeID_1, NodeID_2)));
        else if(!is_endogenous_present_1 && is_endogenous_present_2)
          {
            List_of_Op_RHS.push_back(make_pair(oDivide, make_pair(NodeID_1, (NodeID)NULL) ));
            return(make_pair(1, NodeID_1));
          }
        else if(is_endogenous_present_1 && !is_endogenous_present_2)
          {
            List_of_Op_RHS.push_back(make_pair(oDivide, make_pair(NodeID_2, (NodeID)NULL) ));
            return(make_pair(1, NodeID_2));
          }
        else
          return(make_pair(1, (NodeID)NULL));
        break;
      case oDivide:
        if (!is_endogenous_present_1 && !is_endogenous_present_2)
          return(make_pair(0, datatree.AddDivide(NodeID_1, NodeID_2)));
        else if(!is_endogenous_present_1 && is_endogenous_present_2)
          {
            List_of_Op_RHS.push_back(make_pair(oDivide, make_pair((NodeID)NULL, NodeID_1) ));
            return(make_pair(1, NodeID_1));
          }
        else if(is_endogenous_present_1 && !is_endogenous_present_2)
          {
            List_of_Op_RHS.push_back(make_pair(oTimes, make_pair(NodeID_2, (NodeID)NULL) ));
            return(make_pair(1, NodeID_2));
          }
        else
          return(make_pair(1, (NodeID)NULL));
        break;
      case oPower:
        if (!is_endogenous_present_1 && !is_endogenous_present_2)
          return(make_pair(0, datatree.AddPower(NodeID_1, NodeID_2)));
        else if(is_endogenous_present_1 && !is_endogenous_present_2)
          {
            List_of_Op_RHS.push_back(make_pair(oPower, make_pair(datatree.AddDivide( datatree.AddNumConstant("1"), NodeID_2), (NodeID)NULL) ));
            return(make_pair(1, (NodeID)NULL));
          }
        break;
      case oEqual:
        if (!is_endogenous_present_1 && !is_endogenous_present_2)
          {
              return( make_pair(0,
              datatree.AddEqual(datatree.AddVariable(datatree.symbol_table.getName(datatree.symbol_table.getID(eEndogenous, var_endo)), 0), datatree.AddMinus(NodeID_2, NodeID_1))
              ));
          }
        else if (is_endogenous_present_1 && is_endogenous_present_2)
          {
            return(make_pair(0,
            datatree.AddEqual(datatree.AddVariable(datatree.symbol_table.getName(datatree.symbol_table.getID(eEndogenous, var_endo)), 0), datatree.Zero)
            ));
          }
        else if (!is_endogenous_present_1 && is_endogenous_present_2)
          {
              return(make_pair(0,
              datatree.AddEqual(datatree.AddVariable(datatree.symbol_table.getName(datatree.symbol_table.getID(eEndogenous, var_endo)), 0), /*datatree.AddUMinus(NodeID_1)*/NodeID_1)
              ));
          }
        else if (is_endogenous_present_1 && !is_endogenous_present_2)
          {
              return(make_pair(0,
              datatree.AddEqual(datatree.AddVariable(datatree.symbol_table.getName(datatree.symbol_table.getID(eEndogenous, var_endo)), 0), NodeID_2)
              ));
          }
        break;
      case oMax:
        if (!is_endogenous_present_1 && !is_endogenous_present_2)
          return(make_pair(0, datatree.AddMax(NodeID_1, NodeID_2) ));
        else
          return(make_pair(1, (NodeID)NULL ));
        break;
      case oMin:
        if (!is_endogenous_present_1 && !is_endogenous_present_2)
          return(make_pair(0, datatree.AddMin(NodeID_1, NodeID_2) ));
        else
          return(make_pair(1, (NodeID)NULL ));
        break;
      case oLess:
        if (!is_endogenous_present_1 && !is_endogenous_present_2)
          return(make_pair(0, datatree.AddLess(NodeID_1, NodeID_2) ));
        else
          return(make_pair(1, (NodeID)NULL ));
        break;
      case oGreater:
        if (!is_endogenous_present_1 && !is_endogenous_present_2)
          return(make_pair(0, datatree.AddGreater(NodeID_1, NodeID_2) ));
        else
          return(make_pair(1, (NodeID)NULL ));
        break;
      case oLessEqual:
        if (!is_endogenous_present_1 && !is_endogenous_present_2)
          return(make_pair(0, datatree.AddLessEqual(NodeID_1, NodeID_2) ));
        else
          return(make_pair(1, (NodeID)NULL ));
        break;
      case oGreaterEqual:
        if (!is_endogenous_present_1 && !is_endogenous_present_2)
          return(make_pair(0, datatree.AddGreaterEqual(NodeID_1, NodeID_2) ));
        else
          return(make_pair(1, (NodeID)NULL ));
        break;
      case oEqualEqual:
        if (!is_endogenous_present_1 && !is_endogenous_present_2)
          return(make_pair(0, datatree.AddEqualEqual(NodeID_1, NodeID_2) ));
        else
          return(make_pair(1, (NodeID)NULL ));
        break;
      case oDifferent:
        if (!is_endogenous_present_1 && !is_endogenous_present_2)
          return(make_pair(0, datatree.AddDifferent(NodeID_1, NodeID_2) ));
        else
          return(make_pair(1, (NodeID)NULL ));
        break;
      }
    // Suppress GCC warning
    exit(EXIT_FAILURE);
  }


NodeID
BinaryOpNode::getChainRuleDerivative(int deriv_id, const map<int, NodeID> &recursive_variables)
{
  NodeID darg1 = arg1->getChainRuleDerivative(deriv_id, recursive_variables);
  NodeID darg2 = arg2->getChainRuleDerivative(deriv_id, recursive_variables);
  return composeDerivatives(darg1, darg2);
}

NodeID
BinaryOpNode::toStatic(DataTree &static_datatree) const
  {
    NodeID sarg1 = arg1->toStatic(static_datatree);
    NodeID sarg2 = arg2->toStatic(static_datatree);
    switch (op_code)
      {
      case oPlus:
        return static_datatree.AddPlus(sarg1, sarg2);
      case oMinus:
        return static_datatree.AddMinus(sarg1, sarg2);
      case oTimes:
        return static_datatree.AddTimes(sarg1, sarg2);
      case oDivide:
        return static_datatree.AddDivide(sarg1, sarg2);
      case oPower:
        return static_datatree.AddPower(sarg1, sarg2);
      case oEqual:
        return static_datatree.AddEqual(sarg1, sarg2);
      case oMax:
        return static_datatree.AddMax(sarg1, sarg2);
      case oMin:
        return static_datatree.AddMin(sarg1, sarg2);
      case oLess:
        return static_datatree.AddLess(sarg1, sarg2);
      case oGreater:
        return static_datatree.AddGreater(sarg1, sarg2);
      case oLessEqual:
        return static_datatree.AddLessEqual(sarg1, sarg2);
      case oGreaterEqual:
        return static_datatree.AddGreaterEqual(sarg1, sarg2);
      case oEqualEqual:
        return static_datatree.AddEqualEqual(sarg1, sarg2);
      case oDifferent:
        return static_datatree.AddDifferent(sarg1, sarg2);
      }
    // Suppress GCC warning
    exit(EXIT_FAILURE);
  }


TrinaryOpNode::TrinaryOpNode(DataTree &datatree_arg, const NodeID arg1_arg,
                             TrinaryOpcode op_code_arg, const NodeID arg2_arg, const NodeID arg3_arg) :
    ExprNode(datatree_arg),
    arg1(arg1_arg),
    arg2(arg2_arg),
    arg3(arg3_arg),
    op_code(op_code_arg)
{
  datatree.trinary_op_node_map[make_pair(make_pair(make_pair(arg1, arg2), arg3), op_code)] = this;

  // Non-null derivatives are the union of those of the arguments
  // Compute set union of arg{1,2,3}->non_null_derivatives
  set<int> non_null_derivatives_tmp;
  set_union(arg1->non_null_derivatives.begin(),
            arg1->non_null_derivatives.end(),
            arg2->non_null_derivatives.begin(),
            arg2->non_null_derivatives.end(),
            inserter(non_null_derivatives_tmp, non_null_derivatives_tmp.begin()));
  set_union(non_null_derivatives_tmp.begin(),
            non_null_derivatives_tmp.end(),
            arg3->non_null_derivatives.begin(),
            arg3->non_null_derivatives.end(),
            inserter(non_null_derivatives, non_null_derivatives.begin()));
}

NodeID
TrinaryOpNode::composeDerivatives(NodeID darg1, NodeID darg2, NodeID darg3)
{

  NodeID t11, t12, t13, t14, t15;

  switch (op_code)
    {
    case oNormcdf:
      // normal pdf is inlined in the tree
      NodeID y;
      // sqrt(2*pi)
      t14 = datatree.AddSqrt(datatree.AddTimes(datatree.Two, datatree.Pi));
      // x - mu
      t12 = datatree.AddMinus(arg1,arg2);
      // y = (x-mu)/sigma
      y = datatree.AddDivide(t12,arg3);
      // (x-mu)^2/sigma^2
      t12 = datatree.AddTimes(y,y);
      // -(x-mu)^2/sigma^2
      t13 = datatree.AddUMinus(t12);
      // -((x-mu)^2/sigma^2)/2
      t12 = datatree.AddDivide(t13, datatree.Two);
      // exp(-((x-mu)^2/sigma^2)/2)
      t13 = datatree.AddExp(t12);
      // derivative of a standardized normal
      // t15 = (1/sqrt(2*pi))*exp(-y^2/2)
      t15 = datatree.AddDivide(t13,t14);
      // derivatives thru x
      t11 = datatree.AddDivide(darg1,arg3);
      // derivatives thru mu
      t12 = datatree.AddDivide(darg2,arg3);
      // intermediary sum
      t14 = datatree.AddMinus(t11,t12);
      // derivatives thru sigma
      t11 = datatree.AddDivide(y,arg3);
      t12 = datatree.AddTimes(t11,darg3);
      //intermediary sum
      t11 = datatree.AddMinus(t14,t12);
      // total derivative:
      // (darg1/sigma - darg2/sigma - darg3*(x-mu)/sigma^2) * t15
      // where t15 is the derivative of a standardized normal
      return datatree.AddTimes(t11, t15);
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

NodeID
TrinaryOpNode::computeDerivative(int deriv_id)
{
  NodeID darg1 = arg1->getDerivative(deriv_id);
  NodeID darg2 = arg2->getDerivative(deriv_id);
  NodeID darg3 = arg3->getDerivative(deriv_id);
  return composeDerivatives(darg1, darg2, darg3);
}

int
TrinaryOpNode::precedence(ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms) const
  {
    temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<TrinaryOpNode *>(this));
    // A temporary term behaves as a variable
    if (it != temporary_terms.end())
      return 100;

    switch (op_code)
      {
      case oNormcdf:
        return 100;
      }
    // Suppress GCC warning
    exit(EXIT_FAILURE);
  }

int
TrinaryOpNode::cost(const temporary_terms_type &temporary_terms, bool is_matlab) const
  {
    temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<TrinaryOpNode *>(this));
    // For a temporary term, the cost is null
    if (it != temporary_terms.end())
      return 0;

    int cost = arg1->cost(temporary_terms, is_matlab);
    cost += arg2->cost(temporary_terms, is_matlab);

    if (is_matlab)
      // Cost for Matlab files
      switch (op_code)
        {
        case oNormcdf:
          return cost+1000;
        }
    else
      // Cost for C files
      switch (op_code)
        {
        case oNormcdf:
          return cost+1000;
        }
    // Suppress GCC warning
    exit(EXIT_FAILURE);
  }

void
TrinaryOpNode::computeTemporaryTerms(map<NodeID, int> &reference_count,
                                     temporary_terms_type &temporary_terms,
                                     bool is_matlab) const
  {
    NodeID this2 = const_cast<TrinaryOpNode *>(this);
    map<NodeID, int>::iterator it = reference_count.find(this2);
    if (it == reference_count.end())
      {
        // If this node has never been encountered, set its ref count to one,
        //  and travel through its children
        reference_count[this2] = 1;
        arg1->computeTemporaryTerms(reference_count, temporary_terms, is_matlab);
        arg2->computeTemporaryTerms(reference_count, temporary_terms, is_matlab);
        arg3->computeTemporaryTerms(reference_count, temporary_terms, is_matlab);
      }
    else
      {
        // If the node has already been encountered, increment its ref count
        //  and declare it as a temporary term if it is too costly
        reference_count[this2]++;
        if (reference_count[this2] * cost(temporary_terms, is_matlab) > MIN_COST(is_matlab))
          temporary_terms.insert(this2);
      }
  }

void
TrinaryOpNode::computeTemporaryTerms(map<NodeID, int> &reference_count,
                                     temporary_terms_type &temporary_terms,
                                     map<NodeID, pair<int, int> > &first_occurence,
                                     int Curr_block,
                                     Model_Block *ModelBlock,
                                     int equation,
                                     map_idx_type &map_idx) const
  {
    NodeID this2 = const_cast<TrinaryOpNode *>(this);
    map<NodeID, int>::iterator it = reference_count.find(this2);
    if (it == reference_count.end())
      {
        reference_count[this2] = 1;
        first_occurence[this2] = make_pair(Curr_block,equation);
        arg1->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, Curr_block, ModelBlock, equation, map_idx);
        arg2->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, Curr_block, ModelBlock, equation, map_idx);
        arg3->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, Curr_block, ModelBlock, equation, map_idx);
      }
    else
      {
        reference_count[this2]++;
        if (reference_count[this2] * cost(temporary_terms, false) > MIN_COST_C)
          {
            temporary_terms.insert(this2);
            ModelBlock->Block_List[first_occurence[this2].first].Temporary_Terms_in_Equation[first_occurence[this2].second]->insert(this2);
          }
      }
  }

double
TrinaryOpNode::eval_opcode(double v1, TrinaryOpcode op_code, double v2, double v3) throw (EvalException)
{
  switch (op_code)
    {
    case oNormcdf:
      cerr << "NORMCDF: eval not implemented" << endl;
      exit(EXIT_FAILURE);
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

double
TrinaryOpNode::eval(const eval_context_type &eval_context) const throw (EvalException)
{
  double v1 = arg1->eval(eval_context);
  double v2 = arg2->eval(eval_context);
  double v3 = arg3->eval(eval_context);

  return eval_opcode(v1, op_code, v2, v3);
}

void
TrinaryOpNode::compile(ostream &CompileCode, bool lhs_rhs, const temporary_terms_type &temporary_terms, map_idx_type &map_idx, bool dynamic, bool steady_dynamic) const
  {
    // If current node is a temporary term
    temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<TrinaryOpNode *>(this));
    if (it != temporary_terms.end())
      {
        if(dynamic)
          CompileCode.write(&FLDT, sizeof(FLDT));
        else
          CompileCode.write(&FLDST, sizeof(FLDST));
        int var=map_idx[idx];
        CompileCode.write(reinterpret_cast<char *>(&var), sizeof(var));
        return;
      }
    arg1->compile(CompileCode, lhs_rhs, temporary_terms, map_idx, dynamic, steady_dynamic);
    arg2->compile(CompileCode, lhs_rhs, temporary_terms, map_idx, dynamic, steady_dynamic);
    arg3->compile(CompileCode, lhs_rhs, temporary_terms, map_idx, dynamic, steady_dynamic);
    CompileCode.write(&FBINARY, sizeof(FBINARY));
    TrinaryOpcode op_codel=op_code;
    CompileCode.write(reinterpret_cast<char *>(&op_codel),sizeof(op_codel));
  }

void
TrinaryOpNode::collectTemporary_terms(const temporary_terms_type &temporary_terms, Model_Block *ModelBlock, int Curr_Block) const
  {
    temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<TrinaryOpNode *>(this));
    if (it != temporary_terms.end())
      ModelBlock->Block_List[Curr_Block].Temporary_InUse->insert(idx);
    else
      {
        arg1->collectTemporary_terms(temporary_terms, ModelBlock, Curr_Block);
        arg2->collectTemporary_terms(temporary_terms, ModelBlock, Curr_Block);
        arg3->collectTemporary_terms(temporary_terms, ModelBlock, Curr_Block);
      }
  }


void
TrinaryOpNode::writeOutput(ostream &output, ExprNodeOutputType output_type,
                           const temporary_terms_type &temporary_terms) const
  {
    // TrinaryOpNode not implemented for C output
    assert(!IS_C(output_type));

    // If current node is a temporary term
    temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<TrinaryOpNode *>(this));
    if (it != temporary_terms.end())
      {
        output << "T" << idx;
        return;
      }

    switch (op_code)
      {
      case oNormcdf:
        output << "normcdf(";
        break;
      }
    arg1->writeOutput(output, output_type, temporary_terms);
    output << ",";
    arg2->writeOutput(output, output_type, temporary_terms);
    output << ",";
    arg3->writeOutput(output, output_type, temporary_terms);
    output << ")";
  }

void
TrinaryOpNode::collectVariables(SymbolType type_arg, set<pair<int, int> > &result) const
{
  arg1->collectVariables(type_arg, result);
  arg2->collectVariables(type_arg, result);
  arg3->collectVariables(type_arg, result);
}

pair<int, NodeID>
TrinaryOpNode::normalizeEquation(int var_endo, vector<pair<int, pair<NodeID, NodeID> > > &List_of_Op_RHS) const
  {
    pair<int, NodeID> res = arg1->normalizeEquation(var_endo, List_of_Op_RHS);
    bool is_endogenous_present_1 = res.first;
    NodeID NodeID_1 = res.second;
    res = arg2->normalizeEquation(var_endo, List_of_Op_RHS);
    bool is_endogenous_present_2 = res.first;
    NodeID NodeID_2 = res.second;
    res = arg3->normalizeEquation(var_endo, List_of_Op_RHS);
    bool is_endogenous_present_3 = res.first;
    NodeID NodeID_3 = res.second;
    if (!is_endogenous_present_1 && !is_endogenous_present_2 && !is_endogenous_present_3)
      return(make_pair(0, datatree.AddNormcdf(NodeID_1, NodeID_2, NodeID_3) ));
    else
      return(make_pair(1, (NodeID)NULL ));
  }

NodeID
TrinaryOpNode::getChainRuleDerivative(int deriv_id, const map<int, NodeID> &recursive_variables)
{
  NodeID darg1 = arg1->getChainRuleDerivative(deriv_id, recursive_variables);
  NodeID darg2 = arg2->getChainRuleDerivative(deriv_id, recursive_variables);
  NodeID darg3 = arg3->getChainRuleDerivative(deriv_id, recursive_variables);
  return composeDerivatives(darg1, darg2, darg3);
}

NodeID
TrinaryOpNode::toStatic(DataTree &static_datatree) const
  {
    NodeID sarg1 = arg1->toStatic(static_datatree);
    NodeID sarg2 = arg2->toStatic(static_datatree);
    NodeID sarg3 = arg3->toStatic(static_datatree);
    switch (op_code)
      {
      case oNormcdf:
        return static_datatree.AddNormcdf(sarg1, sarg2, sarg3);
      }
    // Suppress GCC warning
    exit(EXIT_FAILURE);
  }


UnknownFunctionNode::UnknownFunctionNode(DataTree &datatree_arg,
    int symb_id_arg,
    const vector<NodeID> &arguments_arg) :
    ExprNode(datatree_arg),
    symb_id(symb_id_arg),
    arguments(arguments_arg)
{
}

NodeID
UnknownFunctionNode::computeDerivative(int deriv_id)
{
  cerr << "UnknownFunctionNode::computeDerivative: operation impossible!" << endl;
  exit(EXIT_FAILURE);
}

NodeID
UnknownFunctionNode::getChainRuleDerivative(int deriv_id, const map<int, NodeID> &recursive_variables)
{
  cerr << "UnknownFunctionNode::getChainRuleDerivative: operation impossible!" << endl;
  exit(EXIT_FAILURE);
}


void
UnknownFunctionNode::computeTemporaryTerms(map<NodeID, int> &reference_count,
    temporary_terms_type &temporary_terms,
    bool is_matlab) const
  {
    cerr << "UnknownFunctionNode::computeTemporaryTerms: operation impossible!" << endl;
    exit(EXIT_FAILURE);
  }

void UnknownFunctionNode::writeOutput(ostream &output, ExprNodeOutputType output_type,
                                      const temporary_terms_type &temporary_terms) const
  {
    output << datatree.symbol_table.getName(symb_id) << "(";
    for (vector<NodeID>::const_iterator it = arguments.begin();
         it != arguments.end(); it++)
      {
        if (it != arguments.begin())
          output << ",";

        (*it)->writeOutput(output, output_type, temporary_terms);
      }
    output << ")";
  }

void
UnknownFunctionNode::computeTemporaryTerms(map<NodeID, int> &reference_count,
    temporary_terms_type &temporary_terms,
    map<NodeID, pair<int, int> > &first_occurence,
    int Curr_block,
    Model_Block *ModelBlock,
    int equation,
    map_idx_type &map_idx) const
  {
    cerr << "UnknownFunctionNode::computeTemporaryTerms: not implemented" << endl;
    exit(EXIT_FAILURE);
  }

void
UnknownFunctionNode::collectVariables(SymbolType type_arg, set<pair<int, int> > &result) const
{
  for (vector<NodeID>::const_iterator it = arguments.begin();
       it != arguments.end(); it++)
    (*it)->collectVariables(type_arg, result);
}

void
UnknownFunctionNode::collectTemporary_terms(const temporary_terms_type &temporary_terms, Model_Block *ModelBlock, int Curr_Block) const
  {
    temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<UnknownFunctionNode *>(this));
    if (it != temporary_terms.end())
      ModelBlock->Block_List[Curr_Block].Temporary_InUse->insert(idx);
    else
      {
        //arg->collectTemporary_terms(temporary_terms, result);
      }
  }


double
UnknownFunctionNode::eval(const eval_context_type &eval_context) const throw (EvalException)
{
  throw EvalException();
}

void
UnknownFunctionNode::compile(ostream &CompileCode, bool lhs_rhs, const temporary_terms_type &temporary_terms, map_idx_type &map_idx, bool dynamic, bool steady_dynamic) const
  {
    cerr << "UnknownFunctionNode::compile: operation impossible!" << endl;
    exit(EXIT_FAILURE);
  }

pair<int, NodeID>
UnknownFunctionNode::normalizeEquation(int var_endo, vector<pair<int, pair<NodeID, NodeID> > >  &List_of_Op_RHS) const
  {
    vector<pair<bool, NodeID> > V_arguments;
    vector<NodeID> V_NodeID;
    bool present = false;
    for (vector<NodeID>::const_iterator it = arguments.begin();
         it != arguments.end(); it++)
      {
        V_arguments.push_back((*it)->normalizeEquation(var_endo, List_of_Op_RHS));
        present = present || V_arguments[V_arguments.size()-1].first;
        V_NodeID.push_back(V_arguments[V_arguments.size()-1].second);
      }
    if (!present)
      return(make_pair(0, datatree.AddUnknownFunction(datatree.symbol_table.getName(symb_id), V_NodeID)));
    else
      return(make_pair(1, (NodeID)NULL ));
  }

NodeID
UnknownFunctionNode::toStatic(DataTree &static_datatree) const
  {
    vector<NodeID> static_arguments;
    for (vector<NodeID>::const_iterator it = arguments.begin();
         it != arguments.end(); it++)
      static_arguments.push_back((*it)->toStatic(static_datatree));
    return static_datatree.AddUnknownFunction(datatree.symbol_table.getName(symb_id), static_arguments);
  }
