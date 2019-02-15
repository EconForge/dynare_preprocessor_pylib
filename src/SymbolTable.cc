/*
 * Copyright (C) 2003-2019 Dynare Team
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

#include <algorithm>
#include <sstream>
#include <iostream>
#include <cassert>
#include <boost/algorithm/string/replace.hpp>
#include <utility>

#include "SymbolTable.hh"

AuxVarInfo::AuxVarInfo(int symb_id_arg, AuxVarType type_arg, int orig_symb_id_arg, int orig_lead_lag_arg,
                       int equation_number_for_multiplier_arg, int information_set_arg,
                       expr_t expr_node_arg, string unary_op_arg) :
  symb_id{symb_id_arg},
  type{type_arg},
  orig_symb_id{orig_symb_id_arg},
  orig_lead_lag{orig_lead_lag_arg},
  equation_number_for_multiplier{equation_number_for_multiplier_arg},
  information_set{information_set_arg},
  expr_node{expr_node_arg},
  unary_op{move(unary_op_arg)}
{
}

SymbolTable::SymbolTable()
= default;

int
SymbolTable::addSymbol(const string &name, SymbolType type, const string &tex_name, const vector<pair<string, string>> &partition_value) noexcept(false)
{
  if (frozen)
    throw FrozenException();

  if (exists(name))
    {
      if (type_table[getID(name)] == type)
        throw AlreadyDeclaredException(name, true);
      else
        throw AlreadyDeclaredException(name, false);
    }

  string final_tex_name = tex_name;
  if (final_tex_name.empty())
    {
      final_tex_name = name;
      size_t pos = 0;
      while ((pos = final_tex_name.find('_', pos)) != string::npos)
        {
          final_tex_name.insert(pos, "\\");
          pos += 2;
        }
    }

  string final_long_name = name;
  bool non_long_name_partition_exists = false;
  for (auto it : partition_value)
    if (it.first == "long_name")
      final_long_name = it.second;
    else
      non_long_name_partition_exists = true;

  int id = symbol_table.size();

  symbol_table[name] = id;
  type_table.push_back(type);
  name_table.push_back(name);
  tex_name_table.push_back(final_tex_name);
  long_name_table.push_back(final_long_name);
  if (non_long_name_partition_exists)
    {
      map<string, string> pmv;
      for (auto it : partition_value)
        pmv[it.first] = it.second;
      partition_value_map[id] = pmv;
    }
  return id;
}

int
SymbolTable::addSymbol(const string &name, SymbolType type) noexcept(false)
{
  return addSymbol(name, type, "", {});
}

void
SymbolTable::freeze() noexcept(false)
{
  if (frozen)
    throw FrozenException();

  frozen = true;

  for (int i = 0; i < (int) symbol_table.size(); i++)
    {
      int tsi;
      switch (getType(i))
        {
        case SymbolType::endogenous:
          tsi = endo_ids.size();
          endo_ids.push_back(i);
          break;
        case SymbolType::exogenous:
          tsi = exo_ids.size();
          exo_ids.push_back(i);
          break;
        case SymbolType::exogenousDet:
          tsi = exo_det_ids.size();
          exo_det_ids.push_back(i);
          break;
        case SymbolType::parameter:
          tsi = param_ids.size();
          param_ids.push_back(i);
          break;
        default:
          tsi = -1;
          break;
        }
      type_specific_ids.push_back(tsi);
    }
}

void
SymbolTable::unfreeze()
{
  frozen = false;
  endo_ids.clear();
  exo_ids.clear();
  exo_det_ids.clear();
  param_ids.clear();
  type_specific_ids.clear();
}

void
SymbolTable::changeType(int id, SymbolType newtype) noexcept(false)
{
  if (frozen)
    throw FrozenException();

  validateSymbID(id);

  type_table[id] = newtype;
}

int
SymbolTable::getID(SymbolType type, int tsid) const noexcept(false)
{
  if (!frozen)
    throw NotYetFrozenException();

  switch (type)
    {
    case SymbolType::endogenous:
      if (tsid < 0 || tsid >= (int) endo_ids.size())
        throw UnknownTypeSpecificIDException(tsid, type);
      else
        return endo_ids[tsid];
    case SymbolType::exogenous:
      if (tsid < 0 || tsid >= (int) exo_ids.size())
        throw UnknownTypeSpecificIDException(tsid, type);
      else
        return exo_ids[tsid];
    case SymbolType::exogenousDet:
      if (tsid < 0 || tsid >= (int) exo_det_ids.size())
        throw UnknownTypeSpecificIDException(tsid, type);
      else
        return exo_det_ids[tsid];
    case SymbolType::parameter:
      if (tsid < 0 || tsid >= (int) param_ids.size())
        throw UnknownTypeSpecificIDException(tsid, type);
      else
        return param_ids[tsid];
    default:
      throw UnknownTypeSpecificIDException(tsid, type);
    }
}

map<string, map<int, string>>
SymbolTable::getPartitionsForType(SymbolType st) const noexcept(false)
{
  map<string, map<int, string>> partitions;
  for (const auto & it : partition_value_map)
    if (getType(it.first) == st)
      for (auto it1 = it.second.begin();
           it1 != it.second.end(); it1++)
        {
          if (partitions.find(it1->first) == partitions.end())
            partitions[it1->first] = map<int, string> ();
          partitions[it1->first][it.first] = it1->second;
        }
  return partitions;
}

void
SymbolTable::writeOutput(ostream &output) const noexcept(false)
{
  if (!frozen)
    throw NotYetFrozenException();

  if (exo_nbr() > 0)
    {
      output << "M_.exo_names = cell(" << exo_nbr() << ",1);" << endl;
      output << "M_.exo_names_tex = cell(" << exo_nbr() << ",1);" << endl;
      output << "M_.exo_names_long = cell(" << exo_nbr() << ",1);" << endl;
      for (int id = 0; id < exo_nbr(); id++)
        output << "M_.exo_names(" << id+1 << ") = {'" << getName(exo_ids[id]) << "'};" << endl
               << "M_.exo_names_tex(" << id+1 << ") = {'" << getTeXName(exo_ids[id]) << "'};" << endl
               << "M_.exo_names_long(" << id+1 << ") = {'" << getLongName(exo_ids[id]) << "'};" << endl;
      map<string, map<int, string>> partitions = getPartitionsForType(SymbolType::exogenous);
      for (map<string, map<int, string>>::const_iterator it = partitions.begin();
           it != partitions.end(); it++)
        if (it->first != "long_name")
          {
            map<int, string>::const_iterator it1;
            output << "M_.exo_partitions." << it->first << " = { ";
            for (int id = 0; id < exo_nbr(); id++)
              {
                output << "'";
                it1 = it->second.find(exo_ids[id]);
                if (it1 != it->second.end())
                  output << it1->second;
                output << "' ";
              }
            output << "};" << endl;
          }
    }

  if (exo_det_nbr() > 0)
    {
      output << "M_.exo_det_names = cell(" << exo_det_nbr() << ",1);" << endl;
      output << "M_.exo_det_names_tex = cell(" << exo_det_nbr() << ",1);" << endl;
      output << "M_.exo_det_names_long = cell(" << exo_det_nbr() << ",1);" << endl;
      for (int id = 0; id < exo_det_nbr(); id++)
        output << "M_.exo_det_names(" << id+1 << ") = {'" << getName(exo_det_ids[id]) << "'};" << endl
               << "M_.exo_det_names_tex(" << id+1 << ") = {'" << getTeXName(exo_det_ids[id]) << "'};" << endl
               << "M_.exo_det_names_long(" << id+1 << ") = {'" << getLongName(exo_det_ids[id]) << "'};" << endl;
      output << "M_.exo_det_partitions = struct();" << endl;
      map<string, map<int, string>> partitions = getPartitionsForType(SymbolType::exogenousDet);
      for (map<string, map<int, string>>::const_iterator it = partitions.begin();
           it != partitions.end(); it++)
        if (it->first != "long_name")
          {
            map<int, string>::const_iterator it1;
            output << "M_.exo_det_partitions." << it->first << " = { ";
            for (int id = 0; id < exo_det_nbr(); id++)
              {
                output << "'";
                it1 = it->second.find(exo_det_ids[id]);
                if (it1 != it->second.end())
                  output << it1->second;
                output << "' ";
              }
            output << "};" << endl;
          }
    }

  if (endo_nbr() > 0)
    {
      output << "M_.endo_names = cell(" << endo_nbr() << ",1);" << endl;
      output << "M_.endo_names_tex = cell(" << endo_nbr() << ",1);" << endl;
      output << "M_.endo_names_long = cell(" << endo_nbr() << ",1);" << endl;
      for (int id = 0; id < endo_nbr(); id++)
        output << "M_.endo_names(" << id+1 << ") = {'" << getName(endo_ids[id]) << "'};" << endl
               << "M_.endo_names_tex(" << id+1 << ") = {'" << getTeXName(endo_ids[id]) << "'};" << endl
               << "M_.endo_names_long(" << id+1 << ") = {'" << getLongName(endo_ids[id]) << "'};" << endl;
      output << "M_.endo_partitions = struct();" << endl;
      map<string, map<int, string>> partitions = getPartitionsForType(SymbolType::endogenous);
      for (map<string, map<int, string>>::const_iterator it = partitions.begin();
           it != partitions.end(); it++)
        if (it->first != "long_name")
          {
            map<int, string>::const_iterator it1;
            output << "M_.endo_partitions." << it->first << " = { ";
            for (int id = 0; id < endo_nbr(); id++)
              {
                output << "'";
                it1 = it->second.find(endo_ids[id]);
                if (it1 != it->second.end())
                  output << it1->second;
                output << "' ";
              }
            output << "};" << endl;
          }
    }

  if (param_nbr() > 0)
    {
      output << "M_.param_names = cell(" << param_nbr() << ",1);" << endl;
      output << "M_.param_names_tex = cell(" << param_nbr() << ",1);" << endl;
      output << "M_.param_names_long = cell(" << param_nbr() << ",1);" << endl;
      for (int id = 0; id < param_nbr(); id++)
        {
          output << "M_.param_names(" << id+1 << ") = {'" << getName(param_ids[id]) << "'};" << endl
                 << "M_.param_names_tex(" << id+1 << ") = {'" << getTeXName(param_ids[id]) << "'};" << endl
                 << "M_.param_names_long(" << id+1 << ") = {'" << getLongName(param_ids[id]) << "'};" << endl;
          if (getName(param_ids[id]) == "dsge_prior_weight")
            output << "options_.dsge_var = 1;" << endl;
        }
      output << "M_.param_partitions = struct();" << endl;
      map<string, map<int, string>> partitions = getPartitionsForType(SymbolType::parameter);
      for (map<string, map<int, string>>::const_iterator it = partitions.begin();
           it != partitions.end(); it++)
        if (it->first != "long_name")
          {
            map<int, string>::const_iterator it1;
            output << "M_.param_partitions." << it->first << " = { ";
            for (int id = 0; id < param_nbr(); id++)
              {
                output << "'";
                it1 = it->second.find(param_ids[id]);
                if (it1 != it->second.end())
                  output << it1->second;
                output << "' ";
              }
            output << "};" << endl;
          }
    }

  output << "M_.exo_det_nbr = " << exo_det_nbr() << ";" << endl
         << "M_.exo_nbr = " << exo_nbr() << ";" << endl
         << "M_.endo_nbr = " << endo_nbr() << ";" << endl
         << "M_.param_nbr = " << param_nbr() << ";" << endl;

  // Write the auxiliary variable table
  output << "M_.orig_endo_nbr = " << orig_endo_nbr() << ";" << endl;
  if (aux_vars.size() == 0)
    output << "M_.aux_vars = [];" << endl;
  else
    for (int i = 0; i < (int) aux_vars.size(); i++)
      {
        output << "M_.aux_vars(" << i+1 << ").endo_index = " << getTypeSpecificID(aux_vars[i].get_symb_id())+1 << ";" << endl
               << "M_.aux_vars(" << i+1 << ").type = " << aux_vars[i].get_type_id() << ";" << endl;
        switch (aux_vars[i].get_type())
          {
          case AuxVarType::endoLead:
          case AuxVarType::exoLead:
            break;
          case AuxVarType::endoLag:
          case AuxVarType::exoLag:
          case AuxVarType::varModel:
            output << "M_.aux_vars(" << i+1 << ").orig_index = " << getTypeSpecificID(aux_vars[i].get_orig_symb_id())+1 << ";" << endl
                   << "M_.aux_vars(" << i+1 << ").orig_lead_lag = " << aux_vars[i].get_orig_lead_lag() << ";" << endl;
            break;
          case AuxVarType::unaryOp:
            if (aux_vars[i].get_orig_symb_id() >= 0)
              output << "M_.aux_vars(" << i+1 << ").orig_index = " << getTypeSpecificID(aux_vars[i].get_orig_symb_id())+1 << ";" << endl
                     << "M_.aux_vars(" << i+1 << ").orig_lead_lag = " << aux_vars[i].get_orig_lead_lag() << ";" << endl;
            output << "M_.aux_vars(" << i+1 << ").unary_op = '" << aux_vars[i].get_unary_op() << "';" << endl;
            break;
          case AuxVarType::multiplier:
            output << "M_.aux_vars(" << i+1 << ").eq_nbr = " << aux_vars[i].get_equation_number_for_multiplier() + 1 << ";" << endl;
            break;
          case AuxVarType::diffForward:
            output << "M_.aux_vars(" << i+1 << ").orig_index = " << getTypeSpecificID(aux_vars[i].get_orig_symb_id())+1 << ";" << endl;
            break;
          case AuxVarType::expectation:
            output << "M_.aux_vars(" << i+1 << ").orig_expr = '\\mathbb{E}_{t"
                   << (aux_vars[i].get_information_set() < 0 ? "" : "+")
                   << aux_vars[i].get_information_set() << "}(";
            aux_vars[i].get_expr_node()->writeOutput(output, ExprNodeOutputType::latexDynamicModel);
            output << ")';" << endl;
            break;
          case AuxVarType::diff:
          case AuxVarType::diffLag:
          case AuxVarType::diffLead:
            if (aux_vars[i].get_orig_symb_id() >= 0)
              output << "M_.aux_vars(" << i+1 << ").orig_index = " << getTypeSpecificID(aux_vars[i].get_orig_symb_id())+1 << ";" << endl
                     << "M_.aux_vars(" << i+1 << ").orig_lead_lag = " << aux_vars[i].get_orig_lead_lag() << ";" << endl;
            break;
          }
      }

  if (predeterminedNbr() > 0)
    {
      output << "M_.predetermined_variables = [ ";
      for (int predetermined_variable : predetermined_variables)
        output << getTypeSpecificID(predetermined_variable)+1 << " ";
      output << "];" << endl;
    }

  if (observedVariablesNbr() > 0)
    {
      int ic = 1;
      output << "options_.varobs = cell(" << observedVariablesNbr() << ", 1);" << endl;
      for (auto it = varobs.begin();
           it != varobs.end(); it++, ic++)
        output << "options_.varobs(" << ic << ")  = {'" << getName(*it) << "'};" << endl;

      output << "options_.varobs_id = [ ";
      for (int varob : varobs)
        output << getTypeSpecificID(varob)+1 << " ";
      output << " ];"  << endl;
    }

  if (observedExogenousVariablesNbr() > 0)
    {
      int ic = 1;
      output << "options_.varexobs = cell(1);" << endl;
      for (auto it = varexobs.begin();
           it != varexobs.end(); it++, ic++)
        output << "options_.varexobs(" << ic << ")  = {'" << getName(*it) << "'};" << endl;

      output << "options_.varexobs_id = [ ";
      for (int varexob : varexobs)
        output << getTypeSpecificID(varexob)+1 << " ";
      output << " ];"  << endl;
    }
}

void
SymbolTable::writeCOutput(ostream &output) const noexcept(false)
{
  if (!frozen)
    throw NotYetFrozenException();

  output << endl
         << "int exo_nbr = " << exo_nbr() << ";" << endl;
  if (exo_nbr() > 0)
    {
      output << "char *exo_names[" << exo_nbr() << "];" << endl;
      for (int id = 0; id < exo_nbr(); id++)
        output << "exo_names[" << id << "] = \"" << getName(exo_ids[id]) << "\";" << endl;
    }

  output << endl
         << "int exo_det_nbr = " << exo_det_nbr() << ";" << endl;
  if (exo_det_nbr() > 0)
    {
      output << "char *exo_det_names[" << exo_det_nbr() << "];" << endl;
      for (int id = 0; id < exo_det_nbr(); id++)
        output << "exo_det_names[" << id << "] = \"" << getName(exo_det_ids[id]) << "\";" << endl;
    }

  output << endl
         << "int endo_nbr = " << endo_nbr() << ";" << endl;
  if (endo_nbr() > 0)
    {
      output << "char *endo_names[" << endo_nbr() << "];" << endl;
      for (int id = 0; id < endo_nbr(); id++)
        output << "endo_names[" << id << "] = \"" << getName(endo_ids[id]) << "\";" << endl;
    }

  output << endl
         << "int param_nbr = " << param_nbr() << ";" << endl;
  if (param_nbr() > 0)
    {
      output << "char *param_names[" << param_nbr() << "];" << endl;
      for (int id = 0; id < param_nbr(); id++)
        output << "param_names[" << id << "] = \"" << getName(param_ids[id]) << "\";" << endl;
    }

  // Write the auxiliary variable table
  output << "int aux_var_nbr = " << aux_vars.size() << ";" << endl;
  if (aux_vars.size() > 0)
    {
      output << "struct aux_vars_t *av[" << aux_vars.size() << "];" << endl;
      for (int i = 0; i < (int) aux_vars.size(); i++)
        {
          output << "av[" << i << "].endo_index = " << getTypeSpecificID(aux_vars[i].get_symb_id()) << ";" << endl
                 << "av[" << i << "].type = " << aux_vars[i].get_type_id() << ";" << endl;
          switch (aux_vars[i].get_type())
            {
            case AuxVarType::endoLead:
            case AuxVarType::exoLead:
            case AuxVarType::expectation:
            case AuxVarType::multiplier:
            case AuxVarType::diffForward:
              break;
            case AuxVarType::endoLag:
            case AuxVarType::exoLag:
            case AuxVarType::varModel:
              output << "av[" << i << "].orig_index = " << getTypeSpecificID(aux_vars[i].get_orig_symb_id()) << ";" << endl
                     << "av[" << i << "].orig_lead_lag = " << aux_vars[i].get_orig_lead_lag() << ";" << endl;
              break;
            case AuxVarType::unaryOp:
              if (aux_vars[i].get_orig_symb_id() >= 0)
                output << "av[" << i << "].orig_index = " << getTypeSpecificID(aux_vars[i].get_orig_symb_id()) << ";" << endl
                       << "av[" << i << "].orig_lead_lag = " << aux_vars[i].get_orig_lead_lag() << ";" << endl;
              output << "av[" << i << "].unary_op = \"" << aux_vars[i].get_unary_op() << "\";" << endl;
              break;
            case AuxVarType::diff:
            case AuxVarType::diffLag:
            case AuxVarType::diffLead:
              if (aux_vars[i].get_orig_symb_id() >= 0)
                output << "av[" << i << "].orig_index = " << getTypeSpecificID(aux_vars[i].get_orig_symb_id()) << ";" << endl
                       << "av[" << i << "].orig_lead_lag = " << aux_vars[i].get_orig_lead_lag() << ";" << endl;
              break;
            }
        }
    }

  output << "int predeterminedNbr = " << predeterminedNbr() << ";" << endl;
  if (predeterminedNbr() > 0)
    {
      output << "int predetermined_variables[" << predeterminedNbr() << "] = {";
      for (auto it = predetermined_variables.begin();
           it != predetermined_variables.end(); it++)
        {
          if (it != predetermined_variables.begin())
            output << ",";
          output << getTypeSpecificID(*it);
        }
      output << "};" << endl;
    }

  output << "int observedVariablesNbr = " << observedVariablesNbr() << ";" << endl;
  if (observedVariablesNbr() > 0)
    {
      output << "int varobs[" << observedVariablesNbr() << "] = {";
      for (auto it = varobs.begin();
           it != varobs.end(); it++)
        {
          if (it != varobs.begin())
            output << ",";
          output << getTypeSpecificID(*it);
        }
      output << "};" << endl;
    }

  output << "int observedExogenousVariablesNbr = " << observedExogenousVariablesNbr() << ";" << endl;
  if (observedExogenousVariablesNbr() > 0)
    {
      output << "int varexobs[" << observedExogenousVariablesNbr() << "] = {";
      for (auto it = varexobs.begin();
           it != varexobs.end(); it++)
        {
          if (it != varexobs.begin())
            output << ",";
          output << getTypeSpecificID(*it);
        }
      output << "};" << endl;
    }
}

void
SymbolTable::writeCCOutput(ostream &output) const noexcept(false)
{
  if (!frozen)
    throw NotYetFrozenException();

  output << endl
         << "exo_nbr = " << exo_nbr() << ";" << endl;
  for (int id = 0; id < exo_nbr(); id++)
    output << "exo_names[\"" << getName(exo_ids[id]) << "\"] = " << id << ";" << endl;

  output << endl
         << "exo_det_nbr = " << exo_det_nbr() << ";" << endl;
  for (int id = 0; id < exo_det_nbr(); id++)
    output << "exo_det_names[\"" << getName(exo_det_ids[id]) << "\"] = " << id << " ;" << endl;

  output << endl
         << "endo_nbr = " << endo_nbr() << ";" << endl;
  for (int id = 0; id < endo_nbr(); id++)
    output << "endo_names[\"" << getName(endo_ids[id]) << "\"] = " << id << ";" << endl;

  output << endl
         << "param_nbr = " << param_nbr() << ";" << endl;
  for (int id = 0; id < param_nbr(); id++)
    output << "param_names[\"" << getName(param_ids[id]) << "\"] = " << id << ";" << endl;

  // Write the auxiliary variable table
  for (int i = 0; i < (int) aux_vars.size(); i++)
    {
      output << "aux_vars_t av" << i << ";" << endl;
      output << "av" << i << ".endo_index = " << getTypeSpecificID(aux_vars[i].get_symb_id()) << ";" << endl
             << "av" << i << ".type = " << aux_vars[i].get_type_id() << ";" << endl;
      switch (aux_vars[i].get_type())
        {
        case AuxVarType::endoLead:
        case AuxVarType::exoLead:
        case AuxVarType::expectation:
        case AuxVarType::multiplier:
        case AuxVarType::diffForward:
          break;
        case AuxVarType::endoLag:
        case AuxVarType::exoLag:
        case AuxVarType::varModel:
          output << "av" << i << ".orig_index = " << getTypeSpecificID(aux_vars[i].get_orig_symb_id()) << ";" << endl
                 << "av" << i << ".orig_lead_lag = " << aux_vars[i].get_orig_lead_lag() << ";" << endl;
          break;
        case AuxVarType::unaryOp:
          if (aux_vars[i].get_orig_symb_id() >= 0)
            output << "av" << i << ".orig_index = " << getTypeSpecificID(aux_vars[i].get_orig_symb_id()) << ";" << endl
                   << "av" << i << ".orig_lead_lag = " << aux_vars[i].get_orig_lead_lag() << ";" << endl;
          output << "av" << i << ".unary_op = \"" << aux_vars[i].get_unary_op() << "\";" << endl;
          break;
        case AuxVarType::diff:
        case AuxVarType::diffLag:
        case AuxVarType::diffLead:
          if (aux_vars[i].get_orig_symb_id() >= 0)
            output << "av" << i << ".orig_index = " << getTypeSpecificID(aux_vars[i].get_orig_symb_id()) << ";" << endl
                   << "av" << i << ".orig_lead_lag = " << aux_vars[i].get_orig_lead_lag() << ";" << endl;
          break;
        }
      output << "aux_vars.push_back(" << "av" << i << ");" << endl;
    }

  for (int predetermined_variable : predetermined_variables)
    output << "predetermined_variables.push_back(" << getTypeSpecificID(predetermined_variable) << ");" << endl;

  for (int varob : varobs)
    output << "varobs.push_back(" << getTypeSpecificID(varob) << ");" << endl;

  for (int varexob : varexobs)
    output << "varexobs.push_back(" << getTypeSpecificID(varexob) << ");" << endl;
}

int
SymbolTable::addLeadAuxiliaryVarInternal(bool endo, int index, expr_t expr_arg) noexcept(false)
{
  ostringstream varname;
  if (endo)
    varname << "AUX_ENDO_LEAD_";
  else
    varname << "AUX_EXO_LEAD_";
  varname << index;
  int symb_id;
  try
    {
      symb_id = addSymbol(varname.str(), SymbolType::endogenous);
    }
  catch (AlreadyDeclaredException &e)
    {
      cerr << "ERROR: you should rename your variable called " << varname.str() << ", this name is internally used by Dynare" << endl;
      exit(EXIT_FAILURE);
    }

  aux_vars.emplace_back(symb_id, (endo ? AuxVarType::endoLead : AuxVarType::exoLead), 0, 0, 0, 0, expr_arg, "");

  return symb_id;
}

int
SymbolTable::addLagAuxiliaryVarInternal(bool endo, int orig_symb_id, int orig_lead_lag, expr_t expr_arg) noexcept(false)
{
  ostringstream varname;
  if (endo)
    varname << "AUX_ENDO_LAG_";
  else
    varname << "AUX_EXO_LAG_";
  varname << orig_symb_id << "_" << -orig_lead_lag;

  int symb_id;
  try
    {
      symb_id = addSymbol(varname.str(), SymbolType::endogenous);
    }
  catch (AlreadyDeclaredException &e)
    {
      cerr << "ERROR: you should rename your variable called " << varname.str() << ", this name is internally used by Dynare" << endl;
      exit(EXIT_FAILURE);
    }

  aux_vars.emplace_back(symb_id, (endo ? AuxVarType::endoLag : AuxVarType::exoLag), orig_symb_id, orig_lead_lag, 0, 0, expr_arg, "");

  return symb_id;
}

int
SymbolTable::addEndoLeadAuxiliaryVar(int index, expr_t expr_arg) noexcept(false)
{
  return addLeadAuxiliaryVarInternal(true, index, expr_arg);
}

int
SymbolTable::addEndoLagAuxiliaryVar(int orig_symb_id, int orig_lead_lag, expr_t expr_arg) noexcept(false)
{
  return addLagAuxiliaryVarInternal(true, orig_symb_id, orig_lead_lag, expr_arg);
}

int
SymbolTable::addExoLeadAuxiliaryVar(int index, expr_t expr_arg) noexcept(false)
{
  return addLeadAuxiliaryVarInternal(false, index, expr_arg);
}

int
SymbolTable::addExoLagAuxiliaryVar(int orig_symb_id, int orig_lead_lag, expr_t expr_arg) noexcept(false)
{
  return addLagAuxiliaryVarInternal(false, orig_symb_id, orig_lead_lag, expr_arg);
}

int
SymbolTable::addExpectationAuxiliaryVar(int information_set, int index, expr_t expr_arg) noexcept(false)
{
  ostringstream varname;
  int symb_id;

  varname << "AUX_EXPECT_" << (information_set < 0 ? "LAG" : "LEAD") << "_"
          << abs(information_set) << "_" << index;

  try
    {
      symb_id = addSymbol(varname.str(), SymbolType::endogenous);
    }
  catch (AlreadyDeclaredException &e)
    {
      cerr << "ERROR: you should rename your variable called " << varname.str() << ", this name is internally used by Dynare" << endl;
      exit(EXIT_FAILURE);
    }

  aux_vars.emplace_back(symb_id, AuxVarType::expectation, 0, 0, 0, information_set, expr_arg, "");

  return symb_id;
}

int
SymbolTable::addDiffLagAuxiliaryVar(int index, expr_t expr_arg, int orig_symb_id, int orig_lag) noexcept(false)
{
  ostringstream varname;
  int symb_id;

  varname << "AUX_DIFF_LAG_" << index;

  try
    {
      symb_id = addSymbol(varname.str(), SymbolType::endogenous);
    }
  catch (AlreadyDeclaredException &e)
    {
      cerr << "ERROR: you should rename your variable called " << varname.str() << ", this name is internally used by Dynare" << endl;
      exit(EXIT_FAILURE);
    }

  aux_vars.emplace_back(symb_id, AuxVarType::diffLag, orig_symb_id, orig_lag, 0, 0, expr_arg, "");

  return symb_id;
}

int
SymbolTable::addDiffLeadAuxiliaryVar(int index, expr_t expr_arg, int orig_symb_id, int orig_lead) noexcept(false)
{
  ostringstream varname;
  int symb_id;

  varname << "AUX_DIFF_LEAD_" << index;

  try
    {
      symb_id = addSymbol(varname.str(), SymbolType::endogenous);
    }
  catch (AlreadyDeclaredException &e)
    {
      cerr << "ERROR: you should rename your variable called " << varname.str() << ", this name is internally used by Dynare" << endl;
      exit(EXIT_FAILURE);
    }

  aux_vars.emplace_back(symb_id, AuxVarType::diffLead, orig_symb_id, orig_lead, 0, 0, expr_arg, "");

  return symb_id;
}

int
SymbolTable::addDiffAuxiliaryVar(int index, expr_t expr_arg, int orig_symb_id, int orig_lag) noexcept(false)
{
  ostringstream varname;
  int symb_id;

  varname << "AUX_DIFF_" << index;

  try
    {
      symb_id = addSymbol(varname.str(), SymbolType::endogenous);
    }
  catch (AlreadyDeclaredException &e)
    {
      cerr << "ERROR: you should rename your variable called " << varname.str() << ", this name is internally used by Dynare" << endl;
      exit(EXIT_FAILURE);
    }

  aux_vars.emplace_back(symb_id, AuxVarType::diff, orig_symb_id, orig_lag, 0, 0, expr_arg, "");

  return symb_id;
}

int
SymbolTable::addDiffAuxiliaryVar(int index, expr_t expr_arg) noexcept(false)
{
  return addDiffAuxiliaryVar(index, expr_arg, -1, 0);
}

int
SymbolTable::addUnaryOpAuxiliaryVar(int index, expr_t expr_arg, string unary_op, int orig_symb_id, int orig_lag) noexcept(false)
{
  ostringstream varname;
  int symb_id;

  varname << "AUX_UOP_" << index;
  try
    {
      symb_id = addSymbol(varname.str(), SymbolType::endogenous);
    }
  catch (AlreadyDeclaredException &e)
    {
      cerr << "ERROR: you should rename your variable called " << varname.str() << ", this name is internally used by Dynare" << endl;
      exit(EXIT_FAILURE);
    }

  aux_vars.emplace_back(symb_id, AuxVarType::unaryOp, orig_symb_id, orig_lag, 0, 0, expr_arg, unary_op);

  return symb_id;
}

int
SymbolTable::addVarModelEndoLagAuxiliaryVar(int orig_symb_id, int orig_lead_lag, expr_t expr_arg) noexcept(false)
{
  int symb_id;
  ostringstream varname;
  varname << "AUX_VARMODEL_" << orig_symb_id << "_" << -orig_lead_lag;

  try
    {
      symb_id = addSymbol(varname.str(), SymbolType::endogenous);
    }
  catch (AlreadyDeclaredException &e)
    {
      cerr << "ERROR: you should rename your variable called " << varname.str() << ", this name is internally used by Dynare" << endl;
      exit(EXIT_FAILURE);
    }

  aux_vars.emplace_back(symb_id, AuxVarType::varModel, orig_symb_id, orig_lead_lag, 0, 0, expr_arg, "");

  return symb_id;
}

int
SymbolTable::addMultiplierAuxiliaryVar(int index) noexcept(false)
{
  ostringstream varname;
  int symb_id;
  varname << "MULT_" << index+1;

  try
    {
      symb_id = addSymbol(varname.str(), SymbolType::endogenous);
    }
  catch (AlreadyDeclaredException &e)
    {
      cerr << "ERROR: you should rename your variable called " << varname.str() << ", this name is internally used by Dynare" << endl;
      exit(EXIT_FAILURE);
    }

  aux_vars.emplace_back(symb_id, AuxVarType::multiplier, 0, 0, index, 0, nullptr, "");
  return symb_id;
}

int
SymbolTable::addDiffForwardAuxiliaryVar(int orig_symb_id, expr_t expr_arg) noexcept(false)
{
  ostringstream varname;
  int symb_id;
  varname << "AUX_DIFF_FWRD_" << orig_symb_id+1;

  try
    {
      symb_id = addSymbol(varname.str(), SymbolType::endogenous);
    }
  catch (AlreadyDeclaredException &e)
    {
      cerr << "ERROR: you should rename your variable called " << varname.str() << ", this name is internally used by Dynare" << endl;
      exit(EXIT_FAILURE);
    }

  aux_vars.emplace_back(symb_id, AuxVarType::diffForward, orig_symb_id, 0, 0, 0, expr_arg, "");
  return symb_id;
}

int
SymbolTable::searchAuxiliaryVars(int orig_symb_id, int orig_lead_lag) const noexcept(false)
{
  for (const auto & aux_var : aux_vars)
    if ((aux_var.get_type() == AuxVarType::endoLag || aux_var.get_type() == AuxVarType::exoLag)
        && aux_var.get_orig_symb_id() == orig_symb_id && aux_var.get_orig_lead_lag() == orig_lead_lag)
      return aux_var.get_symb_id();
  throw SearchFailedException(orig_symb_id, orig_lead_lag);
}

int
SymbolTable::getOrigSymbIdForAuxVar(int aux_var_symb_id) const noexcept(false)
{
  for (const auto & aux_var : aux_vars)
    if ((aux_var.get_type() == AuxVarType::endoLag
         || aux_var.get_type() == AuxVarType::exoLag
         || aux_var.get_type() == AuxVarType::diff
         || aux_var.get_type() == AuxVarType::diffLag
         || aux_var.get_type() == AuxVarType::diffLead)
        && aux_var.get_symb_id() == aux_var_symb_id)
      return aux_var.get_orig_symb_id();
  throw UnknownSymbolIDException(aux_var_symb_id);
}

int
SymbolTable::getOrigLeadLagForDiffAuxVar(int diff_aux_var_symb_id) const noexcept(false)
{
  int lag = 0;
  for (const auto & aux_var : aux_vars)
    if ((aux_var.get_type() == AuxVarType::diffLag || aux_var.get_type() == AuxVarType::diffLead)
        && aux_var.get_symb_id() == diff_aux_var_symb_id)
      lag += 1 + getOrigLeadLagForDiffAuxVar(aux_var.get_orig_symb_id());
  return lag;
}

int
SymbolTable::getOrigSymbIdForDiffAuxVar(int diff_aux_var_symb_id) const noexcept(false)
{
  int orig_symb_id = -1;
  for (const auto & aux_var : aux_vars)
    if (aux_var.get_symb_id() == diff_aux_var_symb_id)
      if (aux_var.get_type() == AuxVarType::diff)
        orig_symb_id = diff_aux_var_symb_id;
      else if (aux_var.get_type() == AuxVarType::diffLag || aux_var.get_type() == AuxVarType::diffLead)
        orig_symb_id = getOrigSymbIdForDiffAuxVar(aux_var.get_orig_symb_id());
  return orig_symb_id;
}

expr_t
SymbolTable::getAuxiliaryVarsExprNode(int symb_id) const noexcept(false)
// throw exception if it is a Lagrange multiplier
{
  for (const auto & aux_var : aux_vars)
    if (aux_var.get_symb_id() == symb_id)
      {
        expr_t expr_node = aux_var.get_expr_node();
        if (expr_node != nullptr)
          return expr_node;
        else
          throw SearchFailedException(symb_id);
      }
  throw SearchFailedException(symb_id);
}

void
SymbolTable::markPredetermined(int symb_id) noexcept(false)
{
  validateSymbID(symb_id);

  if (frozen)
    throw FrozenException();

  assert(getType(symb_id) == SymbolType::endogenous);

  predetermined_variables.insert(symb_id);
}

bool
SymbolTable::isPredetermined(int symb_id) const noexcept(false)
{
  validateSymbID(symb_id);
  return (predetermined_variables.find(symb_id) != predetermined_variables.end());
}

int
SymbolTable::predeterminedNbr() const
{
  return (predetermined_variables.size());
}

void
SymbolTable::addObservedVariable(int symb_id) noexcept(false)
{
  validateSymbID(symb_id);
  assert(getType(symb_id) == SymbolType::endogenous);
  varobs.push_back(symb_id);
}

int
SymbolTable::observedVariablesNbr() const
{
  return (int) varobs.size();
}

bool
SymbolTable::isObservedVariable(int symb_id) const
{
  return (find(varobs.begin(), varobs.end(), symb_id) != varobs.end());
}

int
SymbolTable::getObservedVariableIndex(int symb_id) const
{
  auto it = find(varobs.begin(), varobs.end(), symb_id);
  assert(it != varobs.end());
  return (int) (it - varobs.begin());
}

void
SymbolTable::addObservedExogenousVariable(int symb_id) noexcept(false)
{
  validateSymbID(symb_id);
  assert(getType(symb_id) != SymbolType::endogenous);
  varexobs.push_back(symb_id);
}

int
SymbolTable::observedExogenousVariablesNbr() const
{
  return (int) varexobs.size();
}

bool
SymbolTable::isObservedExogenousVariable(int symb_id) const
{
  return (find(varexobs.begin(), varexobs.end(), symb_id) != varexobs.end());
}

int
SymbolTable::getObservedExogenousVariableIndex(int symb_id) const
{
  auto it = find(varexobs.begin(), varexobs.end(), symb_id);
  assert(it != varexobs.end());
  return (int) (it - varexobs.begin());
}

vector <int>
SymbolTable::getTrendVarIds() const
{
  vector <int> trendVars;
  for (const auto & it : symbol_table)
    if (getType(it.second) == SymbolType::trend || getType(it.second) == SymbolType::logTrend)
      trendVars.push_back(it.second);
  return trendVars;
}

set<int>
SymbolTable::getExogenous() const
{
  set <int> exogs;
  for (const auto & it : symbol_table)
    if (getType(it.second) == SymbolType::exogenous)
      exogs.insert(it.second);
  return exogs;
}

set<int>
SymbolTable::getObservedExogenous() const
{
  set <int> oexogs;
  for (const auto & it : symbol_table)
    if (getType(it.second) == SymbolType::exogenous)
      if (isObservedExogenousVariable(it.second))
        oexogs.insert(it.second);
  return oexogs;
}

set<int>
SymbolTable::getEndogenous() const
{
  set <int> endogs;
  for (const auto & it : symbol_table)
    if (getType(it.second) == SymbolType::endogenous)
      endogs.insert(it.second);
  return endogs;
}

bool
SymbolTable::isAuxiliaryVariable(int symb_id) const
{
  for (const auto & aux_var : aux_vars)
    if (aux_var.get_symb_id() == symb_id)
      return true;
  return false;
}

bool
SymbolTable::isAuxiliaryVariableButNotMultiplier(int symb_id) const
{
  for (const auto & aux_var : aux_vars)
    if (aux_var.get_symb_id() == symb_id && aux_var.get_type() != AuxVarType::multiplier)
      return true;
  return false;
}

set<int>
SymbolTable::getOrigEndogenous() const
{
  set <int> origendogs;
  for (const auto & it : symbol_table)
    if (getType(it.second) == SymbolType::endogenous && !isAuxiliaryVariable(it.second))
      origendogs.insert(it.second);
  return origendogs;
}

void
SymbolTable::writeJuliaOutput(ostream &output) const noexcept(false)
{
  if (!frozen)
    throw NotYetFrozenException();

  output << "# Endogenous Variables" << endl
         << "model_.endo = [" << endl;
  if (endo_nbr() > 0)
    for (int id = 0; id < endo_nbr(); id++)
      output << "              DynareModel.Endo(\""
             << getName(endo_ids[id]) << "\", raw\""
             << getTeXName(endo_ids[id]) << "\", \""
             << getLongName(endo_ids[id]) << "\")" << endl;
  output << "             ]" << endl;
  output << "model_.endo_nbr = " << endo_nbr() << ";" << endl;

  output << "# Exogenous Variables" << endl
         << "model_.exo = [" << endl;
  if (exo_nbr() > 0)
    for (int id = 0; id < exo_nbr(); id++)
      output << "             DynareModel.Exo(\""
             << getName(exo_ids[id]) << "\", raw\""
             << getTeXName(exo_ids[id]) << "\", \""
             << getLongName(exo_ids[id]) << "\")" << endl;
  output << "            ]" << endl;
  output << "model_.exo_nbr = " << exo_nbr() << ";" << endl;

  if (exo_det_nbr() > 0)
    {
      output << "# Exogenous Deterministic Variables" << endl
             << "model_.exo_det = [" << endl;
      if (exo_det_nbr() > 0)
        for (int id = 0; id < exo_det_nbr(); id++)
          output << "                 DynareModel.ExoDet(\""
                 << getName(exo_det_ids[id]) << "\", raw\""
                 << getTeXName(exo_det_ids[id]) << "\", \""
                 << getLongName(exo_det_ids[id]) << "\")" << endl;
      output << "                ]" << endl;
      output << "model_.exo_det_nbr = " << exo_det_nbr() << ";" << endl;
    }

  output << "# Parameters" << endl
         << "model_.param = [" << endl;
  if (param_nbr() > 0)
    for (int id = 0; id < param_nbr(); id++)
      output << "               DynareModel.Param(\""
             << getName(param_ids[id]) << "\", raw\""
             << getTeXName(param_ids[id]) << "\", \""
             << getLongName(param_ids[id]) << "\")" << endl;
  output << "              ]" << endl;
  output << "model_.param_nbr = " << param_nbr() << ";" << endl;

  output << "model_.orig_endo_nbr = " << orig_endo_nbr() << endl;

  if (aux_vars.size() > 0)
    {
      output << "# Auxiliary Variables" << endl
             << "model_.aux_vars = [" << endl;
      for (const auto & aux_var : aux_vars)
        {
          output << "                   DynareModel.AuxVars("
                 << getTypeSpecificID(aux_var.get_symb_id()) + 1 << ", "
                 << aux_var.get_type_id() << ", ";
          switch (aux_var.get_type())
            {
            case AuxVarType::endoLead:
            case AuxVarType::exoLead:
            case AuxVarType::endoLag:
            case AuxVarType::exoLag:
            case AuxVarType::varModel:
              output << getTypeSpecificID(aux_var.get_orig_symb_id()) + 1 << ", "
                     << aux_var.get_orig_lead_lag() << ", typemin(Int), string(), string()";
              break;
            case AuxVarType::unaryOp:
              if (aux_var.get_orig_symb_id() >= 0)
                output << getTypeSpecificID(aux_var.get_orig_symb_id()) + 1 << ", " << aux_var.get_orig_lead_lag();
              else
                output << "typemin(Int), typemin(Int)";
              output << ", typemin(Int), string(), "
                     << "\"" << aux_var.get_unary_op() << "\"" << endl;
              break;
            case AuxVarType::diff:
            case AuxVarType::diffLag:
            case AuxVarType::diffLead:
              if (aux_var.get_orig_symb_id() >= 0)
                output << getTypeSpecificID(aux_var.get_orig_symb_id()) + 1 << ", "
                       << aux_var.get_orig_lead_lag() << ", typemin(Int), string(), string()";
              break;
            case AuxVarType::multiplier:
              output << "typemin(Int), typemin(Int), " << aux_var.get_equation_number_for_multiplier() + 1
                     << ", string(), string()";
              break;
            case AuxVarType::diffForward:
              output << getTypeSpecificID(aux_var.get_orig_symb_id())+1 << ", typemin(Int), typemin(Int), string(), string()";
              break;
            case AuxVarType::expectation:
              output << "typemin(Int), typemin(Int), typemin(Int), \"\\mathbb{E}_{t"
                     << (aux_var.get_information_set() < 0 ? "" : "+")
                     << aux_var.get_information_set() << "}(";
              aux_var.get_expr_node()->writeOutput(output, ExprNodeOutputType::latexDynamicModel);
              output << ")\"";
              break;
            default:
              output << " typemin(Int), typemin(Int), typemin(Int), string(), string()";
            }
          output << ")" << endl;
        }
      output << "]" << endl;
    }

  if (predeterminedNbr() > 0)
    {
      output << "# Predetermined Variables" << endl
             << "model_.pred_vars = [ " << endl;
      for (int predetermined_variable : predetermined_variables)
        output << "                   DynareModel.PredVars("
               << getTypeSpecificID(predetermined_variable)+1 << ")" << endl;
      output << "                  ]" << endl;
    }

  if (observedVariablesNbr() > 0)
    {
      output << "# Observed Variables" << endl
             << "options_.obs_vars = [" << endl;
      for (int varob : varobs)
        output << "                    DynareModel.ObsVars("
               << getTypeSpecificID(varob)+1 << ")" << endl;
      output << "                   ]" << endl;
    }
}

void
SymbolTable::writeJsonOutput(ostream &output) const
{
  output << "\"endogenous\": ";
  writeJsonVarVector(output, endo_ids);

  output << ", \"exogenous\":";
  writeJsonVarVector(output, exo_ids);

  output << ", \"exogenous_deterministic\": ";
  writeJsonVarVector(output, exo_det_ids);

  output << ", \"parameters\": ";
  writeJsonVarVector(output, param_ids);
}

void
SymbolTable::writeJsonVarVector(ostream &output, const vector<int> &varvec) const
{
  output << "[";
  for (size_t i = 0; i < varvec.size(); i++)
    {
      if (i != 0)
        output << ", ";
      output << "{"
             << "\"name\":\"" << getName(varvec[i]) << "\", "
             << "\"texName\":\"" << boost::replace_all_copy(getTeXName(varvec[i]), "\\", "\\\\") << "\", "
             << "\"longName\":\"" << boost::replace_all_copy(getLongName(varvec[i]), "\\", "\\\\") << "\"}"
             << endl;
    }
  output << "]" << endl;
}
