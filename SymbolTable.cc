/*
 * Copyright (C) 2003-2017 Dynare Team
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

#include "SymbolTable.hh"

AuxVarInfo::AuxVarInfo(int symb_id_arg, aux_var_t type_arg, int orig_symb_id_arg, int orig_lead_lag_arg,
                       int equation_number_for_multiplier_arg, int information_set_arg,
                       expr_t expr_node_arg) :
  symb_id(symb_id_arg),
  type(type_arg),
  orig_symb_id(orig_symb_id_arg),
  orig_lead_lag(orig_lead_lag_arg),
  equation_number_for_multiplier(equation_number_for_multiplier_arg),
  information_set(information_set_arg),
  expr_node(expr_node_arg)
{
}

SymbolTable::SymbolTable() : frozen(false)
{
}

int
SymbolTable::addSymbol(const string &name, SymbolType type, const string &tex_name, const vector<pair<string *, string *> *> *partition_value) throw (AlreadyDeclaredException, FrozenException)
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
  if (partition_value)
    for (vector<pair<string *, string *> *>::const_iterator it = partition_value->begin();
         it != partition_value->end(); it++)
      if (*((*it)->first) == "long_name")
        final_long_name = *((*it)->second);
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
      for (vector<pair<string *, string *> *>::const_iterator it = partition_value->begin();
           it != partition_value->end(); it++)
        pmv[*((*it)->first)] = *((*it)->second);
      partition_value_map[id] = pmv;
    }
  return id;
}

int
SymbolTable::addSymbol(const string &name, SymbolType type) throw (AlreadyDeclaredException, FrozenException)
{
  return addSymbol(name, type, "", NULL);
}

void
SymbolTable::freeze() throw (FrozenException)
{
  if (frozen)
    throw FrozenException();

  frozen = true;

  for (int i = 0; i < symbol_table.size(); i++)
    {
      int tsi;
      switch (getType(i))
        {
        case eEndogenous:
          tsi = endo_ids.size();
          endo_ids.push_back(i);
          break;
        case eExogenous:
          tsi = exo_ids.size();
          exo_ids.push_back(i);
          break;
        case eExogenousDet:
          tsi = exo_det_ids.size();
          exo_det_ids.push_back(i);
          break;
        case eParameter:
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
SymbolTable::changeType(int id, SymbolType newtype) throw (UnknownSymbolIDException, FrozenException)
{
  if (frozen)
    throw FrozenException();

  validateSymbID(id);

  type_table[id] = newtype;
}

int
SymbolTable::getID(SymbolType type, int tsid) const throw (UnknownTypeSpecificIDException, NotYetFrozenException)
{
  if (!frozen)
    throw NotYetFrozenException();

  switch (type)
    {
    case eEndogenous:
      if (tsid < 0 || tsid >= (int) endo_ids.size())
        throw UnknownTypeSpecificIDException(tsid, type);
      else
        return endo_ids[tsid];
    case eExogenous:
      if (tsid < 0 || tsid >= (int) exo_ids.size())
        throw UnknownTypeSpecificIDException(tsid, type);
      else
        return exo_ids[tsid];
    case eExogenousDet:
      if (tsid < 0 || tsid >= (int) exo_det_ids.size())
        throw UnknownTypeSpecificIDException(tsid, type);
      else
        return exo_det_ids[tsid];
    case eParameter:
      if (tsid < 0 || tsid >= (int) param_ids.size())
        throw UnknownTypeSpecificIDException(tsid, type);
      else
        return param_ids[tsid];
    default:
      throw UnknownTypeSpecificIDException(tsid, type);
    }
}

map<string, map<int, string> >
SymbolTable::getPartitionsForType(enum SymbolType st) const throw (UnknownSymbolIDException)
{
  map<string, map<int, string> > partitions;
  for (map<int, map<string, string> >::const_iterator it = partition_value_map.begin();
       it != partition_value_map.end(); it++)
    if (getType(it->first) == st)
      for (map<string, string>::const_iterator it1 = it->second.begin();
           it1 != it->second.end(); it1++)
        {
          if (partitions.find(it1->first) == partitions.end())
            partitions[it1->first] = map<int, string> ();
          partitions[it1->first][it->first] = it1->second;
        }
  return partitions;
}

void
SymbolTable::writeOutput(ostream &output) const throw (NotYetFrozenException)
{
  if (!frozen)
    throw NotYetFrozenException();

  if (exo_nbr() > 0)
    {
      output << "M_.exo_names = '" << getName(exo_ids[0]) << "';" << endl;
      output << "M_.exo_names_tex = '" << getTeXName(exo_ids[0]) << "';" << endl;
      output << "M_.exo_names_long = '" << getLongName(exo_ids[0]) << "';" << endl;
      for (int id = 1; id < exo_nbr(); id++)
        output << "M_.exo_names = char(M_.exo_names, '" << getName(exo_ids[id]) << "');" << endl
               << "M_.exo_names_tex = char(M_.exo_names_tex, '" << getTeXName(exo_ids[id]) << "');" << endl
               << "M_.exo_names_long = char(M_.exo_names_long, '" << getLongName(exo_ids[id]) << "');" << endl;

      map<string, map<int, string> > partitions = getPartitionsForType(eExogenous);
      for (map<string, map<int, string> >::const_iterator it = partitions.begin();
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
      output << "M_.exo_det_names = '" << getName(exo_det_ids[0]) << "';" << endl;
      output << "M_.exo_det_names_tex = '" << getTeXName(exo_det_ids[0]) << "';" << endl;
      output << "M_.exo_det_names_long = '" << getLongName(exo_det_ids[0]) << "';" << endl;
      for (int id = 1; id < exo_det_nbr(); id++)
        output << "M_.exo_det_names = char(M_.exo_det_names, '" << getName(exo_det_ids[id]) << "');" << endl
               << "M_.exo_det_names_tex = char(M_.exo_det_names_tex, '" << getTeXName(exo_det_ids[id]) << "');" << endl
               << "M_.exo_det_names_long = char(M_.exo_det_names_long, '" << getLongName(exo_det_ids[id]) << "');" << endl;

      output << "M_.exo_det_partitions = struct();" << endl;
      map<string, map<int, string> > partitions = getPartitionsForType(eExogenousDet);
      for (map<string, map<int, string> >::const_iterator it = partitions.begin();
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
      output << "M_.endo_names = '" << getName(endo_ids[0]) << "';" << endl;
      output << "M_.endo_names_tex = '" << getTeXName(endo_ids[0]) << "';" << endl;
      output << "M_.endo_names_long = '" << getLongName(endo_ids[0]) << "';" << endl;
      for (int id = 1; id < endo_nbr(); id++)
        output << "M_.endo_names = char(M_.endo_names, '" << getName(endo_ids[id]) << "');" << endl
               << "M_.endo_names_tex = char(M_.endo_names_tex, '" << getTeXName(endo_ids[id]) << "');" << endl
               << "M_.endo_names_long = char(M_.endo_names_long, '" << getLongName(endo_ids[id]) << "');" << endl;

      output << "M_.endo_partitions = struct();" << endl;
      map<string, map<int, string> > partitions = getPartitionsForType(eEndogenous);
      for (map<string, map<int, string> >::const_iterator it = partitions.begin();
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
      output << "M_.param_names = '" << getName(param_ids[0]) << "';" << endl;
      output << "M_.param_names_tex = '" << getTeXName(param_ids[0]) << "';" << endl;
      output << "M_.param_names_long = '" << getLongName(param_ids[0]) << "';" << endl;
      for (int id = 1; id < param_nbr(); id++)
        {
          output << "M_.param_names = char(M_.param_names, '" << getName(param_ids[id]) << "');" << endl
                 << "M_.param_names_tex = char(M_.param_names_tex, '" << getTeXName(param_ids[id]) << "');" << endl
                 << "M_.param_names_long = char(M_.param_names_long, '" << getLongName(param_ids[id]) << "');" << endl;

          if (getName(param_ids[id]) == "dsge_prior_weight")
            output << "options_.dsge_var = 1;" << endl;
        }

      output << "M_.param_partitions = struct();" << endl;
      map<string, map<int, string> > partitions = getPartitionsForType(eParameter);
      for (map<string, map<int, string> >::const_iterator it = partitions.begin();
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
               << "M_.aux_vars(" << i+1 << ").type = " << aux_vars[i].get_type() << ";" << endl;
        switch (aux_vars[i].get_type())
          {
          case avEndoLead:
          case avExoLead:
            break;
          case avEndoLag:
          case avExoLag:
          case avVarModel:
            output << "M_.aux_vars(" << i+1 << ").orig_index = " << getTypeSpecificID(aux_vars[i].get_orig_symb_id())+1 << ";" << endl
                   << "M_.aux_vars(" << i+1 << ").orig_lead_lag = " << aux_vars[i].get_orig_lead_lag() << ";" << endl;
            break;
          case avMultiplier:
            output << "M_.aux_vars(" << i+1 << ").eq_nbr = " << aux_vars[i].get_equation_number_for_multiplier() + 1 << ";" << endl;
            break;
          case avDiffForward:
            output << "M_.aux_vars(" << i+1 << ").orig_index = " << getTypeSpecificID(aux_vars[i].get_orig_symb_id())+1 << ";" << endl;
            break;
          case avExpectation:
            output << "M_.aux_vars(" << i+1 << ").orig_expr = '\\mathbb{E}_{t"
                   << (aux_vars[i].get_information_set() < 0 ? "" : "+")
                   << aux_vars[i].get_information_set() << "}(";
            aux_vars[i].get_expr_node()->writeOutput(output, oLatexDynamicModel);
            output << ")';" << endl;
            break;
          }
      }

  if (predeterminedNbr() > 0)
    {
      output << "M_.predetermined_variables = [ ";
      for (set<int>::const_iterator it = predetermined_variables.begin();
           it != predetermined_variables.end(); it++)
        output << getTypeSpecificID(*it)+1 << " ";
      output << "];" << endl;
    }

  if (observedVariablesNbr() > 0)
    {
      int ic = 1;
      output << "options_.varobs = cell(1);" << endl;
      for (vector<int>::const_iterator it = varobs.begin();
           it != varobs.end(); it++, ic++)
        output << "options_.varobs(" << ic << ")  = {'" << getName(*it) << "'};" << endl;

      output << "options_.varobs_id = [ ";
      for (vector<int>::const_iterator it = varobs.begin();
           it != varobs.end(); it++)
        output << getTypeSpecificID(*it)+1 << " ";
      output << " ];"  << endl;
    }

  if (observedExogenousVariablesNbr() > 0)
    {
      int ic = 1;
      output << "options_.varexobs = cell(1);" << endl;
      for (vector<int>::const_iterator it = varexobs.begin();
           it != varexobs.end(); it++, ic++)
        output << "options_.varexobs(" << ic << ")  = {'" << getName(*it) << "'};" << endl;

      output << "options_.varexobs_id = [ ";
      for (vector<int>::const_iterator it = varexobs.begin();
           it != varexobs.end(); it++)
        output << getTypeSpecificID(*it)+1 << " ";
      output << " ];"  << endl;
    }
}

void
SymbolTable::writeCOutput(ostream &output) const throw (NotYetFrozenException)
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
                 << "av[" << i << "].type = " << aux_vars[i].get_type() << ";" << endl;
          switch (aux_vars[i].get_type())
            {
            case avEndoLead:
            case avExoLead:
            case avExpectation:
            case avMultiplier:
            case avDiffForward:
              break;
            case avEndoLag:
            case avExoLag:
            case avVarModel:
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
      for (set<int>::const_iterator it = predetermined_variables.begin();
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
      for (vector<int>::const_iterator it = varobs.begin();
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
      for (vector<int>::const_iterator it = varexobs.begin();
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
SymbolTable::writeCCOutput(ostream &output) const throw (NotYetFrozenException)
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
             << "av" << i << ".type = " << aux_vars[i].get_type() << ";" << endl;
      switch (aux_vars[i].get_type())
        {
        case avEndoLead:
        case avExoLead:
        case avExpectation:
        case avMultiplier:
        case avDiffForward:
          break;
        case avEndoLag:
        case avExoLag:
        case avVarModel:
          output << "av" << i << ".orig_index = " << getTypeSpecificID(aux_vars[i].get_orig_symb_id()) << ";" << endl
                 << "av" << i << ".orig_lead_lag = " << aux_vars[i].get_orig_lead_lag() << ";" << endl;
          break;
        }
      output << "aux_vars.push_back(" << "av" << i << ");" << endl;
    }

  for (set<int>::const_iterator it = predetermined_variables.begin();
       it != predetermined_variables.end(); it++)
    output << "predetermined_variables.push_back(" << getTypeSpecificID(*it) << ");" << endl;

  for (vector<int>::const_iterator it = varobs.begin();
       it != varobs.end(); it++)
    output << "varobs.push_back(" << getTypeSpecificID(*it) << ");" << endl;

  for (vector<int>::const_iterator it = varexobs.begin();
       it != varexobs.end(); it++)
    output << "varexobs.push_back(" << getTypeSpecificID(*it) << ");" << endl;
}

int
SymbolTable::addLeadAuxiliaryVarInternal(bool endo, int index, expr_t expr_arg) throw (FrozenException)
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
      symb_id = addSymbol(varname.str(), eEndogenous);
    }
  catch (AlreadyDeclaredException &e)
    {
      cerr << "ERROR: you should rename your variable called " << varname.str() << ", this name is internally used by Dynare" << endl;
      exit(EXIT_FAILURE);
    }

  aux_vars.push_back(AuxVarInfo(symb_id, (endo ? avEndoLead : avExoLead), 0, 0, 0, 0, expr_arg));

  return symb_id;
}

int
SymbolTable::addLagAuxiliaryVarInternal(bool endo, int orig_symb_id, int orig_lead_lag, expr_t expr_arg) throw (FrozenException)
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
      symb_id = addSymbol(varname.str(), eEndogenous);
    }
  catch (AlreadyDeclaredException &e)
    {
      cerr << "ERROR: you should rename your variable called " << varname.str() << ", this name is internally used by Dynare" << endl;
      exit(EXIT_FAILURE);
    }

  aux_vars.push_back(AuxVarInfo(symb_id, (endo ? avEndoLag : avExoLag), orig_symb_id, orig_lead_lag, 0, 0, expr_arg));

  return symb_id;
}

int
SymbolTable::addEndoLeadAuxiliaryVar(int index, expr_t expr_arg) throw (FrozenException)
{
  return addLeadAuxiliaryVarInternal(true, index, expr_arg);
}

int
SymbolTable::addEndoLagAuxiliaryVar(int orig_symb_id, int orig_lead_lag, expr_t expr_arg) throw (FrozenException)
{
  return addLagAuxiliaryVarInternal(true, orig_symb_id, orig_lead_lag, expr_arg);
}

int
SymbolTable::addExoLeadAuxiliaryVar(int index, expr_t expr_arg) throw (FrozenException)
{
  return addLeadAuxiliaryVarInternal(false, index, expr_arg);
}

int
SymbolTable::addExoLagAuxiliaryVar(int orig_symb_id, int orig_lead_lag, expr_t expr_arg) throw (FrozenException)
{
  return addLagAuxiliaryVarInternal(false, orig_symb_id, orig_lead_lag, expr_arg);
}

int
SymbolTable::addExpectationAuxiliaryVar(int information_set, int index, expr_t expr_arg) throw (FrozenException)
{
  ostringstream varname;
  int symb_id;

  varname << "AUX_EXPECT_" << (information_set < 0 ? "LAG" : "LEAD") << "_"
          << abs(information_set) << "_" << index;

  try
    {
      symb_id = addSymbol(varname.str(), eEndogenous);
    }
  catch (AlreadyDeclaredException &e)
    {
      cerr << "ERROR: you should rename your variable called " << varname.str() << ", this name is internally used by Dynare" << endl;
      exit(EXIT_FAILURE);
    }

  aux_vars.push_back(AuxVarInfo(symb_id, avExpectation, 0, 0, 0, information_set, expr_arg));

  return symb_id;
}

int
SymbolTable::addVarModelEndoLagAuxiliaryVar(int orig_symb_id, int orig_lead_lag, expr_t expr_arg) throw (AlreadyDeclaredException, FrozenException)
{
  int symb_id;
  ostringstream varname;
  varname << "AUX_VARMODEL_" << orig_symb_id << "_" << -orig_lead_lag;

  try
    {
      symb_id = addSymbol(varname.str(), eEndogenous);
    }
  catch (AlreadyDeclaredException &e)
    {
      cerr << "ERROR: you should rename your variable called " << varname.str() << ", this name is internally used by Dynare" << endl;
      exit(EXIT_FAILURE);
    }

  aux_vars.push_back(AuxVarInfo(symb_id, avVarModel, orig_symb_id, orig_lead_lag, 0, 0, expr_arg));

  return symb_id;
}

int
SymbolTable::addMultiplierAuxiliaryVar(int index) throw (FrozenException)
{
  ostringstream varname;
  int symb_id;
  varname << "MULT_" << index+1;

  try
    {
      symb_id = addSymbol(varname.str(), eEndogenous);
    }
  catch (AlreadyDeclaredException &e)
    {
      cerr << "ERROR: you should rename your variable called " << varname.str() << ", this name is internally used by Dynare" << endl;
      exit(EXIT_FAILURE);
    }

  aux_vars.push_back(AuxVarInfo(symb_id, avMultiplier, 0, 0, index, 0, NULL));
  return symb_id;
}

int
SymbolTable::addDiffForwardAuxiliaryVar(int orig_symb_id, expr_t expr_arg) throw (FrozenException)
{
  ostringstream varname;
  int symb_id;
  varname << "AUX_DIFF_FWRD_" << orig_symb_id+1;

  try
    {
      symb_id = addSymbol(varname.str(), eEndogenous);
    }
  catch (AlreadyDeclaredException &e)
    {
      cerr << "ERROR: you should rename your variable called " << varname.str() << ", this name is internally used by Dynare" << endl;
      exit(EXIT_FAILURE);
    }

  aux_vars.push_back(AuxVarInfo(symb_id, avDiffForward, orig_symb_id, 0, 0, 0, expr_arg));
  return symb_id;
}

int
SymbolTable::searchAuxiliaryVars(int orig_symb_id, int orig_lead_lag) const throw (SearchFailedException)
{
  for (size_t i = 0; i < aux_vars.size(); i++)
    if ((aux_vars[i].get_type() == avEndoLag || aux_vars[i].get_type() == avExoLag)
        && aux_vars[i].get_orig_symb_id() == orig_symb_id && aux_vars[i].get_orig_lead_lag() == orig_lead_lag)
      return aux_vars[i].get_symb_id();
  throw SearchFailedException(orig_symb_id, orig_lead_lag);
}

expr_t
SymbolTable::getAuxiliaryVarsExprNode(int symb_id) const throw (SearchFailedException)
// throw exception if it is a Lagrange multiplier
{
  for (size_t i = 0; i < aux_vars.size(); i++)
    if (aux_vars[i].get_symb_id() == symb_id)
      {
        expr_t expr_node = aux_vars[i].get_expr_node();
        if (expr_node != NULL)
          return expr_node;
        else
          throw SearchFailedException(symb_id);
      }
  throw SearchFailedException(symb_id);
}

void
SymbolTable::markPredetermined(int symb_id) throw (UnknownSymbolIDException, FrozenException)
{
  validateSymbID(symb_id);

  if (frozen)
    throw FrozenException();

  assert(getType(symb_id) == eEndogenous);

  predetermined_variables.insert(symb_id);
}

bool
SymbolTable::isPredetermined(int symb_id) const throw (UnknownSymbolIDException)
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
SymbolTable::addObservedVariable(int symb_id) throw (UnknownSymbolIDException)
{
  validateSymbID(symb_id);
  assert(getType(symb_id) == eEndogenous);
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
  vector<int>::const_iterator it = find(varobs.begin(), varobs.end(), symb_id);
  assert(it != varobs.end());
  return (int) (it - varobs.begin());
}

void
SymbolTable::addObservedExogenousVariable(int symb_id) throw (UnknownSymbolIDException)
{
  validateSymbID(id);
  assert(getType(symb_id) != eEndogenous);
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
  vector<int>::const_iterator it = find(varexobs.begin(), varexobs.end(), symb_id);
  assert(it != varexobs.end());
  return (int) (it - varexobs.begin());
}

vector <int>
SymbolTable::getTrendVarIds() const
{
  vector <int> trendVars;
  for (symbol_table_type::const_iterator it = symbol_table.begin();
       it != symbol_table.end(); it++)
    if (getType(it->second) == eTrend || getType(it->second) == eLogTrend)
      trendVars.push_back(it->second);
  return trendVars;
}

set<int>
SymbolTable::getExogenous() const
{
  set <int> exogs;
  for (symbol_table_type::const_iterator it = symbol_table.begin();
       it != symbol_table.end(); it++)
    if (getType(it->second) == eExogenous)
      exogs.insert(it->second);
  return exogs;
}

set<int>
SymbolTable::getObservedExogenous() const
{
  set <int> oexogs;
  for (symbol_table_type::const_iterator it = symbol_table.begin();
       it != symbol_table.end(); it++)
    if (getType(it->second) == eExogenous)
      if (isObservedExogenousVariable(it->second))
        oexogs.insert(it->second);
  return oexogs;
}

set<int>
SymbolTable::getEndogenous() const
{
  set <int> endogs;
  for (symbol_table_type::const_iterator it = symbol_table.begin();
       it != symbol_table.end(); it++)
    if (getType(it->second) == eEndogenous)
      endogs.insert(it->second);
  return endogs;
}

bool
SymbolTable::isAuxiliaryVariable(int symb_id) const
{
  for (int i = 0; i < (int) aux_vars.size(); i++)
    if (aux_vars[i].get_symb_id() == symb_id)
      return true;
  return false;
}

bool
SymbolTable::isAuxiliaryVariableButNotMultiplier(int symb_id) const
{
  for (int i = 0; i < (int) aux_vars.size(); i++)
    if (aux_vars[i].get_symb_id() == symb_id && aux_vars[i].get_type() != avMultiplier)
      return true;
  return false;
}

set<int>
SymbolTable::getOrigEndogenous() const
{
  set <int> origendogs;
  for (symbol_table_type::const_iterator it = symbol_table.begin();
       it != symbol_table.end(); it++)
    if (getType(it->second) == eEndogenous && !isAuxiliaryVariable(it->second))
      origendogs.insert(it->second);
  return origendogs;
}

void
SymbolTable::writeJuliaOutput(ostream &output) const throw (NotYetFrozenException)
{
  if (!frozen)
    throw NotYetFrozenException();

  output << "# Endogenous Variables" << endl
         << "model_.endo = [" << endl;
  if (endo_nbr() > 0)
    for (int id = 0; id < endo_nbr(); id++)
      output << "              DynareModel.Endo(\""
             << getName(endo_ids[id]) << "\", \""
             << getTeXName(endo_ids[id]) << "\", \""
             << getLongName(endo_ids[id]) << "\")" << endl;
  output << "             ]" << endl;

  output << "# Exogenous Variables" << endl
         << "model_.exo = [" << endl;
  if (exo_nbr() > 0)
    for (int id = 0; id < exo_nbr(); id++)
      output << "             DynareModel.Exo(\""
             << getName(exo_ids[id]) << "\", \""
             << getTeXName(exo_ids[id]) << "\", \""
             << getLongName(exo_ids[id]) << "\")" << endl;
  output << "            ]" << endl;

  if (exo_det_nbr() > 0)
    {
      output << "# Exogenous Deterministic Variables" << endl
             << "model_.exo_det = [" << endl;
      if (exo_det_nbr() > 0)
        for (int id = 0; id < exo_det_nbr(); id++)
          output << "                 DynareModel.ExoDet(\""
                 << getName(exo_det_ids[id]) << "\", \""
                 << getTeXName(exo_det_ids[id]) << "\", \""
                 << getLongName(exo_det_ids[id]) << "\")" << endl;
      output << "                ]" << endl;
    }

  output << "# Parameters" << endl
         << "model_.param = [" << endl;
  if (param_nbr() > 0)
    for (int id = 0; id < param_nbr(); id++)
      output << "               DynareModel.Param(\""
             << getName(param_ids[id]) << "\", \""
             << getTeXName(param_ids[id]) << "\", \""
             << getLongName(param_ids[id]) << "\")" << endl;
  output << "              ]" << endl;

  output << "model_.orig_endo_nbr = " << orig_endo_nbr() << endl;

  if (aux_vars.size() > 0)
    {
      output << "# Auxiliary Variables" << endl
             << "model_.aux_vars = [" << endl;
      for (int i = 0; i < (int) aux_vars.size(); i++)
        {
          output << "                   DynareModel.AuxVars("
                 << getTypeSpecificID(aux_vars[i].get_symb_id()) + 1 << ", "
                 << aux_vars[i].get_type() << ", ";
          switch (aux_vars[i].get_type())
            {
            case avEndoLead:
            case avExoLead:
              break;
            case avEndoLag:
            case avExoLag:
            case avVarModel:
              output << getTypeSpecificID(aux_vars[i].get_orig_symb_id()) + 1 << ", "
                     << aux_vars[i].get_orig_lead_lag() << ", NaN, NaN";
              break;
            case avMultiplier:
              output << "NaN, NaN, " << aux_vars[i].get_equation_number_for_multiplier() + 1
                     << ", NaN";
              break;
            case avDiffForward:
              output << getTypeSpecificID(aux_vars[i].get_orig_symb_id())+1 << ", NaN, ";
              break;
            case avExpectation:
              output << "NaN, NaN, NaN, \"\\mathbb{E}_{t"
                     << (aux_vars[i].get_information_set() < 0 ? "" : "+")
                     << aux_vars[i].get_information_set() << "}(";
              aux_vars[i].get_expr_node()->writeOutput(output, oLatexDynamicModel);
              output << ")\"";
              break;
            }
          output << ")" << endl;
        }
      output << "]" << endl;
    }

  if (predeterminedNbr() > 0)
    {
      output << "# Predetermined Variables" << endl
             << "model_.pred_vars = [ " << endl;
      for (set<int>::const_iterator it = predetermined_variables.begin();
           it != predetermined_variables.end(); it++)
        output << "                   DynareModel.PredVars("
               << getTypeSpecificID(*it)+1 << ")" << endl;
      output << "                  ]" << endl;
    }

  if (observedVariablesNbr() > 0)
    {
      output << "# Observed Variables" << endl
             << "options_.obs_vars = [" << endl;
      for (vector<int>::const_iterator it = varobs.begin();
           it != varobs.end(); it++)
        output << "                    DynareModel.ObsVars("
               << getTypeSpecificID(*it)+1 << ")" << endl;
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
