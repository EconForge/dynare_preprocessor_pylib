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

#include <algorithm>
#include <sstream>
#include <iostream>
#include <cassert>

#include "SymbolTable.hh"

SymbolTable::SymbolTable() : frozen(false), size(0)
{
}

int
SymbolTable::addSymbol(const string &name, SymbolType type, const string &tex_name) throw (AlreadyDeclaredException, FrozenException)
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

  int id = size++;

  symbol_table[name] = id;
  type_table.push_back(type);
  name_table.push_back(name);
  tex_name_table.push_back(tex_name);

  return id;
}

int
SymbolTable::addSymbol(const string &name, SymbolType type) throw (AlreadyDeclaredException, FrozenException)
{
  // Construct "tex_name" by prepending an antislash to all underscores in "name"
  string tex_name = name;
  size_t pos = 0;
  while ((pos = tex_name.find('_', pos)) != string::npos)
    {
      tex_name.insert(pos, "\\");
      pos += 2;
    }
  return addSymbol(name, type, tex_name);
}

void
SymbolTable::freeze() throw (FrozenException)
{
  if (frozen)
    throw FrozenException();

  frozen = true;

  for (int i = 0; i < size; i++)
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
SymbolTable::changeType(int id, SymbolType newtype) throw (UnknownSymbolIDException, FrozenException)
{
  if (frozen)
    throw FrozenException();

  if (id < 0 || id >= size)
    throw UnknownSymbolIDException(id);

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

void
SymbolTable::writeOutput(ostream &output) const throw (NotYetFrozenException)
{
  if (!frozen)
    throw NotYetFrozenException();

  if (exo_nbr() > 0)
    {
      output << "M_.exo_names = '" << getName(exo_ids[0]) << "';" << endl;
      output << "M_.exo_names_tex = '" << getTeXName(exo_ids[0]) << "';" << endl;
      for (int id = 1; id < exo_nbr(); id++)
        {
          output << "M_.exo_names = strvcat(M_.exo_names, '" << getName(exo_ids[id]) << "');" << endl
                 << "M_.exo_names_tex = strvcat(M_.exo_names_tex, '" << getTeXName(exo_ids[id]) << "');" << endl;
        }
    }
  if (exo_det_nbr() > 0)
    {
      output << "M_.exo_det_names = '" << getName(exo_det_ids[0]) << "';" << endl;
      output << "M_.exo_det_names_tex = '" << getTeXName(exo_det_ids[0]) << "';" << endl;
      for (int id = 1; id < exo_det_nbr(); id++)
        {
          output << "M_.exo_det_names = strvcat(M_.exo_det_names, '" << getName(exo_det_ids[id]) << "');" << endl
                 << "M_.exo_det_names_tex = strvcat(M_.exo_det_names_tex, '" << getTeXName(exo_det_ids[id]) << "');" << endl;
        }
    }
  if (endo_nbr() > 0)
    {
      output << "M_.endo_names = '" << getName(endo_ids[0]) << "';" << endl;
      output << "M_.endo_names_tex = '" << getTeXName(endo_ids[0]) << "';" << endl;
      for (int id = 1; id < endo_nbr(); id++)
        {
          output << "M_.endo_names = strvcat(M_.endo_names, '" << getName(endo_ids[id]) << "');" << endl
                 << "M_.endo_names_tex = strvcat(M_.endo_names_tex, '" << getTeXName(endo_ids[id]) << "');" << endl;
        }
    }
  if (param_nbr() > 0)
    {
      output << "M_.param_names = '" << getName(param_ids[0]) << "';" << endl;
      output << "M_.param_names_tex = '" << getTeXName(param_ids[0]) << "';" << endl;
      for (int id = 1; id < param_nbr(); id++)
        {
          output << "M_.param_names = strvcat(M_.param_names, '" << getName(param_ids[id]) << "');" << endl
                 << "M_.param_names_tex = strvcat(M_.param_names_tex, '" << getTeXName(param_ids[id]) << "');" << endl;
        }
    }

  output << "M_.exo_det_nbr = " << exo_det_nbr() << ";" << endl
         << "M_.exo_nbr = " << exo_nbr() << ";" << endl
         << "M_.endo_nbr = " << endo_nbr() << ";" << endl
         << "M_.param_nbr = " << param_nbr() << ";" << endl;

  output << "M_.Sigma_e = zeros(" << exo_nbr() << ", " << exo_nbr() << ");" << endl;

  // Write the auxiliary variable table
  output << "M_.orig_endo_nbr = " << endo_nbr() - aux_vars.size() << ";" << endl;
  if (aux_vars.size() == 0)
    output << "M_.aux_vars = [];" << endl;
  else
    for (int i = 0; i < (int) aux_vars.size(); i++)
      {
        output << "M_.aux_vars(" << i+1 << ").endo_index = " << getTypeSpecificID(aux_vars[i].symb_id)+1 << ";" << endl
               << "M_.aux_vars(" << i+1 << ").type = " << aux_vars[i].type << ";" << endl;
        switch (aux_vars[i].type)
          {
          case avEndoLead:
          case avExoLead:
          case avExpectation:
            break;
          case avEndoLag:
          case avExoLag:
            output << "M_.aux_vars(" << i+1 << ").orig_index = " << getTypeSpecificID(aux_vars[i].orig_symb_id)+1 << ";" << endl
                   << "M_.aux_vars(" << i+1 << ").orig_lead_lag = " << aux_vars[i].orig_lead_lag << ";" << endl;
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
}

int
SymbolTable::addLeadAuxiliaryVarInternal(bool endo, int index) throw (FrozenException)
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

  AuxVarInfo avi;
  avi.symb_id = symb_id;
  avi.type = (endo ? avEndoLead : avExoLead);
  aux_vars.push_back(avi);

  return symb_id;
}

int
SymbolTable::addLagAuxiliaryVarInternal(bool endo, int orig_symb_id, int orig_lead_lag) throw (FrozenException)
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

  AuxVarInfo avi;
  avi.symb_id = symb_id;
  avi.type = (endo ? avEndoLag : avExoLag);
  avi.orig_symb_id = orig_symb_id;
  avi.orig_lead_lag = orig_lead_lag;
  aux_vars.push_back(avi);

  return symb_id;
}

int
SymbolTable::addEndoLeadAuxiliaryVar(int index) throw (FrozenException)
{
  return addLeadAuxiliaryVarInternal(true, index);
}

int
SymbolTable::addEndoLagAuxiliaryVar(int orig_symb_id, int orig_lead_lag) throw (FrozenException)
{
  return addLagAuxiliaryVarInternal(true, orig_symb_id, orig_lead_lag);
}

int
SymbolTable::addExoLeadAuxiliaryVar(int index) throw (FrozenException)
{
  return addLeadAuxiliaryVarInternal(false, index);
}

int
SymbolTable::addExoLagAuxiliaryVar(int orig_symb_id, int orig_lead_lag) throw (FrozenException)
{
  return addLagAuxiliaryVarInternal(false, orig_symb_id, orig_lead_lag);
}

int
SymbolTable::addExpectationAuxiliaryVar(int information_set, int index) throw (FrozenException)
{
  ostringstream varname;
  varname << "AUX_EXPECT_" << (information_set < 0 ? "LAG" : "LEAD") << "_"
          << abs(information_set) << "_" << index;
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

  AuxVarInfo avi;
  avi.symb_id = symb_id;
  avi.type = avExpectation;
  aux_vars.push_back(avi);

  return symb_id;
}

void
SymbolTable::markPredetermined(int symb_id) throw (UnknownSymbolIDException, FrozenException)
{
  if (symb_id < 0 || symb_id >= size)
    throw UnknownSymbolIDException(symb_id);
  if (frozen)
    throw FrozenException();

  assert(getType(symb_id) == eEndogenous);

  predetermined_variables.insert(symb_id);
}

bool
SymbolTable::isPredetermined(int symb_id) const throw (UnknownSymbolIDException)
{
  if (symb_id < 0 || symb_id >= size)
    throw UnknownSymbolIDException(symb_id);

  return (predetermined_variables.find(symb_id) != predetermined_variables.end());
}

int
SymbolTable::predeterminedNbr() const
{
  return (predetermined_variables.size());
}
