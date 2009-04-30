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

#include "SymbolTable.hh"

SymbolTable::SymbolTable() : frozen(false), size(0)
{
}

void
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
}

void
SymbolTable::addSymbol(const string &name, SymbolType type) throw (AlreadyDeclaredException, FrozenException)
{
  // Construct "tex_name" by prepending an antislash to all underscores in "name"
  string tex_name = name;
  size_t pos = 0;
  while((pos = tex_name.find('_', pos)) != string::npos)
    {
      tex_name.insert(pos, "\\");
      pos += 2;
    }
  addSymbol(name, type, tex_name);
}

void
SymbolTable::freeze() throw (FrozenException)
{
  if (frozen)
    throw FrozenException();

  frozen = true;

  for(int i = 0; i < size; i++)
    {
      int tsi;
      switch(getType(i))
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

  switch(type)
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
}
