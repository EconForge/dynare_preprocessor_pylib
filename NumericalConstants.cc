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

#include "NumericalConstants.hh"

NumericalConstants::NumericalConstants()
{
  AddConstant("0");
  AddConstant("1");
}

int
NumericalConstants::AddConstant(const string &iConst)
{
  map<string, int>::iterator iter = numConstantsIndex.find(iConst);
  //cout << "iConst=" << iConst << "\n" ;
  if (iter != numConstantsIndex.end())
    return iter->second;

  if (atof(iConst.c_str()) < 0)
    {
      cerr << "Can't handle a negative constant..!" << endl;
      exit(-1);
    }

  int id = (int) mNumericalConstants.size();
  mNumericalConstants.push_back(iConst);
  numConstantsIndex[iConst] = id;
  return id;
}

string
NumericalConstants::get(int ID) const
{
  if (ID < (int) mNumericalConstants.size())
    return mNumericalConstants[ID];
  else
    {
      cerr << "Unknown constant" << endl;
      exit(-1);
    }
}

double
NumericalConstants::getDouble(int iID) const
{
  return(atof(get(iID).c_str()));
}
