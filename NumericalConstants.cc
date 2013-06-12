/*
 * Copyright (C) 2003-2012 Dynare Team
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
#include <cassert>
#include <cerrno>
#include <cmath>
#include <iostream>

#include "NumericalConstants.hh"

int
NumericalConstants::AddNonNegativeConstant(const string &iConst)
{
  map<string, int>::const_iterator iter = numConstantsIndex.find(iConst);

  if (iter != numConstantsIndex.end())
    return iter->second;

  int id = (int) mNumericalConstants.size();
  mNumericalConstants.push_back(iConst);
  numConstantsIndex[iConst] = id;

  double val = strtod(iConst.c_str(), NULL);

  /* Note that we allow underflows (will be converted to 0) and overflows (will
     be converted to Inf), as MATLAB and Octave do. */

  assert(val >= 0 || isnan(val)); // Check we have a positive constant or a NaN
  double_vals.push_back(val);

  return id;
}

string
NumericalConstants::get(int ID) const
{
  assert(ID >= 0 && ID < (int) mNumericalConstants.size());
  return mNumericalConstants[ID];
}

double
NumericalConstants::getDouble(int ID) const
{
  assert(ID >= 0 && ID < (int) double_vals.size());
  return (double_vals[ID]);
}
