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

#ifndef _NUMERICALCONSTANTS_HH
#define _NUMERICALCONSTANTS_HH

using namespace std;

#include <string>
#include <vector>
#include <map>

//! Handles numerical constants
class NumericalConstants
{
private:
  //! Vector of numerical constants
  vector<string> mNumericalConstants;
  //! Map matching constants to their id
  map<string, int> numConstantsIndex;
public:
  NumericalConstants();
  //! Adds a constant and returns its ID
  int AddConstant(const string &iConst);
  //! Get a constant in string form
  string get(int iID) const;
  //! Get a constant in double form
  double getDouble(int iID) const;
};

#endif
