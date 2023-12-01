/*
 * Copyright © 2003-2023 Dynare Team
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
 * along with Dynare.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef NUMERICAL_CONSTANTS_HH
#define NUMERICAL_CONSTANTS_HH

#include <map>
#include <string>
#include <vector>

using namespace std;

//! Handles non-negative numerical constants
class NumericalConstants
{
private:
  //! Vector of numerical constants
  vector<string> mNumericalConstants;
  //! Double values of these constants
  vector<double> double_vals;
  //! Map matching constants to their id
  map<string, int> numConstantsIndex;

public:
  //! Adds a non-negative constant (possibly Inf or NaN) and returns its ID
  int AddNonNegativeConstant(const string& iConst);
  //! Get a constant in string form
  [[nodiscard]] string get(int ID) const;
  //! Get a constant in double form
  [[nodiscard]] double getDouble(int ID) const;
};

#endif
