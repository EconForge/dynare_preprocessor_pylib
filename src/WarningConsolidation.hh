/*
 * Copyright © 2012-2024 Dynare Team
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

#ifndef WARNING_CONSOLIDATION_HH
#define WARNING_CONSOLIDATION_HH

#include "DynareBisonLocation.hh"

#include <iostream>
#include <ostream>
#include <sstream>
#include <string>

using namespace std;

/* Provide our implementation of operator<< with locations in DynareBison.hh. Note that the
   following is a template specialization of the version provided in DynareBisonLocation.hh.

   Ideally it should go into DynareBisonLocation.hh, but there does not seem to be a way to achieve
   that. */
ostream& operator<<(ostream& stream, const Dynare::location& l);

class WarningConsolidation
{
private:
  const bool no_warn;
  int num_warnings {0};

  // Increases the warning counter by as many newlines as there are in the message
  void incrementWarnings(const string& msg);

public:
  explicit WarningConsolidation(bool no_warn_arg) : no_warn {no_warn_arg}
  {
  }

  // Generic function to print something to the warning stream
  friend WarningConsolidation& operator<<(WarningConsolidation& wcc, auto&& warning);

  /* Print std::endl to the warning stream. Unfortunately, since std::endl is a template of
     functions, it cannot be bound to the universal reference of the generic function, hence the
     need for this specialization. */
  friend WarningConsolidation& operator<<(WarningConsolidation& wcc, ostream& (*pf)(ostream&));

  int
  numWarnings() const
  {
    return num_warnings;
  };
};

WarningConsolidation&
operator<<(WarningConsolidation& wcc, auto&& warning)
{
  if (wcc.no_warn)
    return wcc;

  ostringstream ostr;
  ostr << warning;

  cerr << ostr.str();

  wcc.incrementWarnings(ostr.str());

  return wcc;
}

#endif
