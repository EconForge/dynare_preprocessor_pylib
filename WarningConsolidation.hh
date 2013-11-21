/*
 * Copyright (C) 2012-2013 Dynare Team
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

#ifndef _WARNINGCONSOLIDATION_HH
#define _WARNINGCONSOLIDATION_HH

#include <sstream>
#include <string>
#include "location.hh"

using namespace std;

//! Stores Warnings issued by the Preprocessor
class WarningConsolidation
{
private:
  stringstream warnings;
  bool no_warn;

public:
  WarningConsolidation(bool no_warn_arg) : no_warn(no_warn_arg) { };
  ~WarningConsolidation() { };

  //! Add A Warning to the StringStream
  friend WarningConsolidation& operator<< (WarningConsolidation& wcc, const string &warning);
  friend WarningConsolidation& operator<< (WarningConsolidation& wcc, const Dynare::location &loc);
  friend WarningConsolidation& operator<< (WarningConsolidation& wcc, ostream& (*pf) (ostream&));

  inline void addWarning(const string &w) { warnings << w; };
  inline void addWarning(ostream& (*pf) (ostream&)) { warnings << pf; };

  //! Write Warnings to m file
  void writeOutput(ostream &output) const;
  //! Count warnings
  /*! This is done in a very lousy way, by counting newlines in the
    stringstream... */
  int countWarnings() const;
};

#endif
