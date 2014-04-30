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

#include "WarningConsolidation.hh"
#include <ostream>

WarningConsolidation&
operator<< (WarningConsolidation& wcc, const string &warning)
{
  if (wcc.no_warn)
    return wcc;

  cerr << warning;
  wcc.addWarning(warning);
  return wcc;
};

WarningConsolidation&
operator<< (WarningConsolidation& wcc, const Dynare::location& loc)
{
  if (wcc.no_warn)
    return wcc;

  stringstream ostr;
  Dynare::position last = loc.end - 1;
  ostr << loc.begin;
  if (last.filename
      && (!loc.begin.filename
          || *loc.begin.filename != *last.filename))
    ostr << '-' << last;
  else if (loc.begin.line != last.line)
    ostr << '-' << last.line  << '.' << last.column;
  else if (loc.begin.column != last.column)
    ostr << '-' << last.column;

  cerr << ostr.str();
  wcc.addWarning(ostr.str());
  return wcc;
};

WarningConsolidation&
operator<< (WarningConsolidation& wcc, ostream& (*pf) (ostream&))
{
  if (wcc.no_warn)
    return wcc;

  cerr << pf;
  wcc.addWarning(pf);
  return wcc;
}

void
WarningConsolidation::writeOutput(ostream &output) const
{
  if (warnings.str().empty())
    return;

  output << "disp([char(10) 'Dynare Preprocessor Warning(s) Encountered:']);" << endl;

  bool writedisp = true;
  string warningsstr = warnings.str();
  for (size_t i = 0; i < warningsstr.length(); i++)
    {
      if (writedisp)
        {
          output << "disp('     ";
          writedisp = false;
        }

      if (warningsstr[i] != '\n')
        output << warningsstr[i];
      else
        {
          output << "');" << endl;
          if (i+1 < warningsstr.length())
            writedisp = true;
        }
    }
}

int
WarningConsolidation::countWarnings() const
{
  size_t p = 0;
  int n = 0;
  while ((p = warnings.str().find('\n', p)) != string::npos)
    {
      p++;
      n++;
    }
  return n;
}
