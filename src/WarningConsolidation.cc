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

#include "WarningConsolidation.hh"

ostream&
operator<<(ostream& stream, const Dynare::location& l)
{
  stream << *l.begin.filename << ": line " << l.begin.line;
  if (l.begin.line == l.end.line)
    if (l.begin.column == l.end.column - 1)
      stream << ", col " << l.begin.column;
    else
      stream << ", cols " << l.begin.column << "-" << l.end.column - 1;
  else
    stream << ", col " << l.begin.column << " -"
           << " line " << l.end.line << ", col " << l.end.column - 1;

  return stream;
}

void
WarningConsolidation::incrementWarnings(const string& msg)
{
  size_t p {0};
  while ((p = msg.find('\n', p)) != string::npos)
    {
      p++;
      num_warnings++;
    }
}

WarningConsolidation&
operator<<(WarningConsolidation& wcc, ostream& (*pf)(ostream&))
{
  if (wcc.no_warn)
    return wcc;

  ostringstream ostr;
  ostr << pf;

  cerr << ostr.str();

  wcc.incrementWarnings(ostr.str());

  return wcc;
}
