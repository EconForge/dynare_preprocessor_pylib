/*
 * Copyright (C) 2003-2015 Dynare Team
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

#include <utility>

#include "SigmaeInitialization.hh"

SigmaeStatement::SigmaeStatement(matrix_t matrix_arg) noexcept(false) :
  matrix(move(matrix_arg)),
  matrix_form(determineMatrixForm(matrix))
{
}

SigmaeStatement::MatrixForm
SigmaeStatement::determineMatrixForm(const matrix_t &matrix) noexcept(false)
{
  size_t nbe;
  int inc;
  MatrixForm type;
  // Checking if first or last row has one element.
  if (matrix.front().size() == 1)
    {
      inc = 1;
      nbe = 2;
      type = MatrixForm::lower;
    }
  else if (matrix.back().size() == 1)
    {
      inc = -1;
      nbe = matrix.front().size()-1;
      type = MatrixForm::upper;
    }
  else
    throw MatrixFormException();

  // Checking if matrix is triangular (upper or lower):
  // each row has one element more or less than the previous one
  // and first or last one has one element.
  matrix_t::const_iterator ir;
  for (ir = matrix.begin(), ir++; ir != matrix.end(); ir++, nbe += inc)
    if (ir->size() != nbe)
      throw MatrixFormException();

  return type;
}

void
SigmaeStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  size_t ic, ic1, ir, ir1;

  output << "M_.Sigma_e = [..." << endl;
  for (ir = 0; ir < matrix.size(); ir++)
    {
      for (ic = 0; ic < matrix.size(); ic++)
        {
          if (ic >= ir && matrix_form == MatrixForm::upper)
            {
              ic1 = ic-ir;
              ir1 = ir;
            }
          else if (ic < ir && matrix_form == MatrixForm::upper)
            {
              ic1 = ir-ic;
              ir1 = ic;
            }
          else if (ic > ir && matrix_form == MatrixForm::lower)
            {
              ic1 = ir;
              ir1 = ic;
            }
          else // ic <= ir && matrix_form == eLower
            {
              ic1 = ic;
              ir1 = ir;
            }

          matrix[ir1][ic1]->writeOutput(output);
          output << " ";
        }
      output << ";..." << endl;
    }
  output << "];" << endl;
}
