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

using namespace std;

#include <iostream>

#include "Shocks.hh"
#include "Interface.hh"

AbstractShocksStatement::AbstractShocksStatement(bool mshocks_arg,
                                                  const det_shocks_type &det_shocks_arg,
                                                  const var_and_std_shocks_type &var_shocks_arg,
                                                  const var_and_std_shocks_type &std_shocks_arg,
                                                  const covar_and_corr_shocks_type &covar_shocks_arg,
                                                  const covar_and_corr_shocks_type &corr_shocks_arg,
                                                  const SymbolTable &symbol_table_arg) :
  mshocks(mshocks_arg),
  det_shocks(det_shocks_arg),
  var_shocks(var_shocks_arg),
  std_shocks(std_shocks_arg),
  covar_shocks(covar_shocks_arg),
  corr_shocks(corr_shocks_arg),
  symbol_table(symbol_table_arg)
{
}

void
AbstractShocksStatement::writeDetShocks(ostream &output) const
{
  int exo_det_length = 0;

  for(det_shocks_type::const_iterator it = det_shocks.begin();
      it != det_shocks.end(); it++)
    {
      int id = symbol_table.getID(it->first) + 1;
      bool exo_det = (symbol_table.getType(it->first) == eExogenousDet);
      int set_shocks_index = ((int) mshocks) + 2 * ((int) exo_det);

      for (unsigned int i = 0; i < it->second.size(); i++)
        {
          const int &period1 = it->second[i].period1;
          const int &period2 = it->second[i].period2;
          const NodeID value = it->second[i].value;

          if (period1 == period2)
            {
              output << "set_shocks(" << set_shocks_index << "," << period1
                     << ", " << id << ", ";
              value->writeOutput(output);
              output << ");" << endl;
            }
          else
            {
              output << "set_shocks(" << set_shocks_index << "," << period1
                     << ":" << period2 << ", " << id << ", ";
              value->writeOutput(output);
              output << ");" << endl;
            }

          if (exo_det && (period2 > exo_det_length))
            exo_det_length = period2;
        }
    }
  output << "M_.exo_det_length = " << exo_det_length << ";\n";
}

void
AbstractShocksStatement::writeVarAndStdShocks(ostream &output) const
{
  var_and_std_shocks_type::const_iterator it;

  for(it = var_shocks.begin(); it != var_shocks.end(); it++)
    {
      int id = symbol_table.getID(it->first) + 1;
      const NodeID value = it->second;
      output << "M_.Sigma_e(" << id << ", " << id << ") = ";
      value->writeOutput(output);
      output << ";" << endl;
    }

  for(it = std_shocks.begin(); it != std_shocks.end(); it++)
    {
      int id = symbol_table.getID(it->first) + 1;
      const NodeID value = it->second;
      output << "M_.Sigma_e(" << id << ", " << id << ") = (";
      value->writeOutput(output);
      output << ")^2;" << endl;
    }
}

void
AbstractShocksStatement::writeCovarAndCorrShocks(ostream &output) const
{
  covar_and_corr_shocks_type::const_iterator it;

  for(it = covar_shocks.begin(); it != covar_shocks.end(); it++)
    {
      int id1 = symbol_table.getID(it->first.first) + 1;
      int id2 = symbol_table.getID(it->first.second) + 1;
      const NodeID value = it->second;
      output << "M_.Sigma_e(" << id1 << ", " << id2 << ") = ";
      value->writeOutput(output);
      output << "; M_.Sigma_e(" << id2 << ", " << id1 << ") = M_.Sigma_e("
             << id1 << ", " << id2 << ");\n";
    }

  for(it = corr_shocks.begin(); it != corr_shocks.end(); it++)
    {
      int id1 = symbol_table.getID(it->first.first) + 1;
      int id2 = symbol_table.getID(it->first.second) + 1;
      const NodeID value = it->second;
      output << "M_.Sigma_e(" << id1 << ", " << id2 << ") = ";
      value->writeOutput(output);
      output << "*sqrt(M_.Sigma_e(" << id1 << ", " << id1 << ")*M_.Sigma_e("
             << id2 << ", " << id2 << "); M_.Sigma_e(" << id2 << ", "
             << id1 << ") = M_.Sigma_e(" << id1 << ", " << id2 << ");\n";
    }
}


ShocksStatement::ShocksStatement(const det_shocks_type &det_shocks_arg,
                                 const var_and_std_shocks_type &var_shocks_arg,
                                 const var_and_std_shocks_type &std_shocks_arg,
                                 const covar_and_corr_shocks_type &covar_shocks_arg,
                                 const covar_and_corr_shocks_type &corr_shocks_arg,
                                 const SymbolTable &symbol_table_arg) :
  AbstractShocksStatement(false, det_shocks_arg, var_shocks_arg, std_shocks_arg,
                          covar_shocks_arg, corr_shocks_arg, symbol_table_arg)
{
}

void
ShocksStatement::writeOutput(ostream &output, const string &basename) const
{
  output << interfaces::comment() << endl << interfaces::comment();
  output << "SHOCKS instructions \n";
  output << interfaces::comment() << "\n";

  // Write instruction that initializes a shock
  output << "make_ex_;\n";

  writeDetShocks(output);
  writeVarAndStdShocks(output);
  writeCovarAndCorrShocks(output);
}

MShocksStatement::MShocksStatement(const det_shocks_type &det_shocks_arg,
                                   const var_and_std_shocks_type &var_shocks_arg,
                                   const var_and_std_shocks_type &std_shocks_arg,
                                   const covar_and_corr_shocks_type &covar_shocks_arg,
                                   const covar_and_corr_shocks_type &corr_shocks_arg,
                                   const SymbolTable &symbol_table_arg) :
  AbstractShocksStatement(true, det_shocks_arg, var_shocks_arg, std_shocks_arg,
                          covar_shocks_arg, corr_shocks_arg, symbol_table_arg)
{
}

void
MShocksStatement::writeOutput(ostream &output, const string &basename) const
{
  output << interfaces::comment() << endl << interfaces::comment();
  output << "SHOCKS instructions \n";
  output << interfaces::comment() << "\n";

  // Write instruction that initializes a shock
  output << "make_ex_;\n";

  writeDetShocks(output);
  writeVarAndStdShocks(output);
  writeCovarAndCorrShocks(output);
}
