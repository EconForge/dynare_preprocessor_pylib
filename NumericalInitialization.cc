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

#include "NumericalInitialization.hh"

InitParamStatement::InitParamStatement(const string &param_name_arg,
                                       const NodeID param_value_arg,
                                       const SymbolTable &symbol_table_arg) :
  param_name(param_name_arg),
  param_value(param_value_arg),
  symbol_table(symbol_table_arg)
{
}

void
InitParamStatement::writeOutput(ostream &output, const string &basename) const
{
  int id = symbol_table.getID(param_name) + 1;
  output << "M_.params( " << id << " ) = ";
  param_value->writeOutput(output);
  output << ";" << endl;
  output << param_name << " = M_.params( " << id << " );\n";
}

InitOrEndValStatement::InitOrEndValStatement(const init_values_type &init_values_arg,
                                             const SymbolTable &symbol_table_arg) :
  init_values(init_values_arg),
  symbol_table(symbol_table_arg)
{
}

void
InitOrEndValStatement::writeInitValues(ostream &output) const
{
  for(init_values_type::const_iterator it = init_values.begin();
      it != init_values.end(); it++)
    {
      const string &name = it->first;
      const NodeID expression = it->second;

      SymbolType type = symbol_table.getType(name);
      int id = symbol_table.getID(name) + 1;

      if (type == eEndogenous)
        output << "oo_.steady_state";
      else if (type == eExogenous)
        output << "oo_.exo_steady_state";
      else if (type == eExogenousDet)
        output << "oo_.exo_det_steady_state";

      output << "( " << id << " ) = ";
      expression->writeOutput(output);
      output << ";" << endl;
    }
}

InitValStatement::InitValStatement(const init_values_type &init_values_arg,
                                   const SymbolTable &symbol_table_arg) :
  InitOrEndValStatement(init_values_arg, symbol_table_arg)
{
}

void
InitValStatement::writeOutput(ostream &output, const string &basename) const
{
  output << "%" << endl
         << "% INITVAL instructions" << endl
         << "%" << endl;
  // Writing initval block to set initial values for variables
  output << "options_.initval_file = 0;" << endl
         << "endval_=0;" << endl;

  if (symbol_table.recur_nbr > 0)
    output << "recurs_ = zeros(" << symbol_table.recur_nbr << ", 1);\n";

  writeInitValues(output);

  output << "oo_.endo_simul=[oo_.steady_state*ones(1,M_.maximum_lag)];\n";
  output << "if M_.exo_nbr > 0;\n";
  output << "\too_.exo_simul = [ones(M_.maximum_lag,1)*oo_.exo_steady_state'];\n";
  output <<"end;\n";
  output << "if M_.exo_det_nbr > 0;\n";
  output << "\too_.exo_det_simul = [ones(M_.maximum_lag,1)*oo_.exo_det_steady_state'];\n";
  output <<"end;\n";
}


EndValStatement::EndValStatement(const init_values_type &init_values_arg,
                                 const SymbolTable &symbol_table_arg) :
  InitOrEndValStatement(init_values_arg, symbol_table_arg)
{
}


void
EndValStatement::writeOutput(ostream &output, const string &basename) const
{
  output << "%" << endl
         << "% ENDVAL instructions" << endl
         << "%" << endl;
  // Writing endval block to set terminal values for variables
  output << "ys0_= oo_.steady_state;" << endl
         << "ex0_ = oo_.exo_steady_state;" << endl
         << "recurs0_ = recurs_;" << endl
         << "endval_ = 1;" << endl;

  writeInitValues(output);
}

HistValStatement::HistValStatement(const hist_values_type &hist_values_arg,
                                   const SymbolTable &symbol_table_arg) :
  hist_values(hist_values_arg),
  symbol_table(symbol_table_arg)
{
}

void
HistValStatement::writeOutput(ostream &output, const string &basename) const
{
  output << "%" << endl
         << "% HISTVAL instructions" << endl
         << "%" << endl;

  for(hist_values_type::const_iterator it = hist_values.begin();
      it != hist_values.end(); it++)
    {
      const string &name = it->first.first;
      const int &lag = it->first.second;
      const NodeID expression = it->second;

      SymbolType type = symbol_table.getType(name);
      int id = symbol_table.getID(name) + 1;

      if (type == eEndogenous)
        output << "oo_.endo_simul( " << id << ", M_.maximum_lag + " << lag << ") = ";
      else if (type == eExogenous)
        output << "oo_.exo_simul( M_.maximum_lag + " << lag << ", " << id << " ) = ";
      else if (type != eExogenousDet)
        output << "oo_.exo_det_simul( M_.maximum_lag + " << lag  << ", " << id << " ) = ";

      expression->writeOutput(output);
      output << ";" << endl;
    }
}

InitvalFileStatement::InitvalFileStatement(const string &filename_arg) :
  filename(filename_arg)
{
}

void
InitvalFileStatement::writeOutput(ostream &output, const string &basename) const
{
  output << "%" << endl
         << "% INITVAL_FILE statement" << endl
         << "%" << endl
         << "options_.initval_file = 1;" << endl
         << "initvalf('" << filename << "');" << endl;
}

HomotopyStatement::HomotopyStatement(const homotopy_values_type &homotopy_values_arg,
                                     const SymbolTable &symbol_table_arg) :
  homotopy_values(homotopy_values_arg),
  symbol_table(symbol_table_arg)
{
}

void
HomotopyStatement::writeOutput(ostream &output, const string &basename) const
{
  output << "%" << endl
         << "% HOMOTOPY_SETUP instructions" << endl
         << "%" << endl
         << "options_.homotopy_values = [];" << endl;

  for(homotopy_values_type::const_iterator it = homotopy_values.begin();
      it != homotopy_values.end(); it++)
    {
      const string &name = it->first;
      const NodeID expression1 = it->second.first;
      const NodeID expression2 = it->second.second;

      const SymbolType type = symbol_table.getType(name);
      const int id = symbol_table.getID(name) + 1;

      output << "options_.homotopy_values = vertcat(options_.homotopy_values, [ " << type << ", " << id << ", ";
      if (expression1 != NULL)
        expression1->writeOutput(output);
      else
        output << "NaN";
      output << ", ";
      expression2->writeOutput(output);
      output << "]);" << endl;
    }
}
