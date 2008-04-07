/*
 * Copyright (C) 2006-2008 Dynare Team
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

#include "Statement.hh"

ModFileStructure::ModFileStructure() :
  check_present(false),
  simul_present(false),
  stoch_simul_or_similar_present(false),
  order_option(2)
{
}

Statement::~Statement()
{
}

void
Statement::checkPass(ModFileStructure &mod_file_struct)
{
}

void
Statement::computingPass()
{
}

NativeStatement::NativeStatement(const string &native_statement_arg) :
  native_statement(native_statement_arg)
{
}

void
NativeStatement::writeOutput(ostream &output, const string &basename) const
{
  output << native_statement << endl;
}

void
OptionsList::writeOutput(ostream &output) const
{
  for(num_options_type::const_iterator it = num_options.begin();
      it != num_options.end(); it++)
    output << "options_." << it->first << " = " << it->second << ";" << endl;

  for(paired_num_options_type::const_iterator it = paired_num_options.begin();
      it != paired_num_options.end(); it++)
    output << "options_." << it->first << " = [" << it->second.first << "; "
           << it->second.second << "];" << endl;

  for(string_options_type::const_iterator it = string_options.begin();
      it != string_options.end(); it++)
    output << "options_." << it->first << " = '" << it->second << "';" << endl;

  for(symbol_list_options_type::const_iterator it = symbol_list_options.begin();
      it != symbol_list_options.end(); it++)
    it->second.writeOutput("options_." + it->first, output);
}

void
OptionsList::writeOutput(ostream &output, const string &option_group) const
{
  output << option_group << " = struct();" << endl;

  for(num_options_type::const_iterator it = num_options.begin();
      it != num_options.end(); it++)
    output << option_group << "." << it->first << " = " << it->second << ";" << endl;

  for(paired_num_options_type::const_iterator it = paired_num_options.begin();
      it != paired_num_options.end(); it++)
    output << option_group << "." << it->first << " = [" << it->second.first << "; "
           << it->second.second << "];" << endl;

  for(string_options_type::const_iterator it = string_options.begin();
      it != string_options.end(); it++)
    output << option_group << "." << it->first << " = '" << it->second << "';" << endl;

  for(symbol_list_options_type::const_iterator it = symbol_list_options.begin();
      it != symbol_list_options.end(); it++)
    it->second.writeOutput(option_group + "." + it->first, output);
}

void
OptionsList::clear()
{
  num_options.clear();
  paired_num_options.clear();
  string_options.clear();
  symbol_list_options.clear();
}
