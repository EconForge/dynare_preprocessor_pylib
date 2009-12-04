/*
 * Copyright (C) 2006-2009 Dynare Team
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

#ifndef _STATEMENT_HH
#define _STATEMENT_HH

using namespace std;

#include <ostream>
#include <string>
#include <map>

#include "SymbolList.hh"

class ModFileStructure
{
public:
  ModFileStructure();
  //! Wheter check is present
  bool check_present;
  //! Whether a simul statement is present
  bool simul_present;
  //! Whether a stoch_simul statement is present
  bool stoch_simul_present;
  //! Whether an estimation statement is present
  bool estimation_present;
  //! Whether an osr statement is present
  bool osr_present;
  //! Whether a ramsey_policy statement is present
  bool ramsey_policy_present;
  //! The value of the "order" option of stoch_simul, estimation, osr, ramsey_policy
  //! Derivation order
  /*! First initialized to zero. If user sets order option somewhere in the MOD file, it will be equal to the maximum of order options. Otherwise will default to 2 */
  int order_option;
  //! Whether a bvar_density, bvar_forecast, sbvar, ms_sbvar statement is present
  bool bvar_present;
  //! Whether an svar_identification statement is present
  bool svar_identification_present;
  //! Whether an identification statement is present or the identification option of dynare_sensitivity statement is equal to one
  bool identification_present;
  //! Whether the option partial_information is given to stoch_simul/estimation/osr/ramsey_policy
  bool partial_information;
  //! Whether a shocks or mshocks block is present
  /*! Used for the workaround for trac ticket #35 */
  bool shocks_present;
};

class Statement
{
public:
  virtual ~Statement();
  //! Do some internal check, and fill the ModFileStructure class
  virtual void checkPass(ModFileStructure &mod_file_struct);
  virtual void computingPass();
  //! Write Matlab output code
  /*!
    \param output is the output stream of the main matlab file
    \param basename is the name of the modfile (without extension) which can be used to build auxiliary files
  */
  virtual void writeOutput(ostream &output, const string &basename) const = 0;
};

class NativeStatement : public Statement
{
private:
  const string native_statement;
public:
  NativeStatement(const string &native_statement_arg);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

class OptionsList
{
public:
  typedef map<string, string> num_options_type;
  typedef map<string, pair<string, string> > paired_num_options_type;
  typedef map<string, string> string_options_type;
  typedef map<string, SymbolList> symbol_list_options_type;
  num_options_type num_options;
  paired_num_options_type paired_num_options;
  string_options_type string_options;
  symbol_list_options_type symbol_list_options;
  void writeOutput(ostream &output) const;
  void writeOutput(ostream &output, const string &option_group) const;
  void clear();
};

#endif // ! _STATEMENT_HH
