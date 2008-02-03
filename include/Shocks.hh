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

#ifndef _SHOCKS_HH
#define _SHOCKS_HH

using namespace std;

#include <string>
#include <vector>
#include <map>

#include "Statement.hh"
#include "SymbolTable.hh"
#include "ExprNode.hh"

class AbstractShocksStatement : public Statement
{
public:
  struct DetShockElement
  {
    int period1;
    int period2;
    NodeID value;
  };
  typedef map<string, vector<DetShockElement> > det_shocks_type;
  typedef map<string, NodeID> var_and_std_shocks_type;
  typedef map<pair<string, string>, NodeID> covar_and_corr_shocks_type;
protected:
  //! Is this statement a "mshocks" statement ? (instead of a "shocks" statement)
  const bool mshocks;
  const det_shocks_type det_shocks;
  const var_and_std_shocks_type var_shocks, std_shocks;
  const covar_and_corr_shocks_type covar_shocks, corr_shocks;
  const SymbolTable &symbol_table;
  void writeDetShocks(ostream &output) const;
  void writeVarAndStdShocks(ostream &output) const;
  void writeCovarAndCorrShocks(ostream &output) const;

  AbstractShocksStatement(bool mshocks_arg,
                          const det_shocks_type &det_shocks_arg,
                          const var_and_std_shocks_type &var_shocks_arg,
                          const var_and_std_shocks_type &std_shocks_arg,
                          const covar_and_corr_shocks_type &covar_shocks_arg,
                          const covar_and_corr_shocks_type &corr_shocks_arg,
                          const SymbolTable &symbol_table_arg);
};

class ShocksStatement : public AbstractShocksStatement
{
public:
  ShocksStatement(const det_shocks_type &det_shocks_arg,
                  const var_and_std_shocks_type &var_shocks_arg,
                  const var_and_std_shocks_type &std_shocks_arg,
                  const covar_and_corr_shocks_type &covar_shocks_arg,
                  const covar_and_corr_shocks_type &corr_shocks_arg,
                  const SymbolTable &symbol_table_arg);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

class MShocksStatement : public AbstractShocksStatement
{
public:
  MShocksStatement(const det_shocks_type &det_shocks_arg,
                   const var_and_std_shocks_type &var_shocks_arg,
                   const var_and_std_shocks_type &std_shocks_arg,
                   const covar_and_corr_shocks_type &covar_shocks_arg,
                   const covar_and_corr_shocks_type &corr_shocks_arg,
                   const SymbolTable &symbol_table_arg);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

#endif
