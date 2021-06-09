/*
 * Copyright © 2003-2021 Dynare Team
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

#ifndef _NUMERICALINITIALIZATION_HH
#define _NUMERICALINITIALIZATION_HH

using namespace std;

#include <string>
#include <vector>
#include <map>

#include "SymbolTable.hh"
#include "ExprNode.hh"
#include "Statement.hh"

class InitParamStatement : public Statement
{
private:
  const int symb_id;
  const expr_t param_value;
  const SymbolTable &symbol_table;
public:
  InitParamStatement(int symb_id_arg, const expr_t param_value_arg,
                     const SymbolTable &symbol_table_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
  //! Fill eval context with parameter value
  void fillEvalContext(eval_context_t &eval_context) const;
};

class InitOrEndValStatement : public Statement
{
public:
  /*!
    We use a vector instead of a map, since the order of declaration matters:
    an initialization can depend on a previously initialized variable inside the block
  */
  using init_values_t = vector<pair<int, expr_t>>;
protected:
  const init_values_t init_values;
  const SymbolTable &symbol_table;
  const bool all_values_required;
public:
  InitOrEndValStatement(init_values_t init_values_arg,
                        const SymbolTable &symbol_table_arg,
                        bool all_values_required_arg);
  //! Return set of unused variables by type
  set<int> getUninitializedVariables(SymbolType type);
  //! Fill eval context with variables values
  void fillEvalContext(eval_context_t &eval_context) const;
protected:
  void writeInitValues(ostream &output) const;
  void writeJsonInitValues(ostream &output) const;
};

class InitValStatement : public InitOrEndValStatement
{
public:
  InitValStatement(const init_values_t &init_values_arg,
                   const SymbolTable &symbol_table_arg,
                   bool all_values_required_arg);
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
  //! Writes initializations for oo_.exo_simul and oo_.exo_det_simul
  void writeOutputPostInit(ostream &output) const;
};

class EndValStatement : public InitOrEndValStatement
{
public:
  EndValStatement(const init_values_t &init_values_arg,
                  const SymbolTable &symbol_table_arg,
                  bool all_values_required_arg);
  //! Workaround for trac ticket #35
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class HistValStatement : public Statement
{
public:
  /*!
    Contrary to Initval and Endval, we use a map, since it is impossible to reuse
    a given initialization value in a second initialization inside the block.
    Maps pairs (symbol_id, lag) to expr_t
  */
  using hist_values_t = map<pair<int, int>, expr_t>;
private:
  const hist_values_t hist_values;
  const SymbolTable &symbol_table;
  const bool all_values_required;
public:
  HistValStatement(hist_values_t hist_values_arg,
                   const SymbolTable &symbol_table_arg,
                   bool all_values_required_arg);
  //! Workaround for trac ticket #157
  void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings) override;
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class InitvalFileStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  explicit InitvalFileStatement(OptionsList options_list_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class HistvalFileStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  explicit HistvalFileStatement(OptionsList options_list_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class HomotopyStatement : public Statement
{
public:
  //! Stores the declarations of homotopy_setup
  /*! Order matter so we use a vector. First expr_t can be NULL if no initial value given. */
  using homotopy_values_t = vector<tuple<int, expr_t, expr_t>>;
private:
  const homotopy_values_t homotopy_values;
  const SymbolTable &symbol_table;
public:
  HomotopyStatement(homotopy_values_t homotopy_values_arg,
                    const SymbolTable &symbol_table_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class SaveParamsAndSteadyStateStatement : public Statement
{
private:
  const string filename;
public:
  explicit SaveParamsAndSteadyStateStatement(string filename_arg);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  void writeJsonOutput(ostream &output) const override;
};

class LoadParamsAndSteadyStateStatement : public Statement
{
private:
  const SymbolTable &symbol_table;
  //! Content of the file
  /*! Maps symbol ID to numeric value (stored as string) */
  map<int, string> content;
public:
  LoadParamsAndSteadyStateStatement(const string &filename,
                                    const SymbolTable &symbol_table_arg,
                                    WarningConsolidation &warnings);
  void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const override;
  //! Fill eval context with parameters/variables values
  void fillEvalContext(eval_context_t &eval_context) const;
  void writeJsonOutput(ostream &output) const override;
};

#endif
