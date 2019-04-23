/*
 * Copyright © 2003-2019 Dynare Team
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

#ifndef _SYMBOLTABLE_HH
#define _SYMBOLTABLE_HH

using namespace std;

#include <map>
#include <string>
#include <utility>
#include <vector>
#include <set>
#include <ostream>

#include "CodeInterpreter.hh"
#include "ExprNode.hh"

using expr_t = class ExprNode *;

//! Types of auxiliary variables
enum class AuxVarType
  {
    endoLead = 0,    //!< Substitute for endo leads >= 2
    endoLag = 1,     //!< Substitute for endo lags >= 2
    exoLead = 2,     //!< Substitute for exo leads >= 1
    exoLag = 3,      //!< Substitute for exo lags >= 1
    expectation = 4, //!< Substitute for Expectation Operator
    diffForward = 5, //!< Substitute for the differentiate of a forward variable
    multiplier = 6,  //!< Multipliers for FOC of Ramsey Problem
    varModel = 7,    //!< Variable for var_model with order > abs(min_lag()) present in model
    diff = 8,        //!< Variable for Diff operator
    diffLag = 9,     //!< Variable for timing between Diff operators
    unaryOp = 10,    //!< Variable for allowing the undiff operator to work when diff was taken of unary op, eg diff(log(x))
    diffLead = 11    //!< Variable for timing between Diff operators
  };

//! Information on some auxiliary variables
class AuxVarInfo
{
private:
  int symb_id; //!< Symbol ID of the auxiliary variable
  AuxVarType type; //!< Its type
  int orig_symb_id; //!< Symbol ID of the endo of the original model represented by this aux var. Only used for avEndoLag and avExoLag.
  int orig_lead_lag; //!< Lead/lag of the endo of the original model represented by this aux var. Only used for avEndoLag and avExoLag.
  int equation_number_for_multiplier; //!< Stores the original constraint equation number associated with this aux var. Only used for avMultiplier.
  int information_set; //! Argument of expectation operator. Only used for avExpectation.
  expr_t expr_node; //! Auxiliary variable definition
  string unary_op; //! Used with AuxUnaryOp
public:
  AuxVarInfo(int symb_id_arg, AuxVarType type_arg, int orig_symb_id, int orig_lead_lag, int equation_number_for_multiplier_arg, int information_set_arg, expr_t expr_node_arg, string unary_op_arg);
  int
  get_symb_id() const
  {
    return symb_id;
  };
  AuxVarType
  get_type() const
  {
    return type;
  };
  int
  get_type_id() const
  {
    return static_cast<int>(type);
  }
  int
  get_orig_symb_id() const
  {
    return orig_symb_id;
  };
  int
  get_orig_lead_lag() const
  {
    return orig_lead_lag;
  };
  int
  get_equation_number_for_multiplier() const
  {
    return equation_number_for_multiplier;
  };
  int
  get_information_set() const
  {
    return information_set;
  };
  expr_t
  get_expr_node() const
  {
    return expr_node;
  };
  string
  get_unary_op() const
  {
    return unary_op;
  };
};

//! Stores the symbol table
/*!
  A symbol is given by its name, and is internally represented by a unique integer.

  When method freeze() is called, computes a distinct sequence of IDs for some types
  (endogenous, exogenous, parameters), which are used by the Matlab/Octave functions.
  We call these "type specific IDs".

  Also manages a TeX name for each symbol, which by default is an empty string.
*/
class SymbolTable
{
private:
  //! Has method freeze() been called?
  bool frozen{false};

  using symbol_table_type = map<string, int>;
  //! Maps strings to symbol IDs
  symbol_table_type symbol_table;

  //! Maps IDs to names
  vector<string> name_table;
  //! Maps IDs to TeX names
  vector<string> tex_name_table;
  //! Maps IDs to string names of variables
  vector<string> long_name_table;
  //! Maps IDs to a pair containing the partition and the partition value
  map<int, map<string, string>> partition_value_map;
  //! Maps IDs to types
  vector<SymbolType> type_table;

  //! Maps symbol IDs to type specific IDs
  vector<int> type_specific_ids;

  //! Maps type specific IDs of endogenous to symbol IDs
  vector<int> endo_ids;
  //! Maps type specific IDs of exogenous to symbol IDs
  vector<int> exo_ids;
  //! Maps type specific IDs of exogenous deterministic to symbol IDs
  vector<int> exo_det_ids;
  //! Maps type specific IDs of parameters to symbol IDs
  vector<int> param_ids;
  //! Information about auxiliary variables
  vector<AuxVarInfo> aux_vars;

  //! Stores the predetermined variables (by symbol IDs)
  set<int> predetermined_variables;

  //! Stores the list of observed variables
  vector<int> varobs;

  //! Stores the list of observed exogenous variables
  vector<int> varexobs;

public:
  SymbolTable();
  //! Thrown when trying to access an unknown symbol (by name)
  class UnknownSymbolNameException
  {
  public:
    //! Symbol name
    string name;
    explicit UnknownSymbolNameException(string name_arg) : name{move(name_arg)}
    {
    }
  };
  //! Thrown when trying to access an unknown symbol (by id)
  class UnknownSymbolIDException
  {
  public:
    //! Symbol ID
    int id;
    explicit UnknownSymbolIDException(int id_arg) : id{id_arg}
    {
    }
  };
  //! Thrown when trying to access an unknown type specific ID
  class UnknownTypeSpecificIDException
  {
  public:
    int tsid;
    SymbolType type;
    UnknownTypeSpecificIDException(int tsid_arg, SymbolType type_arg) : tsid{tsid_arg}, type{type_arg}
    {
    }
  };
  //! Thrown when trying to declare a symbol twice
  class AlreadyDeclaredException
  {
  public:
    //! Symbol name
    string name;
    //! Was the previous declaration done with the same symbol type ?
    bool same_type;
    AlreadyDeclaredException(string name_arg, bool same_type_arg) : name{move(name_arg)}, same_type{same_type_arg}
    {
    }
  };
  //! Thrown when table is frozen and trying to modify it
  class FrozenException
  {
  };
  //! Thrown when trying to use the result of freeze() while this method has not yet been called
  class NotYetFrozenException
  {
  };
  //! Thrown when searchAuxiliaryVars() failed
  class SearchFailedException
  {
  public:
    int orig_symb_id, orig_lead_lag, symb_id;
    SearchFailedException(int orig_symb_id_arg, int orig_lead_lag_arg) : orig_symb_id{orig_symb_id_arg},
                                                                         orig_lead_lag{orig_lead_lag_arg}
    {
    }
    explicit SearchFailedException(int symb_id_arg) : symb_id{symb_id_arg}
    {
    }
  };

private:
  //! Factorized code for adding aux lag variables
  int addLagAuxiliaryVarInternal(bool endo, int orig_symb_id, int orig_lead_lag, expr_t arg) noexcept(false);
  //! Factorized code for adding aux lead variables
  int addLeadAuxiliaryVarInternal(bool endo, int index, expr_t arg) noexcept(false);
  //! Factorized code for Json writing
  void writeJsonVarVector(ostream &output, const vector<int> &varvec) const;
  //! Factorized code for asserting that 0 <= symb_id <= symbol_table.size()
  inline void validateSymbID(int symb_id) const noexcept(false);
public:
  //! Add a symbol
  /*! Returns the symbol ID */
  int addSymbol(const string &name, SymbolType type, const string &tex_name, const vector<pair<string, string>> &partition_value) noexcept(false);
  //! Add a symbol without its TeX name (will be equal to its name)
  /*! Returns the symbol ID */
  int addSymbol(const string &name, SymbolType type) noexcept(false);
  //! Adds an auxiliary variable for endogenous with lead >= 2
  /*!
    \param[in] index Used to construct the variable name
    \return the symbol ID of the new symbol */
  int addEndoLeadAuxiliaryVar(int index, expr_t arg) noexcept(false);
  //! Adds an auxiliary variable for endogenous with lag >= 2
  /*!
    \param[in] orig_symb_id symbol ID of the endogenous declared by the user that this new variable will represent
    \param[in] orig_lead_lag lag value such that this new variable will be equivalent to orig_symb_id(orig_lead_lag)
    \return the symbol ID of the new symbol */
  int addEndoLagAuxiliaryVar(int orig_symb_id, int orig_lead_lag, expr_t arg) noexcept(false);
  //! Adds an auxiliary variable for endogenous with lead >= 1
  /*!
    \param[in] index Used to construct the variable name
    \return the symbol ID of the new symbol */
  int addExoLeadAuxiliaryVar(int index, expr_t arg) noexcept(false);
  //! Adds an auxiliary variable for exogenous with lag >= 1
  /*!
    \param[in] orig_symb_id symbol ID of the exogenous declared by the user that this new variable will represent
    \param[in] orig_lead_lag lag value such that this new variable will be equivalent to orig_symb_id(orig_lead_lag)
    \return the symbol ID of the new symbol */
  int addExoLagAuxiliaryVar(int orig_symb_id, int orig_lead_lag, expr_t arg) noexcept(false);
  //! Adds an auxiliary variable for the expectation operator
  /*!
    \param[in] information_set information set (possibly negative) of the expectation operator
    \param[in] index Used to construct the variable name
    \return the symbol ID of the new symbol
  */
  int addExpectationAuxiliaryVar(int information_set, int index, expr_t arg) noexcept(false);
  //! Adds an auxiliary variable for the multiplier for the FOCs of the Ramsey Problem
  /*!
    \param[in] index Used to construct the variable name
    \return the symbol ID of the new symbol
  */
  int addMultiplierAuxiliaryVar(int index) noexcept(false);
  //! Adds an auxiliary variable for the (time) differentiate of a forward var
  /*!
    \param[in] orig_symb_id The symb_id of the forward variable
    \return the symbol ID of the new symbol
  */
  int addDiffForwardAuxiliaryVar(int orig_symb_id, expr_t arg) noexcept(false);
  //! Searches auxiliary variables which are substitutes for a given symbol_id and lead/lag
  /*!
    The search is only performed among auxiliary variables of endo/exo lag.
    \return the symbol ID of the auxiliary variable
    Throws an exception if match not found.
  */
  int searchAuxiliaryVars(int orig_symb_id, int orig_lead_lag) const noexcept(false);
  //! Serches aux_vars for the aux var represented by aux_var_symb_id and returns its associated orig_symb_id
  int getOrigSymbIdForAuxVar(int aux_var_symb_id) const noexcept(false);
  //! Searches for diff aux var and finds the original lag associated with this variable
  int getOrigLeadLagForDiffAuxVar(int diff_aux_var_symb_id) const noexcept(false);
  //! Searches for diff aux var and finds the symb id associated with this variable
  int getOrigSymbIdForDiffAuxVar(int diff_aux_var_symb_id) const noexcept(false);
  //! Adds an auxiliary variable when var_model is used with an order that is greater in absolute value
  //! than the largest lag present in the model.
  int addVarModelEndoLagAuxiliaryVar(int orig_symb_id, int orig_lead_lag, expr_t expr_arg) noexcept(false);
  //! Adds an auxiliary variable when the diff operator is encountered
  int addDiffAuxiliaryVar(int index, expr_t expr_arg) noexcept(false);
  int addDiffAuxiliaryVar(int index, expr_t expr_arg, int orig_symb_id, int orig_lag) noexcept(false);
  //! Takes care of timing between diff statements
  int addDiffLagAuxiliaryVar(int index, expr_t expr_arg, int orig_symb_id, int orig_lag) noexcept(false);
  //! Takes care of timing between diff statements
  int addDiffLeadAuxiliaryVar(int index, expr_t expr_arg, int orig_symb_id, int orig_lead) noexcept(false);
  //! An Auxiliary variable for a unary op
  int addUnaryOpAuxiliaryVar(int index, expr_t expr_arg, string unary_op, int orig_symb_id = -1, int orig_lag = 0) noexcept(false);
  //! Returns the number of auxiliary variables
  int
  AuxVarsSize() const
  {
    return aux_vars.size();
  };
  //! Retruns expr_node for an auxiliary variable
  expr_t getAuxiliaryVarsExprNode(int symb_id) const noexcept(false);
  //! Tests if symbol already exists
  inline bool exists(const string &name) const;
  //! Get symbol name (by ID)
  inline string getName(int id) const noexcept(false);
  //! Get TeX name
  inline string getTeXName(int id) const noexcept(false);
  //! Get long name
  inline string getLongName(int id) const noexcept(false);
  //! Returns true if the partition name is the first encountered for the type of variable represented by id
  bool isFirstOfPartitionForType(int id) const noexcept(false);
  //! Returns a list of partitions and symbols that belong to that partition
  map<string, map<int, string>> getPartitionsForType(SymbolType st) const noexcept(false);
  //! Get type (by ID)
  inline SymbolType getType(int id) const noexcept(false);
  //! Get type (by name)
  inline SymbolType getType(const string &name) const noexcept(false);
  //! Get ID (by name)
  inline int getID(const string &name) const noexcept(false);
  //! Get ID (by type specific ID)
  int getID(SymbolType type, int tsid) const noexcept(false);
  //! Freeze symbol table
  void freeze() noexcept(false);
  //! unreeze symbol table
  //! Used after having written JSON files
  void unfreeze();
  //! Change the type of a symbol
  void changeType(int id, SymbolType newtype) noexcept(false);
  //! Get type specific ID (by symbol ID)
  inline int getTypeSpecificID(int id) const noexcept(false);
  //! Get type specific ID (by symbol name)
  inline int getTypeSpecificID(const string &name) const noexcept(false);
  //! Get number of endogenous variables
  inline int endo_nbr() const noexcept(false);
  //! Get number of exogenous variables
  inline int exo_nbr() const noexcept(false);
  //! Get number of exogenous deterministic variables
  inline int exo_det_nbr() const noexcept(false);
  //! Get number of parameters
  inline int param_nbr() const noexcept(false);
  //! Returns the greatest symbol ID (the smallest is zero)
  inline int maxID();
  //! Get number of user-declared endogenous variables (without the auxiliary variables)
  inline int orig_endo_nbr() const noexcept(false);
  //! Write output of this class
  void writeOutput(ostream &output) const noexcept(false);
  //! Write JSON Output
  void writeJsonOutput(ostream &output) const;
  //! Write Julia output of this class
  void writeJuliaOutput(ostream &output) const noexcept(false);
  //! Write C output of this class
  void writeCOutput(ostream &output) const noexcept(false);
  //! Write CC output of this class
  void writeCCOutput(ostream &output) const noexcept(false);
  //! Mark a symbol as predetermined variable
  void markPredetermined(int symb_id) noexcept(false);
  //! Test if a given symbol is a predetermined variable
  bool isPredetermined(int symb_id) const noexcept(false);
  //! Return the number of predetermined variables
  int predeterminedNbr() const;
  //! Add an observed variable
  void addObservedVariable(int symb_id) noexcept(false);
  //! Return the number of observed variables
  int observedVariablesNbr() const;
  //! Is a given symbol in the set of observed variables
  bool isObservedVariable(int symb_id) const;
  //! Return the index of a given observed variable in the vector of all observed variables
  int getObservedVariableIndex(int symb_id) const;
  //! Add an observed exogenous variable
  void addObservedExogenousVariable(int symb_id) noexcept(false);
  //! Return the number of observed exogenous variables
  int observedExogenousVariablesNbr() const;
  //! Is a given symbol in the set of observed exogenous variables
  bool isObservedExogenousVariable(int symb_id) const;
  //! Return the index of a given observed exogenous variable in the vector of all observed variables
  int getObservedExogenousVariableIndex(int symb_id) const;
  vector <int> getTrendVarIds() const;
  //! Get list of exogenous variables
  set <int> getExogenous() const;
  //! Get list of exogenous variables
  set <int> getObservedExogenous() const;
  //! Get list of endogenous variables
  set <int> getEndogenous() const;
  //! Is a given symbol an auxiliary variable
  bool isAuxiliaryVariable(int symb_id) const;
  //! Is a given symbol an auxiliary variable but not a Lagrange multiplier
  bool isAuxiliaryVariableButNotMultiplier(int symb_id) const;
  //! Is a given symbol a diff, diff lead, or diff lag auxiliary variable
  bool isDiffAuxiliaryVariable(int symb_id) const;
  //! Get list of endogenous variables without aux vars
  set <int> getOrigEndogenous() const;
  //! Returns the original symbol corresponding to this variable
  /* If symb_id is not an auxiliary var, returns symb_id. Otherwise,
     repeatedly call getOrigSymbIDForAuxVar() until an original
     (non-auxiliary) variable is found. */
  int getUltimateOrigSymbID(int symb_id) const;
};

inline void
SymbolTable::validateSymbID(int symb_id) const noexcept(false)
{
  if (symb_id < 0 || symb_id > static_cast<int>(symbol_table.size()))
    throw UnknownSymbolIDException(symb_id);
}

inline bool
SymbolTable::exists(const string &name) const
{
  auto iter = symbol_table.find(name);
  return (iter != symbol_table.end());
}

inline string
SymbolTable::getName(int id) const noexcept(false)
{
  validateSymbID(id);
  return name_table[id];
}

inline string
SymbolTable::getTeXName(int id) const noexcept(false)
{
  validateSymbID(id);
  return tex_name_table[id];
}

inline string
SymbolTable::getLongName(int id) const noexcept(false)
{
  validateSymbID(id);
  return long_name_table[id];
}

inline SymbolType
SymbolTable::getType(int id) const noexcept(false)
{
  validateSymbID(id);
  return type_table[id];
}

inline SymbolType
SymbolTable::getType(const string &name) const noexcept(false)
{
  return getType(getID(name));
}

inline int
SymbolTable::getID(const string &name) const noexcept(false)
{
  auto iter = symbol_table.find(name);
  if (iter != symbol_table.end())
    return iter->second;
  else
    throw UnknownSymbolNameException(name);
}

inline int
SymbolTable::getTypeSpecificID(int id) const noexcept(false)
{
  if (!frozen)
    throw NotYetFrozenException();

  validateSymbID(id);

  return type_specific_ids[id];
}

inline int
SymbolTable::getTypeSpecificID(const string &name) const noexcept(false)
{
  return getTypeSpecificID(getID(name));
}

inline int
SymbolTable::endo_nbr() const noexcept(false)
{
  if (!frozen)
    throw NotYetFrozenException();

  return endo_ids.size();
}

inline int
SymbolTable::exo_nbr() const noexcept(false)
{
  if (!frozen)
    throw NotYetFrozenException();

  return exo_ids.size();
}

inline int
SymbolTable::exo_det_nbr() const noexcept(false)
{
  if (!frozen)
    throw NotYetFrozenException();

  return exo_det_ids.size();
}

inline int
SymbolTable::param_nbr() const noexcept(false)
{
  if (!frozen)
    throw NotYetFrozenException();

  return param_ids.size();
}

inline int
SymbolTable::maxID()
{
  return symbol_table.size() - 1;
}

inline int
SymbolTable::orig_endo_nbr() const noexcept(false)
{
  return (endo_nbr() - aux_vars.size());
}

#endif
