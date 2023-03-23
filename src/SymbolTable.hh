/*
 * Copyright © 2003-2023 Dynare Team
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

#ifndef _SYMBOLTABLE_HH
#define _SYMBOLTABLE_HH

#include <map>
#include <string>
#include <utility>
#include <vector>
#include <set>
#include <ostream>
#include <optional>

#include "CommonEnums.hh"
#include "ExprNode.hh"

using namespace std;

using expr_t = class ExprNode *;

//! Types of auxiliary variables
enum class AuxVarType
  {
   endoLead = 0, //!< Substitute for endo leads >= 2
   endoLag = 1, //!< Substitute for endo lags >= 2
   exoLead = 2, //!< Substitute for exo leads >= 1
   exoLag = 3, //!< Substitute for exo lags >= 1
   expectation = 4, //!< Substitute for Expectation Operator
   diffForward = 5, /* Substitute for the differentiate of a forward variable,
                       for the differentiate_forward_vars option.
                       N.B.: nothing to do with the diff() operator! */
   multiplier = 6, //!< Multipliers for FOC of Ramsey Problem
   logTransform = 7, //!< Log-transformation of a variable declared with “var(log)”
   diff = 8, //!< Variable for Diff operator
   diffLag = 9, //!< Variable for timing between Diff operators (lag)
   unaryOp = 10, //!< Variable for allowing the undiff operator to work when diff was taken of unary op, eg diff(log(x))
   diffLead = 11, //!< Variable for timing between Diff operators (lead)
   pacExpectation = 12, //!< Variable created for the substitution of the pac_expectation operator
   pacTargetNonstationary = 13 //!< Variable created for the substitution of the pac_target_nonstationary operator
  };

//! Information on some auxiliary variables
struct AuxVarInfo
{
  const int symb_id; // Symbol ID of the auxiliary variable
  const AuxVarType type; // Its type
  const optional<int> orig_symb_id; /* Symbol ID of the (only) endo that appears on the RHS of
                                       the definition of this auxvar.
                                       Used by endoLag, exoLag, diffForward, logTransform, diff,
                                       diffLag, diffLead and unaryOp.
                                       For diff and unaryOp, if the argument expression is more complex
                                       than than a simple variable, this value is unset
                                       (hence the need for std::optional). */
  const optional<int> orig_lead_lag; /* Lead/lag of the (only) endo as it appears on the RHS of the
                                        definition of this auxvar. Only set if orig_symb_id is set
                                        (in particular, for diff and unaryOp, unset
                                        if orig_symb_id is unset).
                                        For diff and diffForward, since the definition of the
                                        auxvar is a time difference, the value corresponds to the
                                        time index of the first term of that difference. */
  const int equation_number_for_multiplier; /* Stores the original constraint equation number
                                               associated with this aux var. Only used for
                                               avMultiplier. */
  const int information_set; // Argument of expectation operator. Only used for avExpectation.
  const expr_t expr_node; // Auxiliary variable definition
  const string unary_op; // Used with AuxUnaryOp

  int
  get_type_id() const
  {
    return static_cast<int>(type);
  }
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
  map<int, int> type_specific_ids;

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

  //! Stores the endogenous variables declared with “var(log)”
  set<int> with_log_transform;

public:
  //! Thrown when trying to access an unknown symbol (by name)
  struct UnknownSymbolNameException
  {
    //! Symbol name
    const string name;
  };
  //! Thrown when trying to access an unknown symbol (by id)
  struct UnknownSymbolIDException
  {
    //! Symbol ID
    const int id;
  };
  //! Thrown when trying to access an unknown type specific ID
  struct UnknownTypeSpecificIDException
  {
    const int tsid;
    const SymbolType type;
  };
  /* Thrown when requesting the type specific ID of a symbol which doesn’t
     have one */
  struct NoTypeSpecificIDException
  {
    const int symb_id;
  };
  //! Thrown when trying to declare a symbol twice
  struct AlreadyDeclaredException
  {
    //! Symbol name
    const string name;
    //! Was the previous declaration done with the same symbol type ?
    const bool same_type;
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
  /* Adds an auxiliary variable associated to an endogenous declared with
     “var(log)”.
     – orig_symb_id is the symbol ID of the original variable
     – orig_lead_lag is typically 0
     – expr_arg is typically log(orig_symb_id)
  */
  int addLogTransformAuxiliaryVar(int orig_symb_id, int orig_lead_lag, expr_t expr_arg) noexcept(false);
  //! Adds an auxiliary variable for the (time) differentiate of a forward var
  /*!
    \param[in] orig_symb_id The symb_id of the forward variable
    \return the symbol ID of the new symbol
  */
  int addDiffForwardAuxiliaryVar(int orig_symb_id, int orig_lead_lag, expr_t arg) noexcept(false);
  //! Searches auxiliary variables which are substitutes for a given symbol_id and lead/lag
  /*!
    The search is only performed among auxiliary variables of endo/exo lag.
    \return the symbol ID of the auxiliary variable
    Throws an exception if match not found.
  */
  int searchAuxiliaryVars(int orig_symb_id, int orig_lead_lag) const noexcept(false);
  /* Searches aux_vars for the aux var represented by aux_var_symb_id and
     returns its associated orig_symb_id.
     Throws an UnknownSymbolIDException if there is no orig_symb_id associated to
     this auxvar (either because it’s of the wrong type, or because there is
     no such orig var for this specific auxvar, in case of complex expressions
     in diff or unaryOp). */
  int getOrigSymbIdForAuxVar(int aux_var_symb_id_arg) const noexcept(false);
  /* Unrolls a chain of diffLag or diffLead aux vars until it founds a (regular) diff aux
     var. In other words:
     - if the arg is a (regu) diff aux var, returns the arg
     - if the arg is a diffLag/diffLead, get its orig symb ID, and call the
       method recursively
     - if the arg is something else, throw an UnknownSymbolIDException
       exception
     The 2nd input/output arguments are used to track leads/lags. The 2nd
     output argument is equal to the 2nd input argument, shifted by as many
     lead/lags were encountered in the chain (a diffLag decreases it, a
     diffLead increases it). */
  pair<int, int> unrollDiffLeadLagChain(int symb_id, int lag) const noexcept(false);
  //! Adds an auxiliary variable when the diff operator is encountered
  int addDiffAuxiliaryVar(int index, expr_t expr_arg, optional<int> orig_symb_id = nullopt, optional<int> orig_lag = nullopt) noexcept(false);
  //! Takes care of timing between diff statements
  int addDiffLagAuxiliaryVar(int index, expr_t expr_arg, int orig_symb_id, int orig_lag) noexcept(false);
  //! Takes care of timing between diff statements
  int addDiffLeadAuxiliaryVar(int index, expr_t expr_arg, int orig_symb_id, int orig_lead) noexcept(false);
  //! An Auxiliary variable for a unary op
  int addUnaryOpAuxiliaryVar(int index, expr_t expr_arg, string unary_op, optional<int> orig_symb_id = nullopt, optional<int> orig_lag = nullopt) noexcept(false);
  //! An auxiliary variable for a pac_expectation operator
  int addPacExpectationAuxiliaryVar(const string &name, expr_t expr_arg);
  //! An auxiliary variable for a pac_target_nonstationary operator
  int addPacTargetNonstationaryAuxiliaryVar(const string &name, expr_t expr_arg);
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
  //! Mark a symbol as predetermined variable
  void markPredetermined(int symb_id) noexcept(false);
  //! Mark an endogenous as having been declared with “var(log)”
  void markWithLogTransform(int symb_id) noexcept(false);
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
  //! Is a given symbol a diff, diff lead, or diff lag auxiliary variable
  bool isDiffAuxiliaryVariable(int symb_id) const;
  //! Get list of endogenous variables without aux vars
  set <int> getOrigEndogenous() const;
  //! Returns the original symbol corresponding to this variable
  /* If symb_id has no original variable, returns symb_id. Otherwise,
     repeatedly call getOrigSymbIDForAuxVar() until an original variable is
     found. Note that the result may be an auxiliary variable if the latter has
     no original variable (e.g. aux var for lead, Lagrange Multiplier or diff
     associated to a complex expression). */
  int getUltimateOrigSymbID(int symb_id) const;
  //! If this is a Lagrange multiplier, return its associated equation number; otherwise return nullopt
  optional<int> getEquationNumberForMultiplier(int symb_id) const;
  /* Return all the information about a given auxiliary variable. Throws
     UnknownSymbolIDException if it is not an aux var */
  const AuxVarInfo &getAuxVarInfo(int symb_id) const;
  // Returns the set of all endogenous declared with “var(log)”
  const set<int> &getVariablesWithLogTransform() const;
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
  return symbol_table.contains(name);
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
  if (auto iter = symbol_table.find(name);
      iter != symbol_table.end())
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

  if (auto it = type_specific_ids.find(id);
      it != type_specific_ids.end())
    return it->second;
  else
    throw NoTypeSpecificIDException(id);
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
  return endo_nbr() - aux_vars.size();
}

inline const AuxVarInfo &
SymbolTable::getAuxVarInfo(int symb_id) const
{
  for (const auto &aux_var : aux_vars)
    if (aux_var.symb_id == symb_id)
      return aux_var;
  throw UnknownSymbolIDException(symb_id);
}

#endif
