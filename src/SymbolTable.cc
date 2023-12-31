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

#include <algorithm>
#include <sstream>
#include <iostream>
#include <cassert>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#include <boost/algorithm/string/replace.hpp>
#pragma GCC diagnostic pop
#include <utility>

#include "SymbolTable.hh"

int
SymbolTable::addSymbol(const string &name, SymbolType type, const string &tex_name, const vector<pair<string, string>> &partition_value) noexcept(false)
{
  if (frozen)
    throw FrozenException();

  if (exists(name))
    {
      if (type_table[getID(name)] == type)
        throw AlreadyDeclaredException{name, true};
      else
        throw AlreadyDeclaredException{name, false};
    }

  string final_tex_name = tex_name;
  if (final_tex_name.empty())
    {
      final_tex_name = name;
      size_t pos = 0;
      while ((pos = final_tex_name.find('_', pos)) != string::npos)
        {
          final_tex_name.insert(pos, R"(\)");
          pos += 2;
        }
    }

  string final_long_name = name;
  bool non_long_name_partition_exists = false;
  for (const auto &it : partition_value)
    if (it.first == "long_name")
      final_long_name = it.second;
    else
      non_long_name_partition_exists = true;

  int id = symbol_table.size();

  symbol_table[name] = id;
  type_table.push_back(type);
  name_table.push_back(name);
  tex_name_table.push_back(final_tex_name);
  long_name_table.push_back(final_long_name);
  if (non_long_name_partition_exists)
    {
      map<string, string> pmv;
      for (const auto &it : partition_value)
        pmv[it.first] = it.second;
      partition_value_map[id] = pmv;
    }
  return id;
}

int
SymbolTable::addSymbol(const string &name, SymbolType type) noexcept(false)
{
  return addSymbol(name, type, "", {});
}

void
SymbolTable::freeze() noexcept(false)
{
  if (frozen)
    throw FrozenException();

  frozen = true;

  for (int i = 0; i < static_cast<int>(symbol_table.size()); i++)
    {
      int tsi;
      switch (getType(i))
        {
        case SymbolType::endogenous:
          tsi = endo_ids.size();
          endo_ids.push_back(i);
          break;
        case SymbolType::exogenous:
          tsi = exo_ids.size();
          exo_ids.push_back(i);
          break;
        case SymbolType::exogenousDet:
          tsi = exo_det_ids.size();
          exo_det_ids.push_back(i);
          break;
        case SymbolType::parameter:
          tsi = param_ids.size();
          param_ids.push_back(i);
          break;
        default:
          continue;
        }
      type_specific_ids[i] = tsi;
    }
}

void
SymbolTable::unfreeze()
{
  frozen = false;
  endo_ids.clear();
  exo_ids.clear();
  exo_det_ids.clear();
  param_ids.clear();
  type_specific_ids.clear();
}

void
SymbolTable::changeType(int id, SymbolType newtype) noexcept(false)
{
  if (frozen)
    throw FrozenException();

  validateSymbID(id);

  type_table[id] = newtype;
}

int
SymbolTable::getID(SymbolType type, int tsid) const noexcept(false)
{
  if (!frozen)
    throw NotYetFrozenException();

  switch (type)
    {
    case SymbolType::endogenous:
      if (tsid < 0 || tsid >= static_cast<int>(endo_ids.size()))
        throw UnknownTypeSpecificIDException{tsid, type};
      else
        return endo_ids[tsid];
    case SymbolType::exogenous:
      if (tsid < 0 || tsid >= static_cast<int>(exo_ids.size()))
        throw UnknownTypeSpecificIDException{tsid, type};
      else
        return exo_ids[tsid];
    case SymbolType::exogenousDet:
      if (tsid < 0 || tsid >= static_cast<int>(exo_det_ids.size()))
        throw UnknownTypeSpecificIDException{tsid, type};
      else
        return exo_det_ids[tsid];
    case SymbolType::parameter:
      if (tsid < 0 || tsid >= static_cast<int>(param_ids.size()))
        throw UnknownTypeSpecificIDException{tsid, type};
      else
        return param_ids[tsid];
    default:
      throw UnknownTypeSpecificIDException{tsid, type};
    }
}

map<string, map<int, string>>
SymbolTable::getPartitionsForType(SymbolType st) const noexcept(false)
{
  map<string, map<int, string>> partitions;
  for (const auto &it : partition_value_map)
    if (getType(it.first) == st)
      for (const auto &it1 : it.second)
        partitions[it1.first][it.first] = it1.second;
  return partitions;
}

void
SymbolTable::writeOutput(ostream &output) const noexcept(false)
{
  if (!frozen)
    throw NotYetFrozenException();

  if (exo_nbr() > 0)
    {
      output << "M_.exo_names = cell(" << exo_nbr() << ",1);" << endl;
      output << "M_.exo_names_tex = cell(" << exo_nbr() << ",1);" << endl;
      output << "M_.exo_names_long = cell(" << exo_nbr() << ",1);" << endl;
      for (int id = 0; id < exo_nbr(); id++)
        output << "M_.exo_names(" << id+1 << ") = {'" << getName(exo_ids[id]) << "'};" << endl
               << "M_.exo_names_tex(" << id+1 << ") = {'" << getTeXName(exo_ids[id]) << "'};" << endl
               << "M_.exo_names_long(" << id+1 << ") = {'" << getLongName(exo_ids[id]) << "'};" << endl;
      for (auto &partition : getPartitionsForType(SymbolType::exogenous))
        if (partition.first != "long_name")
          {
            output << "M_.exo_partitions." << partition.first << " = { ";
            for (int id = 0; id < exo_nbr(); id++)
              {
                output << "'";
                if (auto it1 = partition.second.find(exo_ids[id]);
                    it1 != partition.second.end())
                  output << it1->second;
                output << "' ";
              }
            output << "};" << endl;
            if (partition.first == "status")
              output << "M_ = set_observed_exogenous_variables(M_);" << endl;
            if (partition.first == "used")
              output << "M_ = set_exogenous_variables_for_simulation(M_);" << endl;
          }
    }
  else
    {
      output << "M_.exo_names = {};" << endl;
      output << "M_.exo_names_tex = {};" << endl;
      output << "M_.exo_names_long = {};" << endl;
    }

  if (exo_det_nbr() > 0)
    {
      output << "M_.exo_det_names = cell(" << exo_det_nbr() << ",1);" << endl;
      output << "M_.exo_det_names_tex = cell(" << exo_det_nbr() << ",1);" << endl;
      output << "M_.exo_det_names_long = cell(" << exo_det_nbr() << ",1);" << endl;
      for (int id = 0; id < exo_det_nbr(); id++)
        output << "M_.exo_det_names(" << id+1 << ") = {'" << getName(exo_det_ids[id]) << "'};" << endl
               << "M_.exo_det_names_tex(" << id+1 << ") = {'" << getTeXName(exo_det_ids[id]) << "'};" << endl
               << "M_.exo_det_names_long(" << id+1 << ") = {'" << getLongName(exo_det_ids[id]) << "'};" << endl;
      output << "M_.exo_det_partitions = struct();" << endl;
      for (auto &partition : getPartitionsForType(SymbolType::exogenousDet))
        if (partition.first != "long_name")
          {
            output << "M_.exo_det_partitions." << partition.first << " = { ";
            for (int id = 0; id < exo_det_nbr(); id++)
              {
                output << "'";
                if (auto it1 = partition.second.find(exo_det_ids[id]);
                    it1 != partition.second.end())
                  output << it1->second;
                output << "' ";
              }
            output << "};" << endl;
          }
    }

  if (endo_nbr() > 0)
    {
      output << "M_.endo_names = cell(" << endo_nbr() << ",1);" << endl;
      output << "M_.endo_names_tex = cell(" << endo_nbr() << ",1);" << endl;
      output << "M_.endo_names_long = cell(" << endo_nbr() << ",1);" << endl;
      for (int id = 0; id < endo_nbr(); id++)
        output << "M_.endo_names(" << id+1 << ") = {'" << getName(endo_ids[id]) << "'};" << endl
               << "M_.endo_names_tex(" << id+1 << ") = {'" << getTeXName(endo_ids[id]) << "'};" << endl
               << "M_.endo_names_long(" << id+1 << ") = {'" << getLongName(endo_ids[id]) << "'};" << endl;
      output << "M_.endo_partitions = struct();" << endl;
      for (auto &partition : getPartitionsForType(SymbolType::endogenous))
        if (partition.first != "long_name")
          {
            output << "M_.endo_partitions." << partition.first << " = { ";
            for (int id = 0; id < endo_nbr(); id++)
              {
                output << "'";
                if (auto it1 = partition.second.find(endo_ids[id]);
                    it1 != partition.second.end())
                  output << it1->second;
                output << "' ";
              }
            output << "};" << endl;
          }
    }

  if (param_nbr() > 0)
    {
      output << "M_.param_names = cell(" << param_nbr() << ",1);" << endl;
      output << "M_.param_names_tex = cell(" << param_nbr() << ",1);" << endl;
      output << "M_.param_names_long = cell(" << param_nbr() << ",1);" << endl;
      for (int id = 0; id < param_nbr(); id++)
        {
          output << "M_.param_names(" << id+1 << ") = {'" << getName(param_ids[id]) << "'};" << endl
                 << "M_.param_names_tex(" << id+1 << ") = {'" << getTeXName(param_ids[id]) << "'};" << endl
                 << "M_.param_names_long(" << id+1 << ") = {'" << getLongName(param_ids[id]) << "'};" << endl;
          if (getName(param_ids[id]) == "dsge_prior_weight")
            output << "options_.dsge_var = 1;" << endl;
        }
      output << "M_.param_partitions = struct();" << endl;
      for (auto &partition : getPartitionsForType(SymbolType::parameter))
        if (partition.first != "long_name")
          {
            output << "M_.param_partitions." << partition.first << " = { ";
            for (int id = 0; id < param_nbr(); id++)
              {
                output << "'";
                if (auto it1 = partition.second.find(param_ids[id]);
                    it1 != partition.second.end())
                  output << it1->second;
                output << "' ";
              }
            output << "};" << endl;
          }
    }
  else
    {
      output << "M_.param_names = {};" << endl;
      output << "M_.param_names_tex = {};" << endl;
      output << "M_.param_names_long = {};" << endl;
    }

  output << "M_.exo_det_nbr = " << exo_det_nbr() << ";" << endl
         << "M_.exo_nbr = " << exo_nbr() << ";" << endl
         << "M_.endo_nbr = " << endo_nbr() << ";" << endl
         << "M_.param_nbr = " << param_nbr() << ";" << endl;

  // Write the auxiliary variable table
  output << "M_.orig_endo_nbr = " << orig_endo_nbr() << ";" << endl;
  if (aux_vars.size() == 0)
    output << "M_.aux_vars = [];" << endl;
  else
    for (int i = 0; i < static_cast<int>(aux_vars.size()); i++)
      {
        output << "M_.aux_vars(" << i+1 << ").endo_index = " << getTypeSpecificID(aux_vars[i].symb_id)+1 << ";" << endl
               << "M_.aux_vars(" << i+1 << ").type = " << aux_vars[i].get_type_id() << ";" << endl;
        switch (aux_vars[i].type)
          {
          case AuxVarType::endoLead:
          case AuxVarType::exoLead:
          case AuxVarType::expectation:
          case AuxVarType::pacExpectation:
          case AuxVarType::pacTargetNonstationary:
            break;
          case AuxVarType::endoLag:
          case AuxVarType::exoLag:
          case AuxVarType::logTransform:
          case AuxVarType::diffLag:
          case AuxVarType::diffLead:
          case AuxVarType::diffForward:
            output << "M_.aux_vars(" << i+1 << ").orig_index = " << getTypeSpecificID(aux_vars[i].orig_symb_id.value())+1 << ";" << endl
                   << "M_.aux_vars(" << i+1 << ").orig_lead_lag = " << aux_vars[i].orig_lead_lag.value() << ";" << endl;
            break;
          case AuxVarType::unaryOp:
            output << "M_.aux_vars(" << i+1 << ").unary_op = '" << aux_vars[i].unary_op << "';" << endl;
            [[fallthrough]];
          case AuxVarType::diff:
            if (aux_vars[i].orig_symb_id)
              output << "M_.aux_vars(" << i+1 << ").orig_index = " << getTypeSpecificID(*aux_vars[i].orig_symb_id)+1 << ";" << endl
                     << "M_.aux_vars(" << i+1 << ").orig_lead_lag = " << aux_vars[i].orig_lead_lag.value() << ";" << endl;
            break;
          case AuxVarType::multiplier:
            output << "M_.aux_vars(" << i+1 << ").eq_nbr = " << aux_vars[i].equation_number_for_multiplier + 1 << ";" << endl;
            break;
          }

        if (expr_t orig_expr = aux_vars[i].expr_node;
            orig_expr)
          {
            output << "M_.aux_vars(" << i+1 << ").orig_expr = '";
            orig_expr->writeJsonOutput(output, {}, {});
            output << "';" << endl;
          }
      }

  if (predeterminedNbr() > 0)
    {
      output << "M_.predetermined_variables = [ ";
      for (int predetermined_variable : predetermined_variables)
        output << getTypeSpecificID(predetermined_variable)+1 << " ";
      output << "];" << endl;
    }

  if (observedVariablesNbr() > 0)
    {
      output << "options_.varobs = cell(" << observedVariablesNbr() << ", 1);" << endl;
      for (int ic{1};
           int it : varobs)
        output << "options_.varobs(" << ic++ << ")  = {'" << getName(it) << "'};" << endl;

      output << "options_.varobs_id = [ ";
      for (int varob : varobs)
        output << getTypeSpecificID(varob)+1 << " ";
      output << " ];"  << endl;
    }

  if (observedExogenousVariablesNbr() > 0)
    {
      output << "options_.varexobs = cell(1);" << endl;
      for (int ic{1};
           int it : varexobs)
        output << "options_.varexobs(" << ic++ << ")  = {'" << getName(it) << "'};" << endl;

      output << "options_.varexobs_id = [ ";
      for (int varexob : varexobs)
        output << getTypeSpecificID(varexob)+1 << " ";
      output << " ];"  << endl;
    }
}

int
SymbolTable::addLeadAuxiliaryVarInternal(bool endo, int index, expr_t expr_arg) noexcept(false)
{
  string varname{(endo ? "AUX_ENDO_LEAD_" : "AUX_EXO_LEAD_") + to_string(index)};
  int symb_id;
  try
    {
      symb_id = addSymbol(varname, SymbolType::endogenous);
    }
  catch (AlreadyDeclaredException &e)
    {
      cerr << "ERROR: you should rename your variable called " << varname << ", this name is internally used by Dynare" << endl;
      exit(EXIT_FAILURE);
    }

  aux_vars.push_back(AuxVarInfo{symb_id, (endo ? AuxVarType::endoLead : AuxVarType::exoLead), 0, 0, 0, 0, expr_arg, ""});

  return symb_id;
}

int
SymbolTable::addLagAuxiliaryVarInternal(bool endo, int orig_symb_id, int orig_lead_lag, expr_t expr_arg) noexcept(false)
{
  string varname{(endo ? "AUX_ENDO_LAG_" : "AUX_EXO_LAG_") + to_string(orig_symb_id) + "_" + to_string(-orig_lead_lag)};
  int symb_id;
  try
    {
      symb_id = addSymbol(varname, SymbolType::endogenous);
    }
  catch (AlreadyDeclaredException &e)
    {
      cerr << "ERROR: you should rename your variable called " << varname << ", this name is internally used by Dynare" << endl;
      exit(EXIT_FAILURE);
    }

  aux_vars.push_back(AuxVarInfo{symb_id, (endo ? AuxVarType::endoLag : AuxVarType::exoLag), orig_symb_id, orig_lead_lag, 0, 0, expr_arg, ""});

  return symb_id;
}

int
SymbolTable::addEndoLeadAuxiliaryVar(int index, expr_t expr_arg) noexcept(false)
{
  return addLeadAuxiliaryVarInternal(true, index, expr_arg);
}

int
SymbolTable::addEndoLagAuxiliaryVar(int orig_symb_id, int orig_lead_lag, expr_t expr_arg) noexcept(false)
{
  return addLagAuxiliaryVarInternal(true, orig_symb_id, orig_lead_lag, expr_arg);
}

int
SymbolTable::addExoLeadAuxiliaryVar(int index, expr_t expr_arg) noexcept(false)
{
  return addLeadAuxiliaryVarInternal(false, index, expr_arg);
}

int
SymbolTable::addExoLagAuxiliaryVar(int orig_symb_id, int orig_lead_lag, expr_t expr_arg) noexcept(false)
{
  return addLagAuxiliaryVarInternal(false, orig_symb_id, orig_lead_lag, expr_arg);
}

int
SymbolTable::addExpectationAuxiliaryVar(int information_set, int index, expr_t expr_arg) noexcept(false)
{
  string varname{"AUX_EXPECT_"s + (information_set < 0 ? "LAG" : "LEAD") + "_"
    + to_string(abs(information_set)) + "_" + to_string(index)};
  int symb_id;
  try
    {
      symb_id = addSymbol(varname, SymbolType::endogenous);
    }
  catch (AlreadyDeclaredException &e)
    {
      cerr << "ERROR: you should rename your variable called " << varname << ", this name is internally used by Dynare" << endl;
      exit(EXIT_FAILURE);
    }

  aux_vars.push_back(AuxVarInfo{symb_id, AuxVarType::expectation, 0, 0, 0, information_set, expr_arg, ""});

  return symb_id;
}

int
SymbolTable::addLogTransformAuxiliaryVar(int orig_symb_id, int orig_lead_lag, expr_t expr_arg) noexcept(false)
{
  string varname = "LOG_" + getName(orig_symb_id);
  int symb_id;
  try
    {
      symb_id = addSymbol(varname, SymbolType::endogenous);
    }
  catch (AlreadyDeclaredException &e)
    {
      cerr << "ERROR: you should rename your variable called " << varname << ", it conflicts with the auxiliary variable created for representing the log of " << getName(orig_symb_id) << endl;
      exit(EXIT_FAILURE);
    }

  aux_vars.push_back(AuxVarInfo{symb_id, AuxVarType::logTransform, orig_symb_id, orig_lead_lag, 0, 0, expr_arg, ""});

  return symb_id;
}

int
SymbolTable::addDiffLagAuxiliaryVar(int index, expr_t expr_arg, int orig_symb_id, int orig_lag) noexcept(false)
{
  string varname{"AUX_DIFF_LAG_" + to_string(index)};
  int symb_id;
  try
    {
      symb_id = addSymbol(varname, SymbolType::endogenous);
    }
  catch (AlreadyDeclaredException &e)
    {
      cerr << "ERROR: you should rename your variable called " << varname << ", this name is internally used by Dynare" << endl;
      exit(EXIT_FAILURE);
    }

  aux_vars.push_back(AuxVarInfo{symb_id, AuxVarType::diffLag, orig_symb_id, orig_lag, 0, 0, expr_arg, ""});

  return symb_id;
}

int
SymbolTable::addDiffLeadAuxiliaryVar(int index, expr_t expr_arg, int orig_symb_id, int orig_lead) noexcept(false)
{
  string varname{"AUX_DIFF_LEAD_" + to_string(index)};
  int symb_id;
  try
    {
      symb_id = addSymbol(varname, SymbolType::endogenous);
    }
  catch (AlreadyDeclaredException &e)
    {
      cerr << "ERROR: you should rename your variable called " << varname << ", this name is internally used by Dynare" << endl;
      exit(EXIT_FAILURE);
    }

  aux_vars.push_back(AuxVarInfo{symb_id, AuxVarType::diffLead, orig_symb_id, orig_lead, 0, 0, expr_arg, ""});

  return symb_id;
}

int
SymbolTable::addDiffAuxiliaryVar(int index, expr_t expr_arg, optional<int> orig_symb_id, optional<int> orig_lag) noexcept(false)
{
  string varname{"AUX_DIFF_" + to_string(index)};
  int symb_id;
  try
    {
      symb_id = addSymbol(varname, SymbolType::endogenous);
    }
  catch (AlreadyDeclaredException &e)
    {
      cerr << "ERROR: you should rename your variable called " << varname << ", this name is internally used by Dynare" << endl;
      exit(EXIT_FAILURE);
    }

  aux_vars.push_back(AuxVarInfo{symb_id, AuxVarType::diff, move(orig_symb_id), move(orig_lag), 0, 0, expr_arg, ""});

  return symb_id;
}

int
SymbolTable::addUnaryOpAuxiliaryVar(int index, expr_t expr_arg, string unary_op, optional<int> orig_symb_id, optional<int> orig_lag) noexcept(false)
{
  string varname{"AUX_UOP_" + to_string(index)};
  int symb_id;
  try
    {
      symb_id = addSymbol(varname, SymbolType::endogenous);
    }
  catch (AlreadyDeclaredException &e)
    {
      cerr << "ERROR: you should rename your variable called " << varname << ", this name is internally used by Dynare" << endl;
      exit(EXIT_FAILURE);
    }

  aux_vars.push_back(AuxVarInfo{symb_id, AuxVarType::unaryOp, move(orig_symb_id), move(orig_lag), 0, 0, expr_arg, unary_op});

  return symb_id;
}

int
SymbolTable::addMultiplierAuxiliaryVar(int index) noexcept(false)
{
  string varname{"MULT_" + to_string(index+1)};
  int symb_id;
  try
    {
      symb_id = addSymbol(varname, SymbolType::endogenous);
    }
  catch (AlreadyDeclaredException &e)
    {
      cerr << "ERROR: you should rename your variable called " << varname << ", this name is internally used by Dynare" << endl;
      exit(EXIT_FAILURE);
    }

  aux_vars.push_back(AuxVarInfo{symb_id, AuxVarType::multiplier, 0, 0, index, 0, nullptr, ""});
  return symb_id;
}

int
SymbolTable::addDiffForwardAuxiliaryVar(int orig_symb_id, int orig_lead_lag, expr_t expr_arg) noexcept(false)
{
  string varname{"AUX_DIFF_FWRD_" + to_string(orig_symb_id+1)};
  int symb_id;
  try
    {
      symb_id = addSymbol(varname, SymbolType::endogenous);
    }
  catch (AlreadyDeclaredException &e)
    {
      cerr << "ERROR: you should rename your variable called " << varname << ", this name is internally used by Dynare" << endl;
      exit(EXIT_FAILURE);
    }

  aux_vars.push_back(AuxVarInfo{symb_id, AuxVarType::diffForward, orig_symb_id, orig_lead_lag, 0, 0, expr_arg, ""});
  return symb_id;
}

int
SymbolTable::addPacExpectationAuxiliaryVar(const string &name, expr_t expr_arg)
{
  int symb_id;
  try
    {
      symb_id = addSymbol(name, SymbolType::endogenous);
    }
  catch (AlreadyDeclaredException &e)
    {
      cerr << "ERROR: the variable/parameter '" << name << "' conflicts with a variable that will be generated for a 'pac_expectation' expression. Please rename it." << endl;
      exit(EXIT_FAILURE);
    }

  aux_vars.push_back(AuxVarInfo{symb_id, AuxVarType::pacExpectation, 0, 0, 0, 0, expr_arg, ""});
  return symb_id;
}

int
SymbolTable::addPacTargetNonstationaryAuxiliaryVar(const string &name, expr_t expr_arg)
{
  int symb_id;
  try
    {
      symb_id = addSymbol(name, SymbolType::endogenous);
    }
  catch (AlreadyDeclaredException &e)
    {
      cerr << "ERROR: the variable/parameter '" << name << "' conflicts with a variable that will be generated for a 'pac_target_nonstationary' expression. Please rename it." << endl;
      exit(EXIT_FAILURE);
    }

  aux_vars.push_back(AuxVarInfo{symb_id, AuxVarType::pacTargetNonstationary, 0, 0, 0, 0, expr_arg, ""});
  return symb_id;
}

int
SymbolTable::searchAuxiliaryVars(int orig_symb_id, int orig_lead_lag) const noexcept(false)
{
  for (const auto &aux_var : aux_vars)
    if ((aux_var.type == AuxVarType::endoLag || aux_var.type == AuxVarType::exoLag)
        && aux_var.orig_symb_id == orig_symb_id && aux_var.orig_lead_lag == orig_lead_lag)
      return aux_var.symb_id;
  throw SearchFailedException{orig_symb_id, orig_lead_lag};
}

int
SymbolTable::getOrigSymbIdForAuxVar(int aux_var_symb_id_arg) const noexcept(false)
{
  for (const auto &aux_var : aux_vars)
    if ((aux_var.type == AuxVarType::endoLag
         || aux_var.type == AuxVarType::exoLag
         || aux_var.type == AuxVarType::diff
         || aux_var.type == AuxVarType::diffLag
         || aux_var.type == AuxVarType::diffLead
         || aux_var.type == AuxVarType::diffForward
         || aux_var.type == AuxVarType::unaryOp)
        && aux_var.symb_id == aux_var_symb_id_arg)
      {
        if (optional<int> r = aux_var.orig_symb_id; r)
          return *r;
        else
          throw UnknownSymbolIDException{aux_var_symb_id_arg}; // Some diff and unaryOp auxvars have orig_symb_id unset
      }
  throw UnknownSymbolIDException{aux_var_symb_id_arg};
}

pair<int, int>
SymbolTable::unrollDiffLeadLagChain(int symb_id, int lag) const noexcept(false)
{
  for (const auto &aux_var : aux_vars)
    if (aux_var.symb_id == symb_id)
      if (aux_var.type == AuxVarType::diffLag || aux_var.type == AuxVarType::diffLead)
        {
          auto [orig_symb_id, orig_lag] = unrollDiffLeadLagChain(aux_var.orig_symb_id.value(), lag);
          return { orig_symb_id, orig_lag + aux_var.orig_lead_lag.value() };
        }
  return { symb_id, lag };
}

void
SymbolTable::markPredetermined(int symb_id) noexcept(false)
{
  validateSymbID(symb_id);

  if (frozen)
    throw FrozenException();

  assert(getType(symb_id) == SymbolType::endogenous);

  predetermined_variables.insert(symb_id);
}

void
SymbolTable::markWithLogTransform(int symb_id) noexcept(false)
{
  validateSymbID(symb_id);

  if (frozen)
    throw FrozenException();

  assert(getType(symb_id) == SymbolType::endogenous);

  with_log_transform.insert(symb_id);
}

bool
SymbolTable::isPredetermined(int symb_id) const noexcept(false)
{
  validateSymbID(symb_id);
  return predetermined_variables.contains(symb_id);
}

int
SymbolTable::predeterminedNbr() const
{
  return predetermined_variables.size();
}

void
SymbolTable::addObservedVariable(int symb_id) noexcept(false)
{
  validateSymbID(symb_id);
  assert(getType(symb_id) == SymbolType::endogenous);
  varobs.push_back(symb_id);
}

int
SymbolTable::observedVariablesNbr() const
{
  return static_cast<int>(varobs.size());
}

bool
SymbolTable::isObservedVariable(int symb_id) const
{
  return find(varobs.begin(), varobs.end(), symb_id) != varobs.end();
}

int
SymbolTable::getObservedVariableIndex(int symb_id) const
{
  auto it = find(varobs.begin(), varobs.end(), symb_id);
  assert(it != varobs.end());
  return static_cast<int>(it - varobs.begin());
}

void
SymbolTable::addObservedExogenousVariable(int symb_id) noexcept(false)
{
  validateSymbID(symb_id);
  assert(getType(symb_id) != SymbolType::endogenous);
  varexobs.push_back(symb_id);
}

int
SymbolTable::observedExogenousVariablesNbr() const
{
  return static_cast<int>(varexobs.size());
}

bool
SymbolTable::isObservedExogenousVariable(int symb_id) const
{
  return find(varexobs.begin(), varexobs.end(), symb_id) != varexobs.end();
}

int
SymbolTable::getObservedExogenousVariableIndex(int symb_id) const
{
  auto it = find(varexobs.begin(), varexobs.end(), symb_id);
  assert(it != varexobs.end());
  return static_cast<int>(it - varexobs.begin());
}

vector <int>
SymbolTable::getTrendVarIds() const
{
  vector <int> trendVars;
  for (const auto &it : symbol_table)
    if (getType(it.second) == SymbolType::trend || getType(it.second) == SymbolType::logTrend)
      trendVars.push_back(it.second);
  return trendVars;
}

set<int>
SymbolTable::getExogenous() const
{
  set <int> exogs;
  for (const auto &it : symbol_table)
    if (getType(it.second) == SymbolType::exogenous)
      exogs.insert(it.second);
  return exogs;
}

set<int>
SymbolTable::getObservedExogenous() const
{
  set <int> oexogs;
  for (const auto &it : symbol_table)
    if (getType(it.second) == SymbolType::exogenous)
      if (isObservedExogenousVariable(it.second))
        oexogs.insert(it.second);
  return oexogs;
}

set<int>
SymbolTable::getEndogenous() const
{
  set <int> endogs;
  for (const auto &it : symbol_table)
    if (getType(it.second) == SymbolType::endogenous)
      endogs.insert(it.second);
  return endogs;
}

bool
SymbolTable::isAuxiliaryVariable(int symb_id) const
{
  return any_of(aux_vars.begin(), aux_vars.end(), [=](const auto &av) { return av.symb_id == symb_id; });
}

bool
SymbolTable::isDiffAuxiliaryVariable(int symb_id) const
{
  return any_of(aux_vars.begin(), aux_vars.end(),
                [=](const auto &av) { return av.symb_id == symb_id
                    && (av.type == AuxVarType::diff
                        || av.type == AuxVarType::diffLag
                        || av.type == AuxVarType::diffLead); });
}

set<int>
SymbolTable::getOrigEndogenous() const
{
  set <int> origendogs;
  for (const auto &it : symbol_table)
    if (getType(it.second) == SymbolType::endogenous && !isAuxiliaryVariable(it.second))
      origendogs.insert(it.second);
  return origendogs;
}

void
SymbolTable::writeJsonOutput(ostream &output) const
{
  output << R"("endogenous": )";
  writeJsonVarVector(output, endo_ids);

  output << R"(, "exogenous":)";
  writeJsonVarVector(output, exo_ids);

  output << R"(, "exogenous_deterministic": )";
  writeJsonVarVector(output, exo_det_ids);

  output << R"(, "parameters": )";
  writeJsonVarVector(output, param_ids);

  if (observedVariablesNbr() > 0)
    {
      output << R"(, "varobs": [)";
      for (size_t i = 0; i < varobs.size(); i++)
	{
	  if (i != 0)
	    output << ", ";
	  output << R"(")" << getName(varobs[i]) << R"(")";
	}
      output << "]" << endl;
      
      output << R"(, "varobs_ids": [)";
      for (size_t i = 0; i < varobs.size(); i++)
	{
	  if (i != 0)
	    output << ", ";
	  output << getTypeSpecificID(varobs[i])+1;
	}
      output << "]" << endl;
    }

  if (observedExogenousVariablesNbr() > 0)
    {
      output << R"(, "varexobs": [)";
      for (size_t i = 0; i < varexobs.size(); i++)
	{
	  if (i != 0)
	    output << ", ";
	  output << R"(")" << getName(varexobs[i]) << R"(")";
	}
      output << "]" << endl;
      
      output << R"(, "varexobs_ids": [)";
      for (size_t i = 0; i < varexobs.size(); i++)
	{
	  if (i != 0)
	    output << ", ";
	  output << getTypeSpecificID(varexobs[i])+1;
	}
      output << "]" << endl;
    }
  // Write the auxiliary variable table
  output << R"(, "orig_endo_nbr": )" << orig_endo_nbr() << endl;
  if (aux_vars.size() == 0)
    output << R"(, "aux_vars": [])";
  else
    {
      output << R"(, "aux_vars": [)" << endl;
      for (int i = 0; i < static_cast<int>(aux_vars.size()); i++)
	{
	  if (i != 0)
	    output << ", ";
	  output << R"({"endo_index": )" << getTypeSpecificID(aux_vars[i].symb_id)+1
		 << R"(, "type": )" << aux_vars[i].get_type_id();
	  switch (aux_vars[i].type)
	    {
	    case AuxVarType::endoLead:
	    case AuxVarType::exoLead:
	    case AuxVarType::expectation:
	    case AuxVarType::pacExpectation:
            case AuxVarType::pacTargetNonstationary:
	      break;
	    case AuxVarType::endoLag:
	    case AuxVarType::exoLag:
            case AuxVarType::logTransform:
	    case AuxVarType::diffLag:
	    case AuxVarType::diffLead:
	    case AuxVarType::diffForward:
	      output << R"(, "orig_index": )" << getTypeSpecificID(aux_vars[i].orig_symb_id.value())+1
		     << R"(, "orig_lead_lag": )" << aux_vars[i].orig_lead_lag.value();
	      break;
	    case AuxVarType::unaryOp:
              output << R"(, "unary_op": ")" << aux_vars[i].unary_op << R"(")";
              [[fallthrough]];
            case AuxVarType::diff:
	      if (aux_vars[i].orig_symb_id)
		output << R"(, "orig_index": )" << getTypeSpecificID(*aux_vars[i].orig_symb_id)+1
		       << R"(, "orig_lead_lag": )" << aux_vars[i].orig_lead_lag.value();
	      break;
	    case AuxVarType::multiplier:
	      output << R"(, "eq_nbr": )" << aux_vars[i].equation_number_for_multiplier + 1;
	      break;
	    }

	  if (expr_t orig_expr = aux_vars[i].expr_node;
	      orig_expr)
	    {
	      output << R"(, "orig_expr": ")";
	      orig_expr->writeJsonOutput(output, {}, {});
	      output << R"(")";
	    }
	  output << '}' << endl;
	}
      output << "]" << endl;
    }
}

void
SymbolTable::writeJsonVarVector(ostream &output, const vector<int> &varvec) const
{
  output << "[";
  for (size_t i = 0; i < varvec.size(); i++)
    {
      if (i != 0)
        output << ", ";
      output << "{"
             << R"("name":")" << getName(varvec[i]) << R"(", )"
             << R"("texName":")" << boost::replace_all_copy(getTeXName(varvec[i]), R"(\)", R"(\\)") << R"(", )"
             << R"("longName":")" << boost::replace_all_copy(getLongName(varvec[i]), R"(\)", R"(\\)") << R"("})"
             << endl;
    }
  output << "]" << endl;
}

int
SymbolTable::getUltimateOrigSymbID(int symb_id) const
{
  while (isAuxiliaryVariable(symb_id))
    try
      {
        symb_id = getOrigSymbIdForAuxVar(symb_id);
      }
    catch (UnknownSymbolIDException &)
      {
        break;
      }
  return symb_id;
}

optional<int>
SymbolTable::getEquationNumberForMultiplier(int symb_id) const
{
  for (const auto &aux_var : aux_vars)
    if (aux_var.symb_id == symb_id && aux_var.type == AuxVarType::multiplier)
      return aux_var.equation_number_for_multiplier;
  return nullopt;
}

const set<int> &
SymbolTable::getVariablesWithLogTransform() const
{
  return with_log_transform;
}

set<int>
SymbolTable::getLagrangeMultipliers() const
{
  set<int> r;
  for (const auto &aux_var : aux_vars)
    if (aux_var.type == AuxVarType::multiplier)
      r.insert(aux_var.symb_id);
  return r;
}
