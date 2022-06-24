/*
 * Copyright © 2019-2022 Dynare Team
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

#ifndef _ENVIRONMENT_HH
#define _ENVIRONMENT_HH

#include "ForwardDeclarationsAndEnums.hh"

#include <map>
#include <vector>
#include <optional>

namespace macro
{
  class Environment
  {
  private:
    const Environment *parent{nullptr};
    map<string, ExpressionPtr> variables;
    map<string, tuple<FunctionPtr, ExpressionPtr>> functions;
  public:
    Environment() = default;
    Environment(const Environment *parent_arg) : parent{parent_arg} { }
    void define(VariablePtr var, ExpressionPtr value);
    void define(FunctionPtr func, ExpressionPtr value);
    ExpressionPtr getVariable(const string &name) const;
    tuple<FunctionPtr, ExpressionPtr> getFunction(const string &name) const;
    codes::BaseType getType(const string &name) const;
    bool isVariableDefined(const string &name) const noexcept;
    bool isFunctionDefined(const string &name) const noexcept;
    bool isSymbolDefined(const string &name) const noexcept
    {
      return isVariableDefined(name) || isFunctionDefined(name);
    }
    void print(ostream &output, const vector<string> &vars, const optional<int> &line = nullopt, bool save = false) const;
    void printVariable(ostream &output, const string &name, const optional<int> &line, bool save) const;
    void printFunction(ostream &output, const tuple<FunctionPtr, ExpressionPtr> &function, const optional<int> &line, bool save) const;
    size_t
    size() const noexcept
    {
      return variables.size() + functions.size();
    }
    const Environment *
    getGlobalEnv() const noexcept
    {
      return parent == nullptr ? this : parent->getGlobalEnv();
    }
  };
}
#endif
