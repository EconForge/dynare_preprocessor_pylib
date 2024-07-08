/*
 * Copyright © 2019-2024 Dynare Team
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

#ifndef ENVIRONMENT_HH
#define ENVIRONMENT_HH

#include "ForwardDeclarationsAndEnums.hh"

#include <map>
#include <optional>
#include <vector>

namespace macro
{
class Environment
{
private:
  const Environment* parent {nullptr};
  map<string, ExpressionPtr> variables;
  map<string, pair<FunctionPtr, ExpressionPtr>> functions;

public:
  Environment() = default;
  Environment(const Environment* parent_arg) : parent {parent_arg}
  {
  }
  void define(const VariablePtr& var, const ExpressionPtr& value);
  void define(FunctionPtr func, ExpressionPtr value);
  /* The following two functions are not marked [[nodiscard]], because they are used without output
     to check whether they return an exception or not. */
  ExpressionPtr getVariable(const string& name) const; // NOLINT(modernize-use-nodiscard)
  pair<FunctionPtr, ExpressionPtr>                     // NOLINT(modernize-use-nodiscard)
  getFunction(const string& name) const;
  [[nodiscard]] codes::BaseType getType(const string& name) const;
  [[nodiscard]] bool isVariableDefined(const string& name) const noexcept;
  [[nodiscard]] bool isFunctionDefined(const string& name) const noexcept;
  [[nodiscard]] bool
  isSymbolDefined(const string& name) const noexcept
  {
    return isVariableDefined(name) || isFunctionDefined(name);
  }
  void print(ostream& output, const vector<string>& vars, const optional<int>& line = nullopt,
             bool save = false) const;
  void printVariable(ostream& output, const string& name, const optional<int>& line,
                     bool save) const;
  void printFunction(ostream& output, const string& name, const optional<int>& line,
                     bool save) const;
  [[nodiscard]] size_t
  size() const noexcept
  {
    return variables.size() + functions.size();
  }
  [[nodiscard]] const Environment*
  getGlobalEnv() const noexcept
  {
    return parent == nullptr ? this : parent->getGlobalEnv();
  }
};
}
#endif
