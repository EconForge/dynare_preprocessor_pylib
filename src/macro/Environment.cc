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

#include <cassert>

#include "Environment.hh"
#include "Expressions.hh"

using namespace macro;

void
Environment::define(VariablePtr var, ExpressionPtr value)
{
  string name = var->getName();
  if (functions.contains(name))
    throw StackTrace("Variable " + name + " was previously defined as a function");
  variables[move(name)] = value->eval(*this);
}

void
Environment::define(FunctionPtr func, ExpressionPtr value)
{
  string name = func->getName();
  if (variables.contains(name))
    throw StackTrace("Variable " + name + " was previously defined as a variable");
  functions[name] = {move(func), move(value)};
}

ExpressionPtr
Environment::getVariable(const string &name) const
{
  if (auto it = variables.find(name); it != variables.end())
    return it->second;

  if (!parent)
    throw StackTrace("Unknown variable " + name);

  return getGlobalEnv()->getVariable(name);
}

tuple<FunctionPtr, ExpressionPtr>
Environment::getFunction(const string &name) const
{
  if (auto it = functions.find(name); it != functions.end())
    return it->second;

  if (!parent)
    throw StackTrace("Unknown function " + name);

  return parent->getFunction(name);
}

codes::BaseType
Environment::getType(const string &name) const
{
  return getVariable(name)->eval(const_cast<Environment &>(*this))->getType();
}

bool
Environment::isVariableDefined(const string &name) const noexcept
{
  try
    {
      getVariable(name);
      return true;
    }
  catch (StackTrace &ex)
    {
      return false;
    }
}

bool
Environment::isFunctionDefined(const string &name) const noexcept
{
  try
    {
      getFunction(name);
      return true;
    }
  catch (StackTrace &ex)
    {
      return false;
    }
}

void
Environment::print(ostream &output, const vector<string> &vars, const optional<int> &line, bool save) const
{
  if (!save && !variables.empty())
    output << "Macro Variables:" << endl;

  if (vars.empty())
    for (const auto &it : variables)
      printVariable(output, it.first, line, save);
  else
    for (const auto &it : vars)
      if (isVariableDefined(it))
        printVariable(output, it, line, save);

  if (!save && !functions.empty())
    output << "Macro Functions:" << endl;

  if (vars.empty())
    for (const auto &it : functions)
      printFunction(output, it.second, line, save);
  else
    for (const auto &it : vars)
      if (isFunctionDefined(it))
        printFunction(output, functions.at(it), line, save);

  if (parent)
    parent->print(output, vars, line, save);
}

void
Environment::printVariable(ostream &output, const string &name, const optional<int> &line, bool save) const
{
  assert(!save || line);
  output << (save ? "options_.macrovars_line_" + to_string(*line) + "." : "  ")
         << name << " = ";
  getVariable(name)->eval(const_cast<Environment &>(*this))->print(output, save);
  if (save)
    output << ";";
  output << endl;
}

void
Environment::printFunction(ostream &output, const tuple<FunctionPtr, ExpressionPtr> &function, const optional<int> &line, bool save) const
{
  assert(!save || line);
  output << (save ? "options_.macrovars_line_" + to_string(*line) + ".function." : "  ");
  if (save)
    {
      get<0>(function)->printName(output);
      output << " = '";
    }

  get<0>(function)->print(output);
  output << " = ";
  get<1>(function)->print(output);

  if (save)
    output << "';";
  output << endl;
}
