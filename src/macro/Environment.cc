/*
 * Copyright © 2019 Dynare Team
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

#include "Environment.hh"
#include "Expressions.hh"

using namespace macro;

void
Environment::define(VariablePtr var, ExpressionPtr value)
{
  string name = var->getName();
  if (functions.find(name) != functions.end())
    throw StackTrace("Variable " + name + " was previously defined as a function");
  variables[move(name)] = value->eval();
}

void
Environment::define(FunctionPtr func, ExpressionPtr value)
{
  string name = func->getName();
  if (variables.find(name) != variables.end())
    throw StackTrace("Variable " + name + " was previously defined as a variable");
  functions[name] = {move(func), move(value)};
}

ExpressionPtr
Environment::getVariable(const string &name) const
{
  auto it = variables.find(name);
  if (it != variables.end())
    return it->second;

  if (!parent)
    throw StackTrace("Unknown variable " + name);

  return getGlobalEnv()->getVariable(name);
}

tuple<FunctionPtr, ExpressionPtr>
Environment::getFunction(const string &name) const
{
  auto it = functions.find(name);
  if (it != functions.end())
    return it->second;

  if (!parent)
    throw StackTrace("Unknown function " + name);

  return parent->getFunction(name);
}

codes::BaseType
Environment::getType(const string &name)
{
  return getVariable(name)->eval()->getType();
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
Environment::print() const
{
  if (!variables.empty())
    {
      cout << "Macro Variables:" << endl;
      for (auto & it : variables)
        {
          cout << "  " << it.first << " = ";
          getVariable(it.first)->eval()->print(cout);
          cout << endl;
        }
    }

  if (!functions.empty())
    {
      cout << "Macro Functions:" << endl;
      for (auto & it : functions)
        {
          cout << "  ";
          get<0>(it.second)->print(cout);
          cout << " = ";
          get<1>(it.second)->print(cout);
          cout << endl;
        }
    }

  if (parent)
    parent->print();
}
