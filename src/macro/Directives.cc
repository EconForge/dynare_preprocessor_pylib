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

#include "Directives.hh"

using namespace macro;

void
Eval::interpret(ostream &output, bool no_line_macro)
{
  try
    {
      output << expr->eval()->to_string();
    }
  catch (StackTrace &ex)
    {
      ex.push("Evaluation", location);
      error(ex);
    }
  catch (exception &e)
    {
      error(StackTrace("Evaluation", e.what(), location));
    }
}

void
Include::interpret(ostream &output, bool no_line_macro)
{
  try
    {
      StringPtr msp = dynamic_pointer_cast<String>(expr->eval());
      if (!msp)
        throw StackTrace("File name does not evaluate to a string");
      name = *msp;
    }
  catch (StackTrace &ex)
    {
      ex.push("@#include", location);
      error(ex);
    }
  catch (exception &e)
    {
      error(StackTrace("@#include", e.what(), location));
    }
}

void
IncludePath::interpret(ostream &output, bool no_line_macro)
{
  try
    {
      StringPtr msp = dynamic_pointer_cast<String>(expr->eval());
      if (!msp)
        throw StackTrace("File name does not evaluate to a string");
      path = *msp;
    }
  catch (StackTrace &ex)
    {
      ex.push("@#includepath", location);
      error(ex);
    }
  catch (exception &e)
    {
      error(StackTrace("@#includepath", e.what(), location));
    }
}

void
Define::interpret(ostream &output, bool no_line_macro)
{
  try
    {
      if (var)
        env.define(var, value);
      else if (func)
        env.define(func, value);
      else
        throw StackTrace("LHS of can be either a variable or a function");
    }
  catch (StackTrace &ex)
    {
      ex.push("@#define", location);
      error(ex);
    }
  catch (exception &e)
    {
      error(StackTrace("@#define", e.what(), location));
    }
}

void
Echo::interpret(ostream &output, bool no_line_macro)
{
  try
    {
      cout << "@#echo (" << getLocation() << "): " << expr->eval()->to_string() << endl;
    }
  catch (StackTrace &ex)
    {
      ex.push("@#echo", location);
      error(ex);
    }
  catch (exception &e)
    {
      error(StackTrace("@#echo", e.what(), location));
    }
  printEndLineInfo(output, no_line_macro);
}

void
Error::interpret(ostream &output, bool no_line_macro)
{
  try
    {
      throw StackTrace(expr->eval()->to_string());
    }
  catch (StackTrace &ex)
    {
      ex.push("@#error", location);
      error(ex);
    }
  catch (exception &e)
    {
      error(StackTrace("@#error", e.what(), location));
    }
}

void
EchoMacroVars::interpret(ostream &output, bool no_line_macro)
{
  if (save)
    env.print(output, location.begin.line, true);
  else
    env.print(cout);
  printEndLineInfo(output, no_line_macro);
}

void
For::interpret(ostream &output, bool no_line_macro)
{
  ArrayPtr ap;
  try
    {
      ap = dynamic_pointer_cast<Array>(index_vals->eval());
      if (!ap)
        throw StackTrace("The index must loop through an array");
    }
  catch (StackTrace &ex)
    {
      ex.push("@#for", location);
      error(ex);
    }
  catch (exception &e)
    {
      error(StackTrace("@#for", e.what(), location));
    }

  for (size_t i = 0; i < ap->size(); i++)
    {
      if (index_vec.size() == 1)
        env.define(index_vec.at(0), ap->at(i));
      else
        {
          BaseTypePtr btp = dynamic_pointer_cast<BaseType>(ap->at(i));
          if (!btp)
            error(StackTrace("@#for", "Unexpected error encountered in for loop", location));

          if (btp->getType() == codes::BaseType::Tuple)
            {
              TuplePtr mtp = dynamic_pointer_cast<Tuple>(btp);
              if (index_vec.size() != mtp->size())
                error(StackTrace("@#for", "Encountered tuple of size " + to_string(mtp->size())
                                 + " but only have " + to_string(index_vec.size()) + " index variables", location));
              else
                for (size_t j = 0; j < index_vec.size(); j++)
                  env.define(index_vec.at(j), mtp->at(j));
            }
        }

      bool printLine = true;
      for (auto & statement : statements)
        {
          if (printLine)
            {
              statement->printLineInfo(output, no_line_macro);
              printLine = false;
            }
          statement->interpret(output, no_line_macro);
        }
    }
  printEndLineInfo(output, no_line_macro);
}

void
If::interpret(ostream &output, bool no_line_macro)
{
  DoublePtr dp;
  BoolPtr bp;
  try
    {
      dp = dynamic_pointer_cast<Double>(condition->eval());
      bp = dynamic_pointer_cast<Bool>(condition->eval());
      if (!bp && !dp)
        error(StackTrace("@#if", "The condition must evaluate to a boolean or a double", location));
    }
  catch (StackTrace &ex)
    {
      ex.push("@#if", location);
      error(ex);
    }
  catch (exception &e)
    {
      error(StackTrace("@#if", e.what(), location));
    }

  if ((bp && *bp) || (dp && *dp))
    loopIf(output, no_line_macro);
  else
    loopElse(output, no_line_macro);
  printEndLineInfo(output, no_line_macro);
}

void
If::loopIf(ostream &output, bool no_line_macro)
{
  bool printLine = !no_line_macro;
  for (auto & statement : if_statements)
    {
      if (printLine)
        {
          statement->printLineInfo(output, no_line_macro);
          printLine = false;
        }
      statement->interpret(output, no_line_macro);
    }
}

void
If::loopElse(ostream &output, bool no_line_macro)
{
  bool printLine = !no_line_macro;
  for (auto & statement : else_statements)
    {
      if (printLine)
        {
          statement->printLineInfo(output, no_line_macro);
          printLine = false;
        }
      statement->interpret(output, no_line_macro);
    }
}

void
Ifdef::interpret(ostream &output, bool no_line_macro)
{
  BaseTypePtr btnp = dynamic_pointer_cast<BaseType>(condition);
  VariablePtr vp = dynamic_pointer_cast<Variable>(condition);
  if (btnp || (vp && env.isVariableDefined(vp->getName())))
    loopIf(output, no_line_macro);
  else
    loopElse(output, no_line_macro);
  printEndLineInfo(output, no_line_macro);
}

void
Ifndef::interpret(ostream &output, bool no_line_macro)
{
  BaseTypePtr btnp = dynamic_pointer_cast<BaseType>(condition);
  VariablePtr vp = dynamic_pointer_cast<Variable>(condition);
  if (!(btnp || (vp && env.isVariableDefined(vp->getName()))))
    loopIf(output, no_line_macro);
  else
    loopElse(output, no_line_macro);
  printEndLineInfo(output, no_line_macro);
}

