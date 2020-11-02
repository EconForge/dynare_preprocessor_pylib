/*
 * Copyright © 2019-2020 Dynare Team
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
#include "Driver.hh"

#include <fstream>

using namespace macro;

void
Eval::interpret(ostream &output, vector<filesystem::path> &paths)
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
Include::interpret(ostream &output, vector<filesystem::path> &paths)
{
  using namespace filesystem;
  try
    {
      StringPtr msp = dynamic_pointer_cast<String>(expr->eval());
      if (!msp)
        throw StackTrace("File name does not evaluate to a string");
      path filename = msp->to_string();
      ifstream incfile(filename, ios::binary);
      if (incfile.fail())
        {
          for (const auto &dir : paths)
            {
              incfile = ifstream(dir / filename, ios::binary);
              if (incfile.good())
                break;
            }
          if (incfile.fail())
            {
              ostringstream errmsg;
              errmsg << "   * " << current_path().string() << endl;
              for (const auto &dir : paths)
                errmsg << "   * " << absolute(dir).string() << endl;
              error(StackTrace("@#includepath", "Could not open " + filename.string()
                               +". The following directories were searched:\n" + errmsg.str(), location));
            }
        }
      Driver m(env);
      // Calling `string()` method on filename and filename.stem() because of bug in
      // MinGW 8.3.0 that ignores implicit conversion to string from filename::path.
      // Test if bug exists when version of MinGW is upgraded on Debian runners
      m.parse(filename.string(), filename.stem().string(), incfile, false, {}, paths, output);
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
  printLineInfo(output);
}

void
IncludePath::interpret(ostream &output, vector<filesystem::path> &paths)
{
  using namespace filesystem;
  try
    {
      StringPtr msp = dynamic_pointer_cast<String>(expr->eval());
      if (!msp)
        throw StackTrace("File name does not evaluate to a string");
      path ip = static_cast<string>(*msp);
      if (!is_directory(ip))
        throw StackTrace(ip.string() + " does not evaluate to a valid directory");
      if (!exists(ip))
        warning(StackTrace("@#includepath", ip.string() + " does not exist", location));
      paths.emplace_back(ip);
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
Define::interpret(ostream &output, vector<filesystem::path> &paths)
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
Echo::interpret(ostream &output, vector<filesystem::path> &paths)
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
  printEndLineInfo(output);
}

void
Error::interpret(ostream &output, vector<filesystem::path> &paths)
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
EchoMacroVars::interpret(ostream &output, vector<filesystem::path> &paths)
{
  if (save)
    env.print(output, vars, location.begin.line, true);
  else
    env.print(cout, vars);
  printEndLineInfo(output);
}

void
For::interpret(ostream &output, vector<filesystem::path> &paths)
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
      for (const auto &statement : statements)
        {
          if (printLine)
            {
              statement->printLineInfo(output);
              printLine = false;
            }
          statement->interpret(output, paths);
        }
    }
  printEndLineInfo(output);
}

void
If::interpret(ostream &output, vector<filesystem::path> &paths)
{
  for (const auto & [expr, body] : expr_and_body)
    try
      {
        auto tmp = expr->eval();
        RealPtr dp = dynamic_pointer_cast<Real>(tmp);
        BoolPtr bp = dynamic_pointer_cast<Bool>(tmp);
        if (!bp && !dp)
          error(StackTrace("@#if",
                           "The condition must evaluate to a boolean or a double", location));
        if ((bp && *bp) || (dp && *dp))
          {
            interpretBody(body, output, paths);
            break;
          }
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
  printEndLineInfo(output);
}

void
If::interpretBody(const vector<DirectivePtr> &body, ostream &output, vector<filesystem::path> &paths)
{
  bool printLine = true;
  for (const auto &statement : body)
    {
      if (printLine)
        {
          statement->printLineInfo(output);
          printLine = false;
        }
      statement->interpret(output, paths);
    }
}

void
Ifdef::interpret(ostream &output, vector<filesystem::path> &paths)
{
  for (const auto & [expr, body] : expr_and_body)
    if (VariablePtr vp = dynamic_pointer_cast<Variable>(expr);
        dynamic_pointer_cast<BaseType>(expr)
        || (vp && env.isVariableDefined(vp->getName())))
      {
        interpretBody(body, output, paths);
        break;
      }
  printEndLineInfo(output);
}

void
Ifndef::interpret(ostream &output, vector<filesystem::path> &paths)
{
  for (const auto & [expr, body] : expr_and_body)
    if (VariablePtr vp = dynamic_pointer_cast<Variable>(expr);
        dynamic_pointer_cast<BaseType>(expr)
        || (vp && !env.isVariableDefined(vp->getName())))
      {
        interpretBody(body, output, paths);
        break;
      }
  printEndLineInfo(output);
}
