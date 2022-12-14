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

#include "Directives.hh"
#include "Driver.hh"

#include <fstream>
#include <utility>

using namespace macro;

void
Eval::interpret(ostream &output, Environment &env, [[maybe_unused]] vector<filesystem::path> &paths)
{
  try
    {
      output << expr->eval(env)->to_string();
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
Include::interpret(ostream &output, Environment &env, vector<filesystem::path> &paths)
{
  using namespace filesystem;
  try
    {
      StringPtr msp = dynamic_pointer_cast<String>(expr->eval(env));
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
              error(StackTrace("@#include", "Could not open " + filename.string()
                               +". The following directories were searched:\n" + errmsg.str(), location));
            }
        }
      Driver m;
      /* Calling `string()` method on filename and filename.stem() because of
         bug in GCC/MinGW 10.2 (shipped in Debian “Bullseye” 11), that fails
         to accept implicit conversion to string from filename::path. See
         https://en.cppreference.com/w/cpp/filesystem/path/native. */
      m.parse(filename.string(), filename.stem().string(), incfile, false, {}, env, paths, output);
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
  printEndLineInfo(output);
}

void
IncludePath::interpret([[maybe_unused]] ostream &output, Environment &env, vector<filesystem::path> &paths)
{
  using namespace filesystem;
  try
    {
      StringPtr msp = dynamic_pointer_cast<String>(expr->eval(env));
      if (!msp)
        throw StackTrace("File name does not evaluate to a string");
#ifdef _WIN32
      /* Trim trailing slashes and backslashes in the path. This is a
         workaround for a GCC/MinGW bug present in version 10.2
         https://gcc.gnu.org/bugzilla/show_bug.cgi?id=88881, that affects
         gcc-mingw-w64 in Debian “Bullseye” 11. It is fixed in GCC 10.3, and
         thus should be fixed in Debian “Bookworm” 12.
         See Madysson/estimation-codes#11. */
      string ipstr = static_cast<string>(*msp);
      while (ipstr.size() > 1 && (ipstr.back() == '/' || ipstr.back() == '\\'))
        ipstr.pop_back();
      path ip{ipstr};
#else
      path ip = static_cast<string>(*msp);
#endif
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
  printEndLineInfo(output);
}

void
Define::interpret([[maybe_unused]] ostream &output, Environment &env, [[maybe_unused]] vector<filesystem::path> &paths)
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
  printEndLineInfo(output);
}

void
Echo::interpret(ostream &output, Environment &env, [[maybe_unused]] vector<filesystem::path> &paths)
{
  try
    {
      cout << "@#echo (" << getLocation() << "): " << expr->eval(env)->to_string() << endl;
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
Error::interpret([[maybe_unused]] ostream &output, Environment &env, [[maybe_unused]] vector<filesystem::path> &paths)
{
  try
    {
      throw StackTrace(expr->eval(env)->to_string());
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
EchoMacroVars::interpret(ostream &output, Environment &env, [[maybe_unused]] vector<filesystem::path> &paths)
{
  if (save)
    env.print(output, vars, location.begin.line, true);
  else
    env.print(cout, vars);
  printEndLineInfo(output);
}

void
For::interpret(ostream &output, Environment &env, vector<filesystem::path> &paths)
{
  ArrayPtr ap;
  try
    {
      ap = dynamic_pointer_cast<Array>(index_vals->eval(env));
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
          statement->interpret(output, env, paths);
        }
    }
  printEndLineInfo(output);
}

void
If::interpret(ostream &output, Environment &env, vector<filesystem::path> &paths)
{
  for (bool first_clause{true};
       const auto &[expr, body] : expr_and_body)
    try
      {
        if ((ifdef || ifndef) && exchange(first_clause, false))
          {
            VariablePtr vp = dynamic_pointer_cast<Variable>(expr);
            if (!vp)
              error(StackTrace(ifdef ? "@#ifdef" : "@#ifndef",
                               "The condition must be a variable name", location));
            if ((ifdef && env.isVariableDefined(vp->getName()))
                || (ifndef && !env.isVariableDefined(vp->getName())))
              {
                interpretBody(body, output, env, paths);
                break;
              }
          }
        else
          {
            auto tmp = expr->eval(env);
            RealPtr dp = dynamic_pointer_cast<Real>(tmp);
            BoolPtr bp = dynamic_pointer_cast<Bool>(tmp);
            if (!bp && !dp)
              error(StackTrace("@#if",
                               "The condition must evaluate to a boolean or a double", location));
            if ((bp && *bp) || (dp && *dp))
              {
                interpretBody(body, output, env, paths);
                break;
              }
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
If::interpretBody(const vector<DirectivePtr> &body, ostream &output, Environment &env, vector<filesystem::path> &paths)
{
  for (bool printLine{true};
       const auto &statement : body)
    {
      if (exchange(printLine, false))
        statement->printLineInfo(output);
      statement->interpret(output, env, paths);
    }
}
