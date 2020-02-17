/*
 * Copyright © 2006-2020 Dynare Team
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

#include "Statement.hh"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#include <boost/xpressive/xpressive.hpp>
#pragma GCC diagnostic pop
#include <utility>

void
Statement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
}

void
Statement::writeJuliaOutput(ostream &output, const string &basename)
{
}

void
Statement::writeJsonOutput(ostream &output) const
{
}

void
Statement::computingPass()
{
}

NativeStatement::NativeStatement(string native_statement_arg) :
  native_statement{move(native_statement_arg)}
{
}

void
NativeStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  using namespace boost::xpressive;
  string date_regex = R"((-?\d+([YyAa]|[Mm]([1-9]|1[0-2])|[Qq][1-4])))";
  sregex regex_lookbehind = sregex::compile(R"((?<!\$|\d|[a-zA-Z_]|\'))" + date_regex);
  sregex regex_dollar = sregex::compile(R"((\$))"+date_regex);

  string ns = regex_replace(native_statement, regex_lookbehind, "dates('$&')");
  ns = regex_replace(ns, regex_dollar, "$2"); //replace $DATE with DATE
  output << ns << endl;
}

void
NativeStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "native")"
         << R"(, "string": ")" << native_statement << R"(")"
         << "}";
}

VerbatimStatement::VerbatimStatement(string verbatim_statement_arg) :
  verbatim_statement{move(verbatim_statement_arg)}
{
}

void
VerbatimStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << verbatim_statement << endl;
}

void
VerbatimStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "verbatim")"
         << R"(, "string": ")" << verbatim_statement << R"(")"
         << "}";
}

void
OptionsList::writeOutput(ostream &output) const
{
  for (const auto & [name, val] : num_options)
    output << "options_." << name << " = " << val << ";" << endl;

  for (const auto & [name, vals] : paired_num_options)
    output << "options_." << name << " = [" << vals.first << "; "
           << vals.second << "];" << endl;

  for (const auto & [name, val] : string_options)
    output << "options_." << name << " = '" << val << "';" << endl;

  for (const auto & [name, val] : date_options)
    output << "options_." << name << " = " << val << ";" << endl;

  for (const auto & [name, list] : symbol_list_options)
    list.writeOutput("options_." + name, output);

  for (const auto & [name, vals] : vector_int_options)
    {
      output << "options_." << name << " = ";
      if (vals.size() > 1)
        {
          output << "[";
          for (int viit : vals)
            output << viit << ";";
          output << "];" << endl;
        }
      else
        output << vals.front() << ";" << endl;
    }

  for (const auto & [name, vals] : vector_str_options)
    {
      output << "options_." << name << " = ";
      if (vals.size() > 1)
        {
          output << "{";
          for (const auto &viit : vals)
            output << "'" << viit << "';";
          output << "};" << endl;
        }
      else
        output << vals.front() << ";" << endl;
    }
}

void
OptionsList::writeOutput(ostream &output, const string &option_group) const
{
  // Initialize option_group as an empty struct iff the field does not exist!
  if (size_t idx = option_group.find_last_of(".");
      idx != string::npos)
    {
      output << "if ~isfield(" << option_group.substr(0, idx) << ",'" << option_group.substr(idx+1) << "')" << endl;
      output << "    " << option_group << " = struct();" << endl;
      output << "end" << endl;
    }
  else
    output << option_group << " = struct();" << endl;

  for (const auto & [name, val] : num_options)
    output << option_group << "." << name << " = " << val << ";" << endl;

  for (const auto & [name, vals] : paired_num_options)
    output << option_group << "." << name << " = [" << vals.first << "; "
           << vals.second << "];" << endl;

  for (const auto & [name, val] : string_options)
    output << option_group << "." << name << " = '" << val << "';" << endl;

  for (const auto & [name, val] : date_options)
    output << option_group << "." << name << " = " << val << ";" << endl;

  for (const auto & [name, list] : symbol_list_options)
    list.writeOutput(option_group + "." + name, output);

  for (const auto & [name, vals] : vector_int_options)
    {
      output << option_group << "." << name << " = ";
      if (vals.size() > 1)
        {
          output << "[";
          for (int viit : vals)
            output << viit << ";";
          output << "];" << endl;
        }
      else
        output << vals.front() << ";" << endl;
    }

  for (const auto & [name, vals] : vector_str_options)
    {
      output << option_group << "." << name << " = ";
      if (vals.size() > 1)
        {
          output << "{";
          for (const auto &viit : vals)
            output << "'" << viit << "';";
          output << "};" << endl;
        }
      else
        output << vals.front() << ";" << endl;
    }
}

void
OptionsList::writeJsonOutput(ostream &output) const
{
  if (getNumberOfOptions() == 0)
    return;

  output << R"("options": {)";
  for (auto it = num_options.begin();
       it != num_options.end();)
    {
      output << R"(")"<< it->first << R"(": )" << it->second;
      ++it;
      if (it != num_options.end()
          || !(paired_num_options.empty()
               && string_options.empty()
               && date_options.empty()
               && symbol_list_options.empty()
               && vector_int_options.empty()))
        output << ", ";
    }

  for (auto it = paired_num_options.begin();
       it != paired_num_options.end();)
    {
      output << R"(")"<< it->first << R"(": [)" << it->second.first << " " << it->second.second << "]";
      ++it;
      if (it != paired_num_options.end()
          || !(string_options.empty()
               && date_options.empty()
               && symbol_list_options.empty()
               && vector_int_options.empty()))
        output << ", ";
    }

  for (auto it = string_options.begin();
       it != string_options.end();)
    {
      output << R"(")"<< it->first << R"(": ")" << it->second << R"(")";
      ++it;
      if (it != string_options.end()
          || !(date_options.empty()
               && symbol_list_options.empty()
               && vector_int_options.empty()))
        output << ", ";
    }

  for (auto it = date_options.begin();
       it != date_options.end();)
    {
      output << R"(")"<< it->first << R"(": ")" << it->second << R"(")";
      ++it;
      if (it != date_options.end()
          || !(symbol_list_options.empty()
               && vector_int_options.empty()))
        output << ", ";
    }

  for (auto it = symbol_list_options.begin();
       it != symbol_list_options.end();)
    {
      output << R"(")"<< it->first << R"(": {)";
      it->second.writeJsonOutput(output);
      output << "}";
      ++it;
      if (it != symbol_list_options.end()
          || !vector_int_options.empty())
        output << ", ";
    }

  for (auto it = vector_int_options.begin();
       it != vector_int_options.end();)
    {
      output << R"(")"<< it->first << R"(": [)";
      if (it->second.size() > 1)
        {
          for (auto viit = it->second.begin();
               viit != it->second.end();)
            {
              output << *viit;
              ++viit;
              if (viit != it->second.end())
                output << ", ";
            }
        }
      else
        output << it->second.front() << endl;
      output << "]";
      ++it;
      if (it != vector_int_options.end())
        output << ", ";
    }

  for (auto it = vector_str_options.begin();
       it != vector_str_options.end();)
    {
      output << R"(")"<< it->first << R"(": [)";
      if (it->second.size() > 1)
        {
          for (auto viit = it->second.begin();
               viit != it->second.end();)
            {
              output << R"(")" << *viit << R"(")";
              ++viit;
              if (viit != it->second.end())
                output << ", ";
            }
        }
      else
        output << it->second.front() << endl;
      output << "]";
      ++it;
      if (it != vector_str_options.end())
        output << ", ";
    }

  output << "}";
}

void
OptionsList::clear()
{
  num_options.clear();
  paired_num_options.clear();
  string_options.clear();
  date_options.clear();
  symbol_list_options.clear();
  vector_int_options.clear();
  vector_str_options.clear();
}

int
OptionsList::getNumberOfOptions() const
{
  return num_options.size()
    + paired_num_options.size()
    + string_options.size()
    + date_options.size()
    + symbol_list_options.size()
    + vector_int_options.size()
    + vector_str_options.size();
}
