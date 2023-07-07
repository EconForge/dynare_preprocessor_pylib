/*
 * Copyright © 2006-2022 Dynare Team
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

#include "Statement.hh"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#include <boost/xpressive/xpressive.hpp>
#pragma GCC diagnostic pop
#include <utility>

void
Statement::checkPass([[maybe_unused]] ModFileStructure &mod_file_struct,
                     [[maybe_unused]] WarningConsolidation &warnings)
{
}

void
Statement::computingPass([[maybe_unused]] const ModFileStructure &mod_file_struct)
{
}

NativeStatement::NativeStatement(string native_statement_arg) :
  native_statement{move(native_statement_arg)}
{
}

void
NativeStatement::writeOutput(ostream &output, [[maybe_unused]] const string &basename,
                             [[maybe_unused]] bool minimal_workspace) const
{
  using namespace boost::xpressive;
  string date_regex = R"((-?\d+([YyAa]|[Mm]([1-9]|1[0-2])|[Qq][1-4]|[SsHh][1-2])))";
  sregex regex_lookbehind = sregex::compile(R"((?<!\$|\d|[a-zA-Z_]|-|'))" + date_regex);
  sregex regex_dollar = sregex::compile(R"((\$))" + date_regex);

  string ns = regex_replace(native_statement, regex_lookbehind, "dates('$&')");
  ns = regex_replace(ns, regex_dollar, "$2"); //replace $DATE with DATE
  output << ns << endl;
}

void
NativeStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "native")"
         << R"(, "string": ")";

  // A similar code is in VerbatimStatement::writeJsonOutput()
  for (auto ch : native_statement)
    switch (ch)
      {
      case '\b':
	output << R"(\b)";
	break;

      case '\f':
	output << R"(\f)";
	break;

      case '\n':
	output << R"(\n)";
	break;

      case '\r':
	output << R"(\r)";
	break;

      case '\t':
	output << R"(\t)";
	break;

      case '"':
	output << R"(\")";
	break;

      case '\\':
	output << R"(\\)";
        break;

      default:
	output << ch;
	break;
      }

  output << R"("})";
}

VerbatimStatement::VerbatimStatement(string verbatim_statement_arg) :
  verbatim_statement{move(verbatim_statement_arg)}
{
}

void
VerbatimStatement::writeOutput(ostream &output, [[maybe_unused]] const string &basename,
                               [[maybe_unused]] bool minimal_workspace) const
{
  output << verbatim_statement << endl;
}

void
VerbatimStatement::writeJsonOutput(ostream &output) const
{
  output << R"({"statementName": "verbatim")"
         << R"(, "string": ")";

  // A similar code is in NativeStatement::writeJsonOutput()
  for (auto ch : verbatim_statement)
    switch (ch)
      {
      case '\b':
	output << R"(\b)";
	break;

      case '\f':
	output << R"(\f)";
	break;

      case '\n':
	output << R"(\n)";
	break;

      case '\r':
	output << R"(\r)";
	break;

      case '\t':
	output << R"(\t)";
	break;

      case '"':
	output << R"(\")";
	break;

      case '\\':
	output << R"(\\)";
        break;

      default:
	output << ch;
	break;
      }

  output << R"("})";
}

void
OptionsList::writeOutput(ostream &output) const
{
  writeOutputCommon(output, "options_");
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

  writeOutputCommon(output, option_group);
}

void
OptionsList::writeOutputCommon(ostream &output, const string &option_group) const
{
  for (const auto &[name, val] : options)
    std::visit([&, &name = name, &val = val]<class T>(const T &v)
    {
      if constexpr(is_same_v<T, SymbolListVal>)
        v.writeOutput(option_group + "." + name, output);
      else
        {
          output << option_group << "." << name << " = ";
          if constexpr(is_same_v<T, NumVal> || is_same_v<T, DateVal>)
            output << v;
          else if constexpr(is_same_v<T, pair<string, string>>)
            output << '[' << v.first << "; " << v.second << ']';
          else if constexpr(is_same_v<T, StringVal>)
            output << "'" << v << "'";
          else if constexpr(is_same_v<T, vector<int>>)
            {
              if (v.size() > 1)
                {
                  output << '[';
                  for (int it : v)
                    output << it << ";";
                  output << ']';
                }
              else
                output << v.front();
            }
          else if constexpr(is_same_v<T, VecStrVal>)
            {
              if (v.size() > 1)
                {
                  output << '{';
                  for (const auto &it : v)
                    output << "'" << it << "';";
                  output << '}';
                }
              else
                output << v.front();
            }
          else if constexpr(is_same_v<T, VecCellStrVal>)
            {
              /* VecCellStrVal should ideally be merged into VecStrVal.
                 only difference is treatment of v.size==1, where VecStrVal
                 does not add quotes and curly brackets, i.e. allows for type conversion of
                 '2' into the number 2 */
              output << '{';
              for (const auto &it : v)
                output << "'" << it << "';";
              output << '}';
            }
          else if constexpr(is_same_v<T, VecValueVal>)
            {
              /* For historical reason, those vectors are output as row vectors (contrary
                 to vectors of integers which are output as column vectors) */
              output << '[';
              for (const auto &it : v)
                output << it << ',';
              output << ']';
            }
          else if constexpr(is_same_v<T, vector<vector<string>>>)
            {
              // Same remark as for VecValueVal
              output << '{';
              for (const auto &v2 : v)
                {
                  output << '[';
                  for (const auto &it : v2)
                    output << it << ',';
                  output << "], ";
                }
              output << '}';
            }
          else
            static_assert(always_false_v<T>, "Non-exhaustive visitor!");
          output << ";" << endl;
        }
    }, val);
}

void
OptionsList::writeJsonOutput(ostream &output) const
{
  if (empty())
    return;

  output << R"("options": {)";

  for (bool opt_written {false};
       const auto &[name, val] : options)
    {
      if (exchange(opt_written, true))
        output << ", ";
      output << R"(")" << name << R"(": )";
      std::visit([&]<class T>(const T &v)
      {
        if constexpr(is_same_v<T, NumVal>)
          output << v;
        else if constexpr(is_same_v<T, pair<string, string>>)
          output << '[' << v.first << ", " << v.second << ']';
        else if constexpr(is_same_v<T, StringVal> || is_same_v<T, DateVal>)
          output << '"' << v << '"';
        else if constexpr(is_same_v<T, SymbolListVal>)
          {
            output << '{';
            v.writeJsonOutput(output);
            output << '}';
          }
        else if constexpr(is_same_v<T, vector<int>> || is_same_v<T, VecStrVal>
                          || is_same_v<T, VecCellStrVal> || is_same_v<T, VecValueVal>
                          || is_same_v<T, vector<vector<string>>>)
          {
            output << '[';
            for (bool printed_something{false};
                 const auto &it : v)
              {
                if (exchange(printed_something, true))
                  output << ", ";
                if constexpr(is_same_v<T, vector<int>> || is_same_v<T, VecValueVal>)
                  output << it;
                else if constexpr(is_same_v<T, VecStrVal> || is_same_v<T, VecCellStrVal>)
                  output << '"' << it << '"';
                else // vector<vector<string>>
                  {
                    output << '[';
                    for (bool printed_something2{false};
                         const auto &it2 : it)
                      {
                        if (exchange(printed_something2, true))
                          output << ", ";
                        output << it2;
                      }
                    output << ']';
                  }
              }
            output << ']';
          }
        else
          static_assert(always_false_v<T>, "Non-exhaustive visitor!");
      }, val);
    }

  output << "}";
}

void
OptionsList::clear()
{
  options.clear();
}

bool
OptionsList::contains(const string &name) const
{
  return options.contains(name);
}

void
OptionsList::erase(const string &name)
{
  options.erase(name);
}

bool
OptionsList::empty() const
{
  return options.empty();
}
