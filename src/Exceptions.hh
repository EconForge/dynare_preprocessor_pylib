/*
 * Copyright © 2026 Dynare Team
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

#ifndef EXCEPTIONS_HH
#define EXCEPTIONS_HH

#include <filesystem>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace std;

//! Source location representation independent of Bison headers
struct SourceLocation
{
  string filename;
  int begin_line {1};
  int begin_column {1};
  int end_line {1};
  int end_column {1};

  SourceLocation() = default;

  SourceLocation(string fn, int bl, int bc, int el, int ec)
    : filename {move(fn)}, begin_line {bl}, begin_column {bc}, end_line {el}, end_column {ec}
  {
  }

  //! Construct from any Bison-like location (Dynare::parser::location_type or Tokenizer::parser::location_type)
  template <typename LocationType>
  explicit SourceLocation(const LocationType& loc)
    : filename {loc.begin.filename ? *loc.begin.filename : ""},
      begin_line {loc.begin.line},
      begin_column {loc.begin.column},
      end_line {loc.end.line},
      end_column {loc.end.column}
  {
  }

  [[nodiscard]] string
  to_string() const
  {
    ostringstream stream;
    if (!filename.empty())
      stream << filename << ": ";
    stream << "line " << begin_line;
    if (begin_line == end_line)
      {
        if (begin_column == end_column - 1)
          stream << ", col " << begin_column;
        else
          stream << ", cols " << begin_column << "-" << end_column - 1;
      }
    else
      stream << ", col " << begin_column << " -"
             << " line " << end_line << ", col " << end_column - 1;
    return stream.str();
  }
};

//! Base class for all preprocessor exceptions
class DynareException : public runtime_error
{
public:
  explicit DynareException(const string& msg) : runtime_error {msg}
  {
  }
};

//! Exception tied to a specific location in a source file (.mod or included file)
class SourceFileException : public DynareException
{
protected:
  optional<SourceLocation> location;
  string message;

  static string
  formatMessage(const optional<SourceLocation>& loc, const string& msg)
  {
    if (loc && (!loc->filename.empty() || loc->begin_line > 0))
      return "ERROR: " + loc->to_string() + ": " + msg;
    return "ERROR: " + msg;
  }

public:
  template <typename LocationType>
  SourceFileException(const LocationType& loc, string msg)
    : DynareException {formatMessage(SourceLocation {loc}, msg)},
      location {SourceLocation {loc}},
      message {move(msg)}
  {
  }

  SourceFileException(SourceLocation loc, string msg)
    : DynareException {formatMessage(loc, msg)},
      location {move(loc)},
      message {move(msg)}
  {
  }

  explicit SourceFileException(string msg)
    : DynareException {msg.starts_with("ERROR: ") ? msg : "ERROR: " + msg},
      location {nullopt},
      message {move(msg)}
  {
  }

  [[nodiscard]] const optional<SourceLocation>&
  getLocation() const noexcept
  {
    return location;
  }

  [[nodiscard]] const string&
  getMessage() const noexcept
  {
    return message;
  }
};

//! Parser-level exceptions (syntax errors, declaration errors, undeclared variables)
class ParserException : public SourceFileException
{
private:
  vector<pair<string, string>> undeclared_variables;

public:
  using SourceFileException::SourceFileException;

  ParserException(vector<pair<string, string>> undeclared_vars, string formatted_msg)
    : SourceFileException {move(formatted_msg)},
      undeclared_variables {move(undeclared_vars)}
  {
  }

  [[nodiscard]] const vector<pair<string, string>>&
  getUndeclaredVariables() const noexcept
  {
    return undeclared_variables;
  }
};

//! Macroprocessor exceptions
class MacroException : public SourceFileException
{
private:
  string backtrace;

public:
  MacroException(string msg, string bt = "")
    : SourceFileException {move(msg)}, backtrace {move(bt)}
  {
  }

  template <typename LocationType>
  MacroException(const LocationType& loc, string msg, string bt = "")
    : SourceFileException {loc, "in macro-processor: " + msg}, backtrace {move(bt)}
  {
  }

  [[nodiscard]] const string&
  getBacktrace() const noexcept
  {
    return backtrace;
  }
};

//! Semantic exceptions occurring during model checking and transformations
class ModelSemanticException : public DynareException
{
protected:
  optional<int> equation_number;
  optional<int> equation_lineno;
  optional<string> equation_tag;
  optional<string> symbol_name;
  string message;

  static string
  formatMessage(const string& msg,
                const optional<int>& eq_num,
                const optional<int>& eq_line,
                const optional<string>& eq_tag,
                const optional<string>& sym_name)
  {
    ostringstream stream;
    stream << "ERROR: ";
    if (eq_num)
      {
        stream << "in equation " << *eq_num;
        if (eq_line)
          stream << " (line " << *eq_line << ")";
        if (eq_tag)
          stream << " ['" << *eq_tag << "']";
        stream << ": ";
      }
    else if (eq_line)
      stream << "at line " << *eq_line << ": ";

    if (sym_name)
      stream << "symbol '" << *sym_name << "': ";

    stream << msg;
    return stream.str();
  }

public:
  explicit ModelSemanticException(string msg,
                                 optional<int> eq_num = nullopt,
                                 optional<int> eq_line = nullopt,
                                 optional<string> eq_tag = nullopt,
                                 optional<string> sym_name = nullopt)
    : DynareException {formatMessage(msg, eq_num, eq_line, eq_tag, sym_name)},
      equation_number {eq_num},
      equation_lineno {eq_line},
      equation_tag {move(eq_tag)},
      symbol_name {move(sym_name)},
      message {move(msg)}
  {
  }

  [[nodiscard]] const optional<int>&
  getEquationNumber() const noexcept
  {
    return equation_number;
  }

  [[nodiscard]] const optional<int>&
  getEquationLineno() const noexcept
  {
    return equation_lineno;
  }

  [[nodiscard]] const optional<string>&
  getEquationTag() const noexcept
  {
    return equation_tag;
  }

  [[nodiscard]] const optional<string>&
  getSymbolName() const noexcept
  {
    return symbol_name;
  }

  [[nodiscard]] const string&
  getMessage() const noexcept
  {
    return message;
  }
};

//! Specialized equation exception
class EquationException : public ModelSemanticException
{
public:
  using ModelSemanticException::ModelSemanticException;
};

//! Statement and computing task configuration exceptions
class StatementException : public DynareException
{
protected:
  string statement_name;
  optional<string> option_name;
  string message;

  static string
  formatMessage(const string& stmt, const optional<string>& opt, const string& msg)
  {
    ostringstream stream;
    stream << "ERROR: ";
    if (!stmt.empty())
      stream << stmt << ": ";
    if (opt)
      stream << "option '" << *opt << "': ";
    stream << msg;
    return stream.str();
  }

public:
  StatementException(string stmt, string msg, optional<string> opt = nullopt)
    : DynareException {formatMessage(stmt, opt, msg)},
      statement_name {move(stmt)},
      option_name {move(opt)},
      message {move(msg)}
  {
  }

  [[nodiscard]] const string&
  getStatementName() const noexcept
  {
    return statement_name;
  }

  [[nodiscard]] const optional<string>&
  getOptionName() const noexcept
  {
    return option_name;
  }

  [[nodiscard]] const string&
  getMessage() const noexcept
  {
    return message;
  }
};

//! File input / output exception
class FileIOException : public DynareException
{
protected:
  filesystem::path path;
  string action;

public:
  FileIOException(filesystem::path p, string act)
    : DynareException {"ERROR: Can't open file " + p.string() + " for " + act},
      path {move(p)},
      action {move(act)}
  {
  }

  [[nodiscard]] const filesystem::path&
  getPath() const noexcept
  {
    return path;
  }

  [[nodiscard]] const string&
  getAction() const noexcept
  {
    return action;
  }
};

//! Mathematical evaluation exception (division by zero, numerical error)
class MathEvaluationException : public DynareException
{
public:
  using DynareException::DynareException;
};

//! Internal compiler error for impossible cases or unmet assertions
class InternalCompilerException : public DynareException
{
public:
  explicit InternalCompilerException(const string& msg)
    : DynareException {"Internal preprocessor error: " + msg}
  {
  }
};

#endif
