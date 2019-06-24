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

#include "Driver.hh"

#include <fstream>
#include <sstream>

using namespace macro;

void
Driver::parse(const string &file_arg, const string &basename_arg, istream &modfile,
              ostream &output, bool debug, const vector<pair<string, string>> &defines,
              vector<string> paths_arg)
{
#ifdef _WIN32
  string FILESEP = "\\";
#else
  string FILESEP = "/";
#endif

  file = file_arg;
  basename = basename_arg;
  paths = move(paths_arg);

  if (!defines.empty())
    {
      stringstream command_line_defines_with_endl;
      for (auto & define : defines)
        try
          {
            stoi(define.second);
            command_line_defines_with_endl << "@#define " << define.first << " = " << define.second << endl;
          }
        catch (const invalid_argument &)
          {
            if (!define.second.empty() && define.second.at(0) == '[' && define.second.at(define.second.length()-1) == ']')
              // If the input is an array. Issue #1578
              command_line_defines_with_endl << "@#define " << define.first << " = " << define.second << endl;
            else
              command_line_defines_with_endl << "@#define " << define.first << " = \"" << define.second << "\"" << endl;
          }
      Driver m(env, true);
      istream is(command_line_defines_with_endl.rdbuf());
      m.parse("command_line_defines", "command_line_defines", is, output, debug, vector<pair<string, string>>{}, paths);
    }
  stringstream file_with_endl;
  file_with_endl << modfile.rdbuf() << endl;

  lexer = make_unique<TokenizerFlex>(&file_with_endl);
  lexer->set_debug(debug);

  Tokenizer::parser parser(*this);
  parser.set_debug_level(debug);

  // Launch macro-processing
  parser.parse();

  // Interpret parsed statements
  bool printLine = true;
  for (auto & statement : statements)
    {
      if (printLine)
        {
          statement->printLineInfo(output, no_line_macro);
          printLine = false;
        }
      statement->interpret(output, no_line_macro);

      auto ipp = dynamic_pointer_cast<IncludePath>(statement);
      if (ipp)
        {
          auto p = ipp->getPath();
          paths.insert(paths.end(), p.begin(), p.end());
        }

      auto ip = dynamic_pointer_cast<Include>(statement);
      if (ip)
        {
          string filename = ip->getName();
          ifstream incfile(filename, ios::binary);
          if (incfile.fail())
            {
              ostringstream dirs;
              dirs << "." << FILESEP << endl;
              for (const auto & path : paths)
                {
                  string testfile = path + FILESEP + filename;
                  incfile = ifstream(testfile, ios::binary);
                  if (incfile.good())
                    break;
                  dirs << path << endl;
                }
              if (incfile.fail())
                error(statement->getLocation(), "Could not open " + filename +
                      ". The following directories were searched:\n" + dirs.str());
            }

          string basename = filename;
          size_t pos = basename.find_last_of('.');
          if (pos != string::npos)
            basename.erase(pos);

          Driver m(env, no_line_macro);
          m.parse(filename, basename, incfile, output, debug, vector<pair<string, string>>{}, paths);
        }
    }
}

void
Driver::error(const Tokenizer::parser::location_type &location, const string &message) const
{
  cerr << "ERROR in macro-processor: " << location << ": " << message << endl;
  exit(EXIT_FAILURE);
}
