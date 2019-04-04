/*
 * Copyright © 2003-2019 Dynare Team
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

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <regex>
#include <filesystem>

#include <cstdlib>

#ifndef PACKAGE_VERSION
# define PACKAGE_VERSION 4.
#endif

#include <unistd.h>

#include "ParsingDriver.hh"
#include "ExtendedPreprocessorTypes.hh"
#include "ConfigFile.hh"

/* Prototype for second part of main function
   Splitting main() in two parts was necessary because ParsingDriver.h and MacroDriver.h can't be
   included simultaneously (because of Bison limitations).
*/
void main2(stringstream &in, const string &basename, bool debug, bool clear_all, bool clear_global,
           bool no_tmp_terms, bool no_log, bool no_warn, bool warn_uninit, bool console,
           bool nograph, bool nointeractive, bool parallel, const ConfigFile &config_file,
           WarningConsolidation &warnings_arg, bool nostrict, bool stochastic, bool check_model_changes,
           bool minimal_workspace, bool compute_xrefs, FileOutputType output_mode,
           LanguageOutputType lang, int params_derivs_order, bool transform_unary_ops,
           JsonOutputPointType json, JsonFileOutputType json_output_mode, bool onlyjson, bool jsonderivsimple,
           const string &mexext, const filesystem::path &matlabroot,
           const filesystem::path &dynareroot, bool onlymodel);

void main1(const string &filename, const string &basename, istream &modfile, bool debug, bool save_macro, string &save_macro_file,
           bool no_line_macro, bool no_empty_line_macro, const vector<pair<string, string>> &defines, const vector<string> &path, stringstream &macro_output);

void
usage()
{
  cerr << "Dynare usage: dynare mod_file [debug] [noclearall] [onlyclearglobals] [savemacro[=macro_file]] [onlymacro] [nolinemacro] [noemptylinemacro] [notmpterms] [nolog] [warn_uninit]"
       << " [console] [nograph] [nointeractive] [parallel[=cluster_name]] [conffile=parallel_config_path_and_filename] [parallel_slave_open_mode] [parallel_test]"
       << " [-D<variable>[=<value>]] [-I/path] [nostrict] [stochastic] [fast] [minimal_workspace] [compute_xrefs] [output=dynamic|first|second|third] [language=matlab|julia]"
       << " [params_derivs_order=0|1|2] [transform_unary_ops]"
       << " [json=parse|check|transform|compute] [jsonstdout] [onlyjson] [jsonderivsimple] [nopathchange] [nopreprocessoroutput]"
       << " [mexext=<extension>] [matlabroot=<path>] [onlymodel]"
       << endl;
  exit(EXIT_FAILURE);
}

/* Looks for an options list in the first line of the mod file (but rewind
   the input stream afterwards */
vector<string>
parse_options_line(istream &modfile)
{
  vector<string> options;
  string first_nonempty_line;
  regex pat{R"(^\s*//\s*--\+\s*options:([^\+]*)\+--)"};
  smatch matches;

  while (getline(modfile, first_nonempty_line))
    if (first_nonempty_line != "")
      {
        if (regex_search(first_nonempty_line, matches, pat))
          if (matches.size() > 1 && matches[1].matched)
            {
              regex pat2{R"(([^,\s]+))"};
              string s{matches[1]};
              for (sregex_iterator p(s.begin(), s.end(), pat2);
                   p != sregex_iterator{}; ++p)
                options.push_back((*p)[1]);
            }
        break;
      }

  modfile.seekg(0);

  return options;
}

int
main(int argc, char **argv)
{
  /*
    Redirect stderr to stdout.
    Made necessary because MATLAB/Octave can only capture stdout (but not
    stderr), in order to put it in the logfile (see issue #306)
  */
  dup2(STDOUT_FILENO, STDERR_FILENO);

  if (argc < 2)
    {
      cerr << "Missing model file!" << endl;
      usage();
    }

  string filename = argv[1];
  ifstream modfile(filename, ios::binary);
  if (modfile.fail())
    {
      cerr << "ERROR: Could not open file: " << argv[1] << endl;
      exit(EXIT_FAILURE);
    }

  // Create options list, using first line of mod-file and command line
  vector<string> options = parse_options_line(modfile);
  for (int arg = 2; arg < argc; arg++)
    options.emplace_back(argv[arg]);

  // Parse options
  bool clear_all = true;
  bool clear_global = false;
  bool save_macro = false;
  string save_macro_file;
  bool debug = false;
  bool no_tmp_terms = false;
  bool only_macro = false;
  bool no_line_macro = false;
  bool no_empty_line_macro = false;
  bool no_log = false;
  bool no_warn = false;
  int params_derivs_order = 2;
  bool warn_uninit = false;
  bool console = false;
  bool nograph = false;
  bool nointeractive = false;
  string parallel_config_file;
  bool parallel = false;
  string cluster_name;
  bool parallel_slave_open_mode = false;
  bool parallel_test = false;
  bool nostrict = false;
  bool stochastic = false;
  bool check_model_changes = false;
  bool minimal_workspace = false;
  bool compute_xrefs = false;
  bool transform_unary_ops = false;
  vector<pair<string, string>> defines;
  vector<string> path;
  FileOutputType output_mode{FileOutputType::none};
  JsonOutputPointType json{JsonOutputPointType::nojson};
  JsonFileOutputType json_output_mode{JsonFileOutputType::file};
  bool onlyjson = false;
  bool jsonderivsimple = false;
  LanguageOutputType language{LanguageOutputType::matlab};
  string mexext;
  filesystem::path matlabroot;
  filesystem::path dynareroot{argv[0]};
  dynareroot = dynareroot.parent_path();
  dynareroot = dynareroot / ".." / "..";
  bool onlymodel = false;

  for (auto s : options)
    {
      if (s == "debug")
        debug = true;
      else if (s == "noclearall")
        clear_all = false;
      else if (s.substr(0, 19) == "params_derivs_order")
        {
          if (s.length() > 21 || s.at(19) != '='
              || !(s.at(20) == '0' || s.at(20) == '1' || s.at(20) == '2'))
            {
              cerr << "Incorrect syntax for params_derivs_order option" << endl;
              usage();
            }
          params_derivs_order = stoi(s.substr(20));
        }
      else if (s == "onlyclearglobals")
        {
          clear_all = false;
          clear_global = true;
        }
      else if (s == "onlymacro")
        only_macro = true;
      else if (s.substr(0, 9) == "savemacro")
        {
          save_macro = true;
          if (s.length() > 9)
            {
              if (s.length() == 10 || s.at(9) != '=')
                {
                  cerr << "Incorrect syntax for savemacro option" << endl;
                  usage();
                }
              save_macro_file = s.substr(10);
            }
        }
      else if (s == "nolinemacro")
        no_line_macro = true;
      else if (s == "noemptylinemacro")
        no_empty_line_macro = true;
      else if (s == "notmpterms")
        no_tmp_terms = true;
      else if (s == "nolog")
        no_log = true;
      else if (s == "nowarn")
        no_warn = true;
      else if (s == "warn_uninit")
        warn_uninit = true;
      else if (s == "console")
        console = true;
      else if (s == "nograph")
        nograph = true;
      else if (s == "nointeractive")
        nointeractive = true;
      else if (s.substr(0, 8) == "conffile")
        {
          if (s.length() <= 9 || s.at(8) != '=')
            {
              cerr << "Incorrect syntax for conffile option" << endl;
              usage();
            }
          parallel_config_file = s.substr(9);
        }
      else if (s == "parallel_slave_open_mode")
        parallel_slave_open_mode = true;
      else if (s == "parallel_test")
        parallel_test = true;
      else if (s == "nostrict")
        nostrict = true;
      else if (s == "stochastic")
        stochastic = true;
      else if (s == "fast")
        check_model_changes = true;
      else if (s == "minimal_workspace")
        minimal_workspace = true;
      else if (s == "compute_xrefs")
        compute_xrefs = true;
      else if (s == "transform_unary_ops")
        transform_unary_ops = true;
      else if (s.substr(0, 8) == "parallel")
        {
          parallel = true;
          if (s.length() > 8)
            {
              if (s.length() == 9 || s.at(8) != '=')
                {
                  cerr << "Incorrect syntax for parallel option" << endl;
                  usage();
                }
              cluster_name = s.substr(9);
            }
        }
      else if (s.substr(0, 2) == "-D")
        {
          if (s.length() == 2)
            {
              cerr << "Incorrect syntax for command line define: the defined variable "
                   << "must not be separated from -D by whitespace." << endl;
              usage();
            }

          auto equal_index = s.find('=');
          if (equal_index != string::npos)
            defines.emplace_back(s.substr(2, equal_index-2), s.substr(equal_index+1));
          else
            defines.emplace_back(s.substr(2), "1");
        }
      else if (s.substr(0, 2) == "-I")
        {
          if (s.length() == 2)
            {
              cerr << "Incorrect syntax for command line define: the defined variable "
                   << "must not be separated from -I by whitespace." << endl;
              usage();
            }
          path.push_back(s.substr(2));
        }
      else if (s.substr(0, 6) == "output")
        {
          if (s.length() <= 7 || s.at(6) != '=')
            {
              cerr << "Incorrect syntax for output option" << endl;
              usage();
            }

          s.erase(0, 7);

          if (s == "dynamic")
            output_mode = FileOutputType::dynamic;
          else if (s == "first")
            output_mode = FileOutputType::first;
          else if (s == "second")
            output_mode = FileOutputType::second;
          else if (s == "third")
            output_mode = FileOutputType::third;
          else
            {
              cerr << "Incorrect syntax for output option" << endl;
              usage();
            }
        }
      else if (s.substr(0, 8) == "language")
        {
          if (s.length() <= 9 || s.at(8) != '=')
            {
              cerr << "Incorrect syntax for language option" << endl;
              usage();
            }

          s.erase(0, 9);

          if (s == "julia")
            language = LanguageOutputType::julia;
          if (s == "matlab")
            language = LanguageOutputType::matlab;
          else
            {
              cerr << "Incorrect syntax for language option" << endl;
              usage();
            }
        }
      else if (s == "jsonstdout")
        json_output_mode = JsonFileOutputType::standardout;
      else if (s == "onlyjson")
        onlyjson = true;
      else if (s == "nopreprocessoroutput")
        cout.rdbuf(nullptr);
      else if (s == "jsonderivsimple")
        jsonderivsimple = true;
      else if (s.substr(0, 4) == "json")
        {
          if (s.length() <= 5 || s.at(4) != '=')
            {
              cerr << "Incorrect syntax for json option" << endl;
              usage();
            }

          s.erase(0, 5);

          if (s == "parse")
            json = JsonOutputPointType::parsing;
          else if (s == "check")
            json = JsonOutputPointType::checkpass;
          else if (s == "transform")
            json = JsonOutputPointType::transformpass;
          else if (s == "compute")
            json = JsonOutputPointType::computingpass;
          else
            {
              cerr << "Incorrect syntax for json option" << endl;
              usage();
            }
        }
      else if (s.substr(0, 6) == "mexext")
        {
          if (s.length() <= 7 || s.at(6) != '=')
            {
              cerr << "Incorrect syntax for mexext option" << endl;
              usage();
            }
          mexext = s.substr(7);
        }
      else if (s.substr(0, 10) == "matlabroot")
        {
          if (s.length() <= 11 || s.at(10) != '=')
            {
              cerr << "Incorrect syntax for matlabroot option" << endl;
              usage();
            }
          matlabroot = filesystem::path{s.substr(11)};
        }
      else if (s == "onlymodel")
        onlymodel = true;
      else
        {
          cerr << "Unknown option: " << s << endl;
          usage();
        }
    }

  cout << "Starting preprocessing of the model file ..." << endl;

  // Construct basename (i.e. remove file extension if there is one)
  string basename = argv[1];
  size_t pos = basename.find_last_of('.');
  if (pos != string::npos)
    basename.erase(pos);

  WarningConsolidation warnings(no_warn);

  // Process config file
  ConfigFile config_file(parallel, parallel_test, parallel_slave_open_mode, cluster_name);
  config_file.getConfigFileInfo(parallel_config_file);
  config_file.checkPass(warnings);
  config_file.transformPass();

  // If Include option was passed to the [paths] block of the config file, add
  // it to paths before macroprocessing
  vector<string> config_include_paths = config_file.getIncludePaths();
  for (const auto &it : config_include_paths)
    path.push_back(it);

  // Do macro processing
  stringstream macro_output;
  main1(filename, basename, modfile, debug, save_macro, save_macro_file, no_line_macro, no_empty_line_macro,
        defines, path, macro_output);

  if (only_macro)
    return EXIT_SUCCESS;

  // Do the rest
  main2(macro_output, basename, debug, clear_all, clear_global,
        no_tmp_terms, no_log, no_warn, warn_uninit, console, nograph, nointeractive,
        parallel, config_file, warnings, nostrict, stochastic, check_model_changes, minimal_workspace,
        compute_xrefs, output_mode, language, params_derivs_order, transform_unary_ops,
        json, json_output_mode, onlyjson, jsonderivsimple,
        mexext, matlabroot, dynareroot, onlymodel);

  return EXIT_SUCCESS;
}
