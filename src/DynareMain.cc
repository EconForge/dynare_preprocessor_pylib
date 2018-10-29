/*
 * Copyright (C) 2003-2018 Dynare Team
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

#include <cstdlib>
#include <cstring>
#ifndef PACKAGE_VERSION
# define PACKAGE_VERSION 4.
#endif

#include <unistd.h>

#include <boost/filesystem.hpp>

#include "ParsingDriver.hh"
#include "ExtendedPreprocessorTypes.hh"
#include "ConfigFile.hh"

/* Prototype for second part of main function
   Splitting main() in two parts was necessary because ParsingDriver.h and MacroDriver.h can't be
   included simultaneously (because of Bison limitations).
*/
void main2(stringstream &in, string &basename, bool debug, bool clear_all, bool clear_global,
           bool no_tmp_terms, bool no_log, bool no_warn, bool warn_uninit, bool console,
           bool nograph, bool nointeractive, bool parallel, ConfigFile &config_file,
           WarningConsolidation &warnings_arg, bool nostrict, bool stochastic, bool check_model_changes,
           bool minimal_workspace, bool compute_xrefs, FileOutputType output_mode,
           LanguageOutputType lang, int params_derivs_order, bool transform_unary_ops,
           JsonOutputPointType json, JsonFileOutputType json_output_mode, bool onlyjson, bool jsonderivsimple,
           bool nopreprocessoroutput, const string &mexext, const boost::filesystem::path &matlabroot,
           const boost::filesystem::path &dynareroot);

void main1(string &modfile, string &basename, string &modfiletxt, bool debug, bool save_macro, string &save_macro_file,
           bool no_line_macro, bool no_empty_line_macro, map<string, string> &defines, vector<string> &path, stringstream &macro_output);

void
usage()
{
  cerr << "Dynare usage: dynare mod_file [debug] [noclearall] [onlyclearglobals] [savemacro[=macro_file]] [onlymacro] [nolinemacro] [noemptylinemacro] [notmpterms] [nolog] [warn_uninit]"
       << " [console] [nograph] [nointeractive] [parallel[=cluster_name]] [conffile=parallel_config_path_and_filename] [parallel_slave_open_mode] [parallel_test]"
       << " [-D<variable>[=<value>]] [-I/path] [nostrict] [stochastic] [fast] [minimal_workspace] [compute_xrefs] [output=dynamic|first|second|third] [language=julia]"
       << " [params_derivs_order=0|1|2] [transform_unary_ops]"
       << " [json=parse|check|transform|compute] [jsonstdout] [onlyjson] [jsonderivsimple] [nopathchange] [nopreprocessoroutput]"
       << " [mexext=<extension>] [matlabroot=<path>]"
       << endl;
  exit(EXIT_FAILURE);
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
  map<string, string> defines;
  vector<string> path;
  FileOutputType output_mode{FileOutputType::none};
  JsonOutputPointType json{JsonOutputPointType::nojson};
  JsonFileOutputType json_output_mode{JsonFileOutputType::file};
  bool onlyjson = false;
  bool jsonderivsimple = false;
  LanguageOutputType language{LanguageOutputType::matlab};
  bool nopreprocessoroutput = false;
  string mexext;
  boost::filesystem::path matlabroot;
  boost::filesystem::path dynareroot{argv[0]};
  dynareroot = dynareroot.parent_path();
  dynareroot = dynareroot / ".." / "..";

  // Parse options
  for (int arg = 2; arg < argc; arg++)
    {
      if (!strcmp(argv[arg], "debug"))
        debug = true;
      else if (!strcmp(argv[arg], "noclearall"))
        clear_all = false;
      else if (strlen(argv[arg]) >= 19 && !strncmp(argv[arg], "params_derivs_order", 19))
        {
          if (strlen(argv[arg]) >= 22 || argv[arg][19] != '='
              || !(argv[arg][20] == '0' || argv[arg][20] == '1' || argv[arg][20] == '2'))
            {
              cerr << "Incorrect syntax for params_derivs_order option" << endl;
              usage();
            }
          params_derivs_order = atoi(argv[arg] + 20);
        }
      else if (!strcmp(argv[arg], "onlyclearglobals"))
        {
          clear_all = false;
          clear_global = true;
        }
      else if (!strcmp(argv[arg], "onlymacro"))
        only_macro = true;
      else if (strlen(argv[arg]) >= 9 && !strncmp(argv[arg], "savemacro", 9))
        {
          save_macro = true;
          if (strlen(argv[arg]) > 9)
            {
              if (strlen(argv[arg]) == 10 || argv[arg][9] != '=')
                {
                  cerr << "Incorrect syntax for savemacro option" << endl;
                  usage();
                }
              save_macro_file = string(argv[arg] + 10);
            }
        }
      else if (!strcmp(argv[arg], "nolinemacro"))
        no_line_macro = true;
      else if (!strcmp(argv[arg], "noemptylinemacro"))
        no_empty_line_macro = true;
      else if (!strcmp(argv[arg], "notmpterms"))
        no_tmp_terms = true;
      else if (!strcmp(argv[arg], "nolog"))
        no_log = true;
      else if (!strcmp(argv[arg], "nowarn"))
        no_warn = true;
      else if (!strcmp(argv[arg], "warn_uninit"))
        warn_uninit = true;
      else if (!strcmp(argv[arg], "console"))
        console = true;
      else if (!strcmp(argv[arg], "nograph"))
        nograph = true;
      else if (!strcmp(argv[arg], "nointeractive"))
        nointeractive = true;
      else if (strlen(argv[arg]) >= 8 && !strncmp(argv[arg], "conffile", 8))
        {
          if (strlen(argv[arg]) <= 9 || argv[arg][8] != '=')
            {
              cerr << "Incorrect syntax for conffile option" << endl;
              usage();
            }
          parallel_config_file = string(argv[arg] + 9);
        }
      else if (!strcmp(argv[arg], "parallel_slave_open_mode"))
        parallel_slave_open_mode = true;
      else if (!strcmp(argv[arg], "parallel_test"))
        parallel_test = true;
      else if (!strcmp(argv[arg], "nostrict"))
        nostrict = true;
      else if (!strcmp(argv[arg], "stochastic"))
        stochastic = true;
      else if (!strcmp(argv[arg], "fast"))
        check_model_changes = true;
      else if (!strcmp(argv[arg], "minimal_workspace"))
        minimal_workspace = true;
      else if (!strcmp(argv[arg], "compute_xrefs"))
        compute_xrefs = true;
      else if (!strcmp(argv[arg], "transform_unary_ops"))
        transform_unary_ops = true;
      else if (strlen(argv[arg]) >= 8 && !strncmp(argv[arg], "parallel", 8))
        {
          parallel = true;
          if (strlen(argv[arg]) > 8)
            {
              if (strlen(argv[arg]) == 9 || argv[arg][8] != '=')
                {
                  cerr << "Incorrect syntax for parallel option" << endl;
                  usage();
                }
              cluster_name = string(argv[arg] + 9);
            }
        }
      else if (strlen(argv[arg]) >= 2 && !strncmp(argv[arg], "-D", 2))
        {
          if (strlen(argv[arg]) == 2)
            {
              cerr << "Incorrect syntax for command line define: the defined variable "
                   << "must not be separated from -D by whitespace." << endl;
              usage();
            }

          size_t equal_index = string(argv[arg]).find('=');
          if (equal_index != string::npos)
            {
              string key = string(argv[arg]).erase(equal_index).erase(0, 2);
              defines[key] = string(argv[arg]).erase(0, equal_index+1);
            }
          else
            {
              string key = string(argv[arg]).erase(0, 2);
              defines[key] = "1";
            }
        }
      else if (strlen(argv[arg]) >= 2 && !strncmp(argv[arg], "-I", 2))
        {
          if (strlen(argv[arg]) == 2)
            {
              cerr << "Incorrect syntax for command line define: the defined variable "
                   << "must not be separated from -I by whitespace." << endl;
              usage();
            }
          path.push_back(string(argv[arg]).erase(0, 2));
        }
      else if (strlen(argv[arg]) >= 6 && !strncmp(argv[arg], "output", 6))
        {
          if (strlen(argv[arg]) <= 7 || argv[arg][6] != '=')
            {
              cerr << "Incorrect syntax for output option" << endl;
              usage();
            }
          if (strlen(argv[arg]) == 14 && !strncmp(argv[arg] + 7, "dynamic", 7))
            output_mode = FileOutputType::dynamic;
          else if (strlen(argv[arg]) ==  12 && !strncmp(argv[arg] + 7, "first", 5))
            output_mode = FileOutputType::first;
          else if (strlen(argv[arg]) == 13 && !strncmp(argv[arg] + 7, "second", 6))
            output_mode = FileOutputType::second;
          else if (strlen(argv[arg]) == 12 && !strncmp(argv[arg] + 7, "third", 5))
            output_mode = FileOutputType::third;
          else
            {
              cerr << "Incorrect syntax for output option" << endl;
              usage();
            }
        }
      else if (strlen(argv[arg]) >= 8 && !strncmp(argv[arg], "language", 8))
        {
          if (strlen(argv[arg]) <= 9 || argv[arg][8] != '=')
            {
              cerr << "Incorrect syntax for language option" << endl;
              usage();
            }

          if (strlen(argv[arg]) == 14 && !strncmp(argv[arg] + 9, "julia", 5))
            language = LanguageOutputType::julia;
          else
            {
              // we don't want temp terms in external functions (except Julia)
              no_tmp_terms = true;
              if (strlen(argv[arg]) == 13 && !strncmp(argv[arg] + 9, "cuda", 4))
                language = LanguageOutputType::cuda;
              else if (strlen(argv[arg]) == 15 && !strncmp(argv[arg] + 9, "python", 6))
                language = LanguageOutputType::python;
              else
                {
                  cerr << "Incorrect syntax for language option" << endl;
                  usage();
                }
            }
        }
      else if (!strcmp(argv[arg], "jsonstdout"))
        json_output_mode = JsonFileOutputType::standardout;
      else if (!strcmp(argv[arg], "onlyjson"))
        onlyjson = true;
      else if (!strcmp(argv[arg], "nopreprocessoroutput"))
        nopreprocessoroutput = true;
      else if (!strcmp(argv[arg], "jsonderivsimple"))
        jsonderivsimple = true;
      else if (strlen(argv[arg]) >= 4 && !strncmp(argv[arg], "json", 4))
        {
          if (strlen(argv[arg]) <= 5 || argv[arg][4] != '=')
            {
              cerr << "Incorrect syntax for json option" << endl;
              usage();
            }
          if (strlen(argv[arg]) == 10 && !strncmp(argv[arg] + 5, "parse", 5))
            json = JsonOutputPointType::parsing;
          else if (strlen(argv[arg]) ==  10 && !strncmp(argv[arg] + 5, "check", 5))
            json = JsonOutputPointType::checkpass;
          else if (strlen(argv[arg]) == 14 && !strncmp(argv[arg] + 5, "transform", 9))
            json = JsonOutputPointType::transformpass;
          else if (strlen(argv[arg]) == 12 && !strncmp(argv[arg] + 5, "compute", 7))
            json = JsonOutputPointType::computingpass;
          else
            {
              cerr << "Incorrect syntax for json option" << endl;
              usage();
            }
        }
      else if (strlen(argv[arg]) >= 6 && !strncmp(argv[arg], "mexext", 6))
        {
          string s{argv[arg]};
          if (s.length() <= 7 || s.at(6) != '=')
            {
              cerr << "Incorrect syntax for mexext option" << endl;
              usage();
            }
          s.erase(0, 7);
          mexext = s;
        }
      else if (strlen(argv[arg]) >= 10 && !strncmp(argv[arg], "matlabroot", 10))
        {
          string s{argv[arg]};
          if (s.length() <= 11 || s.at(10) != '=')
            {
              cerr << "Incorrect syntax for matlabroot option" << endl;
              usage();
            }
          s.erase(0, 11);
          matlabroot = boost::filesystem::path{s};
        }
      else
        {
          cerr << "Unknown option: " << argv[arg] << endl;
          usage();
        }
    }

  if (!nopreprocessoroutput)
    cout << "Starting preprocessing of the model file ..." << endl;

  // Construct basename (i.e. remove file extension if there is one)
  string basename = argv[1];
  string modfile, modfiletxt;
  size_t fsc = basename.find_first_of(';');
  if (fsc != string::npos)
    {
      // If a semicolon is found in argv[1], treat it as the text of the modfile
      modfile = "mod_file_passed_as_string.mod";
      basename = "mod_file_passed_as_string";
      modfiletxt = argv[1];
    }
  else
    {
      // If a semicolon is NOT found in argv[1], treat it as the name of the modfile
      modfile = argv[1];
      size_t pos = basename.find_last_of('.');
      if (pos != string::npos)
        basename.erase(pos);

      ifstream modfile(argv[1], ios::binary);
      if (modfile.fail())
        {
          cerr << "ERROR: Could not open file: " << argv[1] << endl;
          exit(EXIT_FAILURE);
        }

      stringstream buffer;
      buffer << modfile.rdbuf();
      modfiletxt = buffer.str();
    }

  WarningConsolidation warnings(no_warn);

  // Process config file
  ConfigFile config_file(parallel, parallel_test, parallel_slave_open_mode, cluster_name);
  config_file.getConfigFileInfo(parallel_config_file);
  config_file.checkPass(warnings);
  config_file.transformPass();

  // If Include option was passed to the [paths] block of the config file, add
  // it to paths before macroprocessing
  vector<string> config_include_paths = config_file.getIncludePaths();
  for (vector<string>::const_iterator it = config_include_paths.begin();
       it != config_include_paths.end(); it++)
    path.push_back(*it);

  // Do macro processing
  stringstream macro_output;
  main1(modfile, basename, modfiletxt, debug, save_macro, save_macro_file, no_line_macro, no_empty_line_macro,
        defines, path, macro_output);

  if (only_macro)
    return EXIT_SUCCESS;

  // Do the rest
  main2(macro_output, basename, debug, clear_all, clear_global,
        no_tmp_terms, no_log, no_warn, warn_uninit, console, nograph, nointeractive,
        parallel, config_file, warnings, nostrict, stochastic, check_model_changes, minimal_workspace,
        compute_xrefs, output_mode, language, params_derivs_order, transform_unary_ops,
        json, json_output_mode, onlyjson, jsonderivsimple, nopreprocessoroutput,
        mexext, matlabroot, dynareroot);

  return EXIT_SUCCESS;
}
