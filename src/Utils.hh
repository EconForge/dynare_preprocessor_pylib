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

#ifndef UTILS_HH
#define UTILS_HH

#include <filesystem>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

using namespace std;

/* Writes the contents of “new_contents” to the file “filename”. However, if
   the file already exists and would not be modified by this operation, then do
   nothing. */
void writeToFileIfModified(stringstream& new_contents, const filesystem::path& filename);

/* Equivalent of MATLAB/Octave’s strsplit, except that it ignores empty
   substring components (MATLAB/Octave adds them to the output); in
   particular, returns an empty vector given an empty string. */
vector<string> strsplit(string_view str, char delim);

/*! Takes a MATLAB/Octave package name (possibly with several levels nested using dots),
  and returns the path to the corresponding filesystem directory.
  In practice the package nesting is used for the planner_objective (stored
  inside +objective subdir). */
filesystem::path packageDir(const string_view& package);

/* Hack for removing a directory that is locked by MATLAB/Windows. The
   directory is renamed before being deleted. The renaming must occur in the
   same directory, otherwise it may fail if the destination is not on the
   same filesystem. This technique is applied recursively to subdirectories.
   Works even if the argument does not exist or is not a directory. */
void remove_directory_with_matlab_lock(const filesystem::path& dir);

#endif
