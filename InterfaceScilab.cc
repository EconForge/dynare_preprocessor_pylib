/*
 * Copyright (C) 2006-2008 Dynare Team
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

#include <string>

#include "Interface.hh"

namespace interfaces
{
  string comment()
  {
    return "// ";
  }

  string delete_file(string s)
  {
    return "mdelete " + s;
  }

  string file_exist(string s)
  {
    return "file_exist(" + s + ")";
  }

  string compile(string s)
  {
    return "";
  }

  string function_close()
  {
    return "endfunction\n";
  }

  string function_file_extension()
  {
    return ".sci";
  }

  string strvcat(string s1, string s2)
  {
    return "[" + s1 + ";" + s2 + "]";
  }

  string load_model_function_files(string filename)
  {
    return "getf('" + filename + "_static.sci');\ngetf('" + filename + "_dynamic.sci');\n";
  }
}
