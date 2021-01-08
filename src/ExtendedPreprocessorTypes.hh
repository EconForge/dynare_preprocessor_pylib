/*
 * Copyright © 2014-2021 Dynare Team
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

#ifndef _EXTENDED_PREPROCESSOR_TYPES_HH
#define _EXTENDED_PREPROCESSOR_TYPES_HH

// Values for the “output” option
enum class OutputType
  {
   standard, // Default value, infer the derivation order from .mod file only
   second, // Output at least 2nd dynamic derivatives
   third, // Output at least 3rd dynamic derivatives
  };

// Values for the “language” option
enum class LanguageOutputType
  {
   matlab, // outputs files for MATLAB/Octave processing
   julia, // outputs files for Julia
  };

enum class JsonFileOutputType
  {
   file, // output JSON files to file
   standardout, // output JSON files to stdout
  };

// Values for the “json” option
enum class JsonOutputPointType
  {
   nojson, // don't output JSON
   parsing, // output JSON after the parsing step
   checkpass, // output JSON after the check pass
   transformpass, // output JSON after the transform pass
   computingpass // output JSON after the computing pass
  };

#endif
