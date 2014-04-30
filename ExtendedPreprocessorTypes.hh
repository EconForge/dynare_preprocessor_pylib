/*
 * Copyright (C) 2014 Dynare Team
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

enum FileOutputType
  {
    none,                             // outputs files for Matlab/Octave processing
    dynamic,                          // outputs <fname>_dynamic.* and related files
    first,                            // outputs <fname>_first_derivatives.* and related files 
    second,                           // outputs <fname>_first_derivatives.*, <fname>_second_derivatives.* and related files 
    third,                            // outputs <fname>_first_derivatives.*, <fname>_second_derivatives.*, <fname>_third_derivatives.*  and related files 
  };

enum LanguageOutputType
  {
    matlab,                           // outputs files for Matlab/Octave processing
    c,                                // outputs files for C
    cpp,                              // outputs files for C++
    cuda,                             // outputs files for CUDA (not yet implemented)
    julia,                            // outputs files for Julia (not yet implemented)
    python,                           // outputs files for Python (not yet implemented) (not yet implemented)
  };
#endif
