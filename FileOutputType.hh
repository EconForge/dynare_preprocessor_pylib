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

#ifndef _FILE_OUTPUT_TYPE_HH
#define _FILE_OUTPUT_TYPE_HH

enum FileOutputType
  {
    none,                             // outputs files for Matlab/Octave processing
    dynamic,                          // outputs <fname>_dynamic.cc and related files
    first,                            // outputs <fname>_first_derivatives and related files 
    second,                           // outputs <fname>_first_derivatives, <fname>_second_derivatives.cc and related files 
    third,                            // outputs <fname>_first_derivatives, <fname>_second_derivatives.cc, <fname>_third_derivatives.cc  and related files 
  };
#endif
